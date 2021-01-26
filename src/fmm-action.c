// ---------------------------------------------------------------------------
// Copyright (c) 2014 Bo Zhang, Jingfang Huang, Nikos P. Pitsinis, Xiaobai Sun
//
// This file is part of recFMM
//
// recFMM is free software: you can redistribute it and/or modify it
// under the terms of GNU General Public Licenses as published by the Free
// Software Foundation, either version 3 of the licenses, or any later version.
//
// recFMM is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FINESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details. 
//
// You should have received a copy of the GNU General Public License along with
// recFMM. If not, see <http://www.gnu.org/license/>.
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <assert.h>
#include <math.h>

#include "cilk.h"
#include "fmm-param.h"
#include "fmm-action.h"
#include "fmm-dag.h"

Real_t *sources_sorted; 
Real_t *charges_sorted; 
Real_t *targets_sorted; 
double *potential_sorted; 
double *field_sorted; 

fmm_dag_t *fmm_dag; 
fmm_param_t *fmm_param; 

int part_thres; 

void fmm_compute(const fmm_config_t *fmm_config, fmm_dag_t *_fmm_dag, 
                 const Real_t *sources, const Real_t *charges, 
                 const Real_t *targets, double *potential, double *field) {
  int nsources = fmm_config->nsources; 
  int ntargets = fmm_config->ntargets; 
  int datatype = fmm_config->datatype; 
  int s        = fmm_config->s; 
#if YUKAWA
  double beta  = fmm_config->beta;
#endif 

  sources_sorted = CALLOC(nsources * 3, sizeof(Real_t));
  charges_sorted = CALLOC(nsources, sizeof(Real_t));
  targets_sorted = CALLOC(ntargets * 3, sizeof(Real_t));
  potential_sorted = CALLOC(ntargets, sizeof(double));
  field_sorted = CALLOC(ntargets * 3, sizeof(double)); 

  assert(sources_sorted != NULL);
  assert(charges_sorted != NULL);
  assert(targets_sorted != NULL); 
  assert(field_sorted != NULL); 

  // Construct FMM parameters
  fmm_param = construct_param(fmm_config, _fmm_dag); 

  fmm_dag = _fmm_dag; 

  // Rearrange input data
  CILK_FOR (int i = 0; i < nsources; i++) {
    int j = fmm_dag->mapsrc[i], i3 = i * 3, j3 = j * 3; 
    sources_sorted[i3]     = sources[j3]; 
    sources_sorted[i3 + 1] = sources[j3 + 1]; 
    sources_sorted[i3 + 2] = sources[j3 + 2]; 
    charges_sorted[i]      = charges[j]; 
  }

  CILK_FOR (int i = 0; i < ntargets; i++) {
    int j = fmm_dag->maptar[i], i3 = i * 3, j3 = j * 3; 
    targets_sorted[i3]     = targets[j3]; 
    targets_sorted[i3 + 1] = targets[j3 + 1]; 
    targets_sorted[i3 + 2] = targets[j3 + 2]; 
  }


  // Spawn aggregate operation on the source tree
  aggregate(fmm_dag->source_root); 

  // Spawn disaggregate operation on the target tree
  fmm_dag->target_root->list5[0] = fmm_dag->source_root; 
  fmm_dag->target_root->nlist5 = 1; 
  fmm_dag->target_root->expansion = CALLOC(fmm_param->pgsz, 
                                           sizeof(double complex)); 
  assert(fmm_dag->target_root->expansion != NULL);

  part_thres = s; 
  CILK_FOR (int i = 0; i < 8; i++) {
    fmm_box_t *cbox = fmm_dag->target_root->child[i]; 
    if (cbox != NULL) {
      cbox->expansion = CALLOC(fmm_param->pgsz, sizeof(double complex)); 
      assert(cbox->expansion != NULL); 
      disaggregate(cbox);
    }
  }

  // Rearrange output in original order
  CILK_FOR (int i = 0; i < ntargets; i++) {
    int j = fmm_dag->maptar[i], i3 = i * 3, j3 = j * 3;
    potential[j]  = potential_sorted[i];
    field[j3]     = field_sorted[i3];
    field[j3 + 1] = field_sorted[i3 + 1];
    field[j3 + 2] = field_sorted[i3 + 2];
  }
}

void fmm_cleanup(void) {
  destruct_param(fmm_param);
  FREE(sources_sorted);
  FREE(charges_sorted);
  FREE(targets_sorted);
  FREE(potential_sorted);
  FREE(field_sorted);
}

void aggregate(fmm_box_t *sbox) {
  sbox->expansion = CALLOC(fmm_param->pgsz + fmm_param->nexpmax*6, 
                           sizeof(double complex)); 
  assert(sbox->expansion != NULL);

  if (sbox->nchild) {
    CILK_FOR (int i = 0; i < 8; i++) {
      fmm_box_t *child = sbox->child[i]; 
      if (child) 
        aggregate(child); 
    }
    multipole_to_multipole(sbox); 
  } else {
    source_to_multipole(sbox); 
  }

  multipole_to_exponential(sbox); 
}

void disaggregate(fmm_box_t *tbox) {
  fmm_box_t *parent = tbox->parent; 
  fmm_box_t **plist1 = parent->list1; 
  fmm_box_t **plist5 = parent->list5; 
  int nplist1 = parent->nlist1; 
  int nplist5 = parent->nlist5; 

  fmm_box_t **list1 = tbox->list1; 
  fmm_box_t **list5 = tbox->list5; 
  int nlist1 = 0, nlist5 = 0; 

  // Go over list-5 of the parent box and determine its own list-5
  for (int i = 0; i < nplist5; i++) {
    fmm_box_t *sbox = plist5[i]; 
    for (int j = 0; j < 8; j++) {
      fmm_box_t *cbox = sbox->child[j]; 
      if (cbox != NULL) {
        if (abs(tbox->idx - cbox->idx) <= 1 &&
            abs(tbox->idy - cbox->idy) <= 1 &&
            abs(tbox->idz - cbox->idz) <= 1) 
          list5[nlist5++] = cbox;
      }
    }
  }

  tbox->nlist5 = nlist5; 

  // Go over list-1 of the parent box. If the box is adjacent, then it
  // belongs to coarser level list-1 and the pointer is
  // stored. Otherwise, the box belongs to list-4 and is processed
  // using source-to-local oeprator.
  for (int i = 0; i < nplist1; i++) {
    fmm_box_t *sbox = plist1[i]; 
    if (is_adjacent(sbox, tbox)) {
      list1[nlist1++] = sbox; 
    } else {
      source_to_local(tbox, sbox); 
    }
  }

  // Check if it is possible to prune the subtree rooted at
  // 'tbox'. First, we check if the box is adajcent to any source box
  // by checking if tbox->nlist5 is zero. If not, we check if any of
  // the list-5 box has more than 'part_thres' points. When both
  // conditions are satifised, the partition performed at 'tbox' will
  // be kept. Otherwise, the lower branch will be dropped. 
  if (tbox->npts <= part_thres || tbox->nlist5 == 0) {
    CILK_FOR (int i = 0; i < 8; i++) {
      delete_box(tbox->child[i]); 
      tbox->child[i] = NULL;
    }
    tbox->nchild = 0;
  } else {
    bool remove = true;
    for (int i = 0; i < nlist5; i++) {
      fmm_box_t *sbox = list5[i]; 
      remove &= (sbox->npts <= part_thres); 
    }

    if (remove) {
      CILK_FOR (int i = 0; i < 8; i++) {
        delete_box(tbox->child[i]); 
        tbox->child[i] = NULL;
      } 
      tbox->nchild = 0;
    }
  }

  if (tbox->nchild) {
    // Complete exponential-to-local translation using the
    // merge-and-shift technique. 
    for (int i = 0; i < 8; i++) {
      fmm_box_t *child = tbox->child[i];
      if (child != NULL) {
        child->expansion = CALLOC(fmm_param->pgsz, sizeof(double complex));
        assert(child->expansion != NULL);
      }
    }

    exponential_to_local(tbox); 
    local_to_local(tbox); 

    // Go over list-5 and if any box in that list is a leaf box,
    // relocate it to list-1 and then spawn work on the children of 'tbox'.
    for (int i = 0; i < nlist5; i++) {
      fmm_box_t *sbox = list5[i]; 
      if (sbox->nchild == 0) 
        list1[nlist1++] = sbox; 
    }
    tbox->nlist1 = nlist1; 

    CILK_FOR (int i = 0; i < 8; i++) {
      fmm_box_t *child = tbox->child[i]; 
      if (child != NULL) 
        disaggregate(child); 
    }
  } else {
    // Evaluate the local expansion to get far-field influence. 
    local_to_target(tbox); 
    
    // If tbox has nonempty list-5, go over each entry in the list-5
    // and process them either as a list-1 entry or a list-3 entry. 
    if (nlist5) {
      for (int i = 0; i < nlist5; i++) {
        fmm_box_t *sbox = list5[i]; 
        process_list13(tbox, sbox); 	
      }
    }

    // Process the coarser or same level list-1 boxes
    for (int i = 0; i < nlist1; i++) {
      fmm_box_t *sbox = list1[i]; 
      direct_evaluation(tbox, sbox);
    }
  }	     
}  

void process_list13(fmm_box_t *tbox, fmm_box_t *sbox) {
  if (is_adjacent(tbox, sbox)) {
    if (sbox->nchild) {
      for (int i = 0; i < 8; i++) {
        fmm_box_t *child = sbox->child[i]; 
        if (child != NULL) 
          process_list13(tbox, child);
      }
    } else {
      direct_evaluation(tbox, sbox); 
    }
  } else {
    // This source box is a list-3 box by definition. If 'sbox'
    // contains more than 'pgsz' points, multipole-to-target operator
    // is the more efficient way to process it. Otherwise, process
    // this source box by direct evaluation. 
    if (sbox->npts > fmm_param->pgsz) {
      multipole_to_target(tbox, sbox);
    } else {
      direct_evaluation(tbox, sbox);
    }
  }
}

void make_ulist(int type, fmm_box_t **list, int nlist, 
                int *xoff, int *yoff, double complex *mexpo, int level) {
  int nexpmax = fmm_param->nexpmax; 
  int pgsz    = fmm_param->pgsz; 

  for (int i = 0; i < nexpmax; i++) 
    mexpo[i] = 0;

#if LAPLACE 
  int nexpo = fmm_param->nexptotp; 
  double complex *xs = fmm_param->xs;
  double complex *ys = fmm_param->ys;
#elif YUKAWA
  int nexpo = fmm_param->mnexptotp;
  double complex *xs = &fmm_param->xs[3 * nexpmax * level];
  double complex *ys = &fmm_param->ys[3 * nexpmax * level]; 
#endif 

  if (nlist) {
    for (int i = 0; i < nlist; i++) {
      fmm_box_t *sbox = list[i]; 
      const double complex *temp = &sbox->expansion[pgsz + type * nexpmax]; 
      for (int j = 0; j < nexpo; j++) {
        double complex zmul = 1;
        if (xoff[i] > 0) {
          zmul *= xs[3 * j + xoff[i] - 1];
        } else if (xoff[i] < 0) {
          zmul *= conj(xs[3 * j - xoff[i] - 1]);
        }

        if (yoff[i] > 0) {
          zmul *= ys[3 * j + yoff[i] - 1];
        } else if (yoff[i] < 0) {
          zmul *= conj(ys[3 * j - yoff[i] - 1]);
        }

        mexpo[j] += zmul * temp[j];
      }
    }
  }
}

void make_dlist(int type, fmm_box_t **list, int nlist, 
                int *xoff, int *yoff, double complex *mexpo, int level) {
  int nexpmax = fmm_param->nexpmax;
  int pgsz    = fmm_param->pgsz; 
  
  for (int i = 0; i < nexpmax; i++) 
    mexpo[i] = 0;

#if LAPLACE 
  int nexpo = fmm_param->nexptotp; 
  double complex *xs = fmm_param->xs;
  double complex *ys = fmm_param->ys;
#elif YUKAWA
  int nexpo = fmm_param->mnexptotp;
  double complex *xs = &fmm_param->xs[3 * nexpmax * level];
  double complex *ys = &fmm_param->ys[3 * nexpmax * level]; 
#endif 

  if (nlist) {
    for (int i = 0; i < nlist; i++) {
      fmm_box_t *sbox = list[i]; 
      const double complex *temp = &sbox->expansion[pgsz + type * nexpmax]; 
      for (int j = 0; j < nexpo; j++) {
        double complex zmul = 1;
        if (xoff[i] > 0) {
          zmul *= conj(xs[3 * j + xoff[i] - 1]);
        } else if (xoff[i] < 0) {
          zmul *= xs[3 * j - xoff[i] - 1];
        }

        if (yoff[i] > 0) {
          zmul *= conj(ys[3 * j + yoff[i] - 1]);
        } else if (yoff[i] < 0) {
          zmul *= ys[3 * j - yoff[i] - 1];
        }

        mexpo[j] += zmul * temp[j];
      }
    }
  }
}

void rotz2y(const double complex *multipole, const double *rd,
            double complex *mrotate) {
  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz;

  double complex *mwork = CALLOC(pgsz, sizeof(double complex));
  double complex *ephi = CALLOC(pterms + 1, sizeof(double complex));

  assert(mwork != NULL);
  assert(ephi != NULL);

  ephi[0] = 1.0;
  for (int m =1; m <= pterms; m++)
    ephi[m] = -ephi[m - 1] * _Complex_I;

  for (int m = 0; m <= pterms; m++) {
    int offset = m * (pterms + 1);
    for (int ell = m; ell <= pterms; ell++) {
      int index = offset + ell;
      mwork[index] = ephi[m] * multipole[index];
    }
  }

  for (int m = 0; m <= pterms; m++) {
    int offset = m * (pterms + 1);
    for (int ell = m; ell <= pterms; ell++) {
      int index = ell + offset;
      mrotate[index] = mwork[ell] * rd[ell + (m + pterms) * pgsz];
      for (int mp = 1; mp <= ell; mp++) {
        int index1 = ell + mp * (pterms + 1);
        mrotate[index] +=
          mwork[index1] * rd[ell + mp * (pterms + 1) + (m + pterms) * pgsz] +
          conj(mwork[index1]) * 
          rd[ell + mp * (pterms + 1) + (pterms - m) * pgsz];
      }
    }
  }
  
  FREE(ephi);
  FREE(mwork);
}

void roty2z(const double complex *multipole, const double *rd,
            double complex *mrotate) {
  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz;

  double complex *mwork = CALLOC(pgsz, sizeof(double complex));
  double complex *ephi = CALLOC(1 + pterms, sizeof(double complex));

  ephi[0] = 1.0;
  for (int m = 1; m <= pterms; m++)
    ephi[m] = ephi[m - 1] * _Complex_I;

  for (int m = 0; m <= pterms; m++) {
    int offset = m * (pterms + 1);
    for (int ell = m; ell <= pterms; ell++) {
      int index = ell + offset;
      mwork[index] = multipole[ell] * rd[ell + (m + pterms) * pgsz];
      for (int mp = 1; mp <= ell; mp++) {
        int index1 = ell + mp * (pterms + 1);
        double complex temp = multipole[index1];
        mwork[index] +=
          temp * rd[ell + mp * (pterms + 1) + (m + pterms) * pgsz] +
          conj(temp) * rd[ell + mp * (pterms + 1) + (pterms - m) * pgsz];
      }
    }
  }

  for (int m = 0; m <= pterms; m++) {
    int offset = m * (pterms + 1);
    for (int ell = m; ell <= pterms; ell++) {
      int index = ell + offset;
      mrotate[index] = ephi[m] * mwork[index];
    }
  }

  FREE(ephi);
  FREE(mwork);
}

void rotz2x(const double complex *multipole, const double *rd,
            double complex *mrotate) {
  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz;

  int offset1 = pterms * pgsz;
  for (int m = 0; m <= pterms; m++) {
    int offset2 = m * (pterms + 1);
    int offset3 = m * pgsz + offset1;
    int offset4 = -m * pgsz + offset1;
    for (int ell = m; ell <= pterms; ell++) {
      mrotate[ell + offset2] = multipole[ell] * rd[ell + offset3];
      for (int mp = 1; mp <= ell; mp++) {
        int offset5 = mp * (pterms + 1);
        mrotate[ell + offset2] +=
          multipole[ell + offset5] * rd[ell + offset3 + offset5] +
          conj(multipole[ell + offset5]) * rd[ell + offset4 + offset5];
      }
    }
  }
}

void lgndr(int nmax, double x, double *y) {
  int n;
  n = (nmax + 1) * (nmax + 1);
  for (int m = 0; m < n; m++)
    y[m] = 0.0;

  double u = -sqrt(1 - x * x);
  y[0] = 1;

  y[1] = x * y[0];
  for (int n = 2; n <= nmax; n++)
    y[n] = ((2 * n - 1) * x * y[n - 1] - (n - 1) * y[n - 2]) / n;

  int offset1 = nmax + 2;
  for (int m = 1; m <= nmax - 1; m++) {
    int offset2 = m * offset1;
    y[offset2] = y[offset2 - offset1] * u * (2 * m - 1);
    y[offset2 + 1] = y[offset2] * x * (2 * m + 1);
    for (int n = m + 2; n <= nmax; n++) {
      int offset3 = n + m * (nmax + 1);
      y[offset3] = ((2 * n - 1) * x * y[offset3 - 1] -
                    (n + m - 1) * y[offset3 - 2]) / (n - m);
    }
  }

  y[nmax + nmax * (nmax + 1)] =
    y[nmax - 1 + (nmax - 1) * (nmax + 1)] * u * (2 * nmax - 1);
}

#if LAPLACE

void source_to_multipole(fmm_box_t *sbox) {
  int pgsz      = fmm_param->pgsz; 
  int pterms    = fmm_param->pterms; 
  double *ytopc = fmm_param->ytopc;
  double *scale = fmm_param->scale; 
  int addr      = sbox->addr; 
  int level     = sbox->level; 

  double complex *multipole = &sbox->expansion[0]; 
  Real_t *sources = &sources_sorted[3 * addr]; 
  Real_t *charges = &charges_sorted[addr]; 
  int nsources = sbox->npts; 
  double myscale = scale[level]; 
  double h = fmm_dag->size / (1 << (level + 1)); 
  double center[3]; 
  center[0] = fmm_dag->corner[0] + (2 * sbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * sbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * sbox->idz + 1) * h;

  const double precision = 1e-14;
  double *powers = CALLOC(pterms + 1, sizeof(double));
  double *p = CALLOC(pgsz, sizeof(double));
  double complex *ephi = CALLOC(pterms + 1, sizeof(double complex));
  assert(powers != NULL);
  assert(p != NULL);
  assert(ephi != NULL); 

  for (int i = 0; i < nsources; i++) {
    int i3 = i * 3; 
    double rx = sources[i3]     - center[0];
    double ry = sources[i3 + 1] - center[1];
    double rz = sources[i3 + 2] - center[2];
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 1.0 : rz / d);
    ephi[0] = (proj <= precision * d ? 1.0 : rx / proj + _Complex_I * ry / proj);
    d *= myscale;
    powers[0] = 1.0;
    for (int ell = 1; ell <= pterms; ell++) {
      powers[ell] = powers[ell - 1] * d;
      ephi[ell] = ephi[ell - 1] * ephi[0];
    }

    double charge = charges[i];
    multipole[0] += charge;

    lgndr(pterms, ctheta, p);
    for (int ell = 1; ell <= pterms; ell++) {
      double cp = charge * powers[ell] * p[ell];
      multipole[ell] += cp;
    }

    for (int m = 1; m <= pterms; m++) {
      int offset1 = m * (pterms + 1);
      int offset2 = m * (pterms + 2); 
      for (int ell = m; ell <= pterms; ell++) {
        double cp = charge * powers[ell] * ytopc[ell + offset2] * 
          p[ell + offset1];
        multipole[ell + offset1] += cp * conj(ephi[m - 1]);
      }
    }
  }

  FREE(powers);
  FREE(p);
  FREE(ephi);
}

void multipole_to_multipole(fmm_box_t *sbox) {
  const double complex var[5] =
    {1,-1 + _Complex_I, 1 + _Complex_I, 1 - _Complex_I, -1 - _Complex_I};
  const double arg = sqrt(2)/2.0;
  const int iflu[8] = {3, 4, 2, 1, 3, 4, 2, 1};

  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz;
  double *dc = fmm_param->dc;

  double *powers = CALLOC(pterms + 3, sizeof(double));
  double complex *mpolen = CALLOC(pgsz, sizeof(double complex));
  double complex *marray = CALLOC(pgsz, sizeof(double complex));
  double complex *ephi   = CALLOC(pterms + 3, sizeof(double complex));
  
  assert(powers != NULL);
  assert(mpolen != NULL);
  assert(marray != NULL);
  assert(ephi != NULL); 

  double complex *pmpole = &sbox->expansion[0]; 

  for (int i = 0; i < 8; i++) {
    fmm_box_t *child = sbox->child[i]; 
    if (child != NULL) {
      int ifl = iflu[i]; 
      double *rd = (i < 4 ? fmm_param->rdsq3 : fmm_param->rdmsq3);
      double complex *cmpole = &child->expansion[0]; 

      ephi[0] = 1.0;
      ephi[1] = arg * var[ifl];
      double dd = -sqrt(3) / 2.0;
      powers[0] = 1.0;

      for (int ell = 1; ell <= pterms + 1; ell++) {
        powers[ell] = powers[ell - 1] * dd;
        ephi[ell + 1] = ephi[ell] * ephi[1];
      }

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset; 
          mpolen[index] = conj(ephi[m]) * cmpole[index];
        }
      }

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = offset + ell;
          marray[index] = mpolen[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp++) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }
        }
      }
      
      for (int k = 0; k <= pterms; k++) {
        int offset = k * (pterms + 1);
        for (int j = k; j <= pterms; j++) {
          int index = offset + j;
          mpolen[index] = marray[index];
          for (int ell = 1; ell <= j - k; ell++) {
            int index2 = j - k + ell * (2 * pterms + 1);
            int index3 = j + k + ell * (2 * pterms + 1);
            mpolen[index] += marray[index - ell] * powers[ell] *
              dc[index2] * dc[index3];
          }
        }
      }
      
      for (int m = 0; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = mpolen[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] -= mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }

          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }
        }
      }

      for (int m = 1; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = -mpolen[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }

          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] -= mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }
        }
      }

      powers[0] = 1.0;
      for (int ell = 1; ell <= pterms + 1; ell++)
        powers[ell] = powers[ell - 1] / 2;

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          mpolen[index] = ephi[m] * marray[index] * powers[ell];
        }
      }
      
      for (int m = 0; m < pgsz; m++)
        pmpole[m] += mpolen[m];
    }
  }

  FREE(ephi);
  FREE(powers);
  FREE(mpolen);
  FREE(marray);
}

void multipole_to_exponential(fmm_box_t *sbox) {
  int pgsz        = fmm_param->pgsz;
  int nexpmax     = fmm_param->nexpmax;
  double *rdminus = fmm_param->rdminus; 
  double *rdplus  = fmm_param->rdplus; 

  double complex *mw = CALLOC(pgsz, sizeof(double complex)); 
  double complex *mexpf1 = CALLOC(nexpmax, sizeof(double complex)); 
  double complex *mexpf2 = CALLOC(nexpmax, sizeof(double complex)); 

  assert(mw != NULL);
  assert(mexpf1 != NULL);
  assert(mexpf2 != NULL);

  multipole_to_exponential_p1(&sbox->expansion[0], mexpf1, mexpf2); 
  multipole_to_exponential_p2(mexpf1, &sbox->expansion[pgsz]); 
  multipole_to_exponential_p2(mexpf2, &sbox->expansion[pgsz + nexpmax]); 

  rotz2y(&sbox->expansion[0], rdminus, mw);
  multipole_to_exponential_p1(mw, mexpf1, mexpf2); 
  multipole_to_exponential_p2(mexpf1, &sbox->expansion[pgsz + nexpmax * 2]); 
  multipole_to_exponential_p2(mexpf2, &sbox->expansion[pgsz + nexpmax * 3]);

  rotz2x(&sbox->expansion[0], rdplus, mw); 
  multipole_to_exponential_p1(mw, mexpf1, mexpf2); 
  multipole_to_exponential_p2(mexpf1, &sbox->expansion[pgsz + nexpmax * 4]); 
  multipole_to_exponential_p2(mexpf2, &sbox->expansion[pgsz + nexpmax * 5]); 

  FREE(mw);
  FREE(mexpf1);
  FREE(mexpf2); 
}

void multipole_to_exponential_p1(const double complex *multipole, 
                                 double complex *mexpu, 
                                 double complex *mexpd) {
  int nlambs   = fmm_param->nlambs;
  int *numfour = fmm_param->numfour;
  int pterms   = fmm_param->pterms;
  int pgsz     = fmm_param->pgsz;
  double *rlsc = fmm_param->rlsc;

  int ntot = 0;
  for (int nell = 0; nell < nlambs; nell++) {
    double sgn = -1.0;
    double complex zeyep = 1.0;
    for (int mth = 0; mth <= numfour[nell] - 1; mth++) {
      int ncurrent = ntot + mth;
      double complex ztmp1 = 0.0;
      double complex ztmp2 = 0.0;
      sgn = -sgn;
      int offset = mth * (pterms + 1);
      int offset1 = offset + nell * pgsz;
      for (int nm = mth; nm <= pterms; nm += 2)
        ztmp1 += rlsc[nm + offset1] * multipole[nm + offset];
      for (int nm = mth + 1; nm <= pterms; nm += 2)
        ztmp2 += rlsc[nm + offset1] * multipole[nm + offset];

      mexpu[ncurrent] = (ztmp1 + ztmp2) * zeyep;
      mexpd[ncurrent] = sgn * (ztmp1 - ztmp2) * zeyep;
      zeyep *= _Complex_I;
    }
    ntot += numfour[nell];
  }
}

void multipole_to_exponential_p2(const double complex *mexpf, 
                                 double complex *mexpphys) {
  int nlambs   = fmm_param->nlambs;
  int *numfour = fmm_param->numfour;
  int *numphys = fmm_param->numphys;

  int nftot, nptot, nexte, nexto;
  nftot = 0;
  nptot = 0;
  nexte = 0;
  nexto = 0;

  double complex *fexpe = fmm_param->fexpe;
  double complex *fexpo = fmm_param->fexpo;

  for (int i = 0; i < nlambs; i++) {
    for (int ival = 0; ival < numphys[i] / 2; ival++) {
      mexpphys[nptot + ival] = mexpf[nftot]; 

      for (int nm = 1; nm < numfour[i]; nm += 2) {
        double rt1 = cimag(fexpe[nexte]) * creal(mexpf[nftot + nm]);
        double rt2 = creal(fexpe[nexte]) * cimag(mexpf[nftot + nm]);
        double rtmp = 2 * (rt1 + rt2);

        nexte++;
        mexpphys[nptot + ival] += rtmp * _Complex_I;
      }

      for (int nm = 2; nm < numfour[i]; nm += 2) {
        double rt1 = creal(fexpo[nexto]) * creal(mexpf[nftot + nm]);
        double rt2 = cimag(fexpo[nexto]) * cimag(mexpf[nftot + nm]);
        double rtmp = 2 * (rt1 - rt2);

        nexto++;
        mexpphys[nptot + ival] += rtmp;
      }
    }
    nftot += numfour[i];
    nptot += numphys[i]/2;
  }
}

void exponential_to_local(fmm_box_t *tbox) {
  fmm_box_t *uall[36], *u1234[16], *dall[36], *d5678[16], 
    *nall[24], *n1256[8], *n12[4], *n56[4], 
    *sall[24], *s3478[8], *s34[4], *s78[4],
    *eall[16], *e1357[4], *e13[2], *e57[2], *e1[3], *e3[3], *e5[3], *e7[3], 
    *wall[16], *w2468[4], *w24[2], *w68[2], *w2[3], *w4[3], *w6[3], *w8[3]; 
  
  int nuall, xuall[36], yuall[36];
  int nu1234, x1234[16], y1234[16];
  int ndall, xdall[36], ydall[36];
  int nd5678, x5678[16], y5678[16];
  int nnall, xnall[24], ynall[24];
  int nn1256, x1256[8], y1256[8];
  int nn12, x12[4], y12[4];
  int nn56, x56[4], y56[4];
  int nsall, xsall[24], ysall[24];
  int ns3478, x3478[8], y3478[8];
  int ns34, x34[4], y34[4];
  int ns78, x78[4], y78[4];
  int neall, xeall[16], yeall[16];
  int ne1357, x1357[4], y1357[4];
  int ne13, x13[2], y13[2];
  int ne57, x57[2], y57[2];
  int ne1, x1[3], y1[3];
  int ne3, x3[3], y3[3];
  int ne5, x5[3], y5[3];
  int ne7, x7[3], y7[3];
  int nwall, xwall[16], ywall[16];
  int nw2468, x2468[4], y2468[4];
  int nw24, x24[2], y24[2];
  int nw68, x68[2], y68[2];
  int nw2, x2[3], y2[3];
  int nw4, x4[3], y4[3];
  int nw6, x6[3], y6[3];
  int nw8, x8[3], y8[3];

  build_merged_list2(tbox, uall, &nuall, xuall, yuall, 
                     u1234, &nu1234, x1234, y1234,
                     dall, &ndall, xdall, ydall, 
                     d5678, &nd5678, x5678, y5678,
                     nall, &nnall, xnall, ynall, 
                     n1256, &nn1256, x1256, y1256,
                     n12, &nn12, x12, y12,
                     n56, &nn56, x56, y56,
                     sall, &nsall, xsall, ysall, 
                     s3478, &ns3478, x3478, y3478,
                     s34, &ns34, x34, y34,
                     s78, &ns78, x78, y78,
                     eall, &neall, xeall, yeall, 
                     e1357, &ne1357, x1357, y1357,
                     e13, &ne13, x13, y13,
                     e57, &ne57, x57, y57,
                     e1, &ne1, x1, y1,
                     e3, &ne3, x3, y3,
                     e5, &ne5, x5, y5,
                     e7, &ne7, x7, y7,
                     wall, &nwall, xwall, ywall, 
                     w2468, &nw2468, x2468, y2468,
                     w24, &nw24, x24, y24, 
                     w68, &nw68, x68, y68,
                     w2, &nw2, x2, y2,
                     w4, &nw4, x4, y4,
                     w6, &nw6, x6, y6,
                     w8, &nw8, x8, y8); 

  int nexpmax        = fmm_param->nexpmax; 
  int nexptotp       = fmm_param->nexptotp;
  int pgsz           = fmm_param->pgsz;
  double complex *xs = fmm_param->xs;
  double complex *ys = fmm_param->ys;
  double *zs         = fmm_param->zs; 
  double *rdplus     = fmm_param->rdplus;
  double *rdminus    = fmm_param->rdminus; 
  double *scale      = fmm_param->scale; 
  double myscale = scale[tbox->level + 1]; 

  double complex *work = CALLOC(nexpmax * 16, sizeof(double complex));
  double complex *mw1 = CALLOC(pgsz, sizeof(double complex));
  double complex *mw2 = CALLOC(pgsz, sizeof(double complex));
  double complex *mexpf1 = CALLOC(nexpmax, sizeof(double complex));
  double complex *mexpf2 = CALLOC(nexpmax, sizeof(double complex));
  double complex *temp   = CALLOC(nexpmax, sizeof(double complex));

  assert(work != NULL);
  assert(mw1 != NULL);
  assert(mw2 != NULL);
  assert(mexpf1 != NULL);
  assert(mexpf2 != NULL);
  assert(temp != NULL); 

  // Process z direction exponential expansions
  double complex *mexuall  = work;
  double complex *mexu1234 = work + nexpmax;
  double complex *mexdall  = mexu1234 + nexpmax;
  double complex *mexd5678 = mexdall + nexpmax;

  make_ulist(1, uall, nuall, xuall, yuall, mexuall, 0);
  make_ulist(1, u1234, nu1234, x1234, y1234, mexu1234, 0);
  make_dlist(0, dall, ndall, xdall, ydall, mexdall, 0);
  make_dlist(0, d5678, nd5678, x5678, y5678, mexd5678, 0);

  if (tbox->child[0]) {
    fmm_box_t *child = tbox->child[0]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * myscale;
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[1]) {
    fmm_box_t *child = tbox->child[1]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * conj(xs[3 * j]) * myscale;
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      exponential_to_local_p1(temp, mexpf2);
      iexp2++;
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[2]) {
    fmm_box_t *child = tbox->child[2]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[3]) {
    fmm_box_t *child = tbox->child[3]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[4]) {
    fmm_box_t *child = tbox->child[4]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * myscale;
      iexp2++;
    }

    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }  

  if (tbox->child[5]) {
    fmm_box_t *child = tbox->child[5]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * xs[3 * j] * myscale;
      iexp2++;
    }

    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }
  
  if (tbox->child[6]) {
    fmm_box_t *child = tbox->child[6]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[7]) {
    fmm_box_t *child = tbox->child[7]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  // Process y direction exponential expansions
  double complex *mexnall  = work;
  double complex *mexn1256 = mexnall + nexpmax;
  double complex *mexn12   = mexn1256 + nexpmax;
  double complex *mexn56   = mexn12 + nexpmax;
  double complex *mexsall  = mexn56 + nexpmax;
  double complex *mexs3478 = mexsall + nexpmax;
  double complex *mexs34   = mexs3478 + nexpmax;
  double complex *mexs78   = mexs34 + nexpmax;

  make_ulist(3, nall, nnall, xnall, ynall, mexnall, 0);
  make_ulist(3, n1256, nn1256, x1256, y1256, mexn1256, 0);
  make_ulist(3, n12, nn12, x12, y12, mexn12, 0);
  make_ulist(3, n56, nn56, x56, y56, mexn56, 0);
  make_dlist(2, sall, nsall, xsall, ysall, mexsall, 0);
  make_dlist(2, s3478, ns3478, x3478, y3478, mexs3478, 0);
  make_dlist(2, s34, ns34, x34, y34, mexs34, 0);
  make_dlist(2, s78, ns78, x78, y78, mexs78, 0);

  if (tbox->child[0]) {
    fmm_box_t *child = tbox->child[0];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * myscale;
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * myscale;
      iexp1++;
    }

    if (nn12) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn12[j] * zs[3 * j + 1] * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[1]) {
    fmm_box_t *child = tbox->child[1];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (nn12) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn12[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[2]) {
    fmm_box_t *child = tbox->child[2]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * myscale;
      iexp2++;
    }

    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * myscale;
      iexp2++;
    }

    if (ns34) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs34[j] * zs[3 * j + 1] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[3]) {
    fmm_box_t *child = tbox->child[3];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (ns34) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs34[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[4]) {
    fmm_box_t *child = tbox->child[4];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * conj(xs[3 * j]) * myscale;
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp1++;
    }

    if (nn56) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn56[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      exponential_to_local_p1(temp, mexpf2);
      iexp2++;
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[5]) {
    fmm_box_t *child = tbox->child[5];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
    }

    if (nn56) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn56[j] * zs[3 * j + 1] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[6]) {
    fmm_box_t *child = tbox->child[6];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * xs[3 * j] * myscale;
      iexp2++;
    }

    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp2++;
    }

    if (ns78) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs78[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[7]) {
    fmm_box_t *child = tbox->child[7]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * 
          conj(xs[3 * j] * ys[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (ns78) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs78[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  // Process x direction exponential expansions
  double complex *mexeall  = work;
  double complex *mexe1357 = mexeall + nexpmax;
  double complex *mexe13   = mexe1357 + nexpmax;
  double complex *mexe57   = mexe13 + nexpmax;
  double complex *mexe1    = mexe57 + nexpmax;
  double complex *mexe3    = mexe1 + nexpmax;
  double complex *mexe5    = mexe3 + nexpmax;
  double complex *mexe7    = mexe5 + nexpmax;
  double complex *mexwall  = mexe7 + nexpmax;
  double complex *mexw2468 = mexwall + nexpmax;
  double complex *mexw24   = mexw2468 + nexpmax;
  double complex *mexw68   = mexw24 + nexpmax;
  double complex *mexw2    = mexw68 + nexpmax;
  double complex *mexw4    = mexw2 + nexpmax;
  double complex *mexw6    = mexw4 + nexpmax;
  double complex *mexw8    = mexw6 + nexpmax;

  make_ulist(5, eall, neall, xeall, yeall, mexeall, 0);
  make_ulist(5, e1357, ne1357, x1357, y1357, mexe1357, 0);
  make_ulist(5, e13, ne13, x13, y13, mexe13, 0);
  make_ulist(5, e57, ne57, x57, y57, mexe57, 0);
  make_ulist(5, e1, ne1, x1, y1, mexe1, 0);
  make_ulist(5, e3, ne3, x3, y3, mexe3, 0);
  make_ulist(5, e5, ne5, x5, y5, mexe5, 0);
  make_ulist(5, e7, ne7, x7, y7, mexe7, 0);
  make_dlist(4, wall, nwall, xwall, ywall, mexwall, 0);
  make_dlist(4, w2468, nw2468, x2468, y2468, mexw2468, 0);
  make_dlist(4, w24, nw24, x24, y24, mexw24, 0);
  make_dlist(4, w68, nw68, x68, y68, mexw68, 0);
  make_dlist(4, w2, nw2, x2, y2, mexw2, 0);
  make_dlist(4, w4, nw4, x4, y4, mexw4, 0);
  make_dlist(4, w6, nw6, x6, y6, mexw6, 0);
  make_dlist(4, w8, nw8, x8, y8, mexw8, 0);

  if (tbox->child[0]) {
    fmm_box_t *child = tbox->child[0];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * myscale;
      iexp1++;
    }

    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * myscale;
      iexp1++;
    }

    if (ne13) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe13[j] * zs[3 * j + 1] * myscale;
      iexp1++;
    }

    if (ne1) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1[j] * zs[3 * j + 1] * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[1]) {
    fmm_box_t *child = tbox->child[1]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * myscale;
      iexp2++;
    }

    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * myscale;
      iexp2++;
    }

    if (nw24) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw24[j] * zs[3 * j + 1] * myscale;
      iexp2++;
    }

    if (nw2) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2[j] * zs[3 * j + 1] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[2]) {
    fmm_box_t *child = tbox->child[2];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (ne13) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe13[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (ne3) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe3[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[3]) {
    fmm_box_t *child = tbox->child[3];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * conj(ys[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (nw24) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw24[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (nw4) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw4[j] * zs[3 * j + 1] * ys[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[4]) {
    fmm_box_t *child = tbox->child[4];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * xs[3 * j] * myscale;
      iexp1++;
    }

    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp1++;
    }

    if (ne57) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe57[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp1++;
    }

    if (ne5) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe5[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[5]) {
    fmm_box_t *child = tbox->child[5];
    double complex *local = &child->expansion[0];
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * xs[3 * j] * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * conj(xs[3 * j]) * myscale;
      iexp2++;
    }

    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp2++;
    }

    if (nw68) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw68[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp2++;
    }

    if (nw6) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw6[j] * zs[3 * j + 1] * conj(xs[3 * j]) * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[6]) {
    fmm_box_t *child = tbox->child[6];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * xs[3 * j] * 
          conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * xs[3 * j] * 
          conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (ne57) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe57[j] * zs[3 * j + 1] *  xs[3 * j] * 
          conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (ne7) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe7[j] * zs[3 * j + 1] * xs[3 * j] * 
          conj(ys[3 * j]) * myscale;
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1);

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * 
          conj(xs[3 * j]) * ys[3 * j] * myscale;
      iexp2++;
      exponential_to_local_p1(temp, mexpf2);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[7]) {
    fmm_box_t *child = tbox->child[7]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * xs[3 * j] * 
          conj(ys[3 * j]) * myscale;
      iexp1++;
      exponential_to_local_p1(temp, mexpf1);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * conj(xs[3 * j]) * 
          ys[3 * j] * myscale;
      iexp2++;
    }

    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * conj(xs[3 * j]) * 
          ys[3 * j] * myscale;
      iexp2++;
    }

    if (nw68) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw68[j] * zs[3 * j + 1] * conj(xs[3 * j]) * 
          ys[3 * j] * myscale;
      iexp2++;
    }

    if (nw8) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw8[j] * zs[3 * j + 1] * conj(xs[3 * j]) * 
          ys[3 * j] * myscale;
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  FREE(work);
  FREE(mw1);
  FREE(mw2);
  FREE(mexpf1);
  FREE(mexpf2);
  FREE(temp);
}

void exponential_to_local_p1(const double complex *mexpphys,
                             double complex *mexpf) {
  int nlambs = fmm_param->nlambs;
  int *numfour = fmm_param->numfour;
  int *numphys = fmm_param->numphys;
  double complex *fexpback = fmm_param->fexpback;

  int nftot = 0;
  int nptot = 0;
  int next  = 0;

  for (int i = 0; i < nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha / 2;
    mexpf[nftot] = 0;
    for (int ival = 0; ival < nalpha2; ival++) {
      mexpf[nftot] += 2.0 * creal(mexpphys[nptot + ival]);
    }
    mexpf[nftot] /= nalpha;

    for (int nm = 2; nm < numfour[i]; nm += 2) {
      mexpf[nftot + nm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
        double rtmp = 2 * creal(mexpphys[nptot + ival]);
        mexpf[nftot + nm] += fexpback[next] * rtmp;
        next++;
      }
      mexpf[nftot + nm] /= nalpha;
    }

    for (int nm = 1; nm < numfour[i]; nm += 2) {
      mexpf[nftot + nm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
        double complex ztmp = 2 * cimag(mexpphys[nptot + ival]) * _Complex_I;
        mexpf[nftot + nm] += fexpback[next] * ztmp;
        next++;
      }
      mexpf[nftot + nm] /= nalpha;
    }
    nftot += numfour[i];
    nptot += numphys[i] / 2;
  }
}

void exponential_to_local_p2(int iexpu, const double complex *mexpu,
                             int iexpd, const double complex *mexpd,
                             double complex *local) {

  int pterms = fmm_param->pterms;
  int nlambs = fmm_param->nlambs;
  int nexptot = fmm_param->nexptot;
  int pgsz = fmm_param->pgsz;
  double *whts = fmm_param->whts;
  double *rlams = fmm_param->rlams;
  int *numfour = fmm_param->numfour;
  double *ytopcs = fmm_param->ytopcs;

  double *rlampow = CALLOC(pterms + 1, sizeof(double));
  double complex *zeye = CALLOC(pterms + 1, sizeof(double complex));
  double complex *mexpplus = CALLOC(nexptot, sizeof(double complex));
  double complex *mexpminus = CALLOC(nexptot, sizeof(double complex));

  assert(rlampow != NULL);
  assert(zeye != NULL);
  assert(mexpplus != NULL);
  assert(mexpminus != NULL);

  zeye[0] = 1.0;
  for (int i = 1; i <= pterms; i++)
    zeye[i] = zeye[i - 1] * _Complex_I;

  for (int i = 0; i < pgsz; i++)
    local[i] = 0;

  for (int i = 0; i < nexptot; i++) {
    if (iexpu <= 0) {
      mexpplus[i] = mexpd[i];
      mexpminus[i] = mexpd[i];
    } else if (iexpd <= 0) {
      mexpplus[i] = mexpu[i];
      mexpminus[i] = -mexpu[i];
    } else {
      mexpplus[i] = mexpd[i] + mexpu[i];
      mexpminus[i] = mexpd[i] - mexpu[i];
    }
  }

  int ntot = 0;
  for (int nell = 0; nell < nlambs; nell++) {
    rlampow[0] = whts[nell];
    double rmul = rlams[nell];
    for (int j = 1; j <= pterms; j++)
      rlampow[j] = rlampow[j - 1] * rmul;

    int mmax = numfour[nell]-1;
    for (int mth = 0; mth <= mmax; mth += 2) {
      int offset = mth * (pterms + 1);
      for (int nm = mth; nm <= pterms; nm += 2) {
        int index = offset + nm;
        int ncurrent = ntot + mth;
        rmul = rlampow[nm];
        local[index] += rmul * mexpplus[ncurrent];
      }

      for (int nm = mth + 1; nm <= pterms; nm += 2) {
        int index = offset + nm;
        int ncurrent = ntot + mth;
        rmul = rlampow[nm];
        local[index] += rmul * mexpminus[ncurrent];
      }
    }

    for (int mth = 1; mth <= mmax; mth += 2) {
      int offset = mth * (pterms + 1);
      for (int nm = mth + 1; nm <= pterms; nm += 2) {
        int index = nm + offset;
        int ncurrent = ntot+mth;
        rmul = rlampow[nm];
        local[index] += rmul * mexpplus[ncurrent];
      }

      for (int nm = mth; nm <= pterms; nm += 2) {
        int index = nm + offset;
        int ncurrent = ntot + mth;
        rmul = rlampow[nm];
        local[index] += rmul * mexpminus[ncurrent];
      }
    }
    ntot += numfour[nell];
  }

  for (int mth = 0; mth <= pterms; mth++) {
    int offset1 = mth * (pterms + 1);
    int offset2 = mth * (pterms + 2); 
    for (int nm = mth; nm <= pterms; nm++) {
      int index1 = nm + offset1;
      int index2 = nm + offset2; 
      local[index1] *= zeye[mth] * ytopcs[index2];
    }
  }

  FREE(rlampow);
  FREE(zeye);
  FREE(mexpplus);
  FREE(mexpminus);
}

void local_to_local(fmm_box_t *tbox) {
  const double complex var[5] =
    {1, 1 - _Complex_I, -1 - _Complex_I, -1 + _Complex_I, 1 + _Complex_I};
  const double arg = sqrt(2) / 2.0;
  const int ifld[8] = {1, 2, 4, 3, 1, 2, 4, 3};

  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz;
  double *dc = fmm_param->dc;

  double complex *localn = CALLOC(pgsz, sizeof(double complex));
  double complex *marray = CALLOC(pgsz, sizeof(double complex));
  double complex *ephi = CALLOC(1 + pterms, sizeof(double complex));
  double *powers = CALLOC(1 + pterms, sizeof(double));

  assert(localn != NULL);
  assert(marray != NULL);
  assert(ephi != NULL);
  assert(powers != NULL);

  double complex *plocal = &tbox->expansion[0]; 
  for (int i = 0; i < 8; i++) {
    fmm_box_t *child = tbox->child[i];
    if (child) {
      double *rd = (i < 4 ? fmm_param->rdsq3 : fmm_param->rdmsq3);
      int ifl = ifld[i];
      ephi[0] = 1.0;
      ephi[1] = arg * var[ifl];
      double dd = -sqrt(3) / 4.0;
      powers[0] = 1.0;

      for (int ell = 1; ell <= pterms; ell++)
        powers[ell] = powers[ell - 1] * dd;

      for (int ell = 2; ell <= pterms; ell++)
        ephi[ell] = ephi[ell - 1] * ephi[1];
      
      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          localn[index] = conj(ephi[m]) * plocal[index];
        }
      }
      
      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        int offset1 = (pterms + m) * pgsz;
        int offset2 = (pterms - m) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = localn[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp++) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
        }
      }

      for (int k = 0; k <= pterms; k++) {
        int offset = k * (pterms + 1);
        for (int j = k; j <= pterms; j++) {
          int index = j + offset;
          localn[index] = marray[index];
          for (int ell = 1; ell <= pterms - j; ell++) {
            int index1 = ell + index;
            int index2 = ell + j + k + ell * (2 * pterms + 1);
            int index3 = ell + j - k + ell * (2 * pterms + 1);
            localn[index] += marray[index1] * powers[ell] * 
              dc[index2] * dc[index3];
          }
        }
      }
      
      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        int offset1 = (pterms + m) * pgsz;
        int offset2 = (pterms - m) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = localn[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] -= localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
	    
          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
        }
      }
      
      for (int m = 1; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (pterms + m) * pgsz;
        int offset2 = (pterms - m) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = -localn[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }

          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] -= localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
        }
      }
      
      powers[0] = 1.0;
      for (int ell = 1; ell <= pterms; ell++)
        powers[ell] = powers[ell - 1] / 2;

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          int index = offset + ell;
          localn[index] = ephi[m] * marray[index] * powers[ell];
        }
      }
      
      double complex *clocal = &child->expansion[0]; 
      for (int m = 0; m < pgsz; m++)
        clocal[m] += localn[m];
    }
  }

  FREE(ephi);
  FREE(powers);
  FREE(localn);
  FREE(marray);
}

void local_to_target(fmm_box_t *tbox) {
  int pgsz = fmm_param->pgsz;
  int pterms = fmm_param->pterms;
  double *ytopc = fmm_param->ytopc;
  double *ytopcs = fmm_param->ytopcs;
  double *ytopcsinv = fmm_param->ytopcsinv;
  double *scale = fmm_param->scale; 

  double *p = CALLOC(pgsz, sizeof(double));
  double *powers = CALLOC(pterms + 1, sizeof(double));
  double complex *ephi = CALLOC(pterms + 1, sizeof(double complex));

  assert(p != NULL);
  assert(powers != NULL);
  assert(ephi != NULL);

  double complex *local = &tbox->expansion[0]; 
  double myscale = scale[tbox->level]; 
  double h = fmm_dag->size / (1 << (tbox->level + 1)); 
  double center[3]; 
  center[0] = fmm_dag->corner[0] + (2 * tbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * tbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * tbox->idz + 1) * h;

  int ntargets = tbox->npts; 
  int addr = tbox->addr; 
  const double precision = 1e-14;

  for (int i = 0; i < ntargets; i++) {
    int ptr = addr + i; 
    double *pot = &potential_sorted[ptr]; 
    double *field = &field_sorted[3 * ptr]; 
    double field0, field1, field2, rloc, cp, rpotz = 0.0;
    double complex cpz, zs1 = 0.0, zs2 = 0.0, zs3 = 0.0;
    double rx = targets_sorted[3 * ptr] - center[0];
    double ry = targets_sorted[3 * ptr + 1] - center[1]; 
    double rz = targets_sorted[3 * ptr + 2] - center[2];
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 0.0 : rz / d);
    ephi[0] = (proj <= precision * d ? 1.0 : rx / proj + _Complex_I * ry / proj);
    d *= myscale;
    double dd = d;

    powers[0] = 1.0;
    for (int ell = 1; ell <= pterms; ell++) {
      powers[ell] = dd;
      dd *= d;
      ephi[ell] = ephi[ell - 1] * ephi[0];
    }

    lgndr(pterms, ctheta, p);
    *pot += creal(local[0]);

    field2 = 0.0;
    for (int ell = 1; ell <= pterms; ell++) {
      rloc = creal(local[ell]);
      cp = rloc * powers[ell] * p[ell];
      *pot += cp;
      cp = powers[ell - 1] * p[ell - 1] * ytopcs[ell - 1];
      cpz = local[ell + pterms + 1] * cp * ytopcsinv[ell + pterms + 2];
      zs2 = zs2 + cpz;
      cp = rloc * cp * ytopcsinv[ell];
      field2 += cp;
    }
    
    for (int ell = 1; ell <= pterms; ell++) {
      for (int m = 1; m <= ell; m++) {
        int index = ell + m * (pterms + 1);
        cpz = local[index] * ephi[m - 1];
        rpotz += creal(cpz) * powers[ell] * ytopc[ell + m * (pterms + 2)] * 
          p[index];	
      }

      for (int m = 1; m <= ell - 1; m++) {
        int index1 = ell + m * (pterms + 1);
        int index2 = index1 - 1;
        zs3 += local[index1] * ephi[m - 1] * powers[ell - 1] * 
          ytopc[ell - 1 + m * (pterms + 2)] * p[index2] * 
          ytopcs[ell - 1 + m * (pterms + 2)] * 
          ytopcsinv[ell + m * (pterms + 2)];
      }

      for (int m = 2; m <= ell; m++) {
        int index1 = ell + m * (pterms + 1);
        int index2 = ell - 1 + (m - 1) * (pterms + 1);
        zs2 += local[index1] * ephi[m - 2] * 
          ytopcs[ell - 1 + (m - 1) * (pterms + 2)] * 
          ytopcsinv[ell + m * (pterms + 2)] * powers[ell - 1] * 
          ytopc[ell - 1 + (m - 1) * (pterms + 2)] * p[index2];
      }
      
      for (int m = 0; m <= ell - 2; m++) {
        int index1 = ell + m * (pterms + 1);
        int index2 = ell - 1 + (m + 1) * (pterms + 1);
        zs1 += local[index1] * ephi[m] * 
          ytopcs[ell - 1 + (m + 1) * (pterms + 2)] * 
          ytopcsinv[ell + m * (pterms + 2)] * powers[ell - 1] * 
          ytopc[ell - 1 + (m + 1) * (pterms + 2)] * p[index2];        
      }
    }

    *pot += 2.0 * rpotz;
    field0 = creal(zs2 - zs1);
    field1 = -cimag(zs2 + zs1);
    field2 += 2.0 * creal(zs3);
    
    field[0] += field0 * myscale;
    field[1] += field1 * myscale;
    field[2] -= field2 * myscale;
  }
  FREE(powers);
  FREE(ephi);
  FREE(p);
}

void source_to_local(fmm_box_t *tbox, fmm_box_t *sbox) {
  int pgsz   = fmm_param->pgsz; 
  int pterms = fmm_param->pterms; 
  double *scale = fmm_param->scale; 
  double *ytopc = fmm_param->ytopc; 

  int nsources = sbox->npts; 
  int addr = sbox->addr; 

  double h = fmm_dag->size / (1 << (tbox->level + 1)); 
  double center[3]; 
  center[0] = fmm_dag->corner[0] + (2 * tbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * tbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * tbox->idz + 1) * h;
  double complex *local = &tbox->expansion[0]; 
  double myscale = scale[tbox->level]; 

  double *powers = CALLOC(pterms + 3, sizeof(double));
  double *p = CALLOC(pgsz, sizeof(double));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex));

  assert(powers != NULL);
  assert(p != NULL);
  assert(ephi != NULL); 

  const double precision = 1e-14; 

  for (int i = 0; i < nsources; i++) {
    int ptr = addr + i; 
    double rx = sources_sorted[3 * ptr]     - center[0];
    double ry = sources_sorted[3 * ptr + 1] - center[1]; 
    double rz = sources_sorted[3 * ptr + 2] - center[2]; 
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 1.0 : rz / d);
    ephi[0] = (proj <= precision * d ? 1.0 : rx / proj - _Complex_I * ry / proj);

    d = 1.0/d;
    powers[0] = 1.0;
    powers[1] = d;
    d /= myscale; 

    for (int ell = 2; ell <= pterms + 2; ell++) 
      powers[ell] = powers[ell - 1] * d;

    for (int ell = 1; ell <= pterms + 1; ell++) 
      ephi[ell] = ephi[ell - 1] * ephi[0]; 

    local[0] += charges_sorted[ptr] * powers[1]; 
    
    lgndr(pterms, ctheta, p); 

    for (int ell = 1; ell <= pterms; ell++) 
      local[ell] += charges_sorted[ptr] * p[ell] * powers[ell + 1]; 

    for (int m = 1; m <= pterms; m++) {
      int offset1 = m * (pterms + 1);
      int offset2 = m * (pterms + 2); 
      for (int ell = m; ell <= pterms; ell++) {
        int index1 = offset1 + ell; 
        int index2 = offset2 + ell;
        local[index1] += charges_sorted[ptr] * powers[ell + 1] * 
          ytopc[index2] * p[index1] * ephi[m - 1]; 
      }
    }
  }

  FREE(powers);
  FREE(p);
  FREE(ephi); 
}


void multipole_to_target(fmm_box_t *tbox, fmm_box_t *sbox) {
  int pgsz      = fmm_param->pgsz; 
  int pterms    = fmm_param->pterms; 
  double *scale = fmm_param->scale; 
  double *ytopc = fmm_param->ytopc;
  double *ytopcs = fmm_param->ytopcs; 
  double *ytopcsinv = fmm_param->ytopcsinv; 

  double myscale = scale[sbox->level]; 
  double complex *multipole = &sbox->expansion[0]; 
  double h = fmm_dag->size / (1 << (sbox->level + 1));
  double center[3]; 
  center[0] = fmm_dag->corner[0] + (2 * sbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * sbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * sbox->idz + 1) * h;

  int ntargets = tbox->npts; 
  int addr = tbox->addr; 

  double *p = CALLOC((pterms + 2) * (pterms + 2), sizeof(double));
  double *powers = CALLOC(pterms + 4, sizeof(double));
  double complex *ephi = CALLOC(pterms + 3, sizeof(double complex)); 

  assert(p != NULL);
  assert(powers != NULL);
  assert(ephi != NULL);

  const double precis = 1e-14; 
  
  for (int i = 0; i < ntargets; i++) {
    int ptr = addr + i; 
    double *pot = &potential_sorted[ptr]; 
    double *field = &field_sorted[3 * ptr]; 
    double rpotz = 0.0, field1 = 0.0, field2 = 0.0, field3 = 0.0;
    double complex zs1 = 0.0, zs2 = 0.0, zs3 = 0.0;
    double rx = targets_sorted[3 * ptr]     - center[0];
    double ry = targets_sorted[3 * ptr + 1] - center[1]; 
    double rz = targets_sorted[3 * ptr + 2] - center[2]; 
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj); 
    double d = sqrt(rr); 
    double ctheta = (d <= precis ? 0.0 : rz / d); 
    ephi[0] = (proj <= precis * d ? 1.0 : rx / proj + _Complex_I * ry / proj); 
    d = 1.0 / d; 
    powers[0] = 1.0;
    powers[1] = d; 
    d /= myscale; 

    for (int ell = 1; ell <= pterms + 2; ell++) {
      powers[ell + 1] = powers[ell] * d; 
      ephi[ell] = ephi[ell - 1] * ephi[0]; 
    }

    lgndr(pterms + 1, ctheta, p); 
    
    double rtemp; 

    double rmp = creal(multipole[0]); 
    *pot += rmp * powers[1]; 
    double complex cpz = ephi[0] * rmp * powers[2] * ytopc[1 + pterms + 2] * 
      p[1 + pterms + 2] * ytopcsinv[1 + pterms + 2]; 
    zs1 += cpz; 
    double cp = rmp * powers[2] * p[1] * ytopcs[0] * ytopcsinv[1]; 
    field3 = cp;

    for (int ell = 1; ell <= pterms; ell++) {
      rmp = creal(multipole[ell]); 
      cp = rmp * powers[ell + 1] * p[ell]; 
      *pot += cp; 
      zs1 += ephi[0] * rmp * powers[ell + 2] * ytopcs[ell + 1+ pterms + 2] * 
        p[ell + 1 + pterms+2] * ytopcs[ell] * ytopcsinv[ell + 1 + pterms + 2];
      cpz = multipole[ell + pterms + 1]; 
      rtemp = powers[ell + 2] * p[ell + 1] * ytopcsinv[ell + 1]; 
      zs2 += cpz * rtemp * ytopcs[ell + pterms + 2]; 
      cp = rmp * rtemp * ytopcs[ell]; 
      field3 += cp;
    }

    for (int m = 1; m <= pterms; m++) {
      int offset1 = m * (pterms + 1); 
      int offset2 = m * (pterms + 2);
      int offset5 = (m + 1) * (pterms + 2);
      int offset6 = (m - 1) * (pterms + 2);
      for (int ell = m; ell <= pterms; ell++) {
        int index1 = ell + offset1; 
        int index2 = ell + offset2;
        int index5 = ell + 1 + offset5; 
        int index6 = ell + 1 + offset6;
        cpz = multipole[index1] * powers[ell + 1] * ytopc[index2] * p[index2]; 
        rpotz += creal(cpz * ephi[m - 1]); 
        cpz = multipole[index1] * ytopcs[index2] * powers[ell + 2]; 
        zs1 += cpz * ephi[m] * ytopcsinv[index5] * ytopc[index5] * p[index5];
        if (m > 1) 
          zs2 += cpz * ephi[m - 2] * ytopcsinv[index6] * 
            ytopc[index6] * p[index6]; 
        zs3 += cpz * ephi[m - 1] * ytopc[index2 + 1] * p[index2 + 1] * 
          ytopcsinv[index2 + 1];
      }
    }

    *pot += 2.0 * rpotz;
    field1 = creal(zs2 - zs1);
    field2 = -cimag(zs2 + zs1);
    field3 = field3 + 2.0 * creal(zs3); 
    
    field[0] += field1 * myscale;
    field[1] += field2 * myscale;
    field[2] += field3 * myscale; 
  }
  
  FREE(p);
  FREE(powers);
  FREE(ephi);
}

void direct_evaluation(const fmm_box_t *tbox, const fmm_box_t *sbox) {
  int start1 = tbox->addr, num1 = tbox->npts, end1 = start1 + num1 - 1;
  int start2 = sbox->addr, num2 = sbox->npts, end2 = start2 + num2 - 1;

  for (int i = start1; i <= end1; i++) {
    int i3 = i * 3; 
    double pot = 0, fx = 0, fy = 0, fz = 0; 

#ifdef icc
    #pragma simd reduction(+:pot)
    #pragma simd reduction(+:fx)
    #pragma simd reduction(+:fy)
    #pragma simd reduction(+:fz)
#endif
    for (int j = start2; j <= end2; j++) {
      int j3 = j * 3; 
      double rx = targets_sorted[i3]     - sources_sorted[j3];
      double ry = targets_sorted[i3 + 1] - sources_sorted[j3 + 1];
      double rz = targets_sorted[i3 + 2] - sources_sorted[j3 + 2]; 
      double q = charges_sorted[j]; 
      double rr = rx * rx + ry * ry + rz * rz;
      double rdis = sqrt(rr); 

      if (rr) {
        pot += q / rdis;
        double rmul = q / (rdis * rr);
        fx += rmul * rx;
        fy += rmul * ry;
        fz += rmul * rz;
      }
    }

    potential_sorted[i]   += pot;
    field_sorted[i * 3]     += fx;
    field_sorted[i * 3 + 1] += fy;
    field_sorted[i * 3 + 2] += fz;
  }
}

#elif YUKAWA

void source_to_multipole(fmm_box_t *sbox) {
  int pgsz     = fmm_param->pgsz; 
  int pterms   = fmm_param->pterms; 
  double beta  = fmm_param->beta; 
  double *ytop = fmm_param->ytop; 
  int addr     = sbox->addr; 
  int level    = sbox->level; 

  double complex *multipole = &sbox->expansion[0]; 
  Real_t *sources = &sources_sorted[3 * addr]; 
  Real_t *charges = &charges_sorted[addr]; 
  int nsources = sbox->npts; 
  double myscale = fmm_param->sfactor[level]; 
  double h = fmm_dag->size / (1 << (level + 1)); 
  double center[3]; 
  center[0] = fmm_dag->corner[0] + (2 * sbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * sbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * sbox->idz + 1) * h;

  double *p = CALLOC(pgsz, sizeof(double));
  double *bi = CALLOC(pterms + 2, sizeof(double));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex)); 

  assert(p != NULL);
  assert(bi != NULL);
  assert(ephi != NULL); 

  const double precision=1.0e-14;

  for (int i = 0; i < nsources; i++) {
    int i3 = i * 3; 
    double rx = sources[i3]     - center[0];
    double ry = sources[i3 + 1] - center[1];
    double rz = sources[i3 + 2] - center[2];
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = (d <= precision ? 1.0 : rz / d);
    ephi[0] = (proj <= precision * d ? 1.0 : rx / proj - _Complex_I * ry / proj);
    for (int ell = 1; ell < pterms; ell++) {
      ephi[ell] = ephi[ell-1] * ephi[0];
    }

    double rk = d * beta;
    int ncalc;
    in(myscale, rk, pterms, bi, &ncalc);
    multipole[0] += charges[i] * bi[0];
    lgndr(pterms, ctheta, p);
    for (int ell = 1; ell <= pterms; ell++) {
      double complex cp = charges[i] * p[ell] * bi[ell] * ytop[ell];
      multipole[ell] += cp;
    }

    for (int m = 1; m <= pterms; m++) {
      int offset = m * (pterms + 1);
      for (int ell = m; ell <= pterms; ell++) {
        double cp = charges[i] * bi[ell] * ytop[ell + offset] * p[ell + offset];
        multipole[ell + offset] += cp * ephi[m - 1];
      }
    }
  }
  FREE(p);
  FREE(bi);
  FREE(ephi);
}

void multipole_to_multipole(fmm_box_t *sbox) {
  const double complex var[5] = 
    {0, -1 + _Complex_I, 1 + _Complex_I, 1 - _Complex_I, -1 - _Complex_I};
  const double arg = sqrt(2.0) / 2.0;
  const int yifl[8] = {3, 4, 2, 1, 3, 4, 2, 1}; 

  int level  = sbox->level; 
  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz;
  int dcpgsz = fmm_param->dcpgsz; 
  double *dc = &fmm_param->dcu[dcpgsz * level]; 
  double complex *pmpole = &sbox->expansion[0]; 

  double complex *mpolen = CALLOC(pgsz, sizeof(double complex));
  double complex *marray = CALLOC(pgsz, sizeof(double complex));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex)); 
  assert(mpolen != NULL);
  assert(marray != NULL);
  assert(ephi != NULL); 

  for (int i = 0; i < 8; i++) {
    fmm_box_t *child = sbox->child[i]; 
    if (child != NULL) {
      int ifl = yifl[i];
      double *rd = (i < 4 ? fmm_param->rdsq3 : fmm_param->rdmsq3);
      double complex *cmpole = &child->expansion[0]; 

      ephi[0] = 1.0;
      ephi[1] = arg * var[ifl];
      for (int ell = 1; ell <= pterms; ell++)
        ephi[ell + 1] = ephi[ell] * ephi[1];

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          mpolen[index] = conj(ephi[m]) * cmpole[index];
        }
      }

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = offset + ell;
          marray[index] = mpolen[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp++) {
            int index1 = mp * (pterms + 1) + ell;
            marray[index] += mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }
        }
      }

      for (int k = 0; k<= pterms; k++) {
        int offset = k * (pterms + 1);
        for (int j = k; j <= pterms; j++) {
          int index = offset + j;
          int offset1 = j * (pterms + 1) + k;
          mpolen[index] = 0;
          for (int ell = k; ell <= pterms; ell++)
            mpolen[index] += marray[ell + offset] * dc[offset1 + ell * pgsz];
        }
      }

      for (int m = 0; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = mpolen[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] -= mpolen[index1] * rd[offset1 + index1] +
              conj(mpolen[index1]) * rd[offset2 + index1];
          }
          
          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += mpolen[index1] * rd[offset1 + index1] +
              conj(mpolen[index1]) * rd[offset2 + index1];
          }
        }
      }

      for (int m = 1; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = -mpolen[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] += mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }

          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = ell + mp * (pterms + 1);
            marray[index] -= mpolen[index1] * rd[index1 + offset1] +
              conj(mpolen[index1]) * rd[index1 + offset2];
          }
        }
      }

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          mpolen[index] = ephi[m] * marray[index];
        }
      }
      
      for (int m = 0; m < pgsz; m++)
        pmpole[m] += mpolen[m];
    }
  }
  FREE(mpolen);
  FREE(marray);
  FREE(ephi);
}

void multipole_to_exponential(fmm_box_t *sbox) {
  int level = sbox->level; 

  int pgsz        = fmm_param->pgsz; 
  int nexpmax     = fmm_param->nexpmax; 
  double *rdminus = fmm_param->rdminus;
  double *rdplus  = fmm_param->rdplus; 

  double complex *mw = CALLOC(pgsz, sizeof(double complex));
  double complex *mexpf1 = CALLOC(nexpmax, sizeof(double complex)); 
  double complex *mexpf2 = CALLOC(nexpmax, sizeof(double complex)); 
  assert(mw != NULL);
  assert(mexpf1 != NULL);
  assert(mexpf2 != NULL); 

  multipole_to_exponential_p1(&sbox->expansion[0], mexpf1, mexpf2, level); 
  multipole_to_exponential_p2(mexpf1, &sbox->expansion[pgsz], level); 
  multipole_to_exponential_p2(mexpf2, &sbox->expansion[pgsz + nexpmax], level); 

  rotz2y(&sbox->expansion[0], rdminus, mw); 
  multipole_to_exponential_p1(mw, mexpf1, mexpf2, level);
  multipole_to_exponential_p2(mexpf1, &sbox->expansion[pgsz + nexpmax * 2], 
                              level); 
  multipole_to_exponential_p2(mexpf2, &sbox->expansion[pgsz + nexpmax * 3],
                              level);

  rotz2x(&sbox->expansion[0], rdplus, mw); 
  multipole_to_exponential_p1(mw, mexpf1, mexpf2, level); 
  multipole_to_exponential_p2(mexpf1, &sbox->expansion[pgsz + nexpmax * 4], 
                              level); 
  multipole_to_exponential_p2(mexpf2, &sbox->expansion[pgsz + nexpmax * 5],
                              level);
  FREE(mw);
  FREE(mexpf1);
  FREE(mexpf2); 
}

void multipole_to_exponential_p1(const double complex *multipole, 
                                 double complex *mexpu, 
                                 double complex *mexpd, int level) {
  int nlambs   = fmm_param->nlambs;
  int pterms   = fmm_param->pterms; 
  int pgsz     = fmm_param->pgsz; 
  int *numfour = fmm_param->numfour; 
  double *rlsc = &fmm_param->rlsc[pgsz * nlambs * level];

  int ntot = 0;
  for (int nell = 0; nell < nlambs; nell++) {
    double sgn = -1;
    for (int mth = 0; mth <= numfour[nell] - 1; mth++) {
      int ncurrent = ntot + mth;
      double complex ztmp1 = 0;
      double complex ztmp2 = 0;
      sgn = -sgn;
      int offset = mth * (pterms + 1);
      int offset1 = offset + nell * pgsz;
      for (int nm = mth; nm <= pterms; nm += 2)
        ztmp1 += multipole[nm + offset] * rlsc[nm + offset1];
      for (int nm = mth + 1; nm <= pterms; nm += 2)
        ztmp2 += multipole[nm + offset] * rlsc[nm + offset1];
      mexpu[ncurrent] = ztmp1 + ztmp2;
      mexpd[ncurrent] = sgn * (ztmp1 - ztmp2);
    }
    ntot += numfour[nell];
  }
}


void multipole_to_exponential_p2(const double complex *mexpf, 
                                 double complex *mexpphys, int level) {
  int nlambs = fmm_param->nlambs;
  int *numfour = fmm_param->numfour; 
  int *numphys = fmm_param->numphys; 
  double complex *fexpe = &fmm_param->fexpe[15000 * level]; 
  double complex *fexpo = &fmm_param->fexpo[15000 * level]; 

  int nftot = 0;
  int nptot = 0;
  int nexte = 0;
  int nexto = 0;

  for (int i = 0; i < nlambs; i++) {
    for (int ival = 0; ival < numphys[i]/2; ival++) {
      mexpphys[nptot + ival] = mexpf[nftot];
      double sgn = -2;
      for (int nm = 1; nm < numfour[i]; nm += 2) {
        sgn = -sgn;
        double rtmp = sgn * (fexpe[nexte] * mexpf[nftot + nm]);
        nexte++; 
        mexpphys[nptot + ival] += rtmp * _Complex_I;
      }
      
      sgn = 2;
      for (int nm = 2; nm < numfour[i]; nm += 2) {
        sgn = -sgn;
        double rtmp = sgn * (fexpo[nexto] * mexpf[nftot + nm]);
        nexto++;
        mexpphys[nptot + ival] += rtmp;
      }
    }    
    nftot += numfour[i];
    nptot += numphys[i]/2;
  }
}

void exponential_to_local(fmm_box_t *tbox) {
  fmm_box_t *uall[36], *u1234[16], *dall[36], *d5678[16], 
    *nall[24], *n1256[8], *n12[4], *n56[4], 
    *sall[24], *s3478[8], *s34[4], *s78[4],
    *eall[16], *e1357[4], *e13[2], *e57[2], *e1[3], *e3[3], *e5[3], *e7[3], 
    *wall[16], *w2468[4], *w24[2], *w68[2], *w2[3], *w4[3], *w6[3], *w8[3]; 
  
  int nuall, xuall[36], yuall[36];
  int nu1234, x1234[16], y1234[16];
  int ndall, xdall[36], ydall[36];
  int nd5678, x5678[16], y5678[16];
  int nnall, xnall[24], ynall[24];
  int nn1256, x1256[8], y1256[8];
  int nn12, x12[4], y12[4];
  int nn56, x56[4], y56[4];
  int nsall, xsall[24], ysall[24];
  int ns3478, x3478[8], y3478[8];
  int ns34, x34[4], y34[4];
  int ns78, x78[4], y78[4];
  int neall, xeall[16], yeall[16];
  int ne1357, x1357[4], y1357[4];
  int ne13, x13[2], y13[2];
  int ne57, x57[2], y57[2];
  int ne1, x1[3], y1[3];
  int ne3, x3[3], y3[3];
  int ne5, x5[3], y5[3];
  int ne7, x7[3], y7[3];
  int nwall, xwall[16], ywall[16];
  int nw2468, x2468[4], y2468[4];
  int nw24, x24[2], y24[2];
  int nw68, x68[2], y68[2];
  int nw2, x2[3], y2[3];
  int nw4, x4[3], y4[3];
  int nw6, x6[3], y6[3];
  int nw8, x8[3], y8[3];

  build_merged_list2(tbox, uall, &nuall, xuall, yuall, 
                     u1234, &nu1234, x1234, y1234,
                     dall, &ndall, xdall, ydall, 
                     d5678, &nd5678, x5678, y5678,
                     nall, &nnall, xnall, ynall, 
                     n1256, &nn1256, x1256, y1256,
                     n12, &nn12, x12, y12,
                     n56, &nn56, x56, y56,
                     sall, &nsall, xsall, ysall, 
                     s3478, &ns3478, x3478, y3478,
                     s34, &ns34, x34, y34,
                     s78, &ns78, x78, y78,
                     eall, &neall, xeall, yeall, 
                     e1357, &ne1357, x1357, y1357,
                     e13, &ne13, x13, y13,
                     e57, &ne57, x57, y57,
                     e1, &ne1, x1, y1,
                     e3, &ne3, x3, y3,
                     e5, &ne5, x5, y5,
                     e7, &ne7, x7, y7,
                     wall, &nwall, xwall, ywall, 
                     w2468, &nw2468, x2468, y2468,
                     w24, &nw24, x24, y24, 
                     w68, &nw68, x68, y68,
                     w2, &nw2, x2, y2,
                     w4, &nw4, x4, y4,
                     w6, &nw6, x6, y6,
                     w8, &nw8, x8, y8); 
  
  int level = tbox->level + 1;
  int nexpmax        = fmm_param->nexpmax; 
  int nexptotp       = fmm_param->vnexptotp[level]; 
  int pgsz           = fmm_param->pgsz; 
  double complex *xs = &fmm_param->xs[3 * nexpmax * level];
  double complex *ys = &fmm_param->ys[3 * nexpmax * level]; 
  double *zs         = &fmm_param->zs[3 * nexpmax * level]; 
  double *rdplus     = fmm_param->rdplus;
  double *rdminus    = fmm_param->rdminus;

  double complex *work = CALLOC(nexpmax * 16, sizeof(double complex));
  double complex *mw1 = CALLOC(pgsz, sizeof(double complex));
  double complex *mw2 = CALLOC(pgsz, sizeof(double complex));
  double complex *mexpf1 = CALLOC(nexpmax, sizeof(double complex));
  double complex *mexpf2 = CALLOC(nexpmax, sizeof(double complex));
  double complex *temp   = CALLOC(nexpmax, sizeof(double complex));

  assert(work != NULL);
  assert(mw1 != NULL);
  assert(mw2 != NULL);
  assert(mexpf1 != NULL);
  assert(mexpf2 != NULL);
  assert(temp != NULL); 

  // Process z direction exponential expansions
  double complex *mexuall  = work;
  double complex *mexu1234 = work + nexpmax;
  double complex *mexdall  = mexu1234 + nexpmax;
  double complex *mexd5678 = mexdall + nexpmax;

  make_ulist(1, uall, nuall, xuall, yuall, mexuall, level);
  make_ulist(1, u1234, nu1234, x1234, y1234, mexu1234, level); 
  make_dlist(0, dall, ndall, xdall, ydall, mexdall, level);
  make_dlist(0, d5678, nd5678, x5678, y5678, mexd5678, level);
  
  if (tbox->child[0]) {
    fmm_box_t *child = tbox->child[0];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2];
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1];
      iexp1++;
    }
    
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
    
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[1]) {
    fmm_box_t *child = tbox->child[1]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * conj(xs[3 * j]);
      iexp1++; 
    }
      
    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
      
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * xs[3 * j];
      exponential_to_local_p1(temp, mexpf2, level);
      iexp2++;
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }
  
  if (tbox->child[2]) {
    fmm_box_t *child = tbox->child[2]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * conj(ys[3 * j]);
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
    }
      
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
    
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[3]) {
    fmm_box_t *child = tbox->child[3]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 2] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
    }

    if (nu1234) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexu1234[j] * zs[3 * j + 1] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
    }
      
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);

    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[4]) {
    fmm_box_t *child = tbox->child[4]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;
    
    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1];
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
    
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2];
      iexp2++;
    }
    
    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1];
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }
  
  if (tbox->child[5]) {
    fmm_box_t *child = tbox->child[5]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * xs[3 * j];
      iexp2++;
    }
      
    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * xs[3 * j];
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[6]) {
    fmm_box_t *child = tbox->child[6]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * ys[3 * j];
      iexp2++;
    }
    
    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
    }
      
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  }

  if (tbox->child[7]) {
    fmm_box_t *child = tbox->child[7]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nuall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexuall[j] * zs[3 * j + 1] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (ndall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexdall[j] * zs[3 * j + 2] * xs[3 * j] * ys[3 * j];
      iexp2++;
    }
      
    if (nd5678) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexd5678[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j];
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw1[j];
    }
  } 

  // Process y direction exponential expansions
  double complex *mexnall  = work;
  double complex *mexn1256 = mexnall + nexpmax;
  double complex *mexn12   = mexn1256 + nexpmax;
  double complex *mexn56   = mexn12 + nexpmax;
  double complex *mexsall  = mexn56 + nexpmax;
  double complex *mexs3478 = mexsall + nexpmax;
  double complex *mexs34   = mexs3478 + nexpmax;
  double complex *mexs78   = mexs34 + nexpmax;

  make_ulist(3, nall, nnall, xnall, ynall, mexnall, level);
  make_ulist(3, n1256, nn1256, x1256, y1256, mexn1256, level);
  make_ulist(3, n12, nn12, x12, y12, mexn12, level);
  make_ulist(3, n56, nn56, x56, y56, mexn56, level);
  make_dlist(2, sall, nsall, xsall, ysall, mexsall, level);
  make_dlist(2, s3478, ns3478, x3478, y3478, mexs3478, level);
  make_dlist(2, s34, ns34, x34, y34, mexs34, level);
  make_dlist(2, s78, ns78, x78, y78, mexs78, level);

  if (tbox->child[0]) {
    fmm_box_t *child = tbox->child[0];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0; 

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2];
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1];
      iexp1++;
    }
    
    if (nn12) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn12[j] * zs[3 * j + 1];
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[1]) {
    fmm_box_t *child = tbox->child[1];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * conj(ys[3 * j]);
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
    }

    if (nn12) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn12[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
    
    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[2]) {
    fmm_box_t *child = tbox->child[2]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1];
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2];
      iexp2++;
    }

    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1];
      iexp2++;
    }

    if (ns34) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs34[j] * zs[3 * j + 1];
      iexp2++;
    }
    
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[3]) {
    fmm_box_t *child = tbox->child[3];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
      
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * ys[3 * j];
      iexp2++;
    }
    
    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
    }

    if (ns34) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs34[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
    }

    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[4]) {
    fmm_box_t *child = tbox->child[4];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * conj(xs[3 * j]);
      iexp1++;
    }

    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp1++;
    }

    if (nn56) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn56[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp1++;
    }
    
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
    
    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 1] * xs[3 * j];
      exponential_to_local_p1(temp, mexpf2, level);
      iexp2++;
    }
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }
  
  if (tbox->child[5]) {
    fmm_box_t *child = tbox->child[5];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 2] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
    }
      
    if (nn1256) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn1256[j] * zs[3 * j + 1] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
    }
    
    if (nn56) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexn56[j] * zs[3 * j+1] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
    }

    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
      
    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3*j + 1] * xs[3 * j] * ys[3 * j];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }

    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }
  
  if (tbox->child[6]) {
    fmm_box_t *child = tbox->child[6];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
    
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * xs[3 * j];
      iexp2++;
    }
      
    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * xs[3 * j];
      iexp2++;
    }
    
    if (ns78) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs78[j] * zs[3 * j + 1] * xs[3 * j];
      iexp2++;
    }
    
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  } 

  if (tbox->child[7]) {
    fmm_box_t *child = tbox->child[7]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;
  
    if (nnall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexnall[j] * zs[3 * j + 1] * conj(xs[3 * j] * ys[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);	
    }
      
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nsall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexsall[j] * zs[3 * j + 2] * xs[3 * j] * ys[3 * j];
      iexp2++;
    }
    
    if (ns3478) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs3478[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j];
      iexp2++;
    }
    
    if (ns78) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexs78[j] * zs[3 * j + 1] * xs[3 * j] * ys[3 * j];
      iexp2++;
    }
    
    if (iexp2 > 0)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      roty2z(mw1, rdplus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }
  
  // Process x direction exponential expansions
  double complex *mexeall  = work;
  double complex *mexe1357 = mexeall + nexpmax;
  double complex *mexe13   = mexe1357 + nexpmax;
  double complex *mexe57   = mexe13 + nexpmax;
  double complex *mexe1    = mexe57 + nexpmax;
  double complex *mexe3    = mexe1 + nexpmax;
  double complex *mexe5    = mexe3 + nexpmax;
  double complex *mexe7    = mexe5 + nexpmax;
  double complex *mexwall  = mexe7 + nexpmax;
  double complex *mexw2468 = mexwall + nexpmax;
  double complex *mexw24   = mexw2468 + nexpmax;
  double complex *mexw68   = mexw24 + nexpmax;
  double complex *mexw2    = mexw68 + nexpmax;
  double complex *mexw4    = mexw2 + nexpmax;
  double complex *mexw6    = mexw4 + nexpmax;
  double complex *mexw8    = mexw6 + nexpmax;

  make_ulist(5, eall, neall, xeall, yeall, mexeall, level);
  make_ulist(5, e1357, ne1357, x1357, y1357, mexe1357, level);
  make_ulist(5, e13, ne13, x13, y13, mexe13, level);
  make_ulist(5, e57, ne57, x57, y57, mexe57, level);
  make_ulist(5, e1, ne1, x1, y1, mexe1, level);
  make_ulist(5, e3, ne3, x3, y3, mexe3, level);
  make_ulist(5, e5, ne5, x5, y5, mexe5, level);
  make_ulist(5, e7, ne7, x7, y7, mexe7, level);
  make_dlist(4, wall, nwall, xwall, ywall, mexwall, level);
  make_dlist(4, w2468, nw2468, x2468, y2468, mexw2468, level);
  make_dlist(4, w24, nw24, x24, y24, mexw24, level);
  make_dlist(4, w68, nw68, x68, y68, mexw68, level);
  make_dlist(4, w2, nw2, x2, y2, mexw2, level);
  make_dlist(4, w4, nw4, x4, y4, mexw4, level);
  make_dlist(4, w6, nw6, x6, y6, mexw6, level);
  make_dlist(4, w8, nw8, x8, y8, mexw8, level);

  if (tbox->child[0]) {
    fmm_box_t *child = tbox->child[0];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2];
      iexp1++;
    }
    
    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1];
      iexp1++;
    }
    
    if (ne13) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe13[j] * zs[3 * j + 1];
      iexp1++;
    }
      
    if (ne1) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1[j] * zs[3 * j + 1];
      iexp1++;
    }
    
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
    
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }
  
  if (tbox->child[1]) {
    fmm_box_t *child = tbox->child[1]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1];
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
      
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2];
      iexp2++;
    }
      
    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1];
      iexp2++;
    }
    
    if (nw24) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw24[j] * zs[3 * j + 1];
      iexp2++;
    }
    
    if (nw2) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2[j] * zs[3 * j + 1];
      iexp2++;
    }
    
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[2]) {
    fmm_box_t *child = tbox->child[2];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (ne13) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe13[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (ne3) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe3[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
      
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[3]) {
    fmm_box_t *child = tbox->child[3];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;
   
    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * conj(ys[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
      
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * ys[3 * j];
      iexp2++;
    }
      
    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
    }
      
    if (nw24) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw24[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
    }
    
    if (nw4) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw4[j] * zs[3 * j + 1] * ys[3 * j];
      iexp2++;
    }
    
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[4]) {
    fmm_box_t *child = tbox->child[4];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * xs[3 * j];
      iexp1++;
    }
    
    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * xs[3 * j];
      iexp1++;
    }
    
    if (ne57) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe57[j] * zs[3 * j + 1] * xs[3 * j];
      iexp1++;
    }
    
    if (ne5) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe5[j] * zs[3 * j + 1] * xs[3 * j];
      iexp1++;
    }
      
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
      
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }
  
  if (tbox->child[5]) {
    fmm_box_t *child = tbox->child[5];
    double complex *local = &child->expansion[0];
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * xs[3 * j];
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
      
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
    
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * conj(xs[3 * j]);
      iexp2++;
    }
    
    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp2++;
    }
    
    if (nw68) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw68[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp2++;
    }
      
    if (nw6) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw6[j] * zs[3 * j + 1] * conj(xs[3 * j]);
      iexp2++;
    }
      
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }
  
  if (tbox->child[6]) {
    fmm_box_t *child = tbox->child[6];
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
      
    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 2] * xs[3 * j] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (ne1357) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe1357[j] * zs[3 * j + 1] * xs[3 * j] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (ne57) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe57[j] * zs[3 * j + 1] * xs[3 * j] * conj(ys[3 * j]);
      iexp1++;
    }
      
    if (ne7) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexe7[j] * zs[3 * j + 1] * xs[3 * j] * conj(ys[3 * j]);
      iexp1++;
    }
    
    if (iexp1)
      exponential_to_local_p1(temp, mexpf1, level);
      
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 1] * conj(xs[3 * j]) * ys[3 * j];
      iexp2++;
      exponential_to_local_p1(temp, mexpf2, level);
    }
      
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  if (tbox->child[7]) {
    fmm_box_t *child = tbox->child[7]; 
    double complex *local = &child->expansion[0]; 
    int iexp1 = 0, iexp2 = 0;

    if (neall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexeall[j] * zs[3 * j + 1] * xs[3 * j] * conj(ys[3 * j]);
      iexp1++;
      exponential_to_local_p1(temp, mexpf1, level);
    }
      
    for (int j = 0; j < nexptotp; j++)
      temp[j] = 0;
      
    if (nwall) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] = mexwall[j] * zs[3 * j + 2] * conj(xs[3 * j]) * ys[3 * j];
      iexp2++;
    }
      
    if (nw2468) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw2468[j] * zs[3 * j + 1] * conj(xs[3 * j]) * ys[3 * j];
      iexp2++;
    }
    
    if (nw68) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw68[j] * zs[3 * j + 1] * conj(xs[3 * j]) * ys[3 * j];
      iexp2++;
    }
    
    if (nw8) {
      for (int j = 0; j < nexptotp; j++)
        temp[j] += mexw8[j] * zs[3 * j + 1] * conj(xs[3 * j]) * ys[3 * j];
      iexp2++;
    }
    
    if (iexp2)
      exponential_to_local_p1(temp, mexpf2, level);
    
    if (iexp1 + iexp2) {
      exponential_to_local_p2(iexp2, mexpf2, iexp1, mexpf1, mw1, level);
      rotz2x(mw1, rdminus, mw2);
      for (int j = 0; j < pgsz; j++)
        local[j] += mw2[j];
    }
  }

  FREE(work);
  FREE(mw1);
  FREE(mw2);
  FREE(mexpf1);
  FREE(mexpf2);
  FREE(temp);
}

void exponential_to_local_p1(const double complex *mexpphys,
                             double complex *mexpf, int level) {
  int nlambs = fmm_param->nlambs; 
  int *numfour = fmm_param->numfour; 
  int *numphys = &fmm_param->numphys[nlambs * level]; 
  double complex *fexpback = &fmm_param->fexpback[15000 * level]; 

  int nftot = 0;
  int nptot = 0;
  int next = 0;

  for (int i = 0; i < nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha / 2;
    mexpf[nftot] = 0;
    for (int ival = 0; ival < nalpha2; ival++)
      mexpf[nftot] += 2.0 * creal(mexpphys[nptot + ival]);
    mexpf[nftot] /= nalpha;

    for (int mm = 2; mm < numfour[i]; mm += 2) {
      mexpf[nftot + mm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
        double rtmp = 2 * creal(mexpphys[nptot + ival]);
        mexpf[nftot + mm] += fexpback[next] * rtmp;
        next++;
      }
      mexpf[nftot + mm] /= nalpha;
    }

    for (int mm = 1; mm < numfour[i]; mm += 2) {
      mexpf[nftot + mm] = 0;
      for (int ival = 0; ival < nalpha2; ival++) {
        double complex ztmp = 2 * cimag(mexpphys[nptot + ival]) * _Complex_I;
        mexpf[nftot + mm] += fexpback[next] * ztmp;
        next++;
      }
      mexpf[nftot + mm] /= nalpha;
    }
    nftot += numfour[i];
    nptot += numphys[i] / 2;
  }
}

void exponential_to_local_p2(int iexpu, const double complex *mexpu,
                             int iexpd, const double complex *mexpd,
                             double complex *local, int level) {
  int pterms   = fmm_param->pterms;
  int nlambs   = fmm_param->nlambs;
  int nexpmax  = fmm_param->nexpmax; 
  int pgsz     = fmm_param->pgsz; 
  double beta  = fmm_param->betascale[level]; 
  double *rlsc = &fmm_param->rlsc[pgsz * nlambs * level]; 
  int nexptot  = fmm_param->vnexptot[level]; 
  int *numfour = fmm_param->numfour; 
  double *whts = fmm_param->whts; 
  double *ytop = fmm_param->ytop; 

  double pi = acos(-1); 

  double complex *mexpplus = CALLOC(nexpmax, sizeof(double complex));
  double complex *mexpminus = CALLOC(nexpmax, sizeof(double complex));
  double complex *zeye = CALLOC(pterms + 1, sizeof(double complex));
  assert(mexpplus != NULL);
  assert(mexpminus != NULL);
  assert(zeye != NULL); 

  zeye[0] = 1.0;
  for (int i = 1; i <= pterms; i++)
    zeye[i] = zeye[i - 1] * _Complex_I;

  for (int i = 0; i < pgsz; i++)
    local[i] = 0.0;

  if (iexpu <= 0) {
    for (int nm = 0; nm < nexptot; nm++) {
      mexpplus[nm] = mexpd[nm];
      mexpminus[nm] = mexpd[nm];
    }
  } else if (iexpd <= 0) {
    for (int nm = 0; nm < nexptot; nm++) {
      mexpplus[nm] = mexpu[nm];
      mexpminus[nm] = -mexpu[nm];
    }
  } else {
    for (int nm = 0; nm < nexptot; nm++) {
      mexpplus[nm] = mexpd[nm] + mexpu[nm];
      mexpminus[nm] = mexpd[nm] - mexpu[nm];
    }
  }

  int ntot = 1;

  for (int nl = 0; nl < nlambs; nl++) {
    int mmax = numfour[nl] - 1;
    int offset1 = nl * pgsz;
    for (int mth = 0; mth <= mmax; mth += 2) {
      int ncurrent = ntot + mth - 1;
      int offset = mth * (pterms + 1);
      double complex temp = whts[nl] * mexpplus[ncurrent];

      for (int nm = mth; nm <= pterms; nm += 2) 
        local[nm + offset] += rlsc[nm + offset + offset1] * temp;
      
      temp = whts[nl] * mexpminus[ncurrent];

      for (int nm = mth + 1; nm <= pterms; nm += 2) 
        local[nm + offset] += rlsc[nm + offset + offset1] * temp;
    }

    for (int mth = 1; mth <= mmax; mth += 2) {
      int ncurrent = ntot + mth - 1;
      int offset = mth * (pterms + 1);
      double complex temp = whts[nl] * mexpminus[ncurrent];
      for (int nm = mth; nm <= pterms; nm += 2) 
        local[nm + offset] += rlsc[nm + offset + offset1] * temp;

      temp = whts[nl] * mexpplus[ncurrent];
      for (int nm = mth + 1; nm <= pterms; nm += 2) 
        local[nm + offset] += rlsc[nm + offset + offset1] * temp;
    }
    ntot += numfour[nl];
  }

  double rscale = pi / beta / 2.0; 

  for (int mth = 0; mth <= pterms; mth++) {
    int offset = mth * (pterms + 1);
    double complex temp = zeye[mth] * rscale;
    for (int nm = mth; nm <= pterms; nm++) {
      local[nm + offset] *= temp * ytop[nm + offset];
    }
  }

  FREE(mexpplus);
  FREE(mexpminus);
  FREE(zeye);
}

void local_to_local(fmm_box_t *tbox) {
  const double complex var[5] = 
    {0, 1 - _Complex_I, -1 - _Complex_I, -1 + _Complex_I, 1 + _Complex_I};
  const double arg = sqrt(2.0) / 2.0;
  const int yifl[8] = {3, 4, 2, 1, 3, 4, 2, 1}; 

  int level  = tbox->level; 
  int pterms = fmm_param->pterms;
  int pgsz   = fmm_param->pgsz; 
  int dcpgsz = fmm_param->dcpgsz; 
  double *dc = &fmm_param->dcd[dcpgsz * level]; 
  
  double complex *localn = CALLOC(pgsz, sizeof(double complex));
  double complex *marray = CALLOC(pgsz, sizeof(double complex));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex)); 
  assert(localn != NULL);
  assert(marray != NULL);
  assert(ephi != NULL); 

  double complex *plocal = &tbox->expansion[0]; 

  for (int i = 0; i < 8; i++) {
    fmm_box_t *child = tbox->child[i]; 
    if (child != NULL) {
      int ifl = yifl[i];
      double *rd = (i < 4 ? fmm_param->rdmsq3 : fmm_param->rdsq3);

      ephi[0] = 1.0;
      ephi[1] = arg * var[ifl];
      for (int ell = 1; ell <= pterms; ell++)
        ephi[ell + 1] = ephi[ell] * ephi[1];

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          localn[ell + offset] = conj(ephi[m]) * plocal[ell + offset];
        }
      }

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = localn[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp++) {
            int index1 = ell + mp*(pterms + 1);
            marray[index] += localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
        }
      }
      
      for (int k = 0; k <= pterms; k++) {
        int offset = k * (pterms + 1);
        for (int j = k; j <= pterms; j++) {
          int index = j + offset;
          int offset1 = k + j * (pterms + 1);
          localn[index] = 0;
          for (int ell = k; ell <= pterms; ell++)
            localn[index] += marray[ell + offset] * dc[offset1 + ell * pgsz];
        }
      }

      for (int m = 0; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = ell + offset;
          marray[index] = localn[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = mp * (pterms + 1) + ell;
            marray[index] -= localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
          
          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = mp * (pterms + 1) + ell;
            marray[index] += localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
        }
      }
      
      for (int m = 1; m <= pterms; m += 2) {
        int offset = m * (pterms + 1);
        int offset1 = (m + pterms) * pgsz;
        int offset2 = (-m + pterms) * pgsz;
        for (int ell = m; ell <= pterms; ell++) {
          int index = offset + ell;
          marray[index] = -localn[ell] * rd[ell + offset1];
          for (int mp = 1; mp <= ell; mp += 2) {
            int index1 = mp * (pterms + 1) + ell;
            marray[index] += localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
          
          for (int mp = 2; mp <= ell; mp += 2) {
            int index1 = mp * (pterms + 1) + ell;
            marray[index] -= localn[index1] * rd[index1 + offset1] +
              conj(localn[index1]) * rd[index1 + offset2];
          }
        }
      }

      for (int m = 0; m <= pterms; m++) {
        int offset = m * (pterms + 1);
        for (int ell = m; ell <= pterms; ell++) {
          localn[ell + offset] = ephi[m] * marray[ell + offset];
        }
      }

      double complex *clocal = &child->expansion[0]; 
      for (int m =0; m < pgsz; m++)
        clocal[m] += localn[m];
    }
  }

  FREE(localn);
  FREE(marray);
  FREE(ephi);
}

void local_to_target(fmm_box_t *tbox) {
  int level = tbox->level; 
  int pgsz = fmm_param->pgsz; 
  int pterms = fmm_param->pterms; 
  double myscale = fmm_param->sfactor[level]; 
  double beta = fmm_param->beta; 

  double complex *local = &tbox->expansion[0]; 

  double h = fmm_dag->size / (1 << (level + 1));
  double center[3]; 
  center[0] = fmm_dag->corner[0] + (2 * tbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * tbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * tbox->idz + 1) * h;


  double *p = CALLOC((pterms + 2) * (pterms + 2), sizeof(double));
  double complex *ptemp = CALLOC((pterms + 2) * (pterms + 2), 
                                 sizeof(double complex));
  double *bi = CALLOC(pterms + 2, sizeof(double));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex)); 

  assert(p != NULL);
  assert(ptemp != NULL);
  assert(bi != NULL);
  assert(ephi != NULL); 

  int ntargets = tbox->npts; 
  int addr = tbox->addr; 
  const double precision = 1.0e-14;

  for (int i = 0; i < ntargets; i++) {
    int ptr = addr + i;
    double *pot = &potential_sorted[ptr]; 
    double *field = &field_sorted[3 * ptr]; 
    double complex cpz, fieldtemp[3] = {0};

    double rx = targets_sorted[3 * ptr] - center[0]; 
    double ry = targets_sorted[3 * ptr + 1] - center[1]; 
    double rz = targets_sorted[3 * ptr + 2] - center[2]; 
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj);
    double d = sqrt(rr);

    double ctheta = (d <= precision ? 0.0: rz / d);
    ephi[0] = (proj <= precision*d ? 1.0 : rx / proj + _Complex_I * ry / proj);
    d *= beta;

    for (int ell = 0; ell <= pterms; ell++)
      ephi[ell + 1] = ephi[ell] * ephi[0];

    int ncalc;
    in(myscale, d, pterms + 1, bi, &ncalc);
    lgndr(pterms + 1, ctheta, p);

    *pot += creal(local[0])*bi[0];
    ptemp[0] = bi[0];

    for (int ell = 1; ell <= pterms; ell++) {
      double cp = bi[ell] * p[ell];
      ptemp[ell] = cp;
      *pot += creal(local[ell]) * cp;
    }

    for (int m = 1; m <= pterms; m++) {
      int offset1 = m * (pterms + 2);
      int offset2 = m * (pterms + 1);
      for (int ell = m; ell <= pterms; ell++) {
        ptemp[ell + offset1] = bi[ell] * p[ell + offset1] * ephi[m - 1];
        *pot += 2.0 * creal(local[ell + offset2] * ptemp[ell + offset1]);
      }
    }
    
    ptemp[pterms + 1] = 0;
    for (int m = 1; m <= pterms + 1; m++) {
      int offset1 = m * (pterms + 2);
      ptemp[pterms + 1 + offset1] = 
        bi[pterms + 1] * p[pterms + 1 + offset1] * ephi[m - 1];
    }

    double rpotz = local[0] * myscale;
    cpz = rpotz * ptemp[1 + pterms + 2];
    fieldtemp[0] = conj(cpz);
    fieldtemp[1] = rpotz * ptemp[1];
    fieldtemp[2] = -cpz;
    
    rpotz = local[1] / 3.0;
    cpz = (rpotz * myscale) * ptemp[2 + pterms + 2];
    fieldtemp[0] += conj(cpz);
    fieldtemp[1] += rpotz * (ptemp[0] / myscale + 2.0 * ptemp[2] * myscale);
    fieldtemp[2] -= cpz;

    for (int ell = 2; ell <= pterms; ell++) {
      rpotz = local[ell] / (2 * ell + 1);
      cpz = rpotz * (ptemp[ell + pterms + 1] / myscale - 
                     ptemp[ell + pterms + 3] * myscale);
      fieldtemp[0] -= conj(cpz);
      fieldtemp[1] += rpotz * (creal(ptemp[ell - 1]) * ell / myscale + 
                               (ell + 1) * creal(ptemp[ell + 1]) * myscale);
      fieldtemp[2] += cpz;
    }

    field[0] -= beta * creal(fieldtemp[0] - fieldtemp[2]) / 2;
    field[1] -= beta * creal((fieldtemp[0] + fieldtemp[2]) * _Complex_I ) / 2.0;
    field[2] += beta * creal(fieldtemp[1]);
    
    cpz = local[pterms + 2] / 3.0;
    fieldtemp[0] = cpz * ((creal(ptemp[0]) / myscale - 
                           creal(ptemp[2]) * myscale) * 2.0);
    cpz *= myscale;
    fieldtemp[1] = cpz * ptemp[pterms + 4];
    fieldtemp[2] = cpz * ptemp[2 + 2 * (pterms + 2)];

    for (int m = 1; m <= pterms - 2; m++) {
      int offset1 = m * (pterms + 1);
      int offset2 = (m - 1) * (pterms + 2);
      int offset3 = offset2 + (pterms + 2);
      int offset4 = offset3 + (pterms + 2);
      for (int ell = m + 2; ell <= pterms; ell++) {
        cpz = local[ell + offset1] / (2 * ell + 1);
        fieldtemp[0] += cpz * (ptemp[ell - 1 + offset2] * (ell + m - 1) * 
                               (ell + m) / myscale -
                               ptemp[ell + 1 + offset2] * (ell - m + 1) * 
                               (ell - m + 2) * myscale);
        fieldtemp[1] += cpz * (ptemp[ell - 1 + offset3] * (ell + m) / myscale 
                               + ptemp[ell + 1 + offset3] * (ell - m + 1) * 
                               myscale);
        fieldtemp[2] += cpz * (ptemp[ell - 1 + offset4] / myscale -
                               ptemp[ell + 1 + offset4] * myscale);
      }
    }

    for (int ell = 2; ell <= pterms; ell++) {
      int offset1 = (ell - 1) * (pterms + 1);
      int offset2 = (ell - 2) * (pterms + 2) + ell;
      int offset3 = offset2 + (pterms + 2);
      int offset4 = offset3 + (pterms + 2);
      cpz = local[ell + offset1] / (2 * ell + 1);
      fieldtemp[0] += cpz * (ptemp[-1 + offset2] * (ell + ell - 2) * 
                             (ell + ell - 1) / myscale -
                             ptemp[1 + offset2] * 6 * myscale);
      fieldtemp[1] += cpz * (ptemp[-1 + offset3] * (ell + ell - 1) / myscale +
                             ptemp[1 + offset3] * 2 * myscale);
      fieldtemp[2] -= cpz * ptemp[1 + offset4] * myscale;
    }

    for (int ell = 2; ell <= pterms; ell++) {
      int offset1 = ell * (pterms + 1);
      int offset2 = (ell - 1) * (pterms + 2) + ell;
      int offset3 = offset2 + (pterms + 2);
      int offset4 = offset3 + (pterms + 2);
      cpz = local[ell + offset1] / (2 * ell + 1);
      fieldtemp[0] += cpz * (ptemp[-1 + offset2] * (ell + ell - 1) * 
                             (ell + ell) / myscale -
                             ptemp[1 + offset2] * 2.0 * myscale);
      fieldtemp[1] += cpz * ptemp[1 + offset3] * myscale;
      fieldtemp[2] -= cpz * ptemp[1 + offset4] * myscale;
    }
    
    field[0] -= beta * creal(fieldtemp[0] - fieldtemp[2]);
    field[1] -= beta * creal((fieldtemp[0] + fieldtemp[2]) * _Complex_I);
    field[2] += 2.0 * beta * creal(fieldtemp[1]);
  }

  FREE(p);
  FREE(ptemp);
  FREE(bi);
  FREE(ephi);
}

void source_to_local(fmm_box_t *tbox, fmm_box_t *sbox) {
  int pgsz = fmm_param->pgsz;
  int pterms = fmm_param->pterms; 
  double beta = fmm_param->beta; 
  double *ytop = fmm_param->ytop; 

  int nsources = sbox->npts;
  int addr = sbox->addr; 

  double h = fmm_dag->size / (1 << (tbox->level + 1));
  double center[3];
  center[0] = fmm_dag->corner[0] + (2 * tbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * tbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * tbox->idz + 1) * h;
  double complex *local = &tbox->expansion[0]; 
  double scale = fmm_param->sfactor[tbox->level]; 

  double *p = CALLOC((pterms + 2) * (pterms + 2), sizeof(double));
  double *bk = CALLOC(pterms + 2, sizeof(double));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex));

  assert(p != NULL);
  assert(bk != NULL);
  assert(ephi != NULL); 

  const double precis = 1e-14; 

  for (int i = 0; i < nsources; i++) {
    int ptr = addr + i;
    double rx = sources_sorted[3 * ptr]     - center[0]; 
    double ry = sources_sorted[3 * ptr + 1] - center[1];
    double rz = sources_sorted[3 * ptr + 2] - center[2]; 
    double proj = rx * rx + ry * ry;
    double rr = proj + rz * rz;
    proj = sqrt(proj); 
    double d = sqrt(rr); 

    double ctheta = (d <= precis ? 1.0 : rz / d);
    ephi[0] = (proj <= precis * d ? 1.0 : rx / proj - _Complex_I * ry / proj);

    for (int ell = 1; ell < pterms; ell++) 
      ephi[ell] = ephi[ell - 1] * ephi[0]; 

    double rk = d * beta; 
    int ncalc;
    kn(scale, rk, pterms, bk, &ncalc); 
    lgndr(pterms, ctheta, p); 

    local[0] += charges_sorted[ptr] * bk[0]; 
    
    for (int ell = 1; ell <= pterms; ell++)       
      local[ell] += charges_sorted[ptr] * p[ell] * bk[ell] * ytop[ell]; 

    for (int m = 1; m <= pterms; m++) {
      int offset = m * (pterms + 1);
      for (int ell = m; ell <= pterms; ell++) {
        int index = ell + offset;
        local[index] += charges_sorted[ptr] * bk[ell] * ytop[index] * 
          p[index] * ephi[m - 1]; 
      }
    }
  }

  FREE(p);
  FREE(bk);
  FREE(ephi);
}

void multipole_to_target(fmm_box_t *tbox, fmm_box_t *sbox) {
  int pgsz = fmm_param->pgsz;
  int pterms = fmm_param->pterms; 
  double beta = fmm_param->beta; 

  double complex *multipole = &sbox->expansion[0]; 
  double h = fmm_dag->size / (1 << (sbox->level + 1)); 
  double center[3];
  center[0] = fmm_dag->corner[0] + (2 * sbox->idx + 1) * h;
  center[1] = fmm_dag->corner[1] + (2 * sbox->idy + 1) * h;
  center[2] = fmm_dag->corner[2] + (2 * sbox->idz + 1) * h;
  double scale = fmm_param->scale * fmm_dag->size / (1 << (sbox->level)); 

  double *p = CALLOC((pterms + 2) * (pterms + 2), sizeof(double));
  double complex *ptemp = CALLOC((pterms + 2) * (pterms + 2), 
                                 sizeof(double complex));
  double *bk = CALLOC(pterms + 2, sizeof(double));
  double complex *ephi = CALLOC(pterms + 2, sizeof(double complex));

  assert(p != NULL);
  assert(ptemp != NULL);
  assert(bk != NULL);
  assert(ephi != NULL); 

  int ntargets = tbox->npts; 
  int addr = tbox->addr; 
  const double precis = 1e-14; 

  for (int i = 0; i < ntargets; i++) {
    double complex fieldtemp[3] = {0}; 
    int ptr = addr + i; 
    double *pot = &potential_sorted[ptr]; 
    double *field = &field_sorted[3 * ptr]; 

    double rx = targets_sorted[3 * ptr] - center[0];
    double ry = targets_sorted[3 * ptr + 1] - center[1];
    double rz = targets_sorted[3 * ptr + 2] - center[2]; 
    double proj = rx * rx + ry * ry; 
    double rr = proj + rz * rz; 
    proj = sqrt(proj);
    double d = sqrt(rr); 
    double ctheta = (d <= precis ? 0.0 : rz / d);
    ephi[0] = (proj <= precis * d ? 1.0 : rx / proj + _Complex_I * ry / proj);

    for (int ell = 1; ell <= pterms + 1; ell++) 
      ephi[ell] = ephi[ell - 1] * ephi[0]; 

    double rk = d * beta; 

    int ncalc; 
    kn(scale, rk, pterms + 1, bk, &ncalc); 
    lgndr(pterms + 1, ctheta, p); 

    ptemp[0] = bk[0];
    *pot += multipole[0] * bk[0];

    for (int ell = 1; ell <= pterms; ell++) {
      ptemp[ell] = bk[ell] * p[ell];
      *pot += multipole[ell] * ptemp[ell];
    }

    for (int m = 1; m <= pterms; m++) {
      int offset1 = m * (pterms + 2);
      int offset2 = m * (pterms + 1);
      for (int ell = m; ell <= pterms; ell++) {
        ptemp[ell + offset1] = bk[ell] * p[ell + offset1] * ephi[m - 1];
        *pot += 2.0 * creal(multipole[ell + offset2] * ptemp[ell + offset1] );
      }
    }

    ptemp[pterms + 1] = 0;
    for (int m = 1; m <= pterms + 1; m++) {
      int offset1 = m * (pterms + 2) + pterms + 1;
      ptemp[offset1] = bk[pterms + 1] * p[offset1] * ephi[m - 1];  
    }

    double rpotz = multipole[0] / scale;
    double complex cpz = rpotz * ptemp[1 + pterms + 2]; //ptemp(1,1)
    fieldtemp[0] = -conj(cpz);
    fieldtemp[1] = -rpotz * ptemp[1]; //ptemp(1,0)
    fieldtemp[2] = cpz;

    // contribution from n = 1, m = 0
    rpotz = multipole[1] / 3.0;
    cpz = rpotz / scale * ptemp[2 + pterms + 2]; // ptemp(2,1)
    fieldtemp[0] -= conj(cpz);
    fieldtemp[1]  -= rpotz * (ptemp[0] * scale + 2.0 * ptemp[2] / scale);
    fieldtemp[2] += cpz;

    // contribution from n > 1, m = 0
    for (int ell = 2; ell <= pterms; ell++) {
      rpotz = multipole[ell] / (2 * ell + 1);
      cpz = rpotz * (ptemp[ell - 1 + pterms + 2] * scale - 
                     ptemp[ell + 1 + pterms + 2] / scale);
      fieldtemp[0] += conj(cpz);
      fieldtemp[1] -= rpotz * (creal(ptemp[ell - 1]) * ell * scale + 
                               (ell + 1) * creal(ptemp[ell + 1]) / scale);
      fieldtemp[2] -= cpz;
    }

    // update force 
    field[0] -= beta * creal(fieldtemp[0] - fieldtemp[2]) / 2.0;
    field[1] -= beta * creal((fieldtemp[0] + fieldtemp[2]) * _Complex_I) / 2.0;
    field[2] += beta * creal(fieldtemp[1]);

    // cotribution from n = 1, m = 1
    cpz = multipole[1 + pterms + 1] / 3.0;
    fieldtemp[0] = -cpz * (ptemp[0] * scale - ptemp[2] / scale) * 2.0;
    fieldtemp[1] = -cpz * ptemp[2 + pterms + 2] / scale; //ptemp(2,1)
    fieldtemp[2] = cpz * ptemp[2 + 2 * (pterms + 2)] / scale; //ptemp(2,2)

    // contribution from n = 2, pterms; m = 1, n - 2
    for (int m = 1; m <= pterms - 2; m++) {
      int offset1 = m * (pterms + 1);
      int offset2 = (m - 1) * (pterms + 2);
      int offset3 = offset2 + (pterms + 2);
      int offset4 = offset3 + (pterms + 2);
      for (int ell = m + 2; ell <= pterms; ell++) {
        cpz = multipole[ell + offset1] / (2 * ell + 1);
        fieldtemp[0] -= cpz * (ptemp[ell - 1 + offset2] * (ell + m-1) * 
                               (ell + m) * scale - 
                               ptemp[ell + 1 + offset2] * (ell-m + 1) * 
                               (ell-m + 2) / scale);
        fieldtemp[1] -= cpz * (ptemp[ell - 1 + offset3] * (ell + m) * scale + 
                               ptemp[ell + 1 + offset3] * (ell - m + 1) / scale);
        fieldtemp[2] -= cpz * (ptemp[ell - 1 + offset4] * scale - 
                               ptemp[ell + 1 + offset4] / scale);
      }
    }

    for (int ell = 2; ell <= pterms; ell++) {
      // m = n-1
      cpz = multipole[ell + (ell - 1) * (pterms + 1)] / (2 * ell + 1);
      fieldtemp[0] -= cpz * (ptemp[ell - 1 + (ell - 2) * (pterms + 2)] * 
                             (ell + ell - 2) * (ell + ell - 1) * scale - 
                             ptemp[ell + 1 + (ell - 2) * (pterms + 2)] * 
                             6 / scale );

      fieldtemp[1] -= cpz * (ptemp[ell - 1 + (ell - 1) * (pterms + 2)] * 
                             (ell + ell - 1) * scale + 
                             ptemp[ell + 1 + (ell - 1) * (pterms + 2)] * 
                             2 / scale );

      fieldtemp[2] += cpz * ptemp[ell + 1 + ell * (pterms + 2)] / scale;
    }

    for (int ell = 2; ell <= pterms; ell++) {
      // m = n
      int offset1 = ell * (pterms + 1);
      int offset2 = (ell - 1)*(pterms + 2) + ell;
      int offset3 = offset2 + (pterms + 2);
      int offset4 = offset3 + (pterms + 2);

      cpz = multipole[ell + offset1] / (2 * ell + 1);
      fieldtemp[0] -= cpz *(ptemp[-1  +  offset2] * (ell + ell - 1) * 
                            (ell + ell) * scale - 
                            ptemp[1 + offset2] * 2 / scale);
      fieldtemp[1] -= cpz * ptemp[1 + offset3] / scale;
      fieldtemp[2] += cpz * ptemp[1 + offset4] / scale;
    }

    field[0] -= beta * creal(fieldtemp[0] - fieldtemp[2]);
    field[1] -= beta * creal((fieldtemp[0] + fieldtemp[2]) * _Complex_I);
    field[2] += 2.0 * beta * creal(fieldtemp[1]);
  }

  FREE(p);
  FREE(ptemp);
  FREE(bk);
  FREE(ephi);
}

void direct_evaluation(const fmm_box_t *tbox, const fmm_box_t *sbox) {
  int start1 = tbox->addr, num1 = tbox->npts, end1 = start1 + num1 - 1;
  int start2 = sbox->addr, num2 = sbox->npts, end2 = start2 + num2 - 1;
  double beta = fmm_param->beta; 
  const double pi_2 = acos(0); 

  for (int i = start1; i <= end1; i++) {
    int i3 = i * 3;
    double pot = 0, fx = 0, fy = 0, fz = 0; 

#ifdef icc
    #pragma simd reduction(+:pot)
    #pragma simd reduction(+:fx)
    #pragma simd reduction(+:fy)
    #pragma simd reduction(+:fz)
#endif
    for (int j = start2; j <= end2; j++) {
      int j3 = j * 3; 
      double rx = targets_sorted[i3]     - sources_sorted[j3];
      double ry = targets_sorted[i3 + 1] - sources_sorted[j3 + 1];
      double rz = targets_sorted[i3 + 2] - sources_sorted[j3 + 2]; 
      double q = charges_sorted[j]; 
      double rr = rx * rx + ry * ry + rz * rz;
      double rdis = sqrt(rr) * beta; 
      double expr = cexp(-rdis); 
      if (rdis) {
        double term1 = q / rdis * expr * pi_2; 
        pot += term1;
        double term2 = -term1 * (1 + rdis) / rr;
        fx += rx * term2;
        fy += ry * term2;
        fz += rz * term2;
      }
    }

    potential_sorted[i]  += pot;
    field_sorted[i3]     += fx;
    field_sorted[i3 + 1] += fy;
    field_sorted[i3 + 2] += fz;
  }
}

#endif 
