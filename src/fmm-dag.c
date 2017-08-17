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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "cilk.h"
#include "fmm-dag.h"

const int xoff[] = {0, 1, 0, 1, 0, 1, 0, 1};
const int yoff[] = {0, 0, 1, 1, 0, 0, 1, 1};
const int zoff[] = {0, 0, 0, 0, 1, 1, 1, 1};

int part_thres, nslev; 

fmm_dag_t *fmm_dag; 

int *source_swap, *source_record, *target_swap, *target_record; 

static inline 
void update_list(fmm_box_t **list, int *nlist, int *xoff, int *yoff,
                 fmm_box_t *entry, int ix, int iy) {
  list[*nlist] = entry;
  xoff[*nlist] = ix;
  yoff[*nlist] = iy;
  (*nlist)++;
}

fmm_dag_t *construct_dag(const Real_t *sources, const int nsources, 
                         const Real_t *targets, const int ntargets, 
                         const int s) {
  fmm_dag = CALLOC(1, sizeof(fmm_dag_t)); 
  assert(fmm_dag != NULL);
  
  fmm_dag->nslev = 0; 
  fmm_dag->nsboxes = 0; 
  fmm_dag->ntlev = 0; 
  fmm_dag->ntboxes = 0;

  fmm_dag->mapsrc = CALLOC(nsources, sizeof(int)); 
  fmm_dag->maptar = CALLOC(ntargets, sizeof(int)); 
  source_swap   = CALLOC(nsources, sizeof(int)); 
  source_record = CALLOC(nsources, sizeof(int)); 
  target_swap   = CALLOC(ntargets, sizeof(int)); 
  target_record = CALLOC(ntargets, sizeof(int)); 

  assert(fmm_dag->mapsrc != NULL);
  assert(fmm_dag->maptar != NULL);
  assert(source_swap != NULL);
  assert(source_record != NULL);
  assert(target_swap != NULL);
  assert(target_record != NULL);

#pragma simd
  for (int i = 0; i < nsources; i++) 
    fmm_dag->mapsrc[i] = i; 

#pragma simd
  for (int i = 0; i < ntargets; i++) 
    fmm_dag->maptar[i] = i; 

#ifdef USING_CILK
  CILK_C_REDUCER_MIN(xminR, double, sources[0]); 
  CILK_C_REDUCER_MAX(xmaxR, double, sources[0]);
  CILK_C_REDUCER_MIN(yminR, double, sources[1]);
  CILK_C_REDUCER_MAX(ymaxR, double, sources[1]);
  CILK_C_REDUCER_MIN(zminR, double, sources[2]);
  CILK_C_REDUCER_MAX(zmaxR, double, sources[2]); 
  CILK_C_REGISTER_REDUCER(xminR);
  CILK_C_REGISTER_REDUCER(xmaxR);
  CILK_C_REGISTER_REDUCER(yminR);
  CILK_C_REGISTER_REDUCER(ymaxR);
  CILK_C_REGISTER_REDUCER(zminR);
  CILK_C_REGISTER_REDUCER(zmaxR);

  cilk_for (int i = 1; i <= nsources; i++) {
    int j = 3 * i; 
    REDUCER_VIEW(xminR) = fmin(sources[j],     REDUCER_VIEW(xminR)); 
    REDUCER_VIEW(xmaxR) = fmax(sources[j],     REDUCER_VIEW(xmaxR)); 
    REDUCER_VIEW(yminR) = fmin(sources[j + 1], REDUCER_VIEW(yminR)); 
    REDUCER_VIEW(ymaxR) = fmax(sources[j + 1], REDUCER_VIEW(ymaxR)); 
    REDUCER_VIEW(zminR) = fmin(sources[j + 2], REDUCER_VIEW(zminR)); 
    REDUCER_VIEW(zmaxR) = fmax(sources[j + 2], REDUCER_VIEW(zmaxR)); 
  }

  cilk_for (int i = 1; i <= ntargets; i++) {
    int j = 3 * i; 
    REDUCER_VIEW(xminR) = fmin(targets[j],     REDUCER_VIEW(xminR)); 
    REDUCER_VIEW(xmaxR) = fmax(targets[j],     REDUCER_VIEW(xmaxR)); 
    REDUCER_VIEW(yminR) = fmin(targets[j + 1], REDUCER_VIEW(yminR)); 
    REDUCER_VIEW(ymaxR) = fmax(targets[j + 1], REDUCER_VIEW(ymaxR)); 
    REDUCER_VIEW(zminR) = fmin(targets[j + 2], REDUCER_VIEW(zminR)); 
    REDUCER_VIEW(zmaxR) = fmax(targets[j + 2], REDUCER_VIEW(zmaxR)); 
  }

  double xmin = xminR.value;
  double xmax = xmaxR.value;
  double ymin = yminR.value;
  double ymax = ymaxR.value;
  double zmin = zminR.value;
  double zmax = zmaxR.value; 
#else 
  double xmin = sources[0], xmax = xmin; 
  double ymin = sources[1], ymax = ymin; 
  double zmin = sources[2], zmax = zmin; 

 for (int i = 1; i < nsources; i++) {
    int j = 3 * i; 
    xmin = fmin(sources[j],     xmin);
    xmax = fmax(sources[j],     xmax);
    ymin = fmin(sources[j + 1], ymin);
    ymax = fmax(sources[j + 1], ymax);
    zmin = fmin(sources[j + 2], zmin);
    zmax = fmax(sources[j + 2], zmax);
  }

  for (int i = 0; i < ntargets; i++) {
    int j = 3 * i;
    xmin = fmin(targets[j],     xmin);
    xmax = fmax(targets[j],     xmax);
    ymin = fmin(targets[j + 1], ymin);
    ymax = fmax(targets[j + 1], ymax);
    zmin = fmin(targets[j + 2], zmin);
    zmax = fmax(targets[j + 2], zmax);
  }
#endif

  fmm_dag->size = fmax(fmax(xmax - xmin, ymax - ymin), zmax - zmin); 
  fmm_dag->corner[0] = (xmax + xmin - fmm_dag->size) * 0.5; 
  fmm_dag->corner[1] = (ymax + ymin - fmm_dag->size) * 0.5;
  fmm_dag->corner[2] = (zmax + zmin - fmm_dag->size) * 0.5; 

  fmm_box_t *source_root = CALLOC(1, sizeof(fmm_box_t)); 
  assert(source_root != NULL); 
  source_root->level = 0; 
  source_root->idx = 0; 
  source_root->idy = 0; 
  source_root->idz = 0;
  source_root->npts = nsources; 
  source_root->addr = 0; 
  fmm_dag->source_root = source_root; 

  fmm_box_t *target_root = CALLOC(1, sizeof(fmm_box_t)); 
  assert(target_root != NULL); 
  target_root->level = 0; 
  target_root->idx = 0; 
  target_root->idy = 0; 
  target_root->idz = 0; 
  target_root->npts = ntargets; 
  target_root->addr = 0; 
  fmm_dag->target_root = target_root; 

  part_thres = s; 
  fmm_dag->nsboxes = 1; 
  fmm_dag->ntboxes = 1; 

  CILK_SPAWN partition_box(fmm_dag->source_root, 'S', sources, 
                           &(fmm_dag->nslev), &(fmm_dag->nsboxes)); 
  CILK_SPAWN partition_box(fmm_dag->target_root, 'T', targets, 
                           &(fmm_dag->ntlev), &(fmm_dag->ntboxes)); 
  CILK_SYNC;

  nslev = fmm_dag->nslev; 

  FREE(source_swap);
  FREE(source_record);
  FREE(target_swap);
  FREE(target_record); 

  return fmm_dag; 
}

void partition_box(fmm_box_t *box, const char type, const Real_t *points, 
                   int *level, int *nboxes) {
  int begin, npoints, *imap, *swap, *record; 
 
  begin   = box->addr; 
  npoints = box->npts; 

  if (type == 'S') {
    imap   = &(fmm_dag->mapsrc[begin]);
    swap   = &source_swap[begin]; 
    record = &source_record[begin]; 
  } else {
    imap   = &(fmm_dag->maptar[begin]); 
    swap   = &target_swap[begin]; 
    record = &target_record[begin]; 
  }

  // Compute box center coordinate
  double h  = fmm_dag->size/(1 << (box->level + 1)); 
  double xc = fmm_dag->corner[0] + (2 * box->idx + 1) * h;
  double yc = fmm_dag->corner[1] + (2 * box->idy + 1) * h;
  double zc = fmm_dag->corner[2] + (2 * box->idz + 1) * h;
  
  unsigned child_parts[8] = {0}, addrs[8] = {0}, assigned[8] = {0}; 

  // Determine which box each particle belongs to 
#pragma simd
  for (int i = 0; i < npoints; i++) {
    int j = 3 * imap[i]; 
    int bin = 4 * (points[j + 2] > zc) + 
      2 * (points[j + 1] > yc) + (points[j] > xc); 
    record[i] = bin;
  } 

  for (int i = 0; i < npoints; i++)
    child_parts[record[i]]++; 

  // Compute the offsets
  addrs[1] = addrs[0] + child_parts[0]; 
  addrs[2] = addrs[1] + child_parts[1]; 
  addrs[3] = addrs[2] + child_parts[2]; 
  addrs[4] = addrs[3] + child_parts[3]; 
  addrs[5] = addrs[4] + child_parts[4]; 
  addrs[6] = addrs[5] + child_parts[5]; 
  addrs[7] = addrs[6] + child_parts[6]; 

  // Assign particle to boxes
  for (int i = 0; i < npoints; i++) {
    int bin = record[i]; 
    int offset = addrs[bin] + assigned[bin]++; 
    swap[offset] = imap[i]; 
  }

  // Update this box's portion of global mapping 
#pragma simd
  for (int i = 0; i < npoints; i++) 
    imap[i] = swap[i]; 

  // Create new boxes
  int levels_below[8] = {0}, nboxes_below[8] = {0};
  box->nchild = (child_parts[0] > 0) + (child_parts[1] > 0) + 
    (child_parts[2] > 0) + (child_parts[3] > 0) + (child_parts[4] > 0) + 
    (child_parts[5] > 0) + (child_parts[6] > 0) + (child_parts[7] > 0); 

  CILK_FOR (int i = 0; i < 8; i++) {
    if (child_parts[i] > 0) {
      levels_below[i] = 0; 
      nboxes_below[i] = 1; 
      fmm_box_t *cbox = CALLOC(1, sizeof(fmm_box_t)); 
      assert(cbox != NULL); 
      cbox->level = box->level + 1; 
      cbox->parent = box; 
      cbox->idx = 2 * box->idx + xoff[i]; 
      cbox->idy = 2 * box->idy + yoff[i];
      cbox->idz = 2*box->idz + zoff[i]; 
      cbox->npts = child_parts[i]; 
      cbox->addr = box->addr + addrs[i]; 
      box->child[i] = cbox; 

      if (child_parts[i] > part_thres) 
        partition_box(cbox, type, points, &levels_below[i], &nboxes_below[i]);
    }
  }
  

  *nboxes = *nboxes + nboxes_below[0] + nboxes_below[1] + nboxes_below[2] + 
    nboxes_below[3] + nboxes_below[4] + nboxes_below[5] + nboxes_below[6] + 
    nboxes_below[7]; 

  int maxlevel = levels_below[0]; 
  for (int i = 1; i < 8; i++) 
    maxlevel = (maxlevel > levels_below[i] ? maxlevel : levels_below[i]); 
  
  *level = 1 + maxlevel; 
}

void destruct_dag(fmm_dag_t *fmm_dag) {
  FREE(fmm_dag->mapsrc); 
  FREE(fmm_dag->maptar); 
  CILK_SPAWN delete_box(fmm_dag->source_root); 
  CILK_SPAWN delete_box(fmm_dag->target_root); 
  FREE(fmm_dag); 
  CILK_SYNC; 
}

void delete_box(fmm_box_t *box) {
  if (box != NULL) {
    FREE(box->expansion); 
    if (box->nchild) {
      CILK_FOR (int i = 0; i < 8; i++) 
        delete_box(box->child[i]); 
    }
    FREE(box); 
  } 
}

bool is_adjacent(const fmm_box_t *box1, const fmm_box_t *box2) {
  int dim = 1 << (box2->level - box1->level); 

  return ((box2->idx >= dim*box1->idx - 1) &&
          (box2->idx <= dim*box1->idx + dim) &&
          (box2->idy >= dim*box1->idy - 1) &&
          (box2->idy <= dim*box1->idy + dim) &&
          (box2->idz >= dim*box1->idz - 1) &&
          (box2->idz <= dim*box1->idz + dim)); 
}

void 
build_merged_list2(const fmm_box_t *tbox, 
                   fmm_box_t *uall[36], int *nuall, int *xuall, int *yuall,
                   fmm_box_t *u1234[16], int *nu1234, int *x1234, int *y1234,
                   fmm_box_t *dall[36], int *ndall, int *xdall, int *ydall, 
                   fmm_box_t *d5678[16], int *nd5678, int *x5678, int *y5678,
                   fmm_box_t *nall[24], int *nnall, int *xnall, int *ynall, 
                   fmm_box_t *n1256[8], int *nn1256, int *x1256, int *y1256,
                   fmm_box_t *n12[4], int *nn12, int *x12, int *y12,
                   fmm_box_t *n56[4], int *nn56, int *x56, int *y56, 
                   fmm_box_t *sall[24], int *nsall, int *xsall, int *ysall, 
                   fmm_box_t *s3478[8], int *ns3478, int *x3478, int *y3478,
                   fmm_box_t *s34[4], int *ns34, int *x34, int *y34, 
                   fmm_box_t *s78[4], int *ns78, int *x78, int *y78, 
                   fmm_box_t *eall[16], int *neall, int *xeall, int *yeall, 
                   fmm_box_t *e1357[4], int *ne1357, int *x1357, int *y1357,
                   fmm_box_t *e13[2], int *ne13, int *x13, int *y13,
                   fmm_box_t *e57[2], int *ne57, int *x57, int *y57,
                   fmm_box_t *e1[3], int *ne1, int *x1, int *y1, 
                   fmm_box_t *e3[3], int *ne3, int *x3, int *y3, 
                   fmm_box_t *e5[3], int *ne5, int *x5, int *y5, 
                   fmm_box_t *e7[3], int *ne7, int *x7, int *y7, 
                   fmm_box_t *wall[16], int *nwall, int *xwall, int *ywall, 
                   fmm_box_t *w2468[4], int *nw2468, int *x2468, int *y2468,
                   fmm_box_t *w24[2], int *nw24, int *x24, int *y24,
                   fmm_box_t *w68[2], int *nw68, int *x68, int *y68,
                   fmm_box_t *w2[3], int *nw2, int *x2, int *y2, 
                   fmm_box_t *w4[3], int *nw4, int *x4, int *y4, 
                   fmm_box_t *w6[3], int *nw6, int *x6, int *y6, 
                   fmm_box_t *w8[3], int *nw8, int *x8, int *y8) {
  *nuall = *nu1234 = 0;
  *ndall = *nd5678 = 0;
  *nnall = *nn1256 = *nn12 = *nn56 = 0;
  *nsall = *ns3478 = *ns34 = *ns78 = 0;
  *neall = *ne1357 = *ne13 = *ne57 = *ne1 = *ne3 = *ne5 = *ne7 = 0;
  *nwall = *nw2468 = *nw24 = *nw68 = *nw2 = *nw4 = *nw6 = *nw8 = 0;

  int tidx = tbox->idx; 
  int tidy = tbox->idy; 
  int tidz = tbox->idz; 
  int nlist5 = tbox->nlist5; 

  for (int i = 0; i < nlist5; i++) {
    fmm_box_t *sbox = tbox->list5[i]; 
    int sidx = sbox->idx;
    int sidy = sbox->idy;
    int sidz = sbox->idz; 
    int offset = 9 * (sidz - tidz) + 3 * (sidy - tidy) + sidx - tidx + 13; 

    switch (offset) {
    case 0:
      // [-1][-1][-1], update dall, sall, wall, d5678, s34, w2 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], -2, -2);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], -1,-2);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], -2, -1);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], -1,-1);
      if (sbox->child[4]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[4], -1, -2);
      if (sbox->child[5]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[5], -1, -1);
      if (sbox->child[6]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[6], 1, -1);
      if (sbox->child[7]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[7], -1, -1);
        update_list(s34, ns34, x34, y34, sbox->child[7], -1, -1);
        update_list(w2, nw2, x2, y2, sbox->child[7], 1, -1);
      }
      break;
    case 1://[0][-1][-1], update dall, sall, d5678, s34 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], 0, -2);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], 1, -2);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], 0, -1);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], 1,-1);
      if (sbox->child[4]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[4], -1, 0);
      if (sbox->child[5]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[5], -1, 1);     
      if (sbox->child[6]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[6], 0, -1);
        update_list(s34, ns34, x34, y34, sbox->child[6], -1,0);
      }
      if (sbox->child[7]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[7], 1, -1);
        update_list(s34, ns34, x34, y34, sbox->child[7], -1, 1);
      }
      break;
    case 2: //[1][-1][-1], update dall, sall, eall, d5678, s34, e1 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], 2, -2);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], 3, -2);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], 2, -1);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], 3, -1);
      if (sbox->child[4]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[4], -1, 2);
      if (sbox->child[5]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[5], -1, 3);
      if (sbox->child[6]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[6], 2, -1);
        update_list(s34, ns34, x34, y34, sbox->child[6], -1, 2);
        update_list(e1, ne1, x1, y1, sbox->child[6], 1, -1);
      }
      if (sbox->child[7]) 
        update_list(eall, neall, xeall, yeall, sbox->child[7], 1, -1);
      break;
    case 3: //[-1][0][-1], update dall, wall, d5678, w24 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], -2, 0);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], -1, 0);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], -2, 1);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], -1, 1);
      if (sbox->child[4]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[4], 1, 0);
      if (sbox->child[5]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[5], -1, 0);
        update_list(w24, nw24, x24, y24, sbox->child[5], 1, 0);
      }
      if (sbox->child[6]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[6], 1, 1);
      if (sbox->child[7]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[7], -1, 1);
        update_list(w24, nw24, x24, y24, sbox->child[7], 1, 1);
      }
      break;
    case 4://[0][0][-1], update dall and d5678 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], 0, 0);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], 1, 0);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], 0, 1);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], 1, 1);
      if (sbox->child[4]) 
        update_list(d5678, nd5678, x5678, y5678, sbox->child[4], 0, 0);
      if (sbox->child[5]) 
        update_list(d5678, nd5678, x5678, y5678, sbox->child[5], 1, 0);
      if (sbox->child[6]) 
        update_list(d5678, nd5678, x5678, y5678, sbox->child[6], 0, 1);
      if (sbox->child[7]) 
        update_list(d5678, nd5678, x5678, y5678, sbox->child[7], 1, 1);
      break;
    case 5: // [1][0][-1], update dall, d5678, eall, e13 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], 2, 0);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], 3, 0);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], 2, 1);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], 3, 1);
      if (sbox->child[4]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[4], 2, 0);
        update_list(e13, ne13, x13, y13, sbox->child[4], 1, 0);
      }
      if (sbox->child[5]) 
        update_list(eall, neall, xeall, yeall, sbox->child[5], 1, 0);
      if (sbox->child[6]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[6], 2, 1);
        update_list(e13, ne13, x13, y13, sbox->child[6], 1, 1);
      }
      if (sbox->child[7]) 
        update_list(eall, neall, xeall, yeall, sbox->child[7], 1, 1);
      break;
    case 6://[-1][1][-1], update dall, nall, wall, d5678, n12, w4 list3
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], -2, 2);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], -1, 2);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], -2, 3);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], -1, 3);
      if (sbox->child[4]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[4], 1, 2);
      if (sbox->child[5]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[5], -1, 2);
        update_list(n12, nn12, x12, y12, sbox->child[5], -1, -1);
        update_list(w4, nw4, x4, y4, sbox->child[5], 1, 2);
      }
      if (sbox->child[6]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[6], -1, -2);
      if (sbox->child[7]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[7], -1, -1);
      break;
    case 7: //[0][1][-1], update dallm d5678, nall, n12 lists 
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], 0, 2);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], 1, 2);
      if (sbox->child[2]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[2], 0, 3);
      if (sbox->child[3]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[3], 1, 3);
      if (sbox->child[4]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[4], 0, 2);
        update_list(n12, nn12, x12, y12, sbox->child[4], -1, 0);
      }
      if (sbox->child[5]) {
        update_list(d5678, nd5678,x5678, y5678, sbox->child[5], 1, 2);
        update_list(n12, nn12, x12, y12, sbox->child[5], -1, 1);
      }
      if (sbox->child[6]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[6], -1, 0);
      if (sbox->child[7])
        update_list(nall, nnall, xnall, ynall, sbox->child[7], -1, 1);
      break;
    case 8: //[1][1][-1], update dall, d5678, nall, eall, n12, e3 lists
      if (sbox->child[0]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[0], 2,2);
      if (sbox->child[1]) 
        update_list(dall, ndall, xdall, ydall, sbox->child[1], 3, 2);
      if (sbox->child[2])
        update_list(dall, ndall, xdall, ydall, sbox->child[2], 2, 3);
      if (sbox->child[3])
        update_list(dall, ndall, xdall, ydall, sbox->child[3], 3,3);
      if (sbox->child[4]) {
        update_list(d5678, nd5678, x5678, y5678, sbox->child[4], 2, 2);
        update_list(n12, nn12, x12, y12, sbox->child[4], -1, 2);
        update_list(e3, ne3, x3, y3, sbox->child[4], 1, 2);
      }
      if (sbox->child[5]) 
        update_list(eall, neall, xeall, yeall, sbox->child[5], 1, 2);
      if (sbox->child[6]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[6], -1, 2);
      if (sbox->child[7]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[7], -1, 3);
      break;
    case 9: // [-1][-1][0], update sall, wall, s3478, w2, w6 lists
      if (sbox->child[0]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[0], 0, -2);
      if (sbox->child[1])
        update_list(sall, nsall, xsall, ysall, sbox->child[1], 0, -1);
      if (sbox->child[2])
        update_list(wall, nwall, xwall, ywall, sbox->child[2], 0, -1);
      if (sbox->child[3]) {
        update_list(s3478, ns3478, x3478, y3478, sbox->child[3], 0, -1);
        update_list(w2, nw2, x2, y2, sbox->child[3], 0, -1);
        update_list(w6, nw6, x6, y6, sbox->child[3], 0, -1);
      }
      if (sbox->child[4]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[4], 1, -2);
      if (sbox->child[5])
        update_list(sall, nsall, xsall, ysall, sbox->child[5], 1, -1);
      if (sbox->child[6]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[6], -1, -1);
      if (sbox->child[7]) {
        update_list(s3478, ns3478, x3478, y3478, sbox->child[7], 1, -1);
        update_list(w2, nw2, x2, y2, sbox->child[7], -1, -1);
        update_list(w6, nw6, x6, y6, sbox->child[7], -1, -1);
      }
      break;
    case 10://[0][-1][0], update sall, s3478 lists
      if (sbox->child[0]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[0], 0, 0);
      if (sbox->child[1])
        update_list(sall, nsall, xsall, ysall, sbox->child[1], 0, 1);
      if (sbox->child[2])
        update_list(s3478, ns3478, x3478, y3478, sbox->child[2], 0, 0);
      if (sbox->child[3])
        update_list(s3478, ns3478, x3478, y3478, sbox->child[3], 0, 1);
      if (sbox->child[4]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[4], 1, 0);
      if (sbox->child[5])
        update_list(sall, nsall, xsall, ysall, sbox->child[5], 1, 1);
      if (sbox->child[6])
        update_list(s3478, ns3478, x3478, y3478, sbox->child[6], 1, 0);
      if (sbox->child[7])
        update_list(s3478, ns3478, x3478, y3478, sbox->child[7], 1, 1);
      break;
    case 11: //[1][-1][0], update eall, sall, s3478, e1, e5 lists
      if (sbox->child[0])
        update_list(sall, nsall, xsall, ysall, sbox->child[0], 0, 2);
      if (sbox->child[1])
        update_list(sall, nsall, xsall, ysall, sbox->child[1], 0, 3);
      if (sbox->child[2]) {
        update_list(s3478, ns3478, x3478, y3478, sbox->child[2], 0, 2);
        update_list(e1, ne1, x1, y1, sbox->child[2], 0, -1);
        update_list(e5, ne5, x5, y5, sbox->child[2], 0, -1);
      }
      if (sbox->child[3]) 
        update_list(eall, neall, xeall, yeall, sbox->child[3], 0, -1);
      if (sbox->child[4])
        update_list(sall, nsall, xsall, ysall, sbox->child[4], 1, 2);
      if (sbox->child[5]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[5], 1, 3);
      if (sbox->child[6]) {
        update_list(s3478, ns3478, x3478, y3478, sbox->child[6], 1, 2);
        update_list(e1, ne1, x1, y1, sbox->child[6], -1, -1);
        update_list(e5, ne5, x5, y5, sbox->child[6], -1, -1);
      }
      if (sbox->child[7]) 
        update_list(eall, neall, xeall, yeall, sbox->child[7], -1, -1);
      break;
    case 12: // [-1][0][0], update wall, w2468 lists
      if (sbox->child[0])
        update_list(wall, nwall, xwall, ywall, sbox->child[0], 0, 0);
      if (sbox->child[1])
        update_list(w2468, nw2468, x2468, y2468, sbox->child[1], 0, 0);
      if (sbox->child[2])
        update_list(wall, nwall, xwall, ywall, sbox->child[2], 0, 1);
      if (sbox->child[3])
        update_list(w2468, nw2468, x2468, y2468, sbox->child[3], 0, 1);
      if (sbox->child[4])
        update_list(wall, nwall, xwall, ywall, sbox->child[4], -1, 0);
      if (sbox->child[5])
        update_list(w2468, nw2468, x2468, y2468, sbox->child[5], -1, 0);
      if (sbox->child[6])
        update_list(wall, nwall, xwall, ywall, sbox->child[6], -1, 1);
      if (sbox->child[7])
        update_list(w2468, nw2468, x2468, y2468, sbox->child[7], -1, 1);
      break;
    case 13: //[0][0][0], nothing here
      break;
    case 14: //[1][0][0], update eall, e1357 lists
      if (sbox->child[0])
        update_list(e1357, ne1357, x1357, y1357, sbox->child[0], 0, 0);
      if (sbox->child[1])
        update_list(eall, neall, xeall, yeall, sbox->child[1], 0, 0);
      if (sbox->child[2])
        update_list(e1357, ne1357, x1357, y1357, sbox->child[2], 0, 1);
      if (sbox->child[3])
        update_list(eall, neall, xeall, yeall, sbox->child[3], 0, 1);
      if (sbox->child[4]) 
        update_list(e1357, ne1357, x1357, y1357, sbox->child[4], -1, 0);
      if (sbox->child[5])
        update_list(eall, neall, xeall, yeall, sbox->child[5], -1, 0);
      if (sbox->child[6])
        update_list(e1357, ne1357, x1357, y1357, sbox->child[6], -1, 1);
      if (sbox->child[7])
        update_list(eall, neall, xeall, yeall, sbox->child[7], -1, 1);
      break;
    case 15://[-1][1][0], update wall, nall, n1256, w4, w8 lists
      if (sbox->child[0])
        update_list(wall, nwall, xwall, ywall, sbox->child[0], 0, 2);
      if (sbox->child[1]) {
        update_list(n1256, nn1256, x1256, y1256, sbox->child[1], 0, -1);
        update_list(w4, nw4, x4, y4, sbox->child[1], 0, 2);
        update_list(w8, nw8, x8, y8, sbox->child[1], 0, 2);
      }
      if (sbox->child[2])
        update_list(nall, nnall, xnall, ynall, sbox->child[2], 0, -2);
      if (sbox->child[3])
        update_list(nall, nnall, xnall, ynall, sbox->child[3], 0, -1);
      if (sbox->child[4])
        update_list(wall, nwall, xwall, ywall, sbox->child[4], -1, 2);
      if (sbox->child[5]) {
        update_list(n1256, nn1256, x1256, y1256, sbox->child[5], 1, -1);
        update_list(w4, nw4, x4, y4, sbox->child[5], -1, 2);
        update_list(w8, nw8, x8, y8, sbox->child[5], -1, 2);
      }
      if (sbox->child[6]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[6], 1, -2);
      if (sbox->child[7])
        update_list(nall, nnall, xnall, ynall, sbox->child[7], 1, -1);
      break;
    case 16: //[0][1][0], update nall, n1256 lists
      if (sbox->child[0])
        update_list(n1256, nn1256, x1256, y1256, sbox->child[0], 0, 0);
      if (sbox->child[1])
        update_list(n1256, nn1256, x1256, y1256, sbox->child[1], 0, 1);
      if (sbox->child[2])
        update_list(nall, nnall, xnall, ynall, sbox->child[2], 0, 0);
      if (sbox->child[3])
        update_list(nall, nnall, xnall, ynall, sbox->child[3], 0, 1);
      if (sbox->child[4])
        update_list(n1256, nn1256, x1256, y1256, sbox->child[4], 1, 0);
      if (sbox->child[5])
        update_list(n1256, nn1256, x1256, y1256, sbox->child[5], 1, 1);
      if (sbox->child[6])
        update_list(nall, nnall, xnall, ynall, sbox->child[6], 1, 0);
      if (sbox->child[7])
        update_list(nall, nnall, xnall, ynall, sbox->child[7], 1, 1);
      break;
    case 17: //[1][1][0], update nall, n1256, eall, e3, e7 lists
      if (sbox->child[0]) {
        update_list(n1256, nn1256, x1256, y1256, sbox->child[0], 0, 2);
        update_list(e3, ne3, x3, y3, sbox->child[0], 0, 2);
        update_list(e7, ne7, x7, y7, sbox->child[0], 0, 2);
      }
      if (sbox->child[1]) 
        update_list(eall, neall, xeall, yeall, sbox->child[1], 0, 2);
      if (sbox->child[2])
        update_list(nall, nnall, xnall, ynall, sbox->child[2], 0, 2);
      if (sbox->child[3])
        update_list(nall, nnall, xnall, ynall, sbox->child[3], 0, 3);
      if (sbox->child[4]) {
        update_list(n1256, nn1256, x1256, y1256, sbox->child[4], 1, 2);
        update_list(e3, ne3, x3, y3, sbox->child[4], -1, 2);
        update_list(e7, ne7, x7, y7, sbox->child[4], -1, 2);
      }
      if (sbox->child[5])
        update_list(eall, neall, xeall, yeall, sbox->child[5], -1, 2);
      if (sbox->child[6])
        update_list(nall, nnall, xnall, ynall, sbox->child[6], 1, 2);
      if (sbox->child[7])
        update_list(nall, nnall, xnall, ynall, sbox->child[7], 1, 3);
      break;
    case 18: //[-1][-1][1], update sall, wall, u1234, s78, w6, uall lists
      if (sbox->child[0]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[0], 2, -2);
      if (sbox->child[1]) 
        update_list(sall, nsall, xsall, ysall, sbox->child[1], 2, -1);
      if (sbox->child[2]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[2], -2, -1);
      if (sbox->child[3]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[3], -1,-1);
        update_list(s78, ns78, x78, y78, sbox->child[3], 2, -1);
        update_list(w6, nw6, x6, y6, sbox->child[3], -2, -1);
      }
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], -2, -2);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], -1, -2);
      if (sbox->child[6])
        update_list(uall, nuall, xuall, yuall, sbox->child[6], -2, -1);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], -1, -1);
      break;
    case 19: //[0][-1][1], update sall, u1234, s78, uall lists
      if (sbox->child[0])
        update_list(sall, nsall, xsall, ysall, sbox->child[0], 2, 0);
      if (sbox->child[1])
        update_list(sall, nsall, xsall, ysall, sbox->child[1], 2, 1);
      if (sbox->child[2]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[2], 0, -1);
        update_list(s78, ns78, x78, y78, sbox->child[2], 2, 0);
      }
      if (sbox->child[3]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[3], 1, -1);
        update_list(s78, ns78, x78, y78, sbox->child[3], 2, 1);
      }
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], 0, -2);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], 1, -2);
      if (sbox->child[6])
        update_list(uall, nuall, xuall, yuall, sbox->child[6], 0, -1);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], 1, -1);
      break;
    case 20: // [1][-1][1], update sall, eall, u1234, s78, e5, uall lists
      if (sbox->child[0])
        update_list(sall, nsall, xsall, ysall, sbox->child[0], 2, 2);
      if (sbox->child[1])
        update_list(sall, nsall, xsall, ysall, sbox->child[1], 2, 3);
      if (sbox->child[2]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[2], 2, -1);
        update_list(s78, ns78, x78, y78, sbox->child[2], 2, 2);
        update_list(e5, ne5, x5, y5, sbox->child[2], -2, -1);
      }
      if (sbox->child[3])
        update_list(eall, neall, xeall, yeall, sbox->child[3], -2, -1);
      if (sbox->child[4])
        update_list(uall, nuall, xuall, yuall, sbox->child[4], 2, -2);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], 3, -2);
      if (sbox->child[6])
        update_list(uall, nuall, xuall, yuall, sbox->child[6], 2, -1);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], 3, -1);
      break;
    case 21: // [-1][0][1], update wall, u1234, w68, uall lists
      if (sbox->child[0])
        update_list(wall, nwall, xwall, ywall, sbox->child[0], -2, 0);
      if (sbox->child[1]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[1], -1, 0);
        update_list(w68, nw68, x68, y68, sbox->child[1], -2, 0);
      }
      if (sbox->child[2]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[2], -2, 1);
      if (sbox->child[3]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[3], -1, 1);
        update_list(w68, nw68, x68, y68, sbox->child[3], -2, 1);
      }
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], -2, 0);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], -1, 0);
      if (sbox->child[6])
        update_list(uall, nuall, xuall, yuall, sbox->child[6], -2, 1);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], -1, 1);
      break;
    case 22: //[0][0][1], update u1234, uall lists
      if (sbox->child[0]) 
        update_list(u1234, nu1234, x1234, y1234, sbox->child[0], 0, 0);
      if (sbox->child[1])
        update_list(u1234, nu1234, x1234, y1234, sbox->child[1], 1, 0);
      if (sbox->child[2]) 
        update_list(u1234, nu1234, x1234, y1234, sbox->child[2], 0, 1);
      if (sbox->child[3])
        update_list(u1234, nu1234, x1234, y1234, sbox->child[3], 1, 1);
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], 0, 0);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], 1, 0);
      if (sbox->child[6]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[6], 0, 1);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], 1, 1);
      break;
    case 23: // [1][0][1], update u1234, e57, eall, uall lists
      if (sbox->child[0]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[0], 2, 0);
        update_list(e57, ne57, x57, y57, sbox->child[0], -2, 0);
      }
      if (sbox->child[1])
        update_list(eall, neall, xeall, yeall, sbox->child[1], -2, 0);
      if (sbox->child[2]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[2], 2, 1);
        update_list(e57, ne57, x57, y57, sbox->child[2], -2, 1);
      }
      if (sbox->child[3])
        update_list(eall, neall, xeall, yeall, sbox->child[3], -2, 1);
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], 2, 0);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], 3, 0);
      if (sbox->child[6]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[6], 2, 1);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], 3, 1);
      break;
    case 24: // [-1][1][1], update nall, wall, u1234, n56, w8, uall lists
      if (sbox->child[0]) 
        update_list(wall, nwall, xwall, ywall, sbox->child[0], -2, 2);
      if (sbox->child[1]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[1], -1, 2);
        update_list(n56, nn56, x56, y56, sbox->child[1], 2, -1);
        update_list(w8, nw8, x8, y8, sbox->child[1], -2, 2);
      }
      if (sbox->child[2]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[2], 2, -2);
      if (sbox->child[3])
        update_list(nall, nnall, xnall, ynall, sbox->child[3], 2, -1);
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], -2, 2);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], -1, 2);
      if (sbox->child[6]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[6], -2, 3);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], -1, 3);
      break;
    case 25: // [0][1][1], update u1234, nall, n56, nall lists
      if (sbox->child[0]) {
        update_list(u1234, nu1234, x1234, y1234,  sbox->child[0], 0, 2);
        update_list(n56, nn56, x56, y56,  sbox->child[0], 2, 0);
      }
      if (sbox->child[1]) {
        update_list(u1234, nu1234, x1234, y1234,  sbox->child[1],1, 2);
        update_list(n56, nn56, x56, y56,  sbox->child[1], 2, 1);
      }
      if (sbox->child[2]) 
        update_list(nall, nnall, xnall, ynall,  sbox->child[2], 2, 0);
      if (sbox->child[3])
        update_list(nall, nnall, xnall, ynall,  sbox->child[3], 2, 1);
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], 0, 2);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], 1, 2);
      if (sbox->child[6]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[6], 0, 3);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], 1, 3);
      break;
    case 26: // [1][1][1], update u1234, n56, e7, eall, nall, uall lists 
      if (sbox->child[0]) {
        update_list(u1234, nu1234, x1234, y1234, sbox->child[0], 2, 2);
        update_list(n56, nn56, x56, y56, sbox->child[0], 2, 2);
        update_list(e7, ne7, x7, y7, sbox->child[0], -2, 2);
      }
      if (sbox->child[1]) 
        update_list(eall, neall, xeall, yeall, sbox->child[1], -2, 2);
      if (sbox->child[2]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[2], 2, 2);
      if (sbox->child[3]) 
        update_list(nall, nnall, xnall, ynall, sbox->child[3], 2, 3);
      if (sbox->child[4]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[4], 2, 2);
      if (sbox->child[5])
        update_list(uall, nuall, xuall, yuall, sbox->child[5], 3, 2);
      if (sbox->child[6]) 
        update_list(uall, nuall, xuall, yuall, sbox->child[6], 2, 3);
      if (sbox->child[7])
        update_list(uall, nuall, xuall, yuall, sbox->child[7], 3, 3);
      break;
    default:
      break;
    }
  }
}
