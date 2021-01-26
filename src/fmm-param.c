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
#include <math.h>
#include <assert.h>
#include "cilk.h"
#include "fmm-param.h"
#include "fmm-types.h"

void bnlcft(double *c, int p) {
  for (int n = 0; n <= p; n++)
    c[n] = 1.0;

  for (int m = 1; m <= p; m++) {
    int offset = m * (p + 1);
    int offset1 = offset - p - 1;
    c[m + offset] = 1.0;
    for (int n = m + 1; n <= p; n++) 
      c[n + offset] = c[n - 1 + offset] + c[n - 1 + offset1];
  }

  for (int m = 1; m <= p; m++) {
    int offset = m * (p + 1);
    for (int n = m + 1; n <= p; n++) {
      c[n + offset] = sqrt(c[n + offset]);
    }
  }
}

void vwts(fmm_param_t *fmm_param) {
  int nlambs    = fmm_param->nlambs; 
  double *rlams = fmm_param->rlams; 
  double *whts  = fmm_param->whts; 

  if (nlambs == 9) {
    rlams[0] = 0.99273996739714473469540223504736787e-01;
    rlams[1] = 0.47725674637049431137114652301534079e+00;
    rlams[2] = 0.10553366138218296388373573790886439e+01;
    rlams[3] = 0.17675934335400844688024335482623428e+01;
    rlams[4] = 0.25734262935147067530294862081063911e+01;
    rlams[5] = 0.34482433920158257478760788217186928e+01;
    rlams[6] = 0.43768098355472631055818055756390095e+01;
    rlams[7] = 0.53489575720546005399569367000367492e+01;
    rlams[8] = 0.63576578531337464283978988532908261e+01;
    whts[0]  = 0.24776441819008371281185532097879332e+00;
    whts[1]  = 0.49188566500464336872511239562300034e+00;
    whts[2]  = 0.65378749137677805158830324216978624e+00;
    whts[3]  = 0.76433038408784093054038066838984378e+00;
    whts[4]  = 0.84376180565628111640563702167128213e+00;
    whts[5]  = 0.90445883985098263213586733400006779e+00;
    whts[6]  = 0.95378613136833456653818075210438110e+00;
    whts[7]  = 0.99670261613218547047665651916759089e+00;
    whts[8]  = 0.10429422730252668749528766056755558e+01;
  } else if (nlambs == 18) { 
    rlams[0]  = 0.52788527661177607475107009804560221e-01;
    rlams[1]  = 0.26949859838931256028615734976483509e+00;
    rlams[2]  = 0.63220353174689392083962502510985360e+00;
    rlams[3]  = 0.11130756427760852833586113774799742e+01;
    rlams[4]  = 0.16893949614021379623807206371566281e+01;
    rlams[5]  = 0.23437620046953044905535534780938178e+01;
    rlams[6]  = 0.30626998290780611533534738555317745e+01;
    rlams[7]  = 0.38356294126529686394633245072327554e+01;
    rlams[8]  = 0.46542473432156272750148673367220908e+01;
    rlams[9]  = 0.55120938659358147404532246582675725e+01;
    rlams[10] = 0.64042126837727888499784967279992998e+01;
    rlams[11] = 0.73268800190617540124549122992902994e+01;
    rlams[12] = 0.82774009925823861522076185792684555e+01;
    rlams[13] = 0.92539718060248947750778825138695538e+01;
    rlams[14] = 0.10255602723746401139237605093512684e+02;
    rlams[15] = 0.11282088297877740146191172243561596e+02;
    rlams[16] = 0.12334067909676926788620221486780792e+02;
    rlams[17] = 0.13414920240172401477707353478763252e+02;
    whts[0]   = 0.13438265914335215112096477696468355e+00;
    whts[1]   = 0.29457752727395436487256574764614925e+00;
    whts[2]   = 0.42607819361148618897416895379137713e+00;
    whts[3]   = 0.53189220776549905878027857397682965e+00;
    whts[4]   = 0.61787306245538586857435348065337166e+00;
    whts[5]   = 0.68863156078905074508611505734734237e+00;
    whts[6]   = 0.74749099381426187260757387775811367e+00;
    whts[7]   = 0.79699192718599998208617307682288811e+00;
    whts[8]   = 0.83917454386997591964103548889397644e+00;
    whts[9]   = 0.87570092283745315508980411323136650e+00;
    whts[10]  = 0.90792943590067498593754180546966381e+00;
    whts[11]  = 0.93698393742461816291466902839601971e+00;
    whts[12]  = 0.96382546688788062194674921556725167e+00;
    whts[13]  = 0.98932985769673820186653756536543369e+00;
    whts[14]  = 0.10143828459791703888726033255807124e+01;
    whts[15]  = 0.10400365437416452252250564924906939e+01;
    whts[16]  = 0.10681548926956736522697610780596733e+01;
    whts[17]  = 0.11090758097553685690428437737864442e+01;
  }
}

void numthetahalf(fmm_param_t *fmm_param) {
  int *numfour = fmm_param->numfour; 
  int nlambs   = fmm_param->nlambs; 

  if (nlambs == 9) {
    numfour[0] = 2;
    numfour[1] = 4;
    numfour[2] = 4;
    numfour[3] = 6;
    numfour[4] = 6;
    numfour[5] = 4;
    numfour[6] = 6;
    numfour[7] = 4;
    numfour[8] = 2;
  } else if (nlambs == 18) {
    numfour[0]  = 4;
    numfour[1]  = 6;
    numfour[2]  = 6;
    numfour[3]  = 8;
    numfour[4]  = 8;
    numfour[5]  = 8;
    numfour[6]  = 10;
    numfour[7]  = 10;
    numfour[8]  = 10;
    numfour[9]  = 10;
    numfour[10] = 12;
    numfour[11] = 12;
    numfour[12] = 12;
    numfour[13] = 12;
    numfour[14] = 12;
    numfour[15] = 12;
    numfour[16] = 8;
    numfour[17] = 2;
  }
}

#if LAPLACE

fmm_param_t *construct_param(const fmm_config_t *fmm_config, 
                             const fmm_dag_t *fmm_dag) {
  int accuracy = fmm_config->accuracy; 

  fmm_param_t *param = CALLOC(1, sizeof(fmm_param_t)); 
  assert(param != NULL); 

  if (accuracy == 3) {
    param->pterms = 9;
    param->nlambs = 9;
    param->pgsz = 100;
  } else if (accuracy == 6) {
    param->pterms = 18;    
    param->nlambs = 18;
    param->pgsz = 361;
  }

  int pterms = param->pterms; 
  int nlambs = param->nlambs; 
  int pgsz = param->pgsz; 

  param->numphys = CALLOC(nlambs, sizeof(int)); 
  param->numfour = CALLOC(nlambs, sizeof(int)); 
  param->whts    = CALLOC(nlambs, sizeof(double)); 
  param->rlams   = CALLOC(nlambs, sizeof(double)); 
  param->rdplus  = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->rdminus = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->rdsq3   = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->rdmsq3  = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->dc      = CALLOC((2 * pterms + 1)*(2 * pterms + 1)* (2 * pterms + 1), 
                          sizeof(double)); 
  param->ytopc     = CALLOC((pterms + 2) * (pterms + 2), sizeof(double)); 
  param->ytopcs    = CALLOC((pterms + 2) * (pterms + 2), sizeof(double)); 
  param->ytopcsinv = CALLOC((pterms + 2) * (pterms + 2), sizeof(double)); 
  param->rlsc      = CALLOC(pgsz * nlambs, sizeof(double)); 

  assert(param->numphys != NULL);
  assert(param->numfour != NULL);
  assert(param->whts != NULL);
  assert(param->rlams != NULL);
  assert(param->rdplus != NULL);
  assert(param->rdminus != NULL);
  assert(param->rdsq3 != NULL);
  assert(param->rdmsq3 != NULL);
  assert(param->dc != NULL);
  assert(param->ytopc != NULL);
  assert(param->ytopcs != NULL);
  assert(param->ytopcsinv != NULL);
  assert(param->rlsc != NULL);

  frmini(param); 
  rotgen(param); 
  vwts(param); 
  numthetahalf(param); 
  numthetafour(param); 
  rlscini(param); 

  param->nexptot  = 0;
  param->nthmax   = 0; 
  param->nexptotp = 0; 

  for (int i = 1; i <= nlambs; i++) {
    param->nexptot += param->numfour[i - 1]; 
    if (param->numfour[i - 1] > param->nthmax) 
      param->nthmax = param->numfour[i - 1]; 
    param->nexptotp += param->numphys[i - 1]; 
  }
  
  param->nexptotp *= 0.5; 
  param->nexpmax = (param->nexptot > param->nexptotp ? 
                    param->nexptot : param->nexptotp) + 1;

  param->xs = CALLOC(param->nexpmax * 3, sizeof(double complex)); 
  param->ys = CALLOC(param->nexpmax * 3, sizeof(double complex)); 
  param->zs = CALLOC(param->nexpmax * 3, sizeof(double)); 
  param->fexpe    = CALLOC(15000, sizeof(double complex)); 
  param->fexpo    = CALLOC(15000, sizeof(double complex)); 
  param->fexpback = CALLOC(15000, sizeof(double complex)); 

  assert(param->xs != NULL);
  assert(param->ys != NULL);
  assert(param->zs != NULL);
  assert(param->fexpe != NULL);
  assert(param->fexpo != NULL);
  assert(param->fexpback != NULL);

  mkfexp(param); 
  mkexps(param); 

  int nlevel = (fmm_dag->nslev >= fmm_dag->ntlev ? 
                fmm_dag->nslev : fmm_dag->ntlev);

  param->scale = CALLOC(1 + nlevel, sizeof(double));   
  assert(param->scale != NULL); 

  param->scale[0] = 1 / fmm_dag->size; 
  for (int i = 1; i <= nlevel; i++) 
    param->scale[i] = 2 * param->scale[i - 1]; 

  return param; 
}

void destruct_param(fmm_param_t *fmm_param) {
  FREE(fmm_param->xs); 
  FREE(fmm_param->ys); 
  FREE(fmm_param->zs);
  FREE(fmm_param->fexpe);
  FREE(fmm_param->fexpo);
  FREE(fmm_param->fexpback);
  FREE(fmm_param->numphys);
  FREE(fmm_param->numfour);
  FREE(fmm_param->whts); 
  FREE(fmm_param->rlams);
  FREE(fmm_param->rdplus);
  FREE(fmm_param->rdminus);
  FREE(fmm_param->rdsq3); 
  FREE(fmm_param->rdmsq3); 
  FREE(fmm_param->dc); 
  FREE(fmm_param->ytopc);
  FREE(fmm_param->ytopcs);
  FREE(fmm_param->ytopcsinv);
  FREE(fmm_param->rlsc); 
  FREE(fmm_param->scale); 
  FREE(fmm_param); 
}


void frmini(fmm_param_t *fmm_param) {
  double *ytopc     = fmm_param->ytopc; 
  double *ytopcs    = fmm_param->ytopcs; 
  double *ytopcsinv = fmm_param->ytopcsinv; 
  int pterms        = fmm_param->pterms; 

  double *factorial = CALLOC(3 + 2 * pterms, sizeof(double));
  assert(factorial != NULL); 

  double d = 1.0;
  factorial[0] = d;
  for (int ell = 1; ell <= 2 * pterms + 2; ell++) {
    d *= sqrt(ell);
    factorial[ell] = d;
  }

  ytopcs[0] = 1.0;
  ytopcsinv[0] = 1.0;
  for (int m = 0; m <= pterms + 1; m++) {
    int offset = m * (pterms + 2);
    for (int ell = m; ell <= pterms + 1; ell++) {
      ytopc[ell + offset] = factorial[ell - m] / factorial[ell + m];
      ytopcsinv[ell + offset] = factorial[ell - m] * factorial[ell + m];
      ytopcs[ell + offset] = 1.0 / ytopcsinv[ell + offset];
    }
  }


  FREE(factorial);
}

void rotgen(fmm_param_t *fmm_param) {
  double *dc      = fmm_param->dc; 
  double *rdplus  = fmm_param->rdplus; 
  double *rdminus = fmm_param->rdminus; 
  double *rdsq3   = fmm_param->rdsq3; 
  double *rdmsq3  = fmm_param->rdmsq3; 
  int pterms      = fmm_param->pterms; 
  int pgsz        = fmm_param->pgsz; 

  bnlcft(dc, 2*pterms); 

  double theta = acos(0); 
  fstrtn(pterms, rdplus, dc, theta, pgsz); 
  
  theta = -theta; 
  fstrtn(pterms, rdminus, dc, theta, pgsz); 

  theta = acos(sqrt(3)/3); 
  fstrtn(pterms, rdsq3, dc, theta, pgsz); 

  theta = acos(-sqrt(3)/3); 
  fstrtn(pterms, rdmsq3, dc, theta, pgsz); 
}

void fstrtn(int p, double *d, const double *sqc, double theta, int pgsz) {
  const double precision = 1.0e-19;
  const double ww = sqrt(2)/2;
  double ctheta = cos(theta);
  ctheta = (fabs(ctheta) <= precision ? 0.0 : ctheta);
  double stheta = sin(-theta);
  stheta = (fabs(stheta) <= precision ? 0.0 : stheta);
  double hsthta = ww*stheta;
  double cthtap = ww*(1.0+ctheta);
  double cthtan = -ww*(1.0-ctheta);

  int ij, im, imp; 
  d[p * pgsz] = 1.0;

  for (ij = 1; ij <= p; ij++) {
    for (im = -ij; im <= -1; im++) {
      int index = ij + (im + p)*pgsz; 
      d[index] = -sqc[ij - im + 2 * (1 + 2 * p)] * d[ij-1 + (im + 1 + p) * pgsz];
      if (im > 1 - ij) 
	      d[index] += sqc[ij + im + 2 * (1 + 2 * p)] * 
          d[ij - 1 + (im - 1 + p) * pgsz];
      d[index] *= hsthta;

      if (im > -ij) 
	      d[index] += d[ij - 1 + (im + p) * pgsz] * ctheta *
	        sqc[ij + im + 2 * p + 1] * sqc[ij - im + 2 * p + 1];      
      d[index] /= ij;
    }

    d[ij + p * pgsz] = d[ij - 1 + p * pgsz] * ctheta;

    if (ij > 1) 
      d[ij + p * pgsz] += hsthta * sqc[ij + 2 * (1 + 2 * p)] *
	     (d[ij - 1 + (-1 + p) * pgsz] + d[ij - 1 + (1 + p) * pgsz]) / ij;
    
    for (im = 1; im <= ij; im++) {
      int index = ij + (im + p) * pgsz; 
      d[index] = -sqc[ij + im + 2 * (1 + 2 * p)] * 
        d[ij - 1 + (im - 1 + p) *pgsz];
      if (im < ij-1) 
	      d[index] += sqc[ij - im + 2 * (1 + 2 * p)] *
          d[ij - 1 + (im + 1 + p) * pgsz];
      d[index] *= hsthta;

      if (im < ij) 
	      d[index] += d[ij- 1 + (im + p) * pgsz] * ctheta *
	        sqc[ij + im + 2 * p + 1] * sqc[ij - im + 2 * p + 1];      
      d[index] /= ij;
    }

    for (imp = 1; imp <= ij; imp++) {
      for (im = -ij; im <= -1; im++) {
        int index1 = ij + imp * (p + 1) + (im + p) * pgsz; 
	      int index2 = ij - 1 + (imp - 1) * (p + 1) + (im + p) * pgsz; 
	      d[index1] = d[index2 + pgsz] * cthtan * sqc[ij - im + 2 * (2 * p + 1)]; 
        if (im > 1 - ij) 
          d[index1] -= d[index2 - pgsz] * cthtap * sqc[ij + im + 4 * p + 2]; 

      if (im > -ij) 
	      d[index1] += d[index2] * stheta * sqc[ij + im + 2 * p + 1] *
          sqc[ij - im + 2 * p + 1];
      d[index1] *= ww / sqc[ij + imp + 2 * (2 * p + 1)];
      }      

      int index3 = ij + imp * (p + 1) + p * pgsz; 
      int index4 = ij - 1 + (imp - 1) * (p + 1) + p * pgsz; 
      d[index3] = ij * stheta * d[index4];
      if (ij > 1) 
        d[index3] -= sqc[ij + 2 * (2 * p + 1)] *
          (d[index4 - pgsz] * cthtap + d[index4 + pgsz] * cthtan);
      d[index3] *= ww / sqc[ij + imp + 2 * (2 * p + 1)]; 

      for (im = 1; im <= ij; im++) {
        int index5 = ij + imp * (p + 1) + (im + p) * pgsz; 
        int index6 = ij - 1 + (imp - 1) * (p + 1) + (im + p) * pgsz; 
        d[index5] = d[index6 - pgsz] * cthtap * sqc[ij + im + 2 * (2 * p + 1)]; 
        if (im < ij - 1) 
          d[index5] -= d[index6 + pgsz] * cthtan * sqc[ij - im + 4 * p + 2]; 

        if (im < ij) 
          d[index5] += d[index6] * stheta * sqc[ij + im + 2 * p + 1] *
            sqc[ij - im + 2 * p + 1];
        d[index5] *= ww/sqc[ij + imp + 2 * (2 * p + 1)];
      }
    }
  }
}

void numthetafour(fmm_param_t *fmm_param) {
  int *numphys = fmm_param->numphys; 
  int nlambs   = fmm_param->nlambs; 

  if (nlambs == 9) {
    numphys[0] = 4;
    numphys[1] = 8;
    numphys[2] = 12;
    numphys[3] = 16;
    numphys[4] = 20;
    numphys[5] = 20;
    numphys[6] = 24;
    numphys[7] = 8;
    numphys[8] = 2;
  } else if (nlambs == 18) {
    numphys[0]  = 6;
    numphys[1]  = 8;
    numphys[2]  = 12;
    numphys[3]  = 16;
    numphys[4]  = 20;
    numphys[5]  = 26;
    numphys[6]  = 30;
    numphys[7]  = 34;
    numphys[8]  = 38;
    numphys[9]  = 44;
    numphys[10] = 48;
    numphys[11] = 52;
    numphys[12] = 56;
    numphys[13] = 60;
    numphys[14] = 60;
    numphys[15] = 52;
    numphys[16] = 4;
    numphys[17] = 2;
  }
}


void rlscini(fmm_param_t *fmm_param) {
  int pterms    = fmm_param->pterms; 
  int nlambs    = fmm_param->nlambs; 
  int pgsz      = fmm_param->pgsz; 
  double *rlsc  = fmm_param->rlsc; 
  double *rlams = fmm_param->rlams; 

  double *factorial = CALLOC(2 * pterms + 1, sizeof(double));
  double *rlampow = CALLOC(pterms + 1, sizeof(double));
  assert(factorial != NULL);
  assert(rlampow != NULL);

  factorial[0] = 1;
  for (int i = 1; i <= 2 * pterms; i++)
    factorial[i] = factorial[i-1]*sqrt(i);
 
  for (int nell = 0; nell < nlambs; nell++) {
    double rmul = rlams[nell];
    rlampow[0] = 1;
    for (int j = 1;  j <= pterms; j++)
      rlampow[j] = rlampow[j - 1] * rmul;      
    for (int j = 0; j <= pterms; j++) {
      for (int k = 0; k <= j; k++) {
        rlsc[j + k * (pterms + 1) + nell * pgsz] = rlampow[j] /
          factorial[j - k] / factorial[j + k];
      }
    }    
  }
    
  FREE(factorial);
  FREE(rlampow);
}

void mkfexp(fmm_param_t *fmm_param) {
  int nlambs = fmm_param->nlambs; 
  int *numphys = fmm_param->numphys; 
  int *numfour = fmm_param->numfour; 
  double complex *fexpe = fmm_param->fexpe; 
  double complex *fexpo = fmm_param->fexpo; 
  double complex *fexpback = fmm_param->fexpback; 

  int nexte = 0; 
  int nexto = 0; 
  double m_pi = acos(-1); 

  for (int i = 0; i < nlambs; i++) {
    int nalpha = numphys[i]; 
    int nalpha2 = nalpha / 2; 
    double halpha = 2.0 * m_pi / nalpha; 
    for (int j = 1; j <= nalpha2; j++) {
      double alpha = (j - 1) * halpha; 
      for (int nm = 2; nm <= numfour[i]; nm += 2) {
        fexpe[nexte] = cexp((nm - 1) * _Complex_I * alpha);
        nexte++;
      }

      for (int nm = 3; nm <= numfour[i]; nm += 2) {
        fexpo[nexto] = cexp((nm - 1) * _Complex_I * alpha);
        nexto++;
      }
    }
  }

  int next = 0; 
  for (int i = 0; i < nlambs; i++) {
    int nalpha = numphys[i]; 
    int nalpha2 = nalpha / 2; 
    double halpha = 2.0 * m_pi / nalpha; 
    for (int nm = 3; nm <= numfour[i]; nm += 2) {
      for (int j = 1; j <= nalpha2; j++) {
        double alpha = (j - 1) * halpha; 
        fexpback[next] = cexp(-(nm - 1) * _Complex_I * alpha);
        next++;
      }
    }

    for (int nm = 2; nm <= numfour[i]; nm += 2) {
      for (int j = 1; j <= nalpha2; j++) {
        double alpha = (j - 1) * halpha; 
        fexpback[next] = cexp(-(nm - 1) * _Complex_I * alpha);
        next++;
      }
    }
  }
}

void mkexps(fmm_param_t *fmm_param) {
  int nlambs = fmm_param->nlambs; 
  int *numphys = fmm_param->numphys;
  double *rlams = fmm_param->rlams; 
  double complex *xs = fmm_param->xs; 
  double complex *ys = fmm_param->ys; 
  double *zs = fmm_param->zs; 
  int ntot = 0; 
  double m_pi = acos(-1); 
  for (int nell = 0; nell < nlambs; nell++) {
    double hu = 2.0 * m_pi / numphys[nell]; 
    for (int mth = 0; mth < numphys[nell] / 2; mth++) {
      double u = mth * hu; 
      int ncurrent = 3 * (ntot + mth);
      zs[ncurrent]     = exp(-rlams[nell]);
      zs[ncurrent + 1] = zs[ncurrent] * zs[ncurrent]; 
      zs[ncurrent + 2] = zs[ncurrent] * zs[ncurrent + 1];
      xs[ncurrent]     = cexp(_Complex_I * cos(u) * rlams[nell]);
      xs[ncurrent + 1] = xs[ncurrent] * xs[ncurrent];
      xs[ncurrent + 2] = xs[ncurrent + 1] * xs[ncurrent]; 
      ys[ncurrent]     = cexp(_Complex_I * sin(u) * rlams[nell]);
      ys[ncurrent + 1] = ys[ncurrent] * ys[ncurrent]; 
      ys[ncurrent + 2] = ys[ncurrent + 1] * ys[ncurrent]; 
    }
    ntot += numphys[nell]/2; 
  }
}

#elif YUKAWA

fmm_param_t *construct_param(const fmm_config_t *fmm_config, 
                             const fmm_dag_t *fmm_dag) {
  fmm_param_t *param = CALLOC(1, sizeof(fmm_param_t)); 
  assert(param != NULL); 
  param->beta = fmm_config->beta; 

  int accuracy = fmm_config->accuracy; 
  if (accuracy == 3) {
    param->pterms = 9;
    param->nlambs = 9;
    param->pgsz = 100;
    param->dcpgsz = 6859;
  } else if (accuracy == 6) {
    param->pterms = 18;
    param->nlambs = 18;
    param->pgsz = 361;
    param->dcpgsz = 50653;
  }

  int pterms = param->pterms;
  int nlambs = param->nlambs;
  int pgsz   = param->pgsz; 
  int dcpgsz = param->dcpgsz; 
  int nlevel = (fmm_dag->nslev >= fmm_dag->ntlev ? 
                fmm_dag->nslev : fmm_dag->ntlev);

  param->numfour = CALLOC(nlambs, sizeof(int)); 
  param->numphys = CALLOC(nlambs * (nlevel + 1), sizeof(int)); 
  param->whts = CALLOC(nlambs, sizeof(double));
  param->rlams = CALLOC(nlambs, sizeof(double));
  param->ytop = CALLOC(pgsz, sizeof(double));
  param->rdplus = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->rdminus = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->rdsq3 = CALLOC(pgsz * (2 * pterms + 1), sizeof(double));
  param->rdmsq3 = CALLOC(pgsz * (2 * pterms + 1), sizeof(double)); 

  yhfrmini(param); 
  yhrotgen(param); 
  vwts(param); 
  numthetahalf(param); 
  numthetafour(param, 0); 

  for (int i = 0; i < nlambs; i++) {
    double test1 = param->rlams[i];
    double test2 = sqrt(test1 * test1 + 
                        2.0 * test1 * param->beta * fmm_dag->size);
    int indd = i;
    int mmax = param->numphys[i];
    for (int j = i; j < nlambs; j++) {
      if (test2 <= param->rlams[j]) {
        indd = j;
        break;
      } else {
        mmax = (param->numphys[j] >= mmax ? param->numphys[j] : mmax); 
      }
    }
    param->numphys[i] = (param->numphys[indd] >= mmax ? 
                         param->numphys[indd] : mmax); 
  }
  
  param->nexptot = 0;
  param->nthmax = 0;
  param->mnexptotp = 0;

  for (int i = 0; i < nlambs; i++) {
    param->nexptot += param->numfour[i];
    if (param->numfour[i] > param->nthmax)
      param->nthmax = param->numfour[i];
    param->mnexptotp += param->numphys[i];
  }

  param->mnexptotp /= 2.0;
  param->nexpmax = (param->nexptot >= param->mnexptotp ? 
                    param->nexptot : param->mnexptotp) + 1; 

  param->fexpe = CALLOC(15000 * (nlevel + 1), sizeof(double complex));
  param->fexpo = CALLOC(15000 * (nlevel + 1), sizeof(double complex));
  param->fexpback = CALLOC(15000 * (nlevel + 1), sizeof(double complex));
  param->vnexptot = CALLOC(nlevel + 1, sizeof(int));
  param->vnexptotp = CALLOC(nlevel + 1, sizeof(int));
  param->vnthmax = CALLOC(nlevel + 1, sizeof(int));
  param->dcu = CALLOC(dcpgsz * (nlevel + 1), sizeof(double));
  param->dcd = CALLOC(dcpgsz * (nlevel + 1), sizeof(double));
  param->sfactor = CALLOC(nlevel + 1, sizeof(double));
  //param->sfactor2 = CALLOC(nlevel + 1, sizeof(double));
  param->betascale = CALLOC(nlevel + 1, sizeof(double));
  param->zs = CALLOC(param->nexpmax * 3 * (nlevel + 1), sizeof(double));
  param->rlsc = CALLOC(pgsz * nlambs * (nlevel + 1), sizeof(double));
  param->xs = CALLOC(param->nexpmax * 3 * (nlevel + 1), sizeof(double complex));
  param->ys = CALLOC(param->nexpmax * 3 * (nlevel + 1), sizeof(double complex));

  param->scale = (param->beta * fmm_dag->size > 1.0 ? 
                  1.0 / fmm_dag->size : param->beta);

  for (int lev = 0; lev <= nlevel; lev++) {
    double size = fmm_dag->size / (1 << lev); 
    double r0 = sqrt(3.0) * size / 4;
    param->betascale[lev] = param->beta * size;
    param->sfactor[lev] = param->scale * size; 

    ympshftcoef(param, r0, lev); // dcu at every level
    ylcshftcoef(param, r0, lev); // dcd at every level

    // numphys at each level 
    numthetafour(param, lev); 
    for (int i = 0; i < nlambs; i++) {
      double test1 = param->rlams[i]; 
      double test2 = sqrt(test1 * test1 + test1 * param->sfactor[lev]); 
      int indd = i; 
      int mmax = param->numphys[lev * nlambs + i]; 
      for (int j = i; j < nlambs; j++) {
        if (test2 <= param->rlams[j]) {
          indd = j;
          break;
        } else {
          mmax = (mmax > param->numphys[lev * nlambs + j] ? 
                  mmax : param->numphys[lev * nlambs + j]); 
        }
      }
      param->numphys[lev * nlambs + i] = 
        (mmax > param->numphys[lev * nlambs + indd] ? 
         mmax : param->numphys[lev * nlambs + indd]); 
    }

    param->vnexptot[lev] = 0;
    param->vnthmax[lev] = 0;
    param->vnexptotp[lev] = 0; 

    for (int i = 0; i < nlambs; i++) {
      param->vnexptot[lev] += param->numfour[i];
      if (param->numfour[i] > param->vnthmax[lev])
        param->vnthmax[lev] = param->numfour[i];
      param->vnexptotp[lev] += param->numphys[lev * nlambs + i];
    }

    param->vnexptotp[lev] /= 2.0;

    ymkfexp(param, lev); // fexpe, fexpo, fexpback at every level
    yrlscini(param, lev); // rlsc at every level
    ymkexps(param, lev); // xs, ys, zs at every level
  } 

  return param;
}

void destruct_param(fmm_param_t *fmm_param) {
  FREE(fmm_param->numfour);
  FREE(fmm_param->numphys);
  FREE(fmm_param->whts);
  FREE(fmm_param->rlams);
  FREE(fmm_param->ytop);
  FREE(fmm_param->rdplus);
  FREE(fmm_param->rdminus);
  FREE(fmm_param->rdsq3);
  FREE(fmm_param->rdmsq3);
  FREE(fmm_param->fexpe);
  FREE(fmm_param->fexpo);
  FREE(fmm_param->fexpback);
  FREE(fmm_param->vnexptot);
  FREE(fmm_param->vnexptotp);
  FREE(fmm_param->vnthmax);
  FREE(fmm_param->dcu);
  FREE(fmm_param->dcd);
  FREE(fmm_param->sfactor);
  //FREE(fmm_param->sfactor2);
  FREE(fmm_param->betascale);
  FREE(fmm_param->zs);
  FREE(fmm_param->rlsc);
  FREE(fmm_param->xs);
  FREE(fmm_param->ys);
  FREE(fmm_param); 
}

void yhfrmini(fmm_param_t *fmm_param) {
  int pterms   = fmm_param->pterms; 
  double *ytop = fmm_param->ytop; 

  double *factorial = CALLOC(2 * pterms + 2, sizeof(double));
  assert(factorial != NULL); 

  factorial[0] = 1;
  for (int ell = 1; ell <= 2 * pterms + 1; ell++)
    factorial[ell] = factorial[ell - 1] * ell;

  for (int m = 0; m <= pterms; m++) {
    int offset = m * (pterms + 1);
    for (int ell = m; ell <= pterms; ell++) {
      ytop[ell + offset] = factorial[ell - m] / 
        factorial[ell + m] * (2 * ell + 1);
    }
  }

  FREE(factorial); 
}

void yhrotgen(fmm_param_t *fmm_param) {
  int pterms      = fmm_param->pterms; 
  int pgsz        = fmm_param->pgsz; 
  double *rdplus  = fmm_param->rdplus;
  double *rdminus = fmm_param->rdminus;
  double *rdsq3   = fmm_param->rdsq3;
  double *rdmsq3  = fmm_param->rdmsq3; 
  
  double *carray = CALLOC((4 * pterms + 1) * (4 * pterms + 1), sizeof(double)); 
  assert(carray != NULL); 
  bnlcft(carray, 4 * pterms);

  double theta = acos(0); 

  yhfstrtn(pterms, theta, carray, rdplus, pgsz);

  theta = -theta; 
  yhfstrtn(pterms, theta, carray, rdminus, pgsz);

  theta = acos(sqrt(3) / 3);
  yhfstrtn(pterms, theta, carray, rdsq3, pgsz);

  theta = acos(-sqrt(3) / 3);
  yhfstrtn(pterms, theta, carray, rdmsq3, pgsz);

  FREE(carray); 
}

void yhfstrtn(int p, double theta, const double *sqc, double *d, int pgsz) {
  const double precision = 1.0e-19;
  const double ww = sqrt(2) / 2;
  double ctheta = cos(theta);
  ctheta = (fabs(ctheta) <= precision ? 0.0 : ctheta);
  double stheta = sin(-theta);
  stheta = (fabs(stheta) <= precision ? 0.0 : stheta);
  double hsthta = ww * stheta;
  double cthtap = ww * (1 + ctheta);
  double cthtan = -ww * (1 - ctheta);
  d[p * pgsz] = 1;

  for (int ij = 1; ij <= p; ij++) {
    for (int im = -ij; im <= -1; im++) {
      int index = ij + (im + p) * pgsz;
      d[index] = -sqc[ij - im + 2 * (1 + 4 * p)] * 
        d[ij - 1 + (im + 1 + p) * pgsz];
      if (im > 1 - ij)
        d[index] += sqc[ij + im + 2 * (1 + 4 * p)] * 
          d[ij - 1 + (im - 1 + p) * pgsz];
      d[index] *= hsthta;
      if (im > -ij)
        d[index] += ctheta * sqc[ij + im + 4 * p + 1] * 
          sqc[ij - im + 4 * p + 1] * d[index - 1];
      d[index] /= ij;
    }

    d[ij + p * pgsz] = d[ij - 1 + p * pgsz] * ctheta;

    if (ij > 1)
      d[ij + p * pgsz] += hsthta * sqc[ij + 2 * (1 + 4 * p)] *
        (d[ij - 1 + (-1 + p) * pgsz] + d[ij - 1 + (1 + p) * pgsz]) / ij;

    for (int im = 1; im <= ij; im++) {
      int index = ij + (im + p) * pgsz;
      d[index] -= sqc[ij + im + 2 * (1 + 4 * p)] * 
        d[ij - 1 + (im - 1 + p) * pgsz];
      if (im < ij - 1)
        d[index] += sqc[ij - im + 2 * (1 + 4 * p)] * 
          d[ij -1 + (im + 1 + p) * pgsz];
      d[index] *= hsthta;
      if (im < ij)
        d[index] += ctheta * sqc[ij + im + 4 * p + 1] * 
          sqc[ij - im + 4 * p + 1] * d[index - 1];
      d[index] /= ij;
    }

    for (int imp = 1; imp <= ij; imp++) {
      for (int im = -ij; im <= -1; im++) {
        int index1 = ij + imp * (p + 1) + (im + p) * pgsz;
        int index2 = ij - 1 + (imp - 1) * (p + 1) + (im + p) * pgsz;
        d[index1] = cthtan * sqc[ij - im + 2 * (4 * p + 1)] * d[index2 + pgsz];
        if (im > 1 - ij)
          d[index1] -= d[index2 - pgsz] * cthtap * 
            sqc[ij + im + 2 * (4 * p + 1)];
        if (im > -ij)
          d[index1] += d[index2] * stheta * sqc[ij + im + 4 * p + 1] * 
            sqc[ij - im + 4 * p + 1];
        d[index1] *= ww / sqc[ij + imp + 2 * (4 * p + 1)];
      }

      int index3 = ij + imp * (p + 1) + p * pgsz;
      int index4 = ij - 1 + (imp - 1) * (p + 1) + p * pgsz;
      d[index3] = ij * stheta * d[index4];
      if (ij > 1)
        d[index3] -= sqc[ij + 2 * (4 * p + 1)] *
          (d[index4 - pgsz] * cthtap + d[index4 + pgsz] * cthtan);
      d[index3] *= ww / sqc[ij + imp + 2 * (4 * p + 1)];
      
      for (int im = 1; im <= ij; im++) {
        int index5 = ij + imp * (p + 1) + (im + p) * pgsz;
        int index6 = ij - 1 + (imp - 1) * (p + 1) + (im + p) * pgsz;
        d[index5] = d[index6 - pgsz] * cthtap * sqc[ij + im + 2 * (4 * p + 1)];
        if (im < ij - 1)
          d[index5] -= d[index6 + pgsz] * cthtan * 
            sqc[ij - im + 2 * (4 * p + 1)];
        if (im < ij)
          d[index5] += d[index6] * stheta * sqc[ij + im + 4 * p + 1] *
            sqc[ij - im + 4 * p + 1];
        d[index5] *= ww / sqc[ij + imp + 2 * (4 * p + 1)];
      }
    }
  }
  
  double *fac = CALLOC(2 * p + 1, sizeof(double));
  assert(fac != NULL); 

  fac[0] = 1;
  for (int ij = 1; ij <= 2 * p; ij++)
    fac[ij] = fac[ij - 1] * ij;

  for (int ij = 0; ij <= p; ij++) {
    for (int im = 0; im <= ij; im++) {
      for (int imp = -ij; imp <= ij; imp++) {
        int impabs = fabs(imp);
        d[ij + im * (p + 1) + (imp + p) * pgsz] *=
          sqrt(fac[ij + im] / fac[ij + impabs] * 
               fac[ij - impabs] / fac[ij - im]);
      }
    }
  }
  
  FREE(fac);
}

void numthetafour(fmm_param_t *fmm_param, int level) {
  int nlambs   = fmm_param->nlambs; 
  int *numphys = &(fmm_param->numphys[level * nlambs]); 

  if (nlambs == 9) {
    numphys[0] = 4;
    numphys[1] = 8;
    numphys[2] = 12;
    numphys[3] = 16;
    numphys[4] = 20;
    numphys[5] = 20;
    numphys[6] = 24;
    numphys[7] = 8;
    numphys[8] = 2;
  } else if (nlambs == 18) {
    numphys[0]  = 6;
    numphys[1]  = 8;
    numphys[2]  = 12;
    numphys[3]  = 16;
    numphys[4]  = 20;
    numphys[5]  = 26;
    numphys[6]  = 30;
    numphys[7]  = 34;
    numphys[8]  = 38;
    numphys[9]  = 44;
    numphys[10] = 48;
    numphys[11] = 52;
    numphys[12] = 56;
    numphys[13] = 60;
    numphys[14] = 60;
    numphys[15] = 52;
    numphys[16] = 4;
    numphys[17] = 2;
  }
}

void ympshftcoef(fmm_param_t *fmm_param, double r0, int level) {
  int pgsz     = fmm_param->pgsz; 
  int dcpgsz   = fmm_param->dcpgsz; 
  int pterms   = fmm_param->pterms; 
  double scale = fmm_param->sfactor[level] / 2;
  double *c    = &fmm_param->dcu[dcpgsz*level];
  double beta  = fmm_param->beta; 

  double *fac = CALLOC(2 * pterms + 2, sizeof(double));
  double *bj = CALLOC(2 * pterms + 3, sizeof(double));
  assert(fac != NULL);
  assert(bj != NULL); 

  fac[0] = 1;
  for (int i = 1; i <= 2 * pterms + 1; i++)
    fac[i] = fac[i - 1] * i;

  double r0k = r0 * beta;
  int ncalc;
  in(r0k, r0k, 2 * pterms + 2, bj, &ncalc);

  for (int mnew = 0; mnew <= pterms; mnew++) {
    for (int ellnew = mnew; ellnew <= pterms; ellnew++ ) {
      int offset = ellnew * (pterms + 1);
      for (int nn = mnew; nn <= pterms; nn++) {
        c[mnew + offset + nn * pgsz] = 0;
        int np1 = (nn < ellnew ? nn : ellnew); 
        for (int np = mnew; np <= np1; np++) {
          c[mnew + offset + nn * pgsz] += pow(scale, nn - ellnew) *
            pow(2, -ellnew - np) * pow(-1, ellnew + nn) * (2 * ellnew + 1) *
            fac[ellnew - mnew] / fac[np + mnew] * fac[nn + mnew] * 
            fac[2 * np] / fac[np] / fac[np - mnew] / fac[ellnew - np] / 
            fac[nn - np] *bj[ellnew + nn - np] * pow(r0k, ellnew + nn - 2 * np);
        }
      }
    }
  }
  FREE(fac);
  FREE(bj);
}

void ylcshftcoef(fmm_param_t *fmm_param, double r0, int level) {
  int pgsz     = fmm_param->pgsz; 
  int pterms   = fmm_param->pterms; 
  int dcpgsz   = fmm_param->dcpgsz; 
  double scale = fmm_param->sfactor[level];
  double *c    = &fmm_param->dcd[dcpgsz*level];
  double beta  = fmm_param->beta; 

  double *fac = CALLOC(2 * pterms + 2, sizeof(double));
  double *bj = CALLOC(2 * pterms + 3, sizeof(double));
  assert(fac != NULL);
  assert(bj != NULL); 

  fac[0] = 1;
  for (int i = 1; i <= 2 * pterms + 1; i++)
    fac[i] = fac[i - 1] * i;

  double r0k = r0 * beta;
  int ncalc;
  in(r0k, r0k, 2 * pterms + 2, bj, &ncalc);

  for (int mnew = 0; mnew <= pterms; mnew++) {
    for (int ellnew = mnew; ellnew <= pterms; ellnew++) {
      int offset = ellnew * (pterms + 1);
      for (int nn = mnew; nn <= pterms; nn++) {
        c[mnew + offset + nn * pgsz] = 0;
        int np1 = (nn < ellnew ? nn : ellnew); 
        for (int np = mnew; np <= np1; np++) {
          c[mnew + offset + nn * pgsz] += pow(scale, ellnew - nn) *
            pow(2, -ellnew - np) * (2 * ellnew + 1) * fac[ellnew - mnew] /
            fac[np + mnew] * fac[nn + mnew] * fac[2 * np] / fac[np] / 
            fac[np - mnew] / fac[ellnew - np] / fac[nn - np] * 
            bj[ellnew + nn - np] * pow(r0k, ellnew + nn - 2 * np);
        }
      }
    }
  }
  FREE(fac);
  FREE(bj);
}

void in(double scal, double x, int nb, double *b, int *ncalc) {
  const double ensig = 1.0e-4;
  const double enmten = 1.0e-300;
  double m_pi_2 = acos(0); 

  for (int i = 0; i <= nb; i++)
    b[i] = 0;

  if (x < 0)
    exit(-1);

  if (x <= ensig){
    double xscal = x / scal;
    double term1 = 1.0;
    double term2 = 0.5 * x * x;
    b[0] = term1 * (1.0 + term2 / 3.0);
    for (int i = 1; i <= nb; i++) {
      term1 = term1 * xscal / (2 * i + 1);
      term1 = (term1 <= enmten ? 0.0 : term1);
      b[i] = term1 * (1.0 + term2 / (2 * i + 3));
    }
    *ncalc = nb + 1;
  } else if (x > 1.0e2) {
    for (int i = 0; i <= nb; i++)
      b[i] = 0;
    *ncalc = nb + 1;
  } else {
    double constant = sqrt(m_pi_2 / x);
    double alpha = 0.5;
    int ize = 1;
    int nb1 = nb + 1;
    ribesl_(&x, &alpha, &nb1, &ize, b, ncalc);
    for (int i = 0; i <= nb; i++) {
      b[i] *= constant;
      constant /= scal; 
      constant = (fabs(b[i]) <= enmten ? 0 : constant);
    }
  }
}

void ymkfexp(fmm_param_t *fmm_param, int level) {
  int nlambs = fmm_param->nlambs; 
  int *numfour = fmm_param->numfour; 
  int *numphys = &fmm_param->numphys[nlambs * level]; 
  double complex *fexpe = &fmm_param->fexpe[15000 * level]; 
  double complex *fexpo = &fmm_param->fexpo[15000 * level];
  double complex *fexpback = &fmm_param->fexpback[15000 * level]; 

  int nexte = 0;
  int nexto = 0;
  double pi = acos(-1);

  for (int i = 0; i < nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha / 2;
    double halpha = 2.0 * pi / nalpha;
    for (int j = 1; j <= nalpha2; j++) {
      double alpha = (j - 1) * halpha;
      for (int nm = 2; nm <= numfour[i]; nm += 2) {
        fexpe[nexte] = cexp((nm - 1) * _Complex_I * alpha);
        nexte++;
      }
      
      for (int nm = 3; nm <= numfour[i]; nm += 2) {
        fexpo[nexto] = cexp((nm - 1) * _Complex_I * alpha);
        nexto++;
      }
    }
  }

  int next = 0;
  for (int i = 0; i < nlambs; i++) {
    int nalpha = numphys[i];
    int nalpha2 = nalpha / 2;
    double halpha = 2 * pi / nalpha;
    for (int nm = 3; nm <= numfour[i]; nm += 2) {
      for (int j = 1; j <= nalpha2; j++) {
        double alpha = (j - 1) * halpha;
        fexpback[next] = cexp(-(nm - 1) * _Complex_I * alpha);
        next++;
      }
    }

    for (int nm = 2; nm <= numfour[i]; nm += 2) {
      for (int j = 1; j <= nalpha2; j++) {
        double alpha = (j - 1) * halpha;
        fexpback[next] = cexp(-(nm - 1) * _Complex_I * alpha);
        next++;
      }
    }
  } 
}

void yrlscini(fmm_param_t *fmm_param, int level) {
  int pgsz = fmm_param->pgsz; 
  int nlambs = fmm_param->nlambs; 
  int pterms = fmm_param->pterms; 
  double scale = fmm_param->sfactor[level];
  double beta = fmm_param->betascale[level];
  double *rlsc = &fmm_param->rlsc[level * pgsz * nlambs];

  for (int nell=1; nell <= nlambs; nell++) {
    double u1 = fmm_param->rlams[nell - 1]/beta + 1.0;
    lgndrgt1(scale, pterms, u1, &rlsc[pgsz * (nell - 1)]);
  }
}

void lgndrgt1(double scale, int nmax, double x, double *y) {
  double u = sqrt(x * x - 1) * scale;
  double v = scale * x;
  double w = scale * scale;
  y[0] = 1;
  y[1] = y[0] * v;
  for (int n = 2; n <= nmax; n++)
    y[n] = ((2 * n - 1) * v * y[n - 1] - (n - 1) * w * y[n - 2]) / n;

  for (int m = 1; m < nmax; m++) {
    int offset = m * (nmax + 1);
    y[m + offset] = y[m - 1 + (m - 1) * (nmax + 1)] * u * (2 * m - 1);
    y[m + 1 + offset] = y[m + offset] * (2 * m + 1) * v;
    for (int n = m + 2; n <= nmax; n++)
      y[n + offset] = ((2 * n - 1) * v * y[n - 1 + offset] - 
                       (n + m - 1) * w * y[n - 2 + offset] ) / (n - m);
  }
  y[nmax + nmax * (nmax + 1)] = y[nmax - 1 + (nmax - 1) * (nmax + 1)] * 
    u * (2 * nmax - 1);
}

void ymkexps(fmm_param_t *fmm_param, int level) {
  int nlambs = fmm_param->nlambs; 
  int nexpmax = fmm_param->nexpmax; 
  double beta = fmm_param->betascale[level];
  int *numphys = &fmm_param->numphys[nlambs * level];
  double complex *xs = &fmm_param->xs[nexpmax * 3 * level];
  double complex *ys = &fmm_param->ys[nexpmax * 3 * level];
  double *zs = &fmm_param->zs[nexpmax * 3 * level];

  int ntot = 0;
  double pi = acos(-1); 
  for (int nell = 0; nell < nlambs; nell++) {
    double w1 = fmm_param->rlams[nell] + beta;
    double w2 = sqrt(fmm_param->rlams[nell] * 
                     (fmm_param->rlams[nell] + beta * 2));
    double hu = 2 * pi / numphys[nell];
    for (int mth = 0; mth < numphys[nell]/2; mth++) {
      double u = mth * hu;
      int ncurrent = 3 * (ntot + mth);
      zs[ncurrent]      = exp(-w1);
      zs[ncurrent + 1] = zs[ncurrent] * zs[ncurrent];
      zs[ncurrent + 2] = zs[ncurrent + 1] * zs[ncurrent];
      xs[ncurrent]     = cexp(w2 * cos(u) * _Complex_I);
      xs[ncurrent + 1] = xs[ncurrent] * xs[ncurrent];
      xs[ncurrent + 2] = xs[ncurrent + 1] * xs[ncurrent];
      ys[ncurrent]     = cexp(w2 * sin(u) * _Complex_I);
      ys[ncurrent + 1] = ys[ncurrent] * ys[ncurrent];
      ys[ncurrent + 2] = ys[ncurrent + 1] * ys[ncurrent];
    }
    ntot += numphys[nell] / 2;
  }
}

void kn(double scal, double x, int nb, double *by, int *ncalc) {
  // translated by f2c (version 20090411) and modified by Nikos Pitsianis
  
  /* Initialized data */
  
  static const double xmin = 4.46e-308;
  static const double xinf = 1.79e308;
  static const double xlarge = 1e8;
  
  /* System generated locals */
  int n;
  double d__1;
  
  /* Builtin functions */
  double atan(double), exp(double);
  
  /* Local variables */
  int i;
  double p, u1, u2, ex;
  static const double zero = 0.0;
  static const double halfpi = 1.57079632679; // atan(1.) * 2.;

  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd */
  
  /* purpose: */
  /*  this routine calculates scaled sphecial bessel functions y_n (x) */
  /*  for non-negative argument x, and n=0,...,nb. */
  
  /* on input: */
  /*   scal: the scaling factor. */
  /*   x: working precision non-negative real argument for which */
  /*      y's are to be calculated. */
  /*   nb: integer number of functions to be calculated, nb .gt. 0. */
  /*       the first function calculated is of order 0, and the */
  /*       last is of order nb. */
  /*   by: working precision output vector of length nb+1.  if the */
  /*       routine terminates normally (ncalc=nb+1), the vector by */
  /*       contains the functions y(0,x), ... , y(nb,x), */
  /*       if (0 .lt. ncalc .lt. nb), by(i) contains correct function */
  /*       values for i .le. ncalc, and contains the ratios */
  /*       y(i-1,x)/y(i-2,x) for the rest of the array. */
  /*   ncalc: integer output variable indicating possible errors. */
  /*       before using the vector by, the user should check that */
  /*       ncalc=nb, i.e., all orders have been calculated to */
  /*       the desired accuracy.  see error returns below. */
  
  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd */
  
  /* explanation of machine-dependent constants */
  
  /*   beta   = radix for the floating-point system */
  /*   p      = number of significant base-beta digits in the */
  /*            significand of a floating-point number */
  /*   minexp = smallest representable power of beta */
  /*   maxexp = smallest power of beta that overflows */
  /*   eps    = beta ** (-p) */
  /*   del    = machine number below which sin(x)/x = 1; approximately */
  /*            sqrt(eps). */
  /*   xmin   = smallest acceptable argument for rbesy; approximately */
  /*            max(2*beta**minexp,2/xinf), rounded up */
  /*   xinf   = largest positive machine number; approximately */
  /*            beta**maxexp */
  /*   thresh = lower bound for use of the asymptotic form; approximately */
  /*            aint(-log10(eps/2.0))+1.0 */
  /*   xlarge = upper bound on x; approximately 1/del, because the sine */
  /*            and cosine functions have lost about half of their */
  /*            precision at that point. */
  
  
  /*     approximate values for some important machines are: */
  
  /*                        beta    p     minexp      maxexp      eps */
  
  /*  cray-1        (s.p.)    2    48     -8193        8191    3.55e-15 */
  /*  cyber 180/185 */
  /*    under nos   (s.p.)    2    48      -975        1070    3.55e-15 */
  /*  ieee (ibm/xt, */
  /*    sun, etc.)  (s.p.)    2    24      -126         128    5.96e-8 */
  /*  ieee (ibm/xt, */
  /*    sun, etc.)  (d.p.)    2    53     -1022        1024    1.11d-16 */
  /*  ibm 3033      (d.p.)   16    14       -65          63    1.39d-17 */
  /*  vax           (s.p.)    2    24      -128         127    5.96e-8 */
  /*  vax d-format  (d.p.)    2    56      -128         127    1.39d-17 */
  /*  vax g-format  (d.p.)    2    53     -1024        1023    1.11d-16 */
  
  
  /*                         del      xmin      xinf     thresh  xlarge */
  
  /* cray-1        (s.p.)  5.0e-8  3.67e-2466 5.45e+2465  15.0e0  2.0e7 */
  /* cyber 180/855 */
  /*   under nos   (s.p.)  5.0e-8  6.28e-294  1.26e+322   15.0e0  2.0e7 */
  /* ieee (ibm/xt, */
  /*   sun, etc.)  (s.p.)  1.0e-4  2.36e-38   3.40e+38     8.0e0  1.0e4 */
  /* ieee (ibm/xt, */
  /*   sun, etc.)  (d.p.)  1.0d-8  4.46d-308  1.79d+308   16.0d0  1.0d8 */
  /* ibm 3033      (d.p.)  1.0d-8  2.77d-76   7.23d+75    17.0d0  1.0d8 */
  /* vax           (s.p.)  1.0e-4  1.18e-38   1.70e+38     8.0e0  1.0e4 */
  /* vax d-format  (d.p.)  1.0d-9  1.18d-38   1.70d+38    17.0d0  1.0d9 */
  /* vax g-format  (d.p.)  1.0d-8  2.23d-308  8.98d+307   16.0d0  1.0d8 */
  
  /* ******************************************************************* */
  
  /* error returns */
  
  /*  in case of an error, ncalc .ne. nb, and not all y's are */
  /*  calculated to the desired accuracy. */
  
  /*  ncalc .lt. -1:  an argument is out of range. for example, */
  /*       nb .le. 0, ize is not 1 or 2, or ize=1 and abs(x) .ge. */
  /*       xmax.  in this case, by(1) = 0.0, the remainder of the */
  /*       by-vector is not calculated, and ncalc is set to */
  /*       min0(nb,0)-2  so that ncalc .ne. nb. */
  /*  ncalc = -1:  y(alpha,x) .ge. xinf.  the requested function */
  /*       values are set to 0.0. */
  /*  1 .lt. ncalc .lt. nb: not all requested function values could */
  /*       be calculated accurately.  by(i) contains correct function */
  /*       values for i .le. ncalc, and and the remaining nb-ncalc */
  /*       array elements contain 0.0. */
  

  /* intrinsic functions required are: */
  
  /*     dble, exp, int, max, min, real, sqrt */
  
  
  /* acknowledgement */
  
  /*  this program draws heavily on temme's algol program for y(a,x) */
  /*  and y(a+1,x) and on campbell's programs for y_nu(x).  temme's */
  /*  scheme is used for  x < thresh, and campbell's scheme is used */
  /*  in the asymptotic region.  segments of code from both sources */
  /*  have been translated into fortran 77, merged, and heavily modified. */
  /*  modifications include parameterization of machine dependencies, */
  /*  use of a new approximation for ln(gamma(x)), and built-in */
  /*  protection against over/underflow. */
  
  /* references: "bessel functions j_nu(x) and y_nu(x) of real */
  /*              order and real argument," campbell, j. b., */
  /*              comp. phy. comm. 18, 1979, pp. 133-142. */
  
  /*             "on the numerical evaluation of the ordinary */
  /*              bessel function of the second kind," temme, */
  /*              n. m., j. comput. phys. 21, 1976, pp. 343-350. */
  
  /*  latest modification: march 19, 1990 */
  
  /*  modified by: w. j. cody */
  /*               applied mathematics division */
  /*               argonne national laboratory */
  /*               argonne, il  60439 */
  
  /* ---------------------------------------------------------------------- */
  /* ---------------------------------------------------------------------- */
  /*  machine-dependent constants */
  /* ---------------------------------------------------------------------- */
  /* ---------------------------------------------------------------------- */
  ex = x;
  if (nb >= 0 && x >= xmin && ex < xlarge) {
    
    /* ---------------------------------------------------------------------- */
    /*  now have first one or two y's */
    /* ---------------------------------------------------------------------- */
    p = exp(-ex) / ex * halfpi;
    by[0] = p;
    by[1] = p * scal * (1. / ex + 1.);
    *ncalc = 1;
    u1 = scal / ex;
    u2 = scal * scal;
    n = nb;
    for (i = 2; i <= n; ++i) {
      if ((d__1 = by[i - 1], abs(d__1)) * u1 >= 
          xinf / (double) ((i << 1) - 1)) {
        break;  // goto L450;
      }
      by[i] = (double) ((i << 1) - 1) * u1 * by[i - 1] + u2 * by[i - 2];
      ++(*ncalc);
    }
    // L450:
    n = nb;
    for (i = *ncalc + 1; i <= n; ++i) {
      by[i] = zero;
    }
  } else {
    by[0] = 0.;
    /* Computing MIN */
    n = nb + 1;
    *ncalc = (n <= 0 ? n : 0) - 1;
  }  
}

#endif 


