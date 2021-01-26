// ---------------------------------------------------------------------------
// Copyright (c) 2014 Bo Zhang, Jingfang Huang, Nikos P. Pitsinis, Xiaobai Sun
//
// This file is part of FMM-LAPLACE/YUKAWA
//
// FMM-LAPLACE/YUKAWA is free software: you can redistribute it and/or modify it
// under the terms of GNU General Public Licenses as published by the Free
// Software Foundation, either version 3 of the licenses, or any later version.
//
// FMM-LAPLACE/YUKAWA is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FINESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details. 
//
// You should have received a copy of the GNU General Public License along with
// FMM-LAPLACE/YUKAWA. If not, see <http://www.gnu.org/license/>.
// ----------------------------------------------------------------------------

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include "fmm.h"

static void usage(FILE *stream) {
  fprintf(stream, "Usage: fmm [options]\n"
	  "\t-n, number of source points\n"
	  "\t-m, number of target points\n"
	  "\t-a, accuracy requirement, only 3 and 6 are valid\n"
	  "\t-d, distribution of the particles\n"
	  "\t    type-1 data is uniformly distributed inside a box\n"
	  "\t    type-2 data is uniformly distributed over a sphere\n"
	  "\t-b, parameter beta (>=0) for Yukawa kernel exp(-beta * r) / r\n"
	  "\t-s, maximum number of points allowed per leaf box\n");
}


int main(int argc, char *argv[]) {
  // default test configuration

  fmm_config_t fmm_config = {
#if YUKAWA
    .beta     = 0.1,
#endif
    .nsources = 10000,
    .ntargets = 10000,
    .datatype = 1,
    .accuracy = 3,
    .s        = 40, 
  };

  // Parse command line options
  int opt = 0; 
  while ((opt = getopt(argc, argv, "n:m:a:d:s:b:h")) != -1) {
    switch (opt) {
    case 'n':
      fmm_config.nsources = atoi(optarg);
      break;
    case 'm':
      fmm_config.ntargets = atoi(optarg);
      break;
    case 'a':
      fmm_config.accuracy = atoi(optarg);
      break;
    case 'd':
      fmm_config.datatype = atoi(optarg);
      break;
    case 's':
      fmm_config.s        = atoi(optarg);
      break;
#if YUKAWA
    case 'b':
      fmm_config.beta     = atof(optarg);
      break;
#endif
    case 'h':
      usage(stdout); 
      return 0;
    case '?':
    default:
      usage(stderr); 
      return -1;
    }
  }

  printf("test configuration:\n"
         "  number of sources: %d\n"
         "  number of targets: %d\n"
         "  accuracy requirement: %d\n",
         fmm_config.nsources, fmm_config.ntargets, fmm_config.s);

  // Generate test data according to the configuration
  Real_t *sources = CALLOC(fmm_config.nsources * 3, sizeof(Real_t)); 
  Real_t *charges = CALLOC(fmm_config.nsources, sizeof(Real_t)); 
  Real_t *targets = CALLOC(fmm_config.ntargets * 3, sizeof(Real_t)); 
  double *potential = CALLOC(fmm_config.ntargets, sizeof(double)); 
  double *field = CALLOC(fmm_config.ntargets * 3, sizeof(double)); 

  double pi = acos(-1); 
  if (fmm_config.datatype == 1) {
    for (int i = 0; i < fmm_config.nsources; i++) {
      int j = 3 * i;
      charges[i] = 1.0 * rand() / RAND_MAX - 0.5;
      sources[j] = 1.0 * rand() / RAND_MAX - 0.5;
      sources[j + 1] = 1.0 * rand() / RAND_MAX - 0.5;
      sources[j + 2] = 1.0 * rand() / RAND_MAX - 0.5;
    }
    
    for (int i = 0; i < fmm_config.ntargets; i++) {
      int j = 3 * i;
      targets[j] = 1.0 * rand() / RAND_MAX - 0.5;
      targets[j + 1] = 1.0 * rand() / RAND_MAX - 0.5;
      targets[j + 2] = 1.0 * rand() / RAND_MAX - 0.5;
    }
  } else if (fmm_config.datatype == 2) {
    for (int i = 0; i < fmm_config.nsources; i++) {
      int j = 3 * i;
      charges[i] = 1.0 * rand() / RAND_MAX - 0.5;
      double theta = 1.0 * rand() / RAND_MAX * pi * 0.5;
      double phi = 1.0 * rand() / RAND_MAX * pi * 2;
      sources[j] = sin(theta) * cos(phi);
      sources[j + 1] = sin(theta) * sin(phi);
      sources[j + 2] = cos(theta);
    }

    for (int i = 0; i < fmm_config.ntargets; i++) {
      int j = 3 * i;
      double theta = 1.0 * rand() / RAND_MAX * pi * 0.5;
      double phi = 1.0 * rand() / RAND_MAX * pi;
      targets[j] = sin(theta) * cos(phi);
      targets[j + 1] = sin(theta) * sin(phi);
      targets[j + 2] = cos(theta);
    }
  } else if (fmm_config.datatype == 3) {
    for (int i = 0; i < fmm_config.nsources; i++) {
      int j = 3 * i; 
      charges[i] = 1.0 * rand() / RAND_MAX - 0.5;
      double theta = 1.0 * rand() / RAND_MAX * pi;
      double phi = 1.0 * rand() / RAND_MAX * pi;
      sources[j] = sin(theta) * cos(phi);
      sources[j + 1] = sin(theta) * sin(phi);
      sources[j + 2] = cos(theta);
    }

    for (int i = 0; i < fmm_config.ntargets; i++) {
      int j = 3 * i;
      double theta = 1.0 * rand() / RAND_MAX * pi;
      double phi = (1.0 * rand() / RAND_MAX + 1.0) * pi;
      targets[j] = sin(theta) * cos(phi);
      targets[j + 1] = sin(theta) * sin(phi);
      targets[j + 2] = cos(theta);
    }
  }

  fmm_dag_t *fmm_dag = construct_dag(sources, fmm_config.nsources, 
                                     targets, fmm_config.ntargets, 
                                     fmm_config.s);

  fmm_compute(&fmm_config, fmm_dag, sources, charges, targets, 
              potential, field); 

  // verify numerical accuracy
  int n_verify = (fmm_config.ntargets < 200 ? fmm_config.ntargets : 200); 
  double salg = 0, salg2 = 0, stot = 0, stot2 = 0, errmax = 0;
  double *direct_potential = CALLOC(n_verify, sizeof(double)); 
  double *direct_field     = CALLOC(n_verify * 3, sizeof(double));
  assert(direct_potential != NULL);
  assert(direct_field != NULL); 

#if YUKAWA
  double pi_2 = acos(0); 
  double beta = fmm_config.beta;
#endif 

  CILK_FOR (int i = 0; i < n_verify; i++) {
    int i3 = i * 3;
    const Real_t *t = &targets[i3];
    double pot = 0, fx = 0, fy = 0, fz = 0;
    for (int j = 0; j < fmm_config.nsources; j++) {
      int j3 = j * 3;
      const Real_t *s = &sources[j3];
      const Real_t q = charges[j];
      double rx = t[0] - s[0];
      double ry = t[1] - s[1];
      double rz = t[2] - s[2];
      double rr = rx * rx + ry * ry + rz * rz;

#if LAPLACE 
      double rdis = sqrt(rr);
      if (rr) {
        pot += q / rdis;
        double rmul = q / (rdis * rr);
        fx += rmul * rx;
        fy += rmul * ry;
        fz += rmul * rz;
      }
#elif YUKAWA
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
#endif 
    }

    direct_potential[i]  = pot; 
    direct_field[i3]     = fx;
    direct_field[i3 + 1] = fy;
    direct_field[i3 + 2] = fz;
  }

  for (int i = 0; i < n_verify; i++) {
    int i3 = i * 3; 
    salg += pow(potential[i] - direct_potential[i], 2); 
    stot += pow(direct_potential[i], 2); 
    salg2 += pow(field[i3] - direct_field[i3], 2) +
      pow(field[i3 + 1] - direct_field[i3 + 1], 2) + 
      pow(field[i3 + 2] - direct_field[i3 + 2], 2); 
    stot2 += direct_field[i3]*direct_field[i3] + 
      direct_field[i3 + 1]*direct_field[i3 + 1] + 
      direct_field[i3 + 2]*direct_field[i3 + 2]; 
    errmax = fmax(errmax, fabs(potential[i] - direct_potential[i]));
  }  

  printf("accuracy statistics:\n"
         "  error of potential in L-2 norm: %e\n"
         "  error of potential in L-infty norm: %e\n"
         "  error of field in L-2 norm: %e\n",
         sqrt(salg / stot), errmax, sqrt(salg2 / stot2));

  bool if_pass = (sqrt(salg / stot) <= pow(10.0, -fmm_config.accuracy)) &&
    (sqrt(salg2 / stot2) <= pow(10.0, -fmm_config.accuracy + 1)); 
  printf(if_pass ? "******* pass *******\n" : "******* fail *******\n"); 

  destruct_dag(fmm_dag);
  fmm_cleanup(); 

  FREE(sources);
  FREE(charges);
  FREE(targets);
  FREE(potential);
  FREE(field);
  FREE(direct_potential);
  FREE(direct_field); 

  return 0;
}
