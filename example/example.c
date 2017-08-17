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
#include "fmm.h"

void update_postition(const double *potential, const double *field, 
                      Real_t *sources, Real_t *charges) {
  // implement your favorite time-marching scheme
}

int main(int argc, char *argv[]) {
  fmm_config_t fmm_config = {
    .nsources = 10000,
    .ntargets = 10000,
    .datatype = 1,
    .accuracy = 3,
    .s        = 40
  };

  // Allocate memory to hold particle information and output results
  Real_t *sources = CALLOC(fmm_config.nsources * 3, sizeof(Real_t)); 
  Real_t *charges = CALLOC(fmm_config.nsources, sizeof(Real_t)); 
  Real_t *targets = CALLOC(fmm_config.ntargets * 3, sizeof(Real_t)); 
  double *potential = CALLOC(fmm_config.ntargets, sizeof(double)); 
  double *field = CALLOC(fmm_config.ntargets * 3, sizeof(double)); 
  
  // Generate input data
  double pi = acos(-1);   
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

  // The package can be used in a time-marching or an iterative solver 
  // configuration. In a time-marching scheme, DAG is regenerated in each step. 
  // In an iterative solver configuration, DAG remains the same. 

  // We demonstrate the time-marching configuration here
  double t0 = 0; 
  double tf = 1; 
  double dt = 1; 
  fmm_dag_t *fmm_dag = NULL;

  for (double t = t0; t < tf; t += dt) {
    // destruct the previous DAG if it exists
    if (fmm_dag != NULL)
      destruct_dag(fmm_dag); 

    // partition the source and target ensemble.    
    fmm_dag = construct_dag(sources, fmm_config.nsources, 
                            targets, fmm_config.ntargets, 
                            fmm_config.s); 

    // compute the potential and field information
    fmm_compute(&fmm_config, fmm_dag, sources, charges, targets, 
                potential, field); 

    // update particle position
    update_postition(potential, field, sources, charges); 
  }

  // destruct the final DAG
  destruct_dag(fmm_dag); 

  // clean up the FMM internal variables
  fmm_cleanup(); 
  
  FREE(sources);
  FREE(charges);
  FREE(targets);
  FREE(potential);
  FREE(field);

  return 0;
}
