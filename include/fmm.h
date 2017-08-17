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

/// ----------------------------------------------------------------------------
/// @file fmm.h
/// @author Bo Zhang <zhang416 [at] indiana.edu>
/// @author Jingfang Huang <huang [at] amath.unc.edu>
/// @author Nikos P. Pitsianis <nikos [at] cs.duke.edu>
/// @author Xiaobai Sun <xiaobai [at] cs.duke.edu>
/// @date June 18, 2014
/// @brief recFMM package user interface
/// ----------------------------------------------------------------------------

#pragma once
#ifndef FMM_H
#define FMM_H

#include "fmm-types.h"
#include "cilk.h"

/// ----------------------------------------------------------------------------
/// @brief Generate the source and target trees from input data
/// \param [in] sources The position of source points
/// \param [in] nsources The number of source points
/// \param [in] targets The position of target points
/// \param [in] ntargets The number of target points
/// \param [in] s A parameter controlling whether or not to partition a box
/// \return Pointer to the fmm_dag object
/// ----------------------------------------------------------------------------
fmm_dag_t *construct_dag(const Real_t *sources, const int nsources, 
                         const Real_t *targets, const int ntargets, 
                         const int s); 

/// ----------------------------------------------------------------------------
/// @brief Destruct the DAG
/// \param [in] fmm_dag The dag that will be destructed
/// ----------------------------------------------------------------------------
void destruct_dag(fmm_dag_t *fmm_dag); 

/// ----------------------------------------------------------------------------
/// @brief Compute the potential and force fields according the configuration
/// \param [in] fmm_config Configuration of the computation
/// \param [in] fmm_dag The DAG used for the graph traversal
/// \param [in] sources The position of source points
/// \param [in] charges The charges carried by the source points
/// \param [in] targets The position of target points
/// \param [out] potential The induced potential at each target location
/// \param [out] field The induced force field at each target location
/// ----------------------------------------------------------------------------
void fmm_compute(const fmm_config_t *fmm_config, fmm_dag_t *fmm_dag, 
                 const Real_t *sources, const Real_t *charges, 
                 const Real_t *targets, double *potential, double *field); 

/// ----------------------------------------------------------------------------
/// @brief Cleanup the internal variables used in the FMM computation
/// ----------------------------------------------------------------------------
void fmm_cleanup(void); 

#endif
