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
/// @file fmm-action.h
///
/// @brief Functions traversing the graph and completing FMM computation
///
/// @author Bo Zhang <zhang416 [at] indiana.edu>
/// @author Jingfang Huang <huang [at] amath.unc.edu>
/// @author Nikos P. Pitsianis <nikos [at] cs.duke.edu>
/// @author Xiaobai Sun <xiaobai [at] cs.duke.edu>
/// 
/// ----------------------------------------------------------------------------

#pragma once
#ifndef FMM_ACTION_H
#define FMM_ACTION_H

#include "fmm-types.h"

/// ----------------------------------------------------------------------------
//! @brief Aggregate operation on the source tree
//! \param[in] sbox Pointer to a source box where aggregate operation occurs
/// ----------------------------------------------------------------------------
void aggregate(fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Disaggregate operation on the target tree
/// \param[in] tbox Pointer to a target box where disaggregate operation occurs
/// ----------------------------------------------------------------------------
void disaggregate(fmm_box_t *tbox); 

/// ----------------------------------------------------------------------------
/// @brief Source to multipole translation
/// \param[in] sbox Pointer to the source box where translation occurs
/// ----------------------------------------------------------------------------
void source_to_multipole(fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Multipole to multipole translation
/// \param[in] sbox Pointer to the source box where translation occurs
/// ----------------------------------------------------------------------------
void multipole_to_multipole(fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Multipole to exponential translation
/// \param[in] sbox Pointer to the source box where translation occurs
/// ----------------------------------------------------------------------------
void multipole_to_exponential(fmm_box_t *sbox); 

#if LAPLACE
/// ----------------------------------------------------------------------------
void multipole_to_exponential_p1(const double complex *multipole, 
                                 double complex *mexpu, 
                                 double complex *mexpd);

void multipole_to_exponential_p2(const double complex *mexpf, 
                                 double complex *mexpphys);
#elif YUKAWA
void multipole_to_exponential_p1(const double complex *multipole, 
                                 double complex *mexpu, 
                                 double complex *mexpd, int level);

void multipole_to_exponential_p2(const double complex *mexpf, 
                                 double complex *mexpphys, int level);
#endif 

/// ----------------------------------------------------------------------------
/// @brief Exponential to local translation
/// \param[in] tbox Pointer to a nonleaf target box where translation occurs
/// ----------------------------------------------------------------------------
void exponential_to_local(fmm_box_t *tbox);

#if LAPLACE
void exponential_to_local_p1(const double complex *mexpphys, 
                             double complex *mexpf);

void exponential_to_local_p2(int iexpu, const double complex *mexpu,
                             int iexpd, const double complex *mexpd, 
                             double complex *local); 
#elif YUKAWA
void exponential_to_local_p1(const double complex *mexpphys, 
                             double complex *mexpf, int level);

void exponential_to_local_p2(int iexpu, const double complex *mexpu,
                             int iexpd, const double complex *mexpd, 
                             double complex *local, int level); 
#endif 

/// ----------------------------------------------------------------------------
/// @brief Local to local translation
/// \param[in] tbox Pointer to the target box where translation occurs
/// ----------------------------------------------------------------------------
void local_to_local(fmm_box_t *tbox);

/// ----------------------------------------------------------------------------
/// @brief Local to target translation (evaluating local expansion)
/// \param[in] tbox Pointer to the target box where translation occurs
/// ----------------------------------------------------------------------------
void local_to_target(fmm_box_t *tbox); 

/// ----------------------------------------------------------------------------
/// @brief Source to local translation
/// \param[in] tbox The target box whose local expansion will be updated
/// \param[in] sbox The source box whose points will be used to generate local
///                 expansion of the specified target box
/// ----------------------------------------------------------------------------
void source_to_local(fmm_box_t *tbox, fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Multipole to target translation (evaluating multipole expansion)
/// \param[in] tbox Points in which target box to evaluate the expansion
/// \param[in] sbox Whose multipole expansion will be evaluated
/// ----------------------------------------------------------------------------
void multipole_to_target(fmm_box_t *tbox, fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Process List 1 and List 3 
/// \param[in] tbox Target box whose List 1 and List 3 need to be processed
/// \param[in] sbox Source box which will be processed either as List 1 or List 3
/// ----------------------------------------------------------------------------
void process_list13(fmm_box_t *tbox, fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Pairwise direct interaction between two boxes
///
/// \param[in] tbox Pointer to the target box
/// \param[in] sbox Pointer to the source box
/// ----------------------------------------------------------------------------
void direct_evaluation(const fmm_box_t *tbox, const fmm_box_t *sbox); 

/// ----------------------------------------------------------------------------
/// @brief Evaluates Lengndre polynomial 
/// \param[in] nmax Number of terms in the expansion
/// \param[in] x Where the expansion is evaluated
/// \param[out] y Vector array of output
/// ----------------------------------------------------------------------------
void lgndr(int nmax, double x, double *y); 

void rotz2y(const double complex *multipole, const double *rd, 
            double complex *mrotate); 

void roty2z(const double complex *multipole, const double *rd, 
            double complex *mrotate); 

void rotz2x(const double complex *multipole, const double *rd, 
            double complex *mrotate); 

void make_ulist(int type, fmm_box_t **list, int nlist, 
                int *xoff, int *yoff, double complex *mexpo, int level);

void make_dlist(int type, fmm_box_t **list, int nlist,
                int *xoff, int *yoff, double complex *mexpo, int level); 


#endif 
