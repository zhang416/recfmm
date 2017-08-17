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
/// @file fmm-param.h
/// @author Bo Zhang <zhang416 [at] indiana.edu>
/// @author Jingfang Huang <huang [at] amath.unc.edu>
/// @author Nikos P. Pitsianis <nikos [at] cs.duke.edu>
/// @author Xiaobai Sun <xiaobai [at] cs.duke.edu>
/// @date June 18, 2014
/// @brief Contains functions that compute various translation coefficients
/// ----------------------------------------------------------------------------

#pragma once
#ifndef FMM_PARAM_H
#define FMM_PARAM_H

#include "fmm-types.h"

/// ----------------------------------------------------------------------------
/// @brief Constructs the parameter
/// \param[in] fmm_config FMM configuration
/// \param[in] fmm_dag FMM DAG
/// \return The Parameter for the computation. 
/// ----------------------------------------------------------------------------
fmm_param_t *construct_param(const fmm_config_t *fmm_config, 
                             const fmm_dag_t *fmm_dag); 

/// ----------------------------------------------------------------------------
/// @brief Destructs the parameter
/// ----------------------------------------------------------------------------
void destruct_param(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Computes the binomial coefficents binom(p, n), where n = 0, ..., p
/// ----------------------------------------------------------------------------
void bnlcft(double *c, int p);

/// ----------------------------------------------------------------------------
/// @brief Returns a set of Gaussian nodes and weights for integrating the 
///        functions j0(r * x) * exp(z * x) dx over the range x = 0 to 
///        x = infinity
/// ----------------------------------------------------------------------------
void vwts(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Computes number of Fourier modes needed 
/// ----------------------------------------------------------------------------
void numthetahalf(fmm_param_t *fmm_param);

#if LAPLACE

/// ----------------------------------------------------------------------------
/// @brief Creates factorial scaling factors
/// ----------------------------------------------------------------------------
void frmini(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Precomputes the rotation matrix
/// ----------------------------------------------------------------------------
void rotgen(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Implement the fast version of rotation matrices from recursion
/// ----------------------------------------------------------------------------
void fstrtn(int p, double *d, const double *sqc, double theta, int pgsz);

/// ----------------------------------------------------------------------------
/// @brief Computes number of Fourier modes needed 
/// ----------------------------------------------------------------------------
void numthetafour(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Precomputes coefficients needed by the TME operator
/// ----------------------------------------------------------------------------
void rlscini(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Precomputes exponentials needed in the TME->TEE->TEL operation
/// ----------------------------------------------------------------------------
void mkfexp(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Computes the tables of exponentials needed for translating 
///        exponential representations of harmonic functions
/// ----------------------------------------------------------------------------
void mkexps(fmm_param_t *fmm_param);

#elif YUKAWA

/// ----------------------------------------------------------------------------
/// @brief Creates factorial scaling factors
/// ----------------------------------------------------------------------------
void yhfrmini(fmm_param_t *fmm_param); 

/// ----------------------------------------------------------------------------
/// @brief Precomputes the rotation matrix
/// ----------------------------------------------------------------------------
void yhrotgen(fmm_param_t *fmm_param);

/// ----------------------------------------------------------------------------
/// @brief Implement the fast version of rotation matrices from recursion
/// ----------------------------------------------------------------------------
void yhfstrtn(int p, double theta, const double *sqc, double *d, int pgsz);

/// ----------------------------------------------------------------------------
/// @brief Computes number of Fourier modes needed 
/// ----------------------------------------------------------------------------
void numthetafour(fmm_param_t *fmm_param, int level); 

/// ----------------------------------------------------------------------------
/// @brief Precomputes coefficients for shifting multipole expansion
/// \param[in] r0 The shifting distance
/// \param[in] level Partition level where translation occurs
/// ----------------------------------------------------------------------------
void ympshftcoef(fmm_param_t *fmm_param, double r0, int level);

/// ----------------------------------------------------------------------------
/// @brief Precomptues coefficients for shifting local expansion
/// \param[in] r0 The shifting distance
/// \param[in] level Partition level where translation occurs
/// ----------------------------------------------------------------------------
void ylcshftcoef(fmm_param_t *fmm_param, double r0, int level);

/// ----------------------------------------------------------------------------
/// @brief Computes modified Bessel function i_n(z) = sqrt(pi_2/z)*i_(n+1/2)(z) 
/// \param[in] scal The scaling factor to avoid underflow
/// \param[in] x The parameter for i_n
/// \param[in] nb The number of terms for the subindex n
/// \param[out] b Values of i_n(x) for n = 0, ..., nb
/// \param[out] ncalc Error index needed by the fortran routine ribesl_
/// ----------------------------------------------------------------------------
void in(double scal, double x, int nb, double *b, int *ncalc);

/// ----------------------------------------------------------------------------
/// @brief Precomputes exponentials needed in the TME->TEE->TEL operation
/// ----------------------------------------------------------------------------
void ymkfexp(fmm_param_t *fmm_param, int level); 

/// ----------------------------------------------------------------------------
/// @brief Precomputes coefficients needed by the TME operator
/// ----------------------------------------------------------------------------
void yrlscini(fmm_param_t *fmm_param, int level);

/// ----------------------------------------------------------------------------
/// @brief Computes the Legendre function for x > 1 (scaled version). 
/// \param[in] scale The scaling factor
/// \param[in] nmax The number of terms to compute
/// \param[in] x The value at which the Legendre function to be calculated
/// \param[out] y Vector array of output
/// ----------------------------------------------------------------------------
void lgndrgt1(double scale, int nmax, double x, double *y);

/// ----------------------------------------------------------------------------
/// @brief Computes the tables of exponentials needed for translating 
///        exponential representations of harmonic functions
/// ----------------------------------------------------------------------------
void ymkexps(fmm_param_t *fmm_param, int level); 

/// ----------------------------------------------------------------------------
/// @brief Computes Bessel function i_{n + alpha}(x) for nonnegative x
/// \param[in] x Nonnegative real argument for which the Bessel function will be
///              evaluated
/// \param[in] alpha Fractional part of the order of the Bessel function
/// \param[in] nb Integer number of functions to be calculated
/// \param[in] ize If the value is 1, unscaled Bessel function is calculated. 
///                Otherwise, exponentially scaled Bessel function is calculated.
/// \param[out] b Output vector of length nb + 1
/// \param[out] ncalc Integer indicating possible errors
/// ----------------------------------------------------------------------------
void ribesl_(const double *x, double *alpha, int *nb, int *ize, 
             double *b, int *ncalc); 

/// ----------------------------------------------------------------------------
/// @brief Computes scaled spherical Bessel functions y_n(x) for non-negative 
///        argument x, and n = 0, ..., nb
/// \param[in] scal The scaling factor
/// \param[in] x Working precision nonnegative real argument for which y's are
///              to be calculated. 
/// \param[in] nb Integer number of functions to be calculated. 
/// \param[out] by Working precision output vector of length nb + 1. 
/// \param[out] ncalc Integer indicating possible errors. 
/// ----------------------------------------------------------------------------
void kn(double scal, double x, int nb, double *by, int *ncalc);
#endif

#endif

