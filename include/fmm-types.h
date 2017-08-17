#pragma once
#ifndef FMM_TYPES_H
#define FMM_TYPES_H

#include <complex.h>
#define MAXLEV 128

typedef float real4;
typedef double real8;
typedef real4 Real_t; 

/// ----------------------------------------------------------------------------
/// @brief FMM Box type
/// ----------------------------------------------------------------------------
typedef struct fmm_box_t fmm_box_t; 

struct fmm_box_t {
  int level;  ///< level of the box
  fmm_box_t *parent; ///< parent of the box
  fmm_box_t *child[8]; ///< children of the box
  int nchild; ///< number of children
  int idx; ///< x direction index
  int idy; ///< y direction index
  int idz; ///< z direction index
  int npts; ///< number of points contained 
  int addr; ///< offset to retrieve the first contained point
  fmm_box_t *list1[27]; ///< Coarser or same level List 1 boxes
  fmm_box_t *list5[27]; ///< colleague list
  int nlist1; ///< number of boxes contained in List 1
  int nlist5; ///< number of boxes contained in List 5
  double complex *expansion; ///< multipole/exponential/local expansions
}; 

/// ----------------------------------------------------------------------------
/// @brief FMM DAG type
/// ----------------------------------------------------------------------------
typedef struct {
  int nslev;  ///< level of the source tree
  int nsboxes; ///< total number of boxes on the source tree
  int ntlev; ///< level of the target tree
  int ntboxes; ///< total number of boxes on the target tree
  double size; ///< bounding box dimension
  double corner[3]; ///< corner of the bounding box 
  int *mapsrc; ///< array storing which box each source point is mapped to
  int *maptar; ///< array storing which box each target point is mapped to
  fmm_box_t *source_root; ///< pointer to the root of the source tree
  fmm_box_t *target_root; ///< pointer to the root of the target tree
} fmm_dag_t; 

#if LAPLACE

/// ----------------------------------------------------------------------------
/// @brief FMM Configuration type
/// ----------------------------------------------------------------------------
typedef struct {
  int nsources; ///< number of source points 
  int ntargets; ///< number of target points
  int datatype; ///< type of data to generate
  int accuracy; ///< accuracy requirement
  int s; ///< partition criterion, box with more than s points is partitioned
} fmm_config_t; 

/// ----------------------------------------------------------------------------
/// @brief FMM parameter type
/// ----------------------------------------------------------------------------
typedef struct {
  int pterms; ///< order of the multipole/local expansion
  int nlambs; ///< number of terms in the exponential expansions
  int pgsz; ///< buffer size for holding multipole/local expansion
  int nexptot; ///< total number of exponential expansion terms
  int nthmax; ///< maximum number of Fourier terms in the exponential expansion
  int nexptotp; ///< number of exponential expansions
  int nexpmax; ///< buffer size for holding exponential expansions
  int *numphys; ///< number of modes in the plane-wave expansion
  int *numfour; ///< number of Fourier modes in the expansion
  double *whts; ///< weights for the plane-wave expansion
  double *rlams; ///< nodes for the plane-wave expansion
  double *rdplus; ///< rotation matrix y->z
  double *rdminus; ///< rotation matrix z->x
  double *rdsq3;  /// shifts multipole/local expansion +z direction
  double *rdmsq3; /// shifts multipole/local expansion -z direction
  double *dc; ///< coefficients for local translation along the z-axis
  double *ytopc; ///< precomputed vectors for factorials
  double *ytopcs; ///< precomputed vectors for factorials 
  double *ytopcsinv; ///< precomputed vectors for factorials
  double *rlsc; /// p_n^m for different lambda_k 
  double *zs; ///< TEE operator, z-direction
  double *scale; ///< scaling factor at each level
  double complex *xs; ///< TEE operator, x-direction
  double complex *ys; ///< TEE operator, y-direction
  double complex *fexpe; ///< Coefficients for merging exponentials
  double complex *fexpo; ///< Coefficients for merging exponentials
  double complex *fexpback; ///< Coefficients for merging exponentials
} fmm_param_t; 

#elif YUKAWA

typedef struct {
  int nsources; ///< number of source points 
  int ntargets; ///< number of target points
  int datatype; ///< type of data to generate
  int accuracy; ///< accuracy requirement
  int s; ///< partition criterion, box with more than s points is partitioned
  double beta; ///< parameter of the Yukawa kernel
} fmm_config_t; 

typedef struct {
  int pterms; ///< order of the multipole/local expansion
  int nlambs; ///< number of terms in the exponential expansions
  int pgsz; ///< buffer size for holding multipole/local expansion
  int dcpgsz; ///< buffer size for holding dcu/dcd array
  int *numfour; ///< number of Fourier modes in the expansion 
  int *numphys; ///< number of modes in the plane-wave expansion
  int nexptot; ///< total number of exponential expansion terms
  int nthmax; ///< maximum number of Fourier terms in the exponential expansion
  int mnexptotp; ///< max number of exponential expansions for the worst case
  int nexpmax; ///< buffer size for holding exponential expansions
  int *vnexptot; ///< total number of exponential terms at each level
  int *vnexptotp; ///< number of exponential expansions at each level
  int *vnthmax; ///< maximum number of Fourier terms in the exponential 
                ///< expansion at each level
  double beta; ///< Yukawa kernel parameter
  double *ytop; ///< factorial scaling factors
  double *whts; ///< weights for the plane-wave expansion
  double *rlams; ///< nodes for the plane-wave expansion
  double *rdplus; ///< rotation matrix y->z
  double *rdminus; ///< rotation matrix z->X
  double *rdsq3; /// shifts multipole/local expansion +z direction
  double *rdmsq3; /// shifts multipole/local expansion -z direction
  double *dcu; ///< coefficients for multipole translation along the z-axis
  double *dcd; ///< coefficients for local translation along the z-axis
  double *sfactor;
  double *betascale;
  double *rlsc; /// p_n^m for different lambda_k 
  double scale; 
  double *zs; ///< TEE operator, z-direction
  double complex *xs; ///< TEE operator, x-direction
  double complex *ys; ///< TEE operator, y-direction
  double complex *fexpe; ///< Coefficients for merging exponentials
  double complex *fexpo; ///< Coefficients for merging exponentials
  double complex *fexpback; ///< Coefficients for merging exponentials
} fmm_param_t; 

#endif 

#endif
