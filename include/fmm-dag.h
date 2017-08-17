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
/// @brief Contains interface functions to operate on the DAG
/// ----------------------------------------------------------------------------

#pragma once
#ifndef FMM_DAG_H
#define FMM_DAG_H

#include <stdbool.h>
#include "fmm-types.h"

/// ----------------------------------------------------------------------------
/// @brief Parition a given box
/// \param [in] box Pointer to the box requiring partition
/// \param [in] type Type of the box, source or target
/// \param [in] points Position of the points contained in the box
/// \param [out] level Number of levels below box 
/// \param [out] nboxes Number of boxes below the branch rooted at box
/// ----------------------------------------------------------------------------
void partition_box(fmm_box_t *box, const char type, const Real_t *points, 
		               int *level, int *nboxes); 

/// ----------------------------------------------------------------------------
/// @brief Generate the merged interaction lists 
/// \param [in] tbox A nonleaf target box, the merged interaction lists of whose
///                  child boxes are being built
/// \param [out] uall Pointers of the boxes in the up-all list
/// \param [out] nuall Number of boxes contained in the up-all list
/// \param [out] xuall x-direction offset of boxes in the up-all list
/// \param [out] yuall y-direction offset of boxes in the up-all list
/// \param [out] u1234 Pointers of the boxes in the u1234 list
/// \param [out] nu1234 Number of boxes contained in the u1234 list
/// \param [out] x1234 x-direction offset of boxes in the u1234 list
/// \param [out] y1234 y-direction offset of boxes in the u1234 list
/// \param [out] dall Pointers of the boxes in the down-all list
/// \param [out] ndall Number of boxes contained in the down-all list
/// \param [out] xdall x-direction offset of boxes in the down-all list
/// \param [out] ydall y-direction offset of boxes in the down-all list
/// \param [out] d5678 Pointers of the boxes in the d5678 list
/// \param [out] nd5678 Number of boxes contained in the d5678 list
/// \param [out] x5678 x-direction offset of boxes in the d5678 list
/// \param [out] y5678 y-direction offset of boxes in the d5678 list
/// \param [out] nall Pointers of the boxes in the north-all list
/// \param [out] nnall Number of boxes contained in the north-all list
/// \param [out] xnall x-direction offset of boxes in the north-all list
/// \param [out] ynall y-direction offset of boxes in the north-all list
/// \param [out] n1256 Pointers of the boxes in the north-1256 list
/// \param [out] nn1256 Number of boxes contained in the north-1256 list
/// \param [out] x1256 x-direction offset of boxes in the north-1256 list
/// \param [out] y1256 y-direction offset of boxes in the north-1256 list
/// \param [out] n12 Pointers of the boxes in the north-12 list
/// \param [out] nn12 Number of boxes contained in the north-12 list
/// \param [out] x12 x-direction offset of boxes in the north-12 list
/// \param [out] y12 y-direction offset of boxes in the north-12 list
/// \param [out] n56 Pointers of the boxes in the north-56 list
/// \param [out] nn56 Number of boxes contained in the north-56 list
/// \param [out] x56 x-direction offset of boxes in the north-56 list
/// \param [out] y56 y-direction offset of boxes in the north-56 list
/// \param [out] sall Pointers of the boxes in the south-all list
/// \param [out] nsall Number of boxes contained in the south-all list
/// \param [out] xsall x-direction offset of boxes in the south-all list
/// \param [out] ysall y-direction offset of boxes in the south-all list
/// \param [out] s3478 Pointers of the boxes in the south-3478 list
/// \param [out] ns3478 Number of boxes contained in the south-3478 list
/// \param [out] x3478 x-direction offset of boxes in the south-3478 list
/// \param [out] y3478 y-direction offset of boxes in the south-3478 list
/// \param [out] s34 Pointers of the boxes in the south-34 list
/// \param [out] ns34 Number of boxes contained in the south-34 list
/// \param [out] x34 x-direction offset of boxes in the south-34 list
/// \param [out] y34 y-direction offset of boxes in the south-34 list
/// \param [out] s78 Pointers of the boxes in the south-78 list
/// \param [out] ns78 Number of boxes contained in the south-78 list
/// \param [out] x78 x-direction offset of boxes in the south-78 list
/// \param [out] y78 y-direction offset of boxes in the south-78 list
/// \param [out] eall Pointers of the boxes in the east-all list
/// \param [out] neall Number of boxes contained in the east-all list
/// \param [out] xeall x-direction offset of boxes in the east-all list
/// \param [out] yeall y-direction offset of boxes in the east-all list
/// \param [out] e1357 Pointers of the boxes in the east-1357 list
/// \param [out] ne1357 Number of boxes contained in the east-1357 list
/// \param [out] x1357 x-direction offset of boxes in the east-1357 list
/// \param [out] y1357 y-direction offset of boxes in the east-1357 list
/// \param [out] e13 Pointers of the boxes in the east-13 list
/// \param [out] ne13 Number of boxes contained in the east-13 list
/// \param [out] x13 x-direction offset of boxes in the east-13 list
/// \param [out] y13 y-direction offset of boxes in the east-13 list
/// \param [out] e57 Pointers of the boxes in the east-57 list
/// \param [out] ne57 Number of boxes contained in the east-57 list
/// \param [out] x57 x-direction offset of boxes in the east-57 list
/// \param [out] y57 y-direction offset of boxes in the east-57 list
/// \param [out] e1 Pointers of the boxes in the east-1 list
/// \param [out] ne1 Number of boxes contained in the east-1 list
/// \param [out] x1 x-direction offset of boxes in the east-1 list
/// \param [out] y1 y-direction offset of boxes in the east-1 list
/// \param [out] e3 Pointers of the boxes in the east-3 list
/// \param [out] ne3 Number of boxes contained in the east-3 list
/// \param [out] x3 x-direction offset of boxes in the east-3 list
/// \param [out] y3 y-direction offset of boxes in the east-3 list
/// \param [out] e5 Pointers of the boxes in the east-5 list
/// \param [out] ne5 Number of boxes contained in the east-5 list
/// \param [out] x5 x-direction offset of boxes in the east-5 list
/// \param [out] y5 y-direction offset of boxes in the east-5 list
/// \param [out] e7 Pointers of the boxes in the east-7 list
/// \param [out] ne7 Number of boxes contained in the east-7 list
/// \param [out] x7 x-direction offset of boxes in the east-7 list
/// \param [out] y7 y-direction offset of boxes in the east-7 list
/// \param [out] wall Pointers of the boxes in the west-all list
/// \param [out] nwall Number of boxes contained in the west-all list
/// \param [out] xwall x-direction offset of boxes in the west-all list
/// \param [out] ywall y-direction offset of boxes in the west-all list
/// \param [out] w2468 Pointers of the boxes in the west-2468 list
/// \param [out] nw2468 Number of boxes contained in the west-2468 list
/// \param [out] x2468 x-direction offset of boxes in the west-2468 list
/// \param [out] y2468 y-direction offset of boxes in the west-2468 list
/// \param [out] w24 Pointers of the boxes in the west-24 list
/// \param [out] nw24 Number of boxes contained in the west-24 list
/// \param [out] x24 x-direction offset of boxes in the west-24 list
/// \param [out] y24 y-direction offset of boxes in the west-24 list
/// \param [out] w68 Pointers of the boxes in the west-68 list
/// \param [out] nw68 Number of boxes contained in the west-68 list
/// \param [out] x68 x-direction offset of boxes in the west-68 list
/// \param [out] y68 y-direction offset of boxes in the west-68 list
/// \param [out] w2 Pointers of the boxes in the west-2 list
/// \param [out] nw2 Number of boxes contained in the west-2 list
/// \param [out] x2 x-direction offset of boxes in the west-2 list
/// \param [out] y2 y-direction offset of boxes in the west-2 list
/// \param [out] w4 Pointers of the boxes in the west-4 list
/// \param [out] nw4 Number of boxes contained in the west-4 list
/// \param [out] x4 x-direction offset of boxes in the west-4 list
/// \param [out] y4 y-direction offset of boxes in the west-4 list
/// \param [out] w6 Pointers of the boxes in the west-6 list
/// \param [out] nw6 Number of boxes contained in the west-6 list
/// \param [out] x6 x-direction offset of boxes in the west-6 list
/// \param [out] y6 y-direction offset of boxes in the west-6 list
/// \param [out] w8 Pointers of the boxes in the west-8 list
/// \param [out] nw8 Number of boxes contained in the west-8 list
/// \param [out] x8 x-direction offset of boxes in the west-8 list
/// \param [out] y8 y-direction offset of boxes in the west-8 list
/// ----------------------------------------------------------------------------
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
			             fmm_box_t *w8[3], int *nw8, int *x8, int *y8);

/// ----------------------------------------------------------------------------
/// @brief Delete the branch rooted at the input box
/// \param [in] box Pointer to the box the branch below which will be deleted
/// ----------------------------------------------------------------------------
void delete_box(fmm_box_t *box); 

/// ----------------------------------------------------------------------------
/// @brief Determine whether two boxes are adjacent
/// \param [in] box1 The first box 
/// \param [in] box2 The second box, and its size is no larger than the first.
/// \return true if adjacent, false otherwise. 
/// ----------------------------------------------------------------------------
bool is_adjacent(const fmm_box_t *box1, const fmm_box_t *box2); 
#endif 
