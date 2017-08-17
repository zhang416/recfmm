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
/// @file cilk.h
///
/// @brief Macros to switch between sequential and parallel compilation
///
/// @author Bo Zhang <zhang416 [at] indiana.edu>
/// @author Jingfang Huang <huang [at] amath.unc.edu>
/// @author Nikos P. Pitsianis <nikos [at] cs.duke.edu>
/// @author Xiaobai Sun <xiaobai [at] cs.duke.edu>
///
/// ----------------------------------------------------------------------------

#pragma once
#ifndef CILK_H
#define CILK_H

#ifdef USING_CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_min.h>
#include <cilk/reducer_max.h>

#include "tbb/scalable_allocator.h"

#define CILK_FOR cilk_for
#define CILK_SPAWN cilk_spawn
#define CILK_SYNC cilk_sync
#define CALLOC scalable_calloc
#define FREE scalable_free

#else

#define CILK_FOR for
#define CILK_SPAWN 
#define CILK_SYNC
#define CALLOC calloc
#define FREE free

#endif

#endif 
