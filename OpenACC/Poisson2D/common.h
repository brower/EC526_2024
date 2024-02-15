/* Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef COMMON_H
#define COMMON_H

#include <assert.h>

#ifdef USE_DOUBLE
    typedef double real;
    #define fmaxr fmax
    #define fabsr fabs
    #define expr exp
    #define MPI_REAL_TYPE MPI_DOUBLE
#else
    typedef float real;
    #define fmaxr fmaxf
    #define fabsr fabsf
    #define expr expf
    #define MPI_REAL_TYPE MPI_FLOAT
#endif

typedef struct
{
    int y;
    int x;
} dim2;

#define MAX_MPI_SIZE 16

static dim2 size_to_size2d_map[MAX_MPI_SIZE+1] = { {0,0},
    {1,1}, {2,1}, {3,1}, {2,2},
    {5,1}, {3,2}, {7,1}, {4,2},
    {3,3}, {5,2}, {11,1}, {6,2},
    {13,1}, {7,2}, {5,3}, {4,4}
};

inline int min( int a, int b)
{
    return a < b ? a : b;
}

inline int max( int a, int b)
{
    return a > b ? a : b;
}

double get_time();

void poisson2d_reference( int iter_max, real tol, real* restrict const Aref, real* restrict const Anew, int nx, int ny, const real* restrict const rhs );

int check_results( int ix_start, int ix_end,  int iy_start, int iy_end, real tol, const real* restrict const A, const real* restrict const Aref, int nx );

static dim2 size_to_2Dsize( int size )
{
    assert(size<=MAX_MPI_SIZE);
    return size_to_size2d_map[size];
}

#endif // COMMON_H
