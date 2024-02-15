/* Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
 * Copyright (c) 2017, Andreas Herten/Forschungszentrum Jülich
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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"

// Run with:
// ./poisson2d [NITER [NY [NX]]]
int main(int argc, char** argv)
{
    int ny = 1024;
    int nx = 1024;
    int iter_max = 100;
    const real tol = 1.0e-5;

    if (argc >= 2)
    {
        iter_max = atoi( argv[1] );
    }
    if (argc == 3)
    {
        ny = atoi( argv[2] );
        nx = ny;
    }
    if (argc == 4)
    {
        ny = atoi( argv[2] );
        nx = atoi( argv[3] );
    }

    real* restrict const A    = (real*) malloc(nx*ny*sizeof(real));
    real* restrict const Aref = (real*) malloc(nx*ny*sizeof(real));
    real* restrict const Anew = (real*) malloc(nx*ny*sizeof(real));
    real* restrict const rhs  = (real*) malloc(nx*ny*sizeof(real));
    
    // set rhs
    for (int iy = 1; iy < ny-1; iy++)
    {
        for( int ix = 1; ix < nx-1; ix++ )
        {
            const real x = -1.0 + (2.0*ix/(nx-1));
            const real y = -1.0 + (2.0*iy/(ny-1));
            rhs[iy*nx+ix] = expr(-10.0*(x*x + y*y));
        }
    }

    int ix_start = 1;
    int ix_end   = (nx - 1);

    int iy_start = 1;
    int iy_end   = (ny - 1);

    for( int iy = 0; iy < ny; iy++)
    {
        for( int ix = 0; ix < nx; ix++ )
        {
            A[iy*nx+ix] = 0.0;
        }
    }
    // OpenACC Init / Warm-Up
    #pragma acc init
    
    double start;  // Time measurements
    printf("Jacobi relaxation calculation: max %d iterations on %d x %d mesh\n", iter_max, ny, nx);

    printf("Calculate reference solution and time with serial CPU execution.\n");
    start = get_time();
    poisson2d_reference( iter_max, tol, Aref, Anew, nx, ny, rhs );
    double runtime_ref = get_time() - start;

    printf("Calculate current execution.\n");
    start = get_time();
    int iter  = 0;
    real error = 1.0;

    while ( error > tol && iter < iter_max )
    {
        error = 0.0;

        // Jacobi kernel
        #pragma acc parallel loop reduction(max:error)
        for (int ix = ix_start; ix < ix_end; ix++)
        {
            for (int iy = iy_start; iy < iy_end; iy++)
            {
                Anew[iy*nx+ix] = -0.25 * (rhs[iy*nx+ix] - ( A[iy*nx+ix+1] + A[iy*nx+ix-1]
                                                       + A[(iy-1)*nx+ix] + A[(iy+1)*nx+ix] ));
                error = fmaxr( error, fabsr(Anew[iy*nx+ix]-A[iy*nx+ix]));
            }
        }

        // A <-> Anew
        #pragma acc parallel loop
        for (int iy = iy_start; iy < iy_end; iy++)
        {
            for( int ix = ix_start; ix < ix_end; ix++ )
            {
                A[iy*nx+ix] = Anew[iy*nx+ix];
            }
        }

        //Periodic boundary conditions
        #pragma acc parallel loop
        for (int ix = ix_start; ix < ix_end; ix++)
        {
                A[0*nx+ix]      = A[(ny-2)*nx+ix];
                A[(ny-1)*nx+ix] = A[1*nx+ix];
        }
        #pragma acc parallel loop
        for (int iy = iy_start; iy < iy_end; iy++)
        {
                A[iy*nx+0]      = A[iy*nx+(nx-2)];
                A[iy*nx+(nx-1)] = A[iy*nx+1];
        }
        
        if((iter % 100) == 0) printf("%5d, %0.6f\n", iter, error);
        
        iter++;
    }
    double runtime = get_time() - start;

    int errors = 0;

    if (check_results( ix_start, ix_end, iy_start, iy_end, tol, A, Aref, nx ))
    {
        printf( "%dx%d: Ref: %8.4f s, This: %8.4f s, speedup: %8.2f\n", ny, nx, runtime_ref, runtime, runtime_ref/runtime );
    }
    else
    {
        errors = -1;
    }
    
    free(rhs);
    free(Anew);
    free(Aref);
    free(A);

    return errors;
}
