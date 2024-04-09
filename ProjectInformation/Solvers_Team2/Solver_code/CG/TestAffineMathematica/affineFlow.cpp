/*======================================================================================
  Rich Brower Sat Apr  8 13:59:53 EDT 2017
  Sun Dec 10 21:19:52 EST 2023 Do for flow on fundamental Icosahedral face 

Use Baricentric co-ordinates
better indices. 

L = 4 example old nx = 0,..,L, ny = 0,...,L

Fidual Grid points xi1 = nx/L  xi1 = ny/L  xi1 = nz/L                  
lattice  0 <= i1 <= L   and i1 +i2 + i3 = L
         0 <= i1 <= L 
          0 <= i1 <= L   

Displaed co-ordinate

X[n1][n2] =L *(  xi r1 + xi2 r2 + xi3 r3)     xi1 + xi2  + xi3  = 1 

Flat co-ordinates :  xi1 = nx/L     xi2 = ny/L   = xi3 =  nz/L =1 -  nx/L + ny/L   
 
X[n1][n2] = nx r1 + ny  r2 + nz r3     nx + ny + nz = L

Projected vectors

rvecX[n1][n2] =  R*X[n1][n2]/|X[n1][n2]|

 
                     
                    (2) 
                    *   
                  *    *     i3 = 0
 i1   = 0       *    *    *                  No of Triangles = L*L{xco
             *    *    *    *
           *    *    *    *   * 
    (3)          i2 = 0        (1)   (e.g. i2 = nx , i1 = ny, i3 = nz)

Flart plane vectors 
rvec0[i1][i2] = xi1 r1 +  xi2 r2 +  xi2 r3
Spherical projection
            
nx-->   0    1    2    3    4   =   L  
Ntriagle = L*L

Nsite = (L+1)*(L+2)/2     3*4/2 = 6,  4*5/2 = 10
Edge = 3(L+1)  &  interior  (L+1)*(L+2)/2 -  3(L+1) = (L-4)(L+1)/2

rvec[i2][i1][3] = i1 r1 + i2 r2 + i3 r3  (unit lttice spacing?)

------------------------------------


Gsuss Newton -- See Andreqw Gibiansky!

Minimize    N sites, 2 N Faces and 3 N  edges    N - E + F = 0

E[x] =(1/2) \sum_k ( g_k(x) )^2  \simeq  ( g_k(x0) + dg_k/dx_i dx_i + dg_k/dx_idx_j dx_i dx_j/2 )^2 
 
      =  g_k(x0)^2 + 2 g_k(x0)  (dg_k/dx_i) dx_i +  dg_k/dx_j dg_k/dx_i dx_j  dx_i + 2 g_k(x0)  dg_k/dx_idx_j dx_i dx_j/2

Iterate  0 = g^2_k + 2 g_k A_{ki} dx_i ==>   g_k^2 + 2 g_k A_{ki} dx_i

Look for null state and/or SVD/

This is the same condition with full Hessian: Grad S=x[x] = 0 = \sum_k 2 g_k(x) \Grad G_k = 0 and othogonal projection g_k J_k_i = 0
does not imply g_k = 0 -- just left vector in the null space. Think SVD. There can be flat direction (zero
e.v) for square case or insufficient constratians to imply  g_k = 0 all k.

J_ki =   dg_k/dx_j   Relax normal equation or Psuedo inverse for over constrained.

S[x0 + x] = (J_ki x_i - b_k)^2  by standard CG on normal equations. (J^T J)_ij x_j = (J^T)_ik b_k

x' = x0 + x   f^k(x) \simeq f^k(x_0) +  D[f_k(x0),x_i] \sum_k = J^k_i dx^i - b^k = J x - b

Then shift x

Array layout  ls[y][x][mu] ---xvec[i] with i  mu + Three*x + Three*Lx*y; 

TODO:  Replace Finite Difference by Exact Derivatives.

------------------------------------------------------
Conventions:  Let a =1 so refiment  changes R.

Equlateral triangle  A  = ell^2 sqrt(3)/4 so can define

a^2 A F = 4 Pi R^2 == 
a^2/R^2 = 16 Pi/(sqrt(3) F)   = 8Pi/(sqrt(3) 10 L^2)

Ref "Radial lattice quantization of 3D ϕ 4 field theory." Physical Review D 104.9 (2021): 094502. 

uses a^2/R^2  8Pi/(sqrt(3) 10L^2 ) =  8Pi/(sqrt(3)( 2 +10 L^2))


TODO: Convert arrays to vectors: Why ? Better to pass args? 

std::vector<std::vector<double>> fourier_2pt_sum(k_max + 1);
  for (int k = 0; k <= k_max; k++) {
    fourier_2pt[k].resize(Ny);
    fourier_2pt_sum[k].resize(Ny);
  }

 See old ising_graph for tricks 

rotations about axis n.
https://stackoverflow.com/questions/7582398/rotate-a-vector-about-another-vector


FIRST PASS: Only Area miniminzation and L no factors of 3. 
======================================================================================*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <climits>
//#include "triangle.h"

using namespace std;

#define  L  32
#define TWOPI  6.283185307179586
#define Root2 1.4142135623730951

#define Three 3
#define Debug 0

int printGeometry( double vec[L+1][L+1][Three] ,  FILE *fptlocal);

int main(int argc, char *argv[])
{
 
  double ls[L+1][L+1][Three];
  double xvec[L+1][L+1][Three]; //  vectors on plane at fixed z
  double rvec[L+1][L+1][Three];  // vector on sphere
  int N = L*L;
  
  //Inintial plane
  
  //co-ordinates for iscohedral triangle
  // Cos[2 Pi/3] = -1/2, Sin[ 2 Pi/3] = sqrt(3)/2
  // Nearest Neibor angle is Cos[theta] = 1/Sqrt[5]]
  
  double R = 1.0;
  double a = 1.0;
  double z =   sqrt((7 + 3 * sqrt(5))/8);  // iscosahedron
  // Got to choose frame correctly starting here!  (x,y,z) ord
   // fiducial triangle
  z= sqrt((7 + 3 * sqrt(5))/8);
  double r1[3] =  { sqrt(3)/2.0 , -1/2.0,  z };
  double r2[3] =  {           0,   1.0 ,   z };
  double r3[3] =  { -sqrt(3)/2.0, -1/2.0,  z };

  //Plot xvec = nx r1 + ny r2 + nz r3  

  for(int mu =0; mu < 3; mu++)
    {
      r1[mu] = r1[mu]/sqrt(1 + z*z);
      r2[mu] = r2[mu]/sqrt(1 + z*z);
      r3[mu] = r3[mu]/sqrt(1 + z*z);
    }

 
 
  srand(137); // cout<< "RAND_MAX: " << RAND_MAX << '\n';
  cout << N <<endl;
  
  // Set initial plane fundamental zone:  0 <= ny <= nx <= nz <= L
  
  int nz = 0;

  #if 1  // Orginal code
  for(int nx = 0; nx <= L; nx++)
    {
      for(int ny = 0; ny<= L; ny++)
	{   
	  nz = L -nx - ny; 
	  for(int mu =0; mu < 3; mu++)
	    {
	      xvec[ny][nx][mu] = nx*r1[mu] +  ny*r2[mu] +   nz*r3[mu]; 
	      rvec[ny][nx][mu] = nx*r1[mu] +  ny*r2[mu] +   nz*r3[mu]; 
	      ls[ny][nx][mu]     = a;     
	    }
	}
    }
#endif

  /***************************************
how to permute variabls inside mu loop: Reassign indeces as this.

nz | [ny][nx] & ny | [nx] [nz]   &   nx | [nz ][ny]    rotation,
nz | [nx][ny] & nx | [ny] [nz]   &   ny | [nz ][nx]    reflect rotation


******************************************/  
  
#if 0  //try symmetric form
  for(int nx = 0; nx <= L; nx++)
    {
      for(int ny = 0; ny<= L; ny++)
	{
	  for(int nz = 0; nz<= L; nz++)
	    { 
	      if(nx + ny + nz == L )
		{
		  for(int mu =0; mu < 3; mu++)
		    {
		      xvec[ny][nx][mu] = nx*r1[mu] +  ny*r2[mu] +   nz*r3[mu]; 
		      rvec[ny][nx][mu] = nx*r1[mu] +  ny*r2[mu] +   nz*r3[mu]; 
		      ls[ny][nx][mu]     = a;     
		    }
		}
	    }
	}
    }
#endif
  
  double x_sqr = 1.0;

  // double angle =  TWOPI/64.0;
    double angle = TWOPI/8.0;  // flat

/*==========================

z streth   r2-->  {0, r2[1], w r[2]}   

R^2 = r2[1]^2 + w^2 r[2]^2 = (32 + (47 + 21 Sqrt[5]) w^2)/(79 + 21 Sqrt[5])

==========================*/ 
 #if 1 
  for(int ny = 0; ny<= L; ny++)
    {
      for(int nx = 0; nx <= L; nx++)
	{     nz = L - nx - ny;
	      xvec[ny][nx][0] = Root2*sin(angle)*( nx*r1[0] +  ny*r2[0] +  nz*r3[0] );
	      xvec[ny][nx][1] = Root2*sin(angle)*( nx*r1[1] +  ny*r2[1] +  nz*r3[1] );
	      xvec[ny][nx][2] = Root2*cos(angle)*( nx*r1[2] +  ny*r2[2] +  nz*r3[2] );  
              x_sqr = xvec[ny][nx][0]*xvec[ny][nx][0]
	      	+ xvec[ny][nx][1]*xvec[ny][nx][1] +xvec[ny][nx][2]*xvec[ny][nx][2];	  
	  for(int mu =0; mu < 3; mu++)
	    {	    
	      rvec[ny][nx][mu] = xvec[ny][nx][mu]/sqrt(x_sqr);
	      ls[ny][nx][mu]     = a;     
	    }
	}
    }

 #endif

  #if 1  
  FILE* fptr = NULL;  // C style
  char out_name[64];
  sprintf(out_name,"data/Triangle_%d.dat",L); // filename
  fptr = fopen(out_name,"w");
  if(fptr == NULL)
    {
      printf("Error!");   
      exit(1);             
    }

   for(int ny = 0; ny <= L-1; ny++)
    for(int nx = 0; nx <= L-1; nx++)
      { // if(ny <= 1 - ny - nx)
	{
	  fprintf(fptr,"  %f  %f  %f     %f     %f  ",  xvec[ny][nx][0], xvec[ny][nx][1], 
		  rvec[ny][nx][0],  rvec[ny][nx][1],  rvec[ny][nx][2]);
	  fprintf(fptr,"\n");
	}
      }

   fclose(fptr);
#endif
  
  /* Loop for a hack to experiment with  tracking  evolution of rvec on Mathematica

double r1[3] =  { sqrt(3)/2.0 , -1/2.0,  z };
  double r2[3] =  {           0,   1.0 ,   z };
  double r3[3] =  { -sqrt(3)/2.0, -1/2.0,  z };
  Plot xvec = nx r1 + ny r2 + nz r3  
Simplify[nx {Sqrt[3]/2, -1/2, z} +  ny {0, 1, z} + (L - nx - ny) {Sqrt[3]/2, -1/2, z}]

{1/2 Sqrt[3] (L - ny), 1/2 (-L + 3 ny), L z}
 
  */


  /************* BEGIN EVOLUTION *******************/

  double Rsqr = 1.0;
  double w = 1.0;
  for(int iter =0 ; iter < 4; iter++)  
    {  w = 10*iter + 1.0;
      Rsqr = r2[0]*r2[0] + r2[1]*r2[1] + w*w*r2[2]*r2[2];
  for(int ny = 0; ny<= L; ny++)
    {
      for(int nx = 0; nx <= L; nx++)
	{   
	      nz = L - nx - ny;
	      xvec[ny][nx][2] = w*( nx*r1[2] +  ny*r2[2] +  nz*r3[2] );  
	      x_sqr = (xvec[ny][nx][0]*xvec[ny][nx][0]
		       + xvec[ny][nx][1]*xvec[ny][nx][1] +xvec[ny][nx][2]*xvec[ny][nx][2]);
	       
	  for(int mu =0; mu < 3; mu++)
	    {	    
	      rvec[ny][nx][mu] = sqrt(Rsqr)*xvec[ny][nx][mu]/sqrt(x_sqr);
	      ls[ny][nx][mu]     = a;     
	    }
	}
    }
  
  FILE* fptr = NULL;  // C style
  char out_name[64];
  sprintf(out_name,"data/Triangle_%d_%d.dat",L,iter); // filename
  fptr = fopen(out_name,"w");
  if(fptr == NULL)
    {
      printf("Error!");   
      exit(1);             
    }
  
  //fprintf(fptr,"#  nx   ny   rvec[ny][nx][0] rvec[ny][nx][1]  rvec[ny][nx][2] \n");
  
  for(int ny = 0; ny <= L; ny++)
    for(int nx = 0; nx <= L -ny; nx++)  // Fidual Trangle  
      { 
	{
	  fprintf(fptr,"  %f  %f  %f     %f     %f  ",    xvec[ny][nx][0],  xvec[ny][nx][1],
		  rvec[ny][nx][0],  rvec[ny][nx][1],  rvec[ny][nx][2]);
	  fprintf(fptr,"\n");
	}
      }
  
  fclose(fptr);
}
 /*************  END EVOLUTION *****************/


   
   // test generic file write

 #if 0
  FILE * fptrnew;
   char out_arg[64];
  sprintf(out_arg,"data/TrianglTest_%d.dat",L); // filename
  fptrnew = fopen(out_arg,"w");
  printGeometry(  rvec, fptrnew );
  fclose(fptrnew);
#endif

  #if 0  
  FILE* fptrFull = NULL;  // C style
  char out_nameFull[64];
sprintf(out_nameFull,"data/TriangleFull_%d.dat",L); // filename
  fptrFull = fopen(out_nameFull,"w");
  if(fptrFull == NULL)
    {
      printf("Error!");   
      exit(1);             
    }

   for(int ny = 0; ny <= L-1; ny++)
    for(int nx = 0; nx <= L-1; nx++)
      { // if(ny <= 1 - ny - nx)
	{
	  fprintf(fptrFull,"  %f  %f  %f     %f     %f  ",  xvec[ny][nx][0], xvec[ny][nx][1], 
		  rvec[ny][nx][0],  rvec[ny][nx][1],  rvec[ny][nx][2]);
	  fprintf(fptrFull,"\n");
	}
      }

   fclose(fptrFull);
#endif
  
   return 0;
}

#if 1
int printGeometry( double vec[L+1][L+1][Three] ,  FILE *fptlocal)
{
 
  for(int ny = 0; ny <= L; ny++)
    for(int nx = 0; nx <= L -ny; nx++)
      { 
	{
	  fprintf(fptlocal,"  %f  %f  %f     %f     %f  ",  vec[ny][nx][0], vec[ny][nx][1],
		  vec[ny][nx][0],  vec[ny][nx][1],  vec[ny][nx][2]);
	  fprintf(fptlocal,"\n");
	}
      }
  
  return 0;
}

#endif



