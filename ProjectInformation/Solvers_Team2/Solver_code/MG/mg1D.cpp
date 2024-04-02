/*======================================================================================
Rich Brower April 21, 2020

C program based on simple FMV code of  S. F. McCormick, Multigrid Methods, Frontiers in Applied! Mathematics, vol. 3, SIAM Books, Philadelphia, 1987.

We solve the 1-d Laplace problem:   

 "A phi = b"   with lattice spacing a = h = 1.

As we iterate to reduce "res = b -A phi" we reduce the error:

"phi  = phi_exact + error" 

by  "A error = res"

A error = A (phi_exact - phi)  = b - A phi = residue

(The code just call error  == phi! No need to have two arrays!)

--- Ok here is the equation:

 A phi(x)  =   [2 phi(x) - phi(x+a) - phi(x-a)]/a^2 + m^2 phi(x) = b(x)

in MG we consider latice a --> 2 a --> 4 a ---> a^level    level = 0, 1, ..., nlev

 so we scale by scale(a) = 1/(2 + a^2 m^2)   (a = h)

for the Jacobi iteration:

phi(x) =scale(a) [phi(x+a) + phi(x-a)] + a^2 scale(a)  b(x)

res(x) =  a^2 scale(a)  b(x) - [phi(x) - scale(a) [phi(x+a) + phi(x-a)]


converting error -- residue form we have

err(x) = scale(a) [err(x+a) + err(x-a)] + res(x)  

(See Relax routine below where code calls err = phi! Sorry)

======================================================================================*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
  //#include <stdio.h>
  //#include <stdlib.h>
  //#include <math.h>
  //#include <complex.h>
using namespace std;

FILE * pFile = fopen ("mySolution.txt","w"); 
  
typedef struct{
    int N;
    int Lmax;
    int size[20];
    double a[20];
    double m;
    double scale[20];
  } param_t;

void relax(double *phi, double *res, int lev, int niter, param_t p);
void proj_res(double *res_c, double *rec_f, double *phi_f, int lev,param_t p);
void inter_add(double *phi_f, double *phi_c, int lev,param_t p);
double GetResRoot(double *phi, double *res, int lev, param_t p);

int main()
{  
  double *phi[20], *res[20];
  param_t p;
  int nlev;
  int i,lev;
  
  //set parameters________________________________________
  p.Lmax = 10; // max number of levelsx
  p.N = 2*(int)pow(2.,p.Lmax);  // MUST BE POWER OF  eg. p.Lmax = 11 x==> p.N = 2^{12} = 4048
  p.m = 0.001;
  nlev = 10; // NUMBER OF LEVELS:  nlev = 0 give top level alone
 
  if(nlev  > p.Lmax){ 
    printf("ERROR More levels than available in lattice! \n");
    return 0; }
  
  printf("\n V cycle for  %d lattice with nlev = %d out of max  %d \n \n", p.N, nlev, p.Lmax); 
  
  // initialize arrays__________________________________
  p.size[0] = p.N;
  p.a[0] = 1.0;
  p.scale[0] = 1.0/(2.0 + p.m*p.m);
  
  for(lev = 1;lev< p.Lmax+1; lev++) {
    p.size[lev] = p.size[lev-1]/2;
    p.a[lev] = 2.0 * p.a[lev-1];
    //  p.scale[lev] = 1.0/(2.0 + p.m*p.m*p.a[lev]*p.a[lev]);
    p.scale[lev] = 1.0/(2.0 + p.m*p.m);  // This seem to work better?? Why I am not sure!
  }


  for(lev = 0;lev< p.Lmax+1; lev++)
    {
      phi[lev] = (double *) malloc(p.size[lev] * sizeof(double)); //over ride with error to save space
      res[lev] = (double *) malloc(p.size[lev] * sizeof(double));
      for(i = 0;i< p.size[lev];i++)
       {
	  phi[lev][i] = 0.0;
          res[lev][i] = 0.0;
	};
    }  

   
  res[0][p.N/4] = 1.0*p.scale[0];  //Two unit point middle of fist half  N/2 and anti-source at reflected 
  res[0][3*p.N/4] = - 1.0*p.scale[0];
  
  // iterate to solve_____________________________________
  double resmag = 1.0; // not rescaled.
  int ncycle = 0; 
  int n_per_lev = 10;
  resmag = GetResRoot(phi[0],res[0],0,p);
  printf("At the %d cycle the mag residue is %g \n",ncycle,resmag);
 

  
  while(resmag > 0.000001 && ncycle < 100000 )
     // while(resmag > 0.00001)
    { 
      ncycle +=1; 
      for(lev = 0;lev<nlev; lev++)   //go down 
	{    
	   relax(phi[lev],res[lev],lev, n_per_lev,p); // lev = 1, ..., nlev-1  
	   proj_res(res[lev + 1], res[lev], phi[lev], lev,p);    // res[lev+1] += P^dag res[lev]
	}

      for(lev = nlev;lev >= 0; lev--)  //come up
	{ 
  	  relax(phi[lev],res[lev],lev, n_per_lev,p);   // lev = nlev -1, ... 0;
	  if(lev > 0) inter_add(phi[lev-1], phi[lev], lev, p);   // phi[lev-1] += error = P phi[lev] and set phi[lev] = 0;
	}
      resmag = GetResRoot(phi[0],res[0],0,p);
      printf("At the %d cycle the mag residue is %g \n",ncycle,resmag);
    }


   int L  = p.size[0];
   for(int i = 0; i < p.N; i++)
     {
      double residual = res[0][i]/p.scale[0] - phi[0][i]/p.scale[0]  
	          + (phi[0][(i+1)%L] + phi[0][(i-1+L)%L] );
     fprintf(pFile, " %d     %g     %g      \n ",i,phi[0][i],residual);
     }
   
   fclose (pFile);
   
  return 0;
}

void relax(double *phi, double *res, int lev, int niter, param_t p)
{  
  int i, x,y;
   int L;
   L  = p.size[lev];
  
  for(i=0; i<niter; i++)    
    for(x = 0; x < L; x++)
	{
	  phi[x] = res[x] + p.scale[lev] * (phi[(x+1)%L] + phi[(x-1+L)%L] );			
       	}
  return;    
}

void proj_res(double *res_c, double *res_f, double *phi_f,int lev,param_t p)
{  
  int L, Lc, f_off, c_off, x, y;
  L = p.size[lev];
  double r[L]; // temp residue
  Lc = p.size[lev+1];  // course level
  
  //get residue
  for(x = 0; x< L; x++)
      r[x ] = res_f[x ] -  phi_f[x ] + p.scale[lev]*(phi_f[(x+1)%L] + phi_f[(x-1+L)%L]);
  
  //project residue
  for(x = 0; x< Lc; x++)  res_c[x ] = 0.5*( r[2*x]  + r[(2*x + 1)%L ] );

  return;
}

void inter_add(double *phi_f,double *phi_c,int lev,param_t p)
{  
  int L, Lc, x, y;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1]; 
  
  for(x = 0; x< Lc; x++)
      {
	phi_f[2*x]              += phi_c[x];
	phi_f[(2*x + 1)%L ]     += phi_c[x];
      }
  //set to zero so phi = error 
  for(x = 0; x< Lc; x++) phi_c[x + y*L] = 0.0;
  
  return;
}

#if 0
double GetResRoot(double *phi, double *res, int lev, param_t p)
{ //true residue
  int i, x;
  double residue;
  double ResRoot = 0.0;
  int L;
  if(lev != 0) cout<<"ERROR ON TOP LEVEL RES " << endl;
  L  = p.size[lev];
  for(x = 0; x < L; x++){
    residue = p.a[lev]*p.a[lev]*p.scale[lev]*b[x]-  phi[x]
            +  p.scale[0]*(phi[(x+1)%L] + phi[(x-1+L)%L]) ;
      ResRoot += residue*residue; // true residue
    }

  return sqrt(ResRoot);    
}
#endif

double GetResRoot(double *phi, double *res, int lev, param_t p)
{ //true residue
  int i, x,y;
  double residue;
  double ResRoot = 0.0;
  int L;
  L  = p.size[lev];
  
  for(x = 0; x < L; x++) {
      residue = res[x]/p.scale[lev] - phi[x]/p.scale[lev]  
	+ (phi[(x+1)%L] + phi[(x-1+L)%L] );
      ResRoot += residue*residue; // true residue
    }

  return sqrt(ResRoot);    
}
