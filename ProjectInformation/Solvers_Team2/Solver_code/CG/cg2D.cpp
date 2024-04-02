#include <iostream>
#include <math.h>
#include <stack>
using namespace std;

#define Lx  8// 16 
#define Ly  8 // 16
#define N  Lx*Ly



#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

inline int mod(int x, int n)
{
  return  (n + (x % n))%n;
}


inline int nn(int site,int mu)
{
  int x,y;
  int xp, xm, yp, ym;
  int neighbor;
  
  x = site%Lx; y = site/Lx;
  xp = mod(x+1,Lx); xm =  mod(x-1,Lx); yp  = mod(y+1,Ly); ym =  mod(y-1,Ly);
  
  switch(mu){
  case 0:  neighbor =  xp + Lx*y;break; // x + 1 
  case 1:  neighbor =  x + Lx*yp;  break; // y + 1
  case 2:  neighbor =  xm + Lx*y; break;  // x -1  
  case 3:  neighbor =  x + Lx*ym; break;  // y - 1 
 
  default: neighbor = -1;
  }
  return neighbor;
}
double  CGiter(double * x, double * b, int Niter);
int Avec(double * vec_out, double * vec_in);
double dotA(double * v1, double * v2);
double dot(double * v1,double * v2);
void printArray(double * phi);

int main()
{

  double  phi[N], b[N];

  srand(137);
  // Sample Spin random Configuration
  for(int i = 0; i < N ; i++) {
    rand()%2 == 0 ? phi[i] = 1.0 : phi[i] = -1.0;
    b[i] = 0.0;
  }

  b[Lx/2 + Lx* Ly/2] = 10;

    printArray(phi);
      printArray(b);
    
  Avec(phi, b);
  printArray(phi);
1
   CGiter(phi,  b, 1);

  printArray(phi);
   
  return 0;
}

int Avec(double * vec_out, double * vec_in)
{
  for(int i = 0; i < N; i++)
    {
      vec_out[i] = vec_in[i];
	for(int mu = 0; mu <4; mu++)
	  {
	    vec_out[i] -= 0.25 * vec_in[nn(i, mu)];
	  }
    }   
  return 0;
}

double dot(double *  v1, double  * v2)
{
  double scalar = 0.0;
  
  for(int site = 0; site < N; site++)
    scalar += v1[site]*v2[site];
  
  return scalar;
}



double dotA(double * Av1, double * v2)
{
  double scalar = 0.0;

Avec(Av1, v2);

for(int site = 0; site < N; site++)
  scalar += v2[site]*Av1[site];

return scalar;
}

double  CGiter(double * x, double * b, int Niter)
{
  double r[N], p[N];
  double Ax[N], Ap[N];
  double alpha;
  double rsold, rsnew;
    
  Avec(Ax,x);
  
  for(int i =0; i < N; i++)
    {
      r[i] = b[i] - Ax[i];
      p[i] =r[i];
    };
  
  for(int iter =0; iter < Niter; iter++)
    {  Avec(Ap, p);
      alpha = rsold / dotA(p,Ap);
      for(int i =0 ; i < N ; i++)
	{
	  x[i] = x[i] + alpha * p[i];
	  r[i] = r[i] - alpha * Ap[i];
	};
      rsnew = dot(r,r);
      if (sqrt(rsnew) < 0.0000001) break;
      for(int i =0 ; i < N ; i++)
	{  
	p[i] = r[i] + (rsnew/rsold)* p[i];
	};
      
      rsold = rsnew;
    }

  return rsnew;
  
}


									    
  
void printArray(double * phi)
{
	cout<<"\n--------------------------------------------";
	for(int y = 0; y<Ly; y++)
	  {
	     cout << endl;
	    for(int x= 0 ; x<Lx; x++)
   	      printf(" %10.5f ", phi[x + y* Lx]);    
	  }
	cout<<"\n-------------------------------------------- \n";
}
