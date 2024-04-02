#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>

int Nfft = 64;

using namespace std;

#include "fft.h"

int main()
{
  int N;
  N = Nfft;
  Complex * omega = new Complex[N];
  Complex * F = new Complex[N];
  Complex * Ftilde  = new Complex[N];
  Complex * Fnew = new Complex[N];
  Complex * Fsave = new Complex[N];
  
  makePhase(omega,N);
   cout << std::fixed;
   cout << setprecision(3);
  
  /************  Test SLOW  FT ****************/
  
  for(int x = 0; x < N; x++){
    F[x] = 2* sin( 2.0*PI*x/(double) N) + 4*cos( 2.0*PI*3.0*x/(double) N) - 8*cos( 2.0*PI*4.0*x/(double) N) ;
    Fsave[x] = F[x];
    Ftilde[x] = 0.0;
    Fnew[x] = 0.0;
    cout<<" x = "<< x << "  F =    " <<  F[x]  << endl;
  }
  
  /** MY CONVENTION IS THE OUT ARRAY IS ON THE LEFT and IN NEXT TO IT **/
  cout << endl << "Test  Slow FT " << endl;
  FT(Ftilde, F,omega,N);

 
  for(int k = 0; k < N; k++)
    cout<<" k = "<< k << "  Ftilde =   " <<  Ftilde[k]  << endl;
  
  FTinv(Fnew, Ftilde, omega, N);
  
  for(int x = 0; x < N; x++)
    cout<<  " x = " << setw(4) << x << " Compare newF[x] : " << setw(25) << Fnew[x] <<" with F[x]  "<< F[x] << endl;
  
  /*************    TEST ITERATIVE  FFT **************/
  
  cout << endl << " Test Iterative   FFT " << endl; 
  
  for(int x = 0; x < N; x++){
    F[x] = Fsave[x];
    Ftilde[x] = 0.0;
    Fnew[x] = 0.0;
  }
  
  FFT( Ftilde,  F, omega,  N);
  
  for(int k = 0; k < N; k++)
    cout<<" k = "<< k << "  Ftilde =   " <<  Ftilde[k]  << endl;
  
  FFTinv(Fnew,  Ftilde, omega,  N);

  
  for(int x = 0; x < N; x++)
    cout<<  " x = " << setw(4) << x << " Compare Fnew[x] : " << setw(25) << Fnew[x] <<" with F[x]  "<< F[x] << endl;
  
  /*******************     TEST ITERATIVE  FFT **********************/
  
  cout << endl << " Test  FFTrecursion " << endl; 
  
  for(int x = 0; x < N; x++){
    F[x] = Fsave[x];
    Ftilde[x] = 0.0;
    Fnew[x] = 0.0;
  }

  cout << " Calls at N = " << N << endl;
  FFTrecursion( Ftilde,  F, omega,  N);

  for(int k = 0; k < N; k++)
    cout<<" k = "<< k << "  Ftilde =   " <<  Ftilde[k]  << endl;
  
  FFTrecursion_inv(Fnew,  Ftilde, omega,  N);

  
  for(int x = 0; x < N; x++)
    cout<<  " x = " << setw(4) << x << " Compare Fnew[x] : " << setw(25) << Fnew[x] <<" with F[x]  "<< F[x] << endl;
  
  return 0;
}

