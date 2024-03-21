#include "fft.h"

void makePhase(Complex *omega, int N )
{
  for(int k = 0; k < N; k++)
    omega[k] = exp(2.0*PI*I*(double)k/(double) N);
}

void FT(Complex * Fout, Complex * Fin, Complex * omega, int N)
{
  for(int k = 0; k < N; k++)
    {
      Fout[k] = 0.0;
      for(int x = 0; x < N; x++)
      	Fout[k] += omega[(k*x)%N]*Fin[x];
	
    }
}

void FTinv(Complex * Fout, Complex * Fin, Complex * omega, int N)
{
  for(int x = 0; x < N; x++)
    {
      Fout[x] = 0.0;
      for(int k = 0; k < N; k++)
	Fout[x] +=conj(omega[((k*x))%N])*Fin[k]/(double) N;
    }
}



void BitRevArray(Complex * Frev, Complex * F, int N)
{
  for(int x = 0;x < N; x++)
     F[reverseBits(x, N)] = F[x];
}

unsigned int  reverseBits(unsigned int x, int N) {
        unsigned int xrev = 0;
	unsigned int p = log2(N) ; // p = 4
        unsigned int n;
	unsigned int power = N;
        for(int i = 0;i<p;i++)
	  {
	    n = x%2;  // find lowest bit
	    power /= 2;
 	    xrev += n*power; //  add to highest 2^3
	    x /= 2;
	  }
	return xrev;	  
}


/*********************************************************************

Butterfly Iterative Solution coarsenting to large blocks
Bsize =   2, 4, 8  ...; index  b = 0,..,Bsize- 1  
Nblock  =  N/Bsize  index  nb =0,.., Nblock -1  

******************************************************************/

void FFT( Complex * Fout, Complex * Fin, Complex * omega, int N)
{
  int log2n = log2(N);
  Complex high, low;
  int Bsize = 0;
  int Nblocks = 0;
  
  for(unsigned int x = 0; x < N; x++)
    Fout[reverseBits( x, N)] = Fin[x];

  for(int s = 0; s < log2n; s++) 
    {
      Bsize = pow(2,s+1);                         
      Nblocks  = N/Bsize;                         
        
      for(int nb = 0;nb < Nblocks ; nb++)                      // slow to next Nblocks
	for(int pairs  = 0; pairs <  Bsize/2; pairs++)         // fast index over nomber of  Pairs = Bsize/2
	  {
	    int  k = pairs +  nb * Bsize;                       //low pairs index : {x high = low + Bsize/2
	    low = Fout[k];
	    high =  omega[Nblocks*pairs]*Fout[k +  Bsize/2];    //(omega_Bsize)^pairs = (omega_N^(Nblocks*pairs)
	    Fout[k] = low + high;
	    Fout[k + Bsize/2]  = low -  high;
	  } 	 
    }
}

void FFTinv( Complex * Fout, Complex * Fin, Complex * omega, int N)
{
  int log2n = log2(N);
  Complex high, low;
  int Bsize = 0;
  int Nblocks = 0;
  
  for(unsigned int x = 0; x < N; x++)
    Fout[reverseBits( x, N)] = Fin[x]/(double)N;

  for(int s = 0; s < log2n; s++) 
    {
      Bsize = pow(2,s+1);                         
      Nblocks  = N/Bsize;                         
        
      for(int nb = 0;nb < Nblocks ; nb++)                          // slow to next Nblocks
	for(int pairs  = 0; pairs <  Bsize/2; pairs++)             // fast index over nomber of  Pairs = Bsize/2
	  {
	    int  k = pairs +  nb * Bsize;                          //low pairs index : {x high = low + Bsize/2
	    low = Fout[k];
	    high =  conj(omega[Nblocks*pairs])*Fout[k +  Bsize/2]; //(omega_Bsize)^pairs = (omega_N^(Nblocks*pairs)
	    Fout[k] = low + high;
	    Fout[k + Bsize/2]  = low -  high;
	  } 	 
    }
}

void FFTrecursion(Complex* Fout, Complex * Fin, Complex * omega, int N)
{
  
  cout << "Nfft =  " << Nfft << endl;
   cout << "Enters at  N " << N << endl;
  if(N == 2)
    {
      Fout[0] = Fin[0]  + Fin[1];
      Fout[1] = Fin[0] - Fin[1];
      return;
    }
  
  Complex* Fin_even = new Complex[N/2];
  Complex* Fout_even = new Complex[N/2];
  Complex* Fin_odd = new Complex[N/2];
  Complex* Fout_odd = new Complex[N/2];
  
  for(int i = 0; i < N/2; i++){
    Fin_even[i] = Fin[2*i];
    Fin_odd[i] = Fin[2*i + 1];
  }

  cout << "Calls at twice at N =  "<< N/2 << endl;
  FFTrecursion(Fout_even,Fin_even,omega,N/2);
  FFTrecursion(Fout_odd,Fin_odd,omega,N/2);
 
 int lev  = Nfft/N;
 cout << " at level =  " << lev << " Combine two  N =  " << N/2 << endl;
 
 for(int i = 0; i < N/2; i++){        // N = 4  need omega_4(i) = (omega_16)^4
    Fout[i] = Fout_even[i] + omega[i*lev]*Fout_odd[i]; 
    Fout[i + N/2] =  Fout_even[i] - omega[i*lev]*Fout_odd[i];
  }
  
}

void FFTrecursion_inv(Complex* Fout, Complex * Fin, Complex * omega, int N)
{
  //  cout << "Inters at  N " << N << endl;
  if(N == 2)
    {
      Fout[0] = (Fin[0] + Fin[1])/(double)Nfft;
      Fout[1] = (Fin[0] - Fin[1])/(double)Nfft;
      return;
    }
  
  Complex* Fin_even = new Complex[N/2];
  Complex* Fout_even = new Complex[N/2];
  Complex* Fin_odd = new Complex[N/2];
  Complex* Fout_odd = new Complex[N/2];
  
  for(int i = 0; i < N/2; i++){
    Fin_even[i] = Fin[2*i];
    Fin_odd[i] = Fin[2*i + 1];
  }

  // cout << "Calls at twice at N =  "<< N/2 << endl;
  FFTrecursion_inv(Fout_even,Fin_even,omega,N/2);
  FFTrecursion_inv(Fout_odd,Fin_odd,omega,N/2);
 
 int lev  = Nfft/N;
 //cout << " Combine two  N =  " << N/2 << " at level =  " << lev << endl;
 
  for(int i = 0; i < N/2; i++){
    Fout[i] = Fout_even[i] + conj(omega[i*lev])*Fout_odd[i]; 
    Fout[i + N/2] =  Fout_even[i] - conj(omega[i*lev])*Fout_odd[i];
  }
}
