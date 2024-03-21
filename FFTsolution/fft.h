#ifndef FFT_H
#define FFT_H

#include <iostream>
#include <fstream>


#include <cmath>
#include <complex>
using namespace std;
#define PI 3.141592653589793
#define  I  Complex(0.0, 1.0)
#define  Nfft 16
typedef complex <double> Complex;

void makePhase(Complex *omega, int N);
void FT(Complex * Ftilde, Complex * F, Complex * omega, int N);
void FTinv(Complex * F, Complex * Ftilde, Complex * omega, int N);
void FFT( Complex * Ftilde, Complex * F, Complex * omega, int N);
void FFTinv(Complex * F, Complex * Ftilde, Complex * omega, int N);
void BitRevArray(Complex * Frev, Complex * F, int N);
unsigned int  reverseBits(unsigned int x, int N);
void FFTrecursion(Complex* F, Complex * Ftilde, Complex * omega, int N);
void FFTrecursion_inv(Complex* Ftilde, Complex * F, Complex * omega, int N);

  
#endif 
