#include <iostream>
#include <cmath>

using namespace std;


int getLegendreCoeff(double* a, int n);

double evalPolynomial(double x, double* a, int n);

double evalPolyDeriv(double x, double* a, int n);

int getLegendreZero(double* zero, double* a, int n);



int main(int argc, char** argv)
{

  // Fix the output format.
  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(20);

  // Order of Legendre polynomial
  int order = 0;

  // Space to hold coefficients
  double* a = 0;

  // Describe the problem and prompt the user:
  cout << "What order Legendre polynomial? ";
  cin >> order;

  if (order < 0) {
    cout << "Error! The Legendre polynomial needs to be of positive integer order.\n";
    return -1;
  }

  // Allocate the coefficients
  a = new double[order+1];
  for(int i = 0; i < order+1;i++) a[i] = 0.0;

  // Get the coefficients
  if (getLegendreCoeff(a, order) == -1) {
    cout << "Error! The coefficient array wasn't allocated.\n";
    return -1;
  }

  // Print
  cout << "Coefficients\n"; 
  for (int i = order; i > 0; i--) {
    cout << a[i] << " x^" << i << " +  \n";

  }
  cout << a[0] << " x^0\n";

  // Compute the zeros
  double* zeros = new double[order];
  if (getLegendreZero(zeros, a, order) == -1) {
    cout << "Error! The zeros array wasn't allocated.\n";
    delete[] a;
    return -1;
  }

  // Print
  cout << "Zeros\n";
  for (int i = 0; i < order; i++) {
    cout << zeros[i] << " ";
  }
   cout << "\n";

  // Clean up
  delete[] a;

}

// P_n has  n + 1  coefficents coeffs[n][0], coeffs[n][1], ..,  coeffs[n][n]
// This function's wasteful because every time you run it
// you have to re-compute the entire recursive stack.
// But whatever.
int getLegendreCoeff(double* a, int n)
{
  if (a == 0) { return -1; } // error out if not allocated
  
  // The startup cases
  if (n == 0) { a[0] = 1.0; return 0; }
  if (n == 1) { a[0] = 0.0; a[1] = 1.0; return 0; }

  // Allocate space for each set of coefficients
  //double** coeffs;

  // Allocate things.
   double** coeffs = new double*[n+1];

  // Fill out the startup cases
  coeffs[0] = new double[1];
  coeffs[0][0] = 1.0;
  coeffs[1] = new double[2];
  coeffs[1][0] = 0.0;
  coeffs[1][1] = 1.0;

  // And go go go
  for (int m = 2; m <= n; m++)
  {
    coeffs[m] = new double[m+1];

    // Start filling up the coefficients.
    
  }

  // Copy the last row in
  for (int i = 0; i <= n; i++)
  {
    a[i] = coeffs[n][i];
  }

  // Clean up
  for (int i = 0; i <= n; i++)
  {
    delete[] coeffs[i];
  }
  delete[] coeffs;

  return 0;

}

// evaluate a polynomial of order 'n'
double evalPolynomial(double x, double* a, int n)
{
  if (n == 0) { return a[0]; }
  double sum = a[0];
 

  return sum;
}

// evaluate the derivative of a polynomial of order 'n'
double evalPolyDeriv(double x, double* a, int n)
{
  if (n == 0) { return 0; }
  if (n == 1) { return a[1]; }
  double sum = a[1];


  return sum;
}

int getLegendreZero(double* zero, double* a, int n)
{
  if (zero == 0) { return -1; } // error out if not allocated
  if (a == 0) { return -1; } // error out if not allocated

  for (int i = 1; i <= n; i++)
  {
  
  }

  return 0;

}

