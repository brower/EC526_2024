#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double f(double x)
{
  return sin(x);
}

double derivative(double x)
{
  return cos(x);
}

double forward_diff(double x, double h)
{
  return (f(x+h)-f(x))/h;
}

double backward_diff(double x, double h)
{
 // return the backward difference.
  return 0.0; 
}

double central_diff(double x, double h)
{
  // return the central difference.
  return 0.0;
}

int main(int argc, char** argv)
{
  int i;
  double h;
  const double x = 1.0;

  cout << setprecision(15); // .. 10 points after the decimal.
  
  // Loop over 17 values of h: 1 to 10^(-16).
  h = 1.0;
  for (i = 0; i < 17; i++)
  {
    // Print stuff to graph
  }

  return 0;
}

