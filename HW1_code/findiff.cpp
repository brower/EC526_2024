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
  double x = 1.0;

  // USE PRINTF WHICH HAS BETTER FORMAT CONTROL:
  
  // Format to Plot Data   
  printf("#     h                fwd                bwd                 central               fwd Error             bwd Error          ctl Error \n");
  
  for (h = 1; h >= pow(10, -16); h *= 0.1)  // Loops over 17 decreasing valuse of h
    {
      // Print h, the forward, backward, and central difference of sin(x) at 1 as well as the exact derivative cos(x) at 1.
      double fwd = forward_diff(x, h);
      double bck = backward_diff(x, h);
      double ctl = central_diff(x, h);
      double fwd_error = abs(fwd - derivative(x))/abs(derivative(x));
      double  bck_error = abs(bck - derivative(x))/abs(derivative(x));
      double ctl_error  = abs(ctl - derivative(x))/abs(derivative(x));	
      printf("%20.16f  %20.16f  %20.16f  %20.16f %20.16f %20.16f %20.16f \n", h, fwd, bck,ctl , fwd_error, bck_error, ctl_error );
    }
  
  /************** HELP WITH GNUPLOT  IF YOU USE IT ******************
plot "diff.dat" u (log($1)):(log($5))  
gnuplot> replot "diff.dat" u (log($1)):(log($6))  
gnuplot> replot "diff.dat" u (log($1)):(log($7))  
  ********************************************************/
  


  return 0;
}

