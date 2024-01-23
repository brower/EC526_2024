#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/****************************
Testing accuracy
  float  pi = 3.141592653589793238462643383279502884;  // 8 
  double PI = 3.141592653589793238462643383279502884;  // 16 places
  printf("pi = %10.200f  PI = %5.40f \n", pi, PI);
  **************************/
int main()
{
  float a, b,c;
  a = 1.0/3.0;
  b =  1.0000001;
  c = -1.0000000;

  cout << "Single Results \n";
   printf("  a*(b+c) =  %10.30f  \n",  a*(b+c));
  printf("  a*b+a*c =  %10.30f  \n",  a*b + a*c);
  printf( " a*(b+c) -   (a*b+a*c)  =  %10.30f \n",    a*(b+c) -   (a*b+a*c) );
 
  double A, B,C;
  A = 1.0/3.0;
  B =  1.0000001;
  C = -1.0000000;


  cout << "Double  Results \n";
  printf("  A*(B+C) =  %10.30f  \n",  A*(B+C));
  printf("  A*B+A*C =  %10.30f  \n",  A*B + A*C);
  printf( " A*(B+C) -   A*B+A*C  =  %10.30f \n",   A*B + A*C - (A*B + A*C));

  //cast back to single
 
  float error =  A*(B+C)  - ( A*B+A*C);

  cout<< " Intermediate Double back to single \n";
  
  printf("cast to single  A*(B+C) -   A*B+A*C = %10.30f \n" , error);
  

  
  return 0;
}

// float  pi = 3.141592653589793238462643383279502884;  // 8 
  //  double PI = 3.141592653589793238462643383279502884;  // 16 places
  // printf("pi = %10.200f  PI = %5.40f \n", pi, PI);
