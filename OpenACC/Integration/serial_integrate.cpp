#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>  // for clock_t, clock(), CLOCKS_PER_SEC

/*** 
https://www.techiedelight.com/find-execution-time-c-program/ 
***/

int main(int argc, char** argv){

  // OpenACC initialise 
  //  #pragma acc init

  int N=1000000000;
  double lower = -100;
  double upper = 100;
  double interval = upper - lower;
  double increment = interval/N;  
  double sum = 0.0;

  // to store execution time of code
  double time_spent = 0.0;

  
  //Integrate from -pi to pi.
  //#pragma acc parallel loop reduction(+:sum)

 
  clock_t begin = clock();

  // #pragma acc parallel loop 
  for(int dx=0; dx<N; dx++) {
    
    //The x-axis value
    double x = lower + increment*(dx + 1.0/2);
    sum += increment * exp(-x*x);
    
  }

  clock_t end = clock();

  time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
  printf("Time elpased is %f seconds\n", time_spent);
  
  printf("integral of exp(-x*x) from x=%.2e to %.2e = %.16f, error = %e%%.\n", lower, upper, sum, 100*(sum - sqrt(M_PI))/sqrt(M_PI));
  
  return 0;
}
