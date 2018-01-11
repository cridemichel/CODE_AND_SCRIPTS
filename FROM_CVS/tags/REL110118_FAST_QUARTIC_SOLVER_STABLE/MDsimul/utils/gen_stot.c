#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
const double Stot_a_C0 = 2917.8 - 5898.66 + 12525.1563 + 500.402383,
  T0=5.0, N = 1000.0; 
double stot(double T)
{
  return -1.5*3407.5*(pow(T,-0.4)-pow(T0,-0.4)) + 1.5*N*log(T/T0) + 
    Stot_a_C0;/* + 63.1245*log(T/T0);*/
}
double temp[] = {0.77,0.78,0.79,0.80,0.81,0.82,0.85,0.90,1.0,1.2,
  1.5,2.0,3.0,5.0,-1.0};
int main(int argc, char **argv)
{
  int i=0;
  while (temp[i]>0.0)
    {
      printf("%g %f\n", temp[i], stot(temp[i])); 
      i++;
    }
  return 0;
}
