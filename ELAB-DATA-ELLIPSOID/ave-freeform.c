#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char**argv){

  int N = 0; 
  double x,y;
  double mx = 0., my = 0.;
  double sx = 0., sy = 0.;

  while ( scanf("%lf%lf",&x,&y) != EOF ){
    N++;
    mx += (x-mx)/(double)N;
    my += (y-my)/(double)N;
    sx += (x*x-sx)/(double)N; 
    sy += (y*y-sy)/(double)N;
  }
  
  printf("%lf %lf\n",mx,my);

  return 0;
}
