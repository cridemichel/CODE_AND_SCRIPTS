#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char**argv){

  FILE* fin;

  char name[1000];
  int Ntot; 
  int Nblock;

  int Npoint;

  int i,j,k;
  double x,y;
  double mx,my;
  double sx,sy;

  scanf("%s%d%d",name,&Ntot,&Nblock);
  Npoint = Ntot/Nblock;

  fin=fopen(name,"r");
  for(i=0;i<Npoint;i++){
    fscanf(fin,"%lf%lf",&mx,&my);
    sx=mx*mx; sy=my*my;
    for(j=1;j<Nblock;j++){
      fscanf(fin,"%lf%lf",&x,&y);
      mx += (x-mx)/(double)j;
      my += (y-my)/(double)j;
      sx += (x*x-sx)/(double)j; 
      sy += (y*y-sy)/(double)j;
    }
    printf("%lf %lf\n",mx,my);
  }
  fclose(fin);

  return 0;
}
