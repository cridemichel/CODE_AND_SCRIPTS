#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int Nq;
int Nave;
double* q;
double* Sq;
char nameSq[2000];
FILE* fileSq;


/* Reads Nave=number of files, Nq=numbers of q's in a file;
 * afterwards reads the Nave names of the files to open */

main(){
  int i,j;
  double x,y;

  scanf("%d%d",&Nave,&Nq);
  q=(double*)calloc(Nq,sizeof(double));
  Sq=(double*)calloc(Nq,sizeof(double));

  for(i=0;i<Nave;i++){
    scanf("%s",nameSq); 
    fileSq=fopen(nameSq,"r");
    for(j=0;j<Nq;j++){
      fscanf(fileSq,"%lf%lf",&x,&y);
      q[j]+=x; Sq[j]+=y;
    }
  }

  for(i=0;i<Nq;i++)
    printf("%lf %lf\n",q[i]/Nave,Sq[i]/Nave);

  fclose(fileSq);
}
