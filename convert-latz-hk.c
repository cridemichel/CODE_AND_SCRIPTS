#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif

// format of hk: q 000 200 220 400 420 440 221 421 441 222 422 442 443 444
// where l1 l2 m as m1=m2
// in our language: 
// 1   2    3    4    5    6    7    8    9   10   11   12   13   14   15
// q 0000 2000 2020 4000 4020 4040 2121 4121 4141 2222 4222 4242 4343 4444

double s[14];
int l1[14]={0,2,2,4,4,4,2,4,4,2,4,4,4,4};
int l2[14]={0,0,2,0,2,4,2,2,4,2,2,4,4,4};
int  m[14]={0,0,0,0,0,0,1,1,1,2,2,2,3,4};

int main(int argc, char **argv){

  
  double unit_length;
  double fac_q;
  double fac_Ylm= 1.0;

  char *hk_name=argv[1];
  FILE *hk_file;
  double Phi=atof(argv[2]);
  double X0=atof(argv[3]);
  double rho;

  char line[1000];
  double Q;
  int i,j,k;

  double convert_fs;


  fac_q = pow( X0 , 1./3. );
  rho=6.0*Phi/(M_PI*X0);
  fac_Ylm = rho/(4.0*M_PI);
  hk_file=fopen(hk_name,"r");

#if 0
  convert_fs=pow(3.*Phi/(4.*M_PI),1./3.);
  fac_q /= convert_fs;
#endif


  for(i=0;i<4;i++) fscanf(hk_file,"%s",line);
  while( fgets(line, sizeof(line), hk_file) ){
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	   &Q,&s[0],&s[1],&s[2],&s[3],&s[4],&s[5],&s[6],&s[7],&s[8],
	   &s[9],&s[10],&s[11],&s[12],&s[13]);    

    Q     *= fac_q; 

    if(fabs(Q)>1.0e-10) {

      printf("%lf ",Q);
      for(j=0;j<14;j++){
	s[j] *= fac_Ylm;
#if 0
	if((l1[j]>0)||(l2[j]>0)) s[j]/=M_PI;
#endif
	if((m[j]%2)==1)  s[j] = -s[j];
	if(l1[j]==l2[j]) s[j] += 1.0;
      printf("%lf ",s[j]);
      }
      putchar('\n');

    }//ifQ

  }//while

  fclose(hk_file);

  return 0;
}
