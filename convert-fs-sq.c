#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

int main(int argc, char **argv){
char namefile[1000];
FILE *fp;
double Phi;
int l1,l2,m1,m2;
double deltal1l2;


char line[1000];
double q,sq,err;
double convert_q;
double convert_s;

	sprintf(namefile,"%s",argv[1]);
	fp=fopen(namefile,"r");
	Phi=atof(argv[2]);
	l1=atoi(argv[3]); m1=atoi(argv[4]);
	l2=atoi(argv[5]); m2=atoi(argv[6]);
	
	if(l1==l2) deltal1l2=1.0; 
	else deltal1l2=0.0;

	convert_q=pow(3.*Phi/(4.*M_PI),1./3.);
	if((l1==0)&&(l2==0)) convert_s=1.0;
	else convert_s=M_PI;

	//	fprintf(stderr,"%d %d %d %d : delta = %lf\n",l1,l2,m1,m2,deltal1l2);
		
	while(fscanf(fp,"%lf%lf%lf",&q,&sq,&err)!=EOF)
		printf("%g %g\n", convert_q*q,
		       deltal1l2+convert_s*(sq-deltal1l2));

	fclose(fp);
	
	return 0;
}
