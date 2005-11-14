#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

double** dmatrix( int rows, int cols );
double*  dvector( int rows );



read_cristiano_ellipses( conf_t *conf, int allocate)
{

  char  line[500000];
  char  prefix[100];

  FILE *fp;

  int i,count;
/*   double rx[256],ry[256],rz[256]; */
/*   double R[256][9]; */
/*   double vx[256],vy[256],vz[256]; */
/*   double wx[256],wy[256],wz[256]; */
   
  fp = fopen(conf->nomefile,"r");
  if (fp == NULL) {
    fprintf(stderr, "<%s> open failure\n", conf->nomefile);
    exit(1);
  }


  // read away first block up to @@@
  do fgets(line, sizeof(line), fp);
  while(strncmp(line,"@@@",3)!=0);

  // read away second block up to @@@; get parameters
  do
    {
      fgets(line, sizeof(line), fp);
      if(strncmp(line,"parnum:",7)==0){
	sscanf(line,"%s %d",prefix,&conf->parnum);
	//printf("%s %d\n",prefix,conf->parnum);
      }

      if(strncmp(line,"parnumA:",8)==0){
	sscanf(line,"%s %d",prefix,&conf->parnumA);
	//printf("%s %d\n",prefix,conf->parnumA);
      }
      conf->parnumB = conf->parnum - conf->parnumA;

      if(strncmp(line,"time:",5)==0){
	sscanf(line,"%s %lf",prefix,&conf->time);
	//printf("%s %lf\n",prefix,conf->time);
      }

      if(strncmp(line,"T:",2)==0){
	sscanf(line,"%s %lf",prefix,&conf->kbT);
	conf->kbT /= 2.0; // half for transl, half for rotations...
	//printf("%s %lf\n",prefix,conf->kbT);
      }

      if(strncmp(line,"m:",2)==0){
	sscanf(line,"%s %lf %lf",prefix,&conf->mA,&conf->mB);
	//printf("%s %lf %lf\n",prefix,conf->mA,conf->mB);
      }

      if(strncmp(line,"a:",2)==0){
	sscanf(line,"%s %lf %lf",prefix,&conf->aA,&conf->aB);
	//printf("%s %lf %lf\n",prefix,conf->aA,conf->aB);
      }
     if(strncmp(line,"b:",2)==0){
	sscanf(line,"%s %lf %lf",prefix,&conf->bA,&conf->bB);
	//printf("%s %lf %lf\n",prefix,conf->bA,conf->bB);
      }
      if(strncmp(line,"c:",2)==0){
	sscanf(line,"%s %lf %lf",prefix,&conf->cA,&conf->cB);
	//printf("%s %lf %lf\n",prefix,conf->cA,conf->cB);
      }

      if(strncmp(line,"I:",2)==0){
	sscanf(line,"%s %lf %lf",prefix,&conf->IA,&conf->IB);
	//printf("%s %lf %lf\n",prefix,conf->IA,conf->IB);
      }

     
    }
  while(strncmp(line,"@@@",3)!=0);

  // eventually allocate position, velocities, rotations...
  if(allocate) AllocateConf(conf,conf->parnum);


  // read position rx,ry,rz and rotation matrix R
  for(i=0;i<conf->parnum;i++){
    if( fgets(line, sizeof(line), fp) ) 

      sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &conf->rx[i], &conf->ry[i], &conf->rz[i], 
	     &conf->R[i][0], &conf->R[i][1], &conf->R[i][2], 
	     &conf->R[i][3], &conf->R[i][4], &conf->R[i][5], 
	     &conf->R[i][6], &conf->R[i][7], &conf->R[i][8]
	     );

    else  {
      fprintf(stderr, "read x-y-z-R failure at line %d\n", i);
      exit(2);
    }

  }



  // read velocities vx,vy,vz and angular velocities wx,wy,wz
 for(i=0;i<conf->parnum;i++){
    if( fgets(line, sizeof(line), fp) ) 

      sscanf(line,"%lf %lf %lf %lf %lf %lf",
	     &conf->vx[i], &conf->vy[i], &conf->vz[i], 
	     &conf->wx[i], &conf->wy[i], &conf->wz[i]
	     );

    else  {
      fprintf(stderr, "read v-w failure at line %d\n", i);
      exit(3);
    }

  }


  fgets(line, sizeof(line), fp);
  sscanf(line,"%lf",&conf->box);

  if( fgets(line, sizeof(line), fp) ) {
    fprintf(stderr, "too many lines in this file !!!\n");
    exit(4);
  }

  fclose(fp);

}


AllocateConf(conf_t *conf,int parnum)
{ int i;

 conf->rx  = dvector(parnum);
 conf->ry  = dvector(parnum);
 conf->rz  = dvector(parnum);
 
 conf->R  = dmatrix(parnum,9);
 
 conf->vx = dvector(parnum);
 conf->vy = dvector(parnum);
 conf->vz = dvector(parnum);
 
 conf->wx = dvector(parnum);
 conf->wy = dvector(parnum);
 conf->wz = dvector(parnum);

}

#ifndef JUNK

double* dvector(int N){
  double* pt;
  pt=(double *) malloc(N*sizeof(double));
  if(pt==NULL){
    fprintf(stderr,"Error in allocating vector\n");
    exit(1);
  }
  return pt;
}

double** dmatrix(int N, int M){
  double **pt;
  int i;

  pt=(double**) malloc(N*sizeof(double*));
  if(pt==NULL){
    fprintf(stderr,"Error in allocating matrix\n");
    exit(2);
  }

  for(i=0;i<N;i++){
    pt[i] = (double*) malloc(M*sizeof(double));
    if(pt[i]==NULL){
      fprintf(stderr,"Error in allocating row n. %d\n",i);
      exit(3);
    }
  }

  return pt;
}

#endif

#ifdef JUNK

double *dvector( int rows )
{
   int      i;
   double   *v;


   /* Allocate memory for pointers to rows */
   v = (double *) malloc( rows * sizeof(double) );
   if (!v)
      return( NULL );

   /* Return pointer to rows */
   return( v );
}

double **dmatrix( int rows, int cols )
{
   int      i;
   double   **m;


   /* Allocate memory for pointers to rows */
   m = (double **) malloc( rows * sizeof(double*) );
   if (!m)
      return( NULL );

   /* Pointer to first row allocates memory for box */
   m[0] = (double *) malloc( rows * cols * sizeof(double) );
   if (!m[0])
       return( NULL );

   /* Set pointers to rows */
   for (i = 1; i < rows; i++)
     m[i] = m[i-1] + cols;

   /* Return pointer to array of pointers to rows */
   return( m );
}


#endif
