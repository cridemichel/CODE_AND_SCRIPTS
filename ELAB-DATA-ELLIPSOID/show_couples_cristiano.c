#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define mytrue  (1==1)
#define myfalse (1==0)

int icolor;
char* color[10] = {"red","green","blue","yellow","orange","violet","cyan","magenta","pink","grey"};


int main(int argc, char**argv){

  conf_t conf;
  int allocate = mytrue;

#ifndef Npart
#define Npart 256
#endif

  char  nome_dir_in[1000];
  char  nome_conf_in[1000];


  double box,invbox;
  
  double unit_length, unit_time, fac_r2, fac_time;

  double min_len,max_len;

  int i,j,k; 

  double thetaerr = 1.0e-2*M_PI/180.; // uncertainity in theta 
  double dmin=0.0,dmax=10000;

  double dist,dx,dy,dz,dot_prod;
  double ave_axis[3];
  double soglia=-1000;
  double rcut=2.3;

  // full name of the input conf
  sprintf(conf.nomefile,"%s",argv[1]);
  read_cristiano_ellipses(&conf,allocate);
  
  if(argc>2) { thetaerr = atof(argv[2]); thetaerr*=M_PI/180.;} 
  if(argc>3) dmin = atof(argv[3]); 
  if(argc>4) dmax = atof(argv[4]); 

  if(conf.parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}

  box   = conf.box; 
  printf(".Vol: %lf\n",conf.box*conf.box*conf.box);



  unit_length = 2*pow(conf.aA*conf.bA*conf.cA,1./3.);
  
  fac_r2 = 1./(unit_length*unit_length);
  unit_time = sqrt(conf.kbT/conf.mA)/unit_length;

  //  q_lenght=pow((double)Npart,1./3.)*2.*M_PI/box;
  //  q_lenght=4.*M_PI/box;
  //  printf("%lf %lf\n",pow((double)Npart,1./3.)*2.*M_PI/box,4.*M_PI/box);
  //  exit(0);

  max_len=conf.box/2.; 

  min_len=conf.aA; 
#ifdef MAXMIN
  if(conf.bA>min_len) min_len=conf.bA;
  if(conf.cA>min_len) min_len=conf.cA;
#else
  if(conf.bA<min_len) min_len=conf.bA;
  if(conf.cA<min_len) min_len=conf.cA;
#endif

  fprintf(stderr,"%g %g %g %g %g\n",unit_length,conf.aA,conf.bA,conf.cA);

  icolor=0;
  for(i=1;i<Npart;i++) for(j=0;j<i;j++) {

    dx = conf.rx[i] - conf.rx[j];
    dy = conf.ry[i] - conf.ry[j];
    dz = conf.rz[i] - conf.rz[j];
    /* adjust XYZ displacements so that it lies inside the box */
    /*
      dx=remainder(dx,box);	
      dy=remainder(dy,box);
      dz=remainder(dz,box);
    */
    dist = sqrt( dx*dx + dy*dy + dz*dz );

    if( ( fabs(dz/dist) > cos(thetaerr) )  &&
	(dist < dmax) &&
	(dist > dmin)    ){      

      printf("%lf %lf %lf ",conf.rx[i],conf.ry[i],conf.rz[i]);
      for(k=0;k<9;k++) printf("%lf ",conf.R[i][k]);
      printf("@ %lf %lf %lf C[%s]",conf.aA,conf.bA,conf.cA,color[icolor]);
      putchar('\n');
      printf("%lf %lf %lf ",conf.rx[j],conf.ry[j],conf.rz[j]);
      for(k=0;k<9;k++) printf("%lf ",conf.R[j][k]);
      printf("@ %lf %lf %lf C[%s]",conf.aA,conf.bA,conf.cA,color[icolor]);
      putchar('\n');
      icolor=(icolor+1)%10;

      // plot ball with effective volume
#if 0
      printf("%lf %lf %lf ",conf.rx[i],conf.ry[i],conf.rz[i]);
      for(k=0;k<9;k++) printf("%lf ",conf.R[i][k]);
      printf("@ %lf %lf %lf C[%s]",unit_length/2.,unit_length/2.,unit_length/2.,"grey");
      putchar('\n');
      printf("%lf %lf %lf ",conf.rx[j],conf.ry[j],conf.rz[j]);
      for(k=0;k<9;k++) printf("%lf ",conf.R[j][k]);
      printf("@ %lf %lf %lf C[%s]",unit_length/2.,unit_length/2.,unit_length/2.,"grey");
      putchar('\n');
#endif

  	
    }//if
  }//for i>j


  return 0;

  }


