#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define mytrue  (1==1)
#define myfalse (1==0)



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

  int i,j; 

  // read where is the input conf, which is its name
  scanf("%s %s",nome_dir_in,nome_conf_in);
  // read integer parameters l,m for the spherical harmonic Y_{l,m}


  // full name of the input conf
  sprintf(conf.nomefile,"%s/%s",nome_dir_in,nome_conf_in);
  read_cristiano_ellipses(&conf,allocate);
    

  if(conf.parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}

  box   = conf.box; 

#ifdef SCHILLING
  unit_length = 2.0*conf.aA;
#else
  unit_length = pow(2.0*conf.aA*conf.bA*conf.cA,1./3.);
#endif
  
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

  printf("%lf %lf %5.3lf %5.3lf\n",max_len,min_len,2.*M_PI/max_len,2.*M_PI/min_len);


  return 0;

}


