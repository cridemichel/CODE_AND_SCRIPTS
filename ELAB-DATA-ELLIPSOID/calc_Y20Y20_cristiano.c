#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define MY_TRUE  (1==1)
#define MY_FALSE (1==0)

main(){

  conf_t conf;
  double Y20,Y20Y20;
 
  
#ifndef Npart
#define Npart 256
#endif
  
  
  int  TotalConfs,iConf;
  char nomedir[1000];
  char nome[1000];
  char nomeconf[1000];
  int  allocate;
  double box,invbox;
  
  int i;
  
  double unit_length, unit_time, fac_r2, fac_time;
  
  double cosTheta;
  
  /* initialize averages to zero*/
  Y20Y20 = 0.0;
  
  
  scanf("%d",&TotalConfs); // number of configurations
  scanf("%s",nomedir); // directory where configurations are         
  allocate = MY_TRUE;  // first Blocksize confs must be allocated    
  
  /* I will average each of the BlockSize points over
   * TotalConfs/BlockSize values */		      
  for(iConf=0; iConf<TotalConfs; iConf++){
    
    /* reads the configuration */
    strcpy(conf.nomefile,nomedir); 
    scanf("%s",nomeconf); 
    strcat(conf.nomefile,nomeconf);
    read_cristiano_ellipses(&conf,allocate);
    if(conf.parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}
    
    if(allocate){ // if it is the first time, calculate parameters
      
      /* assume that the size of the box does not change */
      box   = conf.box; 
      
      /* calculate the effective diameter */
      unit_length = 2.0*pow(conf.aA*conf.bA*conf.cA,1./3.);
      
      /* calculate rescaling factor to get things adimensional */
      fac_r2 = 1./(unit_length*unit_length);
      unit_time = sqrt(conf.mA/conf.kbT)*unit_length; 
      fac_time = 1.0/(unit_time*unit_time);
      
      allocate = MY_FALSE;

    }//end calculate parameters
    

    Y20 = 0.0;
    for(i=0;i<Npart;i++) {
      cosTheta = conf.R[i][0]; 
      Y20+=( 3.*cosTheta*cosTheta - 1.) / 2.;
    }
    Y20Y20 += Y20*Y20;
      	
    Y20 = 0.0;
    for(i=0;i<Npart;i++) {
      cosTheta = conf.R[i][1]; 
      Y20+=( 3.*cosTheta*cosTheta - 1.) / 2.;
    }
    Y20Y20 += Y20*Y20;
      	
    Y20 = 0.0;
    for(i=0;i<Npart;i++) {
      cosTheta = conf.R[i][2]; 
      Y20+=( 3.*cosTheta*cosTheta - 1.) / 2.;
    }

    Y20Y20 += Y20*Y20;
      	
    allocate = MY_FALSE; 
    
  }/* iConf - end of reading all */
  
  /* normalize */

  Y20Y20  /= 3.0 * (double)Npart * (double)TotalConfs;

  /* spherical harmonics prefactor */

  Y20Y20  *= 5.0 / (4.0 * M_PI);

  /* rho-rho prefactor */
  Y20Y20  *= 4.0 * M_PI;

  printf("%lf\n",Y20Y20 );
  

  //  return Y20Y20;

}

