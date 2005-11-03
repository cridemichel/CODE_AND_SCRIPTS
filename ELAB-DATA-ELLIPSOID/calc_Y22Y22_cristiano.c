#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define MY_TRUE  (1==1)
#define MY_FALSE (1==0)

main(){

  conf_t conf;
  double Y22re, Y22im;
  double Y22Y22;
 
  
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
  
  double cosTheta,sinsqTheta;
  double sinPhi,cosPhi;
  double sin2Phi,cos2Phi;

  /* initialize averages to zero*/
  Y22Y22 = 0.0;
  
  
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
    
    Y22re = Y22im = 0.0;
    for(i=0;i<Npart;i++) {
      cosTheta = conf.R[i][2];  sinsqTheta = 1 - cosTheta*cosTheta;
      cosPhi =  conf.R[i][0]; sinPhi = conf.R[i][1];
      sin2Phi=2.0*sinPhi*cosPhi; cos2Phi=cosPhi*cosPhi-sinPhi*sinPhi;
      Y22re += sinsqTheta*cos2Phi; Y22im += sinsqTheta*sin2Phi; 
    }
    Y22Y22 += Y22re*Y22re +Y22im*Y22im;
      	
       
  }/* iConf - end of reading all */
  
  /* normalize */

  //  Y22Y22  /= 3.0 * (double)Npart * (double)TotalConfs;
  Y22Y22  /= (double)Npart * (double)TotalConfs;

  /* spherical harmonics prefactor */

  Y22Y22  *= ( 15.0 / (2.0 * M_PI) ) / (4.0 *4.0) ;

  /* rho-rho prefactor */
  Y22Y22  *= 4.0 * M_PI;

  printf("%lf\n",Y22Y22 );
  

  //  return Y22Y22;

}

