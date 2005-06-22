#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define MY_TRUE  (1==1)
#define MY_FALSE (1==0)

main(){
	
  conf_t conf;    // configuration 
  double S[3][3]; // order parameter matrix 
  double norm;
  
#ifndef Npart
#define Npart 256
#endif

  int  TotalConfs,iConf;
  char nomedir[1000];
  char nome[1000];
  char nomeconf[1000];
  int  allocate;
   
  int i,j,k;

  /* initialize averages */
  for(j=0;j<3;j++) for(k=0;k<3;k++) S[j][k]=0.0;
  
  scanf("%d",&TotalConfs); // number of configurations
  scanf("%s",nomedir); // directory where configurations are         
  allocate = MY_TRUE;  // first conf must be allocated  
	      
  for(iConf=0; iConf<(TotalConfs); iConf++){
  	
  	/* reads the first configuration */
    strcpy(conf.nomefile,nomedir); 
    scanf("%s",nomeconf); strcat(conf.nomefile,nomeconf);
    read_cristiano_ellipses(&conf,allocate);
    if(conf.parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}
	
      for(i=0;i<Npart;i++){
      	norm = 
	  	conf.R[i][0] * conf.R[i][0] +
	  	conf.R[i][1] * conf.R[i][1] +
	  	conf.R[i][2] * conf.R[i][2];
	  
		for(j=0;j<3;j++) for(k=0;k<3;k++) 
			S[j][k] += conf.R[i][j]*conf.R[i][k]/norm;
			
      }/* Npart */

    allocate = MY_FALSE; 

  }/* iConf */
  
  /* normalize the result */
  for(j=0;j<3;j++) for(k=0;k<3;k++) S[j][k] /= (double) TotalConfs;
  for(j=0;j<3;j++) for(k=0;k<3;k++) S[j][k] /= (double) Npart;
  
  /* Sij = m_i m_j - delta_ij / 3 */
  for(j=0;j<3;j++) S[j][j] -= 1.0/3.0;
  
  /* write output */
  for(j=0;j<3;j++) printf("%lf %lf %lf\n",S[j][0],S[j][1],S[j][2]);
  
}

