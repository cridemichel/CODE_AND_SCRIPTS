#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define MY_TRUE  (1==1)
#define MY_FALSE (1==0)

main(){

  /* The means square displacement is calculated from an initial configuration 
   * up to the following BlockSize-th configuration */
#ifndef  BlockSize  
#define  BlockSize 332
#endif
  
  int iBlock;
  conf_t conf[BlockSize]; /* i need BlockSize configurations */
  double C1[BlockSize];
  double C2[BlockSize];
  double C3[BlockSize];
  double C4[BlockSize];
  double time[BlockSize];  
  
 
  
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
  for(i=0;i<BlockSize;i++) C4[i]=C3[i]=C2[i]=C1[i]=time[i]=0.0;
  
  
  scanf("%d",&TotalConfs); // number of configurations
  scanf("%s",nomedir); // directory where configurations are         
  allocate = MY_TRUE;  // first Blocksize confs must be allocated    
  
  /* I will average each of the BlockSize points over
   * TotalConfs/BlockSize values */		      
  for(iConf=0; iConf<(TotalConfs/BlockSize); iConf++){
    
    /* reads the first configuration */
    strcpy(conf[0].nomefile,nomedir); 
    scanf("%s",nomeconf); 
    strcat(conf[0].nomefile,nomeconf);
    read_cristiano_ellipses(&conf[0],allocate);
    if(conf[0].parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}
    
    if(allocate){ // if it is the first time, calculate parameters
      
      /* assume that the size of the box does not change */
      box   = conf[0].box; 
      
      /* calculate the effective diameter */
#ifdef SCHILLING 
      unit_length = 2.0*conf[0].aA;
#else
      unit_length = 2.0*pow(conf[0].aA*conf[0].bA*conf[0].cA,1./3.);
#endif
      
      /* calculate rescaling factor to get things adimensional */
      fac_r2 = 1./(unit_length*unit_length);
      unit_time = sqrt(conf[0].mA/conf[0].kbT)*unit_length; 
      fac_time = 1.0/(unit_time*unit_time);
      
    }//end calculate parameters
    
    for(iBlock=1; iBlock<BlockSize; iBlock++) {
      
      /* read the iBlock-th configuration */
      strcpy(conf[iBlock].nomefile,nomedir); 
      scanf("%s",nomeconf); 
      strcat(conf[iBlock].nomefile,nomeconf);
      read_cristiano_ellipses(&conf[iBlock],allocate);
      
      /* average time at iBlock*/
      time[iBlock] += conf[iBlock].time-conf[0].time;
      
      for(i=0;i<Npart;i++){
      	
      	cosTheta = 
	  	conf[iBlock].R[i][0] * conf[0].R[i][0] +
	  	conf[iBlock].R[i][1] * conf[0].R[i][1] +
	  	conf[iBlock].R[i][2] * conf[0].R[i][2];
	
		C1[iBlock] += cosTheta;
		C2[iBlock] += ( 3.*cosTheta*cosTheta - 1.) / 2.;
		C3[iBlock] += ( 5.*cosTheta*cosTheta*cosTheta - 3.*cosTheta) / 2.;
		C4[iBlock] += ( 35.*cosTheta*cosTheta*cosTheta*cosTheta 
				      - 30.*cosTheta*cosTheta + 3. ) / 8.;
				      
		}/* i - end particles */
      
    } /* iBlock - end reading one block */
    
    /* after reading the first BlockSize configuration, I don't need anymore 
     * to allocate the memory for the configurations I am reading */ 
    allocate = MY_FALSE; 
    
  }/* iConf - end of reading all the blocks*/
  
  /* normalize and print the result */
  for(iBlock=1;iBlock<BlockSize;iBlock++){
    
    time[iBlock] /= (double)(TotalConfs/BlockSize);
    C1[iBlock]   /= (double)Npart * (double)(TotalConfs/BlockSize);
    C2[iBlock]  /= (double)Npart * (double)(TotalConfs/BlockSize);
    C3[iBlock]  /= (double)Npart * (double)(TotalConfs/BlockSize);
    C4[iBlock]  /= (double)Npart * (double)(TotalConfs/BlockSize);
    
    /* prints time, correlation of Y1 , correlation of Y2 */
    printf("%lf %lf %lf %lf %lf\n",
	   time[iBlock]*fac_time, 
	   C1[iBlock], C2[iBlock],
	   C3[iBlock], C4[iBlock] );
  }
  
}

