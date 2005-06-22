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
#define  BlockSize 99
#endif

  int iBlock;
  conf_t conf[BlockSize]; /* i need BlockSize configurations */
  double   r2[BlockSize];
  double  dr2[BlockSize];
  double time[BlockSize];
  
#ifndef Npart
#define Npart 256
#endif

  /* temporary vectors to follow the evolution in time of 
   * the XYZ coordinates */
  double rx[Npart],ry[Npart],rz[Npart];

  int  TotalConfs,iConf;
  char nomedir[1000];
  char nome[1000];
  char nomeconf[1000];
  int  allocate;
  double box,invbox;
  
  int i;
  double dx,dy,dz;

  double unit_length, unit_time, fac_r2, fac_time;

  /* initialize averages */
  for(i=0;i<BlockSize;i++) r2[i]=dr2[i]=time[i]=0.;

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
    unit_time = unit_length * sqrt(conf[0].mA/conf[0].kbT); 
    fac_time = 1.0/(unit_time*unit_time);

	/* at the beginning, particles are in zero */
    for(i=0;i<Npart;i++) rx[i]=ry[i]=rz[i]=0.;

    for(iBlock=1; iBlock<BlockSize; iBlock++) {
      
      /* read the iBlock-th configuration */
      strcpy(conf[iBlock].nomefile,nomedir); 
      scanf("%s",nomeconf); 
      strcat(conf[iBlock].nomefile,nomeconf);
      read_cristiano_ellipses(&conf[iBlock],allocate);
      
      /* average time at iBlock*/
      time[iBlock] += conf[iBlock].time-conf[0].time;

      for(i=0;i<Npart;i++){
	
		/* calculate the XYZ displacements between two 
		 * consecutive configurations */
		dx = conf[iBlock].rx[i] - conf[iBlock-1].rx[i];
		dy = conf[iBlock].ry[i] - conf[iBlock-1].ry[i];
		dz = conf[iBlock].rz[i] - conf[iBlock-1].rz[i];
	
		/* adjust XYZ displacements so that it lies inside the box */
		dx=remainder(dx,box);	
		dy=remainder(dy,box);
		dz=remainder(dz,box);
		/* This is reasonable IFF the particle has not moved 
		 * more than one box */
		
		/* calculate the average displacement between two 
		 * consecutive configurations */	
		dr2[iBlock] += dx*dx + dy*dy + dz*dz;
	
		/* adjust the displacement of the i-th particle */
		rx[i] += dx; ry[i] += dy; rz[i] += dz;
		
		/* calculate the average displacement between the initial 
		 * and the iBlock-th configuration*/
		r2[iBlock] += rx[i]*rx[i] +ry[i]*ry[i] +rz[i]*rz[i];
	
      }/* i - end particles */

    } /* iBlock - end reading one block */

	/* after reading the first block, I don't need anymore to allocate 
	 * the memory for the configurations I am reading */ 
    allocate = MY_FALSE; 

  }/* iConf - end of reading all the blocks*/
  
  /* normalize and print the result */
  for(iBlock=1;iBlock<BlockSize;iBlock++){

	time[iBlock] /= (double)(TotalConfs/BlockSize);
    r2[iBlock]   /= (double)Npart * (double)(TotalConfs/BlockSize);
    dr2[iBlock]  /= (double)Npart * (double)(TotalConfs/BlockSize);
    
    /* prints time, m.s.d.(t) , D(t) */
    printf("%lf %lf %lf %lf\n",
	   time[iBlock]*fac_time,
	   r2[iBlock]*fac_r2,
	   r2[iBlock]*fac_r2/(6.*time[iBlock]*fac_time),
	   dr2[iBlock]*fac_r2);
  }

}

