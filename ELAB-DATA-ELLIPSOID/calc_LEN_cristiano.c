#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "struct_conf_cristiano.h"

#define MY_TRUE  (1==1)
#define MY_FALSE (1==0)

#define Nbins 100

main(){
	
  conf_t conf;    // configuration
  int ibin; 
  double hist[Nbins]; // histogram
  double norm[Nbins]; // bin normalization
  double rx,ry,rz,r,dr;
  
  double box, unit_length, cosTheta;
  
  double corr_len;
  
#ifndef Npart
#define Npart 256
#endif

  int  TotalConfs,iConf;
  char nomedir[1000];
  char nome[1000];
  char nomeconf[1000];
  int  allocate;
   
  int i,j,k;

  
  scanf("%d",&TotalConfs); // number of configurations
  scanf("%s",nomedir); // directory where configurations are         
  allocate = MY_TRUE;  // first conf must be allocated  
	      
  for(iConf=0; iConf<(TotalConfs); iConf++){
  	
  	/* reads the first configuration */
    strcpy(conf.nomefile,nomedir); 
    scanf("%s",nomeconf); strcat(conf.nomefile,nomeconf);
    read_cristiano_ellipses(&conf,allocate);
    if(conf.parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}

	box = conf.box; dr = box/(double)Nbins;
#ifdef SCHILLING 
      unit_length = 2.0*conf.aA;
#else
      unit_length = 2.0*pow(conf.aA*conf.bA*conf.cA,1./3.);
#endif


  	/* initialize histograms */
  	for(ibin=0;ibin<Nbins;ibin++) hist[ibin]=norm[ibin]=0.0;


    for(i=1;i<Npart;i++) for(j=0;j<i;j++){

      	cosTheta = 
	  	conf.R[i][0] * conf.R[j][0] +
	  	conf.R[i][1] * conf.R[j][1] +
	  	conf.R[i][2] * conf.R[j][2];
	  	
	  	rx = conf.R[i][0] - conf.R[j][0];
	  	ry = conf.R[i][1] - conf.R[j][1];
	  	rz = conf.R[i][2] - conf.R[j][2];

		rx=remainder(rx,box);	
		ry=remainder(ry,box);
		rz=remainder(rz,box);
	  
		r=sqrt(rx*rx+ry*ry+rz*rz); ibin=r/dr;
		if((ibin>0)&&(ibin<Nbins)) {
			hist[ibin] += cosTheta*cosTheta; 
			norm[ibin]++ ;
		}
		
		
      }/* Npart x Npart */

  	for(ibin=0;ibin<Nbins;ibin++) if(norm[ibin]>0.) 
  		hist[ibin]=(3*hist[ibin]/norm[ibin]-1.)/2.;
  	hist[0] = 1.0;

	corr_len = 0.;
	for(ibin=0;ibin<Nbins;ibin++) corr_len += hist[ibin];
	corr_len *= dr;
	printf("%lf %lf %s\n",corr_len,unit_length,conf.nomefile);
	
    allocate = MY_FALSE; 

  }/* iConf */
    
  
}

