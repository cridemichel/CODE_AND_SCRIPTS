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

  double histY22[Nbins]; // histogram
  double Y2i,Y2j,Y22;
  
  double box, unit_length, cosTheta;
  
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
  
  /* initialize histograms */
  for(ibin=0;ibin<Nbins;ibin++) hist[ibin]=norm[ibin]=0.0;
	      
  for(iConf=0; iConf<(TotalConfs); iConf++){
  	
  	/* reads the first configuration */
    strcpy(conf.nomefile,nomedir); 
    scanf("%s",nomeconf); strcat(conf.nomefile,nomeconf);
    read_cristiano_ellipses(&conf,allocate);
    if(conf.parnum!=Npart){fprintf(stderr,"cazzo!");exit(1);}

	box = conf.box; dr = (box/2.)/(double)Nbins;
#ifdef SCHILLING 
      unit_length = 2.0*conf.aA;
#else
      unit_length = 2.0*pow(conf.aA*conf.bA*conf.cA,1./3.);
#endif
	

    for(i=1;i<Npart;i++) for(j=0;j<i;j++){

      	cosTheta = 
	  	conf.R[i][0] * conf.R[j][0] +
	  	conf.R[i][1] * conf.R[j][1] +
	  	conf.R[i][2] * conf.R[j][2];
	  	
	  	rx = conf.rx[i] - conf.rx[j];
	  	ry = conf.ry[i] - conf.ry[j];
	  	rz = conf.rz[i] - conf.rz[j];

		rx=remainder(rx,box);	
		ry=remainder(ry,box);
		rz=remainder(rz,box);

		Y22 = 0;
		for (k = 0; k < 3; ++k) {
			Y2i=conf.R[i][k];
			Y2j=conf.R[j][k];
			Y2i=(3.*Y2i*Y2i-1.)*0.5;
			Y2j=(3.*Y2j*Y2j-1.)*0.5;
			Y22 += Y2i*Y2j;
		}
		Y22 /=3.;
	  
		r=sqrt(rx*rx+ry*ry+rz*rz); ibin=r/dr;
		if((ibin>0)&&(ibin<Nbins)) {
			hist[ibin] += cosTheta*cosTheta;
			histY22[ibin] += Y22; 
			norm[ibin]++ ;
		}
		
		
      }/* Npart x Npart */

	
    allocate = MY_FALSE; 

  }/* iConf */
    
  for(ibin=0;ibin<Nbins;ibin++) 
  	if(norm[ibin]>0.){
  		hist[ibin]=(3.*hist[ibin]/norm[ibin]-1.)/2.;
  		histY22[ibin]/=norm[ibin];
  	}
  	else
  		histY22[ibin]=hist[ibin]=0.; 	
  
  
  for(ibin=0;ibin<Nbins;ibin++) { 
  	r=(ibin+0.5)*dr; 
    printf("%lf %lf %lf %lf\n",r/unit_length,
    	hist[ibin],histY22[ibin],norm[ibin]);
  }
	
}

