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
  double hist2[Nbins]; // histogram
  double hist4[Nbins];
  double hist6[Nbins];
  double cos2theta, cos4theta;
  double norm[Nbins]; // bin normalization
  double rx,ry,rz,r,dr;

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
  for(ibin=0;ibin<Nbins;ibin++) hist2[ibin]=hist4[ibin]=hist6[ibin]=norm[ibin]=0.0;
	      
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
#ifdef Z_AXIS
	cosTheta = conf.R[i][6] * conf.R[j][6] +
	  	conf.R[i][7] * conf.R[j][7] +
	  	conf.R[i][8] * conf.R[j][8];

#else
      	cosTheta = 
	  	conf.R[i][0] * conf.R[j][0] +
	  	conf.R[i][1] * conf.R[j][1] +
	  	conf.R[i][2] * conf.R[j][2];
#endif	  	
	  	rx = conf.rx[i] - conf.rx[j];
	  	ry = conf.ry[i] - conf.ry[j];
	  	rz = conf.rz[i] - conf.rz[j];

		rx=remainder(rx,box);	
		ry=remainder(ry,box);
		rz=remainder(rz,box);
		r=sqrt(rx*rx+ry*ry+rz*rz); ibin=r/dr;
		if((ibin>0)&&(ibin<Nbins)) {
			cos2theta = cosTheta*cosTheta;
			cos4theta = cos2theta*cos2theta;
			hist2[ibin] += cos2theta;
			hist4[ibin] += cos4theta;
			hist6[ibin] += cos2theta*cos4theta;
			norm[ibin]++ ;
		}
		
		
      }/* Npart x Npart */

	
    allocate = MY_FALSE; 

  }/* iConf */
    
  for(ibin=0;ibin<Nbins;ibin++) 
  	if(norm[ibin]>0.){
		hist2[ibin] /= norm[ibin];
		hist4[ibin] /= norm[ibin];
		hist6[ibin] /= norm[ibin];
		hist6[ibin]= (231.0*hist6[ibin] - 315.0*hist4[ibin] + 105.0*hist2[ibin] - 5.0)/16.0;
		hist4[ibin]= (35.0*hist4[ibin] - 30.0*hist2[ibin] + 3.0)/8.0;
  		hist2[ibin]=(3.*hist2[ibin]-1.)/2.;
  	}
  
  
  for(ibin=0;ibin<Nbins;ibin++) { 
  	r=(ibin+0.5)*dr; if(norm[ibin] > 0.0)
    printf("%lf %lf %lf %lf\n",r/unit_length,
    	hist2[ibin],hist4[ibin],hist6[ibin]);
  }
	
}

