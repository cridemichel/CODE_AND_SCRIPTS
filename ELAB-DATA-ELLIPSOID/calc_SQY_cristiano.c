#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "struct_conf_cristiano.h"

#define MY_TRUE  (1==1)
#define MY_FALSE (1==0)


// number of triples for each q-shell of the mesh
int ntriple[] =
#include "ntripl.incl"

// number of q-shells in the mesh
int Nshell = sizeof(ntriple)/sizeof(int);

// q-mesh
int qmesh[][150][3]=
#include "kmesh.incl"

// Ylm
const double Y00_pref = 0.28209479; // =  sqrt( 1.0/ (4.0*M_PI));
const double Y10_pref = 0.48860251; // =  sqrt( 3.0/ (4.0*M_PI));
const double Y20_pref = 0.31539157; // = sqrt( 5.0/ (4.0*M_PI))/2.0;
double *Ylm;

int main(int argc, char**argv){
  
  
#ifndef  BlockSize  
#define  BlockSize 1
#endif
  
  int iBlock;
  conf_t conf; // for the i/o of a cristiano's file
  int allocate = MY_TRUE;
  
#ifndef Npart
#define Npart 256
#endif
  int ipart;
  
  int  TotalConfs,iConf;
  char nomedir[1000];
  char nomeconf[1000];
  
  double unit_length, unit_time, fac_r2, fac_time;
  
  // for the FT
  int iq;
  double *qave;
  double *Sql0l0y;  
  double qmin;
  double qx,qy,qz;  
  
  // for the fourier transform
  double q_dot_r;
  double cos_qr[Npart],sin_qr[Npart];
  double rhoqYlm_real, rhoqYlm_imag;
  
  // for the euler angles
  double cosTheta;
  
  
  // for the spherical harmonics
  int ll,mm; 
  
  // to average over axis
  double*  axis;
  double   a_dot_x;
   
  
  int i,j;
 
   long int rseed;
#ifdef RSEED
  rseed=RSEED;
#else
  rseed=time(NULL);
#endif
  srand48(rseed); /* initialize random num generator */
  

  
  // read integer parameters l,m for the spherical harmonic Y_{l,m}
  if(argc==3) ll=atoi(argv[1]), mm=atoi(argv[2]);
  else fprintf(stderr,"\nI need ll, mm!\n\n"), exit(1);


  scanf("%d",&TotalConfs); // number of configurations
  scanf("%s",nomedir); // directory where configurations are            
  
  /* initialize averages */
  qave=(double *)malloc(Nshell*sizeof(double));
  Sql0l0y=(double *)malloc(Nshell*sizeof(double));  
  for(iq=0;iq<Nshell;iq++) qave[iq]=Sql0l0y[iq]=0.;
  
  // choose axis for Theta
  axis = (double*) malloc(3*sizeof(double));

  // allocate the Ylm i'll need
  Ylm = (double*) malloc(Npart*sizeof(double));
  
    
  /* I will average each of the BlockSize points over
   * TotalConfs/BlockSize values */		      
  for(iConf=0; iConf<(TotalConfs/BlockSize); iConf++){
    
    for(iBlock=0;iBlock<BlockSize;iBlock++) scanf("%s",nomeconf); 
    
    /* reads the iConf-th configuration */
    strcpy(conf.nomefile,nomedir); 
    strcat(conf.nomefile,nomeconf);
    read_cristiano_ellipses(&conf,allocate); allocate=MY_FALSE;
    if(conf.parnum!=Npart) 
      fprintf(stderr,"\nwrong Npart!\n\n"), exit(1); 
    
    qmin = 2.0 * M_PI / conf.box; // minimum q
    
#ifdef SCHILLING
    unit_length = 2.0*conf.aA;
#else
    unit_length = 2.0*pow(conf.aA*conf.bA*conf.cA,1./3.);
#endif
    fac_r2 = 1./(unit_length*unit_length);
    unit_time = sqrt(conf.mA/conf.kbT)*unit_length;
    fac_time  = 1.0/unit_time;
    
    
    
    for(iq=0;iq<Nshell;iq++){
      
      double mag,norm_axis;
      
      for(j=0;j<ntriple[iq];j++){
	
	qx = qmesh[iq][j][0]*qmin; 
	qy = qmesh[iq][j][1]*qmin; 
	qz = qmesh[iq][j][2]*qmin; 
	qave[iq] += sqrt(qx*qx+qy*qy+qz*qz);
	
	norm_axis = sqrt(qx*qx+qy*qy+qz*qz);
	axis[0] = qx/norm_axis;
	axis[1] = qy/norm_axis;
	axis[2] = qz/norm_axis;

	for(ipart=0;ipart<Npart;ipart++){

	  q_dot_r = ( qx*conf.rx[ipart] + 
		      qy*conf.ry[ipart] + 
		      qz*conf.rz[ipart]);
	  sin_qr[ipart]  =  sin(q_dot_r);
	  cos_qr[ipart]  =  cos(q_dot_r);	  

	  a_dot_x = conf.R[ipart][0]*axis[0] +
	    conf.R[ipart][1]*axis[1] + 
	    conf.R[ipart][2]*axis[2];
	  cosTheta = a_dot_x;
	  // cosTheta = (2.*drand48()-1.);
	
	  if ((ll==0)&&(mm==0)) 
	    Ylm[ipart] = Y00_pref;
	  else if ((ll==1)&&(mm==0))
	    Ylm[ipart] = Y10_pref * cosTheta;
	  else if ((ll==2)&&(mm==0))
	    Ylm[ipart] = Y20_pref * (3.0*cosTheta*cosTheta-1.0);
	  else{
	    fprintf(stderr,"Not yet ready for l=%d, m=%d !!!\n",ll,mm);
	    exit(1);
	  }


	}
	
	
	// calculate rho(q,l,m)
	rhoqYlm_real = 0.0; rhoqYlm_imag = 0.0;
	
	for(ipart=0;ipart<Npart;ipart++){
	  
	  rhoqYlm_real += sin_qr[ipart] * Ylm[ipart] ;
	  rhoqYlm_imag += cos_qr[ipart] * Ylm[ipart] ;
	    
	  }
	  
	  mag  = rhoqYlm_real*rhoqYlm_real 
	       + rhoqYlm_imag*rhoqYlm_imag;
	  Sql0l0y[iq] += mag/(double)Npart;
	  
	
	
      }//ntriple
      
    }//Nshell
    
    
  }//iConf

 // normalize qave, Sql0l0y

  for(iq=0;iq<Nshell;iq++){
  	qave[iq] = qave[iq] * unit_length / (double)ntriple[iq];
  	Sql0l0y[iq]  = Sql0l0y[iq] * 4.0 * M_PI /(double)ntriple[iq];    
  }
  for(iq=0;iq<Nshell;iq++) {
  	qave[iq] /= (double)(TotalConfs/BlockSize);
  	Sql0l0y[iq]  /= (double)(TotalConfs/BlockSize);
  }
  for(iq=0;iq<Nshell;iq++)
  	printf("%lf %lf\n",qave[iq],Sql0l0y[iq]);

  return 0;

}


static int check_open(FILE* file_p, char* nome_file){

  if(file_p==NULL) {
    fprintf(stderr,"Problem opening %s\n",nome_file);
    exit(1);
  }

  return 0;

}
