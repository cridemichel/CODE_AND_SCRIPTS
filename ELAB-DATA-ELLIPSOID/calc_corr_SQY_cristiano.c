#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

// rho(q,t)
double **cosqr_0;
double **cosqr_t;
double **sinqr_0;
double **sinqr_t;

// Ylm(t)
const double Y00_pref = 0.28209479; // =  sqrt( 1.0/ (4.0*M_PI));
const double Y10_pref = 0.48860251; // =  sqrt( 3.0/ (4.0*M_PI));
const double Y20_pref = 0.31539157; // = sqrt( 5.0/ (4.0*M_PI))/2.0;
double **Ylm_0;
double **Ylm_t;


  // to average over q-axis
  double** q_axis;
  int      Nq_axis;
  int      iq_axis;
  double   q_dot_r;
 
  // to average over z-axis
  double** z_axis;
  int      Nz_axis;
  int      iz_axis;
  double   a_dot_x;



/* prototypes */

double **doublematrix(int row, int col);
double  *doublevector(int row);

int mesh2matrix(double **target, int index, int row, int col);
int normalize(double *vector);

int calcY(conf_t *config, double **Ylm, int ll, int mm);
int calcSQ(conf_t *config, double **cosq, double **sinqr, double qmag);

main(int argc, char **argv){

  /* The means square displacement is calculated from an initial configuration 
   * up to the following BlockSize-th configuration */
#ifndef  BlockSize  
#define  BlockSize 200
#endif
  
  int iBlock;
  conf_t conf[BlockSize]; /* i need BlockSize configurations */
  double corrSQYlm[BlockSize];
  double time[BlockSize];  
  
  
#ifndef Npart
#define Npart 256
#endif
int ipart;
  
  int  TotalConfs,iConf;
  char nomedir[1000];
  char nome[1000];
  char nomeconf[1000];
  int  allocate;
  double box,invbox;
  
  int i,j,k;
  
  double unit_length, unit_time, fac_r2, fac_time;
  
  // chosen value of q
  double qmag;

  
  // Euler angles
  double cosTheta;

  // chosen values for the spherical harmonics
  int ll,mm; 
  
  /* read magnitude of the q-vector for the fourier-transform 
   * and the integer parameters l,m for the spherical harmonic Y_{l,m} */
  if(argc==4) qmag=atof(argv[1]), ll=atoi(argv[2]), mm=atoi(argv[3]);  
  else { fprintf(stderr,"\nNeed only q,l,m\n\n");
  		 fprintf(stderr,"q must be in the units of the SQy graphs!\n\n"); 
  		 exit(1);}	
  
  
  // choose z-axis for Theta
  Nz_axis = ntriple[0];
  z_axis = doublematrix(Nz_axis,3);
  mesh2matrix(z_axis,0,Nz_axis,3);
  for(iz_axis=0;iz_axis<Nz_axis;iz_axis++) normalize(z_axis[iz_axis]);
  

  // allocate all the Ylm i need
  Ylm_0 = doublematrix(Npart,Nz_axis);
  Ylm_t = doublematrix(Npart,Nz_axis);
  
  // choose q-axis 
  Nq_axis = ntriple[10];
  q_axis = doublematrix(Nq_axis,3);
  mesh2matrix(q_axis,10,Nq_axis,3);
  for(iq_axis=0;iq_axis<Nq_axis;iq_axis++) normalize(q_axis[iq_axis]);

  // allocate all the rhoq i need
  cosqr_0 = doublematrix(Npart,Nq_axis);
  cosqr_t = doublematrix(Npart,Nq_axis);
  sinqr_0 = doublematrix(Npart,Nq_axis);
  sinqr_t = doublematrix(Npart,Nq_axis);

   
  /* initialize averages to zero*/
  for(i=0;i<BlockSize;i++) corrSQYlm[i]=time[i]=0.;
  
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
      box   = conf[0].box; //assume that the size of the box does not change 
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
    } /* iBlock - end reading one block */
    
    /* after reading the first block, I don't need anymore to allocate 
     * the memory for the configurations I am reading */ 
    allocate = MY_FALSE; 
    
    /* calculate Ylm's for the 0-th conf  */
    calcY(&conf[0],Ylm_0,ll,mm);
    calcSQ(&conf[0],cosqr_0,sinqr_0,qmag/unit_length);  
    
   
    for(iBlock=1; iBlock<BlockSize; iBlock++) { 
      
      calcY(&conf[iBlock],Ylm_t,ll,mm);
      calcSQ(&conf[iBlock],cosqr_t,sinqr_t,qmag/unit_length);
      
      for(ipart=0;ipart<Npart;ipart++) 
      for(iz_axis=0;iz_axis<Nz_axis;iz_axis++)
      for(iq_axis=0;iq_axis<Nq_axis;iq_axis++)
	corrSQYlm[iBlock] += 
	  Ylm_0[ipart][iz_axis]*Ylm_t[ipart][iz_axis]*
	  ( cosqr_0[ipart][iq_axis]*cosqr_t[ipart][iq_axis] + 
	    sinqr_0[ipart][iq_axis]*sinqr_t[ipart][iq_axis] );
 
    } /* iBlock - end averaging one block */
    
  }/* iConf - end of reading all the blocks*/
  
  /* normalize and print the result */
  for(iBlock=1;iBlock<BlockSize;iBlock++){
    
    corrSQYlm[iBlock] /= (double)Npart*(double)Nz_axis*(double)Nq_axis;
    
    time[iBlock]      /= (double)(TotalConfs/BlockSize);
    corrSQYlm[iBlock] /= (double)(TotalConfs/BlockSize);
    
    /* prints time, correlation of SQYlm */
    printf("%lf %lf\n",
	   time[iBlock]*fac_time, corrSQYlm[iBlock]);
  }

  return 0;
}



double *doublevector(int row){
	double *vector;
	
	vector = (double*) malloc( row*sizeof(double) );

    return vector;
}


double **doublematrix(int row, int col){
	double **matrix; int irow;
	
	matrix = (double**) malloc( row*sizeof(double*) );
	for(irow=0;irow<row;irow++)
    	matrix[irow] = (double*) malloc( col*sizeof(double) );
    	
    return matrix;
}

int mesh2matrix(double **target, int index, int row, int col){
	int irow,icol;
	
	for(irow=0;irow<row;irow++)
		for(icol=0;icol<col;icol++)
			target[irow][icol] = qmesh[index][irow][icol];
		
	return 0;
}

int normalize(double *vector){
	double comp, norm=0.0;
    int i;
    
    for(i=0;i<3;i++) comp = vector[i], norm += comp*comp;    
    norm=sqrt(norm);
    for(i=0;i<3;i++) vector[i] /= norm;
	
	return 0;	
}

int calcY(conf_t *config, double **Ylm, int ll, int mm){
	int ipart,iz_axis;
	double cosTheta,a_dot_x;
	
	for(ipart=0;ipart<Npart;ipart++) 
		for(iz_axis=0;iz_axis<Nz_axis;iz_axis++){
      		a_dot_x = config->R[ipart][0]*z_axis[iz_axis][0] + 
      				config->R[ipart][1]*z_axis[iz_axis][1] + 
      				config->R[ipart][2]*z_axis[iz_axis][2];
      		cosTheta = a_dot_x;
      	// cosTheta = (2.*drand48()-1.);
      	if ((ll==0)&&(mm==0)) 
			Ylm[ipart][iz_axis] = Y00_pref;
      	else if ((ll==1)&&(mm==0))
			Ylm[ipart][iz_axis] = Y10_pref * cosTheta;
      	else if ((ll==2)&&(mm==0))
			Ylm[ipart][iz_axis] = Y20_pref * (3.0*cosTheta*cosTheta-1.0);
      	else{
			fprintf(stderr,"Not yet ready for l=%d, m=%d !!!\n",ll,mm);
			exit(1);
      	}
    }

	return 0;
}

int calcSQ(conf_t *config, double **cos_qr, double **sin_qr, double qmag){ 
	
	double q_dot_r;
	int ipart;
	
	for(ipart=0;ipart<Npart;ipart++) 
		for(iq_axis=0;iq_axis<Nq_axis;iq_axis++){
			q_dot_r = qmag * (
				q_axis[iq_axis][0]*config->rx[ipart] + 
				q_axis[iq_axis][1]*config->ry[ipart] + 
				q_axis[iq_axis][2]*config->rz[ipart]);	
			cos_qr[ipart][iq_axis] = cos(q_dot_r);
			sin_qr[ipart][iq_axis] = sin(q_dot_r);
		}
		
	return 0;

}
