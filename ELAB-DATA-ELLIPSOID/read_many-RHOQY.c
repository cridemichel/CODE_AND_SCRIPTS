#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define mytrue  (1==1)
#define myfalse (1==0)

int    N_rhoq;
char   dir_rhoq_in[1000];
char   name_rhoq_in[1000];

double** rho_x;
double** rho_y;
double*  time;
int      N_q;
int      N_part;
double   S_q0;

double dot_prod(int n, double *xi,double *yi,double *xj,double *yj);
int read_rhoq(char* dir, char* name, int i);
double mag_sq(int n, double *xi,double *yi);

int main(int argc, char**argv){
  int i,j;
  
  fscanf(stdin,"%d",&N_rhoq);

  rho_x = (double**) malloc(N_rhoq*sizeof(double*));
  rho_y = (double**) malloc(N_rhoq*sizeof(double*));
  time  = (double*)  malloc(N_rhoq*sizeof(double));

  // read the N_rhoq configurations
  for(i=0;i<N_rhoq;i++){

    if(EOF==fscanf(stdin,"%s %s",dir_rhoq_in,name_rhoq_in)){
      fprintf(stderr,"Error reading name_rhoq_in\n\n");
      exit(1);
    }

    read_rhoq(dir_rhoq_in, name_rhoq_in, i);
    
  }


  S_q0=0;


  for(i=0;i<N_rhoq;i++)    
      S_q0 += 
	mag_sq(N_q, rho_x[i],rho_y[i])/((double)N_q*(double)N_rhoq);

  S_q0 /= (double) N_part;

  printf("%lf\n",S_q0);
  
  return 0;

}



int read_rhoq(char* dir, char* name, int i){
  FILE*  v_in;
  char name_v[2000];
  int j;


  sprintf(name_v,"%s/%s",dir,name);
  v_in=fopen(name_v,"r");
  if(v_in==NULL){
    fprintf(stderr,"Error opening %s\n\n",name_v);
    exit(1);
  }



  if(fscanf(v_in,"%d%d%lf",&N_q,&N_part,&time[i])==EOF){
    fprintf(stderr,"Error reading %s\n\n",name_v);
    exit(1);
  }
  rho_x[i] = (double *) malloc(N_q*sizeof(double));
  rho_y[i] = (double *) malloc(N_q*sizeof(double));
  
  for(j=0;j<N_q;j++)
    if(fscanf(v_in,"%lf %lf",&rho_x[i][j],&rho_y[i][j])==EOF){
      fprintf(stderr,"Error reading %s\n\n",name_v);
      exit(1);
    }


  fclose(v_in);

  return 0;
}



double dot_prod(int n, double *xi,double *yi,double *xj,double *yj){
  double sum = 0.0;
  int k;

  for(k=0;k<n;k++)
    sum += xi[k]*xj[k]+yi[k]*yj[k];

  return sum;

}

double mag_sq(int n, double *xi,double *yi){
  double sum = 0.0;
  int k;

  for(k=0;k<n;k++)
    sum += xi[k]*xi[k]+yi[k]*yi[k];

  return sum;

}
