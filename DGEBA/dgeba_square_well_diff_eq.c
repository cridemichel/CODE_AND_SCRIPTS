/***********************************************************************
*  fick.c                                                              
*  programma per la soluzione dell'equazione di diffusione             
*  in un intervallo illimitato                                          
*  compilare con cc -lm fick.c                                             
*  uscita per grafico 3D gnuplot nel file fick.dat                     
*  comando: splot 'fick.dat'          
* 
*  l'equazione di Fick
* 
* 		dn(x,t)/dt = D grad^2 n(x,t)
* 
*  viene discretizzata nel tempo con un metodo di Eulero
* 
* 		n(x,t+dt) = n(x,t) + D grad^2 n(x,t) dt
* 
*  con la discretizzazione spaziale del Laplaciano 
* 
* 		grad^2 n = [n(x+dx) -2 n(x) + n(x-dx)] / (dx*dx)    
* 
*               
***********************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#define Sqr(x) ((x)*(x))
#if 0
#define dx 1.	     /* passo spaziale */
#define Nx 1000      /* numero di punti spaziali (PARI) */
#define dt 0.5001       /* passo temporale */
#define Nt 10000    /* passi di integrazione da eseguire */
#define D  1.0	     /* coefficiente di diffusione */
#endif
#define width 100       /* larghezza della distribuzione iniziale */

#define maxNx 1000000
double betauI=10.0,  betauO=10.0; 
double wI, wO, concSE, concS;
double D1=1.0, D2=1.0;
int Nt = 10000, Nx=1000, NxL, stpout=1000;
double dx=1.0, dt = 0.5, L1, L2, Dx, nbuf=0.0, nbuf1=0.0;
double n[maxNx][2]; // distribuzione al tempo t e t+dt
void print_usage(void)
{
  printf("square_well_diff_eq <uI> <D1> <L1> <L2> <Nx> <Nt> <dt> [ <stpout> ]\n");
  exit(-1);	
}

int main(int argc, char **argv)
{ 
  int i,j, k;
  double cost1, cost2;     // costante di integrazione
#ifndef UNIDIM
  double cost1b, cost2b;
#endif
  FILE *uscita, *equilib, *kD;

  /* apri il file di uscita */
  uscita=fopen("conc_profile.dat","w+");
  equilib=fopen("n_vs_t.dat", "w+");
  kD = fopen("kD_vs_t.dat", "w+");
  if (argc < 8)
    print_usage();

  betauI = atof(argv[1]);
//  betauO = atof(argv[2]);
  D1 = atof(argv[2]);
//  D2 = atof(argv[4]);
//  Dx = atof(argv[5]);
  L1 = atof(argv[3]);
  L2 = atof(argv[4]);
  Nx = atoi(argv[5]);
  if (Nx > maxNx)
    {
      printf("Troppi punti (maxNx=%d)!\n", maxNx);
      exit(-1);
    }
  //dx = atof(argv[7]);
  Nt = atoi(argv[6]);
  dt = atof(argv[7]);
  if (argc == 8)
    stpout = atoi(argv[8]);

  dx = (L2-L1) / Nx;
  NxL = L1 / dx;
  wI=exp(-betauI)*dx/dt;
 // wO=exp(-betauO)*dx/dt;
  /* la condizione iniziale e` uno scalino centrato in
   * Nx/2 e di larghezza 2 width */
  for(i=0; i<Nx; i++) n[i][0]=0.;  
  for(i=Nx/2-width; i<=Nx/2+width; i++) n[i][0]=10.0;

  /* in 0 le condizioni al bordo sono di tipo "radiation" */
  //for(j=0; j<2; j++) n[0][j] = n[Nx-1][j] = 0.; 
 /* per la stabilita` dell'algoritmo, la costante di 
   * integrazione deve essere << 1*/
  cost1 = D1*dt/dx/dx;
  //cost2 = D2*dt/dx/dx;
#ifndef UNIDIM
  cost1b = dt*D1/dx;
  //cost2b = dt*D2/dx;
#endif
  printf("cost1=%G x=%G wI=%G L1=%f L2=%f NxL=%d\n", cost1, dx, wI, L1, L2, NxL);
  
  
  for(j=0; j<=Nt; j++){  // loop temporale 

    /* n[i][0] contiene la distribuzione al tempo t, 
     * n[i][1] conterra` la distribuzione al tempo t + dt 
     * il loop implementa la discretizzazione di
     * n(x,t+dt) = n(x,t) + D grad^2 n(x,t) dt             */
    
#if 0
#if UNIDIM
    /* radiation boundary condition */
    n[0][0] = n[1][0]*(1.0-dx*wI/D1);
#else
    n[0][0] = n[1][0]*(1.0-dx*wI/D1);
#endif
#endif
    /* reflection boundary condition */
    n[Nx-1][0] = n[Nx-2][0];

    /* radiation boundary condition (in->out)*/
    //ntmp = n[NxL][0];
    
#if 0
    nbuf1 = n[NxL+1][0];
    n[NxL][0] = n[NxL-1][0]*(1.0-dx*wO/D1);
    nbuf = n[NxL][0];
    //n[0][1] = n[1][1]*(dx*wI+D);  
#ifdef UNIDIM
    for(i=1; i < NxL; i++)                
      n[i][1] = n[i][0] 
	+ cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]);
#else
    for(i=1; i < NxL; i++)                
      n[i][1] = n[i][0] 
	+ cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0])
	+ cost1b*(n[i+1][0]-n[i-1][0])/(((double)i)*dx+L1);
#endif
#ifdef UNIDIM
    n[NxL-1][1] += dt*nbuf1*D2/dx/dx;
#else
    n[NxL-1][1] += dt*nbuf1*D2/dx/dx;
#endif
#endif
#ifdef INNER_REFLECTION
    /* adsorbing boundary condition */
    n[0][0] = 0.0;
#else
    n[0][0] = n[1][0]*(1.0-dx*wI/D1);
#endif
#ifdef UNIDIM
    for(i=1; i<(Nx-1); i++)                
      n[i][1] = n[i][0] 
	+ cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]);
#else
    for(i=1; i<(Nx-1); i++)                
      n[i][1] = n[i][0] 
	+ cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0])
	+ cost1b*(n[i+1][0]-n[i-1][0])/(((double)i)*dx+L1);
#endif
#if 0
#ifdef UNIDIM
    n[NxL+1][1] += dt*nbuf*wO;
#else
    n[NxL+1][1] += dt*nbuf*wO;
#endif
#endif
    if((j%stpout==0) || (j==0))
      { // salva ogni 10000 passi 
  	  {
  	    for(i=0 ; i<Nx; i++) // formato 3D per gnuplot 
  	      fprintf(uscita, "%f %f\n", L1+(((double)i)+0.5)*dx, n[i][1]);  
#ifdef GNUPLOT
  	    fprintf(uscita, "\n"); // linea vuota per gnuplot 
#else
  	    fprintf(uscita, "&\n"); // linea vuota per gnuplot 
#endif
  	    fprintf(equilib, "%f %f\n", dt*j, n[NxL-5][1]);
#if 0
	    concSE = 0.0;
	    for (k=1; k < NxL-1; k++)
	      {
		concSE += n[k][0];
	      }
	    concSE /= NxL;
#endif
	    concS = 0.0;
	    for (k=1; k < Nx; k++)
	      {
#ifdef UNIDIM
		concS += n[k][0];
#else
		concS += 4.0*M_PI*n[k][0]*Sqr(k*dx+L1);
#endif
	      }
#ifdef UNIDIM
	    concS /= Nx;
	    if (concS!=0)
	      fprintf(kD, "%G %G\n", dt*j, n[1][0]*wI/concS);
#else
	    concS /= 4.0*M_PI*(L2*Sqr(L2)-L1*Sqr(L1))/3.0;
	    if (concS!=0)
	      fprintf(kD, "%G %G\n", dt*j, 4.0*M_PI*Sqr(L1)*n[1][0]*wI/concS);
#endif
	    //printf("n[NxL]=%G n[NxL-1]=%G\n", nbuf, n[NxL-1][0]);
  	  }
      }
    /* termine di sorgente legato a particelle che rientrano nel dominio */
    //n[Nx-2][1]+=1.0;
    /* copia la soluzione al tempo t+dt in n[x][0] */
    for(i=1; i<(Nx-1); i++) n[i][0]=n[i][1];   
  } 

  fprintf(stderr,"Dati in conc_profile.dat\n");
  fclose(uscita);
  fclose(equilib);
  fclose(kD);
  uscita=fopen("last_profile.dat", "w+");
  for(i=0 ; i<Nx; i++) // formato 3D per gnuplot 
    fprintf(uscita, "%f %f\n", L1+(((double)i)+0.5)*dx, n[i][1]);  
  fclose(uscita);
  return 0;
}
