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
double wI, wO;
double D=1.0;
int Nt = 10000, Nx=1000, NxL, stpout=1000;
double dx=1.0, dt = 0.5, L, Dx, nbuf=0.0, nbuf1=0.0;
double n[maxNx][2]; // distribuzione al tempo t e t+dt
void print_usage(void)
{
  printf("square_well_diff_eq <uI> <uO> <D> <Dx> <L> <Nx> <Nt> <dt> [ <stpout> ]\n");
  exit(-1);	
}

int main(int argc, char **argv)
{ 
  int i,j;
  double cost;     // costante di integrazione
  FILE *uscita, *equilib;

  /* apri il file di uscita */
  uscita=fopen("conc_profile.dat","w+");
  equilib=fopen("n_vs_t.dat", "w+");
  if (argc < 8)
    print_usage();

  betauI = atof(argv[1]);
  betauO = atof(argv[2]);
  D = atof(argv[3]);
  Dx = atof(argv[4]);
  L = atof(argv[5]);
  Nx = atoi(argv[6]);
  if (Nx > maxNx)
    {
      printf("Troppi punti (maxNx=%d)!\n", maxNx);
      exit(-1);
    }
  //dx = atof(argv[7]);
  Nt = atoi(argv[7]);
  dt = atof(argv[8]);
  if (argc == 10)
    stpout = atoi(argv[9]);

  dx = L / Nx;
  NxL = Dx / dx;
  wI=exp(-betauI);
  wO=exp(-betauO);
  /* la condizione iniziale e` uno scalino centrato in
   * Nx/2 e di larghezza 2 width */
  for(i=0; i<Nx; i++) n[i][0]=0.;  
  for(i=Nx/2-width; i<=Nx/2+width; i++) n[i][0]=1.;

  /* in 0 le condizioni al bordo sono di tipo "radiation" */
  //for(j=0; j<2; j++) n[0][j] = n[Nx-1][j] = 0.; 
 /* per la stabilita` dell'algoritmo, la costante di 
   * integrazione deve essere << 1*/
  cost=D*dt/dx/dx;
  
  printf("cost=%G dx=%G wI=%G wI=%G L=%f NxL=%d\n", cost, dx, wI, wO, L, NxL);
  
  for(j=0; j<=Nt; j++){  // loop temporale 

    /* n[i][0] contiene la distribuzione al tempo t, 
     * n[i][1] conterra` la distribuzione al tempo t + dt 
     * il loop implementa la discretizzazione di
     * n(x,t+dt) = n(x,t) + D grad^2 n(x,t) dt             */
    
    /* radiation boundary condition */
    n[0][0] = n[1][0]*(1.0-dx*wI/D);
    /* reflection boundary condition */
    n[Nx-1][0] = n[Nx-2][0];

    /* radiation boundary condition (in->out)*/
    //ntmp = n[NxL][0];
    nbuf = n[NxL][0];
    nbuf1 = n[NxL+1][0];
    n[NxL][0] = n[NxL-1][0]*(1.0-dx*wO/D);
   
    //n[0][1] = n[1][1]*(dx*wI+D);  
    for(i=1; i < NxL; i++)                
      n[i][1] = n[i][0] 
	+ cost*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]);
    n[NxL-1][1] += nbuf1*D/dx;
    n[NxL][0] = 0.0;
    for(i=NxL+1; i<(Nx-1); i++)                
      n[i][1] = n[i][0] 
	+ cost*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]);
    n[NxL+1][1] += nbuf*wO;
    if((j%stpout==0) || (j==0))
      { // salva ogni 10000 passi 
  	  {
  	    for(i=0 ; i<Nx; i++) // formato 3D per gnuplot 
  	      fprintf(uscita, "%f %f\n", (((double)i)+0.5)*dx, n[i][1]);  
#ifdef GNUPLOT
  	    fprintf(uscita, "\n"); // linea vuota per gnuplot 
#else
  	    fprintf(uscita, "&\n"); // linea vuota per gnuplot 
#endif
  	    fprintf(equilib, "%f %f\n", dt*j, n[NxL-5][1]);
  	  }
      }
    n[Nx-2][1]+=1.0;
    /* copia la soluzione al tempo t+dt in n[x][0] */
    for(i=1; i<(Nx-1); i++) n[i][0]=n[i][1];   

  } 

  fprintf(stderr,"Dati in conc_profile.dat\n");
  fclose(uscita);
  fclose(equilib);

  return 0;
}
