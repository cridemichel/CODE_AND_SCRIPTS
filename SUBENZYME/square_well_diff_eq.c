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

#if 0
#define dx 1.	     /* passo spaziale */
#define Nx 1000      /* numero di punti spaziali (PARI) */
#define dt 0.5001       /* passo temporale */
#define Nt 10000    /* passi di integrazione da eseguire */
#define D  1.0	     /* coefficiente di diffusione */
#endif
#define width 100       /* larghezza della distribuzione iniziale */

double betauI=10.0,  betauO=10.0; 
double wI, wO, j0, jDx, jL;
double D=1.0;
int Nt = 10000, Nx=1000;
double dx=1.0, dt = 0.5, L, Dx, NxL;
void print_usage(void)
{
  printf("square_well_diff_eq <uI> <uO> <D> <Dx> <L> <Nx> <dx> <Nt> <dt>\n");
}

int main(int argc, char **argv)
{ 
  int i,j;
  double cost;     // costante di integrazione
  double n[Nx][2]; // distribuzione al tempo t e t+dt
  FILE *uscita;

  /* apri il file di uscita */
  uscita=fopen("conc_profile.dat","w");

  if (argc < 9)
    print_usage();

  betauI = atof(argv[1]);
  betauO = atof(argv[2]);
  D = atof(argv[3]);
  Dx = atof(argv[4]);
  L = atof(argv[5]);
  Nx = atoi(argv[6]);
  dx = atof(argv[7]);
  Nt = atoi(argv[8])
  dt = atof(argv[9]);

  NxL = L / dx;
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
  
  for(j=0; j<=Nt; j++){  // loop temporale 

    /* n[i][0] contiene la distribuzione al tempo t, 
     * n[i][1] conterra` la distribuzione al tempo t + dt 
     * il loop implementa la discretizzazione di
     * n(x,t+dt) = n(x,t) + D grad^2 n(x,t) dt             */
    
    /* radiation boundary condition */
    n[0][0] = n[1][0]*(dx*wI/D+1);
    /* reflection boundary condition */
    n[Nx-1][0] = n[Nx-2][0];

    /* radiation boundary condition (in->out)*/
    n[NxL][0] = n[NxL-1][0]/(dx*wO/D+1);
   
    //n[0][1] = n[1][1]*(dx*wI+D);  
    for(i=1; i<(NxL-1); i++)                
      n[i][1] = n[i][0] 
	+ cost*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]);
   
    for(i=NxL+1; i<(Nx-1); i++)                
      n[i][1] = n[i][0] 
	+ cost*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]);

    if((j%1000==0) || (j==0)){ // salva ogni 10000 passi 
      for(i=0 ; i<Nx; i++) // formato 3D per gnuplot 
	fprintf(uscita, "%f\n", n[i][1]);  
      fprintf(uscita, "\n"); // linea vuota per gnuplot 
    }

    /* copia la soluzione al tempo t+dt in n[x][0] */
    for(i=1; i<(Nx-1); i++) n[i][0]=n[i][1];   

  } 

  fprintf(stderr,"Dati in conc_profile.dat\n");
  fclose(uscita);

  return 0;
}
