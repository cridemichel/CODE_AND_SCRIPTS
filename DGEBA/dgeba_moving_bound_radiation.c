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
double wI, wO, concSE, concS, V0=1.0, delta=0.1;
double D=1.0;
int Nt = 10000, Nx=1000, NxL, stpout=1000;
double dx=1.0, dt = 0.5, L1, L2, Dx, nbuf=0.0, nbuf1=0.0, r;
double n[maxNx][2]; // distribuzione al tempo t e t+dt
void print_usage(void)
{
  printf("square_well_diff_eq <uI> <uO> <D> <delta> <Dx> <L1> <L2> <Nx> <Nt> <dt> [ <stpout> ]\n");
  exit(-1);	
}
double U(double r, double L1, double Dx, double delta)
{
  return V0*(0.5+0.5*tanh(r-(L1+Dx)/delta));
}
double dU(double r, double L1, double Dx, double delta)
{
  double sh;
  sh = 1.0/cosh((r-(L1+Dx))/delta);
  return V0*0.5*Sqr(sh)/delta;
}
double ddU(double r, double L1, double Dx, double delta)
{
  double a, sh;
  a = L1+Dx;
  sh = 1.0/cosh((r-a)/delta);
  return -V0*Sqr(sh)*tanh((r-a)/delta)/Sqr(delta);
}
int main(int argc, char **argv)
{ 
  int i,j, k, rfI;
  double pos, rf, cost0, cost1, cost2;     // costante di integrazione
  FILE *uscita, *equilib, *kD;

  /* apri il file di uscita */
  uscita=fopen("conc_profile.dat","w+");
  equilib=fopen("n_vs_t.dat", "w+");
  kD = fopen("kD_vs_t.dat", "w+");
  if (argc < 8)
    print_usage();

  betauI = atof(argv[1]);
  V0 = atof(argv[2]);
  D = atof(argv[3]);
  delta = atof(argv[4]);
  Dx = atof(argv[5]);
  L1 = atof(argv[6]);
  L2 = atof(argv[7]);
  Nx = atoi(argv[8]);

  if (Nx > maxNx)
    {
      printf("Nx=%d Troppi punti (maxNx=%d)!\n", Nx, maxNx);
      exit(-1);
    }
  //dx = atof(argv[7]);
  Nt = atoi(argv[9]);
  dt = atof(argv[10]);
  if (argc == 12)
    stpout = atoi(argv[11]);

  dx = (L2-L1) / Nx;
  NxL = Dx / dx;
  wI=exp(-betauI)*dx/dt;
  wO=exp(-betauO)*dx/dt;
  /* la condizione iniziale e` uno scalino centrato in
   * Nx/2 e di larghezza 2 width */
  for(i=0; i<Nx; i++) n[i][0]=0.;  
  for(i=Nx/2-width; i<=Nx/2+width; i++) n[i][0]=1.0;

  /* in 0 le condizioni al bordo sono di tipo "radiation" */
  //for(j=0; j<2; j++) n[0][j] = n[Nx-1][j] = 0.; 
 /* per la stabilita` dell'algoritmo, la costante di 
   * integrazione deve essere << 1*/
  cost0 = D*dt;
  cost1 = D*dt/dx/dx;
  cost2 = D*dt/dx/2.0;
  printf("cost1=%G cost2=%G dx=%G wI=%G wI=%G L1=%f L2=%f NxL=%d\n", cost1, cost2, dx, wI, wO, L1, L2, NxL);
  printf("delta=%f Dx=%f V0=%f\n", delta, Dx, V0);
  pos = Dx;
  rfI = pos/dx;
  rf = L1 + pos;
  printf("out flux calculated at r=%f\n", rf);
  for(j=0; j<=Nt; j++){  // loop temporale 

    /* n[i][0] contiene la distribuzione al tempo t, 
     * n[i][1] conterra` la distribuzione al tempo t + dt 
     * il loop implementa la discretizzazione di
     * n(x,t+dt) = n(x,t) + D grad^2 n(x,t) dt             */
    /* radiation boundary conditions */
    n[0][0] = n[1][0]*(1.0-dx*wI/D);
    /* reflection boundary condition */
    n[Nx-1][0] = n[Nx-2][0];

    for(i=1; i<(Nx-1); i++)                
      {
	r = ((double)i)*dx+L1;
	n[i][1] = n[i][0] + cost0*n[i][0]*(dU(r, L1, Dx, delta)*2.0/r+ddU(r, L1, Dx, delta)) +  
	  + cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]) + 
	  + cost2*(dU(r, L1, Dx, delta)+2.0/r)*(n[i+1][0]-n[i-1][0]);
      }
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
	    concS = 0.0;
	    for (k=NxL+1; k < Nx; k++)
	      {
		concS += n[k][0];
	      }
	    concS /= Nx - NxL;

#endif	

	     fprintf(kD, "%G %G %G\n", dt*j, 4.0*M_PI*Sqr(L1)*n[1][0]*wI, 
		    -D*4.0*M_PI*Sqr(rf)*((n[rfI+1][0]-n[rfI-1][0])/2.0/dx + dU(rf, L1, Dx, delta)*n[rfI][0]));
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
