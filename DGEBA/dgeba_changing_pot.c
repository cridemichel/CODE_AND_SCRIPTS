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
int Nt = 10000, Nx=1000, NxL, stpout=1000, stpout2=50000;
double dx=1.0, dt = 0.5, L1, L2, Dx, nbuf=0.0, nbuf1=0.0, r;
double n[maxNx][2]; // distribuzione al tempo t e t+dt
void print_usage(void)
{
  printf("square_well_diff_eq <uI> <D0> <n0> <Df> <L1> <L2> <Nx> <Nt> <dt> [ <stpout> ]\n");
  exit(-1);	
}
double U(double r, double L1, double Dx, double delta)
{
  return V0*(0.5-0.5*tanh(r-(L1+Dx)/delta));
}
double dU(double r, double L1, double Dx, double delta)
{
  double sh;
  sh = 1.0/cosh((r-(L1+Dx))/delta);
  return -V0*0.5*Sqr(sh)/delta;
}
double ddU(double r, double L1, double Dx, double delta)
{
  double a, sh;
  a = L1+Dx;
  sh = 1.0/cosh((r-a)/delta);
  return V0*Sqr(sh)*tanh((r-a)/delta)/Sqr(delta);
}
int main(int argc, char **argv)
{ 
  int i,j, k, i2R;
  double pos, cost0, cost1, cost2;     // costante di integrazione
  double M0, Sd, Df, M[2], D0, gamma, dndr2R, Rt, n0;	
  FILE *uscita, *equilib, *kD;

  /* apri il file di uscita */
  uscita=fopen("conc_profile.dat","w+");
  equilib=fopen("n_vs_t.dat", "w+");
  kD = fopen("kD_vs_t.dat", "w+");
  if (argc < 6)
    print_usage();

  betauI = atof(argv[1]);
  //V0 = atof(argv[2]);
  D0 = atof(argv[2]);
  n0 = atof(argv[3]); 
  delta = atof(argv[4]);
  Dx = atof(argv[5]);
  Df = atof(argv[6]);
  L1 = atof(argv[7]);
  L2 = atof(argv[8]);
  Nx = atoi(argv[9]);

  if (Nx > maxNx)
    {
      printf("Nx=%d Troppi punti (maxNx=%d)!\n", Nx, maxNx);
      exit(-1);
    }
  //dx = atof(argv[7]);
  Nt = atoi(argv[10]);
  dt = atof(argv[11]);
  if (argc == 12)
    stpout = atoi(argv[11]);
  if (argc == 13)
    stpout2 = atoi(argv[12]);


  dx = (L2-L1) / Nx;
  wI=exp(-betauI);
  //wO=exp(-betauO)*dx/dt;
  /* la condizione iniziale e` uno scalino centrato in
   * Nx/2 e di larghezza 2 width */
  for(i=0; i<Nx; i++) n[i][0]=0.;  
  //n0=0.01;
  for(i=0; i<Nx; i++) n[i][0]=n0;

  /* in 0 le condizioni al bordo sono di tipo "radiation" */
  //for(j=0; j<2; j++) n[0][j] = n[Nx-1][j] = 0.; 
 /* per la stabilita` dell'algoritmo, la costante di 
   * integrazione deve essere << 1*/
  printf("dx=%G wI=%G L1=%f L2=%f NxL=%d\n", dx, wI, L1, L2, NxL);
  printf("D0=%f Df=%f\n", D0, Df);
  //Df = 1.9; /* dimensione frattale dei cluster */
  //printf("out flux calculated at r=%f\n", rf);
  M0=1.0;
  M[0] = M0;
  gamma = 1.0/3.0;
  for(j=0; j<=Nt; j++){  // loop temporale 

    /* n[i][0] contiene la distribuzione al tempo t, 
     * n[i][1] conterra` la distribuzione al tempo t + dt 
     * il loop implementa la discretizzazione di
     * n(x,t+dt) = n(x,t) + D grad^2 n(x,t) dt             */
    /* radiation boundary conditions */
    /* reflection boundary condition */
    n[Nx-1][0] = n0*M0/M[0];

    D = 2.0*D0/pow(M[0],gamma);
    cost0 = D*dt;
    cost1 = D*dt/dx/dx;
    cost2 = D*dt/dx/2.0;

    Rt = L1*pow(M[0]/M0,1.0/Df);
    i2R = (2.0*Rt - L1)/dx;
    if (i2R > Nx-1)
      {
	printf("Rt=%f M[0]=%f i2R=%d\n", Rt, M[0], i2R);
	printf("[j=%d] cluster too big, increase L2!\n", j);
	exit(-1);
      }
#ifdef ABSORB
    n[i2R-1][0] = 0;
#else
    //n[i2R-1][0] = n[i2R][0]*(1.0-dx*wI/D);
#endif
    dndr2R = (n[i2R+2][0]-n[i2R][0])/dx/2.0;
    Sd=4.0*M_PI*Sqr(2.0*Rt);
    Dx = 2.0*delta;
    for(i=i2R; i<(Nx-1); i++)                
      {
	r = ((double)i)*dx+L1;
	n[i][1] = n[i][0] + cost0*n[i][0]*(dU(r, Rt, Dx, delta)*2.0/r+ddU(r, Rt, Dx, delta)) +  
	  + cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]) + 
	  + cost2*(dU(r, Rt, Dx, delta)+2.0/r)*(n[i+1][0]-n[i-1][0]) - cost0*n[i][0]*(Sd*dndr2R);
#if 0
	n[i][1] = n[i][0] - cost0*n[i][0]*(Sd*dndr2R) +  
	  + cost1*(n[i+1][0]+n[i-1][0]-2.0*n[i][0]) + 
	  + cost2*(2.0/r)*(n[i+1][0]-n[i-1][0]);
#endif
      }
    M[1] = M[0] + dt*M[0]*D*Sd*(dndr2R + dU(2.0*Rt, L1, Dx, delta)*n[i2R][0]);
    if((j%stpout2==0) || (j==0))
      {
	printf("j=%d Rt=%f\n", j, Rt);
      }
     if((j%stpout==0) || (j==0))
      { // salva ogni 10000 passi 
	//printf("j=%d Rt=%f\n", j, Rt);
#if 0
	for(i=i2R ; i<Nx; i++) // formato 3D per gnuplot 
	  fprintf(uscita, "%f %f\n", (L1+(((double)i)+0.5)*dx)/(2.0*Rt), M[1]*n[i][1]/M0/n0);  
#ifdef GNUPLOT
	fprintf(uscita, "\n"); // linea vuota per gnuplot 
#else
	fprintf(uscita, "&\n"); // linea vuota per gnuplot 
#endif
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

	fprintf(kD, "%G %G %G\n", dt*j, 4.0*M_PI*Sqr(2.0*Rt)*(dndr2R+dU(2.0*Rt, L1, Dx, delta)*n[i2R][0]), Rt);
	//printf("n[NxL]=%G n[NxL-1]=%G\n", nbuf, n[NxL-1][0]);
      }
    /* termine di sorgente legato a particelle che rientrano nel dominio */
    //n[Nx-2][1]+=1.0;
    /* copia la soluzione al tempo t+dt in n[x][0] */
    for(i=1; i<(Nx-1); i++) n[i][0]=n[i][1];   
    M[0] = M[1];
  } 

  fprintf(stderr,"Dati in conc_profile.dat\n");
  fclose(uscita);
  fclose(equilib);
  fclose(kD);
  uscita=fopen("last_profile.dat", "w+");
  for(i=i2R ; i<Nx; i++) // formato 3D per gnuplot 
    fprintf(uscita, "%f %f\n", (L1+(((double)i)+0.5)*dx)/(2.0*Rt), M[1]*n[i][1]/M0/n0);  
  fclose(uscita);
  return 0;
}
