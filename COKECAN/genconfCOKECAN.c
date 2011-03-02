#include <stdlib.h>
#include <stdio.h>
#include <math.h>
double T=1.0, m=1.0, L[3];
int N;
double *rx, *ry, *rz, *vx, *vy, *vz;
double ranf(void);

void FCC(FILE *f)
{
  /*   DESCRIPTION:
       Sets up the alpha fcc lattice for n linear molecules.   
       The simulation box is a unit cube centred at the origin.
       N should be an integer of the form ( 4 * ( Nc ** 3 ) ),
       Where Nc is the number of FCC unit cells in each direction.  
       See figure 5.10 for a diagram of the lattice and a           
       definition of the four orientational sublattices.            
       PRINCIPAL VARIABLES:                                         
       int           Nm                   Number of molecules             
       COORD_TYPE    d                    diatomic molecule length 
       COORD_TYPE    rxCm, ryCm, rzCm     Molecular Center of mass 
                                          positions             
       COORD_TYPE    ex, ey, ez           half of vector joining atom a and b 
                                          in a molecule 
       COORD_TYPE    rRoot3               1.0 / sqrt ( 3.0 ) */
  int cc, deli, Nc, cc1, cc2, type;
  double Cell[3], Cell2[3];
  int i, ix, iy, iz, ii;
  Nc = ceil(  pow( ((double)N), 1.0/3.0 )  );
  fprintf(stderr,"Nc: %d\n", Nc);
  /* Calculate the side of the unit cell */
  cc=cc1=cc2=0;
  L[0]=L[1]=L[2]=Nc*5;

  for (ii=0; ii < 3; ii++)
    {
      Cell[ii] = L[ii] / ((double) Nc); /* unit cell length */
      Cell2[ii] = 0.5 * Cell[ii];              /* half unit cell length */
    }
  /* Construct the lattice from the unit cell */
  /* assumiamo N2 > N1 */
  ii = 0;
  for(iz = 0; iz < Nc; iz++) /* loops over unit cells (that are simply cubes) */ 
    {
      for(iy = 0; iy < Nc; iy++)
	{
	  for(ix = 0; ix < Nc; ix++)
	    {
	      if ( cc >= N ) 
		break;	
	      /* Center of Mass of the actual molecule (m + iref) */
	      rx[cc] = Cell[0] * ((double) ix);
	      ry[cc] = Cell[1] * ((double) iy);
	      rz[cc] = Cell[2] * ((double) iz);
	      rx[cc] = rx[cc] - 0.5 * L[0] + ranf()*1E-7+0.1; 
	      ry[cc] = ry[cc] - 0.5 * L[1] + ranf()*1E-7+0.1;
	      rz[cc] = rz[cc] - 0.5 * L[2] + ranf()*1E-7+0.1;
	      fprintf(f, "%f %f %f 1 0 0 0 1 0 0 0 1 0\n", rx[cc], ry[cc], rz[cc]);
	    }
	}
    }
  
  /* Shift centre of box to the origin */
  
  return;
}
double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return rand() / ( (double) RAND_MAX );
}

double gauss(void)
{
  
  /* 
     Random variate from the standard normal distribution.
     
     The distribution is gaussian with zero mean and unit variance.
     REFERENCE:                                                    
                                                                
     Knuth D, The art of computer programming, (2nd edition        
     Addison-Wesley), 1978                                      
                                                                
     ROUTINE REFERENCED:                                           
                                                                
     COORD_TYPE ranf()                                  
     Returns a uniform random variate on the range zero to one  
  */

  double  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0;i < 12; i++)
    {
      sum = sum + ranf();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}


void vels(FILE *f)
{
  int i;
  double sumx, sumy, sumz, rTemp;
  rTemp = sqrt(T / m); 
  sumx=sumy=sumz=0;
  for (i = 0; i < N; i++)
    {
      vx[i] = rTemp * gauss(); 
      vy[i] = rTemp * gauss();
      vz[i] = rTemp * gauss();
      sumx += vx[i]*m;
      sumy += vy[i]*m;
      sumz += vz[i]*m;
    }
  sumx /= N*m;
  sumy /= N*m;
  sumz /= N*m;
  for (i = 0; i < N; i++)
    {
      vx[i] -= sumx;
      vy[i] -= sumy;
      vz[i] -= sumz;
      fprintf(f, "%.15G %.15G %.15G 0 0 0\n", vx[i], vy[i], vz[i]);
    }
}
int main(int argc, char** argv)
{
  int i, M, i1, i2, i3, ds;
  N=atoi(argv[1]);
  fprintf(stderr,"N=%d\n", N);
  ds=sizeof(double);
  rx = malloc(ds*N);
  ry = malloc(ds*N);
  rz = malloc(ds*N);
  vx = malloc(ds*N);
  vy = malloc(ds*N);
  vz = malloc(ds*N);
  printf("parnum:%d\n", N);
  printf("ntypes: 1\n");        
  printf("ninters: 3\n");
  printf("nintersIJ: 0\n");
  //printf("parnumA: %d\n", N1+N2);
  printf("@@@\n");	
  printf("%d\n", N); 
  printf("2 1 1\n");
  printf("16 2 2\n");
  printf("1 1 1 1 2 0\n");
  printf("2 0\n");
  printf("1.541665 0 0 1.21667\n");
  printf("-1.541665 0 0 1.21667\n");
  printf("0 0 0 0 1 0 0 100000\n");
  printf("0 0 0 1 1 0 0 100000\n");
  printf("0 1 0 1 1 0 0 100000\n");
  printf("@@@\n");	
  FCC(stdout);
  vels(stdout);
  //fprintf(stdout, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  fprintf(stdout, "%.15G\n", L[0]);
  return 0;
}
