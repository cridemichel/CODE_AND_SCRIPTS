#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Sqr(X) ((X)*(X))

double *rx, *ry, *rz, *R[3][3];
double ux, uy, uz, *vx, *vy, *vz, *wx, *wy, *wz;
double mass;
double temp;
double Ix, Iy, Iz;
double k=4.0;
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
/* ------------------------------ */
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;

  /* first row vector */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[0][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[0][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
#if 0
  if (typesArr[0].nspots==3 && type==0)
    {
      for (k=0; k < 3 ; k++)
	u[k] = R[1][k];
      vectProdVec(R[0], u, up);
      /* rotate randomly second axis */
      angle=4.0*acos(0.0)*ranf_vb();
      xx = cos(angle);
      yy = sin(angle);
      for (k=0; k < 3 ; k++)
	R[1][k] = u[k]*xx + up[k]*yy;
      //printf("calc_norm(R[1])=%.15G\n", calc_norm(R[1]));
    }
#endif
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
#if 0
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    Rt[k1][k2]=R[k2][k1];
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    R[k1][k2]=Rt[k1][k2];
#endif

  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
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
double calc_energy(int nummol)
{
  int i, k1;
  double K,wt[3];
  K = 0;
  for (i=0; i < nummol; i++)
    {
	 /* calcola tensore d'inerzia e le matrici delle due quadriche */
	  K += mass*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    K += Sqr(wt[k1])*Ix;
    }
  K *= 0.5;
  return K;
}

void angvel(double *wx, double *wy, double* wz, int Nm)
{
  int i;
  double inert;                 /* momentum of inertia of the molecule */
  double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, Mtot;
  
  Mtot = mass; /* total mass of molecule */

  inert = Ix; /* momentum of inertia */
 
  mean = 3.0*temp / inert;

  for (i = 0; i < Nm; i++)
    {
      xisq = 1.0;
      
      while (xisq >= 1.0)
	{
	  xi1  = ranf() * 2.0 - 1.0;
	  xi2  = ranf() * 2.0 - 1.0;
	  xisq = xi1 * xi1 + xi2 * xi2;
	}
      
      xi = sqrt (fabs(1.0 - xisq));
      ox = 2.0 * xi1 * xi;
      oy = 2.0 * xi2 * xi;
      oz = 1.0 - 2.0 * xisq;
      
      /* Renormalize */
      osq   = ox * ox + oy * oy + oz * oz;
      norm  = sqrt(fabs(osq));
      ox    = ox / norm;
      oy    = oy / norm;
      oz    = oz / norm;
      
      /* Choose the magnitude of the angular velocity
         NOTE: consider that it is an exponential distribution 
	 (i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
      
      osq   = - mean * log(ranf());
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      *wx = ox;
      *wy = oy;
      *wz = oz;
    }
}

void main(int argc, char **argv)
{
  FILE *outf, *inpf;
  int parnum, dummyint, k1, k2, i;
  double K, sf, pi, k, kh, Lx, Ly, Lz, dummydbl, Rt[3][3];
  double rTemp, sumx, sumy, sumz;

  inpf = fopen(argv[1], "r");

  fscanf(inpf, "%d %d %lf %lf %lf %lf\n", &parnum, &dummyint, &Lx, &Ly, &Lz, &dummydbl);

  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  vx = malloc(sizeof(double)*parnum);
  vy = malloc(sizeof(double)*parnum);
  vz = malloc(sizeof(double)*parnum);
  wx = malloc(sizeof(double)*parnum);
  wy = malloc(sizeof(double)*parnum);
  wz = malloc(sizeof(double)*parnum);

  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R[k1][k2] = malloc(sizeof(double)*parnum);
  pi = acos(0.0)*2.0;
  mass = k;
  Ix = Iy = Iz = mass*(Sqr(k)+1.0)*Sqr(0.5)/5.0;
  rTemp = sqrt(temp/mass);
  sumx = sumy = sumz = 0.0;
  for (i=0; i < parnum; i++)
    {
      fscanf(inpf, "%lf %lf %lf %lf %lf %lf\n", &(rx[i]), &(ry[i]), &(rz[i]), &ux, &uy, &uz);
      versor_to_R(ux, uy, uz, Rt);
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	 R[k1][k2][i] = Rt[k1][k2]; 
      vx[i] = rTemp * gauss(); 
      vy[i] = rTemp * gauss();
      vz[i] = rTemp * gauss();
      sumx += vx[i]*mass;
      sumy += vy[i]*mass;
      sumz += vz[i]*mass;
    }
  sumx /= parnum*mass;
  sumy /= parnum*mass;
  sumz /= parnum*mass;
  for (i = 0; i < parnum; i++)
    {
      vx[i] -= sumx;
      vy[i] -= sumy;
      vz[i] -= sumz;
    }	
  K = calc_energy(parnum);
  sf = sqrt( ( (4.0*((double)parnum)-3.0) * temp ) / (2.0*K) );
  for (i = 0; i < parnum; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      wx[i] *= sf;
      wy[i] *= sf;
      wz[i] *= sf;
    } 

  angvel(wx, wy, wz, parnum);

  outf = fopen(argv[2], "w+");
  kh=k*0.5;
  fprintf(outf,"totStep: 50000\n");
  fprintf(outf,"equilibrat: 0\n");
  fprintf(outf,"curStep: 4400\n");
  fprintf(outf,"M: 5\n");
  fprintf(outf,"rcut: 16.766\n");
  fprintf(outf,"P: 1\n");
  fprintf(outf,"ninters: 0\n");
  fprintf(outf,"T: 1\n");
  fprintf(outf,"tol: 0\n");
  fprintf(outf,"time: 0\n");
  fprintf(outf,"Dt: 0.05\n");
  fprintf(outf,"parnum: %d\n", parnum);
  fprintf(outf,"ntypes: 1\n");
  fprintf(outf,"@@@\n");
  fprintf(outf,"%d\n",parnum); 
  fprintf(outf,"%f 0.5 0.5\n", kh); 
  fprintf(outf,"1 1 1\n"); 
  fprintf(outf,"%f %f %f %f 0 0\n", Ix, Iy, Iz, mass); 
  fprintf(outf,"0 0\n");
  fprintf(outf, "@@@\n");
  for (i=0; i < parnum; i++)
    {
      fprintf(outf, "%f %f %f %f %f %f %f %f %f %f %f %f 0\n", rx[i], ry[i], rz[i], R[0][0][i],R[0][1][i],R[0][2][i], 
	      R[1][0][i],R[1][1][i],R[1][2][i], R[2][0][i],R[2][1][i],R[2][2][i]);
    }
  for (i=0; i < parnum; i++)
    {
      fprintf(outf, "%f %f %f %f %f %f\n", vx[i], vy[i], vz[i], wx[i], wy[i], wz[i]);
    }
  fprintf(outf, "%f %f %f\n", Lx, Ly, Lz);
  fclose(outf);
  fclose(inpf);
}
