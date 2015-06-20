//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
#define SYMMETRY
//#define USEGSL
#define GAUSS
#define EULER_ROT
//#define MCGAMMA
#ifdef USEGSL
#include <gsl/gsl_qrng.h>
#endif
#if defined(MPI)
int MPIpid;
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
#endif 
//#define ELEC
//#define ALBERTA
//#define NO_INTERP
double **XI1, **XI2, **XI3, **XI4, **XI5, **XI6;
#ifdef ELEC
double kD, yukcut, yukcutkD, yukcutkDsq;
#endif
#ifdef PARALLEL
struct kDsortS {
int k1;
int k2;
double invkD;
} *kD_sorted;

/* Matrice di Eulero come viene definita nel Goldstein x' = Reul*x dove x' sono 
   le coordinate nel sistema di riferimento del corpo rigido e x in quello
   del laboratorio. */
double Reul[3][3];
void build_euler_matrix(double cosphi, double sinphi, double costheta, double sintheta,
			double cospsi, double sinpsi)
{ /* psi = gamma */
  Reul[0][0] = cospsi*cosphi-costheta*sinphi*sinpsi;
  Reul[0][1] = cospsi*sinphi+costheta*cosphi*sinpsi;
  Reul[0][2] = sinpsi*sintheta;
  Reul[1][0] = -sinpsi*cosphi-costheta*sinphi*cospsi;
  Reul[1][1] = -sinpsi*sinphi+costheta*cosphi*cospsi;
  Reul[1][2] = cospsi*sintheta;
  Reul[2][0] = sintheta*sinphi;
  Reul[2][1] = -sintheta*cosphi;
  Reul[2][2] = costheta;
}
double *cdna_arr, *beta_arr;
int numtemps, numconcs;
double **yukcutkD_arr, **kD_arr, **yuk_corr_fact_arr, **yukcutkDsq_arr, **uel_arr, **vexclel_arr;
double num_kD=0;
double maxyukcutkDsq, maxyukcutkD;
#endif
char dummy1[32], dummy2[32], atname[32], nbname[8];
int type, nat, atnum, nbnum, len;
long long int tot_trials, tt=0, ttini=0;
double L, rx, ry, rz, alpha, dfons_sinth_max, fons_sinth_max, ROMBTOL, Lx, Ly, Lz;
const double thetapts=100000;

int ntheta, nphi, ngamma;
double *xtheta, *xphi, *wtheta, *wphi, *xgamma, *xphi, *wgamma;
#define EPSGAULEG 3.0e-11 
/*EPS is the relative precision.*/
void gauleg(double x1, double x2, double x[], double w[], int n)
  /*Given the lower and upper limits of integration x1 and x2, and given n, this routine returns arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss- Legendre n-point quadrature formula.*/
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  /*
     High precision is a good idea for this rou- tine.
     The roots are symmetric in the interval, so we only have to find half of them.
     Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
     Loop over the desired roots. */
  for (i=1;i<=m;i++) 
    {
      z=cos(3.141592654*(i-0.25)/(n+0.5));
      /*Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.*/
      do {
	p1=1.0;
	p2=0.0;
	for (j=1;j<=n;j++) {
	  p3=p2;
	  p2=p1; p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	}
	/*
	   p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a 
	   standard relation involving also p2, the polynomial of one lower order. 
	 */	   
	pp=n*(z*p1-p2)/(z*z-1.0);
	z1=z;
	z=z1-p1/pp;
      } while (fabs(z-z1) > EPSGAULEG); 
      x[i]=xm-xl*z;
      x[n+1-i]=xm+xl*z; 
      w[i]=2.0*xl/((1.0-z*z)*pp*pp); 
      w[n+1-i]=w[i];
      /*Scale the root to the desired interval, and put in its symmetric counterpart. 
	Compute the weight and its symmetric counterpart.*/
    }
}
#if 0
int j;
float xr,xm,dx,s;
static float x[]={0.0,0.1488743389,0.4333953941,
0.6794095682,0.8650633666,0.9739065285}; static float w[]={0.0,0.2955242247,0.2692667193, 0.2190863625,0.1494513491,0.0666713443};
The abscissas and weights. First value of each array not used.
xm=0.5*(b+a); xr=0.5*(b-a);
s=0;
for (j=1;j<=5;j++) {
dx=xr*x[j];
Will be twice the average value of the function, since the ten weights (five numbers above each used twice) sum to 2.
s += w[j]*((*func)(xm+dx)+(*func)(xm-dx)); }
return s *= xr;
#endif
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
struct DNADallStr {
  double R[3][3];
  double rcm[3];
  double sax[3];
} DNADall[2];

double calcDistBox(void)
{
  double RR, R0, R1, cij[3][3], fabscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k, k1, k2, existsParallelPair = 0;
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  for (k=0; k < 3; k++)
    {
      rA[k] = DNADall[0].rcm[k];
      rB[k] = DNADall[1].rcm[k];
      EA[k] = DNADall[0].sax[k];
      EB[k] = DNADall[1].sax[k];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  AA[k1][k2] = DNADall[0].R[k1][k2];
	  BB[k1][k2] = DNADall[1].R[k1][k2];
	}
    	DD[k1] = rA[k1] - rB[k1];
    }
  /* axis C0+s*A0 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[0][k1] =  scalProd(AA[0], BB[k1]);
      fabscij[0][k1] = fabs(cij[0][k1]);
      if ( fabscij[0][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[0] = scalProd(AA[0],DD);
  RR = fabs(AD[0]);
  R1 = EB[0]*fabscij[0][0]+EB[1]*fabscij[0][1]+EB[2]*fabscij[0][2];
  R01 = EA[0] + R1;
  if ( RR > R01 )
    return 1.0; /* non si intersecano */
  /* axis C0+s*A1 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[1][k1] = scalProd(AA[1],BB[k1]);
      fabscij[1][k1] = fabs(cij[1][k1]);
      if ( fabscij[1][k1] == 1.0  )
	existsParallelPair = 1;
    }
  AD[1] = scalProd(AA[1],DD);
  RR = fabs(AD[1]);
  R1 = EB[0]*fabscij[1][0]+EB[1]*fabscij[1][1]+EB[2]*fabscij[1][2];
  R01 = EA[1] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*A2 */
  for (k1= 0; k1 < 3; k1++)
    {
      cij[2][k1] = scalProd(AA[2], BB[k1]);
      fabscij[2][k1] = fabs(cij[2][k1]);
      if ( fabscij[2][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[2] = scalProd(AA[2],DD);
  RR = fabs(AD[2]);
  R1 = EB[0]*fabscij[2][0]+EB[1]*fabscij[2][1]+EB[2]*fabscij[2][2];
  R01 = EA[2] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*B0 */
  RR = fabs(scalProd(BB[0],DD));
  R0 = EA[0]*fabscij[0][0]+EA[1]*fabscij[1][0]+EA[2]*fabscij[2][0];
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B1 */
  RR = fabs(scalProd(BB[1],DD));
  R0 = EA[0]*fabscij[0][1]+EA[1]*fabscij[1][1]+EA[2]*fabscij[2][1];
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;
  
  /* axis C0+s*B2 */
  RR = fabs(scalProd(BB[2],DD));
  R0 = EA[0]*fabscij[0][2]+EA[1]*fabscij[1][2]+EA[2]*fabscij[2][2];
  R01 = R0 + EB[2];
  if ( RR > R01 )
    return 1.0;

  /* At least one pair of box axes was parallel, therefore the separation is
   * effectively in 2D, i.e. checking the "edge" normals is sufficient for
   * the separation of the boxes. 
   */
  if ( existsParallelPair )
    return -1.0;

  /* axis C0+s*A0xB0 */
  RR = fabs(AD[2]*cij[1][0]-AD[1]*cij[2][0]);
  R0 = EA[1]*fabscij[2][0] + EA[2]*fabscij[1][0];
  R1 = EB[1]*fabscij[0][2] + EB[2]*fabscij[0][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB1 */
  RR = fabs(AD[2]*cij[1][1]-AD[1]*cij[2][1]);
  R0 = EA[1]*fabscij[2][1] + EA[2]*fabscij[1][1];
  R1 = EB[0]*fabscij[0][2] + EB[2]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB2 */
  RR = fabs(AD[2]*cij[1][2]-AD[1]*cij[2][2]);
  R0 = EA[1]*fabscij[2][2] + EA[2]*fabscij[1][2];
  R1 = EB[0]*fabscij[0][1] + EB[1]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB0 */
  RR = fabs(AD[0]*cij[2][0]-AD[2]*cij[0][0]);
  R0 = EA[0]*fabscij[2][0] + EA[2]*fabscij[0][0];
  R1 = EB[1]*fabscij[1][2] + EB[2]*fabscij[1][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB1 */
  RR = fabs(AD[0]*cij[2][1]-AD[2]*cij[0][1]);
  R0 = EA[0]*fabscij[2][1] + EA[2]*fabscij[0][1];
  R1 = EB[0]*fabscij[1][2] + EB[2]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB2 */
  RR = fabs(AD[0]*cij[2][2]-AD[2]*cij[0][2]);
  R0 = EA[0]*fabscij[2][2] + EA[2]*fabscij[0][2];
  R1 = EB[0]*fabscij[1][1] + EB[1]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB0 */
  RR = fabs(AD[1]*cij[0][0]-AD[0]*cij[1][0]);
  R0 = EA[0]*fabscij[1][0] + EA[1]*fabscij[0][0];
  R1 = EB[1]*fabscij[2][2] + EB[2]*fabscij[2][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB1 */
  RR = fabs(AD[1]*cij[0][1]-AD[0]*cij[1][1]);
  R0 = EA[0]*fabscij[1][1] + EA[1]*fabscij[0][1];
  R1 = EB[0]*fabscij[2][2] + EB[2]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB2 */
  RR = fabs(AD[1]*cij[0][2]-AD[0]*cij[1][2]);
  R0 = EA[0]*fabscij[1][2] + EA[1]*fabscij[0][2];
  R1 = EB[0]*fabscij[2][1] + EB[1]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  return -1.0;
}

void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
//#define ALBERTA
char fn[1024];
struct DNA {
  double x;
  double y;
  double z;
  double rad;
#ifdef ELEC
  int atype;
#endif
} *DNAchain;

struct DNA *DNADs[2];
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

#define MC_BENT_DBLCYL

/*
                pdb        radius (angstrom)
    sugar       Xe          3.5        
    phosphate   B           3.0         
    base        Se          4.0    
*/
void body2lab(double xp[3], double x[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
      x[k1] += rO[k1];
    }
}
void print_matrix(double M[3][3], int n)
{
  int k1, k2;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}

void versor_to_R(double ox, double oy, double oz, double gamma, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
#ifdef MC_BENT_DBLCYL
  double Rout[3][3];
  int k1, k2;
#endif
  /* first row vector */
  R[2][0] = ox;
  R[2][1] = oy;
  R[2][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[2][0] && u[1]==R[2][1] && u[2]==R[2][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[2][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[2][k];
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
  vectProdVec(R[1], R[2], u);
 
  for (k=0; k < 3 ; k++)
    R[0][k] = u[k];
#if 0
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    Rt[k1][k2]=R[k2][k1];
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    R[k1][k2]=Rt[k1][k2];
#endif
  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
#ifdef DEBUG
  printf("==============\n");
  print_matrix(R, 3);
  printf("==============\n");
#endif
}
double RMDNA[2][3][3];
#if defined(EULER_ROT)
void place_DNAD(double x, double y, double z, double cosphi12, double sinphi12, double costheta12, 
		double sintheta12, double cosgamma12, double singamma12,
		int which)
{
  FILE *fd;
  char fn[256];
  double xp[3], rO[3], xl[3];
  double R[3][3];
  int i, k1, k2;
  rO[0] = x;
  rO[1] = y;
  rO[2] = z;
  R[0][0] = cosgamma12*costheta12*cosphi12 - singamma12*sinphi12;
  R[0][1] = cosphi12*singamma12 + cosgamma12*costheta12*sinphi12;
  R[0][2] = -cosgamma12*sintheta12;
  R[1][0] = -costheta12*cosphi12*singamma12-cosgamma12*sinphi12;
  R[1][1] = cosgamma12*cosphi12-costheta12*singamma12*sinphi12;
  R[1][2] = singamma12*sintheta12;
  R[2][0] = cosphi12*sintheta12;
  R[2][1] = sintheta12*sinphi12;
  R[2][2] = costheta12; 
#ifdef DEBUG 
  /*printf("costhe=%f sinth=%f cosphi=%f sinphi=%f cosgmma=%f singamma=%f\n", 
	 costheta12, sintheta12, cosphi12, sinphi12, cosgamma12, singamma12);*/
  print_matrix(R,3);
  sprintf(fn, "DNAD%d.mgl", which);
  fd=fopen(fn, "w+");
  fprintf(fd, ".Vol: %f\n", L*L*L);
#endif
  /* ============ */
  for (k1=0; k1 < 3; k1++)
    {
      DNADall[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	DNADall[which].R[k1][k2] = R[k1][k2];
    }
  for (i=0; i < nat; i++)
    {
      xp[0] = DNAchain[i].x;
      xp[1] = DNAchain[i].y;
      xp[2] = DNAchain[i].z;

      body2lab(xp, xl, rO, R);
      DNADs[which][i].x = xl[0];
      DNADs[which][i].y = xl[1];
      DNADs[which][i].z = xl[2];
#ifdef ELEC
      DNADs[which][i].atype = DNAchain[i].atype;
#endif
#ifdef DEBUG
      fprintf(fd,"%f %f %f @ %f\n", xl[0], xl[1], xl[2], DNAchain[i].rad);
#endif
      DNADs[which][i].rad = DNAchain[i].rad;
    }
#ifdef DEBUG
  fclose(fd);
  exit(-1);
#endif
} 
#else
void place_DNAD(double x, double y, double z, double ux, double uy, double uz, double gamma, int which)
{
  double xp[3], rO[3], xl[3];
  double R[3][3];
  int k1, k2;
#ifdef DEBUG
  FILE *fd;
  char fn[128];
#endif
  FILE *f;
  int i; 
  rO[0] = x;
  rO[1] = y;
  rO[2] = z;
  /* build R here from the orientation (ux,uy,uz) */
  versor_to_R(ux, uy, uz, gamma, R);
#ifdef DEBUG 
  sprintf(fn, "DNAD%d.mgl", which);
  fd=fopen(fn, "w+");
  fprintf(fd, ".Vol: %f\n", L*L*L);
#endif
  /* ============ */
  for (k1=0; k1 < 3; k1++)
    {
      DNADall[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	DNADall[which].R[k1][k2] = R[k1][k2];
    }
  for (i=0; i < nat; i++)
    {
      xp[0] = DNAchain[i].x;
      xp[1] = DNAchain[i].y;
      xp[2] = DNAchain[i].z;

      body2lab(xp, xl, rO, R);
      DNADs[which][i].x = xl[0];
      DNADs[which][i].y = xl[1];
      DNADs[which][i].z = xl[2];
#ifdef ELEC
      DNADs[which][i].atype = DNAchain[i].atype;
#endif
#ifdef DEBUG
      fprintf(fd,"%f %f %f @ %f\n", xl[0], xl[1], xl[2], DNAchain[i].rad);
#endif
      DNADs[which][i].rad = DNAchain[i].rad;
    }
#ifdef DEBUG
  fclose(fd);
#endif
}
#endif
/* ============================ >>> ranf <<< =============================== */
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
void init_distbox(void)
{
  int i, k;
  double max_x, max_y, max_z, distx, disty, distz;
  /* maximum distance to z-axis */
  for (i=0; i < nat; i++)
    {
      distx = fabs(DNAchain[i].x) + DNAchain[i].rad;
      disty = fabs(DNAchain[i].y) + DNAchain[i].rad;
      distz = fabs(DNAchain[i].z) + DNAchain[i].rad;
#ifdef ELEC
      if (DNAchain[i].atype==1)/* se si tratta di un P */
	{
#ifdef PARALLEL
	  if (numtemps > 1 || numconcs > 1)
	    {
	      yukcutkD = maxyukcutkD;
	    }
#endif
	  if (yukcutkD*0.5 > DNAchain[i].rad)
	    {
	      distx = fabs(DNAchain[i].x) + yukcutkD*0.5;
	      disty = fabs(DNAchain[i].y) + yukcutkD*0.5;
	      distz = fabs(DNAchain[i].z) + yukcutkD*0.5;
	    }
	}
#endif
      if (i==0 || distx > max_x)
	max_x = distx;
       if (i==0 || disty > max_y)
	max_y = disty;
       if (i==0 || distz > max_z)
	max_z = distz;
    } 
  for (k=0; k < 2; k++)
    {
      DNADall[k].sax[0] = max_x*1.01;
      DNADall[k].sax[1] = max_y*1.01;
      DNADall[k].sax[2] = max_z*1.01;
    }
  //printf("maxx=%f %f\n",DNADall[0].sax[0],DNADall[0].sax[1]);
}
double max3(double a, double b, double c)
{
  double m;
  m = a;
  if (b > m)
    m = b;
  if (c > m)
    m = c;
  return m;
}
double max2(double a, double b)
{
  if (a > b)
    return a;
  else 
    return b;
}
double epsrtbl[21][2]={{0,87.740},{5,85.763},{10,83.832},{15,81.946},{20,80.103},{25,78.304},{30,76.546},{35,74.828},{40,73.151},{45,71.512},{50,69.910},{55,68.345},{60,66.815},{65,65.319},{70,63.857},{75,62.427},{80,61.027},{85,59.659},{90,58.319},{95,57.007},{100,55.720}};
double epsr(double T)
{
  int i;
  i = 0;
  T = T - 273.15;
  if (T <= 0.0 || T >= 100.0)
    {
      printf("Temperature must be between 0 and 100!\n");
      exit(-1);
    }
  for (i=0; i < 21; i++)
    {
      if (T > epsrtbl[i][0])
	{
	  /* linear interpolation */
	  return epsrtbl[i][1]+(T-epsrtbl[i][0])*(epsrtbl[i+1][1]-epsrtbl[i][1])/(epsrtbl[i+1][0]-epsrtbl[i][0]);
	}
    }
}
double integrandv1(double rcmx, double rcmy, double rcmz, double gamma12,
		    double phi12, double theta12)
{
  int i, j;
  double sigsq, distsq, sigijsq, u1z, u2x, u2y, u2z;
  double sintheta12, costheta12, sinphi12, cosphi12, cosgamma12, singamma12;

  costheta12 = cos(theta12);
  sintheta12 = sin(theta12);
  cosphi12 = cos(phi12);
  sinphi12 = sin(phi12);
  cosgamma12 = cos(gamma12);
  singamma12 = sin(gamma12);
  //versor_to_R(u1x, u1y, u1z, gamma1, DNADall[0].R);
  //versor_to_R(u2x, u2y, u2z, gamma2, DNADall[1].R);
#ifdef EULER_ROT
  place_DNAD(rcmx, rcmy, rcmz, cosphi12, sinphi12, costheta12, sintheta12, cosgamma12, singamma12, 1);
#else
  u2x = sintheta12*cosphi12;
  u2y = sintheta12*sinphi12;
  u2z = costheta12;  
  place_DNAD(rcmx, rcmy, rcmz, u2x, u2y, u2z, gamma12, 1);
#endif
  if (calcDistBox() < 0.0)
    {
      for (i=0; i < nat; i++)
	{
	  for (j=0; j < nat; j++)
	    {
	      distsq = Sqr(DNADs[0][i].x-DNADs[1][j].x) + Sqr(DNADs[0][i].y-DNADs[1][j].y) + Sqr(DNADs[0][i].z-DNADs[1][j].z);
	      sigijsq = Sqr(DNADs[0][i].rad + DNADs[1][j].rad);
	      if (distsq < sigijsq)
		{
		  return -1.0;
		}
	    }
	}
    }
  return 0.0;
}

double phi12sav, theta12sav, gamma12sav;
#ifdef QFGAUSS
double rcmxsav, nrcmxsav, rcmysav, nrcmysav, rcmzsav, nrcmzsav;
int nrcmx, nrcmy, nrcmz;
double *xrcmx, *wrcmx, *xrcmy, *wrcmy, *xrcmz, *wrcmz;
#endif
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
#define MAXBIT 30 
#define MAXDIM 6
void sobseq(int *n, double x[])
/*When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[1..n] the next values from n of these sequences. (n must not be changed between initializations.)*/
{
  int j,k,l;
  unsigned long i,im,ipp;
  static double fac;
  static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
  static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
  static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4}; 
  static unsigned long iv[MAXDIM*MAXBIT+1]={0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
  if (*n < 0) 
    { 
      /*Initialize, don’t return a vector. */
      for (k=1;k<=MAXDIM;k++) ix[k]=0;
      in=0;
      if (iv[1] != 1) return;
      fac=1.0/(1L << MAXBIT);
      for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) 
	iu[j] = &iv[k];/* To allow both 1D and 2D addressing.*/
      for (k=1;k<=MAXDIM;k++) 
	{
	  for (j=1;j<=mdeg[k];j++)
	    iu[j][k] <<= (MAXBIT-j); /*Stored values only require normalization.*/
	  for (j=mdeg[k]+1;j<=MAXBIT;j++) 
	    {
	      ipp=ip[k]; i=iu[j-mdeg[k]][k];
	      i ^= (i >> mdeg[k]);
	      for (l=mdeg[k]-1;l>=1;l--) 
		{
		  if (ipp & 1) i ^= iu[j-l][k];
		  ipp >>= 1; 
		}
	      iu[j][k]=i;
	    }
	}
    } 
  else 
    {
      im=in++;
      for (j=1;j<=MAXBIT;j++) {
	if (!(im & 1)) break;
	im >>= 1; }
      if (j > MAXBIT) {
	printf("MAXBIT too small in sobseq");
	exit(-1);
      } 
      im=(j-1)*MAXDIM;
      for (k=1;k<=IMIN(*n,MAXDIM);k++)
	{
	  ix[k] ^= iv[im+k]; 
	  x[k]=ix[k]*fac;
	}
      /*XOR the appropriate direction num- ber into each component of the vector and convert to a floating number.
       */
    }
}
int main(int argc, char**argv)
{
#ifdef QUASIMC
#ifdef USEGSL
  gsl_qrng *qsob;
#endif
#endif
#ifdef MPI
  MPI_Status status;
#endif
#ifdef ELEC
  double uel, beta;
  int interact;
#endif
#if defined(QUASIMC) || defined(QFGAUSS)
  double sv[10];
  int nsv;
#endif
  unsigned int *overlaparr;
  int ircmx, ircmy, ircmz, iphi, itheta, igamma;
  char fn[256];
  int aa, bb;
  double ccc, totfact;
  int cc;
  double gamma1, gamma2;
  FILE *fin, *fout, *f, *fread, *fxi1, *fxi2, *fxi3, *fxi4, *fxi5, *fxi6;
  int ncontrib, k, i, j, overlap, contrib, cont=0, nfrarg;
  long long int fileoutits, outits;
  char fnin[1024],fnout[256];
  double dummydbl, segno, u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, vexclel=0.0, factor, dth, th;
  /* syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */
  if (argc < 5)
    {
      printf("syntax:  build_mesh_for_FG <pdb file> <DNAD length> [outits] [nrcmx] [nrcmy] [nrcmz] [ngamma] [nphi] [ntheta]\n");
      exit(-1);
    }
  strcpy(fnin,argv[1]);
  fin=fopen(fnin,"r");
  len=atoi(argv[2]);
  alpha = atof(argv[3]);
  outits = atoll(argv[4]);
	
  if (argc == 5)
    nrcmx = 20;
  else
    nrcmx = atoi(argv[5]);

  if (argc == 6)
    nrcmy = 20;
  else
    nrcmy = atoi(argv[6]);

  if (argc == 7)
    nrcmz = 20;
  else
    nrcmz = atoi(argv[7]);

  if (argc == 8)
    ngamma = 10;
  else
    ngamma = atoi(argv[8]);

  if (argc == 9)
    nphi = 10;
  else
    nphi = atoi(argv[9]);

  if (argc == 10)
    ntheta = 10;
  else
    ntheta = atoi(argv[10]);

   /*
     rho_salt =2 csalt Nav 1000;
     rho_counter[cdna_]:=(2 cdna)/(660*Dalton);
     (* cdna in mg/ml e csalt Molare (=moli/litro), nu=numero di cariche per unità di lunghezza *)
     InvDebyeScrLen[T_,qdna_,cdna_,qsalt_,csalt_,\[Epsilon]rel_]:= Sqrt[qdna^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) ( 2cdna)/(660*Dalton)+qsalt^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) 2 csalt Nav 1000 ];
     InvDebyeScrLen[300, 2, 200, 2, 1, 20]^-1*10^9

   */
  /* qdna è la carica rilasciata da ogni gruppo fosfato in soluzione (tipicamente=1) */
  /* ELISA: ATOM    39   Xe   G A   14      -5.687  -8.995  37.824 */
  /* ALBERTA: HETATM    1  B            1     -1.067  10.243 -35.117 */
  /* len here is the number of dodecamers, where 70 is the number of atoms per dodecamers
     in our CG model */
  nat = 70*len;
  DNAchain = (struct DNA*) malloc(sizeof(struct DNA)*nat); 
  for (k=0; k < 2; k++)
    DNADs[k] = (struct DNA*) malloc(sizeof(struct DNA)*nat);
  //L = 1.05*3.0*40*len; /* 4 nm is approximately the length of a 12 bp DNAD */ 
  /* read the CG structure */
  cc=0;
  while (!feof(fin))
    {
#ifdef ALBERTA
      fscanf(fin, "%s %d %s %d %lf %lf %lf ", dummy1, &atnum, atname, &nbnum, &rx, &ry, &rz);
#else
      fscanf(fin, "%s %d %s %s %s %d %lf %lf %lf ", dummy1, &atnum, atname, nbname, dummy2, &nbnum, &rx, &ry, &rz);
#endif
      DNAchain[cc].x = rx;
      DNAchain[cc].y = ry;
      DNAchain[cc].z = rz;
      //printf("cc=%d (%f,%f,%f)\n", cc, rx, ry, rz);
#ifdef ALBERTA
      if (!strcmp(atname, "S"))
	{
	  DNAchain[cc].rad = 3.5;
#ifdef ELEC
	  DNAchain[cc].atype = 0;
#endif
	}
      else if (!strcmp(atname, "P"))
	{
	  DNAchain[cc].rad = 3.0;
#ifdef ELEC
	  DNAchain[cc].atype = 1;
#endif
	}
      else if (!strcmp(atname, "B"))
	{
	  DNAchain[cc].rad = 4.0;
#ifdef ELEC
	  DNAchain[cc].atype = 2;
#endif
	}
#else
      if (!strcmp(atname, "Xe"))
	{
	  DNAchain[cc].rad = 3.5;
#ifdef ELEC
	  DNAchain[cc].atype = 0;
#endif
	}
      else if (!strcmp(atname, "B"))
	{
	  DNAchain[cc].rad = 3.0;
#ifdef ELEC
	  DNAchain[cc].atype = 1;
#endif
	}
      else if (!strcmp(atname, "Se"))
	{
	  DNAchain[cc].rad = 4.0;
#ifdef ELEC
	  DNAchain[cc].atype = 2;
#endif
	}
#endif     
      else
	{
	  printf("Unrecognized atom name, exiting...\n");
	  exit(1);
	}
      cc++;
      if (cc >= nat)
	break;
    };
  rcmx=rcmy=rcmz=0.0;
  for (i=0; i < nat; i++)
    {
      rcmx += DNAchain[i].x;
      rcmy += DNAchain[i].y;
      rcmz += DNAchain[i].z;
    }
  rcmx /= (double) nat;
  rcmy /= (double) nat;
  rcmz /= (double) nat;
  for (i=0; i < nat; i++)
    {
      DNAchain[i].x -= rcmx;
      DNAchain[i].y -= rcmy;
      DNAchain[i].z -= rcmz;
    }
  fclose(fin);
  init_distbox();
  L=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*3.0;
  printf("nat=%d L=%f alpha=%f I am going to calculate v%d and I will do %lld trials\n", nat, L, alpha, type, tot_trials);
  printf("box semiaxes=%f %f %f\n", DNADall[0].sax[0], DNADall[0].sax[1], DNADall[0].sax[2]);
#ifdef MPI
  srand48(((int)time(NULL))+my_rank);
#else
  srand48((int)time(NULL));
#endif
  sprintf(fnout, "v%d.dat", type);
  factor=0.0;
#if 0
  dfons_sinth_max=estimate_maximum_dfons(alpha);
  fons_sinth_max=dfons_sinth_max/alpha;
  printf("Estimated maximum of dfons is %f\n", dfons_sinth_max);
#endif
#ifdef QFGAUSS
  printf("Quasi Full Gauss quadrature with nrcmx=%d nrcmy=%d nrcmz=%d points\n", nrcmx, nrcmy, nrcmz);
#endif
#ifdef GAUSS
  printf("Gauss quadrature with %d %d %d points\n", nphi, ntheta, ngamma);
#else
  printf("Romberg method with %.15G tolerance\n", ROMBTOL);
#endif
#ifdef QUASIMC
  printf("I will generate a Quasi Monte Carlo sequence\n");
#endif
  //exit(-1);
  /* avendo diviso l'integrazione in theta negli intervalli [0,pi/2] e [pi/2,pi]
     il fattore si deve ottenere integrando fra 0 e pi/2 */
#if 0
  dth=acos(0.0)/((double)thetapts);
  th=0.0;
  //f = fopen("dfons.dat", "w+");
  for (i=0; i < thetapts; i++)
    {
      factor += 0.5*dth*sin(th)*(dfons(th, alpha) + dfons(th+dth,alpha));
      th += dth;
      //fprintf(f,"%f %.15G\n", th, dfons(th, alpha));
    };
  //fclose(f);
  factor= fabs(factor);
  factor *= 4.0*acos(0.0);
#else
  factor = alpha/2.0;
#endif
  printf("factor=%.15G\n", factor);
  fout = fopen(fnout, "w+");
  fclose(fout);
  Lx=Ly=Lz=L;
  printf("Lx=%f Ly=%f Lz=%f\n", Lx, Ly, Lz);
  printf("type=%d ncontrib=%d\n", type, ncontrib);
#if 0
  rcmysav=0.1;
  rcmxsav=rcmzsav=gamma1sav=gamma2sav=0.0;
  printf("val=%.15G\n", intfunc(0.2, 0.1, 0.2, 0.2));
  exit(-1);
#endif
  //nrfunc = intfunc;
  xtheta = malloc(sizeof(double)*(ntheta+1));
  wtheta = malloc(sizeof(double)*(ntheta+1));
  xphi = malloc(sizeof(double)*(nphi+1));
  wphi = malloc(sizeof(double)*(nphi+1));
  xgamma = malloc(sizeof(double)*(ngamma+1));
  wgamma = malloc(sizeof(double)*(ngamma+1));
  xrcmx = malloc(sizeof(double)*(nrcmx+1));
  wrcmx = malloc(sizeof(double)*(nrcmx+1));
  xrcmy = malloc(sizeof(double)*(nrcmy+1));
  wrcmy = malloc(sizeof(double)*(nrcmy+1));
  xrcmz = malloc(sizeof(double)*(nrcmz+1));
  wrcmz = malloc(sizeof(double)*(nrcmz+1));
  /*printf("x=%.15G %.15G %.15G %.15G %.15G\n w=%.15G %.15G %.15G %.15G %.15G\n",
       (M_PI/2.0-xtheta[1])/(M_PI/2.), xtheta[2]/M_PI, xtheta[3]/M_PI, xtheta[4]/M_PI, ((M_PI/2.)-xtheta[5])/(M_PI/2.0),
       wtheta[1]/(M_PI/2.), wtheta[2]/(M_PI/2.), wtheta[3]/(M_PI/2.), wtheta[4]/(M_PI/2.), wtheta[5]/(M_PI/2.));
  exit(-1);*/
  gauleg(-Lx/2., Lx/2., xrcmx, wrcmx, nrcmx);
  gauleg(-Ly/2., Ly/2., xrcmy, wrcmy, nrcmy);
  gauleg(-Lz/2., Lz/2., xrcmz, wrcmz, nrcmz);
  gauleg(0,2.0*M_PI, xgamma, wgamma ngamma);
  gauleg(0,2.0*M_PI, xphi, wphi nphi);
  gauleg(0, M_PI, xtheta, wtheta, ntheta);
  nsv = -1;
  sobseq(&nsv, sv);
  nsv = 1;
#if 0
  nsv = -1;  
  sobseq(&nsv, sv);
  nsv = 1;
#endif
  
  //printf("XI1[7][8]:%.15G \n", XI1[7][8]);
  /* we use as the reference system the body reference system of first particle */
#ifdef EULER_ROT
  place_DNAD(0.0, 0.0, 0.0, 1., 0., 1., 0., 1., 0., 0);      
#else
  place_DNAD(0.0, 0.0, 0.0, 0., 0., 1., 0., 0);      
#endif
  totfact = 1.0/(2.0*M_PI);
  for (ircmx = 0; ircmx < nrcmx; ircmx++)
    {
      rcmx = xrcmx[ircmx];
      for (ircmy = 0; ircmy < nrcmy; ircmy++)
	{
	  rcmy = xrcmy[ircmy];
	  for (ircmz = 0; itcmz < nrcmz; ircmz++)
	    {
	      rcmz = xrcmz[ircmz];
	      for (igamma = 0; igamma < ngamma; igamma++)
		{
		  gamma = xgamma[igamma];
		  for (iphi = 0; phi < nphi; iphi++)
		    {
		      phi = xphi[iphi]
		      for (itheta = 0; itheta < ntheta; itheta++)
		    	{
			  theta = xtheta[itheta];
		  	  /* place second DNAD randomly */
		  	  if (ircmx > 0 && ircmx % outits == 0)
		  	    {
		  	      printf("building mesh ircmx=%d/%d\n", ircmx, nrcmx)
		  	    }
		  	  if (integrandv1())
		  	    {

		  	    }
		    	}	 
		    }		
		}
	    } 
	}
    }
}
