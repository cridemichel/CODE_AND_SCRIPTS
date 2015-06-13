//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
#define SYMMETRY
#if defined(MPI)
int MPIpid;
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
#endif 
//#define ELEC
//#define ALBERTA
//#define NO_INTERP
#ifdef ELEC
double kD, yukcut, yukcutkD, yukcutkDsq;
#endif
#ifdef PARALLEL
struct kDsortS {
int k1;
int k2;
double invkD;
} *kD_sorted;
double *cdna_arr, *beta_arr;
int numtemps, numconcs;
double **yukcutkD_arr, **kD_arr, **yuk_corr_fact_arr, **yukcutkDsq_arr, **uel_arr, **vexclel_arr;
double num_kD=0;
double maxyukcutkDsq, maxyukcutkD;
#endif
char dummy1[32], dummy2[32], atname[32], nbname[8];
int nat, atnum, nbnum, len;
long long int tot_trials, tt=0, ttini=0;
double L, rx, ry, rz, alpha, dfons_sinth_max, fons_sinth_max;
const double thetapts=100000;

double funcG(double x)
{
  return 1.0;
};

double qgaus(double a, double b)
{
  static const double x[]={0.1488743389816312,0.4333953941292472,
    0.6794095682990244,0.8650633666889845,0.9739065285171717};
  static const double w[]={0.2955242247147529,0.2692667193099963,
    0.2190863625159821,0.1494513491505806,0.0666713443086881};
  double xm, xr, s, dx;
  int j;
  xm=0.5*(b+a);
  xr=0.5*(b-a);
  s=0.;
  for (j=0;j<5;j++) 
    {
      dx=xr*x[j];
      s += w[j]*(funcG(xm+dx)+funcG(xm-dx));
    }
  return s *= xr;
}

#include <math.h> 
#define EPS 1.0e-5
#define JMAX 20 
#define JMAXP (JMAX+1) 
#define K 5
int polinterr=0;
/*Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate; JMAX limits the total number of steps; K is the number of points used in the extrapolation.*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y,
 * and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xai) = yai, 
 * i = 1, . . . , n, then the returned value y = P(x).*/
{ 
  int i,m,ns=1; 
  double den,dif,dift,ho,hp,w;
  double c[K+1], d[K+1];
#if 0
  for (i=0; i < n; i++)
    {
      xa[i+1] = xain[i];
      ya[i+1] = yain[i];
    }
#endif
  dif=fabs(x-xa[1]); 
  //c=vector(n); 
  //d=vector(n); 
  for (i=1;i<=n;i++) 
    { 
      /* Here we find the index ns of the closest table entry,*/
      if ( (dift=fabs(x-xa[i])) < dif) 
	{ 
	  ns=i; 
	  dif=dift;
	} 
      c[i]=ya[i];
      /* and initialize the tableau of c s and d s.*/
      d[i]=ya[i]; 
    } 
  *y=ya[ns--];
  /* This is the initial approximation to y.*/
  for (m=1;m<n;m++) 
    { 
      /* For each column of the tableau,*/
      for (i=1;i<=n-m;i++)
	{
	  /* we loop over the current c s and d s and update them.*/
	  ho=xa[i]-x; 
	  hp=xa[i+m]-x; 
	  w=c[i+1]-d[i]; 
	  if ( (den=ho-hp) == 0.0) 
	    {
	      polinterr=1;
	      return ;
	      //nrerror("Error in routine polint"); 
	    }
	  /* This error can occur only if two input xa s are (to within roundoff)*/
	  den=w/den; d[i]=hp*den; 
	  /*Here the c s and d s are updated. */
	  c[i]=ho*den; 
	} 
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--])); 
      /* After each column in the tableau is completed, we decide which correction, 
       * c or d, we want to add to our accumulating value of y, i.e., which path to take through the tableau 
       * forking up or down. We do this in such a way as to take the most  straight line  route through the 
       * tableau to its apex, updating ns accordingly to keep track of where we are. 
       * This route keeps the partial approximations centered (insofar as possible) on the target x. 
       * The last dy added is thus the error indication. */
    } 
  //free_vector(d); 
  //free_vector(c); 
}
#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
/*This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
as a pointer to the function to be integrated between limits a and b, also input. When called with
n=1, the routine returns the crudest estimate of b f (x)dx. Subsequent calls with n=2,3,... a
(in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.*/
{
  double x,tnm,sum,del; 
  static double s;
  int it,j;
  if (n == 1) 
    {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } 
  else
    {
      for (it=1,j=1;j<n-1;j++) 
	it <<= 1;
      tnm=it;
      del=(b-a)/tnm; /*This is the spacing of the points to be added.*/ 
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) 
	sum += FUNC(x); 
      s=0.5*(s+(b-a)*sum/tnm); /* This replaces s by its refined value. */
      return s;
    } 
}

double qromb(double (*func)(double), double a, double b)
/*Returns the integral of the function func from a to b. Integration is performed by Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.*/
{
  //void nrerror(char error_text[]);
  double ss,dss;
  double s[JMAXP],h[JMAXP+1]; 
  int j;
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd(func,a,b,j); 
    if (j >= K) {
      /* These store the successive trapezoidal approxi- mations and their relative stepsizes.*/
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss)) 
	return ss; 
    }
    h[j+1]=0.25*h[j];
    /*This is a key step: The factor is 0.25 even though the stepsize is decreased by only 0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1), not just a polynomial in h.*/
  }
  printf("Too many steps in routine qromb"); 
  return 0.0; /*Never get here.*/
}

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
#ifdef MC_BENT_DBLCYL
/* apply a random rotation around the supplied axis because 
   bent cylinders do not have azimuthal symmetry */
double thetaGlobalBondangle;
void add_rotation_around_axis(double ox, double oy, double oz, double Rin[3][3], double Rout[3][3], double theta)
{
  double thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random rotation angle between 0 and 2*pi*/
 
  /* set to be used in az. angle distro calculation */
  thetaGlobalBondangle = theta;

  thetaSq=Sqr(theta);
  sinw = sin(theta);
  cosw = (1.0 - cos(theta));
  Omega[0][0] = 0;
  Omega[0][1] = -oz;
  Omega[0][2] = oy;
  Omega[1][0] = oz;
  Omega[1][1] = 0;
  Omega[1][2] = -ox;
  Omega[2][0] = -oy;
  Omega[2][1] = ox;
  Omega[2][2] = 0;
  OmegaSq[0][0] = -Sqr(oy) - Sqr(oz);
  OmegaSq[0][1] = ox*oy;
  OmegaSq[0][2] = ox*oz;
  OmegaSq[1][0] = ox*oy;
  OmegaSq[1][1] = -Sqr(ox) - Sqr(oz);
  OmegaSq[1][2] = oy*oz;
  OmegaSq[2][0] = ox*oz;
  OmegaSq[2][1] = oy*oz;
  OmegaSq[2][2] = -Sqr(ox) - Sqr(oy);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = Rin[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += Rin[k1][k3]*M[k3][k2];
//	  Ro[k1][k2] += M[k1][k3]*Rin[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     Rout[k1][k2] = Ro[k1][k2]; 
}
#endif
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
#ifdef MC_BENT_DBLCYL
  /* add a random rotation around the axis (ox, oy, oz) */
  add_rotation_around_axis(ox, oy, oz, R, Rout, gamma);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R[k1][k2] = Rout[k1][k2];
#endif
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
/* ============================ >>> ranf <<< =============================== */
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}

double fons(double theta, double alpha)
{
  double pi;
  pi = acos(0.0)*2.0;
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return cosh(alpha*cos(theta))*alpha/(4.0*pi*sinh(alpha));
}
/* return an angle theta sampled from an Onsager angular distribution */
double theta_onsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta, dtheta;
  //printf("alpha=%f\n", alpha);
  pi = acos(0.0)*2.0;
  if (first == 1)
    {
      first=0;
#if 0
      f0 = 1.01*fons(0.0,alpha);
#else      
      f0 = 1.01*fons_sinth_max;
#endif
    }
  do 
    {
      /* uniform theta between 0 and pi */
      theta = pi*ranf_vb();
      /* uniform y between 0 and 1 (note that sin(theta) <= 1 for 0 < theta < pi)*/
      y = f0*ranf_vb();
      f = sin(theta)*fons(theta,alpha);
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}
double distro[10000];
const int nfons=100;
void angles_to_R(double *omx, double *omy, double* omz, double alpha)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_onsager(alpha);
  //printf("thos=%f\n", thons);
  distro[(int) (thons/(pi/((double)nfons)))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1;
#if 1 /* along z */
  *omx = verso*sin(thons)*cos(phi);
  *omy = verso*sin(thons)*sin(phi);
  *omz = verso*cos(thons); 
#else /* or along x (but it has to be same of course!) */
  *omy = verso*sin(thons)*cos(phi);
  *omz = verso*sin(thons)*sin(phi);
  *omx = verso*cos(thons); 
#endif
  //printf("norma=%f\n", sqrt(Sqr(*omx)+Sqr(*omy)+Sqr(*omz)));
}

/* first derivative of Onsager distribution */
double dfons(double theta, double alpha)
{
  double pi;
  pi = acos(0.0)*2.0;
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return sinh(alpha*cos(theta))*alpha*alpha/(4.0*pi*sinh(alpha));
}

/* return an angle theta sampled from an Onsager angular distribution */
double theta_donsager(double alpha, int domain)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta, dtheta;
  //printf("alpha=%f\n", alpha);
  pi = acos(0.0)*2.0;
  if (first == 1)
    {
      first=0;
      /* qui va fornito il massimo della funzione nell'intervallo
	 [0,Pi] */
      //f0 = 1.01*alpha*fons(0.0,alpha);
      f0 = 1.01*dfons_sinth_max;
    }
  do 
    {
      /* uniform theta between 0 and pi/2 (domain=0) or pi/2 and pi (domain=1) */
      if (domain==0)
	theta = ranf_vb()*pi/2.;
      else
	theta = (pi/2.)*(1.+ranf_vb());
      /* uniform y between 0 and 1 (note that sin(theta) <= 1 for 0 < theta < pi)*/
      y = f0*ranf_vb();
      f = fabs(sin(theta)*dfons(theta,alpha));
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}

//extern const int nfons;
void orient_donsager(double *omx, double *omy, double* omz, double alpha, int domain)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_donsager(alpha, domain);
  //printf("thos=%f\n", thons);
  //distro[(int) (thons/(pi/((double)nfons)))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1;
#if 1 /* along z */
  *omx = verso*sin(thons)*cos(phi);
  *omy = verso*sin(thons)*sin(phi);
  *omz = verso*cos(thons); 
#else /* or along x (but it has to be same of course!) */
  *omy = verso*sin(thons)*cos(phi);
  *omz = verso*sin(thons)*sin(phi);
  *omx = verso*cos(thons); 
#endif
  //printf("norma=%f\n", sqrt(Sqr(*omx)+Sqr(*omy)+Sqr(*omz)));
}
double estimate_maximum_dfons(double alpha)
{
  double th, dth, maxval, m;
  int i;
  dth=2.0*(acos(0.0))/((double)thetapts);
  th=0.0;
  for (i=0; i < thetapts; i++)
    {
      m=sin(th)*dfons(th,alpha);
      if (i==0 || maxval < m)
	maxval = m;
      th += dth;
      //printf("%f %.15G\n", th, sin(th)*dfons(th, alpha));
    }
  // printf("maxval=%f\n", maxval);
  return maxval;
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
#ifdef ELEC
double delta_rab0=2.0, epsr_prime=1.8, yuk_corr_fact;
double esq_eps, esq_eps_prime; /* = e^2 / (4*pi*epsilon0*epsilon*kB) in J*angstrom */
double esq_eps10, esq_eps_prime10;
const double bmann = 1E-9*0.34/2.0; /* spacing between charged phosphate groups for manning theory */ 
const double Dalton = 1.660538921E-27;
const double kB = 1.3806503E-23, eps0=8.85E-12; /* boltzmann constant */
const double qel = 1.602176565E-19, Nav=6.02214129E23;
const double qdna = 1.0, qsalt = 1.0; /* qsalt è la valenza del sale aggiunto (tipicamente 1 poiché si tratta di NaCl */
double cdna, csalt = 0.0; /* concentrazione del sale aggiunto molare */
double ximanning, deltamann; /* Debye screening length */
/* charge on phosphate groups */
double zeta_a, zeta_b;
double Ucoul(double rab)
{
  //return esq_eps_prime10/rab;
  return esq_eps_prime10*zeta_a*zeta_b/rab;

}
double Uyuk(double rab)
{
#if 0
  printf("qui esq_eps10=%.15G zeta_a=%f exp(-kD*rab)=%.15G rab=%.15G\n", esq_eps10, zeta_a, exp(-kD*rab), rab);
  printf("Uyuk=%.15G\n",esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab );
#endif
  
  return yuk_corr_fact*esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab; 
} 
#ifdef PARALLEL
double Uyuk_arr(double rab)
{
#if 0
  printf("qui esq_eps10=%.15G zeta_a=%f exp(-kD*rab)=%.15G rab=%.15G\n", esq_eps10, zeta_a, exp(-kD*rab), rab);
  printf("Uyuk=%.15G\n",esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab );
#endif
  return yuk_corr_fact*esq_eps10*zeta_a*zeta_b/rab; 
}
double calc_yukawa_arr(int i, int j, double distsq, int *kks)
{
  double ret, rab0, rab, sigab;
  int kk, k1, k2;
  rab = sqrt(distsq);
  sigab = DNADs[0][i].rad + DNADs[1][j].rad;
  rab0 = sigab + delta_rab0; /* we are using Angstrom units here (see pag. S2 in the SI of Frezza Soft Matter 2011) */ 
  *kks=-2;
  if (distsq > maxyukcutkDsq)
    return 0.0;
  for (kk=0; kk < num_kD; kk++)
    {    
      if (rab < rab0)
	{
	  *kks=-1;
	  //printf("interp=%.15G\n",  Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab));
#ifdef NO_INTERP
	  return Ucoul(rab);
#else
	  return Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab);
#endif
#if 0
	  if (isnan(ret))
	    {
	      printf("a=%f rab=%f boh ret=%f Ucoul(rab)=%f Uyuk(rab)=%f\n", zeta_a, rab, ret, Ucoul(rab), Uyuk(rab));
	      exit(-1);
	    }
	  return ret;
#endif    
	}
      /* we set a cutoff for electrostatic interactions */
      else if (rab < yukcut*kD_sorted[kk].invkD)
	{
	  //printf("Yuk=%.15G\n", Uyuk(rab));
	  *kks = kk; 
	  return Uyuk_arr(rab);
	} 
    }
}
#endif
double calc_yukawa(int i, int j, double distsq)
{
  double ret, rab0, rab, sigab;
  rab = sqrt(distsq);
  sigab = DNADs[0][i].rad + DNADs[1][j].rad;
  rab0 = sigab + delta_rab0; /* we are using Angstrom units here (see pag. S2 in the SI of Frezza Soft Matter 2011) */ 
  
  if (rab < rab0)
    {
      //printf("interp=%.15G\n",  Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab));
#ifdef NO_INTERP
      return Ucoul(rab);
#else
      return Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab);
#endif
#if 0
      if (isnan(ret))
	{
	  printf("a=%f rab=%f boh ret=%f Ucoul(rab)=%f Uyuk(rab)=%f\n", zeta_a, rab, ret, Ucoul(rab), Uyuk(rab));
	  exit(-1);
	}
      return ret;
#endif    
    }
  /* we set a cutoff for electrostatic interactions */
  else if (rab < yukcutkD)
    {
      //printf("Yuk=%.15G\n", Uyuk(rab));
      return Uyuk(rab);
    } 
  else return 0.0;
}
#endif
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
#ifdef PARALLEL
int compare_func(const void *aa, const void *bb)
{
  double temp;
  temp = ((struct kDsortS*)aa)->invkD - ((struct kDsortS*)bb)->invkD;
  if (temp < 0)
    return -1;
  else if (temp > 0)
    return 1;
  else
    return 0;
}
#endif
double integrandv1(double rcmx, double rcmy, double rcmz, double theta1, double phi1,
		  double gamma1, double theta2, double phi2, double gamma2, double alpha)
{
  int i, j;
  double sigsq, distsq, sigijsq, u1x, u1y, u1z, u2x, u2y, u2z;
  u1x = sin(theta1)*cos(phi1);
  u1y = sin(theta1)*sin(phi1);
  u1z = cos(theta1); 
  u2x = sin(theta2)*cos(phi2);
  u2y = sin(theta2)*sin(phi2);
  u2z = cos(theta2); 
  //versor_to_R(u1x, u1y, u1z, gamma1, DNADall[0].R);
  //versor_to_R(u2x, u2y, u2z, gamma2, DNADall[1].R);
  place_DNAD(0.0, 0.0, 0.0, u1x, u1y, u1z, gamma1, 0);      
  place_DNAD(rcmx, rcmy, rcmz, u2x, u2y, u2z, gamma2, 1);
  if (calcDistBox() < 0.0)
    {
      for (i=0; i < nat; i++)
	{
	  for (j=0; j < nat; j++)
	    {
	      distsq = Sqr(DNADs[0][i].x-DNADs[1][j].x) + Sqr(DNADs[0][i].y-DNADs[1][j].y) + Sqr(DNADs[0][i].z-DNADs[1][j].z);
	      sigijsq = Sqr(DNADs[0][i].rad + DNADs[1][j].rad);
	      if (distsq < sigsq)
		return u2x*rcmy*sin(theta1)*sin(theta2)*fons(theta1, alpha)*dfons(theta2, alpha);
	    }
	}
    }
  return 0.0;
}
double phi1sav, phi2sav, theta1sav, theta2sav;
double (*nrfunc)(double,double,double,double);
double quad4d(double (*func)(double,double,double,double), 
	      double phi1_1, double phi1_2)
{
  double fphi1(double phi1);
  nrfunc=func;
  return qromb(fphi1,phi1_1,phi1_2);
}
double fphi1(double phi1)
{
  double fphi2(double phi2);
  phi1sav=phi1;
  return qromb(fphi2,0.0,2.0*M_PI); 
}
double fphi2(double phi2) 
{
  double ftheta1(double theta1);
  phi2sav=phi2;
  return qromb(ftheta1,0.0,2.0*M_PI); 
}
double ftheta1(double theta1) 
{
  double ftheta2(double theta2);
  theta1sav=theta1;
  return qromb(ftheta2,0.0,M_PI); 
}
double ftheta2(double theta2) 
{
  return (*nrfunc)(phi1sav,phi2sav,theta1sav,theta2);
}
double rcmxsav, rcmysav, rcmzsav, gamma1sav, gamma2sav, alphasav;
double intfunc(double phi1, double phi2, double theta1, double theta2)
{
  return integrandv1(rcmxsav, rcmysav, rcmzsav, theta1, phi1, gamma1sav, theta2, phi2, gamma2sav, alphasav);
}
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
#ifdef MPI
  MPI_Status status;
#endif
#ifdef ELEC
  double uel, beta;
  int interact;
#endif
  double gamma1, gamma2, Lx, Ly, Lz;
  FILE *fin, *fout, *f, *fread;
#ifdef PARALLEL
  FILE *fp;
  double sigab, rab0, rab0sq, uelcontrib, tempfact;
  int k1, k2, kk;
#endif
  int ncontrib, cc, k, i, j, overlap, type, contrib, cont=0, nfrarg;
  long long int fileoutits, outits;
  char fnin[1024],fnout[256];
  double dummydbl, segno, u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, vexclel=0.0, factor, dth, th;
  /* syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */
#if defined(MPI) 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
  //sprintf(TXT, "rank:%d\n", my_rank);
#endif

  if (argc < 7)
    {
#ifdef ELEC
      printf("syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <fileoutits> [outits] [Temperature (in K)] [DNA concentration in mg/ml] [yukawa cutoff in units of 1/kD] [epsr_prime (1.0-3.0, default=2 ] [delta_rab0 (default=2) ]\n");
#else
      printf("syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <fileoutits> [outits]\n");
#endif
      exit(1);
    }
  strcpy(fnin,argv[1]);
  fin=fopen(fnin,"r");
  len=atoi(argv[2]);
  alpha = atof(argv[3]);
  tot_trials=atoll(argv[4]);
  type = atoi(argv[5]);
  fileoutits = atoll(argv[6]);
  
  if (argc == 7)
    outits=100*fileoutits;
  else
    outits = atoll(argv[7]);

#ifdef ELEC
#ifdef PARALLEL
  if (argc <= 8)
    {
      numtemps=1;
      beta = 1.0;
    }
  else
    {
      if (sscanf(argv[8], "%lf", &dummydbl) < 1)
	{
	  fp=fopen(argv[8],"r");
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      cc++;
	    }
	  beta_arr = malloc(sizeof(double)*cc);
	  rewind(fp);
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      beta_arr[cc] = 1.0/dummydbl;
	      cc++;
	    }
	  fclose(fp);
	  numtemps=cc;
	}
      else
	{
	  numtemps = 1;
	  beta = 1.0/atof(argv[8]);
	} 
    }
  if (argc <= 9)
    {
      numconcs = 1;
      cdna = 600;
    }
  else
    {
      if (!sscanf(argv[9], "%lf", &dummydbl))
	{
	  fp=fopen(argv[9],"r");
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      cc++;
	    }
	  cdna_arr = malloc(sizeof(double)*cc);
	  rewind(fp);
	  numconcs=cc;
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      cdna_arr[cc] = dummydbl;
	      cc++;
	    }
	  fclose(fp);
	}
      else
	{
	  numconcs = 1;
	  cdna = atof(argv[9]);
	}
    }
#else
  if (argc <= 8)
    beta = 1.0;
  else  
    beta = 1.0/atof(argv[8]);

  if (argc <= 9)
    cdna =  600; /* mg/ml */
  else
    cdna = atof(argv[9]);
#endif
  if (argc <= 10)
    yukcut = 2.0;
  else 
    yukcut = atof(argv[10]);

  if (argc <= 11)
    epsr_prime = 2.0;
  else
    epsr_prime = atof(argv[11]);
  if (argc <= 12)
    delta_rab0 = 2.0;
  else
    delta_rab0 = atof(argv[12]);

#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      esq_eps = Sqr(qel)/(4.0*M_PI*eps0)/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
    }
  else
    {
      esq_eps = Sqr(qel)/(4.0*M_PI*eps0*epsr(1.0/beta))/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
    }
#else
  esq_eps = Sqr(qel)/(4.0*M_PI*eps0*epsr(1.0/beta))/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
#endif
  esq_eps10 = esq_eps*1E10;
  esq_eps_prime = Sqr(qel)/(4.0*M_PI*eps0*epsr_prime)/kB;
  esq_eps_prime10 = esq_eps_prime*1E10;
#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      ximanning = esq_eps/bmann;
      deltamann = 1.0/ximanning;
      zeta_a = deltamann;
      zeta_b = deltamann;
      //printf("zeta_a=%f zeta_b:%f\n", zeta_a, zeta_b);
    }
  else
    {
      ximanning = esq_eps*beta/bmann;
      deltamann = 1.0/ximanning;
      zeta_a = deltamann;
      zeta_b = deltamann;
    }
#else
  ximanning = esq_eps*beta/bmann;
  deltamann = 1.0/ximanning;
  zeta_a = deltamann;
  zeta_b = deltamann;
#endif
  /*
     rho_salt =2 csalt Nav 1000;
     rho_counter[cdna_]:=(2 cdna)/(660*Dalton);
     (* cdna in mg/ml e csalt Molare (=moli/litro), nu=numero di cariche per unità di lunghezza *)
     InvDebyeScrLen[T_,qdna_,cdna_,qsalt_,csalt_,\[Epsilon]rel_]:= Sqrt[qdna^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) ( 2cdna)/(660*Dalton)+qsalt^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) 2 csalt Nav 1000 ];
     InvDebyeScrLen[300, 2, 200, 2, 1, 20]^-1*10^9

   */
  /* qdna è la carica rilasciata da ogni gruppo fosfato in soluzione (tipicamente=1) */
#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      kD_arr = malloc(sizeof(double*)*numtemps);
      yukcutkD_arr = malloc(sizeof(double*)*numtemps);
      yukcutkDsq_arr = malloc(sizeof(double*)*numtemps); 
      vexclel_arr = malloc(sizeof(double*)*numtemps); 
      uel_arr = malloc(sizeof(double*)*numtemps);
      yuk_corr_fact_arr = malloc(sizeof(double*)*numtemps); 
      for (k1=0; k1 < numtemps; k1++)
	{
	  kD_arr[k1] = malloc(sizeof(double)*numconcs);
	  yukcutkD_arr[k1] = malloc(sizeof(double)*numconcs);
	  yukcutkDsq_arr[k1] = malloc(sizeof(double)*numconcs);
	  vexclel_arr[k1] = malloc(sizeof(double)*numconcs);
	  uel_arr[k1] = malloc(sizeof(double)*numconcs);
	  yuk_corr_fact_arr[k1] = malloc(sizeof(double)*numconcs); 
	}
      for (k1 = 0; k1 < numtemps; k1++)
	for (k2 = 0; k2 < numconcs; k2++)
	  {
	    kD_arr[k1][k2] = sqrt((4.0*M_PI*esq_eps*(1.0/epsr(1.0/beta_arr[k1])))*beta_arr[k1]*(Sqr(qdna)*2.0*epsr(1.0/beta_arr[k1])*(deltamann/beta_arr[k1])*cdna_arr[k2]*(22.0/24.0)/660.0/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E10;
	    printf("numtemps=%d numconcs=%d kD:%f beta_arr:%f cdna_arr: %f\n", numtemps, numconcs, kD_arr[k1][k2], beta_arr[k1], cdna_arr[k2]);
	    printf("esq_eps: %f qdna=%f deltamann=%f qsalt=%f csalt=%f\n", esq_eps, qdna, deltamann, qsalt, csalt);
	    /* 6.0 Angstrom is the closest distance between phosphate charges */
	    yukcutkD_arr[k1][k2] = yukcut/kD_arr[k1][k2];
	    yukcutkDsq_arr[k1][k2] = Sqr(yukcutkD_arr[k1][k2]);	
	    yuk_corr_fact_arr[k1][k2] = exp(kD_arr[k1][k2]*6.0)/(1.0+kD_arr[k1][k2]*6.0);
	  }
      num_kD = numtemps*numconcs;
      kD_sorted = malloc(sizeof(struct kDsortS)*num_kD);
      cc=0;
      for (k1 = 0; k1 < numtemps; k1++)
	for (k2 = 0; k2 < numconcs; k2++)
	  {
	    kD_sorted[cc].invkD = 1.0/kD_arr[k1][k2];
	    kD_sorted[cc].k1 = k1;
	    kD_sorted[cc].k2 = k2;
	    cc++;
	  }
      qsort(kD_sorted, cc, sizeof(struct kDsortS), compare_func);
      yuk_corr_fact = 1.0;//exp((1.0/kD_sorted[cc-1].invkD)*6.0)/(1.0+(1.0/kD_sorted[cc-1].invkD)*6.0);
      maxyukcutkD = yukcut*kD_sorted[cc-1].invkD;
      maxyukcutkDsq = Sqr(yukcut*kD_sorted[cc-1].invkD);
      printf("min: %f max: %f maxyukcutkD=%f\n",kD_sorted[0].invkD,kD_sorted[cc-1].invkD, sqrt(maxyukcutkDsq));
    }
  else
    {
      kD = sqrt((4.0*M_PI*esq_eps)*beta*(Sqr(qdna)*2.0*deltamann*cdna*(22.0/24.0)/660.0/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E10;
      /* 6.0 Angstrom is the closest distance between phosphate charges */
      yuk_corr_fact = 1.0;//exp(kD*6.0)/(1.0+kD*6.0);
      yukcutkD = yukcut/kD;
      yukcutkDsq = Sqr(yukcutkD);
    }
#else
  kD = sqrt((4.0*M_PI*esq_eps)*beta*(Sqr(qdna)*2.0*deltamann*cdna*(22.0/24.0)/660.0/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E10;
  /* 6.0 Angstrom is the closest distance between phosphate charges */
  yuk_corr_fact = 1.0;//exp(kD*6.0)/(1.0+kD*6.0);
  yukcutkD = yukcut/kD;
  yukcutkDsq = Sqr(yukcutkD);
#endif
#ifdef PARALLEL
  printf("epsr_prime=%f beta=%f deltamanning=%.15G kB=%.15G kD=%.15G (in Angstrom^-1) esq_eps=%.15G esq_eps_prime=%.15G yukcut=%f\n", epsr_prime, beta, deltamann, kB, kD, esq_eps, esq_eps_prime, yukcut);
  printf("yukawa cutoff=%.15G yuk_corr_fact=%.15G\n", yukcutkD, yuk_corr_fact);
#else
  printf("epsr_prime=%f epsr=%f beta=%f deltamanning=%.15G kB=%.15G kD=%.15G (in Angstrom^-1) esq_eps=%.15G esq_eps_prime=%.15G yukcut=%f\n", epsr_prime, epsr(1.0/beta), beta, deltamann, kB, kD, esq_eps, esq_eps_prime, yukcut);
  printf("yukawa cutoff=%.15G yuk_corr_fact=%.15G\n", yukcutkD, yuk_corr_fact);
#endif
#endif
  cont=0;
#ifdef ELEC
  nfrarg = 14;
#else
  nfrarg = 9;
#endif
  if (argc == nfrarg)
    {
      cont=1;
      fread = fopen(argv[nfrarg-1], "r");
      printf("reading file = %s\n", argv[nfrarg-1]);
      while (!feof(fread))
	{
#ifdef ELEC
	  fscanf(fread, "%lld %lf %lf %lf\n", &ttini, &dummydbl, &vexcl, &vexclel);
#else
	  fscanf(fread, "%lld %lf\n", &ttini, &vexcl);
#endif
	}
      fclose(fread);
#ifdef ELEC
      printf("restarting tt=%lld vexcltot=%.15G vexcl=%.15G vexclel=%.15G\n", ttini, vexcl+vexclel, vexcl, vexclel);
#else
      printf("restarting tt=%lld vexcl=%.15G\n", ttini, vexcl);
#endif
    }
  else
    {
      vexcl = 0.0;
#ifdef ELEC
#ifdef PARALLEL
      if (numtemps > 1 || numconcs > 1)
	{
	  for (k1=0; k1 < numtemps; k1++)
	    for (k2=0; k2 < numconcs; k2++)
	      vexclel_arr[k1][k2] = 0.0;
	}
      else
	vexclel = 0.0;
#else
      vexclel = 0.0;
#endif
#endif
      ttini = 0;
    }

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
  dfons_sinth_max=estimate_maximum_dfons(alpha);
  fons_sinth_max=dfons_sinth_max/alpha;
  printf("Estimated maximum of dfons is %f\n", dfons_sinth_max);
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
  nrfunc = intfunc;
  for (tt=ttini+1; tt < tot_trials; tt++)
    {
      /* place second DNAD randomly */
      /* implementare un quasi-MC */
      rcmx = Lx*(drand48()-0.5);
      rcmy = Ly*(drand48()-0.5);
      rcmz = Lz*(drand48()-0.5);
      gamma1 = 2.0*M_PI*drand48();
      gamma2 = 2.0*M_PI*drand48();
      gamma1sav=gamma1;
      gamma2sav=gamma2;
      rcmxsav = rcmx;
      rcmysav = rcmy;
      rcmzsav = rcmz;
      vexcl += quad4d(nrfunc, 0.0, 2.0*M_PI);
      if (tt > 0 && tt % fileoutits == 0)
	{
	  fout = fopen(fnout, "a+");
	  if (type==0)
	    //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	    fprintf(fout,"%lld %.15G\n", tt, Lx*Ly*Lz*vexcl/((double)tt)/1E3);
	  else if (type==1)
	    fprintf(fout,"%lld %.15G\n", tt, (Lx*Ly*Lz*vexcl/((double)tt))*factor/1E4); /* divido per 10^4 per convertire in nm */
	  else
	    fprintf(fout,"%lld %.15G\n", tt, (Lx*Ly*Lz*vexcl/((double)tt))*Sqr(factor)/1E5); /* divido per 10^5 per convertire in nm */
	  fclose(fout);
	}
      if (tt % outits==0)
	printf("trials: %lld/%lld\n", tt, tot_trials);
    } 
}
