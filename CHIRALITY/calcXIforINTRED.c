//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
#define SYMMETRY
#define USEGSL
#ifdef USEGSL
#include <gsl/gsl_qrng.h>
#endif
//#define ELEC
//#define ALBERTA
//#define NO_INTERP
double **XI1, **XI2, **XI3;
char dummy1[32], dummy2[32], atname[32], nbname[8];
int nat, atnum, nbnum, len;
long long int tot_trials, tt=0, ttini=0;
double L, rx, ry, rz, alpha, dfons_sinth_max, fons_sinth_max, ROMBTOL, phi12, theta12;
double costheta12, sintheta12, cosphi12, sinphi12;
const double thetapts=100000;
double gammln(double xx)
  /*Returns the value ln[Γ(xx)] for xx > 0.*/
{
  /*Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.*/
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5}; 
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp); ser=1.000000000190015;
  for (j=0;j<=5;j++) 
    ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x);
}

#if 0
#define EPSJAC 3.0e-14 /*Increase EPS if you don’t have this precision */ 
#define MAXITJAC 10 
void gaujac(double x[], double w[], int n, double alf, double bet)
/*Given alf and bet, the parameters α and β of the Jacobi polynomials, this routine returns arrays x[1..n] and w[1..n] containing the abscissas and weights of the n-point Gauss-Jacobi quadrature formula. The largest abscissa is returned in x[1], the smallest in x[n].*/
{
  double gammln(double xx);
  /*void nrerror(char error_text[]); */
  int i,its,j;
  double alfbet,an,bn,r1,r2,r3;
  double a,b,c,p1,p2,p3,pp,temp,z,z1;
  for (i=1;i<=n;i++) { 
    if (i == 1) { 
      an=alf/n;
      /* High precision is a good idea for this rou- tine.
	 Loop over the desired roots. Initial guess for the largest root.*/
      bn=bet/n; 
      r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n); 
      r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn; 
      z=1.0-r1/r2;
    } 
    else if (i == 2) 
      { /* Initial guess for the second largest root.*/
	r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf)); 
	r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n; 
	r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
	z -= (1.0-z)*r1*r2*r3;
      } 
    else if (i == 3) 
      { /* Initial guess for the third largest root.*/
	r1=(1.67+0.28*alf)/(1.0+0.37*alf); 
	r2=1.0+0.22*(n-8.0)/n; 
	r3=1.0+8.0*bet/((6.28+bet)*n*n);
	z -= (x[1]-z)*r1*r2*r3;
      } 
    else if (i == n-1) 
      { /* Initial guess for the second smallest root. */ 
	r1=(1.0+0.235*bet)/(0.766+0.119*bet); 
	r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0))); 
	r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
	z += (z-x[n-3])*r1*r2*r3;
      } 
    else if (i == n) 
      { /* Initial guess for the smallest root.*/
	r1=(1.0+0.37*bet)/(1.67+0.28*bet); 
	r2=1.0/(1.0+0.22*(n-8.0)/n); 
	r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n)); 
	z += (z-x[n-2])*r1*r2*r3;
      } 
    else 
      { /* Initial guess for the other roots.*/
       	z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
      }
    alfbet=alf+bet;
    for (its=1;its<=MAXITJAC;its++) 
      {
	temp=2.0+alfbet; p1=(alf-bet+temp*z)/2.0; p2=1.0;
	for (j=2;j<=n;j++) {
	  p3=p2;
	  p2=p1;
	  temp=2*j+alfbet;
	  a=2*j*(j+alfbet)*(temp-2.0);
	  b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z); 
	  c=2.0*(j-1+alf)*(j-1+bet)*temp;
	  p1=(b*p2-c*p3)/a;
	} 
	pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z)); 
	/* p1 is now the desired Jacobi polynomial. We next compute pp, its derivative, 
	   by a standard relation involving also p2, the polynomial of one lower order. */
	z1=z;
	z=z1-p1/pp; /*Newton’s formula.*/
	if (fabs(z-z1) <= EPSJAC) break;
      }
    if (its > MAXITJAC)
      {
	printf("too many iterations in gaujac"); 
	exit(-1);
      }
    x[i]=z; /* Store the root and the weight.*/
    w[i]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
	     gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
  }
}
#endif
int ntheta12, nphi12, nphi1, ntheta1, ngamma1;
double *xtheta12, *xphi12, *wtheta12, *wphi12, *xphi1, *xtheta1, *xgamma1, *wtheta1, *wphi1, *wgamma1;
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
/* nangle=0 phi
         =1 theta 
         =2 gamma */
double qgaus(double (*func)(double), double a, double b, double *x, double *w, int np)
{
#if 0
  static const double x[]={0.1488743389816312,0.4333953941292472,
    0.6794095682990244,0.8650633666889845,0.9739065285171717};
  static const double w[]={0.2955242247147529,0.2692667193099963,
    0.2190863625159821,0.1494513491505806,0.0666713443086881};
#endif
  int j;
  double s;
  s=0.0;
  for (j=1;j<=np;j++) 
    {
      s += w[j]*func(x[j]);
    }
  return s;
}
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}

void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
//#define ALBERTA
char fn[1024];
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

/* ============================ >>> ranf <<< =============================== */
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
double fonsfact, dfonsfact; 

double fons(double costheta, double alpha)
{
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return cosh(alpha*costheta);
}
/* first derivative of Onsager distribution */
double dfons(double costheta, double alpha)
{
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return sinh(alpha*costheta);
}
double func_u2z(double costheta1, double sintheta1, double cosgamma1, double singamma1)
{
  return costheta1*costheta12-cosgamma1*cosphi12*sintheta1*sintheta12+singamma1*sintheta1*sintheta12*sinphi12;
}
double commonfunc_v1(double cosphi1, double sinphi1, double costheta1, double sintheta1,
		     double cosgamma1, double singamma1)
{
  return (costheta12*cosphi1*sintheta1 + cosphi12*sintheta12*(cosgamma1*costheta1*cosphi1 - singamma1*sinphi1) +
    sintheta12*sinphi12*(-costheta1*cosphi1*singamma1 - cosgamma1*sinphi1))*fons(costheta1,alpha)*dfons(func_u2z(costheta1, sintheta1, cosgamma1, singamma1),alpha)*sintheta1*sintheta12;
}
double integrandXI1_v1(double phi1, double theta1, double gamma1)
{
  double cosphi1, sinphi1, cosgamma1, singamma1, costheta1, sintheta1;
  cosphi1=cos(phi1);
  sinphi1=sin(phi1);
  cosgamma1=cos(gamma1);
  singamma1=sin(gamma1);
  costheta1=cos(theta1);
  sintheta1=sin(theta1);
  return commonfunc_v1(cosphi1, sinphi1, costheta1, sintheta1, cosgamma1, singamma1)
    *(cosphi1*singamma1+cosgamma1*costheta1*sinphi1);
}
double integrandXI2_v1(double phi1, double theta1, double gamma1)
{
  double cosphi1, sinphi1, cosgamma1, singamma1, costheta1, sintheta1;
  cosphi1=cos(phi1);
  sinphi1=sin(phi1);
  cosgamma1=cos(gamma1);
  singamma1=sin(gamma1);
  costheta1=cos(theta1);
  sintheta1=sin(theta1);
  return commonfunc_v1(cosphi1, sinphi1, costheta1, sintheta1, cosgamma1, singamma1)*(cosgamma1*cosphi1-costheta1*singamma1*sinphi1);
}
double integrandXI3_v1(double phi1, double theta1, double gamma1)
{
  double cosphi1, sinphi1, cosgamma1, singamma1, costheta1, sintheta1;
  cosphi1=cos(phi1);
  sinphi1=sin(phi1);
  cosgamma1=cos(gamma1);
  singamma1=sin(gamma1);
  costheta1=cos(theta1);
  sintheta1=sin(theta1);
  return commonfunc_v1(cosphi1, sinphi1, costheta1, sintheta1, cosgamma1, singamma1)*(sintheta1*sinphi1);
}

double phi1sav, theta1sav;
double (*nrfunc)(double,double,double);
double ftheta1(double theta1);
double fphi1(double phi1);
double fgamma1(double gamma1);
double quad3d(double (*func)(double,double,double), 
	      double phi1_1, double phi1_2)
{
  nrfunc=func;
  return qgaus(fphi1,phi1_1,phi1_2, xphi1,wphi1,nphi1);
}
double fphi1(double phi1) 
{
  phi1sav=phi1;
  return qgaus(ftheta1,0.0,M_PI, xtheta1, wtheta1, ntheta1); 
}
double ftheta1(double theta1) 
{
  theta1sav=theta1;
  return qgaus(fgamma1,0.0,2.0*M_PI, xgamma1, wgamma1, ngamma1); 
}
double fgamma1(double gamma1)
{
  return (*nrfunc)(phi1sav,theta1sav,gamma1);
}

int main(int argc, char**argv)
{
  char fn[256], fxi1n[256], fxi2n[256], fxi3n[256];
  int aa, bb;
  double cc, xi1, xi2, xi3, totfact;
  double gamma1, gamma2, Lx, Ly, Lz;
  FILE *fin, *fout, *f, *fread, *fxi1, *fxi2, *fxi3;
  int k, i, j, overlap, type;
  long long int fileoutits, outits;
  char fnin[1024],fnout[256];
  double dummydbl, segno, u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, vexclel=0.0, factor, dth, th;
  /* syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */

  if (argc < 5)
    {
      printf("syntax:  calcXIfirINTRED <alpha> <type:0=v0, 1=v1, 2=v2> <nphi> <ntheta>\n");
      exit(1);
    }
  alpha = atof(argv[1]);
  type = atoi(argv[2]);
  nphi12 = atoi(argv[3]);
  ntheta12 = atoi(argv[4]);
  printf("alpha=%f I am going to calculate XI for v%d", alpha, type);
  sprintf(fxi1n, "XI1_v%d.dat", type);
  sprintf(fxi2n, "XI2_v%d.dat", type);
  sprintf(fxi3n, "XI3_v%d.dat", type);
  fxi1 = fopen(fxi1n, "w+");
  fxi2 = fopen(fxi2n, "w+");
  fxi3 = fopen(fxi3n, "w+");
  printf("Gauss quadrature for %d %d points\n", nphi12, ntheta12);
 
  ntheta1 = nphi1 = ngamma1 = 50; 
  xtheta12 = malloc(sizeof(double)*(ntheta12+1));
  xphi12 = malloc(sizeof(double)*(nphi12+1));
  wtheta12 = malloc(sizeof(double)*(ntheta12+1));
  wphi12 = malloc(sizeof(double)*(nphi12+1));

  xtheta1 = malloc(sizeof(double)*(ntheta1+1));
  xphi1 = malloc(sizeof(double)*(nphi1+1));
  xgamma1 = malloc(sizeof(double)*(ngamma1+1));
  wtheta1 = malloc(sizeof(double)*(ntheta1+1));
  wphi1 = malloc(sizeof(double)*(nphi1+1));
  wgamma1 = malloc(sizeof(double)*(ngamma1+1));

  gauleg(0.0, 2.0*M_PI, xphi1, wphi1, nphi1);
  gauleg(0.0, M_PI, xtheta1, wtheta1, ntheta1);
  gauleg(0.0, 2.0*M_PI, xgamma1, wgamma1, ngamma1);

  gauleg(0.0, 2.0*M_PI, xphi12, wphi12, nphi12);
  gauleg(0.0, M_PI, xtheta12, wtheta12, ntheta12);
  
#if 0
  XI1=malloc(sizeof(double)*nphi12);
  XI2=malloc(sizeof(double)*nphi12);
  XI3=malloc(sizeof(double)*nphi12);
  for (i=0; i < nphi12; i++)
    {
      XI1[i] = malloc(sizeof(double)*ntheta12);
      XI2[i] = malloc(sizeof(double)*ntheta12);
      XI3[i] = malloc(sizeof(double)*ntheta12);
    }
#endif
  fprintf(fxi1,"%.15G %d %d\n", alpha, nphi12, ntheta12);
  fprintf(fxi2,"%.15G %d %d\n", alpha, nphi12, ntheta12);
  fprintf(fxi3,"%.15G %d %d\n", alpha, nphi12, ntheta12);

  /* we use as the reference system the body reference system of first particle */
  fonsfact= alpha/(4.0*M_PI*sinh(alpha));
  dfonsfact = alpha*alpha/(4.0*M_PI*sinh(alpha));
  printf("alpha=%f factors=%.15G %.15G sinh(alpha)=%f\n", alpha, fonsfact, dfonsfact, sinh(alpha));
  totfact=fonsfact*dfonsfact/2.0/M_PI;
  for (i=1; i <= nphi12; i++)
    {
      for (j=1; j <= ntheta12; j++)
	{ 
	  phi12 = xphi12[i];
	  theta12 = xtheta12[j];
	  costheta12 = cos(theta12);
	  sintheta12 = sin(theta12);
	  cosphi12 = cos(phi12);
	  sinphi12 = sin(phi12);
	  /* XI1 */
	  xi1=totfact*quad3d(integrandXI1_v1, 0., 2.0*M_PI);
	  /* XI2 */
	  xi2=totfact*quad3d(integrandXI2_v1, 0., 2.0*M_PI);
	  /* XI3 */
	  xi3=totfact*quad3d(integrandXI3_v1, 0., 2.0*M_PI);

	  fprintf(fxi1, "%.15G ", xi1);
	  fprintf(fxi2, "%.15G ", xi2);
	  fprintf(fxi3, "%.15G ", xi3);
	}
      if (i <= nphi12)
	{
	  fprintf(fxi1, "\n");
	  fprintf(fxi2, "\n");
	  fprintf(fxi3, "\n");
	}
    }
  fclose(fxi1);
  fclose(fxi2);
  fclose(fxi3);
}
