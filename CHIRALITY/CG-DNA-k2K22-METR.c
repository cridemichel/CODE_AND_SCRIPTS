//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
#define SYMMETRY
#define ALBERTA
#ifdef ELEC
double kD, yukcut, yukcutkD, yukcutkDsq;
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

double calcDistBox(double shift[3])
{
  double RR, R0, R1, cij[3][3], fabscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k, k1, k2, existsParallelPair = 0;
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  for (k=0; k < 3; k++)
    {
      rA[k] = DNADall[0].rcm[k];
      rB[k] = DNADall[1].rcm[k] + shift[k];
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

char dummy1[32], dummy2[32], atname[32], nbname[8];
int nat, atnum, nbnum, len;
double deltaMC, dthetaMC;
long long int tot_trials, tt=0, ttini=0, tramoveMC=0, rotmoveMC=0;
double L, rx, ry, rz, alpha, dfons_sinth_max, fons_sinth_max;
const double thetapts=100000;
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
void add_rotation_around_axis(double ox, double oy, double oz, double Rin[3][3], double Rout[3][3])
{
  double theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random rotation angle between 0 and 2*pi*/
  theta = 4.0*acos(0.0)*drand48();
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
void orient(double *omx, double *omy, double* omz)
{
  int i;
  //double inert;                 /* momentum of inertia of the molecule */
  //double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, osq, norm;
  
  //Mtot = m; /* total mass of molecule */

  //inert = I; /* momentum of inertia */
 
  //mean = 3.0*temp / inert;

  xisq = 1.0;

  while (xisq >= 1.0)
    {
      xi1  = drand48() * 2.0 - 1.0;
      xi2  = drand48() * 2.0 - 1.0;
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

  *omx = ox;
  *omy = oy;
  *omz = oz; 
  //distro[(int) (acos(oz)/(pi/1000.0))] += 1.0;
}

void versor_to_R(double ox, double oy, double oz, double R[3][3])
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
  add_rotation_around_axis(ox, oy, oz, R, Rout);
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
void place_DNAD(int which)
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
  rO[0] = DNADall[which].rcm[0];
  rO[1] = DNADall[which].rcm[1];
  rO[2] = DNADall[which].rcm[2];
  /* build R here from the orientation (ux,uy,uz) */
  //versor_to_R(ux, uy, uz, R);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	R[k1][k2] = DNADall[which].R[k1][k2];
      }
#ifdef DEBUG 
  sprintf(fn, "DNAD%d.mgl", which);
  fd=fopen(fn, "w+");
  fprintf(fd, ".Vol: %f\n", L*L*L);
#endif
  /* ============ */
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
double min(double a, double b)
{
  return (a < b)?a:b;
}
double distro[10000];
const int nfons=100;
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

void orient_onsager(double *omx, double *omy, double* omz, double alpha)
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
/* return an angle theta sampled from an Onsager angular distribution */
double acc_onsager(double alpha, double theta_old, double theta_new)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  double thr;
  
  thr = min(1.0,sin(theta_new)*fons(theta_new,alpha)/(sin(theta_old)*fons(theta_old,alpha)));
  if (drand48() < thr)
    return 1;
  else 
    return 0;
    
}
double acc_donsager(double alpha, double theta_old, double theta_new)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  double thr;
  
  thr= min(1.0,sin(theta_new)*dfons(theta_new,alpha)/(sin(theta_old)*dfons(theta_old,alpha)));
  if (drand48() < thr)
    return 1;
  else 
    return 0;
    
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
  //return esq_eps_prime10*zeta_a*zeta_b/rab;
  return esq_eps_prime10*zeta_a*zeta_b/rab;

}
double Uyuk(double rab)
{
  return esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab;
} 

double calc_yukawa(int i, int j, double distsq)
{
  double ret, rab0, rab, sigab;
  rab = sqrt(distsq);
  sigab = DNADs[0][i].rad + DNADs[1][j].rad;
  rab0 = sigab + 2; /* we are using Angstrom units here (see pag. S2 in the SI of Frezza Soft Matter 2011) */ 
  
  if (rab < rab0)
    {
      return Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab);
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
void tra_move(int ip)
{
  double dx, dy, dz;
  dx = deltaMC*(drand48()-0.5);
  dy = deltaMC*(drand48()-0.5);
  dz = deltaMC*(drand48()-0.5);
  DNADall[ip].rcm[0] += dx;
  DNADall[ip].rcm[1] += dy;
  DNADall[ip].rcm[2] += dz;
  tramoveMC++; 
}
void rot_move(int ip)
{
  double theta, thetaSq, sinw, cosw;
  double ox, oy, oz, OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random orientation */
#ifdef MC_ALT_ROT
  double thor, xp[3], xl[3];
  thor = 4.0*acos(0.0)*ranf(); /* random angle between 0 and 2*pi */
  xp[0] = 0.0;
  xp[1] = cos(thor);
  xp[2] = sin(thor);
  body2labR(ip, xp, xl, NULL, DNADall[which].R); 
  ox = xl[0];
  oy = xl[1];
  oz = xl[2];
#else
  orient(&ox,&oy,&oz);
    //remove_parall(ip, &ox, &oy, &oz);
#endif
  /* pick a random rotation angle */
  theta= dthetaMC*(drand48()-0.5);
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
	Ro[k1][k2] = DNADall[ip].R[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += DNADall[ip].R[k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     DNADall[ip].R[k1][k2] = Ro[k1][k2]; 
  rotmoveMC++;
}
int random_move(int ip)
{
  double p;
  //printf("random move ip=%d\n", ip);
  p=drand48();
  if (p <= 0.5)
   {
     tra_move(ip);
     return 0;
   }
  else
    {
      rot_move(ip);
      return 1;
    } 
}
double areoverlapping(double shift[3])
{
  int i, j, overlap=0;
  double distsq, sigijsq;
  for (i=0; i < nat; i++)
    {
      for (j=0; j < nat; j++)
	{
	  distsq = Sqr(DNADs[0][i].x-DNADs[1][j].x-shift[0]) + Sqr(DNADs[0][i].y-DNADs[1][j].y-shift[1]) + Sqr(DNADs[0][i].z-DNADs[1][j].z-shift[2]);
	  sigijsq = Sqr(DNADs[0][i].rad + DNADs[1][j].rad);
	  if (distsq < sigijsq)
	    {
	      overlap=1;
	      break;
	    }
#ifdef ELEC
	  /* if they are both phosphate groups we need to calculate electrostatic interaction here */
	  if (DNADs[0][i].atype==1 && DNADs[1][j].atype==1 && distsq < yukcutkDsq)
	    {
#if 0
	      if (distsq < Sqr(yukcut/kD))
			      printf("tt=%lld boh... dist=%f sigij=%f yukcut/kD=%f\n", tt, sqrt(distsq), sqrt(sigijsq), yukcut/kD);
#endif
#if 0
	      if (calc_yukawa(i,j,distsq) < 0.0)
		printf("tt=%lld boh... dist=%f sigij=%f yukcut/kD=%f yuk=%f\n", tt, sqrt(distsq), sqrt(sigijsq), yukcut/kD, calc_yukawa(i,j,distsq));
#endif
	      uel += calc_yukawa(i, j, distsq); 

	    }
#endif
	}
      if (overlap)
	break;
    }
  if (overlap)
    return 1;
  else
    return 0;
}
double calc_theta(int which)
{
  return acos(DNADall[which].R[2][2]);
}
double Rold[3][3], rold[3]; 
void store_state(int ip)
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      rold[k1] = DNADall[ip].rcm[k1];
      for (k2=0; k2 < 3; k2++)
    	Rold[k1][k2] = DNADall[ip].R[k1][k2];
    }
}
void restore_state(int ip)
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      DNADall[ip].rcm[k1] = rold[k1];
      for (k2=0; k2 < 3; k2++)
	DNADall[ip].R[k1][k2] = Rold[k1][k2];
    }
}
int main(int argc, char**argv)
{
#ifdef ELEC
  double uel, beta;
  int interact;
#endif
  int ip, k1, k2, reject, movtype;
  double Lx, Ly, Lz, theta1, theta2, theta_old, theta_new, shift[3];
  FILE *fin, *fout, *f, *fread;
  int ncontrib, cc, k, i, j, overlap, type, contrib, cont=0, nfrarg;
  long long int fileoutits, outits;
  char fnin[1024],fnout[256];
  double dummydbl, segno, u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, vexclel=0.0, factor, dth, th;
  /* syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */
  if (argc < 7)
    {
#ifdef ELEC
      printf("syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <fileoutits> [outits] [Temperature (in K)] [DNA concentration in mg/ml] [yukawa cutoff in units of 1/kD]\n");
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
  if (argc <= 8)
    beta = 1.0;
  else  
    beta = 1.0/atof(argv[8]);

  if (argc <= 9)
    cdna =  600; /* mg/ml */
  else
    cdna = atof(argv[9]);
  
  if (argc <= 10)
    yukcut = 2.0;
  else 
    yukcut = atof(argv[10]);
 
  esq_eps = Sqr(qel)/(4.0*M_PI*eps0*80.1)/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
  esq_eps_prime = Sqr(qel)/(4.0*M_PI*eps0*2.0)/kB;
  esq_eps10 = esq_eps10*1E10;
  esq_eps_prime10 = esq_eps_prime*1E10;
  ximanning = esq_eps*beta/bmann;
  deltamann = 1.0/ximanning;
  zeta_a = deltamann;
  zeta_b = deltamann;
 /*
     rho_salt =2 csalt Nav 1000;
     rho_counter[cdna_]:=(2 cdna)/(660*Dalton);
     (* cdna in mg/ml e csalt Molare (=moli/litro), nu=numero di cariche per unità di lunghezza *)
     InvDebyeScrLen[T_,qdna_,cdna_,qsalt_,csalt_,\[Epsilon]rel_]:= Sqrt[qdna^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) ( 2cdna)/(660*Dalton)+qsalt^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) 2 csalt Nav 1000 ];
     InvDebyeScrLen[300, 2, 200, 2, 1, 20]^-1*10^9

  */
  /* qdna è la carica rilasciata da ogni grupppo fosfato in soluzione (tipicamente=1) */
  kD = sqrt((4.0*M_PI*esq_eps)*beta*(Sqr(qdna)*2.0*cdna*(22.0/24.0)/660.0/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E10;
  yukcutkD = yukcut/kD;
  yukcutkDsq = Sqr(yukcutkD);
  printf("beta=%f deltamanning=%.15G kB=%.15G kD=%.15G (in Angstrom^-1) esq_eps=%.15G esq_eps_prime=%.15G yukcut=%f\n", beta, deltamann, kB, kD, esq_eps, esq_eps_prime, yukcut);
  printf("yukawa cutoff=%.15G\n", yukcutkD);
#endif
  cont=0;
#ifdef ELEC
  nfrarg = 12;
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
      printf("restarting tt=%lld vexcltot=%.15G vexcl=%.15G vexclel=%.15G\n", ttini, vexcl+vexclel, vexclel);
#else
      printf("restarting tt=%lld vexcl=%.15G\n", ttini, vexcl);
#endif
    }
  else
    {
      vexcl = 0.0;
#ifdef ELEC
      vexclel = 0.0;
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
  srand48((int)time(NULL));
  sprintf(fnout, "v%d.dat", type);
  factor=0.0;
  dfons_sinth_max=estimate_maximum_dfons(alpha);
  fons_sinth_max=dfons_sinth_max/alpha;
  printf("Estimated maximum of dfons is %f\n", dfons_sinth_max);
  factor = alpha/2.0;
  printf("factor=%.15G\n", factor);
#if 0
  if (cont)
  { 
#ifdef ELEC
      if (type == 0)
	vexclel *= 1E3*((double)ttini)/(L*L*L);
      else if (type==1)
	vexclel *= 1E4*((double)ttini)/(L*L*L);
      else
	vexclel *= 1E5*((double)ttini)/(L*L*L);
#endif    
      if (type == 0)
	vexcl *= 1E3*((double)ttini)/(L*L*L);
      else if (type==1)
	vexcl *= 1E4*((double)ttini)/(L*L*L);
      else
	vexcl *= 1E5*((double)ttini)/(L*L*L);
    }
  if (!cont)
    {
      fout = fopen(fnout, "w+");
      fclose(fout);
    }
#endif  
   Lx=Ly=Lz=L;
   printf("Lx=%f Ly=%f Lz=%f\n", Lx, Ly, Lz);

   for (k1=0; k1 < 3; k1++)
     {
       DNADall[0].rcm[k1] = DNADall[0].sax[0]/1.001;
       DNADall[1].rcm[k1] = -DNADall[1].sax[0]/1.001;
       for (k2=0; k2 < 3; k2++)
	 {
	   DNADall[0].R[k1][k2] = (k1==k2)?1:0;
	   DNADall[1].R[k1][k2] = (k1==k2)?1:0;
	 } 
     }
   for (tt=ttini; tt < tot_trials; tt++)
    {
      ip=(int) 2.0*drand48();
      store_state(ip);
      theta_old=calc_theta(ip);
      movtype=random_move(ip);
      theta_new=calc_theta(ip); 
      reject=0;
      /* apply periodic boundary conditions */ 
      if (movtype==0)
	{
	   if (DNADall[ip].rcm[0] > Lx*0.5)
	     DNADall[ip].rcm[0] -= Lx;
	   if (DNADall[ip].rcm[0] < -Lx*0.5)
	     DNADall[ip].rcm[0] += Lx;
	   if (DNADall[ip].rcm[1] > Ly*0.5)
	     DNADall[ip].rcm[1] -= Ly;
	   if (DNADall[ip].rcm[1] < -Ly*0.5)
	     DNADall[ip].rcm[1] += Ly;
	   if (DNADall[ip].rcm[2] > Lz*0.5)
	     DNADall[ip].rcm[2] -= Lz;
	   if (DNADall[ip].rcm[2] < -Lz*0.5)
	     DNADall[ip].rcm[2] += Lz;
	}
      shift[0] = Lx*rint((DNADall[0].rcm[0]-DNADall[1].rcm[0])/Lx);
      shift[1] = Ly*rint((DNADall[0].rcm[1]-DNADall[1].rcm[1])/Ly);
      shift[2] = Lz*rint((DNADall[0].rcm[2]-DNADall[1].rcm[2])/Lz);

      if (calcDistBox(shift) > 0.0)
	reject=1;
      else
	{
	  if (type==0)
	    {
	      if (!acc_onsager(alpha, theta_old, theta_new) || !acc_onsager(alpha, theta_old, theta_new))
		reject=1;
	    }
	  else if (type==1)
	    {
	      if (!acc_onsager(alpha, theta_old, theta_new) || !acc_donsager(alpha, theta_old, theta_new))
		reject=1;
	    }  
	  else
	    {
	      if (!acc_donsager(alpha, theta_old, theta_new) || !acc_donsager(alpha, theta_old, theta_new))
		reject=1;
	    }
	}
      if (reject)
	{
	  /* reject move */
	  restore_state(ip);
	  continue;
	}
      place_DNAD(ip);      
      if (areoverlapping(shift))
	{
	  if (type==0)
	    vexcl += 1.0;
	  else if (type==1)
	    vexcl += u2x*rcmy; /* questo '-' rende negativa la k2 e viene dalla derivata della funzione di Onsager! */
	  else 
	    vexcl += -u1x*u2x*rcmy*rcmy;
	}
      if (tt > 0 && tt % fileoutits == 0)
	{
	 fout = fopen(fnout, "a+");
#ifdef ELEC
	 if (type==0)
	   //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	   fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, Lx*Ly*Lz*(vexcl+vexclel)/((double)tt)/1E3, Lx*Ly*Lz*vexcl/((double)tt)/1E3,
		   Lx*Ly*Lz*vexclel/((double)tt)/1E3);
	 else if (type==1)
	   fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel)/((double)tt))*factor/1E4,
		   (Lx*Ly*Lz*vexcl/((double)tt))*factor/1E4,
		   (Lx*Ly*Lz*vexclel/((double)tt))*factor/1E4); /* divido per 10^4 per convertire in nm */
	 else
	   fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel)/((double)tt))*Sqr(factor)/1E5,
		   (Lx*Ly*Lz*vexcl/((double)tt))*Sqr(factor)/1E5,
		   (Lx*Ly*Lz*vexclel/((double)tt))*Sqr(factor)/1E5); /* divido per 10^5 per convertire in nm */
#else 
	 if (type==0)
	   //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	   fprintf(fout,"%lld %.15G\n", tt, vexcl/((double)tt)/1E3);
	 else if (type==1)
	   fprintf(fout,"%lld %.15G\n", tt, vexcl/((double)tt)/1E4); /* divido per 10^4 per convertire in nm */
	 else
	   fprintf(fout,"%lld %.15G\n", tt, vexcl/((double)tt)/1E5); /* divido per 10^5 per convertire in nm */
#endif
	 fclose(fout);
	}
    }
  if (tt % outits==0)
    printf("trials: %lld/%lld\n", tt, tot_trials);
}
