#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_sf_ellint.h>
// il quasi MC è molto più lento del MC nella generazione dei numeri casuali quindi meglio non usarlo qui
//#define SOBOL_LL
#ifdef SOBOL_LL
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))


/*#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))*/
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
#endif
#ifdef SOBOL_LL
double sv[10];
int nsv;
#endif
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
#ifdef SOBOL_LL
  sobseq(&nsv, sv);
  return sv[4];  
#else
  return drand48();
#endif
}
double fons(double theta, double alpha)
{
  double pi;
  //pi = acos(0.0)*2.0;
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return cosh(alpha*cos(theta))*alpha/(4.0*M_PI*sinh(alpha));
}


double theta_onsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta, dtheta;
  //printf("alpha=%f\n", alpha);
  //pi = acos(0.0)*2.0;
  if (first == 1)
    {
      first=0;
      f0 = 1.01*fons(0.0,alpha);
    }

  do 
    {
      /* uniform theta between 0 and pi */
      theta = M_PI*ranf_vb();
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
  //distro[(int) (thons/(pi/((double)nfons)))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1.0;
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
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
double u1[3], u2[3], alpha;
char fn[1024];
int main(int argc, char** argv)
{
  long long int tt = 0, maxtrials, outstps;
  double aini, aend, da, integAn=0.0, integkn=0.0, sp, fact, avgkn, avgAn, avgold, singamma;
  FILE *f;
  fact = M_PI/4.0; 
  if (argc==1)
    {
      printf("calc_ons_int <alphaini> <alphaend> <deltaalpha> <maxtrials> <outits>\n");
      exit(-1);
    }
  aini = atof(argv[1]);
  aend = atof(argv[2]);
  da   = atof(argv[3]);
  maxtrials = atoll(argv[4]);
  if (argc > 5)
    outstps = atoll(argv[5]);
  else
    outstps = 500000;
  printf("aini=%f aend=%f da=%f totsteps=%lld outits=%lld\n", aini, aend, da, maxtrials, outstps);
#ifdef SOBOL_LL
  nsv = -1;  
  sobseq(&nsv, sv);
  nsv = 4;
#endif
  //printf("EllipticE(%f)=%.15G\n", 0.9455, gsl_sf_ellint_Ecomp(sqrt(0.9455), GSL_PREC_DOUBLE));
  //exit(-1);
  f = fopen("An.dat", "w+");
  fprintf(f, "{ ");
  fclose(f);
  f = fopen("kn.dat", "w+");
  fprintf(f, "{ ");
  fclose(f);
  for (alpha=aini; alpha <= aend; alpha+=da)
    {
      sprintf(fn, "onsint_alpha_%f.dat", alpha);
      f = fopen(fn, "w+");
      tt = 0;
      while (tt < maxtrials)
	{
	  orient_onsager(&(u1[0]), &(u1[1]), &(u1[2]), alpha);
	  orient_onsager(&(u2[0]), &(u2[1]), &(u2[2]), alpha);
	  sp = scalProd(u1,u2);
	  singamma = sin(acos(sp));
	  integAn += singamma;
	  /* NOTA: ci va la radice poiché in mathematica l'integrale ellittico è definito in funzione di m=k^2 
	   * dove k è l'argomento dell'integrale ellittico nelle gsl */
	  integkn += 1.0 + fabs(sp) + (4.0/M_PI)*gsl_sf_ellint_Ecomp(sqrt(singamma), GSL_PREC_DOUBLE);
	  //printf("gamma=%f u1=%f %f %f u2=%f %f %f\n", 180.0*acos(sp)/M_PI,
	  //     u1[0], u1[1], u1[2], u2[0], u2[1], u2[2]);
	  //printf("gamma=%f\n", 180.0*acos(sp)/M_PI);
	  if (tt % outstps == 0 && tt != 0) 
	    {
	      avgAn = fact*integAn/((double)tt);
	      avgkn = fact*integkn/((double)tt);
	      printf("An=%.15G kn=%.15G\n", avgAn, avgkn);
	      fprintf(f, "%lld %.15G %.15G\n", tt, avgAn, avgkn);
	    }	  
	  tt++;
	}
      fclose(f);

      f = fopen ("An.dat", "a+");
      if (alpha==aini)
	fprintf(f,"{%.15G,%.15G}", alpha, fact*integAn/((double)tt));
      else
	fprintf(f,",{%.15G,%.15G}", alpha, fact*integAn/((double)tt));
      fclose(f);
      f = fopen ("kn.dat", "a+");
      if (alpha==aini)
	fprintf(f,"{%.15G,%.15G}", alpha, fact*integkn/((double)tt));
      else
	fprintf(f,",{%.15G,%.15G}", alpha, fact*integkn/((double)tt));
      fclose(f);
    }
  f = fopen("An.dat", "a+");
  fprintf(f, " };\n");
  fclose(f);
  f = fopen("kn.dat", "a+");
  fprintf(f, " };\n");
  fclose(f);

}
