#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_sf_ellint.h>
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
int main(int argc, char** argv)
{
  long long int tt = 0, maxtrials, outstps;
  double integAn=0.0, integkn=0.0, sp, fact, avgkn, avgAn, avgold, singamma;
  FILE *f;
  fact = M_PI/4.0; 
  maxtrials = atoll(argv[1]);
  alpha = atof(argv[2]);
  if (argc > 3)
    outstps = atoll(argv[3]);
  else
    outstps = 100000;
  printf("alpha=%f totsteps=%lld outits=%lld\n", alpha, maxtrials, outstps);
  f = fopen("onsint.dat", "w+");

  printf("EllipticE(%f)=%.15G\n", 0.9455, gsl_sf_ellint_Ecomp(sqrt(0.9455), GSL_PREC_DOUBLE));
  exit(-1);
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
  printf("final An=%.15G kn=%.15G\n", fact*integAn/((double)tt), fact*integkn/((double)tt));
}
