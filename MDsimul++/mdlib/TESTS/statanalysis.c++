#include<stdlib.h>
#include<stdio.h>
#include "../rpoly.H"
#include <complex>
extern int perm[24][4];
int perm[24][4]={
      {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1}, 
      {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0}, 
      {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
      {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};


double ranf(void)
{
  return drand48();
}

void sort_sol_opt(pvector<complex<double>,5>& sol, pvector<complex<double>,5>& exsol)
{
  int k1, k2, k1min;
  double v, vmin;
  pvector<complex<double>,5> solt;
  for (k1=0; k1 < 24; k1++)
    {
      v = 0;
      for (k2=0; k2 < 5; k2++)
	{
	  v += (exsol[k2]==complex<double>(0,0))?abs(sol[perm[k1][k2]]-exsol[k2]):abs((sol[perm[k1][k2]]-exsol[k2])/exsol[k2]);
	}
      if (k1==0 || v < vmin)
	{
	  k1min=k1;
	  vmin = v;
	}
    } 
  for (k2=0; k2 < 5; k2++)
    solt[k2] = sol[k2];

  for (k2=0; k2 < 5; k2++)
    sol[k2] = solt[perm[k1min][k2]];
}
#define MAXSOLV 2
#define PEPTS 100
int cmplxreal=0, restart, dojust=-1;
double cumPEall[MAXSOLV][PEPTS], PEall[MAXSOLV][PEPTS];
pvector<complex<double>,5> csolall[MAXSOLV];
pvector<complex<double>,5> exsol;
char algs[2][64] = {"HQR", "OQS"};
char *ic2algo(int ic)
{
  if (ic > 2)
    {
      exit(-1);
    }
  return algs[ic];
}

void print_legend(FILE *f)
{
  int ic;
  if (dojust < 0)
    {
      for (ic=0; ic < 8; ic++)
	fprintf(f,"@    s%d legend \"%s\"\n", ic, ic2algo(ic));
    }
  else
    fprintf(f,"@    s%d legend \"%s\"\n", dojust, ic2algo(dojust));
}
int maxic=2, icref;
char fname[256];
void save_PE(long long int numtrials, int numpts, double dlogdE, double logdEmin)
{
  FILE *f;
  int k, kk, ic;
  for (ic=0; ic < maxic; ic++)
    {
     if (dojust >= 0 && ic != dojust)
	continue;
      sprintf(fname,"PE-%s.dat", ic2algo(ic));
      f = fopen(fname, "w+");
       for (k=0; k < numpts; k++)
	{
	  if (PEall[ic][k]==0) 
	    continue;
	  fprintf(f, "%.32G %.32G\n", k*dlogdE+logdEmin, PEall[ic][k]/((double)numtrials)/4.);
	}
      fclose(f);
    }
  for (ic=0; ic < maxic; ic++)
    {
      if (dojust >= 0 && ic != dojust)
	continue;
      for (k=0; k < numpts; k++)
	{
	  cumPEall[ic][k] = 0.0;
	  for (kk=k; kk < numpts; kk++) 
	    {
	      cumPEall[ic][k] += PEall[ic][kk]/((double)numtrials)/4.0;
	    }
	}
    }
  for (ic=0; ic < maxic; ic++)
    {
      if (dojust >= 0 && ic != dojust)
	continue;
      
      sprintf(fname,"cumPE-%s.dat", ic2algo(ic));
      f = fopen(fname, "w+");
      for (k=0; k < numpts; k++)
	{
	  if (cumPEall[ic][k]==0)
	    continue;
	  fprintf(f, "%.32G %.32G\n", k*dlogdE+logdEmin, cumPEall[ic][k]);
	}
      fclose(f);
    }
}

int main(int argc, char **argv)
{
  complex<long double> x1c, x2c, x3c, x4c, x5c; 
  double logdE, dlogdE, logdEmax, logdEmin, sig, sig2;
  pvector<double,6> c;
  double dE;
  long double x1,y1;
  long long int numtrials, its, numout, itsI;
  int numpts, ilogdE;
  int num, k, k2, ic=0, okHQR, okHQRL;
  rpoly<double,5> oqs;
  rpoly<double,5,true> hqr;
  srand48(4242);
  
  for (ic=0; ic < 8; ic++)
    for (k=0; k < PEPTS; k++)
      PEall[ic][k] = 0;

  sig = 1.0;
  sig2= 1.0;
  logdEmax=10.0;
  logdEmin=-22.0;
  numpts = PEPTS; 
  dlogdE = (logdEmax -logdEmin)/numpts;

  if (argc>=2)
    numtrials=atoll(argv[1]);
  else 
    numtrials=1000000000;

  restart = 0;
  itsI = 0;
  if (argc>=3)
    numout=atoll(argv[2]);
  else
    numout=100;
  if (argc>=4)
    cmplxreal = atoi(argv[3]);
  if (cmplxreal < 0 || cmplxreal > 5)
    {
      printf("cmplxreal must be between 0 and 5!\n");
      exit(-1);
    }
  if (cmplxreal==3)
    {
      sig = 1.0;
      sig2= 1E6;
      cmplxreal=1;
    }
  else if (cmplxreal==4)
    {
      sig = 1E6;
      sig2 = 1E6;
      cmplxreal = 2;
    } 
  if (argc  >= 5)
    dojust=atoi(argv[4]);
  if (dojust > 8)
    {
      printf("which test should I have to perform?!?\n Last arg is too big!\n");
      exit(-1);
    }
  if (numtrials < 0)
    {
      printf("number of trials must be a positive integer!\n");
      exit(-1);
    }
  else
    {
      printf("numtrials=%lld numout=%lld cmplxreal=%d dojust=%d\n", 
	     numtrials, numout, cmplxreal, dojust);
    }

  x1c=x2c=x3c=x4c=x5c=0;
  for (its=itsI; its < numtrials; its++)
    {
      if (its > 0 && (its % (numtrials/numout) == 0))
	{
          if (cmplxreal == 0)
	    printf("[SAMPLE A sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
          else if (cmplxreal==1)
	    printf("[SAMPLE B sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
	  else if (cmplxreal==2)
	    printf("[SAMPLE C sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
          else if (cmplxreal==3)
	    printf("[SAMPLE D sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
	  else if (cmplxreal==4)
            printf("[SAMPLE E sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
          else if (cmplxreal==5)
            printf("[SAMPLE F sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
	  save_PE(its, numpts, dlogdE, logdEmin);
	}
      /* generate 4 random roots */
      if (cmplxreal==2) /* 4 complex */
	{
	  x1 = sig2*(ranf()-0.5);
	  y1 = sig2*(ranf()-0.5);
	  x1c = complex<long double>(x1, y1);
	  x2c = complex<long double>(x1,-y1);
	  x1 = sig*(ranf()-0.5);
	  y1 = sig*(ranf()-0.5);
	  x3c = complex<long double>(x1,y1);
	  x4c = complex<long double>(x1,-y1);
	}
      else if (cmplxreal==1) /* two complex two real */
	{
	  x1 = sig2*(ranf()-0.5);
	  y1 = sig2*(ranf()-0.5);
	  x1c = complex<long double>(x1,y1);
	  x2c = complex<long double>(x1,-y1);
	  x1 = sig*(ranf()-0.5);
	  y1 = sig*(ranf()-0.5);
	  x3c = x1;
	  x4c = y1;
	}
      else if (cmplxreal==0)/* four real */
	{
	  x1c = sig*(ranf()-0.5);
	  x2c = sig*(ranf()-0.5);
	  x3c = sig*(ranf()-0.5);
	  x4c = sig*(ranf()-0.5);
	}
      
      if (cmplxreal == 5)
	{
	  c[5]=1.0;
          c[4]= ranf()-0.5;
	  c[3]=ranf()-0.5;
	  c[2]=ranf()-0.5;
	  c[1]=ranf()-0.5;
	  c[0]=ranf()-0.5;
	}
      else
	{
          c[5] = 1.0;
          c[4] = (-(x1c+x2c+x3c+x4c+x5c)).real();
          c[3] = (x1c*x2c + x1c*x3c + x2c*x3c + x1c*x4c + x2c*x4c + x3c*x4c + x1c*x5c + 
                  x2c*x5c + x3c*x5c + x4c*x5c).real(); 
          c[2] = (-(x1c*x2c*x3c) - x1c*x2c*x4c - x1c*x3c*x4c - x2c*x3c*x4c - x1c*x2c*x5c - x1c*x3c*x5c - x2c*x3c*x5c - 
                  x1c*x4c*x5c - x2c*x4c*x5c - x3c*x4c*x5c).real();
          c[1] =(x1c*x2c*x3c*x4c + x1c*x2c*x3c*x5c + x1c*x2c*x4c*x5c + 
                 x1c*x3c*x4c*x5c + x2c*x3c*x4c*x5c).real();
          c[0] =(-x1c*x2c*x3c*x4c*x5c).real(); 
	  exsol[0] = x1c;
	  exsol[1] = x2c;
	  exsol[2] = x3c;
	  exsol[3] = x4c;
          exsol[4] = x5c;
        }

      ic = 0;
      if (dojust==-1 || dojust == ic)
        {
          hqr.set_coeff(c);
          hqr.find_roots(csolall[ic]);
        }
      ic++;	
      if (dojust==-1 || dojust == ic)
	{
          oqs.set_coeff(c);
          oqs.find_roots(csolall[ic]);
	}
      for (ic = 0; ic < maxic; ic++)
	{
	  if (dojust==-1 || dojust==ic)
	    sort_sol_opt(csolall[ic], exsol);
	  if (ic==0 && !okHQR)
	    continue;
	}
      for (ic=0; ic < maxic; ic++)
	{
	  if (dojust >= 0 && dojust!=ic)
	    continue;
	  if (ic==0 && !okHQR)
	    continue;
	  for (k=0; k < 5; k++)
	    {	
	     dE = (exsol[k]!=complex<double>(0,0))?abs((csolall[ic][k] - exsol[k])/exsol[k]):
	       abs(csolall[ic][k] - exsol[k]); 
	      if (dE > 0.0)
		{
		  logdE=log10(dE)-logdEmin;
		  ilogdE=(int)(logdE/dlogdE);
		  if (ilogdE >= 0 && ilogdE < numpts)
		    {
		      (PEall[ic][ilogdE])++;
		    }
		}
	    }
	}
    }
  save_PE(numtrials, numpts, dlogdE, logdEmin);
  printf("Finished\n");
  exit(0);
}
