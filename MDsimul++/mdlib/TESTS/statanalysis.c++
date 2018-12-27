#include<stdlib.h>
#include<stdio.h>
#include "../rpoly.H"
#include <complex>
#define AURENTZ
#define MPC_MP
#define NDEG 15
#ifdef AURENTZ
extern "C" {
extern void damvw_(int *, double*, double *, double *, int *, int*);
extern void balance_(int *N, double *c, int *nnew, double *cn, double *alpha);
}
void wrap_balance(int n, pvector<double,NDEG+1>& c, double *alpha)
{
  double caur[NDEG], caurn[NDEG];
  int i, nnew; 
  //  for (i=0; i <= NDEG; i++)
  //printf("c[%d]=%.15G\n", i, c[i]);
  for (i=0; i < NDEG; i++)
    {
      caur[i] = c[NDEG-i-1]/c[NDEG];
      //printf("caur[%d]=%.15G\n", i, caur[i]);
    }
  balance_(&n, caur, &nnew, caurn, alpha);
#if 1
  for (i=0; i < NDEG; i++)
    {
      c[i] = caurn[i];
      //printf("roots= %.15G %.15G\n", rroots[i], iroots[i]);
    }
  c[NDEG]=1.0;
#endif
  //for (i=0; i <= NDEG; i++)
    //printf("dopo c[%d]=%.15G\n", i, c[i]);
 
}
void wrap_damvw(int N, pvector<double,NDEG+1> c, pvector<complex<double>,NDEG>& roots, int balanced)
{
  double rroots[NDEG], iroots[NDEG];
  int i, flag, iterations[NDEG];
  double caur[NDEG];
  int n=N;
  //for (i=0; i <= PDEG; i++)
    //printf("c[%d]=%.15G\n", i, c[i]);
  if (balanced==0)
    {
      for (i=0; i < NDEG; i++)
        {
          caur[i] = c[NDEG-i-1]/c[NDEG];
        }
    }
  else
    {
      for (i=0; i < NDEG; i++)
        {
          //cout << "caur=" << caur[i] << "\n";
          caur[i] = c[i];
        }
    }
  damvw_(&n,caur,rroots,iroots,iterations,&flag);
  for (i=0; i < NDEG; i++)
    {
      roots[i] = complex<double>(rroots[i],iroots[i]);
      //printf("roots= %.15G %.15G\n", rroots[i], iroots[i]);
    }
} 
#endif
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl = number<cpp_bin_float<200>>;
using cmplx = cpp_complex<200>;
using pdbl=vldbl;
using pcmplx=cmplx;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<gmp_float<200>>;
using cmplx=complex<numty>;
using pdbl=vldbl;
using pcmplx=cmplx;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<mpfr_float_backend<200>>;
using cmplx=number<mpc_complex_backend<200>>;
using pdbl=double;
using pcmplx=complex<double>;
#else
using vldbl=long double;
using cmplx=complex<vldbl>;
using pdbl=double;
using pcmplx=complex<pdbl>;
#endif

void calc_coeff(pvector<pdbl,NDEG+1>& co, pvector<cmplx,NDEG> er)
{
  vldbl rr[NDEG], ir[NDEG], c[NDEG+1], alpha, beta, zero;
  int ii, jj;

  zero = 0.0;
  for (ii=0; ii < NDEG; ii++)
    {
      rr[ii] = er[ii].real();
      ir[ii] = er[ii].imag();
      c[ii]  = 0.0;
    }
  c[NDEG]=1.0;
  ii=0;
  
  while (ii < NDEG)
    { 
      if (ir[ii] == zero) 
        {
          alpha = -rr[ii];
          for (jj=ii; jj >= 0; jj--)
            {         
              //do jj=ii,1,-1
              if (jj==0)
                c[jj] = c[jj] + alpha;
              else
                c[jj] = c[jj] + alpha*c[jj-1];
            }
          ii=ii+1;
        }
      else
        {
          alpha = -rr[ii]*2.0;
          beta = rr[ii]*rr[ii] + ir[ii]*ir[ii];
          for (jj=ii+1; jj >= 0; jj--)
            { 
              //cout << "jj=" << jj << "\n";
              //do jj=ii+1,1,-1
              if (jj == 1)
                {
                  c[jj] = c[jj] + alpha*c[jj-1] + beta;
                }
              else if (jj == 0) 
                {
                  c[jj] = c[jj] + alpha;
                }
              else 
                c[jj] = c[jj] + alpha*c[jj-1] + beta*c[jj-2];
            }
          ii=ii+2;
        }
    }
  for (ii=0; ii < NDEG; ii++)
     co[ii] = pdbl(c[NDEG-ii-1]);
  co[NDEG]=1.0;
}
double ranf(void)
{
  return drand48();
}
void sort_sol_opt(pvector<pcmplx,NDEG>& csol, pvector<cmplx,NDEG>& exsol, vldbl allrelerr[])
{
  int k1, k2, k2min;
  int perm[NDEG];
  vldbl relerr, relerrmin, relerrmax;
  cmplx diff, solt[NDEG];
  bool used_exsol[NDEG];
  for (k1=0; k1 < NDEG; k1++)
    used_exsol[k1]=false;
  for (k1=0; k1 < NDEG; k1++)
    {
      bool ini = true;
      for (k2=0; k2 < NDEG; k2++)
        {
          if (used_exsol[k2]==true)
            continue;
          diff = cmplx(csol[k1]) - exsol[k2];
          relerr = (exsol[k2]==cmplx(0.0,0.0))?abs(diff):abs(diff/exsol[k2]);
          if (ini==true || relerr <= relerrmin)
           {
             ini=false;
             k2min=k2;
             relerrmin = relerr;
           } 
        }
      perm[k1] = k2min;
      //cout << "perm[" << k1 << "]=" << k2min << "\n";
      allrelerr[k2min] = relerrmin;
      used_exsol[k2min]=true;
    }

  for (k1=0; k1 < NDEG; k1++)
    solt[k1] = cmplx(csol[k1]);

  for (k1=0; k1 < NDEG; k1++)
    csol[perm[k1]] = pcmplx(solt[k1]);
}
#if 0
void sort_sol_opt(pvector<complex<double>,5>& sol, pvector<complex<double>,5>& exsol)
{
  int k1, k2, k1min;
  double v, vmin;
  pvector<complex<double>,5> solt;
  for (k1=0; k1 < 120; k1++)
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
#endif
#define MAXSOLV 2
#define PEPTS 100
int cmplxreal=0, restart, dojust=-1;
double cumPEall[MAXSOLV][PEPTS], PEall[MAXSOLV][PEPTS];
pvector<pcmplx,NDEG> csolall[MAXSOLV];
char algs[3][64] = {"OQS", "HQR", "AUR"};
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
      for (ic=0; ic < 2; ic++)
	fprintf(f,"@    s%d legend \"%s\"\n", ic, ic2algo(ic));
    }
  else
    fprintf(f,"@    s%d legend \"%s\"\n", dojust, ic2algo(dojust));
}
int maxic=2, icref;
char fname[256];
void save_PE(long long int numtrials, int numpts, vldbl dlogdE, vldbl logdEmin)
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
	  fprintf(f, "%.32G %.32G\n", k*dlogdE+logdEmin, double(PEall[ic][k]/((double)numtrials)/4.));
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
	  fprintf(f, "%.32G %.32G\n", double(k*dlogdE+logdEmin), double(cumPEall[ic][k]));
	}
      fclose(f);
    }
}

int main(int argc, char **argv)
{
  vldbl logdE, dlogdE, logdEmax, logdEmin, sig, sig2;
  vldbl dE, x1, y1;
  pvector<cmplx,NDEG> exsol;
  pvector<pcmplx,NDEG> csol;
  pvector<pdbl,NDEG+1> co;
  long long int numtrials, its=0, numout, itsI;
  int numpts, ilogdE;
  int k, ic=0, i;
  rpoly<pdbl,NDEG> oqs;
  rpoly<pdbl,NDEG,true> hqr;
  vldbl allrelerr[NDEG];
  srand48(42);
  
  for (ic=0; ic < 2; ic++)
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
  if (dojust > 2)
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

  for (its=itsI; its < numtrials; its++)
    {
      if (its > 0 && (its % (numtrials/numout) == 0))
	{
          if (cmplxreal == 0)
	    printf("[SAMPLE A sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (cmplxreal==1)
	    printf("[SAMPLE B sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  else if (cmplxreal==2)
	    printf("[SAMPLE C sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (cmplxreal==3)
	    printf("[SAMPLE D sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  else if (cmplxreal==4)
            printf("[SAMPLE E sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (cmplxreal==5)
            printf("[SAMPLE F sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  save_PE(its, numpts, dlogdE, logdEmin);
	}
      /* generate 4 random roots */
      if (cmplxreal==2) /* 4 complex */
	{
          for (i=0; i < NDEG; i++)
            {
              x1 = sig2*(ranf()-0.5);
              y1 = sig2*(ranf()-0.5);
              exsol[i] = cmplx(x1, y1);
              if (i+1==NDEG)
                break;
              exsol[i+1]=cmplx(x1,-y1);
            }
          if (NDEG%2==1)
            {
              exsol[NDEG] = sig2*(ranf()-0.5);
            }
        }
      else if (cmplxreal==1) /* two complex two real */
	{
          for (i=0; i < NDEG; i++)
            { 
              x1 = sig2*(ranf()-0.5);
              y1 = sig2*(ranf()-0.5);
              exsol[i] = cmplx(x1,y1);
              if (i+1==NDEG) 
                break;
              exsol[i+1] = cmplx(x1,-y1);
              if (i+2==NDEG) 
                break;
              exsol[i+2] = cmplx(sig*(ranf()-0.5),0);
              if (i+3==NDEG)
                break;
              exsol[i+3] = cmplx(sig*(ranf()-0.5),0);
            }
	}
      else if (cmplxreal==0)/* four real */
	{
          for (i=0; i < NDEG; i++)
            exsol[i] = cmplx(sig*(ranf()-0.5),0.0);
	}
      if (cmplxreal == 5)
	{
	  co[5]=1.0;
          for (i=0; i < NDEG; i++)
            co[i]= ranf()-0.5;
	}
      else
	{
          calc_coeff(co, exsol);
        }

      ic = 0;
      if (dojust==-1 || dojust == ic)
	{
          oqs.set_coeff(co);
          oqs.find_roots(csolall[ic]);
	}

      ic++;	
      if (dojust==-1 || dojust == ic)
        {
          hqr.set_coeff(co);
          hqr.find_roots(csolall[ic]);
        }
#ifdef AURENTZ
      ic++;	
      if (dojust==-1 || dojust == ic)
	{
          //  oqs.set_coeff(c);
          //  oqs.find_roots(csolall[ic]);
	}

#endif
      for (ic = 0; ic < maxic; ic++)
	{
	  if (dojust==-1 || dojust==ic)
	    sort_sol_opt(csolall[ic], exsol, allrelerr);
	}
      for (ic=0; ic < maxic; ic++)
	{
	  if (dojust >= 0 && dojust!=ic)
	    continue;
	  for (k=0; k < NDEG; k++)
	    {	
	     dE = allrelerr[k]; 
               //(exsol[k]!=complex<double>(0,0))?abs((csolall[ic][k] - exsol[k])/exsol[k]):
	       //abs(csolall[ic][k] - exsol[k]); 
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
