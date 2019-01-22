#include<stdlib.h>
#include<stdio.h>
#define CPOLY
#ifdef CPOLY
#include "../cpoly.H"
#else
#include "../rpoly.H"
#endif
#include <complex>
//#define AURENTZ
//#define MPC_MP
#ifndef NDEG
#define NDEG 6
#endif
//#define STATIC

#ifdef MPSOLVE
extern "C" {
  void mpsolve_wrap(int N);
  void mpsolve_init();
  void mpsolve_free(); 
}
void mpsolve_init_cpp(int n, pvector<double,NDEG+1>& c, double *alpha)
{
  pdbl ca[NDEG+1];
}
void mpsolve_cpp(int n, pvector<double,NDEG+1>& c, double *alpha)
{
  pdbl ca[NDEG+1];


}
void mpsolve_free_cpp()
{

}
#endif
#ifdef AURENTZ
int balance=0;
extern "C" {
extern void damvw_(int *, double*, double *, double *, int *, int*);
extern void balance_(int *N, double *c, int *nnew, double *cn, double *alpha);
}
void wrap_balance(int n, pvector<double,NDEG+1>& c, double *alpha)
{
  double caur[NDEG], caurn[NDEG];
  int i, nnew; 
  //for (i=0; i <= NDEG; i++)
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
  damvw_(&N,caur,rroots,iroots,iterations,&flag);
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
using vldbl = number<cpp_bin_float<500>>;
using cmplx = cpp_complex<500>;
using pdbl=vldbl;
using pcmplx=cmplx;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<gmp_float<500>>;
using cmplx=complex<numty>;
using pdbl=vldbl;
using pcmplx=cmplx;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<mpfr_float_backend<50>>;
using cmplx=number<mpc_complex_backend<50>>;
using pdbl=vldbl;
using pcmplx=cmplx;
#else
using vldbl=long double;
using cmplx=complex<vldbl>;
using pdbl=double;
using pcmplx=complex<pdbl>;
#endif
#if 0
void calc_coeff(pvector<pdbl,NDEG+1>& c, pvector<cmplx,NDEG> er)
{
  cmplx x1c,x2c,x3c,x4c,x5c,x6c;
  x1c=er[0];
  x2c=er[1];
  x3c=er[2];
  x4c=er[3];
  x5c=er[4];
  x6c=er[5];

  c[6] = 1.0;
  c[5] = real(-x1c - x2c - x3c - x4c - x5c - x6c);
  c[4] = real(x1c*x2c + x1c*x3c + x2c*x3c + x1c*x4c + x2c*x4c + x3c*x4c + x1c*x5c + x2c*x5c + 
                x3c*x5c + x4c*x5c + x1c*x6c + x2c*x6c + x3c*x6c + x4c*x6c + x5c*x6c); 
  c[3] = real(-(x1c*x2c*x3c) - x1c*x2c*x4c - x1c*x3c*x4c - x2c*x3c*x4c - x1c*x2c*x5c - 
                x1c*x3c*x5c - x2c*x3c*x5c - x1c*x4c*x5c - x2c*x4c*x5c - x3c*x4c*x5c - 
                x1c*x2c*x6c - x1c*x3c*x6c - x2c*x3c*x6c - x1c*x4c*x6c - x2c*x4c*x6c - 
                x3c*x4c*x6c - x1c*x5c*x6c - x2c*x5c*x6c - x3c*x5c*x6c - x4c*x5c*x6c);
  c[2] = real(x1c*x2c*x3c*x4c + x1c*x2c*x3c*x5c + x1c*x2c*x4c*x5c + x1c*x3c*x4c*x5c + 
                x2c*x3c*x4c*x5c + x1c*x2c*x3c*x6c + x1c*x2c*x4c*x6c + x1c*x3c*x4c*x6c + 
                x2c*x3c*x4c*x6c + x1c*x2c*x5c*x6c + x1c*x3c*x5c*x6c + x2c*x3c*x5c*x6c + 
                x1c*x4c*x5c*x6c + x2c*x4c*x5c*x6c + x3c*x4c*x5c*x6c);
  c[1] = real((-(x1c*x2c*x3c*x4c*x5c) - x1c*x2c*x3c*x4c*x6c - x1c*x2c*x3c*x5c*x6c - 
                 x1c*x2c*x4c*x5c*x6c - x1c*x3c*x4c*x5c*x6c - x2c*x3c*x4c*x5c*x6c)); 
  c[0] = real(x1c*x2c*x3c*x4c*x5c*x6c);
}
#else
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
          alpha = -rr[ii]*vldbl(2.0L);
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
#endif
double ranf(void)
{
  return drand48();
}
#if 1
void sort_sol_opt(pvector<pcmplx,NDEG>& csol, pvector<cmplx,NDEG>& exsol, vldbl allrelerr[])
{
  int k1, k2, k2min;
  int perm[NDEG];
  vldbl relerr, relerrmin;
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
#if 0
      if (k2min==80)
        {
          cout << "k1=" << k1 << "\n";
          cout << "csol=" << csol[k1] << "\n";
          cout << "exsol="<< exsol[k2min]<< "\n";
        }
#endif
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
#else
void sort_sol_opt(pvector<pcmplx,NDEG>& sol, pvector<cmplx,NDEG>& exsol, vldbl *allrelerr)
{
  int k1, k2, k1min;
  vldbl v, vmin;
  cmplx diff;
  pvector<cmplx,NDEG> solt;
  for (k1=0; k1 < 720; k1++)
    {
      v = 0;
      for (k2=0; k2 < NDEG; k2++)
	{
          if (exsol[k2]==cmplx(0,0))
            v += abs(cmplx(sol[perm[k1][k2]])-exsol[k2]);
          else
            v +=abs((cmplx(sol[perm[k1][k2]])-exsol[k2])/exsol[k2]);
        }
      if (k1==0 || v < vmin)
	{
	  k1min=k1;
	  vmin = v;
	}
    } 
  for (k2=0; k2 < NDEG; k2++)
    solt[k2] = cmplx(sol[k2]);

  for (k2=0; k2 < NDEG; k2++)
    sol[k2] = pcmplx(solt[perm[k1min][k2]]);
  for (k2=0; k2 < NDEG; k2++)
    {
      diff = cmplx(sol[k2]) - exsol[k2];
      if (exsol[k2]==cmplx(0.0,0.0))
        allrelerr[k2] = abs(diff);
      else
        allrelerr[k2] = abs(diff/exsol[k2]);
    }
}
#endif
#define MAXSOLV 10
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
int maxic=3, icref;
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
	  fprintf(f, "%.32G %.32G\n", double(vldbl(k)*dlogdE+logdEmin), double(PEall[ic][k]/((double)numtrials)/double(NDEG)));
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
	      cumPEall[ic][k] += PEall[ic][kk]/((double)numtrials)/((double)NDEG);
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

vldbl allrelerr[10][NDEG];
int main(int argc, char **argv)
{
  pdbl  sig, sig2, x1, y1;
  vldbl logdE, dlogdE, logdEmax, logdEmin;
  vldbl dE;
  pvector<cmplx,NDEG> exsol;
  pvector<pcmplx,NDEG> csol;
  pvector<pdbl,NDEG+1> co;
  long long int numtrials, its=0, numout, itsI;
  int numpts, ilogdE;
  int k, ic=0, i;
#ifdef STATIC
  rpoly<pdbl,NDEG> oqs;
  rpoly<pdbl,NDEG,true> hqr;
#else
#ifdef CPOLY
  pvector<pcmplx,-1> csold(NDEG);
  pvector<pcmplx,-1> cod(NDEG+1);
  cpoly<pcmplx,-1,pdbl> oqs(NDEG);
  cpoly<pcmplx,-1,pdbl> hqr(NDEG);
#else
  pvector<pcmplx,-1> csold(NDEG);
  pvector<pdbl,-1> cod(NDEG+1);
  rpoly<pdbl,-1,false,pcmplx> oqs(NDEG);
  rpoly<pdbl,-1,true,pcmplx> hqr(NDEG);
#endif
#endif
  srand48(4242);
  
  for (ic=0; ic < 2; ic++)
    for (k=0; k < PEPTS; k++)
      PEall[ic][k] = 0;

  sig = 1.0;
  sig2= 1.0;
  logdEmax=10.0;
  logdEmin=((int)log10(numeric_limits<vldbl>::epsilon()))-6;
  //cout << "logdEmin=" << logdEmin << "\n";
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
  //oqs.set_output_prec(1E-40);
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
      if (cmplxreal==2) /* all complex */
	{
          for (i=0; i < 2*(NDEG/2); i+=2)
            {
              x1 = sig*(ranf()-0.5);
              y1 = sig*(ranf()-0.5);
              exsol[i] = cmplx(x1, y1);
              exsol[i+1]=cmplx(x1,-y1);
            }
          if (NDEG%2==1)
            {
              exsol[NDEG] = sig2*(ranf()-0.5);
            }
        }
      else if (cmplxreal==1) /* half complex half real */
	{
          
          for (i=0; i < NDEG/2; i=i+2)
            { 
              x1 = sig2*(ranf()-0.5);
              y1 = sig2*(ranf()-0.5);
              exsol[i] = cmplx(x1,y1);
              exsol[i+1] = cmplx(x1,-y1);
            }

          for (; i < NDEG; i++)
            {
              exsol[i] = cmplx(sig*(ranf()-0.5),0);
            }
	}
      else if (cmplxreal==0)/* all real */
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
#ifdef STATIC
          oqs.set_coeff(co);
          oqs.find_roots(csolall[ic]);
#else
#ifdef CPOLY
          for (i=0; i <= NDEG; i++)
            cod[i] = pcmplx(co[i],0.0);
#else
          for (i=0; i <= NDEG; i++)
            cod[i] = co[i];
#endif
          oqs.set_coeff(cod);
          oqs.find_roots(csold);
          //oqs.find_roots_laguerre(csold);
          //oqs.zroots(csold,true);
          for (i=0; i < NDEG; i++)
            csolall[ic][i] = csold[i];
#endif
        }

      ic++;	
      if (dojust==-1 || dojust == ic)
        {
#ifdef STATIC
          hqr.set_coeff(co);
          hqr.find_roots(csolall[ic]);
#else
#ifdef CPOLY
          for (i=0; i <= NDEG; i++)
            cod[i] = pcmplx(co[i],0.0);
#else
           for (i=0; i <= NDEG; i++)
            cod[i] = co[i];
#endif
          hqr.set_coeff(cod);
          hqr.find_roots(csold);
          for (i=0; i < NDEG; i++)
            csolall[ic][i] = csold[i];
#endif
        }
#ifdef AURENTZ
      ic++;	
      if (dojust==-1 || dojust == ic)
	{
          //pvector<pdbl,NDEG+1> ca;
          //pvector<pcmplx,NDEG> rootsa;
          pdbl alpha;
          //for (i=0; i < NDEG+1; i++)
            //ca[i]=co[i]; 
          if (balance==1)
            wrap_balance(NDEG,co, &alpha);
          else
            alpha=1.0L;
          wrap_damvw(NDEG,co, csolall[ic], balance);
          for (i=0; i < NDEG; i++)
            csolall[ic][i]=alpha*csolall[ic][i];
          //  oqs.set_coeff(c);
          //  oqs.find_roots(csolall[ic]);
	}

#endif
      //csolall[2].show("boh");
      //exsol.show("exsol");
      for (ic = 0; ic < maxic; ic++)
	{
	  if (dojust==-1 || dojust==ic)
	    sort_sol_opt(csolall[ic], exsol, allrelerr[ic]);
	}
      for (ic=0; ic < maxic; ic++)
	{
	  if (dojust >= 0 && dojust!=ic)
	    continue;
	  for (k=0; k < NDEG; k++)
	    {	
	     dE = allrelerr[ic][k]; 
             //(exsol[k]!=complex<double>(0,0))?abs((csolall[ic][k] - exsol[k])/exsol[k]):
	     //abs(csolall[ic][k] - exsol[k]); 
             if (dE > 0.0)
		{
		  logdE=log10(dE)-logdEmin;
                  if (log10(dE) >=10)
                    {
                      cout << "exsol[" << k << "]=" << exsol[k] << " csol=" << csolall[ic][k] << "\n";
                      cout << "allrelel=" << allrelerr[ic][k] << "\n";
                      cout << "dE=" << dE << "\n";
                      exit(1);
                    }
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
