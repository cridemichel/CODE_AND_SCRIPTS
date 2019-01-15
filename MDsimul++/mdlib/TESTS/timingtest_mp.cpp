#include <stdlib.h>
#include <stdio.h>
//#include "pmatrix.H"
#define CPOLY
#define AURENTZ
#define BACKSTAB
#ifdef CPOLY
#include "./cpoly.H"
#else
#include "./rpoly.H"
#endif
//#include<complex>
#ifndef NDEG
#define NDEG 6
#endif
#ifdef AURENTZ
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


//#define STATIC
// fino a N=40 la versione statica è più veloce
//By defining both numty and cmplx and supplying cmplx to rpoly class as template parameter (see below)
//one can use: 
//1) mpfr+mpc by defining both numty and cmplx 
//2) cpp_bin_float + cpp_complex
//otherwise one can use mpfr, cpp_bin_float or gmp_float
//with complex<numty> where numty can be mpfr, cpp_bin_float or gmp_float
//The most efficient solution is 1)

double gauss(void)
{
  double  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0; i < 12; i++)
    {
      sum = sum + drand48();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;

}

#define MPC_MP
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using numty = number<cpp_bin_float<100>>;
using cmplx = cpp_complex<100>;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using numty=number<gmp_float<50>>;
using cmplx=complex<numty>;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using numty=number<mpfr_float_backend<16>>;
using cmplx=number<mpc_complex_backend<16>>;
#ifdef BACKSTAB
using bsdbl=number<mpfr_float_backend<500>>;
using bscmplx=number<mpc_complex_backend<500>>;
#endif
#elif defined(FL128)
#include<complex>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using numty=float128;
using cmplx=complex128;
#else
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
#ifdef BACKSTAB
using bsdbl=number<mpfr_float_backend<200>>;
using bscmplx=number<mpc_complex_backend<200>>;
#endif
using numty=double;
using cmplx=complex<numty>;
#endif
#ifdef BACKSTAB
#ifdef CPOLY
void calc_coeff(bscmplx co[], bscmplx er[])
{
  bscmplx c[NDEG+1], alpha;

  int ii, jj;

  for (ii=0; ii < NDEG; ii++)
    {
      c[ii]  = 0.0;
    }
  c[NDEG]=1.0;
  ii=0;
  
  while (ii < NDEG)
    { 
      alpha = -er[ii];
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
  for (ii=0; ii < NDEG; ii++)
     co[ii] = c[NDEG-ii-1];
  co[NDEG]=1.0;
}
#else
void calc_coeff(bsdbl co[], bscmplx er[])
{
  bscmplx c[NDEG+1], alpha;

  int ii, jj;

  for (ii=0; ii < NDEG; ii++)
    {
      c[ii]  = 0.0;
    }
  c[NDEG]=1.0;
  ii=0;
  
  while (ii < NDEG)
    { 
      alpha = -er[ii];
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
  for (ii=0; ii < NDEG; ii++)
     co[ii] = real(c[NDEG-ii-1]);
  co[NDEG]=1.0;
}
#endif
// in questa routine si assume che la radice è reale solo se la parte immaginaria è 0
// ma le radici calcolate possono avere una parte immaginaria piccola 
void calc_coeff_ideal(bsdbl co[], bscmplx er[])
{
  bsdbl rr[NDEG], ir[NDEG], c[NDEG+1], alpha, beta, zero;
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
     co[ii] = c[NDEG-ii-1];
  co[NDEG]=1.0;
}
#ifdef STATIC
#define NN NDEG
#else
#define NN -1
#endif
#ifdef CPOLY
bsdbl calc_backward_err(pvector<cmplx,NN>& roots,  pvector<cmplx,NN>& c)
{
  int i;
  vector<cmplx> rv;
  bscmplx *er = new bscmplx[NDEG];
  bscmplx *cbs = new bscmplx[NDEG+1];
  rv.resize(NDEG);
  for (i=0; i < NDEG; i++)
    {
      rv[i] = cmplx(roots[i]);
    }
  std::sort(rv.begin(),rv.end(),[&](cmplx a, cmplx b)->bool {return real(a) < real(b);}); 

  for (i=0; i < NDEG; i++)
    {
      er[i] = bscmplx(rv[i]);
      //cout << setprecision(16) <<"roots #" << i << " =" << er[i]<< "\n";
    }  
  calc_coeff(cbs, er);
  bsdbl err, errmax;
  for (i=0; i < NDEG+1; i++)
    {
      if (abs(real(c[i]))==0)
        err=abs(cbs[i] - bscmplx(c[i]));
      else
        err=abs((cbs[i] - bsdbl(real(c[i])))/bsdbl(real(c[i])));
      if (err > 100)
        cout << "i=" << i << " cbs=" << cbs[i] << " c=" << c[i] << "\n";
      if (i==0 || err > errmax)
        errmax=err;
   }
#if 0
  bsdbl sum;
  for (i=0; i < NDEG+1; i++)
    sum+=abs(c[i])*abs(c[i]);
  sum=sqrt(sum);
  //cout << "backward error=" << errmax << "\n";
  return errmax/sum;
#else
  return errmax;
#endif
}
#else
bsdbl calc_backward_err(pvector<cmplx,NN> roots,  pvector<numty,NN> c)
{
  int i;
  bscmplx *er = new bscmplx[NDEG];
  bsdbl *cbs = new bsdbl[NDEG+1];
  vector<cmplx> rv;
   
  bsdbl err, errmax;
  rv.resize(NDEG);
  for (i=0; i < NDEG; i++)
    {
      rv[i] = cmplx(roots[i]);
    }
  std::sort(rv.begin(),rv.end(),[&](cmplx a, cmplx b)->bool {return real(a) < real(b);}); 
  for (i=0; i < NDEG; i++)
    {
      //cout << setprecision(16) <<"roots #" << i << " =" << rv[i]<< "\n";
      er[i] = bscmplx(rv[i]); 
    }
  calc_coeff(cbs, er);
  for (i=0; i < NDEG+1; i++)
    {
      if (abs(c[i])==0)
        err=abs(cbs[i] - bsdbl(c[i]));
      else
        err=abs((cbs[i] - bsdbl(c[i]))/bsdbl(c[i]));
      if (i==0 || err > errmax)
        errmax=err;
   }
  //cout << "backward error=" << errmax << "\n";
  return errmax;
}
#endif
#endif
using namespace std;
int main(int argc, char* argv[])
{
  //numty::default_precision(unsigned(50));
  //cmplx::default_precision(unsigned(50));
#ifdef BACKSTAB
  bsdbl berr, berrmax;
  int poltype, m;
#endif
  int i, caso, j, maxiter;
#ifdef STATIC
  rpoly<numty,NDEG,false,cmplx> rp;
  rpoly<numty,NDEG,true,cmplx> rphqr;
  pvector<numty,NDEG+1> c;
  pvector<cmplx,NDEG> roots;
#else
#ifdef CPOLY
  cpoly<cmplx,-1,numty> rp(NDEG);
  cpoly<cmplx,-1,numty> rphqr(NDEG);
  pvector<cmplx,-1> c(NDEG+1);
  pvector<cmplx,-1> roots(NDEG);
#else
  rpoly<numty,-1,false,cmplx> rp(NDEG);
  rpoly<numty,-1,true,cmplx> rphqr(NDEG);
  pvector<numty,-1> c(NDEG+1);
  pvector<cmplx,-1> roots(NDEG);
#endif
#endif
#ifdef AURENTZ
  pvector<double,NDEG+1> ca;
  pvector<complex<double>,NDEG> rootsa;
  int balance=1;
  double alpha;
#endif
  //cmplx x1c, x2c, x3c, x4c, x5c, x6c;
#if 0
c << -0.2269860014469,0.106758402093,-0.02494545844908,0.08966693274224,-0.2613448138378,0.07302387365935,0.1020656417214,0.5145630427322,-0.2961782004896,0.6045784994495,0.07598402268157,0.5646850212872,-0.4609870294448,0.5138572317693,0.3120856456867,0.4563159692286,-0.1124231835099,-0.1705281517195,1.064042629743,0.0663882813567,-0.2026415551499,-0.1901851699741,1.123244995215,0.1317247039573,-0.5084796650622,0.4040805572591,0.8429166795897,0.3121074800984,-0.5361179219351,0.7010215597226,0.5584933046676,0.7533943656272,-0.3972609283811,0.7823959052639,0.604244653732,0.6033151610707,-0.1446604425831,0.4195759742675,1.28744638082,-0.3579521715032,0.7225961738922,-0.1533940066523,1.824029509441,-0.7882870995137,0.6704658200805,0.3968920262438,0.9526890443073,0.09151351510636,-0.627859854702,1.702549871028,0.04725650650517,0.9733198321932,-0.9317249099251,1.459727055467,0.8950835960221,-0.5571578417624,0.6524509030881,0.1825026629376,1.581461868482,-1.447104688938,1.507724529246,-0.4688142341557,1.497523120474,-1.574147291393,0.9029853827528,0.1654747506775,0.4153270868669,-0.4820269951578,-0.3575370185607,1.720501320067,-1.541642561857,1.254911564578,-1.460040768286,2.359610705454,-1.583213062129,0.9807226180953,-0.2439465124685,0.2985870719795,0.9638778804117,-1.988926658131,2.550919395608,-2.170807129614,3.090072930263,-2.83422936892,2.937493913292,-1.232399518157,0.9818935581534,0.565301752405,-0.4606020556921,2.196114430632,-2.395567443036,3.878564226627,-2.807328514903,3.28794873823,-2.659010318394,2.794837141484,-0.5206632037172,0.230691343815,0.5826038352876,-0.3653458678977,2.470979219813,-2.10254322298,2.499255355467,-1.520801284527,2.703763017745,-1.744208039905,1.704652838219,0.4263066052074,0.6216913434341,0.8929244166327,-1.066214842247,3.540422974982,-1.74198313961,2.573665231147,-1.824631544958,3.290070090854,-0.6006656827238,0.7970419705946,0.426977593818,0.6435003086037,1.635603249666,-1.478287370301,2.411691578045,-0.4687210636769,2.228956555628,-1.39496781422,1.527822936876,0.7844391908241,0.7257456300956,0.161152885238,0.3872368656017,2.305495951409,-0.756420500514,1.623087831223,-0.3004932419363,2.627095003688,-0.8761196363053,1.278794383368,0.05451113004596,1.695596107032,0.6616560838342,-0.7495896024239,2.314533555958,-0.03954626389995,2.102287367432,-1.748718692833,3.064203764808,-0.08282908821097,1.481874750239,-0.9769302026346,1.438015309603,1.521177779804,-0.825462606153,0.9799047244901,-0.09806184624147,2.682161007069,-1.496014489989,0.9806768903669,0.4315459002696,1.369133142378,-0.2323470138455,-0.1471718041068,1.797199600697,-0.2804711337398,1.454603096959,-1.324269003413,2.496879145034,-0.8576884559593,1.756893012034,-1.204440707589,1.957245485921,0.4359519489942,0.1524483258536,0.9938083997065,-0.222843159256,2.050109758177,-0.9659765396799,1.57347589888,-0.8677435713656,1.994576675547,-0.6758898706807,1.178342861448,0.3979020199079,0.4889109595084,0.7202793687941,-0.5132990829363,1.483899533789,-0.7709366790157,1.197427587561,-1.069346256398,1.571690312777,-0.5263989751231,0.51063763032,-0.3351629272919,1;
      rp.set_coeff(c);
      //rp.zroots(roots,false);
      //rp.show();
      rp.find_roots(roots);
      //roots.show();  
      exit(-1); 
#endif


#if 0
  rp.allocate(NDEG);
  c.allocate(NDEG+1);
  roots.allocate(NDEG);
#endif
#if 1
#if 0
  cq[5] = 1.0;
  cq[4] = -2.0;
  cq[3] = 5.0;
  cq[2] = -10.0;
  cq[1] = 4.0;
  cq[0] = -8.0;
  qp.set_coeff(cq);
  qp.show("quintic=");
  qp.find_roots(qroots);
  qroots.show("quintic roots=");

  exit(-1);
#endif
#if 0
  c[0] = 274;
  c[1] = -450;
  c[2] = 255;
  c[3] = -60;
  c[4] = 5;
  rp.set_coeff(c);
  rp.find_roots(roots);
  roots.show("boh"); 
#endif
#if 0
  x1c = 1;
  x2c = 2;
  x3c = 8;
  x4c = 9;
  x5c = 10;
  x6c = 100;
  cout << "ntype size=" << sizeof(numty) << "\n";
  c[6] = 1.0;
  c[5] = (-x1c - x2c - x3c - x4c - x5c - x6c).real();
  c[4] = (x1c*x2c + x1c*x3c + x2c*x3c + x1c*x4c + x2c*x4c + x3c*x4c + x1c*x5c + x2c*x5c + x3c*x5c + x4c*x5c + 
   x1c*x6c + x2c*x6c + x3c*x6c + x4c*x6c + x5c*x6c).real(); 
  c[3] = (-(x1c*x2c*x3c) - x1c*x2c*x4c - x1c*x3c*x4c - x2c*x3c*x4c - x1c*x2c*x5c - x1c*x3c*x5c - x2c*x3c*x5c - 
   x1c*x4c*x5c - x2c*x4c*x5c - x3c*x4c*x5c - x1c*x2c*x6c - x1c*x3c*x6c - x2c*x3c*x6c - x1c*x4c*x6c - 
   x2c*x4c*x6c - x3c*x4c*x6c - x1c*x5c*x6c - x2c*x5c*x6c - x3c*x5c*x6c - x4c*x5c*x6c).real();
  c[2] = (x1c*x2c*x3c*x4c + x1c*x2c*x3c*x5c + x1c*x2c*x4c*x5c + x1c*x3c*x4c*x5c + x2c*x3c*x4c*x5c + 
   x1c*x2c*x3c*x6c + x1c*x2c*x4c*x6c + x1c*x3c*x4c*x6c + x2c*x3c*x4c*x6c + x1c*x2c*x5c*x6c + 
   x1c*x3c*x5c*x6c + x2c*x3c*x5c*x6c + x1c*x4c*x5c*x6c + x2c*x4c*x5c*x6c + x3c*x4c*x5c*x6c).real();
  c[1] = (-(x1c*x2c*x3c*x4c*x5c) - x1c*x2c*x3c*x4c*x6c - x1c*x2c*x3c*x5c*x6c - x1c*x2c*x4c*x5c*x6c - 
   x1c*x3c*x4c*x5c*x6c - x2c*x3c*x4c*x5c*x6c).real(); 
  c[0] = (x1c*x2c*x3c*x4c*x5c*x6c).real();
#else
  //c << -1685011.48498632,162908947.86425, 464224691.56007,20202439788.8921,31022318528.7487, -352127.445371458,1;
  // c << 1.0,-1685011.4849863185081630945205688,162908947.864249706268310546875,3464224691.560069561004638671875,20202439788.892055511474609375,31022318528.748729705810546875,-352127.44537145795766264200210571,1;
#endif
#if 0
  rp.set_coeff(c);
  rp.show();
  rp.find_roots(roots);
  roots.show("ows=");  
  rphqr.set_coeff(c);
  rphqr.find_roots(roots);
  roots.show("hqr=");
  exit(-1);
#endif
#if 0
  roots.show("oqs roots=");
  
  rphqr.set_coeff(c);
  rphqr.show();
  rphqr.find_roots(roots);
  roots.show("hqr roots=");
  //cout << "p(x)=" << rp.evaldpoly(1.2) << "\n";
#endif
  //exit(-1);
#endif
#if 0
      bsdbl berrmaxmin;
      int smin;
  for (int s=0; s < 5000; s++)
#endif
  srand48(244);// 244 best
  numty sig=1.0;
  if (argc>=2)
    caso = atoi(argv[1]);
  else
    caso = 0;
  if (caso < 0 || caso > 2)
    {
      cout << "for now you can test only OQS (0) and HQR (1)\n";
      exit(-1);
    }
  if (argc>=3)
    {
      maxiter = atoi(argv[2]);
    }
  else
    maxiter = 1000000;
#ifdef BACKSTAB
  if (argc >= 4)
    poltype = atoi(argv[3]);
  else
    poltype=0;
#endif
  c[NDEG]=1.0;
#if defined(MPC_MP) || defined(CPP_MP) || defined(FL128) || defined(GMP_MP)
  // set precision equal to precision of input coefficients (i.e. double epsilon)
  //rp.set_output_prec(numeric_limits<double>::epsilon());
  rp.set_output_prec(1E-16);
  //rp.set_output_prec(numeric_limits<numty>::epsilon());
  //cout << "qui\n" << " eps=" << numeric_limits<double>::epsilon() << "\n" ;
#else
  //rp.set_output_prec(1E-12);
  //cout << "outputprec=" << rp.get_output_prec() << "\n";
#endif
#ifdef BACKSTAB
  if (poltype !=0 )
    maxiter=1;
#endif
  for (int i=0; i < maxiter; i++)
    {
#ifdef BACKSTAB
#ifdef CPOLY
      for (j=0; j < NDEG; j++)
        c[j]=cmplx(2.0*sig*(drand48()-0.5),0.0);
#else
      if (poltype==0)
        {
          for (j=0; j < NDEG; j++)
            c[j]=sig*(drand48()-0.5);
        }
      //from Aurentz p1, p2 and p3 (see pag. 963 of paper 2015) 
      else if (poltype==1)
        {
          if (NDEG%2==1)
            {
              cout << "You have to set an even polynomial degree for this test\n";
              exit(1);
            }
          for (j=0; j < NDEG; j++)
            c[j] = 0.0;
          c[0] = 1.0;
          c[NDEG/2] = NDEG/(NDEG+1.0)+(NDEG+1.0)/NDEG;
          c[NDEG]=1.0;
        }
      else if (poltype==2)
        {
          if (NDEG%2==1)
            {
              cout << "You have to set an even polynomial degree for this test\n";
              exit(1);
            }
          m=NDEG/2;
          c[m]=m+1.0;
          for (j=0; j < m-1; j++)
            c[j] = m+j;
          for (j=2*m; j >= m+1; j--)
            c[j] = 3*m-j;
          for (j=0; j <= 2*m; j++)
            c[j] /= m;
        }
      else if (poltype==3)
        {
          m=NDEG-1;
          numty lambda=0.9;
          for (j=0; j <= NDEG; j++)
           c[j] = 0;
          c[m+1] = 1.0-lambda;
          c[m] = -(lambda+1.0);
          c[1] = lambda+1.0;
          c[0] = -(1.0-lambda); 
          for (j=1; j <= NDEG; j++)
            c[j] /= c[NDEG];
          c[NDEG]=1.0;
        }
      else if (poltype==4)
        {
          m=NDEG-1;
          numty lambda=0.99;
          for (j=0; j <= NDEG; j++)
           c[j] = 0;
          c[m+1] = 1.0-lambda;
          c[m] = -(lambda+1.0);
          c[1] = lambda+1.0;
          c[0] = -(1.0-lambda);
          for (j=1; j <= NDEG; j++)
            c[j] /= c[NDEG];
          c[NDEG]=1.0;
         }
#endif
#else
#ifdef CPOLY
      for (j=0; j < NDEG; j++)
        c[j]=cmplx(2.0*sig*(drand48()-0.5),0.0);
#else

      for (j=0; j < NDEG; j++)
        c[j]=sig*(drand48()-0.5);
#endif
#endif
      //c << -0.45979133123592,0.016854936689896,0.16284047520515,0.096589904467187,0.49233387428065,0.070947441607572,0.29190001457468,-0.12865167346755,-0.24557117982394,-0.044901894571694,-0.31653706131302,-0.17726595898969,0.3202022896471,0.45510206357928,-0.4201483108709,-0.49967657235846,0.40552818974112,-0.042645962421801,0.086290448952195,0.28437788338989;
      //c << -0.34335362518832,0.3227620422213,-0.33523991717434,0.30603232236776,-0.0004913632795045,0.11895207913477,0.1106512587181,-0.12790756907631,0.48389077297934,-0.0030871046557621,-0.27676517724439,0.21125418248601,0.4490961224484,0.033645447195426,0.27715120795188,-0.012734115584344,0.072796330215358,-0.10955456879709,0.11290799827972,0.083943354978793, 1.0;
      if (caso==0)
        {
          rp.set_coeff(c);
          //rp.zroots(roots, false);
          //rp.show("f(x)=");
          //c.show("coeff=");
          //rp.find_roots(roots, true);
          rp.find_roots_cam(roots);
          roots.set_show_digits(50);
          //roots.show("roots");
          //for (i=0; i < NDEG; i++)
            //cout << setprecision(16) << "roots #" << i << " =" << roots[i]<< "\n";

          //rp.show("poly=");
        }
      else if (caso==1)
        {
          rphqr.set_coeff(c);
          rphqr.find_roots(roots);
          //for (i=0; i < NDEG; i++)
            //cout << setprecision(16) <<"roots #" << i << " =" << roots[i]<< "\n";
        }
#ifdef AURENTZ
      else if (caso==2)
        {
          for (j=0; j < NDEG+1; j++)
            {
              ca[j]=double(real(c[j]));
              //cout << setprecision(16) <<"coeff #" << j << "=" << ca[j] << "\n";
            } 
          if (balance==1)
            wrap_balance(NDEG,ca, &alpha);
          else
            alpha=1.0;
          wrap_damvw(NDEG,ca, rootsa, balance);
          for (j=0; j < NDEG; j++)
            {
              roots[j]=cmplx(alpha*rootsa[j]); 
              //cout << "roots # " << j << " = " << roots[j] << "\n"; 
            }
        }
#endif
      else 
        {
          cout << "Which case?!?\n";
          exit(-1);
        }
#ifdef BACKSTAB
      berr=calc_backward_err(roots, c);
      if (i==0 || berr > berrmax)
        berrmax=berr;
#if 1
      if (berrmax > 100)
        {
          cout << "?!?!?!?\n";
          //rp.show("poly=");
          //for (i=0; i < NDEG; i++)
           // cout << c[i] << ",";
#if 0
          for (i=0; i < NDEG+1; i++)
            {
              cout << setprecision(6) << roots[i] << "\n";
              //cout << setprecision(80) << real(c[i]) << "\n";
            }
#endif
          exit(-1);
        }
#endif
#endif
    }
#if 0
  if (s==0 || berrmax < berrmaxmin) 
    {
      smin=s;
      berrmaxmin=berrmax;
    }
#endif
#ifdef BACKSTAB
  cout << "Max backward error=" << berrmax << "\n";
#endif
#if 0
  printf("root 0: %.15LG %.15LG\n", roots[0].real(), roots[0].imag());
  printf("root 1: %.15LG %.15LG\n", roots[1].real(), roots[1].imag());
  printf("root 2: %.15LG %.15LG\n", roots[2].real(), roots[2].imag());
  printf("root 3: %.15LG %.15LG\n", roots[3].real(), roots[3].imag());
  std::cout.precision(16);
  cout << abs((roots[3]-((cmplx)x2c))/((cmplx)x2c))<<  setprecision(16) << "\n";
  cout << setprecision(16) << abs((roots[0]-((cmplx)x3c))/((cmplx)x3c))<<"\n"; 
  cout << setprecision(16) << abs((roots[2]-((cmplx)x1c))/((cmplx)x1c))<<"\n";
  cout << setprecision(16) << abs((roots[1]-((cmplx)x4c))/((cmplx)x4c))<<"\n";
  roots.show();
#endif

  return 0;
}
