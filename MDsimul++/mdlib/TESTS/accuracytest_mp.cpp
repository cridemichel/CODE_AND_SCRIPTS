#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
//#define CPOLY
#ifdef CPOLY
#include "./cpoly.H"
#else
#include "./rpoly.H"
#endif
#include<complex>
#include<list>
#include<string>
#include <iomanip>
static int balance=0;
using namespace std;
//#define AURENTZ
#define MPC_MP
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
//using boost::multiprecision::polar;
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
using cmplx=complex<vldbl>;
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
//template <int N, int digits=200>
//using rpolymp = rpoly<number<mpfr_float_backend<digits>>,N,false,number<mpc_complex_backend<digits>>>;
#define Complex(x,y) cmplx((x),(y))
bool allreal=false, doswap=false;
#ifndef CASO
#define CASO 4
#endif
#undef M_PI
#define M_PI 3.1415926535897932384626433832795029L
using numty = vldbl;
cmplx er[1000];
void calc_coeff(vldbl c[], cmplx er[]);
numty gauss(void)
{
  numty  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  numty sum, r, r2;
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

void calc_coeff_dep_on_case(vldbl c[], cmplx er[])
{
  int i;
#if CASO==1
  // wilkinson
  //
#define NDEG 10
  //static cmplx er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
#elif CASO==2
  // wilkinson
  //
#define NDEG 15
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);

#elif CASO==3
#define NDEG 20
  //static cmplx er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
#elif CASO==4
#define NDEG 20
  allreal=true;
  er[0] = vldbl(-2.1L);
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]+vldbl(0.2L);
    }
  calc_coeff(c, er);
#elif CASO==5
#define NDEG 10
  //static cmplx er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = vldbl(1.0L)/vldbl(i);
    }
  calc_coeff(c, er);
  //r = er;
#elif CASO==6
#define NDEG 15
  //static cmplx er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = vldbl(1.0L)/vldbl(i);
    }
  calc_coeff(c, er);
  //r = er;
#elif CASO==7
#define NDEG 20
  //static cmplx er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = cmplx(1.0)/cmplx(i);
    }
  calc_coeff(c, er);
  //r = er;
#elif CASO==8
#define NDEG 20
  //static cmplx er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = vldbl(1.0L)/pow(vldbl(2),NDEG/2-i);
    }
  calc_coeff(c, er);
  //r = er;
#elif CASO==9
#define NDEG 20
  //static cmplx er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = 1.0L/pow(((vldbl)2),NDEG/2-i)-3.0L;
    }
  calc_coeff(c, er);
  //r = er;
#elif CASO==10
#define NDEG 20
#if 0
  static cmplx er[NDEG];

  for (i=0; i < NDEG; i++)
    { 
      er[i] = cos(M_PI*(2.0L*(i+1L)-1L)/40.0L);
    }
#endif
  allreal=true;
#if 1
  doswap=true;
  static complex<long double> erl[NDEG]=
    {-0.98883082622512854506974288293400861L - 
      0.14904226617617444692935471527721756L*1il, \
        -0.98883082622512854506974288293400861L + 
        0.14904226617617444692935471527721756L*1il, \
        -0.90096886790241912623610231950744505L - 
        0.43388373911755812047576833284835875L*1il, \
        -0.90096886790241912623610231950744505L + 
        0.43388373911755812047576833284835875L*1il, \
        -0.73305187182982632852243148927067191L - 
        0.68017273777091939018735870103374024L*1il, \
        -0.73305187182982632852243148927067191L + 
        0.68017273777091939018735870103374024L*1il, \
        -0.50000000000000000000000000000000000L - 
        0.86602540378443864676372317075293618L*1il, \
        -0.50000000000000000000000000000000000L + 
        0.86602540378443864676372317075293618L*1il, \
        -0.22252093395631440428890256449679476L - 
        0.97492791218182360701813168299393122L*1il, \
        -0.22252093395631440428890256449679476L + 
        0.97492791218182360701813168299393122L*1il, 
      0.07473009358642425429093974573476665L - 
        0.99720379718118014822502987087811927L*1il, 
      0.07473009358642425429093974573476665L + 
        0.99720379718118014822502987087811927L*1il, 
      0.36534102436639501454473799892976880L - 
        0.93087374864420425563779924195127531L*1il, 
      0.36534102436639501454473799892976880L + 
        0.93087374864420425563779924195127531L*1il, 
      0.62348980185873353052500488400423981L - 
        0.78183148246802980870844452667405775L*1il, 
      0.62348980185873353052500488400423981L + 
        0.78183148246802980870844452667405775L*1il, 
      0.82623877431599487194516257377267840L - 
        0.56332005806362202774926153802976051L*1il, 
      0.82623877431599487194516257377267840L + 
        0.56332005806362202774926153802976051L*1il, 
      0.95557280578614073281133405376746667L - 
        0.29475517441090421683077298196019097L*1il, 
      0.95557280578614073281133405376746667L + 
        0.29475517441090421683077298196019097L*1il};  
#endif
  for (i=0; i < NDEG; i++)
    er[i]=cmplx(erl[i]);

  calc_coeff(c, er);
#elif CASO==11
#define NDEG 20
  static complex<long double> erl[NDEG]=
    {-0.98883082622512854506974288293400861L- 
      0.14904226617617444692935471527721756L*1il, 
      -0.98883082622512854506974288293400861L+ 
        0.14904226617617444692935471527721756L*1il, 
      -0.90096886790241912623610231950744505L- 
        0.43388373911755812047576833284835875L*1il, 
      -0.90096886790241912623610231950744505L+ 
        0.43388373911755812047576833284835875L*1il, 
      -0.73305187182982632852243148927067191L- 
        0.68017273777091939018735870103374024L*1il, 
      -0.73305187182982632852243148927067191L+ 
        0.68017273777091939018735870103374024L*1il, 
      -0.50000000000000000000000000000000000L- 
        0.86602540378443864676372317075293618L*1il, 
      -0.50000000000000000000000000000000000L+ 
        0.86602540378443864676372317075293618L*1il, 
      -0.22252093395631440428890256449679476L- 
        0.97492791218182360701813168299393122L*1il, 
      -0.22252093395631440428890256449679476L+ 
        0.97492791218182360701813168299393122L*1il, 
      0.07473009358642425429093974573476665L- 
        0.99720379718118014822502987087811927L*1il, 
      0.07473009358642425429093974573476665L+ 
        0.99720379718118014822502987087811927L*1il, 
      0.36534102436639501454473799892976880L- 
        0.93087374864420425563779924195127531L*1il, 
      0.36534102436639501454473799892976880L+ 
        0.93087374864420425563779924195127531L*1il, 
      0.62348980185873353052500488400423981L- 
        0.78183148246802980870844452667405775L*1il, 
      0.62348980185873353052500488400423981L+ 
        0.78183148246802980870844452667405775L*1il, 
      0.82623877431599487194516257377267840L- 
        0.56332005806362202774926153802976051L*1il, 
      0.82623877431599487194516257377267840L+ 
        0.56332005806362202774926153802976051L*1il, 
      0.95557280578614073281133405376746667L- 
        0.29475517441090421683077298196019097L*1il, 
      0.95557280578614073281133405376746667L+ 
        0.29475517441090421683077298196019097L*1il};
  doswap=true;
  allreal=true;
  for (i=0; i < NDEG; i++)
    er[i]=cmplx(erl[i]);
  
  calc_coeff(c, er);
#elif CASO==12
#define NDEG 24
  doswap=true;
  allreal=true;
 
  static cmplx erl[NDEG];

  erl[0]=Complex(-3.52E2L, 0);
  erl[1]=Complex(-3.52E2L, 0);
  erl[2]=Complex(-2.8371450777E2L, -2.9920517772E2L);
  erl[3]=Complex(-2.8371450777E2L,  2.9920517772E2L);
  erl[4]=Complex(-2.7867414048E2L,  6.1005469197E2L);
  erl[5]=Complex(-2.7867414048E2L, -6.1005469197E2L);
  erl[6]=Complex(-2.74892372E2L, 0L);
  erl[7]=Complex(-2.014171531E2L, 0L);
  erl[8]=Complex(-1.255366582E2L, 0L);
  erl[9]=Complex(-9.599999999E1L, 0L);
  erl[10]=Complex(-8.8692435121E1L,  5.5009607430E2L);
  erl[11]=Complex(-8.869243512E1L, -5.5009607430E2L);
  erl[12]=Complex(-1.6000000000E1L, 0L);
  erl[13]=Complex( 8.23178509855E1L, 0L);
  erl[14]=Complex( 8.8692435121E1L, -5.50096074303E2L);
  erl[15]=Complex( 8.8692435121E1L,  5.5009607430E2L);
  erl[16]=Complex( 1.9293739373E2L,  1.60865921259E3L);
  erl[17]=Complex( 1.929373937E2L, -1.6086592125E3L);
  erl[18]=Complex( 2.0141715312E2L, 0L);
  erl[19]=Complex( 2.7489237213E2L, 0L);
  erl[20]=Complex( 7.52E2L, 0L);
  erl[21]=Complex( 7.52E2L, 0L);
  erl[22]=Complex( 9.1106065E2L,  1.5722L);
  erl[23]=Complex( 9.1106065E2L, -1.5722L);
  static vldbl cs[NDEG+1];
  cs[0]=-54765291428198020791747503747742749163073958404455022926495744.L;
  cs[1]=-4052135566767965847649766745769409681058667331648450681896960.L;
  cs[2]=-31969984081155943263834965670035075493639295858977076674560.L;
  cs[3]=575060225471570237690073740639182419333523437771848417280.L;
  cs[4]=7337981286595499156409929740830030318565357725459415040.L;
  cs[5]=6611223380089859336490797585290455483968982077145088.L;
  cs[6]=-195514288747757987122118583800597358656801082441728.L;
  cs[7]=-726907419403715013562762609680450059293446635520.L;
  cs[8]=197178719520196724204974332265013056299335680.L;
  cs[9]=5968852409133617129605588058090797893943296.L;
  cs[10]=16576506891508825500182005531742679597056.L;
  cs[11]=23375026506968330494765978581548924928.L;
  cs[12]=2206941937668751746514177591607296.L;
  cs[13]=-75617855277818001758431020580864.L;
  cs[14]=-204797687173976372829472423936.L;
  cs[15]= -143150263927579584306872320.L;
  cs[16]=  20214880144364480233472.L;
  cs[17]=  453786251090072698880.L;
  cs[18]=  1265052493274939392.L;
  cs[19]= -968887355572224.L;
  cs[20]=  1015406084096.L;
  cs[21]= -3949133824.L;
  cs[22]=  3284992.L;
  cs[23]= -1728.L;

  cs[24]=1.0L;

  for (i=0; i < NDEG+1; i++)
    {
      c[i] = cs[i];
    }
  for (i=0; i < NDEG; i++)
    er[i]=erl[i];
  

#elif CASO==13


#elif CASO==14


#elif CASO==15
// Noferini
#define NDEG 12
  //roots and coefficients were calculated by Wolfram Mathematica with a precision of 1000 digits
  //print*, "Vanni Noferini's example degree 12 or 35"
  for ( i=0; i < NDEG; i++)
    {
      er[i]=gauss();
    }
  er[0] *= 1E9;
  er[1] *=1E12;

  calc_coeff(c, er);

#if 0
  for (i=0; i < NDEG+1; i++)
    {
      c[i] = cs[i];
    }
  
#endif
#elif CASO==16
// Noferini
#define NDEG 35
  //roots and coefficients were calculated by Wolfram Mathematica with a precision of 1000 digits
  //print*, "Vanni Noferini's example degree 12 or 35"
  for ( i=0; i < NDEG; i++)
    {
      er[i]=gauss();
    }
  er[0] *= 1E9;
  er[1] *=1E12;

  calc_coeff(c, er);


#define NDEG 35
 
#elif CASO==21
#define NDEG 10
  allreal=true;
  er[0]=vldbl(0.1);
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/vldbl(10.0L);
    }
  calc_coeff(c, er);
#if 0
  for (i=0; i < NDEG; i++)
    cout << "[CALC] c[" << i << "]=" << c[i] << setprecision(20) << "\n";
  c[0] = 1.00000000000000053791679056187E-55;
  c[1] = 1.11111111100000055584244679722E-45;
  c[2] = 1.12233445443322162787800313739E-36;
  c[3] = 1.12345790111098789834840761654E-28;
  c[4] = 1.12357014577977569473069186308E-21;
  c[5] = 1.12358025801221010756462539505E-15;
  c[6] = 1.12357014577977572218174173241E-10;
  c[7] = 1.12345790111098777138700591932E-6;
  c[8] = 0.0011223344544332213099796513589;
  c[9] = 0.111111111100000006790544659907;
  c[10] = 1.0;
#endif
  //for (i=0; i < NDEG; i++)
  //cout << " c[" << i << "]=" << c[i] << setprecision(20) << "\n";

#elif CASO==22
#define NDEG 20
  allreal=true;
  er[0]=vldbl(0.1);
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/vldbl(10.0L);
    }
  calc_coeff(c, er);
#if 0
  for (i=0; i < NDEG; i++)
    cout << "[CALC] c[" << i << "]=" << c[i] << setprecision(20) << "\n";
  c[0] = 1.00000000000000053791679056187E-55;
  c[1] = 1.11111111100000055584244679722E-45;
  c[2] = 1.12233445443322162787800313739E-36;
  c[3] = 1.12345790111098789834840761654E-28;
  c[4] = 1.12357014577977569473069186308E-21;
  c[5] = 1.12358025801221010756462539505E-15;
  c[6] = 1.12357014577977572218174173241E-10;
  c[7] = 1.12345790111098777138700591932E-6;
  c[8] = 0.0011223344544332213099796513589;
  c[9] = 0.111111111100000006790544659907;
  c[10] = 1.0;
#endif
  //for (i=0; i < NDEG; i++)
  //cout << " c[" << i << "]=" << c[i] << setprecision(20) << "\n";
#elif CASO==26
#define NDEG 60
  //allreal=true;
  int ii;
  int m=15;
#ifdef MPC_MP
  vldbl pi = boost::math::constants::pi<vldbl>();//2.0*acos(vldbl(0.0));
#else
  vldbl pi=M_PI;
#endif
  for  (ii=-m+1; ii <= 0; ii++)
    er[m-1+ii]= vldbl(0.9L)*exp(cmplx(0,vldbl(ii)*pi/vldbl(2.0L)/vldbl(m)));
  for  (ii=1; ii <= m; ii++)
    er[m-1+ii] = conj(er[m-ii]);
  for  (ii=m+1; ii <= 2*m; ii++)
    er[m-1+ii] = exp(cmplx(0,vldbl(ii)*pi/vldbl(2.0L)/vldbl(m)));
  for  (ii=2*m+1; ii <= 3*m; ii++)
    er[m-1+ii] = conj(er[5*m-ii]);
  calc_coeff(c, er);
#if 0
  for (i=0; i < NDEG; i++)
    cout << "[CALC] c[" << i << "]=" << c[i] << setprecision(20) << "\n";
  c[0] = 1.00000000000000053791679056187E-55;
  c[1] = 1.11111111100000055584244679722E-45;
  c[2] = 1.12233445443322162787800313739E-36;
  c[3] = 1.12345790111098789834840761654E-28;
  c[4] = 1.12357014577977569473069186308E-21;
  c[5] = 1.12358025801221010756462539505E-15;
  c[6] = 1.12357014577977572218174173241E-10;
  c[7] = 1.12345790111098777138700591932E-6;
  c[8] = 0.0011223344544332213099796513589;
  c[9] = 0.111111111100000006790544659907;
  c[10] = 1.0;
#endif
  //for (i=0; i < NDEG; i++)
  //cout << " c[" << i << "]=" << c[i] << setprecision(20) << "\n";

#elif CASO==31
  //Polynomials with few very clustered roots.
  //Kameny  
#define NDEG 9
  numty K = 1E50;
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[0] = 9.0;
  c[2] = -6.0*K*K;
  c[4] = K*K*K*K;
  c[9] = K*K; 
  //r = NULL;
#elif CASO==32
#define NDEG 10
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[NDEG] = 1.0;
  c[0]=1.0;
  c[1] = -300.0;
  c[2] = 30000.0;
  c[3] = -1E6;
  //r=NULL;
#elif CASO==5
#define NDEG 11
  static numty ct ={-1E22, 2E21, -1E20, 0, 0, 0, 0, 0, 0, 0, 1.0};
  for (i=0; i < 11; i++)
    c[i] = ct[i]; 
  r=NULL;
#elif CASO==33
#define NDEG 20
  static cmplx er[NDEG];
  er[0]=1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/10.0L;
    }
  calc_coeff(c, er);
  *r = er;
#endif
}
#ifdef AURENTZ
extern "C" {
extern void damvw_(int *, double*, double *, double *, int *, int*);
extern void balance_(int *N, double *c, int *nnew, double *cn, double *alpha);
}
void wrap_balance(int n, pvector<double,NDEG+1>& c, double *alpha)
{
  double caur[NDEG], caurn[NDEG];
  int i, nnew; 
  for (i=0; i <= NDEG; i++)
    printf("c[%d]=%.15G\n", i, c[i]);
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
  for (i=0; i <= NDEG; i++)
    printf("dopo c[%d]=%.15G\n", i, c[i]);
 
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

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
void sort_sol_opt(cmplx csol[], cmplx exsol[], vldbl allrelerr[])
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
          diff = csol[k1] - exsol[k2];
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
    solt[k1] = csol[k1];

  for (k1=0; k1 < NDEG; k1++)
    csol[perm[k1]] = solt[k1];
}
numty print_accuracy_at(char *str, cmplx csol[], cmplx exsol[], vldbl allrelerr[])
{
  /* we follow FLocke here */
  int k1;
  vldbl relerrmax;
  for (k1=0; k1 < NDEG; k1++)
    {
      if (k1==0 || allrelerr[k1] > relerrmax)
        {
          relerrmax=allrelerr[k1];
        }
    }

  //printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
  cout << setprecision(16) << "[" << str << "relative accuracy=" <<  relerrmax << "\n";
  return relerrmax;
}

void print_roots(char *str, cmplx er[], cmplx cr[], vldbl allrelerr[])
{
  printf("CASE %s\n", str);
  for (auto i=0; i < NDEG; i++)
    {
      cout << setprecision(16) << "root #" << i << " EX: "<< er[i] << " C:" << cr[i];
      cout << setprecision(16) << " [ eps: " << allrelerr[i] << " ]\n"; 
    }
}

void print( list<int> l){
    for(list<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}


void calc_coeff(vldbl co[], cmplx er[])
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
     co[ii] = c[NDEG-ii-1];
  co[NDEG]=1.0;
}
#if 0
void subset(list<list<int>>& L, int arr[], int size, int left, int index, list<int> &l)
{
    if(left==0)
      {
        L.push_back(l);
        //print(L.back());
        return;
      }
    for(int i=index; i<size;i++)
      {
        l.push_back(arr[i]);
        subset(L,arr,size, left-1,i+1,l);
        l.pop_back();
      }
} 

void calc_coeff(vldbl c[], cmplx er[])
{
  int i, j;
  std::list<std::list<int>> subsets;// list containing all subsets of a given length
  cmplx cc[NDEG+1];
  cmplx term, segno;
  // list with all list of subsets
  list<list<int>> lt;   
  // a subset of integers
  list<int> l;
  int array[NDEG];
  for (j=0; j < NDEG; j++)
    array[j] = j;
  cc[NDEG] = 1.0L;
  for (i=0; i < NDEG; i++)
    {
      // all subsets of NDEG-i elements are stored in list lt
      subset(lt,array,NDEG,NDEG-i,0,l);
      cc[i] = 0L;
      // loop over all subsets
      for(list<list<int>>::iterator it=lt.begin(); it!=lt.end(); ++it)
        {
          term=1.0;
          //print(*it);
          // loop over one subset
          for(list<int>::iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
            {
              term *= er[*it2]; 
              //cout << "er[" << *it2 << "]=" << er[*it2] << "\n";
            }
          if (NDEG % 2==1)
            {
              if (i%2==0)
                segno=-1;
              else
                segno=1;
            }
          else
            {
              if (i%2==0)
                segno=1;
              else
                segno=-1;
            }
          cc[i] += term*segno;
        }
      lt.clear();
    }
  for (i=0; i < NDEG+1; i++)
    {
      c[i] = real(cc[i]);
    }
}
#endif
#if 0
void print_backward_err(char *str, vldbl c[], cmplx cr[])
{
  /* we follow FLocke here */
  int k1;
  vldbl relerr, relerrmax, cc[NDEG+1];
  calc_coeff(cc,  cr);
  for (k1=0; k1 < NDEG+1; k1++)
    {
      relerr=abs((c[k1]==vldbl(0.0))?(cc[k1] - c[k1]):(cc[k1] - c[k1])/c[k1]); 
      cout << "k1=" << k1 << " relerr" << relerr << " " << c[k1] <<" " << cc[k1] << "\n" ;
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((c[k1]==0)?(cc[k1] - c[k1]):(cc[k1] - c[k1])/c[k1]); 
        }
    }
  cout << "[" << str << "relative accuracy=" << relerrmax << "\n";
}
#endif
//#define STATIC
int main(int argc, char *argv[])
{
#ifdef STATIC
  pvector<pcmplx,NDEG> roots;
#else
  pvector<pcmplx> roots(NDEG);
#endif
  char testo2[256];
  cmplx cr[NDEG];
#ifdef STATIC
  pvector<pdbl,NDEG+1> c(NDEG+1);
#else
#ifdef CPOLY
  pvector<pcmplx> c(NDEG+1);
#else
  pvector<pdbl> c(NDEG+1);
#endif
#endif 
  int algo, i;
  numty ca[NDEG+1];
  numty allrelerr[NDEG];
  srand48(0);
  if (argc == 2)
    {
      algo = atoi(argv[1]);
    }
  else
    {
      algo = 1;
    }
  if (algo < 0 || algo > 2)
    {
      printf("algorithm must be between 1 and XX algo=%d\n", algo);
      exit(-1);
    }
  calc_coeff_dep_on_case(ca, er);
#ifdef CPOLY
  for (i=0; i < NDEG+1; i++)
    c[i]=pcmplx(pdbl(ca[i]),0.0);
#else
  for (i=0; i < NDEG+1; i++)
    c[i]=pdbl(ca[i]);
#endif
  //c.show("boh");
  //cout << "coeff=" << c[NDEG] << "\n";
  if (algo==0)
    {
#ifdef STATIC
      rpoly<pdbl,NDEG,false, pcmplx> rp;
#else
#ifdef CPOLY
      cpoly<pcmplx,-1,pdbl> rp(NDEG);
#else
      rpoly<pdbl,-1,false,pcmplx> rp(NDEG);
#endif
#endif
      rp.set_coeff(c);
      rp.show();
      //rp.zroots(roots,false);
      rp.find_roots(roots,false);
      //roots.show();  
      sprintf(testo2, "OPS");
    }
  else if (algo==1)
    {
#ifdef STATIC
      rpoly<pdbl,NDEG,true,pcmplx> rphqr;
#else
#ifdef CPOLY
      cpoly<pcmplx,-1,pdbl> rphqr(NDEG);
#else
      rpoly<pdbl,-1,true,pcmplx> rphqr;
#endif
#endif
      rphqr.set_coeff(c);
      //rphqr.show("p(x)=");
      rphqr.find_roots(roots);
      sprintf(testo2, "HQR");
      //roots.show();  
    }
#ifdef AURENTZ
  else if (algo==2)
    {
      pvector<double,NDEG+1> ca;
      pvector<complex<double>,NDEG> rootsa;
      double alpha;
      for (i=0; i < NDEG+1; i++)
       ca[i]=c[i]; 
      if (balance==1)
        wrap_balance(NDEG,ca, &alpha);
      else
        alpha=1.0;
      wrap_damvw(NDEG,ca, rootsa, balance);
      for (i=0; i < NDEG; i++)
        roots[i]=alpha*rootsa[i]; 
      sprintf(testo2, "AUR");
      //roots.show();  
    }
#endif
  for (i=0; i < NDEG; i++)
    cr[i] = cmplx(roots[i]);
  // sort roots and calculate relative error
#if 0
  if (allreal==true)
    {
      struct rerot { vldbl re; vldbl im;};
      std::array<struct rerot,NDEG> rero;
      for (i=0; i < NDEG; i++)
        {
          rero[i].re = er[i].real();
          rero[i].im = er[i].imag();
        }
      std::sort(rero.begin(),rero.end(), [&] (struct rerot a, struct rerot b)-> bool {return a.re < b.re;});
      if (doswap==true)
        {
          for (i=0; i+1 < NDEG; i+=2)
            {
              if (rero[i].im > rero[i+1].im)
                swap(rero[i],rero[i+1]);
            }
        }
      for (i=0; i < NDEG; i++)
        {
          er[i] = rero[i].re + 1il*rero[i].im;
          cout << "EX root= " << er[i] << "\n";
        }
      for (i=0; i < NDEG; i++)
        {
          rero[i].re = cr[i].real();
          rero[i].im = cr[i].imag();
        }
      std::sort(rero.begin(),rero.end(), [&] (struct rerot a, struct rerot b)-> bool {return a.re < b.re;});
      if (doswap==true)
        {
          for (i=0; i+1 < NDEG; i+=2)
            {
              if (rero[i].im > rero[i+1].im)
                swap(rero[i],rero[i+1]);
            }
        }
      for (i=0; i < NDEG; i++)
        {
          cr[i] = rero[i].re + 1il*rero[i].im;
          cout << "CALC root= " << cr[i] << "\n";
        }
    }
  else
    {
      sort_sol_opt(cr, er);
    };
#endif
  sort_sol_opt(cr, er, allrelerr);
  print_roots(testo2, er, cr, allrelerr);
  cout << "Forward relarive error:\n";
  print_accuracy_at(testo2, cr, er, allrelerr);
  //cout << "Backward relarive error:\n";
  //print_backward_err(testo2, ca, cr); 
  return 0;
}
