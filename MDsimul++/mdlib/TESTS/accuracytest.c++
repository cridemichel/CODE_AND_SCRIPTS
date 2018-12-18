#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
#include<list>
#include<string>
bool allreal=false, doswap=false;
#ifndef CASO
#define CASO 4
#endif
#undef M_PI
#define M_PI 3.1415926535897932384626433832795029L
using numty = double;
void calc_coeff(long double* c, complex<long double> er[]);
void calc_coeff_dep_on_case(long double* c, complex<long double> **r)
{
  int i;
#if CASO==1
// wilkinson
//
#define NDEG 10
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==2
// wilkinson
//
#define NDEG 15
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==3
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==4
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  er[0] = -2.1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]+0.2L;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==5
#define NDEG 10
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = 1.0L/i;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==6
#define NDEG 15
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = 1.0L/i;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==7
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = 1.0L/i;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==8
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = 1.0L/pow(2,NDEG/2-i);
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==9
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = 1.0L/pow(2,NDEG/2-i)-3.0L;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==10
#define NDEG 20
#if 1
  static complex <long double> er[NDEG];

  for (i=0; i < NDEG; i++)
    { 
      er[i] = cos(M_PI*(2.0L*(i+1L)-1L)/40.0L);
    }
#endif
  calc_coeff(c, er);
  allreal=true;
#if 0
  static complex <long double> er[NDEG]=
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
  *r = er;
#elif CASO==11
#define NDEG 20
  static complex <long double> er[NDEG]={-0.9888308262251285450697429 - 
  0.1490422661761744469293547*1i, -0.9888308262251285450697429 + 
  0.1490422661761744469293547*1i, -0.9009688679024191262361023 - 
  0.4338837391175581204757683*1i, -0.9009688679024191262361023 + 
  0.4338837391175581204757683*1i, -0.7330518718298263285224315 - 
  0.6801727377709193901873587*1i, -0.7330518718298263285224315 + 
  0.6801727377709193901873587*1i, -0.5000000000000000000000000 - 
  0.8660254037844386467637232*1i, -0.5000000000000000000000000 + 
  0.8660254037844386467637232*1i, -0.2225209339563144042889026 - 
  0.9749279121818236070181317*1i, -0.2225209339563144042889026 + 
  0.9749279121818236070181317*1i, 
 0.0747300935864242542909397 - 0.9972037971811801482250299*1i, 
 0.0747300935864242542909397 + 0.9972037971811801482250299*1i, 
 0.3653410243663950145447380 - 0.9308737486442042556377992*1i, 
 0.3653410243663950145447380 + 0.9308737486442042556377992*1i, 
 0.6234898018587335305250049 - 0.7818314824680298087084445*1i, 
 0.6234898018587335305250049 + 0.7818314824680298087084445*1i, 
 0.8262387743159948719451626 - 0.5633200580636220277492615*1i, 
 0.8262387743159948719451626 + 0.5633200580636220277492615*1i, 
 0.9555728057861407328113341 - 0.2947551744109042168307730*1i, 
 0.9555728057861407328113341 + 0.2947551744109042168307730*1i};
  doswap=true;
  allreal=true;
  calc_coeff(c, er);
  *r = er;
#elif CASO==30
#define NDEG 20
  allreal=true;
  er[0]=0.1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/(10.0L);
    }
  calc_coeff(c, er);
  *r = er;
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
  r = NULL;
#elif CASO==32
#define NDEG 10
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[NDEG] = 1.0;
  c[0]=1.0;
  c[1] = -300.0;
  c[2] = 30000.0;
  c[3] = -1E6;
  r=NULL;
#elif CASO==5
#define NDEG 11
  static numty ct ={-1E22, 2E21, -1E20, 0, 0, 0, 0, 0, 0, 0, 1.0};
  for (i=0; i < 11; i++)
   c[i] = ct[i]; 
  r=NULL;
#elif CASO==33
#define NDEG 20
  static complex <long double> er[NDEG];
  er[0]=1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/10.0L;
    }
  calc_coeff(c, er);
  *r = er;
#endif
} 
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
void sort_sol_opt(complex<long double> sol[NDEG], complex<long double>* exsol)
{
  int k1, k2, kk;
  double v, vmin;
  complex<numty> solt[NDEG];
  int perm[NDEG];
  int perm_min[NDEG];
  bool ini=true;
  for (k1=0; k1 < NDEG; k1++)
    {
      perm[k1] = k1;
    }
  do {
    //  std::cout << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';
    //cout << perm << "\n"; 
    v = 0;
    for (k2=0; k2 < NDEG; k2++)
      {
        v += (exsol[k2]==complex<long double>(0,0))?abs(sol[perm[k2]]-exsol[k2]):abs((sol[perm[k2]]-exsol[k2])/exsol[k2]);
      }
    if (ini==true || v < vmin)
      {
        ini=false;
        for (kk=0; kk < NDEG; kk++)
          perm_min[kk]=perm[kk];
        vmin = v;
      }
  } 
  while (std::next_permutation(perm,perm+NDEG));

  for (k2=0; k2 < NDEG; k2++)
    solt[k2] = sol[k2];

  for (k2=0; k2 < NDEG; k2++)
    sol[k2] = solt[perm_min[k2]];
}
numty print_accuracy_at(char *str, complex<long double>* csol, complex<long double> *exsol)
{
  /* we follow FLocke here */
  int k1;
  numty relerr, relerrmax;
  for (k1=0; k1 < NDEG; k1++)
    {
      relerr=abs((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16G\n", str, relerrmax);
  return relerrmax;
}

void print_roots(char *str, complex<long double> *er)
{
  printf("CASE %s\n", str);
  for (auto i=0; i < NDEG; i++)
    cout << "root #" << i << " "<< er[i] << setprecision(16) << "\n";
}

void print( list<int> l){
    for(list<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}

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

void calc_coeff(long double* c, complex<long double> er[NDEG])
{
  int i, j;
  std::list<std::list<int>> subsets;// list containing all subsets of a given length
  complex<long double> cc[NDEG+1];
  complex<long double> term, segno;
  // list with all list of subsets
  list<list<int>> lt;   
  // a subset of integers
  list<int> l;
  int array[NDEG];
  for (j=0; j < NDEG; j++)
    array[j] = j;
  cc[NDEG] = 1.0;
  for (i=0; i < NDEG; i++)
    {
      // all subsets of NDEG-i elements are stored in list lt
      subset(lt,array,NDEG,NDEG-i,0,l);
      cc[i] = 0;
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
void print_backward_err(char *str, long double* c, complex<long double> *cr)
{
  /* we follow FLocke here */
  int k1;
  long double relerr, relerrmax, cc[NDEG+1];
  calc_coeff(cc,  cr);
  for (k1=0; k1 < NDEG+1; k1++)
    {
      relerr=abs((c[k1]==0)?(cc[k1] - c[k1]):(cc[k1] - c[k1])/c[k1]); 
      cout << "k1=" << k1 << " relerr" << relerr << " " << c[k1] <<" " << cc[k1] << "\n" ;
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((c[k1]==0)?(cc[k1] - c[k1]):(cc[k1] - c[k1])/c[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
}
int main(int argc, char *argv[])
{
  pvector<complex<numty>,NDEG> roots;
  char testo2[256];
  complex<long double> cr[NDEG], *er;
  pvector<numty,NDEG+1> c;
  int algo, i;
  long double ca[NDEG+1];
  if (argc == 2)
    {
      algo = atoi(argv[1]);
    }
  else
    {
      algo = 1;
    }
  if (algo < 0 || algo > 1)
    {
      printf("algorithm must be between 1 and XX algo=%d\n", algo);
      exit(-1);
    }
  calc_coeff_dep_on_case(ca, &er);
  char testo[] = "Case CASO";
  if (er!=nullptr)
    print_roots(testo, er); 
  for (i=0; i < NDEG+1; i++)
    c[i]=ca[i];
  if (algo==0)
    {
      rpoly<numty,NDEG> rp;
      rp.set_coeff(c);
      rp.show();
      rp.find_roots(roots);
      roots.show();  
      sprintf(testo2, "OPS");
    }
  else if (algo==1)
    {
      rpoly<numty,NDEG,true> rphqr;
      rphqr.set_coeff(c);
      rphqr.show("p(x)=");
      rphqr.find_roots(roots);
      sprintf(testo2, "HQR");
      //roots.show();  
    }
  for (i=0; i < NDEG; i++)
    cr[i] = roots[i];
  // sort roots and calculate relative error
  if (allreal==true)
    {
      struct rerot { long double re; long double im;};
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
  cout << "Forward relarive error:\n";
  print_accuracy_at(testo2, cr, er);
  //cout << "Backward relarive error:\n";
  //print_backward_err(testo2, ca, cr); 
  return 0;
}
