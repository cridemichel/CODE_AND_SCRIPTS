#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
//#define KAMENY
//#define MIGNOTTE
#define WILK2
#ifdef KAMENY
#define NDEG 9
#elif defined(MIGNOTTE)
#define NDEG 6
#elif defined(WILKINSON)
#define NDEG 10
#elif defined(WILK2)
#define NDEG 10
#else
#endif
int main()
{
  using numty=double;
  rpoly<numty,NDEG> rp;
  rpoly<numty,NDEG,true> rphqr;
  pvector<numty, NDEG+1> c;
  pvector<complex<numty>, NDEG> roots;
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
   //c << -1685011.4849863185081630945205688,162908947.864249706268310546875,3464224691.560069561004638671875,20202439788.892055511474609375,31022318528.748729705810546875,-352127.44537145795766264200210571,1;
#ifdef KAMENY
//Polynomials with few very clustered roots.
//Kameny  
  numty K = 1E50;
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[0] = 9.0;
  c[2] = -6.0*K*K;
  c[4] = K*K*K*K;
  c[9] = K*K; 
#elif defined(MIGNOTTE)
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[NDEG] = 1.0;
  c[0]=1.0;
  c[1] = -300.0;
  c[2] = 30000.0;
  c[3] = -1E6;
#elif defined(WILKINSON)
c << 3628800, - 10628640, +12753576, - 8409500, 3416930, -902055, 157773, -18150,1320,-55, 1.0;
#elif defined(WILK2)
c << -1E22, 2E21, -1E20, 0, 0, 0, 0, 0, 0, 0, 1.0;
#else

#endif
#endif
  rp.set_coeff(c);
  rp.show();
  rp.find_roots(roots);
  roots.show();  
  exit(-1);
  return 0;
}
