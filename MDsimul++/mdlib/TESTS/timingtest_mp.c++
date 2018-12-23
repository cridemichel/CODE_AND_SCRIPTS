#include <stdlib.h>
#include <stdio.h>
#include "pmatrix_mp.H"
#include "./rpoly_mp.H"
//#include<complex>
#ifndef NDEG
#define NDEG 6
#endif
//#define STATIC
// fino a N=40 la versione statica è più veloce
using numty=double;
int main(int argc, char* argv[])
{
  int caso, j, maxiter;
#ifdef STATIC
  rpoly<numty,NDEG> rp;
  rpoly<numty,NDEG,true> rphqr;
  pvector<numty,NDEG+1> c;
  pvector<cmplx,NDEG> roots;
#else
  rpoly<numty,-1> rp(NDEG);
  rpoly<numty,-1,true> rphqr(NDEG);
  pvector<numty,-1> c(NDEG+1);
  pvector<cmplx,-1> roots(NDEG);
#endif
  //cmplx x1c, x2c, x3c, x4c, x5c, x6c;
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
  srand48(4242);
  numty sig=1E10;
  if (argc>=2)
    caso = atoi(argv[1]);
  else
    caso = 0;
  if (caso < 0 || caso > 1)
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
  for (int i=0; i < maxiter; i++)
    {
      for (j=0; j < NDEG+1; j++)
        c[j]=sig*(drand48()-0.5);
      if (caso==0)
        {
          rp.set_coeff(c);
          //rp.zroots(roots, true);
          rp.find_roots(roots, false);
        }
      else
        {
          rphqr.set_coeff(c);
          rphqr.find_roots(roots);
        }
    }
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
