#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
int main()
{
  using numty=double;
  rpoly<numty,5> rp;
  rpoly<numty,5,true> rphqr;
  pvector<numty, 6> c;
  complex<numty> x1c, x2c, x3c, x4c, x5c;
  pvector<complex<numty>, 5> roots;
#if 0
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
  x1c = 1;
  x2c = 1;
  x3c = 9.9999;
  x4c = 9.9999;
  x5c = 10.001;
  cout << "ntype size=" << sizeof(numty) << "\n";
  c[5] = 1.0;
  c[4] = (-(x1c+x2c+x3c+x4c+x5c)).real();
  c[3] = (x1c*x2c + x1c*x3c + x2c*x3c + x1c*x4c + x2c*x4c + x3c*x4c + x1c*x5c + 
 x2c*x5c + x3c*x5c + x4c*x5c).real(); 
  c[2] = (-(x1c*x2c*x3c) - x1c*x2c*x4c - x1c*x3c*x4c - x2c*x3c*x4c - x1c*x2c*x5c - x1c*x3c*x5c - x2c*x3c*x5c - 
   x1c*x4c*x5c - x2c*x4c*x5c - x3c*x4c*x5c).real();
  c[1] =(x1c*x2c*x3c*x4c + x1c*x2c*x3c*x5c + x1c*x2c*x4c*x5c + 
   x1c*x3c*x4c*x5c + x2c*x3c*x4c*x5c).real();
  c[0] =(-x1c*x2c*x3c*x4c*x5c).real(); 

  rp.set_coeff(c);
  rp.show();
  rp.find_roots(roots);
      
  roots.show("oqs roots=");
  
  rphqr.set_coeff(c);
  rphqr.show();
  rphqr.find_roots(roots);
  roots.show("hqr roots=");
  //cout << "p(x)=" << rp.evaldpoly(1.2) << "\n";
  exit(-1);
#endif
  srand48(4242);
  double sig=1E30;
  for (int i=0; i < 10000000; i++)
    {
      c[5]=sig*(drand48()-0.5);
      c[4]=sig*(drand48()-0.5);
      c[3]=sig*(drand48()-0.5);
      c[2]=sig*(drand48()-0.5);
      c[1]=sig*(drand48()-0.5);
      c[0]=sig*(drand48()-0.5);
#if 0
      rp.set_coeff(c);
      rp.find_roots(roots);
#else
      rphqr.set_coeff(c);
      rphqr.find_roots(roots);
#endif
    }
#if 0
  printf("root 0: %.15LG %.15LG\n", roots[0].real(), roots[0].imag());
  printf("root 1: %.15LG %.15LG\n", roots[1].real(), roots[1].imag());
  printf("root 2: %.15LG %.15LG\n", roots[2].real(), roots[2].imag());
  printf("root 3: %.15LG %.15LG\n", roots[3].real(), roots[3].imag());
  std::cout.precision(16);
  cout << abs((roots[3]-((complex<numty>)x2c))/((complex<numty>)x2c))<<  setprecision(16) << "\n";
  cout << setprecision(16) << abs((roots[0]-((complex<numty>)x3c))/((complex<numty>)x3c))<<"\n"; 
  cout << setprecision(16) << abs((roots[2]-((complex<numty>)x1c))/((complex<numty>)x1c))<<"\n";
  cout << setprecision(16) << abs((roots[1]-((complex<numty>)x4c))/((complex<numty>)x4c))<<"\n";
  roots.show();
#endif
  return 0;
}
