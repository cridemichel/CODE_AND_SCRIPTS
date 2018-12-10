#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
int main()
{
  using numty=double;
  rpoly<numty,4> rp;
  rpoly<numty,5> qp;
  pvector<numty, 6> cq;
  pvector<numty, 5> c;
  complex<numty> x1c, x2c, x3c, x4c;
  pvector<complex<numty>, 4> roots;
  pvector<complex<numty>, 5> qroots;

#if 1
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
  c[0] = 274;
  c[1] = -450;
  c[2] = 255;
  c[3] = -60;
  c[4] = 5;
  rp.set_coeff(c);
  rp.find_roots(roots);
  roots.show("boh"); 
  x1c = 1;
  x2c = 2;
  x3c = 3;
  x4c = 4;
  cout << "ntype size=" << sizeof(numty) << "\n";
  c[4] = 1.0;
  c[3] = (-(x1c+x2c+x3c+x4c)).real();
  c[2] = (x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c).real(); 
  c[1] = (-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c)).real();
  c[0] = (x1c*x2c*x3c*x4c).real();
  rp.set_coeff(c);
  rp.show();
  //cout << "p(x)=" << rp.evaldpoly(1.2) << "\n";
  exit(-1);
  for (int i=0; i < 100000000; i++)
    {
      c[4]=1.0;
      c[3]=drand48()-0.5;
      c[2]=drand48()-0.5;
      c[1]=drand48()-0.5;
      c[0]=drand48()-0.5;
      rp.set_coeff(c);
      rp.find_roots(roots);
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
