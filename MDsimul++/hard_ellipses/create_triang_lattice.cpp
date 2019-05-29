#include "../mdlib/pvector.H"
#include "./hardell.H"
using ntype = long double;
int main(int argc, char **argv)
{ 
  ntype dx, dy, sig, exfact, a, b;
  pvector<ntype,2> u1, u2, Ls, L;
  vector<hardell<ntype>> he; 
  ntype fact, phim;
  int N=5000;
  //ntype x, y;
  int i=0, ix, iy, maxix, maxiy;
  // assum pars.a > pars.b 

  a=0.55;
  b=0.5;
  he.resize(N);
  maxix = maxiy = ceil(sqrt(((ntype) (N/2)))); 
  if (argc >= 2)
    maxix = atoi(argv[1]);
  if (argc >= 3)
    maxiy = atoi(argv[2]);
  cerr << "maxix=" << maxix << " maxiy=" << maxiy << "\n";
  u1 << 1,0;
  u2 << 0,sqrt(3.0);
  sig = 2.0*b;
  exfact = 2.0*a;
  dx = sig;
  dy = sig;
  phim = M_PI/(sqrt(3.0)*2.0);
  cerr << setprecision(16) << "phi max=" << phim << "\n";
  //Ls[0] = sqrt(N*M_PI*a*b/phim);
  //Ls[1] = sqrt(N*M_PI*a*b/phim);
  //cerr << "Ls = " << Ls[0] << " " << Ls[1] << "\n";
  //x=y=0;
  cerr << "dx=" << dx << " dy=" << dy << "\n";
  cerr << "a=" << a << " b=" << b << "\n";
  //pars.L.show("box");

  pvector <ntype,2> offset;
  offset << 0,0;
  for (auto nl=0; nl < 2; nl++)
    {
      if (nl==1)
       offset << 0.5, sqrt(3.0)/2.0; 
      for (ix = 0; ix < maxix; ix++)
        {
          for (iy = 0; iy < maxiy; iy++)
            {
              if (nl==1 && iy==maxiy-1 && ix==maxix-1)
                break;
              //he[i].r.show("r");
              //cout << "ix=" << ix << " iy=" << iy << "\n";
              he[i].r = offset + (ix*dx)*u1 + (iy*dy)*u2;
              he[i].r[0] *= exfact;
              he[i].na << 1,0;
              he[i].nb << 0,1;
              he[i].a = a;
              he[i].b = b;
              i++;
            }
        }
    }
  N=i;
  L << 0,0;
  for (i=0; i < N; i++)
    {
      if (i==0 || he[i].r[0]+a > L[0])
        L[0]=he[i].r[0]+a;
      if (i==0 || he[i].r[1]+b > L[1])
        L[1]=he[i].r[1]+b;
    }
  cerr<< "L= " << L[0] << " " << L[1] << "\n";
  auto phi=N*M_PI*a*b/(L[0]*L[1]);

  fact=1.0;
  cerr << "phi=" << phi << " N=" << N << "\n";
  for (i=0; i < N; i++)
    { 
      he[i].r -= 0.5L*L;
      he[i].r *=fact;
      //he[i].show("mah");
    }
  L *= fact;
  cout << L[0] << " " << L[1] << "\n";
  for (i=0; i < N; i++)
    { 
      cout << he[i].r[0]  << " " << he[i].r[1] << " " << he[i].na[0] << " " <<
        he[i].na[1] << "\n";
    } 
  return 0;
}
