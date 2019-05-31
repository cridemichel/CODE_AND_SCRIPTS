#include "./mcsim_hardell.H"
using ntype = long double;
int main(int argc, char **argv)
{ 
  double fx=1.0, fy=1.0;
  mcsim<long double> heMC;
  bool scalax=false;
  if (argc==1)
    {
      cout << "expand <factorx> <factory> <scale-semiaxes>\n";
      exit(1);
    }
  if (argc >= 2)
    fx = atof(argv[1]);
  if (argc >= 3)
    fy = atof(argv[2]);
  if (argc >= 4)
    scalax=(atoi(argv[3])==1)?true:false;

  heMC.expand(fx, fy, true, scalax);
  cout << "Saved conf: cnf-expanded\n";
  return 0;
}
