#include "./mcsim_hardell.H"
using ntype = long double;
int main(int argc, char **argv)
{ 
  double fx=1.0, fy=1.0;
  mcsim<double> heMC;

  if (argc==1)
    {
      cout << "expand <factorx> <factory>\n";
      exit(1);
    }
  if (argc >= 2)
    fx = atof(argv[1]);
  if (argc >= 3)
    fy = atof(argv[2]);

  heMC.expand(fx, fy, true, true);
  cout << "Saved conf: cnf-expanded\n";
  return 0;
}
