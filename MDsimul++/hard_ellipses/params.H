#include "../mdlib/pvector.H"
#include<string>
#include<iostream>
#include<fstream>
using namespace std;
template<class ntype>
class params
{
  ifstream parfile;
  string filename;
  void fopen(const char *fn)
    {
      filename.assign(fn);
      parfile.open(filename.c_str());
      if (parfile.fail())
        {
          cout << "Parameters file named '" << fn << "' not found, exiting...\n";
          exit(0);
        }
    }
  void fclose(void)
    {
      parfile.close();
    }
public:
  string iniconf; // configurazione iniziale da leggere
  long long int steps;
  int simtype;
  int seed;// -1 means random seed 
  ntype a, b;
  int parnum; // number of particles
  int& N=parnum;// reference declarations are aliases!
  pvector<ntype,2> L; // box
  pvector<ntype,2> M; // linked list cells 
  ntype rcut;
  ntype P, T, deltra, delrot, vmax;
  params()
    {
      N=1000;
      L[0]=L[1]=10.0;
      simtype = 0; // 0=NTV; 1=NPT; 2=gran canonico
      seed=0;
      T=1.0;
      P=10.0;
      vmax=0.00001;
      delrot = 0.3;
      deltra = 0.1;
  }

  void read(const char *fn)
    {
      // read simulation parameters first
      fopen(fn);
      parfile >> simtype;
      parfile >> iniconf;
      parfile >> steps;
      parfile >> seed;
      parfile >> a;// a semi-axis
      parfile >> b;// b semi-axis
      parfile >> N;
      parfile >> L[0];
      parfile >> L[1];
      parfile >> P;
      parfile >> T;
      parfile >> deltra;
      parfile >> delrot;
      parfile >> vmax;
      fclose();
    }
  
};

template <class ntype> 
class simstate 
{
public:
  long long int trarej=0, rotrej=0, totrej=0, totmoves=0, tramoves=0, rotmoves=0;
  long long int curstep;
  //ntype deltra, delrot, vmax;
  void write_state(void)
    {
      // write full state to ascii file
    }
  void read_state(void)
    {
      // write full state to ascii file
    }
  simstate ()
    {
    }
};


