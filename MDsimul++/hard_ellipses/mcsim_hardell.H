#include "../mdlib/pvector.H"
#include<string>
#include<iostream>
#include<fstream>
#include "./hardell.H"
#include "../mdlib/boxes.H"
template<class ntype>
class params
{
  ifstream parfile;
  string filename;
  void fopen(const char *fn)
    {
      filename.assign(fn);
      parfile.open(filename.c_str());
    }
  void fclose(void)
    {
      parfile.close();
    }
public:
  long long int steps;
  int simtype;
  int seed;// -1 means random seed 
  ntype a, b;
  int N; // number of particles
  pvector<ntype,2> L; // box
  pvector<ntype,2> M; // linked list cells 
  params()
    {
      N=1000;
      L[0]=L[1]=L[2]=10.0;
      simtype = 0; // 0=NTV; 1=NPT; 2=gran canonico
      seed=0;
    }

  void read(const char *fn)
    {
      // read simulation parameters first
      fopen(fn);
      parfile >> simtype;
      parfile >> steps;
      parfile >> seed;
      parfile >> a;// a semi-axis
      parfile >> b;// b semi-axis
      parfile >> N;
      parfile >> L[0];
      parfile >> L[1];
      fclose();
    }
};

template <class ntype>
class mcsim {
  ofstream snapfile, restorefile;
  vector<rectangle<double>> rect; 
  vector<hardell<double>> he;
public:
  params<ntype> pars;
  void init(void)
    {
      rect.resize(pars.N);
      he.resize(pars.N); 

      if (pars.seed > 0)
        srand48(pars.seed);
      else
        srand48(time(NULL));

      for(auto &value: he) 
        {
          value.a = pars.a;
          value.b = pars.b;
        }
      for(auto &value: rect) 
        {
          value.sax[0] = pars.a;
          value.sax[1] = pars.b;
        } 
#if 0 // move elsewhere
      rect.r = he.r;
      rectA.R.set_row(0,A.na);
      rectA.R.set_row(1,A.nb);
#endif
 
    }
  void run(void)
    {


    }
  void restart(void)
    {

    }
  void createconf(int type)
    { // type = 0 SC  1 = TRIANGULAR
      switch (type)
        {
        case 0:
          break;
        default:
          break;
        }
    }
  void read_pars(const char *fn)
    {
      pars.read(fn);
    }

};
