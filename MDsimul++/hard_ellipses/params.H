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
  long long int savemeas, savesnap;
  int curStep;
  ntype targetAccept, targetAcceptVol;
  long long int adjsteps, outsteps, resetaccept, resetacceptVol; 
  char fnt[256];
  params()
    {
      N=1000;
      L[0]=L[1]=10.0;
      simtype = 0; // 0=NTV; 1=NPT; 2=gran canonico
      seed=0;
      T=1.0;
      P=10.0;
      vmax=0.001;
      delrot = 0.3;
      deltra = 0.1;
      savemeas=50000;
      savesnap=0;
      targetAcceptVol=0.5;
      targetAccept=0.5;
      adjsteps = 0;
      resetaccept = 10000;
      resetacceptVol = 10000;
   }
  void read(const char *fn)
    {
      std::size_t found;
      fstream pft;
      string line;
      fopen(fn);
      sprintf(fnt, "%s.tmp", fn);
      pft.open(fnt, ios::trunc|ios::out|ios::in);
      if (pft.is_open())
        {
          while ( getline (parfile,line) )
            {
              found=line.find_first_of("#");
              if (found==std::string::npos)
                {
                  pft << line << "\n";
                }
            }
        }
      else
        {
          cout << "[ERROR] Error in opening file " << fnt << "\n";
          exit(0);
        }
      fclose();
      pft.seekg(ios::beg);
      pft >> simtype;
      pft >> iniconf;
      pft >> steps;
      pft >> seed;
      pft >> a;// a semi-axis
      pft >> b;// b semi-axis
      pft >> N;
      pft >> L[0];
      pft >> L[1];
      pft >> P;
      pft >> T;
      pft >> deltra;
      pft >> delrot;
      pft >> vmax;
      pft >> savemeas;
      pft >> savesnap;
      pft >> adjsteps;
      pft >> outsteps;
      pft >> resetaccept;
      pft >> resetacceptVol;
      pft.close();
      remove(fnt);// remove temporary file
    }
};

template <class ntype> 
class simstate 
{
public:
  long long int trarej=0, rotrej=0, totrej=0, totmoves=0, tramoves=0, rotmoves=0;
  long long int volmoves=0, volrej=0;
  long long int curStep, adjsteps;
  //long long int resetaccept=0, resetacceptVol=0;

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


