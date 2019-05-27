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
  long long int savemeas, savesnap, saverest;
  long long int curStep;
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
      saverest=0;
      targetAcceptVol=0.5;
      targetAccept=0.5;
      adjsteps = 0;
      resetaccept = 10000;
      resetacceptVol = 10000;
      curStep=0;
   }
  void write(fstream& o)
    {
      o << simtype << " ";
      o << iniconf << " ";
      o << steps << " ";
      o << seed << " ";
      o << a << " ";// a semi-axis
      o << b << " ";// b semi-axis
      o << N << " ";
      o << L[0] << " ";
      o << L[1] << "\n";
      o << P << " ";
      o << T << " ";
      o << deltra << " ";
      o << delrot << " ";
      o << vmax << "\n";
      o << savemeas << " ";
      o << savesnap << " ";
      o << saverest << "\n";
      o << adjsteps << " ";
      o << outsteps << "\n";
      o << resetaccept << " ";
      o << resetacceptVol << "\n";
    }
  void read(fstream& i)
    {
      i >> simtype;
      i >> iniconf;
      i >> steps;
      i >> seed;
      i >> a;// a semi-axis
      i >> b;// b semi-axis
      i >> N;
      i >> L[0];
      i >> L[1];
      i >> P;
      i >> T;
      i >> deltra;
      i >> delrot;
      i >> vmax;
      i >> savemeas;
      i >> savesnap;
      i >> saverest;
      i >> adjsteps;
      i >> outsteps;
      i >> resetaccept;
      i >> resetacceptVol;
    }
  void readpars(const char *fn)
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
      read(pft); 
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
  long long int curStep;
  //long long int resetaccept=0, resetacceptVol=0;

  //ntype deltra, delrot, vmax;
  void write(fstream& o)
    {
      // write full state to ascii file
      o << trarej << " " << rotrej << " " << totrej << "\n";
      o << totmoves << " " << tramoves << " " << rotmoves << "\n";
      o << volmoves << " " << volrej << "\n";
      o << curStep << " " << "\n";
    }

  void read(fstream& i)
    {
      // read full state from ascii file
      i >> trarej;
      i >> rotrej;
      i >> totrej;
      i >> totmoves;
      i >> tramoves;
      i >> rotmoves;
      i >> volmoves;
      i >> volrej;
      i >> curStep;
    }
  simstate ()
    {
    }
};


