#include "../mdlib/pvector.H"
#include<vector>
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
  ntype zeta;
  int susnmin, susnmax, npav, nexc;
  ntype maxmsd;
  char fnt[256];
  ntype targetPhi;
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
      zeta=0;
      susnmin=1000;
      susnmax=1001;
      nexc=10;
      npav=1000;
      curStep=0;
      targetPhi=-1; //<= 0 means disabled
      maxmsd = -1;// equilibration
    }
  void write(fstream& o)
    {
      o << simtype << " ";
      o << iniconf << " ";
      o << steps << " ";
      o << curStep << " ";
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
      o << resetacceptVol << " ";
      o << targetPhi << " ";
      o << maxmsd << "\n";
      o << zeta << " ";
      o << susnmin << " ";
      o << susnmax << " ";
      o << npav << " ";
      o << nexc << "\n";
    }
  void read(fstream& i)
    {
      i >> simtype;
      i >> iniconf;
      i >> steps;
      i >> curStep;
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
      i >> targetPhi;
      i >> maxmsd;
      i >> zeta, susnmin, susnmax, npav, nexc;
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
  long long int volmoves=0, volrej=0, excmove=0, excrej=0;
  //long long int curStep;
  unsigned int which=0;
  vector<long long int> sushisto;
  params<ntype> *ppars;
  //long long int resetaccept=0, resetacceptVol=0;

  //ntype deltra, delrot, vmax;
  void write(fstream& o)
    {
      // write full state to ascii file
      o << trarej << " " << rotrej << " " << totrej << "\n";
      o << totmoves << " " << tramoves << " " << rotmoves << "\n";
      o << volmoves << " " << volrej << "\n";
      //o << curStep << " ";
      o << which << "\n";
      o << excmove << " " << excrej << "\n";
      for (auto &e: sushisto)
        o << e << " ";
      o << "\n";
    }

  void read(fstream& i)
    {
      init_vecs();
      // read full state from ascii file
      i >> trarej;
      i >> rotrej;
      i >> totrej;
      i >> totmoves;
      i >> tramoves;
      i >> rotmoves;
      i >> volmoves;
      i >> volrej;
      //i >> curStep;
      i >> which;
      i >> excmove >> excrej;
      for (auto &e: sushisto)
        i >> e;
    }
  void init_vecs(void)
    {
      sushisto.resize(ppars->susnmax-ppars->susnmin+1);
      for (auto &e: sushisto)
        e = 0;
    }
  simstate ()
    {
    }
};


