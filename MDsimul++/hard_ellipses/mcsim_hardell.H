#include "../mdlib/pvector.H"
#include<string>
#include<iostream>
#include<fstream>
#include "./hardell.H"
#include "../mdlib/boxes.H"
#include "../mdlib/linked_cell_lists.H"
#include <time.h>
#include "./params.H"
template <class ntype>
class mcsim {
  ofstream snapfile, restorefile;
  ifstream inicfile;
  vector<rectangle<ntype>> rect; 
  vector<hardell<ntype>> he;
  linked_cell_lists<ntype,hardell<ntype>> heLL;
  void init_LL(void)
    {
      heLL.init(&he,pars.L);
    }
  void build_LL(void)
    {
      heLL.build();
    }
public:
  const char *parsfile="hepars.asc";
  params<ntype> pars;
  mcsim()
    {
      read_pars(parsfile);
      init();
    }
  mcsim(const char* fn)
    {
      read_pars(fn);
      init();
    }

  void init(void)
    {
      rect.resize(pars.N);
      he.resize(pars.N); 
      init_LL();
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
    }
  void run(void)
    {
      long long int tt;
      for (tt=0; tt < pars.steps; tt++)
        {
          //cout << "nothing to do\n";
        } 
    }
  void restart(void)
    {
      build_LL();
    }
  void readconf(string& icn)
    {
      inicfile.open(icn.c_str());
      // read ascii iniconf here 
      inicfile.close();
    }
  void update_all_rect(void)
    {
      for (auto i=0; i < pars.parnum; i++)
        {
          update_rect(i);
        }
    }
  void update_rect(int i)
    {
      rect[i].r = he[i].r;
      rect[i].R.set_row(0,he[i].na);
      rect[i].R.set_row(1,he[i].nb);
    }
  void createconf(int type)
    { // type = 0 SC  1 = TRIANGULAR
      switch (type)
        {
        case 0:
            {
              int i=0, ix, iy, maxix, maxiy;
              maxix = maxiy = ceil(sqrt(((ntype) pars.N))); 
              for (ix = 0; ix < maxix; ix++)
                for (iy = 0; iy < maxiy; iy++)
                  {
                    he[i].r << 0,0;
                    he[i].na << 1,0;
                    he[i].nb << 0,1;
                    update_rect(i);
                    i++;
                    if (i==pars.N)
                      break;
                  }
              break;
            }
        case 1:
            {
              // triangular to be implemented
              break;
            }
        default:
            {

            }
        }
      build_LL();
    }
  void read_pars(const char *fn)
    {
      if (fn==nullptr)
        pars.read(parsfile);
      else
        pars.read(fn);
    }
};
