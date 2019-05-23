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
  ntype ranf(void)
    {
      return drand48();
    }
  ofstream snapfile, restorefile;
  ifstream inicfile;
  vector<rectangle<ntype>> rect; 
  vector<hardell<ntype>> he;
  linked_cell_lists_2d<hardell<ntype>,rectangle<ntype>> heLL;
  pvector<ntype,2> rold;
  pvector<ntype,2> naold, nbold;
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
  void init_LL(void)
    {
      heLL.init(&he,&rect,pars.L,pars.rcut);
    }
  void build_LL(void)
    {
      heLL.build();
    }
  void read_pars(const char *fn)
    {
      if (fn==nullptr)
        pars.read(parsfile);
      else
        pars.read(fn);
    }
  void pbc(int ip)
    {
      pvector<ntype,2> L2, delr;
      L2 = 0.5*pars.L;
      delr << pars.L[0]*rint(he[ip].r[0]/pars.L[0]), pars.L[1]*rint(he[ip].r[1]/pars.L[1]);
      he[ip].r = he[ip].r - delr;
    }
  int random_move(int ip) // return move type 0=tra 1=rot
    {
      ntype p;
      p=ranf();
      if (p <= 0.5)
        {
          he[ip].tra_move(pars.deltra);
          update_rect(ip);
          return 0;
        }
      else
        {
          he[ip].rot_move(pars.delrot);
          update_rect(ip);
          return 1;
        } 
    }
  void store_coord(int ip)
    {
      rold = he[ip].r;
      naold = he[ip].na;
      nbold = he[ip].nb;
    }
  void restore_coord(int ip)
    {
      he[ip].r = rold;
      he[ip].na = naold;
      he[ip].nb = nbold;
    }
  void mcmotion(int ip)
    {
      store_coord(ip);
      random_move(ip);
      // check overlap by using LL here!
      pbc(ip);
      heLL.update_LL(ip);
      if (heLL.check_overlap(ip)==true)
        {
          restore_coord(ip);
          heLL.update_LL(ip);
          // reject move
        }
    }
  void move_box(void)
    {
      //NPT volume move
      int i, ii;
      pvector<ntype,2> L2;
      ntype nn;
      ntype vo, lnvn, vn, Lfact, enn, eno, arg, delv=0.0;
      vo = pars.L[0]*pars.L[1];
      lnvn = log(vo) + (ranf()-0.5)*pars.vmax;
      vn = exp(lnvn);
      //cout << "vo=" << vo << " vn=" << vn << " vmax=" << state.vmax << "\n";
      Lfact = pow(vn/vo,1.0/2.0);
      //cout << "Lfact=" << Lfact << "\n";
      pars.L *= Lfact;
      L2 = 0.5*pars.L;
      //cout << "Lnew=" << pars.L[0] << " " << pars.L[1] << "\n";
      for (i=0; i < pars.parnum; i++)
        {
          he[i].r *= Lfact;
          update_rect(i);
          pbc(i);
        }
      //cout << "qui1\n";
      heLL.update_numcells(pars.L);
      //cout << "qui2\n";
#ifdef MC_STORELL
      store_ll_mc();
#endif
      heLL.build();
      for (i=0; i < pars.parnum; i++)
        {
          //cout << "i=" << i << "\n"; 
          //he[i].r.show("r");
          if (heLL.check_overlap(i)==true)
            {
              /* move rejected restore old positions */
              pars.L /= Lfact;
              L2 = 0.5*pars.L;
              for (ii=0; ii < pars.parnum; ii++)
                {
                  he[ii].r /= Lfact;
                  update_rect(ii);
                  pbc(ii);
                }
              //volrejMC++;
              heLL.update_numcells(pars.L);
#ifdef MC_STORELL
              restore_ll_mc();
#else
              heLL.build();
#endif
              return;
            }
        }
      arg = -(1.0/pars.T)*(pars.P*(vn-vo)-(pars.parnum+1)*log(vn/vo)*pars.T);
      if (ranf() > exp(arg))
        {
          /* move rejected restore old positions */
          pars.L /= Lfact;
          L2 = 0.5*pars.L;

          for (i=0; i < pars.parnum; i++)
            {
              he[i].r /= Lfact;
              pbc(i);
            }
          //volrejMC++;
          heLL.update_numcells(pars.L);
#ifdef MC_STORELL
          restore_ll_mc();
#else
          build_LL();
#endif
        }
    }
  void save_snap(int tt)
    {
      if (tt%5000==0)
        {

        }
    }
  void save_restore(int tt)
    {
      if (tt%5000==0)
        {

        }
    }
  void init_measures(void)
    {
      ofstream of;
      of.open("v.dat", ios::trunc);
      of.close();
    }
  void save_measures(int tt)
    {
      ofstream of;
      of.open("v.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      if (tt!=0 && tt%5000==0)
        {
          of << tt << " " << pars.L[0]*pars.L[1] << "\n";
        }
      of.close();
    }
public:
  const char *parsfile="hepars.asc";
  params<ntype> pars;
  simstate<ntype> state;
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
      pars.rcut = 2.0*max(pars.a,pars.b);
      init_LL();
      cout << "deltra=" << pars.deltra << " delrot=" << pars.delrot << " vmax=" << pars.vmax << "\n";
      cout << "P= " << pars.P << " T=" << pars.T << "\n";
      init_measures();
    }
  void run(void)
    {
      long long int tt;
      int ntot, ran, i;
      for (tt=0; tt < pars.steps; tt++)
        {
          ntot=pars.parnum;
          if (pars.simtype==1) //NPT
            ntot+=1;
          heLL.set_ini_numcells();
          for (i=0; i < ntot; i++)
            {
              //cout << "step #" << tt << " i=" << i << "\n";
              ran = (int) (ntot*ranf());
              if (ran < pars.parnum)
                mcmotion(ran);
              else 
                {
                  move_box();
                }
            } 
          // save snapshots and backup file 
          save_snap(tt);
          save_restore(tt);
          save_measures(tt);
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
 
  void createconf(int type)
    { // type = 0 SC  1 = TRIANGULAR
      switch (type)
        {
        case 0:
            {
              ntype dx, dy, x, y;
              int i=0, ix, iy, maxix, maxiy;
              maxix = maxiy = ceil(sqrt(((ntype) pars.N)))+1; 
              dx = pars.L[0]/((ntype)maxix);
              dy = pars.L[1]/((ntype)maxiy);
              x=y=0;
              cout << "dx=" << dx << " dy=" << dy << "\n";
              for (ix = 0; ix < maxix; ix++)
                {
                  for (iy = 0; iy < maxiy; iy++)
                    {
                      he[i].r << x-pars.L[0]*0.5,y-pars.L[1]*0.5;
                      //he[i].r.show("r");
                      he[i].na << 1,0;
                      he[i].nb << 0,1;
                      update_rect(i);
                      i++;
                      //cout << "i=" << i << "pars.N="<< pars.N <<"\n";
                      if (i==pars.N)
                        {
                          build_LL();
                          return;
                        }
                      y += dy;
                    }
                  y=0;
                  x += dx;
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
    }
  
};
