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
  linked_cell_lists_2d<hardell<ntype>,rectangle<ntype>> ellLL;
  pvector<ntype,2> rold;
  pvector<ntype,2> naold, nbold;
  ntype acceptance, traaccept, rotaccept, volaccept;
  void update_all_rect(void)
    {
      for (auto i=0; i < pars.parnum; i++)
        {
          update_rect(i);
        }
    }
  void calc_acceptance_and_adjust(void)
    {
      ntype avL;
      acceptance=((ntype)(state.totmoves-state.totrej))/state.totmoves;
      traaccept = ((ntype)(state.tramoves-state.trarej))/state.tramoves;
      rotaccept = ((ntype)(state.rotmoves-state.rotrej))/state.rotmoves; 
      volaccept = ((ntype)(state.volmoves-state.volrej))/state.volmoves;
      if (pars.adjsteps < 0 || pars.curStep <= pars.adjsteps)
        { 
          if (pars.targetAccept > 0.0 && state.curStep % state.resetaccept==0)
            {
              if (traaccept > pars.targetAccept)
                pars.deltra *= 1.1;
              else
                pars.deltra /= 1.1;
#if 0
              if (state.restrmove==0)
                {
                  if (rotaccept > pars.targetAccept)
                    pars.delrot *= 1.1;
                  else
                    pars.deltor /= 1.1;
                }
#endif
              if (pars.deltra > (avL=pow(pars.L[0]*pars.L[1],1.0/2.0))*0.1)
                pars.deltra = avL*0.1;
              if (pars.delrot > 3.14)
                pars.delrot = 3.14;
            }

          if (pars.targetAcceptVol > 0.0 && state.volmoves > 0     
              && (state.curStep % state.resetacceptVol==0))
            {
              if (volaccept > pars.targetAcceptVol)
                pars.vmax *= 1.1;
              else
                pars.vmax /= 1.1;
            }
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
      ellLL.init(&he,&rect,pars.L,pars.rcut);
    }
  void build_LL(void)
    {
      ellLL.build();
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
      state.totmoves++;
      if (p <= 0.5)
        {
          he[ip].tra_move(pars.deltra);
          update_rect(ip);
          state.tramoves++; 
          return 0;
        }
      else
        {
          he[ip].rot_move(pars.delrot);
          update_rect(ip);
          state.rotmoves++;
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
      auto movetype=random_move(ip);
      // check overlap by using LL here!
      pbc(ip);
      ellLL.update_LL(ip);
      if (ellLL.check_overlap(ip)==true)
        {
          restore_coord(ip);
          ellLL.update_LL(ip);
          if (movetype==0)
            state.trarej++;
          else if (movetype==1)
            state.rotrej++;
          state.totrej++;
          // reject move
        }
    }
  void move_box(void)
    {
      //NPT volume move
      int i, ii, dir=0;
      pvector<ntype,2> L2;
      ntype vo, lnvn, vn, Lfact, arg;
      state.volmoves++;
      if (pars.simtype==2)
        {
          dir=3*ranf();
        }
      vo = pars.L[0]*pars.L[1];
      lnvn = log(vo) + (ranf()-0.5)*pars.vmax;
      vn = exp(lnvn);
      //cout << "vo=" << vo << " vn=" << vn << " vmax=" << state.vmax << "\n";
      if (pars.simtype==2)
        {
          Lfact = vn/vo;
          pars.L[dir] *= Lfact;
          L2[dir] = 0.5*pars.L[dir];
        }
      else
        {
          Lfact = pow(vn/vo,1.0/2.0);
          //cout << "Lfact=" << Lfact << "\n";
          pars.L *= Lfact;
          L2 = 0.5*pars.L;
        }
      //ellLL.store();
      //cout << "Lnew=" << pars.L[0] << " " << pars.L[1] << "\n";
      for (i=0; i < pars.parnum; i++)
        {
          if (pars.simtype==2)
            {
              he[i].r[dir] *= Lfact;
            }
          else
            {
              he[i].r *= Lfact;
            }
          update_rect(i);
          pbc(i);
        }
      ellLL.update_numcells(pars.L);
      ellLL.build();
      for (i=0; i < pars.parnum; i++)
        {
          if (ellLL.check_overlap(i)==true)
            {
              /* move rejected restore old positions */
              if (pars.simtype==2)
                {
                  pars.L[dir] /= Lfact;
                  L2[dir] = 0.5*pars.L[dir];
                }
              else
                {
                  pars.L /= Lfact;
                  L2 = 0.5*pars.L;
                }
              for (ii=0; ii < pars.parnum; ii++)
                {
                  if (pars.simtype==2)
                    he[ii].r[dir] /= Lfact;
                  else 
                    he[ii].r /= Lfact;
                  update_rect(ii);
                  pbc(ii);
                }
              //volrejMC++;
              ellLL.update_numcells(pars.L);
              build_LL();
              state.volrej++;
              return;
            }
        }
      arg = -(1.0/pars.T)*(pars.P*(vn-vo)-(pars.parnum+1)*log(vn/vo)*pars.T);
      if (ranf() > exp(arg))
        {
          /* move rejected restore old positions */
          if (pars.simtype==2)
            {
              pars.L[dir] /= Lfact;
              L2[dir] = 0.5*pars.L[dir];
            }
          else
            {
              pars.L /= Lfact;
              L2 = 0.5*pars.L;
            }
          for (i=0; i < pars.parnum; i++)
            {
              if (pars.simtype==2)
                he[i].r[dir] /= Lfact;
              else
                he[i].r /= Lfact;
              pbc(i);
            }
          //volrejMC++;
          ellLL.update_numcells(pars.L);
          build_LL();
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
      of.open("phi.dat", ios::trunc);
      of.close();
    }
  void save_measures(int tt)
    {
      ofstream of;
      of.open("v.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      if (tt!=0 && tt%pars.savemeas==0)
        {
          of << tt << " " << setprecision(16) << pars.L[0]*pars.L[1] << "\n";
        }
      of.close();
      of.open("phi.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      if (tt!=0 && tt%pars.savemeas==0)
        {
          of << tt << " " << setprecision(16) << pars.parnum*M_PI*pars.a*pars.b/(pars.L[0]*pars.L[1]) << "\n";
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
          state.curStep=tt;
          ntot=pars.parnum;
          if (pars.simtype==1) //NPT
            ntot+=1;
          ellLL.set_ini_numcells();
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
          calc_acceptance_and_adjust();
          // save snapshots and backup file 
          save_snap(tt);
          save_restore(tt);
          save_measures(tt);
          //output
          outinfo(tt);
        }
    }
  void outinfo(long long int tt)
    {
      if (tt> 0 && tt % pars.outsteps == 0)
        {
          cout << "steps=" << tt << "\n";
          cout << "acceptance tra=" << traaccept << " rot=" << rotaccept << "\n";
          if (pars.simtype==1 || pars.simtype==2)
            cout << "vol acceptance=" << volaccept << "\n";
 
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
              //pars.L.show("box");
              for (ix = 0; ix < maxix; ix++)
                {
                  for (iy = 0; iy < maxiy; iy++)
                    {
                      he[i].r << x-pars.L[0]*0.5,y-pars.L[1]*0.5;
                      //he[i].r.show("r");
                      he[i].na.random_orient();
                      he[i].nb = he[i].na.orto();
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
