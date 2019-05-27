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
  fstream snapfile;
  fstream restorefile;
  fstream inicfile;
  char sfn[512], rfn[512];
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
      if (pars.simtype==1 || pars.simtype==2)
        volaccept = ((ntype)(state.volmoves-state.volrej))/state.volmoves;
      if (pars.adjsteps < 0 || pars.curStep <= pars.adjsteps)
        { 
          if (pars.targetAccept > 0.0 && state.curStep % pars.resetaccept==0)
            {
              if (traaccept > pars.targetAccept)
                pars.deltra *= 1.1;
              else
                pars.deltra /= 1.1;
              if (rotaccept > pars.targetAccept)
                pars.delrot *= 1.1;
              else
                pars.delrot /= 1.1;
              if (pars.deltra > (avL=pow(pars.L[0]*pars.L[1],1.0/2.0))*0.1)
                pars.deltra = avL*0.1;
              if (pars.delrot > 3.14)
                pars.delrot = 3.14;
              state.totrej=state.trarej=state.rotrej=state.totmoves=state.rotmoves=state.tramoves=0;
            }

          if ((pars.simtype==1 || pars.simtype==2) && pars.targetAcceptVol > 0.0 && state.volmoves > 0     
              && (state.curStep % pars.resetacceptVol==0))
            {
              if (volaccept > pars.targetAcceptVol)
                pars.vmax *= 1.1;
              else
                pars.vmax /= 1.1;
              state.volmoves=state.volrej=0;
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
        pars.readpars(parsfilestr);
      else
        pars.readpars(fn);
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
  void check_overlap_debug(void) 
    { 
      int i;
      for (i=0; i < pars.parnum; i++)
        {
          if (ellLL.check_overlap(i)==true)
            {
              cout << "step #" <<  pars.curStep << "boh overlap i=" << i << "\n";
              exit(1);
            }
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
          Lfact = pow(vn/vo,0.5);
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
      //cout << "Lfact=" << Lfact << "\n";

      for (i=0; i < pars.parnum; i++)
        {
          if (ellLL.check_overlap(i,-2)==true)
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
  void read_coords(fstream& i)
    {
      for (auto ip=0; ip < pars.parnum; ip++)
        {
          // save axis a
          i >> he[ip].r[0];
          i >> he[ip].r[1];
          i >> he[ip].na[0];
          i >> he[ip].na[1];
          he[ip].nb = he[ip].na.orto();
          he[ip].a = pars.a;
          he[ip].b = pars.b;
          update_rect(ip);
        }
    }
  void write_coords(fstream& o)
    {
      for (auto i=0; i < pars.parnum; i++)
        {
          // save axis a
          o << setprecision(16) << he[i].r[0] << " " << he[i].r[1] << " " << he[i].na[0] << " " << he[i].na[1] << "\n";
        }
    } 
  void save_snap_str(const char *s)
    {
      snapfile.open(s, ios::out|ios::trunc);
      write_coords(snapfile);
      snapfile.close();
    }
  void save_snap(long long int tt)
    {
      if (pars.savesnap > 0 && tt > 0 && tt % pars.savesnap == 0)
        {
          sprintf(sfn, "cnf-%lld", tt);
          save_snap_str(sfn); 
        }
    }
  void save_restart(long long int tt)
    {
      static unsigned int w=0;
      state.curStep = tt;
      if (tt > 0 && pars.saverest > 0 && tt % pars.saverest == 0)
        {
          sprintf(rfn, "restart-%u", w);
          restorefile.open(rfn,ios::out|ios::trunc);
          pars.write(restorefile);
          state.write(restorefile);
          write_coords(restorefile);
          restorefile.close();
          w = (~w)&1;
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
  void save_measures(long long int tt)
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
          of << tt << " " << setprecision(16) << calc_phi() << "\n";
        }
      of.close();
    }
  void outinfo(long long int tt)
    {
      if (tt> 0 && tt % pars.outsteps == 0)
        {
          cout << "steps=" << tt << "\n";
          cout << "tra acceptance=" << traaccept << " [ deltra=" << pars.deltra << " ]\n";
          cout << "rot acceptance=" << rotaccept << " [ delrot=" << pars.delrot << " ]\n";
          if (pars.simtype==1 || pars.simtype==2||state.volmoves > 0)
            cout << "vol acceptance=" << volaccept << " [vmax=" << pars.vmax << " ]\n";
 
        } 
    }
  void init(void)
    {
      rect.resize(pars.N);
      he.resize(pars.N);
      if (pars.seed >= 0)
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
      iniStep=0;
    }

public:
  long long int iniStep;
  const char *parsfilestr="hepars.asc";
  params<ntype> pars;
  simstate<ntype> state;
  mcsim()
    {
      read_pars(parsfilestr);
      init();
    }
  mcsim(const char* fn)
    {
      read_pars(fn);
      init();
    }
  ntype calc_phi()
    {
      return pars.parnum*M_PI*pars.a*pars.b/(pars.L[0]*pars.L[1]) ;
    }

  void run(void)
    {
      long long int tt;
      int ntot, ran, i;
      //cout << "saverest=" << pars.saverest << "\n";
      for (tt=iniStep; tt < pars.steps; tt++)
        {
          state.curStep=tt;
          ntot=pars.parnum;
          if (pars.simtype==1) //NPT
            ntot+=1;
          ellLL.set_ini_numcells();

          cout << ">>>tt=" << tt << "\n";
          //exit(1);
          for (i=0; i < ntot; i++)
            {
              cout << "step #" << tt << " i=" << i << "ntot=" << ntot << "\n";
              ran = (int) (ntot*ranf());
              if (ran < pars.parnum)
                {
                  mcmotion(ran);
                }
              else 
                {
                  move_box();
                }
            } 
          //check_overlap_debug();
          calc_acceptance_and_adjust();
          // save snapshots and backup file 
          save_snap(tt);
          save_restart(tt);
          save_measures(tt);
          //output
          outinfo(tt);
        }

      save_snap_str("cnf-final");
    }
  
  // read restart file
  void restart(char *rfn)
    {
      restorefile.open(rfn);
      pars.read(restorefile);
      state.read(restorefile);
      pars.curStep = state.curStep;
      read_coords(restorefile);
      restorefile.close();
      //cout << "state.curStep=" << state.curStep << "\n";
      cout << "semi-axes=" << pars.a << " " << pars.b << "\n";
      iniStep = pars.curStep;
      build_LL();
    }
  void readconf(string& icn)
    {
      inicfile.open(icn.c_str(),ios::in);
      // read ascii iniconf here 
      read_coords(inicfile);
      inicfile.close();
      build_LL();
    }
 
  void createconf(int type)
    { // type = 0 SC  1 = TRIANGULAR
      switch (type)
        {
        case 0:
            {
              ntype dx, dy, x, y;
              int i=0, ix, iy, maxix, maxiy;
              maxix = maxiy = ceil(sqrt(((ntype) pars.N))); 
              cout << "maxix=" << maxix << " maxiy=" << maxiy << "\n";
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
                      he[i].a = pars.a;
                      he[i].b = pars.b;
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
