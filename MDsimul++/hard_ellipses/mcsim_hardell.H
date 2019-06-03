#include "../mdlib/pvector.H"
#include<string>
#include<iostream>
#include<fstream>
#include "./hardell.H"
#include "../mdlib/boxes.H"
#include "../mdlib/linked_cell_lists.H"
#include <time.h>
#include "./params.H"
/* N.B. per il restart basta prendere o il file di restart ascii più recente
 * oppure quello che ha il numero di campi maggiori calcolati così
 * NF=`awk '{numf+=NF} END {print numf}'`*/
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
  vector<pvector<ntype,2>> r0;
  vector<ntype> dth;
  vector<pvector<int,2>> cross;
  linked_cell_lists_2d<hardell<ntype>,rectangle<ntype>> ellLL;
  pvector<ntype,2> rold;
  pvector<ntype,2> naold, nbold;
  ntype acceptance, traaccept, rotaccept, volaccept, msd, dth_i;
  pvector<int,2> crossed;
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
      if (pars.adjsteps >= 0 || pars.curStep <= pars.adjsteps)
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
      resize_vec_N();
      init_LL();
    }
  void adjust_cross(int ip)
    {
      cross[ip] += crossed;
    }
  void pbc(int ip)
    {
      pvector<ntype,2> delr;
      //L2 = 0.5*pars.L;
      delr << pars.L[0]*rint(he[ip].r[0]/pars.L[0]), pars.L[1]*rint(he[ip].r[1]/pars.L[1]);
      he[ip].r = he[ip].r - delr;
      //delr.show("delr");
      crossed << 0,0;
      if (delr[0]!=0)
        crossed[0] = (int)copysign(1.0,delr[0]);
      if (delr[1]!=0)
        crossed[1] = (int)copysign(1.0,delr[1]);
    }
  void update_rotmsd(int ip)
    {
      dth[ip] += dth_i;
    }
  int random_move(int ip) // return move type 0=tra 1=rot
    {
      ntype p;
      p=ranf();
      state.totmoves++;
      if (p <= 0.5)
        {
          he[ip].tra_move(pars.deltra);
          pbc(ip);
          update_rect(ip);
          state.tramoves++; 
          return 0;
        }
      else
        {
          dth_i=he[ip].rot_move(pars.delrot);
          //cout << "dth=" << dth_i << "\n";
          pbc(ip);
          update_rect(ip);
          state.rotmoves++;
          //update_rotmsd(ip);    
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
      //cout << "attempt to move ip=" << ip << "\n";
      bool reject=false;
      store_coord(ip);
#if 0
      if (ellLL.check_overlap(ip)==true)
        {
          cout << "prima boh ip=" << ip << "\n";
          exit(0);
        }
#endif

#if 0
      he[ip].r.show("prima r=");
      he[ip].shift.show("prima shift="); 
      he[ip].na.show("na=");
      he[ip].nb.show("na=");
#endif       
      auto movetype=random_move(ip);
      // check overlap by using LL here!
      //pbc(ip);
      ellLL.update_LL(ip);
      if (ellLL.check_overlap(ip)==true)
        {
          reject=true;
          restore_coord(ip);
          update_rect(ip);
          ellLL.update_LL(ip);
          //cout << "incell=" << ellLL.inCell[0][ip] << " " << ellLL.inCell[1][ip] << "\n";
          if (movetype==0)
            state.trarej++;
          else if (movetype==1)
            state.rotrej++;
          state.totrej++;
#if 0
          if (ellLL.check_overlap(ip)==true)
            {
              he[ip].r.show("dopo r="); 
              he[ip].shift.show("dopo shift="); 
              he[ip].na.show("na=");
              he[ip].nb.show("na=");
              cout << "boh ip=" << ip << "\n";
              exit(0);
            }
#endif
          // reject move
        }
      if (!reject)
        {
          adjust_cross(ip);
          update_rotmsd(ip);
        } 
      //cout << "attempt to move ip=" << ip << "END\n";
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
#if 0
      for (ii=0; ii < pars.parnum; ii++)
        {
          if (ellLL.check_overlap(ii)==true)
            {
              cout << "prima boh ip=" << ii << "\n";
              exit(0);
            }
        }
#endif
  
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
#if 0
          if (Lfact < 1.0)
            return;
#endif
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
          pbc(i);
          update_rect(i);
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
                  pbc(ii);
                  update_rect(ii);
                }
              //volrejMC++;
              ellLL.update_numcells(pars.L);
              build_LL();
#if 0
              for (ii=0; ii < pars.parnum; ii++)
                {
                  if (ellLL.check_overlap(ii)==true)
                    {
                      cout << "Lfact=" << Lfact << "\n";
                      cout << "dopo boh ip=" << ii << "\n";
                      exit(0);
                    }
                }
#endif
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
  void assign_semiaxes(void)
    {
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
  void read_coords(fstream& i)
    {
      i >> pars.N;
      i >> pars.L[0];
      i >> pars.L[1];
      i >> pars.a;
      i >> pars.b;
      pars.rcut = 2.0*max(pars.a,pars.b);
      resize_vec_N();
      init_LL();
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
          rect[ip].sax[0] = pars.a;
          rect[ip].sax[1] = pars.b;
          update_rect(ip);
        }
    }
  void write_coords(fstream& o)
    {
      o << pars.N << " ";
      o << setprecision(20) << pars.L[0] << " " << pars.L[1] << " " << pars.a << " " << pars.b << "\n";
      for (auto i=0; i < pars.parnum; i++)
        {
          // save axis a
          o << setprecision(20) << he[i].r[0] << " " << he[i].r[1] << " " << he[i].na[0] << " " << he[i].na[1] << "\n";
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
  void write_msd(fstream &o)
    {
      for (auto i=0; i < pars.N; i++)
        {
          o << setprecision(20) << r0[i][0] << " " << r0[i][1] << " " << cross[i][0] << " " 
            << cross[i][1] << " " << dth[i] <<"\n";
        }
    }
  void read_msd(fstream &is)
    {
      for (auto i=0; i < pars.N; i++)
        {
          is >> r0[i][0];
          is >> r0[i][1];
          is >> cross[i][0];
          is >> cross[i][1];
          is >> dth[i];
        }
    }
  void save_restart(long long int tt)
    {
      static unsigned int w=0;
      state.curStep = tt;
      if (tt != iniStep && pars.saverest > 0 && tt % pars.saverest == 0)
        {
          sprintf(rfn, "restart-%u", w);
          restorefile.open(rfn,ios::out|ios::trunc);
          pars.write(restorefile);
          state.write(restorefile);
          write_coords(restorefile);
          write_msd(restorefile);
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
      of.open("msd.dat", ios::trunc);
      of.close();
      of.open("rotmsd.dat", ios::trunc);
      of.close();
    }
  void save_measures(long long int tt)
    {
      ofstream of;
      of.open("v.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      if (tt!=iniStep && tt%pars.savemeas==0)
        {
          of << tt << " " << setprecision(16) << pars.L[0]*pars.L[1] << "\n";
        }
      of.close();
      of.open("phi.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      if (tt!=iniStep && tt%pars.savemeas==0)
        {
          of << tt << " " << setprecision(16) << calc_phi() << "\n";
        }
      of.close();
      of.open("msd.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      msd=0;
      pvector<ntype,2> dr;

      for (auto i=0; i < pars.N; i++)
        {
          dr[0] = ntype(cross[i][0])*pars.L[0];
          dr[1] = ntype(cross[i][1])*pars.L[1];
          dr += he[i].r - r0[i];
          msd+=dr*dr;
        }
      msd /= ntype(pars.N);
      if (tt!=iniStep && tt%pars.savemeas==0)
        {
          of << tt << " " << setprecision(16) << msd << "\n";
        }
      of.close();
      of.open("rotmsd.dat",ios::app);
      // save observables, e.g. energy, volume, pressure, etc...
      ntype dthtot=0;
      for (auto i=0; i < pars.N; i++)
        {
          dthtot+=dth[i]*dth[i];
        }
      dthtot /= ntype(pars.N);
      if (tt!=iniStep && tt%pars.savemeas==0)
        {
          of << tt << " " << setprecision(16) << dthtot << "\n";
        }
      of.close();
      if (tt!=iniStep && tt%pars.savemeas==0)
       calc_hex_ordpar(3*pars.b);
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
  void resize_vec_N(void)
    {
      int Nold = rect.size();
      if (Nold < pars.N)
        {
          rect.resize(pars.N);
          he.resize(pars.N);
          r0.resize(pars.N);
          cross.resize(pars.N);
          dth.resize(pars.N);
        }
    }
  void init(void)
    {
      if (pars.seed >= 0)
        srand48(pars.seed);
      else
        srand48(time(NULL));
      //cout << "deltra=" << pars.deltra << " delrot=" << pars.delrot << " vmax=" << pars.vmax << "\n";
      //cout << "P= " << pars.P << " T=" << pars.T << "\n";
      pars.rcut = 2.0*max(pars.a,pars.b);
      //resize_vec_N();
      //init_LL();
      iniStep=0;
      msd=0;
    }

public:
  long long int iniStep;
  int boot_type;//0=create conf 1=restart 2=read conf
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

  void check_overlap_debug(void) 
    { 
      int i;
      cout << "Checking overlap\n";
      for (i=0; i < pars.parnum; i++)
        {
          if (ellLL.check_overlap(i,-2,true)==true)
            {
              cout << "step #" <<  pars.curStep << " boh overlap i=" << i << "\n";
              save_snap_str("cnf-final");
              exit(1);
            }
        }
    }
  void prepare_observables(void)
    {
      if (boot_type!=1) // 1=restart
        {
          for (auto i=0; i < pars.N; i++)
            {
              r0[i] = he[i].r;
              cross[i] << 0,0;
              //cerr << "r(" << i << ")=" << r0[i][0] << " " << r0[i][1] << "\n"; 
            }
        }
#if 0
      else
        {
          for (auto i=0; i < pars.N; i++)
            {
              cerr << "r(" << i << ")=" << r0[i][0] << " " << r0[i][1] << "\n"; 
            } 
        }
#endif
    }
  void run(void)
    {
      long long int tt;
      int ntot, ran, i;
      //cout << "saverest=" << pars.saverest << "\n";
      cout << "N=" << pars.parnum << " a=" << pars.a << " b=" << pars.b << "\n";
      prepare_observables();
      for (tt=iniStep; tt < pars.steps; tt++)
        {
          state.curStep=tt;
          ntot=pars.parnum;
          if (pars.simtype==1) //NPT
            ntot+=1;
          ellLL.set_ini_numcells();

          //cout << ">>>tt=" << tt << "\n";
          //exit(1);
          for (i=0; i < ntot; i++)
            {
              //cout << "step #" << tt << " i=" << i << "ntot=" << ntot << "\n";
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
          if (pars.targetPhi > 0 && abs(calc_phi()-pars.targetPhi) < 1E-7)
            break;
          if (pars.maxmsd > 0 && msd > pars.maxmsd)
            {
              cout << "Equilibration reached (msd=" << msd << ")\n";
              break;
            } 
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
      read_msd(restorefile);
      restorefile.close();
      //cout << "state.curStep=" << state.curStep << "\n";
      cout << "semi-axes=" << pars.a << " " << pars.b << "\n";
      iniStep = pars.curStep;
      build_LL();
      boot_type=1; 
    }
  void readconf(string& icn)
    {
      init_measures();
      inicfile.open(icn.c_str(),ios::in);
      // read ascii iniconf here 
      read_coords(inicfile);
      inicfile.close();
      build_LL();
      boot_type=2;
    }
  void place_on_square_lattice(bool maxphi=true)
    { 
      ntype dx, dy;
      pvector<ntype,2> u1, u2, fact;
      //ntype x, y;
      int i=0, ix, iy, maxix, maxiy;
      maxix = maxiy = ceil(sqrt(((ntype) pars.N))); 
      cout << "maxix=" << maxix << " maxiy=" << maxiy << "\n";
      u1 << 1,0;
      u2 << 0,1;
      fact << pars.a,pars.b;
      fact *= 2.0;
      if (maxphi==true)
        {
          pars.L[0] = maxix*fact[0];
          pars.L[1] = maxiy*fact[1];
          dx = fact[0];
          dy = fact[1];
        }
      else
        {
          dx = pars.L[0]/((ntype)maxix);
          dy = pars.L[1]/((ntype)maxiy);
        }
      //x=y=0;
      cout << "dx=" << dx << " dy=" << dy << "\n";
      //pars.L.show("box");
      for (ix = 0; ix < maxix; ix++)
        {
          for (iy = 0; iy < maxiy; iy++)
            {
              //he[i].r.show("r");
              he[i].r = (ix*dx)*u1 + (iy*dy)*u2 - 0.5*pars.L;
              if (maxphi==false)
                {
                  //he[i].r << x-pars.L[0]*0.5,y-pars.L[1]*0.5;
                  he[i].na.random_orient();
                  he[i].nb = he[i].na.orto();
                }
              else
                {
                  he[i].na << 1,0;
                  he[i].nb << 0,1;
                }
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
              //y += dy;
            }
          //y=0;
          //x += dx;
        }
    }
  void expand(ntype fx, ntype fy, bool save=true, bool updsemiax=false)
    {
      ellLL.set_ini_numcells();
      
      readconf(pars.iniconf);
      cout << "L=" << pars.L[0] << " " << pars.L[1] << "\n";
      pars.L[0] *= fx;
      pars.L[1] *= fy;
      cout << "expanding conf " << pars.iniconf << " by factors " << fx << " " << fy << "\n";
      for (auto i=0; i < pars.N; i++)
        {
          he[i].r[0] *= fx;
          he[i].r[1] *= fy;
          update_rect(i);
        }
      if (updsemiax)
        {
          pars.a *= fx;
          pars.b *= fy;
        }
      ellLL.update_numcells(pars.L);
      build_LL();
      if (save)
        save_snap_str("cnf-expanded");
    }
  void create_triang_lattice(int maxix=0, int maxiy=0, bool save=true)
    {
      ntype dx, dy, sig, exfact;
      pvector<ntype,2> u1, u2, Ls, L;
      ntype fact, phim, K;
      bool fine=false;
      //ntype x, y;
      int i=0, ix, iy;
      // assum pars.a > pars.b 

      //ellLL.set_ini_numcells();
      K=sqrt(3.0);
      if (maxix > 0 && maxiy > 0)
        {
          pars.N=maxix*maxiy*2;
        }
      else
        {
          //estimate maxix and maxiy from N
          maxiy = rint(sqrt(ntype(pars.N)*(pars.a/(K*pars.b))/2.0));
          maxix = rint(ntype(maxiy)*K*pars.b/pars.a);
        }
      resize_vec_N();
      if (pars.N > maxix*maxiy*2)
        {
          cerr << "ERROR N > available lattice sites\n";
          exit(1); 
        }
      Ls << maxix*pars.a*2.0,K*maxiy*pars.b*2.0;
      cout << ">>>Lest=" << Ls[0] << " " << Ls[1] << "\n";
      cerr << "maxix=" << maxix << " maxiy=" << maxiy << "\n";
      u1 << 1,0;
      u2 << 0,sqrt(3.0);
      if (pars.a < pars.b)
        {
          cout << "[ERROR a < b] a must be greater than b!\n";
          exit(1);
        }
      sig = 2.0*pars.b;
      exfact = 2.0*pars.a;
      dx = sig;
      dy = sig;
      phim = M_PI/(sqrt(3.0)*2.0);
      cerr << setprecision(16) << "phi max=" << phim << "\n";
      cerr << "dx=" << dx << " dy=" << dy << "\n";
      cerr << "a=" << pars.a << " b=" << pars.b << "\n";

      pvector <ntype,2> offset;
      offset << 0,0;
      for (auto nl=0; nl < 2 && fine==false; nl++)
        {
          if (nl==1)
            offset << 0.5*sig, 0.5*sig*sqrt(3.0); 
          for (ix = 0; ix < maxix && fine==false; ix++)
            {
              for (iy = 0; iy < maxiy && fine==false; iy++)
                {
                  he[i].r = offset + (ix*dx)*u1 + (iy*dy)*u2;
                  he[i].r[0] *= exfact;
                  he[i].na << 1,0;
                  he[i].nb << 0,1;
                  he[i].a = pars.a;
                  he[i].b = pars.b;
                  i++;
                  if (i >= pars.N)
                    { 
                      fine = true;
                    }
                }
            }
        }
      if (pars.N < maxix*maxiy*2)
        cout << "[WARNING] number of particles (N=" << pars.N << ") is not optimal being less than maxix*maxiy*2=" << maxix*maxiy*2 << "\n";
      pars.N=i;
      pars.L << 0,0;
      for (i=0; i < pars.N; i++)
        {
          if (i==0 || he[i].r[0]+pars.a > pars.L[0])
            pars.L[0]=he[i].r[0]+pars.a;
          if (i==0 || he[i].r[1]+pars.b*K > pars.L[1])
            pars.L[1]=he[i].r[1]+pars.b*K;
        }
      cerr<< "L= " << pars.L[0] << " " << pars.L[1] << "\n";

      fact=ntype(1.0+1E-8);
      for (i=0; i < pars.N; i++)
        { 
          he[i].r -= ntype(0.5)*pars.L;
          he[i].r *= fact;
          update_rect(i);
        }
      pars.L *= fact;
      auto phi=calc_phi();
      cerr << "phi=" << phi << " N=" << pars.N << "\n";
      if (save)
        save_snap_str("cnf-trianglattice");
      //ellLL.update_numcells(pars.L);
      init_LL();
      build_LL();
    }
#if 0
  void place_on_triang_lattice(bool randorient=false)
    { 
      ntype dx, dy, sig,exfact;
      pvector<ntype,2> u1, u2, Ls;
      ntype fact, phim;
      //ntype x, y;
      bool fine=false;
      int i=0, ix, iy, maxix, maxiy;
      // assum pars.a > pars.b 
      maxix = maxiy = ceil(sqrt(((ntype) pars.N))); 
      cout << "maxix=" << maxix << " maxiy=" << maxiy << "\n";
      u1 << 1,0;
      u2 << -0.5,sqrt(3.0)/2.0;
      sig = 2.0*pars.b;
      exfact = 2.0*pars.a;
      dx = sig;
      dy = sig;
      phim = M_PI/(sqrt(3.0)*2.0);
      cout << setprecision(16) << "phi max=" << phim << "\n";
      Ls[0] = sqrt(pars.N*M_PI*pars.a*pars.b/phim);
      Ls[1] = sqrt(pars.N*M_PI*pars.a*pars.b/phim);
      cout << "Ls = " << Ls[0] << " " << Ls[1] << "\n";
      //x=y=0;
      cout << "dx=" << dx << " dy=" << dy << "\n";
      cout << "a=" << pars.a << " b=" << pars.b << "\n";
      //pars.L.show("box");

      ellLL.set_ini_numcells();
      for (ix = -3*maxix; ix < 3*maxix; ix++)
        {
          for (iy = -3*maxiy; iy < 3*maxiy; iy++)
            {
              //he[i].r.show("r");
              //cout << "ix=" << ix << " iy=" << iy << "\n";
              he[i].r = (ix*dx)*u1 + (iy*dy)*u2;
              he[i].r[0] *= exfact;
              if (he[i].r[0] <  0 || he[i].r[0] >= Ls[0]
                  || he[i].r[1] < 0 || he[i].r[1] >= Ls[1])
                {
                  continue;
                }

              //cerr << he[i].r[0] << " " << he[i].r[1] << "\n";
              if (randorient==true)
                {
                  //he[i].r << x-pars.L[0]*0.5,y-pars.L[1]*0.5;
                  he[i].na.random_orient();
                  he[i].nb = he[i].na.orto();
                }
              else
                {
                  //he[i].r[maxax] *= fact;
                  he[i].na << 1,0;
                  he[i].nb << 0,1;
                }
              he[i].a = pars.a;
              he[i].b = pars.b;
              i++;
              //cout << "i=" << i << " pars.N="<< pars.N <<"\n";
              if (i==pars.N)
                {
                  fine=true;
                  break;
                }
              //y += dy;
            }
          //y=0;
          //x += dx;
          if (fine==true)
            break;
        }
      if (i < pars.N)
        {
          cout << "Placed only " << i << " particles out of " << pars.N << "\n";
          exit(1); 
        }
      pars.L << 0,0;
      for (i=0; i < pars.N; i++)
        {
          if (i==0 || he[i].r[0]+2.0*pars.a > pars.L[0])
            pars.L[0]=he[i].r[0]+2.0*pars.a;
          if (i==0 || he[i].r[1]+2.0*pars.b*sqrt(3.0) > pars.L[1])
            pars.L[1]=he[i].r[1]+2.0*pars.b*sqrt(3.0);
        }
      cout << "L= " << pars.L[0] << " " << pars.L[1] << "\n";
      ntype phiT = 0.875;
      auto phi=calc_phi();

      fact = sqrt(phi/phiT); 
      cout << "phi=" << phi << " fact=" << fact << "\n"; 
      for (i=0; i < pars.N; i++)
        { 
          he[i].r -= 0.5*pars.L;
          he[i].r *=fact;
          update_rect(i);
          //he[i].show("mah");
        }
      pars.L *= fact;
      //cout << "quiqui\n";
      ellLL.update_numcells(pars.L);
      build_LL();
      //init_LL();
    }
#endif
  void createconf(int type)
    { // type = 0 SC  1 = TRIANGULAR
      resize_vec_N();
      init_LL();
      init_measures();
      boot_type=0;
      switch (type)
        {
        case 0:
            {
              place_on_square_lattice(false);  
              check_overlap_debug();
              break;
            }
        case 1:
            {
              create_triang_lattice();
              check_overlap_debug();
              cout << "no overlap in initial conf\n";
              break;
            }
        default:
            {

            }
        }
    }
  void calc_hex_ordpar(ntype rcut)
    {
      int i;
      vector<int> lnn;
      ntype theta;
      complex<ntype> psi;
      psi=0;
      for (i=0; i < pars.N; i++)
        {
          ellLL.get_neighbours(i);
          for (auto & el: ellLL.lista_vicini)
            {
              auto d = he[i].r - (he[el].r + he[el].shift);
              if (d*d < rcut*rcut)
                {
                  //cout << "found i=" << el << " dd=" << d*d << " rcut*rcut=" << rcut*rcut << "\n";
                  theta=0.0;
                  psi += complex<ntype>(cos(6.0*theta),sin(6.0*theta));
                }
            }
          //cout << "fine\n";
        }
      psi /= (ntype)pars.N;
      return abs(psi);
    }
  void calc_transl_order(void)
    {
      const pvector<ntype,2> G;
      complex<ntype> psi;
      psi=0;
      G << 2.0*M_PI/pars.L[0], 2.0*M_PI/pars.L[1];
      for (auto i=0; i < pars.N; i++)
        {
          psi += complex<ntype>(cos(G*he[i].r),sin(G*he[i].r));
        }
      psi /= (ntype)pars.N;
      return abs(psi);
    }
};
