#include "pmatrix.H"
#include "./linked_cell_lists.H"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
//#include "./phardcyl.H"
typedef long long int llint;
// patch will have a standard definition and a vector of patches
// will be a member of any particle object
template <class ntype>
class patch
{
public:
  pvector<ntype,3> r;
  int pidx;
};

 
#if 0
template <class ntype> class extrapars
{
public:
  extrapars()
    {
    }
  ~extrapars()
    {

    }
};
#endif
// particles are patchy hard cylinders here!
//template<class ntype> 
//using particle=phardcyl<ntype>;

template <class ntype, class particle> class simstate
{
  //extrapars<ntype> epars;
  int curstep, maxstep;
  int ensemble;
  llint excmoveMC, volmoveMC, totmovesMC;
public:
  llint getcurstep(void)
    {
      return curstep;
    }
  llint getfinalstep(void)
    {
      return maxstep;
    }
  simstate()
    {
      curstep=1;
    }
  ~simstate()
    {

    }
};

template <class ntype, class particle, class hardobj, class params> 
class sim
{
  simstate<ntype,particle> state;
  params pars;
  //extrapars<ntype> epars; 
  /* auxiliary vectors */
  vector<particle> particles;
  vector<patch<ntype>*> patches;
  vector<hardobj*> hardobjs;
  int NP; //num. of tot. patches 
  int NHO;//num. of tot. hard objects
  linked_cell_list<ntype,vector<patch<ntype>*>> lclpa;
  linked_cell_list<ntype,vector<hardobj*>> lclho;

public:
  int N;  //actual number of simulated particles 
  int Nalloc;// number of allocated particles (Nalloc >= N)
  void run(void);
  llint getinistep(void)
    {
      return state.getcurstep();
    }
  llint getfinalstep(void)
    {
      return state.finalstep();
    }
  ntype ranf(void)
    {
      return ((ntype)drand48());
    }
  void seed(int s)
    {
      srand48((long)s);
    }
  virtual void onestep(void);
  void saveobserv(void)
    {}
  void terminate(void)
    {
      // save last conf here
      exit(0);
    }
  void savesnap(int type=0);
  void get_patches()
    {
      for (auto &p: state.particles) 
        {
          for (auto &pa: state.particles.patches)
            {
              patches.push_back(&pa);
            }
        }
    }
  void get_hard_objs()
    {
      for (auto &p: state.particles) 
        {
          for (auto &ho: state.particles.hardobjs)
            hardobjs.push_back(&p);
        }
    }  
  int get_num_patches()
    {
      int totHO=0;
      for (auto &p: state.particles) 
        totHO+=p.get_num_hardobjs();
      return totHO;
    }
  int get_num_hard_objs()
    {
      int totNP=0;
      for (auto &p: state.particles) 
        totNP+=p.get_num_patches();
      return totNP;
    }
  void resize(int n)
    {
      // number of particles must be set elsewhere
      // this is for SUS simulations where size of allocagted vectors is conveniently greater than actual number of particles
      state.particles.resize(n);
      patches.resize(NP=get_num_patches());
      hardobjs.resize(NHO=get_num_hard_objs());
    }
  sim() = delete;
  sim(int n)
    {
      N=n;
      state.particles.resize(n);
      NP=get_num_patches();
      NHO=get_hard_objs();
      patches.resize(NP);
      hardobjs.resize(NHO);
    }
  void init_lcls(void)
    {
      lclpa.init_lcls(pars.Lbox, pars.rcutpa, patches);
      lclho.init_lcls(pars.Lbox, pars.rcutho, hardobjs);
    }
  void init(void)
    {
      seed(42);
      // place all particles if a configuration is not provided
      init_lcls();
    }
};

template <class ntype, class particle, class hardobj, class params>
void sim<ntype,particle,hardobj,params>::run(void)
{
  int i;
  for (i=getinistep(); i < getfinalstep(); i++)
    {
      onestep();
      saveobserv();
      savesnap();
    }
}


