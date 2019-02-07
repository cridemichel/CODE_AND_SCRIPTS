#include "simulation.H"
template <class ntype, class particle, class hardobj, class params> 
class simulmc: public sim<ntype,particle,hardobj,params>
{
  using sim = sim<ntype,particle,hardobj,params>;
  void mcexc(int& err) 
    {err=0;}
  void mcoutin(ntype beta, ntype pbias)
    {

    }
  void mcinout(ntype beta, ntype pbias)
    {

    }
  void move_box_xyz(int& err)
    {

    }
  void move_box(int& err)
    {}
  virtual void onestep(void);

 int mcmotion(void); 
public:

simulmc() = delete;
simulmc(int n): sim(n)
  {

  } 
};

template <class ntype, class particle, class hardobj, class params>
int simulmc<ntype,particle,hardobj,params>::mcmotion(void)
{
  int ip, dorej, movetype, err;
  double enn, eno;
  if (N==0)
    return 0;
  ip = N*ranf();
  /* qui basta calcolare l'energia della particella che sto muovendo */
  eno = calcpotene_GC(ip);
  //printf("eno=%f\n", eno);
  store_coord(ip);
  movetype=random_move(ip);

  pbc(ip);
  update_LL(ip);
 
  state.totmovesMC++;
  /* overlapMC() aggiorna anche i bond */
  //err=0;
  dorej = overlapMC(ip, &err);
  if (!dorej)
    {
#ifdef MC_STOREBONDS
      store_bonds_mc(ip);
#endif
      update_bonds_MC(ip);

      /* qui basta calcolare l'energia della particella che sto muovendo */
      enn=calcpotene_GC(ip);
      if (enn <= eno)
	{
	  //printf("ene=%f (old=%f)\n", enn, eno);
	  //	  if (abs(enn-eno) >=1 )
	  //	    printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	  dorej=0;
	}
      else
	{
	  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
	    dorej=0;
	  else
	    {
	      dorej=2;
	    }
	}
    }

  if (dorej != 0)
    {
      /* move rejected */
      totrejMC++;
      if (movetype==0)
	trarejMC++;
      else 
	rotrejMC++;
      //printf("restoring coords\n");
      if(err)
	{
	  cout << "[random_move] NR failed...I rejected this trial move...\n";
	  err=0;
	}
      restore_coord(ip);
      //rebuildLinkedList();
      update_LL(ip);
      if (dorej==2)
	{
#ifdef MC_STOREBONDS
	  remove_allbonds_ij(ip, -2);
	  restore_bonds_mc(ip);
#else
	  update_bonds_MC(ip);
#endif
	}
    }
  if (OprogStatus.useNNL && dorej==0 )
    {
      overestimate_of_displ[ip] += displMC; 
    }
  return movetype;
}
template <class ntype, class particle, class hardobj, class params>
void simulmc<ntype,particle,hardobj,params>::onestep(void)
{
  //if (state.particles.size()==0)
   // return 0;
  int i, movetype, err;
  int deln, ran, ntot;
  for (i=0; i < N; i++)
    {
      if (pars.ensemble==0 && deln==0)
	ran = 0;
      else
	ran = ntot*ranf();
      if (pars.ensemble == 4)
	{
	  movetype = mcmotion();
	}	
      else if (state.ensemble==2 && ran >= state.npav)
	{
	  if (ran >= pars.npav + pars.nexc)
	    {
	      if (ranf() < state.pbias)
		mcoutin(1.0/pars.T,state.pbias);
	      else
		mcinout(1.0/pars.T,state.pbias);
	    }
	  else
	    {
	      mcexc(&err);
	      movetype=4; /* 4 = insert/remove particle */
	      if(err)
		{
		  cout << "[mcexc] NR failed...I rejected this trial move...\n";
		  err=0;
		}
	      state.excmoveMC++;
	    }
	}
      else if ((state.ensemble==1 || state.ensemble==3||state.ensemble==4)
	  && ran >= N)
	{
	  //err=0;
	  if (ran > N)
	    {
	      if (ranf() < state.pbias)
		mcoutin(1.0/pars.T,state.pbias);
	      else
		mcinout(1.0/pars.T,state.pbias);
	    }
	  else
	    {
	      //check_bonds_mc("boh");
	      if (state.ensemble==3)
		move_box_xyz(&err);
	      else 
		move_box(&err);
	      movetype=3; /* 0 = tra; 1 = rot 2 = tra and rot; 3 = box */
	      if(err)
		{
		  cout << "[move_box] NR failed...I rejected this trial move...\n";
		  err=0;
		}
	      state.volmoveMC++;
	    }
	}
      else
	{
	  if (state.ensembleMC==0 && ran >= N)
	    {
	      if (ranf() < state.pbias)
		mcoutin(1.0/pars.T,state.pbias);
	      else
		mcinout(1.0/pars.T,state.pbias);
	    }
	  else
	    {
	      movetype=mcmotion();
	    }
	}

    }
}
