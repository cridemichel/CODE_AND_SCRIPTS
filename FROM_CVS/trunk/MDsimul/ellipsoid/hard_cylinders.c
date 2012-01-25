#ifdef MC_HC
#include<mdsimul.h>
const double saxfactMC[3]={0.85,0.68,0.68};
#ifdef MC_QUASI_CUBE
const double saxfactMC_QC[3]={0.832,0.832,0.832};
#endif
const int nfons=200;
extern void init_rng(int mdseed, int mpi, int my_rank);
#ifdef MC_SIMUL
#ifdef MC_STORELL
int *cellListMC;
#endif
#ifdef MC_STOREBONDS
#ifdef MD_LL_BONDS
long long int **bondsMC;
int *numbondsMC;
#else
int *numbondsMC, **bondsMC;
#endif
#endif
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG10(x) 
#define MD_DEBUG11(x) 
#define MD_DEBUG15(x) 
#define MD_DEBUG20(x)  
#define MD_DEBUG31(x) 
#define MD_DEBUG32(x) 
#define MD_DEBUG33(x) 
#define MD_DEBUG34(x) 
#define MD_DEBUG35(x)  
#define MD_DEBUG36(x) 
#define MD_DEBUG38(x) 
#define MD_DEBUG39(x) 
#define MD_DEBUG40(x) 
extern double calc_norm(double *vec);
extern void vectProdVec(double *A, double *B, double *C);
extern void print_matrix(double **M, int n);
extern void update_MSDrot(int i);
extern void update_MSD(int i);
#ifdef MD_SUPERELLIPSOID
extern int is_superellipse(int i);
extern void fdjacSE(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3]);
#endif
#ifdef MD_ASYM_ITENS
extern void calc_omega(int i, double *wwx, double *wwy, double *wwz);
extern void calc_angmom(int i, double **I);
extern void upd_refsysM(int i);
#endif
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
#ifdef MD_CALENDAR_HYBRID
extern int *linearLists;
extern int numevPQ, totevHQ, overevHQ;
#endif
#ifdef MD_SPOT_GLOBAL_ALLOC
extern double **ratA, **ratB;
extern double **ratAll;
#endif
extern void ProcessCellCrossingMLL(void);
extern void PredictEventMLL(int na, int nb);
extern void PredictEventMLL_NLL(void);
extern double *mapBheightFlex, *mapBhinFlex, *mapBhoutFlex, *mapSigmaFlex; 

extern double DphiSqA, DphiSqB, DrSqTotA, DrSqTotB;
extern double minaxA, minaxB, minaxAB;
extern int do_check_negpairs;
#ifdef MD_ASYM_ITENS
extern double **Ia, **Ib, **invIa, **invIb, **Iatmp, **Ibtmp;
#else
extern double Ia, Ib, invIa, invIb;
#endif
#ifdef MD_PATCHY_HE
extern void bumpSP(int i, int j, int ata, int atb, double* W, int bt);
extern void assign_bond_mapping(int i, int j);
#endif
#ifdef MD_ASYM_ITENS
extern double *phi0, *psi0, *costheta0, *sintheta0, **REt, **RE0, *angM, ***RM, **REtA, **REtB, **Rdot;
extern double cosEulAng[2][3], sinEulAng[2][3];
#endif
#ifdef MD_PATCHY_HE
extern struct LastBumpS *lastbump;
extern void check_all_bonds(void);
#else
extern int *lastbump;
#endif
#ifdef MD_SPHERICAL_WALL
extern int sphWall, sphWallOuter;
#endif
extern double *axa, *axb, *axc;
#ifdef EDHE_FLEX
extern double *a0I;
#endif
extern int *scdone;
extern double *maxax;
extern double calcDistNegNeighPlane(double t, double t1, int i, double *r1, double *r2, double *vecgsup, int calcguess, int calcgradandpoint, int *err, int nplane);
void calc_energy(char *msg);
#ifdef EDHE_FLEX
extern int *is_a_sphere_NNL;
#endif
extern double min3(double a, double b, double c);
extern double min(double a, double b);
extern double max3(double a, double b, double c);
extern double *lastupdNNL, *totDistDispl;
double rA[3], rB[3];
/* Routines for LU decomposition from Numerical Recipe online */
extern void ludcmpR(double **a, int* indx, double* d, int n);
extern void lubksbR(double **a, int* indx, double *b, int n);
extern void InvMatrix(double **a, double **b, int NB);
extern double invaSq[2], invbSq[2], invcSq[2];
extern double rxC, ryC, rzC, trefG;
extern int SolveLineq (double **a, double *x, int n); 
extern int calcdist_retcheck;
extern void comvel_brown (COORD_TYPE temp, COORD_TYPE *m);
extern void InitEventList (void);
#ifdef MD_HSVISCO
extern void calcT(void);
#endif
extern void writeAsciiPars(FILE* fs, struct pascii strutt[]);
extern void writeAllCor(FILE* fs, int saveAll);
extern struct nebrTabStruct *nebrTab;
extern double nextNNLrebuild;
extern void rebuildNNL(void);
extern void updrebuildNNL(int na);
extern void PredictEventNNL(int na, int nb);
extern void updAllNNL();
#ifdef MD_PATCHY_HE
extern int isSymItens(int i);
extern int locate_contactSP(int i, int j, double shift[3], double t1, double t2, double *evtime, int *ata, int *atb, int *collCode);
extern void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent);
extern int *mapbondsa;
extern int *mapbondsb;
#endif
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern long long accngA, accngB;
extern void print_matrix(double **M, int n);
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern char TXT[MSG_LEN];
extern void resetCM(int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void kinet(int Nm, COORD_TYPE* velx, COORD_TYPE* vely, 
		  COORD_TYPE* velz, COORD_TYPE VOL1);
extern void ScheduleEvent (int idA, int idB, double tEvent);
extern void NextEvent (void);
extern void distanza(int ia, int ib);
#ifdef MD_LXYZ
extern double invL[3];
#else
extern double invL;
#endif
extern double pi, Vz;
#ifdef MD_LXYZ
extern double L2[3];
#else
extern double L2;
#endif 
#ifdef MD_GRAVITY
extern double Lz2;
#endif
extern double W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB;
#ifdef MD_HSVISCO
extern double  DQxxOld, DQyyOld, DQzzOld, DQxyOld, DQyzOld, DQzxOld, DQxxOldKin, 
       DQyyOldKin, DQzzOldKin, DQxxOldHS, DQyyOldHS, DQzzOldHS, DQxxOldST, DQyyOldST, DQzzOldST,
       PxxKin, PyyKin, PzzKin, PxxHS, PyyHS, PzzHS, PxxST, PyyST, PzzST;
#endif
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
/* used by linked list routines */
#ifdef MD_GRAVITY
extern double g2, mgA, mgB;
#endif
extern double *lastcol;
extern double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
extern int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
    cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
#ifdef MD_EDHEFLEX_OPTNNL
extern int *inCell_NNL[3], *cellList_NNL;
extern double *rxNNL, *ryNNL, *rzNNL;
#endif
extern int evIdA, evIdB;
extern int parnumA, parnumB;
#ifdef MD_PATCHY_HE
extern int evIdC, evIdD, evIdE;
extern double *treeRxC, *treeRyC, *treeRzC;
#ifdef MD_LL_BONDS
extern long long int *bondscache, **bonds;
extern int *numbonds;
#else
extern int *bondscache, *numbonds, **bonds;
#endif
#endif
extern void newtDist(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3]);
extern void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[], double xb2[], 
	   int *nb);
extern void newt(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3]);
extern void rebuildCalendar(void);
extern void R2u(void);
extern void store_bump(int i, int j);
#endif
#if (defined(MC_SIMUL) || defined(MD_STANDALONE)) && 1
double scalProd(double *A, double *B);
double toteneini=0.0;
long long int ttini=0;
int covrestart = 0;
const int nmboxMC=5;
double totdist=0.0, distcc=0.0;
static double maxarg1,maxarg2;

#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

int check_convergence(double Told[3], double Tnew[3])
{
  double test=0.0;
  int i;
  for (i=0;i<3;i++) 
    {
      temp=(fabs(Tnew[i]-Told[i]))/FMAX(fabs(Tnew[i]),1.0); 
      //temp=(fabs(x[i]-xold[i]))/fabs(x[i]); 
      if (temp > test) 
	test=temp; 
    }
  if (test < 1.0E-8)
    return 1;
  else 
    return 0;
}

double calcDistNegHC(int i, int j, double shift[3])
{
  const int MAX_ITERATIONS = 1000;
  int it;
  double ViVj[3], lambdai, lambdaj;
  double sp, Q1, Q2, normPiDi, normPjDj, normN, L, D, DiN, DjN, niN[3], njN[3], Djni, Djnj;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DiCi[3], DiCini, normDiCi, DjCi[3], normDjCi;
  double PiDi[3], PjDj[3], Ai[3], Tj[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3], TipCjnj, TimCjnj;
  double Aj[3], AjDini, AjDinivec[3], AjDi[3], Tip[3], Tim[3], TipCj[3], TimCj[3], Dini;
  double DiCj[3], normDiCj, DiCjnj, Uj[3], DiUj[3], normDiUj, DiUjnj;
  double Tim_perp[3], Tip_perp[3], Tim_para[3], Tip_para[3], normTim_perp, DjCini;
  double Tjm_perp[3], Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp;
  double TiOld[3], TiNew[3], TiNewCj[3], TiNewCjnj;	
  double normCiCj;	
  int kk, j1, j2;

  for (kk=0; kk < 3; kk++)
    {
      ni[kk] = R[i][0][kk];
      nj[kk] = R[j][0][kk];
    }
  Ci[0] = rx[i];
  Ci[1] = ry[i];
  Ci[2] = rz[i]; 
  Cj[0] = rx[j] + shift[0];
  Cj[1] = ry[j] + shift[1];
  Cj[2] = rz[j] + shift[2]; 
  L = 2.0*typesArr[typeOfPart[i]].sax[0];
  D = 2.0*typesArr[typeOfPart[i]].sax[1];
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }

  for (kk=0; kk < 3; kk++)
    {
      /* center of masses of disks */
      Di[0][kk]=Ci[kk]+0.5*L*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*L*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*L*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*L*nj[kk];
    }

  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==ni[1] && ni[2]==ni[2])
    {
      /* special case of collinear cylinders */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= L)
	    return -1;
	  else
	    return 1;
	}

      /* parallel disks */
      for (j1=0; j1 < 2; j1++)
	for (j2=j1; j2 < 2; j2++)
	  {
	    sp=0.0;
	    for (kk=0; kk < 3; kk++)
	      {
		VV[kk] = Di[j1][kk]-Dj[j2][kk];
		sp += ni[kk]*VV[kk];
	      }
	    if (sp == 0 && calc_norm(VV) < D)
	      return -1;
	  }
    }
  else 
    {
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    vectProdVec(ni, nj, N);
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    vectProdVec(ni,N,niN);
	    vectProdVec(nj,N,njN);
	    normN=calcnorm(N);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/Sqr(normN);
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/Sqr(normN);
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calcnorm(PiDi);
	    normPjDj = calcnorm(PjDj);
	    if (normPiDi <= 0.5*D && normPjDj <= 0.5*D)
	      {
		Q1 = sqrt(Sqr(D)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calcnorm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  return -1;
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  for (j2=0; j2 < 2; j2++)
    {
      for (kk=0; kk < 3; kk++)
	DjCi[kk] = Dj[j2][kk] - Ci[kk];
      normDjCi = calcnorm(DjCi);
      DjCini = scalProd(DjCi,ni);
      for (kk=0; kk < 3; kk++)
	{
	  Ui[kk] = Ci[kk] + DjCini*ni[kk];
	  DjUi[kk] = Dj[j2][kk] - Ui[kk];
	}

      DjUini = scalProd(DjUi,ni);
      normDjUi = calcnorm(DjUi);

      if (normDjUi <= D*0.5 && DjUini > L*0.5)
	return -1;

      for (kk=0; kk < 3; kk++)
	{
	  Ai[kk] = Ci[kk];
	}
      for (it = 0; it < MAX_ITERATIONS; it++);
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      AiDj[kk] = Ai[kk] - Dj[j2][kk];
	    }
	  AiDjnj = scalProd(AiDj,nj);
	  vectProdVec(AiDj,nj,AiDjnjvec);
	  for (kk=0; kk < 3; kk++)
	    VV[kk] =  0.5*D*(Ai[kk]-Dj[j2][kk]-AiDjnj*nj[kk])/calcnorm(AiDjnjvec);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tjp[kk] = Dj[j2][kk] + VV[kk];
	      Tjm[kk] = Dj[j2][kk] - VV[kk];
	      TjpCi[kk] = Tjp[kk] - Ci[kk];
	      TjmCi[kk] = Tjm[kk] - Ci[kk];
	    }
	  TjpCini = scalProd(TjpCi,ni);  
	  TjmCini = scalProd(TjmCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tjp_perp[kk] = TjpCi[kk]-TjpCini*ni[kk];
	      Tjp_para[kk] = TjpCini*ni[kk];
	      Tjm_perp[kk] = TjmCi[kk]-TjmCini*ni[kk];
	      Tjm_para[kk] = TjmCini*ni[kk];
	    } 
	  normTjm_perp = calcnorm(Tjp_perp);
	  for (kk=0; kk < 3; kk++)
	    TjOld[kk] = TjNew[kk];
	  if (calcnorm(Tjm_perp) < calcnorm(Tjp_perp))
	    {
	      for (kk=0; kk < 3; kk++)
		TjNew[kk] = Tjm[kk];
	    }	  
	  else
	    {
	      for (kk=0; kk < 3; kk++)
		TjNew[kk] = Tjp[kk];
	    }

	  for (kk=0; kk < 3; kk++)
	    TjNewCi[kk] = TjNew[kk] - Ci[kk];
	  TjNewCini = scalProd(TjNewCi,ni);
	  Ai[kk] = TjNewCini*ni[kk] + Ci[kk]; 
	  if ( it > 0 && check_convergence(TjOld,TjNew) ) 
	    break;
	} 
      if ( (calcnorm(Tjp_para) <= L*0.5 && calcnorm(Tjp_perp) >= D*0.5)||
	   (calcnorm(Tjm_para) <= L*0.5 && calcnorm(Tjm_perp) >= D*0.5) )
	return -1;
    }

  for (j1=0; j1 < 2; j1++)
    {
      for (kk=0; kk < 3; kk++)
	DiCj[kk] = Di[j1][kk] - Cj[kk];
      normDiCj = calcnorm(DiCj);
      DiCjnj = scalProd(DiCj,nj);
      for (kk=0; kk < 3; kk++)
	{
	  Uj[kk] = Cj[kk] + DiCjnj*nj[kk];
	  DiUj[kk] = Di[j1][kk] - Uj[kk];
	}

      DiUjnj = scalProd(DiUj,nj);
      normDiUj = calcnorm(DiUj);

      if (normDiUj <= D*0.5 && DiUjnj > L*0.5)
	return -1;

      for (kk=0; kk < 3; kk++)
	{
	  Aj[kk] = Cj[kk];
	}
      for (it = 0; it < MAX_ITERATIONS; it++);
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      AjDi[kk] = Aj[kk] - Di[j1][kk];
	    }
	  AjDini = scalProd(AjDi,ni);
	  vectProdVec(AjDi,nj,AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    VV[kk] =  0.5*D*(Aj[kk]-Di[j1][kk]-AjDini*ni[kk])/calcnorm(AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tip[kk] = Di[j2][kk] + VV[kk];
	      Tim[kk] = Di[j2][kk] - VV[kk];
	      TipCj[kk] = Tip[kk] - Cj[kk];
	      TimCj[kk] = Tim[kk] - Cj[kk];
	    }
	  TipCjnj = scalProd(TipCj,nj);  
	  TimCjnj = scalProd(TimCj,nj);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tip_perp[kk] = TipCj[kk]-TipCjnj*nj[kk];
	      Tip_para[kk] = TipCjnj*nj[kk];
	      Tim_perp[kk] = TimCj[kk]-TimCjnj*nj[kk];
	      Tim_para[kk] = TimCjnj*nj[kk];
	    } 
	  normTim_perp = calcnorm(Tip_perp);
	  for (kk=0; kk < 3; kk++)
	    TiOld[kk] = TiNew[kk];
	  if (calcnorm(Tim_perp) < calcnorm(Tip_perp))
	    {
	      for (kk=0; kk < 3; kk++)
		TiNew[kk] = Tim[kk];
	    }	  
	  else
	    {
	      for (kk=0; kk < 3; kk++)
		TiNew[kk] = Tip[kk];
	    }

	  for (kk=0; kk < 3; kk++)
	    TiNewCj[kk] = TiNew[kk] - Cj[kk];
	  TiNewCjnj = scalProd(TiNewCj,nj);
	  Aj[kk] = TiNewCjnj*nj[kk] + Cj[kk]; 
	  if ( it > 0 && check_convergence(TiOld,TiNew) ) 
	    break;
	} 
      if ( (calcnorm(Tip_para) <= L*0.5 && calcnorm(Tip_perp) >= D*0.5)||
	   (calcnorm(Tim_para) <= L*0.5 && calcnorm(Tim_perp) >= D*0.5) )
	return -1;

    }
  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  detA = Sqr(ninj)-1;

  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calcnorm(ViVj) < D && fabs(lambdai) < 0.5*L && fabs(lambdaj) < 0.5*L)
    return -1;

  return 1;
}
#endif
#endif
