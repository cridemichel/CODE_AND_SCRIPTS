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
struct mboxstr 
{
  double sa[3];
  double dr[3];
} **mbox;
double toteneini=0.0;
long long int ttini=0;
int covrestart = 0;
const int nmboxMC=5;
double totdist=0.0, distcc=0.0;
void build_parallelepipeds(void)
{
  double sa[3], dx;
  int tt, kk, k1, k2;

  mbox = malloc(sizeof(struct mboxstr*)*Oparams.ntypes);
  /* 2 multibox per tipo */
  for (tt=0; tt < Oparams.ntypes; tt++)
    mbox[tt] = malloc(sizeof(struct mboxstr)*nmboxMC); 
  for (tt=0; tt < Oparams.ntypes; tt++)
    {
      for (kk = 0; kk < 3; kk++)
	sa[kk] = typesArr[tt].sax[kk]; 
      /* secondo set di parallelepipedi: 
	 ne considero 2 che sono lungo x lunghi quanto il precedente ma che vanno a rimepire
	 i "buchi" lungo y e z. Per ora li costruisco a mano "ad hoc" */
      mbox[tt][0].dr[0]=0.0;
      mbox[tt][0].dr[1]=0.0;
      //mbox[tt][1].dr[2]=(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][0].dr[2]=0.0;
      mbox[tt][0].sa[0]=saxfactMC[0]*sa[0];
      mbox[tt][0].sa[1]=0.48*sa[1]; 
      mbox[tt][0].sa[2]=0.83*sa[2];

      mbox[tt][1].dr[0]=0;
      mbox[tt][1].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][1].dr[2]=0.0;
      mbox[tt][1].sa[0]=saxfactMC[0]*sa[0];
      mbox[tt][1].sa[1]=0.83*sa[1]; 
      mbox[tt][1].sa[2]=0.48*sa[2];

      /* con questo ultimo multibox si evitano configurazioni molto overlapped
	 con asse x quasi parallelo nel caso di grandi elongazioni */
      mbox[tt][2].dr[0]=0;
      mbox[tt][2].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][2].dr[2]=0.0;
      mbox[tt][2].sa[0]=0.94*sa[0];
      mbox[tt][2].sa[1]=0.55*sa[1]; 
      mbox[tt][2].sa[2]=0.55*sa[2];
     
      mbox[tt][3].dr[0]=0;
      mbox[tt][3].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][3].dr[2]=0.0;
      mbox[tt][3].sa[0]=0.99*sa[0];
      mbox[tt][3].sa[1]=0.27*sa[1]; 
      mbox[tt][3].sa[2]=0.27*sa[2];

#if 1
      mbox[tt][4].dr[0]=0;
      mbox[tt][4].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][4].dr[2]=0.0;
      mbox[tt][4].sa[0]=0.999*sa[0];
      mbox[tt][4].sa[1]=0.08*sa[1]; 
      mbox[tt][4].sa[2]=0.08*sa[2];
#endif
    }
}
#endif

#ifdef MC_SIMUL
extern double *overestimate_of_displ;
#ifdef MC_GRANDCAN
int allocnpGC;
#endif
double calcDistBox(int i, int j, double rbi[3], double rbj[3], double saxi[3], double saxj[3], double shift[3])
{
  double RR, R0, R1, cij[3][3], fabscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k, k1, k2, existsParallelPair = 0;
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  for (k=0; k < 3; k++)
    {
      rA[k] = rbi[k];
      rB[k] = rbj[k]+shift[k];
      EA[k] = saxi[k];
      EB[k] = saxj[k];
    }
#if 0
  /* riportare qua anche l'analogo routin per sfere se servirà */
  if (is_a_sphere_NNL[i] && is_a_sphere_NNL[j])
    return calcDistNegNNLoverlapPlaneHS(i, j, rA, rB);
#endif

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  AA[k1][k2] = R[i][k1][k2];
	  BB[k1][k2] = R[j][k1][k2];
	}
    	DD[k1] = rA[k1] - rB[k1];
    }
  /* axis C0+s*A0 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[0][k1] =  scalProd(AA[0], BB[k1]);
      fabscij[0][k1] = fabs(cij[0][k1]);
      if ( fabscij[0][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[0] = scalProd(AA[0],DD);
  RR = fabs(AD[0]);
  R1 = EB[0]*fabscij[0][0]+EB[1]*fabscij[0][1]+EB[2]*fabscij[0][2];
  R01 = EA[0] + R1;
  if ( RR > R01 )
    return 1.0; /* non si intersecano */
  /* axis C0+s*A1 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[1][k1] = scalProd(AA[1],BB[k1]);
      fabscij[1][k1] = fabs(cij[1][k1]);
      if ( fabscij[1][k1] == 1.0  )
	existsParallelPair = 1;
    }
  AD[1] = scalProd(AA[1],DD);
  RR = fabs(AD[1]);
  R1 = EB[0]*fabscij[1][0]+EB[1]*fabscij[1][1]+EB[2]*fabscij[1][2];
  R01 = EA[1] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*A2 */
  for (k1= 0; k1 < 3; k1++)
    {
      cij[2][k1] = scalProd(AA[2], BB[k1]);
      fabscij[2][k1] = fabs(cij[2][k1]);
      if ( fabscij[2][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[2] = scalProd(AA[2],DD);
  RR = fabs(AD[2]);
  R1 = EB[0]*fabscij[2][0]+EB[1]*fabscij[2][1]+EB[2]*fabscij[2][2];
  R01 = EA[2] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*B0 */
  RR = fabs(scalProd(BB[0],DD));
  R0 = EA[0]*fabscij[0][0]+EA[1]*fabscij[1][0]+EA[2]*fabscij[2][0];
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B1 */
  RR = fabs(scalProd(BB[1],DD));
  R0 = EA[0]*fabscij[0][1]+EA[1]*fabscij[1][1]+EA[2]*fabscij[2][1];
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;
  
  /* axis C0+s*B2 */
  RR = fabs(scalProd(BB[2],DD));
  R0 = EA[0]*fabscij[0][2]+EA[1]*fabscij[1][2]+EA[2]*fabscij[2][2];
  R01 = R0 + EB[2];
  if ( RR > R01 )
    return 1.0;

  /* At least one pair of box axes was parallel, therefore the separation is
   * effectively in 2D, i.e. checking the "edge" normals is sufficient for
   * the separation of the boxes. 
   */
  if ( existsParallelPair )
    return -1.0;

  /* axis C0+s*A0xB0 */
  RR = fabs(AD[2]*cij[1][0]-AD[1]*cij[2][0]);
  R0 = EA[1]*fabscij[2][0] + EA[2]*fabscij[1][0];
  R1 = EB[1]*fabscij[0][2] + EB[2]*fabscij[0][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB1 */
  RR = fabs(AD[2]*cij[1][1]-AD[1]*cij[2][1]);
  R0 = EA[1]*fabscij[2][1] + EA[2]*fabscij[1][1];
  R1 = EB[0]*fabscij[0][2] + EB[2]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB2 */
  RR = fabs(AD[2]*cij[1][2]-AD[1]*cij[2][2]);
  R0 = EA[1]*fabscij[2][2] + EA[2]*fabscij[1][2];
  R1 = EB[0]*fabscij[0][1] + EB[1]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB0 */
  RR = fabs(AD[0]*cij[2][0]-AD[2]*cij[0][0]);
  R0 = EA[0]*fabscij[2][0] + EA[2]*fabscij[0][0];
  R1 = EB[1]*fabscij[1][2] + EB[2]*fabscij[1][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB1 */
  RR = fabs(AD[0]*cij[2][1]-AD[2]*cij[0][1]);
  R0 = EA[0]*fabscij[2][1] + EA[2]*fabscij[0][1];
  R1 = EB[0]*fabscij[1][2] + EB[2]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB2 */
  RR = fabs(AD[0]*cij[2][2]-AD[2]*cij[0][2]);
  R0 = EA[0]*fabscij[2][2] + EA[2]*fabscij[0][2];
  R1 = EB[0]*fabscij[1][1] + EB[1]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB0 */
  RR = fabs(AD[1]*cij[0][0]-AD[0]*cij[1][0]);
  R0 = EA[0]*fabscij[1][0] + EA[1]*fabscij[0][0];
  R1 = EB[1]*fabscij[2][2] + EB[2]*fabscij[2][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB1 */
  RR = fabs(AD[1]*cij[0][1]-AD[0]*cij[1][1]);
  R0 = EA[0]*fabscij[1][1] + EA[1]*fabscij[0][1];
  R1 = EB[0]*fabscij[2][2] + EB[2]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB2 */
  RR = fabs(AD[1]*cij[0][2]-AD[0]*cij[1][2]);
  R0 = EA[0]*fabscij[1][2] + EA[1]*fabscij[0][2];
  R1 = EB[0]*fabscij[2][1] + EB[1]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  return -1.0;
}

/* MONTE CARLO CODE START HERE */
double rxold, ryold,rzold, Rold[3][3];
void set_semiaxes_vb_mc(int ii, double fx, double fy, double fz)
{
  nebrTab[ii].axa = fx;//typesArr[typeOfPart[ii]].sax[0]+typesArr[typeOfPart[ii]].spots[0].sigma*0.5;
  nebrTab[ii].axb = fy;//typesArr[typeOfPart[ii]].sax[1];
  nebrTab[ii].axc = fz;//typesArr[typeOfPart[ii]].sax[2];
}

void store_coord(int ip)
{
  int k1, k2;
  rxold=rx[ip];
  ryold=ry[ip];
  rzold=rz[ip];
  for (k1=0; k1<3; k1++)
    for (k2=0; k2<3; k2++)
      Rold[k1][k2]=R[ip][k1][k2]; 
}
void restore_coord(int ip)
{
  int k1, k2;
  rx[ip]=rxold;
  ry[ip]=ryold;
  rz[ip]=rzold;
  for (k1=0; k1<3; k1++)
    for (k2=0; k2<3; k2++)
      R[ip][k1][k2]=Rold[k1][k2]; 
}
void rebuild_nnl_MC(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    overestimate_of_displ[i]=0.0;
  rebuildNNL();
  OprogStatus.lastNNLrebuildMC = Oparams.curStep;
}
extern double ranf(void);
extern void orient(double *ox, double *oy, double *oz);
long long int rotmoveMC=0, tramoveMC=0, totmovesMC=0, trarejMC=0, rotrejMC=0, totrejMC=0, volrejMC=0, volmoveMC=0, excrejMC=0, excmoveMC=0;
double displMC;
void tra_move(int ip)
{
  double dx, dy, dz;
  dx = OprogStatus.deltaMC*(ranf()-0.5);
  dy = OprogStatus.deltaMC*(ranf()-0.5);
  dz = OprogStatus.deltaMC*(ranf()-0.5);
  rx[ip]+= dx;
  ry[ip]+= dy;
  rz[ip]+= dz;
  if (OprogStatus.useNNL)
    {
      displMC = sqrt(Sqr(dx)+Sqr(dy)+Sqr(dz))*1.001;
    }
  tramoveMC++; 
}
extern double calc_norm(double *v);
extern double scalProd(double *, double *);
void remove_parall(int ip, double *ox, double *oy, double *oz)
{
  double sp, vp[3], v[3], nn;
  int k;
  v[0]=*ox;
  v[1]=*oy;
  v[2]=*oz;
  sp = 0;
  for (k=0; k < 3; k++)
   {
     vp[k] = R[ip][0][k];
     sp += vp[k]*v[k];
   }
  
  for (k=0; k < 3; k++)
    {
      v[k] -= vp[k]*sp;
    }
  nn = calc_norm(v);
  for (k=0; k < 3; k++)
    {
      v[k] /= nn;
    } 
  *ox = v[0];
  *oy = v[1];
  *oz = v[2];
  //printf("DOPO o=%f %f %f\n", *ox,*oy,*oz);
  //printf("sp=%f\n", scalProd(v,vp));
}
void rot_move(int ip)
{
  double theta, thetaSq, sinw, cosw;
  double ox, oy, oz, OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random orientation */
  orient(&ox,&oy,&oz);
  remove_parall(ip, &ox, &oy, &oz);
  /* pick a random rotation angle */
  theta= OprogStatus.dthetaMC*(ranf()-0.5);
  if (OprogStatus.useNNL)
    displMC = abs(theta)*0.5001*maxax[ip]; 
  thetaSq=Sqr(theta);
  sinw = sin(theta);
  cosw = (1.0 - cos(theta));
  Omega[0][0] = 0;
  Omega[0][1] = -oz;
  Omega[0][2] = oy;
  Omega[1][0] = oz;
  Omega[1][1] = 0;
  Omega[1][2] = -ox;
  Omega[2][0] = -oy;
  Omega[2][1] = ox;
  Omega[2][2] = 0;
  OmegaSq[0][0] = -Sqr(oy) - Sqr(oz);
  OmegaSq[0][1] = ox*oy;
  OmegaSq[0][2] = ox*oz;
  OmegaSq[1][0] = ox*oy;
  OmegaSq[1][1] = -Sqr(ox) - Sqr(oz);
  OmegaSq[1][2] = oy*oz;
  OmegaSq[2][0] = ox*oz;
  OmegaSq[2][1] = oy*oz;
  OmegaSq[2][2] = -Sqr(ox) - Sqr(oy);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = R[ip][k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += R[ip][k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     R[ip][k1][k2] = Ro[k1][k2]; 
  rotmoveMC++;
}
int random_move(int ip)
{
  double p;
  //printf("random move ip=%d\n", ip);
  p=ranf();
  if (OprogStatus.restrmove==1)
    p=0.0;
  
  if (p <= 0.5)
   {
     tra_move(ip);
     return 0;
   }
  else
    {
      rot_move(ip);
      return 1;
    } 
}
extern double calcDistNegNNLoverlapPlane(double t, double t1, int i, int j, double shift[3]);
extern double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
double overlap_using_multibox(int i, int j, double shift[3])
{
  /* per ora ci sono due livelli di multibox quindi nmb=0 o 1,
     a livello 0 c'è un solo multibox a livello 1 ce ne sono 4 */
  int k1, k2, nmbox, typei, typej, nmboxi, nmboxj;
  double rA[3], rB[3];
  double d0;  

  typei=typeOfPart[i];
  typej=typeOfPart[j];
  nmboxi=nmboxMC;
  nmboxj=nmboxMC;
  /* N.B. non sommo il  mbox[typei][k1].dr[0] poichè so che è zero, se così
     non fosse bisogna ricordarsi di spostare tali assegnazioni nel doppio loop sotto e 
     che il displacemente dr è definito nel sistema
     di riferimento del corpo rigido quindi va trasformato nel laboratorio moltiplicandolo per la
     trasposta di R[i][...][...] (usare la funzione body2labR(...)) */
	
  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  rB[0] = rx[j];
  rB[1] = ry[j];
  rB[2] = rz[j];

  for (k1=0; k1 < nmboxi; k1++)
    {
      for (k2=0; k2 < nmboxj; k2++)
	{
#if 0
	  set_semiaxes_vb_mc(i, mbox[typei][k1].sa[0],
			        mbox[typei][k1].sa[1],
	    		        mbox[typei][k1].sa[2]);
	  set_semiaxes_vb_mc(j, mbox[typej][k2].sa[0],
			        mbox[typej][k2].sa[1],
	    		        mbox[typej][k2].sa[2]);
#endif	  
	  if ((d0=calcDistBox(i, j, rA, rB, mbox[typei][k1].sa, mbox[typej][k2].sa, shift))<0)
	    return d0;
	}
    }
  return 1.0;
}
#undef DEBUG_HCMC
#undef MC_OF_BOXES
double check_overlap(int i, int j, double shift[3], int *errchk)
{
  int k, k1, k2;
  double vecg[8], vecgNeg[8], rA[3], rB[3], saxi[3], saxj[3];
  double d, d0, r1[3], r2[3], alpha; 
  *errchk=0;
  OprogStatus.optnnl = 0;

#ifdef DEBUG_HCMC
  d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
  printf("d=%f\n",d);
  exit(-1);
#endif
  if (are_spheres(i,j))
    {
      OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
      calcdist_retcheck = 0;
      d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
      *errchk = calcdist_retcheck;
      return d;
    }

#if 0
  nebrTab[i].r[0] = rx[i];
  nebrTab[i].r[1] = ry[i];
  nebrTab[i].r[2] = rz[i];
  nebrTab[j].r[0] = rx[j];
  nebrTab[j].r[1] = ry[j];
  nebrTab[j].r[2] = rz[j];
  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  nebrTab[i].R[k1][k2] = R[i][k1][k2];
	  nebrTab[j].R[k1][k2] = R[j][k1][k2];
	}
    }
#endif
#if 0
#if 1
  set_semiaxes_vb_mc(i, 1.01*(typesArr[typeOfPart[i]].sax[0]),
		  1.01*(typesArr[typeOfPart[i]].sax[1]), 
		  1.01*(typesArr[typeOfPart[i]].sax[2]));
  set_semiaxes_vb_mc(j, 1.01*(typesArr[typeOfPart[j]].sax[0]),
		  1.01*(typesArr[typeOfPart[j]].sax[1]), 
		  1.01*(typesArr[typeOfPart[j]].sax[2]));
#else
 set_semiaxes_vb(1.01*(typesArr[typeOfPart[0]].sax[0]),
		  1.01*(typesArr[typeOfPart[0]].sax[1]), 
		  1.01*(typesArr[typeOfPart[0]].sax[2]));
#endif
#endif 
 for (k=0; k < 3; k++)
   {
#ifdef MC_OF_BOXES
     saxi[k] = typesArr[typeOfPart[i]].sax[k];
     saxj[k] = typesArr[typeOfPart[j]].sax[k];
#else
     saxi[k] = 1.01* typesArr[typeOfPart[i]].sax[k];
     saxj[k] = 1.01* typesArr[typeOfPart[j]].sax[k];
#endif
}  
 rA[0] = rx[i];
 rA[1] = ry[i];
 rA[2] = rz[i];
 rB[0] = rx[j];
 rB[1] = ry[j];
 rB[2] = rz[j];
 d0 = calcDistBox(i, j, rA, rB, saxi, saxj, shift);
 //d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
  /* se d0 è positiva vuol dire che i due parallelepipedi non s'intersecano */
  if (d0 > 0.0)
    {
      return 1.0;
    }
#ifdef MC_OF_BOXES
   else
     return -1;
#endif
#if 0
#if 1
  set_semiaxes_vb_mc(i, saxfactMC[0]*typesArr[typeOfPart[i]].sax[0],
		  saxfactMC[1]*typesArr[typeOfPart[i]].sax[1], 
		  saxfactMC[2]*typesArr[typeOfPart[i]].sax[2]);
  set_semiaxes_vb_mc(j, saxfactMC[0]*typesArr[typeOfPart[j]].sax[0],
		  saxfactMC[1]*typesArr[typeOfPart[j]].sax[1], 
		  saxfactMC[2]*typesArr[typeOfPart[j]].sax[2]);
#else
  set_semiaxes_vb(saxfactMC[0]*typesArr[typeOfPart[0]].sax[0],
		  saxfactMC[1]*typesArr[typeOfPart[0]].sax[1], 
		  saxfactMC[2]*typesArr[typeOfPart[0]].sax[2]);
#endif
#endif
  for (k=0; k < 3; k++)
    {
#ifdef MC_QUASI_CUBE
      saxi[k] = saxfactMC_QC[k]*typesArr[typeOfPart[i]].sax[k];
      saxj[k] = saxfactMC_QC[k]*typesArr[typeOfPart[j]].sax[k];
#else
      saxi[k] = saxfactMC[k]*typesArr[typeOfPart[i]].sax[k];
      saxj[k] = saxfactMC[k]*typesArr[typeOfPart[j]].sax[k];
#endif    
}  
 d0 = calcDistBox(i, j, rA, rB, saxi, saxj, shift);
 //d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
  /* se d0 è positiva vuol dire che i due parallelepipedi non s'intersecano */
  if (d0 < 0.0)
    {
      return -1.0;
    }

#if 0
  d0=overlap_using_multibox(i, j, shift);
  if (d0 < 0)
    {
      return d0;
    }
#endif
  OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
  calcdist_retcheck = 0;
  
  d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
  *errchk = calcdist_retcheck;
#ifndef MC_QUASI_CUBE
  if (*errchk)
    {
      d0=overlap_using_multibox(i, j, shift);
      if (d0 > 0)
	printf("I used multibox routine and d0=%f\n", d0);
      //printf("errcheck=%d\n", *errchk);
      if (d0 < 0)
	{
	  *errchk=0;
	  return d0;
	}
      return -1.0;
    }
#endif
  //printf("QUI d=%f\n", d);
  return d;
}
extern void find_bonds_one_NLL(int i);
#ifdef MC_QUASI_CUBE
double calcdistsaQC(double ra[3], double rb[3], double u1a[3], double u2a[3], double u1b[3], double u2b[3], double shift[3], int *errchk)
{
  /* N.B. u1, u2 ed u3 sono i vettori del sistema di riferimento solidale con il corpo rigido 
     espressi nel riferminto del laboratorio */
  int k1, k2;
  double nn, vecg[8], vecgNeg[8], Rla[3][3], Rlb[3][3];
  double d, r1[3], r2[3], alpha, d0;
  rx[0] = ra[0];
  ry[0] = ra[1];
  rz[0] = ra[2];
  rx[1] = rb[0];
  ry[1] = rb[1];
  rz[1] = rb[2];
  for (k1 = 0; k1 < 3; k1++)
    {
      Rla[0][k1] = u1a[k1];
      Rlb[0][k1] = u1b[k1];
      Rla[1][k1] = u2b[k1];
      Rlb[1][k1] = u2b[k1];
    }
  
  vectProdVec(u1a, u2a, Rla[2]);
  vectProdVec(u1b, u2b, Rlb[2]);
  
  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  R[0][k1][k2] = Rla[k1][k2];
	  R[1][k1][k2] = Rlb[k1][k2];
	}
    }

  OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
  calcdist_retcheck = 0;
  d = check_overlap(0, 1, shift, errchk);
  return d;
}
#endif
int overlapMC_NNL(int na, int *err)
{
  int i, signDir[NDIM]={0,0,0}, evCode, k, n;
  double vecg[5], shift[3], t1, t2, t, tm[NDIM];
  double sigSq, tInt, d, b, vv, dv[3], dr[3], distSq;
  int overlap;
#ifdef MD_PATCHY_HE
  int nb, ac, bc, collCode, collCodeOld, acHC, bcHC;
  double evtime, evtimeHC;
#ifdef MD_EDHEFLEX_WALL
  int nplane=-1;
#endif
#endif
#ifdef MD_ABSORPTION
  int hwcell;
#endif
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
      PredictEventNNL_MLL(na, nb);
      return;
    }
#endif
  nb=-1;
  for (i=0; i < nebrTab[na].len; i++)
    {
      n = nebrTab[na].list[i]; 
      if (!(n != na && n!=nb && (nb >= -1 || n < na)))
	continue;
#ifdef MD_LXYZ
      shift[0] = L[0]*rint((rx[na]-rx[n])/L[0]);
      shift[1] = L[1]*rint((ry[na]-ry[n])/L[1]);
#ifdef MD_EDHEFLEX_WALL
      if (!OprogStatus.hardwall)
	shift[2] = L[2]*rint((rz[na]-rz[n])/L[2]);
      else
	shift[2] = 0.0;
#else
      shift[2] = L[2]*rint((rz[na]-rz[n])/L[2]);
#endif
#else
      shift[0] = L*rint((rx[na]-rx[n])/L);
      shift[1] = L*rint((ry[na]-ry[n])/L);
#ifdef MD_EDHEFLEX_WALL
      if (!OprogStatus.hardwall)
	shift[2] = L*rint((rz[na]-rz[n])/L);
      else
	shift[2] = 0.0;
#else
      shift[2] = L*rint((rz[na]-rz[n])/L);
#endif
#endif
      if (check_overlap(na, n, shift, err)<0.0)
	{
	  //printf("checking i=%d j=%d: ", na, n);
	  //printf("overlap!\n");
	  return 1;
	}
    }
  return 0;
}
int overlapMC_LL(int ip, int *err)
{
  int kk, nb, k, iZ, jZ, iX, jX, iY, jY, n, na;
  int cellRangeT[6];
  double shift[3];
  na=ip;
  nb=-1;
  for (kk = 0;  kk < 3; kk++)
    {
      cellRange[2*kk]   = - 1;
      cellRange[2*kk+1] =   1;
    }

  for (k = 0; k < 6; k++) cellRangeT[k] = cellRange[k];

  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
      if (jZ == -1) 
	{
	  jZ = cellsz - 1;    
#ifdef MD_LXYZ
	  shift[2] = - L[2];
#else
	  shift[2] = - L;
#endif
	} 
      else if (jZ == cellsz) 
	{
	  jZ = 0;    
#ifdef MD_LXYZ
	  shift[2] = L[2];
#else
	  shift[2] = L;
#endif
	}
      
      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	{
	  jY = inCell[1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cellsy - 1;    
#ifdef MD_LXYZ
	      shift[1] = -L[1];
#else
	      shift[1] = -L;
#endif
	    } 
	  else if (jY == cellsy) 
	    {
	      jY = 0;    
#ifdef MD_LXYZ
	      shift[1] = L[1];
#else
	      shift[1] = L;
#endif
	    }
	  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
	    {
	      jX = inCell[0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cellsx - 1;    
#ifdef MD_LXYZ
		  shift[0] = - L[0];
#else
		  shift[0] = - L;
#endif
		} 
	      else if (jX == cellsx) 
		{
		  jX = 0;   
#ifdef MD_LXYZ
		  shift[0] = L[0];
#else
		  shift[0] = L;
#endif
		}
	      n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	      for (n = cellList[n]; n > -1; n = cellList[n]) 
		{
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
#if defined(EDHE_FLEX) && 0
		      if (!may_interact_all(na, n))
			continue;
#endif

		      if (check_overlap(na, n, shift, err)<0.0)
			{
			  //printf("checking i=%d j=%d: ", na, n);
			  //printf("overlap!\n");
		  	  return 1;
			}
		    }
		} 
	    }
	}
    }
  return 0;
}
int overlapMC(int ip, int *err)
{
  *err=0;
  if (OprogStatus.useNNL)
    return overlapMC_NNL(ip, err);
  else
    return overlapMC_LL(ip, err);
}
void remove_from_current_cell(int i)
{
  int n;
  n = (inCell[2][i] * cellsy + inCell[1][i])*cellsx + inCell[0][i]
    + Oparams.parnum;
  
  while (cellList[n] != i) 
    n = cellList[n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellList[n] = cellList[i];
}

void insert_in_new_cell(int i)
{
  int n;
  n = (inCell[2][i] * cellsy + inCell[1][i])*cellsx + 
    inCell[0][i] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellList[i] = cellList[n];
  cellList[n] = i;
}
void update_LL(int n)
{
  int cox, coy, coz, cx, cy, cz;

  cox=inCell[0][n];
  coy=inCell[1][n];
  coz=inCell[2][n];
#ifdef MD_LXYZ
  cx =  (rx[n] + 0.5*L[0]) * cellsx / L[0];
  cy =  (ry[n] + 0.5*L[1]) * cellsy / L[1];
  cz =  (rz[n] + 0.5*L[2]) * cellsz / L[2];
#else
  cx =  (rx[n] + 0.5*L) * cellsx / L;
  cy =  (ry[n] + 0.5*L) * cellsy / L;
#ifdef MD_GRAVITY
  cz =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
  cz =  (rz[n] + 0.5*L)  * cellsz / L;
#endif
#endif
  if (cx!=cox || cy!=coy || cz!=coz)
    {
      remove_from_current_cell(n);
      inCell[0][n] = cx;
      inCell[1][n] = cy;
      inCell[2][n] = cz;
      insert_in_new_cell(n);
    }
}

void pbc(int ip)
{
  double L2[3], Ll[3];
  
#ifdef MD_LXYZ
  L2[0] = L[0]*0.5;
  L2[1] = L[1]*0.5;
  L2[2] = L[2]*0.5;
  Ll[0] = L[0];
  Ll[1] = L[1];
  Ll[2] = L[2];
#else
  L2[0] = L2[1] = L2[2] = L*0.5;
  Ll[0] = Ll[1] = Ll[2] = L;
#endif
  if (rx[ip] > L2[0])
    {
      rx[ip] -= L[0];
    }
  else if (rx[ip] < -L2[0])
    {    
      rx[ip] += L[0];
    }
  if (ry[ip] > L2[1])
    {
      ry[ip] -= L[1];
    }
  else if (ry[ip] < -L2[1])
    {
      ry[ip] += L[1];
    }
  if (rz[ip] > L2[2])
    {
      rz[ip] -= L[2];
    }
  else if (rz[ip] < -L2[2])
    {
      rz[ip] += L[2];
    }
}
void update_bonds_MC(int ip);
extern void find_bonds_flex_all(void);
extern double calcpotene(void);
extern void rebuildLinkedList(void);
int cxini=-1, cyini=-1, czini=-1;
void update_numcells(void)
{
#ifdef MD_LXYZ
  cellsx = L[0] / Oparams.rcut;
  cellsy = L[1] / Oparams.rcut;
  cellsz = L[2] / Oparams.rcut;
#else
  cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
  cellsz = L / Oparams.rcut;
#endif
  if (cellsx > cxini || cellsy > cyini || cellsz > czini)
    {
      cellList = realloc(cellList, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
    } 
}
void find_bonds_flex_NNL(void);
#ifdef MC_GRANDCAN
extern int *listtmp;
void remove_from_nnl_MC(int ip)
{ 
  int k1, k2, len, p;

  for (k1=0; k1 < nebrTab[ip].len; k1++)
    {
      p = nebrTab[ip].list[k1];
      len = nebrTab[p].len;
      memcpy(listtmp, nebrTab[p].list, nebrTab[p].len*sizeof(int));
      nebrTab[p].len=0;
      for (k2=0; k2 < len; k2++)
	{
	  if (nebrTab[p].list[k2]!=ip)
	    {
	      nebrTab[p].list[nebrTab[p].len] = listtmp[k2];
	      nebrTab[p].len++;
	    }
	}
    }
}
void remove_bonds_GC(int ip);
#ifdef MCGC_OPTLLREBUILD
void adjLinkedListRemove(void)
{
  int k;
  for (k = Oparams.parnum; k < cellsx*cellsy*cellsz + Oparams.parnum; k++)
    {
      cellList[k] = cellList[k+1];
    }

}
void adjLinkedListInsert(void)
{
  int k;
  for (k = cellsx*cellsy*cellsz + Oparams.parnum - 1; k >= Oparams.parnum; k--)
    {
      cellList[k] = cellList[k-1];
    }
}
#endif
void update_bonds_GC(int oldi, int newi)
{
  int kk;
#ifdef MD_LL_BONDS
  int i, nb;
  long long int aa, bb, ii, jj, jj2;
#else
  int i, nb, ii, jj, aa, bb, jj2;
#endif
 
  for (i=0; i < Oparams.parnum; i++)
    {
      nb = numbonds[i];
      if (i==newi)
	continue;
      for (kk=0; kk < nb; kk++)
	{
	  jj = bonds[i][kk] / (NANA);
	  jj2 = bonds[i][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (jj==oldi)
	    {
	      //printf("QUIIIIII\n");
#ifdef MD_LL_BONDS
	      bonds[i][kk] = newi*(((long long int)NA)*NA)+aa*((long long int)NA)+bb;
#else
	      bonds[i][kk] = newi*(NANA)+aa*NA+bb;
#endif
	    }
	}
    }

}
void remove_par_GC(int ip)
{
  int lp, k, k1, k2;

  lp = Oparams.parnum-1;
  
  if (lp==ip)
    {
      //printf("1) [#%d] i=%d Oparams.parnum=%d\n QUI\n", Oparams.curStep, ip,  Oparams.parnum);
      remove_bonds_GC(ip);

#ifdef MCGC_OPTLLREBUILD
      remove_from_current_cell(ip);
#endif
      if (OprogStatus.useNNL)
	{
	  /* rimuove ip da tutte le NNL del sistema */
	  remove_from_nnl_MC(ip);
	}
      Oparams.parnum--;
      typeNP[0]--;
#ifdef MCGC_OPTLLREBUILD
      /* N.B. essendo cambiato Oparams.parnum le linked list vanno aggiustate (poichè
	 le linked lists head stanno da Oparams.parnum compreso in poi quindi vanno shiftate, 
	 in questo caso vanno spostate di uno a sinistra) */
      adjLinkedListRemove();
#else
      rebuildLinkedList(); 
#endif
      return;
    }

  //printf("2) [#%d] i=%d Oparams.parnum=%d\n QUI\n", Oparams.curStep, ip,  Oparams.parnum);
  /* copy all data from particle lp to ip */
  remove_bonds_GC(ip);
  //remove_from_current_cell(ip);
  rx[ip] = rx[lp];
  ry[ip] = ry[lp];
  rz[ip] = rz[lp];
  typeOfPart[ip] = typeOfPart[lp];
  is_a_sphere_NNL[ip] = is_a_sphere_NNL[lp];
  numbonds[ip] = numbonds[lp];
#ifdef MD_LL_BONDS
  memcpy( (void*) bonds[ip], (void*) bonds[lp], sizeof(long long int)*numbonds[lp]);
#else
  memcpy( (void*) bonds[ip], (void*) bonds[lp], sizeof(int)*numbonds[lp]);
#endif
  axa[ip] = axa[lp];
  axb[ip] = axb[lp];
  axc[ip] = axb[lp];
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	R[ip][k1][k2] = R[lp][k1][k2];
      }
  if (OprogStatus.useNNL)
    {
      /* rimuove ip da tutte le NNL del sistema */
      remove_from_nnl_MC(ip);
      nebrTab[ip].len = nebrTab[lp].len;
      memcpy(nebrTab[ip].list, nebrTab[lp].list, sizeof(int)*OprogStatus.nebrTabFac);
    }
  maxax[ip] = maxax[lp]; 
#ifdef MCGC_OPTLLREBUILD
  remove_from_current_cell(ip);
  remove_from_current_cell(lp);
#endif
  for (k=0; k < 3; k++)
    inCell[k][ip] = inCell[k][lp];
  typeNP[0]--;
  Oparams.parnum--;
  update_bonds_GC(lp, ip);
#ifdef MCGC_OPTLLREBUILD
  /* N.B. essendo cambiato Oparams.parnum le linked list vanno aggiustate (poichè
     le linked lists head stanno da Oparams.parnum compreso in poi quindi vanno shiftate, 
     in questo caso vanno spostate di uno a sinistra) */
  adjLinkedListRemove();
  insert_in_new_cell(ip);
#else
  rebuildLinkedList(); 
#endif
}
int dyn_realloc_oprog(int np)
{
  int i;  
  void *last_ptr;
#ifdef MD_CALC_DPP
  OprogStatus.len = sizeof(double)*22*np;
#else
  OprogStatus.len = sizeof(double)*10*np;
#endif
  //printf("DYNAMIC ALLOCATION of %d bytes\n", OprogStatus.len);
  OprogStatus.ptr = realloc(OprogStatus.ptr, OprogStatus.len);
  last_ptr = OprogStatus.ptr;
#ifdef MD_CALC_DPP
  OprogStatus.sumdx = (double*)last_ptr;
  OprogStatus.sumdy = OprogStatus.sumdx + np;
  OprogStatus.sumdz = OprogStatus.sumdy + np;
  OprogStatus.lastu1x = OprogStatus.sumdz + np;
  OprogStatus.lastu1y = OprogStatus.lastu1x + np;
  OprogStatus.lastu1z = OprogStatus.lastu1y + np;
  OprogStatus.lastu2x = OprogStatus.lastu1z + np;
  OprogStatus.lastu2y = OprogStatus.lastu2x + np;
  OprogStatus.lastu2z = OprogStatus.lastu2y + np;
  OprogStatus.lastu3x = OprogStatus.lastu2z + np;
  OprogStatus.lastu3y = OprogStatus.lastu3x + np;
  OprogStatus.lastu3z = OprogStatus.lastu3y + np;
  last_ptr = (void*) (OprogStatus.lastu3z + np);
#endif
  OprogStatus.sumox = (double*)last_ptr;
  OprogStatus.sumoy = OprogStatus.sumox + np;
  OprogStatus.sumoz = OprogStatus.sumoy + np;
  OprogStatus.lastcolltime = OprogStatus.sumoz + np;
  OprogStatus.rxCMi = OprogStatus.lastcolltime + np;
  OprogStatus.ryCMi = OprogStatus.rxCMi + np;
  OprogStatus.rzCMi = OprogStatus.ryCMi + np;
  OprogStatus.DR = realloc(OprogStatus.DR, sizeof(double*)*np);
  for (i=0; i < np; i++)
    {
      OprogStatus.DR[i] = OprogStatus.rzCMi + np + i*3;
    }
#if 0
  OprogStatus.vcmx0 = OprogStatus.DR[np-1] + 3;
  OprogStatus.vcmy0 = OprogStatus.vcmx0 + np;
  OprogStatus.vcmz0 = OprogStatus.vcmy0 + np;
#endif
  OprogStatus.set_dyn_ascii();
  return OprogStatus.len;
}
extern double **matrix(int n, int m);
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);
extern double *max_step_MC;
void check_alloc_GC(void)
{
  int size, allocnpGCold, k, i;
  double *rt[3];
  if (Oparams.parnum+1 > allocnpGC)
    {
      allocnpGCold=allocnpGC;
      allocnpGC = (int) (1.2*((double)allocnpGC));
      if (allocnpGC==allocnpGCold)
	allocnpGC *= 2;
      printf("new allocnpGC=%d old=%d\n", allocnpGC, allocnpGCold);
      rt[0] = malloc(sizeof(double)*Oparams.parnum);
      rt[1] = malloc(sizeof(double)*Oparams.parnum);
      rt[2] = malloc(sizeof(double)*Oparams.parnum);
      is_a_sphere_NNL = realloc(is_a_sphere_NNL, sizeof(int)*allocnpGC);
      for (i=0; i < Oparams.parnum; i++)
	{
	  rt[0][i] = rx[i];
	  rt[1][i] = ry[i];
	  rt[2][i] = rz[i];
	}	  
      free(rx);
      AllocCoord(allocnpGC*sizeof(double), ALLOC_LIST, NULL); 
      for (i=0; i < Oparams.parnum; i++)
	{
	  rx[i] = rt[0][i];
	  ry[i] = rt[1][i];
	  rz[i] = rt[2][i];
	}
      numbonds = realloc(numbonds,allocnpGC*sizeof(int));
#ifdef MC_STOREBONDS
      numbondsMC = realloc(numbondsMC,allocnpGC*sizeof(int));
#endif
#ifdef MD_LL_BONDS
      bonds = realloc(bonds, allocnpGC*sizeof(long long int*));
#ifdef MC_STOREBONDS
      bondsMC = realloc(bondsMC, allocnpGC*sizeof(long long int*));
#endif
      size = sizeof(long long int);
#else
      bonds = realloc(bonds, allocnpGC*sizeof(int *));
#ifdef MC_STOREBONDS
      bondsMC = realloc(bondsMC, allocnpGC*sizeof(int *));
#endif
      size = sizeof(int);
#endif
      for (k=allocnpGCold; k < allocnpGC; k++)
	{
  	  bonds[k] = malloc(size*OprogStatus.maxbonds);
#ifdef MC_STOREBONDS
	  bondsMC[k] = malloc(size*OprogStatus.maxbonds);
#endif
	}
      typeOfPart = realloc(typeOfPart, sizeof(int)*allocnpGC);
      R = realloc(R,sizeof(double**)*allocnpGC);
#ifdef MD_MATRIX_CONTIGOUS
      //printf("prima: "); 
      //print_matrix(R[10], 3);

      R[0] = malloc(sizeof(double*)*3);
      R[0][0] = realloc(R[k][0], sizeof(double)*allocnpGC*9);
      R[0][1] = R[0][0] + 3;
      R[0][2] = R[0][1] + 3;
      for (k=1; k < allocnpGC; k++)
	{
	  R[k] = malloc(sizeof(double*)*3);
	  R[k][0] = R[k-1][2] + 3;
	  R[k][1] = R[k][0] + 3;
	  R[k][2] = R[k][1] + 3;
	}
#else
      for (k=allocnpGCold; k < allocnpGC; k++)
	{
    	  R[k] = matrix(3, 3);
	}
#endif
      //printf("dopo: "); 
      //print_matrix(R[10], 3);
      maxax = realloc(maxax, allocnpGC*sizeof(double));
      axa   = realloc(axa,   allocnpGC*sizeof(double));
      axb   = realloc(axb,   allocnpGC*sizeof(double));
      axc   = realloc(axc,   allocnpGC*sizeof(double));
      for (k=allocnpGCold; k < allocnpGC; k++)
	{
	  /* per ora presuppone tutte le particelle uguali */
	  maxax[k] = maxax[0];
	}
      if (OprogStatus.useNNL)
	{
	  nebrTab = realloc(nebrTab, sizeof(struct nebrTabStruct)*allocnpGC);
	  for (k=allocnpGCold; k <  allocnpGC; k++)
	    {
	      nebrTab[k].len = 0;
	      nebrTab[k].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	      //nebrTab[i].shift = matrix(OprogStatus.nebrTabFac, 3);
	      /* tanto shift non viene usata per cui è inutile sprecare memoria*/ 
	      nebrTab[k].shift = NULL;
	      nebrTab[k].R = matrix(3, 3);
	    }
	  overestimate_of_displ = realloc(overestimate_of_displ, sizeof(double)*allocnpGC);
	  max_step_MC = realloc(max_step_MC, sizeof(double)*allocnpGC);
	}
      cellList = realloc(cellList, sizeof(int)*(cellsx*cellsy*cellsz+allocnpGC));
      inCell[0] = realloc(inCell[0],sizeof(int)*allocnpGC);
      inCell[1] = realloc(inCell[1],sizeof(int)*allocnpGC);
      inCell[2] = realloc(inCell[2],sizeof(int)*allocnpGC);
      dyn_realloc_oprog(allocnpGC);
    }
}
extern void remove_bond(int na, int n, int a, int b);

void remove_bonds_GC(int ip)
{
  int kk;
#ifdef MD_LL_BONDS
  int i, nb;
  long long int aa, bb, ii, jj, jj2;
#else
  int i, nb, ii, jj, aa, bb, jj2;
#endif
  nb = numbonds[ip];
  for (kk=0; kk < nb; kk++)
    {
      jj = bonds[ip][kk] / (NANA);
      jj2 = bonds[ip][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      remove_bond(jj, ip, bb, aa);
    }
}
extern void BuildNNL(int na);
extern void nextNNLupdate(int na);

void build_one_nnl_GC(int ip)
{
  int p, k1;
  nextNNLupdate(ip);
  BuildNNL(ip);
  /* and now add ip to neighbor lists of existing particles */
  for (k1=0; k1 < nebrTab[ip].len; k1++)
    {
      p=nebrTab[ip].list[k1];
      nebrTab[p].list[nebrTab[p].len] = ip;
      nebrTab[p].len++;
    }
}
extern void versor_to_R(double ox, double oy, double oz, double R[3][3]);
extern void find_bonds_one(int i);
double calc_maxstep_MC(int i);
#ifdef MCGC_OPTLLREBUILD
void assign_cell_GC(int np)
{
#ifdef MD_LXYZ
  inCell[0][np] =  (rx[np] + L2[0]) * cellsx / L[0];
  inCell[1][np] =  (ry[np] + L2[1]) * cellsy / L[1];
  inCell[2][np] =  (rz[np] + L2[2]) * cellsz / L[2];
#else
  inCell[0][np] =  (rx[np] + L2) * cellsx / L;
  inCell[1][np] =  (ry[np] + L2) * cellsy / L;
#ifdef MD_GRAVITY
  inCell[2][np] =  (rz[np] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
  inCell[2][np] =  (rz[np] + L2)  * cellsz / L;
#endif
#endif
}
#endif
int insert_particle_GC(void)
{
  int np, k1, k2;
  double ox, oy, oz, Rl[3][3];
  np = Oparams.parnum;
  /* controlla se c'è abbastanza allocazione se no espandere
     tutto ciò che va espanso */
  check_alloc_GC();
  /* choose a random position to insert particle */
#ifdef MD_LXYZ
  rx[np] = L[0]*(ranf()-0.5);
  ry[np] = L[1]*(ranf()-0.5);
  rz[np] = L[2]*(ranf()-0.5);
#else
  rx[np] = L*(ranf()-0.5);
  ry[np] = L*(ranf()-0.5);
  rz[np] = L*(ranf()-0.5);
#endif
  if (OprogStatus.restrmove==0)
    {
      orient(&ox,&oy,&oz);
      versor_to_R(ox, oy, oz, Rl);
    }
  else
    {
      Rl[0][0]=Rl[1][1]=Rl[2][2]=1.0;
      Rl[0][1]=Rl[0][2]=Rl[1][0]=Rl[1][2]=Rl[2][0]=Rl[2][1]=0.0;
    }
  /* per ora assumiamo un solo tipo di particelle nel GC */
  numbonds[np] = 0;
  typeOfPart[np]=0;
  vx[np]=vy[np]=vz[np]=wx[np]=wy[np]=wz[np]=Mx[np]=My[np]=Mz[np]=0.0;
  is_a_sphere_NNL[np] = is_a_sphere_NNL[0]; 
  maxax[np] = maxax[0];
  axa[np]=axa[0];
  axb[np]=axb[0];
  axc[np]=axc[0];
  //printf("np=%d v=%f %f %f\n", np, vx[np], vy[np], vz[np]);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	R[np][k1][k2] = Rl[k1][k2];
      }
  //update_LL(np);
  typeNP[0]++;
  Oparams.parnum++;
#ifdef MCGC_OPTLLREBUILD
  /* N.B. essendo cambiato Oparams.parnum le linked list vanno aggiustate (poichè
     le linked lists head stanno da Oparams.parnum compreso in poi  vanno shiftate, 
     in questo caso vanno spostate di uno a destra) */
  adjLinkedListInsert();
  assign_cell_GC(np);
  insert_in_new_cell(np);
#else
  rebuildLinkedList();
#endif
  if (OprogStatus.useNNL)
    {
      build_one_nnl_GC(np);
      overestimate_of_displ[np]=0.0;
      max_step_MC[np] = calc_maxstep_MC(np);
    }
  return np;
}
void find_bonds_GC(int ip);
extern int is_in_ranges(int A, int B, int nr, rangeStruct* r);
double calcpotene_GC(int ip)
{
  double Epot; 
  int na;
#ifdef EDHE_FLEX
#ifdef MD_LL_BONDS
  long long int jj, jj2, aa, bb;
  int kk, kk2;
#else
  int kk2, jj, kk, jj2, aa, bb;
#endif
#endif
  Epot = 0;
#ifdef EDHE_FLEX
  na=ip;
  for (kk=0; kk < numbonds[na]; kk++)
    {
      jj = bonds[na][kk]/(NANA);
      jj2 = bonds[na][kk]%(NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      //printf("numbonds[%d]=%d aa=%lld bb=%lld ninters:%d\n", na, numbonds[na], aa, bb, Oparams.ninters);
      for (kk2 = 0; kk2 < Oparams.ninters; kk2++)
	{
	  //printf("type[%d]=%d type[%lld]=%d \n", na, typeOfPart[na], jj, typeOfPart[jj]);
	  if ( (is_in_ranges(typeOfPart[na], intersArr[kk2].type1, intersArr[kk2].nr1, intersArr[kk2].r1) && 
		is_in_ranges(typeOfPart[jj], intersArr[kk2].type2, intersArr[kk2].nr2, intersArr[kk2].r2) &&
		intersArr[kk2].spot1 == aa-1 && intersArr[kk2].spot2 == bb-1) || 
	       (is_in_ranges(typeOfPart[jj], intersArr[kk2].type1, intersArr[kk2].nr1, intersArr[kk2].r1) && 
		is_in_ranges(typeOfPart[na], intersArr[kk2].type2, intersArr[kk2].nr2, intersArr[kk2].r2) &&
		intersArr[kk2].spot1 == bb-1 && intersArr[kk2].spot2 == aa-1) )  
	    {
	      Epot -= intersArr[kk2].bheight;
#ifdef MD_SPHERICAL_WALL
	      /* 14/07/08 NOTA: i muri sferici hanno sempre numbonds[sphWall/sphWallOuter]=0 quindi
		 si considerano doppie le interazioni tra una certa particella ed un muro sferico,
		 poiché il potenziale viene calcolato considerando che per ogni interazione ce n'è una simmetrica */ 
	      if (jj==sphWall || jj==sphWallOuter)
		Epot -= intersArr[kk2].bheight;
#endif

	    }		 
	}
#ifdef MD_SWDUMBBELL
      Epot+=swdb_adjust_Epot(na);
#endif
      if (Oparams.nintersIJ > 0)
	{
	  for (kk2 = 0; kk2 < Oparams.nintersIJ; kk2++)
	    {
	      if ( (is_in_ranges(na, intersArrIJ[kk2].i, intersArrIJ[kk2].nr1, intersArrIJ[kk2].r1) && 
		    is_in_ranges(jj, intersArrIJ[kk2].j, intersArrIJ[kk2].nr2, intersArrIJ[kk2].r2) &&
		    intersArrIJ[kk2].spot1 == aa-1 && intersArrIJ[kk2].spot2 == bb-1) || 
		   (is_in_ranges(jj, intersArrIJ[kk2].i, intersArrIJ[kk2].nr1, intersArrIJ[kk2].r1) && 
		    is_in_ranges(na, intersArrIJ[kk2].j, intersArrIJ[kk2].nr2, intersArrIJ[kk2].r2) &&
		    intersArrIJ[kk2].spot1 == bb-1 && intersArrIJ[kk2].spot2 == aa-1) )  
		{
		  Epot -= intersArrIJ[kk2].bheight;
		}
	    }	
	}
    }
  //Epot -= numbonds[na];
#else
  Epot -= numbonds[na];
#endif
  return Epot;
}

void mcexc(int *ierr)
{
  int o, np;
  double enn, eno, arg, vol;
  double xn[3]; 
#ifdef MD_LXYZ
  vol = L[0]*L[1]*L[2];
#else
  vol = L*L*L;
#endif
  //printf("INIZIO MCEXC\n"); 
  /* grand canonical ensemble */
  if (ranf() < 0.5)
    {
      if (Oparams.parnum==0)
	{
#ifdef MC_SUS	  
	  if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
	    OprogStatus.sushisto[Oparams.parnum-OprogStatus.susnmin]++;
#endif
	  return;
	}
      o = Oparams.parnum*ranf();
      eno=calcpotene_GC(o);
      arg = Oparams.parnum*exp((1.0/Oparams.T)*eno)/(OprogStatus.zetaMC*vol);
      //printf("arg=%.15G zetaMC=%.15G \n", arg, OprogStatus.zetaMC);
      if (ranf() < arg)
	{
	  //printf("removing #%d\n", o);
#ifdef MC_SUS	  
	  if ((OprogStatus.susnmin < 0 || OprogStatus.susnmax <= 0) || Oparams.parnum > OprogStatus.susnmin)
	    {
	      remove_par_GC(o);
	    }

#else
	  remove_par_GC(o);
#endif
	}
#ifdef MC_SUS
      if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
	OprogStatus.sushisto[Oparams.parnum-OprogStatus.susnmin]++;
#endif
      //printf("FINE MCEXC\n"); 
    }
  else
    {
      /* nella seguente routine deve aggiungere la particella nelle LL e nelle NNL se usate */
      //printf("Inserting #%d\n", Oparams.parnum);
      np=insert_particle_GC();
      //printf("FINE-2 MCEXC\n"); 
      if (overlapMC(np, ierr))
	{
	  /* reject insertion */
	  //remove_from_current_cell(Oparams.parnum-1);
	  //printf("overlap Insertion rejected #%d\n", Oparams.parnum-1);
#ifdef MCGC_OPTLLREBUILD
	  remove_from_current_cell(np);
#endif
	  Oparams.parnum--;
	  typeNP[0]--;
#ifdef MCGC_OPTLLREBUILD
	  //printf("Celllist[parnum]=%d\n", cellList[Oparams.parnum+1]);
	  adjLinkedListRemove();
	  //printf("Celllist[parnum+1]=%d\n", cellList[Oparams.parnum]);
#else
	  rebuildLinkedList();
#endif
	  //printf("FINE-3 MCEXC\n"); 
	  /* rimuove ip da tutte le NNL del sistema */
	  if (OprogStatus.useNNL)
	    remove_from_nnl_MC(np);
	  excrejMC++;
#ifdef MC_SUS	  
	  if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
	    OprogStatus.sushisto[Oparams.parnum-OprogStatus.susnmin]++;
#endif
	  return;
	}	
      find_bonds_GC(np);
      enn=calcpotene_GC(np);
      //printf("enn=%.15G\n", enn);
      arg = OprogStatus.zetaMC*vol*exp(-(1.0/Oparams.T)*enn)/Oparams.parnum;

      if (ranf() >= arg
#ifdef MC_SUS
	  || (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0 && Oparams.parnum > OprogStatus.susnmax)
#endif
	  )
	{
	  //printf("Insertion rejected #%d\n", np);
	  /* insertion rejected */
	  //remove_from_current_cell(Oparams.parnum-1);
	  remove_bonds_GC(np);
#ifdef MCGC_OPTLLREBUILD
	  remove_from_current_cell(np);
#endif
	  Oparams.parnum--;
	  typeNP[0]--;
#ifdef MCGC_OPTLLREBUILD
	  adjLinkedListRemove();
#else
	  rebuildLinkedList();
#endif
	  /* rimuove ip da tutte le NNL del sistema */
	  if (OprogStatus.useNNL)
	    remove_from_nnl_MC(np);
	  excrejMC++;
#ifdef MC_SUS	  
	  if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
	    OprogStatus.sushisto[Oparams.parnum-OprogStatus.susnmin]++;
#endif
	  return;
	}
#ifdef MC_SUS
      if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
	OprogStatus.sushisto[Oparams.parnum-OprogStatus.susnmin]++;
#endif
    }
  //printf("MCEXC FINE\n");
}
#endif
#ifdef MC_STORELL
/* queste routine vengono eventualmente usate (#ifdef MC_STORELL)
   in move_box() per il ripristino delle LL se la mossa di volume viene rifiutata */
void store_ll_mc(void)
{
  int k;
  for (k=0; k < Oparams.parnum + cellsx*cellsy*cellsz; k++)
    {
      cellListMC[k] = cellList[k];
    }
}
void restore_ll_mc(void)
{
  int k;
  for (k=0; k < Oparams.parnum + cellsx*cellsy*cellsz; k++)
    {
      cellList[k] = cellListMC[k];
    }
}
#endif
#ifdef MC_STOREBONDS
void store_bonds_mc(int ip)
{
#ifdef MD_LL_BONDS
  int kk;
  long long int jj, jj2, aa, bb;
#else
  int jj, jj2, kk, nn, aa, bb;
#endif  
  int i, k;
  /* ip=-1 means all bonds */
  if (ip!=-1)
    {
      numbondsMC[ip] = numbonds[ip];
      for (k=0; k < numbonds[ip]; k++)
	{
  	  bondsMC[ip][k] = bonds[ip][k]; 
	  /* ...and store all bonds of bonded particles */
	  jj = bonds[ip][k] / (NANA);
	  //jj2 = bonds[ip][kk] % (NANA);
	  //aa = jj2 / NA;
	  //bb = jj2 % NA;
	  numbondsMC[jj] = numbonds[jj];
	  for (kk=0; kk < numbonds[jj]; kk++)
	    {
	      bondsMC[jj][kk] = bonds[jj][kk];
	    } 

	}
    }
  else
    {
      for (i=0; i < Oparams.parnum; i++)
	{
	  numbondsMC[i] = numbonds[i];
	  for (k=0; k < numbonds[i]; k++)
	    bondsMC[i][k] = bonds[i][k]; 
	}
    }
}
void restore_bonds_mc(int ip)
{
  int i, k;
#ifdef MD_LL_BONDS
  int nb, nn, kk;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, nn, aa, bb;
#endif  

  /* ip=-1 means all bonds */
  if (ip!=-1)
    {
      numbonds[ip] = numbondsMC[ip];
      for (k=0; k < numbondsMC[ip]; k++)
	{
	  bonds[ip][k] = bondsMC[ip][k]; 
	  /* ...and restore all bonds of bonded particles */
	  jj = bonds[ip][k] / (NANA);
	  //jj2 = bonds[ip][kk] % (NANA);
	  //aa = jj2 / NA;
	  //bb = jj2 % NA;
	  numbonds[jj] = numbondsMC[jj];
	  for (kk=0; kk < numbondsMC[jj]; kk++)
	    {
	      bonds[jj][kk] = bondsMC[jj][kk];
	    } 
	}
    }
  else
    {
      for (i=0; i < Oparams.parnum; i++)
	{
	  numbonds[i] = numbondsMC[i];
	  for (k=0; k < numbondsMC[i]; k++)
	    bonds[i][k] = bondsMC[i][k]; 

	}
    }
}
#endif
void move_box(int *ierr)
{
  int i, ii;
  double nn;
#if 1
  double vo, lnvn, vn, Lfact, enn, eno, arg, delv=0.0;
  //printf("moving box\n");
#ifdef MD_LXYZ
  vo = L[0]*L[1]*L[2];
#else
  vo = L*L*L;
#endif
  eno = calcpotene();
  lnvn = log(vo) + (ranf()-0.5)*OprogStatus.vmax;
  vn = exp(lnvn);
  Lfact = pow(vn/vo,1.0/3.0);
  if (OprogStatus.useNNL)
    delv = (Lfact - 1.0); /* FINISH HERE */
 
//	Lfact=1;
#ifdef MD_LXYZ
  L[0] *= Lfact;
  L[1] *= Lfact;
  L[2] *= Lfact;
  L2[0] = L[0]*0.5;
  L2[1] = L[1]*0.5;
  L2[2] = L[2]*0.5;
#else
  L *= Lfact;
  L2 = L*0.5;
#endif

  for (i=0; i < Oparams.parnum; i++)
    {
      rx[i] *= Lfact;
      ry[i] *= Lfact;
      rz[i] *= Lfact; 
      pbc(i);
    }
  update_numcells();
#ifdef MC_STORELL
  store_ll_mc();
#endif
  rebuildLinkedList();
  for (i=0; i < Oparams.parnum; i++)
    {
      if (overlapMC(i, ierr))
	{
	  /* move rejected restore old positions */
#ifdef MD_LXYZ
	  L[0] /= Lfact;
	  L[1] /= Lfact;
	  L[2] /= Lfact;
	  L2[0] = L[0]*0.5;
	  L2[1] = L[1]*0.5;
	  L2[2] = L[2]*0.5;
#else
	  L /= Lfact;
	  L2 = L*0.5;
#endif
	  for (ii=0; ii < Oparams.parnum; ii++)
	    {
	      rx[ii] /= Lfact;
	      ry[ii] /= Lfact;
	      rz[ii] /= Lfact; 
	      pbc(ii);
	    }
	  volrejMC++;
	  update_numcells();
#ifdef MC_STORELL
	  restore_ll_mc();
#else
	  rebuildLinkedList();
#endif
	  return;
	}
    }
#ifdef MC_STOREBONDS
  store_bonds_mc(-1);
#endif
  /* update all bonds with new positions */
  for (i=0; i < Oparams.parnum; i++)
    numbonds[i] = 0;
  if (OprogStatus.useNNL)
    find_bonds_flex_NNL();
  else
    find_bonds_flex_all();
  enn = calcpotene();
  //printf("enn-eno/T=%f\n", (enn-eno)/Oparams.T);
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-(Oparams.parnum+1)*log(vn/vo)*Oparams.T);
  if (ranf() > exp(arg))
    {
      /* move rejected restore old positions */
#ifdef MD_LXYZ
      L[0] /= Lfact;
      L[1] /= Lfact;
      L[2] /= Lfact;
      L2[0] = L[0]*0.5;
      L2[1] = L[1]*0.5;
      L2[2] = L[2]*0.5;
#else
      L /= Lfact;
      L2 = L*0.5;
#endif
      
      for (i=0; i < Oparams.parnum; i++)
	{
	  rx[i] /= Lfact;
	  ry[i] /= Lfact;
	  rz[i] /= Lfact; 
	  pbc(i);
	}
      volrejMC++;
      update_numcells();
#ifdef MC_STORELL
      restore_ll_mc();
#else
      rebuildLinkedList();
#endif
      /* restore all bonds*/
#ifdef MC_STOREBONDS
      restore_bonds_mc(-1);
#else
      for (i=0; i < Oparams.parnum; i++)
	numbonds[i] = 0;
      if (OprogStatus.useNNL)
	find_bonds_flex_NNL();
      else
	find_bonds_flex_all();
#endif
      return;
    }
#endif
  if (OprogStatus.useNNL)
    {
      /* move accepted update displacements */
      for (i=0; i < Oparams.parnum; i++)
	{
	  nn = sqrt(Sqr(rx[i])+Sqr(ry[i])+Sqr(rz[i]));
	  overestimate_of_displ[i] += nn*delv; 
	}
    }
}
extern void remove_bond(int , int , int , int);
void find_bonds_one(int ip);
void update_bonds_MC(int ip)
{
#ifdef MD_LL_BONDS
  int nb, nn, kk;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, nn, aa, bb;
#endif  
  nb = numbonds[ip];
  for (kk = 0; kk < nb; kk++)
    {
      jj = bonds[ip][kk] / (NANA);
      jj2 = bonds[ip][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      remove_bond(jj, ip, bb, aa);
    }
  numbonds[ip] = 0;
  if (OprogStatus.useNNL)
    find_bonds_one_NLL(ip);
  else
    find_bonds_one(ip);
}
void find_bonds_GC(int ip)
{
  numbonds[ip] = 0;
  if (OprogStatus.useNNL)
    find_bonds_one_NLL(ip);
  else
    find_bonds_one(ip);
}
double calc_maxstep_MC(int i);

double get_max_displ_MC(void)
{
  int i;
  double max=0.0, displ, ms;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (OprogStatus.ensembleMC==1)
	{
	  /* nel caso di simulazione NPT il passo massimo
	     dipende dalla posizione della particella e quindi
	     va ricalcolata ogni volta qua */
	  ms = calc_maxstep_MC(i);
	}
      else
	ms = max_step_MC[i];
      //printf("maxax=%f i=%d ms=%.15G ov=%.15G\n", maxax[i], i, ms, overestimate_of_displ[i]);
      displ = overestimate_of_displ[i]+ms;
      if (i==0 || displ > max)
	{
	  max = displ;
	}
    }
  return max;
}
double calc_maxstep_MC(int i)
{
  double vo, d1, d2, d3, d4, nn;
  d1 = 0.5*OprogStatus.deltaMC;
  d2 = OprogStatus.dthetaMC*0.25001*maxax[i];
  if (OprogStatus.ensembleMC==0)
    return (d1>d2)?d1:d2;
  else 
    {
#ifdef MD_LXYZ
      vo = L[0]*L[1]*L[2];
#else
      vo = L*L*L;
#endif
      d3 = (d1>d2)?d1:d2;
      nn = sqrt(rx[i]+ry[i]+rz[i]);
      d4 = nn*(pow(1.0 + exp(0.5*OprogStatus.vmax)/vo,1.0/3.0)-1.0);	
      return (d4>d3)?d4:d3;
    }
}
#if 1
void calc_all_maxsteps(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      max_step_MC[i] = calc_maxstep_MC(i);
    }

}
#endif

int do_nnl_rebuild(void)
{
  int i=0; 
  double displ; 
  displ = get_max_displ_MC();
  /* NOTA: ms è il massimo displacemente finora mentre 
   max_step_MC[i] è il displacement massimo che ci potrà essere successiva mossa */
  if (displ > OprogStatus.rNebrShell)
    return 1;
  else
    return 0;
}
extern void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3]);
extern double **RtA;
extern double ranf_vb(void);
extern double *dists;
extern int bound(int na, int n, int a, int b);
extern void add_bond(int na, int n, int a, int b);

extern void orient_onsager(double *omx, double *omy, double* omz, double alpha);
extern int nbondsFlex;
extern int *mapbondsaFlex, *mapbondsbFlex; 
extern double calcDistNegSP(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, double *dists, int bondpair);

double find_bonds_covadd(int i, int j)
{
  int nn,  amin, bmin, nbonds;
  double shift[3], dist;
#ifdef MD_LXYZ
  shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
  shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
  shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
  shift[0] = L*rint((rx[i]-rx[j])/L);
  shift[1] = L*rint((ry[i]-ry[j])/L);
  shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
  assign_bond_mapping(i, j);  
  dist = calcDistNegSP(0.0, 0.0, i, j, shift, &amin, &bmin, dists, -1);
  
  nbonds = nbondsFlex;
  for (nn=0; nn < nbonds; nn++)
    {
      if (dists[nn]<0.0 && !bound(i, j, mapbondsaFlex[nn], mapbondsbFlex[nn]))
	{
	  add_bond(i, j, mapbondsaFlex[nn], mapbondsbFlex[nn]);
	  add_bond(j, i, mapbondsbFlex[nn], mapbondsaFlex[nn]);
	}
    }

  //if (dist < 0)
    //printf("inside dist=%.15G\n", dist);
  return dist;
  
}
double find_bonds_fake(int i, int j);
#ifdef MD_SPOT_GLOBAL_ALLOC
extern void BuildAtomPos(int i, double *rO, double **R, double **rat);
#else
extern void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3]);
#endif
int is_bonded_mc(int ip, int numb);
int mcinAVB(int i, int j, int dist_type, double alpha, int *merr)
{
  /* dist_type=0 -> isotropic
     dist_type=1 -> onsager */
  static double calls=0.0, tottrials=0.0;
#ifndef MD_SPOT_GLOBAL_ALLOC
  double ratAll[NA][3]; 
#endif
  const int maxtrials=1000000;
  double rA[3], rat[3], norm, sax, cc[3], ene;
  double shift[3], Rl[3][3], vv[3];
  double ox, oy, oz, d, dx, dy, dz;
  int nb, ierr, bonded, k1, k2, trials, nbold;
  /* N.B. questa routine dipende dalla geometria degli spot, in particolare il volume
     in cui inserire va scelto in relazione alla geometria adottata, per ora qui assumiamo
     due spot sulle basi di uno pseudo-cilindro di elongazione data elongazione (sax, vedi sotto) */
  /* place particle i bonded to bond nb of particle j */
  rA[0] = rx[j];
  rA[1] = ry[j];
  rA[2] = rz[j];
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  RtA[k1][k2] = R[j][k1][k2];
	}
    }
  BuildAtomPos(j, rA, RtA, ratAll);
  sax = typesArr[typeOfPart[j]].sax[0];
  bonded=0;
  trials=0;

  *merr=0;
  /* choose randmoly a spot */
  /* ================================================= */
   while (!bonded)
    {
      nb = ranf()*typesArr[typeOfPart[j]].nspots;
      /* we assume here that at most there is one bond per spot */
      //if (is_bonded_mc(j, nb))
	//return 0;

      for (k1=0; k1 < 3; k1++)
	vv[k1] = ratAll[nb+1][k1] - rA[k1];
      norm = calc_norm(vv);
      for (k1=0; k1 < 3; k1++)
	vv[k1] /= norm;
      for (k1=0; k1 < 3; k1++)
	cc[k1] = rA[k1] + vv[k1]*sax*2.0;
#if 1
      /* chose a random position inside a sphere or cube (removing do...while loop)*/
      do {
	dx = 2.0*(ranf_vb()-0.5);
	dy = 2.0*(ranf_vb()-0.5);
	dz = 2.0*(ranf_vb()-0.5);
      } 
      while (dx*dx+dy*dy+dz*dz > 1);
      rx[i] = cc[0]+dx*sax;
      ry[i] = cc[1]+dy*sax;
      rz[i] = cc[2]+dz*sax;
      pbc(i);
#else
      /* chose a random position inside a sphere (this is of course more efficient than
       using a cube as above) */
      orient(&ox, &oy, &oz);
      /* elevando alla 1/3 ci assicuriamo che la distribuzione è uniforme nella sfera di raggio sax */
      norm = pow(ranf_vb(),1.0/3.0)*sax;
      rx[i] = cc[0]+ox*norm;
      ry[i] = cc[1]+oy*norm;
      rz[i] = cc[2]+oz*norm;
      typeOfPart[i]=0;
      pbc(i);
#endif
      if (dist_type==4)
	{
	  ox = 1;
	  oy = 0;
	  oz = 0;
	}
      else if (dist_type==0)
	orient(&ox, &oy, &oz);
      else 
	orient_onsager(&ox, &oy, &oz, alpha);
      versor_to_R(ox, oy, oz, Rl);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      R[i][k1][k2] = Rl[k1][k2]; 
	    }
	}
      ene = find_bonds_fake(i, j);
      if (ene < 0)
	{
#ifdef MD_LXYZ
	  shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
	  shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
	  shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
	  shift[0] = L*rint((rx[i]-rx[j])/L);
	  shift[1] = L*rint((ry[i]-ry[j])/L);
	  shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
	  if ((d=check_overlap(i, j, shift, &ierr))>0.0)
	    {
	      if (ierr==0)
		bonded=1;
	      else 
		{
		  printf("[mcin] NR failure\n");
		}
	    }
	}

      //printf("i=%d j=%d ene=%f\n", i, j, ene);
      trials++;
      if (trials > maxtrials)
	{
	  *merr=1;
	  /* mcinAVB failed */
	  exit(-1);
	  //return 0;
	}
    }
  tottrials += trials;
  calls += 1.0;
  //printf("trials=%d nb=%d avgtrials=%.15G\n", trials,nb, tottrials/calls);
  return nb;
}
void mcin(int i, int j, int nb, int dist_type, double alpha, int *merr, int fake)
{
  /* dist_type=0 -> isotropic
     dist_type=1 -> onsager */
  const int maxtrials=1000000;
  double bondlen, dist=0.0, rA[3], rat[3], norm, sax, cc[3], ene;
  double shift[3], Rl[3][3], vv[3];
  double ox, oy, oz, d, dx, dy, dz;
  int ierr, bonded, k1, k2, trials, nbold;
  static double tottrials=0, calls=0;
  /* N.B. questa routine dipende dalla geometria degli spot, in particolare il volume
     in cui inserire va scelto in relazione alla geometria adottata, per ora qui assumiamo
     due spot sulle basi di uno pseudo-cilindro di elongazione data elongazione (sax, vedi sotto) */
  /* place particle i bonded to bond nb of particle j */
  rA[0] = rx[j];
  rA[1] = ry[j];
  rA[2] = rz[j];
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  RtA[k1][k2] = R[j][k1][k2];
	}
    }
  BuildAtomPosAt(j, nb+1, rA, RtA, rat);
  for (k1=0; k1 < 3; k1++)
    vv[k1] = rat[k1] - rA[k1];
  norm = calc_norm(vv);
  for (k1=0; k1 < 3; k1++)
    vv[k1] /=norm;
  sax = typesArr[typeOfPart[j]].sax[0];
  for (k1=0; k1 < 3; k1++)
    cc[k1] = rA[k1] + vv[k1]*sax*2.0;
  
  assign_bond_mapping(i,j);
  /* N.B. here we assume one bond per pair */
  if (are_spheres(i,j))
    bondlen = mapSigmaFlex[0];
  else 
    bondlen = 1.05*2.0*(norm+0.5*mapSigmaFlex[0]-sax); /* overestimate */
  //printf("bondlen=%.15G norm=%.15G mapSismaFlex=%.15G\n", bondlen, norm, mapSigmaFlex[0]);
  bonded=0;
  trials=0;
  if (!fake)
    numbonds[i]=0;
  *merr=0;
  while (!bonded)
    {
#if 1
      if (are_spheres(i,j))
	{
	  dx = 2.0*(ranf_vb()-0.5);
	  dy = 2.0*(ranf_vb()-0.5);
      	  dz = 2.0*(ranf_vb()-0.5);
	}
      else
	{
	  /* chose a random position inside a sphere or cube (removing do...while loop)*/
	  do {
	    dx = 2.0*(ranf_vb()-0.5);
	    dy = 2.0*(ranf_vb()-0.5);
	    dz = 2.0*(ranf_vb()-0.5);
	  } 
	  while (dx*dx+dy*dy+dz*dz > 1);
	}
      rx[i] = cc[0]+dx*sax;
      ry[i] = cc[1]+dy*sax;
      rz[i] = cc[2]+dz*sax;
#if 1
      if (Sqr(rx[j]-rx[i])+Sqr(ry[j]-ry[i])+Sqr(rz[j]-rz[i]) > Sqr(2.0*sax+bondlen))
	{
	  trials++;
  	  continue;
	}
#endif
      pbc(i);
#else
      /* chose a random position inside a sphere (this is of course more efficient than
       using a cube as above) */
      orient(&ox, &oy, &oz);
      /* elevando alla 1/3 ci assicuriamo che la distribuzione è uniforme nella sfera di raggio sax */
      norm = pow(ranf_vb(),1.0/3.0)*sax;
      rx[i] = cc[0]+ox*norm;
      ry[i] = cc[1]+oy*norm;
      rz[i] = cc[2]+oz*norm;
      typeOfPart[i]=0;
      pbc(i);
#endif
      if (dist_type==4)
	{
	  ox = 1;
	  oy = 0;
	  oz = 0;
	}
      else if (dist_type==0)
	orient(&ox, &oy, &oz);
      else 
	orient_onsager(&ox, &oy, &oz, alpha);
      versor_to_R(ox, oy, oz, Rl);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      R[i][k1][k2] = Rl[k1][k2]; 
	    }
	}
#if 0
      store_bonds_mc(-1);
#else
      //store_bonds_mc(j);
      nbold = numbonds[j];
#endif
      if (fake)
	{
	  ene = find_bonds_fake(i, j);
	}
      else
	{
	  dist=find_bonds_covadd(i, j);
	  ene = calcpotene_GC(i);
	}
#if 0
      find_bonds_GC(i);
      if (numbonds[i] > 1)
	{
	  restore_bonds_mc(-1);
	  numbonds[i]=0;
	  continue;
	}
#endif
      if (ene < 0)
	{
#ifdef MD_LXYZ
	  shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
	  shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
	  shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
	  shift[0] = L*rint((rx[i]-rx[j])/L);
	  shift[1] = L*rint((ry[i]-ry[j])/L);
	  shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
	  if ((d=check_overlap(i, j, shift, &ierr))>0.0)
	    {
	      if (ierr==0)
		bonded=1;
	      else 
		{
		  printf("[mcin] NR failure\n");
#if 0
		  restore_bonds_mc(-1);
#else
		  if (!fake)
		    {
		      numbonds[j]=nbold;
		    }
		      //restore_bonds_mc(j);
#endif
		  if (!fake)
		    numbonds[i]=0;
		}
	    }
	  else
	    {
#if 0
	      restore_bonds_mc(-1);
#else
	      if (!fake)
		numbonds[j]=nbold;
	      //restore_bonds_mc(j);
#endif
	      if (!fake)
		numbonds[i]=0;
#if 0
	      if (ierr!=0)
		{
		  printf("u%d.u%d=%.10G\n", i, j, R[i][0][0]*R[j][0][0]+R[i][0][1]*R[j][0][1]+R[i][0][2]*R[j][0][2]);
		  printf("distRM=%.15G\n", sqrt(Sqr(rx[i]-rx[j]-shift[0])+Sqr(ry[i]-ry[j]-shift[1])+Sqr(rz[i]-rz[j]-shift[2])));
		  printf("[mcin] NR failure\n");
		  store_bump(i, j);
		}
#endif
	    }
	 // printf("d=%.15G trials #%d\n", d, trials);
	}

      //printf("i=%d j=%d ene=%f\n", i, j, ene);
      trials++;
      if (trials > maxtrials)
	{
	  *merr=1;
	  return;
	}
    }
  tottrials += trials;
  calls += 1.0;
  if (dist < 0.0)
    {
      totdist += dist;
      distcc += 1.0;
    }

  //printf("trials=%d avg_trials=%.15G\n", trials, tottrials/calls);
}
/* bonding moves */
int are_bonded(int i, int j)
{
  int kk;
#ifdef MD_LL_BONDS
  int nb;
  long long int aa, bb, ii, jj, jj2;
#else
  int nb, ii, jj, aa, bb, jj2;
#endif
  for (kk=0; kk < numbonds[i]; kk++)
    {
      jj = bonds[i][kk] / (NANA);
#if 0
      jj2 = bonds[ip][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
#endif
      if (jj==j)
	return 1;
    }
  return 0;


}
extern double vnonbond;
int is_bonded_mc(int ip, int numb);
void remove_allbonds_ij(int ip, int j);

void mcoutin(double beta, double pbias)
{
  int j, i, l, nout, ierr, k1, k2, dorej, nb, nbc;
  double rcn, re, epotenoldj, epotenoldi, ideltae, epotennewi;
  int nbj, done, nin, accetto;
  /* select a molecule j=0...Oparams.parnum-1 */
  j=(int) (ranf()*Oparams.parnum);
  epotenoldj=calcpotene_GC(j);  
  nbj = numbonds[j];
  //printf("[mcoutin]-enter\n");
  if (numbonds[j]==Oparams.parnum)
    return;
  if (numbonds[j]==typesArr[typeOfPart[j]].nspots)
    return; /* tutti i siti sono occupati (assumendo che ogni sito è legato al massimo con un altro sito )*/
  /* sceglie una particella esterna */
  done = 0;
  while (!done)
    {    
      i = ranf()*Oparams.parnum;
      if (i==j)
	continue;
      /* controlla che i non sia già bondata a j */ 
      if (are_bonded(i, j))
	{
	  continue;
	}
      else break;
    }
  /*numero di particelle esterne a Vb */
  nout = Oparams.parnum-1-nbj;
  /*numero di particelle interne a Vb */
  nin = nbj;
  store_coord(i);
  epotenoldi=calcpotene_GC(i);
#ifdef MC_STOREBONDS
  store_bonds_mc(i);
#endif
  /* cerca un bond libero. NOTA: se fossero possibili più legami per spot allora 
     bisogna usare una routine che inserisca a caso fra tutti i bond, ossia bisogna usare una routine
     come mcinAVB che sceglie ad ogni tentativo un bond a caso. Infatti mcin inserisce nel bond prescelto 
     (a caso). */
#ifndef MC_USE_MCINAVB
  nb = (int)(ranf_vb()*typesArr[typeOfPart[j]].nspots);
  /* N.B. qui assumo che al massimo ogni patch abbia un solo legame! */
  if (is_bonded_mc(j, nb))
    return;
#endif
  /* inserisce la particella la particella legata a j con il bond nb appena scelto.
     Notare che passando 1 (=fake) come ultimo argomento mcin cambia la posizione e l'orientazione
     della particella ma non aggiorna i legami di i e j.*/
#ifdef MC_USE_MCINAVB
  nbc=mcinAVB(i, j, 0, -1, &ierr); 
  if (is_bonded_mc(j, nbc))
    {
      restore_coord(i);
      update_LL(i);
      if (OprogStatus.useNNL)
	{
	  remove_from_nnl_MC(i);
	  build_one_nnl_GC(i);
	  overestimate_of_displ[i]=0.0;
	  max_step_MC[i] = calc_maxstep_MC(i);
	}
      return;
    }
#else
  mcin(i, j, nb, 0, -1, &ierr, 1);
#endif
  update_LL(i);
  if (OprogStatus.useNNL)
    {
      remove_from_nnl_MC(i);
      build_one_nnl_GC(i);
      overestimate_of_displ[i]=0.0;
      max_step_MC[i] = calc_maxstep_MC(i);
    }

  dorej = overlapMC(i, &ierr);
  if (dorej)
    {
      //remove_allbonds_ij(i, -2);
      //restore_bonds_mc(i);
      //numbonds[j] = nbj;
      restore_coord(i);
      update_LL(i);
      if (OprogStatus.useNNL)
	{
	  remove_from_nnl_MC(i);
	  build_one_nnl_GC(i);
	  overestimate_of_displ[i]=0.0;
	  max_step_MC[i] = calc_maxstep_MC(i);
	}
      //printf("[mcoutin]-exit\n");
      return;
    }

#if 0
  remove_allbonds_ij(i, -2);
  find_bonds_GC(i);
#else
  update_bonds_MC(i);
#endif
  epotennewi = calcpotene_GC(i);
  ideltae=epotennewi-epotenoldi;
  re = ((1.0-pbias)/pbias)*(OprogStatus.vbond/vnonbond)*exp(-beta*ideltae);
  //printf("vbond=%.15G vnonbond: %.15G\n", OprogStatus.vbond, vnonbond);
  re = re*((double) nout)/(nin+1.0); 
  rcn = ranf();
  accetto = (rcn < re)?1:0;
  if (!accetto)
    {
      restore_coord(i);
      update_LL(i);
#ifdef MC_STOREBONDS
      remove_allbonds_ij(i, -2);
      restore_bonds_mc(i);
#else
      update_bonds_MC(i);
#endif
      //numbonds[j] = nbj;
      if (OprogStatus.useNNL)
	{
	  remove_from_nnl_MC(i);
	  build_one_nnl_GC(i);
	  overestimate_of_displ[i]=0.0;
	  max_step_MC[i] = calc_maxstep_MC(i);
	}
      //printf("[mcoutin]-exit\n");
      return;
    }
}
void remove_allbonds_ij(int ip, int j)
{
  int kk;
#ifdef MD_LL_BONDS
  int i, nb;
  long long int aa, bb, ii, jj, jj2;
#else
  int i, nb, ii, jj, aa, bb, jj2;
#endif
  /* NOTA: se j > 0 rimuove tutti i legami trai i e j
           se j==-1 rimuove tutti i legami di i (anche dalle particelle legate)
	   se j==-2 rimuove tutti i legami di i solo dalle particelle legate ad i e non da i stessa
   */
  nb = numbonds[ip];
  for (kk=0; kk < nb; kk++)
    {
      jj = bonds[ip][kk] / (NANA);
      
      if (j < 0 || jj==j)
	{
	  jj2 = bonds[ip][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (j==-1|| j > 0)
	    remove_bond(ip, jj, aa, bb);
	  remove_bond(jj, ip, bb, aa);
	}
    }
}
double find_bonds_fake(int i, int j)
{
  int nn,  amin, bmin, nbonds;
  double ene, shift[3], dist;
#ifdef MD_LXYZ
  shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
  shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
  shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
  shift[0] = L*rint((rx[i]-rx[j])/L);
  shift[1] = L*rint((ry[i]-ry[j])/L);
  shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
  assign_bond_mapping(i, j);  
  dist = calcDistNegSP(0.0, 0.0, i, j, shift, &amin, &bmin, dists, -1);
  nbonds = nbondsFlex;
  ene=0;
  for (nn=0; nn < nbonds; nn++)
    {
      if (dists[nn]<0.0)
	{
	  ene -= 1.0;
	}
    }
  return ene;

}
void mcout(int i, int j, int nb)
{
  double ene, ox, oy, oz, Rl[3][3];
  int k1, k2, bonded=1, nbjold;

  //store_bonds_mc(i);
  //store_bonds_mc(j);

  //numbonds[i]=0;
  //numbonds[j]=0;
  while (bonded)
    {
#ifdef MD_LXYZ
      rx[i] = L[0]*(ranf_vb()-0.5);
      ry[i] = L[1]*(ranf_vb()-0.5); 
      rz[i] = L[2]*(ranf_vb()-0.5); 
#else
      rx[i] = L*(ranf_vb()-0.5);
      ry[i] = L*(ranf_vb()-0.5); 
      rz[i] = L*(ranf_vb()-0.5); 
#endif
      orient(&ox, &oy, &oz);
      versor_to_R(ox, oy, oz, Rl);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      R[i][k1][k2] = Rl[k1][k2]; 
	    }
	}
      ene=find_bonds_fake(i, j);
      if (ene==0)
	{
	  bonded=0;
	}
    }
  //restore_bonds_mc(j);
  //numbonds[i]=0;
  /* rimuove da j tutti i bond tra i e j */
  //remove_all_bonds(j, i);
}
void mcinout(double beta, double pbias)
{
  int i, l, j, nbt;
  int dorej, accetto, ierr, nout, nin;
  double re, rcn, epotenoldj, epotenoldi, epotennewi, ideltae;
  j = (int) (ranf()*Oparams.parnum);
  if (j >= Oparams.parnum)
    j=Oparams.parnum-1;
  epotenoldj = calcpotene_GC(j);
  nbt = numbonds[j];
  if (nbt==0)
    return; /* nessuna particella legata esci */
  /* N.B. sceglie a caso un legame da rompere, 
     assumendo una sola particella per legame */
  l = ranf()*nbt;
  i = bonds[j][l] / (NANA);  
  nout = Oparams.parnum-nbt-1;
  nin = nbt;
  epotenoldi = calcpotene_GC(i);
  store_coord(i);
#ifdef MC_STOREBONDS
  store_bonds_mc(i);
  //store_bonds_mc(j);
#endif
  /* N.B.: notare che mcout posiziona la particella ma non aggiorna i legami */
  mcout(i, j, l);
  /* rimuove dai bond di j tutti i bond tra j ed i */
  //printf("boh numbonds[%d]=%d numbonds[%d]=%d\n", i, numbonds[i], j, numbonds[j]);
  update_LL(i);
  //printf("[mcinout]-enter\n");
  if (OprogStatus.useNNL)
    {
      /* rimuove i da tutte le NNL del sistema */
      remove_from_nnl_MC(i);
      build_one_nnl_GC(i);
      overestimate_of_displ[i]=0.0;
      max_step_MC[i] = calc_maxstep_MC(i);
    }
  dorej = overlapMC(i, &ierr);
  if (dorej)
    {
      /* rimuove tutti i legami di i */
      restore_coord(i);
      update_LL(i);
      //numbonds[j] = nbt;
      if (OprogStatus.useNNL)
	{
	  /* rimuove i da tutte le NNL del sistema */
	  remove_from_nnl_MC(i);
	  build_one_nnl_GC(i);
	  overestimate_of_displ[i]=0.0;
	  max_step_MC[i] = calc_maxstep_MC(i);
	}
      //printf("[mcinout]-exit\n");
      return;
    }

  update_bonds_MC(i);
  epotennewi = calcpotene_GC(i);
 
  ideltae=epotennewi-epotenoldi;
  re = (pbias/(1.0-pbias))*(vnonbond/OprogStatus.vbond)*exp(-beta*ideltae);
  re = re*((double) nin)/(nout+1.0); 
  rcn = ranf();
  accetto = (rcn < re)?1:0;
  if (!accetto)
    {
      restore_coord(i);
      update_LL(i);
#ifdef MC_STOREBONDS
      remove_allbonds_ij(i, -2);
      restore_bonds_mc(i);
      //restore_bonds_mc(j);
#else
      update_bonds_MC(i);
      //update_bonds_MC(j);
#endif
      //numbonds[j] = nbt;
      if (OprogStatus.useNNL)
	{
	  /* rimuove i da tutte le NNL del sistema */
	  remove_from_nnl_MC(i);
	  build_one_nnl_GC(i);
	  overestimate_of_displ[i]=0.0;
	  max_step_MC[i] = calc_maxstep_MC(i);
	}
      //printf("[mcinout]-exit\n");
      return;
    }
}
int is_bonded_mc(int ip, int numb)
{
  int kk;
#ifdef MD_LL_BONDS
  int i, nb;
  long long int aa, bb, ii, jj, jj2;
#else
  int i, nb, ii, jj, aa, bb, jj2;
#endif
  for (kk=0; kk < numbonds[ip]; kk++)
    {
      jj = bonds[ip][kk] / (NANA);
      jj2 = bonds[ip][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (aa==numb+1)
	return 1;
    }
  return 0;
}
#ifdef MC_CALC_COVADD
extern const char sepStr[];
void save_conf_mc(int i)
{
  char fileop2[512], fileop[512], fileop3[512];
  FILE* bf;
  sprintf(fileop2 ,"Store-%d-%d", 
	  i, 0);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
  for (i=0; i < Oparams.parnum; i++)
    {
      update_MSDrot(i);
#ifdef MD_CALC_DPP
      update_MSD(i);
      store_last_u(i);
#endif
      OprogStatus.lastcolltime[i] = Oparams.time;
    }
  R2u();
#ifdef MD_SAVE_DISTANCE
  if (mgl_mode==1)
    {
      FILE *f;
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
      double dists[MD_PBONDS], d;
      int i;
#endif
      f = fopen("distance.dat", "a");
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
      d=calcJustDistNegSP(Oparams.time, 0, 1, dists);
      fprintf(f, "%.15G %.15G ", Oparams.time + OprogStatus.refTime, calcJustDistNeg(Oparams.time, 0, 1));
#if 1
      for (i=0; i < MD_PBONDS; i++)
	{
	  fprintf(f, "%.15G ", dists[i]);
	}
#else
      fprintf(f, "%.15G ", dists[1]);
      printf("mapbonds[0]:%d mapbondsb[1]:%d\n", mapbondsa[0], mapbondsb[1]);
#endif
      fprintf(f, "\n");
#else
      fprintf(f, "%.15G %.15G\n", Oparams.time + OprogStatus.refTime, calcJustDistNeg(Oparams.time, 0, 1));
#endif
      fclose(f);
    }
#endif

  if (mgl_mode==0)
    {
      writeAsciiPars(bf, opro_ascii);
      fprintf(bf, sepStr);
      writeAsciiPars(bf, opar_ascii);
      fprintf(bf, sepStr);
    }	      
  MD_DEBUG(printf("[Store event]: %.15G JJ=%d KK=%d\n", Oparams.time, OprogStatus.JJ, OprogStatus.KK));
  //fprintf(bf, ".semiAxes: %f %f %f, %f %f %f\n",
  //	  Oparams.a[0], Oparams.b[0], Oparams.c[0],
  //  Oparams.a[1], Oparams.b[1], Oparams.c[1]);
  writeAllCor(bf, 0);
  fclose(bf);
  if (mgl_mode==0)
    {
#ifdef MPI
#ifdef MD_MAC
      sprintf(fileop3, "/usr/bin/gzip -f %s_R%d", fileop, my_rank);
#else
      sprintf(fileop3, "/bin/gzip -f %s_R%d", fileop, my_rank);
#endif
#else 
#ifdef MD_MAC
      sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
#else
      sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
#endif
#ifndef MD_NO_SYSTEM
      system(fileop3);
#endif	    
    }
}

void build_ll(int *pl)
{
  int i, i0=-1, iold, in0, in1;
  //printf("build ll\n");
  for (i=0; i < Oparams.parnum; i++)
    {
      if ( numbonds[i]==1) 
	{
  	  i0 = i;
	  break;
	}
    }
  for (i=0; i < Oparams.parnum+1; i++)
    pl[i] = -1;
  pl[Oparams.parnum] = i0; /* list head */
  iold = i0;
  i =  bonds[i0][0] / (NANA);
  while (1)
    {
      pl[i] = pl[Oparams.parnum]; 
      pl[Oparams.parnum] = i;
      if (numbonds[i]==1)
	break;
      in0 = bonds[i][0]/(NANA);
      in1 = bonds[i][1]/(NANA);
      if (in1==iold)
	{
	  iold = i;
	  i = in0;
	}
      else
	{
	  iold = i;
	  i = in1;
	}
    }

  //printf("done build ll\n");
}
void accum_persist_len_mono_orient(int *parlist, double *pl, double *cc)
{
#ifdef MD_LL_BONDS
  int nb;
  long long int aa, bb, ii, jj, jj2, a1=-1, a2=-1;
#else
  int nb, ii, jj, aa, bb, jj2, a1=-1, a2=-1;
#endif
  double Ri[3], Rj[3], norm;
  double rA[3], r2[3], r1[3];	
  int k, a, i, j, NP, np, in, jn, k1, k2, kk, ip, jp;
  double costh2, normi, normj;
#ifndef MD_SPOT_GLOBAL_ALLOC
  double ratAll[NA][3]; 
#endif
 
  NP=Oparams.parnum;
  build_ll(parlist);
  i = Oparams.parnum;
  //printf("INIZIO\n");

  ip = -1;
  while ((i=parlist[i])!=-1) 
    {
      // printf("i=%d\n", i);
      j=i;
      np=0;
      in = parlist[i];
      if (in==-1)
	break;
      rA[0] = rx[i];
      rA[1] = ry[i];
      rA[2] = rz[i];
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      RtA[k1][k2] = R[i][k1][k2];
	    }
    	}
      BuildAtomPos(i, rA, RtA, ratAll);
      for (kk = 0; kk < numbonds[i]; kk++)
	{
	  jj = bonds[i][kk] / (NANA);
    	  jj2 = bonds[i][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  /* we assume two bonds max here */
      	  if (in!=-1 && jj==in)
	    {
	      a2 = aa;
	      a1 = (aa==1)?2:1;
	    }
	  else if (ip == jj)
	    {
	      a1 = aa;
	      a2 = (aa==1)?2:1; 
	    }
	  /* we assume two bonds max here */
	}
      for (k1=0; k1 < 3; k1++)
	{
	  r2[k1] = ratAll[a2][k1];
	  r1[k1] = ratAll[a1][k1];
	} 
      for (k1=0; k1 < 3; k1++)
	Ri[k1] = r2[k1]-r1[k1];
      norm=calc_norm(Ri);
      for (k1=0; k1 < 3; k1++)
	Ri[k1] /= norm; 
      jp = j;
      while ((j=parlist[j])!=-1)
	{
	  //printf("i=%d j=%d\n", i, j);
	  costh2 = 0.0;
	  np++;
	  jn = parlist[j];
	  if (jn==-1)
	    break;
	  rA[0] = rx[j];
	  rA[1] = ry[j];
	  rA[2] = rz[j];
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      for (k2=0; k2 < 3; k2++)
		{
		  RtA[k1][k2] = R[j][k1][k2];
		}
	    }
	  BuildAtomPos(j, rA, RtA, ratAll);

	  for (kk = 0; kk < numbonds[j]; kk++)
	    {
	      jj = bonds[j][kk] / (NANA);
	      jj2 = bonds[j][kk] % (NANA);
	      aa = jj2 / NA;
	      bb = jj2 % NA;
	      /* we assume two bonds max here (aa=1,2) */
	      if (jn!=-1 && jj==jn)
		{
		  a2 = aa;
		  a1 = (aa==1)?2:1;
		}
	      else if (jp == jj)
		{
		  a1 = aa;
		  a2 = (aa==1)?2:1; 
		}
	    }
	  for (k1=0; k1 < 3; k1++)
    	    {
	      r2[k1] = ratAll[a2][k1];
	      r1[k1] = ratAll[a1][k1];
	    } 
	  for (k1=0; k1 < 3; k1++)
	    Rj[k1] = r2[k1] - r1[k1];
	  norm=calc_norm(Rj);
	  for (k1=0; k1 < 3; k1++)
	    Rj[k1] /= norm; 
     
	  for (a = 0; a < 3; a++)
	    {
	      costh2 += Ri[a]*Rj[a];
	    }
	  pl[np] += costh2;
	  //printf("np=%d costh2=%.15G pl=%f\n", np, costh2, pl[np]);
	  cc[np] += 1.0;
	  jp = j;
	}
      ip = i;
    }
}


void accum_persist_len(int *parlist, double *pl, double *cc)
{
  double Ri[3], Rj[3];
  int a, i, j, NP, np, in, jn;
  double costh2, normi, normj;

  NP=Oparams.parnum;
  build_ll(parlist);
  i = Oparams.parnum;
  //printf("INIZIO\n");
  while ((i=parlist[i])!=-1) 
    {
  //    printf("i=%d\n", i);
      j=i;
      np=0;
      in = parlist[i];
      if (in==-1)
	break;
      Ri[0] = rx[in] - rx[j];
      Ri[1] = ry[in] - ry[j];
      Ri[2] = rz[in] - rz[j];
      normi=calc_norm(Ri);
      for (a=0; a < 3; a++)
	{
	  Ri[a] /= normi;
	} 
			 

      while ((j=parlist[j])!=-1)
	{
	  //printf("i=%d j=%d\n", i, j);
	  costh2 = 0.0;
	  np++;
	  jn = parlist[j];
	  if (jn==-1)
	    break;
	  Rj[0] = rx[jn] - rx[j];
	  Rj[1] = ry[jn] - ry[j];
	  Rj[2] = rz[jn] - rz[j];
	  normj=calc_norm(Rj);
	  for (a=0; a < 3; a++)
	    {
	      Rj[a] /= normj;
	    } 
	
	  //printf("j=%d bonds[%d]/NANA=%lld %lld\n", j, i, bonds[i][0]/NANA, bonds[i][1]/NANA);
	  for (a = 0; a < 3; a++)
    	    {
	      costh2 += Ri[a]*Rj[a];
	    }
	  pl[np] += costh2;
	  //printf("np=%d costh2=%.15G pl=%f\n", np, costh2, pl[np]);
	  cc[np] += 1.0;
	}
    }
}
void calc_persistence_length_mc(long long int maxtrials, int outits, int size1)
{
  int abort=0, *parlist, i, j, nb, k1, k2, ierr, merr, jj;
  long long int tt;
  double dist, *pl, *cc, shift[3];
  FILE *fi, *f;
  /* first particle is always in the center of the box with the same orientation */
  printf("calculating persistence length\n");
  cc = malloc(sizeof(double)*Oparams.parnum);
  pl = malloc(sizeof(double)*Oparams.parnum);
  parlist = malloc(sizeof(int)*(Oparams.parnum+1));
  for (i=0; i < Oparams.parnum; i++)
    {
      cc[i] = pl[i] = 0;
      numbonds[i] = 0;
    }
  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
      }

  f = fopen("perslenlast.dat","w+");
  fclose(f);
  for (tt=0; tt < maxtrials; tt++)
    {
      abort=0;
      if (tt%outits==0)
	{
	  
	  printf("tt=%lld\n", tt);
	  if (distcc > 0)
	    printf("average bond distance=%.15G\n", totdist/distcc); 
	  if (tt!=0 && cc[Oparams.parnum-2]!=0)
	    {
	      f = fopen("perslenlast.dat", "a");
	      fprintf(f, "%lld %.15G\n", tt, pl[Oparams.parnum-2]/cc[Oparams.parnum-2]);
	      fclose(f);
	      fi = fopen("persist.dat","w+");
	      pl[0]=1.0;
	      for (i=(size1!=0)?0:1; i < Oparams.parnum-1; i++)
		{
		  if (size1!=0 && (cc[i]!=0 || i==0))
		    {
		      if (i==0)
			fprintf(fi, "0 1\n");
		      else
			fprintf(fi, "%d %.15G\n", i, pl[i]/cc[i]);
		    }
		  else
	      	    { 
	    	      if ((cc[i]!=0 && cc[i-1]!=0)||i==1)
	    		{
	    		  if (i==1)
	    		    fprintf(fi, "0 1\n");
	    		  else
	    		    fprintf(fi, "%d %.15G\n", i-1, 0.5*(pl[i]/cc[i]+pl[i-1]/cc[i-1]));
	    		}
		    }
		}
	      fclose(fi);	
	      sync();
	    }
	}
      for (i=0; i < Oparams.parnum; i++)
	numbonds[i]=0;
      for (i=1; i < Oparams.parnum; i++)
	{
	  while (1)
	    {
	      nb = (int)(ranf_vb()*2.0);
	      j = (int) (ranf_vb()*i);
#if 0
	      if (numbonds[j]==2)
		{
		  printf("j=%d has 2 bonds!\n", j);
		  exit(-1);
		}
#endif
	      if (is_bonded_mc(j, nb))
		continue;
	      else
		break;
	    }
	  mcin(i, j, nb, 0, 0.0, &merr, 0);
	  if (merr!=0)
	    {
	      printf("[mcin] attempt to add a bonded particle failed!\n");
	      /* mc in ha fallito scarta tale conf */
	      abort=1;
	      break;
	    }
#if 1
	  /* qui controlla che non ci siano overlap con le particelle
	     già inserite e che mcin abbia formato un solo legame */
	  for (jj=0; jj < i; jj++)
	    {
	      if (jj==j)
		continue;
	      find_bonds_covadd(i, jj);
      	    }
	  if (numbonds[i] > 1)
	    {
#if 0
	      restore_bonds_mc(-1);
	      numbonds[i]=0;
	      i--;
	      continue;
#else
	      /* scarta tale configurazione "chiusa" (loop) */
	      abort=1;
	      break;
#endif
	    }
	  for (jj=0; jj < i; jj++)
	    {
	      if (jj==j) 
		continue;
#ifdef MD_LXYZ
	      shift[0] = L[0]*rint((rx[i]-rx[jj])/L[0]);
	      shift[1] = L[1]*rint((ry[i]-ry[jj])/L[1]);
	      shift[2] = L[2]*rint((rz[i]-rz[jj])/L[2]);
#else
	      shift[0] = L*rint((rx[i]-rx[jj])/L);
	      shift[1] = L*rint((ry[i]-ry[jj])/L);
	      shift[2] = L*rint((rz[i]-rz[jj])/L);
#endif
	      if (check_overlap(i, jj, shift, &ierr) < 0.0)
		{
#if 0
		  restore_bonds_mc(-1);
		  numbonds[i]=0;
		  i--;
#endif
		  /* scarta la conf self-overlapping */
		  abort=1; 
		  break;
		}
	    }
#endif
	  if (abort)
	    break;
	}
#if 0
      save_conf_mc(tt);
#else
      if (!abort)
	{
	  /* N.B. size1 viene solo usato per scegliere il modo in cui viene calcolata 
	     la funzione di correlazione dei bond */
	  if (size1!=0)
	    accum_persist_len_mono_orient(parlist, pl, cc);
	  else
	    accum_persist_len(parlist, pl, cc);
	}
#endif
    }
  fi = fopen("persist.dat","w+");
  
  for (i=(size1==0)?1:0; i < Oparams.parnum-1; i++)
    {
      /* la media della funzione di correlazione serve per avere una stima 
	 più vicina alla funzione di correlazione delle orientazioni dei monomeri (il -1 serve per iniziare
	 correttamente da s=0), ossia se voglio avere <u(0)*u(s)> dove u è l'orientazione di un monomero allora
	 <bu(0)*(bu(s)+bu(s-1))*0.5> dove ora bu è la direzione di calcolata usando i centri di 
	 massa dei monomeri e bu(-1)=0 */
      if (size1!=0)
	{
	  if (i==0)
	    fprintf(fi, "0 1\n");
	  else
	    fprintf(fi, "%d %.15G\n", i, pl[i]/cc[i]);
	}
      else
	{
	  if (i==1)
	    fprintf(fi, "0 1\n");
	  else
	    fprintf(fi, "%d %.15G\n", i-1, 0.5*(pl[i]/cc[i]+pl[i-1]/cc[i-1]));
	}
    }
  fclose(fi);	
}
int check_overlap_all(int i0, int i_ini, int i_fin);
int insert_remaining_particles(void)
{
  int i, j, tt2, maxattempts = 1000000000;
  long long int totattempts=0;
  int k1, k2;
  double ox, oy, oz, Rl[3][3];
  /* insert all remaining size1-1  particles randomly */
  for (i=2; i < Oparams.parnum; i++)
    {
      for (tt2=0; tt2 < maxattempts; tt2++)
	{
	  orient(&ox, &oy, &oz);
	  versor_to_R(ox, oy, oz, Rl);
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		R[i][k1][k2] = Rl[k1][k2];
	      }
	  numbonds[i] = 0;
#ifdef MD_LXYZ
	  rx[i] = (ranf_vb()-0.5)*L[0];
	  ry[i] = (ranf_vb()-0.5)*L[1];
	  rz[i] = (ranf_vb()-0.5)*L[2];	
#else
	  rx[i] = (ranf_vb()-0.5)*L;
	  ry[i] = (ranf_vb()-0.5)*L;
	  rz[i] = (ranf_vb()-0.5)*L;	
#endif
	  if (!check_overlap_all(i, 1, i-1))
	    break;
	}
      totattempts+=tt2+1;
      if (tt2==maxattempts)
	return 0;/* 0 =failed */
    }
  //printf("average attempts = %f\n", ((double)totattempts)/(Oparams.parnum-2));
  return 1;
}
int check_overlap_all(int i0, int i_ini, int i_fin)
{
  double shift[3];
  int i, ierr=0;
  for (i=i_ini; i <= i_fin; i++)
    {
#ifdef MD_LXYZ
      shift[0] = L[0]*rint((rx[i0]-rx[i])/L[0]);
      shift[1] = L[1]*rint((ry[i0]-ry[i])/L[1]);
      shift[2] = L[2]*rint((rz[i0]-rz[i])/L[2]);
#else
      shift[0] = L*rint((rx[i0]-rx[i])/L);
      shift[1] = L*rint((ry[i0]-ry[i])/L);
      shift[2] = L*rint((rz[i0]-rz[i])/L);
#endif

      if (check_overlap(i0, i, shift, &ierr) < 0.0)
	{
	  if (ierr==0)
	    {
	      return 1;
	    }
	} 
    }
  return 0;
}
#undef MC_VB_PLOT
void calc_bonding_volume_mc(long long int maxtrials, int outits, int type, double alpha)
{
  FILE *f;	
  int k1, k2, ierr, i;
  long long int tt;
  double deldistcc=0.0, deltotdist=0.0, dist, fact, shift[3], Lb, totene=0.0, ox, oy, oz, Rl[3][3];
#ifdef MC_VB_PLOT
  FILE *vbf;
#endif
  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
#ifdef MC_VB_PLOT
  vbf = fopen("bonding-volume.dat", "w+");
#endif
  printf("calc vbonding MC\n");
  
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
      }
  f = fopen("vbonding.dat","w+");
  fclose(f);
  totene = toteneini;
  for (tt=ttini; tt < maxtrials; tt++)
    {
      if (tt%outits==0 && tt!=0)
	{
	  printf("tt=%lld\n", tt); 
	  if (distcc > 0)
	    {
	      printf("Bonding distance=%.15G\n", totdist/distcc);
	    }
	  f=fopen("vbonding.dat", "a");
	  fact=Oparams.parnum-1;
#ifdef MD_LXYZ
	  fprintf(f, "%lld %.15G %G\n", tt, (totene/((double)tt))*(L[0]*L[1]*L[2])/Sqr(typesArr[0].nspots)/fact, totene);
#else
	  fprintf(f, "%lld %.15G %G\n", tt, (totene/((double)tt))*(L*L*L)/Sqr(typesArr[0].nspots)/fact, totene);
#endif
	  sync();
	  fclose(f);
	}
      numbonds[1] = numbonds[0] = 0;
#ifdef MD_LXYZ
      rx[1] = (ranf_vb()-0.5)*L[0];
      ry[1] = (ranf_vb()-0.5)*L[1];
      rz[1] = (ranf_vb()-0.5)*L[2];	
#else
      rx[1] = (ranf_vb()-0.5)*L;
      ry[1] = (ranf_vb()-0.5)*L;
      rz[1] = (ranf_vb()-0.5)*L;	
#endif
      if (type==5)
	{
	  orient_onsager(&ox,&oy,&oz,alpha);
	  versor_to_R(ox, oy, oz, Rl);
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      for (k2=0; k2 < 3; k2++)
		{
		  R[0][k1][k2] = Rl[k1][k2]; 
		}
	    }
	}

      if (type==5)
	orient_onsager(&ox,&oy,&oy,alpha);
      else
	orient(&ox, &oy, &oz);
      versor_to_R(ox, oy, oz, Rl);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      R[1][k1][k2] = Rl[k1][k2]; 
	    }
	}
      if (Oparams.parnum > 2)
	{
	  if (!insert_remaining_particles())
	    {
	      /* se arriva qui vuold dire che non è riuscito a creare
		 la conf alla volume di size1 particelle */
	      printf("Failed to create background configuration\n");
	      exit(-1);
	    }
	}

      if (Oparams.parnum > 2)
	{
	  deltotdist=0.0;
	  deldistcc=0.0;
	  for (i=1; i < Oparams.parnum; i++) 
	    {
	      dist=find_bonds_covadd(0, i);
	      if (dist < 0.0)
		{
		  deltotdist += dist;
		  deldistcc += 1.0;
		}
	    }
	}
      else
	  find_bonds_covadd(0, 1);
      if (numbonds[0] > 0)
	{	
	  if (Oparams.parnum > 2)
	    {
	      if (!check_overlap_all(0, 1, Oparams.parnum-1))
		{
		  totene += numbonds[0];		  
		  totdist += deltotdist;
		  distcc += deldistcc;
		}
	    }
	  else
	    {
#ifdef MD_LXYZ
	      shift[0] = L[0]*rint((rx[0]-rx[1])/L[0]);
	      shift[1] = L[1]*rint((ry[0]-ry[1])/L[1]);
	      shift[2] = L[2]*rint((rz[0]-rz[1])/L[2]);
#else
    	      shift[0] = L*rint((rx[0]-rx[1])/L);
    	      shift[1] = L*rint((ry[0]-ry[1])/L);
    	      shift[2] = L*rint((rz[0]-rz[1])/L);
#endif

	      if (check_overlap(0,1, shift, &ierr) > 0.0)
		{
		  if (ierr==0)
		    {
#ifdef MC_VB_PLOT
		      fprintf(vbf, "%.15G %.15G %.15G\n", rx[1], ry[1], rz[1]);
#endif
		      totene += 1.0;
		    }
		} 
	    }
	}
#if 0
      if (tt%500000==0)
	{
	  printf("tt=%lld\n", tt); 
	}
#endif
    }
  if (Oparams.parnum > 2)
    {
#ifdef MD_LXYZ
      printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)tt))*(L[0]*L[1]*L[2])/Sqr(typesArr[0].nspots)/Oparams.parnum-1, totene);
#else
      printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)tt))*(L*L*L)/Sqr(typesArr[0].nspots)/Oparams.parnum-1, totene);
#endif
    }
  else
    {
#ifdef MD_LXYZ
      printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)tt))*(L[0]*L[1]*L[2])/Sqr(typesArr[0].nspots), totene);
#else
      printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)tt))*(L*L*L)/Sqr(typesArr[0].nspots), totene);
#endif
    }
  //fclose(f);
#ifdef MC_VB_PLOT
  fclose(vbf);
#endif
 
}
extern double *distro;
int check_self_overlap(int i0, int i)
{
  int j;
  double shift[3];	
  int overlap, ierr;
  overlap = 0;
  for (j=i0; j < i; j++)
    {
#ifdef MD_LXYZ
      shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
      shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
      shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
      shift[0] = L*rint((rx[i]-rx[j])/L);
      shift[1] = L*rint((ry[i]-ry[j])/L);
      shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
      ierr=0;
      if (check_overlap(i, j, shift, &ierr)<0.0)
	{
	  overlap=1;
	  //printf("i=%d j=%d overlap!\n", i, j);
	  //printf("shift=%f %f %f\n", shift[0], shift[1], shift[2]);
	  //printf("r=%f %f %f  j %f %f %f\n", rx[i], ry[i], rz[i], rx[j], ry[j], rz[j]);
	  //printf("self overlap!\n");
	  return 1;
	}
    }
   return 0; 
}
extern double fons(double theta, double alpha);
void calc_cov_additive(void)
{
  FILE *fi, *f=NULL;
  double cov, Lb, totene = 0.0, alpha, shift[3];
  int bt=0, merr=0, i, j=-1, selfoverlap=0, size1, size2, nb, k1, k2, overlap=0, ierr, type, outits;
  long long int maxtrials, tt;
  double ox, oy, oz, Rl[3][3];
  fi = fopen("covmc.conf", "r");
  fscanf(fi, "%lld %d %d %d ", &maxtrials, &type, &size1, &outits);
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i]=0;
    }
  
  printf("ene iniziale=%f\n", calcpotene());
  if (type==1||type==5)
    fscanf(fi, " %lf ", &alpha);
  /* if it is a restart put info initial values here */
  covrestart = 0;
  if (!feof(fi))
    {
      covrestart=1;
      fscanf(fi, "%lf %lld ", &toteneini, &ttini);
    }
   
  /* type = 0 -> covolume 
     type = 1 -> covolume nematic
     type = 2 -> persistence length
     type = 3 -> bonding volume
     type = 4 -> covolume if perfectly aligned (i.e. alpha -> infinity)
     type = 5 -> bonding volume using onsager trial function
   */
  if (type==1||type==5)
    {
      const int n=1000;
      distro=malloc(sizeof(double)*n);
      for (i=0; i < n; i++)
	distro[i] = 0.0;
    }
  
  if (size1 >= Oparams.parnum)
    {
      printf("size1=%d must be less than parnum=%d\n", size1, Oparams.parnum);
      exit(-1);
    } 
  fclose(fi);
  init_rng(-1, 0, -1);
  OprogStatus.optnnl = 0;
  tt = ttini;
  totene = toteneini;
  size2 = Oparams.parnum-size1;
  if (type==3||type==5)
    {
      calc_bonding_volume_mc(maxtrials, outits, type, alpha);
      exit(-1);
    }
  if (type==2)
    {
      calc_persistence_length_mc(maxtrials, outits, size1);
      exit(-1);
    }
  if (type==4)
    {
      ox=1;
      oy=0;
      oz=0; 
      versor_to_R(ox, oy, oz, Rl);
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    R[0][k1][k2] = Rl[k1][k2];
	  }
    }
  else
    {
      /* first particle is always in the center of the box with the same orientation */
      rx[0] = 0;
      ry[0] = 0;
      rz[0] = 0;
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    R[0][k1][k2] = (k1==k2)?1:0;
	  }
    }
  if (type==0||type==4)
    f = fopen("covolume.dat","w+");
  else if (type==1)
    f = fopen("covolume-nem.dat","w+");
  fclose(f);
    
  while (tt < maxtrials) 
    {
      merr=ierr=0;
      selfoverlap=0;
      
      if (type==1)
	{
	  /* first particle is always in the center of the box with the same orientation */
	  rx[0] = 0;
	  ry[0] = 0;
	  rz[0] = 0;
	  orient_onsager(&ox, &oy, &oz, alpha); 
	  versor_to_R(ox, oy, oz, Rl);
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		R[0][k1][k2] = Rl[k1][k2];
	      }
	}
      /* place first cluster */
      if (tt%outits==0)
	{
	  if (tt!=0)
	    {
#ifdef MD_LXYZ
	      cov = (totene/((double)tt))*(L[0]*L[1]*L[2]);
	    
#else
	      cov = (totene/((double)tt))*(L*L*L);
#endif
	      if (type==0||type==4)
		f=fopen("covolume.dat", "a");
	      else
		f=fopen("covolume-nem.dat", "a");
	      printf("co-volume=%.10f (totene=%f/%lld)\n", cov, totene, tt);
	      fprintf(f, "%lld %.15G %.15G\n", tt, cov, totene);
	      fclose(f);
	      sync();
	    }
#if 0
     	  save_conf_mc(tt); 
#endif
	}
      for (i=0; i < Oparams.parnum; i++)
	{
	  numbonds[i] = 0;
	}
      for (i=1; i < size1; i++)
	{
	  bt = 0;
	  while (1)
	    {
	      nb = (int)(ranf_vb()*2.0);
	      j = (int) (ranf_vb()*i);
#if 1
	      if (bt > Oparams.parnum*1000)
		{
		  printf("insertion cluster #1: maximum number of iteration reached\n");
    		  if (fabs(calcpotene()+i*2.0) < 1E-10)
		    {
		      save_conf_mc(0);
		      printf("During insertion of #1 cluster: every particle has 2 bonds! i=%d\n", i);
		      exit(-1);
		    }
		}
#endif
	      if (is_bonded_mc(j, nb))
		continue;
	      else
		break;
	      bt++;
	    }
	  mcin(i, j, nb, type, alpha, &merr, 0);
	  if (merr!=0)
	    {
	      save_conf_mc(0);
	      printf("[mcin] attemp to add a bonded particle failed!\n");
	      break;
	    }
	  if (check_self_overlap(0, i))
	    {
	      selfoverlap = 1;
	      break;
	    }
	  /* N.B. per ora non controlla il self-overlap della catena 
	     e la formazione dopo mcin di legami multipli poiché
	     si presuppone che al massimo stiamo considerando dimeri */
	}

      if (selfoverlap)
	{
	  tt++;
	  continue;
	}
      /* place second cluster */
      overlap=0;
      for (i=size1; i < Oparams.parnum; i++)
	{
	  if (i==size1)
	    {
#ifdef MD_LXYZ
	      rx[i] = L[0]*(ranf_vb()-0.5);
	      ry[i] = L[1]*(ranf_vb()-0.5); 
	      rz[i] = L[2]*(ranf_vb()-0.5); 
#else
    	      rx[i] = L*(ranf_vb()-0.5);
    	      ry[i] = L*(ranf_vb()-0.5); 
    	      rz[i] = L*(ranf_vb()-0.5); 
#endif
	      if (type==4)
		{
		  ox = 1;
		  oy = 0;
		  oz = 0;
		}
	      else if (type==1)
		orient_onsager(&ox, &oy, &oz, alpha);
	      else
		orient(&ox, &oy, &oz);
    	      versor_to_R(ox, oy, oz, Rl);
    	      for (k1=0; k1 < 3; k1++)
    		for (k2=0; k2 < 3; k2++)
    		  R[i][k1][k2] = Rl[k1][k2];
	    }
	  else
	    {
	      bt = 0;
	      while (1)
		{
		  nb = (int)(ranf_vb()*2.0);
		  j = ((int) (ranf_vb()*(i-size1)))+size1;

#if 1
	    	  if (bt > Oparams.parnum*1000)
	    	    {
		      printf("insertion cluster #2: maximum number of iteration reached\n");
		      if (fabs(calcpotene()+(i-size1)*2.0) < 1E-10)
			{
			  save_conf_mc(0);
			  printf("During insertion of #2 cluster: every particle has 2 bonds! i=%d\n", i);
			  exit(-1);
			}
	    	    }
#endif
		  if (is_bonded_mc(j, nb))
		    continue;
		  else
		    break;
		  bt++;
		}
	      /* mette la particella i legata a j con posizione ed orientazione a caso */
	      //printf("i=%d j=%d size1=%d size2=%d\n", i, j, size1, size2);
	      mcin(i, j, nb, type, alpha, &merr, 0);
	      if (merr!=0)
		{
		  save_conf_mc(0);
		  printf("[mcin] attemp to add a bonded particle failed!\n");
		  break;
		}
	      if (check_self_overlap(size1, i))
		{
		  selfoverlap = 1;
		  break;
		}

	      /* N.B. per ora non controlla il self-overlap della catena 
		 e la formazione dopo mcin di legami multipli poiché
		 si presuppone che al massimo stiamo considerando dimeri */
	    }
	  if (selfoverlap||merr)
	    {
	      tt++;
	      continue;
	    }
	  overlap = 0;
	  for (j=0; j < size1; j++)
	    {
#ifdef MD_LXYZ
	      shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
	      shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
	      shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
	      shift[0] = L*rint((rx[i]-rx[j])/L);
	      shift[1] = L*rint((ry[i]-ry[j])/L);
	      shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
	      ierr=0;
	      if (check_overlap(i, j, shift, &ierr)<0.0)
		{
		  overlap=1;
		  //printf("i=%d j=%d overlap!\n", i, j);
		  //printf("shift=%f %f %f\n", shift[0], shift[1], shift[2]);
		  //printf("r=%f %f %f  j %f %f %f\n", rx[i], ry[i], rz[i], rx[j], ry[j], rz[j]);
		  break;
		}
	    }
	  if (overlap)
	    break;
	}
      if (merr)
	{
	  tt++;
	  continue;
	}
      if (overlap)// && ierr==0)
	totene += 1.0;

      if (ierr!=0)
	{
	  printf("COV main loop NR failure\n");
	  printf("i=%d j=%d scalp=%.15G\n", i, j, R[i][0][0]*R[j][0][0]+R[i][0][1]*R[j][0][1]+R[i][0][2]*R[j][0][2]);
	  store_bump(i,j);
	}
      tt++;
   }
  if (type==1)
    {
      double norm, dtheta, pi;
      FILE *F;
      pi=2.0*acos(0.0);
      norm=0.0;
      dtheta = pi/nfons;
      for (i=0; i < nfons; i++)
	norm += distro[i]*2*pi*dtheta;
      for (i=0; i < nfons; i++)
	distro[i]/= norm;

      F=fopen("fons.dat", "w");  
      for (i=0; i < nfons; i++)
	fprintf(F, "%.15G %.15G\n",(i+0.5)*dtheta,distro[i]); 
      fclose(F);

      dtheta=pi/1000;
      F=fopen("fonsExact.dat", "w");  
      for (i=0; i < 1000; i++)
	fprintf(F, "%.15G %.15G\n",i*dtheta,sin(i*dtheta)*fons(i*dtheta,alpha)); 
      fclose(F);
    }

  //fclose(f);	
#ifdef MD_LXYZ
  printf("co-volume=%.10f (totene=%f)\n", (totene/((double)tt))*(L[0]*L[1]*L[2]), totene);
  printf("%.15G\n",(totene/((double)tt))*(L[0]*L[1]*L[2]));
#else 
  printf("co-volume=%.10f (totene=%f)\n", (totene/((double)tt))*(L*L*L), totene);
  printf("%.15G\n",(totene/((double)tt))*(L*L*L));
#endif
  exit(-1);
}
#endif
#if 0
double calcfLab(int i, double *x, double *rA, double **Ri);
int check_overlp_in_calcdist(double *x, double *fx, double *gx, int iA, int iB)
{
  double rA[3], rB[3], rmid[3], shift[3];
  int kk;
  double A, B;
  return 0;
  rA[0] = rx[iA];
  rA[1] = ry[iA];
  rA[2] = rz[iA];
#if 0
  if (!OprogStatus.dist5)
    rB=&(x[3]);
#endif
#ifdef MD_LXYZ
  shift[0] = L[0]*rint((rx[iA]-rx[iB])/L[0]);
  shift[1] = L[1]*rint((ry[iA]-ry[iB])/L[1]);
  shift[2] = L[2]*rint((rz[iA]-rz[iB])/L[2]);
#else
  shift[0] = L*rint((rx[iA]-rx[iB])/L);
  shift[1] = L*rint((ry[iA]-ry[iB])/L);
  shift[2] = L*rint((rz[iA]-rz[iB])/L);
#endif
  rB[0] = rx[iB]+shift[0];
  rB[1] = ry[iB]+shift[1];
  rB[2] = rz[iB]+shift[2];

  for (kk=0; kk < 3; kk++)
    {
      if (OprogStatus.dist5)
	rmid[kk] = x[kk] - x[4]*fx[kk]*0.5;
      else
	rmid[kk] = x[kk] - x[7]*fx[kk]*0.5;
    }
  //printf("rmid=%f %f %f rA=%f %f %f x[]=%.15G %.15G\n", rmid[0], rmid[1], rmid[2], rA[0], rA[1], rA[2],x[6], x[7]);
  //printf("iA=%d iB=%d\n", iA, iB);
  A=calcfLab(iA, rmid, rA, R[iA]);
  B=calcfLab(iB, rmid, rB, R[iB]);
  
  if (A<=0 && B<=0)
    {
      printf("rmid=%f %f %f rA=%f %f %f x[]=%.15G %.15G\n", rmid[0], rmid[1], rmid[2], rA[0], rA[1], rA[2],x[6], x[7]);
      printf("beccato! A=%.15G B=%.15G\n", A, B);
      return 1;
    }
  else
    return 0;
}
#endif
int mcmotion(void)
{
  int ip, dorej, movetype, err;
  double enn, eno;
  ip = Oparams.parnum*ranf();
  /* qui basta calcolare l'energia della particella che sto muovendo */
#if 0
  eno = calcpotene();
#else
  eno = calcpotene_GC(ip);
#endif
  store_coord(ip);
  movetype=random_move(ip);
  pbc(ip);
  update_LL(ip);
  //rebuildLinkedList();
  //printf("i=%d\n", i);
  totmovesMC++;
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
#if 0
      enn=calcpotene();
#else
      enn=calcpotene_GC(ip);
#endif
      if (enn <= eno)
	{
	  //	  if (abs(enn-eno) >=1 )
	  //	    printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	  dorej=0;
	}
      else
	{
	  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
	    dorej=0;
	  else
	    dorej=2;
	  // if (dorej==0)
	  // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
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
	  printf("[random_move] NR failed...I rejected this trial move...\n");
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
      //printf("restoring finished\n");
      //rebuildLinkedList();
    }
  if (OprogStatus.useNNL && dorej==0 )
    {
      overestimate_of_displ[ip] += displMC; 
    }
  return movetype;
}
void move(void)
{
  double acceptance, traaccept, ene, eno, rotaccept, volaccept=0.0;
#ifdef MD_LXYZ
  double avL;
#endif
  int ran, movetype, i, ip, err=0, dorej, enn, ntot, deln;
  /* 1 passo monte carlo = num. particelle tentativi di mossa */
  //printf("Doing MC step #%d\n", Oparams.curStep);
#if 1 
#ifdef MC_CALC_COVADD
  calc_cov_additive();
#endif
  if (OprogStatus.useNNL && do_nnl_rebuild())
    {
      //printf("Step #%d Rebuilding NNL\n", Oparams.curStep);
      for (i=0; i < Oparams.parnum; i++)
	overestimate_of_displ[i]=0.0;
      rebuildNNL();
      OprogStatus.lastNNLrebuildMC = Oparams.curStep;
    }
#endif
  if (cxini==-1)
    {
      cxini=cellsx;
      cyini=cellsy;
      czini=cellsz;
    }
  if (OprogStatus.nvbbias > 0)
    deln=OprogStatus.nvbbias;
  else
    deln=0;
  //printf("deln=%d\n", deln);
  if (OprogStatus.ensembleMC==0)
    ntot = Oparams.parnum+deln;
  else if (OprogStatus.ensembleMC==1)
    ntot = Oparams.parnum+1+deln;
#ifdef MC_GRANDCAN
  else if (OprogStatus.ensembleMC==2)
    ntot = OprogStatus.npav+OprogStatus.nexc+deln;
#endif
  else
    ntot=Oparams.parnum;
  for (i=0; i < ntot; i++)
    {
      if (OprogStatus.ensembleMC==0 && deln==0)
	ran = 0;
      else
	ran = ntot*ranf();

#ifdef MC_GRANDCAN
      if (OprogStatus.ensembleMC==2 && ran >= OprogStatus.npav)
	{
	  if (ran >= OprogStatus.npav + OprogStatus.nexc)
	    {
	      if (ranf() < OprogStatus.pbias)
		mcoutin(1.0/Oparams.T,OprogStatus.pbias);
	      else
		mcinout(1.0/Oparams.T,OprogStatus.pbias);
	    }
	  else
	    {
	      mcexc(&err);
	      movetype=4; /* 4 = insert/remove particle */
	      if(err)
		{
		  printf("[mcexc] NR failed...I rejected this trial move...\n");
		  err=0;
		}
	      excmoveMC++;
	    }
	}
      else
#endif
      if (OprogStatus.ensembleMC==1 && ran>=Oparams.parnum)
	{
	  //err=0;
	  if (ran > Oparams.parnum)
	    {
	      if (ranf() < OprogStatus.pbias)
		mcoutin(1.0/Oparams.T,OprogStatus.pbias);
	      else
		mcinout(1.0/Oparams.T,OprogStatus.pbias);
	    }
	  else
	    {
	      move_box(&err);
	      movetype=3; /* 0 = tra; 1 = rot 2 = tra and rot; 3 = box */
	      if(err)
		{
		  printf("[move_box] NR failed...I rejected this trial move...\n");
		  err=0;
		}
	      volmoveMC++;
	    }
	}
      else
	{
	  if (OprogStatus.ensembleMC==0 && ran >= Oparams.parnum)
	    {
	      if (ranf() < OprogStatus.pbias)
		mcoutin(1.0/Oparams.T,OprogStatus.pbias);
	      else
		mcinout(1.0/Oparams.T,OprogStatus.pbias);
	    }
	  else
	    movetype=mcmotion();
	}
      //printf("done\n");
    }
  if (OprogStatus.targetAccept > 0.0 && Oparams.curStep % OprogStatus.resetaccept==0)
    {
      acceptance=((double)(totmovesMC-totrejMC))/totmovesMC;
      traaccept = ((double)(tramoveMC-trarejMC))/tramoveMC;
      rotaccept = ((double)(rotmoveMC-rotrejMC))/rotmoveMC; 
 
      if (traaccept > OprogStatus.targetAccept)
	OprogStatus.deltaMC *= 1.1;
      else
	OprogStatus.deltaMC /= 1.1;
      if (OprogStatus.restrmove==0)
	{
	  if (rotaccept > OprogStatus.targetAccept)
	    OprogStatus.dthetaMC *= 1.1;
	  else
	    OprogStatus.dthetaMC /= 1.1;
	}
#ifdef MD_LXYZ
      if (OprogStatus.deltaMC > (avL=pow(L[0]*L[1]*L[2],1.0/3.0))*0.1)
	OprogStatus.deltaMC = avL*0.1;
#else
      if (OprogStatus.deltaMC > L*0.1)
	OprogStatus.deltaMC = L*0.1;
#endif
      if (OprogStatus.dthetaMC > 3.14)
	OprogStatus.dthetaMC = 3.14;
      /*NOTA: questo credo che si possa eliminare */
      if (OprogStatus.useNNL)
	calc_all_maxsteps();
    }
  if (OprogStatus.targetAcceptVol > 0.0 && OprogStatus.ensembleMC==1 && volmoveMC > 0 
      && (Oparams.curStep % OprogStatus.resetacceptVol==0))
    {
      volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
      //printf("sono qui volaccept=%.15G\n", volaccept);
      if (volaccept > OprogStatus.targetAcceptVol)
	OprogStatus.vmax *= 1.1;
      else
	OprogStatus.vmax /= 1.1;
    }
  if (Oparams.curStep%OprogStatus.outMC==0)
    {
      if (OprogStatus.ensembleMC==1 && volmoveMC > 0)
	volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
      //totmoves=((long long int)Oparams.parnum*(long long int)Oparams.curStep);
      acceptance=((double)(totmovesMC-totrejMC))/totmovesMC;
      traaccept = ((double)(tramoveMC-trarejMC))/tramoveMC;
      if (rotmoveMC!=0)
	rotaccept = ((double)(rotmoveMC-rotrejMC))/rotmoveMC; 
      else
	rotaccept = -1.0;
      printf("MC Step #%d pressure=%f temperature=%f Npar=%d\n", Oparams.curStep, Oparams.P, Oparams.T, Oparams.parnum);
      printf("Acceptance=%.15G (tra=%.15G rot=%.15G) deltaMC=%.15G dthetaMC=%.15G\n", acceptance, traaccept, 
	     rotaccept, OprogStatus.deltaMC, OprogStatus.dthetaMC);
      printf("rotmoveMC:%lld rotrefMC: %lld cells= %d %d %d\n", rotmoveMC, rotrejMC, cellsx, cellsy, cellsz);
      if (OprogStatus.ensembleMC==1 && volmoveMC > 0)
	printf("Volume moves acceptance = %.15G vmax = %.15G\n", volaccept, OprogStatus.vmax);
    }
  if ((Oparams.curStep % OprogStatus.resetacceptVol == 0) && OprogStatus.ensembleMC==1)
    volmoveMC=volrejMC=0;
  if (Oparams.curStep % OprogStatus.resetaccept==0)
    {
      totmovesMC=totrejMC=0;
      tramoveMC=trarejMC=0;
      rotmoveMC=rotrejMC=0;
    }
  if (Oparams.curStep==Oparams.totStep)
    {
      if (OprogStatus.targetAccept > 0.0)
	{
	  printf("deltrafin %.15G\n", OprogStatus.deltaMC);
	  printf("delrotfin %.15G\n", OprogStatus.dthetaMC);
	}
      if (OprogStatus.targetAcceptVol > 0.0 && OprogStatus.ensembleMC==1)
	printf("delvolfin %.15G\n", OprogStatus.vmax);
    }
}
#endif
