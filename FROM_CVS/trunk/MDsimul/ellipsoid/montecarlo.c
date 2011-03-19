#include<mdsimul.h>
const double saxfactMC[3]={0.85,0.68,0.68};
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
#endif
extern void ProcessCellCrossingMLL(void);
extern void PredictEventMLL(int na, int nb);
extern void PredictEventMLL_NLL(void);

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
extern void newtDistNeg(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3],int);
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
const int nmboxMC=2;

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
double check_overlap(int i, int j, double shift[3], int *errchk)
{
  int k, k1, k2;
  double vecg[8], vecgNeg[8], rA[3], rB[3], saxi[3], saxj[3];
  double d, d0, r1[3], r2[3], alpha; 
  OprogStatus.optnnl = 0;
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
     saxi[k] = 1.01* typesArr[typeOfPart[i]].sax[k];
     saxj[k] = 1.01* typesArr[typeOfPart[j]].sax[k];
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
      saxi[k] = saxfactMC[k]*typesArr[typeOfPart[i]].sax[k];
      saxj[k] = saxfactMC[k]*typesArr[typeOfPart[j]].sax[k];
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
  if (*errchk)
    {
      d0=overlap_using_multibox(i, j, shift);
      printf("I used multibox routine and d0=%f\n", d0);
      if (d0 < 0)
	{
	  *errchk=0;
	  return d0;
	}
      return -1.0;
    }
  //printf("QUI d=%f\n", d);
  return d;
}
extern void find_bonds_one_NLL(int i);

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
void remove_par_GC(int ip)
{
  int lp, k, k1, k2;
  lp = Oparams.parnum-1;
  /* copy all data from particle lp to ip */
  remove_bonds_GC(ip);
  //remove_from_current_cell(ip);
  rx[ip] = rx[lp];
  ry[ip] = ry[lp];
  rz[ip] = rz[lp];
  typeOfPart[ip] = typeOfPart[lp];
  is_a_sphere_NNL[ip] = is_a_sphere_NNL[lp];
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
  typeNP[0]++;
  Oparams.parnum--;
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
  for (kk=0; kk < numbonds[ip]; kk++)
    {
      jj = bonds[ip][kk] / (NANA);
      jj2 = bonds[ip][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      remove_bond(jj2,bb,jj,aa);
    }
}
extern void BuildNNL(int na);
extern void nextNNLupdate(int na);

void build_one_nnl_GC(int ip)
{
  int p, k1;
  nextNNLupdate(ip);
  BuildNNL(ip);
  /* and now add ip to neighbor lists of exisisting particles */
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
  orient(&ox,&oy,&oz);
  versor_to_R(ox, oy, oz, Rl);
  /* per ora assumiamo un solo tipo di particelle nel GC */
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
 
  /* grand canonical ensemble */
  if (ranf() < 0.5)
    {
      if (Oparams.parnum==0)
	return;
      o = Oparams.parnum*ranf();
      eno=calcpotene_GC(o);
      arg = Oparams.parnum*exp((1.0/Oparams.T)*eno)/(OprogStatus.zetaMC*vol);
      //printf("arg=%.15G zetaMC=%.15G \n", arg, OprogStatus.zetaMC);
      if (ranf() < arg)
	{
	  //printf("removing #%d\n", o);
	  remove_par_GC(o);
	}
    }
  else
    {
      /* nella seguente routine deve aggiungere la particella nelle LL e nelle NNL se usate */
      //printf("Inserting #%d\n", Oparams.parnum);
      np=insert_particle_GC();
      if (overlapMC(np, ierr))
	{
	  /* reject insertion */
	  //remove_from_current_cell(Oparams.parnum-1);
	  //printf("overlap Insertion rejected #%d\n", Oparams.parnum-1);
#ifdef MCGC_OPTLLREBUILD
	  remove_from_current_cell(np);
#endif
	  Oparams.parnum--;
#ifdef MCGC_OPTLLREBUILD
	  //printf("Celllist[parnum]=%d\n", cellList[Oparams.parnum+1]);
	  adjLinkedListRemove();
	  //printf("Celllist[parnum+1]=%d\n", cellList[Oparams.parnum]);
#else
	  rebuildLinkedList();
#endif
	  /* rimuove ip da tutte le NNL del sistema */
	  if (OprogStatus.useNNL)
	    remove_from_nnl_MC(np);
	  excrejMC++;
	  return;
	}	
      find_bonds_GC(np);
      enn=calcpotene_GC(np);
      //printf("enn=%.15G\n", enn);
      arg = OprogStatus.zetaMC*vol*exp(-(1.0/Oparams.T)*enn)/Oparams.parnum;
      if (ranf() >= arg)
	{
	  //printf("Insertion rejected #%d\n", np);
	  /* insertion rejected */
	  //remove_from_current_cell(Oparams.parnum-1);
	  remove_bonds_GC(np);
#ifdef MCGC_OPTLLREBUILD
	  remove_from_current_cell(np);
#endif
	  Oparams.parnum--;
#ifdef MCGC_OPTLLREBUILD
	  adjLinkedListRemove();
#else
	  rebuildLinkedList();
#endif
	  /* rimuove ip da tutte le NNL del sistema */
	  if (OprogStatus.useNNL)
	    remove_from_nnl_MC(np);
	  excrejMC++;
	  return;
	}
    }
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

extern int nbondsFlex;
extern int *mapbondsaFlex, *mapbondsbFlex; 
extern double calcDistNegSP(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, double *dists, int bondpair);

void find_bonds_covadd(int i, int j)
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
}
void mcin(int i, int j, int nb)
{
  double rA[3], rat[3], norm, sax, cc[3], ene;
  double shift[3], Rl[3][3], vv[3];
  double ox, oy, oz, d;
  int ierr, bonded, k1, k2, trials;
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
  bonded=0;
  trials=0;

  numbonds[i]=0;
  while (!bonded)
    {
#if 0
      /* chose a random position inside a cube*/
      rx[i] = cc[0]+2.0*sax*(ranf_vb()-0.5);
      ry[i] = cc[1]+2.0*sax*(ranf_vb()-0.5);
      rz[i] = cc[2]+2.0*sax*(ranf_vb()-0.5);
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
      orient(&ox, &oy, &oz);
      versor_to_R(ox, oy, oz, Rl);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      R[i][k1][k2] = Rl[k1][k2]; 
	    }
	}
      store_bonds_mc(j);
      find_bonds_covadd(i, j);
      if ((ene=calcpotene_GC(i)) < 0)
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
	      bonded=1;
	    }
	  else
	    {
	      restore_bonds_mc(j);
	      numbonds[i]=0;
	    }
	 // printf("d=%.15G trials #%d\n", d, trials);
	}

      //printf("i=%d j=%d ene=%f\n", i, j, ene);
      trials++;
    }
  //printf("trials=%d\n", trials);
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
void calc_persistence_length_mc(int maxtrials)
{
  int i, j, nb, k1, k2, tt;
  /* first particle is always in the center of the box with the same orientation */
  printf("calculating persistence length\n");

  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
      }


  for (tt=0; tt < maxtrials; tt++)
    {
      if (tt%100==0)
	{
	  printf("tt=%d\n", tt); 
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
	  mcin(i, j, nb);
	}
      save_conf_mc(tt);
    }
}
void calc_cov_additive(void)
{
  FILE *fi;
  double Lb, totene = 0.0, alpha, shift[3];
  int i, j, size1, size2, nb, tt, k1, k2, overlap=0, ierr, type;
  long long int maxtrials;
  double ox, oy, oz, Rl[3][3];
  fi = fopen("covmc.conf", "r");
  fscanf(fi, "%lld %d %d ", &maxtrials, &type, &size1);
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i]=0;
    }
  printf("ene iniziale=%f\n", calcpotene());
  if (type==1)
    fscanf(fi, " %lf ", &alpha);
  /* type = 0 -> covolume 
     type = 1 -> covolume nematic
     type = 2 -> persistence length
     type = 3 -> bonding volume
   */
  
  if (size1 >= Oparams.parnum)
    {
      printf("size1=%d must be less than parnum=%d\n", size1, Oparams.parnum);
      exit(-1);
    } 
  fclose(fi);
#ifdef MD_MAC
  srandomdev();
#else
  srandom((int)(time(NULL)));
#endif
  OprogStatus.optnnl = 0;
  tt=0;
  size2 = Oparams.parnum-size1;
  if (type==2)
    {
      calc_persistence_length_mc(maxtrials);
      exit(-1);
    }
  /* first particle is always in the center of the box with the same orientation */
  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
      }

  while (tt < maxtrials) 
    {
      /* place first cluster */
      if (tt%1000==0)
	{
	  printf("tt=%d\n", tt);
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
	  while (1)
	    {
	      nb = (int)(ranf_vb()*2.0);
	      j = (int) (ranf_vb()*i);
	      if (numbonds[j]==2)
		{
		  printf("j=%d has 2 bonds!\n", j);
		  exit(-1);
		}
	      if (is_bonded_mc(j, nb))
		continue;
	      else
		break;
	    }
	  mcin(i, j, nb);
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
    	      orient(&ox, &oy, &oz);
    	      versor_to_R(ox, oy, oz, Rl);
    	      for (k1=0; k1 < 3; k1++)
    		for (k2=0; k2 < 3; k2++)
    		  R[i][k1][k2] = Rl[k1][k2];
	    }
	  else
	    {
	      while (1)
		{
		  nb = (int)(ranf_vb()*2.0);
		  j = ((int) (ranf_vb()*(i-size1)))+size1;
	    	  if (numbonds[j]==2)
	    	    {
	    	      printf("j=%d has 2 bonds!\n", j);
	    	      exit(-1);
	    	    }
		  if (is_bonded_mc(j, nb))
		    continue;
		  else
		    break;
		}
	      /* mette la particella i legata a j con posizione ed orientazione a caso */
	      //printf("i=%d j=%d size1=%d size2=%d\n", i, j, size1, size2);
	      mcin(i, j, nb);
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
      if (overlap)
	totene += 1.0;
      tt++;
   }
#ifdef MD_LXYZ
  Lb = L[0];
#else
  Lb = L;
#endif
 
  printf("co-volume=%.10f (totene=%f)\n", (totene/((double)tt))*(Lb*Lb*Lb), totene);
  printf("%.15G\n",(totene/((double)tt))*(Lb*Lb*Lb));
}
#endif
void move(void)
{
  double acceptance, traaccept, ene, eno, rotaccept, volaccept=0.0;
#ifdef MD_LXYZ
  double avL;
#endif
  int ran, movetype, i,ip, err=0, dorej, enn;
  /* 1 passo monte carlo = num. particelle tentativi di mossa */
  //printf("Doing MC step #%d\n", Oparams.curStep);
#if 1 
#ifdef MC_CALC_COVADD
  calc_cov_additive();
  exit(-1);
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
  for (i=0; i < Oparams.parnum; i++)
    {
      /* pick a particle at random */
#ifdef MC_GRANDCAN
      if (OprogStatus.ensembleMC==2) /* grand canonical */
	ran = (OprogStatus.npav+OprogStatus.nexc)*ranf();
      else
#endif
      if (OprogStatus.ensembleMC==1)
	ran=(Oparams.parnum+1)*ranf();
      else 
	ran = 0;
#ifdef MC_GRANDCAN
      if (OprogStatus.ensembleMC==2 && ran >= OprogStatus.npav)
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
      else
#endif
      if (OprogStatus.ensembleMC==1 && ran==Oparams.parnum)
	{
	  //err=0;
	  move_box(&err);
	  movetype=3; /* 0 = tra; 1 = rot 2 = tra and rot; 3 = box */
    	  if(err)
    	    {
	      printf("[move_box] NR failed...I rejected this trial move...\n");
	      err=0;
	    }
	  volmoveMC++;
	}
      else
	{
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
	}
      //printf("done\n");
    }
  if (Oparams.curStep%OprogStatus.outMC==0)
    {
      //totmoves=((long long int)Oparams.parnum*(long long int)Oparams.curStep);
      acceptance=((double)(totmovesMC-totrejMC))/totmovesMC;
      traaccept = ((double)(tramoveMC-trarejMC))/tramoveMC;
      rotaccept = ((double)(rotmoveMC-rotrejMC))/rotmoveMC; 
      if (OprogStatus.ensembleMC==1 && volmoveMC > 0)
	volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
      if (OprogStatus.targetAccept > 0.0 && Oparams.curStep % OprogStatus.resetaccept==0)
	{
	  if (traaccept > OprogStatus.targetAccept)
	    OprogStatus.deltaMC *= 1.1;
	  else
	    OprogStatus.deltaMC /= 1.1;
	  if (rotaccept > OprogStatus.targetAccept)
	    OprogStatus.dthetaMC *= 1.1;
	  else
	    OprogStatus.dthetaMC /= 1.1;
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
      if (OprogStatus.targetAcceptVol > 0.0 && OprogStatus.ensembleMC==1
	  && (Oparams.curStep % OprogStatus.resetacceptVol==0))
	{
	  //printf("sono qui volaccept=%.15G\n", volaccept);
	  if (volaccept > OprogStatus.targetAcceptVol)
	    OprogStatus.vmax *= 1.1;
	  else
	    OprogStatus.vmax /= 1.1;
	}
      printf("MC Step #%d pressure=%f temperature=%f Npar=%d\n", Oparams.curStep, Oparams.P, Oparams.T, Oparams.parnum);
      printf("Acceptance=%.15G (tra=%.15G rot=%.15G) deltaMC=%.15G dthetaMC=%.15G\n", acceptance, traaccept, 
	     rotaccept, OprogStatus.deltaMC, OprogStatus.dthetaMC);
      printf("rotmoveMC:%lld rotrefMC: %lld cells= %d %d %d\n", rotmoveMC, rotrejMC, cellsx, cellsy, cellsz);
      if (OprogStatus.ensembleMC==1 && volmoveMC>0)
	printf("Volume moves acceptance = %.15G vmax = %.15G\n", volaccept, OprogStatus.vmax);
      if ((Oparams.curStep % OprogStatus.resetacceptVol == 0) && OprogStatus.ensembleMC==1)
	volmoveMC=volrejMC=0;
      if (Oparams.curStep % OprogStatus.resetaccept==0)
	{
	  totmovesMC=totrejMC=0;
	  tramoveMC=trarejMC=0;
	  rotmoveMC=rotrejMC=0;
	}

    }
}
#endif
