#ifdef MC_HC
#undef DEBUG_HCMC
#undef MC_HC_SPHERO_OPT
#include<mdsimul.h>
#include<complex.h>
extern const double saxfactMC[3];
#ifdef MC_QUASI_CUBE
extern const double saxfactMC_QC[3];
#endif
extern const int nfons;
extern void init_rng(int mdseed, int mpi, int my_rank);
#ifdef MC_SIMUL
#ifdef MC_STORELL
int *cellListMC;
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
extern double rA[3], rB[3];
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
extern double toteneini;
extern long long int ttini;
extern int covrestart;
extern const int nmboxMC;
extern double totdist, distcc;
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
  if (test < 1.0E-14)
    {
      //printf("convergence reached! test=%.15G\n", test);
      return 1;
    }
  else 
    return 0;
}
double totitsHC = 0.0;
double numcallsHC = 0.0;
#ifdef DEBUG_HCMC
extern int dostorebump;
#endif
extern double calc_norm(double *vec);
extern void versor_to_R(double ox, double oy, double oz, double R[3][3]);
void body2labHC(int i, double xp[3], double x[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      /* NOTE: k2 starts from 1 because xp[0] = 0.0 see function find_initial_guess() below */
      for (k2=1; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
      x[k1] += rO[k1];
    }
}
#ifdef HC_ALGO_OPT
#define MESH_PTS 10000
double meshptsGbl;
struct brentOpt 
{
  double th; 
  double UipPjp[3];
  double dUipPjp[3];
  double normUipPjp;
  double normUipPjpSq;
  double Uip[3];
  double Pjp[3];
  double Pip[3];
  double PjCi[3];
  double PjPi[3];
  double lambda;
  double sinth;
  double costh;
  double D2sinth;
  double D2costh;
  double minPgbl[3];
  int id; /* 0 = not yet calculated; 1 = calculated by drimdisk; 2 = calculated by rimdisk */
} brentmsg;

double CipGbl[3], nipGbl[3], Dgbl, minPgbl[3];
void calcrimdisk(double th);
#if 0
double find_initial_guess_bracket2(double *thg, int meshpts)
{
  static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3], dist, mindist;
  static double *tharr;
  static struct brentOpt *mesh;
  FILE* f;
  int k1, k2, nn;
  /* bracketing */
  if (firstcall || thg==NULL)
    {
      mesh = malloc(sizeof(struct brentOpt)*meshpts);
      firstcall=0;
      dth = 2.0*M_PI/((double)meshpts);

      th=0.0;
      for (nn=0; nn < meshpts; nn++)
	{
	  calcrimdisk(th);
	  mesh[nn].sinth = brentmsg.sinth;
	  mesh[nn].costh = brentmsg.costh;
	  mesh[nn].id = 3; 
	  for (k1=0; k1 < 3; k1++)
	    {
	      mesh[nn].Pjp[k1] = brentmsg.Pjp[k1];
	    } 
	  th += dth;
	}
      if (thg==NULL)
	return -1.0;
    }
  dth = 2.0*M_PI/meshpts;
  th = 0;
  mindist = -1;
  if (meshpts > 256)
    f=fopen("dtheta.dat", "w+");
  for (k1 = 0; k1 < meshpts; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	brentmsg.PjCi[k2] = mesh[k1].Pjp[k2] - CipGbl[k2]; 
      brentmsg.lambda = scalProd(brentmsg.PjCi,nipGbl);
      for (k2=0; k2 < 3; k2++)
	{
	  brentmsg.Uip[k2] = CipGbl[k2] + brentmsg.lambda*nipGbl[k2];
	  brentmsg.UipPjp[k2] = brentmsg.Uip[k2] - mesh[k1].Pjp[k2];
	  //PjPi[k1] = Pjp[k1] - (CipGbl[k1] + lambda*nipGbl[k1]);
	}
      dist = calc_norm(brentmsg.UipPjp);
      if (meshpts > 256)
	fprintf(f, "%.15G %.15G\n", th, dist);
      //printf("bracket 2 (%f,%d)=%.15G\n", th, k1, dist);
      if (k1==0 || dist < mindist)
	{
	  mindist = dist;
	  *thg = th;
	}
      th+=dth;
    }
  if (meshpts > 256)
  fclose(f);  
  return mindist;
}
#endif
struct maxminS { 
  double max;
  double min;
  int nnmax;
  int nnmin;
  double thgmax;
  double thgmin;
} maxmin;

struct iniguessStruct { 
  int num;
  int multmin;
  double Pg[3][3];
  double thg[3];
  double dist[3];
  int which;
} iniguess;
int check_convergence_vec(double Told[3], double Tnew[3], double tol)
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
  if (test < tol)
    {
      //printf("convergence reached! test=%.15G\n", test);
      return 1;
    }
  else 
    return 0;
}

double rimdiskfunc(double th);
double drimdiskfunc(double th);

void find_initial_linecircle(double diam, double Cip[3], double nip[3], struct iniguessStruct *ig)
{
  //static int firstcall=1;
  double m, q, xp[3], Ui[3], UiPj[3], dist, maxdist;
  double normCipP, CipP[3], ab2, asq1, nipP[3], PgCip[3], sqrtdelta, PgOld[3];
  double Pp[3];
  //static double *tharr;
  //static struct brentOpt *mesh;
  double sp, D2, lambda, aR, bR, delta, norm;
  int k1, k2, nn, numsol=0, iter;
  /* bracketing */
  maxdist = -1;

  normCipP = sqrt(Sqr(Cip[1])+Sqr(Cip[2])); 
  D2 = diam*0.5;
  CipP[0]= 0.0;
  CipP[1]=Cip[1];
  CipP[2]=Cip[2];
  nipP[0]= 0.0;
  nipP[1]=nip[1];
  nipP[2]=nip[2];
 
  /* cerco l'intersezione della retta proiettata sul piano
   * del disco con la circonferenza */
  ig->Pg[0][0]=ig->Pg[1][0]=0.0;

  if (nipGbl[1]==0.0)
    {
      if (fabs(CipP[1]) == D2)
	{
	  ig->num = 1;
	  ig->Pg[0][1] = CipP[1];
	  ig->Pg[0][2]=0.0;
	}
      else if (fabs(CipP[1]) < D2)
	{
	  numsol = 2;
	  ig->Pg[0][1]=ig->Pg[1][1]=CipP[1];
	  ig->Pg[0][2]= sqrt(Sqr(D2)-Sqr(Cip[1]));
	  ig->Pg[1][2]= -ig->Pg[0][2];
	}
      else 
	{
	  numsol = 1;
	  ig->Pg[0][1] = CipP[1];
	  ig->Pg[0][2]=0.0;
	}
    }
  else 
    {
      m = nipP[2]/nipP[1];
      q = CipP[2]-m*CipP[1];
      
      delta = 4.0*(Sqr(m*q) - (Sqr(m)+1.0)*(Sqr(q)-Sqr(D2)));
#if 1
      printf("m=%.15G q=%.15G delta=%.15G\n", m, q, delta);
      printf("Ci=%f %f %f ni=%f %f %f\n", Cip[0], Cip[1], Cip[2], nip[0], nip[1], nip[2]); 
#endif
      if (delta==0)
	{
	  ig->num = 1;
	  ig->Pg[0][1] = -m*q/(Sqr(m)+1.0);  
	  ig->Pg[0][2] = m*ig->Pg[0][1]+q;
	}
      else if (delta > 0.)
	{
	  printf("qui>>\n");
	  ig->num = 2;
	  ab2 = -2.0*m*q;
	  asq1 = 2.0*(Sqr(m)+1);
	  sqrtdelta = sqrt(delta);
	  ig->Pg[0][1] = (ab2 + sqrtdelta)/asq1; 
	  ig->Pg[0][2] = m*ig->Pg[0][1]+q;
	  ig->Pg[1][1] = (ab2 - sqrtdelta)/asq1;
	  ig->Pg[1][2] = m*ig->Pg[1][1]+q;
	}
      else 
	{
	  /* prende il punto sulla circonferenza che minimizza la distanza 
	   con la retta proiettata sul piano del disco */
	  ig->num = 1;
	  sp = scalProd(CipP,nipP); 
	  for (k1=0; k1 < 3; k1++)
	    ig->Pg[0][k1] = CipP[k1] - sp*nipP[k1];
	  norm = calc_norm(ig->Pg[0]);
	  for (k1=0; k1 < 3; k1++)
	   ig->Pg[0][k1] = D2*ig->Pg[0][k1]/norm; 
	}
    }
  /* further refine the guessed points */
#if 1
  for (k2 = 0; k2 < ig->num; k2++)
    {
      printf("Prima sol # %d %.15G %15G %.15G\n", k2, ig->Pg[k2][0], ig->Pg[k2][1],ig->Pg[k2][2]);
      printf("norm=%f\n", calc_norm(ig->Pg[k2]));
    }
#endif 
 #if 1
  for (k2 = 0; k2 < ig->num; k2++)
    { 
#if 0
      for (k1=0; k1 < 3; k1++)
	{
	  //PgOld[k1] = ig->Pg[k2][k1];
	  //PgCip[k1] = ig->Pg[k2][k1];
	  PgCip[k1] = ig->Pg[k2][k1] - Cip[k1];
	}
      sp = scalProd(PgCip,nip);
      for (k1=0; k1 < 3; k1++)
	{
	  PgCip[k1] = Cip[k1] + sp*nip[k1]; 
	}
      PgCip[0] = 0.0;
      norm = calc_norm(PgCip);

      for (k1=0; k1 < 3; k1++)
	ig->Pg[k2][k1] = D2*PgCip[k1]/norm;	
      if (check_convergence_vec(PgOld, ig->Pg[k2]))
	break;
#endif
      if (ig->Pg[k2][1]>=D2)
	ig->thg[k2] = 0.0;
      else if (ig->Pg[k2][1] <= -D2)
	ig->thg[k2] = M_PI;
      else if (ig->Pg[k2][2] >= 0.0)
	ig->thg[k2] = acos(ig->Pg[k2][1]);
      else
	ig->thg[k2] = 2.0*M_PI-acos(ig->Pg[k2][1]);
      //printf("thg[%d]=%.15G\n", k2, ig->thg[k2]);
      //iter++;
    }
#endif
  if (ig->num==2)
    ig->multmin = 1;
  else
    ig->multmin = 0;
  /* further refine the guessed points */
  brentmsg.id = 0;
  ig->dist[0]=rimdiskfunc(ig->thg[0]);
  brentmsg.id = 0;
  ig->dist[1]=rimdiskfunc(ig->thg[1]);
  if (ig->num==1)
    {
      ig->which = 0;
    }
  else if (ig->dist[0] <= ig->dist[1])
    {
      //printf("1) qui?!? numsol=%d thg=%.15G\n", ig->num, ig->thg[0]);
      ig->which = 0;
    }
  else
    {
      //printf("2) qui?!? numsol=%d thg=%.15G\n", ig->num, ig->thg[1]);
      ig->which = 1;
    }

  /* infine considera l'intersezione della retta con il piano in cui giace il disco */
  if (nip[0]!=0.0)
    {
      lambda = -Cip[0]/nip[0];
      Pp[0]=  0.0;
      Pp[1] = Cip[1]+lambda*nip[1];
      Pp[2] = Cip[2]+lambda*nip[2];
      norm = sqrt(Sqr(Pp[1]) + Sqr(Pp[2]));      
      for (k1=0; k1 < 3; k1++)
	Pp[k1] = D2*Pp[k1]/norm;
#if MC_RIMDISK_NORMSQ
      dist = Sqr(norm-D2);      
#else
      dist = norm-D2;      
#endif
      if (dist < ig->dist[ig->which])
	{
    	  for (k1=0; k1 < 3; k1++)
	    ig->Pg[ig->num][k1] = Pp[k1];
	  ig->which = 0;
	  if (ig->Pg[ig->num][1]>=D2)
	    ig->thg[ig->num] = 0.0;
	  else if (ig->Pg[ig->num][1] <= -D2)
	    ig->thg[ig->num] = M_PI;
	  else if (ig->Pg[ig->num][2] >= 0.0)
	    ig->thg[ig->num] = acos(ig->Pg[ig->num][1]);
	  else
	    ig->thg[ig->num] = 2.0*M_PI-acos(ig->Pg[ig->num][1]);
	  ig->which = ig->num; 
	  ig->num++;
	  //brentmsg.id = 0;
	  //printf("3) qui?!? dfx=%15G\n", drimdiskfunc(ig->thg[0]));
	}
    }
  /* trova la distanza tra Pp e la circonferenza */
#if 0
  for (k2 = 0; k2 < ig->num; k2++)
    {
      printf("Dopo sol # %d %.15G %15G %.15G\n", k2, ig->Pg[k2][0], ig->Pg[k2][1],ig->Pg[k2][2]);
    } 
  printf("ig->num=%d\n", ig->num);
#endif
}
double find_initial_guess_proj(double diam, double Cip[3], double nip[3], struct iniguessStruct *ig)
{
  int k1, k2, iter;
  const int MAXNUMITER = 100;
  double PgOld[3], D2, sp, PgCip[3], norm, Cold[3];

  D2 = diam*0.5;
  for (k1=0; k1 < 3; k1++)
    {
      ig->Pg[0][k1] = 0.0;
      Cold[k1] = Cip[k1];
    }
  for (iter=0; iter < MAXNUMITER; iter++)
    {
      for (k1=0; k1 < 3; k1++)
	PgOld[k1] = ig->Pg[0][k1];
      for (k1=0; k1 < 3; k1++)
	{
	  //PgOld[k1] = ig->Pg[k2][k1];
	  //PgCip[k1] = ig->Pg[k2][k1];
	  PgCip[k1] = Cip[k1] - ig->Pg[0][k1];
	}
      sp = scalProd(PgCip,nip);
      for (k1=0; k1 < 3; k1++)
	{
	  PgCip[k1] = Cip[k1] - sp*nip[k1]; 
	}
      PgCip[0] = 0.0;
      norm = calc_norm(PgCip);

      Cold[k1] = PgCip[k1];
      
      for (k1=0; k1 < 3; k1++)
	ig->Pg[0][k1] = D2*PgCip[k1]/norm;	
      if (check_convergence_vec(PgOld, ig->Pg[0], 1.0E-14))
	break;
    }
  if (ig->Pg[0][1]>=D2)
    ig->thg[0] = 0.0;
  else if (ig->Pg[0][1] <= -D2)
    ig->thg[0] = M_PI;
  else if (ig->Pg[0][2] >= 0.0)
    ig->thg[0] = acos(ig->Pg[0][1]);
  else
    ig->thg[0] = 2.0*M_PI-acos(ig->Pg[0][1]);
  printf("iter=%d guessed=%.15G\n", iter, ig->thg[0]);
}
double find_initial_guess_bracket(double *thg, int meshpts, struct brentOpt* mesh, int *nnini, struct maxminS *maxmin)
{
  //static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3], dist, maxdist;
  //static double *tharr;
  //static struct brentOpt *mesh;
  int k1, k2, nn;
  /* bracketing */
  dth = 2.0*M_PI/meshpts;
  th = 0;
  maxdist = -1;
  for (k1 = 0; k1 < meshpts; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	mesh[k1].PjCi[k2] = mesh[k1].Pjp[k2] - CipGbl[k2]; 
      mesh[k1].lambda = scalProd(mesh[k1].PjCi,nipGbl);
      for (k2=0; k2 < 3; k2++)
	{
	  mesh[k1].Uip[k2] = CipGbl[k2] + mesh[k1].lambda*nipGbl[k2];
	  mesh[k1].UipPjp[k2] = mesh[k1].Uip[k2] - mesh[k1].Pjp[k2];
	  //PjPi[k1] = Pjp[k1] - (CipGbl[k1] + lambda*nipGbl[k1]);
	}
      dist = mesh[k1].normUipPjp = calc_norm(mesh[k1].UipPjp);
      //printf("bracket 1 (%f,%d)=%.15G\n", th, k1, dist);
      if (k1==0 || dist > maxmin->max)
	{
	  maxmin->max = dist;
	  *thg=maxmin->thgmax = th;
	  *nnini=maxmin->nnmax = k1;
	}
      if (k1==0 || dist < maxmin->min)
	{
	  maxmin->min = dist;
	  maxmin->thgmin = th;
	  maxmin->nnmin = k1;
	}
      th+=dth;
    }
  return maxmin->max;
}
double find_initial_guess_opt(double *Aj, double Ci[3], double ni[3], double Dj[3], double nj[3], double D, double *thmin)
{
  const int meshpts = MESH_PTS;
  double Pj[3], Rj[3][3], AiCi[3];
  int kk, k1, k2, nn;
  static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3];
  static double *tharr;
  static double **mesh; /* {{1,0},{0.707106781186547, 0.707106781186547},{0,1},
      {-0.707106781186547,0.707106781186547},{-1,0},{-0.707106781186547,-0.707106781186547},
      {0,-1},{0.707106781186547,-0.707106781186547}};*/
  double PjCini, PjCi[3], normPjCi, d, mindist=-1.0; 
  versor_to_R(nj[0],nj[1],nj[2], Rj); 
#if 1
  if (firstcall)
    {
      mesh = malloc(sizeof(double*)*meshpts);
      tharr =malloc(sizeof(double)*meshpts); 
      for (nn=0; nn < meshpts; nn++)
	mesh[nn] = malloc(sizeof(double)*3);
      firstcall=0;
      dth = acos(0)*4.0/((double)meshpts);

      th=0.0;
      for (nn=0; nn < meshpts; nn++)
	{
	  mesh[nn][0] = cos(th);
	  mesh[nn][1] = sin(th);
	  tharr[nn] = th;
	  th += dth;
	}
    }
#endif
  for (nn=0; nn < meshpts; nn++)
    {
      //xp[0] = 0.0;
      xp[1] = D*0.5*mesh[nn][0];
      xp[2] = D*0.5*mesh[nn][1];
      body2labHC(0, xp, Pj, Dj, Rj);    
//	printf("xp=%f\n", xp[0]);
      for (kk=0; kk < 3; kk++)
	PjCi[kk] = Pj[kk] - Ci[kk];
      //normPjCi = calc_norm(PjCi);
      PjCini = scalProd(PjCi,ni);
      for (kk=0; kk < 3; kk++)
	{
	  Ui[kk] = Ci[kk] + PjCini*ni[kk];
	  UiPj[kk] = Ui[kk]-Pj[kk];
	}
      if ((d=calc_norm(UiPj)) < mindist || nn==0)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      Aj[kk] = Pj[kk];
    	    }
	  *thmin = tharr[nn];
	  mindist=d;
	  //printf("nn=%d mindist=%.15G d=%.15G\n", nn, mindist, d);
	  //printf("Ui=%f %f %f Pi=%f %f %f\n", Ui[0],Ui[1], Ui[2], Pj[0], Pj[1], Pj[2]);
	}
    }
  //printf("mindist=%f thmin=%f\n", mindist, *thmin);
  return mindist;
  //printf("done\n");
#if 0
   for (kk=0; kk < 3; kk++)
     AiCi[kk]  = Ai[kk] - Ci[kk]; 
  printf("norm AiCi=%.15G sp=%.15G\n", calc_norm(AiCi), scalProd(AiCi,ni)/calc_norm(AiCi));
  for (kk=0; kk < 3; kk++)
    AiCi[kk]  = Pj[kk] - Dj[kk]; 

  printf("norm AiCi=%.15G sp=%.15G\n", calc_norm(AiCi), scalProd(AiCi,nj));
#endif 
}
#endif
void find_initial_guess(double *Ai, double Ci[3], double ni[3], double Dj[3], double nj[3], double D)
{
  const int meshpts = 8;
  double Pj[3], Rj[3][3], AiCi[3];
  int kk, k1, k2, nn;
  static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3];
  static double **mesh; /* {{1,0},{0.707106781186547, 0.707106781186547},{0,1},
      {-0.707106781186547,0.707106781186547},{-1,0},{-0.707106781186547,-0.707106781186547},
      {0,-1},{0.707106781186547,-0.707106781186547}};*/
  double PjCini, PjCi[3], normPjCi, d, mindist=-1.0; 
  versor_to_R(nj[0],nj[1],nj[2], Rj); 
#if 1
  if (firstcall)
    {
      mesh = malloc(sizeof(double*)*meshpts);
      for (nn=0; nn < meshpts; nn++)
	mesh[nn] = malloc(sizeof(double)*3);
      firstcall=0;
      dth = acos(0)*4.0/((double)meshpts);

      th=0.0;
      for (nn=0; nn < meshpts; nn++)
	{
	  mesh[nn][0] = cos(th);
	  mesh[nn][1] = sin(th);
	  th += dth;
	}
    }
#endif
  for (nn=0; nn < meshpts; nn++)
    {
      //xp[0] = 0.0;
      xp[1] = D*0.5*mesh[nn][0];
      xp[2] = D*0.5*mesh[nn][1];
      body2labHC(0, xp, Pj, Dj, Rj);    
      for (kk=0; kk < 3; kk++)
	PjCi[kk] = Pj[kk] - Ci[kk];
      //normPjCi = calc_norm(PjCi);
      PjCini = scalProd(PjCi,ni);
      for (kk=0; kk < 3; kk++)
	{
	  Ui[kk] = Ci[kk] + PjCini*ni[kk];
	  UiPj[kk] = Ui[kk]-Pj[kk];
	}
      if ((d=calc_norm(UiPj)) < mindist || nn==0)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      Ai[kk] = Ui[kk];
    	    }
	  mindist=d;
	  //printf("nn=%d mindist=%.15G d=%.15G\n", nn, mindist, d);
	  //printf("Ui=%f %f %f Pi=%f %f %f\n", Ui[0],Ui[1], Ui[2], Pj[0], Pj[1], Pj[2]);
	}
    }
  //printf("done\n");
#if 0
   for (kk=0; kk < 3; kk++)
     AiCi[kk]  = Ai[kk] - Ci[kk]; 
  printf("norm AiCi=%.15G sp=%.15G\n", calc_norm(AiCi), scalProd(AiCi,ni)/calc_norm(AiCi));
  for (kk=0; kk < 3; kk++)
    AiCi[kk]  = Pj[kk] - Dj[kk]; 

  printf("norm AiCi=%.15G sp=%.15G\n", calc_norm(AiCi), scalProd(AiCi,nj));
#endif 
}

#ifdef MC_HC_SPHERO_OPT
double check_spherocyl(double CiCj[3], double D, double L, double Di[2][3], double *Ci, double *ni, double Dj[2][3], double *Cj, double *nj, int *rim);
#endif
double calcDistNegHCdiff(int i, int j, double shift[3], int* retchk)
{
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int it, k2;
  double normNSq, ViVj[3], lambdai, lambdaj, Li, Diami, Lj, Diamj; 
  double LiTmp, LjTmp, DiamiTmp, DiamjTmp;
  double sp, Q1, Q2, normPiDi, normPjDj, normN, DiN, DjN, niN[3], njN[3], Djni, Djnj;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DiCi[3], DiCini, normDiCi, DjCi[3], normDjCi;
  double PiDi[3], PjDj[3], Ai[3], Tj[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3], TipCjnj, TimCjnj;
  double Aj[3], AjDini, AjDinivec[3], AjDi[3], Tip[3], Tim[3], TipCj[3], TimCj[3], Dini;
  double DiCj[3], normDiCj, DiCjnj, Uj[3], DiUj[3], normDiUj, DiUjnj;
  double Tim_perp[3], Tip_perp[3], Tim_para[3], Tip_para[3], normTim_perp, DjCini;
  double Tjm_perp[3], Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp;
  double TiOld[3], TiNew[3], TiNewCj[3], TiNewCjnj, Tjpara, Tjperp[3];	
  double normCiCj;	
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3];
  int kk, j1, j2;
  *retchk = 0; 

  // return calcDistNegHCsame(i, j, shift, retchk);
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
  Li = 2.0*typesArr[typeOfPart[i]].sax[0];
  Diami = 2.0*typesArr[typeOfPart[i]].sax[1];
  Lj = 2.0*typesArr[typeOfPart[j]].sax[0];
  Diamj = 2.0*typesArr[typeOfPart[j]].sax[1];

  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }

  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*Li*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*Li*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*Lj*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*Lj*nj[kk];
    }
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= 0.5*(Li+Lj))
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
	    if (sp == 0 && calc_norm(VV) < 0.5*(Diami+Diamj))
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
#ifdef DEBUG_HCMC
	    printf("Di=%f %f %f\n", Di[j1][0], Di[j1][1], Di[j1][2]);
	    printf("Dj=%f %f %f\n", Dj[j2][0], Dj[j2][1], Dj[j2][2]);
	    printf("normPiDi: %f normPjDj=%f\n", normPiDi, normPjDj);
	    printf("0.5*Diami=%f 0.5*Diamj=%f\n", 0.5*Diami, 0.5*Diamj);
#endif
	    if (normPiDi <= 0.5*Diami && normPjDj <= 0.5*Diamj)
	      {
		Q1 = sqrt(Sqr(Diami)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(Diamj)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  /* =================================== >>> Part A <<< ========================= */
  for (j1=0; j1 < 2; j1++)
    {

      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      DiamiTmp = Diami;
	      DiamjTmp = Diamj;
	      LiTmp = Li;
	      LjTmp = Lj;
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	      Diami = Diamj;
	      Diamj = DiamiTmp;
	      Li = Lj;
	      Lj = LiTmp;
	    }
	}
      for (j2=0; j2 < 2; j2++)
	{
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);

	  if (normDjUi > 0.5*(Diami+Diamj))
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Uj| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

	   */
	  if (normDjUi < Diami*0.5 && fabs(DjCini) > Li*0.5)
	    continue;

	  if (normDjUi < Diami*0.5 && fabs(DjCini) <= Li*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }
#if 1
	  find_initial_guess(Ai, Ci, ni, Dj[j2], nj, Diamj);
#else
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
	  for (it = 0; it < MAX_ITERATIONS; it++)
	    {
	      for (kk=0; kk < 3; kk++)
		{
		  AiDj[kk] = Ai[kk] - Dj[j2][kk];
		}
	      AiDjnj = scalProd(AiDj,nj);
	      vectProdVec(AiDj,nj,AiDjnjvec);
	      for (kk=0; kk < 3; kk++)
		VV[kk] =  0.5*Diamj*(AiDj[kk]-AiDjnj*nj[kk])/calc_norm(AiDjnjvec);
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
	      normTjm_perp = calc_norm(Tjp_perp);
	      for (kk=0; kk < 3; kk++)
		TjOld[kk] = TjNew[kk];
	      if (calc_norm(Tjm_perp) < calc_norm(Tjp_perp))
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

#ifdef DEBUG_HCMC
	      printf("j1=%d A it=%d Aiold=%.15G %.15G %.15G\n", j1, it, Ai[0], Ai[1], Ai[2]);
#endif
	      for (kk=0; kk < 3; kk++)
		Ai[kk] = TjNewCini*ni[kk] + Ci[kk]; 
#ifdef DEBUG_HCMC
	      printf("A it=%d Ainew=%.15G %.15G %.15G TjNewCini=%.15G\n", it, Ai[0], Ai[1], Ai[2], TjNewCini);
	      printf("A Ci=%.15G %.15G %.15G\n", Ci[0], Ci[1], Ci[2]);
	      printf("A ni=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
#endif
	      if ( it > 0 && check_convergence(TjOld,TjNew) ) 
		break;
	    }
	  totitsHC += it;
#ifdef DEBUG_HCMC
	  printf("A #1 number of iterations=%d Tjold=%.15G %.15G %.15G Tjnew=%.15G %.15G %.15G\n",it, 
		 TjOld[0], TjOld[1], TjOld[2], TjNew[0], TjNew[1], TjNew[2]);
#endif
	  if (it >= MAX_ITERATIONS)
	    {
	      printf("MAX ITERATIONS REACHED in A!\n");
	      *retchk=1;
	      return -1;
	    }
	  if ( (calc_norm(Tjp_para) <= Li*0.5 && calc_norm(Tjp_perp) <= Diami*0.5)||
	       (calc_norm(Tjm_para) <= Li*0.5 && calc_norm(Tjm_perp) <= Diami*0.5) )
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #2 disk-rim\n");
#endif	   
	      return -1;
	    }
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	      Diami = DiamiTmp;
	      Diamj = DiamjTmp;
	      Li = LiTmp;
	      Lj = LjTmp;
	    }
	}

    }
  /* =================================== >>> Part B <<< ========================= */
  numcallsHC += 4.0; 

  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < 0.5*(Diami+Diamj) && fabs(lambdai) < 0.5*Li && fabs(lambdaj) < 0.5*Lj)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 1;
}
#ifdef HC_ALGO_OPT
extern double zbrent(double (*func)(double), double x1, double x2, double tol);
extern double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
extern double dbrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double), double tol, double *xmin);
double signGbl;
void calcrimdisk(double th)
{
  int k1, k2;
  double D2;

  D2 = 0.5*Dgbl;
  brentmsg.th = th;
  brentmsg.costh=cos(th);
  brentmsg.sinth=sin(th);
  brentmsg.D2costh=D2*brentmsg.costh;
  brentmsg.D2sinth=D2*brentmsg.sinth;
  brentmsg.Pjp[0] = 0.0;
  brentmsg.Pjp[1] = brentmsg.D2costh;
  brentmsg.Pjp[2] = brentmsg.D2sinth;

  for (k1=0; k1 < 3; k1++)
    brentmsg.PjCi[k1] = brentmsg.Pjp[k1] - CipGbl[k1]; 
  brentmsg.lambda = scalProd(brentmsg.PjCi,nipGbl);
  for (k1=0; k1 < 3; k1++)
    {
      brentmsg.Uip[k1] = CipGbl[k1] + brentmsg.lambda*nipGbl[k1];
      brentmsg.UipPjp[k1] = brentmsg.Uip[k1] - brentmsg.Pjp[k1];

      //PjPi[k1] = Pjp[k1] - (CipGbl[k1] + lambda*nipGbl[k1]);
      brentmsg.minPgbl[k1] = brentmsg.Pjp[k1]; 
    }
#ifdef MC_RIMDISK_NORMSQ
  brentmsg.normUipPjpSq = Sqr(brentmsg.UipPjp[0])+Sqr(brentmsg.UipPjp[1])+Sqr(brentmsg.UipPjp[2]);
#else
  brentmsg.normUipPjp = calc_norm(brentmsg.UipPjp);
#endif
}
void calcdrimdisk(double th)
{
  double fact;
  calcrimdisk(th);
  fact = -nipGbl[1]*brentmsg.D2sinth+nipGbl[2]*brentmsg.D2costh;
  brentmsg.dUipPjp[0] = nipGbl[0]*fact;
  brentmsg.dUipPjp[1] = nipGbl[1]*fact+brentmsg.D2sinth;
  brentmsg.dUipPjp[2] = nipGbl[2]*fact-brentmsg.D2costh;
}
double ddrimdiskfunc(double th)
{
  /* i è il rim e j il disco */
  //double lambda, Pjp[3], Pip[3], D2, tj[3], PjCi[3], PjPi[3];
  //double UipPjp[3], Uip[3], fact, dUipPjp[3];
  double dfact, ddUipPjp[3];
  int k1, k2;
  //if (!(th==brentmsg.th) || brentmsg.id==0)
    //{
      calcdrimdisk(th);
      brentmsg.id = 3;
   // }

  dfact = -nipGbl[1]*brentmsg.D2costh-nipGbl[2]*brentmsg.D2sinth;
  ddUipPjp[0] = nipGbl[0]*dfact;
  ddUipPjp[1] = nipGbl[1]*dfact+brentmsg.D2costh;
  ddUipPjp[2] = nipGbl[2]*dfact+brentmsg.D2sinth;
#ifdef MC_RIMDISK_NORMSQ
  return signGbl*(2.0*(scalProd(brentmsg.dUipPjp,brentmsg.dUipPjp)+scalProd(ddUipPjp,brentmsg.UipPjp)));
#else
  return signGbl*((scalProd(ddUipPjp,brentmsg.UipPjp)+scalProd(brentmsg.dUipPjp,brentmsg.dUipPjp))/sqrt(brentmsg.normUipPjp)-Sqr(scalProd(brentmsg.dUipPjp,brentmsg.UipPjp))/(brentmsg.normUipPjp,1.5));
#endif
  //tj[0] = -D2*sinth;
  //tj[1] = D2*costh;
  //tj[2] = 0.0;
 //return scalProd(PjPi,tj);
}
double drimdiskfunc(double th)
{
  /* i è il rim e j il disco */
  //double lambda, Pjp[3], Pip[3], D2, tj[3], PjCi[3], PjPi[3];
  //double UipPjp[3], Uip[3], fact, dUipPjp[3];
  double dUipPjp[3], fact, D2;
  int k1, k2;
  D2 = 0.5*Dgbl;

  //if (!(th==brentmsg.th) || brentmsg.id==0)
    //{
      calcdrimdisk(th);
      brentmsg.id = 2;
   // }
#ifdef MC_RIMDISK_NORMSQ
  return signGbl*2.0*scalProd(brentmsg.dUipPjp,brentmsg.UipPjp);
#else
  return signGbl*scalProd(brentmsg.dUipPjp,brentmsg.UipPjp)/sqrt(brentmsg.normUipPjp);
#endif
  //tj[0] = -D2*sinth;
  //tj[1] = D2*costh;
  //tj[2] = 0.0;
  //return scalProd(PjPi,tj);
}

double rimdiskfunc(double th)
{
  /* i è il rim e j il disco */
  //double lambda, Pjp[3], Pip[3], D2, tj[3], PjCi[3], PjPi[3];
  //double UipPjp[3], Uip[3];
  int k1;
  //double sinth, costh;

  //if (!(th == brentmsg.th) || brentmsg.id==0)
    //{
      calcrimdisk(th);
      brentmsg.id = 1;
    //}
  for (k1=0; k1 < 3; k1++)
    minPgbl[k1] = brentmsg.minPgbl[k1];
#ifdef MC_RIMDISK_NORMSQ
  return signGbl*brentmsg.normUipPjpSq;
#else
  return signGbl*brentmsg.normUipPjp;
#endif
  //tj[0] = -D2*sinth;
  //tj[1] = D2*costh;
  //tj[2] = 0.0;
  //return scalProd(PjPi,tj);
}
void versor_to_R_opt(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
#ifdef MC_BENT_DBLCYL
  double Rout[3][3];
  int k1, k2;
#endif
  /* first row vector */
  R[2][0] = ox;
  R[2][1] = oy;
  R[2][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[2][0] && u[1]==R[2][1] && u[2]==R[2][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[2][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[2][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
  /* third row vector */
  vectProdVec(R[1], R[2], u);
 
  for (k=0; k < 3 ; k++)
    R[0][k] = u[k];
  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
}
double calcDistNegHCdiffbrent(int i, int j, double shift[3], int* retchk);
extern double newton1D(double ax, double (*f)(double), double (*df)(double), double (*ddf)(double), double tol, double *xmin);
void solve_quadratic(double coeff[3], int *numsol, double *sol)
{
  double delta, a2inv, sqrtd;
  delta = Sqr(coeff[1]) - 4.0*coeff[2]*coeff[0];
  if (delta > 0.0)
    {
      sqrtd = sqrt(delta);
      a2inv = 1.0/(2.0*coeff[2]);
      sol[0] = (-coeff[1]+sqrtd)*a2inv;
      sol[1] = (-coeff[1]-sqrtd)*a2inv; 
      *numsol = 2;
    } 
  else if (delta == 0)
    {
      sol[0] = -coeff[1]/(2.0*coeff[2]);
      *numsol = 1;
    }
  else
    {
      *numsol = 0;
    }
}
void solve_cubic(double *coeff, int *numsol, double sol[3], int justone)
{
  const double sqrt3=sqrt(3.0);
  /* solution from Abramovitz */
  double s1, s2, q, r, q3, r2, a2, a1, a0, a2sq, H;
  double complex cs1, cs2, s1ps2, s1ms2;
  if (coeff[3]==0)
    {
      printf("orca troia...coeff[3] è zero!\n");
      exit(-1);
    }
  a2 = coeff[2]/coeff[3];
  a1 = coeff[1]/coeff[3];
  a0 = coeff[0]/coeff[3];
  a2sq=Sqr(a2);
  q = (1.0/3.0)*a1-(1.0/9.0)*a2sq;
  r = (1.0/6.0)*(a1*a2 - 3.0*a0)-a2sq*a2/27.0; 
  q3 = Sqr(q)*q;
  r2 = Sqr(r);
  H = r2 + q3;
  //printf("H=%.15G\n", H);
  if (H <= 0.0)
    {
      cs1 = cpow(r + csqrt(H),1.0/3.0);
      cs2 = cpow(r - csqrt(H),1.0/3.0);
      sol[0] = creal(cs1 + cs2) - a2/3.0;
      if (!(justone==1))
	{
	  /* per il calcolo degli zeri di una quartica basta uno zero reale qualsiasi quindi 
	   * con questo flag si puo' ottimizzare calcolandone solo uno */
	  s1ps2 = (cs1+cs2)/2.0; 
	  s1ms2 = (cs1-cs2)/2.0;
	  sol[1] = creal(-s1ps2 - a2/3.0 + I*sqrt3*s1ms2);
	  sol[2] = creal(-s1ps2 - a2/3.0 - I*sqrt3*s1ms2);
	}
      *numsol = 3;
    }
  else
    {
      s1 = cbrt(r + sqrt(H));
      s2 = cbrt(r - sqrt(H));
      sol[0] = s1 + s2 - a2/3.0;
      *numsol = 1;
    }
}
void solve_fourth_deg(double *coeff, int *numsol, double sol[4])
{
  /* solution from H.E. Salzer, "A Note on Solution of Quartic Equations" Am. Math Society Proceedings, 279-281 (1959) */ 
  const double ROUNDOFFERR = 1E-12;
  double x1, A, B, C, D, solc[3], cb[4], m, n;
  double Asq, alpha, beta, gamma, delta, rho, mp, np, m2;
  int nsc;
  if (coeff[4]==0)
    {
      printf("[ERROR: solve_fourth_deg] coeff[4] must be different from 0!\n");
      exit(-1);
    }
  *numsol = 0;
  A = coeff[3]/coeff[4];
  B = coeff[2]/coeff[4];
  C = coeff[1]/coeff[4];
  D = coeff[0]/coeff[4];
  Asq = Sqr(A);
  cb[3] = 1.0;
  cb[2] = -B;
  cb[1] = A*C-4*D;
  cb[0] = D*(4.0*B - Asq)-Sqr(C);
  solve_cubic(cb, &nsc, solc, 1);
  //printf("x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", cb[2], cb[1], cb[0]);
  x1=solc[0];
  m2 = Asq/4.0-B+x1;
  //printf("solc=%.15G numsol=%d m2=%.15G\n", solc[0], nsc, m2);
  
  if (m2 > 0)
    {
      //printf("qui\n");
      m = sqrt(m2);
      n = (A*x1 - 2.0*C)/(4.0*m);
    }
  else if (m2==0.0)
    {
      m = 0.0;
      n = sqrt(Sqr(x1)/4.0 - D);
    }
  else 
    {
      mp = sqrt(-m2);
      np = (A*x1 - 2.0*C)/(4.0*mp);
      alpha=0.5*Asq-x1-B;
      beta = 4.0*np-A*mp;
      rho = sqrt(Sqr(alpha)+Sqr(beta));
      gamma = sqrt((alpha + rho)/2.0); /* notare che alpha + rho è sempre positivo */
      if (gamma == 0.0)
	{
	  delta = sqrt(-alpha); 
	}
      else
	delta  = beta/gamma/2.0;
      if (fabs(mp + delta) < ROUNDOFFERR) /* in questo caso la parte immaginaria è zero e ho due soluzioni reali coincidenti */
	{
	  *numsol = 1;
	  sol[0] = (-A/2.0 + gamma)/2.0;
	}
      if (fabs(mp - delta) < ROUNDOFFERR)
	{
	  sol[*numsol] = (-A/2.0 - gamma)/2.0;
	  *numsol += 1;
	}
      return;
    }
  alpha=0.5*Asq-x1-B;
  beta=4*n-A*m;
  //printf("alpha+beta=%.15G alpha-beta=%.15G gamma=%.15G\n", alpha+beta, alpha-beta,sqrt(alpha+beta));
  //printf("A=%.15G m=%.15G\n", A, m);
  if (alpha+beta > 0)
    {
      /* pair of real solutions */
      gamma = sqrt(alpha+beta);
      sol[0]=-0.25*A+0.5*m+0.5*gamma;
      sol[1]=-0.25*A+0.5*m-0.5*gamma;
      //printf("A=%.15G\n",-0.25*A+0.5*m+0.5*gamma);
      *numsol=2;
    }
  else if (fabs(alpha+beta) < ROUNDOFFERR)
    /* nel caso che alpha+beta è minore di zero ma di poco (entro gli errori di roundoff) 
     * allora significa che abbiamo zeri degeneri e quindi faccio il calcolo con i numeri complessi
     * e se la parte immaginaria è piccola allora considero la soluzione reale */
    {
      sol[0]=-0.25*A+0.5*m;
      sol[1]=-0.25*A+0.5*m;
      *numsol=2;
    }

  if (alpha-beta > 0)
    {    
      /* another pair or real solutions */
      delta = sqrt(alpha-beta);
      sol[*numsol]=-0.25*A-0.5*m+0.5*delta;
      sol[*numsol+1]=-0.25*A-0.5*m-0.5*delta;
      *numsol+=2;
    }
  else if (fabs(alpha - beta) < ROUNDOFFERR)
    {
      sol[*numsol]=-0.25*A+0.5*m;
      sol[*numsol+1]=-0.25*A+0.5*m;
      *numsol+=2;
    }
}

double calcDistNegHCbrent(int i, int j, double shift[3], int* retchk)
{
  static int firstcall=1;
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  static struct brentOpt *mesh;
  int it, kk1, kk2, k2, k1, nz, nl, nn, nng;
  double th, dth, normNSq, ViVj[3], lambdai, lambdaj, Rl[3][3], PjPi[3], PjCi[3], D2, thg, Pjp[3], PiCi[3], lambda, dist, maxmind[2];
  double sp, Q1, Q2, normPiDi, normPjDj, normN, L, D, DiN, DjN, niN[3], njN[3], Djni, Djnj, assex[3], nEy[3], nEz[3];
  double dthg, distleft, distcenter, distright, mindistb, maxdist, semmaxE, semminE, sp1, sp2;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DiCi[3], DiCini, normDiCi, DjCi[3], normDjCi;
  double PiDi[3], PjDj[3], Ai[3], Tj[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3], TipCjnj, TimCjnj;
  double Aj[3], AjDini, AjDinivec[3], AjDi[3], Tip[3], Tim[3], TipCj[3], TimCj[3], Dini;
  double DiCj[3], normDiCj, DiCjnj, Uj[3], DiUj[3], normDiUj, DiUjnj;
  double Tim_perp[3], Tip_perp[3], Tim_para[3], Tip_para[3], normTim_perp, DjCini;
  double Tjm_perp[3], Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp, Tj_para, Tj_perp[3];
  double TiOld[3], TiNew[3], TiNewCj[3], TiNewCjnj, nip[3], Cip[3], Aip[3];	
  double normCiCj, thL, thR, solarr[4][3], coeff[5], solec[4][2], solcc[2][2], solqua[4], solcub[2];	
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3], mindist, PminCip[3], mindistL, mindistR, PminCipL[3], PminCipR[3];
  int kk, j1, j2, numsol;
  
  /* if we have two cylinder with different L or D use calcDistNegHCdiff() function
   * which is able to handle this! */
  if (typesArr[typeOfPart[i]].sax[0]!=typesArr[typeOfPart[j]].sax[0]
      || typesArr[typeOfPart[i]].sax[1] != typesArr[typeOfPart[j]].sax[1])
    return calcDistNegHCdiffbrent(i, j, shift, retchk);

  *retchk = 0; 

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
  meshptsGbl = MESH_PTS;
  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*L*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*L*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*L*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*L*nj[kk];
    }
#ifdef MC_HC_SPHERO_OPT
  if ((sphov=check_spherocyl(CiCj, D, L, Di, Ci, ni, Dj, Cj, nj, &rim)) > 0.0)
    return 1;
#endif
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
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
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
	    if (normPiDi <= 0.5*D && normPjDj <= 0.5*D)
	      {
		Q1 = sqrt(Sqr(D)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  /* =================================== >>> Part A <<< ========================= */
  for (j1=0; j1 < 2; j1++)
    {
      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	    }
	}
      for (j2=0; j2 < 2; j2++)
	{
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);
#if 0
	  if (dostorebump)
	    {
	      printf("normDjUi=%.15G DjUini=%.15G\n", normDjUi, DjUini);
	      printf("Ci=%f %f %f Dj=%f %f %f\n", Ci[0], Ci[1], Ci[2], Dj[0], Dj[1], Dj[2]);
	      printf("DjUi=%.15G %.15G %.15G\n", DjUi[0], DjUi[1], DjUi[2]); 
	      printf("Uj=%.15G %.15G %.15G\n", Ui[0], Ui[1], Ui[2]); 
	      printf("nj=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
	      printf("DjCini= %.15G\n", DjCini);
	    }
#endif 
	  if (normDjUi > D)
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Ui| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

	   */
	  if (normDjUi < D*0.5 && fabs(DjCini) > L*0.5)
	    continue;

	  if (normDjUi < D*0.5 && fabs(DjCini) <= L*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }
#if 1
	  //mindist=find_initial_guess_opt(Ai, Ci, ni, Dj[j2], nj, D, &thg);

	  //printf("Ai-Dj=%f\n", sqrt(Sqr(Ai[0]-Dj[j2][0]) + Sqr(Ai[1]-Dj[j2][1]) +Sqr(Ai[2]-Dj[j2][2])));
	  //printf("Ai=%f %f %f\n", Ai[0], Ai[1], Ai[2]);
#else
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
	  versor_to_R(nj[0], nj[1], nj[2], Rl);
	  for (kk1=0; kk1 < 3; kk1++)
	    {
	      nip[kk1] = 0;
	      //Aip[kk1] = 0;
	      Cip[kk1] = 0;
	      for (kk2=0; kk2 < 3; kk2++)
		{
		  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
		  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[j2][kk2]);
		  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
		} 
	    }
	  //printf("norm Aip=%f\n", calc_norm(Aip));
	  //printf("NormAip=%f\n", sqrt(Sqr(Aip[0])+Sqr(Aip[1]))/(D/2));
	  //printf("thgmin found=%f\n", thg);
#if 0
	  if (Aip[0] >= D/2.)
	    thg = 0;
	  else if (Aip[0] <= -D/2.0)
	    thg = M_PI;
	  else if (Aip[1] < 0.0)
	    thg = 2.0*M_PI-acos(2.0*Aip[0]/D);
	  else
	    thg = acos(2.0*Aip[0]/D);
	  printf("thgcalc=%f\n", thg);
#endif
#if 0
	    {
	      double PP[3];
	      for (kk1=0; kk1 < 3; kk1++)
		{
		  PP[kk1] = 0;
		  for (kk2=0; kk2 < 3; kk2++)
		    {
		      PP[kk1] += Rl[kk2][kk1]*Aip[kk2];
		    } 
		}
	      PP[0] += Dj[j2][0];
	      PP[1] += Dj[j2][1];
	      PP[2] += Dj[j2][2];
	      printf("PP=%f %f %f\n", PP[0], PP[1], PP[2]);
	    }
#endif
	  //printf("Ai=%f %f %f Dj=%f %f %f\n", Ai[0], Ai[1], Ai[2], Dj[j2][0], Dj[j2][1], Dj[j2][2]);
	  //printf("thg=%f Aip=%f %f %f D=%f\n", thg, Aip[0], Aip[1], Aip[2], D);
	  //	  for (kk1=0; kk1 < 3; kk1++)

	  for (kk1=0; kk1 < 3; kk1++)
	    {
	      CipGbl[kk1] = Cip[kk1];
	      nipGbl[kk1] = nip[kk1];
	    }
	  Dgbl = D;
	  /* se l'asse del rim è parallelo al piano del disco bisogna considerare un caso a parte */
	  semminE=D*0.5;
	  if (nip[0]==0.0)
	    {
	      /* devo considerare l'intersezione di due rette (che si ottengono intersecando 
		 il rim con il piano del disco) con il disco */	 
	    }
	  else
	    {
	      /* determina l'ellisse che si ottiene intersecando il rim con il piano del disco */ 
	      lambda =  -Cip[0]/nip[0];
	      /* centro dell'ellisse */
	      rE[0] = 0.0;
	      rE[1] = Cip[1]+lambda*nip[1];
	      rE[2] = Cip[2]+lambda*nip[2];
	      /* versore semiasse minore */
	      if (np[1]==0.0)
		{
		  nEy[0]=0.0;
		  nEy[1]=1.0;
		  nEy[2]=0.0;
		  nEz[0]=0.0;
		  nEz[1]=0.0;
		  nEz[2]=1.0;
		}
	      else
		{
		  nEy[0]=0.0;
		  nEy[2]=1.0/sqrt(1.0+Sqr(nip[2]/nip[1]));
		  nEy[1]=nEy[2]*nip[2]/nip[1];
		  assex[0]=1.0;
		  assex[1]=assex[2]=0.0;
		  vectProdVec(assex,nEy,nEz);
		}
	    }
	  semmaxE=D*0.5/fabs(scalProd(nEy,nip));
	  /* determino le coordinate del centro del cerchio rispetto al riferimento dell'ellisse */
	  rC[0] = 0.0;
	  rC[1] = -rE[1];
	  rC[2] = -rE[2];
	  sp1 = scalProd(rC,nEy);
	  sp2 = scalProd(rC,nEz);
	  rC[1] = sp1;
	  rC[2] = sp2;
	  /* ora trovo l'intersezione dell'ellisse con il cerchio risolvendo l'equazione di quarto grado */
	  /* prima calcolo i coefficienti del polinomio */

	  /* coeff è un array di 5 elementi ossia a,b,c,d,e (coeff. del polinomio c0+c1*x+c2*x^2... )
	   * solarr un array con le numsol soluzioni 
	   * */
	   if (seminE==semaxE)
	    {
	      /* se a=b si ha un equazione quadratica poiché si tratta di due circonferenze */
	      double a,a2,a4,b4,R2,xC,yC,xC2,yC2;
	      a=semminE;
	      a2=Sqr(a);
	      a4=Sqr(a2);
	      R2=Sqr(D2);
	      xC=rC[1];
	      yC=rC[2];
	      xC2=Sqr(xC);
	      yC2=Sqr(yC);
	      if (xC!=0)
		{
		  coeff[2] = 1.0 + yC2/xC2;
		  coeff[1] = -yC - (a2*yC)/xC2 + (R2*yC)/xC2 - yC2*yC/xC2;
		  coeff[0] = -(a2/2.0) - R2/2.0 + a4/(4.0*xC2) - (a2*R2)/(2.0*xC2) + 
		    R2*R2/(4.0*xC2) + xC2/4.0 + yC2/2.0 + (a2*yC2)/(2.0*xC2) - (R2*yC2)/(2.0*xC2) 
		    + yC*yC2/(4.0*xC2);
		  /* sto risolvendo in y */
		  solve_quadratic(coeff, &numsol, solcub);
		  /* assegno solcc e calcolo x */
		  for (kk1=0; kk1 < numsol; kk1++)
		    {
		      solcc[kk1][0] = (a2 - R2 + xC2 - 2.0*solcub[kk1]*yC + yC2)/(2.0*xC);
		      solcc[kk1][1] = solcub[kk1];
		    }
		  /* torno al sistema di coordinate del disco */
		  for (kk1=0; kk1 < numsol; kk1++)
		    ellips2disk(solcc[kk1],solarr[kk1]);
		}
	      else if (yC!=0)
		{
		  coeff[2] = 1.0 + xC2/yC2;
		  coeff[1] = -xC - (a2*xC)/yC2 + (R2*xC)/yC2 - xC2*xC/yC2;
		  coeff[0] = -(a2/2.0) - R2/2.0 + xC2/2.0 + a4/(4.0*yC2) - 
		    (a2*R2)/(2.0*yC2) + R2*R2/(4.0*yC2) + 
		    (a2*xC2)/(2.0*yC2) - (R2*xC2)/(2.0*yC2) + xC2*xC2/(4.0*yC2) + yC2/4.0;
		  /* sto risolvendo in x */
		  solve_quadratic(coeff, &numsol, solcub);
		  /* assegno solcc e calcolo y */
		   for (kk1=0; kk1 < numsol; kk1++)
		    {
		      solcc[kk1][0] = solcub[kk1];
		      solcc[kk1][1] = (a2 - R2 - 2.0*solcub[kk1]*xC + xC2 + yC2)/(2.0*yC) ;
		    }
		   /* torno al sistema di coordinate del disco */
		   for (kk1=0; kk1 < numsol; kk1++)
		     ellips2disk(solcc[kk1],solarr[kk1]);
		}
	      else
		{
		  /* se semminE==D2 allora ho infinite soluzioni */

		}
	    }
	  else
	    {
	      double a,b,b2,a2,a4,b4,R2,xC,yC,xC2,yC2, sqA, sqB, sqC, sqD;
	      a=semminE;
	      b=semmaxE;
	      a2=Sqr(a);
	      b2=Sqr(b);
	      a4=Sqr(a2);
	      b4=Sqr(b2);
	      R2=Sqr(D2);
	      xC=rC[1];
	      yC=rC[2];
	      xC2=Sqr(xC);
	      yC2=Sqr(yC);
	      /* notare che coeff[4]=0 solo se a=b che però è un caso a parte! */
	      if (xC!=0)
		{
		  coeff[4] =1.0/(4.0*xC2) + a4/(4.0*b4*xC2) - a2/(2.0*b2*xC2); 
		  coeff[3] = -(yC/xC2) + (a2*yC)/(b2*xC2);
		  coeff[2] = 1.0/2.0 + a2/(2.0*b2) + a2/(2.0*xC2) - a4/(2.0*b2*xC2) - R2/(2.0*xC2) + 
		    (a2*R2)/(2.0*b^2*xC2) + (3.0*yC2)/(2.0*xC2) - (a2*yC2)/(2.0*b2*xC2);
		  coeff[1] = -yC - (a2*yC)/xC2 + (R2*yC)/xC2 - yC2*yC/xC2;
		  coeff[0] = -(a2/2.0) - R2/2.0 + a4/(4.0*xC2) - (a2*R2)/(2.0*xC2) +
		    R2*R2/(4.0*xC2) + xC2/4.0+ yC2/2.0 + (a2*yC2)/(2.0*xC2) - (R2*yC2)/(2.0*xC2) + 
		    yC2*yC2/(4.0*xC2) ;  
		  /* qui risolvo in y */
		  solve_fourth_deg(coeff, &numsol, solqua);
		  /* ora assegno a solec[][] e calcolo x */
		  for (kk1=0; kk1 < numsol; kk1++)
		    {
		      solec[kk1][0] = (a2*b2 - b2*R2 + b2*xC2 + (-a2 + b2)*Sqr(solqua[0]) - 
				       2.0*b2*solqua[0]*yC + b2*yC2)/(2.0*b2*xC) ;
		      solec[kk1][1] = solqua[kk1];
		    }
		  /* torno al sistema di coordinate del disco */
		  for (kk1=0; kk1 < numsol; kk1++)
		    ellips2disk(solec[kk1],solarr[kk1]);
		}
	      else if (yC!=0)
		{
		  coeff[4] = 1.0/(4.0*yC2) - b2/(2.0*a2*yC2) + b4/(4.0*a4*yC2);
		  coeff[3] = -(xC/yC2) + (b2*xC)/(a2*yC2);
		  coeff[2] = 1.0/2.0 + b2/(2.0*a2) + b2/(2.0*yC2) - b4/(2.0*a2*yC2) - R2/(2.0*yC2) + 
		    (b2*R2)/(2.0*a2*yC2) + (3.0*xC2)/(2.0*yC2) - (b2*xC2)/(2.0*a2*yC2);
		  coeff[1] = -xC - (b2*xC)/yC2 + (R2*xC)/yC2 - xC2*xC/yC2;
		  coeff[0] = -(b2/2.0) - R2/2.0 + xC2/2.0 + b4/(4.0*yC2) - (b2*R2)/(2.0*yC2) + R2*R2/(4.0*yC2) + 
		    (b2*xC2)/(2.0*yC2) - (R2*xC2)/(2.0*yC2) + xC2*xC2/(4.0*yC2) + yC2/4.0; 
		  /* qui risolvo in x */
		  solve_fourth_deg(coeff, &numsol, solqua);
		  /* ora assegno solec[][] e calcolo y */
		  for (kk1=0; kk1 < numsol; kk1++)
		    {
		      solec[kk1][0] = solqua[kk1];
		      solec[kk1][1] = (a2*b2 - a2*R2 + (a2 - b2)*Sqr(solqua[kk1]) -
				       2.0*a2*solqua[kk1]*xC + a2*xC2 + a2*yC2)/(2.0*a2*yC);
		    }
		  /* torno al sistema di coordinate del disco */
		  for (kk1=0; kk1 < numsol; kk1++)
		    ellips2disk(solec[kk1],solarr[kk1]);
		}
	      else
		{
		  if (a <= D2 && b >= D2)
		    {
		      sqA = R2-a2;
		      sqB = 1.0-a2/b2;
		      sqC = a2*(b2-R2);
		      sqD = b2 - a2;
		      solec[0][0] = solec[1][0] = -sqrt(sC/sqD);
	    	      solec[0][1] = -(sqrt(sqA/sqB));
    		      solec[1][1] = -solec[0][1];
		      solec[2][0] = solec[3][0] = -solec[0][0];
		      solec[2][1] = solec[0][1];
		      solec[3][1] = solec[1][1];
		      numsol = 4;
		    }
		  else
		    {
		      numsol = 0;
		    }
		}
	      /* in questo caso ho due circonferenze e quindi al più due soluzioni ossia ho un'equazione quadratica */
	    }
	  /* verifico che le soluzioni siano nella parte di rim che fa parte del cilindro */


	  /* =========================================================================== */
	  //printf("mindist=%f mindst from rimdisk=%.15G\n", mindist, rimdiskfunc(thg));
	  //if (fabs(rimdiskfunc(thg)-mindist) > 1E-7)
	  //exit(-1);
#if 0
	  if (firstcall)
	    {
	      mesh = malloc(sizeof(struct brentOpt)*MESH_PTS);
	      firstcall=0;
	      dth = 2.0*M_PI/((double)MESH_PTS);

	      th=0.0;
	      for (nn=0; nn < MESH_PTS; nn++)
		{
		  calcrimdisk(th);
		  mesh[nn].sinth = brentmsg.sinth;
		  mesh[nn].costh = brentmsg.costh;
		  mesh[nn].id = 3; 
		  for (k1=0; k1 < 3; k1++)
		    {
		      mesh[nn].Pjp[k1] = brentmsg.Pjp[k1];
		    } 
		  th += dth;
		}
	    }
#endif
#if 0
  	  brentmsg.id = 0;
	  maxdist=find_initial_guess_bracket(&thg, MESH_PTS, mesh, &nng, &maxmin);
	  //maxdist=find_initial_guess_bracket(&thg, MESH_PTS, mesh, &nng, &maxmin);
	  find_initial_linecircle(Dgbl, CipGbl, nipGbl, &iniguess);
	  //find_initial_guess_proj(Dgbl, CipGbl, nipGbl, &iniguess);
#endif
       	  //printf("ax=%f bx(mindist)=%f cx=%f\n", rimdiskfunc(thg-2.0*M_PI/MESH_PTS), rimdiskfunc(thg), rimdiskfunc(thg+2.0*M_PI/MESH_PTS));
#if 0
	  signGbl=-1.0;
	  maxdist = -newton1D(thg, rimdiskfunc, drimdiskfunc, ddrimdiskfunc, 1E-7, &th);
#endif

#if 0
	  if (maxdist < Sqr(maxmin.max))
	    {
	      printf("maxdist find=%.15G mindist=%.15G\n", maxmin.max,maxmin.min);	
	      printf("maxdist=%.15G\n", maxdist);
	      exit(-1);
	    }	
	  for (nn=nng; ;nn++)
	    {
	      nnL = nn-1;
	      if (nnL < MESH_PTS)
	      if mesh[nn]   
	    }
#endif	
#if 0 
	  /* at most one has two minima and two maxima... */
	  signGbl= 1.0;
	  /* left */
	  //steepest1D(th-2.0*M_PI/MESH_PTS, rimdiskfunc, drimdiskfunc, ddrimdiskfunc, 1E-3, &thg);
	  thg = iniguess.thg[iniguess.which];
	  for (k1=0; k1 < iniguess.num; k1++)
	    printf("iniguess[%d]:%.15G\n", k1, iniguess.thg[k1]);
	  printf("best guess[%d]=%.15G\n", iniguess.which, iniguess.thg[iniguess.which]);
	  printf("multmin=%d which=%d\n", iniguess.multmin, iniguess.which);
  	  brentmsg.id = 0;
	  if (iniguess.multmin == 0)
	    {
	      dthg = fabs(iniguess.thg[0]-iniguess.thg[1]);
	      mindistL = dbrent(thg-dthg, thg, thg+dthg, rimdiskfunc, drimdiskfunc, 1.0E-14, &th);
	    }
	  else if (iniguess.multmin==1 && iniguess.which==2)
	    {
	      dthg = fabs(iniguess.thg[0]-iniguess.thg[2]);
	      mindistL = dbrent(thg-dthg, thg, thg+dthg, rimdiskfunc, drimdiskfunc, 1.0E-14, &th);
	    }
	  else
	    {
	      dthg = fabs(iniguess.thg[0]-iniguess.thg[1])*0.5;
	      thg = iniguess.thg[0];
	      brentmsg.id = 0;
	      mindistL = dbrent(thg-dthg, thg, thg+dthg, rimdiskfunc, drimdiskfunc, 1.0E-14, &thL);
	      thg = iniguess.thg[1];
	      brentmsg.id = 0;
	      mindistR = dbrent(thg-dthg, thg, thg+dthg, rimdiskfunc, drimdiskfunc, 1.0E-14, &thR);
	      if (fabs(mindistR - mindistL) > 1E-7)
		{
		  printf("mindistL=%.15G th=%.15G mindistR=%.15G th=%.15G\n", mindistL, thL, mindistR, thR);
		  mindistL = -1;
		}
	      //mindistL = newton1D(thg, rimdiskfunc, drimdiskfunc, ddrimdiskfunc, 1E-14, &th);
	    }
#if 1
	  if (mindistL < 0)
	    {
	      FILE *f;
	      printf(">>>>>>>>>>>>>>>><<<<<<<<<<<<<<\n");
	      f=fopen("dtheta.dat", "w+");
	      th = 0;
	      for (k1 = 0; k1 < MESH_PTS; k1++)
		{
		  for (k2=0; k2 < 3; k2++)
		    brentmsg.PjCi[k2] = mesh[k1].Pjp[k2] - CipGbl[k2]; 
		  brentmsg.lambda = scalProd(brentmsg.PjCi,nipGbl);
		  for (k2=0; k2 < 3; k2++)
		    {
		      brentmsg.Uip[k2] = CipGbl[k2] + brentmsg.lambda*nipGbl[k2];
		      brentmsg.UipPjp[k2] = brentmsg.Uip[k2] - mesh[k1].Pjp[k2];
		      //PjPi[k1] = Pjp[k1] - (CipGbl[k1] + lambda*nipGbl[k1]);
		    }
		  dist = calc_norm(brentmsg.UipPjp);
		  fprintf(f, "%.15G %.15G\n", th, Sqr(dist));
		  //printf("bracket 2 (%f,%d)=%.15G\n", th, k1, dist);
		  if (k1==0 || dist < mindist)
		    {
		      mindist = dist;
		      thg = th;
		    }
		  th+=dth;
		}
	      fclose(f); 
	      exit(-1);
	    }
#endif 
	  //printf("maxdist=%f th=%f mindistL=%f\n", maxdist, thg, mindistL);
	  for (k1=0; k1 < 3; k1++)
	    {
	      PminCip[k1] = minPgbl[k1] - Cip[k1];
	    }
#if 0
	  /* right */
	  //steepest1D(th+2.0*M_PI/MESH_PTS, rimdiskfunc, drimdiskfunc, ddrimdiskfunc, 1E-3, &thg);
	  if (iniguess.num == 2)
	    {
	      brentmsg.id = 0;
	      thg = iniguess.thg[1];
	      mindistR = newton1D(thg, rimdiskfunc, drimdiskfunc, ddrimdiskfunc, 1E-10, &th);
#if 1
	      if (mindistR < 0)
		{
		  FILE *f;
		  printf("RRRRRRRR>>>>>>>>>>>>>>>><<<<<<<<<<<<<<\n");
		  f=fopen("dtheta.dat", "w+");
		  th = 0;
		  for (k1 = 0; k1 < MESH_PTS; k1++)
		    {
		      for (k2=0; k2 < 3; k2++)
			brentmsg.PjCi[k2] = mesh[k1].Pjp[k2] - CipGbl[k2]; 
		      brentmsg.lambda = scalProd(brentmsg.PjCi,nipGbl);
		      for (k2=0; k2 < 3; k2++)
			{
			  brentmsg.Uip[k2] = CipGbl[k2] + brentmsg.lambda*nipGbl[k2];
			  brentmsg.UipPjp[k2] = brentmsg.Uip[k2] - mesh[k1].Pjp[k2];
			  //PjPi[k1] = Pjp[k1] - (CipGbl[k1] + lambda*nipGbl[k1]);
			}
		      dist = calc_norm(brentmsg.UipPjp);
		      fprintf(f, "%.15G %.15G\n", th, Sqr(dist));
		      //printf("bracket 2 (%f,%d)=%.15G\n", th, k1, dist);
		      if (k1==0 || dist < mindist)
			{
			  mindist = dist;
			  thg = th;
			}
		      th+=dth;
		    }
		  fclose(f); 
		  exit(-1);
		}
#endif
	      for (k1=0; k1 < 3; k1++)
		{
		  PminCipR[k1] = minPgbl[k1] - Cip[k1];
		}
	      if (mindistL < mindistR)	
		{
		  for (k1=0; k1 < 3; k1++)
		    {
		      PminCip[k1] = PminCipL[k1];
		    }
		}
	      else
		{
		  for (k1=0; k1 < 3; k1++)
		    {
		      PminCip[k1] = PminCipR[k1];
		    }
		}
	    }
	  else
	    {
    	      for (k1=0; k1 < 3; k1++)
		{
		  PminCip[k1] = PminCipL[k1];
		}
	    }
#endif
	  //printf("mindist=%f dist=%.15G step=%d\n", mindist, dist, Oparams.curStep);
#if 0
	  brentmsg.id = 0;
	  dist=brent(thg-2.0*M_PI/MESH_PTS, thg, thg+2.0*M_PI/MESH_PTS, rimdiskfunc, 1.0E-14, &th);
	  //printf("mindist=%f dist=%.15G step=%d\n", mindist, dist, Oparams.curStep);
#endif
	  Tj_para = scalProd(PminCip,nip);
	  for (k1=0; k1 < 3; k1++)
	    Tj_perp[k1] = PminCip[k1] - Tj_para*nip[k1];
#if 0
	    {
	      int kkk;
	      const int mp=16;
	      double thgold, mindist2, dist2, minPgbl2[3], thold;
	      //printf("rimdisk4(%f)=%.15G\n", th, rimdiskfunc(th));
	      brentmsg.id = 0;
	      thgold = thg;
	      mindist2=find_initial_guess_bracket2(&thg, mp);
	      for (kkk=0; kkk < 3; kkk++)
		minPgbl2[kkk] = minPgbl[kkk];
	      brentmsg.id = 0;
	      thold = th;
	      dist2=dbrent(thg-2.0*M_PI/mp, thg, thg+2.0*M_PI/mp, rimdiskfunc, drimdiskfunc, 1.0E-14, &th);
	      printf("th=%.15G thold=%.15G\n", th, thold);
	      printf("x1 x2 x3=%f %f %f 2 x1 x2 x3=%f %f %f\n", thg-2.0*M_PI/mp, thg, thg+2.0*M_PI/mp, thgold-2.0*M_PI/MESH_PTS, thgold, thgold+2.0*M_PI/MESH_PTS);
	      //brentmsg.id = 0;
	      //dist2=brent(thg-2.0*M_PI/MESH_PTS, thg, thg+2.0*M_PI/MESH_PTS, rimdiskfunc, 1.0E-14, &th);
	      //printf("rimdisk8(%f)=%.15G\n", th, rimdiskfunc(th));
	      if (fabs(minPgbl[0]-minPgbl2[0])> 1E-7 ||
		 fabs(minPgbl[1]-minPgbl2[1])> 1E-7 ||
		 fabs(minPgbl[2]-minPgbl2[2])> 1E-7)
		{
		  brentmsg.id=0;
		  find_initial_guess_bracket2(NULL, 1000000);
		  printf("thg=%.15G thgold=%.15G\n", thg, thgold);
		  printf("rimdisk(%f)=%.15G\n", th, rimdiskfunc(th));
		  printf(">>>> dist=%.15G\n",find_initial_guess_bracket2(&thg, 1000000));
		  printf("nip=%.15G %.15G %.15G\n", nip[0], nip[1], nip[2]);
		  printf("dist tramite rimdisk(%f)=%.15G\n", thg, rimdiskfunc(thg));
		  printf("Ci=%f %f %f Dj=%f %f %f\n", Cip[0], Cip[1], Cip[2], Dj[j2][0], Dj[j2][1], Dj[j2][2]);
		  printf("calc_norm Ci=%.15G\n", calc_norm(Cip));
		  printf("P8=%.15G %.15G %.15G P4=%.15G %.15G %.15G\n", minPgbl2[0],minPgbl2[1],minPgbl2[2],minPgbl[0],minPgbl[1],minPgbl[2]);
		  printf("normPgbl=%.15G normPgbl2=%.15G\n", calc_norm(minPgbl), calc_norm(minPgbl2));
		  printf("dist=%.15G dist2=%.15G\n", dist, dist2);
		  printf("mindist=%.15G mindist2=%.15G\n", mindist, mindist2);
		  exit(-1);
		}
	    }
#endif
#if 0
	  if (rimdiskfunc(thg) > rimdiskfunc(thg+2.0*M_PI/MESH_PTS)
	       || rimdiskfunc(thg) > rimdiskfunc(thg-2.0*M_PI/MESH_PTS))
	    {
	      printf("boh...\n");
	      exit(-1);
	    }
	  
#endif
	  if ( (fabs(Tj_para) <= L*0.5 && calc_norm(Tj_perp) <= D*0.5))
	    {
	      return -1;
	    }
#endif
	  //printf("dist=%f th=%.15G (min=%f max=%f)\n", dist, th, thg-M_PI/MESH_PTS, thg+M_PI/MESH_PTS);
	  //if (dist < D2)
	  //return -1;
#if 0
	  D2 = 0.5*D;
	  Pjp[0] = D2*cos(th);
	  Pjp[1] = D2*sin(th);
	  Pjp[2] = 0.0;
	  for (k1=0; k1 < 3; k1++)
	    PjCi[k1] = Pjp[k1] - Cip[k1]; 
	  lambda = scalProd(PjCi,nipGbl);
	  for (k1=0; k1 < 3; k1++)
	    {
	      PjPi[k1] = Pjp[k1] - (Cip[k1] + lambda*nip[k1]);
	    }
	  if (calc_norm(PjPi) < D2)
	    return -1;
#endif
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	    }
	}

    }
  /* =================================== >>> Part B <<< ========================= */
#if 0
  for (j1=0; j1 < 2; j1++)
    {
      for (kk=0; kk < 3; kk++)
	DiCj[kk] = Di[j1][kk] - Cj[kk];
      normDiCj = calc_norm(DiCj);
      DiCjnj = scalProd(DiCj,nj);
      for (kk=0; kk < 3; kk++)
	{
	  Uj[kk] = Cj[kk] + DiCjnj*nj[kk];
	  DiUj[kk] = Di[j1][kk] - Uj[kk];
	}

      DiUjnj = scalProd(DiUj,nj);
      normDiUj = calc_norm(DiUj);
#ifdef DEBUG_HCMC
      if (dostorebump)
	{
	  printf("B normDiUj=%.15G DiUjnj=%.15G\n", normDiUj, DiUjnj);
	  printf("B Cj=%f %f %f Di=%f %f %f\n", Cj[0], Cj[1], Cj[2], Di[j1][0], Di[j1][1], Di[j1][2]);
	  printf("B DiUj=%.15G %.15G %.15G\n", DiUj[0], DiUj[1], DiUj[2]); 
	  printf("B Uj=%.15G %.15G %.15G\n", Uj[0], Uj[1], Uj[2]); 
	  printf("B nj=%.15G %.15G %.15G\n", nj[0], nj[1], nj[2]);
	  printf("DjCini= %.15G\n", DjCini);
	}
#endif 
 
      if (normDiUj > D)
	continue;

      if (normDiUj < D*0.5 && fabs(DiCjnj) > L*0.5)
	continue;

      if (normDiUj < D*0.5 && fabs(DiCjnj) <= L*0.5)
	{
#ifdef DEBUG_HCMC
	  if (dostorebump)
	    printf("B #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	  return -1;
	}      
      for (kk=0; kk < 3; kk++)
	{
	  Aj[kk] = Cj[kk];
	}
      for (it = 0; it < MAX_ITERATIONS; it++)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      AjDi[kk] = Aj[kk] - Di[j1][kk];
	    }
	  AjDini = scalProd(AjDi,ni);
	  vectProdVec(AjDi,ni,AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    VV[kk] =  0.5*D*(AjDi[kk]-AjDini*ni[kk])/calc_norm(AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tip[kk] = Di[j1][kk] + VV[kk];
	      Tim[kk] = Di[j1][kk] - VV[kk];
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
	  normTim_perp = calc_norm(Tip_perp);
	  for (kk=0; kk < 3; kk++)
	    TiOld[kk] = TiNew[kk];
	  if (calc_norm(Tim_perp) < calc_norm(Tip_perp))
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
#ifdef DEBUG_HCMC
	  printf("B it=%d Ajold=%.15G %.15G %.15G\n", it, Aj[0], Aj[1], Aj[2]);
#endif
	  for (kk=0; kk < 3; kk++)
	    Aj[kk] = TiNewCjnj*nj[kk] + Cj[kk]; 
#ifdef DEBUG_HCMC
	  printf("B it=%d Ajnew=%.15G %.15G %.15G TiNewCjnj=%.15G\n", it, Aj[0], Aj[1], Aj[2], TiNewCjnj);
	  printf("B Ci=%.15G %.15G %.15G\n", Cj[0], Cj[1], Cj[2]);
	  printf("B ni=%.15G %.15G %.15G\n", nj[0], nj[1], nj[2]);
	  printf("B #1 number of iterations=%d Tiold=%.15G %.15G %.15G Tinew=%.15G %.15G %.15G\n",it, 
	     TiOld[0], TiOld[1], TiOld[2], TiNew[0], TiNew[1], TiNew[2]);

#endif
	
	  if ( it > 0 && check_convergence(TiOld,TiNew) ) 
	    {
	      break;
	    }
	} 
      totitsHC += it;
#ifdef DEBUG_HCMC
      printf("B #1 number of iterations=%d Tiold=%.15G %.15G %.15G Tinew=%.15G %.15G %.15G\n",it, 
	     TiOld[0], TiOld[1], TiOld[2], TiNew[0], TiNew[1], TiNew[2]);
#endif
 
      if (it >= MAX_ITERATIONS)
       	{
 	  printf("MAX ITERATIONS REACHED IN B\n");
	  *retchk=1;
#ifdef DEBUG_HCMC
	  //exit(-1);
#endif
 	  return -1;
  	}
      
     // printf("#2 number of iterations=%d\n",it);
      if ( (calc_norm(Tip_para) <= L*0.5 && calc_norm(Tip_perp) <= D*0.5)||
	   (calc_norm(Tim_para) <= L*0.5 && calc_norm(Tim_perp) <= D*0.5) )
	{
#ifdef DEBUG_HCMC
	  if (dostorebump)
	    printf("B #2 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	  return -1;
	}
    }
#endif
  numcallsHC += 4.0; 

  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < D && fabs(lambdai) < 0.5*L && fabs(lambdaj) < 0.5*L)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 1;
}
double calcDistNegHCdiffbrent(int i, int j, double shift[3], int* retchk)
{
  /* NOTA 291117: va ancora testata! */
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int it, k2, k1, kk1, kk2, nl, nn, nz;
  static struct brentOpt *mesh;
  double normNSq, ViVj[3], lambdai, lambdaj, Li, Diami, Lj, Diamj, dist, mindist, Tj_para, Tj_perp[3]; 
  double LiTmp, LjTmp, DiamiTmp, DiamjTmp, nip[3], Cip[3], th, thg, PminCip[3], Rl[3][3];
  double sp, Q1, Q2, normPiDi, normPjDj, normN, DiN, DjN, niN[3], njN[3], Djni, Djnj;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DiCi[3], DiCini, normDiCi, DjCi[3], normDjCi;
  double PiDi[3], PjDj[3], Ai[3], Tj[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3], TipCjnj, TimCjnj;
  double Aj[3], AjDini, AjDinivec[3], AjDi[3], Tip[3], Tim[3], TipCj[3], TimCj[3], Dini;
  double DiCj[3], normDiCj, DiCjnj, Uj[3], DiUj[3], normDiUj, DiUjnj;
  double Tim_perp[3], Tip_perp[3], Tim_para[3], Tip_para[3], normTim_perp, DjCini;
  double Tjm_perp[3], Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp;
  double TiOld[3], TiNew[3], TiNewCj[3], TiNewCjnj, Tjpara, Tjperp[3];	
  double normCiCj;	
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3];
  int kk, j1, j2;
  *retchk = 0; 

  // return calcDistNegHCsame(i, j, shift, retchk);
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
  Li = 2.0*typesArr[typeOfPart[i]].sax[0];
  Diami = 2.0*typesArr[typeOfPart[i]].sax[1];
  Lj = 2.0*typesArr[typeOfPart[j]].sax[0];
  Diamj = 2.0*typesArr[typeOfPart[j]].sax[1];

  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }

  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*Li*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*Li*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*Lj*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*Lj*nj[kk];
    }
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= 0.5*(Li+Lj))
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
	    if (sp == 0 && calc_norm(VV) < 0.5*(Diami+Diamj))
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
#ifdef DEBUG_HCMC
	    printf("Di=%f %f %f\n", Di[j1][0], Di[j1][1], Di[j1][2]);
	    printf("Dj=%f %f %f\n", Dj[j2][0], Dj[j2][1], Dj[j2][2]);
	    printf("normPiDi: %f normPjDj=%f\n", normPiDi, normPjDj);
	    printf("0.5*Diami=%f 0.5*Diamj=%f\n", 0.5*Diami, 0.5*Diamj);
#endif
	    if (normPiDi <= 0.5*Diami && normPjDj <= 0.5*Diamj)
	      {
		Q1 = sqrt(Sqr(Diami)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(Diamj)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  /* =================================== >>> Part A <<< ========================= */
  for (j1=0; j1 < 2; j1++)
    {

      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      DiamiTmp = Diami;
	      DiamjTmp = Diamj;
	      LiTmp = Li;
	      LjTmp = Lj;
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	      Diami = Diamj;
	      Diamj = DiamiTmp;
	      Li = Lj;
	      Lj = LiTmp;
	    }
	}
      for (j2=0; j2 < 2; j2++)
	{
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);

	  if (normDjUi > 0.5*(Diami+Diamj))
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Uj| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

	   */
	  if (normDjUi < Diami*0.5 && fabs(DjCini) > Li*0.5)
	    continue;

	  if (normDjUi < Diami*0.5 && fabs(DjCini) <= Li*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }
#if 1
	  //find_initial_guess(Ai, Ci, ni, Dj[j2], nj, Diamj);
#else
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
	  //mindist=find_initial_guess_opt(Ai, Ci, ni, Dj[j2], nj, Diamj, &thg);
	  versor_to_R(nj[0], nj[1], nj[2], Rl);
	  for (kk1=0; kk1 < 3; kk1++)
	    {
	      nip[kk1] = 0;
	      //Aip[kk1] = 0;
	      Cip[kk1] = 0;
	      for (kk2=0; kk2 < 3; kk2++)
		{
		  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
		  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[j2][kk2]);
		  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
		} 
	    }

	  for (kk1=0; kk1 < 3; kk1++)
	    {
	      CipGbl[kk1] = Cip[kk1];
	      nipGbl[kk1] = nip[kk1];
	    }
	  Dgbl = Diamj;
	  brentmsg.id = 0;
#if 0
	  /* bracketing */
	  dth = 2.0*M_PI/MESH_PTS;
	  th = 0;
	  mindist = -1;
	  for (k1 = 0; k1 < MESH_PTS; k1++)
	    {
	      dist = rimdiskfunc(th);
	      if (k1==0 || dist < mindist)
		{
		  mindist = dist;
		  thg = th;
		}
	      th+=dth;
	    }
#endif
       	  //printf("ax=%f bx(mindist)=%f cx=%f\n", rimdiskfunc(thg-2.0*M_PI/MESH_PTS), rimdiskfunc(thg), rimdiskfunc(thg+2.0*M_PI/MESH_PTS));
	  //mindist=find_initial_guess_bracket(&thg, MESH_PTS, mesh);

	  dist=dbrent(thg-2.0*M_PI/MESH_PTS, thg, thg+2.0*M_PI/MESH_PTS, rimdiskfunc, drimdiskfunc, 1.0E-14, &th);
	  //dist=brent(thg-2.0*M_PI/MESH_PTS, thg, thg+2.0*M_PI/MESH_PTS, rimdiskfunc, 1.0E-7, &th);
	  for (k1=0; k1 < 3; k1++)
	    {
	      PminCip[k1] = minPgbl[k1] - Cip[k1];
	    }
	  Tj_para = scalProd(PminCip,nip);
	  for (k1=0; k1 < 3; k1++)
	    Tj_perp[k1] = PminCip[k1] - Tj_para*nip[k1];
	  if ( (fabs(Tj_para) <= Li*0.5 && calc_norm(Tj_perp) <= Diami*0.5))
	    {
	      return -1;
	    }
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	      Diami = DiamiTmp;
	      Diamj = DiamjTmp;
	      Li = LiTmp;
	      Lj = LjTmp;
	    }
	}

    }
  /* =================================== >>> Part B <<< ========================= */
  numcallsHC += 4.0; 

  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < 0.5*(Diami+Diamj) && fabs(lambdai) < 0.5*Li && fabs(lambdaj) < 0.5*Lj)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 1;
}
#endif
double calcDistNegHC(int i, int j, double shift[3], int* retchk)
{
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int it, k2;
  double normNSq, ViVj[3], lambdai, lambdaj;
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
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3];
  int kk, j1, j2;

#ifdef HC_ALGO_OPT
  return calcDistNegHCbrent(i, j, shift, retchk);
#endif
  /* if we have two cylinder with different L or D use calcDistNegHCdiff() function
   * which is able to handle this! */
  if (typesArr[typeOfPart[i]].sax[0]!=typesArr[typeOfPart[j]].sax[0]
      || typesArr[typeOfPart[i]].sax[1] != typesArr[typeOfPart[j]].sax[1])
    return calcDistNegHCdiff(i, j, shift, retchk);

  *retchk = 0; 

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
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*L*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*L*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*L*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*L*nj[kk];
    }
#ifdef MC_HC_SPHERO_OPT
  if ((sphov=check_spherocyl(CiCj, D, L, Di, Ci, ni, Dj, Cj, nj, &rim)) > 0.0)
    return 1;
#endif
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
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
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
	    if (normPiDi <= 0.5*D && normPjDj <= 0.5*D)
	      {
		Q1 = sqrt(Sqr(D)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  /* =================================== >>> Part A <<< ========================= */
  for (j1=0; j1 < 2; j1++)
    {
      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	    }
	}
      for (j2=0; j2 < 2; j2++)
	{
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);
#if 0
	  if (dostorebump)
	    {
	      printf("normDjUi=%.15G DjUini=%.15G\n", normDjUi, DjUini);
	      printf("Ci=%f %f %f Dj=%f %f %f\n", Ci[0], Ci[1], Ci[2], Dj[0], Dj[1], Dj[2]);
	      printf("DjUi=%.15G %.15G %.15G\n", DjUi[0], DjUi[1], DjUi[2]); 
	      printf("Uj=%.15G %.15G %.15G\n", Ui[0], Ui[1], Ui[2]); 
	      printf("nj=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
	      printf("DjCini= %.15G\n", DjCini);
	    }
#endif 
	  if (normDjUi > D)
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Ui| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

	   */
	  if (normDjUi < D*0.5 && fabs(DjCini) > L*0.5)
	    continue;

	  if (normDjUi < D*0.5 && fabs(DjCini) <= L*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }
#if 1
	  find_initial_guess(Ai, Ci, ni, Dj[j2], nj, D);

#else
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
	  for (it = 0; it < MAX_ITERATIONS; it++)
	    {
	      for (kk=0; kk < 3; kk++)
		{
		  AiDj[kk] = Ai[kk] - Dj[j2][kk];
		}
	      AiDjnj = scalProd(AiDj,nj);
	      vectProdVec(AiDj,nj,AiDjnjvec);
	      for (kk=0; kk < 3; kk++)
		VV[kk] =  0.5*D*(AiDj[kk]-AiDjnj*nj[kk])/calc_norm(AiDjnjvec);
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
	      normTjm_perp = calc_norm(Tjp_perp);
	      for (kk=0; kk < 3; kk++)
		TjOld[kk] = TjNew[kk];
	      if (calc_norm(Tjm_perp) < calc_norm(Tjp_perp))
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

#ifdef DEBUG_HCMC
	      printf("j1=%d A it=%d Aiold=%.15G %.15G %.15G\n", j1, it, Ai[0], Ai[1], Ai[2]);
#endif
	      for (kk=0; kk < 3; kk++)
		Ai[kk] = TjNewCini*ni[kk] + Ci[kk]; 
#ifdef DEBUG_HCMC
	      printf("A it=%d Ainew=%.15G %.15G %.15G TjNewCini=%.15G\n", it, Ai[0], Ai[1], Ai[2], TjNewCini);
	      printf("A Ci=%.15G %.15G %.15G\n", Ci[0], Ci[1], Ci[2]);
	      printf("A ni=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
#endif
	      if ( it > 0 && check_convergence(TjOld,TjNew) ) 
		break;
	    }
	  totitsHC += it;
#ifdef DEBUG_HCMC
	  printf("A #1 number of iterations=%d Tjold=%.15G %.15G %.15G Tjnew=%.15G %.15G %.15G\n",it, 
		 TjOld[0], TjOld[1], TjOld[2], TjNew[0], TjNew[1], TjNew[2]);
#endif
	  if (it >= MAX_ITERATIONS)
	    {
	      printf("MAX ITERATIONS REACHED in A!\n");
	      *retchk=1;
	      return -1;
	    }
	  if ( (calc_norm(Tjp_para) <= L*0.5 && calc_norm(Tjp_perp) <= D*0.5)||
	       (calc_norm(Tjm_para) <= L*0.5 && calc_norm(Tjm_perp) <= D*0.5) )
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #2 disk-rim\n");
#endif	   
	      return -1;
	    }
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	    }
	}

    }
  /* =================================== >>> Part B <<< ========================= */
#if 0
  for (j1=0; j1 < 2; j1++)
    {
      for (kk=0; kk < 3; kk++)
	DiCj[kk] = Di[j1][kk] - Cj[kk];
      normDiCj = calc_norm(DiCj);
      DiCjnj = scalProd(DiCj,nj);
      for (kk=0; kk < 3; kk++)
	{
	  Uj[kk] = Cj[kk] + DiCjnj*nj[kk];
	  DiUj[kk] = Di[j1][kk] - Uj[kk];
	}

      DiUjnj = scalProd(DiUj,nj);
      normDiUj = calc_norm(DiUj);
#ifdef DEBUG_HCMC
      if (dostorebump)
	{
	  printf("B normDiUj=%.15G DiUjnj=%.15G\n", normDiUj, DiUjnj);
	  printf("B Cj=%f %f %f Di=%f %f %f\n", Cj[0], Cj[1], Cj[2], Di[j1][0], Di[j1][1], Di[j1][2]);
	  printf("B DiUj=%.15G %.15G %.15G\n", DiUj[0], DiUj[1], DiUj[2]); 
	  printf("B Uj=%.15G %.15G %.15G\n", Uj[0], Uj[1], Uj[2]); 
	  printf("B nj=%.15G %.15G %.15G\n", nj[0], nj[1], nj[2]);
	  printf("DjCini= %.15G\n", DjCini);
	}
#endif 
 
      if (normDiUj > D)
	continue;

      if (normDiUj < D*0.5 && fabs(DiCjnj) > L*0.5)
	continue;

      if (normDiUj < D*0.5 && fabs(DiCjnj) <= L*0.5)
	{
#ifdef DEBUG_HCMC
	  if (dostorebump)
	    printf("B #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	  return -1;
	}      
      for (kk=0; kk < 3; kk++)
	{
	  Aj[kk] = Cj[kk];
	}
      for (it = 0; it < MAX_ITERATIONS; it++)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      AjDi[kk] = Aj[kk] - Di[j1][kk];
	    }
	  AjDini = scalProd(AjDi,ni);
	  vectProdVec(AjDi,ni,AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    VV[kk] =  0.5*D*(AjDi[kk]-AjDini*ni[kk])/calc_norm(AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tip[kk] = Di[j1][kk] + VV[kk];
	      Tim[kk] = Di[j1][kk] - VV[kk];
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
	  normTim_perp = calc_norm(Tip_perp);
	  for (kk=0; kk < 3; kk++)
	    TiOld[kk] = TiNew[kk];
	  if (calc_norm(Tim_perp) < calc_norm(Tip_perp))
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
#ifdef DEBUG_HCMC
	  printf("B it=%d Ajold=%.15G %.15G %.15G\n", it, Aj[0], Aj[1], Aj[2]);
#endif
	  for (kk=0; kk < 3; kk++)
	    Aj[kk] = TiNewCjnj*nj[kk] + Cj[kk]; 
#ifdef DEBUG_HCMC
	  printf("B it=%d Ajnew=%.15G %.15G %.15G TiNewCjnj=%.15G\n", it, Aj[0], Aj[1], Aj[2], TiNewCjnj);
	  printf("B Ci=%.15G %.15G %.15G\n", Cj[0], Cj[1], Cj[2]);
	  printf("B ni=%.15G %.15G %.15G\n", nj[0], nj[1], nj[2]);
	  printf("B #1 number of iterations=%d Tiold=%.15G %.15G %.15G Tinew=%.15G %.15G %.15G\n",it, 
	     TiOld[0], TiOld[1], TiOld[2], TiNew[0], TiNew[1], TiNew[2]);

#endif
	
	  if ( it > 0 && check_convergence(TiOld,TiNew) ) 
	    {
	      break;
	    }
	} 
      totitsHC += it;
#ifdef DEBUG_HCMC
      printf("B #1 number of iterations=%d Tiold=%.15G %.15G %.15G Tinew=%.15G %.15G %.15G\n",it, 
	     TiOld[0], TiOld[1], TiOld[2], TiNew[0], TiNew[1], TiNew[2]);
#endif
 
      if (it >= MAX_ITERATIONS)
       	{
 	  printf("MAX ITERATIONS REACHED IN B\n");
	  *retchk=1;
#ifdef DEBUG_HCMC
	  //exit(-1);
#endif
 	  return -1;
  	}
      
     // printf("#2 number of iterations=%d\n",it);
      if ( (calc_norm(Tip_para) <= L*0.5 && calc_norm(Tip_perp) <= D*0.5)||
	   (calc_norm(Tim_para) <= L*0.5 && calc_norm(Tim_perp) <= D*0.5) )
	{
#ifdef DEBUG_HCMC
	  if (dostorebump)
	    printf("B #2 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	  return -1;
	}
    }
#endif
  numcallsHC += 4.0; 

  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < D && fabs(lambdai) < 0.5*L && fabs(lambdaj) < 0.5*L)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 1;
}
#ifdef MC_HC_SPHERO_OPT
/*
 Revision of
 Carlos Vega & Santiago Lago
 Computers Chem. 18, 55-59, 1994

 Subrutine to evaluate the shortest distance between two rods of
 different length

 The original code did not give the symmetry property of the distance for almost parallel rods.
 The coordinates of the centers of the rods should be given in a periodic system

 r1,r2: centers of rods
 w1,w2: unit orientation vectors of rods
 lh1,lh2: halves of the length of rods
 Lv.x,Lv.y,Lv.z the edges of the periodic simulation cell
*/

//----------------- VECTOR operations: -----------------------------------------------------


#define VECT_COMMA ,
#define VECT_PAR (
#define VECT_PSEQ(_,SEP) (_ x)) SEP (_ y)) SEP (_ z))

#define VECT_COMP(x) .x
#define VECT_OP(A,COMP,OP,x) A COMP(x) OP
#define VECT_A_OP_B(A,OP,B,x) VECT_OP(A,VECT_COMP,OP,x) VECT_OP(B,VECT_COMP,,x)

#define VECT_OSEQ_(A,OP,B,SEP,_) \
 VECT_PSEQ(VECT_A_OP_B VECT_PAR A VECT_COMMA OP VECT_COMMA B VECT_COMMA,SEP##_)

#define VECT_OSEQ(A,OP,B,SEP) VECT_OSEQ_(A,OP,B,SEP,)
#define VECT_PROD(A,B) VECT_OSEQ(A,*,B,+)  /* product of A and B */
#define VECT_NORM2(A) VECT_PROD(A,A)  /* square of the norm of A */

#define VECT_OLIST(A,OP,B) VECT_OSEQ_(A,OP,B,VECT_COMMA,) /* (A.x OP B.x), ... */

#define VECT_SEQ(V,SEP) V(x) SEP V(y) SEP V(z)  /* because of the single macro expansion */
#define VECT_LIST(V) VECT_SEQ(V,VECT_COMMA)  /* V(x), ... */

typedef struct { double VECT_LIST(); } coo_t;

//---------------------------------------------------------------------------------------


coo_t Lv;


// Minimum distance in the periodic system:

//#define MIN_RIJ(x) ( FX= fabs(rij.x),(FX<Lv.x-FX)?rij.x:(rij.x-((rij.x >0)?Lv.x:-Lv.x) ) )
#define MIN_RIJ(x) (rij.x)

#define PW2(x) (x*x)

static inline double sign(double a,double b) { return a= fabs(a),(b<0)?-a:a; }


//---------------- Distance of two rods: -------------------------------------

double dist2_rods(coo_t r1, coo_t r2, coo_t w1, coo_t w2,double lh1,double lh2)
{
 coo_t rij= { VECT_OLIST(r2,-,r1) };
 register double FX;
 coo_t min_rij= { VECT_LIST(MIN_RIJ) };
 double
  xla,xmu,
  rr= VECT_NORM2(min_rij),
  rw1= VECT_PROD(min_rij,w1),
  rw2= VECT_PROD(min_rij,w2),
  w1w2= VECT_PROD(w1,w2),
  cc= 1-PW2(w1w2);

// Checking whether the rods are or not parallel:
// The original code is modified to have symmetry:

 if(cc<1e-15) {
  if(rw1 && rw2) {
   xla= rw1/2;
   xmu= -rw2/2;
  }
  else return rr;
 }

 else {

// Step 1

  xla= (rw1-w1w2*rw2)/cc;
  xmu= (-rw2+w1w2*rw1)/cc;
 }

// Step 2

if( fabs(xla)>lh1 || fabs(xmu)>lh2 ) {

// Step 3 - 7

  if(fabs(xla)-lh1>fabs(xmu)-lh2) {
   xla= sign(lh1,xla);
   xmu= xla*w1w2-rw2;
   if( fabs(xmu)>lh2 ) xmu= sign(lh2,xmu);
  }
  else {
   xmu= sign(lh2,xmu);
   xla= xmu*w1w2+rw1;
   if( fabs(xla)>lh1 ) xla= sign(lh1,xla);
  }
 }

// Step 8

 return rr+PW2(xla)+PW2(xmu) + 2*(xmu*rw2 -xla*(rw1+xmu*w1w2));
}
double check_spherocyl(double CiCj[3], double D, double Lc, double Di[2][3], double *Ci, double *ni, double Dj[2][3], double *Cj, double *nj, int *rim)
{
  coo_t r1, r2, w1, w2;
  double sum, d, normDiCj, normDjCi, DiCj[3], DjCi[3], Ui[3], Uj[3], DjUi[3], DiUj[3], DjCini, DiCjnj;
  int kk, j1, j2;

  r1.x = Ci[0];
  r1.y = Ci[1];
  r1.z = Ci[2];
  r2.x = Cj[0];
  r2.y = Cj[1];
  r2.z = Cj[2];
  w1.x = ni[0];
  w1.y = ni[1];
  w1.z = ni[2];
  w2.x = nj[0];
  w2.y = nj[1];
  w2.z = nj[2];

#ifdef MD_LXYZ
  Lv.x = L[0];
  Lv.y = L[1];
  Lv.z = L[2];
#else
  Lv.x = Lv.y = Lv.z = L;
#endif

  for (j1=0; j1 < 2; j1++)
    for (j2=0; j2 < 2; j2++)
      {
	sum=0.0;
	for (kk=0; kk < 3; kk++)
	  sum += Sqr(Di[j1][kk]-Dj[j2][kk]);
	if (sum < Sqr(D))
	  {
	    //printf("qui -1\n");
	    return -1;
	  }
      }

  for (j2=0; j2 < 2; j2++)
    {
      for (kk=0; kk < 3; kk++)
	DjCi[kk] = Dj[j2][kk] - Ci[kk];
      normDjCi = calc_norm(DjCi);
      DjCini = scalProd(DjCi,ni);
      for (kk=0; kk < 3; kk++)
	{
	  Ui[kk] = Ci[kk] + DjCini*ni[kk];
	  DjUi[kk] = Dj[j2][kk] - Ui[kk];
	}
      if (calc_norm(DjUi) < D && fabs(DjCini) <= Lc*0.5)
	{
	  //printf("qui0\n");
	  return -1;
	}
    }
  for (j1=0; j1 < 2; j1++)
    {
      for (kk=0; kk < 3; kk++)
	DiCj[kk] = Di[j1][kk] - Cj[kk];
      normDiCj = calc_norm(DiCj);
      DiCjnj = scalProd(DiCj,nj);
      for (kk=0; kk < 3; kk++)
	{
	  Uj[kk] = Cj[kk] + DiCjnj*nj[kk];
	  DiUj[kk] = Di[j1][kk] - Uj[kk];
	}
      if (calc_norm(DiUj) < D && fabs(DiCjnj) <= Lc*0.5)
	{
	  //printf("SC qui1\n");
	  return -1;
	}
    }

  *rim = 1;

  if ((d=dist2_rods(r1, r2, w1, w2, Lc*0.5, Lc*0.5)) <= Sqr(D)) 
    {
      *rim=-1;
      return -1;
    }

  
  return 1;
}
#endif
#endif
#endif
