#include<mdsimul.h>
//#define DEBUG_HCMC
#ifdef MC_BOND_POS
void update_spot_pos(int i, double dx, double dy, double dz);
void update_spot_pos_dir(int i, int dir, double fact);
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
extern int rejectMove, checkMoveKF;
#endif
#ifdef MC_BOND_POS
extern int *inCellBP[3], *cellListBP, cellsxBP, cellsyBP, cellszBP;
extern int totspots;
#ifdef MC_BOUNDING_SPHERES
extern int *inCellBS[3], *cellListBS, cellsxBS, cellsyBS, cellszBS;
extern int totspotsBS;
#ifdef MC_OPT_CHECK_BS
extern int *numcheckBS, **checkBS;
#ifdef MC_OPT_ALMARZA
extern int *numcheckBS_MC, **checkBS_MC;
extern int *numcheckBS_P, **checkBS_P;
#endif
#else
extern int *checkBS;
#endif
#endif
#endif
#ifdef MC_BOND_POS
extern struct n2sp_struct *n2sp_map;
extern int **sp2n_map;
#endif
#ifdef MD_SPOT_GLOBAL_ALLOC
extern void BuildAtomPos(int i, double *rO, double **R, double **rat);
#else
extern void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3]);
#endif
#if defined(MC_HC)
const double saxfactMC[3]={0.999,0.7071,0.7071};
#elif defined(MC_SWELL)
const double saxfactMC[3]={0.85,0.37249,0.37249};
#else
const double saxfactMC[3]={0.85,0.68,0.68};
#endif
#ifdef MC_QUASI_CUBE
const double saxfactMC_QC[3]={0.832,0.832,0.832};
#endif
const int nfons=200;
extern void init_rng(int mdseed, int mpi, int my_rank);
#ifdef MC_SIMUL
int dostorebump=0;
#ifdef MC_STORELL
int *cellListMC;
int cellsxMC, cellsyMC, cellszMC;
#endif
#ifdef MC_STOREBONDS
#ifdef MD_LL_BONDS
long long int **bondsMC;
int *numbondsMC;
#else
int *numbondsMC, **bondsMC;
#endif
#endif
#if defined(MC_CLUSTER_NPT) || defined(MC_CLUSTER_MOVE) || defined(MC_NPT_XYZ)
extern int *color, *color_dup, *clsdim, *nbcls, *clsarr, *firstofcls, numOfClusters;  
extern double *Dxpar, *Dypar, *Dzpar, *Dxcls, *Dycls, *Dzcls, *clsCoM[3];
#ifdef MC_ALMARZA
extern int **bondsYN, **bondsPYN;
#ifdef MC_OPT_ALMARZA
extern int **bondsYN_BS;
#endif
#endif
#ifdef MC_NEW_PERC
#ifdef MD_LL_BONDS
extern long long int **bondsP;
extern int *numbondsP;
#else
extern int *numbondsP, **bondsP;
#endif
extern int *refind_bondsP;
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
#ifdef MC_CLUSTER_MOVE
extern double ***RoldAll;
#endif
#ifdef MC_BOND_POS
double **bpos[3], **bposold[3];
#endif
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
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
extern long long int *bondscache2;
#endif
extern int *numbonds;
#else
extern int *bondscache, *numbonds, **bonds;
#if defined(MC_KERN_FRENKEL)||defined(MC_GAPDNA)
extern int *bondscache2;
#endif
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
#if defined(MC_CLUSTER_NPT) || defined(MC_CLUSTER_MOVE) || defined(MC_NPT_XYZ)
#if 0
struct cluster_sort_struct { 
  int dim;
  int color;
};
struct cluster_sort_struct *cluster_sort;
#endif
int compare_func (const void *aa, const void *bb)
{
  int *a, *b;
  int temp;
  a = (int*) aa;
  b = (int*) bb;
  temp = *b - *a;
  if (temp < 0)
    return 1;
  else if (temp > 0)
    return -1;
  else
    return 0;
}

void change_all_colors(int NP, int* color, int colorsrc, int colordst)
{
  int ii;
  for (ii = 0; ii < NP; ii++)
    {
      if (color[ii] == colorsrc)
	color[ii] = colordst;
    }
}

int findmaxColor(int NP, int *color)
{
  int i, maxc=-1;
  for (i = 0; i < NP; i++) 
    {
      if (color[i] > maxc)
	maxc = color[i];
    }
  return maxc;
}

#if 0
struct freecolstk
{
  int idx;
  int *arr;
  int ncls;
};
typedef struct freecolstk freecolT;
freecolT fcstack;

void init_freecolor(freecolT *stack, int nitems)
{
  if (stack->arr==NULL)
    stack->arr = malloc(sizeof(int)*nitems);
  stack->idx=0;
}
void push_freecolor(freecolT *stack, int col)
{
  int oldcol, ii;

  /* NOTA stack->idx da l'indice dello stack che è libero */
  stack->arr[stack->idx] = col;
  if (stack->idx > 0)
    {
      ii=stack->idx;
      while (stack->arr[ii] > stack->arr[ii-1])
	{
	  /* keep array ordered */
	  //printf("qui?!? col=%d stack=%d -2=%d\n", stack->arr[ii], stack->arr[ii-1], stack->arr[ii-2]);
	  /* scambia il valore attuale con quello precedente più piccolo */
	  oldcol = stack->arr[ii-1];
	  stack->arr[ii-1] = stack->arr[ii];
	  stack->arr[ii] = oldcol;
	  ii--;
	}
    }
  (stack->idx)++; 
  //printf("push idx=%d\n", stack->idx);
}

int pop_freecolor(freecolT *stack)
{
  int val;
  /* NOTA stack->idx da l'indice dello stack che è libero */
  (stack->idx)--;
  val = stack->arr[stack->idx];
  //printf("pop idx=%d\n", stack->idx);
  return val;
}
#endif

#ifdef MC_NEW_PERC
const int images_array[8][3]={{0,0,0},
{1,0,0},{0,1,0},{0,0,1},{1,1,1},{1,1,0},{0,1,1},{1,0,1}};
void choose_image(int img, int *dix, int *diy, int *diz)
{
  *dix = images_array[img][0];
  *diy = images_array[img][1];
  *diz = images_array[img][2];
}

void clone_bonds(void)
{
  int i, k;
#ifdef MD_LL_BONDS
  int nb, nn;
  long long int jj, jj2, aa, bb, bo;
#else
  int nb, jj, jj2, kk, nn, aa, bb, bo;
#endif
  int numimg; 
  for (i=0; i < Oparams.parnum*8; i++)
    {
      numbondsP[i] = numbonds[i%Oparams.parnum];
      numimg = i / Oparams.parnum;
      for (k=0; k < numbonds[i%Oparams.parnum]; k++)
	{
	  jj = bonds[i%Oparams.parnum][k] / (NANA);
	  jj2 = bonds[i%Oparams.parnum][k] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;

#ifdef MD_LL_BONDS
	  bondsP[i][k] = (jj+numimg*Oparams.parnum)*(((long long int)NA)*NA)+aa*((long long int)NA)+bb;
#else
	  bondsP[i][k] = (jj+numimg*Oparams.parnum)*(NANA)+aa*NA+bb;
#endif
	}
    }
}
extern int *colorP, *color_dupP;
extern int cellsxP, cellsyP, cellszP;
extern int *inCellP[3], cellsxP, cellsyP, cellsyP, *cellListP;
void build_linked_list_perc(void)
{
  int kk, img, j, n, np, dix, diy, diz;
  double L2big[3], Lbig[3]; 
  for (kk=0; kk < 3; kk++) 
    {
      Lbig[kk] = L[kk]*2.0;
      L2big[kk] = 0.5 * Lbig[kk];
    }
  for (j = 0; j < cellsxP*cellsyP*cellszP + Oparams.parnum*8; j++)
    cellListP[j] = -1;
  //printf("cells=%d %d %d Lbig=%.15G L=%.15G\n", cellsx, cellsy, cellsz, Lbig, L);
  for (n = 0; n < Oparams.parnum*8; n++)
    {
      if (!refind_bondsP[n%Oparams.parnum])
	continue;
      img = n / Oparams.parnum;
      choose_image(img, &dix, &diy, &diz);
      np = n % 8;
      inCellP[0][n] =  (rx[np] + dix*L[0] + L2[0]) * cellsxP / Lbig[0];
      inCellP[1][n] =  (ry[np] + diy*L[1] + L2[1]) * cellsyP / Lbig[1];
      inCellP[2][n] =  (rz[np] + diz*L[2] + L2[2]) * cellszP / Lbig[2];
      if (inCellP[0][n] == cellsxP)
	inCellP[0][n]=cellsxP-1;
      if (inCellP[1][n] == cellsyP)
	inCellP[1][n]=cellsyP-1;
      if (inCellP[2][n] == cellszP)
	inCellP[2][n]=cellszP-1;
      if (inCellP[0][n] == -1)
	inCellP[0][n]=0;
      if (inCellP[1][n] == -1)
	inCellP[1][n]=0;
      if (inCellP[2][n] == -1)
	inCellP[2][n]=0;
      j = (inCellP[2][n]*cellsyP + inCellP[1][n])*cellsxP + 
	inCellP[0][n] + Oparams.parnum*8;
      cellListP[n] = cellListP[j];
      cellListP[j] = n;
    }
}
int boundP(int na, int n, int a, int b)
{
  int i;
  for (i = 0; i < numbondsP[na]; i++)
#ifdef MD_LL_BONDS
    //if (na==1841 && n==8965 && a==2 && b==2)
    //	  printf("[BOUND] TUTTI (%d,%d)-(%d,%d) bonds=%lld\n", na, a, n, b, bonds[na][i]);

    if (bondsP[na][i] == n*(((long long int)NA)*NA)+a*((long long int)NA)+b)
      {
	//printf("[BOUND] (%d,%d)-(%d,%d) bonds=%lld\n", na, a, n, b, bonds[na][i]);
       	return 1;
      }
#else
    if (bondsP[na][i] == n*(NANA)+a*NA+b)
      return 1;
#endif
  return 0;
}
void add_bondP(int na, int n, int a, int b)
{
  if (boundP(na, n, a, b))
    {
      //printf("il bond (%d,%d),(%d,%d) esiste gia'!\n", na, a, n, b);
      return;
    }
#ifdef MD_LL_BONDS
  bondsP[na][numbondsP[na]] = n*(((long long int)NA)*NA)+a*((long long int)NA)+b;
#else
  bondsP[na][numbondsP[na]] = n*(NANA)+a*NA+b;
#endif
  numbondsP[na]++;
  //check_bonds_size(na); 
}
extern int *mapbondsaFlex, *mapbondsbFlex, nbondsFlex;
extern double *dists;
extern void remove_bond(int na, int n, int a, int b);
extern double calcDistNegSP(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, 
		   double *dists, int bondpair);
void store_bonds_mc(int ip);
void restore_bonds_mc(int ip);
int get_image(int i, int j, int k)
{
  int ii;
  for (ii=0; ii < 8; ii++)
    {
      if ((images_array[ii][0]==i)&&
	  (images_array[ii][1]==j)&&
	  (images_array[ii][2]==k))
	return ii;
    }
}
int is_percolating(int ncls)
{
  /* nel caso di catene rigide: cluster di lunghezza l 
     è percolante se e solo se il numero totale di bond è pari a 2*l */
#ifdef MD_LL_BONDS
  int nb, nn;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, nn, aa, bb;
#endif  
  int NP, kk, i, j, k, percola, i2, j2, a, nc, imgi2, imgj2, nbonds, amin, bmin, curcolor, nclsP;
  double dist, dx, dy, dz, Dx, Dy, Dz, shift[3];
  int check, jX, jY, jZ, iX, iY, iZ, jold, dix, diy, diz, djx, djy, djz, ip;
  double Drx, Dry, Drz, imgx, imgy, imgz, Lbig[3];
  int torefind, numcolors, numimg_i, numimg_jj;

  check=0;
  for (nc=0; nc < ncls; nc++)
    {
      if (nbcls[nc] == clsdim[nc]*2)
	{
  	  check=1;
	} 
    }
  //printf("check:%d\n", check);
  if (!check)
    return 0;
  //store_bonds_mc(-1);
  torefind=0; 
  for (i=0; i < Oparams.parnum*8; i++)
    {
      numbondsP[i] = numbonds[i%Oparams.parnum];
      numimg_i = i / Oparams.parnum;
      for (k=0; k < numbonds[i%Oparams.parnum]; k++)
	{
	  jj = bonds[i%Oparams.parnum][k] / (NANA);
	  jj2 = bonds[i%Oparams.parnum][k] % (NANA);
	  Drx = L[0]*rint((rx[i%Oparams.parnum]-rx[jj])/L[0]);
	  Dry = L[1]*rint((ry[i%Oparams.parnum]-ry[jj])/L[1]);
	  Drz = L[2]*rint((rz[i%Oparams.parnum]-rz[jj])/L[2]); 
	  if (fabs(Drx) > L[0]*0.5)
	    djx = 1;
	  else
	    djx = 0;
	  if (fabs(Dry) > L[1]*0.5)
	    djy = 1;
	  else
	    djy = 0;
	  if (fabs(Drz) > L[2]*0.5)
	    djz = 1;
	  else 
	    djz = 0;
#if 0
  	  if (dix || diy || diz)
	    printf("di=%d %d %d\n", dix, diy, diz);
#endif
	  choose_image(numimg_i, &dix, &diy, &diz);
	  numimg_jj = get_image((dix+djx)%2, (diy+djy)%2, (diz+djz)%2);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
#if 0
#ifdef MC_ALMARZA
	  if (aa > 1 && bb > 1)
	    bondsPYN[i][k]=1;
	  else
	    bondsPYN[i][k]=0;
#endif
#endif
#ifdef MD_LL_BONDS
	  bondsP[i][k] = (jj+numimg_jj*Oparams.parnum)*(((long long int)NA)*NA)+aa*((long long int)NA)+bb;
#else
	  bondsP[i][k] = (jj+numimg_jj*Oparams.parnum)*(NANA)+aa*NA+bb;
#endif
	}
#ifdef MC_OPT_ALMARZA
      numcheckBS_P[i] = numcheckBS[i%Oparams.parnum];
      for (k=0; k < numcheckBS[i%Oparams.parnum]; k++)
	{
	  jj = checkBS[i%Oparams.parnum][k];
	  Drx = L[0]*rint((rx[i%Oparams.parnum]-rx[jj])/L[0]);
	  Dry = L[1]*rint((ry[i%Oparams.parnum]-ry[jj])/L[1]);
	  Drz = L[2]*rint((rz[i%Oparams.parnum]-rz[jj])/L[2]); 
	  if (fabs(Drx) > L[0]*0.5)
	    djx = 1;
	  else
	    djx = 0;
	  if (fabs(Dry) > L[1]*0.5)
	    djy = 1;
	  else
	    djy = 0;
	  if (fabs(Drz) > L[2]*0.5)
	    djz = 1;
	  else 
	    djz = 0;
	  choose_image(numimg_i, &dix, &diy, &diz);
	  numimg_jj = get_image((dix+djx)%2, (diy+djy)%2, (diz+djz)%2);
	  checkBS_P[i][k] = jj+numimg_jj*Oparams.parnum;
	}	
#endif
    }


#if 0
  for (i=0; i < Oparams.parnum; i++)
    {
      refind_bondsP[i] = 0;
      for (k=0; k < numbonds[i]; k++)
	{
	  jj = bondsMC[i][k] / (NANA);
	  jj2 = bondsMC[i][k] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  Drx = L[0]*rint((rx[i]-rx[jj])/L[0]);
	  Dry = L[1]*rint((ry[i]-ry[jj])/L[1]);
	  Drz = L[2]*rint((rz[i]-rz[jj])/L[2]); 
	  if (fabs(Drx) > L[0]*0.5 || fabs(Dry) > L[1]*0.5 || fabs(Drz) > L[2]*0.5)
	    {
	      refind_bondsP[i] = 1;
	      torefind++;
	      remove_bond(i, jj, aa, bb);
	    }
	}
    }
  //printf("refind bonds: %d\n", torefind);
  clone_bonds();
  for (kk=0; kk < 3; kk++) 
    {
      Lbig[kk] = L[kk]*2.0;
      //L2big[kk] = 0.5 * Lbig[kk];
    }

  build_linked_list_perc();
  curcolor = 0;
  //printf("na=%d clsdim=%d\n", na, cluster_sort[nc].dim);
  for (i2 = 0; i2 < Oparams.parnum*8; i2++)
    {
      if (!refind_bondsP[i2 % Oparams.parnum])
	continue;
      //printf("nc=%d na*NUMREP=%d i2=%d\n",nc, na*NUMREP, i2);
      //printf("curcolor:%d\n", curcolor);
      i = i2 % Oparams.parnum;//dupcluster[i2];
      for (iZ = -1; iZ <= 1; iZ++) 
	{
	  jZ = inCellP[2][i2] + iZ;    
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellszP - 1;    
	    } 
	  else if (jZ == cellszP) 
	    {
	      jZ = 0;    
	    }
	  for (iY = -1; iY <= 1; iY ++) 
	    {
	      jY = inCellP[1][i2] + iY;    
	      if (jY == -1) 
		{
		  jY = cellsyP - 1;    
		} 
	      else if (jY == cellsyP) 
		{
		  jY = 0;    
		}
	      for (iX = -1; iX <= +1; iX ++) 
		{
		  jX = inCellP[0][i2] + iX;    
		  if (jX == -1) 
		    {
		      jX = cellsxP - 1;    
		    } 
		  else if (jX == cellsxP) 
		    {
		      jX = 0;   
		    }
		  j2 = (jZ *cellsyP + jY) * cellsxP + jX + Oparams.parnum*8;
		  for (j2 = cellListP[j2]; j2 > -1; j2 = cellListP[j2]) 
		    {
		      //if (color[j] != cluster_sort[nc].color)
		      //continue;
		      //coppie++;
		      if (!refind_bondsP[j2 % Oparams.parnum])
			continue;
		      j = j2 % Oparams.parnum;//dupcluster[j2];
		      if (j2 <= i2) 
			continue;
		      //dix = diy = diz = 0;
		      //djx = djy = djz = 0;
		      imgi2 = i2 / Oparams.parnum;
		      imgj2 = j2 / Oparams.parnum;
		      //printf("i2=%d j2=%d imgi2=%d imgj2=%d i=%d j=%d\n", i2, j2, imgi2, imgj2, i, j);
		      choose_image(imgi2, &dix, &diy, &diz);
		      choose_image(imgj2, &djx, &djy, &djz);
		      dx = L[0]*(dix-djx);
		      dy = L[1]*(diy-djy);
		      dz = L[2]*(diz-djz);
		      Dx = rx[i] - rx[j] + dx;
		      Dy = ry[i] - ry[j] + dy;
		      Dz = rz[i] - rz[j] + dz;
		      imgx = -Lbig[0]*rint(Dx/Lbig[0]);
		      imgy = -Lbig[1]*rint(Dy/Lbig[1]);
		      imgz = -Lbig[2]*rint(Dz/Lbig[2]);
		      shift[0] = -(imgx + dx);
		      shift[1] = -(imgy + dy);
		      shift[2] = -(imgz + dx);
		      assign_bond_mapping(i,j);
		      dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
		      nbonds = nbondsFlex;

		      //if (dix!=0||diy!=0||diz!=0||djx!=0||djy!=0||djz!=0)
		      //printf("(%d,%d,%d)-(%d,%d,%d)\n", dix, diy, diz, djx, djy, djz);
		      for (nn=0; nn < nbonds; nn++)
			{
			  if (dists[nn] < 0.0 && !boundP(i2, j2, mapbondsaFlex[nn], mapbondsbFlex[nn]))
			    {
			      if ((i/OprogStatus.polylen != j/OprogStatus.polylen) && 
				  mapbondsaFlex[nn] == 1 && mapbondsbFlex[nn] == 1)
				{
				  //rejectMove = 1;
				  continue;
				}
			      add_bondP(i2, j2, mapbondsaFlex[nn], mapbondsbFlex[nn]);
			      add_bondP(j2, i2, mapbondsbFlex[nn], mapbondsaFlex[nn]);
	    		    }
			}
		    }
		}
	    }
	}
    }
#endif
  /* build clusters of replicated systems */
  NP = Oparams.parnum*8;
  curcolor=0;
  numcolors = 0;
  //init_freecolor(&fcstack, NP);
  for (i=0; i < NP; i++)
    {
      colorP[i] = -1;
    }
 
  for (i=0; i < NP; i++)
    {
      if (colorP[i] == -1)
	{
	  colorP[i] = curcolor;
	  numcolors++;
	  //printf("pop i=%d idx=%d col=%d\n", i, fcstack.idx+1, color[i]);
	}
      //printf("numbonds[%d]=%d\n", i, numbonds[i]);	      
      for (j=0; j < numbondsP[i]; j++)
	{
	  jj = bondsP[i][j] / (NANA);
#ifdef MC_ALMARZA
#if 1
	  if (bondsYN[i%Oparams.parnum][jj%Oparams.parnum]==0)
	    continue;
#else
	  if (bondsPYN[i][j]==0)
	    continue;
#endif
#endif
	  //printf("i=%d jj=%d\n", i, jj);
	  if (colorP[jj] == -1)
	    colorP[jj] = colorP[i];
	  else
	    {
	      if (colorP[i] < colorP[jj])
		{
		  //printf("1) color[%d]=%d to color[%d]=%d\n", jj, color[jj], i, color[i]);
		  //printf("push 1) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], jj, fcstack.idx, color[jj]);
		  //push_freecolor(&fcstack, color[jj]);
		  change_all_colors(NP, colorP, colorP[jj], colorP[i]);
		  numcolors--;
		}	
	      else if (colorP[i] > colorP[jj])
		{
		  //printf("2) color[%d]=%d to color[%d]=%d\n", i, color[i], jj, color[jj]);
		 // printf("push 2) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], fcstack.idx, jj, color[jj]);
		  //push_freecolor(&fcstack, color[i]);
		  change_all_colors(NP, colorP, colorP[i], colorP[jj]);
		  numcolors--;
		}
	    }
	}
#if defined(MC_ALMARZA) && defined(MC_OPT_ALMARZA)
      //printf("numbonds[%d]=%d\n", i, numbonds[i]);	      
      for (j=0; j < numcheckBS_P[i]; j++)
	{
	  jj = checkBS_P[i][j];
	  if (bondsYN_BS[i%Oparams.parnum][jj%Oparams.parnum]==0)
	    continue;
	  //printf("i=%d jj=%d\n", i, jj);
	  if (colorP[jj] == -1)
	    colorP[jj] = colorP[i];
	  else
	    {
	      if (colorP[i] < colorP[jj])
		{
		  //printf("1) color[%d]=%d to color[%d]=%d\n", jj, color[jj], i, color[i]);
		  //printf("push 1) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], jj, fcstack.idx, color[jj]);
		  //push_freecolor(&fcstack, color[jj]);
		  change_all_colors(NP, colorP, colorP[jj], colorP[i]);
		  numcolors--;
		}	
	      else if (colorP[i] > colorP[jj])
		{
		  //printf("2) color[%d]=%d to color[%d]=%d\n", i, color[i], jj, color[jj]);
		 // printf("push 2) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], fcstack.idx, jj, color[jj]);
		  //push_freecolor(&fcstack, color[i]);
		  change_all_colors(NP, colorP, colorP[i], colorP[jj]);
		  numcolors--;
		}
	    }
	}
#endif

      curcolor = findmaxColor(NP, colorP)+1;
      //printf("curcolor=%d\n", curcolor);
    }
  nclsP = numcolors;
   
  /* considera la particelle singole come cluster da 1 */
  for (i = 0; i < NP; i++)
    {
      if (colorP[i]==-1)
	{	    
	  colorP[i] = curcolor;//pop_freecolor(&fcstack);
	  curcolor++;
	  numcolors++;
	}
    }
#if 0
  memcpy(color_dupP,colorP,sizeof(int)*NP);
  qsort(color_dupP,NP,sizeof(int),compare_func);
  nclsP = 1;
  if (color_dupP[0]!=0)
    change_all_colors(NP, colorP, color_dupP[0], 0);
  for (i = 1; i < NP; i++)
    {
      //clsdim[color_dup[i]]++;
      if (color_dupP[i]!=color_dupP[i-1])
	{
	  if (color_dupP[i]!=nclsP)
	    change_all_colors(NP, colorP, color_dupP[i], nclsP);
	  nclsP++;
	}
    } 
  if (nclsP!=numcolors)
    {
      printf("bohboh nclsP=%d numcolors=%d\n", nclsP, numcolors);
    }
#endif
  //printf("ncls2=%d\n", ncNV2);
  if (nclsP < 8*ncls)
    percola = 1;
  else
    percola = 0;
#if 0
  if (percola)
    {
      printf("percolating!\n");
      printf("ncls=%d nclsP=%d\n", ncls, nclsP);
    }
  else
    {
      printf("*NOT* percolating\n");
    }
#endif
  //restore_bonds_mc(-1);

  return percola;
}

int is_cls_percolating(int nc)
{
  /* NOTA: nel caso della cluster move non si usano le BS (come nel caso dell'algoritmo clusterNPT)
   * ma solo i bond dovuti alle patch */
  /* nel caso di catene rigide: cluster di lunghezza l 
     è percolante se e solo se il numero totale di bond è pari a 2*l */
#ifdef MD_LL_BONDS
  int nb, nn;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, nn, aa, bb;
#endif  
  int NP, kk, i, j, k, percola, i2, j2, a, imgi2, imgj2, nbonds, amin, bmin, curcolor, nclsP;
  double dist, dx, dy, dz, Dx, Dy, Dz, shift[3];
  int check, jX, jY, jZ, iX, iY, iZ, jold, dix, diy, diz, djx, djy, djz, ip;
  double Drx, Dry, Drz, imgx, imgy, imgz, Lbig[3];
  int numcolors, numimg_i, numimg_jj, ni;

  if (nbcls[nc] < clsdim[nc]*2)
    {
      return 0;
    } 
    
  //store_bonds_mc(-1);
  //torefind=0; 
  for (i=0; i < Oparams.parnum*8; i++)
    {
      if (color[i%Oparams.parnum]!=nc)
	continue;
      //printf("i=%d\n", i) ;
      numbondsP[i] = numbonds[i%Oparams.parnum];
      numimg_i = i / Oparams.parnum;
      for (k=0; k < numbonds[i%Oparams.parnum]; k++)
	{
	  jj = bonds[i%Oparams.parnum][k] / (NANA);
	  jj2 = bonds[i%Oparams.parnum][k] % (NANA);
	  Drx = L[0]*rint((rx[i%Oparams.parnum]-rx[jj])/L[0]);
	  Dry = L[1]*rint((ry[i%Oparams.parnum]-ry[jj])/L[1]);
	  Drz = L[2]*rint((rz[i%Oparams.parnum]-rz[jj])/L[2]); 
	  if (fabs(Drx) > L[0]*0.5)
	    djx = 1;
	  else
	    djx = 0;
	  if (fabs(Dry) > L[1]*0.5)
	    djy = 1;
	  else
	    djy = 0;
	  if (fabs(Drz) > L[2]*0.5)
	    djz = 1;
	  else 
	    djz = 0;
	  choose_image(numimg_i, &dix, &diy, &diz);
	  numimg_jj = get_image((dix+djx)%2, (diy+djy)%2, (diz+djz)%2);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
#ifdef MD_LL_BONDS
	  bondsP[i][k] = (jj+numimg_jj*Oparams.parnum)*(((long long int)NA)*NA)+aa*((long long int)NA)+bb;
#else
	  bondsP[i][k] = (jj+numimg_jj*Oparams.parnum)*(NANA)+aa*NA+bb;
#endif
	}
    }

#if 0
  for (i=0; i < Oparams.parnum; i++)
    {
      refind_bondsP[i] = 0;
      for (im=0; im < 8; im++)
	{
	  numbondsP[i+im*Oparams.parnum] = numbonds[i];
	}	  
      numimg = i / Oparams.parnum;
      for (k=0; k < numbonds[i]; k++)
	{
	  jj = bonds[i][k] / (NANA);
	  jj2 = bonds[i][k] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  Drx = L[0]*rint((rx[i]-rx[jj])/L[0]);
	  Dry = L[1]*rint((ry[i]-ry[jj])/L[1]);
	  Drz = L[2]*rint((rz[i]-rz[jj])/L[2]); 
	  if (fabs(Drx) > L[0]*0.5 || fabs(Dry) > L[1]*0.5 || fabs(Drz) > L[2]*0.5)
	    {
	      refind_bondsP[i]=1;
	      remove_bond(i, jj, aa, bb);
	      torefind++;
	    }
	}
    }
#endif
  //printf("refind bonds: %d\n", torefind);
#if 0
  clone_bonds();
  for (kk=0; kk < 3; kk++) 
    {
      Lbig[kk] = L[kk]*2.0;
      //L2big[kk] = 0.5 * Lbig[kk];
    }


  build_linked_list_perc();
  curcolor = 0;
  //printf("na=%d clsdim=%d\n", na, cluster_sort[nc].dim);
  for (i2 = 0; i2 < Oparams.parnum*8; i2++)
    {
      if (!refind_bondsP[i2 % Oparams.parnum])
	continue;
      //printf("nc=%d na*NUMREP=%d i2=%d\n",nc, na*NUMREP, i2);
      //printf("curcolor:%d\n", curcolor);
      i = i2 % Oparams.parnum;//dupcluster[i2];
      for (iZ = -1; iZ <= 1; iZ++) 
	{
	  jZ = inCellP[2][i2] + iZ;    
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellszP - 1;    
	    } 
	  else if (jZ == cellszP) 
	    {
	      jZ = 0;    
	    }
	  for (iY = -1; iY <= 1; iY ++) 
	    {
	      jY = inCellP[1][i2] + iY;    
	      if (jY == -1) 
		{
		  jY = cellsyP - 1;    
		} 
	      else if (jY == cellsyP) 
		{
		  jY = 0;    
		}
	      for (iX = -1; iX <= +1; iX ++) 
		{
		  jX = inCellP[0][i2] + iX;    
		  if (jX == -1) 
		    {
		      jX = cellsxP - 1;    
		    } 
		  else if (jX == cellsxP) 
		    {
		      jX = 0;   
		    }
		  j2 = (jZ *cellsyP + jY) * cellsxP + jX + Oparams.parnum*8;
		  for (j2 = cellListP[j2]; j2 > -1; j2 = cellListP[j2]) 
		    {
		      //if (color[j] != cluster_sort[nc].color)
		      //continue;
		      //coppie++;
		      if (!refind_bondsP[j2 % Oparams.parnum])
			continue;
		      j = j2 % Oparams.parnum;//dupcluster[j2];
		      if (j2 <= i2) 
			continue;
		      //dix = diy = diz = 0;
		      //djx = djy = djz = 0;
		      imgi2 = i2 / Oparams.parnum;
		      imgj2 = j2 / Oparams.parnum;
		      //printf("i2=%d j2=%d imgi2=%d imgj2=%d i=%d j=%d\n", i2, j2, imgi2, imgj2, i, j);
		      choose_image(imgi2, &dix, &diy, &diz);
		      choose_image(imgj2, &djx, &djy, &djz);
		      dx = L[0]*(dix-djx);
		      dy = L[1]*(diy-djy);
		      dz = L[2]*(diz-djz);
		      Dx = rx[i] - rx[j] + dx;
		      Dy = ry[i] - ry[j] + dy;
		      Dz = rz[i] - rz[j] + dz;
		      imgx = -Lbig[0]*rint(Dx/Lbig[0]);
		      imgy = -Lbig[1]*rint(Dy/Lbig[1]);
		      imgz = -Lbig[2]*rint(Dz/Lbig[2]);
		      shift[0] = -(imgx + dx);
		      shift[1] = -(imgy + dy);
		      shift[2] = -(imgz + dx);
		      assign_bond_mapping(i,j);
		      dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
		      nbonds = nbondsFlex;

		      //if (dix!=0||diy!=0||diz!=0||djx!=0||djy!=0||djz!=0)
		      //printf("(%d,%d,%d)-(%d,%d,%d)\n", dix, diy, diz, djx, djy, djz);
		      for (nn=0; nn < nbonds; nn++)
			{
			  if (dists[nn] < 0.0 && !boundP(i2, j2, mapbondsaFlex[nn], mapbondsbFlex[nn]))
			    {
			      if ((i/OprogStatus.polylen != j/OprogStatus.polylen) && 
				  mapbondsaFlex[nn] == 1 && mapbondsbFlex[nn] == 1)
				{
				  //rejectMove = 1;
				  continue;
				}
			      add_bondP(i2, j2, mapbondsaFlex[nn], mapbondsbFlex[nn]);
			      add_bondP(j2, i2, mapbondsbFlex[nn], mapbondsaFlex[nn]);
	    		    }
			}
		    }
		}
	    }
	}
    }
#endif
  /* build clusters of replicated systems */
  NP = Oparams.parnum*8;
  curcolor=0;
  numcolors = 0;
  //init_freecolor(&fcstack, NP);
  for (i=0; i < NP; i++)
    {
      colorP[i] = -1;
    }
 
  for (i=0; i < NP; i++)
    {
      if (color[i%Oparams.parnum]!=nc)
	continue;
      //printf("i=%d i mod pn=%d\n", i, i%Oparams.parnum);
      if (colorP[i] == -1)
	{
	  colorP[i] = curcolor;
	  numcolors++;
	  //printf("pop i=%d idx=%d col=%d\n", i, fcstack.idx+1, color[i]);
	}
      //printf("numbonds[%d]=%d\n", i, numbonds[i]);	      
      for (j=0; j < numbondsP[i]; j++)
	{
	  jj = bondsP[i][j] / (NANA);
#ifdef MC_ALMARZA
	  if (bondsYN[i%Oparams.parnum][jj%Oparams.parnum]==0)
	    continue;
#endif
	  //printf("i=%d jj=%d\n", i, jj);
	  if (colorP[jj] == -1)
	    colorP[jj] = colorP[i];
	  else
	    {
	      if (colorP[i] < colorP[jj])
		{
		  //printf("1) color[%d]=%d to color[%d]=%d\n", jj, color[jj], i, color[i]);
		  //printf("push 1) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], jj, fcstack.idx, color[jj]);
		  //push_freecolor(&fcstack, color[jj]);
		  change_all_colors(NP, colorP, colorP[jj], colorP[i]);
		  numcolors--;
		}	
	      else if (colorP[i] > colorP[jj])
		{
		  //printf("2) color[%d]=%d to color[%d]=%d\n", i, color[i], jj, color[jj]);
		 // printf("push 2) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], fcstack.idx, jj, color[jj]);
		  //push_freecolor(&fcstack, color[i]);
		  change_all_colors(NP, colorP, colorP[i], colorP[jj]);
		  numcolors--;
		}
	    }
	}
      curcolor = findmaxColor(NP, colorP)+1;
      //printf("curcolor=%d\n", curcolor);
    }
  nclsP = numcolors;
   
  /* considera la particelle singole come cluster da 1 */
  for (i = 0; i < NP; i++)
    {
      if (colorP[i]==-1 && color[i%Oparams.parnum]==nc)
	{	    
	  colorP[i] = curcolor;//pop_freecolor(&fcstack);
	  curcolor++;
	  numcolors++;
	}
    }
#if 0
  memcpy(color_dupP,colorP,sizeof(int)*NP);
  qsort(color_dupP,NP,sizeof(int),compare_func);
  nclsP = 1;
  if (color_dupP[0]!=0)
    change_all_colors(NP, colorP, color_dupP[0], 0);
  for (i = 1; i < NP; i++)
    {
      //clsdim[color_dup[i]]++;
      if (color_dupP[i]!=color_dupP[i-1])
	{
	  if (color_dupP[i]!=nclsP)
	    change_all_colors(NP, colorP, color_dupP[i], nclsP);
	  nclsP++;
	}
    } 
  if (nclsP!=numcolors)
    {
      printf("bohboh nclsP=%d numcolors=%d\n", nclsP, numcolors);
    }
#endif
  //printf("ncls2=%d\n", ncNV2);
#if 1
  if (nclsP < 8)
    percola = 1;
  else
    percola = 0;
#endif
#if 0
  if (percola)
    {
      printf("percolating!\n");
      printf("clsmove nclsP=%d\n",nclsP);
    }
  else
    {
      printf("*NOT* percolating\n");
    }
#endif
  //restore_bonds_mc(-1);

  return percola;
}
#else
int is_cls_percolating(int nc)
{
  if (nbcls[nc] == clsdim[nc]*2)
    return 1;
}
int is_percolating(int ncls)
{
  int nc;
  for (nc=0; nc < ncls; nc++)
    {
      if (is_cls_percolating(nc))
	{
  	  return 1;
	} 
    }
  return 0;
}
#endif
extern int inside_cluster_move;

void build_clusters(int *Ncls, int *percolating, int check_perc)
{
  int NP, i, j, ncls, nc, a, curcolor, maxc=-1;
  long long int jj;
  int freecolor, cidx, first;
  NP = Oparams.parnum;
  curcolor=0;

  //init_freecolor(&fcstack, NP);
  for (i=0; i < NP; i++)
    {
      color[i] = -1;
    }
#if 0
  for (i=NP-1; i >= 0; i--)
    {
      push_freecolor(&fcstack, i);
      //printf("i=%d stackidx=%d fcstack=%d\n", i, fcstack.idx-1, fcstack.arr[fcstack.idx-1]);
    }
#endif
  for (i=0; i < NP; i++)
    {
      if (color[i] == -1)
	{
	  color[i] = curcolor;
	  //printf("pop i=%d idx=%d col=%d\n", i, fcstack.idx+1, color[i]);
	}
      //printf("numbonds[%d]=%d\n", i, numbonds[i]);	      
      for (j=0; j < numbonds[i]; j++)
	{
	  jj = bonds[i][j] / (NANA);
#ifdef MC_ALMARZA
#if 0
	  if (bondsYN[i][j]==-1)
	    {
	      printf("bond -1 nooooo!\n");
	      exit(-1);
	    }
#endif
	  if (bondsYN[i][j]==0)
	    continue;	    
#endif
	  //printf("i=%d jj=%d\n", i, jj);
	  if (color[jj] == -1)
	    color[jj] = color[i];
	  else
	    {
	      if (color[i] < color[jj])
		{
		  //printf("1) color[%d]=%d to color[%d]=%d\n", jj, color[jj], i, color[i]);
		  //printf("push 1) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], jj, fcstack.idx, color[jj]);
		  //push_freecolor(&fcstack, color[jj]);
		  change_all_colors(NP, color, color[jj], color[i]);
		}	
	      else if (color[i] > color[jj])
		{
		  //printf("2) color[%d]=%d to color[%d]=%d\n", i, color[i], jj, color[jj]);
		 // printf("push 2) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], fcstack.idx, jj, color[jj]);
		  //push_freecolor(&fcstack, color[i]);
		  change_all_colors(NP, color, color[i], color[jj]);
		}
	    }
	}
#if defined(MC_ALMARZA) && defined(MC_OPT_ALMARZA)
      if (!inside_cluster_move)
	{
	  for (j=0; j < numcheckBS[i]; j++)
	    {
	      jj = checkBS[i][j];
	      if (bondsYN_BS[i][j]==0 || (color[i]!=-1 && color[i] == color[jj]))
		continue;	    
	      //printf("i=%d jj=%d\n", i, jj);
	      if (color[jj] == -1)
		color[jj] = color[i];
	      else
		{
		  if (color[i] < color[jj])
		    {
		      //printf("1) color[%d]=%d to color[%d]=%d\n", jj, color[jj], i, color[i]);
		      //printf("push 1) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], jj, fcstack.idx, color[jj]);
		      //push_freecolor(&fcstack, color[jj]);
		      change_all_colors(NP, color, color[jj], color[i]);
		    }	
		  else if (color[i] > color[jj])
		    {
		      //printf("2) color[%d]=%d to color[%d]=%d\n", i, color[i], jj, color[jj]);
		      // printf("push 2) color[%d]=%d idx=%d col[jj=%d]=%d\n", i, color[i], fcstack.idx, jj, color[jj]);
		      //push_freecolor(&fcstack, color[i]);
		      change_all_colors(NP, color, color[i], color[jj]);
		    }
		}
	    }
	}
#endif
      curcolor = findmaxColor(NP, color)+1;
      //printf("curcolor=%d\n", curcolor);
    }
  /* considera la particelle singole come cluster da 1 */
  for (i = 0; i < NP; i++)
    {
      if (color[i]==-1)
	{	    
	  color[i] = curcolor;//pop_freecolor(&fcstack);
	  curcolor++;
	}
    }
  memcpy(color_dup,color,sizeof(int)*NP);
  qsort(color_dup,NP,sizeof(int),compare_func);
#if 0
  for (i = 0; i < NP; i++)
    {
      printf("color_dup[%d]=%d\n",i, color_dup[i]);
    }
#endif
  ncls = 1;
  if (color_dup[0]!=0)
    change_all_colors(NP, color, color_dup[0], 0);
  for (i = 1; i < NP; i++)
    {
      //clsdim[color_dup[i]]++;
      if (color_dup[i]!=color_dup[i-1])
	{
	  if (color_dup[i]!=ncls)
	    change_all_colors(NP, color, color_dup[i], ncls);
	  ncls++;
	}
    } 
#if 0
  for (i=0; i < Oparams.parnum; i++)
    printf("color[%d]=%d\n", i, color[i]);
  exit(-1);
#endif
  //printf("ncls = %d\n", ncls);
  for (nc = 0; nc < ncls; nc++)
    {
      clsdim[nc] = 0; 
      /* nbcls è il numero di bond nel cluster color[a] */
      nbcls[nc] = 0;
    }
  for (nc = 0; nc < ncls; nc++)
    {
      for (a = 0; a < NP; a++)
	{
	  if (color[a] == nc)
	    {
	      clsdim[color[a]]++;
	      //clscol[nc] = color[a];
	      nbcls[color[a]] += numbonds[a];
	    }
	}
    }
  cidx=0;
  for (nc=0; nc < ncls; nc++)
    {
      first=1;
      for (a=0; a < NP; a++)
	{
	  if (color[a]==nc)
	    {
	      clsarr[cidx] = a;
	      if (first)
		firstofcls[nc] = cidx; 
	      cidx++;
	      first = 0;
	    }
	}
    }
#if 0
  //for (i = 0; i < NP; i++)
    //printf("clsarr[%d]=%d\n", i, clsarr[i]);
  saveBakAscii("first.cnf");
  for (nc=0; nc < ncls; nc++)
    printf("cls dim[%d]=%d firstofcls=%d\n", nc, clsdim[nc], firstofcls[nc]);
  printf("ncls=%d\n", ncls);
  exit(-1);
#endif
  //printf("ncls=%d\n", ncls);
  if (check_perc)
    *percolating = is_percolating(ncls);
  *Ncls = ncls;
  //exit(-1);
}
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
#ifdef MC_SWELL
      mbox[tt][0].dr[0]=0.0;
      mbox[tt][0].dr[1]=0.0;
      //mbox[tt][1].dr[2]=(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][0].dr[2]=0.0;
      mbox[tt][0].sa[0]=0.85*sa[0];
      mbox[tt][0].sa[1]=0.48*sa[1]; 
      mbox[tt][0].sa[2]=0.21702*sa[2];

      mbox[tt][1].dr[0]=0;
      mbox[tt][1].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][1].dr[2]=0.0;
      mbox[tt][1].sa[0]=0.85*sa[0];
      mbox[tt][1].sa[1]=0.21702*sa[1]; 
      mbox[tt][1].sa[2]=0.48*sa[2];

      /* con questo ultimo multibox si evitano configurazioni molto overlapped
	 con asse x quasi parallelo nel caso di grandi elongazioni */
      mbox[tt][2].dr[0]=0;
      mbox[tt][2].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][2].dr[2]=0.0;
      mbox[tt][2].sa[0]=0.94*sa[0];
      mbox[tt][2].sa[1]=0.24124*sa[1]; 
      mbox[tt][2].sa[2]=0.24124*sa[2];
     
      mbox[tt][3].dr[0]=0;
      mbox[tt][3].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][3].dr[2]=0.0;
      mbox[tt][3].sa[0]=0.99*sa[0];
      mbox[tt][3].sa[1]=0.099749*sa[1]; 
      mbox[tt][3].sa[2]=0.099749*sa[2];

#if 1
      mbox[tt][4].dr[0]=0;
      mbox[tt][4].dr[1]=0;
      //mbox[tt][1].dr[2]=-(saxfactMC[2]+0.1*0.5)*sa[2];
      mbox[tt][4].dr[2]=0.0;
      mbox[tt][4].sa[0]=0.999*sa[0];
      mbox[tt][4].sa[1]=0.031614*sa[1]; 
      mbox[tt][4].sa[2]=0.031614*sa[2];
#endif
#else
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
#ifdef MC_BOND_POS
  int extraspots = 0;
#endif
  rxold=rx[ip];
  ryold=ry[ip];
  rzold=rz[ip];
  for (k1=0; k1<3; k1++)
    for (k2=0; k2<3; k2++)
      Rold[k1][k2]=R[ip][k1][k2]; 
#ifdef MC_BOND_POS
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < typesArr[typeOfPart[ip]].nspots + extraspots; k2++)
      {
	bposold[k1][ip][k2] = bpos[k1][ip][k2];
      } 
#endif
}
void restore_coord(int ip)
{
#ifdef MC_BOND_POS
  int extraspots = 0;
#endif
  int k1, k2;
  rx[ip]=rxold;
  ry[ip]=ryold;
  rz[ip]=rzold;
  for (k1=0; k1<3; k1++)
    for (k2=0; k2<3; k2++)
      R[ip][k1][k2]=Rold[k1][k2]; 
#ifdef MC_BOND_POS
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < typesArr[typeOfPart[ip]].nspots + extraspots; k2++)
      {
	bpos[k1][ip][k2] = bposold[k1][ip][k2];
      } 
#endif
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
#ifdef MC_CLUSTER_MOVE
long long int totclsrejMC=0, totclsmovesMC=0, rotclsmoveMC=0, traclsmoveMC=0, rotclsrejMC=0, traclsrejMC=0;
#endif
#ifdef MC_BIGROT_MOVE
long long int totbigrotmovMC=0, totbigrotmovrejMC=0;
#endif
void tra_move(int ip)
{
  double dx, dy, dz;
#ifdef MC_BOND_POS
  int k1, extraspots=0;
#endif
  dx = OprogStatus.deltaMC*(ranf()-0.5);
  dy = OprogStatus.deltaMC*(ranf()-0.5);
  dz = OprogStatus.deltaMC*(ranf()-0.5);
  rx[ip]+= dx;
  ry[ip]+= dy;
  rz[ip]+= dz;
#ifdef MC_BOND_POS
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
  for (k1=0; k1 < typesArr[typeOfPart[ip]].nspots + extraspots; k1++)
    {
	bpos[0][ip][k1] += dx;
	bpos[1][ip][k1] += dy;
	bpos[2][ip][k1] += dz;
    }
#endif
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
void rot_move(int ip, int flip)
{
#ifdef MC_BOND_POS
  int extraspots = 0;
#endif
  double theta, thetaSq, sinw, cosw;
  double ox, oy, oz, OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random orientation */
#ifdef MC_ALT_ROT
  double thor, xp[3], xl[3];
  thor = 4.0*acos(0.0)*ranf(); /* random angle between 0 and 2*pi */
  xp[0] = 0.0;
  xp[1] = cos(thor);
  xp[2] = sin(thor);
  body2labR(ip, xp, xl, NULL, R[ip]); 
  ox = xl[0];
  oy = xl[1];
  oz = xl[2];
#else
  orient(&ox,&oy,&oz);
#ifdef MC_BENT_DBLCYL
  if (typesArr[typeOfPart[ip]].nhardobjs == 0)
    remove_parall(ip, &ox, &oy, &oz);
#else
#if !defined(MC_HELIX) && !defined(MC_KERN_FRENKEL) && !defined(MC_SWELL)
  remove_parall(ip, &ox, &oy, &oz);
#endif
#endif
#endif
  /* pick a random rotation angle */
#ifdef MC_FLIP_MOVE
  if (flip==1)
    theta = acos(0.0); /* 90 degree rotation */
  else
    theta= OprogStatus.dthetaMC*(ranf()-0.5);
#else
  theta= OprogStatus.dthetaMC*(ranf()-0.5);
#endif
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
#ifdef MC_BOND_POS
  rA[0] = rx[ip];
  rA[1] = ry[ip];
  rA[2] = rz[ip];
  BuildAtomPos(ip, rA, R[ip], ratA);
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
  for (k1 = 0; k1 < typesArr[typeOfPart[ip]].nspots + extraspots; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	bpos[k2][ip][k1] = ratA[k1+1][k2];
    }
#endif
  rotmoveMC++;
}

int random_move(int ip)
{
  double p;
  //printf("random move ip=%d\n", ip);
  p=ranf();
  if (OprogStatus.restrmove==1|| OprogStatus.restrmove==2
#ifdef MC_RESTR_MATRIX
      || OprogStatus.restrmove==3
#endif
      )
    p=0.0;
  
  if (p <= 0.5)
   {
     tra_move(ip);
     return 0;
   }
  else
    {
#ifdef MC_FLIP_MOVE
      /* flip_prob must be between 0 and 1 
	 (in Blaak, Frenkel and Mulder, J. Chem. Phys., Vol. 110, (1999) 
	 they suggest 0.1=1/10) */
      if (OprogStatus.flip_prob > 0.0 && ranf()< OprogStatus.flip_prob)
	rot_move(ip, 1);
      else
	rot_move(ip, 0);
#else
      rot_move(ip, 0);
#endif
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
#undef MC_OF_BOXES
extern int are_spheres(int i, int j);
double check_overlap(int i, int j, double shift[3], int *errchk);

#ifdef MC_BENT_DBLCYL
void save_pos_R(int i, double xi[3], double Ri[3][3])
{
  int a, b;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      Ri[a][b] = R[i][a][b];
  xi[0] = rx[i]; 
  xi[1] = ry[i];
  xi[2] = rz[i];
}

void load_pos_R(int i, double xi[3], double Ri[3][3])
{
  int a, b;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      R[i][a][b] = Ri[a][b];
  rx[i] = xi[0]; 
  ry[i] = xi[1];
  rz[i] = xi[2];
}
void set_pos_R_ho(int i, int a)
{
  int k1, k2, k3;
  double Rho[3][3], xl[3], rA[3];
  /* NOTA: La posizione dell'hard object è riferita al laboratory
     reference system e così pure l'orientazione */
  for (k1 = 0; k1 < 3; k1++)
    {
      typesArr[typeOfPart[i]].sax[k1] = typesArr[typeOfPart[i]].hardobjs[a].sax[k1]; 
    }
  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  
  body2lab(i, typesArr[typeOfPart[i]].hardobjs[a].x, xl, rA, R[i]);
  rx[i] = xl[0];
  ry[i] = xl[1];
  rz[i] = xl[2];
  /* i campi hardobjs[].n[] vengono usati come orientazioni degli HC */
  /* multiply Rho=R[i]*hardobjs[a].R */
  /* la matrice di rotazione del hard object è una rotazione ulteriore 
     rispetto al body reference system */
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  Rho[k1][k2] = 0.0;//R[i][k1][k2];
#if 1
	  //printf("(%d,%d) %f\n", k1, k2,typesArr[typeOfPart[i]].hardobjs[a].R[k1][k2]);
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      /* matrix multiplication: riga * colonna */
	      Rho[k1][k2] += typesArr[typeOfPart[i]].hardobjs[a].R[k1][k3]*R[i][k3][k2];
	    }  
#endif
	}
    }
  for (k1 =0; k1 < 3; k1++)
    {
      for (k2 =0; k2 < 3; k2++)
	{
	  R[i][k1][k2] = Rho[k1][k2];
	}
    }
  /* first apply rotatation then transform the traslation which is given in the body 
     reference system to laboratory coordinates */	  
  /* hard object position is relative to the body reference system */
#if 0
  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  
  body2lab(i, typesArr[typeOfPart[i]].hardobjs[a].x, xl, rA, R[i]);
  rx[i] = xl[0];
  ry[i] = xl[1];
  rz[i] = xl[2];
#endif
}
#endif
#ifdef MC_PERWER
void tRDiagRpw(int i, double M[3][3], double D[3], double **Ri)
{
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = D[0];
  Di[1][1] = D[1];
  Di[2][2] = D[2];
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
void xlambda(double lambda, double rA[3], double A[3][3], double rB[3], double B[3][3], double x[3])
{
  double lamA[3][3], onemlamB[3][3], ABL[3][3], invABL[3][3];
  double x1[3], x2[3], x3[3], detinvABL;
  int k1, k2;
  /* calcola xlambda, vedi L. Paramonov and S. N. Yaliraki J. Chem. Phys. 123, 194111 (2005) */
  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  lamA[k1][k2] = lambda*A[k1][k2];
	  onemlamB[k1][k2] = (1.0-lambda)*B[k1][k2];
 	  ABL[k1][k2] = lamA[k1][k2] + onemlamB[k1][k2];
	}
    }
  for (k1=0; k1 < 3; k1++)
    {
      x1[k1]=0;
      x2[k1]=0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  x1[k1] += lamA[k1][k2]*rA[k2];
	  x2[k1] += onemlamB[k1][k2]*rB[k2];
	}
      x3[k1] = x1[k1] + x2[k1];
    }
  detinvABL=-ABL[0][2]*ABL[1][1]*ABL[2][0] + ABL[0][1]*ABL[1][2]*ABL[2][0] + 
    ABL[0][2]*ABL[1][0]*ABL[2][1] - ABL[0][0]*ABL[1][2]*ABL[2][1] - 
    ABL[0][1]*ABL[1][0]*ABL[2][2] + ABL[0][0]*ABL[1][1]*ABL[2][2]; 

  invABL[0][0] = -ABL[1][2]*ABL[2][1] + ABL[1][1]*ABL[2][2];
  invABL[0][1] =  ABL[0][2]*ABL[2][1] - ABL[0][1]*ABL[2][2];
  invABL[0][2] = -ABL[0][2]*ABL[1][1] + ABL[0][1]*ABL[1][2];
  invABL[1][0] =  ABL[1][2]*ABL[2][0] - ABL[1][0]*ABL[2][2]; /* a12 a20 - a10 a22 */
  invABL[1][1] = -ABL[0][2]*ABL[2][0] + ABL[0][0]*ABL[2][2]; /* -a02 a20 + a00 a22 */
  invABL[1][2] =  ABL[0][2]*ABL[1][0] - ABL[0][0]*ABL[1][2]; /* a02 a10 - a00 a12 */
  invABL[2][0] = -ABL[1][1]*ABL[2][0] + ABL[1][0]*ABL[2][1]; /* -a11 a20 + a10 a21 */
  invABL[2][1] =  ABL[0][1]*ABL[2][0] - ABL[0][0]*ABL[2][1]; /* a01 a20 - a00 a21 */
  invABL[2][2] = -ABL[0][1]*ABL[1][0] + ABL[0][0]*ABL[1][1]; /* -a01 a10 + a00 a11 */

  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      invABL[k1][k2] /= detinvABL;

  for (k1 = 0; k1 < 3; k1++)
    {
      x[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  x[k1] += invABL[k1][k2]*x3[k2];
	}
    }
}

double Slam(double lambda, double rA[3], double A[3][3], double rB[3], double B[3][3])
{
  int k1, k2;
  double xlam[3], fA[3], fB[3], SA, SB;

  xlambda(lambda, rA, A, rB, B, xlam);

  for (k1=0; k1 < 3; k1++)
    {
      fA[k1] = 0;
      fB[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  fA[k1] += A[k1][k2]*(xlam[k2]-rA[k2]);
	  fB[k1] += B[k1][k2]*(xlam[k2]-rB[k2]);
	}
    }

  SA = SB = 0.0;
  for (k1=0; k1 < 3; k1++)
    {
      SA += lambda*(xlam[k1]-rA[k1])*fA[k1];
      SB += (1.0-lambda)*(xlam[k1]-rB[k1])*fB[k1];
    }
  /* ho messo un - così la funzione ha un minimo invece
     che un massimo e questo minimo viene trovato dalla funzione brentPW */
  return -(SA+SB);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
int brentPWTooManyIter=0;
double brentPW(double ax, double bx, double cx, double tol, double *xmin, double rA[3], double A[3][3], double rB[3], double B[3][3])
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
 * (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)),
 * this routine isolates the minimum to a fractional precision of about tol using Brent's
 * method. The abscissa of the minimum is returned as xmin, and the minimum function value 
 * is returned as brent, the returned function value. */
{ 
  int iter, ITMAXBR=100;
  const double CGOLD=0.3819660;
  const double ZEPSBR=1E-20;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0, fuold;
  brentPWTooManyIter=0;
  /* This will be the distance moved on the step before last.*/
  a=(ax < cx ? ax : cx); /*a and b must be in ascending order, 
			   but input abscissas need not be.*/
  b=(ax > cx ? ax : cx);
  x=w=v=bx; /*Initializations...*/
  fw=fv=fx=Slam(x, rA, A, rB, B); 
  if (fw < -1.0)
    {
      /* non-overlap! */
      *xmin=x;
      return -100.0;
    }
  fuold = fv;
  for (iter=1;iter<=ITMAXBR;iter++)
    { 
      /*Main program loop.*/
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPSBR); 
      if (fabs(x-xm) <= (tol2-0.5*(b-a)))
	{ /*Test for done here.*/
	  *xmin=x;
	  return fx;
	} 
      if (fabs(e) > tol1) 
	{ /*Construct a trial parabolic fit.*/
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if (q > 0.0)
	    p = -p; 
	  q=fabs(q);
	  etemp=e; 
	  e=d; 
	  if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=CGOLD*(e=(x >= xm ? a-x : b-x)); 
	    /*The above conditions determine the acceptability of the parabolic fit.
	     * Here we take the golden section step into the larger of the two segments.*/
	  else
	    {
	      d=p/q; /* Take the parabolic step.*/
	      u=x+d; 
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x); 
	    }
	}
      else
	{
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	} 
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=Slam(u, rA, A, rB, B); /*This is the one function evaluation per iteration.*/
      if (fu < -1.0)
	{
	  /* non overlap! */
	  *xmin=x;
	  return -100.0;
	}
#if 0
      if (2.0*fabs(fuold-fu) <= tol*(fabs(fuold)+fabs(fu)+ZEPSBR)) 
	{ 
	  *xmin=u;
	  return fu;
	}
#endif
      fuold = fu;//
      if (fu <= fx)
	{ /*Now decide what to do with our function evaluation.*/
	  if (u >= x) 
	    a=x;
	  else
	    b=x;
	  SHFT(v,w,x,u); /* Housekeeping follows:*/
	  SHFT(fv,fw,fx,fu); 
	} 
      else
	{ 
	  if (u < x) 
	    a=u; 
	  else 
	    b=u; 
	  if (fu <= fw || w == x)
	    {
	      v=w; w=u; fv=fw; fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w)
	    { 
	      v=u; fv=fu;
	    }
	} /* Done with housekeeping. Back for another iteration.*/
    }
  printf("Too many iterations in brent!\n");
  brentPWTooManyIter=1;
  //nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}

double check_overlap_pw(int i, int j, double shift[3])
{
  const double tolPW=1.0E-12;
  double res, A[3][3], B[3][3], xmin; 
  int k1, k2;
  double  DA[3], DB[3], rA[3], rB[3];
  int typei, typej;

  typei = typeOfPart[i];
  typej = typeOfPart[j];

  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];

  rB[0] = rx[j]+shift[0];
  rB[1] = ry[j]+shift[1];
  rB[2] = rz[j]+shift[2];

  for (k1=0; k1 < 3; k1++)
    {
      DA[k1]= 1.0/Sqr(typesArr[typei].sax[k1]);
      DB[k1]= 1.0/Sqr(typesArr[typej].sax[k1]);
    }
  tRDiagRpw(i, A, DA, R[i]);
  tRDiagRpw(j, B, DB, R[j]);

  res =  - brentPW(0, 0.5, 1.0, tolPW, &xmin, rA, A, rB, B);
  if (brentPWTooManyIter)
    {
      printf("res=%f xmin=%f\n", res, xmin);
      exit(-1);
    }
  //printf("res=%f\n", res);
  return res - 1.0;
}
#endif
double check_overlap_ij(int i, int j, double shift[3], int *errchk)
{
  int k, k1, k2, kk;
  double tpo, vecg[8], vecgNeg[8], rA[3], rB[3], saxi[3], saxj[3];
  double daSq, drSq, d, d0, r1[3], r2[3], alpha; 
  *errchk=0;
  OprogStatus.optnnl = 0;

#if defined(DEBUG_HCMC) && 0
  d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
  //printf("d=%f\n",d);
  return d;
  //exit(-1);
#endif
  if (are_spheres(i,j))
    {
      tpo = OprogStatus.targetPhi;
      OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
      calcdist_retcheck = 0;
      d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
      *errchk = calcdist_retcheck;
      OprogStatus.targetPhi = tpo;
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
#if 1
  /* if bounding spheres do not overlap parallelepiped will not overlap */
  daSq = drSq = 0.0;
  for (kk=0; kk < 3; kk++)
    {
      daSq += Sqr(typesArr[typeOfPart[i]].sax[kk]+typesArr[typeOfPart[j]].sax[kk]);
      drSq += Sqr(rA[kk]-rB[kk]-shift[kk]);
    } 

  //printf("daSq=%f\n", sqrt(daSq));
  if (drSq > daSq)
    return 1.0;
#endif 
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
#ifndef MC_PERWER
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
#endif
  tpo = OprogStatus.targetPhi;
  OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
  calcdist_retcheck = 0;
  d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
  *errchk = calcdist_retcheck;
  OprogStatus.targetPhi = tpo;
#ifndef MC_QUASI_CUBE
  if (*errchk)
    {
      //store_bump(0,1);
      //exit(-1);
      //printf("parnum=%d\n", Oparams.parnum);
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
#ifdef MC_HELIX
#define nmaxNxi 1000
double xihel[nmaxNxi], xhel[nmaxNxi][3], xhelA[nmaxNxi][3], xhelB[nmaxNxi][3];
void build_helix(void);
void mgl_helix(FILE* fs, int i, char *col)
{
  double rA[3], rB[3];
  int Nxi, jj, j1, j2;
  static double sigSq;
  static int first=1;

  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  if (first)
    {
      first = 0;
      build_helix();
    }
  for (jj=0; jj < OprogStatus.Nxi; jj++)
    {
      body2lab(i, xhel[jj], xhelA[jj], rA, R[i]);
      fprintf(fs, "%.15G %.15G %.15G @ %.8G C[%s]\n",xhelA[jj][0], xhelA[jj][1], xhelA[jj][2],
	      OprogStatus.sighelix/2.0, col);
    }
}
void build_helix(void)
{
  double temp, pi, npitch, deltaxi, length_eucl, radius;
  int jj, Nxi, kk;
  double xcm[3];

  radius = OprogStatus.radhelix; /* x,y perpendicular to helix axis in body reference frame */
  length_eucl = OprogStatus.lenhelix; /* helix axis along z in body reference frame */
  Nxi = OprogStatus.Nxi;
  pi = acos(0.0)*2.0;
  npitch = length_eucl/OprogStatus.pitch;
  temp=npitch*sqrt(Sqr(radius)+Sqr(OprogStatus.pitch/(2.0*pi)));
  deltaxi=2.0*pi*npitch/((double)(Nxi-1));
  // length_eucl=OprogStatus.npitch*OprogStatus.pitch;
  for (jj=0; jj < Nxi; jj++)
    {
      xihel[jj]=((double)jj)*deltaxi;
      xhel[jj][0]=radius*cos(xihel[jj]);
      xhel[jj][1]=radius*sin(xihel[jj]);
      xhel[jj][2]=OprogStatus.pitch*xihel[jj]/(2.0*pi);
    }
  xcm[0]=xcm[1]=xcm[2]=0.0;
  for (jj=0; jj < Nxi; jj++)
    {
      for (kk=0; kk < 3; kk++)
	xcm[kk] += xhel[jj][kk]; 
    } 

  for (kk=0; kk < 3; kk++)
    xcm[kk] /= Nxi;
  for (jj=0; jj < Nxi; jj++)
    {
      for (kk=0; kk < 3; kk++)
	xhel[jj][kk] -= xcm[kk];
    }   
}
double check_overlap_helices(int i, int j, double shift[3])
{
  double rA[3], rB[3];
  int jj, j1, j2;
  static double sigSq;
  static int first=1;

  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  rB[0] = rx[j] + shift[0];
  rB[1] = ry[j] + shift[1];
  rB[2] = rz[j] + shift[2];
  if (first)
    {
      first = 0;
      sigSq=Sqr(OprogStatus.sighelix);
      build_helix();
    }
   for (jj=0; jj < OprogStatus.Nxi; jj++)
     {
       body2lab(i, xhel[jj], xhelA[jj], rA, R[i]);
       body2lab(j, xhel[jj], xhelB[jj], rB, R[j]);
     }
   for (j1=0; j1 < OprogStatus.Nxi; j1++)
     for (j2=0; j2 < OprogStatus.Nxi; j2++)
       {
	 if (Sqr(xhelA[j1][0]-xhelB[j2][0])+Sqr(xhelA[j1][1]-xhelB[j2][1])+
	     Sqr(xhelA[j1][2]-xhelB[j2][2]) < sigSq)
	   return -1.0;
       }

   return 1.0;
}
#endif
#ifdef MC_BENT_DBLCYL
double check_overlap_bent_dblcyl(int i, int j, double shift[3], int *errchk)
{
  int a, b;
  double Ri[3][3], Rj[3][3], xi[3], xj[3], d=0.0;
  /* NOTE: check overlaps only for all pairs of hard objects,
     hence in the configuration file two hard cylinders have to be provided
     with given position and orientation */

  for (a=0; a < typesArr[typeOfPart[i]].nhardobjs; a++)
    for (b=0; b < typesArr[typeOfPart[j]].nhardobjs; b++)
      {
	save_pos_R(i, xi, Ri);
	save_pos_R(j, xj, Rj);
	set_pos_R_ho(i, a);
	set_pos_R_ho(j, b);
	d=check_overlap_ij(i, j, shift, errchk);
	load_pos_R(i, xi, Ri);
	load_pos_R(j, xj, Rj);
	if (d < 0)
	  return d;
      }
  return d;
}
#endif
double check_overlap(int i, int j, double shift[3], int *errchk)
{
#if 0
  int k, k1, k2, kk;
  double vecg[8], vecgNeg[8], rA[3], rB[3], saxi[3], saxj[3];
  double daSq, drSq, d, d0, r1[3], r2[3], alpha; 
#endif
#ifdef MC_BENT_DBLCYL
  //printf("nhard %d %d\n", typesArr[typeOfPart[i]].nhardobjs, typesArr[typeOfPart[i]].nhardobjs);
  if (typesArr[typeOfPart[i]].nhardobjs > 0 && typesArr[typeOfPart[j]].nhardobjs > 0)
    return check_overlap_bent_dblcyl(i, j, shift, errchk);
  else
    return check_overlap_ij(i, j, shift, errchk);
#else
  return check_overlap_ij(i, j, shift, errchk);
#endif
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
#if defined(MC_CLUSTER_MOVE)|| defined(MC_CLUSTER_NPT)|| defined(MC_NPT_XYZ)
int clsNPT=0;
#endif
#ifdef MC_BOUNDING_SPHERES
#ifdef MC_OPT_CHECK_BS
int func_checkBS(int i, int j)
{
  int kk;
  for (kk=0; kk < numcheckBS[i]; kk++)
    {
      if (checkBS[i][kk] == j)
       	return 1;
    }
  return 0;
}
#endif
void find_part_with_BS(int i)
{
  int kk, nb, k, iZ, jZ, iX, jX, iY, jY, n, na, ns, ns2, j;
  int cellRangeT[6];
  double shift[3], sigmabs, distSq;
  double rspA[3], rspB[3];
  
  nb=-1;
  for (j=0; j < Oparams.parnum; j++)
    {
#ifdef MC_OPT_CHECK_BS
      numcheckBS[j] = 0;
#else
      checkBS[j] = 0;
#endif
    }
  for (ns = typesArr[typeOfPart[i]].nspots; ns < typesArr[typeOfPart[i]].nspotsBS + typesArr[typeOfPart[i]].nspots; ns++)
    {
      na = sp2n_map[i][ns] - totspots;
      for (k = 0; k < 6; k++) cellRangeT[k] = cellRange[k];
      for (kk = 0;  kk < 3; kk++)
	{
	  cellRange[2*kk]   = - 1;
	  cellRange[2*kk+1] =   1;
	}

      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCellBS[2][na] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellszBS - 1;    
#ifdef MD_LXYZ
	      shift[2] = - L[2];
#else
	      shift[2] = - L;
#endif
	    } 
	  else if (jZ == cellszBS) 
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
	      jY = inCellBS[1][na] + iY;    
	      shift[1] = 0.0;
	      if (jY == -1) 
		{
		  jY = cellsyBS - 1;    
#ifdef MD_LXYZ
		  shift[1] = -L[1];
#else
		  shift[1] = -L;
#endif
		} 
	      else if (jY == cellsyBS) 
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
		  jX = inCellBS[0][na] + iX;    
		  shift[0] = 0.0;
		  if (jX == -1) 
		    {
		      jX = cellsxBS - 1;    
#ifdef MD_LXYZ
		      shift[0] = - L[0];
#else
		      shift[0] = - L;
#endif
		    } 
		  else if (jX == cellsxBS) 
		    {
		      jX = 0;   
#ifdef MD_LXYZ
		      shift[0] = L[0];
#else
		      shift[0] = L;
#endif
		    }
		  n = (jZ *cellsyBS + jY) * cellsxBS + jX + totspotsBS;

		  for (n = cellListBS[n]; n > -1; n = cellListBS[n]) 
		    {
		      if (n != na && n != nb && (nb >= -1 || n < na)) 
			{
			  //assign_bond_mapping(i,j);
			  j = n2sp_map[n+totspots].i;
#ifdef MC_OPT_CHECK_BS
			  if (func_checkBS(i,j))
			    continue;
#else
			  if (checkBS[j])
			    continue;
#endif
			  ns2 = n2sp_map[n+totspots].ns;
#if 0
			  if (ns2 < 2)
			    {
			      printf("we have a problem\n");
			      exit(-1);
			    }
#endif
			  if (j==i)
			    continue;
			  for (kk=0; kk < 3; kk++)
			    shift[kk] = L[kk]*rint((bpos[kk][i][ns]-bpos[kk][j][ns2])/L[kk]); 

			  for (kk=0; kk < 3; kk++)
			    {
			      rspA[kk] = bpos[kk][i][ns];
			      rspB[kk] = bpos[kk][j][ns2] + shift[kk];
			    }

			  distSq = Sqr(rspA[0]-rspB[0])+Sqr(rspA[1]-rspB[1])+Sqr(rspA[2]-rspB[2]);
			  sigmabs = 0.5*(typesArr[typeOfPart[i]].spots[ns].sigma + typesArr[typeOfPart[j]].spots[ns2].sigma);
#if 0
			  if ((i==2328 && j==2330)||(i==2330 && j==2328))
			    {
			      printf("[%d-%d] distSq=%f sigmaSq=%f ns=%d ns2=%d\n", i, j, distSq, Sqr(sigmabs), ns, ns2);
			    }
#endif	
			  if (distSq <= Sqr(sigmabs))
			    {
#if MC_OPT_CHECK_BS
			      checkBS[i][numcheckBS[i]]=j;
			      numcheckBS[i]++;
#else
			      checkBS[j] = 1; 
#endif
			    }
			}
		    } 
		}
	    }
	}
    }
}
#endif
#ifdef MC_NEW_NNL_CHECK
int mc_check_NNL_escape(int ip)
{
  int ns, kk;
  /* ux[0,1] sono i versori associati alle due facce del parallelepipedo perpendicolare all'asse x
   * stesso discorso per uy e uz 
   * px[0,1] sono due punti che giacciono sui suddetti piani rispettivamente */
  double dpos[3], ux[2][3], uy[2][3], uz[2][3], px[2][3], py[2][3], pz[2][3], dist, sigmah;

  for (kk=0; kk < 3; kk++)
    {
      ux[0][kk] = nebrTab[ip].R[0][kk];
      ux[1][kk] = -nebrTab[ip].R[0][kk];
      uy[0][kk] = nebrTab[ip].R[1][kk];
      uy[1][kk] = -nebrTab[ip].R[1][kk];
      uz[0][kk] = nebrTab[ip].R[2][kk];
      uz[1][kk] = -nebrTab[ip].R[2][kk];
      px[0][kk] = nebrTab[ip].r[kk] + ux[0][kk]*nebrTab[ip].axa;
      px[1][kk] = nebrTab[ip].r[kk] + ux[1][kk]*nebrTab[ip].axa;
      py[0][kk] = nebrTab[ip].r[kk] + uy[0][kk]*nebrTab[ip].axb;
      py[1][kk] = nebrTab[ip].r[kk] + uy[1][kk]*nebrTab[ip].axb;
      pz[0][kk] = nebrTab[ip].r[kk] + uz[0][kk]*nebrTab[ip].axc;
      pz[1][kk] = nebrTab[ip].r[kk] + uz[1][kk]*nebrTab[ip].axc;
    }	
  //printf("px[0]=%f %f %f nebrTab=%f %f %f sax=%f %f %f\n", px[0][0], px[0][1], px[0][2], nebrTab[ip].r[0],
  //	 nebrTab[ip].r[1], nebrTab[ip].r[2], nebrTab[ip].axa, nebrTab[ip].axb, nebrTab[ip].axc);
  for (ns=typesArr[typeOfPart[ip]].nspots; ns < typesArr[typeOfPart[ip]].nspotsBS + typesArr[typeOfPart[ip]].nspots; ns++)
    {
      //printf("ns=%d r=%f %f %f\n", ns, bpos[0][ip][ns],bpos[1][ip][ns],bpos[2][ip][ns]);
      sigmah = 0.5*typesArr[typeOfPart[ip]].spots[ns].sigma;
      /* check all 6 faces of the parallelepiped if distance
       * between any of these faces and the bounding spheres is less than 0 then
       * rebuild NNL! */
      for (kk=0; kk < 3; kk++)
	dpos[kk] = bpos[kk][ip][ns] - px[0][kk];
      dist=-scalProd(dpos, ux[0])-sigmah;
      if (dist < 0.0)
	return 1;
      for (kk=0; kk < 3; kk++)
	dpos[kk] = bpos[kk][ip][ns] - px[1][kk];
      dist=-scalProd(dpos, ux[1])-sigmah;
      if (dist < 0.0)
	return 1;
      for (kk=0; kk < 3; kk++)
	dpos[kk] = bpos[kk][ip][ns] - py[0][kk];
      dist=-scalProd(dpos, uy[0])-sigmah;
      if (dist < 0.0)
	return 1;
      for (kk=0; kk < 3; kk++)
	dpos[kk] = bpos[kk][ip][ns] - py[1][kk];
      dist=-scalProd(dpos, uy[1])-sigmah;
      if (dist < 0.0)
	return 1;
      for (kk=0; kk < 3; kk++)
	dpos[kk] = bpos[kk][ip][ns] - pz[0][kk];
      dist=-scalProd(dpos, uz[0])-sigmah;
      if (dist < 0.0)
	return 1;
      for (kk=0; kk < 3; kk++)
	dpos[kk] = bpos[kk][ip][ns] - pz[1][kk];
      dist=-scalProd(dpos, uz[1])-sigmah;
      if (dist < 0.0)
	return 1;
    }
  return 0;
}
#endif
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
#ifdef MC_BOUNDING_SPHERES
  find_part_with_BS(ip);
#endif
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
#ifdef MC_BOUNDING_SPHERES
#if 1
		      if (
#ifdef MC_OPT_CHECK_BS
			  !func_checkBS(na, n)
#else
			  !checkBS[n]
#endif
			 )
			{
#if 0
			  if (check_overlap(na, n, shift, err)<0.0)
			    {
			      printf("we have a big problem!\n");
			      store_bump(na, n);
			      exit(-1);
			    }			    
#endif
			  continue;
			}
#endif
#endif
#if defined(EDHE_FLEX) && 0
		      if (!may_interact_all(na, n))
			continue;
#endif

#if defined(MC_CLUSTER_MOVE) || defined(MC_CLUSTER_NPT) || defined(MC_NPT_XYZ)
		      /* nel caso di cluster move particelle appartenenti allo stesso cluster
			 non possono overlapparsi dopo la cluster move */
		      if (clsNPT==1 && color[na] == color[n])
			{
			  continue;
			}
#endif
		      if (check_overlap(na, n, shift, err)<0.0)
			{
#ifdef DEBUG_HCMC
			  printf("checking i=%d j=%d: ", na, n);
			  printf("overlap!\n");
			  if (dostorebump)
			    store_bump(na, n);
#endif
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
#ifdef MC_BOUNDING_SPHERES
void remove_from_current_cell_BS(int i)
{
  int n;
  n = (inCellBS[2][i] * cellsyBS + inCellBS[1][i])*cellsxBS + inCellBS[0][i]
    + totspotsBS;
  
  while (cellListBS[n] != i) 
    n = cellListBS[n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellListBS[n] = cellListBS[i];
}
void insert_in_new_cell_BS(int i)
{
  int n;
  n = (inCellBS[2][i] * cellsyBS + inCellBS[1][i])*cellsxBS + 
    inCellBS[0][i] + totspotsBS;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellListBS[i] = cellListBS[n];
  cellListBS[n] = i;
}


#endif
#ifdef MC_BOND_POS
void remove_from_current_cell_BP(int i)
{
  int n;
  n = (inCellBP[2][i] * cellsyBP + inCellBP[1][i])*cellsxBP + inCellBP[0][i]
    + totspots;
  
  while (cellListBP[n] != i) 
    n = cellListBP[n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellListBP[n] = cellListBP[i];
}

void insert_in_new_cell_BP(int i)
{
  int n;
  n = (inCellBP[2][i] * cellsyBP + inCellBP[1][i])*cellsxBP + 
    inCellBP[0][i] + totspots;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellListBP[i] = cellListBP[n];
  cellListBP[n] = i;
}

#endif
#if defined(MC_BOND_POS) || defined(MC_BOUNDING_SPHERES)
void remove_from_current_cell_all_spots(int n)
{
  int ns, k1;
  for (k1 = 0; k1 < typesArr[typeOfPart[n]].nspots; k1++)
    {
      ns = sp2n_map[n][k1];
      remove_from_current_cell_BP(ns);
    } 

#ifdef MC_BOUNDING_SPHERES
  for (k1 = typesArr[typeOfPart[n]].nspots; k1 < typesArr[typeOfPart[n]].nspots + typesArr[typeOfPart[n]].nspotsBS; k1++)
    {
      ns = sp2n_map[n][k1] - totspots;
      remove_from_current_cell_BS(ns);
    } 
#endif
}
#endif
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
#if defined(MC_BOND_POS) || defined(MC_BOUNDING_SPHERES)
void insert_in_new_cell_all_spots(int n)
{
  int k1, ns;
  for (k1 = 0; k1 < typesArr[typeOfPart[n]].nspots; k1++)
    {
      ns = sp2n_map[n][k1];
      insert_in_new_cell_BP(ns);
    }
#ifdef MC_BOUNDING_SPHERES
  for (k1 = typesArr[typeOfPart[n]].nspots; k1 < typesArr[typeOfPart[n]].nspots + typesArr[typeOfPart[n]].nspotsBS; k1++)
    {
      ns = sp2n_map[n][k1] - totspots;
      insert_in_new_cell_BS(ns);
    } 
#endif
}
#endif
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
#ifdef MC_BOND_POS
  int k1, ns, kk;
  double bp[3];
#endif
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
#ifdef MC_BOND_POS
  for (k1 = 0; k1 < typesArr[typeOfPart[n]].nspots; k1++)
    {
      ns = sp2n_map[n][k1];
      cox=inCellBP[0][ns];
      coy=inCellBP[1][ns];
      coz=inCellBP[2][ns];
      /* consider spot in the first box */
      for (kk=0; kk < 3; kk++)
	bp[kk] = bpos[kk][n][k1] - L[kk]*rint(bpos[kk][n][k1]/L[kk]);

      cx =  (bp[0] + 0.5*L[0]) * cellsxBP / L[0];
      cy =  (bp[1] + 0.5*L[1]) * cellsyBP / L[1];
      cz =  (bp[2] + 0.5*L[2]) * cellszBP / L[2];
      if (cx!=cox || cy!=coy || cz!=coz)
	{
	  remove_from_current_cell_BP(ns);
	  inCellBP[0][ns] = cx;
	  inCellBP[1][ns] = cy;
	  inCellBP[2][ns] = cz;
	  insert_in_new_cell_BP(ns);
	}
    }
#endif
#ifdef MC_BOUNDING_SPHERES
  for (k1 = typesArr[typeOfPart[n]].nspots; k1 < typesArr[typeOfPart[n]].nspots + typesArr[typeOfPart[n]].nspotsBS; k1++)
    {
      ns = sp2n_map[n][k1] - totspots;
      cox=inCellBS[0][ns];
      coy=inCellBS[1][ns];
      coz=inCellBS[2][ns];
      /* consider spot in the first box */
      for (kk=0; kk < 3; kk++)
	bp[kk] = bpos[kk][n][k1] - L[kk]*rint(bpos[kk][n][k1]/L[kk]);

      cx =  (bp[0] + 0.5*L[0]) * cellsxBS / L[0];
      cy =  (bp[1] + 0.5*L[1]) * cellsyBS / L[1];
      cz =  (bp[2] + 0.5*L[2]) * cellszBS / L[2];
      if (cx!=cox || cy!=coy || cz!=coz)
	{
	  remove_from_current_cell_BS(ns);
	  inCellBS[0][ns] = cx;
	  inCellBS[1][ns] = cy;
	  inCellBS[2][ns] = cz;
	  insert_in_new_cell_BS(ns);
	}
    }
#endif
}
#if defined(MC_CLUSTER_NPT) || defined(MC_NPT_XYZ)
void pbccls(int ip)
{
  double L2[3], Ll[3];
  double Dx, Dy, Dz;  
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
  Dx =  L[0]*rint(rx[ip]/L[0]); 
  Dy =  L[1]*rint(ry[ip]/L[1]);
  Dz =  L[2]*rint(rz[ip]/L[2]);
  Dxpar[ip] += Dx;
  Dypar[ip] += Dy;
  Dzpar[ip] += Dz;
  rx[ip] -= Dx; 
  ry[ip] -= Dy;
  rz[ip] -= Dz;
#ifdef MC_BOND_POS
  update_spot_pos(ip, -Dx, -Dy, -Dz);
#endif
}
#endif
#if 1 
void pbc(int ip)
{
  double L2[3], Ll[3];
#ifdef MC_BOND_POS
  int k1, extraspots = 0;
  double Dx, Dy, Dz;
#endif  
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
#ifdef MC_BOND_POS
  Dx = L[0]*rint(rx[ip]/L[0]);
  Dy = L[1]*rint(ry[ip]/L[1]);
  Dz = L[2]*rint(rz[ip]/L[2]);
  rx[ip] -= Dx;
  ry[ip] -= Dy;
  rz[ip] -= Dz;
  if (OprogStatus.useNNL)
    {
      nebrTab[ip].r[0] -=  Dx; 
      nebrTab[ip].r[1] -=  Dy; 
      nebrTab[ip].r[2] -=  Dz;
    }
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
  for (k1=0; k1 < typesArr[typeOfPart[ip]].nspots + extraspots; k1++)
    {
      bpos[0][ip][k1] -= Dx;
      bpos[1][ip][k1] -= Dy;
      bpos[2][ip][k1] -= Dz;
    }
#else
  rx[ip] -= L[0]*rint(rx[ip]/L[0]); 
  ry[ip] -= L[1]*rint(ry[ip]/L[1]);
  rz[ip] -= L[2]*rint(rz[ip]/L[2]);
  if (OprogStatus.useNNL)
    {
      nebrTab[ip].r[0] -=  L[0]*rint(rx[ip]/L[0]); 
      nebrTab[ip].r[1] -=  L[1]*rint(rx[ip]/L[1]); 
      nebrTab[ip].r[2] -=  L[2]*rint(rx[ip]/L[2]);
    }
#endif
}
#else
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
#endif
void update_bonds_MC(int ip);
extern void find_bonds_flex_all(void);
extern double calcpotene(void);
extern void rebuildLinkedList(void);
int cxini=-1, cyini=-1, czini=-1;
#ifdef MC_BOND_POS
int cxiniBP=-1, cyiniBP=-1, cziniBP=-1;
#endif
#ifdef MC_BOUNDING_SPHERES
int cxiniBS=-1, cyiniBS=-1, cziniBS=-1;
#endif
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
  if (cellsx*cellsy*cellsz > cxini*cyini*czini)
    {
      cellList = realloc(cellList, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
    } 
#ifdef MC_BOND_POS
#ifdef MD_LXYZ
  cellsxBP = L[0] / Oparams.rcutBP;
  cellsyBP = L[1] / Oparams.rcutBP;
  cellszBP = L[2] / Oparams.rcutBP;
#else
  cellsxBP = L / Oparams.rcutBP;
  cellsyBP = L / Oparams.rcutBP;
  cellszBP = L / Oparams.rcutBP;
#endif
  if (cellsxBP*cellsyBP*cellszBP > cxiniBP*cyiniBP*cziniBP)
    {
      cellListBP = realloc(cellListBP, sizeof(int)*(cellsxBP*cellsyBP*cellszBP+totspots));
    } 
#ifdef MC_BOUNDING_SPHERES
#ifdef MD_LXYZ
  cellsxBS = L[0] / Oparams.rcutBS;
  cellsyBS = L[1] / Oparams.rcutBS;
  cellszBS = L[2] / Oparams.rcutBS;
#else
  cellsxBS = L / Oparams.rcutBS;
  cellsyBS = L / Oparams.rcutBS;
  cellszBS = L / Oparams.rcutBS;
#endif
  if (cellsxBS*cellsyBS*cellszBS > cxiniBS*cyiniBS*cziniBS)
    {
      cellListBS = realloc(cellListBS, sizeof(int)*(cellsxBS*cellsyBS*cellszBS+totspotsBS));
    } 
#endif
#endif
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
#if defined(MC_BOND_POS) || defined(MC_BOUNDING_SPHERES)
      remove_from_current_cell_all_spots(ip);
#endif
#endif
      if (OprogStatus.useNNL)
	{
	  /* rimuove ip da tutte le NNL del sistema */
	  remove_from_nnl_MC(ip);
	}
      Oparams.parnum--;
#ifdef MC_BOND_POS
      totspots -= typesArr[0].nspots;
#ifdef MC_BOUNDING_SPHERES
      totspotsBS -= typesArr[0].nspotsBS;
#endif
#endif
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
#if defined(MC_BOND_POS) || defined(MC_BOUNDING_SPHERES)
  remove_from_current_cell_all_spots(ip);
  remove_from_current_cell_all_spots(lp);
#endif
#endif
  for (k=0; k < 3; k++)
    inCell[k][ip] = inCell[k][lp];
  typeNP[0]--;
  Oparams.parnum--;
#ifdef MC_BOND_POS
  totspots -= typesArr[0].nspots;
#ifdef MC_BOUNDING_SPHERES
  totspotsBS -= typesArr[0].nspotsBS;
#endif
#endif
  update_bonds_GC(lp, ip);
#ifdef MCGC_OPTLLREBUILD
  /* N.B. essendo cambiato Oparams.parnum le linked list vanno aggiustate (poichè
     le linked lists head stanno da Oparams.parnum compreso in poi quindi vanno shiftate, 
     in questo caso vanno spostate di uno a sinistra) */
  adjLinkedListRemove();
  insert_in_new_cell(ip);
#if defined(MC_BOND_POS) || defined(MC_BOUNDING_SPHERES)
  insert_in_new_cell_all_spots(ip);
#endif
#else
  rebuildLinkedList(); 
#endif
}
#ifdef MD_DYNAMIC_OPROG
int dyn_realloc_oprog(int np)
{
#ifdef MC_SUS
  double *sushisto;
#endif
  double **DRold;
  double *rxCMiold, *ryCMiold, *rzCMiold;
  int i, aa;  
  void *last_ptr;
#ifdef MD_CALC_DPP
  OprogStatus.len = sizeof(double)*22*np;
#else
  OprogStatus.len = sizeof(double)*10*np;
#endif
   
  DRold = malloc(sizeof(double*)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    {
      DRold[i] = malloc(sizeof(double)*3);
    }
  rxCMiold = malloc(sizeof(double)*Oparams.parnum);
  ryCMiold = malloc(sizeof(double)*Oparams.parnum);
  rzCMiold = malloc(sizeof(double)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    {
      rxCMiold[i] = OprogStatus.rxCMi[i];
      ryCMiold[i] = OprogStatus.ryCMi[i];
      rzCMiold[i] = OprogStatus.rzCMi[i];
      for (aa=0; aa < 3; aa++)
	DRold[i][aa] = OprogStatus.DR[i][aa];
    } 
#ifdef MC_SUS
  if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
    {
      sushisto = malloc(sizeof(double)*(OprogStatus.susnmax-OprogStatus.susnmin+1));
      /* salva l'istogramma prima della realloc */
      for (i=0; i < OprogStatus.susnmax-OprogStatus.susnmin+1; i++)
	sushisto[i] = OprogStatus.sushisto[i];
    }
  if (OprogStatus.susnmax > 0 && OprogStatus.susnmin >= 0) 
    OprogStatus.len += sizeof(double)*(OprogStatus.susnmax-OprogStatus.susnmin+1);
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
#ifdef MC_SIMUL
  for (i=0; i < np; i++)
    {
      if (i < Oparams.parnum)
	{
	  OprogStatus.DR[i][0] = DRold[i][0];
	  OprogStatus.DR[i][1] = DRold[i][1];
	  OprogStatus.DR[i][2] = DRold[i][2];
	  OprogStatus.rxCMi[i] = rxCMiold[i];
	  OprogStatus.ryCMi[i] = ryCMiold[i];
	  OprogStatus.rzCMi[i] = rzCMiold[i];
	}
      else
	{
	  OprogStatus.DR[i][0] = 0;
	  OprogStatus.DR[i][1] = 0;
	  OprogStatus.DR[i][2] = 0;
	  OprogStatus.rxCMi[i] = 0;
	  OprogStatus.ryCMi[i] = 0;
	  OprogStatus.rzCMi[i] = 0;
	}
#ifdef MD_CALC_DPP
      OprogStatus.sumdx[i] = 0;
      OprogStatus.sumdy[i] = 0;
      OprogStatus.sumdz[i] = 0;
      OprogStatus.lastu1x[i] = 0;
      OprogStatus.lastu1y[i] = 0;
      OprogStatus.lastu1z[i] = 0;
      OprogStatus.lastu2x[i] = 0;
      OprogStatus.lastu2y[i] = 0;
      OprogStatus.lastu2z[i] = 0;
      OprogStatus.lastu3x[i] = 0;
      OprogStatus.lastu3y[i] = 0;
      OprogStatus.lastu3z[i] = 0;
#endif
      OprogStatus.lastcolltime[i] = 0;
      OprogStatus.sumox[i] = 0;
      OprogStatus.sumoy[i] = 0;
      OprogStatus.sumoz[i] = 0;
    } 
  free(rxCMiold);
  free(ryCMiold);
  free(rzCMiold);
  for (i=0; i < Oparams.parnum; i++)
    {
      free(DRold[i]);
    };
  free(DRold);
#endif
#ifdef MC_SUS
  /* reinizializza tutto fregandosene per ora */
 if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
   {
     if (np==0)
       OprogStatus.sushisto = OprogStatus.ptr;
     else
       OprogStatus.sushisto = OprogStatus.DR[np-1]+3;
   }
 /* ripristina l'istogramma */
 if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0)
   {
     for (i=0; i < OprogStatus.susnmax-OprogStatus.susnmin+1; i++)
       OprogStatus.sushisto[i] = sushisto[i];
     free(sushisto);
   }
#endif
  OprogStatus.set_dyn_ascii();
  return OprogStatus.len;
}
#endif
extern double **matrix(int n, int m);
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);
extern double *max_step_MC;
void check_alloc_GC(void)
{
  int size, allocnpGCold, k, i;
  double *rt[3];
#ifdef MC_BOND_POS
  int extraspots, kk, k1;
#endif
  if (Oparams.parnum+1 > allocnpGC)
    {
      allocnpGCold=allocnpGC;
      allocnpGC = (int) (1.2*((double)allocnpGC));
      if (allocnpGC==allocnpGCold)
	allocnpGC *= 2;
      if (allocnpGC == 0)
	allocnpGC = 1;
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
      /* un array contiguo facendo la free di rx dealloca tutte la variabili
	 in  ALLOC_LIST */
      free(rx);
      AllocCoord(allocnpGC*sizeof(double), ALLOC_LIST, NULL); 
      for (i=0; i < allocnpGC; i++)
	{
	  if (i < Oparams.parnum)
	    {
	      rx[i] = rt[0][i];
	      ry[i] = rt[1][i];
	      rz[i] = rt[2][i];
	    }
	  else
	    {
	      rx[i] = 0;
	      ry[i] = 0;
	      rz[i] = 0;
	    }
#ifndef MD_ASYM_ITENS
	  lastcol[i] = 0;
#endif
	  uxx[i] = 0;
	  uxy[i] = 0;
	  uxz[i] = 0;
	  uyx[i] = 0;
	  uyy[i] = 0;
	  uyz[i] = 0;
	  uzx[i] = 0;
	  uzy[i] = 0;
	  uzz[i] = 0;
	  wx[i] = 0;
	  wy[i] = 0;
	  wz[i] = 0;
	  Mx[i] = 0;
	  My[i] = 0;
	  Mz[i] = 0;
	  vx[i] = 0;
	  vy[i] = 0;
	  vz[i] = 0;
#ifdef MD_POLYDISP
	  axaP[i] = 0;
	  axbP[i] = 0;
	  axcP[i] = 0;
#endif
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
#ifdef MC_BOND_POS
      /* per continuiamo ad assumere un solo tipo di particelle nel GC */
      bpos[0] = malloc(sizeof(double*)*allocnpGC);
      bpos[1] = malloc(sizeof(double*)*allocnpGC);
      bpos[2] = malloc(sizeof(double*)*allocnpGC); 
      bposold[0] = malloc(sizeof(double*)*allocnpGC);
      bposold[1] = malloc(sizeof(double*)*allocnpGC);
      bposold[2] = malloc(sizeof(double*)*allocnpGC); 
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[0].nspotsBS;
#else
      extraspots = 0;
#endif
      bpos[0] = realloc(bpos[0],sizeof(double*)*allocnpGC*(typesArr[0].nspots+extraspots));
      bpos[1] = realloc(bpos[1],sizeof(double*)*allocnpGC*(typesArr[0].nspots+extraspots));
      bpos[2] = realloc(bpos[2],sizeof(double*)*allocnpGC*(typesArr[0].nspots+extraspots)); 
      bposold[0] = realloc(bposold[0],sizeof(double*)*allocnpGC*(typesArr[0].nspots+extraspots));
      bposold[1] = realloc(bposold[1],sizeof(double*)*allocnpGC*(typesArr[0].nspots+extraspots));
      bposold[2] = realloc(bposold[2],sizeof(double*)*allocnpGC*(typesArr[0].nspots+extraspots)); 
      for (i=1; i < allocnpGC; i++)
	{
	  for (kk=0; kk < 3; kk++)
	    {
#ifdef MC_BOUNDING_SPHERES
	      bpos[kk][i] = bpos[kk][i-1] + typesArr[typeOfPart[i-1]].nspots + typesArr[typeOfPart[i-1]].nspotsBS;
	      bposold[kk][i] = bposold[kk][i-1] + typesArr[typeOfPart[i-1]].nspots + typesArr[typeOfPart[i-1]].nspotsBS;
#else
    	      bpos[kk][i] = bpos[kk][i-1] + typesArr[typeOfPart[i-1]].nspots;
    	      bposold[kk][i] = bposold[kk][i-1] + typesArr[typeOfPart[i-1]].nspots;
#endif
	    }
	}
      cellListBP = realloc(cellListBP,sizeof(int)*(cellsxBP*cellsyBP*cellszBP+allocnpGC*typesArr[0].nspots));
      inCellBP[0] = realloc(inCellBP[0],sizeof(int)*allocnpGC*typesArr[0].nspots);
      inCellBP[1] = realloc(inCellBP[1],sizeof(int)*allocnpGC*typesArr[0].nspots);
      inCellBP[2] = realloc(inCellBP[2],sizeof(int)*allocnpGC*typesArr[0].nspots);
#ifdef MC_BOUNDING_SPHERES
      cellListBS = realloc(cellListBS,sizeof(int)*(cellsxBS*cellsyBS*cellszBS+allocnpGC*extraspots));
      inCellBS[0] = realloc(inCellBS[0],sizeof(int)*allocnpGC*extraspots);
      inCellBS[1] = realloc(inCellBS[1],sizeof(int)*allocnpGC*extraspots);
      inCellBS[2] = realloc(inCellBS[2],sizeof(int)*allocnpGC*extraspots);
      checkBS = realloc(checkBS,sizeof(int)*allocnpGC);
#endif
      n2sp_map = realloc(n2sp_map,sizeof(struct n2sp_struct)*allocnpGC*(typesArr[0].nspots+extraspots));
      sp2n_map = realloc(sp2n_map,sizeof(int*)*allocnpGC);
      sp2n_map[0] = realloc(sp2n_map[0],sizeof(int)*allocnpGC*(typesArr[0].nspots+extraspots));
      for (i=1; i < allocnpGC; i++)
	{
#ifdef MC_BOUNDING_SPHERES
	  sp2n_map[i] = sp2n_map[i-1] + typesArr[typeOfPart[i-1]].nspots + typesArr[typeOfPart[i-1]].nspotsBS;
#else  
	  sp2n_map[i] = sp2n_map[i-1] + typesArr[typeOfPart[i-1]].nspots;
#endif
	}
#endif 
#ifdef MD_DYNAMIC_OPROG
      dyn_realloc_oprog(allocnpGC);
#endif
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
#ifdef MC_BOND_POS
  int ns, k1;
#endif
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
#ifdef MC_BOND_POS
  for (k1 = 0; k1 < typesArr[typeOfPart[np]].nspots; k1++)
    {
      ns = sp2n_map[np][k1];
      inCellBP[0][ns] =  (bpos[0][np][k1] + L2[0]) * cellsxBP / L[0];
      inCellBP[1][ns] =  (bpos[1][np][k1] + L2[1]) * cellsyBP / L[1];
      inCellBP[2][ns] =  (bpos[2][np][k1] + L2[2]) * cellszBP / L[2];
    }
#endif
#ifdef MC_BOUNDING_SPHERES
  for (k1 = typesArr[typeOfPart[np]].nspots; k1 < typesArr[typeOfPart[np]].nspots + typesArr[typeOfPart[np]].nspotsBP; k1++)
    {
      ns = sp2n_map[n][k1] - totspots;
      inCellBS[0][ns] =  (bpos[0][np][k1] + L2[0]) * cellsxBS / L[0];
      inCellBS[1][ns] =  (bpos[1][np][k1] + L2[1]) * cellsyBS / L[1];
      inCellBS[2][ns] =  (bpos[2][np][k1] + L2[2]) * cellszBS / L[2];
    } 
#endif
}
#endif

#ifdef MC_RESTR_MATRIX
extern double restrMatrix[3][3];
#endif
extern double eval_max_dist_for_spots(int pt);
void calc_maxax(void)
{
  double MAXAX, maxSpots;
  int i, a;
  MAXAX = 0.0;

  for (i = 0; i < Oparams.parnum; i++)
    {
      maxax[i] = 0.0;
#if defined(MD_POLYDISP) 
      if (axaP[i] > maxax[i])
	maxax[i] = axaP[i];
      if (axbP[i] > maxax[i])
	maxax[i] = axbP[i];
      if (axcP[i] > maxax[i])
	maxax[i] = axcP[i];
#elif defined(EDHE_FLEX)
#ifdef MD_SUPERELLIPSOID
#if 0
      if (typesArr[typeOfPart[i]].sax[0] > maxax[i])
	maxax[i] = typesArr[typeOfPart[i]].sax[0];
      if (typesArr[typeOfPart[i]].sax[1] > maxax[i])
	maxax[i] = typesArr[typeOfPart[i]].sax[1];
      if (typesArr[typeOfPart[i]].sax[2] > maxax[i])
	maxax[i] = typesArr[typeOfPart[i]].sax[2];
#else 
      /* per i superellissoidi il centroide non puo' essere uguale al raggio maggiore,
	 così lo sceglo pari alla metà della diagonale maggiore del parallelepipedo
	 con i lati pari al doppio dei semi-assi */
      maxax[i] = sqrt(Sqr(typesArr[typeOfPart[i]].sax[0])+Sqr(typesArr[typeOfPart[i]].sax[1])+
	Sqr(typesArr[typeOfPart[i]].sax[2]));
#endif
#else
      maxax[i] = sqrt(Sqr(typesArr[typeOfPart[i]].sax[0])+Sqr(typesArr[typeOfPart[i]].sax[1])+
	Sqr(typesArr[typeOfPart[i]].sax[2]));
#endif
#ifdef MC_SIMUL
      /* se è una sfera setta maxax pari al raggio della stessa */
      if (typesArr[typeOfPart[i]].sax[0] == typesArr[typeOfPart[i]].sax[1] && 
	    typesArr[typeOfPart[i]].sax[1] == typesArr[typeOfPart[i]].sax[2]) 
	{
#ifdef MD_SUPERELLIPSOID
    	  if (!is_superellipse(i))
    	    {
	      maxax[i] = typesArr[typeOfPart[i]].sax[0];
    	    }
#else
	  maxax[i] = typesArr[typeOfPart[i]].sax[0];
#endif
	}
#endif
      maxSpots = eval_max_dist_for_spots(typeOfPart[i]);
      if (maxSpots > maxax[i])
	maxax[i] = maxSpots;
#ifdef MC_SIMUL
      maxax[i] *= 1.0001;
#endif
      //printf("maxax[%d]:%f maxSpots:%f\n", i, 2.0*maxax[i], 2.0*maxSpots);
#else
      a=(i<Oparams.parnumA)?0:1;
      if (Oparams.a[a] > maxax[i])
	maxax[i] = Oparams.a[a];
      if (Oparams.b[a] > maxax[i])
	maxax[i] = Oparams.b[a];
      if (Oparams.c[a] > maxax[i])
	maxax[i] = Oparams.c[a];
#endif
      //printf("distSPA=%.15G distSPB=%.15G\n", distSPA, distSPB);
      maxax[i] *= 2.0;

      if (maxax[i] > MAXAX)
	MAXAX = maxax[i];
      //printf("maxax aft[%d]: %.15G\n", i, maxax[i]);
    }

}
void calc_ax(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      axa[i] = typesArr[typeOfPart[i]].sax[0];
      axb[i] = typesArr[typeOfPart[i]].sax[1];
      axc[i] = typesArr[typeOfPart[i]].sax[2];
    }
}
extern void find_spheres_NNL(void);
int insert_particle_GC(void)
{
#ifdef MC_BOND_POS
  int extraspots = 0;
#endif
  int np, k1, k2;
  double ox, oy, oz, Rl[3][3];
  np = Oparams.parnum;
  /* controlla se c'è abbastanza allocazione se no espandere
     tutto ciò che va espanso */
  check_alloc_GC();
  /* choose a random position to insert particle */
  if (OprogStatus.restrmove == 2)
    {
      if (Oparams.parnum < OprogStatus.susnmax/2)
	{
	  /* insert nematic particles in half box z > 0 */
#ifdef MD_LXYZ
    	  rx[np] = L[0]*(ranf()-0.5);
	  ry[np] = L[1]*(ranf()-0.5);
	  rz[np] = L[2]*(ranf()*0.5);
#else
	  rx[np] = L*(ranf()-0.5);
	  ry[np] = L*(ranf()-0.5);
	  rz[np] = L*(ranf()*0.5);
#endif
	}
      else
	{
	  /* insert isotropic particles in half box z < 0 */
#ifdef MD_LXYZ
	  rx[np] = L[0]*(ranf()-0.5);
	  ry[np] = L[1]*(ranf()-0.5);
	  rz[np] = -L[2]*(ranf()*0.5);
#else
	  rx[np] = L*(ranf()-0.5);
	  ry[np] = L*(ranf()-0.5);
	  rz[np] = -L*(ranf()*0.5);
#endif
	}
    }
  else 
    {
#ifdef MD_LXYZ
      rx[np] = L[0]*(ranf()-0.5);
      ry[np] = L[1]*(ranf()-0.5);
      rz[np] = L[2]*(ranf()-0.5);
#else
      rx[np] = L*(ranf()-0.5);
      ry[np] = L*(ranf()-0.5);
      rz[np] = L*(ranf()-0.5);
#endif
    }
  if (OprogStatus.restrmove==0)
    {
      orient(&ox,&oy,&oz);
      versor_to_R(ox, oy, oz, Rl);
    }
  else
    {
      orient(&ox,&oy,&oz);
      versor_to_R(ox, oy, oz, Rl);
      if (OprogStatus.restrmove == 2)
	{
	  if (Oparams.parnum < OprogStatus.susnmax/2)
	    {
	      Rl[0][0]=Rl[1][1]=Rl[2][2]=1.0;
	      Rl[0][1]=Rl[0][2]=Rl[1][0]=Rl[1][2]=Rl[2][0]=Rl[2][1]=0.0;
	    }
	}
#ifdef MC_RESTR_MATRIX
      else if (OprogStatus.restrmove == 3)
	{
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		Rl[k1][k2] = restrMatrix[k1][k2];
		//printf("Rl[%d][%d]=%f\n", k1, k2, restrMatrix[k1][k2]);
	      }
	}
#endif
      else
	{
	  Rl[0][0]=Rl[1][1]=Rl[2][2]=1.0;
	  Rl[0][1]=Rl[0][2]=Rl[1][0]=Rl[1][2]=Rl[2][0]=Rl[2][1]=0.0;
	}
    }
  /* per ora assumiamo un solo tipo di particelle nel GC */
  numbonds[np] = 0;
  typeOfPart[np]=0;
  vx[np]=vy[np]=vz[np]=wx[np]=wy[np]=wz[np]=Mx[np]=My[np]=Mz[np]=0.0;

  if (Oparams.parnum==1)
    {
      /* se parnum=1 vuol dire che prima c'erano 0 particelle quindi vanno calcolati i valori in 0 per
	 axa, axb, axc e maxax */
      find_spheres_NNL();
      calc_maxax();
      calc_ax();
    }
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
#ifdef MC_BOND_POS
  totspots += typesArr[0].nspots;
#ifdef MC_BOUNDING_SPHERES
  totspotsBS += typesArr[0].nspotsBS;
#endif
#endif
#ifdef MC_BOND_POS
  rA[0] = rx[np];
  rA[1] = ry[np];
  rA[2] = rz[np];
  BuildAtomPos(np, rA, R[np], ratA);
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[np]].nspotsBS;
#endif
  for (k1 = 0; k1 < typesArr[typeOfPart[np]].nspots + extraspots; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	bpos[k2][np][k1] = ratA[k1+1][k2];
    }
#endif
 
#ifdef MCGC_OPTLLREBUILD
  /* N.B. essendo cambiato Oparams.parnum le linked list vanno aggiustate (poichè
     le linked lists head stanno da Oparams.parnum compreso in poi  vanno shiftate, 
     in questo caso vanno spostate di uno a destra) */
  adjLinkedListInsert();
  assign_cell_GC(np);
  insert_in_new_cell(np);
#if defined(MC_BOND_POS) || defined(MC_BOUNDING_SPHERES)
  insert_in_new_cell_all_spots(np);
#endif
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
#ifdef MC_AMYLOID_FIBRILS
double calc_el_ij(int i, int j)
{
  double norm, theta, sp, uxi[3], uyi[3], uzi[3], uxj[3], uyj[3], uzj[3];
  int kk;
  for (kk=0; kk < 3; kk++)
    {
      /* permanent patches are along x-axis */
      uxi[kk] = R[i][0][kk];
      uyi[kk] = R[i][1][kk];
      uxj[kk] = R[j][0][kk];
      uyj[kk] = R[j][1][kk];
    }
#if 0
  sp=scalProd(uyj,uxi);
  for (kk=0; kk < 3; kk++)
    uyj[kk] = uyj[kk] - sp*uxi[kk];
  norm = calc_norm(uyj);
  for (kk=0; kk < 3; kk++)
    uyj[kk] /= norm;
#endif  
  sp = scalProd(uyi,uyj);
  /* theta in gradi */
  theta = fabs(90.0*acos(sp)/acos(0.0));
  return OprogStatus.tors_k*Sqr(theta - OprogStatus.tors_theta0);
}
double calc_elastic_torsional_energy(int ip)
{
  int nextip, previp, nn, ncov;
  double elene; 
#ifdef MD_LL_BONDS
  int i, nb;
  long long int aa, bb, ii, jj, jj2;
  long long int covmonomers[2];
#else
  int i, nb, ii, jj, aa, bb, jj2;
  int covmonomers[2];
#endif
  ncov = 0;
  elene = 0;
  /* fine previsous (if any) and next (if any) particles which
     are covalently bonded */
  for (nn = 0; nn < numbonds[ip]; nn++)
    {
      jj = bonds[ip][nn] / (NANA);
      jj2 = bonds[ip][nn] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (aa == 1 || aa == 2)
	{
	  //printf("ip=%d numbonds[ip]=%d jj=%d ncov=%d\n", ip, numbonds[ip], jj, ncov);
	  covmonomers[ncov] = jj;
	  ncov++;
	}
    }
  //printf("covmonomers: %lld %lld\n", covmonomers[0], covmonomers[1]);
  for (nn = 0; nn < ncov; nn++)
    {
      elene += calc_el_ij(ip, ((int)covmonomers[nn]));
    }
  //printf("elene(i=%d)=%.15G\n", ip, elene);
  return elene;
}
#endif
#ifdef MC_HYDROPHOBIC_INT
extern double **eneij;
#endif
#ifdef ALIGN_POT
double calc_align_pot(int ip)
{
  return OprogStatus.Ual*(1.0-Sqr(R[ip][0][OprogStatus.alignaxis]));
}
#endif
#ifdef GAPDNA_BENDING_ENERGY
double calc_bending_energy(int i)
{
  double upot, norm, theta, sp, uxi[3], uxj[3];
  int j, kk;
  /* i cilindri che rappresentano le due stick sono numerati 
   * in maniera contigua, da cui la condizione di seguito */
  if ((i%OprogStatus.polylen)==0)
    j = i+1;
  else
    j = i-1;
  for (kk=0; kk < 3; kk++)
    {
      /* permanent patches are along x-axis */
      uxi[kk] = R[i][0][kk];
      uxj[kk] = R[j][0][kk];
    }
  sp = scalProd(uxi,uxj);
  upot = OprogStatus.kbend*Sqr((1.0+sp)/2.0);
  //printf("sp=%f upot=%f\n", sp, upot);
  /* theta in gradi */
  return upot;
}
#endif
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
#ifdef MC_HYDROPHOBIC_INT
	      Epot -= eneij[ip][jj]*intersArr[kk2].bheight;
#else
	      //if (aa==2 && bb==2)
		//printf("height=%f\n", intersArr[kk2].bheight);
	      Epot -= intersArr[kk2].bheight;
#endif
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
#ifdef GAPDNA_BENDING_ENERGY
  Epot += calc_bending_energy(ip);
#endif
#ifdef MC_AMYLOID_FIBRILS
  /* add torsional elastic energy here! */
  Epot += calc_elastic_torsional_energy(ip);
#endif
#ifdef ALIGN_POT
  Epot += calc_align_pot(ip);
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
      /* zetaMC= exp(mu/kBT)/lambda^3 (vedi Frenkel pag. 132) */
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
      if (OprogStatus.susnmin==-1 && OprogStatus.susnmax > 0 && Oparams.parnum >= OprogStatus.susnmax)
	return;
#if defined(MC_SUS) && defined(MC_SWELL)	  
      if (OprogStatus.susnmin >= 0 && OprogStatus.susnmax > 0 &&  Oparams.parnum >= OprogStatus.susnmax)
	{
	  OprogStatus.sushisto[Oparams.parnum-OprogStatus.susnmin]++;
	  return;
	}
#endif

      np=insert_particle_GC();
      //printf("FINE-2 MCEXC\n");//// 

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
#ifdef MC_BOND_POS
    	  totspots -= typesArr[0].nspots;
#ifdef MC_BOUNDING_SPHERES
	  totspotsBS -= typesArr[0].nspotsBS;
#endif
#endif
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
#ifdef MC_BOND_POS
    	  totspots -= typesArr[0].nspots;
#ifdef MC_BOUNDING_SPHERES
	  totspotsBS -= typesArr[0].nspotsBS;
#endif
#endif
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
  printf("Current code for storing/restoring linked lists in volume change move is broken!\n");
  printf("Please compile with -UMC_STORELL flag\n");
  exit(-1);
#if defined(MC_CLUSTER_NPT) || defined(MC_NPT_XYZ)
  cellsxMC = cellsx;
  cellsyMC = cellsy;
  cellszMC = cellsz;
#endif 

  for (k=0; k < Oparams.parnum + cellsx*cellsy*cellsz; k++)
    {
      cellListMC[k] = cellList[k];
    }
}
void restore_ll_mc(void)
{
  int k;
#if defined(MC_CLUSTER_NPT) || defined(MC_NPT_XYZ)
  cellsx = cellsxMC;
  cellsy = cellsyMC;
  cellsz = cellszMC; 
#endif
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
	    {
  	      bondsMC[i][k] = bonds[i][k]; 
	    }
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
#if 1
int check_bonds_mc(char *txt)
{
  int i;
  double enoo, eno;
  enoo=calcpotene();

  for (i=0; i < Oparams.parnum; i++)
    numbonds[i] = 0;

  find_bonds_flex_all();
  eno=calcpotene();
  if (eno!=enoo)
    {
      if (txt) 
	printf("%s\n", txt);
      printf("[step #%d] boh eno=%f enoo=%f\n", Oparams.curStep, eno, enoo);
      return 1;
      //exit(-1);
    }
  return 0;
}
#endif

#ifdef MC_BOND_POS
void update_spot_pos_dir(int i, int dir, double del)
{
  int k2;
  int extraspots = 0;

#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[i]].nspotsBS;
#endif
  for (k2=0; k2 < typesArr[typeOfPart[i]].nspots + extraspots; k2++)
    {
      bpos[dir][i][k2] += del;
    }
}
void update_spot_pos(int i, double dx, double dy, double dz)
{
  int k2;
  int extraspots = 0;

#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[i]].nspotsBS;
#endif
  for (k2=0; k2 < typesArr[typeOfPart[i]].nspots + extraspots; k2++)
    {
      bpos[0][i][k2] += dx;
      bpos[1][i][k2] += dy;
      bpos[2][i][k2] += dz;
    }
}
#endif
#ifdef MC_STORE_ALL_COORDS
void store_all_coords(void)
{
  int i, k1, k2;
#ifdef MC_BOND_POS
  int extraspots = 0;
#endif
  memcpy(vx, rx, sizeof(double)*Oparams.parnum);
  memcpy(vy, ry, sizeof(double)*Oparams.parnum);
  memcpy(vz, rz, sizeof(double)*Oparams.parnum);
#ifdef MC_BOND_POS
  for (i=0; i < Oparams.parnum; i++)
    {
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[typeOfPart[i]].nspotsBS;
#endif
      for (k1=0; k1 < 3; k1++)
    	for (k2=0; k2 < typesArr[typeOfPart[i]].nspots + extraspots; k2++)
  	  {
  	    bposold[k1][i][k2] = bpos[k1][i][k2];
  	  } 
    }
#endif
}
void restore_all_coords(void)
{
  int k1, k2, i;
#ifdef MC_BOND_POS
  int extraspots = 0;
#endif
  memcpy(rx, vx, sizeof(double)*Oparams.parnum);
  memcpy(ry, vy, sizeof(double)*Oparams.parnum);
  memcpy(rz, vz, sizeof(double)*Oparams.parnum);
#ifdef MC_BOND_POS
  for (i=0; i < Oparams.parnum; i++)
    {
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[typeOfPart[i]].nspotsBS;
#endif
      for (k1=0; k1 < 3; k1++)
    	for (k2=0; k2 < typesArr[typeOfPart[i]].nspots + extraspots; k2++)
  	  {
  	    bpos[k1][i][k2] = bposold[k1][i][k2];
  	  } 
    }
#endif
}
#endif
#if defined(MC_NPT_XYZ) 
#define MC_CLSNPT_LOG

#ifdef MC_ALMARZA
double calc_almarza_prob(void)
{
  int i, k;
  long long int jj, jj2, aa, bb; 
  double totprod=1.0, prod=1.0, thr;

  thr = OprogStatus.almarza_thr;
  for (i=0; i < Oparams.parnum; i++)
    {
      prod=1.0;
      for (k=0; k < numbondsMC[i]; k++)
	{
	  jj = bondsMC[i][k] / (NANA);
	  jj2 = bondsMC[i][k] & (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (color[i] != color[jj] && i < jj)
	    {
	      //printf("qui 1-thr=%f prod=%f i=%d\n", 1.0-thr, prod, i);
  	      prod *= 1.0/(1.0-thr);
	    }
	}
      /* costruisce i cluster della trial configuration */
      for (k=0; k < numbonds[i]; k++)
	{
	  jj = bonds[i][k] / (NANA);
	  jj2 = bonds[i][k] & (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (color[i]!=color[jj] && i < jj)
	    prod *= (1.0-thr);

	}
#ifdef MC_OPT_ALMARZA
      if (!inside_cluster_move)
	{
	  for (k=0; k < numcheckBS_MC[i]; k++)
	    {
	      if (color[i] != color[checkBS_MC[i][k]] && i < checkBS_MC[i][k])
		prod *= 1.0/(1.0-thr);
	    }
	  for (k=0; k < numcheckBS[i]; k++)
	    {
	      if (color[i] != color[checkBS[i][k]] && i < checkBS[i][k])
		prod *= (1.0-thr);
	    }
	}
#endif
      totprod *= prod;
   }
#if 0
  if (totprod!=1.0)
    printf("totprod=%f\n", totprod);
#endif
  return totprod;
}

void assign_bonds_almarza(void)
{
  int i, k, yn, k2;
  long long int jj, jj2, jj3, aa, bb, aa3, bb3;
  double thr, xi;
  thr = OprogStatus.almarza_thr;//exp(-1.0/Oparams.T);// poi diventerà un parametro da settare in nel file .par tramite OprogStatus
  //printf("bondsYN[0][0]=%d\n", bondsYN[0][0]);
  for (i=0; i < Oparams.parnum; i++)
    {
      for (k=0; k < numbonds[i]; k++)
	{
	  bondsYN[i][k] = -1;
	}
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      for (k=0; k < numbonds[i]; k++)
	{
	  xi = ranf();
	  jj = bonds[i][k] / (NANA);
	  jj2 = bonds[i][k] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (bondsYN[i][k]==-1) 
	    {
	      if (aa == 2 && bb == 2)
	      	{ 
	  	  yn = (ranf() > thr)?0:1;
		  bondsYN[i][k] = yn;  
		  for (k2 = 0; k2 < numbonds[jj]; k2++)
	  	    {
	  	      jj2 = bonds[jj][k2] / (NANA); 
		      jj3 = bonds[jj][k2] % (NANA);
		      aa3 = jj3 / NA;
		      bb3 = jj3 % NA;
		      if (jj2 == i && aa3==2 && bb3==2)
			{
			  bondsYN[jj][k2] = yn;  
			}
		    }
	    	}
	      else if (aa == 1 && bb == 1)
		{
		  //printf("1!!! jj2=%lld aa=%lld bb==%lld i=%d k=%d\n",jj2, aa, bb, i, k); 
		  bondsYN[i][k] = 1; /* for permanent bonds always yes! */
		  for (k2 = 0; k2 < numbonds[jj]; k2++)
		    {
		      jj2 = bonds[jj][k2] / (NANA); 
		      jj3 = bonds[jj][k2] % (NANA);
		      aa3 = jj3 / NA;
		      bb3 = jj3 % NA;
	  	      if (jj2 == i && aa3==1 && bb3==1)
			bondsYN[jj][k2] = 1;  
		    }
		}
	      else
		{ printf("aa=%lld bb=%lld BOH\n", aa, bb); 
		  exit(-1);
		}
	    }
	}
#ifdef MC_OPT_ALMARZA
      for (k = 0; k < numcheckBS[i]; k++)
	{
	  if (bondsYN_BS[i][k] == -1)
	    {
	      yn = (ranf() > thr)?0:1;
	      bondsYN_BS[i][k] = yn;  
	    }
	  jj2 = checkBS[i][k];
	  for (k2=0; k2 < numcheckBS[jj2]; k2++)
	    {
	      if (checkBS[jj2][k2] == i)
		bondsYN_BS[jj2][k2] = yn;
	    }
	}
#endif
    }
#if 0
  for (i=0; i < Oparams.parnum; i++)
    {
      for (k=0; k < numbonds[i]; k++)
	{
	  if (bondsYN[i][k] == -1)
	    {
	      printf("i=%d k=%d numbonds[%d]=%d\n", i, k, i, numbonds[i]);
	      exit(-1);
	    }
	}
    }
#endif
}
#endif
void move_box_cluster_xyz(int *ierr)
{
  int i0, i, ii, k, nc, ncls=0, percolating=0, np_in_cls, np, in0, in1, iold, kk;
  double nn, distance, lnvn;
  double vo, vn, Lfact, enn, eno, arg, delv=0.0, lastrx=0, lastry=0, lastrz=0, Dx, Dy, Dz;
  int dir;
  double *rr[3], del[3];
  double omegaratio;
#ifdef MC_BOND_POS
  int k1, k2;
#endif
#ifdef MC_STORE_ALL_COORDS
#ifdef MD_LXYZ
  double Lold[3];
#else
  double Lold;
#endif
#endif
  //printf("moving box\n");
#ifdef MD_LXYZ
  vo = L[0]*L[1]*L[2];
#else
  vo = L*L*L;
#endif
  rr[0] = rx;
  rr[1] = ry;
  rr[2] = rz;
  dir = 3*ranf();
  eno = calcpotene();

#ifdef MC_CLSNPT_LOG
  lnvn = log(vo) + (ranf()-0.5)*OprogStatus.vmax;
  vn = exp(lnvn);
#else
  vn = vo + (2.0*ranf()-1.0)*OprogStatus.vmax;
#endif

  Lfact = vn/vo;

  if (OprogStatus.useNNL)
    delv = (Lfact - 1.0); /* FINISH HERE */
  //printf("Lfact=%.15G vmax=%f vn=%f vo=%f\n", Lfact, OprogStatus.vmax, vn, vo); 
  // Lfact=1;
  
#ifdef MC_ALMARZA
#ifdef MC_OPT_ALMARZA
  for (i=0; i < Oparams.parnum; i++)
    find_part_with_BS(i);
  for (i=0; i < Oparams.parnum; i++)
    {
      numcheckBS_MC[i] = numcheckBS[i];
      for (k=0; k < numcheckBS[i]; k++)
	{
	  checkBS_MC[i][k] = checkBS[i][k];
	}
    }
#endif
  assign_bonds_almarza();
#endif

  build_clusters(&ncls, &percolating, 1);
  numOfClusters = ncls;

  //printf("ncls=%d percolating=%d\n", ncls, percolating);
  /* calcola i centri di massa dei cluster */
#if 1
  if (percolating)
    {
      //printf("is percolating\n");
      volrejMC++;
      return;
    }
#endif
#ifdef MC_STORE_ALL_COORDS
  /* I use velocity arrays to store coordinates, in MC simulations they are unused! */
  store_all_coords();
#endif
  for (nc=0; nc < ncls; nc++)
    {
      for (k=0; k < 3; k++)
	clsCoM[k][nc] = 0.0;
      /* if it is percolating all particles have 2 bonds, 
	 hence we have to initialize i0 */
      i0 =  clsarr[firstofcls[nc]];
      for (np=0; np < clsdim[nc]; np++)
	{
	  i =  clsarr[firstofcls[nc]+np];
	  if ( numbonds[i]==1) 
	    {
	      i0 = i;
	      //printf("qui i0=%d\n", i0);
	      break;
	    }
	}

      np_in_cls=1;
      i =  bonds[i0][0] / (NANA);
      iold=i0;	
      clsCoM[0][nc] = rx[i0];
      clsCoM[1][nc] = ry[i0];
      clsCoM[2][nc] = rz[i0]; 
#if !defined(MC_STORE_ALL_COORDS)
      Dxpar[i0]=Dypar[i0]=Dzpar[i0]=0;
#endif
      /* if particles has no bond we have to check... */ 
      while (numbonds[i0]>0)
	{
	  /* rebuild the cluster as a whole */
	  lastrx = rx[iold];
	  lastry = ry[iold];
	  lastrz = rz[iold];
#ifdef MC_STORE_ALL_COORDS
	  del[0] = L[0]*rint((rx[i]-lastrx)/L[0]);
	  del[1] = L[1]*rint((ry[i]-lastry)/L[1]); 
	  del[2] = L[2]*rint((rz[i]-lastrz)/L[2]);
	  rx[i] -= del[0];
	  ry[i] -= del[1];
	  rz[i] -= del[2];
#ifdef MC_BOND_POS
	  update_spot_pos(i, -del[0], -del[1], -del[2]);
#endif
#else
	  Dxpar[i] = L[0]*rint((rx[i]-lastrx)/L[0]);
	  Dypar[i] = L[1]*rint((ry[i]-lastry)/L[1]);
	  Dzpar[i] = L[2]*rint((rz[i]-lastrz)/L[2]); 
	  rx[i] -= Dxpar[i];
	  ry[i] -= Dypar[i];
	  rz[i] -= Dzpar[i];
#ifdef MC_BOND_POS
	  update_spot_pos(i, -Dxpar[i], -Dypar[i], -Dzpar[i]);
#endif
#endif
	  clsCoM[0][nc] += rx[i];
	  clsCoM[1][nc] += ry[i];
	  clsCoM[2][nc] += rz[i];
	  np_in_cls++;
#if 0
	  if (numbonds[i] > 2)
	    {
	      printf("more than 2 bonds (%d) for particle %d!\n", numbonds[i], i);
	    }
#endif
	  if(numbonds[i]==1 || np_in_cls == clsdim[nc])
	    {
	      //printf("numbonds[i=%d]=%d numbonds[i0=%d]=%d", i, numbonds[i], i0, numbonds[i0]);
	      break;
	    }	   
	  /* qui assumo di avere a che fare con particelle bifunzionali */
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
#if 0
      if (np_in_cls != clsdim[nc])
	{
	  printf("we got a problem np_in_cls=%d  but clsdim[%d]=%d\n", np_in_cls, nc, clsdim[nc]);
	  printf("percolating=%d\n", percolating);
	  exit(-1);
	}
	if (np_in_cls!=clsdim[nc])
	  {
	    printf("boh...\n");
	    exit(-1);
	  }
#endif
    }
#ifdef MC_STORE_ALL_COORDS
  for (kk=0; kk < 3; kk++)
    Lold[kk] = L[kk];
#endif
  L[dir] *= Lfact;
  L2[dir] = L[dir]*0.5;

  for (nc=0; nc < ncls; nc++)
    {
      for (k=0; k < 3; k++)
	clsCoM[k][nc] /= clsdim[nc];
      Dxcls[nc] = Dycls[nc] = Dzcls[nc] = 0;
      switch (dir) {
      case 0:
	Dxcls[nc]=(Lfact-1.0)*clsCoM[0][nc];
	break;
      case 1:
	Dycls[nc]=(Lfact-1.0)*clsCoM[1][nc];
	break;
      case 2:
	Dzcls[nc]=(Lfact-1.0)*clsCoM[2][nc];
	break;
      }
    }
#if 0
  for (nc=0; nc < ncls; nc++)
    {
      for (np=0; np < clsdim[nc]; np++)
	{
	  i = clsarr[firstofcls[nc]+np];
	  rx[i] += Dxcls[color[i]];
	  ry[i] += Dycls[color[i]];
	  rz[i] += Dzcls[color[i]]; 
	  pbc(i);	
	}
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    {  
      rx[i] += Dxcls[color[i]];
      ry[i] += Dycls[color[i]];
      rz[i] += Dzcls[color[i]];
#ifdef MC_BOND_POS
      update_spot_pos(i, Dxcls[color[i]],Dycls[color[i]],Dzcls[color[i]]); 
#endif
#ifdef MC_STORE_ALL_COORDS
      pbc(i);
#else
      pbccls(i);
#endif
    }
#endif
#ifdef MC_STORELL
  store_ll_mc();
#endif
  update_numcells();
  rebuildLinkedList();
  for (i=0; i < Oparams.parnum; i++)
    {
      clsNPT=1;
      if (overlapMC(i, ierr))
	{
	  //printf("overlap di %d Lfact=%.15G\n", i, Lfact);
	  /* move rejected restore old positions */
#ifdef MC_STORE_ALL_COORDS
	  for (kk=0; kk < 3; kk++)
	    {
	      L[kk] = Lold[kk];
	      L2[kk] = 0.5*Lold[kk];
	    }
#else
	  L[dir] /= Lfact;
	  L2[dir] = L[dir]*0.5;
#endif
#if !defined(MC_STORE_ALL_COORDS)
	  for (ii=0; ii < Oparams.parnum; ii++)
	    {
	      rx[ii] -= Dxcls[color[ii]];
	      ry[ii] -= Dycls[color[ii]];
	      rz[ii] -= Dzcls[color[ii]];
	      rx[ii] += Dxpar[ii];
	      ry[ii] += Dypar[ii];
	      rz[ii] += Dzpar[ii]; 
#ifdef MC_BOND_POS
	      update_spot_pos(i, Dxpar[color[ii]]-Dxcls[color[ii]], Dypar[color[ii]]-Dycls[color[ii]],
			      Dzpar[color[ii]]-Dzcls[color[ii]]);
#endif
	      pbc(ii);
	    }
#else
	  restore_all_coords();
#endif

	  volrejMC++;
#ifdef MC_STORELL
	  restore_ll_mc();
#else
	  update_numcells();
	  rebuildLinkedList();
#endif
	  clsNPT=0;
	  return;
	}
      clsNPT=0;
    }
#ifndef MC_OPT_CLSNPT
#ifdef MC_STOREBONDS
  store_bonds_mc(-1);
#endif
#endif
  /* update all bonds with new positions for calculating new energy */
#ifndef MC_OPT_CLSNPT
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
    }
#endif
#ifdef MC_OPT_CLSNPT
  clsNPT=1;
#endif

  if (OprogStatus.useNNL)
    find_bonds_flex_NNL();
  else
    find_bonds_flex_all();

#ifdef MC_OPT_CLSNPT
  if (clsNPT==2)
    enn=eno-1;
  else
    enn=eno;
  clsNPT=0;
#else
  enn = calcpotene();
#endif
#if 0
  if (enn > eno)
    {
#if 0
      int *clsarr_bak, NP;
      int *clsdim_bak, kk;
      NP=Oparams.parnum;
      clsarr_bak = malloc(sizeof(int)*NP);
      clsdim_bak = malloc(sizeof(int)*NP);
      memcpy(clsarr_bak,clsarr,sizeof(int)*NP);
      memcpy(clsdim_bak,clsdim,sizeof(int)*NP);
      printf("CRITICAL: a cluster breaks, exiting... percolating=%d\n", percolating);
      printf("enn=%f eno=%f\n", enn, eno);
      build_clusters(&ncls, &percolating);
      printf("ncls=%d\n", ncls);
      exit(-1);
#endif
    }
#endif
  //printf("enn-eno/T=%f\n", (enn-eno)/Oparams.T);

#ifdef MC_CLSNPT_LOG
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-(ncls+1.0)*log(vn/vo)*Oparams.T);
#else
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-ncls*log(vn/vo)*Oparams.T);
#endif
#ifdef MC_ALMARZA
  omegaratio=calc_almarza_prob();
#if 0
  if (omegaratio != 1.)
    printf("omegaratio:%G enn-eno=%f\n", omegaratio, enn-eno);
#endif
#endif
  /* enn < eno means that a new bonds form, actually we have to reject such move */
  if (
#ifdef MC_ALMARZA
      ranf() > exp(arg)*omegaratio 
#else 
      ranf() > exp(arg) || enn != eno
#endif
      )
    {
      /* move rejected restore old positions */
#ifdef MC_STORE_ALL_COORDS
      for (kk=0; kk < 3; kk++)
	{
	  L[kk] = Lold[kk];
	  L2[kk] = 0.5*Lold[kk];
	}
#else
      L[dir] /= Lfact;
      L2[dir] = L[dir]*0.5;
#endif
#if !defined(MC_STORE_ALL_COORDS)
      for (i=0; i < Oparams.parnum; i++)
	{
	  rx[i] -= Dxcls[color[i]];
	  ry[i] -= Dycls[color[i]];
	  rz[i] -= Dzcls[color[i]];
	  rx[i] += Dxpar[i];
	  ry[i] += Dypar[i];
	  rz[i] += Dzpar[i]; 
#ifdef MC_BOND_POS
  	  update_spot_pos(i,  Dxpar[i]-Dxcls[color[i]], Dxpar[i]-Dxcls[color[i]], Dxpar[i]-Dxcls[color[i]]);
#endif 
	  pbc(i);
	}
#else
      restore_all_coords();
#endif
      volrejMC++;
#ifdef MC_STORELL
      restore_ll_mc();
#else
      update_numcells();
      rebuildLinkedList();
#endif
      /* restore all bonds*/
#ifndef MC_OPT_CLSNPT
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
#endif
      return;
    }
  if (OprogStatus.useNNL)
    {
      /* move accepted update displacements */
      for (i=0; i < Oparams.parnum; i++)
	{
	  nn = sqrt(Sqr(rx[i])+Sqr(ry[i])+Sqr(rz[i]));
	  overestimate_of_displ[i] += nn*delv; 
	}
    }
  //exit(-1);
}
#endif
#ifdef MC_CLUSTER_NPT
#define MC_CLSNPT_LOG
void move_box_cluster(int *ierr)
{
  int i0, i, ii, k, nc, ncls=0, percolating=0, np_in_cls, np, in0, in1, iold, kk;
  double nn, distance, lnvn;
  double del[3];
  double vo, vn, Lfact, enn, eno, arg, delv=0.0, lastrx=0, lastry=0, lastrz=0, Dx, Dy, Dz;
#ifdef MC_STORE_ALL_COORDS
#ifdef MD_LXYZ
  double Lold[3];
#else
  double Lold;
#endif
#endif
  //printf("moving box\n");
#ifdef MD_LXYZ
  vo = L[0]*L[1]*L[2];
#else
  vo = L*L*L;
#endif

  eno = calcpotene();
      
#ifdef MC_CLSNPT_LOG
  lnvn = log(vo) + (ranf()-0.5)*OprogStatus.vmax;
  vn = exp(lnvn);
#else
  vn = vo + (2.0*ranf()-1.0)*OprogStatus.vmax;
#endif
  Lfact = pow(vn/vo,1.0/3.0);
  if (OprogStatus.useNNL)
    delv = (Lfact - 1.0); /* FINISH HERE */
  //printf("Lfact=%.15G vmax=%f vn=%f vo=%f\n", Lfact, OprogStatus.vmax, vn, vo); 
  // Lfact=1;
#ifdef MC_ALMARZA
#ifdef MC_OPT_ALMARZA
  for (i=0; i < Oparams.parnum; i++)
    find_part_with_BS(i);
  for (i=0; i < Oparams.parnum; i++)
    {
      numcheckBS_MC[i] = numcheckBS[i];
      for (k=0; k < numcheckBS[i]; k++)
	{
	  checkBS_MC[i][k] = checkBS[i][k];
	}
    }
#endif 
  assign_bonds_almarza();
#endif
  build_clusters(&ncls, &percolating, 1);
  //printf("ncls=%d percolating=%d\n", ncls, percolating);
  /* calcola i centri di massa dei cluster */
  if (percolating)
    {
      volrejMC++;
      return;
    }
#ifdef MC_STORE_ALL_COORDS
  /* I use velocity arrays to store coordinates, in MC simulations they are unused! */
  store_all_coords();
#endif
  for (nc=0; nc < ncls; nc++)
    {
      for (k=0; k < 3; k++)
	clsCoM[k][nc] = 0.0;
      /* if it is percolating all particles have 2 bonds, 
	 hence we have to initialize i0 */
      i0 =  clsarr[firstofcls[nc]];
      for (np=0; np < clsdim[nc]; np++)
    	{
	  i =  clsarr[firstofcls[nc]+np];
	  if ( numbonds[i]==1) 
	    {
	      i0 = i;
	      //printf("qui i0=%d\n", i0);
	      break;
	    }
	}

      np_in_cls=1;
      i =  bonds[i0][0] / (NANA);
      iold=i0;	
      clsCoM[0][nc] = rx[i0];
      clsCoM[1][nc] = ry[i0];
      clsCoM[2][nc] = rz[i0]; 
#if !defined(MC_STORE_ALL_COORDS)
      Dxpar[i0]=Dypar[i0]=Dzpar[i0]=0;
#endif
      /* if particles has no bond we have to check... */ 
      while (numbonds[i0]>0)
	{
	  /* rebuild the cluster as a whole */
	  lastrx = rx[iold];
	  lastry = ry[iold];
	  lastrz = rz[iold];
#ifdef MC_STORE_ALL_COORDS
	  del[0] = L[0]*rint((rx[i]-lastrx)/L[0]);
	  del[1] = L[1]*rint((ry[i]-lastry)/L[1]); 
	  del[2] = L[2]*rint((rz[i]-lastrz)/L[2]);
	  rx[i] -= del[0];
	  ry[i] -= del[1];
	  rz[i] -= del[2];
#ifdef MC_BOND_POS
	  update_spot_pos(i, -del[0], -del[1], -del[2]);
#endif
#else
#ifdef MD_LXYZ
	  Dxpar[i] = L[0]*rint((rx[i]-lastrx)/L[0]);
	  Dypar[i] = L[1]*rint((ry[i]-lastry)/L[1]);
	  Dzpar[i] = L[2]*rint((rz[i]-lastrz)/L[2]); 
#else
	  Dxpar[i] = L*rint((rx[i]-lastrx)/L);
	  Dypar[i] = L*rint((ry[i]-lastry)/L);
	  Dzpar[i] = L*rint((rz[i]-lastrz)/L); 
#endif
    	  rx[i] -= Dxpar[i];
	  ry[i] -= Dypar[i];
	  rz[i] -= Dzpar[i];
#ifdef MC_BOND_POS
	  update_spot_pos(i, -Dxpar[i], -Dypar[i], -Dzpar[i]);
#endif
#endif
	  clsCoM[0][nc] += rx[i];
	  clsCoM[1][nc] += ry[i];
	  clsCoM[2][nc] += rz[i];
	  np_in_cls++;
#if 0
	  if (numbonds[i] > 2)
	    {
	      printf("more than 2 bonds (%d) for particle %d!\n", numbonds[i], i);
	    }
#endif
	  if(numbonds[i]==1 || np_in_cls == clsdim[nc])
	    {
	      //printf("numbonds[i=%d]=%d numbonds[i0=%d]=%d", i, numbonds[i], i0, numbonds[i0]);
	      break;
	    }	   
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
#if 0
      if (np_in_cls != clsdim[nc])
	{
	  printf("we got a problem np_in_cls=%d  but clsdim[%d]=%d\n", np_in_cls, nc, clsdim[nc]);
	  printf("percolating=%d\n", percolating);
	  exit(-1);
	}
#endif
    }
#ifdef MC_STORE_ALL_COORDS
#ifdef MD_LXYZ
  for (kk=0; kk < 3; kk++)
    Lold[kk] = L[kk];
#else
  Lold = L;
#endif
#endif
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

  for (nc=0; nc < ncls; nc++)
    {
      for (k=0; k < 3; k++)
	clsCoM[k][nc] /= clsdim[nc];
      Dxcls[nc]=(Lfact-1.0)*clsCoM[0][nc];
      Dycls[nc]=(Lfact-1.0)*clsCoM[1][nc];
      Dzcls[nc]=(Lfact-1.0)*clsCoM[2][nc];
    }
#if 0
  for (nc=0; nc < ncls; nc++)
    {
      for (np=0; np < clsdim[nc]; np++)
	{
	  i = clsarr[firstofcls[nc]+np];
	  rx[i] += Dxcls[color[i]];
	  ry[i] += Dycls[color[i]];
	  rz[i] += Dzcls[color[i]]; 
	  pbc(i);	
	}
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    {  
      rx[i] += Dxcls[color[i]];
      ry[i] += Dycls[color[i]];
      rz[i] += Dzcls[color[i]];
#ifdef MC_BOND_POS
      update_spot_pos(i, Dxcls[color[i]],Dycls[color[i]],Dzcls[color[i]]); 
#endif
#ifdef MC_STORE_ALL_COORDS
      pbc(i);
#else
      pbccls(i);
#endif
    }
#endif
#ifdef MC_STORELL
  store_ll_mc();
#endif
  update_numcells();
  rebuildLinkedList();
 
  for (i=0; i < Oparams.parnum; i++)
    {
      clsNPT=1;
      if (overlapMC(i, ierr))
	{
	  //printf("overlap di %d Lfact=%.15G\n", i, Lfact);
	  /* move rejected restore old positions */
#ifdef MC_STORE_ALL_COORDS
#ifdef MD_LXYZ
	  for (kk=0; kk < 3; kk++)
	    {
	      L[kk] = Lold[kk];
	      L2[kk] = 0.5*Lold[kk];
	    }
#else
	  L = Lold;
	  L2 = L*0.5;
#endif
#else
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
#endif
#if !defined(MC_STORE_ALL_COORDS)
	  for (ii=0; ii < Oparams.parnum; ii++)
	    {
	      rx[ii] -= Dxcls[color[ii]];
	      ry[ii] -= Dycls[color[ii]];
	      rz[ii] -= Dzcls[color[ii]];
	      
	      rx[ii] += Dxpar[ii];
	      ry[ii] += Dypar[ii];
	      rz[ii] += Dzpar[ii]; 
#ifdef MC_BOND_POS
	      update_spot_pos(i, Dxpar[color[ii]]-Dxcls[color[ii]], Dypar[color[ii]]-Dycls[color[ii]],
			      Dzpar[color[ii]]-Dzcls[color[ii]]);
#endif

	      pbc(ii);
	    }
#else
	restore_all_coords();
#endif

	  volrejMC++;
#ifdef MC_STORELL
	  restore_ll_mc();
#else
	  update_numcells();
	  rebuildLinkedList();
#endif
	  clsNPT=0;
	  return;
	}
      clsNPT=0;
    }
#ifndef MC_OPT_CLSNPT
#ifdef MC_STOREBONDS
  store_bonds_mc(-1);
#endif
#endif
  /* update all bonds with new positions */
#ifndef MC_OPT_CLSNPT
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
    }
#endif
#ifdef MC_OPT_CLSNPT
  clsNPT=1;
#endif
  if (OprogStatus.useNNL)
    find_bonds_flex_NNL();
  else
    find_bonds_flex_all();
#ifdef MC_OPT_CLSNPT
  if (clsNPT==2)
    enn=eno-1;
  else
    enn=eno;
  clsNPT=0;
#else
  enn = calcpotene();
#endif
#if 0
  if (enn > eno)
    {
      int *clsarr_bak, NP;
      int *clsdim_bak, kk;
      NP=Oparams.parnum;
      clsarr_bak = malloc(sizeof(int)*NP);
      clsdim_bak = malloc(sizeof(int)*NP);
      memcpy(clsarr_bak,clsarr,sizeof(int)*NP);
      memcpy(clsdim_bak,clsdim,sizeof(int)*NP);
      printf("CRITICAL: a cluster breaks, exiting... percolating=%d\n", percolating);
      printf("enn=%f eno=%f\n", enn, eno);
      build_clusters(&ncls, &percolating, 1);
      printf("ncls=%d\n", ncls);
      exit(-1);
    }
#endif
  //printf("enn-eno/T=%f\n", (enn-eno)/Oparams.T);
#ifdef MC_CLSNPT_LOG
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-(ncls+1.0)*log(vn/vo)*Oparams.T);
#else
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-ncls*log(vn/vo)*Oparams.T);
#endif
#ifdef MC_ALMARZA
  omegaratio=calc_almarza_prob();
#endif

  /* enn < eno means that a new bonds form, actually we have to reject such move */
  if (
#ifdef MC_ALMARZA
      ranf() > exp(arg)*omegaratio 
#else 
      ranf() > exp(arg) || enn != eno
#endif
      )
    {
      /* move rejected restore old positions */
#ifdef MC_STORE_ALL_COORDS
#ifdef MD_LXYZ
	  for (kk=0; kk < 3; kk++)
	    {
	      L[kk] = Lold[kk];
	      L2[kk] = 0.5*Lold[kk];
	    }
#else
	  L = Lold;
	  L2 = L*0.5;
#endif
#else
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
#endif
#if !defined(MC_STORE_ALL_COORDS)
      for (i=0; i < Oparams.parnum; i++)
	{
      	  rx[i] -= Dxcls[color[i]];
	  ry[i] -= Dycls[color[i]];
	  rz[i] -= Dzcls[color[i]];
	  rx[i] += Dxpar[i];
	  ry[i] += Dypar[i];
	  rz[i] += Dzpar[i]; 
#ifdef MC_BOND_POS
  	  update_spot_pos(i,  Dxpar[i]-Dxcls[color[i]], Dxpar[i]-Dxcls[color[i]], Dxpar[i]-Dxcls[color[i]]);
#endif 
	  pbc(i);
	}
#else
      restore_all_coords();
#endif
      volrejMC++;
#ifdef MC_STORELL
      restore_ll_mc();
#else
      update_numcells();
      rebuildLinkedList();
#endif
      /* restore all bonds*/
#ifndef MC_OPT_CLSNPT
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
#endif
      return;
    }
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
#endif
#ifdef MC_FREEZE_BONDS
extern int refFB, fakeFB;
#endif
#ifdef MC_NPT_XYZ
void move_box_xyz(int *ierr)
{
  int k1, extraspots, k2;
#ifdef MC_BOND_POS
  double oldrr, del;
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
#ifdef MD_LL_BONDS
  int nb, kk;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, aa, bb;
#endif  
#else
  int kk;
#endif
  int i, ii, dir;
  double nn;
#if 1
  double vo, lnvn, vn, Lfact, enn, eno, arg, delv=0.0;
#ifdef MC_STORE_ALL_COORDS
  double Lold[3];
#endif
  double *rr[3];
  //printf("moving box\n");
#ifdef MD_LXYZ
  vo = L[0]*L[1]*L[2];
#else
  vo = L*L*L;
#endif

  rr[0] = rx;
  rr[1] = ry;
  rr[2] = rz;
  dir = 3*ranf();
  eno = calcpotene();
  lnvn = log(vo) + (ranf()-0.5)*OprogStatus.vmax;
  vn = exp(lnvn);
  Lfact = vn/vo;
  if (OprogStatus.useNNL)
    delv = (Lfact - 1.0); /* FINISH HERE */
#ifdef MC_STORE_ALL_COORDS
  for (kk=0; kk < 3; kk++)
    Lold[kk] = L[kk];
  store_all_coords();
#endif
 //	Lfact=1;
  L[dir] *= Lfact;
  L2[dir] = L[dir]*0.5;
 for (i=0; i < Oparams.parnum; i++)
    {
#ifdef MC_BOND_POS
      oldrr = rr[dir][i];
#endif
      rr[dir][i] *= Lfact;
#ifdef MC_BOND_POS
      del = rr[dir][i] - oldrr;
      update_spot_pos_dir(i, dir, del);
#endif
      
#if 0
      rA[0] = rx[i];
      rA[1] = ry[i];
      rA[2] = rz[i];

      BuildAtomPos(i, rA, R[i], ratA);
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[typeOfPart[i]].nspotsBS;
#endif
      for (k1 = 0; k1 < typesArr[typeOfPart[i]].nspots + extraspots; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    bpos[k2][i][k1] = ratA[k1+1][k2];
	}
#endif
      pbc(i);
    }
  update_numcells();
#ifdef MC_STORELL
  store_ll_mc();
#endif
  rebuildLinkedList();
  for (i=0; i < Oparams.parnum; i++)
    {
#if !defined(MC_BENT_DBLCYL)
      if (Lfact > 1.0)
	break;
#endif
      if (overlapMC(i, ierr))
	{
	  //printf("overlap di %d Lfact=%.15G\n", i, Lfact);
	  /* move rejected restore old positions */
#ifdef MC_STORE_ALL_COORDS
	  for (kk=0; kk < 3; kk++)
	    {
	      L[kk] = Lold[kk];
	      L2[kk] = Lold[kk]*0.5;
	    }
#else
	  L[dir] /= Lfact;
	  L2[dir] = L[dir]*0.5;

#endif	 
#ifndef MC_STORE_ALL_COORDS
	  for (ii=0; ii < Oparams.parnum; ii++)
	    {
#ifdef MC_BOND_POS
	      oldrr = rr[dir][ii];
#endif
	      rr[dir][ii] /= Lfact;
#ifdef MC_BOND_POS
	      del = rr[dir][ii] - oldrr;
	      update_spot_pos_dir(ii, dir, del);
#endif
 	      pbc(ii);
	    }
#else
       	  restore_all_coords();
#endif
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
#ifdef MC_FREEZE_BONDS
  if (!OprogStatus.freezebonds)
    {
#ifdef MC_STOREBONDS
      store_bonds_mc(-1);
#endif
      /* update all bonds with new positions */
      for (i=0; i < Oparams.parnum; i++)
	numbonds[i] = 0;
    }
#else
#ifdef MC_STOREBONDS
  store_bonds_mc(-1);
#endif
  /* update all bonds with new positions */
#ifdef MC_KERN_FRENKEL
  for (i=0; i < Oparams.parnum; i++)
    {
      nb = numbonds[i];
      for (kk = 0; kk < nb; kk++)
	{
	  jj = bonds[i][kk] / (NANA);
	  jj2 = bonds[i][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (aa == 3 && bb == 3)
    	    {
	      //remove_bond(jj, i, bb, aa);
	      remove_bond(i, jj, aa, bb);
	    }
	}
    }
#elif defined(MC_GAPDNA)
  for (i=0; i < Oparams.parnum; i++)
    {
      nb = numbonds[i];
      for (kk = 0; kk < nb; kk++)
	{
	  jj = bonds[i][kk] / (NANA);
	  jj2 = bonds[i][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (aa > 1 && bb > 1)
    	    {
	      //remove_bond(jj, i, bb, aa);
	      remove_bond(i, jj, aa, bb);
	    }
	}
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    numbonds[i] = 0;
#endif
#endif  
#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds)
    { 
      refFB=0;
      fakeFB=1;
    }
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
  rejectMove = 0;
  checkMoveKF=1;
#endif

  if (OprogStatus.useNNL)
    find_bonds_flex_NNL();
  else
    find_bonds_flex_all();

#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
  checkMoveKF = 0;
#endif	

#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds==1)
    fakeFB=0;
#endif
  enn = calcpotene();
#if 0
  if (enn > -840)
    {
      printf("buuuummm\n");
      exit(-1);
    }
#endif
  //printf("enn-eno/T=%f\n", (enn-eno)/Oparams.T);
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-(Oparams.parnum+1)*log(vn/vo)*Oparams.T);
#ifdef MC_FREEZE_BONDS
  if (ranf() > exp(arg)|| (OprogStatus.freezebonds && refFB==1))
#elif defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
  if (ranf() > exp(arg) || rejectMove)
#else
  if (ranf() > exp(arg))
#endif
    {
     
#ifdef MC_STORE_ALL_COORDS
      for (kk=0; kk < 3; kk++)
	{
	  L[kk] = Lold[kk];
	  L2[kk] = Lold[kk]*0.5;
	}
#else
      /* move rejected restore old positions */
      L[dir] /= Lfact;
      L2[dir] = L[dir]*0.5;
#endif 
#ifndef MC_STORE_ALL_COORDS
      //printf("qui\n");
      for (i=0; i < Oparams.parnum; i++)
	{
#ifdef MC_BOND_POS
      	  oldrr = rr[dir][i];
#endif
	  rr[dir][i] /= Lfact;
#ifdef MC_BOND_POS
	  del = rr[dir][i] - oldrr;
    	  update_spot_pos_dir(i, dir, del);
#endif
 	  pbc(i);
	}
#else
      restore_all_coords();
#endif
      volrejMC++;
      update_numcells();
#ifdef MC_STORELL
      restore_ll_mc();
#else
      rebuildLinkedList();
#endif
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds)
	refFB=0;
      if (!OprogStatus.freezebonds)
	{
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
  	}
#if 0
      else
   	{
 	  for (i=0; i < Oparams.parnum; i++)
 	    numbonds[i] = 0;
 	  if (OprogStatus.useNNL)
	    find_bonds_flex_NNL();
 	  else
 	    find_bonds_flex_all();
  	}
#endif
#else
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
#endif
      return;
    }
  else
    {
#if 0
      //if (Oparams.curStep==46)
      printf("step=%d volume change accepted rejectMove=%d energy=%f!\n", Oparams.curStep, rejectMove, calcpotene());
      for (i=0;i< Oparams.parnum;i++)
	numbonds[i]=0;
      find_bonds_flex_all();
      printf("step=%d volume change accepted rejectMove=%d energy=%f!\n", Oparams.curStep, rejectMove, calcpotene());
#endif
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
#endif

void move_box(int *ierr)
{
#ifdef MC_BOND_POS
  double del[3], oldrr[3];
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
#ifdef MD_LL_BONDS
  int nb, kk;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, aa, bb;
#endif  
#endif
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
#ifdef MC_BOND_POS
      oldrr[0] = rx[i];
      oldrr[1] = ry[i];
      oldrr[2] = rz[i];
#endif
      rx[i] *= Lfact;
      ry[i] *= Lfact;
      rz[i] *= Lfact; 
#ifdef MC_BOND_POS
      del[0] = rx[i] - oldrr[0];
      del[1] = ry[i] - oldrr[1];
      del[2] = rz[i] - oldrr[2];
      update_spot_pos(i, del[0], del[1], del[2]);
#endif
      pbc(i);
    }
  update_numcells();
#ifdef MC_STORELL
  store_ll_mc();
#endif
  rebuildLinkedList();
  for (i=0; i < Oparams.parnum; i++)
    {
#if !defined(MC_BENT_DBLCYL)
      if (Lfact > 1.0)
	break;
#endif
      if (overlapMC(i, ierr))
	{
	  //printf("overlap di %d Lfact=%.15G\n", i, Lfact);
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
#ifdef MC_BOND_POS
	      oldrr[0] = rx[ii];
	      oldrr[1] = ry[ii];
	      oldrr[2] = rz[ii];
#endif
	      rx[ii] /= Lfact;
	      ry[ii] /= Lfact;
	      rz[ii] /= Lfact; 
#ifdef MC_BOND_POS
	      del[0] = rx[ii] - oldrr[0];
	      del[1] = ry[ii] - oldrr[1];
	      del[2] = rz[ii] - oldrr[2];
	      update_spot_pos(ii, del[0], del[1], del[2]);
#endif
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
#ifdef MC_FREEZE_BONDS
  if (!OprogStatus.freezebonds)
    {
#ifdef MC_STOREBONDS
      store_bonds_mc(-1);
#endif
      /* update all bonds with new positions */
      for (i=0; i < Oparams.parnum; i++)
	numbonds[i] = 0;
    }
#else
#ifdef MC_STOREBONDS
  store_bonds_mc(-1);
#endif
  /* update all bonds with new positions */
#ifdef MC_KERN_FRENKEL
  for (i=0; i < Oparams.parnum; i++)
    {
      nb = numbonds[i];
      for (kk = 0; kk < nb; kk++)
	{
	  jj = bonds[i][kk] / (NANA);
	  jj2 = bonds[i][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (aa == 3 && bb == 3)
    	    {
	      //remove_bond(jj, i, bb, aa);
	      remove_bond(i, jj, aa, bb);
	    }
	}
    }
#elif defined(MC_GAPDNA)
  for (i=0; i < Oparams.parnum; i++)
    {
      nb = numbonds[i];
      for (kk = 0; kk < nb; kk++)
	{
	  jj = bonds[i][kk] / (NANA);
	  jj2 = bonds[i][kk] % (NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  if (aa > 1 && bb > 1)
    	    {
	      //remove_bond(jj, i, bb, aa);
	      remove_bond(i, jj, aa, bb);
	    }
	}
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    numbonds[i] = 0;
#endif
#endif  
#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds)
    { 
      refFB=0;
      fakeFB=1;
    }
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      rejectMove = 0;
      checkMoveKF=1;
#endif

  if (OprogStatus.useNNL)
    find_bonds_flex_NNL();
  else
    find_bonds_flex_all();

#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
  checkMoveKF = 0;
#endif	

#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds==1)
    fakeFB=0;
#endif
  enn = calcpotene();
  //printf("enn-eno/T=%f\n", (enn-eno)/Oparams.T);
  arg = -(1.0/Oparams.T)*((enn-eno)+Oparams.P*(vn-vo)-(Oparams.parnum+1)*log(vn/vo)*Oparams.T);
#ifdef MC_FREEZE_BONDS
  if (ranf() > exp(arg)|| (OprogStatus.freezebonds && refFB==1))
#elif defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
  if (ranf() > exp(arg) || rejectMove)
#else
  if (ranf() > exp(arg))
#endif
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
#ifdef MC_BOND_POS
	  oldrr[0] = rx[i];
	  oldrr[1] = ry[i];
	  oldrr[2] = rz[i];
#endif
	  rx[i] /= Lfact;
	  ry[i] /= Lfact;
	  rz[i] /= Lfact; 
#ifdef MC_BOND_POS
    	  del[0] = rx[i] - oldrr[0];
	  del[1] = ry[i] - oldrr[1];
	  del[2] = rz[i] - oldrr[2];
	  update_spot_pos(i, del[0], del[1], del[2]);
#endif
	  pbc(i);
	}
      volrejMC++;
      update_numcells();
#ifdef MC_STORELL
      restore_ll_mc();
#else
      rebuildLinkedList();
#endif
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds)
	refFB=0;
      if (!OprogStatus.freezebonds)
	{
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
  	}
#if 0
      else
   	{
 	  for (i=0; i < Oparams.parnum; i++)
 	    numbonds[i] = 0;
 	  if (OprogStatus.useNNL)
	    find_bonds_flex_NNL();
 	  else
 	    find_bonds_flex_all();
  	}
#endif
#else
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
#if defined(MC_KERN_FRENKEL)||defined(MC_GAPDNA)
extern int checkMoveKF;
#endif
void update_bonds_MC(int ip)
{
#ifdef MD_LL_BONDS
  int nb, nn, kk;
  long long int jj, jj2, aa, bb;
#else
  int nb, jj, jj2, kk, nn, aa, bb;
#endif  

#ifdef MC_FREEZE_BONDS
  if (!OprogStatus.freezebonds)
    {
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
    }
#else
  nb = numbonds[ip];

#ifdef MC_KERN_FRENKEL

#ifdef MD_LL_BONDS
  memcpy(bondscache2, bonds[ip], sizeof(long long int)*numbonds[ip]);
#else
  memcpy(bondscache2, bonds[ip], sizeof(int)*numbonds[ip]);
#endif
  for (kk = 0; kk < nb; kk++)
    {
      jj = bondscache2[kk] / (NANA);
      jj2 = bondscache2[kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (checkMoveKF==1)
	{
#ifdef MC_AMYLOID_FIBRILS
	  /* 04/12/14 BUG FIX: all spots but first two are kern frenkel ones! */
	  if (aa >= 3 && bb >= 3)
	    {
	      remove_bond(jj, ip, bb, aa);
	      remove_bond(ip, jj, aa, bb);
	    }

#else
	  if (aa == 3 && bb == 3)
	    {
	      remove_bond(jj, ip, bb, aa);
	      remove_bond(ip, jj, aa, bb);
	    }
#endif
	}
      else
	remove_bond(jj, ip, bb, aa);
    }
#elif defined(MC_GAPDNA)
#ifdef MD_LL_BONDS
  memcpy(bondscache2, bonds[ip], sizeof(long long int)*numbonds[ip]);
#else
  memcpy(bondscache2, bonds[ip], sizeof(int)*numbonds[ip]);
#endif
  for (kk = 0; kk < nb; kk++)
    {
      jj = bondscache2[kk] / (NANA);
      jj2 = bondscache2[kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      //printf("========>aa=%lld bb=%lld\n", aa, bb);
      if (checkMoveKF==1)
	{
	  /* 04/12/14 BUG FIX: all spots but first two are kern frenkel ones! */
	  if (aa > 1 && bb > 1)
	    {
	      remove_bond(jj, ip, bb, aa);
	      remove_bond(ip, jj, aa, bb);
	    }
	}
      else
	remove_bond(jj, ip, bb, aa);
    }
#else
  for (kk = 0; kk < nb; kk++)
    {
      jj = bonds[ip][kk] / (NANA);
      jj2 = bonds[ip][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      remove_bond(jj, ip, bb, aa);     
    }
#endif
    
#ifdef MC_KERN_FRENKEL
  if (!checkMoveKF)
    numbonds[ip] = 0;
#elif defined(MC_GAPDNA)
  if (!checkMoveKF)
    numbonds[ip] = 0;
#else
  numbonds[ip] = 0;
#endif
#endif

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
#ifdef MC_CLUSTER_NPT
      if (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3)
#elif defined(MC_NPT_XYZ)
      if (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4)
#else
      if (OprogStatus.ensembleMC==1)
#endif
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
#ifdef MC_KERN_FRENKEL
double find_bonds_covadd_kf(int i, int j)
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
      /* consider only kerne-frenkel bonds */
      if (dists[nn] < 0.0 && mapbondsaFlex[nn] == 3 && mapbondsbFlex[nn] == 3)
	{
	  add_bond(i, j, mapbondsaFlex[nn], mapbondsbFlex[nn]);
	  add_bond(j, i, mapbondsbFlex[nn], mapbondsaFlex[nn]);
	}
    }

  //if (dist < 0)
    //printf("inside dist=%.15G\n", dist);
  return dist;
  
}

#endif
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
double find_bonds_fake(int i, int j, int *nbf);
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
  int nb, ierr, bonded, k1, k2, trials, nbold, nbf;
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
      else if (dist_type==0||dist_type==6)
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
      ene = find_bonds_fake(i, j, &nbf);
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
#if defined(MC_HC) || defined(MC_BIFUNC_SPHERES)
int check_bond_added(int j, int nb)
{
#ifdef MD_LL_BONDS
  long long int aa, bb, ii, jj, jj2;
#else
  int ii, jj, aa, bb, jj2;
#endif
  /* controlla che il bond aggiunto sia effettivamente il bond nb della particella j */
  jj = bonds[j][numbonds[j]-1] / (NANA);
  jj2 = bonds[j][numbonds[j]-1] % (NANA);
  aa = jj2 / NA;
  bb = jj2 % NA;
  if (aa-1!=nb)
    {
      return 0;
    }
  return 1;
}
#endif

#ifdef MC_BENT_DBLCYL
void addRestrMatrix(double Rl[3][3]);
#endif
void save_conf_mc(int i, int ii);
#ifdef MC_ELASTIC_CONSTANTS
double *ec_segno, *ec_ux, *ec_uy, *ec_uz;
void orient_donsager(double *ox, double *oy, double *oz, double alpha, double *segno); 
int cur_ii, cur_jj;
#endif
void mcin(int i, int j, int nb, int dist_type, double alpha, int *merr, int fake)
{
  /* dist_type=0 -> isotropic
     dist_type=1 -> onsager */
#ifdef MC_BENT_DBLCYL
  double *spXYZ=NULL;
  double splab[3], spR[3], rO[3];
#endif
  double sphrad, rmin1, rmin2, rmin, rmax, drSq;
  const int maxtrials=1000000;
  double bondlen, dist=0.0, rA[3], rat[3], norm, normB, sax, cc[3], ene;
#ifdef MCIN_OPT
#ifndef MD_SPOT_GLOBAL_ALLOC
  double ratAll[NA][3];
#endif
  double rB[3], normo;
  int nbB;
#endif
  double shift[3], Rl[3][3], vv[3];
  double ox, oy, oz, d, dx, dy, dz;
  int nbf, ierr, bonded, k1, k2, trials, nbold;
  static double tottrials=0, calls=0;
  /* N.B. questa routine dipende dalla geometria degli spot, in particolare il volume
     in cui inserire va scelto in relazione alla geometria adottata, per ora qui assumiamo
     due spot sulle basi di uno pseudo-cilindro di elongazione data elongazione (sax, vedi sotto) */
  /* place particle i bonded to bond nb of particle j */
  /* N.B. Questa routine presuppone di avere due spot lungo l'asse di simmetria della particella
    e che tutte le particelle siano uguali */

  rA[0] = rx[j];
  rA[1] = ry[j];
  rA[2] = rz[j];
#ifdef MCIN_OPT
  rB[0] = rx[i];
  rB[1] = ry[i];
  rB[2] = rz[i];
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  RtA[k1][k2] = R[j][k1][k2];
#ifdef MCIN_OPT
	  RtB[k1][k2] = R[i][k1][k2];
#endif
	}
    }
#if 0
#ifdef MCIN_OPT
  BuildAtomPos(i, rB, RtB, ratAll);
#endif
#endif
  BuildAtomPosAt(j, nb+1, rA, RtA, rat);
  for (k1=0; k1 < 3; k1++)
    vv[k1] = rat[k1] - rA[k1];
  norm = calc_norm(vv);
  for (k1=0; k1 < 3; k1++)
    vv[k1] /=norm;
  assign_bond_mapping(i,j);
  if (are_spheres(i, j))
    {
      sax = (norm+mapSigmaFlex[0])*1.05;
      for (k1=0; k1 < 3; k1++)
	cc[k1] = rat[k1];
    }
  else
    {
      sax = typesArr[typeOfPart[j]].sax[0];
      for (k1=0; k1 < 3; k1++)
	cc[k1] = rA[k1] + vv[k1]*sax*2.0;
    }
  /* N.B. here we assume one bond per pair */
  if (are_spheres(i,j))
    bondlen = mapSigmaFlex[0];
  else 
    {
      bondlen = 1.05*2.0*(norm+0.5*mapSigmaFlex[0]-sax); /* overestimate */
    }
  //printf("bondlen=%.15G norm=%.15G mapSismaFlex=%.15G\n", bondlen, norm, mapSigmaFlex[0]);
  bonded=0;
  trials=0;
  if (!fake)
    numbonds[i]=0;
  *merr=0;
  while (!bonded)
    {
#ifdef MCIN_OPT
      //nbB = ranf()*typesArr[typeOfPart[i]].nspots;
      nbB = (ranf()>0.5)?1:-1;
#if 0
      for (k1=0; k1 < 3; k1++)
	vv[k1] = ratAll[nbB+1][k1] - rB[k1];
      normB = calc_norm(vv);
#endif
#endif
#if 1
#ifdef MCIN_OPT
      do {
	dx = 2.0*(ranf_vb()-0.5);
	dy = 2.0*(ranf_vb()-0.5);
	dz = 2.0*(ranf_vb()-0.5);
      } 
      while (dx*dx+dy*dy+dz*dz > 1);
#ifdef MC_BIFUNC_SPHERES
      dx = dx*OprogStatus.distKF;
      dy = dy*OprogStatus.distKF;
      dz = dz*OprogStatus.distKF;
#else
      dx = dx*mapSigmaFlex[0];
      dy = dy*mapSigmaFlex[0];
      dz = dz*mapSigmaFlex[0];
#endif
      //printf("nbB=%d dx=%f %f %f\n",nbB, dx, dy, dz);
#else
#if 0
	/* N.B. qui si assume che gli spot siano uguali per tutte le particelle, che siano lungo x 
	   e che siano simmetrici rispetto al centro di massa della particella */

      if (are_spheres(i,j))
	{
      	  rmax = (norm+mapSigmaFlex[0])*1.01;
	  //rmin1 = 0.5*(typesArr[typeOfPart[i]].sax[0]+typesArr[typeOfPart[j]].sax[0])*0.99;
       	  rmin = (norm-mapSigmaFlex[0])*0.99;
	  //rmin = max(rmin1, rmin2);
	  //printf("rmin=%.15G rmax=%.15G\n", rmin, rmax);
	  do {
	    dx = 2.0*(ranf_vb()-0.5);
	    dy = 2.0*(ranf_vb()-0.5);
	    dz = 2.0*(ranf_vb()-0.5);
	    //printf("sphrad=%.15G\n", sphrad);
	    dx *= rmax;
	    dy *= rmax;
	    dz *= rmax;
	    drSq = Sqr(dx)+Sqr(dy)+Sqr(dz);
	  }
	  while (drSq > Sqr(rmax) || drSq < Sqr(rmin));
	  
	}
      else
	{
	  rmax = (norm+mapSigmaFlex[0])*1.01;
	  //rmin = (typesArr[typeOfPart[i]].sax[0]+typesArr[typeOfPart[j]].sax[0])*0.99;
       	  rmin = (norm-mapSigmaFlex[0])*0.99;
	  //rmin = max(rmin1, rmin2);
	  do {
	    dx = 2.0*(ranf_vb()-0.5);
	    dy = 2.0*(ranf_vb()-0.5);
	    dz = 2.0*(ranf_vb()-0.5);
	    //printf("sphrad=%.15G\n", sphrad);
	    dx *= rmax;
	    dy *= rmax;
	    dz *= rmax;
	    drSq = Sqr(dx)+Sqr(dy)+Sqr(dz);
	  }
	  while (drSq > Sqr(rmax) || drSq < Sqr(rmin));
	}
      rx[i] = rat[0]+dx;
      ry[i] = rat[1]+dy;
      rz[i] = rat[2]+dz;
#else
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
#endif
#endif
#if 1
#if !defined(MC_HC) 
      if (!are_spheres(i,j) && Sqr(rx[j]-rx[i])+Sqr(ry[j]-ry[i])+Sqr(rz[j]-rz[i]) > Sqr(2.0*sax+bondlen))
	{
	  trials++;
  	  continue;
	}
#endif
#endif

#ifndef MCIN_OPT
#ifndef MC_CALC_COVADD
      pbc(i);
#endif
#endif
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
      else if (dist_type==0||dist_type==6)
	{
  	  orient(&ox, &oy, &oz);
	}
#ifdef MC_ELASTIC_CONSTANTS
      else if (dist_type == 11 || dist_type == 12 || dist_type == 13)
	{
	  if (i == cur_ii || i == cur_jj)
	    {
	      orient_donsager(&ox, &oy, &oz, alpha, &(ec_segno[i]));
	      ec_ux[i] = ox;
	      ec_uy[i] = oy;
	      ec_uz[i] = oz;
	    }
	  else
	    {
	      orient_onsager(&ox, &oy, &oz, alpha);
	      ec_segno[i] = 1.0;
	      ec_ux[i] = ox;
	      ec_uy[i] = oy;
	      ec_uz[i] = oz;
	    }
	}
#endif
      else 
	{
  	  //orient(&ox, &oy, &oz);
	  orient_onsager(&ox, &oy, &oz, alpha);
	}
#ifdef MCIN_OPT
#ifndef MC_BENT_DBLCYL
#ifdef MC_BIFUNC_SPHERES
      rx[i] = rA[0] + dx;
      ry[i] = rA[1] + dy;
      rz[i] = rA[2] + dz;
#else
      rx[i] = dx + ox*norm + rat[0];
      ry[i] = dy + oy*norm + rat[1];
      rz[i] = dz + oz*norm + rat[2];
#endif
#if 0
      ox = rat[0] + dx - rx[i];
      oy = rat[1] + dy - ry[i];
      oz = rat[2] + dz - rz[i];
      normo = sqrt(Sqr(ox)+Sqr(oy)+Sqr(oz));
      ox /= normo;
      oy /= normo;
      oz /= normo;
#else
      /* nbB = +1 o -1 così viene scelta a caso uno dei due spot */
      ox = nbB*ox;
      oy = nbB*oy;
      oz = nbB*oz;
#endif
#ifndef MC_CALC_COVADD
      pbc(i);
#endif
#endif
#endif
      versor_to_R(ox, oy, oz, Rl);
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
      if (dist_type != 0 && dist_type != 4 && dist_type !=6)
	{
  	  addRestrMatrix(Rl);
	}
#endif
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      R[i][k1][k2] = Rl[k1][k2]; 
	    }
	}
#ifdef MC_BENT_DBLCYL 
      /* pick one spot randomly */
      nbB =(ranf()>0.5)?0:1;
      spXYZ = typesArr[typeOfPart[i]].spots[nbB].x;
      /* we want the position of the spot bonded to spot of particle j */
      body2labR(i, spXYZ, spR, NULL, R[i]);
      splab[0] = dx + rat[0];
      splab[1] = dy + rat[1];
      splab[2] = dz + rat[2];
      /* center of mass position rcm has to be consistent with the position of the spot 
	 i.e  spot_position_Lab = R*spot_position_Body + rcm => rcm = spot_position_Lab - R*spot_position_Body  */
      rx[i] = splab[0] - spR[0];
      ry[i] = splab[1] - spR[1];
      rz[i] = splab[2] - spR[2];
      pbc(i);
#endif
#if 0
      store_bonds_mc(-1);
#else
      //store_bonds_mc(j);
      nbold = numbonds[j];
#endif

      if (fake)
	{
#if !defined(MCIN_OPT) || defined(MC_BIFUNC_SPHERES)
	  ene = find_bonds_fake(i, j, &nbf);
#else
	  /* N.B. the optimized version always create a bonded conf
	     we have only to check for overlaps between hard cores */
#if 0
	  ene = find_bonds_fake(i, j, &nbf);
	  if (ene >= 0)
	    exit(-1);
	  else
	    printf("ene=%f nbB=%d\n", ene,nbB);
#else
	  /* con il metodo ottimizzato il legame si forma sempre quindi non
	     si devono controllare i legami */
	  ene = -1;
#endif
#endif
	}
      else
	{
#if defined(MCIN_OPT) && !defined(MC_BIFUNC_SPHERES)
	  /* con il metodo ottimizzato creiamo un legame tra gli spot
          (i, nbB) e (j, nb) */
#if defined(MC_BENT_DBLCYL)
	  ene = -1;
	  add_bond(i, j, nbB+1, nb+1);
      	  add_bond(j, i, nb+1, nbB+1);

#if 0
	  if (numbonds[i] > 1)
	    printf("BOH numbonds=%d nbold=%d\n", numbonds[i], nbold);
	  dist=find_bonds_fake(i, j, &nbf);
	  if (dist > 0)
	    {
    	      printf("[WARNING] MCIN FAILED PARTICLES NOT BONDED!\n");
	    }
#endif	    
#else
	  dist=find_bonds_covadd(i, j);
	  ene = calcpotene_GC(i);
	  if (ene >= 0)
	    {
	      printf("[WARNING] MCIN FAILED PARTICLES NOT BONDED!\n");
	    }
#endif
#else
	  dist=find_bonds_covadd(i, j);
	  ene = calcpotene_GC(i);
	  if (ene >= 0)
	    {
#ifndef MC_BIFUNC_SPHERES
	      printf("[WARNING] MCIN FAILED PARTICLES NOT BONDED!\n");
#endif
	      //exit(-1);
#if 0
	      rO[0] = rx[i];
	      rO[1] = ry[i];
	      rO[2] = rz[i];
	      body2lab(i, spXYZ, spR, rO, R[i]);
	      printf("spot i=%f %f %f spot j= %f %f %f dist=%f\n", rat[0], rat[1], rat[2], spR[0], spR[1], spR[2], sqrt(Sqr(rat[0]-spR[0])+
															Sqr(rat[1]-spR[1])+Sqr(rat[2]-spR[2])) );
	      printf("dx=%f dy=%f dz=%f\n", dx, dy, dz);
	      printf("qui ene=%f\n", ene);
#endif
	    }
#endif
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
#if defined(MC_HC) || defined(MC_BIFUNC_SPHERES) 
	if (bonded)
	  {
	    if (!fake)
	      {
		/* 16/11/13: nel metodo ottimizzato formo un legame sicuramente tra (j,nb) e uno dei due spot di i 
		   quindi non si deve controllare nulla */
#if !defined(MCIN_OPT) || defined(MC_BIFUNC_SPHERES)
		if (!check_bond_added(j, nb))
		  {
		    numbonds[j] = nbold;
		    numbonds[i] = 0;
		    bonded=0;
		    //printf("BOH\n");
		  }
#endif
#if !defined(MCIN_OPT) || defined(MC_BIFUNC_SPHERES)
    		else
     		  {
		    /* se il bond trovato non è quello prescelto allora continua a provare */
#ifndef MC_BIFUNC_SPHERES
		    if (nbf != nb)
    		      bonded=0;	
#endif
     		  } 
#endif
	      }
	  }
#endif 
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
#if 0
// DEBUGGING STUFF
   if (Oparams.curStep%10==0)
    printf("trials=%d avg_trials=%.15G\n", trials, tottrials/calls);
#endif
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
    return; /* tutti i siti sono occupati ( assumendo che ogni sito è legato al massimo con un altro sito )*/
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
  /* inserisce la particella legata a j con il bond nb appena scelto.
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
  if (OprogStatus.ensembleMC==1 || OprogStatus.ensembleMC ==3)
    {
#ifdef MD_LXYZ
      vnonbond = L[0]*L[1]*L[2] - OprogStatus.vbond;
#else
      vnonbond = L*L*L - OprogStatus.vbond;
#endif
    }
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
double find_bonds_fake(int i, int j, int *nbf)
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
	  *nbf = mapbondsbFlex[nn]-1;
	}
    }
  return ene;

}
void mcout(int i, int j, int nb)
{
  double ene, ox, oy, oz, Rl[3][3];
  int k1, k2, bonded=1, nbjold;
  int nbf;
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
      ene=find_bonds_fake(i, j, &nbf);
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
  if (OprogStatus.ensembleMC==1 || OprogStatus.ensembleMC ==3)
    {
#ifdef MD_LXYZ
      vnonbond = L[0]*L[1]*L[2] - OprogStatus.vbond;
#else
      vnonbond = L*L*L - OprogStatus.vbond;
#endif
    }
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
void save_conf_mc(int i, int ii)
{
  char fileop2[512], fileop[512], fileop3[512];
  FILE* bf;
  sprintf(fileop2 ,"Store-%d-%d", 
	  i, ii);
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
#ifdef MC_BENT_DBLCYL
void adjust_r1r2(int i, double rcom[3], double r1[3], double r2[3])
{
  double dr1[3], dr2[3], norm1, norm2;
  int k;
  for (k=0; k < 3; k++)
    {
      dr1[k] = r1[k] - rA[k];
      dr2 [k] = r2[k] - rA[k]; 
    }
  norm1= calc_norm(dr1);
  norm2= calc_norm(dr2);
  for (k=0; k < 3; k++)
    {
      dr1[k] *= typesArr[typeOfPart[i]].sax[0]/norm1;
      dr2[k] *= typesArr[typeOfPart[i]].sax[0]/norm2;
    }
  for (k=0; k < 3; k++)
    {
      dr1[k] += rcom[k];
      dr2[k] += rcom[k];
    }
}
void accum_persist_len_mono_orient_base_base(int *parlist, double *pl, double *cc)
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
      adjust_r1r2(i, rA, r1, r2);  
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

	  adjust_r1r2(j, rA, r1, r2);  
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
#endif

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
	      /* we assume here that bonds 0 and 1 are permanent ones! */
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
		  //save_conf_mc(tt, 0);
		  break;
		}
	    }
#endif
	  if (abort)
	    break;
	}
#if 0
	if (!abort)
	  save_conf_mc(tt, 0);
#else
      if (!abort)
	{
	  /* N.B. size1 viene solo usato per scegliere il modo in cui viene calcolata 
	     la funzione di correlazione dei bond */
	  if (size1!=0)
	    {
#ifdef MC_BENT_DBLCYL
	      if (size1==2)
		accum_persist_len_mono_orient_base_base(parlist, pl, cc);
	      else
		accum_persist_len_mono_orient(parlist, pl, cc);
#else
	      accum_persist_len_mono_orient(parlist, pl, cc);
#endif
	    }
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
#if defined(MC_KERN_FRENKEL) && !defined(MC_BIFUNC_SPHERES)
int insert_two_polymers_kf(int size1, int size2)
{
  int i, j, nb, bt;
  /* insert firs cluster */
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
		  save_conf_mc(0, 0);
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
	  save_conf_mc(0, 0);
	  printf("[mcin] attempt to add a bonded particle failed!\n");
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
		      save_conf_mc(0, 0);
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
	      save_conf_mc(0, 0);
	      printf("[mcin] attempt to add a bonded particle failed!\n");
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
	  return 0;
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

  /* 0 means failed, i.e. overlap */
  if (selfoverlap||merr||overlap)
    return 0;
  else
    return 1;
}

void calc_bonding_volume_kf(long long int maxtrials, int outits, int type, double alpha)
{
  /* BONDING VOLUME BETWEEN TWO SEMI-FLEXIBLE POLYMERS */
  FILE *f;	
  int k1, k2, ierr, i;
  long long int tt;
  double deldistcc=0.0, deltotdist=0.0, dist, fact, shift[3], Lb, totbonds=0.0, totene=0.0, ox, oy, oz, Rl[3][3];
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
#ifdef MD_LXYZ
	  fprintf(f, "%lld %.15G %G\n", tt, (totene/((double)tt))*(L[0]*L[1]*L[2])/Sqr(size1), totene);
#else
	  fprintf(f, "%lld %.15G %G\n", tt, (totene/((double)tt))*(L*L*L)/Sqr(size1), totene);
#endif
	  sync();
	  fclose(f);
	}
      for (i=0; i < Oparams.parnum; i++)
	{
	  numbonds[i] = 0;
	}
      if (insert_two_polymers_kf(size1, size2)==0)
	{
	  tt++;
	  continue;
	}
      for (i=0; i < Oparams.parnum; i++)
	{
	  numbonds[i] = 0;
	}
      deltotdist=0.0;
      deldistcc=0.0;
      for (i=0; i < size1; i++) 
	{
	  for (j=i+1; j < Oparams.parnum; j++) 
	    dist=find_bonds_covadd_kf(i, j);
	  if (dist < 0.0)
	    {
	      deltotdist += dist;
    	      deldistcc += 1.0;
	    }
	}
      totbonds=0;
      for (i=0; i < Oparams.parnum; i++)
	totbonds+=numbonds[i]; 
      
      if (totbonds > 0)
	{	
    	  totene += totbonds/2.0;		  
	  totdist += deltotdist;
	  distcc += deldistcc;
	}
#if 0
      if (tt%500000==0)
	{
	  printf("tt=%lld\n", tt); 
	}
#endif
      /* we are calculating the avaerage bonding volume per sphere */
#ifdef MD_LXYZ
      printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)tt))*(L[0]*L[1]*L[2])/Sqr(size1), totene);
#else
      printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)tt))*(L*L*L)/Sqr(size1), totene);
#endif
      //fclose(f);
#ifdef MC_VB_PLOT
      fclose(vbf);
#endif
    }
}
#endif
void calc_bonding_volume_mc(long long int maxtrials, int outits, int type, double alpha)
{
  FILE *f;	
  int k1, k2, ierr, i;
  long long int tt;
  double Pi,deldistcc=0.0, deltotdist=0.0, dist, fact, shift[3], Lb, totene=0.0, ox, oy, oz, Rl[3][3];
#ifdef MC_VB_PLOT
  FILE *vbf;
#endif
  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
#ifdef MC_VB_PLOT
  vbf = fopen("bonding-volume.dat", "w+");
#endif
  if (type==9)
    printf("calc B2 MC\n");
  else
    printf("calc vbonding MC!!\n");
  Pi = 2.0*acos(0.0);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
      }
  if (type==9)
    f = fopen("B2.dat","w+");
  else
    f = fopen("vbonding.dat","w+");
  fclose(f);
  printf("toteneini=%f\n", toteneini);
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
	  if (type==9)
	      f=fopen("B2.dat", "a");
	  else
	      f=fopen("vbonding.dat", "a");
	  fact=Oparams.parnum-1;
	 if (type==9)
	   {
#ifdef MD_LXYZ
	     fprintf(f, "%lld %.15G %G\n", tt, 0.5*(totene/((double)tt))*(L[0]*L[1]*L[2])/fact, totene);
#else
   	     fprintf(f, "%lld %.15G %G\n", tt, 0.5*(totene/((double)tt))*(L*L*L)/fact, totene);
#endif
	   }
	 else
	   {
#ifdef MD_LXYZ
	     fprintf(f, "%lld %.15G %G\n", tt, (totene/((double)tt))*(L[0]*L[1]*L[2])/Sqr(typesArr[0].nspots)/fact, totene);
#else
   	     fprintf(f, "%lld %.15G %G\n", tt, (totene/((double)tt))*(L*L*L)/Sqr(typesArr[0].nspots)/fact, totene);
#endif
	   }
	  sync();
	  fclose(f);
	}
      /* nel seguito deve inserire due polimeri di uguale lunghezza */
      
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
	    	      if (type==9)
	    		{
	    		  totene += 1.0-exp(1.0/Oparams.T);
	    		}
	    	      else
	    		{
			  totene += 1.0;
			}
		    }
		} 
	      else
		{
		  if (type==9)
		    totene += 1.0;
		}
	    }
	}
      else
	{
	  if (type==9 && (check_overlap(0,1, shift, &ierr) <= 0.0))
	    totene += 1.0;
      	}
	      
#if 0
      if (tt%500000==0)
	{
	  printf("tt=%lld\n", tt); 
	}
#endif
    }
  if (type==9)
    {
#ifdef MD_LXYZ
      printf("B2=%.10f (totene=%f)\n", 0.5*(totene/((double)tt))*(L[0]*L[1]*L[2]), totene);
#else
      printf("B2=%.10f (totene=%f)\n", 0.5*(totene/((double)tt))*(L*L*L), totene);
#endif
    }
  else
    {
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
    }
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
#if 0
      	  printf("i=%d j=%d overlap!\n", i, j);
	  printf("shift=%f %f %f\n", shift[0], shift[1], shift[2]);
	  printf("r=%f %f %f  j %f %f %f\n", rx[i], ry[i], rz[i], rx[j], ry[j], rz[j]);
	  printf("self overlap!\n");
#endif
	  return 1;
	}
    }
   return 0; 
}
#define MC_COV_SIG
#ifdef MC_COV_SIG
double check_over(int size1, int size2, double xi, double rhx, double rhy, double rhz)
{
  int i, j, overlap, ierr;  
  double shift[3], posx, posy, posz, rCMx, rCMy, rCMz;
  /* check overlaps */
  rCMx = rCMy = rCMz = 0;
  for (i=size1; i < Oparams.parnum; i++)
    {
      rCMx += rx[i];
      rCMy += ry[i];
      rCMz += rz[i];
    }	
  rCMx /= size2;
  rCMy /= size2;
  rCMz /= size2;

  posx = rhx*xi;
  posy = rhy*xi;
  posz = rhz*xi;

  for (i=size1; i < Oparams.parnum; i++)
    {
      rx[i] = (rx[i] - rCMx) + posx;
      ry[i] = (ry[i] - rCMy) + posy;
      rz[i] = (rz[i] - rCMz) + posz;
    }	
  for (i=size1; i < Oparams.parnum; i++)
    { 
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
  if (overlap)
    return -1;
  else
    return 1;
}
double do_bisection(double rmin, double rmax, int size1, int size2, double rhx, double rhy, double rhz, int *found)
{
  /* if odd = 1 look for odd solutions */
  const int maxit=500;
  int cc;
  double rmid; // I am looking for E such that Psi(xmax)=0
  double k1,k2,k3; // values of Psi(xmax) for Emin,Emax,E
  double tol=1.0e-6; // tolerance in finding a root

  k1 = check_over(size1, size2, rmin, rhx, rhy, rhz);
  k2 = check_over(size1, size2, rmax, rhx, rhy, rhz);
  if( k1*k2 > 0.0 )
    {
      //fprintf(stderr,"Don't know if there is  a zero in [%g,%g]\n",Emin,Emax);
      //fprintf(stderr, "k1=%-15G k2=%.15G\n", k1, k2);
      //fflush(stderr); 
      printf("no bracketing in bisect.\n");
      *found = 0;
      return 0.0; 
    }

  cc=0;
  do
    {
      rmid = 0.5*(rmin+rmax); k3 = check_over(size1, size2, rmid, rhx, rhy, rhz);
      if( k3*k2 < 0.0 ) rmin=rmid; 
      else rmax=rmid;
      cc++;    
      if (cc > maxit)
	{
	  printf("Too many iterations in bisections\n");
	  *found = 0;
	  return 0;
	}
    } 
  while(fabs(rmax-rmin)>tol);

  //printf("Bisection finished after %d iterations\n", cc);
  *found=1;
  return rmid;
}
#ifdef MC_ELASTIC_CONSTANTS
double dfons_sinth_max, fons_sinth_max;
/* return an angle theta sampled from an Onsager angular distribution */
double dfons(double theta, double alpha);
double theta_donsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta, dtheta;
  //printf("alpha=%f\n", alpha);
  pi = acos(0.0)*2.0;
  if (first == 1)
    {
      first=0;
      /* qui va fornito il massimo della funzione nell'intervallo
	 [0,Pi] */
      //f0 = 1.01*alpha*fons(0.0,alpha);
      f0 = 1.01*dfons_sinth_max;
    }
  do 
    {
      /* uniform theta between 0 and pi */
      theta = ranf_vb()*pi;
      /* uniform y between 0 and 1 (note that sin(theta) <= 1 for 0 < theta < pi)*/
      y = f0*ranf_vb();
      f = fabs(sin(theta)*dfons(theta,alpha));
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}

/* first derivative of Onsager distribution */
double dfons(double theta, double alpha)
{
  double pi;
  pi = acos(0.0)*2.0;
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return sinh(alpha*cos(theta))*alpha*alpha/(4.0*pi*sinh(alpha));
}

void orient_donsager(double *omx, double *omy, double* omz, double alpha, double *segno)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_donsager(alpha);
  if (thons < pi*0.5)
    *segno = 1.0;
  else
    *segno = -1.0;
  //printf("thos=%f\n", thons);
  //distro[(int) (thons/(pi/((double)nfons)))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1;
#if 1 /* along z */
  *omx = verso*sin(thons)*cos(phi);
  *omy = verso*sin(thons)*sin(phi);
  *omz = verso*cos(thons); 
#else /* or along x (but it has to be same of course!) */
  *omy = verso*sin(thons)*cos(phi);
  *omz = verso*sin(thons)*sin(phi);
  *omx = verso*cos(thons); 
#endif
  //printf("norma=%f\n", sqrt(Sqr(*omx)+Sqr(*omy)+Sqr(*omz)));
}
const double thetapts=100000;

double estimate_maximum_dfons(double alpha)
{
  double th, dth, maxval, m;
  int i;
  dth=2.0*(acos(0.0))/((double)thetapts);
  th=0.0;
  for (i=0; i < thetapts; i++)
    {
      m=sin(th)*dfons(th,alpha);
      if (i==0 || maxval < m)
	maxval = m;
      th += dth;
      //printf("%f %.15G\n", th, sin(th)*dfons(th, alpha));
    }
  // printf("maxval=%f\n", maxval);
  return maxval;
}
void calc_com_cls(double Rcm1[3], double Rcm2[3], int size1)
{
  int i0, i, kk, in0, in1, iold;
  double lastrx, lastry, lastrz;
  Rcm1[0]=Rcm1[1]=Rcm1[2]=0;
  Rcm2[0]=Rcm2[1]=Rcm2[2]=0;
  
#if 0
  for (i = 0; i < size1; i++)
    if (numbonds[i]==1)
      {
	i0 = i;
	break;
      }

  i =  bonds[i0][0] / (NANA);
  iold=i0;
  Rcm1[0] = rx[i0];
  Rcm1[1] = ry[i0];
  Rcm1[2] = rz[i0];
   /* if particles has no bond we have to check... */ 
  while (numbonds[i0]>0)
    {
      /* rebuild the cluster as a whole */
      lastrx = rx[iold];
      lastry = ry[iold];
      lastrz = rz[iold];
      rx[i] -= L[0]*rint((rx[i]-lastrx)/L[0]);
      ry[i] -= L[1]*rint((ry[i]-lastry)/L[1]);
      rz[i] -= L[2]*rint((rz[i]-lastrz)/L[2]); 
      Rcm1[0] += rx[i];
      Rcm1[1] += ry[i];
      Rcm1[2] += rz[i];
#if 0
      if (numbonds[i] > 2)
	{
	  printf("more than 2 bonds (%d) for particle %d!\n", numbonds[i], i);
	}
#endif
      if (numbonds[i]==1)
	{
	  //printf("numbonds[i=%d]=%d numbonds[i0=%d]=%d", i, numbonds[i], i0, numbonds[i0]);
	  break;
	}	   
      /* qui assumo che le particelle siano bifunzionali */
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
  for (i = size1; i < Oparams.parnum; i++)
    if (numbonds[i]==1)
      {
	i0 = i;
	break;
      }
  i =  bonds[i0][0] / (NANA);
  iold=i0;
  Rcm2[0] = rx[i0];
  Rcm2[1] = ry[i0];
  Rcm2[2] = rz[i0];
   /* if particles has no bond we have to check... */ 
  while (numbonds[i0]>0)
    {
      /* rebuild the cluster as a whole */
      lastrx = rx[iold];
      lastry = ry[iold];
      lastrz = rz[iold];
      rx[i] -= L[0]*rint((rx[i]-lastrx)/L[0]);
      ry[i] -= L[1]*rint((ry[i]-lastry)/L[1]);
      rz[i] -= L[2]*rint((rz[i]-lastrz)/L[2]); 
      Rcm2[0] += rx[i];
      Rcm2[1] += ry[i];
      Rcm2[2] += rz[i];
#if 0
      if (numbonds[i] > 2)
	{
	  printf("more than 2 bonds (%d) for particle %d!\n", numbonds[i], i);
	}
#endif
      if (numbonds[i]==1)
	{
	  //printf("numbonds[i=%d]=%d numbonds[i0=%d]=%d", i, numbonds[i], i0, numbonds[i0]);
	  break;
	}	   
      /* qui assumo che le particelle siano bifunzionali */
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
#endif
#if 1
  for (i=0; i < size1; i++)
    {
      Rcm1[0]+=rx[i];
      Rcm1[1]+=ry[i];
      Rcm1[2]+=rz[i];
      Rcm2[0]+=rx[i+size1];
      Rcm2[1]+=ry[i+size1];
      Rcm2[2]+=rz[i+size1];
    }
#endif
  for (kk=0; kk < 3; kk++)
    {
      Rcm1[kk] /= ((double)size1);
      Rcm2[kk] /= ((double)size1);
      //Rcm1[kk] -= L[kk]*rint(Rcm1[kk]/L[kk]);
      //Rcm2[kk] -= L[kk]*rint(Rcm2[kk]/L[kk]);
    }
  
}
void calc_elastic_constants(int type, double alpha, int maxtrials, int outits, int size1)
{
  int tt, k1, k2, i, j, kk, merr, selfoverlap=0, overlap, bt, nb, ierr;
  int ii, jj;
  double delrxii, delryii, delrzii, delrxjj, delryjj, delrzjj;
  const int nn=1000;
   FILE *f;
  double cov, totene, shift[3], segno;
  double ox, oy, oz, Rl[3][3], Rcm1[3], Rcm2[3], Rcm[3], fact1, fact2; 
  if (type == 11)
    f = fopen("v11.dat","w+");
  else if (type==12)
    f = fopen("v22.dat","w+");
  else
    f = fopen("v33.dat","w+");

  distro=malloc(sizeof(double)*nn);
  for (i=0; i < nn; i++)
    distro[i] = 0.0;
  
  fclose(f);
  ec_segno = malloc(sizeof(double)*Oparams.parnum);
  ec_ux =    malloc(sizeof(double)*Oparams.parnum);
  ec_uy =    malloc(sizeof(double)*Oparams.parnum);
  ec_uz =    malloc(sizeof(double)*Oparams.parnum);
  dfons_sinth_max=estimate_maximum_dfons(alpha);
  fons_sinth_max=dfons_sinth_max/alpha;
  printf("Estimated maximum of dfons is %f maxtrials=%d alpha=%f\n", dfons_sinth_max, maxtrials, alpha);
  tt=0;
  totene=0.0;
  while (tt < maxtrials) 
    {
#if 0
      for (ii=0; ii < Oparams.parnum; ii++)
	for (jj=size1; jj < Oparams.parnum; jj++)
	  { 
	    if (ii==jj)
	      continue;


	    if (ii >= size1 && jj >= size1)
	      segno=1.0;
	    else 
	      segno=-1.0;
#endif
	    cur_ii=-1; 
	    cur_jj=-1;
	    while (cur_ii==cur_jj)
	      {
		cur_ii = (int)(ranf()*Oparams.parnum);
#if 1
		cur_jj = (int)(ranf()*Oparams.parnum);
#else
		cur_jj = size1+((int)(ranf()*size1));
#endif
		//	printf("qui cur_ii=%d cur_jj=%d\n", cur_ii, cur_jj);
	      }
	    //printf("qui cur_ii=%d cur_jj=%d\n", cur_ii, cur_jj);
#if 1
	    segno = 1.0;
#else
	    if (cur_ii >= size1 && cur_jj >= size1)
	      segno=1.0;
	    else 
	      segno=-1.0;
#endif
	    segno *= 0.5*(Oparams.parnum-1)*Oparams.parnum;
	    merr=ierr=0;
	    selfoverlap = overlap = 0;
	    //cur_ii = ii;
	    //cur_jj = jj;
	    /* first particle is always in the center of the box with the same orientation */
	    rx[0] = 0;
	    ry[0] = 0;
	    rz[0] = 0;
	    if (cur_ii==0)
	      orient_donsager(&ox, &oy, &oz, alpha, &(ec_segno[0])); 
	    else
	      {
		orient_onsager(&ox, &oy, &oz, alpha);
		ec_segno[0] = 1.0;
	      }
	    ec_ux[0] = ox;
	    ec_uy[0] = oy;
	    ec_uz[0] = oz;
	    //ox = 0; oy=0; oz=1;
	    //printf("alpha=%f ox=%f oy=%f oz=%f norm=%f\n", alpha, ox, oy, oz, sqrt(Sqr(ox)+ Sqr(oy)+Sqr(oz)));
	    versor_to_R(ox, oy, oz, Rl);
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
	    addRestrMatrix(Rl);
#endif
	    for (k1=0; k1 < 3; k1++)
	      for (k2=0; k2 < 3; k2++)
		{
		  R[0][k1][k2] = Rl[k1][k2];
		}
	    /* place first cluster */
	    if (tt%outits==0)// && ii==0 && jj==size1)
	      {
		if (tt!=0)
		  {
#ifdef MD_LXYZ
		    cov = (totene/((double)tt))*(L[0]*L[1]*L[2]);

#else
		    cov = (totene/((double)tt))*(L*L*L);
#endif
		    cov *= Sqr(alpha);
		    if (type == 11)
		      f = fopen("v11.dat","a");
		    else if (type==12)
		      f = fopen("v22.dat","a");
		    else
		      f = fopen("v33.dat","a");

		    if (type == 11)
		      printf("v11=%.10f (totene=%f/%lld)\n", cov, totene, tt);
		    else if (type==12)
		      printf("v22=%.10f (totene=%f/%lld)\n", cov, totene, tt);
		    else
		      printf("v33=%.10f (totene=%f/%lld)\n", cov, totene, tt);

		    fprintf(f, "%lld %.15G %.15G\n", tt, cov, totene);
		    fclose(f);
		    sync();
		  }
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
			    save_conf_mc(0, 0);
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
		    save_conf_mc(0, 0);
		    printf("[mcin covolume] attempt to add a bonded particle failed!\n");
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

	    if (selfoverlap||merr)
	      {
		//tt++;
		continue;
	      }
	    /* place second cluster */
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
		    /* NOTA 28/06/17:
		     * ora devo generalre una distribuzione che è la derivata della funzione di Onsager
		     * ma devo anche aggiustare i segni negativi che spuntano fuori! */
		    if (i == cur_jj)
		      {
			orient_donsager(&ox, &oy, &oz, alpha, &(ec_segno[i]));
		      }
		    else
		      {
			orient_onsager(&ox, &oy, &oz, alpha);
			ec_segno[i] = 1.0;
		      }
		    ec_ux[i] = ox;
		    ec_uy[i] = oy;
		    ec_uz[i] = oz;

		    versor_to_R(ox, oy, oz, Rl);
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
		    addRestrMatrix(Rl);
#endif
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
			//printf("i=%d j=%d size1=%d\n", i, j, size1);

#if 1
			if (bt > Oparams.parnum*1000)
			  {
			    printf("insertion cluster #2: maximum number of iteration reached\n");
			    if (fabs(calcpotene()+(i-size1)*2.0) < 1E-10)
			      {
				save_conf_mc(0, 0);
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
			save_conf_mc(0, 0);
			printf("[mcin] attempt to add a bonded particle failed!\n");
			break;
		      }
		    if (check_self_overlap(size1, i))
		      {
#if 0
			printf("BOH i=%d j=%d\n", i, j);
			save_conf_mc(tt, 0); 
#endif
			selfoverlap = 1;
			continue;
			//break;
		      }

		    /* N.B. per ora non controlla il self-overlap della catena 
		       e la formazione dopo mcin di legami multipli poiché
		       si presuppone che al massimo stiamo considerando dimeri */
		  }
#if 0
		if (selfoverlap)
		  {
		    break;
		    //tt++;
		    //continue;
		  }
#endif
		//overlap = 0;
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
			//break;
		      }
		  }
	      }
#if 0
	    if (selfoverlap||merr)
	      {
		//tt++;
#if 0	  
		save_conf_mc(tt, 0); 
#endif
		break;
	      }
#endif
	    /* qui va modificato perché ora calcoliamo le costanti elastiche */
	    if (overlap && selfoverlap==0 && merr==0 && ierr==0)
	      {
		//printf("qui\n");
		calc_com_cls(Rcm1, Rcm2, size1);
		for (kk=0; kk < 3; kk++)
		  {
		    Rcm[kk] = Rcm2[kk] - Rcm1[kk];
		  }
#if 1
		if (type==11) // K11
		  {
		    delrxii = rx[cur_ii]-Rcm1[0];
		    delryii = ry[cur_ii]-Rcm1[1];
		    delrxjj = rx[cur_jj]-Rcm1[0];
		    delryjj = ry[cur_jj]-Rcm1[1];
		    fact1 = ec_ux[cur_ii]*ec_ux[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj]*
		      (-2.0*delrxii*delrxjj + delrxii*delrxii + delrxjj*delrxjj);
		    fact2 = ec_uy[cur_ii]*ec_uy[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj]*
		      (-2.0*delryii*delryjj + delryii*delryii + delryjj*delryjj);
		    //(ry[cur_ii]-Rcm1[0])*(ry[cur_jj]-Rcm1[0]);
		    //printf("Rcm %f\n", Rcm[0]);
		    //printf("segno %f %f\n", ec_segno[cur_ii], ec_segno[cur_jj]);
		    totene += segno*(fact1+fact2)*0.5;
		  }
		else if (type==12) // K22
		  {
		    delrxii = rx[cur_ii]-Rcm1[0];
		    delryii = ry[cur_ii]-Rcm1[1];
		    delrxjj = rx[cur_jj]-Rcm1[0];
		    delryjj = ry[cur_jj]-Rcm1[1];
		    fact1 = ec_ux[cur_ii]*ec_ux[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    fact2 = ec_uy[cur_ii]*ec_uy[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    fact1 *= (-2.0*delryii*delryjj + delryii*delryii + delryjj*delryjj);
		    fact2 *= (-2.0*delrxii*delrxjj + delrxii*delrxii + delrxjj*delrxjj);
		    totene += segno*(fact1+fact2)*0.5;
		  }
		else if (type==13) // K33
		  {
		    delrzii = rz[cur_ii]-Rcm1[2];
		    delrzjj = rz[cur_jj]-Rcm1[2];
		    fact1 = ec_ux[cur_ii]*ec_ux[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    fact2 = ec_uy[cur_ii]*ec_uy[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    totene += segno*(fact1+fact2)*0.5*(-2.0*delrzii*delrzjj + delrzii*delrzii + delrzjj*delrzjj);
		  }
	      }
#else
		if (type==11) // K11
		  {
		    fact1 = ec_ux[cur_ii]*ec_ux[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj]*Rcm[0]*Rcm[0];
		    fact2 = ec_uy[cur_ii]*ec_uy[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj]*Rcm[1]*Rcm[1];
		    //printf("Rcm %f\n", Rcm[0]);
		    //printf("segno %f %f\n", ec_segno[cur_ii], ec_segno[cur_jj]);
		    totene += segno*(fact1+fact2)*0.5;
		  }
		else if (type==12) // K22
		  {
		    fact1 = ec_ux[cur_ii]*ec_ux[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    fact2 = ec_uy[cur_ii]*ec_uy[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    fact1 *= Rcm[1]*Rcm[1];
		    fact2 *= Rcm[0]*Rcm[0];
		    totene += segno*(fact1+fact2)*0.5;
		  }
		else if (type==13) // K33
		  {
		    fact1 = ec_ux[cur_ii]*ec_ux[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    fact2 = ec_uy[cur_ii]*ec_uy[cur_jj]*ec_segno[cur_ii]*ec_segno[cur_jj];
		    totene += segno*(fact1+fact2)*0.5*Rcm[2]*Rcm[2];
		  }
	      }
#endif
#if 0
	    if (tt!=0 && tt%outits==0)
	      save_conf_mc(tt, 0); 
#endif

	    if (ierr!=0)
	      {
		printf("COV main loop NR failure\n");
		printf("i=%d j=%d scalp=%.15G\n", i, j, R[i][0][0]*R[j][0][0]+R[i][0][1]*R[j][0][1]+R[i][0][2]*R[j][0][2]);
		store_bump(i,j);
	      }
#if 0
	  }
#endif
      tt++;
    }

}
#endif
void calc_sigma_parsons(int size1, int size2, double alpha, int type, int outits, int maxtrials)
{
  FILE *f;	 
  long long int tt=0;
  int ttt;
  double A, Aold;
  double Rl[3][3], posx, posy, posz, xi, xiold, sigma, sumsigma=0.0;
  int found_contact, selfoverlap, k1, k2, i, j;
  double rCMx, rCMy, rCMz, ox, oy, oz, rhx, rhy, rhz, rmin, rmax, rmid;
  int merr, ierr, bt=0, nb, found, overlap;

  f = fopen("sigma_parsons.dat","w+");
  fclose(f);
  
  while (tt < maxtrials) 
    {
      merr=ierr=0;
      selfoverlap=0;
      //printf("parsons tt=%d\n", tt);
      /* first particle is always in the center of the box with the same orientation */
      rx[0] = 0;
      ry[0] = 0;
      rz[0] = 0;
      if (type==7)
	orient_onsager(&ox, &oy, &oz, alpha); 
      else
	orient(&ox, &oy, &oz); 
      versor_to_R(ox, oy, oz, Rl);
      orient(&rhx, &rhy, &rhz); 
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    R[0][k1][k2] = Rl[k1][k2];
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
		      save_conf_mc(0, 0);
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
	      save_conf_mc(0, 0);
	      printf("[mcin] attempt to add a bonded particle failed!\n");
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
      /* place first cluster with the center of mass in the center
	 of the box */
      rCMx = rCMy = rCMz = 0;
      for (i=0; i < size1; i++)
	{
	  rCMx += rx[i];
	  rCMy += ry[i];
	  rCMz += rz[i];
	}	
      rCMx /= size1;
      rCMy /= size1;
      rCMz /= size1;
      for (i=0; i < size1; i++)
	{
     	  rx[i] = (rx[i] - rCMx);
	  ry[i] = (ry[i] - rCMy);
	  rz[i] = (rz[i] - rCMz);
	}	
      if (selfoverlap)
	{
	  printf("self overlap!\n");
	}
      /* place second cluster */
      overlap=0;
      for (i=size1; i < Oparams.parnum; i++)
	{
	  if (i==size1)
	    {
#if 0
#ifdef MD_LXYZ
	      rx[i] = L[0]*(ranf_vb()-0.5);
	      ry[i] = L[1]*(ranf_vb()-0.5); 
	      rz[i] = L[2]*(ranf_vb()-0.5); 
#else
    	      rx[i] = L*(ranf_vb()-0.5);
    	      ry[i] = L*(ranf_vb()-0.5); 
    	      rz[i] = L*(ranf_vb()-0.5); 
#endif
#endif
	      rx[i] = 0;
	      ry[i] = 0;
	      rz[i] = 0;
    
	      if (type==7)
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
			  save_conf_mc(0, 0);
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
		  save_conf_mc(0, 0);
		  printf("[mcin] attempt to add a bonded particle failed!\n");
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
	      printf("self overlap!\n");
	    }	  
	}

      /* place first cluster */
      if (tt%outits==0)
	{
	  if (tt!=0)
	    {
	      f=fopen("sigma_parsons.dat", "a");
	      printf("sigma=%.10f (tries=%lld)\n", 4.0*acos(0.0)*2*sumsigma / tt, tt);
	      fprintf(f, "%lld %.15G\n", tt, 4.0*acos(0.0)*2*sumsigma / tt);
	      fclose(f);
	      sync();
	    }
#if 0
     	  save_conf_mc(tt, 0); 
#endif
	}
      // tt++;
      // continue;
      xi = Oparams.parnum*(typesArr[typeOfPart[0]].sax[0]+mapSigmaFlex[0])*2.0;
      xiold = xi - 0.5*typesArr[typeOfPart[0]].sax[0]; 

      Aold = A = 0.0;
      found_contact=0;
      ttt=0;
      do
	{
	  //printf("found=%d xi=%.15G dxi=%.15G\n", found_contact, xi, typesArr[typeOfPart[0]].sax[0]);
	  if (xi < 0.0)
	    {
	      printf("contact bracketing not found...boh!\n");
	      //exit(-1);
	      break;
	    }
	 
	  Aold = A;
 	  if ((A=check_over(size1, size2, xi, rhx, rhy, rhz)) < 0.0)
	    {
	      //printf("A=%f Aold=%f\n", A, Aold);
	      sigma=do_bisection(xiold, xi, size1, size2, rhx, rhy, rhz, &found);
	      found_contact = 1;
	      break;
	    }
#if 0
	  save_conf_mc(tt, ttt);
#endif
	  xiold = xi;
	  xi -= 0.5*typesArr[typeOfPart[0]].sax[0]; 
	  ttt++;
	}
      while (!found_contact);
      tt++;
      sumsigma += sigma*sigma*sigma/3.0; 
   }

  printf("sigma=%.15G\n", 4.0*acos(0.0)*2*sumsigma / tt);
}
#endif
extern double fons(double theta, double alpha);
#ifdef MC_BENT_DBLCYL

extern double thetaGlobalBondangle;
void calc_stddev_angle(double alpha, long long int maxtrials, int outits, int size1)
{
  FILE *f, *f1;	
  int k1, k2, ierr, i, abort=0, nb, j, merr, jj;
  long long int tt;
  double deldistcc=0.0, deltotdist=0.0, dist, fact, shift[3], Lb, totene=0.0, ox, oy, oz, Rl[3][3];
  double *ad, deltheta, avgangle=0.0, stddevangle=0.0, sum=0.0, norm=0.0;
  int num_bins = 100;
  if (size1 > 0)
    num_bins=size1;
  printf("calculating bonding angle\n");
  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
  ox = 0; 
  oy = 0;
  oz = 1;
#if 0
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
      }
#endif
//   versor_to_R(ox, oy, oz, Rl);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Rl[k1][k2] = restrMatrix[k1][k2];
	//printf("Rl[%d][%d]=%f\n", k1, k2, restrMatrix[k1][k2]);
      }

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  R[0][k1][k2] = Rl[k1][k2]; 
	}
    }
  f = fopen("bondingangle.dat","w+");
  fclose(f);
  totene = toteneini;
  if (Oparams.parnum > 2 )
    {
      printf("For calculating the bonding azimuthal angle 2 particles are enough!\n");
      exit(-1);
    }
  //cc = malloc(sizeof(double)*num_bins);
  ad = malloc(sizeof(double)*(num_bins+1));
  for (i=0; i < num_bins; i++)
    {
      ad[i] = 0;
    }
  deltheta = 4.0*acos(0.0)/100.0;
  for (tt=0; tt < maxtrials; tt++)
    {
      abort=0;
      if (tt%outits==0)
	{
	  printf("trials=%lld\n", tt);
	  if (tt!=0)
	    {
	      f = fopen("bondingangle.dat", "a");
	      fprintf(f, "%lld %.15G %.15G\n", tt, avgangle/((double)tt), 
		      stddevangle/((double)tt)-Sqr(avgangle/((double)tt)));
	      fclose(f);
	    }
	  if (tt!=0)
	    {
	      f1 = fopen("bonddistro.dat", "w+");
	      norm=0.0;
	      for (i=0; i < num_bins; i++)
		norm += ad[i]*deltheta;
	      for (i=0; i < num_bins; i++)
		fprintf(f1,"%.15G %.15G\n", i*deltheta, ad[i]/norm);
	      fclose(f1);
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
	  mcin(i, j, nb, 1, alpha, &merr, 0);
	  if (merr!=0)
	    {
	      printf("[mcin] attempt to add a bonded particle failed!\n");
	      /* mc in ha fallito scarta tale conf */
	      abort=1;
	      break;
	    }
	  /* qui controlla che non ci siano overlap con le particelle
	     già inserite e che mcin abbia formato un solo legame */
#if 0
	  for (jj=0; jj < i; jj++)
	    {
	      if (jj==j)
		continue;
	      find_bonds_covadd(i, jj);
      	    }
#endif
#if 0
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
#endif
#if 0
#ifdef MD_LXYZ
	  shift[0] = L[0]*rint((rx[0]-rx[1])/L[0]);
	  shift[1] = L[1]*rint((ry[0]-ry[1])/L[1]);
	  shift[2] = L[2]*rint((rz[0]-rz[1])/L[2]);
#else
	  shift[0] = L*rint((rx[0]-rx[1])/L);
	  shift[1] = L*rint((ry[0]-ry[1])/L);
	  shift[2] = L*rint((rz[0]-rz[1])/L);
#endif
	  if (check_overlap(0, 1, shift, &ierr) < 0.0)
	    {
#if 0
	      restore_bonds_mc(-1);
	      numbonds[i]=0;
	      i--;
#endif
	      printf("qui\n");
	      /* scarta la conf self-overlapping */
	      abort=1; 
	    }
#endif
	  if (abort)
	    break;
	}
#if 0
      save_conf_mc(tt, 0);
#else
      if (!abort)
	{
	  /* N.B. size1 viene solo usato per scegliere il modo in cui viene calcolata 
	     la funzione di correlazione dei bond */
	  /* >>>calc_angle_avg_stddev <<< */
	  //printf("numbonds[...]=%d %d\n", numbonds[0], numbonds[1]);
	  i = (int) (thetaGlobalBondangle/deltheta);
	  ad[i]++;
	  //save_conf_mc(tt,0);
	  avgangle += thetaGlobalBondangle;
	  stddevangle += Sqr(thetaGlobalBondangle);
	}
#endif
    }
}
#endif
#ifdef MC_BENT_DBLCYL
void addRestrMatrix(double Rl[3][3])
{
  int k1, k2, k3;
  double Rll[3][3];
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rll[k1][k2] = 0;
	for (k3 = 0; k3 < 3; k3++)
	  //Ro[k1][k2] += M[k1][k3]*R[i][k3][k2];
	  Rll[k1][k2] += restrMatrix[k1][k3]*Rl[k3][k2];
	  //Rll[k1][k2] += Rl[k1][k3]*restrMatrix[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      Rl[k1][k2] = Rll[k1][k2];
}
#endif
#ifdef MC_SMECTIC_VEXCL
void calc_cov_additive(void)
{
  FILE *fi, *f=NULL;
  double cov, Lb, totene = 0.0, alpha, shift[3];
  int bt=0, merr=0, i, j=-1, selfoverlap=0, size1, size2, nb, k1, k2, overlap=0, ierr, type, outits;
  long long int maxtrials, tt;
  double ox, oy, oz, Rl[3][3];

  fi = fopen("covmc.conf", "r");
  fscanf(fi, "%lld %d %d %d ", &maxtrials, &type, &size1, &outits);

  printf("type=%d\n", type);
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i]=0;
    }
  
  printf("ene iniziale=%f\n", calcpotene());
  fscanf(fi, " %lf ", &alpha);
  printf("type=%d alpha=%.15G\n", type, alpha);
    
  /* if it is a restart put info initial values here */
  covrestart = 0;
  if (!feof(fi))
    {
      covrestart=1;
      fscanf(fi, "%lf %lld ", &toteneini, &ttini);
    }
  /* type = 0 -> covolume smectic unfolded 
     type = 1 -> covolume smectic folded
   */
  const int n=1000;
  distro=malloc(sizeof(double)*n);
  for (i=0; i < n; i++)
    distro[i] = 0.0;
  
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
  /* first particle is always in the center of the box with the same orientation */
  rx[0] = 0;
  ry[0] = 0;
  rz[0] = 0;
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	R[0][k1][k2] = (k1==k2)?1:0;
      }

  if (type==0)
    f = fopen("covolume-sm-unfold.dat","w+");
  else if (type==1)
    f = fopen("covolume-sm-fold.dat","w+");
  fclose(f);
    
  while (tt < maxtrials) 
    {
      merr=ierr=0;
      selfoverlap=0;

      /* first particle is always in the center of the box with the same orientation */
      rx[0] = 0;
      ry[0] = 0;
      rz[0] = 0;
      orient_onsager(&ox, &oy, &oz, alpha); 
      //ox = 0; oy=0; oz=1;
      //printf("alpha=%f ox=%f oy=%f oz=%f norm=%f\n", alpha, ox, oy, oz, sqrt(Sqr(ox)+ Sqr(oy)+Sqr(oz)));
      versor_to_R(ox, oy, oz, Rl);
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
      addRestrMatrix(Rl);
#endif
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    R[0][k1][k2] = Rl[k1][k2];
	  }
      /* place first cluster */
      if (tt%outits==0)
	{
	  if (tt!=0)
	    {
#ifdef MD_LXYZ
	      cov = (totene/((double)tt))*(L[1]*L[2]);
	    
#else
	      cov = (totene/((double)tt))*(L*L);
#endif
	      if (type==0)
		f=fopen("covolume-sm-unfold.dat", "a");
	      else
		f=fopen("covolume-sm-fold.dat", "a");
	      printf("co-volume=%.10f (totene=%f/%lld)\n", cov, totene, tt);
	      fprintf(f, "%lld %.15G %.15G\n", tt, cov, totene);
	      fclose(f);
	      sync();
	    }
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
		      save_conf_mc(0, 0);
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
	      save_conf_mc(0, 0);
	      printf("[mcin covolume] attempt to add a bonded particle failed!\n");
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
	      /* since it is smectic we put second particle onto plane x=0 */
	      rx[i] = 0.0;//L[0]*(ranf_vb()-0.5);//L[0]*(ranf_vb()-0.5);
	      ry[i] = L[1]*(ranf_vb()-0.5); 
	      rz[i] = L[2]*(ranf_vb()-0.5); 
#else
	      /* since it is smectic we put second particle onto plane x=0 */
    	      rx[i] = 0.0;//L*(ranf_vb()-0.5);
    	      ry[i] = L*(ranf_vb()-0.5); 
    	      rz[i] = L*(ranf_vb()-0.5); 
#endif
	      orient_onsager(&ox, &oy, &oz, alpha);
    	      versor_to_R(ox, oy, oz, Rl);
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
    	      addRestrMatrix(Rl);
#endif
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
		  //printf("i=%d j=%d size1=%d\n", i, j, size1);

#if 1
	    	  if (bt > Oparams.parnum*1000)
	    	    {
		      printf("insertion cluster #2: maximum number of iteration reached\n");
		      if (fabs(calcpotene()+(i-size1)*2.0) < 1E-10)
			{
			  save_conf_mc(0, 0);
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
		  save_conf_mc(0, 0);
		  printf("[mcin] attempt to add a bonded particle failed!\n");
		  break;
		}
	      if (check_self_overlap(size1, i))
		{
#if 0
		  printf("BOH i=%d j=%d\n", i, j);
		  save_conf_mc(tt, 0); 
#endif
		  selfoverlap = 1;
		  break;
		}

	      /* N.B. per ora non controlla il self-overlap della catena 
		 e la formazione dopo mcin di legami multipli poiché
		 si presuppone che al massimo stiamo considerando dimeri */
	    }
	  if (selfoverlap||merr)
	    {
	      break;
	      //tt++;
	      //continue;
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

      if (selfoverlap||merr)
	{
	  tt++;
#if 0	  
	  save_conf_mc(tt, 0); 
#endif
	  continue;
	}
      if (overlap)// && ierr==0)
	totene += 1.0;
#if 0
      if (tt!=0 && tt%outits==0)
	save_conf_mc(tt, 0); 
#endif

      if (ierr!=0)
	{
	  printf("COV main loop NR failure\n");
	  printf("i=%d j=%d scalp=%.15G\n", i, j, R[i][0][0]*R[j][0][0]+R[i][0][1]*R[j][0][1]+R[i][0][2]*R[j][0][2]);
	  store_bump(i,j);
	}
      tt++;
   }
  if (type==1||type==0)
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
  printf("co-volume=%.10f (totene=%f)\n", (totene/((double)tt))*(L[1]*L[2]), totene);
  printf("%.15G\n",(totene/((double)tt))*(L[1]*L[2]));
#else 
  printf("co-volume=%.10f (totene=%f)\n", (totene/((double)tt))*(L*L), totene);
  printf("%.15G\n",(totene/((double)tt))*(L*L));
#endif
  exit(-1);
}

#else
void calc_cov_additive(void)
{
  FILE *fi, *f=NULL;
  double cov, Lb, totene = 0.0, alpha, shift[3];
  int bt=0, merr=0, i, j=-1, selfoverlap=0, size1, size2, nb, k1, k2, overlap=0, ierr, type, outits;
  long long int maxtrials, tt;
  double ox, oy, oz, Rl[3][3];

  fi = fopen("covmc.conf", "r");
  fscanf(fi, "%lld %d %d %d ", &maxtrials, &type, &size1, &outits);

  printf("type=%d\n", type);
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i]=0;
    }
  
  printf("ene iniziale=%f\n", calcpotene());
  if (type==1||type==5||type==8||type==11||type==12||type==13)
    {
      fscanf(fi, " %lf ", &alpha);
      printf("type=%d alpha=%.15G\n", type, alpha);
    }
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
     type = 8 -> standard deviation of bonding angle in the nematic phase
     type = 9 -> B2
     type = 10 -> K11
     type = 11 -> K22
     type = 12 -> K33
   */
 if (type==1||type==5||type==7||type==8)
    {
      const int n=1000;
      distro=malloc(sizeof(double)*n);
      for (i=0; i < n; i++)
	distro[i] = 0.0;
    }
  
  if (size1 >= Oparams.parnum && !type==8)
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
#ifdef MC_COV_SIG
  if (type==6||type==7)
    {
      calc_sigma_parsons(size1, size2, alpha, type, outits, maxtrials);
      exit(-1);
    }
#endif
#ifdef MC_BENT_DBLCYL
  if (type == 8) 
    {
      calc_stddev_angle(alpha, maxtrials, outits, size1);
      exit(-1);
    }
#endif
#ifdef MC_ELASTIC_CONSTANTS
  if (type == 11 || type == 12 || type == 13)
    {
      printf("type=%d alpha=%.15G size1=%d outits=%d\n", type, alpha, size1, outits);
      calc_elastic_constants(type, alpha, maxtrials, outits, size1);
      exit(-1);
    }
#endif
  if (type==3||type==5||type==9)
    {
#if defined(MC_KERN_FRENKEL) && !defined(MC_BIFUNC_SPHERES)
	  /* calc_poly_bonding_volume */
      calc_bonding_volume_kf(maxtrials, outits, type, alpha);
      exit(-1);
#else
      calc_bonding_volume_mc(maxtrials, outits, type, alpha);
      exit(-1);
#endif
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
	  //ox = 0; oy=0; oz=1;
	  //printf("alpha=%f ox=%f oy=%f oz=%f norm=%f\n", alpha, ox, oy, oz, sqrt(Sqr(ox)+ Sqr(oy)+Sqr(oz)));
	  versor_to_R(ox, oy, oz, Rl);
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
	  addRestrMatrix(Rl);
#endif
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
		      save_conf_mc(0, 0);
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
	      save_conf_mc(0, 0);
	      printf("[mcin covolume] attempt to add a bonded particle failed!\n");
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
#if defined(MC_CALC_COVADD) && defined(MC_BENT_DBLCYL)
    	      addRestrMatrix(Rl);
#endif
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
		  //printf("i=%d j=%d size1=%d\n", i, j, size1);

#if 1
	    	  if (bt > Oparams.parnum*1000)
	    	    {
		      printf("insertion cluster #2: maximum number of iteration reached\n");
		      if (fabs(calcpotene()+(i-size1)*2.0) < 1E-10)
			{
			  save_conf_mc(0, 0);
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
		  save_conf_mc(0, 0);
		  printf("[mcin] attempt to add a bonded particle failed!\n");
		  break;
		}
	      if (check_self_overlap(size1, i))
		{
#if 0
		  printf("BOH i=%d j=%d\n", i, j);
		  save_conf_mc(tt, 0); 
#endif
		  selfoverlap = 1;
		  break;
		}

	      /* N.B. per ora non controlla il self-overlap della catena 
		 e la formazione dopo mcin di legami multipli poiché
		 si presuppone che al massimo stiamo considerando dimeri */
	    }
	  if (selfoverlap||merr)
	    {
	      break;
	      //tt++;
	      //continue;
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

      if (selfoverlap||merr)
	{
	  tt++;
#if 0	  
	  save_conf_mc(tt, 0); 
#endif
	  continue;
	}
      if (overlap)// && ierr==0)
	totene += 1.0;
#if 0
      if (tt!=0 && tt%outits==0)
	save_conf_mc(tt, 0); 
#endif

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
  if (Oparams.parnum==0)
    return 0;
  ip = Oparams.parnum*ranf();
  /* qui basta calcolare l'energia della particella che sto muovendo */
  //return;
#if 0
  eno = calcpotene();
#else
  eno = calcpotene_GC(ip);
#endif
#ifdef ALIGN_POT
      /* correct eno here adding external field contrib */
#endif
  store_coord(ip);
  movetype=random_move(ip);

#ifdef MC_NEW_NNL_CHECK
  if (OprogStatus.useNNL && mc_check_NNL_escape(ip))
    {
      //printf("step # %d rebuilding NNL\n", Oparams.curStep);
      rebuildNNL();
    }
#endif
  pbc(ip);
  update_LL(ip);
 
  //rebuildLinkedList();
  //printf("i=%d\n", i);
  totmovesMC++;
  /* overlapMC() aggiorna anche i bond */
  //err=0;
#ifdef MC_FREEZE_BONDS
  //check_all_bonds();
#endif
  dorej = overlapMC(ip, &err);

#if 0
  if (OprogStatus.deltaMC==0 && OprogStatus.dthetaMC ==0)
    {
      if (dorej)
	printf("overlap!\n");
    }
#endif
  if (!dorej)
    {
#ifdef MC_STOREBONDS
#ifdef MC_FREEZE_BONDS
      if (!OprogStatus.freezebonds)
#endif
	store_bonds_mc(ip);
#endif
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds==1)
	{
	  refFB=0;
	  fakeFB=1;
	}
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      rejectMove = 0;
      checkMoveKF=1;
#endif
      update_bonds_MC(ip);

#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      checkMoveKF = 0;
#endif	
      /* qui basta calcolare l'energia della particella che sto muovendo */
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds==1)
	fakeFB=0;
#endif
#if 0
      enn=calcpotene();
#else
      enn=calcpotene_GC(ip);
#endif
#ifdef ALIGN_POT
      /* correct enn here adding external field contrib */
#endif
#if 0
      if (enn > eno) 
	printf("BOH rejectMove=%d\n", rejectMove);
#endif
#ifdef MC_FREEZE_BONDS
 //exit(-1);
      if (OprogStatus.freezebonds)
	{
	  if (refFB==1)
	    {
	      dorej=2;
	    }
	  else
	    dorej=0;
	}
      else
	{
     	  if (enn <= eno)
    	    {
    	      dorej=0;
#ifdef MC_NVE
	      if (OprogStatus.ensembleMC==4)
		OprogStatus.Ed += (enn-eno);
#endif
    	    }
	  else
	    {
#ifdef MC_NVE
	      if (OprogStatus.ensembleMC==4)
		{ 
		  if (OprogStatus.Ed >= (enn-eno))
		    {
		      OprogStatuts.Ed -= (enn-eno);
		      dorej = 0;
		    }
		  else 
		    dorej = 2;
		}
	      else
		{
		  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		    dorej=0;
		  else
		    dorej=2;
		}
#else
	      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		dorej=0;
	      else
		dorej=2;
#endif
	      // if (dorej==0)
	      // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	    }
	}
#else
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      if (rejectMove==1)
	{
	  dorej=2;
	}
      else
#endif
      if (enn <= eno)
	{
	  //printf("ene=%f (old=%f)\n", enn, eno);
	  //	  if (abs(enn-eno) >=1 )
	  //	    printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	  dorej=0;
#ifdef MC_NVE
	  if (OprogStatus.ensembleMC==4)
	    {
	      //printf("OprogStatus.Ed=%f step=%d\n", OprogStatus.Ed, Oparams.curStep);
	      OprogStatus.Ed += (eno-enn);
	    }
#endif
	}
      else
	{
#ifdef MC_NVE
	  if (OprogStatus.ensembleMC==4)
	    { 
	      if (OprogStatus.Ed >= (enn-eno))
		{
		  OprogStatus.Ed -= (enn-eno);
		  dorej = 0;
		}
	      else 
		dorej = 2;
	    }
	  else
	    {
	      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		dorej=0;
	      else
		dorej=2;
	    }
#else
	  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
	    dorej=0;
	  else
	    {
	      dorej=2;
	    }
#endif
	  // if (dorej==0)
	  // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	}
#endif
    }

#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds)
    refFB=0;
#endif
#if 0
#ifdef MC_KERN_FRENKEL
  if (dorej && abs(enn-eno) >= 1)
    {
      printf("ene=%f (old=%f)\n", enn, eno);
    }
#endif
#endif
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
#ifdef MC_FREEZE_BONDS
      if (dorej==2 && !OprogStatus.freezebonds)
#else
      if (dorej==2)
#endif
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
extern double totitsHC, numcallsHC;
extern double calc_phi(void);
#ifdef MC_CLUSTER_MOVE
void tra_cls_move(int nc)
{
  int ip, np;
  double dx, dy, dz;
  dx = OprogStatus.delTclsMC*(ranf()-0.5);
  dy = OprogStatus.delTclsMC*(ranf()-0.5);
  dz = OprogStatus.delTclsMC*(ranf()-0.5);
  //  printf("dr=%f %f %f\n", dx, dy, dz);
  for (np=0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np]; 
      rx[ip] += dx;
      ry[ip] += dy;
      rz[ip] += dz;
#ifdef MC_BOND_POS
      update_spot_pos(ip,  dx, dy, dz);
#endif
    }
  if (OprogStatus.useNNL)
    {
      displMC = sqrt(Sqr(dx)+Sqr(dy)+Sqr(dz))*1.001;
    }
  traclsmoveMC++; 
}
#if 0
int is_cls_percolating(int nc)
{
#ifdef MC_GAPDNA
  return 0;
#endif
  if (nbcls[nc]==clsdim[nc]*2)
    return 1;
  else 
    return 0;
}
#endif
void rot_cls_move(int nc, int flip)
{
#ifdef MC_BOND_POS
  int extraspots=0; 
#endif
   double theta, thetaSq, sinw, cosw;
  double ox, oy, oz, OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  double rhb[3], rhbn[3];
  int k1, k2, k3, ip, np;
  /* pick a random orientation */
  orient(&ox,&oy,&oz);
  /* pick a random rotation angle */
#ifdef MC_FLIP_MOVE
  if (flip==1)
    theta = acos(0.0); /* 90 degree rotation */
  else
    theta= OprogStatus.delRclsMC*(ranf()-0.5);
#else
  theta= OprogStatus.delRclsMC*(ranf()-0.5);
#endif
  if (OprogStatus.useNNL)
    {
      for (np=0; np < clsdim[nc]; np++)
	{ 
	  ip = clsarr[firstofcls[nc]+np];
	  displMC = abs(theta)*0.5001*maxax[ip];
	} 
    }
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
  for (np=0; np < clsdim[nc]; np++)
    {
      /* rotate all particles within cluster */
      ip = clsarr[firstofcls[nc]+np];
      rhb[0] = rx[ip]-clsCoM[0][nc];
      rhb[1] = ry[ip]-clsCoM[1][nc];
      rhb[2] = rz[ip]-clsCoM[2][nc];
      /* note that I+M is the orthogonal matrix associated to the rotation around axiz (ox, oy, oz)
	 (see my J. Comp. Phys. 2010 on cond-mat pag. 4)
       */
      for (k1 = 0; k1 < 3; k1++)
	{
	  rhbn[k1] = 0;
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      if (k1==k2)
		rhbn[k1] += (1.0+M[k2][k1])*rhb[k2];
	      else
		rhbn[k1] += M[k2][k1]*rhb[k2];
	    }
	}
       //    printf("r[%d]=%f %f %f rhb=%f %f %f\n", ip, rx[ip], ry[ip], rz[ip], rhbn[0], rhbn[1], rhbn[2]);
#if 1
      rx[ip] = rhbn[0] + clsCoM[0][nc];
      ry[ip] = rhbn[1] + clsCoM[1][nc];
      rz[ip] = rhbn[2] + clsCoM[2][nc];
#endif
#if 1
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
#if 0
      printf("{");
      for (k1 = 0; k1 < 3; k1++)
	{
	  printf("{");
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      printf("%.15G", M[k1][k2]);
	      if (k2 < 3 - 1)
		printf(", ");
	    }
	  printf("}");
      if (k1 < 3-1)
	printf(",\n");
	}
      printf("}\n");
#endif
#endif
#ifdef MC_BOND_POS
      rA[0] = rx[ip];
      rA[1] = ry[ip];
      rA[2] = rz[ip];
      BuildAtomPos(ip, rA, R[ip], ratA);
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
      for (k1 = 0; k1 < typesArr[typeOfPart[ip]].nspots + extraspots; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    bpos[k2][ip][k1] = ratA[k1+1][k2];
	}
#endif

    }
  rotclsmoveMC++;
}

int random_cls_move(int nc)
{
  double p;
  //printf("random move ip=%d\n", ip);
  p=ranf();
  if (OprogStatus.restrmove==1|| OprogStatus.restrmove==2
#ifdef MC_RESTR_MATRIX
      || OprogStatus.restrmove==3
#endif
     )
  p=0.0;
  if (p <= 0.5)
    {
      tra_cls_move(nc);
      return 0;
    }
  else
    {
#ifdef MC_FLIP_MOVE
      /* flip_prob must be between 0 and 1 
	 (in Blaak, Frenkel and Mulder, J. Chem. Phys., Vol. 110, (1999) 
	 they suggest 0.1=1/10) */
      if (OprogStatus.flip_prob > 0.0 && ranf()< OprogStatus.flip_prob)
	rot_cls_move(nc, 1);
      else
	rot_cls_move(nc, 0);
#else
      rot_cls_move(nc, 0);
#endif
      return 1;
    } 
}
void find_clsbonds_MC(int nc)
{
  int ip, np;
  for (np=0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      if (OprogStatus.useNNL)
	find_bonds_one_NLL(ip);
      else
	find_bonds_one(ip);
      if (clsNPT==2)
	break;
    }
}

void store_cls_coord(int nc)
{
#ifdef MC_BOND_POS
  int extraspots=0;
#endif
  int np, ip, k1, k2;
  for (np = 0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      vx[ip]=rx[ip];
      vy[ip]=ry[ip];
      vz[ip]=rz[ip];
      for (k1=0; k1<3; k1++)
	for (k2=0; k2<3; k2++)
	  RoldAll[ip][k1][k2]=R[ip][k1][k2];
#ifdef MC_BOND_POS
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < typesArr[typeOfPart[ip]].nspots + extraspots; k2++)
	  {
	    bposold[k1][ip][k2] = bpos[k1][ip][k2];
	  } 
#endif
    } 
}
void restore_cls_coord(int nc)
{
#ifdef MC_BOND_POS
  int extraspots=0;
#endif
  int np, ip, k1, k2;

  for (np = 0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      rx[ip]=vx[ip];
      ry[ip]=vy[ip];
      rz[ip]=vz[ip];
      for (k1=0; k1<3; k1++)
	for (k2=0; k2<3; k2++)
	  R[ip][k1][k2]=RoldAll[ip][k1][k2]; 
#ifdef MC_BOND_POS
#ifdef MC_BOUNDING_SPHERES
      extraspots = typesArr[typeOfPart[ip]].nspotsBS;
#endif
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < typesArr[typeOfPart[ip]].nspots + extraspots; k2++)
	  {
	    bpos[k1][ip][k2] = bposold[k1][ip][k2];
	  } 
#endif
    }
}
void store_particles_cls(int nc)
{
  char fileop2[512], fileop[512];
  FILE *bf;
  int i, na, np;
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
#ifdef MD_PATCHY_HE
  int nn;
#ifndef MD_SPOT_GLOBAL_ALLOC
  double ratA[NA][3], ratB[NA][3];
#endif
#endif
  double Drx, Dry, Drz, RCMx, RCMy, RCMz;
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G @ %.15G %.15G C[%s]\n";
  sprintf(fileop2 ,"StoreBump-cls-%d", nc);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
  R2u();

#if 1
  fprintf(bf, ".Vol: %f\n", Oparams.rcut*Oparams.rcut*Oparams.rcut);
#else
  fprintf(bf, ".Vol: %f\n", L[0]*L[1]*L[2]);
#endif
  MD_DEBUG(printf("[Store bump]: %.15G\n", Oparams.time));
  RCMx = clsCoM[0][nc];
  RCMy = clsCoM[1][nc];
  RCMz = clsCoM[2][nc];
  for (np = 0; np < clsdim[nc]; np++)
    {
      i = clsarr[firstofcls[nc]+np]; 
      for (kk=0; kk < 3; kk++)
	{
	  axi[kk] = typesArr[typeOfPart[i]].sax[kk];
	}
      fprintf(bf, tipodat2,rx[i]-RCMx, ry[i]-RCMy, rz[i]-RCMz, uxx[i], uxy[i], uxz[i],  axi[1], 2*axi[0], "red");
      //writeAllCor(bf);
#ifdef MD_PATCHY_HE
      rA[0] = rx[i]-RCMx;
      rA[1] = ry[i]-RCMy;
      rA[2] = rz[i]-RCMz;
      BuildAtomPos(i, rA, R[i], ratA);
#ifdef EDHE_FLEX
      for (nn = 1; nn < typesArr[typeOfPart[i]].nspots+1; nn++)
	fprintf(bf,"%.15f %.15f %.15f @ %.15G C[orange]\n", 
		ratA[nn][0], ratA[nn][1], ratA[nn][2], typesArr[typeOfPart[i]].spots[nn-1].sigma*0.5);
#else
      for (nn = 1; nn < ((i < Oparams.parnumA)?MD_STSPOTS_A+1:MD_STSPOTS_B+1); nn++)
	fprintf(bf,"%.15f %.15f %.15f @ %.15G C[orange]\n", 
		ratA[nn][0], ratA[nn][1], ratA[nn][2], Oparams.sigmaSticky*0.5);
#endif
#endif
    }
  fclose(bf);
}
#ifdef MC_GAPDNA
int cluster_is_not_a_dimer(int nc)
{
  int kk, i, par1, par2, cc;
  long long int jj, jj2, aa, bb;

  if (clsdim[nc] != 4)
    return 1;

  i = clsarr[firstofcls[nc]];
  par1 = -1;
  par2 = -1;
  /* find reversible bond of first particle */
  for (kk=0; kk < numbonds[i]; kk++)
    {
      jj = bonds[i][kk] / (NANA);
      jj2 = bonds[i][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (aa > 1)
	{
	  par1=jj;
	  break;
	}
    }
  /* se i è pari allora i+1 è l'altro cilindro del gapped DNA 
     altrimenti è i-1. E' così per costruzione della conf iniziale. */
  if (i % 2 == 0)
    i += 1;
  else
    i -= 1;

  for (kk=0; kk < numbonds[i]; kk++)
    {
      jj = bonds[i][kk] / (NANA);
      jj2 = bonds[i][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (aa > 1)
	{
	  par2=jj;
	  break;
	}
    }
  /* se una delle due particelle non è legata
     in maniera reversibile allora non abbiamo un dimero */
  if (par1 == -1 || par2 == -1)
    return 1;

  /* se il cluster è formato da 4 particelle e due legate covalentemente
     hanno entrambi legami reversibili allora non puo' che essere un dimero */
  return 0;
}
#endif
#ifdef MC_GAPDNA
int cls_move_bonds=0;
#endif
int inside_cluster_move=0;

int cluster_move(void)
{
  int ncls, percolating, nc, totbonds;
  int iold, k, i, np, ip, dorej, movetype, err, in0, in1, np_in_cls, i0, ret;
  double enn, eno;
  double lastrx, lastry, lastrz;
#ifdef MC_ALMARZA
  double omegaratio;
  assign_bonds_almarza();
#endif
  inside_cluster_move=1;
  build_clusters(&ncls, &percolating, 0);
  inside_cluster_move=0;
  /* pick randomly a cluster */
  nc = ranf()*ncls;
  /* discard move if the cluster is percolating */ 
  //printf("clsdim[%d]=%d\n", nc, clsdim[nc]);
 
  //printf("CLUSTER MOVE dim=%d\n", clsdim[nc]); 
#if 1 // FIX THIS!!!
  if (is_cls_percolating(nc))
    {
      return -1;
    }
#endif
  eno = calcpotene();

#if 0
#ifdef MC_GAPDNA
  if (cluster_is_not_a_dimer(nc))
    {
      /* NOTA: le cluster move le applico solo a dimeri! */
      return -1;
    }
#endif
#endif
  /* qui basta calcolare l'energia della particella che sto muovendo */
  store_cls_coord(nc);
  //printf("BEG eno=%f\n", eno);
  for (k=0; k < 3; k++)
    clsCoM[k][nc] = 0.0;
  /* if it is percolating all particles have 2 bonds, 
     hence we have to initialize i0 */
  i0 =  clsarr[firstofcls[nc]];
  for (np=0; np < clsdim[nc]; np++)
    {
      i =  clsarr[firstofcls[nc]+np];
      if ( numbonds[i]==1) 
	{
	  i0 = i;
	  //printf("qui i0=%d\n", i0);
	  break;
	}
    }

  np_in_cls=1;
  i =  bonds[i0][0] / (NANA);
  iold=i0;	
  clsCoM[0][nc] = rx[i0];
  clsCoM[1][nc] = ry[i0];
  clsCoM[2][nc] = rz[i0]; 
  /* if particles has no bond we have to check... */ 
  while (numbonds[i0]>0)
    {
      /* rebuild the cluster as a whole */
      lastrx = rx[iold];
      lastry = ry[iold];
      lastrz = rz[iold];
      rx[i] -= L[0]*rint((rx[i]-lastrx)/L[0]);
      ry[i] -= L[1]*rint((ry[i]-lastry)/L[1]);
      rz[i] -= L[2]*rint((rz[i]-lastrz)/L[2]); 
      clsCoM[0][nc] += rx[i];
      clsCoM[1][nc] += ry[i];
      clsCoM[2][nc] += rz[i];
      np_in_cls++;
#if 0
      if (numbonds[i] > 2)
	{
	  printf("more than 2 bonds (%d) for particle %d!\n", numbonds[i], i);
	}
#endif
      if (numbonds[i]==1 || np_in_cls == clsdim[nc])
	{
	  //printf("numbonds[i=%d]=%d numbonds[i0=%d]=%d", i, numbonds[i], i0, numbonds[i0]);
	  break;
	}	   
      /* qui assumo che le particelle siano bifunzionali */
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
#if 0
  if (np_in_cls != clsdim[nc])
    {
      printf("we got a problem np_in_cls=%d  but clsdim[%d]=%d\n", np_in_cls, nc, clsdim[nc]);
      printf("percolating=%d\n", percolating);
      exit(-1);
    }
#endif
  for (k=0; k < 3; k++)
    clsCoM[k][nc] /= clsdim[nc];
  movetype=random_cls_move(nc);
  
#if 0
for (np=0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      printf("PRIMAnumbonds[%d]=%d\n", ip, numbonds[ip]);
    }
#endif
  for (np=0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      pbc(ip);
      update_LL(ip);
    };
  //rebuildLinkedList();
#if 0
  ret=check_bonds_mc("boh");
  for (np=0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      printf("DOPOnumbonds[%d]=%d\n", ip, numbonds[ip]);
    } 
  if (ret)
    {
      store_particles_cls(nc);
      exit(-1);
    }
#endif
  //printf("i=%d\n", i);
  totclsmovesMC++;
  /* overlapMC() aggiorna anche i bond */
  //err=0;
  dorej=0;
  //printf("clsdim=%d\n", clsdim[nc]);
  for (np=0; np < clsdim[nc]; np++)
    {
      ip = clsarr[firstofcls[nc]+np];
      clsNPT=1;
      dorej = overlapMC(ip, &err);
      clsNPT=0;
      if (dorej!=0)
	break;
    }
  /* NOTA: se nel caso dei gapped DNA si muovono solo i dimeri
     che non hanno bond liberi per formare ulteriori legami, 
     il seguente check diventa inutile e si può commentare */
  if (!dorej)
    {
#if 0
      totbonds = 0;
      store_bonds_mc(-1);
      for (np=0; np < clsdim[nc]; np++)
	{
	  ip = clsarr[firstofcls[nc]+np];
	  totbonds += numbonds[ip];
	}
      for (i=0; i < Oparams.parnum; i++)
	numbonds[i] = 0;
#else
#ifndef MC_ALMARZA
      clsNPT=1;
#endif
#endif
      /* NOTA: nel caso di GAPDNA la seguente routine calcola 
	 i bond per tutte le particelle che appartengono al cluster,
	 così verifico se il numero di bond è cambiato.
	 A seguite di una cluster move si possono formare legami o anche rompere
	 se si è mosso un cluster percolante.
       */
      find_clsbonds_MC(nc);
#if 0
      cls_move_bonds=0;
      for (np=0; np < clsdim[nc]; np++)
	{
	  ip = clsarr[firstofcls[nc]+np];
 	  cls_move_bonds += numbonds[ip]; 
	}
      if (cls_move_bonds != totbonds)
	{
	  dorej=1;
	  //printf("qui cls_move_bonds=%d totbonds=%d\n", cls_move_bonds, totbonds);
	}
      restore_bonds_mc(-1);
#else
#ifndef MC_ALMARZA
      if (clsNPT==2)
	enn=eno-1;
      else
	enn=eno;
      clsNPT=0;
      if (enn != eno)
	{
	  dorej=1;
	}
#else
      inside_cluster_move = 1;
      omegaratio=calc_almarza_prob();
      inside_cluster_move = 0;
      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno))*omegaratio)
	{
	  dorej=0;
	}
      else
	{
	  dorej=2;
	}
#endif
#endif
    }
  if (dorej != 0)
    {
      /* move rejected */
      totclsrejMC++;
      if (movetype==0)
	traclsrejMC++;
      else 
	rotclsrejMC++;
      if(err)
	{
	  printf("[random_move] NR failed...I rejected this trial move...\n");
	  err=0;
	}
      restore_cls_coord(nc);
      //rebuildLinkedList();
      for (np=0; np < clsdim[nc]; np++)
	{
	  ip = clsarr[firstofcls[nc]+np];
	  update_LL(ip);
	}
    }
  else
    {
   //   printf("step=%d cluster move accettata\n", Oparams.curStep);
    }
  if (OprogStatus.useNNL && dorej==0 )
    {
      for (np = 0; np < clsdim[nc]; np++)
	{	  
	  ip =  clsarr[firstofcls[nc]+np];
	  overestimate_of_displ[ip] += displMC; 
	}
    }
  return movetype;
}
#endif
#ifdef MC_BIGROT_MOVE
/* da straight (i.e. fully folded or unfolded) a qualsiasi orientazione */
#ifdef MC_BIGROT_BIASED
/* da qualsiasi orientazione a fully folded o fully unfolded */
extern void orient_biased(double *omx, double *omy, double* omz, double refax[3], int type);
int bigrot_move_outin(int bmtype) // bmtype=1 => out->in ; bmtype=0 => in->out
{
  int np, kk, k1, k2, dorej, err, extraspots, j, i;
  double ox, oy, oz, dist, enn, eno, fact, cos0;
  double rc[3], Rl[3][3], refaxi[3], refaxj[3];
  double *spXYZ=NULL, costheta;//, oldbene;
  /* first pick a particle randomly */
  if (Oparams.parnum==0)
    return 0;
  np = ranf()*Oparams.parnum;
  eno = calcpotene_GC(np);
  //oldbene= calc_bending_energy(np);
  /* now set center of rotation */
  store_coord(np);
#ifdef MC_BOND_POS
  for (kk=0; kk < 3; kk++)
    {
      rc[kk] = bpos[kk][np][0];
      ratA[1][kk] = rc[kk];
    }
#else
  rA[0] = rx[np];
  rA[1] = ry[np];
  rA[2] = rz[np];
  BuildAtomPos(np, rA, R[np], ratA);
  for (kk=0; kk < 3; kk++)
    {
      rc[kk] = ratA[1][kk];
    }
#endif
  i=np;
  if ((i%OprogStatus.polylen)==0)
    j = i+1;
  else
    j = i-1;
  for (kk=0; kk < 3; kk++)
    {
      /* l'asse di riferimento è quello della molecola a cui è quello del cilindro 
       * a cui è legato in maniera permanente.*/
      refaxj[kk] = R[j][0][kk];
      refaxi[kk] = R[i][0][kk];
    }
#if 1
  if (bmtype==1 || bmtype==2)
    {
      /* se la mossa è out->in (o out -> out) e non è out non fare nulla */
      costheta = scalProd(refaxi, refaxj);
      if (fabs(costheta) > cos(M_PI*OprogStatus.bigrotTheta0/180.0))
	{
	  //printf("i=%d j=%d non è out!!!!!!!!!!!!!! costheta=%.15G (%f)\n", i, j, costheta,cos(M_PI*OprogStatus.bigrotTheta0/180.0));
	  return 0;
	}
    }
  else
    {
      /* se la mossa è in->out (o in->in) e non è in non fare nulla */
      costheta = scalProd(refaxi, refaxj);
      if (fabs(costheta) <= cos(M_PI*OprogStatus.bigrotTheta0/180.0))
	{
	  //printf("i=%d j=%d non è out!!!!!!!!!!!!!! costheta=%.15G (%f)\n", i, j, costheta,cos(M_PI*OprogStatus.bigrotTheta0/180.0));
	  return 0;
	}
    }
#endif
  totbigrotmovMC++;
  //if (bmtype==1 || bmtype==2)
  //printf("bmtype=%d ( 0=in -> out ; 2 = out->out 1 = out-> in; 3 = in -> in)\n", bmtype);
  if (bmtype==0 || bmtype==2)// 0=in -> out ; 2 = out->out
    orient_biased(&ox,&oy,&oz, refaxj, 1);
  else //if (bmtype==1 || bmtype==3) // 1 = out-> in; 3 = in -> in
    orient_biased(&ox,&oy,&oz, refaxj, 0);
  spXYZ = typesArr[typeOfPart[np]].spots[0].x;
  dist=calc_norm(spXYZ);
  rx[np] = rc[0] + ox*dist;
  ry[np] = rc[1] + oy*dist;
  rz[np] = rc[2] + oz*dist;
 
  ox = -ox;
  oy = -oy;
  oz = -oz;
  versor_to_R(ox, oy, oz, Rl);
  for (k1 = 0; k1 < 3; k1++)
   for (k2 = 0; k2 < 3; k2++)
     {
       R[np][k1][k2] = Rl[k1][k2];
     }
#ifdef MC_BOND_POS
  rA[0] = rx[np];
  rA[1] = ry[np];
  rA[2] = rz[np];
  BuildAtomPos(np, rA, R[np], ratA);
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[np]].nspotsBS;
#endif
  for (k1 = 0; k1 < typesArr[typeOfPart[np]].nspots + extraspots; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	bpos[k2][np][k1] = ratA[k1+1][k2];
    }
#endif

#ifdef MC_NEW_NNL_CHECK
  if (OprogStatus.useNNL && mc_check_NNL_escape(np))
    {
      //printf("step # %d rebuilding NNL\n", Oparams.curStep);
      rebuildNNL();
    }
#endif
  pbc(np);
  update_LL(np);
 
  dorej = overlapMC(np, &err);

  if (!dorej)
    {
#ifdef MC_STOREBONDS
#ifdef MC_FREEZE_BONDS
      if (!OprogStatus.freezebonds)
#endif
	store_bonds_mc(np);
#endif
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds==1)
	{
	  refFB=0;
	  fakeFB=1;
	}
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      rejectMove = 0;
      checkMoveKF=1;
#endif
      update_bonds_MC(np);

#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      checkMoveKF = 0;
#endif	
      /* qui basta calcolare l'energia della particella che sto muovendo */
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds==1)
	fakeFB=0;
#endif
#if 0
      enn=calcpotene();
#else
      enn=calcpotene_GC(np);
#endif
      cos0 =cos(M_PI*OprogStatus.bigrotTheta0/180.0);
      if (bmtype==1) // out->in
	{
	  fact = ((1.0-cos0)/cos0)*((1.0 - OprogStatus.bigrotbias)/OprogStatus.bigrotbias);
	  //printf("1 (out->in) fact =%.15G pi=%.15G cos0=%.15G exp=%.15G new bend ene=%.15G (old:%.15G)\n", fact, M_PI, cos0,exp(-(1.0/Oparams.T)*(enn-eno)), calc_bending_energy(np), oldbene);
	}
      else if (bmtype==0) // in -> out
	{
	  //cos0 =cos(M_PI*OprogStatus.bigrotTheta0/180.0);
	  fact = cos0*OprogStatus.bigrotbias/((1.0-cos0)*(1.0-OprogStatus.bigrotbias));
	  //printf("0 (in->out) fact =%.15G pi=%.15G cos0=%.15G exp=%.15G new bend energy=%.15G (old: %.15G)\n", fact, M_PI, cos0,exp(-(1.0/Oparams.T)*(enn-eno)),calc_bending_energy(np), oldbene);
	}
      else if (bmtype==3) // in->in
	{
	  fact = 1.0;
	}
      else if (bmtype==2) // out -> out
	{
	  fact = 1.0;
	}
#if 0
      if (enn > eno) 
	printf("BOH rejectMove=%d\n", rejectMove);
#endif
#ifdef MC_FREEZE_BONDS
 //exit(-1);
      if (OprogStatus.freezebonds)
	{
	  if (refFB==1)
	    {
	      dorej=2;
	    }
	  else
	    dorej=0;
	}
      else
	{
     	  if (enn <= eno)
    	    {
    	      dorej=0;
#ifdef MC_NVE
	      if (OprogStatus.ensembleMC==4)
		OprogStatus.Ed += (enn-eno);
#endif
    	    }
	  else
	    {
#ifdef MC_NVE
	      if (OprogStatus.ensembleMC==4)
		{ 
		  if (OprogStatus.Ed >= (enn-eno))
		    {
		      OprogStatuts.Ed -= (enn-eno);
		      dorej = 0;
		    }
		  else 
		    dorej = 2;
		}
	      else
		{
		  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		    dorej=0;
		  else
		    dorej=2;
		}
#else
	      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		dorej=0;
	      else
		dorej=2;
#endif
	      // if (dorej==0)
	      // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	    }
	}
#else

#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      if (rejectMove==1)
	{
	  dorej=2;
	}
      else
#endif
      if (fact >= 1.0 && enn <= eno)
	{
	  //printf("qui adso!!!\n");
	  //printf("ene=%f (old=%f)\n", enn, eno);
	  //	  if (abs(enn-eno) >=1 )
	  //	    printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	  dorej=0;
#ifdef MC_NVE
	  if (OprogStatus.ensembleMC==4)
	    {
	      //printf("OprogStatus.Ed=%f step=%d\n", OprogStatus.Ed, Oparams.curStep);
	      OprogStatus.Ed += (eno-enn);
	    }
#endif
	}
      else
	{
#ifdef MC_NVE
	  if (OprogStatus.ensembleMC==4)
	    { 
	      if (OprogStatus.Ed >= (enn-eno))
		{
		  OprogStatus.Ed -= (enn-eno);
		  dorej = 0;
		}
	      else 
		dorej = 2;
	    }
	  else
	    {
	      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		dorej=0;
	      else
		dorej=2;
	    }
#else
#if 0
	  if (log(ranf()) < log(fact)-(1.0/Oparams.T)*(enn-eno))
	    dorej=0;
	  else
	    {
	      dorej=2;
	    }
#else 
	  if (ranf() < fact*exp(-(1.0/Oparams.T)*(enn-eno)))
	    dorej=0;
	  else
	    {
	      dorej=2;
	    }
#endif
	  // if (dorej==0)
	  // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	}
#endif
    }

#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds)
    refFB=0;
#endif
#if 0
#ifdef MC_KERN_FRENKEL
  if (dorej && abs(enn-eno) >= 1)
    {
      printf("ene=%f (old=%f)\n", enn, eno);
    }
#endif
#endif
  if (dorej != 0)
    {
      totbigrotmovrejMC++;
      /* move rejected */
#if 0
      totrejMC++;
      if (movetype==0)
	trarejMC++;
      else 
	rotrejMC++;
#endif
      //printf("restoring coords\n");
      if(err)
	{
	  printf("[random_move] NR failed...I rejected this trial move...\n");
	  err=0;
	}
      restore_coord(np);
      //rebuildLinkedList();
      update_LL(np);
#ifdef MC_FREEZE_BONDS
      if (dorej==2 && !OprogStatus.freezebonds)
#else
      if (dorej==2)
#endif
	{
#ifdef MC_STOREBONDS
	  remove_allbonds_ij(np, -2);
	  restore_bonds_mc(np);
#else
	  update_bonds_MC(np);
#endif
	}
      //printf("restoring finished\n");
      //rebuildLinkedList();
#endif
    }
#if 0
  else
    {
      printf("BIGROT move accepted!\n");
    }
#endif
  if (OprogStatus.useNNL && dorej==0 )
    {
      overestimate_of_displ[np] += displMC; 
    }
  //printf("accettata!\n");
  if (bmtype==1)
    return 1;
  else
    return 2;
}
#endif
int bigrot_move(void)
{
  int np, kk, k1, k2, dorej, err, extraspots, j;
  double ox, oy, oz, dist, enn, eno;
  double rc[3], Rl[3][3];
  double *spXYZ=NULL;
  /* first pick a particle randomly */
  totbigrotmovMC++;
  if (Oparams.parnum==0)
    return 0;
  np = ranf()*Oparams.parnum;
  eno = calcpotene_GC(np);
  /* now set center of rotation */
  store_coord(np);
#ifdef MC_BOND_POS
  for (kk=0; kk < 3; kk++)
    {
      rc[kk] = bpos[kk][np][0];
      ratA[1][kk] = rc[kk];
    }
#else
  rA[0] = rx[np];
  rA[1] = ry[np];
  rA[2] = rz[np];
  BuildAtomPos(np, rA, R[np], ratA);
  for (kk=0; kk < 3; kk++)
    {
      rc[kk] = ratA[1][kk];
    }
#endif
  orient(&ox,&oy,&oz);
  spXYZ = typesArr[typeOfPart[np]].spots[0].x;
  dist=calc_norm(spXYZ);
  rx[np] = rc[0] + ox*dist;
  ry[np] = rc[1] + oy*dist;
  rz[np] = rc[2] + oz*dist;
 
  ox = -ox;
  oy = -oy;
  oz = -oz;
  versor_to_R(ox, oy, oz, Rl);
  for (k1 = 0; k1 < 3; k1++)
   for (k2 = 0; k2 < 3; k2++)
     {
       R[np][k1][k2] = Rl[k1][k2];
     }
#ifdef MC_BOND_POS
  rA[0] = rx[np];
  rA[1] = ry[np];
  rA[2] = rz[np];
  BuildAtomPos(np, rA, R[np], ratA);
#ifdef MC_BOUNDING_SPHERES
  extraspots = typesArr[typeOfPart[np]].nspotsBS;
#endif
  for (k1 = 0; k1 < typesArr[typeOfPart[np]].nspots + extraspots; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	bpos[k2][np][k1] = ratA[k1+1][k2];
    }
#endif

#ifdef MC_NEW_NNL_CHECK
  if (OprogStatus.useNNL && mc_check_NNL_escape(np))
    {
      //printf("step # %d rebuilding NNL\n", Oparams.curStep);
      rebuildNNL();
    }
#endif
  pbc(np);
  update_LL(np);
 
  dorej = overlapMC(np, &err);

  if (!dorej)
    {
#ifdef MC_STOREBONDS
#ifdef MC_FREEZE_BONDS
      if (!OprogStatus.freezebonds)
#endif
	store_bonds_mc(np);
#endif
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds==1)
	{
	  refFB=0;
	  fakeFB=1;
	}
#endif
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      rejectMove = 0;
      checkMoveKF=1;
#endif
      update_bonds_MC(np);

#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      checkMoveKF = 0;
#endif	
      /* qui basta calcolare l'energia della particella che sto muovendo */
#ifdef MC_FREEZE_BONDS
      if (OprogStatus.freezebonds==1)
	fakeFB=0;
#endif
#if 0
      enn=calcpotene();
#else
      enn=calcpotene_GC(np);
#endif
#if 0
      if (enn > eno) 
	printf("BOH rejectMove=%d\n", rejectMove);
#endif
#ifdef MC_FREEZE_BONDS
 //exit(-1);
      if (OprogStatus.freezebonds)
	{
	  if (refFB==1)
	    {
	      dorej=2;
	    }
	  else
	    dorej=0;
	}
      else
	{
     	  if (enn <= eno)
    	    {
    	      dorej=0;
#ifdef MC_NVE
	      if (OprogStatus.ensembleMC==4)
		OprogStatus.Ed += (enn-eno);
#endif
    	    }
	  else
	    {
#ifdef MC_NVE
	      if (OprogStatus.ensembleMC==4)
		{ 
		  if (OprogStatus.Ed >= (enn-eno))
		    {
		      OprogStatuts.Ed -= (enn-eno);
		      dorej = 0;
		    }
		  else 
		    dorej = 2;
		}
	      else
		{
		  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		    dorej=0;
		  else
		    dorej=2;
		}
#else
	      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		dorej=0;
	      else
		dorej=2;
#endif
	      // if (dorej==0)
	      // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	    }
	}
#else
#if defined(MC_KERN_FRENKEL) || defined(MC_GAPDNA)
      if (rejectMove==1)
	{
	  dorej=2;
	}
      else
#endif
      if (enn <= eno)
	{
	  //printf("ene=%f (old=%f)\n", enn, eno);
	  //	  if (abs(enn-eno) >=1 )
	  //	    printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	  dorej=0;
#ifdef MC_NVE
	  if (OprogStatus.ensembleMC==4)
	    {
	      //printf("OprogStatus.Ed=%f step=%d\n", OprogStatus.Ed, Oparams.curStep);
	      OprogStatus.Ed += (eno-enn);
	    }
#endif
	}
      else
	{
#ifdef MC_NVE
	  if (OprogStatus.ensembleMC==4)
	    { 
	      if (OprogStatus.Ed >= (enn-eno))
		{
		  OprogStatus.Ed -= (enn-eno);
		  dorej = 0;
		}
	      else 
		dorej = 2;
	    }
	  else
	    {
	      if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
		dorej=0;
	      else
		dorej=2;
	    }
#else
	  if (ranf() < exp(-(1.0/Oparams.T)*(enn-eno)))
	    dorej=0;
	  else
	    {
	      dorej=2;
	    }
#endif
	  // if (dorej==0)
	  // printf("accetto la mossa energetica enn-eno=%.15G\n", enn-eno);
	}
#endif
    }

#ifdef MC_FREEZE_BONDS
  if (OprogStatus.freezebonds)
    refFB=0;
#endif
#if 0
#ifdef MC_KERN_FRENKEL
  if (dorej && abs(enn-eno) >= 1)
    {
      printf("ene=%f (old=%f)\n", enn, eno);
    }
#endif
#endif
  if (dorej != 0)
    {
      totbigrotmovrejMC++;
      /* move rejected */
#if 0
      totrejMC++;
      if (movetype==0)
	trarejMC++;
      else 
	rotrejMC++;
#endif
      //printf("restoring coords\n");
      if(err)
	{
	  printf("[random_move] NR failed...I rejected this trial move...\n");
	  err=0;
	}
      restore_coord(np);
      //rebuildLinkedList();
      update_LL(np);
#ifdef MC_FREEZE_BONDS
      if (dorej==2 && !OprogStatus.freezebonds)
#else
      if (dorej==2)
#endif
	{
#ifdef MC_STOREBONDS
	  remove_allbonds_ij(np, -2);
	  restore_bonds_mc(np);
#else
	  update_bonds_MC(np);
#endif
	}
      //printf("restoring finished\n");
      //rebuildLinkedList();
    }
#if 0
  else
    {
      printf("BIGROT move accepted!\n");
    }
#endif
  if (OprogStatus.useNNL && dorej==0 )
    {
      overestimate_of_displ[np] += displMC; 
    }
  //printf("accettata!\n");
  return 3;
}
#endif
void set_ini_numcells(void)
{
  cxini=cellsx;
  cyini=cellsy;
  czini=cellsz;

#ifdef MC_BOND_POS
  cxiniBP=cellsxBP;
  cyiniBP=cellsyBP;
  cziniBP=cellszBP;
    
#endif
#ifdef MC_BOUNDING_SPHERES
  cxiniBS=cellsxBS;
  cyiniBS=cellsyBS;
  cziniBS=cellszBS;
#endif
}
void move(void)
{
  double acceptance, traaccept, ene, eno, rotaccept, volaccept=0.0, volfrac;
#ifdef MC_BIGROT_BIASED
  double pbr;
#endif
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
  
#if 0
    {
      int ncls, perc, nc, ii;
      FILE *f;
      check_bonds_mc("boh");
      build_clusters(&ncls, &perc);
      f=fopen("clusters.txt", "w+");
      fprintf(f,"ncls=%d perc=%d\n", ncls, perc);
      for (nc=0; nc < ncls; nc++)
	{
	  for (ii=0; ii < clsdim[nc]; ii++)
	    fprintf(f,"%d ", clsarr[firstofcls[nc]+ii]);
	  fprintf(f,"\n");
	}
      fclose(f);
      exit(-1);
    }
#endif
#if 0 || defined(DEBUG_HCMC)
    {int overlap, ierr; 

      overlap=0;
      ierr=0;
      dostorebump=1;
      //printf("checking overlaps\n");
  for (i=0; i < Oparams.parnum; i++)
    {
      if (overlapMC(i, &ierr))
	{ 
	  overlap=1;
	  printf("PRIMA overlap di %d\n", i);
	  exit(-1);
	}
    }
  //store_bump(0, 1);

  //printf("...done overlap=%d\n", overlap);
  dostorebump=0;
  if (overlap)
    {
      printf("MC step=%d\n", Oparams.curStep);
      exit(-1);
    }
    }
#endif

#ifndef MC_NEW_NNL_CHECK
  if (OprogStatus.useNNL && do_nnl_rebuild())
    {
      //printf("Step #%d Rebuilding NNL\n", Oparams.curStep);
      for (i=0; i < Oparams.parnum; i++)
	overestimate_of_displ[i]=0.0;
      rebuildNNL();
      OprogStatus.lastNNLrebuildMC = Oparams.curStep;
    }
#endif
#endif
  set_ini_numcells();
  if (OprogStatus.nvbbias > 0)
    deln=OprogStatus.nvbbias;
  else
    deln=0;
#ifdef MC_NVE
  if (OprogStatus.ensembleMC==4)
    deln=0;
#endif
  //printf("deln=%d\n", deln);
  if (OprogStatus.ensembleMC==0
#ifdef MC_NVE
      || OprogStatus.ensembleMC==4
#endif
      )
    ntot = Oparams.parnum+deln;
#ifdef MC_CLUSTER_NPT
  else if (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3)
    ntot = Oparams.parnum+1+deln;
#elif defined(MC_NPT_XYZ)
  else if (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4)
    {
      ntot = Oparams.parnum+1+deln;
     // printf("ensemble=%d ntot=%d\n", OprogStatus.ensembleMC, ntot);
    }
#else
  else if (OprogStatus.ensembleMC==1)
    ntot = Oparams.parnum+1+deln;
#endif
#ifdef MC_GRANDCAN
  else if (OprogStatus.ensembleMC==2)
    ntot = OprogStatus.npav+OprogStatus.nexc+deln;
#endif
  else
    ntot=Oparams.parnum;
  //check_all_bonds();
  for (i=0; i < ntot; i++)
    {
#ifdef MC_BIGROT_MOVE
      if (OprogStatus.bigrotmov > 0.0 && ranf() < OprogStatus.bigrotmov)
	{
#ifdef MC_BIGROT_BIASED
	  if (OprogStatus.bigrotTheta0 < 0.0)
    	    movetype = bigrot_move();
	  else
	    {
	      pbr = ranf();
	      //printf("pbr=%f\n", pbr);
	      /* with equal probability select a in->out/out->in or in->in or out->out move*/
	      if (pbr < (1.0/3.0))
		{
		  //movetype = bigrot_move_outin(0);
		  if (ranf() < OprogStatus.bigrotbias)
		    {
		      movetype = bigrot_move_outin(1);// out -> in (in = theta < theta0)
		    }
		  else
		    {
		      movetype = bigrot_move_outin(0); //in -> out
		    }
		}
#if 1
	      else if (pbr > (2.0/3.0)) // out -> out
		movetype = bigrot_move_outin(2);
	      else 
		movetype = bigrot_move_outin(3); // in -> in	
#endif	
	    }
#else
	  movetype = bigrot_move();
#endif
	  /* se la mossa bigrot non è stata effettuata allora tenta una mossa normale */
	  //if (movetype != 0)
	  continue;
	}
#endif
#ifdef MC_CLUSTER_MOVE
      if (OprogStatus.clsmovprob > 0.0 && ranf() < OprogStatus.clsmovprob)
	{
	  //printf("MCstep=%d i=%d\n", Oparams.curStep, i);
	  movetype=cluster_move();
	  continue;
	}
#endif
      if (OprogStatus.ensembleMC==0 && deln==0)
	ran = 0;
      else
	ran = ntot*ranf();
#if 0
      /* FOR DEBUG */
      ran = Oparams.parnum;
#endif
#ifdef MC_NVE
      if (OprogStatus.ensembleMC == 4)
	{
	  movetype = mcmotion();
	}	
#endif
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
      if ((OprogStatus.ensembleMC==1 
#ifdef MC_CLUSTER_NPT
	   || OprogStatus.ensembleMC==3
#elif defined(MC_NPT_XYZ)
	   || OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4
#endif
	   )
	  && ran>=Oparams.parnum)
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
#ifdef MC_CLUSTER_NPT
	      //check_bonds_mc("boh");
	      if  (OprogStatus.ensembleMC==3)
		move_box_cluster(&err);
	      else
		move_box(&err);
#elif defined(MC_NPT_XYZ)
	      if (OprogStatus.ensembleMC == 4)
		{
		  move_box_cluster_xyz(&err);
		}
	      else if (OprogStatus.ensembleMC==3)
		move_box_xyz(&err);
	      else 
		move_box(&err);
#else
	      move_box(&err);
#endif
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
	    {
	      movetype=mcmotion();
	    }
	}
      //printf("done\n");
    }
  if (OprogStatus.adjstepsMC < 0 || Oparams.curStep <= OprogStatus.adjstepsMC)
    {
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

#ifdef MC_CLUSTER_MOVE
	  acceptance=((double)(totclsmovesMC-totclsrejMC))/totclsmovesMC;
	  traaccept = ((double)(traclsmoveMC-traclsrejMC))/traclsmoveMC;
	  rotaccept = ((double)(rotmoveMC-rotrejMC))/rotmoveMC; 
	  if (traaccept > OprogStatus.targetAccept)
	    OprogStatus.delTclsMC *= 1.1;
	  else
	    OprogStatus.delTclsMC /= 1.1;
	  if (OprogStatus.restrmove==0)
	    {
	      if (rotaccept > OprogStatus.targetAccept)
		OprogStatus.delRclsMC *= 1.1;
	      else
		OprogStatus.delRclsMC /= 1.1;
	    }
#ifdef MD_LXYZ
	  if (OprogStatus.delTclsMC > (avL=pow(L[0]*L[1]*L[2],1.0/3.0))*0.1)
	    OprogStatus.delTclsMC = avL*0.1;
#else
	  if (OprogStatus.delTclsMC > L*0.1)
	    OprogStatus.delTclsMC = L*0.1;
#endif
	  if (OprogStatus.delRclsMC > 3.14)
	    OprogStatus.delRclsMC = 3.14;
#endif
	}
    }
  if (OprogStatus.adjstepsMC < 0 || Oparams.curStep <= OprogStatus.adjstepsMC)
    { 
#ifdef MC_CLUSTER_NPT
      if (OprogStatus.targetAcceptVol > 0.0 && (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3) && volmoveMC > 0 
	  && (Oparams.curStep % OprogStatus.resetacceptVol==0))
#elif defined(MC_NPT_XYZ)
      if (OprogStatus.targetAcceptVol > 0.0 && (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3 || OprogStatus.ensembleMC==4) && volmoveMC > 0 
          && (Oparams.curStep % OprogStatus.resetacceptVol==0))
#else
      if (OprogStatus.targetAcceptVol > 0.0 && OprogStatus.ensembleMC==1 && volmoveMC > 0     
  	  && (Oparams.curStep % OprogStatus.resetacceptVol==0))
#endif
	  {
	    volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
	    //printf("sono qui volaccept=%.15G\n", volaccept);
	    if (volaccept > OprogStatus.targetAcceptVol)
	      OprogStatus.vmax *= 1.1;
	    else
	      OprogStatus.vmax /= 1.1;
	  }
   }
  if (Oparams.curStep%OprogStatus.outMC==0)
    {
      if (OprogStatus.targetPhiMC > 0.0)
	printf("Current Phi=%.12G\n", calc_phi());
#ifdef MC_CLUSTER_NPT
      if ((OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3) && volmoveMC > 0)
	volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
#elif defined(MC_NPT_XYZ)
      if ((OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4) && volmoveMC > 0)
	volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
#else
      if (OprogStatus.ensembleMC==1 && volmoveMC > 0)
	volaccept = ((double)(volmoveMC-volrejMC))/volmoveMC;
#endif
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
#ifdef MC_CLUSTER_NPT
      if ((OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3) && volmoveMC > 0)
	printf("Volume moves acceptance = %.15G vmax = %.15G\n", volaccept, OprogStatus.vmax);
#elif defined(MC_NPT_XYZ)
      if ((OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4) && volmoveMC > 0)
	printf("Volume moves acceptance = %.15G vmax = %.15G\n", volaccept, OprogStatus.vmax);
#else
      if (OprogStatus.ensembleMC==1 && volmoveMC > 0)
	printf("Volume moves acceptance = %.15G vmax = %.15G\n", volaccept, OprogStatus.vmax);
#endif
#ifdef MC_CLUSTER_MOVE
      if (OprogStatus.clsmovprob > 0.0)
	{
	  if (totclsmovesMC > 0)
	    acceptance = ((double)(totclsmovesMC-totclsrejMC))/totclsmovesMC;
	  else
	    acceptance = 0;
	  if (traclsmoveMC > 0)
	    traaccept = ((double)(traclsmoveMC-traclsrejMC))/traclsmoveMC;
	  else 
	    traaccept = -1;
	  if (rotclsmoveMC > 0)
	    rotaccept = ((double)(rotclsmoveMC-rotclsrejMC))/rotclsmoveMC;
	  else 
	    rotaccept = -1;	  

	  printf("Cluster move acceptance: %.15G (tra=%.15G rot=%.15G) delTcls=%.12G delRcls=%.12G\n", acceptance, traaccept, rotaccept, OprogStatus.delTclsMC, OprogStatus.delRclsMC);
	}
#endif
#ifdef MC_BIGROT_MOVE
      if (OprogStatus.bigrotmov > 0.0)
	{
	  if (totbigrotmovMC > 0)
	    acceptance = ((double)(totbigrotmovMC - totbigrotmovrejMC))/totbigrotmovMC;
	  else
	    acceptance = 0;
	  printf("Bigrotmov move acceptance: %.15G (accettate=%lld)\n", acceptance, totbigrotmovMC - totbigrotmovrejMC);
	}

#endif
#if defined(MC_HC) && !defined(MC_HELIX)
      printf("Average iterations in case A.2:%G\n", totitsHC/numcallsHC);
#endif
    }
  if (OprogStatus.adjstepsMC < 0 || Oparams.curStep <= OprogStatus.adjstepsMC)
    {
#ifdef MC_CLUSTER_NPT
      if ((Oparams.curStep % OprogStatus.resetacceptVol == 0) && (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3))
	volmoveMC=volrejMC=0;
#elif defined(MC_NPT_XYZ)
      if ((Oparams.curStep % OprogStatus.resetacceptVol == 0) && (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4))
	volmoveMC=volrejMC=0;
#else
      if ((Oparams.curStep % OprogStatus.resetacceptVol == 0) && OprogStatus.ensembleMC==1)
	volmoveMC=volrejMC=0;
#endif
     if (Oparams.curStep % OprogStatus.resetaccept==0)
	{
#ifdef MC_CLUSTER_MOVE
    	  totclsmovesMC = totclsrejMC = 0;
	  traclsmoveMC = traclsrejMC = 0;
	  rotclsmoveMC = rotclsrejMC = 0;
#endif
#ifdef MC_BIGROT_MOVE
	  totbigrotmovMC=totbigrotmovrejMC=0;
#endif
 	  totmovesMC=totrejMC=0;
	  tramoveMC=trarejMC=0;
	  rotmoveMC=rotrejMC=0;
	}
    }
  if (Oparams.curStep==Oparams.totStep)
    {
      if (OprogStatus.targetAccept > 0.0)
	{
	  printf("deltrafin %.15G\n", OprogStatus.deltaMC);
	  printf("delrotfin %.15G\n", OprogStatus.dthetaMC);
	}
#ifdef MC_CLUSTER_NPT
      if (OprogStatus.targetAcceptVol > 0.0 && (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3))
	printf("delvolfin %.15G\n", OprogStatus.vmax);
#elif defined(MC_NPT_XYZ)
      if (OprogStatus.targetAcceptVol > 0.0 && (OprogStatus.ensembleMC==1||OprogStatus.ensembleMC==3||OprogStatus.ensembleMC==4))
	printf("delvolfin %.15G\n", OprogStatus.vmax);
#else
      if (OprogStatus.targetAcceptVol > 0.0 && OprogStatus.ensembleMC==1)
	printf("delvolfin %.15G\n", OprogStatus.vmax);
#endif
    }
#ifdef MC_SUS
  /* If susnmin == -1 and susnmax > 0 then susnmax is used to end the simulations: this is useful
     if we want to use GC MC to grow a system */
  if (OprogStatus.susnmin==-1 && OprogStatus.susnmax > 0 && Oparams.parnum >= OprogStatus.susnmax)
    ENDSIM=1;
#endif
#if 1
  if (OprogStatus.targetPhiMC > 0.0)
    { 
      volfrac = calc_phi();
      if (fabs((volfrac - OprogStatus.targetPhiMC)/OprogStatus.targetPhiMC) < OprogStatus.phitol)
	{
	  ENDSIM=1;
	} 
    }
#endif
}
#endif
