const double saxfactMC[3]={0.85,0.68,0.68};
#ifdef MC_SIMUL
#include<mdsimul.h>
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
/* MONTE CARLO CODE START HERE */
double rxold, ryold,rzold, Rold[3][3];
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
extern double ranf(void);
extern void orient(double *ox, double *oy, double *oz);
long long int rotmoveMC=0, tramoveMC=0, totmovesMC=0, trarejMC=0, rotrejMC=0, totrejMC=0, volrejMC=0, volmoveMC=0;
void tra_move(int ip)
{
  rx[ip]+=OprogStatus.deltaMC*ranf(); 
  ry[ip]+=OprogStatus.deltaMC*ranf();
  rz[ip]+=OprogStatus.deltaMC*ranf();
  tramoveMC++; 
}
void rot_move(int ip)
{
  double theta, thetaSq, sinw, cosw;
  double ox, oy, oz, OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random orientation */
  orient(&ox,&oy,&oz);
  /* pick a random rotation angle */
  theta= OprogStatus.dthetaMC*ranf();
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
extern void set_semiaxes_vb(double fx, double fy, double fz);
extern double calcDistNegNNLoverlapPlane(double t, double t1, int i, int j, double shift[3]);
extern double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);

double check_overlap(int i, int j, double shift[3], int *errchk)
{
  int k1, k2;
  double vecg[8], vecgNeg[8];
  double d, d0, r1[3], r2[3], alpha; 
  OprogStatus.optnnl = 0;

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
  set_semiaxes_vb(1.01*(typesArr[typeOfPart[0]].sax[0]),
		  1.01*(typesArr[typeOfPart[0]].sax[1]), 
		  1.01*(typesArr[typeOfPart[0]].sax[2]));
  
  d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
  /* se d0 è positiva vuol dire che i due parallelepipedi non s'intersecano */
  if (d0 > 0.0)
    {
      return 1.0;
    }

  set_semiaxes_vb(saxfactMC[0]*typesArr[typeOfPart[0]].sax[0],
		  saxfactMC[1]*typesArr[typeOfPart[0]].sax[1], 
		  saxfactMC[2]*typesArr[typeOfPart[0]].sax[2]);

  d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
  /* se d0 è positiva vuol dire che i due parallelepipedi non s'intersecano */
  if (d0 < 0.0)
    {
      return -1.0;
    }
  OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
  calcdist_retcheck = 0;
  
  d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
  *errchk = calcdist_retcheck;
  if (*errchk)
    {
      return -1.0;
    }
  //printf("QUI d=%f\n", d);
  return d;
}
int overlapMC(int ip, int *err)
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
void remove_from_current_cell(int i)
{
  int n;
  n = (inCell[2][i] * cellsy + inCell[1][i] )*cellsx + inCell[0][i]
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
    rx[ip] -= L[0];
  else if (rx[ip] < -L2[0])
    rx[ip] += L[0];
  if (ry[ip] > L2[1])
    ry[ip] -= L[1];
  else if (ry[ip] < -L2[1])
    ry[ip] += L[1];
  if (rz[ip] > L2[2])
    rz[ip] -= L[2];
  else if (rz[ip] < -L2[2])
    rz[ip] += L[2];
}
void update_bonds_MC(int ip);
extern void find_bonds_flex_all(void);
extern double calcpotene(void);
extern void rebuildLinkedList(void);
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
}
void move_box(int *ierr)
{
  int i, ii;
#if 1
  double vo, lnvn, vn, Lfact, enn, eno, arg;
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
#ifdef MD_LXYZ
  L[0] *= Lfact;
  L[1] *= Lfact;
  L[2] *= Lfact;
#else
  L *= Lfact;
#endif
  for (i=0; i < Oparams.parnum; i++)
    {
      rx[i] *= Lfact;
      ry[i] *= Lfact;
      rz[i] *= Lfact; 
      pbc(i);
    }
  update_numcells();
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
#else
	  L /= Lfact;
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
	  rebuildLinkedList();
	  return;
	}
    }
  /* update all bonds with new positions */
  for (i=0; i < Oparams.parnum; i++)
    numbonds[i] = 0;
  find_bonds_flex_all();
  enn = calcpotene();
  arg = -(1.0/Oparams.T)*(enn-eno)+Oparams.P*(vn-vo) - (Oparams.parnum+1)*log(vn/vo)*Oparams.T;
  if (ranf() > exp(arg))
    {
      /* move rejected restore old positions */
#ifdef MD_LXYZ
      L[0] /= Lfact;
      L[1] /= Lfact;
      L[2] /= Lfact;
#else
      L /= Lfact;
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
      rebuildLinkedList();
      /* restore all bonds*/
      for (i=0; i < Oparams.parnum; i++)
	numbonds[i] = 0;
      find_bonds_flex_all();
    }
#endif
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
  find_bonds_one(ip);
}
void move(void)
{
  double acceptance, traaccept, ene, eno, rotaccept, volaccept;
  int ran, movetype, i,ip, err, dorej, enn;
  /* 1 passo monte carlo = num. particelle tentativi di mossa */
  //printf("Doing MC step #%d\n", Oparams.curStep);
  for (i=0; i < Oparams.parnum; i++)
    {
      /* pick a particle at random */
      if (OprogStatus.ensembleMC==1)
	ran=(Oparams.parnum+1)*ranf();
      else 
	ran = 0;
      if (ran==Oparams.parnum)
	{
	  move_box(&err);
	  movetype=3; /* 0 = tra; 1 = rot 2 = tra and rot; 3 = box */
	  volmoveMC++;
	}
      else
	{
	  ip = Oparams.parnum*ranf();
	  eno = calcpotene();
	  store_coord(ip);
	  movetype=random_move(ip);
	  pbc(ip);
	  update_LL(ip);
	  //rebuildLinkedList();
	  //printf("i=%d\n", i);
	  totmovesMC++;
	  /* overlapMC() aggiorna anche i bond */
	  dorej = overlapMC(ip, &err);
	  if (!dorej)
	    {
	      update_bonds_MC(ip);
	      enn=calcpotene();
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
		  printf("NR failed...I rejected this trial move...\n");
		}
	      restore_coord(ip);
	      //rebuildLinkedList();
	      update_LL(ip);
	      if (dorej==2)
		update_bonds_MC(ip);
	      //printf("restoring finished\n");
	      //rebuildLinkedList();
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
      printf("MC Step #%d pressure=%f temperature=%f\n", Oparams.curStep, Oparams.P, Oparams.T);
      printf("Acceptance=%.15G (tra=%.15G rot=%.15G) deltaMC=%.15G dthetaMC=%.15G\n", acceptance, traaccept, 
	     rotaccept, OprogStatus.deltaMC, OprogStatus.dthetaMC);
      printf("rotmoveMC:%d rotrefMC: %d cells= %d %d %d\n", rotmoveMC, rotrejMC, cellsx, cellsy, cellsz);
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
