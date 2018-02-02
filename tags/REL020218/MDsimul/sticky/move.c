#include<mdsimul.h>
//#define MD_BIGDT_REBUILD
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG09(x) 
#define MD_DEBUG10(x)  
#define MD_DEBUG11(x) 
#define MD_DEBUG15(x) 
#define MD_DEBUG20(x) 
#define MD_DEBUG29(x) 
#define MD_DEBUG30(x) 
#define MD_DEBUG36(x)  
#define MD_NEGPAIRS
#define MD_NO_STRICT_CHECK
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
#ifdef MD_CALENDAR_HYBRID
extern int *linearLists;
extern int numevPQ, overevHQ, totevHQ;
#endif
#ifdef MD_ASYM_ITENS
double **Ia, **Ib, **invIa, **invIb;
#else
double Ia, Ib, invIa, invIb;
#endif
#ifdef MD_GRAVITY
extern double g2;
double Lz2;
#endif

double *lastcol;
void AdjustLastcol(void)
{
  int i;
  for (i=0; i <  Oparams.parnum; i++)
    lastcol[i] -= Oparams.time;
}

struct LastBumpS *lastbump;
extern double *axa, *axb, *axc;
extern int *scdone;
extern double *maxax;
/* Routines for LU decomposition from Numerical Recipe online */
void ludcmpR(double **a, int* indx, double* d, int n);
void lubksbR(double **a, int* indx, double *b, int n);
void InvMatrix(double **a, double **b, int NB);
double min(double a, double b);
double zbrent(double (*func)(double), double x1, double x2, double tol);
extern double invaSq[2], invbSq[2], invcSq[2];
double rxC, ryC, rzC;
extern int SolveLineq (double **a, double *x, int n); 
int calcdist_retcheck;
double rA[3], rB[3];
int polinterr, polinterrRyck;
int do_check_negpairs = 0;
int set_pbonds(int i, int j)
{
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    return MD_PBONDS_AA;
  else if (i >= Oparams.parnumA && j >= Oparams.parnumA)
    return MD_PBONDS_BB;
  else
    return MD_PBONDS_AB;
#elif defined(MD_AB41)
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    return MD_PBONDS_AA;
#if 0
  else if (i >= Oparams.parnumA && j >= Oparams.parnumA)
    return MD_PBONDS_BB;
#endif  
  else
    return MD_PBONDS_AB;
#else
return MD_PBONDS;
#endif
#else
return MD_PBONDS;
#endif
}
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
int mapbondsaAB[MD_PBONDS_AB]={1,1,2,2,3,3};
int mapbondsbAB[MD_PBONDS_AB]={1,2,1,2,1,2};
int mapbondsaAA[MD_PBONDS_AA]={1,1,2,2};
int mapbondsbAA[MD_PBONDS_AA]={1,2,1,2};
int mapbondsaBB[MD_PBONDS_BB]={1,1,1,2,2,2,3,3,3};
int mapbondsbBB[MD_PBONDS_BB]={1,2,3,1,2,3,1,2,3};
#elif defined(MD_AB41)
int mapbondsaAB[MD_PBONDS_AB]={1,2,3,4};
int mapbondsbAB[MD_PBONDS_AB]={1,1,1,1};
int mapbondsaAA[MD_PBONDS_AA]={1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
int mapbondsbAA[MD_PBONDS_AA]={1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4};
#else
int mapbondsaSiO[MD_PBONDS]={1,1,2,2,3,3,4,4};
int mapbondsbSiO[MD_PBONDS]={1,2,1,2,1,2,1,2};
#endif
int *mapbondsa;
int *mapbondsb;
extern int *crossevtodel;
#else
int mapbondsa[MD_PBONDS]={1,1,2,2,3,3,4,4};
int mapbondsb[MD_PBONDS]={3,4,3,4,1,2,1,2};
#endif
long long int itsF=0, timesF=0, itsS=0, timesS=0, numcoll=0;
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern double accngA, accngB;
void ScheduleEventBarr (int idA, int idB, int idata, int atb, int idcollcode, double tEvent);
double calcDistNeg(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, double dists[MD_PBONDS], int bondpair);
void comvel_brown (COORD_TYPE temp, COORD_TYPE *m);
void remove_bond(int na, int n, int a, int b);
void add_bond(int na, int n, int a, int b);
void writeAsciiPars(FILE* fs, struct pascii strutt[]);
void writeAllCor(FILE* fs);
void InitEventList (void);
int bound(int na, int n, int a, int b);
extern void check_all_bonds(void);
#ifdef MD_HSVISCO
void calcT(void);
#endif
void print_matrix(double **M, int n)
{
  int k1, k2;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}
void print_matrixArr(double M[3][3])
{
  int k1, k2;
  const int n = 3;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}
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
void distanza(int ia, int ib);
double pi, invL, L2, Vz;   
double W, K, T1xx, T1yy, T1zz,
       T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz, 
       Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB, DQxxOld, 
       DQyyOld, DQzzOld, DQxyOld, DQyzOld, DQzxOld, DQxxOldKin, 
       DQyyOldKin, DQzzOldKin, DQxxOldHS, DQyyOldHS, DQzzOldHS, DQxxOldST, DQyyOldST, DQzzOldST,
       PxxKin, PyyKin, PzzKin, PxxHS, PyyHS, PzzHS, PxxST, PyyST, PzzST;
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
double DrSq = 0.0; 
const double timbig = 1E12;
/* used by linked list routines */
double *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
int  **tree, cellRange[2*NDIM], initUcellx, initUcelly, initUcellz;
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
int *inCell[2][3], *cellList[4], cellsx[4], cellsy[4], cellsz[4];
#else
int *inCell[3], *cellList, cellsx, cellsy, cellsz;
#endif
int evIdA, evIdB, parnumB, parnumA, evIdD, evIdE;
int evIdC;
extern int *bondscache, *numbonds, **bonds;
void newtDist(double x[], int n, int *check, 
	      void (*vecfunc)(int, double [], double [], int, int, double []),
	      int iA, int iB, double shift[3]);
void newtDistNeg(double x[], int n, int *check, 
		 void (*vecfunc)(int, double [], double [], int, int, double []),
		 int iA, int iB, double shift[3]);
void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[], double xb2[], 
	   int *nb);
void newt(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3]);
void rebuildCalendar(void);
void R2u(void);
void store_bump(int i, int j);
#ifdef MD_ROTDIFF_MIS
double calc_norm(double *vec);
double scalProd(double *A, double *B);
void vectProdVec(double *A, double *B, double *C);
void build_orient_matrix(int i, double Rt[3][3])
{
  double norm, u1[3], u2[3], u3[3], wsz;
  int a;
#ifdef MD_THREESPOTS
  if (i < Oparams.parnumA)
    {
      for (a=0; a < 3; a++)
	{
  	  u3[a]=R[i][a][2];
	  u2[a]=R[i][a][1];
	}
      norm = calc_norm(u3);
      for (a=0; a < 3; a++)
	u3[a] /= norm;
      wsz = scalProd(u2, u3);
      for (a=0; a < 3; a++)
	u2[a] = u2[a]-u3[a]*wsz;
      norm=calc_norm(u2);
      for (a=0; a < 3; a++)
	u2[a] /= norm;
      vectProdVec(u2, u3, u1);
      for (a = 0; a < 3; a++)
	{
	  Rt[0][a] = u1[a];
	  Rt[1][a] = u2[a];
	  Rt[2][a] = u3[a];
	}
    }
  else
    {
      /* TODO: in tale caso ho tre vettori che vanno ortonormalizzati */
    }
#else
  /* TODO: in tale caso ho tre vettori che vanno ortonormalizzati */
  /* caso dell'acqua e della silica */
  for (a=0; a < 3; a++)
    {
      u3[a]=R[i][a][2];      
      u2[a]=R[i][a][1];
    }
  norm = calc_norm(u3);
  for (a=0; a < 3; a++)
    u3[a] /= norm;
  wsz = scalProd(u2, u3);
  for (a=0; a < 3; a++)
    u2[a] = u2[a]-u3[a]*wsz;
  norm=calc_norm(u2);
  for (a=0; a < 3; a++)
    u2[a] /= norm;
  vectProdVec(u2, u3, u1);
  for (a = 0; a < 3; a++)
    {
      Rt[0][a] = u1[a];
      Rt[1][a] = u2[a];
      Rt[2][a] = u3[a];
    }
#endif
}
void update_MSDrot(int i)
{
  double ti;
  double Rt[3][3];
  ti = Oparams.time - OprogStatus.lastcolltime[i];
  /* sumox, sumoy e sumoz sono gli integrali nel tempo delle componenti della velocità
   * angolare lungo gli assi dell'ellissoide */
  /* ERRORE: la matrice d'orientazione va costruita poiché nel presente codice 
     R[][] contiene la posizione di 3 spots (o meno)!!!  */
#if 1
  build_orient_matrix(i, Rt);	
  OprogStatus.sumox[i] += (wx[i]*Rt[0][0]+wy[i]*Rt[0][1]+wz[i]*Rt[0][2])*ti;
  OprogStatus.sumoy[i] += (wx[i]*Rt[1][0]+wy[i]*Rt[1][1]+wz[i]*Rt[1][2])*ti;
  OprogStatus.sumoz[i] += (wx[i]*Rt[2][0]+wy[i]*Rt[2][1]+wz[i]*Rt[2][2])*ti;
#else
  OprogStatus.sumox[i] += wx[i]*ti;
  OprogStatus.sumoy[i] += wy[i]*ti;
  OprogStatus.sumoz[i] += wz[i]*ti;
#endif
}
#endif
/* ========================== >>> scalCor <<< ============================= */
void scalCor(int Nm)
{ 
  int i;
  double DRx, DRy, DRz, invL;

  invL = 1.0 / L;
  /* Reduced particles to first box */
  for(i=0; i < Oparams.parnum; i++)
    {
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * rx[i]);
      DRy = - L * rint(invL * ry[i]);
      DRz = - L * rint(invL * rz[i]);
      rx[i] += DRx;
      ry[i] += DRy;
      rz[i] += DRz;
    }
}

double calc_norm(double *vec);

#if 0
double max_ax(int i)
{
  double ma;
  ma = 0;
  if (axa[i]>ma)
    ma = axa[i];
  if (axb[i]>ma)
    ma = axb[i];
  if (axc[i]>ma)
    ma = axc[i];
  return ma;
}
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
void rebuild_linked_list()
{
  double L2;
  int j, n, nl, iA, nc;
  L2 = 0.5 * L;
  for (nl = 0; nl < 4; nl++)
    {
      cellsx[nl] = L / Oparams.rcut[nl];
      cellsy[nl] = L / Oparams.rcut[nl];
#ifdef MD_GRAVITY
      cellsz[nl] = (Lz+OprogStatus.extraLz) / Oparams.rcut[nl];
#else
      cellsz[nl] = L / Oparams.rcut[nl];
#endif 
    }
  for (nl = 0; nl < 4; nl++)
    {
      for (j = 0; j < cellsx[nl]*cellsy[nl]*cellsz[nl] + Oparams.parnum; j++)
	cellList[nl][j] = -1;
    }

  for (nc = 0; nc < 2; nc++)
    {
      /* rebuild event calendar */
      for (n = 0; n < Oparams.parnum; n++)
	{
	  iA = (n < Oparams.parnumA)?0:1;
	  if (iA==0 && nc == 0)
	    nl = 0;
	  else if (iA == 1 && nc == 0)
	    nl = 1;
	  else if (iA == 0 && nc == 1)
	    nl = 2;
	  else
	    nl = 3;
	  inCell[nc][0][n] =  (rx[n] + L2) * cellsx[nl] / L;
	  inCell[nc][1][n] =  (ry[n] + L2) * cellsy[nl] / L;
#ifdef MD_GRAVITY
	  inCell[nc][2][n] =  (rz[n] + Lz2) * cellsz[nl] / (Lz+OprogStatus.extraLz);
#else
	  inCell[nc][2][n] =  (rz[n] + L2)  * cellsz[nl] / L;
#endif
	  j = (inCell[nc][2][n]*cellsy[nl] + inCell[nc][1][n])*cellsx[nl] + 
	    inCell[nc][0][n] + Oparams.parnum;
	  cellList[nl][n] = cellList[nl][j];
	  cellList[nl][j] = n;
	}
    }
}
#else
void rebuild_linked_list()
{
  double L2;
  int j, n;
  L2 = 0.5 * L;
  cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
  cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
  cellsz = L / Oparams.rcut;
#endif 
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList[j] = -1;

  /* rebuild event calendar */
  for (n = 0; n < Oparams.parnum; n++)
    {
      inCell[0][n] =  (rx[n] + L2) * cellsx / L;
      inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
      inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
      inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
#endif
      j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + Oparams.parnum;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}
#endif
#ifdef MD_GROWTH_CODE
double calcDistNegHS(double t, double t1, int i, int j, double shift[3], double *r1, double *r2)
{
  double ti, rAB[3], rABn[3], norm, rA[3], rB[3], sigma, distSq;
  double rad1, rad2;
  int kk;
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  distSq = 0;
  
  for (kk = 0; kk < 3; kk++)
    rAB[kk] = rA[kk] - rB[kk];
  norm = calc_norm(rAB);
  //printf("ti=%.15G rA=%f %f %f rB=%f %f %f\n", ti,rA[0], rA[1], rA[2], rB[0], rB[1], rB[2]);
  //printf("norm=%.15G rABn=%.15G %.15G %.15G sax=%.15G i=%d j=%d\n", norm, rAB[0],rAB[1],rAB[2], typesArr[typeOfPart[i]].sax[0],
  //i, j);
  for (kk=0; kk < 3; kk++)
    distSq += Sqr(rAB[kk]);
  for (kk = 0; kk < 3; kk++)
    rABn[kk] = rAB[kk]/norm;
  /* 04/05/2010: questa funzione viene chiamata solo durante la crescita quindi
     questo if è inutile, tuttavia non si sa mai che serva altrove... */
  if (OprogStatus.targetPhi > 0.0)
    {
      rad1 = axa[i];
      rad2 = axa[j];
    }
  else
    {
      if (i < Oparams.parnumA)
	rad1 = Oparams.sigma[0][0]*0.5;
      else
	rad1 = Oparams.sigma[1][1]*0.5;
      if (j < Oparams.parnumA)
	rad2 = Oparams.sigma[0][0]*0.5;
      else
	rad2 = Oparams.sigma[1][1]*0.5;
    }
  for (kk = 0; kk < 3; kk++)
    {
      r1[kk] = rA[kk] - rABn[kk]*rad1;
      r2[kk] = rB[kk] + rABn[kk]*rad2;
    }
#if 0
  for (kk = 0; kk < 3; kk++)
    rAB[kk] = r1[kk] - r2[kk];

  printf("rAB norm = %.15G\n", calc_norm(rAB));
#endif
  sigma = rad1 + rad2; 
  //printf("sigma=%.15G\n", sigma);
  return  sqrt(distSq) - sigma;
}
#ifdef MD_USE_SINGLE_LL
double get_min_dist (int na, int *jmin, double *rCmin, double *rDmin, double *shiftmin) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double distMin=1E10,dist,vecg[8], alpha, shift[3];
  /*double cells[NDIM];*/
  int kk;
  double r1[3], r2[3];
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
 /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  calcdist_retcheck = 0;
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
		  if (n!=na) 
		    {
		      dist = calcDistNegHS(Oparams.time, 0.0, na, n, shift, r1, r2);
		      if (calcdist_retcheck)
			continue;
#if 0
		      if ((na==125||na==15) && (n==15||n==125))
			printf("$$$$ dist: %.12G\n", dist);
#endif
		      if (*jmin == -1 || dist<distMin)
			{
			  distMin = dist;
			  for (kk = 0; kk < 3; kk++)
			    {
			      rCmin[kk] = r1[kk];
			      rDmin[kk] = r2[kk];
			      shiftmin[kk] = shift[kk];
			    }
			  *jmin = n;
			}
		    }
		} 
	    }
	}
    }
  return distMin;
}
#else
/* NOTA 06/05/2010: voidLL[0] = 0 se la lista AA o BB ha almeno una particelle
   mentre voidLL[1] = 0 se la lista AB contiene almeno una particella.
   Lascio tale array poiché potrebbe tornare anche utile, chissà... */
int voidLL[2];
double get_min_dist (int na, int *jmin, double *rCmin, double *rDmin, double *shiftmin) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  int iA, nl_ignore;
  double sigSq=0.0, dr[NDIM], dv[NDIM], shift[NDIM],  
	 b, d, t, tInt, vv, distSq, t1=0.0, t2=0.0, evtime=0, evtimeHC;
  int overlap, ac, bc, acHC, collCodeOld, nc;
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cells[NDIM];*/
  double dist, r1[3], r2[3], distMin=1E10;
  int collCode, kk, nl;
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
#if 0
  int evCode;
#endif
  MD_DEBUG29(printf("PredictEvent: %d,%d\n", na, nb));
  MD_DEBUG(calc_energy("PredEv"));
  /* NOTA: le linked list sono tre:
   *  0 = lista dell'interazione AA
   *  1 = lista dell'interazione BB
   *  2 = lista dell'interazione AB 
   *  inoltre inCell[ncel][dir][na] contiene 
   *  la cella a cui appartiene la particella A 
   *  Per ogni particella ci sono due insiemi di celle
   *  ncel = 0 celle relative all'interazione della particella na
   *  con particelle della stessa specie.
   *  ncel = 1 celle relative all'interazione della particella na 
   *  con particelle della stessa specie. */
  iA = (na<Oparams.parnumA)?0:1;
  nl_ignore = (na<Oparams.parnumA)?1:0;
  voidLL[0] = voidLL[1] = 0;
  for (nl=0; nl < 4; nl++)
    { 
      if (nl==nl_ignore || nl==iA+2)
	continue;


      if (nl < 2)
	{
	  nc = 0;
	}
      else
	{
	  nc = 1;
	}
      for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCell[nc][2][na] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellsz[nl] - 1;    
	      shift[2] = - L;
	    } 
	  else if (jZ == cellsz[nl]) 
	    {
	      jZ = 0;    
	      shift[2] = L;
	    }
	  for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	    {
	      jY = inCell[nc][1][na] + iY;    
	      shift[1] = 0.0;
	      if (jY == -1) 
		{
		  jY = cellsy[nl] - 1;    
		  shift[1] = -L;
		} 
	      else if (jY == cellsy[nl]) 
		{
		  jY = 0;    
		  shift[1] = L;
		}
	      for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
		{
		  jX = inCell[nc][0][na] + iX;    
		  shift[0] = 0.0;
		  if (jX == -1) 
		    {
		      jX = cellsx[nl] - 1;    
		      shift[0] = - L;
		    } 
		  else if (jX == cellsx[nl]) 
		    {
		      jX = 0;   
		      shift[0] = L;
		    }
		  n = (jZ *cellsy[nl] + jY) * cellsx[nl] + jX + Oparams.parnum;
		  for (n = cellList[nl][n]; n > -1; n = cellList[nl][n]) 
		    {
		      //printf("nl=%d cellList[nl=%d][n=%d]:%d na=%d\n", nl, nl, n, cellList[nl][n], na);
		      if (n != na) 
			{
			  dist = calcDistNegHS(Oparams.time, 0.0, na, n, shift, r1, r2);
			  if (*jmin == -1 || dist<distMin)
			    {
			      distMin = dist;
			      for (kk = 0; kk < 3; kk++)
				{
				  rCmin[kk] = r1[kk];
				  rDmin[kk] = r2[kk];
				  shiftmin[kk] = shift[kk];
				}
			      *jmin = n;
			    }
			  voidLL[nc] = 1; 
			}
		    } 
		}
	    }
	}
    }
  return distMin;

}
#endif

#ifdef MD_USE_SINGLE_LL
double rcutIni;
#else
double rcutIni[4];
#endif
double calc_phi(void)
{
  double N = 0;
  //const double pi = acos(0)*2;
  int i ;
  /* Anche se la miscela è non additiva la volume fraction è calcolata usando i diametri di ogni specie
     e ignorando il diametro per l'interazione AB */
  if (OprogStatus.targetPhi > 0.0)
    {
      for (i=0; i < Oparams.parnum; i++)
	{
	  N += axa[i]*axa[i]*axa[i];
	}
    }
  else
    {
      N  = Oparams.parnumA*pow(Oparams.sigma[0][0]*0.5,3);
      N += (Oparams.parnum-Oparams.parnumA)*pow(Oparams.sigma[1][1]*0.5,3);
    }
  N *= 4.0*pi/3.0;
#ifdef MD_LXYZ
  return N / (L[0]*L[1]*L[2]);
#else
  return N / (L*L*L);
#endif
}
double aL, bL, cL;
#ifdef MD_USE_SINGLE_LL
double rcutL;
#else
double rcutL[4];
#endif
void store_values(int i)
{
  int a;
#ifdef MD_USE_SINGLE_LL
  rcutL = Oparams.rcut;
#else
  for (a=0; a < 4; a++)
    rcutL[a] = Oparams.rcut[a];
#endif
  aL = axa[i];
#if 0
  bL = axb[i];
  cL = axc[i];
#endif
}
void restore_values(int i)
{
  int a;
#ifdef MD_USE_SINGLE_LL
  Oparams.rcut= rcutL;
#else
  for (a=0; a < 4; a++)
      Oparams.rcut[a] = rcutL[a];
#endif
  axa[i] = aL;
#if 0
  axb[i] = bL;
  axc[i] = cL;
#endif
}
void scale_coords(double sf)
{
  int i; 
#ifdef MD_LXYZ
  for (i=0; i < 3; i++)
    L[i] *= sf;
#else
  L *= sf;
#endif
  for (i = 0; i < Oparams.parnum; i++)
    {
      rx[i] *= sf;
      ry[i] *= sf;
      rz[i] *= sf;
    }
}
double max3(double a, double b, double c);
#if !defined(MD_USE_SINGLE_LL)
void calc_rcut(void)
{
  int i;
  double rc, mA=0.0, mB=0.0;
  /* rcut[2] = rcut[3] = max_i{axb[i]}+max_j{axb[i]}, 
     dove i è di tipo A e j di tipo B,
     mentre rcut[0] = max_i{axa[i]} con i di tipo A 
     rcut[1] = max_j{axa[j]} con j di tipo B
     chiaramente queste ultime sono sovrastime 
   */ 
  for (i=0; i < Oparams.parnumA; i++)
    {
      if (i==0)
	mA = axa[i];
      else
	{
	  if (axa[i] > mA)
	    mA = axa[i];
	}
    }

  if (mA*2.0 > Oparams.rcut[0])
    Oparams.rcut[0] = mA*2.0*1.01;
  
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      if (i==Oparams.parnumA)
	mB = axa[i];
      else
	{
	  if (axa[i] > mB)
	    mB = axa[i];
	}
    }

  if (mB*2.0 > Oparams.rcut[1])
    Oparams.rcut[1] = mB*2.0*1.01;

  if ((mA+mB) > Oparams.rcut[2])
    Oparams.rcut[2]=Oparams.rcut[3]=(mA+mB)*1.01;
}
#endif
double scale_axes(int i, double d, double rA[3], double rC[3], double rB[3], double rD[3], 
		      double shift[3], double scalfact, double *factor, int j)
{
  int kk;
  int ii, typei, a;
  double nnlfact, rmix;
  double C, Ccur, F, phi, fact, rAC[3], rBD[3], fact1, fact2;
  double boxdiag, factNNL=1.0;
  static double phi0;
  static int first=1;
#ifdef MD_LXYZ
  for (kk=0; kk < 3; kk++)
    L2[kk] = 0.5 * L[kk];
#else
  L2 = 0.5 * L;
#endif
#if !defined(MD_OPT_SCALEPHI)
  phi = calc_phi();
#endif
  if (first)
    {
      phi0 = ((double)Oparams.parnumA)*pow(Oparams.sigma[0][0]*0.5,3);
      phi0 +=((double)Oparams.parnum-Oparams.parnumA)*pow(Oparams.sigma[1][1]*0.5,3);
      phi0 *= 4.0*pi/3.0;
#ifdef MD_LXYZ
      phi0 /= L[0]*L[1]*L[2];
#else
      phi0 /= L*L*L;
#endif
      first = 0;
    }
  //printf("phi0=%.15G\n", phi0);
  C = cbrt(OprogStatus.targetPhi/phi0);

  if (i < Oparams.parnumA)
    Ccur = axa[i]*2.0/Oparams.sigma[0][0]; 
  else
    Ccur = axa[i]*2.0/Oparams.sigma[1][1]; 
  F = C / Ccur;
  if (j != -1)
    {
      for (kk=0; kk < 3; kk++)
	{
	  rAC[kk] = rA[kk] - rC[kk];
#ifdef MD_LXYZ
	  if (fabs (rAC[kk]) > L2[kk])
	    rAC[kk] -= SignR(L[kk], rAC[kk]);
#else
	  if (fabs (rAC[kk]) > L2)
	    rAC[kk] -= SignR(L, rAC[kk]);
#endif
	}
      for (kk=0; kk < 3; kk++)
	{
	  rBD[kk] = rB[kk] - rD[kk];
#ifdef MD_LXYZ
	  if (fabs (rBD[kk]) > L2[kk])
	    rBD[kk] -= SignR(L[kk], rBD[kk]);
#else
	  if (fabs (rBD[kk]) > L2)
	    rBD[kk] -= SignR(L, rBD[kk]);
#endif
	}
    }
  /* 0.99 serve per evitare che si tocchino */
  if (F < 1)
    fact = F;
  else
    {
     /*
	 if (OprogStatus.useNNL && maxax[i] / maxaxNNL > 1.0)
	 fact1 = 1.0 + scalfact*(maxaxNNL / maxaxStore - 1.0);
	 else
	 */	 
      if (j != -1)
	{
	  //printf("===> d=%.15G norm(rAC):%.15G scalfact=%.15G\n", d, calc_norm(rAC), scalfact);
	  fact1 = 1 + scalfact*(d / (calc_norm(rAC)));//+calc_norm(rBD)));
	}
      else
	{
#if 0
	  /* NOTA 29/04/2010: prima era così il che non aveva molto senso poiché
	     se non si usano le NNL non serve limitare fact1, 
	     comunque testare tale cambiamento. */
	    fact1 = 1.0 + 0.99*(factNNL - 1.0);
#else
	    /* prova il massimo se non va bene tanto poi riduce la crescita fino ad evitare overlap */
	    fact1 = F;
#endif
	}

      fact2 = F;
      if (fact2 < fact1)
	fact = fact2;
      else
	fact = fact1;
      //printf("i=%d j=%d fact=%.15G fact1=%.15G, fact2=%.15G factNNL=%.15G d=%.15G rAC=%.15G\n", i, j, fact, fact1, fact2, factNNL,d, calc_norm(rAC));

    }
 // printf("phi=%f fact1=%.8G fact2=%.8G scaling factor: %.8G\n", phi, fact1, fact2, fact);
  axa[i] *= fact;
  //axb[i] *= fact;
  *factor = fact;
#if !defined(MD_USE_SINGLE_LL)
  /* rcut[2] e rcut[3] sono le dimensioni delle celle per 
     le LL relative alle interazioni AB, per cui deve essere:
     rcut[2] > max_i{axa[i]}+max_i{axa[j]} (1)
     con i è una particella di tipo A e j di tipo B.
     Data la (1) allora se axa[i] > max_i{axa[i]} allora 
     il nuovo rcut sarà:
     rcut_new[2] = axa[i] + max_i{axa[j]}
   */
  calc_rcut();
#if 0
  Oparams.rcut[2] = rmix;
  Oparams.rcut[3] = rmix;
    
  if (i < Oparams.parnumA)
    {

      if (axa[i]*2.0 > Oparams.rcut[0])
	{
	  Oparams.rcut[0] = axa[i]*2.02;
	}
    }
  else 
    {
      if (axa[i]*2.0 > Oparams.rcut[1]) 
	{
	  Oparams.rcut[1] = axa[i]*2.0*1.01;
	}
    }
#endif
#else
  if (axa[i]*2.0 > Oparams.rcut)
    Oparams.rcut = 2.0*axa[i]*1.01;
#endif
#if !defined(MD_OPT_SCALEPHI)
  return calc_phi();
#else
  return 0.0;
#endif
}
int check_type_done(int t)
{
  int i, cc=0;
  if (t == 0)
    {
      for (i = 0; i < Oparams.parnumA; i++)
	{
	  if (scdone[i]==1)
	    cc++;
	}
      if (cc==Oparams.parnumA)
	return 1;
      else
	return 0;
    }
  else
    {
      for (i = Oparams.parnumA; i < Oparams.parnum; i++)
	{
	  if (scdone[i]==1)
	    cc++;
	}
      if (cc==Oparams.parnum-Oparams.parnumA)
	return 1;
      else
	return 0;
    }

}
double check_dist_min(int i, char *msg)
{
  int j;
  double distMin=1E60, rC[3], rD[3], shift[3];
  
  j = -1;
  distMin = get_min_dist(i, &j, rC, rD, shift);

  if (msg)
    {
      printf("[check_dist_min] %s distMin: %.12G i=%d j=%d\n", msg, distMin, i, j);
      //printf("rcut=%.15G %.15G %.15G\n", Oparams.rcut[0], Oparams.rcut[1], Oparams.rcut[2]);
    }
  return distMin;
} 
extern double *a0I;

double calc_safe_factor(void)
{
  int i, iMin;
  double axaMin;
  axaMin = axa[0];
  iMin = 0;
  for (i=1; i < Oparams.parnum; i++)
    {
      if (axa[i] < axaMin)
	{
	  axaMin = axa[i];
	  iMin = i;
	}
    }
  return a0I[iMin]/axaMin;
}

void scale_Phi(void)
{
#ifdef MD_SCALEPHI_STAGES
  static int curType=0;
#endif
  int a;
  int i, j, imin, kk, its, done=0;
  static int first = 1;
#if 0
  static double a0I;
#endif  
  static double target;
  double distMinT, distMin=1E60, rAmin[3], rBmin[3], rC[3]={0,0,0}, 
	 rD[3]={0,0,0};
  double shift[3], phi, scalfact, factor, axai;
  if (OprogStatus.targetPhi <= 0)
    return;
  phi=calc_phi();
  printf("Scaling Axes actual phi is %.15G\n", phi);
  if (first)
    {
      first = 0;
#if 1
      for (i = 0; i < Oparams.parnumA; i++)
	{
	  a0I[i] = Oparams.sigma[0][0]*0.5;
	}
      for (i = Oparams.parnumA; i < Oparams.parnum; i++)
	{
	  a0I[i] = Oparams.sigma[1][1]*0.5;
	}
#else
      a0I = Oparams.sigma[0][0]*0.5;
#endif
#if defined(MD_USE_SINGLE_LL)
      rcutIni = Oparams.rcut;
#else
      for (a=0; a < 4; a++)
       	rcutIni[a] = Oparams.rcut[a];
#endif
      target = cbrt(OprogStatus.targetPhi/calc_phi());
    }
  //UpdateSystem();   
#ifdef MD_LXYZ
  for (kk=0; kk < 3; kk++)
    L2[kk] = 0.5 * L[kk];
#else
  L2 = 0.5 * L;
#endif
  /* get the minimum distance in the system */
  //phi = calc_phi();
  for (kk = 0;  kk < 3; kk++)
    {
      cellRange[2*kk]   = - 1;
      cellRange[2*kk+1] =   1;
    }
  imin = -1;
#ifdef MD_SCALEPHI_STAGES
  if (OprogStatus.growthType == 1)
    {
      if (check_type_done(curType))
	{
	  curType++;
	}	
      printf("[GROWTH IN STAGES]: Actual growing type is: %d\n", curType);
    }
#endif
  for (i = 0; i < Oparams.parnum; i++)
    {
#ifdef MD_SCALEPHI_STAGES
      if (OprogStatus.growthType == 1)
	{
	  if (curType==0 && i >= Oparams.parnumA)
	    {
	      if (scdone[i]==1)
		done++;
	      continue;
	    }
	  else if (curType == 1 && i < Oparams.parnumA)
	    {
	      if (scdone[i]==1)
		done++;
	      continue;
	    }
	}
#endif
      j = -1;
      if (scdone[i]==1)
	{
	  done++;
	  continue;
	}
      /* NOTA 16/04/2010 */
      distMin = get_min_dist(i, &j, rC, rD, shift);
      //printf("distMin[%d]:%.15G\n", i, distMin);
      if (calcdist_retcheck)
	continue;
      if (j == -1)
	continue;
      //printf("i=%d j=%d distmin=%.10G\n", i, j, distMin);
      rAmin[0] = rx[i];
      rAmin[1] = ry[i];
      rAmin[2] = rz[i];
      if (j>-1)
	{
	  rBmin[0] = rx[j];
	  rBmin[1] = ry[j];
	  rBmin[2] = rz[j];
	}
      else
	{
	  rBmin[0] = 0.0;
	  rBmin[1] = 0.0;
	  rBmin[2] = 0.0;
	}
      scalfact = OprogStatus.scalfact;
      store_values(i);
      if (distMin < OprogStatus.minDist)// || fabs(distMin)<1E-10)//OprogStatus.epsd/10.0)
	continue;
      //printf("===> j=%d rA=(%f,%f,%f) rC=(%f,%f,%f)\n", j, rA[0], rA[1], rA[2], rC[0], rC[1], rC[2]);
      //printf("===> i=%d\n", i);
      phi = scale_axes(i, distMin, rAmin, rC, rBmin, rD, shift, scalfact, &factor, j);
      rebuild_linked_list();
      distMinT = check_dist_min(i, NULL);
      //printf("i=%d jmin=%d distMinT[%d]=%.15G\n", i, j, i, distMinT);
      //printf("axb[%d]=%.15g axb[%d]=%.15G\n", i, axb[i], j, axb[j]);
      if (calcdist_retcheck)
	{
	  restore_values(i);
	  rebuild_linked_list();
	  continue;
	}
      its = 0;
      while ( j != -1 && (distMinT < 0 || 
	     (distMinT > distMin && factor > 1.0) ||
	     (distMinT < distMin && factor < 1.0)) )
	{
	  restore_values(i);
	  phi = scale_axes(i, distMin, rAmin, rC, rBmin, rD, shift, scalfact, &factor, j);
  	  rebuild_linked_list();
  	  distMinT = check_dist_min(i, "Alla fine di calc_Phi()");
	  //printf("t=%.8G cellsx: %d rcut: %.8G imin=%d jmin=%d distMinT= %.15G phi=%.8G\n", 
	  //     Oparams.time, cellsx, Oparams.rcut, imin, jmin, distMinT, phi);
	  scalfact *= OprogStatus.reducefact;
	  its++;
	}

#if 0
      if (dist < distMin)
	{
	  distMin = dist;
	  jmin = j;
	  imin = i;
	  for (kk=0; kk < 3; kk++)
	    {
	      rCmin[kk] = rC[kk];
	      rDmin[kk] = rD[kk];
	      shiftmin[kk] = shift[kk];
	    }
	}
#endif
      if (i < Oparams.parnumA)
	axai = Oparams.sigma[0][0]*0.5;
      else
	axai = Oparams.sigma[1][1]*0.5;
#ifdef MD_OPT_SCALEPHI
      if (fabs(axa[i] / axai / target - 1.0) < OprogStatus.axestol)
#else
      if (fabs(axa[i] / axai - target) < OprogStatus.axestol)
#endif	
	{
	  done++;
	  scdone[i] = 1;
	  if (done == Oparams.parnum)
	    break;
	  continue;
	}
     
    }

  //check_alldist_min("DOPO");
#if 0
  if (phi > 0.7)
    {
      for (i=0; i < Oparams.parnum; i++)
	if (axa[i]>3.0||axb[i]>3.0||axc[i]>3.0)
	  printf("%d-(%f,%f,%f) ", i, axa[i], axb[i], axc[i]);
    }
#endif
#if defined(MD_OPT_SCALEPHI) 
  phi = calc_phi();
#endif
  printf("Scaled axes succesfully phi=%.8G\n", phi);
  //check_dist_min("PRIMA");
  rebuild_linked_list();
  rebuildCalendar();
  if (OprogStatus.intervalSum > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
#ifdef MD_DOUBLE_DT
  if (OprogStatus.brownian)
    ScheduleEvent(-1, ATOM_LIMIT+12,OprogStatus.nextDtR);
#endif

  printf("Scaled successfully %d/%d ellipsoids \n", done, Oparams.parnum);
#ifdef MD_OPT_SCALEPHI
  if (done == Oparams.parnum || fabs(phi / OprogStatus.targetPhi - 1.0)<OprogStatus.phitol)
#else
  if (done == Oparams.parnum || fabs(phi - OprogStatus.targetPhi)<OprogStatus.phitol)
#endif
    {
      printf("GROWTH DONE\n");
      R2u();
#if 0
      for (i=0; i < Oparams.parnum; i++)
	{
	  printf("i=%d axa=%.15G axb=%15G axc=%.15G\n", i, axa[i], axb[i], axc[i]);
	  distMin=check_dist_min(i, NULL);
	  if (distMin < 0.0)
	    {
	      printf("Prima distanza negativa per i=%d = %.15G\n", i, distMin);
	    }
	}
#endif
      ENDSIM = 1;
      /* riduce gli ellissoidi alle dimensioni iniziali e comprime il volume */
#if 1
      factor = calc_safe_factor();
#else
      factor = a0I/axa[0];
#endif
#ifdef MD_USE_SINGLE_LL
      Oparams.rcut = rcutIni;
#else
      for (a=0; a < 4; a++)
	Oparams.rcut[a] = rcutIni[a];
#endif
#if 0
      Oparams.a[0] *= factor;
      Oparams.b[0] *= factor;
      Oparams.c[0] *= factor;
      Oparams.a[1] *= factor;
      Oparams.b[1] *= factor;
      Oparams.c[1] *= factor;
#endif
      scale_coords(factor);
    }
#if 0
  for (i=0; i < Oparams.parnum; i++)
    {
      if (ENDSIM)
	OprogStatus.targetPhi = 0.0;
      distMin=check_dist_min(i, NULL);
      if (distMin < 0.0)
	{
	  printf("distanza negativa per i=%d = %.15G\n", i, distMin);
	}
    }
#endif
}

#endif
#ifdef MD_CALENDAR_HYBRID
void adjust_HQ_params(void);
#endif
void outputSummary(void)
{
  /* mettere qualcosa qui */
  if (timesS>0)
    printf("Average iterations in locate_contact: %.6G\n", ((double)itsS)/timesS);
  if (timesF>0)
    printf("Average iterations in search_contact_faster: %.6G\n",  ((double)itsF)/timesF);
  if (callsfrprmn>0)
    {
      printf("percentage of convergence in frprmn: %.6G\n", ((double) callsok)*100.0/callsfrprmn);
      printf("avg its in frprmn: %.10G\n", ((double) itsfrprmn)/callsfrprmn);
      printf("avg percentage of gradient: %.10G\n", (accngA+accngB)/callsfrprmn/2.0);
    }
  if (callsprojonto>0)
    printf("Average iterations in projonto: %.8G\n", ((double) itsprojonto)/callsprojonto);
  printf("Number of collisions: %lld\n", numcoll);
#ifdef MD_CALENDAR_HYBRID
#if 1
  if (OprogStatus.adjustHQ)
    {
      //printf("Average number of events in binary tree: %d (Actual is %d)\n", sumnumevPQ/callsAdjsuHQ, numevPQ);
      adjust_HQ_params();
    }
#endif
  printf("Actual number of events in binary tree: %d (of %d, overflows: %d)\n", numevPQ, totevHQ, overevHQ);
  //if (OprogStatus.curIndex > 0)
    //printf("lists[%d-1]=%d\n", OprogStatus.curIndex,linearLists[OprogStatus.curIndex-1]);
#endif
  if (OprogStatus.checkGrazing)
    check_all_bonds();
  //scale_Phi();
#ifdef MD_GROWTH_CODE
  if (!ENDSIM)
    scale_Phi();
#endif
#ifdef MD_GRAVITY
  printf("K= %.15f V=%.15f T=%.15f Vz: %f\n", K, V, 
	 (2.0*K/(3.0*Oparams.parnum-3.0)), Vz);
#else
#if 0
  printf("time=%.15f\n", OprogStatus.time);
  printf("K= %.15f T=%.15f\n", K, 
	 (2.0*K/(3.0*Oparams.parnum-3.0)));
#endif
#endif
#if 0
  f = fopenMPI(MD_HD_MIS "T.dat", "a");
  if (OprogStatus.brownian==1)
    fprintf(f, "%.15f %.15f\n", Oparams.time, (2.0*K/(3.0*Oparams.parnum)));
  else
    fprintf(f, "%.15f %.15f\n", Oparams.time, (2.0*K/(3.0*Oparams.parnum-3.0)));
  fclose(f);
#endif
}

void scalevels(double temp, double K)
{
  int i; 
  double sf;
#ifdef MD_SPOT_OFF
  sf = sqrt( ( (3.0*((double)Oparams.parnum)-3.0) * temp ) / (2.0*K) );
#else
#ifdef MD_THREESPOTS
  double dogTot, dogTra, dogRot;
  dogTra = 3.0*Oparams.parnum-3;
  dogRot = 2.0*Oparams.parnumA+3.0*(Oparams.parnum-Oparams.parnumA);
  dogTot = dogTra + dogRot; 
  sf = sqrt( dogTot * temp / (2.0*K) );
#else
  sf = sqrt( ( (6.0*((double)Oparams.parnum)-3.0) * temp ) / (2.0*K) );
#endif
#endif
  for (i = 0; i < Oparams.parnumA; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      wx[i] *= sf;
      wy[i] *= sf;
      wz[i] *= sf;
      /* scala anche i tempi di collisione! */
    } 
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      wx[i] *= sf;
      wy[i] *= sf;
      wz[i] *= sf;
      /* scala anche i tempi di collisione! */
    } 
  MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}
/* ============================ >>> updateQ <<< =========================== */
void updateDQ(COORD_TYPE dt)
{
  /* Here we use molecular pressure tensor */
  OprogStatus.DQxy += dt*Pxy;
  OprogStatus.DQyz += dt*Pyz;
  OprogStatus.DQzx += dt*Pzx;

  /*         _ t
	     |
	     DQab =  |  Pab dt
	     |
	     - 0

	     con a, b = x,y,z 
	     */
}

/* ============================= >>> updatePE <<< ========================= */
void updatePE(int Nm)
{
  int iE;
  double ENmin, ENmax;

  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  iE = (int) ((V - ((double) Nm)*ENmin) / 
	      ( ((double) Nm) * ((double) ENmax - ENmin)) * 
	      ((double) PE_POINTS));
  /*  sprintf(TXT, "energia: %f iE:%d\n", EE[my_rank], iE);
      mdPrintf(ALL, TXT, NULL);
      printf("PEPOINTS: %d Vc:%f\n", PE_POINTS, Vc);*/
  if ( (iE >= 0) && (iE < PE_POINTS)) 
    {
      ++(OprogStatus.PE[iE]);
    }

}

extern double WLJ;
#if 1
void check (int *overlap, double *K, double *V)
{
  *overlap = 0;
}
#else
void check (int *overlap, double *K, double *V)
{
  /* *******************************************************************
   ** TESTS FOR PAIR OVERLAPS AND CALCULATES KINETIC ENERGY.        **
   *******************************************************************
   */
  int  i, j;
  double rxi, ryi, rzi, rxij, ryij, rzij, rijSq, sigSq, rij;
  double tol=OprogStatus.overlaptol;

  printf("overlaptol: %f\n", OprogStatus.overlaptol);
  *overlap = 0;
  *V = 0.0;
  *K = 0.0;   
  for (i = 0; i < Oparams.parnum-1; i++)
    {     
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i];
      for (j = i+1; j < Oparams.parnum; j++)
	{
	  rxij = rxi - rx[j];
	  ryij = ryi - ry[j];
	  rzij = rzi - rz[j];
	  rxij = rxij - L*rint(invL*rxij);
	  ryij = ryij - L*rint(invL*ryij);
#if !defined(MD_GRAVITY)
	  rzij = rzij - L*rint(invL*rzij);
#endif
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  if (i < parnumA && j < parnumA)
	    sigSq = Sqr(Oparams.sigma[0][0]);
	  else if (i >= parnumA && j >= parnumA)
	    sigSq = Sqr(Oparams.sigma[1][1]);
	  else
	    sigSq = Sqr(Oparams.sigma[0][1]);
	  if ( rijSq < sigSq ) 
	    {
	      rij = sqrt(rijSq / sigSq);
	      if ( ( 1.0 - rij ) > tol ) 
		{
		  *overlap = 1;
#if 0
		  printf("#%d,#%d,rij/sigma =%f15.8 | ", i, j, rij);
		  printf("(%f,%f,%f)-(%f,%f,%f)\n", rxi, ryi, rzi, 
			 rx[j], ry[j], rz[j]);
#endif

		}
	    }
	}
    }
}
#endif
/* effettua la seguente moltiplicazione tra matrici:
 * 
 *                 | a 0 0 |               
 *  Trasposta(R) * | 0 b 0 | * R  dove R = {matrix unitaria relativa all'orientazione dell'ell.}
 *                 | 0 0 c |             
 * */
void tRDiagR(int i, double **M, double a, double b, double c, double **Ri)
{
  int na;
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  na = (i < Oparams.parnumA)?0:1;
  Di[0][0] = a;
  Di[1][1] = b;
  Di[2][2] = c;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  MD_DEBUG2(printf("a=%f b=%f c=%f Di=\n", a, b, c));
  MD_DEBUG2(print_matrixArr(Di));
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
  MD_DEBUG2(print_matrixArr(Rtmp));
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
void RDiagtR(int i, double **M, double a, double b, double c, double **Ri)
{
  int na;
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  na = (i < Oparams.parnumA)?0:1;
  Di[0][0] = a;
  Di[1][1] = b;
  Di[2][2] = c;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  MD_DEBUG2(printf("a=%f b=%f c=%f Di=\n", a, b, c));
  MD_DEBUG2(print_matrixArr(Di));
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k2][k3];
	  }
      }
  MD_DEBUG2(print_matrixArr(Rtmp));
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k1][k3]*Rtmp[k3][k2];
	  }
      }
}

void check_contact(int i, int j, double** Xa, double **Xb, double *rAC, double *rBC)
{
  int k1, k2;
  double f, g, del[3];
  f = g = -1; 
  for (k1=0; k1 < 3; k1++)
    del[0] = -L*rint((rAC[k1] - rBC[k1])/L);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	f += rAC[k1]*Xa[k1][k2]*rAC[k2];
	g += rBC[k1]*Xb[k1][k2]*rBC[k2];
      } 
  MD_DEBUG(printf("f(rC)=%.15G g(rC)=%.15G\n", f, g));
  //if (fabs(f) > 1E-5||fabs(g) > 1E-5)
  // exit(-1);
}
extern double **matrix(int n, int m);
extern void free_matrix(double **M, int n);
double calcDist(double t, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
		double *vecgsup, int calcguess);
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[]);
void core_bump(int i, int j, double *W, double sigSq)
{
  double rxij, ryij, rzij, factor;
  double delvx, delvy, delvz, invmi, invmj, denom;
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, taus, DTxx, DTyy, DTzz;
#endif

  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;
  denom = invmi + invmj; 
  
  rxij = rx[i] - rx[j];
  if (fabs (rxij) > L2)
    rxij = rxij - SignR(L, rxij);
  ryij = ry[i] - ry[j];
  if (fabs (ryij) > L2)
    ryij = ryij - SignR(L, ryij);
  rzij = rz[i] - rz[j];
#ifndef MD_GRAVITY
  if (fabs (rzij) > L2)
    rzij = rzij - SignR(L, rzij);
#endif
  factor = ( rxij * ( vx[i] - vx[j] ) +
	     ryij * ( vy[i] - vy[j] ) +
	     rzij * ( vz[i] - vz[j] ) ) / sigSq;
#if 0
  printf("prefcat=%.15G rxij * ( vx[i] - vx[j] ):%.15G %.15G %.15G\n", factor,rxij * ( vx[i] - vx[j] ) ,
       	    ryij * ( vy[i] - vy[j] ),  rzij * ( vz[i] - vz[j] ) );
#endif
  //if (!((Oparams.time - lastcol[i] < 1E-5)||
//	(Oparams.time - lastcol[j] < 1E-5)))
  //  factor *= (1+0.9)*0.5;

  factor *= 2.0 / denom;
#if 0
 if (factor > 0.00001)
   { 
     printf("v[i]=%f %f %f sig=%f\n", vx[i], vy[i], vz[i], sigSq);
     printf("v[j]=%f %f %f sig=%f\n", vx[j], vy[j], vz[j], sigSq);
     printf("rij=%f %f %f\n", rxij, ryij, rzij); 
     printf("factor=%.15G denom=%.15G\n", factor, denom);
   }
#endif
  /* Dissipation */
  delvx = - factor * rxij;
  delvy = - factor * ryij;
  delvz = - factor * rzij;
#ifdef MD_HSVISCO
  DTxy = delvx*delvy*invmi + vx[i]*delvy + delvx*vy[i];
  DTxy += delvx*delvy*invmj - vx[j]*delvy - delvx*vy[j]; 
  DTyz = delvy*delvz*invmi + vy[i]*delvz + delvy*vz[i];
  DTyz += delvy*delvz*invmj - vy[j]*delvz - delvy*vz[j];
  DTzx = delvz*delvx*invmi + vz[i]*delvx + delvz*vx[i];
  DTzx += delvz*delvx*invmj - vz[j]*delvx - delvz*vx[j];

  DTxx = delvx*delvx*invmi + vx[i]*delvx + delvx*vx[i];
  DTxx += delvx*delvx*invmj - vx[j]*delvx - delvx*vx[j]; 
  DTyy = delvy*delvy*invmi + vy[i]*delvy + delvy*vy[i];
  DTyy += delvy*delvy*invmj - vy[j]*delvy - delvy*vy[j];
  DTzz = delvz*delvz*invmi + vz[i]*delvz + delvz*vz[i];
  DTzz += delvz*delvz*invmj - vz[j]*delvz - delvz*vz[j];
#endif

  vx[i] = vx[i] + delvx*invmi;
  vx[j] = vx[j] - delvx*invmj;
  vy[i] = vy[i] + delvy*invmi;
  vy[j] = vy[j] - delvy*invmj;
  vz[i] = vz[i] + delvz*invmi;
  vz[j] = vz[j] - delvz*invmj;
#ifdef MD_HSVISCO 
  if (OprogStatus.lastcoll!=-1)
    {
      taus = Oparams.time - OprogStatus.lastcoll;
      OprogStatus.DQTxy += OprogStatus.Txy*taus; 
      OprogStatus.DQTyz += OprogStatus.Tyz*taus;
      OprogStatus.DQTzx += OprogStatus.Tzx*taus;

      OprogStatus.DQTxx += OprogStatus.Txx*taus; 
      OprogStatus.DQTyy += OprogStatus.Tyy*taus;
      OprogStatus.DQTzz += OprogStatus.Tzz*taus;
      //taus = Oparams.time - OprogStatus.lastcoll;
      //printf("DQT= %f %f %f\n", OprogStatus.DQTxy, OprogStatus.DQTyz, OprogStatus.DQTzx);
      OprogStatus.DQWxy += rxij*delvy;
      OprogStatus.DQWyz += ryij*delvz;
      OprogStatus.DQWzx += rzij*delvx;

      OprogStatus.DQWxx += rxij*delvx;
      OprogStatus.DQWyy += ryij*delvy;
      OprogStatus.DQWzz += rzij*delvz;
      OprogStatus.DQWxxHS += rxij*delvx;
      OprogStatus.DQWyyHS += ryij*delvy;
      OprogStatus.DQWzzHS += rzij*delvz;
      //printf("DQW= %f %f %f\n", OprogStatus.DQWxy, OprogStatus.DQWyz, OprogStatus.DQWzx);
    }
  OprogStatus.Txy += DTxy; 
  OprogStatus.Tyz += DTyz;
  OprogStatus.Tzx += DTzx;
  OprogStatus.Txx += DTxx; 
  OprogStatus.Tyy += DTyy;
  OprogStatus.Tzz += DTzz;
#endif

  /* TO CHECK: il viriale ha senso solo se non c'è la gravità */
  //*W = delvx * rxij + delvy * ryij + delvz * rzij;
}
double calcpotene(void);
void calc_energy(char *msg);
#ifdef MD_STORE_BONDS
void save_all_bonds(int i, int oldnb)
{
  FILE *fnb;
  char fname[128];
  int save_all, ii;
  sprintf(fname, "numbonds-%d", i);
  fnb = fopen(fname, "a");
  if (numbonds[i] == 0 || (numbonds[i] == 1 && oldnb == 0) )
    save_all = 1;
  else
    save_all = 0;
#ifdef MD_BIG_DT
  fprintf(fnb, "%.15G %d %d\n", Oparams.time+OprogStatus.refTime, numbonds[i], save_all);
#else
  fprintf(fnb, "%.15G %d %d\n", Oparams.time, numbonds[i], save_all);
#endif
  if (save_all)
    {
      UpdateSystem();
      for (ii = 0; ii < Oparams.parnum; ii++)
	{
	  fprintf(fnb, "%.15G %.15G %.15G\n", rx[ii]+L*OprogStatus.DR[ii][0], ry[ii]+L*OprogStatus.DR[ii][1],
		  rz[ii]+L*OprogStatus.DR[ii][2]);
	}
    }
  else
    {
      fprintf(fnb, "%.15G %.15G %.15G\n", rx[i]+L*OprogStatus.DR[i][0], ry[i]+L*OprogStatus.DR[i][1],
	      rz[i]+L*OprogStatus.DR[i][2]);
    }
  fclose(fnb);
}
#endif
void bump (int i, int j, int ata, int atb, double* W, int bt)
{
  /*
   *******************************************************************
   ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
   **                                                               **
   ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
   ** THE ROUTINE ALSO COMPUTES COLLISIONAL VIRIAL W.               **
   *******************************************************************
   */
  /* NOTA: Controllare che inizializzare factor a 0 è corretto! */
  double factor=0, invmi, invmj, sigmai, mredl;
  double delpx, delpy, delpz, wrx, wry, wrz, rACn[3], rBCn[3];
  double rAB[3], rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double ratA[3], ratB[3], norm[3];
  double bheight;
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, taus, DTxx, DTyy, DTzz;
#endif
#ifdef MD_STORE_BONDS
  int oldnbi, oldnbj;
#endif
  double denom, rCx, rCy, rCz, nrAB, Dr;
#ifndef MD_ASYM_ITENS
  double factorinvIa, factorinvIb;
#endif
  //double shift[3];
  double sigmaSticky;
  int na, a, kk;
#if 1
#ifdef MD_GROWTH_CODE
  if (OprogStatus.targetPhi > 0.0)
    {
      sigmai = axa[i]+axa[j];
    }
  else
    {
      if (i < Oparams.parnumA && j < Oparams.parnumA)
	{
	  /* qui si assume che ci possano essere due specie e che i diametri degli atomi
	   * componenti la molecola possano avere diametri diversi */ 
	  sigmai = Oparams.sigma[0][0];
	}
      else if (i >= parnumA &&  j >= Oparams.parnumA)
	{
	  sigmai = Oparams.sigma[1][1];
	}
      else
	{
	  sigmai = Oparams.sigma[0][1];
	}
    }
#else
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    {
      /* qui si assume che ci possano essere due specie e che i diametri degli atomi
       * componenti la molecola possano avere diametri diversi */ 
      sigmai = Oparams.sigma[0][0];
    }
  else if (i >= parnumA &&  j >= Oparams.parnumA)
    {
      sigmai = Oparams.sigma[1][1];
    }
  else
    {
      sigmai = Oparams.sigma[0][1];
    }
#endif
#endif
#ifdef MD_AB41
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    sigmaSticky = Oparams.sigmaStickyAA;
  else 
    sigmaSticky = Oparams.sigmaStickyAB;
#else
  sigmaSticky = Oparams.sigmaSticky;
#endif
  /*printf("(i:%d,j:%d sigSq:%f\n", i, j, sigSq);*/
  /*printf("mredl: %f\n", mredl);*/
  //MD_DEBUG(calc_energy("dentro bump1"));
  numcoll++;
  //printf("BUMP %d-%d btnone=%d  bthc=%d bt=%d time=%.15G\n", i, j,bt==MD_EVENT_NONE, bt==MD_CORE_BARRIER, bt, Oparams.time);
#if 1
#if 0
  calc_energy(NULL);
  printf("ene prima=%.10G enepot=%f K=%f\n", calcpotene()+K, calcpotene(), K);
#endif
  if (bt == MD_CORE_BARRIER)
    {
      /* do a normal collison between hard spheres 
       * (or whatever kind of collision between core objects
       * (ellipsoids as well...)*/
      MD_DEBUG36(printf("i=%d j=%d time=%.15G HARD CORE BUMP\n", i, j, Oparams.time));
      //printf("HC ene prima=%.10G enepot=%f K=%f\n", calcpotene()+K, calcpotene(), K);
      core_bump(i, j, W, Sqr(sigmai));
      //printf("HC ene dopo=%.10G enepot=%f K=%f\n", calcpotene()+K, calcpotene(), K);
      //check_all_bonds();
      MD_DEBUG10(printf(">>>>>>>>>>collCode: %d\n", bt));
      MD_DEBUG30(printf("time=%.15G collision type= %d %d-%d %d-%d ata=%d atb=%d\n",Oparams.time, bt, i, j, j, i, ata, atb));
#if 0
      calc_energy(NULL);
      printf("ene dopo=%.10G enepot=%f K=%f\n", calcpotene()+K, calcpotene(), K);
#endif 
      return;
    }
#endif
#if 0
  if ((i < Oparams.parnumA && j >= Oparams.parnumA) ||
      (i >= Oparams.parnumA && j < Oparams.parnumA))
    {
      printf("ok\n");
      exit(-1);
    }
#endif
  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  rB[0] = rx[j];
  rB[1] = ry[j];
  rB[2] = rz[j];
#if 1
  for (a=0; a < 3; a++)
    {
#ifdef MD_GRAVITY
      if (a==2)
       	continue;
#endif

      Dr = rA[a]-rB[a];
      if (fabs(Dr) > L2)
	{
	  
	  rB[a] += SignR(L, Dr);
//	  shift[a] =SignR(L, Dr);
//	  printf("in if L2: %f a=%d shift: %f Dr:%f\n",L2, a, SignR(L, Dr), Dr);
	}
  //   else
//	shift[a] = 0.0;
    }
#endif
  MD_DEBUG30(printf("[STICKY bump] t=%f contact point: %f,%f,%f \n", Oparams.time, rxC, ryC, rzC));
  /* qui calcolo il punto di contatto */
  MD_DEBUG30(printf("STICKY i=%d j=%d ata: %d atb: %d\n", i, j, ata, atb));
  MD_DEBUG20(printf("rA %f %f %f\n", rA[0], rA[1], rA[2]));
  MD_DEBUG20(printf("rB %f %f %f\n", rB[0], rB[1], rB[2]));
  BuildAtomPosAt(i, ata, rA, R[i], ratA);
  BuildAtomPosAt(j, atb, rB, R[j], ratB);
  MD_DEBUG20(printf("ratA %f %f %f\n", ratA[0], ratA[1], ratA[2]));
  MD_DEBUG20(printf("ratB %f %f %f\n", ratB[0], ratB[1], ratB[2]));
  for (kk = 0; kk < 3; kk++)
    rAB[kk] = ratA[kk] - ratB[kk];
  /* reduce to minimum image rAB[]!! */
#if 0
  calc_energy(NULL); 
  Kold = K;
  Vold = calcpotene();
  Eold = K + Vold;
  //printf("E PRIMA= %.15f\n", E);
#endif
  nrAB = calc_norm(rAB);
  MD_DEBUG(printf("sigmaSticky= %.15G norm rAB: %.15G\n", Oparams.sigmaSticky, calc_norm(rAB)));
#if 0
  if (fabs(nrAB - Oparams.sigmaSticky)> 1E-4)
    {
      printf("distance (%d,%d)-(%d,%d): %.15G)between sticky point is wrong!", 
	     i, j, ata, atb, nrAB);
      exit(-1);
    }
  for (kk=0; kk < 3; kk++)
    r12[kk] = rA[kk]-rB[kk];
  if (fabs(calc_norm(r12) - (Oparams.sigma[0][0]+Oparams.sigmaSticky))> Oparams.sigmaSticky)
    {
      printf("distance (%d,%d)-(%d,%d): %.15G sigma+sigsticky=%.15G between oxygens molecules is wrong!\n", 
	     i, j, ata, atb, calc_norm(r12), Oparams.sigma[0][0]+Oparams.sigmaSticky);
      for (kk=0; kk < 3; kk++)
	r12[kk] = ratA[kk]-rA[kk];
      printf("dist atA-A:%.15G\n", calc_norm(r12));
      for (kk=0; kk < 3; kk++)
	r12[kk] = ratB[kk]-rB[kk];
      printf("dist atB-B:%.15G\n", calc_norm(r12));
      for (kk=0; kk < 3; kk++)
	r12[kk] = ratA[kk]-ratB[kk];
      printf("dist atA-atB:%.15G\n", calc_norm(r12));
      exit(-1);
    }
	
#endif
  for (kk = 0; kk < 3; kk++)
    rAB[kk] /= nrAB;
  /* controllare con cura la scelta dei parametri relativi ai diametri delle sferette
   * e alle larghezze delle buche dei potenziali a buca quadrata */
  MD_DEBUG20(printf("coll code: %d\n", bt));

  rCx = ratA[0] - rAB[0]*sigmaSticky*0.5;
  rCy = ratA[1] - rAB[1]*sigmaSticky*0.5;
  rCz = ratA[2] - rAB[2]*sigmaSticky*0.5;

  rAC[0] = rA[0] - rCx;
  rAC[1] = rA[1] - rCy;
  rAC[2] = rA[2] - rCz;
 
  rBC[0] = rB[0] - rCx;
  rBC[1] = rB[1] - rCy;
  rBC[2] = rB[2] - rCz;

  //printf("contact point %f %f %f\n", rCx, rCy, rCz);
#if 0
  for (a=0; a < 3; a++)
    if (fabs (rAC[a]) > L2)
      rAC[a] -= SignR(L, rAC[a]);
#endif
#if 0
  r12[0] = ratA[0]-rCx;
  r12[1] = ratA[1]-rCy;
  r12[2] = ratA[2]-rCz;
  printf("ratAC: %.15G\n", calc_norm(r12));
  r12[0] = ratB[0]-rCx;
  r12[1] = ratB[1]-rCy;
  r12[2] = ratB[2]-rCz;
  printf("ratBC: %.15G\n",calc_norm(r12));
#endif
#if 0
    {
      int kk, amin, bmin;
      double Dr, dist, distSq, dists[MD_PBONDS], d1;
      /* calcola sigmaSq[][]!!! */
      distSq = 0;
      for (kk=0; kk < 3; kk++)
	distSq += Sqr(ratA[kk]-ratB[kk]);
      dist =  sqrt(distSq)-Oparams.sigmaSticky;
      //if (fabs(dist) > 1E-6)
      printf("[BUMP] %d-%d %d-%d type=%d ata=%d atb=%d t=%.15G dist= %.15G\n", i, j, j, i, bt,ata, atb, Oparams.time, dist);
      for (kk=0; kk < 3; kk++)
    	{
	  Dr = rA[kk]-rB[kk];
	  if (fabs(Dr) > L2)
	    shift[kk] = SignR(L, Dr);
	  printf("shift[%d]: %f\n", kk, shift);
	}
      printf("time-atomtime[%d]:%.15G time-atomtime[%d]:%.15G\n", i, Oparams.time - 
      atomTime[i], j, Oparams.time - atomTime[j]);
      d1 = calcDistNeg(Oparams.time, i, j, shift, &amin, &bmin, dists);
      for (kk=0; kk < MD_PBONDS; kk++)
  	printf("time=%.15G i=%d j=%d d1=%.15G dists[%d]:%.15G\n", Oparams.time, i, j, d1, kk, dists[kk]);
  }
#endif
#if 0
  for (a=0; a < 3; a++)
    {
      MD_DEBUG(printf("P rBC[%d]:%.15f ", a, rBC[a]));
      if (fabs (rBC[a]) > L2)
	rBC[a] -= SignR(L, rBC[a]);
      MD_DEBUG(printf("D rBC[%d]:%.15f ", a, rBC[a]));
    }
#endif
  MD_DEBUG(printf("\n"));
  /* calcola tensore d'inerzia e le matrici delle due quadriche */
  na = (i < Oparams.parnumA)?0:1;
#if 0
  if (OprogStatus.targetPhi > 0)
    {
      /* scalare tutti i raggi qui */
    }
#endif
#ifdef MD_ASYM_ITENS
  RDiagtR(i, Ia, ItensD[na][0], ItensD[na][1], ItensD[na][2], R[i]);
#else
  Ia = Oparams.I[na];
#endif
  na = (j < Oparams.parnumA)?0:1;
#if 0
  if (OprogStatus.targetPhi > 0)
    {
      /* scalare tutti i raggi qui */
    }
#endif
#ifdef MD_ASYM_ITENS
  RDiagtR(j, Ib, ItensD[na][0], ItensD[na][1], ItensD[na][2], R[j]);
#else
  Ib = Oparams.I[na];
#endif
  //MD_DEBUG(calc_energy("dentro bump2"));
  MD_DEBUG(check_contact(evIdA, evIdB, Xa, Xb, rAC, rBC));

  //MD_DEBUG(calc_energy("dentro bump3"));
  //calc_energy(NULL);
  //printf("ene prima=%.10G enepot=%f K=%f\n", calcpotene()+K, calcpotene(), K);
  /* calcola le matrici inverse del tensore d'inerzia */
#ifdef MD_ASYM_ITENS
  InvMatrix(Ia, invIa, 3);
  InvMatrix(Ib, invIb, 3);
#else
  invIa = 1/Ia;
  invIb = 1/Ib;
  MD_DEBUG20(printf("Ia=%f Ib=%f\n", Ia, Ib));
#endif
  for (a=0; a < 3; a++)
    norm[a] = rAB[a];

  MD_DEBUG(printf("CYL %f %f %f %f %f %f\n", rCx, rCy, rCz, norm[0], norm[1], norm[2]));
  /* calcola le velocità nel punto di contatto */
  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[0] = vx[i] + wrx;
  vCA[1] = vy[i] + wry;
  vCA[2] = vz[i] + wrz;
  vectProd(wx[j], wy[j], wz[j], -rBC[0], -rBC[1], -rBC[2], &wrx, &wry, &wrz);
  vCB[0] = vx[j] + wrx;
  vCB[1] = vy[j] + wry;
  vCB[2] = vz[j] + wrz;

  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;

  denom = invmi + invmj; 
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  MD_DEBUG20(printf("[bump] before bump vc=%.15G\n", vc));
#if 0
    {
      double sp=0;
      for (a=0; a < 3; a++)
	{
	  sp += -rAC[a]*norm[a];	  
	} 
      sp = 0;
      for (a=0; a < 3; a++)
	{
	  sp += -rBC[a]*norm2[a];	  
	} 
    }
#endif
#if 1
  if ((vc > 0 && bt == MD_OUTIN_BARRIER) ||
      (vc < 0 && bt == MD_INOUT_BARRIER))// && fabs(vc) > 1E-10)
    {
      MD_DEBUG(printf("norm = (%f,%f,%f)\n", norm[0], norm[1],norm[2]));
      MD_DEBUG(printf("vel  = (%f,%f,%f)\n", vx[i], vy[i], vz[i]));
      MD_DEBUG(printf("i=%d r = (%f,%f,%f)\n", i, rx[i], ry[i], rz[i]));
      printf("[WARNING] maybe second collision has been wrongly predicted\n");
      printf("relative velocity (vc=%.15G) at contact point is negative! I ignore this event...\n", vc);
      return;
    }
#endif
  vectProd(rAC[0], rAC[1], rAC[2], norm[0], norm[1], norm[2], &rACn[0], &rACn[1], &rACn[2]);
#ifdef MD_ASYM_ITENS 
  for (a=0; a < 3; a++)
    {
      rnI[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  rnI[a] += invIa[a][b]*rACn[b]; 
	}
    }
  for (a = 0; a < 3; a++)
    denom += rnI[a]*rACn[a];
#else
  for (a = 0; a < 3; a++)
    denom += invIa*Sqr(rACn[a]);
#endif
  vectProd(rBC[0], rBC[1], rBC[2], norm[0], norm[1], norm[2], &rBCn[0], &rBCn[1], &rBCn[2]);
#ifdef MD_ASYM_ITENS  
  for (a=0; a < 3; a++)
    {
      rnI[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  rnI[a] += invIb[a][b]*rBCn[b]; 
	}
    }
  for (a = 0; a < 3; a++)
    denom += rnI[a]*rBCn[a];
#else
  for (a = 0; a < 3; a++)
    denom += invIb*Sqr(rBCn[a]);
#endif
  /* Nel caso di gravita' e' intuile implementare il TC-model di Luding
   * per evitare il collasso inelastico.
   * Gli urti in tale caso sono tutti elastici. */ 
  /* SQUARE WELL: modify here */
  //factor =2.0*vc/denom;
  /* NOTA + o -????*/
#ifdef MD_STORE_BONDS
  oldnbi = numbonds[i];
  oldnbj = numbonds[j];
#endif
  mredl = 1.0 / denom;
#ifdef MD_AB41
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    bheight = Oparams.bheightAA;
  else
    bheight = Oparams.bheightAB;
#else
  bheight = Oparams.bheight;
#endif
  switch (bt)
    {
      /* N.B.
       * Notare che Oparams.bheight è la profondità della buca ed 
       * è una quantità positiva!!*/
      /* b = vc*r ma ora b -> vc riscrivere correttamente le seguenti equazioni!! */
    case MD_CORE_BARRIER:
      factor = -2.0*vc;
      factor *= mredl;
      break;
    case MD_INOUT_BARRIER:
#ifdef MD_INF_BARRIER
      factor = -2.0*vc;
#else
      if (Sqr(vc) < 2.0*bheight/mredl)
	{
#if 1
	  MD_DEBUG36(printf("t=%.15G vc=%.15G NOT ESCAPEING collType: %d d=%.15G\n", Oparams.time,
		 vc,  bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
#endif
	  factor = -2.0*vc;
	}
      else
	{
#if 1
	  MD_DEBUG36(printf("t=%.15G vc=%.15G ESCAPING collType: %d d=%.15G\n", Oparams.time, vc, bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
#endif
	  factor = -vc + sqrt(Sqr(vc) - 2.0*bheight/mredl);
	  //oldbond = bound(i, j, ata, atb);
	  //printf("qui\n");
	  remove_bond(i, j, ata, atb);
	  remove_bond(j, i, atb, ata);
#ifdef MD_STORE_BONDS
	  save_all_bonds(i, oldnbi);
	  save_all_bonds(j, oldnbj);
#endif
	}
#endif
      factor *= mredl;
#if 0
      if (fabs(distSq - sigDeltaSq)>1E-12)    
	printf("[bump]dist:%.20f\n",sqrt(distSq));
#endif
      break;
    case MD_OUTIN_BARRIER:
#ifdef MD_INF_BARRIER
      factor = -2.0*vc;
#else
#if 0
      if ((i==78 && j==98)|| (i==98 && j==78))
	printf("ENTERING 78-98 collType: %d\n", bt);
#endif
      add_bond(i, j, ata, atb);
      add_bond(j, i, atb, ata);
#ifdef MD_STORE_BONDS
      save_all_bonds(i, oldnbi);
      save_all_bonds(j, oldnbj);
#endif
      factor = -vc - sqrt(Sqr(vc) + 2.0*bheight/mredl);
#endif     
      MD_DEBUG36(printf("delta= %f height: %f mredl=%f\n", 
		      Sqr(vc) + 2.0*bheight/mredl, bheight, mredl));
      factor *= mredl;
      break;
    }

  MD_DEBUG(printf("factor=%f denom=%f\n", factor, denom));
  delpx = factor * norm[0];
  delpy = factor * norm[1];
  delpz = factor * norm[2];
#if 0
  ene= (Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])+
	Sqr(vx[j])+Sqr(vy[j])+Sqr(vz[j])); 
#endif
#ifdef MD_HSVISCO
  DTxy = delpx*delpy*invmi + vx[i]*delpy + delpx*vy[i];
  DTxy += delpx*delpy*invmj - vx[j]*delpy - delpx*vy[j]; 
  DTyz = delpy*delpz*invmi + vy[i]*delpz + delpy*vz[i];
  DTyz += delpy*delpz*invmj - vy[j]*delpz - delpy*vz[j];
  DTzx = delpz*delpx*invmi + vz[i]*delpx + delpz*vx[i];
  DTzx += delpz*delpx*invmj - vz[j]*delpx - delpz*vx[j];

  DTxx = delpx*delpx*invmi + vx[i]*delpx + delpx*vx[i];
  DTxx += delpx*delpx*invmj - vx[j]*delpx - delpx*vx[j]; 
  DTyy = delpy*delpy*invmi + vy[i]*delpy + delpy*vy[i];
  DTyy += delpy*delpy*invmj - vy[j]*delpy - delpy*vy[j];
  DTzz = delpz*delpz*invmi + vz[i]*delpz + delpz*vz[i];
  DTzz += delpz*delpz*invmj - vz[j]*delpz - delpz*vz[j];
#endif
  vx[i] = vx[i] + delpx*invmi;
  vx[j] = vx[j] - delpx*invmj;
  vy[i] = vy[i] + delpy*invmi;
  vy[j] = vy[j] - delpy*invmj;
  vz[i] = vz[i] + delpz*invmi;
  vz[j] = vz[j] - delpz*invmj;
#ifdef MD_HSVISCO 
  if (OprogStatus.lastcoll!=-1)
    {
      taus = Oparams.time - OprogStatus.lastcoll;
      OprogStatus.DQTxy += OprogStatus.Txy*taus; 
      OprogStatus.DQTyz += OprogStatus.Tyz*taus;
      OprogStatus.DQTzx += OprogStatus.Tzx*taus;

      OprogStatus.DQTxx += OprogStatus.Txx*taus; 
      OprogStatus.DQTyy += OprogStatus.Tyy*taus;
      OprogStatus.DQTzz += OprogStatus.Tzz*taus;
      //taus = Oparams.time - OprogStatus.lastcoll;
      //printf("DQT= %f %f %f\n", OprogStatus.DQTxy, OprogStatus.DQTyz, OprogStatus.DQTzx);
      OprogStatus.DQWxy += (rA[0]-rB[0])*delpy;
      OprogStatus.DQWyz += (rA[1]-rB[1])*delpz;
      OprogStatus.DQWzx += (rA[2]-rB[2])*delpx;

      OprogStatus.DQWxx += (rA[0]-rB[0])*delpx;
      OprogStatus.DQWyy += (rA[1]-rB[1])*delpy;
      OprogStatus.DQWzz += (rA[2]-rB[2])*delpz;
      OprogStatus.DQWxxST += (rA[0]-rB[0])*delpx;
      OprogStatus.DQWyyST += (rA[1]-rB[1])*delpy;
      OprogStatus.DQWzzST += (rA[2]-rB[2])*delpz;
      //printf("DQW= %f %f %f\n", OprogStatus.DQWxy, OprogStatus.DQWyz, OprogStatus.DQWzx);
    }
  OprogStatus.Txy += DTxy; 
  OprogStatus.Tyz += DTyz;
  OprogStatus.Tzx += DTzx;
  OprogStatus.Txx += DTxx; 
  OprogStatus.Tyy += DTyy;
  OprogStatus.Tzz += DTzz;
#endif

  MD_DEBUG(printf("delp=(%f,%f,%f)\n", delpx, delpy, delpz));
#ifdef MD_ROTDIFF_MIS
  update_MSDrot(i);
  update_MSDrot(j);
#endif
#ifdef MD_ASYM_ITENS
  for (a=0; a < 3; a++)
    {
      wx[i] += factor*invIa[0][a]*rACn[a];
      wx[j] -= factor*invIb[0][a]*rBCn[a];
      wy[i] += factor*invIa[1][a]*rACn[a];
      wy[j] -= factor*invIb[1][a]*rBCn[a];
      wz[i] += factor*invIa[2][a]*rACn[a];
      wz[j] -= factor*invIb[2][a]*rBCn[a];
    }
#else
  MD_DEBUG(printf(">>>>>>>>>>collCode: %d\n", bt));
  MD_DEBUG(printf("numbonds[0]:%d numbonds[1]:%d\n", numbonds[0], numbonds[1]));
  factorinvIa = -factor*invIa;
  factorinvIb = -factor*invIb;
  wx[i] += factorinvIa*rACn[0];
  wx[j] -= factorinvIb*rBCn[0];
  wy[i] += factorinvIa*rACn[1];
  wy[j] -= factorinvIb*rBCn[1];
  wz[i] += factorinvIa*rACn[2];
  wz[j] -= factorinvIb*rBCn[2];
#if 0
  if (isnan(wx[i]))
    {
      printf("factor: %f denom: %f vc: %f invIa: %f\n", factor, denom, vc, invIa);
      printf("bt: %d\n", bt);
      exit(-1);
    }
#endif
#endif
  //check_all_bonds();
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, vx[i],vy[i],vz[i], j, vx[j],vy[j],vz[j]));
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, wx[i],wy[i],wz[i], j, wx[j],wy[j],wz[j]));
  //printf("rC-ratA: %f\n", sqrt(Sqr(ratA[0]-rCx)+Sqr(ratA[1]-rCy)+Sqr(ratA[2]-rCz)));
  //printf("rC-ratB: %f\n", sqrt(Sqr(ratB[0]-rCx)+Sqr(ratB[1]-rCy)+Sqr(ratB[2]-rCz)));
  //calc_energy(NULL);
  //printf("ene dopo=%.10G enepot=%f K=%f\n", calcpotene()+K, calcpotene(), K);
#if 0
  calc_energy(NULL); 
  V = calcpotene();
  E = K + V;
  //printf("E DOPO= %.15f\n", E);
  if (fabs(E-Eold) > 1E-3)
    {
      printf("K-Kold: %f\n", K-Kold);
      printf("V-Vold: %f\n", V-Vold);
      printf("rC-ratA: %f\n", sqrt(Sqr(ratA[0]-rCx)+Sqr(ratA[1]-rCy)+Sqr(ratA[2]-rCz)));
      printf("rC-rA: %f\n",   sqrt(Sqr(rA[0]-rCx)+Sqr(rA[1]-rCy)+Sqr(rA[2]-rCz)));
      for (kk=0; kk < 3; kk++)
	r12[kk] = rA[kk]-rB[kk];
      printf("rAC: %f rBC=%f\n", calc_norm(rAC), calc_norm(rBC));
      printf("rAB: %f\n", calc_norm(rAB));
      printf("factor=%.15f\n", factor);
      printf("distance (%d,%d)-(%d,%d): %.15G sigma+sigsticky=%.15G between oxygens molecules is wrong!\n", 
	     i, ata, j, atb, calc_norm(r12), Oparams.sigma[0][0]+Oparams.sigmaSticky);
      for (kk=0; kk < 3; kk++)
	r12[kk] = ratA[kk]-rA[kk];
      printf("dist atA-A:%.15G\n", calc_norm(r12));
      for (kk=0; kk < 3; kk++)
	r12[kk] = ratB[kk]-rB[kk];
      printf("dist atB-B:%.15G\n", calc_norm(r12));
      for (kk=0; kk < 3; kk++)
	r12[kk] = ratA[kk]-ratB[kk];
      printf("dist atA-atB:%.15G\n", calc_norm(r12));
      printf("E-Eold= %.15G bt=%d\n", E-Eold,  bt);
      printf("bound:%d oldbond:%d\n", bound(i, j, ata, atb), oldbond);
      //exit(-1);
    }
#endif
 //check_bonds("CHECKING", i, j, ata, atb, 0);
  /* TO CHECK: il viriale ha senso solo se non c'è la gravità */
#if 0
  *W = delpx * rxij + delpy * ryij + delpz * rzij;
#endif
#if 0
  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[0] = vx[i] + wrx;
  vCA[1] = vy[i] + wry;
  vCA[2] = vz[i] + wrz;
  vectProd(wx[j], wy[j], wz[j], -rBC[0], -rBC[1], -rBC[2], &wrx, &wry, &wrz);
  vCB[0] = vx[j] + wrx;
  vCB[1] = vy[j] + wry;
  vCB[2] = vz[j] + wrz;

  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;

  printf("vcprima=%f\n", vc);
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  printf("vcdopo=%f\n", vc);
    {
      double d, shift[3], dists[MD_PBONDS];
      int i, amin, bmin;
      shift[0] = L*rint((rx[0]-rx[1])/L);
      shift[1] = L*rint((ry[0]-ry[1])/L);
      shift[2] = L*rint((rz[0]-rz[1])/L);
      d = calcDistNeg(Oparams.time+0.1, 0, 1, shift, &amin, &bmin, dists);
      for (i=0; i < MD_PBONDS; i++)
	{
	  printf("t=%.15G dists[%d]: %.15G\n",Oparams.time,i,dists[i]);
	}
    }
#endif
}

void calcObserv(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  double Drx, Dry, Drz;
  int i;
#ifdef MPI
  int equilib, sumeq;
#endif

  DrSqTot = 0.0;
  K = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      K += Oparams.m[0]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
    }
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      K += Oparams.m[1]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
    }
  K *= 0.5;
  for(i=0; i < Oparams.parnumA; i++)
    {
      Drx = rx[i] - OprogStatus.rxCMi[i] + L*OprogStatus.DR[i][0]; 
      Dry = ry[i] - OprogStatus.ryCMi[i] + L*OprogStatus.DR[i][1];
      Drz = rz[i] - OprogStatus.rzCMi[i] + L*OprogStatus.DR[i][2];
      DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
    }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
#ifdef MD_BIG_DT
  if (Oparams.time + OprogStatus.refTime > 0)
    Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time + OprogStatus.refTime) *
			 ((double) Oparams.parnumA ) );   
  else 
    Dtrans = 0;
#else
  if (Oparams.time>0)
    Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time) *
			 ((double) Oparams.parnumA ) );   
  else 
    Dtrans = 0;
#endif
  DrSqTot /= ((double) Oparams.parnumA);
  if (OprogStatus.eqlevel > 0.0)
    {
      if (!OprogStatus.equilibrated)
	OprogStatus.equilibrated = (DrSqTot>OprogStatus.eqlevel*(2.0*maxax[0])?1:0);
#ifdef MPI
      equilib = OprogStatus.equilibrated;
      MPI_Allgather(&equilib, 1, MPI_INT, 
		    equilibrated, 1, MPI_INT, MPI_COMM_WORLD);
      sumeq = 0;
      for (i=0; i < numOfProcs; i++)
	sumeq += equilibrated[i];
      /* se sumeq = numOfProcs vuol dire che tutti i processi sono
       * equilibrati quindi la simulazione può terminare */
      if (sumeq == numOfProcs)
	{
	  mdPrintf(ALL,"All systems reached equilibrium, simulation completed", 
		   NULL);
	  ENDSIM = 1;
	}
#else
      mdPrintf(ALL, "All systems reached equilibrium, simulation completed",
	       NULL);
      ENDSIM = OprogStatus.equilibrated;
#endif
    }
#if 0
  if (Oparams.time>0)
    {
      f = fopenMPI(MD_HD_MIS "D.dat", "a");
      fprintf(f, "%.15f %.15f\n", Oparams.time,  Dtrans);
      fclose(f);
    }
#endif
}
double *treeTime;
void adjust_norm(double **R)
{
  int k1, k2; 
  double n[3];
  for (k1 = 0; k1 < 3; k1++)
    {
      n[k1]=0;
      for(k2 = 0; k2 < 3; k2++)
	n[k1] += Sqr(R[k2][k1]);
      n[k1] = sqrt(n[k1]);
      if (fabs((n[k1])-1.0)>1E-10)
	{
	  MD_DEBUG(printf("Adjusting norm of orientations time=%.15f\n", Oparams.time));
	  MD_DEBUG(printf("delta = %.15f\n", fabs(n[k1]-1.0)));
	  for(k2 = 0; k2 < 3; k2++)
	    R[k2][k1] /= n[k1];
	}
    }
}

void UpdateAtom(int i)
{
  double ti;
  double wSq, w, sinw, cosw;
  double Omega[3][3], OmegaSq[3][3], Rtmp[3][3], M[3][3];
  int k1, k2, k3;
  ti = Oparams.time - atomTime[i];
  rx[i] += vx[i]*ti;
  ry[i] += vy[i]*ti;
#ifdef MD_GRAVITY
  rz[i] += vz[i]*ti - g2*Sqr(ti);
  vz[i] += -Oparams.ggrav*ti;
#else
  rz[i] += vz[i]*ti;
#endif
  /* ...and now orientations */
#if 0 && defined(MD_GRAVITY)
  if (i==24)
    {
      printf("[UPDATE ATOM] time=%.15G inCells[%d]=%d %d %d\n", Oparams.time, i, inCell[0][0][i],inCell[0][1][i],inCell[0][2][i]);
      printf("relpos=%.15G vz=%.15G\n", rz[i]+Lz2-0.5, vz[i]);
      printf("BAHBAH real cells: %d %d %d\n", (int)((rx[i] + L2) * cellsx[0] / L),
	     (int)((ry[i] + L2) * cellsy[0] / L), (int)((rz[i] + Lz2)  * cellsz[0] / (Lz+OprogStatus.extraLz)));
      printf("[END UPD ATOM]\n");
    }
#endif
 
#ifndef MD_SPOT_OFF
  wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
  w = sqrt(wSq);
  if (w != 0.0) 
    {
#if 0
      if (fabs(w*ti) < 1E-8)
	{
	  sinw = ti*(1-Sqr(w*ti)/6.0);	  
	  cosw = Sqr(ti)*(0.5 - Sqr(w*ti)/24.0);
	}
#endif
      sinw = sin(w*ti)/w;
      cosw = (1.0 - cos(w*ti))/wSq;
      Omega[0][0] = 0;
      Omega[0][1] = -wz[i];
      Omega[0][2] = wy[i];
      Omega[1][0] = wz[i];
      Omega[1][1] = 0;
      Omega[1][2] = -wx[i];
      Omega[2][0] = -wy[i];
      Omega[2][1] = wx[i];
      Omega[2][2] = 0;
      OmegaSq[0][0] = -Sqr(wy[i]) - Sqr(wz[i]);
      OmegaSq[0][1] = wx[i]*wy[i];
      OmegaSq[0][2] = wx[i]*wz[i];
      OmegaSq[1][0] = wx[i]*wy[i];
      OmegaSq[1][1] = -Sqr(wx[i]) - Sqr(wz[i]);
      OmegaSq[1][2] = wy[i]*wz[i];
      OmegaSq[2][0] = wx[i]*wz[i];
      OmegaSq[2][1] = wy[i]*wz[i];
      OmegaSq[2][2] = -Sqr(wx[i]) - Sqr(wy[i]);

      for (k1 = 0; k1 < 3; k1++)
	{

	  for (k2 = 0; k2 < 3; k2++)
	    {
	      //Omega[k1][k2] = -Omega[k1][k2];
	      Rtmp[k1][k2] = R[i][k1][k2];
	      M[k1][k2] = sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
#if 0
	      if (k1==k2)
		M[k1][k1] += 1.0;
#endif
	    }
	}
#if 0
#if MD_DEBUG(x)==x
      w1[0] = wx[i];
      w1[1] = wy[i];
      w1[2] = wz[i];
      for (k1 = 0; k1 < 3; k1++)
	{
	  w2[k1] = 0;
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      w2[k1] += (sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2])*w1[k2];   
	    }	 

	}
      if (fabs(w2[0])>1E-12 || fabs(w2[1])>1E-12
	  || fabs(w2[2])>1E-12)
	printf("ti=%.15G w=%.15G cosw=%.15G sinw=%.15G w = (%.15f,%.15f,%.15f) Exp(Omega t) w = (%.15f,%.15f,%.15f)\n",
	       ti, w, cosw, sinw, w1[0],w1[1],w1[2],w2[0],w2[1],w2[2]);	
#endif
#endif
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
#if 0
	    R[i][k1][k2] = 0.0;
#endif
	    for (k3 = 0; k3 < 3; k3++)
	      //  R[i][k1][k2] += M[k1][k3]*Rtmp[k3][k2];
	      R[i][k1][k2] += M[k1][k3]*Rtmp[k3][k2];
	  }
      //adjust_norm(R[i]);
#if 0
      if (rz[i]+Lz*0.5-Oparams.sigma/2.0 < 0. && OprogStatus.quenchend > 0.0)
	{
	  int no;
	  printf("rz[i](t-ti):%.15f rz[i]:%.15f ti:%.15f", 
		 rz[i] - vz[i]*ti + g2*Sqr(ti), rz[i],
		 ti);
	  printf("vz[i]:%.15f\n", vz[i]);
	  printf("i=%d SOTTO MURO PT rz+Lz*0.5=%.30f\n", i, rz[i]+Lz*0.5-Oparams.sigma/2.0);
	  printf("*******************************************************************\n");
	  if (vz[i] < 0 && evIdA != i && evIdB != ATOM_LIMIT +4)
	    {

	      exit(-1);
	    }
	}
#endif
    }
#endif
  atomTime[i] = Oparams.time;
}
void UpdateSystem(void)
{
  int i;
  /* porta tutte le particelle allo stesso tempo */
  for (i=0; i < Oparams.parnum; i++)
    {
      UpdateAtom(i);
#if 0
      if (ry[i] > L2 || ry[i] < -L2)
	{
	  printf("Porca Troia!!!\n");
	  exit(-1);
	}
#endif
    }
}
void check_bonds(char* msg, int i, int j, int ata, int atb, int yesexit)
{
  int a, b, B1;
  for (a = 0; a < numbonds[i]-1; a++)
    {
      B1 = bonds[i][a];
      
      for (b = a+1;  b < numbonds[i]; b++)
	{
	  if (B1 == bonds[i][b])
	    {
	      printf("Due bond uguali!!\n");
	      printf("bond=%d\n", B1);
	      printf("[%s] i=%d j=%d ata=%d atb=%d\n", msg, i, j, ata, atb);
	      if (yesexit)
		exit(-1);
	    }
	}
    }
}
void remove_bond(int na, int n, int a, int b)
{
  int i, nb, ii, jj, aa, bb, jj2;
  nb = numbonds[na];
  if (!nb)
    return;
  ii = 0;
  memcpy(bondscache, bonds[na], sizeof(int)*numbonds[na]);
  /* bonds[i] = j*(NA*NA) + a * NA + b 
   * dove b è l'atomo di j */
  for (i = 0; i < nb; i++)
    {
      jj = bondscache[i] / (NA*NA);
      jj2 = bondscache[i] % (NA*NA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (jj != n || aa != a || bb != b)
	{
	  bonds[na][ii++] = bondscache[i];
	} 
      else
	numbonds[na]--;
    }
  
  if (nb==numbonds[na])
    printf("nessun bond rimosso fra %d,%d\n", n, na);
#if 0
  if (abs(nb - numbonds[na])==2)
    printf("rimossi due bond boh...\n");
#endif
}

int bound(int na, int n, int a, int b);
void add_bond(int na, int n, int a, int b)
{
  if (bound(na, n, a, b))
    {
      printf("il bond %d,%d eiste già!\n", na, n);
      return;
    }
  bonds[na][numbonds[na]] = n*(NA*NA)+a*NA+b;
  numbonds[na]++;
#if 0
  if (numbonds[na]>4)
    {
      printf(">>>>>>>>>>> numbonds[%d]=%d\n", na, numbonds[na]);
      exit(-1);
    }
#endif
}

int bound(int na, int n, int a, int b)
{
  int i;
  for (i = 0; i < numbonds[na]; i++)
    if (bonds[na][i] == n*(NA*NA)+a*NA+b)
      return 1;
  return 0;
}
/* array con le posizioni degli atomi nel riferimento del corpo rigido 
 * nel caso dell'acqua i siti idrogeno ed elettroni sono disposti su 
 * di un tetraedro */
double rat_body[NA][3] = {{0,0,0},{0,0,1},{1,-1,0},{0,1,0},{0,0,1}};
extern void vectProdVec(double *A, double *B, double *C);
#if 1
#ifdef MD_SILICA
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3])
{
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  double r1[3], r2[3], r3[3], nr;
  double radius; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */
  /* WARNING: se le particelle che interagiscono hanno diametro diverso da sigma[0][1] qui va cambiato!!!! */ 
#ifdef MD_AB41
  if (i < Oparams.parnumA)
    radius = Oparams.sigma[0][0] / 2.0;
  else
    radius = Oparams.sigma[1][1] / 2.0;
  //printf("i=%d sigma=%f\n", i, radius*2.0);
#else
  radius = Oparams.sigma[0][1] / 2.0;
#endif
  if (ata == 0)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk];
      //printf("%f %f %f @ 0.5 C[red]\n", rat[0], rat[1], rat[2]);
    }
  else if (ata <= 3)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] + R[kk][ata-1]; 
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
      //printf("ata=%d %f %f %f @ 0.075 C[blue]\n", ata, R[0][ata-1], R[1][ata-1], R[2][ata-1]);
    }
  else
    {
      /* l'atomo restante è un electron site */
      for (kk = 0; kk < 3; kk++)
	{
	  r1[kk] = R[kk][1]-R[kk][0];
	  r2[kk] = R[kk][2]-R[kk][0];
	}
      vectProdVec(r1, r2, r3);
      nr = calc_norm(r3);
      for (kk = 0; kk < 3; kk++)
	r3[kk] *= radius/nr;
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] - r3[kk]; 
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
    }

}
#else
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3])
{
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  double r1[3], r2[3], r3[3], nr, fact;
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */

  fact = MD_DIST_HYDROSITES/MD_DIST_ELECTSITES;
  if (ata == 0)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk];
      //printf("%f %f %f @ 0.5 C[red]\n", rat[0], rat[1], rat[2]);
    }
  else if (ata <= 3)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] + R[kk][ata-1]; 
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
      //printf("ata=%d %f %f %f @ 0.075 C[blue]\n", ata, R[0][ata-1], R[1][ata-1], R[2][ata-1]);
    }
  else
    {
      /* l'atomo restante è un electron site */
      for (kk = 0; kk < 3; kk++)
	{
	  r1[kk] = R[kk][1]-R[kk][0];
	  r2[kk] = fact*R[kk][2]-R[kk][0];
	}
      vectProdVec(r1, r2, r3);
      nr = calc_norm(r3);
      for (kk = 0; kk < 3; kk++)
	r3[kk] *= MD_DIST_ELECTSITES/nr;
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] - r3[kk]; 
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
    }
}
#endif
void BuildAtomPos(int i, double *rO, double **R, double rat[5][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a, NUMAT;
  /* l'atomo zero si suppone nell'origine */
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
  if (i >= Oparams.parnumA)
    NUMAT = 4;
  else
    NUMAT = 3;
#elif defined(MD_AB41)
  if (i >= Oparams.parnumA)
    NUMAT = 2;
  else
    NUMAT = 5;
#else
  if (i >= Oparams.parnumA)
    NUMAT = 5;
  else
    NUMAT = 3;
#endif
#else
  NUMAT = 5;
#endif
  for (a=0; a < NUMAT; a++)
    BuildAtomPosAt(i, a, rO, R, rat[a]);

  //if (i==215)
    //printf("dist=%.15G\n", sqrt( Sqr(rat[1][0]-rat[4][0]) + Sqr(rat[1][1]-rat[4][1]) + Sqr(rat[1][2]-rat[4][2])));

}
#else
void body2lab(int i, double xp[], double x[], double *rO, double **R);
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[])
{
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  /* l'atomo zero si suppone nell'origine */
  if (ata == 0)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk];
    }
  else
    {
      body2lab(i, rat_body[ata], rat, rO, R);
    }
}
void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int kk, a;
  /* l'atomo zero si suppone nell'origine */
  for (kk = 0; kk < 3; kk++)
    rat[0][kk] = rO[kk];
  for (a=1; a < 5; a++)
    body2lab(i, rat_body[a], rat[a], rO, R);
}
#endif
void UpdateOrient(int i, double ti, double **Ro, double Omega[3][3], int bondpair)
{ 
  double wSq, w, OmegaSq[3][3], M[3][3];
  double sqrwx, sqrwy, sqrwz, sinw, cosw;
  int k1, k2, k3;
  sqrwx = Sqr(wx[i]);
  sqrwy = Sqr(wy[i]);
  sqrwz = Sqr(wz[i]);
  wSq = sqrwx+sqrwy+sqrwz;
  w = sqrt(wSq);
  if (w != 0.0 && ti != 0.0) 
    {
#if 0
      if (fabs(w*ti) < 1E-8)
	{
	  sinw = ti*(1-Sqr(w*ti)/6.0);	  
	  cosw = Sqr(ti)*(0.5 - Sqr(w*ti)/24.0);
	}
      else 
	{
	  sinw = sin(w*ti)/w;
	  cosw = (1.0 - cos(w*ti))/wSq;
	}
#endif
      sinw = sin(w*ti)/w;
      cosw = (1.0 - cos(w*ti))/wSq;
      Omega[0][0] = 0;
      Omega[0][1] = -wz[i];
      Omega[0][2] = wy[i];
      Omega[1][0] = wz[i];
      Omega[1][1] = 0;
      Omega[1][2] = -wx[i];
      Omega[2][0] = -wy[i];
      Omega[2][1] = wx[i];
      Omega[2][2] = 0;
      OmegaSq[0][0] = -sqrwy - sqrwz;
      OmegaSq[0][1] = wx[i]*wy[i];
      OmegaSq[0][2] = wx[i]*wz[i];
      OmegaSq[1][0] = OmegaSq[0][1];//wx[i]*wy[i];
      OmegaSq[1][1] = -sqrwx - sqrwz;
      OmegaSq[1][2] = wy[i]*wz[i];
      OmegaSq[2][0] = OmegaSq[0][2];//wx[i]*wz[i];
      OmegaSq[2][1] = OmegaSq[1][2];//wy[i]*wz[i];
      OmegaSq[2][2] = -sqrwx - sqrwy;

      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      //Omega[k1][k2] = -Omega[k1][k2];
	      M[k1][k2] = sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
#ifdef MD_USE_CBLAS
	      Rtmp2[k1][k2] = Rtmp[k1][k2] = R[i][k1][k2];
#endif
#if 0
	      if (k1==k2)
		M[k1][k1] += 1.0;
#endif
	    }
	}
#ifdef MD_USE_CBLAS
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  3, 3, 3, 1.0, &Rtmp[0][0],
		  3, &M[0][0], 3,
		  1.0, &Rtmp2[0][0], 3);
#if 1
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  Ro[k1][k2] = Rtmp2[k1][k2];
#endif
#else 
      for (k2 = 0; k2 < 3; k2++)
	{
	  /* se bondpair = 4 ci servono tutte e tre gli altri atomi per costruirlo
	   * quindi dobbiamo ruotarli tutti e 3! 
	   * NOTA: quindi bondpair = sticky point (bondpair={1,2,3,4}) */
	  if (bondpair != -1 && bondpair != 4 && (k2+1 != bondpair))
	    {
	      //printf("qui in UpdateOrient\n");
	      continue;
	    }
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      Ro[k1][k2] = R[i][k1][k2];
	      for (k3 = 0; k3 < 3; k3++)
		Ro[k1][k2] += M[k1][k3]*R[i][k3][k2];
	      //Ro[k1][k2] += R[i][k1][k3]*M[k3][k2];
	    }
	}
#endif
      //adjust_norm(Ro);
    }
  else
    {
      Omega[0][0] = 0;
      Omega[0][1] = 0;
      Omega[0][2] = 0;
      Omega[1][0] = 0;
      Omega[1][1] = 0;
      Omega[1][2] = 0;
      Omega[2][0] = 0;
      Omega[2][1] = 0;
      Omega[2][2] = 0;
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    Ro[k1][k2] = R[i][k1][k2];
	  }
    }
}

double calcDistNegOne(double t, double t1, int i, int j, int nn, double shift[3]);
extern double **matrix(int n, int m);
extern void free_matrix(double **M, int n);
int ibr, jbr, nnbr; 
double shiftbr[3], trefbr;
#if 1
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
void assign_bond_mapping(int i, int j)
{
  /* NOTA: l'interazione bonded è solo tra Si e O 
   * i <  Oparams.parnumA => O
   * i >=  Oparams.parnumA => Si */
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    {
      mapbondsa = mapbondsaAA;
      mapbondsb = mapbondsbAA;
    }
  else if (i >= Oparams.parnumA && j >= Oparams.parnumA)
    {
      mapbondsa = mapbondsaBB;
      mapbondsb = mapbondsbBB;

    }
  else if (i < Oparams.parnumA && j >= Oparams.parnumA)
    {
      mapbondsa = mapbondsbAB;
      mapbondsb = mapbondsaAB;
    }
  else
    {
      mapbondsa = mapbondsaAB;
      mapbondsb = mapbondsbAB;
    }
}
#elif defined(MD_AB41)
void assign_bond_mapping(int i, int j)
{
  /* NOTA: l'interazione bonded è solo tra Si e O 
   * i <  Oparams.parnumA => O
   * i >=  Oparams.parnumA => Si */
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    {
      mapbondsa = mapbondsaAA;
      mapbondsb = mapbondsbAA;
    }
  else if (i < Oparams.parnumA && j >= Oparams.parnumA)
    {
      mapbondsa = mapbondsaAB;
      mapbondsb = mapbondsbAB;
    }
  else
    {
      mapbondsa = mapbondsbAB;
      mapbondsb = mapbondsaAB;
    }
}
#else
void assign_bond_mapping(int i, int j)
{
  /* NOTA: l'interazione bonded è solo tra Si e O 
   * i <  Oparams.parnumA => O
   * i >=  Oparams.parnumA => Si */
  if (i < Oparams.parnumA)
    {
      mapbondsa = mapbondsbSiO;
      mapbondsb = mapbondsaSiO;
    }
  else
    {
      mapbondsa = mapbondsaSiO;
      mapbondsb = mapbondsbSiO;
    }
}
#endif
#endif
#endif
double funcs2beZeroed(double x, double tref, int i, int j, int nn, double shift[3])
{
#if 0
  int na, ata, atb; 
  double  rA[3], rB[3], ti;
  double Omega[3][3];
  /* x = (r, alpha, t) */ 
#ifdef MD_SILICA
  assign_bond_mapping(i, j);
#endif
  ata = mapbondsa[nn];
  atb = mapbondsb[nn];
  ti = x - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
  UpdateOrient(i, ti, Rt, Omega);
  na = (i < Oparams.parnumA)?0:1;
  if (OprogStatus.targetPhi > 0)
    {
      /* per ora la crescita nei primitive models non è implementata */
    }
  //tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = x - atomTime[j];
  MD_DEBUG(printf("atomTime[%d]:%.15f\n", j, atomTime[j]));
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, Rt, Omega);
  na = (j < Oparams.parnumA)?0:1;
  if (OprogStatus.targetPhi > 0)
    {
      /* per ora la crescita nei primitive models non è implementata */
    }
  //tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
#endif
  return calcDistNegOne(x, trefbr, i, j, nn, shift);
}
double  funcs2beZeroedBrent(double x)
{
  return funcs2beZeroed(x, trefbr, ibr, jbr, nnbr, shiftbr); 
}
double tdist;
double rA[3], rB[3];

double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
extern int check_point(char* msg, double *p, double *rc, double **XX);
extern void distconjgrad(int i, int j, double shift[3], double *vecg, double lambda, int halfspring);
extern int maxitsRyck;
extern double sigmaSqSticky;
double calcDistNegOne(double t, double t1, int i, int j, int nn, double shift[3])
{
  double distSq, ti;
  double ratA[NA][3], ratB[NA][3];
  int kk;
  double Omega[3][3];
  double sigmaSticky;
  int na;
#ifdef MD_SILICA
  assign_bond_mapping(i, j);
#endif
#ifdef MD_AB41
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    sigmaSticky = Oparams.sigmaStickyAA;
  else
    sigmaSticky = Oparams.sigmaStickyAB;
#else
  sigmaSticky = Oparams.sigmaSticky;
#endif
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  MD_DEBUG20(printf("BRENT nn=%d\n", nn));
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
  UpdateOrient(i, ti, RtA, Omega, mapbondsa[nn]);
  /* calcola le posizioni nel laboratorio degli atomi della molecola */
  BuildAtomPos(i, rA, RtA, ratA);
  na = (i < Oparams.parnumA)?0:1;
#if 0
  if (OprogStatus.targetPhi > 0)
    {
      /* qui deve scalare i raggi degli atomi che compongono la molecola */
    }
#endif
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, RtB, Omega, mapbondsb[nn]);
  na = (j < Oparams.parnumA)?0:1;
  BuildAtomPos(j, rB, RtB, ratB);
#if 0
  if (OprogStatus.targetPhi > 0)
    {
      /* qui deve scalare i raggi degli atomi che compongono la molecola */
    }
#endif
  /* calcola sigmaSq[][]!!! */
  distSq = 0;
  for (kk=0; kk < 3; kk++)
    distSq += Sqr(ratA[mapbondsa[nn]][kk]-ratB[mapbondsb[nn]][kk]);
  MD_DEBUG20(printf("dist= %.15G\n", sqrt(distSq)-sigmaSticky));
  return sqrt(distSq) - sigmaSticky;
#if 0
  if (firstdist || fabs(dist) < fabs(distmin))
    {
      firstdist = 0;
      distmin = dist;
      *amin = mapbondsa[nn];
      *bmin = mapbondsb[nn];
    }
#endif
}

/* N.B. per la silica tale routine va cambiata! */
double calcDistNeg(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, 
		   double dists[MD_PBONDS], int bondpair)
{
  double distmin, distSq, ti;
  double ratA[NA][3], ratB[NA][3], dist;
  int firstdist = 1, nn, kk;
  double Omega[3][3];
  double sigmaSticky;
  int na, npbonds;
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
#ifdef MD_AB41
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    sigmaSticky = Oparams.sigmaStickyAA;
  else
    sigmaSticky = Oparams.sigmaStickyAB;
#else
  sigmaSticky = Oparams.sigmaSticky;
#endif
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
  UpdateOrient(i, ti, RtA, Omega, (bondpair==-1)?-1:mapbondsa[bondpair]);
  /* calcola le posizioni nel laboratorio degli atomi della molecola */
  BuildAtomPos(i, rA, RtA, ratA);
  na = (i < Oparams.parnumA)?0:1;
#if 0
  if (OprogStatus.targetPhi > 0)
    {
      /* qui deve scalare i raggi degli atomi che compongono la molecola */
    }
#endif
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
 
  UpdateOrient(j, ti, RtB, Omega, (bondpair==-1)?-1:mapbondsb[bondpair]);
  na = (j < Oparams.parnumA)?0:1;
  BuildAtomPos(j, rB, RtB, ratB);
#if 0
  if (OprogStatus.targetPhi > 0)
    {
      /* qui deve scalare i raggi degli atomi che compongono la molecola */
    }
#endif
  /* calcola sigmaSq[][]!!! */
  distmin = 0;
  npbonds = set_pbonds(i, j);
  for (nn = 0; nn < npbonds; nn++)
    {
     if (bondpair != -1 && bondpair != nn)
	{
	  //printf("qui in calcDistNeg\n");
	  continue;
	}
      distSq = 0;
      for (kk=0; kk < 3; kk++)
	distSq += Sqr(ratA[mapbondsa[nn]][kk]-ratB[mapbondsb[nn]][kk]);
      dists[nn] = dist = sqrt(distSq) - sigmaSticky;
      if (firstdist || fabs(dist) < fabs(distmin))
	{
	  firstdist = 0;
	  distmin = dist;
	  *amin = mapbondsa[nn];
	  *bmin = mapbondsb[nn];
	}
    }
  return distmin;
}

void rebuildCalendar(void);
int vc_is_pos(int i, int j, double rCx, double rCy, double rCz,
	      double t)
{ 
  double rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double norm[3], wrx, wry, wrz, OmegaA[3][3], OmegaB[3][3];
  double modn;
  int na, a, b;
  MD_DEBUG(printf("[bump] t=%f contact point: %f,%f,%f \n", Oparams.time, rxC, ryC, rzC));
  rAC[0] = rx[i] + vx[i]*(t - atomTime[i]) - rCx;
  rAC[1] = ry[i] + vy[i]*(t - atomTime[i]) - rCy;
  rAC[2] = rz[i] + vz[i]*(t - atomTime[i]) - rCz;
#if 1
  for (a=0; a < 3; a++)
    if (fabs (rAC[a]) > L2)
      rAC[a] -= SignR(L, rAC[a]);
#endif
  rBC[0] = rx[j] + vx[j]*(t-atomTime[j]) - rCx;
  rBC[1] = ry[j] + vy[j]*(t-atomTime[j]) - rCy;
  rBC[2] = rz[j] + vz[j]*(t-atomTime[j]) - rCz;
#if 1
  for (a=0; a < 3; a++)
    {
      MD_DEBUG(printf("P rBC[%d]:%.15f ", a, rBC[a]));
      if (fabs (rBC[a]) > L2)
	rBC[a] -= SignR(L, rBC[a]);
      MD_DEBUG(printf("D rBC[%d]:%.15f ", a, rBC[a]));
    }
  MD_DEBUG(printf("\n"));
#endif 
  /* calcola tensore d'inerzia e le matrici delle due quadriche */
  na = (i < Oparams.parnumA)?0:1;

  UpdateOrient(i, t-atomTime[i], RA, OmegaA, -1);
#if 0
  /* durante la crescita gli spot vengono disattivati
     quindi questo if non serve */
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], RA);

  na = (j < Oparams.parnumA)?0:1;
  UpdateOrient(j, t-atomTime[j], RB, OmegaB, -1);
#if 0
  /* durante la crescita gli spot vengono disattivati
     quindi questo if non serve */
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
#endif
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], RB);
  for (a=0; a < 3; a++)
    {
      norm[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  norm[a] += -Xa[a][b]*rAC[b];
	}
    }
  modn = 0.0;
  for (a = 0; a < 3; a++)
    modn += Sqr(norm[a]);
  modn = sqrt(modn);
  for (a=0; a < 3; a++)
    norm[a] /= modn;
  /* calcola le velocità nel punto di contatto */
  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[0] = vx[i] + wrx;
  vCA[1] = vy[i] + wry;
  vCA[2] = vz[i] + wrz;
  vectProd(wx[j], wy[j], wz[j], -rBC[0], -rBC[1], -rBC[2], &wrx, &wry, &wrz);
  vCB[0] = vx[j] + wrx;
  vCB[1] = vy[j] + wry;
  vCB[2] = vz[j] + wrz;
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  MD_DEBUG(printf("VCPOS vc=%.15f\n", vc));
  return (vc > 0);
}

void evolveVec(int i, double ti, double *vecout, double *vecin)
{
  double wSq, w, OmegaSq[3][3], M[3][3], Omega[3][3];
  double sinw, cosw;
  int k1, k2;
  wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
  w = sqrt(wSq);
  if (w != 0.0) 
    {
      sinw = sin(w*ti)/w;
      cosw = (1.0 - cos(w*ti))/wSq;
      Omega[0][0] = 0;
      Omega[0][1] = -wz[i];
      Omega[0][2] = wy[i];
      Omega[1][0] = wz[i];
      Omega[1][1] = 0;
      Omega[1][2] = -wx[i];
      Omega[2][0] = -wy[i];
      Omega[2][1] = wx[i];
      Omega[2][2] = 0;
      OmegaSq[0][0] = -Sqr(wy[i]) - Sqr(wz[i]);
      OmegaSq[0][1] = wx[i]*wy[i];
      OmegaSq[0][2] = wx[i]*wz[i];
      OmegaSq[1][0] = wx[i]*wy[i];
      OmegaSq[1][1] = -Sqr(wx[i]) - Sqr(wz[i]);
      OmegaSq[1][2] = wy[i]*wz[i];
      OmegaSq[2][0] = wx[i]*wz[i];
      OmegaSq[2][1] = wy[i]*wz[i];
      OmegaSq[2][2] = -Sqr(wx[i]) - Sqr(wy[i]);

      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      //Omega[k1][k2] = -Omega[k1][k2];
	      M[k1][k2] = sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	    }
	}

      for (k1 = 0; k1 < 3; k1++)
	{
	  vecout[k1] = vecin[k1];
	  for (k2 = 0; k2 < 3; k2++)
	    vecout[k1] += M[k1][k2]*vecin[k1];
	}
    }
  else
    {
      for (k1 = 0; k1 < 3; k1++)
	{
	  vecout[k1] = vecin[k1];
	}
    }

}
#undef MD_DDOT_OPT
double calcvecF(int i, int j, double t, double *r1, double *r2, double* ddot, double shift[3])
{
  int kk;
  double rcat[3], rdbt[3], wra[3], wrb[3];
#ifdef MD_DDOT_OPT
  double normrcd, r12[3], sp;
#endif
  //evolveVec(i, t-t1, rcat, r1);
  //evolveVec(j, t-t1, rdbt, r2);
  rcat[0] = r1[0] - (rx[i] + vx[i]*(t-atomTime[i])); 
  rcat[1] = r1[1] - (ry[i] + vy[i]*(t-atomTime[i]));
  rcat[2] = r1[2] - (rz[i] + vz[i]*(t-atomTime[i]));
  rdbt[0] = r2[0] - (rx[j] + vx[j]*(t-atomTime[j]))-shift[0]; 
  rdbt[1] = r2[1] - (ry[j] + vy[j]*(t-atomTime[j]))-shift[1];
  rdbt[2] = r2[2] - (rz[j] + vz[j]*(t-atomTime[j]))-shift[2];
  ddot[0] = vx[i] - vx[j];
  ddot[1] = vy[i] - vy[j];
  ddot[2] = vz[i] - vz[j];
  vectProd(wx[i], wy[i], wz[i], rcat[0], rcat[1], rcat[2], &wra[0], &wra[1], &wra[2]);
  vectProd(wx[j], wy[j], wz[j], rdbt[0], rdbt[1], rdbt[2], &wrb[0], &wrb[1], &wrb[2]);
  for (kk=0; kk < 3; kk++)
    ddot[kk] += wra[kk] - wrb[kk];

#ifdef MD_DDOT_OPT
  for (kk=0; kk < 3; kk++)
    r12[kk] = r1[kk] - r2[kk];
  normrcd = 0;
  for (kk=0; kk < 3; kk++)
    normrcd += Sqr(r12[kk]);
  normrcd = sqrt(normrcd);

  for (kk=0; kk < 3; kk++)
    r12[kk] /= normrcd;
  sp = 0;
  for (kk=0; kk < 3; kk++)
    sp += r12[kk] * ddot[kk];
  return fabs(sp);
#else
  return calc_norm(ddot);
#endif  
}


const COORD_TYPE bc1 = 14.0/45.0, bc2 = 64.0/45.0, bc3 = 24.0/45.0;
/* =========================== >>> BodeTerm <<< ============================*/
double BodeTerm(double dt, double* fi)
{
  return dt * (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	       bc1 * fi[4]);
}
int refine_contact(int i, int j, double tref, double t1, double t2, int nn, double shift[3], double *troot)
{
  int kk;//, retcheck;

  polinterr=0;
  //newt(vecg, 5, &retcheck, funcs2beZeroed, i, j, shift); 
  ibr = i;
  jbr = j;
  nnbr = nn;
  trefbr = tref;
  for (kk=0; kk < 3; kk++)
    shiftbr[kk] = shift[kk];
  *troot=zbrent(funcs2beZeroedBrent, t1, t2, 1E-16);
  *troot += tref;
  if (polinterr==1)
    {
      MD_DEBUG10(printf("newt did not find any contact point!\n"));
      return 0;
    }
  else
    {
      return 1; 
    }
}
/* 0 = atomo grosso 
 * 1,2 = hydrogen sites
 * 3,4 = electon sites */
/* in tale modello non c'è interazione fra due
 * hydrogen sites o due electron sites */
/* N.B. l'urto 0-0 è tra due sfere dure nei primitive model 
 * dell'acqua e della silica, quindi lo tratto a parte.
 * Inoltre non c'è interazione tra un atomo grosso e un atomo sticky. */
int check_cross(double distsOld[MD_PBONDS], double dists[MD_PBONDS], 
		int crossed[MD_PBONDS], int bondpair, int npbonds)
{
  int nn;
  int retcross = 0;
 
  for (nn = 0; nn < npbonds; nn++)
    {
      crossed[nn] = MD_EVENT_NONE;
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (dists[nn]*distsOld[nn] < 0)
	{
	  if (distsOld[nn] > 0)
	    crossed[nn] = MD_OUTIN_BARRIER; 
	  else
	    crossed[nn] = MD_INOUT_BARRIER;
	  retcross = 1;
	}
    }
  return retcross;
}
int check_cross_strictcheck(int i, int j, double distsOld[MD_PBONDS], double dists[MD_PBONDS], 
		int crossed[MD_PBONDS], int bondpair)
{
  int nn;
  int retcross = 0, npbonds;

  npbonds = set_pbonds(i, j);
  for (nn = 0; nn < npbonds; nn++)
    {
      crossed[nn] = MD_EVENT_NONE;
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (dists[nn] > 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	{
	  crossed[nn] = MD_INOUT_BARRIER; 
	  retcross = 1;
	}	  
      if (dists[nn] < 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	{
	  crossed[nn] = MD_OUTIN_BARRIER;
	  retcross = 1;
	}
    }
  return retcross;
}
int check_cross_strictcheck_sf(int i, int j, double distsOld[MD_PBONDS], double dists[MD_PBONDS], int crossed[MD_PBONDS], int bondpair)
{
  int nn;
  int retcross = 0, npbonds;
  npbonds = set_pbonds(i, j);
  for (nn = 0; nn < npbonds; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (fabs(dists[nn]) < OprogStatus.epsdFast)
	return 1;
      if (dists[nn] > 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
      if (dists[nn] < 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
    }
  return retcross;
}


int get_dists_tocheck(double distsOld[], double dists[], int tocheck[], int dorefine[],
		      int bondpair, int npbonds)
{
  int nn;
  int rettochk = 0;
  for (nn = 0; nn < npbonds; nn++)
    {
      tocheck[nn] = 0;
      if ( dists[nn]*distsOld[nn] > 0.0 &&
	  fabs(dists[nn]) < OprogStatus.epsd && fabs(distsOld[nn]) < OprogStatus.epsd &&
	  dorefine[nn] == MD_EVENT_NONE && (bondpair== -1 || bondpair == nn))
	{
	  tocheck[nn] = 1; 
	  rettochk++;
	}
    }
  return rettochk;
}
double get_max_deldist(double distsOld[MD_PBONDS], double dists[MD_PBONDS], int bondpair, int npbonds)
{
  int nn, first = 1;
  double maxdd=0.0, dd;
  for (nn = 0; nn < npbonds; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
      dd = fabs(dists[nn]-distsOld[nn]);
      if (first || dd > maxdd)
	{
	  first = 0;
	  maxdd = dd;
	}
    }
  return maxdd;
}

void assign_dists(double a[], double b[])
{
  memcpy(b, a, MD_PBONDS*sizeof(double));
}
#define MD_OPTDDIST
/* NOTA: tale stima ottimizzata della maggiorazione per la velocità di variazione della distanza
 * sembra corretta, fare comunque dei test.*/
#if 1 
double eval_maxddist(int i, int j, int bondpair, double t1, double *maxddotOpt)
{
  double ti, rA[3], rB[3], Omega[3][3], ratA[NA][3], ratB[NA][3], wri[3], wrj[3], nwri, nwrj,
	 r12i[3], r12j[3];//, maxddotOpt[MD_PBONDS];
  double maxddot=0.0, nr12i, nr12j;
  double sig2, factori, factorj;
  int nn, kk, npbonds;
  ti = t1 - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
  UpdateOrient(i, ti, RtA, Omega, (bondpair==-1)?-1:mapbondsa[bondpair]);
  BuildAtomPos(i, rA, RtA, ratA);
  ti = t1 - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti;
  rB[1] = ry[j] + vy[j]*ti;
  rB[2] = rz[j] + vz[j]*ti;
  /* ...and now orientations */
  UpdateOrient(j, ti, RtB, Omega, (bondpair==-1)?-1:mapbondsb[bondpair]);
  BuildAtomPos(j, rB, RtB, ratB);
  npbonds = set_pbonds(i, j);
  for (nn = 0; nn < npbonds; nn++)
    {
      for (kk = 0; kk < 3; kk++)
	{
	  r12i[kk] = (ratA[mapbondsa[nn]][kk]-rA[kk]);
  	  r12j[kk] = (ratB[mapbondsb[nn]][kk]-rB[kk]);	  
	}
      nr12i = calc_norm(r12i);
      nr12j = calc_norm(r12j);
#ifdef MD_AB41
      if (i < Oparams.parnumA && j < Oparams.parnumA)
	sig2 = 0.5*Oparams.sigmaStickyAA;
      else
	sig2 = 0.5*Oparams.sigmaStickyAB;
#else
      sig2 = 0.5*Oparams.sigmaSticky;
#endif
      nr12i += sig2;
      nr12j += sig2;
      factori = nr12i + OprogStatus.epsd;
      factorj = nr12j + OprogStatus.epsd;
      maxddotOpt[nn] = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
	sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori + 
	sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;
      if (OprogStatus.assumeOneBond && nn==bondpair)
	{
	  maxddot = maxddotOpt[nn];
	  return maxddot;
	}
      else
	{
	  if (nn==0 || maxddotOpt[nn] > maxddot)
	    maxddot = maxddotOpt[nn];
	}
    }
  return maxddot;
}
#else
double eval_maxddist(int i, int j, int bondpair, double t1, double *maxddotOpt)
{
  double ti, rA[3], rB[3], Omega[3][3], ratA[NA][3], ratB[NA][3], wri[3], wrj[3], nwri, nwrj,
	 r12i[3], r12j[3];//, maxddotOpt[MD_PBONDS];
  double maxddot=0.0, nr12i, nr12j;
#if 1
  double sig2;
#endif
  int nn, kk, npbonds;
  ti = t1 - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
  UpdateOrient(i, ti, RtA, Omega, (bondpair==-1)?-1:mapbondsa[bondpair]);
  BuildAtomPos(i, rA, RtA, ratA);
  ti = t1 - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti;
  rB[1] = ry[j] + vy[j]*ti;
  rB[2] = rz[j] + vz[j]*ti;
  /* ...and now orientations */
  UpdateOrient(j, ti, RtB, Omega, (bondpair==-1)?-1:mapbondsb[bondpair]);
  BuildAtomPos(j, rB, RtB, ratB);
  npbonds = set_pbonds(i, j);
  for (nn = 0; nn < npbonds; nn++)
    {
      for (kk = 0; kk < 3; kk++)
	{
	  r12i[kk] = (ratA[mapbondsa[nn]][kk]-rA[kk]);
  	  r12j[kk] = (ratB[mapbondsb[nn]][kk]-rB[kk]);	  
	}
      nr12i = calc_norm(r12i);
      nr12j = calc_norm(r12j);
#if 0 
#ifdef MD_AB41
      if (i < Oparams.parnumA && j < Oparams.parnumA)
	sig2 = 0.5*Oparams.sigmaStickyAA;
      else
	sig2 = 0.5*Oparams.sigmaStickyAB;
#else
      sig2 = 0.5*Oparams.sigmaSticky;
#endif
      nr12i += sig2;
      nr12j += sig2;
#endif
      for (kk = 0; kk < 3; kk++)
	{
	  //printf("nr12i=%.15G nr12j=%.15G\n", nr12i, nr12j);
	  r12i[kk] *= (nr12i+OprogStatus.epsd)/nr12i;
	  r12j[kk] *= (nr12j+OprogStatus.epsd)/nr12j;
	}	  
      vectProd(wx[i], wy[i], wz[i], r12i[0], r12i[1], r12i[2], &wri[0], &wri[1], &wri[2]);
      nwri = calc_norm(wri);
      vectProd(wx[j], wy[j], wz[j], r12j[0], r12j[1], r12j[2], &wrj[0], &wrj[1], &wrj[2]);
      nwrj = calc_norm(wrj);
      maxddotOpt[nn] = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
	nwri + nwrj;
      if (OprogStatus.assumeOneBond && nn==bondpair)
	{
	  maxddot = maxddotOpt[nn];
	  return maxddot;
	}
      else
	{
	  if (nn==0 || maxddotOpt[nn] > maxddot)
	    maxddot = maxddotOpt[nn];
	}
    }
  return maxddot;
}
#endif
void calc_delt(double *maxddoti, double *delt, double *dists, int bondpair, int npbonds)
{
  int nn;
  double dt;
  for (nn = 0; nn < npbonds; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
      dt = fabs(dists[nn]) / maxddoti[nn];
      //printf("nn=%d dt=%.15G delt=%.15G dists=%.15G maxddoti=%15G\n", nn, dt, *delt, dists[nn], maxddoti[nn]);
      if (nn==0 || bondpair != -1 || dt < (*delt))
	*delt = dt;
    }
  //printf("I chose dt=%.15G\n", *delt);
}

int search_contact_faster(int i, int j, double *shift, double *t, double t1, double t2, double epsd, double *d1, double epsdFast, double dists[MD_PBONDS], int bondpair, double maxddot, double *maxddoti)
{
  /* NOTA: 
   * MAXOPTITS è il numero massimo di iterazioni al di sopra del quale esce */
  double told, delt, distsOld[MD_PBONDS];
  const int MAXOPTITS = 500;
  int its=0, amin, bmin, crossed[MD_PBONDS], npbonds; 
  /* estimate of maximum rate of change for d */
#if 0
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]*0.5
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j]*0.5;
#endif
  *d1 = calcDistNeg(*t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
  MD_DEBUG30(printf("[IN SEARCH CONTACT FASTER]*d1=%.15G t=%.15G\n", *d1, *t));
  timesF++;
  MD_DEBUG10(printf("Pri distances between %d-%d d1=%.12G epsd*epsdTimes:%f\n", i, j, *d1, epsdFast));
  told = *t;
  delt = OprogStatus.h;
  npbonds = set_pbonds(i, j);
  if (fabs(*d1) < epsdFast)
    {
      assign_dists(distsOld, dists);
      return 0;
    }
  while (fabs(*d1) > epsdFast && its < MAXOPTITS)
    {
      if (maxddot*(t2-(t1+*t)) < fabs(*d1)-OprogStatus.epsd)
	return 1;
#ifdef MD_OPTDDIST
      calc_delt(maxddoti, &delt, distsOld, bondpair, npbonds);
#else
      delt = fabs(*d1) / maxddot;
#endif
#if 0
      /* CALCOLARE normddot usando la distanza vecchia come in locate_contact */
      normddot = calcvecF(i, j, *t, r1, r2, ddot, shift);
      //printf("normddot: %.15G\n", epsd/normddot);
      /* check for convergence */

      if (normddot!=0 && delt < (epsd / normddot))
	{
	  MD_DEBUG10(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
#endif
      //printf("[SEARCH_CONTACT_FASTER] t=%.15G maxddot=%.15G\n", *t, eval_maxddist(i, j, bondpair, *t));
      *t += delt;
      *d1 = calcDistNeg(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
      if (check_cross(distsOld, dists, crossed, bondpair, npbonds))
	{
	  /* go back! */
	  MD_DEBUG30(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNeg(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 0;
	}
#if 0
      if (check_cross_strictcheck_sf(i,j,distsOld, dists, crossed, bondpair))
	{
	  /* go back! */
	  MD_DEBUG30(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNeg(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 0;
	}
#endif
      if (*t+t1 > t2)
	{
	  *t = told;
	  MD_DEBUG30(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("convergence t>t2\n"));
	  *d1 = calcDistNeg(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 1;
	}
      told = *t;
      assign_dists(dists, distsOld);
      its++;
      itsF++;
    }

  MD_DEBUG10(printf("max iterations %d iterations reached t=%f t2=%f\n", its, *t, t2));
  return 0;
}
extern double **Aip;
extern void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double xa[3], ya[3];
double distfunc(double x)
{
  double dy, y;
  polint(xa, ya, 3, x, &y, &dy);
  if (polinterr==1)
    return 0.0;
  if (dy > OprogStatus.epsd)
    {
      if (OprogStatus.phitol <= 0)
	printf("dy=%.15G\n", dy);
      polinterr = 1;
    }
  else 
    polinterr = 0;
  return y;
}
#if 0
double choose_seg(int i, int j, double tref, int nn, double shift[3], double xa[3], double ya[3],
		  double *tmin, double *dmin)
{
  double xal[3], yal[3], xar[3], yar[3], tminr, tminl, dminr, dminl;
  double A, B, C;
  xal[0] = xa[0];
  yal[0] = ya[0];
  xal[1] = (xal[0]+xal[1])*0.5;
  yal[1] = calcDistNegOne(xal[1], tref, i, j, nn, shift);
  xal[2] = xa[1];
  yal[2] = ya[1];
  A = xal[2]*(yal[0]-yal[1]);
  B = xal[0]*(yal[1]-yal[2]);
  C = xal[1]*(yal[2]-yal[0]);
  tminl = (xal[2]*A+xal[0]*B+xal[1]*C)/(A+B+C)/2.0;
  dminl = calcDistNegOne(tminl, tref, i, j, nn, shift);
  xar[0] = xa[1];
  yar[0] = ya[1];
  xar[1] = (xal[1]+xal[2])*0.5;
  yar[1] = calcDistNegOne(xar[1], tref, i, j, nn, shift);
  xar[2] = xa[2];
  yar[2] = ya[2];
  A = xar[2]*(yar[0]-yar[1]);
  B = xar[0]*(yar[1]-yar[2]);
  C = xar[1]*(yar[2]-yar[0]);
  tminr = (xar[2]*A+xar[0]*B+xar[1]*C)/(A+B+C)/2.0;
  dminr = calcDistNegOne(tminr, tref, i, j, nn, shift);
  if ((tminl < xal[0] || tminl > xal[2]) && (tminr < xar[0] || tminr > xar[2]))
   return 1;
  else if (tminl > xal[0] && tminl < xal[2])
    {
      xa[0] = xal[0];
      ya[0] = yal[0];
      xa[1] = xal[1];
      ya[1] = yal[1];
      xa[2] = xal[2];
      ya[2] = yal[2];
      *tmin = tminl;
      *dmin = dminl;
      return 0;
    }
  else
    {
      xa[0] = xar[0];
      ya[0] = yar[0];
      xa[1] = xar[1];
      ya[1] = yar[1];
      xa[2] = xar[2];
      ya[2] = yar[2];
      *tmin = tminr;
      *dmin = dminr;
      return 0;
    }
  if (fabs(dminl) < fabs(dminr))
    {
      xa[0] = xal[0];
      ya[0] = yal[0];
      xa[1] = xal[1];
      ya[1] = yal[1];
      xa[2] = xal[2];
      ya[2] = yal[2];
      *tmin = tminl;
      *dmin = dminl;
      return 0;
    }
  else
    {
      xa[0] = xar[0];
      ya[0] = yar[0];
      xa[1] = xar[1];
      ya[1] = yar[1];
      xa[2] = xar[2];
      ya[2] = yar[2];
      *tmin = tminr;
      *dmin = dminr;
      return 0;
    }
} 
#endif
int interpol(int i, int j, int nn, 
	     double tref, double t, double delt, double d1, double d2,
	     double *tmin, double shift[3], int ignoresignchg)
{
  double d3, A, dmin;
  /* NOTA: dists di seguito può non essere usata? controllare!*/
  d3 = calcDistNegOne(t+delt*0.5, tref, i, j, nn, shift);
  xa[0] = 0;
  ya[0] = d1;
  xa[1] = delt*0.5;
  ya[1] = d3;
  xa[2] = delt;
  ya[2] = d2;
  if (ya[0]-ya[1] == 0.0)
    {
      *tmin = t + delt*0.25;
    }
  else if (ya[2]-ya[0] ==0.0)
    {
      *tmin = t + delt*0.5;
    }
  else
    {      
      A = (ya[2]-ya[0])/(ya[0]-ya[1]);
      *tmin = t + 0.5*delt*((1.0 + A * 0.25)/( 1.0 + A * 0.5));
    }
  dmin = calcDistNegOne(*tmin, tref, i, j, nn, shift);
#if 0
  printf("delt=%.15G *tmin=%.15G BAH=%.15G\n", delt, *tmin, 0.5*delt*((1.0 + A * 0.25)/( 1.0 + A * 0.5)));
  printf("A=%.15G i=%d j=%d\n", A, i, j);
  printf("{{%.15G,%.15G},{%.15G,%.15G},{%.15G,%.15G}} - *tmin=%.15G,%.15G\n", t, d1, t+delt*0.5,d3, t+delt, d2, *tmin, dmin);
#endif
  if (*tmin < t+delt && *tmin > t)
    {
      *tmin += tref;
      if (!ignoresignchg)
	{
	  if (d1*dmin < 0.0)
	    return 0;
	  /* we can call grazing_try_harder() (i.e. brent) only if dmin is less than d1 and d2 (bracketing condition) */
	  else if (fabs(dmin) < fabs(d1) && fabs(dmin) < fabs(d2))
	    return 2;
	  else
	    return 1;
	}
      else
	return 0;
    }
#if 0
  if (OprogStatus.tryharder==0)
    return 1;
  if (*tmin > t + delt || *tmin < t) 
    return 1;
  /* inizia la dicotomia */
  for (its = 0; its < OprogStatus.tryharder; its++)
    {
      xa[1] = *tmin;
      ya[1] = dmin;
      if (choose_seg(i, j, tref, nn, shift, xa, ya, &tminnew, &dminnew))
	return 1;
      if (fabs(dminnew) > fabs(dmin))
	{
 	  if (dmin*d1 < 0.0)
	    {
	      printf("beccato provando duramente :) its=%d\n", its);
	      *tmin += tref;
	      return 0;
	    }
	  return 1;
	}
      dmin = dminnew;
      *tmin = tminnew;
    }
#endif
   // printf("t=%.15G tmin=%.15G t+delt=%.15G dmin=%.15G d1=%.15G d2=%.15G d3=%.15G\n", t, *tmin, t+delt, dmin, d1, d3, d2);
  return 1;
}

int valid_collision(int i, int j, int ata, int atb, int collCode)
{
  MD_DEBUG30(printf("lastbump[i=%d].mol=%d lastbump[j=%d].mol=%d lastbump[i=%d].at=%d lastbump[j=%d].at=%d\n",
	i, lastbump[i].mol, j, lastbump[j].mol, i, lastbump[i].at, j, lastbump[j].at));
  MD_DEBUG30(printf("collCode=%d ata=%d atb=%d bound: %d\n", collCode, ata, atb,  bound(i, j, ata, atb)));
  
  if ((collCode==MD_INOUT_BARRIER && !bound(i, j, ata, atb)) ||
      (collCode==MD_OUTIN_BARRIER && bound(i, j, ata, atb)) ) 
    return 0; 
  else
    return 1;
}
int get_bonded(int i, int j)
{
  int nb, jj, jj2, kk, nn, aa, bb, npbonds;
  nb = numbonds[i];
  if (!OprogStatus.assumeOneBond)
    return -1;
  npbonds = set_pbonds(i, j);
  for (kk = 0; kk < nb; kk++)
    {
      jj = bonds[i][kk] / (NA*NA);
      jj2 = bonds[i][kk] % (NA*NA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (jj == j)
	{
	  for (nn=0; nn < npbonds; nn++)
	    if (mapbondsa[nn]==aa && mapbondsb[nn]==bb)
	      return nn;
	} 
    }
  return -1; 

}

int check_negpairs(int *negpairs, int bondpair, int i, int j)
{
  int nn, sum, npbonds;
  sum = 0;
//  if (lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == 0 
  //    && lastbump[j].at == 0)
    //return 2;
  npbonds = set_pbonds(i, j);
  for (nn = 0; nn < npbonds; nn++)
    {
      negpairs[nn] = 0;
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (!(lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == mapbondsa[nn]
	&& lastbump[j].at == mapbondsb[nn]))
	continue;
      negpairs[nn] = 1;
      return nn+1;
#if 0
      if (bound(i, j, mapbondsa[nn], mapbondsb[nn]) && dists[nn] > 0.0)
	negpairs[nn] = 1;
      else if (!bound(i, j, mapbondsa[nn], mapbondsb[nn]) && dists[nn] < 0.0)
	negpairs[nn] = 1;
      sum += negpairs[nn];
#endif
      //printf("bondpair: %d dists[%d]:%.15G\n", bondpair, nn, dists[nn]);
    }
  return 0;
}

int delt_is_too_big_hc(int i, int j, int bondpair, double *dists)
{
  int nn, npbonds;
  npbonds = set_pbonds(i, j);
  for (nn=0; nn < npbonds; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
   if (dists[nn] > 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
      return 1;
    if (dists[nn] < 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
      return 1;
    }
  return 0;
}
int delt_is_too_big(int i, int j, int bondpair, double *dists, double *distsOld,
		    int *negpairs)
{
  int nn, npbonds;
  npbonds = set_pbonds(i, j);
  for (nn=0; nn < npbonds; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (!negpairs[nn])
	continue;
#if 0
      if (!(lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == mapbondsa[nn]
	&& lastbump[j].at == mapbondsb[nn]))
	continue;
#endif
#if 0
      if (distsOld[nn] > 0.0 && dists[nn] > 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
      if (distsOld[nn] < 0.0 && dists[nn] < 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
#endif
      /* N.B. distsOld[nn] non va controllato poiché dopo un urto cmq ci deve essere un estremo
       * di distanza dmin t.c. dmin > 0 se !bound(i,j..) o dmin < 0 se bound(i,j...) */
      if (dists[nn] >= 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1; 
      if (dists[nn] <= 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
    }
  return 0;
}
extern double max(double a, double b);
#define MD_BASIC_DT
double brentSign, brent_tref, shiftBrent[3];
int iBrent, jBrent, nnBrent;
double distSPbrent(double t)
{
  return brentSign*calcDistNegOne(t, brent_tref, iBrent, jBrent, nnBrent, shiftBrent);
}
#define MD_BRENT_TOL 1E-15
extern double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
int grazing_try_harder(int i, int j, int nn, double tref, double t1, double delt, double d1, double d2, double shift[3], double *troot, double *dmin)
{
  int a;
#ifndef MD_GRAZING_TRYHARDER
  return 0;
#endif
  printf("[grazing_try_harder] i=%d j=%d nn=%d time=%.15G\n", i, j, nn, tref+t1);
  /* Brent looks always for a minimum hence we have to change sign
     if grazing occurrs coming from negative distances */
  brentSign = ((d1>0.0)?1.0:-1.0);
  brent_tref = tref;
  iBrent = i;
  jBrent = j;
  nnBrent = nn;
  for (a=0; a < 3; a++)
    shiftBrent[a] = shift[a];
  /* use brent to find the exact minimum */
  *dmin = brent(t1, t1+delt*0.5, t1+delt, distSPbrent, MD_BRENT_TOL, troot);
  *dmin *= brentSign;
  if (*troot >= t1 && *troot <= t1+delt && *dmin*d1 < 0.0)
    {
      /* found a crossing! */
      *troot += tref;
      return 1;
    }
  return 0;/* no collision found */
}

int locate_contact(int i, int j, double shift[3], double t1, double t2, 
		   double *evtime, int *ata, int *atb, int *collCode)
{
  const double minh = 1E-20;
  double h, d, dold, t2arr[MD_PBONDS], t, dists[MD_PBONDS], distsOld[MD_PBONDS]; 
  double maxddot, delt, troot, tmin, tini; //distsOld2[MD_PBONDS];
  int retip;
  double dmin, deltini;
  //const int MAXOPTITS = 4;
  int bondpair, itstb, adjt1=0, npbonds;
  int its, foundrc;
#if 0
  //const int MAXITS = 100;
  //int itsRef;
  //int goback;
  double t1ini, delthc;
#endif
#ifndef MD_BASIC_DT
  double dold2, deltth, normddot, distsOld2[MD_PBONDS], deldist;
#endif
  double maxddoti[MD_PBONDS], epsd, epsdFast, epsdFastR, epsdMax; 
  int tocheck[MD_PBONDS], dorefine[MD_PBONDS], ntc, ncr, nn, gotcoll, amin, bmin,
      crossed[MD_PBONDS], firstaftsf;
#ifdef MD_NEGPAIRS
  int negpairs[MD_PBONDS], sumnegpairs;
#endif
  const double GOLD= 1.618034;
  epsd = OprogStatus.epsd;
  epsdFast = OprogStatus.epsdFast;
  epsdFastR= OprogStatus.epsdFastR;
  epsdMax = OprogStatus.epsdMax;
#ifdef MD_SPOT_OFF
  return 0;
#endif
#ifdef MD_GROWTH_CODE
  /* switch off spots during growth */
  if (OprogStatus.targetPhi > 0.0)
    return 0;
#endif
#ifdef MD_SILICA
  assign_bond_mapping(i, j);
#endif
  npbonds = set_pbonds(i, j);
  /* NOTA: 
   * - epsd è di quanto varia d ad ogni iterazione e quindi determina il grado di accuratezza
   * con cui viene individuato il punto di contatto. In generale se due ellissoidi si "spizzicano"
   * ad una distanza minore di epsd tali urti non vengono rilevati 
   * - epsdTimes*epsd è la soglia sotto la quale la ricerca veloce fatta in search_contact_faster termina 
   * - epsdTimesIsteresi*epsd è la sogli al di sopra della quale viene di nuovo usata search_contact_faster
   *   Tale valore è bene che sia abbastanza grande di modo che non faccia continue ricerche veloci e lente 
   *   che rallentano moltissimo. Infatti tale ricerca veloce serve solo nel caso in cui due ellissoidi si 
   *   sfiorano per poi allontanrsi. 
   */
  t = 0.0;
  bondpair = get_bonded(i, j);
#ifdef MD_OPTDDIST
  maxddot = eval_maxddist(i, j, bondpair, t1, maxddoti);
#else
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]*0.5
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j]*0.5;
  //printf("wx:%f maxax: %f,%f\n", sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i])), maxax[i], maxax[j]);
#endif

  //printf("[LOCATE_CONTACT INIZIO] t=%.15G maxddot=%.15G\n", t1, eval_maxddist(i, j, bondpair, t1));
  MD_DEBUG30(printf("[locate_contact] %d-%d t1=%f t2=%f shift=(%f,%f,%f)\n", i,j,t1, t2, shift[0], shift[1], shift[2]));
  h = OprogStatus.h; /* last resort time increment */
  if (*collCode!=MD_EVENT_NONE)
    {
      if (t2 > *evtime)
	t2 = *evtime+1E-6;
    }
  delt = h;
  MD_DEBUG(printf("QUIIII collCode=%d\n", *collCode));
#ifndef MD_NEGPAIRS
  /* NOTA: le strategie per evitare problemi dopo una collisione sono due:
   * 1) andare avanti nel tempo finché la distanza non è corretta.
   * 2) fare un passo ed eventualmente ridurlo finchè la distanza non è corretta.
   */
  df = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
  for (nn=0; nn < npbonds; nn++)
    {
      if (!do_check_negpairs)
	break;
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (!(lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == mapbondsa[nn]
		       && lastbump[j].at == mapbondsb[nn]))
	continue;
      if (dists[nn] > 0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	{
	  its = 0;
	  df = dists[nn];
	  while (df > 0 && its < MAXITS)
	    {
	      t += h;
	      its++;
	      if (t + t1 > t2)
		return 0;
	      df = calcDistNegOne(t, t1, i, j, nn, shift);
	    }
	}
      if (dists[nn] < 0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	{
	  its = 0;
	  df = dists[nn];
	  while (df < 0 && its < MAXITS)
	    {
	      t += h;
	      its++;
	      if (t + t1 > t2)
		return 0;
	      df = calcDistNegOne(t, t1, i, j, nn, shift);
	    }
	}
    }
#endif
  MD_DEBUG30(printf("[BEFORE SEARCH CONTACT FASTER]Dopo distances between %d-%d t=%.15G t2=%.15G\n", i, j, t, t2));
#ifdef MD_NEGPAIRS
  /* tale check ha senso solo appena dopo un urto fra i e j */
  if (do_check_negpairs)
    sumnegpairs = check_negpairs(negpairs, bondpair, i, j); 
  else
    sumnegpairs = 0;
#endif
#if 0
#ifdef MD_NEGPAIRS
  /* NOTA: inizia poco prima di t1 se l'ultimo urto tra le molecole è stata una collisione delle sfere dure 
   * per evitare problemi legati al fatto che il punto iniziale in tale caso puo' essere molto a ridosso 
   * del crossing (ciò accade nella regione intorno all'intersezione della sticky spheres con la sfera dura. */
  d = calcDistNeg(0, t1, i, j, shift, &amin, &bmin, dists, bondpair);
  adjt1 = 0;  
  if (lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == 0 
      && lastbump[j].at == 0 && fabs(d) < epsd)
    {
      delthc = max(epsdFast/maxddot,OprogStatus.h);
      //printf("INIZIO t1=%.18G delthc=%.18G d=%.20G\n", t1, delthc, d);
      adjt1 = 1;
      t1ini = t1;
      t1 -= delthc;
      d = calcDistNeg(0, t1, i, j, shift, &amin, &bmin, dists, bondpair);
      while (delt_is_too_big_hc(i, j, bondpair, dists) && 
	     delthc > minh)
    	{
	  delthc /= GOLD; 
	  t1 = t1ini - delthc;
	  d = calcDistNeg(0, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  //printf("d=%.15G t1=%.15G delthc=%.15G tini=%.15G\n", d, t1, delthc, t1ini);
	}
	
    }
#endif
#endif 
#if 0
  dold=calcDistNeg(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
  d=calcDistNeg(t+1E-20, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
#endif
  MD_DEBUG30(
  if (sumnegpairs)
    printf("_======> d=%.30G\n", calcDistNeg(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair));
    )
  if (search_contact_faster(i, j, shift, &t, t1, t2, epsd, &d, epsdFast, dists, bondpair, maxddot, maxddoti))
    {
      return 0;  
    }
  timesS++;
  MD_DEBUG30(printf("[AFTER SEARCH CONTACT FASTER]Dopo distances between %d-%d d1=%.12G\n", i, j, d));
   
#if 0
  /* N.B. prova per vedere il minimo delt */
  dold = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
  d = calcDistNeg(t+1E-18, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
  printf("d:%.15G dold: %.15G d-dold=%.15G\n", d, dold, d-dold);
#endif
  MD_DEBUG(printf(">>>>d:%f\n", d));
  foundrc = 0;
#if 1
  assign_dists(dists, distsOld);
  dold = d;
#else
  dold = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
#endif
  firstaftsf = 1;
  its = 0;
  while (t+t1 < t2)
    {
#ifdef MD_BASIC_DT
      tini = t;
      delt = epsd/maxddot;
      t += delt;
      while ((d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair))==0.0)
	{
	  delt *= GOLD;
	  t = tini + delt;
	}
#else
#if 0
      deldist = get_max_deldist(distsOld2, distsOld);
      normddot = fabs(deldist)/delt;
      /* NOTA: forse qui si potrebbe anche usare sempre delt = epsd/maxddot */
      if (normddot!=0)
	{
	  delt = epsd/normddot;
	}
      else
	delt = h;
      if (fabs(dold) < epsd)
	{
	  delt = epsd / maxddot;
	}
#else
      if (!firstaftsf)
	{
	  deldist = get_max_deldist(distsOld2, distsOld, bondpair, npbonds);
	  normddot = fabs(deldist)/delt;
	  /* NOTA: forse qui si potrebbe anche usare sempre delt = epsd/maxddot */
	  if (normddot!=0)
	    {
	      delt = epsd/normddot;
	    }
	  else
	    {
	      //if (fabs(t)<1E-15)
	      delt = epsd/maxddot;
	      //delt = h;
	      //else
	      //delt = t*h;
	    }

	  if (fabs(dold) < epsd)
	    {
	      delt = epsd / maxddot;
	    }
	}
      else
	{
	  delt = h;
	  firstaftsf = 0;
	  //t += delt;
	  dold2 = calcDistNeg(t-delt, t1, i, j, shift, &amin, &bmin, distsOld2, bondpair);
	  MD_DEBUG30(printf("==========>>>>> t=%.15G t2=%.15G\n", t, t2));
	  //assign_dists(distsOld,  distsOld2);
	  //assign_dists(dists, distsOld);
	  continue;
	}
#endif
#if 0
      /* se all'inizio c'erano sticky spots che si overlappavano finché le distanze non sono corrette
       * usa il passo minimo (dell'ordine della precisione di macchina) 
       * NOTA: evita di fare più di MAXITS passi "minimi" */
      if (its < MAXITS && use_min_delt(negpairs, distsOld, bondpair))
	{
	  its++;
	  delt = sh;
	}
#endif
      MD_DEBUG30(printf("delt: %f epsd/maxddot:%f h*t:%f maxddot:%f\n", delt, epsd/maxddot,h*t,maxddot));
      ///delt = h;///
#if 0
      if (sumnegpairs)
	delt = 1E-6;
#endif
      tini = t;
      t += delt;
      //printf("normddot= %.15G t=%.15G delt=%.15G maxddot: %.15G t*h=%.15G\n", 
      //   normddot, t, delt, maxddot, t*h);
      //printf("normddot=%f dt=%.15G\n",normddot, epsd/normddot); 
      //dold2 = dold;
      d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);

      deldist = get_max_deldist(distsOld, dists, bondpair, npbonds);
      if (deldist > epsdMax)
	{
	  /* se la variazione di d è eccessiva 
	   * cerca di correggere il passo per ottenere un valore
	   * più vicino a epsd*/
	  t -= delt;
	  //delt = d2old / maxddot;
	  delt = epsd/maxddot;
	  /* NOTE: prob. la seguente condizione si puo' rimuovere 
	   * o cambiare in > */
#if 0
	  deltth = h;
	  if (delt < deltth)
	    {
	      delt = deltth;
	    }
#endif
	  ///delt = h;///
	  t += delt; 
	  //t += delt*epsd/fabs(d2-d2old);
	  itsS++;
	  d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	}
#endif
#ifdef MD_NEGPAIRS
      itstb = 0;
      /* NOTA: se la distanza tra due sticky spheres è positiva a t (per errori numerici 
       * accade spesso) e t+delt allora delt è troppo grande e qui lo riduce fino ad un 
       * valore accettabile. */
      if (sumnegpairs)// && !firstaftsf)
	{
#if 0
	  while (fabs(dold-d)<1E-12)
	    {
	      delt *= GOLD;
	      t = tini + delt;
	      d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	    }
#endif
	  if (delt_is_too_big(i, j, bondpair, dists, distsOld, negpairs))
	    {
	      //printf("delt=%.15G\n", delt);
	      deltini = delt;
	      if (!interpol(i, j, sumnegpairs-1, t1, tini, delt, distsOld[sumnegpairs-1], 
			    dists[sumnegpairs-1], &tmin, shift, 1))
		{
		  //printf("qui\n");
		  tmin -= t1;
		  //printf("i=%d j=%d QUIIIIIIIIIIIIII delt=%.15G t1=%.15G\n", i,j,delt, t1);
		  delt = tmin - tini;
		  t = tmin;
		  //printf(">>>> QUI delt=%.15G tmin=%.15G tini=%.15G tmin-t=%.15G\n", delt, tmin, tini, tmin-t);
		  d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
		}
#if 1
	      /* if, for any reason, the minimum just found does not exist, lay outside [tini, tini+delt]
		 or breaks bond, try harder using follwing loop, alternatively we could use brent minimization
		 algorithm */
	      if (delt_is_too_big(i, j, bondpair, dists, distsOld, negpairs))
	      //else 
		{
		  delt = deltini;
		  t = tini;
		  printf("[INFO] using old goldenfactor method to reduce delt\n");
		  MD_DEBUG30(exit(-1));
		  /*printf("[INFO] using old goldenfactor method to reduce delt\n");
		  printf("tini=%.15G tmin=%.15G t+delt=%.15G sumnegpairs=%d delt=%.15G d1=%.15G d2=%.15G\n",
			 tini, tmin, tini+delt, sumnegpairs, delt, distsOld[sumnegpairs-1],
			 dists[sumnegpairs-1]);
		  exit(-1);*/	  
		  while (delt_is_too_big(i, j, bondpair, dists, distsOld, negpairs) && 
			 delt > minh)
		    {
		      delt /= GOLD; 
		      t = tini + delt;
		      d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
		      itstb++;
		      if (!interpol(i, j, sumnegpairs-1, t1, tini, delt, distsOld[sumnegpairs-1], 
				    dists[sumnegpairs-1], &tmin, shift, 1))
			{
			  tmin -= t1;
			  delt = tmin - tini;
			  t = tmin;
			  d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
			  break;
			}
		    }
		}
#endif
	    }
	  sumnegpairs = 0;
	}
#endif
      MD_DEBUG30(printf(">>>>> d = %.15G\n", d));
      for (nn=0; nn < npbonds; nn++)
	dorefine[nn] = MD_EVENT_NONE;
      //ncr=check_cross_strictcheck(i, j, distsOld, dists, crossed, bondpair);
      ncr=check_cross(distsOld, dists, crossed, bondpair, npbonds);
      /* N.B. crossed[] e tocheck[] sono array relativi agli 8 possibili tipi di attraversamento fra gli atomi
       * sticky */
      for (nn = 0; nn < npbonds; nn++)
	{
	  t2arr[nn] = t; 

	  dorefine[nn] = MD_EVENT_NONE;
	  if (crossed[nn]!=MD_EVENT_NONE)
	    {
	      /* se dorefine è 2 vuol dire che due superfici si sono
	       * attraversate */
	      if (valid_collision(i, j, mapbondsa[nn], mapbondsb[nn], crossed[nn]))
		{
		  MD_DEBUG30(printf("type: %d i=%d j=%d ata=%d atb=%d bound:%d\n", crossed[nn], i, j, mapbondsa[nn],
				    mapbondsb[nn], bound(i, j, mapbondsa[nn], mapbondsb[nn])));
#if 0
		  if (collCode==MD_INOUT_BARRIER && lastbump[i].mol==j && lastbump[j].mol==i
		      && lastbump[i].at == ata && lastbump[j].at==atb// && bound(i, j, ata, atb) 
		      && fabs(lastcol[i]-t)< 1E-14) 
		    {
		      //printf("qui\n");

		    }
		  else
#endif
		    dorefine[nn] = crossed[nn];
		  //printf("CROSSING dorefine[%d]:%d\n", nn, dorefine[nn]);
		}
	    }
	}

#define MD_INTERPOL
#ifdef MD_INTERPOL
      ntc = get_dists_tocheck(distsOld, dists, tocheck, dorefine, bondpair, npbonds);
      for (nn = 0; nn < npbonds; nn++)
	{
	  if (tocheck[nn])
	    {
	      //printf("tocheck[%d]:%d\n", nn, tocheck[nn]);
	      if ((retip=interpol(i, j, nn, t1, t-delt, delt, distsOld[nn], dists[nn], 
			   &troot, shift, 0)))
		{
		  if (retip==1)
		    dorefine[nn] = MD_EVENT_NONE;
		  else
		    {
		      /* interpolSP ha trovato un minimo ma non c'è stato cambio di segno */
		      if (grazing_try_harder(i, j, nn, t1, t-delt, delt, distsOld[nn], dists[nn], shift, &troot, &dmin))
			{
			  if (distsOld[nn] > 0.0)
			    dorefine[nn] = MD_OUTIN_BARRIER;
			  else
			    dorefine[nn] = MD_INOUT_BARRIER;
			  if (!valid_collision(i, j, mapbondsa[nn], mapbondsb[nn], crossed[nn]))
			    dorefine[nn] = MD_EVENT_NONE;
			  else
			    t2arr[nn] = troot - t1;
			}
		      else
			dorefine[nn] = MD_EVENT_NONE;
		    }

		}
	      else 
		{
		  if (distsOld[nn] > 0.0)
		    dorefine[nn] = MD_OUTIN_BARRIER;
		  else
		    dorefine[nn] = MD_INOUT_BARRIER;
		  if (!valid_collision(i, j, mapbondsa[nn], mapbondsb[nn], crossed[nn]))
		    dorefine[nn] = MD_EVENT_NONE;
		  else
		    t2arr[nn] = troot - t1;
		}
	    }
	}
#endif
      tmin = 0;
      gotcoll = 0;
      for (nn = 0; nn < npbonds; nn++)
	{
	  if (dorefine[nn]!=MD_EVENT_NONE)
	    {
	      MD_DEBUG30(printf("REFINE dorefine[%d]:%d\n", nn, dorefine[nn]));
	      MD_DEBUG30(printf("t-delt=%.15G t2arr[%d]=%.15G\n", t-delt, nn, t2arr[nn]));
	      if (refine_contact(i, j, t1, t-delt, t2arr[nn], nn, shift, &troot))
		{
		  //printf("[locate_contact] Adding collision between %d-%d\n", i, j);
		  MD_DEBUG30(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
		  MD_DEBUG30(printf("[locate_contact] t=%.15G nn=%d\n", t, nn));
		  MD_DEBUG(printf("[locate_contact] its: %d\n", its));
		  /* se il legame già c'è e con l'urto si forma tale legame allora
		   * scarta tale urto */
		  if ((adjt1 && troot < Oparams.time) || troot > t2 || troot < t1)
#if 0
		      || 
		      (lastbump[i].mol == j && lastbump[j].mol==i && 
		       lastbump[i].at == mapbondsa[nn]
		       && lastbump[j].at == mapbondsb[nn] && fabs(troot - lastcol[i]) < 1E-14))
#endif
		    {
#if 0
		      if (lastbump[i].mol == j && lastbump[j].mol==i && 
			  lastbump[i].at == mapbondsa[nn]
			  && lastbump[j].at == mapbondsb[nn] && fabs(troot - lastcol[i]) < 1E-14)
			{
			  printf("dold[%d]:%.15G d[%d]:%.15G\n", nn, distsOld[nn], nn, dists[nn]);
			  printf("state: %d collision type: %d\n", 
				 bound(i,j,mapbondsa[nn],mapbondsb[nn]), dorefine[nn]);
			  printf("lastcol[%d]: %.15G troot: %.15G\n", i, lastcol[i], troot);
			  printf("fabs(lastcol[i]-troot):%.15G\n",  fabs(troot - lastcol[i]));
			  printf("t1: %.15G t2: %.15G\n", t1, t2);
			  printf("delt: %.15G\n", delt);
			}
#endif
		      //gotcoll = -1;
		      continue;
		    }
		  else
		    {
		      if (*collCode == MD_EVENT_NONE || troot < *evtime)
			{
			  gotcoll = 1;
			  *ata = mapbondsa[nn];
			  *atb = mapbondsb[nn];
			  *evtime = troot;
			  *collCode = dorefine[nn]; 
			}
		      continue;
		    }
		}
	      else 
		{
		  MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
#ifdef MD_INTERPOL
		  if (!tocheck[nn])
#endif
		    mdPrintf(ALL,"[locate_contact] can't find contact point!\n",NULL);
		  /* Se refine_contact fallisce deve cmq continuare a cercare 
		   * non ha senso smettere...almeno credo */
		  //gotcoll = -1;
		  continue;
#if 0
		  if (dorefine[nn] == 2)
		    {
		      MD_DEBUG10(printf("t=%.15G d2 < 0 and I did not find contact point, boh...\n",t));
		      MD_DEBUG10(printf("d1: %.15G d2: %.15G\n", d1, d2));
		      MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
		      gotcoll = -1;
		      continue;
		    }
		  else
		    {
		      printf("dorefine[%d]:%d\n", nn, dorefine[nn]);
		      printf("d: %.15G dold: %.15G\n", d, dold); 
		    }
#endif
		}
	    }
	}
      if (gotcoll == 1)
	return 1;
      else if (gotcoll == -1)
	return 0;
      if (fabs(d) > epsdFastR)
	{
	  if (search_contact_faster(i, j, shift, &t, t1, t2, epsd, &d, epsdFast, dists, bondpair,
				    maxddot, maxddoti))
	    {
	      MD_DEBUG30(printf("[search contact faster locate_contact] d: %.15G\n", d));
	      return 0;
	    }
#if 1
	  dold = d;
	  assign_dists(dists, distsOld);
#else
	  dold = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
#endif
	  firstaftsf = 1;
	  its++;
	  //itsS++;
	  continue;
	}
      dold = d;
      MD_DEBUG30(printf("==========>>>>> t=%.15G t2=%.15G\n", t, t2));
#ifndef MD_BASIC_DT
      assign_dists(distsOld,  distsOld2);
#endif
      assign_dists(dists, distsOld);
      its++;
      itsS++;
    }
  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return 0;
}

#define EPS 1e-4
#if 0
double estimate_tmin(double t, int na, int nb)
{
  int cellRangeT[2 * NDIM];
  int iX, iY, iZ, jX, jY, jZ, k, n;
  double te, d, tmin=-1, shiftd[3], shift[3], r1[3], r2[3], alpha, vecg[8], maxddot;

  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
#ifndef MD_GRAVITY
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
      if (jZ == -1) 
	{
	  jZ = cellsz - 1;    
	  shift[2] = - L;
	} 
      else if (jZ == cellsz) 
	{
	  jZ = 0;    
	  shift[2] = L;
	}
#endif
      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	{
	  jY = inCell[1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cellsy - 1;    
	      shift[1] = -L;
	    } 
	  else if (jY == cellsy) 
	    {
	      jY = 0;    
	      shift[1] = L;
	    }
	  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
	    {
	      jX = inCell[0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cellsx - 1;    
		  shift[0] = - L;
		} 
	      else if (jX == cellsx) 
		{
		  jX = 0;   
		  shift[0] = L;
		}
	      n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	      for (n = cellList[n]; n > -1; n = cellList[n]) 
		{
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
		      d = calcDistNeg(t, na, n, shift, r1, r2, &alpha, vecg, 1);
		      maxddot = sqrt(Sqr(vx[na]-vx[n])+Sqr(vy[na]-vy[n])+Sqr(vz[na]-vz[n])) +
			sqrt(Sqr(wx[na])+Sqr(wy[na])+Sqr(wz[na]))*maxax[na]
			+ sqrt(Sqr(wx[n])+Sqr(wy[n])+Sqr(wz[n]))*maxax[n];

		      if (d>0)// && n!=lastbump[na] && lastbump[n]!=na)
			{
			  te = t + d / maxddot;
			  if (tmin==-1 ||t < tmin)
			    tmin = te;
			}
		    }
		} 
	    }
	}
    }
  return tmin;
}
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
//extern int currentIndex;
void PredictCellCross(int na, int nc)
{
  int ignorecross[3], k, evCode, signDir[NDIM]={0,0,0}, iA, nl;
  double tm[3], sigma;
#ifdef MD_GRAVITY
  double Lzx, cells[NDIM]; 
  double h1, h2, hh1;
#endif

  if (nc==1 && (Oparams.parnumA==Oparams.parnum || Oparams.parnumA==0))
    return;
  iA = (na < Oparams.parnumA)?0:1;
  ignorecross[0] = ignorecross[1] = ignorecross[2] = 1;
  if (iA == 0 && nc == 0)
    nl = 0;
  else if (iA == 0 && nc == 1)
    nl = 2;
  else if (iA == 1 && nc == 0)
    nl = 1;
  else
    nl = 3;
#if 0
  printf("[PredictCellCross ]time=%f nl=%d nc=%d n=%d inCell: %d %d %d cells: %d %d %d\n",
	 Oparams.time, nl, nc, na, inCell[nc][0][na], inCell[nc][1][na], inCell[nc][2][na],
	 cellsx[nl], cellsy[nl], cellsz[nl]);
#endif
#ifdef MD_GRAVITY
  Lzx = OprogStatus.extraLz + Lz;
  //ignorecross[2]=0;
  /* NOTA: Il muro inferiore è posto a Lz / 2 cioè sul fondo del box,
   * in questo modo la base di ogni sfera poggia esattamente 
   * sul fondo della scatola */
  if (nc==0 && inCell[nc][2][na] == 0)
    {
      hh1 =  vz[na] * vz[na] + 2.0 * Oparams.ggrav * (rz[na] + Lz2);
      if (OprogStatus.targetPhi > 0.0)
	{
	  sigma = 2.0*axa[na];
	}
      else
	{
    	  if (na < parnumA)
	    sigma = Oparams.sigma[0][0];
	  else if (na >= parnumA)
	    sigma = Oparams.sigma[1][1];
	}
      h1 = hh1 -  Oparams.ggrav * sigma;
    }
  else
    {
      h1 = hh1 = vz[na] * vz[na] + 2.0 * Oparams.ggrav *
	(rz[na] + Lz2 - inCell[nc][2][na] * (Lzx) / cellsz[nl]);
    }
  if (vz[na] > 0.0) 
    {
      /* h1 è il Discriminante dell'intersezione con la faccia 
       * inferiore della cella lungo z, h2 di quella superiore,
       * per cui se vz > 0 e h2 > 0 allora si la particella attraverserà
       * la faccia superiore poiché h2 < h1 => t_h2 < t_h1
       * (si noti che la soluzione con tempo positivo se vz > 0 è quella
       * con il + nella formula per la risoluz. di un'eq. di secondo grado */
      h2 = hh1 - 2.0 * Oparams.ggrav * Lzx  / cellsz[nl];
      if (h2 > 0.0) 
	{
       	  tm[2] =  (vz[na] - sqrt (h2)) / Oparams.ggrav;
	  signDir[2] = 0;/* signDir = 0 vuol dire che la direzione è
			    positiva (faccia superiore in questo caso) */
	  if (nc==1 && inCell[nc][2][na]==cellsz[nl]-1)
	    {
	      tm[2] = timbig;
	      ignorecross[2] = 1;
	    }
	  else
	    ignorecross[2] = 0;
	} 
      else 
	{
	  tm[2] =  (vz[na] + sqrt (h1)) / Oparams.ggrav;
	  signDir[2] = 1;/* direzione negativa (faccia inferiore in questo caso) */
	  if (nc==1 && inCell[nc][2][na]==0)
	    {
	      tm[2] = timbig;
	      ignorecross[2] = 1;
	    }
	  else
	    ignorecross[2] = 0;

	}
    } 
  else 
    {
      tm[2] =  (vz[na] + sqrt (h1)) / Oparams.ggrav;
      signDir[2] = 1;
      if (nc==1 && inCell[nc][2][na]==0)
	{
	  ignorecross[2] = 1;
	  tm[2] = timbig;
	}
      else
	ignorecross[2] = 0;
    }
#else
  if (vz[na] != 0.0) 
    {
      if (vz[na] > 0.0) 
	{
	  signDir[2] = 0;/* direzione positiva */
	  if (nc==1 && inCell[nc][2][na]==cellsz[nl]-1)
	    ignorecross[2] = 1;
	  else
	    ignorecross[2] = 0;
	}
      else 
	{
	  signDir[2] = 1;/* direzione negativa */
	  if (nc==1 && inCell[nc][2][na]==0)
	    ignorecross[2] = 1;
	  else
	    ignorecross[2] = 0;
	}
      if (ignorecross[2])
	tm[2] = timbig;
      else
	tm[2] = ((inCell[nc][2][na] + 1 - signDir[2]) * L /
  		 cellsz[nl] - rz[na] - L2) / vz[na];
    } 
  else 
    tm[2] = timbig;
  /* end forcefield[k] != 0*/
#endif
  if (vx[na] != 0.0) 
    {
      if (vx[na] > 0.0) 
	{
	  if (nc==1 && inCell[nc][0][na]==cellsx[nl]-1)
	    ignorecross[0] = 1;
	  else
	    ignorecross[0] = 0;
	  signDir[0] = 0;/* direzione positiva */
	}
      else
	{
	  if (nc==1 && inCell[nc][0][na]==0)
	    ignorecross[0] = 1;
	  else
	    ignorecross[0] = 0;
	  signDir[0] = 1;/* direzione negativa */
	}
      if (ignorecross[0])
	tm[0] = timbig;
      else
	tm[0] = ((inCell[nc][0][na] + 1 - signDir[0]) * L /
	      	 cellsx[nl] - rx[na] - L2) / vx[na];
    } 
  else 
    tm[0] = timbig;

  if (vy[na] != 0.) 
    {
      if (vy[na] > 0.) 
	{
	  if (nc==1 && inCell[nc][1][na]==cellsy[nl]-1)
	    ignorecross[1] = 1;
	  else
	    ignorecross[1] = 0;
	  signDir[1] = 0;
	}
      else
	{
	  if (nc==1 && inCell[nc][1][na]==0)
	    ignorecross[1] = 1;
	  else
	    ignorecross[1] = 0;
	  signDir[1] = 1;
	}
      if (ignorecross[1])
	tm[1] = timbig;
      else
	tm[1] = ((inCell[nc][1][na] + 1 - signDir[1]) * L /
	      	 cellsy[nl] - ry[na] - L2) / vy[na];
    } 
  else 
    tm[1] = timbig;
  /* ====== */
  /* Find minimum time */
  k = -1; /* giusto per dare un valore ed evitare una warning */
  if (tm[1] <= tm[2]) 
    {
      if (tm[0] <= tm[1]) k = 0;
      else k = 1;
    } 
  else
    {
      if (tm[0] <= tm[2]) 
	k = 0;
      else 
	k = 2;
    }
  /* Se un errore numerico fa si che tm[k] < 0 allora lo poniamo uguale a 0
   * (ved. articolo Lubachevsky) */
#if 1
  if (tm[k]<0.0)
    {
      printf(">>>[PRED. CELL CROSS NEG. TIME] tm[%d]: %.15G time=%.15G\n", k, tm[k], Oparams.time);
      tm[k] = 0.0;
      printf("Cells(%d,%d,%d)\n", inCell[nc][0][na], inCell[nc][1][na], 
	     inCell[nc][2][na]);
#if 0
#ifdef MD_GRAVITY
      printf("na=%d real cells: %d %d %d\n", na, (int)((rx[na] + L2) * cellsx[nl] / L),
	     (int)((ry[na] + L2) * cellsy[nl] / L), (int)((rz[na] + Lz2)  * cellsz[nl] / (Lz+OprogStatus.extraLz)));
      printf("vz[na]=%.15G pos=%f %f %f h1=%.15G sqrt(h1)=%.15G relpos=%.15G h2=%.15G Lzx=%.15G\n", vz[na], rx[na], ry[na], rz[na], h1, sqrt(h1), rz[na]-0.5+Lz2, h2, Lzx);

      printf("<<< ===== \n");
#else 
      printf("real cells: %d %d %d\n", (int)((rx[na] + L2) * cellsx[nl] / L),
	     (int)((ry[na] + L2) * cellsy[nl] / L), (int)((rz[na] + L2)  * cellsz[nl] / L));
#endif
#endif
#if 0
      printf("idA=%d idB=%d treeQIndex[%d]=%d treeStatus[]=%d\n ", treeIdA[na+1], treeIdB[na+1], na+1, treeQIndex[na+1], treeStatus[na+1]);
      printf("currentIndex=%d\n", OprogStatus.curIndex);
      printf("nc=%d na=%d nl=%d\n",nc,na,nl);
      printf("tm[%d]<0 step %lld na=%d\n", k, (long long int)Oparams.curStep, na);
      printf("Cells(%d,%d,%d)\n", inCell[nc][0][na], inCell[nc][1][na], 
	     inCell[nc][2][na]);
      printf("cells= (%d,%d,%d)\n", cellsx[nl], cellsy[nl], cellsz[nl]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
      exit(-1);
      /*tm[k] = 0.0;*/
#endif
    }
#endif
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;// + 3*nc;
#ifdef MD_GRAVITY
  if (k == 2 && nc == 0 && inCell[0][2][na] == 0 && signDir[2] == 1) 
    {
      evCode = 4;/* sarebbe 2*k con k=2 (z) e per me vuol dire urto con parete in basso
		    che è anche l'unica nel presente caso */
      MD_DEBUG2(printf("wall!!! (evIdA: %d)\n", na));
#if 0
      if (na==24)
	{
	  printf("[PRED COLL WALL ]coll wall na=%d time=%.15G\n", na, Oparams.time+tm[k]);
	  printf("<<< ==== \n");
	}
#endif
    }
#endif
#if 0
  if (na==24)
    {
      printf("[PRED CELL CROSS] coll wall na=%d time=%.15G\n", na, Oparams.time+tm[k]);
      printf("<<< ==== \n");
    }
#endif
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 
  MD_DEBUG15(printf("schedule event [WallCrossing](%d,%d) tm[%d]: %.8G\n", 
		    na, ATOM_LIMIT+evCode, k, tm[k]));

  if (!ignorecross[k])
    {
#if 0
      if (na==24 && nc==1 && inCell[nc][2][na]==0 && signDir[2]==1 && k==2)
	printf("????????ma che cazzo succede? evCode=%d k=%d\n", evCode, k);
#endif	
      ScheduleEventBarr (na, ATOM_LIMIT + evCode, nc, 0, MD_EVENT_NONE, Oparams.time + tm[k]);
    }
  //printf("===>crossevtodel[%d]:%d\n", na, crossevtodel[na]);
  //printf("schedule event [WallCrossing](%d,%d) tm[%d]: %.16G time=%.15G evCode:%d\n", 
//	 na, ATOM_LIMIT+evCode, k, tm[k], tm[k]+Oparams.time, evCode);

}
int sticky_bump(int n, int na, int nl)
{
#ifdef MD_SPOT_OFF
  return 0;
#endif
#ifdef MD_GROWTH_CODE
  /* gli spot sono disattivi durante la crescita */	
  if (OprogStatus.targetPhi > 0.0)
    return 0;
#endif
#ifdef MD_THREESPOTS
  return 1;
#elif defined(MD_AB41)
  /* le particelle B (nl=1 i.e. interazione B-B) interagiscono solo come HS*/
  if (nl==0||nl==2||nl==3)
    return 1;
  else
    return 0;
#else
  if (nl==2||nl==3)
    return 1;
  else 
    return 0;
#endif
}
void PredictColl (int na, int nb, int nl) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq=0.0, dr[NDIM], dv[NDIM], shift[NDIM],  
	 b, d, t, tInt, vv, distSq, t1=0.0, t2=0.0, evtime=0, evtimeHC;
  int overlap, ac, bc, acHC, collCodeOld, nc;
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cells[NDIM];*/
  int collCode;
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
#if 0
  int evCode;
#endif
  MD_DEBUG29(printf("PredictEvent: %d,%d\n", na, nb));
  MD_DEBUG(calc_energy("PredEv"));
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  if (nl < 2)
    {
      nc = 0;
    }
  else
    {
      nc = 1;
    }
  /*
  printf("[PredictEvent ]nl=%d nc=%d n=%d inCell: %d %d %d cells: %d %d %d\n",
	 nl, nc, na, inCell[nc][0][na], inCell[nc][1][na], inCell[nc][2][na],
	 cellsx[nl], cellsy[nl], cellsz[nl]);
  */
    /* NOTA: le linked list sono tre:
   *  0 = lista dell'interazione AA
   *  1 = lista dell'interazione BB
   *  2 = lista dell'interazione AB 
   *  inoltre inCell[ncel][dir][na] contiene 
   *  la cella a cui appartiene la particella A 
   *  Per ogni particella ci sono due insiemi di celle
   *  ncel = 0 celle relative all'interazione della particella na
   *  con particelle della stessa specie.
   *  ncel = 1 celle relative all'interazione della particella na 
   *  con particelle della stessa specie. */
 
#if 0
  {
      int ii,nn,np;    
      for (np=0; np < Oparams.parnum; np++) 
	{
	  nn = (inCell[0][2][np] *cellsy[1] + inCell[0][1][np]) * cellsx[1] + inCell[0][0][np] + Oparams.parnum;
	  if (cellList[1][nn]!=-1)
	    { printf("lista nl=1 non vuota per n=%d\n",n); exit(-1);}
	  nn = (inCell[1][2][np] *cellsy[3] + inCell[1][1][np]) * cellsx[3] + inCell[1][0][np] + Oparams.parnum;
	  if (cellList[3][nn]!=-1)
	    {printf("lista nl=3 non vuota per n=%d\n",n); exit(-1);}
	}
  }
#endif
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_GRAVITY
      /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (inCell[nc][2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
  if (inCell[nc][2][na] + cellRangeT[2 * 2 + 1] == cellsz[nl]) cellRangeT[2 * 2 + 1] = 0;
#endif
 
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[nc][2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
#ifndef MD_GRAVITY
      /* no periodic boundary conditions under gravity along z axis */
      if (jZ == -1) 
	{
	  jZ = cellsz[nl] - 1;    
	  shift[2] = - L;
	} 
      else if (jZ == cellsz[nl]) 
	{
	  jZ = 0;    
	  shift[2] = L;
	}
#endif
      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	{
	  jY = inCell[nc][1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cellsy[nl] - 1;    
	      shift[1] = -L;
	    } 
	  else if (jY == cellsy[nl]) 
	    {
	      jY = 0;    
	      shift[1] = L;
	    }
	  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
	    {
	      jX = inCell[nc][0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cellsx[nl] - 1;    
		  shift[0] = - L;
		} 
	      else if (jX == cellsx[nl]) 
		{
		  jX = 0;   
		  shift[0] = L;
		}
	      n = (jZ *cellsy[nl] + jY) * cellsx[nl] + jX + Oparams.parnum;
	      for (n = cellList[nl][n]; n > -1; n = cellList[nl][n]) 
		{
		  //printf("nl=%d cellList[nl=%d][n=%d]:%d na=%d\n", nl, nl, n, cellList[nl][n], na);
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
		      collCode = MD_EVENT_NONE;
		      tInt = Oparams.time - atomTime[n];
	    	      dr[0] = rx[na] - (rx[n] + vx[n] * tInt) - shift[0];	  
    		      dv[0] = vx[na] - vx[n];
		      dr[1] = ry[na] - (ry[n] + vy[n] * tInt) - shift[1];
		      dv[1] = vy[na] - vy[n];
#ifdef MD_GRAVITY
		      dr[2] = rz[na] - 
	    		(rz[n] + (vz[n] - 0.5 * Oparams.ggrav * tInt) * tInt) - shift[2];
    		      dv[2] = vz[na] - (vz[n] - Oparams.ggrav * tInt);

#else
		      dr[2] = rz[na] - (rz[n] + vz[n] * tInt) - shift[2];
		      dv[2] = vz[na] - vz[n];
#endif		      
	    	      b = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
#ifndef MD_SPOT_OFF
    		      distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
    		      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
#endif
		      /* N.B. 2 e 3 sono le liste per urti tra specie diverse */
		      if (sticky_bump(n,na,nl))
			{
			  /* NOTA 040510: notare che qui non c'entra se stiamo facendo una crescita 
			     poiché gli spot vengono disattivati durante la crescita */	
			  /* maxax[...] è il diametro dei centroidi dei due tipi
			   * di ellissoidi */
		      	  if (na < parnumA && n < parnumA)
	    		    sigSq = Sqr(maxax[na]);
    			  else if (na >= parnumA && n >= parnumA)
			    sigSq = Sqr(maxax[na]);
			  else
			    sigSq = Sqr((maxax[n]+maxax[na])*0.5);
			  MD_DEBUG2(printf("sigSq: %f\n", sigSq));
 			  d = Sqr (b) - vv * (distSq - sigSq);

			  collCode = MD_EVENT_NONE;
			  if (d < 0 || (b > 0.0 && distSq > sigSq)) 
			    {
			      /* i centroidi non collidono per cui non ci può essere
			       * nessun urto sotto tali condizioni */
			      continue;
			    }
			  MD_DEBUG(printf("PREDICTING na=%d n=%d\n", na , n));
			  if (vv==0.0)
			    {
			      if (distSq >= sigSq)
				continue;
			      /* la vel relativa è zero e i centroidi non si overlappano quindi
			       * non si possono urtare! */
			      t1 = t = 0;
			      t2 = 10.0;/* anche se sono fermi l'uno rispetto all'altro possono 
					   urtare ruotando */
			    }
			  else if (distSq >= sigSq)
			    {
			      t = t1 = - (sqrt (d) + b) / vv;
			      t2 = (sqrt (d) - b) / vv;
			      MD_DEBUG29(printf("NOT OVERLAP t1=%.15G t2=%.15G\n", t1, t2));
			      overlap = 0;
			    }
			  else 
			    {
			      MD_DEBUG29(printf("Centroids overlap!\n"));
			      t2 = t = (sqrt (d) - b) / vv;
			      t1 = 0;//-OprogStatus.h;
			      overlap = 1;
			      MD_DEBUG(printf("altro d=%f t=%.15f\n", d, (-sqrt (d) - b) / vv));
			      MD_DEBUG(printf("vv=%f dv[0]:%f\n", vv, dv[0]));
			    }
			  //printf("na=%d j=%d type=%d t1=%.15G\n", na, lastbump[na].mol,
			  //     lastbump[na].type, 
			  //   t1);
			  MD_DEBUG(printf("t=%f curtime: %f b=%f d=%f\n", t, Oparams.time, b ,d));
			  MD_DEBUG(printf("dr=(%f,%f,%f) sigSq: %f", dr[0], dr[1], dr[2], sigSq));
			  //t += Oparams.time; 
			  t2 += Oparams.time;
			  t1 += Oparams.time;
			  //printf("t1=%.15G t2=%.15G\n",t1,t2);
			  /* calcola cmq l'urto fra le due core spheres */
			}
#ifdef MD_GROWTH_CODE
		      if (OprogStatus.targetPhi > 0.0)
			{
			  /* N.B. nel caso della crescita per ora si assume che il diametro misto sia additivo
			     cioè sigma_ij = (sigma_i + sigma_j)*0.5 per ogni possibile coppia i,j */	  
			  sigSq = Sqr(axa[na]+axa[n]);
			}
		      else
			{
			  if (na < parnumA && n < parnumA)
			    sigSq = Sqr(Oparams.sigma[0][0]);
			  else if (na >= parnumA && n >= parnumA)
			    sigSq = Sqr(Oparams.sigma[1][1]);
			  else
			    sigSq = Sqr(Oparams.sigma[0][1]);
			}
#else
		      if (na < parnumA && n < parnumA)
			sigSq = Sqr(Oparams.sigma[0][0]);
		      else if (na >= parnumA && n >= parnumA)
			sigSq = Sqr(Oparams.sigma[1][1]);
		      else
			sigSq = Sqr(Oparams.sigma[0][1]);
#endif
#if 0
		      if ((distSq=Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2])) < 0.99*sigSq)
	    		{
    			  printf("[WARNING] Significant overlap: %.15G (%d-%d)\n", distSq/sigSq, na, n);
			}
#endif

		      if (b < 0.0) 
			{
#ifdef MD_SPOT_OFF
			  distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
#if 0
			  if (distSq < 0.99*sigSq)
			    {
			      printf("[WARNING] Significant overlap: %.15G (%d-%d)\n", distSq/sigSq, na, n);
			    }
#endif
			  vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
#endif
#if 1
			  d = Sqr (b) - vv * (distSq - sigSq);
#else
			  d = Sqr (b) - vv * 
			    (Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]) - sigSq);
#endif
			  if (d >= 0.) 
			    {
			      t = - (sqrt (d) + b) / vv;
			      if (t < 0)
				{
#if 0
				  printf("time:%.15f tInt:%.15f\n", Oparams.time,
					 tInt);
				  printf("t = %.15G\n", t);
				  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
							      Sqr(dr[2]))-1.0 );
				  printf("STEP: %lld\n", (long long int)Oparams.curStep);
				  printf("atomTime: %.10f \n", atomTime[n]);
				  printf("n:%d na:%d\n", n, na);
				  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
				  t = 0;
				}
			      collCode = MD_CORE_BARRIER;
			      evtime = Oparams.time + t;
			      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
			      //printf("============>CORE COLLISION vv=%f d=%f evtime=%.15G\n", vv, d, evtime);
			    } 
			}
		      collCodeOld = collCode;
		      evtimeHC = evtime;
		      //printf("))))))))))evtime=%f collCode:%d\n", evtime, collCode);
		      acHC = ac = 0; 
		      acHC = bc = 0;
		      //calcDist(Oparams.time, na, n, shift, r1, r2);
		      //continue;
		      //exit(-1);
#if 1
		      /* gli sticky spots esistono solo nell'interazione AB */
		      if (sticky_bump(n,na,nl))
			{
			  if (!locate_contact(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			    {
			      if (collCode == MD_EVENT_NONE)
				continue;
			    }
			}
		      else if (collCode == MD_EVENT_NONE)
			continue;
#else
		      if (collCode == MD_EVENT_NONE)
			continue;
#endif
		      MD_DEBUG29(printf("evtime=%.15G evtimeHC=%.15G\n",evtime,evtimeHC));
#if 0
		      if (collCode!=MD_CORE_BARRIER && collCodeOld==MD_CORE_BARRIER &&
			  fabs(evtime - evtimeHC)<1E-5)
			{
			  ac = bc = 0;
			  evtime = evtimeHC;
			  collCode = collCodeOld;
			}
#endif
		      MD_DEBUG(printf("A x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)-B x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)",
				      rx[na], ry[na], rz[na], vx[na], vy[na], vz[na],
				      rx[n], ry[n], rz[n], vx[n], vy[n], vz[n]));
		      t = evtime;
#if 1
		      if (t < Oparams.time)
			{
#if 0
			  printf("time:%.15f tInt:%.15f\n", Oparams.time,
				 tInt);
			  printf("t = %.15G\n", t);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
						      Sqr(dr[2]))-1.0 );
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = Oparams.time;
			}
#endif
		      /* il tempo restituito da newt() è già un tempo assoluto */
		      MD_DEBUG30(printf("EVENT TIME=%.15G\n",evtime));
		      MD_DEBUG30(printf("time: %f Adding collision %d-%d\n", Oparams.time, na, n));
		      MD_DEBUG30(printf("t=%.15G %d-%d %d-%d ac=%d bc=%d type=%d\n",
					t, na,n,n,na,ac,bc, collCode));
		      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
		      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
		    }
		} 
	    }
	}
    }
}
#ifdef MD_GRAVITY
void ProcessCollWall(void)
{
  double sigma;
  int iA, nl_ignore, k, nc, nl;
  /* Dissipation */
  MD_DEBUG2(printf("Collision with wall evIdA: %d vz: %.15f\n", evIdA, vz[evIdA]));
  UpdateAtom(evIdA);
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }

  vz[evIdA] = -vz[evIdA];
  lastcol[evIdA] = Oparams.time;
  iA = (evIdA<Oparams.parnumA)?0:1;
  nl_ignore = (evIdA<Oparams.parnumA)?1:0;
  for (nc = 0; nc < 2; nc++)
    {
      PredictCellCross(evIdA, nc);
    }
  for (nl = 0; nl < 4; nl++)
    {
      if (nl==nl_ignore || nl==iA+2)
	continue;
      PredictColl(evIdA, -1, nl);
    }
}

#endif

#else
#ifdef MD_GRAVITY
void PredictEvent (int na, int nb);
void ProcessCollWall(void)
{
  double sigma;
  int k;
  /* Dissipation */
  MD_DEBUG2(printf("Collision with wall evIdA: %d vz: %.15f\n", evIdA, vz[evIdA]));
  UpdateAtom(evIdA);
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  vz[evIdA] = -vz[evIdA];
  lastcol[evIdA] = Oparams.time;
  PredictEvent(evIdA, -1);
  PredictEvent(evIdB, evIdA);
}

#endif

void PredictEvent (int na, int nb) 
{
#ifdef MD_GRAVITY
  double Lzx, h1, h2, sig, hh1;
#endif
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq=0.0, dr[NDIM], dv[NDIM], shift[NDIM], tm[NDIM],  
	 b, d, t, tInt, vv, distSq, t1=0.0, t2=0.0, evtime=0, evtimeHC;
  int overlap, ac, bc, acHC, collCodeOld;
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cells[NDIM];*/
  int collCode;
  int cellRangeT[2 * NDIM], signDir[NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
#if !defined(MD_SILICA) || defined(MD_USE_SINGLE_LL)
  int evCode;
#endif
  MD_DEBUG29(printf("PredictEvent: %d,%d\n", na, nb));
  MD_DEBUG(calc_energy("PredEv"));
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
#ifdef MD_GRAVITY
  Lzx = OprogStatus.extraLz + Lz;
  if (na < parnumA)
    sig = Oparams.sigma[0][0];
  else
    sig = Oparams.sigma[1][1];
  /* NOTA: Il muro inferiore è posto a Lz / 2 cioè sul fondo del box,
   * in questo modo la base di ogni sfera poggia esattamente 
   * sul fondo della scatola */
  if (inCell[2][na] == 0)
    {
      hh1 =  vz[na] * vz[na] + 2.0 * Oparams.ggrav *
	(rz[na] + Lz2);
      h1 = hh1 -  Oparams.ggrav * sig;
    }
  else
    h1 = hh1 = vz[na] * vz[na] + 2.0 * Oparams.ggrav *
      (rz[na] + Lz2 - inCell[2][na] * (Lzx) / cellsz);
  if (vz[na] > 0.0) 
    {
      /* h1 è il Discriminante dell'intersezione con la faccia 
       * inferiore della cella lungo z, h2 di quella superiore,
       * per cui se vz > 0 e h2 > 0 allora si la particella attraverserà
       * la faccia superiore poiché h2 < h1 => t_h2 < t_h1
       * (si noti che la soluzione con tempo positivo se vz > 0 è quella
       * con il + nella formula per la risoluz. di un'eq. di secondo grado */
      h2 = hh1 - 2.0 * Oparams.ggrav * Lzx  / cellsz;
      if (h2 > 0.0) 
	{
	  tm[2] =  (vz[na] - sqrt (h2)) / Oparams.ggrav;
	  signDir[2] = 0;/* signDir = 0 vuol dire che la direzione è
	  		    positiva (faccia superiore in questo caso) */
	} 
      else 
	{
	  tm[2] =  (vz[na] + sqrt (h1)) / Oparams.ggrav;
	  signDir[2] = 1;/* direzione negativa (faccia inferiore in questo caso) */
	}
    } 
  else 
    {
      tm[2] =  (vz[na] + sqrt (h1)) / Oparams.ggrav;
      signDir[2] = 1;
    }

#else
  if (vz[na] != 0.0) 
    {
      if (vz[na] > 0.0) 
	signDir[2] = 0;/* direzione positiva */
      else 
	signDir[2] = 1;/* direzione negativa */
      tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L /
	       cellsz - rz[na] - L2) / vz[na];
    } 
  else 
    tm[2] = timbig;
  /* end forcefield[k] != 0*/
#endif
  if (vx[na] != 0.0) 
    {
      if (vx[na] > 0.0) 
	signDir[0] = 0;/* direzione positiva */
      else 
	signDir[0] = 1;/* direzione negativa */
      tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L /
	       cellsx - rx[na] - L2) / vx[na];
    } 
  else 
    tm[0] = timbig;

  if (vy[na] != 0.) 
    {
      if (vy[na] > 0.) 
	signDir[1] = 0;
      else 
	signDir[1] = 1;
      tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L /
	       cellsy - ry[na] - L2) / vy[na];
    } 
  else 
    tm[1] = timbig;
  /* ====== */
  /* Find minimum time */
  k = -1; /* giusto per dare un valore ed evitare una warning */
  if (tm[1] <= tm[2]) 
    {
      if (tm[0] <= tm[1]) k = 0;
      else k = 1;
    } 
  else
    {
      if (tm[0] <= tm[2]) 
	k = 0;
      else 
	k = 2;
    }
  /* Se un errore numerico fa si che tm[k] < 0 allora lo poniamo uguale a 0
   * (ved. articolo Lubachevsky) */
#if 1
  if (tm[k]<0)
    {
      tm[k] = 0.0;
#if 0
      printf("tm[%d]<0 step %lld na=%d\n", k, (long long int)Oparams.curStep, na);
      printf("Cells(%d,%d,%d)\n", inCell[0][na], inCell[1][na], inCell[2][na]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
      /*exit(-1);*/
      /*tm[k] = 0.0;*/
#endif
    }
#endif
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 
  MD_DEBUG15(printf("schedule event [WallCrossing](%d,%d) tm[%d]: %.8G\n", 
		    na, ATOM_LIMIT+evCode, k, tm[k]));

#ifdef MD_GRAVITY
  if (k == 2 && inCell[2][na] == 0 && signDir[2] == 1) 
    {
      evCode = 4;/* sarebbe 2*k con k=2 (z) e per me vuol dire urto con parete in basso
  		    che è anche l'unica nel presente caso */
      MD_DEBUG2(printf("wall!!! (evIdA: %d)\n", na));
    }
#endif

  ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]); 
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_GRAVITY
      /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (inCell[2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
  if (inCell[2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
#endif

#if 0
  /* no collisions at all! */
  return;
#endif
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
#ifndef MD_GRAVITY
      if (jZ == -1) 
	{
	  jZ = cellsz - 1;    
	  shift[2] = - L;
	} 
      else if (jZ == cellsz) 
	{
	  jZ = 0;    
	  shift[2] = L;
	}
#endif
      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	{
	  jY = inCell[1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cellsy - 1;    
	      shift[1] = -L;
	    } 
	  else if (jY == cellsy) 
	    {
	      jY = 0;    
	      shift[1] = L;
	    }
	  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
	    {
	      jX = inCell[0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cellsx - 1;    
		  shift[0] = - L;
		} 
	      else if (jX == cellsx) 
		{
		  jX = 0;   
		  shift[0] = L;
		}
	      n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	      for (n = cellList[n]; n > -1; n = cellList[n]) 
		{
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
		      collCode = MD_EVENT_NONE;
		      MD_DEBUG2(printf("sigSq: %f\n", sigSq));
		      tInt = Oparams.time - atomTime[n];
		      dr[0] = rx[na] - (rx[n] + vx[n] * tInt) - shift[0];	  
		      dv[0] = vx[na] - vx[n];
		      dr[1] = ry[na] - (ry[n] + vy[n] * tInt) - shift[1];
		      dv[1] = vy[na] - vy[n];
#ifdef MD_GRAVITY
		      dr[2] = rz[na] - 
			(rz[n] + (vz[n] - 0.5 * Oparams.ggrav * tInt) * tInt) - shift[2];
		      dv[2] = vz[na] - (vz[n] - Oparams.ggrav * tInt);
#else
		      dr[2] = rz[na] - (rz[n] + vz[n] * tInt) - shift[2];
		      dv[2] = vz[na] - vz[n];
#endif
		      b = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
#ifndef MD_SPOT_OFF
	    	      distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
    		      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
#endif
#ifndef MD_SPOT_OFF
		      /* N.B. 12/05/2010: se gli spot sono disattivi non ha senso usare i centroidi */
#ifdef MD_GROWTH_CODE 
	            if (OprogStatus.targetPhi <= 0.0)
			{
#endif
    			  /* maxax[...] è il diametro dei centroidi dei due tipi
    			   * di ellissoidi */
    			  if (na < parnumA && n < parnumA)
    			    sigSq = Sqr(maxax[na]);
    			  else if (na >= parnumA && n >= parnumA)
    			    sigSq = Sqr(maxax[na]);
    			  else
    			    sigSq = Sqr((maxax[n]+maxax[na])*0.5);

			  d = Sqr (b) - vv * (distSq - sigSq);
			  if (d < 0 || (b > 0.0 && distSq > sigSq)) 
			    {
			      /* i centroidi non collidono per cui non ci può essere
			       * nessun urto sotto tali condizioni */
			      continue;
			    }
			  MD_DEBUG(printf("PREDICTING na=%d n=%d\n", na , n));
			  if (vv==0.0)
			    {
			      if (distSq >= sigSq)
				continue;
			      /* la vel relativa è zero e i centroidi non si overlappano quindi
			       * non si possono urtare! */
			      t1 = t = 0;
			      t2 = 10.0;/* anche se sono fermi l'uno rispetto all'altro possono 
					   urtare ruotando */
			    }
			  else if (distSq >= sigSq)
			    {
			      t = t1 = - (sqrt (d) + b) / vv;
			      t2 = (sqrt (d) - b) / vv;
			      MD_DEBUG29(printf("NOT OVERLAP t1=%.15G t2=%.15G\n", t1, t2));
			      overlap = 0;
			    }
			  else 
			    {
			      MD_DEBUG29(printf("Centroids overlap!\n"));
			      t2 = t = (sqrt (d) - b) / vv;
			      t1 = 0;//-OprogStatus.h;
			      overlap = 1;
			      MD_DEBUG(printf("altro d=%f t=%.15f\n", d, (-sqrt (d) - b) / vv));
			      MD_DEBUG(printf("vv=%f dv[0]:%f\n", vv, dv[0]));
			    }
			  //printf("na=%d j=%d type=%d t1=%.15G\n", na, lastbump[na].mol,
			  //     lastbump[na].type, 
			  //   t1);
			  MD_DEBUG(printf("t=%f curtime: %f b=%f d=%f\n", t, Oparams.time, b ,d));
			  MD_DEBUG(printf("dr=(%f,%f,%f) sigSq: %f", dr[0], dr[1], dr[2], sigSq));
			  //t += Oparams.time; 
			  t2 += Oparams.time;
			  t1 += Oparams.time;
#ifdef MD_GROWTH_CODE
			}
#endif
#endif
		      //printf("t1=%.15G t2=%.15G\n",t1,t2);
		      /* calcola cmq l'urto fra le due core spheres */
#ifdef MD_GROWTH_CODE
		      if (OprogStatus.targetPhi > 0.0)
			{
			  /* N.B. nel caso della crescita per ora si assume che il diametro misto sia additivo
			     cioè sigma_ij = (sigma_i + sigma_j)*0.5 per ogni possibile coppia i,j */	  
			  sigSq = Sqr(axa[na]+axa[n]);
			}
		      else
			{
			  if (na < parnumA && n < parnumA)
			    sigSq = Sqr(Oparams.sigma[0][0]);
			  else if (na >= parnumA && n >= parnumA)
			    sigSq = Sqr(Oparams.sigma[1][1]);
			  else
			    sigSq = Sqr(Oparams.sigma[0][1]);
			}
#else
		      if (na < parnumA && n < parnumA)
			sigSq = Sqr(Oparams.sigma[0][0]);
		      else if (na >= parnumA && n >= parnumA)
			sigSq = Sqr(Oparams.sigma[1][1]);
		      else
			sigSq = Sqr(Oparams.sigma[0][1]);
#endif
		      if (b < 0.0) 
			{
#ifdef MD_SPOT_OFF
			  distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
			  vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
#endif
#if 1
			  d = Sqr (b) - vv * (distSq - sigSq);
#else
			  d = Sqr (b) - vv * 
			    (Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]) - sigSq);
#endif
			  if (d >= 0.) 
			    {
			      t = - (sqrt (d) + b) / vv;
			      if (t < 0)
				{
#if 0
				  printf("time:%.15f tInt:%.15f\n", Oparams.time,
					 tInt);
				  printf("t = %.15G\n", t);
				  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
							      Sqr(dr[2]))-1.0 );
				  printf("STEP: %lld\n", (long long int)Oparams.curStep);
				  printf("atomTime: %.10f \n", atomTime[n]);
				  printf("n:%d na:%d\n", n, na);
				  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
				  t = 0;
				}
			      collCode = MD_CORE_BARRIER;
			      evtime = Oparams.time + t;
			      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
			    } 
			}
		      collCodeOld = collCode;
		      evtimeHC = evtime;
		      //printf("))))))))))evtime=%f collCode:%d\n", evtime, collCode);
		      acHC = ac = 0; 
		      acHC = bc = 0;
		      //calcDist(Oparams.time, na, n, shift, r1, r2);
		      //continue;
		      //exit(-1);
#if 1
		      if (!locate_contact(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			{
			  if (collCode == MD_EVENT_NONE)
			    continue;
			}
#else
		      if (collCode == MD_EVENT_NONE)
			continue;
#endif
		      MD_DEBUG29(printf("evtime=%.15G evtimeHC=%.15G\n",evtime,evtimeHC));
#if 0
		      if (collCode!=MD_CORE_BARRIER && collCodeOld==MD_CORE_BARRIER &&
			  fabs(evtime - evtimeHC)<1E-5)
			{
			  ac = bc = 0;
			  evtime = evtimeHC;
			  collCode = collCodeOld;
			}
#endif
		      MD_DEBUG(printf("A x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)-B x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)",
				      rx[na], ry[na], rz[na], vx[na], vy[na], vz[na],
				      rx[n], ry[n], rz[n], vx[n], vy[n], vz[n]));
		      t = evtime;
#if 1
		      if (t < Oparams.time)
			{
#if 0
			  printf("time:%.15f tInt:%.15f\n", Oparams.time,
				 tInt);
			  printf("t = %.15G\n", t);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
						      Sqr(dr[2]))-1.0 );
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = Oparams.time;
			}
#endif
		      /* il tempo restituito da newt() è già un tempo assoluto */
		      MD_DEBUG30(printf("EVENT TIME=%.15G\n",evtime));
		      MD_DEBUG30(printf("time: %f Adding collision %d-%d\n", Oparams.time, na, n));
		      MD_DEBUG30(printf("t=%.15G %d-%d %d-%d ac=%d bc=%d type=%d\n",
					t, na,n,n,na,ac,bc, collCode));
		      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
		      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
		    }
		} 
	    }
	}
    }
}
#endif
void calc_energynew(char *msg)
{
  int i, k1;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#endif
  K = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (i<Oparams.parnumA)
	{
	  /* calcola tensore d'inerzia e le matrici delle due quadriche */
	  K += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    {
#ifdef MD_ASYM_ITENS
	      K += wt[k1]*ItensD[0][k1]*wt[k1];
#else
	      K += Oparams.I[0]*Sqr(wt[k1]);
#endif
	    }
	}
      else
	{
	  K += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    {
#ifdef MD_ASYM_ITENS
	      K += wt[k1]*ItensD[1][k1]*wt[k1];
#else
	      K += Oparams.I[1]*Sqr(wt[k1]);
#endif
	    }
	}
    }
  K *= 0.5;
#ifdef MD_ASYM_ITENS
  free_matrix(Ia,3);
  free_matrix(Ib,3);
#endif
  //printf("[%s] Kinetic Energy: %f\n", msg, K);
}
void calc_energy_i(char *msg, int i)
{
  int k1;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#endif
  K = 0;
  if (i<Oparams.parnumA)
    {
      /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
      RDiagtR(i, Ia, ItensD[0][0], ItensD[0][1], ItensD[0][2], R[i]);
#endif
      K += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
      wt[0] = wx[i];
      wt[1] = wy[i];
      wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    K += wt[k1]*Ia[k1][k2]*wt[k2];
	  }
#else
      for (k1=0; k1 < 3; k1++)
	K += Sqr(wt[k1])*Oparams.I[0];
#endif
    }
  else
    {
#ifdef MD_ASYM_ITENS
      RDiagtR(i, Ib, ItensD[1][0], ItensD[1][1], ItensD[1][2], R[i]);
#endif
      K += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
      wt[0] = wx[i];
      wt[1] = wy[i];
      wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    K += wt[k1]*Ib[k1][k2]*wt[k2];
	  }
#else
      for (k1=0; k1 < 3; k1++)
	K += Sqr(wt[k1])*Oparams.I[1];
#endif
    }
  K *= 0.5;
#ifdef MD_ASYM_ITENS
  free_matrix(Ia,3);
  free_matrix(Ib,3);
#endif
  //printf("[%s] Kinetic Energy of %d: %f\n", msg, i, K);

}

void calc_energy(char *msg)
{
  int i, k1;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#endif

  K = Ktra = Krot = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (i < Oparams.parnumA)
	{
	  /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
	  RDiagtR(i, Ia, ItensD[0][0], ItensD[0][1], ItensD[0][2], R[i]);
#endif
	  Ktra += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		K += wt[k1]*Ia[k1][k2]*wt[k2];
	      }
#else
	  for (k1=0; k1 < 3; k1++)
	    Krot += Sqr(wt[k1])*Oparams.I[0];
#endif
	}
      else
	{
#ifdef MD_ASYM_ITENS
	  RDiagtR(i, Ib, ItensD[1][0], ItensD[1][1], ItensD[1][2], R[i]);
#endif
	  Ktra += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		K += wt[k1]*Ib[k1][k2]*wt[k2];
	      }
#else
	  for (k1=0; k1 < 3; k1++)
	    Krot += Sqr(wt[k1])*Oparams.I[1];
#endif
	}
    }
  Ktra *= 0.5;
  Krot *= 0.5;
  K = Ktra + Krot;
#ifdef MD_ASYM_ITENS
  free_matrix(Ia,3);
  free_matrix(Ib,3);
#endif
  if (msg)
    printf("[%s] Kinetic Energy: %f\n", msg, K);
}

void store_bump(int i, int j)
{
  char fileop2[512], fileop[512];
  int ii,a;
  FILE *bf;
  double rat[5][3], rO[3], **Rl;
#if 0
  const char tipodat2[]= "%.15G %.15G %.15G\n";
#endif
#ifdef MD_BIG_DT
  sprintf(fileop2 ,"StoreBump-%d-%d-t%.8f", i, j, Oparams.time + OprogStatus.refTime);
#else
  sprintf(fileop2 ,"StoreBump-%d-%d-t%.8f", i, j, Oparams.time);
#endif
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
  R2u();
  fprintf(bf, ".Vol: %f\n", L*L*L);
  MD_DEBUG(printf("[Store bump]: %.15G\n", Oparams.time));
  Rl = matrix(3,3);
  for (ii = 0; ii < Oparams.parnum; ii++)
    {
      if (ii==i || ii==j)
	{
	  rO[0] = rx[ii];
	  rO[1] = ry[ii];
	  rO[2] = rz[ii];

	  Rl[0][0] = uxx[ii];
	  Rl[0][1] = uxy[ii];
	  Rl[0][2] = uxz[ii];
	  Rl[1][0] = uyx[ii];
	  Rl[1][1] = uyy[ii];
	  Rl[1][2] = uyz[ii];
	  Rl[2][0] = uzx[ii];
	  Rl[2][1] = uzy[ii];
	  Rl[2][2] = uzz[ii];

	  BuildAtomPos(ii, rO, Rl, rat);
	  /* write coords */
	  for (a = 0; a < 5; a++)
	    {
	      if (a == 0)
		{
		  fprintf(bf, "%.15G %.15G %.15G @ %f C[red]\n", rat[a][0], rat[a][1], rat[a][2],
			  Oparams.sigma[0][0]/2.0);
		}
	      else if (a < 3)
		{
#ifdef MD_AB41
		  if (ii < Oparams.parnumA)
		    fprintf(bf, "%.15G %.15G %.15G @ %f C[blue]\n", rat[a][0], rat[a][1], rat[a][2],
		  	    Oparams.sigmaStickyAA/2.0);
		  else
		    fprintf(bf, "%.15G %.15G %.15G @ %f C[blue]\n", rat[a][0], rat[a][1], rat[a][2],
			  Oparams.sigmaStickyAB/2.0);

#else
		  fprintf(bf, "%.15G %.15G %.15G @ %f C[blue]\n", rat[a][0], rat[a][1], rat[a][2],
			  Oparams.sigmaSticky/2.0);
#endif
		}
	      else if (a < 5)
		{
#ifdef MD_AB41
		  if (ii < Oparams.parnumA)
		    fprintf(bf, "%.15G %.15G %.15G @ %f C[green]\n", rat[a][0], rat[a][1], rat[a][2],
			    Oparams.sigmaStickyAA/2.0);
		  else
		    fprintf(bf, "%.15G %.15G %.15G @ %f C[green]\n", rat[a][0], rat[a][1], rat[a][2],
			    Oparams.sigmaStickyAB/2.0);

#else
		  fprintf(bf, "%.15G %.15G %.15G @ %f C[green]\n", rat[a][0], rat[a][1], rat[a][2],
			  Oparams.sigmaSticky/2.0);
#endif
		}
	    }
	}
    }
  //writeAllCor(bf);
  fprintf(bf,"%.15f %.15f %.15f @ 0.1 C[green]\n", rxC, ryC, rzC);
  fclose(bf);
  free_matrix(Rl,3);
}

extern void delete_events(int evIdA);
#ifdef MD_SURV_PROB
extern double *sp_firstcolltime, sp_start_time;
extern int *sp_has_collided, sp_equilib, *sp_coll_type;
extern int sp_tot_collisions;
extern void sp_reset_fct(void);
extern double gauss(void);
extern void save_sp(void);
extern double ranf(void);
void sp_update_cal(int ip)
{
#if 0
  if (ip!=evIdA && ip!=evIdB)
    PredictEvent(ip, -1);
#else
  UpdateSystem();
  rebuildCalendar();
  if (OprogStatus.intervalSum > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
#endif
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
}
#endif

void ProcessCollision(void)
{
  int k;
#ifdef MD_SURV_PROB
  int rebuilt_cal;
#endif

#ifdef MD_SILICA
  int nl, nl_ignore, iA, iB, nc;
#endif
  UpdateAtom(evIdA);
  UpdateAtom(evIdB);
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
#ifdef MD_SURV_PROB
  rebuilt_cal=0;
  if (!sp_equilib)  
    {
      int i,ip, ii, cur_trap=-1, cur_part=-1;
      double rTemp, mass;
#if 0
      if (Oparams.ntypes==2)
	{
	  /* Particelle di tipo 0 = assorbite
	     Particelle di tipo 1 = trappole, quindi se typeNP[1] == 1 allora TARGET problem
	     altrimenti se typeNP[0]=1 allora TRAPPING problem 
	     le TRAPPOLE sono IMMOBILI, dopo ogni assorbimento si fa andare la simulazione
	     per un tempo pari a OprogStatus.spdeltat per randomizzare il sistema */
	  if (typeNP[1]==1)
	    {
	      /* Una sola trappola fissa (le trappole sono fisse) => TARGET problem */
	      if (typeOfPart[evIdA]==1 || typeOfPart[evIdB]==1)
		{
	    	  save_sp();
		  sp_reset_fct();
		  //sp_start_time = Oparams.time;
		  sp_start_time = Oparams.time;
		  sp_equilib = 1;
		  if (typeOfPart[evIdA] == 1)
		    cur_trap=evIdA;
		  else
		    cur_trap=evIdB;
		  do
		    ip = (int) (Oparams.parnum*ranf());
		  while (typeOfPart[ip]==1 || ip==evIdA||ip==evIdB);
		 /* la particella scelta diventa una trappola e la trappola corrente
		    diventa una particella normale (cioè di tipo 0)*/ 
		  typeOfPart[ip] = 1;
		  typeOfPart[cur_trap] = 0;
		  if (ip!=evIdA && ip!=evIdB)
		    UpdateAtom(ip);
		  vx[ip] = 0.0; vy[ip]=0.0; vz[ip]=0.0;
		  sp_update_cal(ip);
		  rebuilt_cal=1;
		}
	    }
	  else
	    {
	      /* Oparams.parnum-1 trappole fisse => TRAPPING */
	      if (typeOfPart[evIdA]==0 || typeOfPart[evIdB]==0)
		{
		  save_sp();
		  sp_reset_fct();
		  //sp_start_time = Oparams.time;
		  if (typeOfPart[evIdA] == 0)
		    cur_part=evIdA;
		  else
		    cur_part=evIdB;
		  do
		    ip = (int) (Oparams.parnum*ranf());
		  while (ip==evIdA||ip==evIdB);
		 /* la particella scelta diventa una trappola e la trappola corrente
		    diventa una particella normale (cioè di tipo 0)*/ 
		  typeOfPart[ip] = 0;
		  typeOfPart[cur_part] = 1;
		  printf("sp_equilib=%d time=%.15G cur_part=%d ip=%d\n", sp_equilib, Oparams.time, cur_part, ip);
		  sp_start_time=Oparams.time;
		  sp_equilib = 1;
		  //if (ip!=evIdA && ip!=evIdB)
		    //UpdateAtom(ip);
		 //vx[cur_part] = 0.0; vy[cur_part]=0.0; vz[cur_part]=0.0;
		  typesArr[1].brownian = 1;
#if 0
		  for (i=0; i < Oparams.parnum; i++)
		    {
		      mass = typesArr[typeOfPart[i]].m;
		      rTemp = sqrt(Oparams.T / mass);  
		      vx[i] = rTemp * gauss(); 
		      vy[i] = rTemp * gauss();
		      vz[i] = rTemp * gauss();
		    }
		  sp_update_cal(cur_part);
		  rebuilt_cal=1;
#endif		  
		}
	    }
	}
#endif
      //printf("time=%f sp_start_time=%f sp_equilib=%d spdelta=%f\n", Oparams.time, sp_start_time, sp_equilib, OprogStatus.spdeltat);
      if (!sp_has_collided[evIdA])
	{
	  sp_firstcolltime[evIdA] = Oparams.time - sp_start_time;
	  /* sp_coll_type = 0->AA 1->BB 2->AB */
	  if (evIdA < Oparams.parnumA && evIdB < Oparams.parnumA)
	    {
	      sp_coll_type[evIdA] = 0;
	    }
	  else if (evIdA >= Oparams.parnumA && evIdB >= Oparams.parnumA)
	    {
	      sp_coll_type[evIdA] = 1;
	    }
	  else 
	    {
	      sp_coll_type[evIdA] = 2;
	    }

 	  sp_has_collided[evIdA] = 1;
	  sp_tot_collisions++;
	}
      if (!sp_has_collided[evIdB])
	{
	  sp_firstcolltime[evIdB] = Oparams.time - sp_start_time;
	  sp_has_collided[evIdB] = 1;
	  if (evIdA < Oparams.parnumA && evIdB < Oparams.parnumA)
	    {
	      sp_coll_type[evIdB] = 0;
	    }
	  else if (evIdA >= Oparams.parnumA && evIdB >= Oparams.parnumA)
	    {
	      sp_coll_type[evIdB] = 1;
	    }
	  else 
	    {
	      sp_coll_type[evIdB] = 2;
	    }
	  sp_tot_collisions++;
	}

      if (sp_tot_collisions == Oparams.parnum)
	{
	  save_sp();
	  sp_reset_fct();
	  sp_start_time = Oparams.time;
	  sp_equilib = 1;
	}
    }
  else 
    {
      double rTemp, mass;
      if (fabs(Oparams.time - sp_start_time) >= OprogStatus.spdeltat)
	{
	  int i;
  	  sp_equilib = 0;
	  sp_start_time = Oparams.time;
  	}
    }
#endif

  MD_DEBUG10(calc_energy("prima"));
  /* i primi due bit sono il tipo di event (uscita buca, entrata buca, collisione con core 
   * mentre nei bit restanti c'e' la particella con cui tale evento e' avvenuto */
  //printf("evIdC: %d evIdD: %d\n", evIdC, evIdD);
  bump(evIdA, evIdB, evIdC, evIdD, &W, evIdE);
  MD_DEBUG10(calc_energy("dopo"));
  MD_DEBUG(store_bump(evIdA, evIdB));
  //ENDSIM=1;
  /*printf("qui time: %.15f\n", Oparams.time);*/
#ifdef MD_ROTDIFF_MIS
  OprogStatus.lastcolltime[evIdA] =  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
  OprogStatus.lastcolltime[evIdB] = lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
#else
  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
#endif
  lastbump[evIdA].mol = evIdB;
  lastbump[evIdB].mol = evIdA;
  lastbump[evIdA].at = evIdC;
  lastbump[evIdB].at = evIdD;
  lastbump[evIdA].type = evIdE;
  lastbump[evIdB].type = evIdE;
#ifdef MD_HSVISCO
  OprogStatus.lastcoll = Oparams.time;
#endif
#ifdef MD_SURV_PROB
  if (rebuilt_cal)
    {
      do_check_negpairs=0;
      return;
    }
#endif 
  do_check_negpairs = 1;
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  iA = (evIdA<Oparams.parnumA)?0:1;
  nl_ignore = (evIdA<Oparams.parnumA)?1:0;
  for (nc = 0; nc < 2; nc++)
    {
      PredictCellCross(evIdA, nc);
      PredictCellCross(evIdB, nc);
    }
  for (nl = 0; nl < 4; nl++)
    {
      if (nl==nl_ignore || nl==iA+2)
	continue;
      PredictColl(evIdA, -1, nl);
    }
  iB = (evIdB<Oparams.parnumA)?0:1;
  nl_ignore = (evIdB<Oparams.parnumA)?1:0;
  for (nl = 0; nl < 4; nl++)
    {
      if (nl==nl_ignore || nl==iB+2)
	continue;
      PredictColl(evIdB, evIdA, nl);
    }
  
#else
  PredictEvent(evIdA, -1);
  PredictEvent(evIdB, evIdA);
#endif
  do_check_negpairs = 0;

}
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
void docellcross2(int k, double velk, int cellsk, int nc)
{
  if (velk > 0.0)
    {
      inCell[nc][k][evIdA] = inCell[nc][k][evIdA] + 1;
      cellRange[2 * k] = 1;
      if (inCell[nc][k][evIdA] == cellsk) 
	inCell[nc][k][evIdA] = 0;
    }
  else
    { 
      cellRange[2 * k + 1] = -1;
      inCell[nc][k][evIdA] = inCell[nc][k][evIdA] - 1;
      if (inCell[nc][k][evIdA] == -1) 
	inCell[nc][k][evIdA] = cellsk - 1;
    }
}
void docellcross(int k, double velk, double *rkptr, int cellsk, int nc)
{
  if (velk > 0.0)
    {
      inCell[nc][k][evIdA] = inCell[nc][k][evIdA] + 1;
      cellRange[2 * k] = 1;
#ifdef MD_GRAVITY
      if (k!=2 && inCell[nc][k][evIdA] == cellsk) 
	{
	  inCell[nc][k][evIdA] = 0;
	  *rkptr = -L2;
	  OprogStatus.DR[evIdA][k]++;
	}

#else
      if (inCell[nc][k][evIdA] == cellsk) 
	{
	  inCell[nc][k][evIdA] = 0;
	  *rkptr = -L2;
	  OprogStatus.DR[evIdA][k]++;
	}
#endif
    }
  else
    { 
      cellRange[2 * k + 1] = -1;
      inCell[nc][k][evIdA] = inCell[nc][k][evIdA] - 1;
#ifdef MD_GRAVITY
      if (k!=2 && inCell[nc][k][evIdA] == -1) 
	{
	  inCell[nc][k][evIdA] = cellsk - 1;
	  *rkptr = L2;
	  OprogStatus.DR[evIdA][k]--;
	}
#else
      if (inCell[nc][k][evIdA] == -1) 
	{
	  inCell[nc][k][evIdA] = cellsk - 1;
	  *rkptr = L2;
	  OprogStatus.DR[evIdA][k]--;
	}
#endif
    }
}
#else
void docellcross(int k, double velk, double *rkptr, int cellsk)
{
#if 0
  if (inCell[0][evIdA]+1> cellsx ||inCell[1][evIdA]+1> cellsy||inCell[2][evIdA]+1> cellsz) 
    {printf("PRIMAin cell cross ?!?\n");
      printf("velk: %f (%d,%d,%d) (%d,%d,%d) k=%d cellsk:%d\n",velk,  cellsx , cellsy,cellsz,
	     inCell[0][evIdA],inCell[1][evIdA], inCell[2][evIdA], k, cellsk );}
#endif
      if (velk > 0.0)
	{
	  inCell[k][evIdA] = inCell[k][evIdA] + 1;
	  cellRange[2 * k] = 1;
	  if (inCell[k][evIdA] == cellsk) 
	    {
	      inCell[k][evIdA] = 0;
#ifdef MD_GRAVITY
    	      if (k==2)
    		{
    		  printf("Un particella ha superato la massima altezza consentita (%.15f)\n",
    			 Lz2 + OprogStatus.extraLz);
    		  printf("Aumentare il parametro extraLz e rilanciare la simulazione\n");
    		  exit(-1);
		}
#endif
	      *rkptr = -L2;
	      OprogStatus.DR[evIdA][k]++;
	    }
	}
      else
	{ 
	  cellRange[2 * k + 1] = -1;
	  inCell[k][evIdA] = inCell[k][evIdA] - 1;
	  if (inCell[k][evIdA] == -1) 
	    {
	      inCell[k][evIdA] = cellsk - 1;
	      *rkptr = L2;
	      OprogStatus.DR[evIdA][k]--;
	    }
	}
#if 0
      if (inCell[0][evIdA]> cellsx ||inCell[1][evIdA]> cellsy||inCell[2][evIdA]> cellsz) 
	{printf("in cell cross ?!?\n");
	  printf("velk: %f(%d,%d,%d) (%d,%d,%d) k=%d cellsk:%d\n",  velk,cellsx , cellsy,cellsz,
		 inCell[0][evIdA],inCell[1][evIdA], inCell[2][evIdA], k, cellsk );}
#endif
}
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
int check_boxwall(int k, int nc, int nl)
{
  int cellsk=0;
  double vel=0.0;
#ifdef MD_GRAVITY
  if (k==2)
    return 0;
#endif
  switch (k)
    {
    case 0:
      cellsk = cellsx[nl];
      vel = vx[evIdA];
      break;
    case 1:
      cellsk = cellsy[nl];
      vel = vy[evIdA];
      break;
    case 2:
      cellsk = cellsz[nl];
      vel = vz[evIdA];
      break;
    }
  //printf("CHECK BOXWALL k=%d inCell[%d][%d][%d]:%d\n", k, nc, k, evIdA, inCell[nc][k][evIdA]);
  if ((vel < 0 && inCell[nc][k][evIdA]==0) || (vel > 0 && inCell[nc][k][evIdA]==cellsk-1))
    return 1;
  else 
    return 0;
}
extern void DeleteEvent(int id);
void ProcessCellCrossing(void)
{
  int kk, n, iA, k;
  int nc, boxwall, nlcoll, nlcross, nc_bw, nlcross_bw, nlcoll_bw;
#if 0
  if (evIdA==24)
    printf("[PROCESS CELL CROSSIN] time=%.15G\n", Oparams.time);
#endif
  UpdateAtom(evIdA);
  kk = evIdB - 100 - ATOM_LIMIT; 
  /* trattandosi di due specie qui c'è un due */
  nc = evIdC;
  //printf("time=%.15G Proc CellCrossing nc=%d kk=%d evIdA=%d\n", Oparams.time, nc, kk, evIdA);
  iA = (evIdA < Oparams.parnumA)?0:1;
  if (iA == 0 && nc == 0)
    {
      nlcoll = 0;
      nlcross = 0;
      nlcoll_bw = 3;
      nlcross_bw = 2;
      nc_bw = 1;
    }
  else if (iA == 1 && nc == 0)
    { 
      nlcoll = 1;
      nlcross = 1;
      nlcoll_bw = 2;
      nlcross_bw = 3;
      nc_bw = 1; 
    }
  else if (iA == 0 && nc == 1)
    {
      nlcoll = 3;
      nlcross = 2;
      nlcoll_bw = 0;
      nlcross_bw = 0;
      nc_bw = 0;
    }
  else /* iA == 1 && nc == 1 */
    {
      nlcoll = 2;
      nlcross = 3;
      nlcoll_bw = 1;
      nlcross_bw = 1;
      nc_bw = 0;
    }
  
  boxwall = check_boxwall(kk, nc, nlcross);
  /* questa condizione non si dovrebbe *mai* verificare, quindi 
   * in linea di principio le due linee che seguono potrebbero anche essere eliminate */
  if (nc==1 && boxwall)
    {
      printf("nc=1 and boxwall!!!! <===!!!\n");
      return;
    }
  //printf("ProcellCellCross nl=%d nc=%d k=%d\n", nl, nc, k);
  /* NOTA: cellList[i] con 0 < i < Oparams.parnum è la cella in cui si trova la particella
   * i-esima mentre cellList[j] con 
   * Oparams.parnum <= j < cellsx*cellsy*cellsz+Oparams.parnum
   * è la prima particella che si trova nella cella j-esima
   */
  n = (inCell[nc][2][evIdA] * cellsy[nlcross] + inCell[nc][1][evIdA])*cellsx[nlcross] + 
    inCell[nc][0][evIdA] + Oparams.parnum;
#if 0
  printf("nc=%d n=%d cellList[%d][%d]:%d\n",nc, n, nlcross, n, cellList[nlcross][n]);
  printf("vel=(%f,%f,%f) inCell= %d %d %d\n", vx[evIdA], vy[evIdA], vz[evIdA], inCell[nc][0][evIdA],inCell[nc][1][evIdA], inCell[nc][2][evIdA]);
#endif
  while (cellList[nlcross][n] != evIdA) 
    n = cellList[nlcross][n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellList[nlcross][n] = cellList[nlcross][evIdA];
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }

  if (boxwall)
    {
      //printf("BOXWALL nc=%d nc2=%d nl=%d nl2=%d evIdA=%d time=%.15G\n", nc, nc2, nl, nl2, evIdA, Oparams.time);
      n = (inCell[nc_bw][2][evIdA] * cellsy[nlcross_bw] + inCell[nc_bw][1][evIdA])*cellsx[nlcross_bw] + 
	inCell[nc_bw][0][evIdA]
	+ Oparams.parnum;
      while (cellList[nlcross_bw][n] != evIdA) 
	n = cellList[nlcross_bw][n];
      /* Eliminazione di evIdA dalla lista della cella n-esima della lista nl2 */
      cellList[nlcross_bw][n] = cellList[nlcross_bw][evIdA];
    }
  switch (kk)
    {
    case 0: 
      docellcross(0, vx[evIdA], &(rx[evIdA]), cellsx[nlcross], nc);
      break;
    case 1: 
      docellcross(1, vy[evIdA], &(ry[evIdA]), cellsy[nlcross], nc);
      break;
    case 2:
      docellcross(2, vz[evIdA], &(rz[evIdA]), cellsz[nlcross], nc);
      break;
    }
  PredictCellCross(evIdA, nc);
  PredictColl(evIdA, evIdB, nlcoll);
  n = (inCell[nc][2][evIdA] * cellsy[nlcross] + inCell[nc][1][evIdA])*cellsx[nlcross] + 
    inCell[nc][0][evIdA] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellList[nlcross][evIdA] = cellList[nlcross][n];
  cellList[nlcross][n] = evIdA;
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
#if 0
  printf("DOPO boxwall=%d nc=%d n=%d cellList[%d][%d]:%d\n",boxwall, nc, n, nlcross, n, cellList[nlcross][n]);
  printf("DOPO vel=(%f,%f,%f) inCell= %d %d %d\n", vx[evIdA], vy[evIdA], vz[evIdA], inCell[nc][0][evIdA],inCell[nc][1][evIdA], inCell[nc][2][evIdA]);
#endif
  if (boxwall)
    {
      switch (kk)
	{
	case 0: 
	  docellcross2(0, vx[evIdA], cellsx[nlcross_bw], nc_bw);
	  break;
	case 1: 
	  docellcross2(1, vy[evIdA], cellsy[nlcross_bw], nc_bw);
	  break;
      	case 2:
	  docellcross2(2, vz[evIdA], cellsz[nlcross_bw], nc_bw);
	  break;
	}
      if (crossevtodel[evIdA]!=-1)
	{
	  //printf("DELETING CROSS EVENT evIdA=%d\n", evIdA);
	  DeleteEvent(crossevtodel[evIdA]);
	  crossevtodel[evIdA] = -1;
	}
      PredictCellCross(evIdA, nc_bw);
      PredictColl(evIdA, evIdB, nlcoll_bw);
      n = (inCell[nc_bw][2][evIdA] * cellsy[nlcross_bw] + inCell[nc_bw][1][evIdA])*cellsx[nlcross_bw] + 
	inCell[nc_bw][0][evIdA] + Oparams.parnum;
      /* Inserimento di evIdA nella nuova cella (head) */
      cellList[nlcross_bw][evIdA] = cellList[nlcross_bw][n];
      cellList[nlcross_bw][n] = evIdA;
    }
}
#else
void ProcessCellCrossing(void)
{
  int k, n;
  UpdateAtom(evIdA);
  /* NOTA: cellList[i] con 0 < i < Oparams.parnum è la cella in cui si trova la particella
   * i-esima mentre cellList[j] con 
   * Oparams.parnum <= j < cellsx*cellsy*cellsz+Oparams.parnum
   * è la prima particella che si trova nella cella j-esima
   */
  n = (inCell[2][evIdA] * cellsy + inCell[1][evIdA] )*cellsx + inCell[0][evIdA]
    + Oparams.parnum;

  while (cellList[n] != evIdA) 
    n = cellList[n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellList[n] = cellList[evIdA];
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  k = evIdB - 100 - ATOM_LIMIT; 
  switch (k)
    {
    case 0: 
      docellcross(0, vx[evIdA], &(rx[evIdA]), cellsx);
      break;
    case 1: 
      docellcross(1, vy[evIdA], &(ry[evIdA]), cellsy);
      break;
    case 2:
      docellcross(2, vz[evIdA], &(rz[evIdA]), cellsz);
      break;
    }
  PredictEvent(evIdA, evIdB);
  n = (inCell[2][evIdA] * cellsy + inCell[1][evIdA])*cellsx + 
    inCell[0][evIdA] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellList[evIdA] = cellList[n];
  cellList[n] = evIdA;
}
#endif
void velsBrown(double T)
{
  comvel_brown(T, Oparams.m); 
}
#if 0
void rebuildLinkedList(void)
{
  int j, n;
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList[j] = -1;
  /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
  for (n = 0; n < Oparams.parnum; n++)
    {
      atomTime[n] = Oparams.time;
      inCell[0][n] =  (rx[n] + L2) * cellsx / L;
      inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
      inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
      inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
#endif
      j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + Oparams.parnum;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
void rebuildCalendar(void)
{
  int k, n, nl, nl_ignore, iA, nc, i;


  InitEventList();
  for (i=0; i < Oparams.parnum; i++)
    {
      crossevtodel[i] = -1;
    } 

  for (k = 0;  k < 3; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (n = 0; n < Oparams.parnum; n++)
    {
      iA = (n<Oparams.parnumA)?0:1;
      nl_ignore = (n<Oparams.parnumA)?1:0;
      for (nc = 0; nc < 2; nc++)
	PredictCellCross(n, nc);
      for (nl = 0; nl < 4; nl++)
	{
	  if (nl == nl_ignore || nl == iA + 2)
	    continue;
	  PredictColl(n, -2, nl); 
	}
    }
}
#else
void rebuildCalendar(void)
{
  int k, n;

  InitEventList();
  for (k = 0;  k < 3; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (n = 0; n < Oparams.parnum; n++)
    PredictEvent(n, -2); 
}
#endif
void distanza(int ia, int ib)
{
  double dx, dy, dz;
  dx = rx[ia]-rx[ib];
  dy = ry[ia]-ry[ib];
  dz = rz[ia]-rz[ib];
  dx = dx - L*rint(dx/L);
  dy = dx - L*rint(dy/L);
  //printf("dist(%d,%d): %f\n", ia, ib, sqrt(Sqr(dx)+Sqr(dy)+Sqr(dz)));
}
//void rebuildLinkedList(void);
#ifdef MD_BIG_DT
void timeshift_variables(void)
{
  int i;
  Oparams.time = 0.0;
  if (OprogStatus.rescaleTime > 0.0)
    OprogStatus.nextcheckTime -= OprogStatus.bigDt;
  OprogStatus.nextDt -= OprogStatus.bigDt;
#ifdef MD_DOUBLE_DT
  if (OprogStatus.brownian)
    OprogStatus.nextDtR -= OprogStatus.bigDt;
#endif
  if (OprogStatus.intervalSum > 0.0)
    OprogStatus.nextSumTime -= OprogStatus.bigDt;
  //nextStoreTime viene calcolato opportunamente ogni volta quindi non va shiftato
  if (OprogStatus.storerate > 0.0)
    OprogStatus.nextStoreTime -= OprogStatus.bigDt;
  for (i = 0; i < Oparams.parnum; i++)
    {
      /* NOTA 10/05/2010: se si definisce MD_BIDGT_REBUILD ad ogni passo bigDt il calendario viene
	 ricostruito completamente predicendo di nuovo tutti gli eventi (è l'extrema ratio se 
	 l'altra soluzione attualmente di default non dovesse funzionare) */
#ifdef MD_CALENDAR_HYBRID
#ifdef MD_BIGDT_REBUILD
      /* nel caso del calendario ibrido (e se si definisce MD_BIGDT_REBUILD) si fa un UpdateSystem()
	 cosicché tutto gli atomi vengono portati al tempo attuale che corrisponde 
	 ad un passo bigDt e quindi si puo' tranquillamente porre a 0 il tempo di ogni atomo */

      atomTime[i] = 0.0;
#else
      atomTime[i] -= OprogStatus.bigDt;
#endif
#else
      atomTime[i] -= OprogStatus.bigDt;
#endif
      lastcol[i] -= OprogStatus.bigDt;
#ifdef MD_ROTDIFF_MIS
      OprogStatus.lastcolltime[i] -= OprogStatus.bigDt;
#endif
    }
#ifdef MD_SURV_PROB
      sp_firstcolltime[i] -= OprogStatus.bigDt;
      sp_start_time -= OprogStatus.bigDt;
#endif
#ifdef MD_HSVISCO
  OprogStatus.lastcoll -= OprogStatus.bigDt;
#endif

}
#ifdef MD_CALENDAR_HYBRID
extern int insertInEventQ(int p);
//extern int currentIndex;
void insertInCircularLists(int idNew);
void initHQlist(void);
#endif
void timeshift_calendar(void)
{
  int poolSize, id, idB, idA;
#ifdef MD_CALENDAR_HYBRID
  int k, e, i, j;
#endif
  poolSize = Oparams.parnum*OprogStatus.eventMult;
  /* parte da 1 perché tree[0] è solo l'inzio dell'albero e non un evento */
#ifndef MD_CALENDAR_HYBRID
  for (id=1; id < poolSize; id++) 
    {
      if (treeUp[id] != -1)
	{
	  treeTime[id] -= OprogStatus.bigDt;
	  
	}
    }
#else
  for (id=1; id < poolSize; id++) 
    {
      if (treeStatus[id] != 0) /* 1 or 2 mean "node not free", i.e. used */ 
	{
	  treeTime[id] -= OprogStatus.bigDt;
	}
    }
  /* svuota l'albero binario */
  treeLeft[0] = treeRight[0] = -1;
  /*treeIdA[0] = 2*Oparams.parnum + 1;  questa non serve poiche' il nodo libero rimane invariato */
  /* azzera le linked lists che vanno ripopolate */
  for (i=0; i < OprogStatus.nlistsHQ+1; i++)
    linearLists[i] = -1;
  OprogStatus.baseIndex = 0; 
  OprogStatus.curIndex = 0;
  for (id=1; id < poolSize; id++) 
    {
      if (treeStatus[id] != 0) /* 1 or 2 mean "node not free", i.e. used */ 
	{
	  /* ripopola le linked lists e l'albero binario con i nodi preesistenti */
       	  insertInEventQ(id);
	}
    }
#endif
#if 0
      if (treeTime[id] < 0.0)
	treeTime[id] = 0.0;
#endif
}
#endif
#ifdef MD_SAVEFRA
void save_fra(void)
{
  char fileop[1024], fileop2[1024], fileop3[1024];
  double rat[5][3];
  double rcm[3];	
  int i, ret, rerr;
  FILE* f;
  sprintf(fileop2 ,"Cnf-%d-%d", 
	  OprogStatus.KK, OprogStatus.JJ);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (f = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Error saving store file!\n", NULL);
      exit(-1);
    }
  rerr = 0;
#ifdef MD_BIG_DT
  ret = fprintf(f, "%lld 0 %d %d 0\n", 
		(long long int)((Oparams.time+OprogStatus.refTime)*1000.0/Oparams.Dt), 
		Oparams.parnum, Oparams.parnum-Oparams.parnumA);
#else
  ret = fprintf(f, "%lld 0 %d %d 0\n", (long long int)(Oparams.time*1000.0/Oparams.Dt), 
		Oparams.parnum, Oparams.parnum-Oparams.parnumA) 
#endif
  rerr |= (ret < 0)?1:0;  
  ret = fprintf(f, "%.15G %.15G %.15G 0 0 %.15G\n", L, L, L, Oparams.Dt/1000.0);
  rerr |= (ret < 0)?1:0;  
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      //printf("i=%d\n",i);
      rcm[0] = rx[i]+L*OprogStatus.DR[i][0];
      rcm[1] = ry[i]+L*OprogStatus.DR[i][1];
      rcm[2] = rz[i]+L*OprogStatus.DR[i][2];
      BuildAtomPos(i, rcm, R[i], rat);
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[1][0], rat[1][1], rat[1][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[2][0], rat[2][1], rat[2][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[0][0], rat[0][1], rat[0][2]);
      rerr |= (ret < 0)?1:0;  
    }
  for (i = 0; i < Oparams.parnumA; i++)
    {
      rcm[0] = rx[i]+L*OprogStatus.DR[i][0];
      rcm[1] = ry[i]+L*OprogStatus.DR[i][1];
      rcm[2] = rz[i]+L*OprogStatus.DR[i][2];
      BuildAtomPos(i, rcm, R[i], rat);
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[1][0], rat[1][1], rat[1][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[2][0], rat[2][1], rat[2][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[0][0], rat[0][1], rat[0][2]);
      rerr |= (ret < 0)?1:0;  
    }
#ifdef MD_SAVEFRA_COMPRESSED
#ifdef MD_MAC
  sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
#else
  sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
#ifndef MD_NO_SYSTEM
  system(fileop3);
#endif
#endif
  if (rerr)
    error_on_writing(f, fileop, "save_fra", "fprintf");
  fclose(f);
}
#endif
#if defined(MD_SILICA) && defined(MD_USE_SINGLE_LL)
void rebuild_bonds(void)
{
  int i, j, k, nn;
  int cellRangeT[2 * NDIM], amin, bmin, npbonds;
  double shift[3], drx, dry, drz, dist, dists[MD_PBONDS];
  int jZ, iZ, iY, jY, iX, jX, n, na, aa, bb;
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
    }
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  for (na = 0; na < Oparams.parnum; na++)
    {
      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCell[2][na] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellsz - 1;    
	      shift[2] = - L;
	    } 
	  else if (jZ == cellsz) 
	    {
	      jZ = 0;    
	      shift[2] = L;
	    }
	  for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	    {
	      jY = inCell[1][na] + iY;    
	      shift[1] = 0.0;
	      if (jY == -1) 
		{
		  jY = cellsy - 1;    
		  shift[1] = -L;
		} 
	      else if (jY == cellsy) 
		{
		  jY = 0;    
		  shift[1] = L;
		}
	      for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
		{
		  jX = inCell[0][na] + iX;    
		  shift[0] = 0.0;
		  if (jX == -1) 
		    {
		      jX = cellsx - 1;    
		      shift[0] = - L;
		    } 
		  else if (jX == cellsx) 
		    {
		      jX = 0;   
		      shift[0] = L;
		    }
		  n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
		  for (n = cellList[n]; n > -1; n = cellList[n]) 
		    {
		      if (n <= na)
			continue;
		      /* l'interazione bonded è solo tra Si e O!! */
		      i=na; 
		      j=n;
		      drx = rx[i] - rx[j];
		      shift[0] = L*rint(drx/L);
		      dry = ry[i] - ry[j];
		      shift[1] = L*rint(dry/L);
		      drz = rz[i] - rz[j]; 
		      shift[2] = L*rint(drz/L);
		      assign_bond_mapping(i, j);
		      dist = calcDistNeg(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
		      npbonds = set_pbonds(i, j);
		      for (nn=0; nn < npbonds; nn++)
			{
			  if (dists[nn]<0.0)
			    {
			      //printf("(%d,%d)-(%d,%d)\n", i, mapbondsa[nn], j, mapbondsb[nn]);
			      aa = mapbondsa[nn];
			      bb = mapbondsb[nn];
			      add_bond(i, j, aa, bb);
			      add_bond(j, i, bb, aa);
			    }
			}
		    }
		}
	    }
	}
    }
}
#endif

void save_coordtmp_ascii(unsigned char hdWhich)
{
  char fileop[1024], fileop2[1024];
  FILE *bf;
  char fileop3[1024];
  const char sepStr[] = "@@@\n";
  int i;	
  sprintf(fileop2 ,"COORD_TMP_ASCII%d", (int)hdWhich);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "fopen: error in save_coordtmp_ascii!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
#ifdef MD_ROTDIFF_MIS
  for (i=0; i < Oparams.parnum; i++)
    {
      update_MSDrot(i);
      OprogStatus.lastcolltime[i] = Oparams.time;
    }
#endif
  R2u();
  writeAsciiPars(bf, opro_ascii);
  fprintf(bf, sepStr);
  writeAsciiPars(bf, opar_ascii);
  fprintf(bf, sepStr);
  writeAllCor(bf);
  fclose(bf);
#ifdef MPI
#ifdef MD_MAC
  sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
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
extern void angvel(void);
#ifdef MD_BIG_DT
void rebuid_all_events(void)
{
  rebuildCalendar();
  if (OprogStatus.intervalSum > 0.0)	
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
  if (OprogStatus.scalevel)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  else
    OprogStatus.scalevel = 0;
}
#endif
#ifdef MD_SAVE_REALLY_ALL
void saveTreeBondsLL(char *fn)
{
  int bf, segsize;
  int a, i, nc, nl; 
  bf = creat(fn, 0666);
  //printf("time=%.15G step:%d treeRight[0]:%d\n", Oparams.time, Oparams.curStep, treeRight[0]);
  segsize = Oparams.parnum*OprogStatus.eventMult;
  for (a = 0; a < 12; a++)
    write(bf, tree[a], sizeof(int)*segsize);
  write(bf, treeTime, sizeof(double)*segsize);
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  for (nl = 0; nl < 4; nl++)
    write(bf, cellList[nl], sizeof(int)*(cellsx[nl]*cellsy[nl]*cellsz[nl]+Oparams.parnum));
  for (nc = 0; nc < 2; nc++)
    {
      for(a = 0; a < 3; a++)
	write(bf, inCell[nc][a], sizeof(int)*Oparams.parnum); 
    }
#else
  write(bf, cellList, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  for(a = 0; a < 3; a++)
    write(bf, inCell[a], sizeof(int)*Oparams.parnum); 
#endif
  write(bf, numbonds, sizeof(int)*Oparams.parnum);
  for (i = 0; i < Oparams.parnum; i++)
    {
      write(bf, bonds[i], sizeof(int)*OprogStatus.maxbonds);
    } 
  write(bf, atomTime, sizeof(double)*Oparams.parnum);
  write(bf, lastbump, sizeof(struct LastBumpS)*Oparams.parnum);
  write(bf, lastcol,  sizeof(double)*Oparams.parnum);
  write(bf, crossevtodel, sizeof(int)*Oparams.parnum);
  close(bf);
}
void readTreeBondsLL(char *fn)
{
  int bf, segsize, a, nc, nl, i;
 
  printf("Reading binary tree: %s ", fn); 
  bf = open(fn,  EXIT | O_RDONLY, 0666);
  segsize = Oparams.parnum*OprogStatus.eventMult;
  for (a = 0; a < 12; a++)
    read(bf, tree[a], sizeof(int)*segsize);
  read(bf, treeTime, sizeof(double)*segsize);
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  for (nl = 0; nl < 4; nl++)
    read(bf, cellList[nl], sizeof(int)*(cellsx[nl]*cellsy[nl]*cellsz[nl]+Oparams.parnum));
  for (nc = 0; nc < 2; nc++)
    {
      for(a = 0; a < 3; a++)
	read(bf, inCell[nc][a], sizeof(int)*Oparams.parnum); 
    }
#else
  read(bf, cellList, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  for(a = 0; a < 3; a++)
    read(bf, inCell[a], sizeof(int)*Oparams.parnum); 
#endif
  read(bf, numbonds, sizeof(int)*Oparams.parnum);
  for (i = 0; i < Oparams.parnum; i++)
    {
      read(bf, bonds[i], sizeof(int)*OprogStatus.maxbonds);
    } 
  read(bf, atomTime, sizeof(double)*Oparams.parnum);
  read(bf, lastbump, sizeof(struct LastBumpS)*Oparams.parnum);
  read(bf, lastcol,  sizeof(double)*Oparams.parnum);
  read(bf, crossevtodel, sizeof(int)*Oparams.parnum);
  close(bf);
  //printf("time=%.15G step:%d treeRight[0]:%d\n", Oparams.time, Oparams.curStep, treeRight[0]);
  printf("...done\n");
}
void backup_all(int jj, int kk)
{
  char fn[256];
  sprintf(fn, "BinaryBak-%d-%d",jj,kk);
  saveBak(fn);
  sprintf(fn, "SaveTree-%d-%d",jj,kk);
  saveTreeBondsLL(fn);
}
void readBinBak(char *fn)
{
  int bf;
  printf("Reading binary restart file: %s ", fn);
  bf = open(fn,  EXIT | O_RDONLY, 0666);
  readBak(bf);
  close(bf);
  printf("...done\n");
}
#endif
/* ============================ >>> move<<< =================================*/
#ifdef MD_CALENDAR_HYBRID
void rebuildCalend_TS(void)
{
  int j, k, n, nl, nc, iA, nl_ignore, i;

  /* for safety reset linked lists */
  for (i=0; i < Oparams.parnum; i++)
    crossevtodel[i] = -1;
  rebuild_linked_list();
  OprogStatus.baseIndex = 0;
  OprogStatus.curIndex = 0;
  InitEventList();
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  rebuildCalendar();
  if (OprogStatus.intervalSum > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_DOUBLE_DT
  if (OprogStatus.brownian)
    ScheduleEvent(-1, ATOM_LIMIT+12,OprogStatus.nextDtR);
#endif
}
#endif
void move(void)
{
#ifdef MD_SAVE_REALLY_ALL
  int jj, kk;
#endif
  char fileop[1024], fileop2[1024];
#ifndef MD_STOREMGL
  char fileop3[1024];
#endif
  FILE *bf;
#ifndef MD_STOREMGL
  const char sepStr[] = "@@@\n";
#endif
  int innerstep=0;
  int i;
  /* Zero all components of pressure tensor */
#if 0
  Wxy = 0.0;
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;
#endif
  /* get next event */
  while (!ENDSIM)
    {
      innerstep++;
      NextEvent();
      /* Descrizione Eventi:
       * 0 <= evIdB < ATOM_LIMIT: 
       *        urto fra evIdA e evIdB 
       *
       * ATOM_LIMIT <= evIdB <= ATOM_LIMIT + 5:
       *        ATOM_LIMIT   -> urto con parete lungo x nella direzione negativa
       *        ATOM_LIMIT+1 -> urto con parete lungo x nella direzione positiva
       *
       * ATOM_LIMIT + 5 < evIdB < ATOM_LIMT + 100:
       *        eventi per usi vari (misure, output e altro)
       *
       * evIdB >= ATOM_LIMIT+100: 
       *        attraversamento della cella (cell-crossing) */        
      /* PROCESS EVENTS */ 
      if (evIdB < ATOM_LIMIT)
	{
	  MD_DEBUG(printf("collision (evIdA: %d evIdB:%d)\n", evIdA, evIdB));
	  ProcessCollision();
	  OprogStatus.collCount++;
	}
#ifdef MD_GRAVITY
      else if (evIdB >= ATOM_LIMIT + 100 || evIdB < ATOM_LIMIT + NDIM * 2)
	{
	  if (evIdB < ATOM_LIMIT + NDIM * 2)
	    ProcessCollWall();
	  else
	    ProcessCellCrossing();
	  OprogStatus.crossCount++;
	}
#else
      else if ( evIdB >= ATOM_LIMIT + 100 )
	{
	  ProcessCellCrossing();
	  OprogStatus.crossCount++;
	}
#endif
      /* ATOM_LIMIT +6 <= evIdB < ATOM_LIMIT+100 eventi che si possono usare 
       * liberamente */
      else if (evIdB == ATOM_LIMIT + 7)
	{
	  UpdateSystem();
	  OprogStatus.nextSumTime += OprogStatus.intervalSum;
	  ScheduleEvent(-1, ATOM_LIMIT + 7, OprogStatus.nextSumTime);
	  /*calcObserv();*/
	  outputSummary(); 
	}
      else if (evIdB == ATOM_LIMIT + 8)
	{
#ifdef MD_SAVE_REALLY_ALL
	  jj=OprogStatus.JJ;
	  kk=OprogStatus.KK;
#endif
#ifndef MD_SAVEFRA
	  sprintf(fileop2 ,"Store-%d-%d", 
		  OprogStatus.KK, OprogStatus.JJ);
	  /* store conf */
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	  if ( (bf = fopenMPI(fileop, "w")) == NULL)
	    {
	      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	      exit(-1);
	    }
#endif
	  UpdateSystem();
#ifdef MD_ROTDIFF_MIS
	  for (i=0; i < Oparams.parnum; i++)
	    {
	      update_MSDrot(i);
	      OprogStatus.lastcolltime[i] = Oparams.time;
	    }
#endif
	  R2u();
#ifdef MD_SAVEFRA
	  save_fra();
#else
#ifndef MD_STOREMGL
	  writeAsciiPars(bf, opro_ascii);
	  fprintf(bf, sepStr);
	  writeAsciiPars(bf, opar_ascii);
	  fprintf(bf, sepStr);
#endif
	  MD_DEBUG(printf("[Store event]: %.15G JJ=%d KK=%d\n", Oparams.time, OprogStatus.JJ, OprogStatus.KK));
#ifdef MD_STOREMGL
	  fprintf(bf, ".Vol: %f\n", L*L*L);
#endif
	  writeAllCor(bf);
	  fclose(bf);
#ifndef MD_STOREMGL
#ifdef MPI
#ifdef MD_MAC
	  sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
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
#endif
#endif
	  OprogStatus.JJ++;
	  if (OprogStatus.JJ == OprogStatus.NN)
	    {
	      OprogStatus.JJ = 0;
	      OprogStatus.KK++;
	    }
	  OprogStatus.nextStoreTime = OprogStatus.storerate *
	    (pow(OprogStatus.base,OprogStatus.NN)*OprogStatus.KK+pow(OprogStatus.base,OprogStatus.JJ));
#ifdef MD_BIG_DT
	  OprogStatus.nextStoreTime -= OprogStatus.refTime;
#endif
	  ScheduleEvent(-1, ATOM_LIMIT + 8, OprogStatus.nextStoreTime);
#ifdef MD_SAVE_REALLY_ALL
	  if (OprogStatus.saveReallyAll)
	    backup_all(jj, kk);
#endif
	}
      else if (evIdB == ATOM_LIMIT + 10)
	{
	  UpdateSystem();
#if defined(MD_ROTDIFF_MIS) && !defined(MD_DOUBLE_DT)
	  for (i=0; i < Oparams.parnum; i++)
	    {
	      update_MSDrot(i);
	      OprogStatus.lastcolltime[i] = Oparams.time;
	    }
#endif
	  R2u();
#if 0
	    {
	      double d, shift[3], dists[MD_PBONDS];
	      int i, amin, bmin;
	      shift[0] = L*rint((rx[0]-rx[1])/L);
	      shift[1] = L*rint((ry[0]-ry[1])/L);
	      shift[2] = L*rint((rz[0]-rz[1])/L);
	      d = calcDistNeg(Oparams.time, 0, 1, shift, &amin, &bmin, dists);
	      for (i=0; i < MD_PBONDS; i++)
		{
		  printf("t=%.15G dists[%d]: %.15G\n",Oparams.time,i,dists[i]);
		}
	    }
#endif
#if 0
	    {
	      static double shift[3] = {0,0,0}, vecg[8], vecgNeg[8];
	      double d,r1[3], r2[3], alpha;
	      FILE* f;
	      static int first = 1;
	      shift[0] = L*rint((rx[0]-rx[1])/L);
	      shift[1] = L*rint((ry[0]-ry[1])/L);
	      shift[2] = L*rint((rz[0]-rz[1])/L);
	      MD_DEBUG(printf("[EVENT10] shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]));
#if 1
	      d=calcDist(Oparams.time, 0, 1, shift, r1, r2, &alpha, vecg, 1);
	      if (first)
		f = fopen("distPos.dat","w");
	      else
		f = fopen("distPos.dat","a");
	      fprintf(f,"%.15G %.15G %.15G %.15G %.15G %.15G\n", Oparams.time, d,vecg[0],vecg[1],vecg[2],vecg[4]);
	      fclose(f);
#endif
	      d=calcDistNeg(Oparams.time, 0, 1, shift, r1, r2, &alpha, vecgNeg, 1);
	      if (first)
		f = fopen("distNeg.dat","w");
	      else
		f = fopen("distNeg.dat","a");
	      fprintf(f,"%.15G %.15G %.15G %.15G %.15G %.15G\n", Oparams.time, d,vecgNeg[0],vecgNeg[1],vecgNeg[2],vecgNeg[4]);
	      fclose(f);
	      if (first == 1)
		first = 0;
	    }
#endif
#if 0
	  if (sp_equilib && Oparams.ntypes==2 && typeNP[0]==1 && fabs(Oparams.time - sp_start_time) >= OprogStatus.spdeltat)
	    {
	      sp_equilib=0;
	      sp_start_time = Oparams.time;
	      /* TRAPPING */
	      typesArr[1].brownian = 0;
#if 1
	      for (i=0; i < Oparams.parnum; i++)
		{
		  if (typeOfPart[i]==1)
		    {
		      vx[i] = vy[i] = vz[i] = 0;
		      //typesArr[1].m = 1E200;
		    }
		}
#endif
	    }
#endif
	  if (OprogStatus.brownian)
	    {
#ifdef MD_HSVISCO
	      double taus, Vol;
	      if (OprogStatus.lastcoll!=-1)
		{
		  /* notare che nel caso di dinamica browniana
		   * lastcoll è in generale l'ultima collisione o tra due particelle
		   * o tra le particelle e il fluido (reset delle velocità)*/
		  taus = Oparams.time - OprogStatus.lastcoll; 
  		  OprogStatus.DQTxy += taus * OprogStatus.Txy;
		  OprogStatus.DQTyz += taus * OprogStatus.Tyz;
		  OprogStatus.DQTzx += taus * OprogStatus.Tzx;
		  OprogStatus.DQTxx += taus * OprogStatus.Txx;
		  OprogStatus.DQTyy += taus * OprogStatus.Tyy;
		  OprogStatus.DQTzz += taus * OprogStatus.Tzz;

		  DQxyOld = OprogStatus.DQxy;
		  DQyzOld = OprogStatus.DQyz;
		  DQzxOld = OprogStatus.DQzx;
		  DQxxOld = OprogStatus.DQxx;
		  DQyyOld = OprogStatus.DQyy;
		  DQzzOld = OprogStatus.DQzz;
		  OprogStatus.DQxy = OprogStatus.DQTxy + OprogStatus.DQWxy;
		  OprogStatus.DQyz = OprogStatus.DQTyz + OprogStatus.DQWyz;
		  OprogStatus.DQzx = OprogStatus.DQTzx + OprogStatus.DQWzx;
		  OprogStatus.DQxx = OprogStatus.DQTxx + OprogStatus.DQWxx;
		  OprogStatus.DQyy = OprogStatus.DQTyy + OprogStatus.DQWyy;
		  OprogStatus.DQzz = OprogStatus.DQTzz + OprogStatus.DQWzz;
		  Vol = L*L*L;
		  OprogStatus.DQxy /= Vol;
		  OprogStatus.DQyz /= Vol;
		  OprogStatus.DQzx /= Vol;
		  OprogStatus.DQxx /= Vol;
		  OprogStatus.DQyy /= Vol;
		  OprogStatus.DQzz /= Vol;

		  Pxy = (OprogStatus.DQxy - DQxyOld)/Oparams.Dt;
		  Pyz = (OprogStatus.DQyz - DQyzOld)/Oparams.Dt;
		  Pzx = (OprogStatus.DQzx - DQzxOld)/Oparams.Dt;
		  Pxx = (OprogStatus.DQxx - DQxxOld)/Oparams.Dt;
		  Pyy = (OprogStatus.DQyy - DQyyOld)/Oparams.Dt;
		  Pzz = (OprogStatus.DQzz - DQzzOld)/Oparams.Dt;
		  PxxKin = (OprogStatus.DQTxx - DQxxOldKin)/Oparams.Dt;
		  PyyKin = (OprogStatus.DQTyy - DQyyOldKin)/Oparams.Dt;
		  PzzKin = (OprogStatus.DQTzz - DQzzOldKin)/Oparams.Dt;
		  PxxHS = (OprogStatus.DQWxxHS - DQxxOldHS)/Oparams.Dt;
		  PyyHS = (OprogStatus.DQWyyHS - DQyyOldHS)/Oparams.Dt;
		  PzzHS = (OprogStatus.DQWzzHS - DQzzOldHS)/Oparams.Dt;
		  PxxST = (OprogStatus.DQWxxST - DQxxOldST)/Oparams.Dt;
		  PyyST = (OprogStatus.DQWyyST - DQyyOldST)/Oparams.Dt;
		  PzzST = (OprogStatus.DQWzzST - DQzzOldST)/Oparams.Dt;
		  press = (Pxx+Pyy+Pzz)/3.0;
		  pressKin = (PxxKin+PyyKin+PzzKin)/3.0/Vol;
		  pressHS  = (PxxHS+PyyHS+PzzHS)/3.0/Vol;
		  pressST  = (PxxST+PyyST+PzzST)/3.0/Vol; 
		  OprogStatus.lastcoll = Oparams.time;
		}
#endif
	      velsBrown(Oparams.T);
#ifdef MD_HSVISCO
	      calcT();
#endif
#if defined(MD_SILICA) && defined(MD_USE_SINGLE_LL)
	      //rebuild_linked_list();
	      //rebuild_bonds();
#endif
	      rebuildCalendar();
#ifdef MD_BIG_DT
	      if (OprogStatus.bigDt > 0.0)
		ScheduleEvent(-1,ATOM_LIMIT + 11, OprogStatus.bigDt);
#endif
	      if (OprogStatus.intervalSum > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      if (OprogStatus.scalevel)
		ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
#ifdef MD_DOUBLE_DT
	      ScheduleEvent(-1, ATOM_LIMIT+12, OprogStatus.nextDtR);
#endif
	    }
#ifdef MD_HSVISCO
	  else
	    {
	      double Vol;
	      DQxyOld = OprogStatus.DQxy;
	      DQyzOld = OprogStatus.DQyz;
	      DQzxOld = OprogStatus.DQzx;
	      DQxxOld = OprogStatus.DQxx;
	      DQyyOld = OprogStatus.DQyy;
	      DQzzOld = OprogStatus.DQzz;
	      
	      OprogStatus.DQxy = OprogStatus.DQTxy + OprogStatus.DQWxy;
	      OprogStatus.DQyz = OprogStatus.DQTyz + OprogStatus.DQWyz;
	      OprogStatus.DQzx = OprogStatus.DQTzx + OprogStatus.DQWzx;
	      OprogStatus.DQxx = OprogStatus.DQTxx + OprogStatus.DQWxx;
	      OprogStatus.DQyy = OprogStatus.DQTyy + OprogStatus.DQWyy;
	      OprogStatus.DQzz = OprogStatus.DQTzz + OprogStatus.DQWzz;
	      //printf("DQTxx: %.15G DQWxx:%.15G\n", OprogStatus.DQTxx, OprogStatus.DQWxx);
	      Vol = L*L*L;
	      OprogStatus.DQxy /= Vol;
	      OprogStatus.DQyz /= Vol;
	      OprogStatus.DQzx /= Vol;
	      OprogStatus.DQxx /= Vol;
	      OprogStatus.DQyy /= Vol;
	      OprogStatus.DQzz /= Vol;

    	      Pxy = (OprogStatus.DQxy - DQxyOld)/Oparams.Dt;
	      Pyz = (OprogStatus.DQyz - DQyzOld)/Oparams.Dt;
	      Pzx = (OprogStatus.DQzx - DQzxOld)/Oparams.Dt;
	      Pxx = (OprogStatus.DQxx - DQxxOld)/Oparams.Dt;
	      Pyy = (OprogStatus.DQyy - DQyyOld)/Oparams.Dt;
	      Pzz = (OprogStatus.DQzz - DQzzOld)/Oparams.Dt;
	      PxxKin = (OprogStatus.DQTxx - DQxxOldKin)/Oparams.Dt;
    	      PyyKin = (OprogStatus.DQTyy - DQyyOldKin)/Oparams.Dt;
	      PzzKin = (OprogStatus.DQTzz - DQzzOldKin)/Oparams.Dt;
	      PxxHS = (OprogStatus.DQWxxHS - DQxxOldHS)/Oparams.Dt;
	      PyyHS = (OprogStatus.DQWyyHS - DQyyOldHS)/Oparams.Dt;
	      PzzHS = (OprogStatus.DQWzzHS - DQzzOldHS)/Oparams.Dt;
	      PxxST = (OprogStatus.DQWxxST - DQxxOldST)/Oparams.Dt;
	      PyyST = (OprogStatus.DQWyyST - DQyyOldST)/Oparams.Dt;
	      PzzST = (OprogStatus.DQWzzST - DQzzOldST)/Oparams.Dt;
	      press = (Pxx+Pyy+Pzz)/3.0;
    	      pressKin = (PxxKin+PyyKin+PzzKin)/3.0/Vol;
	      pressHS  = (PxxHS+PyyHS+PzzHS)/3.0/Vol;
	      pressST  = (PxxST+PyyST+PzzST)/3.0/Vol; 
	      //printf("STEP #%d DQ= %f %f %f\n", Oparams.curStep, OprogStatus.DQxy, OprogStatus.DQyz, OprogStatus.DQzx);
	    }
	  DQxxOldHS = OprogStatus.DQWxxHS;
      	  DQyyOldHS = OprogStatus.DQWyyHS;
	  DQzzOldHS = OprogStatus.DQWzzHS;
	  DQxxOldST = OprogStatus.DQWxxST;
	  DQyyOldST = OprogStatus.DQWyyST;
	  DQzzOldST = OprogStatus.DQWzzST;
	  DQxxOldKin = OprogStatus.DQTxx;
      	  DQyyOldKin = OprogStatus.DQTyy;
	  DQzzOldKin = OprogStatus.DQTzz;
#endif
	  OprogStatus.nextDt += Oparams.Dt;
	  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	  break;
	}
#ifdef MD_BIG_DT
      else if (evIdB == ATOM_LIMIT + 11)
	{
	  //UpdateSystem();
	  //printf("BIG DT!\n");
#ifndef MD_CALENDAR_HYBRID
	  timeshift_calendar();
#else
#ifdef MD_BIGDT_REBUILD
	  UpdateSystem();
#else
	  //UpdateSystem();
	  /* N.B. se si pone atomTime[i] = 0 in timeshif_variables() bisogna scommentare questa riga */
	  timeshift_calendar();
#endif
#endif
	  timeshift_variables();
	  OprogStatus.refTime += OprogStatus.bigDt;
#if defined(MD_CALENDAR_HYBRID) && defined(MD_BIGDT_REBUILD)
	  /* reinitialize linke lists and predict all event again */
	  rebuildCalend_TS();
#endif
       	  ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
	  //printf("FINE BIG DT\n");
	}
#endif
#ifdef MD_DOUBLE_DT
      else if (evIdB == ATOM_LIMIT + 12)
	{
	  UpdateSystem();
#if defined(MD_ROTDIFF_MIS) 
	  for (i=0; i < Oparams.parnum; i++)
	    {
	      update_MSDrot(i);
	      OprogStatus.lastcolltime[i] = Oparams.time;
	    }
#endif
#ifndef MD_SPOT_OFF
	  angvel();
#endif
	  rebuildCalendar();		  
#ifdef MD_BIG_DT
	  if (OprogStatus.bigDt > 0.0)
	    ScheduleEvent(-1,ATOM_LIMIT + 11, OprogStatus.bigDt);
#endif
	  if (OprogStatus.intervalSum > 0.0)
	    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	  if (OprogStatus.storerate > 0.0)
	    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	  if (OprogStatus.scalevel)
	    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	  OprogStatus.nextDtR += Oparams.DtR;
	  ScheduleEvent(-1, ATOM_LIMIT+12,OprogStatus.nextDtR);
	  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	}	
#endif
      else if (evIdB == ATOM_LIMIT + 9)
	{
	  if (OprogStatus.scalevel)
	    {
	      UpdateSystem();
	      OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	      MD_DEBUG2(printf("[TAPTAU < 0] SCALVEL #%lld Vz: %.15f\n", 
			       (long long int)Oparams.curStep,Vz));
	      calc_energy(NULL);
#if 0
	      K = 0.0;
	      for (i = 0; i < Oparams.parnumA; i++)
		{
		  K += Oparams.m[0]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
		}
	      for (i = Oparams.parnumA; i < Oparams.parnum; i++)
		{
		  K += Oparams.m[1]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
		}
	      K *= 0.5;
#endif
	      scalevels(Oparams.T, K);
	      rebuildCalendar();
	      if (OprogStatus.intervalSum > 0.0)	
		ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	      if (OprogStatus.scalevel)
		ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	      else
		OprogStatus.scalevel = 0;
#ifdef MD_BIG_DT
	      if (OprogStatus.bigDt > 0.0)
		ScheduleEvent(-1,ATOM_LIMIT + 11, OprogStatus.bigDt);
#endif
	    }
	}
      /* termina la simulazione se si supera un certo MSD totale */
      if ( ( ( OprogStatus.tmsd2endA > 0.0 &&
	   DrSqTotA > Sqr(OprogStatus.tmsd2endA) ) || OprogStatus.tmsd2endA <= 0.0 ) &&
	   ( ( OprogStatus.tmsd2endB > 0.0  && 
	    DrSqTotB > Sqr(OprogStatus.tmsd2endB)) || OprogStatus.tmsd2endB <= 0.0 ) &&
	   (OprogStatus.tmsd2endA > 0.0 || OprogStatus.tmsd2endB > 0.0) )
	{
	  printf("[MSDcheck] steps %d time %.15G\n", Oparams.curStep, Oparams.time);
	  ENDSIM=1;
	}
      /* termina la simulazione dopo un certo numero di collisioni
	 OprogStatus.maxcoll > 0
	 */
      if (maxcoll > 0)
	{
	  if (numcoll >= maxcoll)
	    ENDSIM=1;
	}
#ifdef MD_BIG_DT
      if (OprogStatus.endtime > 0 && Oparams.time + OprogStatus.refTime > OprogStatus.endtime)
	ENDSIM = 1;
#else
      if (OprogStatus.endtime > 0 && Oparams.time > OprogStatus.endtime)
	ENDSIM = 1;
#endif
      if (ENDSIM)
	{
	  outputSummary();
	}

#if 0
      if (Oparams.curStep == Oparams.totStep)
	{
	  printf(" **** end of dynamics **** \n");
	  printf(" final colliding pair (%d,%d)\n", i, j);
	  /* ** CHECK FOR PARTICLE OVERLAPS ** */
	  UpdateSystem();
	  /*check (Oparams.sigma, &overlap, &K, &V);*/
	  if ( overlap ) 
	    {
	      printf(" particle overlap in final configuration\n");
	    }
	}
#endif
#if 0
      if ( ( (OprogStatus.CMreset > 0) &&
	     /* 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) */
	   (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	     /* 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) */
	  (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
	  resetCM(Oparams.parnum);

      /* Update the integral of the pressure tensor */
      updateDQ(tij);
      updatePE(Oparams.parnum);
#endif
    }
}
