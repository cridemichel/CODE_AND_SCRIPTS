#include<mdsimul.h>
#define SIMUL
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
void update_MSDrot(int i);
void update_MSD(int i);
#ifdef MD_ASYM_ITENS
void calc_omega(int i);
void calc_angmom(int i, double **I);
extern void upd_refsysM(int i);
#endif
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
extern double DphiSqA, DphiSqB, DrSqTotA, DrSqTotB;
double minaxA, minaxB, minaxAB;
int do_check_negpairs = 0;
#ifdef MD_ASYM_ITENS
double **Ia, **Ib, **invIa, **invIb, **Iatmp, **Ibtmp;
#else
double Ia, Ib, invIa, invIb;
#endif
#ifdef MD_PATCHY_HE
void bumpSP(int i, int j, int ata, int atb, double* W, int bt);
extern void assign_bond_mapping(int i, int j);
#endif
#ifdef MD_ASYM_ITENS
extern double *phi0, *psi0, *costheta0, *sintheta0, **REt, **RE0, *angM, ***RM, **REtA, **REtB, **Rdot;
extern double cosEulAng[2][3], sinEulAng[2][3];
#endif
#ifdef MD_PATCHY_HE
struct LastBumpS *lastbump;
extern void check_all_bonds(void);
#else
int *lastbump;
#endif
extern double *axa, *axb, *axc;
extern int *scdone;
extern double *maxax;
extern double calcDistNegNeighPlane(double t, double t1, int i, double *r1, double *r2, double *vecgsup, int calcguess, int calcgradandpoint, int *err, int nplane);
void calc_energy(char *msg);

extern double min3(double a, double b, double c);
extern double min(double a, double b);
extern double max3(double a, double b, double c);
extern double *lastupdNNL, *totDistDispl;
double rA[3], rB[3];
/* Routines for LU decomposition from Numerical Recipe online */
void ludcmpR(double **a, int* indx, double* d, int n);
void lubksbR(double **a, int* indx, double *b, int n);
void InvMatrix(double **a, double **b, int NB);
extern double invaSq[2], invbSq[2], invcSq[2];
double rxC, ryC, rzC, trefG;
extern int SolveLineq (double **a, double *x, int n); 
int calcdist_retcheck;
void comvel_brown (COORD_TYPE temp, COORD_TYPE *m);
void InitEventList (void);
#ifdef MD_HSVISCO
void calcT(void);
#endif
void writeAsciiPars(FILE* fs, struct pascii strutt[]);
extern void writeAllCor(FILE* fs, int saveAll);
extern struct nebrTabStruct *nebrTab;
double nextNNLrebuild;
extern void rebuildNNL(void);
extern void updrebuildNNL(int na);
extern void PredictEventNNL(int na, int nb);
extern void updAllNNL();
#ifdef MD_PATCHY_HE
extern int locate_contactSP(int i, int j, double shift[3], double t1, double t2, double *evtime, int *ata, int *atb, int *collCode);
extern void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent);
extern int *mapbondsa;
extern int *mapbondsb;
#endif
long long int itsF=0, timesF=0, itsS=0, timesS=0, numcoll=0, itsFNL=0, timesFNL=0, 
     timesSNL=0, itsSNL=0, numcalldist=0, numdisttryagain=0;
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern long long accngA, accngB;
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
double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
#ifdef MD_CALC_DPP
void store_last_u(int i)
{
  OprogStatus.lastu1x[i] = R[i][0][0];
  OprogStatus.lastu1y[i] = R[i][0][1];
  OprogStatus.lastu1z[i] = R[i][0][2];
  OprogStatus.lastu2x[i] = R[i][1][0];
  OprogStatus.lastu2y[i] = R[i][1][1];
  OprogStatus.lastu2z[i] = R[i][1][2];
  OprogStatus.lastu3x[i] = R[i][2][0];
  OprogStatus.lastu3y[i] = R[i][2][1];
  OprogStatus.lastu3z[i] = R[i][2][2];
}
#endif
double calc_dist_ij(int i, int j, double t)
{
  static double shift[3] = {0,0,0}, vecgNeg[8];
  double d,r1[3], r2[3], alpha;
  shift[0] = L*rint((rx[i]-rx[j])/L);
  shift[1] = L*rint((ry[i]-ry[j])/L);
  shift[2] = L*rint((rz[i]-rz[j])/L);
  d=calcDistNeg(t, 0, i, j, shift, r1, r2, &alpha, vecgNeg, 1);
  return d;
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
#ifdef MD_GRAVITY
int checkz(char *msg)
{
  int ii;
  double sig;
  for (ii=0; ii < Oparams.parnum; ii++)
    {
      if (ii < Oparams.parnumA)
	sig = Oparams.sigma[0][0];
      else
        sig = Oparams.sigma[1][1];
      if ((rz[ii]+Lz*0.5-sig/2.0)<0. &&  
	fabs(rz[ii]+Lz*0.5-sig/2.0) > 0.1)
	{
	  printf("[%s]ii=%d diff:%.15f\n", msg,
		 ii, rz[ii]+Lz*0.5-sig/2.0);
	  printf("sotto more in checkz!!!\n");
	  return 1;
	}
    }
  return 0;
}
#endif
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
#ifdef MD_GRAVITY
double Lz2;
#endif
double W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB;
#ifdef MD_HSVISCO
double  DQxxOld, DQyyOld, DQzzOld, DQxyOld, DQyzOld, DQzxOld, DQxxOldKin, 
       DQyyOldKin, DQzzOldKin, DQxxOldHS, DQyyOldHS, DQzzOldHS, DQxxOldST, DQyyOldST, DQzzOldST,
       PxxKin, PyyKin, PzzKin, PxxHS, PyyHS, PzzHS, PxxST, PyyST, PzzST;
#endif
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
double DrSq = 0.0; 
const double timbig = 1E12;
/* used by linked list routines */
#ifdef MD_GRAVITY
extern double g2, mgA, mgB;
#endif
double *lastcol;
double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
  cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
#ifdef MD_EDHEFLEX_OPTNNL
extern int *inCell_NNL[3], *cellList_NNL;
extern double *rxNNL, *ryNNL, *rzNNL;
#endif
int evIdA, evIdB, parnumB, parnumA;
#ifdef MD_PATCHY_HE
int evIdC, evIdD, evIdE;
extern double *treeRxC, *treeRyC, *treeRzC;
extern int *bondscache, *numbonds, **bonds;
#endif
void newtDist(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3]);
void newtDistNeg(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3],int);
void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[], double xb2[], 
	   int *nb);
void newt(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3]);
void rebuildCalendar(void);
void R2u(void);
void store_bump(int i, int j);
void check_shift(int i, int j, double *shift)
{
  double drx, dry, drz;
  if (cellsx <= 2)
    {
      drx = rx[i] - rx[j];
      shift[0] = L*rint(drx/L);
    }
  if (cellsy <= 2)
    {
      dry = ry[i] - ry[j];
      shift[1] = L*rint(dry/L);
    }
  if (cellsz <=2)
    {
      drz = rz[i] - rz[j]; 
      shift[2] = L*rint(drz/L);
    }
}
#ifdef EDHE_FLEX
extern int *is_a_sphere_NNL;
void check_inf_mass(int typei, int typej, int *infMass_i, int *infMass_j);
void saveFullStore(char* fname)
{
  char fileop3[1024], fileop2[512], fileop[512];
  int ii, i, stripStoreBak, globSaveAllBak;
  FILE *bf;
  const char sepStr[] = "@@@\n";

  sprintf(fileop2 ,fname);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii finale!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
  for (i=0; i < Oparams.parnum; i++)
    {
#ifdef MD_ASYM_ITENS
      calc_omega(i);
#endif
      update_MSDrot(i);
#ifdef MD_CALC_DPP
      update_MSD(i);
      store_last_u(i);
#endif
    }
  R2u();
  stripStoreBak = OprogStatus.stripStore;
  OprogStatus.stripStore = 0;
  globSaveAllBak = globSaveAll;
  globSaveAll = 1;
   if (mgl_mode==0)
    {
      writeAsciiPars(bf, opro_ascii);
      fprintf(bf, sepStr);
      writeAsciiPars(bf, opar_ascii);
      fprintf(bf, sepStr);
   }	      
  writeAllCor(bf, 1);
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
      system(fileop3);
    }
  OprogStatus.stripStore = stripStoreBak;
  globSaveAll = globSaveAllBak;
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
#ifdef MD_GRAVITY
void calcKVz(void)
{
  int i;
  double dd;
  Vz = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      Vz += vz[i];
    }
  Vz *= Oparams.m[0];
  dd = 0.0;
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += vz[i];
    }
  dd *= Oparams.m[1];
  Vz += dd;
  Vz /= Mtot;
  K = 0.0;
  for (i=0; i < Oparams.parnumA; i++)
    {
      K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
    }
  K *= Oparams.m[0]*0.5;
  dd = 0.0;
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
    }
  dd *= Oparams.m[1]*0.5;
  K += dd;
}
#endif
double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);

double get_min_dist_NNL (int na, int *jmin, double *rCmin, double *rDmin, double *shiftmin) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double distMin=1E10,dist,vecg[8], alpha, shift[3];
  /*double cells[NDIM];*/
  int kk, i;
  double r1[3], r2[3];
  int n;
 /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  MD_DEBUG32(printf("nebrTab[%d].len=%d\n", na, nebrTab[na].len));
  calcdist_retcheck = 0;
  for (i=0; i < nebrTab[na].len; i++)
    {
      n = nebrTab[na].list[i]; 
      if (!(n != na))
	continue;
      shift[0] = L*rint((rx[na]-rx[n])/L);
      shift[1] = L*rint((ry[na]-ry[n])/L);
      shift[2] = L*rint((rz[na]-rz[n])/L);

      if (n!=na) 
	{
	  dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
	  if (calcdist_retcheck)
	    continue;
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
  return distMin;
}

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
		  if (n!=na) 
		    {
		      dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
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

double calc_phi(void)
{
  double N = 0;
  //const double pi = acos(0)*2;
  int i ;
#ifdef EDHE_FLEX
  int typei;
  for (i=0; i < Oparams.parnum; i++)
    {
      typei = typeOfPart[i];
      N += typesArr[typei].sax[0]*typesArr[typei].sax[1]*typesArr[typei].sax[2];
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    {
      N += axa[i]*axb[i]*axc[i];
    }
#endif
  N *= 4.0*pi/3.0;
  return N / (L*L*L);
}
double calc_norm(double *vec);
void calc_ellips_norms(double *rAC, double *rBD, double *norm, double *norm2)
{ 
  int a, b;
  double modn;
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

  for (a=0; a < 3; a++)
    {
      norm2[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  norm2[a] += -Xb[a][b]*rBD[b];
	}
    }
  modn = 0.0;
  for (a = 0; a < 3; a++)
    modn += Sqr(norm2[a]);
  modn = sqrt(modn);
  for (a=0; a < 3; a++)
    norm2[a] /= modn;
}

double rcutL, aL, bL, cL;
void store_values(int i)
{
  rcutL = Oparams.rcut;
  aL = axa[i];
  bL = axb[i];
  cL = axc[i];
}
void restore_values(int i)
{
  Oparams.rcut = rcutL;
  axa[i] = aL;
  axb[i] = bL;
  axc[i] = cL;
}
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
void scale_coords(double sf)
{
  int i; 
  L *= sf;
  for (i = 0; i < Oparams.parnum; i++)
    {
      rx[i] *= sf;
      ry[i] *= sf;
      rz[i] *= sf;
    }
}
double max3(double a, double b, double c);

double scale_axes(int i, double d, double rA[3], double rC[3], double rB[3], double rD[3], 
		      double shift[3], double scalfact, double *factor, int j)
{
  int kk;
#ifndef EDHE_FLEX
  int ii;
#endif
  double nnlfact;
  double C, Ccur, F, phi0, phi, fact, L2, rAC[3], rBD[3], fact1, fact2, fact3;
  double boxdiag, factNNL=1.0;
  L2 = 0.5 * L;
  phi = calc_phi();
#ifdef MD_POLYDISP
  phi0 = 0.0;
  for (ii = 0; ii < Oparams.parnum; ii++)
    phi0 += axaP[ii]*axbP[ii]*axcP[ii];
#else
  phi0 = ((double)Oparams.parnumA)*Oparams.a[0]*Oparams.b[0]*Oparams.c[0];
  phi0 +=((double)Oparams.parnum-Oparams.parnumA)*Oparams.a[1]*Oparams.b[1]*Oparams.c[1];
#endif
  phi0 *= 4.0*pi/3.0;
  phi0 /= L*L*L;
  C = cbrt(OprogStatus.targetPhi/phi0);
#ifdef MD_POLYDISP
  Ccur = axa[i] / axaP[i];
#else
  if (i < Oparams.parnumA)
    Ccur = axa[i]/Oparams.a[0]; 
  else
    Ccur = axa[i]/Oparams.a[1]; 
#endif
  F = C / Ccur;
  if (j != -1)
    {
      for (kk=0; kk < 3; kk++)
	{
	  rAC[kk] = rA[kk] - rC[kk];
	  if (fabs (rAC[kk]) > L2)
	    rAC[kk] -= SignR(L, rAC[kk]);
	}
      for (kk=0; kk < 3; kk++)
	{
	  rBD[kk] = rB[kk] - rD[kk];
	  if (fabs (rBD[kk]) > L2)
	    rBD[kk] -= SignR(L, rBD[kk]);
	}
    }
  /* 0.99 serve per evitare che si tocchino */
  if (F < 1)
    fact = F;
  else
    {
      if (OprogStatus.useNNL)
	{
	  //maxaxNNL = 2.0*max3(nebrTab[i].axa, nebrTab[i].axb,nebrTab[i].axc);
	  //maxaxStore = 2.0*max3(aL,bL,cL);
	  factNNL = min3(nebrTab[i].axa/aL, nebrTab[i].axb/bL, nebrTab[i].axc/cL); 
	}
     /*
	 if (OprogStatus.useNNL && maxax[i] / maxaxNNL > 1.0)
	 fact1 = 1.0 + scalfact*(maxaxNNL / maxaxStore - 1.0);
	 else
	 */	 
      if (j != -1)
	{
	  //printf("===> d=%.15G norm(rAC):%.15G scalfact=%.15G\n", d, calc_norm(rAC), scalfact);
	  fact1 = 1 + scalfact*(d / (calc_norm(rAC)));//+calc_norm(rBD)));
	  if (OprogStatus.useNNL)
	    {
	      fact3 =  1.0 + 0.99*(factNNL - 1.0);
	      /* nella crescita l'ellissoide non deve uscire dal suo NNL box */ 
	      if (fact3 < fact1)
		fact1 = fact3;
	    }
	}
      else
	fact1 = 1.0 + 0.99*(factNNL - 1.0);
      fact2 = F;
      if (fact2 < fact1)
	fact = fact2;
      else
	fact = fact1;

    }
  //printf("phi=%f fact1=%.8G fact2=%.8G scaling factor: %.8G\n", phi, fact1, fact2, fact);
  axa[i] *= fact;
  axb[i] *= fact;
  axc[i] *= fact;
  maxax[i] *= fact;
  *factor = fact;
  if (!OprogStatus.useNNL)
    {
      if (2.0*max_ax(i) > Oparams.rcut)
	{
	  Oparams.rcut = 2.0*max_ax(i)*1.01;
	}
    }
  else
    {
#ifdef MD_POLYDISP
      nnlfact = axa[i]/axaP[i];
#else
      nnlfact = axa[i]/Oparams.a[i<Oparams.parnumA?0:1];
#endif
      boxdiag = 2.0*sqrt(Sqr(axa[i]+OprogStatus.rNebrShell*nnlfact)+
			 Sqr(axb[i]+OprogStatus.rNebrShell*nnlfact)+
			 Sqr(axc[i]+OprogStatus.rNebrShell*nnlfact)); 
      if ( boxdiag > Oparams.rcut)
	Oparams.rcut = 1.01*boxdiag;
    }
  return calc_phi();
}
#ifdef MD_EDHEFLEX_OPTNNL
void rebuild_linked_list_NNL()
{
  /* N.B. Se si usano le NNL ottimizzate il centro di massa geometrico delle NNL non coincide con il 
     centro di massa degli ellissoidi, in tale caso quindi le linked lists degli ellissoidi vengono usate 
     per avere una sovrastima in caso di urti con la parete dura e per mantenere gli ellissoidi nella
     first box. */
  double L2, invL;
  int j, n;
  L2 = 0.5 * L;
  invL = 1.0/L;
  cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
  cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
  cellsz = L / Oparams.rcut;
#endif 
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList_NNL[j] = -1;
  /* NOTA: rcut delle LL per le NNL ###guale a quello delle LL per gli ellissoidi
     ma quest'ultimo potrebbe essere scelto ad hoc. Inoltre le LL per gli ellissoidi 
     vengono solo usate per avere una upper limit per il tempo di collisione contro la pareti dure */  
  /* rebuild event calendar */
 
  for (n = 0; n < Oparams.parnum; n++)
    {
      /* reduce to first box */
      rxNNL[n] = nebrTab[n].r[0] - L*rint(nebrTab[n].r[0]*invL);
      ryNNL[n] = nebrTab[n].r[1] - L*rint(nebrTab[n].r[1]*invL);
      rzNNL[n] = nebrTab[n].r[2] - L*rint(nebrTab[n].r[2]*invL);
      inCell_NNL[0][n] =  (rxNNL[n] + L2) * cellsx / L;
      inCell_NNL[1][n] =  (ryNNL[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
      inCell_NNL[2][n] =  (rzNNL[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
      inCell_NNL[2][n] =  (rzNNL[n] + L2)  * cellsz / L;
#endif
      j = (inCell_NNL[2][n]*cellsy + inCell_NNL[1][n])*cellsx + 
       inCell_NNL[0][n] + Oparams.parnum;
      cellList_NNL[n] = cellList_NNL[j];
      cellList_NNL[j] = n;
    }
}
#endif
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
double check_dist_min(int i, char *msg)
{
  int j;
  double distMin=1E60, rC[3], rD[3], shift[3];
  
  j = -1;
  if (OprogStatus.useNNL)
    distMin = get_min_dist_NNL(i, &j, rC, rD, shift);
  else
    distMin = get_min_dist(i, &j, rC, rD, shift);

  if (msg)
    printf("[check_dist_min] %s distMin: %.12G\n", msg, distMin);

  return distMin;
} 
double check_alldist_min(char *msg)
{
  int j, i;
  double distMin=1E60, dist;
  double rC[3], rD[3], shift[3];
  for (i=0; i < Oparams.parnum; i++)
    {
      j = -1;
      if (OprogStatus.useNNL)
	dist = get_min_dist_NNL(i, &j, rC, rD, shift);
      else
	dist = get_min_dist(i, &j, rC, rD, shift);
      if (j != -1 && dist < distMin)
	distMin = dist;
    }
  if (msg)
    printf("[dist all] %s: %.10G\n", msg, distMin);
  return distMin;
  
}
double rcutIni;
void scale_Phi(void)
{
  int i, j, imin, kk, its, done=0;
  static int first = 1;
  static double a0I, target;
  double distMinT, distMin=1E60, rAmin[3], rBmin[3], rC[3]={0,0,0}, 
	 rD[3]={0,0,0};
  double L2, shift[3], phi, scalfact, factor, axai;
  if (OprogStatus.targetPhi <= 0)
    return;
  
  if (OprogStatus.useNNL)
    rebuildNNL();
  phi=calc_phi();
  if (first)
    {
      first = 0;
#ifdef MD_POLYDISP
      a0I = axaP[0];
#else
      a0I = Oparams.a[0];
#endif
      rcutIni = Oparams.rcut;
      target = cbrt(OprogStatus.targetPhi/calc_phi());
    }
  //UpdateSystem();   
  L2 = 0.5 * L;
  /* get the minimum distance in the system */
  phi = calc_phi();
  for (kk = 0;  kk < 3; kk++)
    {
      cellRange[2*kk]   = - 1;
      cellRange[2*kk+1] =   1;
    }
  imin = -1;
  for (i = 0; i < Oparams.parnum; i++)
    {
      j = -1;
      if (scdone[i]==1)
	{
	  done++;
	  continue;
	}
      if (OprogStatus.useNNL)
	distMin = get_min_dist_NNL(i, &j, rC, rD, shift);
      else
	distMin = get_min_dist(i, &j, rC, rD, shift);
      if (calcdist_retcheck)
	continue;
      if (j == -1 && !OprogStatus.useNNL)
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
      phi = scale_axes(i, distMin, rAmin, rC, rBmin, rD, shift, scalfact, &factor, j);
      rebuild_linked_list();
      distMinT = check_dist_min(i, NULL);
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
#ifdef MD_POLYDISP
      axai = axaP[i];
#else
      if (i < Oparams.parnumA)
	axai = Oparams.a[0];
      else
	axai = Oparams.a[1];
#endif
      if (fabs(axa[i] / axai - target) < OprogStatus.axestol)
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
  printf("Scaled axes succesfully phi=%.8G\n", phi);
#if 0
  if (fabs(phi - OprogStatus.targetPhi)<1E-8)
    {
      for (i=0; i < Oparams.parnum; i++)
	printf("axes of %d (%.15f,%.15f,%.15f)\n", i, axa[i], axb[i], axc[i]);
    }
#endif
  //check_dist_min("PRIMA");
  rebuild_linked_list();
  if (OprogStatus.useNNL)
    rebuildNNL();
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
  printf("Scaled successfully %d/%d ellipsoids \n", done, Oparams.parnum);
  if (done == Oparams.parnum || fabs(phi - OprogStatus.targetPhi)<OprogStatus.phitol)
    {
      R2u();
      ENDSIM = 1;
      /* riduce gli ellissoidi alle dimensioni iniziali e comprime il volume */
      factor = a0I/axa[0];
      Oparams.rcut = rcutIni;
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
}

long long itsNRdist=0, callsdistNR=0;
void outputSummary(void)
{
#if 0
  FILE *f;
  int i;
#endif
  /* mettere qualcosa qui */

  if (callsdistNR)
    printf("Average iteration in Newton-Raphson for distance: %.6G\n", 
	   ((double)itsNRdist)/callsdistNR);
  if (timesS>0)
    printf("Average iterations in locate_contact: %.6G\n", ((double)itsS)/timesS);
  if (timesF>0)
    printf("Average iterations in search_contact_faster: %.6G\n",  ((double)itsF)/timesF);
  if (callsfrprmn>0)
    {
      printf("percentage of convergence in frprmn: %.6G%%\n", ((double) callsok)*100.0/callsfrprmn);
      printf("avg its in frprmn: %.10G\n", ((double) itsfrprmn)/callsfrprmn);
      printf("avg percentage of gradient: %.10G%%\n", 50.0*((double)(accngA+accngB))/callsfrprmn);
    }
  if (timesFNL>0)
    printf("Average iterations in search_contact_faster_neigh: %.6G\n",  ((double)itsFNL)/timesFNL);
  if (timesSNL>0)
    printf("Average iterations in locate_contact_neigh: %.6G\n", ((double)itsSNL)/timesSNL);
  if (callsprojonto>0)
    printf("Average iterations in projonto: %.8G\n", ((double) itsprojonto)/callsprojonto);
  if (numcalldist && OprogStatus.SDmethod)
    printf("Percentage of failed dist=%.6f%%\n", 100.0*(((double) numdisttryagain) / numcalldist));
  printf("Number of collisions: %lld\n", numcoll);
#ifdef MD_PATCHY_HE 
  if (OprogStatus.checkGrazing)
    check_all_bonds();
#endif
  if (!ENDSIM)
    scale_Phi();
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
#ifdef MD_GRAVITY
  f = fopenMPI(MD_HD_MIS "Vz2.dat", "a");
  fprintf(f, "%.15f %.15f\n", Oparams.time, Sqr(Vz));
  fclose(f);
#endif
}
#ifdef EDHE_FLEX
int dofTot;
int all_spots_on_symaxis(int sa, int pt)
{
  int sp;
  int axA, axB;
  axA = (sa+1) % 3;
  axB = (axA+1) % 3;
  MD_DEBUG35(sprintf("pt=%d sa=%d axA=%d axB=%d ", pt, sa, axA, axB));
  for (sp=0; sp < typesArr[pt].nspots; sp++)
    {
      /* N.B. x[2] is the z-axis! */
      if (typesArr[pt].spots[sp].x[axA]!=0.0 || typesArr[pt].spots[sp].x[axB]!=0.0) 
	{
	  MD_DEBUG35(printf("all spots are *not* on sym axis\n"));
	  return 0;
	}
    }
  MD_DEBUG35(printf("all spots *are* on sym axis\n"));
  return 1;
}
int all_spots_in_CoM(int pt)
{
  int sp;  
  for (sp=0; sp < typesArr[pt].nspots; sp++)
    {
      /* N.B. x[2] is the z-axis! */
      if (typesArr[pt].spots[sp].x[0]!=0.0 || typesArr[pt].spots[sp].x[1]!=0.0 ||  
	  typesArr[pt].spots[sp].x[2]!=0.0 ) 
	return 0;
    }
  return 1;	
}
int two_axes_are_equal(int pt)
{
  /* tale funzione restituisce l'asse di simmetria dell'ellissoide (0=x, 1=y, 2=z) */ 
  if (typesArr[pt].sax[0] == typesArr[pt].sax[1])
    return 2;
  if (typesArr[pt].sax[0] == typesArr[pt].sax[2])
    return 1;
  if (typesArr[pt].sax[1] == typesArr[pt].sax[2])
    return 0;
  return -1;
}
int get_dof_flex(int filter)
{
  int pt, dofOfType, dofTot, sa;
  dofTot = 0;
  for (pt = 0; pt < Oparams.ntypes; pt++)
    {
      if (filter != 0 && typesArr[pt].brownian != filter)
	continue;
      /* Sphere */
      if (typesArr[pt].sax[0] == typesArr[pt].sax[1] &&
	  typesArr[pt].sax[1] == typesArr[pt].sax[2])
	{
	  /* sfere con o senza sticky spots */
	  if (typesArr[pt].nspots == 0 || (typesArr[pt].nspots!=0 && all_spots_in_CoM(pt)))	  
	    dofOfType = 3;
	  else
	    dofOfType = 5;
	}
      else if (typesArr[pt].nspots == 0)
	{
	  if (two_axes_are_equal(pt)!=-1)
	    dofOfType = 5;
	  else
	    /* ellissoide senza sticky spots */
	    dofOfType = 6;
   	}
      else
	{
	  /* loop over all spots to see whether they are along z-axis or not */
	  sa = two_axes_are_equal(pt);
	  if (sa == -1)
	    dofOfType = 6;
	  else
	    {
	      if (all_spots_on_symaxis(sa, pt))
		{
		  dofOfType = 5;
		}
	      else
		{
		  dofOfType = 6;
		}
	    }
	}
      MD_DEBUG36(printf("pt=%d dofOfType=%d filter=%d brown=%d ntypes=%d\n", pt, dofOfType, filter, typesArr[pt].brownian, Oparams.ntypes));
      dofTot += dofOfType*typeNP[pt];
    }
  /* il centro di massa dell'anticorpo è fermo */
  return dofTot;
}
#endif
#ifdef MD_GRAVITY
void scalevels(double temp, double K, double Vz)
{
  int i; 
  double sf, VVx, VVy, VVz, ddx, ddy, ddz;
  sf = sqrt( ( (3.0*((double)Oparams.parnum)-3.0) * temp ) / (2.0*K) );
  VVx = VVy = VVz = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] = (vz[i] - Vz)*sf + Vz;
      VVx += vx[i];
      VVy += vy[i];
      VVz += vz[i];
      /* scala anche i tempi di collisione! */
    } 
  VVx *= Oparams.m[0];
  VVy *= Oparams.m[0];
  VVz *= Oparams.m[0];
  ddx = ddy = ddz = 0.0;
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] = (vz[i] - Vz)*sf + Vz;
      ddx += vx[i];
      ddy += vy[i];
      ddz += vz[i];
      /* scala anche i tempi di collisione! */
    } 
  ddx *= Oparams.m[1];
  ddy *= Oparams.m[1];
  ddz *= Oparams.m[1];
  VVx += ddx;
  VVy += ddy;
  VVz += ddz;
  VVx = VVx / Mtot;
  VVy = VVy / Mtot;
  VVz = VVz / Mtot;
  for (i = 0; i < Oparams.parnum; i++)
    {
      vx[i] -= VVx;
      vy[i] -= VVy;
      vz[i] = vz[i] - VVz + Vz;
    }
  MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}
#else
#ifdef EDHE_FLEX
void scalevels(double temp, double K)
{
  int i; 
  double sf;
  double dof;
#ifdef EDHE_FLEX
  double Ti;
  static int first = 1;
#endif 
#ifdef MD_FOUR_BEADS
  dof = ((double)Oparams.parnum)*6.0-6.0;
#else
  dof = get_dof_flex(2) - OprogStatus.frozenDOF;
#endif
#ifdef EDHE_FLEX
  /* quench per il folding */
  if (OprogStatus.scalevel == 2)
    {
      if (first) 
	{
	  first = 0;
	  Ti = Oparams.T;
	}
      else
	Ti = 2.0*K/dof;
      temp = Oparams.T = OprogStatus.xi*OprogStatus.Tf + (1.0-OprogStatus.xi)*Ti;
    } 
#endif
  sf = sqrt( ( dof * temp ) / (2.0*K) );
  //printf("dof=%f temp=%.15G sf=%.15G\n", dof, temp, sf );
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (typesArr[typeOfPart[i]].brownian!=2)
	continue;
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      wx[i] *= sf;
      wy[i] *= sf;
      wz[i] *= sf;
#ifdef MD_ASYM_ITENS
      /* N.B. notare che il reference system con l'asse-z parallelo ad M non 
       * cambia poiché la direzione di M non cambia! */
      angM[i] *= sf;
      Mx[i] *= sf;
      My[i] *= sf;
      Mz[i] *= sf;
#endif
      /* scala anche i tempi di collisione! */
    } 
   MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}
#else
void scalevels(double temp, double K)
{
  int i; 
  double sf;
#ifndef EDHE_FLEX
  double dof;
#endif
#ifdef EDHE_FLEX
  sf = sqrt( ( dofTot * temp ) / (2.0*K) );
#else
  dof = OprogStatus.dofA*((double)Oparams.parnumA) + 
    OprogStatus.dofB*((double) (Oparams.parnum-Oparams.parnumA));
  sf = sqrt( ( (dof - 3.0) * temp ) / (2.0*K) );
#endif
  for (i = 0; i < Oparams.parnumA; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      wx[i] *= sf;
      wy[i] *= sf;
      wz[i] *= sf;
#ifdef MD_ASYM_ITENS
      /* N.B. notare che il reference system con l'asse-z parallelo ad M non 
       * cambia poiché la direzione di M non cambia! */
      angM[i] *= sf;
      Mx[i] *= sf;
      My[i] *= sf;
      Mz[i] *= sf;
#endif
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
#ifdef MD_ASYM_ITENS
      angM[i] *= sf;
      Mx[i] *= sf;
      My[i] *= sf;
      Mz[i] *= sf;
#endif
      /* scala anche i tempi di collisione! */
    } 
   MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}
#endif
#endif
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
#if 1
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
  MD_DEBUG11(printf("f(rC)=%.15G g(rC)=%.15G\n", f, g));
  //if (fabs(f) > 1E-5||fabs(g) > 1E-5)
   // exit(-1);
}
extern double **matrix(int n, int m);
extern void free_matrix(double **M, int n);
double calcDist(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
void update_MSDrot(int i)
{
  double ti;
  ti = Oparams.time - OprogStatus.lastcolltime[i];
  /* sumox, sumoy e sumoz sono gli integrali nel tempo delle componenti della velocità
   * angolare lungo gli assi dell'ellissoide */
  OprogStatus.sumox[i] += (wx[i]*R[i][0][0]+wy[i]*R[i][0][1]+wz[i]*R[i][0][2])*ti;
  OprogStatus.sumoy[i] += (wx[i]*R[i][1][0]+wy[i]*R[i][1][1]+wz[i]*R[i][1][2])*ti;
  OprogStatus.sumoz[i] += (wx[i]*R[i][2][0]+wy[i]*R[i][2][1]+wz[i]*R[i][2][2])*ti;
}
#ifdef MD_CALC_DPP
extern double scalProd(double *A, double *B);
extern void vectProdVec(double *A, double *B, double *C);

void update_MSD(int i)
{
  int a, b;
  double ti;
  double lu[3][3], nw[3],un, normw, dsum[3][3], nvecu[3], vel[3];
  ti = Oparams.time - OprogStatus.lastcolltime[i];
  /* sumox, sumoy e sumoz sono gli integrali nel tempo delle componenti della velocità
   * angolare lungo gli assi dell'ellissoide */
  lu[0][0] = OprogStatus.lastu1x[i];
  lu[0][1] = OprogStatus.lastu1y[i];
  lu[0][2] = OprogStatus.lastu1z[i];
  lu[1][0] = OprogStatus.lastu2x[i];
  lu[1][1] = OprogStatus.lastu2y[i];
  lu[1][2] = OprogStatus.lastu2z[i];
  lu[2][0] = OprogStatus.lastu3x[i];
  lu[2][1] = OprogStatus.lastu3y[i];
  lu[2][2] = OprogStatus.lastu3z[i];
  nw[0] = wx[i];
  nw[1] = wy[i];
  nw[2] = wz[i];
  vel[0] = vx[i];
  vel[1] = vy[i];
  vel[2] = vz[i];
  normw = calc_norm(nw);
  if (normw==0.0)
    {
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  dsum[a][b] = ti*lu[a][b];
    }
  else
    {
      for (a=0; a < 3; a++)
	nw[a] /= normw; 
      for (b=0; b < 3; b++)
	{
	  un = scalProd(nw, lu[b]);
	  for (a=0; a < 3; a++)
	    dsum[b][a] = ti*un*nw[a];
	  vectProdVec(nw, lu[b], nvecu);
	  for (a=0; a < 3; a++)
	    dsum[b][a] += (sin(normw*ti)/normw)*(lu[b][a]-un*nw[a]) - ((cos(normw*ti)-1.0)/normw)*nvecu[a];
	}
    }
  OprogStatus.sumdx[i] += scalProd(dsum[0],vel);
  OprogStatus.sumdy[i] += scalProd(dsum[1],vel);
  OprogStatus.sumdz[i] += scalProd(dsum[2],vel);
}
#endif
#ifdef MD_ASYM_ITENS
extern void calc_angmom(int i, double **I);
extern void upd_refsysM(int i);
#endif
#ifdef EDHE_FLEX
void bumpHS(int i, int j, double *W)
{
  double rxij, ryij, rzij, factor;
  double delvx, delvy, delvz, invmi, invmj, denom;
  double sigSq;
  int typei, typej;
#ifdef MD_HANDLE_INFMASS
  int infMass_i=0, infMass_j=0;
#endif
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, taus, DTxx, DTyy, DTzz;
#endif

  typei = typeOfPart[i];
  typej = typeOfPart[j];

  numcoll++;
#ifdef MD_HANDLE_INFMASS
  check_inf_mass(typei, typej, &infMass_i, &infMass_j);
  MD_DEBUG32(printf("[BUMPHS] numcoll:%lld\n", numcoll));
  if (infMass_i)
    invmi = 0.0;
  else
    invmi = 1.0/typesArr[typei].m;
  if (infMass_j)
    invmj = 0.0;
  else
    invmj = 1.0/typesArr[typej].m; 
#else
  invmi = 1.0/typesArr[typei].m;
  invmj = 1.0/typesArr[typej].m; 
#endif 
  denom = invmi + invmj; 
  sigSq = Sqr(typesArr[typei].sax[0]+typesArr[typej].sax[0]); 
  rxij = rx[i] - rx[j];
  if (fabs (rxij) > L2)
    rxij = rxij - SignR(L, rxij);
  ryij = ry[i] - ry[j];
  if (fabs (ryij) > L2)
    ryij = ryij - SignR(L, ryij);
  rzij = rz[i] - rz[j];
  if (fabs (rzij) > L2)
    rzij = rzij - SignR(L, rzij);
  factor = ( rxij * ( vx[i] - vx[j] ) +
	     ryij * ( vy[i] - vy[j] ) +
	     rzij * ( vz[i] - vz[j] ) ) / sigSq;
  factor *= 2.0 / denom;
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
  update_MSDrot(i);
  update_MSDrot(j);
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
  /* prob. quanto segue non serve */
#if 0
#ifdef MD_ASYM_ITENS
  calc_angmom(i, Ia);
  upd_refsysM(i);
  calc_angmom(j, Ib);
  upd_refsysM(j);
#endif
#endif 
}
int is_sphere(int i)
{
  int type1;
  type1 = typeOfPart[i];
  if (typesArr[type1].sax[0] == typesArr[type1].sax[1] &&
      typesArr[type1].sax[1] == typesArr[type1].sax[2])
    return 1;
  else 
    return 0;
}
int are_spheres(int i, int j)
{
  int type1, type2;
  type1 = typeOfPart[i];
  type2 = typeOfPart[j];
  if (typesArr[type1].sax[0] == typesArr[type1].sax[1] &&
      typesArr[type1].sax[1] == typesArr[type1].sax[2] &&
      typesArr[type2].sax[0] == typesArr[type2].sax[1] &&
      typesArr[type2].sax[1] == typesArr[type2].sax[2])  
    return 1;
  else 
    return 0;
}
#endif
#if defined(EDHE_FLEX) && defined(MD_HANDLE_INFMASS)
void check_inf_mass_itens(int typei, int typej, int *infMass_i, int *infMass_j, int *infItens_i, int *infItens_j)
{
#if 0
  *infMass_i = *infMass_j = *infItens_i = *infItens_j = 0;
  return;
#endif
  if (typesArr[typei].I[0] > MD_INF_ITENS||
      typesArr[typei].I[1] > MD_INF_ITENS||
      typesArr[typei].I[2] > MD_INF_ITENS)
    *infItens_i = 1;
  else
    *infItens_i = 0;
  if (typesArr[typej].I[0] > MD_INF_ITENS||
      typesArr[typej].I[1] > MD_INF_ITENS||
      typesArr[typej].I[2] > MD_INF_ITENS)
    *infItens_j = 1;
  else
    *infItens_j = 0;
  if (typesArr[typei].m > MD_INF_MASS)
    *infMass_i = 1;
  else 
    *infMass_i = 0;
  if (typesArr[typej].m > MD_INF_MASS)
    *infMass_j = 1;
  else 
    *infMass_j = 0;
}
void check_inf_mass(int typei, int typej, int *infMass_i, int *infMass_j)
{
#if 0
  *infMass_i = *infMass_j = 0;
  return;
#endif
  if (typesArr[typei].m > MD_INF_MASS)
    *infMass_i = 1;
  else 
    *infMass_i = 0;
  if (typesArr[typej].m > MD_INF_MASS)
    *infMass_j = 1;
  else 
    *infMass_j = 0;
}
#endif
#ifdef EDHE_FLEX
extern void set_angmom_to_zero(int i);
#endif
void bump (int i, int j, double rCx, double rCy, double rCz, double* W)
{
  /*
   *******************************************************************
   ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
   **                                                               **
   ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
   ** THE ROUTINE ALSO COMPUTES COLLISIONAL VIRIAL W.               **
   *******************************************************************
   */
  double factor, invmi, invmj;
  double delpx, delpy, delpz, wrx, wry, wrz, rACn[3], rBCn[3];
  double rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double norm[3];
  double modn, denom;
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, taus, DTxx, DTyy, DTzz ;
  double rxij, ryij, rzij, Dr;
#endif
#ifndef MD_ASYM_ITENS
  double factorinvIa, factorinvIb;
#endif
#ifdef MD_ASYM_ITENS
  int k1,k2;
  double rnI[3];
  double Mvec[3], omega[3];
#endif
#ifdef EDHE_FLEX
  int typei, typej;
  int infMass_j=0, infMass_i=0, infItens_i=0, infItens_j=0;
#endif
  int na, a, b;
#ifdef EDHE_FLEX
  if (are_spheres(i, j))
    {
      numcoll++;
      bumpHS(i, j, W);
      return;
    } 
#endif

  MD_DEBUG36(calc_energy("dentro bump1"));
  numcoll++;
  MD_DEBUG36(printf("[BUMP] numcoll:%lld\n", numcoll));
  MD_DEBUG36(printf("i=%d j=%d [bump] t=%f contact point: %f,%f,%f \n", i, j, Oparams.time, rxC, ryC, rzC));
  rAC[0] = rx[i] - rCx;
  rAC[1] = ry[i] - rCy;
  rAC[2] = rz[i] - rCz;
#if 1
  for (a=0; a < 3; a++)
    if (fabs (rAC[a]) > L2)
      rAC[a] -= SignR(L, rAC[a]);
#endif
  rBC[0] = rx[j] - rCx;
  rBC[1] = ry[j] - rCy;
  rBC[2] = rz[j] - rCz;
#if 0
    {
      double shift[3], r1[3], r2[3], alpha, vecgd[8], r12[3];
      shift[0] = L*rint((rx[i]-rx[j])/L);
      shift[1] = L*rint((ry[i]-ry[j])/L);
      shift[2] = L*rint((rz[i]-rz[j])/L);
      printf("shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]):
      printf("[bump] distance between %d-%d: %.15f\n", i, j, calcDistNeg(Oparams.time, i, j, shift, r1, r2, &alpha, vecgd, 1));
    }
#endif 
  for (a=0; a < 3; a++)
    {
      MD_DEBUG(printf("P rBC[%d]:%.15f ", a, rBC[a]));
      if (fabs (rBC[a]) > L2)
	rBC[a] -= SignR(L, rBC[a]);
      MD_DEBUG(printf("D rBC[%d]:%.15f ", a, rBC[a]));
    }
  MD_DEBUG(printf("\n"));
  /* calcola tensore d'inerzia e le matrici delle due quadriche */
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[i];
  typej = typeOfPart[j];
#ifdef MD_HANDLE_INFMASS
  check_inf_mass_itens(typei, typej, &infMass_i, &infMass_j, &infItens_i, &infItens_j);
#endif
  invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    { 
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
  else
    { 
      invaSq[na] = 1/Sqr(axaP[i]);
      invbSq[na] = 1/Sqr(axbP[i]);
      invcSq[na] = 1/Sqr(axcP[i]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], R[i]);
#ifdef MD_ASYM_ITENS
#ifdef EDHE_FLEX
  tRDiagR(i, Ia, typesArr[typei].I[0], typesArr[typei].I[1], typesArr[typei].I[2], R[i]);
#else
  tRDiagR(i, Ia, Oparams.I[na][0], Oparams.I[na][1], Oparams.I[na][2], R[i]);
#endif
#else
  Ia = Oparams.I[na];
#endif
  na = (j < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  invaSq[na] = 1/Sqr(typesArr[typej].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typej].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typej].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[j]);
      invbSq[na] = 1/Sqr(axbP[j]);
      invcSq[na] = 1/Sqr(axcP[j]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
#endif
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], R[j]);
#ifdef MD_ASYM_ITENS
#ifdef EDHE_FLEX
  tRDiagR(j, Ib, typesArr[typej].I[0], typesArr[typej].I[1], typesArr[typej].I[2], R[j]);
#else
  tRDiagR(j, Ib, Oparams.I[na][0], Oparams.I[na][1], Oparams.I[na][2], R[j]);
#endif
#else
  Ib = Oparams.I[na];
#endif
  //MD_DEBUG(calc_energy("dentro bump2"));
  MD_DEBUG11(check_contact(evIdA, evIdB, Xa, Xb, rAC, rBC));
  MD_DEBUG36(calc_energy("dentro bump3"));
  /* calcola le matrici inverse del tensore d'inerzia */
#ifdef MD_ASYM_ITENS
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Iatmp[k1][k2] = Ia[k1][k2];
	Ibtmp[k1][k2] = Ib[k1][k2];
      } 
#if defined(EDHE_FLEX) 
  if (!infItens_i && !is_a_sphere_NNL[i])
    {
      InvMatrix(Iatmp, invIa, 3);
      Mvec[0] = Mx[i];
      Mvec[1] = My[i];
      Mvec[2] = Mz[i];
      for (k1 = 0; k1 < 3; k1++)
	{
	  omega[k1] = 0.0;
	  for (k2 = 0; k2 < 3; k2++)
	    omega[k1] += invIa[k1][k2]*Mvec[k2]; 
	}
      wx[i] = omega[0];
      wy[i] = omega[1];
      wz[i] = omega[2];
    }
  else
    {
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  invIa[k1][k2] = 0.0;
      wx[i] = 0.0;
      wy[i] = 0.0;
      wz[i] = 0.0;
    }
  if (!infItens_j && !is_a_sphere_NNL[j])
    {
      InvMatrix(Ibtmp, invIb, 3);
      Mvec[0] = Mx[j];
      Mvec[1] = My[j];
      Mvec[2] = Mz[j];
      for (k1 = 0; k1 < 3; k1++)
	{
	  omega[k1] = 0.0;
	  for (k2 = 0; k2 < 3; k2++)
	    omega[k1] += invIb[k1][k2]*Mvec[k2]; 
	}
      wx[j] = omega[0];
      wy[j] = omega[1];
      wz[j] = omega[2];
    }  
  else
    {
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  invIb[k1][k2] = 0.0;
      wx[j] = 0.0;
      wy[j] = 0.0;
      wz[j] = 0.0;  
    }
#else
  InvMatrix(Iatmp, invIa, 3);
  InvMatrix(Ibtmp, invIb, 3);
  Mvec[0] = Mx[i];
  Mvec[1] = My[i];
  Mvec[2] = Mz[i];
  for (k1 = 0; k1 < 3; k1++)
    {
      omega[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	omega[k1] += invIa[k1][k2]*Mvec[k2]; 
    }
  wx[i] = omega[0];
  wy[i] = omega[1];
  wz[i] = omega[2];
  Mvec[0] = Mx[j];
  Mvec[1] = My[j];
  Mvec[2] = Mz[j];
  for (k1 = 0; k1 < 3; k1++)
    {
      omega[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	omega[k1] += invIb[k1][k2]*Mvec[k2]; 
    }
  wx[j] = omega[0];
  wy[j] = omega[1];
  wz[j] = omega[2];
#endif
#else
#if defined(EDHE_FLEX) && defined(MD_HANDLE_INFMASS)
  if (!infItens_i)
    invIa = 1/Ia;
  else
    invIa = 0.0;
  if (!infItens_j)
    invIb = 1/Ib;
  else
    invIb = 0.0;
#else
  invIa = 1/Ia;
  invIb = 1/Ib;
  MD_DEBUG(printf("Ia=%f Ib=%f\n", Ia, Ib));
#endif
#endif
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
#if 0 
  for (a=0; a < 3; a++)
    {
      norm2[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  norm2[a] += -Xb[a][b]*rBC[b];
	}
    }
  modn = 0.0;
  for (a = 0; a < 3; a++)
    modn += Sqr(norm2[a]);
  modn = sqrt(modn);
  for (a=0; a < 3; a++)
    norm2[a] /= modn;
#endif
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
#ifdef EDHE_FLEX
  invmi = 1.0/typesArr[typei].m;
  invmj = 1.0/typesArr[typej].m; 
#else 
  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;
#endif 
#if defined(EDHE_FLEX) && defined(MD_HANDLE_INFMASS)
  if (infMass_i)
    invmi = 0.0;
  if (infMass_j)
    invmj = 0.0;
#endif
  denom = invmi + invmj; 
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  MD_DEBUG(printf("[bump] before bump vc=%.15G\n", vc));
  //printf("[bump] before bump vc=%.15G\n", vc);
  if (vc < 0)// && fabs(vc) > 1E-10)
    {
      MD_DEBUG(printf("norm = (%f,%f,%f)\n", norm[0], norm[1],norm[2]));
      MD_DEBUG(printf("vel  = (%f,%f,%f)\n", vx[i], vy[i], vz[i]));
      MD_DEBUG(printf("i=%d r = (%f,%f,%f)\n", i, rx[i], ry[i], rz[i]));
      printf("[ERROR t=%.15G] maybe second collision has been wrongly predicted %d-%d\n",Oparams.time,i,j);
      printf("relative velocity (vc=%.15G) at contact point is negative! I ignore this event...\n", vc);
      MD_DEBUG(store_bump(i,j));
      MD_DEBUG(exit(-1));
      return;
    }
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
#ifdef EDHE_FLEX
#if 0
  if (is_a_sphere_NNL[i])
    for (a = 0; a < 3; a++)
      rBCn[a] = 0.0;
  if (is_a_sphere_NNL[j])
    for (a = 0; a < 3; a++)
      rACn[a] = 0.0;
#endif
#endif
#if 0
  if (i==2 && j==259)
    printf("rBCn[]=%.15G %.15G %.15G w=%.15G %.15G %.15G rnI=%.15G %.15G %.15G\n", rBCn[0], rBCn[1], rBCn[2],
	   wx[j], wy[j], wz[j], rnI[0], rnI[1], rnI[2]);
#endif
#ifdef MD_GRAVITY
  factor =2.0*vc/denom;
  /* Dissipation */
#if 0
  /* se si vuole avere una dissipazione nell'urto questo va sistemato!!! */
  if (!((Oparams.time - lastcol[i] < OprogStatus.tc)||
  	(Oparams.time - lastcol[j] < OprogStatus.tc)))
    factor *= mredl*(1+Oparams.partDiss);
#endif
#else
  /* SQUARE WELL: modify here */
  factor =2.0*vc/denom;
#ifdef MD_INELASTIC
  if (!((Oparams.time - lastcol[i] < OprogStatus.tc)||
  	(Oparams.time - lastcol[j] < OprogStatus.tc)))
    factor *= (1+Oparams.partDiss)/2.0;
#endif
  MD_DEBUG(printf("factor=%f denom=%f\n", factor, denom));
#endif
  delpx = - factor * norm[0];
  delpy = - factor * norm[1];
  delpz = - factor * norm[2];
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
  MD_DEBUG(printf("delp=(%f,%f,%f)\n", delpx, delpy, delpz));

  update_MSDrot(i);
  update_MSDrot(j);
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
  factorinvIa = factor*invIa;
  factorinvIb = factor*invIb;
  wx[i] += factorinvIa*rACn[0];
  wx[j] -= factorinvIb*rBCn[0];
  wy[i] += factorinvIa*rACn[1];
  wy[j] -= factorinvIb*rBCn[1];
  wz[i] += factorinvIa*rACn[2];
  wz[j] -= factorinvIb*rBCn[2];
#endif
#if 0
/* se si azzerano invIa o invIb e le vel. ang. nel caso si tratti di oggetti sferici questo codice è ridondante */
#ifdef EDHE_FLEX
  if (is_a_sphere_NNL[i])
    set_angmom_to_zero(i);
  if (is_a_sphere_NNL[j])
    set_angmom_to_zero(j);
#endif
#endif
#ifdef MD_ASYM_ITENS
#if defined(EDHE_FLEX) 
  if (!infItens_i && !is_a_sphere_NNL[i])
    {
      calc_angmom(i, Ia);
      upd_refsysM(i);
    }
  if (!infItens_j && !is_a_sphere_NNL[j])
    {
      calc_angmom(j, Ib);
      upd_refsysM(j);
    }
#else
  calc_angmom(i, Ia);
  upd_refsysM(i);
  calc_angmom(j, Ib);
  upd_refsysM(j);
#endif
#endif
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, vx[i],vy[i],vz[i], j, vx[j],vy[j],vz[j]));
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  
		  i, wx[i],wy[i],wz[i], j, wx[j],wy[j],wz[j]));
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
      Dr = L*rint((rx[i]-rx[j])/L);
      rxij = (rx[i]-rx[j]) - Dr;
      Dr = L*rint((ry[i]-ry[j])/L);
      ryij = (ry[i]-ry[j]) - Dr;
      Dr = L*rint((rz[i]-rz[j])/L);
      rzij = (rz[i]-rz[j]) - Dr;
      OprogStatus.DQWxy += rxij*delpy;
      OprogStatus.DQWyz += ryij*delpz;
      OprogStatus.DQWzx += rzij*delpx;

      OprogStatus.DQWxx += rxij*delpx;
      OprogStatus.DQWyy += ryij*delpy;
      OprogStatus.DQWzz += rzij*delpz;
      OprogStatus.DQWxxHS += rxij*delpx;
      OprogStatus.DQWyyHS += ryij*delpy;
      OprogStatus.DQWzzHS += rzij*delpz;
    }
  OprogStatus.Txy += DTxy; 
  OprogStatus.Tyz += DTyz;
  OprogStatus.Tzx += DTzx;
  OprogStatus.Txx += DTxx; 
  OprogStatus.Tyy += DTyy;
  OprogStatus.Tzz += DTzz;
#endif
}
#ifdef MD_GRAVITY
void calccmz(void)
{
  int i;
  double dd;
  rcmz = 0.0;
  
  for(i = 0; i < Oparams.parnumA; i++)
    {
      rcmz += rz[i];
    }
  rcmz *= Oparams.m[0];
  dd = 0.0;
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += rz[i];
    }
  dd *= Oparams.m[1];
  rcmz += dd;
  rcmz /= Mtot;

}
void calcRho(void)
{
  int jZ, jY, jX, n, npart, jZmax;
  double rhohcp, hhcp, Lz2, sig;
  double dia;
  Lz2 = Lz*0.5;
  /* hhcp è l'altezza delle particelle con diametro piu' piccolo 
   * se fossero close-packed */
  dia = (Oparams.sigma[0][0] < Oparams.sigma[1][1])?Oparams.sigma[0][0]:Oparams.sigma[1][1];
  rhohcp =0.7405*24/(4*pi*dia*dia*dia) ; 
  hhcp = Oparams.parnum/(rhohcp*Sqr(L));

  if (OprogStatus.rhobh <= 0)
    hhcp = Oparams.parnum/(rhohcp*Sqr(L))+OprogStatus.rhobh;
  else
    hhcp = (Oparams.parnum/(rhohcp*Sqr(L)))*OprogStatus.rhobh;
  /* Se rhobh > 0 allora l'altezza per il calcolo della densità è:
   * h_closepacking * rhobh */
  jZmax = (int) ((Lz / (OprogStatus.extraLz + Lz)) * cellsz)+1;  
  jZmax = (jZmax < cellsz)?jZmax:cellsz;
  MD_DEBUG2(printf("cellsz: %d iZmax: %d hhcp: %.15f\n", cellsz, jZmax, hhcp));
  npart = 0;
  for (jZ = 0; jZ < jZmax; jZ++)
    for (jX = 0; jX < cellsx; jX++)
      for (jY = 0; jY < cellsy; jY++)
	{
	  n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	  for (n = cellList[n]; n > -1; n = cellList[n]) 
	    {
	      if (n < parnumA)
		sig = Oparams.sigma[0][0];
	      else
		sig = Oparams.sigma[1][1];
	      if (rz[n] + Lz2 + sig*0.5 < hhcp)
		npart++;
	    }
	}
  rho = ((double)npart)/(L*L*hhcp);
}

void save_rho(void)
{
  FILE *f;
  f = fopenMPI(MD_HD_MIS "rho.dat", "a");
  fprintf(f, "%d %.15f\n", OprogStatus.numquench, rho);
  fclose(f);
}
void save_rzcm(void)
{
  FILE* f;
  f = fopenMPI(MD_HD_MIS "rcmz.dat", "a");
  fprintf(f, "%d %.15f\n", OprogStatus.numquench, rcmz);
  fclose(f);
}
#endif
#ifdef MD_GRAVITY
void calcObserv(void)
{
  int i;
  double dd1, dd;
  double dd2;
  K = 0.0;
  V = 0.0;
  Vz = 0.0;
  /* Bisogna considerare le velocità rispetto al centro di massa! */
  for (i=0; i < Oparams.parnumA; i++)
    {
      Vz += vz[i]; 
    }
  Vz *= Oparams.m[0];
  dd = 0.0;
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += vz[i];
    }
  dd *= Oparams.m[1];
  Vz += dd;
  Vz /= (Oparams.parnumA*Oparams.m[0] + (Oparams.parnum-Oparams.parnumA)*Oparams.m[1]);
  calccmz();

  for (i = 0; i < Oparams.parnumA; i++)
    {
      K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
      V += rz[i];
    }
  K *= Oparams.m[0] * 0.5;
  V *= mgA;
  dd1 = 0.0;
  dd2 = 0.0;
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd1 += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
      dd2 += rz[i];
    }
  dd1 *= Oparams.m[1] * 0.5;
  dd2 *= mgB;
  K += dd1;
  V += dd2;
}
#else
void calcObserv(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  FILE *f;
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
  if (Oparams.time + OprogStatus.refTime > 0.0)
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
  if (Oparams.time>0)
    {
      f = fopenMPI(MD_HD_MIS "D.dat", "a");
#ifdef MD_BIG_DT
      fprintf(f, "%.15f %.15f\n", Oparams.time + OprogStatus.refTime,  Dtrans);
#else
      fprintf(f, "%.15f %.15f\n", Oparams.time,  Dtrans);
#endif
      fclose(f);
    }
}
#endif
extern double *treeTime;
#ifdef MD_ASYM_ITENS
void adjust_norm(double **R)
{
  int k1, k2; 
  double n[3];
  for (k1 = 0; k1 < 3; k1++)
    {
      n[k1]=0;
      for(k2 = 0; k2 < 3; k2++)
	n[k1] += Sqr(R[k1][k2]);
      n[k1] = sqrt(n[k1]);
      if (fabs((n[k1])-1.0)>1E-10)
	{
	  MD_DEBUG(printf("Adjusting norm of orientations time=%.15f\n", Oparams.time));
	  MD_DEBUG(printf("delta = %.15f\n", fabs(n[k1]-1.0)));
	  for(k2 = 0; k2 < 3; k2++)
	    R[k1][k2] /= n[k1];
	}
    }
}
#else
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
#endif
#ifdef MD_HANDLE_INFMASS
int is_infinite_Itens(int i);
#endif
#ifdef MD_ASYM_ITENS
void symtop_evolve_orient(int i, double ti, double **Ro, double **REt, double cosea[3], double sinea[3], double *phi, double *psi);
void UpdateAtom(int i)
{
  double ti, phi, psi;
  int k1, k2;


  ti = Oparams.time - atomTime[i];
  
  rx[i] += vx[i]*ti;
  ry[i] += vy[i]*ti;
#if defined(MD_GRAVITY)
  rz[i] += vz[i]*ti - g2*Sqr(ti);
  vz[i] += -Oparams.ggrav*ti;
#else
  rz[i] += vz[i]*ti;
#endif
#ifdef MD_HANDLE_INFMASS
  if (is_infinite_Itens(i))
    {
      atomTime[i] = Oparams.time;
      return;
    }
#endif
  symtop_evolve_orient(i, ti, RA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R[i][k1][k2] = RA[k1][k2];
  phi0[i] = phi;
  psi0[i] = psi;
  //adjust_norm(R[i]);
  atomTime[i] = Oparams.time;
}
#else
void UpdateAtom(int i)
{
  double ti;
  double wSq, w, sinw, cosw;
  double Omega[3][3], OmegaSq[3][3], Rtmp[3][3], M[3][3];
  int k1, k2, k3;
  ti = Oparams.time - atomTime[i];
  
  rx[i] += vx[i]*ti;
  ry[i] += vy[i]*ti;
#if defined(MD_GRAVITY)
  rz[i] += vz[i]*ti - g2*Sqr(ti);
  vz[i] += -Oparams.ggrav*ti;
#else
  rz[i] += vz[i]*ti;
#endif
  /* ...and now orientations */
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
	      Omega[k1][k2] = -Omega[k1][k2];
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
	        R[i][k1][k2] += Rtmp[k1][k3]*M[k3][k2];
	  }
      adjust_norm(R[i]);
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
  atomTime[i] = Oparams.time;
}
#endif
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

#ifdef MD_ASYM_ITENS
extern double pi;
void calc_euler_angles(int i, double **M, double *phi, double *theta, double *psi)
{
  double sintheta;
  *theta = acos(M[2][2]);
#if 0
  if (i==117)
    {
      print_matrix(M, 3);
    } 
#endif
  if (*theta == 0.0)
    {
      *phi = *psi = 0;
      return;
    }
  sintheta = sin(*theta);
  /*
   *phi = acos(-M[2][1]/sintheta);
   */
  if (M[2][1] == 0.0)
    {
      if (M[2][0]/sintheta < 0.0)
	*phi = 1.5*pi;
      else
	*phi = 0.5*pi;
    }
  else
    {
#if 1
      *phi = atan2(M[2][0],-M[2][1]);
#else
      *phi = atan(-M[2][0]/M[2][1]);
      if (M[2][1]/sintheta > 0.0)
	*phi = pi + *phi;
#endif
    }
  //*psi = acos(-M[1][2]/sintheta);
  if (M[1][2] == 0.0)
    {
      if (M[0][2]/sintheta < 0.0)
	*psi = 1.5*pi;
      else
	*psi = 0.5*pi;
    }
  else
    {
#if 1
      *psi = atan2(M[0][2],M[1][2]);
#else
      *psi = atan(M[0][2]/M[1][2]);
      if (M[1][2]/sintheta < 0.0)
	*psi = pi + *psi;
#endif
    }
  //printf("psi=%.15G theta=%.15G phi=%.15G\n", *psi, *theta, *phi);
#if 0
  printf("*psi: %.15G psi0: %.15G sintheta: %.15G M[1][2]:%.15G\n", *psi, psi0[i], sintheta, M[1][2]);
  printf("M[1][2]:%.15G sintheta:%.15G\n", M[1][2], sintheta);
  printf("M[1][2]/sintheta=%.15G", -M[1][2]/sintheta);
#endif
}
/* matrice di eulero 
 * a_(00)	=	cos(psi)cos(phi)-cos(theta)sin(phi)sin(psi)	(6)
 * a_(01)	=	cos(psi)sin(phi)+cos(theta)cos(phi)sin(psi)	(7)
 * a_(02)	=	sin(psi)sin(theta)	(8)
 * a_(10)	=	-sin(psi)cos(phi)-cos(theta)sin(phi)cos(psi)	(9)
 * a_(11)	=	-sin(psi)sin(phi)+cos(theta)cos(phi)cos(psi)	(10)
 * a_(12)	=	cos(psi)sin(theta)	(11)
 * a_(20)	=	sin(theta)sin(phi)	(12)
 * a_(21)	=	-sin(theta)cos(phi)	(13)
 * a_(22)	=	cos(theta)	(14)
 * */

void build_euler_matrix(double cosphi, double sinphi, double costheta, double sintheta,
			double cospsi, double sinpsi, double **Reul)
{
  Reul[0][0] = cospsi*cosphi-costheta*sinphi*sinpsi;
  Reul[0][1] = cospsi*sinphi+costheta*cosphi*sinpsi;
  Reul[0][2] = sinpsi*sintheta;
  Reul[1][0] = -sinpsi*cosphi-costheta*sinphi*cospsi;
  Reul[1][1] = -sinpsi*sinphi+costheta*cosphi*cospsi;
  Reul[1][2] = cospsi*sintheta;
  Reul[2][0] = sintheta*sinphi;
  Reul[2][1] = -sintheta*cosphi;
  Reul[2][2] = costheta;
}
void evolve_euler_angles_symtop(int i, double ti, double *phi, double *psi)
{
  /* N.B. per la trottola simmetrica theta è costante */
  double I1, I3, invI1;
#ifdef EDHE_FLEX
  I1 = typesArr[typeOfPart[i]].I[0];
  I3 = typesArr[typeOfPart[i]].I[2]; 
#else
  I1 = (i < Oparams.parnumA) ? Oparams.I[0][0]:Oparams.I[1][0];
  I3 = (i < Oparams.parnumA) ? Oparams.I[0][2]:Oparams.I[1][2]; 
#endif
  invI1 = 1.0/I1;
  /* see Landau - Mechanics */
  *phi = invI1 * angM[i] * ti + phi0[i];
  *psi = ti * angM[i] * costheta0[i] * (I1 - I3) / (I3*I1) + psi0[i];
#if 0
  printf("*phi=%.15G *psi=%.15G ti=%.15G angM[%d]:%.15G I1: %.15G I2:%.15G\n", *phi, *psi, ti, i, angM[i], I1, I3);
  printf("costheta0[]:%.15G psi0[]:%.15G phi0[]=%.15G Delta=%.15G\n", costheta0[i], psi0[i], phi0[i], invI1*angM[i]*ti);
#endif
}
#ifdef EDHE_FLEX
int is_infinite_Itens(int i)
{
  int typei;
  typei = typeOfPart[i];
  if (typesArr[typei].I[0] > MD_INF_ITENS ||
      typesArr[typei].I[1] > MD_INF_ITENS ||
      typesArr[typei].I[2] > MD_INF_ITENS)
    return 1;
  else
    return 0;
}
int is_infinite_mass(int i)
{
  int typei;
  typei = typeOfPart[i];
  if (typesArr[typei].m > MD_INF_MASS)
    return 1;
  else
    return 0;
}
#endif
void UpdateOrient(int i, double ti, double **Ro, double Omega[3][3]);
void symtop_evolve_orient(int i, double ti, double **Ro, double **REt, double cosea[3], double sinea[3], double *phir, double *psir)
{
  double phi, psi, cospsi, sinpsi, cosphi, sinphi;
  int k1, k2, k3;
  //wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
#if 1
#ifdef EDHE_FLEX
  int typei;
  double It;
  typei=typeOfPart[i];
  if ((It=typesArr[typei].I[0])==typesArr[typei].I[1] &&
      typesArr[typei].I[1]==typesArr[typei].I[2])
    {
      double Omega[3][3];
      UpdateOrient(i, ti, Ro, Omega);
      return;
    }
#endif
#endif
#ifdef MD_HANDLE_INFMASS
  if (is_infinite_Itens(i))
    {    
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    Ro[k1][k2] = R[i][k1][k2];
	  }
      cosea[0] = cos(phi0[i]);
      cosea[1] = costheta0[i];
      cosea[2] = cos(psi0[i]);
      sinea[0] = sin(phi0[i]);
      sinea[1] = sintheta0[i];
      sinea[2] = sin(psi0[i]);
      *psir = psi0[i];
      *phir = phi0[i];
      return;
    }
#endif
#if 1
   if (ti == 0.0 || angM[i] == 0.0)
    {
      //printf("phi0=%.15G costheta0=%.15G sintheta0=%.15G psi0=%.15G\n", phi0[i], costheta0[i], sintheta0[i], psi0[i]);
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    Ro[k1][k2] = R[i][k1][k2];
	  }
      cosea[0] = cos(phi0[i]);
      cosea[1] = costheta0[i];
      cosea[2] = cos(psi0[i]);
      sinea[0] = sin(phi0[i]);
      sinea[1] = sintheta0[i];
      sinea[2] = sin(psi0[i]);
      *psir = psi0[i];
      *phir = phi0[i];
      return;
    }
#endif
  evolve_euler_angles_symtop(i, ti, &phi, &psi);
#if 0
  if (isnan(phi) || isnan(psi))
    {
      printf("NAN\n");
      exit(-1);
    }
#endif 
  cosphi = cos(phi);
  sinphi = sin(phi);
  cospsi = cos(psi);
  sinpsi = sin(psi);
  cosea[0] = cosphi;
  cosea[1] = costheta0[i];
  cosea[2] = cospsi;
  sinea[0] = sinphi;
  sinea[1] = sintheta0[i];
  sinea[2] = sinpsi;
  //printf("costheta: %.15G sintheta: %.15G\n", costheta0[i], sintheta0[i]);
  //printf("cosphi: %.15G sinphi: %.15G cospsi: %.15G sinpsi:%.15G\n", cosphi, sinphi, cospsi, sinpsi);
  build_euler_matrix(cosphi, sinphi, costheta0[i], sintheta0[i],
		     cospsi, sinpsi, REt);
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro [k1][k2] = 0.0;
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += REt[k1][k3]*RM[i][k3][k2];
      }
  adjust_norm(Ro);
  *psir = psi;
  *phir = phi;
}
#endif
void UpdateOrient(int i, double ti, double **Ro, double Omega[3][3])
{ 
  double wSq, w, OmegaSq[3][3], M[3][3];
  double sinw, cosw;
  int k1, k2, k3;
#ifdef MD_HANDLE_INFMASS
  if (is_infinite_Itens(i))
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
      return;
    }
#endif
  wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
  w = sqrt(wSq);
  if (ti != 0.0 && w != 0.0) 
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
	      M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
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
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    Ro[k1][k2] = R[i][k1][k2];
	    for (k3 = 0; k3 < 3; k3++)
	      //Ro[k1][k2] += M[k1][k3]*R[i][k3][k2];
	        Ro[k1][k2] += R[i][k1][k3]*M[k3][k2];
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

#ifndef MD_APPROX_JACOB
extern double **matrix(int n, int m);
extern void free_matrix(double **M, int n);
#ifdef MD_ASYM_ITENS
void calc_Rdot(int i, double cosea[3], double sinea[3], double **Ro)
{
  int k1,k2,k3;
  /* cosea[] = {0:cos(phi),1:cos(theta),2:cos(psi)}
   * sinea[] = {0:sin(phi),1:sin(theta),2:sin(psi)}*/
  /* here we assume d/dt[theta(t)] = 0 */
  double A, B, I1, I3;
  double costh, sinth, cospsi, sinphi, cosphi, sinpsi;
  double sinphisinpsi, cosphicospsi, sinphicospsi, cosphisinpsi;
  double cosphisinth;
  cosphi = cosea[0];
  sinphi = sinea[0];
  costh = cosea[1];
  sinth = sinea[1];
  cospsi = cosea[2];
  sinpsi = sinea[2];
  sinphisinpsi = sinphi*sinpsi;
  cosphicospsi = cosphi*cospsi;
  sinphicospsi = sinphi*cospsi;
  cosphisinpsi = cosphi*sinpsi;
  cosphisinth = cosphi*sinth;
#ifdef EDHE_FLEX
  I1 = typesArr[typeOfPart[i]].I[0];
  I3 = typesArr[typeOfPart[i]].I[2];
#else
  I1 = Oparams.I[i<Oparams.parnumA?0:1][0];
  I3 = Oparams.I[i<Oparams.parnumA?0:1][2];
#endif
  A = angM[i]/I1;
  B = angM[i]*cosea[1]*(I1 - I3)/(I1*I3);
  REt[0][0] = -B*cosphisinpsi - A*sinphicospsi - costh*( A*cosphisinpsi + B*sinphicospsi);
  REt[0][1] = -B*sinphisinpsi + A*cosphicospsi + costh*( B*cosphicospsi - A*sinphisinpsi);
  REt[0][2] =  B*cospsi*sinth;
  REt[1][0] = -B*cosphicospsi + A*sinphisinpsi - costh*( A*cosphicospsi - B*sinphisinpsi);
  REt[1][1] = -B*sinphicospsi - A*cosphisinpsi + costh*(-A*sinphicospsi - B*cosphisinpsi);
  REt[1][2] = -B*sinpsi*sinth;
  REt[2][0] =  A*cosphi*sinth;
  REt[2][1] =  A*sinphi*sinth;
  REt[2][2] =  0;
  for (k1=0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = 0.0;
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += REt[k1][k3] * RM[i][k3][k2];	
      }
}
/* N.B. questa va riscritta per la trottola simmetrica! */
void calcFxtFt(int i, double x[3], double **RM, double cosea[3], double sinea[3], double **X,
	       double D[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft)
{
  double tRDRdot[3][3], tRdotDR[3][3], DR[3][3]; 
  double DtX[3][3], dx[3];
  int k1, k2, k3;
  calc_Rdot(i, cosea, sinea, Rdot);
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  DR[k1][k2] = D[k1][k1]*R[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRDRdot[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (DR[k3][k1] == 0.0 || Rdot[k3][k2] == 0.0)
		continue;
	      tRDRdot[k1][k2] += DR[k3][k1]*Rdot[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRdotDR[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (Rdot[k3][k1] == 0.0 || DR[k3][k2] == 0.0)
		continue;
	      tRdotDR[k1][k2] += Rdot[k3][k1]*DR[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      DtX[k1][k2] = tRdotDR[k1][k2] + tRDRdot[k1][k2];
  for (k1 = 0; k1 < 3; k1++)
    {
      Fxt[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	Fxt[k1] += DtX[k1][k2]*(x[k2]-pos[k2]) - X[k1][k2]*vel[k2]; 
      Fxt[k1] *= 2.0;
     } 
   *Ft = 0;
   for (k1 = 0; k1 < 3; k1++)
     dx[k1] = x[k1]-pos[k1];
   for (k1 = 0; k1 < 3; k1++)
     {
       for (k2 = 0; k2 < 3; k2++)
	 {
	   *Ft += -vel[k1]*X[k1][k2]*dx[k2]+dx[k1]*DtX[k1][k2]*dx[k2]-dx[k1]*X[k1][k2]*vel[k2];
	 }
     }
}
#else
void calcFxtFt(double x[3], double **X,
	       double D[3][3], double Omega[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft)
{
  double OmegaX[3][3], XOmega[3][3]; 
  double DtX[3][3], dx[3];
  int k1, k2, k3;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  XOmega[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (X[k1][k3] == 0.0 || Omega[k3][k2] == 0.0)
		continue;
	      XOmega[k1][k2] += X[k1][k3]*Omega[k3][k2];
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  OmegaX[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (X[k3][k2] == 0.0 || Omega[k3][k1] == 0.0)
		continue;
	      OmegaX[k1][k2] += Omega[k1][k3]*X[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      DtX[k1][k2] = OmegaX[k1][k2] - XOmega[k1][k2];
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  Mtmp[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    Mtmp[k1][k2] += sumDOmega[k1][k3]*R[k3][k2]; 
	}
    }
   for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  DtX[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    DtX[k1][k2] += R[k3][k1]*Mtmp[k3][k2]; 
	}
    }
#endif
   for (k1 = 0; k1 < 3; k1++)
     {
       Fxt[k1] = 0;
       for (k2 = 0; k2 < 3; k2++)
	 Fxt[k1] += DtX[k1][k2]*(x[k2]-pos[k2]) - X[k1][k2]*vel[k2]; 
       Fxt[k1] *= 2.0;
     } 
   *Ft = 0;
   for (k1 = 0; k1 < 3; k1++)
     dx[k1] = x[k1]-pos[k1];
   for (k1 = 0; k1 < 3; k1++)
     {
       for (k2 = 0; k2 < 3; k2++)
	 {
	   *Ft += -vel[k1]*X[k1][k2]*dx[k2]+dx[k1]*DtX[k1][k2]*dx[k2]-dx[k1]*X[k1][k2]*vel[k2];
	 }
     }
}
#endif
//#define MD_GLOBALNR
#undef MD_GLOBALNR2
void fdjacGuess(int n, double x[], double fvec[], double **df, 
	   void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  double fx[3], gx[3], tmp;
  int k1, k2;

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*(Xa[k1][k2] + Sqr(x[3])*Xb[k1][k2]);
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]);
	}
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = 2.0*x[3]*gx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1] - gx[k1];
      printf("[fdjacGuess] fx[%d]:%.15G gx=%.15G\n",k1 , fx[k1],gx[k1]);
    } 

  df[3][3] = 0.0;
#ifndef MD_GLOBALNR2
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = 0.0;
  tmp = 0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      tmp += (x[k1]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0 - (0.5*tmp-1.0);
#endif
}
/* funzione che calcola lo Jacobiano */
void fdjac(int n, double x[], double fvec[], double **df, 
	   void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  int na; 
#ifdef EDHE_FLEX
  int typei, typej;
#endif
  double  rA[3], rB[3], ti, vA[3], vB[3];
#ifndef MD_ASYM_ITENS
  double OmegaA[3][3], OmegaB[3][3];
#endif
  double DA[3][3], DB[3][3], fx[3], gx[3];
  double Fxt[3], Gxt[3], Ft, Gt;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  int k1, k2;
  ti = x[4] + (trefG - atomTime[iA]);
  rA[0] = rx[iA] + vx[iA]*ti;
  rA[1] = ry[iA] + vy[iA]*ti;
  rA[2] = rz[iA] + vz[iA]*ti;
  vA[0] = vx[iA];
  vA[1] = vy[iA];
  vA[2] = vz[iA];
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(iA, ti, RA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(iA, ti, RA, OmegaA);
#endif
  MD_DEBUG2(printf("i=%d ti=%f", iA, ti));
  MD_DEBUG2(print_matrix(RA, 3));
  na = (iA < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[iA];
  invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[iA]);
      invbSq[na] = 1/Sqr(axb[iA]);
      invcSq[na] = 1/Sqr(axc[iA]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[iA]);
      invbSq[na] = 1/Sqr(axbP[iA]);
      invcSq[na] = 1/Sqr(axcP[iA]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[iA]);
      invbSq[na] = 1/Sqr(axb[iA]);
      invcSq[na] = 1/Sqr(axc[iA]);
    }
#endif
  tRDiagR(iA, Xa, invaSq[na], invbSq[na], invcSq[na], RA);
  MD_DEBUG2(printf("invabc: (%f,%f,%f)\n", invaSq[na], invbSq[na], invcSq[na]));
  MD_DEBUG2(print_matrix(Xa, 3));
  DA[0][1] = DA[0][2] = DA[1][0] = DA[1][2] = DA[2][0] = DA[2][1] = 0.0;
  DA[0][0] = invaSq[na];
  DA[1][1] = invbSq[na];
  DA[2][2] = invcSq[na];
  ti = x[4] + (trefG - atomTime[iB]);
  rB[0] = rx[iB] + vx[iB]*ti + shift[0];
  rB[1] = ry[iB] + vy[iB]*ti + shift[1];
  rB[2] = rz[iB] + vz[iB]*ti + shift[2];
  vB[0] = vx[iB];
  vB[1] = vy[iB];
  vB[2] = vz[iB];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(iB, ti, RB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  UpdateOrient(iB, ti, RB, OmegaB);
#endif
  na = (iB < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typej = typeOfPart[iB];
  invaSq[na] = 1/Sqr(typesArr[typej].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typej].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typej].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[iB]);
      invbSq[na] = 1/Sqr(axb[iB]);
      invcSq[na] = 1/Sqr(axc[iB]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[iB]);
      invbSq[na] = 1/Sqr(axbP[iB]);
      invcSq[na] = 1/Sqr(axcP[iB]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[iB]);
      invbSq[na] = 1/Sqr(axb[iB]);
      invcSq[na] = 1/Sqr(axc[iB]);
    }
#endif
  tRDiagR(iB, Xb, invaSq[na], invbSq[na], invcSq[na], RB);
  DB[0][1] = DB[0][2] = DB[1][0] = DB[1][2] = DB[2][0] = DB[2][1] = 0.0;
  DB[0][0] = invaSq[na];
  DB[1][1] = invbSq[na];
  DB[2][2] = invcSq[na];
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*(Xa[k1][k2] + Sqr(x[3])*Xb[k1][k2]);
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]);
	}
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[3][k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[3][k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]); 
#endif
      df[3][k1] = fx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[4][k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[4][k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]); 
#endif
      df[4][k1] = gx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[k1][3] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[k1][3] += 4.0*x[3]*Xb[k1][k2]*(x[k2]-rB[k2]); 
#endif
      df[k1][3] = 2.0*x[3]*gx[k1];
    } 
  df[3][3] = 0.0;
  df[4][3] = 0.0;
#ifdef MD_ASYM_ITENS
  calcFxtFt(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fx, Fxt, &Ft);
  calcFxtFt(iB, x, RM[iB], cosEulAng[1], sinEulAng[1], Xb, DB, RB, rB, vB, gx, Gxt, &Gt);
#else
  calcFxtFt(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
  calcFxtFt(x, Xb, DB, OmegaB, RB, rB, vB, gx, Gxt, &Gt);
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      //df[k1][4] = 0;
      //for (k2 = 0; k2 < 3; k2++)
      df[k1][4] = Fxt[k1]+Sqr(x[3])*Gxt[k1]; 
    } 
 df[3][4] = Ft;
 df[4][4] = Gt;
#ifndef MD_GLOBALNR
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
 MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
#endif
}
#endif
void upd2tGuess(int i, int j, double shift[3], double tGuess)
{
  double ti;
  int na;
#ifdef EDHE_FLEX
  int typei, typej;
#endif
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
#ifdef MD_ASYM_ITENS
  double cosea[3], sinea[3], phi, psi;
#endif
  ti = tGuess - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, Rt, REt, cosea, sinea, &phi, &psi);
#else
  UpdateOrient(i, ti, Rt, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[i];
  invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[i]);
      invbSq[na] = 1/Sqr(axbP[i]);
      invcSq[na] = 1/Sqr(axcP[i]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = tGuess - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, Rt, REt, cosea, sinea, &phi, &psi);
#else
  UpdateOrient(j, ti, Rt, Omega);
#endif
  na = (j < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typej = typeOfPart[j];
  invaSq[na] = 1/Sqr(typesArr[typej].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typej].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typej].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[j]);
      invbSq[na] = 1/Sqr(axbP[j]);
      invcSq[na] = 1/Sqr(axcP[j]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
#endif
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);

}
void funcs2beZeroedGuess(int n, double x[], double fvecG[], int i, int j, double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3], tmp;
  /* x = (r, alpha, t) */ 
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvecG[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
#if 0
  fvec[3] = -1.0;
  fvec[4] = -1.0;
#endif
#if 0
  MD_DEBUG(printf("fx+Sqr(alpha)*gx=(%f,%f,%f) fx=(%f,%f,%f) gx=(%f,%f,%f)\n", fvec[0], fvec[1], fvec[2],
	 fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
#endif
  fvecG[3] = 0.0;
  tmp = 0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvecG[3] += (x[k1]-rA[k1])*fx[k1];
      tmp += (x[k1]-rB[k1])*gx[k1];
    }
  fvecG[3] = 0.5*fvecG[3]-1.0 - (0.5*tmp-1.0);
}

void funcs2beZeroed(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int na, k1, k2; 
  double  rA[3], rB[3], ti;
  double fx[3], gx[3];
#ifdef EDHE_FLEX
  int typei, typej;
#endif
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
#ifdef MD_ASYM_ITENS
  double cosea[3], sinea[3], phi, psi;
#endif
  /* x = (r, alpha, t) */ 
  
  ti = x[4] + (trefG - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, Rt, REt, cosea, sinea, &phi, &psi);
#else
  UpdateOrient(i, ti, Rt, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[i];
  invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[i]);
      invbSq[na] = 1/Sqr(axbP[i]);
      invcSq[na] = 1/Sqr(axcP[i]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = x[4] + (trefG - atomTime[j]);
  MD_DEBUG(printf("x[4]:%.15f atomTime[%d]:%.15f\n",x[4], j, atomTime[j]));
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(j, ti, Rt, REt, cosea, sinea, &phi, &psi);
#else
  UpdateOrient(j, ti, Rt, Omega);
#endif
  na = (j < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typej = typeOfPart[j];
  invaSq[na] = 1/Sqr(typesArr[typej].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typej].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typej].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[j]);
      invbSq[na] = 1/Sqr(axbP[j]);
      invcSq[na] = 1/Sqr(axcP[j]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
#endif
 tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
      MD_DEBUG(printf("[FUNC2BEZ]x[%d]:%.15f rA[%d]:%f fx:%.15f\n", k1, x[k1], k1, rA[k1],fx[k1]));

    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2] - rB[k2]);
      MD_DEBUG(printf("[FUNC2BEZ]x[%d]:%.15f rB[%d]:%f gx:%.15f\n", k1, x[k1], k1, rB[k1],gx[k1]));
    }

  MD_DEBUG(print_matrix(Xb,3));
  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      fvec[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fvec[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]) + 2.0*Sqr(x[3])*Xb[k1][k2]*(x[k2] - rB[k2]);
#endif
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
#if 0
  fvec[3] = -1.0;
  fvec[4] = -1.0;
#endif
#if 0
  MD_DEBUG(printf("fx+Sqr(alpha)*gx=(%f,%f,%f) fx=(%f,%f,%f) gx=(%f,%f,%f)\n", fvec[0], fvec[1], fvec[2],
	 fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
#endif
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      for (k2 = 0; k2 < 3; k2++)
	{
	  fvec[3] += (x[k1]-rA[k1])*Xa[k1][k2]*(x[k2]-rA[k2]);
	  fvec[4] += (x[k1]-rB[k1])*Xb[k1][k2]*(x[k2]-rB[k2]);
	}
#endif
#if 1
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1]-rB[k1])*gx[k1];
#endif
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;
  MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
}
extern double max3(double a, double b, double c);
double max(double a, double b);

double scalProd(double *A, double *B);

double rA[3], rB[3];
//#define MD_GLOBALNRD
int fdjac_disterr;
void fdjacDistNeg5SE(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int, int, double []), 
		   int iA, int iB, double shift[3], double *fx, double *gx);

void fdjacDistNeg5(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3], double *fx, double *gx)
{
  double rDC[3], rD[3], A[3][3], b[3], c[3];
  int k1, k2;
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(iA) || is_superellipse(iB))
    {
      fdjacDistNeg5SE(n, x, fvec, df, vecfunc, iA, iB, shift, fx, gx);
      return;
    }
#endif
 
#ifdef MD_USE_CBLAS
  double XaL[3][3], XbL[3][3];
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	XaL[k1][k2] = Xa[k1][k2];
	XbL[k1][k2] = Xb[k1][k2];
	A[k1][k2] = 2.0*Xb[k1][k2]*Sqr(x[3]);
      }
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      3, 3, 3, 4.0*Sqr(x[3])*x[4], &XaL[0][0],
	      3, &XbL[0][0], 3,
	      1.0, &A[0][0], 3);
#else
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
#if 0
	  A[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    A[k1][k2] += Xb[k1][k3]*Xa[k3][k2];
#else
	  A[k1][k2] = XbXa[k1][k2];
#endif
	  A[k1][k2] *= 4.0*Sqr(x[3])*x[4];
	  A[k1][k2] += 2.0*Xb[k1][k2]*Sqr(x[3]);
	}
    }	
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2] + A[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	}
      rD[k1] = x[k1] + fx[k1]*x[4];
    } 
  //printf("rC: %f %f %f rD: %f %f %f\n", x[0], x[1], x[2], rD[0], rD[1], rD[2]);
  //printf("fx: %f %f %f x[4]: %f\n", fx[0], fx[1], fx[2], x[4]);
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  gx[k1] += 2.0*Xb[k1][k2]*(rD[k2]-rB[k2]);
	}
    }
#if 1
  if (OprogStatus.SDmethod==2 || OprogStatus.SDmethod==3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = rD[k1] - x[k1];
#ifdef EDHE_FLEX
      for (kk=0; kk < 3; kk++)
	{
	  axi[kk] = typesArr[typeOfPart[iA]].sax[kk];
	  axj[kk] = typesArr[typeOfPart[iB]].sax[kk];
	}
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axi[0],axi[1],axi[2])+max3(axj[0],axj[1],axj[2])))
	{
	  fdjac_disterr = 1;	
	}

#else
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axa[iA],axb[iA],axc[iA])+max3(axa[iB],axb[iB],axc[iB])))
	{
	  fdjac_disterr = 1;	
	}
#endif
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      b[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  b[k1] += Xb[k1][k2]*fx[k2];
	}
      b[k1] *= 2.0*Sqr(x[3]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      c[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  c[k1] += gx[k2]*Xa[k2][k1];
	}
      c[k1] *= 2.0*x[4];
      c[k1] += gx[k1];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  df[3][3] = 0.0;
  df[3][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = c[k1];
    } 
  df[4][3] = 0.0;
  df[4][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    df[4][4] += gx[k1]*fx[k1];
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = 2.0*x[3]*gx[k1];
      df[k1][4] = b[k1];
    } 

#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
 //print_matrix(df, 5);
 //printf("fx: %f %f %f gx: %f %f %f\n", fx[0]/calc_norm(fx), fx[1]/calc_norm(fx), 
//	fx[2]/calc_norm(fx), gx[0]/calc_norm(gx), gx[1]/calc_norm(gx), gx[2]/calc_norm(gx));
 //printf("F2BZdistNeg5 fvec (%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]);
#endif
}
extern void fdjacDistNegSE(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3], double *fx, double *gx);

void fdjacDistNeg(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3], double *fx, double *gx)
{
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
  double rDC[3];
  int k1, k2;
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(iA) || is_superellipse(iB))
    {
      fdjacDistNegSE(n, x, fvec, df, vecfunc, iA, iB, shift, fx, gx);
      return;
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = 2.0*Sqr(x[6])*Xb[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3]-rB[k2]);
	}
    }
#if 1
  if (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = x[k1+3] - x[k1];
#ifdef EDHE_FLEX
      for (kk=0; kk < 3; kk++)
	{
	  axi[kk] = typesArr[typeOfPart[iA]].sax[kk];
	  axj[kk] = typesArr[typeOfPart[iB]].sax[kk];
	}
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axi[0],axi[1],axi[2])+max3(axj[0],axj[1],axj[2])))
	{
	  fdjac_disterr = 1;	
	}
#else
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axa[iA],axb[iA],axc[iA])+max3(axa[iB],axb[iB],axc[iB])))
	{
	  fdjac_disterr = 1;	
	}
#endif 
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  for (k1 = 0; k1 < 5; k1++)
    {
      df[3][k1+3] = 0;
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = 0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1+3] = gx[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = 2.0*x[6]*gx[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1 + 2.0*x[7]*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 2.0*x[7]*Xa[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = fx[k1];
#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
}

/* funzione che calcola lo Jacobiano */
void fdjacDist(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3])
{
  double fx[3], gx[3];
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = 2.0*Sqr(x[6])*Xb[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3]-rB[k2]);
	}
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  for (k1 = 0; k1 < 5; k1++)
    {
      df[3][k1+3] = 0;
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = 0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1+3] = gx[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = 2.0*x[6]*gx[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1 + 2.0*Sqr(x[7])*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 2.0*Sqr(x[7])*Xa[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = 2.0*x[7]*fx[k1];
#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*Sqr(x[7]); 
  //MD_DEBUG(printf("F2BZdist fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
}
extern void funcs2beZeroedDistNegSE(int n, double x[], double fvec[], int i, int j, double shift[3]);

void funcs2beZeroedDistNeg(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3];
  /* x = (r, alpha, t) */ 

#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif

#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(i) || is_superellipse(j))
    {
      funcs2beZeroedDistNegSE(n, x, fvec, i, j, shift);
      return;
    }
#endif
 
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
extern void funcs2beZeroedDistNeg5SE(int n, double x[], double fvec[], int i, int j, double shift[3]);

void funcs2beZeroedDistNeg5(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3], rD[3];
  /* x = (r, alpha, t) */ 
  
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(i) || is_superellipse(j))
    {
      funcs2beZeroedDistNeg5SE(n, x, fvec, i, j, shift);
      return;
    }
#endif
 
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
	}
      rD[k1] = x[k1] + fx[k1]*x[4];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(rD[k2] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
void funcs2beZeroedDist(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3];
  /* x = (r, alpha, t) */ 
  
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*Sqr(x[7]); 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
void calc_intersec(double *rB, double *rA, double **Xa, double* rI)
{
  double A, B=0.0, C=0.0, D=0.0, tt=0.0;
  double rBA[3];
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = rB[k1] - rA[k1];
  MD_DEBUG(printf("rBA=(%f,%f,%f)\n", rBA[0], rBA[1], rBA[2])); 
  MD_DEBUG(printf("rB= (%f,%f,%f)\n", rB[0], rB[1], rB[2]));
  A = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	MD_DEBUG(printf("Xa[%d][%d]%f\n", k1, k2,Xa[k1][k2]);)
    	A += rBA[k1]*Xa[k1][k2]*rBA[k2];
      }
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {   
      XarA[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	XarA[k1] += Xa[k1][k2]*rA[k2];
    }
  B = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    B += rBA[k1]*XarA[k1];
  B *= 2.0;
  C = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    C +=  rA[k1]*XarA[k1];
  C = C - 1.0;
  D = Sqr(B) - 4.0*A*C;
  if (D < 0 || A == 0)
    {
      printf("[calc_intersec] Serious problem guessing distance, aborting...\n");
      printf("tt = %f D=%f A=%f B=%f C=%f\n", tt, D, A, B, C);
      printf("distance: %f\n", sqrt(Sqr(rBA[0])+Sqr(rBA[1])+Sqr(rBA[2])));
      print_matrix(Xa,3);
      print_matrix(Xb,3);
      exit(-1);
    }
  tt = (-B + sqrt(D))/(2.0*A); 
#endif
  if (A <= 0)
    {
      printf("[calc_intersec] Serious problem guessing distance, aborting...\n");
      printf("tt = %f D=%f A=%f B=%f C=%f\n", tt, D, A, B, C);
      printf("distance: %f\n", sqrt(Sqr(rBA[0])+Sqr(rBA[1])+Sqr(rBA[2])));
      print_matrix(Xa,3);
      print_matrix(Xb,3);
      exit(-1);
    }
  tt = sqrt(1 / A); 
  for (k1 = 0; k1 < 3; k1++)
    {
      rI[k1] = rA[k1] + tt*rBA[k1];  
    }
}

void guess_dist(int i, int j, 
		double *rA, double *rB, double **Xa, double **Xb, double *rC, double *rD,
		double **RA, double **RB)
{
  double gradA[3], gradB[3], gradaxA[3], gradaxB[3], dA[3], dB[3];
  int k1, n;
  double saA[3], saB[3];
#ifdef EDHE_FLEX
  int typei, typej;
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  saA[0] = typesArr[typei].sax[0];
  saA[1] = typesArr[typei].sax[1];
  saA[2] = typesArr[typei].sax[2];
  saB[0] = typesArr[typej].sax[0];
  saB[1] = typesArr[typej].sax[1];
  saB[2] = typesArr[typej].sax[2];
#else
  saA[0] = axa[i];
  saA[1] = axb[i];
  saA[2] = axc[i];
  saB[0] = axa[j];
  saB[1] = axb[j];
  saB[2] = axc[j];
#endif
  //printf("axes[%d]=%f,%f,%f axes[%d]=%f,%f,%f\n", i, axa[i], axb[i], axc[i], j, axa[j], axb[j], axc[j]);
  for (k1 = 0; k1 < 3; k1++)
    {
      gradA[k1] =  (rB[k1]-rA[k1]);
      gradB[k1] = -(rB[k1]-rA[k1]);
    }
  for (n = 0; n < 3; n++)
    {
      gradaxA[n] = 0;
      gradaxB[n] = 0;
      for (k1 = 0; k1 < 3; k1++) 
	{
	  gradaxA[n] += gradA[k1]*RA[n][k1];
	  gradaxB[n] += gradB[k1]*RB[n][k1];
	}
    }
#if 0 
  gAn = calc_norm(gradaxA);
  gBn = calc_norm(gradaxB);
  for (n = 0; n < 3; n++)
    {
      gradaxA[n] /= gAn;
      gradaxB[n] /= gBn;
    }
#endif
  for (k1=0; k1 < 3; k1++)
    {
      dA[k1] = rA[k1];
      dB[k1] = rB[k1];
      for (n=0; n < 3;n++)
	{
	  dA[k1] += gradaxA[n]*RA[n][k1]*saA[n]/2.0; 
	  dB[k1] += gradaxB[n]*RB[n][k1]*saB[n]/2.0;
	}
    }
  calc_intersec(dA, rA, Xa, rC);
  calc_intersec(dB, rB, Xb, rD);
}

void calc_grad(double *rC, double *rA, double **Xa, double *grad)
{
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      grad[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	grad[k1] += 2.0*Xa[k1][k2]*(rC[k2]-rA[k2]); 
    }
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
extern int check_point(char* msg, double *p, double *rc, double **XX);
extern void distSD(int i, int j, double shift[3], double *vecg, double lambda, int halfspring);
extern void distconjgrad(int i, int j, double shift[3], double *vec);
extern int maxitsRyck;
extern double scalProd(double *A, double *B);
extern int is_superellipse(int i);
extern int newtDistNegSE(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3], int tryagain);

double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
     		double *vecgsup, int calcguess)
{
  /* SDmethod=1 usa la riduzione del passo nello Steepest Descent (SD) e applica lo SD sempre
     SDmethod=2 non usa la riduzione del passo e applica lo SD solo se il calcolo della distanza fallisce 
     SDmethod=3 usa la riduzione del passo e applica lo SD solo se il calcolo della distanza fallisce 
     SDmethod=4 non usa la riduzione del passo e applica lo SD sempre */
  double vecg[8], rC[3], rD[3], rDC[3], r12[3], vecgcg[6], fx[3];
  double ti, segno, segno2;
  double g1=0.0, g2=0.0, SP, nrDC, vecnf[3], nvecnf;
  int retcheck, tryagain = 0;
#ifdef EDHE_FLEX
  int typei, typej;
  double axaiF, axbiF, axciF;
  double axajF, axbjF, axcjF;
#endif
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
  double nf, ng, gradf[3], gradg[3];
  int k1, k2, na, k3;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  ti = t + (t1 - atomTime[i]);
  MD_DEBUG50(printf("[calcDistNeg] - BEGIN");)
#ifdef EDHE_FLEX
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  axaiF = typesArr[typei].sax[0];
  axbiF = typesArr[typei].sax[1];
  axciF = typesArr[typei].sax[2];
  minaxA = min3(axaiF,axbiF,axciF);
  axajF = typesArr[typej].sax[0];
  axbjF = typesArr[typej].sax[1];
  axcjF = typesArr[typej].sax[2];
  minaxB = min3(axajF,axbjF,axcjF);
#else
  minaxA = min3(axa[i],axb[i],axc[i]);
  minaxB = min3(axa[j],axb[j],axc[j]);
#endif
  minaxAB = min(minaxA,minaxB);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */

#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, RtA, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[i];
  invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[i]);
      invbSq[na] = 1/Sqr(axbP[i]);
      invcSq[na] = 1/Sqr(axcP[i]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], RtA);

  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(j, ti, RtB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  UpdateOrient(j, ti, RtB, Omega);
#endif
  na = (j < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typej = typeOfPart[j];
  invaSq[na] = 1.0/Sqr(typesArr[typej].sax[0]);
  MD_DEBUG50(printf("BOH semiasse a del tipo %d: %f\n",typei, typesArr[typej].sax[0]));
  invbSq[na] = 1.0/Sqr(typesArr[typej].sax[1]);
  invcSq[na] = 1.0/Sqr(typesArr[typej].sax[2]);
  MD_DEBUG50(printf("sax of %d %f %f %f %f %f %f\n", i, typesArr[typei].sax[0],typesArr[typei].sax[1],
	 typesArr[typei].sax[2], invaSq[na], invbSq[na], invbSq[na]));
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[j]);
      invbSq[na] = 1/Sqr(axbP[j]);
      invcSq[na] = 1/Sqr(axcP[j]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
#endif
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], RtB);

  if (OprogStatus.dist5)
    {
      for (k1 = 0; k1 < 3; k1++)
    	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      XbXa[k1][k2] = 0;
	      for (k3 = 0; k3 < 3; k3++)
		XbXa[k1][k2] += Xb[k1][k3]*Xa[k3][k2];
	    }
	}
    }
  //printf(">>>>>>> BOHHHHHHHHHHHHHHHHH\n");
#ifdef MC_SIMUL
  /* N.B. se si fa un Monte Carlo di ellissoidi queste condizioni evitano problemi nel calcolo 
   * della distanza qualora ci sia un overlap tale che il centro dell'ellissoide A sia all'interno
   * di B o viceversa. 
   * In tale caso infatti la mossa monte carlo va rigettata e quindi opportunamente la routine ritorna
   * una distanza negativa. */
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (rB[k1]-rA[k1])*Xa[k1][k2]*(rB[k2]-rA[k2]); 
  if (segno < 0)
    return -1.0;
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (rA[k1]-rB[k1])*Xb[k1][k2]*(rA[k2]-rB[k2]); 
  if (segno < 0)
    return -1.0;
#endif
retry:
  if (OprogStatus.forceguess)
    calcguess = 1;
  if (calcguess || tryagain)
    {
      if (OprogStatus.guessDistOpt==1)
	{
	  guess_dist(i, j, rA, rB, Xa, Xb, rC, rD, RtA, RtB);
#if 0
	  for(k1=0; k1 < 3; k1++)
	    r12[k1] = rC[k1]-rD[k1]; 
	  printf("PRIMA PRIMA dist=%.15f\n",calc_norm(r12));
#endif
	}
      else
	{
	  calc_intersec(rB, rA, Xa, rC);
	  calc_intersec(rA, rB, Xb, rD);
	}
#if 1
      for(k1=0; k1 < 3; k1++)
	r12[k1] = rC[k1]-rD[k1]; 
      if ((OprogStatus.SDmethod==1 || OprogStatus.SDmethod==4) || tryagain)
	{
	  for (k1=0; k1 < 3; k1++)
	    {
	      vecgcg[k1] = rC[k1];
	      vecgcg[k1+3] = rD[k1];
	    }
	  //check_point("calcDistNeg", &vecgcg[3], rB, Xb);
#if 0
	  for(k1=0; k1 < 3; k1++)
	    r12[k1] = rC[k1]-rD[k1]; 
	  printf("PRIMA dist=%.15f\n",calc_norm(r12));
	  //printf("distVera=%.15f\n", calcDist(t, i, j, shift, r1, r2, alpha, vecgsup, 1));
#endif
	  distSD(i, j, shift, vecgcg, OprogStatus.springkSD, 1);
	  for (k1=0; k1 < 3; k1++)
	    {
	      rC[k1] = vecgcg[k1];
	      rD[k1] = vecgcg[k1+3];
	    }	 
#endif
#if 0
	{
	  double dist, distVera;
	  for(k1=0; k1 < 3; k1++)
	    r12[k1] = rC[k1]-rD[k1]; 
	  dist = calc_norm(r12);
	  printf("DOPO2 dist: %.15G\n", dist);
	      //printf("dist=%.15f\n",calc_norm(r12));
	}
#endif
	}
      MD_DEBUG50(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		      rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
      calc_grad(rC, rA, Xa, gradf);
      calc_grad(rD, rB, Xb, gradg);
      MD_DEBUG50(printf("gradf=(%f,%f,%f) gradg=(%f,%f,%f)\n",
		      gradf[0], gradf[1], gradf[2], gradg[0], gradg[1], gradg[2]));
      nf = calc_norm(gradf);
      ng = calc_norm(gradg);
      if (OprogStatus.dist5)
	vecg[3] = sqrt(nf/ng);
      else
	vecg[6] = sqrt(nf/ng);
      for (k1=0; k1 < 3; k1++)
	{
	  vecg[k1] = rC[k1];
	  if (!OprogStatus.dist5)
	    vecg[k1+3] = rD[k1];
	  else
	    vecg[k1+5] = rD[k1];
	  rDC[k1] = rD[k1] - rC[k1];
	}
      if (OprogStatus.epsdGDO > 0.0)
	{
	  g1 = calc_norm(rDC)/nf;
	  nrDC = calc_norm(rDC);
	  SP = scalProd(rDC,gradf)/nf;
	  for (k1=0; k1 < 3; k1++)
	    {
	      vecnf[k1] = rDC[k1] - SP*gradf[k1]; 
	    }
	  nvecnf = calc_norm(vecnf);
	  if ( nvecnf > 0.0)
#ifdef EDHE_FLEX
	    g2 = OprogStatus.epsdGDO*min3(axajF,axbjF,axcjF)/calc_norm(vecnf); 
#else
	    g2 = OprogStatus.epsdGDO*min3(axa[j],axb[j],axc[j])/calc_norm(vecnf); 
#endif
	  else 
	    g2 = g1;
	}	  
      if (OprogStatus.dist5)
	{
	  /* N.B. se il guess dei due punti è abbastanza accurato, ossia se si è usato lo SD
	   * è sicuramente meglio dare una stima di vecg[4] altrimenti settarlo a 0 
	   * è un buon guess. */
#if 1
	  if ((OprogStatus.SDmethod==1||OprogStatus.SDmethod==4) || tryagain)
	    {
	      if (scalProd(gradf, rDC) < 0.0)
		vecg[4] = 0.0;
	      else
		vecg[4] = calc_norm(rDC)/nf;  
	    }
	  else
	    {
	      if (OprogStatus.epsdGDO > 0.0)
		{
		  if (scalProd(gradf, rDC) < 0.0)
		    vecg[4] = 0.0;
		  else
		    vecg[4] = min(g1,g2);
		}
	      else
		vecg[4] = 0.0;
	    }
#endif
	//   vecg[4] = 0.0;
	}
      else
	{
#if 1
	  if ((OprogStatus.SDmethod==1 || OprogStatus.SDmethod==4) || tryagain)
	    {
	      if (scalProd(gradf, rDC) < 0.0)
		vecg[7] = 0.0;
	      else
		vecg[7] = calc_norm(rDC)/nf;  
	    }
	  else
	    {
	      if (OprogStatus.epsdGDO > 0.0)
		{
		  if (scalProd(gradf, rDC) < 0.0)
		    vecg[7] = 0.0;
		  else
		    vecg[7] = min(g1,g2);
		}
	      else
		vecg[7] = 0.0;
	    }
#endif
	  //vecg[7] = 0.0;
	}
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  
  //printf(">>>>>>>>>>>> INIZIO NR\n");
  MD_DEBUG(printf(">>>>>>> alpha: %f beta: %f\n", vecg[6], vecg[7]));
  if (OprogStatus.dist5)
    {
      double vecg8[8];
      if ((tryagain && OprogStatus.dist8stps > 0) || OprogStatus.dist8stps < 0)
	{
	  OprogStatus.dist5 = 0;
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      vecg8[k1] = vecg[k1];
	      vecg8[k1+3] = vecg[k1+5];
	    }
	  vecg8[6] = vecg[3];
	  vecg8[7] = vecg[4];
	  newtDistNeg(vecg8, 8, &retcheck, funcs2beZeroedDistNeg, i, j, shift, 0); 
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      vecg[k1] = vecg8[k1];
	      vecg[k1+5] = vecg8[k1+3];
	    }
	  vecg[3] = vecg8[6];
	  vecg[4] = vecg8[7];
	  OprogStatus.dist5 = 1;
	}
#ifdef MD_SUPERELLIPSOID
      if (is_superellipse(i) || is_superellipse(j))
	newtDistNegSE(vecg, 5, &retcheck, funcs2beZeroedDistNeg5, i, j, shift, tryagain); 
      else
	newtDistNeg(vecg, 5, &retcheck, funcs2beZeroedDistNeg5, i, j, shift, tryagain); 
#else
      newtDistNeg(vecg, 5, &retcheck, funcs2beZeroedDistNeg5, i, j, shift, tryagain); 
#endif
    }
  else 
    newtDistNeg(vecg, 8, &retcheck, funcs2beZeroedDistNeg, i, j, shift, tryagain); 
  numcalldist++;
  if (retcheck != 0)
    {
      if (tryagain && OprogStatus.targetPhi <= 0.0)
	{
	  printf("[ERROR] I'm sorry but I can't really calculate distance\n");
	  exit(-1);
     	}
	 
      if (!tryagain && (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod==3))
	{
	  numdisttryagain++;
	  tryagain = 1; 
	  goto retry;
	}
     if (OprogStatus.targetPhi > 0)
	{
	  calcdist_retcheck=1;
	  return 0.0;
	} 
      if (calcguess==0)
	{
	  calcguess=2;
	  goto retry;
	} 
      printf("I couldn't calculate distance between %d and %d\n, exiting....\n", i, j);
      Oparams.time = t + t1;
      store_bump(i, j);
      exit(-1);
    }
  if (!OprogStatus.dist5)
    {
      for (k1 = 0; k1 < 8; k1++)
	{
	  vecgsup[k1] = vecg[k1]; 
	}
    }
  else
    {
      for (k1 = 0; k1 < 5; k1++)
	{
	  vecgsup[k1] = vecg[k1]; 
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      if (!OprogStatus.dist5)
      	r2[k1] = vecg[k1+3];
#if 0
      if (OprogStatus.dist5 && tryagain)
	r2[k1] = vecg[k1+5];
#endif
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      if (OprogStatus.dist5)// && !tryagain)
	{
    	  fx[k1] = 0;
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      fx[k1] += 2.0*Xa[k1][k2]*(r1[k2]-rA[k2]);
	    }
	  r2[k1] = r1[k1] + fx[k1]*vecg[4];
	  vecgsup[k1+5] = r2[k1];
	}
      r12[k1] = r1[k1] - r2[k1];
    } 
  
  *alpha = vecg[3];
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (r2[k1]-rA[k1])*Xa[k1][k2]*(r2[k2]-rA[k2]); 
#if 0
  if (!OprogStatus.dist5)
    {
      if (segno*vecg[7]<0 && fabs(segno*vecg[7])>3E-8)
	{
	  if (OprogStatus.targetPhi>0)
	    {
	      calcdist_retcheck = 1;
	      return 0.0;
	    }
	  printf("segno: %.8G vecg[7]: %.8G dist=%.15G\n", segno, vecg[4], calc_norm(r12));
	  exit(-1);
	  //return calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
	  //exit(-1);
	}
    }
#endif

  MD_DEBUG50(printf("[calcDistNeg] - END");)
  if (OprogStatus.dist5)
    {
      if (segno*vecg[4]<0 && fabs(segno*vecg[4])>3E-8)
	{
	  if (tryagain && OprogStatus.targetPhi <= 0.0)
    	    {
	      printf("[ERROR] I'm sorry but I can't really calculate distance\n");
	      exit(-1);
	    }
	  if (!tryagain && (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod==3))
	    {
	      numdisttryagain++;
	      tryagain = 1; 
	      goto retry;
	    }
	  if (OprogStatus.targetPhi>0)
	    {
	      calcdist_retcheck = 1;
	      return 0.0;
	    } 
	  printf("segno: %.8G vecg[7]: %.8G dist=%.15G\n", segno, vecg[4], calc_norm(r12));
	  return calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
	  //exit(-1);
	}
    }
  else
    {
      if (segno*vecg[7]<0 && fabs(segno*vecg[7])>3E-8)
	{
	  
      	  if (tryagain && OprogStatus.targetPhi <= 0)
	    {
	      printf("[ERROR] I'm sorry but I can't really calculate distance\n");
	      exit(-1);
	    } 
	  if (!tryagain && ( OprogStatus.SDmethod==2 || OprogStatus.SDmethod==3 ))
	    {
	      numdisttryagain++;
	      tryagain = 1; 
	      goto retry;
	    }
	  if (OprogStatus.targetPhi>0)
	    {
	      calcdist_retcheck = 1;
	      return 0.0;
	    }
	  printf("segno: %.8G vecg[7]: %.8G dist=%.15G\n", segno, vecg[7], calc_norm(r12));
	  return calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
	}
    }
#ifdef MD_PATCHY_HE
  segno2 = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno2 += (r1[k1]-rB[k1])*Xb[k1][k2]*(r1[k2]-rB[k2]); 
  if (segno2*segno < 0.0 && fabs(segno*segno2) > 3E-8)
    {
      if (tryagain && OprogStatus.targetPhi <= 0)
	{
	  printf("[ERROR segno*segno2] I'm sorry but I can't really calculate distance\n");
	  exit(-1);
	} 
      if (!tryagain && ( OprogStatus.SDmethod==2 || OprogStatus.SDmethod==3 ))
	{
	  numdisttryagain++;
	  tryagain = 1; 
	  goto retry;
	}
      if (OprogStatus.targetPhi>0)
	{
	  calcdist_retcheck = 1;
	  return 0.0;
	}
      printf("segno: %.8G segno2: %.15G dist=%.15G\n", segno, segno2, calc_norm(r12));
      return calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
    }
#endif
  if (segno > 0)
    return calc_norm(r12);
  else
    return -calc_norm(r12);
}

double calcDist(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess)
{
  double vecg[8], rC[3], rD[3], rDC[3], r12[3];
  double ti;
  int retcheck;
  double Omega[3][3], nf, ng, gradf[3], gradg[3];
  int k1, na;
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
  UpdateOrient(i, ti, Rt, Omega);
  na = (i < Oparams.parnumA)?0:1;
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, Rt, Omega);
  na = (j < Oparams.parnumA)?0:1;
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
  if (calcguess)
    {
      calc_intersec(rB, rA, Xa, rC);
      calc_intersec(rA, rB, Xb, rD);
      MD_DEBUG(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		      rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
      calc_grad(rC, rA, Xa, gradf);
      calc_grad(rD, rB, Xb, gradg);
      MD_DEBUG(printf("gradf=(%f,%f,%f) gradg=(%f,%f,%f)\n",
		      gradf[0], gradf[1], gradf[2], gradg[0], gradg[1], gradg[2]));
      nf = calc_norm(gradf);
      ng = calc_norm(gradg);
      
      vecg[6] = sqrt(nf/ng);
      for (k1=0; k1 < 3; k1++)
	{
	  vecg[k1] = rC[k1];
	  vecg[k1+3] = rD[k1];
	  rDC[k1] = rD[k1] - rC[k1];
	}
      vecg[7] = sqrt(calc_norm(rDC)/nf);  
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  MD_DEBUG(printf("alpha: %f beta: %f\n", vecg[6], vecg[7]));
  newtDist(vecg, 8, &retcheck, funcs2beZeroedDist, i, j, shift); 
  if (retcheck != 0)
    {
      printf("[calcDist] I couldn't calculate distance between %d and %d\n, exiting....\n", i, j);
      //Oparams.time = t;
      //store_bump(i, j);
      //exit(-1);
    }
  for (k1 = 0; k1 < 8; k1++)
    {
      vecgsup[k1] = vecg[k1]; 
    }  
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      r2[k1] = vecg[k1+3];
      r12[k1] = r1[k1] - r2[k1];
    }
  *alpha = vecg[6];
#if 0
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (r2[k1]-rA[k1])*Xa[k1][k2]*(r2[k2]-rA[k2]); 
#endif
#if 1
  return calc_norm(r12);
#if 0
  if (segno > 0)
    return calc_norm(r12);
  else
    return -calc_norm(r12);
#endif
#else
#if 0
#if MD_DEBUG(x) == x
  for (k1 = 0; k1 < 3; k1++)
    rDC[k1] = r1[k1] - r2[k1];
    {
      FILE *f;
      f =fopen("dist.dat","a");
      fprintf(f,"%.15f %.15f\n", Oparams.time, (segno>0)?calc_norm(rDC):-calc_norm(rDC));  
      fclose(f);
    }
  return calc_norm(rDC);
#endif
#endif
#endif
}
void rebuildCalendar(void);
int vc_is_pos(int i, int j, double rCx, double rCy, double rCz,
	       double t)
{ 
  double rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double norm[3], wrx, wry, wrz, OmegaA[3][3], OmegaB[3][3];
#ifdef EDHE_FLEX
  int typei, typej;
#endif
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
  
  UpdateOrient(i, t-atomTime[i], RA, OmegaA);

#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[i];
  invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif  
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], RA);
  na = (j < Oparams.parnumA)?0:1;
  UpdateOrient(j, t-atomTime[j], RB, OmegaB);
#ifdef EDHE_FLEX
  na = 0;
  typej = typeOfPart[j];
  invaSq[na] = 1/Sqr(typesArr[typej].sax[0]);
  invbSq[na] = 1/Sqr(typesArr[typej].sax[1]);
  invcSq[na] = 1/Sqr(typesArr[typej].sax[2]);
#else
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
double calcvecF(int i, int j, double t, double t1, double *r1, double *r2, double* ddot, double shift[3])
{
  int kk;
  double rcat[3], rdbt[3], wra[3], wrb[3], dti, dtj;
#ifdef MD_DDOT_OPT
  double normrcd, r12[3], sp;
#endif
  //evolveVec(i, t-t1, rcat, r1);
  //evolveVec(j, t-t1, rdbt, r2);
  dti = t1 - atomTime[i];
  dtj = t1 - atomTime[j];
  rcat[0] = r1[0] - (rx[i] + vx[i]*(t+dti)); 
  rcat[1] = r1[1] - (ry[i] + vy[i]*(t+dti));
  rcat[2] = r1[2] - (rz[i] + vz[i]*(t+dti));
  rdbt[0] = r2[0] - (rx[j] + vx[j]*(t+dtj))-shift[0]; 
  rdbt[1] = r2[1] - (ry[j] + vy[j]*(t+dtj))-shift[1];
  rdbt[2] = r2[2] - (rz[j] + vz[j]*(t+dtj))-shift[2];
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
int refine_contact(int i, int j, double t1, double t, double vecgd[8], double shift[3],double  vecg[5])
{
  int kk, retcheck;

  for (kk = 0; kk < 3; kk++)
    {
      if (OprogStatus.dist5)
	vecg[kk] = (vecgd[kk]+vecgd[kk+5])*0.5; 
      else
	vecg[kk] = (vecgd[kk]+vecgd[kk+3])*0.5; 
    }
  if (OprogStatus.dist5)
    vecg[3] = vecgd[3];
  else
    vecg[3] = vecgd[6];
  vecg[4] = t-t1;
  trefG = t1;
  newt(vecg, 5, &retcheck, funcs2beZeroed, i, j, shift); 
  vecg[4] += t1;
  if (retcheck==2)
    {
      MD_DEBUG31(printf("newt did not find any contact point!\n"));
      return 0;
    }
#if 0
  else if (vecg[4] < Oparams.time ||
	   fabs(vecg[4] - Oparams.time)<1E-12 )
    {
      /* se l'urto è nel passato chiaramente va scartato
       * tuttavia se t è minore di zero per errori di roundoff? */
      /* Notare che i centroidi si possono overlappare e quindi t può
       * essere tranquillamente negativo */
      MD_DEBUG(printf("i=%d j=%d <<< vecg[4]=%.15f time:%.15f\n",
		      i, j, vecg[4], Oparams.time));
      return 0;
    }
#endif
  else
    {
#if 0
      if (!vc_is_pos(i, j, vecg[0], vecg[1], vecg[2], vecg[4]))
	{
	  MD_DEBUG(printf("vc is positive!\n"));
	  MD_DEBUG(printf("t=%.15f collision predicted %d-%d\n",
			  Oparams.time, i, j));
	  return 0;
	}
#endif
      return 1; 
    }
}
#ifdef MD_ASYM_ITENS
double calcopt_maxddot(int i, int j, double *r1 , double *r2, double factori, double factorj)
{
  int kk;
  double Iamin, Ibmin;
  double dd[3], ndd;
#ifdef EDHE_FLEX
  int typei, typej;
#else
  int na;
#endif
  for (kk=0; kk < 3; kk++)
    dd[kk] = r2[kk] - r1[kk];
  ndd = calc_norm(dd);
#ifdef EDHE_FLEX
  typei = typeOfPart[i];
  Iamin = min(typesArr[typei].I[0], typesArr[typei].I[1]);
  typej = typeOfPart[j];
  Ibmin = min(typesArr[typej].I[0], typesArr[typej].I[0]);
#else
  na = i<Oparams.parnumA?0:1;
  Iamin = min(Oparams.I[na][0],Oparams.I[na][2]);
  na = j<Oparams.parnumA?0:1;
  Ibmin = min(Oparams.I[na][0],Oparams.I[na][2]);
#endif
  for (kk=0; kk < 3; kk++)
    dd[kk] /= ndd;
  return (vx[i]-vx[j])*dd[0]+(vy[i]-vy[j])*dd[1]+(vz[i]-vz[j])*dd[2] +
	angM[i]*factori/Iamin + angM[j]*factorj/Ibmin;
}
#else
double calcopt_maxddot(int i, int j, double *r1 , double *r2, double factori, double factorj)
{
  int kk;
  double dd[3], ndd;
  for (kk=0; kk < 3; kk++)
    dd[kk] = r2[kk] - r1[kk];
  ndd = calc_norm(dd);
  for (kk=0; kk < 3; kk++)
    dd[kk] /= ndd;
  return (vx[i]-vx[j])*dd[0]+(vy[i]-vy[j])*dd[1]+(vz[i]-vz[j])*dd[2] +
	sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori
	+ sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;
}
#endif
#ifdef MD_ASYM_ITENS
double calc_maxddot(int i, int j);
#endif
int search_contact_faster(int i, int j, double *shift, double *t, double t1, double t2, double *vecgd, double epsd, double *d1, double epsdFast, double *r1, double *r2)
{
  /* NOTA: 
   * MAXOPTITS è il numero massimo di iterazioni al di sopra del quale esce */
  double maxddot, told, delt, normddot, ddot[3];
  const int MAXOPTITS = 500;
  double alpha;
#ifndef MD_ASYM_ITENS
  double factori, factorj;
#endif
  int its=0; 
    
  /* estimate of maximum rate of change for d */
#ifdef MD_ASYM_ITENS
  maxddot = calc_maxddot(i, j);
#else
  factori = 0.5*maxax[i]+OprogStatus.epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
  factorj = 0.5*maxax[j]+OprogStatus.epsd;//sqrt(Sqr(axa[j])+Sqr(axb[j])+Sqr(axc[j]));
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;
#endif
  *d1 = calcDistNeg(*t, t1, i, j, shift, r1, r2, &alpha, vecgd, 1);
  timesF++;
  MD_DEBUG31(printf("Pri distances between %d-%d d1=%.12G epsd*epsdTimes:%f\n", i, j, *d1, epsdFast));
  told = *t;
  while (*d1 > epsdFast && its < MAXOPTITS)
    {
#if 1
      if (*t+t1 < t2 && (t2 - (*t + t1))*maxddot < fabs(*d1) - OprogStatus.epsd)
	{
	  return 1;
	}
#endif
#if 0
      factori = 0.5*maxax[i]+OprogStatus.epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
      factorj = 0.5*maxax[j]+OprogStatus.epsd;//sqrt(Sqr(axa[j])+Sqr(axb[j])+Sqr(axc[j]));
      maxddot = calcopt_maxddot(i, j, r1, r2, factori, factorj);
      if (maxddot < 0)
	return 1;
#endif 
      delt = *d1  / maxddot;
#if defined(MD_BASIC_DT) || defined (EDHE_FLEX)
      if (delt < (epsd / maxddot))
	{
	  MD_DEBUG10(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
#else
      normddot = calcvecF(i, j, *t, t1, r1, r2, ddot, shift);
      /* check for convergence */
     
      if (normddot!=0 && delt < (epsd / normddot))
	{
	  MD_DEBUG10(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
#endif
      *t += delt;
#if 1
      if (*t + t1 > t2)
	{
	  *t = told;
	  MD_DEBUG31(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG31(printf("convergence t>t2\n"));
	  *d1 = calcDistNeg(*t, t1, i, j, shift, r1, r2, &alpha, vecgd, 1);
	  return 1;
	}
#endif
      *d1 = calcDistNeg(*t, t1, i, j, shift, r1, r2, &alpha, vecgd, 1);
      if (*d1 < 0)
	{
	  /* go back! */
	  MD_DEBUG31(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG31(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNeg(*t, t1, i, j, shift, r1, r2, &alpha, vecgd, 1);
	  return 0;
	}
      told = *t;
      its++;
      itsF++;
    }

  MD_DEBUG31(printf("max iterations %d iterations reached t=%f t2=%f\n", its, *t, t2));
  return 0;
}
extern double **Aip;
extern void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double zbrent(double (*func)(double), double x1, double x2, double tol);
int polinterr, polinterrRyck;
double xa[3], ya[3];
double sogliaErr_zbrent;
double distfunc(double x)
{
  double dy, y;
  polint(xa, ya, 3, x, &y, &dy);
  if (polinterr==1)
    return 0.0;
  if (dy > sogliaErr_zbrent)
    {
      if (OprogStatus.phitol <= 0)
	printf("dy=%.15G\n", dy);
      polinterr = 1;
    }
  else 
    polinterr = 0;
  return y;
}
int interpolSNP(int i, int j, double tref, double t, double delt, double d1, double d2, double *tmin, double* vecg, double shift[3])
{
  double d3, A;
  double r1[3], r2[3], alpha;
  double dmin;
  d3 = calcDistNeg(t+delt*0.5, tref, i, j, shift, r1, r2, &alpha, vecg, 0);
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
  if (*tmin < t+delt && *tmin > t)
    {
      dmin = calcDistNeg(*tmin, tref, i, j, shift, r1, r2, &alpha, vecg, 0);
      *tmin += tref;
      return 0;
    }
  return 1;
}
int interpol(int i, int j, double tref, double t, double delt, double d1, double d2, double *troot, double* vecg, double shift[3], int bracketing)
{
  int nb;
  double d3, t1, t2, A;
  double r1[3], r2[3], alpha, xb1[2], xb2[2];
  double tmin, dmin;
  d3 = calcDistNeg(t+delt*0.5, tref, i, j, shift, r1, r2, &alpha, vecg, 0);
#if 0
  if (d1 > OprogStatus.epsd)
    {
      printf("d1=%.15G t=%.15G\n", d1, t);
      exit(-1);
    }
#endif
  xa[0] = t;
  ya[0] = d1;
  xa[1] = t+delt*0.5;
  ya[1] = d3;
  xa[2] = t+delt;
  ya[2] = d2;
  //printf("(%.8f,%.8f) (%.8f,%.8f) (%.8f,%.8f)\n", t, d1, t+delt, d2, t+delt*0.5, d3);
  //printf("{%.8f %.8f %.8f}\n", bip[0], bip[1], bip[2]);
  //print_matrix(Aip, 3);

  //printf("d1:%.10f d2: %.10f\n", d1, d2); 
  //printf("polint1: %.10f polint2: %.10f\n",distfunc(t), distfunc(t+delt));
  polinterr = 0;
  sogliaErr_zbrent = OprogStatus.epsd;
  if (!bracketing)
    {
      t1 = t;
      t2 = t+delt;
    }
  else
    {
      if (OprogStatus.zbrakn==0)
	{
	  if (ya[0]-ya[1] == 0.0)
	    {
	      tmin = t + delt*0.25;
	    }
	  else if (ya[2]-ya[0] ==0.0)
	    {
	      tmin = t + delt*0.5;
	    }
	  else
	    {      
	      A = (ya[2]-ya[0])/(ya[0]-ya[1]);
	      tmin = t + 0.5*delt*((1.0 + A * 0.25)/( 1.0 + A * 0.5));
	    }
	  if (tmin < t+delt && tmin > t)
	    {
	      dmin = calcDistNeg(tmin, tref, i, j, shift, r1, r2, &alpha, vecg, 0);
	      if (d1*dmin < 0.0)
		{
		  t2 = tmin;
		  t1 = t;
		}
	      else
		return 1;
	    }
	  else
	    {
	      return 1;
	    }
	}
      else
	{
	  t1 = t;
	  t2 = t+delt;
	  nb = 1;
	  zbrak(distfunc, t1, t2, OprogStatus.zbrakn, xb1, xb2, &nb);
	  //printf("nb=%d t2=%.15G tmin=%.15G dmin=%.15G d1=%.15G\n", nb, t2, tmin, dmin, d1);
	  //printf("(%.15G,%.15G)-(%.15G,%.15G)-(%.15G,%.15G)\n", t, d1, t+delt*0.5, d3, t+delt, d2);
	  if (nb==0 || polinterr==1)
	    {
	      return 1;
	    }
	  t1 = xb1[0];
	  t2 = xb2[0];
	}
    }
  if (polinterr)
    return 1;
  if (OprogStatus.zbrentTol <= 0.0)
    {
      *troot = tref + (t1+t2)*0.5;
      return 0;
    }
  *troot=zbrent(distfunc, t1, t2, OprogStatus.zbrentTol);
  if (polinterr)
    {
      printf("[interpol] bracketing: %d polinterr=%d t1=%.15G t2=%.15G\n", bracketing,polinterr, t1, t2);
      printf("d: %.15G,%.15G,%.15G\n", d1, d3, d2);
      printf("t: %.15G,%.15G,%.15G\n", t, t+delt*0.5, t+delt);
      printf("distfunc(t1)=%.10G distfunc(t2)=%.10G\n", distfunc(t), distfunc(t+delt));
      return 1;
    }
  if ((*troot < t && fabs(*troot-t)>3E-8) || (*troot > t+delt && fabs(*troot - (t+delt))>3E-8))
    {
      if (OprogStatus.zbrakn > 0)
	printf("[interpol] brack: %d xb1: %.10G xb2: %.10G\n", bracketing, xb1[0], xb2[0]);
      printf("*troot: %.15G t=%.15G t+delt:%.15G\n", *troot, t, t+delt);
      printf("d1=%.10G d2=%.10G d3:%.10G\n", d1, d2, d3);
      printf("distfunc(t1)=%.10G distfunc(t2)=%.10G\n", distfunc(t), distfunc(t+delt));
      printf("distfunc(t+delt*0.5)=%.10G\n", distfunc(t+delt*0.5));
      return 1;
    }
  calcDistNeg(*troot, tref, i, j, shift, r1, r2, &alpha, vecg, 0);
  //printf("t=%.8G t+delt=%.8G troot=%.8G\n", t, t+delt, *troot);
  *troot += tref;
  return 0;
}
#ifdef MD_ASYM_ITENS
double calc_maxddot(int i, int j)
{
#if 0
  int na;
  double Iamin, Ibmin;
#endif
  double factori, factorj;
  factori = 0.5*maxax[i]+OprogStatus.epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
  factorj = 0.5*maxax[j]+OprogStatus.epsd;//sqrt(Sqr(axa[j])+Sqr(axb[j])+Sqr(axc[j]));
  MD_DEBUG20(printf("[calc_maxddot] maxax[%d]:%.15G maxax[%d]:%.15G\n", i, maxax[i], j, maxax[j]);)
#if 0
  /* N.B. nel caso della trottola simmetrica il modulo di w comunque non varia (anche 
   * se w ruota, ossia precede, intorno al vettore momento angolare) nel 
   * tempo per questo la maggiorazione rimane inalterata */
  na = i<Oparams.parnumA?0:1;
  Iamin = min(Oparams.I[na][0],Oparams.I[na][2]);
  na = j<Oparams.parnumA?0:1;
  Ibmin = min(Oparams.I[na][0],Oparams.I[na][2]);
  return sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    angM[i]*factori/Iamin + angM[j]*factorj/Ibmin;
#else
  return sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori + 
    sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;
#endif
}
#endif
#ifdef EDHE_FLEX
int locate_contact_HS(int i, int j, double shift[3], double t1, double t2, double vecg[5])
{
  double d, t, evtime;
  int typei, typej;
  double sigSq, b, dr[3], dv[3], tInt, vv;
  int collCode;
  /* calcola cmq l'urto fra le due core spheres */
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  sigSq = Sqr(typesArr[typei].sax[0]+typesArr[typej].sax[0]); 
  tInt = Oparams.time - atomTime[j];
  dr[0] = rx[i] - (rx[j] + vx[j] * tInt) - shift[0];	  
  dv[0] = vx[i] - vx[j];
  dr[1] = ry[i] - (ry[j] + vy[j] * tInt) - shift[1];
  dv[1] = vy[i] - vy[j];
  dr[2] = rz[i] - (rz[j] + vz[j] * tInt) - shift[2];
  dv[2] = vz[i] - vz[j];
  b = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
  collCode = MD_EVENT_NONE;
  evtime = timbig;
  if (b < 0.0) 
    {
      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
      d = Sqr (b) - vv * 
	(Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]) - sigSq);
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
  if (collCode != MD_EVENT_NONE && evtime > t1 && evtime < t2)
    {
      /* N.B. nel caso di urto fra sfere dure le coordinate del punto di contatto non servono,
	 per questo qui vengono messe a 0 
       */
      vecg[0] = 0;
      vecg[1] = 0;
      vecg[2] = 0;
      vecg[3] = 0;
      vecg[4] = evtime;
      return 1;  
    }
  else 
    return 0;
}
#endif
int locate_contact(int i, int j, double shift[3], double t1, double t2, double vecg[5])
{
  double h, d, dold, alpha, vecgd[8], vecgdold[8], t, r1[3], r2[3]; 
  double maxddot, delt, troot, vecgroot[8], tini, tmin;
#ifndef MD_BASIC_DT
  double normddot, ddot[3], dold2, vecgdold2[8];
#endif
  //const int MAXOPTITS = 4;
  double epsd, epsdFast, epsdFastR, epsdMax;
#ifndef MD_ASYM_ITENS
  double factori, factorj; 
#endif
  int dorefine, sumnegpairs=0;
  int its, foundrc, kk;
#ifdef EDHE_FLEX
  if (typesArr[typeOfPart[i]].ignoreCore || typesArr[typeOfPart[j]].ignoreCore)
    return 0;
  if (are_spheres(i, j))
    {
      return locate_contact_HS(i, j, shift, t1, t2, vecg);
    }
#endif
  epsd = OprogStatus.epsd;
  epsdFast = OprogStatus.epsdFast;
  epsdFastR= OprogStatus.epsdFastR;
  epsdMax = OprogStatus.epsdMax;

  /* NOTA: 
   * - epsd è di quanto varia d ad ogni iterazione e quindi determina il grado di accuratezza
   * con cui viene individuato il punto di contatto. In generale se due ellissoidi si "spizzicano"
   * ad una distanza minore di epsd tali urti non vengono rilevati 
   * - epsdTimes*epsd è la soglia sotto la quale la ricerca veloce fatta in search_contact_faster termina 
   * - epsdTimesIsteresi*epsd è la soglia al di sopra della quale viene di nuovo usata search_contact_faster
   *   Tale valore è bene che sia abbastanza grande di modo che non faccia continue ricerche veloci e lente 
   *   che rallentano moltissimo. Infatti tale ricerca veloce serve solo nel caso in cui due ellissoidi si 
   *   sfiorano per poi allontanarsi. 
   */
  t = 0;//t1;

#ifdef MD_ASYM_ITENS
  maxddot = calc_maxddot(i, j);  
#else
  factori = 0.5*maxax[i]+OprogStatus.epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
  factorj = 0.5*maxax[j]+OprogStatus.epsd;//sqrt(Sqr(axa[j])+Sqr(axb[j])+Sqr(axc[j]));
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;
#endif
  MD_DEBUG10(printf("[locate_contact] %d-%d t1=%f t2=%f shift=(%f,%f,%f)\n", i,j,t1, t2, shift[0], shift[1], shift[2]));
  h = OprogStatus.h; /* last resort time increment */
#if 0
  if (lastbump[i]==j && lastbump[j]==i)
    {
      t += h;
      MD_DEBUG(printf("last collision was between %d-%d\n",i,j));
      MD_DEBUG(printf("atomTime[%d]:%.15G atomTime[%d]:%.15G\n", i, atomTime[i], j, atomTime[j])); 
    }
#endif
#if 0
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j];
  d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 1);
  MD_DEBUG(printf("Pri distances between %d-%d d1=%.12G", i, j, d1));
  told = t;
  its = 0;
  while (d1 > 10*epsd && its < MAXOPTITS)
    {
      delt = d1 / maxddot;
      t += delt;
      d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 1);
      if (d1 < 0)
	{
	  /* go back! */
	  t = told;	  
	  d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 1);
	}
      told = t;
      its++;
    }
#else
  //d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 1);
  timesS++;
#endif
  MD_DEBUG(printf("Dopo distances between %d-%d d1=%.12G", i, j, d));
#if 1
  MD_DEBUG20(printf("[LOCATE_CONTACT] INIZIO\n"));
  d = calcDistNeg(t, t1, i, j, shift, r1, r2, &alpha, vecgd, 1);
#if 1
#if 0
  maxddot = calcopt_maxddot(i, j, r1, r2, factori, factorj);
  if (maxddot < 0)
    return 0;
#endif
  if ((t2 - t1)*maxddot < fabs(d) - OprogStatus.epsd)
    {
      return 0;
    }
#endif
#ifdef MD_PATCHY_HE
  if (lastbump[j].mol==i && lastbump[j].at==0 && lastbump[i].mol==j && lastbump[i].at==0)
#else
  if (lastbump[j]==i && lastbump[i]==j)
#endif
    {
      MD_DEBUG11(printf("last collision was between (%d-%d) d=%.15G\n", i, j, d));
      if (d < 0 && fabs(d) > epsd)
	{
	  printf("i=%d j=%d BOH d=%.15G\n", i, j, d);
	  exit(-1);
	}
      sumnegpairs = 1;
#if 0
      its = 0;	
      while (d < 0)
	{
	  //printf("===>its=%d d=%.15G d=%.15G\n", its, fabs(d)/maxddot, d);
	  t += h;
	  its++;
	  if (t + t1 > t2)
	    return 0;
	  d = calcDistNeg(t, t1, i, j, shift, r1, r2, &alpha, vecgd, 1);
	}
#endif
    }
  else if (d<0&&fabs(d)>OprogStatus.epsd)
    {
      printf("[WARNING] t=%.10G d=%.15G < 0 i=%d j=%d\n", t+t1, d, i, j);
      printf("[WARNING] Some collision has been missed, ellipsoid may overlap!\n");
      store_bump(i, j);
      return 0;
    }
#endif
  MD_DEBUG31(printf("[LOCATE_CONTACT] prima search contact faster\n"));
  if (search_contact_faster(i, j, shift, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2))
    return 0; 
  MD_DEBUG31(printf("[LOCATE_CONTACT]>>>>d:%f t=%.15G\n", d,t));
  foundrc = 0;
#if 0
  if (d1 < 0)
    {
      if (refine_contact(i, j, t, vecgd, shift, vecg))
	{
	  MD_DEBUG(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
	  return 1;
	}
      else 
	{
	  MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
	  if (d2 < 0)
	    {
	      MD_DEBUG(printf("d2 < 0 and I did not find contact point, boh...\n"));
	      return 0;
	    }
	}
    }
#endif
  for (kk = 0; kk < 8; kk++)
    vecgdold[kk] = vecgd[kk];
  dold = d;
  its = 0;
  while (t + t1 < t2)
    {
#ifdef MD_BASIC_DT
      delt = epsd/maxddot;
      tini = t;
      t += delt;
      d = calcDistNeg(t, t1, i, j, shift, r1, r2, &alpha, vecgd, 0);
#else
#if 1
      normddot = calcvecF(i, j, t, t1, r1, r2, ddot, shift);
      if (normddot!=0)
	delt = epsd/normddot;
#endif
      else
	delt = epsd / maxddot;
	//delt = h;
      if (dold < epsd)
	delt = epsd / maxddot;
      tini = t;
      t += delt;
      for (kk = 0; kk < 8; kk++)
	vecgdold2[kk] = vecgd[kk];
      dold2 = dold;
      d = calcDistNeg(t, t1, i, j, shift, r1, r2, &alpha, vecgd, 0);
      if (fabs(d-dold2) > epsdMax)
	{
	  /* se la variazione di d è eccessiva 
	   * cerca di correggere il passo per ottenere un valore
	   * più vicino a epsd*/
	  //printf("P delt: %.15G d2-d2o:%.15G d2:%.15G d2o:%.15G\n", delt, fabs(d2-d2old), d2, d2old);
	  t -= delt;
	  //delt = d2old / maxddot;
	  delt = epsd / maxddot;
#if 0
	  if (delt < h)
	    delt = h;
#endif
	  t += delt; 
	  //t += delt*epsd/fabs(d2-d2old);
	  itsS++;
	  d = calcDistNeg(t, t1, i, j, shift, r1, r2, &alpha, vecgdold2, 0);
	  for (kk = 0; kk < 8; kk++)
	    vecgd[kk] = vecgdold2[kk];
	  //printf("D delt: %.15G d2-d2o:%.15G d2:%.15G d2o:%.15G\n", delt*epsd/fabs(d2-d2old), fabs(d2-d2old), d2, d2old);
	}
#endif
      MD_DEBUG(printf(">>>> t = %f d1:%f d2:%f d1-d2:%.15G\n", t, d1, d2, fabs(d1-d2)));
      if (sumnegpairs)
	{
	  MD_DEBUG(printf("sumnnegpairs d=%.15G\n", d));
  	  if (d <= 0.0)
	    {
	      if(!interpolSNP(i, j, t1, tini, delt, dold, d, &tmin, vecgd, shift))
		{
		  tmin -= t1;
		  delt = tmin - tini;
		  t = tmin;
		  d = calcDistNeg(t, t1, i, j, shift, r1, r2, &alpha, vecgd, 0);
		  MD_DEBUG(printf("qui tmin = %.15G d=%.15G delt=%.15G\n", tmin, d, delt));
		}
	    }
	  sumnegpairs = 0;	  
	}
      dorefine = 0;      
      if (dold > 0 && d < 0)
	{
#if 0
	  if (d2 <0)
	    {
	      t -= h;
	      d2 = calcDist(t, i, j, shift, r1, r2, &alpha, vecgd, 0);
	    }
#endif
       	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
#ifndef MD_NOINTERPOL  
	 if (interpol(i, j, t1, t-delt, delt, dold, d, &troot, vecgroot, shift, 0))
#endif
	    {
	      /* vecgd2 è vecgd al tempo t-delt */
	      for (kk=0; kk < 8; kk++)
		vecgroot[kk] = vecgdold[kk];
	      /* forse è meglio scegliere il valore di t più grande per ridurre il rischio
	       * di eventi coincidenti! */
	      /* VECCHIA SOLUZIONE: troot = t + t1 - delt;*/
	      troot = t + t1;
	    }
	  dorefine = 1;
	}
      else if (dold < OprogStatus.epsd && d < OprogStatus.epsd)
	{
#ifndef MD_NOINTERPOL
	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
	  
	  if (interpol(i, j, t1, t-delt, delt, dold, d, &troot, vecgroot, shift, 1))
	    dorefine = 0;
	  else 
	    dorefine = 1;
#endif
	}
      if (dorefine)
	{
	  if (refine_contact(i, j, t1, troot, vecgroot, shift, vecg))
	    {
	      MD_DEBUG31(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
	      MD_DEBUG31(printf("collision will occur at time %.15G\n", vecg[4])); 
	      MD_DEBUG31(printf("[locate_contact] its: %d\n", its));
#ifdef MD_PATCHY_HE
	      if (vecg[4]>t2 || vecg[4]<t1)
#if 0
		|| 
		  (lastbump[i].mol==j && lastbump[i].at==0 && lastbump[j].mol==i && lastbump[j].at==0
		   && fabs(vecg[4] - lastcol[i])<1E-15))
#endif
		  // && !vc_is_pos(i, j, vecg[0], vecg[1], vecg[2], vecg[4]))
		return 0;
	      else
		return 1;

#else
	      if (vecg[4]>t2 || vecg[4]<t1)
#if 0
		|| 
		  (lastbump[i] == j && lastbump[j]==i && fabs(vecg[4] - lastcol[i])<1E-14))
		  // && !vc_is_pos(i, j, vecg[0], vecg[1], vecg[2], vecg[4]))
#endif
		  return 0;
	      else
		{
		  return 1;
		}
#endif
	    }
	  else 
	    {
	      MD_DEBUG31(printf("[locate_contact] can't find contact point!\n"));
	      MD_DEBUG31(printf("troot=%.15G t1=%.15G delt=%.15G t=%.15G troot-t=%.15G\n", troot, t1, delt, t, troot-t-t1));
	      if (d < 0)
		{
		  MD_DEBUG31(printf("t=%.15G d2 < 0 and I did not find contact point, boh...\n",t));
		  //MD_DEBUG10(printf("d1: %.15G d2: %.15G\n", d1, d2));
		  MD_DEBUG31(printf("[locate_contact] its: %d d=%.15G\n", its, d));
		  MD_DEBUG31(printf("dold=%.15G d=%.15G\n", dold, d));
		  return 0;
		  //if (lastbump[i] == j && lastbump[j]==i )
		   // return 0;
		  //exit(-1);
		}
	      else
		{
		  printf("d: %.15G dold: %.15G\n", d, dold); 
		}
		//continue;
	      
	    }
	}
#if 1
      if (d > epsdFastR)
	{
	  if (search_contact_faster(i, j, shift, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2))
	    {
	      MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
	      return 0;
	    }
	  for (kk = 0; kk < 8; kk++)
	    vecgdold[kk] = vecgd[kk];
	  dold = d;
	  its++;
	  //itsS++;
	  continue;
	}
#endif
 
      dold = d;
      for (kk = 0; kk < 8; kk++)
	vecgdold[kk] = vecgd[kk];
      its++;
      itsS++;
    }
  MD_DEBUG(  
  if (foundrc==0)
    printf("%d-%d t=%.12G > t2=%.12G I did not find any contact point!\n", i, j, t, t2);
  );

  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return foundrc;
}

#define EPS 1e-4
#if 1
double estimate_tmin(double t, int na, int nb)
{
  int cellRangeT[2 * NDIM];
  int iX, iY, iZ, jX, jY, jZ, k, n;
  double te, d, tmin=-1, shift[3], r1[3], r2[3], alpha, vecg[8], maxddot;

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
		      d = calcDistNeg(t, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
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
#ifdef EDHE_FLEX
extern void check_these_bonds(int i, int j, double *shift, double t);
#endif
#ifdef MD_EDHEFLEX_WALL
#ifdef MD_ABSORPTION
void calc_grad_and_point_plane_hwsemiperm(int i, double *grad, double *point)
{
  /* nplane = 0 -> -L/2, nplane = 1 -> L/2 */  
  int kk;
  double del=0.0, segno;
  MD_DEBUG33(printf("[PRIMA] del=%.15G grad=%f %f %f point=%f %f %f\n", del, grad[0], grad[1], grad[2], 
		    point[0], point[1], point[2]));

  for (kk=0; kk < 3; kk++)
    {
      grad[kk] = (kk==2)?1:0;
    }
  del = L/2.0-OprogStatus.bufHeight;
  /* NOTA: epsdNL+epsd viene usato come buffer per evitare problemi numerici 
   * nell'update delle NNL. */
  //del -= OprogStatus.epsdNL+OprogStatus.epsd;
  point[2] = del;
  /* N.B. aggiusta il punto in modo che sia allineato con la particella lungo x e y.
     Notare che l'unica cosa che conta è che point[] appartenga al piano in questione.*/
  point[0] = rx[i];
  point[1] = ry[i];
  MD_DEBUG35(printf("[DOPO] del=%.15G grad=%f %f %f point=%f %f %f\n", del, grad[0], grad[1], grad[2], 
		    point[0], point[1], point[2]));
}
#endif
void calc_grad_and_point_plane_hwbump(int i, double *grad, double *point, int nplane)
{
  /* nplane = 0 -> -L/2, nplane = 1 -> L/2 */ 
  int kk;
  double del=0.0, segno;
  MD_DEBUG33(printf("[PRIMA] del=%.15G grad=%f %f %f point=%f %f %f\n", del, grad[0], grad[1], grad[2], 
		    point[0], point[1], point[2]));

  for (kk=0; kk < 3; kk++)
    {
      grad[kk] = (kk==2)?1:0;
    }
  del = L/2.0;
  /* NOTA: epsdNL+epsd viene usato come buffer per evitare problemi numerici 
   * nell'update delle NNL. */
  //del -= OprogStatus.epsdNL+OprogStatus.epsd;
  for (kk=0; kk < 3; kk++)
    {
      if (nplane % 2 == 0)
	segno = -1.0;
      else
	segno = 1.0;
      grad[kk] *= segno;
      /* rB[] (i.e. nebrTab[i].r[]) è un punto appartenente al piano */
      point[kk] = (kk==2)?del*grad[kk]:0; 
    }
  /* N.B. aggiusta il punto in modo che sia allineato con la particella lungo x e y.
     Notare che l'unica cosa che conta è che point[] appartenga al piano in questione.*/
  point[0] = rx[i];
  point[1] = ry[i];
  MD_DEBUG35(printf("[DOPO] del=%.15G grad=%f %f %f point=%f %f %f\n", del, grad[0], grad[1], grad[2], 
		    point[0], point[1], point[2]));
}
int globalHW = 0;
void bumpHW(int i, int nplane, double rCx, double rCy, double rCz, double *W)
{
  double factor, invmi;
  double wrx, wry, wrz, rACn[3] ;
  double rAC[3], vCA[3], vc;
  double norm[3];
  double denom;
  double invaSq, invbSq, invcSq;
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, taus, DTxx, DTyy, DTzz ;
  double rxij, ryij, rzij, Dr;
#endif
#ifndef MD_ASYM_ITENS
  double factorinvIa;
#endif
  int k;
#ifdef MD_ASYM_ITENS
  int k1,k2;
  double rnI[3];
  double Mvec[3], omega[3];
#endif
  int typei;
  int a, b;

  MD_DEBUG38(printf("Collision with wall evIdA: %d vz: %.15f rz:%.15G evIdB=%d\n", i, vz[i], rz[i], evIdB-ATOM_LIMIT));
  if (is_sphere(i))
    {
      update_MSDrot(i);
      vz[i] = -vz[i];
#if 0
#ifdef MD_ASYM_ITENS
      calc_angmom(i, Ia);
      upd_refsysM(i);
#endif
#endif
      return;
    }
  rAC[0] = rx[i] - rCx;
  rAC[1] = ry[i] - rCy;
  rAC[2] = rz[i] - rCz;
  
  MD_DEBUG35(printf("i=%d nplane=%d rAC=%f %f %fi rC=%f %f %f\n", i, nplane, rAC[0], rAC[1], rAC[2], rCx, rCy, rCz));
  MD_DEBUG35(printf("r(%d)=%f %f %f time=%.15G\n", i, rx[i], ry[i], rz[i], Oparams.time));
  typei = typeOfPart[i];
  invaSq = 1/Sqr(typesArr[typei].sax[0]);
  invbSq = 1/Sqr(typesArr[typei].sax[1]);
  invcSq = 1/Sqr(typesArr[typei].sax[2]);
  tRDiagR(i, Xa, invaSq, invbSq, invcSq, R[i]);
#ifdef MD_ASYM_ITENS
  tRDiagR(i, Ia, typesArr[typei].I[0], typesArr[typei].I[1], typesArr[typei].I[2], R[i]);
#else
  Ia = Oparams.I[na];
#endif
#ifdef MD_ASYM_ITENS
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Iatmp[k1][k2] = Ia[k1][k2];
      } 
  InvMatrix(Iatmp, invIa, 3);
  Mvec[0] = Mx[i];
  Mvec[1] = My[i];
  Mvec[2] = Mz[i];
  for (k1 = 0; k1 < 3; k1++)
    {
      omega[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	omega[k1] += invIa[k1][k2]*Mvec[k2]; 
    }
  wx[i] = omega[0];
  wy[i] = omega[1];
  wz[i] = omega[2];
#else
  invIa = 1/Ia;
  MD_DEBUG(printf("Ia=%f Ib=%f\n", Ia, Ib));
#endif
  for (k=0; k<3; k++)
    norm[k] = (k==2)?(2.0*nplane-1.0):0;

  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[2] = vz[i] + wrz;
  invmi = 1.0/typesArr[typei].m;
  denom = invmi;
  vc = vCA[2]*norm[2];
  rACn[0] = rAC[1]*norm[2];
  rACn[1] = -rAC[0]*norm[2];

#ifdef MD_ASYM_ITENS 
  for (a=0; a < 3; a++)
    {
      rnI[a] = 0;
      for (b = 0; b < 2; b++)
	{
	  rnI[a] += invIa[a][b]*rACn[b]; 
	}
    }
  for (a = 0; a < 2; a++)
    denom += rnI[a]*rACn[a];
#else
  for (a = 0; a < 2; a++)
    denom += invIa*Sqr(rACn[a]);
#endif
  factor = 2.0*vc/denom;
  vz[i] = vz[i] - invmi*norm[2]*factor;
  update_MSDrot(i);
  if (rAC[0]!=0.0 || rAC[1]!=0.0)
    { 
      MD_DEBUG35(printf("QUI factor=%.15G wrz=%.15G\n", factor,wrz));
       //rACn[2] = 0.0;
      //vectProd(rAC[0], rAC[1], rAC[2], norm[0], norm[1], norm[2], &rACn[0], &rACn[1], &rACn[2]);
#ifdef MD_ASYM_ITENS
      /* rACn[2]=0 per quello a < 2*/
      for (a=0; a < 2; a++)
	{
	  wx[i] += factor*invIa[0][a]*rACn[a];
	  wy[i] += factor*invIa[1][a]*rACn[a];
	  wz[i] += factor*invIa[2][a]*rACn[a];
	}
#else
      factorinvIa = factor*invIa;
      wx[i] += factorinvIa*rACn[0];
      wy[i] += factorinvIa*rACn[1];
      /* N.B. rACn[2] = 0 */
      //wz[i] += factorinvIa*rACn[2];
#endif
#ifdef MD_ASYM_ITENS
      calc_angmom(i, Ia);
      upd_refsysM(i);
#endif
    }
}
void PredictEvent (int na, int nb); 

void ProcessWallColl(void)
{
  int k;
  MD_DEBUG35(printf("PRIr(%d)=%f %f %f time=%.15G\n", evIdA, rx[evIdA], ry[evIdA], rz[evIdA], Oparams.time));
  UpdateAtom(evIdA);
  MD_DEBUG35(printf("DOPr(%d)=%f %f %f time=%.15G\n", evIdA, rx[evIdA], ry[evIdA], rz[evIdA], Oparams.time));
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  MD_DEBUG36(calc_energy("prima"));
  MD_DEBUG38(printf("[BUMP WALLCOLL] t=%.15G i=%d at=%d j=%d at=%d collCode=%d\n", 
		    Oparams.time,evIdA,evIdC, evIdB, evIdD, evIdE)); 

  MD_DEBUG38( printf("i=%d pos=%f %f %f vel=%f %f %f\n", evIdA, rx[evIdA], ry[evIdA], rz[evIdA], vx[evIdA], vy[evIdA], vz[evIdA]));
  bumpHW(evIdA, evIdB-ATOM_LIMIT, rxC, ryC, rzC, &W);
  MD_DEBUG36(calc_energy("dopo"));
  MD_DEBUG(store_bump(evIdA, evIdB));
  //ENDSIM=1;
  /*printf("qui time: %.15f\n", Oparams.time);*/
  OprogStatus.lastcolltime[evIdA] = lastcol[evIdA] = Oparams.time;
#ifdef MD_CALC_DPP
  store_last_u(evIdA);
#endif
  lastbump[evIdA].mol=evIdB-ATOM_LIMIT;
  lastbump[evIdA].at = evIdC;
  //printf("lastbump[%d].at=%d lastbump[%d].at=%d\n", evIdA, lastbump[evIdA].at, evIdB, lastbump[evIdB].at);
  lastbump[evIdA].type = evIdE;
#ifdef MD_HSVISCO
  OprogStatus.lastcoll = Oparams.time;
#endif
  if (OprogStatus.useNNL)
    {
      /* ricalcola i tempi di collisione con la NL */
      updrebuildNNL(evIdA);
      PredictEventNNL(evIdA, -1);
    }
  else
    {
      PredictEvent(evIdA, -1);
    }
}
extern int locate_contact_neigh_plane(int i, double vecg[5], int nplane, double tsup);
extern double gradplane[3], rB[3];

extern void calc_grad_and_point_plane(int i, double *grad, double *point, int nplane);
extern int locate_contact_neigh_plane_HS_one(int i, double *evtime, double t2);

int locateHardWall(int na, int nplane, double tsup, double vecg[5], int ghw)
{
  double delt, tone;
  if ((is_infinite_Itens(na) && is_infinite_mass(na)) ||
      (is_infinite_mass(na) && is_a_sphere_NNL[na]))
    {
      return 0;
    }
  globalHW = ghw;
#ifdef EDHE_FLEX
  /* ottimizzazione nel caso di sfere che urtano il muro duro */
  if (is_sphere(na))
    {
      calc_grad_and_point_plane(na, gradplane, rB, nplane);
      //printf("type=%d ghw=%d na=%d (%.15G %.15G %.15G)\n", typeOfPart[na], ghw, na, rx[na], ry[na], rz[na]);

      //if (ghw==2)
	//tsup = timbig;
      if(!locate_contact_neigh_plane_HS_one(na, &vecg[4], tsup))
	{
	  globalHW = 0;
	  //printf("NO time one=%.15G\n", vecg[4]);
	  return 0;
	}

      globalHW=0;
      delt = vecg[4]-Oparams.time;
      vecg[0] = rx[na]+vx[na]*delt;
      vecg[1] = ry[na]+vy[na]*delt;
      if (nplane==1)
	vecg[2] = rz[na]+vz[na]*delt+typesArr[typeOfPart[na]].sax[0];
      else
	vecg[2] = rz[na]+vz[na]*delt-typesArr[typeOfPart[na]].sax[0];
 
      if (vecg[4] < tsup)
	return 1;
      return 0;
    }
#endif

  MD_DEBUG33(printf("inCell=%d pos of %d=%f %f %f\n", inCell[2][na], na, rx[na], ry[na], rz[na]));
  if (!locate_contact_neigh_plane(na, vecg, nplane, tsup))
    {
      //printf("NO locate neigh full %.15G\n", vecg[4]);
      globalHW = 0;
      return 0;
    }
  globalHW = 0;
  if (vecg[4] < tsup)
    {
#if 0
      printf("OK locate neigh full %.15G\n", vecg[4]);
      if (fabs(tone-vecg[4]) > 1E-10)
	exit(-1);
#endif
      return 1;
    }
  return 0;
}
#endif
#ifdef EDHE_FLEX
extern int may_interact_all(int i, int j);
extern int *is_a_sphere_NNL;
#endif
void PredictEvent (int na, int nb) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq, dr[NDIM], dv[NDIM], shift[NDIM], tm[NDIM], 
	 b, d, t, tInt, vv, distSq, t1, t2;
  int overlap;
#ifdef MD_ABSORPTION
  int hwcell;
#endif
#ifdef MD_PATCHY_HE
  int ac, bc, collCode, collCodeOld, acHC, bcHC;
  double evtime, evtimeHC;
  int nplane=-1;
#endif
#ifdef EDHE_FLEX
  double cctime;
#endif
  double vecg[5];
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cells[NDIM];*/
#ifdef MD_GRAVITY
  double Lzx, h1, h2, sig, hh1;
#endif
  int cellRangeT[2 * NDIM], signDir[NDIM]={0,0,0}, evCode,
  iX, iY, iZ, jX, jY, jZ, k, n;

  MD_DEBUG38(printf("PredictEvent: %d,%d\n", na, nb));
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
#if 0
  if (evIdB < ATOM_LIMIT)
    {
  if (rx[na] < inCell[0][na]*L/cellsx - L2 || rx[na] > (inCell[0][na]+1)*L/cellsx-L2)
    {
      printf("ERRORE nelle CELLE LUNGO X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("na=%d (%f,%f,%f) cells(%d,%d,%d)\n",
	     na, rx[na], ry[na], rz[na], inCell[0][na], inCell[1][na], inCell[2][na]);
      exit(-1);
    }
   else if (ry[na] < inCell[1][na]*L/cellsy -L2 || ry[na] > (inCell[1][na]+1)*L/cellsy-L2)
     {
       printf("ERRORE nelle CELLE LUNGO Y!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
       printf("na=%d (%f,%f,%f) cells(%d,%d,%d)\n",
 	      na, rx[na], ry[na], rz[na], inCell[0][na], inCell[1][na], inCell[2][na]);
       exit(-1);
     }
   else if (rz[na] < inCell[2][na]*(Lz+OprogStatus.extraLz)/cellsz - Lz2 || rz[na] > (inCell[2][na]+1)*(Lz+OprogStatus.extraLz)/cellsz -Lz2) 
     {
       printf("ERRORE nelle CELLE LUNGO Z!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
       printf("na=%d (%f,%f,%f) cells(%d,%d,%d)\n",
	      na, rx[na], ry[na], rz[na], inCell[0][na], inCell[1][na], inCell[2][na]);
       printf("cellz0: %f cellz1: %f\n", 
	      inCell[2][na]*(Lz+OprogStatus.extraLz)/cellsz -Lz2,
	     (inCell[2][na]+1)*(Lz+OprogStatus.extraLz)/cellsz -Lz2 );
       printf("evIdA: %d evIdB: %d\n", evIdA, evIdB);
       printf("vz[na]:%f\n", vz[na]);
       exit(-1);
     } } 
#endif
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
#ifdef MD_EDHEFLEX_WALL
      if (OprogStatus.hardwall && ((signDir[2]==0 && inCell[2][na]==cellsz-1) || (signDir[2]==1 && inCell[2][na]==0)))
	tm[2] = timbig;
      else
	tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L /
		 cellsz - rz[na] - L2) / vz[na];
#else
      tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L /
	       cellsz - rz[na] - L2) / vz[na];
#endif
    } 
  else 
    tm[2] = timbig;
#endif
  /* end forcefield[k] != 0*/
  
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
  if (tm[1] <= tm[2]) {
    if (tm[0] <= tm[1]) k = 0;
    else k = 1;
  } else {
    if (tm[0] <= tm[2]) k = 0;
    else k = 2;
  }
  /* Se un errore numerico fa si che tm[k] < 0 allora lo poniamo uguale a 0
   * (ved. articolo Lubachevsky) */
#if 1
  if (tm[k]<0)
    {
#if 1
      printf("tm[%d]<0=%.15G step %lld na=%d\n", k, tm[k], (long long int)Oparams.curStep, na);
#ifdef MD_GRAVITY
      printf("rz:%f diff:%f\n", rz[na], rz[na]+Lz2);
      printf("h1:%f hh1:%f vz:%f cellz:%d\n", h1, hh1, vz[na], inCell[2][na]);
#endif
      printf("Cells(%d,%d,%d)\n", inCell[0][na], inCell[1][na], inCell[2][na]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
      /*exit(-1);*/
      /*tm[k] = 0.0;*/
#endif
      tm[k] = 0.0;
    }
#endif
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 
#ifdef MD_GRAVITY
  if (k == 2 && inCell[2][na] == 0 && signDir[2] == 1) 
    {
      evCode = 4;/* sarebbe 2*k con k=2 (z) e per me vuol dire urto con parete in basso
  		    che è anche l'unica nel presente caso */
      MD_DEBUG2(printf("wall!!! (evIdA: %d)\n", na));
    }
#endif
  MD_DEBUG15(printf("schedule event [WallCrossing](%d,%d) tm[%d]: %.8G\n", 
		    na, ATOM_LIMIT+evCode, k, tm[k]));
#ifndef MD_EDHEFLEX_WALL
  ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]);
#else
  cctime = Oparams.time+tm[k];
#endif
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#if (defined(MD_GRAVITY) || defined(MD_EDHEFLEX_WALL))
  /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (OprogStatus.hardwall)
    {
#if defined(MD_ABSORPTION) 
      if (OprogStatus.bufHeight > 0.0)
	{
	  if (vz[na] != 0.0)
	    {
	      hwcell = (L-OprogStatus.bufHeight)*cellsz/L;
#if 1
	      if (hwcell-inCell[2][na] < 2)
		{
		  /* the semi-permeable plane is just one (nplane=0) */
		  if (locateHardWall(na, 0, cctime, vecg, 2))
		    {
		      rxC = vecg[0];
		      ryC = vecg[1];
		      rzC = vecg[2];
		      MD_DEBUG38(printf("SEMIPERM Located Contact with WALL rC=%f %f %f time=%.15G i=%d\n", rxC, ryC, rzC, vecg[4], na));
		      MD_DEBUG38(printf("r=%f %f %f\n", rx[na], ry[na], rz[na]));
		      ScheduleEventBarr (na, ATOM_LIMIT+50, 0, 0, MD_WALL, vecg[4]);
		    }
		}
#endif
	    }
	}
#endif

      if (inCell[2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (inCell[2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
      if (inCell[2][na] == 0)
	nplane = 0;
      else if (inCell[2][na] == cellsz-1)
	nplane = 1;
      if (nplane!=-1 && locateHardWall(na, nplane, cctime, vecg, 1))
	{
	  rxC = vecg[0];
	  ryC = vecg[1];
	  rzC = vecg[2];
	  MD_DEBUG38(printf("scheduled collision with wall na=%d nplane=%d evtime=%.15G\n", na, nplane, vecg[4]));
	  ScheduleEventBarr (na, ATOM_LIMIT + nplane, 0, 0, MD_WALL, vecg[4]);
	}
      else
	{
	  MD_DEBUG38(printf("scheduled cell-crossing with wall na=%d evtime=%.15G\n", na, cctime));
	  ScheduleEvent (na, ATOM_LIMIT + evCode, cctime);
	}
    }
  else
    ScheduleEvent (na, ATOM_LIMIT + evCode, cctime);
#endif

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
#ifdef EDHE_FLEX
		      if (!may_interact_all(na, n))
			continue;
#endif
	
		      /* maxax[...] è il diametro dei centroidi dei due tipi
		       * di ellissoidi */
		      if (OprogStatus.targetPhi > 0)
			{
			  sigSq = Sqr(max_ax(na)+max_ax(n)+OprogStatus.epsd);
			}
		      else
		       {
#if defined(MD_POLYDISP) || defined(EDHE_FLEX)
			 sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#else
			 if (na < parnumA && n < parnumA)
			   sigSq = Sqr(maxax[na]+OprogStatus.epsd);
			 else if (na >= parnumA && n >= parnumA)
			   sigSq = Sqr(maxax[na]+OprogStatus.epsd);
			 else
			   sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#endif
		       }
		      MD_DEBUG20(printf("sigSq: %f maxax[%d]:%.15G maxax[%d]:%.15G\n", sigSq, na, maxax[na],
					n, maxax[n]));
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
		      distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
		      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
	    	      d = Sqr (b) - vv * (distSq - sigSq);
		
		      if (d < 0 || (b > 0.0 && distSq > sigSq)) 
			{
			  /* i centroidi non collidono per cui non ci può essere
			   * nessun urto sotto tali condizioni */
			  continue;
			}
		      MD_DEBUG20(printf("PREDICTING na=%d n=%d\n", na , n));
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
			  overlap = 0;
			  MD_DEBUG20(printf("qui..boh sig:%.15G dist=%.15G\n", sqrt(sigSq), sqrt(distSq)));
			}
		      else 
			{
			  MD_DEBUG(printf("Centroids overlap!\n"));
			  t2 = t = (sqrt (d) - b) / vv;
			  t1 = 0.0; 
			  overlap = 1;
			  MD_DEBUG(printf("altro d=%f t=%.15f\n", d, (-sqrt (d) - b) / vv));
			  MD_DEBUG(printf("vv=%f dv[0]:%f\n", vv, dv[0]));
			}
		      MD_DEBUG(printf("t=%f curtime: %f b=%f d=%f\n", t, Oparams.time, b ,d));
		      MD_DEBUG(printf("dr=(%f,%f,%f) sigSq: %f", dr[0], dr[1], dr[2], sigSq));
		      //t += Oparams.time; 
		      t2 += Oparams.time;
		      t1 += Oparams.time;
#if 0 
		      calcDist(Oparams.time, na, n, shift, r1, r2);
		      continue;
		      exit(-1);
#endif
#ifdef MD_PATCHY_HE
		      evtime = t2;
		      collCode = MD_EVENT_NONE;
		      rxC = ryC = rzC = 0.0;
		      MD_DEBUG20(printf("time=%.15G t1=%.15G t2=%.15G maxax[%d]:%.15G maxax[%d]:%.15G\n", 
					Oparams.time, t1, t2, na, maxax[na], n, maxax[n]));
		      MD_DEBUG31(printf("t1=%.15G t2=%.15G\n", t1, t2));
		      collCodeOld = collCode;
		      evtimeHC = evtime;
		      acHC = ac = 0;
		      bcHC = bc = 0;
#ifdef EDHE_FLEX
		      if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			{
			  collCode = MD_EVENT_NONE;
			}
		      MD_DEBUG20(if (collCode!=MD_EVENT_NONE) printf("na=%d ac=%d n=%d bc=%d\n", na, ac, n, bc));
		      MD_DEBUG20(if (collCode!=MD_EVENT_NONE) check_these_bonds(na, n, shift, evtime));
#else
		      if (OprogStatus.targetPhi <=0 && ((na < Oparams.parnumA && n >= Oparams.parnumA)|| 
							(na >= Oparams.parnumA && n < Oparams.parnumA)))
			{
			  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			    {
			      collCode = MD_EVENT_NONE;
			    }
			}
#endif
		      if (collCode!=MD_EVENT_NONE)
			t2 = evtime+1E-7;
		      if (locate_contact(na, n, shift, t1, t2, vecg))
			{
			  if (collCode == MD_EVENT_NONE || (collCode!=MD_EVENT_NONE && vecg[4] <= evtime))
			    {
			      collCode = MD_CORE_BARRIER;
			      ac = bc = 0;
			      evtime = vecg[4];
			      rxC = vecg[0];
			      ryC = vecg[1];
			      rzC = vecg[2];
			    }
			}
		      else
			{
			  if (collCode == MD_EVENT_NONE)
			    continue;
			}

		      t = evtime;
		      MD_DEBUG31(printf("(%d,%d)-(%d,%d) troot=%.15G\n", na, ac, n, bc, evtime));
		      MD_DEBUG31(printf("evtimeHC=%.15G",evtimeHC));
		      MD_DEBUG31(printf("dist(t=%.15G,%d-%d):%.15G dist(t=%.15G)=%.15G\n", evtimeHC, na, n, 
					calc_dist_ij(na, n, evtimeHC), evtime, calc_dist_ij(na, n, evtime)));
#else
		      if (!locate_contact(na, n, shift, t1, t2, vecg))
		      	continue;
		      rxC = vecg[0];
		      ryC = vecg[1];
		      rzC = vecg[2];
		      MD_DEBUG(printf("A x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)-B x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)",
				      rx[na], ry[na], rz[na], vx[na], vy[na], vz[na],
				      rx[n], ry[n], rz[n], vx[n], vy[n], vz[n]));
		      t = vecg[4];

#endif
#ifdef MD_PATCHY_HE
		      if (t < Oparams.time)
			{
#if 1
			  printf("time:%.15f tInt:%.15f\n", Oparams.time,
				 tInt);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
	     					      Sqr(dr[2]))-1.0 );
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = Oparams.time;
			}

#else
		      if (t < 0)
			{
#if 1
			  printf("time:%.15f tInt:%.15f\n", Oparams.time,
				 tInt);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
	     					      Sqr(dr[2]))-1.0 );
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = 0;
			}
#endif
		      /* il tempo restituito da newt() è già un tempo assoluto */
#ifdef MD_PATCHY_HE
		      MD_DEBUG38(printf("time: %f Adding collision (%d,%d)-(%d,%d)\n", t, na, ac, n, bc));
		      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
#else
		      ScheduleEvent (na, n, t);
#endif
		      MD_DEBUG20(printf("schedule event [collision](%d,%d)-(%d,%d) collCode=%d\n", na, ac, n, bc, collCode));
		    }
		} 
	    }
	}
    }
}

#ifdef MD_GRAVITY
void ProcessCollWall(void)
{
  /* Dissipation */
  MD_DEBUG2(printf("Collision with wall evIdA: %d vz: %.15f\n", evIdA, vz[evIdA]));
#if 0
  printf("timeNow-lastcol[%d]:%.15f\n", evIdA, Oparams.time - lastcol[evIdA]);
  printf("walldiss:%.15f vz:%.10f\n", Oparams.wallDiss, vz[evIdA]);
  printf("timeNow:%.10f lastcol:%.10f\n", Oparams.time, lastcol[evIdA]);
#endif
  if (Oparams.time - lastcol[evIdA] < OprogStatus.tc)
    {
      vz[evIdA] = -vz[evIdA];
    }
  else
    {
      vz[evIdA] = -Oparams.wallDiss*vz[evIdA];
    }	 
  lastcol[evIdA] = Oparams.time;
}
#endif
void calc_energynew(char *msg)
{
  int i, k1;
  double wt[3];
#ifdef EDHE_FLEX
  int typei;
#endif
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
#ifdef EDHE_FLEX
	      typei = typeOfPart[i];
      	      K += wt[k1]*typesArr[typei].I[k1]*wt[k1];
#else
      	      K += wt[k1]*Oparams.I[0][k1]*wt[k1];
#endif
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
#ifdef EDHE_FLEX
	      typei = typeOfPart[i];
      	      K += wt[k1]*typesArr[typei].I[k1]*wt[k1];
#else
	      K += wt[k1]*Oparams.I[1][k1]*wt[k1];
#endif
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
  printf("[%s] Kinetic Energy: %f\n", msg, K);
}
void calc_energy_i(char *msg, int i)
{
  int k1;
  double wt[3];
#ifdef EDHE_FLEX
  int typei;
#endif
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  int k2;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#endif
  K = 0;
  if (i<Oparams.parnumA)
    {
      /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
#ifdef EDHE_FLEX
      typei = typeOfPart[i];
      RDiagtR(i, Ia, typesArr[typei].I[0], typesArr[typei].I[1], typesArr[typei].I[2], R[i]);
#else
      RDiagtR(i, Ia, Oparams.I[0][0], Oparams.I[0][1], Oparams.I[0][2], R[i]);
#endif
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
#ifdef EDHE_FLEX
      typei = typeOfPart[i];
      RDiagtR(i, Ia, typesArr[typei].I[0], typesArr[typei].I[1], typesArr[typei].I[2], R[i]);
#else
      RDiagtR(i, Ib, Oparams.I[1][0], Oparams.I[1][1], Oparams.I[1][2], R[i]);
#endif
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
  printf("[%s] Kinetic Energy of %d: %f\n", msg, i, K);
  
}
#ifdef MD_ASYM_ITENS
void calc_omega(int i)
{
  double Mvec[3], omega[3];
  int k1, k2, na;
#ifdef EDHE_FLEX
  int typei;
#endif
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  typei = typeOfPart[i];	
  if (is_infinite_Itens(i))
    {
      wx[i] = wy[i] = wz[i] = 0.0;
      return;     
    }
  tRDiagR(i, Ia, typesArr[typei].I[0], typesArr[typei].I[1], typesArr[typei].I[2], R[i]);
#else
  tRDiagR(i, Ia, Oparams.I[na][0], Oparams.I[na][1], Oparams.I[na][2], R[i]);
#endif
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Iatmp[k1][k2] = Ia[k1][k2];
      }
  InvMatrix(Iatmp, invIa, 3);
  Mvec[0] = Mx[i];
  Mvec[1] = My[i];
  Mvec[2] = Mz[i];
  for (k1 = 0; k1 < 3; k1++)
    {
      omega[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	omega[k1] += invIa[k1][k2]*Mvec[k2]; 
    }
  wx[i] = omega[0];
  wy[i] = omega[1];
  wz[i] = omega[2];
  //printf("i=%d M=%.15G %.15G %.15G w[%d]=%.15G %.15G %.15G\n", i, Mx[i], My[i], Mz[i], i, wx[i], wy[i], wz[i]);
}
#endif
#ifdef MD_PATCHY_HE
extern double calcpotene(void);
#endif
double Krot, Ktra;
#ifdef EDHE_FLEX
void calc_energy_filtered(int filter)
{
  int i, k1;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double wtp[3];
  int k2;
  //double **Ia, **Ib;
  //Ia = matrix(3,3); 
  //Ib = matrix(3,3);
#endif
  K = Ktra = Krot = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (filter!=0 && typesArr[typeOfPart[i]].brownian!=filter)
	continue;
      /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
      //RDiagtR(i, Ia, Oparams.I[0][0], Oparams.I[0][1], Oparams.I[0][2], R[i]);
#endif
      Ktra += typesArr[typeOfPart[i]].m*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
#ifdef MD_ASYM_ITENS
      calc_omega(i);
#endif
      wt[0] = wx[i];
      wt[1] = wy[i];
      wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
      for (k1=0; k1 < 3; k1++)
	{
	  wtp[k1] = 0.0;
	  for (k2=0; k2 < 3; k2++)
	    wtp[k1] += R[i][k1][k2]*wt[k2];
	}
      //printf("calcnorm wt: %.15G wtp:%.15G\n", calc_norm(wt), calc_norm(wtp));
      for (k1=0; k1 < 3; k1++)
	{
	  Krot += Sqr(wtp[k1])*typesArr[typeOfPart[i]].I[k1];
	  //printf("I[%d][%d]=%.15G wt[%d]:%.15G wtp[%d]:%.15G\n", 0, k1, Oparams.I[0][k1],
	  //     k1, wt[k1], k1, wtp[k1]);
	}
#else
      for (k1=0; k1 < 3; k1++)
	Krot += Sqr(wt[k1])*typesArr[typeOfPart[i]].I[k1];
#endif
    }
  Ktra *= 0.5;
  Krot *= 0.5;
  K = Ktra + Krot;
}
#endif

void calc_energy(char *msg)
{
  int i, k1;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double wtp[3];
  int k2;
  //double **Ia, **Ib;
  //Ia = matrix(3,3); 
  //Ib = matrix(3,3);
#endif

  K = Ktra = Krot = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (i<Oparams.parnumA)
	{
	  /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
	  //RDiagtR(i, Ia, Oparams.I[0][0], Oparams.I[0][1], Oparams.I[0][2], R[i]);
#endif
#ifdef EDHE_FLEX
	  Ktra += typesArr[typeOfPart[i]].m*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
#else
	  Ktra += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
#endif
#ifdef MD_ASYM_ITENS
	  calc_omega(i);
#endif
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
	  for (k1=0; k1 < 3; k1++)
	    {
	      wtp[k1] = 0.0;
	      for (k2=0; k2 < 3; k2++)
		wtp[k1] += R[i][k1][k2]*wt[k2];
	    }
	  //printf("calcnorm wt: %.15G wtp:%.15G\n", calc_norm(wt), calc_norm(wtp));
	  for (k1=0; k1 < 3; k1++)
	    {
#ifdef EDHE_FLEX
	      Krot += Sqr(wtp[k1])*typesArr[typeOfPart[i]].I[k1];
#else
	      Krot += Sqr(wtp[k1])*Oparams.I[0][k1];
#endif
	      //printf("I[%d][%d]=%.15G wt[%d]:%.15G wtp[%d]:%.15G\n", 0, k1, Oparams.I[0][k1],
		//     k1, wt[k1], k1, wtp[k1]);
	    }
#else
	  for (k1=0; k1 < 3; k1++)
	    {
#ifdef EDHE_FLEX
	      Krot += Sqr(wt[k1])*typesArr[typeOfPart[i]].I[k1];
#else
	      Krot += Sqr(wt[k1])*Oparams.I[0];
#endif
	    }
#endif
	}
      else
	{
#ifdef MD_ASYM_ITENS
	  //RDiagtR(i, Ib, Oparams.I[1][0], Oparams.I[1][1], Oparams.I[1][2], R[i]);
#endif
	  Ktra += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
#ifdef MD_ASYM_ITENS
	  calc_omega(i);
#endif
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
	  for (k1=0; k1 < 3; k1++)
	    {
	      wtp[k1] = 0.0;
	      for (k2=0; k2 < 3; k2++)
		wtp[k1] += R[i][k1][k2]*wt[k2];
	    }
	  for (k1=0; k1 < 3; k1++)
	    Krot += Sqr(wtp[k1])*Oparams.I[1][k1];
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
  //free_matrix(Ia,3);
  //free_matrix(Ib,3);
#endif
  if (msg)
    {
      printf("[%s] Kinetic Energy: %.15G\n", msg, K);
#ifdef MD_PATCHY_HE
      printf("tot energy: %.15G\n", K+calcpotene());
#endif
    }
}
void store_bump_neigh(int i, double *r1, double *r2)
{
  char fileop2[512], fileop[512];
  int ii;
  FILE *bf;
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G @ %.15G %.15G %.15G C[%s]\n";
#ifdef EDHE_FLEX
  int kk;
  double axi[3];
#endif
#ifdef MD_BIG_DT
  sprintf(fileop2 ,"StoreBumpNeigh-%d-t%.8f", i, Oparams.time + OprogStatus.refTime);
#else
  sprintf(fileop2 ,"StoreBumpNeigh-%d-t%.8f", i, Oparams.time);
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
  for (ii = 0; ii < Oparams.parnum; ii++)
    {
      if (ii==i)
	{
#ifdef EDHE_FLEX
	  for (kk=0; kk < 3; kk++)
	    axi[kk] = typesArr[typeOfPart[i]].sax[kk];
	  fprintf(bf, tipodat2,rx[ii], ry[ii], rz[ii], uxx[ii], uxy[ii], uxz[ii], uyx[ii], uyy[ii], 
	  	  uyz[ii], uzx[ii], uzy[ii], uzz[ii], axi[0], axi[1], axi[2], "red");

#else
	  fprintf(bf, tipodat2,rx[ii], ry[ii], rz[ii], uxx[ii], uxy[ii], uxz[ii], uyx[ii], uyy[ii], 
	  	  uyz[ii], uzx[ii], uzy[ii], uzz[ii], axa[i], axb[i], axc[i], "red");
#endif
	  fprintf(bf, tipodat2,nebrTab[i].r[0], nebrTab[i].r[1], nebrTab[i].r[2], nebrTab[i].R[0][0], nebrTab[i].R[0][1], 
		  nebrTab[i].R[0][2], nebrTab[i].R[1][0], nebrTab[i].R[1][1], nebrTab[i].R[1][2], 
	  	  nebrTab[i].R[2][0], nebrTab[i].R[2][1], nebrTab[i].R[2][2], 
		  nebrTab[i].axa, nebrTab[i].axb, 
		  nebrTab[i].axc, "green");


	}
    }
  //writeAllCor(bf);
  fprintf(bf,"%.15f %.15f %.15f @ 0.2 C[blue]\n", r1[0], r1[1], r1[2]);
  fprintf(bf,"%.15f %.15f %.15f @ 0.2 C[blue]\n", r2[0], r2[1], r2[2]);
  fclose(bf);

}
#ifdef MD_PATCHY_HE
extern void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3]);
#endif
void store_bump(int i, int j)
{
  char fileop2[512], fileop[512];
  FILE *bf;
  int na;
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
#ifdef MD_PATCHY_HE
  int nn;
  double ratA[NA][3], ratB[NA][3];
#endif
  double Drx, Dry, Drz, RCMx, RCMy, RCMz;
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G @ %.15G %.15G %.15G C[%s]\n";
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
  fprintf(bf, ".Vol: %f\n", Oparams.rcut*Oparams.rcut*Oparams.rcut);
  MD_DEBUG(printf("[Store bump]: %.15G\n", Oparams.time));
  Drx = L*rint((rx[i]-rx[j])/L);
  Dry = L*rint((ry[i]-ry[j])/L);
  Drz = L*rint((rz[i]-rz[j])/L);

  RCMx = (rx[i]+rx[j]+Drx)*0.5;
  RCMy = (ry[i]+ry[j]+Dry)*0.5;
  RCMz = (rz[i]+rz[j]+Drz)*0.5;
#ifdef EDHE_FLEX
  for (kk=0; kk < 3; kk++)
    {
      axi[kk] = typesArr[typeOfPart[i]].sax[kk];
      axj[kk] = typesArr[typeOfPart[j]].sax[kk];
    }
  fprintf(bf, tipodat2,rx[i]-RCMx, ry[i]-RCMy, rz[i]-RCMz, uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	  uyz[i], uzx[i], uzy[i], uzz[i], axi[0], axi[1], axi[2], "red");
  fprintf(bf, tipodat2,rx[j]+Drx-RCMx, ry[j]+Dry-RCMy, rz[j]+Drz-RCMz, uxx[j], uxy[j], uxz[j], uyx[j], uyy[j], 
	  uyz[j], uzx[j], uzy[j], uzz[j], axj[0], axj[1], axj[2], "blue");

#elif defined(MD_POLYDISP)
  fprintf(bf, tipodat2,rx[i]-RCMx, ry[i]-RCMy, rz[i]-RCMz, uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	  uyz[i], uzx[i], uzy[i], uzz[i], axaP[i], axbP[i], axcP[i], "red");
  fprintf(bf, tipodat2,rx[j]+Drx-RCMx, ry[j]+Dry-RCMy, rz[j]+Drz-RCMz, uxx[j], uxy[j], uxz[j], uyx[j], uyy[j], 
	  uyz[j], uzx[j], uzy[j], uzz[j], axaP[j], axbP[j], axcP[j], "blue");
#else
  na = (i < Oparams.parnumA)?0:1;
  fprintf(bf, tipodat2,rx[i]-RCMx, ry[i]-RCMy, rz[i]-RCMz, uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	  uyz[i], uzx[i], uzy[i], uzz[i], Oparams.a[na], Oparams.b[na], Oparams.c[na], "red");
  na = (j < Oparams.parnumA)?0:1;
  fprintf(bf, tipodat2,rx[j]+Drx-RCMx, ry[j]+Dry-RCMy, rz[j]+Drz-RCMz, uxx[j], uxy[j], uxz[j], uyx[j], uyy[j], 
	  uyz[j], uzx[j], uzy[j], uzz[j], Oparams.a[na], Oparams.b[na], Oparams.c[na], "blue");
#endif
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
  rB[0] = rx[j]-RCMx+Drx;
  rB[1] = ry[j]-RCMy+Dry;
  rB[2] = rz[j]-RCMz+Drz;
  BuildAtomPos(j, rB, R[j], ratB);
#ifdef EDHE_FLEX
  for (nn = 1; nn < typesArr[typeOfPart[j]].nspots+1; nn++)
    fprintf(bf,"%.15f %.15f %.15f @ %.15G C[orange]\n", 
	    ratA[nn][0], ratA[nn][1], ratA[nn][2], typesArr[typeOfPart[j]].spots[nn-1].sigma*0.5);
#else
  for (nn = 1; nn < ((j < Oparams.parnumA)?MD_STSPOTS_A+1:MD_STSPOTS_B+1); nn++)
    fprintf(bf,"%.15f %.15f %.15f @ %.15G C[brown]\n",
	    ratB[nn][0], ratB[nn][1], ratB[nn][2], Oparams.sigmaSticky*0.5);
#endif
#endif
  fprintf(bf,"%.15f %.15f %.15f @ 0.1 C[green]\n", rxC, ryC, rzC);
  fclose(bf);

}

extern void delete_events(int evIdA);
void ProcessCollision(void)
{
  int k;
  UpdateAtom(evIdA);
  UpdateAtom(evIdB);
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  MD_DEBUG10(calc_energy("prima"));
  MD_DEBUG38(printf("[BUMP] t=%.15G i=%d at=%d j=%d at=%d collCode=%d\n", 
		    Oparams.time,evIdA,evIdC, evIdB, evIdD, evIdE)); 

#ifdef MD_PATCHY_HE
  /* i primi due bit sono il tipo di event (uscit buca, entrata buca, collisione con core 
   * mentre nei bit restanti c'e' la particella con cui tale evento e' avvenuto */
  bumpSP(evIdA, evIdB, evIdC, evIdD, &W, evIdE);
#else
  bump(evIdA, evIdB, rxC, ryC, rzC, &W);
#endif
  MD_DEBUG10(calc_energy("dopo"));
  MD_DEBUG(store_bump(evIdA, evIdB));
  //ENDSIM=1;
  /*printf("qui time: %.15f\n", Oparams.time);*/
#ifdef MD_GRAVITY
  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
#else
  OprogStatus.lastcolltime[evIdA] = OprogStatus.lastcolltime[evIdB] = 
    lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
  do_check_negpairs = 1;
#ifdef MD_CALC_DPP
  store_last_u(evIdA);
  store_last_u(evIdB);
#endif
#ifdef MD_PATCHY_HE
  lastbump[evIdA].mol=evIdB;
  lastbump[evIdA].at = evIdC;
  lastbump[evIdB].mol=evIdA;
  lastbump[evIdB].at = evIdD;
  //printf("lastbump[%d].at=%d lastbump[%d].at=%d\n", evIdA, lastbump[evIdA].at, evIdB, lastbump[evIdB].at);
  lastbump[evIdA].type = evIdE;
  lastbump[evIdB].type = evIdE;
#else
  lastbump[evIdA]=evIdB;
  lastbump[evIdB]=evIdA;
#endif
#endif
#ifdef MD_HSVISCO
  OprogStatus.lastcoll = Oparams.time;
#endif
  if (OprogStatus.useNNL)
    {
      /* ricalcola i tempi di collisione con la NL */
      updrebuildNNL(evIdA);
      updrebuildNNL(evIdB);
      PredictEventNNL(evIdA, -1);
      PredictEventNNL(evIdB, evIdA);
    }
  else
    {
      PredictEvent(evIdA, -1);
      PredictEvent(evIdB, evIdA);
    }
  do_check_negpairs = 0;
}
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
	  if (OprogStatus.useNNL)
	    nebrTab[evIdA].r[k] -= L;
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
	  if (OprogStatus.useNNL)
	    nebrTab[evIdA].r[k] += L;
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
void ProcessCellCrossing(void)
{
#ifdef MD_GRAVITY
  int j; 
#endif
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
  MD_DEBUG34(printf("OLD cellCrossing evIdA=%d k=%d inCells=%d %d %d\n", evIdA, k, inCell[0][evIdA], inCell[1][evIdA], inCell[2][evIdA]));
#ifdef MD_GRAVITY
  j = evIdB - ATOM_LIMIT;
  if (j >= 100)
    {
      k = j - 100; 
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
    }
  else
    {
      k = j / 2;
      cellRange[j] = 0;
      ProcessCollWall();
    }
#else
  k = evIdB - 100 - ATOM_LIMIT; 
#if 0
  if (inCell[0][evIdA]> cellsx ||inCell[1][evIdA]> cellsy||inCell[2][evIdA]> cellsz) 
  printf("Cells(%d,%d,%d)\n", inCell[0][evIdA],inCell[1][evIdA],inCell[2][evIdA]);
#endif
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
#endif
  /* NOTA: ogni cella delle linked list deve poter contenere il parallelepipedo delle NNL
   * ed inoltre tali celle servono solo per costruire le NNL e non più per predire gli 
   * urti fra gli ellissoidi. */
  MD_DEBUG32(printf("i=%d PROCESS CELL CROSSING\n", evIdA));
  MD_DEBUG34(printf("NEW cellCrossing evIdA=%d k=%d inCells=%d %d %d\n", evIdA, k, inCell[0][evIdA], inCell[1][evIdA], inCell[2][evIdA]));
  if (OprogStatus.useNNL)
    PredictEventNNL(evIdA, evIdB);
  else
    PredictEvent(evIdA, evIdB);
  n = (inCell[2][evIdA] * cellsy + inCell[1][evIdA])*cellsx + 
    inCell[0][evIdA] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellList[evIdA] = cellList[n];
  cellList[n] = evIdA;
}
void velsBrown(double T)
{
  comvel_brown(T, Oparams.m); 
}

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
    {
      if (OprogStatus.useNNL)
	PredictEventNNL(n, -2); 
      else
	PredictEvent(n, -2); 
    }
}
void distanza(int ia, int ib)
{
  double dx, dy, dz;
  dx = rx[ia]-rx[ib];
  dy = ry[ia]-ry[ib];
  dz = rz[ia]-rz[ib];
  dx = dx - L*rint(dx/L);
  dy = dx - L*rint(dy/L);
  printf("dist(%d,%d): %f\n", ia, ib, sqrt(Sqr(dx)+Sqr(dy)+Sqr(dz)));
}
void rebuildLinkedList(void);
#ifdef MD_BIG_DT
void timeshift_variables(void)
{
  int i;
  if (OprogStatus.scalevel)
    OprogStatus.nextcheckTime -= OprogStatus.bigDt;
  OprogStatus.nextDt -= OprogStatus.bigDt;
  if (OprogStatus.intervalSum > 0.0)
    OprogStatus.nextSumTime -= OprogStatus.bigDt;
  //nextStoreTime viene calcolato opportunamente ogni volta quindi non va shiftato
  if (OprogStatus.storerate > 0.0)
    OprogStatus.nextStoreTime -= OprogStatus.bigDt;
  nextNNLrebuild -= OprogStatus.bigDt;
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (OprogStatus.useNNL)
	nebrTab[i].nexttime -= OprogStatus.bigDt;
      atomTime[i] -= OprogStatus.bigDt;
      lastcol[i] -= OprogStatus.bigDt;
      OprogStatus.lastcolltime[i] -= OprogStatus.bigDt;
#ifdef MD_HSVISCO
      OprogStatus.lastcoll -= OprogStatus.bigDt;
#endif
    }
}
void timeshift_calendar(void)
{
  int poolSize, id;
  poolSize = Oparams.parnum*OprogStatus.eventMult;
  /* parte da 1 perché tree[0] è solo l'inzio dell'albero e non un evento */
  for (id=1; id < poolSize; id++) 
    {
      if (treeUp[id] != -1)
	treeTime[id] -= OprogStatus.bigDt;
    } 
}
#endif
#ifdef EDHE_FLEX
extern int getnumbonds(int np, interStruct ts);
int first=1;
#endif
#ifdef EDHE_FLEX
int termination(void)
{
  int nb;
  interStruct ts;
  ts.type1 = 0;
  ts.type2 = 5;
  ts.spot1 = 4;
  ts.spot2 = 0;
  nb = getnumbonds(0,ts);
  ts.type1 = 1;
  ts.type2 = 5;
  ts.spot1 = 5;
  ts.spot2 = 0;
  nb += getnumbonds(1,ts);			
  if ((nb==1 || nb==2) && first)
    {
      double ti = 0;
      FILE* f;
      first = 0;
      ti = Oparams.time;	
#ifdef 	MD_BIG_DT
      ti += OprogStatus.refTime;	
#endif
      printf("1 BOND time = %.15G\n", ti);
      f=fopen("onebond.dat","w");
      fprintf(f,"%.15G\n", ti);
      fclose(f);
      return 0;
    }	 
  if (nb==2)
    {
      double ti = 0;
      FILE* f;
      if (first == 1)
	first = 0;
      ti = Oparams.time;	
#ifdef 	MD_BIG_DT
      ti += OprogStatus.refTime;	
#endif
      printf("2 BONDS time = %.15G ... FINISHED !\n", ti);
      ENDSIM=1;
      f=fopen("twobonds.dat","w");
      fprintf(f,"%.15G\n", ti);
      fclose(f);
      return 1;
    }
  return 0; 
}	
#endif

/* ============================ >>> move<<< =================================*/
void move(void)
{
  char fileop[1024], fileop2[1024]; 
  char fileop3[1024];
  FILE *bf;
  const char sepStr[] = "@@@\n";
  int i, innerstep=0;
#ifdef MD_GRAVITY
  int ii;
  double rzmax, zfact;
#endif
  int k, n;
  double timeold;
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
#ifdef EDHE_FLEX 
      if (termination())
	{
	  ENDSIM=1;
	  break;
	}
#endif
      innerstep++;
      if (OprogStatus.useNNL)
	timeold = Oparams.time;
      NextEvent();
      /* l'evento di ricostruzione della NNL è mantenuto fuori dal calendario degli eventi per
       * semplicità */
      if (OprogStatus.useNNL)
	{
	  if (Oparams.time >= nextNNLrebuild)
	    {
	      Oparams.time = nextNNLrebuild;
	      InitEventList();
	      rebuildNNL();
	      for (k = 0;  k < 3; k++)
		{
		  cellRange[2*k]   = - 1;
		  cellRange[2*k+1] =   1;
		}
	      for (n = 0; n < Oparams.parnum; n++)
		{
		  PredictEventNNL(n, -2); 
		}
	      //rebuildCalendar();
	      if (OprogStatus.intervalSum > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
	    	ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
	      if (OprogStatus.bigDt > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
	      if (OprogStatus.scalevel)
		ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	      continue;
#if 0
	      NextEvent();
	      if (Oparams.time >= nextNNLrebuild)
		{
		  printf("Le NNL devono essere aggiornat troppo frequentemente!\n");
		  exit(-1);
		}
#endif
	    }
	}
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
#if defined(MD_ABSORPTION) && 1
      else if (evIdB == ATOM_LIMIT+50)
	{
	  ProcessWallColl();
	}
#endif
#if defined(MD_GRAVITY) || defined(MD_EDHEFLEX_WALL)
      else if (evIdB >= ATOM_LIMIT + 100 || evIdB < ATOM_LIMIT + NDIM * 2)
	{
	  if (evIdE == MD_WALL)
	    ProcessWallColl();
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
	  sprintf(fileop2 ,"Store-%d-%d", 
      		  OprogStatus.KK, OprogStatus.JJ);
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
#ifdef MD_ASYM_ITENS
	      calc_omega(i);
#endif
	      update_MSDrot(i);
#ifdef MD_CALC_DPP
	      update_MSD(i);
	      store_last_u(i);
#endif
	      OprogStatus.lastcolltime[i] = Oparams.time;
	    }
	  R2u();
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
    	      system(fileop3);
	    }
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
	}
     else if (evIdB == ATOM_LIMIT + 10)
	{
	  UpdateSystem();
	  R2u();
	  if (Oparams.curStep == Oparams.totStep)
	    {
#ifdef EDHE_FLEX
	      //if (!globSaveAll || OprogStatus.stripStore)
	      saveFullStore("StoreFinal");
#endif
	      outputSummary();
	    }

#if 0
	    {
	      static double shift[3] = {0,0,0}, vecg[8], vecgNeg[8];
	      double d,r1[3], r2[3], alpha;
	      int distfail;
	      FILE* f;
	      static int first = 1;
	      shift[0] = L*rint((rx[0]-rx[1])/L);
	      shift[1] = L*rint((ry[0]-ry[1])/L);
	      shift[2] = L*rint((rz[0]-rz[1])/L);
	      MD_DEBUG(printf("[EVENT10] shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]));
#if 0
	      d=calcDist(Oparams.time, 0, 1, shift, r1, r2, &alpha, vecg, 1);
	      if (first)
		f = fopen("distPos.dat","w");
	      else
		f = fopen("distPos.dat","a");
	      fprintf(f,"%.15G %.15G %.15G %.15G %.15G %.15G\n", Oparams.time, d,vecg[0],vecg[1],vecg[2],vecg[4]);
	      fclose(f);
#endif
	      d=calcDistNeg(Oparams.time, 0, 0, 1, shift, r1, r2, &alpha, vecgNeg, 1);
	      //d=calcDistNegNeighPlane(0, Oparams.time, 2, r1, r2, vecgNeg, 0, 1, &distfail, 2);
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
	  for (i=0; i < Oparams.parnum; i++)
	    {
	      update_MSDrot(i);
#ifdef MD_CALC_DPP
	      update_MSD(i);
	      store_last_u(i);
#endif
	      OprogStatus.lastcolltime[i] = Oparams.time;
	    }
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
	      if (OprogStatus.useNNL)
		rebuildNNL();
	      rebuildCalendar();
	      if (OprogStatus.intervalSum > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      if (OprogStatus.scalevel)
		ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
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
	  timeshift_calendar();
	  timeshift_variables();
	  OprogStatus.refTime += OprogStatus.bigDt;
	  ScheduleEvent(-1, ATOM_LIMIT + 11, OprogStatus.bigDt);
	}
#endif
#if 0 && defined(MD_NNL)
      else if (evIdB == ATOM_LIMIT + 11)
	{
	  UpdateSystem();
	  for (i=0; i < Oparams.parnum; i++)
	    {
	      BuildNNL(i);
	      if (i==0 || nebrTab[i].nexttime < nltime)
		nltime = nebrTab[i].nexttime;
	    }
	  printf("rebuilt NNL next time=%.15G\n", nltime);
	  /* next complete update */
	  ScheduleEvent(-1, ATOM_LIMIT + 11, nltime); 
	}
#endif
#ifdef MD_GRAVITY
      else if (evIdB == ATOM_LIMIT + 9)
	{
	  UpdateSystem();
	  if (OprogStatus.taptau > 0.0)
	    {
	      if (OprogStatus.quenchend < 0.0)
		{
#if 0
		  if ((V - Vold)/V < OprogStatus.quenchtol)
		    {
		      printf("QUENCH DONE! %d\n", Oparams.curStep);
		      /* se l'energia potenziale è ormai stabile considera il quench finito */
		      OprogStatus.quenchend = Oparams.curStep;
		    }
#else
		  calcKVz();
		  OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
		  if ( (2.0*K/(3.0*((double)Oparams.parnum)-3.0)) < 
		       OprogStatus.quenchtol)
		    {
#endif
		      printf("QUENCH DONE! %lld\n", (long long int) Oparams.curStep);
		      OprogStatus.numquench++;
		      /* calcola e salva le misure quando finisce il quench */
		      calcRho();
		      save_rho();
		      calccmz();
		      save_rzcm();
		      OprogStatus.quenchend = Oparams.time;
		      comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
#if 1
		      calcKVz();
		      scalevels(Oparams.T, K, Vz);
#endif
		      MD_DEBUG3(printf("rzmax:%f\n", rzmax));
		      rzmax = -Lz2;
		      for (ii=0; ii < Oparams.parnum; ii++)
			{
			  if (rz[ii] > rzmax)
			    rzmax = rz[ii];
			}
		      if (Lz / (rzmax+Lz2) < OprogStatus.expandFact)
			zfact = Lz/(rzmax+Lz2);
		      else
			zfact = OprogStatus.expandFact;
		      for (ii=0; ii < Oparams.parnum; ii++)
			{
			  rz[ii] = zfact*(rz[ii]+Lz2)-Lz2;
			  rz[ii] += OprogStatus.rzup;
			  vz[ii] += OprogStatus.vztap; 
			}
		      rebuildLinkedList();
		      MD_DEBUG3(distanza(996, 798));
		      if (OprogStatus.useNNL)
			rebuildNNL();
		      rebuildCalendar();
		      if (OprogStatus.storerate > 0.0)
			ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
		      if (OprogStatus.intervalSum > 0.0)
			ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
		      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
		      if (OprogStatus.bigDt > 0.0)
			ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
		    }
		}
	      else if ((Oparams.time - OprogStatus.quenchend)  < OprogStatus.taptau)
		{
		  /* se scalevelsteps = 0 allora scala ogni passo se si sta facendo il 
		     tapping */
		  OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
		  calcKVz();
		  MD_DEBUG4(printf("SCALVEL #%lld Vz: %.15f\n", (long long int) Oparams.curStep,Vz));
		  scalevels(Oparams.T, K, Vz);
		  rebuildLinkedList();
		  if (OprogStatus.useNNL)
		    rebuildNNL();
		  rebuildCalendar();
		  if (OprogStatus.storerate > 0.0)
		    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
		  if (OprogStatus.intervalSum > 0.0)
		    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
		  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
		  if (OprogStatus.bigDt > 0.0)
		    ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
		}
	      else
		{
		  /* start quench (-1  significa che il quench è iniziato) */
		  OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
		  OprogStatus.quenchend = -1;
		}
	      ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	    }
	  else if (OprogStatus.scalevel)
	    {
	      OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	      calcKVz();
	      MD_DEBUG2(printf("[TAPTAU < 0] SCALVEL #%lld Vz: %.15f\n", (long long int)Oparams.curStep,Vz));
	      scalevels(Oparams.T, K, Vz);
	      if (OprogStatus.useNNL)
		updAllNNL();
	      rebuildCalendar();
	      if (OprogStatus.intervalSum > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
	      if (OprogStatus.bigDt > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
	    }
#if 0
	  else if (2.0*K/(3.0*Oparams.parnum-3.0)>Oparams.T)
	    {
	      UpdateSystem();
	      calcKVz();
	      scalevels(Oparams.T, K, Vz);
	    }
#endif
	}
      if (OprogStatus.maxquench && OprogStatus.numquench == OprogStatus.maxquench)
	ENDSIM = 1;
#else
      else if (evIdB == ATOM_LIMIT + 9)
	{
	  if (OprogStatus.scalevel)
	    {
	      UpdateSystem();
	      OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	      MD_DEBUG2(printf("[TAPTAU < 0] SCALVEL #%lld Vz: %.15f\n", 
			       (long long int)Oparams.curStep,Vz));
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
#ifdef EDHE_FLEX
	      calc_energy_filtered(2);
#else
	      calc_energy(NULL);
#endif
	      scalevels(Oparams.T, K);
	      if (OprogStatus.useNNL)
		updAllNNL();
	      rebuildCalendar();
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	      if (OprogStatus.rescaleTime > 0)
		ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	      else
		OprogStatus.scalevel = 0;
#ifdef MD_BIG_DT
	      if (OprogStatus.bigDt > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
	    }
	}
#endif
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
      if (mgl_mode==2)
	{
	  ENDSIM=1;
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

  if ((OprogStatus.rmsd2end > 0.0 && OprogStatus.tmsd2end > 0.0 &&
       DphiSqA > Sqr(OprogStatus.rmsd2end) 
       && DrSqTotA > Sqr(OprogStatus.tmsd2end)) 
      || (OprogStatus.rmsd2end > 0.0  && OprogStatus.tmsd2end <= 0.0 
	  && DphiSqA > Sqr(OprogStatus.rmsd2end))
      || (OprogStatus.tmsd2end > 0.0  && OprogStatus.rmsd2end <= 0.0 
	  && DrSqTotA > Sqr(OprogStatus.tmsd2end))
      )
    {
      printf("[MSDcheck] steps %d time %.15G\n", Oparams.curStep, Oparams.time);
      ENDSIM=1;
    }
  if (ENDSIM)
    {
#ifdef EDHE_FLEX
      //if (!globSaveAll || OprogStatus.stripStore)
      saveFullStore("StoreFinal");
#endif
      R2u();
    }
}
