#ifdef MD_PATCHY_HE
#include<mdsimul.h>
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG10(x)  
#define MD_DEBUG11(x) 
#define MD_DEBUG15(x) 
#define MD_DEBUG20(x) 
#define MD_DEBUG29(x) 
#define MD_DEBUG30(x)  //qui 
#define MD_DEBUG31(x)  //qui 
#define MD_DEBUG32(x) 
#define MD_NEGPAIRS
#define MD_NO_STRICT_CHECK
#define MD_OPTDDIST
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **REtA, **REtB;
extern double cosEulAng[2][3], sinEulAng[2][3];
extern long long int itsFNL, timesFNL, timesSNL, itsSNL;
#ifdef MD_ASYM_ITENS
extern double **Ia, **Ib, **invIa, **invIb, **Iatmp, **Ibtmp, *angM;
#else
extern double Ia, Ib, invIa, invIb;
#endif
extern double gradplane[3];
struct LastBumpS *lastbump;
extern double *axa, *axb, *axc;
extern int *scdone;
extern double *maxax;
/* Routines for LU decomposition from Numerical Recipe online */
void ludcmpR(double **a, int* indx, double* d, int n);
void lubksbR(double **a, int* indx, double *b, int n);
extern void update_MSDrot(int i);
extern void InvMatrix(double **a, double **b, int NB);
extern void calc_angmom(int i, double **I);
double min(double a, double b);
extern void upd_refsysM(int i);
extern double calc_maxddot(int i, int j);
double zbrent(double (*func)(double), double x1, double x2, double tol);
extern double invaSq[2], invbSq[2], invcSq[2];
extern double rxC, ryC, rzC;
extern int SolveLineq (double **a, double *x, int n); 
extern double calc_norm(double *vec);
int calcdist_retcheck;
extern double rA[3], rB[3];
extern double gradplane_all[6][3], rBall[6][3];
extern int polinterr, polinterrRyck;
/* *** change here if you change the number sticky spots *** */
int mapbondsaAB[MD_PBONDS]={1,1,2,2,3,3,4,4,5,5};
int mapbondsbAB[MD_PBONDS]={1,2,1,2,1,2,1,2,1,2};
int *mapbondsa;
int *mapbondsb;
/* ------------------------------------------------------------ */
extern void bump (int i, int j, double rCx, double rCy, double rCz, double* W);
extern int *crossevtodel;
extern long long int itsF, timesF, itsS, timesS, numcoll;
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern double accngA, accngB;
void ScheduleEventBarr (int idA, int idB, int idata, int atb, int idcollcode, double tEvent);
double calcDistNeg(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, double dists[MD_PBONDS], int bondpair);
extern void symtop_evolve_orient(int i, double ti, double **Ro, double **REt, double cosea[3], double sinea[3], double *phir, double *psir);

void comvel_brown (COORD_TYPE temp, COORD_TYPE *m);
void remove_bond(int na, int n, int a, int b);
void add_bond(int na, int n, int a, int b);
void writeAsciiPars(FILE* fs, struct pascii strutt[]);
void writeAllCor(FILE* fs);
void InitEventList (void);
extern void UpdateOrient(int i, double ti, double **Ro, double Omega[3][3]);
int bound(int na, int n, int a, int b);
#ifdef MD_HSVISCO
void calcT(void);
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
double calcDistNegOneSP(double t, double t1, int i, int j, int nn, double shift[3]);
extern double pi, invL, L2, Vz;   
extern double W, K, T1xx, T1yy, T1zz,
       T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz, 
       Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB, DQxxOld, 
       DQyyOld, DQzzOld, DQxyOld, DQyzOld, DQzxOld, DQxxOldKin, 
       DQyyOldKin, DQzzOldKin, DQxxOldHS, DQyyOldHS, DQzzOldHS, DQxxOldST, DQyyOldST, DQzzOldST,
       PxxKin, PyyKin, PzzKin, PxxHS, PyyHS, PzzHS, PxxST, PyyST, PzzST;
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
//extern double DrSq = 0.0; 
//extern const double timbig = 1E12;
/* used by linked list routines */
extern double *lastcol;
extern double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
extern int  **tree, cellRange[2*NDIM], initUcellx, initUcelly, initUcellz;
extern int *inCell[3], *cellList, cellsx, cellsy, cellsz;
extern int evIdA, evIdB, parnumB, parnumA, evIdD, evIdE;
extern int evIdC;
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
double rcutL, aL, bL, cL;
extern double max_ax(int i);
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[]);
#define MD_SP_DELR 0.0
double spApos[MD_STSPOTS_A][3] = {{MD_SP_DELR, 0.54, 0.0},{MD_SP_DELR, 0.54, 3.14159},{MD_SP_DELR, 2.60159,0.0},
    {MD_SP_DELR, 2.60159, 3.14159},{MD_SP_DELR, 1.5708, 0.0}};
double spBpos[MD_STSPOTS_B][3] = {{MD_SP_DELR, 0.0, 0.0},{MD_SP_DELR, 3.14159, 0.0}};

double spXYZ_A[MD_STSPOTS_A][3];
double spXYZ_B[MD_STSPOTS_B][3];
void build_atom_positions(void)
{
 /* N.B. le coordinate spXpos sono del tipo (Dr, theta, phi),
  * dove se Dr=0 la sfera sticky viene posizionata esattamente in 
  * maniera tangente e theta (0 <= theta <= Pi) e phi (0 <= phi < 2Pi)
  * sono gli angoli in coordinate sferiche che individuano il punto di contatto
  * tra sticky sphere ed ellissoide.
  * Tale routine converte le coordinate spXpos in coordinate cartesiane 
  * riferite al riferimento del corpo rigido. */  
  int kk, k1, aa;
  double x,y,z, grad[3], ng, dd[3];

  spApos[0][1] = Oparams.theta;
  spApos[1][1] = Oparams.theta;
  spApos[2][1] = pi - Oparams.theta;
  spApos[3][1] = pi - Oparams.theta;
  spApos[0][0] = Oparams.Dr;
  spApos[1][0] = Oparams.Dr;
  spApos[2][0] = Oparams.Dr;
  spApos[3][0] = Oparams.Dr;
  spApos[4][0] = Oparams.Dr;
  spBpos[0][0] = Oparams.Dr;
  spBpos[1][0] = Oparams.Dr;
  for (k1 = 0; k1 < MD_STSPOTS_A; k1++)
    {
      x = Oparams.a[0]*cos(spApos[k1][2])*sin(spApos[k1][1]);
      y = Oparams.b[0]*sin(spApos[k1][2])*sin(spApos[k1][1]);
      z = Oparams.c[0]*cos(spApos[k1][1]);
      //printf("xyz=%f %f %f\n", x, y, z);
      grad[0] = 2.0 * x / Sqr(Oparams.a[0]);
      grad[1] = 2.0 * y / Sqr(Oparams.b[0]);
      grad[2] = 2.0 * z / Sqr(Oparams.c[0]);
      ng = calc_norm(grad);
      for (aa = 0; aa < 3; aa++)
	grad[aa] /= ng;
      spXYZ_A[k1][0] = x + grad[0]*(Oparams.sigmaSticky*0.5 + spApos[k1][0]);
      spXYZ_A[k1][1] = y + grad[1]*(Oparams.sigmaSticky*0.5 + spApos[k1][0]);
      spXYZ_A[k1][2] = z + grad[2]*(Oparams.sigmaSticky*0.5 + spApos[k1][0]);
	
	      //printf("k1=%d %f %f %f \n", k1,  spXYZ_A[k1][0] ,    spXYZ_A[k1][1] ,  spXYZ_A[k1][2]  );
    }
  for (kk=0; kk < 3; kk++)
    dd[kk] = spXYZ_A[0][kk] - spXYZ_A[1][kk];
  printf("Molecule A distance between Atoms 0 and 1: %.15G\n", calc_norm(dd));;
  for (kk=0; kk < 3; kk++)
    dd[kk] = spXYZ_A[2][kk] - spXYZ_A[3][kk];
  printf("Molecule A distance between Atoms 2 and 3: %.15G\n", calc_norm(dd));;

  for (k1 = 0; k1 < MD_STSPOTS_B; k1++)
    {
      x = Oparams.a[1]*cos(spBpos[k1][2])*sin(spBpos[k1][1]);
      y = Oparams.b[1]*sin(spBpos[k1][2])*sin(spBpos[k1][1]);
      z = Oparams.c[1]*cos(spBpos[k1][1]);
      grad[0] = 2.0 * x / Sqr(Oparams.a[1]);
      grad[1] = 2.0 * y / Sqr(Oparams.b[1]);
      grad[2] = 2.0 * z / Sqr(Oparams.c[1]);
      ng = calc_norm(grad);
      for (aa = 0; aa < 3; aa++)
	grad[aa] /= ng;
      spXYZ_B[k1][0] = x + grad[0]*(Oparams.sigmaSticky*0.5 + spBpos[k1][0]);
      spXYZ_B[k1][1] = y + grad[1]*(Oparams.sigmaSticky*0.5 + spBpos[k1][0]);
      spXYZ_B[k1][2] = z + grad[2]*(Oparams.sigmaSticky*0.5 + spBpos[k1][0]) ;
    }
}
extern void tRDiagR(int i, double **M, double a, double b, double c, double **Ri);
extern void calc_energy(char *msg);
extern void print_matrix(double **M, int n);
int getnumbonds(int np, int at)
{
  int kk, jj, jj2, aa, bb, nb;
  nb=0;
  for (kk = 0; kk < numbonds[np]; kk++)
    {
      jj = bonds[np][kk] / (NA*NA);
      jj2 = bonds[np][kk] % (NA*NA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (aa == at)
	nb++;
    }
  return nb;
}
int one_is_bonded(int i, int a, int j, int b)
{
  /* per ora è disbilitato */
  //return 0;
  if (getnumbonds(i,a) >= Oparams.nmax || getnumbonds(j,b) >= Oparams.nmax)
    return 1;
  else
    return 0;
}
void bumpSP(int i, int j, int ata, int atb, double* W, int bt)
{
  /* NOTA: Controllare che inizializzare factor a 0 è corretto! */
  double factor=0, invmi, invmj, mredl;
  double delpx, delpy, delpz, wrx, wry, wrz, rACn[3], rBCn[3];
  double rAB[3], rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double ratA[3], ratB[3], norm[3];
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, taus, DTxx, DTyy, DTzz;
#endif
  double denom, rCx, rCy, rCz, nrAB, Dr;
#ifndef MD_ASYM_ITENS
  double factorinvIa, factorinvIb;
#endif
#ifdef MD_ASYM_ITENS
  int k1,k2, b;
  double rnI[3];
  double Mvec[3], omega[3];
#endif
  int na, a, kk;
#if 0
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
  numcoll++;
  //printf("collision code: %d (%d,%d)\n", bt, i, j);
  MD_DEBUG31(calc_energy("PRIMA"));
  if (bt == MD_CORE_BARRIER)
    {
      bump(i, j, rxC, ryC, rzC, W);
      MD_DEBUG31(calc_energy("DOPO HARD COLL"));
      MD_DEBUG10(printf(">>>>>>>>>>collCode: %d\n", bt));
      MD_DEBUG30(printf("time=%.15G collision type= %d %d-%d %d-%d ata=%d atb=%d\n",Oparams.time, bt, i, j, j, i, ata, atb));
      return;
    }
  rA[0] = rx[i];
  rA[1] = ry[i];
  rA[2] = rz[i];
  rB[0] = rx[j];
  rB[1] = ry[j];
  rB[2] = rz[j];
  for (a=0; a < 3; a++)
    {
      Dr = rA[a]-rB[a];
      if (fabs(Dr) > L2)
	{
	  rB[a] += SignR(L, Dr);
	}
    }
  MD_DEBUG10(printf("[bump] t=%f contact point: %f,%f,%f \n", Oparams.time, rxC, ryC, rzC));
  /* qui calcolo il punto di contatto */
  MD_DEBUG20(printf("ata: %d atb: %d\n", ata, atb));
  MD_DEBUG20(printf("rA %f %f %f\n", rA[0], rA[1], rA[2]));
  MD_DEBUG20(printf("rB %f %f %f\n", rB[0], rB[1], rB[2]));
  BuildAtomPosAt(i, ata, rA, R[i], ratA);
  BuildAtomPosAt(j, atb, rB, R[j], ratB);
  //printf("ata:%d atb:%d\n", ata, atb);
  MD_DEBUG20(printf("ratA %f %f %f\n", ratA[0], ratA[1], ratA[2]));
  MD_DEBUG20(printf("ratB %f %f %f\n", ratB[0], ratB[1], ratB[2]));
  for (kk = 0; kk < 3; kk++)
    rAB[kk] = ratA[kk] - ratB[kk];
  /* reduce to minimum image rAB[]!! */
 nrAB = calc_norm(rAB);
 MD_DEBUG(printf("sigmaSticky= %.15G norm rAB: %.15G\n", Oparams.sigmaSticky, calc_norm(rAB)));
 for (kk = 0; kk < 3; kk++)
    rAB[kk] /= nrAB;
  /* controllare con cura la scelta dei parametri relativi ai diametri delle sferette
   * e alle larghezze delle buche dei potenziali a buca quadrata */
  MD_DEBUG20(printf("coll code: %d\n", bt));
  rCx = ratA[0] - rAB[0]*Oparams.sigmaSticky*0.5;
  rCy = ratA[1] - rAB[1]*Oparams.sigmaSticky*0.5;
  rCz = ratA[2] - rAB[2]*Oparams.sigmaSticky*0.5;
  rAC[0] = rA[0] - rCx;
  rAC[1] = rA[1] - rCy;
  rAC[2] = rA[2] - rCz;
 
  rBC[0] = rB[0] - rCx;
  rBC[1] = rB[1] - rCy;
  rBC[2] = rB[2] - rCz;
  MD_DEBUG(printf("\n"));
  /* calcola tensore d'inerzia e le matrici delle due quadriche */
  na = (i < Oparams.parnumA)?0:1;
  if (OprogStatus.targetPhi > 0)
    {
      /* scalare tutti i raggi qui */
    }
#ifdef MD_ASYM_ITENS
  tRDiagR(i, Ia, Oparams.I[na][0], Oparams.I[na][1], Oparams.I[na][2], R[i]);
#else
  Ia = Oparams.I[na];
#endif
  na = (j < Oparams.parnumA)?0:1;
  if (OprogStatus.targetPhi > 0)
    {
      /* scalare tutti i raggi qui */
    }
#ifdef MD_ASYM_ITENS
  tRDiagR(j, Ib, Oparams.I[na][0], Oparams.I[na][1], Oparams.I[na][2], R[j]);
#else
  Ib = Oparams.I[na];
#endif
  MD_DEBUG(check_contact(evIdA, evIdB, Xa, Xb, rAC, rBC));
  /* calcola le matrici inverse del tensore d'inerzia */
#ifdef MD_ASYM_ITENS
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Iatmp[k1][k2] = Ia[k1][k2];
	Ibtmp[k1][k2] = Ib[k1][k2];
      } 
  InvMatrix(Iatmp, invIa, 3);
  InvMatrix(Ibtmp, invIb, 3);
  //printf("invIa=\n");
  //print_matrix(invIa, 3);
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

  
  invmi = (i<Oparams.parnumA)?1/Oparams.m[0]:1/Oparams.m[1];
  invmj = (j<Oparams.parnumA)?1/Oparams.m[0]:1/Oparams.m[1];

  denom = invmi + invmj; 
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  MD_DEBUG20(printf("[bump] before bump vc=%.15G\n", vc));
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
  mredl = 1.0 / denom;
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
      if (Sqr(vc) < 2.0*(Oparams.bheight+Oparams.bhout)/mredl)
	{
	  MD_DEBUG31(printf("MD_INOUT_BARRIER (%d,%d)-(%d,%d) t=%.15G vc=%.15G NOT ESCAPEING collType: %d d=%.15G\n",  i, ata, j, atb, 
			    Oparams.time, vc,  bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
	  factor = -2.0*vc;
	}
      else
	{
	  MD_DEBUG31(printf("_MD_INOUT_BARRIER (%d-%d)-(%d,%d) t=%.15G vc=%.15G ESCAPING collType: %d d=%.15G\n", i, ata, j, atb, Oparams.time, vc, bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
	  factor = -vc + sqrt(Sqr(vc) - 2.0*Oparams.bheight/mredl);
	  remove_bond(i, j, ata, atb);
	  remove_bond(j, i, atb, ata);
	}
#if 0
      printf("factor: %.15G mredl=%.15G vc=%.15G sqrt(XXX=%.15G) bheight: %.15G\n",
		 factor, mredl, vc, sqrt(Sqr(vc) - 2.0*Oparams.bheight/mredl), Oparams.bheight);
#endif	 
      factor *= mredl;
      break;
    case MD_OUTIN_BARRIER:
      if (one_is_bonded(i, ata, j, atb) || (Oparams.bhin >= 0.0 && Sqr(vc) < 2.0*Oparams.bhin/mredl))
	{
	  MD_DEBUG31(printf("MD_INOUT_BARRIER (%d,%d)-(%d,%d) t=%.15G vc=%.15G NOT ESCAPEING collType: %d d=%.15G\n",  i, ata, j, atb, 
			    Oparams.time, vc,  bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
	  factor = -2.0*vc;
	}
      else
	{
	  add_bond(i, j, ata, atb);
	  add_bond(j, i, atb, ata);
	  factor = -vc - sqrt(Sqr(vc) + 2.0*Oparams.bheight/mredl);
	}
	MD_DEBUG31(printf("[MD_OUTIN_BARRIER] (%d,%d)-(%d,%d)  delta= %f height: %f mredl=%f\n", 
		      i, ata, j, atb, Sqr(vc) + 2.0*Oparams.bheight/mredl, Oparams.bheight, mredl));
#if 0
	{ double dist;
	  double shift[3]={0,0,0};
	  printf("mapbondsa[1]:%d mapbondsb[1]:%d\n", mapbondsa[1], mapbondsb[1]);
	  dist = calcDistNegOneSP(Oparams.time, 0.0, i, j, 7, shift);
	  printf("dists[7]=%.15G\n", dist);
	  dist = calcDistNegOneSP(Oparams.time, 0.0, i, j, 1, shift);
	  printf("dists[1]:%.15G\n", dist);
	}
#endif
      factor *= mredl;
      break;
    }
  MD_DEBUG(printf("factor=%f denom=%f\n", factor, denom));
  delpx = factor * norm[0];
  delpy = factor * norm[1];
  delpz = factor * norm[2];
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
      OprogStatus.DQWxy += (rA[0]-rB[0])*delpy;
      OprogStatus.DQWyz += (rA[1]-rB[1])*delpz;
      OprogStatus.DQWzx += (rA[2]-rB[2])*delpx;
      OprogStatus.DQWxx += (rA[0]-rB[0])*delpx;
      OprogStatus.DQWyy += (rA[1]-rB[1])*delpy;
      OprogStatus.DQWzz += (rA[2]-rB[2])*delpz;
      OprogStatus.DQWxxST += (rA[0]-rB[0])*delpx;
      OprogStatus.DQWyyST += (rA[1]-rB[1])*delpy;
      OprogStatus.DQWzzST += (rA[2]-rB[2])*delpz;
    }
  OprogStatus.Txy += DTxy; 
  OprogStatus.Tyz += DTyz;
  OprogStatus.Tzx += DTzx;
  OprogStatus.Txx += DTxx; 
  OprogStatus.Tyy += DTyy;
  OprogStatus.Tzz += DTzz;
#endif
  MD_DEBUG(printf("delp=(%f,%f,%f)\n", delpx, delpy, delpz));
#ifdef MD_ASYM_ITENS
  factor = -factor;
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
#endif
#ifdef MD_ASYM_ITENS
  calc_angmom(i, Ia);
  upd_refsysM(i);
  calc_angmom(j, Ib);
  upd_refsysM(j);
#endif
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, vx[i],vy[i],vz[i], j, vx[j],vy[j],vz[j]));
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, wx[i],wy[i],wz[i], j, wx[j],wy[j],wz[j]));

  MD_DEBUG31(calc_energy("DOPO"));
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
void assign_bond_mapping(int i, int j)
{
  /* NOTA: l'interazione bonded è solo tra Si e O 
   * i <  Oparams.parnumA => O
   * i >=  Oparams.parnumA => Si */
  if (i < Oparams.parnumA)
    {
      mapbondsa = mapbondsaAB;
      mapbondsb = mapbondsbAB;
    }
  else 
    {
      mapbondsb = mapbondsaAB;
      mapbondsa = mapbondsbAB;
    }
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
  MD_DEBUG31(printf("numbonds[%d]=%d bonds[][numbonds-1]:%d a=%d b=%d\n", na, numbonds[na],bonds[na][numbonds[na]-1],
  a, b));
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
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3])
{
  /* QUESTA VA RISCRITTA PER GLI ELLISSOIDI STICKY!!! */
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int k1, k2;
  double *spXYZ=NULL;
  //double radius; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */

  if (ata > 0)
    {
      if (i < Oparams.parnumA)
	spXYZ = spXYZ_A[ata-1];
      else  
	spXYZ = spXYZ_B[ata-1];
    }
  //radius = Oparams.sigma[0][1] / 2.0;
  if (ata == 0)
    {
      for (k1 = 0; k1 < 3; k1++)
	rat[k1] = rO[k1];
      //printf("%f %f %f @ 0.5 C[red]\n", rat[0], rat[1], rat[2]);
    }
  else 
    {
      for (k1 = 0; k1 < 3; k1++)
	{ 
	  rat[k1] = rO[k1];
	  for (k2 = 0; k2 < 3; k2++)
	    rat[k1] += R[k2][k1]*spXYZ[k2]; 
	}
      //printf("ata= %d rat= %f %f %f\n", ata, rat[0], rat[1], rat[2]);
      //printf("rO = %f %f %f \n", rO[0], rO[1], rO[2]);
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
      //printf("ata=%d %f %f %f @ 0.075 C[blue]\n", ata, R[0][ata-1], R[1][ata-1], R[2][ata-1]);
    }
  
}
void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a;
  /* l'atomo zero si suppone nell'origine */
  if (i < Oparams.parnumA)
    {
      for (a=0; a < MD_STSPOTS_A+1; a++)
	BuildAtomPosAt(i, a, rO, R, rat[a]);
    }
  else
    {
      for (a=0; a < MD_STSPOTS_B+1; a++)
	BuildAtomPosAt(i, a, rO, R, rat[a]);
    }
}
int ibr, jbr, nnbr; 
double shiftbr[3], trefbr;
double calcDistNegOneSP(double t, double t1, int i, int j, int nn, double shift[3]);

double funcs2beZeroedSP(double x, double tref, int i, int j, int nn, double shift[3])
{
  return calcDistNegOneSP(x, trefbr, i, j, nn, shift);
}

double  funcs2beZeroedBrent(double x)
{
  return funcs2beZeroedSP(x, trefbr, ibr, jbr, nnbr, shiftbr); 
}

extern double sigmaSqSticky;
double calcDistNegOneSP(double t, double t1, int i, int j, int nn, double shift[3])
{
  double distSq, ti;
  double ratA[NA][3], ratB[NA][3];
  int kk;
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
  int na;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  assign_bond_mapping(i, j);
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  MD_DEBUG20(printf("BRENT nn=%d\n", nn));
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  //UpdateOrient(i, ti, RtA, Omega, mapbondsa[nn]);
  UpdateOrient(i, ti, RtA, Omega);
#endif
  /* calcola le posizioni nel laboratorio degli atomi della molecola */
  BuildAtomPos(i, rA, RtA, ratA);
  na = (i < Oparams.parnumA)?0:1;
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(j, ti, RtB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  //UpdateOrient(j, ti, RtB, Omega, mapbondsb[nn]);
  UpdateOrient(j, ti, RtB, Omega);
#endif
  na = (j < Oparams.parnumA)?0:1;
  BuildAtomPos(j, rB, RtB, ratB);
  /* calcola sigmaSq[][]!!! */
  distSq = 0;
  for (kk=0; kk < 3; kk++)
    distSq += Sqr(ratA[mapbondsa[nn]][kk]-ratB[mapbondsb[nn]][kk]);
  MD_DEBUG20(printf("dist= %.15G\n", sqrt(distSq)-Oparams.sigmaSticky));
  return sqrt(distSq) - Oparams.sigmaSticky;
}

/* N.B. per la silica tale routine va cambiata! */
double calcDistNegSP(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, 
		   double dists[MD_PBONDS], int bondpair)
{
  double distmin, distSq, ti;
  double ratA[NA][3], ratB[NA][3], dist;
  int firstdist = 1, nn, kk;
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
  int na;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  //UpdateOrient(i, ti, RtA, Omega, (bondpair==-1)?-1:mapbondsa[bondpair]);
  UpdateOrient(i, ti, RtA, Omega);
#endif
  /* calcola le posizioni nel laboratorio degli atomi della molecola */
  BuildAtomPos(i, rA, RtA, ratA);
  na = (i < Oparams.parnumA)?0:1;
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(j, ti, RtB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  //UpdateOrient(j, ti, RtB, Omega, (bondpair==-1)?-1:mapbondsb[bondpair]);
  UpdateOrient(j, ti, RtB, Omega);
#endif
  na = (j < Oparams.parnumA)?0:1;
  BuildAtomPos(j, rB, RtB, ratB);
  /* calcola sigmaSq[][]!!! */
  distmin = 0;
  for (nn = 0; nn < MD_PBONDS; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	{
	  //printf("qui in calcDistNeg\n");
	  continue;
	}
      distSq = 0;
      for (kk=0; kk < 3; kk++)
	distSq += Sqr(ratA[mapbondsa[nn]][kk]-ratB[mapbondsb[nn]][kk]);
      dists[nn] = dist = sqrt(distSq) - Oparams.sigmaSticky;
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
int refine_contactSP(int i, int j, double tref, double t1, double t2, int nn, double shift[3], double *troot)
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
int check_crossSP(double distsOld[MD_PBONDS], double dists[MD_PBONDS], 
		int crossed[MD_PBONDS], int bondpair)
{
  int nn;
  int retcross = 0;
  for (nn = 0; nn < MD_PBONDS; nn++)
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
int get_dists_tocheckSP(double distsOld[], double dists[], int tocheck[], int dorefine[],
		      int bondpair)
{
  int nn;
  int rettochk = 0;
  for (nn = 0; nn < MD_PBONDS; nn++)
    {
      tocheck[nn] = 0;
      if ( dists[nn]*distsOld[nn] > 0.0 &&
	  fabs(dists[nn]) < OprogStatus.epsdSP && fabs(distsOld[nn]) < OprogStatus.epsdSP &&
	  dorefine[nn] == MD_EVENT_NONE && (bondpair== -1 || bondpair == nn))
	{
	  tocheck[nn] = 1; 
	  rettochk++;
	}
    }
  return rettochk;
}
double get_max_deldistSP(double distsOld[MD_PBONDS], double dists[MD_PBONDS], int bondpair)
{
  int nn, first = 1;
  double maxdd=0.0, dd;
  for (nn = 0; nn < MD_PBONDS; nn++)
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

void assign_distsSP(double a[], double b[])
{
  memcpy(b, a, MD_PBONDS*sizeof(double));
}
/* NOTA: tale stima ottimizzata della maggiorazione per la velocità di variazione della distanza
 * sembra corretta, fare comunque dei test.*/
double eval_maxddistSP(int i, int j, int bondpair, double t1, double *maxddotOpt)
{
  double ti, rA[3], rB[3], ratA[NA][3], ratB[NA][3], wri[3], wrj[3], nwri, nwrj,
	 r12i[3], r12j[3];//, maxddotOpt[MD_PBONDS];
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  double maxddot=0.0, nr12i, nr12j;
  int nn, kk;
  ti = t1 - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  //UpdateOrient(i, ti, RtA, Omega, (bondpair==-1)?-1:mapbondsa[bondpair]);
  UpdateOrient(i, ti, RtA, Omega);
#endif
  BuildAtomPos(i, rA, RtA, ratA);
  ti = t1 - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti;
  rB[1] = ry[j] + vy[j]*ti;
  rB[2] = rz[j] + vz[j]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(j, ti, RtB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  //UpdateOrient(j, ti, RtB, Omega, (bondpair==-1)?-1:mapbondsb[bondpair]);
  UpdateOrient(j, ti, RtB, Omega);
#endif
  BuildAtomPos(j, rB, RtB, ratB);
  for (nn = 0; nn < MD_PBONDS; nn++)
    {
      for (kk = 0; kk < 3; kk++)
	{
	  r12i[kk] = (ratA[mapbondsa[nn]][kk]-rA[kk]);
  	  r12j[kk] = (ratB[mapbondsb[nn]][kk]-rB[kk]);	  
	}
      nr12i = calc_norm(r12i);
      nr12j = calc_norm(r12j);
      for (kk = 0; kk < 3; kk++)
	{
	  r12i[kk] *= (nr12i+OprogStatus.epsdSP)/nr12i;
	  r12j[kk] *= (nr12j+OprogStatus.epsdSP)/nr12j;
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
void calc_deltSP(double *maxddoti, double *delt, double *dists, int bondpair)
{
  int nn;
  double dt;
  for (nn = 0; nn < MD_PBONDS; nn++)
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

int search_contact_fasterSP(int i, int j, double *shift, double *t, double t1, double t2, double epsd, double *d1, double epsdFast, double dists[MD_PBONDS], int bondpair, double maxddot, double *maxddoti)
{
  /* NOTA: 
   * MAXOPTITS è il numero massimo di iterazioni al di sopra del quale esce */
  double told, delt, distsOld[MD_PBONDS];
  const int MAXOPTITS = 500;
  int its=0, amin, bmin, crossed[MD_PBONDS]; 
  /* estimate of maximum rate of change for d */
#if 0
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]*0.5
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j]*0.5;
#endif
  *d1 = calcDistNegSP(*t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
  MD_DEBUG30(printf("[IN SEARCH CONTACT FASTER]*d1=%.15G t=%.15G bondpair=%d\n", *d1, *t, bondpair));
  timesF++;
  MD_DEBUG10(printf("Pri distances between %d-%d d1=%.12G epsd*epsdTimes:%f\n", i, j, *d1, epsdFast));
  told = *t;
  delt = OprogStatus.h;
  if (fabs(*d1) < epsdFast)
    {
      assign_distsSP(distsOld, dists);
      return 0;
    }
  while (fabs(*d1) > epsdFast && its < MAXOPTITS)
    {
      if (maxddot*(t2-(t1+*t)) < fabs(*d1)-OprogStatus.epsdSP)
	{
	  MD_DEBUG30(printf("SP maxddot*(t2-(t1+*t)=%.15G fabs(*d1)=%.15G\n", maxddot*(t2-(t1+*t)),fabs(*d1)));
	  MD_DEBUG30(printf("SP maxddot=%.15G t2-(t1+*t)=%.15G\n",maxddot, t2-(t1+*t)));
	  MD_DEBUG30(printf("SP t2=%.15G t1=%.15G *t=%.15G\n", t2, t1, *t));
	  return 1;
	}
#ifdef MD_OPTDDIST
      calc_deltSP(maxddoti, &delt, distsOld, bondpair);
#else
      delt = fabs(*d1) / maxddot;
#endif
      *t += delt;
      *d1 = calcDistNegSP(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
      if (check_crossSP(distsOld, dists, crossed, bondpair))
	{
	  /* go back! */
	  MD_DEBUG30(printf("SP d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("SP d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNegSP(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 0;
	}
      if (*t+t1 > t2)
	{
	  *t = told;
	  MD_DEBUG30(printf("SP t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("SP convergence t>t2\n"));
	  *d1 = calcDistNegSP(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 1;
	}
      told = *t;
      assign_distsSP(dists, distsOld);
      its++;
      itsF++;
    }

  MD_DEBUG30(printf("SP max iterations %d iterations reached t=%f t2=%f\n", its, *t, t2));
  return 0;
}
extern double **Aip;
extern void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double xa[3], ya[3];
double distfuncSP(double x)
{
  double dy, y;
  polint(xa, ya, 3, x, &y, &dy);
  if (polinterr==1)
    return 0.0;
  if (dy > OprogStatus.epsdSP)
    {
      if (OprogStatus.phitol <= 0)
	printf("dy=%.15G\n", dy);
      polinterr = 1;
    }
  else 
    polinterr = 0;
  return y;
}
int interpolSP(int i, int j, int nn, 
	     double tref, double t, double delt, double d1, double d2,
	     double *tmin, double shift[3], int ignoresignchg)
{
  double d3, A, dmin;
  /* NOTA: dists di seguito può non essere usata? controllare!*/
  d3 = calcDistNegOneSP(t+delt*0.5, tref, i, j, nn, shift);
  xa[0] = 0;
  ya[0] = d1;
  xa[1] = delt*0.5;
  ya[1] = d3;
  xa[2] = delt;
  ya[2] = d2;
#if 0
  A = xa[2]*(ya[0]-ya[1]);
  B = xa[0]*(ya[1]-ya[2]);
  C = xa[1]*(ya[2]-ya[0]);
  *tmin = (xa[2]*A+xa[0]*B+xa[1]*C)/(A+B+C)/2.0;
#endif
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
  dmin = calcDistNegOneSP(*tmin, tref, i, j, nn, shift);
  if (*tmin < t+delt && *tmin > t && (d1*dmin < 0.0 || ignoresignchg) )
    {
      *tmin += tref;
      return 0;
    }
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
  int nb, jj, jj2, kk, nn, aa, bb;
  nb = numbonds[i];
  if (!OprogStatus.assumeOneBond)
    return -1;
  for (kk = 0; kk < nb; kk++)
    {
      jj = bonds[i][kk] / (NA*NA);
      jj2 = bonds[i][kk] % (NA*NA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if (jj == j)
	{
	  for (nn=0; nn < MD_PBONDS; nn++)
	    if (mapbondsa[nn]==aa && mapbondsb[nn]==bb)
	      return nn;
	} 
    }
  return -1; 

}

int check_negpairs(int *negpairs, int bondpair, int i, int j)
{
  int nn, sum;
  sum = 0;
//  if (lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == 0 
  //    && lastbump[j].at == 0)
    //return 2;
  for (nn = 0; nn < MD_PBONDS; nn++)
    {
      negpairs[nn] = 0;
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (!(lastbump[i].mol == j && lastbump[j].mol==i && lastbump[i].at == mapbondsa[nn]
	&& lastbump[j].at == mapbondsb[nn]))
	continue;
      negpairs[nn] = 1;
      return nn+1;
      //sum += 1;
    }
  return 0;
  //return (sum > 0)?1:0;
}

int delt_is_too_big_hc(int i, int j, int bondpair, double *dists, double *distsOld)
{
  int nn;
  for (nn=0; nn < MD_PBONDS; nn++)
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
  int nn;
  for (nn=0; nn < MD_PBONDS; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (!negpairs[nn])
	continue;
      /* N.B. distsOld[nn] non va controllato poiché dopo un urto cmq ci deve essere un estremo
       * di distanza dmin t.c. dmin > 0 se !bound(i,j..) o dmin < 0 se bound(i,j...) */
      if (dists[nn] >= 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
      if (dists[nn] <= 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
	return 1;
    }
  return 0;
}
#ifdef MD_ASYM_ITENS
double calc_maxddotSP(int i, int j, double *maxddoti)
{
  int kk;
#if 0
  int na;
  double Iamin, Ibmin;
#endif
  double factori, factorj, maxddot=0.0;

  for (kk=0; kk < MD_PBONDS; kk++)
    {
      if (i < Oparams.parnumA)
	factori = calc_norm(spXYZ_A[mapbondsa[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
      else
	factori = calc_norm(spXYZ_B[mapbondsa[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
      if (j < Oparams.parnumA)
	factorj = calc_norm(spXYZ_A[mapbondsb[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
      else
	factorj = calc_norm(spXYZ_B[mapbondsb[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
#if 0
	factori = 0.5*maxax[i]+OprogStatus.epsdSP;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
	factorj = 0.5*maxax[j]+OprogStatus.epsdSP;//sqrt(Sqr(axa[j])+Sqr(axb[j])+Sqr(axc[j]));
#endif
#if 0
	na = i<Oparams.parnumA?0:1;
	Iamin = min(Oparams.I[na][0],Oparams.I[na][2]);
	na = j<Oparams.parnumA?0:1;
	Ibmin = min(Oparams.I[na][0],Oparams.I[na][2]);
	maxddoti[kk] = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
	  angM[i]*factori/Iamin + angM[j]*factorj/Ibmin;
#else
	maxddoti[kk] = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
	  sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori + 
	  sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;

#endif
	if (kk==0 || maxddoti[kk] > maxddot)
	  {
	    maxddot = maxddoti[kk];
	  }
    }
  return maxddot;
}
#else
double calc_maxddotSP(int i, int j, double *maxddoti)
{
  int kk;
  double maxddot=0.0;
  double factori, factorj;

  for (kk=0; kk < MD_PBONDS; kk++)
    {
#if 0
      if (i < Oparams.parnumA)
	printf("kk=%d mapbondsa[kk]=%d mapbondsb[kk]:%d norm: %.15G w = (%.15G, %.15G, %.15G)\n",kk, mapbondsa[kk], mapbondsb[kk],
	       calc_norm(spXYZ_A[mapbondsa[kk]-1]), wx[i], wy[i], wz[i]);
      else
	printf("kk=%d mapbondsa[kk]=%d mapbondsb[kk]:%d norm: %.15G w = (%.15G, %.15G, %.15G)\n",
	       kk, mapbondsa[kk], mapbondsb[kk],
	       calc_norm(spXYZ_B[mapbondsa[kk]-1]), wx[i], wy[i], wz[i]);
#endif
      if (i < Oparams.parnumA)
	factori = calc_norm(spXYZ_A[mapbondsa[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
      else
	factori = calc_norm(spXYZ_B[mapbondsa[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
      if (j < Oparams.parnumA)
	factorj = calc_norm(spXYZ_A[mapbondsb[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
      else
	factorj = calc_norm(spXYZ_B[mapbondsb[kk]-1]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
#if 0
      factori = 0.5*maxax[i]+OprogStatus.epsdSP;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
      factorj = 0.5*maxax[j]+OprogStatus.epsdSP;//sqrt(Sqr(axa[j])+Sqr(axb[j])+Sqr(axc[j]));
#endif
      maxddoti[kk] = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j]))+
	sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori + 
	sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*factorj;  
      
      if (kk==0 || maxddoti[kk] > maxddot)
	{
	  maxddot = maxddoti[kk];
	}
      //printf("maxddot:%.15G maxddoti:%.15G\n", maxddot, maxddoti[kk]);
    }
  return maxddot;
}
#endif

int locate_contactSP(int i, int j, double shift[3], double t1, double t2, 
		   double *evtime, int *ata, int *atb, int *collCode)
{
  const double minh = 1E-20;
  double h, d, dold, t2arr[MD_PBONDS], t, dists[MD_PBONDS], distsOld[MD_PBONDS];
  double maxddot, delt, troot, tmin, tini; //distsOld2[MD_PBONDS];
#ifndef MD_BASIC_DT
  double deldist, normddot, dold2, distsOld2[MD_PBONDS]; 
#endif
  //const int MAXOPTITS = 4;
  int bondpair, itstb;
  int its, foundrc;
  double maxddoti[MD_PBONDS], epsd, epsdFast, epsdFastR, epsdMax; 
  int tocheck[MD_PBONDS], dorefine[MD_PBONDS], ntc, ncr, nn, gotcoll, amin, bmin,
      crossed[MD_PBONDS], firstaftsf;
#ifdef MD_NEGPAIRS
  int negpairs[MD_PBONDS], sumnegpairs;
#endif
  const double GOLD= 1.618034;
  epsd = OprogStatus.epsdSP;
  epsdFast = OprogStatus.epsdFastSP;
  epsdFastR= OprogStatus.epsdFastSP;
  epsdMax = OprogStatus.epsdSP;
  assign_bond_mapping(i, j);
  bondpair = get_bonded(i, j);
  t = 0.0;
#ifdef MD_OPTDDIST
#ifndef MD_ASYM_ITENS
  maxddot = eval_maxddistSP(i, j, bondpair, t1, maxddoti);
#else
  maxddot = calc_maxddotSP(i, j, maxddoti);
#endif
#else
  //maxddot = calc_maxddotSP(i, j, maxddoti);
#ifndef MD_ASYM_ITENS
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]*0.5
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j]*0.5;
#else
  maxddot = calc_maxddot(i, j);
#endif
#endif
  MD_DEBUG10(printf("[locate_contact] %d-%d t1=%f t2=%f shift=(%f,%f,%f)\n", i,j,t1, t2, shift[0], shift[1], shift[2]));
  h = OprogStatus.h; /* last resort time increment */
  if (*collCode!=MD_EVENT_NONE)
    {
      if (t2 > *evtime)
	t2 = *evtime+1E-7;
    }
  delt = h;
  MD_DEBUG(printf("QUIIII collCode=%d\n", *collCode));
#ifndef MD_NEGPAIRS
  /* NOTA: le strategie per evitare problemi dopo una collisione sono due:
   * 1) andare avanti nel tempo finché la distanza non è corretta.
   * 2) fare un passo ed eventualmente ridurlo finchè la distanza non è corretta.
   */
  df = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
  for (nn=0; nn < MD_PBONDS; nn++)
    {
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
	      df = calcDistNegOneSP(t, t1, i, j, nn, shift);
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
	      df = calcDistNegOneSP(t, t1, i, j, nn, shift);
	    }
	}
    }
#endif
  MD_DEBUG30(printf("[BEFORE SEARCH CONTACT FASTER_SP]Dopo distances between %d-%d t=%.15G t2=%.15G\n", i, j, t, t2));
#ifdef MD_NEGPAIRS
  sumnegpairs = check_negpairs(negpairs, bondpair, i, j); 
#endif
  if (search_contact_fasterSP(i, j, shift, &t, t1, t2, epsd, &d, epsdFast, dists, bondpair, maxddot, maxddoti))
    {
      return 0;  
    }
  timesS++;
  MD_DEBUG30(printf("[AFTER SEARCH CONTACT FASTER_SP]Dopo distances between %d-%d d1=%.12G\n", i, j, d));

  MD_DEBUG(printf(">>>>d:%f\n", d));
  foundrc = 0;
#if 1
  assign_distsSP(dists, distsOld);
  dold = d;
#else
  dold = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
#endif
  firstaftsf = 1;
  its = 0;
  while (t+t1 < t2)
    {
#ifdef MD_BASIC_DT
      delt = epsd/maxddot;
      tini = t;
      t += delt;
      d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
#else
      if (!firstaftsf)
	{
	  deldist = get_max_deldistSP(distsOld2, distsOld, bondpair);
	  normddot = fabs(deldist)/delt;
	  /* NOTA: forse qui si potrebbe anche usare sempre delt = epsd/maxddot */
	  if (normddot!=0)
	    {
	      delt = epsd/normddot;
	    }
	  else
	    {
	      delt = epsd/maxddot;
	      //delt = h;
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
	  dold2 = calcDistNegSP(t-delt, t1, i, j, shift, &amin, &bmin, distsOld2, bondpair);
	  MD_DEBUG30(printf("SP ==========>>>>> t=%.15G t2=%.15G\n", t, t2));
	  continue;
	}
      MD_DEBUG30(printf("SP delt: %f epsd/maxddot:%f h*t:%f maxddot:%f\n", delt, epsd/maxddot,h*t,maxddot));
      tini = t;
      t += delt;
      d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);

      deldist = get_max_deldistSP(distsOld, dists, bondpair);
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
	  t += delt; 
	  itsS++;
	  d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	}
#endif
#ifdef MD_NEGPAIRS
      itstb = 0;
      /* NOTA: se la distanza tra due sticky spheres è positiva a t (per errori numerici 
       * accade spesso) e t+delt allora delt è troppo grande e qui lo riduce fino ad un 
       * valore accettabile. */
      if (sumnegpairs)// && !firstaftsf)
	{
	  if (delt_is_too_big(i, j, bondpair, dists, distsOld, negpairs))
	    {
	      if(!interpolSP(i, j, sumnegpairs-1, t1, tini, delt, distsOld[sumnegpairs-1], 
			     dists[sumnegpairs-1], &tmin, shift, 1))
		{
		  //printf("qui\n");
		  tmin -= t1;
		  delt = tmin - tini;
		  t = tmin;
		  d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
		}
#if 1
	      else 
		{
		  while (delt_is_too_big(i, j, bondpair, dists, distsOld, negpairs) && 
			 delt > minh)
		    {
		      delt /= GOLD; 
		      t = tini + delt;
		      d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
		      itstb++;
		      if (!interpolSP(i, j, sumnegpairs-1, t1, tini, delt, distsOld[sumnegpairs-1], 
				      dists[sumnegpairs-1], &tmin, shift, 1))
			{
			  tmin -= t1;
			  delt = tmin - tini;
			  t = tmin;
			  d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
			  break;
			}
		    }
		}
#endif
	    } 
	  sumnegpairs = 0;
	}

#if 0
	  
	  while (delt_is_too_big(i, j, bondpair, dists, distsOld, negpairs) && 
		 delt > minh)
	    {
	      delt /= GOLD; 
	      t = tini + delt;
	      d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	      itstb++;
	    }
	  sumnegpairs = 0;
	}
#endif
#endif
      MD_DEBUG30(printf(">>>>> d = %.15G\n", d));
      for (nn=0; nn < MD_PBONDS; nn++)
	dorefine[nn] = MD_EVENT_NONE;
      ncr=check_crossSP(distsOld, dists, crossed, bondpair);
      /* N.B. crossed[] e tocheck[] sono array relativi agli 8 possibili tipi di attraversamento fra gli atomi
       * sticky */
      for (nn = 0; nn < MD_PBONDS; nn++)
	{
	  t2arr[nn] = t; 

	  dorefine[nn] = MD_EVENT_NONE;
	  if (crossed[nn]!=MD_EVENT_NONE)
	    {
	      /* se dorefine è 2 vuol dire che due superfici si sono
	       * attraversate */
	      if (valid_collision(i, j, mapbondsa[nn], mapbondsb[nn], crossed[nn]))
		{
		  MD_DEBUG30(printf("SP type: %d i=%d j=%d ata=%d atb=%d bound:%d\n", crossed[nn], i, j, mapbondsa[nn],
			 mapbondsb[nn], bound(i, j, mapbondsa[nn], mapbondsb[nn])));
		    dorefine[nn] = crossed[nn];
		}
	    }
	}

#define MD_INTERPOL
#ifdef MD_INTERPOL
      ntc = get_dists_tocheckSP(distsOld, dists, tocheck, dorefine, bondpair);
      for (nn = 0; nn < MD_PBONDS; nn++)
	{
	  if (tocheck[nn])
	    {
	      //printf("tocheck[%d]:%d\n", nn, tocheck[nn]);
	      if (interpolSP(i, j, nn, t1, t-delt, delt, distsOld[nn], dists[nn], 
			   &troot, shift, 0))
		dorefine[nn] = MD_EVENT_NONE;
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
      for (nn = 0; nn < MD_PBONDS; nn++)
	{
	  if (dorefine[nn]!=MD_EVENT_NONE)
	    {
	      MD_DEBUG30(printf("REFINE dorefine[%d]:%d\n", nn, dorefine[nn]));
	      if (refine_contactSP(i, j, t1, t-delt, t2arr[nn], nn, shift, &troot))
		{
		  //printf("[locate_contact] Adding collision between %d-%d\n", i, j);
		  MD_DEBUG30(printf("[locate_contact_sp] Adding collision between %d-%d\n", i, j));
		  MD_DEBUG30(printf("[locate_contact_sp] t=%.15G nn=%d\n", t, nn));
		  MD_DEBUG30(printf("[locate_contact_sp] troot=%.15G\n", troot));
		  MD_DEBUG(printf("[locate_contact_sp] its: %d\n", its));
		  /* se il legame già c'è e con l'urto si forma tale legame allora
		   * scarta tale urto */
		  if (troot > t2 || troot < t1
#if 1
		  || 
		      (lastbump[i].mol == j && lastbump[j].mol==i && 
		       lastbump[i].at == mapbondsa[nn]
		       && lastbump[j].at == mapbondsb[nn] && fabs(troot - lastcol[i]) < 1E-14))
#endif
		    {
		     MD_DEBUG31(printf("SP lastbump[%d].mol=%d lastbump[%d].at=%d lastbump[%d].mol=%d lastbump[%d].at=%d\n", i, lastbump[i].mol, i, lastbump[i].at, j, lastbump[j].mol, j, lastbump[j].at));
		      //gotcoll = -1;
		      continue;
		    }
		  else
		    {
		      gotcoll = 1;
		      MD_DEBUG31(printf("SP *evtime=%.15G troot=%.15G troot-*evtime:%.15G\n", *evtime, troot, 
					troot-*evtime));
		      if (*collCode == MD_EVENT_NONE || troot < *evtime)
			{
			  *ata = mapbondsa[nn];
			  *atb = mapbondsb[nn];
			  *evtime = troot;
			  *collCode = dorefine[nn]; 
			  MD_DEBUG31(printf("SP ok scheduling collision between %d-%d nn=%d\n", i, j, nn));
			  MD_DEBUG31(printf("SP collcode=%d bound(i, j, nn):%d\n", *collCode, 
					    bound(i, j, mapbondsa[nn], mapbondsb[nn])));
			}
		      continue;
		    }
		}
	      else 
		{
		  MD_DEBUG(printf("[locate_contactSP] can't find contact point!\n"));
#ifdef MD_INTERPOL
		  if (!tocheck[nn])
#endif
		  mdPrintf(ALL,"[locate_contactSP] can't find contact point!\n",NULL);
		  /* Se refine_contact fallisce deve cmq continuare a cercare 
		   * non ha senso smettere...almeno credo */
		  //gotcoll = -1;
		  continue;
		}
	    }
	}
      if (gotcoll == 1)
	return 1;
      else if (gotcoll == -1)
	return 0;
      if (fabs(d) > epsdFastR)
	{
	  if (search_contact_fasterSP(i, j, shift, &t, t1, t2, epsd, &d, epsdFast, dists, bondpair,
				    maxddot, maxddoti))
	    {
	      MD_DEBUG30(printf("[search contact faster locate_contact_sp] d: %.15G\n", d));
	      return 0;
	    }
#if 1
	  dold = d;
	  assign_distsSP(dists, distsOld);
#else
	  dold = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
#endif
	  firstaftsf = 1;
	  its++;
	  //itsS++;
	  continue;
	}
      dold = d;
      MD_DEBUG30(printf("SP ==========>>>>> t=%.15G t2=%.15G\n", t, t2));
#ifndef MD_BASIC_DT
      assign_distsSP(distsOld,  distsOld2);
#endif
      assign_distsSP(dists, distsOld);
      its++;
      itsS++;
    }
  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return 0;
}
/* -------- >>> neighbour list stuff <<< --------- */
double get_max_deldist_sp(int nsp, double distsOld[6][NA], double dists[6][NA])
{
  int nn, nn2, first = 1;
  double maxdd=0.0, dd;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < nsp; nn2++)
	{

	  dd = fabs(dists[nn][nn2]-distsOld[nn][nn2]);
	  if (first || dd > maxdd)
	    {
	      first = 0;
	      maxdd = dd;
	    }
	}
    }
  return maxdd;
}

double calcDistNegOneNNL_sp(double t, double t1, int i, int nn);
int interpolNeighPlane_sp(int i, double tref, double t, double delt, double d1, double d2, double *tmin, int nplane, int nn)
{
  double d3, A, dmin;
  /* NOTA: dists di seguito può non essere usata? controllare!*/
  d3 = calcDistNegOneNNL_sp(t+delt*0.5, tref, i, nn);
#if 0
  xa[0] = t;
  ya[0] = d1;
  xa[1] = t+delt*0.5;
  ya[1] = d3;
  xa[2] = t+delt;
  ya[2] = d2;
  A = xa[2]*(ya[0]-ya[1]);
  B = xa[0]*(ya[1]-ya[2]);
  C = xa[1]*(ya[2]-ya[0]);
  *tmin = (xa[2]*A+xa[0]*B+xa[1]*C)/(A+B+C)/2.0;
#endif
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
  else if (ya[2]-ya[0] == 0.0)
    {
      *tmin = t + delt*0.5;
    }
  else 
    {
      A = (ya[2]-ya[0])/(ya[0]-ya[1]);
      *tmin = t + 0.5*delt*((1.0 + A * 0.25)/( 1.0 + A * 0.5));
    }
  dmin = calcDistNegOneNNL_sp(*tmin, tref, i, nn);
  if (*tmin < t+delt && *tmin > t && d1*dmin < 0.0)
    {
      *tmin += tref;
      return 0;
    }
  return 1;
}

int check_cross_sp(int nsp, double distsOld[6][NA], double dists[6][NA], int crossed[6][NA])
{
  int nn, nn2;
  int retcross = 0;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < nsp; nn2++)
	{

	  crossed[nn][nn2] = 0;
	  //printf("dists[%d]=%.15G distsOld[%d]:%.15G\n", nn, dists[nn], nn, distsOld[nn]);
	  if (dists[nn][nn2] < 0.0  && distsOld[nn][nn2] > 0.0)
	    {
	      crossed[nn][nn2] = 1;
	      //printf("CROSSED[%d]:%d\n", nn, crossed[nn]);
	      retcross = 1;
	    }
	}
    }
  return retcross;
}

int check_cross_scf_sp(int NSP, double distsOld[6][NA], double dists[6][NA], int crossed[6][NA])
{
  int nn, nn2;
  int retcross = 0;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < NSP; nn2++)
	{
	  crossed[nn][nn2] = 0;
	  //printf("dists[%d]=%.15G distsOld[%d]:%.15G\n", nn, dists[nn], nn, distsOld[nn]);
	  if ((fabs(dists[nn][nn2]) < 1E-15 || dists[nn][nn2] < 0.0) && distsOld[nn][nn2] > 0.0)
	    {
	      crossed[nn][nn2] = 1;
	      retcross = 1;
	    }
	}
    }
  return retcross;
}
void assign_dists_sp(int nsp, double a[6][NA], double b[6][NA])
{
  int k1, k2;
  //memcpy(b, a, nsp*6*sizeof(double));
  for (k1=0; k1 < 6; k1++)
    for (k2 = 0; k2 < nsp; k2++)
      b[k1][k2] = a[k1][k2];
}
int get_dists_tocheck_sp(int nsp, double distsOld[6][NA], double dists[6][NA], int tocheck[6][NA], int dorefine[6][NA])
{
  int nn, nn2;
  int rettochk = 0;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < nsp; nn2++)
	{
	  tocheck[nn][nn2] = 0;
	  if (dists[nn][nn2] < OprogStatus.epsdSPNL && distsOld[nn][nn2] < OprogStatus.epsdSPNL &&
	      dorefine[nn][nn2] == 0)
	    {
	      tocheck[nn][nn2] = 1; 
	      rettochk++;
	    }
	}
    }
  return rettochk;
}
double calcDistNegOneNNL_sp_norient(double t, double t1, int i, int nn, double ratA[NA][3])
{
  double dist;
  int kk;
  dist = 0;
  for (kk=0; kk < 3; kk++)
    dist += -(ratA[nn+1][kk]-rB[kk])*gradplane[kk];
  MD_DEBUG32(printf("DIST NOORIENT nn=%d t=%.15G dist=%.15G\n", nn, t+t1, dist- Oparams.sigmaSticky*0.5));
  return dist - Oparams.sigmaSticky*0.5;
}

double calcDistNegNeighPlaneAll_sp(int nsp, double t, double t1, int i, double dists[6][NA])
{
  int nn, kk, nn2;
  double dmin=0.0;
  double ti;
  double ratA[NA][3];
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  MD_DEBUG20(printf("BRENT nn=%d\n", nn));
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  //UpdateOrient(i, ti, RtA, Omega, mapbondsa[nn]);
  UpdateOrient(i, ti, RtA, Omega);
#endif
  /* calcola le posizioni nel laboratorio degli atomi della molecola */
  BuildAtomPos(i, rA, RtA, ratA);

  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < nsp; nn2++)
	{
	  for (kk = 0; kk < 3; kk++)
	    {
	      gradplane[kk] = gradplane_all[nn][kk];
	      rB[kk] = rBall[nn][kk];
	    }
	  dists[nn][nn2] = calcDistNegOneNNL_sp_norient(t, t1, i, nn2, ratA);
	  //printf("dist[%d]:%.15G\n", nn, dists[nn]);
	  if ((nn==0 && nn2==0) || dists[nn][nn2] < dmin)
	    dmin = dists[nn][nn2];
	}
    }
  return dmin; 
}
int check_distance_sp(int nsp, double maxddoti[6][NA], double dists[6][NA], double t1, double t2, double t)
{
  int nn, cc, nn2;
  cc = 0;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < nsp; nn2++)
	{
	  if ((t2 - (t1 + t))*maxddoti[nn][nn2] <  dists[nn][nn2] - OprogStatus.epsdSPNL)
	    cc++;
	}
    }
  if (cc == 6*nsp)
    return 1;
  else
    return 0;
  //printf("I chose dt=%.15G\n", *delt);
}
void calc_delt_sp(int nsp, double maxddoti[6][NA], double *delt, double dists[6][NA])
{
  int nn, nn2;
  double dt;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < nsp; nn2++)
	{
	  dt = fabs(dists[nn][nn2]) / maxddoti[nn][nn2];
	  //printf("nn=%d dt=%.15G delt=%.15G dists=%.15G maxddoti=%15G\n", nn, dt, *delt, dists[nn], maxddoti[nn]);
	  if ((nn==0 && nn2 ==0) || dt < (*delt))
	    *delt = dt;
	}
    }
  //printf("I chose dt=%.15G\n", *delt);
}
extern double max(double a, double b);
extern const double mddotfact;
void adjust_maxddoti_sp(int i, int NSP, double *maxddot, double maxddotiLC[6][NA], double maxddoti[6][NA])
{
  double K = 1.0;
  int a, b;
#ifdef MD_ASYM_ITENS
  if (Mx[i] == 0.0 && My[i] == 0.0 && Mz[i] == 0.0)
    K = mddotfact;
#else
  if (wx[i] == 0.0 && wy[i] == 0.0 && wz[i] == 0.0)
    K = mddotfact;
#endif
  *maxddot *= K;
  for (a = 0; a < 6; a++)
    for (b = 0; b < NSP; b++)
      maxddoti[a][b] = K*maxddotiLC[a][b];

}
int search_contact_faster_neigh_plane_all_sp(int i, double *t, double t1, double t2, 
					  double epsd, double *d1, double epsdFast, 
					  double dists[6][NA], double maxddotiLC[6][NA], double maxddot)
{
  double told, delt=1E-15, distsOld[6][NA];
  const double GOLD= 1.618034;
  const int MAXOPTITS = 500;
  double maxddoti[6][NA];
  int its=0, crossed[6][NA], itsf, NSP; 

  if (i < Oparams.parnumA)
    NSP = MD_STSPOTS_A;
  else
    NSP = MD_STSPOTS_B;
  *d1 = calcDistNegNeighPlaneAll_sp(NSP, *t, t1, i, distsOld);
#if 0
  if ((t2-t1)*maxddot < *d1 - OprogStatus.epsdSPNL)
    return 1;
#endif
  timesFNL++;
  told = *t;
  if (fabs(*d1) < epsdFast)
    {
      assign_dists_sp(NSP, distsOld, dists);
      return 0;
    }
  adjust_maxddoti_sp(i, NSP, &maxddot, maxddotiLC, maxddoti);
  while (fabs(*d1) > epsdFast && its < MAXOPTITS)
    {
#if 1
      calc_delt_sp(NSP, maxddoti, &delt, distsOld);
      if (check_distance_sp(NSP, maxddoti, dists, t1, t2, *t))
	return 1;
#else
      delt = fabs(*d1) / maxddot;
      if (*t + t1 < t2 && (t2 - (*t + t1))*maxddot < fabs(*d1) - OprogStatus.epsdSPNL)
	return 1;
#endif
      *t += delt;
      *d1 = calcDistNegNeighPlaneAll_sp(NSP, *t, t1, i, dists);
#if 0
      if (its > 100 && its%10 == 0)
	{
	  printf("NNL SEARCH CONTACT FASTER t=%.15G its=%d\n", *t+t1, its);
	}
#endif
      //printf("d=%.15G t=%.15G\n", *d1, *t+t1);
#if 1
      itsf = 0;
      while (check_cross_sp(NSP, distsOld, dists, crossed))
	{
	  /* reduce step size */
	  if (itsf == 0 && delt - OprogStatus.h > 0)
	    delt -= max(OprogStatus.h, OprogStatus.h*delt);
	  else
	    delt /= GOLD;
	  *t = told + delt;
	  *d1 = calcDistNegNeighPlaneAll_sp(NSP, *t, t1, i, dists);
	  itsf++;	
	  if (itsf > 100)
	    {
	      printf("*d1=%.15G too many times calculation of distance failed!\n", *d1);
	      printf("aborting...\n");
	      exit(-1);
	    }
	}
#else
     if (check_cross_sp(NSP, distsOld, dists, crossed))
       {
	 /* go back! */
	 MD_DEBUG30(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	 MD_DEBUG30(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	 *t = told;	  
	 *d1 = calcDistNegNeighPlaneAll_sp(NSP, *t, t1, i, dists);
	 return 0;
       }
#endif
      if (*t+t1 > t2)
	{
	  *t = told;
	  MD_DEBUG30(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("convergence t>t2\n"));
	  *d1 = calcDistNegNeighPlaneAll_sp(NSP, *t, t1, i, dists);
	  return 1;
	}
      told = *t;
      assign_dists_sp(NSP, dists, distsOld);
      its++;
      itsFNL++;
    }
  return 0;
}
double calcDistNegOneNNL_sp(double t, double t1, int i, int nn)
{
  double dist, ti;
  double ratA[NA][3];
  int kk;
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  MD_DEBUG20(printf("BRENT nn=%d\n", nn));
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  //UpdateOrient(i, ti, RtA, Omega, mapbondsa[nn]);
  UpdateOrient(i, ti, RtA, Omega);
#endif
  /* calcola le posizioni nel laboratorio degli atomi della molecola */
  BuildAtomPos(i, rA, RtA, ratA);
  /* calcola sigmaSq[][]!!! */
  dist = 0;
  for (kk=0; kk < 3; kk++)
    dist += -(ratA[nn+1][kk]-rB[kk])*gradplane[kk];
  MD_DEBUG32(printf("BRENT nn=%d t=%.15G dist= %.15G\n", nn,
		    t1+t,  dist - Oparams.sigmaSticky*0.5));
  return dist - Oparams.sigmaSticky*0.5;
}
double funcs2beZeroedNNL_sp(double x, double tref, int i, int nn)
{
  return calcDistNegOneNNL_sp(x, trefbr, i, nn);
}

double  funcs2beZeroedBrentNNL_sp(double x)
{
  return funcs2beZeroedNNL_sp(x, trefbr, ibr, nnbr); 
}

int refine_contact_neigh_plane_sp(int i, double tref, double t1, double t2, double *troot,
			       int nplane, int nsp)
{
  polinterr=0;
  //newt(vecg, 5, &retcheck, funcs2beZeroed, i, j, shift); 
  ibr = i;
  nnbr = nsp;
  trefbr = tref;
  *troot=zbrent(funcs2beZeroedBrentNNL_sp, t1, t2, 1E-16);
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
extern const double timbig;
extern void calc_grad_and_point_plane_all(int i, double gradplaneALL[6][3], double rBALL[6][3]);
extern void assign_plane(int nn);
#ifdef MD_ASYM_ITENS
double calc_maxddot_nnl_sp(int i, int nn, double *gradplane)
{
#if 0
  int na;
  double Iamin;
#endif
  double factori;
#if 0
  factori = 0.5*maxax[i]+OprogStatus.epsdSP;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
#else
  if (i < Oparams.parnumA)
    factori = calc_norm(spXYZ_A[nn]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
  else
    factori = calc_norm(spXYZ_B[nn]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
#endif
#if 0
  na = i<Oparams.parnumA?0:1;
  Iamin = min(Oparams.I[na][0],Oparams.I[na][2]);
  return fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
     angM[i]*factori/Iamin;
#else
  return fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
}
#else
double calc_maxddot_nnl_sp(int i, int nn, double *gradplane)
{
  double factori;
#if 0
  factori = 0.5*maxax[i]+OprogStatus.epsdSP;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
#else
   if (i < Oparams.parnumA)
    factori = calc_norm(spXYZ_A[nn]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
  else
    factori = calc_norm(spXYZ_B[nn]) + 0.5*Oparams.sigmaSticky + OprogStatus.epsdSP;
#endif
  return fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
}
#endif


int locate_contact_neigh_plane_parall_sp(int i, double *evtime, double t2)
{
  /* const double minh = 1E-14;*/
  double h, d, dold, t2arr[6][NA], t, dists[6][NA], distsOld[6][NA]; 
  double maxddot, delt, troot, tini, maxddoti[6][NA];
#ifndef MD_BASIC_DT
  double distsOld2[6][NA], dold2, normddot, deldist;
#endif
  int firstev, nn2;
  /*
  const int MAXITS = 100;
  const double EPS=3E-8;*/ 
  /* per calcolare derivate numeriche questo è il magic number in doppia precisione (vedi Num. Rec.)*/
  int its, foundrc, NSP;
  double t1, epsd, epsdFast, epsdFastR, epsdMax; 
  int tocheck[6][NA], dorefine[6][NA], ntc, ncr, nn, gotcoll, crossed[6][NA], firstaftsf;
  epsd = OprogStatus.epsdSPNL;
  epsdFast = OprogStatus.epsdFastSPNL;
  epsdFastR= OprogStatus.epsdFastSPNL;
  epsdMax = OprogStatus.epsdSPNL;
  t = 0;//t1;
  t1 = Oparams.time;
  //t2 = timbig;
  if (i < Oparams.parnumA)
    NSP = MD_STSPOTS_A;
  else
    NSP = MD_STSPOTS_B;
  calc_grad_and_point_plane_all(i, gradplane_all, rBall);
  //factori = 0.5*maxax[i]+OprogStatus.epsdSPNL;
  maxddot = 0.0;
  for (nn = 0; nn < 6; nn++)
    {
      for (nn2 = 0; nn2 < NSP; nn2++)
	{
	  maxddoti[nn][nn2] = calc_maxddot_nnl_sp(i, nn2, gradplane_all[nn]);
	  //printf("maxddoti[%d][%d]=%.15G\n", nn, nn2, maxddoti[nn][nn2]);
#if 0
	  maxddoti[nn][nn2] = fabs(vx[i]*gradplane_all[nn][0]+vy[i]*gradplane_all[nn][1]+vz[i]*gradplane_all[nn][2])+
	    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
	  if ((nn==0 && nn2==0) || maxddoti[nn][nn2] > maxddot)
	    maxddot = maxddoti[nn][nn2];
	  //printf("nn=%d maxddoti=%.15G\n", nn, maxddoti);
	}
    }
  h = OprogStatus.h; /* last resort time increment */
  delt = h;
  
  if (search_contact_faster_neigh_plane_all_sp(i, &t, t1, t2, epsd, &d, epsdFast, 
					       dists, maxddoti, maxddot))
    {
      return 0;  
    }

  MD_DEBUG30(printf("t=%.15G d=%.15G\n", t, d));
  timesSNL++;
  foundrc = 0;
  assign_dists_sp(NSP, dists, distsOld);
  dold = d;
  firstaftsf = 1;
  its = 0;
  while (t+t1 < t2)
    {
#if 0
      if (its > 500 && its%10 == 0)
	printf("[LOCATE_CONTACT NNL] i=%d its=%d t=%.15G d=%.15G\n", i, its, t+t1, d);
#endif
      //normddot = calcvecF(i, j, t, r1, r2, ddot, shift);
#ifdef MD_BASIC_DT
      delt = epsd/maxddot;
      tini = t;
      t += delt;
      d = calcDistNegNeighPlaneAll_sp(NSP, t, t1, i, dists);
#else
      if (!firstaftsf)
	{
	  deldist = get_max_deldist_sp(NSP, distsOld2, distsOld);
	  normddot = fabs(deldist)/delt;
	  /* NOTA: forse qui si potrebbe anche usare sempre delt = epsd/maxddot */
	  if (normddot!=0)
	    delt = epsd/normddot;
	  else
	    delt = epsd/maxddot;
	    //delt = h;
	  if (fabs(dold) < epsd)
	    delt = epsd / maxddot;
	}
      else
	{
	  delt = h;//EPS*fabs(t);
	  firstaftsf = 0;
	  dold2 = calcDistNegNeighPlaneAll_sp(NSP, t-delt, t1, i, distsOld2);
	  continue;
	}
      tini = t;
      t += delt;
      d = calcDistNegNeighPlaneAll_sp(NSP, t, t1, i, dists);
      MD_DEBUG30(printf("t=%.15G d=%.15G\n", t, d));
      deldist = get_max_deldist_sp(NSP, distsOld, dists);
      if (deldist > epsdMax)
	{
	  /* se la variazione di d è eccessiva 
	   * cerca di correggere il passo per ottenere un valore
	   * più vicino a epsd*/
	  t -= delt;
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
	  t += delt; 
	  itsSNL++;
	  d = calcDistNegNeighPlaneAll_sp(NSP, t, t1, i, dists);
	  //assign_vec(vecgdold2, vecgd);
	}
#endif
      for (nn=0; nn < 6; nn++)
	for (nn2=0; nn2 < NSP; nn2++)
	  dorefine[nn][nn2] = 0;
      ncr=check_cross_sp(NSP, distsOld, dists, crossed);
      for (nn = 0; nn < 6; nn++)
	{
	  for (nn2 = 0; nn2 < NSP; nn2++)
	    {
	      t2arr[nn][nn2] = t; 
	      dorefine[nn][nn2] = crossed[nn][nn2];
#if 0
	      if (crossed[nn]!=0)
		{
		  if (distsOld[nn] > 0 && dists[nn] < 0)
		    {
		      dorefine[nn] = 1;
		    }
		}
#endif
	    }
	}
#if 1
      ntc = get_dists_tocheck_sp(NSP, distsOld, dists, tocheck, dorefine);
      for (nn = 0; nn < 6; nn++)
	{
	  for (nn2 = 0; nn2 < NSP; nn2++)
	    {
	      if (tocheck[nn][nn2])
		{
		  if (interpolNeighPlane_sp(i, t1, t-delt, delt, distsOld[nn][nn2], dists[nn][nn2], 
					    &troot, nn, nn2))
		    {
		      dorefine[nn][nn2] = 0;
		    }
		  else 
		    {
		      MD_DEBUG32(printf("qui-1 t-delt=%.15G t=%.15G t2arr=%.15G\n", t-delt,t, troot));
		      dorefine[nn][nn2] = 1;
		      t2arr[nn][nn2] = troot-t1;
		    }
		}
	      else if (dorefine[nn][nn2])
		{
		  t2arr[nn][nn2] = t;
		}
	    }
	}
#endif
      gotcoll = 0;
      firstev = 1;
      for (nn = 0; nn < 6; nn++)
	{
	  for (nn2 = 0; nn2 < NSP; nn2++)
	    {
	      if (dorefine[nn][nn2]!=0)
		{
		  assign_plane(nn);
		  MD_DEBUG32(printf("t1=%.15G t2=%.15G delt=%.15G dists[%d][%d]: %.15G distsOld[%d][%d]:%.15G\n", t1+t-delt, t1+t2arr[nn][nn2], delt, nn, nn2, dists[nn][nn2], nn, nn2, distsOld[nn][nn2]));
		  MD_DEBUG32(printf("d(%.15G)=%.15G d(%.15G)=%.15G\n",
				    t-delt, calcDistNegOneNNL_sp(t-delt, t1, i, nn2),
				    t2arr[nn][nn2], calcDistNegOneNNL_sp(t2arr[nn][nn2],
									 t1, i, nn2)));
		  if (refine_contact_neigh_plane_sp(i, t1, t-delt, t2arr[nn][nn2], &troot, nn, nn2))
		    {
		      //printf("[locate_contact] Adding collision for ellips. N. %d t=%.15G t1=%.15G t2=%.15G\n", i,
		      //	 vecg[4], t1 , t2);
		      MD_DEBUG(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
		      MD_DEBUG(printf("[locate_contact] t=%.15G nn=%d\n", t, nn));
		      MD_DEBUG(printf("[locate_contact] its: %d\n", its));
		      /* se il legame già c'è e con l'urto si forma tale legame allora
		       * scarta tale urto */
		      if (troot > t2 || troot < t1)
			{
			  continue;
			}
		      else
			{
			  gotcoll = 1;

			  if (firstev || troot < *evtime)
			    {
			      firstev = 0;
			      *evtime = troot;
			    }
			  //printf("QUI\n");
			  if (nn==5)
			    return 1;
			  else
			    continue;
			}
		    }
		  else 
		    {
		      MD_DEBUG(printf("[locate_contact_sp] can't find contact point!\n"));
#ifdef MD_INTERPOL
		      if (!tocheck[nn][nn2])
#endif
			mdPrintf(ALL,"[locate_contact_nnl_sp] can't find contact point!\n",NULL);

		      MD_DEBUG32(printf("t1=%.15G t2=%.15G delt=%.15G dists[%d][%d]: %.15G distsOld[%d][%d]:%.15G\n", t-delt, t2arr[nn][nn2], delt, nn, nn2, dists[nn][nn2], nn, nn2, distsOld[nn][nn2]));
		      /* Se refine_contact fallisce deve cmq continuare a cercare 
		       * non ha senso smettere...almeno credo */
		      //gotcoll = -1;
		      continue;
		    }
		}
	    }
	}
      if (gotcoll == 1)
	return 1;
      if (fabs(d) > epsdFastR)
	{
	  if (search_contact_faster_neigh_plane_all_sp(i, &t, t1, t2, epsd, &d, epsdFast, 
						       dists, maxddoti, maxddot))
	    {
	      MD_DEBUG30(printf("[search contact faster locate_contact] d: %.15G\n", d));
	      return 0;
	    }
	  dold = d;
	  assign_dists_sp(NSP, dists, distsOld);
	  firstaftsf = 1;
	  its++;
	  continue;
	}
      dold = d;
#ifndef MD_BASIC_DT
      assign_dists_sp(NSP, distsOld,  distsOld2);
#endif
      assign_dists_sp(NSP, dists, distsOld);
      its++;
      itsSNL++;
    }
  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return 0;
}
#endif
