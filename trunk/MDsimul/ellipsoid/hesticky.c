#include<mdsimul.h>
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG(x) 
#define MD_DEBUG10(x)  
#define MD_DEBUG11(x) 
#define MD_DEBUG15(x) 
#define MD_DEBUG20(x) 
#define MD_DEBUG29(x) 
#define MD_DEBUG30(x) 
#define MD_NEGPAIRS
#define MD_NO_STRICT_CHECK
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
#ifdef MD_ASYM_ITENS
double **Ia, **Ib, **invIa, **invIb;
#else
double Ia, Ib, invIa, invIb;
#endif
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
/* *** change here if you change the number sticky spots *** */
int mapbondsaAA[MD_PBONDS]={1,1,2,2,3,3,4,4};
int mapbondsbAA[MD_PBONDS]={1,2,1,2,1,2,1,2};
int mapbondsaAB[MD_PBONDS]={1,1,2,2,3,3,4,4};
int mapbondsbAB[MD_PBONDS]={1,2,1,2,1,2,1,2};
int mapbondsaBB[MD_PBONDS]={1,1,2,2,3,3,4,4};
int mapbondsbBB[MD_PBONDS]={1,2,1,2,1,2,1,2};
int *mapbondsa;
int *mapbondsb;
/* ------------------------------------------------------------ */
extern int *crossevtodel;
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
double *lastcol;
double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC � la coordinata del punto di contatto */
int  **tree, cellRange[2*NDIM], initUcellx, initUcelly, initUcellz;
int *inCell[3], *cellList, cellsx, cellsy, cellsz;
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
double calc_norm(double *vec);
double rcutL, aL, bL, cL;
extern double max_ax(int i);
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[]);

void bumpSP(int i, int j, int ata, int atb, double* W, int bt)
{
  /* NOTA: Controllare che inizializzare factor a 0 � corretto! */
  double factor=0, invmi, invmj, sigmai, mredl, shift[3];
  double delpx, delpy, delpz, wrx, wry, wrz, rACn[3], rBCn[3], r12[3];
  double rAB[3], rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double ratA[3], ratB[3], norm[3];
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, Txyold, Tyzold, Tzxold, taus, 
	  DTxx, DTyy, DTzz, Txxold, Tyyold, Tzzold;
#endif

  double denom, rCx, rCy, rCz, nrAB, Dr, V, Eold, E, Vold, Kold;
#ifndef MD_ASYM_ITENS
  double factorinvIa, factorinvIb;
#endif
  int na, a, kk, oldbond=-1;
#if 1
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
  if (bt == MD_CORE_BARRIER)
    {
      bump(i, j, W, Sqr(sigmai));
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
  RDiagtR(i, Ia, ItensD[na][0], ItensD[na][1], ItensD[na][2], R[i]);
#else
  Ia = Oparams.I[na];
#endif
  na = (j < Oparams.parnumA)?0:1;
  if (OprogStatus.targetPhi > 0)
    {
      /* scalare tutti i raggi qui */
    }
#ifdef MD_ASYM_ITENS
  RDiagtR(j, Ib, ItensD[na][0], ItensD[na][1], ItensD[na][2], R[j]);
#else
  Ib = Oparams.I[na];
#endif
  MD_DEBUG(check_contact(evIdA, evIdB, Xa, Xb, rAC, rBC));
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
  /* calcola le velocit� nel punto di contatto */
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
       * Notare che Oparams.bheight � la profondit� della buca ed 
       * � una quantit� positiva!!*/
      /* b = vc*r ma ora b -> vc riscrivere correttamente le seguenti equazioni!! */
    case MD_CORE_BARRIER:
      factor = -2.0*vc;
      factor *= mredl;
      break;
    case MD_INOUT_BARRIER:
      if (Sqr(vc) < 2.0*Oparams.bheight/mredl)
	{
	  MD_DEBUG30(printf("t=%.15G vc=%.15G NOT ESCAPEING collType: %d d=%.15G\n", Oparams.time,
		 vc,  bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
	  factor = -2.0*vc;
	}
      else
	{
	  MD_DEBUG(printf("t=%.15G vc=%.15G ESCAPING collType: %d d=%.15G\n", Oparams.time, vc, bt,
		 sqrt(Sqr(ratA[0]-ratB[0])+Sqr(ratA[1]-ratB[1])+Sqr(ratA[2]-ratB[2]))));
	  factor = -vc + sqrt(Sqr(vc) - 2.0*Oparams.bheight/mredl);
	  remove_bond(i, j, ata, atb);
	  remove_bond(j, i, atb, ata);
	}
      factor *= mredl;
      break;
    case MD_OUTIN_BARRIER:
      add_bond(i, j, ata, atb);
      add_bond(j, i, atb, ata);
      factor = -vc - sqrt(Sqr(vc) + 2.0*Oparams.bheight/mredl);
      MD_DEBUG(printf("delta= %f height: %f mredl=%f\n", 
		      Sqr(vc) + 2.0*Oparams.bheight/mredl, Oparams.bheight, mredl));
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
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, vx[i],vy[i],vz[i], j, vx[j],vy[j],vz[j]));
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, wx[i],wy[i],wz[i], j, wx[j],wy[j],wz[j]));
}
void check_bonds(char* msg, int i, int j, int ata, int atb, int yesexit)
{
  int a, b, B1, B2;
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
   * dove b � l'atomo di j */
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
      printf("il bond %d,%d eiste gi�!\n", na, n);
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
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3])
{
  /* QUESTA VA RISCRITTA PER GLI ELLISSOIDI STICKY!!! */
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  double r1[3], r2[3], r3[3], nr, fact;
  double radius; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */

  radius = Oparams.sigma[0][1] / 2.0;
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
      /* l'atomo restante � un electron site */
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
void BuildAtomPos(int i, double *rO, double **R, double rat[5][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a;
  /* l'atomo zero si suppone nell'origine */
  for (a=0; a < 5; a++)
    BuildAtomPosAt(i, a, rO, R, rat[a]);
}

void assign_bond_mapping(int i, int j)
{
  /* NOTA: l'interazione bonded � solo tra Si e O 
   * i <  Oparams.parnumA => O
   * i >=  Oparams.parnumA => Si */
  if (i < Oparams.parnumA && j < Oparams.parnumA)
    {
      mapbondsa = mapbondsbAA;
      mapbondsb = mapbondsaAA;
    }
  else if (i > Oparams.parnumA && j > Oparams.parnumA) 
    {
      mapbondsa = mapbondsaBB;
      mapbondsb = mapbondsbBB;
    }
  else if (i < Oparams.parnumA && j > Oparams.parnumA)
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
int ibr, jbr, nnbr; 
double shiftbr[3], trefbr;

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
  double Omega[3][3];
  int na;
  assign_bond_mapping(i, j);
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
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, RtB, Omega, mapbondsb[nn]);
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
  double Omega[3][3];
  int na;
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
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
  ti = t + (t1 - atomTime[j]);
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
 
  UpdateOrient(j, ti, RtB, Omega, (bondpair==-1)?-1:mapbondsb[bondpair]);
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
	  fabs(dists[nn]) < OprogStatus.epsd && fabs(distsOld[nn]) < OprogStatus.epsd &&
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
#define MD_OPTDDIST
/* NOTA: tale stima ottimizzata della maggiorazione per la velocit� di variazione della distanza
 * sembra corretta, fare comunque dei test.*/
double eval_maxddistSP(int i, int j, int bondpair, double t1, double *maxddotOpt)
{
  double ti, rA[3], rB[3], Omega[3][3], ratA[NA][3], ratB[NA][3], wri[3], wrj[3], nwri, nwrj,
	 r12i[3], r12j[3];//, maxddotOpt[MD_PBONDS];
  double maxddot=0.0, nr12i, nr12j;
  int nn, kk;
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
   * MAXOPTITS � il numero massimo di iterazioni al di sopra del quale esce */
  double told, delt, distsOld[MD_PBONDS];
  const int MAXOPTITS = 500;
  int nn, its=0, amin, bmin, crossed[MD_PBONDS]; 
  /* estimate of maximum rate of change for d */
#if 0
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]*0.5
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j]*0.5;
#endif
  *d1 = calcDistNegSP(*t, t1, i, j, shift, &amin, &bmin, distsOld, bondpair);
  MD_DEBUG30(printf("[IN SEARCH CONTACT FASTER]*d1=%.15G t=%.15G\n", *d1, *t));
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
      if (maxddot*(t2-(t1+*t)) < fabs(*d1)-OprogStatus.epsd)
	return 1;
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
	  MD_DEBUG30(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNeg(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 0;
	}
      if (*t+t1 > t2)
	{
	  *t = told;
	  MD_DEBUG30(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("convergence t>t2\n"));
	  *d1 = calcDistNegSP(*t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	  return 1;
	}
      told = *t;
      assign_distsSP(dists, distsOld);
      its++;
      itsF++;
    }

  MD_DEBUG10(printf("max iterations %d iterations reached t=%f t2=%f\n", its, *t, t2));
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
int interpolSP(int i, int j, int nn, 
	     double tref, double t, double delt, double d1, double d2,
	     double *tmin, double shift[3])
{
  int its, nb, amin, bmin;
  double d3, A, B, C, dmin, dminnew, tminnew;
  double dists[MD_PBONDS];
  /* NOTA: dists di seguito pu� non essere usata? controllare!*/
  d3 = calcDistNegOne(t+delt*0.5, tref, i, j, nn, shift);
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
  dmin = calcDistNegOneSP(*tmin, tref, i, j, nn, shift);
  if (*tmin < t+delt && *tmin > t && d1*dmin < 0.0)
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
      sum += 1;
    }
  return (sum > 0)?1:0;
}

int delt_is_too_big_hc(int i, int j, int bondpair, double *dists, double *distsOld)
{
  int nn, retval=0;
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
  int nn, retval=0;
  for (nn=0; nn < MD_PBONDS; nn++)
    {
      if (bondpair != -1 && bondpair != nn)
	continue;
      if (!negpairs[nn])
	continue;
    if (dists[nn] > 0.0 && bound(i,j,mapbondsa[nn],mapbondsb[nn]))
      return 1;
    if (dists[nn] < 0.0 && !bound(i,j,mapbondsa[nn],mapbondsb[nn]))
      return 1;
    }
  return 0;
}
int locate_contactSP(int i, int j, double shift[3], double t1, double t2, 
		   double *evtime, int *ata, int *atb, int *collCode)
{
  const double minh = 1E-14;
  double h, d, dold, dold2, t2arr[MD_PBONDS], t, dists[MD_PBONDS], distsOld[MD_PBONDS],
	 distsOld2[MD_PBONDS], deltth; 
  double normddot, maxddot, delt, troot, tmin, tini; //distsOld2[MD_PBONDS];
  //const int MAXOPTITS = 4;
  int bondpair, itstb;
  const int MAXITS = 100;
  int itsRef;
  int its, foundrc, goback;
  double t1ini, delthc;
  double maxddoti[MD_PBONDS], epsd, epsdFast, epsdFastR, epsdMax, deldist, df; 
  int kk,tocheck[MD_PBONDS], dorefine[MD_PBONDS], ntc, ncr, nn, gotcoll, amin, bmin,
      crossed[MD_PBONDS], firstaftsf;
#ifdef MD_NEGPAIRS
  int negpairs[MD_PBONDS], sumnegpairs;
#endif
  const double GOLD= 1.618034;
  epsd = OprogStatus.epsd;
  epsdFast = OprogStatus.epsdFast;
  epsdFastR= OprogStatus.epsdFastR;
  epsdMax = OprogStatus.epsdMax;
  assign_bond_mapping(i, j);
  bondpair = get_bonded(i, j);
#ifdef MD_OPTDDIST
  maxddot = eval_maxddistSP(i, j, bondpair, t1, maxddoti);
#else
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i]*0.5
    + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j]*0.5;
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
   * 1) andare avanti nel tempo finch� la distanza non � corretta.
   * 2) fare un passo ed eventualmente ridurlo finch� la distanza non � corretta.
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
  MD_DEBUG30(printf("[BEFORE SEARCH CONTACT FASTER]Dopo distances between %d-%d t=%.15G t2=%.15G\n", i, j, t, t2));
#ifdef MD_NEGPAIRS
  sumnegpairs = check_negpairs(negpairs, bondpair, i, j); 
#endif
  if (search_contact_fasterSP(i, j, shift, &t, t1, t2, epsd, &d, epsdFast, dists, bondpair, maxddot, maxddoti))
    {
      return 0;  
    }
  timesS++;
  MD_DEBUG30(printf("[AFTER SEARCH CONTACT FASTER]Dopo distances between %d-%d d1=%.12G\n", i, j, d));

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
	      delt = h;
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
	  MD_DEBUG30(printf("==========>>>>> t=%.15G t2=%.15G\n", t, t2));
	  continue;
	}
      MD_DEBUG30(printf("delt: %f epsd/maxddot:%f h*t:%f maxddot:%f\n", delt, epsd/maxddot,h*t,maxddot));
      tini = t;
      t += delt;
      d = calcDistNegSP(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);

      deldist = get_max_deldistSP(distsOld, dists, bondpair);
      if (deldist > epsdMax)
	{
	  /* se la variazione di d � eccessiva 
	   * cerca di correggere il passo per ottenere un valore
	   * pi� vicino a epsd*/
	  t -= delt;
	  //delt = d2old / maxddot;
	  delt = epsd/maxddot;
	  /* NOTE: prob. la seguente condizione si puo' rimuovere 
	   * o cambiare in > */
	  deltth = h;
	  if (delt < deltth)
	    {
	      delt = deltth;
	    }
	  t += delt; 
	  itsS++;
	  d = calcDistNeg(t, t1, i, j, shift, &amin, &bmin, dists, bondpair);
	}
#ifdef MD_NEGPAIRS
      itstb = 0;
      /* NOTA: se la distanza tra due sticky spheres � positiva a t (per errori numerici 
       * accade spesso) e t+delt allora delt � troppo grande e qui lo riduce fino ad un 
       * valore accettabile. */
      if (sumnegpairs)// && !firstaftsf)
	{
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
	      /* se dorefine � 2 vuol dire che due superfici si sono
	       * attraversate */
	      if (valid_collision(i, j, mapbondsa[nn], mapbondsb[nn], crossed[nn]))
		{
		  MD_DEBUG30(printf("type: %d i=%d j=%d ata=%d atb=%d bound:%d\n", crossed[nn], i, j, mapbondsa[nn],
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
			   &troot, shift))
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
		    t2arr[nn] = troot;
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
		  MD_DEBUG30(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
		  MD_DEBUG30(printf("[locate_contact] t=%.15G nn=%d\n", t, nn));
		  MD_DEBUG(printf("[locate_contact] its: %d\n", its));
		  /* se il legame gi� c'� e con l'urto si forma tale legame allora
		   * scarta tale urto */
		  if (troot > t2 || troot < t1
#if 1
		  || 
		      (lastbump[i].mol == j && lastbump[j].mol==i && 
		       lastbump[i].at == mapbondsa[nn]
		       && lastbump[j].at == mapbondsb[nn] && fabs(troot - lastcol[i]) < 1E-15))
#endif
		    {
		      //gotcoll = -1;
		      continue;
		    }
		  else
		    {
		      gotcoll = 1;
		      if (*collCode == MD_EVENT_NONE || troot < *evtime)
			{
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
	      MD_DEBUG30(printf("[search contact faster locate_contact] d: %.15G\n", d));
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
      MD_DEBUG30(printf("==========>>>>> t=%.15G t2=%.15G\n", t, t2));
      assign_distsSP(distsOld,  distsOld2);
      assign_distsSP(dists, distsOld);
      its++;
      itsS++;
    }
  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return 0;
}


