#include<mdsimul.h>
#define MD_NNLPLANES
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG10(x)  
#define MD_DEBUG11(x) 
#define MD_DEBUG15(x) 
#define MD_DEBUG20(x) 
#define MD_DEBUG30(x)
#define MD_DEBUG31(x) 
#define MD_DEBUG32(x) 
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern const double timbig;
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
extern double rA[3], rB[3];
extern double rxC, ryC, rzC;
extern double pi, invL, L2, Vz; 
#ifdef MD_ASYM_ITENS
extern double **Ia, **Ib, **invIa, **invIb;
#else
extern double Ia, Ib, invIa, invIb;
#endif
extern double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
extern int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
  cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
extern int evIdA, evIdB, parnumB, parnumA;
extern double *axa, *axb, *axc;
extern int *scdone;
extern double *maxax;
extern double xa[3], ya[3];
extern int polinterr, polinterrRyck;
extern double *lastupdNNL, *totDistDispl;
double gradplane[3];
/* Routines for LU decomposition from Numerical Recipe online */
#ifdef MD_ASYM_ITENS
extern double *phi0, *psi0, *costheta0, *sintheta0, **REt, **RE0, *angM, ***RM, **REtA, **REtBi, **Rdot, **REtA, **REtB;
extern double cosEulAng[2][3], sinEulAng[2][3];
void symtop_evolve_orient(int i, double ti, double **Ro, double **REt, double cosea[3], double sinea[3], double *phir, double *psir);
#endif

void ludcmpR(double **a, int* indx, double* d, int n);
void lubksbR(double **a, int* indx, double *b, int n);
void InvMatrix(double **a, double **b, int NB);
extern double invaSq[2], invbSq[2], invcSq[2];
extern double rxC, ryC, rzC, trefG;
extern int SolveLineq (double **a, double *x, int n); 
extern int calcdist_retcheck;
void comvel_brown (COORD_TYPE temp, COORD_TYPE *m);
void InitEventList (void);
void writeAsciiPars(FILE* fs, struct pascii strutt[]);
void writeAllCor(FILE* fs);
extern struct nebrTabStruct *nebrTab;
double nextNNLrebuild;
extern void UpdateSystem(void);
extern void UpdateOrient(int i, double ti, double **Ro, double Omega[3][3]);
extern void nextNNLupdate(int na);
extern void BuildNNL(int na);
extern long long int itsF, timesF, itsS, timesS, numcoll, itsFNL, itsSNL, timesSNL, timesFNL;
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern double accngA, accngB;
extern void tRDiagR(int i, double **M, double a, double b, double c, double **Ri);
extern double min(double a, double b);
#ifdef MD_ASYM_ITENS
extern void calcFxtFt(int i, double x[3], double **RM, double cosea[3], double sinea[3], double **X,
	       double D[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft);
#else
extern void calcFxtFt(double x[3], double **X,
	       double D[3][3], double Omega[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft);
#endif
double calcDistNegNNLoverlapPlane(double t, double t1, int i, int j, double shift[3]);
extern double calc_norm(double *vec);
extern void funcs2beZeroedDistNeg(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroedDistNeg5(int n, double x[], double fvec[], int i, int j, double shift[3]);

extern void calc_intersec_neigh(double *rB, double *rA, double **Xa, double* rI, double alpha);
extern double scalProd(double *A, double *B);
extern double distfunc(double x);
extern double min3(double a, double b, double c);
extern double max3(double a, double b, double c);
extern void distSD(int i, int j, double shift[3], double *vecg, double lambda, int halfspring);
extern void calc_grad(double *rC, double *rA, double **Xa, double *grad);
extern void calc_intersec(double *rB, double *rA, double **Xa, double* rI);
extern double min(double a, double b);
extern void newtDistNegNeigh(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int),
	  int iA);
extern void newtDistNeg(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3]);
extern void store_bump(int i, int j);
extern void store_bump_neigh(int i, double *r1, double *r2);
extern double calcDist(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
extern void print_matrix(double **M, int n);
extern void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[], double xb2[], 
	   int *nb);
extern double zbrent(double (*func)(double), double x1, double x2, double tol);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void guess_dist(int i, int j, 
		double *rA, double *rB, double **Xa, double **Xb, double *rC, double *rD,
		double **RA, double **RB);
extern int search_contact_faster_neigh(int i, double *t, double t1, double t2, 
				double *vecgd, double epsd, double *d1, double epsdFast,
				double *r1, double *r2);
extern double max_ax(int i);
extern int locate_contact(int i, int j, double shift[3], double t1, double t2, double vecg[5]);
extern void ScheduleEvent (int idA, int idB, double tEvent);
extern int refine_contact_neigh(int i, double t1, double t, double vecgd[8], double  vecg[5]);
extern void newtNeigh(double x[], int n, int *check, 
  void (*vecfunc)(int, double [], double [], int),
  int iA);
extern void newtDistNegNeighPlane(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int),
	  int iA);
extern double max(double a, double b);
#ifdef MD_ASYM_ITENS
double calc_maxddot_nnl(int i, double *gradplane)
{
#if 0
  int na;
  double Iamin;
#endif
  double factori;
  factori = 0.5*maxax[i]+OprogStatus.epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
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
#endif

void calc_grad_and_point_plane(int i, double *grad, double *point, int nplane)
{
  int kk;
  double del=0.0, segno;
  for (kk=0; kk < 3; kk++)
    {
      /* NOTA: controllare che non si debbano scambiare kk e nplane/2 */ 
      grad[kk] = nebrTab[i].R[nplane/2][kk];
    }
  switch (nplane/2)
    {
    case 0:
      del = nebrTab[i].axa;	
      break;
    case 1:
      del = nebrTab[i].axb;	
      break;
    case 2:
      del = nebrTab[i].axc;	
      break;
    }

  /* NOTA: epsdNL+epsd viene usato come buffer per evitare problemi numerici 
   * nell'update delle NNL. */
  del -= OprogStatus.epsdNL+OprogStatus.epsd;
  for (kk=0; kk < 3; kk++)
    {
      if (nplane % 2 == 0)
	segno = 1;
      else
	segno = -1;
      grad[kk] *= segno;
      /* rB[] (i.e. nebrTab[i].r[]) è un punto appartenente al piano */
      point[kk] = nebrTab[i].r[kk] + del*grad[kk]; 
    }
}
void growth_rebuildNNL(int i)
{
  int ii, n;
  nextNNLupdate(i);
  BuildNNL(i);
  if (nebrTab[i].nexttime < nextNNLrebuild)
    nextNNLrebuild = nebrTab[i].nexttime;
  for (ii=0; ii < nebrTab[i].len; ii++)
    {
      n = nebrTab[i].list[ii]; 
      BuildNNL(n);
    }
}

void rebuildNNL(void)
{
  int i;
  double nltime=timbig;
  UpdateSystem();
  if (OprogStatus.useNNL==2)
    printf("Rebuilding NNL t=%.15G\n", Oparams.time);
  for (i=0; i < Oparams.parnum; i++)
    {
      nextNNLupdate(i);
      if (i==0 || nebrTab[i].nexttime < nltime)
	nltime = nebrTab[i].nexttime;
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      BuildNNL(i);
    }
  /* next complete update */
  nextNNLrebuild = nltime;
  if (OprogStatus.useNNL==2)
    printf("nextNNLrebuild=%.15G\n", nextNNLrebuild);
  //ScheduleEvent(-1, ATOM_LIMIT + 11, nltime); 
}

void fdjacNeighPlane(int n, double x[], double fvec[], double **df, 
		     void (*vecfunc)(int, double [], double [], int), int iA)
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  double  rA[3], ti, vA[3], vB[3], OmegaB[3][3];
  double DA[3][3], fx[3], invaSqN, invbSqN, invcSqN;
  double Fxt[3], Ft;
#ifdef MD_ASYM_ITENS
  double psi, phi;
#else
  double OmegaA[3][3];
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
  invaSqN = 1/Sqr(axa[iA]);
  invbSqN = 1/Sqr(axb[iA]);
  invcSqN = 1/Sqr(axc[iA]);

  tRDiagR(iA, Xa, invaSqN, invbSqN, invcSqN, RA);
  MD_DEBUG2(printf("invabc: (%f,%f,%f)\n", invaSq, invbSq, invcSq));
  MD_DEBUG2(print_matrix(Xa, 3));
  DA[0][1] = DA[0][2] = DA[1][0] = DA[1][2] = DA[2][0] = DA[2][1] = 0.0;
  DA[0][0] = invaSqN;
  DA[1][1] = invbSqN;
  DA[2][2] = invcSqN;

  /*N.B. l'ellissoide B in tale caso non evolve! */
  //rB[0] = nebrTab[iA].r[0];
  //rB[1] = nebrTab[iA].r[1];
  //rB[2] = nebrTab[iA].r[2];

  vB[0] = 0.0;
  vB[1] = 0.0;
  vB[2] = 0.0;
  OmegaB[0][0] = 0;
  OmegaB[0][1] = 0;
  OmegaB[0][2] = 0;
  OmegaB[1][0] = 0;
  OmegaB[1][1] = 0;
  OmegaB[1][2] = 0;
  OmegaB[2][0] = 0;
  OmegaB[2][1] = 0;
  OmegaB[2][2] = 0;
   
#if 0
  invaSqN = 1.0/Sqr(nebrTab[iA].axa);
  invbSqN = 1.0/Sqr(nebrTab[iA].axb);
  invcSqN = 1.0/Sqr(nebrTab[iA].axc);
  tRDiagR(iA, Xb, invaSqN, invbSqN, invcSqN, nebrTab[iA].R);
  DB[0][1] = DB[0][2] = DB[1][0] = DB[1][2] = DB[2][0] = DB[2][1] = 0.0;
  DB[0][0] = invaSqN;
  DB[1][1] = invbSqN;
  DB[2][2] = invcSqN;
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];// - Sqr(x[3])*Xb[k1][k2]);
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      //gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  //gx[k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]);
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
      df[4][k1] = gradplane[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[k1][3] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[k1][3] += 4.0*x[3]*Xb[k1][k2]*(x[k2]-rB[k2]); 
#endif
      df[k1][3] = -2.0*x[3]*gradplane[k1];
    } 
  df[3][3] = 0.0;
  df[4][3] = 0.0;
#ifdef MD_ASYM_ITENS
  calcFxtFt(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fx, Fxt, &Ft);
#else
  calcFxtFt(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
#endif
  //calcFxtFt(x, Xb, DB, OmegaB, RB, rB, vB, gx, Gxt, &Gt);
  for (k1 = 0; k1 < 3; k1++)
    {
      //df[k1][4] = 0;
      //for (k2 = 0; k2 < 3; k2++)
      df[k1][4] = Fxt[k1];//-Sqr(x[3])*Gxt[k1]; 
    } 
 df[3][4] = Ft;
 df[4][4] = 0.0;//Gt;
#ifndef MD_GLOBALNRNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1]-rB[k1])*gradplane[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 //fvec[4] = 0.5*fvec[4]-1.0;
 MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
#endif
}
/* funzione che calcola lo Jacobiano */
/* N.B. In tale caso l'ellissoide B è semplicemente l'ellissoide A riscalato di un opportuno 
 * fattore */
void fdjacNeigh(int n, double x[], double fvec[], double **df, 
		void (*vecfunc)(int, double [], double [], int), int iA)
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  double  rA[3], rB[3], ti, vA[3], vB[3], OmegaB[3][3];
  double DA[3][3], DB[3][3], fx[3], gx[3], invaSqN, invbSqN, invcSqN;
  double Fxt[3], Ft;
  int k1, k2;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#else
  double OmegaA[3][3];
#endif
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
  invaSqN = 1/Sqr(axa[iA]);
  invbSqN = 1/Sqr(axb[iA]);
  invcSqN = 1/Sqr(axc[iA]);

  tRDiagR(iA, Xa, invaSqN, invbSqN, invcSqN, RA);
  MD_DEBUG2(printf("invabc: (%f,%f,%f)\n", invaSq, invbSq, invcSq));
  MD_DEBUG2(print_matrix(Xa, 3));
  DA[0][1] = DA[0][2] = DA[1][0] = DA[1][2] = DA[2][0] = DA[2][1] = 0.0;
  DA[0][0] = invaSqN;
  DA[1][1] = invbSqN;
  DA[2][2] = invcSqN;

  /*N.B. l'ellissoide B in tale caso non evolve! */
  rB[0] = nebrTab[iA].r[0];
  rB[1] = nebrTab[iA].r[1];
  rB[2] = nebrTab[iA].r[2];
  vB[0] = 0.0;
  vB[1] = 0.0;
  vB[2] = 0.0;
  OmegaB[0][0] = 0;
  OmegaB[0][1] = 0;
  OmegaB[0][2] = 0;
  OmegaB[1][0] = 0;
  OmegaB[1][1] = 0;
  OmegaB[1][2] = 0;
  OmegaB[2][0] = 0;
  OmegaB[2][1] = 0;
  OmegaB[2][2] = 0;
   
  invaSqN = 1.0/Sqr(nebrTab[iA].axa);
  invbSqN = 1.0/Sqr(nebrTab[iA].axb);
  invcSqN = 1.0/Sqr(nebrTab[iA].axc);

  tRDiagR(iA, Xb, invaSqN, invbSqN, invcSqN, nebrTab[iA].R);
  DB[0][1] = DB[0][2] = DB[1][0] = DB[1][2] = DB[2][0] = DB[2][1] = 0.0;
  DB[0][0] = invaSqN;
  DB[1][1] = invbSqN;
  DB[2][2] = invcSqN;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*(Xa[k1][k2] - Sqr(x[3])*Xb[k1][k2]);
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
      df[k1][3] = -2.0*x[3]*gx[k1];
    } 
  df[3][3] = 0.0;
  df[4][3] = 0.0;
#ifdef MD_ASYM_ITENS
  calcFxtFt(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fx, Fxt, &Ft);
#else
  calcFxtFt(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
#endif
  //calcFxtFt(x, Xb, DB, OmegaB, RB, rB, vB, gx, Gxt, &Gt);
  for (k1 = 0; k1 < 3; k1++)
    {
      //df[k1][4] = 0;
      //for (k2 = 0; k2 < 3; k2++)
      df[k1][4] = Fxt[k1];//-Sqr(x[3])*Gxt[k1]; 
    } 
 df[3][4] = Ft;
 df[4][4] = 0.0;//Gt;
#ifndef MD_GLOBALNRNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gx[k1];
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
void funcs2beZeroedNeighPlane(int n, double x[], double fvec[], int i)
{
  int na, k1, k2; 
  double  rA[3], ti;
  double fx[3];
  double invaSqN, invbSqN, invcSqN;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#else
  double Omega[3][3];
#endif
  /* x = (r, alpha, t) */ 
  ti = x[4] + (trefG - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, Rt, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, Rt, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
  tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, Rt);

  /* il secondo ellissoide resta fermo al tempo iniziale */
  MD_DEBUG(printf("x[4]:%.15f atomTime[%d]:%.15f\n",x[4], j, atomTime[j]));
  /* il secondo ellissoide non deve evolvere nel tempo */
#if 0
  rB[0] = nebrTab[i].r[0];
  rB[1] = nebrTab[i].r[1];
  rB[2] = nebrTab[i].r[2];
  invaSqN = 1.0/Sqr(nebrTab[i].axa);
  invbSqN = 1.0/Sqr(nebrTab[i].axb);
  invcSqN = 1.0/Sqr(nebrTab[i].axc);
  tRDiagR(i, Xb, invaSqN, invbSqN, invcSqN, nebrTab[i].R);
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
      	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
      MD_DEBUG(printf("[FUNC2BEZ]x[%d]:%.15f rA[%d]:%f fx:%.15f\n", k1, x[k1], k1, rA[k1],fx[k1]));
      
   }
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2] - rB[k2]);
      MD_DEBUG(printf("[FUNC2BEZ]x[%d]:%.15f rB[%d]:%f gx:%.15f\n", k1, x[k1], k1, rB[k1],gx[k1]));
    }
#endif
  MD_DEBUG(print_matrix(Xb,3));
  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      fvec[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fvec[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]) + 2.0*Sqr(x[3])*Xb[k1][k2]*(x[k2] - rB[k2]);
#endif
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
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
      fvec[4] += (x[k1]-rB[k1])*gradplane[k1];
#endif
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  //fvec[4] = 0.5*fvec[4]-1.0;
  MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
}
void funcs2beZeroedNeigh(int n, double x[], double fvec[], int i)
{
  int na, k1, k2; 
  double  rA[3], rB[3], ti;
  double fx[3], gx[3];
  double invaSqN, invbSqN, invcSqN;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#else
  double Omega[3][3];
#endif
  /* x = (r, alpha, t) */ 
  ti = x[4] + (trefG - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, Rt, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, Rt, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
  tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, Rt);

  /* il secondo ellissoide resta fermo al tempo iniziale */
  MD_DEBUG(printf("x[4]:%.15f atomTime[%d]:%.15f\n",x[4], j, atomTime[j]));
  /* il secondo ellissoide non deve evolvere nel tempo */
  rB[0] = nebrTab[i].r[0];
  rB[1] = nebrTab[i].r[1];
  rB[2] = nebrTab[i].r[2];
  invaSqN = 1.0/Sqr(nebrTab[i].axa);
  invbSqN = 1.0/Sqr(nebrTab[i].axb);
  invcSqN = 1.0/Sqr(nebrTab[i].axc);
  tRDiagR(i, Xb, invaSqN, invbSqN, invcSqN, nebrTab[i].R);
  
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
      fvec[k1] = fx[k1] - Sqr(x[3])*gx[k1];
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
extern double gradplane[3];
void funcs2beZeroedDistNegNeighPlane5(int n, double x[], double fvec[], int i)
{
  int k1, k2; 
  double fx[3], rD[3];
  /* x = (r, alpha, t) */ 
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
	}
      rD[k1] = x[k1] + gradplane[k1]*x[4];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gradplane[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}

void fdjacDistNegNeighPlane5(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3], rD[3];
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	}
    }
  /* calc fx*/
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	}
      rD[k1] = x[k1] + gradplane[k1]*x[4];
    } 
  //printf("rC: %f %f %f rD: %f %f %f\n", x[0], x[1], x[2], rD[0], rD[1], rD[2]);
  //printf("fx: %f %f %f x[4]: %f\n", fx[0], fx[1], fx[2], x[4]);
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  df[3][3] = 0.0;
  df[3][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = gradplane[k1];
    } 
  df[4][3] = 0.0;
  df[4][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    df[4][4] += gradplane[k1]*gradplane[k1];
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = -2.0*x[3]*gradplane[k1];
      df[k1][4] = 0.0;
    } 

#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gradplane[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
#endif
}

void fdjacDistNegNeighPlane(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3];
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = 0;
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      //gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  //gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3]-rB[k2]);
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
      df[4][k1+3] = gradplane[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = -2.0*x[6]*gradplane[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1.0; //+ 2.0*x[7]*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 0.0;//2.0*x[7]*Xa[k1][k2];
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
    df[k1+5][7] = gradplane[k1];//fx[k1];
#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gradplane[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 //fvec[4] = 0.5*fvec[4]-1.0;
 /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
}
void fdjacDistNegNeigh(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3], gx[3];
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = -2.0*Sqr(x[6])*Xb[k1][k2];
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
      df[k1][6] = -2.0*x[6]*gx[k1];
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
#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gx[k1];
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
void funcs2beZeroedDistNegNeighPlane(int n, double x[], double fvec[], int i)
{
  int k1, k2; 
  double fx[3];
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
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3] - rB[k2]);
    }
#endif
   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gradplane[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  //fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//fx[k1]*x[7]; 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
void funcs2beZeroedDistNegNeigh(int n, double x[], double fvec[], int i)
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
      fvec[k1] = fx[k1] - Sqr(x[6])*gx[k1];
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
void calc_intersec_neigh_plane(double *rA, double *rB, double **Xa, double *grad, double* rC, double* rD)
{
  double A, B=0.0, C=0.0, D=0.0, tt=0.0;
  double rBA[3], rBAgrad;
  int k1, k2;
  /* la direzione è quella perpendicolare al piano */
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = grad[k1];
  MD_DEBUG(printf("rBA=(%f,%f,%f)\n", rBA[0], rBA[1], rBA[2])); 
  MD_DEBUG(printf("rB= (%f,%f,%f)\n", rB[0], rB[1], rB[2]));
  A = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      A += rBA[k1]*Xa[k1][k2]*rBA[k2];
 
  if (A <= 0)
    {
      printf("NNL [calc_intersec] Serious problem guessing distance, aborting...\n");
      printf("tt = %f D=%f A=%f B=%f C=%f\n", tt, D, A, B, C);
      printf("grad = (%.10f, %.10f, %.10f\n", grad[0], grad[1], grad[2]);
      printf("distance: %f\n", sqrt(Sqr(rBA[0])+Sqr(rBA[1])+Sqr(rBA[2])));
      print_matrix(Xa,3);
      print_matrix(Xb,3);
      exit(-1);
    }
  tt = sqrt(1 / A); 
  for (k1 = 0; k1 < 3; k1++)
    {
      rC[k1] = rA[k1] + tt*rBA[k1];  
    }

  /* ...e ora calcoliamo rD (guess sul piano) */
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = rB[k1] - rA[k1];
  rBAgrad = scalProd(rBA, grad);
  for (k1=0; k1 < 3; k1++)
    rD[k1] = rA[k1] + rBAgrad*grad[k1];
}
void guess_distNeigh_plane(int i, 
		double *rA, double *rB, double **Xa, double *grad, double *rC, double *rD,
		double **RA)
{
  double gradA[3], gradaxA[3], dA[3], dB[3];
  int k1, n;
  double saA[3], saB[3], sp;

  //printf("===============>SONO QUI\n");
  saA[0] = axa[i];
  saA[1] = axb[i];
  saA[2] = axc[i];
  saB[0] = nebrTab[i].axa;
  saB[1] = nebrTab[i].axb;
  saB[2] = nebrTab[i].axb;
  for (k1 = 0; k1 < 3; k1++)
    gradA[k1] =  grad[k1];
  for (n = 0; n < 3; n++)
    {
      gradaxA[n] = 0;
      for (k1 = 0; k1 < 3; k1++) 
	gradaxA[n] += gradA[k1]*RA[n][k1];
    }
  for (k1=0; k1 < 3; k1++)
    {
      dA[k1] = rA[k1];
      for (n=0; n < 3;n++)
	dA[k1] += gradaxA[n]*RA[n][k1]*saA[n]/2.0; 
    }
  calc_intersec_neigh(dA, rA, Xa, rC, 1);
  for (k1=0; k1 < 3; k1++)
    dB[k1] = rB[k1] - rC[k1];
  sp = scalProd(dB, grad);
  for (k1=0; k1 < 3; k1++)
    rD[k1] = rC[k1] + sp*grad[k1];
}
void calc_intersec_neigh(double *rB, double *rA, double **Xa, double* rI, double alpha)
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
      A += rBA[k1]*Xa[k1][k2]*rBA[k2];
 
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
      rI[k1] = rA[k1] + alpha*tt*rBA[k1];  
    }
}
void guess_distNeigh(int i, 
		double *rA, double *rB, double **Xa, double **Xb, double *rC, double *rD,
		double **RA, double **RB)
{
  double gradA[3], gradB[3], gradaxA[3], gradaxB[3], dA[3], dB[3];
  int k1, n;
  double saA[3], saB[3];

  //printf("===============>SONO QUI\n");
  saA[0] = axa[i];
  saA[1] = axb[i];
  saA[2] = axc[i];
  saB[0] = nebrTab[i].axa;
  saB[1] = nebrTab[i].axb;
  saB[2] = nebrTab[i].axb;
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
  calc_intersec_neigh(dA, rA, Xa, rC, -1);
  calc_intersec_neigh(dB, rB, Xb, rD, 1);
}
double calcDistNegNeighPlane(double t, double t1, int i, double *r1, double *r2, double *vecgsup, int calcguess, int calcgradandpoint, int *err, int nplane)

{
  /* NOTA: nplane = {0...7} e indica il piano rispetto al quale dobbiamo calcolare la distanza */
  double vecg[8], rC[3], rD[3], rDC[3], r12[3], invaSqN, invbSqN, invcSqN;
  double ti, segno;
  int retcheck;
  double nf, ng, gradf[3];
#ifdef MD_ASYM_ITENS
  double psi, phi;
#else
  double Omega[3][3];
#endif
  int k1;
  MD_DEBUG20(printf("t=%f tai=%f i=%d\n", t, t+t1-atomTime[i],i));
  MD_DEBUG20(printf("v = (%f,%f,%f)\n", vx[i], vy[i], vz[i]));
  *err = 0;
  ti = t + (t1 - atomTime[i]);
  //printf("t1-atomTime[%d]:%.15G\n", i, t1-atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG20(printf("AAAA ti= %.15G rA (%.15G,%.15G,%.15G)\n", ti, rA[0], rA[1], rA[2]));
  MD_DEBUG20(printf("AAAA t1=%.15G atomTime[%d]=%.15G\n",t1,i,atomTime[i]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, RtA, Omega);
#endif
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);

  tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, RtA);
  //printf("ti= %.15G rNebrShell: %f\n", ti, OprogStatus.rNebrShell);
  ti = 0.0;
  MD_DEBUG20(printf("BBBB ti= %.15G rB (%.15G,%.15G,%.15G)\n", ti, rB[0], rB[1], rB[2]));
  /* NOTA: dato l'ellissoide e la sua neighbour list a t=0 bisognerebbe stimare con esattezza 
   * la loro distanza e restituirla di seguito */
  /* NOTA2: per ora fa in ogni caso il guess ma si potrebbe facilmente fare in modo
   * di usare il vecchio vecg. */
  if (calcgradandpoint)
    calc_grad_and_point_plane(i, gradplane, rB, nplane);
  
  if (OprogStatus.guessDistOpt==1)
    {
      guess_distNeigh_plane(i, rA, rB, Xa, gradplane, rC, rD, RtA);
    }
  else
    {
      calc_intersec_neigh_plane(rA, rB, Xa, gradplane, rC, rD);
    }
  for(k1=0; k1 < 3; k1++)
    r12[k1] = rC[k1]-rD[k1]; 
  MD_DEBUG(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		  rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
  calc_grad(rC, rA, Xa, gradf);
  MD_DEBUG(printf("gradf=(%f,%f,%f) gradplane=(%f,%f,%f)\n",
		  gradf[0], gradf[1], gradf[2], gradplane[0], gradplane[1], gradplane[2]));
  nf = calc_norm(gradf);
  ng = calc_norm(gradplane);
  if (OprogStatus.dist5NL)
    vecg[3] = sqrt(nf/ng);
  else
    vecg[6] = sqrt(nf/ng);
  for (k1=0; k1 < 3; k1++)
    {
      vecg[k1] = rC[k1];
      if (!OprogStatus.dist5NL)
	vecg[k1+3] = rD[k1];
      rDC[k1] = rD[k1] - rC[k1];
    }

  if (OprogStatus.dist5NL)
    vecg[4] = 0.0;
  /*vecg[4] = calc_norm(rDC);
   * QUESTO GUESS POTREBBE ESSERE MIGLIORE ANCHE SE IN PRATICA SEMBRA
   * LO STESSO. */
  else
    vecg[7] = 0.0;
  MD_DEBUG(printf("alpha: %f beta: %f\n", vecg[6], vecg[7]));
  if (OprogStatus.dist5NL)
    newtDistNegNeighPlane(vecg, 5, &retcheck, funcs2beZeroedDistNegNeighPlane5, i); 
  else
    newtDistNegNeighPlane(vecg, 8, &retcheck, funcs2beZeroedDistNegNeighPlane, i); 
  if (retcheck != 0)
    {
      if (OprogStatus.targetPhi>0)
	{
	  calcdist_retcheck=1;
	  return 0.0;
	}
      printf("[NNL] I couldn't calculate distance between %d and its NL, calcguess=%d, exiting....\n", i, calcguess);
      exit(-1);
    }
  if (OprogStatus.dist5NL)
    {
      for (k1 = 0; k1 < 5; k1++)
	{
	  vecgsup[k1] = vecg[k1]; 
	}
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	{
	  vecgsup[k1] = vecg[k1]; 
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      if (!OprogStatus.dist5NL)
	r2[k1] = vecg[k1+3];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      if (OprogStatus.dist5NL)
	{
	  r2[k1] = r1[k1] + gradplane[k1]*vecg[4];
	  vecgsup[k1+5] = r2[k1];
	}
      r12[k1] = r1[k1] - r2[k1];
    } 
  if (OprogStatus.dist5NL)
    segno = vecg[4];
  else
    segno = vecg[7];
  if (segno > 0)
    {
      //printf("t=%.15G distanza: %.15G\n", t, calc_norm(r12));
      return calc_norm(r12);
    }
  else
    {
      //printf("t=%.15G distanza: %.15G\n", t, -calc_norm(r12));
      return -calc_norm(r12);
    }
}
double calcDistNegNeigh(double t, double t1, int i, double *r1, double *r2, double *vecgsup, int calcguess, int ignorefail, int *err)

{
  double vecg[8], rC[3], rD[3], rDC[3], r12[3], vecgcg[6], invaSqN, invbSqN, invcSqN;
  double shift[3] = {0.0, 0.0, 0.0};
  double ti, segno;
  int retcheck;
  double nf, ng, gradf[3], gradg[3];
#ifdef MD_ASYM_ITENS
  double psi, phi;
#else
  double Omega[3][3]; 
#endif
  int k1;
  MD_DEBUG20(printf("t=%f tai=%f i=%d\n", t, t+t1-atomTime[i],i));
  MD_DEBUG20(printf("v = (%f,%f,%f)\n", vx[i], vy[i], vz[i]));
  *err = 0;
  ti = t + (t1 - atomTime[i]);
  //printf("t1-atomTime[%d]:%.15G\n", i, t1-atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG20(printf("AAAA ti= %.15G rA (%.15G,%.15G,%.15G)\n", ti, rA[0], rA[1], rA[2]));
  MD_DEBUG20(printf("AAAA t1=%.15G atomTime[%d]=%.15G\n",t1,i,atomTime[i]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, RtA, Omega);
#endif
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);

  tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, RtA);
  //printf("ti= %.15G rNebrShell: %f\n", ti, OprogStatus.rNebrShell);
  ti = 0.0;
  rB[0] = nebrTab[i].r[0];
  rB[1] = nebrTab[i].r[1];
  rB[2] = nebrTab[i].r[2];
  MD_DEBUG20(printf("BBBB ti= %.15G rB (%.15G,%.15G,%.15G)\n", ti, rB[0], rB[1], rB[2]));
  /* NOTA: dato l'ellissoide e la sua neighbour list a t=0 bisognerebbe stimare con esattezza 
   * la loro distanza e restituirla di seguito */

  if ((Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]))==0.0 || (rA[0]==rB[0] && rA[1]==rB[1] && rA[2]==rB[2]))
    {
      //firstDist = 1;
      return OprogStatus.rNebrShell*(min3(axa[i],axb[i],axc[i])/max3(axa[i],axb[i],axc[i]));
    }
  invaSqN = 1.0/Sqr(nebrTab[i].axa);
  invbSqN = 1.0/Sqr(nebrTab[i].axb);
  invcSqN = 1.0/Sqr(nebrTab[i].axc);
  tRDiagR(i, Xb, invaSqN, invbSqN, invcSqN, nebrTab[i].R);
retryneigh:
  //printf("<<<<<<<CALCGUESS=%d forceguess=%d>>>>>>>>\n", calcguess, OprogStatus.forceguess);
  if (OprogStatus.forceguess)
    calcguess = 1;
  if (calcguess)
    {
#if 0
      if (firstDist)
	{
	  double rBtmp[3]={1.0,0.0,0.0}, rAtmp[3]={-1.0,0.0,0.0};
	  calc_intersec_neigh(rBtmp, rAtmp, Xa, rC, -1);
	  calc_intersec_neigh(rAtmp, rBtmp, Xb, rD, 1);
	}
      else
#endif
      if (OprogStatus.guessDistOpt==1)
	guess_distNeigh(i, rA, rB, Xa, Xb, rC, rD, RtA, nebrTab[i].R);
      else
	{
	  calc_intersec_neigh(rB, rA, Xa, rC, -1);
	  calc_intersec_neigh(rA, rB, Xb, rD, 1);
	}
     for(k1=0; k1 < 3; k1++)
	r12[k1] = rC[k1]-rD[k1]; 
      if ((OprogStatus.SDmethod==1 || OprogStatus.SDmethod==4) && OprogStatus.stepSDA>0 && OprogStatus.stepSDB>0)
	{
	  for (k1=0; k1 < 3; k1++)
	    {
	      vecgcg[k1] = rC[k1];
	      vecgcg[k1+3] = rD[k1];
	    }
	  distSD(i, i, shift, vecgcg, OprogStatus.springkSD, 1);
	  for (k1=0; k1 < 3; k1++)
	    {
	      rC[k1] = vecgcg[k1];
	      rD[k1] = vecgcg[k1+3];
	    }	 
	}
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
#if 0
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
	    g2 = OprogStatus.epsdGDO*min3(axa[j],axb[j],axc[j])/calc_norm(vecnf); 
	  else 
	    g2 = g1;
	}	  
#endif
      if (OprogStatus.SDmethod==1 || OprogStatus.SDmethod==4)
	{
	  if (scalProd(gradf, rDC) < 0.0)
	    vecg[7] = 0.0;
	  else
	    vecg[7] = calc_norm(rDC)/nf;  
	}
      else
	{
#if 0	
	  if (OprogStatus.epsdGDO > 0.0)
	    {
	      if (scalProd(gradf, rDC) < 0.0)
		vecg[7] = 0.0;
	      else
		vecg[7] = min(g1,g2);
	    }
#endif
	  vecg[7] = 0.0;
	}
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  MD_DEBUG(printf("alpha: %f beta: %f\n", vecg[6], vecg[7]));
  newtDistNegNeigh(vecg, 8, &retcheck, funcs2beZeroedDistNegNeigh, i); 
  if (retcheck != 0)
    {
      if (OprogStatus.targetPhi>0)
	{
	  calcdist_retcheck=1;
	  return 0.0;
	}
      printf("[NNL] I couldn't calculate distance between %d and its NL, calcguess=%d, exiting....\n", i, calcguess);
      if (calcguess==0)
	{
	  calcguess=2;
	  goto retryneigh;
	} 
      Oparams.time = t + t1;
      //store_bump(i, j);
      if (ignorefail)
	{
	  *err = 1;
	  return 0.0;
	}
      else
	exit(-1);
    }
  for (k1 = 0; k1 < 8; k1++)
    {
      vecgsup[k1] = vecg[k1]; 
    }  
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      r2[k1] = vecg[k1+3];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      r12[k1] = r1[k1] - r2[k1];
    } 
  segno = vecg[7];
#if 1
    {
      int kk;
      double rCA[3], rDB[3], grA[3], grB[3];
      for (kk=0; kk < 3; kk++)
	{
	  rCA[kk] = rC[kk] - rA[kk];
	  rDB[kk] = rD[kk] - rB[kk];
	}
      calc_grad(r1, rA, Xa, grA);
      calc_grad(r2, rB, Xb, grB);
      printf("rC = (%.12f,%.12f,%.12f) norm(rCA)=%f rD = (%.12f, %.12f, %.12f) norm(rDB)=%.14f scalprod: %f\n"
	     , rC[0]-rA[0], rC[1]-rA[1], rC[2]-rA[2], calc_norm(rCA),  
	     rD[0]-rB[0], rD[1]-rB[1], rD[2]-rB[2], calc_norm(rDB), scalProd(grA, grB));
    }
#endif
 
#if 0
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (r2[k1]-rA[k1])*Xa[k1][k2]*(r2[k2]-rA[k2]);
#endif
#if 0
  printf("=====> r1=%f %f %f r2=%f %f %f segno=%.15G\n", r1[0], r1[1], r1[2], r2[0],
	 r2[1], r2[2], segno);
  if (t+t1 > 30.0)//(fabs(segno) > 0.5 && segno < 0.0)
    {
      double check;
      Oparams.time = t+t1;
      store_bump_neigh(i, r1, r2);
      check = -1;
      /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++) 
	  check += (r2[k1]-rB[k1])*Xb[k1][k2]*(r2[k2]-rB[k2]);
      printf("===>check=%.15G\n", check);
      //exit(-1);	
    }
#endif
#if 0
  if (segno*vecg[7]<0 && fabs(segno*vecg[7])>3E-8)
    {
#if 0
      if (OprogStatus.targetPhi>0)
	{
	  calcdist_retcheck = 1;
	  return 0.0;
	}
#endif
      printf("segno: %.8G vecg[7]: %.8G dist=%.15G\n", segno, vecg[7], calc_norm(r12));
      exit(-1);
      //return calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
    }
#endif
  if (segno > 0)
    {
      printf("t=%.15G distanza: %.15G\n", t, calc_norm(r12));
      return calc_norm(r12);
    }
  else
    {
      printf("t=%.15G distanza: %.15G\n", t, -calc_norm(r12));
      return -calc_norm(r12);
    }
}
double calcDistNegNNLoverlapPlane(double t, double t1, int i, int j, double shift[3])
{
  double RR, R0, R1, cij[3][3], fabscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k1, k2, existsParallelPair = 0;
#if 0
  double distBuf;
  distBuf = OprogStatus.epsd;
#endif
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  rA[0] = nebrTab[i].r[0];
  rA[1] = nebrTab[i].r[1];
  rA[2] = nebrTab[i].r[2];
  rB[0] = nebrTab[j].r[0] + shift[0];
  rB[1] = nebrTab[j].r[1] + shift[1];
  rB[2] = nebrTab[j].r[2] + shift[2];
 
  EA[0] = nebrTab[i].axa;
  EA[1] = nebrTab[i].axb;
  EA[2] = nebrTab[i].axc;
  EB[0] = nebrTab[j].axa;
  EB[1] = nebrTab[j].axb;
  EB[2] = nebrTab[j].axc;
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      EA[k1] += distBuf;
      EB[k1] += distBuf;
    }
#endif
  /* verificare che AA e BB sono effettivamente gli assi principali degli ellissoidi */
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  AA[k1][k2] = nebrTab[i].R[k1][k2];
	  BB[k1][k2] = nebrTab[j].R[k1][k2];
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
double calcDistNegNNLoverlap(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha,

     		double *vecgsup, int calcguess)
{
  double vecg[8], rC[3], rD[3], rDC[3], r12[3], vecgcg[6], fx[3];
  double ti, segno;
  double g1=0.0, g2=0.0, SP, nrDC, vecnf[3], nvecnf;
  int retcheck;
  double nf, ng, gradf[3], gradg[3], invaSq, invbSq, invcSq;
  int k1, k2;
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  ti = nebrTab[i].time - Oparams.time;
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
  invaSq = 1/Sqr(nebrTab[i].axa);
  invbSq = 1/Sqr(nebrTab[i].axb);
  invcSq = 1/Sqr(nebrTab[i].axc);
  tRDiagR(i, Xa, invaSq, invbSq, invcSq, nebrTab[i].R);

  ti = nebrTab[j].time - Oparams.time;
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  invaSq = 1/Sqr(nebrTab[j].axa);
  invbSq = 1/Sqr(nebrTab[j].axb);
  invcSq = 1/Sqr(nebrTab[j].axc);
    
  tRDiagR(j, Xb, invaSq, invbSq, invcSq, nebrTab[j].R);
retryoverlap:
  if (OprogStatus.forceguess)
    calcguess = 1;
  if (calcguess)
    {
      if (OprogStatus.guessDistOpt==1)
	guess_dist(i, j, rA, rB, Xa, Xb, rC, rD, nebrTab[i].R, nebrTab[j].R);
      else
	{
	  calc_intersec(rB, rA, Xa, rC);
	  calc_intersec(rA, rB, Xb, rD);
	}
#if 1

      for(k1=0; k1 < 3; k1++)
	r12[k1] = rC[k1]-rD[k1]; 

      if ((OprogStatus.SDmethod==1||OprogStatus.SDmethod==4) && OprogStatus.stepSDA>0 && OprogStatus.stepSDB>0)
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
	  //distSD(i, j, shift, vecgcg, 100, 1); 
	  //distSD(i, j, shift, vecgcg, 1/OprogStatus.tolSD, 1);
	  //guessdistByMesh(i, j, shift, vecgcg);
	  distSD(i, j, shift, vecgcg, OprogStatus.springkSD, 1);
#if 0
	  if (maxitsRyck)
	    printf("distVera=%.15G\n", calcDist(t, i, j, shift, r1, r2, alpha, vecgsup, 1));
#endif
	  for (k1=0; k1 < 3; k1++)
	    {
	      rC[k1] = vecgcg[k1];
	      rD[k1] = vecgcg[k1+3];
	    }	 
#endif
	}
#if 0
	{
	  double dist, distVera;
	  for(k1=0; k1 < 3; k1++)
	    r12[k1] = rC[k1]-rD[k1]; 
	  dist = calc_norm(r12);
	  printf("dist: %.15G\n", dist);
	      //printf("dist=%.15f\n",calc_norm(r12));
	}
#endif
      //exit(-1);
      MD_DEBUG(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		      rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
      calc_grad(rC, rA, Xa, gradf);
      calc_grad(rD, rB, Xb, gradg);
      MD_DEBUG(printf("gradf=(%f,%f,%f) gradg=(%f,%f,%f)\n",
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
	    g2 = OprogStatus.epsdGDO*min3(axa[j],axb[j],axc[j])/calc_norm(vecnf); 
	  else 
	    g2 = g1;
	}	  
      if (OprogStatus.dist5)
	{
	  if (OprogStatus.SDmethod==1||OprogStatus.SDmethod==4)
	    if (scalProd(gradf, rDC) < 0.0)
	      vecg[4] = 0.0;
	    else
	      vecg[4] = calc_norm(rDC)/nf;  
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
	}
      else
	{
	  if (OprogStatus.SDmethod==1||OprogStatus.SDmethod==4)
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
	}
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  MD_DEBUG(printf(">>>>>>> alpha: %f beta: %f\n", vecg[6], vecg[7]));
  //printf(">>>>>>> MAHHH alpha: %f beta: %f\n", vecg[6], vecg[7]);
  if (OprogStatus.dist5)
    newtDistNeg(vecg, 5, &retcheck, funcs2beZeroedDistNeg5, i, j, shift); 
  else
    newtDistNeg(vecg, 8, &retcheck, funcs2beZeroedDistNeg, i, j, shift); 
  if (retcheck != 0)
    {
      if (OprogStatus.targetPhi>0)
	{
	  calcdist_retcheck=1;
	  return 0.0;
	}
      printf("I couldn't calculate distance between %d and %d\n, exiting....\n", i, j);
      if (calcguess==0)
	{
	  calcguess=2;
	  goto retryoverlap;
	} 
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
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      if (OprogStatus.dist5)
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
#if 1
  if (OprogStatus.dist5)
    {
      if (segno*vecg[4]<0 && fabs(segno*vecg[4])>3E-8)
	{
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
	  if (OprogStatus.targetPhi>0)
	    {
	      calcdist_retcheck = 1;
	      return 0.0;
	    }
	  printf("segno: %.8G vecg[7]: %.8G dist=%.15G\n", segno, vecg[7], calc_norm(r12));
	  return calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
	  //exit(-1);
	}
    }
#endif
  if (segno > 0)
    return calc_norm(r12);
  else
    return -calc_norm(r12);
}

double gradplane_all[6][3], rBall[6][3];
void assign_plane(int nn)
{
  int kk;
  for (kk = 0; kk < 3; kk++)
    {
      gradplane[kk] = gradplane_all[nn][kk];
      rB[kk] = rBall[nn][kk];
    }
}
extern double sogliaErr_zbrent;
int interpolNeighPlane(int i, double tref, double t, double delt, double d1, double d2, double *troot, double* vecg, int bracketing, int nplane)
{
  int nb, distfail;
  double d3, t1, t2;
  double r1[3], r2[3], xb1[2], xb2[2];
  if (OprogStatus.paralNNL)
    assign_plane(nplane);
  sogliaErr_zbrent = OprogStatus.epsdNL;
  d3 = calcDistNegNeighPlane(t+delt*0.5, tref, i, r1, r2, vecg, 0, 0, &distfail, nplane);
  xa[0] = t;
  ya[0] = d1;
  xa[1] = t+delt*0.5;
  ya[1] = d3;
  xa[2] = t+delt;
  ya[2] = d2;
  polinterr = 0;
  if (!bracketing)
    {
      t1 = t;
      t2 = t+delt;
    }
  else
    {
      t1 = t;
      t2 = t+delt;
      nb = 1;
      zbrak(distfunc, t1, t2, OprogStatus.zbrakn, xb1, xb2, &nb);
      if (nb==0 || polinterr==1)
	{
	  return 1;
	}
      t1 = xb1[0];
      t2 = xb2[0];
    }
  if (polinterr)
    return 1;
  *troot=zbrent(distfunc, t1, t2, OprogStatus.zbrentTol);
  if (polinterr)
    {
      printf("[interpol_neigh] bracketing: %d polinterr=%d t1=%.15G t2=%.15G\n", bracketing,polinterr, t1, t2);
      printf("d: %.15G,%.15G,%.15G\n", d1, d3, d2);
      printf("tref=%.15G t: %.15G,%.15G,%.15G\n", tref, t, t+delt*0.5, t+delt);
      printf("distfunc(t1)=%.10G distfunc(t2)=%.10G\n", distfunc(t), distfunc(t+delt));
      return 1;
    }
  if ((*troot < t && fabs(*troot-t)>3E-8) || (*troot > t+delt && fabs(*troot - (t+delt))>3E-8))
    {
      printf("[interpol neigh] brack: %d xb1: %.10G xb2: %.10G\n", bracketing, xb1[0], xb2[0]);
      printf("*troot: %.15G t=%.15G t+delt:%.15G\n", *troot, t, t+delt);
      printf("d1=%.10G d2=%.10G d3:%.10G\n", d1, d2, d3);
      printf("distfunc(t1)=%.10G distfunc(t2)=%.10G\n", distfunc(t), distfunc(t+delt));
      printf("distfunc(t+delt*0.5)=%.10G\n", distfunc(t+delt*0.5));
      return 1;
    }
  calcDistNegNeighPlane(*troot, tref, i, r1, r2, vecg, 0, 0, &distfail, nplane);
  *troot += tref;
  return 0;
}
int interpolNeigh(int i, double tref, double t, double delt, double d1, double d2, double *troot, double* vecg, int bracketing)
{
  int nb, distfail;
  double d3, t1, t2;
  double r1[3], r2[3], xb1[2], xb2[2];
  d3 = calcDistNegNeigh(t+delt*0.5, tref, i, r1, r2, vecg, 0, 0, &distfail);
  sogliaErr_zbrent = OprogStatus.epsdNL;
  xa[0] = t;
  ya[0] = d1;
  xa[1] = t+delt*0.5;
  ya[1] = d3;
  xa[2] = t+delt;
  ya[2] = d2;
  polinterr = 0;
  if (!bracketing)
    {
      t1 = t;
      t2 = t+delt;
    }
  else
    {
      t1 = t;
      t2 = t+delt;
      nb = 1;
      zbrak(distfunc, t1, t2, OprogStatus.zbrakn, xb1, xb2, &nb);
      if (nb==0 || polinterr==1)
	{
	  return 1;
	}
      t1 = xb1[0];
      t2 = xb2[0];
    }
  if (polinterr)
    return 1;
  *troot=zbrent(distfunc, t1, t2, OprogStatus.zbrentTol);
  if (polinterr)
    {
      printf("bracketing: %d polinterr=%d t1=%.15G t2=%.15G\n", bracketing,polinterr, t1, t2);
      printf("d: %.15G,%.15G,%.15G\n", d1, d3, d2);
      printf("t: %.15G,%.15G,%.15G\n", t, t+delt*0.5, t+delt);
      printf("distfunc(t1)=%.10G distfunc(t2)=%.10G\n", distfunc(t), distfunc(t+delt));
      return 1;
    }
  if ((*troot < t && fabs(*troot-t)>3E-8) || (*troot > t+delt && fabs(*troot - (t+delt))>3E-8))
    {
      printf("brack: %d xb1: %.10G xb2: %.10G\n", bracketing, xb1[0], xb2[0]);
      printf("*troot: %.15G t=%.15G t+delt:%.15G\n", *troot, t, t+delt);
      printf("d1=%.10G d2=%.10G d3:%.10G\n", d1, d2, d3);
      printf("distfunc(t1)=%.10G distfunc(t2)=%.10G\n", distfunc(t), distfunc(t+delt));
      printf("distfunc(t+delt*0.5)=%.10G\n", distfunc(t+delt*0.5));
      return 1;
    }
  calcDistNegNeigh(*troot, tref, i, r1, r2, vecg, 0, 0, &distfail);
  *troot += tref;
  return 0;
}
double calcvecFNeigh(int i, double t, double t1, double* ddot, double *r1)
{
  int kk;
  double rcat[3], wra[3], dti;
  dti = t1 - atomTime[i];
  rcat[0] = r1[0] - (rx[i] + vx[i]*(t+dti)); 
  rcat[1] = r1[1] - (ry[i] + vy[i]*(t+dti));
  rcat[2] = r1[2] - (rz[i] + vz[i]*(t+dti));
  ddot[0] = vx[i];
  ddot[1] = vy[i];
  ddot[2] = vz[i];
  vectProd(wx[i], wy[i], wz[i], rcat[0], rcat[1], rcat[2], &wra[0], &wra[1], &wra[2]);
  for (kk=0; kk < 3; kk++)
    ddot[kk] += wra[kk];
  //printf("^^^ sp = %.15G norm= %.15G\n", fabs(scalProd(ddot, gradplane)), calc_norm(ddot));
#if 0
  return calc_norm(ddot);
#else
  /* WARNING:  quella che conta è la velocità rispetto alla perpendicolare al piano, 
   * poiché la distanza tra piano ed ellisse è sempre perpendicolare al piano...tuttavia
   * verificare e testare tale maggiorazione !!! */
  return fabs(scalProd(ddot, gradplane));
#endif
}

int search_contact_faster_neigh_plane(int i, double *t, double t1, double t2, 
				double *vecgd, double epsd, double *d1, double epsdFast,
				double *r1, double *r2, int nplane)
{
  /* NOTA: 
   * MAXOPTITS è il numero massimo di iterazioni al di sopra del quale esce */
  double maxddot, told, delt, normddot, ddot[3];
  const int MAXOPTITS = 500;
  double factori;
  int its=0, distfailed, itsf=0; 
  const double GOLD= 1.618034;  
  factori = 0.5*maxax[i]+OprogStatus.epsdNL;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));

  /* estimate of maximum rate of change for d */
#if 0
  maxddot = sqrt(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
#else
  /* WARNING: questa maggiorazione si applica al caso specifico dell'urto di un ellissoide con un piano,
   * è molto migliore della precedente ma va testata accuratamente! */
#ifdef MD_ASYM_ITENS
  maxddot = calc_maxddot_nnl(i, gradplane);
#else
  maxddot = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
#endif
  //printf("SCALPROD = %.15G\n", vx[0]*gradplane[0]+vy[1]*gradplane[1]+vz[2]*gradplane[2]);
  *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed, nplane);
  timesFNL++;
  MD_DEBUG20(printf("Pri distances between %d d1=%.12G epsd*epsdTimes:%f\n", i, *d1, epsdFast));
  //printf("[SEARCH_CONTACT_FASTER] t=%.15G ellips N. %d d=%.15G\n", *t, i, *d1); 
  while (*d1 > epsdFast && its < MAXOPTITS)
    {
      told = *t;
      delt = *d1 / maxddot;
      normddot = calcvecFNeigh(i, *t, t1, ddot, r1);
      //printf("normddot: %.15G\n", epsd/normddot);
      /* check for convergence */
     
      if (normddot!=0 && delt < (epsd / normddot))
	{
	  MD_DEBUG20(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
      *t += delt;
#if 1
      if (*t + t1 > t2)
	{
	  *t = told;
	  MD_DEBUG20(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG20(printf("convergence t>t2\n"));
	  *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed, nplane);
	  return 1;
	}
#endif
      *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 1, &distfailed, nplane);
      //printf("NNL LOOP SEARCH CONTACT FASTER i=%d *d1=%.15G *t=%.15G\n", i, *d1, *t+t1);
      /* NOTA: nel caso di urto di un ellissoide con la sua neighbour list 
       * a t=0 i due ellissoidi hanno stesso centro e assi principali paralleli
       * e questo implica che ci possono essere più soluzioni possibili per le eq.
       * che definiscono la distanza.
       * Questo implica che la distanza che ottengo è una sovrastima di quella vera e quindi 
       * con il loop che segue compenso questo problema riducendo il passo fino a che 
       * non arrivo ad una distanza positiva con il passo veloce */
#if 1
      itsf = 0;
      while (*d1 < 0 || distfailed)
	{
	  /* reduce step size */
	  if (itsf == 0 && delt - OprogStatus.h > 0)
	    delt -= max(OprogStatus.h,delt*OprogStatus.h);
	  else
	    delt /= GOLD;
	  
	  *t = told + delt;
	  *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 1, &distfailed, nplane);
	  //printf("itsf=%d BUBU SEARCH_CONTACT_FASTER_NEIGH *d1=%.15G\n",itsf,*d1);
	  itsf++;	
	  if (itsf > 100)
	    {
	      printf("*d1=%.15G too many times calculation of distance failed!\n", *d1);
	      printf("aborting...\n");
	      exit(-1);
	    }
	}
#else
      if (*d1 < 0)
	{
	  /* go back! */
	  MD_DEBUG31(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG31(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed, nplane);
	  return 0;
	}
#endif
      told = *t;
      its++;
      itsFNL++;
    }

  MD_DEBUG20(printf("max iterations %d iterations reached t=%f t2=%f\n", its, *t, t2));
  return 0;

}
int refine_contact_neigh_plane(int i, double t1, double t, double vecgd[8], double  vecg[5],
			       int nplane)
{
  int kk, retcheck;

  for (kk = 0; kk < 3; kk++)
    {
      if (OprogStatus.dist5NL)
	vecg[kk] = (vecgd[kk]+vecgd[kk+5])*0.5; 
      else
	vecg[kk] = (vecgd[kk]+vecgd[kk+3])*0.5; 
    }
  if (OprogStatus.dist5NL)
    vecg[3] = vecgd[3];
  else
    vecg[3] = vecgd[6];
  vecg[4] = t-t1;
  trefG = t1;
  newtNeigh(vecg, 5, &retcheck, funcs2beZeroedNeighPlane, i); 
  vecg[4] += t1;
  if (retcheck==2)
    {
      MD_DEBUG10(printf("newtNeigh did not find any contact point!\n"));
      return 0;
    }
  else
    {
      return 1; 
    }
}
int bracket_neigh(int i, double *t1, double *t2, int nplane)
{
  double dd, veln, vel[3], r1[3], dr[3], fabsveln;
  int kk;
  vel[0] = vx[i];
  vel[1] = vy[i];
  vel[2] = vz[i];
  veln = scalProd(vel, gradplane);
  
  r1[0] = rx[i];
  r1[1] = ry[i];
  r1[2] = rz[i];
  for (kk=0; kk < 3; kk++)
    dr[kk] = r1[kk] - rB[kk]; 
  dd = fabs(scalProd(dr, gradplane));
  if (veln == 0.0)
    return 0;
  if (veln < 0 && dd > maxax[i])
    return 0;
  fabsveln = fabs(veln);
  if (dd > maxax[i])
    *t1 = (dd - maxax[i]) / fabsveln;
  else
    *t1 = 0.0;
  *t2 = (dd + maxax[i]) / fabsveln;
  if (*t2 < *t1)
    {
      printf("problema nel bracketing per i=%d, t1=%.15G t2=%.15G\n", i, *t1, *t2);
      exit(-1);
    }
  return 1;
}
int check_cross_scf(double distsOld[6], double dists[6], int crossed[6])
{
  int nn;
  int retcross = 0;
  for (nn = 0; nn < 6; nn++)
    {
      crossed[nn] = 0;
      //printf("dists[%d]=%.15G distsOld[%d]:%.15G\n", nn, dists[nn], nn, distsOld[nn]);
      if (fabs(dists[nn]) < 1E-14 && distsOld[nn] > 0.0)
	{
	  crossed[nn] = 1;
	  retcross = 1;
	}
    }
  return retcross;
}

int check_cross(double distsOld[6], double dists[6], int crossed[6])
{
  int nn;
  int retcross = 0;
  for (nn = 0; nn < 6; nn++)
    {
      crossed[nn] = 0;
      //printf("dists[%d]=%.15G distsOld[%d]:%.15G\n", nn, dists[nn], nn, distsOld[nn]);
      if (dists[nn] < 0.0  && distsOld[nn] > 0.0)
	{
	  crossed[nn] = 1;
	  //printf("CROSSED[%d]:%d\n", nn, crossed[nn]);
	  retcross = 1;
	}
    }
  return retcross;
}
double get_max_deldist(double distsOld[6], double dists[6])
{
  int nn, first = 1;
  double maxdd=0.0, dd;
  for (nn = 0; nn < 6; nn++)
    {
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
  memcpy(b, a, 6*sizeof(double));
}
int get_dists_tocheck(double distsOld[6], double dists[6], int tocheck[6], int dorefine[6])
{
  int nn;
  int rettochk = 0;
  for (nn = 0; nn < 6; nn++)
    {
      tocheck[nn] = 0;
      if (dists[nn] < OprogStatus.epsdNL && distsOld[nn] < OprogStatus.epsdNL &&
	  dorefine[nn] == 0)
	{
	  tocheck[nn] = 1; 
	  rettochk++;
	}
    }
  return rettochk;
}
double calcDistNegNeighPlaneAll(double t, double tref, int i, double dists[6], double vecgd[6][8], int calcguess)
{
  int nn, err, kk;
  double r1[3], r2[3], dmin=0.0;
  for (nn = 0; nn < 6; nn++)
    {
      for (kk = 0; kk < 3; kk++)
	{
	  gradplane[kk] = gradplane_all[nn][kk];
	  rB[kk] = rBall[nn][kk];
	}
      dists[nn] = calcDistNegNeighPlane(t, tref, i, r1, r2, vecgd[nn], calcguess, 0, &err, nn);
      //printf("dist[%d]:%.15G\n", nn, dists[nn]);
      if (nn==0 || dists[nn] < dmin)
	dmin = dists[nn];
    }
  return dmin; 
}
int check_distance(double maxddoti[6], double dists[6], double t1, double t2, double t)
{
  int nn, cc;
  cc = 0;
  for (nn = 0; nn < 6; nn++)
    {
      if ((t2 - (t1 + t))*maxddoti[nn] <  dists[nn] - OprogStatus.epsd)
	cc++;
    }
  if (cc == 6)
    return 1;
  else
    return 0;
  //printf("I chose dt=%.15G\n", *delt);
}
void calc_delt(double maxddoti[6], double *delt, double dists[6])
{
  int nn;
  double dt;
  for (nn = 0; nn < 6; nn++)
    {
      dt = fabs(dists[nn]) / maxddoti[nn];
      //printf("nn=%d dt=%.15G delt=%.15G dists=%.15G maxddoti=%15G\n", nn, dt, *delt, dists[nn], maxddoti[nn]);
      if (nn==0 || dt < (*delt))
	*delt = dt;
    }
  //printf("I chose dt=%.15G\n", *delt);
}
int search_contact_faster_neigh_plane_all(int i, double *t, double t1, double t2, 
					  double vecgd[6][8], double epsd, 
					  double *d1, double epsdFast, 
					  double dists[6], double maxddoti[6], double maxddot)
{
  double told, delt=1E-15, distsOld[6];
  const double GOLD= 1.618034;
  const int MAXOPTITS = 500;
  int its=0, crossed[6], itsf; 
  *d1 = calcDistNegNeighPlaneAll(*t, t1, i, distsOld, vecgd, 1);
#if 0
  if ((t2-t1)*maxddot < *d1 - OprogStatus.epsd)
    return 1;
#endif
  timesFNL++;
  told = *t;
  if (fabs(*d1) < epsdFast)
    {
      assign_dists(distsOld, dists);
      return 0;
    }
  while (fabs(*d1) > epsdFast && its < MAXOPTITS)
    {
#if 1
      calc_delt(maxddoti, &delt, distsOld);
      if (check_distance(maxddoti, dists, t1, t2, *t))
	return 1;
#else
      delt = fabs(*d1) / maxddot;
      if (*t + t1 < t2 && (t2 - (*t + t1))*maxddot < fabs(*d1) - OprogStatus.epsd)
	return 1;
#endif
      *t += delt;
      *d1 = calcDistNegNeighPlaneAll(*t, t1, i, dists, vecgd, 1);
#if 0
      if (its > 100 && its%10 == 0)
	{
	  printf("NNL SEARCH CONTACT FASTER t=%.15G its=%d\n", *t+t1, its);
	}
#endif
      //printf("d=%.15G t=%.15G\n", *d1, *t+t1);
#if 1
      itsf = 0;
      while (check_cross_scf(distsOld, dists, crossed))
	{
	  /* reduce step size */
	  if (itsf == 0 && delt - OprogStatus.h > 0)
	    delt -= max(OprogStatus.h, OprogStatus.h*delt);
	  else
	    delt /= GOLD;
	  *t = told + delt;
	  *d1 = calcDistNegNeighPlaneAll(*t, t1, i, dists, vecgd, 1);
	  itsf++;	
	  if (itsf > 100)
	    {
	      printf("*d1=%.15G too many times calculation of distance failed!\n", *d1);
	      printf("aborting...\n");
	      exit(-1);
	    }
	}
#else
     if (check_cross(distsOld, dists, crossed))
       {
	 /* go back! */
	 MD_DEBUG30(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	 MD_DEBUG30(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	 *t = told;	  
	 *d1 = calcDistNegNeighPlaneAll(*t, t1, i, dists, vecgd, 1);
	 return 0;
       }
#endif
      if (*t+t1 > t2)
	{
	  *t = told;
	  MD_DEBUG30(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG30(printf("convergence t>t2\n"));
	  *d1 = calcDistNegNeighPlaneAll(*t, t1, i, dists, vecgd, 1);
	  return 1;
	}
      told = *t;
      assign_dists(dists, distsOld);
      its++;
      itsFNL++;
    }
  return 0;
}
void assign_vec(double vecsource[6][8], double vecdest[6][8])
{
  int nn, kk;
  for (nn = 0; nn < 6; nn++)
    for (kk = 0; kk < 8; kk++)
      {
	vecdest[nn][kk] = vecsource[nn][kk];
      }
}
   	   
void calc_grad_and_point_plane_all(int i, double gradplaneALL[6][3], double rBALL[6][3])
{
  int nn;
  for (nn = 0; nn < 6; nn++)
    {
      calc_grad_and_point_plane(i, gradplaneALL[nn], rBALL[nn], nn);
    }
}
int locate_contact_neigh_plane_parall(int i, double *evtime, double t2)
{
  /* const double minh = 1E-14;*/
  double h, d, dold, t2arr[6], t, dists[6], distsOld[6], 
	 vecg[5], vecgroot[6][8], vecgd[6][8], vecgdold[6][8], factori; 
#ifndef MD_BASIC_DT
  double normddot, distsOld2[6], vecgdold2[6][8], dold2, deldist;
#endif
  double maxddot, delt, troot, tini, maxddoti[6];
  int firstev;
  /*
  const int MAXITS = 100;
  const double EPS=3E-8;*/ 
  /* per calcolare derivate numeriche questo è il magic number in doppia precisione (vedi Num. Rec.)*/
  int its, foundrc;
  double t1, epsd, epsdFast, epsdFastR, epsdMax; 
  int kk,tocheck[6], dorefine[6], ntc, ncr, nn, gotcoll, crossed[6], firstaftsf;
  epsd = OprogStatus.epsdNL;
  epsdFast = OprogStatus.epsdFastNL;
  epsdFastR= OprogStatus.epsdFastRNL;
  epsdMax = OprogStatus.epsdMaxNL;
  t = 0;//t1;
  t1 = Oparams.time;
  //t2 = timbig;
  calc_grad_and_point_plane_all(i, gradplane_all, rBall);
  factori = 0.5*maxax[i]+OprogStatus.epsdNL;
  maxddot = 0.0;
  for (nn = 0; nn < 6; nn++)
    {
#ifdef MD_ASYM_ITENS
      maxddoti[nn] = calc_maxddot_nnl(i, gradplane_all[nn]);
#else
      maxddoti[nn] = fabs(vx[i]*gradplane_all[nn][0]+vy[i]*gradplane_all[nn][1]+vz[i]*gradplane_all[nn][2])+
	sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
      if (nn==0 || maxddoti[nn] > maxddot)
	maxddot = maxddoti[nn];
      //printf("nn=%d maxddoti=%.15G\n", nn, maxddoti);
    }
  h = OprogStatus.h; /* last resort time increment */
  delt = h;
  
  if (search_contact_faster_neigh_plane_all(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, 
					       dists, maxddoti, maxddot))
    {
      return 0;  
    }
  assign_vec(vecgd, vecgdold);/* assegna a vecgdold vecgd */
  timesSNL++;
  foundrc = 0;
  assign_dists(dists, distsOld);
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
      d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgd, 0);
#else
      if (!firstaftsf)
	{
	  deldist = get_max_deldist(distsOld2, distsOld);
	  normddot = fabs(deldist)/delt;
	  /* NOTA: forse qui si potrebbe anche usare sempre delt = epsd/maxddot */
	  if (normddot!=0)
	    delt = epsd/normddot;
	  else
	    delt = epsd/maxddot;
	  //  delt = h;
	  if (fabs(dold) < epsd)
	    delt = epsd / maxddot;
	}
      else
	{
	  delt = h;//EPS*fabs(t);
	  firstaftsf = 0;
	  dold2 = calcDistNegNeighPlaneAll(t-delt, t1, i, distsOld2, vecgdold2, 0);
	  continue;
	}
      tini = t;
      t += delt;
      d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgd, 0);
      deldist = get_max_deldist(distsOld, dists);
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
	  d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgdold2, 0);
	  assign_vec(vecgdold2, vecgd);
	}
#endif
      for (nn=0; nn < 6; nn++)
	dorefine[nn] = 0;
      ncr=check_cross(distsOld, dists, crossed);
      for (nn = 0; nn < 6; nn++)
	{
	  t2arr[nn] = t + t1; 
	  dorefine[nn] = crossed[nn];
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
      assign_vec(vecgd, vecgroot);
#if 1
      ntc = get_dists_tocheck(distsOld, dists, tocheck, dorefine);
		
      //assign_vec(vecgd, vecgroot);
      for (nn = 0; nn < 6; nn++)
	{
	  if (tocheck[nn])
	    {
	      for (kk=0; kk < 8; kk++)
    		vecgroot[nn][kk] = vecgd[nn][kk];
	      if (interpolNeighPlane(i, t1, t-delt, delt, distsOld[nn], dists[nn], 
				     &troot, vecgroot[nn],  1, nn))
		{
		  dorefine[nn] = 0;
		}
	      else 
		{
	       	  dorefine[nn] = 1;
      		  t2arr[nn] = troot;
		}
	    }
	  else if (dorefine[nn])
	    {
  	      for (kk=0; kk < 8; kk++)
		vecgroot[nn][kk] = vecgd[nn][kk];
	      if (interpolNeighPlane(i, t1, t-delt, delt, distsOld[nn], dists[nn], 
				     &troot, vecgroot[nn],  0, nn))
		{
		  for (kk=0; kk < 8; kk++)
		    vecgroot[nn][kk] = vecgdold[nn][kk];
		  t2arr[nn] = t1 + t - delt;
		}
	      else
		{
  		  t2arr[nn] = troot;
		}
	    }
	}
#endif
      gotcoll = 0;
      firstev = 1;
      for (nn = 0; nn < 6; nn++)
	{
	  if (dorefine[nn]!=0)
	    {
	      assign_plane(nn);
	      //printf("nn=%d dists[%d]: %.15G distsOld[%d]:%.15G\n", nn, nn, dists[nn], nn, distsOld[nn])
	      if (refine_contact_neigh_plane(i, t1, t2arr[nn], vecgroot[nn], vecg, nn))
		{
		  //printf("[locate_contact] Adding collision for ellips. N. %d t=%.15G t1=%.15G t2=%.15G\n", i,
		//	 vecg[4], t1 , t2);
		  MD_DEBUG30(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
		  MD_DEBUG30(printf("[locate_contact] t=%.15G nn=%d\n", t, nn));
		  MD_DEBUG(printf("[locate_contact] its: %d\n", its));
		  /* se il legame già c'è e con l'urto si forma tale legame allora
		   * scarta tale urto */
		  if (vecg[4] > t2 || vecg[4] < t1)
		    {
		      continue;
		    }
		  else
		    {
		      gotcoll = 1;

		      if (firstev || vecg[4] < *evtime)
			{
			  firstev = 0;
			  *evtime = vecg[4];
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
		  MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
#ifdef MD_INTERPOL
		  if (!tocheck[nn])
#endif
		  mdPrintf(ALL,"[locate_contact_nnl] can't find contact point!\n",NULL);
		  /* Se refine_contact fallisce deve cmq continuare a cercare 
		   * non ha senso smettere...almeno credo */
		  //gotcoll = -1;
		  continue;
		}
	    }
	}
      if (gotcoll == 1)
	return 1;
      if (fabs(d) > epsdFastR)
	{
	  if (search_contact_faster_neigh_plane_all(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, 
						       dists, maxddoti, maxddot))
	    {
	      MD_DEBUG30(printf("[search contact faster locate_contact] d: %.15G\n", d));
	      return 0;
	    }
	  dold = d;
	  assign_dists(dists, distsOld);
	  assign_vec(vecgd, vecgdold);
	  firstaftsf = 1;
	  its++;
	  continue;
	}
      dold = d;
      assign_vec(vecgd, vecgdold);
#ifndef MD_BASIC_DT
      assign_dists(distsOld,  distsOld2);
#endif
      assign_dists(dists, distsOld);
      its++;
      itsSNL++;
    }
  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return 0;
}
int dist_too_big(int i, double t, double t1)
{
  double DR[3], rA[3], ti;
  int kk;
  ti = t + (t1 - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  for (kk = 0; kk < 3; kk++)
    {
      DR[kk] = rA[kk] - rB[kk];
    }
  if (calc_norm(DR) > 0.5*maxax[i])
    return 1;
  else 
    return 0;
}
int locate_contact_neigh_plane(int i, double vecg[5], int nplane, double tsup)
{
  double h, d, dold, vecgd[8], vecgdold[8], t, r1[3], r2[3]; 
  double dtmp, t1, t2, maxddot, delt, troot, vecgroot[8];
  //const int MAXOPTITS = 4;
#ifndef MD_BASIC_DT
  double ddot[3], dold2, vecgdold2[8], normddot;
#endif
#ifndef MD_ASYM_ITENS
  double factori;
#endif
  double epsd, epsdFast, epsdFastR, epsdMax; 
  int dorefine, distfail;
  int its, foundrc, kk;
  epsd = OprogStatus.epsdNL;
  epsdFast = OprogStatus.epsdFastNL;
  epsdFastR= OprogStatus.epsdFastRNL;
  epsdMax = OprogStatus.epsdMaxNL;
  /* NOTA: implementare le varie funzioni _neigh (search_contact_faster_neigh, ecc.)
   * in tali funzioni la particella j non è altro che un ellissoide più grande di i
   * con lo stesso centro e immobile */
  t = 0.0;//Oparams.time;
  timesSNL++;
  calc_grad_and_point_plane(i, gradplane, rB, nplane);
  dtmp = calcDistNegNeighPlane(t, Oparams.time, i, r1, r2, vecgd, 0, 0, &distfail, nplane);
  if (dtmp < 0)
    {
      printf("La distanza fra l'ellissoide N. %d e il piano %d è negativa d=%.15G\n", i, nplane, dtmp);
      printf("nexttime[%d]: %.15G\n", i, nebrTab[i].nexttime);
      exit(-1);
    }
  if (!bracket_neigh(i, &t1, &t2, nplane))
    {
      //printf("NOT BRACK NEIGH\n");
      return 0;
    }
  t1 += Oparams.time;	
  if (tsup < t2)
    t2 = tsup;
  else
    t2 += Oparams.time;
  //printf("LOCATE_CONTACT_NNL nplane=%d grad=%.8f %.8f %.8f  rB=%.8f %.8f %.8f t1=%.8f t2=%.8f tsup=%.8f maxax[%d]=%f\n", nplane, 
  //	 gradplane[0], gradplane[1], gradplane[2], rB[0], rB[1], rB[2], t1, t2, tsup, i, maxax[i]);
#ifdef MD_ASYM_ITENS
  maxddot = calc_maxddot_nnl(i, gradplane);
#else
  factori = 0.5*maxax[i]+OprogStatus.epsdNL;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
  maxddot = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
  h = OprogStatus.h; /* last resort time increment */
  if (search_contact_faster_neigh_plane(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2, nplane))
    return 0;  
  MD_DEBUG(printf(">>>>d:%f\n", d));
  foundrc = 0;
  for (kk = 0; kk < 8; kk++)
    vecgdold[kk] = vecgd[kk];
  dold = d;
  its = 0;
  while (t + t1 < t2 || !dist_too_big(i,t,t1))
    {
      MD_DEBUG31(printf("LOC CONT rB = %.15G %.15G %.15G\n", rB[0], rB[1], rB[2]));
#ifdef MD_BASIC_DT
      delt = epsd / maxddot;
      t += delt;
      d = calcDistNegNeighPlane(t, t1, i, r1, r2, vecgd, 0, 0, &distfail, nplane);
#else
      normddot = calcvecFNeigh(i, t, t1, ddot, r1);
      if (normddot!=0)
	delt = epsd / normddot;
      else
	delt = epsd / maxddot;
	//delt = h;
      if (dold < epsd)
	delt = epsd / maxddot;
      t += delt;
      //printf("normddot=%f dt=%.15G\n",normddot, epsd/normddot); 
      for (kk = 0; kk < 8; kk++)
	vecgdold2[kk] = vecgd[kk];
      dold2 = dold;
      //printf("NNL [LOCATE_CONTACT] >>>>>>>>>>>><<<<<<<<<<\n"); 
      d = calcDistNegNeighPlane(t, t1, i, r1, r2, vecgd, 0, 0, &distfail, nplane);
      MD_DEBUG31(printf("NNL [LOCATE_CONTACT] delt=%.15G (%.15G,maxddot=%.10G,normddot=%.15G) t=%.15G ellips N. %d d=%.15G dold=%.15G its=%d t1=%.15G t2=%.15G vparall=%.15G\n", 
      delt, epsd/maxddot, maxddot, normddot, t, i, d, dold, its, t1, t2, fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2]))); 
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
	  itsSNL++;
	  d = calcDistNegNeighPlane(t, t1, i, r1, r2, vecgdold2, 0, 0, &distfail, nplane);
	  for (kk = 0; kk < 8; kk++)
	    vecgd[kk] = vecgdold2[kk];
	  //printf("D delt: %.15G d2-d2o:%.15G d2:%.15G d2o:%.15G\n", delt*epsd/fabs(d2-d2old), fabs(d2-d2old), d2, d2old);
	}
#endif
#if 1
      if (d > epsdFastR)
	{
	  if (search_contact_faster_neigh_plane(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2, nplane))
	    {
	      MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
	      return 0;
	    }
	  for (kk = 0; kk < 8; kk++)
	    vecgdold[kk] = vecgd[kk];
	  dold = d;
	  its++;
	  itsSNL++;
	  continue;
	}
#endif
      MD_DEBUG(printf(">>>> t = %f d1:%f d2:%f d1-d2:%.15G\n", t, d1, d2, fabs(d1-d2)));
      dorefine = 0;      
      if (dold > 0 && d < 0)
	{
	  //printf(">>>>>>>>>QUI\n");
       	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
#ifndef MD_NOINTERPOL  
	 if (interpolNeighPlane(i, t1, t-delt, delt, dold, d, &troot, vecgroot, 0, nplane))
#endif
	    {
	      /* vecgd2 è vecgd al tempo t-delt */
	      for (kk=0; kk < 8; kk++)
		vecgroot[kk] = vecgdold[kk];
	      troot = t1 + t - delt;
	    }
	  dorefine = 1;
	}
      else if (dold < OprogStatus.epsdNL && d < OprogStatus.epsdNL)
	{
#ifndef MD_NOINTERPOL
	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
	  
	  if (interpolNeighPlane(i, t1, t-delt, delt, dold, d, &troot, vecgroot, 1, nplane))
	    dorefine = 0;
	  else 
	    dorefine = 1;
#endif
	}
      if (dorefine)
	{
	  MD_DEBUG31(printf("NNL REFINING CONTACT t=%.15G troot=%.15G t1=%.15G\n", t, troot, t1));
	  if (refine_contact_neigh_plane(i, t1, troot, vecgroot, vecg, nplane))
	    {
	      MD_DEBUG31(printf("[locate_contact] Adding collision between %d\n", i));
	      MD_DEBUG31(printf("collision will occur at time %.15G\n", vecg[4])); 
	      MD_DEBUG31(printf("[locate_contact] its: %d\n", its));
	      if (vecg[4]>t2 || vecg[4]<t1)
		return 0;
	      else
		{
#if 0
		  if (i==2)
		    printf("[i=2] nexttime = %.15G\n", vecg[4]);
#endif		    
		      
		  return 1;
		}
	    }
	  else 
	    {
	      MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
	      if (d < 0)
		{
		  MD_DEBUG10(printf("t=%.15G d2 < 0 and I did not find contact point, boh...\n",t));
		  MD_DEBUG10(printf("d1: %.15G d2: %.15G\n", d1, d2));
		  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
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
      dold = d;
      for (kk = 0; kk < 8; kk++)
	vecgdold[kk] = vecgd[kk];
      its++;
      itsSNL++;
    }
  MD_DEBUG(  
  if (foundrc==0)
    printf("%d-%d t=%.12G > t2=%.12G I did not find any contact point!\n", i, j, t, t2);
  );

  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return foundrc;
}
int search_contact_faster_neigh(int i, double *t, double t1, double t2, 
				double *vecgd, double epsd, double *d1, double epsdFast,
				double *r1, double *r2)
{
  /* NOTA: 
   * MAXOPTITS è il numero massimo di iterazioni al di sopra del quale esce */
  double maxddot, told, delt, normddot, ddot[3];
  const int MAXOPTITS = 500;
  double factori;
  int its=0, distfailed, itsf=0; 
  const double GOLD= 1.618034;  
  factori = 0.5*maxax[i]+OprogStatus.epsdNL;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));

  /* estimate of maximum rate of change for d */
  maxddot = sqrt(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
  *d1 = calcDistNegNeigh(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed);
  timesFNL++;
  MD_DEBUG20(printf("Pri distances between %d d1=%.12G epsd*epsdTimes:%f\n", i, *d1, epsdFast));
  printf("[SEARCH_CONTACT_FASTER] t=%.15G ellips N. %d d=%.15G\n", *t, i, *d1); 
  while (*d1 > epsdFast && its < MAXOPTITS)
    {
      told = *t;
      delt = *d1 / maxddot;
      MD_DEBUG20(printf("SEARCH_CONTACT_FASTER delt=%.15G maxddot: %.15G\n", delt, maxddot));
      normddot = calcvecFNeigh(i, *t, t1, ddot, r1);
      //printf("normddot: %.15G\n", epsd/normddot);
      /* check for convergence */
     
      if (normddot!=0 && delt < (epsd / normddot))
	{
	  MD_DEBUG20(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
      *t += delt;
#if 1
      if (*t + t1 > t2)
	{
	  *t = told;
	  MD_DEBUG20(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG20(printf("convergence t>t2\n"));
	  *d1 = calcDistNegNeigh(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed);
	  return 1;
	}
#endif
      *d1 = calcDistNegNeigh(*t, t1, i, r1, r2, vecgd, 1, 1, &distfailed);
      //printf("NNL LOOP SEARCH CONTACT FASTER i=%d *d1=%.15G *t=%.15G\n", i, *d1, *t+t1);
      /* NOTA: nel caso di urto di un ellissoide con la sua neighbour list 
       * a t=0 i due ellissoidi hanno stesso centro e assi principali paralleli
       * e questo implica che ci possono essere più soluzioni possibili per le eq.
       * che definiscono la distanza.
       * Questo implica che la distanza che ottengo è una sovrastima di quella vera e quindi 
       * con il loop che segue compenso questo problema riducendo il passo fino a che 
       * non arrivo ad una distanza positiva con il passo veloce */
#if 1
      itsf = 0;
      while (*d1 < 0 || distfailed)
	{
	  /* reduce step size */
	  delt /= GOLD;
	  *t = told + delt;
	  *d1 = calcDistNegNeigh(*t, t1, i, r1, r2, vecgd, 1, 1, &distfailed);
	  printf("itsf=%d BUBU SEARCH_CONTACT_FASTER_NEIGH *d1=%.15G\n",itsf,*d1);
	  itsf++;	
	  if (itsf > 100)
	    {
	      printf("*d1=%.15G too many times calculation of distance failed!\n", *d1);
	      printf("aborting...\n");
	      exit(-1);
	    }
	}
#else
      if (*d1 < 0)
	{
	  /* go back! */
	  MD_DEBUG20(printf("d1<0 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG20(printf("d1 negative in %d iterations d1= %.15f\n", its, *d1));
	  *t = told;	  
	  *d1 = calcDistNegNeigh(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed);
	  return 0;
	}
#endif
      told = *t;
      its++;
      itsFNL++;
    }

  MD_DEBUG20(printf("max iterations %d iterations reached t=%f t2=%f\n", its, *t, t2));
  return 0;

}
int refine_contact_neigh(int i, double t1, double t, double vecgd[8], double  vecg[5])
{
  int kk, retcheck;

  for (kk = 0; kk < 3; kk++)
    vecg[kk] = (vecgd[kk]+vecgd[kk+3])*0.5; 
    
  vecg[3] = vecgd[6];
  vecg[4] = t-t1;
  trefG = t1;
  newtNeigh(vecg, 5, &retcheck, funcs2beZeroedNeigh, i); 
  vecg[4] += t1;
  if (retcheck==2)
    {
      MD_DEBUG10(printf("newtNeigh did not find any contact point!\n"));
      return 0;
    }
  else
    {
      return 1; 
    }
}
int locate_contact_neigh(int i, double vecg[5])
{
  double h, d, dold, dold2, vecgdold2[8], vecgd[8], vecgdold[8], t, r1[3], r2[3]; 
  double normddot, ddot[3], t1, t2, maxddot, delt, troot, vecgroot[8];
  //const int MAXOPTITS = 4;
  double epsd, epsdFast, epsdFastR, epsdMax, factori; 
  int dorefine, distfail;
  int its, foundrc, kk;
  epsd = OprogStatus.epsdNL;
  epsdFast = OprogStatus.epsdFastNL;
  epsdFastR= OprogStatus.epsdFastRNL;
  epsdMax = OprogStatus.epsdMaxNL;
  /* NOTA: implementare le varie funzioni _neigh (search_contact_faster_neigh, ecc.)
   * in tali funzioni la particella j non è altro che un ellissoide più grande di i
   * con lo stesso centro e immobile */
  t = 0.0;//Oparams.time;
  t1 = Oparams.time;	
  t2 = timbig;
  
  factori = 0.5*maxax[i]+OprogStatus.epsdNL;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
  maxddot = sqrt(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
  h = OprogStatus.h; /* last resort time increment */
  t += h;
  if (search_contact_faster_neigh(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2))
    return 0;  
  MD_DEBUG(printf(">>>>d:%f\n", d));
  foundrc = 0;
  for (kk = 0; kk < 8; kk++)
    vecgdold[kk] = vecgd[kk];
  dold = d;
  its = 0;
  while (t + t1 < t2)
    {
      normddot = calcvecFNeigh(i, t, t1, ddot, r1);
      if (normddot!=0)
	delt = epsd / normddot;
      else
	delt = h;
      if (dold < epsd)
	delt = epsd / maxddot;
      t += delt;
      //printf("normddot=%f dt=%.15G\n",normddot, epsd/normddot); 
      for (kk = 0; kk < 8; kk++)
	vecgdold2[kk] = vecgd[kk];
      dold2 = dold;
      printf("[LOCATE_CONTACT] >>>>>>>>>>>><<<<<<<<<<\n"); 
      d = calcDistNegNeigh(t, t1, i, r1, r2, vecgd, 0, 0, &distfail);
      printf("[LOCATE_CONTACT] t=%.15G ellips N. %d d=%.15G dold=%.15G its=%lld\n", t, i, d, dold, itsSNL); 
      if (fabs(d-dold2) > epsdMax)
	{
	  /* se la variazione di d è eccessiva 
	   * cerca di correggere il passo per ottenere un valore
	   * più vicino a epsd*/
	  //printf("P delt: %.15G d2-d2o:%.15G d2:%.15G d2o:%.15G\n", delt, fabs(d2-d2old), d2, d2old);
	  t -= delt;
	  //delt = d2old / maxddot;
	  delt = epsd / maxddot;
	  if (delt < h)
	    delt = h;
	  t += delt; 
	  //t += delt*epsd/fabs(d2-d2old);
	  itsSNL++;
	  d = calcDistNegNeigh(t, t1, i, r1, r2, vecgdold2, 0, 0, &distfail);
	  for (kk = 0; kk < 8; kk++)
	    vecgd[kk] = vecgdold2[kk];
	  //printf("D delt: %.15G d2-d2o:%.15G d2:%.15G d2o:%.15G\n", delt*epsd/fabs(d2-d2old), fabs(d2-d2old), d2, d2old);
	}
#if 1
      if (d > epsdFastR)
	{
	  if (search_contact_faster_neigh(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2))
	    {
	      MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
	      return 0;
	    }
	  for (kk = 0; kk < 8; kk++)
	    vecgdold[kk] = vecgd[kk];
	  dold = d;
	  its++;
	  //itsSNL++;
	  continue;
	}
#endif
      MD_DEBUG(printf(">>>> t = %f d1:%f d2:%f d1-d2:%.15G\n", t, d1, d2, fabs(d1-d2)));
      dorefine = 0;      
      if (dold > 0 && d < 0)
	{
	  printf(">>>>>>>>>QUI\n");
       	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
#ifndef MD_NOINTERPOL  
	 if (interpolNeigh(i, t1, t-delt, delt, dold, d, &troot, vecgroot, 0))
#endif
	    {
	      /* vecgd2 è vecgd al tempo t-delt */
	      for (kk=0; kk < 8; kk++)
		vecgroot[kk] = vecgdold[kk];
	      troot = t1 + t - delt;
	    }
	  dorefine = 1;
	}
      else if (dold < OprogStatus.epsdNL && d < OprogStatus.epsdNL)
	{
#ifndef MD_NOINTERPOL
	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
	  
	  if (interpolNeigh(i, t1, t-delt, delt, dold, d, &troot, vecgroot, 1))
	    dorefine = 0;
	  else 
	    dorefine = 1;
#endif
	}
      if (dorefine)
	{
	  printf("REFINING CONTACT t=%.15G\n", t);
	  if (refine_contact_neigh(i, t1, troot, vecgroot, vecg))
	    {
	      MD_DEBUG20(printf("[locate_contact] Adding collision between %d\n", i));
	      MD_DEBUG20(printf("collision will occur at time %.15G\n", vecg[4])); 
	      MD_DEBUG20(printf("[locate_contact] its: %d\n", its));
	      if (vecg[4]>t2 || vecg[4]<t1)
		return 0;
	      else
		return 1;
	    }
	  else 
	    {
	      MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
	      if (d < 0)
		{
		  MD_DEBUG10(printf("t=%.15G d2 < 0 and I did not find contact point, boh...\n",t));
		  MD_DEBUG10(printf("d1: %.15G d2: %.15G\n", d1, d2));
		  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
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
      dold = d;
      for (kk = 0; kk < 8; kk++)
	vecgdold[kk] = vecgd[kk];
      its++;
      itsSNL++;
    }
  MD_DEBUG(  
  if (foundrc==0)
    printf("%d-%d t=%.12G > t2=%.12G I did not find any contact point!\n", i, j, t, t2);
  );

  MD_DEBUG10(printf("[locate_contact] its: %d\n", its));
  return foundrc;
}
extern double max(double a, double b);
#ifdef MD_PATCHY_HE
extern int locate_contactSP(int i, int j, double shift[3], double t1, double t2, double *evtime, int *ata, int *atb, int *collCode);
extern void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent);
#endif
void PredictEventNNL(int na, int nb) 
{
  int i, signDir[NDIM]={0,0,0}, evCode, k, n;
  double vecg[5], shift[3], t1, t2, t, tm[NDIM];
  double sigSq, tInt, d, b, vv, dv[3], dr[3], distSq;
  int overlap;
#ifdef MD_PATCHY_HE
  int ac, bc, collCode, collCodeOld, acHC, bcHC;
  double evtime, evtimeHC;
#endif
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
#if 1
  if (tm[k]<0)
    {
      tm[k] = 0.0;
#if 1
      printf("tm[%d]<0 step %lld na=%d\n", k, (long long int)Oparams.curStep, na);
      printf("Cells(%d,%d,%d)\n", inCell[0][na], inCell[1][na], inCell[2][na]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
#endif
    }
#endif
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 
  ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]);
  /* NOTA: nel caso di attraversamento di una cella non deve predire le collisioni */
  if (nb >= ATOM_LIMIT)
    return;
  MD_DEBUG32(printf("nebrTab[%d].len=%d\n", na, nebrTab[na].len));
  for (i=0; i < nebrTab[na].len; i++)
    {
      n = nebrTab[na].list[i]; 
#if 0
      if (na==35 || na==16)
	printf("[PredictEventNNL] na=%d n=%d nb=%d\n",na,  n, nb);
#endif
      if (!(n != na && n!=nb && (nb >= -1 || n < na)))
	continue;
     // for (kk=0; kk < 3; kk++)
     //	shift[kk] = nebrTab[na].shift[i][kk];
      shift[0] = L*rint((rx[na]-rx[n])/L);
      shift[1] = L*rint((ry[na]-ry[n])/L);
      shift[2] = L*rint((rz[na]-rz[n])/L);

      /* maxax[...] è il diametro dei centroidi dei due tipi
       * di ellissoidi */
      if (OprogStatus.targetPhi > 0)
	{
	  sigSq = Sqr(max_ax(na)+max_ax(n)+OprogStatus.epsd);
	}
      else
	{
	  if (na < parnumA && n < parnumA)
	    sigSq = Sqr(maxax[na]+OprogStatus.epsd);
	  else if (na >= parnumA && n >= parnumA)
	    sigSq = Sqr(maxax[na]+OprogStatus.epsd);
	  else
	    sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
	}
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
      distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
      d = Sqr (b) - vv * (distSq - sigSq);
      if (d < 0 || (b > 0.0 && distSq > sigSq)) 
	{
	  /* i centroidi non collidono per cui non ci può essere
	   * nessun urto sotto tali condizioni */
	  continue;
	}
      MD_DEBUG32(printf("PREDICTING na=%d n=%d\n", na , n));
      if (vv==0.0)
	{
	  if (distSq >= sigSq)
	    {
	      continue;
	    }
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
      tnnl = min(nebrTab[na].nexttime,nebrTab[n].nexttime);
      if (tnnl < t2)
	t2 = tnnl;
#else      
      /* WARNING: OprogStatus.h è un buffer di sicurezza */
      if (nextNNLrebuild < t2)
	t2 = nextNNLrebuild + OprogStatus.h;
#endif
      // t1 = Oparams.time;
      // t2 = nebrTab[na].nexttime;//,nebrTab[n].nexttime);

      MD_DEBUG32(printf("nexttime[%d]:%.15G\n", n, nebrTab[n].nexttime));
      MD_DEBUG32(printf("locating contact between %d and %d t1=%.15G t2=%.15G\n", na, n, t1, t2));
#ifdef MD_PATCHY_HE
      evtime = t2;
      collCode = MD_EVENT_NONE;
      rxC = ryC = rzC = 0.0;
      MD_DEBUG31(printf("t1=%.15G t2=%.15G\n", t1, t2));
      collCodeOld = collCode;
      evtimeHC = evtime;
      acHC = ac = 0;
      bcHC = bc = 0;
      if (OprogStatus.targetPhi <=0 && ((na < Oparams.parnumA && n >= Oparams.parnumA)|| 
					(na >= Oparams.parnumA && n < Oparams.parnumA)))
	{
	  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
	    {
	      collCode = MD_EVENT_NONE;
	    }
	}
      if (collCode!=MD_EVENT_NONE)
	t2 = evtime+1E-7;
      if (locate_contact(na, n, shift, t1, t2, vecg))
	{
	  if (collCode == MD_EVENT_NONE || (collCode!=MD_EVENT_NONE && vecg[4] <= evtime))
	    {
	      collCode = MD_CORE_BARRIER;
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
#else
      if (!locate_contact(na, n, shift, t1, t2, vecg))
	{
	  continue;
	}
      rxC = vecg[0];
      ryC = vecg[1];
      rzC = vecg[2];
      t = vecg[4];
#endif
      MD_DEBUG32(printf("Scheduling collision between %d and %d at t=%.15G\n", na, n, t));
#ifdef MD_PATCHY_HE
      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
#else
      ScheduleEvent (na, n, t);
#endif
    }
}
#ifdef MD_PATCHY_HE
extern int locate_contact_neigh_plane_parall_sp(int i, double *evtime, double t2);
#endif
void updrebuildNNL(int na)
{
  /* qui ricalcola solo il tempo di collisione dell'ellisoide na-esimo con 
   * la sua neighbour list */
  double vecg[5];
  double nnltime1, nnltime2;
  int ip;
#ifdef MD_NNLPLANES
#ifdef MD_PATCHY_HE
  if (OprogStatus.targetPhi <= 0.0)
    {
      if (!locate_contact_neigh_plane_parall_sp(na, &nnltime1, timbig))
	{
	  printf("[ERROR] failed to find escape time for sticky spots\n");
	  exit(-1);
	}
    }
  else 
    nnltime1 = timbig;
  MD_DEBUG32(printf("sptime: %.15G nexttime=%.15G\n", sptime, nebrTab[na].nexttime));
#else
  nnltime1 = timbig;
#endif
  nnltime2 = timbig;
  if (OprogStatus.paralNNL)
    {
      if (!locate_contact_neigh_plane_parall(na, &nnltime2, nnltime1))
	{
#ifndef MD_PATCHY_HE
	  printf("[ERROR] failed to find escape time for ellipsoid N. %d\n", na);
	  exit(-1);
#endif
	}
    }
 else
   {
     nnltime2 = nnltime1;
     for (ip = 0; ip < 6; ip++)
       {
	 if (!locate_contact_neigh_plane(na, vecg, ip, nnltime1))
	   continue;
	 if (vecg[4] < nnltime2)
	   nnltime2 = vecg[4];
       }
   }
 nebrTab[na].nexttime = min(nnltime1, nnltime2);
#else
  if (!locate_contact_neigh(na, vecg))
    nebrTab[na].nexttime = timbig;
  else
    nebrTab[na].nexttime = vecg[4];
#endif
#if 0
    {
      double gradA[3], gradB[3], sp;
      calc_grad(vecg, rA, Xa, gradA);
      calc_grad(vecg, rB, Xb, gradB);
      sp = scalProd(gradA, gradB);
      if (sp < 0)
	{
	  printf("BOH i prodotti scalari sono diretti in maniera opposta!\n");
	  exit(-1);
	}
      else 
	printf("sp > 0!!!!!!!!!\n");
    }
#endif
 
  //printf("updneigh REBUILD i=%d t=%.15G\n", na, nebrTab[na].nexttime);
  if (nebrTab[na].nexttime < nextNNLrebuild)
    nextNNLrebuild = nebrTab[na].nexttime;
}
void updAllNNL()
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    updrebuildNNL(i);
}
void nextNNLupdate(int na)
{
  int i1, i2, ip;
  double DelDist, nnlfact;
  const double distBuf = 0.1;
  double vecg[5];
  double nnltime1, nnltime2;
#ifdef MD_ASYM_ITENS
  double psi, phi;
#else
  double Omega[3][3];
#endif
#ifndef MD_NNLPLANES
  nebrTab[na].axa = OprogStatus.rNebrShell*axa[na];
  nebrTab[na].axb = OprogStatus.rNebrShell*axb[na];
  nebrTab[na].axc = OprogStatus.rNebrShell*axc[na];
  DelDist = max3(nebrTab[na].axa,nebrTab[na].axb,nebrTab[na].axc) -
    max3(axa[na],axb[na],axc[na]);
#else
  if (OprogStatus.targetPhi > 0.0)
    {
      nnlfact = axa[na]/Oparams.a[na<Oparams.parnumA?0:1];
      nebrTab[na].axa = nnlfact*OprogStatus.rNebrShell+axa[na];
      nebrTab[na].axb = nnlfact*OprogStatus.rNebrShell+axb[na];
      nebrTab[na].axc = nnlfact*OprogStatus.rNebrShell+axc[na];
      DelDist = nnlfact*OprogStatus.rNebrShell;
    }
  else
    {
      nebrTab[na].axa = OprogStatus.rNebrShell+axa[na];
      nebrTab[na].axb = OprogStatus.rNebrShell+axb[na];
      nebrTab[na].axc = OprogStatus.rNebrShell+axc[na];
      DelDist = OprogStatus.rNebrShell;
    }
#endif
  DelDist += distBuf;
  MD_DEBUG31(printf("DelDist=%.15G\n", DelDist));
  nebrTab[na].r[0] = rx[na];
  nebrTab[na].r[1] = ry[na];
  nebrTab[na].r[2] = rz[na];
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(na, 0, RtB, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(na, 0, RtB, Omega);
#endif
  for (i1 = 0; i1 < 3; i1++)
    for (i2 = 0; i2 < 3; i2++)
      nebrTab[na].R[i1][i2] = RtB[i1][i2];
  /* calcola il tempo a cui si deve ricostruire la NNL */
  MD_DEBUG31(printf("BUILDING NNL FOR i=%d\n",na));
#ifdef MD_NNLPLANES
#ifdef MD_PATCHY_HE
  if (OprogStatus.targetPhi <= 0.0)
    {
      if (!locate_contact_neigh_plane_parall_sp(na, &nnltime1, timbig))
	{
	  printf("[ERROR] failed to find escape time for sticky spots\n");
	  exit(-1);
	}
    }
  else
    nnltime1 = timbig;
  MD_DEBUG32(printf("[nextNNLupdate] sptime: %.15G nexttime=%.15G\n", sptime, nebrTab[na].nexttime));
#else
  nnltime1 = timbig; 
#endif
  nnltime2 = timbig;
  if (OprogStatus.paralNNL)
    {
      if (!locate_contact_neigh_plane_parall(na, &nnltime2, nnltime1))
	{
#ifndef MD_PATCHY_HE
	  printf("[ERROR] failed to find escape time for ellipsoid N. %d\n", na);
	  exit(-1);
#endif
	}
    }
  else
    {
      for (ip = 0; ip < 6; ip++)
       	{
 	  if (!locate_contact_neigh_plane(na, vecg, ip, nnltime1))
 	    continue;
 	  if (vecg[4] < nnltime2)
 	    nnltime2 = vecg[4];
  	}
    }
  nebrTab[na].nexttime = min(nnltime1, nnltime2);
  //printf(">> nexttime=%.15G\n", nebrTab[na].nexttime);
#else
  if (!locate_contact_neigh(na, vecg))
    nebrTab[na].nexttime = timbig;
  else
    nebrTab[na].nexttime = vecg[4];
#endif
#if 0
    {
      double gradA[3], gradB[3], sp;
      calc_grad(vecg, rA, Xa, gradA);
      calc_grad(vecg, rB, Xb, gradB);
      sp = scalProd(gradA, gradB);
      if (sp < 0)
	{
	  printf("[BUILDNNL] BOH i prodotti scalari sono diretti in maniera opposta!\n");
	  exit(-1);
	}
      else 
	printf("[BUILD NNL] sp > 0!!!!!!!!!\n");
    }
#endif
  
#if 0
  factori = 0.5*maxax[na];
  maxddot = sqrt(Sqr(vx[na])+Sqr(vy[na])+Sqr(vz[na])) +
    sqrt(Sqr(wx[na])+Sqr(wy[na])+Sqr(wz[na]))*factori;
  nebrTab[na].nexttime = Oparams.time+0.5*OprogStatus.rNebrShell/maxddot;
  //ScheduleEvent(na, ATOM_LIMIT + 11, vecg[4]); 
#endif
  nebrTab[na].len=0;
  nebrTab[na].time = Oparams.time;
}
void BuildNNL(int na) 
{
  double shift[NDIM];
  int kk;
  double dist;
#ifndef MD_NNLPLANES
  double vecgsup[8], alpha;
#endif
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cels[NDIM];*/
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
  nebrTab[na].len = 0;
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (kk=0; kk < 3; kk++)
    shift[kk] = 0;
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodic boundary condition along z if gravitational
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
		  if (n != na)// && n != nb && (nb >= -1 || n < na)) 
		    {
      		      //dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
#ifdef MD_NNLPLANES
		      dist = calcDistNegNNLoverlapPlane(Oparams.time, 0.0, na, n, shift); 
#else
		      dist = calcDistNegNNLoverlap(Oparams.time, 0.0, na, n, shift,
						   r1, r2, &alpha, vecgsup, 1); 
#endif
		      /* 0.1 è un buffer per evitare problemi, deve essere un parametro 
		       * in OprogStatus */
#if 0
		      if (n==132 && na==106)
			printf("beccato!! dist=%f\n", calcDistNegNNLoverlapPlane(Oparams.time, 0.0, 132, 106, shift));
#endif
		      if (dist < 0)
			{
			  MD_DEBUG32(printf("Adding ellipsoid N. %d to NNL of %d\n", n, na));
#if 0
			  if (n==132 && na==106)
			    { printf("bahbah (106-132): %f\n", calcDistNegNNLoverlapPlane(Oparams.time, 0.0, 106, 132, shift)); 

			      printf("bahbah (132-106): %f\n", calcDistNegNNLoverlapPlane(Oparams.time, 0.0, 132, 106, shift)); 
			      printf("shift= (%f,%f,%f)\n", shift[0], shift[1], shift[2]);
}
#endif
			  nebrTab[na].list[nebrTab[na].len] = n;
			  //for (kk=0; kk < 3; kk++)
			    //nebrTab[na].shift[nebrTab[na].len][kk] = shift[kk];
			  nebrTab[na].len++;
			}
		    }
		} 
	    }
	}
    }
}
