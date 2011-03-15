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
#define MD_DEBUG33(x) 
#define MD_DEBUG34(x) 
#define MD_DEBUG35(x)  
#define MD_DEBUG36(x) 
#define MD_DEBUG37(x) 
#define MD_DEBUG38(x) 
#define MD_DEBUG39(x) 
#ifdef MD_ABSORPTION
extern int *listtmp;
#endif
#ifdef EDHE_FLEX
extern int is_sphere(int i);
extern int *is_a_sphere_NNL;
#endif
#if defined(EDHE_FLEX) || defined(MD_PATCHY_HE)
extern int isSymItens(int i);
#endif
#ifdef MD_SUPERELLIPSOID
int is_superellipse(int i);
#endif
#ifdef EDHE_FLEX
extern int is_infinite_Itens(int i);
extern int is_infinite_mass(int i);
#endif
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
#ifdef MD_HE_PARALL
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
#endif
extern const double timbig;
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
extern double rA[3], rB[3];
extern double rxC, ryC, rzC;
#ifdef MD_LXYZ
extern double pi, invL[3], L2[3], Vz; 
#else
extern double pi, invL, L2, Vz; 
#endif
#ifdef MD_ASYM_ITENS
extern double **Ia, **Ib, **invIa, **invIb;
#else
extern double Ia, Ib, invIa, invIb;
#endif
extern double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
extern int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
  cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
#ifdef MD_EDHEFLEX_OPTNNL
extern int *inCell_NNL[3], *cellList_NNL;
extern double *rxNNL, *ryNNL, *rzNNL;
#endif
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
void BuildNNL(int na);
extern long long int itsF, timesF, itsS, timesS, numcoll, itsFNL, itsSNL, timesSNL, timesFNL;
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern double accngA, accngB;
extern void tRDiagR(int i, double **M, double a, double b, double c, double **Ri);
extern double min(double a, double b);
extern void calcFxtFt(int i, double x[3], double **RM, double cosea[3], double sinea[3], double **X,
	       double D[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft);
extern void calcFxtFtSym(double x[3], double **X,
	       double D[3][3], double Omega[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft);
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
double calc_maxddot_nnl(int i, double *gradplane, double epsd)
{
#if 0
  int na;
  double Iamin;
#endif
  double factori;
  double A, B;
  factori = 0.5*maxax[i]+epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
#if 0
  na = i<Oparams.parnumA?0:1;
  Iamin = min(Oparams.I[na][0],Oparams.I[na][2]);
  return fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
     angM[i]*factori/Iamin;
#else
  A = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2]);
  B = sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
  /* 19/05/2010: se la velocità angolare è nulla (entro gli errori di macchina) 
     allora A è esattamente la velocità d'avvicinamento alla parete e non 
     più una sovrastima per quello moltiplico A per ADJ_MAXDDOT_NL!
     Notare che se vale la seguente condizione A + B = A quindi bisogna correggere A */
  if (B / A < MAXDDOT_NNL_THR)
    return A*ADJ_MAXDDOT_NL + B;
  else
    return A + B;
#endif
}
#endif
#ifdef MD_EDHEFLEX_WALL
extern int globalHW;
extern void calc_grad_and_point_plane_hwbump(int i, double *grad, double *point, int nplane);
extern int locateHardWall(int na, int nplane, double tsup, double vecg[5], int ghw);
#endif
#ifdef MD_ABSORPTION
#ifdef MD_EDHEFLEX_ISOCUBE
void calc_grad_and_point_plane_hwsemiperm(int i, double *grad, double *point, int nplane);
#else
void calc_grad_and_point_plane_hwsemiperm(int i, double *grad, double *point);
#endif
#endif
void calc_grad_and_point_plane(int i, double *grad, double *point, int nplane)
{
  int kk;
  double del=0.0, segno;
#ifdef MD_EDHEFLEX_WALL
  if (globalHW==1)
    {
      calc_grad_and_point_plane_hwbump(i, grad, point, nplane);
      return;
    }
#ifdef MD_ABSORPTION
  if (globalHW==2)
    {
#ifdef MD_EDHEFLEX_ISOCUBE
      calc_grad_and_point_plane_hwsemiperm(i, grad, point, nplane);
#else
      calc_grad_and_point_plane_hwsemiperm(i, grad, point);
#endif
      return;
    }
#endif
#endif
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
  MD_DEBUG39(printf("i=%d lati parallel %f %f %f\n", i, nebrTab[i].axa, nebrTab[i].axb, nebrTab[i].axc));
  /* NOTA: epsdNL+epsd viene usato come buffer per evitare problemi numerici 
   * nell'update delle NNL. */
  /* 19/05/2010 prima usavo OprogStatus.epsd (non epsdFastNL) ma se la SQ o l'HE 
     non ha velocità angolare in locate_contact_neigh_plane_parall_sp() il passo delt=epsd/maxddot 
     porta lo spot esattamente sul piano creando problemi. */
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
  MD_DEBUG39(printf("i=%d del=%f point=%f %f %f grad =%f %f %f\n", i, del, point[0], point[1], point[2], grad[0], grad[1], grad[2] ));
  MD_DEBUG39(printf("nplane=%d pos=%f %f %f\n", nplane, rx[i], ry[i], rz[i]));
}
#ifdef MD_EDHEFLEX_OPTNNL
extern void rebuild_linked_list_NNL();
#endif
void growth_rebuildNNL(int i)
{
  int ii, n;
  nextNNLupdate(i);
#ifdef MD_EDHEFLEX_OPTNNL
  /* il centro di massa del parallelepipedo può non coincidere con quello dell'ellissoide
     in tal caso per cui le linked lists vanno generate ad hoc */
  if (OprogStatus.optnnl)
    rebuild_linked_list_NNL();
#endif
  BuildNNL(i);
  if (nebrTab[i].nexttime < nextNNLrebuild)
    nextNNLrebuild = nebrTab[i].nexttime;
  for (ii=0; ii < nebrTab[i].len; ii++)
    {
      n = nebrTab[i].list[ii]; 
      BuildNNL(n);
    }
}
double n1_bak, n2_bak, n3_bak;
#ifdef EDHE_FLEX
void switch_to_HE(int i)
{
  int t;
  t = typeOfPart[i];
  n1_bak = typesArr[t].n[0];
  n2_bak = typesArr[t].n[1];
  n3_bak = typesArr[t].n[2];
  typesArr[t].n[0] = 2.0;
  typesArr[t].n[1] = 2.0;
  typesArr[t].n[2] = 2.0;
  typesArr[t].sax[0] *= 2.0;
  typesArr[t].sax[1] *= 2.0;
  typesArr[t].sax[2] *= 2.0;
}
void back_to_SQ(int i)
{
  int t;
  t = typeOfPart[i];
  typesArr[t].n[0] = n1_bak;
  typesArr[t].n[1] = n2_bak;
  typesArr[t].n[2] = n3_bak;
  typesArr[t].sax[0] /= 2.0;
  typesArr[t].sax[1] /= 2.0;
  typesArr[t].sax[2] /= 2.0;
}
#endif
#ifdef MD_EDHEFLEX_OPTNNL
extern void rebuild_linked_list_NNL(void);
#endif
#ifdef MD_GHOST_IGG
void update_ghost_status(void)
{
  int i, curnigg, a, notransition, np;

  /* ghostsim=3 means that all IgG are ghost forever (i.e. forever in state 3) */
  if (Oparams.ghostsim==3)
    return;
  np = 0;
  for (a=0; a < 4; a++)
    {
      np += typeNP[a];
    }
  for (i = 0; i < np; i++)
    {
      curnigg = ghostInfoArr[i].iggnum;
      /* 3 is an intermediate state going from 1 to 2, hence
       when all bonds get broken we have following transitions:
       2 -> 3 -> 1 */
      if (ghostInfoArr[i].ghost_status == 3)
	{
	  notransition = 0;
	  for (a = 0; a < nebrTab[i].len; a++)
	    {
	      if ( ghostInfoArr[nebrTab[i].list[a]].iggnum != -1 &&
		   ghostInfoArr[nebrTab[i].list[a]].iggnum != curnigg )
		notransition = 1;
	    }
	  /* if notransition == 1 it means that there are possible overlap 
	     in the transition from 3 (not bound but possibly overlapping) to 1 (bulk) */ 
	  if (notransition == 0)
	    {
	      //printf("transition i=%d 3->1\n", i);
	      ghostInfoArr[i].ghost_status = 1;
	    }
	}
    }
}
#endif
#ifdef MD_MULTIPLE_LL
extern void BuildAllNNL_MLL(void);
#endif
void rebuildNNL(void)
{
  int i;
  double nltime=timbig;
#ifndef MC_SIMUL
  UpdateSystem();
  if (OprogStatus.useNNL==2||OprogStatus.useNNL==4)
    printf("Rebuilding NNL t=%.15G numcoll=%lld\n", Oparams.time, numcoll);
  for (i=0; i < Oparams.parnum; i++)
    {
      //printf("Updating i=%d\n", i);
      nextNNLupdate(i);
      if (i==0 || nebrTab[i].nexttime < nltime)
	nltime = nebrTab[i].nexttime;
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    {
      //printf("Updating i=%d\n", i);
      nextNNLupdate(i);
    }

  if (OprogStatus.useNNL==2||OprogStatus.useNNL==4)
    printf("Rebuilding NNL step=%d\n", Oparams.curStep);
#endif
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
      BuildAllNNL_MLL();
    }
  else
    {
#ifdef MD_EDHEFLEX_OPTNNL
      /* il centro di massa del parallelepipedo può non coincidere con quello dell'ellissoide
	 in tal caso per cui le linked lists vanno generate ad hoc */
      if (OprogStatus.optnnl)
	rebuild_linked_list_NNL();
#endif
      for (i=0; i < Oparams.parnum; i++)
	{
	  BuildNNL(i);
	}
    }
#else
#ifdef MD_EDHEFLEX_OPTNNL
  /* il centro di massa del parallelepipedo può non coincidere con quello dell'ellissoide
     in tal caso per cui le linked lists vanno generate ad hoc */
  if (OprogStatus.optnnl)
    rebuild_linked_list_NNL();
#endif
  for (i=0; i < Oparams.parnum; i++)
    {
      BuildNNL(i);
    }
#endif
#ifdef MD_GHOST_IGG
  /* if ghostsim=1/2 allow transition from 3 to 1 */	
  if (Oparams.ghostsim)
    update_ghost_status();
#endif
 
#ifndef MC_SIMUL
  /* next complete update */
  nextNNLrebuild = nltime;
  if (OprogStatus.useNNL==2||OprogStatus.useNNL==4)
    printf("nextNNLrebuild=%.15G\n", nextNNLrebuild);
#endif
  //ScheduleEvent(-1, ATOM_LIMIT + 11, nltime); 
}
#ifdef MD_SUPERELLIPSOID
extern void funcs2beZeroedNeighPlaneSE(int n, double x[], double fvec[], int i);
extern void fdjacNeighPlaneSE(int n, double x[], double fvec[], double **df, 
		     void (*vecfunc)(int, double [], double [], int), int iA);
#endif

void fdjacNeighPlane(int n, double x[], double fvec[], double **df, 
		     void (*vecfunc)(int, double [], double [], int), int iA)
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  double  rA[3], ti, vA[3], vB[3], OmegaB[3][3];
  double DA[3][3], fx[3], invaSqN, invbSqN, invcSqN;
  double Fxt[3], Ft;
#ifdef EDHE_FLEX
  int typei;
#endif
  double psi, phi;
  double OmegaA[3][3];
  int k1, k2;
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(iA))
    {
      fdjacNeighPlaneSE(n, x, fvec, df, vecfunc, iA);
      return;
    }
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
  if (isSymItens(iA))
    UpdateOrient(iA, ti, RA, OmegaA);
  else
    symtop_evolve_orient(iA, ti, RA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(iA, ti, RA, OmegaA);
#endif
  MD_DEBUG2(printf("i=%d ti=%f", iA, ti));
  MD_DEBUG2(print_matrix(RA, 3));
#ifdef EDHE_FLEX
  typei = typeOfPart[iA];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1/Sqr(axa[iA]);
      invbSqN = 1/Sqr(axb[iA]);
      invcSqN = 1/Sqr(axc[iA]);
    }
  else
    {
      invaSqN = 1/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1/Sqr(axa[iA]);
  invbSqN = 1/Sqr(axb[iA]);
  invcSqN = 1/Sqr(axc[iA]);
#endif
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
  if (isSymItens(iA))
    calcFxtFtSym(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
  else
    calcFxtFt(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fx, Fxt, &Ft);
#else
  calcFxtFtSym(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
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
#ifdef EDHE_FLEX
  int typei, typej;
#endif
  double phi, psi;
  double OmegaA[3][3];
  ti = x[4] + (trefG - atomTime[iA]);
  rA[0] = rx[iA] + vx[iA]*ti;
  rA[1] = ry[iA] + vy[iA]*ti;
  rA[2] = rz[iA] + vz[iA]*ti;
  vA[0] = vx[iA];
  vA[1] = vy[iA];
  vA[2] = vz[iA];
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  if (isSymItens(iA))
    UpdateOrient(iA, ti, RA, OmegaA);
  else
    symtop_evolve_orient(iA, ti, RA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(iA, ti, RA, OmegaA);
#endif
  MD_DEBUG2(printf("i=%d ti=%f", iA, ti));
  MD_DEBUG2(print_matrix(RA, 3));
#ifdef EDHE_FLEX
  typei = typeOfPart[iA];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1/Sqr(axa[iA]);
      invbSqN = 1/Sqr(axb[iA]);
      invcSqN = 1/Sqr(axc[iA]);
    }
  else
    {
      invaSqN = 1/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1/Sqr(axa[iA]);
  invbSqN = 1/Sqr(axb[iA]);
  invcSqN = 1/Sqr(axc[iA]);
#endif
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
  if (isSymItens(iA))
    calcFxtFtSym(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
  else
    calcFxtFt(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fx, Fxt, &Ft);
#else
  calcFxtFtSym(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
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
#ifdef EDHE_FLEX
  int typei;
#endif
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(i))
    {
      funcs2beZeroedNeighPlaneSE(n, x, fvec, i);
      return;
    }
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
#ifdef EDHE_FLEX
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1/Sqr(axa[i]);
      invbSqN = 1/Sqr(axb[i]);
      invcSqN = 1/Sqr(axc[i]);
    }
  else
    {
      invaSqN = 1/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
#endif
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
#ifdef EDHE_FLEX
  int typei;
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
#ifdef EDHE_FLEX
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1/Sqr(axa[i]);
      invbSqN = 1/Sqr(axb[i]);
      invcSqN = 1/Sqr(axc[i]);
    }
  else
    {
      invaSqN = 1/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
  tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, Rt);
#endif
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
#ifdef MD_SUPERELLIPSOID
extern void funcs2beZeroedDistNegNeighPlane5SE(int n, double x[], double fvec[], int i);
#endif
void funcs2beZeroedDistNegNeighPlane5(int n, double x[], double fvec[], int i)
{
  int k1, k2; 
  double fx[3], rD[3];
  /* x = (r, alpha, t) */ 
  
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(i))
    {
      funcs2beZeroedDistNegNeighPlane5SE(n, x, fvec, i);
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
extern void fdjacDistNegNeighPlane5SE(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int), int iA);

void fdjacDistNegNeighPlane5(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3], rD[3];
  int k1, k2;
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(iA))
    {
      fdjacDistNegNeighPlane5SE(n, x, fvec, df, vecfunc, iA);
      return;
    }
#endif
 
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
extern void fdjacDistNegNeighPlaneSE(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA);

void fdjacDistNegNeighPlane(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3];
  int k1, k2;
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(iA))
    {
      fdjacDistNegNeighPlaneSE(n, x, fvec, df, vecfunc, iA);
      return;
    }
#endif
 
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
extern void funcs2beZeroedDistNegNeighPlaneSE(int n, double x[], double fvec[], int i);

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
  
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(i))
    {
      funcs2beZeroedDistNegNeighPlaneSE(n, x, fvec, i);
      return;
    }
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
#ifdef MD_SUPERELLIPSOID
extern int calc_intersecSE(int i, double *rB, double *rA, double **Ri, double* rI);

extern double calcfLab(int i, double *x, double *rA, double **Ri);
extern long long numcoll;
void calc_intersec_neigh_planeSE(int i, double *rA, double *rB, double **Rt, double *grad, double* rC, double* rD)
{
  double rAA[3], rBA[3], rBAgrad;
  int k, k1;
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = rB[k1] - rA[k1];
  /* notare che grad è un vettore unitario */
  rBAgrad = scalProd(rBA, grad);
  //printf("rBAgrad=%f\n", rBAgrad);
  /* N.B. 06/03/09: il parallelepipedo usato per il calcolo dell'escape time è più piccolo 
     di quello reale di  OprogStatus.epsdNL+OprogStatus.epsd (vedi nnl.c:calc_grad_and_point_plane())
     per cui scegliendo come RD l'intersezione della normale al piano con il parallelepipedo
     reale (sommando la quantità testé detta) ci si assicura che rD sia esterno al SE */
  for (k1=0; k1 < 3; k1++)
    rD[k1] = rA[k1] + (rBAgrad+OprogStatus.epsdNL+OprogStatus.epsd)*grad[k1];

#if 0
  printf("A qui rA=%f %f %f f(rD)=%f\n", rA[0], rA[1], rA[2], calcfLab(i, rD, rA, Rt) );
  printf("rD=%f %f %f\n", rD[0], rD[1], rD[2]);
#endif
  /* N.B. 06/03/09: il punto rD deve essere sempre esterno al SE */
  calc_intersecSE(i, rD, rA, Rt, rC);
  /* il guess iniziale per rD va calcolato però con il piano attuale */
  for (k1=0; k1 < 3; k1++)
    rD[k1] = rA[k1] + rBAgrad*grad[k1];
#if 0
    {double dd[3], dd1[3], dd2[3], dd3[3]; int kk;
      for (kk=0; kk < 3; kk++)
	{
	  dd1[kk] = rC[kk]-rA[kk];
	  dd2[kk] = rD[kk]-rA[kk];
	  dd3[kk] = rD[kk]-rB[kk];
	  dd[kk]=rC[kk]-rD[kk];
	}
      printf("numcoll=%lld norm rCD=%.15G rCA.rDA:%.15G\n", numcoll, calc_norm(dd), scalProd(dd1,dd2));
      printf("rB.grad=%.15G\n", scalProd(dd3,grad));
      printf("rB=%f %f %f rA=%f %f %f\n", rB[0], rB[1], rB[2], rA[0], rA[1], rA[2]);
    }
#endif
  //printf("B qui\n");
  /* ...e ora calcoliamo rD (guess sul piano) */
#if 0
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = rB[k1] - rA[k1];
  /* notare che grad è un vettore unitario */
  rBAgrad = scalProd(rBA, grad);
  for (k1=0; k1 < 3; k1++)
    rD[k1] = rA[k1] + rBAgrad*grad[k1];
#endif
} 
#endif
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
      printf("NNL2 [calc_intersec] Serious problem guessing distance, aborting...\n");
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
  MD_DEBUG35(if (globalHW && Oparams.curStep > 5100) printf("[GuessSimple] rC=%f %f %f rD=%f %f %f\n", rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
}
#ifdef MD_SUPERELLIPSOID
extern void calc_intersec_neighSE(int i, double *rB, double *rA, double **Ri, double* rI);
#endif
void guess_distNeigh_plane(int i, 
		double *rA, double *rB, double **Xa, double *grad, double *rC, double *rD,
		double **RA)
{
  double gradA[3], gradaxA[3], dA[3], dB[3];
  int k1, n;
  double saA[3], sp;
#ifdef EDHE_FLEX
  int typei;
#endif
#ifdef MD_SUPERELLIPSOID
  double DdA[3], sfA, nDdA;
#endif

  //printf("===============>SONO QUI\n");
#ifdef EDHE_FLEX
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      saA[0] = axa[i];
      saA[1] = axb[i];
      saA[2] = axc[i];
    }
  else
    {
      saA[0] = typesArr[typei].sax[0];
      saA[1] = typesArr[typei].sax[1];
      saA[2] = typesArr[typei].sax[2];
    }
#else
  saA[0] = axa[i];
  saA[1] = axb[i];
  saA[2] = axc[i];
#endif
  for (k1 = 0; k1 < 3; k1++)
    gradA[k1] =  grad[k1];
  for (n = 0; n < 3; n++)
    {
      gradaxA[n] = 0;
      for (k1 = 0; k1 < 3; k1++) 
	gradaxA[n] += gradA[k1]*RA[n][k1];
    }
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(i))
    {
      sfA = 0.0;	
      for (k1=0; k1 < 3; k1++)
	{
	  sfA += Sqr(saA[k1]); 
	}
      sfA = sqrt(sfA)+OprogStatus.epsdNL+OprogStatus.epsd;
      for (k1=0; k1 < 3; k1++)
	{
	  dA[k1] = rA[k1];
	  DdA[k1] = 0;
	  for (n=0; n < 3;n++)
	    DdA[k1] += gradaxA[n]*RA[n][k1]*saA[n]/2.0; 
	}
      nDdA = calc_norm(DdA);
      for (k1=0; k1 < 3; k1++)
	{
	  dA[k1] = sfA*DdA[k1]/nDdA + dA[k1]; 
	}

      calc_intersec_neighSE(i, dA, rA, RA, rC);
    }
  else
    {
      for (k1=0; k1 < 3; k1++)
	{
	  dA[k1] = rA[k1];
	  for (n=0; n < 3;n++)
    	    dA[k1] += gradaxA[n]*RA[n][k1]*saA[n]/2.0; 
	}

      calc_intersec_neigh(dA, rA, Xa, rC, 1);
    }
#else
  for (k1=0; k1 < 3; k1++)
    {
      dA[k1] = rA[k1];
      for (n=0; n < 3;n++)
	dA[k1] += gradaxA[n]*RA[n][k1]*saA[n]/2.0; 
    }

  calc_intersec_neigh(dA, rA, Xa, rC, 1);
#endif
  for (k1=0; k1 < 3; k1++)
    dB[k1] = rB[k1] - rC[k1];
  sp = scalProd(dB, grad);
  for (k1=0; k1 < 3; k1++)
    rD[k1] = rC[k1] + sp*grad[k1];
  MD_DEBUG35(if (0 && globalHW && Oparams.curStep >= 5100) printf("[GuessOpt] rC=%f %f %f rD=%f %f %f\n", rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
  MD_DEBUG35(if (0 && globalHW && Oparams.curStep >= 5100) printf("[GuessOpt] grad=%f %f %f rA=%f %f %f rB=%f %f %f\n",
							     grad[0], grad[1], grad[2], rA[0], rA[1], rA[2], rB[0], rB[1], rB[2]));
}
#ifdef MD_SUPERELLIPSOID
extern int calc_intersecSE(int i, double *rB, double *rA, double **Ri, double* rI);
void calc_intersec_neighSE(int i, double *rB, double *rA, double **Ri, double* rI)
{
  calc_intersecSE(i, rB, rA, Ri, rI);
}
#endif
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
      printf("NNL [calc_intersec] Serious problem guessing distance, aborting...\n");
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
#ifdef EDHE_FLEX
  int typei;
#endif
  //printf("===============>SONO QUI\n");
#ifdef EDHE_FLEX
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      saA[0] = axa[i];
      saA[1] = axb[i];
      saA[2] = axc[i];
    }
  else
    {
      saA[0] = typesArr[typei].sax[0];
      saA[1] = typesArr[typei].sax[1];
      saA[2] = typesArr[typei].sax[2];
    }
#else
  saA[0] = axa[i];
  saA[1] = axb[i];
  saA[2] = axc[i];
#endif
  saB[0] = nebrTab[i].axa;
  saB[1] = nebrTab[i].axb;
  saB[2] = nebrTab[i].axc;
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
#ifdef MD_SUPERELLIPSOID
extern double calcf(double *x,int i);
extern void calcfxLabSE(int i, double *x, double *r, double **Ri, double fx[3]);
extern void distSD_NNL(int i, double *vecg, double lambda, int halfspring);
#endif
double calcDistNegNeighPlane(double t, double t1, int i, double *r1, double *r2, double *vecgsup, int calcguess, int calcgradandpoint, int *err, int nplane)

{
  /* NOTA: nplane = {0...7} e indica il piano rispetto al quale dobbiamo calcolare la distanza */
  double vecg[8], rC[3], rD[3], rDC[3], r12[3], invaSqN, invbSqN, invcSqN;
  double ti, segno;
  int retcheck;
  double nf, ng, gradf[3];
#ifdef MD_SUPERELLIPSOID
  double vecgcg[6], rCD[3];
  int tryagain=0;
#endif
#ifdef MD_ASYM_ITENS
  double psi, phi;
#else
  double Omega[3][3];
#endif
#ifdef EDHE_FLEX
  int typei;
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
#if defined(EDHE_FLEX) && defined (MD_EDHEFLEX_WALL)
  if (globalHW)
    {
      rB[0] = rA[0];
      rB[1] = rA[1];
    }
#endif
  MD_DEBUG20(printf("AAAA ti= %.15G rA (%.15G,%.15G,%.15G)\n", ti, rA[0], rA[1], rA[2]));
  MD_DEBUG20(printf("AAAA t1=%.15G atomTime[%d]=%.15G\n",t1,i,atomTime[i]));
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, RtA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, RtA, Omega);
#endif
#ifdef EDHE_FLEX
  //printf("INIZIO CALC DIST NEIGH PLABE t=%.15G\n", t);
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1.0/Sqr(axa[i]);
      invbSqN = 1.0/Sqr(axb[i]);
      invcSqN = 1.0/Sqr(axc[i]);
    }
  else
    {
      invaSqN = 1.0/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1.0/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1.0/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
#endif
#ifndef MD_SUPERELLIPSOID
  tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, RtA);
#else
  if (!is_superellipse(i))
    tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, RtA);
#endif
  //printf("ti= %.15G rNebrShell: %f\n", ti, OprogStatus.rNebrShell);
  ti = 0.0;
  MD_DEBUG20(printf("BBBB ti= %.15G rB (%.15G,%.15G,%.15G)\n", ti, rB[0], rB[1], rB[2]));
  /* NOTA: dato l'ellissoide e la sua neighbour list a t=0 bisognerebbe stimare con esattezza 
   * la loro distanza e restituirla di seguito */
  /* NOTA2: per ora fa in ogni caso il guess ma si potrebbe facilmente fare in modo
   * di usare il vecchio vecg. */
#if 0
  if (calcgradandpoint)
    calc_grad_and_point_plane(i, gradplane, rB, nplane);
#endif
retry:
  MD_DEBUG34(printf("PRI gradf=(%f,%f,%f) gradplane=(%f,%f,%f)\n",
		  gradf[0], gradf[1], gradf[2], gradplane[0], gradplane[1], gradplane[2]));
  if (OprogStatus.guessDistOpt==1)
    {
      guess_distNeigh_plane(i, rA, rB, Xa, gradplane, rC, rD, RtA);
    }
  else
    {
#ifdef MD_SUPERELLIPSOID
      if (is_superellipse(i))
	calc_intersec_neigh_planeSE(i, rA, rB, RtA, gradplane, rC, rD);
      else
	calc_intersec_neigh_plane(rA, rB, Xa, gradplane, rC, rD);
#else
      calc_intersec_neigh_plane(rA, rB, Xa, gradplane, rC, rD);
#endif
    }
#if 0
    {int uu;
      double aa[3], bb[3];
      for (uu=0; uu < 3; uu++)
	{aa[uu] = rC[uu]-rD[uu];
	  bb[uu]=rA[uu]-rB[uu];
	}
      printf("rCD=%.15G rAB=%.15G\n", calc_norm(aa), calc_norm(bb));
    }
#endif
  //printf("i=%d rA=%f %f %f rB=%f %f %f\n", i, rA[0], rA[1], rA[2], rB[0], rB[1], rB[2]);
  MD_DEBUG30(if (globalHW && Oparams.curStep >= 5100) printf("[GuessOpt] t=%.15G\n",t+t1));
#if 0
  /* check that point rC lies on the surface of object A (SE or HE) */
    {
      double xpA[3];
      lab2body(i, rC, xpA, rA, RtA);
      printf("after guess f(rC)=%.15G\n", calcf(xpA, i));
    }
#endif
  for(k1=0; k1 < 3; k1++)
    r12[k1] = rC[k1]-rD[k1]; 
  MD_DEBUG34(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		  rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
  //printf("nplane=%d rC=(%f,%f,%f) rD=(%f,%f,%f) r12=%.15G\n", nplane,
//		  rC[0], rC[1], rC[2], rD[0], rD[1], rD[2], calc_norm(r12));
#ifdef MD_SUPERELLIPSOID
#if 0
  if ((OprogStatus.SDmethod==1 || OprogStatus.SDmethod==4) || tryagain)
    {
      for (k1=0; k1 < 3; k1++)
	{
	  vecgcg[k1] = rC[k1];
	  vecgcg[k1+3] = rD[k1];
	}
      for (k1=0; k1 < 3; k1++)
	rCD[k1]=rC[k1]-rD[k1];
      //printf("PRIMA rC=%f %f %f rD=%f %f %f rCD=%.15G\n", rC[0], rC[1], rC[2], rD[0], rD[1], rD[2], calc_norm(rCD));
      distSD_NNL(i, vecgcg, OprogStatus.springkSD, 1);
      for (k1=0; k1 < 3; k1++)
	{
	  rC[k1] = vecgcg[k1];
	  rD[k1] = vecgcg[k1+3];
	}	
      for (k1=0; k1 < 3; k1++)
	rCD[k1]=rC[k1]-rD[k1];
      //printf("DOPO rC=%f %f %f rD=%f %f %fi rCD=%.15G\n", rC[0], rC[1], rC[2], rD[0], rD[1], rD[2], calc_norm(rCD));

    }
  //printf("DOPO SD nplane=%d rC=(%f,%f,%f) rD=(%f,%f,%f) r12=%.15G\n", nplane,
//		  rC[0], rC[1], rC[2], rD[0], rD[1], rD[2], calc_norm(rCD));

#endif
  if (is_superellipse(i))
    calcfxLabSE(i, rC, rA, RtA, gradf);
  else
    calc_grad(rC, rA, Xa, gradf);
#else
  calc_grad(rC, rA, Xa, gradf);
#endif
  MD_DEBUG34(printf("DOPO gradf=(%f,%f,%f) gradplane=(%f,%f,%f)\n",
		  gradf[0], gradf[1], gradf[2], gradplane[0], gradplane[1], gradplane[2]));
  nf = calc_norm(gradf);
  ng = calc_norm(gradplane);
  //printf("gradf.gradplane=%.15G\n", scalProd(gradf, gradplane)/(nf*ng));
#if 0
  printf ("gradf=%.15G %.15G %.15G gradg=%.15G %.15G %.15G\n", gradf[0], gradf[1], gradf[2],
	  gradplane[0], gradplane[1], gradplane[2]);
  printf("gradf.gradg=%.15G\n", scalProd(gradf, gradplane)/nf/ng);
#endif
  if (OprogStatus.dist5NL)
    vecg[3] = sqrt(nf/ng);
  else
    vecg[6] = sqrt(nf/ng);
  //printf("alpha guess is %.15G\n", vecg[6]);
  for (k1=0; k1 < 3; k1++)
    {
      vecg[k1] = rC[k1];
      if (!OprogStatus.dist5NL)
	vecg[k1+3] = rD[k1];
      rDC[k1] = rD[k1] - rC[k1];
    }

  //printf("boh[%d]=%.15G\n",i, scalProd(rDC,gradplane));
  if (OprogStatus.dist5NL)
    {
#ifdef MD_SUPERELLIPSOID
      if (is_superellipse(i))
	{
	  if (scalProd(rDC,gradplane) >= 0)
	    vecg[4] = calc_norm(rDC)/ng;
	  else
	    vecg[4] = -calc_norm(rDC)/ng;
	}
      else
	vecg[4] = 0.0;
#else
      vecg[4] = 0.0;
#endif
  /*vecg[4] = calc_norm(rDC);
   * QUESTO GUESS POTREBBE ESSERE MIGLIORE ANCHE SE IN PRATICA SEMBRA
   * LO STESSO. */
    }
  else
    {
#ifdef MD_SUPERELLIPSOID
      if (is_superellipse(i))
	{
	  if (scalProd(rDC,gradplane) >= 0)
	    vecg[7] = calc_norm(rDC)/ng;
	  else
	    vecg[7] = -calc_norm(rDC)/ng;
	  //printf("beta guess is %.15G\n", vecg[7]);
#ifdef MD_FDJAC_SYM
	  if (scalProd(rDC,gradplane) >= 0)
	    vecg[6] = calc_norm(rDC)/nf;
	  else
	    vecg[6] = -calc_norm(rDC)/nf;
#endif
	}
      else
	vecg[7] = 0.0;
#else
      vecg[7] = 0.0;
#endif
    }
#ifdef MD_SUPERELLIPSOID
#if 0
  if (!calcguess)
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
#endif
#endif
  //if (isnan(vecg[1]))
   //   printf("BEFORE vecg: %f %f %f %f %f %f\n", vecg[0], vecg[1], vecg[2], vecg[3], vecg[4], vecg[5]);
  MD_DEBUG(printf("alpha: %f beta: %f\n", vecg[6], vecg[7]));
  if (OprogStatus.dist5NL)
    newtDistNegNeighPlane(vecg, 5, &retcheck, funcs2beZeroedDistNegNeighPlane5, i); 
  else
    newtDistNegNeighPlane(vecg, 8, &retcheck, funcs2beZeroedDistNegNeighPlane, i); 

  //if (isnan(vecg[1]))
    //  printf("AFTER vecg: %f %f %f %f %f %f\n", vecg[0], vecg[1], vecg[2], vecg[3], vecg[4], vecg[5]);
#ifdef MD_SUPERELLIPSOID
  if (retcheck != 0)
    {
#if 0
      if (!tryagain && (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod==3))
	{
	  tryagain=1;
	  goto retry;
	}
#endif
      if (OprogStatus.targetPhi>0)
	{
	  calcdist_retcheck=1;
	  return 0.0;
	}  
      printf("[NNL] I couldn't calculate distance between %d and its NL, calcguess=%d, exiting....\n", i, calcguess);
      //printf("vec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", vecg[0], vecg[1], vecg[2], vecg[3], vecg[4],
      //     vecg[5], vecg[6], vecg[7]);
      exit(-1);
    }
#else
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
#endif
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
      //printf("[DISTNEIGHPLANE] t=%.15G distanza: %.15G\n", t, calc_norm(r12));
      return calc_norm(r12);
    }
  else
    {
      //printf("[DISTNEIGHPLANE] t=%.15G distanza: %.15G\n", t, -calc_norm(r12));
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
#ifdef EDHE_FLEX
  int typei;
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
#ifdef EDHE_FLEX
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1.0/Sqr(axa[i]);
      invbSqN = 1.0/Sqr(axb[i]);
      invcSqN = 1.0/Sqr(axc[i]);
    }
  else
    { 
      invaSqN = 1.0/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1.0/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1.0/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
#endif
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
#ifdef EDHE_FLEX
      if (OprogStatus.targetPhi > 0.0)
	return OprogStatus.rNebrShell*(min3(axa[i],axb[i],axc[i])/max3(axa[i],axb[i],axc[i]));
      else
  	return OprogStatus.rNebrShell*(min3(typesArr[typei].sax[0],typesArr[typei].sax[1],typesArr[typei].sax[2])
				     /max3(typesArr[typei].sax[0],typesArr[typei].sax[1],typesArr[typei].sax[2]));
#else
      return OprogStatus.rNebrShell*(min3(axa[i],axb[i],axc[i])/max3(axa[i],axb[i],axc[i]));
#endif
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
#ifdef EDHE_FLEX
extern int are_spheres(int i, int j);

double calcDistNegNNLoverlapPlaneHS(int i, int j, double rA[3], double rB[3])
{
  double RR, R01, DD[3];
  double EA[3], EB[3];
  int k1;
  EA[0] = nebrTab[i].axa;
  EA[1] = nebrTab[i].axb;
  EA[2] = nebrTab[i].axc;
  EB[0] = nebrTab[j].axa;
  EB[1] = nebrTab[j].axb;
  EB[2] = nebrTab[j].axc;
  for (k1 = 0; k1 < 3; k1++)
    {
      DD[k1] = rA[k1] - rB[k1];
    }
  /* axis C0+s*A0 */
  RR = fabs(DD[0]);
  R01 = EA[0] + EB[0];
  if ( RR > R01 )
    return 1.0; /* non si intersecano */
  /* axis C0+s*A1 */
  RR = fabs(DD[1]);
  R01 = EA[1] + EB[1];
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*A2 */
  RR = fabs(DD[2]);
  R01 = EA[2] + EB[2];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B0 */
  RR = fabs(DD[0]);
  R01 = EA[0] + EB[0];
  if ( RR > R01 )
    return 1.0;

  return -1.0;
}
#endif
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
#ifdef MD_EDHEFLEX_OPTNNL
  if (OprogStatus.optnnl)
    {
      rA[0] = rxNNL[i];
      rA[1] = ryNNL[i];
      rA[2] = rzNNL[i];
      rB[0] = rxNNL[j] + shift[0];
      rB[1] = ryNNL[j] + shift[1];
      rB[2] = rzNNL[j] + shift[2];
    }
  else
    {
      rA[0] = nebrTab[i].r[0];
      rA[1] = nebrTab[i].r[1];
      rA[2] = nebrTab[i].r[2];
      rB[0] = nebrTab[j].r[0] + shift[0];
      rB[1] = nebrTab[j].r[1] + shift[1];
      rB[2] = nebrTab[j].r[2] + shift[2];
    }
#else
  rA[0] = nebrTab[i].r[0];
  rA[1] = nebrTab[i].r[1];
  rA[2] = nebrTab[i].r[2];
  rB[0] = nebrTab[j].r[0] + shift[0];
  rB[1] = nebrTab[j].r[1] + shift[1];
  rB[2] = nebrTab[j].r[2] + shift[2];
#endif
#ifdef EDHE_FLEX
  if (is_a_sphere_NNL[i] && is_a_sphere_NNL[j])
    return calcDistNegNNLoverlapPlaneHS(i, j, rA, rB);
#endif

#if 0
  for (k1=0; k1<3; k1++)
    { 
      if (fabs(rA[k1]-rB[k1]) > L)
	{
	  printf("boh... shift=%f %f %f\n", shift[0], shift[1], shift[2]);
	  printf("rA=%f %f %f\n", rA[0], rA[1], rA[2]);
	  printf("rB=%f %f %f\n", nebrTab[j].r[0], rB[1], rB[2]);
	  exit(-1);
	}
    }
#endif
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
#ifdef EDHE_FLEX
  int typei, typej;
  double axaF, axbF, axcF;
#endif
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
	    {
#ifdef EDHE_FLEX
	      typej = typeOfPart[j];
	      if (OprogStatus.targetPhi > 0.0)
		{
		  axaF = axa[j];
		  axbF = axb[j];
		  axcF = axc[j];
		}
	      else
		{
		  axaF = typesArr[typej].sax[0];
		  axbF = typesArr[typej].sax[1];
		  axcF = typesArr[typej].sax[2];
		} 
	      g2 = OprogStatus.epsdGDO*min3(axaF,axbF,axcF)/calc_norm(vecnf); 
#else
  	      g2 = OprogStatus.epsdGDO*min3(axa[j],axb[j],axc[j])/calc_norm(vecnf); 
#endif
	    }
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

double brent_tref_plane, brentSign_plane, *brent_vecg_plane;
int iBrent_plane, brent_nplane;
double distbrent_plane(double t)
{
  double r1[3], r2[3];
  int distfail;
  return brentSign_plane*calcDistNegNeighPlane(t, brent_tref_plane, iBrent_plane, r1, r2, brent_vecg_plane, 
					       0, 0, &distfail, brent_nplane);
}

#define MD_BRENT_TOL 5E-16
extern int brentTooManyIter;
extern double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
int grazing_try_harder_plane(int i, double tref, double t1, double delt, double d1, double d2, double *vecg, double *troot, double *dmin, int nplane)
{
  int a;
#ifndef MD_GRAZING_TRYHARDER
  return 0;
#endif
  printf("[grazing_try_harder_plane] i=%d time=%.15G\n", i, tref+t1);
  /* Brent looks always for a minimum hence we have to change sign
     if grazing occurrs coming from negative distances */
  brentSign_plane = ((d1>0.0)?1.0:-1.0);
  brent_tref_plane = tref;
  brent_vecg_plane = vecg;
  iBrent_plane = i;
  brent_nplane= nplane;
  /* use brent to find the exact minimum */
  *dmin = brent(t1, t1+delt*0.5, t1+delt, distbrent_plane, MD_BRENT_TOL, troot);
  *dmin *= brentSign_plane;
  //printf("DOPO dmin(%.15G)=%.15G\n", *troot, *dmin);
  //printf("try harder cross dmin=%.15G\n", *dmin);
  if (!brentTooManyIter && *troot >= t1 && *troot <= t1+delt && *dmin*d1 < 0.0)
    {
      /* found a crossing! */
      /* 19/06/08 NOTE: note that we need times relative to tref, i.e.
	 we do not need to add tref to the solution *troot here!*/
      return 1;
    }
  return 0;/* no collision found */
}
int interpolNeighPlaneSNP(int i, double tref, double t, double delt, double d1, double d2, double *tmin, double* vecg,  int nplane)
{
  double d3, A;
  double r1[3], r2[3];
  int distfail;
  double dmin;
  d3 = calcDistNegNeighPlane(t+delt*0.5, tref, i, r1, r2, vecg, 0, 0, &distfail, nplane);
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
      //dmin = calcDistNeg(*tmin, tref, i, j, shift, r1, r2, &alpha, vecg, 0);
      dmin = calcDistNegNeighPlane(t+delt*0.5, tref, i, r1, r2, vecg, 0, 0, &distfail, nplane);
      *tmin += tref;
      return 0;
    }
  return 1;
}
#ifdef MD_GRAZING_TRYHARDER
extern double trefGT, *vecgGT;
extern int iGT;
int nplaneGT;
static double distfuncGT(double t)
{
  double r1[3], r2[3];
  int  distfail;
  return calcDistNegNeighPlane(t, trefGT, iGT, r1, r2, vecgGT, 
			       0, 0, &distfail, nplaneGT);
}
#endif

int interpolNeighPlane(int i, double tref, double t, double delt, double d1, double d2, double *troot, double* vecg, int bracketing, int nplane)
{
  int triedharder=0, a;
  int nb, distfail;
  double d3, t1, t2, A;
  double r1[3], r2[3], xb1[2], xb2[2];
  double tmin, dmin;
#ifdef MD_EDHEFLEX_WALL
  if (OprogStatus.paralNNL && !globalHW)
    assign_plane(nplane);
#else
  if (OprogStatus.paralNNL)
    assign_plane(nplane);
#endif
#ifdef MD_EDHEFLEX_WALL
  if (globalHW)
    sogliaErr_zbrent = OprogStatus.epsd;
  else
    sogliaErr_zbrent = OprogStatus.epsdNL;
#else
  sogliaErr_zbrent = OprogStatus.epsdNL;
#endif
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
      if (OprogStatus.zbrakn==0)
	{
	  if (ya[0]-ya[1] == 0.0)
	    {
	      tmin = t + delt*0.25;
	    }
	  else if (ya[2]-ya[0] == 0.0)
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
	      dmin = calcDistNegNeighPlane(tmin, tref, i, r1, r2, vecg, 0, 0, &distfail, nplane);
	      if (d1*dmin < 0.0)
		{
		  //printf("brcketing done d1=%.15G dmin=%.15G t1=%.15G t2=%.15G\n");
		  t2 = tmin;
		  t1 = t;
		}
	      else if (fabs(dmin) < fabs(d1) && fabs(dmin) < fabs(d2))
		{
		  /* differently from interpolSP we call grazing_try_harder_HE() func
		     from here because for sumnegpairs case there is a dedicated 
		     function called interpolSNP */
		  //printf("PRIMA delt=%.15G\n", delt);
		  //printf("PRIMA d1(%.15G)=%.15G dmin(%.15G) = %.15G d2(%.15G)=%.15G\n", t, d1, tmin, dmin, t+delt, d2);
		  /* call grazing_try_harder_plane(...) only if this routine is called from locateHardWall()
		   (i.e. globalHW != 0), se rimuovo la condizione su globalHW, grazing_try_harder_plane()
		   viene usato anche nella predizione dell'uscita dal bounding box delle NNL  */
#if defined(EDHE_FLEX) && defined(MD_EDHEFLEX_WALL)
		  //printf ("globalHW=%d dmin=%.15G d1=%.15G d2=%.15G\n", globalHW, dmin, d1, d2);
		  if (globalHW && grazing_try_harder_plane(i, tref, t, delt, d1, d2, vecg, &tmin, &dmin, nplane))
		    {
		      /* in grazing_try_harderHE() already checks for d1*dmin < 0.0 */
		      triedharder=1;
		      t2 = tmin;
	    	      t1 = t;
		    }
		  else
		    return 1;
#else
		  return 1;
#endif		
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
#ifdef MD_GRAZING_TRYHARDER
  /* se grazing_try_harderHE() ha trovato un minimo allora 
     non si deve usare in zbrent() il polinomio interpolante di secondo grado 
     (chiamando distfunc), che ha un minimo ma non negativo (notare che proprio
     sotto questa condizione si chiama grazing_try_harderHE()),
     ma la distanza esatta tra i e il piano nplane calcolata da distfuncGT.
     */
  if (triedharder)
   {
      trefGT = tref;
      vecgGT = vecg;
      iGT = i;
      nplaneGT=nplane;
      *troot=zbrent(distfuncGT, t1, t2, OprogStatus.zbrentTol);
    }
  else
    {
      *troot=zbrent(distfunc, t1, t2, OprogStatus.zbrentTol);
    }
#else
  *troot=zbrent(distfunc, t1, t2, OprogStatus.zbrentTol);
#endif
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
  /* N.B. 21/07/08: InterpolNeigh è la funzione per l'interpolazione per le NNL ellissoidali che non vengono però usate! */
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
  MD_DEBUG35(printf("^^^ sp = %.15G norm= %.15G\n", fabs(scalProd(ddot, gradplane)), calc_norm(ddot)));
  MD_DEBUG35(printf("w=%.15G %.15G %.15G rcat=%.15G %.15G %.15G v=%f %f %f wra=%f %f %f\n", wx[i], wy[i], wz[i],
		    rcat[0], rcat[1], rcat[2], vx[i], vy[i], vz[i], wra[0],wra[1], wra[2]));
#if 0
  return calc_norm(ddot);
#else
  /* WARNING:  quella che conta è la velocità rispetto alla perpendicolare al piano, 
   * poiché la distanza tra piano ed ellisse è sempre perpendicolare al piano...tuttavia
   * verificare e testare tale maggiorazione !!! */
  return fabs(scalProd(ddot, gradplane));
#endif
}
void adjust_maxddot(int i, double *maxddot);

int search_contact_faster_neigh_plane(int i, double *t, double t1, double t2, 
				double *vecgd, double epsd, double *d1, double epsdFast,
				double *r1, double *r2, int nplane)
{
  /* NOTA: 
   * MAXOPTITS è il numero massimo di iterazioni al di sopra del quale esce */
  double maxddot, told, delt, normddot, ddot[3];
  const int MAXOPTITS = 500;
  double factori, A, B;
  int its=0, distfailed, itsf=0; 
  const double GOLD= 1.618034;  
  factori = 0.5*maxax[i]+epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));

  /* estimate of maximum rate of change for d */
#if 0
  maxddot = sqrt(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
#else
  /* WARNING: questa maggiorazione si applica al caso specifico dell'urto di un ellissoide con un piano,
   * è molto migliore della precedente ma va testata accuratamente! */
#ifdef MD_ASYM_ITENS
  maxddot = calc_maxddot_nnl(i, gradplane, epsd);
#else
  A = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2]);
  B = sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
  if (B/A < MAXDDOT_NNL_THR)
    maxddot = A*ADJ_MAXDDOT_NL + B;
  else
    maxddot = A + B;
 // maxddot = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
   // sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
#endif
  adjust_maxddot(i, &maxddot);
  //printf("SCALPROD = %.15G\n", vx[0]*gradplane[0]+vy[1]*gradplane[1]+vz[2]*gradplane[2]);
  *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed, nplane);
  timesFNL++;
  MD_DEBUG35(printf("Pri distances between %d d1=%.12G epsd*epsdTimes:%f\n", i, *d1, epsdFast));
  //printf("[SEARCH_CONTACT_FASTER] t=%.15G ellips N. %d d=%.15G\n", *t, i, *d1); 
  while (*d1 > epsdFast && its < MAXOPTITS)
    {
      told = *t;
      delt = *d1 / maxddot;
      //printf("normddot: %.15G\n", epsd/normddot);
      /* check for convergence */
#if defined(EDHE_FLEX) || defined(MD_BASIC_DT)
      if (delt < (epsd / maxddot))
	{
	  MD_DEBUG35(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
      *t += delt;
#else 
      normddot = calcvecFNeigh(i, *t, t1, ddot, r1);
      if (normddot!=0 && delt < (epsd / normddot))
	{
	  MD_DEBUG35(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
      *t += delt;
#endif
#if 1
      if (*t + t1 > t2)
	{
	  *t = told;
	  MD_DEBUG35(printf("t>t2 %d iterations reached t=%f t2=%f\n", its, *t, t2));
	  MD_DEBUG35(printf("convergence t>t2\n"));
	  *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 0, &distfailed, nplane);
	  return 1;
	}
#endif
      *d1 = calcDistNegNeighPlane(*t, t1, i, r1, r2, vecgd, 1, 1, &distfailed, nplane);
      MD_DEBUG35(printf("Inside Loop distances between %d d1=%.12G\n", i, *d1));
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
      while (*d1 <= 0 || distfailed)
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
      //printf("nplane=%d no contact point found\n", nplane);
      MD_DEBUG35(printf("newtNeigh did not find any contact point!\n"));
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
  //if (veln == 0.0)
    //return 0;
  if (veln <= 0 && dd > maxax[i])
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
      if ((fabs(dists[nn]) < 1E-14 || dists[nn] < 0.0)  && distsOld[nn] > 0.0)
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
      /* N.B. 21/07/08: questa routine viene usata in locate_contact_neigh_plane_parall() quindi non 
	 per il calcolo dell'urto con una muro, per cui si puo' usare OprogStatus.epsdNK */
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
      //printf("rB=%f %f %f vecgd: %f %f %f %f %f %f\n", rB[0], rB[1], rB[2], vecgd[nn][0], vecgd[nn][1], vecgd[nn][2], vecgd[nn][3], vecgd[nn][4], vecgd[nn][5]);
#if 0
      if (i==121)
	printf("nn=%d calcguess=%d\n", nn, calcguess);
#endif
      dists[nn] = calcDistNegNeighPlane(t, tref, i, r1, r2, vecgd[nn], calcguess, 0, &err, nn);
      //printf("NNL i=%d dist[%d]:%.15G\n", i, nn, dists[nn]);
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
      //if (dt < 1E-14)
	//printf("nn=%d dt=%.15G delt=%.15G dists=%.15G maxddoti=%15G\n", nn, dt, *delt, dists[nn], maxddoti[nn]);
      if (nn==0 || dt < (*delt))
	*delt = dt;
    }
  //printf("I chose dt=%.15G\n", *delt);
}
const double mddotfact = 1.001;
void adjust_maxddot(int i, double *maxddot)
{
  double K = 1.0;
#ifdef EDHE_FLEX
  int typei;
  double axaL, axbL, axcL;
#endif
#ifdef MD_ASYM_ITENS
  if (Mx[i] == 0.0 && My[i] == 0.0 && Mz[i] == 0.0)
    K = mddotfact;
#else
  if (wx[i] == 0.0 && wy[i] == 0.0 && wz[i] == 0.0)
    K = mddotfact;
#endif
#ifdef MD_POLYDISP
  if (axaP[i] == axbP[i] && axbP[i] == axcP[i])
    K = mddotfact;
#else
#ifdef EDHE_FLEX
  typei = typeOfPart[i];
  if (OprogStatus.targetPhi > 0.0)
    {
      axaL = axa[i];
      axbL = axb[i];
      axcL = axc[i];
    }
  else
    {
      axaL = typesArr[typei].sax[0];
      axbL = typesArr[typei].sax[1];
      axcL = typesArr[typei].sax[2];
    }
  if (axaL == axbL && axbL == axcL)
    K = mddotfact;
#else
  if (i < Oparams.parnumA)
    {
      if (Oparams.a[0] == Oparams.b[0] && Oparams.b[0] == Oparams.c[0])
	K = mddotfact;
    }
  else
    {
      if (Oparams.a[1] == Oparams.b[1] && Oparams.b[1] == Oparams.c[1])
	K = mddotfact;
    }
#endif
#endif
  *maxddot *= K;
}
void adjust_maxddoti(int i, double *maxddot, double maxddotiLC[6], double maxddoti[6])
{
  double K = 1.0;
  int a;
#ifdef EDHE_FLEX
  int typei;
  double axaL, axbL, axcL;
#endif
#ifdef MD_ASYM_ITENS
  if (Mx[i] == 0.0 && My[i] == 0.0 && Mz[i] == 0.0)
    K = mddotfact;
#else
  if (wx[i] == 0.0 && wy[i] == 0.0 && wz[i] == 0.0)
    K = mddotfact;
#endif
#ifdef MD_POLYDISP
  if (axaP[i] == axbP[i] && axbP[i] == axcP[i])
    K = mddotfact;
#else
#ifdef EDHE_FLEX
  typei = typeOfPart[i];
  if (OprogStatus.targetPhi > 0.0)
    {
      axaL = axa[i];
      axbL = axb[i];
      axcL = axc[i];
    }
  else
    {
      axaL = typesArr[typei].sax[0];
      axbL = typesArr[typei].sax[1];
      axcL = typesArr[typei].sax[2];
    }
  if (axaL == axbL && axbL == axcL)
    K = mddotfact;
#else
  if (i < Oparams.parnumA)
    {
      if (Oparams.a[0] == Oparams.b[0] && Oparams.b[0] == Oparams.c[0])
	K = mddotfact;
    }
  else
    {
      if (Oparams.a[1] == Oparams.b[1] && Oparams.b[1] == Oparams.c[1])
	K = mddotfact;
    }
#endif
#endif
  *maxddot *= K;
  for (a = 0; a < 6; a++)
    maxddoti[a] = K*maxddotiLC[a];
}
int search_contact_faster_neigh_plane_all(int i, double *t, double t1, double t2, 
					  double vecgd[6][8], double epsd, 
					  double *d1, double epsdFast, 
					  double dists[6], double maxddotiLC[6], double maxddot)
{
  double told, delt=1E-15, distsOld[6];
  const double GOLD= 1.618034;
  const int MAXOPTITS = 500;
  double maxddoti[6];
  int its=0, crossed[6], itsf; 

  *d1 = calcDistNegNeighPlaneAll(*t, t1, i, distsOld, vecgd, 1);
  //printf("SCF NNL *d1=%.15G t=%.15G\n", *d1, *t);
  MD_DEBUG36(printf("[SEARCH_CONTACT_FASTER_NNL_PARALL]t=%.15G d=%.15G\n", *t, *d1));
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
  adjust_maxddoti(i, &maxddot, maxddotiLC, maxddoti);

  while (fabs(*d1) > epsdFast && its < MAXOPTITS)
    {
#if 1
      calc_delt(maxddoti, &delt, distsOld);
      /* 10/06/2010 a check_distance prima passavo erroneamente dists e non distsOld 
	 ma dists poteva non essere inizializzato!! grazie valgrind! :-)*/
      if (check_distance(maxddoti, distsOld, t1, t2, *t))
	{
	  return 1;
	}
#else
      delt = fabs(*d1) / maxddot;
      if (*t + t1 < t2 && (t2 - (*t + t1))*maxddot < fabs(*d1) - OprogStatus.epsd)
	return 1;
#endif
      *t += delt;

      //printf("PRIMA i=%d SCF dist=%.15G\n",i,*d1);
      *d1 = calcDistNegNeighPlaneAll(*t, t1, i, dists, vecgd, 1);
#if 0
      if (its > 100 && its%10 == 0)
	{
	  printf("NNL SEARCH CONTACT FASTER t=%.15G its=%d\n", *t+t1, its);
	}
#endif
      //if (delt < 1E-12)
	//printf("maxddot=%.15G delt=%.15G d=%.15G t=%.15G\n", maxddot, delt, *d1, *t+t1);

      //printf("DOPO i=%d SCF dist=%.15G\n",i,*d1);
#if 1
      itsf = 0;
      while (check_cross(distsOld, dists, crossed)||(*d1==0.0))
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
     if (check_cross(distsOld, dists, crossed)||(*d1==0.0))
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
#ifdef EDHE_FLEX
int locate_contact_neigh_plane_HS_one(int i, double *evtime, double t2)
{
  int nn, typei;
  double t1, b, dist, dr[3], colltime=0.0, dv[3];
  typei = typeOfPart[i];

  t1 = Oparams.time;
  dr[0] = rx[i] - rB[0];
  dr[1] = ry[i] - rB[1];
  dr[2] = rz[i] - rB[2];  
#ifdef MD_ABSORPTION
  if (globalHW==2 && typeOfPart[i]==2)
    {
      //printf("qui globalHW=%d type=%d\n", globalHW, typeOfPart[i]);
      if (dr[2]+typesArr[typeOfPart[i]].sax[0] < 0.0)
	{
	  MD_DEBUG38(printf("time=%.15G i=%d switched to type 1 rz[]=%f\n",Oparams.time, i, rz[i])); 
	  typeOfPart[i]=1;
	  //return 0;
	}
      else
	{
	  return 0;
	}
    }
#endif

  dv[0] = vx[i];
  dv[1] = vy[i];
  dv[2] = vz[i];
  //printf("calc_norm(dr)=%.15G v=%.15G %.15G %.15G grad=%.15G %.15G %.15G\n", calc_norm(dr), vx[i], vy[i], vz[i],
  //     gradplane_all[nn][0], gradplane_all[nn][1], gradplane_all[nn][2]);
  /* N.B. controllare che il gradiente sia a norma unitaria e che sia uscente rispetto 
     al parallelepipedo delle NNL! */
  if (OprogStatus.targetPhi > 0.0)
    {
      dist = fabs(scalProd(dr, gradplane)) - axa[i];
    }
  else
    {
      dist = fabs(scalProd(dr, gradplane)) - typesArr[typei].sax[0];
    }
  b = scalProd(dv, gradplane);
  if (b < 0)
    return 0;
#if 0
  if (typeOfPart[i]==1 && (rz[i]-L*0.5-OprogStatus.bufHeight+typesArr[typei].sax[0] > 0.0 || 
			   rz[i]-0.5 < -L*0.5))
    {

      printf("i=%d dist=%.15G b=%.15G typei=%d incellz=%d rz=%.15G\n", i, dist, b, typeOfPart[i], inCell[2][i], rz[i]);
      printf("gradplane=%f %f %f point=%f %f %f dr=%f %f %f ghw=%d\n", gradplane[0], gradplane[1], gradplane[2], rB[0],
	     rB[1], rB[2], dr[0], dr[1], dr[2], globalHW);\
      exit(-1);
    }
#endif  
  colltime = dist/b+Oparams.time;
  if (colltime > t1 && colltime < t2)
    {
      *evtime = colltime;	
      return 1;
    }
  MD_DEBUG36(printf("t1=%.15G t2=%.15G colltime=%.15G\n", t1, t2, *evtime));
  return 0;
}

int locate_contact_neigh_plane_HS(int i, double *evtime, double t2)
{
  int nn, typei, first=1;
  double t1, b, dist, dr[3], ti, t, colltime=0.0, dv[3];
  typei = typeOfPart[i];
  t1 = Oparams.time;

  for (nn = 0; nn < 6; nn++)
    {
      dr[0] = rx[i] - rBall[nn][0];
      dr[1] = ry[i] - rBall[nn][1];
      dr[2] = rz[i] - rBall[nn][2];  
      dv[0] = vx[i];
      dv[1] = vy[i];
      dv[2] = vz[i];
      //printf("calc_norm(dr)=%.15G v=%.15G %.15G %.15G grad=%.15G %.15G %.15G\n", calc_norm(dr), vx[i], vy[i], vz[i],
	//     gradplane_all[nn][0], gradplane_all[nn][1], gradplane_all[nn][2]);
      /* N.B. controllare che il gradiente sia a norma unitaria e che sia uscente rispetto 
	 al parallelepipedo delle NNL! */
      if (OprogStatus.targetPhi > 0.0)
	{
	  dist = fabs(scalProd(dr, gradplane_all[nn])) - axa[i];
	}
      else
	{
	  dist = fabs(scalProd(dr, gradplane_all[nn])) - typesArr[typei].sax[0];
	}
      b = scalProd(dv, gradplane_all[nn]);
      if (b < 0)
	continue;
      //printf("dist=%.15G b=%.15G\n", dist, b);
      colltime = dist/b+Oparams.time;
      if (colltime > t1 && colltime < t2)
	{
	  if (colltime < *evtime || first)
	    {
	      first = 0;  
	      *evtime = colltime;	
	    }
	}
    }	
  MD_DEBUG36(printf("t1=%.15G t2=%.15G colltime=%.15G\n", t1, t2, *evtime));
  if (first)
    return 0;
  else
    return 1;
}
#endif
int locate_contact_neigh_plane_parall(int i, double *evtime, double t2)
{
  /* const double minh = 1E-14;*/
  double h, d, dold, t2arr[6], t, dists[6], distsOld[6], 
	 vecg[5], vecgroot[6][8], vecgd[6][8], vecgdold[6][8], factori, A, B; 
#ifndef MD_BASIC_DT
  double normddot, distsOld2[6], vecgdold2[6][8], dold2, deldist;
#endif
  double maxddot, delt, troot, tini, maxddoti[6];
  int firstev;
#ifdef MD_PARANOID_CHECKS
  double deltini;
#endif
  /*
  const int MAXITS = 100;
  const double EPS=3E-8;*/ 
  /* per calcolare derivate numeriche questo è il magic number in doppia precisione (vedi Num. Rec.)*/
  const double GOLD= 1.618034;
  int its, foundrc;
  double t1, epsd, epsdFast, epsdFastR, epsdMax; 
  int kk,tocheck[6], dorefine[6], ntc, ncr, nn, gotcoll, crossed[6], firstaftsf;
#ifdef EDHE_FLEX
  if (typesArr[typeOfPart[i]].ignoreCore)
    return 0;
#endif
  epsd = OprogStatus.epsdNL;
  epsdFast = OprogStatus.epsdFastNL;
  epsdFastR= OprogStatus.epsdFastRNL;
  epsdMax = OprogStatus.epsdMaxNL;
  t = 0;//t1;
  t1 = Oparams.time;
  //t2 = timbig;
  //printf("inizio locate contact \n");
#ifdef MD_HANDLE_INFMASS
  if (is_infinite_Itens(i) && is_infinite_mass(i))
    {
      *evtime = timbig;
      return 1;
    }
  if (is_infinite_mass(i) && is_sphere(i))
    {
      *evtime = timbig;
      return 1;    
    }
#endif
  calc_grad_and_point_plane_all(i, gradplane_all, rBall);
  /* la collisione di una sfera con i vari piani si puo' calcolare 
     velocemente senza passare per il codice che segue */
#ifdef EDHE_FLEX
  if (is_sphere(i))// nel caso di superellissoidi deve anche valere che gli esponenti
    // siano tutti (fix this bug but do some testing before)
    {
#ifdef MD_SUPERELLIPSOID
	if (!is_superellipse(i))
#endif
	  //printf("qui\n");
	  return locate_contact_neigh_plane_HS(i, evtime, t2);
      MD_DEBUG37(printf("HS evtime=%.15G\n", *evtime));
    }
#endif
  factori = 0.5*maxax[i]+epsd;
  maxddot = 0.0;
  for (nn = 0; nn < 6; nn++)
    {
#ifdef MD_ASYM_ITENS
      maxddoti[nn] = calc_maxddot_nnl(i, gradplane_all[nn], epsd);
#else
      A = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2]);
      B = sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
      /* 19/05/2010: se vale questa condizione A + B = A quindi bisogna correggere A */
      if (B / A < MAXDDOT_NNL_THR)
	maxddoti[nn] = A*ADJ_MAXDDOT_NL + B;
      else
	maxddoti[nn] = A + B;
      //maxddoti[nn] = fabs(vx[i]*gradplane_all[nn][0]+vy[i]*gradplane_all[nn][1]+vz[i]*gradplane_all[nn][2])+
	//sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
      if (nn==0 || maxddoti[nn] > maxddot)
	maxddot = maxddoti[nn];
      //printf("nn=%d maxddoti=%.15G\n", nn, maxddoti);
    }
  h = OprogStatus.h; /* last resort time increment */
  delt = h;
 

  MD_DEBUG36(printf("[LOCATE_CONTACT_PARALL_NNL] BEGIN t=%.15G\n", t)); 
  //printf("STARTING i=%d\n",i);
  if (search_contact_faster_neigh_plane_all(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, 
					       dists, maxddoti, maxddot))
    {
      return 0;  
    }
  MD_DEBUG36(printf("[LOCATE_CONTACT_PARALL_NNL] dopo primo search_contact_faster\n")); 
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
      //printf("prima i=%d d=%.15G\n", i, d);
      while ((d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgd, 0))==0.0)
	{

	  //printf("d=%.15G QUIIIIIIIIIIIQUOOO\n",d);
	  delt *= GOLD;
	  t = tini + delt;
	}
      //printf("NNLlocate t=%.15G d=%.15G QUIIIIIIIIIIIQUOOO\n",t, d);
#ifdef MD_PARANOID_CHECKS
      if ((fabs(d-dold)==0.0))
	{
	  /* first we reduce delt to check d=dold is just random */
	  deltini = delt;
	  delt /= GOLD;
	  t = tini + delt;
	  d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgd, 0);
	  if (fabs(d-dold)!=0.0)
	    {
	      delt = deltini;
	      t = tini + delt;
	      d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgd, 0);
	    } 
	  while (fabs(d-dold)==0.0)
	    {
	      delt *= GOLD;
	      t = tini + delt;
	      d = calcDistNegNeighPlaneAll(t, t1, i, dists, vecgd, 0);
	    }	
	}
#endif
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
      MD_DEBUG36(printf("[LOCATE_CONTACT_PARALL_NNL] dentro loop t=%.15G d=%.15G\n", t, d)); 
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
	      //printf("TO CHECK\n");
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
		  MD_DEBUG36(printf("[LOCATE_CONTACT_PARALL_NNL] do refine t1=%.15G t=%.15G vecg[4]:%.15G\n", t, t2arr[nn], vecg[4])); 
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
			{
			  MD_DEBUG37(printf(">>>1 evtime=%.15G\n", *evtime));
			  return 1;
			}
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
	{
	  //printf("i=%d BECCATA\n",i);
	  MD_DEBUG37(printf(">>>2 evtime=%.15G\n", *evtime));
	  return 1;
	}
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
#ifdef MD_PATCHY_HE
extern struct LastBumpS *lastbump;
#endif
int locate_contact_neigh_plane(int i, double vecg[5], int nplane, double tsup)
{
  double h, d, dold, vecgd[8], vecgdold[8], t, r1[3], r2[3]; 
  double dtmp, t1, t2, maxddot, delt, troot, vecgroot[8], A, B;
  //const int MAXOPTITS = 4;
  const double GOLD= 1.618034;
#ifndef MD_BASIC_DT
  double ddot[3], dold2, vecgdold2[8], normddot;
#endif
#ifndef MD_ASYM_ITENS
  double factori;
#endif
  double epsd, epsdFast, epsdFastR, epsdMax; 
  int dorefine, distfail;
  int its, foundrc, kk;
  double tini, tmin;
  int sumnegpairs=0;
#ifdef EDHE_FLEX
  if (typesArr[typeOfPart[i]].ignoreCore==2)
    return 0;
#endif
#if defined(EDHE_FLEX) && defined(MD_EDHEFLEX_WALL)
  if (globalHW)
    {
      /* 09/07/08: if we are looking for a collision with a wall
	 be more accurate than in the case of NNL, 
	 grazing collision can be a problem in this case,
	 se in futuro volessi usare grazing_try_harder_plane() anche per le NNL
	 allora questo if va rimosso e epsdFastNL, epsdFastRNL, ecc. vanno tutti 
	 settati nel file .par con valori intorno a 1E-4/1E-5 
       */
      epsd = OprogStatus.epsdPlane;
      epsdFast = OprogStatus.epsdFastPlane;
      epsdFastR= OprogStatus.epsdFastPlane;
      epsdMax = OprogStatus.epsdPlane;
    }
  else
    {
      epsd = OprogStatus.epsdNL;
      epsdFast = OprogStatus.epsdFastNL;
      epsdFastR= OprogStatus.epsdFastRNL;
      epsdMax = OprogStatus.epsdMaxNL;
    }
#else
  epsd = OprogStatus.epsdNL;
  epsdFast = OprogStatus.epsdFastNL;
  epsdFastR= OprogStatus.epsdFastRNL;
  epsdMax = OprogStatus.epsdMaxNL;
#endif
  /* NOTA: implementare le varie funzioni _neigh (search_contact_faster_neigh, ecc.)
   * in tali funzioni la particella j non è altro che un ellissoide più grande di i
   * con lo stesso centro e immobile */
  t = 0.0;//Oparams.time;
  timesSNL++;
  calc_grad_and_point_plane(i, gradplane, rB, nplane);
  dtmp = calcDistNegNeighPlane(t, Oparams.time, i, r1, r2, vecgd, 0, 0, &distfail, nplane);
#ifdef MD_ABSORPTION
  if (globalHW==2 && typeOfPart[i]==2)
    {
      if (dtmp > 0.0)
	{
	  typeOfPart[i]=1;
	  return 0;
	}
      else
	{
	  return 0;
	}
    }
#endif
  if (dtmp < 0)
    {
#ifdef MD_EDHEFLEX_WALL
      /* N.B. nel caso di urto con muro la distanza può essere negativa all'inizio! */
      if (!globalHW)
	{
	  printf("La distanza fra l'ellissoide N. %d e il piano %d e' negativa d=%.15G\n", i, nplane, dtmp);
	  printf("nexttime[%d]: %.15G\n", i, nebrTab[i].nexttime);
	  exit(-1);
	}
      else if (!(lastbump[i].mol==nplane && lastbump[i].type==MD_WALL))
	{
	  printf("[WARNING] distance (d=%.15G) between particle N. %d and plane %d along z is negative!\n", dtmp, i, nplane);
	}
#else
      printf("La distanza fra l'ellissoide N. %d e il piano %d e' negativa d=%.15G\n", i, nplane, dtmp);
      printf("nexttime[%d]: %.15G\n", i, nebrTab[i].nexttime);
      exit(-1);
#endif
    }
  if (!bracket_neigh(i, &t1, &t2, nplane))
    {
      //printf("NOT BRACK NEIGH\n");
      return 0;
    }
  t1 += Oparams.time;	
  if (t1 > tsup)
    return 0;
  if (tsup < Oparams.time+t2)
    t2 = tsup;
  else
    t2 += Oparams.time;
  MD_DEBUG35(if (globalHW) printf("[GuessOpt] t1=%.15Gi t2=%.15G\n",t1, t2));
  //printf("t1=%.15G t2=%.15G time=%.15G tsup=%.15G\n", t1, t2, Oparams.time, tsup);
  //printf("LOCATE_CONTACT_NNL nplane=%d grad=%.8f %.8f %.8f  rB=%.8f %.8f %.8f t1=%.8f t2=%.8f tsup=%.8f maxax[%d]=%f\n", nplane, 
  //	 gradplane[0], gradplane[1], gradplane[2], rB[0], rB[1], rB[2], t1, t2, tsup, i, maxax[i]);
#ifdef MD_ASYM_ITENS
  maxddot = calc_maxddot_nnl(i, gradplane, epsd);
#else
  factori = 0.5*maxax[i]+epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
  A = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2]);
  B = sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;
  if (B / A < MAXDDOT_NNL_THR)
    maxddot =  A*ADJ_MAXDDOT_NL + B;
  else
    maxddot =  A + B;
  //maxddot = fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2])+
    //sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*factori;  
#endif
  h = OprogStatus.h; /* last resort time increment */
#ifdef MD_EDHEFLEX_WALL
  if (globalHW && lastbump[i].mol==nplane && lastbump[i].type==MD_WALL)
    {
      sumnegpairs=1;
    }
#endif
  if (search_contact_faster_neigh_plane(i, &t, t1, t2, vecgd, epsd, &d, epsdFast, r1, r2, nplane))
    return 0;  
  MD_DEBUG35(printf(">>>>d:%f\n", d));
  foundrc = 0;
  for (kk = 0; kk < 8; kk++)
    vecgdold[kk] = vecgd[kk];
  dold = d;
  its = 0;
  /* dist_too_big seems to be an unsafe optimization */
  while (t + t1 < t2)// || !dist_too_big(i,t,t1))
    {
      MD_DEBUG31(printf("LOC CONT rB = %.15G %.15G %.15G\n", rB[0], rB[1], rB[2]));
#if defined(MD_BASIC_DT) 
      delt = epsd / maxddot;
      tini = t;
      t += delt;
#ifdef MD_EDHEFLEX_WALL
      if (globalHW)
	{
	  while ((d = calcDistNegNeighPlane(t, t1, i, r1, r2, vecgd, 0, 0, &distfail, nplane))==0.0)
	    {
	      delt *= GOLD;
	      t = tini + delt;
	    }
	}
      else
#endif
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
      MD_DEBUG35(printf("NNL [LOCATE_CONTACT] delt=%.15G (%.15G,maxddot=%.10G,normddot=%.15G) t=%.15G ellips N. %d d=%.15G dold=%.15G its=%d t1=%.15G t2=%.15G vparall=%.15G\n", 
      delt, epsd/maxddot, maxddot, normddot, t, i, d, dold, its, t1, t2, fabs(vx[i]*gradplane[0]+vy[i]*gradplane[1]+vz[i]*gradplane[2]))); 
      if (fabs(d-dold2) > epsdMax)
	{
	  /* se la variazione di d è eccessiva 
	   * cerca di correggere il passo per ottenere un valore
	   * più vicino a epsd*/
	  MD_DEBUG35(printf("P delt: %.15G d-d2o:%.15G d2o:%.15G\n", delt, fabs(d-dold2), dold2));
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
	  MD_DEBUG35(printf("new dist d=%.15G delt=%.15G\n", d, delt));
	  for (kk = 0; kk < 8; kk++)
	    vecgd[kk] = vecgdold2[kk];
	  //printf("D delt: %.15G d2-d2o:%.15G d2:%.15G d2o:%.15G\n", delt*epsd/fabs(d2-d2old), fabs(d2-d2old), d2, d2old);
	}
#endif
#ifdef MD_EDHEFLEX_WALL
      if (globalHW && sumnegpairs)
	{
	  MD_DEBUG(printf("sumnnegpairs d=%.15G\n", d));
	  if (d <= 0.0)
	    {
	      if (!interpolNeighPlaneSNP(i, t1, t-delt, delt, dold, d, &tmin, vecgd, nplane))	
		  {
		    tmin -= t1;
		    delt = tmin - tini;
		    t = tmin;
		    d = calcDistNegNeighPlane(t, t1, i, r1, r2, vecgd, 0, 0, &distfail, nplane);
		  }
	    }
	  sumnegpairs = 0;
	}
#endif
      MD_DEBUG(printf(">>>> t = %f d1:%f d2:%f d1-d2:%.15G\n", t, d1, d2, fabs(d1-d2)));
      dorefine = 0;      
      if (dold > 0 && d < 0)
	{
	  MD_DEBUG35(printf(">>>>>>>>>QUI\n"));
       	  for (kk=0; kk < 8; kk++)
	    vecgroot[kk] = vecgd[kk];
#ifndef MD_NOINTERPOL  
	 if (interpolNeighPlane(i, t1, t-delt, delt, dold, d, &troot, vecgroot, 0, nplane))
#endif
	    {
	      /* vecgd2 è vecgd al tempo t-delt */
	      for (kk=0; kk < 8; kk++)
		vecgroot[kk] = vecgdold[kk];
	      /* VECCHIA SOLUZIONE: troot = t + t1 - delt;*/
	      troot = t + t1;
	    }
	  dorefine = 1;
	}
      else if (dold < epsd && d < epsd)
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
	  MD_DEBUG35(printf("NNL REFINING CONTACT t=%.15G troot=%.15G t1=%.15G\n", t, troot, t1));
	  if (refine_contact_neigh_plane(i, t1, troot, vecgroot, vecg, nplane))
	    {
	      MD_DEBUG35(printf("[locate_contact] Adding collision between %d\n", i));
	      MD_DEBUG35(printf("collision will occur at time %.15G\n", vecg[4])); 
	      MD_DEBUG35(printf("[locate_contact] its: %d\n", its));
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
	      MD_DEBUG35(printf("[locate_contact] can't find contact point!\n"));
	      if (d < 0)
		{
		  MD_DEBUG35(printf("t=%.15G d2 < 0 and I did not find contact point, boh...\n",t));
		  MD_DEBUG35(printf("[locate_contact] its: %d\n", its));
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
	  MD_DEBUG35(printf("[d>epsdFastR] d=%.15G\n", d));
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
  factori = 0.5*maxax[i]+epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));

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
#if defined(EDHE_FLEX) || defined(MD_BASIC_DT)
      /* check for convergence */
      if (delt < (epsd / maxddot))
	{
	  MD_DEBUG20(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
#else
      normddot = calcvecFNeigh(i, *t, t1, ddot, r1);
      //printf("normddot: %.15G\n", epsd/normddot);
      /* check for convergence */
     
      if (normddot!=0 && delt < (epsd / normddot))
	{
	  MD_DEBUG20(printf("convergence reached in %d iterations\n", its));
	  return 0;
	}
#endif
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
  
  factori = 0.5*maxax[i]+epsd;//sqrt(Sqr(axa[i])+Sqr(axb[i])+Sqr(axc[i]));
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
      else if (dold < epsd && d < epsd)
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
#ifdef MD_HE_PARALL
extern MPI_Datatype Particletype;
extern MPI_Datatype Eventtype;
#endif
#ifdef MD_HE_PARALL
void find_contact_parall(int na, int n, parall_event_struct *parall_event)
{
  double shift[3];
  double sigSq, tInt, d, b, vv, dv[3], dr[3], distSq, t, vecg[5];
  double t1, t2;
  int overlap;
#ifdef MD_PATCHY_HE
  double evtimeHC, evtime;
  int ac, bc, acHC, bcHC, collCode, collCodeOld;
#endif
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
  /* maxax[...] è il diametro dei centroidi dei due tipi
   * di ellissoidi */
  if (OprogStatus.targetPhi > 0)
    {
      sigSq = Sqr(max_ax(na)+max_ax(n)+OprogStatus.epsd);
    }
  else
    {
#ifdef MD_POLYDISP
      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#else
#ifdef EDHE_FLEX
      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#else
      if (na < parnumA && n < parnumA)
	sigSq = Sqr(maxax[na]+OprogStatus.epsd);
      else if (na >= parnumA && n >= parnumA)
	sigSq = Sqr(maxax[na]+OprogStatus.epsd);
      else
	sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#endif
#endif
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
      parall_event->t = -timbig;
      return;
    }
  MD_DEBUG32(printf("PREDICTING na=%d n=%d\n", na , n));
  if (vv==0.0)
    {
      if (distSq >= sigSq)
	{
	  parall_event->t = -timbig;
	  return;	 
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
#ifdef EDHE_FLEX
  if (OprogStatus.targetPhi <=0)
    {
      if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
	{
	  collCode = MD_EVENT_NONE;
	}
    }
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
	  evtime = vecg[4];
	  rxC = vecg[0];
	  ryC = vecg[1];
	  rzC = vecg[2];
	}
    }
  else
    {
      if (collCode == MD_EVENT_NONE)
	{
	  parall_event->t = -timbig;
	  return;
	}
    }

  t = evtime;
#else
  if (!locate_contact(na, n, shift, t1, t2, vecg))
    {
      parall_event->t = -timbig;
      return;
    }
  rxC = vecg[0];
  ryC = vecg[1];
  rzC = vecg[2];
  t = vecg[4];
#endif
  parall_event->t = t;
  parall_event->a = na;
  parall_event->b = n;
  parall_event->rC[0] = rxC;
  parall_event->rC[1] = ryC;
  parall_event->rC[2] = rzC;
  //printf("scheduling t=%.15G na=%d n=%d rC=%.15G %.15G %.15G\n", 
//	 t, na, n, rxC, ryC, rzC);
#ifdef MD_PATCHY_HE
  parall_event->sp[0] = ac;
  parall_event->sp[1] = bc;
  parall_event->sp[2] = collCode;
#endif
}
MPI_Status parall_status;
#ifdef MD_PATCHY_HE
extern struct LastBumpS *lastbump;
#else
extern int *lastbump;
#endif
extern double *lastcol;
void parall_slave_get_data(parall_pair_struct *parall_pair)
{
  int i, j, a, b;
  i = parall_pair->p[0];
  j = parall_pair->p[1];
  rx[i] = parall_pair->pos[0];
  ry[i] = parall_pair->pos[1];
  rz[i] = parall_pair->pos[2];
  rx[j] = parall_pair->pos[3];
  ry[j] = parall_pair->pos[4];
  rz[j] = parall_pair->pos[5];
  for (a = 0; a < 3; a++)
    {
      for (b = 0; b < 3; b++)
	{
	  R[i][a][b] = parall_pair->R[a*3+b];
	  R[j][a][b] = parall_pair->R[a*3+b+9]; 
	}
    }
  vx[i] = parall_pair->vels[0];
  vy[i] = parall_pair->vels[1];
  vz[i] = parall_pair->vels[2];
  vx[j] = parall_pair->vels[3];
  vy[j] = parall_pair->vels[4];
  vz[j] = parall_pair->vels[5];
  wx[i] = parall_pair->vels[6];
  wy[i] = parall_pair->vels[7];
  wz[i] = parall_pair->vels[8];
  wx[j] = parall_pair->vels[9];
  wy[j] = parall_pair->vels[10];
  wz[j] = parall_pair->vels[11];
  axa[i] = parall_pair->axes[0];
  axb[i] = parall_pair->axes[1];
  axc[i] = parall_pair->axes[2];
  axa[j] = parall_pair->axes[3];
  axb[j] = parall_pair->axes[4];
  axc[j] = parall_pair->axes[5];
  maxax[i] = parall_pair->axes[6];
  maxax[j] = parall_pair->axes[7];
  inCell[0][i] = parall_pair->cells[0];
  inCell[1][i] = parall_pair->cells[1];
  inCell[2][i] = parall_pair->cells[2];
  inCell[0][j] = parall_pair->cells[3];
  inCell[1][j] = parall_pair->cells[4];
  inCell[2][j] = parall_pair->cells[5];
#ifdef MD_PATCHY_HE
  lastbump[i].mol = parall_pair->lastbump[0];
  lastbump[i].at = parall_pair->lastbump[1];
  lastbump[i].type = parall_pair->lastbump[2];
  lastbump[j].mol = parall_pair->lastbump[3];
  lastbump[j].at = parall_pair->lastbump[4];
  lastbump[j].type = parall_pair->lastbump[5];
#else
  lastbump[i] = parall_pair->lastbump[0];
  lastbump[i] = parall_pair->lastbump[1];
#endif
  Oparams.time = parall_pair->time[0];
  nextNNLrebuild = parall_pair->time[1];
  lastcol[i] = parall_pair->time[2];
  lastcol[j] = parall_pair->time[3];
  atomTime[i] = parall_pair->atomTime[0];
  atomTime[j] = parall_pair->atomTime[1];
#ifdef MD_ASYM_ITENS
  angM[i] = parall_pair->angM[0];
  angM[j] = parall_pair->angM[1];
  sintheta0[i] = parall_pair->sintheta0[0];
  sintheta0[j] = parall_pair->sintheta0[1];
  costheta0[i] = parall_pair->costheta0[0];
  costheta0[j] = parall_pair->costheta0[1];
  phi0[i] = parall_pair->phi0[0];
  phi0[j] = parall_pair->phi0[1];
  psi0[i] = parall_pair->psi0[0];
  psi0[j] = parall_pair->psi0[1];
  for (a = 0; a < 3; a++)
    {
      for (b = 0; b < 3; b++)
	{
	  RM[i][a][b] = parall_pair->RM[a*3+b];
	  RM[j][a][b] = parall_pair->RM[a*3+b+9]; 
	}
    }
#endif
}
void parall_set_data_terminate(parall_pair_struct *parall_pair)
{
  parall_pair->p[0] = -1;
  parall_pair->p[1] = -1;
}
void parall_set_data(int i, int j, parall_pair_struct *parall_pair)
{
  int a, b;
  parall_pair->p[0] = i;
  parall_pair->p[1] = j;
  parall_pair->pos[0] = rx[i];
  parall_pair->pos[1] = ry[i];
  parall_pair->pos[2] = rz[i];
  parall_pair->pos[3] = rx[j];
  parall_pair->pos[4] = ry[j];
  parall_pair->pos[5] = rz[j];
  for (a = 0; a < 3; a++)
    {
      for (b = 0; b < 3; b++)
	{
	  parall_pair->R[a*3+b] = R[i][a][b];
	  parall_pair->R[a*3+b+9] = R[j][a][b]; 
	}
    }
  parall_pair->vels[0] = vx[i];
  parall_pair->vels[1] = vy[i];
  parall_pair->vels[2] = vz[i];
  parall_pair->vels[3] = vx[j];
  parall_pair->vels[4] = vy[j];
  parall_pair->vels[5] = vz[j];
  parall_pair->vels[6] = wx[i];
  parall_pair->vels[7] = wy[i];
  parall_pair->vels[8] = wz[i];
  parall_pair->vels[9] = wx[j];
  parall_pair->vels[10] = wy[j];
  parall_pair->vels[11] = wz[j];
  parall_pair->axes[0] = axa[i];
  parall_pair->axes[1] = axb[i];
  parall_pair->axes[2] = axc[i];
  parall_pair->axes[3] = axa[j];
  parall_pair->axes[4] = axb[j];
  parall_pair->axes[5] = axc[j];
  parall_pair->axes[6] = maxax[i];
  parall_pair->axes[7] = maxax[j];
  parall_pair->cells[0] = inCell[0][i];
  parall_pair->cells[1] = inCell[1][i];
  parall_pair->cells[2] = inCell[2][i];
  parall_pair->cells[3] = inCell[0][j];
  parall_pair->cells[4] = inCell[1][j];
  parall_pair->cells[5] = inCell[2][j];
#ifdef MD_PATCHY_HE
  parall_pair->lastbump[0] = lastbump[i].mol;
  parall_pair->lastbump[1] = lastbump[i].at;
  parall_pair->lastbump[2] = lastbump[i].type;
  parall_pair->lastbump[3] = lastbump[j].mol;
  parall_pair->lastbump[4] = lastbump[j].at;
  parall_pair->lastbump[5] = lastbump[j].type;
#else
  parall_pair->lastbump[0] = lastbump[i];
  parall_pair->lastbump[1] = lastbump[j];
#endif
  parall_pair->time[0] = Oparams.time;
  parall_pair->time[1] = nextNNLrebuild;
  parall_pair->time[2] = lastcol[i];
  parall_pair->time[3] = lastcol[j];
  parall_pair->atomTime[0] = atomTime[i];
  parall_pair->atomTime[1] = atomTime[j];
#ifdef MD_ASYM_ITENS
  parall_pair->angM[0] = angM[i];
  parall_pair->angM[1] = angM[j];
  parall_pair->sintheta0[0] = sintheta0[i];
  parall_pair->sintheta0[1] = sintheta0[j];
  parall_pair->costheta0[0] = costheta0[i];
  parall_pair->costheta0[1] = costheta0[j];
  parall_pair->phi0[0] = phi0[i];
  parall_pair->phi0[1] = phi0[j];
  parall_pair->psi0[0] = psi0[i];
  parall_pair->psi0[1] = psi0[j];
  for (a = 0; a < 3; a++)
    {
      for (b = 0; b < 3; b++)
	{
	  parall_pair->RM[a*3+b] = RM[i][a][b];
	  parall_pair->RM[a*3+b+9] = RM[j][a][b]; 
	}
    }

#endif
}
void parall_get_data_and_schedule(parall_event_struct parall_event)
{
  int na, n;
#ifdef MD_PATCHY_HE
  int collCode, ac, bc;
#endif
  double t;
  t = parall_event.t;
  if (t == - timbig)
    return;
  na = parall_event.a;
  n = parall_event.b;
  rxC = parall_event.rC[0];
  ryC = parall_event.rC[1];
  rzC = parall_event.rC[2];
#ifdef MD_PATCHY_HE
  ac = parall_event.sp[0];
  bc = parall_event.sp[1];
  collCode = parall_event.sp[2];
  ScheduleEventBarr (na, n,  ac, bc, collCode, t);
#else
  ScheduleEvent (na, n, t);
#endif
}
/* PredicEventNNL parallelizzata */
const int iwtagEvent = 1, iwtagPair = 2;

void PredictEventNNL_PARALL(int na, int nb) 
{
  int i, signDir[NDIM]={0,0,0}, evCode, k, n;
  double  tm[NDIM];
  char msgtype;
  int npr, iriceve, ndone=0;
#ifdef MD_HE_PARALL
  int missing, num_work_request, njob, njob2i[512];
#endif

 if (vz[na] != 0.0) 
    {
      if (vz[na] > 0.0) 
	signDir[2] = 0;/* direzione positiva */
      else 
	signDir[2] = 1;/* direzione negativa */
#ifdef MD_LXYZ
#ifdef MD_EDHEFLEX_WALL
      if (OprogStatus.hardwall && ((signDir[2]==0 && inCell[2][na]==cellsz-1) || (signDir[2]==1 && inCell[2][na]==0)))
	tm[2] = timbig;
      else
	tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L[2] /
		 cellsz - rz[na] - L2[2]) / vz[na];
#else
      tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L[2] /
	       cellsz - rz[na] - L2[2]) / vz[na];
#endif
#else
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
#endif    
    } 
  else 
    tm[2] = timbig;
  if (vx[na] != 0.0) 
    {
#ifdef MD_LXYZ
      if (vx[na] > 0.0) 
	signDir[0] = 0;/* direzione positiva */
      else 
	signDir[0] = 1;/* direzione negativa */
      tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L[0] /
	       cellsx - rx[na] - L2[0]) / vx[na];
#else
      if (vx[na] > 0.0) 
	signDir[0] = 0;/* direzione positiva */
      else 
	signDir[0] = 1;/* direzione negativa */
      tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L /
	       cellsx - rx[na] - L2) / vx[na];
#endif
    } 
  else 
    tm[0] = timbig;
  
  if (vy[na] != 0.) 
    {
#ifdef MD_LXYZ
      if (vy[na] > 0.) 
	signDir[1] = 0;
      else 
	signDir[1] = 1;
      tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L[1] /
	       cellsy - ry[na] - L2[1]) / vy[na];
#else
      if (vy[na] > 0.) 
	signDir[1] = 0;
      else 
	signDir[1] = 1;
      tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L /
	       cellsy - ry[na] - L2) / vy[na];
#endif
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
  /* NOTA: nel caso di attraversamento di una cella non deve predire le collisioni (visto che stiamo in PredictEventNNL
     e si stanno usando le neighborurs lists */
  if (nb >= ATOM_LIMIT)
    return;
  MD_DEBUG32(printf("nebrTab[%d].len=%d\n", na, nebrTab[na].len));
  if (my_rank == 0)
    {
      num_work_request = 0;
      for (i=0; i < nebrTab[na].len; i++)
	{
	  n = nebrTab[na].list[i]; 
	  if (!(n != na && n!=nb && (nb >= -1 || n < na)))
	    continue;
	  if (numOfProcs == 1)
	    {
	      parall_set_data(na, n, &parall_pair);
	      find_contact_parall(parall_pair.p[0], parall_pair.p[1], &parall_event);
	      parall_get_data_and_schedule(parall_event);
	    }
	  //printf("to locate: n=%d\n", n);
	  njob2i[num_work_request] = i;
	  num_work_request++; 
	}	
      if (numOfProcs == 1)
	return;
      if (num_work_request == 0)
	{
	  for(npr = 1; npr < numOfProcs; npr++)
	    {
	      //printf("sending termination to %d\n", npr);
	      msgtype = 'D';
	      MPI_Send(&msgtype, 1, MPI_CHAR, npr, iwtagPair, MPI_COMM_WORLD); 
	      //MPI_Send(&parall_pair, 1, Particletype, npr, iwtagPair, MPI_COMM_WORLD);
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	  return;
	}
      /* master process (rank=0) distributes jobs */
      for (njob=1; njob < numOfProcs; njob++)
	{
	  if (njob > num_work_request)
       	    break;
	  n = nebrTab[na].list[njob2i[njob-1]]; 
	  //printf(">>> num_work_request: %d numofproc:%d sending to %d\n", num_work_request, numOfProcs, njob);
	  //printf("na=%d sending n=%d to %d\n",na, n, njob);
	  parall_set_data(na ,n, &parall_pair);
	  //printf("sending to %d\n", njob);
	  msgtype = 'F';
	  MPI_Send(&msgtype, 1, MPI_CHAR, njob, iwtagPair, MPI_COMM_WORLD); 
	  MPI_Send(&parall_pair, 1, Particletype, njob, iwtagPair, MPI_COMM_WORLD);
	  //printf("sent!!!!\n");
	}    
      //printf("????????????????\n");
      //printf("num_work_request=%d numofproc=%d\n", num_work_request, numOfProcs);
      for (njob=numOfProcs; njob <= num_work_request; njob++)
	 {
	   //printf("my_rank=%d njob = %d\n", my_rank, njob);
	   n = nebrTab[na].list[njob2i[njob-1]]; 
	   MPI_Recv(&parall_event, 1, Eventtype, MPI_ANY_SOURCE, iwtagEvent, MPI_COMM_WORLD, &parall_status);
	   ndone++;
	   //printf("§§§§§§ ndone=%d\n", ndone);
	   iriceve = parall_status.MPI_SOURCE;
	   //printf(">>>na=%d sending n=%d to %d\n",na, n, iriceve);
	   parall_get_data_and_schedule(parall_event);
	   parall_set_data(na, n, &parall_pair);
 	   msgtype = 'F';
 	   MPI_Send(&msgtype, 1, MPI_CHAR, iriceve, iwtagPair, MPI_COMM_WORLD); 
	   MPI_Send(&parall_pair, 1, Particletype, iriceve, iwtagPair, MPI_COMM_WORLD);
	 }
      for (missing = ndone+1; missing <= num_work_request; missing++)
	{
	  //printf("my_rank=%d receiving missing = %d\n", my_rank, missing);
	  MPI_Recv(&parall_event, 1, Eventtype, MPI_ANY_SOURCE, iwtagEvent, MPI_COMM_WORLD, &parall_status);
	  iriceve = parall_status.MPI_SOURCE;
	  parall_get_data_and_schedule(parall_event);
	}
      //printf("mmm\n");
      for(npr = 1; npr < numOfProcs; npr++)
	{
	  //printf("RANK0000 terminating process=%d\n", npr);
	  //parall_set_data_terminate(&parall_pair);
  	  msgtype = 'D';
	  MPI_Send(&msgtype, 1, MPI_CHAR, npr, iwtagPair, MPI_COMM_WORLD); 
  	  //MPI_Send(&parall_pair, 1, Particletype, npr, iwtagPair, MPI_COMM_WORLD);
	}
    }
      
  MPI_Barrier(MPI_COMM_WORLD);
  //printf("FINISHED rank=%d\n", my_rank);
}
#else
#ifdef MD_SPHERICAL_WALL
extern int sphWall, sphWallOuter;
extern void locate_spherical_wall(int na, int outer);
#endif
extern int use_bounding_spheres(int na, int n);
#ifdef MD_MULTIPLE_LL
void PredictEventNNL_MLL(int na, int nb);
#endif
#if ED_PARALL_DD
extern int dd_is_virtual(int i);
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
#ifdef MD_EDHEFLEX_WALL
  int nplane=-1;
#endif
#endif
#ifdef MD_ABSORPTION
  int hwcell;
#endif
#ifdef ED_PARALL_DD
  if (dd_is_virtual(na))
    return;
#endif
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
      PredictEventNNL_MLL(na, nb);
      return;
    }
#endif
#ifdef MD_SPHERICAL_WALL
  if (na==sphWall|| nb==sphWall)
    return;
  if (na==sphWallOuter|| nb==sphWallOuter)
    return;
#endif
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=0;
#endif

  if (vz[na] != 0.0) 
    {
      if (vz[na] > 0.0) 
	signDir[2] = 0;/* direzione positiva */
      else 
	signDir[2] = 1;/* direzione negativa */
#ifdef MD_LXYZ
#if defined(MD_EDHEFLEX_WALL) 
      if (OprogStatus.hardwall && ((signDir[2]==0 && inCell[2][na]==cellsz-1) || (signDir[2]==1 && inCell[2][na]==0)))
	tm[2] = timbig;
      else
	tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L[2] /
		 cellsz - rz[na] - L2[2]) / vz[na];
#else
      tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L[2]/
	       cellsz - rz[na] - L2[2]) / vz[na];
#endif
#else
#if defined(MD_EDHEFLEX_WALL) 
      if (OprogStatus.hardwall && ((signDir[2]==0 && inCell[2][na]==cellsz-1) || (signDir[2]==1 && inCell[2][na]==0)))
	tm[2] = timbig;
      else
	tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L /
		 cellsz - rz[na] - L2) / vz[na];
#else
      tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L /
	       cellsz - rz[na] - L2) / vz[na];
#endif
#endif
    } 
  else 
    tm[2] = timbig;
  if (vx[na] != 0.0) 
    {
      if (vx[na] > 0.0) 
	signDir[0] = 0;/* direzione positiva */
      else 
	signDir[0] = 1;/* direzione negativa */
#ifdef MD_LXYZ
      tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L[0] /
	       cellsx - rx[na] - L2[0]) / vx[na];
#else
      tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L /
	       cellsx - rx[na] - L2) / vx[na];
#endif    
    } 
  else 
    tm[0] = timbig;
  
  if (vy[na] != 0.) 
    {
      if (vy[na] > 0.) 
	signDir[1] = 0;
      else 
	signDir[1] = 1;
#ifdef MD_LXYZ
      tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L[1] /
	       cellsy - ry[na] - L2[1]) / vy[na];
#else
      tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L /
	       cellsy - ry[na] - L2) / vy[na];
#endif
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
#ifdef MD_EDHEFLEX_WALL
  if (OprogStatus.hardwall)
    {
#if defined(MD_ABSORPTION) 
      if (OprogStatus.bufHeight > 0.0)
	{
	  if (vz[na] != 0.0)
	    {
#ifdef MD_LXYZ
	      hwcell = (L[2]-OprogStatus.bufHeight)*cellsz/L[2];
#else
	      hwcell = (L-OprogStatus.bufHeight)*cellsz/L;
#endif
#if 1
	      if (hwcell-inCell[2][na] < 2)
		{
		  /* the semi-permeable plane is just one (nplane=0) */
		  if (locateHardWall(na, 0, Oparams.time+tm[k], vecg, 2))
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
      if (inCell[2][na] == 0)
	nplane = 0;
      else if (inCell[2][na] == cellsz-1)
	nplane = 1;
      if (nplane!=-1 && locateHardWall(na, nplane, Oparams.time+tm[k], vecg, 1))
	{
	  rxC = vecg[0];
	  ryC = vecg[1];
	  rzC = vecg[2];
	  MD_DEBUG35(printf("Located Contact with WALL rC=%f %f %f time=%.15G i=%d\n", rxC, ryC, rzC, vecg[4], na));
	  MD_DEBUG35(printf("r=%f %f %f\n", rx[na], ry[na], rz[na]));
	  ScheduleEventBarr (na, ATOM_LIMIT + nplane, 0, 0, MD_WALL, vecg[4]);
	}
      else 
	ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]);
    }
  else
    ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]);
#else
  ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]);
#endif

#ifdef MD_SPHERICAL_WALL
  locate_spherical_wall(na, 0);
  locate_spherical_wall(na, 1);
#endif 
  /* NOTA: nel caso di attraversamento di una cella non deve predire le collisioni (visto che in tal caso stiamo 
     usando le NNL */

  if (nb >= ATOM_LIMIT+2*NDIM)
    {
 #ifdef MD_SOLVENT_NOHW
      if (typeOfPart[na]==4)
	OprogStatus.hardwall=1;
#endif
      return;
    }
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
#ifdef ED_PARALL_DD
      /* avoid to predict twice collision time of updated partcles */

      if (doing_rollback_load && rb_since_load_changed[n] && n >= na)
	continue;
#endif
      //
      // for (kk=0; kk < 3; kk++)
      //	shift[kk] = nebrTab[na].shift[i][kk];
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
      /* maxax[...] è il diametro dei centroidi dei due tipi
       * di ellissoidi */
      if (use_bounding_spheres(na, n))
	{
	  if (OprogStatus.targetPhi > 0)
	    {
	      sigSq = Sqr(max_ax(na)+max_ax(n)+OprogStatus.epsd);
	    }
	  else
	    {
#ifdef MD_POLYDISP
	      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#else
#ifdef EDHE_FLEX
	      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
	      //printf("max=(%d)%.15G (%d)%.15G\n", n, maxax[n], na, maxax[na]);
#else
	      if (na < parnumA && n < parnumA)
		sigSq = Sqr(maxax[na]+OprogStatus.epsd);
	      else if (na >= parnumA && n >= parnumA)
		sigSq = Sqr(maxax[na]+OprogStatus.epsd);
	      else
		sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#endif
#endif
	    }
	  MD_DEBUG32(printf("n=%d na=%d maxaxes=%f %f sigSq: %f\n", n, na, maxax[n], maxax[na], sigSq));
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
	}
      else
	{
	  t1 = 0.0;
	  t2 = timbig;
	}     
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
#ifdef EDHE_FLEX
      if (OprogStatus.targetPhi <=0 && typesArr[typeOfPart[na]].nspots > 0 && typesArr[typeOfPart[n]].nspots > 0)
	{
	  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
	    {
	      collCode = MD_EVENT_NONE;
	    }
	}
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
#ifdef MD_SAVE_DISTANCE
	  printf("found collision exiting...\n");
	  //exit(-1);
#endif
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
#ifdef MD_PATCHY_HE
      MD_DEBUG32(printf("Scheduling collision between %d and %d ac=%d bc=%d at t=%.15G\n", na, n, ac, bc, t));
      //printf("Scheduling collision between %d and %d ac=%d bc=%d at t=%.15G\n", na, n, ac, bc, t);
      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
#else
      ScheduleEvent (na, n, t);
#endif
    }
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=1;
#endif
}
#endif
#ifdef MD_PATCHY_HE
extern int locate_contact_neigh_plane_parall_sp(int i, double *evtime, double t2);
#endif
#ifdef MD_HANDLE_INFMASS
extern int is_infinite_mass(int i);
#endif
#ifdef EDHE_FLEX
extern int *is_a_sphere_NNL;
int locateNNLSP;
#endif
void updrebuildNNL(int na)
{
  /* qui ricalcola solo il tempo di collisione dell'ellisoide na-esimo con 
   * la sua neighbour list */
  double vecg[5];
  double nnltime1, nnltime2;
  int ip;
#ifdef MD_SPHERICAL_WALL
  if (na==sphWall || na==sphWallOuter)
    return;
#endif
#ifdef MD_NNLPLANES
#ifdef MD_PATCHY_HE
  nnltime1 = timbig;
  MD_DEBUG33(printf("qui\n"));
  MD_DEBUG33(printf("posEll=%f %f %f posNNL %f %f %f \n", rx[na], ry[na], rz[na], nebrTab[na].r[0],nebrTab[na].r[1],
		    nebrTab[na].r[2]));
#ifdef EDHE_FLEX
  if (is_infinite_Itens(na) && is_infinite_mass(na))
    {
      return;
    }
  if (is_infinite_mass(na) && is_a_sphere_NNL[na])
    {
      return;    
    }
#endif
#ifdef EDHE_FLEX
  if ((OprogStatus.targetPhi <= 0.0 && typesArr[typeOfPart[na]].nspots > 0) || 
      (OprogStatus.useNNL==3 || OprogStatus.useNNL==4))
    {
      locateNNLSP=1;
      if (!locate_contact_neigh_plane_parall_sp(na, &nnltime1, timbig))
	{
#ifdef MD_HANDLE_INFMASS
	  if (is_infinite_mass(na))
	    nnltime1 = timbig; 
	  else
	    {
	      printf("[ERROR] failed to find escape time for sticky spots na=%d\n", na);
	      printf("i=%d nspots=%d\n", na,  typesArr[typeOfPart[na]].nspots);
	      exit(-1);
	    }
#else
	  printf("[ERROR] failed to find escape time for sticky spots na=%d\n", na);
	  printf("i=%d nspots=%d\n", na,  typesArr[typeOfPart[na]].nspots);
	  exit(-1);
#endif
	}
      locateNNLSP=0;
    }
  else 
    nnltime1 = timbig;
#else
  if (OprogStatus.targetPhi <= 0.0)
    {
      if (!locate_contact_neigh_plane_parall_sp(na, &nnltime1, timbig))
	{
	  printf("[ERROR] failed to find escape time for sticky spots na=%d\n", na);
	  exit(-1);
	}
    }
  else 
    nnltime1 = timbig;

#endif
  MD_DEBUG32(printf("nexttime=%.15G\n", nebrTab[na].nexttime));
#else
  nnltime1 = timbig;
#endif
  nnltime2 = timbig;
  if (OprogStatus.useNNL==1 || OprogStatus.useNNL==2)
    {
      if (OprogStatus.paralNNL)
	{
	  //switch_to_HE(na);
	  if (!locate_contact_neigh_plane_parall(na, &nnltime2, nnltime1))
	    {
#ifndef MD_PATCHY_HE
	      printf("[ERROR] failed to find escape time for ellipsoid N. %d\n", na);
	      exit(-1);
#endif
	    }
	  //back_to_SQ(na);
	}
      else
	{
	  nnltime2 = nnltime1;
#ifdef EDHE_FLEX
	  if (is_sphere(na))
	    {
	      calc_grad_and_point_plane_all(na, gradplane_all, rBall);
	      if (!locate_contact_neigh_plane_HS(na, &nnltime2, nnltime1))
		{
		  /* do nothing */
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
#else
	  for (ip = 0; ip < 6; ip++)
	    {
	      if (!locate_contact_neigh_plane(na, vecg, ip, nnltime1))
		continue;
	      if (vecg[4] < nnltime2)
		nnltime2 = vecg[4];
	    }
#endif
	}
    }
 if (OprogStatus.useNNL==1 || OprogStatus.useNNL == 2)
   nebrTab[na].nexttime = min(nnltime1, nnltime2);
 else
   nebrTab[na].nexttime = nnltime1;
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
 
  MD_DEBUG32(printf("updneigh REBUILD i=%d t=%.15G nextNNLrebuild=%.15G\n", na, nebrTab[na].nexttime, nextNNLrebuild));
  if (nebrTab[na].nexttime < nextNNLrebuild)
    nextNNLrebuild = nebrTab[na].nexttime;
}
void updAllNNL()
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    updrebuildNNL(i);
}
#ifdef EDHE_FLEX
void body2labR(int i, double xp[], double x[], double *rO, double **R);
#endif
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
#ifdef EDHE_FLEX
  double xl[3];
  int typena;
#endif
  MD_DEBUG33(printf("nextNNLupdate...\n"));
  MD_DEBUG33(printf("posEll=%f %f %f posNNL %f %f %f \n", rx[na], ry[na], rz[na], nebrTab[na].r[0],nebrTab[na].r[1],
		    nebrTab[na].r[2]));
#ifndef MD_NNLPLANES
#ifdef EDHE_FLEX
  typena = typeOfPart[na];
  if (OprogStatus.targetPhi > 0.0)
    {
      nebrTab[na].axa = OprogStatus.rNebrShell*axa[na];
      nebrTab[na].axb = OprogStatus.rNebrShell*axb[na];
      nebrTab[na].axc = OprogStatus.rNebrShell*axc[na];
    }
  else
    {
      nebrTab[na].axa = OprogStatus.rNebrShell*typesArr[typena].sax[0];
      nebrTab[na].axb = OprogStatus.rNebrShell*typesArr[typena].sax[1];
      nebrTab[na].axc = OprogStatus.rNebrShell*typesArr[typena].sax[2];
    }
#else
  nebrTab[na].axa = OprogStatus.rNebrShell*axa[na];
  nebrTab[na].axb = OprogStatus.rNebrShell*axb[na];
  nebrTab[na].axc = OprogStatus.rNebrShell*axc[na];
#endif
#ifdef EDHE_FLEX
  if (OprogStatus.targetPhi > 0.0)
    {
      DelDist = max3(nebrTab[na].axa,nebrTab[na].axb,nebrTab[na].axc) -
	max3(axa[na],axb[na],axc[na]);
    }
  else
    {
      DelDist = max3(nebrTab[na].axa,nebrTab[na].axb,nebrTab[na].axc) -
	max3(typesArr[typena].sax[0],typesArr[typena].sax[1],typesArr[typena].sax[2]);
    }
#else
  DelDist = max3(nebrTab[na].axa,nebrTab[na].axb,nebrTab[na].axc) -
    max3(axa[na],axb[na],axc[na]);
#endif
#else
  if (OprogStatus.targetPhi > 0.0)
    {
#ifdef EDHE_FLEX
      nnlfact = axa[na]/typesArr[typeOfPart[na]].sax[0];
#else
#ifdef MD_POLYDISP
      nnlfact = axa[na]/axaP[na];
#else
      nnlfact = axa[na]/Oparams.a[na<Oparams.parnumA?0:1];
#endif
#endif
      /* NOTA 12/04/10:
	 ppsax sono in generale i semiassi per le NNL ottimizzate, ossia che tengono conto
	 anche degli spot, tuttavia durante la crescita gli spot vengono disattivati,
	 quindi si possono usare tranquillamente i semiassi degli ellissoidi polidispersi
	 duraten la crescita (axa[...],axb[...],axc[...])*/
      nebrTab[na].axa = nnlfact*OprogStatus.rNebrShell+axa[na];
      nebrTab[na].axb = nnlfact*OprogStatus.rNebrShell+axb[na];
      nebrTab[na].axc = nnlfact*OprogStatus.rNebrShell+axc[na];
      DelDist = nnlfact*OprogStatus.rNebrShell;
    }
  else
    {
#ifdef EDHE_FLEX
      typena = typeOfPart[na]; 
      nebrTab[na].axa = OprogStatus.rNebrShell+typesArr[typena].ppsax[0];
      nebrTab[na].axb = OprogStatus.rNebrShell+typesArr[typena].ppsax[1];
      nebrTab[na].axc = OprogStatus.rNebrShell+typesArr[typena].ppsax[2];
#else
      nebrTab[na].axa = OprogStatus.rNebrShell+axa[na];
      nebrTab[na].axb = OprogStatus.rNebrShell+axb[na];
      nebrTab[na].axc = OprogStatus.rNebrShell+axc[na];
#endif
      DelDist = OprogStatus.rNebrShell;
    }
#endif
  DelDist += distBuf;
  MD_DEBUG31(printf("DelDist=%.15G\n", DelDist));
#ifdef EDHE_FLEX
  typena = typeOfPart[na]; 
#endif
#ifndef EDHE_FLEX
  nebrTab[na].r[0] = rx[na];
  nebrTab[na].r[1] = ry[na];
  nebrTab[na].r[2] = rz[na];
#endif
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(na, 0, RtB, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(na, 0, RtB, Omega);
#endif
  for (i1 = 0; i1 < 3; i1++)
    for (i2 = 0; i2 < 3; i2++)
      nebrTab[na].R[i1][i2] = RtB[i1][i2];
  /* calcola il tempo a cui si deve ricostruire la NNL */
#ifdef EDHE_FLEX
#ifdef MD_EDHEFLEX_OPTNNL
  if (OprogStatus.optnnl)
    {
      body2labR(na, typesArr[typena].ppr, xl, NULL, nebrTab[na].R);
      nebrTab[na].r[0] = rx[na]+xl[0];
      nebrTab[na].r[1] = ry[na]+xl[1];
      nebrTab[na].r[2] = rz[na]+xl[2];
    }
  else
    {
      nebrTab[na].r[0] = rx[na];
      nebrTab[na].r[1] = ry[na];
      nebrTab[na].r[2] = rz[na];
    }
  //if (Oparams.curStep > 2690 && na <= 4)
  //printf("r[%d]=%f %f %f nebrTab.r=%f %f %f xl=%f %f %f\n", na, rx[na], ry[na], rz[na], nebrTab[na].r[0], nebrTab[na].r[1], nebrTab[na].r[2], xl[0], xl[1], xl[2]);
#else
  nebrTab[na].r[0] = rx[na];
  nebrTab[na].r[1] = ry[na];
  nebrTab[na].r[2] = rz[na];
#endif
#endif
  MD_DEBUG31(printf("BUILDING NNL FOR i=%d\n",na));
#ifdef MC_SIMUL
  nebrTab[na].len=0; /* WARNING: check this for MC sims! */
  return;
#endif
#ifdef MD_NNLPLANES
#ifdef MD_PATCHY_HE
#ifdef EDHE_FLEX
  nnltime1 = timbig;
  if ((OprogStatus.targetPhi <= 0.0 && typesArr[typeOfPart[na]].nspots > 0) ||
      (OprogStatus.useNNL==3 || OprogStatus.useNNL==4))
    {
      locateNNLSP = 1;
      if (!locate_contact_neigh_plane_parall_sp(na, &nnltime1, timbig))
	{
	  printf("[ERROR] failed to find escape time for sticky spots\n");
	  printf("typeOfPart[%d]=%d\n", na, typeOfPart[na]);
	  exit(-1);
	}
      locateNNLSP = 0;
    }
  else
    nnltime1 = timbig;
#else
  nnltime1 = timbig;
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

#endif
  MD_DEBUG32(printf("[nextNNLupdate] nexttime=%.15G\n", nebrTab[na].nexttime));
#else
  nnltime1 = timbig; 
#endif
  nnltime2 = timbig;
  if (OprogStatus.useNNL==1 || OprogStatus.useNNL==2)
    {
      if (OprogStatus.paralNNL)
	{
	  //switch_to_HE(na);
	  if (!locate_contact_neigh_plane_parall(na, &nnltime2, nnltime1))
	    {
#ifndef MD_PATCHY_HE
	      printf("[ERROR] failed to find escape time for ellipsoid N. %d\n", na);
	      exit(-1);
#endif
	    }
	  //back_to_SQ(na);
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
    }
  if (OprogStatus.useNNL==1 || OprogStatus.useNNL==2)
    nebrTab[na].nexttime = min(nnltime1, nnltime2);
  else
    nebrTab[na].nexttime = nnltime1;
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
#ifdef EDHE_FLEX
extern int is_in_ranges(int A, int B, int nr, rangeStruct* r);
int may_interact_core(int i, int j)
{
  int typei, typej;
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  if (!typesArr[typei].ignoreCore && 
      !typesArr[typej].ignoreCore)
   return 1;
 else
   return 0; 
}
#ifdef EDHE_FLEX
extern int *nbondsFlexS;
extern int get_linked_list_type(int type1, int type2);
#endif
int may_interact_spots(int i, int j)
{
  int type1, type2, ni, nl, numll;
  type1 = typeOfPart[i];
  type2 = typeOfPart[j];
  if (typesArr[type1].nspots == 0 || typesArr[type2].nspots == 0)
    return 0;
#ifdef EDHE_FLEX
  if (OprogStatus.optbm)
    {
      nl = get_linked_list_type(type1, type2);
      if (nbondsFlexS[nl] > 0)
	return 1;
    }
#endif
  for (ni = 0; ni < Oparams.ninters; ni++)
    {
      if (is_in_ranges(type1, intersArr[ni].type1, intersArr[ni].nr1, intersArr[ni].r1) && 
	  is_in_ranges(type2, intersArr[ni].type2, intersArr[ni].nr2, intersArr[ni].r2))
	{
	  return 1;
	}	
      else if (is_in_ranges(type2, intersArr[ni].type1, intersArr[ni].nr1, intersArr[ni].r1) && 
	       is_in_ranges(type1, intersArr[ni].type2, intersArr[ni].nr2, intersArr[ni].r2))
	{
	  return 1;
	}
    }
  if (Oparams.nintersIJ > 0)
    {
      for (ni = 0; ni < Oparams.nintersIJ; ni++)
	{
	  if (is_in_ranges(i, intersArrIJ[ni].i, intersArrIJ[ni].nr1, intersArrIJ[ni].r1) && 
	      is_in_ranges(j, intersArrIJ[ni].j, intersArrIJ[ni].nr2, intersArrIJ[ni].r2))
	    {
	      return 1;
	    }	
	  else if (is_in_ranges(j, intersArrIJ[ni].i, intersArrIJ[ni].nr1, intersArrIJ[ni].r1) && 
		   is_in_ranges(i, intersArrIJ[ni].j, intersArrIJ[ni].nr2, intersArrIJ[ni].r2))
	    {
	      return 1;
	    }
	}
    }
  return 0;
}
int may_interact_all(int i, int j)
{
  /* infinite mass spheres can not interact! */
  if (is_infinite_mass(i) && is_a_sphere_NNL[i] &&
      is_infinite_mass(j) && is_a_sphere_NNL[j])
    return 0;
  if (may_interact_core(i, j))
    return 1;
  if (may_interact_spots(i, j))
    return 1;
  return 0;
}
#endif
void check_nnl_size(int na)
{
  int i;
  double dblsize;
  int nnlsize;

  if (nebrTab[na].len < OprogStatus.nebrTabFac - 1)
    return;
  
  dblsize = (double)OprogStatus.nebrTabFac;
  while (((int)dblsize)==OprogStatus.nebrTabFac)
    {
      dblsize *= 1.10; /* try 10% increment */
    }

  OprogStatus.nebrTabFac = (int) dblsize;
  for (i=0; i < Oparams.parnum; i++)
    {
      nnlsize = OprogStatus.nebrTabFac*sizeof(int); 
      nebrTab[i].list = (int*)realloc(nebrTab[i].list, nnlsize);
#ifdef MD_ABSORPTION
      listtmp = (int*) realloc(listtmp, nnlsize);
#endif
      if (nebrTab[i].list==NULL)
	{
	  printf("[CRITICAL ERROR] neighbor lists size exceed memory allocation of %d elements\n", OprogStatus.nebrTabFac);
	  printf("and realloc() to enlarge nebrTabFac array failed\n");
	  exit(-1);
	}
    }
  printf("[INFO]  neighbor lists exceed allocated memory I have just enlarged it with realloc()\n");
  printf("new size (nebrTabFac) is %d\n", OprogStatus.nebrTabFac);
}
#ifdef MD_MULTIPLE_LL
extern void BuildAllNNL_MLL_one(int i);
#endif
void BuildNNL(int na) 
{
  double shift[NDIM];
  int kk;
  double dist;
  int *inCellL[3], *cellListL;
#ifndef MD_NNLPLANES
  double vecgsup[8], alpha;
#endif
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cels[NDIM];*/
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
  nebrTab[na].len = 0;
  MD_DEBUG32(printf("Building NNL...\n"));
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
      BuildAllNNL_MLL_one(na);
      return;
    }
#endif
#ifdef MD_SPHERICAL_WALL
  if (na==sphWall || na==sphWallOuter)
    return;
#endif

#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=0;
#endif
 
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (kk=0; kk < 3; kk++)
    shift[kk] = 0;
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_EDHEFLEX_OPTNNL
  if (OprogStatus.optnnl)
    {
      for (k=0; k < 3; k++)
	inCellL[k] = inCell_NNL[k];
      cellListL = cellList_NNL;
    }
  else
    {
      for (k=0; k < 3; k++)
	inCellL[k] = inCell[k];
      cellListL = cellList;
    }
#else
  for (k=0; k < 3; k++)
    inCellL[k] = inCell[k];
  cellListL = cellList;
#endif
#if defined(MD_EDHEFLEX_WALL)
  /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (OprogStatus.hardwall)
    {
      if (inCellL[2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (inCellL[2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
    }
#endif
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCellL[2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodic boundary condition along z if gravitational
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
	  jY = inCellL[1][na] + iY;    
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
	      jX = inCellL[0][na] + iX;    
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
	      for (n = cellListL[n]; n > -1; n = cellListL[n]) 
		{
		  if (n != na)// && n != nb && (nb >= -1 || n < na)) 
		    {
#ifdef EDHE_FLEX
		      if (!may_interact_all(na, n))
			continue;
#endif
		      //assign_bond_mapping(na, n);
		      //dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
#ifdef MD_NNLPLANES
		      dist = calcDistNegNNLoverlapPlane(Oparams.time, 0.0, na, n, shift); 
#else
		      dist = calcDistNegNNLoverlap(Oparams.time, 0.0, na, n, shift,
						   r1, r2, &alpha, vecgsup, 1); 
#endif
		      /* 0.1 è un buffer per evitare problemi, deve essere un parametro 
		       * in OprogStatus */
		      if (dist < 0)
			{
			  MD_DEBUG32(printf("Adding ellipsoid N. %d to NNL of %d\n", n, na));
 			  nebrTab[na].list[nebrTab[na].len] = n;
			  //for (kk=0; kk < 3; kk++)
			  //nebrTab[na].shift[nebrTab[na].len][kk] = shift[kk];
			  nebrTab[na].len++;
			  check_nnl_size(na);
			}
		    }
		} 
	    }
	}
    }

#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=1;
#endif
}
