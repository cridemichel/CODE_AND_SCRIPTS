#ifdef MC_HC
#undef DEBUG_HCMC
#undef MC_HC_SPHERO_OPT
#include<mdsimul.h>
#undef MD_MAC
#ifdef MD_MAC
#include <Accelerate/Accelerate.h>
#endif
#include<float.h>
#include<complex.h>
#ifndef CMPLX
#define CMPLX(x,y) (x)+I*(y)
#endif
#ifndef CMPLXL
#define CMPLXL(x,y) (x)+I*(y)
#endif

/* NOTA SU RISOLUZIONE QUARTICA
 * la routine gsl è un 15% più lenta della routine hqr() di Numerical Recipe.
 * La routine di Numerical Recipe sembra essere più accurata di quella delle gsl.
 * laguerre va come la routine gsl ma sembra molto inaccurate, mentre le lapack sono un po più lente.*/
//#define MC_QUART_LONG_DOUBLE
//#define POLY_SOLVE_GSL
//#define USE_LAPACK
//#define USE_LAGUERRE
//#define MC_DEBUG_HCALGO
//#define MC_EXCHG_QUART_SOL
//#define MC_QUART_VERBOSE
#define MC_QUART_HYBRID
#define FAST_QUARTIC_SOLVER
#define LDLT_OPT_CUBIC
//#define USE_CUBIC_HANDLE_INF_ONLY

/* NOTA 21/02/18: per quanto riguarda le simulazioni queste due define incidono pochissimo e danno più robustezza
 * volendo quindi possono essere attivate senza problemi */
#define LDLT_USENR /* newton-raphson per il calcolo di c e a in LDLT */
#define LDLT_REFINE_PHI_WITH_NR  /* newtown-raphson per il refinement di phi0 in LDLT */
//#define LDLT_LOSSOFSIG_FIX
//#define LDLT_USENRCMPLX
//#define MC_QUART_USE_ANALYTIC
//#include <gsl/gsl_poly.h>
//#include <gsl/gsl_errno.h>
void versor_to_R_alt(double *Ci, double *ni, double *Dj, double *nj, double R[3][3], double D);
double test_overlap_parall_cyll(long double *Ci, long double *ni, long double *Dj, long double *nj, long double Li, 
				long double Diami, long double Diamj);
double test_overlap_parall_cyl(double *Ci, double *ni, double *Dj, double *nj, double Li, 
	       		       double Diami, double Diamj);
void solve_numrecl(long double coeff[5], int *numrealsol, long double rsol[4], int *ok, int m);
void solve_numrec(double coeff[5], int *numrealsol, double rsol[4], int *ok, int m);

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
#define SIGNL(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#ifndef MAX
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
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
extern long double calc_norml(long double *vec);
extern double calc_normsq(double *vec);
extern void vectProdVec(double *A, double *B, double *C);
extern void vectProdVecl(long double *A, long double *B, long double *C);
extern void print_matrix(double **M, int n);
extern void update_MSDrot(int i);
extern void update_MSD(int i);
#ifdef MC_DEBUG_HCALGO
int iGbl, jGbl;
#endif
#ifdef MD_SUPERELLIPSOID
extern int is_superellipse(int i);
extern void fdjacSE(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3]);
#endif
#ifdef MD_ASYM_ITENS
extern void calc_omega(int i, double *wwx, double *wwy, double *wwz);
extern void calc_angmom(int i, double **);
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
extern double Vz;
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
long double scalProdl(long double *A, long double *B);

extern double toteneini;
extern long long int ttini;
extern int covrestart;
extern const int nmboxMC;
extern double totdist, distcc;
static double maxarg1,maxarg2;
/* ==================== FLOCKE ======================= */
#if 1
void quadraticRoots (double q1, double q0, int *nReal, complex double root[2]);

double Max2(double a, double b)
{
  if (a >= b)
    return a;
  else
    return b;
}
double Max3(double a, double b, double c)
{
  double t;
  t = Max2(a,b);
  return Max2(t,c);
}
double Max4(double a, double b, double c, double d)
{
  double t;
  t = Max3(a,b,c);
  return Max2(t,d);
}
 double Min2(double a, double b)
{
  if (a <= b)
    return a;
  else
    return b;
}
 double Min3(double a, double b, double c)
{
  double t;
  t = Min2(a,b);
  return Min2(t,c);
}
 double Min4(double a, double b, double c, double d)
{
  double t;
  t = Min3(a,b,c);
  return Min2(t,d);
}
//
// CUBIC POLYNOMIAL ROOT SOLVER
//
// SYNOPSIS
//
//  call cubicRoots (real,              intent (in)  :: c2,
//                   real,              intent (in)  :: c1,
//                   real,              intent (in)  :: c0,
//                   integer,           intent (out) :: nReal,
//                   real,              intent (out) :: root (1:3,1:2),
//                   logical, optional, intent (in)  :: printInfo)
//
// DESCRIPTION
//
//  Calculates all real + complex roots of the cubic polynomial:
//
//                 x^3 + c2 * x^2 + c1 * x + c0
//
//  The first real root (which always exists) is obtained using an optimized
//  Newton-Raphson scheme. The other remaining roots are obtained through
//  composite deflation into a quadratic. An option for printing detailed info
//  about the intermediate stages in solving the cubic is available.
//
//  The cubic root solver can handle any size of cubic coefficients and there is
//  no danger of overflow due to proper rescaling of the cubic polynomial.
//
//  The order of the roots is as follows:
//
//        1) For real roots, the order is according to their algebraic value
//           on the number scale (largest positive first, largest negative last).
//
//        2) Since there can be only one complex conjugate pair root, no order
//           is necessary.
//
//        3) All real roots preceede the complex ones.
//
// ARGUMENTS
//
//  c2         : coefficient of x^2 term
//  c1         : coefficient of x term
//  c0         : independent coefficient
//  nReal      : number of real roots found
//  root (n,1) : real part of n-th root
//  root (n,2) : imaginary part of n-th root
//  printInfo  : if given and true, detailed info will be printed about intermediate stages
//
// NOTES
//
//***

void cubicRoots (double c2, double c1, double c0, int *nReal, complex double root[3])
{
  int bisection, converged;
  int cubicType, deflateCase, oscillate;

  enum costanti {allzero = 0, linear = 1, quadratic = 2, general   = 3};
  const int Re = 1, Im = 2;

  double a0, a1, a2, a=0, b=0, c=0, k, s, t, u, x, y, z, xShift;

  double macheps =2.2204460492503131E-16, one27th = 1.0 / 27.0;
  double two27th = 2.0 / 27.0, third   = 1.0 /  3.0;
  double p1 = 1.09574,q1 = 3.23900E-1,r1=3.23900e-1;      //
  double s1 = 9.57439E-2;      //

  double p3 = 1.14413;         //
  double q3 = 2.75509E-1;      // Newton-Raphson coeffs for class 3
  double r3 = 4.45578E-1;      //
  double s3 = 2.59342E-2;      //

  double q4 = 7.71845E-1;      // Newton-Raphson coeffs for class 4
  double s4 = 2.28155E-1;      //

  double p51 = 8.78558E-1;     //
  double p52 = 1.92823E-1;     //
  double p53 = 1.19748;        //
  double p54 = 3.45219E-1;     //
  double q51 = 5.71888E-1;     //
  double q52 = 5.66324E-1;     //
  double q53 = 2.83772E-1;     // Newton-Raphson coeffs for class 5 and 6
  double q54 = 4.01231E-1;     //
  double r51 = 7.11154E-1;     //
  double r52 = 5.05734E-1;     //
  double r53 = 8.37476E-1;     //
  double r54 = 2.07216E-1;     //
  double s51 = 3.22313E-1;     //
  double s52 = 2.64881E-1;     //
  double s53 = 3.56228E-1;     //
  double s54 = 4.45532E-3;     //
  //
  //
  //     ...Start.
  //
  //
  //
  //
  //     ...Handle special cases.
  //
  //            1) all terms zero
  //            2) only quadratic term is nonzero -> linear equation.
  //            3) only independent term is zero -> quadratic equation.
  //
  //
  if (c0 == 0.0 && c1 == 0.0 && c2 == 0.0) {

    cubicType = allzero;

  } else if (c0 == 0.0 && c1 == 0.0) {

    k  = 1.0;
    a2 = c2;

    cubicType = linear;

  } else if (c0 == 0.0) {

    k  = 1.0;
    a2 = c2;
    a1 = c1;

    cubicType = quadratic;

  } else {
    //
    //
    //     ...The general case. Rescale cubic polynomial, such that largest fabsolute coefficient
    //        is (exactly!) equal to 1. Honor the presence of a special cubic case that might have
    //        been obtained during the rescaling process (due to underflow in the coefficients).
    //
    //
    x = fabs (c2);
    y = sqrt (fabs (c1));
    z = pow(fabs (c0), third);
    u = Max3(x,y,z);

    if (u == x) {

      k  = 1.0 / x;
      a2 = copysign (1.0 , c2);
      a1 = (c1 * k) * k;
      a0 = ((c0 * k) * k) * k;

    } else if (u == y) {

      k  = 1.0 / y;
      a2 = c2 * k;
      a1 = copysign (1.0 , c1);
      a0 = ((c0 * k) * k) * k;

    } else {

      k  = 1.0 / z;
      a2 = c2 * k;
      a1 = (c1 * k) * k;
      a0 = copysign (1.0 , c0);

    }

    k = 1.0 / k;

    if (a0 == 0.0 && a1 == 0.0 && a2 == 0.0) {
      cubicType = allzero;
    } else if (a0 == 0.0 && a1 == 0.0) {
      cubicType = linear;
    } else if (a0 == 0.0) {
      cubicType = quadratic;
    } else {
      cubicType = general;
    }

  }
  //
  //
  //     ...Select the case.
  //
  //        1) Only zero roots.
  //
  //
  switch (cubicType)
    {
    case allzero:

      *nReal = 3;

      root[0]=root[1]=root[2] = 0.0+I*0.0;
      //root (:,Im) = 0.0;
      break;
      //
      //
      //     ...2) The linear equation case -> additional 2 zeros.
      //
      //
    case linear:

      x = - a2 * k;

	*nReal = 3;

	root [0] = Max2 (0.0, x)+I*0.0;
	root [1] = 0.0+I*0.0;
	root [2] = Min2 (0.0, x)+I*0.0;
	//
	//
	//     ...3) The quadratic equation case -> additional 1 zero.
	//
	//
	break;
    case quadratic:

      quadraticRoots (a2, a1, nReal, root);

      if (*nReal == 2) {

	x = creal(root[0]) * k;         // real roots of quadratic are ordered x >= y
	y = creal(root[1]) * k;

	*nReal = 3;

	root[0] = Max2 (x, 0.0)+I*0.0;
	root[1] = Max2 (y, Min2 (x, 0.0))+I*0.0;
	root[2] = Min2 (y, 0.0)+I*0.0;
	//root (:,Im) = 0.0;

      } else {

	*nReal = 1;

	root [2] = root [1] * k;
	root [1] = root [0] * k;
	root [0] = 0.0+I*0.0;
	//root (3,Im) = root (2,Im) * k;
	//root (2,Im) = root (1,Im) * k;
	//root (1,Im) = 0.0;

      }
      break;
      //
      //
      //     ...3) The general cubic case. Set the best Newton-Raphson root estimates for the cubic.
      //           The easiest and most robust conditions are checked first. The most complicated
      //           ones are last and only done when fabsolutely necessary.
      //
      //
    case general:

      if (a0 == 1.0) {

	x = - p1 + q1 * a1 - a2 * (r1 - s1 * a1);

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a0 == - 1.0) {

	x = p1 - q1 * a1 - a2 * (r1 - s1 * a1);

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a1 == 1.0) {

	if (a0 > 0.0) {
	  x = a0 * (- q4 - s4 * a2);
	} else {
	  x = a0 * (- q4 + s4 * a2);
	}

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a1 == - 1.0) {

	y = - two27th;
	y = y * a2;
	y = y * a2 - third;
	y = y * a2;

	if (a0 < y) {
	  x = + p3 - q3 * a0 - a2 * (r3 + s3 * a0);               // + guess
	} else {
	  x = - p3 - q3 * a0 - a2 * (r3 - s3 * a0);               // - guess
	}

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a2 == 1.0) {

	b = a1 - third;
	c = a0 - one27th;

	if (fabs (b) < macheps && fabs (c) < macheps) {        // triple -1/3 root

	  x = - third * k;

	  *nReal = 3;

	  root[0]=root[1]=root[2] = x + I*0.0;
	  return;

	} else {

	  y = third * a1 - two27th;

	  if (a1 <= third) {
	    if (a0 > y) {
	      x = - p51 - q51 * a0 + a1 * (r51 - s51 * a0);   // - guess
	    } else {
	      x = + p52 - q52 * a0 - a1 * (r52 + s52 * a0);   // + guess
	    }
	  } else {
	    if (a0 > y) {
	      x = - p53 - q53 * a0 + a1 * (r53 - s53 * a0);   // <-1/3 guess
	    } else {
	      x = + p54 - q54 * a0 - a1 * (r54 + s54 * a0);   // >-1/3 guess
	    }
	  }

	  if (fabs (b) < 1.e-2 && fabs (c) < 1.e-2) {  // use shifted root
	    c = - third * b + c;
	    if (fabs (c) < macheps) c = 0.0;                  // prevent random noise
	    a = 0.0;
	    xShift = third;
	    x = x + xShift;
	  } else {
	    a = a2;
	    b = a1;
	    c = a0;
	    xShift = 0.0;
	  }

	}

      } else if (a2 == - 1.0) {

	b = a1 - third;
	c = a0 + one27th;

	if (fabs (b) < macheps && fabs (c) < macheps) {        // triple 1/3 root

	  x = third * k;

	  *nReal = 3;

	  root [0]=root[1]=root[2] = x + I*0.0;
	  //root (:,Im) = 0.0;

	  return;

	} else {

	  y = two27th - third * a1;

	  if (a1 <= third) {
	    if (a0 < y) {
	      x = + p51 - q51 * a0 - a1 * (r51 + s51 * a0);   // +1 guess
	    } else {
	      x = - p52 - q52 * a0 + a1 * (r52 - s52 * a0);   // -1 guess
	    }
	  } else {
	    if (a0 < y) {
	      x = + p53 - q53 * a0 - a1 * (r53 + s53 * a0);   // >1/3 guess
	    } else {
	      x = - p54 - q54 * a0 + a1 * (r54 - s54 * a0);   // <1/3 guess
	    }
	  }

	  if (fabs (b) < 1.e-2 && fabs (c) < 1.e-2) {  // use shifted root
	    c = third * b + c;
	    if (fabs (c) < macheps) c = 0.0;                  // prevent random noise
	    a = 0.0;
	    xShift = - third;
	    x = x + xShift;
	  } else {
	    a = a2;
	    b = a1;
	    c = a0;
	    xShift = 0.0;
	  }

	}

      }
      //
      //
      //     ...Perform Newton/Bisection iterations on x^3 + ax^2 + bx + c.
      //
      //
      z = x + a;
      y = x + z;
      z = z * x + b;
      y = y * x + z;       // C'(x)
      z = z * x + c;       // C(x)
      t = z;               // save C(x) for sign comparison
      x = x - z / y;       // 1st improved root

      oscillate = 0;
      bisection = 0;
      converged = 0;

      while (!converged && !bisection)    // Newton-Raphson iterates
	{
	  z = x + a;
	  y = x + z;
	  z = z * x + b;
	  y = y * x + z;
	  z = z * x + c;

	  if (z * t < 0.0) {                       // does Newton start oscillating ?
	    if (z < 0.0) {
	      oscillate = oscillate + 1;              // increment oscillation counter
	      s = x;                                  // save lower bisection bound
	    } else {
	      u = x;                                  // save upper bisection bound
	    }
	    t = z;                                      // save current C(x)
	  }

	  y = z / y;                                      // Newton correction
	  x = x - y;                                      // new Newton root

	  bisection = (oscillate > 2)?1:0;                      // activate bisection
	  converged = (fabs (y) <= fabs (x) * macheps)?1:0;       // Newton convergence indicator


	}

      if (bisection) {

	t = u - s;                                     // initial bisection interval
	while (fabs (t) > fabs (x) * macheps)        // bisection iterates
	  {
	    z = x + a;                                  //
	    z = z * x + b;                              // C (x)
	    z = z * x + c;                              //

	    if (z < 0.0) {                       //
	      s = x;                                  //
	    } else {                                       // keep bracket on root
	      u = x;                                  //
	    }                                     //

	    t = 0.5 * (u - s);                       // new bisection interval
	    x = s + t;                                  // new bisection root


	  }
      }


      x = x - xShift;                                   // unshift root
      //
      //
      //     ...Forward / backward deflate rescaled cubic (if needed) to check for other real roots.
      //        The deflation analysis is performed on the rescaled cubic. The actual deflation must
      //        be performed on the original cubic, not the rescaled one. Otherwise deflation errors
      //        will be enhanced when undoing the rescaling on the extra roots.
      //
      //
      z = fabs (x);
      s = fabs (a2);
      t = fabs (a1);
      u = fabs (a0);

      y = z * Max2 (s,z);           // take maximum between |x^2|,|a2 * x|

      deflateCase = 1;             // up to now, the maximum is |x^3| or |a2 * x^2|

      if (y < t) {             // check maximum between |x^2|,|a2 * x|,|a1|
	y = t * z;               // the maximum is |a1 * x|
	deflateCase = 2;         // up to now, the maximum is |a1 * x|
      } else {
	y = y * z;               // the maximum is |x^3| or |a2 * x^2|
      }

      if (y < u) {             // check maximum between |x^3|,|a2 * x^2|,|a1 * x|,|a0|
	deflateCase = 3;         // the maximum is |a0|
      }

      y = x * k;                   // real root of original cubic

      switch (deflateCase)
	{
	case 1:
	  x = 1.0 / y;
	  t = - c0 * x;              // t -> backward deflation on unscaled cubic
	  s = (t - c1) * x;          // s -> backward deflation on unscaled cubic
	  break;	
	case 2:
	  s = c2 + y;                // s ->  forward deflation on unscaled cubic
	  t = - c0 / y;              // t -> backward deflation on unscaled cubic
	  break;
	case 3:
	  s = c2 + y;                // s ->  forward deflation on unscaled cubic
	  t = c1 + s * y;            // t ->  forward deflation on unscaled cubic
	}

      quadraticRoots (s, t, nReal, root);

      if (*nReal == 2) {

	x = creal(root[0]);         // real roots of quadratic are ordered x >= z
	z = creal(root[1]);         // use 'z', because 'y' is original cubic real root


	*nReal = 3;

	root [0] = Max2 (x, y)+I*0.0;
	root [1] = Max2 (z, Min2 (x, y))+I*0.0;
	root [2] = Min2 (z, y)+I*0.0;
	//root (:,Im) = 0.0;

      } else {

	*nReal = 1;

	root [2] = root[1];
	root [1] = root[0];
	root [0] = y + I*0.0;
	//root (3,Im) = root (2,Im);
	//root (2,Im) = root (1,Im)
	  //root (1,Im) = 0.0;

      }
    }
  // end select
  //
  //
  //     ...Ready!
  //
  //
  return;
}
//end subroutine cubicRoots



//-----------------------------------------------------------------------------------
//
// QUADRATIC POLYNOMIAL ROOT SOLVER
//
// SYNOPSIS
//
//  call quadraticRoots (real,    intent (in)  :: q1,
//                       real,    intent (in)  :: q0,
//                       integer, intent (out) :: nReal,
//                       real,    intent (out) :: root (1:2,1:2))
//
// DESCRIPTION
//
//  Calculates all real + complex roots of the quadratic polynomial:
//
//                 x^2 + q1 * x + q0
//
//  The code checks internally, if rescaling of the coefficients is needed to
//  avoid overflow.
//
//  The order of the roots is as follows:
//
//        1) For real roots, the order is according to their algebraic value
//           on the number scale (largest positive first, largest negative last).
//
//        2) Since there can be only one complex conjugate pair root, no order
//           is necessary.
//
// ARGUMENTS
//
//  q1         : coefficient of x term
//  q0         : independent coefficient
//  nReal      : number of real roots found
//  root (n,1) : real part of n-th root
//  root (n,2) : imaginary part of n-th root
//
// NOTES
//
//***

void quadraticRoots (double q1, double q0, int *nReal, complex double root[2])
{

  int rescale;

  double a0, a1;
  double k, x, y, z;

  const int Re = 1, Im = 2;

  const double LPN = 1.7976931348623157E+308;   // the (L)argest (P)ositive (N)umber
  const double sqrtLPN = sqrt (LPN);      // and the square root of it
//
//
//     ...Handle special cases.
//
//
  if (q0 == 0.0 && q1 == 0.0) {

      *nReal = 2;

      root [0] = root[1] = 0.0+I*0.0;
      //root (:,Im) = 0.0;

  } else if (q0 == 0.0) {

      *nReal = 2;

      root [0] = Max2 (0.0, - q1)+I*0.0;
      root [1] = Min2 (0.0, - q1)+I*0.0;

  } else if (q1 == 0.0) {

      x = sqrt (fabs (q0));

      if (q0 < 0.0) {

          *nReal = 2;

          root [0] = x+I*0.0;
          root [1] = - x+I*0.0;
	  // root (:,Im) = 0.0;

      } else {

          *nReal = 0;

          root [0] = 0.0+I*x;
          root [1] = 0.0 - x*I;

      }

  } else {
//
//
//     ...The general case. Do rescaling, if either squaring of q1/2 or evaluation of
//        (q1/2)^2 - q0 will lead to overflow. This is better than to have the solver
//        crashed. Note, that rescaling might lead to loss of accuracy, so we only
//        invoke it when fabsolutely necessary.
//
//
      rescale = (q1 > (sqrtLPN + sqrtLPN))?1:0;     // this detects overflow of (q1/2)^2

      if (!rescale) {
           x = q1 * 0.5;                      // we are sure here that x*x will not overflow
           rescale = (q0 < (x * x - LPN));      // this detects overflow of (q1/2)^2 - q0
      }

      if (rescale) {

          x = fabs (q1);
          y = sqrt (fabs (q0));

          if (x > y) {
              k  = x;
              z  = 1.0 / x;
              a1 = copysign (1.0 , q1);
              a0 = (q0 * z) * z;
          } else {
              k  = y;
              a1 = q1 / y;
              a0 = copysign (1.0 , q0);
          }

      } else {
          a1 = q1;
          a0 = q0;
      }
//
//
//     ...Determine the roots of the quadratic. Note, that either a1 or a0 might
//        have become equal to zero due to underflow. But both cannot be zero.
//
//
      x = a1 * 0.5;
      y = x * x - a0;

      if (y >= 0.0) {

          y = sqrt (y);

          if (x > 0.0) {
              y = - x - y;
          } else {
              y = - x + y;
          }

          if (rescale) {
              y = y * k;                     // very important to convert to original
              z = q0 / y;                    // root first, otherwise complete loss of
          } else {                              // root due to possible a0 = 0 underflow
              z = a0 / y;
          }

          *nReal = 2;

          root [0] = Max2 (y,z)+I*0.0;           // 1st real root of x^2 + a1 * x + a0
          root [1] = Min2 (y,z)+I*0.0;           // 2nd real root of x^2 + a1 * x + a0

      } else {

          y = sqrt (- y);

          *nReal = 0;

          root [0] = - x+I*y;
          root [1] = - x - y*I;

          if (rescale) {
              root[0] = root[0] * k;
	      root[1] = root[1] * k;
          }

      }

  }
//
//
//     ...Ready!
//
//
  return;
}




//-----------------------------------------------------------------------------------
//
// QUARTIC POLYNOMIAL ROOT SOLVER
//
// SYNOPSIS
//
//  call quarticRoots (real,              intent (in)  :: q3,
//                     real,              intent (in)  :: q2,
//                     real,              intent (in)  :: q1,
//                     real,              intent (in)  :: q0,
//                     integer,           intent (out) :: nReal,
//                     real,              intent (out) :: root (1:4,1:2),
//                     logical, optional, intent (in)  :: printInfo)
//
// DESCRIPTION
//
//  Calculates all real + complex roots of the quartic polynomial:
//
//                 x^4 + q3 * x^3 + q2 * x^2 + q1 * x + q0
//
//  An option for printing detailed info about the intermediate stages in solving
//  the quartic is available. This enables a detailed check in case something went
//  wrong and the roots obtained are not proper.
//
//  The quartic root solver can handle any size of quartic coefficients and there is
//  no danger of overflow, due to proper rescaling of the quartic polynomial.
//
//  The order of the roots is as follows:
//
//        1) For real roots, the order is according to their algebraic value
//           on the number scale (largest positive first, largest negative last).
//
//        2) For complex conjugate pair roots, the order is according to the
//           algebraic value of their real parts (largest positive first). If
//           the real parts are equal, the order is according to the algebraic
//           value of their imaginary parts (largest first).
//
//        3) All real roots preceede the complex ones.
//
// ARGUMENTS
//
//  q3         : coefficient of x^3 term
//  q2         : coefficient of x^2 term
//  q1         : coefficient of x term
//  q0         : independent coefficient
//  nReal      : number of real roots found
//  root (n,1) : real part of n-th root
//  root (n,2) : imaginary part of n-th root
//  printInfo  : if given and true, detailed info will be printed about intermediate stages
//
// NOTES
//
//***
// VERSIONE AGGIORNATA AL 2018 più veloce e robusta

//end subroutine quarticRoots
void quarticRoots (double cc[5], int *nReal, complex double root[4])
{
  int bisection;
  int converged;
  //int doPrint;
  int iterate;
  int minimum;
  int overshoot;
  int dsignflip;
  int notZero;

  int deflateCase;
  //int oscillate;
  int quarticType;

  const int Re = 1, Im = 2;
  double q1, q2, q3, q0;
  //const int biquadratic = 2, cubic = 3, general = 4; 
  enum casi {biquadratic=2, cubic=3, general=4 };
  double a0, a1, a2, a3;
  double a, b, c, d, k, s, t, u, x, y, z;

  const double macheps = 2.2204460492503131E-16; 
  const double third   = 1.0 / 3.0;
  //
  //
  //     ...Start.
  //
  //
#if 0
  if (present (printInfo)) {
    doPrint = printInfo
  } else {
    doPrint = 0
  }
#endif

#if 0
  if (doPrint) {
    write (*,wpformat) ' initial quartic q3    = ',q3
      write (*,wpformat) ' initial quartic q2    = ',q2
      write (*,wpformat) ' initial quartic q1    = ',q1
      write (*,wpformat) ' initial quartic q0    = ',q0
      write (*,wpformat) ' ------------------------------------------------'
  }
#endif
  //
  //
  //     ...Handle special cases. Since the cubic solver handles all its
  //        special cases by itself, we need to check only for two cases:
  //
  //            1) independent term is zero -> solve cubic and include
  //               the zero root
  //
  //            2) the biquadratic case.
  //
  //
  q3=cc[3]/cc[4];
  q2=cc[2]/cc[4];
  q1=cc[1]/cc[4];
  q0=cc[0]/cc[4];
  if (q0 == 0.0) {

    k  = 1.0;
    a3 = q3;
    a2 = q2;
    a1 = q1;

    quarticType = cubic;

  } else if (q3 == 0.0 && q1 == 0.0) {

    k  = 1.0;
    a2 = q2;
    a0 = q0;

    quarticType = biquadratic;

  } else {
    //
    //
    //     ...The general case. Rescale quartic polynomial, such that largest fabsolute coefficient
    //        is (exactly!) equal to 1. Honor the presence of a special quartic case that might have
    //        been obtained during the rescaling process (due to underflow in the coefficients).
    //
    //
    s = fabs (q3);
    t = sqrt (fabs (q2));
    u = pow(fabs (q1), third);
    x = sqrt (sqrt (fabs (q0)));
    y = Max4 (s,t,u,x);

    if (y == s) {

      k  = 1.0 / s;
      a3 = copysign (1.0 , q3);
      a2 = (q2 * k) * k;
      a1 = ((q1 * k) * k) * k;
      a0 = (((q0 * k) * k) * k) * k;

    } else if (y == t) {

      k  = 1.0 / t;
      a3 = q3 * k;
      a2 = copysign (1.0 , q2);
      a1 = ((q1 * k) * k) * k;
      a0 = (((q0 * k) * k) * k) * k;

    } else if (y == u) {

      k  = 1.0 / u;
      a3 = q3 * k;
      a2 = (q2 * k) * k;
      a1 = copysign (1.0 , q1);
      a0 = (((q0 * k) * k) * k) * k;

    } else {

      k  = 1.0 / x;
      a3 = q3 * k;
      a2 = (q2 * k) * k;
      a1 = ((q1 * k) * k) * k;
      a0 = copysign (1.0 , q0);

    }

    k = 1.0 / k;

#if 0
    if (doPrint) {
      write (*,wpformat) ' rescaling factor      = ',k
	write (*,wpformat) ' ------------------------------------------------'
	write (*,wpformat) ' rescaled quartic q3   = ',a3
	write (*,wpformat) ' rescaled quartic q2   = ',a2
	write (*,wpformat) ' rescaled quartic q1   = ',a1
	write (*,wpformat) ' rescaled quartic q0   = ',a0
	write (*,wpformat) ' ------------------------------------------------'
    }
#endif
    if (a0 == 0.0) {
      quarticType = cubic;
    } else if (a3 == 0.0 && a1 == 0.0) {
      quarticType = biquadratic;
    } else {
      quarticType = general;
    }

  }
  //
  //
  //     ...Select the case.
  //
  //        1) The quartic with independent term = 0 -> solve cubic and add a zero root.
  //
  //
  switch (quarticType)
    {
    case cubic:

      cubicRoots (a3, a2, a1, nReal, root);

      if (*nReal == 3) {

	x = creal(root[0]) * k;       // real roots of cubic are ordered x >= y >= z
	y = creal(root[1])* k;
	z = creal(root[2]) * k;

	*nReal = 4;

	root [0] = Max2 (x, 0.0) + I*0.0;
	root [1] = Max2 (y, Min2 (x, 0.0))+I*0.0;
	root [2] = Max2 (z, Min2 (y, 0.0))+I*0.0;
	root [3] = Min2 (z, 0.0)+I*0.0;
	//root (:,Im) = 0.0;

      } else {                          // there is only one real cubic root here

	x = creal(root[0]) * k;

	*nReal = 2;

	root [3] = root [2] * k;
	root [2] = root [1] * k;
	root [1] = Min2 (x, 0.0)+I*0.0;
	root [0] = Max2 (x, 0.0)+I*0.0;

      }
      break;
      //
      //
      //     ...2) The quartic with x^3 and x terms = 0 -> solve biquadratic.
      //
      //
      //
    case biquadratic:

      quadraticRoots (q2, q0, nReal, root);

	if (*nReal == 2) {

	  x = creal(root [0]);         // real roots of quadratic are ordered x >= y
	  y = creal(root [1]);

	    if (y >= 0.0) {

	      x = sqrt (x) * k;
      	      y = sqrt (y) * k;

      	      *nReal = 4;

      	      root [0] = x+I*0.0;
	      root [1] = y+I*0.0;
	      root [2] = - y+I*0.0;
	      root [3] = - x+I*0.0;

	    } else if (x >= 0.0 && y < 0.0) {

	      x = sqrt (x)       * k;
      	      y = sqrt (fabs (y)) * k;

      	      *nReal = 2;

      	      root [0] = x+I*0.0;
	      root [1] = - x+I*0.0;
	      root [2] = 0.0+I*y;
	      root [3] = 0.0-I*y;

	    } else if (x < 0.0) {

	      x = sqrt (fabs (x)) * k;
      	      y = sqrt (fabs (y)) * k;

      	      *nReal = 0;

      	      root [0] = 0+I*y;
	      root [1] = 0+I*x;
	      root [2] = 0- I*x;
	      root [3] = 0- I*y;

	    }

	} else {          // complex conjugate pair biquadratic roots x +/- iy.

	  x = creal(root [0]) * 0.5;
	  y = cimag(root[0]) * 0.5;
	  z = sqrt (x * x + y * y);
	  y = sqrt (z - x) * k;
	  x = sqrt (z + x) * k;

	  *nReal = 0;

	  root [0] = x+I*y;
	  root [1] = x-I*y;
	  root [2] = - x+I*y;
	  root [3] = - x-I*y;

	}
      break;
      //
      //
      //     ...3) The general quartic case. Search for stationary points. Set the first
      //           derivative polynomial (cubic) equal to zero and find its roots.
      //           Check, if any minimum point of Q(x) is below zero, in which case we
      //           must have real roots for Q(x). Hunt down only the real root, which
      //           will potentially converge fastest during Newton iterates. The remaining
      //           roots will be determined by deflation Q(x) -> cubic.
      //
      //           The best roots for the Newton iterations are the two on the opposite
      //           ends, i.e. those closest to the +2 and -2. Which of these two roots
      //           to take, depends on the location of the Q(x) minima x = s and x = u,
      //           with s > u. There are three cases:
      //
      //              1) both Q(s) and Q(u) < 0
      //                 ----------------------
      //
      //                 The best root is the one that corresponds to the lowest of
      //                 these minima. If Q(s) is lowest -> start Newton from +2
      //                 downwards (or zero, if s < 0 and a0 > 0). If Q(u) is lowest
      //                 -> start Newton from -2 upwards (or zero, if u > 0 and a0 > 0).
      //
      //              2) only Q(s) < 0
      //                 -------------
      //
      //                 With both sides +2 and -2 possible as a Newton starting point,
      //                 we have to avoid the area in the Q(x) graph, where inflection
      //                 points are present. Solving Q''(x) = 0, leads to solutions
      //                 x = -a3/4 +/- discriminant, i.e. they are centered around -a3/4.
      //                 Since both inflection points must be either on the r.h.s or l.h.s.
      //                 from x = s, a simple test where s is in relation to -a3/4 allows
      //                 us to avoid the inflection point area.
      //
      //              3) only Q(u) < 0
      //                 -------------
      //
      //                 Same of what has been said under 2) but with x = u.
      //
      //
    case general:

      x = 0.75 * a3;
      y = 0.50 * a2;
      z = 0.25 * a1;

      cubicRoots (x, y, z, nReal, root);

      s = creal(root[0]);        // Q'(x) root s (real for sure)
      x = s + a3;
      x = x * s + a2;
      x = x * s + a1;
      x = x * s + a0;         // Q(s)

      y = 1.0;             // dual info: Q'(x) has more real roots, and if so, is Q(u) < 0 ? 

      if (*nReal > 1) {
	u = creal(root [2]);   // Q'(x) root u
	y = u + a3;
	y = y * u + a2;
	y = y * u + a1;
	y = y * u + a0;     // Q(u)
      }

#if 0
      if (doPrint) {
	write (*,wpformat) ' dQ(x)/dx root s       = ',s
	  write (*,wpformat) ' Q(s)                  = ',x
	  write (*,wpformat) ' dQ(x)/dx root u       = ',u
	  write (*,wpformat) ' Q(u)                  = ',y
	  write (*,wpformat) ' ------------------------------------------------'
      }
#endif
      if (x == 0.0 && y == 0.0) 
	{
	  if (fabs(s) > fabs(u))
	    x=s;
	  else
	    x=u;
	  *nReal = 1;
	  iterate = 0;
	} 
      else if (x==0)
	{
	  x=s;
	  *nReal = 1;
	  iterate = 0;
	}
      else if (y==0)
	{
	  x = u;
	  *nReal = 1;
	  iterate = 0;
	}
      else if (x < 0.0 && y < 0.0)
	{
	  if (s < 0.0 && a0 > 0.0) 
	    x = 0.0;
          else
	    x = 2.0;
          
          if (u > 0.0 && a0 > 0.0) 
              y = 0.0;
          else
              y = -2.0;

          a = x + a3;
          b = x + a;
          a = a * x + a2;
          b = b * x + a;
          a = a * x + a1;
          b = b * x + a;     // b = Q'(x)

          c = y + a3;
          d = y + c;
          c = c * y + a2;
          d = d * y + c;
          c = c * y + a1;
          d = d * y + c;     //! d = Q'(y)

	  if (fabs (b) > fabs (d))   
	    {                           // if Q'(y) < Q'(x),
              x = y;                     // take root u for Newton iterations
              s = u;                     // save for lower bisecion bound just in case
	    }

          *nReal = 1;
          iterate = 1;
	}
      else if (x < 0.0) {

	if (s < - a3 * 0.25) {
	  if (s > 0.0 && a0 > 0.0) {
	    x = 0.0;
	  } else {
	    x = - 2.0;
	  }
	} else {
	  if (s < 0.0 && a0 > 0.0) {
	    x = 0.0;
	  } else {
	    x = 2.0;
	  }
	}

	*nReal = 1;
	iterate = 1;

      } else if (y < 0.0) {

	if (u < - a3 * 0.25) {
	  if (u > 0.0 && a0 > 0.0) {
	    x =  0.0;
	  } else {
	    x = -2.0;
	  }
	} else {
	  if (u < 0.0 && a0 > 0.0) {
	    x = 0.0;
	  } else {
	    x = 2.0;
	  }
	}
	s = u;
	*nReal = 1;
	iterate = 1;
      } else {
	*nReal = 0;
      }
      //
      //
      //     ...Do all necessary Newton iterations. In case we have more than 2 oscillations,
      //        exit the Newton iterations and switch to bisection. Note, that from the
      //        definition of the Newton starting point, we always have Q(x) > 0 and Q'(x)
      //        starts (-ve/+ve) for the (-2/+2) starting points and (increase/decrease) smoothly
      //        and staying (< 0 / > 0). In practice, for extremely shallow Q(x) curves near the
      //        root, the Newton procedure can overshoot slightly due to rounding errors when
      //        approaching the root. The result are tiny oscillations around the root. If such
      //        a situation happens, the Newton iterations are abandoned after 3 oscillations
      //        and further location of the root is done using bisection starting with the
      //        oscillation brackets.
      //
      //
      if (*nReal > 0) {
	if (iterate)
	  {
	    y = x + a3;                                     //
	    z = x + y;                                     //
	    y = y * x + a2;                                 // y = Q(x)
	    z = z * x + y;                                  //
	    y = y * x + a1;                                 // z = Q'(x)
	    z = z * x + y;                                  //
	    y = y * x + a0;                                 //
	    t = z;
	    y = y / z;
	    x = x - y;

	    dsignflip = 0;
	    overshoot = 0;
	    bisection = 0;
	    converged = 0;

	    while (!converged && !bisection)    // Newton-Raphson iterates
	      {
		y = x + a3;                                      //
		z = x + y;                                      //
		y = y * x + a2;                                 // y = Q(x)
		z = z * x + y;                                  //
		y = y * x + a1;                                 // z = Q'(x)
		z = z * x + y;                                  //
		y = y * x + a0;                                 //

		if (y < 0.0)    
		  {                                            // does Newton start overshooting ?
		    overshoot = overshoot + 1;                 // increment overshoot counter
		    s = x;
		  }                                            // save lower bisection bound
		else
		  u = x;                                      // save upper bisection bound

		if (z * t < 0.0)  
		  {  // does Q'(x) have a sign flip ?
		    dsignflip = dsignflip + 1;                  // increment sign flip counter
		    t = z;
		  }  // save Q'(x) for next sign check

		if (z == 0.0)                           // safeguard against accidental
		  {
		    bisection = 1;                        // Q'(x) = 0 due to roundoff
		    break;
		  } // errors -> activate bisection
		// with current bracket [s,u]

		y = y / z;                                      // Newton correction
		x = x - y;                                      // new Newton root

		bisection = (overshoot > 2 || dsignflip > 2)?1:0;   // activate bisection
		converged = (fabs (y) <= fabs (x) * macheps)?1:0;   // Newton convergence indicator
	      }

	    if (bisection) 
	      {
		t = u - s;                                     // initial bisection interval
		while (fabs (t) > fabs (x) * macheps)        // bisection iterates
		  {
		    y = x + a3;                                 //
		    y = y * x + a2;                             // y = Q(x)
		    y = y * x + a1;                             //
		    y = y * x + a0;                             //

		    if (y < 0.0)                       //
		      s = x;                                  //
		    else                                       // keep bracket on root
		      u = x;                                  //

		    t = 0.5 * (u - s);                       // new bisection interval
		    x = s + t;                                  // new bisection root
		  }
	      }
	  }

	//
	//
	//     ...Find remaining roots -> reduce to cubic. The reduction to a cubic polynomial
	//        is done using composite deflation to minimize rounding errors. Also, while
	//        the composite deflation analysis is done on the reduced quartic, the actual
	//        deflation is being performed on the original quartic again to avoid enhanced
	//        propagation of root errors.
	//
	//
	z = fabs (x);            //
	a = fabs (a3);           //
	b = fabs (a2);           // prepare for composite deflation
	c = fabs (a1);           //
	d = fabs (a0);           //

	y = z * Max2 (a,z);      // take maximum between |x^2|,|a3 * x|

	deflateCase = 1;        // up to now, the maximum is |x^4| or |a3 * x^3|

	if (y < b) {        // check maximum between |x^2|,|a3 * x|,|a2|
	  y = b * z;          // the maximum is |a2| -> form |a2 * x|
	  deflateCase = 2;    // up to now, the maximum is |a2 * x^2|
	} else {
	  y = y * z;          // the maximum is |x^3| or |a3 * x^2|
	}

	if (y < c) {        // check maximum between |x^3|,|a3 * x^2|,|a2 * x|,|a1|
	  y = c * z;          // the maximum is |a1| -> form |a1 * x|
	  deflateCase = 3;    // up to now, the maximum is |a1 * x|
	} else {
	  y = y * z;          // the maximum is |x^4|,|a3 * x^3| or |a2 * x^2|
	}

	if (y < d) {        // check maximum between |x^4|,|a3 * x^3|,|a2 * x^2|,|a1 * x|,|a0|
	  deflateCase = 4;    // the maximum is |a0|
	}

	x = x * k;              // 1st real root of original Q(x)

	switch (deflateCase)
	  {
	  case 1:
	    z = 1.0 / x;
	    u = - q0 * z;         // u -> backward deflation on original Q(x)
	    t = (u - q1) * z;     // t -> backward deflation on original Q(x)
	    s = (t - q2) * z;     // s -> backward deflation on original Q(x)
	    break;
	  case 2:
	    z = 1.0 / x;
	    u = - q0 * z;         // u -> backward deflation on original Q(x)
	    t = (u - q1) * z;     // t -> backward deflation on original Q(x)
	    s = q3 + x;           // s ->  forward deflation on original Q(x)
	    break;
	  case 3:
	    s = q3 + x;           // s ->  forward deflation on original Q(x)
	    t = q2 + s * x;       // t ->  forward deflation on original Q(x)
	    u = - q0 / x;         // u -> backward deflation on original Q(x)
	    break;
	  case 4:
	    s = q3 + x;           // s ->  forward deflation on original Q(x)
	    t = q2 + s * x;       // t ->  forward deflation on original Q(x)
	    u = q1 + t * x;       // u ->  forward deflation on original Q(x)
	  }

	cubicRoots (s, t, u, nReal, root);

	  if (*nReal == 3) {

	    s = creal(root [0]);    //
	    t = creal(root [1]);    // real roots of cubic are ordered s >= t >= u
	    u = creal(root [2]);    //

	    root [0] = Max2 (s, x)+I*0.0;
	    root [1] = Max2 (t, Min2 (s, x))+I*0.0;
	    root [2] = Max2 (u, Min2 (t, x))+I*0.0;
	    root [3] = Min2 (u, x)+I*0.0;

	    *nReal = 4;

	  } else {                   // there is only one real cubic root here

	    s = creal(root [0]);

	    root [3] = root [2];
	    root [2] = root [1];
	    root [1] = Min2 (s, x)+I*0.0;
	    root [0] = Max2 (s, x)+I*0.0;

	    *nReal = 2;

	  }

      } else {
	//
	//
	//     ...If no real roots have been found by now, only complex roots are possible.
	//        Find real parts of roots first, followed by imaginary components.
	//
	//
	s = a3 * 0.5;
	t =  s * s - a2;
	u =  s * t + a1;                   // value of Q'(-a3/4) at stationary point -a3/4


	notZero = (fabs (u) >= macheps)?1:0;    // H(-a3/4) is considered > 0 at stationary point

	if (a3 != 0.0) {
	  s = a1 / a3;
	  minimum = (a0 > (s * s))?1:0;                            // H''(-a3/4) > 0 -> minimum
	} else {
	  minimum = ((4 * a0) > (a2 * a2))?1:0;                      // H''(-a3/4) > 0 -> minimum
	}

	iterate = (notZero || (!notZero && minimum))?1:0;

	if (iterate) {

	  x = copysign (2.0,a3);                              // initial root -> target = smaller mag root

	  overshoot = 0;
	  bisection = 0;
	  converged = 0;

	  while (!converged && !bisection)    // Newton-Raphson iterates
	    {
	      a = x + a3;                                     //
	      b = x + a;                                      // a = Q(x)
	      c = x + b;                                      //
	      d = x + c;                                      // b = Q'(x)
	      a = a * x + a2;                                 //
	      b = b * x + a;                                  // c = Q''(x) / 2
	      c = c * x + b;                                  //
	      a = a * x + a1;                                 // d = Q'''(x) / 6
	      b = b * x + a;                                  //
	      a = a * x + a0;                                 //
	      y = a * d * d - b * c * d + b * b;              // y = H(x), usually < 0
	      z = 2 * d * (4 * a - b * d - c * c);            // z = H'(x)
	     
	      if (y > 0.0)                            // does Newton start oscillating ?
		{
		  overshoot = overshoot + 1;                  // increment oscillation counter
		  s = x;
		}
	      else
		{
		  u = x;	  
		}
	      if (z == 0.0)
		{
		  converged = 1;
		  bisection = 0;
		  break;
		}
	      y = y / z;                                      // Newton correction
	      x = x - y;                                      // new Newton root

	      bisection = (overshoot > 2)?1:0;                      // activate bisection
	      converged = (fabs (y) <= fabs (x) * macheps)?1:0;       // Newton convergence criterion

	      // if (doPrint) write (*,wpformat) ' Newton H(x) root      = ',x

	    }

	  if (bisection) {

	    t = u - s;                                     // initial bisection interval
	    while (fabs (t) > fabs (x) * macheps)        // bisection iterates
	      {
		a = x + a3;                                 //
		b = x + a;                                  // a = Q(x)
		c = x + b;                                  //
		d = x + c;                                  // b = Q'(x)
		a = a * x + a2;                             //
		b = b * x + a;                              // c = Q''(x) / 2
		c = c * x + b;                              //
		a = a * x + a1;                             // d = Q'''(x) / 6
		b = b * x + a;                              //
		a = a * x + a0;                             //
		y = a * d * d - b * c * d + b * b;          // y = H(x)

		if (y > 0.0) {                       //
		  s = x;                                  //
		} else {                                       // keep bracket on root
		  u = x;                                  //
		}                                     //

		t = 0.5 * (u - s);                       // new bisection interval
		x = s + t;                                  // new bisection root

	      }

	  }


	  a = x * k;                                         // 1st real component -> a
	  b = - 0.5 * q3 - a;                             // 2nd real component -> b
	  c = a * a;
	  d = b * b;

	  x = 4 * a + q3;                                    // Q'''(a)
	  y = x + q3 + q3;                                   //
	  y = y * a + q2 + q2;                               // Q'(a)
	  y = y * a + q1;                                    //
	  y = y / x;                                         // Q'(a) / Q'''(a)
	  s = c + y;                                         // magnitude^2 of (a + iy) root
	  x = 4 * b + q3;                                    // Q'''(b)
	  z = x + q3 + q3;                                   //
	  z = z * b + q2 + q2;                               // Q'(b)
	  z = z * b + q1;                                    //
	  z = z / x;                                         // Q'(b) / Q'''(b)
	  t = d + z;                                         // magnitude^2 of (b + iz) root

	  if (s > t) 
	    {                                   // minimize imaginary error
	      y = Max2(y, 0.0);                           // ensure >= 0 for sqrt
      	      d = Max2(q0 / s - d, 0.0);                  // small component using Vieta
	      c = sqrt(y);                                  // 1st imaginary component -> c
	      d = sqrt(d);
	    }// 2nd imaginary component -> d
	  else
	    {
	      c = Max2(q0 / t - c, 0.0);                  // small component using Vieta
   	      z = Max2(z, 0.0);                           // ensure >= 0 for sqrt
	      c = sqrt(c);                                  // 1st imaginary component -> c
	      d = sqrt(z);                                  // 2nd imaginary component -> d
	    }

	} else {                                                  // no bisection -> real components equal

	  a = - 0.25 * q3;                                // 1st real component -> a
	  b = a;                                             // 2nd real component -> b = a

	  x = a + q3;                                        //
	  x = x * a + q2;                                    // Q(a)
	  x = x * a + q1;                                    //
	  x = x * a + q0;                                    //
	  y = - 0.1875 * q3 * q3 + 0.5 * q2;           // Q''(a) / 2
	  z = Max2 (y * y - x, 0.0);                       // force discriminant to be >= 0
	  z = sqrt (z);                                      // square root of discriminant
	  y = y + copysign (z,y);                                // larger magnitude root
	  if (y == 0.0)                                    // guard against larger magnitude root = 0
	    x = 0.0;                                       // in this case smaller magnitude root must be 0
	  else
	    x = x / y;                                     // smaller magnitude root
	  c = Max2 (y, 0.0);                               // ensure root of biquadratic > 0
	  d = Max2 (x, 0.0);                               // ensure root of biquadratic > 0
	  c = sqrt (c);                                      // large magnitude imaginary component
	  d = sqrt (d);                                      // small magnitude imaginary component

	}

	if (a > b) {

	  root [0] = a+I*c;
	  root [1] = a-I*c;
	  root [2] = b+I*d;
	  root [3] = b-I*d;

	} else if (a < b) {

	  root [0] = b+I*d;
	  root [1] = b-I*d;
	  root [2] = a+I*c;
	  root [3] = a-I*c;

	} else {
	  root [0] = a+I*c;
	  root [1] = a-I*c;
	  root [2] = a+I*d;
	  root [3] = a-I*d;

	}

      }    // # of real roots 'if'
    }
  //
  //
  //     ...Ready!
  //
  //
  return;
}
//end subroutine quarticRoots
//

#endif




/* ====================================================== */
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
#ifdef HC_ALGO_OPT
double numcallsRR=0.0, numcallsRD=0.0, numcallsDD=0.0;
#endif
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
double calc_distance(double *A, double *B)
{
  int k;
  double d[3];
  for (k=0; k < 3; k++)
    {
      d[k] = A[k] - B[k];
    }
  return calc_norm(d);
}
#ifdef HC_ALGO_OPT
#define MESH_PTS 8 
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
  double m, q, dist, maxdist;
  double normCipP, CipP[3], ab2, asq1, nipP[3], sqrtdelta;
  double Pp[3];
  //static double *tharr;
  //static struct brentOpt *mesh;
  double sp, D2, lambda, delta, norm;
  int k1, k2, numsol=0;
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
  int k1, iter;
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
  return 0.0;
}
double find_initial_guess_bracket(double *thg, int meshpts, struct brentOpt* mesh, int *nnini, struct maxminS *maxmin)
{
  //static int firstcall=1;
  double th, dth, dist, maxdist;
  //static double *tharr;
  //static struct brentOpt *mesh;
  int k1, k2;
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
  double Pj[3], Rj[3][3];
  int kk, nn;
  static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3];
  static double *tharr;
  static double **mesh; /* {{1,0},{0.707106781186547, 0.707106781186547},{0,1},
      {-0.707106781186547,0.707106781186547},{-1,0},{-0.707106781186547,-0.707106781186547},
      {0,-1},{0.707106781186547,-0.707106781186547}};*/
  double PjCini, PjCi[3], d, mindist=-1.0; 
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
void find_initial_guess_simpler(double *Ai, double Ci[3], double ni[3], double Dj[3], double nj[3], double D)
{
  int kk1, kk2;
  double sp, dsc[3], norm;
  
  for (kk2=0; kk2 < 3; kk2++)
    dsc[kk2] = Ci[kk2] - Dj[kk2]; 
  sp = scalProd(dsc, ni);
  for (kk1=0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] - sp*ni[kk1];
} 
void find_initial_guess(double *Ai, double Ci[3], double ni[3], double Dj[3], double nj[3], double D)
{
  const int meshpts = 8;
  double Pj[3], Rj[3][3];
  int kk, nn;
  static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3];
  static double **mesh; /* {{1,0},{0.707106781186547, 0.707106781186547},{0,1},
      {-0.707106781186547,0.707106781186547},{-1,0},{-0.707106781186547,-0.707106781186547},
      {0,-1},{0.707106781186547,-0.707106781186547}};*/
  double PjCini, PjCi[3], d, mindist=-1.0; 
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
  double sp, Q1, Q2, normPiDi, normPjDj, normN, DiN, DjN, niN[3], njN[3], Djnj;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DjCi[3], normDjCi;
  double PiDi[3], PjDj[3], Ai[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3];
  double Dini;
  double DjCini;
  double Tjm_perp[3], Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp;
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
#if 0
	  find_initial_guess(Ai, Ci, ni, Dj[j2], nj, Diamj);
#else
	  find_initial_guess_simpler(Ai, Ci, ni, Dj[j2], nj, Diamj);
#if 0
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
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
  int k1;
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
  double fact, D2;
  int k2;
  D2 = 0.5*Dgbl;

  if (!(th==brentmsg.th) || brentmsg.id==0)
    {
      calcdrimdisk(th);
      brentmsg.id = 2;
    }
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

  if (!(th == brentmsg.th) || brentmsg.id==0)
    {
      calcrimdisk(th);
      brentmsg.id = 1;
    }
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
  double u[3], sp, norm;
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
//double calcDistNegHCdiffbrent(int i, int j, double shift[3], int* retchk);
extern double newton1D(double ax, double (*f)(double), double (*df)(double), double (*ddf)(double), double tol, double *xmin);
void solve_quadraticl(long double coeff[3], int *numsol, long double *sol)
{
  long double delta, a2inv, sqrtd;
  delta = Sqr(coeff[1]) - 4.0*coeff[2]*coeff[0];
  if (delta > 0.0)
    {
      sqrtd = sqrtl(delta);
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
void solve_quadratic_cmplx(double coeff[3], complex double *sol)
{
  /* numeric error safe version of solve_quadratic from Numerical Recipe */
  double delta, a, b, c, q;
  complex double cq;
  a = coeff[2];
  b = coeff[1];
  c = coeff[0];
  delta = Sqr(b) - 4.0*a*c;
  if (delta > 0.0)
    {
      q = -0.5*(b+copysign(1.0,b)*sqrt(delta));
      sol[0] = q/a;
      sol[1] = c/q;
    } 
  else if (delta == 0)
    {
      sol[0] = sol[1] = -b/(2.0*a);
    }
  else
    {
      cq = -0.5*(b+copysign(1.0,b)*csqrt(delta));
      sol[0] = cq/a;
      sol[1] = c/cq; 
    }
}

void solve_quadratic(double coeff[3], int *numsol, double *sol);
void csolve_cubic(double *coeff, double complex sol[3])
{
  const double sqrt3=sqrt(3.0);
  /* solution from Abramovitz */
  double q, r, q3, r2, a2, a1, a0, a2sq, H;
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
  cs1 = cpow(r + csqrt(H),1.0/3.0);
  cs2 = cpow(r - csqrt(H),1.0/3.0);
  sol[0] = cs1 + cs2 - a2/3.0;
  /* per il calcolo degli zeri di una quartica basta uno zero reale qualsiasi quindi 
   * con questo flag si puo' ottimizzare calcolandone solo uno */
  s1ps2 = (cs1+cs2)/2.0; 
  s1ms2 = (cs1-cs2)/2.0;
  sol[1] = -s1ps2 - a2/3.0 + I*sqrt3*s1ms2;
  sol[2] = -s1ps2 - a2/3.0 - I*sqrt3*s1ms2;
}
#if 1
void solve_cubicl(long double *coeff, int *numsol, long double sol[3])
{
  int ok;
  const long double sqrt3=sqrtl(3.0);
  /* solution from Abramovitz */
  long double s1, s2, q, r, q3, r2, a2, a1, a0, a2sq, H;
  long double complex cs1, cs2, s1ps2, s1ms2;
  if (coeff[3] == 0.0)
    {
      printf("[WARNING] fallback to quadratic from cubic\n");
      solve_numrecl(coeff, numsol, sol, &ok, 2);
      // solve_quadraticl(coeff, numsol, sol);
      return;
    }
  solve_numrecl(coeff, numsol, sol, &ok, 3);
  return;
  //printf("qu?????????????????????????????\n");
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
      cs1 = cpowl(r + csqrtl(H),1.0/3.0);
      cs2 = cpowl(r - csqrtl(H),1.0/3.0);
      sol[0] = creall(cs1 + cs2) - a2/3.0;
      /* per il calcolo degli zeri di una quartica basta uno zero reale qualsiasi quindi 
       * con questo flag si puo' ottimizzare calcolandone solo uno */
      s1ps2 = (cs1+cs2)/2.0; 
      s1ms2 = (cs1-cs2)/2.0;
      sol[1] = creall(-s1ps2 - a2/3.0 + I*sqrt3*s1ms2);
      sol[2] = creall(-s1ps2 - a2/3.0 - I*sqrt3*s1ms2);
      *numsol = 3;
    }
  else
    {
      s1 = cbrtl(r + sqrtl(H));
      s2 = cbrtl(r - sqrtl(H));
      sol[0] = s1 + s2 - a2/3.0;
      *numsol = 1;
    }
}
#ifdef LDLT_OPT_CUBIC
void solve_cubic_analytic_depressed_handle_inf(double b, double c, double *sol)
{
  double Q, R, theta, A, B;//Q3, R2;
  //const double sqrt32=sqrt(3)/2.0;
  double QR, QRSQ, KK, sqrtQ, RQ;
  //printf("qu?????????????????????????????\n");
  Q = -b/3.0;
  R = 0.5*c;
  //printf("Q=%.15G R=%.16G\n", Q, R);
  if (R==0)
    {
      //sol[0]=0;
      if (b <= 0)
	{
	  *sol=sqrt(-b);
	  //sol[2]=-sqrt(-b);
	}
      else
	{
	  *sol=0;
	  //sol[1]=I*sqrt(b);
	  //sol[2]=-I*sqrt(b);
	}
      return;
    }
  
  if (fabs(Q) < fabs(R))
    {
      //printf("Q < R\n");
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      //printf("Q >= R\n");
      RQ = R/Q;
      KK = copysign(1.0,Q)*(RQ*RQ/Q-1.0);
      //printf("RQ=%.15G sq=%.16G\n", RQ, RQ*RQ/Q-1.0);
      //printf("qui KK=%.15G Q=%.15G\n", KK, Q);
    }

  if (KK < 0.0)
    {
      //printf("KK < 0\n");
      sqrtQ=sqrt(Q);
      theta = acos((R/fabs(Q))/sqrtQ);
      //printf("theta=%.15G arg=%.15G sqrtQ=%.15G\n", theta, R/fabs(Q), sqrtQ);
#if 1
      if (theta < M_PI/2) 
	*sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
	*sol = -2.0*sqrtQ*cos((theta+2.0*M_PI)/3.0);
#else
      sol[0] = -2.0*sqrtQ*cos(theta/3.0);
      // to solve quartic we need only one real root 
      sol[1] = -2.0*sqrtQ*cos((theta+2.0*M_PI)/3.0);
      sol[2] = -2.0*sqrtQ*cos((theta-2.0*M_PI)/3.0);
#endif
    }
  else
    {
      //printf("KK >= 0\n");
      if (fabs(Q) < fabs(R))
	A = -copysign(1.0,R)*pow(fabs(R)*(1.0+sqrt(KK)),1.0/3.0);
      else
	{
	  A = -copysign(1.0,R)*pow(fabs(R)+sqrt(fabs(Q))*fabs(Q)*sqrt(KK),1.0/3.0);
	  //printf("A=%.15G arg=%.15G\n", A, sqrt(fabs(Q))*Q*sqrt(KK));
	}
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      *sol = A+B;
#if 0
      sol[1] = -0.5*(A+B)+I*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-I*sqrt32*(A-B);
#endif
    }
}

void solve_cubic_analytic_depressed(double b, double c, double *sol)
{
  double Q, R, theta, Q3, R2, A, B, sqrtQ;
  //const double sqrt32=sqrt(3)/2.0;
  //printf("qu?????????????????????????????\n");
  Q = -b/3.0;
  R = 0.5*c;
  /* 
     if (fabs(Q) > 1E102 o fabs(R) > 1E154)
     si ottengono degli inf quindi meglio usare la versione handle_inf
     printf("Q=%.16G R=%.16G\n", Q, R);
  */
  if (fabs(Q) > 1E102 || fabs(R) > 1E154)
    {
      solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  //sol[1]=sol[2]=0.0;
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sqrtQ=-2.0*sqrt(Q);
#if 1
      /* si puo' determinare facilmente la soluizioni più grande poiché basta studiare l'argomento del coseno */
      //t3 = (theta-2.0*M_PI)/3.0;
      if (theta < M_PI/2) 
	*sol = sqrtQ*cos(theta/3.0);
      else 
	*sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
      /* N.B. la terza soluzione è sempre più piccola delle altre due nell'intervallo
	 0 e Pi che è il range di angoli restituito da acos(...)
	 */ 
#else
      sol[0] = sqrtQ*cos(theta/3.0);
      // to solve quartic we need only one real root 
      sol[1] = sqrtQ*cos((theta+2.0*M_PI)/3.0);
      sol[2] = sqrtQ*cos((theta-2.0*M_PI)/3.0);
#endif
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      *sol = A+B; /* this is always largest root even if A=B */
#if 0
      sol[1] = -0.5*(A+B)+I*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-I*sqrt32*(A-B);
#endif
	
    }
}

#else
void solve_cubic_analytic_depressed_handle_inf(double b, double c, complex double sol[3])
{
  double Q, R, theta, A, B;//Q3, R2;
  const double sqrt32=sqrt(3)/2.0;
  double QR, QRSQ, KK, sqrtQ, RQ;
  //printf("qu?????????????????????????????\n");
  Q = -b/3.0;
  R = 0.5*c;
  //printf("Q=%.15G R=%.16G\n", Q, R);
  if (R==0)
    {
      sol[0]=0;
      if (b <= 0)
	{
	  sol[1]=sqrt(-b);
	  sol[2]=-sqrt(-b);
	}
      else
	{
	  sol[1]=I*sqrt(b);
	  sol[2]=-I*sqrt(b);
	}
      return;
    }
  
  if (fabs(Q) < fabs(R))
    {
      //printf("Q < R\n");
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      //printf("Q >= R\n");
      RQ = R/Q;
      KK = copysign(1.0,Q)*(RQ*RQ/Q-1.0);
      //printf("RQ=%.15G sq=%.16G\n", RQ, RQ*RQ/Q-1.0);
      //printf("qui KK=%.15G Q=%.15G\n", KK, Q);
    }

  if (KK < 0.0)
    {
      //printf("KK < 0\n");
      sqrtQ=sqrt(Q);
      theta = acos((R/fabs(Q))/sqrtQ);
      //printf("theta=%.15G arg=%.15G sqrtQ=%.15G\n", theta, R/fabs(Q), sqrtQ);
      sol[0] = -2.0*sqrtQ*cos(theta/3.0);
      // to solve quartic we need only one real root 
      sol[1] = -2.0*sqrtQ*cos((theta+2.0*M_PI)/3.0);
      sol[2] = -2.0*sqrtQ*cos((theta-2.0*M_PI)/3.0);
    }
  else
    {
      //printf("KK >= 0\n");
      if (fabs(Q) < fabs(R))
	A = -copysign(1.0,R)*pow(fabs(R)*(1.0+sqrt(KK)),1.0/3.0);
      else
	{
	  A = -copysign(1.0,R)*pow(fabs(R)+sqrt(fabs(Q))*fabs(Q)*sqrt(KK),1.0/3.0);
	  //printf("A=%.15G arg=%.15G\n", A, sqrt(fabs(Q))*Q*sqrt(KK));
	}
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = A+B;
      sol[1] = -0.5*(A+B)+I*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-I*sqrt32*(A-B);
    }
}

void solve_cubic_analytic_depressed(double b, double c, complex double sol[3])
{
  double Q, R, theta, Q3, R2, A, B, sqrtQ;
  const double sqrt32=sqrt(3)/2.0;
  //printf("qu?????????????????????????????\n");
  Q = -b/3.0;
  R = 0.5*c;
  if (fabs(Q) > 1E102 || fabs(R) > 1E154)
    {
      solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sqrtQ= -2.0*sqrt(Q);
      sol[0] = sqrtQ*cos(theta/3.0);
      // to solve quartic we need only one real root 
      sol[1] = sqrtQ*cos((theta+2.0*M_PI)/3.0);
      sol[2] = sqrtQ*cos((theta-2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = A+B;
      sol[1] = -0.5*(A+B)+I*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-I*sqrt32*(A-B);
    }
}
#endif
void solve_cubic_analytic(double *coeff, complex double sol[3])
{
  double a, b, c, Q, R, theta, Q3, R2, A, B;
  const double sqrt32=sqrt(3)/2.0;
  //printf("qu?????????????????????????????\n");
  a = coeff[2]/coeff[3];
  b = coeff[1]/coeff[3];
  c = coeff[0]/coeff[3];
  Q = (Sqr(a) - 3.0*b)/9.0;
  R = (2.0*Sqr(a)*a - 9.0*a*b + 27.0*c)/54.0;

  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sol[0] = -2.0*sqrt(Q)*cos(theta/3.0)- a/3.0;
      // to solve quartic we need only in real root 
      sol[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a/3.0;
      sol[2] = -2.0*sqrt(Q)*cos((theta-2.0*M_PI)/3.0) - a/3.0;
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = (A+B) - a/3.0;
      sol[1] = -0.5*(A+B)-a/3.0+I*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-a/3.0-I*sqrt32*(A-B);
    }
}

void solve_cubic(double *coeff, int *numsol, double sol[3])
{
  //const double sqrt3=sqrt(3.0);
  int ok;
  /* solution from Abramovitz */
  if (coeff[3] == 0.0)
    {
      printf("[WARNING] fallback to quadratic from cubic\n");
      //solve_quadratic(coeff, numsol, sol);
      solve_numrec(coeff, numsol, sol, &ok, 2);
      return;
    }
  /* la soluzione numnerica è comunque più accurata di quella analitica */
  solve_numrec(coeff, numsol, sol, &ok, 3);
}
#endif
double solcgbl[3];
int nscgbl;
#ifdef POLY_SOLVE_GSL
int solve_gslpoly(double c[5], int *numrealsol, double rsol[4])
{
  static gsl_poly_complex_workspace * w;
  double z[8];
  static int first=1;
  int k, status;
  if (first)
    {
      w  = gsl_poly_complex_workspace_alloc (5);
      first=0;
    }
 status=gsl_poly_complex_solve (c, 5, w, z);
 if (status!=GSL_SUCCESS)
   printf("[WARNING gsl_poly_complex_solve] probably not converging!\n");
 *numrealsol=0;
 for (k=0; k < 4; k++)
   {
     //printf("csol[%d]=%.15G+%.15G I\n", k, creal(csol[k]), cimag(csol[k]));
     if (z[2*k+1] == 0.0)
       {
	 //printf("cimag(csol[%d])=%.15G\n", k, cimag(csol[k]));
	 rsol[*numrealsol] = z[2*k];
	 (*numrealsol)++;
       }
   }
  //gsl_poly_complex_workspace_free (w);

  return 0;
}
#endif
void csolve_quartic_ferrari_cmplx(double *coeff, complex double sol[4])
{
  /* questa soluzione di fatto è quella di abramovitz */
  double a43, a2a3a4, a44, a4, a3, a2, a1, a0, a32, a12, a42, a3a4;
  double cb[4], p, q, r, m, cq[3];
  double complex solc[3], solq[2];
  complex double sm, A, B, C, Dp, Dm;
  const double sqrt2=1.4142135623730950488016887242097; 
  int k, mzero;
  a4 = coeff[4];
  a3 = coeff[3];
  a2 = coeff[2];
  a1 = coeff[1];
  a0 = coeff[0];
  if (a3==0 && a2==0 && a1==0 && a0 ==0)
    {
      sol[0] = 0.0;
      sol[1] = 0.0;
      sol[2] = 0.0;
      sol[3] = 0.0;
      return;
    }
  a32 = Sqr(a3);
  a12 = Sqr(a1);
  a42 = Sqr(a4);
  a3a4=-0.25*a3/a4;
  a43 = a42*a4;
  a2a3a4 = a2*a3*a4;
  a44 = a42*a42;
  p = (8.0*a2*a4-3.0*a32)/8.0/a42;
  q = (a32*a3 - 4.0*a2a3a4 + 8.0*a1*a42)/8.0/(a43);
  r = (-3.0*a32*a32+256.0*a0*a43-64.0*a1*a3*a42+16*a3*a2a3a4)/256.0/(a44);
  cb[3] = 8.0;
  cb[2] = 8.0*p;
  cb[1] = 2*Sqr(p)-8.0*r;
  cb[0] = -Sqr(q);
  solve_cubic_analytic(cb, solc);
  mzero=1;
  for (k=0; k < 3; k++)
    {
      if (cimag(solc[k])==0)
	{
	  if (creal(solc[k])!=0)
	    m = creal(solc[k]);
	  mzero=0;
	  break;
	}
    }
  if (mzero)
    {
      /* hence q=0 and quartic is a biquadratic */
      cq[2] = 1.0;
      cq[1] = p;
      cq[0] = r;
      solve_quadratic_cmplx(cq, solq);
      sol[0] = csqrt(solq[0]);
      sol[1] = -csqrt(solq[0]);
      sol[2] = csqrt(solq[1]);
      sol[3] = -csqrt(solq[1]);
      return;
    }
  sm = csqrt(m);
  A = 0.5*sqrt2*sm;
  B = sqrt2*q/sm;
  C = 2.0*(p+m);
  Dp = 0.5*csqrt(-(C + B));
  Dm = 0.5*csqrt(-(C - B));
  sol[0] = a3a4 + A + Dp;
  sol[1] = a3a4 + A - Dp;
  sol[2] = a3a4 - A + Dm;
  sol[3] = a3a4 - A - Dm;
}

void csolve_quartic_abramovitz_cmplx(double *coeff, complex double sol[4])
{
  double a3, a2, a1, a0, a32, a12;
  double lambda;
  double cb[4], y1;
  double complex solc[3], R, D, E, A, B;
  a3 = coeff[3]/coeff[4];
  a2 = coeff[2]/coeff[4];
  a1 = coeff[1]/coeff[4];
  a0 = coeff[0]/coeff[4];
  if (a3==0 && a2==0 && a1==0 && a0 ==0)
    {
      sol[0] = 0.0;
      sol[1] = 0.0;
      sol[2] = 0.0;
      sol[3] = 0.0;
      return;
    }

  a32 = Sqr(a3);
  a12 = Sqr(a1);
  lambda = 0.25*a3;
  cb[3] = 1.0;
  cb[2] = -a2;
  cb[1] = a1*a3-4.0*a0;
  cb[0] = 4.0*a2*a0-a12-a32*a0;
  solve_cubic_analytic(cb, solc);
  y1 = creal(solc[0]);
  R = csqrt(0.25*a32-a2+y1);
  if (R==0)
    {
      A = 0.75*a32-2.0*a2;
      B = 2.0*csqrt(Sqr(y1)-4*a0);
    }
  else
    {
      A = 0.75*a32-R*R-2.0*a2;
      B = 0.25*(4.0*a3*a2-8.0*a1-a32*a3)/R;
    }  
  D = csqrt(A+B);
  E = csqrt(A-B); 
  sol[0] = -0.25*a3 + 0.5*R + 0.5*D;
  sol[1] = -0.25*a3 + 0.5*R - 0.5*D;
  sol[2] = -0.25*a3 - 0.5*R + 0.5*E;
  sol[3] = -0.25*a3 - 0.5*R - 0.5*E;
}
void csolve_quartic_abramovitz(double *coeff, int *numrealsol, double rsol[4])
{
  int k;
  complex double csol[4];
  csolve_quartic_abramovitz_cmplx(coeff, csol);
  /* restituisce solo le soluzioni reali */
  *numrealsol=0;
  for (k=0; k < 4; k++)
    {
      //printf("csol[%d]=%.15G+%.15G I\n", k, creal(csol[k]), cimag(csol[k]));
      if (cimag(csol[k]) == 0.0)
	{
	  //printf("cimag(csol[%d])=%.15G\n", k, cimag(csol[k]));
	  rsol[*numrealsol] = creal(csol[k]);
	  (*numrealsol)++;
	}
    }
}
void solve_fourth_deg_cmplx(double *coeff, complex double sol[4])
{
  /* solution from H.E. Salzer, "A Note on Solution of Quartic Equations" Am. Math Society Proceedings, 279-281 (1959) */ 
  //const double ROUNDOFFERR = 1E-20;
  double x1, A, B, C, D, cb[4], m, n;
  double Asq, alpha, beta, rho, mp, np, m2;
  complex double gamma, delta, solc[3];
  A = coeff[3]/coeff[4];
  B = coeff[2]/coeff[4];
  C = coeff[1]/coeff[4];
  D = coeff[0]/coeff[4];
  if (A==0 && B==0 && C==0 && D ==0)
    {
      sol[0] = 0.0;
      sol[1] = 0.0;
      sol[2] = 0.0;
      sol[3] = 0.0;
      return;
    }

  Asq = Sqr(A);
  cb[3] = 1.0;
  cb[2] = -B;
  cb[1] = A*C-4*D;
  cb[0] = D*(4.0*B - Asq)-Sqr(C);
  solve_cubic_analytic(cb, solc);
  //printf("x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", cb[2], cb[1], cb[0]);
  x1 = creal(solc[0]);
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
  else // m is imaginary 
    {
      mp = sqrt(-m2);
      np = (A*x1 - 2.0*C)/(4.0*mp);
      alpha=0.5*Asq-x1-B;
      beta = 4.0*np-A*mp;
      rho = sqrt(Sqr(alpha)+Sqr(beta));
      gamma = csqrt((alpha + rho)/2.0); /* notare che alpha + rho è sempre positivo */
      if (gamma == 0.0)
	{
	  delta = csqrt(-alpha); 
	}
      else
	delta  = beta/gamma/2.0;
      sol[0] = (-0.5*A + gamma + I*(mp+delta))*0.5;
      sol[1] = conj(sol[0]);
      sol[2] = (-0.5*A - gamma + I*(mp-delta))*0.5;
      sol[3] = conj(sol[2]);
      return;
    }
  // m is real 
  alpha=0.5*Asq-x1-B;
  beta=4*n-A*m;
  gamma = csqrt(alpha+beta);
  sol[0]=-0.25*A+0.5*m+0.5*gamma;
  sol[1]=-0.25*A+0.5*m-0.5*gamma;
  delta = csqrt(alpha-beta);
  sol[2]=-0.25*A-0.5*m+0.5*delta;
  sol[3]=-0.25*A-0.5*m-0.5*delta;
}
void solve_fourth_deg(double *coeff, int *numsol, double sol[4])
{
  /* solution from H.E. Salzer, "A Note on Solution of Quartic Equations" Am. Math Society Proceedings, 279-281 (1959) */ 
  //const double ROUNDOFFERR = 1E-20;
  double x1, A, B, C, D, solc[3], cb[4], m, n;
  double Asq, alpha, beta, gamma, delta, rho, mp, np, m2;
  int nsc, kk;
#if 1
  if (coeff[4]==0)
    {
      //printf("[WARNING: solve_fourth_deg] coeff[4] must be different from 0!\n");
      solve_cubic(coeff, numsol, sol);
      return;
      //exit(-1);
    }
#endif
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
  solve_cubic(cb, &nsc, solc);
  //printf("x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", cb[2], cb[1], cb[0]);
#if 1
  x1 = solc[0];
#else
  if (nsc==1)
    {
      x1=solc[0];
    }
  else
    {
      dd0 = fabs(solc[0]-B);
      dd1 = fabs(solc[1]-B);
      dd2 = fabs(solc[2]-B);
      if (dd0 < dd1 && dd0 < dd2)
	x1 = solc[0];
      else if (dd1 < dd0 && dd1 < dd2)
	x1 = solc[1];
      else 
	x1 = solc[2];
    }
  //if (fabs(coeff[4]) < 0.01)
    //printf("solc[]=%.15G %.15G %.15G \n", solc[0], solc[1], solc[2]);
#endif
  m2 = Asq/4.0-B+x1;
  //printf("solc=%.15G numsol=%d m2=%.15G\n", solc[0], nsc, m2);
  nscgbl=nsc;
  for (kk=0; kk < nsc; kk++)
   solcgbl[kk] = solc[kk]; 
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
      if (mp + delta==0) /* in questo caso la parte immaginaria è zero e ho due soluzioni reali coincidenti */
	{
	  *numsol = 1;
	  sol[0] = (-A/2.0 + gamma)/2.0;
	}
      if (mp - delta == 0)
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
#if 0
  else if (fabs(alpha+beta) < ROUNDOFFERR)
    /* nel caso che alpha+beta è minore di zero ma di poco (entro gli errori di roundoff) 
     * allora significa che abbiamo zeri degeneri e quindi faccio il calcolo con i numeri complessi
     * e se la parte immaginaria è piccola allora considero la soluzione reale */
    {
      sol[0]=-0.25*A+0.5*m;
      sol[1]=-0.25*A+0.5*m;
      *numsol=2;
    }
#endif
  if (alpha-beta > 0)
    {    
      /* another pair or real solutions */
      delta = sqrt(alpha-beta);
      sol[*numsol]=-0.25*A-0.5*m+0.5*delta;
      sol[*numsol+1]=-0.25*A-0.5*m-0.5*delta;
      *numsol+=2;
    }
#if 0
  else if (fabs(alpha - beta) < ROUNDOFFERR)
    {
      sol[*numsol]=-0.25*A+0.5*m;
      sol[*numsol+1]=-0.25*A+0.5*m;
      *numsol+=2;
    }
#endif
}
void ellips2disk(double *solE, double *solD, double xEr, double yEr, double a, double b)
{
  solD[0] = 0;
  solD[1] = xEr + solE[0];
  solD[2] = yEr + solE[1];

  solD[1]*=a;
  solD[2]*=b;
  /* here we are again in the ppp reference system */
}

void projectback(double solarr[3], double rD[3], double nD[3])
{
  solarr[0] =  (nD[0]*rD[0] + nD[1]*rD[1] + nD[2]*rD[2] - nD[1]*solarr[1] - nD[2]*solarr[2])/nD[0]; 
}
void print_vec(char *lbl, double *V)
{
  printf("%s(%.15G,%.15G,%.15G)\n", lbl, V[0], V[1], V[2]);
}
double perpcompsq(double *V, double *C, double *n)
{
  int kk2;
  double dsc[3], sp, dscperp[3];
  for (kk2=0; kk2 < 3; kk2++)
    dsc[kk2] = V[kk2] - C[kk2]; 
  sp = scalProd(dsc, n);
  for (kk2=0; kk2 < 3; kk2++)
    dscperp[kk2] = dsc[kk2]-sp*n[kk2];
  return calc_normsq(dscperp);
}
double perpcomp(double *V, double *C, double *n)
{
  int kk2;
  double dsc[3], sp, dscperp[3];
  for (kk2=0; kk2 < 3; kk2++)
    dsc[kk2] = V[kk2] - C[kk2]; 
  sp = scalProd(dsc, n);
  for (kk2=0; kk2 < 3; kk2++)
    dscperp[kk2] = dsc[kk2]-sp*n[kk2];
  return calc_norm(dscperp);
}
long double perpcompl(long double *V, long double *C, long double *n)
{
  int kk2;
  long double dsc[3], sp, dscperp[3];
  for (kk2=0; kk2 < 3; kk2++)
    dsc[kk2] = V[kk2] - C[kk2]; 
  sp = scalProdl(dsc, n);
  for (kk2=0; kk2 < 3; kk2++)
    dscperp[kk2] = dsc[kk2]-sp*n[kk2];
  return calc_norml(dscperp);
}
double PowerM(double x, int n)
{
  double xsq;
  if (n==2)
    return x*x;
  else if (n==3)
    {
      xsq = x*x;
      return xsq*x;
    }
  else if (n==4)
    {
      xsq = x*x;
      return xsq*xsq;
    }
  else if (n==5)
    {
      xsq = x*x;
      return xsq*xsq*x;
    }
  else if (n==6)
    {
      xsq = x*x;
      return xsq*xsq*xsq;
    }
  else
    {
      printf("n too big max 6\n");
      exit(-1);
    }
}
extern int compare_func (const void *aa, const void *bb);
#if 0
int compare_func (const void *aa, const void *bb)
{
  double ai, bi, temp;
  ai = *((double *)aa);
  bi = *((double *)bb);
  if (ai > bi)
    return 1;
  else if (ai==bi)
    return 0;
  else
    return -1;
}
#endif
double polyquart(double x, double *coeff)
{
  double x4, x2, x3;
  x2=Sqr(x);
  x3=x2*x;
  x4=x2*x2;
  return coeff[4]*x4+coeff[3]*x3+coeff[2]*x2+coeff[1]*x+coeff[0];
}
void polyquartd(double x, double *coeff, double *fx, double *dfx)
{
  double x4, x2, x3;
  x2=Sqr(x);
  x3=x2*x;
  x4=x2*x2;
  *fx =coeff[4]*x4+coeff[3]*x3+coeff[2]*x2+coeff[1]*x+coeff[0];
  *dfx=4.0*coeff[4]*x3+3.0*coeff[3]*x2+2.0*coeff[2]*x+coeff[1];
}
double rtsafe(double c[5], double xg, double x1, double x2, double  xacc, int guess)
{
  /* xg is the initial guess and x1, x2 must bracket the solution */
#if 0
  p=c[n]*x+c[n-1];
  p1=c[n];
  for(i=n-2;i>=0;i--) {
    p1=p+p1*x;
    p=c[i]+p*x;
  }
  if (p1 == 0.0) throw("derivative should not vanish"); 
  x -= p/p1;
#endif
  const int MAXIT=100; //Maximum allowed number of iterations. Doub xh,xl;
  double fl, fh, xl, xh;
  double rts, dx, dxold, df, f, temp;
  int j;
  fl=polyquart(x1, c);
  fh=polyquart(x2, c);
#if 1
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) 
    {
      printf("[WARNING] Root must be bracketed in rtsafe\n");
      printf("fl=%.15G fh=%.15G\n", fl, fh);
      printf("xg=%.15G x1=%.15G x2=%.15G\n", xg, x1, x2);
      return xg;
    }
#endif
  if (fl == 0.0) 
    return x1;
  if (fh == 0.0) 
    return x2;
  if (fl < 0.0) 
    {
      xl=x1;
      xh=x2;
    } 
  else 
    {
      xh=x1;
      xl=x2; 
    }
  if (guess)
    rts = xg;
  else
    rts = 0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  polyquartd(rts,c,&f,&df);
  for (j=0;j<MAXIT;j++) 
    {
      //Orient the search so that f .xl/ < 0.
      //Initialize the guess for root, the “stepsize before last,” and the last step.
      //Loop over allowed iterations.
      if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) 
	  || (fabs(2.0*f) > fabs(dxold*df)) || (df==0)) 
	{ 
	  //Bisect if Newton out of range, 
	  //or not decreasing fast enough.
	  dxold=dx;
	  dx=0.5*(xh-xl);
	  rts=xl+dx;
	  if (xl == rts) 
	    return rts;
	} 
      else 
	{
	  dxold=dx;
	  dx=f/df;
	  temp=rts;
	  rts -= dx;
	  if (temp == rts) 
	    return rts;
	}
      if (fabs(dx) < xacc) 
	return rts;
      polyquartd(rts, c, &f, &df);
      //The one new function evaluation per iteration.
      if (f < 0.0) //Maintain the bracket on the root.
	xl=rts; 
      else
	xh=rts;
    }

  printf("[WARNING] Maximum number of iterations exceeded in rtsafe\n");
  return rts;
}
double diskdiskdiff(double *D, double *L, double Di[2][3], double Ci[3], double ni[3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, j2, kk; 
  double sp, normCiCj, CiCj[3], VV[3], Q1;
  double DiN, DjN, niN[3], njN[3], Djni, Djnj, assex[3], Dini, Pi[3], N[3], Pj[3], normN;
  double normNSq, PiPj[3], normPiPj, Q2, PjDj[3], normPiDi, normPjDj, PiDi[3];
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= 0.5*(L[0]+L[1]))
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
	    if (sp == 0 && calc_norm(VV) < 0.5*(D[0]+D[1]))
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
	    if (normPiDi <= 0.5*D[0] && normPjDj <= 0.5*D[1])
	      {
		Q1 = sqrt(Sqr(D[0])/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D[1])/4.0-Sqr(normPjDj));
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
  return 0;
}
double diskdisk(double D, double L, double Di[2][3], double Ci[3], double ni[3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, j2, kk; 
  double sp, normCiCj, CiCj[3], VV[3], Q1;
  double DiN, DjN, niN[3], njN[3], Djnj, Dini, Pi[3], N[3], Pj[3], normN;
  double normNSq, PiPj[3], normPiPj, Q2, PjDj[3], normPiDi, normPjDj, PiDi[3];
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  numcallsDD++;
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
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
#if 0
	    sp=0;
	    for (kk=0; kk < 3; kk++)
	      {
		sp  += Sqr(Di[j1][kk]-Dj[j2][kk]);
	      }	  
	    if (sp > D*D)
	      continue;
#endif
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
  return 0;
}
void balancel(long double a[4][4], int n)
{
  const long double RADIX=FLT_RADIX;// numeric_limits<Doub>::radix;
  int i, j;
  long double scale[4]={1.0,1.0,1.0,1.0};
  int done=0;
  long double r, c, g, f, s, sqrdx=RADIX*RADIX;
  while (!done) 
    {
      done=1;
      for (i=0;i<n;i++) 
	{
	  //Calculate row and column norms.
	  //If both are nonzero,
	  //find the integer power of the machine radix that comes closest to balancing the matrix.
	  r=0.0;
	  c=0.0;
	  for (j=0;j<n;j++)
	    if (j != i) 
	      {
		c += fabsl(a[j][i]);
		r += fabsl(a[i][j]);
	      }
	  if (c != 0.0 && r != 0.0) 
	    {
	      g=r/RADIX;
	      f=1.0;
	      s=c+r;
	      while (c<g) {
		f *= RADIX;
		c *= sqrdx;
	      }
	      g=r*RADIX;
	      while (c>g) 
		{
		  f /= RADIX;
		  c /= sqrdx; 
		}
	      if ((c+r)/f < 0.95*s) 
		{
		  done=0;
		  g=1.0/f;
		  scale[i] *= f;
		  for (j=0;j<n;j++) a[i][j] *= g; //Apply similarity transformation
		  for (j=0;j<n;j++) a[j][i] *= f;
		}
	    }
	}
    }
}
void hqrl(long double a[4][4], complex long double wri[4], int *ok, int n)
{
  int nn,m,l,k,j,its,i,mmin;
  long double z,y,x,w,v,u,t,s,r,q,p, anorm=0.0;
  const int MAXITS = 120;
  const long double EPS=1E-30;//3E-16;//numeric_limits<Doub >::epsilon();
  //const int n=4;
  for (i=0;i<n;i++)
    //Compute matrix no rm for possible use in lo- cating single small sub diagonal element.
    for (j=MAX(i-1,0);j<n;j++)
      anorm += fabsl(a[i][j]);
  nn=n-1;
  t=0.0;
  *ok = 1;
  //Gets changed only by an exceptional shift.
  while (nn >= 0) 
    {
      //Begin search for next eigenvalue.
      its=0;
      do 
	{
	  for (l=nn;l>0;l--)
	    {
	      //Begin iteration: look for single small sub di- agonal element.
	      s=fabsl(a[l-1][l-1])+fabsl(a[l][l]);
	      if (s == 0.0)
		s=anorm;
	
	      if (fabsl(a[l][l-1]) <= EPS*s)
		{
		  a[l][l-1] = 0.0;
		  break;
		}
	    }
	  x=a[nn][nn];
	  if (l == nn)
	    {
	      //One root found.  
	      wri[nn--]=x+t;
	    } 
	  else
	    {
	      y=a[nn-1][nn-1];
	      w=a[nn][nn-1]*a[nn-1][nn];
	      if (l == nn-1)
		{
		  //Two roots found...
		  p=0.5*(y-x);
		  q=p*p+w;
		  z=sqrtl(fabsl(q));
		  x += t;
		  if (q >= 0.0)
		    {
		      //...a real pair.
		      z=p+SIGNL(z,p);
		      wri[nn-1]=wri[nn]=x+z;
		      if (z != 0.0)
			wri[nn]=x-w/z;
		    } 
		  else
		    {
		      //...a complex pair.
		      wri[nn]=CMPLXL(x+p,-z);
		      wri[nn-1]=conjl(wri[nn]);
		    }
		  nn -= 2;
		} 
	      else
		{
		  //No roots found.  Continue iteration.
		  if (its == MAXITS)// il valore era 30 ma l'ho aumentato a 200 altrimenti non ce la fa...boh
		    {
		      printf("Too many iterations in hqr");
		      *ok = 0;
		      return;
		      //exit(-1);
		    }
		  /* NOTA: con questa modifica (ossia ho messo its % 10 invece di its==10 e its==20 come in Numerical Recipe) 
		   * mutuata dalle GSL l'algoritmo diventa un multishift e non si rischia
		   * che raggiunga MAXITS */
		  if (its % 10 == 0 && its > 0)
		    {
		      //Form exceptional shift.
		      t += x;
		      for (i=0;i<nn+1;i++)
			a[i][i] -= x;
		      s=fabsl(a[nn][nn-1])+fabsl(a[nn-1][nn-2]);
		      y=x=0.75*s;
		      w = -0.4375*s*s;
		    }
		  ++its;
		  for (m=nn-2;m>=l;m--)
		    {
		      //Form shift and then look for 2 consecutive small sub- diagonal elements.
		      z=a[m][m];
		      r=x-z;
		      s=y-z;
		      p=(r*s-w)/a[m+1][m]+a[m][m+1];
		      //Equation (W ebnote 16.21).
		      q=a[m+1][m+1]-z-r-s;
		      r=a[m+2][m+1];
		      s=fabsl(p)+fabsl(q)+fabsl(r);
		      //Scale to prevent over flow or under flow.
		      p /= s;
		      q /= s;
		      r /= s;
		      if (m == l) 
			break;
		      u=fabsl(a[m][m-1])*(fabsl(q)+fabsl(r));
		      v=fabsl(p)*(fabsl(a[m-1][m-1])+fabsl(z)+fabsl(a[m+1][m+1]));
		      /*  if (a1 * (fabs (q) + fabs (r)) <= GSL_DBL_EPSILON * fabs (p) * (a2 + a3))
			  break; */
 
		      if (u <= EPS*v)
			break;
		      //Equation (W ebnote 16.24).
		    }
		  for (i=m;i<nn-1;i++)
		    {
		      a[i+2][i]=0.0;
		      if (i != m) a[i+2][i-1]=0.0;
		    }
		  for (k=m;k<nn;k++)
		    {
		      //Double QR step on rows l to nn and columns m to nn .
		      if (k != m) 
			{
			  p=a[k][k-1];
			  //Begin setup of Householder vector.
			  q=a[k+1][k-1];
			  r=0.0;
			  if (k+1 != nn) 
			    r=a[k+2][k-1];
			  if ((x=fabsl(p)+fabsl(q)+fabsl(r)) != 0.0)
			    {
			      p /= x;
			      //Scale to prevent over flow or under flow.
			      q /= x;
			      r /= x;
			    }
			}
		      if ((s=SIGNL(sqrtl(p*p+q*q+r*r),p)) != 0.0)
			{
			  if (k == m) 
			    {
			      if (l != m)
				a[k][k-1] = -a[k][k-1];
			    } 
			  else
			    a[k][k-1] = -s*x;
			  p += s;
			  //Equations (Webnote 16.22).
			  x=p/s;
			  y=q/s;
			  z=r/s;
			  q /= p;
			  r /= p;
			  for (j=k;j<nn+1;j++)
			    {
			      //Row mo di cation.
			      p=a[k][j]+q*a[k+1][j];
			      if (k+1 != nn)
				{
				  p += r*a[k+2][j];
				  a[k+2][j] -= p*z;
				}
			      a[k+1][j] -= p*y;
			      a[k][j] -= p*x;
			    }
			  mmin = nn < k+3 ? nn : k+3;
			  for (i=l;i<mmin+1;i++)
			    {
			      //Column modification.
			      p=x*a[i][k]+y*a[i][k+1];
			      if (k+1 != nn) {
				p += z*a[i][k+2];
				a[i][k+2]
				  -= p*r;
			      }
			      a[i][k+1] -= p*q;
			      a[i][k] -= p;
			    }
			}
		    }
		}
	    }
	} 
      while (l+1 < nn);
    }
}

extern void balance(double a[4][4], int n);
#if 0
void balance(double a[4][4], int n)
{
  const double RADIX=FLT_RADIX;// numeric_limits<Doub>::radix;
  int i, j;
  double scale[4]={1.0,1.0,1.0,1.0};
  int done=0;
  double r, c, g, f, s, sqrdx=RADIX*RADIX;
  //const int n=4;
  while (!done) 
    {
      done=1;
      for (i=0;i<n;i++) 
	{
	  //Calculate row and column norms.
	  //If both are nonzero,
	  //find the integer power of the machine radix that comes closest to balancing the matrix.
	  r=0.0;
	  c=0.0;
	  for (j=0;j<n;j++)
	    if (j != i) 
	      {
		c += fabs(a[j][i]);
		r += fabs(a[i][j]);
	      }
	  if (c != 0.0 && r != 0.0) 
	    {
	      g=r/RADIX;
	      f=1.0;
	      s=c+r;
	      while (c<g) {
		f *= RADIX;
		c *= sqrdx;
	      }
	      g=r*RADIX;
	      while (c>g) 
		{
		  f /= RADIX;
		  c /= sqrdx; 
		}
	      if ((c+r)/f < 0.95*s) 
		{
		  done=0;
		  g=1.0/f;
		  scale[i] *= f;
		  for (j=0;j<n;j++) a[i][j] *= g; //Apply similarity transformation
		  for (j=0;j<n;j++) a[j][i] *= f;
		}
	    }
	}
    }
}
#endif
extern void hqr(double a[4][4], complex double wri[4], int *ok, int n);
#if 0
void hqr(double a[4][4], complex double wri[4], int *ok, int n)
{
  int nn,m,l,k,j,its,i,mmin;
  double z,y,x,w,v,u,t,s,r,q,p, anorm=0.0;
  const int MAXITS = 120;
  const double EPS=2.2204460492503131E-16;//3E-16;//numeric_limits<Doub >::epsilon();
  //const int n=4;
  for (i=0;i<n;i++)
    //Compute matrix no rm for possible use in lo- cating single small sub diagonal element.
    for (j=MAX(i-1,0);j<n;j++)
      anorm += fabs(a[i][j]);
  nn=n-1;
  t=0.0;
  *ok = 1;
  //Gets changed only by an exceptional shift.
  while (nn >= 0) 
    {
      //Begin search for next eigenvalue.
      its=0;
      do 
	{
	  for (l=nn;l>0;l--)
	    {
	      //Begin iteration: look for single small sub di- agonal element.
	      s=fabs(a[l-1][l-1])+fabs(a[l][l]);
	      if (s == 0.0)
		s=anorm;
	
	      if (fabs(a[l][l-1]) <= EPS*s)
		{
		  a[l][l-1] = 0.0;
		  break;
		}
	    }
	  x=a[nn][nn];
	  if (l == nn)
	    {
	      //One root found.  
	      wri[nn--]=x+t;
	    } 
	  else
	    {
	      y=a[nn-1][nn-1];
	      w=a[nn][nn-1]*a[nn-1][nn];
	      if (l == nn-1)
		{
		  //Two roots found...
		  p=0.5*(y-x);
		  q=p*p+w;
		  z=sqrt(fabs(q));
		  x += t;
		  if (q >= 0.0)
		    {
		      //...a real pair.
		      z=p+SIGN(z,p);
		      wri[nn-1]=wri[nn]=x+z;
		      if (z != 0.0)
			wri[nn]=x-w/z;
		    } 
		  else
		    {
		      //...a complex pair.
		      wri[nn]=CMPLX(x+p,-z);
		      wri[nn-1]=conj(wri[nn]);
		    }
		  nn -= 2;
		} 
	      else
		{
		  //No roots found.  Continue iteration.
		  if (its == MAXITS)// il valore era 30 ma l'ho aumentato a 200 altrimenti non ce la fa...boh
		    {
		      printf("Too many iterations in hqr");
		      *ok = 0;
		      return;
		      //exit(-1);
		    }
		  /* NOTA: con questa modifica (ossia ho messo its % 10 invece di its==10 e its==20 come in Numerical Recipe) 
		   * mutuata dalle GSL l'algoritmo diventa un multishift e non si rischia
		   * che raggiunga MAXITS */
		  if (its % 10 == 0 && its > 0)
		    {
		      //Form exceptional shift.
		      t += x;
		      for (i=0;i<nn+1;i++)
			a[i][i] -= x;
		      s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
		      y=x=0.75*s;
		      w = -0.4375*s*s;
		    }
		  ++its;
		  for (m=nn-2;m>=l;m--)
		    {
		      //Form shift and then look for 2 consecutive small sub- diagonal elements.
		      z=a[m][m];
		      r=x-z;
		      s=y-z;
		      p=(r*s-w)/a[m+1][m]+a[m][m+1];
		      //Equation (W ebnote 16.21).
		      q=a[m+1][m+1]-z-r-s;
		      r=a[m+2][m+1];
		      s=fabs(p)+fabs(q)+fabs(r);
		      //Scale to prevent over flow or under flow.
		      p /= s;
		      q /= s;
		      r /= s;
		      if (m == l) 
			break;
		      u=fabs(a[m][m-1])*(fabs(q)+ fabs(r));
		      v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
		      /*  if (a1 * (fabs (q) + fabs (r)) <= GSL_DBL_EPSILON * fabs (p) * (a2 + a3))
			  break; */
 
		      if (u <= EPS*v)
			break;
		      //Equation (W ebnote 16.24).
		    }
		  for (i=m;i<nn-1;i++)
		    {
		      a[i+2][i]=0.0;
		      if (i != m) a[i+2][i-1]=0.0;
		    }
		  for (k=m;k<nn;k++)
		    {
		      //Double QR step on rows l to nn and columns m to nn .
		      if (k != m) 
			{
			  p=a[k][k-1];
			  //Begin setup of Householder vector.
			  q=a[k+1][k-1];
			  r=0.0;
			  if (k+1 != nn) 
			    r=a[k+2][k-1];
			  if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0)
			    {
			      p /= x;
			      //Scale to prevent over flow or under flow.
			      q /= x;
			      r /= x;
			    }
			}
		      if ((s=SIGN(sqrt(p*p+q*q+ r*r),p)) != 0.0)
			{
			  if (k == m) 
			    {
			      if (l != m)
				a[k][k-1] = -a[k][k-1];
			    } 
			  else
			    a[k][k-1] = -s*x;
			  p += s;
			  //Equations (Webnote 16.22).
			  x=p/s;
			  y=q/s;
			  z=r/s;
			  q /= p;
			  r /= p;
			  for (j=k;j<nn+1;j++)
			    {
			      //Row mo di cation.
			      p=a[k][j]+q*a[k+1][j];
			      if (k+1 != nn)
				{
				  p += r*a[k+2][j];
				  a[k+2][j] -= p*z;
				}
			      a[k+1][j] -= p*y;
			      a[k][j] -= p*x;
			    }
			  mmin = nn < k+3 ? nn : k+3;
			  for (i=l;i<mmin+1;i++)
			    {
			      //Column modification.
			      p=x*a[i][k]+y*a[i][k+1 ];
			      if (k+1 != nn) {
				p += z*a[i][k+2];
				a[i][k+2]
				  -= p*r;
			      }
			      a[i][k+1] -= p*q;
			      a[i][k] -= p;
			    }
			}
		    }
		}
	    }
	} 
      while (l+1 < nn);
    }
}
#endif
#define MR 10
#define MT 8
void laguer(complex double a[], int m, complex double *x, int *its)
/*m i Given the degree m and the m+1 complex coefficients a[0..m] of the polynomial i=0 a[i]x ,
and given a complex value x, this routine improves x by Laguerre’s method until it converges, within the achievable roundoff limit, to a root of the given polynomial. The number of iterations taken is returned as its.*/
{
  double EPSS=2.2204460492503131E-16;
  int iter,j;
  int MAXIT=MR*MT;
  double abx,abp,abm,err;
  complex double dx,x1,b,d,f,g,h,sq,gp,gm,g2;
  static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0}; //Fractions used to break a limit cycle.
  for (iter=1;iter<=MAXIT;iter++) 
    { 
      *its=iter;
      b=a[m];
      err=cabs(b); 
      d=f=CMPLX(0.0,0.0); 
      abx=cabs(*x);
      for (j=m-1;j>=0;j--) 
	{
  	  f=*x*f+d; 
	  d=*x*d+b; 
	  b=*x*b+a[j]; 
	  err=cabs(b)+abx*err;
	  //Loop over iterations up to allowed maximum.
	  //Efficient computation of the polynomial and its first two derivatives. f stores P′′/2.
	}
      err *= EPSS;
      //Estimate of roundoff error in evaluating polynomial.
      if (cabs(b) <= err) return; //We are on the root.
      g=d/b; //The generic case: use Laguerre’s formula. 
      g2=g*g;
      h=g2-2.0*(f/b);
      sq=csqrt((double) (m-1)*((((double) m)*h)-g2)); 
      gp=g+sq;
      gm=g-sq;
      abp=cabs(gp);
      abm=cabs(gm);
      if (abp < abm) gp=gm;
      dx=FMAX(abp,abm) > 0.0 ? CMPLX(((double) m),0.0)/gp
	   : (1.0+abx)*CMPLX(cos((double)iter),sin((double)iter)); 
      x1=*x-dx;
      if (*x == x1) 
	return; //Converged.
      if (iter % MT) 
	*x=x1;
      else 
	*x=*x-(frac[iter/MT]*dx);
      //Every so often we take a fractional step, to break any limit cycle (itself a rare occur- rence).
    }
  printf("too many iterations in laguer");
  //Very unusual — can occur only for complex roots. Try a different starting guess for the root.
  return;
}

void zroots(complex double *a, int m, complex double *radici, int polish, double EPS)
/* Copy of coefficients for successive deflation. Loop over each root to be found.
Start at zero to favor convergence to small-
est remaining root, and find the root. Forward deflation.
 Given the degree m and the m+1 complex coefficients a[0..m] of the polynomial
this routine successively calls laguer and finds all m complex roots in roots[1..m]. The boolean variable polish should be input as true (1) if polishing (also by Laguerre’s method) is desired, false (0) if the roots will be subsequently polished by other means.*/
{
  //void laguer(fcomplex a[], int m, fcomplex *x, int *its); 
  int i,its,j,jj;
  complex double roots[5];
  complex double x,b,c,ad[20];
  for (j=0;j<=m;j++) ad[j]=a[j]; 
  for (j=m;j>=1;j--) 
    {
      x=CMPLX(0.0,0.0);
      laguer(ad,j,&x,&its);
      if (fabs(cimag(x)) <= 2.0*EPS*fabs(creal(x))) 
	x=CMPLX(creal(x),0.0); 
      roots[j]=x;
      b=ad[j];
      for (jj=j-1;jj>=0;jj--) {
	c=ad[jj];
	ad[jj]=b; 
	b=x*b+c;
      }
    }
  if (polish)
    for (j=1;j<=m;j++)
      laguer(a,m,&roots[j],&its); 
  for (j=2;j<=m;j++) 
    {
      x=roots[j];
      for (i=j-1;i>=1;i--) 
	{
	  if (creal(roots[i]) <= creal(x)) break;
	  roots[i+1]=roots[i];
	}
	roots[i+1]=x;
    }
  for (j=0; j < m; j++)
    {
      radici[j] = roots[j+1];
    }
}

#ifdef USE_LAPACK
#ifndef MD_MAC
extern void dgebal_(char *, int *, double *, int *, int *, int *, double *, int *);
extern void dhseqr_(char *, char*, int *, int *, int *, double *, int *, double *, double *, double *,
		    int *, double *, int *, int *);
#endif
/* int dhseqr_(char *__job, char *__compz, __CLPK_integer *__n,
        __CLPK_integer *__ilo, __CLPK_integer *__ihi, __CLPK_doublereal *__h__,
        __CLPK_integer *__ldh, __CLPK_doublereal *__wr, __CLPK_doublereal *__wi,
        __CLPK_doublereal *__z__, __CLPK_integer *__ldz,
        __CLPK_doublereal *__work, __CLPK_integer *__lwork,
        __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
        __IPHONE_4_0); */

void wrap_dgebal(double a[4][4], int *ilo, int *ihi, int n, int *ok)
{
  char JOB='S', COMPZ='N';
  double AT[4*4], z[4*4];
  int i, j, c1, c2, k, N;
  for (i=0; i<n; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<n; j++) AT[j+n*i]=a[j][i];		
    }						
  c1 = 4;
  c2 = 4;
  N=n;
 //dgesv_(&c1, &c2, AT, &c1, pivot, &x, &c1, ok); 
  dgebal_(&JOB, &N, AT, &c1, ilo, ihi, z, ok);      
  for (i=0; i<n; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<n; j++) a[j][i] = AT[j+n*i];		
    }	
}
void wrap_hseqr(double a[4][4], double *wr, double *wi, int n, int ilo, int ihi, int *ok)
{
  char JOB='E', COMPZ='N';
  double AT[4*4], z[4*4],  work[4];
  int i, j, c1, c2, k, pivot[4], lwork, N;
  for (i=0; i<n; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<n; j++) AT[j+n*i]=a[j][i];		
    }						
  c1 = 4;
  c2 = 4;
  lwork=4;
  N=n;
  dhseqr_(&JOB, &COMPZ, &N, &ilo, &ihi, AT, &c1, wr, wi, z, &c2, work, &lwork, ok);      
}
void wrap_dgeev(double a[4][4], double *wr, double *wi, int n, int *ok)
{
  char JOB='N', COMPZ='N';
  double AT[4*4], z[4*4], z2[4*4], work[40], work2[4];
  int i, j, c1, c2, k, lwork, N;
  for (i=0; i<n; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<n; j++) AT[j+n*i]=a[j][i];		
    }						
  c1 = 4;
  c2 = 4;
  lwork=40;
  N=n;
  dgeev_(&JOB, &COMPZ, &N, AT, &c1, wr, wi, z, &c2, z2, &c2, work, &lwork, ok);      
}

#endif
void QRfactorizationl(long double hess[4][4], complex long double sol[4], int *ok, int n)
{
  //int ok;
  /* pagina 615 Num. Rec. */  
  /* ora funziona si doveva solo aumentare il numero massimo d'iterazioni */
  balancel(hess, n);
  hqrl(hess, sol, ok, n);
}
extern void QRfactorization(double hess[4][4], complex double sol[4], int *ok, int n);
#if 0
void QRfactorization(double hess[4][4], complex double sol[4], int *ok, int n)
{
  //int ok;
  double zr[4], zi[4];
  int ilo, ihi,k;
  /* pagina 615 Num. Rec. */  
#ifdef USE_LAPACK
#if 1
  wrap_dgebal(hess, &ilo, &ihi, n, ok);
  wrap_hseqr(hess, zr, zi, n, ilo, ihi, ok);
  *ok = 1;
#else
  wrap_dgeev(hess, zr, zi, 4, &ok);
#endif
  //printf("hess=%.15G %.15G %.15G %.15G\n", hess[0][0], hess[0][1], hess[0][2], hess[0][3]);
  for (k=0; k < n; k++)
    {
      //printf("PRIMA QRDECOMP csol[k=%d]=%.15G+%.15G I\n", k, zr[k], zi[k]);
      sol[k] =CMPLX(zr[k],zi[k]);
    }
#else
  /* ora funziona si doveva solo aumentare il numero massimo d'iterazioni */
  balance(hess, n);
  hqr(hess, sol, ok, n);
  //sort(); 
#endif
}
#endif
int tinyimagGBL=0;

void solve_numrecl(long double coeff[5], int *numrealsol, long double rsol[4], int *ok, int m)
{
//void zrhqr(VecDoub_I &a, VecComplex_O &rt) Pm i
  /*Find all the roots of a polynomial with real coefficients, a4*x^4+a3*x^3+a2*x^2+a1*x+a0, 
   * given the coefficients a[0..m]. The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pagina 497 Num. Rec. */
  complex long double csol[4]; 
  //const int m=4;
  const double TINYEPS=1E-7;
  long double hess[4][4];
  int j, k, smallimag;
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorizationl(hess, csol, ok, m);
  
  //for (j=0;j<m;j++)
    //rt[j]=h.wri[j];
  *numrealsol=0;
  //printf("{%.16G,%.16G,%.16G,%.16G,%.16G}\n", coeff[0],coeff[1], coeff[2], coeff[3], coeff[4]);
  smallimag=0;
  for (k=0; k < m; k++)
   {
     //printf("QRDECOMP csol[%d]=%.15G+%.15G I\n", k, creal(csol[k]), cimag(csol[k]));
     if (cimagl(csol[k]) == 0.0)
       {
	 //printf("cimag(csol[%d])=%.15G\n", k, cimag(csol[k]));
	 rsol[*numrealsol] = creall(csol[k]);
	 (*numrealsol)++;
       }
     else if (fabsl(cimagl(csol[k])) < TINYEPS*fabsl(creall(csol[k])))
       {
	 smallimag=1;
       }
   }
  if (smallimag==1 && *numrealsol==0)
    {
      tinyimagGBL=1;
    }
  else
    tinyimagGBL=0;
  //gsl_poly_complex_workspace_free (w);
}
#if 0
void solve_numrec(double coeff[5], int *numrealsol, double rsol[4], int *ok, int m)
{
//void zrhqr(VecDoub_I &a, VecComplex_O &rt) Pm i
  /*Find all the roots of a polynomial with real coefficients, a4*x^4+a3*x^3+a2*x^2+a1*x+a0, 
   * given the coefficients a[0..m]. The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pagina 497 Num. Rec. */
  complex double csol[4], cc[5]; 
  const double TINYEPS=1E-7;
  double hess[4][4];
  int j, k, smallimag;
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
#ifdef USE_LAGUERRE
  for(k=0; k < m+1; k++)
    cc[k]=CMPLX(coeff[k],0.0);
  zroots(cc, csol, 1, m, 1E-14);
#else
  QRfactorization(hess, csol, ok, m);
#endif
  //for (j=0;j<m;j++)
    //rt[j]=h.wri[j];
  *numrealsol=0;
  smallimag=0;
  //printf("{%.16G,%.16G,%.16G,%.16G,%.16G}\n", coeff[0],coeff[1], coeff[2], coeff[3], coeff[4]);
  for (k=0; k < m; k++)
   {
#ifdef MC_DEBUG_HCALGO
     //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
       //printf("QRDECOMP csol[%d]=%.15G+%.15G I\n", k, creal(csol[k]), cimag(csol[k]));
#endif
     if (cimag(csol[k]) == 0.0)
       {
	 //printf("cimag(csol[%d])=%.15G\n", k, cimag(csol[k]));
	 rsol[*numrealsol] = creal(csol[k]);
	 (*numrealsol)++;
       }
     else if (fabs(cimag(csol[k])) < TINYEPS*fabs(creal(csol[k])))
       {
	 smallimag=1;
       }
   }
  if (smallimag==1 && *numrealsol==0)
    tinyimagGBL=1;
  else
    {
      tinyimagGBL=0;
    }
  //gsl_poly_complex_workspace_free (w);
}
#endif
/* NOTA: nell'articolo dice di usare un buffer di 4 quindi in tutto 5 elementi, tuttavia
 * stando ai risultati mostrati sembra che abbia usato un buffer più grande che in effetti
 * garantisce una migliore convergenza dei casi ciclici */
#define FQSBUFLEN 12
int eps_identical(double *eps)
{
  int j;
  for (j=1; j < FQSBUFLEN; j++)
    {
      /* if actual eps is identical to one of the four previous ones then terminate! */
      if (eps[0]==eps[j])
	return 1;
    }
  // if it gets here they are all equal!
  return 0;
}
/* ============================= FAST QUARTIC SOLVER STROBACH ===================================== */
#ifdef FAST_QUARTIC_SOLVER
void two_lin_eqs(double fmat[2][2],double evec[2], double dalf[2])
{ 
  double r,x1,x2,c,s;

  r=sqrt(Sqr(fmat[0][0])+Sqr(fmat[1][0]));
  if(r>0.0)
    {
      c=fmat[0][0]/r;
      s=fmat[1][0]/r;
    }
  else
    {
      c=1.0;
      s=0.0;
    }

  x1=c*fmat[0][0]+s*fmat[1][0];
  x2=-s*fmat[0][0]+c*fmat[1][0];
  fmat[0][0]=x1;
  fmat[1][0]=x2;

  x1=c*fmat[0][1]+s*fmat[1][1];
  x2=-s*fmat[0][1]+c*fmat[1][1];
  fmat[0][1]=x1;
  fmat[1][1]=x2;      

  x1=c*evec[0]+s*evec[1];
  x2=-s*evec[0]+c*evec[1];
  evec[0]=x1;
  evec[1]=x2; 

  if(fmat[1][1]==0.0)
    dalf[1]=0.0;
  else
    dalf[1]=evec[1]/fmat[1][1];

  if(fmat[0][0]==0.0)
    dalf[0]=0.0;
  else
    dalf[0]=(evec[0]-fmat[0][1]*dalf[1])/fmat[0][0];
}  

void quadratic(double aa,double bb,double cc, double dd, double a, double b, complex double roots[2])
{ 
  double diskr,div,zmax,zmin,evec[2], fmat[2][2],dpar[2],at,bt,err,errt;
  int iter; 

  //-------------------------------- parameter backward correction step:       

  evec[0]=bb*b-b*b-a*b*aa+a*a*b-dd;            // equation (5.18)
  evec[1]=cc*b-b*b*aa+b*b*a-a*dd;              // equation (5.18)
  err=fabs(evec[0])+fabs(evec[1]);                 // equation (5.23)

  if(err!=0.0)    // terminate
    {      
      for (iter=0; iter < 8; iter++) 
	{ 
	  fmat[0][0]=-b*aa+2*a*b;                         // equation (5.19)
	  fmat[0][1]=bb-2*b-a*aa+a*a;                    // equation (5.19)
	  fmat[1][0]=b*b-dd;                             // equation (5.19)
	  fmat[1][1]=cc-2*b*aa+2*b*a;                     // equation (5.19)

	  evec[0]=-evec[0];
	  evec[1]=-evec[1];

	  two_lin_eqs(fmat,evec,dpar);              // equation (5.20)

	  at=a;
	  bt=b;

	  a=a+dpar[0];                                   // equation (5.21)
	  b=b+dpar[1];                                   // equation (5.22)

	  evec[0]=bb*b-b*b-a*b*aa+a*a*b-dd;            // equation (5.18)
	  evec[1]=cc*b-b*b*aa+b*b*a-a*dd;              // equation (5.18)

	  errt=err;
	  err=fabs(evec[0])+fabs(evec[1]);                 // equation (5.23) 

	  if(err==0.0)
	    break;	// terminate

	  if(err>=errt)   // terminate without parameter update
	    {
	      a=at;
	      b=bt;  
	      break;
	    }

	}
    } 

  //---------------------------------------- solve a quadratic equation:     

  diskr=a*a-4*b;   

  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
	zmin=0.0;
      else
	zmin=b/zmax;

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      roots[0]=CMPLX(-a/2,sqrt(-diskr)/2);
      roots[1]=CMPLX(-a/2,-sqrt(-diskr)/2);      
    }   
} 
      
/* qui ottendo c e a con il metodo dei minimi quadrati (i calcoli sono più semplici di quelli di Strobach
 * poiché faccio un semplice least-square fit mentre lui fa un weighted least square fit) */
void ac_fit_alt(double a, double b, double c, double d, double *AQ, double *BQ, double *CQ, double *DQ)
{
  double aq, bq, cq, dq;

  aq = *AQ;
  bq = *BQ;
  cq = *CQ;
  dq = *DQ; 

  if (fabs(*AQ) > fabs(*CQ))
    cq = (-aq + a - aq*bq + aq*b + bq*c - aq*dq - aq*bq*dq)/(1 + aq*aq + bq*bq);
  else
    aq = (a - cq - bq*cq + b*cq - cq*dq - bq*cq*dq + c*dq)/(1 + cq*cq + dq*dq); 
  *AQ=aq;
  *BQ=bq;
  *CQ=cq;
  *DQ=dq;
}    

void  ac_fit(double a, double b, double c, double d, double *AQ, double *BQ, double *CQ, double *DQ)
{
  double cmat[3][2],x1,x2,r,cr,sr; 
  double aq, bq, cq, dq;
  int j;

  aq = *AQ;
  bq = *BQ;
  cq = *CQ;
  dq = *DQ;
  if (fabs(aq) > fabs(cq))
    {                    // equation (5.12)
      cmat[0][0]=aq;
      cmat[0][1]=b-bq-dq;      
      cmat[1][0]=bq;
      cmat[1][1]=c-aq*dq;      
      cmat[2][0]=1.0;
      cmat[2][1]=a-aq;

      x1=cmat[0][0];
      x2=cmat[2][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2)>fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1; 
	    }
	}     

      for (j=0; j < 2; j++)
	{
	  x1=cr*cmat[0][j]+sr*cmat[2][j];
	  x2=-sr*cmat[0][j]+cr*cmat[2][j];
	  cmat[0][j]=x1;
	  cmat[2][j]=x2;
	}

      x1=cmat[0][0];
      x2=cmat[1][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2) > fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1; 
	    }
	}

      for (j=0; j < 2; j++) 
	{
	  x1=cr*cmat[0][j]+sr*cmat[1][j];
	  x2=-sr*cmat[0][j]+cr*cmat[1][j];
	  cmat[0][j]=x1;
	  cmat[1][j]=x2;
	}

      cq=cmat[0][1]/cmat[0][0];
    } 
  else                                          // equation (5.13)
    {
      cmat[0][0]=cq;
      cmat[0][1]=b-bq-dq;      
      cmat[1][0]=dq;
      cmat[1][1]=c-bq*cq;      
      cmat[2][0]=1.0;
      cmat[2][1]=a-cq;

      x1=cmat[0][0];
      x2=cmat[2][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2) > fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1;
	    } 
	}     

      for (j=0; j < 2; j++)
	{
	  x1=cr*cmat[0][j]+sr*cmat[2][j];
	  x2=-sr*cmat[0][j]+cr*cmat[2][j];
	  cmat[0][j]=x1;
	  cmat[2][j]=x2;
	}

      x1=cmat[0][0];
      x2=cmat[1][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2) > fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1; 
	    }
	}   

      for (j=0; j < 2; j++) 
	{
	  x1=cr*cmat[0][j]+sr*cmat[1][j];
	  x2=-sr*cmat[0][j]+cr*cmat[1][j];
	  cmat[0][j]=x1;
	  cmat[1][j]=x2;
	}

      aq=cmat[0][1]/cmat[0][0];
    } 
  *AQ=aq;
  *BQ=bq;
  *CQ=cq;
  *DQ=dq;
  //endif   

  //      return
  //      end
}
void r_quadratic(double a,double b,double c, double *xmax,double *xmin,int *iflag) // real solution: iflag=1
                                                    // no real sol:   iflag=0
{
  double diskr,nenn;

  diskr=b*b-4*a*c;

  if(diskr < 0.0)
    {
      *iflag=0;
      *xmax=-b/(2*a);
      *xmin=-b/(2*a);
    }
  else
    {
      *iflag=1;

      if(b > 0.0)
	nenn=-b-sqrt(diskr);
      else
	nenn=-b+sqrt(diskr);

      if(nenn==0.0 || a==0.0)
	{	
	  *xmax=0.0;
	  *xmin=0.0;
	}
      else
	{
	  *xmax=nenn/(2*a);
	  *xmin=2*c/nenn;
	}
    }
}
//====================================================================

void depressed_cubic_root(double g, double h, double *root)
{
  double x0,xt0,xr,xmax,xmin,a,b,c,xabs,oldabs;
  int iter, iflag;    
  xabs=0.0;

  if(h==0.0)
    {
      x0=0.0;
    }
  else
    {
      //-------------------------------------------------------------------- 
      if(h<0.0)
	{

	  xr=sqrt(-h);                                  // equation (3.20)

	  if(g > xr)
	    x0=-h/g;                                       // equation (3.26)
	  else if(g < -xr*xr)
	    x0=sqrt(-g);                                  // equation (3.22)
	  else
	    x0=xr;                                         // equation (3.24)

	}
      else if(h > 0.0)
	{
	  xr=sqrt(h);                                   // equation (3.21)

	  if(g > xr)
	    x0=-h/g;                                       // equation (3.26)
	  else if(g < -xr*xr)
	    x0=-sqrt(-g);                                 // equation (3.23)
	  else
	    x0=-xr;                                        // equation (3.25)

	}

      //-------------------------------      
      //do iter=1,8
      for (iter=0; iter < 8; iter++)	
	{
	  a=x0*x0;
	  b=-h;
	  c=g*x0*x0+2*h*x0;

	  r_quadratic(a,b,c,&xmax,&xmin,&iflag);      // equation (3.10)

	  x0=xmin;                                       // equation (3.11)
	  xt0=x0;

	  a=2*x0;
	  b=g-x0*x0;
	  c=h;

	  r_quadratic(a,b,c,&xmax,&xmin,&iflag);       // equation (3.15)

	  if(h < 0.0)
	    {     // take the positive solution
	      if(xmax>=0.0)
		x0=xmax;                                       // equation (3.16)
	      else
		x0=xmin;                                       // equation (3.16)

	    }

	  if(h > 0.0)
	    {     
	      // take the negative solution 
	      if(xmax <= 0.0)
		x0=xmax;                                       // equation (3.17)
	      else
		x0=xmin;                                       // equation (3.17)

	    }       

	  oldabs=xabs;
	  xabs=fabs(xt0-x0);
	  if(iter > 0 && xabs == oldabs)   // terminate
	    { 
	      //fine = 1;
	      //printf("qui1 x0=%f\n", x0);
	      break; //goto 100
	    }
	  if(iter > 0  && xabs == 0.0)      // terminate
	    { 
	      //printf("qui2\n");
	      //fine = 1;
	      break;
	      //goto 100
	    }     

	}//enddo

      *root=x0;

      //--------------------------------------------------------------------      
    }
  //--------------------------------------------------------------------      

  //100   continue

    //*root=x0    
} 
void  cubic_B_shift(double a, double b, double c, double d, double *phi0)
{
  double g,h,gg,hh,aq,bq,cq,dq,s,diskr;

  // !------------------------------------------------------- the B-shift:

  diskr=9*a*a-24*b;                    //         ! equation (3.58)
      
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
	s=-2*b/(3*a+diskr);                            // equation (3.58)
      else
	s=-2*b/(3*a-diskr);                            // equation (3.58)
	  
    }
  else
    {      
      s=-a/4;                                    // equations (3.59-60)
    }
      
  // !--------------------------- the shift transformation (Horner forms):

  aq=a+4*s;                                      // equation (3.45)
  bq=b+3*s*(a+2*s);                              // equation (3.46)
  cq=c+s*(2*b+s*(3*a+4*s));                      // equation (3.47)
  dq=d+s*(c+s*(b+s*(a+s)));                      // equation (3.48)     

  gg=bq*bq/9;
  hh=aq*cq;      
  g=hh-4*dq-3*gg;                                // equation (3.49)
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq;         // equation (3.50)

//------------------ call the parabola/reciprocal intersection method:      
      
  depressed_cubic_root(g,h,phi0);     
      
// -------------------------------------------------------------------- 
      
} 
double calc_res(double d2, double l2, double del2, double dml3l3, double d2eq46)
{
  return fabs(d2eq46-d2)+fabs(del2-2.0*l2*d2)+fabs(dml3l3-d2*l2*l2);
}
void solve_numrecl_cmplx(long double *coeff, complex long double *csol, int m)
{
//void zrhqr(VecDoub_I &a, VecComplex_O &rt) Pm i
  /*Find all the roots of a polynomial with real coefficients, a4*x^4+a3*x^3+a2*x^2+a1*x+a0, 
   * given the coefficients a[0..m]. The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pagina 497 Num. Rec. */
  //const int m=4;
  //const double TINYEPS=1E-7;
  long double hess[4][4]; //coeff[5];
  int j, k, ok;
  //for (k=0; k < 5; k++)
    //coeff[k] = (long double) coeffA[k];
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorizationl(hess, csol, &ok, m);
}

/* ============================================================================================ */  
void myquadratic(double aa,double bb,double cc, double dd, double a, double b, complex double roots[2])
{ 
  // double MACHEPS=2.2204460492503131E-16;
  double sqrtd, diskr,div,zmax,zmin;
  int iter; 
#ifdef LDLT_USENR
  int k;
  double tt[9], ttmax;
#endif
  //-------------------------------- parameter backward correction step:       
#if 0
  /* newtown-raphson must be improved by using bisection if needed (see Numerical Recipat algo NR safe with bisection) */
  evec[0]=bb*b-b*b-a*b*aa+a*a*b-dd;            // equation (5.18)
  evec[1]=cc*b-b*b*aa+b*b*a-a*dd;              // equation (5.18)
  err=fabs(evec[0])+fabs(evec[1]);                 // equation (5.23)
  //printf("BEGIN err=%.16G a=%.16G b=%.15G\n", err, a, b);

  tt[0]=fabs(bb*b);
  tt[1]=fabs(b*b);
  tt[2]=fabs(a*b*aa);
  tt[3]=fabs(a*a*b);
  tt[4]=fabs(dd);
  tt[5]=fabs(cc*b);
  tt[6]=fabs(b*b*aa);
  tt[7]=fabs(b*b*a);
  tt[8]=fabs(a*dd);
  for (k=0; k < 9; k++)
    if (k==0 || tt[k] > ttmax)
      ttmax = tt[k];
  if (err > 2.0*MACHEPS*ttmax)    
    if (err!=0)
    {      
      //errold=err;
      //printf("PRIMA err=%.16G\n", err);
      //for (iter=0; iter < 8; iter++) 
      for (iter=0; iter < 8; iter++) 
	{ 
#if 0
	  if (err <= MACHEPS*(fabs(bb*b)+fabs(b*b)+fabs(a*b*aa)+fabs(a*a*b)+fabs(dd)+fabs(cc*b)+fabs(b*b*aa)+
			      fabs(b*b*a)+fabs(a*dd)))
	    break;
#endif
	  fmat[0][0]=-b*aa+2*a*b;                         // equation (5.19)
	  fmat[0][1]=bb-2*b-a*aa+a*a;                    // equation (5.19)
	  fmat[1][0]=b*b-dd;                             // equation (5.19)
	  fmat[1][1]=cc-2*b*aa+2*b*a;                     // equation (5.19)
	  /* usuale criterio di convergenza di un NR */
	  //if (err <= MACHEPS)
	    //break;
	  evec[0]=-evec[0];
	  evec[1]=-evec[1];

	  two_lin_eqs(fmat,evec,dpar);              // equation (5.20)
	  at=a;
	  bt=b;
	  a=a+dpar[0];                                   // equation (5.21)
	  b=b+dpar[1];                                   // equation (5.22)
	  //printf("iter=%d a=%.15G b=%.16G evec=%.16G %.16G\n", iter, a, b, evec[0], evec[1]);
  	  evec[0]=bb*b-b*b-a*b*aa+a*a*b-dd;            // equation (5.18)
	  evec[1]=cc*b-b*b*aa+b*b*a-a*dd;              // equation (5.18)
	  errt=err;
	  err=fabs(evec[0])+fabs(evec[1]);                 // equation (5.23) 
  	  if(err == 0.0)
	    break;	// terminate
	  if(err>=errt)   // terminate without parameter update
	    {
	      //printf("TOO LARGE QUI\n");
	      a=at;
	      b=bt;  
	      break;
	    }
  	  //if (fabs(dpar[0]) + fabs(dpar[1]) <= MACHEPS*(fabs(a)+fabs(b)))
	    //break;
	}
  //printf("BEGIN err=%.16G a=%.16G b=%.15G\n", err, a, b);
      //if (iter > 8)
	//printf("FINAL err=%.16G a=%.16G b=%.15G errt=%.15G iter=%d\n", err, a, b, errt, iter);
      
      //if (err < errold)
      //printf("DOPO err=%.16G\n", err);
    } 
  //---------------------------------------- solve a quadratic equation:     
  //if (iter>4 && err > 3E-16)
   // printf("iter=%d max iterations reached err=%.16G!\n", iter, err);
#endif

  diskr=a*a-4*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
	zmin=0.0;
      else
	zmin=b/zmax;

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=CMPLX(-a/2,sqrtd/2);
      roots[1]=CMPLX(-a/2,-sqrtd/2);      
    }   
}
#ifdef LDLT_OPT_CUBIC
void solve_cubic_analytic_depressed(double b, double c, double *sol);
void solve_cubic_analytic_depressed_handle_inf(double b, double c, double *sol);

#else
void solve_cubic_analytic_depressed(double b, double c, complex double sol[3]);
void solve_cubic_analytic_depressed_handle_inf(double b, double c, complex double sol[3]);
#endif
void  mycubic_B_shift(double a, double b, double c, double d, double *phi0)
{
  double g,h,gg,hh,aq,bq,cq,dq,s,diskr, coeff[4], err;
  int nreal,k, trovato=0, k2;
  double rmax;
  complex double radici[3];
  complex long double radicil[3];
  long double coeffl[4];
#ifdef LDLT_REFINE_PHI_WITH_NR
  double MACHEPS=2.2204460492503131E-16;
  double maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  int iter;
#endif
  double kd0, kd1, kd2; 
  //------------------ call the parabola/reciprocal intersection method:      
  // !------------------------------------------------------- the B-shift:
#if 1
  diskr=9*a*a-24*b;                    //         ! equation (3.58)

  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
#if 0
      /* b + 3 a s +  6 s^2 = 0*/ 	
      /* if we use standard formula for solving a quadratic
       * the most accurate root is the largest one */
      if (a > 0)
	{
	  s = (-3*a - diskr)/12.0; 
	}
      else
	{
	  s = (-3*a + diskr)/12.0;
	}
#else

      if(a > 0.0)
	s=-2*b/(3*a+diskr);                            // equation (3.58)
      else
	s=-2*b/(3*a-diskr);                            // equation (3.58)
#endif
    }
  else
    {      
      s=-a/4;                                    // equations (3.59-60)
    }
  // !--------------------------- the shift transformation (Horner forms):
#if 0
    {
      double sold;
      sold=s;
      s=0; 
      aq=a+4*s;                                      // equation (3.45)
      bq=b+3*s*(a+2*s);                              // equation (3.46)
      cq=c+s*(2*b+s*(3*a+4*s));                      // equation (3.47)
      dq=d+s*(c+s*(b+s*(a+s)));                      // equation (3.48)     
      printf("[s=0] (8*dq+hh-2*gg)*bq/3=%.16G cq*cq=%16G dq*aq*aq=%.16G\n", (8*dq+hh-2*gg)*bq/3, cq*cq, dq*aq*aq);
      s=sold;
    }
#endif
  //s=0;
#if 0
  aq=a+4*s;                                      // equation (3.45)
  bq =b + 3*a*s + 6*s*s;
  cq = c + 2*b*s + 3*a*s*s + 4*s*s*s; 
  dq = d + c*s + b*s*s + a*s*s*s + s*s*s*s;
#else
  aq=a+4*s;                                      // equation (3.45)
  bq=b+3*s*(a+2*s);                              // equation (3.46)
  cq=c+s*(2*b+s*(3*a+4*s));                      // equation (3.47)
  dq=d+s*(c+s*(b+s*(a+s)));                      // equation (3.48)     
#endif
  //printf("a=%.15G ap=%.15G b=%.15G bp=%.15G c=%.15G cp=%.15G d=%.15G dp=%.15G\n", a, aq, b, bq, c, cq, d, dq);
  gg=bq*bq/9;
  hh=aq*cq;      
  g=hh-4*dq-3*gg;                                // equation (3.49)
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq;         // equation (3.50)
  //printf("g=%.15G h=%.15G\n", g, h);
  //printf("(8*dq+hh-2*gg)*bq/3=%.16G cq*cq=%16G dq*aq*aq=%.16G\n", (8*dq+hh-2*gg)*bq/3, cq*cq, dq*aq*aq);
#if 0
  if (1)
    {
      kd1=fabs(h);
      kd0 = cbrt(kd1);
      kd2= kd0*kd0;
      gg=(bq/kd1)*bq/9;
      hh=aq*(cq/kd1);      
      h = (8*dq/kd1+hh-2*gg)*bq/3-cq*(cq/kd1)-aq*aq*(dq/kd1);
      //h=1.0;
      hh=aq*(cq/kd2);
      gg=(bq/kd2)*bq/9.0;
      g =hh-4*dq/kd2-3*gg;
      //printf("h=%.15G g=%.15G\n", g, h);
    }
#endif
#if 0
  coeff[3]=1.0;
  coeff[2]=0.0;
  coeff[1]=g;
  coeff[0]=h;
#endif
#else
  coeff[3]=1.0;
  coeff[2]=3.0*aq/4.0;
  coeff[1]= bq/2.0;
  coeff[0]= cq/4.0;
  solve_cubic_analytic(coeff,radici);
  s = creal(radici[0]);
  /* provo ad azzerare C anziché B comunque non sembra dare risultati migliori rispetto ad azzerare B
   * ed è un po' più lento  */
  aq=a+4*s;                                      // equation (3.45)
  bq=b+3*s*(a+2*s);                              // equation (3.46)
  cq=c+s*(2*b+s*(3*a+4*s));                      // equation (3.47)
  dq=d+s*(c+s*(b+s*(a+s)));                      // equation (3.48)     
  gg=bq*bq/9;
  hh=aq*cq; 
  g=hh-4*dq-3*gg;                                // equation (3.49)
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq;         // equation (3.50)

  coeff[3]=1.0;
  coeff[2]=0.0;
  coeff[1]=g;
  coeff[0]=h; 
#endif
  /* NOTA 09/02/18 La soluzione analitica sembra accurata quanto il metodo usato da Strobach ma è più veloce,
   * In alternativa si può usare l'algoritmo di Flocke per le cubiche */
  //cubicRoots(0,g,h,&nreal,radici);
  //solve_cubic_analytic(coeff,radici);
#if 1
  //cubicRoots(0,g,h,&nreal,radici);
#ifdef LDLT_OPT_CUBIC
#ifdef USE_CUBIC_HANDLE_INF_ONLY
  solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
#else
  solve_cubic_analytic_depressed(g, h, &rmax);
  for (k=0; k < 3; k++)
    {
      if (isnan(rmax) || isinf(rmax))
	{
	  solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
	  //cubicRoots(0,g,h,&nreal,radici);
	  break;
	}
    }
#endif
#else
#ifdef USE_CUBIC_HANDLE_INF_ONLY
  solve_cubic_analytic_depressed_handle_inf(g, h, radici);
#else
  solve_cubic_analytic_depressed(g, h, radici);
#if 1
  for (k=0; k < 3; k++)
    {
      if (isnan(creal(radici[k])) || isnan(cimag(radici[k])) 
	  || isinf(creal(radici[k])) || isinf(cimag(radici[k])) )
	{
	  //printf("QUIQUIIIIIIIIIIIIIIIIIIIIII\n");
#if 0
	  for (k2=0; k2 < 4; k2++)
	    coeffl[k2] = coeff[k2];
	  solve_numrecl(coeffl, radicil, 3);
	  for (k2=0; k2 < 3; k2++)
	    radici[k2] = radicil[k2];
#endif
	  solve_cubic_analytic_depressed_handle_inf(g, h, radici);
	  //cubicRoots(0,g,h,&nreal,radici);
	  break;
	}
    }
#endif
#endif
  for (k=0; k < 3; k++)
    {
      //printf("rad[%d]=%.15G+I*%.15G\n", k, creal(radici[k]), cimag(radici[k]));
      if (cimag(radici[k])==0)
	{
	  if (!trovato || fabs(creal(radici[k])) > fabs(rmax) )
	    {
	      trovato=1;
	      rmax=creal(radici[k]);
	    }
	}
    }
#endif
#ifdef LDLT_REFINE_PHI_WITH_NR
  x = rmax;
  xsq=x*x;
  xxx=x*xsq;
  gx=g*x;
  f = x*(xsq + g) + h;
  if (fabs(xxx) > fabs(gx))
    maxtt = fabs(xxx);
  else
    maxtt = fabs(gx);
  if (fabs(h) > maxtt)
    maxtt = fabs(h);

  if (fabs(f) > MACHEPS*maxtt)
    {
      for (iter=0; iter < 8; iter++)
	{   
	  df =  3.0*xsq + g;
	  if (df==0)
	    {
	      break;
	    }
	  xold = x;
	  x += -f/df;
	  fold = f;
	  xsq = x*x;
	  f = x*(xsq + g) + h;
	  //f = x*x*x + g*x + h;
	  if (f==0)
	    {
	      break;
	    } 
	  //printf("iter %d f=%.25G\n", iter, fabs(f));
       	  if (fabs(f) >= fabs(fold))
	    {
	      x = xold;
	      break;
	    }
    	}
    }
  *phi0 = x;
#else
  *phi0 = rmax;
#endif
  //printf("phi=%.1G\n", *phi0);
#else
  depressed_cubic_root(g,h,phi0);
#endif


  /* il Newton raphson sembra che non serva a niente */
  // printf("PRIMA phi0=%.16G\n", *phi0);
  //cubicNR(phi0, g, h);

  //printf("err3=%.20G err2=%.15G phi0=%.15G\n", targetfunc(4, coeff, *phi0), err, *phi0);
  // -------------------------------------------------------------------- 
}
double calc_err_ldlt(double a, double b, double c, double d, double d2, double l1, double l2, double l3)
{
  /* d3 + d2 l2^2 + d1 l3^2 + (2 d2 l2 + 2 d1 l1 l3) z + (d2 + d1 l1^2 + 2 d1 l3) z^2 + 
     2 d1 l1 z^3 + d1 z^4 */
  /* d3 + d2 l2^2 + d1 l3^2 + (2 d2 l2 + 2 d1 l1 l3) z + (d2 + d1 l1^2 + 2 d1 l3) z^2 + 
     2 d1 l1 z^3 + d1 z^4 */
  double sum;
#if 0
  return max3((b==0)?fabs(d2 + l1*l1 + 2.0*l3):fabs(((d2 + l1*l1 + 2.0*l3)-b)/b),
              (c==0)?fabs(2.0*d2*l2 + 2.0*l1*l3):fabs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c),
             (d==0)?fabs(d2*l2*l2 + l3*l3):fabs(((d2*l2*l2 + l3*l3)-d)/d));
#else
  //return fabs((d2 + l1*l1 + 2.0*l3)-b)+fabs((2.0*d2*l2 + 2.0*l1*l3)-c)+fabs((d2*l2*l2 + l3*l3)-d);
  sum =  (b==0)?fabs(d2 + l1*l1 + 2.0*l3):fabs(((d2 + l1*l1 + 2.0*l3)-b)/b);
  sum += (c==0)?fabs(2.0*d2*l2 + 2.0*l1*l3):fabs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  sum += (d==0)?fabs(d2*l2*l2 + l3*l3):fabs(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
#endif
}
double calc_err_abc(double a, double b, double c, double aq, double bq, double cq, double dq)
{
  double sum;

  //return fabs((bq*cq + aq*dq) - c) + fabs((bq + aq*cq + dq) - b) + fabs((aq + cq) - a);
#if 1
  sum = (c==0)?fabs(bq*cq + aq*dq):fabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabs(bq + aq*cq + dq):fabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabs(aq + cq):fabs(((aq + cq) - a)/a);
  //printf("abc SUM=%.15G\n", sum);
  return sum;
#endif
}
#if 1
//#define LU_DECOMP
void qludcmp(double a[4][4], int indx[4], double* d, int *ok)
{
  /* A[i][j] = Aij 
   * A x = b  
   * per semplicità nel seguito si assume che l'ordine della matrice è 3 */
  int i,imax=-1,j,k;
  const int n=4;
  const double TINY=1E-40;
  double big,dum,sum,temp; 
  double vv[8]; /* vv stores the implicit scaling of each row.*/
  /*vv = vector(1,n);*/
  *d=1.0; /* No row interchanges yet. */
  *ok = 0;
  for (i=0;i<n;i++) 
    { 
      /* Loop over rows to get the implicit scaling information.*/ 
      big=0.0; 
      for (j=0;j<n;j++)
	{
	  if ((temp=fabs(a[i][j])) > big) big=temp; 
	}
      if (big == 0.0)
	{
	  //ERROR: Singular matrix in routine ludcmp 
	  *ok = 1;
	  return;
	}
      /* No nonzero largest element. */
      vv[i]=1.0/big; /* Save the scaling.*/
    } 
  for (j=0;j<n;j++) 
    { /* This is the loop over columns of Crout s method.*/
      for (i=0;i<j;i++) 
	{ 
	  /* This is equation (2.3.12) except for i = j. */
	  sum=a[i][j]; 
	  for (k=0;k<i;k++) 
	    sum -= a[i][k]*a[k][j]; 
	  a[i][j]=sum; 
	} 
      big=0.0; /* Initialize for the search for largest pivot element. */ 
      for (i=j;i<n;i++) 
	{ 
	  /* This is i = j of equation (2.3.12) and i = j+1. . .N of equation (2.3.13).*/
	  sum=a[i][j]; 
	  for (k=0;k<j;k++)
	    sum -= a[i][k]*a[k][j]; 
	    a[i][j]=sum; 
	    if ( (dum=vv[i]*fabs(sum)) >= big) 
	      { 
		/* Is the  gure of merit for the pivot better than the best so far? */
		big=dum; imax=i; 
	      } 
	} 
      if (j != imax) 
	{ 
	  /* Do we need to interchange rows? */
	  for (k=0;k<n;k++) 
	    { 
	      /* Yes, do so...*/ 
	      dum=a[imax][k]; 
	      a[imax][k]=a[j][k]; 
	      a[j][k]=dum; 
	    } 
	  *d = -(*d); 
	  /* ...and change the parity of d. */ 
	  vv[imax]=vv[j]; 
	  /* Also interchange the scale factor.*/ 
	} 
      indx[j]=imax; 
      if (a[j][j] == 0.0) 
	a[j][j]=TINY; 
      /* If the pivot element is zero the matrix is singular 
       * (at least to the precision of the algorithm). 
       * For some applications on singular matrices, 
       * it is desirable to substitute TINY for zero. */ 
      if (j != n) 
	{ 
	  /* Now,  nally, divide by the pivot element.*/
	  dum=1.0/(a[j][j]); 
	  for (i=j+1;i<n;i++) a[i][j] *= dum; 
	} 
    } 
  /* Go back for the next column in the reduction.*/
  /*free_vector(vv,1,n); */
}

void qlubksb(double a[4][4], int indx[4], double b[4])
{ 
  int i,ii=0,ip,j; 
  double sum; 
  const int n=4;
  for (i=0;i<n;i++) 
    { 
      /* When ii is set to a positive value, it will become the index of the  
       * rst nonvanishing element of b. Wenow do the forward substitution,
       * equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go. */
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i]; 
      if (ii>-1) 
	for (j=ii;j<=i-1;j++) 
	  sum -= a[i][j]*b[j]; 
      else if (sum) 
	ii=i; 
      /* A nonzero element was encountered, so from now on we will have to do 
       * the sums in the loop above. */ 
      b[i]=sum; 
    } 
  for (i=n-1;i>=0;i--) 
    { 
      /* Now we do the backsubstitution, equation (2.3.7).*/
      sum=b[i]; 
      for (j=i+1;j<n;j++) 
	sum -= a[i][j]*b[j]; b[i]=sum/a[i][i]; 
      /* Store a component of the solution vector X. */ 
    } /* All done! */
}
void NRabcdCmplx(double a, double b, double c, double d, 
		 double complex *AQ, double complex *BQ, double complex *CQ, double complex *DQ)
{
  int iter, k1, k2;
  complex double x02, xold[4], x[4], dx[4], det, aq, bq, cq, dq, Jinv[4][4], fvec[4];
  double errf, errfold;
#ifdef LU_DECOMP
  complex double fjac[4][4], p[4];
  int indx[4], ok;
  complex double dlu;
#endif
  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
  /* risolvo con un NR il sistema:
   *
   * bq*dq - d = 0
   * bq*cq + aq*dq - c = 0
   * bq + aq*cq + dq - b = 0
   * aq + cq - a = 0
   *
   * che si ottiene moltiplicando i due polinomi quadratici
   *
   * (x^2 + aq*x + bq)*(x^2 + cq*x + dq)
   *
   * ed equagualiando i coefficienti che si ottengono con quelli della quartica (cioè a,b,c e d)
   * */

  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  for (k1=0; k1 < 4; k1++)
    {
      errf += cabs(fvec[k1]);
    }
  if (errf==0)
    return;
  //printf("abcd=%f %f %f %f\n", x[0], x[1], x[2], x[3]);
  for (iter = 0; iter < 8; iter++)
    {
#ifdef LU_DECOMP
      fjac[0][0]=0;
      fjac[0][1]=x[3];
      fjac[0][2]=0;
      fjac[0][3]=x[1];
      fjac[1][0]=x[3];
      fjac[1][1]=x[2];
      fjac[1][2]=x[1];
      fjac[1][3]=x[0];
      fjac[2][0]=x[2];
      fjac[2][1]=1;
      fjac[2][2]=x[0];
      fjac[2][3]=1;
      fjac[3][0]=1;
      fjac[3][1]=0;
      fjac[3][2]=1;
      fjac[3][3]=0;	
      qludcmp(fjac,indx,&dlu, &ok); /* Solve linear equations by LU decomposition.*/
      if (ok != 0)
	break;
      for (k1=0; k1 < 4; k1++)
	dx[k1] = fvec[k1];
      qlubksb(fjac,indx, dx);
      //printf("iter=%d dx=%.16G %.16G %.16G %.16G\n", iter, dx[0], dx[1], dx[2], dx[3]);
      for (k1=0; k1 < 4; k1++)
      	xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1];
	}
#else
#if 1
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - 2.0*x[3]) + x[3]*(x[0]*x02 + x[3]);
#else
      det = x[1]*x[1] + x[1]*(-x[0]*x[2] + x[2]*x[2] - 2.0*x[3]) + x[3]*(x[0]*x[0] - x[0]*x[2] + x[3]);
#endif
      if (det==0.0)
	break;
#if 1
      Jinv[0][0] = x02;
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = -x[1]*Jinv[0][1] - x[0]*Jinv[0][2]; 
      Jinv[1][0] = x[0]*Jinv[0][0] + Jinv[0][1];
      Jinv[1][1] = -x[1]*Jinv[0][0];
      Jinv[1][2] = -x[1]*Jinv[0][1];   
      Jinv[1][3] = -x[1]*Jinv[0][2];
      Jinv[2][0] = -Jinv[0][0];
      Jinv[2][1] = -Jinv[0][1];
      Jinv[2][2] = -Jinv[0][2];
      Jinv[2][3] = Jinv[0][2]*x[2] + Jinv[0][1]*x[3];
      Jinv[3][0] = -x[2]*Jinv[0][0] - Jinv[0][1];
      Jinv[3][1] = Jinv[0][0]*x[3];
      Jinv[3][2] = x[3]*Jinv[0][1];
      Jinv[3][3] = x[3]*Jinv[0][2];
#else
      Jinv[0][0] = x[0] - x[2];
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = x[1]*x[1] + x[0]*x[0]*x[3] - x[1]*(x[0]*x[2] + x[3]); 
      Jinv[1][0] = x[0]*x[0] - x[1] - x[0]*x[2] + x[3];
      Jinv[1][1] = x[1]*(-x[0] + x[2]);
      Jinv[1][2] = x[1]*(x[1] - x[3]);   
      Jinv[1][3] = x[1]*(-x[1]*x[2] + x[0]*x[3]);
      Jinv[2][0] =-x[0] + x[2];
      Jinv[2][1] = x[1] - x[3];
      Jinv[2][2] = -x[1]*x[2] + x[0]*x[3];
      Jinv[2][3] = x[1]*x[2]*x[2] - (x[1] + x[0]*x[2])*x[3] + x[3]*x[3];
      Jinv[3][0] = x[1] - x[0]*x[2] + x[2]*x[2] - x[3];
      Jinv[3][1] = (x[0] - x[2])*x[3];
      Jinv[3][2] = x[3]*(-x[1] + x[3]);
      Jinv[3][3] = x[3]*(x[1]*x[2] - x[0]*x[3]);
#endif
      for (k1=0; k1 < 4; k1++)
	{
	  dx[k1] = 0;
	  for (k2=0; k2 < 4; k2++)
	    dx[k1] += Jinv[k1][k2]*fvec[k2];
	}
      for (k1=0; k1 < 4; k1++)
      	xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1]/det;
	}

#endif
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 

      errfold = errf;
      errf=0;
      for (k1=0; k1 < 4; k1++)
	{
	  errf += cabs(fvec[k1]);
	}
      
      //printf("iter=%d errf=%.15G\n", iter, errf);
      //printf("iter=%d abcd=%f %f %f %f errf=%.16G errfold=%.16G\n", iter, x[0], x[1], x[2], x[3], errf, errfold);
      if (errf==0)
	break;
      if (errf >= errfold)
	{
	  for (k1=0; k1 < 4; k1++)
	    x[k1] = xold[k1];
	  break;
	}
    }

  *AQ=x[0];
  *BQ=x[1];
  *CQ=x[2];
  *DQ=x[3];
}
void NRabcd(double a, double b, double c, double d, double *AQ, double *BQ, double *CQ, double *DQ)
{
  int iter, k1, k2;
  double x02, errf, errfold, xold[4], x[4], dx[4], det, aq, bq, cq, dq, Jinv[4][4], fvec[4];
#ifdef LU_DECOMP
  double fjac[4][4], p[4];
  int indx[4], ok;
  double dlu;
#endif
  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
  /* risolvo con un NR il sistema:
   *
   * bq*dq - d = 0
   * bq*cq + aq*dq - c = 0
   * bq + aq*cq + dq - b = 0
   * aq + cq - a = 0
   *
   * che si ottiene moltiplicando i due polinomi quadratici
   *
   * (x^2 + aq*x + bq)*(x^2 + cq*x + dq)
   *
   * ed equagualiando i coefficienti che si ottengono con quelli della quartica (cioè a,b,c e d)
   * */

  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  for (k1=0; k1 < 4; k1++)
    {
      errf += fabs(fvec[k1]);
    }
  if (errf==0)
    return;
  //printf("abcd=%f %f %f %f\n", x[0], x[1], x[2], x[3]);
  for (iter = 0; iter < 8; iter++)
    {
#ifdef LU_DECOMP
      fjac[0][0]=0;
      fjac[0][1]=x[3];
      fjac[0][2]=0;
      fjac[0][3]=x[1];
      fjac[1][0]=x[3];
      fjac[1][1]=x[2];
      fjac[1][2]=x[1];
      fjac[1][3]=x[0];
      fjac[2][0]=x[2];
      fjac[2][1]=1;
      fjac[2][2]=x[0];
      fjac[2][3]=1;
      fjac[3][0]=1;
      fjac[3][1]=0;
      fjac[3][2]=1;
      fjac[3][3]=0;	
      qludcmp(fjac,indx,&dlu, &ok); /* Solve linear equations by LU decomposition.*/
      if (ok != 0)
	break;
      for (k1=0; k1 < 4; k1++)
	dx[k1] = fvec[k1];
      qlubksb(fjac,indx, dx);
      //printf("iter=%d dx=%.16G %.16G %.16G %.16G\n", iter, dx[0], dx[1], dx[2], dx[3]);
      for (k1=0; k1 < 4; k1++)
      	xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1];
	}
#else
#if 1
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - 2.0*x[3]) + x[3]*(x[0]*x02 + x[3]);
#else
      det = x[1]*x[1] + x[1]*(-x[0]*x[2] + x[2]*x[2] - 2.0*x[3]) + x[3]*(x[0]*x[0] - x[0]*x[2] + x[3]);
#endif
      if (det==0.0)
	break;
#if 1
      Jinv[0][0] = x02;
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = -x[1]*Jinv[0][1] - x[0]*Jinv[0][2]; 
      Jinv[1][0] = x[0]*Jinv[0][0] + Jinv[0][1];
      Jinv[1][1] = -x[1]*Jinv[0][0];
      Jinv[1][2] = -x[1]*Jinv[0][1];   
      Jinv[1][3] = -x[1]*Jinv[0][2];
      Jinv[2][0] = -Jinv[0][0];
      Jinv[2][1] = -Jinv[0][1];
      Jinv[2][2] = -Jinv[0][2];
      Jinv[2][3] = Jinv[0][2]*x[2] + Jinv[0][1]*x[3];
      Jinv[3][0] = -x[2]*Jinv[0][0] - Jinv[0][1];
      Jinv[3][1] = Jinv[0][0]*x[3];
      Jinv[3][2] = x[3]*Jinv[0][1];
      Jinv[3][3] = x[3]*Jinv[0][2];
#else
      Jinv[0][0] = x[0] - x[2];
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = x[1]*x[1] + x[0]*x[0]*x[3] - x[1]*(x[0]*x[2] + x[3]); 
      Jinv[1][0] = x[0]*x[0] - x[1] - x[0]*x[2] + x[3];
      Jinv[1][1] = x[1]*(-x[0] + x[2]);
      Jinv[1][2] = x[1]*(x[1] - x[3]);   
      Jinv[1][3] = x[1]*(-x[1]*x[2] + x[0]*x[3]);
      Jinv[2][0] =-x[0] + x[2];
      Jinv[2][1] = x[1] - x[3];
      Jinv[2][2] = -x[1]*x[2] + x[0]*x[3];
      Jinv[2][3] = x[1]*x[2]*x[2] - (x[1] + x[0]*x[2])*x[3] + x[3]*x[3];
      Jinv[3][0] = x[1] - x[0]*x[2] + x[2]*x[2] - x[3];
      Jinv[3][1] = (x[0] - x[2])*x[3];
      Jinv[3][2] = x[3]*(-x[1] + x[3]);
      Jinv[3][3] = x[3]*(x[1]*x[2] - x[0]*x[3]);
#endif
      for (k1=0; k1 < 4; k1++)
	{
	  dx[k1] = 0;
	  for (k2=0; k2 < 4; k2++)
	    dx[k1] += Jinv[k1][k2]*fvec[k2];
	}
      for (k1=0; k1 < 4; k1++)
      	xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1]/det;
	}

#endif
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 

      errfold = errf;
      errf=0;
      for (k1=0; k1 < 4; k1++)
	{
	  errf += fabs(fvec[k1]);
	}
      //printf("iter=%d abcd=%f %f %f %f errf=%.16G errfold=%.16G\n", iter, x[0], x[1], x[2], x[3], errf, errfold);
      if (errf==0)
	break;
      if (errf >= errfold)
	{
	  for (k1=0; k1 < 4; k1++)
	    x[k1] = xold[k1];
	  break;
	}
    }

  *AQ=x[0];
  *BQ=x[1];
  *CQ=x[2];
  *DQ=x[3];
}
#else
void NRabcd(double a, double b, double c, double d, double *AQ, double *BQ, double *CQ, double *DQ)
{
  int iter, k1, k2;
  double errf, errfold, xold[4], x[4], dx[4], det, aq, bq, cq, dq, Jinv[4][4], fvec[4];
  double  x00, x11, x01, x02, x12, x22, x03, x13, x23, x33; 

  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
  /* risolvo con un NR il sistema:
   *
   * bq*dq - d = 0
   * bq*cq + aq*dq - c = 0
   * bq + aq*cq + dq - b = 0
   * aq + cq - a = 0
   *
   * che si ottiene moltiplicando i due polinomi quadratici
   *
   * (x^2 + aq*x + bq)*(x^2 + cq*x + dq)
   *
   * ed equagualiando i coefficienti che si ottengono con quelli della quartica (cioè a,b,c e d)
   * */

  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  for (k1=0; k1 < 4; k1++)
    {
      errf += fabs(fvec[k1]);
    }
  if (errf==0)
    return;
  //printf("abcd=%f %f %f %f\n", x[0], x[1], x[2], x[3]);
  for (iter = 0; iter < 8; iter++)
    {
      x00 = x[0]*x[0];
      x11 = x[1]*x[1];
      x12 = x[1]*x[2];
      x01 = x[0]*x[1];
      x02 = x[0]*x[2];
      x22 = x[2]*x[2];
      x03 = x[0]*x[3];
      x13 = x[1]*x[3];
      x23 = x[2]*x[3];
      x33 = x[3]*x[3];

      det = x11 + x[1]*(-x02 + x22 - 2.0*x[3]) + x[3]*(x00 - x02 + x[3]);
      if (det==0.0)
	break;
      Jinv[0][0] = x[0] - x[2];
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x12 - x03;
      Jinv[0][3] = x11 + x[0]*x03 - x[1]*x02 - x13; 
      Jinv[1][0] = x00 - x[1] - x02 + x[3];
      Jinv[1][1] = -x01 + x12;// x[1]*(-x[0] + x[2]);
      Jinv[1][2] =  x11 - x13; //x[1]*(x[1] - x[3]);   
      Jinv[1][3] = x[1]*x12 + x[1]*x03;//x[1]*(-x[1]*x[2] + x[0]*x[3]);
      Jinv[2][0] =-x[0] + x[2];
      Jinv[2][1] = x[1] - x[3];
      Jinv[2][2] = -x12 + x03;
      Jinv[2][3] = x[1]*x22 - x13 - x02*x[3] + x33;//x[1]*x22 - (x[1] + x02)*x[3] + x33;
      Jinv[3][0] = x[1] - x02 + x22 - x[3];
      Jinv[3][1] = x03 - x23;//(x[0] - x[2])*x[3];
      Jinv[3][2] = -x13 + x33;//x[3]*(-x[1] + x[3]);
      Jinv[3][3] = x[3]*(x12 - x03);

      for (k1=0; k1 < 4; k1++)
	{
	  dx[k1] = 0;
	  for (k2=0; k2 < 4; k2++)
	    dx[k1] += Jinv[k1][k2]*fvec[k2];
	}
     for (k1=0; k1 < 4; k1++)
       xold[k1] = x[k1];
 
      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1]/det;
	}
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 

      errfold = errf;
      errf=0;
      for (k1=0; k1 < 4; k1++)
	{
	  errf += fabs(fvec[k1]);
	}
      if (errf==0)
	break;
      //printf("iter=%d abcd=%f %f %f %f errf=%.16G errfold=%.16G\n", iter, x[0], x[1], x[2], x[3], errf, errfold);
      if (errf >= errfold)
	{
	  for (k1=0; k1 < 4; k1++)
	    x[k1] = xold[k1];
	  break;
	}
    }

  *AQ=x[0];
  *BQ=x[1];
  *CQ=x[2];
  *DQ=x[3];
}
#endif
double calc_err_abcd_cmplx(double a, double b, double c, double d, 
			   complex double aq, complex double bq, complex double cq, complex double dq)
{
  double sum;
  sum = (d==0)?cabs(bq*dq):cabs((bq*dq-d)/d);
  sum += (c==0)?cabs(bq*cq + aq*dq):cabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?cabs(bq + aq*cq + dq):cabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?cabs(aq + cq):cabs(((aq + cq) - a)/a);
  return sum;
}
double calc_err_abcd(double a, double b, double c, double d, double aq, double bq, double cq, double dq)
{
  double sum;
  sum = (d==0)?fabs(bq*dq):fabs((bq*dq-d)/d);
  sum += (c==0)?fabs(bq*cq + aq*dq):fabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabs(bq + aq*cq + dq):fabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabs(aq + cq):fabs(((aq + cq) - a)/a);
  return sum;
}

void LDLT_quartic(double coeff[5], complex double roots[4])      
{
  int realcase[2], whichcase;
  double err0, err1;
  double aq1, bq1, cq1, dq1; 
  complex double acx1, bcx1, ccx1, dcx1;
  double l2m[12], d2m[12], res[12], resmin;
  double d2eq46, bl311, dml3l3;
  int k1, k, kmin, nsol; 
  double a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3];
  double gamma,del2;
  const double macheps =2.2204460492503131E-16;
  double cbq[3];
  complex double sbq[2]; 
  double cubc[4];
  double complex acx,bcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
  double ssd;
  complex long double rri, rmri;
  complex double ccx, dcx;
  //int d2z=0;
  //----------------------------- calculate the antidiagonal shift phi0:

  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
#if 1
  /* special cases to handle */
  if (a==0 && b==0 && c==0 && d==0)
    {
      roots[0]=roots[1]=roots[2]=roots[3]=0;
      return;
    }
  else if (a==0 && b==0 && c==0)
    {
      if (d < 0.0)
	{
	  ssd = pow(-d,0.25);
	  roots[0]=ssd+I*0.0;
	  roots[1]=-ssd+I*0.0;
	  roots[2]=0.0+I*ssd;
	  roots[3]=0.0-I*ssd;
	}
      else
	{
	  //printf("qui\n");
	  ssd = pow(d,0.25);
	  rri = cosl(M_PI/4.0)+sinl(M_PI/4.0)*I;
	  rmri=I*(cosl(M_PI/4.0)+sinl(M_PI/4.0)*I);
	  roots[0]=rri*ssd;
	  roots[1]=-rri*ssd;
	  roots[2]=rmri*ssd;
	  roots[3]=-rmri*ssd;
	}
      return;
    }
  else if (a==0 && b==0 && d==0) 
    {
      cubc[3]=1.0;
      cubc[2]=0.0;
      cubc[1]=0.0;
      cubc[0]=c;
      //printf("QUIIIIIIIIIII\n");
      //csolve_cubic(cubc, roots);
      solve_cubic_analytic(cubc, roots);
      //cubicRoots(0, 0, c, &nreal, roots);
      roots[3]=0.0+I*0.0;
      return;
    }
  else if (b==0 && c==0 && d==0)
    {
      roots[0]=roots[1]=roots[2]=0.0+I*0.0;
      roots[3]=-a+I*0.0;
      return;
    }
#else
  if ((a==0 && b==0 && c==0 && d==0)||
      (a==0 && b==0 && c==0)||
      (a==0 && b==0 && d==0)|| 
      (b==0 && c==0 && d==0))
    {
      /* temporary fix */
      quarticRoots(coeff, &nreal, roots);
      return;
    }
#endif
  mycubic_B_shift(a,b,c,d,&phi0);
  l1=a/2;                                         // equation (4.2)
  l3=b/6+phi0/2;                                  // equation (4.3)
  del2=c-a*l3;                                    // equation (4.10) 
  nsol=0;

  d2eq46 =2.*b/3.-phi0-l1*l1;
  dml3l3 = d-l3*l3;
  //printf("dml3l3=%.15G\n", dml3l3);
  bl311 = d2eq46;//b-2*l3-l1*l1;
  if (bl311!=0.0)
    {
      d2m[nsol] = d2eq46;  // eq. (4.6)
      l2m[nsol] = del2/(2.0*d2m[nsol]);   // eq. (4.12)
      //res[nsol] = fabs(del2*l2m[nsol]-2.0*dml3l3);
      //res[nsol] = fabs(d2m[nsol]*l2m[nsol]*l2m[nsol]-dml3l3);
      res[nsol] = calc_err_ldlt(a,b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }
  if (del2!=0)
    {
      l2m[nsol]=2*dml3l3/del2;
      if (l2m[nsol]!=0)
	{
  	  d2m[nsol]=del2/(2*l2m[nsol]);
	  //res[nsol] = fabs(d2m[nsol]-bl311); // nel calcolo della soluzione non uso la (4.6)
	  res[nsol] = calc_err_ldlt(a,b,c,d,d2m[nsol], l1, l2m[nsol], l3);
	  nsol++;
	}

      d2m[nsol] = d2eq46;
      l2m[nsol] = 2.0*dml3l3/del2;
      //res[nsol] = fabs(2.0*d2m[nsol]*l2m[nsol] - del2); 
      res[nsol] = calc_err_ldlt(a,b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }

  if (nsol==0)
    {
      l2=d2=0.0;
    }
  else
    {
      for (k1=0; k1 < nsol; k1++)
	{
	  if (k1==0 || res[k1] < resmin)
	    {
	      resmin = res[k1];
	      kmin = k1;	
	    }
	}
      d2 = d2m[kmin];
      l2 = l2m[kmin];
    }
  if(a==0.0 && c==0.0)  // handle a bi-quadratic equation
    {
      /* temporary fix */
     /* if I choose a set of root such as 1,-1,0.0001,-0.0001 d2=0
      * algorithm fails, hence I decided to make flocke algo handle  
      * biquadratics */
      /* solve directly */
      cbq[2]=1.0;
      cbq[1]=b;	
      cbq[0]=d;
      solve_quadratic_cmplx(cbq,sbq);
      roots[0] = csqrt(sbq[0]);
      roots[1] = -csqrt(sbq[0]);
      roots[2] = csqrt(sbq[1]);
      roots[3] = -csqrt(sbq[1]);
      return;
    }

  whichcase = 0; // d2 != 0
  if(d2 < 0.0) 
    {
      gamma=sqrt(-d2);                               
      aq=l1+gamma;                                  
      bq=l3+gamma*l2;                              

      cq=l1-gamma;                                
      dq=l3-gamma*l2;                            
      if(fabs(dq) < fabs(bq))
	dq=d/bq;                                
      else if(fabs(dq) > fabs(bq))
	bq=d/dq;                               

      if (fabs(aq) < fabs(cq))
	{
	  nsol=0;
	  if (dq !=0)
	    {
	      aqv[nsol] = (c - bq*cq)/dq;
	      errv[nsol]=calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	      nsol++;
	    }
	  if (cq != 0) 
	    {
	      aqv[nsol] = (b - dq - bq)/cq;
	      errv[nsol] = calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	      nsol++;
	    }
	  aqv[nsol] = a - cq;
	  errv[nsol] = calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	  nsol++;
	  for (k=0; k < nsol; k++)
	    {
	      if (k==0 || errv[k] < errmin)
		{
		  kmin = k;
		  errmin = errv[k];
		}
	    }
	  aq = aqv[kmin];
	}
      else 
	{
	  nsol = 0;
	  if (bq != 0)
	    { 
	      cqv[nsol] = (c - aq*dq)/bq;
	      errv[nsol] = calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	      nsol++;
	    }
	  if (aq != 0)
	    {
	      cqv[nsol] = (b - bq - dq)/aq;
	      errv[nsol] = calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	      nsol++;
	    }
	  cqv[nsol] = a - aq;
	  errv[nsol] = calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	  nsol++;	  
	  for (k=0; k < nsol; k++)
	    {
	      if (k==0 || errv[k] < errmin)
		{
		  kmin = k;
		  errmin = errv[k];
		}
	    }
	  cq = cqv[kmin];
	}
      realcase[0]=1;
    }
  else if (d2 > 0)   // (d2.gt.0.0) assume a complex quadratic decomposition: 
    {
      // --------------------------------------------------------------------       
      gamma=sqrt(d2); 

      acx=CMPLX(l1,gamma);  
      bcx=CMPLX(l3,gamma*l2);
      ccx = conj(acx);
      dcx = conj(bcx);
      realcase[0] = 0; // complex 
    }
  else 
    realcase[0]=-1; // d2=0
  
  if (fabs(d2) <= macheps*max3(fabs(2.*b/3.), fabs(phi0), l1*l1)) // if d2 is nearly zero it is a special case
    {
      d3 = d - l3*l3;
      if (realcase[0]==1)
	err0 = calc_err_abcd(a, b, c, d, aq, bq, cq, dq);
      else if (realcase[0]==0)
	err0 = calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);

      if (d3 <= 0)
	{
	  realcase[1] = 1;
	  aq1 = l1;
	  bq1 = l3 + sqrt(-d3);
	  cq1 = l1;
	  dq1 = l3 - sqrt(-d3);
	  if(fabs(dq1) < fabs(bq1))
	    dq1=d/bq1;                                        // equation (5.9)
	  else if(fabs(dq1) > fabs(bq1))
	    bq1=d/dq1;                                       // equation (5.10)
	  err1 = calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1);
	}
      else /* caso complesso */
	{
	  printf("qui1 d2=%.16G d3=%.16G phi0=%.16G\n", d2, d-l3*l3, phi0);
	  realcase[1] = 0;
	  //printf("XX c[5]={%.32G,%.32G,%.32G,%.32G,%.32G}\n", d, c, b, a, 1.0);
	  acx1 = l1;
	  bcx1 = l3 + I*sqrt(d3);
	  ccx1 = l1;
	  dcx1 = conj(bcx1);
	  err1 = calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1); 
	}
      if (realcase[0]==-1 || err1 < err0)
	{
          whichcase=1; // d2 = 0
	  if (realcase[1]==1)
	    {
	      aq = aq1;
	      bq = bq1;
	      cq = cq1;
	      dq = dq1;
	    }
	  else
	    {
	      acx = acx1;
	      bcx = bcx1;
	      ccx = ccx1;
	      dcx = dcx1;
	    }
	}
    }
  if (realcase[whichcase]==1)
    {
#ifdef LDLT_USENR
      NRabcd(a,b,c,d,&aq,&bq,&cq,&dq);      
#endif
      myquadratic(a,b,c,d,aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      myquadratic(a,b,c,d,cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];
    }
  else
    {
      /* complex coefficients of p1 and p2 */
      if (whichcase==0) // d2!=0
	{
#ifdef LDLT_USENRCMPLX
	  NRabcdCmplx(a, b, c, d, &acx, &bcx, &ccx, &dcx);
#endif
	  cdiskr=acx*acx/4-bcx;               

	  zx1=-acx/2+csqrt(cdiskr);
	  zx2=-acx/2-csqrt(cdiskr);

	  if(cabs(zx1) > cabs(zx2))
	    zxmax=zx1;
	  else
	    zxmax=zx2;

	  zxmin=bcx/zxmax;        

	  roots[0]=zxmin;
	  roots[1]=conj(zxmin);
	  roots[2]=zxmax;
	  roots[3]=conj(zxmax);
	}
      else // d2 ~ 0
	{
#ifdef USENRCMPLX
	  NRabcdCmplx(a, b, c, d, &acx, &bcx, &ccx, &dcx);
#endif
	  cdiskr=csqrt(acx*acx-4.0*bcx);
	  zx1 = -0.5*(acx+cdiskr);
	  zx2 = -0.5*(acx-cdiskr);
	  if (cabs(zx1) > cabs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = bcx/zxmax;
	  roots[0] = zxmax;
	  roots[1] = zxmin;
	  cdiskr=csqrt(ccx*ccx-4.0*dcx);
	  zx1 = -0.5*(ccx+cdiskr);
	  zx2 = -0.5*(ccx-cdiskr);
	  if (cabs(zx1) > cabs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = dcx/zxmax;
	  roots[2]= zxmax;
	  roots[3]= zxmin;
	}
    }
  //--------------------------------------------------------------------      
}

/* 11/01/18 NOTA: dai test che ho effettuato fast quartic solver (FQS)è circa 3 ordini di grandezza più
 * accurato dell'algoritm hqr() nel trovare gli zeri della quartica. I test li ho fatto calcolando
 * con quale accuratezza la quartica fa zero con le soluzioni reali trovate e FQS è sotto 5E-16 mentre
 * hqr() sotto 5E-13 */
//#define DEBUG_BACKWARD_OPT 
void backward_optimizer(double *alpha, double *beta, double *gamma, double *delta, double a, double b, double c, double d, int *kchosen)
{
#ifdef DEBUG_BACKWARD_OPT
  static long int nmax=0, totcall=0;
#endif
  //double FASTQSEPS=2.2204460492503131E-16;
  double e1[2], e2[2], e3[2], e4[2];
  double U23[2], U33[2], L43[2], U44[2], x1[2], x2[2], x3[2], x4[2], y1[2], y2[2], y3[2], y4[2];
  double eps[2][FQSBUFLEN];
  const int MAXITS=16; 
  //16 è il valore usato nell'articolo Journal of Computational and Applied Mathematics 234 (2010) 3007–3024
  int k, j, its;
#ifdef DEBUG_BACKWARD_OPT
  totcall++;
#endif
  for (k=0; k < FQSBUFLEN; k++)
    {
      eps[0][k] = eps[1][k] = 0;
    }

  for (k=0; k < 2; k++)
    {
      e1[k] = a - alpha[k] - gamma[k];
      e2[k] = b - beta[k] - alpha[k]*gamma[k] - delta[k];
      e3[k] = c - beta[k]*gamma[k] - alpha[k]*delta[k];
      e4[k] = d - beta[k]*delta[k];
    }
  for (its=0; its < MAXITS; its++)
    {
      for (k=0; k < 2; k++)
	{
	  U23[k] = alpha[k] - gamma[k];
	  U33[k] = beta[k] - delta[k] - gamma[k]*U23[k];
	  L43[k] = -delta[k]*U23[k]/U33[k];
	  U44[k] = beta[k] - delta[k] - L43[k]*U23[k];
	  x1[k] = e1[k];
	  x2[k] = e2[k] - gamma[k]*x1[k];
	  x3[k] = e3[k] - delta[k]*x1[k] - gamma[k]*x2[k];
	  x4[k] = e4[k] - delta[k]*x2[k] - L43[k]*x3[k];
	  y4[k] = x4[k]/U44[k];
	  y3[k] = (x3[k]-U23[k]*y4[k])/U33[k];
	  y2[k] = x2[k] - U23[k]*y3[k] - y4[k];
	  y1[k] = x1[k] - y3[k];
	  alpha[k] = alpha[k] + y1[k];
	  beta[k] = beta[k] + y2[k];
	  gamma[k] = gamma[k] + y3[k];
	  delta[k] = delta[k] + y4[k];
	  e1[k] = a - alpha[k] - gamma[k];
	  e2[k] = b - beta[k] - alpha[k]*gamma[k] - delta[k];
	  e3[k] = c - beta[k]*gamma[k] - alpha[k]*delta[k];
	  e4[k] = d - beta[k]*delta[k];
	  // 0 is the latest one
	  /* shift epsilon's */ 
	  for (j=FQSBUFLEN-1; j > 0; j--)
	    eps[k][j] = eps[k][j-1];
	  eps[k][0] = fabs(e1[k])+fabs(e2[k])+fabs(e3[k])+fabs(e3[k]);
	  // convergence
	  if (eps[k][0] == 0.0)
	    {
	      *kchosen=k;
	      //printf("qui1\n");
	      return;
	    }
	  // cyclic condition
	  else if (eps_identical(eps[k]))
	    {
	      //printf("qui2\n");
	      *kchosen=k;
	      return;
	    }
	}
    }
  //printf("max iterations reached eps=%.16G %.16G\n", eps[0][0], eps[1][0]);
#ifdef DEBUG_BACKWARD_OPT
  nmax++;
  if (nmax % 100==0 && totcall > 0)
    {
      printf("nmax=%ld frac=%f\n", nmax, ((double)nmax)/((double)totcall));
      printf("eps=%.16G %.16G\n", eps[0][0], eps[1][0]);
    }
#endif
  if (eps[0][0] < eps[1][0])
    {
      *kchosen = 0;
    }
  else
    {
      *kchosen = 1;
    }
}
int error_handler1(double a, double b, double c, double d, int *numsol, double solqua[4])
{
  double sq[2], cq[3], alpha, beta, eps1, eps2;
  int k, nsq;
  alpha = a*0.5;
  beta = (b - Sqr(alpha))*0.5;
  eps1 = c - 2.0*(alpha)*(beta);
  eps2 = d - Sqr(beta);
  if (eps1==0 && eps2==0)
    {
      /* now we solve the quadratic equation providing the roots of the original quartic */
      cq[2] = 1.0;
      cq[1] = alpha;
      cq[0] = beta;
      solve_quadratic(cq, &nsq, sq);
      *numsol=0;
      for (k=0; k < nsq; k++)
	{
	  solqua[*numsol] = sq[k];
	  (*numsol)++;
	}
      return 1;
    }
  return 0; 
}
int error_handler2(double a, double b , double c, double d, int *numsol, double solqua[4])
{
  double cq[3], x1[2], x2[2], eps1, eps2;
  int nsq, k;
  cq[2] = 1.0;
  cq[1] = a*0.5;
  cq[0] = b/6.0;
  solve_quadratic(cq, &nsq, x1);
  if (nsq==0)
    {
      *numsol=0;
      return 0;
    } 
  x2[0] = -a - 3.0*x1[0];
  x2[1] = -a - 3.0*x1[1];
     
  for (k=0; k < 2; k++)
    {
      eps1 = c + Sqr(x1[k])*(x1[k]+3.0*x2[k]);
      eps2 = d - Sqr(x1[k])*x1[k]*x2[k];
      if (eps1==0 && eps2==0)
	{
	  *numsol=2;
	  solqua[0] = x1[k];
	  solqua[1] = x2[k];	  
	  return 1;
	}
      
    }
  return 0;
   
}
int fast_solver_cmp_func(const void* aa, const void *bb)
{
  complex double ai, bi;
  ai = *((complex double *)aa);
  bi = *((complex double *)bb);
  if (cabs(ai) > cabs(bi))
    return -1;
  else if (cabs(ai)==cabs(bi))
    return 0;
  else
    return 1;

}
void initial_guess_fast_quart_solver(double *alpha, double *beta, double *gamma, double *delta, double a, double b, double c, double d)
{
  int k;
  complex double csol[4];
  double coeff[5];
  double phi1, phi2, c1, c2, L1, L2, L3, y1, y2;
  coeff[4] = 1.0;
  coeff[3] = a;
  coeff[2] = b;
  coeff[1] = c;
  coeff[0] = d;
  //printf("coeff=%f %f %f %f", a, b, c, d);
#if 1
  csolve_quartic_abramovitz_cmplx(coeff, csol);
  //csolve_quartic_ferrari_cmplx(coeff, csol);
#else
  //solve_fourth_deg_cmplx(coeff, csol);
#endif
  qsort(csol, 4, sizeof(complex double), fast_solver_cmp_func);
  alpha[0] = -creal(csol[0]+csol[1]);
  beta[0] = creal(csol[0]*csol[1]);
  alpha[1] = -creal(csol[1]+csol[2]);
  beta[1] = creal(csol[1]*csol[2]);
  for (k=0; k < 2; k++)
    {
      phi1 = 1.0 + Sqr(alpha[k])+Sqr(beta[k]);
      phi2 = alpha[k]*(1.0+beta[k]);
      c1 = a - alpha[k] + alpha[k]*(b-beta[k])+beta[k]*c;
      c2 = b - beta[k] + alpha[k]*c + beta[k]*d;
      L1 = sqrt(phi1);
      L3 = phi2/L1;
      L2 = sqrt(phi1-phi2*phi2/phi1);
      y1 = c1/L1;
      y2 = (c2 - y1*L3)/L2;
      delta[k] = y2/L2;
      gamma[k] = (y1 - delta[k]*L3)/L1;
    }
#if 0
  for (k=0; k < 4; k++)
    {
      printf("csol=%.15G + (%.15G)*I (norm=%.15G)\n", creal(csol[k]), cimag(csol[k]), cabs(csol[k]));
    }
  printf("\n");
  exit(-1);
#endif
}
void fast_quartic_solver(double coeff[5], int *numsol, double solqua[4])
{
  double alpha[2], beta[2], gamma[2], delta[2];
  double a, b, c, d;
  double cq[3], sq[2];
  int k, nsq, setchosen; 
 
  return;
  a = coeff[3]/coeff[4];
  b = coeff[2]/coeff[4];
  c = coeff[1]/coeff[4];
  d = coeff[0]/coeff[4];


  initial_guess_fast_quart_solver(alpha, beta, gamma, delta, a, b, c, d);

  if (error_handler1(a, b, c, d, numsol, solqua))
    return;

  if (error_handler2(a, b, c, d, numsol, solqua))
    return;

  backward_optimizer(alpha, beta, gamma, delta, a, b, c, d, &setchosen);
  /* now we solve the two quadratic equation providing the four roots of the original quartic */
  cq[2] = 1.0;
  cq[1] = alpha[setchosen];
  cq[0] = beta[setchosen];
  solve_quadratic(cq, &nsq, sq);
  *numsol=0;
  for (k=0; k < nsq; k++)
    {
      solqua[*numsol] = sq[k];
      (*numsol)++;
    }
  cq[2] = 1.0;
  cq[1] = gamma[setchosen];
  cq[0] = delta[setchosen];
  solve_quadratic(cq, &nsq, sq);
  for (k=0; k < nsq; k++)
    {
      solqua[*numsol] = sq[k];
      (*numsol)++;
    }
}
#endif
#if 1
void solve_quarticl(long double coeff[5], int *numsol, long double solqua[4])
{
  int ok;
  long double EPS=2.2204460492503131E-32;

  if (coeff[4] == 0)
    {
      solve_cubicl(coeff, numsol, solqua);
      return;
    }
#if 0
  if (fabsl(coeff[4]) < EPS*fabsl(coeff[2]))
    {
      printf("[WARNING] fallback to quadratic from quartic\n");
      solve_numrecl(coeff, numsol, solqua, &ok, 2);
      //solve_cubicl(coeff, numsol, solqua);
      return;
    }
  else if (coeff[4]==0.0)
    {
      printf("[WARNING] fallback to cubic from quartic\n");
      /* N.B. una cubica non puo' essere, i primi due coefficienti devono essere
       * simultaneamente nulli  */
      solve_numrecl(coeff, numsol, solqua, &ok, 3);
      //solve_cubicl(coeff, numsol, solqua);
      return;
    }
#endif
  solve_numrecl(coeff, numsol, solqua, &ok, 4);
}
#else
void solve_quarticl(long double coeff[5], int *numsol, long double solqua[4])
{
  int ok, k;
  double coeffd[5], solquad[4];
  long double EPS=2.2204460492503131E-32;
  for (k=0; k < 5; k++)
    coeffd[k] = (double) coeff[k];

  if (fabsl(coeff[4]) < EPS*fabsl(coeff[2]))
    {
      printf("[WARNING] fallback to quadratic from quartic\n");
      solve_numrec(coeffd, numsol, solquad, &ok, 2);
      for (k=0; k < *numsol; k++)
	solqua[k] = (long double)solquad[k];

      //solve_cubicl(coeff, numsol, solqua);
      return;
    }
#if 0
  else if (coeff[4]==0.0)
    {
      printf("[WARNING] fallback to cubic from quartic\n");
      /* N.B. una cubica non puo' essere, i primi due coefficienti devono essere
       * simultaneamente nulli  */
      solve_numrecl(coeff, numsol, solqua, &ok, 3);
      //solve_cubicl(coeff, numsol, solqua);
      return;
    }
#endif
  solve_numrec(coeffd, numsol, solquad, &ok, 4);
  for (k=0; k < *numsol; k++)
    solqua[k] = (long double)solquad[k];
}
#endif
#ifdef FAST_QUARTIC_SOLVER
complex double quarticfunc(double c[5], double complex sol)
{
  complex double sum=0.0, psol;
  int k;
  
  sum=c[0];
  psol=1.0+I*0.0;
  for (k=1; k < 5; k++)
    {
      psol*=sol;
      sum+=c[k]*psol;
    }
  return sum/cabs(sol);
}

void wrap_LDLT_quartic(double coeff[5], int *numsol, double solqua[4])
{
  //double ERRTHR=1E-10;
  double complex csol[4];
  double sum, sum2;
  long double complex csoll[4];
  long double coeffl[5];
  int k, planB=0;//, printsol;
  //LDLT_quartic(coeff, csol);
  oqs_quartic_solver(coeff,csol); // questo è il quartic solver ultimo sottomesso in ACM trans. math. softw.
  //quarticRoots(coeff, numsol, csol);
  *numsol=0;

#if 0
  sum=0;
  for (k=0; k < 4; k++)
    {
      sum += cabs(quarticfunc(coeff, csol[k]));
    }
#endif

#if 0
  for (k=0; k < 4; k++)
    {
      if (cabsl(csol[k]) > 1E6)
	printf("%.15G+I*(%.15G)\n", creal(csol[k]), cimag(csol[k]));
    }
#endif
  for (k=0; k < 4; k++)
    {
      if (isnan(creal(csol[k]))||isnan(cimag(csol[k]))|| 
	  isinf(creal(csol[k]))||isinf(cimag(csol[k])))
	planB=1;
    }
 // if error is too large use quadruple precision HQR routine as a safe (although much slower) backup */ 
  if (planB)
    {
      //printf("[USING HQRL] BEFORE sum=%.15G\n", sum);
      for (k=0; k < 4; k++)
	{
	  printf("BEFORE csol=%.16G + I*(%.16G)\n", creal(csol[k]), cimag(csol[k]));
	}
      for (k=0; k < 5; k++)
	coeffl[k]=coeff[k];
      solve_numrecl_cmplx(coeffl, csoll, 4);
      for (k=0; k < 4; k++)
	{
	  csol[k] = csoll[k];
	  printf("AFTER csol=%.16G + I*(%.16G)\n", creal(csol[k]), cimag(csol[k]));
	}
      sum2=0;
      for (k=0; k < 4; k++)
	{
	  sum2 += cabs(quarticfunc(coeff, csol[k]));
	}
      printf("[USING HQRL] AFTER sum=%.20G\n", sum2);
    }
#if 0
  printsol=0;
  for (k=0; k < 4; k++)
    if (fabs(cimag(csol[k])) > 5E5)
      {
	printsol=1;
      }
  if (printsol)
    printf("STAMPO SOL:\n");
#endif
  for (k=0; k < 4; k++)
    {
#if 0
      if (printsol)
	{
	  printf("sol=%.15G+(%.15G)*I\n", creal(csol[k]), cimag(csol[k]));
	}
#endif
      if (cimag(csol[k])==0)
	{
	  solqua[*numsol] = creal(csol[k]);
	  (*numsol)++;
	}
    }
}
#endif
void solve_quartic(double coeff[5], int *numsol, double solqua[4])
{
  int ok, k;
  //double EPS=2.2204460492503131E-16;
  double target;
  if (coeff[4]==0)
    {
      solve_cubic(coeff, numsol, solqua);
      return;
    }
#if 0
  if (fabs(coeff[4]) < EPS*fabs(coeff[2]))
    {
      printf("[WARNING] fallback to quadratic from quartic\n");
      //solve_quadratic(coeff, numsol, solqua);
      solve_numrec(coeff, numsol, solqua, &ok, 2);
      return;
    }
  else if (fabs(coeff[4]) < EPS)
    {
      printf("[WARNING] fallback to cubic from quartic\n");
      //solve_quadratic(coeff, numsol, solqua);
      solve_numrec(coeff, numsol, solqua, &ok, 3);
      return;
    }
#endif
#ifdef FAST_QUARTIC_SOLVER
#if 0
    {
      double solquaQR[4], solquaFA[4], solquaQR2[4];
      int numsolQR;
      solve_numrec(coeff, &numsolQR, solquaQR, &ok, 4);
      
#if 0
      printf("BEGIN\n");
      for (k=0; k < *numsol; k++)
	{
	  printf("QR sol[%d]=%.16G\n", k, solqua[k]);
	}
#endif
      //solve_numrec(coeff, numsol, solqua, &ok, 4);
#endif
      wrap_LDLT_quartic(coeff, numsol, solqua);   
      //fast_quartic_solver(coeff, numsol, solqua);
#elif defined(POLY_SOLVE_GSL)
  solve_gslpoly(coeff, numsol, solqua);
#else
#ifdef MC_QUART_USE_ANALYTIC
  csolve_quartic_abramovitz(coeff, numsol, solqua);
#else
  solve_numrec(coeff, numsol, solqua, &ok, 4);
  if (!ok)
    {
      /* if hqr() fails fallback to GSL */
      solve_gslpoly(coeff, numsol, solqua);
    }
#endif
#endif
}
void newt2Dquartic(double c[6], double sol[3], double D2)
{
  double M[2][2], invM[2][2], detM;
  double x[2], dx[2], f[2], fini[2], D2sq;
  int i, j, it;
  D2sq = Sqr(D2);
  x[0] = sol[1]/D2;
  x[1] = sol[2]/D2;

#if 0
  c[0] /= D2sq;
  c[1] /= D2sq;
  c[2] /= D2sq;
  c[4] /= D2;
  c[5] /= D2;
#endif
  fini[0]= f[0] = x[0]*x[0] + x[1]*x[1] - 1.0;
  fini[1]= f[1] = c[0]*x[0]*x[0] + c[1]*x[1]*x[1] + c[2]*x[0]*x[1] + c[3] + c[4]*x[0] + c[5]*x[1];

  for (it = 0; it < 3; it++)
    {
      M[0][0] = 2.0*x[0];
      M[0][1] = 2.0*x[1];
      M[1][0] = c[4] + 2*c[0]*x[0];
      M[1][1] = c[5] + 2*c[1]*x[1];

      detM = -M[0][1]*M[1][0] + M[0][0]*M[1][1];
      if (detM==0)
	return;
      invM[0][0] =  M[1][1]/detM;
      invM[0][1] = -M[0][1]/detM;
      invM[1][0] = -M[1][0]/detM;
      invM[1][1] =  M[0][0]/detM;
      f[0] = x[0]*x[0] + x[1]*x[1] - 1.0;
      f[1] = c[0]*x[0]*x[0] + c[1]*x[1]*x[1] + c[2]*x[0]*x[1] + c[3] + c[4]*x[0] + c[5]*x[1];

      for (i = 0; i < 2; i++)
	{ 
	  dx[i] = 0;
	  for (j = 0; j < 2; j++)
	    {
	      dx[i] += -invM[i][j]*f[j];
	    }
	}
      //printf("dx=%.18G %.18G x=%.18G %.18G\n", dx[0], dx[1], x[0], x[1]);
      for (i = 0; i < 2; i++)
	x[i] += dx[i];
    }
  f[0] = x[0]*x[0] + x[1]*x[1] - 1.0;
  f[1] = c[0]*x[0]*x[0] + c[1]*x[1]*x[1] + c[2]*x[0]*x[1] + c[3] + c[4]*x[0] + c[5]*x[1];

  //printf("fold=%.18G %.18G fnew=%.18G %.18G\n", fini[0], fini[1], f[0], f[1]);
  if (fabs(f[0]) < fabs(fini[0]) && fabs(f[1]) < fabs(fini[1]))
    {
      sol[1] = x[0]*D2;
      sol[2] = x[1]*D2;
      //printf("fold=%.18G %.18G fnew=%.18G %.18G\n", fini[0], fini[1], f[0], f[1]);
   //   printf("fold=%.15G %.15G fnew=%.18G %.18G\n", fini[0], fini[1], f[0], f[1]);
    }
}
double rimdiskone_ibarradiff(double Diami, double Diamj, double Li, double Lj, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, it, k;
  const int MAX_ITERATIONS=1000;
  double Ai[3], AiDj[3], Tnew[3], Told[3], VV[3], AiDjnj, dscpara[3];
  double dscperp[3], ragg, ragg2, TnCi[3]; 
  for (kk1 = 0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
  //printf("distance=%.15G\n", calc_distance(Ai, Dj[j2]));
  for (it = 0; it < MAX_ITERATIONS; it++)
    {
      for(k=0;k<3;k++)
	Told[k] = Tnew[k];
      for (kk1=0; kk1 < 3; kk1++)
	AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
      AiDjnj = scalProd(AiDj, nj);
      for (kk1=0; kk1 < 3; kk1++)
	{
	  VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
	}
      for (kk1=0; kk1 < 3; kk1++)
	dscpara[kk1] = dscperp[kk1] - Dj[kk1];
      ragg = calc_norm(VV);

      for(k=0;k<3;k++)
	{
	  VV[k] = VV[k]/ragg;
	  Tnew[k] = Dj[k] + VV[k]*Diamj*0.5;
	  TnCi[k] = Tnew[k]-Ci[k];
	}

      ragg = scalProd(TnCi,ni);
      for (k=0;k<3;k++)
	Ai[k] = Ci[k] + ragg*ni[k];	
#if 0
      for (k=0;k<3;k++)
	Bi[k] = Tnew[k]-Ci[k]-ragg*ni[k];

      ragg2 = calc_norm(Bi);

      if ((fabs(ragg) < L*0.5) && ((ragg2) < D*0.5))
	return -1;
#endif
      if ( it > 0 && check_convergence(Told,Tnew) ) 
	break;
    } 
  for(k=0;k<3;k++)
    TnCi[k] = Tnew[k]-Ci[k];

  ragg = scalProd(TnCi,ni);

  for (k=0;k<3;k++)
    Ai[k] = Tnew[k]-Ci[k]-ragg*ni[k];

  ragg2 = calc_norm(Ai);

  if ((fabs(ragg) < Li*0.5) && ((ragg2) < Diami*0.5))
    return -1;

  return 1;
}
double rimdiskone_ibarra(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, it, k;
  const int MAX_ITERATIONS=1000;
  double Bi[3], Ai[3], AiDj[3], Tnew[3], Told[3], VV[3], AiDjnj, dsc[3];//dscpara[3], dsc[3];
  double ragg, ragg2, TnCi[3]; 
  for (kk1 = 0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
  //printf("distance=%.15G\n", calc_distance(Ai, Dj[j2]));
  for (it = 0; it < MAX_ITERATIONS; it++)
    {
      for(k=0;k<3;k++)
	Told[k] = Tnew[k];
      for (kk1=0; kk1 < 3; kk1++)
	AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
      AiDjnj = scalProd(AiDj, nj);
      for (kk1=0; kk1 < 3; kk1++)
	{
	  VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
	}
      //    for (kk1=0; kk1 < 3; kk1++)
      //	dscpara[kk1] = dscperp[kk1] - Dj[kk1];
      ragg = calc_norm(VV);

      for(k=0;k<3;k++)
	{
	  VV[k] = VV[k]/ragg;
	  Tnew[k] = Dj[k] + VV[k]*D*0.5;
	  TnCi[k] = Tnew[k]-Ci[k];
	}

      ragg = scalProd(TnCi,ni);
      for (k=0;k<3;k++)
	Ai[k] = Ci[k] + ragg*ni[k];	
#if 0
      for (k=0;k<3;k++)
	Bi[k] = Tnew[k]-Ci[k]-ragg*ni[k];

      ragg2 = calc_norm(Bi);

      if ((fabs(ragg) < L*0.5) && ((ragg2) < D*0.5))
	return -1;
#endif
      if ( it > 0 && check_convergence(Told,Tnew) ) 
	break;
    } 
  for(k=0;k<3;k++)
    TnCi[k] = Tnew[k]-Ci[k];

  ragg = scalProd(TnCi,ni);

  for (k=0;k<3;k++)
    Ai[k] = Tnew[k]-Ci[k]-ragg*ni[k];

  ragg2 = calc_norm(Ai);

  if ((fabs(ragg) < L*0.5) && ((ragg2) < D*0.5))
    return -1;

  return 1;
}
void discard_spuriousl(long double *solqua, int *numsol)
{
  /* each solution x must be such that |x| <= 1 */
  int k, nsol;
  const long double EPS=1E-9;
  long double solL[4];
  nsol=0;
  for (k=0; k < *numsol; k++)
    solL[k] = solqua[k];
  for (k=0; k < *numsol; k++)
    {
      if (fabsl(solL[k]) < 1.0+EPS)
	{
	  solqua[nsol] = solL[k];
	  nsol++;
	}
    }
  //printf("numsol=%d nsol=%d\n", *numsol, nsol);
  *numsol = nsol;
}
int test_for_fallbackldiff(long double *P, long double *Cip, long double *nip, long double Di2, long double Dj2, long double *diff)
{
  const long double DIST_THR=5E-15;
  long double diff1, diff2;
  diff1=fabsl(perpcompl(P, Cip, nip)-Di2); // qui D2 è il diametro del rim
  diff2=fabsl(sqrtl(Sqr(P[1])+Sqr(P[2]))-Dj2);// qui D2 è il diametro del disco

  *diff=diff1+diff2;
  if (diff1 > DIST_THR*Di2 || diff2 > DIST_THR*Dj2)
    return 1;
  else
    return 0;
#if 0
  if ((*diff=fabsl(perpcompl(P, Cip, nip)-D2)) > DIST_THR*D2)
    return 1;
  else 
    return 0;
#endif
}
int test_for_fallbackl(long double *P, long double *Cip, long double *nip, long double D2, long double *diff)
{
  const long double DIST_THR=5E-15;
  long double diff1, diff2;
  diff1=fabsl(perpcompl(P, Cip, nip)-D2); // qui D2 è il diametro del rim
  diff2=fabsl(sqrtl(Sqr(P[1])+Sqr(P[2]))-D2);// qui D2 è il diametro del disco

  *diff=diff1+diff2;
  if (diff1 > DIST_THR*D2 || diff2 > DIST_THR*D2)
    return 1;
  else
    return 0;
#if 0
  if ((*diff=fabsl(perpcompl(P, Cip, nip)-D2)) > DIST_THR*D2)
    return 1;
  else 
    return 0;
#endif
}

#if 0
int test_solution_xyl(long double *sol)
{
  const double EPS=2.2204460492503131E-16;
  if (fabsl(1.0 - (Sqr(sol[0]) + Sqr(sol[1]))) > EPS)
     return 1;
   else 
     return 0;
}

int test_solution_xy(double *sol, double *diff)
{
  const double EPS=1E-7;
  //printf("diff=%.16G\n", fabs(1.0 - (Sqr(sol[0]) + Sqr(sol[1]))));
  if ((*diff=fabs(1.0 - (Sqr(sol[0]) + Sqr(sol[1])))) > EPS)
     return 1;
   else 
     return 0;
}
#endif
void discard_spurious(double *solqua, int *numsol)
{
  /* each solution x must be such that |x| <= 1 */
  int k, nsol;
  const double EPS=5E-8; // it should be at least not greater than DIST_THR below
  double solL[4];
  nsol=0;
  for (k=0; k < *numsol; k++)
    solL[k] = solqua[k];
  for (k=0; k < *numsol; k++)
    {
      if (fabs(solL[k]) < 1.0+EPS)
	{
	  solqua[nsol] = solL[k];
	  nsol++;
	}
    }
#if 0
  if (numsol==4)
    printf("numsol=%d nsol=%d\n", *numsol, nsol);
#endif
  *numsol = nsol;
}
/* nel caso di cilindri diversi bisogna passare i diametri del disco e del rim e non solo D2*/
int test_for_fallbackdiff(double *P, double *Cip, double *nip, double Di2, double Dj2, double *diff)
{
  const double DIST_THR=5.0E-12;
  double diff1, diff2;
  diff1=fabs(perpcomp(P, Cip, nip)-Di2); // qui D2 è il diametro del rim
  diff2=fabs(sqrt(Sqr(P[1])+Sqr(P[2]))-Dj2);// qui D2 è il diametro del disco

  *diff=diff1+diff2;
  if (diff1 > DIST_THR*Di2 || diff2 > DIST_THR*Dj2)
    return 1;
  else
    return 0;
#if 0
  const double DIST_THR=5E-8;
  if ((*diff=fabs(perpcomp(P, Cip, nip)-D2)) > DIST_THR*D2)
    return 1;
  else 
    return 0;
#endif
}


int test_for_fallback(double *P, double *Cip, double *nip, double D2, double *diff)
{
#ifdef MC_QUART_USE_ANALYTIC
  const double DIST_THR=1E-4;
#else
#ifdef MC_QUART_HYBRID
#ifdef FAST_QUARTIC_SOLVER
  const double DIST_THR=5E-13;// old value 5.0E-13 -- NOTA 22/02/18 con il nuovo quartic solver ci si può spingere a 5E-14 volendo!
  // NOTA 19/04/19: 5E-12 sarebbe comunque un valore più sicuro specialmente quando si fanno calcoli di volume escluso 
#else
  const double DIST_THR=5E-12;
#endif
#else
  const double DIST_THR=5E-12;
#endif
#endif
  double diff1, diff2;
  diff1=fabs(perpcomp(P, Cip, nip)-D2); // qui D2 è il diametro del rim
  diff2=fabs(sqrt(Sqr(P[1])+Sqr(P[2]))-D2);// qui D2 è il diametro del disco

  *diff=diff1+diff2;
  if (diff1 > DIST_THR*D2 || diff2 > DIST_THR*D2)
    return 1;
  else
    return 0;
#if 0
  const double DIST_THR=5E-8;
  if ((*diff=fabs(perpcomp(P, Cip, nip)-D2)) > DIST_THR*D2)
    return 1;
  else 
    return 0;
#endif
}
extern void versor_to_Rl(long double ox, long double oy, long double oz, long double R[3][3]);

/* NOTA: anche in rimdisonediff test_for_fallback va riscritta poichè il test ha bisogno sia di Diami che di Diamj */
double rimdiskoneldiff(double Diamis, double Diamjs, double Lis, double Ljs, double Cis[3], double nis[3], double Djs[3], double njs[3], double DjCinis)
{
  long double Diami, Diamj, Li, Lj, Ci[3], ni[3], nj[3], DjCini, Dj[3], temp;
  int kk1, kk2, numsol[2], nsc, fallback, solset;
  long double diff[2][4], maxdiff[2], sumdiff[2], tmp;
//  double coeffs[5], solquas[4];
  long double sp, coeff[5],solarr[2][4][3], solec[4][2], solqua[4];
  long double dsc[3], c0, c1, c2, c3, c02, c12, c22,  coeffEr[6], norm, c32, c42, c52, c4, c5;  
  long double Cip[3], nip[3];
  long double nip02,nip12,nip22,nip03,nip13,nip23,nip04,nip14,nip24,Cip02,Cip12,Cip22;
  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  long double Di2sq, Dj2sq, Di2, Dj2, Cip0, Cip1, Cip2, nip0, nip1 , nip2, Rl[3][3]; 
/* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */


  DjCini = (long double)DjCinis;
  Diami      = (long double)Diamis;
  Diamj      = (long double)Diamjs;

  Li      = (long double)Lis;
  Lj      = (long double)Ljs;
  for (kk1=0; kk1 < 3; kk1++)
    {
      ni[kk1] = ((long double)nis[kk1]);
      nj[kk1] = ((long double)njs[kk1]);
      Ci[kk1] = ((long double)Cis[kk1]);
      Dj[kk1] = ((long double)Djs[kk1]);
    } 
  Di2 = Diami*0.5;
  Dj2 = Diamj*0.5;
  Di2sq = Sqr(Di2); 
  Dj2sq = Sqr(Dj2);
  /* mi metto nel riferimento del disco (p) */
  versor_to_Rl(nj[0], nj[1], nj[2], Rl);
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[kk1] = 0;
      //Aip[kk1] = 0;
      Cip[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
#if 0
  printf("Rl=%.15LG %.15LG %.15LG\n", Rl[0][0], Rl[0][1], Rl[0][2]);
  printf("Rl=%.15LG %.15LG %.15LG\n", Rl[1][0], Rl[1][1], Rl[1][2]);
  printf("Rl=%.15LG %.15LG %.15LG\n", Rl[2][0], Rl[2][1], Rl[2][2]);
#endif
  norm = calc_norml(nip);
  nip0 = nip[0]/norm;
  nip1 = nip[1]/norm;
  nip2 = nip[2]/norm;
  //printf("norm=%.16G\n", sqrt(nip0*nip0+nip1*nip1+nip2*nip2));
  Cip0 = Cip[0];
  Cip1 = Cip[1];
  Cip2 = Cip[2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);   
#if 0
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
 /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23);
  /* N.B. sistemare qui Dj2 e Di2 */ 
  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(Di2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
 
  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
 
  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;
 
#endif
  /* check ellipse */
#if 0
  {
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
    double cq[3], x, lam, solq[2], p[3];
    int numsol;
    lam = -Cip0/nip0;
    x = Cip1 + lam*nip1;
    cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
    cq[1] = coeffEr[5]+coeffEr[2]*x;
    cq[2] = coeffEr[1];
    solve_quadratic(cq, &numsol, solq);
    p[0] = 0.0;
    p[1] = x;
    p[2] = solq[0];
    if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
      {
	printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      	printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
	printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
      }
  }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= Dj2sq;
  coeffEr[1] *= Dj2sq; 
  coeffEr[2] *= Dj2sq;
  coeffEr[4] *= Dj2;
  coeffEr[5] *= Dj2;
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#if 0
  for (kk1=0; kk1 < 5; kk1++)
    coeffs[kk1] = ((double)coeff[kk1]);
#endif
  solve_quarticl(coeff, &(numsol[0]), solqua);
  discard_spuriousl(solqua, &(numsol[0]));

#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    solqua[kk1] = (long double)solquas[kk1];
#endif
  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
#if 0
      c0l = c0;
      c1l = c1;
      c2l = c2;
      c3l = c3;
      c4l = c4;
      c5l = c5;
      solqual=solqua[kk1];
      templ = c4l + c2l*solqual;
      solec[kk1][0]=((double)((-c0l - c3l - c5l*solqual + (c0l - c1l)*Sqr(solqual))/templ));
      solec[kk1][1] = solqua[kk1];
#else
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
      if (temp==0)// || test_solution_xyl(solec[kk1]))
	fallback = 1;
    }
#if 0
      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
	     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
	     coeff[1]*solqua[kk1]+coeff[0]);
      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
#if 0
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      solve_quartic(coeff, &numsol, solqua);
      for (kk1=0; kk1 < numsol; kk1++)
	{
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	}
    }
#endif
  sumdiff[0]=maxdiff[0]=0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = Dj2*solec[kk1][0];
      solarr[0][kk1][2] = Dj2*solec[kk1][1];
      if (test_for_fallbackldiff(solarr[0][kk1], Cip, nip, Di2, Dj2, &(diff[0][kk1])))
	{
	  fallback=1;
	}	
      sumdiff[0] += diff[0][kk1];  
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];
      //ellips2disk(solec[kk1], solarr[kk1], 0, 0, D2, D2);
    }
  if (tinyimagGBL)
    {
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
#if 1
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#if 0
      for (kk1=0; kk1 < 5; kk1++)
	coeffs[kk1] = (double)coeff[kk1];
#endif
      solve_quarticl(coeff, &(numsol[1]), solqua);
      discard_spuriousl(solqua, &(numsol[1]));

#if 0
      for (kk1=0; kk1 < numsol; kk1++)
	solqua[kk1] = (long double)solquas[kk1];
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0
    	  c0l = c0;
	  c1l = c1;
	  c2l = c2;
	  c3l = c3;
	  c4l = c4;
	  c5l = c5;
	  solqual=solqua[kk1];
	  templ = c5l + c2l*solqual;
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = ((double)((-c1l - c3l - c4l*solqual + (c1l - c0l)*Sqr(solqual))/templ)); 
#else
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif
	}
      sumdiff[1]=maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = Dj2*solec[kk1][0];
	  solarr[1][kk1][2] = Dj2*solec[kk1][1];
	  test_for_fallbackldiff(solarr[1][kk1], Cip, nip, Di2, Dj2, &(diff[1][kk1]));
	  //ellips2disk(solec[kk1], solarr[kk1], 0, 0, D2, D2);
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1]; 
	}
      if (fallback==2)
	solset=1;
      else
	{
    	  if (maxdiff[1] < maxdiff[0])
	  //if (sumdiff[1] < sumdiff[0])
    	    solset = 1;
	  else 
	    solset = 0;
	}
    }
#endif

#if 0
  construct_inner_points(solarr, Ci, ni, Dj, nj, D);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
#if 0
      printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
      printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 0
      if (calc_normsq(solarr[kk1])-Sqr(D2) > NEWT_THR)
	{
	  newt2Dquartic(coeffEr, solarr[kk1], D2);
	}
#endif
#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallbackldiff(solarr[solset][kk1], Cip, nip, Di2, Dj2, &tmp)) 
	{
	  printf("distanza punto-centro disk: %.15LG\n", calc_norml(solarr[solset][kk1]));
#if 1
	  printf("distanza punto-centro disksq: %.15LG D2^2=%.15LG\n", calc_norml(solarr[solset][kk1]), Sqr(Dj2));
	  printf("BOH2BOH2 perpcom=%.15LG\n", perpcompl(solarr[solset][kk1], Cip, nip));
	  printf("Cip1=%15LG Cip2=%.15LG\n", Cip[1], Cip[2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  //print_vec("ni=",ni);
	  //print_vec("nj=",nj);
	  printf("c02=%.15LG c0=%.15LG c1=%.15LG c12=%.15LG c22=%.15LG\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15LG c5=%.15LG\n", c4, c5);
	  printf("solec[%d]=%.15LG\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16LG %.16LG %.16LG %.16LG %.16LG %.16LG\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
#endif
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15LG\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15LG\n", scalProdl(ni,nj));
	  printf("(%.15LG)*x^4+(%.15LG)*x^3+(%.15LG)*x^2+(%.15LG)*x+(%.15LG)\n", coeff[4], coeff[3], coeff[2], 
		 coeff[1], coeff[0]);
	  printf("{%.15LG,%.15LG,%.15LG,%.15LG,%.15LG}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15LG\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15LG\n", temp);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProdl(dsc, nip);
      if (fabsl(sp) < Li*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
void versor_to_R_altl(long double *Ci, long double *ni, long double *Dj, long double *nj, long double R[3][3], long double D);

double rimdiskonel(double Ds, double Ls, double Cis[3], double nis[3], double Djs[3], double njs[3], double DjCinis)
{
  long double D, L, Ci[3], ni[3], nj[3], DjCini, Dj[3], temp;
  int kk1, kk2, numsol[2], fallback, solset;
  long double diff[2][4], maxdiff[2], sumdiff[2], tmp;
//  double coeffs[5], solquas[4];
  long double sp, coeff[5],solarr[2][4][3], solec[4][2], solqua[4];
  long double dsc[3], c0, c1, c2, c3, c02, c12, c22,  coeffEr[6], norm, c32, c42, c52, c4, c5;  
  long double Cip[3], nip[3];
  long double nip02,nip12,nip22,nip03,nip13,nip23,nip04,nip14,nip24,Cip02,Cip12,Cip22;
  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  long double D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2, Rl[3][3]; 
/* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */


  DjCini = (long double)DjCinis;
  D      = (long double)Ds;
  L      = (long double)Ls;
  for (kk1=0; kk1 < 3; kk1++)
    {
      ni[kk1] = ((long double)nis[kk1]);
      nj[kk1] = ((long double)njs[kk1]);
      Ci[kk1] = ((long double)Cis[kk1]);
      Dj[kk1] = ((long double)Djs[kk1]);
    } 
  D2 = D*0.5; 
  D2sq = Sqr(D2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_Rl(nj[0], nj[1], nj[2], Rl);
#else
  versor_to_R_altl(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[kk1] = 0;
      //Aip[kk1] = 0;
      Cip[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
#if 0
  printf("Rl=%.15LG %.15LG %.15LG\n", Rl[0][0], Rl[0][1], Rl[0][2]);
  printf("Rl=%.15LG %.15LG %.15LG\n", Rl[1][0], Rl[1][1], Rl[1][2]);
  printf("Rl=%.15LG %.15LG %.15LG\n", Rl[2][0], Rl[2][1], Rl[2][2]);
#endif
  norm = calc_norml(nip);
  nip0 = nip[0]/norm;
  nip1 = nip[1]/norm;
  nip2 = nip[2]/norm;
  //printf("norm=%.16G\n", sqrt(nip0*nip0+nip1*nip1+nip2*nip2));
  Cip0 = Cip[0];
  Cip1 = Cip[1];
  Cip2 = Cip[2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);   
#if 0
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
 /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23);
 
  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
 
  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
 
  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;
 
#endif
  /* check ellipse */
#if 0
  {
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
    double cq[3], x, lam, solq[2], p[3];
    int numsol;
    lam = -Cip0/nip0;
    x = Cip1 + lam*nip1;
    cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
    cq[1] = coeffEr[5]+coeffEr[2]*x;
    cq[2] = coeffEr[1];
    solve_quadratic(cq, &numsol, solq);
    p[0] = 0.0;
    p[1] = x;
    p[2] = solq[0];
    if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
      {
	printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      	printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
	printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
      }
  }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2sq;
  coeffEr[1] *= D2sq; 
  coeffEr[2] *= D2sq;
  coeffEr[4] *= D2;
  coeffEr[5] *= D2;
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#if 0
  for (kk1=0; kk1 < 5; kk1++)
    coeffs[kk1] = ((double)coeff[kk1]);
#endif
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyll(Ci, ni, Dj, nj, L, D, D);
    }
  else
    solve_quarticl(coeff, &(numsol[0]), solqua);
  discard_spuriousl(solqua, &(numsol[0]));

#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    solqua[kk1] = (long double)solquas[kk1];
#endif
  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
#if 0
      c0l = c0;
      c1l = c1;
      c2l = c2;
      c3l = c3;
      c4l = c4;
      c5l = c5;
      solqual=solqua[kk1];
      templ = c4l + c2l*solqual;
      solec[kk1][0]=((double)((-c0l - c3l - c5l*solqual + (c0l - c1l)*Sqr(solqual))/templ));
      solec[kk1][1] = solqua[kk1];
#else
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
      if (temp==0)// || test_solution_xyl(solec[kk1]))
	fallback = 1;
    }
#if 0
      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
	     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
	     coeff[1]*solqua[kk1]+coeff[0]);
      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
#if 0
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      solve_quartic(coeff, &numsol, solqua);
      for (kk1=0; kk1 < numsol; kk1++)
	{
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	}
    }
#endif
  sumdiff[0]=maxdiff[0]=0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2*solec[kk1][0];
      solarr[0][kk1][2] = D2*solec[kk1][1];
      if (test_for_fallbackl(solarr[0][kk1], Cip, nip, D2, &(diff[0][kk1])))
	{
	  fallback=1;
	}	
      sumdiff[0] += diff[0][kk1];  
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];
      //ellips2disk(solec[kk1], solarr[kk1], 0, 0, D2, D2);
    }
  if (tinyimagGBL)
    {
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
#if 1
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#if 0
      for (kk1=0; kk1 < 5; kk1++)
	coeffs[kk1] = (double)coeff[kk1];
#endif
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyll(Ci, ni, Dj, nj, L, D, D);
	}
      else
	solve_quarticl(coeff, &(numsol[1]), solqua);
      discard_spuriousl(solqua, &(numsol[1]));

#if 0
      for (kk1=0; kk1 < numsol; kk1++)
	solqua[kk1] = (long double)solquas[kk1];
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0
    	  c0l = c0;
	  c1l = c1;
	  c2l = c2;
	  c3l = c3;
	  c4l = c4;
	  c5l = c5;
	  solqual=solqua[kk1];
	  templ = c5l + c2l*solqual;
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = ((double)((-c1l - c3l - c4l*solqual + (c1l - c0l)*Sqr(solqual))/templ)); 
#else
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif
	}
      sumdiff[1]=maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2*solec[kk1][0];
	  solarr[1][kk1][2] = D2*solec[kk1][1];
	  test_for_fallbackl(solarr[1][kk1], Cip, nip, D2, &(diff[1][kk1]));
	  //ellips2disk(solec[kk1], solarr[kk1], 0, 0, D2, D2);
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1]; 
	}
      if (fallback==2)
	solset=1;
      else
	{
    	  if (maxdiff[1] < maxdiff[0])
	  //if (sumdiff[1] < sumdiff[0])
    	    solset = 1;
	  else 
	    solset = 0;
	}
    }
#endif

#if 0
  construct_inner_points(solarr, Ci, ni, Dj, nj, D);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
#if 0
      printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
      printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 0
      if (calc_normsq(solarr[kk1])-Sqr(D2) > NEWT_THR)
	{
	  newt2Dquartic(coeffEr, solarr[kk1], D2);
	}
#endif
#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallbackl(solarr[solset][kk1], Cip, nip, D2, &tmp)) 
	{
	  printf("distanza punto-centro disk: %.15LG\n", calc_norml(solarr[solset][kk1]));
#if 1
	  printf("distanza punto-centro disksq: %.15LG D2^2=%.15LG\n", calc_norml(solarr[solset][kk1]), Sqr(D2));
	  printf("BOH2BOH2 perpcom=%.15LG\n", perpcompl(solarr[solset][kk1], Cip, nip));
	  printf("Cip1=%15LG Cip2=%.15LG\n", Cip[1], Cip[2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  //print_vec("ni=",ni);
	  //print_vec("nj=",nj);
	  printf("c02=%.15LG c0=%.15LG c1=%.15LG c12=%.15LG c22=%.15LG\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15LG c5=%.15LG\n", c4, c5);
	  printf("solec[%d]=%.15LG\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16LG %.16LG %.16LG %.16LG %.16LG %.16LG\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
#endif
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15LG\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15LG\n", scalProdl(ni,nj));
	  printf("(%.15LG)*x^4+(%.15LG)*x^3+(%.15LG)*x^2+(%.15LG)*x+(%.15LG)\n", coeff[4], coeff[3], coeff[2], 
		 coeff[1], coeff[0]);
	  printf("{%.15LG,%.15LG,%.15LG,%.15LG,%.15LG}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15LG\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15LG\n", temp);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProdl(dsc, nip);
      if (fabsl(sp) < L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
double rimdiskonediff(double Diami, double Diamj, double Li, double Lj, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], fallback, solset;
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4];
  double dsc[3], c0, c1, c2, c3, c02, c12, c22,  coeffEr[6], c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2];
  double Cip[3], nip[3], norm, Rl[3][3];
  double nip02,nip12,nip22,nip03,nip13,nip23,nip04,nip14,nip24,Cip02,Cip12,Cip22, temp;
  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double Di2sq, Dj2sq, Di2, Dj2, Cip0, Cip1, Cip2, nip0, nip1 , nip2; 
/* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  Dj2 = Diamj*0.5; 
  Di2 = Diami*0.5;
  Dj2sq = Sqr(Dj2);
  Di2sq = Sqr(Di2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  versor_to_R_alt(Ci, ni, Dj, nj, Rl, Diamj); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[kk1] = 0;
      //Aip[kk1] = 0;
      Cip[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip);
  nip0 = nip[0]/norm;
  nip1 = nip[1]/norm;
  nip2 = nip[2]/norm;
  Cip0 = Cip[0];
  Cip1 = Cip[1];
  Cip2 = Cip[2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);   
#if 1
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  /* sistemare quest con i giusti Dj2sq e Di2sq */
  coeffEr[3] = Cip02 + Cip12 + Cip22 - Di2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
 /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23);
 
  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
 
  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
 
  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;
 
#endif
  /* check ellipse */
#if 0
  {
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
    double cq[3], x, lam, solq[2], p[3];
    int numsol;
    lam = -Cip0/nip0;
    x = Cip1 + lam*nip1;
    cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
    cq[1] = coeffEr[5]+coeffEr[2]*x;
    cq[2] = coeffEr[1];
    solve_quadratic(cq, &numsol, solq);
    p[0] = 0.0;
    p[1] = x;
    p[2] = solq[0];
    if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
      {
	printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      	printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
	printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
      }
  }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= Dj2sq;
  coeffEr[1] *= Dj2sq; 
  coeffEr[2] *= Dj2sq;
  coeffEr[4] *= Dj2;
  coeffEr[5] *= Dj2;
  //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
  if (coeff[4]==0)
    {
      /* N.B. 08/01/18 forse così è troppo restrittiva e dovrò cambiare in una condizione del tipo 
       * fabs(coeff[4]) < EPSILON, devo fare delle prove per stabilirlo... */
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, Li, Diami, Diamj);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
#if 0
  if (numsol==1)
    {
      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
      printf("sol=%.15G\n", solqua[0]);
      printf("BOH\n");
    }
#endif
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
#if 0
      c0l = c0;
      c1l = c1;
      c2l = c2;
      c3l = c3;
      c4l = c4;
      c5l = c5;
      solqual=solqua[kk1];
      templ = c4l + c2l*solqual;
      solec[kk1][0]=((double)((-c0l - c3l - c5l*solqual + (c0l - c1l)*Sqr(solqual))/templ));
      solec[kk1][1] = solqua[kk1];
#else
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
      //printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
#if 0
      if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	{
  	  printf("solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	  printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
	  printf("numeratore=%.16G\n", -c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]));
	}
#endif
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
#if 0
      if (test_solution_xy(solec[kk1], &(diffxy[0][kk1])))
	  fallback=1;
#endif
      if (temp==0.0) 
	{
	  fallback=1;
	}
    }
#if 0
      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
	     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
	     coeff[1]*solqua[kk1]+coeff[0]);
      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
#if 0
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      solve_quartic(coeff, &numsol, solqua);
      for (kk1=0; kk1 < numsol; kk1++)
	{
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	}
    }
#endif
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = Dj2*solec[kk1][0];
      solarr[0][kk1][2] = Dj2*solec[kk1][1];
#if 1
      if (test_for_fallbackdiff(solarr[0][kk1], Cip, nip, Di2, Dj2, &(diff[0][kk1])))
	{
	  fallback=1;
#if 0
	  if (numsol==4)
	    {
	      printf("%d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		     coeff[1]*solqua[kk1]+coeff[0]);
	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[0][kk1]);
	      printf(">>>> =====================\n");
	    }
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
#endif
    }
  if (tinyimagGBL)
    {
      //printf("BOHHHH\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
#if 1
  if (fallback)
    {
      //printf("falling back\n");
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, Li, Diami, Diamj);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));

      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0
    	  c0l = c0;
	  c1l = c1;
	  c2l = c2;
	  c3l = c3;
	  c4l = c4;
	  c5l = c5;
	  solqual=solqua[kk1];
	  templ = c5l + c2l*solqual;
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = ((double)((-c1l - c3l - c4l*solqual + (c1l - c0l)*Sqr(solqual))/templ)); 
#else
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = Dj2*solec[kk1][0];
	  solarr[1][kk1][2] = Dj2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
#if 1
	  test_for_fallbackdiff(solarr[1][kk1], Cip, nip, Di2, Dj2, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  

#endif
#if 0
  	  if (numsol==4)
  	    {
  	      printf("FALLBACK %d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
  	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
  	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
  	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
  	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
  	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
  		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
  		     coeff[1]*solqua[kk1]+coeff[0]);
  	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[1][kk1]);
	      printf(">>>> =====================\n");
	  }
#endif
	
	}
#if 1
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	  //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
#if 0
	  if (fallback==3 && solset != 1)
	    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
	}
#endif
    }
#endif
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
#if 0
  construct_inner_points(solarr, Ci, ni, Dj, nj, D);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
#if 0
      printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
      printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 0
      if (calc_normsq(solarr[kk1])-Sqr(D2) > NEWT_THR)
	{
	  newt2Dquartic(coeffEr, solarr[kk1], D2);
	}
#endif
#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallbackdiff(solarr[solset][kk1], Cip, nip, Di2, Dj2, &tmp)) 
	{
	  printf("# %d ===================== <<<< \n", kk1);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
#if 1
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(Dj2));
	  printf("BOH2BOH2 perpcom=%.15G\n", perpcomp(solarr[solset][kk1], Cip, nip));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[1], Cip[2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
#endif
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProd(dsc, nip);
      if (fabs(sp) < Li*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
void versor_to_R_altl(long double *Ci, long double *ni, long double *Dj, long double *nj, long double R[3][3], long double D)
{
  int k, kk1, kk;
  long double u[3]; 
  long double normDjCi, DjCi[3], DjCini, Ai[3], AiDjnj, AiDj[3], VV[3], ragg;
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

#if 1
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  normDjCi = calc_norml(DjCi);
  DjCini = scalProdl(DjCi,ni);

  for (kk1 = 0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
  for (kk1=0; kk1 < 3; kk1++)
    AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
  AiDjnj = scalProdl(AiDj, nj);
  for (kk1=0; kk1 < 3; kk1++)
    {
      VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
    }
  //for (kk1=0; kk1 < 3; kk1++)
    //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
  ragg = calc_norml(VV);

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k]/ragg;
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
#if 0
  ragg = scalProd(TnCi,ni);
  for (k=0;k<3;k++)
    Ai[k] = Ci[k] + ragg*ni[k];
#endif
#else

 for (k=0; k < 3; k++)
    dsc[k] = Ci[k] - Dj[k]; 
  sp = scalProd(dsc, nj);
  for (k=0; k < 3; k++)
    R[1][k] = dsc[k] - sp*nj[k];
#endif  
  //printf("scalProd=%.15G\n", scalProd(R[1],R[0]));
  vectProdVecl(R[0], R[1], u);
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
extern double ranf(void);
void rotate_axes_on_plane(double RR[3][3])
{
  double Rin[3][3];
  int k1, k2, k3;
  double ox, oy, oz, theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];

  ox = RR[0][0];
  oy = RR[0][1];
  oz = RR[0][2];
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Rin[k1][k2]=RR[k1][k2];
      }
  //NOTA 15/01/2018
  //se uso ranf altero la sequenza casuale e perdo il confronto la simulazione da 50x10^6 già fatta con il metodo di Alberto 
  //theta = (ranf()>0.5?1.:-1.)*M_PI/4.0;
  theta = M_PI/4.0;
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
	Ro[k1][k2] = Rin[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += Rin[k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     RR[k1][k2] = Ro[k1][k2];
}
void versor_to_R_alt_fb(double *Ci, double *ni, double *Dj, double *nj, double R[3][3], double D, double *Tj, int MAXITS)
{
  int k, kk1, kk, its;
  double u[3]; 
  double normDjCi, DjCi[3], DjCini, Ai[3], AiDjnj, AiDj[3], VV[3], ragg;
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

  /* N.B. Questa viene utilizzata in fallback e faccio un ulteriore step dell'algorito di Ibarra semplificato per
   * determinare l'asse y del riferimenti del disco (l'asse x è l'asse perpendicolare
   * al disco e l'asse z si ottiene con il prodotto vettore dell'asse x e y) */
#if 1
  for (its = 0; its < MAXITS; its++)
    {
      for (kk=0; kk < 3; kk++)
	DjCi[kk] = Dj[kk] - Ci[kk];
      normDjCi = calc_norm(DjCi);
      DjCini = scalProd(DjCi,ni);

      for (kk1 = 0; kk1 < 3; kk1++)
	Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
      for (kk1=0; kk1 < 3; kk1++)
	AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
      AiDjnj = scalProd(AiDj, nj);
      for (kk1=0; kk1 < 3; kk1++)
	{
	  VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
	}
      //for (kk1=0; kk1 < 3; kk1++)
      //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
      ragg = calc_norm(VV);
      for(k=0;k<3;k++)
	{
	  VV[k] /= ragg;
	  Tj[k] = D*0.5*VV[k];

	}
    }

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k];
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
#if 0
      ragg = scalProd(TnCi,ni);
      for (k=0;k<3;k++)
	Ai[k] = Ci[k] + ragg*ni[k];
#endif
#else

      for (k=0; k < 3; k++)
	dsc[k] = Ci[k] - Dj[k]; 
      sp = scalProd(dsc, nj);
      for (k=0; k < 3; k++)
	R[1][k] = dsc[k] - sp*nj[k];
#endif  
  //printf("scalProd=%.15G\n", scalProd(R[1],R[0]));
  vectProdVec(R[0], R[1], u);
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
void versor_to_R_alt(double *Ci, double *ni, double *Dj, double *nj, double R[3][3], double D)
{
  int k, kk1, kk;
  double u[3]; 
  double normDjCi, DjCi[3], DjCini, Ai[3], AiDjnj, AiDj[3], VV[3], ragg;
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

  /* N.B. qui faccio uno step dell'algorito di Ibarra semplificato per
   * determinare l'asse y del riferimenti del disco (l'asse x è l'asse perpendicolare
   * al disco e l'asse z si ottiene con il prodotto vettore dell'asse x e y) */
#if 1
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  normDjCi = calc_norm(DjCi);
  DjCini = scalProd(DjCi,ni);

  for (kk1 = 0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
  for (kk1=0; kk1 < 3; kk1++)
    AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
  AiDjnj = scalProd(AiDj, nj);
  for (kk1=0; kk1 < 3; kk1++)
    {
      VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
    }
  //for (kk1=0; kk1 < 3; kk1++)
    //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
  ragg = calc_norm(VV);

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k]/ragg;
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
#if 0
  ragg = scalProd(TnCi,ni);
  for (k=0;k<3;k++)
    Ai[k] = Ci[k] + ragg*ni[k];
#endif
#else

 for (k=0; k < 3; k++)
    dsc[k] = Ci[k] - Dj[k]; 
  sp = scalProd(dsc, nj);
  for (k=0; k < 3; k++)
    R[1][k] = dsc[k] - sp*nj[k];
#endif  
  //printf("scalProd=%.15G\n", scalProd(R[1],R[0]));
  vectProdVec(R[0], R[1], u);
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
double test_overlap_parall_cyll(long double *Ci, long double *ni, long double *Dj, long double *nj, long double Li, 
				long double Diami, long double Diamj)
{
  int kk;
  long double DjCi[3], DjCini, Ui[3], DjUi[3], normDjUi;
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  //normDjCi = calc_norm(DjCi);
  DjCini = scalProdl(DjCi,ni);

  for (kk=0; kk < 3; kk++)
    {
      Ui[kk] = Ci[kk] + DjCini*ni[kk];
      DjUi[kk] = Dj[kk] - Ui[kk];
    }
  normDjUi = calc_norml(DjUi);

  if (normDjUi <= 0.5*(Diamj+Diami) && fabsl(DjCini) <= Li*0.5)
    return -1.0;
  else
    return 1.0;
}
double test_overlap_parall_cyl(double *Ci, double *ni, double *Dj, double *nj, double Li, double Diami, double Diamj)
{
  int kk;
  double DjCi[3], DjCini, Ui[3], DjUi[3], normDjUi;
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  //normDjCi = calc_norm(DjCi);
  DjCini = scalProd(DjCi,ni);

  for (kk=0; kk < 3; kk++)
    {
      Ui[kk] = Ci[kk] + DjCini*ni[kk];
      DjUi[kk] = Dj[kk] - Ui[kk];
    }
  normDjUi = calc_norm(DjUi);

  if (normDjUi <= 0.5*(Diamj+Diami) && fabs(DjCini) <= Li*0.5)
    return -1.0;
  else
    return 1.0;
}
/* HYBRID version where for fallback a second ibarra iteration is performed to determine the disk reference system axes  
*/
#ifdef LDLT_LOSSOFSIG_FIX
/* error compensated sums */
double loft[100];
double KahanSum(double input[], int n)
{    
  double sum = 0.0;
  double c = 0.0;
  double y, t;  // A running compensation for lost low-order bits.
  int i;
  for (i = 0; i < n; i++)
    {
      y = input[i] - c; // So far, so good: c is zero.
      t = sum + y;      // Alas, sum is big, y small, so low-order digits of y are lost.
      c = (t - sum) - y;// (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
      sum = t;          // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    }                   // Next time around, the lost low part will be added to y in a fresh attempt.
  return sum;
}

double NeumaierSum(double input[], int n)
{  
  double sum, c, t;
  int i;
  sum = input[0];
  c = 0.0;                 // A running compensation for lost low-order bits.
  for(i = 1; i < n; i++)
    {
      t = sum + input[i];
      if (fabs(sum) >= fabs(input[i]))
	c += (sum - t) + input[i]; // If sum is bigger, low-order digits of input[i] are lost.
      else
	c += (input[i] - t) + sum; // Else low-order digits of sum are lost
      sum = t;
    }
  return sum + c;              // Correction only applied once in the very end
}
#endif
double rimdiskone_hybrid_diff(double Diami, double Li, double Diamj, double Lj, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], nsc, fallback, solset;
#ifdef MC_QUART_VERBOSE
  static long int numfb=0;
#endif
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4], uy[3];
  double dsc[3], c0, c1, c2, c3, c02, c12, c22,  coeffEr[6], c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2];
  double Cip[2][3], nip[2][3], norm, Rl[3][3];
  double nip02,nip12,nip22,Cip02,Cip12,Cip22, temp;
  double omnip02, omnip12, omnip22;
  //double nip04,nip14,nip24,nip03,nip13,nip23;

  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double D2isq, D2i, D2jsq, D2j, Cip0, Cip1, Cip2, nip0, nip1 , nip2, nip1nip2, nip0nip2, nip0nip1; 
  /* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  D2j = Diamj*0.5; 
  D2jsq = Sqr(D2j);
  D2i = Diami*0.5; 
  D2isq = Sqr(D2i);

  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  for (kk1=0; kk1 < 3; kk1++)
    uy[kk1]=Dj[kk1];
  versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, Diamj, uy, 1); 

  //versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[0][kk1] = 0;
      //Aip[kk1] = 0;
      Cip[0][kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[0][kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[0][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip[0]);
  nip0 = nip[0][0]/norm;
  nip1 = nip[0][1]/norm;
  nip2 = nip[0][2]/norm;
  Cip0 = Cip[0][0];
  Cip1 = Cip[0][1];
  Cip2 = Cip[0][2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
#if 0
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
#endif
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);
#if 1
  /* with some simplifications we save a bunch of FLOPS... */
  omnip02 = 1.0 - nip02;
  omnip12 = 1.0 - nip12;
  omnip22 = 1.0 - nip22;
  nip1nip2 = nip1*nip2;
  nip0nip2 = nip0*nip2;
  nip0nip1 = nip0*nip1;
  coeffEr[0] = omnip12;
  coeffEr[1] = omnip22;
  coeffEr[2] = -2.0*nip1nip2;  
  coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
								    Cip1*Cip2*nip1nip2) - D2jsq;
  coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
  coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  
#elif 0
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
				    2*nip1*nip23);

  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;

  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
  /* check ellipse */
#if 0
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  double cq[3], x, lam, solq[2], p[3];
  int numsol;
  lam = -Cip0/nip0;
  x = Cip1 + lam*nip1;
  cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
  cq[1] = coeffEr[5]+coeffEr[2]*x;
  cq[2] = coeffEr[1];
  solve_quadratic(cq, &numsol, solq);
  p[0] = 0.0;
  p[1] = x;
  p[2] = solq[0];
  if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
    {
      printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
      printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
    }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2jsq;
  coeffEr[1] *= D2jsq; 
  coeffEr[2] *= D2jsq;
  coeffEr[4] *= D2j;
  coeffEr[5] *= D2j;
  //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
#ifndef MC_EXCHG_QUART_SOL
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
  coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
  coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
  coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, Li, Diami, Diamj);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
#if 0
  if (numsol==1)
    {
      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
      printf("sol=%.15G\n", solqua[0]);
      printf("BOH\n");
    }
#endif
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
#if 0
      temp = c5 + c2*solqua[kk1];
      solec[kk1][0] = solqua[kk1];
      solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif     
      //printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
#if 0
      if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	{
	  printf("solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	  printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
	  printf("numeratore=%.16G\n", -c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]));
	}
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
      if (temp==0.0) 
	{
	  fallback=1;
	}
    }
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2j*solec[kk1][0];
      solarr[0][kk1][2] = D2j*solec[kk1][1];
      if (test_for_fallbackdiff(solarr[0][kk1], Cip[0], nip[0], D2i, D2j, &(diff[0][kk1])))
	{
	  fallback=1;
#if 0
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[0][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[0][kk1], Cip[0], nip[0]));

	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
    }
  if (tinyimagGBL)
    {
      //printf("qui\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
  if (fallback)
    {
#if defined(MC_QUART_VERBOSE) && 0
      printf("numsol=%d,", numsol[0]);
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  printf("sol=%.16G ", solqua[kk1]);
	}
      printf("\n");
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  temp = c4 + c2*solqua[kk1];
	  printf("temp[%d]=%.16G\n", kk1, temp);
	}
      printf("c2=%.16G c4=%.16G c5=%.16G\n", c2, c4, c5);
      store_bump(iGbl,jGbl);
#endif

#if 0
      for (kk1=0; kk1 < 3; kk1++)
	uy[kk1] = D2*Rl[1][kk1];
      versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 2); 
#else
      rotate_axes_on_plane(Rl);
#endif
      for (kk1=0; kk1 < 3; kk1++)
	{
	  nip[1][kk1] = 0;
	  //Aip[kk1] = 0;
	  Cip[1][kk1] = 0;
	  for (kk2=0; kk2 < 3; kk2++)
	    {
	      nip[1][kk1] += Rl[kk1][kk2]*ni[kk2];
	      Cip[1][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	      //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	    } 
	}
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
      norm = calc_norm(nip[1]);
      nip0 = nip[1][0]/norm;
      nip1 = nip[1][1]/norm;
      nip2 = nip[1][2]/norm;
      Cip0 = Cip[1][0];
      Cip1 = Cip[1][1];
      Cip2 = Cip[1][2];
      nip02=Sqr(nip0);
      nip12=Sqr(nip1);
      nip22=Sqr(nip2);
      //nip04=Sqr(nip02);
      //nip14=Sqr(nip12);
      //nip24=Sqr(nip22);
      //nip03=nip02*nip0;
      //nip13=nip12*nip1;
      //nip23=nip22*nip2;
      Cip02=Sqr(Cip0);
      Cip12=Sqr(Cip1);
      Cip22=Sqr(Cip2);   
#if 1
      omnip02 = 1.0 - nip02;
      omnip12 = 1.0 - nip12;
      omnip22 = 1.0 - nip22;
      nip1nip2 = nip1*nip2;
      nip0nip2 = nip0*nip2;
      nip0nip1 = nip0*nip1;
      coeffEr[0] = omnip12;
      coeffEr[1] = omnip22;
      coeffEr[2] = -2.0*nip1nip2;  
      coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
									Cip1*Cip2*nip1nip2) - D2jsq;
      coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
      coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  

#elif 0
      coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
	nip12*nip22;
      coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
	nip12*nip22 + nip24;
      coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
	2*nip1*nip23;
      coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
	2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
	2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
	4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
	2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
	Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
	2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
      coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
	4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
	2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;
      coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
	4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
	- 2*Cip2*nip24;
#else
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

      coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
      coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
      coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
					2*nip1*nip23);

      coeffEr[3] = 
	(- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
	 Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
	 + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
	 Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
	+ Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
	2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
	+ Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

      coeffEr[4] =
	(4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
	 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
	 - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
	- 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;

      coeffEr[5] = 
	(4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
	 - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
      /* check ellipse */

      /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
      coeffEr[0] *= D2jsq;
      coeffEr[1] *= D2jsq; 
      coeffEr[2] *= D2jsq;
      coeffEr[4] *= D2j;
      coeffEr[5] *= D2j;
      //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
      c0 = coeffEr[0];
      c1 = coeffEr[1];
      c2 = coeffEr[2];
      c3 = coeffEr[3];
      c4 = coeffEr[4];
      c5 = coeffEr[5];
      c02 = Sqr(c0);
      c12 = Sqr(c1);
      c22 = Sqr(c2);
      c32 = Sqr(c3);
      c42 = Sqr(c4);
      c52 = Sqr(c5);

#if 0
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, Li, Diami, Diamj);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));
#ifdef MC_QUART_VERBOSE
      printf("falling back [#%ld] type=%d numsol=%d %d\n", numfb++,fallback, numsol[0], numsol[1]);
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{

#if 0 
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif
	  temp = c4 + c2*solqua[kk1];
	  //printf("tempnew=%.16G\n", temp);
	  if (temp==0)
	    {
	      printf("[WARNING] temp is 0 in fallback hybrid numsol=%d %d\n", numsol[0], numsol[1]);
	    }
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2j*solec[kk1][0];
	  solarr[1][kk1][2] = D2j*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
	  test_for_fallbackdiff(solarr[1][kk1], Cip[1], nip[1], D2i, D2j, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  
	}
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	    //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
	}
    }
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[solset][kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if !defined(MC_ELASTIC_CONSTANTS) && !defined(MC_ELCONST_MC)
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallbackdiff(solarr[solset][kk1], Cip[solset], nip[solset], D2i, D2j, &tmp)) 
	{
	  printf("# %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip[solset], nip[solset]));

	  if (kk1 < numsol[1-solset])
	    {
	      printf("DISCARDED SOLSET [%d]\n", 1-solset);
	      printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[1-solset][kk1]));
	      printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip[1-solset], nip[1-solset]));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[solset][1], Cip[solset][2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
#endif
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProd(dsc, nip[solset]);
      if (fabs(sp) < Li*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
double rimdiskone_hybrid(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], fallback, solset;
#ifdef MC_QUART_VERBOSE
  static long int numfb=0;
#endif
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4], uy[3];
  double dsc[3], c0, c1, c2, c3, c02, c12, c22, coeffEr[6], c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2];
  double Cip[2][3], nip[2][3], norm, Rl[3][3];
  double nip02,nip12,nip22,Cip02,Cip12,Cip22, temp;
  double omnip02, omnip12, omnip22;
  //double nip04,nip14,nip24,nip03,nip13,nip23;

  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2, nip1nip2, nip0nip2, nip0nip1; 
  /* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  D2 = D*0.5; 
  D2sq = Sqr(D2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  for (kk1=0; kk1 < 3; kk1++)
    uy[kk1]=Dj[kk1];
  versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 1); 

  //versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[0][kk1] = 0;
      //Aip[kk1] = 0;
      Cip[0][kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[0][kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[0][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip[0]);
  nip0 = nip[0][0]/norm;
  nip1 = nip[0][1]/norm;
  nip2 = nip[0][2]/norm;
  Cip0 = Cip[0][0];
  Cip1 = Cip[0][1];
  Cip2 = Cip[0][2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
#if 0
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
#endif
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);
#if 1
  /* with some simplifications we save a bunch of FLOPS... */
  omnip02 = 1.0 - nip02;
  omnip12 = 1.0 - nip12;
  omnip22 = 1.0 - nip22;
  nip1nip2 = nip1*nip2;
  nip0nip2 = nip0*nip2;
  nip0nip1 = nip0*nip1;
#ifdef LDLT_LOSSOFSIG_FIX
  coeffEr[0] = omnip12;
  coeffEr[1] = omnip22;
  coeffEr[2] = -2.0*nip1nip2;
  loft[0] =  Cip02*omnip02;
  loft[1] = Cip12*omnip12;
  loft[2] = Cip22*omnip22; 
  loft[3] = -2.0*Cip0*Cip1*nip0nip1; 
  loft[4] = -2.0*Cip0*Cip2*nip0nip2;
  loft[5] = -2.0* Cip1*Cip2*nip1nip2;
  loft[6] = -D2sq;
  coeffEr[3]= NeumaierSum(loft, 7);
  loft[0] = Cip2*nip1nip2;
  loft[1] = Cip0*nip0nip1; 
  loft[2] = -Cip1*omnip12;
  coeffEr[4] = 2.0*NeumaierSum(loft,3);
  loft[0] = Cip0*nip0nip2; 
  loft[1] = Cip1*nip1nip2;
  loft[2] = - Cip2*omnip22;
  coeffEr[5] = 2.0*NeumaierSum(loft,3);
#else
  coeffEr[0] = omnip12;
  coeffEr[1] = omnip22;
  coeffEr[2] = -2.0*nip1nip2;  
  coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
								    Cip1*Cip2*nip1nip2) - D2sq;
  coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
  coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  
#endif
#elif 0
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
				    2*nip1*nip23);

  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;

  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
  /* check ellipse */
#if 0
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  double cq[3], x, lam, solq[2], p[3];
  int numsol;
  lam = -Cip0/nip0;
  x = Cip1 + lam*nip1;
  cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
  cq[1] = coeffEr[5]+coeffEr[2]*x;
  cq[2] = coeffEr[1];
  solve_quadratic(cq, &numsol, solq);
  p[0] = 0.0;
  p[1] = x;
  p[2] = solq[0];
  if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
    {
      printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
      printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
    }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2sq;
  coeffEr[1] *= D2sq; 
  coeffEr[2] *= D2sq;
  coeffEr[4] *= D2;
  coeffEr[5] *= D2;
  //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
#ifndef MC_EXCHG_QUART_SOL
#ifdef LDLT_LOSSOFSIG_FIX
  loft[0] = c02; 
  loft[1] = -2*c0*c1;
  loft[2] = c12;
  loft[3] = c22;
  coeff[4] = NeumaierSum(loft,4);
  loft[0] = 2*c2*c4;
  loft[1] = - 2*c0*c5;
  loft[2] = 2*c1*c5;
  coeff[3] = NeumaierSum(loft,3);
  loft[0] = -2*c02;
  loft[1] = 2*c0*c1;
  loft[2] = -c22;
  loft[3] = - 2*c0*c3;
  loft[4] = 2*c1*c3;
  loft[5] = c42; 
  loft[6] = c52;
  coeff[2] = NeumaierSum(loft,7);
  loft[0] = -2*c2*c4; 
  loft[1] = 2*c0*c5;
  loft[2] = 2*c3*c5;
  coeff[1] = NeumaierSum(loft,3);
  loft[0] = c02;
  loft[1] = 2*c0*c3;
  loft[2] = c32;
  loft[3] = -c42;
  coeff[0] = NeumaierSum(loft,4);
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#endif
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
  coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
  coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
  coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
#if 0
  if (numsol==1)
    {
      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
      printf("sol=%.15G\n", solqua[0]);
      printf("BOH\n");
    }
#endif
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      temp = c4 + c2*solqua[kk1];
#ifdef LDLT_LOSSOFSIG_FIX
      loft[0] = -c0;
      loft[1] = -c3;
      loft[2] = -c5*solqua[kk1];
      loft[3] = c0*Sqr(solqua[kk1]);
      loft[4] =-c1*Sqr(solqua[kk1]);
      solec[kk1][0] = NeumaierSum(loft,5)/temp;
#else
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
#endif
      solec[kk1][1] = solqua[kk1];
#if 0
      temp = c5 + c2*solqua[kk1];
      solec[kk1][0] = solqua[kk1];
      solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif     
      //printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
#if 0
      if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	{
	  printf("solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	  printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
	  printf("numeratore=%.16G\n", -c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]));
	}
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
      if (temp==0.0) 
	{
	  fallback=1;
	}
    }
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2*solec[kk1][0];
      solarr[0][kk1][2] = D2*solec[kk1][1];
      if (test_for_fallback(solarr[0][kk1], Cip[0], nip[0], D2, &(diff[0][kk1])))
	{
	  fallback=1;
#if 0
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[0][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[0][kk1], Cip[0], nip[0]));

	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
    }
  if (tinyimagGBL)
    {
      //printf("qui\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
  if (fallback)
    {
#if defined(MC_QUART_VERBOSE) && 0
      printf("numsol=%d,", numsol[0]);
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  printf("sol=%.16G ", solqua[kk1]);
	}
      printf("\n");
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  temp = c4 + c2*solqua[kk1];
	  printf("temp[%d]=%.16G\n", kk1, temp);
	}
      printf("c2=%.16G c4=%.16G c5=%.16G\n", c2, c4, c5);
      store_bump(iGbl,jGbl);
#endif

#if 0
      for (kk1=0; kk1 < 3; kk1++)
	uy[kk1] = D2*Rl[1][kk1];
      versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 2); 
#else
      rotate_axes_on_plane(Rl);
#endif
      for (kk1=0; kk1 < 3; kk1++)
	{
	  nip[1][kk1] = 0;
	  //Aip[kk1] = 0;
	  Cip[1][kk1] = 0;
	  for (kk2=0; kk2 < 3; kk2++)
	    {
	      nip[1][kk1] += Rl[kk1][kk2]*ni[kk2];
	      Cip[1][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	      //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	    } 
	}
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
      norm = calc_norm(nip[1]);
      nip0 = nip[1][0]/norm;
      nip1 = nip[1][1]/norm;
      nip2 = nip[1][2]/norm;
      Cip0 = Cip[1][0];
      Cip1 = Cip[1][1];
      Cip2 = Cip[1][2];
      nip02=Sqr(nip0);
      nip12=Sqr(nip1);
      nip22=Sqr(nip2);
      //nip04=Sqr(nip02);
      //nip14=Sqr(nip12);
      //nip24=Sqr(nip22);
      //nip03=nip02*nip0;
      //nip13=nip12*nip1;
      //nip23=nip22*nip2;
      Cip02=Sqr(Cip0);
      Cip12=Sqr(Cip1);
      Cip22=Sqr(Cip2);   
#if 1
      omnip02 = 1.0 - nip02;
      omnip12 = 1.0 - nip12;
      omnip22 = 1.0 - nip22;
      nip1nip2 = nip1*nip2;
      nip0nip2 = nip0*nip2;
      nip0nip1 = nip0*nip1;
      coeffEr[0] = omnip12;
      coeffEr[1] = omnip22;
      coeffEr[2] = -2.0*nip1nip2;  
      coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
									Cip1*Cip2*nip1nip2) - D2sq;
      coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
      coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  

#elif 0
      coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
	nip12*nip22;
      coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
	nip12*nip22 + nip24;
      coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
	2*nip1*nip23;
      coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
	2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
	2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
	4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
	2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
	Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
	2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
      coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
	4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
	2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;
      coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
	4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
	- 2*Cip2*nip24;
#else
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

      coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
      coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
      coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
					2*nip1*nip23);

      coeffEr[3] = 
	(- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
	 Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
	 + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
	 Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
	+ Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
	2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
	+ Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

      coeffEr[4] =
	(4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
	 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
	 - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
	- 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;

      coeffEr[5] = 
	(4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
	 - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
      /* check ellipse */

      /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
      coeffEr[0] *= D2sq;
      coeffEr[1] *= D2sq; 
      coeffEr[2] *= D2sq;
      coeffEr[4] *= D2;
      coeffEr[5] *= D2;
      //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
      c0 = coeffEr[0];
      c1 = coeffEr[1];
      c2 = coeffEr[2];
      c3 = coeffEr[3];
      c4 = coeffEr[4];
      c5 = coeffEr[5];
      c02 = Sqr(c0);
      c12 = Sqr(c1);
      c22 = Sqr(c2);
      c32 = Sqr(c3);
      c42 = Sqr(c4);
      c52 = Sqr(c5);

#if 0
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));
#ifdef MC_QUART_VERBOSE
      printf("falling back [#%ld] type=%d numsol=%d %d\n", numfb++,fallback, numsol[0], numsol[1]);
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{

#if 0 
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif
	  temp = c4 + c2*solqua[kk1];
	  //printf("tempnew=%.16G\n", temp);
	  if (temp==0)
	    {
	      printf("[WARNING] temp is 0 in fallback hybrid numsol=%d %d\n", numsol[0], numsol[1]);
	    }
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2*solec[kk1][0];
	  solarr[1][kk1][2] = D2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
	  test_for_fallback(solarr[1][kk1], Cip[1], nip[1], D2, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  
	}
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	    //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
	}
    }
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[solset][kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if !defined(MC_ELASTIC_CONSTANTS) && !defined(MC_ELCONST_MC)
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallback(solarr[solset][kk1], Cip[solset], nip[solset], D2, &tmp)) 
	{
	  printf("# %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip[solset], nip[solset]));

	  if (kk1 < numsol[1-solset])
	    {
	      printf("DISCARDED SOLSET [%d]\n", 1-solset);
	      printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[1-solset][kk1]));
	      printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip[1-solset], nip[1-solset]));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[solset][1], Cip[solset][2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
#endif
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProd(dsc, nip[solset]);
      if (fabs(sp) < L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
double rimdiskone_solvxy(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], fallback, solset;
#ifdef MC_QUART_VERBOSE
  static long int numfb=0;
#endif
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4];
  double dsc[3], c0, c1, c2, c3, c02, c12, c22, coeffEr[6], c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2];
  double Cip[3], nip[3], norm, Rl[3][3];
  double nip02,nip12,nip22,nip03,nip13,nip23,nip04,nip14,nip24,Cip02,Cip12,Cip22, temp;
  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2; 
/* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  D2 = D*0.5; 
  D2sq = Sqr(D2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[kk1] = 0;
      //Aip[kk1] = 0;
      Cip[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip);
  nip0 = nip[0]/norm;
  nip1 = nip[1]/norm;
  nip2 = nip[2]/norm;
  Cip0 = Cip[0];
  Cip1 = Cip[1];
  Cip2 = Cip[2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);   
#if 1
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
 /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23);
 
  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
 
  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
 
  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;
 
#endif
  /* check ellipse */
#if 0
  {
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
    double cq[3], x, lam, solq[2], p[3];
    int numsol;
    lam = -Cip0/nip0;
    x = Cip1 + lam*nip1;
    cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
    cq[1] = coeffEr[5]+coeffEr[2]*x;
    cq[2] = coeffEr[1];
    solve_quadratic(cq, &numsol, solq);
    p[0] = 0.0;
    p[1] = x;
    p[2] = solq[0];
    if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
      {
	printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      	printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
	printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
      }
  }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2sq;
  coeffEr[1] *= D2sq; 
  coeffEr[2] *= D2sq;
  coeffEr[4] *= D2;
  coeffEr[5] *= D2;
  //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
#ifndef MC_EXCHG_QUART_SOL
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
  coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
  coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
  coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
#if 0
  if (numsol==1)
    {
      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
      printf("sol=%.15G\n", solqua[0]);
      printf("BOH\n");
    }
#endif
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
#if 0
      c0l = c0;
      c1l = c1;
      c2l = c2;
      c3l = c3;
      c4l = c4;
      c5l = c5;
      solqual=solqua[kk1];
      templ = c4l + c2l*solqual;
      solec[kk1][0]=((double)((-c0l - c3l - c5l*solqual + (c0l - c1l)*Sqr(solqual))/templ));
      solec[kk1][1] = solqua[kk1];
#else
#ifndef MC_EXCHG_QUART_SOL
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
#else
      temp = c5 + c2*solqua[kk1];
      solec[kk1][0] = solqua[kk1];
      solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif     
 //printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
#if 0
      if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	{
  	  printf("solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	  printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
	  printf("numeratore=%.16G\n", -c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]));
	}
#endif
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
#if 0
      if (test_solution_xy(solec[kk1], &(diffxy[0][kk1])))
	  fallback=1;
#endif
      if (temp==0.0) 
	{
	  printf("[WARNING] temp is 0 in hybrid numsol=%d\n", numsol[0]);
	  fallback=1;
	}
    }
#if 0
      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
	     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
	     coeff[1]*solqua[kk1]+coeff[0]);
      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
#if 0
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      solve_quartic(coeff, &numsol, solqua);
      for (kk1=0; kk1 < numsol; kk1++)
	{
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	}
    }
#endif
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2*solec[kk1][0];
      solarr[0][kk1][2] = D2*solec[kk1][1];
#if 1
      if (test_for_fallback(solarr[0][kk1], Cip, nip, D2, &(diff[0][kk1])))
	{
	  fallback=1;
#if 0
	  if (numsol==4)
	    {
	      printf("%d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		     coeff[1]*solqua[kk1]+coeff[0]);
	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[0][kk1]);
	      printf(">>>> =====================\n");
	    }
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
#endif
    }
  if (tinyimagGBL)
    {
      //printf("BOHHHH\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
#if 1
  if (fallback)
    {
#if defined(MC_QUART_VERBOSE) && 0
      printf("numsol=%d,", numsol[0]);
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  printf("sol=%.16G ", solqua[kk1]);
	}
      printf("\n");
      printf("temp=%.16G\n", temp);
      store_bump(iGbl,jGbl);
#endif

#ifndef MC_EXCHG_QUART_SOL
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#else
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#endif
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));
#ifdef MC_QUART_VERBOSE
      printf("falling back [#%ld] type=%d numsol=%d %d\n", numfb++,fallback, numsol[0], numsol[1]);
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0
    	  c0l = c0;
	  c1l = c1;
	  c2l = c2;
	  c3l = c3;
	  c4l = c4;
	  c5l = c5;
	  solqual=solqua[kk1];
	  templ = c5l + c2l*solqual;
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = ((double)((-c1l - c3l - c4l*solqual + (c1l - c0l)*Sqr(solqual))/templ)); 
#else
#ifndef MC_EXCHG_QUART_SOL
	  temp = c5 + c2*solqua[kk1];
	  //printf("tempnew=%.16G\n", temp);
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#else
    	  temp = c4 + c2*solqua[kk1];
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
#endif
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2*solec[kk1][0];
	  solarr[1][kk1][2] = D2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
#if 1
	  test_for_fallback(solarr[1][kk1], Cip, nip, D2, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  

#endif
#if 0
  	  if (numsol==4)
  	    {
  	      printf("FALLBACK %d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
  	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
  	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
  	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
  	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
  	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
  		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
  		     coeff[1]*solqua[kk1]+coeff[0]);
  	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[1][kk1]);
	      printf(">>>> =====================\n");
	  }
#endif
	
	}
#if 1
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	  //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
#if 0
	  if (fallback==3 && solset != 1)
	    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
	}
#endif
    }
#endif
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
#if 0
  construct_inner_points(solarr, Ci, ni, Dj, nj, D);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
#if 0
      printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
      printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 0
      if (calc_normsq(solarr[kk1])-Sqr(D2) > NEWT_THR)
	{
	  newt2Dquartic(coeffEr, solarr[kk1], D2);
	}
#endif
#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallback(solarr[solset][kk1], Cip, nip, D2, &tmp)) 
	{
	  printf("# %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip, nip));
	  
	  if (kk1 < numsol[1-solset])
	    {
	      printf("DISCARDED SOLSET [%d]\n", 1-solset);
    	      printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[1-solset][kk1]));
    	      printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip, nip));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[1], Cip[2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.14G\n", scalProd(ni,nj));
	  printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
#endif
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProd(dsc, nip);
      if (fabs(sp) < L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
double rimdiskone(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
#ifdef MC_QUART_HYBRID
  return rimdiskone_hybrid(D, L, Ci, ni, Dj, nj, DjCini);
#else
  return rimdiskone_solvxy(D, L, Ci, ni, Dj, nj, DjCini);
#endif

}
double rimdiskdiff(double *D, double *L, double Ci[3], double ni[3], double Di[2][3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, kk, j2, k2, ignore[2];
  double DjUini, DjTmp[2][3], DjCi[3], DjUi[3], niTmp[3], njTmp[3], perpdist[2];
  double normDjUi, normDjCi, DjCini, Ui[3], CiTmp[3];
  double LiTmp, LjTmp, DiamiTmp, DiamjTmp, Li, Lj, Diami, Diamj;
  Diami=D[0];
  Diamj=D[1];
  Li=L[0];
  Lj=L[1];
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
      
      for (k2=0; k2 < 2; k2++)
	{
  	  perpdist[k2]=perpcomp(Dj[k2], Ci, ni);
	  ignore[k2] = 0;
	}

      if (perpdist[0] < perpdist[1])
	{ 
    	  ignore[1] = 1;
	}
      else
	{
	  ignore[0] = 1;
	}

      for (j2=0; j2 < 2; j2++)
	{
	  if (ignore[j2])
	    continue;
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
	  if (normDjUi > 0.5*(Diami+Diamj))
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Ui| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

*/
#ifndef MC_IBARRA_SIMPLER
	  /* se sono quasi paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) < 1.0E-8)
	    {
	      if (normDjUi <= 0.5*(Diami+Diamj) && fabs(DjCini) <= Li*0.5)
		return -1;
	      else
		continue;
	    }
#else
	  /* se sono quasi paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) < 3E-16)
	    {
	      if (normDjUi <= 0.5*(Diami+Diamj) && fabs(DjCini) <= Li*0.5)
		return -1;
	      else
		continue;
	    }
#endif
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

#ifdef MC_IBARRA_SIMPLER
  	  if (rimdiskone_ibarradiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;
#else
#ifdef MC_QUART_LONG_DOUBLE
#ifdef MC_DEBUG_HCALGO
	    {
	      double alg1, alg2; 
	      alg1 = rimdiskone_ibarradiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      alg2 = rimdiskoneldiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
		  //printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
	      if (alg1!=alg2)
		{
		  store_bump(iGbl, jGbl);
		  printf("Discrepancy between i=%d and j=%d!!\n", iGbl, jGbl);
		  printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
		  //saveCorAscii();
		  sprintf(fileop2 ,"coord-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveCoord(fileop2);
		  sprintf(fileop2 ,"coorbakascii-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveBakAscii(fileop2);
		  //exit(-1);
		}
	      if (alg2 < 0)
		return -1;
	    }
#else
	  /* N.B. NON ANCORA IMPLEMENTATA */
	  if (rimdiskonediffl(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;
#endif
#else
#ifdef MC_DEBUG_HCALGO
	    {
	      double alg1, alg2; 
	      alg1 = rimdiskone_ibarradiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      alg2 = rimdiskonediff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
		  //printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
	      if (alg1!=alg2)
		{
		  store_bump(iGbl, jGbl);
		  printf("Discrepancy between i=%d and j=%d!!\n", iGbl, jGbl);
		  printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
		  //saveCorAscii();
		  sprintf(fileop2 ,"coord-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveCoord(fileop2);
		  sprintf(fileop2 ,"coorbakascii-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveBakAscii(fileop2);
		  //exit(-1);
		}
	      if (alg2 < 0)
		return -1;
	    }
#else
#ifdef MC_QUART_HYBRID
         if (rimdiskone_hybrid_diff(Diami, Li, Diamj, Lj, Ci, ni, Dj, nj, DjCini) < 0.0)
           return -1;
#else
         if (rimdiskonediff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
           return -1;
#endif
#endif
#endif
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
  return 0;
}
double rimdisk(double D, double L, double Ci[3], double ni[3], double Di[2][3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, kk, j2, k2, ignore[2];
  double DjUini, DjTmp[2][3], DjCi[3], DjUi[3], niTmp[3], njTmp[3], perpdist[2];
  double normDjUi, normDjCi, DjCini, Ui[3], CiTmp[3];

  numcallsRD++;

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
      
      for (k2=0; k2 < 2; k2++)
	{
  	  perpdist[k2]=perpcomp(Dj[k2], Ci, ni);
	  ignore[k2] = 0;
	}

      if (perpdist[0] < perpdist[1])
	{ 
    	  ignore[1] = 1;
	}
      else
	{
	  ignore[0] = 1;
	}

      for (j2=0; j2 < 2; j2++)
	{
	  if (ignore[j2])
	    continue;
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
#ifndef MC_IBARRA_SIMPLER
	  /* se sono paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) == 0)
	  //if (1.0-fabs(scalProd(ni,nj)) < 1.0E-8)
	    {
	      if (normDjUi <= D && fabs(DjCini) <= L*0.5)
		return -1;
	      else
		continue;
	    }
#else
	  /* se sono paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) == 0)
	  //if (1.0-fabs(scalProd(ni,nj)) < 3E-16)
	    {
	      if (normDjUi <= D && fabs(DjCini) <= L*0.5)
		return -1;
	      else
		continue;
	    }
#endif
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

#ifdef MC_IBARRA_SIMPLER
  	  if (rimdiskone_ibarra(D, L, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;
#else
#ifdef MC_QUART_LONG_DOUBLE
#ifdef MC_DEBUG_HCALGO
	    {
	      double alg1, alg2; 
	      alg1 = rimdiskone_ibarra(D, L, Ci, ni, Dj[j2], nj, DjCini);
	      alg2 = rimdiskonel(D, L, Ci, ni, Dj[j2], nj, DjCini);
	      //alg1=alg2;
	      //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
		  //printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
	      if (alg1!=alg2)
		{
		  store_bump(iGbl, jGbl);
		  printf("Discrepancy between i=%d and j=%d!!\n", iGbl, jGbl);
		  printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
		  //saveCorAscii();
		  sprintf(fileop2 ,"coord-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveCoord(fileop2);
		  sprintf(fileop2 ,"coorbakascii-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveBakAscii(fileop2);
		  //exit(-1);
		}
	      if (alg2 < 0)
		return -1;
	    }
#else
	  if (rimdiskonel(D, L, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;
#endif
#else
#ifdef MC_DEBUG_HCALGO
	    {
	      double alg1, alg2; 
	      alg1 = rimdiskone_ibarra(D, L, Ci, ni, Dj[j2], nj, DjCini);
	      alg2 = rimdiskone(D, L, Ci, ni, Dj[j2], nj, DjCini);
	      //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
		  //printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
	      if (alg1!=alg2)
		{
		  store_bump(iGbl, jGbl);
		  printf("Discrepancy between i=%d and j=%d!!\n", iGbl, jGbl);
		  printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
		  //saveCorAscii();
		  sprintf(fileop2 ,"coord-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveCoord(fileop2);
		  sprintf(fileop2 ,"coorbakascii-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveBakAscii(fileop2);
		  //exit(-1);
		}
	      if (alg2 < 0)
		return -1;
	    }
#else
	  if (rimdiskone(D, L, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;

#endif
#endif
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
  return 0;
}
double rimrimdiff(double *D, double *L, double Ci[3],double ni[3], double Cj[3], double nj[3])
{
  int kk;
  double ViVj[3], lambdai, lambdaj, ninj;
  double CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3]; 
  /* case A.3 rim-rim overlap */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
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
  if (calc_norm(ViVj) < 0.5*(D[0]+D[1]) && fabs(lambdai) < 0.5*L[0] && fabs(lambdaj) < 0.5*L[1])
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 0;
}
double rimrim(double D, double L, double Ci[3],double ni[3], double Cj[3], double nj[3])
{
  int kk;
  double ViVj[3], lambdai, lambdaj, ninj;
  double CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3]; 
  /* case A.3 rim-rim overlap */
  numcallsRR++;
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
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
  if ( //Sqr(ViVj[0])+Sqr(ViVj[1])+Sqr(ViVj[2])<Sqr(D) 
       calc_norm(ViVj) < D 
       && fabs(lambdai) < 0.5*L && fabs(lambdaj) < 0.5*L)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 0;
}
#if 0
double calcDistNegHCbrentDiff(int i, int j, double shift[3], int* retchk)
{
  static int firstcall=1;
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int kk;
  double Ci[3], Cj[3], L, D, Di[2][3], Dj[2][3], ni[3], nj[3], ret;
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
#ifdef MC_DEBUG_HCALGO
  iGbl=i;
  jGbl=j;
#endif
  if (L >= D) // prolate
    {
      if ((ret=rimrim(D, L, Ci, ni, Cj, nj)) != 0.0)
	return ret;

      if ((ret=diskdisk(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimdisk(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
	return ret;
     }
  else // oblate
    {
      if ((ret=diskdisk(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimrim(D, L, Ci, ni, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimdisk(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
	return ret;
    }
  /* case A.2 overlap of rim and disk */
  /* =================================== >>> Part A <<< ========================= */
 /* =================================== >>> Part B <<< ========================= */
  numcallsHC += 4.0; 
  return 1;
}
#endif
double calcDistNegHCoptdiff(int i, int j, double shift[3], int* retchk)
{
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int kk;
  double Ci[3], Cj[3], L[2], D[2], Di[2][3], Dj[2][3], ni[3], nj[3], ret;

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
  L[0] = 2.0*typesArr[typeOfPart[i]].sax[0];
  D[0] = 2.0*typesArr[typeOfPart[i]].sax[1];
  L[1] = 2.0*typesArr[typeOfPart[j]].sax[0];
  D[1] = 2.0*typesArr[typeOfPart[j]].sax[1];
  
  meshptsGbl = MESH_PTS;
  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*L[0]*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*L[0]*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*L[1]*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*L[1]*nj[kk];
    }
#if 0
#ifdef MC_HC_SPHERO_OPT
  if ((sphov=check_spherocyl(CiCj, D, L, Di, Ci, ni, Dj, Cj, nj, &rim)) > 0.0)
    return 1;
#endif
#endif
#ifdef MC_DEBUG_HCALGO
  iGbl=i;
  jGbl=j;
#endif
  if ((ret=rimrimdiff(D, L, Ci, ni, Cj, nj)) != 0.0)
    return ret;

  if ((ret=diskdiskdiff(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
    return ret;

  if ((ret=rimdiskdiff(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
    return ret;
  /* case A.2 overlap of rim and disk */
  /* =================================== >>> Part A <<< ========================= */
 /* =================================== >>> Part B <<< ========================= */
  numcallsHC += 4.0; 
  return 1;
}
double calcDistNegHCopt(int i, int j, double shift[3], int* retchk)
{
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int kk;
  double Ci[3], Cj[3], L, D, Di[2][3], Dj[2][3], ni[3], nj[3], ret;

  if (typesArr[typeOfPart[i]].sax[0]!=typesArr[typeOfPart[j]].sax[0]
      || typesArr[typeOfPart[i]].sax[1] != typesArr[typeOfPart[j]].sax[1])
    return calcDistNegHCoptdiff(i, j, shift, retchk);

  *retchk = 0; 

  numcallsHC += 1.0; 
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
#ifdef MC_DEBUG_HCALGO
  iGbl=i;
  jGbl=j;
#endif
  if (L >= D) // prolate
    {
      if ((ret=diskdisk(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimrim(D, L, Ci, ni, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimdisk(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
	return ret;
     }
  else // oblate
    {
      if ((ret=diskdisk(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimrim(D, L, Ci, ni, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimdisk(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
	return ret;
    }
  /* case A.2 overlap of rim and disk */
  /* =================================== >>> Part A <<< ========================= */
 /* =================================== >>> Part B <<< ========================= */
  return 1;
}
#if 0
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
#endif
#ifdef MC_SPHEROCYL
double check_spherocyl(double CiCj[3], double D, double Lc, double Di[2][3], double *Ci, double *ni, double Dj[2][3], double *Cj, double *nj, int *rim);

double calcDistSpheroCyl(int i, int j, double shift[3])
{
  int k, kk, rim;
  double dist, CiCj[3], D, L, Di[2][3], Ci[3], ni[3], Dj[2][3], Cj[3], nj[3];

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

  dist = check_spherocyl(CiCj, D, L, Di, Ci, ni, Dj, Cj, nj, &rim);
  //printf("dist=%.15G\n", dist);
  return dist;
}
#endif
double calcDistNegHC(int i, int j, double shift[3], int* retchk)
{
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int it, k2, kk1, kk2, k;
  double ragg, ragg2, norm, normNSq, ViVj[3], lambdai, lambdaj, Aiold[3], TnCi[3], Tnew[3], Told[3];
  double sp, Q1, Q2, normPiDi, normPjDj, normN, L, D, DiN, DjN, niN[3], njN[3], Djni, Djnj, AiOld[3];
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
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3], dsc[3], dscperp[3], dscpara[3];
  int kk, j1, j2;

#ifdef HC_ALGO_OPT
   if (typesArr[typeOfPart[i]].sax[0]!=typesArr[typeOfPart[j]].sax[0]
      || typesArr[typeOfPart[i]].sax[1] != typesArr[typeOfPart[j]].sax[1])
     return calcDistNegHCoptdiff(i, j, shift, retchk);
   else
     return calcDistNegHCopt(i, j, shift, retchk);
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
#if 0
	  find_initial_guess(Ai, Ci, ni, Dj[j2], nj, D);

#else
	  find_initial_guess_simpler(Ai, Ci, ni, Dj[j2], nj, D);
 //	  for (kk=0; kk < 3; kk++)
//	    {
	      //Ai[kk] = Ci[kk];
//	      Ai[kk] = Ui[kk];  
//	    }
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
	      //printf("rimdisk overlap found i=%d j=%d\n", i, j);
	      return -1;
	    }
	  totitsHC += it;
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
#if defined(MC_HC_SPHERO_OPT) || defined(MC_SPHEROCYL)
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
#if 1
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

 if(cc<1e-6) {
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
#endif
double rimrim_sphcyl(double D, double L, double Ci[3],double ni[3], double Cj[3], double nj[3])
{
  int kk;
  double ViVj[3], lambdai, lambdaj, ninj;
  double CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3]; 
  /* case A.3 rim-rim overlap */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
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
double check_spherocyl(double CiCj[3], double D, double Lc, double Di[2][3], double *Ci, double *ni, double Dj[2][3], double *Cj, double *nj, int *rim)
{
  coo_t r1, r2, w1, w2;
  double sum, d, normDiCj, normDjCi, DiCj[3], DjCi[3], Ui[3], Uj[3], DjUi[3], DiUj[3], DjCini, DiCjnj;
  int kk, j1, j2;

#if 1
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
#endif
#if 0
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
#if 1
  if  (rimrim_sphcyl(D, Lc, Ci, ni, Cj, nj) < 0)
    {
      //printf("BOHBOH\n");
      return -1;
    } 
#else
  if ((d=dist2_rods(r1, r2, w1, w2, Lc*0.5, Lc*0.5)) <= Sqr(D)) 
    {
      *rim=-1;
      return -1;
    }
#endif
#else
  // dist2_rods da la distanza tra due segmenti quindi è quello che serve per capire se gli sferocilindri
  // si sovrappongono o meno
  if ((d=dist2_rods(r1, r2, w1, w2, Lc*0.5, Lc*0.5)) <= Sqr(D)) 
    {
      *rim=-1;
      return -1;
    }
#endif
  
  return 1;
}
#endif
#endif
#endif
/* === discareded stuff if calcDistNegHCdiff */
#if 0
	  if (scalProd(ni,nj)==1)
	    {
	      printf("per ora esco...\n");
	      exit(-1);
	    }
	  else
	    {
	      /* =================================================================================== */
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
#if 1
	      /* mi metto nel riferimento del disco (p) */
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
#endif
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
#if 0

	      for (kk1=0; kk1 < 3; kk1++)
		{
		  CipGbl[kk1] = Cip[kk1];
		  nipGbl[kk1] = nip[kk1];
		}
	      Dgbl = D;
#endif
	      D2 = D*0.5;
	      /* calcolo i semiassi ed il centro dell'ellisse che si ottiene tagliando il rim con il piano del disco */
	      lambda = -Cip[0]/nip[0];
	      /* centro dell'ellisse del rim */
	      rErp[0] = 0;
	      rErp[1] = Cip[1]+lambda*nip[1];
	      rErp[2] = Cip[2]+lambda*nip[2];

	      nErxp[0] = 1.0;
	      nErxp[1] = 0.0;
	      nErxp[2] = 0.0;
	      nEryp[0]=0.0;
	      nEryp[2]=1.0/sqrt(1.0+Sqr(nip[2]/nip[1]));
	      nEryp[1]=-nEryp[2]*nip[2]/nip[1];
	      vectProdVec(nErxp,nEryp,nErzp);
	      aErcut=D2;	
	      delta = 1.0-Sqr(scalProd(nErzp,nip));
	      if (delta < 0)// se è < 0 è soltanto per errori di roundoff 
		delta = 0.0;
	      bErcut=D2/sqrt(delta);
#if 0
	      printf("scalprod nErx.nip=%.15G\n", scalProd(nErzp,nip)); 
	      printf("nip=%.15G %.15G %.15G\n", nip[0], nip[1], nip[2]);
	      printf("Erzp=%.15G %.15G %.15G\n", nErzp[0], nErzp[1], nErzp[2]);
#endif
	      /* torno al riferimento del laboratorio */
	      for (kk1=0; kk1 < 3; kk1++)
		{
		  nErcutx[kk1] = nErcuty[kk1]=nErcutz[kk1]=0.0;
		  //Aip[kk1] = 0;
		  rErcut[kk1] = Dj[j2][kk1];
		  for (kk2=0; kk2 < 3; kk2++)
		    {
		      nErcutx[kk1] += Rl[kk2][kk1]*nErxp[kk2];
		      nErcuty[kk1] += Rl[kk2][kk1]*nEryp[kk2];
		      nErcutz[kk1] += Rl[kk2][kk1]*nErzp[kk2];
		      rErcut[kk1]  += Rl[kk2][kk1]*rErp[kk2];
		      //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
		    } 
		}
	      /* >>> scelgo un piano che sia la "media" dei piano del disco e quello perperdincolare al rim <<< */
	      for (kk1=0; kk1 < 3; kk1++)
		{
#if 0
		  Cpl[kk1] = Dj[j2][kk1];
		  npl[kk1] = nj[kk1];
#else
		  Cpl[kk1] = (Dj[j2][kk1] + Ci[kk1])*0.5;
		  if (scalProd(ni,nj) < 0.0)
		    npl[kk1] = (-ni[kk1]+nj[kk1])*0.5;
		  else
		    npl[kk1] = (ni[kk1]+nj[kk1])*0.5;
#endif
		}
	      norm = calc_norm(npl);
	      for (kk1=0; kk1 < 3; kk1++)
		npl[kk1] /= norm;
	      
	      /* passo al sistema di riferimento del piano medio (l'asse perpendicolare è x) (p) */
	      versor_to_R(npl[0], npl[1], npl[2], Rl);
	      for (kk1=0; kk1 < 3; kk1++)
		{
		  nplpx[kk1] = Rl[0][kk1];
		  nplpy[kk1] = Rl[1][kk1];
		  nplpz[kk1] = Rl[2][kk1];
		  nip[kk1] = njp[kk1]=0;
		  Cip[kk1] = Cjp[kk1]=0;
		  rErcutp[kk1]=0;
		  nErcutxp[kk1] = 0;
		  nErcutyp[kk1] = 0;
		  nErcutzp[kk1] = 0;
		  for (kk2=0; kk2 < 3; kk2++)
		    {
		      /* rim */
		      nip[kk1] += Rl[kk1][kk2]*ni[kk2];
		      Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Cpl[kk2]);
		      /* disk */
		      njp[kk1] += Rl[kk1][kk2]*nj[kk2];
		      Cjp[kk1] += Rl[kk1][kk2]*(Dj[j2][kk2]-Cpl[kk2]);
		      rErcutp[kk1] += Rl[kk1][kk2]*(rErcut[kk2]-Cpl[kk2]);
		      nErcutxp[kk1] += Rl[kk1][kk2]*nErcutx[kk2];
		      nErcutyp[kk1] += Rl[kk1][kk2]*nErcuty[kk2];
		      nErcutzp[kk1] += Rl[kk1][kk2]*nErcutz[kk2];
		    } 
		}
#if 0
		{
		  double p[3];
		  for (kk1=0; kk1 < 3; kk1++)
		    p[kk1] = rErcutp[kk1]+bErcut*nErcutzp[kk1];
		  printf("rErp=%f %f %f\n", rErcutp[0], rErcutp[1], rErcutp[2]);
		  //projectback(p, Cjp, njp);
		  printf("PRIMA DISTANCE=%.15G\n", perpcomp(p, Cip, nip));
		}
#endif

      
	      /* >> rim ellipse (ottenuta proiettando l'ellisse ottenuta in precedenza sul piano) <<< */
	      /* notare che prendendo il piano medio non potrà mai essere ni.np = 0 */
	      //lambda = -Cip[0]/nip[0];
	      /* centro dell'ellisse del rim */
      	      rErp[0] = 0;
	      rErp[1] = rErcutp[1];//Cip[1]+lambda*nip[1];
	      rErp[2] = rErcutp[2];//Cip[2]+lambda*nip[2];
	      aErcut2 = Sqr(aErcut);
	      bErcut2 = Sqr(bErcut);


#if 0
	      nErcutyp12 = Sqr(nErcutyp[1]);
	      nErcutyp22 = Sqr(nErcutyp[2]);
	      nErcutzp12 = Sqr(nErcutzp[1]);
	      nErcutzp22 = Sqr(nErcutzp[2]);
#endif
	      fact = 1.0/(-nErcutyp[2]*nErcutzp[1] + nErcutyp[1]*nErcutzp[2]);
	      ia00=nErcutzp[2]*fact;
	      ia01=-nErcutzp[1]*fact;
	      ia10=-nErcutyp[2]*fact;
	      ia11=nErcutyp[1]*fact;
	      ia002=Sqr(ia00);
	      ia102=Sqr(ia10);
	      ia112=Sqr(ia11);
	      ia012=Sqr(ia01);
	      /* e ora basta determinare autovalori e autovettori di M */
	      m00 =ia002/aErcut2 + ia102/bErcut2;
	      m01 = (ia00*ia01)/aErcut2 + (ia10*ia11)/bErcut2;
	      m10 = (ia00*ia01)/aErcut2 + (ia10*ia11)/bErcut2;
	      m11 = ia012/aErcut2 + ia112/bErcut2;
	      //printf("aErcut2=%f bErcut2=%f m00=%.15G m01=%.15G m10=%.15G, m11=%.15G\n", aErcut2, bErcut2,
		//     m00, m01, m10, m11);
	      m002 = Sqr(m00);
	      m112 = Sqr(m11);
	      nErxp[0] = 1.0;
	      nErxp[1] = 0.0;
	      nErxp[2] = 0.0;
	      delta = m002 + 4*m01*m10 - 2*m00*m11 + m112;
	      //printf("BOH ANCORA PRIMA delta=%.15G m10=%.15G\n", delta, m10);
	      //printf("m002=%.15G m01=%.15G, m10=%.15G m00=%.15G m11=%.15G, m112=%.15G\n",
		//     m002, m01, m10, m00, m11, m112);
	      if (delta < 0)
		{
		  printf("Huston abbiamo un problema...\n");
		  printf("m01=%.15G m10=%.15G m11=%.15G m01=%.15G\n", m01, m10, m11, m01);
		  printf("ia00=%.15G ia01=%.15G ia10=%.15G ia11=%.15G\n", ia00, ia01, ia10, ia11);
		  printf("delta=%.15G\n", delta);
		  printf("aErcut=%f bErcut=%f\n", aErcut, bErcut);
		  printf("BOH=%.15G\n",1./sqrt(1.0-Sqr(scalProd(nErzp,nip))));
		  printf("ni.nj=%.15G\n", scalProd(ni,nj));
		  printf("1/fact=%.15G\n",(-nErcutyp[2]*nErcutzp[1] + nErcutyp[1]*nErcutzp[2]));	
		  exit(-1);
		} 
	      invm10 = -1.0/2.0/m10;
	      AA0 = -m00+m11;
	      BB0 = sqrt(delta);
	      //printf("BOH PRIMA delta=%.15G AA0=%.15G BB0=%.15G\n", delta, AA0, BB0);
	      AA = AA0*invm10;
	      BB = BB0*invm10;
	      nEryp[0] = 0.0;
	      nEryp[1] =AA+BB;
	      nEryp[2] = 1.0;
	      norm = calc_norm(nEryp);
	      for (kk1=0; kk1 < 3; kk1++)
	       nEryp[kk1]/=norm;
	      nErzp[0] = 0.0;
	      nErzp[1] =AA-BB;
	      nErzp[2] = 1.0;
	      norm = calc_norm(nErzp);
	      for (kk1=0; kk1 < 3; kk1++)
	       nErzp[kk1]/=norm;
	      //printf("BOH DOPO delta=%.15G AA0=%.15G BB0=%.15G\n", delta, AA0, BB0);
	      AA0=m00+m11;
	      ev0 = 0.5*(AA0-BB0);
	      ev1 = 0.5*(AA0+BB0);
	      aEr = 1.0/sqrt(ev0);
	      bEr = 1.0/sqrt(ev1);   
	      //printf("fact=%.15G aErcut2=%.15G bErcut2=%.15G\n", fact, aErcut2, bErcut2);
	      //printf("ev0=%f ev1=%f semi-axes=%f %f\n", ev0, ev1, 1.0/sqrt(ev0), 1.0/sqrt(ev1));
#if 0
	      aEr=D2;	
	      bEr=D2/sqrt(1.0-Sqr(scalProd(nErzp,nip)));
	       nEryp[0]=0.0;
	      nEryp[2]=1.0/sqrt(1.0+Sqr(nrcutzp[2]/nrcutzp[1]));
	      nEryp[1]=-nEryp[2]*nrcutzp[2]/nrcutzp[1];
	      vectProdVec(nErxp,nEryp,nErzp);
#endif
	      /*
	       * c'è un modo per calcolare i semiassi e i vettori degli assi che ho trovato in rete (da verificare) 
	       */

	      //spy=scalProd(nEryp,nrcutyp);
	      //spz=scalProd(nErzp,nrcutzp);
	      //aEr = sqrt(1.0/(Sqr(spy/aErcut)+Sqr(spz/bErcut)));
	      //bEr;
#if 0
	      printf("?!?rErp=%f %f %f nErzp= %f %f %f\n", rErp[0], rErp[1], rErp[2], nErzp[0], nErzp[1], nErzp[2]);
	      printf("?!?rEdp=%f %f %f nEdzp= %f %f %f\n", rEdp[0], rEdp[1], rEdp[2], nEdzp[0], nEdzp[1], nEdzp[2]);
	      printf("ni=%f %f f%f nj=%f %f %f\n", ni[0], ni[1], ni[2], nj[0], nj[1], nj[2]);
	      printf("nip=%f %f f%f njp=%f %f %f\n", nip[0], nip[1], nip[2], njp[0], njp[1], njp[2]);
#endif
	     /* >>> disk ellipse (ottenuta proiettando il disco sul piano) <<< */
	      //lambda = -Cjp[0]/njp[0];
	      /* centro dell'ellisse del disk*/
	      rEdp[0] = 0;
	      rEdp[1] = Cjp[1];//+lambda*njp[1];
	      rEdp[2] = Cjp[2];//+lambda*njp[2];

	      nEdxp[0] = 1.0;
	      nEdxp[1]=nEdxp[2] = 0.0;
	      nEdyp[0]=0.0;
	      nEdyp[2]=1.0/sqrt(1.0+Sqr(njp[2]/njp[1]));
	      nEdyp[1]=-nEdyp[2]*njp[2]/njp[1];
	      vectProdVec(nEdxp,nEdyp,nEdzp);
	      aEd=D2;	
	      bEd=D2*scalProd(nEdxp,njp);//D2/sqrt(1.0-Sqr(scalProd(nEdzp,njp)));
	      
#if 0
	      printf("semiaxes cut ellipse: %.15G %.15G\n", aErcut, bErcut);
	      printf("ni= %.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
	      printf("semiaxes of projectd rim ellipse: %.15G %.15G\n", aEr, bEr);
	      printf("semiaxes of projectd disk ellipse: %.15G %,.15G\n", aEd, bEd);
#endif
#if 0
		{
		  double p[3];
		  for (kk1=0; kk1 < 3; kk1++)
		    p[kk1] = rErp[kk1]+bEr*nErzp[kk1];
		  printf("rErp=%f %f %f\n", rErp[0], rErp[1], rErp[2]);
		  //projectback(p, Cjp, njp);
		  printf("DOPO DISTANCE=%.15G\n", perpcomp(p, Cip, nip));
		}
#endif


	      /* passo al sistema di riferimento dell'ellisse del disk (pp) */	  
	      //printf("distance rEdp-rErp=%.15G\n",calc_distance(rErp, rEdp));
	      for (kk1=0; kk1 < 3; kk1++)
		{
		  Rl[0][kk1] = nEdxp[kk1];
		  Rl[1][kk1] = nEdyp[kk1];
		  Rl[2][kk1] = nEdzp[kk1];
		}
	      for (kk1=0; kk1 < 3; kk1++)
		{
		  nErxpp[kk1] = nErypp[kk1]=nErzpp[kk1]=0;
		  Cipp[kk1] = Cjpp[kk1]=0;
		  nipp[kk1] = njpp[kk1]=0;
		  rErpp[kk1] = 0;
		  for (kk2=0; kk2 < 3; kk2++)
		    {
		      /* rim ellipse */
		      nErxpp[kk1] += Rl[kk1][kk2]*nErxp[kk2];
		      nErypp[kk1] += Rl[kk1][kk2]*nEryp[kk2];
		      nErzpp[kk1] += Rl[kk1][kk2]*nErzp[kk2];
		      rErpp[kk1] +=  Rl[kk1][kk2]*(rErp[kk2]-rEdp[kk2]);
		      /* cylinder */
		      nipp[kk1] += Rl[kk1][kk2]*nip[kk2];
		      Cipp[kk1] += Rl[kk1][kk2]*(Cip[kk2]-rEdp[kk2]);
		      /* disk */
		      njpp[kk1] += Rl[kk1][kk2]*njp[kk2];
		      Cjpp[kk1] += Rl[kk1][kk2]*(Cjp[kk2]-rEdp[kk2]);
		    } 
		}

	      //printf("rErpp=%f %f %f nErzpp= %f %f %f\n", rErpp[0], rErpp[1], rErpp[2], nErzpp[0], nErzpp[1], nErzpp[2]);
	      //printf("Dnorm rEpp=%.15G\n",calc_norm(rErpp));
	      //printf(">>>>nxp=%f %f %f\n", nErxp[0], nErxp[1], nErxp[2]);
	      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
	      nEry1sq=Sqr(nErypp[1]);
	      nEry2sq=Sqr(nErypp[2]);
	      aErsq = Sqr(aEr);
	      bErsq = Sqr(bEr);
	      nErz1sq=Sqr(nErzpp[1]);
	      nErz2sq=Sqr(nErzpp[2]);
	      rErpp1sq=Sqr(rErpp[1]);
	      rErpp2sq=Sqr(rErpp[2]);
	      //printf("rErpp=%f %f %f nErzpp= %f %f %f\n", rErpp[0], rErpp[1], rErpp[2], nErzpp[0], nErzpp[1], nErzpp[2]);
	      coeffEr[0] = nEry1sq/aErsq + nErz1sq/bErsq;
	      coeffEr[1] = nEry2sq/aErsq + nErz2sq/bErsq;
	      coeffEr[2] = (2.0*nErypp[1]*nErypp[2])/aErsq + (2.0*nErzpp[1]*nErzpp[2])/bErsq;
	      coeffEr[3] = -1.0 + (nEry1sq*rErpp1sq)/aErsq + (nErz1sq*rErpp1sq)/bErsq + 
		(2.0*nErypp[1]*nErypp[2]*rErpp[1]*rErpp[2])/aErsq + (2.0*nErzpp[1]*nErzpp[2]*rErpp[1]*rErpp[2])/bErsq + 
		(nEry2sq*rErpp2sq)/aErsq + (nErz2sq*rErpp2sq)/bErsq;
	      coeffEr[4] = -((2.0*nEry1sq*rErpp[1])/aErsq) - (2.0*nErz1sq*rErpp[1])/bErsq - 
		(2.0*nErypp[1]*nErypp[2]*rErpp[2])/aErsq - (2.0*nErzpp[1]*nErzpp[2]*rErpp[2])/bErsq;
	      coeffEr[5] = -((2.0*nErypp[1]*nErypp[2]*rErpp[1])/aErsq) - (2.0*nErzpp[1]*nErzpp[2]*rErpp[1])/bErsq - 
		(2.0*nEry2sq*rErpp[2])/aErsq - (2.0*nErz2sq*rErpp[2])/bErsq;

	      /* ora faccio un'affinità (x'=x/aEd; y'=y/bEd) per far diventare la prima ellisse un cerchio (e trasformo conseguentemente
	       * i coefficienti coeffEr[...] dell'ellisse del rim) */ 
	      coeffEr[0] *= Sqr(aEd);
	      coeffEr[1] *= Sqr(bEd); 
	      coeffEr[2] *= aEd*bEd;
	      coeffEr[4] *= aEd;
	      coeffEr[5] *= bEd;
	      /* mi metto nel centro dell'ellisse del rim (vedi pagina wolfram su ellisse)
	       * x' = x - xEr; y' = y - yEr dove (xEr, yEr) è il centro dell'ellisse ottenuta dopo aver
	       * applicato l'affinità */
#if 0
	      delta = 4.0*coeffEr[0]*coeffEr[1] - Sqr(coeffEr[2]);
	      xEr = -((2.0*coeffEr[1]*coeffEr[4] - coeffEr[2]*coeffEr[5])/delta);
	      yEr = (coeffEr[2]*coeffEr[4] - 2.0*coeffEr[0]*coeffEr[5])/delta; 
	      /* calcolo i nuovi coeffEricienti a seguito della traslazione (notare che c0, c1 e c2 non cambiano per traslazione) */
	      coeffEr[3] = -((Sqr(coeffEr[2])*coeffEr[3] + coeffEr[1]*Sqr(coeffEr[4]) - coeffEr[2]*coeffEr[4]*coeffEr[5] 
			      + coeffEr[0]*(-4.0*coeffEr[1]*coeffEr[3] + Sqr(coeffEr[5])))/delta);
	      coeffEr[4] = 0.0;
	      coeffEr[5] = 0.0;
#endif
#if 0
		{
		  double p[3],xx, cq[3], solqa[2];
		  int numsolL;
		  /* prendo un punto sull'ellisse del rim e verfico che il punto
		   * appartiene al rim */
		  xx= D2/2.;

		  printf("xx=%f semi-axes=%f %f\n", xx, aErcut, bErcut);
		  cq[2] = coeffEr[1];
		  cq[1] = coeffEr[2]*xx;
		  cq[0] = coeffEr[3]+coeffEr[0]*Sqr(xx);
		  solve_quadratic(cq, &numsolL, solqa);
		  printf("QQQnumsol=%d\n", numsol);
		  p[1] = (xEr+xx)*aEd;
		  p[2] = (yEr+solqa[0])*bEd;
		  p[0] = 0.0;

		  printf("rErp=%f %f %f\n", rErpp[0], rErpp[1], rErpp[2]);
		  //projectback(p, Cjp, njp);
		  printf("!!!!DOPO DISTANCE=%.15G\n", perpcomp(p, Cipp, nipp));
		}
#endif


	      c0 = coeffEr[0];
	      c1 = coeffEr[1];
	      c2 = coeffEr[2];
	      c3 = coeffEr[3];
#if 1 
	      c4 = coeffEr[4];
	      c5 = coeffEr[5];
#endif
	      c02 = Sqr(c0);
	      c12 = Sqr(c1);
	      c22 = Sqr(c2);
	      c32 = Sqr(c3);
	      c42 = Sqr(c4);
	      c52 = Sqr(c5);
#if 1
	      xC=yC=0;
	      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
	      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
	      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
	      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
	      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#else
	      xC = -xEr;
	      yC = -yEr;
	      xC2=Sqr(xC);
	      yC2=Sqr(yC);
	      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
	      coeff[3] =2*c0*c2*xC + 2*c1*c2*xC - 4*c02*yC + 4*c0*c1*yC - 2*c22*yC;
	      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + 2*c02*xC2 + 2*c0*c1*xC2 + 
		c22*xC2 - 4*c0*c2*xC*yC + 6*c02*yC2 - 2*c0*c1*yC2 + c22*yC2;
	      coeff[1] = -2*c0*c2*xC + 2*c2*c3*xC + 2*c0*c2*xC2*xC + 4*c02*yC + 4*c0*c3*yC - 
		4*c02*xC2*yC + 2*c0*c2*xC*yC2 - 4*c02*yC2*yC;
	      coeff[0] = c02 + 2*c0*c3 + Sqr(c3) - 2*c02*xC2 + 2*c0*c3*xC2 + c02*xC2*xC2 
		-2*c02*yC2 - 2*c0*c3*yC2 + 2*c02*xC2*yC2 + c02*yC2*yC2;
#endif
	      solve_fourth_deg(coeff, &numsol, solqua);
	      /* ora assegno a solec[][] e calcolo x */
	      
	      //printf("numsol=%d\n", numsol);
	      for (kk1=0; kk1 < numsol; kk1++)
		{
#if 1
		  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/(c4 + c2*solqua[kk1]);
#else

		  solec[kk1][0] = (-c0 - c3 + c0*xC2 + (c0 - c1)*Sqr(solqua[kk1]) - 
				   2*c0*solqua[kk1]*yC + c0*yC2)/(2*c0*xC + c2*solqua[kk1]);
#endif
		  solec[kk1][1] = solqua[kk1];
#if 0
		  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
			 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
			 coeff[1]*solqua[kk1]+coeff[0]);
		  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
		  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
		}
	      for (kk1=0; kk1 < numsol; kk1++)
		{
		  ellips2disk(solec[kk1], solarr[kk1], 0, 0, aEd, bEd);/* torno al sistema di riferimento pp dell'ellisse del disco */
		  //printf("norm solarr[%d]=%.15G solec[]=%.15G\n", kk1, calc_norm(solarr[kk1]), solec[kk1][1]);
		  //printf("solarr=%f %f %f\n", solarr[kk1][0],solarr[kk1][1], solarr[kk1][2]);
		  projectback(solarr[kk1], Cjpp, njpp) /* ottengo i punti sul disco che corrispondonoa alle intersezioni calcolate */;
		  //printf("DISTDISTDIST*****=%.15G\n", perpcomp(solarr[kk1],Cipp,nipp));
		  //printf("DIST SOL-DISK=%.15G\n", calc_distance(solarr[kk1], Cjpp));
		}
#if 0
	      printf("cos theta=%.15G acos=%.15G\n", scalProd(nEz, nip), 180.0*acos(scalProd(nEz,nip))/M_PI);
	      printf("nEz=%.15G %.15G %.15G\n", nEz[0], nEz[1], nEz[2]);
	      printf("nEy=%.15G %.15G %.15G\n", nEy[0], nEy[1], nEy[2]);
	      printf("calcnorm nEy=%.15G\n", calc_norm(nEy));
#endif	     
	      //printf("semiassi=%f %f\n", semminE, semmaxE);
	      /* determino le coordinate del centro del cerchio rispetto al riferimento dell'ellisse */
#if 0
	      rC[0] = 0.0;
	      rC[1] = -rE[1];
	      rC[2] = -rE[2];
	      sp1 = scalProd(rC,nEy);
	      sp2 = scalProd(rC,nEz);
	      rC[1] = sp1;
	      rC[2] = sp2;
#endif
#if 0
	      printf("vers ell= nEy %f %f %f nEz %f %f %f\n", nEy[0], nEy[1], nEy[2], nEz[0], nEz[1], nEz[2]);
	      printf("prima rE=%f %f %f\n", rE[0], rE[1], rE[2]);
		{
		  double t[3];
		  for (kk=0; kk < 3; kk++)
		    {
		      t[kk] = nEz[kk]*semmaxE;
		    }
		  sp = scalProd(t, nip);
		  for (kk=0; kk < 3; kk++)
		    {
		      t[kk] -= nip[kk]*sp;
		    }
		  printf(">>>> semiax=%.15G semmin=%.15G norm=%.15G\n", semmaxE, semminE, calc_norm(t));
		}
#endif
#if 0
	      /* ora trovo l'intersezione dell'ellisse con il cerchio risolvendo l'equazione di quarto grado */
	      /* prima calcolo i coefficienti del polinomio */

	      /* coeff è un array di 5 elementi ossia a,b,c,d,e (coeff. del polinomio c0+c1*x+c2*x^2... )
	       * solarr un array con le numsol soluzioni 
	       * */
	      docirc=0;
	      a=aErp;
	      b=bErp;
	      a2=Sqr(a);
	      b2=Sqr(b);
	      a4=Sqr(a2);
	      b4=Sqr(b2);
	      R2=Sqr(D2); 
	      //xC=rC[1];
	      //yC=rC[2];
	      xC = -rErpp[1];
	      yC = -rErpp[2];

	      if (xC!=0)
		{
		  coeff[4] = 1.0/(4.0*xC2) + a4/(4.0*b4*xC2) - a2/(2.0*b2*xC2); 
		  if (coeff[4]==0)
		    docirc=1;
		}
	      else if (yC!=0)
		{
		  coeff[4] = 1.0/(4.0*yC2) - b2/(2.0*a2*yC2) + b4/(4.0*a4*yC2);
		  if (coeff[4]==0)
		    docirc=1;
		}
	      else if (semminE <= D2 && semmaxE >= D2)
		{
		  sqB = 1.0-a2/b2;
		  if (sqB==0)
		    docirc=1;
		}

	      if (docirc) 
		/* equivale a semminE=semmaxE ma assicura che non si tenti di risolvere una quartica con
		   coefficiente quartico nullo */
		{
		  /* se a=b si ha un equazione quadratica poiché si tratta di due circonferenze */
		  //double a,a2,a4,b4,R2,xC,yC,xC2,yC2;
		  //printf("doing circle\n");
#if 0
		  a=semminE;
		  a2=Sqr(a);
		  a4=Sqr(a2);
		  R2=Sqr(D2);
		  xC=rC[1];
		  yC=rC[2];
		  xC2=Sqr(xC);
		  yC2=Sqr(yC);
#endif
		  if (xC!=0)
		    {
		      coeff[2] = 1.0 + yC2/xC2;
		      coeff[1] = -yC - (a2*yC)/xC2 + (R2*yC)/xC2 - yC2*yC/xC2;
		      coeff[0] = -(a2/2.0) - R2/2.0 + a4/(4.0*xC2) - (a2*R2)/(2.0*xC2) + 
			R2*R2/(4.0*xC2) + xC2/4.0 + yC2/2.0 + (a2*yC2)/(2.0*xC2) - (R2*yC2)/(2.0*xC2) 
			+ yC*yC2/(4.0*xC2);
		      /* sto risolvendo in y */
		      solve_quadratic(coeff, &numsol, solquad);
		      /* assegno solcc e calcolo x */
		      for (kk1=0; kk1 < numsol; kk1++)
			{
			  solcc[kk1][0] = (a2 - R2 + xC2 - 2.0*solquad[kk1]*yC + yC2)/(2.0*xC);
			  solcc[kk1][1] = solquad[kk1];
			}
		      /* torno al sistema di coordinate dell'ellisse del disco e individuo i punti
		       * sul disco che metto in solarr */
		      for (kk1=0; kk1 < numsol; kk1++)
			ellips2disk(solcc[kk1],solarr[kk1], rErpp, nErypp, nErzpp, aEd, bEd);
		    }
		  else if (yC!=0)
		    {
		      coeff[2] = 1.0 + xC2/yC2;
		      coeff[1] = -xC - (a2*xC)/yC2 + (R2*xC)/yC2 - xC2*xC/yC2;
		      coeff[0] = -(a2/2.0) - R2/2.0 + xC2/2.0 + a4/(4.0*yC2) - 
			(a2*R2)/(2.0*yC2) + R2*R2/(4.0*yC2) + 
			(a2*xC2)/(2.0*yC2) - (R2*xC2)/(2.0*yC2) + xC2*xC2/(4.0*yC2) + yC2/4.0;
		      /* sto risolvendo in x */
		      solve_quadratic(coeff, &numsol, solquad);
		      /* assegno solcc e calcolo y */
		      for (kk1=0; kk1 < numsol; kk1++)
			{
			  solcc[kk1][0] = solquad[kk1];
			  solcc[kk1][1] = (a2 - R2 - 2.0*solquad[kk1]*xC + xC2 + yC2)/(2.0*yC) ;
			}
		      /* torno al sistema di coordinate del disco */
		      for (kk1=0; kk1 < numsol; kk1++)
			ellips2disk(solcc[kk1],solarr[kk1], rErpp, nErypp, nErzpp, aEd, bEd);
		    }
		  else
		    {
		      /* se semminE==D2 allora ho infinite soluzioni */

		    }
		}
	      else
		{
#if 0
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
#endif
		  /* notare che coeff[4]=0 solo se a=b che però è un caso a parte! */
		  if (xC!=0)
		    {
		      //coeff[4] =1.0/(4.0*xC2) + a4/(4.0*b4*xC2) - a2/(2.0*b2*xC2); 
#if 0
		      if (coeff[4]==0)
			{
			  printf("xCnot0 xC2=%.15G a4=%.15G b4=%.15G a2=%.15G b2=%.15G\n", xC2, a4, b4, a2, b2);
			  printf("coeff[4]=%.15G\n", coeff[4]);
			  printf("coeff[3]=%.15G\n", coeff[3]);
			  printf("ni=%.15G %.15G %.15G\n", nip[0], nip[1], nip[2]);
			  printf("a4/b4=%.15G a2/b2 %.15G xC2=%.15G\n", a4/b4, a2/b2, xC2);
			}
#endif
		      coeff[3] = -(yC/xC2) + (a2*yC)/(b2*xC2);
		      coeff[2] = 1.0/2.0 + a2/(2.0*b2) + a2/(2.0*xC2) - a4/(2.0*b2*xC2) - R2/(2.0*xC2) + 
			(a2*R2)/(2.0*b2*xC2) + (3.0*yC2)/(2.0*xC2) - (a2*yC2)/(2.0*b2*xC2);
		      coeff[1] = -yC - (a2*yC)/xC2 + (R2*yC)/xC2 - yC2*yC/xC2;
		      coeff[0] = -(a2/2.0) - R2/2.0 + a4/(4.0*xC2) - (a2*R2)/(2.0*xC2) +
			R2*R2/(4.0*xC2) + xC2/4.0+ yC2/2.0 + (a2*yC2)/(2.0*xC2) - (R2*yC2)/(2.0*xC2) + 
			yC2*yC2/(4.0*xC2) ;  
		      /* qui risolvo in y */
		      solve_fourth_deg(coeff, &numsol, solqua);
		      /* ora assegno a solec[][] e calcolo x */
		      for (kk1=0; kk1 < numsol; kk1++)
			{
			  solec[kk1][0] = (a2*b2 - b2*R2 + b2*xC2 + (-a2 + b2)*Sqr(solqua[kk1]) - 
					   2.0*b2*solqua[kk1]*yC + b2*yC2)/(2.0*b2*xC) ;
			  solec[kk1][1] = solqua[kk1];
#if 0
			  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
				 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
				 coeff[1]*solqua[kk1]+coeff[0]);
			  printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
			}
		      /* torno al sistema di coordinate del disco */
#if 0
		      printf("dopo rE=%f %f %f\n", rE[0], rE[1], rE[2]);
		      printf("perp rE=%.15G\n", perpcomp(rE, Cip, nip));
#endif
		      for (kk1=0; kk1 < numsol; kk1++)
			{
			  ellips2disk(solec[kk1],solarr[kk1], rErpp, nErypp, nErzpp, aEd, bEd);
			  if (fabs(perpcompsq(solarr[kk1], Cip, nip)-D2) > 1E-4)
			    printf("B2 perpcom=%.15G semmaxE=%.15G\n", perpcomp(solarr[kk1], Cip, nip),semmaxE);

#if 0
			    {
			      double VV[3];
			      double y, x;
			      x = sqrt(1.0-Sqr(solec[kk1][1]/semmaxE))*0.55;
			      for (kk2=0; kk2 < 3; kk2++)
				VV[kk2] = rE[kk2]+x*nEy[kk2]+solec[kk1][1]*nEz[kk2];
			      printf("ellips(solsol)=%.15G\n", Sqr(0.5/a)+Sqr(y/b)-1.0);
			      printf("BOH x=0.5 y=%.15G perpcom=%.15G\n", y, perpcomp(VV, Cip, nip)); 
			      printf("BOH2              perpcom=%.15G\n", perpcomp(solarr[kk1], Cip, nip));
			    }
#endif
			}
		    }
		  else if (yC!=0)
		    {
		      //printf("QUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
		      //coeff[4] = 1.0/(4.0*yC2) - b2/(2.0*a2*yC2) + b4/(4.0*a4*yC2);
#if 0
		      if (coeff[4]==0)
			{
			  printf("yCnot0 xC2=%.15G a4=%.15G b4=%.15G a2=%.15G b2=%.15G\n", xC2, a4, b4, a2, b2);
			}
#endif
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
			ellips2disk(solec[kk1],solarr[kk1], rErpp, nErypp, nErzpp, aEd, bEd);
		    }
		  else
		    {
		      if (a <= D2 && b >= D2)
			{
			  sqA = R2-a2;
			  //sqB = 1.0-a2/b2;
			  sqC = a2*(b2-R2);
			  sqD = b2 - a2;
			  solec[0][0] = solec[1][0] = -sqrt(sqC/sqD);
			  solec[0][1] = -(sqrt(sqA/sqB));
			  solec[1][1] = -solec[0][1];
			  solec[2][0] = solec[3][0] = -solec[0][0];
			  solec[2][1] = solec[0][1];
			  solec[3][1] = solec[1][1];
			  numsol = 4;
			  /* torno al sistema di coordinate del disco */
			  for (kk1=0; kk1 < numsol; kk1++)
			    ellips2disk(solec[kk1],solarr[kk1], rErpp, nErypp, nErzpp, aEd, bEd);
			}
		      else
			{
			  /* se non ci sono soluzioni verificare che il disco o l'ellisse non sia l'uno dentro l'altro */
			  numsol = 0;
			}
		    }
		  /* in questo caso ho due circonferenze e quindi al più due soluzioni ossia ho un'equazione quadratica */
		}
#endif
	      /* verifico che le soluzioni siano nella parte di rim che fa parte del cilindro */
	      for (kk1=0; kk1 < numsol; kk1++)
		{
#if 0
		  printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
		  printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
		  for (kk2=0; kk2 < 3; kk2++)
		    {
		      dsc[kk2] = solarr[kk1][kk2] - Cipp[kk2];
		    }
		  //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));
		  if (fabs(perpcomp(solarr[kk1], Cipp, nipp)-D2) > 1E-7)
		    {
		      printf("BOH2BOH2 perpcom=%.15G\n", perpcomp(solarr[kk1], Cipp, nipp));
		      printf("distanza punto-centro disk: %.15G\n", calc_distance(solarr[kk1], Cjpp));
		      printf("semi axes=%f %f\n", aEr, bEr);
		      printf("rcut axes=%f %f\n", aErcut, bErcut);
		      print_vec("ni=",ni);
		      print_vec("nj=",nj);
		      print_vec("npl=", npl);
		      printf("ni.nj=%.15G\n", scalProd(ni,nj));
		      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
		      for (kk2=0; kk2 < nscgbl; kk2++)
			printf("solc[%d]=%.15G\n", kk2, solcgbl[kk2]);
		      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
    		      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
    			     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
    			     coeff[1]*solqua[kk1]+coeff[0]);
    		      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
    		      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
		    }
		  sp = scalProd(dsc, nipp);
		  if (fabs(sp) < L*0.5)
		    {
		      return -1;
		    }
		}
	    }
 	  /* ========================= */

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
#endif

