/* ================= >>>  SIMULATION DEPENDENT HEADER FILE <<< ================
   This contains file contains the follow simulation dependent data types 
   and instances:
    - instance of some structures ( params, filenames and singlePar) 
    - number macros ( NUM_MISURE, BAK_STEPS, STA_STEPS ) 
    - list macros ( that is coordinates list: SAVE_LIST, ALLOC_LIST,
                    DECL_LIST ) */

/* ================== >>> PROGRAM DEFINES(CUSTOMIZE!) <<< ===================*/
/* flags per varie ottimizzazioni del codice Monte Carlo */
//#define MC_RESTR_MATRIX
//#define MD_DGEBA_NMAX
#define MC_NEW_PERC
#ifndef MC_ALMARZA
#define MC_OPT_CLSNPT
#else
#undef MC_OPT_CLSNPT
#endif
#undef MC_STORELL
#define MC_STOREBONDS
#define MC_STORE_ALL_COORDS
#define MCGC_OPTLLREBUILD
#ifdef MC_BOND_POS
/* questa ottimizzazione per il GC richiederebbe un po' di tempo per essere adattata
 * a BONDPOS per cui per ora evito di usarla */
#undef MCGC_OPTLLREBUILD
#endif 
#define MD_RAND48
#define MD_HARDSPHERES
//#undef MD_SAVE_DISTANCE
#ifdef EDHE_FLEX
#ifndef MD_ASYM_ITENS
#define MD_ASYM_ITENS
#endif
#ifndef MD_PATCHY_HE
#define MD_PATCHY_HE
#endif
/* EDHE_FLEX e MD_POLYDISP sono incompatibili */
#undef MD_POLYDISP
#endif
#ifdef MD_HE_PARALL
#include <mpi.h>
#endif
#define MDSIMUL "/home/demichel/shared/simul/mdsimul"
/* this is the executable, you must change this to your directory */
#ifdef MD_USE_CBLAS
#include <cblas.h>
#endif
#define XTERM   "/usr/X11R6/bin/nxterm"
#undef UPDATE_SYSTEM
void UpdateSystem(void);
#define UPDATE_SYSTEM UpdateSystem();
#ifdef MD_GRAVITY
#undef ADJUST_LASTCOL 
#define ADJUST_LASTCOL AdjustLastcol();
#endif
#define MD_HOME "./"
#define MD_SIMDAT MD_HOME ""
#define MD_HD_TMP MD_SIMDAT ""
/* directory to store temporary files */ 
#define Power(x,y) pow((x),(y))
#define MD_HD_MIS MD_SIMDAT "" 
/* directory to store measures files */
#if defined(MD_POLYDISP) && defined(EDHE_FLEX)
#error "-DMD_POLYDISP is not compatible with -DEDHE_FLEX!"
#endif
//#undef MD_BASIC_DT
#define MD_TAPE_TMP "/iomega/mdtmp/"
/* directory on Tape to store some temporary files (restore files and 
   measures files)*/

#define MD_TAPE_MIS "/iomega/measures/"
/* directory on tape to store measures */
#ifdef EDHE_FLEX
#define MD_COORDTMP_ASCII(x) save_coordtmp_ascii(x)
#endif
/* 16/1/1998 ADD: <XVA defines> */
#define MD_HD_XVA "/work2/xvafiles/"
#define MD_TAPE_XVA "/iomega/xvafiles/"
#ifdef EDHE_FLEX
#define MD_EDHEFLEX_OPTNNL
#define MD_OPT_MULTLL
#endif
#define MD_LL_BONDS
#define NUM_MISURE 30         /* maximum number of measure done during simulation */
#ifdef MDLLINT
#undef  MDINT
#define MDINT long long int
#define MDINTFMT "%lld"
#else
#undef  MDINT
#define MDINT int
#define MDINTFMT "%d"
#endif
#define MAXDDOT_NNL_THR 1.0E-8
#define ADJ_MAXDDOT_NL 1.001
#ifdef EDHE_FLEX
/* NOTA 16/04/2010: numero di spot per l'implementazione delle NNL tramite spots */
#define MD_SPNNL_NUMSP 8
#define MD_SPNNL_SIGMA 0.0
#define PAR2SAVE_LEN 1024
#define MD_INF_ITENS 1E199
#define MD_INF_MASS  1E199
#define MD_HANDLE_INFMASS
#ifdef MD_GHOST_IGG
typedef struct {
  int iggnum;    /* IgG a cui appartiene l'ellissoide */
  int ghost_status; /* 1 (bulk) o 2 (legto) */
} ghostInfo;
#endif
#ifdef MC_SIMUL
#define MD_STANDALONE
#endif
#ifdef MD_STANDALONE
#define  MD_CALC_VBONDING
#endif
typedef struct {
  double x[3];
  double sigma;
  int same;
} spotStruct;
/* come possono esserci n sticky spots, 
   si puo' anche pensare di mettere n corpi rigidi 
   di forma super-ellissoidale 
 */

#ifdef MC_BOND_POS
struct n2sp_struct 
{
  int i;
  int ns;
};
#endif

typedef struct {
  double x[3];
  double sax[3];
  double ppsax[3]; /* semi-lati di un parallelepipedo che circoscrive l'ellissoide con tutti i suoi sticky spots */
  double ppr[3];
  double n[3];	
#if 1
  double R[3][3]; /* orientazione relativa al sistema di riferimento del corpo rigido */
#endif
} hardobjsStruct;

typedef struct 
{
  double sax[3];
  double ppsax[3]; /* semi-lati di un parallelepipedo che circoscrive l'ellissoide con tutti i suoi sticky spots */
  double ppr[3];
  double n[3];/*super-ellipsoids double parameters (as in Povray, i.e. only n[0] and n[1] are used, n[2] is a place holder 
		for now) */
  double m;
  double I[3];
#if 0
  double xoff[3];
#endif
  int ignoreCore;/* if 1 ignore collisions of the core object */
  int brownian;
  int nspots;
#ifdef MC_BOUNDING_SPHERES
  int nspotsBS;
  double bsdiam;
#endif
  spotStruct* spots; 
  int nhardobjs; /* one can add other sub-objects to build a super-object made of 
		   several super-ellipsoids with their spots */
#ifdef MD_MULTIPLE_LL
  double rcutFact;
#endif
  hardobjsStruct* hardobjs;
} partType;

typedef struct 
{
  int min;
  int max;
} rangeStruct;

typedef struct 
{
  int type1; /* se type1 == -2 usa i range r1 */
  int spot1; 
  int type2; /* se type2 == -2 usa i range r2 */
  int spot2;
  double bheight;
  double bhin;
  double bhout;
  int nmax;
  int nr1; /* numero di range */
  rangeStruct* r1;
  int nr2; /* numero di range */
  rangeStruct* r2;
} interStruct;

typedef struct 
{
  int i; /* se i==-2 usa i range r1 */
  int spot1;
  int j; /* se j==-2 usa i range r2 */
  int spot2;
  double bheight;
  double bhin;
  double bhout;
  int nr1; /* numero di range */
  rangeStruct* r1;
  int nr2; /* numero di range */
  rangeStruct* r2;

} interStructIJ;
#endif
#ifdef MD_HE_PARALL
typedef struct 
{
  int p[2];
  double pos[6];
  double R[18];
  double vels[12];
  double axes[8];
  int cells[6];
  int lastbump[6];
  double time[4];
  double atomTime[2];
#ifdef MD_ASYM_ITENS
  double angM[2];
  double sintheta0[2];
  double costheta0[2];
  double phi0[2];
  double psi0[2];
  double RM[18];
#endif
} parall_pair_struct;
typedef struct 
{
  double t;
  double rC[3];
  int a;
  int b;
#ifdef MD_PATCHY_HE
  int sp[3];
#endif
} parall_event_struct;
#ifdef MAIN
parall_pair_struct parall_pair;
parall_event_struct parall_event;
#else
extern parall_pair_struct parall_pair;
extern parall_event_struct parall_event;
#endif
void md_mpi_init(int *argc, char***argv);
void md_mpi_finalize(void);
#undef MD_EXT_INIT
#undef MD_EXT_END
#define MD_EXT_INIT(X,Y) md_mpi_init(X,Y)
#define MD_EXT_END() md_mpi_finalize()
#endif
#define NDIM 3
#define MD_DEBUG(X)  
#define MD_DEBUG2(X)     
#define MD_DEBUG3(X) 
#define MD_DEBUG4(X) 
#define MD_DEBUG50(X) 
/* ========================================================================= */

/* ====================== >>> SIMULATION DEFINES <<< ========================*/
#ifdef ED_PARALL_DD
enum {DD_COLLISION=0, DD_CELLCROSSING};
#endif
#ifdef MD_EDHEFLEX_WALL
enum {MD_CORE_BARRIER=0,MD_INOUT_BARRIER,MD_OUTIN_BARRIER,MD_EVENT_NONE,MD_WALL};
#else
enum {MD_CORE_BARRIER=0,MD_INOUT_BARRIER,MD_OUTIN_BARRIER,MD_EVENT_NONE};
#endif
#ifndef LLONG_MAX 
#define LLONG_MAX 0x7fffffffffffffffLL
#endif
#ifndef INT_MAX
#define INT_MAX 2147483647
#endif
#ifdef MD_LL_BONDS
#define MAX_ALLOWED_INT LLONG_MAX 
#else
#define MAX_ALLOWED_INT INT_MAX
#endif
#define C_T COORD_TYPE
#define NK 10000
#ifndef EDHE_FLEX
#define NA 10 /* number of atoms for each molecule (particle) */
#else
#ifdef MD_SPHERICAL_WALL
#define MD_SPOT_GLOBAL_ALLOC
#if defined(MD_SPOT_GLOBAL_ALLOC) && defined(MD_LL_BONDS)
#define NA 1000000 /* number of atoms for each molecule (particle) */
#else
#define NA 1000
#endif
#else
#define MD_SPOT_GLOBAL_ALLOC 
/* 21/04/2010: con questa define (MD_SPOT_GLOBAL_ALLOC) non ci sono più allocazioni nello stack (locali) 
   dipendenti da NA tale soluzione è più sicura se si dovesse avere la necessità 
   di simulare tantissimi spot per particellea. Tale eventualità si puo' porre 
   per "complessi" di ellissoidi che richiedano ad es. poliedri con molte facce 
   per le NNL. */
#if defined(MD_SPOT_GLOBAL_ALLOC) && defined(MD_LL_BONDS)
#define NA 1000000
#else
#define NA 100//50000 /* number of atoms for each molecule (particle) */
#endif
#endif
#endif
#ifdef MD_LL_BONDS
#define NANA (((long long int)NA)*((long long int)NA))
#else
#define NANA (NA*NA)
#endif
#define MAXPAR 100000      /* maximum number of simulated particles */
#ifdef MC_SUS
#define MAXSUSWINDOW 1000
#endif
#ifdef MD_PATCHY_HE
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2
#define MD_PBONDS 10
#endif
#define NUM_PAR 1000   /* Number of particles for the simulation */
#define NUMK 99    /* number of k-points in which we must  calculate the 
		       structure factor */ 
#define MAXBIN 1000  /* Number of radius in which calculate the radial 
		       distribution function */
#define NUMV 150
#define GSRMAX 2.0 /* Maximum r for which we calculate the van Hove function*/
#define GSPOINT 30  /* Pont for which we plot the van Hove function */

#define FSPOINT 30    /* Point for which we plot the self-part of the 
			 intermediate scattering function */
#define FSKMAX 15.0 /* see the maximum k for the static structure factor */

#define PE_POINTS 300

/* ========================================================================= */

/* ======= >>> YOU MUST CHANGE THE SAVE_LIST AND ALLOC_LIST <<< ========== */

/* Follows a list of pointer to COORD_TYPE, you must think these objects like
   array of COORD_TYPE, in this implementation you can allocate whatever
   COORD_TYPE array you want, each sizeof(COORD_TYPE)*<particles number> bytes
   long 
   WARNING: you can't allocate more than 120 COORD_TYPE array.
   TODO: - Implement a better way to assign shared memory pointers, avoiding
           to use one "shared" pointer for each COORD_TYPE array ( for example
	   it could allocate shared memory in blocks of 16Mb, 
	   filling each one ).
	 - Implement doubly dimensioned array as a definition apart.
*/
#ifdef MD_POLYDISP
#define MD_SAVE_EXTRAS , axaP, axbP, axcP
#define MD_ALLOC_EXTRAS , &axaP, &axbP, &axcP
#define MD_DECL_EXTRAS  , *axaP, *axbP,  *axcP
#else
#define MD_SAVE_EXTRAS
#define MD_ALLOC_EXTRAS
#define MD_DECL_EXTRAS
#endif
#ifdef MD_ASYM_ITENS
#ifdef MD_GRAVITY
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, Mx, My, Mz, lastcol MD_SAVE_EXTRAS
#else
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, Mx, My, Mz MD_SAVE_EXTRAS
#endif
#else
#ifdef MD_GRAVITY
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, vx, vy, vz, wx, wy, wz, lastcol MD_SAVE_EXTRAS
#else
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, vx, vy, vz, wx, wy, wz MD_SAVE_EXTRAS
#endif
#endif
#undef  EXT_SLST
#ifdef MD_LXYZ
#define EXT_SLST  &L[0], &L[1], &L[2]
#else
#ifdef MD_GRAVITY
#define EXT_SLST  &L, &Lz
#else
#define EXT_SLST  &L
#endif
#endif
/* Reduced list of variables to save on xva file (tape file) */ 
#ifdef MD_GRAVITY
#define XVA_LIST rx, ry, rz, vx, vy, vz, lastcol
#define XVA_NUM 7 /* this is the number of vars in XVA_LIST <--- SET THIS!!!*/
#else
#define XVA_LIST rx, ry, rz, vx, vy, vz
#define XVA_NUM 6 /* this is the number of vars in XVA_LIST <--- SET THIS!!!*/
#endif

#ifdef MD_GRAVITY
#define XVA_ALST &rx, &ry, &rz, &vx, &vy, &vz, &lastcol
#define XVA_DLST *rx, *ry, *rz, *vx, *vy, *vz, *lastcol
#else
#define XVA_ALST &rx, &ry, &rz, &vx, &vy, &vz
#define XVA_DLST *rx, *ry, *rz, *vx, *vy, *vz
#endif
/* Follows a list of pointer to pointer to COORDTYPE.
   For ex., let ax1 be of type COORD_TYPE*, that is a pointer to COORD_TYPE, 
   then if we take its address by the '&' operator we obtain a pointer 
   to a pointer.
   You must understand it considering that taking an address of a pointer,
   we obtain a pointer (address) to a pointer.
   This list is used by AllocCoord() to allocate shared memory for these 
   arrays (to see how AllocCoord() works see AllocCoord() code in 
   mdarray.c file).
*/
#ifdef MD_ASYM_ITENS
#ifdef MD_GRAVITY
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &Mx, &My, &Mz, &lastcol MD_ALLOC_EXTRAS
#else
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &Mx, &My, &Mz MD_ALLOC_EXTRAS
#endif
#else
#ifdef MD_GRAVITY
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &lastcol MD_ALLOC_EXTRAS
#else
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz MD_ALLOC_EXTRAS
#endif
#endif
/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#ifdef MD_ASYM_ITENS
#ifdef MD_GRAVITY
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *Mx, *My, *Mz, *lastcol MD_DECL_EXTRAS
#else
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *Mx, *My, *Mz MD_DECL_EXTRAS
#endif
#else
#ifdef MD_GRAVITY
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *lastcol MD_DECL_EXTRAS
#else
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz MD_DECL_EXTRAS
#endif
#endif				   
#undef EXT_DLST
#ifdef MD_LXYZ
#define EXT_DLST  L[3] 
#else
#ifdef MD_GRAVITY
#define EXT_DLST  L, Lz 
#else
#define EXT_DLST  L 
#endif
#endif
#ifdef ED_PARALL_DD
#define treeLeftBZ   treeBZ[0]
#define treeRightBZ  treeBZ[1]
#define treeUpBZ     treeBZ[2]
#define treeIdABZ    treeBZ[3]
#endif
#if 1
#define treeLeft   tree[0]
#define treeRight  tree[1]
#define treeUp     tree[2]
#define treeCircAL tree[3]
#define treeCircBL tree[4]
#define treeCircAR tree[5]
#define treeCircBR tree[6]
#define treeIdA    tree[7]
#define treeIdB    tree[8]
#define treeStatus tree[9]
#ifdef MD_PATCHY_HE
#define treeIdC    tree[10]
#define treeIdD    tree[11]
#define treeIdE    tree[12]
#ifdef MD_CALENDAR_HYBRID
#define treeNext   tree[13]
#define treePrev   tree[14]
#define treeQIndex tree[15]
#endif
#else
#ifdef MD_CALENDAR_HYBRID
#define treeNext   tree[10]
#define treePrev   tree[11]
#define treeQIndex tree[12]
#endif
#endif
#endif
#if 0
#ifdef MD_CALENDAR_HYBRID
#define nlists 1000000
#define scale 200000
#endif
#endif
#ifdef MD_PATCHY_HE
struct LastBumpS 
{
  int mol;
  int at;
  int type;
};
#endif
#define ATOM_LIMIT 10000000
#ifdef POVRAY
#define POVCLEN 128
#endif
struct nebrTabStruct 
{
  int *list;       /* ellissoidi nella NNL */
  double **shift;
  int len;         /* numero di ellissoidi nella NNL */
  double time;     /* tempo a cui è stata costruita la NNL */
  double nexttime; /* tempo a cui la NNL deve essere ricostruita */
  double **R;
  double r[3];
  double axa;
  double axb;
  double axc;
};
/* ======================== >>> struct progStatus <<< =======================*/
struct progStatus
{
  /* DESCRIPTION:
     This structure, together with params structure is the header of restore
     files (see TECH_INFO).
     It contains datas used by the program and it is saved 
     into the restore files ( see TECH_INFO file) */
  char xvafile[NAME_LENGTH]; /* file containing positions, velocities and 
				accelerations, saved at regular intervals
			       (this is the "tape file" of Allen-Tildesley) */
  char inifile[NAME_LENGTH];
  int iniFormat; /* 0 = binary 1 = ascii 2 = both */
  int endFormat; /* 0 = binary 1 = ascii 2 = both */
#ifdef EDHE_FLEX
#ifdef POVRAY
  double povtransmit;
  char povpcol[POVCLEN];
  char povcolA[POVCLEN];
  char povcolB[POVCLEN];
  double povcx;// camera position
  double povcy;
  double povcz;
  double povlx;// light source position
  double povly;
  double povlz;
#endif
  int optbm;
#ifdef MD_SCALEPHI_STAGES
  int growthType;
#endif
  int frozenDOF;
#endif
#ifdef MD_EDHEFLEX_WALL
  int hardwall;
#endif
#ifdef MD_BILOG
  double basew;
  int lastbilogsaved;
#endif
#ifdef MC_ELCONST_MC
  double lp;
  int calcvexcl;
  int eqstps;
  int polylen;
  double alpha;
  int curi[2];
  double totene[3];
  double tottrials;
#endif
  char endfile[NAME_LENGTH];
  /* stringa contenente il nome del file in cui
     si salvano le coordinate delle particelle alla fine della simulazione.
     Se il crash del sistema avviene proprio durante la copia in tale file
     allora le coordinate delle particelle saranno nei file temporanei
     (ved. prima) */
  char dataFiles[NUM_MISURE][NAME_LENGTH];
  /* array di stringhe con i nomi dei files su cui memorizzare i valori 
     calcolati durante la simulazione (ad es.coeff. di diff.).
     Tale array va allocato nel file incluso alla fine che dipende
     dalla simulazione */
  int tapeTimes;    /* every 'tapeTimes * Omeasure[].saveSteps' 
				steps save on Tape measures ans 
				restore files, if '0' don't use Tape */ 
  MDINT xvaSteps;     /* steps between two tape file savings */
  MDINT bakSteps;    /* steps between two savings of restore files on HD*/

  MDINT bakStepsAscii; 
  MDINT staSteps;     /* steps after which must save sim_stat structure 
				into the STATUS_FILE */
  MDINT measSteps[NUM_MISURE];/*steps after which save every measure */
  
  MDINT measCalc[NUM_MISURE]; /*steps between two measure calculation */
				        
  MDINT initCalc[NUM_MISURE];
  /* first step to begin to calculate a certain measure, by default 
     this quantity is 1, but it is possible to change this value in the 
     parameter file */

  MDINT initStep[NUM_MISURE]; 
  /* first step to begin to save a certain measure, by default this quantity 
     is 1, but it is possible to change this value in the parameter file*/

  /* =============== >>> PUT HERE YOUR STATUS FILEDS <<< =================== 
     For example accumalators (see Allen - Tildesley)*/

  COORD_TYPE sumEta; /* accumulators for obtaining the mean value of eta */
  
  /* Accumulators for the integral of angular velocity */
#ifdef MD_HSVISCO
  double DQxx;
  double DQyy;
  double DQzz;
  double DQTxy;
  double DQTyz;
  double DQTzx;
  double DQTxx;
  double DQTyy;
  double DQTzz;
  double DQWxy;
  double DQWyz;
  double DQWzx;
  double DQWxx;
  double DQWyy;
  double DQWzz;
  double DQWxxHS;
  double DQWyyHS;
  double DQWzzHS;
  double DQWxxST;
  double DQWyyST;
  double DQWzzST;
  double Txy;
  double Tyz;
  double Tzx;
  double Txx;
  double Tyy;
  double Tzz;
#endif
  COORD_TYPE DQxy;
  COORD_TYPE DQyz;
  COORD_TYPE DQzx;

  COORD_TYPE PxyArr[5];
  COORD_TYPE PyzArr[5];
  COORD_TYPE PzxArr[5];
#ifdef MD_DYNAMIC_OPROG
#ifdef MD_CALC_DPP
  double *sumdx;
  double *sumdy;
  double *sumdz;
  double *lastu1x;
  double *lastu1y;
  double *lastu1z;
  double *lastu2x;
  double *lastu2y;
  double *lastu2z;
  double *lastu3x;
  double *lastu3y;
  double *lastu3z;
#endif
  void *ptr;
  int len;
  int (*dyn_alloc_oprog)(void);
  void (*set_dyn_ascii)(void);
  double *sumox;
  double *sumoy;
  double *sumoz;
  double *lastcolltime;
#if 0
  double *vcmx0;
  double *vcmy0;
  double *vcmz0;
#endif  
  double *rxCMi; /* initial coordinates of center of mass */
  double *ryCMi; /* MAXPAR is the maximum number of particles */
  double *rzCMi;
  double **DR;
#ifdef MC_SUS
  double *sushisto;
#endif
#endif
#ifdef MC_BOUNDING_SPHERES
  int useboundsph;
#endif
#ifndef MD_DYNAMIC_OPROG
  double sumox[MAXPAR];
  double sumoy[MAXPAR];
  double sumoz[MAXPAR];
#ifdef MD_CALC_DPP
  double sumdx[MAXPAR];
  double sumdy[MAXPAR];
  double sumdz[MAXPAR];
  double lastu1x[MAXPAR];
  double lastu1y[MAXPAR];
  double lastu1z[MAXPAR];
  double lastu2x[MAXPAR];
  double lastu2y[MAXPAR];
  double lastu2z[MAXPAR];
  double lastu3x[MAXPAR];
  double lastu3y[MAXPAR];
  double lastu3z[MAXPAR];
#endif
  double lastcolltime[MAXPAR];
#endif
#ifdef MD_CALENDAR_HYBRID
  int nlistsHQ;
  double scaleHQ;
  int adjustHQ;
  double baseIndex;
  int curIndex;
  int overthrHQ;
#endif
  double springkSD;
  int SDmethod;
  double toldxNR;
  /* toldx da usare come secondo tentativo se fallisce la prima volta il calcolo della distanza 
     con SDmethod=2 */
  double toldxNRta;
  double tolAngNR;
  double stepSDA;
  double stepSDB;
  int maxitsSD;
  double tolSD;
  double tolSDlong;
  double tolSDconstr;
  double tolSDgrad;
  double tolAngSD;
  double epsdSD;
  double epsdGDO;
  /* Accumulator for the radial distribution function */
  int hist[MAXBIN];
  int equilibrated;
  double eqlevel;
  int zbrakn;
  int n1;
  int n2;
  double zbrentTol;
  /* Accumulator for the static structure factor */ 
  COORD_TYPE sumS[NUMK];

  /* Accumulator for the MAXWELL-BOLTZMANN distribution */ 
  int histMB[NUMV];

  COORD_TYPE sumTemp;
  COORD_TYPE sumPress;
  
#ifndef MD_DYNAMIC_OPROG
#if 0
  COORD_TYPE vcmx0[MAXPAR];
  COORD_TYPE vcmy0[MAXPAR];
  COORD_TYPE vcmz0[MAXPAR];
#endif
  COORD_TYPE rxCMi[MAXPAR]; /* initial coordinates of center of mass */
  COORD_TYPE ryCMi[MAXPAR]; /* MAXPAR is the maximum number of particles */
  COORD_TYPE rzCMi[MAXPAR];
  COORD_TYPE DR[MAXPAR][3];
#ifdef MC_SUS
  double sushisto[MAXSUSWINDOW];
#endif
#endif
  COORD_TYPE W;

  int savedXva; 
  int CMreset;
  int mdseed;
  double lastcoll;
  double rNebrShell;   /* Dr of shell of neighbour list shell see Rapaport pag. 53 */
  int nebrTabFac;
  int useNNL;
  int dist5;
  int dist8stps;
  int dist5NL;
  int paralNNL;
  COORD_TYPE tolT;
  int ipart;
  int HNBOX;
  int avngS;
  int avnggr;
  int avngTemp;
  int avngPress;
  int avngMB;
  char nRun[32];
  double rescaleTime;
  int scalevel;
  double endtime;
  double scaleVolTime;
  int brownian;
  double h;
#ifdef MD_EDHEFLEX_WALL
  double epsdPlane;
  double epsdFastPlane;
#endif
  double epsd;
  double epsdFast;
  double epsdFastR;
  double epsdMax;
  double epsdNL;
  double epsdFastNL;
  double epsdFastRNL;
  double epsdMaxNL;
#ifdef MD_PATCHY_HE
  int autocat;
  /* parametri per la reazione autocatalitica: k_c(p) = k0 + k1 * p^m */
  double k0;
  double k1;
  double mac;
  double epsdSP;
  double epsdFastSP;
  double epsdSPNL;
  double epsdFastSPNL;
#endif
#ifdef MD_ABSORPTION
  double bufHeight;
#ifdef MD_SPHERICAL_WALL
  int halfsolidangle;
#endif
#endif
  int dofA;
  int dofB;
  int guessDistOpt;
  int forceguess;
  double targetPhi;
#ifdef MC_RESTR_MATRIX
  char restrMatrix[NAME_LENGTH];
#endif
#ifdef MC_KERN_FRENKEL
  int polylen;
  double costhKF;
  double distKF;
#endif
#ifdef MC_GAPDNA
  int polylen;
#ifdef GAPDNA_BENDING_ENERGY
  double kbend;
#endif
#endif
#ifdef MC_AMYLOID_FIBRILS
  double tors_k;
  double tors_theta0;
#endif
#ifdef MC_SIMUL
  double targetPhiMC;
#endif
#ifdef MD_POLYDISP
#ifdef MD_POLYDISP_XYZ
  double polydispX;
  double polydispY;
  double polydispZ;
#else
  double polydisp;
#endif
  double polycutoff;
#endif
  double scalfact;
  double reducefact;
  double phitol;
  double axestol;
  double minDist;
  double rmsd2end;
  double tmsd2end;
#ifdef MD_GRAVITY
  double taptau;
  int tapampl;
  double checkquenchTime;
#endif
  double nextcheckTime;
  double nextSumTime;
  double nextDt;
#ifdef MD_BIG_DT
  double refTime;
  double bigDt;
#endif
  /* questi servono per salvare le conf usando la stessa formula di Giuseppe */
  double nextStoreTime;
  int KK;
  int JJ;
  double storerate;
#ifdef MD_PATCHY_HE
  int checkGrazing;
  int maxbonds;
  int assumeOneBond;
#endif
#ifdef MD_GRAVITY
  int numquench;
  int maxquench;
  double quenchtol;
  double rhobh;
  double extraLz;
  double vztap;
  double rzup;
  double expandFact;
  double quenchend;
  double tc;
  double lastV;
  double accV;
  int wallcollCount;
  double accrcmz;
#endif
  double intervalSum;
  int eventMult;
#if defined(MD_INELASTIC) || (MD_GRAVITY)
  double tc;
#endif
  /* ADD 13/4/2000 
     Logarithmic saving implemented for xva file */
  int xvaSaveMode;/* 0 = linear 1 = semilog 2 = bilog (not impl. yet) */
  double overlaptol;
  /* AGGIUNTI 24/04/01 */
  int bakSaveMode;/* save mode for ascii backup */
  char tmpPath[NAME_LENGTH];
  char misPath[NAME_LENGTH];
  /* =================================================== */
  /* questi servono per il salvataggio bilog */
  double base;    /* We save at base^^NN step */
  int NN;         /* Logatithmic block length */
  double fstps;         /* There are KK block each base^NN long */
  int eventCount;
  int collCount;
  int crossCount;
  int PE[PE_POINTS];
  double ENmin;
  double ENmax;
#ifdef EDHE_FLEX
  int n_gauleg; /* ordine */
  char par2save[PAR2SAVE_LEN];
  int stripStore;
  double Tf; /* temperatura finale per il quench */
  double xi; /* rate di riduzione durante il quench per il folding */
#endif
#ifdef ALIGN_POT
  double Ual;
  int alignaxis; /*0=x 1=y 2=z*/
#endif
#ifdef MD_SUBENZYME
  double rateSE[10];
  int SEreaction;
  double SEp;
#endif
#ifdef MD_EDHEFLEX_OPTNNL
  int optnnl;
#endif
#if defined(MD_RABBIT) || defined(MD_NANOBODY)
  double first_time;
  double time_limit;
  double rate[10];
  double rhozBinSize;
#endif
#ifdef MD_PROTEIN_DESIGN
  char nativeConf[NAME_LENGTH];
#endif
#ifdef MD_MULTIPLE_LL
  int multipleLL;
  double rcutfactMLL;
#endif
  /* ======================================================================= */
#ifdef MC_SIMUL
#ifdef MC_FREEZE_BONDS
  int freezebonds;
#endif
#ifdef MC_NVE
  double Ed;
#endif
#ifdef MC_CLUSTER_MOVE
  double clsmovprob;
  double delTclsMC;
  double delRclsMC;
#endif
#ifdef MC_BIGROT_MOVE
  double bigrotmov;
#ifdef MC_BIGROT_BIASED
  double bigrotbias;
  double bigrotTheta0;
#endif
#endif
  int restrmove; /* 0=no restriction (default) 1=fixed rot */
#ifdef MC_GRANDCAN
  double zetaMC;
  int npav;
  int nexc;
#ifdef MC_SUS
  int susnmax;
  int susnmin;
#endif
#if defined(MC_SWHC) || defined(MC_SWELL)
  double deltasw[2];
  int constDelta;
#endif
#ifdef MC_HYDROPHOBIC_INT
  int maxtrialsH;
#endif
#ifdef MC_HELIX
  double pitch;
  int Nxi;
  double sighelix;
  double lenhelix;
  double radhelix;
#endif
#endif
#ifdef MC_FLIP_MOVE
  double flip_prob;
#endif
  int nvbbias;
  double pbias;
  double lastNNLrebuildMC;
  double vbond;
  double targetAccept;
  double targetAcceptVol;
  double dthetaMC;
  double deltaMC;
  int ensembleMC;
#ifdef MC_ALMARZA
  double almarza_thr;
#endif
  double vmax;
  int outMC;
  int adjstepsMC;
  int resetacceptVol;
  int resetaccept;
#endif
#ifdef MD_SURV_PROB
  double spdeltat;
#endif
};

/* ======================== >>> filenames instance <<< ======================= 
   MAIN is a macro defined only in mdsimul.c, so there we put the instance
   of the object, and elsewhere we put only extern declarations */
#ifdef MAIN
struct progStatus OprogStatus;
#else
extern struct progStatus OprogStatus;
#endif

/* ========================== >>> struct params <<< =========================*/
struct params
{
  /* DESCRIPTION:
     This is the header put at the top of a coordinate file. 
     It is strongly dependent upon simulation, besides to make 
     it more general as possible, it has some unused attribs */

  /* =================== >>> DON'T TOUCH THESE !!!! <<< =================== */
  int parnum;        	/* total number of particles */
  int parnumA;          /* number of particles A */
  MDINT totStep;	/* temporal step number that simulation 
				   must do */
  MDINT curStep;	/* current step of simulation */
  double time;
  double Dt;
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  double wallDiss;             /* dissipazione negli urti contro il muro */
  double partDiss;             /*dissipazione negli urti fra particelle */
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m[2];             /* atoms masses */
#ifdef EDHE_FLEX
  int ntypes;
  int ninters;
  int saveBonds;
  int maxbondsSaved;
  int nintersIJ;
#endif
#ifdef MD_GHOST_IGG
  int ghostsim;
#endif
double a[2];
  double b[2];
  double c[2];
  double rcut;
#ifdef MC_BOND_POS
  double rcutBP;
#endif
#ifdef MC_BOUNDING_SPHERES
  double rcutBS;
#endif
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
#ifndef MD_ASYM_ITENS
  double I[2];
#else
  double I[2][3];
#endif
#ifdef MD_PATCHY_HE
  int nmax;
  double Dr;
  double theta;  
  double sigmaSticky; /* ampiezza della buca */
  double bheight;
  double bhin;
  double bhout;
#endif
#ifdef MD_GRAVITY
  double ggrav;
#endif
  COORD_TYPE tol;               /* Tolerance of the shake algoritm used 
				   by RATTLE */
  /*======================================================================= */
};

/* ======================== >>> params instance <<< ========================
   MAIN is a macro defined only in mdsimul.c, so there we put the instance
   of the object, and elsewhere we put only extern declarations */
#ifdef MAIN
struct params Oparams;
#else
extern struct params Oparams;
#endif

#define OS(_A_) OprogStatus._A_
#define OP(_A_) Oparams._A_

/* ======================== >>> opro_ascii <<< =========================== */
#ifdef MAIN
struct pascii opro_ascii[] =
{
  /*
  {ID, campo di OprogStatus, Qty, formato per printf        } */
  {"xvafile",      OS(xvafile),                    1,     NAME_LENGTH, "%s"},
  {"inifile",      OS(inifile),                    1,     NAME_LENGTH, "%s"},
  {"endFile",      OS(endfile),                     1,   NAME_LENGTH,   "%s"},
  {"iniFormat",    &OS(iniFormat),                  1,               1, "%d"},
  {"endFormat",    &OS(endFormat),                  1,               1, "%d"},
#ifdef EDHE_FLEX
#ifdef POVRAY
  {"povpcol",      OS(povpcol),                     1,   POVCLEN,   "%s"},
  {"povcolA",      OS(povcolA),                     1,   POVCLEN,   "%s"},
  {"povcolB",      OS(povcolB),                     1,   POVCLEN,   "%s"},
  {"povtransmit",  &OS(povtransmit),                1,              1, "%.10G"},
  {"povcx",         &OS(povcx),                1,              1, "%.10G"},
  {"povcy",         &OS(povcy),                1,              1, "%.10G"},
  {"povcz",         &OS(povcz),                1,              1, "%.10G"},
  {"povlx",         &OS(povlx),                1,              1, "%.10G"},
  {"povly",         &OS(povly),                1,              1, "%.10G"},
  {"povlz",        &OS(povlz),                1,              1, "%.10G"},
#endif
  {"optbm",        &OS(optbm),                      1,               1, "%d"},
#ifdef MD_SCALEPHI_STAGES
  {"growthType",      &OS(growthType),                1,       1,  "%d" },
#endif
  {"frozenDOF",      &OS(frozenDOF),                1,       1,  "%d" },
#endif
#ifdef MD_EDHEFLEX_WALL
  {"hardwall",     &OS(hardwall),                   1,               1, "%d"},
#endif
  {"dataFiles",    OS(dataFiles),         NUM_MISURE,     NAME_LENGTH, "%s"},
  {"xvaSteps",     &OS(xvaSteps),                  1,               1, MDINTFMT},
  {"bakSteps",     &OS(bakSteps),                  1,               1, MDINTFMT},
  {"bakStepsAscii",&OS(bakStepsAscii),             1,               1, MDINTFMT},
  {"staSteps",     &OS(staSteps),                  1,               1, MDINTFMT},
  {"measSteps",    OS(measSteps),         NUM_MISURE,               1, MDINTFMT},
  {"measCalc",     OS(measCalc),          NUM_MISURE,               1, MDINTFMT},
  {"initCalc",     OS(initCalc),          NUM_MISURE,               1, MDINTFMT},
  {"initStep",     OS(initStep),          NUM_MISURE,               1, MDINTFMT},
  {"sumEta",       &OS(sumEta),                     1,              1, "%.10G"},
  {"DQxy",         &OS(DQxy),                       1,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                       1,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                       1,              1, "%.10G"},
#ifdef MD_HSVISCO
  {"DQxx",         &OS(DQxx),                       1,              1, "%.10G"},
  {"DQyy",         &OS(DQyy),                       1,              1, "%.10G"},
  {"DQzz",         &OS(DQzz),                       1,              1, "%.10G"},
  {"DQTxy",        &OS(DQTxy),                       1,              1, "%.10G"},
  {"DQTyz",        &OS(DQTyz),                       1,              1, "%.10G"},
  {"DQTzx",        &OS(DQTzx),                       1,              1, "%.10G"},
  {"DQTxx",        &OS(DQTxx),                       1,              1, "%.10G"},
  {"DQTyy",        &OS(DQTyy),                       1,              1, "%.10G"},
  {"DQTzz",        &OS(DQTzz),                       1,              1, "%.10G"},
  {"DQWxy",        &OS(DQWxy),                       1,              1, "%.10G"},
  {"DQWyz",        &OS(DQWyz),                       1,              1, "%.10G"},
  {"DQWzx",        &OS(DQWzx),                       1,              1, "%.10G"},
  {"DQWxx",        &OS(DQWxx),                       1,              1, "%.10G"},
  {"DQWyy",        &OS(DQWyy),                       1,              1, "%.10G"},
  {"DQWzz",        &OS(DQWzz),                       1,              1, "%.10G"},
  {"DQWxxHS",      &OS(DQWxxHS),                       1,              1, "%.10G"},
  {"DQWyyHS",      &OS(DQWyyHS),                       1,              1, "%.10G"},
  {"DQWzzHS",      &OS(DQWzzHS),                       1,              1, "%.10G"},
  {"DQWxxST",      &OS(DQWxxST),                       1,              1, "%.10G"},
  {"DQWyyST",      &OS(DQWyyST),                       1,              1, "%.10G"},
  {"DQWzzST",      &OS(DQWzzST),                       1,              1, "%.10G"},
  {"Txy",          &OS(Txy),                       1,              1, "%.10G"},
  {"Tyz",          &OS(Tyz),                       1,              1, "%.10G"},
  {"Tzx",          &OS(Tzx),                       1,              1, "%.10G"},
  {"Txx",          &OS(Txx),                       1,              1, "%.10G"},
  {"Tyy",          &OS(Tyy),                       1,              1, "%.10G"},
  {"Tzz",          &OS(Tzz),                       1,              1, "%.10G"},
  {"lastcoll",     &OS(lastcoll),                   1,              1, "%.15G"},
#endif
  {"PxyArr",       OS(PxyArr),                      5,              1, "%.10G"},
  {"PyzArr",       OS(PyzArr),                      5,              1, "%.10G"},
  {"PzxArr",       OS(PzxArr),                      5,              1, "%.10G"},
#ifdef MD_DYNAMIC_OPROG
  {"sumox",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"sumoy",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"sumoz",        NULL,                       -MAXPAR,        1, "%.15G"},
#else
  {"sumox",        OS(sumox),                       -MAXPAR,        1, "%.15G"},
  {"sumoy",        OS(sumoy),                       -MAXPAR,        1, "%.15G"},
  {"sumoz",        OS(sumoz),                       -MAXPAR,        1, "%.15G"},
#endif
#ifdef MD_CALENDAR_HYBRID
  {"nlistsHQ",     &OS(nlistsHQ),                    1,   1,  "%d"},
  {"scaleHQ",      &OS(scaleHQ),                    1,   1,  "%.15G"},
  {"adjustHQ",     &OS(adjustHQ),                   1,   1,  "%d"},
  {"overthrHQ",    &OS(overthrHQ),                  1,   1,  "%d"},
#if 0
  /* questi non serve salvarli poichè li inizializza ogni volta in InitEventList() */
  {"baseIndex",    &OS(baseIndex),                  1,   1,  "%.15G"},
  {"curIndex", &OS(curIndex),               1,   1,  "%.15G"},
#endif
#endif
#ifndef MD_DYNAMIC_OPROG
#ifdef MD_CALC_DPP
  {"sumdx",        OS(sumdx),                       -MAXPAR,        1, "%.15G"},
  {"sumdy",        OS(sumdy),                       -MAXPAR,        1, "%.15G"},
  {"sumdz",        OS(sumdz),                       -MAXPAR,        1, "%.15G"},
  {"lastu1x",       OS(lastu1x),                      -MAXPAR,        1, "%.15G"},
  {"lastu1y",       OS(lastu1y),                      -MAXPAR,        1, "%.15G"},
  {"lastu1z",       OS(lastu1z),                      -MAXPAR,        1, "%.15G"},
  {"lastu2x",       OS(lastu2x),                      -MAXPAR,        1, "%.15G"},
  {"lastu2y",       OS(lastu2y),                      -MAXPAR,        1, "%.15G"},
  {"lastu2z",       OS(lastu2z),                      -MAXPAR,        1, "%.15G"},
  {"lastu3x",       OS(lastu3x),                      -MAXPAR,        1, "%.15G"},
  {"lastu3y",       OS(lastu3y),                      -MAXPAR,        1, "%.15G"},
  {"lastu3z",       OS(lastu3z),                      -MAXPAR,        1, "%.15G"},
#endif
#endif
#ifdef MD_DYNAMIC_OPROG
  {"rxCMi",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"ryCMi",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"rzCMi",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"DR",           NULL,                       -MAXPAR,        3, "%.15G"}, 
  {"lastcolltime", NULL,                       -MAXPAR,        1, "%.15G"},    
#ifdef MD_CALC_DPP
  {"sumdx",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"sumdy",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"sumdz",        NULL,                       -MAXPAR,        1, "%.15G"},
  {"lastu1x",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu1y",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu1z",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu2x",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu2y",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu2z",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu3x",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu3y",      NULL,                      -MAXPAR,        1, "%.15G"},
  {"lastu3z",      NULL,                      -MAXPAR,        1, "%.15G"},
#endif
#ifdef MC_SUS
  {"sushisto",     NULL,                      -MAXSUSWINDOW,        1, "%.15G"},
#endif
#else
  {"rxCMi",        OS(rxCMi),                       -MAXPAR,        1, "%.15G"},
  {"ryCMi",        OS(ryCMi),                       -MAXPAR,        1, "%.15G"},
  {"rzCMi",        OS(rzCMi),                       -MAXPAR,        1, "%.15G"},
  {"DR",           OS(DR),                          -MAXPAR,        3, "%.15G"}, 
  {"lastcolltime", OS(lastcolltime),                -MAXPAR,        1, "%.15G"},
#ifdef MC_SUS
  {"sushisto",     OS(sushisto),                      -MAXSUSWINDOW,       1, "%.15G"},
#endif
#endif
  {"hist",         OS(hist),                  MAXBIN,               1, "%d"},
  {"sumS",         OS(sumS),                    NUMK,               1, "%.6G"},
  {"histMB",       OS(histMB),                  NUMV,               1, "%d"},
  {"sumTemp",      &OS(sumTemp),                    1,              1, "%.6G"},
  {"sumPress",     &OS(sumPress),                   1,              1, "%.6G"},
#ifdef MC_ELCONST_MC
  {"alpha",    &OS(alpha),   1, 1, "%.10G"},
  {"lp",       &OS(lp),      1, 1, "%.10G"},
  {"curi",    OS(curi),    2, 1, "%d"}, 
  {"totene",  &OS(totene), 3, 1, "%.15G"},
  {"tottrials", &OS(tottrials), 1, 1, "%.15G"},
  {"polylen", &OS(polylen), 1, 1, "%d"},
  {"eqstps",  &OS(eqstps),  1, 1, "%d"},
  {"calcvexcl", &OS(calcvexcl), 1, 1, "%d"},
#endif
  //  {"sumMedia",     &OS(sumMedia),                   1,   1, "%.6f"},
  //{"sumVx",        OS(sumVx),                    MAXPAR,  1, "%.10f"},
  //{"sumVy",        OS(sumVy),                    MAXPAR,  1, "%.10f"},
  //{"sumVz",        OS(sumVz),                    MAXPAR,  1, "%.10f"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
#ifdef MD_POLYDISP  
#ifdef MD_POLYDISP_XYZ
  {"polydispX",     &OS(polydispX),                  1,          1, "%.15G"}, 
  {"polydispY",     &OS(polydispY),                  1,          1, "%.15G"},
  {"polydispZ",     &OS(polydispZ),                  1,          1, "%.15G"},
#else
  {"polydisp",     &OS(polydisp),                  1,          1, "%.15G"}, 
#endif
  {"polycutoff",   &OS(polycutoff),                1,          1, "%.8G"},
#endif
#ifdef MC_KERN_FRENKEL
  {"costhKF",   &OS(costhKF),                       1,         1, "%.8G"},
  {"distKF" ,   &OS(distKF),                        1,         1, "%.8G"},
  {"polylen",   &OS(polylen),                       1,         1, "%d"},
#endif
#ifdef MC_GAPDNA
  {"polylen",   &OS(polylen),                       1,         1, "%d"},
#ifdef GAPDNA_BENDING_ENERGY
  {"kbend",     &OS(kbend),                         1,         1, "%.8G"},
#endif
#endif
#ifdef MC_AMYLOID_FIBRILS
  {"tors_theta0",   &OS(tors_theta0),               1,         1, "%.8G"},
  {"tors_k",        &OS(tors_k),                    1,         1, "%.8G"},
#endif
#ifdef MC_SIMUL
  {"targetPhiMC",  &OS(targetPhiMC),               1,          1, "%.12G"},
#endif
  {"targetPhi",    &OS(targetPhi),                 1,          1, "%.12G"},
#ifdef MC_RESTR_MATRIX
  {"restrmat",   &OS(restrMatrix),                 1,          NAME_LENGTH,  "%s"},
#endif
  {"savedXva",     &OS(savedXva),                   1,   1,   "%d"},
  {"CMreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.14G"},
  {"tolT",         &OS(tolT),                       1,   1, "%.8G"},
  {"useNNL",       &OS(useNNL),                     1,   1, "%d"},
  {"dist5",        &OS(dist5),                      1,   1, "%d"},
  {"dist8stps",    &OS(dist8stps),                  1,   1, "%d"},   
  {"dist5NL",      &OS(dist5NL),                    1,   1, "%d"},
  {"paralNNL",     &OS(paralNNL),                   1,   1, "%d"},
#ifdef MD_GRAVITY
  {"tc",           &OS(tc),                          1,              1, "%.15G"},
  {"quenchtol",    &OS(quenchtol),                  1,   1, "%.10G"},
  {"rhobh",        &OS(rhobh),                  1,   1, "%.10G"},
  {"quenchend",    &OS(quenchend),                  1,   1, "%f"},
  {"taptau",       &OS(taptau),                1,   1, "%f"},
  {"vztap",        &OS(vztap),                 1,   1, "%.15G"},
  {"rzup",         &OS(rzup),                  1,   1, "%.15G"},
  {"expandFact",   &OS(expandFact),            1,   1, "%.15G"},
  {"numquench",    &OS(numquench),             1,   1, "%d"},
  {"maxquench",    &OS(maxquench),             1,   1, "%d"},
  {"extraLz",      &OS(extraLz),                    1,   1, "%.15G"},
  {"checkquechTime",&OS(checkquenchTime),           1,  1,    "%.15G"},
#endif
  {"scalevel",     &OS(scalevel),              1,   1, "%d"},
  {"equilibrated", &OS(equilibrated),          1,   1, "%d"},
#ifdef MD_EDHEFLEX_WALL
  {"epsdPlane",    &OS(epsdPlane),             1,   1, "%.12G"},
  {"epsdFastPlabe", &OS(epsdFastPlane),        1,   1, "%.12G"},
#endif
  {"epsd",         &OS(epsd),                  1,   1, "%.12G"},
  {"epsdNL",       &OS(epsdNL),                1,   1, "%.12G"},
  {"epsdSD",       &OS(epsdSD),                1,   1, "%.12G"},
  {"epsdGDO",      &OS(epsdGDO),               1,   1, "%.12G"},
  {"h",            &OS(h),                     1,   1, "%.15G"},
  {"epsdFast",     &OS(epsdFast),              1,   1, "%.12G"},
  {"epsdFastR",    &OS(epsdFastR),             1,   1, "%.12G"},
  {"epsdMax",      &OS(epsdMax),               1,   1, "%.12G"},
  {"epsdFastNL",     &OS(epsdFastNL),              1,   1, "%.12G"},
  {"epsdFastRNL",    &OS(epsdFastRNL),             1,   1, "%.12G"},
  {"epsdMaxNL",      &OS(epsdMaxNL),               1,   1, "%.12G"},
#ifdef MD_PATCHY_HE
  {"epsdSP",         &OS(epsdSP),                  1,   1, "%.12G"},
  {"epsdFastSP",     &OS(epsdFastSP),              1,   1, "%.12G"},
  {"epsdSPNL",       &OS(epsdSPNL),                1,   1, "%.12G"},
  {"epsdFastSPNL",   &OS(epsdFastSPNL),            1,   1, "%.12G"},
  {"autocat",        &OS(autocat),                 1,   1, "%d"},
  {"k0",             &OS(k0),                      1,   1, "%.12G"},
  {"k1",             &OS(k1),                      1,   1, "%.12G"},
  {"mac",            &OS(mac),                     1,   1, "%.12G"},
#endif
  {"guessDistOpt", &OS(guessDistOpt),          1,   1, "%d"},
  {"springkSD",    &OS(springkSD),              1,   1, "%.12G"},
  {"SDmethod",     &OS(SDmethod),               1,   1, "%d"},
  {"stepSDA",       &OS(stepSDA),                1,   1, "%.12G"},
  {"toldxNR",       &OS(toldxNR),                1,   1, "%.15G"}, 
  {"toldxNRta",     &OS(toldxNRta),              1,   1, "%.15G"},
  {"tolAngNR",      &OS(tolAngNR),               1,   1, "%.15G"},
  {"stepSDB",       &OS(stepSDB),                1,   1, "%.12G"},
  {"maxitsSD",     &OS(maxitsSD),              1,   1, "%d"},
  {"tolSD",        &OS(tolSD),                 1,   1, "%.15G"},         
  {"tolSDlong",    &OS(tolSDlong),             1,   1, "%.15G"},
  {"tolSDconstr",  &OS(tolSDconstr),           1,   1, "%.15G"},
  {"tolSDgrad",    &OS(tolSDgrad),             1,   1, "%.15G"},
  {"tolAngSD",     &OS(tolAngSD),              1,   1, "%.15G"},
  {"forceguess",   &OS(forceguess),            1,   1, "%d"},
  {"zbrakn",       &OS(zbrakn),              1,   1,  "%d"},
  {"zbrentTol",    &OS(zbrentTol),           1,   1,  "%.15G"},
  {"scalfact",     &OS(scalfact),              1,   1, "%.12G"},
  {"reducefact",   &OS(reducefact),            1,   1, "%.12G"},
  {"rescaleTime",  &OS(rescaleTime),                1,  1,    "%.10G"},
  {"phitol",       &OS(phitol),                     1,  1,    "%.14G"},
  {"axestol",      &OS(axestol),                    1,  1,    "%.14G"},
  {"minDist",      &OS(minDist),                    1,  1,    "%.14G"},
  {"rmsd2end",     &OS(rmsd2end),                   1,  1,    "%.6G"},
  {"tmsd2end",     &OS(tmsd2end),                   1,  1,    "%.6G"},
  {"endtime",      &OS(endtime),                    1,  1,    "%.15G"},
  {"nextcheckTime",&OS(nextcheckTime),              1,  1,    "%.15G"},
  {"nextSumTime"  ,&OS(nextSumTime),                1,  1,    "%.15G"},
  {"nextDt",       &OS(nextDt),                     1,  1,    "%.15G"},
#ifdef MD_BIG_DT
  {"refTime",      &OS(refTime),                    1,  1,    "%.15G"},
  {"bigDt",        &OS(bigDt),                      1,  1,    "%.15G"},
#endif
  {"eqlevel",     &OS(eqlevel),                    1,  1,    "%.12G"},
  {"eventMult",    &OS(eventMult),                  1,   1,  "%d"},  
#if defined(MD_INELASTIC) || defined(MD_GRAVITY)
  {"tc",      &OS(tc),                          1,  1,    "%.15G"},
#endif
  {"overlaptol"   ,&OS(overlaptol),                 1,   1, "%f"},
  {"ipart",        &OS(ipart),                      1,   1, "%d"},
  {"brownian",     &OS(brownian),                   1,   1, "%d"},
  {"HNBOX",        &OS(HNBOX),                      1,   1,  "%d"},
  {"avngTemp",     &OS(avngTemp),                   1,   1,  "%d"},
  {"avngPress",    &OS(avngPress),                  1,   1,  "%d"},
  {"avnggr",       &OS(avnggr),                     1,  1,  "%d"},
  {"avngS",        &OS(avngS),                      1,   1,  "%d"},
  //{"avgMedia",     &OS(avgMedia),                   1,     1, "%d"},
  {"xvaSavedMode", &OS(xvaSaveMode),                1,  1,    "%d"},
  {"bakSavedMode", &OS(bakSaveMode),                1,  1,    "%d"},
  {"intervalSum"   ,&OS(intervalSum),               1,  1,    "%.10G"}, 
  {"nextStoreTime", &OS(nextStoreTime),             1, 1,     "%.10G"},
  {"storerate",     &OS(storerate),                 1, 1,     "%.10G"},
  {"KK",            &OS(KK),                        1, 1,     "%d"},
  {"JJ",            &OS(JJ),                        1, 1,     "%d"},
  {"tmpPath",      OS(tmpPath),                     1,  NAME_LENGTH, "%s"},
  {"misPath",      OS(misPath),                     1,  NAME_LENGTH, "%s"},
#ifdef MD_ABSORPTION
  {"bufHeight",    &OS(bufHeight),                 1, 1,     "%.10G"},
#ifdef MD_SPHERICAL_WALL
  {"halfsolidangle",    &OS(halfsolidangle),        1, 1,     "%d"},
#endif
#endif
  {"dofA",         &OS(dofA),                        1,  1, "%d"},
  {"dofB",         &OS(dofB),                        1,  1, "%d"},
  {"base",         &OS(base),                       1,  1, "%.6G"},
#ifdef MD_MULTIPLE_LL
  {"multipleLL",   &OS(multipleLL),                 1,  1, "%d"},
  {"rcutfactMLL",  &OS(rcutfactMLL),                1,  1, "%.15G"},
#endif
#ifdef MD_BILOG
  {"basew",        &OS(basew),                      1,  1, "%.6G"},
  {"lastbilogsaved",&OS(lastbilogsaved),            1,  1, "%d"},
#endif
#ifdef MD_PATCHY_HE
  {"assumeOneBond",     &OS(assumeOneBond),               1,   1, "%d"},
  {"checkGrazing",      &OS(checkGrazing),                1,   1, "%d"},
  {"maxbonds",          &OS(maxbonds),                    1,   1, "%d"}, 
#endif
  {"NN",           &OS(NN),                         1,  1,   "%d"},
  {"fstps",        &OS(fstps),                      1,  1,   "%.15G"},
  {"nRun",         OS(nRun),                        1,  32,   "%s"},
  {"ENmin",        &OS(ENmin),                      1,  1,  "%.6G"},
  {"ENmax",        &OS(ENmax),                      1,  1,  "%.6G"},
  {"PE",           OS(PE),                 PE_POINTS,    1,  "%d"},
#ifdef EDHE_FLEX
  {"par2save",     &OS(par2save),          1, PAR2SAVE_LEN,  "%s"},
  {"n_gauleg",     &OS(n_gauleg),          1, 1,             "%d"},
  {"stripStore",   &OS(stripStore),        1, 1,             "%d"},
  {"Tf",           &OS(Tf),                1, 1,             "%.15G"},
  {"xi",           &OS(xi),                1, 1,             "%.15G"},
#endif
#ifdef ALIGN_POT
  {"Ual",     &OS(Ual),   1, 1, "%f"},
  {"alignaxis", &OS(alignaxis), 1, 1, "%d"},
#endif
#ifdef MD_SUBENZYME
  {"SEreaction",  &OS(SEreaction),          1, 1, "%d"},
  {"SEp",         &OS(SEp),                 1, 1, "%f"},
#endif
#ifdef MD_EDHEFLEX_OPTNNL
  {"optnnl",       &OS(optnnl),            1, 1,              "%d"},
#endif
#ifdef MC_BOUNDING_SPHERES
  {"useboundsph",   &OS(useboundsph), 1, 1, "%d"},
#endif
#if defined(MD_RABBIT) || defined(MD_NANOBODY)
  {"time_limit",           &OS(time_limit),                1, 1,             "%.15G"},
  {"first_time",           &OS(first_time),                1, 1,             "%.15G"},
  {"rhozBinSize",          &OS(rhozBinSize),               1, 1,             "%.15G"},
  {"rate",                 &OS(rate),                      10, 1,            "%f"},
#endif
#ifdef MD_SUBENZYME
  {"rateSE",               &OS(rateSE),                    10, 1,            "%f"},
#endif
#ifdef MD_PROTEIN_DESIGN
  {"nativeConf",           &OS(nativeConf),                  1, NAME_LENGTH, "%s"},
#endif
#ifdef MC_SIMUL
  {"restrmove",      &OS(restrmove),                     1,  1, "%d"},
#ifdef MC_FREEZE_BONDS
  {"freezebonds",  &OS(freezebonds),                     1, 1, "%d"},
#endif
#ifdef MC_NVE
  {"Ed", &OS(Ed), 1, 1, "%.12G"},
#endif
#ifdef MC_CLUSTER_MOVE
  {"clsmovprob",      &OS(clsmovprob),                  1, 1, "%.12G"},
  {"delTclsMC",       &OS(delTclsMC),                   1, 1, "%.12G"},
  {"delRclsMC",       &OS(delRclsMC),                   1, 1, "%.12G"},
#endif
#ifdef MC_BIGROT_MOVE
  {"bigrotmov",       &OS(bigrotmov),                   1, 1, "%.12G"},
#ifdef MC_BIGROT_BIASED
  {"bigrotbias",       &OS(bigrotbias),                   1, 1, "%.12G"},
  {"bigrotTheta0",     &OS(bigrotTheta0),                   1, 1, "%.12G"},
#endif
#endif
#ifdef MC_GRANDCAN
  {"zetaMC",     &OS(zetaMC),                            1,  1, "%.12G"},
  {"npav" ,     &OS(npav),                             1,  1, "%d"},
  {"nexc" ,     &OS(nexc),                             1,  1, "%d"},
#ifdef MC_SUS
  {"susnmin",     &OS(susnmin),                    1, 1, "%d"},
  {"susnmax",     &OS(susnmax),                    1, 1, "%d"},
#endif
#if defined(MC_SWHC) || defined(MC_SWELL)
  {"deltasw",         OS(deltasw),                   2,               1, "%.12G"},
  {"constDelta",      &OS(constDelta),                1,               1, "%d"},
#endif
#if defined(MC_HYDROPHOBIC_INT) 
  {"maxtrialsH",         &OS(maxtrialsH),                   1,               1, "%d"},
#endif
#ifdef MC_HELIX
  {"Nxi", &OS(Nxi),  1, 1, "%d"},
  {"pitch", &OS(pitch), 1, 1, "%.8G"},
  {"sighelix", &OS(sighelix), 1, 1, "%.8G"},
  {"lenhelix", &OS(lenhelix), 1, 1, "%.8G"},
  {"radhelix", &OS(radhelix), 1, 1, "%.8G"},
#endif
#endif
#ifdef MC_FLIP_MOVE
  {"flip_prob",    &OS(flip_prob),                1, 1,   "%.12G"},
#endif
  {"nvbbias",   &OS(nvbbias),                           1, 1,  "%d"},
  {"vbond",     &OS(vbond),                            1,  1, "%.15G"},
  {"pbias",     &OS(pbias),                            1,  1, "%.15G"},
  {"lastNNLrebuildMC", &OS(lastNNLrebuildMC),                1, 1, "%d"},
  {"targetAccept", &OS(targetAccept),                        1, 1, "%.15G"},
  {"targetAcceptvol", &OS(targetAcceptVol),                  1, 1, "%.15G"},
  {"ensembleMC",   &OS(ensembleMC),                        1,  1, "%d"},
#ifdef MC_ALMARZA
  {"almarza_thr",        &OS(almarza_thr),                1,  1, "%.12G"},        
#endif
  {"dthetaMC",     &OS(dthetaMC),                            1,  1, "%.12G"},
  {"deltaMC" ,     &OS(deltaMC),                             1,  1, "%.12G"},
  {"vmax",        &OS(vmax),                                1,  1, "%.12G"},        
  {"resetaccept", &OS(resetaccept),                         1, 1, "%d"},
  {"resetacceptVol", &OS(resetacceptVol),                   1, 1, "%d"},
  {"outMC",          &OS(outMC),                            1, 1, "%d"},
  {"adjstepsMC",     &OS(adjstepsMC),                       1, 1, "%d"},
#endif
#ifdef MD_SURV_PROB
  {"spdeltat",          &OS(spdeltat),                       1,   1, "%.15G"},
#endif
  {"", NULL, 0, 0, ""}
};
#else
extern struct pascii opro_ascii[];
#endif

#ifdef MAIN
/* ========================= >>> opar_ascii <<< =========================== */
struct pascii opar_ascii[]=
{
  {"parnum",            &OP(parnum),                      1,     1, "%d"},
#ifndef EDHE_FLEX
  {"parnumA",           &OP(parnumA),                      1,     1, "%d"},
#endif
  {"totStep",           &OP(totStep),                     1,   1,  MDINTFMT},
  {"time",              &OP(time),                        1,   1, "%.15G"},
  {"curStep",           &OP(curStep),                     1,   1,  MDINTFMT},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
#ifndef EDHE_FLEX
  {"m",                 OP(m),                            2,   1, "%.6G"},
  {"a",                 OP(a),                             2,   1, "%.8G"},
  {"b",                 OP(b),                             2,   1, "%.8G"},
  {"c",                 OP(c),                             2,   1, "%.8G"},
#endif
#ifdef MD_GHOST_IGG
  {"ghostsim",          &OP(ghostsim),                   1,               1, "%d"},
#endif
  {"rcut",              &OP(rcut),                        1,   1, "%.10G"},
#ifdef MC_BOND_POS
  {"rcutBP",              &OP(rcutBP),                        1,   1, "%.10G"},
#endif
  {"equilibrat",        &OP(equilibrat),                  1,   1,   "%d"},
  {"Dt",                &OP(Dt),                          1,   1, "%.15G"},
#ifndef EDHE_FLEX
#ifndef MD_ASYM_ITENS
  {"I",                 OP(I),                             2,   1, "%.8G"},
#else
  {"I",                 OP(I),                             2,   3, "%.8G"},
#endif
#endif
#ifdef MD_GRAVITY
  {"wallDiss",          &OP(wallDiss),                    1,   1,   "%f"},
  //{"partDiss",          &OP(partDiss),                    1,   1,   "%f"},
  {"ggrav",             &OP(ggrav),                       1,   1,   "%f"},
#endif
#if defined(MD_INELASTIC) || defined(MD_GRAVITY)
  {"partDiss",          &OP(partDiss),                    1,   1,   "%f"},
#endif
  {"M",                 &OP(M),                           1,   1,   "%d"},
  {"tol",               &OP(tol),                         1,   1, "%.15G"},
#ifdef EDHE_FLEX
  {"ninters",       &OP(ninters),                         1,  1,    "%d"},
  {"nintersIJ",     &OP(nintersIJ),                       1,  1,    "%d"},
  {"ntypes",        &OP(ntypes),                          1,  1,    "%d"},
  {"saveBonds",      &OP(saveBonds),                      1,  1,    "%d"},
  {"maxbondsSaved",  &OP(maxbondsSaved),                  1,  1,    "%d"},
#endif
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
  {"sigmaSticky",       &OP(sigmaSticky),                       1,   1, "%.15G"},
  {"bheight",           &OP(bheight),                     1,   1, "%.15G"},
  {"bhin",               &OP(bhin),                         1,   1, "%.15G"},
  {"bhout",              &OP(bhout),                         1,   1, "%.15G"},
  {"Dr",                 &OP(Dr),                            1,   1, "%.15G"},
  {"theta",              &OP(theta),                         1,   1, "%.15G"},
  {"nmax",               &OP(nmax),                          1,   1, "%d"},
#endif
 {"", NULL, 0, 0, ""}
};
#else
extern struct pascii opar_ascii[];
#endif

/* ------------------------ PARTICLES COORDINATES -------------------- */

/* MAIN is a macro defined only in mdsimul.c, so there we put the instance
of the object, and elsewhere we put only extern declarations */
#ifdef MAIN
COORD_TYPE DECL_LIST;
COORD_TYPE EXT_DLST;
int maxcoll=-1;
#ifdef EDHE_FLEX
int MD_DECL_INT;
partType* typesArr;
interStruct* intersArr;
int *typeNP, *typeOfPart;
interStructIJ* intersArrIJ;
#endif
#else
extern int maxcoll;
extern COORD_TYPE DECL_LIST; 
extern COORD_TYPE EXT_DLST;
#ifdef EDHE_FLEX
extern int MD_DECL_INT;
extern partType* typesArr; /* array con tutti i tipi presenti nella simulazione */
extern interStruct* intersArr; /* array di strutture contenente tutte le interazioni */
extern interStructIJ* intersArrIJ;
extern int *typeNP, *typeOfPart; /* array contentente il numero di particelle di ogni specie */
#endif
#endif 

/* ======================== >>> singlePar Array <<< ========================*/

/* MAIN is a macro defined only in mdsimul.c, so there we put the instance
of the object, and elsewhere we put only extern declarations */
#ifdef MAIN
struct singlePar OsinglePar[] = { 
  /* =================== >>> DON'T TOUCH THESE !!!! <<< =================== 
     The first string of each array elemen ( singlePas.parName ), is the
     string you can use in the parameters file to set its value.
     For example considering the definition below of the array me could in the
     parameters file the line:
        parnum: 100 
     It then put at begin of the simulation the value 10 in Opramas.parnum/
     Every parameter could be of this types:
     - INT = integer ( Linux-> 4 bytes signed integer )
     - STR = string NAME_LENGTH bytes long (see mdsimul.h)
     - CT = COORD_TYPE 
     WARNING: the params types should be consistent with previous params 
       structure declaration.
  */
  {"parnum" ,    &Oparams.parnum,             INT},
#ifndef EDHE_FLEX
  {"parnumA" ,   &Oparams.parnumA,            INT},
#endif
#ifdef EDHE_FLEX
  {"ntypes",     &Oparams.ntypes,             INT},
  {"ninters",    &Oparams.ninters,            INT},
  {"nintersIJ",  &Oparams.nintersIJ,          INT},
#endif
  {"stepnum",    &Oparams.totStep,            LLINT},
  {"inifile" ,   &OprogStatus.inifile,        STR},
  {"endfile" ,   &OprogStatus.endfile,        STR},
  {"xvafile" ,   &OprogStatus.xvafile,        STR},
  {"inistep" ,   &Oparams.curStep,            LLINT},
  {"endFormat",  &OprogStatus.endFormat,      INT},
#ifdef MC_SIMUL
  {"restrmove",  &OprogStatus.restrmove,    INT},
#ifdef MC_FREEZE_BONDS
  {"freezebonds", &OprogStatus.freezebonds,  INT},
#endif
#ifdef MC_NVE
  {"Ed",         &OprogStatus.Ed,  CT},
#endif
#ifdef MC_CLUSTER_MOVE
  {"clsmovprob", &OprogStatus.clsmovprob,   CT},
  {"delTclsMC",  &OprogStatus.delTclsMC,    CT},
  {"delRclsMC",  &OprogStatus.delRclsMC,    CT},
#endif
#ifdef MC_BIGROT_MOVE
  {"bigrotmov", &OprogStatus.bigrotmov, CT},
#ifdef MC_BIGROT_BIASED
  {"bigrotbias", &OprogStatus.bigrotbias, CT},
  {"bigrotTheta0", &OprogStatus.bigrotTheta0, CT},
#endif
#endif
#ifdef MC_ELCONST_MC
  {"alpha", &OprogStatus.alpha, CT},
  {"lp",    &OprogStatus.lp,    CT},
  {"polylen", &OprogStatus.polylen, INT},
  {"eqstps",  &OprogStatus.eqstps, INT},
  {"calcvexcl", &OprogStatus.calcvexcl, INT},
#endif
#ifdef MC_GRANDCAN
  {"zetaMC",   &OprogStatus.zetaMC,         CT},
  {"npav",  &OprogStatus.npav,        INT},
  {"nexc",&OprogStatus.nexc,        INT},
#ifdef MC_SUS
  {"susnmin", &OprogStatus.susnmin, INT},
  {"susnmax", &OprogStatus.susnmax, INT},
#endif
#if defined(MC_SWHC) || defined(MC_SWELL)
  {"deltaswL2", &OprogStatus.deltasw[0], CT},
  {"deltaswD2", &OprogStatus.deltasw[1], CT},
  {"constDelta", &OprogStatus.constDelta, INT},
#endif
#ifdef MC_HYDROPHOBIC_INT
  {"maxtrialsH", &OprogStatus.maxtrialsH, INT},
#endif
#ifdef MC_HELIX
  {"Nxi", &OprogStatus.Nxi, INT},
  {"pitch", &OprogStatus.pitch, CT},
  {"sighelix", &OprogStatus.sighelix, CT},
  {"radhelix", &OprogStatus.radhelix, CT},
  {"lenhelix", &OprogStatus.lenhelix, CT},
#endif
#endif
#ifdef MC_FLIP_MOVE
  {"flip_prob",  &OprogStatus.flip_prob,     CT},
#endif
  {"nvbbias",  &OprogStatus.nvbbias,          INT},
  {"pbias",    &OprogStatus.pbias,            CT},
  {"deltaMC",   &OprogStatus.deltaMC,         CT},
  {"dthetaMC",  &OprogStatus.dthetaMC,        CT},
  {"ensembleMC",&OprogStatus.ensembleMC,      INT},
#ifdef MC_ALMARZA
  {"almarza_thr", &OprogStatus.almarza_thr,    CT},
#endif
  {"vmax",      &OprogStatus.vmax,            CT},
  {"resetaccept",  &OprogStatus.resetaccept,    INT},
  {"resetacceptVol", &OprogStatus.resetacceptVol, INT},
  {"outMC",          &OprogStatus.outMC,      INT},
  {"adjstepsMC",     &OprogStatus.adjstepsMC, INT},
  {"targetAccept",   &OprogStatus.targetAccept,           CT},
  {"targetAcceptVol",&OprogStatus.targetAcceptVol,        CT},
  {"vbond",          &OprogStatus.vbond,                  CT},
#endif
#ifdef EDHE_FLEX
#ifdef POVRAY
  {"povpcol",  &OprogStatus.povpcol, STR},
  {"povcolA",  &OprogStatus.povcolA, STR},
  {"povcolB",  &OprogStatus.povcolB, STR},
  {"povtransmit", &OprogStatus.povtransmit, CT},
  {"povcx"      , &OprogStatus.povcx, CT},
  {"povcy"      , &OprogStatus.povcy, CT},
  {"povcz"      , &OprogStatus.povcz, CT},
  {"povlx"      , &OprogStatus.povlx, CT},
  {"povly"      , &OprogStatus.povly, CT},
  {"povlz"      , &OprogStatus.povlz, CT},
#endif
  {"optbm",      &OprogStatus.optbm,          INT},
#ifdef MD_SCALEPHI_STAGES
  {"growthType",       &OprogStatus.growthType,  INT},
#endif
  {"frozenDOF",       &OprogStatus.frozenDOF,  INT},
#endif
#ifdef MD_ABSORPTION
  {"bufHeight", &OprogStatus.bufHeight, CT},
#ifdef MD_SPHERICAL_WALL
  {"halfsolidangle", &OprogStatus.halfsolidangle, INT},
#endif
#endif
#ifdef MD_EDHEFLEX_WALL
  {"hardwall",   &OprogStatus.hardwall,       INT},
#endif
#ifdef MD_GHOST_IGG
  {"ghostsim",   &Oparams.ghostsim,       INT},
#endif
  {"iniFormat",  &OprogStatus.iniFormat,      INT},
  {"tapeTimes",  &OprogStatus.tapeTimes,      LLINT},
  {"bakSteps",   &OprogStatus.bakSteps,       LLINT},
  {"bakStepsAscii", &OprogStatus.bakStepsAscii,LLINT},
  {"staSteps",   &OprogStatus.staSteps,       LLINT},
  {"xvaSteps",   &OprogStatus.xvaSteps,       LLINT},  
  {"energySteps",&OprogStatus.measSteps[0],   LLINT},
  {"energyCalc", &OprogStatus.measCalc[0],    LLINT},
  {"energyName", &OprogStatus.dataFiles[0],   STR},
  {"energyBegin",&OprogStatus.initStep[0],   LLINT},
  {"DtrSteps",   &OprogStatus.measSteps[1],   LLINT}, /* steps between measure
						       savings */
  {"DtrCalc",    &OprogStatus.measCalc[1],    LLINT}, /* steps between measure
						       calculation */
  {"DtrName",    &OprogStatus.dataFiles[1],   STR},
  {"brownian",   &OprogStatus.brownian,       INT},
  {"tempSteps",  &OprogStatus.measSteps[2],   LLINT},
  {"tempCalc",   &OprogStatus.measCalc[2],    LLINT},
  {"tempName",   &OprogStatus.dataFiles[2],   STR},
  {"grSteps",    &OprogStatus.measSteps[3],   LLINT},
  {"grCalc",     &OprogStatus.measCalc[3],    LLINT},
  {"grName",     &OprogStatus.dataFiles[3],   STR},
  {"MBSteps",    &OprogStatus.measSteps[4],   LLINT},
  {"MBCalc",     &OprogStatus.measCalc[4],    LLINT},
  {"MBName",     &OprogStatus.dataFiles[4],   STR},
  {"SSteps",     &OprogStatus.measSteps[5],   LLINT},
  {"SCalc",      &OprogStatus.measCalc[5],    LLINT},
  {"SName",      &OprogStatus.dataFiles[5],   STR},
  {"pressSteps", &OprogStatus.measSteps[6],   LLINT},
  {"pressCalc",  &OprogStatus.measCalc[6],    LLINT},   
  {"pressName",  &OprogStatus.dataFiles[6],   STR},
  {"PtensSteps", &OprogStatus.measSteps[7],  LLINT},
  {"PtensCalc",  &OprogStatus.measCalc[7],   LLINT},   
  {"PtensName",  &OprogStatus.dataFiles[7],  STR},
  {"DQtensSteps",&OprogStatus.measSteps[8],  LLINT},
  {"DQtensCalc", &OprogStatus.measCalc[8],   LLINT},
  {"DQtensName", &OprogStatus.dataFiles[8],  STR},
  {"VSteps",  &OprogStatus.measSteps[9],   LLINT},
  {"VCalc",   &OprogStatus.measCalc[9],    LLINT},   
  {"VName",   &OprogStatus.dataFiles[9],   STR},
  {"rotMSDSteps",  &OprogStatus.measSteps[10],   LLINT},
  {"rotMSDCalc",   &OprogStatus.measCalc[10],    LLINT},   
  {"rotMSDName",   &OprogStatus.dataFiles[10],   STR},
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
  {"useNNL"    , &OprogStatus.useNNL,         INT},
  {"dist5",    &OprogStatus.dist5,          INT},
  {"dist8stps",&OprogStatus.dist8stps,      INT},
  {"dist5NL",    &OprogStatus.dist5NL,        INT},
  {"paralNNL",   &OprogStatus.paralNNL,       INT},
  {"overlaptol", &OprogStatus.overlaptol,     CT},
  {"intervalSum", &OprogStatus.intervalSum,   CT},
  {"rescaleTime", &OprogStatus.rescaleTime,   CT},
  {"storerate",     &OprogStatus.storerate,     CT},
  {"scalevel",   &OprogStatus.scalevel,       INT},
  {"endtime",    &OprogStatus.endtime,        CT},
  {"Dt",         &Oparams.Dt,                 CT},
#ifdef MD_BIG_DT
  {"bigDt",      &OprogStatus.bigDt,           CT},
#endif
#ifdef MD_EDHEFLEX_WALL
  {"epsdPlane",       &OprogStatus.epsdPlane,      CT},
  {"epsdFastPlane",   &OprogStatus.epsdFastPlane,  CT},
#endif
  {"epsd",       &OprogStatus.epsd,           CT},
  {"epsdNL",     &OprogStatus.epsdNL,         CT},
  {"epsdSD",     &OprogStatus.epsdSD,         CT},
  {"epsdGDO",    &OprogStatus.epsdGDO,        CT},
  {"h",          &OprogStatus.h,              CT},
  {"epsdFast",   &OprogStatus.epsdFast,       CT},
  {"epsdFastR",  &OprogStatus.epsdFastR,      CT},
  {"epsdMax",    &OprogStatus.epsdMax,        CT},
  {"epsdFastNL",   &OprogStatus.epsdFastNL,       CT},
  {"epsdFastRNL",  &OprogStatus.epsdFastRNL,      CT},
  {"epsdMaxNL",    &OprogStatus.epsdMaxNL,        CT},
#ifdef MD_PATCHY_HE
  {"epsdSP",       &OprogStatus.epsdSP,           CT},
  {"epsdFastSP",   &OprogStatus.epsdFastSP,       CT},
  {"epsdSPNL",     &OprogStatus.epsdSPNL,         CT},
  {"epsdFastSPNL", &OprogStatus.epsdFastSPNL,     CT},
  {"autocat",      &OprogStatus.autocat,          INT},
  {"k0",           &OprogStatus.k0,               CT},
  {"k1",           &OprogStatus.k1,               CT},
  {"mac",          &OprogStatus.mac,              CT},
#endif
  {"guessDistOpt",&OprogStatus.guessDistOpt,  INT},
  {"tolSD",      &OprogStatus.tolSD,          CT},
  {"tolSDlong",  &OprogStatus.tolSDlong,      CT},
  {"tolSDconstr",&OprogStatus.tolSDconstr,    CT},
  {"tolSDgrad",  &OprogStatus.tolSDgrad,      CT},
  {"tolAngSD",   &OprogStatus.tolAngSD,       CT},
  {"forceguess", &OprogStatus.forceguess,     INT},
  {"springkSD",  &OprogStatus.springkSD,    CT},
  {"SDmethod",    &OprogStatus.SDmethod,     INT},
  {"toldxNR",    &OprogStatus.toldxNR,  CT},
  {"toldxNRta",  &OprogStatus.toldxNRta,  CT},
  {"tolAngNR",   &OprogStatus.tolAngNR, CT},
#ifdef EDHE_FLEX
  {"stepSD",      &OprogStatus.stepSDA,         CT},
#else
  {"stepSDA",     &OprogStatus.stepSDA,         CT},
  {"stepSDB",     &OprogStatus.stepSDB,         CT},
#endif
  {"maxitsSD",   &OprogStatus.maxitsSD,       INT},
  {"zbrakn",     &OprogStatus.zbrakn,         INT},
  {"zbrentTol",  &OprogStatus.zbrentTol,      CT},
  {"scalfact",   &OprogStatus.scalfact,       CT},
  {"reducefact", &OprogStatus.reducefact,     CT},
  {"phitol",      &OprogStatus.phitol,        CT},
#ifdef MD_CALENDAR_HYBRID
  {"scaleHQ",  &OprogStatus.scaleHQ,          CT},
  {"nlistsHQ", &OprogStatus.nlistsHQ,         INT},
  {"adjustHQ", &OprogStatus.adjustHQ,         INT},
  {"overthrHQ",&OprogStatus.overthrHQ,        INT},
#endif
  {"axestol",     &OprogStatus.axestol,       CT},
  {"minDist",     &OprogStatus.minDist,       CT},
  {"tmsd2end",    &OprogStatus.tmsd2end,      CT},
  {"rmsd2end",    &OprogStatus.rmsd2end,      CT},
#ifdef MD_GRAVITY
  {"taptau",     &OprogStatus.taptau,         CT},
  {"quenchtol",  &OprogStatus.quenchtol,      CT},
  {"checkQuench", &OprogStatus.checkquenchTime, CT},
  {"Lz",         &Lz,                       CT},
  {"extraLz",    &OprogStatus.extraLz,      CT},
  {"vztap",      &OprogStatus.vztap,        CT},
  {"rzup",       &OprogStatus.rzup,         CT},
  {"rhobh",      &OprogStatus.rhobh,        CT},
  {"expandFact",  &OprogStatus.expandFact,   CT},
#endif
  {"W",          &OprogStatus.W,              CT},
  {"P",          &Oparams.P,                  CT},
#ifdef MD_LXYZ
  {"Lx",          &L[0],                     CT},
  {"Ly",          &L[1],                     CT},
  {"Lz",          &L[2],                     CT},
#else
  {"L",          &L,                        CT},
#endif
#ifdef MD_PATCHY_HE
  {"sigmaSticky", &Oparams.sigmaSticky,     CT},
  {"bheight",    &Oparams.bheight,          CT},
  {"bhin",         &Oparams.bhin,              CT},
  {"bhout",       &Oparams.bhout,             CT},
  {"assumeOneBond", &OprogStatus.assumeOneBond, INT},
  {"checkGrazing",  &OprogStatus.checkGrazing, INT},
  {"maxbonds",      &OprogStatus.maxbonds,     INT},
#ifdef EDHE_FLEX
  {"saveBonds",     &Oparams.saveBonds,        INT},
#endif
  {"nmax",          &Oparams.nmax,         INT},
  {"Dr",            &Oparams.Dr,           CT},
  {"theta",         &Oparams.theta,        CT},
#endif
  {"avngTemp",   &OprogStatus.avngTemp,       INT},
  {"avngPress",  &OprogStatus.avngPress,      INT},
  {"avngS",      &OprogStatus.avngS,          INT},
  {"avnggr",     &OprogStatus.avnggr,         INT},
  {"avngMB",     &OprogStatus.avngMB,         INT},
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
#ifndef EDHE_FLEX
  {"A0",      &Oparams.a[0],      CT},
  {"A1",      &Oparams.a[1],      CT},
  {"B0",      &Oparams.b[0],      CT},
  {"B1",      &Oparams.b[1],      CT},
  {"C0",      &Oparams.c[0],      CT},
  {"C1",      &Oparams.c[1],      CT},
#endif
  {"targetPhi", &OprogStatus.targetPhi, CT},
#ifdef MC_RESTR_MATRIX
  {"restrmat",  &OprogStatus.restrMatrix, STR},
#endif
#ifdef MC_KERN_FRENKEL
  {"costhKF",   &OprogStatus.costhKF, CT},
  {"distKF",    &OprogStatus.distKF,  CT},
  {"polylen",   &OprogStatus.polylen, INT},
#endif
#ifdef MC_GAPDNA
  {"polylen",   &OprogStatus.polylen, INT},
#ifdef GAPDNA_BENDING_ENERGY
  {"kbend",     &OprogStatus.kbend,   CT},
#endif
#endif
#ifdef MC_AMYLOID_FIBRILS
  {"tors_theta0", &OprogStatus.tors_theta0,  CT},
  {"tors_k",       &OprogStatus.tors_k,       CT},
#endif
#ifdef MC_SIMUL
  {"targetPhiMC", &OprogStatus.targetPhiMC, CT},
#endif
#ifdef MD_POLYDISP
#ifdef MD_POLYDISP_XYZ
  {"polydispX",  &OprogStatus.polydispX, CT},  
  {"polydispY",  &OprogStatus.polydispY, CT},
  {"polydispZ",  &OprogStatus.polydispZ, CT},
#else
  {"polydisp",  &OprogStatus.polydisp, CT},  
#endif
  {"polycutoff",&OprogStatus.polycutoff, CT},
#endif
#ifndef EDHE_FLEX
#ifndef MD_ASYM_ITENS
  {"Ia",      &Oparams.I[0],      CT},
  {"Ib",      &Oparams.I[1],      CT},
#else
  {"I1a",      &Oparams.I[0][0],      CT},
  {"I3a",      &Oparams.I[0][2],      CT},
  {"I1b",      &Oparams.I[1][0],      CT},
  {"I3b",      &Oparams.I[1][2],      CT},
#endif
#endif
#ifdef MD_GRAVITY
  {"ggrav",      &Oparams.ggrav,            CT},
#endif
#ifndef EDHE_FLEX
  {"mass0",       &Oparams.m[0],                CT},
  {"mass1",       &Oparams.m[1],                CT},
#endif
#ifdef MD_GRAVITY
  {"wallDiss",   &Oparams.wallDiss,         CT},
  {"quenchend",  &OprogStatus.quenchend,    CT},
  {"maxquench",  &OprogStatus.maxquench,    INT},
#endif
#if defined(MD_INELASTIC) || defined(MD_GRAVITY)
  {"partDiss",   &Oparams.partDiss,         CT},
  {"tc",         &OprogStatus.tc,           CT},
#endif
  {"eventMult",  &OprogStatus.eventMult,    INT},
  {"rcut",       &Oparams.rcut,             CT},
#ifdef MC_BOND_POS
  {"rcutBP",       &Oparams.rcutBP,         CT},
#endif
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"eqlevel",    &OprogStatus.eqlevel,       CT},
  {"temperat",   &Oparams.T,                CT},
  {"tol",        &Oparams.tol,              CT},
  {"seed",       &mdseed,                   INT},
  {"dofA",       &OprogStatus.dofA,         INT},
  {"dofB",       &OprogStatus.dofB,         INT},
  /* parametri per scegliere il tipo di salvataggio del file xva
     (lineare o semilog) */
  {"xvaSaveMode",&OprogStatus.xvaSaveMode,  INT},
  {"bakSaveMode", &OprogStatus.bakSaveMode, INT},
  {"tmpPath",    OprogStatus.tmpPath,        STR},
  {"misPath",    OprogStatus.misPath,       STR},
  {"base",       &OprogStatus.base,         CT},
#ifdef MD_MULTIPLE_LL
  {"multipleLL", &OprogStatus.multipleLL,   INT},
  {"rcutfactMLL", &OprogStatus.rcutfactMLL, CT},
#endif
  {"NN",         &OprogStatus.NN,           INT},
#ifdef MD_BILOG
  {"basew",     &OprogStatus.basew,         CT},
#endif
#ifdef EDHE_FLEX
  {"par2save",   &OprogStatus.par2save,    STR},
  {"stripStore",  &OprogStatus.stripStore,  INT},
  {"Tf",          &OprogStatus.Tf,          CT},
  {"xi",          &OprogStatus.xi,          CT},
  {"n_gauleg",    &OprogStatus.n_gauleg,    INT},
#endif
#ifdef ALIGN_POT
  {"Ual",          &OprogStatus.Ual,        CT},
  {"alignaxis",    &OprogStatus.alignaxis,  INT},
#endif
#ifdef MD_SUBENZYME
  {"SEreaction",   &OprogStatus.SEreaction, INT},
  {"SEp",          &OprogStatus.SEp,        CT},
#endif
#ifdef MC_BOUNDING_SPHERES
  {"useboundsph",   &OprogStatus.useboundsph, INT},
#endif
#ifdef MD_EDHEFLEX_OPTNNL
  {"optnnl",      &OprogStatus.optnnl,      INT},
#endif
  {"ENmin",      &OprogStatus.ENmin,        CT},
  {"ENmax",      &OprogStatus.ENmax,        CT},
  {"nRun",       &OprogStatus.nRun,         STR},
#if defined(MD_RABBIT) || defined(MD_NANOBODY)
  {"time_limit",  &OprogStatus.time_limit,   CT},
  {"rhozBinSize", &OprogStatus.rhozBinSize,  CT},
#endif
  {"maxcoll",     &maxcoll,                 INT},
#ifdef MD_SURV_PROB
  {"spdeltat",    &OprogStatus.spdeltat,    CT},
#endif
  /* ======================================================================= */

  {"", NULL, 0} /* end of list, don't touch this !!! */
};
#else 
extern struct singlePar OsinglePar[];
#endif 
/* ------------------------------ MEASURES --------------------------------*/
 /* initialize all Omeasure objects, this objects are used to store measures
    on disk and are analogous to OsinglePar objects */

/* MEASURE VARIABLE DECLARATION AND ALL OTHER VARIBALES  YOU NEED 
   Put here all the variables for your measure calculation 
   and if you want a variable on disk, that is if you want to make a measure
   put pointer to this variable int the Omeasure structure.
   Note that you must write this variable in two copies, with one precedeed 
   by the specifier 'extern' */

/* TODO:
   These could be put in OprogStatus so the user should give them in the
   parameter file */

#define VBEG 0.0   /*Interval of velocities for the maxwellian */ 
#define VEND 6.5 
#define KBEG 0.63
#define KEND 25.0
#define RBEG 0.2
#define REND 5.65
#define HNBOX 0
#define NUMK2AV 150 /* Number of k on the same modulus over which we must 
		       perform a mean */

#ifdef MAIN
#ifdef MD_GHOST_IGG
ghostInfo *ghostInfoArr=NULL;
#endif
COORD_TYPE E, Dtrans, temp, S[NUMK], dummy, eta, gr[MAXBIN], invs, press,
  press_m, press_at, rcmz, rho, ItensD[2][3], pressST, pressHS, pressKin;
COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, Aa, V, DrSqTot, temp_transl, DphiSq;
int globSaveAll=1, MB[NUMV];
#else 
#ifdef MD_GHOST_IGG
extern ghostInfo *ghostInfoArr;
#endif
extern COORD_TYPE E, Dtrans, temp, S[NUMK], dummy, eta, gr[MAXBIN], invs,
  press, press_m, press_at, temp_transl, rcmz, rho;
extern COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, V, Aa, DrSqTot,
  DphiSq, ItensD[2][3], DphiSq, pressST, pressHS, pressKin;
extern int MB[NUMV];
extern int globSaveAll;
#endif

/* ============= >>> PUT HERE MEASURING FUNCTION PROTOTYPES <<< ============*/

void energy(void);
void calcV(void);
void transDiff(void);
void temperat(void);
void structFacts(void);
void maxwell(void);
void rotDiff(void);
void viscosity(void);
void radDens(void);
void Ptensor(void);
void DQtensor(void);
void vanHoveSelf(void);
void intermScatt(void);
void gaussapprox(void);
void calcrotMSD(void);
/* =========================================================================*/
#ifdef MAIN
struct measure Omeasure[NUM_MISURE]=
{
  /* DESCRIPTION:
     {<1>, <2>, <3>} where:
     <1> = pointer to the single measure, that could a variable, an array, 
     struct or whatever you want. 
     <2> = size of each measure ( variable, struct, array, etc.)
     <3> = pointer to the function that performs the calculation of the 
           measure ( NULL means that there is no need for this function )*/
  /* ==================== >>> PUT HERE YOUR MEASURES <<< =================== */
  
  {&E,      sizeof(COORD_TYPE),       energy      },
  {&Dtrans, sizeof(COORD_TYPE),       transDiff   },
  {&temp,   sizeof(COORD_TYPE),       temperat    },
  {gr,      sizeof(COORD_TYPE)*MAXBIN,radDens     },
  {MB,      sizeof(int)*NUMV,         maxwell     },
  {S,       sizeof(COORD_TYPE)*NUMK,  structFacts },
  {&press,  sizeof(COORD_TYPE),       NULL        },
  {Ptens,   sizeof(COORD_TYPE)*3,     Ptensor     },
  {DQtens,  sizeof(COORD_TYPE)*3,     DQtensor    },
  {&V,      sizeof(double),           calcV       },
  {&DphiSq, sizeof(double),           calcrotMSD    },
  /* ======================================================================= */
  {NULL, 0, NULL}                   /* End of list don't touch this */
};
#else 
extern struct measure Omeasure[];
#endif

/*============================= >>> measHead <<< ==========================*/
struct measHead
{
  /* DESCRIPTION:
     This the Header put at the top of each mesure file ( see TECH_INFO ) */
  /* ADD 14/09/2000 */
  double dt;
  double T;
  double Vol;
  int N;
  MDINT saveSteps;
  int size;  /* size in bytes of each measure */
};

/* ===================== >>> measHead Instanziation <<< =====================*/
#ifdef MAIN
struct measHead OmeasHead[NUM_MISURE];
#else 
extern struct measHead OmeasHead[NUM_MISURE];
#endif 


/* ======================= >>> xvaHead structure <<< ========================*/
struct xvaHead 
{
  /* DESCRIPTION:
     This strcture is put at begin of the tape file (xva file) */
  int size ;     /* size of each xva savings = 
			     <number coordinate in XVA_LIST> * SEGSIZE */
  int mode;
  MDINT saveSteps; /* save between two tape savings */
  int NN;
  double T;
  double Vol;
  double base;
  double dt;     /* dt of each steps */
  int parnum;
};

/* ===================== >>> xvaHead Instanziation <<< =====================*/
#ifdef MAIN
struct xvaHead OxvaHead;
#else 
extern struct xvaHead OxvaHead;
#endif 


/* ==================== >>> convStruct Instantiation <<< ====================*/
#ifdef MD2ASCII         /* this macro is defined only in md2ascii.c */

/* ======================= >>> VARIABLES <<< ===============================
   THIS ARE NEEDED BY YOUR CONVERTERS (see md2ascii.c for examples of 
   converters)*/
char precision[64];         
/* precision used to output floating point numbers (it is a string containing
   a number)*/
void* mis;
/* pointer to a buffer to store every measure read from the input file */

/* =============== >>> USER DEFINED CONVERTERS CODE <<< =====================*/
void writeVec(FILE* afs, MDINT step, int size);/*Predefined*/
void writeFp(FILE* afs, MDINT step, int size); /*Predefined*/
void writeInt(FILE* afs, MDINT step, int size);/*Predefined*/
void writeS(FILE*afs, MDINT step, int size);
void writeg(FILE* afs, MDINT step, int size);
void writeGs(FILE* afs, MDINT step, int size);
void writeMB(FILE* afs, MDINT step, int size);

/* ==========================================================================*/
struct convStruct OconvStruct[]=
{
  /* Predefined converter, their code is in md2ascii.c file (Don't touch!) */
  {"float",  writeFp  }, 
  {"vector", writeVec },
  {"int",    writeInt },
  {"double", writeFp  },
  {"S",      writeS   },
  {"g",      writeg   },
  {"MB",     writeMB   },

/* ==================>>> PUT HERE YOUR CONVERTERS!!! <<< ==================
     For ex. actually it is possible to write in the command line:
     'md2ascii -t float -o measure.xmgr measure0.0'
     The previuos command means: 'use the float converter to convert the 
     measure file measure0.0 and put the result in measure.xmgr' */
  

  /* ========================================================================*/
  {NULL, NULL} /* End of list, don't touch this line */
};
#endif

/* ====================== >>> move() AND SUBROUTINES <<< ====================*/

void maps(void);

/* ==========================================================================*/
