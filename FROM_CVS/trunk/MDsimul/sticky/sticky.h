/* ================= >>>  SIMULATION DEPENDENT HEADER FILE <<< ================
   This contains file contains the follow simulation dependent data types 
   and instances:
    - instance of some structures ( params, filenames and singlePar) 
    - number macros ( NUM_MISURE, BAK_STEPS, STA_STEPS ) 
    - list macros ( that is coordinates list: SAVE_LIST, ALLOC_LIST,
                    DECL_LIST ) */

/* ================== >>> PROGRAM DEFINES(CUSTOMIZE!) <<< ===================*/
#define MD_HARDSPHERES
//#define MD_USE_SINGLE_LL
#define MDSIMUL "/home/demichel/shared/simul/mdsimul"
/* this is the executable, you must change this to your directory */
#ifdef MD_USE_CBLAS
#include <cblas.h>
#endif
#define XTERM   "/usr/X11R6/bin/nxterm"
#define MD_COORDTMP_ASCII(x) save_coordtmp_ascii(x)
void save_coordtmp_ascii(unsigned char);
extern void UpdateSystem(void);
#undef UPDATE_SYSTEM
#define UPDATE_SYSTEM UpdateSystem();
#ifdef MD_GRAVITY
#undef ADJUST_LASTCOL 
#define ADJUST_LASTCOL AdjustLastcol();
#endif
#define MD_HOME "./"
#define MD_SIMDAT MD_HOME ""
#define MD_HD_TMP MD_SIMDAT ""
/* directory to store temporary files */ 

#define MD_HD_MIS MD_SIMDAT "" 
/* directory to store measures files */

#define MD_TAPE_TMP "/iomega/mdtmp/"
/* directory on Tape to store some temporary files (restore files and 
   measures files)*/

#define MD_TAPE_MIS "/iomega/measures/"
/* directory on tape to store measures */

/* 16/1/1998 ADD: <XVA defines> */
#define MD_HD_XVA "/work2/xvafiles/"
#define MD_TAPE_XVA "/iomega/xvafiles/"


#define NUM_MISURE 30         /* maximum number of
				 measure done during simulation */
#ifdef MDLLINT
#undef  MDINT
#define MDINT long long int
#define MDINTFMT "%lld"
#else
#undef  MDINT
#define MDINT int
#define MDINTFMT "%d"
#endif

#define NDIM 3
#define MD_DEBUG(X)  
#define MD_DEBUG2(X)  
#define MD_DEBUG3(X) 
#define MD_DEBUG4(X) 
/* ========================================================================= */

/* ====================== >>> SIMULATION DEFINES <<< ========================*/
enum {MD_CORE_BARRIER=0,MD_INOUT_BARRIER,MD_OUTIN_BARRIER,MD_EVENT_NONE};

#define C_T COORD_TYPE
#define NK 10000
#define NA 5 /* number of atoms for each molecule (particle) */
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
#define MD_PBONDS_AA 4 /* questo è il max num di bonds possibili fra due molecole */
#define MD_PBONDS_BB 9
#define MD_PBONDS_AB 6
#define MD_PBONDS 9
#elif defined(MD_AB41)
#define MD_PBONDS_AA 16 /* questo è il max num di bonds possibili fra due molecole */
#define MD_PBONDS_AB 4
#define MD_PBONDS 16
#else
#define MD_PBONDS 8 /* questo è il max num di bonds possibili fra due molecole */
#endif
/* chiedere a Francesco!! */
#else
#define MD_PBONDS 8
#endif
#define MD_DIST_ELECTSITES 0.45
#define MD_DIST_HYDROSITES 0.5
#define MAXPAR 200000      /* maximum number of simulated particles */

#define NUM_PAR 2000   /* Number of particles for the simulation */
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
#ifdef MD_GRAVITY
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, vx, vy, vz, wx, wy, wz, lastcol
#else
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, vx, vy, vz, wx, wy, wz
#endif
#undef  EXT_SLST
#ifdef MD_GRAVITY
#define EXT_SLST  &L, &Lz
#else
#define EXT_SLST  &L
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
#ifdef MD_GRAVITY
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &lastcol
#else
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz 
#endif

/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#ifdef MD_GRAVITY
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *lastcol
#else
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz
#endif
				   
#undef EXT_DLST
#ifdef MD_GRAVITY
#define EXT_DLST  L, Lz 
#else
#define EXT_DLST  L
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
#define treeIdC    tree[9]
#define treeIdD    tree[10]
#define treeIdE    tree[11]
#endif
#ifdef MD_CALENDAR_HYBRID
#define treeNext   tree[12]
#define treePrev   tree[13]
#define treeQIndex tree[14]
/* serve per realizzare il timeshift */
#define treeStatus tree[15] 
//#define nlists 1000000
//#define scale 200000
#if 0
typedef struct
{
  int next;
  int previous;
  //int p;
  //int q;
  //int c;
  double t;
  //unsigned int qMarker;
  int qIndex;
  //statusType status;
  /* 06/05/2010: my fields */
  int evId;
  int evIdA;
  int evIdB;
  int idata;
  int idatb;
  int idcollcode;
} eventQEntry;
#endif
#endif
struct LastBumpS 
{
  int mol;
  int at;
  int type;
};
#ifdef MD_ROTDIFF_MIS
void calcrotMSD(void);
#endif
#define ATOM_LIMIT 10000000

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
#ifdef MD_SCALEPHI_STAGES
  int growthType;
#endif

#ifdef MD_BILOG
  double basew;
  int lastbilogsaved;
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
  COORD_TYPE DQxy;
  COORD_TYPE DQyz;
  COORD_TYPE DQzx;
  double tmsd2endA;
  double tmsd2endB; 
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
  double lastcoll;
#endif
  COORD_TYPE PxyArr[5];
  COORD_TYPE PyzArr[5];
  COORD_TYPE PzxArr[5];
#ifdef MD_DYNAMIC_OPROG
  void *ptr;
  int len;
  int (*dyn_alloc_oprog)(void);
  void (*set_dyn_ascii)(void);
#ifdef MD_ROTDIFF_MIS
  double *sumox;
  double *sumoy;
  double *sumoz;
  double *lastcolltime;
#endif
#if 0
  double *vcmx0;
  double *vcmy0;
  double *vcmz0;
#endif
  double *rxCMi; /* initial coordinates of center of mass */
  double *ryCMi; /* MAXPAR is the maximum number of particles */
  double *rzCMi;
  double **DR;
#endif

#ifdef MD_ROTDIFF_MIS
#ifndef MD_DYNAMIC_OPROG
  double sumox[MAXPAR];
  double sumoy[MAXPAR];
  double sumoz[MAXPAR];
  double lastcolltime[MAXPAR];
#endif
#endif
#ifdef MD_CALENDAR_HYBRID
  int nlistsHQ;
  double scaleHQ;
  int adjustHQ;
  double baseIndex;
  int curIndex;
  int overthrHQ;
#endif
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
#endif
  COORD_TYPE W;

  int savedXva; 
  int CMreset;
  int mdseed;
  int nebrTabFac;                /* How much storage sould be provided for 
                                   the neighbour list (see Rapaport pag.53
				   for details )*/
  COORD_TYPE rNebrShell;   /* Dr of shell of neighbour list shell see Rapaport pag. 53 */

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
  int checkGrazing;
  //int tryharder;
  double h;
  double epsd;
  double epsdFast;
  double epsdFastR;
  double epsdMax;
  int assumeOneBond;
  int forceguess;
  double targetPhi;
  double scalfact;
  double reducefact;
  double phitol;
  double minDist;
  double axestol;

#if defined(MD_GRAVITY)
#if 0
  double taptau;
  int tapampl;
  double checkquenchTime;
#endif
  double extraLz;
#endif
  double nextcheckTime;
  double nextSumTime;
  double nextDt;
#ifdef MD_DOUBLE_DT
  double nextDtR;
#endif
#ifdef MD_BIG_DT
  double refTime;
  double bigDt;
#endif
  /* questi servono per salvare le conf usando la stessa formula di Giuseppe */
  double nextStoreTime;
  int KK;
  int JJ;
  double storerate;
  int maxbonds;
  double intervalSum;
  int eventMult;
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
#ifdef MD_SAVE_REALLY_ALL
  int saveReallyAll;
  int readBinTree;
  char iniTree[NAME_LENGTH];
  char iniBak[NAME_LENGTH];
#endif
  /* ======================================================================= */
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
#ifdef MD_DOUBLE_DT
  double DtR;
#endif
  double Dt;
#ifdef MD_GRAVITY
  double ggrav;
#endif
 
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  double wallDiss;             /* dissipazione negli urti contro il muro */
  double partDiss;             /*dissipazione negli urti fra particelle */
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m[2];             /* atoms masses */
  double sigma[2][2];
#ifdef MD_AB41
  double sigmaStickyAA;
  double sigmaStickyAB;
#else
  double sigmaSticky;
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  double rcut[4];
#else
  double rcut;
#endif

  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
#ifdef MD_AB41
  double bheightAA;
  double bheightAB;
#else
  double bheight;
#endif
#ifndef MD_ASYM_ITENS
  double I[2];
#endif
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

#ifdef MD_DYNAMIC_OPROG
#define DOS(_A_) OprogStatusDyn._A_
#endif
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
#ifdef MD_SCALEPHI_STAGES
  {"growthType",      &OS(growthType),                1,       1,  "%d" },
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
  {"tmsd2endA",    &OS(tmsd2endA),                  1,              1, "%.10G"},
  {"tmsd2endB",    &OS(tmsd2endB),                  1,              1, "%.10G"},
  {"DQxy",         &OS(DQxy),                       1,              1, "%.10G"},
  {"DQyz",         &OS(DQyz),                       1,              1, "%.10G"},
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

#ifdef MD_ROTDIFF_MIS
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
  /* non serve salvarli in quanto vengono sempre inizializzati */
  {"baseIndex",    &OS(baseIndex),                  1,   1,  "%.15G"},
  {"curIndex", &OS(curIndex),               1,   1,  "%d"},
#endif
#endif
  {"hist",         OS(hist),                  MAXBIN,               1, "%d"},
  {"sumS",         OS(sumS),                    NUMK,               1, "%.6G"},
#if 1
#ifdef MD_DYNAMIC_OPROG
  /* qui vengo inizializzati con NULL poiché il puntatore verrà poi assegnato correttamente
     dal metodo OprogStatus.set_dyn_ascii() */
  {"DR",           NULL,                  -MAXPAR,              3, "%.10G"},
  {"rxCMi",        NULL,                  -MAXPAR,              1, "%.10G"},
  {"ryCMi",        NULL,                  -MAXPAR,              1, "%.10G"},
  {"rzCMi",        NULL,                  -MAXPAR,              1, "%.10G"},
#else
  {"DR",           OS(DR),                     -MAXPAR,              3, "%.10G"},
  {"rxCMi",        OS(rxCMi),                  -MAXPAR,              1, "%.10G"},
  {"ryCMi",        OS(ryCMi),                  -MAXPAR,              1, "%.10G"},
  {"rzCMi",        OS(rzCMi),                  -MAXPAR,              1, "%.10G"},
#endif
#endif
  {"histMB",       OS(histMB),                  NUMV,               1, "%d"},
  {"sumTemp",      &OS(sumTemp),                    1,              1, "%.6G"},
  {"sumPress",     &OS(sumPress),                   1,              1, "%.6G"},
  //{"sumMedia",     &OS(sumMedia),                   1,   1, "%.6f"},
  //{"sumVx",        OS(sumVx),                    MAXPAR,  1, "%.10f"},
  //{"sumVy",        OS(sumVy),                    MAXPAR,  1, "%.10f"},
  //{"sumVz",        OS(sumVz),                    MAXPAR,  1, "%.10f"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
  {"savedXva",     &OS(savedXva),                   1,   1,   "%d"},
  {"CMreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.6G"},
  {"tolT",         &OS(tolT),                       1,   1, "%.8G"},
#if defined(MD_GRAVITY)
#if 0
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
  {"checkquechTime",&OS(checkquenchTime),           1,  1,    "%.15G"},
#endif  
  {"extraLz",      &OS(extraLz),                    1,   1, "%.15G"},
#endif
  {"scalevel",     &OS(scalevel),              1,   1, "%d"},
  {"equilibrated", &OS(equilibrated),          1,   1, "%d"},
  {"assumeOneBond",&OS(assumeOneBond),         1,   1, "%d"},
  {"epsd",         &OS(epsd),                  1,   1, "%.12G"},
  {"h",            &OS(h),                     1,   1, "%.15G"}, 
  {"epsdFast",     &OS(epsdFast),              1,   1, "%.12G"},
  {"epsdFastR",    &OS(epsdFastR),             1,   1, "%.12G"},
  {"epsdMax",      &OS(epsdMax),               1,   1, "%.12G"},
  //{"guessDistOpt", &OS(guessDistOpt),          1,   1, "%d"},
  {"forceguess",   &OS(forceguess),            1,   1, "%d"},
  {"zbrakn",       &OS(zbrakn),              1,   1,  "%d"},
  {"zbrentTol",    &OS(zbrentTol),           1,   1,  "%.15G"},
  {"scalfact",     &OS(scalfact),              1,   1, "%.12G"},
  {"reducefact",   &OS(reducefact),            1,   1, "%.12G"},
  {"rescaleTime",  &OS(rescaleTime),                1,  1,    "%.10G"},
  {"phitol",       &OS(phitol),                     1,  1,    "%.14G"},
  {"axestol",      &OS(axestol),                    1,  1,    "%.14G"},
  {"minDist",      &OS(minDist),                    1,  1,    "%.14G"},
  {"endtime",      &OS(endtime),                    1,  1,    "%.15G"},
  {"nextcheckTime",&OS(nextcheckTime),              1,  1,    "%.15G"},
  {"nextSumTime"  ,&OS(nextSumTime),                1,  1,    "%.15G"},
#ifdef MD_DOUBLE_DT
  {"nextDtR",       &OS(nextDtR),                     1,  1,    "%.15G"},
#endif
  {"nextDt",       &OS(nextDt),                     1,  1,    "%.15G"},
#ifdef MD_BIG_DT
  {"refTime",      &OS(refTime),                    1,  1,    "%.15G"},
  {"bigDt",        &OS(bigDt),                      1,  1,    "%.15G"},
#endif
  {"eqlevel",     &OS(eqlevel),                    1,  1,    "%.12G"},
  {"eventMult",    &OS(eventMult),                  1,   1,  "%d"},  
  {"overlaptol"   ,&OS(overlaptol),                 1,   1, "%f"},
  {"ipart",        &OS(ipart),                      1,   1, "%d"},
  {"brownian",     &OS(brownian),                   1,   1, "%d"},
  {"checkGrazing", &OS(checkGrazing),               1,   1, "%d"},
  //{"tryharder",    &OS(tryharder),                  1,   1, "%d"},
  {"HNBOX",        &OS(HNBOX),                      1,   1,  "%d"},
  {"avngTemp",     &OS(avngTemp),                   1,   1,  "%d"},
  {"avngPress",    &OS(avngPress),                  1,   1,  "%d"},
  {"avnggr",       &OS(avnggr),                     1,  1,  "%d"},
  {"avngS",        &OS(avngS),                      1,   1,  "%d"},
  //{"avgMedia",     &OS(avgMedia),                   1,     1, "%d"},
  {"xvaSavedMode", &OS(xvaSaveMode),                1,  1,    "%d"},
  {"bakSavedMode", &OS(bakSaveMode),                1,  1,    "%d"},
  {"intervalSum"   ,&OS(intervalSum),               1,  1,    "%.10G"}, 
  {"nextStoreTime", &OS(nextStoreTime),             1, 1,     "%.15G"},
  {"storerate",     &OS(storerate),                 1, 1,     "%.10G"},
  {"KK",            &OS(KK),                        1, 1,     "%d"},
  {"JJ",            &OS(JJ),                        1, 1,     "%d"},
  {"tmpPath",      OS(tmpPath),                     1,  NAME_LENGTH, "%s"},
  {"misPath",      OS(misPath),                     1,  NAME_LENGTH, "%s"},
  {"base",         &OS(base),                       1,  1, "%.6G"},
#ifdef MD_BILOG
  {"basew",        &OS(basew),                      1,  1, "%.6G"},
  {"lastbilogsaved",&OS(lastbilogsaved),            1,  1, "%d"},
#endif
  {"NN",           &OS(NN),                         1,  1,   "%d"},
  {"fstps",        &OS(fstps),                      1,  1,   "%.15G"},
  {"nRun",         OS(nRun),                        1,  32,   "%s"},
  {"ENmin",        &OS(ENmin),                      1,  1,  "%.6G"},
  {"ENmax",        &OS(ENmax),                      1,  1,  "%.6G"},
  {"PE",           OS(PE),                 PE_POINTS,    1,  "%d"},
  {"targetPhi",   &OS(targetPhi),        1, 1, "%.15G"},
  {"maxbonds",      &OS(maxbonds),       1, 1, "%d"},
  {"eventCount",    &OS(eventCount),     1, 1, "%d"},
  {"collCount",     &OS(collCount),      1, 1, "%d"},
  {"crossCount",    &OS(crossCount),     1, 1, "%d"},
  {"scaleVolTime",  &OS(scaleVolTime),   1, 1, "%.15G"},
#ifdef MD_SAVE_REALLY_ALL
  {"saveReallyAll",  &OS(saveReallyAll),  1, 1, "%d"},
  {"readBinTree",    &OS(readBinTree),    1, 1, "%d"},
  {"iniTree",        &OS(iniTree),         1, NAME_LENGTH, "%s"},
  {"iniBak",         &OS(iniBak),          1, NAME_LENGTH, "%s"},
#endif
#ifdef MD_SURV_PROB
  {"spdeltat",          &OS(spdeltat),  1,   1, "%.15G"},
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
  {"parnumA",           &OP(parnumA),                      1,     1, "%d"},
  {"totStep",           &OP(totStep),                     1,   1,  MDINTFMT},
  {"time",              &OP(time),                        1,   1, "%.15G"},
  {"curStep",           &OP(curStep),                     1,   1,  MDINTFMT},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
  {"m",                OP(m),                            2,   1, "%.10G"},
  {"sigma",           OP(sigma),                        2,   2, "%.10G"},
#ifdef MD_AB41
  {"sigmaStickyAA",     &OP(sigmaStickyAA),                  1,   1, "%.10G"},
  {"sigmaStickyAB",     &OP(sigmaStickyAB),                  1,   1, "%.10G"},
#else
  {"sigmaSticky",     &OP(sigmaSticky),                  1,   1, "%.10G"},
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  {"rcut",              OP(rcut),                        3,   1, "%.10G"},
#else
  {"rcut",              &OP(rcut),                        1,   1, "%.10G"},
#endif
#ifdef MD_GRAVITY
  {"ggrav",             &OP(ggrav),                       1,   1,   "%f"},
#endif
  {"equilibrat",        &OP(equilibrat),                  1,   1,   "%d"},
#ifdef MD_DOUBLE_DT
  {"DtR",                &OP(DtR),                          1,   1, "%.15G"},
#endif
  {"Dt",                &OP(Dt),                          1,   1, "%.15G"},
#ifndef MD_ASYM_ITENS
  {"I",                 OP(I),                             2,   1, "%.8G"},
#endif
  {"M",                 &OP(M),                           1,   1,   "%d"},
#ifdef MD_AB41
  {"bheightAA",           &OP(bheightAA),                     1,   1, "%.15G"},
  {"bheightAB",           &OP(bheightAB),                     1,   1, "%.15G"},
#else
  {"bheight",           &OP(bheight),                     1,   1, "%.15G"},
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
int maxcoll=-1;
COORD_TYPE DECL_LIST;
COORD_TYPE EXT_DLST;
#else
extern int maxcoll;
extern COORD_TYPE DECL_LIST; 
extern COORD_TYPE EXT_DLST;
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
  {"parnumA" ,   &Oparams.parnumA,            INT},
  {"stepnum",    &Oparams.totStep,            LLINT},
  {"inifile" ,   &OprogStatus.inifile,        STR},
  {"endfile" ,   &OprogStatus.endfile,        STR},
  {"xvafile" ,   &OprogStatus.xvafile,        STR},
  {"inistep" ,   &Oparams.curStep,            LLINT},
  {"endFormat",  &OprogStatus.endFormat,      INT},
#ifdef MD_SCALEPHI_STAGES
  {"growthType",       &OprogStatus.growthType,  INT},
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
  {"tmsd2endA",  &OprogStatus.tmsd2endA,      CT},
  {"tmsd2endB",  &OprogStatus.tmsd2endB,      CT},
  {"brownian",   &OprogStatus.brownian,       INT},
  {"checkGrazing",  &OprogStatus.checkGrazing, INT},
  //{"tryharder",     &OprogStatus.tryharder,    INT},
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
#ifdef MD_ROTDIFF_MIS
  {"rotMSDSteps",  &OprogStatus.measSteps[10],   LLINT},
  {"rotMSDCalc",   &OprogStatus.measCalc[10],    LLINT},   
  {"rotMSDName",   &OprogStatus.dataFiles[10],   STR},
#endif
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"overlaptol", &OprogStatus.overlaptol,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
  {"intervalSum", &OprogStatus.intervalSum,   CT},
  {"rescaleTime", &OprogStatus.rescaleTime,   CT},
  {"storerate",     &OprogStatus.storerate,     CT},
  {"scalevel",   &OprogStatus.scalevel,       INT},
  {"endtime",    &OprogStatus.endtime,        CT},
#ifdef MD_DOUBLE_DT
  {"DtR",         &Oparams.DtR,                 CT},
#endif
  {"Dt",         &Oparams.Dt,                 CT},
#ifdef MD_BIG_DT
  {"bigDt",      &OprogStatus.bigDt,           CT},
#endif
  {"epsd",       &OprogStatus.epsd,           CT},
  {"h"   ,     &OprogStatus.h,              CT},
  {"epsdFast",   &OprogStatus.epsdFast,       CT},
  {"epsdFastR",  &OprogStatus.epsdFastR,      CT},
  {"epsdMax",    &OprogStatus.epsdMax,        CT},
  //{"guessDistOpt",&OprogStatus.guessDistOpt,  INT},
  {"zbrakn",     &OprogStatus.zbrakn,         INT},
  {"assumeOneBond", &OprogStatus.assumeOneBond, INT},
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
  {"W",          &OprogStatus.W,              CT},
  {"P",          &Oparams.P,                  CT},
  {"L",          &L,                          CT},
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
  {"sigmaAA", &Oparams.sigma[0][0], CT},
  {"sigmaBB",  &Oparams.sigma[1][1], CT},
  {"sigmaAB", &Oparams.sigma[0][1], CT},
  {"massA",   &Oparams.m[1], CT},
  {"massB",    &Oparams.m[0], CT},
  {"bheight",  &Oparams.bheight,            CT},
#ifndef MD_ASYM_ITENS
  {"IA",      &Oparams.I[1],      CT},
  {"IB",       &Oparams.I[0],      CT},
#endif
#elif defined(MD_AB41)
  {"sigmaAA", &Oparams.sigma[0][0], CT},
  {"sigmaBB",  &Oparams.sigma[1][1], CT},
  {"sigmaAB", &Oparams.sigma[0][1], CT},
  {"massA",   &Oparams.m[1], CT},
  {"massB",    &Oparams.m[0], CT},
  {"bheightAA",  &Oparams.bheightAA,  CT},
  {"bheightAB",  &Oparams.bheightAB,  CT},
#ifndef MD_ASYM_ITENS
  {"IA",      &Oparams.I[1],      CT},
  {"IB",       &Oparams.I[0],      CT},
#endif
#else
  {"sigmaSiSi", &Oparams.sigma[1][1], CT},
  {"sigmaOO",  &Oparams.sigma[0][0], CT},
  {"sigmaSiO", &Oparams.sigma[0][1], CT},
  {"massSi",   &Oparams.m[1], CT},
  {"massO",    &Oparams.m[0], CT},
  {"bheight",  &Oparams.bheight,            CT},
#ifndef MD_ASYM_ITENS
  {"ISi",      &Oparams.I[1],      CT},
  {"IO",       &Oparams.I[0],      CT},
#endif
#endif
#else
  /* sigma[0][0] = atomo grande  
   * bheight = barriere dell'interazione a buca quadrata tra H e Elettrone */
  {"sigma",    &Oparams.sigma[0][0],          CT},
  {"mass",     &Oparams.m[0], CT},
  {"bheight",  &Oparams.bheight,            CT},
  {"I",        &Oparams.I[0],      CT},
#endif
#ifdef MD_AB41
  {"sigmaStickyAA",&Oparams.sigmaStickyAA,        CT},
  {"sigmaStickyAB",&Oparams.sigmaStickyAB,        CT},
#else
  {"sigmaSticky",&Oparams.sigmaSticky,        CT},
#endif
  {"avngTemp",   &OprogStatus.avngTemp,       INT},
  {"avngPress",  &OprogStatus.avngPress,      INT},
  {"avngS",      &OprogStatus.avngS,          INT},
  {"avnggr",     &OprogStatus.avnggr,         INT},
  {"avngMB",     &OprogStatus.avngMB,         INT},
#ifdef MD_GRAVITY
  {"ggrav",      &Oparams.ggrav,            CT},
  {"Lz",         &Lz,                       CT},
  {"extraLz",    &OprogStatus.extraLz,      CT},
#endif
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
  {"targetPhi", &OprogStatus.targetPhi, CT},
  {"eventMult",  &OprogStatus.eventMult,    INT},
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
#if defined(MD_THREESPOTS) || defined(MD_AB41)
  {"rcutAA",       &Oparams.rcut[0],             CT},
  {"rcutBB",      &Oparams.rcut[1],             CT},
  {"rcutAB",      &Oparams.rcut[2],             CT},
#else
  {"rcutOO",       &Oparams.rcut[0],             CT},
  {"rcutSiSi",     &Oparams.rcut[1],             CT},
  {"rcutSiO",      &Oparams.rcut[2],             CT},
#endif
#else
  {"rcut",       &Oparams.rcut,             CT},
#endif
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"eqlevel",    &OprogStatus.eqlevel,       CT},
  {"temperat",   &Oparams.T,                CT},
  //{"tol",        &Oparams.tol,              CT},
  {"seed",       &mdseed,                   INT},
  /* parametri per scegliere il tipo di salvataggio del file xva
     (lineare o semilog) */
  {"xvaSaveMode",&OprogStatus.xvaSaveMode,  INT},
  {"bakSaveMode", &OprogStatus.bakSaveMode, INT},
  {"tmpPath",    OprogStatus.tmpPath,        STR},
  {"misPath",    OprogStatus.misPath,       STR},
  {"base",       &OprogStatus.base,         CT},
  {"NN",         &OprogStatus.NN,           INT},
#ifdef MD_BILOG
  {"basew",     &OprogStatus.basew,         CT},
#endif
  {"ENmin",      &OprogStatus.ENmin,        CT},
  {"ENmax",      &OprogStatus.ENmax,        CT},
  {"nRun",       &OprogStatus.nRun,         STR},
#ifdef MD_SAVE_REALLY_ALL
  {"saveReallyAll",  &OprogStatus.saveReallyAll,  INT},
  {"readBinTree",    &OprogStatus.readBinTree,    INT},
  {"iniTree",        &OprogStatus.iniTree,        STR},
  {"iniBak",         &OprogStatus.iniBak,         STR},
#endif
  {"maxcoll",    &maxcoll,                        INT},  
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
COORD_TYPE E, Dtrans, temp, S[NUMK], dummy, eta, gr[MAXBIN], invs, press,
  press_m, press_at, rcmz, rho, ItensD[2][3], pressST, pressHS, pressKin;
COORD_TYPE Ptens[3], DQtens[3], DphiSq, DphiSqA, DphiSqB, Ktra, Krot, sqrtdr2, Aa, V, DrSqTot, temp_transl, DrSqTotA, DrSqTotB;
int MB[NUMV];
#else 
extern COORD_TYPE E, Dtrans, temp, S[NUMK], dummy, eta, gr[MAXBIN], invs,
  press, press_m, press_at, temp_transl, rcmz, rho, pressST, pressHS, pressKin;
extern COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, V, Aa, DrSqTot, Krot, Ktra,
  DphiSq, ItensD[2][3], DphiSqA, DphiSqB, DrSqTotA, DrSqTotB;
extern int MB[NUMV];
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
void calcpress(void);

/* =========================================================================*/
#ifdef MAIN
struct measure Omeasure[NUM_MISURE]=
{
  /* DESCRIPTION:
     {<1>, <2>, <3>} where:
     <1> = pointer to the single measure, that could a variable, an array, 
     struct or whatevere you want. 
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
  {&press,  sizeof(COORD_TYPE),       calcpress   },
  {Ptens,   sizeof(COORD_TYPE)*3,     Ptensor     },
  {DQtens,  sizeof(COORD_TYPE)*3,     DQtensor    },
  {&V,      sizeof(double),          calcV       },
#ifdef MD_ROTDIFF_MIS
  {&DphiSq, sizeof(double),           calcrotMSD    },
#endif
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
