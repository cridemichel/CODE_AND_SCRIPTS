/* ================= >>>  SIMULATION DEPENDENT HEADER FILE <<< ================
   This contains file contains the follow simulation dependent data types 
   and instances:
    - instance of some structures ( params, filenames and singlePar) 
    - number macros ( NUM_MISURE, BAK_STEPS, STA_STEPS ) 
    - list macros ( that is coordinates list: SAVE_LIST, ALLOC_LIST,
                    DECL_LIST ) */

/* ================== >>> PROGRAM DEFINES(CUSTOMIZE!) <<< ===================*/
#define MD_HARDSPHERES

#define MDSIMUL "/home/demichel/shared/simul/mdsimul"
/* this is the executable, you must change this to your directory */
#ifdef MD_USE_CBLAS
#include <cblas.h>
#endif
#define XTERM   "/usr/X11R6/bin/nxterm"
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
#define NA 6 /* number of atoms for each molecule (particle) */

#define MAXPAR 5000      /* maximum number of simulated particles */
#ifdef MD_PATCHY_HE
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2
#define MD_PBONDS 10
#endif
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
#ifdef MD_ASYM_ITENS
#ifdef MD_GRAVITY
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, Mx, My, Mz, lastcol
#else
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, Mx, My, Mz
#endif
#else
#ifdef MD_GRAVITY
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, vx, vy, vz, wx, wy, wz, lastcol
#else
#define SAVE_LIST rx, ry, rz, vx, vy, vz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, vx, vy, vz, wx, wy, wz
#endif
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
#ifdef MD_ASYM_ITENS
#ifdef MD_GRAVITY
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &Mx, &My, &Mz, &lastcol
#else
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &Mx, &My, &Mz 
#endif
#else
#ifdef MD_GRAVITY
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz, &lastcol
#else
#define ALLOC_LIST  &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &uzx, &uzy, &uzz, &vx, &vy, &vz, &wx, &wy, &wz 
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
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *Mx, *My, *Mz, *lastcol
#else
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *Mx, *My, *Mz
#endif
#else
#ifdef MD_GRAVITY
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz, *lastcol
#else
#define DECL_LIST   *rx, *ry, *rz, *uxx, *uxy, *uxz, *uyx, *uyy, *uyz, *uzx, *uzy, *uzz, *vx, *vy, *vz, *wx, *wy, *wz
#endif
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
#ifdef MD_PATCHY_HE
#define treeIdC    tree[9]
#define treeIdD    tree[10]
#define treeIdE    tree[11]
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
  double sumox[MAXPAR];
  double sumoy[MAXPAR];
  double sumoz[MAXPAR];
  double lastcolltime[MAXPAR];
  double springkSD;
  int SDmethod;
  double toldxNR;
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
  
  COORD_TYPE vcmx0[MAXPAR];
  COORD_TYPE vcmy0[MAXPAR];
  COORD_TYPE vcmz0[MAXPAR];
  
  COORD_TYPE rxCMi[MAXPAR]; /* initial coordinates of center of mass */
  COORD_TYPE ryCMi[MAXPAR]; /* MAXPAR is the maximum number of particles */
  COORD_TYPE rzCMi[MAXPAR];
  COORD_TYPE DR[MAXPAR][3];
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
  double epsd;
  double epsdFast;
  double epsdFastR;
  double epsdMax;
  double epsdNL;
  double epsdFastNL;
  double epsdFastRNL;
  double epsdMaxNL;
#ifdef MD_PATCHY_HE
  double epsdSP;
  double epsdFastSP;
  double epsdSPNL;
  double epsdFastSPNL;
#endif
  int guessDistOpt;
  int forceguess;
  double targetPhi;
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

  /* ======================================================================= */
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
  
  double a[2];
  double b[2];
  double c[2];
  double rcut;
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   

#ifndef MD_ASYM_ITENS
  double I[2];
#else
  double I[2][3];
#endif
#ifdef MD_PATCHY_HE
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
#ifndef MAIN
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
  {"sumox",        OS(sumox),                       -MAXPAR,        1, "%.15G"},
  {"sumoy",        OS(sumoy),                       -MAXPAR,        1, "%.15G"},
  {"sumoz",        OS(sumoz),                       -MAXPAR,        1, "%.15G"},
  {"rxCMi",        OS(rxCMi),                       -MAXPAR,        1, "%.15G"},
  {"ryCMi",        OS(ryCMi),                       -MAXPAR,        1, "%.15G"},
  {"rzCMi",        OS(rzCMi),                       -MAXPAR,        1, "%.15G"},
  {"DR",           OS(DR),                          -MAXPAR,        3, "%.15G"}, 
  {"hist",         OS(hist),                  MAXBIN,               1, "%d"},
  {"sumS",         OS(sumS),                    NUMK,               1, "%.6G"},
  {"histMB",       OS(histMB),                  NUMV,               1, "%d"},
  {"sumTemp",      &OS(sumTemp),                    1,              1, "%.6G"},
  {"sumPress",     &OS(sumPress),                   1,              1, "%.6G"},
  //  {"sumMedia",     &OS(sumMedia),                   1,   1, "%.6f"},
  //{"sumVx",        OS(sumVx),                    MAXPAR,  1, "%.10f"},
  //{"sumVy",        OS(sumVy),                    MAXPAR,  1, "%.10f"},
  //{"sumVz",        OS(sumVz),                    MAXPAR,  1, "%.10f"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
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
  {"epsdSP",         &OS(epsdSP),                  1,   1, "%.12G"},
  {"epsdFastSP",     &OS(epsdFastSP),              1,   1, "%.12G"},
  {"epsdSPNL",       &OS(epsdSPNL),                1,   1, "%.12G"},
  {"epsdFastSPNL",   &OS(epsdFastSPNL),            1,   1, "%.12G"},
  {"guessDistOpt", &OS(guessDistOpt),          1,   1, "%d"},
  {"springkSD",    &OS(springkSD),              1,   1, "%.12G"},
  {"SDmethod",     &OS(SDmethod),               1,   1, "%d"},
  {"stepSDA",       &OS(stepSDA),                1,   1, "%.12G"},
  {"toldxNR",       &OS(toldxNR),                1,   1, "%.15G"}, 
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
  {"zbrentTol",    &OS(zbrentTol),           1,   1,  ".15G"},
  {"scalfact",     &OS(scalfact),              1,   1, ".12G"},
  {"reducefact",   &OS(reducefact),            1,   1, ".12G"},
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
  {"base",         &OS(base),                       1,  1, "%.6G"},
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
  {"m",                 OP(m),                            2,   1, "%.6G"},
  {"a",                 OP(a),                             2,   1, "%.8G"},
  {"b",                 OP(b),                             2,   1, "%.8G"},
  {"c",                 OP(c),                             2,   1, "%.8G"},
  {"rcut",              &OP(rcut),                        1,   1, "%.10G"},
  {"equilibrat",        &OP(equilibrat),                  1,   1,   "%d"},
  {"Dt",                &OP(Dt),                          1,   1, "%.15G"},
#ifndef MD_ASYM_ITENS
  {"I",                 OP(I),                             2,   1, "%.8G"},
#else
  {"I",                 OP(I),                             2,   3, "%.8G"},
#endif
#ifdef MD_GRAVITY
  {"wallDiss",          &OP(wallDiss),                    1,   1,   "%f"},
  {"partDiss",          &OP(partDiss),                    1,   1,   "%f"},
  {"ggrav",             &OP(ggrav),                       1,   1,   "%f"},
#endif
  {"M",                 &OP(M),                           1,   1,   "%d"},
  {"tol",               &OP(tol),                         1,   1, "%.15G"},
#ifdef MD_PATCHY_HE
  {"sigmaSticky",       &OP(sigmaSticky),                       1,   1, "%.15G"},
  {"bheight",           &OP(bheight),                     1,   1, "%.15G"},
  {"bhin",               &OP(bhin),                         1,   1, "%.15G"},
  {"bhout",              &OP(bhout),                         1,   1, "%.15G"},
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
#else
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
  {"epsdSP",       &OprogStatus.epsdSP,           CT},
  {"epsdFastSP",   &OprogStatus.epsdFastSP,       CT},
  {"epsdSPNL",     &OprogStatus.epsdSPNL,         CT},
  {"epsdFastSPNL", &OprogStatus.epsdFastSPNL,     CT},
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
  {"tolAngNR",   &OprogStatus.tolAngNR, CT},
  {"stepSDA",     &OprogStatus.stepSDA,         CT},
  {"stepSDB",     &OprogStatus.stepSDB,         CT},
  {"maxitsSD",   &OprogStatus.maxitsSD,       INT},
  {"zbrakn",     &OprogStatus.zbrakn,         INT},
  {"zbrentTol",  &OprogStatus.zbrentTol,      CT},
  {"scalfact",   &OprogStatus.scalfact,       CT},
  {"reducefact", &OprogStatus.reducefact,     CT},
  {"phitol",      &OprogStatus.phitol,        CT},
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
  {"L",          &L,                        CT},
#ifdef MD_PATCHY_HE
  {"sigmaSticky", &Oparams.sigmaSticky,     CT},
  {"bheight",    &Oparams.bheight,          CT},
  {"bhin",         &Oparams.bhin,              CT},
  {"bhout",       &Oparams.bhout,             CT},
  {"assumeOneBond", &OprogStatus.assumeOneBond, INT},
  {"checkGrazing",  &OprogStatus.checkGrazing, INT},
  {"maxbonds",      &OprogStatus.maxbonds,     INT},
#endif
  {"avngTemp",   &OprogStatus.avngTemp,       INT},
  {"avngPress",  &OprogStatus.avngPress,      INT},
  {"avngS",      &OprogStatus.avngS,          INT},
  {"avnggr",     &OprogStatus.avnggr,         INT},
  {"avngMB",     &OprogStatus.avngMB,         INT},
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
  {"A0",      &Oparams.a[0],      CT},
  {"A1",      &Oparams.a[1],      CT},
  {"B0",      &Oparams.b[0],      CT},
  {"B1",      &Oparams.b[1],      CT},
  {"C0",      &Oparams.c[0],      CT},
  {"C1",      &Oparams.c[1],      CT},
  {"targetPhi", &OprogStatus.targetPhi, CT},
#ifndef MD_ASYM_ITENS
  {"Ia",      &Oparams.I[0],      CT},
  {"Ib",      &Oparams.I[1],      CT},
#else
  {"I1a",      &Oparams.I[0][0],      CT},
  {"I3a",      &Oparams.I[0][2],      CT},
  {"I1b",      &Oparams.I[1][0],      CT},
  {"I3b",      &Oparams.I[1][2],      CT},
#endif
#ifdef MD_GRAVITY
  {"ggrav",      &Oparams.ggrav,            CT},
#endif
  {"mass0",       &Oparams.m[0],                CT},
  {"mass1",       &Oparams.m[1],                CT},
#ifdef MD_GRAVITY
  {"wallDiss",   &Oparams.wallDiss,         CT},
  {"partDiss",   &Oparams.partDiss,         CT},
  {"quenchend",  &OprogStatus.quenchend,    CT},
  {"maxquench",  &OprogStatus.maxquench,    INT},
  {"tc",         &OprogStatus.tc,           CT},
#endif
  {"eventMult",  &OprogStatus.eventMult,    INT},
  {"rcut",       &Oparams.rcut,             CT},
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"eqlevel",    &OprogStatus.eqlevel,       CT},
  {"temperat",   &Oparams.T,                CT},
  {"tol",        &Oparams.tol,              CT},
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
COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, Aa, V, DrSqTot, temp_transl, DphiSq;
int MB[NUMV];
#else 
extern COORD_TYPE E, Dtrans, temp, S[NUMK], dummy, eta, gr[MAXBIN], invs,
  press, press_m, press_at, temp_transl, rcmz, rho;
extern COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, V, Aa, DrSqTot,
  DphiSq, ItensD[2][3], DphiSq, pressST, pressHS, pressKin;
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
