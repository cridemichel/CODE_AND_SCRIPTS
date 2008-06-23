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

#define XTERM   "/usr/X11R6/bin/nxterm"
#undef UPDATE_SYSTEM
#define UPDATE_SYSTEM UpdateSystem();
#undef ADJUST_LASTCOL 
#define ADJUST_LASTCOL AdjustLastcol();

#define MD_HOME "./"
#ifdef MD_LOADMESH
#define MD_MESHDIR MD_HOME "/simdat"
#endif

#define MD_HD_TMP MD_HOME ""
/* directory to store temporary files */ 

#define MD_HD_MIS MD_HOME "" 
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

#define NDIM 3
#define MD_DEBUG(X)  
#define MD_DEBUG2(X)    
#define MD_DEBUG3(X) 
#define MD_DEBUG4(X) 
/* ========================================================================= */

/* ====================== >>> SIMULATION DEFINES <<< ========================*/

#define C_T COORD_TYPE
#define NK 10000
#define NA 1 /* number of atoms for each molecule (particle) */

#define MAXPAR 10000      /* maximum number of simulated particles */

#define NUM_PAR 500   /* Number of particles for the simulation */
#define NUMK 99    /* number of k-points in which we must  calculate the 
		       structure factor */ 
#define MAXBIN 300  /* Number of radius in which calculate the radial 
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
#if defined(MD_POLYDISP)
#define SAVE_LIST rx, ry, rz, vx, vy, vz, radii, lastcol
#elif defined(MD_FULL_LANG) || defined(MD_MICRO_LANG)
#define SAVE_LIST rx, ry, rz, vx, vy, vz, lastcol
#else
#define SAVE_LIST rx, ry, rz, vx, vy, vz, lastcol
#endif
#undef  EXT_SLST
#define EXT_SLST  &L, &Lz

/* Reduced list of variables to save on xva file (tape file) */ 
#define XVA_LIST rx, ry, rz, vx, vy, vz, lastcol

#define XVA_NUM 6 /* this is the number of vars in XVA_LIST <--- SET THIS!!!*/

#define XVA_ALST &rx, &ry, &rz, &vx, &vy, &vz, &lastcol

#define XVA_DLST *rx, *ry, *rz, *vx, *vy, *vz, *lastcol
 
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
#if defined(MD_POLYDISP)
#define ALLOC_LIST  &rx, &ry, &rz, &vx, &vy, &vz, &radii, &lastcol
#elif defined(MD_FULL_LANG) || defined (MD_MICRO_LANG)
#define ALLOC_LIST  &rx, &ry, &rz, &vx, &vy, &vz, &lastcol
#else
#define ALLOC_LIST  &rx, &ry, &rz, &vx, &vy, &vz, &lastcol
#endif
/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#if defined(MD_POLYDISP)
#define DECL_LIST   *rx, *ry, *rz, *vx, *vy, *vz, *radii, *lastcol
#elif defined(MD_FULL_LANG) || defined(MD_MICRO_LANG)
#define DECL_LIST   *rx, *ry, *rz, *vx, *vy, *vz, *lastcol
#else
#define DECL_LIST   *rx, *ry, *rz, *vx, *vy, *vz, *lastcol
#endif
				   
#undef EXT_DLST
#define EXT_DLST  L, Lz 

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

  int xvaSteps;     /* steps between two tape file savings */
  int bakSteps;    /* steps between two savings of restore files on HD*/

  int bakStepsAscii; 
  int staSteps;     /* steps after which must save sim_stat structure 
				into the STATUS_FILE */
  int measSteps[NUM_MISURE];/*steps after which save every measure */
  
  int measCalc[NUM_MISURE]; /*steps between two measure calculation */
				        
  int initCalc[NUM_MISURE];
  /* first step to begin to calculate a certain measure, by default 
     this quantity is 1, but it is possible to change this value in the 
     parameter file */

  int initStep[NUM_MISURE]; 
  /* first step to begin to save a certain measure, by default this quantity 
     is 1, but it is possible to change this value in the parameter file*/

  /* =============== >>> PUT HERE YOUR STATUS FILEDS <<< =================== 
     For example accumalators (see Allen - Tildesley)*/

  COORD_TYPE sumEta; /* accumulators for obtaining the mean value of eta */
  
  /* Accumulators for the integral of angular velocity */
  COORD_TYPE DQxy;
  COORD_TYPE DQyz;
  COORD_TYPE DQzx;
#ifdef MD_HSVISCO
  double DQW;
  double DQTxy;
  double DQTyz;
  double DQTzx;
  double DQWxy;
  double DQWyz;
  double DQWzx;
  double Txy;
  double Tyz;
  double Tzx;
  double lastcoll;
#endif

  COORD_TYPE PxyArr[5];
  COORD_TYPE PyzArr[5];
  COORD_TYPE PzxArr[5];
  
  /* Accumulator for the radial distribution function */
  int hist[MAXBIN];
  
  /* Accumulator for the static structure factor */ 
  COORD_TYPE sumS[NUMK];

  /* Accumulator for the MAXWELL-BOLTZMANN distribution */ 
  int histMB[NUMV];

  COORD_TYPE sumTemp;
  COORD_TYPE sumPress;
  
  COORD_TYPE vcmx0[MAXPAR];
  COORD_TYPE vcmy0[MAXPAR];
  COORD_TYPE vcmz0[MAXPAR];
  //double radii[MAXPAR];
  /*double lastcol[MAXPAR];
    double atomTime[MAXPAR];*/

  COORD_TYPE rxCMi[MAXPAR]; /* initial coordinates of center of mass */
  COORD_TYPE ryCMi[MAXPAR]; /* MAXPAR is the maximum number of particles */
  COORD_TYPE rzCMi[MAXPAR];
  COORD_TYPE DR[MAXPAR][3];

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
  double taptau;
  int tapampl;
  int scalevel;
  int brownian;
  double targetPhi;
  double polydisp;
  double polycutoff;
  double phitol;
  double axestol;
  double scalfact;
  double checkquenchTime;
  double rescaleTime;
  double nextcheckTime;
  double nextSumTime;
  double nextDt;
  double nextStoreTime;
  int KK;
  int JJ;
  double storerate;
  int numquench;
  int maxquench;
  double rhobh;
  double intervalSum;
  double quenchtol;
  int eventMult;
  double extraLz;
  double vztap;
  double rzup;
  double expandFact;
  double quenchend;
  /* ADD 13/4/2000 
     Logarithmic saving implemented for xva file */
  int xvaSaveMode;/* 0 = linear 1 = semilog 2 = bilog (not impl. yet) */
  double accV;
  double overlaptol;
  double tc;
  double lastV;
  /* AGGIUNTI 24/04/01 */
  int bakSaveMode;/* save mode for ascii backup */
  char tmpPath[NAME_LENGTH];
  char misPath[NAME_LENGTH];
  /* =================================================== */
#ifdef MD_BILOG
  double logblockbase;
  double logblock;
#endif
  double base;    /* We save at base^^NN step */
  int NN;         /* Logatithmic block length */
  double fstps;         /* There are KK block each base^NN long */

  double accrcmz;
  double time;
  int eventCount;
  long long int collCount;
  long long int crossCount;
  int wallcollCount;
  int PE[PE_POINTS];
  double ENmin;
  double ENmax;
#ifdef MD_POLYDISP
  int ext_radii;
  char radiiFile[NAME_LENGTH];
#endif
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
  int parnum;        	/* particles number */
  int totStep;	/* temporal step number that simulation 
				   must do */
  double Dt;
  int curStep;	/* current step of simulation */
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  double wallDiss;             /* dissipazione negli urti contro il muro */
  double partDiss;             /*dissipazione negli urti fra particelle */
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m;             /* atoms masses */
#if defined(MD_FPBROWNIAN) || defined(MD_FULL_LANG) || defined(MD_MICRO_LANG)
  COORD_TYPE xi;            /* Fokker-Planck damping xi for Brownian dyn. */
#endif
  COORD_TYPE sigma;     /* pair potential length parameters */
  double rcut;
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   

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
  {"dataFiles",    OS(dataFiles),         NUM_MISURE,     NAME_LENGTH, "%s"},
  {"xvaSteps",     &OS(xvaSteps),                  1,               1, "%d"},
  {"bakSteps",     &OS(bakSteps),                  1,               1, "%d"},
  {"bakStepsAscii",&OS(bakStepsAscii),             1,               1, "%d"},
  {"staSteps",     &OS(staSteps),                  1,               1, "%d"},
  {"mdseed",       &OS(mdseed),                    1,               1, "%d"},
  {"measSteps",    OS(measSteps),         NUM_MISURE,               1, "%d"},
  {"measCalc",     OS(measCalc),          NUM_MISURE,               1, "%d"},
  {"initCalc",     OS(initCalc),          NUM_MISURE,               1, "%d"},
  {"initStep",     OS(initStep),          NUM_MISURE,               1, "%d"},
  {"sumEta",       &OS(sumEta),                     1,              1, "%.10G"},
  {"DQxy",         &OS(DQxy),                       1,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                       1,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                       1,              1, "%.10G"},
  {"PxyArr",       OS(PxyArr),                      5,              1, "%.10G"},
  {"PyzArr",       OS(PyzArr),                      5,              1, "%.10G"},
  {"PzxArr",       OS(PzxArr),                      5,              1, "%.10G"},
  {"rxCMi",        OS(rxCMi),                       -MAXPAR,        1, "%.15G"},
  {"ryCMi",        OS(ryCMi),                       -MAXPAR,        1, "%.15G"},
  {"rzCMi",        OS(rzCMi),                       -MAXPAR,        1, "%.15G"},
  {"DR",           OS(DR),                          -MAXPAR,        3, "%.15G"}, 
  {"hist",         OS(hist),                  MAXBIN,               1, "%d"},
  {"time",         &OS(time),                        1,              1, "%.15G"},
  {"tc",           &OS(tc),                          1,              1, "%.15G"},
  {"sumS",         OS(sumS),                    NUMK,               1, "%.6G"},
  {"histMB",       OS(histMB),                  NUMV,               1, "%d"},
  {"sumTemp",      &OS(sumTemp),                    1,              1, "%.6G"},
  {"sumPress",     &OS(sumPress),                   1,              1, "%.6G"},
  //  {"sumMedia",     &OS(sumMedia),                   1,   1, "%.6f"},
  //{"sumVx",        OS(sumVx),                    MAXPAR,  1, "%.10f"},
  //{"sumVy",        OS(sumVy),                    MAXPAR,  1, "%.10f"},
  //{"sumVz",        OS(sumVz),                    MAXPAR,  1, "%.10f"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
  //{"radii",        OS(radii),                      -MAXPAR,    1, "%.15G"},
  {"targetPhi",    &OS(targetPhi),                 1,          1, "%.12G"},
  {"polydisp",     &OS(polydisp),                  1,          1, "%.15G"}, 
  {"polycutoff",   &OS(polycutoff),                1,          1, "%.8G"},
  {"savedXva",     &OS(savedXva),                   1,   1,   "%d"},
  {"Cmreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.6G"},
  {"tolT",         &OS(tolT),                       1,   1, "%.8G"},
  {"quenchtol",    &OS(quenchtol),                  1,   1, "%.10G"},
  {"rhobh",        &OS(rhobh),                  1,   1, "%.10G"},
  {"quenchend",    &OS(quenchend),                  1,   1, "%f"},
  {"taptau",       &OS(taptau),                1,   1, "%f"},
  {"vztap",        &OS(vztap),                 1,   1, "%.15G"},
  {"rzup",         &OS(rzup),                  1,   1, "%.15G"},
  {"expandFact",   &OS(expandFact),            1,   1, "%.15G"},
  {"numquench",    &OS(numquench),             1,   1, "%d"},
  {"maxquench",    &OS(maxquench),             1,   1, "%d"},
  {"scalevel",     &OS(scalevel),              1,   1, "%d"},
  {"phitol",       &OS(phitol),                     1,  1,    "%.14G"},
  {"scalfact",     &OS(scalfact),             1, 1,           "%.14G"},
  {"axestol",      &OS(axestol),                    1,  1,    "%.14G"},
  {"extraLz",      &OS(extraLz),                    1,   1, "%.15G"},
  {"eventMult",    &OS(eventMult),                  1,   1,  "%d"},  
  {"overlaptol"   ,&OS(overlaptol),                 1,   1, "%f"},
  {"ipart",        &OS(ipart),                      1,   1, "%d"},
  {"HNBOX",        &OS(HNBOX),                      1,   1,  "%d"},
  {"avngTemp",     &OS(avngTemp),                   1,   1,  "%d"},
  {"avngPress",    &OS(avngPress),                  1,   1,  "%d"},
  {"avnggr",       &OS(avnggr),                     1,  1,  "%d"},
  {"avngS",        &OS(avngS),                      1,   1,  "%d"},
  //{"avgMedia",     &OS(avgMedia),                   1,     1, "%d"},
  {"xvaSavedMode", &OS(xvaSaveMode),                1,  1,    "%d"},
  {"bakSavedMode", &OS(bakSaveMode),                1,  1,    "%d"},
  {"rescaleTime",  &OS(rescaleTime),                1,  1,    "%.10G"},
  {"checkquenchTime",&OS(checkquenchTime),           1,  1,    "%.10G"},
  {"intervalSum"   ,&OS(intervalSum),               1,  1,    "%.10G"}, 
  {"nextStoreTime", &OS(nextStoreTime),             1, 1,     "%.10G"},
  {"storerate",     &OS(storerate),                 1, 1,     "%.10G"},
  {"KK",            &OS(KK),                        1, 1,     "%d"},
  {"JJ",            &OS(JJ),                        1, 1,     "%d"},
  {"nextcheckTime",&OS(nextcheckTime),              1,  1,    "%.15G"},
  {"nextSumTime"  ,&OS(nextSumTime),                1,  1,    "%.15G"},
  {"nextDt",       &OS(nextDt),                     1,  1,    "%.15G"},
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
#ifdef MD_POLYDISP
  {"ext_radii",           &OS(ext_radii),          1, 1, "%d"},
  {"radiiFile",           OS(radiiFile),          1, NAME_LENGTH, "%s"},
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
  {"totStep",           &OP(totStep),                     1,   1,  "%d"},
  {"curStep",           &OP(curStep),                     1,   1,  "%d"},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
  {"m",                 &OP(m),                           1,   1, "%.6G"},
#if defined(MD_FPBROWNIAN) || defined(MD_FULL_LANG) || defined (MD_MICRO_LANG)
  {"xi",                &OP(xi),                          1,   1, "%.6G"},
#endif
  {"sigma",             &OP(sigma),                       1,   1, "%.8G"},
  {"rcut",              &OP(rcut),                        1,   1, "%.10G"},
  {"equilibrat",        &OP(equilibrat),                  1,   1,   "%d"},
  {"wallDiss",          &OP(wallDiss),                    1,   1,   "%f"},
  {"partDiss",          &OP(partDiss),                    1,   1,   "%f"},
#ifdef MD_GRAVITY
  {"ggrav",             &OP(ggrav),                       1,   1,   "%f"},
#endif
  {"M",                 &OP(M),                           1,   1,   "%d"},
  {"tol",               &OP(tol),                         1,   1, "%.15G"},
  {"Dt",                &OP(Dt),                          1,   1, "%.15G"},
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
  {"stepnum",    &Oparams.totStep,            INT},
  {"inifile" ,   &OprogStatus.inifile,        STR},
  {"endfile" ,   &OprogStatus.endfile,        STR},
  {"xvafile" ,   &OprogStatus.xvafile,        STR},
  {"inistep" ,   &Oparams.curStep,            INT},
  {"endFormat",  &OprogStatus.endFormat,      INT},
  {"iniFormat",  &OprogStatus.iniFormat,      INT},
  {"tapeTimes",  &OprogStatus.tapeTimes,      INT},
  {"bakSteps",   &OprogStatus.bakSteps,       INT},
  {"bakStepsAscii", &OprogStatus.bakStepsAscii,       INT},
  {"staSteps",   &OprogStatus.staSteps,       INT},
  {"xvaSteps",   &OprogStatus.xvaSteps,       INT},  
  {"energySteps",&OprogStatus.measSteps[0],   INT},
  {"energyCalc", &OprogStatus.measCalc[0],    INT},
  {"energyName", &OprogStatus.dataFiles[0],   STR},
  {"energyBegin",&OprogStatus.initStep[0],   INT},
  {"DtrSteps",   &OprogStatus.measSteps[1],   INT}, /* steps between measure
						       savings */
  {"DtrCalc",    &OprogStatus.measCalc[1],    INT}, /* steps between measure
						       calculation */
  {"DtrName",    &OprogStatus.dataFiles[1],   STR},
  {"brownian",   &OprogStatus.brownian,       INT},
  {"tempSteps",  &OprogStatus.measSteps[2],   INT},
  {"tempCalc",   &OprogStatus.measCalc[2],    INT},
  {"tempName",   &OprogStatus.dataFiles[2],   STR},
  {"grSteps",    &OprogStatus.measSteps[3],   INT},
  {"grCalc",     &OprogStatus.measCalc[3],    INT},
  {"grName",     &OprogStatus.dataFiles[3],   STR},
  {"MBSteps",    &OprogStatus.measSteps[4],   INT},
  {"MBCalc",     &OprogStatus.measCalc[4],    INT},
  {"MBName",     &OprogStatus.dataFiles[4],   STR},
  {"SSteps",     &OprogStatus.measSteps[5],   INT},
  {"SCalc",      &OprogStatus.measCalc[5],    INT},
  {"SName",      &OprogStatus.dataFiles[5],   STR},
  {"pressSteps", &OprogStatus.measSteps[6],   INT},
  {"pressCalc",  &OprogStatus.measCalc[6],    INT},   
  {"pressName",  &OprogStatus.dataFiles[6],   STR},
  {"PtensSteps", &OprogStatus.measSteps[7],  INT},
  {"PtensCalc",  &OprogStatus.measCalc[7],   INT},   
  {"PtensName",  &OprogStatus.dataFiles[7],  STR},
  {"DQtensSteps",&OprogStatus.measSteps[8],  INT},
  {"DQtensCalc", &OprogStatus.measCalc[8],   INT},
  {"DQtensName", &OprogStatus.dataFiles[8],  STR},
  {"VSteps",  &OprogStatus.measSteps[9],   INT},
  {"VCalc",   &OprogStatus.measCalc[9],    INT},   
  {"VName",   &OprogStatus.dataFiles[9],   STR},
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"overlaptol", &OprogStatus.overlaptol,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
  {"taptau",     &OprogStatus.taptau,         CT},
  {"quenchtol",  &OprogStatus.quenchtol,      CT},
  {"intervalSum", &OprogStatus.intervalSum,   CT},
  {"storerate",     &OprogStatus.storerate,     CT},
  {"checkQuench", &OprogStatus.checkquenchTime, CT},
  {"rescaleTime", &OprogStatus.rescaleTime,   CT},
  {"scalevel",   &OprogStatus.scalevel,       INT},
  {"Dt",         &Oparams.Dt,                 CT},
  {"W",          &OprogStatus.W,              CT},
  {"P",          &Oparams.P,                  CT},
  {"L",          &L,                        CT},
  {"Lz",         &Lz,                       CT},
  {"extraLz",    &OprogStatus.extraLz,      CT},
  {"vztap",      &OprogStatus.vztap,        CT},
  {"rzup",       &OprogStatus.rzup,         CT},
  {"rhobh",      &OprogStatus.rhobh,        CT},
  {"expandFact",  &OprogStatus.expandFact,   CT},
  {"avngTemp",   &OprogStatus.avngTemp,       INT},
  {"avngPress",  &OprogStatus.avngPress,      INT},
  {"avngS",      &OprogStatus.avngS,          INT},
  {"avnggr",     &OprogStatus.avnggr,         INT},
  {"avngMB",     &OprogStatus.avngMB,         INT},
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
  {"sigma",      &Oparams.sigma,            CT},
#ifdef MD_GRAVITY
  {"ggrav",      &Oparams.ggrav,            CT},
#endif
  {"mass",       &Oparams.m,                CT},
#if defined(MD_FPBROWNIAN) || defined(MD_FULL_LANG) || defined(MD_MICRO_LANG)
  {"xi",         &Oparams.xi,               CT},
#endif
  {"targetPhi", &OprogStatus.targetPhi, CT},
  {"polydisp",  &OprogStatus.polydisp, CT},  
  {"polycutoff",&OprogStatus.polycutoff, CT},
  {"phitol",      &OprogStatus.phitol,        CT},
  {"scalfact",    &OprogStatus.scalfact,      CT},
  {"axestol",     &OprogStatus.axestol,       CT},
  {"wallDiss",   &Oparams.wallDiss,         CT},
  {"partDiss",   &Oparams.partDiss,         CT},
  {"quenchend",  &OprogStatus.quenchend,    CT},
  {"maxquench",  &OprogStatus.maxquench,    INT},
  {"tc",         &OprogStatus.tc,           CT},
  {"eventMult",  &OprogStatus.eventMult,    INT},
  {"rcut",       &Oparams.rcut,             CT},
  {"equilibrat", &Oparams.equilibrat,       INT},
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
#ifdef MD_BILOG
  {"basew",     &OprogStatus.basew,         CT},
#endif
  {"NN",         &OprogStatus.NN,           INT},
  {"ENmin",      &OprogStatus.ENmin,        CT},
  {"ENmax",      &OprogStatus.ENmax,        CT},
  {"nRun",       &OprogStatus.nRun,         STR},
#ifdef MD_POLYDISP
  {"ext_radii",  &OprogStatus.ext_radii,    INT},
  {"radiiFile",  OprogStatus.radiiFile,    STR},
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
  press_m, press_at, rcmz, rho;
COORD_TYPE Ptens[3], DQtens[3], 
  sqrtdr2, Aa, V, DrSqTot, temp_transl;
int MB[NUMV];
#else 
extern COORD_TYPE E, Dtrans, temp, S[NUMK], dummy, eta, gr[MAXBIN], invs,
  press, press_m, press_at, temp_transl, rcmz, rho;
extern COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, V, Aa, DrSqTot,
  DphiSq;
extern int MB[NUMV];
#endif

/* ============= >>> PUT HERE MEASURING FUNCTION PROTOTYPES <<< ============*/

void energy(void);
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
  {&press,  sizeof(COORD_TYPE),       NULL        },
  {Ptens,   sizeof(COORD_TYPE)*3,     Ptensor     },
  {DQtens,  sizeof(COORD_TYPE)*3,     DQtensor    },
  {&V,       sizeof(double),           NULL       },
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
  int saveSteps;
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
  int saveSteps; /* save between two tape savings */
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
void writeVec(FILE* afs, int step, int size);/*Predefined*/
void writeFp(FILE* afs, int step, int size); /*Predefined*/
void writeInt(FILE* afs, int step, int size);/*Predefined*/
void writeS(FILE*afs, int step, int size);
void writeg(FILE* afs, int step, int size);
void writeGs(FILE* afs, int step, int size);
void writeMB(FILE* afs, int step, int size);

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
