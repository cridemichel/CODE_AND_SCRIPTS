/* ================= >>>  SIMULATION DEPENDENT HEADER FILE <<< ================
   This contains file contains the follow simulation dependent data types 
   and instances:
    - instance of some structures ( params, filenames and singlePar) 
    - number macros ( NUM_MISURE, BAK_STEPS, STA_STEPS ) 
    - list macros ( that is coordinates list: SAVE_LIST, ALLOC_LIST,
                    DECL_LIST ) */

/* ================== >>> PROGRAM DEFINES(CUSTOMIZE!) <<< ===================*/

#define MDSIMUL "/home/demichel/shared/simul/mdsimul"
/* this is the executable, you must change this to your directory */

#define XTERM   "/usr/X11R6/bin/nxterm"
#define MD_HOME "/Users/demichel"
/* directory to store temporary files */ 
#define MD_SIMDAT MD_HOME "/simdat"

#ifdef MD_LOADMESH
#define MD_MESHDIR MD_HOME "/simdat"
#endif

#define MD_HD_TMP MD_SIMDAT "/mdtmp/"
/* directory to store temporary files */ 

#define MD_HD_MIS MD_SIMDAT "/measures/" 
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

/* ========================================================================= */

/* ====================== >>> SIMULATION DEFINES <<< ========================*/

#define NA 20 /* number of atoms for each molecule (particle) */

#define MAXPAR 15000      /* maximum number of simulated particles */

#define NUM_PAR 500   /* Number of particles for the simulation */
#define NUMK 100    /* number of k-points in which we must  calculate the 
		       structure factor */ 
#define MAXBIN 150  /* Number of radius in which calculate the radial 
		       distribution function */
#define NUMV 150
#define GSRMAX 2.0 /* Maximum r for which we calculate the van Hove function*/
#define GSPOINT 30  /* Pont for which we plot the van Hove function */
#define GSANGPOINT 30 /* Point for which we plot the Angular Van Hove 
			 function */
#define FSPOINT 30    /* Point for which we plot the self-part of the 
			 intermediate scattering function */
#define FSKMAX 15.0 /* see the maximum k for the static structure factor */
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
#define MD_ALLOC_POLY
#ifdef MD_RESPA
#define SAVE_LIST rx, NA, ry, NA, rz, NA,\
                  vx, NA, vy, NA, vz, NA,\
                  Fx, NA, Fy, NA, Fz, NA,\
                  vxo1, NA, vyo1, NA, vzo1, NA,\
                  vxo2, NA, vyo2, NA, vzo2, NA 
		  
		  //FxLong, NA, FyLong, NA, FzLong, NA,
#else
#define SAVE_LIST rx, NA, ry, NA, rz, NA,\
                  vx, NA, vy, NA, vz, NA,\
                  Fx, NA, Fy, NA, Fz, NA,\
                  vxo1, NA, vyo1, NA, vzo1, NA,\
                  vxo2, NA, vyo2, NA, vzo2, NA 
#endif
#undef  EXT_SLST
#define EXT_SLST  &s, &s1, &s2, &Vol, &Vol1, &Vol2, &Vol1o1, &s1o1, &Vol1o2,\
                  &s1o2

/* Reduced list of variables to save on xva file (tape file) */ 
#define XVA_LIST rx, NA, ry, NA, rz, NA,\
                 vx, NA, vy, NA, vz, NA

#define XVA_NUM 15 /* this is the number of vars in XVA_LIST <--- SET THIS!!!*/

#define XVA_ALST &rx[0], NA, &ry[0], NA, &rz[0], NA,\
                 &vx[0], NA, &vy[0], NA, &vz[0], NA 

#define XVA_DLST *rx[NA], *ry[NA], *rz[NA],\
                 *vx[NA], *vy[NA], *vz[NA]
 
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
#ifdef MD_RESPA
#define ALLOC_LIST  &rx[0], NA, &ry[0], NA, &rz[0], NA,\
                    &vx[0], NA, &vy[0], NA, &vz[0], NA,\
                    &Fx[0], NA, &Fy[0], NA, &Fz[0], NA,\
		    &FxLong[0], NA, &FyLong[0], NA, &FzLong[0], NA,\
                    &vxt[0], NA, &vyt[0], NA, &vzt[0], NA,\
                    &FxLJ[0], NA, &FyLJ[0], NA, &FzLJ[0], NA,\
                    &vxo1[0], NA, &vyo1[0], NA, &vzo1[0], NA,\
                    &vxo2[0], NA, &vyo2[0], NA, &vzo2[0], NA,\
		    &vxg[0], NA,  &vyg[0], NA,  &vzg[0], NA,\
		    &vxt2[0], NA,  &vyt2[0], NA,  &vzt2[0], NA,\
		    &Rx, 1, &Ry, 1, &Rz, 1, &px[0], NA, &py[0], NA, &pz[0], NA
#else
#define ALLOC_LIST  &rx[0], NA, &ry[0], NA, &rz[0], NA,\
                    &vx[0], NA, &vy[0], NA, &vz[0], NA,\
                    &Fx[0], NA, &Fy[0], NA, &Fz[0], NA,\
                    &vxt[0], NA, &vyt[0], NA, &vzt[0], NA,\
                    &FxLJ[0], NA, &FyLJ[0], NA, &FzLJ[0], NA,\
                    &vxo1[0], NA, &vyo1[0], NA, &vzo1[0], NA,\
                    &vxo2[0], NA, &vyo2[0], NA, &vzo2[0], NA,\
		    &vxg[0], NA,  &vyg[0], NA,  &vzg[0], NA,\
		    &vxt2[0], NA,  &vyt2[0], NA,  &vzt2[0], NA,\
		    &Rx, 1, &Ry, 1, &Rz, 1
#endif
/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#ifdef MD_RESPA
#define DECL_LIST   *rx[NA], *ry[NA], *rz[NA],\
                    *vx[NA], *vy[NA], *vz[NA],\
                    *vxt[NA], *vyt[NA], *vzt[NA],\
                    *Fx[NA], *Fy[NA], *Fz[NA],\
		    *FxLong[NA], *FyLong[NA], *FzLong[NA],\
                    *FxLJ[NA], *FyLJ[NA], *FzLJ[NA],\
                    *vxo1[NA], *vyo1[NA], *vzo1[NA],\
                    *vxg[NA], *vyg[NA], *vzg[NA],\
                    *vxt2[NA], *vyt2[NA], *vzt2[NA],\
                    *vxo2[NA], *vyo2[NA], *vzo2[NA],\
                    *Rx, *Ry, *Rz, *px[NA], *py[NA], *pz[NA]         
#else
#define DECL_LIST   *rx[NA], *ry[NA], *rz[NA],\
                    *vx[NA], *vy[NA], *vz[NA],\
                    *vxt[NA], *vyt[NA], *vzt[NA],\
                    *Fx[NA], *Fy[NA], *Fz[NA],\
                    *FxLJ[NA], *FyLJ[NA], *FzLJ[NA],\
                    *vxo1[NA], *vyo1[NA], *vzo1[NA],\
                    *vxg[NA], *vyg[NA], *vzg[NA],\
                    *vxt2[NA], *vyt2[NA], *vzt2[NA],\
                    *vxo2[NA], *vyo2[NA], *vzo2[NA],\
                    *Rx, *Ry, *Rz        
#endif				   
#undef EXT_DLST
#define EXT_DLST    s, s1, s2, Vol, Vol1, Vol2, Vol1o1, s1o1, Vol1o2, s1o2

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
  char endfile[NAME_LENGTH];
  int endFormat; /* 0 = binary 1 = ascii 2 = both */
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
  char tmpPath[NAME_LENGTH];
  char misPath[NAME_LENGTH];
 
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
  int mdseed; 
  /* Accumulators for calculating the rotational diffusion coefficent */
  COORD_TYPE sumox[MAXPAR];
  COORD_TYPE sumoy[MAXPAR];
  COORD_TYPE sumoz[MAXPAR];
#if 0 
  COORD_TYPE vcmx0[MAXPAR];
  COORD_TYPE vcmy0[MAXPAR];
  COORD_TYPE vcmz0[MAXPAR];
  COORD_TYPE ux0[MAXPAR];
  COORD_TYPE uy0[MAXPAR];
  COORD_TYPE uz0[MAXPAR];
  COORD_TYPE ox0[MAXPAR];
  COORD_TYPE oy0[MAXPAR];
  COORD_TYPE oz0[MAXPAR];
#endif  
  COORD_TYPE rxCMi[MAXPAR]; /* initial coordinates of center of mass */
  COORD_TYPE ryCMi[MAXPAR]; /* MAXPAR is the maximum number of particles */
  COORD_TYPE rzCMi[MAXPAR];
  double DR[MAXPAR][3];
  int ipart;            /* Particle to follow in its motion inside the box */
  COORD_TYPE Q;
  COORD_TYPE W;
  int Nose; /* if Nose=1 use Nose method, otherwise not */
  int sResetSteps; /* Steps at which reset s to 1 */
  int savedXva;
  int CMreset;
#ifdef MD_RESPA
  int nebrTabFacLong;                /* How much storage sould be provided for 
     					the neighbour list (see Rapaport pag.53
     					for details )*/
  COORD_TYPE rNebrShellLong;         /* = Dr see Rapaport pag. 53 */
  int nrespa;                    /* numero di iterazioni del reference system */
  double rcutInner;
  double lambda;
  //int keepInvMat;
#endif
#ifdef MD_RAPACONSTR
  double rapatol;
#endif
  int nebrTabFac;                /* How much storage sould be provided for 
                                   the neighbour list (see Rapaport pag.53
				   for details )*/
  COORD_TYPE rNebrShell;         /* = Dr see Rapaport pag. 53 */

  int noLinkedList;              /* If true use neighbour list method without
				    linked list */
  COORD_TYPE avVol;
  COORD_TYPE avs;
  COORD_TYPE avVol1;
  COORD_TYPE tolVol;
  COORD_TYPE tols;
  COORD_TYPE tolVol1;
  int snapSteps;
  int grow;
  double growth;
  int snapmode;
  int HNBOX;
  int avngS;
  int avnggr;
  int avngEta;
  int avngTemp;
  int avngPress;
  int avngMB;
  int bakSaveMode;/* save mode for ascii backup */
  int chkeqstps; /* ogni chkeqstps controlla se il sistema è equilibrato */
  int equilibrated; /* 1 = system equilibrated */
  char nRun[132];

  /* ADD 13/4/2000 
     Logarithmic saving implemented for xva file */
  int xvaSaveMode;/* 0 = linear 1 = semilog 2 = bilog (not impl. yet) */
  double base;    /* We save at base^^NN step */
  int NN;         /* Logatithmic block length */
  double fstps;         /* There are KK block each base^NN long */
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
  COORD_TYPE steplength;		/* temporal step length */
  int totStep;	/* temporal step number that simulation 
				   must do */
  int curStep;	/* current step of simulation */
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
#if 1
  COORD_TYPE d;                 /* distance between atoms */
#endif
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m[NA];                 /* atom mass */
  int nsites;
  COORD_TYPE rcut;              /* cutoff for the pair potential */ 
  int NN;
  int MM;
  double epsilon;               /* dieletric constant */
  double sigma;
#if 0
  double chsi;                  /* friction coefficient (brownian motion) */
#endif
#ifdef MD_FENE
  double kfe;
  double R0;
#endif
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
  
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
#ifdef MAIN
double mass;
#else
extern double mass;
#endif
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
  {"steplength", &Oparams.steplength,         CT},
  {"stepnum",    &Oparams.totStep,            INT},
  {"inifile" ,   &OprogStatus.inifile,        STR},
  {"endfile" ,   &OprogStatus.endfile,        STR},
  {"xvafile" ,   &OprogStatus.xvafile,        STR},
  {"meas0name",  &OprogStatus.dataFiles[0],   STR}, /* files to store 
						       measures */
  {"inistep" ,   &Oparams.curStep,            INT},
  {"tapeTimes",  &OprogStatus.tapeTimes,      INT},
  {"endFormat",  &OprogStatus.endFormat,      INT},
  {"xvaSaveMode",&OprogStatus.xvaSaveMode,  INT},
  {"bakSaveMode", &OprogStatus.bakSaveMode, INT},
  {"bakSteps",   &OprogStatus.bakSteps,       INT},
  {"bakStepsAscii",&OprogStatus.bakStepsAscii, INT},
  {"staSteps",   &OprogStatus.staSteps,       INT},
  {"tmpPath",    OprogStatus.tmpPath,       STR},
  {"misPath",    OprogStatus.misPath,       STR},
  {"xvaSteps",   &OprogStatus.xvaSteps,       INT},  
  {"energySteps",&OprogStatus.measSteps[0],   INT},
  {"energyCalc", &OprogStatus.measCalc[0],    INT},
  {"energyName", &OprogStatus.dataFiles[0],   STR},
  {"DtrSteps",   &OprogStatus.measSteps[1],   INT}, /* steps between measure
						       savings */
  {"DtrCalc",    &OprogStatus.measCalc[1],    INT}, /* steps between measure
						       calculation */
  {"DtrName",    &OprogStatus.dataFiles[1],   STR},
  {"tempSteps",  &OprogStatus.measSteps[2],   INT},
  {"tempCalc",   &OprogStatus.measCalc[2],    INT},
  {"tempName",   &OprogStatus.dataFiles[2],   STR},
  {"grSteps",    &OprogStatus.measSteps[3],   INT},
  {"grCalc",     &OprogStatus.measCalc[3],    INT},
  {"grName",     &OprogStatus.dataFiles[3],   STR},
  {"MBSteps",    &OprogStatus.measSteps[4],   INT},
  {"MBCalc",     &OprogStatus.measCalc[4],    INT},
  {"MBName",     &OprogStatus.dataFiles[4],   STR},
  {"DrotSteps",  &OprogStatus.measSteps[5],   INT},
  {"DrotCalc",   &OprogStatus.measCalc[5],    INT},
  {"DrotName",   &OprogStatus.dataFiles[5],   STR},
  {"etaSteps",   &OprogStatus.measSteps[6],   INT},
  {"etaCalc",    &OprogStatus.measCalc[6],    INT},
  {"etaName",    &OprogStatus.dataFiles[6],   STR},
  {"SSteps",     &OprogStatus.measSteps[7],   INT},
  {"SCalc",      &OprogStatus.measCalc[7],    INT},
  {"SName",      &OprogStatus.dataFiles[7],   STR},
  {"pressSteps", &OprogStatus.measSteps[8],   INT},
  {"pressCalc",  &OprogStatus.measCalc[8],    INT},   
  {"pressName",  &OprogStatus.dataFiles[8],   STR},
  {"VolSteps",   &OprogStatus.measSteps[9],   INT},
  {"VolCalc",    &OprogStatus.measCalc[9],    INT},
  {"VolName",    &OprogStatus.dataFiles[9],   STR},
  {"sSteps",     &OprogStatus.measSteps[10],  INT},
  {"sName",      &OprogStatus.dataFiles[10],  STR},
  {"sCalc",      &OprogStatus.measCalc[10],   INT},
  {"PtensSteps", &OprogStatus.measSteps[11],  INT},
  {"PtensCalc",  &OprogStatus.measCalc[11],   INT},   
  {"PtensName",  &OprogStatus.dataFiles[11],  STR},
  {"DQtensSteps",&OprogStatus.measSteps[12],  INT},
  {"DQtensCalc", &OprogStatus.measCalc[12],   INT},
  {"DQtensName", &OprogStatus.dataFiles[12],  STR},
  {"DphiSqName", &OprogStatus.dataFiles[13],  STR},
  {"DphiSqSteps",&OprogStatus.measSteps[13],  INT},
  {"VSteps", &OprogStatus.measSteps[14],      INT},
  {"VName",  &OprogStatus.dataFiles[14],       STR},
  {"VCalc",  &OprogStatus.measCalc[14],        INT},
#if 0
  {"C1Name",     &OprogStatus.dataFiles[14],  STR},
  {"C1Steps",    &OprogStatus.measSteps[14],  INT},
  {"C2Name",     &OprogStatus.dataFiles[15],  STR},
  {"C2Steps",    &OprogStatus.measSteps[15],  INT},
  {"C3Name",     &OprogStatus.dataFiles[16],  STR},
  {"C3Steps",    &OprogStatus.measSteps[16],  INT},
  {"C4Name",     &OprogStatus.dataFiles[17],  STR},
  {"C4Steps",    &OprogStatus.measSteps[17],  INT},
  {"psi1Name",   &OprogStatus.dataFiles[18],  STR},
  {"psi1Steps",  &OprogStatus.measSteps[18],  INT},
  {"psi2Name",   &OprogStatus.dataFiles[19],  STR},
  {"psi2Steps",  &OprogStatus.measSteps[19],  INT},
  {"velacfName", &OprogStatus.dataFiles[20],  STR},
  {"velacfSteps",&OprogStatus.measSteps[20],  INT},
#endif
  {"ipart",      &OprogStatus.ipart,          INT},
#if 0
  {"nonGaussName",&OprogStatus.dataFiles[21], STR},
  {"nonGaussSteps",&OprogStatus.measSteps[21],INT},
  {"DrSqName",   &OprogStatus.dataFiles[22],  STR},
  {"DrSqSteps",  &OprogStatus.measSteps[22],  INT},
  {"sqrtdr2Name",   &OprogStatus.dataFiles[23],  STR},
  {"sqrtdr2Steps",  &OprogStatus.measSteps[23],  INT},
  {"GsName",        &OprogStatus.dataFiles[24],  STR},
  {"GsSteps",       &OprogStatus.measSteps[24],  INT},
  {"GsAngName",     &OprogStatus.dataFiles[25],  STR},
  {"GsAngSteps",    &OprogStatus.measSteps[25],  INT},
  {"FsName",        &OprogStatus.dataFiles[26],  STR},
  {"FsSteps",       &OprogStatus.measSteps[26],  INT},
  {"GsGsgName",     &OprogStatus.dataFiles[27],  STR},
  {"GsGsgSteps",    &OprogStatus.measSteps[27],  INT},
#endif
  {"initStep6",  &OprogStatus.initStep[6],    INT},
  {"initCalc6",  &OprogStatus.initCalc[6],    INT}, 
  {"Q",          &OprogStatus.Q,              CT},
  {"grow",       &OprogStatus.grow,           INT},
  {"growth",     &OprogStatus.growth,         CT},
  {"Nose",       &OprogStatus.Nose,           INT},
  {"sResetSteps",&OprogStatus.sResetSteps,    INT},
  {"snapSteps",  &OprogStatus.snapSteps,      INT},
  {"snapmode",   &OprogStatus.snapmode,       INT},
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"noLinkedList",&OprogStatus.noLinkedList,  INT},
#ifdef MD_RESPA
  {"rNebrShellLong", &OprogStatus.rNebrShellLong,     CT},
  {"nebrTabFacLong", &OprogStatus.nebrTabFacLong,     INT},
#endif
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
#ifdef MD_RESPA
  {"rNebrShellLong", &OprogStatus.rNebrShellLong,     CT},
  {"nebrTabFacLong", &OprogStatus.nebrTabFacLong,     INT},
#endif
#ifdef MD_RAPACONSTR
  {"rapatol",       &OprogStatus.rapatol,     CT},
#endif
  {"W",          &OprogStatus.W,              CT},
  {"P",          &Oparams.P,                  CT},
#if 0
  {"nsites",     &Oparams.nsites,            INT},
#endif
  {"Vol",        &Vol,                        CT},
  {"avVol",      &OprogStatus.avVol,          CT},
  {"avs",        &OprogStatus.avs,            CT},
  {"tolVol",     &OprogStatus.tolVol,         CT},
  {"tols",       &OprogStatus.tols,           CT},
  {"tolVol1",    &OprogStatus.tolVol1,        CT},
  {"avVol1",     &OprogStatus.avVol1,         CT},
  {"avngEta",    &OprogStatus.avngEta,        INT},
  {"avngTemp",   &OprogStatus.avngTemp,       INT},
  {"avngPress",  &OprogStatus.avngPress,      INT},
  {"avngS",      &OprogStatus.avngS,          INT},
  {"avnggr",     &OprogStatus.avnggr,         INT},
  {"avngMB",     &OprogStatus.avngMB,         INT},
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
#if 0
  {"zplus",      &Oparams.zplus,      CT},
  {"zminus",     &Oparams.zminus,      CT},
  {"nplus",      &Oparams.nplus,      CT},
  {"nminus",     &Oparams.nminus,     CT},
#endif
  {"seed",       &OprogStatus.mdseed,                   INT},
  {"m",          &Oparams.m[0],             CT},
#ifdef MD_FENE
  {"kfe",        &Oparams.kfe,              CT},
  {"R0",         &Oparams.R0,               CT},
#endif
#if 0
  {"chsi",       &Oparams.chsi,          CT},
#endif
  {"cellNum",    &Oparams.M,                INT},
  {"epsilon",    &Oparams.epsilon,      CT},
  {"sigma",      &Oparams.sigma,        CT},
  {"d",          &Oparams.d,            CT},
  {"rcut",       &Oparams.rcut,             CT},
#ifdef MD_RESPA
  {"rcutInner",  &OprogStatus.rcutInner,            CT},
  {"nrespa",     &OprogStatus.nrespa,      INT},
  {"lambda",     &OprogStatus.lambda,      CT},
  //{"keepInvMat", &OprogStatus.keepInvMat,  INT},
#endif
#if 0
  {"atomsDist",  &Oparams.d,                CT},
#endif
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"temperat",   &Oparams.T,                CT},
  {"tol",        &Oparams.tol,              CT},
    /* parametri per scegliere il tipo di salvataggio del file xva
     (lineare o semilog) */
  {"xvaSaveMode",&OprogStatus.xvaSaveMode,  INT},
  {"base",       &OprogStatus.base,         CT},
  {"NN",         &OprogStatus.NN,           INT},
  {"N",          &Oparams.NN,               INT},
  {"M",          &Oparams.MM,               INT},
  /*Parametri per attivare il parallel tempering */
  {"nRun",       OprogStatus.nRun,         STR},
  /* ======================================================================= */
  {"", NULL, 0} /* end of list, don't touch this !!! */
};
#else 
extern struct singlePar OsinglePar[];
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
  {"endFile",      OS(endfile),                     1,    NAME_LENGTH,   "%s"},
  {"endFormat",    &OS(endFormat),                  1,               1, "%d"},
  {"dataFiles",    OS(dataFiles),         NUM_MISURE,     NAME_LENGTH, "%s"},
  {"xvaSteps",     &OS(xvaSteps),                  1,               1, "%d"},
  {"bakSteps",     &OS(bakSteps),                  1,               1, "%d"},
  {"bakStepsAscii",&OS(bakStepsAscii),             1,               1, "%d"},
  {"staSteps",     &OS(staSteps),                  1,               1, "%d"},
  {"measSteps",    OS(measSteps),         NUM_MISURE,               1, "%d"},
  {"measCalc",     OS(measCalc),          NUM_MISURE,               1, "%d"},
  {"initCalc",     OS(initCalc),          NUM_MISURE,               1, "%d"},
  {"tmpPath",      OS(tmpPath),                     1,  NAME_LENGTH, "%s"},
  {"misPath",      OS(misPath),                     1,  NAME_LENGTH, "%s"},
  {"initStep",     OS(initStep),          NUM_MISURE,               1, "%d"},
  {"sumEta",       &OS(sumEta),                     1,              1, "%.10G"},
  {"DQxy",         &OS(DQxy),                       1,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                       1,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                       1,              1, "%.10G"},
  {"PxyArr",       OS(PxyArr),                      5,              1, "%.10G"},
  {"PyzArr",       OS(PyzArr),                      5,              1, "%.10G"},
  {"PzxArr",       OS(PzxArr),                      5,              1, "%.10G"},
  {"hist",         OS(hist),                  MAXBIN,               1, "%d"},
  {"sumS",         OS(sumS),                    NUMK,               1, "%.6G"},
  {"histMB",       OS(histMB),                  NUMV,               1, "%d"},
  {"sumTemp",      &OS(sumTemp),                    1,              1, "%.6G"},
  //mancano le rxi per i coeff. di diffusione!!!!
  {"sumPress",     &OS(sumPress),                   1,              1, "%.6G"},
  {"sumox",        OS(sumox),          MAXPAR,               1, "%.15G"},
  {"sumoy",        OS(sumoy),          MAXPAR,               1, "%.15G"},
  {"sumoz",        OS(sumoz),          MAXPAR,               1, "%.15G"},
  {"rxCMi",        OS(rxCMi),          MAXPAR,               1, "%.15G"},
  {"ryCMi",        OS(ryCMi),          MAXPAR,               1, "%.15G"},
  {"rzCMi",        OS(rzCMi),          MAXPAR,               1, "%.15G"},
/// {"DR",           OS(DR),             MAXPAR,               3, "%.15G"},
#if 0 
  {"vcmx0",        OS(vcmx0),          MAXPAR,               1, "%.15G"},
  {"vcmy0",        OS(vcmy0),          MAXPAR,               1, "%.15G"},
  {"vcmz0",        OS(vcmz0),          MAXPAR,               1, "%.15G"},
#endif
  {"Q",            &OS(Q),                          1,              1 , "%.6G"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
  {"Nose",         &OS(Nose),                       1,   1,  "%d"},
  //  {"zeroall_s",    &OS(zeroall_s),                  1,   1,  "%d"},
  {"sResetSteps",  &OS(sResetSteps),                1,   1,  "%d"},
  {"snapSteps",    &OS(snapSteps),                  1,   1,  "%d"},
  {"Cmreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.6G"},
#ifdef MD_RESPA
  {"nebrTabFacLong",   &OS(nebrTabFacLong),                 1,   1,   "%d"},
  {"rNebrShellLong",   &OS(rNebrShellLong),                 1,   1, "%.6G"},
  {"lambda",           &OS(lambda),                         1,   1, "%.10G"},
  {"nrespa",           &OS(nrespa),                         1,   1,  "%d"},
  {"rcutInner",        &OS(rcutInner),                      1,   1, "%.10G"},
  //{"keepInvMat",       &OS(keepInvMat),                     1,   1,  "%d"},
#endif
#ifdef MD_RAPACONSTR
  {"rapatol",     &OS(rapatol),                      1,  1, "%.15G"},
#endif
  {"noLinkedList", &OS(noLinkedList),               1,   1,  "%d"},
  {"avVol",          &OS(avVol),                    1,   1, "%.8G"},
  {"avVol1",         &OS(avVol1),                   1,   1, "%.8G"},
  {"tolVol",         &OS(tolVol),                   1,   1, "%.8G"},
  {"tolVol1",        &OS(tolVol1),                  1,   1, "%.8G"},
  {"avs",          &OS(avs),                        1,   1, "%.8G"},
  {"tols",         &OS(tols),                       1,   1, "%.8G"},
  {"HNBOX",        &OS(HNBOX),                      1,   1,  "%d"},
  {"avngTemp",     &OS(avngTemp),                   1,   1,  "%d"},
  {"avngPress",    &OS(avngPress),                  1,   1,  "%d"},
  {"avngEta",      &OS(avngEta),                    1,   1,  "%d"},
  {"avnggr",       &OS(avnggr),                      1,  1,  "%d"},
  {"avngS",        &OS(avngS),                      1,   1,  "%d"},
  //  {"avgMedia",     &OS(avgMedia),                   1,     1, "%d"},
  {"xvaSavedMode", &OS(xvaSaveMode),                1,  1,    "%d"},
  {"bakSaveMode",  &OS(bakSaveMode),                1,  1,    "%d"},
  {"chkeqstps",    &OS(chkeqstps),                  1,  1,    "%d" },
  {"equilibrated", &OS(equilibrated),               1,  1,    "%d"},
  {"bakSavedMode", &OS(bakSaveMode),                1,  1,    "%d"},
  {"base",         &OS(base),                       1,  1, "%.6G"},
  {"NN",           &OS(NN),                         1,  1,   "%d"},
  {"fstps",        &OS(fstps),                      1,  1,   "%.15G"},
  {"nRun",         OS(nRun),                       1,  1,   "%s"},
  {"", NULL, 0, 0, ""}
};
#else
extern struct pascii opro_ascii[];
#endif

#ifdef MAIN
/* ========================= >>> opar_ascii <<< =========================== */
struct pascii opar_ascii[]=
{
  {"parnum",            &OP(parnum),                       1,   1, "%d"},
  {"steplength",        &OP(steplength),                  1,   1, "%.6G"},
  {"nsites",            &OP(nsites),                     1,   1,  "%d"},
  {"totStep",           &OP(totStep),                     1,   1,  "%d"},
  {"curStep",           &OP(curStep),                     1,   1,  "%d"},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
  {"m",                 OP(m),                            3,   1, "%.6G"},
  {"NN",                &OP(NN),                          1,   1,  "%d"},
  {"MM",                &OP(MM),                          1,   1,  "%d"},
  {"rcut",              &OP(rcut),                        1,   1, "%.8G"},
  {"epsilon",           &OP(epsilon),                     1,   1, "%.10G"},
  {"sigma",             &OP(sigma),                       1,   1, "%.8G"},
  {"d",                 &OP(d),                           1,   1, "%.15G"},
#if 0
  {"chsi",              &OP(chsi),                        1,   1, "%.10G"},
#endif
#if 0
  {"epsilon",           OP(epsab),                        2,   2, "%.8G"},
#endif
  {"equilibrat",         &OP(equilibrat),                  1,   1,   "%d"},
  {"M",                 &OP(M),                           1,   1,  "%d"},
  {"tol",               &OP(tol),                         1,   1, "%.15G"},
  {"", NULL, 0, 0, ""}
};
#else
extern struct pascii opar_ascii[];
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
#define REND 57.0
#define HNBOX 0
#define NUMK2AV 150 /* Number of k on the same modulus over which we must 
		       perform a mean */

#ifdef MAIN
#ifdef MD_RESPA
double VLong, VcLong;
#endif
COORD_TYPE Vfe, Vc, V, E, Dtrans, temp, S[NUMK], dummy, eta, Drot, gr[MAXBIN], invs, press,
  gr23[MAXBIN], gr33[MAXBIN], press_m, press_at;
COORD_TYPE Ptens[3], DQtens[3], C1c, C2c, C3c, C4c, velc, Gs[GSPOINT], 
  psi1c, psi2c, sqrtdr2, Aa, DrSqTot, DphiSq, GsAng[GSANGPOINT], Fs[FSPOINT],
  GsGsg[GSPOINT], temp_transl;
#ifdef MD_RESPA_NPT
double Ps, Pv;
#endif
int MB[NUMV];
COORD_TYPE *Dphix, *Dphiy, *Dphiz;/* Time integrals of angulars velocity
				     components */
#else 
#ifdef MD_RESPA
extern double VLong, VcLong;
#endif
extern COORD_TYPE Vfe, Vc, V, E, Dtrans, temp, S[NUMK], dummy, eta, Drot, gr[MAXBIN], 
  gr23[MAXBIN], gr33[MAXBIN], invs,
  press, press_m, press_at, C1c, C2c, C3c, C4c, velc, psi1c, psi2c, 
  Gs[GSPOINT], 
  GsAng[GSANGPOINT], Fs[FSPOINT], GsGsg[GSPOINT], temp_transl;
#ifdef MD_RESPA_NPT
extern double Ps, Pv;
#endif
extern COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, Aa, DrSqTot,
  DphiSq;
extern int MB[NUMV];
extern COORD_TYPE *Dphix, *Dphiy, *Dphiz;/* Time integrals of angulars velocity
					    components */
#endif

/* ============= >>> PUT HERE MEASURING FUNCTION PROTOTYPES <<< ============*/

void energy(void);
void transDiff(void);
void temperat(void);
void SkCM(void);
void maxwell(void);
void rotDiff(void);
void viscosity(void);
void radDens(void);
void Ptensor(void);
void DQtensor(void);
void corrFuncs(void);
void vanHoveSelf(void);
void vanHoveAng(void);
void intermScatt(void);
void gaussapprox(void);
void radDens23(void);
void radDens33(void);
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
  {gr,    sizeof(COORD_TYPE)*MAXBIN,  radDens     },
  {MB,      sizeof(int)*NUMV,         maxwell     },
  {&Drot,   sizeof(COORD_TYPE),       rotDiff     },
  {&eta,    sizeof(COORD_TYPE),       viscosity   },
  {S,       sizeof(COORD_TYPE)*NUMK,  SkCM        },
  {&press,  sizeof(COORD_TYPE),       NULL        },
  {&Vol,    sizeof(COORD_TYPE),       NULL        },
  {&s,      sizeof(COORD_TYPE),       NULL        },
  {Ptens,   sizeof(COORD_TYPE)*3,     Ptensor     },
  {DQtens,  sizeof(COORD_TYPE)*3,     DQtensor    },
  {&DphiSq, sizeof(COORD_TYPE),       NULL        },
  {&Vc,     sizeof(COORD_TYPE),       NULL},
  {gr23,    sizeof(COORD_TYPE)*MAXBIN,radDens23   },
  {gr33,    sizeof(COORD_TYPE)*MAXBIN,radDens33   },
/* {&C1c,    sizeof(COORD_TYPE),       corrFuncs  },
  {&C2c,    sizeof(COORD_TYPE),       NULL       },
  {&C3c,    sizeof(COORD_TYPE),       NULL       },
  {&C4c,    sizeof(COORD_TYPE),       NULL       },
  {&psi1c,  sizeof(COORD_TYPE),       NULL       },
  {&psi2c,  sizeof(COORD_TYPE),       NULL       },
  {&velc,   sizeof(COORD_TYPE),       NULL       },
  {&Aa,     sizeof(COORD_TYPE),       NULL       },
  {&DrSqTot,sizeof(COORD_TYPE),       NULL       },
  {&sqrtdr2, sizeof(COORD_TYPE),      NULL       },
  {Gs,      sizeof(COORD_TYPE)*GSPOINT,    vanHoveSelf},
  {GsAng,   sizeof(COORD_TYPE)*GSANGPOINT, vanHoveAng},
  {Fs,      sizeof(COORD_TYPE)*FSPOINT,    intermScatt},
  {GsGsg,   sizeof(COORD_TYPE)*GSPOINT,    gaussapprox},*/
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
  int saveSteps;
  int size;  /* size in bytes of each measure */
  /* ADD 14/09/200 */
  double dt;
  double T;
  double Vol;
  int N;

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
  int parnum;
  int saveSteps; /* save between two tape savings */
  /* ADd 14/09/2000 */
  int mode;
  int NN;
  double T;
  double Vol;
  double base;
  double dt;     /* dt of each steps */
  int N;

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
void writeGsAng(FILE* afs, int step, int size);
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
/* {"Gs",     writeGs  },
  {"GsAng",  writeGsAng},*/\
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
