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

#define MD_HD_TMP "/home/demichel/simdat/mdtmp/"
/* directory to store temporary files */ 

#define MD_HD_MIS "/home/demichel/simdat/measures/" 
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

#define C_T COORD_TYPE

/* Numero di sottosistemi */
#define MAX_M 16

#define NK 10000
#define NA 1 /* number of atoms for each molecule (particle) */

#define MAXPAR 3000      /* maximum number of simulated particles */

#define NUM_PAR 500   /* Number of particles for the simulation */
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

/* ADDED  22/09/2000: PARALLEL TEMPERING  */
#define PE_POINTS 300 /* Punti della distribuzione P(E,Beta) */
#define SOGLIA_PTEQ 0.15

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
#define SAVE_LIST rx, ry, rz, vx, vy, vz,\
                  Fx, Fy, Fz, vxo1, vyo1, vzo1,\
                  vxo2, vyo2, vzo2
#undef  EXT_SLST
#define EXT_SLST  &Vol, &Vol1, &Vol2, &Vol1o1, &s1o1, &Vol1o2

#undef EXT_PT_SLST
#define EXT_PT_SLST s, s1, s2, s1o1, s1o2  

/* Reduced list of variables to save on xva file (tape file) */ 
#define XVA_LIST rx, ry, rz, vx, vy, vz

#define XVA_NUM 6 /* this is the number of vars in XVA_LIST <--- SET THIS!!!*/

#define XVA_ALST &rx, &ry, &rz, &vx, &vy, &vz

/* MODIFIED 25/09/2000: Parallel tempering 
   aggiunta una * perche' ora abbiamo M sottosistemi
 */
#define XVA_DLST **rx, **ry, **rz, **vx, **vy, **vz

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
#define ALLOC_LIST  &rx, &ry, &rz, &vx, &vy, &vz,\
                    &Fx, &Fy, &Fz, &vxt, &vyt, &vzt,\
                    &FxLJ, &FyLJ, &FzLJ,\
                    &vxo1, &vyo1, &vzo1,\
                    &vxo2, &vyo2, &vzo2,\
		    &vxg,  &vyg,  &vzg,\
	            &vxt2,  &vyt2,  &vzt2
                   
/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#define DECL_LIST   **rx, **ry, **rz, **vx, **vy, **vz, **vxt, **vyt, **vzt,\
                    **Fx, **Fy, **Fz, **FxLJ, **FyLJ, **FzLJ,\
                    **vxo1, **vyo1, **vzo1, **vxg, **vyg, **vzg,\
                    **vxt2, **vyt2, **vzt2, **vxo2, **vyo2, **vzo2
				   
#undef EXT_DLST
#define EXT_DLST  Vol, Vol1, Vol2, Vol1o1, Vol1o2
/*questi sono array con un numero di elementi pari ai sottosistemi simulati */

#undef EXT_PT_DLST
#define EXT_PT_DLST s[MAX_M], s1[MAX_M], s2[MAX_M], s1o1[MAX_M], s1o2[MAX_M]

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
  int iniFormat;
  int endFormat;

  int nRun;
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

  COORD_TYPE sumEta[MAX_M]; /* accumulators for obtaining the mean value of eta */
  
  /* Accumulators for the integral of angular velocity */
  COORD_TYPE DQxy[MAX_M];
  COORD_TYPE DQyz[MAX_M];
  COORD_TYPE DQzx[MAX_M];

  COORD_TYPE PxyArr[MAX_M][5];
  COORD_TYPE PyzArr[MAX_M][5];
  COORD_TYPE PzxArr[MAX_M][5];
  
  /* Accumulator for the radial distribution function */
  int hist[MAX_M][MAXBIN];
  
  /* Accumulator for the static structure factor */ 
  COORD_TYPE sumS[MAX_M][NUMK];

  /* Accumulator for the MAXWELL-BOLTZMANN distribution */ 
  int histMB[NUMV];

  COORD_TYPE sumTemp;
  COORD_TYPE sumPress[MAX_M];
  C_T sumMedia;
  
  //COORD_TYPE vcmx0[MAXPAR];
  //COORD_TYPE vcmy0[MAXPAR];
  //COORD_TYPE vcmz0[MAXPAR];
  
  //COORD_TYPE rxCMi[MAXPAR]; /* initial coordinates of center of mass */
  //COORD_TYPE ryCMi[MAXPAR]; /* MAXPAR is the maximum number of particles */
  //COORD_TYPE rzCMi[MAXPAR];
  /* Di ogni particella mi memorizzo lo spostamento ad una certa 
     temperatura indicata dal primo indice */
  /* Parallel Tempering */
  C_T sumVx[MAX_M][MAXPAR];
  C_T sumVy[MAX_M][MAXPAR];
  C_T sumVz[MAX_M][MAXPAR];

  COORD_TYPE Q;
  COORD_TYPE W;

  C_T alpha;
  C_T S0;
  C_T kmax;
  C_T Dk;

  int savedXva; 
  int Nose; /* if Nose=1 use Nose method, otherwise not */
  int zeroall_s;
  int sResetSteps; /* Steps at which reset s to 1 */
  int CMreset;
  int nebrTabFac;                /* How much storage sould be provided for 
                                   the neighbour list (see Rapaport pag.53
				   for details )*/
  COORD_TYPE rNebrShell;         /* = Dr see Rapaport pag. 53 */

  int noLinkedList;             /* If true use neighbour list method without
				    linked list */
  int AntiCry;
  COORD_TYPE avs;
  COORD_TYPE tols;
  COORD_TYPE tolT;
  int ipart;
  int HNBOX;
  int avngTemp;
  int avngPress;
  int avnggr;
  int avngS;

  /* 12/10/2000 ADDED: Parallel Tempering */
  int avgMedia;
  /* ADD 13/4/2000 
     Logarithmic saving implemented for xva file */
  int xvaSaveMode;/* 0 = linear 1 = semilog 2 = bilog (not impl. yet) */
  int bakSaveMode;/* save mode for ascii backup */
  char tmpPath[NAME_LENGTH];
  char misPath[NAME_LENGTH];

  double base;    /* We save at base^^NN step */
  int NN;         /* Logatithmic block length */
  int fstps;         /* There are KK block each base^NN long */

  /* ADDED [Parallel Tempering] 14/09/2000 */
  int DtPT; /* 1 = Parallel tempering ON 0 = parallel tempering OFF */
  int DtPTEq;
  int srescale;
  /* lambda g relativo al processo (ogni processo ha il suo) */
  double E;
  int PT; /* = 1 vuol dire che si deve usare il Parallel Tempering */
  int PTequilibrated;
  C_T ENmin;
  C_T ENmax;
  C_T sogliaPEij;
  /* distribuzioni delle energie */
  int PE[MAX_M][PE_POINTS];
  char savedSnaps[MAX_M][MAX_M];
  /* il primo indice e' la replica e il secondo la temperatura 
     cioe' savedSnaps[0][2] ci dice se la replica 0 ha salvato 
     la temperatura 2 */
  int scambi[MAX_M];
  int attempted[MAX_M];
  int linearDiff;
  int RW; 
  /* all'inizio e' pari 0 poi quando il sistema a piu' bassa
     temperatura arriva alla temper. piu' alta diventa 1 e quando
     alla fine torna alla temperatura piu' bassa diventa 2 
     in tal caso il sistema si puo' considerare equilibrato
  */
  int refReplica; 
  /* E' la replica di piu' bassa  ossia quella che all'inizio ha
     temperatura piu' bassa */
  int NUM_PEij; 
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

  C_T Vol;                      /* Volume */
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m;             /* atoms masses */
  COORD_TYPE rcut;              /* cutoff for the pair potential */ 
  COORD_TYPE sigma;     /* pair potential length parameters */
  COORD_TYPE epsilon;     /* pair potential energy parameters */
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
  
  COORD_TYPE tol;               /* Tolerance of the shake algoritm used 
				   by RATTLE */
  /* ADDED 26/09/2000 [Parallel Tempering] */
  /* Lambda asseganti all'inizio */
  double lambda0[MAX_M]; 
  /* puntatori a lambda0 che evelvono nel tempo secondo una Markov Chain
     (ved. art. di Kob e Yamamoto */
  int lambdat[MAX_M];
  int PTM; /* numero di sottosistemi per il Parallel Tempering */
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
  {"xvaSteps",     &OS(xvaSteps),                  1,               1, "%d"},
  {"bakSteps",     &OS(bakSteps),                  1,               1, "%d"},
  {"bakStepsAscii",&OS(bakStepsAscii),             1,               1, "%d"},
  {"staSteps",     &OS(staSteps),                  1,               1, "%d"},
  {"measSteps",    OS(measSteps),         NUM_MISURE,               1, "%d"},
  {"measCalc",     OS(measCalc),          NUM_MISURE,               1, "%d"},
  {"initCalc",     OS(initCalc),          NUM_MISURE,               1, "%d"},
  {"initStep",     OS(initStep),          NUM_MISURE,               1, "%d"},
  {"sumEta",       &OS(sumEta),                MAX_M,              1, "%.10G"},
  {"DQxy",         &OS(DQxy),                  MAX_M,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                  MAX_M,              1, "%.10G"},
  {"DQzx",         &OS(DQzx),                  MAX_M,              1, "%.10G"},
  {"PxyArr",       OS(PxyArr),                      MAX_M,         5, "%.10G"},
  {"PyzArr",       OS(PyzArr),                      MAX_M,         5, "%.10G"},
  {"PzxArr",       OS(PzxArr),                      MAX_M,         5, "%.10G"},
  {"hist",         OS(hist),                  MAX_M,          MAXBIN,   "%d"},
  {"sumS",         OS(sumS),                    MAX_M,     NUMK,        "%.6G"},
  {"histMB",       OS(histMB),                  NUMV,               1, "%d"},
  {"sumTemp",      &OS(sumTemp),                    1,              1, "%.6G"},
  {"sumPress",     &OS(sumPress),                   MAX_M,              1, "%.6G"},
  {"sumMedia",     &OS(sumMedia),                   1,   1, "%.6f"},
  {"sumVx",        OS(sumVx),                    MAX_M, MAXPAR, "%.10f"},
  {"sumVy",        OS(sumVy),                    MAX_M, MAXPAR,  "%.10f"},
  {"sumVz",        OS(sumVz),                    MAX_M, MAXPAR,  "%.10f"},
  {"Q",            &OS(Q),                          1,              1 , "%.6G"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
  {"alpha",        &OS(alpha),                      1,              1, "%.6G"},
  {"S0",           &OS(S0),                         1,              1, "%.6G"},
  {"kmax",         &OS(kmax),                       1,   1, "%.6G"},
  {"Dk",           &OS(Dk),                         1,   1, "%.6G"},
  {"savedXva",     &OS(savedXva),                   1,   1,   "%d"},
  {"Nose",         &OS(Nose),                       1,   1,  "%d"},
  {"zeroall_s",    &OS(zeroall_s),                  1,   1,  "%d"},
  {"sResetSteps",  &OS(sResetSteps),                1,   1,  "%d"},
  {"Cmreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.6G"},
  {"noLinkedList", &OS(noLinkedList),               1,   1,  "%d"},
  {"AntiCry",      &OS(AntiCry),                    1,   1,   "%d"},
  {"avs",          &OS(avs),                        1,   1, "%.8G"},
  {"tols",         &OS(tols),                       1,   1, "%.8G"},
  {"tolT",         &OS(tolT),                       1,   1, "%.8G"},
  {"ipart",        &OS(ipart),                      1,   1, "%d"},
  {"HNBOX",        &OS(HNBOX),                      1,   1,  "%d"},
  {"avngTemp",     &OS(avngTemp),                   1,   1,  "%d"},
  {"avngPress",    &OS(avngPress),                  1,   1,  "%d"},
  {"avnggr",       &OS(avnggr),                     1,  1,  "%d"},
  {"avngS",        &OS(avngS),                      1,   1,  "%d"},
  {"avgMedia",     &OS(avgMedia),                   1,     1, "%d"},
  {"xvaSavedMode", &OS(xvaSaveMode),                1,  1,    "%d"},
  {"bakSavedMode", &OS(bakSaveMode),                1,  1,    "%d"},
  {"tmpPath",      OS(tmpPath),                     1,  NAME_LENGTH, "%s"},
  {"misPath",      OS(misPath),                     1,  NAME_LENGTH, "%s"},
  {"base",         &OS(base),                       1,  1, "%.6G"},
  {"NN",           &OS(NN),                         1,  1,   "%d"},
  {"fstps",        &OS(fstps),                      1,  1,   "%d"},
  {"nRun",         &OS(nRun),                       1,  1,   "%d"},
    {"srescale",     &OS(srescale),                   1, 1,    "%d"},
  {"PT",           &OS(PT),                         1,  1,   "%d"},
  {"PTequilibrated",&OS(PTequilibrated),             1,   1,  "%d"},
  {"ENmin",        &OS(ENmin),                      1,  1, "%.6f"},
  {"ENmax",        &OS(ENmax),                      1,  1, "%.6f"},
  {"sogliaPEij",   &OS(sogliaPEij),                 1,  1, "%.10f"},
  {"PE",           OS(PE),                     MAX_M, PE_POINTS,  "%d"},
  {"savedSnaps",   OS(savedSnaps),             MAX_M, MAX_M,      "%d"},
  {"scambi",       OS(scambi),                 MAX_M,     1, "%d"},
  {"attempted",    OS(attempted),              MAX_M,     1, "%d"},
  {"linearDiff",   &OS(linearDiff),                 1,    1, "%d"},
  {"RW",           &OS(RW),                         1,    1,  "%d"},
  {"refReplica",   &OS(refReplica),                 1,    1,  "%d"},
  {"NUM_PEij",     &OS(NUM_PEij),                   1,    1,   "%d"},
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
  {"steplength",        &OP(steplength),                 1,    1, "%.6G"},
  {"totStep",           &OP(totStep),                     1,   1,  "%d"},
  {"curStep",           &OP(curStep),                     1,   1,  "%d"},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
  {"m",                 &OP(m),                           1,   1, "%.6G"},
  {"rcut",              &OP(rcut),                        1,   1, "%.8G"},
  {"sigma",             &OP(sigma),                       1,   1, "%.8G"},
  {"epsilon",           &OP(epsilon),                     1,   1, "%.8G"},
  {"equilibrat",         &OP(equilibrat),                  1,   1,   "%d"},
  {"M",                 &OP(M),                           1,   1,  "%d"},
  {"tol",               &OP(tol),                         1,   1, "%.15G"},
  {"lambda0",           OP(lambda0),                  MAX_M,  "%.15f"},
  {"lambdat",           OP(lambdat),                  MAX_M,     "%d"},
  {"PTM",               &OP(PTM),                         1,     "%d"},
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
COORD_TYPE EXT_PT_DLST;
/* ADDED 25/09/2000: Parallel Tempering */
char PT_temps[MAX_M*4];
#else
extern COORD_TYPE DECL_LIST; 
extern COORD_TYPE EXT_DLST;
extern COORD_TYPE EXT_PT_DLST;
extern char PT_temps[MAX_M*4];
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
  {"steplength", &Oparams.steplength,         CT},
  {"stepnum",    &Oparams.totStep,            INT},
  {"inifile" ,   &OprogStatus.inifile,        STR},
  {"endfile" ,   &OprogStatus.endfile,        STR},
  {"xvafile" ,   &OprogStatus.xvafile,        STR},
  {"inistep" ,   &Oparams.curStep,            INT},
  {"tapeTimes",  &OprogStatus.tapeTimes,      INT},
  {"bakSteps",   &OprogStatus.bakSteps,       INT},
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
  {"tempSteps",  &OprogStatus.measSteps[2],   INT},
  {"tempCalc",   &OprogStatus.measCalc[2],    INT},
  {"tempName",   &OprogStatus.dataFiles[2],   STR},
  {"sSteps",     &OprogStatus.measSteps[3],  INT},
  {"sName",      &OprogStatus.dataFiles[3],  STR},
  {"sCalc",      &OprogStatus.measCalc[3],   INT},
  /* Parallel Tempering */
  {"PESteps",    &OprogStatus.measSteps[4],  INT},
  {"PEName",     &OprogStatus.dataFiles[4],  STR},
  {"PECalc",     &OprogStatus.measCalc[4],   INT},
  {"PEBegin",    &OprogStatus.initStep[4],   INT},
  {"PEijSteps",    &OprogStatus.measSteps[5],  INT},
  {"PEijName",     &OprogStatus.dataFiles[5],  STR},
  {"PEijCalc",     &OprogStatus.measCalc[5],   INT},
  {"PEijBegin",    &OprogStatus.initStep[5],   INT},
  {"mediaSteps",    &OprogStatus.measSteps[6],  INT},
  {"mediaName",     &OprogStatus.dataFiles[6],  STR},
  {"mediaCalc",     &OprogStatus.measCalc[6],   INT},
  {"lambdaSteps",    &OprogStatus.measSteps[7],  INT},
  {"lambdaName",     &OprogStatus.dataFiles[7],  STR},
  {"lambdaCalc",     &OprogStatus.measCalc[7],   INT},
  {"exchgSteps",    &OprogStatus.measSteps[8],  INT},
  {"exchgName",     &OprogStatus.dataFiles[8],  STR},
  {"exchgCalc",     &OprogStatus.measCalc[8],   INT},
  {"Q",          &OprogStatus.Q,              CT},
  {"zeroall_s",  &OprogStatus.zeroall_s,      INT},
  {"Nose",       &OprogStatus.Nose,           INT},
  {"sResetSteps",&OprogStatus.sResetSteps,    INT},
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"noLinkedList",&OprogStatus.noLinkedList,  INT},
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
  {"W",          &OprogStatus.W,              CT},
  {"P",          &Oparams.P,                  CT},
  {"Vol",        &Vol,                        CT},
  {"tolT",       &OprogStatus.tolT,           CT},
  {"avs",        &OprogStatus.avs,            CT},
  {"tols",       &OprogStatus.tols,           CT},
  {"avngTemp",   &OprogStatus.avngTemp,       INT},
  {"avgMedia",   &OprogStatus.avgMedia,       INT},
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
  {"S0",         &OprogStatus.S0,            CT},
  {"AntiCry",    &OprogStatus.AntiCry,      INT},
  {"kmax",       &OprogStatus.kmax,          CT},
  {"alpha",      &OprogStatus.alpha,         CT},
  {"Dk",         &OprogStatus.Dk,            CT},
  {"sigma",      &Oparams.sigma,             CT},
  {"mass",       &Oparams.m,                 CT},
  {"cellNum",    &Oparams.M,                INT},
  {"epsilon",    &Oparams.epsilon,           CT},
  {"rcut",       &Oparams.rcut,              CT},
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"temperat",   &Oparams.T,                 CT},
  {"tol",        &Oparams.tol,               CT},
  /* parametri per scegliere il tipo di salvataggio del file xva
     (lineare o semilog) */
  {"xvaSaveMode",&OprogStatus.xvaSaveMode,  INT},
  {"bakSaveMode", &OprogStatus.bakSaveMode, INT},
  {"tmpPath",    OprogStatus.tmpPath,       STR},
  {"misPath",    OprogStatus.misPath,       STR},
  {"base",       &OprogStatus.base,          CT},
  {"NN",         &OprogStatus.NN,           INT},
  /*Parametri per attivare il parallel tempering */
  /* DtPt e' in stepsa!!!!!!!!!!!!!!!!!!!! */
  /* Parallel Tempering */
  {"DtPT",       &OprogStatus.DtPT,         INT},
  {"DtPTEq",     &OprogStatus.DtPTEq,       INT},
  {"PTtemps",    PT_temps,                  STR},
  {"PTM",        &Oparams.PTM,              INT},
  {"PT",         &OprogStatus.PT,           INT},
  {"ENmin",      &OprogStatus.ENmin,         CT},
  {"ENmax",      &OprogStatus.ENmax,         CT},
  {"srescale",   &OprogStatus.srescale,     INT},
  {"linearDiff", &OprogStatus.linearDiff,   INT},
  {"sogliaPEij", &OprogStatus.sogliaPEij,    CT},
  {"RW",         &OprogStatus.RW,           INT},
  {"NUM_PEij",   &OprogStatus.NUM_PEij,     INT},
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
COORD_TYPE Etot, E[MAX_M], Dtrans[MAX_M], temp, S[NUMK], dummy, eta, gr[MAXBIN], 
  invs, press[MAX_M], diffMedia, probSc[MAX_M], 
  press_m[MAX_M], press_at[MAX_M], PE[MAX_M*PE_POINTS], PEijM[2*MAX_M*PE_POINTS];
COORD_TYPE Ptens[3], DQtens[3], 
  sqrtdr2, Aa, DrSqTot, temp_transl, smed;
int MB[NUMV];
#else 
extern COORD_TYPE Etot, E[MAX_M], Dtrans[MAX_M], temp, S[NUMK], 
  dummy, eta, probSc[MAX_M], 
  gr[MAXBIN], invs, diffMedia,
  press[MAX_M], press_m[MAX_M], press_at[MAX_M], temp_transl;
extern COORD_TYPE Ptens[3], DQtens[3], sqrtdr2, Aa, DrSqTot,
  DphiSq, PE[MAX_M*PE_POINTS], smed, PEijM[2*MAX_M*PE_POINTS];
extern int MB[NUMV];
#endif

/* ============= >>> PUT HERE MEASURING FUNCTION PROTOTYPES <<< ============*/

void energy(void);
void transDiff(void);
void temperat(void);
/* ADDED Parallel Tempering */
void probE(void);
void misuras(void);
void media(void);
void probScMis(void);

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
  
  {&Etot,   sizeof(COORD_TYPE),            energy    },
  {Dtrans,  sizeof(COORD_TYPE)*MAX_M,      transDiff },
  {&temp,   sizeof(COORD_TYPE),            temperat  },
  {&smed,   sizeof(COORD_TYPE),            misuras   },
  /* ADDED Parallel Tempering */
  {PE,      sizeof(C_T)*MAX_M*PE_POINTS,   probE     },
  {PEijM,   2*sizeof(C_T)*MAX_M*PE_POINTS, NULL      },
  {&diffMedia, sizeof(C_T),                media     },
  {Oparams.lambdat,  sizeof(C_T)*MAX_M,    NULL      },
  {probSc,  sizeof(C_T)*MAX_M,             probScMis },
  /* ==================================================================== */
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
  int PTM;   /* numero di repliche */
  int saveSteps;
  int size;  /* size in bytes of each measure */
  double ENmin;
  double ENmax;
  int PEpoints;
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
  /* ADDED Parallel Tempering */
  int PTM;  /* numero di repliche */
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
/* ADDED Parallel Tempering */
void writePE(FILE* afs, int step, int size);
void writePEij(FILE* afs, int step, int size);
/* ==========================================================================*/
struct convStruct OconvStruct[]=
{
  /* Predefined converter, their code is in md2ascii.c file (Don't touch!) */
  {"float",    writeFp  }, 
  {"vector",   writeVec },
  {"int",      writeInt },
  {"double",   writeFp  },
  {"S",        writeS   },
  {"g",        writeg   },
  {"MB",       writeMB  },
  /* ADDED 24/10/2000: Parallel Tempering */
  {"PE",       writePE  },
  {"PEij",   writePEij},
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
