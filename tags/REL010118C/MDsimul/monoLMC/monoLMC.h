/* ================= >>>  SIMULATION DEPENDENT HEADER FILE <<< ================
   This contains file contains the follow simulation dependent data types 
   and instances:
    - instance of some structures ( params, filenames and singlePar) 
    - number macros ( NUM_MISURE, BAK_STEPS, STA_STEPS ) 
    - list macros ( that is coordinates list: SAVE_LIST, ALLOC_LIST,
                    DECL_LIST ) */

/* ================== >>> PROGRAM DEFINES(CUSTOMIZE!) <<< ===================*/

#undef PRESSURE
#define MLMC

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

#define NA 1 /* types of atoms in the mixture */

#define MAXPAR 1500      /* maximum number of simulated particles */

#define NUM_PAR 500   /* Number of particles for the simulation */
#define NUMK 500    /* number of k-points in which we must  calculate the 
		       structure factor */ 
#define MAXBIN 150  /* Number of radius in which calculate the radial 
		       distribution function */

#define NUMV 150
/* ========================================================================= */

/* ======= >>> YOU MUST CHANGE THE SAVE_LIST AND ALLOC_LIST <<< ========== */

/* Follows a list of pointer to COORD_TYPE, you must think these objects like
   array of COORD_TYPE, in this implementation you can allocate whatever
   COORD_TYPE array you want, each sizeof(COORD_TYPE)*<particles number> bytes
   long 
*/
#define SAVE_LIST rx, ry, rz
            

#undef  EXT_SLST
#define EXT_SLST  &Vol, &Vol1, &Vol2, &Vol1o1, &Vol1o2

/* Reduced list of variables to save on xva file (tape file) */ 
#define XVA_LIST rx, ry, rz


#define XVA_NUM 3 /* this is the number of vars in XVA_LISTA(B) 
		     <--- SET THIS!!! */

#define XVA_ALST &rx, &ry, &rz
          
#define XVA_DLST *rx, *ry, *rz
           
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
#define ALLOC_LIST &rx, &ry, &rz

              
/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#define DECL_LIST   *rx, *ry, *rz            				   

#undef EXT_DLST
#define EXT_DLST Vol, Vol1, Vol2, Vol1o1, Vol1o2

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
  int endFormat;

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

  int rxi[MAXPAR];
  int ryi[MAXPAR];
  int rzi[MAXPAR];
  
  /* Accumulator for the static structure factor */ 
  COORD_TYPE sumS[NUMK];

  /* Accumulators for calculating the rotational diffusion coefficent */
  int savedXva; 

  int mosseAccettate;
  int mosseTentate;
  int rateCheck;
  int CMreset;
  int nebrTabFac;                /* How much storage sould be provided for 
                                   the neighbour list (see Rapaport pag.53
				   for details )*/
  COORD_TYPE rNebrShell;         /* = Dr see Rapaport pag. 53 */

  /* Accumulator for the MAXWELL-BOLTZMANN distribution */ 
  int hist[MAXBIN];

  COORD_TYPE avVol;
  COORD_TYPE avVol1;
  COORD_TYPE tolVol;
  COORD_TYPE tolVol1;
  int HNBOX;
  int avngS;
  int avnggr;
  int avngTemp;
  int avngPress;
  /* ADD 13/4/2000 
     Logarithmic saving implemented for xva file */
  int xvaSaveMode;/* 0 = linear 1 = semilog 2 = bilog (not impl. yet) */
  int bakSaveMode;/* save mode for ascii backup */
  char tmpPath[NAME_LENGTH];
  char misPath[NAME_LENGTH];

  double alpha;
  double S0;
  double kmax;
  double Dk;
  int AntiCry;
  int chiamate; /* ogni quanti passi bisogna aggiornare le S(k)-neighbour list */
  double DS0;/* questa serve per costruire le neighbour list per S(k) */

  double base;    /* We save at base^^NN step */
  int NN;         /* Logatithmic block length */
  int fstps;         /* There are KK block each base^NN long */
  char nRun[32];
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
  int parnum;        	/* particles number of
				 species 0 (A) and 1(B)*/
  COORD_TYPE steplength;		/* temporal step length */
  int totStep;	/* temporal step number that simulation 
				   must do */
  int curStep;	/* current step of simulation */
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m;             /* atoms masses */
  COORD_TYPE rcut;              /* cutoff for the pair potential */ 
  COORD_TYPE sigma;     /* pair potential length parameters */
  COORD_TYPE epsilon;     /* pair potential energy parameters */
  int equilibrat;               /* != 0 if equilibrating */
  
  /* questi li ho messi per poter leggere i .cor file prodotti 
     dalla simulazione MC */
  int lattice_M;                /* Numero di celle su di un asse */
  double lattice_a;             /* passo del reticolo */
  int num_neigh;                
  /* Ogni volta che tenta una mossa nella MC deve scegliere dei vicini
     questo parametro indica quanti bisogna considerarne lungo un asse 
     (1,2,..) */
  
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

/* MAIN is a macro defined only in mdsimul.c, so there we put the instance
of the object, and elsewhere we put only extern declarations */
#ifdef MAIN
int DECL_LIST;
COORD_TYPE EXT_DLST;
#else
extern int DECL_LIST;
extern COORD_TYPE EXT_DLST;
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
  {"initStep",     OS(initStep),          NUM_MISURE,               1, "%d"},
  {"hist",         OS(hist),                  MAXBIN,               1, "%d"},
  {"sumS",         OS(sumS),                    NUMK,               1, "%.6G"},
  //mancano le rxi per i coeff. di diffusione!!!!
  {"mosseAccettate",&OS(mosseAccettate),            1,   1,  "%d"},
  {"mosseTentate",  &OS(mosseTentate),              1,   1,  "%d"},
  {"rateCheck",     &OS(rateCheck),                 1,   1,  "%d"},
  {"Cmreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.6G"},
  {"HNBOX",        &OS(HNBOX),                      1,   1,  "%d"},
  {"avngTemp",     &OS(avngTemp),                   1,   1,  "%d"},
  {"avngPress",    &OS(avngPress),                  1,   1,  "%d"},
  {"avnggr",       &OS(avnggr),                      1,  1,  "%d"},
  {"avngS",        &OS(avngS),                      1,   1,  "%d"},
  //  {"avgMedia",     &OS(avgMedia),                   1,     1, "%d"},
  {"xvaSavedMode", &OS(xvaSaveMode),                1,  1,    "%d"},
  {"bakSavedMode", &OS(bakSaveMode),                1,  1,    "%d"},
  {"tmpPath",      OS(tmpPath),                     1,  NAME_LENGTH, "%s"},
  {"misPath",      OS(misPath),                     1,  NAME_LENGTH, "%s"},
  {"base",         &OS(base),                       1,  1, "%.6G"},
  {"NN",           &OS(NN),                         1,  1,   "%d"},
  {"fstps",        &OS(fstps),                      1,  1,   "%d"},
  {"nRun",         OS(nRun),                       1,  32,   "%s"},
  {"alpha",        &OS(alpha),                      1,   1, "%.6G"},
  {"DS0",           &OS(DS0),                         1,   1, "%.6G"},
  {"S0",           &OS(S0),                         1,   1, "%.6G"},
  {"kmax",         &OS(kmax),                       1,   1, "%.6G"},
  {"Dk",           &OS(Dk),                         1,   1, "%.6G"},
  {"chiamate",     &OS(chiamate),                   1,   1,  "%d"},
  {"", NULL, 0, 0, ""}
};
#else
extern struct pascii opro_ascii[];
#endif

#ifdef MAIN
/* ========================= >>> opar_ascii <<< =========================== */
struct pascii opar_ascii[]=
{
  {"parnum",            &OP(parnum),                      1,   1, "%d"},
  {"steplenght",        &OP(steplength),                  1,   1, "%.6G"},
  {"totStep",           &OP(totStep),                     1,   1,  "%d"},
  {"curStep",           &OP(curStep),                     1,   1,  "%d"},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
  {"m",                 &OP(m),                           1,   1, "%.6G"},
  {"rcut",              &OP(rcut),                        1,   1, "%.8G"},
  {"sigma",             &OP(sigma),                       1,   1, "%.8G"},
  {"epsilon",           &OP(epsilon),                     1,   1, "%.8G"},
  {"equlibrat",         &OP(equilibrat),                  1,   1,   "%d"},
  {"lattice_M",         &OP(lattice_M),                   1,   1, "%d"},
  {"lattice_a",         &OP(lattice_a),                   1,   1,  "%.15G"},
  {"", NULL, 0, 0, ""}
};
#else
extern struct pascii opar_ascii[];
#endif



/* ------------------------ PARTICLES COORDINATES -------------------- */


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
  {"bakStepsAscii",&OprogStatus.bakStepsAscii, INT},
  {"staSteps",   &OprogStatus.staSteps,       INT},
  {"xvaSteps",   &OprogStatus.xvaSteps,       INT},  
  {"energySteps",&OprogStatus.measSteps[0],   INT},
  {"energyCalc", &OprogStatus.measCalc[0],    INT},
  {"energyName", &OprogStatus.dataFiles[0],   STR},
  {"DtrSteps",   &OprogStatus.measSteps[1],   INT}, /* steps between measure
						       savings */
  {"DtrCalc",    &OprogStatus.measCalc[1],    INT}, /* steps between measure
						       calculation */
  {"DtrName",    &OprogStatus.dataFiles[1],   STR},
  {"grSteps",    &OprogStatus.measSteps[2],   INT},
  {"grInitStep", &OprogStatus.initStep[2],    INT},
  {"grCalc",     &OprogStatus.measCalc[2],    INT},
  {"grName",     &OprogStatus.dataFiles[2],   STR},

  {"SSteps",     &OprogStatus.measSteps[3],   INT},
  {"SCalc",      &OprogStatus.measCalc[3],    INT},
  {"SName",      &OprogStatus.dataFiles[3],   STR},
  {"VcSteps",    &OprogStatus.measSteps[4],   INT},
  {"VcCalc",     &OprogStatus.measCalc[4],    INT},
  {"VcName",     &OprogStatus.dataFiles[4],   STR},
  {"rateCheck",  &OprogStatus.rateCheck,      INT},
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
  {"Vol",        &Vol,                        CT},
  {"avVol",      &OprogStatus.avVol,          CT},
  {"tolVol",     &OprogStatus.tolVol,         CT},
  {"tolVol1",    &OprogStatus.tolVol1,        CT},
  {"avVol1",     &OprogStatus.avVol1,         CT},
  {"avngS",      &OprogStatus.avngS,          INT},
  {"avnggr",     &OprogStatus.avnggr,         INT},
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */
  {"DS0",        &OprogStatus.DS0,            CT},
  {"S0",         &OprogStatus.S0,            CT},
  {"AntiCry",    &OprogStatus.AntiCry,      INT},
  {"chiamate",   &OprogStatus.chiamate,      INT},
  {"kmax",       &OprogStatus.kmax,          CT},
  {"alpha",      &OprogStatus.alpha,         CT},
  {"Dk",         &OprogStatus.Dk,            CT},
  {"sigma",      &Oparams.sigma,           CT},
  {"mass",       &Oparams.m,             CT},
  {"epsilon",    &Oparams.epsilon,          CT},
  {"rcut",       &Oparams.rcut,             CT},
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"temperat",   &Oparams.T,                CT},
  {"lattice_M",  &Oparams.lattice_M,        INT},
  /* parametri per scegliere il tipo di salvataggio del file xva
     (lineare o semilog) */
  {"xvaSaveMode",&OprogStatus.xvaSaveMode,  INT},
  {"bakSaveMode", &OprogStatus.bakSaveMode, INT},
  {"tmpPath",    OprogStatus.tmpPath,       STR},
  {"misPath",    OprogStatus.misPath,       STR},
  {"base",       &OprogStatus.base,         CT},
  {"NN",         &OprogStatus.NN,           INT},
  {"nRun",       OprogStatus.nRun,         STR},
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
#define REND 3.7
#define HNBOX 0
#define NUMK2AV 150 /* Number of k on the same modulus over which we must 
		       perform a mean */

#ifdef MAIN
COORD_TYPE V, W, VLJ, WLJ, VAC, WAC, V, E, Dtrans, temp, S[NUMK], dummy, gr[MAXBIN], press;
#else 
extern COORD_TYPE V, W, WLJ, VLJ, VAC, WAC, V, E, Dtrans, temp, S[NUMK], dummy, gr[MAXBIN], press;
#endif

/* ============= >>> PUT HERE MEASURING FUNCTION PROTOTYPES <<< ============*/

void energy(void);
void transDiff(void);
void structFacts(void);
void radDens(void);

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
  {gr,      sizeof(COORD_TYPE)*MAXBIN,radDens     },
  {S,       sizeof(COORD_TYPE)*NUMK,  structFacts },
  {&V,     sizeof(COORD_TYPE),       NULL        },
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
  int parnum;
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
  /* ADDED 14/09/2000 */
  int mode;
  int NN;
  double T;
  double Vol;
  double base;
  double dt;     /* dt of each steps */
  //int N;
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
