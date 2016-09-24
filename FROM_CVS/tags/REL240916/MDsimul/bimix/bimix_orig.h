/* ================= >>>  SIMULATION DEPENDENT HEADER FILE <<< ================
   This contains file contains the follow simulation dependent data types 
   and instances:
    - instance of some structures ( params, filenames and singlePar) 
    - number macros ( NUM_MISURE, BAK_STEPS, STA_STEPS ) 
    - list macros ( that is coordinates list: SAVE_LIST, ALLOC_LIST,
                    DECL_LIST ) */

/* ================== >>> PROGRAM DEFINES(CUSTOMIZE!) <<< ===================*/

#define BIMIX

#define MDSIMUL "/home/demichel/shared/simul/mdsimul"
/* this is the executable, you must change this to your directory */

#define XTERM   "/usr/X11R6/bin/nxterm"

#define MD_HOME "/home/demichel"
/* directory to store temporary files */ 
#define MD_SIMDAT MD_HOME "/simdat"

#ifdef MD_LOADMESH
#define MD_MESHDIR MD_HOME "/simdat"
#endif

#define MD_HD_TMP MD_SIMDAT "/mdtmp/"
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

#define NA 2 /* types of atoms in the mixture */

#define MAXPAR 1500      /* maximum number of simulated particles */

#define NUM_PAR 500   /* Number of particles for the simulation */
#define NUMK 149    /* number of k-points in which we must  calculate the 
		       structure factor */ 
#define MAXBIN 1000  /* Number of radius in which calculate the radial 
		       distribution function */
#define NUMV 150
#define GSRMAX 2.0 /* Maximum r for which we calculate the van Hove function*/
#define GSPOINT 30  /* Pont for which we plot the van Hove function */

#define FSPOINT 30    /* Point for which we plot the self-part of the 
			 intermediate scattering function */
#define FSKMAX 15.0 /* see the maximum k for the static structure factor */
/* ========================================================================= */

/* ======= >>> YOU MUST CHANGE THE SAVE_LIST AND ALLOC_LIST <<< ========== */

/* Follows a list of pointer to COORD_TYPE, you must think these objects like
   array of COORD_TYPE, in this implementation you can allocate whatever
   COORD_TYPE array you want, each sizeof(COORD_TYPE)*<particles number> bytes
   long 
*/
#define SAVE_LISTA rx[0], ry[0], rz[0],\
                  vx[0], vy[0], vz[0],\
                  Fx[0], Fy[0], Fz[0],\
                  vxo1[0], vyo1[0], vzo1[0],\
                  vxo2[0], vyo2[0], vzo2[0]
#define SAVE_LISTB rx[1], ry[1], rz[1],\
                   vx[1], vy[1], vz[1],\
                   Fx[1], Fy[1], Fz[1],\
                   vxo1[1], vyo1[1], vzo1[1],\
                   vxo2[1], vyo2[1], vzo2[1]

#undef  EXT_SLST
#define EXT_SLST  &s, &s1, &s2, &Vol, &Vol1, &Vol2, &Vol1o1, &s1o1, &Vol1o2,\
                  &s1o2

/* Reduced list of variables to save on xva file (tape file) */ 
#define XVA_LISTA rx[0], ry[0], rz[0],\
                  vx[0], vy[0], vz[0]

#define XVA_LISTB rx[1], ry[1], rz[1],\
                  vx[1], vy[1], vz[1]          

#define XVA_NUM 6 /* this is the number of vars in XVA_LISTA(B) 
		     <--- SET THIS!!!*/

#define XVA_ALSTA &rx[0], &ry[0], &rz[0],\
                 &vx[0], &vy[0], &vz[0]

#define XVA_ALSTB &rx[1], &ry[1], &rz[1],\
                  &vx[1], &vy[1], &vz[1]
          
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
#if defined(TAPPING) || defined(BASINS)
#define ALLOC_LISTA  &rx[0], &ry[0], &rz[0],\
                    &vx[0], &vy[0], &vz[0],\
                    &Fx[0], &Fy[0], &Fz[0],\
                    &vxt[0], &vyt[0], &vzt[0],\
                    &FxLJ[0], &FyLJ[0], &FzLJ[0],\
                    &vxo1[0], &vyo1[0], &vzo1[0],\
                    &vxo2[0], &vyo2[0], &vzo2[0],\
                    &vxg[0],  &vyg[0],  &vzg[0],\
                    &vxt2[0],  &vyt2[0],  &vzt2[0],\
		    &xicomx[0], &xicomy[0], &xicomz[0],\
		    &pcomx[0], &pcomy[0], &pcomz[0],\
		    &xix[0], &xiy[0], &xiz[0],\
		    &Gx[0], &Gy[0], &Gz[0],\
		    &Hx[0], &Hy[0], &Hz[0]

#define ALLOC_LISTB &rx[1], &ry[1], &rz[1],\
                    &vx[1], &vy[1], &vz[1],\
                    &Fx[1], &Fy[1], &Fz[1],\
                    &vxt[1], &vyt[1], &vzt[1],\
                    &FxLJ[1], &FyLJ[1], &FzLJ[1],\
                    &vxo1[1], &vyo1[1], &vzo1[1],\
                    &vxo2[1], &vyo2[1], &vzo2[1],\
                    &vxg[1],  &vyg[1],  &vzg[1],\
                    &vxt2[1],  &vyt2[1],  &vzt2[1],\
		    &xicomx[1], &xicomy[1], &xicomz[1],\
		    &pcomx[1], &pcomy[1], &pcomz[1],\
		    &xix[1], &xiy[1], &xiz[1],\
		    &Gx[1], &Gy[1], &Gz[1],\
		    &Hx[1], &Hy[1], &Hz[1]

#else
#define ALLOC_LISTA  &rx[0], &ry[0], &rz[0],\
                    &vx[0], &vy[0], &vz[0],\
                    &Fx[0], &Fy[0], &Fz[0],\
                    &vxt[0], &vyt[0], &vzt[0],\
                    &FxLJ[0], &FyLJ[0], &FzLJ[0],\
                    &vxo1[0], &vyo1[0], &vzo1[0],\
                    &vxo2[0], &vyo2[0], &vzo2[0],\
                    &vxg[0],  &vyg[0],  &vzg[0],\
                    &vxt2[0],  &vyt2[0],  &vzt2[0]
#define ALLOC_LISTB &rx[1], &ry[1], &rz[1],\
                    &vx[1], &vy[1], &vz[1],\
                    &Fx[1], &Fy[1], &Fz[1],\
                    &vxt[1], &vyt[1], &vzt[1],\
                    &FxLJ[1], &FyLJ[1], &FzLJ[1],\
                    &vxo1[1], &vyo1[1], &vzo1[1],\
                    &vxo2[1], &vyo2[1], &vzo2[1],\
                    &vxg[1],  &vyg[1],  &vzg[1],\
                    &vxt2[1],  &vyt2[1],  &vzt2[1]
#endif              
/* this is used to declare the particle variables ( see below ) 
   NOTE: rx[0][2] means the x-coordinate of the first atoms in the second 
   molecules (particle).
   To remember this note that when we must take a coordinate, going 
   from right to left, first we choose the molucule, then the atom and 
   finally the coordinate, for example consider the position: 
   coordinate(rx, ry, rz) <- atom <- molecule*/
#if defined(TAPPING) || defined(BASINS)  
#define DECL_LIST   *rx[NA], *ry[NA], *rz[NA],\
		    *vx[NA], *vy[NA], *vz[NA],\
                    *vxt[NA], *vyt[NA], *vzt[NA],\
                    *Fx[NA], *Fy[NA], *Fz[NA],\
                    *FxLJ[NA], *FyLJ[NA], *FzLJ[NA],\
                    *vxo1[NA], *vyo1[NA], *vzo1[NA],\
                    *vxg[NA], *vyg[NA], *vzg[NA],\
                    *vxt2[NA], *vyt2[NA], *vzt2[NA],\
                    *vxo2[NA], *vyo2[NA], *vzo2[NA],\
		    *xicomx[NA], *xicomy[NA], *xicomz[NA],\
		    *pcomx[NA], *pcomy[NA], *pcomz[NA],\
		    *xix[NA], *xiy[NA], *xiz[NA],\
		    *Gx[NA], *Gy[NA], *Gz[NA],\
		    *Hx[NA], *Hy[NA], *Hz[NA]
	    
#else
#define DECL_LIST   *rx[NA], *ry[NA], *rz[NA],\
                    *vx[NA], *vy[NA], *vz[NA],\
                    *vxt[NA], *vyt[NA], *vzt[NA],\
                    *Fx[NA], *Fy[NA], *Fz[NA],\
                    *FxLJ[NA], *FyLJ[NA], *FzLJ[NA],\
                    *vxo1[NA], *vyo1[NA], *vzo1[NA],\
                    *vxg[NA], *vyg[NA], *vzg[NA],\
                    *vxt2[NA], *vyt2[NA], *vzt2[NA],\
                    *vxo2[NA], *vyo2[NA], *vzo2[NA]
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

  int endFormat; /* 0 = binary 1 = ascii 2 = both */

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

  COORD_TYPE PxyArr[5];
  COORD_TYPE PyzArr[5];
  COORD_TYPE PzxArr[5];

  double rxi[NA][MAXPAR];
  double ryi[NA][MAXPAR];
  double rzi[NA][MAXPAR];
  /* Accumulator for the radial distribution function */
  int hist[MAXBIN];
  
  /* Accumulator for the static structure factor */ 
  COORD_TYPE sumS[NUMK];

  /* Accumulator for the MAXWELL-BOLTZMANN distribution */ 
  int histMB[NUMV];

  COORD_TYPE sumTemp;
  COORD_TYPE sumPress;
  
  /* Accumulators for calculating the rotational diffusion coefficent */
  int savedXva; 
  COORD_TYPE Q;
  COORD_TYPE W;
  int Nose; /* if Nose=1 use Nose method, otherwise not */
  int sResetSteps; /* Steps at which reset s to 1 */
  int CMreset;
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
  int HNBOX;
  int avngS;
  int avnggr;
  int avngEta;
  int avngTemp;
  int avngPress;
  int avngMB;
#ifdef TAPPING
  int tapping;
  int taptau;
  double taptol;
#endif  
  /* ADD 13/4/2000 
     Logarithmic saving implemented for xva file */
  int xvaSaveMode;/* 0 = linear 1 = semilog 2 = bilog (not impl. yet) */
  int bakSaveMode;/* save mode for ascii backup */
  char tmpPath[NAME_LENGTH];
  char misPath[NAME_LENGTH];

  double base;    /* We save at base^^NN step */
  int NN;         /* Logatithmic block length */
  double fstps;         /* There are KK block each base^NN long */
  char nRun[132];
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
  int parnum[NA];        	/* particles number of
				 species 0 (A) and 1(B)*/
  COORD_TYPE steplength;		/* temporal step length */
  int totStep;	/* temporal step number that simulation 
				   must do */
  int curStep;	/* current step of simulation */
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  COORD_TYPE d;                 /* distance between atoms */
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m[NA];             /* atoms masses */
  COORD_TYPE rcut;              /* cutoff for the pair potential */ 
  COORD_TYPE sigab[NA][NA];     /* pair potential length parameters */
  COORD_TYPE epsab[NA][NA];     /* pair potential energy parameters */
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
#ifdef SOFT_SPHERE
  int PP; 
#endif

  COORD_TYPE tol;               /* Tolerance of the shake algoritm used 
				   by RATTLE */
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
  {"Q",            &OS(Q),                          1,              1 , "%.6G"},
  {"W",            &OS(W),                          1,              1, "%.6G"},
  {"Nose",         &OS(Nose),                       1,   1,  "%d"},
  //  {"zeroall_s",    &OS(zeroall_s),                  1,   1,  "%d"},
  {"sResetSteps",  &OS(sResetSteps),                1,   1,  "%d"},
  {"Cmreset",      &OS(CMreset),                    1,   1,  "%d"},
  {"nebrTabFac",   &OS(nebrTabFac),                 1,   1,   "%d"},
  {"rNebrShell",   &OS(rNebrShell),                 1,   1, "%.6G"},
  {"noLinkedList", &OS(noLinkedList),               1,   1,  "%d"},
  {"avs",          &OS(avs),                        1,   1, "%.8G"},
  {"tols",         &OS(tols),                       1,   1, "%.8G"},
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
  {"fstps",        &OS(fstps),                      1,  1,   "%.15G"},
  {"nRun",         OS(nRun),                       1,  1,   "%s"},
#ifdef TAPPING
  {"taptau",       &OS(taptau),                    1,  1,   "%d"},
  {"tapping",      &OS(tapping),                   1,  1,   "%d"},
  {"taptol",       &OS(taptol),                    1,  1,   "%.20G"},
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
  {"parnum",            OP(parnum),                       2,   1, "%d"},
  {"steplength",        &OP(steplength),                  1,   1, "%.6G"},
  {"totStep",           &OP(totStep),                     1,   1,  "%d"},
  {"curStep",           &OP(curStep),                     1,   1,  "%d"},
  {"P",                 &OP(P),                           1,   1, "%.6G"},
  {"T",                 &OP(T),                           1,   1, "%.6G"},
  {"m",                 OP(m),                            2,   1, "%.6G"},
#ifdef SOFT_SPHERE
  {"PP",                &OP(PP),                          1,   1,  "%d"},
#endif
  {"rcut",              &OP(rcut),                        1,   1, "%.8G"},
  {"sigma",             OP(sigab),                        2,   2, "%.8G"},
  {"epsilon",           OP(epsab),                        2,   2, "%.8G"},
  {"equilibrat",         &OP(equilibrat),                  1,   1,   "%d"},
  {"M",                 &OP(M),                           1,   1,  "%d"},
  {"tol",               &OP(tol),                         1,   1, "%.15G"},
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
  {"parnumA" ,    &Oparams.parnum[0],             INT},
  {"parnumB" ,    &Oparams.parnum[1],             INT},
  {"steplength", &Oparams.steplength,         CT},
  {"stepnum",    &Oparams.totStep,            INT},
  {"inifile" ,   &OprogStatus.inifile,        STR},
  {"endfile" ,   &OprogStatus.endfile,        STR},
  {"endFormat",  &OprogStatus.endFormat,      INT},
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
  {"VolSteps",   &OprogStatus.measSteps[7],   INT},
  {"VolCalc",    &OprogStatus.measCalc[7],    INT},
  {"VolName",    &OprogStatus.dataFiles[7],   STR},
  {"sSteps",     &OprogStatus.measSteps[8],  INT},
  {"sName",      &OprogStatus.dataFiles[8],  STR},
  {"sCalc",      &OprogStatus.measCalc[8],   INT},
  {"PtensSteps", &OprogStatus.measSteps[9],  INT},
  {"PtensCalc",  &OprogStatus.measCalc[9],   INT},   
  {"PtensName",  &OprogStatus.dataFiles[9],  STR},
  {"DQtensSteps",&OprogStatus.measSteps[10],  INT},
  {"DQtensCalc", &OprogStatus.measCalc[10],   INT},
  {"DQtensName", &OprogStatus.dataFiles[10],  STR},
  {"VcSteps",&OprogStatus.measSteps[11],   INT},
  {"VcCalc", &OprogStatus.measCalc[11],    INT},
  {"VcName", &OprogStatus.dataFiles[11],   STR},
  {"PexSteps",  &OprogStatus.measSteps[12], INT},
  {"PexCalc",    &OprogStatus.measCalc[12], INT},
  {"PexName",    &OprogStatus.dataFiles[12], STR},
  {"PexInit",    &OprogStatus.initStep[12],  INT}, 
  {"Q",          &OprogStatus.Q,              CT},
  {"Nose",       &OprogStatus.Nose,           INT},
  {"sResetSteps",&OprogStatus.sResetSteps,    INT},
  {"CMreset",    &OprogStatus.CMreset,        INT},
  {"noLinkedList",&OprogStatus.noLinkedList,  INT},
  {"rNebrShell", &OprogStatus.rNebrShell,     CT},
  {"nebrTabFac", &OprogStatus.nebrTabFac,     INT},
  {"W",          &OprogStatus.W,              CT},
  {"P",          &Oparams.P,                  CT},
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
  {"sigaa",      &Oparams.sigab[0][0],      CT},
  {"sigbb",      &Oparams.sigab[1][1],      CT},
  {"sigab",      &Oparams.sigab[0][1],      CT},
  {"mass0",      &Oparams.m[0],             CT},
  {"mass1",      &Oparams.m[1],             CT},
  {"cellNum",    &Oparams.M,                INT},
  {"epsaa",      &Oparams.epsab[0][0],      CT},
  {"epsbb",      &Oparams.epsab[1][1],      CT},
  {"epsab",      &Oparams.epsab[0][1],      CT},
#ifdef SOFT_SPHERE
  {"PP",         &Oparams.PP,               INT},
#endif 
#ifdef TAPPING
  {"tapping",    &OprogStatus.tapping,      INT},
  {"taptau",     &OprogStatus.taptau,       INT},
  {"taptol",     &OprogStatus.taptol,       CT},
#endif
  {"rcut",       &Oparams.rcut,             CT},
  {"atomsDist",  &Oparams.d,                CT},
  {"equilibrat", &Oparams.equilibrat,       INT},
  {"temperat",   &Oparams.T,                CT},
  {"tol",        &Oparams.tol,              CT},
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
#define REND 5.65
#define HNBOX 0
#define NUMK2AV 150 /* Number of k on the same modulus over which we must 
		       perform a mean */

#ifdef MAIN
COORD_TYPE Vc, E, Dtrans[NA], temp, S[NUMK], dummy, eta, Drot, gr[MAXBIN], invs, press, Pex;
COORD_TYPE Ptens[3], DQtens[3], DrSqTot[NA];
int MB[NUMV];
#else 
extern COORD_TYPE Vc, E, Dtrans[NA], temp, S[NUMK], dummy, eta, Drot, gr[MAXBIN], invs, press, Pex;
extern COORD_TYPE Ptens[3], DQtens[3], DrSqTot[NA];
extern int MB[NUMV];
#endif

/* ============= >>> PUT HERE MEASURING FUNCTION PROTOTYPES <<< ============*/
#if !defined(BASINS)
void energy(void);
void transDiff(void);
void temperat(void);
void structFacts(void);
void maxwell(void);
void rotDiff(void);
void radDens(void);
void Ptensor(void);
void DQtensor(void);
void press_ex(void);

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
  {&Dtrans[0],  sizeof(COORD_TYPE),   transDiff   },
  {&temp,   sizeof(COORD_TYPE),       temperat    },
  {gr,      sizeof(COORD_TYPE)*MAXBIN,radDens     },
  {MB,      sizeof(int)*NUMV,         maxwell     },
  {S,       sizeof(COORD_TYPE)*NUMK,  structFacts },
  {&press,  sizeof(COORD_TYPE),       NULL        },
  {&Vol,    sizeof(COORD_TYPE),       NULL        },
  {&s,      sizeof(COORD_TYPE),       NULL        },
  {Ptens,   sizeof(COORD_TYPE)*3,     Ptensor     },
  {DQtens,  sizeof(COORD_TYPE)*3,     DQtensor    },
  {&Vc,     sizeof(COORD_TYPE),       NULL        },
  {&Pex,    sizeof(double),           press_ex    },
  /* ======================================================================= */
  {NULL, 0, NULL}                   /* End of list don't touch this */
};
#else 
extern struct measure Omeasure[];
#endif
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
  int parnum[NA];
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
#if !defined(BASINS)
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
#endif
/* ====================== >>> move() AND SUBROUTINES <<< ====================*/

void maps(void);

/* ==========================================================================*/
