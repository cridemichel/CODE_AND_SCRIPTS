#undef MAIN
#include<mdsimul.h>

/* WARNING: Extern variables must be of the same type of the real definition */

/* ====================== >>> Descriptors and streams <<< ================== */

extern FILE* output;

/* ======================= >>> Simulation vars <<< ==========================*/

extern unsigned char BAK, STA;       /* global switches for HD(0/1) */
extern unsigned char BAKT;           /* global switch for tape (0/1) */

int logBlock;

extern int pid;
extern int SEGSIZE;
/* ======================== >>> Structures <<< ============================= */

extern struct simStat OsimStat;
extern char msgStrA[MSG_LEN], msgStrB[MSG_LEN], msgStrC[MSG_LEN];
extern char paramFile[NAME_LENGTH];
extern char TXT[MSG_LEN], TXTA[10][MSG_LEN];

extern int my_rank, printerRank, numOfProcs;
/* ADDED 23/04/2001 */
int iniCorFormat = 0, iniBakFormat = 0; /* 0 = binary 1 = ascii */
char iniBakFile[NAME_LENGTH];

extern void AllocCoord(int PTM, int size, COORD_TYPE*** pointer, ...);
extern void readBakAscii(char* );
extern void readCorAscii(char* );

/* ============================= >>> Continue <<< ===========================*/
void Continue(void)
{
  /* Load restore file using an algorithm explained in the TECH_INFO
     file */

  int i;
   /* Opens the restore file that contains all informations
     to restart the simulation (see. TECH_INFO file for details).
     At least one file is not corrupted so if it fails to open 
     one file, it tries the second, for explanations on how it decides
     the file to use for restoring see TECH_INFO.
  */

  newSim = 0; /* New = 0 means: not a new simulation, i.e. continuing a
		 previously interrupted one */
  /* <------------------------------------------------------- OPEN LOG FILE */
  openLog("a"); /* append: continue to write on previously created log file */

  /* ======================= >>> CHOOSE RESTORE FILE <<< ===================*/
 
  /* NOTE: the restore file on disk is the only  used to restart simulation ,
     in fact the restore files on Tape should be used manually copying them
     on temporary dir on disk.
     So if a failure occur and both restore files are currupted on this 
     you must restart the simulation manually */

  mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	"Trying to restart simulation",
	"Choosing restore file on Hard Disk...",
	NULL);

  if (iniBakFormat == 0)
    {
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    "Trying to restart simulation",
	    "Choosing restore file on Hard Disk...",
	    NULL);
      /* -1 means : 'both restore files corrupted' */
      if (chooseRestore( absTmpHD(BAK_FILE_NAME), &BAK, readBak) == -1)
	{
	  mdPrintfSpcWr(ALL, "Simulation aborted.");
	  mdPrintfSpcWr(ALL, "Take restore files from Tape to restart");
	  exit(-1);
	}
    }
  else
    {
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    "Trying to restart simulation",
	    "Reading ascii restore file on Hard Disk...",
	    NULL);
      readBakAscii(iniBakFile);
    }



  /* readBak reads the restore file loading structure and coordinates into
     the memory (it allocates also shared memory to do this) */
  
  if (OprogStatus.tapeTimes != 0) /* Tape used? */
    {
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    "Choosing restore file on Tape...",
	    NULL);

      /* If Y choose next restore file to save on Tape.
         -1 = both restore files corrupted */
      if (chooseRestore( absTmpTape(BAK_FILE_NAME), &BAKT, readFile)
	   == -1)
	 {
	   mdPrintfSpcWr(ALL, "Both restore files on Tape corrupted");
	 }
      /* readFile reads the file byte by byte, checking... */
    }
  /* choose the good restore file in the Hard Disk temporary directory for 
     simulation ( MD_HD_TMP, see mdsimdep.h ) */
  
  /* ======================= >>> CHOOSE STATUS FILE <<< ====================*/
  if (OprogStatus.staSteps != 0)
    chooseStatus( absTmpHD(STATUS_FILE_NAME), NULL);
  /* second arg not used actually !!!! */
 
  //printf("xvasteps: %d\n", OprogStatus.xvaSteps);

  /* ======================= >>> CHOOSE MEASURE FILE <<< ===================*/
  for (i=0; Omeasure[i].buf != NULL; ++i)
    {
      /* 26/1/1998 ADD: */
      if (OprogStatus.measSteps[i] == 0) continue;
	      
      sprintf(msgStrA, "Choosing measure file %s on Hard Disk...",
	      OprogStatus.dataFiles[i]);
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    msgStrA,
	    NULL);
      /* If the measure was saved choose the good file */
      chooseMeasure( absMisHD(OprogStatus.dataFiles[i]), readFile );
      /* the file passed as argument is with path (absolute) */

      if (OprogStatus.tapeTimes == 0) /* Tape used ? */
	{ 
	  continue; /* continue to loop over measure list, ignoring 
		       the rest */
	}
      
      sprintf(msgStrA, "Choosing measure file %s on Tape...",
	      OprogStatus.dataFiles[i]);
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    msgStrA,
	    NULL);

      /* If the Tape is used make the same for it, that is choose the good
	 file and copy it onto the corrupted one*/
      if (OprogStatus.measSteps[i] != 0)
	chooseMeasure( absMisTape(OprogStatus.dataFiles[i]), readFile );
      /* the file passed as argument is with path (absolute)*/
      
    }
  
  /* ================ >>> Initialize OmeasHead structure <<< =============== 
     ( probably this is not necessary because the measure file was yet 
       opened and the new measures are appended to it )*/
  for (i = 0; Omeasure[i].buf != NULL; ++i)    
    {
      /* in the mean while initialize OmeasHead[] array */
      OmeasHead[i].saveSteps = OprogStatus.measSteps[i];
      OmeasHead[i].size = Omeasure[i].size;
      OmeasHead[i].PTM = Oparams.PTM;
      OmeasHead[i].N = Oparams.parnum;
      OmeasHead[i].PEpoints = PE_POINTS;
      OmeasHead[i].ENmin = OprogStatus.ENmin;
      OmeasHead[i].ENmax = OprogStatus.ENmax;
      OmeasHead[i].dt = Oparams.steplength;
      OmeasHead[i].T = Oparams.T;
      OmeasHead[i].Vol = Vol;
    }     
  /* ============== >>> Initialize OxvaHead structure <<< ===================*/
  //printf("SEGSIZE: %d\n", SEGSIZE);

  OxvaHead.saveSteps = OprogStatus.xvaSteps;
  OxvaHead.size = SEGSIZE * XVA_NUM; 
  /* ADD 13/09/2000: */
  OxvaHead.mode = OprogStatus.xvaSaveMode;
  OxvaHead.NN = OprogStatus.NN;
  OxvaHead.base = OprogStatus.base;
  OxvaHead.dt = Oparams.steplength;
  OxvaHead.parnum = Oparams.parnum;
  OxvaHead.T = Oparams.T;
  OxvaHead.Vol = Vol;
  logBlock = pow(OprogStatus.base, (double) OprogStatus.NN);

  /* avendo letto il file di backup vol e' gia' inizializzata
     correttamente */

  /* Manca da inizializzare il volume!!! */

  /* ---------------- */
  /* XVA_NUM is the number of variables in the XVA_LIST list (see mdsimdep.h)*/
}

/* ============================== >>> Parsing <<< ===========================*/
void Parsing(char stringA[NAME_LENGTH], char stringB[NAME_LENGTH])
{
  int i;
  /* check if 'stra', read from the params file, is one of the strings in the 
     singlePar array, if N => ERROR */ 
  for  (i=0; OsinglePar[i].ptr != NULL; ++i) /* parname=NULL menas END */
    {
      if (!strcmp(OsinglePar[i].parName, stringA))
	{
	  switch (OsinglePar[i].type)
	    {
	    case STR : 
	      strcpy((char *) OsinglePar[i].ptr, stringB);
	      break;
	    case INT :
	      *(int*) OsinglePar[i].ptr = atoi(stringB);
	      break;
	    case CT :
	      ( *((COORD_TYPE*)OsinglePar[i].ptr) ) = 
		(COORD_TYPE) atof(stringB); 
	      break;
	    default: 
	      mdMsg(ALL, NOSYS, "Parsing", "ERROR", NULL,
		    "Not valid parameter type in singlePar array(in mdsimdep.h).",
		    NULL);
	      
	      exit(-1);
	    }
	  /* successfully read prameter */
	  return;
	}
    }
   /* no one paremeter defined in the singlePar array matches the read one =>
     ERROR */
  sprintf(msgStrA, "Paramater %s is not valid", stringA);
  mdMsg(ALL, NOSYS, "Parsing", "ERROR", NULL,
	msgStrA,
	NULL);
  exit(-1);
}

void parallel_tempering_init(void);

/* =========================== >>> scanFile <<< ============================*/
void scanFile(char* argom)
{
  /* DESCRIPTION:
     scan the parameter file determining setting variuos parameters, the 
     scan proceeds as follows:
     1) Determine the inital coordinate file searching for 'inifile' paramters 
     
     2) Read coordinate file (therefore setting Oparams structure ) 
     
     3) Read particles number, that is 'parnum' parameter    
   
     4) Allocate shared memory for coordinates 
     
     5) Read all the other parameters values 
  */
  
  char str1[255],str2[255];            /* used to scan parameters file */
  FILE* pfs;
  int cfd;
  int parNumYN = 0, initYN = 0;
  int PTMYN = 0;
  char line[1024];

  if ( (pfs = fopen(argom, "r")) == NULL )
    {
      /* legge i parametri della simulazione dal file indicato
         come argomento e mette i valori nella struttura 
         Params */
      sprintf(msgStrA, "Error during opening file %s for reading.", argom);
      mdMsg(ALL, errno, "Init", "ERROR", "fopen",
	    msgStrA,
	    NULL);
      exit(-1);
    }
  
  /* Get particles number */ 
  while (!feof(pfs))
    {
      /* si aspetta un sintassi del tipo <parametro>:<valore> 
	 GLI SPAZI VENGONO IGNORATI <---------------------- !!!! */
      fscanf(pfs, "%[^\n] ", line);
      
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n# ] ", str1, str2) < 2)
	continue;
      
      if ( !strcmp(str1, "parnum") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum = atoi(str2);
	} 
     
      if ( !strcmp(str1, "PTM") && 
	   strcmp(str2, "*") )
	{
	  PTMYN = 1; /* parnum parameter found */
	  Oparams.PTM = atoi(str2);
	} 
    
    }  
  
  rewind(pfs);
  /* Then it determines initial coordinates file */
  while (!feof(pfs))
    {
      fscanf(pfs, "%[^\n] ", line);

      //printf("line: %s\n", line);
      /* si aspetta un sintassi del tipo <parametro>:<valore> 
         GLI SPAZI VENGONO IGNORATI <---------------------- !!!! */
      //fscanf(pfs, " %[^: ] : %[^'\n' ]", str1, str2);
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n# ] ", str1, str2) < 2)
	continue;
    
      if (!strcmp(str1, "inifile"))
	{
	  strcpy(OprogStatus.inifile, str2);
	}
    }

  /* Se il campo inifile di OprogStatus contiene un asterisco (*) 
     allora il set di coordinate iniziale deve essere generato 
     in qualche modo ( ved. dopo o SIMUL_INFO ) */
  if (!strcmp(OprogStatus.inifile, "*"))
    {
      mdMsg(ALL, NOSYS, "Init", "NOTICE", NULL,
	    "Generating initial coordinates...",
	    NULL);
      /* allocate  coordinates needed by simulation (see. mdsimul_p) */
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
      
      /* SEGSIZE is the size in bytes of an array of coordinates */
      
      /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
	 all addresses of the coordinates declared in the simulaiton
	 (see that file) */
      printf ("Oparams.PTM:%d\n", Oparams.PTM);
      
      AllocCoord(Oparams.PTM, SEGSIZE, ALLOC_LIST,
		 NULL);
      initYN = 1; /* Yes, initialize coordinates after reading 
		     all parameters */
    }
  else
    {
      if (parNumYN || PTMYN) /* Is parnum parameter specified? */
	{
	  /* if you specify an initial file particles number is not needed, 
	     because the particle number is set in coordinate file */
	  mdMsg(ALL, NOSYS, "Init", "WARNING", NULL,
		"You specify an initial coordinates file, so you have not to give",
		"particles number or the number of systems, I ignore parnum parameter",
		NULL
		);
	}
      
      sprintf(msgStrA, "Using file %s for initial coordinates",
	      OprogStatus.inifile);
      mdMsg(ALL, NOSYS, "Init", "NOTICE", NULL,
	    msgStrA,
	    NULL);
      
      cfd = mdOpen(OprogStatus.inifile, "Init", 
		   NULL, EXIT | O_RDONLY );
      
      /* MODIFICA 24/04/01 */
      if (iniCorFormat == 0)
	{
	  
	  /* carica le coordinate dal file 'OprogStatus.inifile' */
	  readCoord(cfd);
	}
      else
	{
	  readCorAscii(OprogStatus.inifile);
	}
      
      close(cfd);
    }  

 /* go to the begin of file pfs */
  rewind(pfs);

  /* scanning the file for the other params */
  while (!feof(pfs))
    {
      /* The syntax must be <parameter>:<value> 
         if <value>  is a '*' then it use the default value or 
         the value loaded from the coordinates file inifile, if
         specified. 
         SPACES ARE IGNORED <---------------------- !!!! */
      fscanf(pfs, "%[^\n] ", line);
      
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n#] ", str1, str2) < 2)
	continue;
    
      //fscanf(pfs, " %[^: ] : %[^'\n' ]", str1, str2);
      /* analyzes the parameter just read 
         This function depends strongly upon simulation */
      Parsing(str1, str2);
    }
  fclose(pfs);	  /* Aggiungere controllo della chiusura <-----!!! */
  
  if (initYN) /* If not initial coords file was supplied initialize them */ 
    {
      parallel_tempering_init();
      initCoord();
    }
}

/* ============================= >>> Newsimul <<< ===========================*/
void Newsimul (char *argom)
{
  /* New simulation: read parameter file specified as argument and initialize
   simulation structure consequently */
  struct dirent *POdirent;
  DIR *rfdir;                          /* pointer to a directory */
  int car;
  int i,ii;                            /* counters */
  
  newSim = 1;
    
  /* <------------------------------------------------ INITIALIZE STRUCTURES */
  initBefore();		/* initialize Oparams with default settings */

  /* <---------------------------------------------------------DELETE TMP FILES
     delete all files created by a previously interrepted simulation,
     if there was one. (restore files deleted at exit of this procedure) */
  
  
  delDoubleHD(STATUS_FILE_NAME);
  
  delTmpHD(CHK_LOG);       /* delete log file for chk_status executable*/
  delTmpHD(MDS_LOG);       /* delete log file for simulation */
  delTmpHD(CF);
  /* <-----------------------------------------------------   OPEN LOG FILE  */
  openLog("w");          /* if exist, truncate it */

  /* check if exist mesure directory */
  existDir(MD_HD_MIS);

  /* check if exist mesure directory on tape*/
  existDir(MD_TAPE_MIS);

  /* check if exist temporary directory on tape*/
  existDir(MD_TAPE_TMP);

 
  /* check if exist temporary directory on HD */
  existDir(MD_HD_TMP);
  
  /* check that no restore files  exist, it would mean that there is 
     a previously interrupted simulation */
  rfdir = opendir(MD_HD_TMP);

  while ( (POdirent = readdir(rfdir)) )
    {
      if ( fileExists(POdirent, appSw(BAK_FILE_NAME, 0)) || 
	   fileExists(POdirent, appSw(BAK_FILE_NAME, 1)) 
	   )
	{
#ifdef MPI	  
	  if (my_rank == MPI_OUT_PROCESS)
	  {
#endif
	    //printf("my_rank = %d\n", my_rank);
	    printf("Init [WARNING]:");
	    printf("   At least one restore file exist!\n");
	    printf("   Probably a previuos interrupted simulation was interrupted.\n");
	    printf("   contiuing anyway? (loosing them)? (y/n) [n] ");
#ifdef MPI
	  }
#endif 
	  car = MPIgetchar();
	  if (car != 'y')
	    {
#ifdef MPI 
#ifdef PRINTER_PROC
	      MPI_Send("fine", 5, MPI_CHAR, 
			printerRank, 
			MD_MPI_END, MPI_COMM_WORLD);
#endif   
#endif
	      exit (-1);
	    }
	  break;
	}
    }
  closedir(rfdir);
  
  mdMsg(ALL, NOSYS, "Init", "NOTICE", NULL,
	"Beginning a new simulation...",
	NULL);
  /* legge i parametri della simulazione dal file passato 
     come argomento (ved. inizio) e setta i corrispondenti campi
     delle strutture Oparams e OprogStatus*/


  scanFile(argom); /* argom is the name of the parameter file */


  if (OprogStatus.tapeTimes > 1 && OprogStatus.xvaSaveMode != 0)
    {
#ifdef MPI
      if (my_rank == 0)
	{
#endif 
	  printf("tapeTimes must be equal to 1 (or 0) if xvaSaveMode is not 0!\n");
	  exit(-1);
#ifdef MPI
	}
#endif 
    }

  /* initialize BAK */
  BAK = 0;
  /* the first restore file saved is BAK_FILE_NAME0 , 
     it is not important which file is the first file when beginning 
     a new smiulation */
  STA = 0;  
  /* similarly it is not important which status file to use 
   at begin of a new simulation */
 
  BAKT = 0; 
  /* this is the switch for tape doubleBuf savings */

  delDoubleHD(BAK_FILE_NAME);  
  /* delete restore files on HD<------ DELETE RESTORE FILES*/
  
  delDoubleTape(BAK_FILE_NAME);
  /* delete restore files on Tape */
  
  /* 26/1/1998 ADD: Initialize the OprogStatus.measCalc[] array, so that if
     it is not specified in the parameter file the Calc steps then they are
     assumed equal to the Save Steps */
  for(i=0; i < NUM_MISURE; ++i)
    {
      if (OprogStatus.measCalc[i] == -1)
	{
	  OprogStatus.measCalc[i] = OprogStatus.measSteps[i];
	}
    }
  
  /* INITIALIZE OmeasHead structure */
  for (i = 0; Omeasure[i].buf != NULL; ++i)  
    {  
      /* This needs an explanation: Omeasure structure is a struct to contain 
	 various information on the measure asavings and calculations 
         ( pointer to measure datas on runtime, etc. ), instead
	 OmeasHead is structure that contains information only on the kind
	 of measure (size of each measure, steps between savings, etc.) .
	 OmeasHead is the analogous of Oparams, while Omeasure is the anlogous
	 of singlePar array.*/
      OmeasHead[i].saveSteps = OprogStatus.measSteps[i];
      OmeasHead[i].size = Omeasure[i].size;
      /* ADDED 28/09/2000 */
      OmeasHead[i].PTM = Oparams.PTM;
      OmeasHead[i].N = Oparams.parnum;
      OmeasHead[i].PEpoints = PE_POINTS;
      OmeasHead[i].ENmin = OprogStatus.ENmin;
      OmeasHead[i].ENmax = OprogStatus.ENmax;
      OmeasHead[i].dt = Oparams.steplength;
      OmeasHead[i].T = Oparams.T;
      OmeasHead[i].Vol = Vol;
      /* RIEMPIRE QUI ANCHE GLI ALTRI CAMPI VED. il file simdep.h!!!!*/
    }
  /* NOTICE: CREATE MEASURE FILE FOR FUTURE SAVINGS AND PUT THE HEADER 
     AT THE TOP OF EACH MEASURE FILE */
  for (i = 0; Omeasure[i].buf != NULL; ++i)
    {
      if (OprogStatus.measSteps[i] == 0) continue;
     // printf("Steps[%d]:%d\n", i, OprogStatus.measSteps[i]);
      for ( ii=0; ii<2; ++ii) /* ii=0 and ii=1 */
	{
	  creatWithHead( appSw(absMisHD(OprogStatus.dataFiles[i]), ii), 
			 "Init", "Unable to create measure file on Hard Disk",
			 EXIT, sizeof(struct measHead), &OmeasHead[i]);
	  if (OprogStatus.tapeTimes == 0) continue;      
	  /* don't use tape if zero ! */  
	  
	  /* ======== >>> Temporary change of OmeasHead struct <<< =========
	     the restore files on zip are saved every
	     'TAPE_TIMES * OmeasHead[].saveSteps' steps, so we must put
	     this value in the header of measures file on Tape.
	     oldSaveSteps is used to restore the original value later */
	  
	  pushSteps(&OmeasHead[i].saveSteps); 
	  /* store saveSteps field value and load in saveSteps the right 
	     value for tape measure savings */
	  creatWithHead( appSw(absMisTape(OprogStatus.dataFiles[i]), ii),
			 "Init", "Unable to create measure file on Hard Disk",
			 EXIT, sizeof(struct measHead), &OmeasHead[i]);
	  
	  /* =====>>> restore saveSteps field of OmeasHead struct <<< =======*/
	  popSteps(&OmeasHead[i].saveSteps);
	}
    }
  /* ============== >>> Initialize OxvaHead structure <<< ===================*/
  
  OxvaHead.saveSteps = OprogStatus.xvaSteps;
  OxvaHead.size = SEGSIZE * XVA_NUM; /* XVA_NUM is the number of variables 
					in the XVA_LIST list (see mdsimdep.h)*/
  /* ADD 13/09/2000: */
  OxvaHead.mode = OprogStatus.xvaSaveMode;
  OxvaHead.NN = OprogStatus.NN;
  OxvaHead.base = OprogStatus.base;
  OxvaHead.dt = Oparams.steplength;
  OxvaHead.parnum = Oparams.parnum;
  OxvaHead.T = Oparams.T;
  OxvaHead.Vol = Vol;
  OprogStatus.savedXva = -1;
  OprogStatus.fstps = 1;
  logBlock = pow(OprogStatus.base, (double) OprogStatus.NN);

  /* La prima volta che va in saveXva savedXva diventa 0 e cosi'
     deve essere !!! */
  /* CREATE TAPE FILE (XVA FILE) ON HARD DISK */
  
  /* 16/1/1998 CHG: absMis... to absXva */
  if ( OprogStatus.xvaSteps == 0 ) return;
  
  creatWithHead(absXvaHD(OprogStatus.xvafile), "Init", NULL, EXIT,
		sizeof(struct xvaHead), &OxvaHead);

  if ( OprogStatus.tapeTimes == 0) return;
  
  /* CREATE TAPE FILE ON TAPE IF NEEDED */

  pushSteps(&OxvaHead.saveSteps);
  creatWithHead(absXvaTape(OprogStatus.xvafile), "Init", NULL, EXIT, 
		sizeof(struct xvaHead), &OxvaHead);
  popSteps(&OxvaHead.saveSteps);
}
 
/* ========================= >>> args <<< =================================*/
int argsMd(int argc, char** argv)  
{
  int car;
  int i;

  /* You must specify -c or -f <filename> as argument */ 
  if ( argc == 1 ) 
	{
	  sprintf(TXTA[0], "What should I do? (No Arguments)\n");
	  sprintf(TXTA[1], "Try 'mdsimul -help' for valid arguments\n");
	  mdPrintfR0(STD, TXTA[0], TXTA[1]);
	  exit(-1);
	}
    if ((!strcmp(argv[1], "-continue_ascii")) ||
      (!strcmp (argv[1], "-ca")))
    {
      if (argc == 2)
	{
	  sprintf(TXT, "ERROR: -ca needs a ascii restart file name.\n");
	  mdPrintfR0(STD, TXT);
	  exit(-1);
	}
      iniBakFormat = 1; /* 1 = ascii */
      strcpy(paramFile, argv[2]); /* set iniBakFile string to the name */ 
      return 0;          /* 0 = continue a previuosly interrupted simulation */
    }

  else if (!strcmp(argv[1], "-fa"))
    {
      if (argc == 2)
	{
	  sprintf(TXT, "ERROR: -fa needs an ascii coordinate file name.\n");
	  mdPrintfR0(STD, TXT);
	  exit(-1);
	}
      iniCorFormat = 1; /* 1 = ascii */
      strcpy(paramFile, argv[2]); /* set paramFile string to the name 
				     of the parameter file */
      return 1;                   /* 1 = new simulation */
    }

  if ((!strcmp(argv[1], "-continue")) ||
      (!strcmp (argv[1], "-c")))
    {
      iniBakFormat = 0;
      return 0;          /* 0 = continue a previuosly interrupted simulation */
    }
  else if (!strcmp(argv[1], "-f"))
    {
      if (argc == 2)
	{
	  sprintf(TXT, "ERROR: -f needs a coordinate file name.\n");
	  mdPrintfR0(STD, TXT);
	  exit(-1);
	}
      strcpy(paramFile, argv[2]); /* set paramFile string to the name 
				     of the parameter file */
      return 1;                   /* 1 = new simulation */
    }
  else if (!strcmp(argv[1], "-del"))
    {
      sprintf(TXT, "Delete all temporary files(restore,status,check)? [n] ");
      mdPrintfR0(STD, TXT);
      car = MPIgetchar();
      if (car != 'y') exit(-1);
      
      delDoubleHD(BAK_FILE_NAME);
      if (OprogStatus.tapeTimes != 0) delDoubleTape(BAK_FILE_NAME);
      delDoubleHD(STATUS_FILE_NAME);
      delTmpHD(CF);
      delTmpHD(CHK_LOG);
      delTmpHD(MDS_LOG);
      /*  measure files not touched */
      exit(-1); 
    }
  else if (!strcmp(argv[1], "-dsf"))
    {
      sprintf(TXT, "Delete status files? [n] ");
      mdPrintfR0(STD, TXT);
      car = MPIgetchar();
      if (car != 'y') exit(-1);
      delDoubleHD(STATUS_FILE_NAME);
      exit(-1);
    }
  else if (!strcmp(argv[1], "-help"))   
    {
      sprintf(TXTA[0], "Syntax: mdsimul <option>\n");
      sprintf(TXTA[1], "Where <option> is one of the following:\n");
      sprintf(TXTA[2], "-c : continue a previuosly interrupted simulation\n");
      sprintf(TXTA[3], "-f <file> : begin a new simulation using <file>\n");
      sprintf(TXTA[4], "            to get simulation parameters (see SIMUL_INFO file)\n");
      sprintf(TXTA[5], "-del : delete all simulation temporary files\n");
      sprintf(TXTA[6], "-dst : delete only status files\n");
      sprintf(TXTA[7], "-list: list all valid parameters in file specified with -f\n");
      mdPrintfR0(STD, TXTA[0], TXTA[1], TXTA[2], TXTA[3], TXTA[4], TXTA[5],
		 TXTA[6], TXTA[7]);
      exit (-1);
    }
  else if (!strcmp(argv[1], "-list")) /* list all valid parameters in params 
					 file */
    {
      mdPrintfR0(STD, "Valid parameters in file specified with -f option are:\n");
      for (i=0; OsinglePar[i].ptr != NULL; ++i)
	{
	  sprintf(TXT, "%-15s ",OsinglePar[i].parName);
	  mdPrintfR0(STD, TXT);
	  if (  ( (i+1) % 4 ) == 0  ) 
	    {
	      /* every 4 params go to newline */
	      mdPrintfR0(STD, "\n");
	    }
	}  
      if (  ( (i+1) % 4 ) != 0  ) /* If there isn't a return print it */
	{
	  mdPrintfR0(STD, "\n");
	}
      exit(-1);
    }
  else 
    {
      sprintf(TXTA[0], "ERROR: Invalid argument!\n");
      sprintf(TXTA[1], "Try mdsimul -help for valid arguments\n");
      mdPrintfR0(STD, TXTA[0], TXTA[1]);
      exit (-1);
    }
  
  fflush(stdout);		/* flush output buffer (printf) */
}
extern int my_rank;
/* =========================== >> SuspendFather <<< =========================*/
void SuspendJob(int sigid)
{
#ifdef MPI
  if (my_rank == MPI_OUT_PROCESS)
#endif
    printf("[STEP %d] Simulation suspended, unlink status file\n", Oparams.curStep);
  fflush(stdout);
  /* delete both status file, so that chk_status can't restart simulation */
  delDoubleHD(STATUS_FILE_NAME);
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  exit(-1);
#endif
  //MPI_Abort(MPI_COMM_WORLD, -1);
  //  exit(-1);           /* on exit calls endf, see below */
}


/* ============================ >>> endFather <<< ===========================*/
void endJob(int ec, void* arg)
{
  /* This function is called at exit of father ( see mdsimul.c at the begin 
     of the FATHER code ) */ 
  //fclose(output); /*close log stream */
#ifdef MPI
 
  if (ec != -1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
    }
  else
    {
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
#else
  exit(-1);
#endif
}

/* ========================== >>> initBefore <<< ============================*/
void initBefore(void)
{
  /* DESCRIPTION:
     Initialize all program structure, that is :
     OsimStat, Oparams, OprogStatus.
      === >>> This procedure is called only by newsimul() <<< === */
  int jj;
  char tmpstr[NAME_LENGTH];
  /* =============== >>> Initialize Oparams structure <<< ===================
  */
  Oparams.parnum = 100;
  Oparams.steplength = 0.000001;	      /* 1 microsecondo */
  Oparams.totStep = 1000;
  Oparams.curStep = 1;

  /* ================ >>> Initialize OsimStat structure <<< ================ */
  OsimStat.running = 1;                       /* 1 = true = running 
						 (NOT USED!!!) */
  OsimStat.curStep = 1;

  /* ================= >>> Initialize progStatus struct <<< =================*/

  strcpy(OprogStatus.inifile, "*");	      /* random particle initialization
						 by default */
  strcpy(OprogStatus.endfile, "endfile.coor");/* defualt name for file 
					         containing final coords */

  strcpy(OprogStatus.xvafile, "tapefile.xva");/* default name for tape file */
  
  OprogStatus.tapeTimes = 0;   /* see mdsimdep.h */
   
  OprogStatus.bakSteps = 40;   /* steps between two restore file savings */ 
  
  OprogStatus.staSteps = 20;   /* steps between two status file savings */

  OprogStatus.xvaSteps = 0;    /* steps between two tape file savings, 0 means
				 no tape file */ 
  OprogStatus.fstps = 1;
  OprogStatus.xvaSaveMode = 0; /* default is linear mode */
  OprogStatus.NN = 0;         /* log block lenght is base^NN */
  OprogStatus.base = 1.3;      /* default base is 1.3 */
  
  strcpy(OprogStatus.tmpPath, "*");/* "*" = use macro for path */
  strcpy(OprogStatus.misPath, "*");/* "*" = use macro for path */

  OprogStatus.bakSaveMode = 0; /* default is linear mode */
  OprogStatus.bakStepsAscii = 0;

  strcpy(PT_temps, ""); 
   /* Files containing measure are named by default: 
     Misura1- Misura2- etc. */
  for (jj = 0; jj < NUM_MISURE; ++jj)
    {
      OprogStatus.measSteps[jj] = 0; /* By default don't
					save any measure (these should be
					set in params file) */
      OprogStatus.measCalc[jj] = -1; 
      /* NOTE:
	 -1 => OprogStatus.measCalc[i] = OprogStatus.measSteps[i] if 
	 the user doesn't specify the calculation interval in the parameters
	 file */
      
      OprogStatus.initStep[jj] = 1;
      OprogStatus.initCalc[jj] = 1;
      strcpy(OprogStatus.dataFiles[jj], "measure");
      sprintf(tmpstr, "%d", jj);
      strcat(OprogStatus.dataFiles[jj], tmpstr);
      strcat(OprogStatus.dataFiles[jj], "-");
    }
  usrInitBef();                                /* user initializations */
}

/* ================== >>> parallel_tempering_init <<< =================== */
void parallel_tempering_init(void)
{
  /* La sintassi del parametro e' la seguente:
     <T1>-<T2>
     e allora assegna i lambda in modo da avere numOfProcs temperature
     equispaziate a partire da T1 fino a T2 
    
     OPPURE
     
     T1,T2, ...T<numOfProcs>
     cioe' si danno tutte le temperature esplicitamente
 
  */
  char* strN;
  double Tini, Tfin, DT, Ti, Lini, Lfin, DL, Li;
  char ptt[MAXPROCS*4]; 
  int i, PTM, ss;
  if (strlen(PT_temps) == 0) 
    return;
  //printf("PT string parsing PT_temps: %s\n", PT_temps);
  strcpy(ptt, PT_temps);
  strN = strtok(ptt, ",");
  PTM = Oparams.PTM;

  //printf("QUA!!! PT init\n");
  if (PTM == 1)
    {
      printf("[ERRORE] Parallel tempering with only one system?!??!?!\n");
      exit (-1);
    }
  if (!strcmp(strN,PT_temps))
    {
      if (sscanf(PT_temps, "L[%lf - %lf]", &Tini, &Tfin) == 2)
	{
	  /* Costruisce una sequenza lineare nei lambda a partire 
	     dalle temperature */
	  Lfin = Oparams.T / Tfin;
	  Lini = Oparams.T / Tini;
	  DL = (Lfin - Lini) / ((double)PTM - 1);
	  Li = Lini;
	  for (i = 0; i < PTM; i++)
	    {
	      /* le temperature a cui equilibreranno
		 i vari sistemi sono:
		 T_i = T_0 / lambda_i
		 (ved. articolo di Kob 
		 "Equilibrating Glassy Systems with Parallel Tempering"
	      */
	      Oparams.lambda0[i] = Li;
	      printf("T[%d]: %f lambda: %f\n", i, 
		     Oparams.T / Li, Oparams.lambda0[i]);
	      Li += DL;
	    }
	}
      else if ((sscanf(PT_temps,"T[%lf - %lf]", &Tini, &Tfin) == 2) ||
	       (sscanf(PT_temps,"%lf - %lf", &Tini, &Tfin) == 2))
	{
	  /* Costruisce una sequenza lineare di temperature */
	  DT = (Tfin - Tini) / ((double)PTM - 1);
	  Ti = Tini;
	  for (i = 0; i < PTM; i++)
	    {
	      /* le temperature a cui equilibreranno
		 i vari sistemi sono:
		 T_i = T_0 / lambda_i
		 (ved. articolo di Kob 
		 "Equilibrating Glassy Systems with Parallel Tempering"
	      */
	      Oparams.lambda0[i] = Oparams.T / Ti;
	      printf("T[%d]:%f lambda: %f\n", i, Ti, Oparams.lambda0[i]);
	      Ti += DT;
	    }
	}
      else
	{
	  printf("[ERROR] Syntax error in parallel tempering string!!!\n");
	  exit(-1);
	}
      //printf("ECCOMI QUA Tfin: %f Tini: %f\n", Tini, Tfin);
     
    }
  else
    {
      i = 0;
      sscanf(strN, "%lf", &Ti);
     
      Oparams.lambda0[0] = Oparams.T / Ti;
    
      for(i=1; (strN = strtok(NULL, ","))!=NULL; i++)
	{
	  sscanf(strN, "%lf", &Ti);
	  Oparams.lambda0[i] = Oparams.T / Ti;
	}
      if ((strN == NULL) && (i != Oparams.PTM))
	{
	  printf("[ERROR] numeber of lambdas must be equal to numOfProcs!!!\n");
	  exit(-1);
	}
    }
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      Oparams.lambdat[ss] = ss;
    }
}
