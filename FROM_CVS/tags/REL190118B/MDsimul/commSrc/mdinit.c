#undef MAIN
#include<mdsimul.h>

/* WARNING: Extern variables must be of the same type of the real definition */

/* ====================== >>> Descriptors and streams <<< ================== */

extern FILE* output;

/* ======================= >>> Simulation vars <<< ==========================*/

extern unsigned char BAK, STA;       /* global switches for HD(0/1) */
extern unsigned char BAKT;           /* global switch for tape (0/1) */

extern int logBlock;

extern int pid;
extern int SEGSIZE;
/* ======================== >>> Structures <<< ============================= */
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);
extern struct simStat OsimStat;
extern char msgStrA[MSG_LEN], msgStrB[MSG_LEN], msgStrC[MSG_LEN];
extern char paramFile[NAME_LENGTH];
extern char TXT[MSG_LEN], TXTA[10][MSG_LEN];
#ifdef MPI
extern int my_rank;
#else
int my_rank;
#endif
extern int printerRank, numOfProcs;
/* ADDED 23/04/2001 */
int iniCorFormat = 0, iniBakFormat = 0; /* 0 = binary 1 = ascii */
char iniBakFile[NAME_LENGTH];

void reread(char* );
extern void readBakAscii(char* );
extern void readCorAscii(char* );

char bakList[512][NAME_LENGTH];
char corList[512][NAME_LENGTH];
char iniBakListFile[NAME_LENGTH];

/* ADDED [reread] 07/05/01 */
int rereadbool = 0;
void init_rng(int mdseed, int mpi, int my_rank);

/* ============================= >>> Continue <<< ===========================*/
void Continue(void)
{
  /* Load restore file using an algorithm explained in the TECH_INFO
     file */
  FILE* il;
  char firest[NAME_LENGTH];
  int i, rnk; 
  char perc[NAME_LENGTH];

   /* Opens the restore file that contains all informations
     to restart the simulation (see. TECH_INFO file for details).
     At least one file is not corrupted so if it fails to open 
     one file, it tries the second, for explanations on how it decides
     the file to use for restoring see TECH_INFO.
  */

  newSim = 0; /* New = 0 means: not a new simulation, i.e. continuing a
		 previously interrupted one */
#if !defined(MPI)
  printf("setting seed:%d time:%d\n", OprogStatus.mdseed, (int)time(NULL));
  init_rng(OprogStatus.mdseed, 0, -1);
#endif  
  /* 30/10/2005: aggiunta l'inizializzazione anche quando si continua!
   *             Verificare che funzioni! */  
  initBefore();		/* initialize Oparams with default settings */

  /* <------------------------------------------------------- OPEN LOG FILE */
  openLog("a"); /* append: continue to write on previously created log file */

  /* ======================= >>> CHOOSE RESTORE FILE <<< ===================*/
 
  /* NOTE: the restore file on disk is the only  used to restart simulation ,
     in fact the restore files on Tape should be used manually copying them
     on temporary dir on disk.
     So if a failure occur and both restore files are currupted on this 
     you must restart the simulation manually */

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
      if (BAK==0)
	strcpy(inifile_for_mgl,"COORD_TMP1.mgl");
      else
	strcpy(inifile_for_mgl,"COORD_TMP0.mgl");
    }
  else if (iniBakFormat == 2)
    {
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    "Trying to restart simulation",
	    "Reading ascii restore file list on Hard Disk...",
	    NULL);
      if ((il = fopen(iniBakListFile, "r")) == NULL)
	{
	  mdPrintf(ALL, "Error opening bak list file\n", NULL);
	  exit(-1);
	}
      
      i = 0;
      /* La prima linea deve contenere il percorso */
      fscanf(il, "%s\n", perc);
      while (!feof(il))
	{
	  /* mette il nome del file letto nella posizione relativa al rango */
	  fscanf(il, "%[^_]_R%d\n", firest, &rnk);
	  printf("1:%s 2: %d\n", firest, rnk);
	  sprintf(bakList[rnk], "%s%s_R%d", perc, firest, rnk);
	  i++;
	}
      fclose(il);
      readBakAscii(bakList[my_rank]);
      strcpy(inifile_for_mgl, bakList[my_rank]);

#if 0     
      for (i = 0; i < Oparams.PTM; i++)
	{ 
	  printf("[%d] %s\n", my_rank, bakList[my_rank]);
	}
#endif
    }
  else if (iniBakFormat == 1)
    {
      mdMsg(ALL, NOSYS, "Restore", "NOTICE", NULL,
	    "Trying to restart simulation",
	    "Reading ascii restore file on Hard Disk...",
	    NULL);
#ifdef MPI 
     strcpy(firest, iniBakFile);
     sprintf(iniBakFile, "%s_R%d", firest, my_rank);
#endif
     readBakAscii(iniBakFile);
     strcpy(inifile_for_mgl, iniBakFile);
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
  /* AGGIUNTA CONDIZIONE 24/04/01 */
  if (OprogStatus.staSteps != 0)
    chooseStatus( absTmpHD(STATUS_FILE_NAME), NULL);
  /* second arg not used actually !!!! */
  
  /* second arg not used actually !!!! */
  if (rereadbool)
    {
      //printf("[%d] QUI!!!!!!%s\n", my_rank, paramFile);
      reread(paramFile);
    }

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
    }     
  /* ============== >>> Initialize OxvaHead structure <<< ===================*/
  //printf("SEGSIZE: %d\n", SEGSIZE);

#ifdef BIMIX
  SEGSIZE=(Oparams.parnum[0]+Oparams.parnum[1])*sizeof(double);
#else
  SEGSIZE=Oparams.parnum*sizeof(double);
#endif
  OxvaHead.saveSteps = OprogStatus.xvaSteps;
  OxvaHead.size = SEGSIZE * XVA_NUM; 
  /* ADD 13/09/2000: */
  OxvaHead.mode = OprogStatus.xvaSaveMode;
  OxvaHead.NN = OprogStatus.NN;
  OxvaHead.base = OprogStatus.base;
#if !defined(MD_HARDSPHERES) && !defined(MC_SIMUL)
  OxvaHead.dt = Oparams.steplength;
#endif
#ifdef BIMIX
    OxvaHead.parnum[0] = Oparams.parnum[0];
    OxvaHead.parnum[1] = Oparams.parnum[1];
#else
    OxvaHead.parnum = Oparams.parnum;
#endif
  OxvaHead.T = Oparams.T;
#if !defined(MD_HARDSPHERES) && !defined(MC_SIMUL)
  OxvaHead.Vol = Vol;
#endif
  logBlock = (int) rint(pow(OprogStatus.base, (double) OprogStatus.NN+1));

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
#ifndef MDLLINT
	    case LLINT:
#endif
	    case INT :
	      *(int*) OsinglePar[i].ptr = atoi(stringB);
	      break;
#ifdef MDLLINT
	    case LLINT :
	      *(long long int*) OsinglePar[i].ptr = atoll(stringB);
	      break;
#endif
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


/* =========================== >>> reread <<< ============================*/
void reread(char* argom)
{
  /* DESCRIPTION:
     legge tutti i parametri scartando ovviamente del numero di
     particelle e il file di coordinate iniziali
  */
  
  char str1[255],str2[255];            /* used to scan parameters file */
  FILE* pfs;
  int parNumYN = 0;
#if 0  
  initYN = 0;
#endif 
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
#ifdef BIMIX      
      if ( !strcmp(str1, "parnumA") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum[0] = atoi(str2);
	}
      if ( !strcmp(str1, "parnumB") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum[1] = atoi(str2);
	}
#ifdef LMC
      if ( !strcmp(str1, "lattice_M") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1;
	  Oparams.lattice_M = atoi(str2);
	}
      if ( !strcmp(str1, "Vol") && 
	   strcmp(str2, "*") )
	{
	  Vol = atof(str2);
	}

#endif
#else
      if ( !strcmp(str1, "parnum") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum = atoi(str2);
	} 
#endif
	} 
    }  
  
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
#if defined(BIMIX) && defined(LMC)
  Oparams.lattice_a = cbrt(Vol) / ((double)Oparams.lattice_M);
#endif
}

/* =========================== >>> scanFile <<< ============================*/
void getseed(char* argom)
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

      if (!strcmp(str1, "seed")) 
	{
	  OprogStatus.mdseed = atoi(str2);
	  break;
	}
    }  
  
  rewind(pfs);
  fclose(pfs);
}
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

#ifdef BIMIX      
      if ( !strcmp(str1, "parnumA") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum[0] = atoi(str2);
	}
      if ( !strcmp(str1, "parnumB") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum[1] = atoi(str2);
	}
#ifdef LMC
      if ( !strcmp(str1, "lattice_M") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1;
	  Oparams.lattice_M = atoi(str2);
	}
      if ( !strcmp(str1, "Vol") && 
	   strcmp(str2, "*") )
	{
	  Vol = atof(str2);
	}

#endif
#else
      if ( !strcmp(str1, "parnum") && 
	   strcmp(str2, "*") )
	{
	  parNumYN = 1; /* parnum parameter found */
	  Oparams.parnum = atoi(str2);
	} 
#endif
      
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
	  strcpy(inifile_for_mgl, str2);
	}
    }

  /* Se il campo inifile di OprogStatus contiene un asterisco (*) 
     allora il set di coordinate iniziale deve essere generato 
     in qualche modo ( ved. dopo o SIMUL_INFO ) */
  if (!strcmp(OprogStatus.inifile, "*"))
    {
#ifdef EDHE_FLEX
      printf("You have to supply an inifile if EDHE_FLEX macro is defined!\n");
      exit(-1);
#endif
      mdMsg(ALL, NOSYS, "Init", "NOTICE", NULL,
	    "Generating initial coordinates...",
	    NULL);
      /* allocate  coordinates needed by simulation (see. mdsimul_p) */
      
#ifdef BIMIX
#ifdef LMC
      SEGSIZE = sizeof(int) * Oparams.parnum[0];
      AllocCoord(SEGSIZE, ALLOC_LISTA,
		 NULL);
      SEGSIZE = sizeof(int) * Oparams.parnum[1];
      AllocCoord(SEGSIZE, ALLOC_LISTB,
		 NULL);
#else
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[0];
      AllocCoord(SEGSIZE, ALLOC_LISTA,
		 NULL);
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[1];
      AllocCoord(SEGSIZE, ALLOC_LISTB,
		 NULL);
#endif
#else
#ifdef MD_ALLOC_POLY
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
      /* SEGSIZE is the size in bytes of an array of coordinates */
      
      /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
	 all addresses of the coordinates declared in the simulaiton
	 (see that file) */
      AllocCoordPoly(SEGSIZE, ALLOC_LIST,
		 NULL);
#else
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
      /* SEGSIZE is the size in bytes of an array of coordinates */
      
      /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
	 all addresses of the coordinates declared in the simulaiton
	 (see that file) */
      AllocCoord(SEGSIZE, ALLOC_LIST,
		 NULL);
#endif
#endif      
      initYN = 1; /* Yes, initialize coordinates after reading 
		     all parameters */
    }
  else
    {
      if (parNumYN) /* Is parnum parameter specified? */
	{
	  /* if you specify an initial file particles number is not needed, 
	     because the particle number is set in coordinate file */
	  mdMsg(ALL, NOSYS, "Init", "WARNING", NULL,
		"You specify an initial coordinates file, so you have not to give",
		"particles number, I ignore parnum parameter",
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
      
      if (sscanf(line, "%[^:# ] : %[^\n# ] ", str1, str2) < 2)
	continue;
    
      //fscanf(pfs, " %[^: ] : %[^'\n' ]", str1, str2);
      /* analyzes the parameter just read 
         This function depends strongly upon simulation */
      Parsing(str1, str2);
    }
  fclose(pfs);	  /* Aggiungere controllo della chiusura <-----!!! */
  
#if defined(BIMIX) && defined(LMC)
  Oparams.lattice_a = cbrt(Vol) / ((double)Oparams.lattice_M);
#endif  
  if (initYN) /* If not initial coords file was supplied initialize them */ 
    {
      initCoord();
    }
}
char rng_array[256];
void init_rng(int mdseed, int mpi, int my_rank)
{
  int fact, k;
  long rngseed;

  if (mpi==0)
    fact = 1;
  else
    fact = my_rank + 1;
 if (mdseed>=0) 
    {
#if defined(MD_RAND48)      
      printf("[RAND48] Initializing RNG with prefixed seed...\n");
      srand48(fact*mdseed);
#elif defined(MD_RANDOM)
#if 1
      printf("[RANDOM] Initializing RNG with prefixed seed...\n");
      initstate(fact*mdseed, rng_array, 128);
      setstate(rng_array);
#else
      srandom(fact*mdseed);
#endif
#else
      printf("[RAND] Initializing RNG with prefixed seed...\n");
      srand(fact*mdseed);
#endif
    }
  else
    {
#if defined(MD_RAND48)
      rngseed=fact*((int)time(NULL));
      printf("[RAND48] Initializing RNG with random seed (=%ld)...\n", rngseed);
      srand48(rngseed);
#elif defined(MD_RANDOM)
#if 1
      rngseed=fact*((int)time(NULL));
      printf("[RANDOM] Initializing RNG with random seed (=%ld)...\n",rngseed);
	/* se 8 < size < 32 allora random usa un algoritmo simile a drand48 mentre 
	  se size > 32 ne usa uno differente, il default è size=128.
	  Dai test effettuati al 19/05/11 sembrerebbe che o la distribuzione non è uniforme o comunque 
	  che numeri successivi siano un po' correlati se si usano 256 byte, da 128 bytes in giù le cose 
	  migliorano molto. */ 
      initstate(rngseed, rng_array, 128);
      setstate(rng_array);
#else
      srandom(fact*((int)time(NULL)));
#endif
#else
      rngseed=fact*((int)time(NULL));
      printf("[RAND] Initializing RNG with random seed (=%ld)...\n", rngseed);
      srand(rngseed);
#endif
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
  
  output = stdout;
  newSim = 1;
#if !defined(MPI)
  getseed(argom);
  printf("setting seed:%d time:%d\n", OprogStatus.mdseed, (int)time(NULL));
  init_rng(OprogStatus.mdseed, 0, -1); 
#if 0
  if (OprogStatus.mdseed>=0) 
    srand(OprogStatus.mdseed);
  else
    srand(((int)time(NULL)));
#endif
#endif

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
#if !defined(MD_HARDSPHERES) && !defined(MC_SIMUL)
  OxvaHead.dt = Oparams.steplength;
#endif
#ifdef BIMIX
  OxvaHead.parnum[0] = Oparams.parnum[0];
  OxvaHead.parnum[1] = Oparams.parnum[1];
#else
  OxvaHead.parnum = Oparams.parnum;
#endif
  OxvaHead.T = Oparams.T;
#if !defined(MD_HARDSPHERES) && !defined(MC_SIMUL)
  OxvaHead.Vol = Vol;
#endif
  OprogStatus.savedXva = -1;
  OprogStatus.fstps = 1;
  logBlock = (int) rint(pow(OprogStatus.base, (double) OprogStatus.NN+1));
  printf("logBlock: %d base:%f NN:%d\n", logBlock, OprogStatus.base, OprogStatus.NN);
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
void print_usage(void)
{
  printf("Syntax: mdsimul <option>\n");
  printf("Where <option> is one of the following:\n");
  printf("-c : continue a previuosly interrupted simulation\n");
  printf("-ca <ascii_store_file : continue using the ascii file <ascii_store_file>\n");
  printf("-f <file> : begin a new simulation using <file>\n");
  printf("            to get simulation parameters\n");
  printf("-fa <file>:  same as -f but starting from an ascii file\n");
  printf("-del : delete all simulation temporary files\n");
  printf("-dsf : delete only status files\n");
  printf("-list: list all valid parameters in file specified with -f\n");
  printf("-creread <param_file>: continue re-reading the parameter file <param_file>\n");
  printf("-mgl <mgl_mode>: <mgl_mode> = 0 -> save normal ascii file\n");
  printf("                              1 -> save ascii files in mgl format\n");
  printf("                              2 -> save ascii file in mgl format and exit\n");
  exit (-1);
}
#ifdef MD_PROTEIN_DESIGN
char nativeConf[512];
int nativeConfYN=0;
#endif
/* ========================= >>> args <<< =================================*/
int argsMd(int argc, char** argv)  
{
  int car, ret=-1, i, cc=1;
  /* You must specify -c or -f <filename> as argument */ 
  if ( argc == 1 ) 
    {
      sprintf(TXTA[0], "What should I do? (No Arguments)\n");
      sprintf(TXTA[1], "Try 'mdsimul -help' for valid arguments\n");
      mdPrintfR0(STD, TXTA[0], TXTA[1], NULL);
      exit(-1);
    }
  while (cc < argc)
    {
      if ((!strcmp(argv[cc], "-continue")) ||
	  (!strcmp (argv[cc], "-c")))
	{
	  iniBakFormat = 0; /* 0 = binary */
	  //strcpy(paramFile, argv[2]); /* set paramFile string to the name */

	  ret=0;          /* 0 = continue a previuosly interrupted simulation */
	}
      /* ADDED 04/05/2001 */
      else if ((!strcmp(argv[cc], "-continue_ascii_list")) ||
	       (!strcmp (argv[cc], "-cal")))
	{
	  cc++;
	  if (cc == argc)
	    {
	      sprintf(TXT, "ERROR: -cal needs a file name.\n");
	      mdPrintfR0(STD, TXT, NULL);
	      exit(-1);
	    }
	  iniBakFormat = 2; /* 2 = lista di ascii file */
	  strcpy(iniBakListFile, argv[cc]); /* set iniBakFile string to the name */ 
	  ret=0;         /* 0 = continue a previuosly interrupted simulation */
	}
      else if ((!strcmp(argv[cc], "-mgl_mode")) ||
	       (!strcmp (argv[cc], "-mgl")))
	{
	  cc++;
	  if (cc==argc)
	    {
	      printf("You should supply the mgl mode (0,1,2)\n");
	      printf("0 = Store files saved not in mgl format [Default].\n");
	      printf("1 = Store files saved in mgl format.\n");
	      printf("2 = save configuration in mgl format and exit\n");
	      exit(-1);
	    }
	  mgl_mode = atoi(argv[cc]);
	} 
      else if ((!strcmp(argv[cc], "-continue_ascii")) ||
	       (!strcmp (argv[cc], "-ca")))
	{
	  cc++;
	  if (argc == cc)
	    {
	      sprintf(TXT, "ERROR: -ca needs a ascii restart file name.\n");
	      mdPrintfR0(STD, TXT, NULL);
	      exit(-1);
	    }
	  iniBakFormat = 1; /* 1 = ascii */
	  strcpy(iniBakFile, argv[cc]); /* set iniBakFile string to the name */ 
	  strcpy(inifile_for_mgl, argv[cc]);
	  ret=0;          /* 0 = continue a previuosly interrupted simulation */
	}
#ifdef MD_PROTEIN_DESIGN
      else if (!strcmp(argv[cc], "-nc"))
	{
	  cc++;
	  if (argc == cc)
	    {
	      sprintf(TXT, "ERROR: -nc needs an ascii coordinate file name.\n");
	      mdPrintfR0(STD, TXT, NULL);
	      exit(-1);
	    }
	  nativeConfYN=1;
	  strcpy(nativeConf, argv[cc]); /* set paramFile string to the name 
					  of the parameter file */
	}
#endif
      else if (!strcmp(argv[cc], "-fa"))
	{
	  cc++;
	  if (argc == cc)
	    {
	      sprintf(TXT, "ERROR: -fa needs an ascii coordinate file name.\n");
	      mdPrintfR0(STD, TXT, NULL);
	      exit(-1);
	    }
	  iniCorFormat = 1; /* 1 = ascii */
	  strcpy(paramFile, argv[cc]); /* set paramFile string to the name 
					  of the parameter file */
	  ret=1;                   /* 1 = new simulation */
	}
      else if (!strcmp(argv[cc], "-f"))
	{
	  cc++;
	  if (argc == cc)
	    {
	      sprintf(TXT, "ERROR: -f needs a coordinate file name.\n");
	      mdPrintfR0(STD, TXT, NULL);
	      exit(-1);
	    }
	  iniCorFormat = 0; /* 0 = binary */
	  strcpy(paramFile, argv[cc]); /* set paramFile string to the name 
					  of the parameter file */
	  ret = 1;                   /* 1 = new simulation */
	}
      else if (!strcmp(argv[cc], "-del"))
	{
	  sprintf(TXT, "Delete all temporary files(restore,status,check)? [n] ");
	  mdPrintfR0(STD, TXT, NULL);
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
      else if (!strcmp(argv[cc], "-dsf"))
	{
	  sprintf(TXT, "Delete status files? [n] ");
	  mdPrintfR0(STD, TXT, NULL);
	  car = MPIgetchar();
	  if (car != 'y') exit(-1);
	  delDoubleHD(STATUS_FILE_NAME);
	  exit(-1);
	}
      /* ADDED 22/03/2001: Rilegge i parametri */	
      else if (!strcmp(argv[cc], "-creread"))
	{
	  cc++;
	  if (argc == cc)
	    {
	      sprintf(TXT, "ERROR: -creread needs a paramter file name.\n");
	      mdPrintfR0(STD, TXT, NULL);
	      exit(-1);
	    }
	  strcpy(paramFile, argv[cc]);
	  mdPrintfR0(STD, "Re-reading parameter file...", NULL);
	  rereadbool = 1;
	  ret=0; /* -reread = -c rileggi il file di parametri */
	}
      else if (!strcmp(argv[cc], "-help"))   
	{
	  print_usage();	  
	}
      else if (!strcmp(argv[cc], "-list")) /* list all valid parameters in params 
					      file */
	{
	  mdPrintfR0(STD, "Valid parameters in file specified with -f option are:\n", NULL);
	  for (i=0; OsinglePar[i].ptr != NULL; ++i)
	    {
	      sprintf(TXT, "%-15s ",OsinglePar[i].parName);
	      mdPrintfR0(STD, TXT, NULL);
	      if (  ( (i+1) % 4 ) == 0  ) 
		{
		  /* every 4 params go to newline */
		  mdPrintfR0(STD, "\n", NULL);
		}
	    }  
	  if (  ( (i+1) % 4 ) != 0  ) /* If there isn't a return print it */
	    {
	      mdPrintfR0(STD, "\n", NULL);
	    }
	  exit(-1);
	}
      else if (cc == argc)
	print_usage();
      else 
	{
	  sprintf(TXTA[0], "ERROR: Invalid argument!\n");
	  sprintf(TXTA[1], "Try mdsimul -help for valid arguments\n");
	  mdPrintfR0(STD, TXTA[0], TXTA[1], NULL);
	  exit (-1);
	}
      cc++;
    }
  if (ret==-1)
    {
      printf("ERRORE: invalid arguments!\n");
      printf("Try mdsimul -help for valid arguments\n");
      exit(-1);
    }
  fflush(stdout);		/* flush output buffer (printf) */
  return ret;
}
extern int my_rank;
/* =========================== >> SuspendFather <<< =========================*/
void SuspendJob(int sigid)
{
#ifdef MPI
  if (my_rank == MPI_OUT_PROCESS)
#endif
    printf("[STEP %lld] Simulation suspended, unlink status file\n", 
	   (long long int)Oparams.curStep);
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

#if defined(CRAY) || defined(ORIGIN)
void endJob(void)
{}
#else
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
#endif
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
#ifdef BIMIX
  Oparams.parnum[0] = 100;
  Oparams.parnum[1] = 100;
#else
  Oparams.parnum = 100;
#endif
#if !defined(MD_HARDSPHERES) && !defined(MC_SIMUL)
  Oparams.steplength = 0.000001;	      /* 1 microsecondo */
#endif
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
  
  strcpy(OprogStatus.tmpPath, "*");/* "*" = use macro for path */
  strcpy(OprogStatus.misPath, "*");/* "*" = use macro for path */

  OprogStatus.tapeTimes = 0;   /* see mdsimdep.h */
   
  OprogStatus.bakSteps = 10000;   /* steps between two restore file savings */ 
  
  
  OprogStatus.staSteps = 0;   /* steps between two status file savings */

  OprogStatus.xvaSteps = 0;    /* steps between two tape file savings, 0 means
				 no tape file */ 
  OprogStatus.fstps = 1;
#ifdef MD_BILOG
  OprogStatus.basew = 1.3;
  OprogStatus.lastbilogsaved = 0;
#endif
  OprogStatus.xvaSaveMode = 0; /* default is linear mode */
  OprogStatus.bakSaveMode = 0; /* default is linear mode */
  OprogStatus.bakStepsAscii = 0;

  OprogStatus.NN = 0;         /* log block lenght is base^NN */
  OprogStatus.base = 1.3;      /* default base is 1.3 */

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

