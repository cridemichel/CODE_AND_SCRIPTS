/*      $Id: mdsimul.c,v 1.20 2011-04-05 07:48:25 demichel Exp $     */
#ifndef lint
//static char vcid[] = "$Id: mdsimul.c,v 1.20 2011-04-05 07:48:25 demichel Exp $";
#endif /* lint */
/* Sintassi: mdsimul -f <nomefile> 
   dove <nomefile> e' il nome del file contenente i parametri della 
   simulazione.
   Tale file deve essere della forma:
   <PARAM>:<VALORE> , dove PARAM puo'assumere uno dei valori contenuti 
   nell'array pars_arr (ved. dopo).
   Se tale opzione viene (cioe'  si lancia semplicemente mdsimul ) allora 
   si assume che si tratta di una simulazione interrotta e viene caricata
   direttamente la struttura simul_pars contenuta nel file PARS_FILE */
#define MAIN
#include<mdsimul.h>
#ifdef MPI
int *equilibrated;
#endif
struct simStat OsimStat;	     /* global istance of sim_status 
					structure */
/* ----------------------------------------------------------------------*/

int SEGSIZE;         /* arrays dimension in bytes 
				 ( set in NewSimul() or Continue(),
				 see mdinit.c file )*/
unsigned char BAK, STA;       /* each is a switch (0/1) which determine 
				 the next file to save to on HD 
				 ("double unsigned") */
unsigned char BAKT;           /* switch used for tape (analogous to BAK for 
		 		 Tape)*/ 
				 
extern char TXT[MSG_LEN];
/* shared array of two integers for active lock/unlock (semaphores) */

/* log stream */
extern FILE* output;

/* strings to store messages (use calling mdMsg()) */
extern char msgStrA[MSG_LEN], msgStrB[MSG_LEN], msgStrC[MSG_LEN]; 

/* variables that counts the number of accesses to function 
   doubleBufBak and doubleSaveMeasure respectively */
int bakTimes = 0, measureTimes = 0, xvaTimes = 0;
char paramFile[NAME_LENGTH];

int ENDSIM = 0;

extern int printerRank;

#if defined(MPI)
int MPIpid;
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
#endif 


/* ============================= >>> commMD <<< ======================= */
#if !defined(CRAY) && !defined(ORIGIN)
inline
#endif
extern void saveBakAscii(char *fn);
extern void saveCorAscii(void);
extern int chkBakAsciiSteps(void);
extern void doubleSaveMeasure(int PN, int misNum, int* times, int tapeTimes);
void commMD(void)
{
  int i;   
  time_t tempo;         /* returned by time(NULL) */
#ifdef MDLLINT
  long long int msteps;
#else
  int msteps;
#endif
#ifdef MPI
  mdMsg(ALL, NOSYS, "Begin", "NOTICE", NULL,
	    "Farm MPI simulation",
	    NULL);
#else
  mdMsg(ALL, NOSYS, "Begin", "NOTICE", NULL,
	    "Single Job simulation",
	NULL);
#endif
  
  tempo = time(NULL);
  sprintf(msgStrA, "INITIAL TIME: %s", ctime(&tempo));
  mdMsg(ALL,NOSYS, "Begin", "NOTICE", NULL,
	msgStrA,
	NULL);
      
  //printf("ECCOCI QUA RANK %d\n", my_rank);
  //fflush(stdout);
  /* <---------------------------------------------- NEW STEP (FATHER) */
  while ( (Oparams.curStep <= Oparams.totStep) && (ENDSIM == 0) )
    {
      //printf("rank[%d] vx[10]: %f\n", my_rank, vx[10]);
      //ProcSync0();

      move(); /* syncronization is inside move() (see mdfuncSD.c )*/
      
      /* <------------------------------------------- MEASURES (FATHER)
	 ACTUALLY NOT PARALLELIZED !!! */
      
      calc(); 

      /* --------------------------------------------------------------- */
      
      /* <--------------------------------------------------- SAVE MEASURE 
	 Every 'OprogStatus.measSteps[i]' steps save measure structure 'i' 
	 into the file Ofilenames.data_files[i] (see TECH_INFO file) */
      
      for (i = 0; Omeasure[i].buf != NULL; ++i) 
	/* NULL in the buf pointer of measure sruct means end of list */
	{
	  /* ADDED 09/04/2001: Se negativi glil steps allora salva sopra */
#ifdef MDLLINT
	  msteps = llabs(OprogStatus.measSteps[i]);
#else
	  msteps = abs(OprogStatus.measSteps[i]);
#endif
	  if ( (msteps != 0) &&
	       ((Oparams.curStep % msteps) == 0) &&
	       ((Oparams.curStep / msteps) >=
		OprogStatus.initStep[i]))

	    /* It begins to save from step: 
	       Oprogstatus.initStep[i] * OprogStatus.measStep[i], 
	       that is the initial step is given in unit of 
	       OprogStatus.measStep[i] */
	    
	    {
	      //printf("curStep: %d initStep[%d]: %d", Oparams.curStep, 
	      //     i, OprogStatus.initStep[i]);

	      /* makes two savings of the i-th measure, the second 
		 arguments indicate which measure to save first (must 
		 be a pointer), the third arguments must be a pointer 
		 to a global variable that store the number of access 
		 to the function, the fourth  arguments says to save 
		 on zip every 'OprogStatus.measSteps[i] * <value>' steps, 
		 where value is the third arg value. (if value is 0, 
		 don't save on zip) */
	      doubleSaveMeasure(OprogStatus.measSteps[i], 
				i, &measureTimes, OprogStatus.tapeTimes);
#if 0
	      printf("saving measure n.%d step: %d [%d] name: %s\n", i, Oparams.curStep,
		     OprogStatus.measSteps[i], OprogStatus.dataFiles[i]);
#endif
	      /* save measures on Tape and on HD doublely
		 i -> measure number 
		 measureTimes -> number of access to doubleSavemeasure
		 OprogSatus.tapeTimes -> multiplicative factor for 
		 Tape savings */ 
	      
	      //sprintf(TXT, "measure %s saved\n", 
	      //OprogStatus.dataFiles[i]);
	      //mdPrintf(STD, TXT, NULL);
	    }
	}
      /* this message  will be removed <------------!!!!!!!!!! */
      //mdMsg(STD, NOSYS, NULL, "NOTICE", NULL,
      //	"All measures saved",
      //	NULL);
      
      /* <-------------------------------------------- SAVE XVA FILE 
	 NOTE: xva file is the file containing positions, velocities and
	 accelerations save every 'OprogStatus.xvaSteps' steps*/
      if (chkXvaSteps())
	//((Oparams.curStep % OprogStatus.xvaSteps) == 0)  )
	{
	  UPDATE_SYSTEM;
	  saveXva(&xvaTimes, OprogStatus.tapeTimes); 
	  /* don't use double save in this case, anyway save also on tape*/
	  mdMsg(STD, NOSYS, NULL, "NOTICE", NULL, /* only on screen */
		"Saved on tape file",
		NULL);
	}

      /* <--------------------------------------------- SAVE RESTORE FILES
	 save restore files every OprogStatus.bakSteps steps 
	 (see header file) */
      if ( (OprogStatus.bakSteps != 0) &&
	   ((Oparams.curStep % OprogStatus.bakSteps) == 0) )
	{
	
	  UPDATE_SYSTEM;
	  doubleBufBak(&BAK, &BAKT, &bakTimes, OprogStatus.tapeTimes);	
	  /* save restore datas on Tape and on HD, BAK is the switch
	     for HD and BAKT is the switch for Tape 
	     bakTimes -> number of access to doubleBufBak
	     OprogStatus.tapeTimes -> multiplicative factor for 
		                          savings on Tape 
	  */
	  mdMsg(STD, NOSYS,NULL, "NOTICE", NULL,
		"Restore file saved",
		NULL);
	}

      /* ADDED 09/04/2001: restart file in ascii format */
      if (chkBakAsciiSteps())
	{
	  UPDATE_SYSTEM;
	  saveBakAscii(NULL);	
	  /* save restore datas on Tape and on HD, BAK is the switch
	     for HD and BAKT is the switch for Tape 
	     bakTimes -> number of access to doubleBufBak
	     OprogStatus.tapeTimes -> multiplicative factor for 
		                          savings on Tape 
	  */
	  mdMsg(STD, NOSYS,NULL, "NOTICE", NULL,
		"ASCII Restore file saved",
		NULL);
	}

      /* <--------------------------------------------------- SAVE STATUS 
	 Every OprogStatus.staSteps save simStat structure into the file 
	 STATUS_FILE (see TECH_INFO file) */
      if ( (OprogStatus.staSteps != 0) &&
	   ((Oparams.curStep % OprogStatus.staSteps) == 0) )
	{
	  doubleBufStatus();  
	  /* save staus only on HD */
	  mdMsg(STD, NOSYS, NULL, "NOTICE", NULL,
		"Status file saved",
		NULL);
	}
      
      /* incrementa il contatore degli step temporali */
      ++Oparams.curStep;
    }


  Oparams.curStep--;

  /* <-------------------------------------------- MAIN LOOP END(FATHER) */
  tempo = time(NULL);
  sprintf(msgStrA, "FINAL TIME: %s", ctime (&tempo));

  mdMsg(ALL, NOSYS, "End", "NOTICE", NULL,
	msgStrA,
	NULL);
  
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /* <----------------------------------------------------- SAVE ENDFILE */

  if (mgl_mode==2)
    {
      UPDATE_SYSTEM;
      ADJUST_LASTCOL;
      saveCorAscii(); 
    }
  else
    {
      if (OprogStatus.endFormat == 0 || OprogStatus.endFormat == 2)
	{
	  UPDATE_SYSTEM;
	  ADJUST_LASTCOL;
	  saveCoord(absTmpHD(OprogStatus.endfile)); 
	}
      if (OprogStatus.endFormat == 1 || OprogStatus.endFormat == 2)
	{
	  UPDATE_SYSTEM;
	  ADJUST_LASTCOL;
	  saveCorAscii(); 
	  printf("salvo le coordinate in formato ascii\n");
	}
    }
  /* take coords from restore files on HD and write them to the
	 coordinate file specified as argument */	
  
  /* <------------------------------------------------- DELETE TMP FILES */
  
  /* delete files STATUS_FILE_NAME + '0' and + '1' */
  delDoubleHD(STATUS_FILE_NAME);      
  /* This file is deleted before because if the system crash occurs
     just after this unlinking, chkStatus will not try to restart 
     simulation (not being the STATUS_FILE). */
  
  /* se e' stato regolarmente creato e scritto 
     il file chiamato Oparams.endfile allora cancella i 
     files  temporanei di ripristino */
  
  /* delete restore files on the Hard Disk, that is 
     BAK_FILE_NAME + '0' and + '1' */
  delDoubleHD(BAK_FILE_NAME);    
  
  /* delete restore files on tape */ 
  if (OprogStatus.tapeTimes !=0) delDoubleTape(BAK_FILE_NAME);
  
  /* NOTA: se accadono degli errori durante l'unlink i files non vengono
     cancellati !!! */
  
  delTmpHD(CF); /* remove check file if exist, because 
		   at this point it is no more useful */
  
  /* <------------------------------------------------REMOVE SHM END SEM */
  
  /* all shared memory segment are removed by delAllShm 
     (see endFather in mdinit.c file) */  
  
  mdMsg(ALL, NOSYS, "End", "NOTICE", NULL,
	"Simulation successfully terminated",
	NULL);
} /* <------------------------------------------------------- END FATHER */
#ifdef MD_BILOG
#ifdef MDLLINT
long long int *bilog_arr;
#else
int *bilog_arr;
#endif
int compare_int(const void *a, const void *b)
{
  int n1, n2;
#ifdef MDLLINT
  n1 = *((long long int*)a);
  n2 = *((long long int*)b);
#else
  n1 = *((int*)a);
  n2 = *((int*)b);
#endif
  if (n1 < n2)
    return -1;
  else if (n1 == n2)
    return 0;
  else
    return 1;
}
#ifdef MDLLINT
#define MD_INT long long int
#define MD_INTP long long int *
#define MD_INTSTR "%lld"
#else
#define MD_INT int
#define MD_INTP int *
#define MD_INTSTR "%d"
#endif
inline MD_INT bilog_step(double fstpsw, double fstps)
{
  return ((MD_INT) (rint(fstps) - 1.0 + rint(fstpsw)));
}
void build_bilog_arr(void)
{
  double fstpsw=1.0, fstps=1.0;
  int maxpts, maxptw=80, maxpt=80; 
  int exist, end=0, npts, nptst, nptstw, i; 
  int tots;
  tots = Oparams.totStep?Oparams.totStep:2000000; 
  maxptw = (int)(log10((double)tots)/log10(OprogStatus.basew));
  maxpt =  1 + (int)(log10((double)tots)/log10(OprogStatus.base));
  maxpts = maxptw*maxpt;
  printf("Building bilog array maxptw: %d maxpt: %d base: %f basew: %f\n", maxptw, maxpt,
	 OprogStatus.base, OprogStatus.basew);
  bilog_arr = malloc(sizeof(MD_INT)*maxpts);
  npts = 0;
  for (i=0; i < maxpts; i++)
    bilog_arr[i] = 0;
  nptstw = 0;
  while (!end)
    {
      nptst = 0;
      fstps = 1.0;
      while (!end)
	{
	  exist = 0;
	  for (i=0; i < npts; i++)
	    {
		if (bilog_arr[i]==bilog_step(fstpsw, fstps))
	  	  {
	  	    exist = 1;
	  	    break;
	  	  }
	    }
	  if (!exist)
	    {
      	      bilog_arr[npts] = bilog_step(fstpsw, fstps);
	      npts++;
	      if (npts == maxpts)
		{
		  end = 1;
		  break;
		}
	    }
	  fstps = fstps * OprogStatus.base;
	  nptst++; 
	  if (bilog_step(fstpsw, fstps) >= tots)
	    break;
  	}
      fstpsw = fstpsw * OprogStatus.basew;
      nptstw++;
      if ((MD_INT) rint(fstpsw) >= tots)
	end=1;
      /*printf("tots: %d nptstw: %d nptst:%d npts: %d fstpw: %f fstps: %f\n",tots, nptstw, nptst, npts, fstpsw, fstps);*/
    }
  qsort((void*)bilog_arr, npts, sizeof(MD_INT), compare_int);
#if 0
  for (i=0; i < npts; i++)
    {
      printf("bilog_arr[%d]=" MD_INTSTR "\n", i, bilog_arr[i]);
    }
  exit(-1);
#endif
}
#endif
extern void scanFile(char* argom);
/* =============================== >>> MAIN <<< ============================*/
extern void init_rng(int mdseed, int mpi, int my_rank);
#if !defined(MD_STANDALONE) || defined(MC_SIMUL) 
int main(int argc, char *argv[])
{
#if defined(MPI)  
  int mpiStatus, i;
  MPI_Status status;
  char txt[1];
#endif
  /* check if program is already running 
     19/3/99 REMOVED this check */
  
  /*if (progExist("mdsimul"))
    {
    printf("Program already running!\n");
    exit(-1);
    }*/      

#if defined(MPI) 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);

  init_rng(-1, 1, my_rank);
#if 0
#if defined(MD_RAND48)
  srand48(my_rank*((long int)time(NULL)));
#elif defined(MD_RANDOM)
#if 1
  for (k=0; k < 256; k++)
    ss[k] = (char) ((int)time(NULL))%256;
  initstate(my_rank*((int)time(NULL)), ss, 256);
#else
  srandom(my_rank*((int)time(NULL)));
#endif
#else
  srand(my_rank*((int)time(NULL)));
#endif

  //srand(my_rank*((int)time(NULL)));
#endif
  equilibrated = malloc(sizeof(int)*numOfProcs);
  for (i=0; i < numOfProcs; i++)
    equilibrated[i] = 0;
#endif
  MD_EXT_INIT(&argc, &argv);
  // printf("PID CHILD: %d\n",pid);
  /* ------------------>>> function to call at exit <<<-------------- */ 
#if !defined(CRAY) && !defined(ORIGIN) && !defined(MD_MAC)
  on_exit(endJob, NULL);
#endif
  /* ---------------------------------------------------------------- */
  
  /* --------------------->>> signal handling<<<--------------------- */
  signal(SIGINT, SuspendJob); /* CtrlC */
  signal(SIGTSTP, SuspendJob);/* CtrlZ */
  /* ------------------------------------------------------------------*/

#ifdef MPI
#if defined(PRINTER_PROC)
  printerRank = numOfProcs - 1;
  numOfProcs -= 1;
  if (my_rank == printerRank)
    {
      _printing();
      /* Termina regolarmente */
      exit(0);
    }
#endif
#endif

  //nice(-20);  /* this works only if the superuser starts the simulation */
   
  if ( argsMd(argc, argv) == 1) /* args parse commandline arguments  */
    /* args = 1 means:  'begin  a new simulation' */
    {
      Newsimul(paramFile);
    }
  else /* otherwise continue a previuosly interrupted simulation */
    {
      /* inside Continue() BAK is initialized in a way that the next 
	 saving is in a restore file different to the file used to 
	 restart simulation.
	 This is because if another system crash occurs too closely 
	 it can damage the unique good file and the simulation 
	 couldn't restart again. 
	 Similarly STA is initialized so that the next status file 
	 saved is not the corrupted one. */
      Continue();
      /* increments the step counter, so that the next step to save
	 is not this one ( see main loop ) */
      ++Oparams.curStep;
    }       
  usrInitAft(); /* user initialization */
#ifdef MD_BILOG
  build_bilog_arr();
#endif  
  delTmpHD(CF); /* Delete the check file */
  
  ENDSIM = 0;
  /* This variable could be set to 1 inside loop to ende the simulation */

  
  commMD();

#ifdef MD_HE_PARALL
  MD_EXT_END();
#endif
#if defined(MPI)
#if defined(PRINTER_PROC)
  MPI_Ssend("fine", 5, MPI_CHAR, 
	   printerRank, MD_MPI_END, MPI_COMM_WORLD);

#endif
  /* Termina regolarmente */
  exit(0);
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Finalize();
#endif 
  return 0;
}
#endif
