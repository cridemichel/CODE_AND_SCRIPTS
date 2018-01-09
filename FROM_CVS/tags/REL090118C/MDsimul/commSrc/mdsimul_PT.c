/*      $Id: mdsimul_PT.c,v 1.1.1.1 2004-02-28 10:53:56 demichel Exp $     */
#ifndef lint
static char vcid[] = "$Id: mdsimul_PT.c,v 1.1.1.1 2004-02-28 10:53:56 demichel Exp $";
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


struct simStat OsimStat;	     /* global istance of sim_status 
					structure */
/* ----------------------------------------------------------------------*/
int SEGSIZE;         /* arrays dimension in bytes 
				 ( set in NewSimul() or Continue(),
				 see mdinit.c file )*/
unsigned char BAK, STA;       /* each is a switch (0/1) which determine 
				 the next file to save to on HD 
				 ("double unsigned") */
char BAKT;           /* switch used for tape (analogous to BAK for 
			Tape)*/ 
				 
extern char TXT[MSG_LEN];
/* shared array of two integers for active lock/unlock (semaphores) */

/* log stream */
FILE* output;

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

extern void parallel_tempering_init(void);
extern int doubleSaveMeasure(int PN, int misNum, int* times, int tapeTimes);

extern void saveBakAscii(char *);
extern void saveCorAscii(void);

/* ============================= >>> commMD <<< ======================= */
void commMD(void)
{
  int i;   
  int stato;            /* set by wait */
  time_t tempo;         /* returned by time(NULL) */
  int msteps;
  //ProcSync0();  
  
#ifdef MPI
  mdMsg(ALL, NOSYS, "Begin", "NOTICE", NULL,
	    "MPI simulation",
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
      //mdPrintf(ALL, "fin qui OK!\n", NULL);

      //ProcSync0(); /* <-----NOT NECESSARY !!!!--------SYNCRONIZATION */
      
      /* <------------------------------------------- MEASURES (FATHER)
	 ACTUALLY NOT PARALLELIZED !!! */
      calc(); 
      
      /* --------------------------------------------------------------- */
      
      //ProcSync0(); /* <------------------------------ SYNCRONIZATION */
      
      /* <--------------------------------------------------- SAVE MEASURE 
	 Every 'OprogStatus.measSteps[i]' steps save measure structure 'i' 
	 into the file Ofilenames.data_files[i] (see TECH_INFO file) */
      
      for (i = 0; Omeasure[i].buf != NULL; ++i) 
	/* NULL in the buf pointer of measure sruct means end of list */
	{
	  msteps = abs(OprogStatus.measSteps[i]);
	  
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
	      //mdPrintf(ALL, "misura salvata!!\n", NULL);
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
	{
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
      /* ADDED 02/04/2001: restart file in ascii format */
      if (chkBakAsciiSteps())
	{
	  
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
  

  /* <----------------------------------------------------- SAVE ENDFILE */
     
  if (OprogStatus.endFormat == 0 || OprogStatus.endFormat == 2)
    {
      saveCoord(absTmpHD(OprogStatus.endfile)); 
    }
  if (OprogStatus.endFormat == 1 || OprogStatus.endFormat == 2)
    {
      saveCorAscii(); 
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


/* =============================== >>> MAIN <<< ============================*/
void main(int argc, char *argv[])
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
  srand(((int)time(NULL))+my_rank);
  //sprintf(TXT, "rank:%d\n", my_rank);
#endif
    // printf("PID CHILD: %d\n",pid);
  /* ------------------>>> function to call at exit <<<-------------- */ 
#if defined(ORIGIN) || defined(CRAY)
  //atexit(endJob);
#else
  on_exit(endJob, NULL);
#endif
  /* ---------------------------------------------------------------- */
  
  /* --------------------->>> signal handling<<<--------------------- */
  signal(SIGINT, SuspendJob); /* CtrlC */
  signal(SIGTSTP, SuspendJob);/* CtrlZ */
  /* ------------------------------------------------------------------*/

#ifdef MPI
#if defined(PRINTER_PROC)
  printf("PRINTER PROC ACTIVE!!!\n");
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
      Newsimul(paramFile);/* parsFile is the params file passed as argument */
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
  /* ADDED 15/09/2000: Parallel Tempering Initialization */
  //parallel_tempering_init();
  
  usrInitAft(); /* user initialization */
  
  delTmpHD(CF); /* Delete the check file */
  
  ENDSIM = 0;
  /* This variable could be set to 1 inside loop to end the simulation */

  
  commMD();

#if defined(MPI)
#if defined(PRINTER_PROC)
  MPI_Ssend("fine", 5, MPI_CHAR, 
	   printerRank, MD_MPI_END, MPI_COMM_WORLD);

#endif
  /* Termina regolarmente */
#if defined(ORIGIN) || defined(CRAY)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else  
  exit(0);
#endif
#endif 
}
