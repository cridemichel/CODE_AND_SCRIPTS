#undef MAIN
#include<mdsimul.h>

/* GENERAL NOTE VERY IMPORTANT
   NEVER USE STRING OPERATORS appSw, abs... MORE THEN ONE TIMES AS ARGS OF 
   sum(char *, ...), BECAUSE ALL THE OPERATOR PUT THE RESULT IN THE SAME 
   MEMORY AREA */
#ifdef MDLLINT
long long int logBlock;
#else
int logBlock;
#endif
#ifdef MPI
extern int my_rank;
extern int numOfProcs;
#endif

/* ======================== >>> absFileNameHD <<< ===========================*/
char* absTmpAsciiHD(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory MD_HD_TMP with 'fileSrc'
     and return the resulting.
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];

  strcpy(absFile, OprogStatus.tmpPath);

  /* MD_TMP is the directory on Hard Diskc ontaining simulation temporary 
     files */
  
  strcat(absFile, fileSrc);  
  return absFile;
}

/* ======================== >>> absFileNameHD <<< ===========================*/
char* absTmpHD(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory MD_HD_TMP with 'fileSrc'
     and return the resulting.
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];
  
  strcpy(absFile, MD_HD_TMP); 

  /* MD_TMP is the directory on Hard Diskc ontaining simulation temporary 
     files */
  
  strcat(absFile, fileSrc);  
  return absFile;
}

/* ======================== >>> absFileNameTape <<< =========================*/
char* absTmpTape(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory MD_TAPE_TMP with 'fileSrc'
     and return the resulting string. 
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];
  strcpy(absFile, MD_TAPE_TMP); 
  /* MD_TAPE_TMP is the directory on the tape containing simulation temporary 
     files */
  strcat(absFile, fileSrc);  
  return absFile;
}

/* ======================== >>> absFileMisTape <<< =========================*/
char* absMisTape(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory on Tape 
     MD_TAPE_MIS with 'fileSrc' and return resulting string.
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];
  
  strcpy(absFile, MD_TAPE_MIS); 
  /* MD_TAPE_MIS is the directory on the tape containing 
     measures */
  strcat(absFile, fileSrc);  
  return absFile;
}

/* 16/1/1998 ADD: absXvaTape() and absXvaHD() */
/* ======================== >>> absFileMisTape <<< =========================*/
char* absXvaTape(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory on Tape 
     MD_TAPE_XVA with 'fileSrc' and return resulting string.
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];
  
  strcpy(absFile, MD_TAPE_XVA); 
  /* MD_TAPE_MIS is the directory on the tape containing 
     xva files */
  strcat(absFile, fileSrc);  
  return absFile;
}

/* ======================== >>> absFileMisTape <<< =========================*/
char* absXvaHD(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory on Tape 
     MD_HD_XVA with 'fileSrc' and return resulting string.
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];
  
  strcpy(absFile, MD_HD_XVA); 
  /* MD_HD_XVA is the directory on the Hard Disk containing 
     xva files */
  strcat(absFile, fileSrc);  
  return absFile;
}

/* ======================== >>> absFileMisTape <<< =========================*/
char* absMisHD(const char* fileSrc)
{
  /* Concatenate the simulation tenporary directory MD_HD_MIS with 'fileSrc'
     and return the resulting string.
     You can think to this like an operator that returns absolute files
     taking relative ones */
  static char absFile[NAME_LENGTH];
  if (!strcmp(OprogStatus.misPath,"*"))
    strcpy(absFile, MD_HD_MIS); 
  else
    strcpy(absFile, OprogStatus.misPath); 
  
  /* MD_HD_MIS is the directory on the hard disk  containing 
     simulation measures files  */
  strcat(absFile, fileSrc);  
  return absFile;
}

/* =========================== >>> appSw <<< =============================*/
char* appSw(char* srcStr, const unsigned char which)
{
  /* DESCRIPTION:
     append to 'fileName' the value of 'which', and returns 
     the resulting string 
     You can think to this like an operator that returns a files  with the
     switch number appended. */
  char tmpStr[10];
  static char resultStr[NAME_LENGTH];
  strcpy(resultStr, srcStr);
  sprintf(tmpStr, "%d", which);
  strcat(resultStr, tmpStr);
  return resultStr;
}

/* ============================= >>> MPIunlink <<< ======================= */
void MPIunlink(char * stri)
{
#ifdef MPI
  char fileName[NAME_LENGTH];
  sprintf(fileName, "%s_R%d", stri, my_rank);
  unlink(fileName);
#else 
  unlink(stri);
#endif  
}


/* ===================== >>> doubleSaveMeasure <<< ==========================*/
unsigned char sw(unsigned char onOff)
{
  /* DESCRIPTION:
     If onOff = 0 return 1 
     If onOff = 1 retrun 0, that is switch onOff */
  return ~(onOff) & 1;
}

/* ======================== >>> delBakFilesHD <<< ==========================*/
void delDoubleHD(char *fileName)
{
  /* DESCRIPTION:
     If 'fileName' is 'pippo' this procedure deletes 'pippo1' and 'pippo2'
     on the Hard Disk temporary directory 
  */
  MPIunlink(  appSw( absTmpHD(fileName), 0)  );
  MPIunlink(  appSw( absTmpHD(fileName), 1)  );
}


/* ======================== >>> delDoubleTape <<< =========================*/
void delDoubleTape(char *fileName)
{
  /* DESCRIPTION:
     If fileName is 'pippo' this procedure deletes 'pippo1' and 'pippo2'
     on tape temporary directory */
  MPIunlink(  appSw( absTmpTape(fileName), 0)  );
  MPIunlink(  appSw( absTmpTape(fileName), 1)  );
}

/* ========================== >>> delTmpHD <<< =============================*/
void delTmpHD(char* relName)
{
  /* DESCRIPTION:
     This procedure deletes the file 'fileName'in simulation temporary 
     directory MD_HD_TMP (HARD DISK) */
  MPIunlink( absTmpHD(relName) );
}

/* ========================== >>> delTmpTape <<< ============================*/
void delTmpTape(char* relName)
{
  /* DESCRIPTION:
     This procedure deletes the file 'fileName'in simulation temporary 
     directory MD_TAPE_TMP (TAPE) */
  MPIunlink( absTmpTape(relName) );
}

/* ============================ >>> sum  <<< =============================== */
char* sum(char *str, ...)
{
  /* DESCRIPTION:
     concatenate the string passed as arguments and return the resulting 
     string */
  char* sptr;              /* pointer to the next arg, returned by va_arg()*/
  va_list ap;
  static char retStr[255]; /* resulting string */
  va_start(ap, str);
  
  strcpy(retStr, str);     /* at least str must be present */
  /* NOTE: the list of string must end  with a NULL !!!!! */
  while ( (sptr = va_arg(ap, char*)) != NULL )
    { 
      strcat(retStr, sptr);/* concatenate variuos args */
    }
  va_end(ap);
  return retStr;
}

/* ========================= >>> simulExist <<< ===========================*/
int progExist(char *name)
{
  /* Check if the program named 'name' is already running, 
     using the special proc files in the /proc directory 
     'name' is the name used by command line to execute the program 
     >>> RETURN VALUE:
     If the program exist return true (!=0), false (0) otherwise. */
  
  /* Every process has a directory in /proc, the name of which is its PID.
     If for ex. the process pid is 800, then we have the following files:
     '/proc/800/cmdline' : command used to execute the process 
     '/proc/800/mem'
     '/proc/800/status'
     etc.
     'status' is a file of this kind:
     Name: <program name>
     State: <program state (Running,sleeping,etc.) >
     .
     .
  */
  
  FILE* pfs;                  /* process file descriptor */
  char tstri1[NAME_LENGTH];   /* fictitiuos string */ 
  char procName[NAME_LENGTH]; /* name of process found in the status file */
  struct dirent *POdirent;    /* structure returned by readdir() */
  int procPid,curPid;         /* Pid of process found in status file (procPid) 
				 and pid of actual process (curPid) */
  DIR* procDir;               /* pointer to e directory returned by opendir */
  procDir = opendir("/proc"); /* open the dir /prov */
  
  curPid = getpid();
  //printf("Actual pid: %d\n", curPid);
  if (procDir == NULL)
    {
      perror("opendir");
      printf("I can't open proc directory");
      exit(-1);
    }
  while ( (POdirent = readdir(procDir)) ) /* read the next directory item 
					     ( NULL if EOF ) */
    {
      //printf("inside while: %s\n", POdirent->d_name);
      /* if the directory is a number then it is a process directory,
	 with inside a status file to check */
      if (atoi(POdirent->d_name) > 0) 
	{
	  strcpy(tstri1, "/proc/");
	  strcat(tstri1, POdirent->d_name);
	  strcat(tstri1, "/status");
	  /* tstri is now the full name (with path) of the status proc file,
	     for ex. /proc/820/status if 820 is the process PID */ 
	  if ( (pfs = fopen(tstri1, "r")) != NULL )
	  {  
	    fscanf(pfs,"Name: %s\n", procName);
	    
	    //fgets(tstri1,65535,pfs); /* skip the State line */
	    fscanf(pfs, "State: %s %s\n", tstri1, tstri1);
	    fscanf(pfs, "Pid: %s\n", tstri1);
	    procPid = atoi(tstri1);
	    /* take the first three lines of status proc file chars and put 
	       them in tstri, we use only Name and Pid */
	
	    /* The question is: exist a process named 'name' with different 
	       pid of this process? if Y then the process is already 
	       running */
	    if ( (!strcmp(procName,name)) && 
		 (procPid != curPid) )
	      {
		
		closedir(procDir);
		return 1;             /* 1 = program is already running */
	      }
	    fclose(pfs);
	  }
	}
    }
  closedir(procDir);  
  return 0;                     /* 0 = program does not already exist */
}

/* ============================== >>> copy <<< ==============================*/
int copy(char *src, char* dest) /* src and dest must be absolute */
{
  /* DESCRIPTION:
     This function copy 'src' to 'dest' and if 'dest' exists truncate it or
     create it otherwise */
  
  int i;           /* counter */
  int srcd, destd; /*descriptors associated to files */
  char buf;        /* one byte buffer to store byte read by 'read' system call
		    */ 
  int readSrc;     /* return value of reading of src file */
  if ( (srcd = open(src, O_RDONLY)) == -1 )
    {
      return -1;  /* -1 = ERROR */
    }
  if ( (destd = creat(dest, 0666)) == -1)
    {
      return -1;
    }
  /* not very efficient, but doesn't matter ... */
  for (i=1; (readSrc = read(srcd, &buf, 1)) > 0; ++i)
    {
      if ( write(destd, &buf, 1) == -1 )
	{
	  return -1;
	}
    } 
  
  /* If it was exit from loop because of error => complain and exit */ 
  if (readSrc == -1)
    {  
      return -1;
    }
  
  /* close both files */
  if ( close(destd) == -1) 
    {  
      return -1;
    }
  
  /* if an  error closing src doesn't matter */
  close(srcd);
  sync();
  return 0;   /* 0 means SUCCESS */
}
#ifdef MD_BILOG
extern int *bilog_arr;
#endif
/* =========================== >>> chkBakAsciiSteps <<< ==================*/
int chkBakAsciiSteps(void)
{
  double base; 
#ifdef MDLLINT
  long long int oldfstps, retval;// *timeout; 
  long long int cslb;
#else
  int oldfstps, retval;// *timeout; 
  int cslb;
#endif
  base = OprogStatus.base;
  //timeout = &OprogStatus.timeout;
  switch (OprogStatus.bakSaveMode){ 
  case 0:
    /* Linear saving of xva file */
    if (OprogStatus.bakStepsAscii == 0)
      return 0;
    retval =  ((Oparams.curStep % OprogStatus.bakStepsAscii) == 0);
    break;
  case 1:
    /* Semilog saving */
    //*timeout = OprogStatus.KK * logBlock + OprogStatus.fastSteps;
    if (OprogStatus.NN == 0)
      return 0;

    //printf("logBlock:%d\n", logBlock);
    retval = 0;
    cslb = Oparams.curStep % logBlock;
  #ifdef MDLLINT
    if ( (cslb == ((long long int)rint(OprogStatus.fstps))))
#else
    if ( (cslb == ((int)rint(OprogStatus.fstps))))
#endif
     {
       if (cslb == 0)
	 retval = 0;
       else
	 retval = 1;
	/*printf("[%d]  cs mod lb: %d fsteps: %f logblock: %d\n", Oparams.curStep,
	  Oparams.curStep % logBlock,
	  OprogStatus.fstps, logBlock);*/

	//OprogStatus.fstps = (int) (base * 
	//((double) OprogStatus.fstps));
#ifdef MDLLINT
	oldfstps = (long long int)rint(OprogStatus.fstps);
	while (oldfstps == (long long int)rint(OprogStatus.fstps))
	  OprogStatus.fstps = base * OprogStatus.fstps;
#else
	oldfstps = (int)rint(OprogStatus.fstps);
	while (oldfstps == (int)rint(OprogStatus.fstps))
	  OprogStatus.fstps = base * OprogStatus.fstps;
#endif
#ifdef MPI
	if(my_rank == 0)
#endif
#ifdef MDLLINT
	printf("[%lld] fstps: %.6f\n", Oparams.curStep, OprogStatus.fstps);
	if ( ((long long int)rint(OprogStatus.fstps)) == logBlock )
	  {
	    OprogStatus.fstps = 1;
	    printf("[%lld] fstps: %.6f\n\n", Oparams.curStep, OprogStatus.fstps);
	  }
#else
	 printf("[%d] fstps: %.6f\n", Oparams.curStep, OprogStatus.fstps);
	 if ( ((int)rint(OprogStatus.fstps)) == logBlock )
	   {
	     OprogStatus.fstps = 1;
	     printf("[%d] fstps: %.6f\n\n", Oparams.curStep, OprogStatus.fstps);
	  }
#endif
     }
    break;
  case 2:
#if !defined(MD_BILOG)
    /* Bilog saving */
    printf("add necessary fields to OprogStatus struct and define macro MD_BILOG!!!\n");
    exit(-1);
#else
    /* OprogStatus.fstps viene usata per memorizzare l'indice nell'array bilog_arr
     * del prossimo step a cui si salvare */
    if (bilog_arr[OprogStatus.lastbilogsaved] && Oparams.curStep == bilog_arr[OprogStatus.lastbilogsaved])
      {
	OprogStatus.lastbilogsaved++;
	return 1;
      }
    return 0;
#endif
    break;
  default:
    printf("ERROR: Invalid save mode!!!!\n");
    exit(-1);
  }
  return retval;
}


/* =========================== >>> chkXvaStep <<< ========================== */
int chkXvaSteps(void)
{
  double base; 
#ifdef MDLLINT
  long long int cslb, oldfstps, retval;// *timeout; 
#else
  int cslb, oldfstps, retval;// *timeout; 
#endif
  base = OprogStatus.base;
  //timeout = &OprogStatus.timeout;
  switch (OprogStatus.xvaSaveMode){ 
  case 0:
    /* Linear saving of xva file */
    if (OprogStatus.xvaSteps == 0)
      return 0;
    retval =  ((Oparams.curStep % OprogStatus.xvaSteps) == 0);
    break;
  case 1:
    /* Semilog saving */
    //*timeout = OprogStatus.KK * logBlock + OprogStatus.fastSteps;
    if (OprogStatus.NN == 0)
      return 0;
    cslb = Oparams.curStep % logBlock;
    retval = 0;
#ifdef MDLLINT
    if ( (cslb == ((long long int)OprogStatus.fstps)) || (cslb == 0))
#else
    if ( (cslb == ((int)OprogStatus.fstps)) || (cslb == 0))
#endif
      {
	retval = 1;
	
	//printf("fstps: %.6f logBlock: %d\n", OprogStatus.fstps, logBlock);

	//printf("cs:%d  cs mod lb: %d fsteps: %f logblock: %d\n", Oparams.curStep,
	//	 Oparams.curStep % logBlaock,
	//	 OprogStatus.fstps, logBlock);

	//OprogStatus.fstps = (int) (base * 
	//((double) OprogStatus.fstps));
#ifdef MDLLINT
	oldfstps = (long long int)OprogStatus.fstps;
	while (oldfstps == (long long int)OprogStatus.fstps)
	  OprogStatus.fstps = base * OprogStatus.fstps;
	if ( ((long long int)OprogStatus.fstps) > logBlock )
	  {
	    OprogStatus.fstps = 1;
	  }
#else
	oldfstps = (int)OprogStatus.fstps;
	while (oldfstps == (int)OprogStatus.fstps)
	  OprogStatus.fstps = base * OprogStatus.fstps;
	if ( ((int)OprogStatus.fstps) > logBlock )
	  {
	    OprogStatus.fstps = 1;
	  }
#endif
      }
    break;
  case 2:
    /* Bilog saving */
#if !defined(MD_BILOG)
    /* Bilog saving */
    printf("add necessary fields to OprogStatus struct and define macro MD_BILOG!!!\n");
    exit(-1);
#else
    /* OprogStatus.fstps viene usata per memorizzare l'indice nell'array bilog_arr
     * del prossimo step a cui si salvare */
    if (bilog_arr[OprogStatus.lastbilogsaved] && Oparams.curStep == bilog_arr[OprogStatus.lastbilogsaved])
      {
	OprogStatus.lastbilogsaved++;
	return 1;
      }
    return 0;

#endif
    break;
  default:
    printf("ERROR: Invalid save mode!!!!\n");
    exit(-1);
  }
  return retval;
}


