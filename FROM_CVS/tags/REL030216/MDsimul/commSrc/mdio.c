#undef MAIN
#include<mdsimul.h>
/* string used to send  messages */

char msgStrA[MSG_LEN],msgStrB[MSG_LEN],msgStrC[MSG_LEN]; 
FILE* output;

#if defined(MPI)
int printerRank;
int my_rank, numOfProcs;
MPI_Status status;
#endif
/* ============================ >>> openMPI <<< ============================ */
int openNewMPI(char* fn, int how, int perm)
{
  /* DESCRIPTION:
     Open an existant file.
     if MPI macro is defined open a file name adding _R<rank>, that is 
     if <fn> is the file name passed by fn, then the file opened by this 
     function is: <fn>_R<rank>, where <rank> is the rank of the process */
#ifdef MPI
  char fileName[MAX_LENGTH];
  char tmpStr[MAX_LENGTH];
  int fd;
  strcpy(tmpStr, fn);
  sprintf(fileName, "%s_R%d", tmpStr, my_rank);
  fd = open(fileName, how, perm);
  return fd;
#else 
  return open(fn, how, perm);
#endif
}


/* ============================ >>> openMPI <<< ============================ */
int openMPI(char* fn, int how)
{
  /* DESCRIPTION:
     if MPI macro is defined open a file name adding _R<rank>, that is 
     if <fn> is the file name passed by fn, then the file opened by this 
     function is: <fn>_R<rank>, where <rank> is the rank of the process */
#ifdef MPI
  char fileName[MAX_LENGTH];
  char tmpStr[MAX_LENGTH];
  int fd;
  strcpy(tmpStr, fn);
  sprintf(fileName, "%s_R%d", tmpStr, my_rank);
  fd = open(fileName, how);
  return fd;
#else 
  return open(fn, how);
#endif
}

/* ============================ >>> openMPI <<< ============================ */
FILE* fopenMPI(char* fn, char* how)
{
  /* DESCRIPTION:
     if MPI macro is defined open a file name adding _R<rank>, that is 
     if <fn> is the file name passed by fn, then the file opened by this 
     function is: <fn>_R<rank>, where <rank> is the rank of the process */
#ifdef MPI
  char fileName[MAX_LENGTH];
  char tmpStr[MAX_LENGTH];
  FILE* fs;
  strcpy(tmpStr, fn);
  sprintf(fileName, "%s_R%d", tmpStr, my_rank);
  fs = fopen(fileName, how);
  return fs;
#else 
  return fopen(fn, how);
#endif
}
#if 0
/* valgrind non supporta creat e quindi bisogna usare questo wrapper */
#define creat(A,B) open(A,O_CREAT|O_WRONLY|O_TRUNC,B)
#endif
/* ============================ >>> openMPI <<< ============================ */
int creatMPI(char* fn, int how)
{
  /* DESCRIPTION:
     if MPI macro is defined open a file name adding _R<rank>, that is 
     if <fn> is the file name passed by fn, then the file opened by this 
     function is: <fn>_R<rank>, where <rank> is the rank of the process */
#ifdef MPI
  char fileName[MAX_LENGTH];
  char tmpStr[MAX_LENGTH];
  int fd;
  strcpy(tmpStr, fn);
  sprintf(fileName, "%s_R%d", tmpStr, my_rank);
  fd = creat(fileName, how);
  return fd;
#else 
  return creat(fn, how);
#endif
}

/* ============================= >>> mdCreat <<< ============================*/
int mdCreat(char* fileName, char* when, char* errMsg, int mode) 
{
  /* DESCRIPTION:
     'Filename' is the name of the file to create and mode could be EXIT or 
     CONT, that are predefined macros.
     If mode == EXIT, if an error occurs then it complains and exit, instead
     if mode == CONT it complains only. 
     errMsg is an extra error message it prints when an error occurs
     when = NULL => print the actual step 
          = "..." => print the string indicating the phase of the simulation */
  int fd;
  if (  ( fd = creatMPI(fileName, 0666) ) 
       == -1  )
    {
      sprintf(msgStrA,"Can't open file %s for writing",
	      fileName);
      mdMsg(ALL, errno, when, "ERROR", "open",
	    msgStrA,
	    errMsg,
	    NULL);
      
      if (mode == EXIT) /* exit if the user want this */
	{
	  exit(-1);
	}
    }
  return fd; /* return the file descriptor or -1 if an error occure and 
		mode == CONT */ 
}

/* ========================= >>> ceratWithHead <<< =========================*/
int creatWithHead(char* fileName, char* when, char* errMsg, int mode,
		   int headSize, void* header)
{
 /* DESCRIPTION:
    create a file fileName, putting the header '*header' at the begin,
    you must supply header length as third argument and a pointer to it in as 
    second argument.
    - when can NULL = print current step on error or a string that specify
      the phase of simulation.
    - errMsg is an extra message to print on error, actually not used.
    - mode = EXIT => complain and exit on error.
        "  = CONT => complain but don't exit on error.
    
 */
  int fd;
  if (  ( fd = creatMPI(fileName,
		     0666) ) == -1  )
    {
      sprintf(msgStrA, "Unable to create %s file", fileName);
      mdMsg(ALL, errno, when, "ERROR", "creat",
	    msgStrA,
	    NULL);
      if (mode == EXIT)
        {
	  exit(-1);
	}
      return -1;
    }
  /* header of mesure file consist of the interval in steps between
     two measure savings and in the size of each measure in bytes */
  else if ( write(fd, header, 
		  headSize) == -1 )
    {
      sprintf(msgStrA, "Unable to write the %s file", fileName);
      mdMsg(ALL, errno, when, "ERROR", "write",
	    msgStrA,
	    NULL);
      if (mode == EXIT)
	{
	  exit(-1);
	}
      return -1;
    }
  else if (close(fd) == -1)
    {
      sprintf(msgStrA, "Unable to close %s file", fileName);
      mdMsg(ALL, errno, when, "ERROR", "close",
	    msgStrA,
	    NULL);
      if (mode == EXIT)
	{
	  exit(-1);
	}    
      return -1;
    }

  //printf("savesteps: %d\n", ((struct measHead*) header)->saveSteps);
    
  return fd; /* actually the returned value is not used in the program */
 }

/* =========================== >>> mdOpen <<< ===============================*/
int mdOpen( char* fileName, char* when, char* errMsg, int mode)
{
  /* DESCRIPTION:
     mode is the usual open mode except that you can also specify the EXIT
     flag or-ed with the usual flags ( actually EXIT = 0900000000 )*/
  int fd;
  int exitFlag;

  exitFlag = mode & EXIT; /* exitFlag = EXIT or 0 in this way */
  mode &= ~EXIT; /* set to zero the EXIT flag not used by open system call */ 

  
  if (  ( fd = openNewMPI(fileName, mode, 0666) ) 
       == -1  )
    {
      sprintf(msgStrA,"Can't open file %s for writing",
	      fileName);
      mdMsg(ALL, errno, when, "ERROR", "open",
	    msgStrA,
	    errMsg,
	    NULL);
      
      if (exitFlag == EXIT) /* exit if the user want this */
	{
	  exit(-1);
	}
    }
  return fd; /* return the file descriptor or -1 if an error occure and 
		mode == CONT */ 

} 

/* ============================ >>> readArr <<< =============================*/
int  mdRead(int fdes, char* when, char* errMsg, int mode, 
	    int size, void *pointer)
{
  /* 
     DESCRIPTION:
     - errMsg is an extra error message print when an error occurs
     - when could be NULL => it print the actual step or a string that
       must indicate the phase of the simulation. 
     - size is the nnumber of bytes to read an pointer is the location where
       put the datas.
     - If mode == EXIT, if an error occurs then it complains and exit, instead
       if mode == CONT it complains only.
  */
  
  int br; /* bytes read by 'read' system call */ 
  /* legge size bytes dallo stream fdes , e li mette nel buffer puntato
     da 'pointer' */
  
  /* controlla il buon esito della lettura, se fallisce esce, poiche' 
     la simulazione sarebbe compromessa */
  br = read(fdes, pointer, size); 
  
  if (br == -1) 
    {
      mdMsg(ALL, errno, when, "ERROR", "read",
	    errMsg,
	    NULL);
      if ( mode == EXIT )
	{
	  exit(-1);
	}
      return -1;			/* -1 = ERROR */
    }
  else if (br < size) 
    {
      mdMsg(ALL, NOSYS, when, "ERROR", NULL,
	    errMsg, 
	    "Too few bytes read",
	    NULL);
      if (mode == EXIT)
	{
	  exit(-1);
	}
      
      return -1;               /* -1 = ERROR */
    }
  return 0;		       /* 0 = OK */
}

/* ============================= >>> writeArr <<< ===========================*/
int mdWrite(int fdes, char* when, char* errMsg, int mode,
	     int size, void *punta_arr)
{
  /* simile a mdRead solo che scrive un'buffer di lunghezza 
     size ( ved. mdRead() ) */
  //printf("mdWrite fdes: %d\n", fdes);
  /* controlla il buon esito della scrittura, se fallisce esce, poiche' 
     la simulazione sarebbe compromessa */
  if (write(fdes, punta_arr, size) == -1)
    {
      mdMsg(ALL, errno, when, "ERROR", "write",
	    errMsg,
	    "Try to restart with: mdsimul -continue",
	    NULL);
     if (mode == EXIT) 
       {
	 exit(-1);
       }
     return -1; /* -1 = ERROR */
    }
return 0; /* 0 = OK */
}

/* ======================== >>> mdClose <<< ================================ */
int mdClose(int fd, char* when, char* errMsg, int mode)
{
  int ret;
  
  if ( (ret = close(fd)) == -1)
    {
      mdMsg(ALL, errno, when, "ERROR", "close",
	    "Can't close file opened for writing",
	    errMsg,
	    NULL);
      if ( mode == EXIT ) 
	{
	  exit(-1);
	}
    }
  return ret; /* 0 = OK, -1 = ERROR (the same as close system call) */
}
#ifdef MD_ALLOC_POLY
/* ======================== >>> loadSegs  <<< ===============================*/
int readSegsPoly(int fdes, char* when, char *errMsg, int mode,
	     int size, void** pointer, ...)
{
  /* Put in each segments pointed by pointers (*pointer) passed as arguments 
     'size' bytes read from the file, whose descriptor is fdes.
     The list of pointer must end with a NULL and you must supply at least one
     pointer.
     - 'errMsg' is an extra message to print on error and 'when' could be 
     NULL or a string, if NULL => print current step , if a string print it. 
     This procedure is very usefule to load coordinates arrays.*/
  void* sptr;
  int i, num;
  unsigned char rerr=0;
  va_list ap;
  va_start(ap, pointer);
  
  if (pointer == NULL) return -rerr; 
  /* If the first pointer read is NULl then read nothing */
  num =  va_arg(ap, int);
  /* Load first pointerof the list */
  for ( i = 0; i < num; i++)
    {
      rerr |= -mdRead(fdes, when, errMsg, mode, size, ((COORD_TYPE**)pointer)[i]);
    }
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      /* if mode == EXIT the if an error occurs, it exits */
      num =  va_arg(ap, int);
      for ( i = 0; i < num; i++)
    	{
	  rerr |= -mdRead(fdes, when, errMsg, mode, size, ((COORD_TYPE**)sptr)[i]);
	  /* mdRead = -1 = ERROR or
	     "    =  0 = OK then rerr = 1 => at least one error occurred */
	}
    }
  va_end(ap);
  return -rerr; /* -1 = ERROR , 0 = OK */
}

#endif
/* ======================== >>> loadSegs  <<< ===============================*/
int readSegs(int fdes, char* when, char *errMsg, int mode,
	     int size, void* pointer, ...)
{
  /* Put in each segments pointed by pointers (*pointer) passed as arguments 
     'size' bytes read from the file, whose descriptor is fdes.
     The list of pointer must end with a NULL and you must supply at least one
     pointer.
     - 'errMsg' is an extra message to print on error and 'when' could be 
     NULL or a string, if NULL => print current step , if a string print it. 
     This procedure is very usefule to load coordinates arrays.*/
  void* sptr;
  unsigned char rerr=0;
  va_list ap;
  va_start(ap, pointer);
  
  if (pointer == NULL) return -rerr; 
  /* If the first pointer read is NULl then read nothing */

  /* Load first pointerof the list */
  rerr |= -mdRead(fdes, when, errMsg, mode, size, pointer);
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      /* if mode == EXIT the if an error occurs, it exits */
      rerr |= -mdRead(fdes, when, errMsg, mode, size, sptr);
      /* mdRead = -1 = ERROR or
	   "    =  0 = OK then rerr = 1 => at least one error occurred */
    }
  va_end(ap);
  return -rerr; /* -1 = ERROR , 0 = OK */
}
#ifdef MD_ALLOC_POLY
/* ======================== >>> saveSegs  <<< ===============================*/
void writeSegsPoly(int fdes, char* when, char* errMsg, int mode, int size,
	       void* pointer, ...)
{
  /* Put each segments 'size' bytes long and pointed by pointers 
     (*pointer) passed as arguments into the file, whose descriptor is fdes.
     The list of pointer must end with a NULL and you must supply at least one
     pointer.
     This procedure is very usefule to save coordinates arrays.
     For explanation of 'errMsg' and 'when' see 'mdWrite' or 'mdRead' */
 
  /* DEBUG: int kk=0; */
  void* sptr;
  va_list ap;
  int i, num;
  va_start(ap, pointer);
  if (pointer == NULL) return; 
  /* writes nothing if the first pointer read is NULL */

  //printf("writing %d\n", fdes);
  num =  va_arg(ap, int);
  /* Load first pointerof the list */
  for ( i = 0; i < num; i++)
    {
      mdWrite(fdes, when, errMsg, mode, size, ((double**)pointer)[i]); 
    }
  /* at least one pointer should be present */
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      /* if mode == EXIT if an error occurs then exit */
      num =  va_arg(ap, int);
      for ( i = 0; i < num; i++)
	{
  	  mdWrite(fdes, when, errMsg, mode, size, ((double**)sptr)[i]);
	}
    }
  /* DEBUG:
     printf("N. %d scritture, scritti %d bytes \n",kk+1,kk*4000+4000);
     printf("Size of params %d\n",sizeof(struct params)); */
  va_end(ap);
}


#endif

/* ======================== >>> saveSegs  <<< ===============================*/
void writeSegs(int fdes, char* when, char* errMsg, int mode, int size,
	       void* pointer, ...)
{
  /* Put each segments 'size' bytes long and pointed by pointers 
     (*pointer) passed as arguments into the file, whose descriptor is fdes.
     The list of pointer must end with a NULL and you must supply at least one
     pointer.
     This procedure is very usefule to save coordinates arrays.
     For explanation of 'errMsg' and 'when' see 'mdWrite' or 'mdRead' */
 
  /* DEBUG: int kk=0; */
  void* sptr;
  va_list ap;
  va_start(ap, pointer);
  if (pointer == NULL) return; 
  /* writes nothing if the first pointer read is NULL */

  //printf("writing %d\n", fdes);

  mdWrite(fdes, when, errMsg, mode, size, pointer); 
  /* at least one pointer should be present */
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      /* if mode == EXIT if an error occurs then exit */
      mdWrite(fdes, when, errMsg, mode, size, sptr);
    }
  /* DEBUG:
     printf("N. %d scritture, scritti %d bytes \n",kk+1,kk*4000+4000);
     printf("Size of params %d\n",sizeof(struct params)); */
  va_end(ap);
}


/* ============================ >>> _mdPrintfNr <<< =======================*/
void _mdPrintf(int mode, char *text)
{
  if (  ( (mode == LOG) || (mode == ALL) ) && 
	(output != stdout)  ) 
    {
      fprintf(output, text);
      fflush(output);
    }
  
  if ( (mode == STD) ||
       (mode == ALL) )
    {
#if defined(MPI)
      if (my_rank == MPI_OUT_PROCESS)
	{
	  printf(text);
	  fflush(stdout);
	}

#else
        printf(text);
	fflush(stdout);
#endif	
    }
}

/* =========================== >>> _printing <<< ===================*/
void _printing(void)
{
  /* Routine used by the printing process, obtained by forking the process
     of rank MPI_OUT_PROCESS, which is tipically 0 */
#if defined(MPI) && defined(PRINTER_PROC)
  char txt[MSG_LEN];
  int end = 0;
  MPI_Status status;
  int endedProcs = 0; 
  
  while (!end) 
    {
      /* Print all incoming messages */
      MPI_Recv(txt, MSG_LEN, MPI_CHAR, MPI_ANY_SOURCE, 
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      
      if (status.MPI_TAG == MD_MPI_PRINT)
	{
	  printf(txt); 
	  fflush(stdout);
	}
      else if (status.MPI_TAG == MD_MPI_END)
	{
	  ++endedProcs;
	  if (endedProcs == numOfProcs)
	    end = 1;
	}
 
    }
#endif	 
}

/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintf(int mode, char *text, ...)
{
/* print text on:
   - stdout if mode = STD 
   - output if mode = LOG
   - output and stdout if mode = ALL 
   without putting a /n after text
   NOTE: STD, LOG and ALL are defined in mdsimul.h */ 
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;
  
  va_start(ap, text);

  /* Put "[R<my_rank>]" at the begin of the string */
#ifdef MPI
  sprintf(tstri, "[R%d]%s", my_rank, text);
#else 
  strcpy(tstri, text);
#endif

  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      //printf("qui!\n");

      strcat(tstri, sptr);
    }
#if defined(MPI) && defined(PRINTER_PROC)
  if ( (mode == ALL || mode == STD) && (my_rank != MPI_OUT_PROCESS) )
    {
       MPI_Send(tstri, strlen(tstri)+1, MPI_CHAR, 
      	       printerRank, 
             MD_MPI_PRINT, MPI_COMM_WORLD);
    } 
#endif 
  _mdPrintf(mode, tstri);
  va_end(ap);    
    
}


/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintfR0(int mode, char *text, ...)
{
/* print text on:
   - stdout if mode = STD 
   - output if mode = LOG
   - output and stdout if mode = ALL 
   without putting a /n after text
   NOTE: STD, LOG and ALL are defined in mdsimul.h */ 
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;

#if defined(MPI)
  if (my_rank != MPI_OUT_PROCESS)
    return;
#endif
  va_start(ap, text);

  /* Put "[R<my_rank>]" at the begin of the string */
#ifdef MPI
  sprintf(tstri, "[R%d]%s", my_rank, text);
#else
  strcpy(tstri, text);
#endif

  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      //printf("qui!\n");

      strcat(tstri, sptr);
    }
#if defined(MPI) && defined(PRINTER_PROC)
  if ( (mode == ALL || mode == STD) && (my_rank != MPI_OUT_PROCESS) )
    {
      MPI_Send(tstri, MSG_LEN, MPI_CHAR, 
	       printerRank, 
	       MD_MPI_PRINT, MPI_COMM_WORLD);
    }
#endif 
  _mdPrintf(mode, tstri);
  va_end(ap);    

}


/* ============================= >>> mdPrintf <<< ===========================*/
void mdPrintfWr(int mode, char* text)
{
  /* print text with a return at the end */

  char tstri[MSG_LEN];
  strcpy(tstri, text); 
  strcat(tstri, "\n");
  mdPrintf(mode, tstri, NULL);
}

/* ============================ >>> mdPrintfSpc <<< =========================*/
void mdPrintfSpcWr(int mode, char* text)
{
  /* put spaces before text to make it more readable */
  char tstri[MSG_LEN];
  strcpy(tstri, "  ");
  strcat(tstri, text);
  mdPrintfWr(mode, tstri);
}

/* ========================== >>> mdPrintfSpcNr <<< =========================*/
void mdPrintfSpc(int mode, char* text)
{
  /* like mdprintfSpc but without the return at the end */
  char tstri[MSG_LEN];
  strcpy(tstri, "  ");
  strcat(tstri, text);
  mdPrintf(mode, tstri, NULL);
}

/* ============================= >>> mdMsg <<< ==============================*/
void mdMsg(int mode, int errnum, char *when, char *errType, char* sysCall,
	   char* text, ...)
{
  /* procedure to write e message, where:
     - mode could be ALL, STD, or LOG (see mdPrintf)
     - errnum must be errno if system call fails, otherwise
     errnumm must be NOSYS ( defined in mdsimul.h )
     - when could be NULL => print STEP <number>
     or a string like "init" or "restore" to specify when the message 
     occurr
     - errType is the error type, that is: "CRITICAL ERROR" or "ERROR"
     - sysCall is the sysCall that fails (used only if errnum!= NOSYS)
     - text is a string containing a text the message, put on the second line. 
     the other pars (of variable number) are other strings.
     Every additional string is put on a separate line 
     EXAMPLE:
     mdMsg(ALL,NOSYS,NULL,"CRITICAL ERROR",NULL,"Unable to...",NULL)
     print if current step is 10:
     STEP 10 [CRITICAL ERROR]: 
     Unable to... */

  char firstLine[MSG_LEN]; /* first line of the message to print */
  char tmps[MSG_LEN];      /* fictitiuos string */
  char *sptr;
  va_list ap;
 
  va_start(ap, text);

   if (when != NULL) 
    {
      strcpy(firstLine, when);
    }
  else 
    {
#ifdef MDLLINT
      sprintf(firstLine, "STEP %lld", (long long int)Oparams.curStep);
#else
      sprintf(firstLine, "STEP %d", (int)Oparams.curStep);
#endif
   }
  sprintf(tmps, " [%s]: ", errType);
  strcat(firstLine, tmps);
  if ( (errnum !=  NOSYS) && (sysCall != NULL) )
    {
      strcat(firstLine, sysCall);
      strcat(firstLine, " - ");
      
      strcat(firstLine, strerror(errnum));/*sys_errlist[errnum]);*/
      /* this array contains the system error
	 and errno indexed the last occurred */
    } 
  mdPrintfWr(mode, firstLine); /* write the builded line */
  /* the list of messages must terminate with NULL */ 
  if (text != NULL)
    {
      mdPrintfSpcWr(mode, text);
    }
  else {
    return;
  }
  while ( ( sptr = va_arg(ap, char *)) != NULL)
    {
      mdPrintfSpcWr(mode, sptr); /* write every additional string on a line
			      apart */
    }
  va_end(ap);
}

/* ============================= >>> MPIgetchar <<< =======================*/
int MPIgetchar(void)
{
#ifdef MPI
  int car;
  int dr;
  if (my_rank == MPI_OUT_PROCESS)
    {
      car =  getchar();
      for (dr = 0; dr < numOfProcs; dr++)
	if (dr != MPI_OUT_PROCESS)
	  {
	    MPI_Send(&car, 1, MPI_INT, dr,
		     MPI_GC, MPI_COMM_WORLD);
	  }
      return car;
    }
  else
   {
     MPI_Recv(&car, 1, MPI_INT, MPI_OUT_PROCESS, MPI_GC, 
	      MPI_COMM_WORLD, &status);
     return car;
   }
  
  //MPI_Barrier(MPI_COMM_WORLD); /* Synchronize all processes */
#else
  return getchar();
#endif
}
/* ============================== >>> openLog <<< ===========================*/
void openLog(char* mode)
{
  char realMode[5];
  struct stat OstatLog; 

  if (stat(absTmpHD(MDS_LOG),&OstatLog)==-1)
    {
      printf("not existing log file, skipping...\n");
      strcpy(realMode, mode);
    }
  else
    {
      /* If the log file is too long and it is request to 
	 append datas truncate it */
      if  ( (OstatLog.st_size > MAXSIZE) &&
	    !strcmp(mode,"a") )
	{
	  printf("%s log file too long, truncated to zero length.\n", 
		 absTmpHD(MDS_LOG));
	  strcpy(realMode,"w"); /* truncate */
	}  
      else 
	{
	  /* in the usual case realMode = mode */
	  strcpy(realMode, mode);
	}
    }
  /* Open log file, mode can be "w","a",etc. ( see fopen manual page ) */
  if (  ( output = fopenMPI(absTmpHD(MDS_LOG), realMode) ) 
	== NULL  )                          /* TRUNCATE LOG FILE IF EXISTS */
    {
      output = stdout; /* in this case stderr unmodified */
    }
}

/* ============================= >>> opendirMPI <<< ======================= */
int fileExists(struct dirent* d, const char* fn)
{
#if defined(MPI) 
  char fileName[MAX_LENGTH];
  char tmpStr[MAX_LENGTH];
  strcpy(tmpStr, fn);
  sprintf(fileName, "%s_R%d", tmpStr, my_rank);
  return !strcmp(d->d_name, fileName);
#else
  return !strcmp(d->d_name, fn);
#endif 
}

/*================================ >>> existDir <<< ========================*/ 
void existDir(char* dirName)
{
  /* DESCRIPTION:
     check if exist directory dirName */
  DIR* rfdir;
  rfdir = opendir(MD_HD_MIS);
  if (rfdir == NULL)
    {
      sprintf(msgStrA, "Not existing directory: %s", MD_HD_MIS);
      mdMsg(ALL, errno, "Init", "CRITICAL ERROR", "opendir",
	    msgStrA,
	    "create it with mkdir",
	    NULL);
      exit(-1);
    }
  closedir(rfdir);
}

void error_on_writing(FILE* f, char* filename, char *where, char *command)
{
  mdMsg(ALL, errno, where, "CRITICAL ERROR", command, NULL);
  fclose(f);
  remove(filename);
  exit(-1);
}

