#undef MAIN
#include<mdsimul.h>

/* string used to send  messages */
extern char msgStrA[MSG_LEN],msgStrB[MSG_LEN],msgStrC[MSG_LEN]; 

extern struct simStat OsimStat;

extern int my_rank;
extern unsigned char BAK, STA;
extern unsigned char BAKT;            /* global switch for restore on tape*/
extern int whichCorrupted(char* absFileName, int (*readFunc)(int), unsigned char* Pnewer);

/* ========================== >>> chooseMeasure <<< =========================*/
void chooseMeasure(char* absFile, int (*readFunc)(int))  /* file -> absolute */
{
  /* DESCRIPTION:
     Choose the good measure 'file':
     if both measure 'file0' and 'file1' were good, then don't touch 
     files, otherwise copy the good file onto the corrupted one, using
     the copy() function 
     readFunc is a pointer to the function to use to check datas, reading it.
     NOTE: file must an absolute measure file name */

  int wc;             /* return value of whichCorrupted function */
  char src[NAME_LENGTH], dest[NAME_LENGTH];
  unsigned char newer;/* flag, if 1 means that whichGood has choosen then
			 newer file */
  if (  ( wc = whichCorrupted(absFile, readFunc, &newer) ) == -1  )
	{  
	  sprintf(msgStrA, "The measure in the file %s is lost!",
		  absFile);
	  mdMsg(ALL, NOSYS, "Restore", "ALERT", NULL,
		msgStrA,
		NULL);
	}
  /* DEBUGGING */ 
  //printf("newer: %d\n",newer);
  //printf("wc = %d sw(wc) = %d\n", wc, sw(wc));
  
  /* if 'whichCorrupted' has choose the newer file don't copy one file 
     onto the other one */
  if (newer == 1) /* If the two files are identical except for date 
		     don't touch them */
    {
      return;
    }
  
  /* If whichCorrupted has choose a corrupted file delete that one 
     (wc is the corrupted measure file), then build a command of this 
     form:
     /bin/cp <good measure file> <corrupted one>, 
     that is copy the good file onto the good one*/
#ifdef MPI
   sprintf(src, "%s_R%d", appSw(absFile, sw(wc)), my_rank);
   sprintf(dest, "%s_R%d", appSw(absFile, wc), my_rank);
#else
   strcpy( src, appSw(absFile, sw(wc)) );
   /* sw(wc) = good measure file */
   strcpy( dest, appSw(absFile, wc) );
#endif
   /* wc is the corrupted file */
  if ( copy(src, dest) == -1 )
    {
      sprintf(msgStrA, "Unable to copy %s in %s\n", src, dest);
      mdMsg(ALL, errno, "init", "WARNING", "open",
	    msgStrA,
	    NULL);
    }
}

/* ========================== >>> chooseREstore <<< =========================*/
int chooseRestore(char* absFile, unsigned char* PbakSwitch, 
		   int readFunc(int)) 
{
  /* DESCRIPTION: 
     choose the next restore file to write to ( the switch 
     pointed by PbakSwitch contains the next restore file ).
     - readFunc is a pointer to the function used to read the file 
     RETURN VALUE: -1 if both files corrupted 
                    0 ptherwise */ 

  int wc;
  unsigned char newer;/* flag, if 1 means that whichGoos has choosen then
			 newer file */
  if ( (wc = whichCorrupted(absFile, readFunc, &newer)) == -1)
    /* file->absolute */
    {  
      *PbakSwitch = 1; /* choose arbitrarly */
      return -1; /* -1 = both files corrupted */
    }
  *PbakSwitch = (unsigned char) wc; 
  /* set the switch pointed by PbakSwitch, that is the next restore file to 
     save in, to the corrupted file.
     Actually the pointed switch are BAK and BAKT */ 
  /* DEBUGGING !!!!!!!!!!*/
 // printf("next restore file to write on : %d\n", BAK);
  return 0; /* 0 = ok et least one file good */
}

/* ========================== >>> chooseStatus <<< ========================*/
void chooseStatus(char* absFile, int (*readFunc)(int))
{
  /* DESCRIPTION:
     Choose the good status file, and set STA so that the next write is on 
     status file corrupted
     absFile is the absolute file name of the status file, with path but 
     without a switch appended (that 0 or 1 at the end of the string)
     >>> readFunc pointer ACTUALLT NOT USED !!!!!!!!!!!*/ 
  int sf0, sf1, sfr0, 
    sfr1, sfc0, sfc1; /* (see later) */
  int err0, err1;     /* used to store errno file after opening of 
			status file */

  mdMsg(ALL, NOSYS, "restore", "NOTICE", NULL,
	"Choosing STATUS FILE...",
	NULL);

  sf0 = openMPI( appSw(absFile, 0),
	      O_RDONLY );
  err0 = errno;
  

  
  sf1 = openMPI( appSw(absFile, 1),
	      O_RDONLY );
  err1 = errno;
  
  sfr0 = read(sf0, &OsimStat, sizeof(struct simStat));
  sfr1 = read(sf1, &OsimStat, sizeof(struct simStat));
  sfc0 = close(sf0);
  sfc1 = close(sf1);
  
  if (  (sf0 == -1) && (err0 == ENOENT) && 
	(sf1 == -1) && (err1 == ENOENT)  )
    { 
       mdMsg(ALL, NOSYS, NULL, "WARNING", NULL,
	     "No status file found, probably user has interrupted simulation",
	     NULL); 
       STA = 0; /* choose a value for STA, it doesn't matter which */
    }
  
  /* This condition says: if I'can't open OR read OR close both status 
     file, they are both corrupted */
  
  else if (  ( (sf0 == -1) || (sfr0 == -1) || (sfc0 == -1) ) &&
	     ( (sf1 == -1) || (sfr1 == -1) || (sfc1 == -1) )  )
    {
      mdMsg(ALL, NOSYS, NULL, "WARNING", NULL,
	    "Both status file corrupted, actually if a system crash occured",
	    "I couldn't restart.",
	    NULL);
      /* if both files are corrupted, but one file cannot be closed only,
         then choose the other file as next status file */
      if ( (sf0 != -1) && (sfr0 != -1) && (sfc0 == -1) )
	{
	  STA = 1;
	}      
      else
	{
	  STA = 0;
	}
    }
  /* ----------- CHOOSE THE CORRUPTED ONE ( see TECH_INFO) */
  else if ( (sf0 == -1) || (sfr0 == -1) || (sfc0 == -1) )  /* sf0 is corrupted?
							      if y use that */
    {
      STA = 0;
    }
  else
    {
      STA = 1;
    }
 
  /* DEBUGGING !!!!!!!!!!*/
  //printf("Status file to write on: %d\n", STA);

}

/* =========================== >>> TryOlderFile <<< =========================*/
int  TryOlderFile(int bf1, int (*readFunc)(int), 
		  int which)
{
  /* DESCRIPTION:
     bf1 is the descriptor of the older file,
     readFunc is apointer to a function able to read from bf1 all the datas,
     and which is the actual choosen file, that is the file choose before
     calling this function.
     RETURNN VALUE: -1 if the older file is also corrupted,
                     the corrupted file otherwise, that is the newer file */
  
  /* Most recent file corrupted, try the other ...  */
  if ((*readFunc)(bf1) == -1)
    {				/* ERROR READING ? */
      mdMsg(ALL,NOSYS, "restore", "CRITICAL ERROR",NULL,
	    "Unable to read both files",
	    NULL);
     return -1;
    }
  /* oldest file successfully read, try to close it */
  else if (close(bf1) == -1)
    {				/* ERROR CLOSING ? */
      mdMsg(ALL, errno, NULL, "CRITICAL ERROR", "close",
	    "Error closing oldest restore file",
	    NULL);
      return -1;
    }
  return ~which & 1;	/* 1 => 0 or 0 => 1 , i.e. next file 
			   is bf0 (newer) , that is the corrupted one, in fact
			   previuosly we set 'which' to the older file */
}

/* ========================= >>> whichCorrupted <<< =========================*/
int whichCorrupted(char* absFileName, int (*readFunc)(int), 
		   unsigned char* Pnewer)
{
  /* DESCRIPTION:
     If a file named pippo is doubleSave() or doubleBuffer(), on the disk 
     you have two files named 'pippo0' and 'pippo1'.
     This procedure returns 1 if the good file is 'pippo0' or 0 if 
     the good file is 'pippo1', that is the corrupted one.
     - readFunc is a pointer to a function that is able to read the files 
     'pippo0' and 'pippo1', and that returns 0 if the read was made 
     successfully or -1 otherwise.
     It is very important to note that the datas are read once from the file,
     so by readFunc you can also use whichFileGood to store them into memory.
     This is done for example when choosing restore file.
     NOTE: the first arg must an absolut file name, that is a filename with
     path. */ 
  
  struct stat Ostat0, Ostat1;  /* structures for file infos (see 'man stat') */
  /* Load restore file using an algorithm explained in the TECH_INFO
     file */
  int bf, bf0, bf1 ;   /* restore file descriptors (see later) */
  int sr1,sr2;         /* return value of the stat system call */   
  
  int which;           /* integer containing the good file */
  char fileA[NAME_LENGTH], fileB[NAME_LENGTH];
#ifdef MPI
  char fA[NAME_LENGTH], fB[NAME_LENGTH];
#endif
  *Pnewer = 0; /* initially the flag is 0 that is : 
		  'one file at least corrupted' */

  /* build absolute names of the two measure files used to realize 
     double buffer */
  strcpy( fileA,
	  appSw(absFileName, 0) ); /* absFileName must be absolute */
  
  strcpy( fileB,
	  appSw(absFileName, 1) );
  
#ifdef MPI
  sprintf(fA, "%s_R%d", fileA, my_rank);
  sprintf(fB, "%s_R%d", fileB, my_rank);
  sr1 = stat(fA, &Ostat0);
  sr2 = stat(fB, &Ostat1);
#else
  sr1 = stat(fileA, &Ostat0);
  sr2 = stat(fileB, &Ostat1);
#endif

  bf0 = openMPI(fileA, O_RDONLY);
  bf1 = openMPI(fileB, O_RDONLY);
  
  if ( (bf0 == -1) && 
       (bf1 == -1) )
    {
      mdMsg(ALL, errno, NULL, "CRITICAL ERROR", "open",
	    "Unable to read both files",
	    NULL);
      return -1;
    }
  /* Chooses the file longer, this criterion works because if a 
     problem (system crash) occurs when writing the last file
     his length is not right,i.e. shorter than the other file
     (see TECH_INFO for details).
     NOTE: if the two files have same length see below. */

  if ( (sr1 == -1) && 
       (sr2 == -1) )
    {
      perror("boh");
      mdMsg(ALL, NOSYS, "Choose", "CRITICAL ERROR", NULL, 
	    "Both stat calls returns -1 (error).",
	    "and both files was opened regularly, I can't decide which",
	    "file is good, try manually...",
	    NULL);
      return -1;
    }
  /* IS ONE FILE NOT OPENED ? */
  if ( (bf0 == -1) || 
       (bf1 == -1) )
    {
      /* if bf0 not opened use bf1 */
      if (bf0 == -1)
	{
	  bf = bf1;
	  which = 0;	 /* next file to write on is the
			    corrupted one, i.e. 'fileA' */
	}
      /* if bf1 not opened use bf0 */
      else
	{
	  bf = bf0;
	  which = 1;	 /* next restore to write on is the corrupted one,
			    i.e. 'fileB' */
	}
      if ((*readFunc)(bf) == -1)
	{
	  mdMsg(ALL, NOSYS, "Restore", "CRITICAL ERROR", NULL,
		"There is only one file and I'm unable to read it",
		NULL);
	  return -1;
	}
      if (close(bf) == -1)
	{
	  mdMsg(ALL, errno, "Restore", "CRITICAL ERROR", "close",
		"There is only one file and I'm unable to close it",
		NULL);
	  return -1;
	}
    }
  /* both files regularly opened => IS ONE FILE LONGER ? */
  else if (Ostat0.st_size != Ostat1.st_size)
    {
      /* bf1 longer than bf0 ? => interchange bf0 and bf1,
         so bf0 is the longest file in anycase */
      if (Ostat0.st_size < Ostat1.st_size)
	{
	  bf = bf0;
	  bf0 = bf1;
	  bf1 = bf;
	  which = 0;	 /* next file to write on is the corrupted one,
			    i.e. 'fileA' */
	}
      else
	{
	  which = 1;	 /* next file to write on is the corrupted one,
			    i.e. 'fileB' in this case */
	}
      /* bf0 is always the longest file (see previous if statement) */
      if ((*readFunc)(bf0) == -1)
	{
	  mdMsg(ALL, NOSYS, "Restore", "CRITICAL ERROR", NULL,
		"Longest file not readable",
		"Other file probably corrupted",
		NULL);
	  return -1;
	}
      if (close(bf0) == -1)
	{
	  mdMsg(ALL, errno, "Restore", "CRITICAL ERROR", "close", NULL,
		"I can't close longest file",
		"Other file probably corrupted",
		NULL);
	  return -1;
	}
      /* AGGIUNTO 07/05/2001 */
      //close(bf1);
    }
  /* if the two files are both regularly opened and have same length 
     than try the most recent one, i.e. the file most recently
     modified --- IS ONE FILE MOST RECENT ? */
  else
    {
      /* if bf1 is more recent than bf0 then interchange  bf0 and bf1,
         i.e. new bf1 = old bf0 and  new bf0 = old bf1.
         In this way bf0 is always the most recent file. 
         NOTE1: difftime(t1,t2) = t1 - t2 in seconds 
         NOTE2: which is set to point to the older file <------- !!!!!! */
      if (difftime(Ostat1.st_mtime, Ostat0.st_mtime) > 0)
	{
	  bf = bf0;
	  bf0 = bf1;
	  bf1 = bf;
	  which = 0;       /* next file to write on is 'fileA',
			      i.e. older one */
	}
      else
	{
	  which = 1;       /* next restore file to write on is 'fileB',
			      i.e. older one */
	}
      /* bf0 is the most recent file 
         ------ LOAD MOST RECENT FILE */
      if ((*readFunc)(bf0) == -1)
	{			/* ----- ERROR READING ? */
	  mdMsg(ALL, NOSYS, "Choose", "WARNING", NULL,
		"Most recent file not readable so I try the older one.",
		NULL);
	  /* here which is older restore file */	
	  which = TryOlderFile(bf1, readFunc, which);
	  printf("Older file is ok, most recent one it is not, delete most recent file %s%d,", BAK_FILE_NAME, which);
	  printf("then restart simulation with -c option\n");
	  /* esco poichè avendo chiamato a questo punto readFunc due volte è meglio uscire, in quanto 
	     ha allocato tutto due volte! */
	  exit(-1);
	  /* AGGIUNTO 07/05/2001 */
	  //close(bf0);
	}

      /* Most recent restore file successfully read, try to close it */
      else if (close(bf0) == -1)
	{			/* ERROR CLOSING ? */
	  mdMsg(ALL, errno, "Restore", "WARNING", "close",
		"Error closing most recent measure file",
		NULL);
	  which = TryOlderFile(bf1, readFunc, which);
	  printf("Older file is ok, most recent one it is not, delete most recent file %s%d,", BAK_FILE_NAME, which);
	  printf("then restart simulation with -c option\n");
	  /* esco poichè avendo chiamato a questo punto readFunc due volte è meglio uscire, in quanto 
	     ha allocato tutto due volte! */
	  exit(-1);
	}
      
      /* If reach this point the newer file was regularly opened, read and 
	 closed */
      *Pnewer = 1; /* the return value of this procedure refers to the 
		     newer file */
    }

  if (bf0!=-1)
    close(bf0);
  if (bf1!=-1)
    close(bf1);
  return which; /* 'which' is the file corrupted (0/1) */

}
