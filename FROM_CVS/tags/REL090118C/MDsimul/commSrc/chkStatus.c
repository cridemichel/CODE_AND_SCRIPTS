#include<mdsimul.h>
struct progStatus OprogStatus;

FILE *chkout, *output;

void chkEnd(void)
{
  fflush(output);
  fclose(output);
}

void main(int argc, char** argv) 
{
  struct simStat OsimStat0, OsimStat1, OsimStat;
  struct chkStr OchkStr;
  struct stat Ostatc; 
  int fd0,fdr0,fdc0,fd1,fdr1,fdc1,cf,sf;
  int err0,err1;                    /* errnos of status files open */  
  char stri[NAME_LENGTH];    
  //char argBash[NAME_LENGTH];

  /* absolute file names */
  char chkLogFile[NAME_LENGTH];
  char chkFile[NAME_LENGTH];
  char SF0[NAME_LENGTH];
  char SF1[NAME_LENGTH];
  
  atexit(chkEnd);
  /* build absolute names for files used by this program */
  
  strcpy(chkLogFile, MD_HD_TMP);
  strcat(chkLogFile, CHK_LOG); /* chkLogFile is an absolute log file name */
  strcpy(chkFile, MD_HD_TMP);
  strcat(chkFile, CF);         /* chkFile is an absolute check file name */
  
  stat(chkLogFile,&Ostatc); 
 
  /* SF0 and SF1 cotains absolute name for status files */
  strcpy(SF0, MD_HD_TMP);
  strcat(SF0, STATUS_FILE_NAME);
  sprintf(stri, "%d", 0);
  strcat(SF0, stri);
 
  strcpy(SF1, MD_HD_TMP);
  strcat(SF1, STATUS_FILE_NAME);
  sprintf(stri, "%d", 1);
  strcat(SF1, stri);
  
  /* check if the log file is too long, if it is so truncate it */
  if (Ostatc.st_size > MAXSIZE) 
    {
      strcpy(stri,"w"); /* truncate */
    }  
  else 
    {
      strcpy(stri,"a"); /* append */
    }
  
  //printf("log file: %s\n",chkLogFile);
  if ( (output = fopenMPI(chkLogFile, stri)) == NULL ) 
    {
      output = stdout;
    }
  
  /* Check if already exists CHK_FILE, if Y then load chkStr structure
     fromt it and the status from the status file specified in the 
     that structure, if N read status file and store chkStr in the 
     CHK_FILE */ 
  if ( ( cf = openMPI(chkFile, O_RDONLY)) != -1 ) 
    {
      
      if ( read(cf, &OchkStr, sizeof(struct chkStr)) == -1) 
	{
	  /* writes on output stream: it is a log file named chkFile 
	     (see below) */
	  fprintf(output,">>> Error reading CHK_FILE, retrying later...\n");
	  exit(-1);
	}
      if ( close(cf) == -1 ) 
	{
	  fprintf(output,">>> Error closing CHK_FILE, retrying later...\n");
	  exit(-1);
	}
      //printf("Current step read CHK_FILE: %d\n",OchkStr.curStep);
      if ( OchkStr.whichSf == 0 ) 
	{
	  strcpy(stri,SF0);
	}      
      else if ( OchkStr.whichSf == 1 ) 
	{
	  strcpy(stri,SF1);
	}      
      else 
	{
	  fprintf(output,">>> UNUSUAL ERROR, which_sf attribute of chk_str struct not 0 nor 1\n");
	  fprintf(output,"    I assume which_sf = 1\n");
             
	  strcpy(stri, SF1);
	}
      /* If the status file specified in the check file is not "openable" or 
	 unreadable exit, it will retry later */
      if ( ( sf = openMPI(stri,O_RDONLY)) == -1 ) 
	{
	  /* if somthing goes wrong opening, reading or closing, it deletes
	     the CHK_FILE, so that is made another trial ab initio */
	  unlink(chkFile);
	  fprintf(output,">>> Error opening status file, retrying later...\n");
	  exit(-1); 
	}
	
      if ( read(sf,&OsimStat,sizeof(OsimStat)) == -1 ) 
	{
	  unlink(chkFile); /* see above (open) */
	  fprintf(output,">>> Error reading status file, retrying later...\n");
	  exit(-1);
	}
	
      if ( close(sf) == -1 ) 
	{
	  unlink(chkFile); /* see above (open) */
	  fprintf(output,">>> Error closing status file, retrying later...\n");
	  exit(-1);
	}
      /* The simulation is not going ? if N execute 'mdsimul -c', i.e.
	 restart simulation */
      if ( OsimStat.curStep == OchkStr.curStep ) 
	{
	  unlink(chkFile);
#ifdef XDISP
	  if ( execlp(XTERM,
		      "nxterm","-display", ":0.0",
		      "-e", MDSIMUL, "-c", NULL)== -1 )
	    fprintf(output,">>> CRITICAL ERROR: I can't execute mdsimul!\n"); 
	    //printf("versione X\n");
#else
	    if ( execlp(MDSIMUL,
			"mdsimul","-c",NULL)== -1 )
	      fprintf(output,">>> CRITICAL ERROR: I can't execute mdsimul!\n");
#endif
	}
      else 
	{
	  /* If the file exists but the simulation seems to be going then 
	     unlink the 'chkFile' */
	  unlink(chkFile);
	  exit(-1);
	}
    } 
  /* If the file exists but it is not possible to open it, then exit */
  else if ( errno != ENOENT) 
    {
      fprintf(output,">>> CHK_FILE exists but I can't read, retrying later...\n");
      unlink(chkFile);
      exit(-1);
    }
  /* the file doesn't exist create CHK_FILE and exit */
  else 
    {
      /* if the file CHK_FILE doesn't exist then create it and strore in 
	 the curStep attribute */
      /* open, read and close both status file */
      fd0 = openMPI(SF0,O_RDONLY); 
      err0 = errno;                       /* err0 contains error id */
      fdr0 = read(fd0,&OsimStat0,sizeof(struct simStat));
      fdc0 = close(fd0);
      fd1 = openMPI(SF1,O_RDONLY); 
      err1 = errno;                       /* err1 contains error id */
      fdr1 = read(fd1,&OsimStat1,sizeof(struct simStat));
      fdc1 = close(fd1);
      if ( (fd0 == -1) && (err0 == ENOENT) &&  /* ENOENT = no such file */
	   (fd1 == -1) && (err1 == ENOENT) )
	{
	  exit(-1); /* exit without complaining if status files don't exist */
	}
      
      /* If both status files are corrupted ( not opened or not readable ) 
	 then exit. */
      if (  ( ( fd0 == -1) || (fdr0 == -1) ) &&
	    ( (fd1 == -1)  || (fdr1 == -1) )  )
	{
	  fprintf(output,">>> Error: both status file corrupted\n");
	  fprintf(output,"    or no interrupted simulation.\n");
	  exit(-1);
	}
      /* take the good structure */
      if ( (fd0 == -1) || (fdr0 == -1) ) 
	{
	  OchkStr.curStep = OsimStat1.curStep;
	  OchkStr.whichSf = 1;
	}
      else 
	{ 
	  OchkStr.curStep = OsimStat0.curStep;
	  OchkStr.whichSf = 0;
	}
      if ( (cf = creatMPI(chkFile,0666)) == -1 ) 
	{
	  fprintf(output,">>> Error opening CHK_FILE, retrying later...\n");
	  exit(-1);
	}
      if (write(cf,&OchkStr,sizeof(struct chkStr)) == -1) 
	{
	  unlink(chkFile);
	  fprintf(output,">>> Error writing CHK_FILE, retrying later...\n");
	  exit(-1);
	}
      if  (close(cf) == -1) 
	{
	  unlink(chkFile);
	  fprintf(output,">>> Error closing CHK_FILE, retrying later...");
	  exit(-1);
	}
    }
}

