#undef MAIN
#include<mdsimul.h>
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   +                   THIS MODULES CONTAINS A LOT OF                        +
   +                     SIMULATION DEPENDENT CODE                           +
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

/* Extern variables must be of the same type of the real definition */

/* =========================== >>> STRUCTURES <<< ===========================*/

extern struct simStat OsimStat;

/* ==========================================================================*/

extern int SEGSIZE;        
/* length in bytes of an array */

/* string used to send  messages */
extern char msgStrA[MSG_LEN],msgStrB[MSG_LEN],msgStrC[MSG_LEN]; 
								   
extern unsigned char BAK, STA, BAKT;

extern int *rxS[NA], *ryS[NA], *rzS[NA];

#ifdef LMC
#define CORTYPE int
#else
#define CORTYPE double
#endif

/* AGGIUNTO 23/04/01 */
extern void writeAllCor(FILE* );
extern void readAllCor(FILE* );
const char sepStr[] = "@@@\n";

/* ============================ >>> AllocCoord <<< ==========================*/
void AllocCoord(int size, CORTYPE** pointer, ...)
{
  /* Allocate contigously memory */
  va_list ap;
  CORTYPE** ptrs[MAXVARS];
  int totBytes, totArr; 
  int i;
  va_start(ap, pointer);
  
  ptrs[0] = pointer;
  totBytes = size;
  for (i=1; (ptrs[i] = va_arg(ap, CORTYPE**)) != 0; ++i)
    { 
      totBytes += size; /* total bytes to allocate */
    }
  totArr = i; /* total number of arrays to allocate */
  
  *ptrs[0] = malloc(totBytes);
  
  // printf("totBytes = %d size: %d\n", totBytes, size);
  for (i = 1; i < totArr; ++i)
    {
      *ptrs[i] =(CORTYPE*)((char*) *ptrs[i-1] + size);  
      //intf("pointer value ptrs[%d]:%d\n", i, (int)*ptrs[i]);
      //intf("*sptr relative %d\n\n",(int) *ptrs[i] - (int)*ptrs[0]);
      //intf("totArr: %d\n", totArr);
    }	
  va_end(ap);
}

/* ==================== >>> AllocCoordFrag <<< =====================*/
void AllocCoordFrag(int size, COORD_TYPE** pointer, ...)
{
  va_list ap;
  CORTYPE** sptr;
  va_start(ap, pointer);
  *pointer = malloc(size);
  while((sptr = va_arg(ap, CORTYPE**)) != NULL)
    {
      *sptr = malloc(size); 
    }
  va_end(ap);
}


/* ======================== >>> Allocint <<< ==============================*/
void AllocInt(int num, int **ptr, ...)
{
  /* DESCRIPTION:
     'mode' could be SHR = shared arrays or NOT_SHR = not shared  
     The other arguments are of the form;
     <num>, <ptr>, <num1>, <ptr1>, ...
     and it allocate an array of <num> integer for the pointer '*<ptr>'
     (note that <ptr> is the address of the pointer not the pointer itself)
     NOTE: max you can allocate arrays for 16Mb
  */
  va_list ap;
  int** ptrs[MAXVARS];
  int arrLen[MAXVARS], totBytes, totArr; 
  int i;
  va_start(ap, ptr);
  
  ptrs[0] = ptr;
  arrLen[0] = num;   /* number of elements for each array */
  totBytes = num * sizeof(int);
  for (i=1; (arrLen[i] = va_arg(ap, int)) != 0; ++i)
    { 
      totBytes += arrLen[i] * sizeof(int); /* total bytes to allocate */
      ptrs[i]= va_arg(ap, int**);  
    }
  totArr = i; /* total number of arrays to allocate */
  
  *ptrs[0] = malloc(totBytes);
  for (i = 1; i < totArr; ++i)
    {
      *ptrs[i] = *ptrs[i-1] + arrLen[i-1];  
      
      //printf("pointer value ptrs[%d]:%d, arrLen: %d\n", i, (int)*ptrs[i], arrLen[i]);
      //printf("*sptr relative %d\n\n",(int) *ptrs[i] - (int)*ptrs[0]);
      //printf("totArr: %d\n", totArr);
    }	
  va_end(ap); 
}

/* ============================ >>> AllocCT <<< =========================== */
void AllocCT(int num, COORD_TYPE **ptr, ...)
{
  /* DESCRIPTION:
     Arguments are of the form;
     <num>, <ptr>, <num1>, <ptr1>, ...
     and it allocates an array of <num> 'COORD_TYPE' numbers for 
     the pointer '*<ptr>'
     (note that <ptr> is the address of the pointer not the pointer itself)
     NOTE: max you can allocate arrays for 16Mb
  */
  va_list ap;
  COORD_TYPE** ptrs[MAXVARS];
  int arrLen[MAXVARS], totBytes, totArr; 
  int i;
  va_start(ap, ptr);
  
  ptrs[0] = ptr;
  arrLen[0] = num;   /* number of elements for each array */
  totBytes = num * sizeof(COORD_TYPE);   
  for (i=1; (arrLen[i] = va_arg(ap, int)) != 0; ++i)
    { 
      totBytes += arrLen[i] * sizeof(COORD_TYPE); /*total bytes to allocate*/
      ptrs[i] = va_arg(ap, COORD_TYPE**);  
    }
  totArr = i; /* total number of arrays to allocate */
 
  *ptrs[0] = malloc(totBytes);
    
  // DEBUGGING 
  //printf("Tot bytes: %d, totArr: %d\n", totBytes, totArr);
  
  for (i = 1; i < totArr; ++i)
    {
      *ptrs[i] = *ptrs[i-1] + arrLen[i-1];  
    }	
  va_end(ap);
}

/* ==================== >>> DeleteCoord(NOT USED!!!) <<< ====================*/
void DeleteCoord(COORD_TYPE **pointer, ...)
{
  /* ======== >>> See AllocCoord() for an explanation. <<<=======
     Remove all shared memory allocated for coordinates 
     WARNING: If it deallocs an uninitialized  pointer  we obtain a
     'segmentation fault' ! */
  COORD_TYPE** sptr;
  va_list ap;
  va_start(ap,pointer);
  free((void *) *pointer);
  while ( (sptr = va_arg(ap, COORD_TYPE**)) != NULL )
    {
      /* sptr contains the address of a COORD_TYPE* variable, passed as 
	 arguments, so if for example sptr = &ax1 then 
	 *sptr = *(&ax1), that is ax1 itself!
	 So the pointer passed to freeshmem() is the one stored in
	 ax1, and this is just what we want.*/ 
      free((void *) *sptr);
    }
  va_end(ap);
}

/* Local Memory allocation function */

/* Allocate a vectore of COORD_TYPE */
COORD_TYPE* AllocVecR(int size)
{
  return (COORD_TYPE*) malloc (size * sizeof(COORD_TYPE));
}

void FreeVecR(COORD_TYPE* v)
{
  free(v);
}

/* Allocate a vectore of int */
int* AllocVecI(int size)
{
  return (int*) malloc (size * sizeof(int));
}

void FreeVecI(int* v)
{
  free(v);
}

/* Allocate memory for a matrix of COORD_TYPE */
COORD_TYPE** AllocMatR(int size1, int size2)
{
  COORD_TYPE** v;
  int k;
  v = (COORD_TYPE**) malloc(size1 * sizeof(COORD_TYPE*));
  v[0] = (COORD_TYPE*) malloc(size1 * size2 * sizeof(COORD_TYPE));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}

void FreeMatR(COORD_TYPE** v)
{
  free(v[0]);
  free(v);
}

/* Allocate memory for a matrix of integers */
int** AllocMatI(int size1, int size2)
{
  int** v;
  int k;
  v = (int**) malloc(size1 * sizeof(COORD_TYPE*));
  v[0] = (int*) malloc(size1 * size2 * sizeof(COORD_TYPE));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}
void FreeMatI(int** v)
{
  free(v[0]);
  free(v);
}

/* =========================== >>> setToZero <<< ========================== */
void setToZero(int parnum, CORTYPE* ptr, ...)
{
  /* DESCRIPTION:
     This procedure set to zero all the coordinates, passed as arguments */
  va_list ap;
  CORTYPE* sptr;
  int i;
  
  va_start(ap, ptr);
  for(i = 0; i < parnum; ++i) ptr[i]=0.0;
  while ( (sptr = va_arg(ap, CORTYPE*)) != NULL)
    {
       for(i = 0; i < parnum; ++i) sptr[i] = 0.0;
    }
  va_end(ap);
}

/* ========================= >>> ReadCoord <<< =============================*/
int readCoord(int cfd)
{
  /* DESCRIPTION:
     Read structure params and coordinates from file whose file descriptor is 
     cfd.
     NOTE: nel salvataggio rispettare l'ordine di caricamento <----!!!!!
     Attualmente non viene effettuato nessun controllo sugli fread <----!!!!
   */
  unsigned char rerr = 0;

  mdRead( cfd, "Init", "I can't read params structure", EXIT,
	  sizeof(struct params), &Oparams);

  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];
 /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  AllocCoord(SEGSIZE, ALLOC_LISTA,
	     NULL);

  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];
  AllocCoord(SEGSIZE, ALLOC_LISTB,
	     NULL);

  /* loads all arrays from the file associated with the fdes descriptor */
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];
  rerr |= -readSegs(cfd, "Init", "Error reading coordinates", CONT,
		    SEGSIZE, SAVE_LISTA,
		    NULL);    /* NULL means: 'no more pointers to load' */
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];
  rerr |= -readSegs(cfd, "Init", "Error reading coordinates", CONT,
		    SEGSIZE, SAVE_LISTB,
		    NULL);

  rerr |= -readSegs(cfd, "Init", "Error reading extra coordinates", CONT,
		    sizeof(COORD_TYPE), EXT_SLST,
		    NULL);    /* NULL means: 'no more pointers to load' */
  return -rerr;
  /* -1 = ERROR that        
     0 = OK */
}

/* =========================== >>> saveCoord <<< ===========================*/
void saveCoord(char* fileName)
{
  /* DESCRIPTION:
     save last coordinates on file fileName */
 
  int cfd; /* descriptor of coordinate file named fileName */ 
 
  cfd = mdCreat(fileName, "End",
		"Final coordinates are in the restore files on HD and Tape", 
		EXIT);  
   
  mdWrite(cfd, NULL, 
	  "Error writing the params struct.", EXIT,
	  sizeof(struct params), &Oparams);

  /* writes all arrays to the disk physically making a sync() */
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];
  writeSegs(cfd, "End", "Error writing final coordinates", EXIT,
	    SEGSIZE, SAVE_LISTA,
	    NULL);
  
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];
  writeSegs(cfd, "End", "Error writing final coordinates", EXIT,
	    SEGSIZE, SAVE_LISTB,
	    NULL);
  /* NOTE:
     Actually the EXT_SLST should be a list of COORD_TYPE variables to save */
  writeSegs(cfd, "End", "Error writing file coordinates", EXIT, 
	    sizeof(COORD_TYPE), EXT_SLST,
	    NULL);

  mdClose(cfd, "End", 
	  "Final coordinates are in restore files on HD and Tape", EXIT);

  sync(); /* <-------------------------------------------SYNC() !!!!!!!!!*/
  
  mdMsg(ALL, NOSYS, "End", "NOTICE", NULL,
	"Final coordinates file successfully saved",
	NULL);
}
  
/* ======================== >>> ReadBak <<< ================================*/ 
int readBak(int bfd)
{	
  int br; /* bytes read by 'read' system call */
  unsigned char rerr = 0;
  /* Read structures params and filenames and coordinates from restore file
     whose file descriptor is 'bfd'. 
     NOTE: nel salvataggio rispettare l'ordine di caricamento <---------!!!!!
     Attualmente non viene effettuato nessun controllo sugli fread <----!!!!!*/
  
  /* Ensure to start reading from begin of file */
  lseek(bfd, 0, SEEK_SET);
  
  /* reads params structure ( see TECH_INFO for details) */
  
  br = mdRead(bfd , "Restore", "Error reading the params struct.",
	      CONT, sizeof(struct params), &Oparams);
  if (br == -1) 
    {
      return 1; /* 1 = ERROR for readBak */
    }

  /* NOW we know the particles number so we can allocate shared memory, to
     get their initial coordinates */
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];
  /* SEGSIZE is the size in bytes of an array of coordinates */

  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */

  AllocCoord(SEGSIZE, ALLOC_LISTA,
	     NULL);
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];

  AllocCoord(SEGSIZE, ALLOC_LISTB,
	     NULL);
  
  /* reads filenames structure ( see TECH_INFO for details) */
  br = mdRead(bfd, "Restore",  "Error reading the program status.", CONT, 
	      sizeof(struct progStatus), &OprogStatus);
  
 if (br == -1)
    { 
      return 1; /* 1 = ERROR */
    }
 
 /* loads all arrays from the file associated with the fdes descriptor */

 SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];
 rerr |= -readSegs(bfd, "Restore", "Error reading restore file", CONT, 
		   SEGSIZE, SAVE_LISTA,
		   NULL);   /* NULL means: 'no more pointers to load' */
 SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];
 rerr |= -readSegs(bfd, "Restore", "Error reading restore file", CONT, 
		   SEGSIZE, SAVE_LISTB,
		   NULL);   /* NULL means: 'no more pointers to load' */

 rerr |= -readSegs(bfd, "Restore", 
		   "Error reading extra coords from restore file", 
		   CONT, sizeof(COORD_TYPE), EXT_SLST,
		   NULL);   /* NULL means: 'no more pointers to load' */
 return rerr; /* 1 = ERROR , 0 = OK */
}
/* ============================== >>> Parsing <<< ===========================*/
void asciiParsing(struct pascii strutt[], 
		  char* stringA, char* stringB)
{
  int i, n, e;
  double* bd;
  int *bi;
  char *bc, *subs;
  char *bs;
  
  for  (i=0; strutt[i].ptr != NULL; ++i) /* parname=NULL menas END */
    {
      if (!strcmp(strutt[i].parName, stringA))
	{
	  if (strchr(strutt[i].type,'f') || strchr(strutt[i].type, 'e')
	      || strchr(strutt[i].type, 'E') || 
	      strchr(strutt[i].type, 'g') || strchr(strutt[i].type, 'G'))
	    {
	      bd = (double *) strutt[i].ptr;
	      subs = strtok(stringB, " "); //lista separata da spazi
	      n = 0;
	      e = 0;
	      while (subs != NULL)
		{
		  sscanf(subs, "%lf", bd + n*strutt[i].block_length + e);
		  subs = strtok(NULL, " ");
		  e++;
		  if (e == strutt[i].block_length)
		    {
		      e = 0;
		      n++;
		    }
		}
	    }
	  else if (strchr(strutt[i].type, 'd'))
	    {
	      bi = (int *) strutt[i].ptr;
	      subs = strtok(stringB, " "); //lista separata da spazi
	      n = 0;
	      e = 0;
	      while (subs != NULL)
		{
		  sscanf(subs, "%d",  bi + n*strutt[i].block_length + e);
		  subs = strtok(NULL, " ");
		  e++;
		  if (e == strutt[i].block_length)
		    {
		      e = 0;
		      n++;
		    }
		}
	    }
	  else if (strchr(strutt[i].type, 's'))
	    {
	      bs = (char *) strutt[i].ptr;
	      subs = strtok(stringB, " "); //lista separata da spazi
	      n = 0;
	      while (subs != NULL)
		{
		  sscanf(subs, "%s", bs + n * strutt[i].block_length);
		  subs = strtok(NULL, " ");
		  n++;
		}
	    }
	  else if (strchr(strutt[i].type, 'c'))
	    {
	      bc = (char *) strutt[i].ptr;
	      subs = strtok(stringB, " "); //lista separata da spazi
	      n = 0;
	      e = 0;
	      while (subs != NULL)
		{
		  sscanf(subs, "%c", bc + n * strutt[i].block_length + e);
		  subs = strtok(NULL, " ");
		  e++;
		  if (e == strutt[i].block_length)
		    {
		      e = 0;
		      n++;
		    }
		  
		}
	      
	    }
	  
	  else 
	    { 
	      mdMsg(ALL, NOSYS, "Parsing", "ERROR", NULL,
		    "Not valid parameter type in singlePar array(in mdsimdep.h).",
		    NULL);
	    }
	  //exit(-1);
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
  //exit(-1);
}

/* ======================== >>> readAsciiPars <<< ========================= */
void readAsciiPars(FILE* pfs, struct pascii strutt[])
{
  char line[256000+NAME_LENGTH];
  char str1[NAME_LENGTH], str2[256000];

  /* scanning the file for the other params */
  while (!feof(pfs))
    {
      /* The syntax must be <parameter>:<value> 
         if <value>  is a '*' then it use the default value or 
         the value loaded from the coordinates file inifile, if
         specified. 
         SPACES ARE IGNORED <---------------------- !!!! */
      fscanf(pfs, "%[^\n] ", line);
      
      if (!strcmp(line, "@@@"))
	break;
      //printf("line: %s\n", line);
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n#] ", str1, str2) < 2)
	continue;
      
      //fscanf(pfs, " %[^: ] : %[^'\n' ]", str1, str2);
      /* analyzes the parameter just read 
         This function depends strongly upon simulation */
      asciiParsing(strutt, str1, str2);
    }
}

/* ======================== >>> readAsciiPars <<< ========================= */
void writeAsciiPars(FILE* fs, struct pascii strutt[])
{
  int i, n, e;
  double* bd;
  int *bi;
  char *bc;
  char *bs;
  
  for  (i=0; strutt[i].ptr != NULL; ++i) /* parname=NULL menas END */
    {
      fprintf(fs, "%s: ", strutt[i].parName);
      if (strchr(strutt[i].type,'f') || strchr(strutt[i].type, 'e')
	  || strchr(strutt[i].type, 'E') || strchr(strutt[i].type, 'g')
	  || strchr(strutt[i].type, 'G'))
	{
	  bd = (double *) strutt[i].ptr;
			
	  for (n = 0; n < strutt[i].qty; n++)
	    {
	      for (e = 0; e < strutt[i].block_length; e++)
		{
		  fprintf(fs, strutt[i].type, *(bd + n*strutt[i].block_length + e));
		  fprintf(fs, " ");
		}
	    }
	}
      else if (strchr(strutt[i].type, 'd'))
	{
	  bi = (int *) strutt[i].ptr;
			
	  for (n = 0; n < strutt[i].qty; n++)
	    {
	      for (e = 0; e < strutt[i].block_length; e++)
		{
		  fprintf(fs, strutt[i].type, *(bi + n*strutt[i].block_length + e));
		  fprintf(fs, " ");
		}
	    }
	  
	}
      else if (strchr(strutt[i].type, 's'))
	{
	  bs = (char *) strutt[i].ptr;
			
	  for (n = 0; n < strutt[i].qty; n++)
	    {
	      fprintf(fs, strutt[i].type, bs + n*strutt[i].block_length);
	      fprintf(fs, " ");
	    }
	    
	}
      else if (strchr(strutt[i].type, 'c'))
	{
	  bc = (char *) strutt[i].ptr;
	  
	  for (n = 0; n < strutt[i].qty; n++)
	    {
	      for (e = 0; e < strutt[i].block_length; e++)
		{
		  fprintf(fs, strutt[i].type, *(bc + n*strutt[i].block_length + e));
		  fprintf(fs, " ");
		}
	    }
	}
      fprintf(fs, "\n");
    }
}

/* ADDED 19/04/2001 */
/* ========================= >>> SaveCorAscii <<< ========================== */
void saveCorAscii(void)
{
  FILE *bf;
  char fn[NAME_LENGTH];
  sync();

  sprintf(fn ,"CorT%.6G-%s.%d", 
	  Oparams.T, 
	  OprogStatus.nRun, Oparams.curStep);

  if ( (bf = fopenMPI(absTmpAsciiHD(fn), "w")) == NULL )
    {
      sprintf(msgStrA, "Problem opening for writing cor file %s ", fn);
      mdMsg(ALL, NOSYS, "saveCorAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
    }

  writeAsciiPars(bf, opar_ascii);

  fprintf(bf, sepStr);
  writeAllCor(bf);  
  
  fclose(bf);
  sync();
}

/* ========================= >>> readCorAscii <<< ========================== */
void readCorAscii(char *fn)
{
  FILE* fs; 
  int SEGSIZE;

  if ((fs = fopenMPI(fn, "r")) == NULL)
    {
      sprintf(msgStrA, "Problem openingn restart file %s ", fn);
      mdMsg(ALL, NOSYS, "ReadBakAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
    }

  readAsciiPars(fs, opar_ascii);
  
  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[0];
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  AllocCoord(SEGSIZE, ALLOC_LISTA, NULL);
    /* allocate  coordinates needed by simulation (see. mdsimul_p) */

  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[1];
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  AllocCoord(SEGSIZE, ALLOC_LISTB, NULL);

  
  readAllCor(fs);
  
  fclose(fs);
}

/* ========================= >>> reaBakAscii <<< =========================== */
void readBakAscii(char* fn)
{
  FILE* fs; 
  if ((fs = fopen(fn, "r")) == NULL)
    {
      sprintf(msgStrA, "Problem opening restart file %s ", fn);
      mdMsg(ALL, NOSYS, "ReadBakAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
    }
  readAsciiPars(fs, opro_ascii);

  readAsciiPars(fs, opar_ascii);

  /* Entrambe queste macro sono definite nel file mono_DPT.h */
  /* read up to coordinates begin */

  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[0];
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  
  AllocCoord(SEGSIZE, ALLOC_LISTA, NULL);
  
  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[1];
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
	 all addresses of the coordinates declared in the simulaiton
	 (see that file) */
  AllocCoord(SEGSIZE, ALLOC_LISTB, NULL);

  readAllCor(fs);
  
  fclose(fs);
}

/* ADDED 04/04/2001 */
/* =========================== >>> SaveBakAscii <<< ========================*/
void saveBakAscii(char* fn)
{
  /* Write restore file, it choose not the last file written (see TECH_INFO) */
  FILE* bf; /* restore file descriptor */
  char fileop[1024];
  //double T;
  //int i;
  /* <---------------------------------------------------- OPEN RESTORE FILE */
  
  /* 6/5/99 ADD
     Before saving the restore do a sync to update all the measure files!!!*/


  sync();

  if (fn == NULL)
    {
      sprintf(fileop ,"CnfT%.6G-%s-%d", 
	      Oparams.T, 
	      OprogStatus.nRun, Oparams.curStep);
      if ( (bf = fopenMPI(absTmpAsciiHD(fileop), "w")) == NULL)
	{
	  mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	  exit(-1);
	}

    }
  else
    {
      strcpy(fileop, fn);
      if ( (bf = fopen(fileop, "w")) == NULL)
	{
	  mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	  exit(-1);
	}
    }

  /* Questo e' l'header nel formato ASCII!!!!!! */
  /*
    fprintf(bf, "%d %d %.5f %.5f\n", Oparams.curStep, Oparams.parnum, 
    Vol, Oparams.T);*/

  
  writeAsciiPars(bf, opro_ascii);
  fprintf(bf, sepStr);
  writeAsciiPars(bf, opar_ascii);
  fprintf(bf, sepStr);

  /* Questa deve essere definita nei file dipendenti dalla simulazione */  
  writeAllCor(bf);
  
  fclose(bf);
 
  sync();/* <--------------------------------------------------------  SYNC */
}

/* =========================== >>> SaveBak <<< =============================*/
void saveBak(char *fileName)
{
  /* Write restore file, it choose not the last file written (see TECH_INFO) */
  int bf; /* restore file descriptor */

  /* <---------------------------------------------------- OPEN RESTORE FILE */
  
  /* 6/5/99 ADD
     Before saving the restore do a sync to update all the measure files!!!*/
  sync();
  bf = mdCreat(fileName, NULL, 
	       "Error during opening restore file for writing.", EXIT);

  /* <--------------------------------------------------------- WRITE DATAS */
  
  mdWrite(bf, NULL, "Error writing the params struct.", EXIT,
	  sizeof(struct params), &Oparams);

  /* write the progStatus struct (see TECH_INFO file) */
  mdWrite(bf, NULL, "Error writing the program status.", EXIT,
	  sizeof (struct progStatus), &OprogStatus);

  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];
  /* writes all arrays to the disk physically making a sync() */
  writeSegs(bf, NULL, "Error writing coordinates array on restore file", EXIT,
	    SEGSIZE, SAVE_LISTA,
	    NULL);             /* NULL means: 'no more pointers to load' */  
  
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];
  writeSegs(bf, NULL, "Error writing coordinates array on restore file", EXIT,
	    SEGSIZE, SAVE_LISTB,
	    NULL);
  /* NOTE:
     Actually the EXT_SLST should be a list of COORD_TYPE variables to save */
  writeSegs(bf, NULL, "Error writing file coordinates on restore file", EXIT, 
	    sizeof(COORD_TYPE), EXT_SLST,
	    NULL);

  mdClose(bf, NULL, "Error closing the restore file", EXIT);

  /* Save on disk the restore file */
  sync();/* <--------------------------------------------------------  SYNC */
}
/* ADDED 09/04/2001: salvataggio di file di restart ascii  */


/* ========================= >>> doubleBufferBak <<< =======================*/
void doubleBufBak(unsigned char* hdWhich, unsigned char* tapeWhich,
		     int* times, int tapeTimes)
{
  /* DESCRIPTION:
     Hard Disk:
      This procedure save the restore datas on file named 
      'COORD_TMP_NAME' + '*hdWhich', that is if COORD_TIME_NAMe = pippo 
      and *hdWhich = 1, then the filename is pippo1.
     Tape:
      Analogous to Hard Disk but it use *tapeWhich switch.
     - times is a pointer to a counter that counts the number of accesses 
     to this procedure; if tapeTimes is != 0 then every 
     (tapeTimes * 2 + 1) * BAK_STEP  steps save on Tape.
     
  */ 
     
    ++(*times);                   /* increment access coounter by 1 */
  /* absTmpHD operator build an absolute name for restore file, 
     that is path + BAK_FILE_NAME ( see mdsimul.h ), while 
     appSw operator append 0 or 1 to the name */
  saveBak(  appSw( absTmpHD(BAK_FILE_NAME), *hdWhich )  );
 
  *hdWhich = sw(*hdWhich);	/* switch buffer for next saving on disk*/
  if ( (*times == OprogStatus.tapeTimes) 
       && (OprogStatus.tapeTimes != 0) )
    {
      *times = 0;               /* reset the access coounter */
      
      /* see above */
      saveBak(  appSw( absTmpTape(BAK_FILE_NAME), *tapeWhich )  );
      *tapeWhich = sw(*tapeWhich);	
      /* switch buffer for next saving on tape*/
    }
}

/* ========================== >>> saveStatus <<< ============================*/
void doubleBufStatus()
{
  /* DESCRIPTION:
     Save status file using STA global switch to ensure it is not the last
     status file saved, in fact at the end of the procedure the 'switch is 
     switched' ( 0->1, 1->0) */ 
 
  int sf; /* sf : file STATUS_FILE (stato della simulaz.) */
  
  /* By appSW and absTmpHD build absolute name for the next status 
     file to save in, for example:
     if STA = 1, STATUS_FILE_NAME = "STATUS_FILE", then 
     file = "/<temporary directory on HD>/STATUS_FILE1" */
  sf = mdCreat( appSw(absTmpHD(STATUS_FILE_NAME), STA), NULL,
		"Simulation continues but it couldn't restart automatically.",
		CONT); 
  
  OsimStat.curStep = Oparams.curStep;
  
  mdWrite(sf, NULL, "Error during writing STATUS_FILE", CONT,
	  sizeof(struct simStat), &OsimStat);
  
  mdClose(sf, NULL, "STATUS_FILE not closed", CONT);
   sync(); /* <---------------------------------       SYNC !!!!!!!!!!!! */
  STA = sw(STA);	/* 1 => 0 or 0 => 1 */  

} 

/* ========================= >>> readFile <<< =========================*/
int readFile(int fdes)
{
  /* DESCRIPTION:
     Read  a file byte by byte, checking it.
     If some problem occurrs then a -1 is reported */
  char byte;
  int result;
 
  /* try read up to the end of file */
  while( (result = read(fdes, &byte, 1)) > 0 );
   
  if (result == -1) 
    {
      return -1; /* if read returns -1 means ERROR */ 
    }
  else 
    {
      return 0;  /* if read return 0 = end of file, that is file read 
		    successfully */
    }
}

/* ========================== >>> saveMeasure <<< ===========================*/
int saveMeasure(int PN, char* fileName, int misNum, char* msgOpen, char* msgWrite, 
		 char* msgClose)
{
  /* DESCRIPTION: 
     This porcedure appends the current measure to the measure file 
     - 'fileName' is the absolute name of the measure file to save to 
     -  misNum is the measure number 
     - 'msgOpen' is the message to print if an open error occurrs
     - 'msgWrite' is the message to print if a write error occurs 
     - 'msgClose' is the message to print if a close error occurrs 
     RETURN VALUE:
      -1 on open/write error
      -2 on close error            
       0 otherwise */
  
  int savedMeasures;
  int mfd1;         /* file desciptor of the file fileName */
  int retValue = 0; /* value to return */ 
  /* save first masure on first datas file */
  if ( (mfd1 = mdOpen(fileName, NULL, msgOpen, CONT | O_WRONLY)) == -1)
    {
      return -1;
    }
  
  /* it is important that the next save is on the measure corrisponding
     to the actual step, so we put the following lseek().
     In fact if we simply append the next measure to the file imagine the 
     following situation: 
     - STEP 10 -> restore file saved
     - STEP 15 -> measure saved
     - STEP 16 system crash 
     When we restart we restart from STEP 10 and we would append STEP 11 
     measure to STEP 15 measure, that is wrong */
  
  /* calculate the number of measure saved, accordingly to the current 
     step number */
  if (OprogStatus.initStep[misNum] == 0)
    savedMeasures = Oparams.curStep / OmeasHead[misNum].saveSteps - 1;
  else
    savedMeasures = Oparams.curStep / OmeasHead[misNum].saveSteps
     - OprogStatus.initStep[misNum];
  
  //  printf("savedMeasure: %d\n", savedMeasures);
  /* -1 because the current measure is not yet saved */

    /* 18/10/2000 MODIFICATO */
  if (PN > 0)
    {
      /* position just after the 'savedMeasure' measures saved */ 
      lseek(mfd1, sizeof(struct measHead) + 
	    OmeasHead[misNum].size * savedMeasures, SEEK_SET);
      /* see TECH_INFO for measure file structure */
    }
  else
    {
      //printf("all'inizio!!!\n");
      lseek(mfd1, sizeof(struct measHead), SEEK_SET);
    }


  /* append to the file the current measure 
     - Omeasure[].buf is a pointer to the current measure
     - Omeasure[].size is the size of the measure */
  if ( mdWrite(mfd1, NULL, msgWrite, CONT,
	       Omeasure[misNum].size, Omeasure[misNum].buf) 
       == -1)
    {
      retValue = -1;
    }

  /* close the file */
  else if (mdClose(mfd1, NULL, msgClose, CONT) == -1)
    {
      retValue = -2; /* -2 = close error */ 
    }

  /* 6/5/99 We wliminate the sync here because it is enough that it does
     the sync when the backup file is saved */
  //sync(); 
  
  //printf("measure saved\n");
  return retValue; /* -1 or -2 or 0 (see above) */ 
}

/* ===================== >>> doubleSaveMeasure <<< ==========================*/
void doubleSaveMeasure(int PN, int misNum, int* times, 
		       int tapeTimes)
{
  /* DESCRIPTION:
     - misNum is the measue number ( 0 <= misNum < NUM_MISURE)
     - whichFirst is the first measure file to save to 
     - if tapeTimes != 0 then it saves on 'tapeDir' (generally a tape)
       every 'tapeTimes * measure.saveSteps' steps 
     - times is a pointer to integer that counts the numeber of access
       to this function */
  char fileA[NAME_LENGTH], fileB[NAME_LENGTH];  /*absolute measure file names*/
  int retValueA, retValueB;                       
  /*  values returned by 'saveMeasure' of 'strTmpA' and 'strTmpB' file */
  
  /* increment the access counter by 1 */
  ++(*times);
  
  /* Build measure file names with absolute path, fileA is the name of the 
     file indicated in whichFirst, fileB the other one */
  
  /* Build fileA, that is an absolute filename
     NOTE: measures are in a directory different to the 
     temporary files dir */
  strcpy(  fileA,
	   appSw( absMisHD(OprogStatus.dataFiles[misNum]), 
		  0 )  ); 

  /* Build fileB, that is the other measure file
     NOTE: sw() return its argument switched */
  strcpy(  fileB,
	   appSw( absMisHD(OprogStatus.dataFiles[misNum]), 
		  1 )  ); 
   
  //printf("0- OmeasHead :%d\n",OmeasHead[misNum].saveSteps);
 
  /* save current measure on both files on hard disk */
  retValueA = saveMeasure(PN, fileA, misNum,
			  "First measure file not openable",
			  "First measure file not writeable",
			  "First measure file not close, anyway datas probably not corrupted");
  retValueB = saveMeasure(PN, fileB, misNum,
			  "Second measure file not openable",
			  "Second measure file not writeable",
			  "Second measure file not close, anyway datas probably not corrupted");
  if ( (retValueA == -1) && (retValueB == -1) ) 
    {
      sprintf(msgStrA, "Measure corrisponding to file %s lost", 
	      OprogStatus.dataFiles[misNum]);
      mdMsg(ALL, NOSYS, NULL, "ALERT", 
	    msgStrA,
	    "I continue anyway...",
	    NULL);
    }
  else if ( (retValueA == -2) && (retValueB == -2) ) 
    {
      sprintf(msgStrA, "Both measure files %s and %s not regularly closed", 
	      fileA, fileB);
      mdMsg(ALL, NOSYS, NULL, "WARNING", NULL,
	    msgStrA,
	    "Datas probably not corrupted",
	    NULL);
    }

  /* ================== >>> save on a backup unit <<< ===================== 
   if tapeTimes!=0, it saves on the media specified by tapeDir*/
  
  /* WARNING: if the number of accesses are tapeTimes save on the tape,we
     doesn't care of saving order for zip */
  if ( (*times == OprogStatus.tapeTimes) && (OprogStatus.tapeTimes != 0) )
    {
      
     
      *times = 0; /* reset the access counter */
   
      pushSteps(&OmeasHead[misNum].saveSteps); 
      /* store Omeasure[misNum].saveSteps value and store in this 
	 fiels the right value  for Tape measure savings*/
     
      //printf("saving on tape measure %f\n", *(int *)Omeasure[0].buf);
      /* build absolute file names for measure files on tape */
      strcpy(  fileA,
	       appSw( absMisTape(OprogStatus.dataFiles[misNum]), 
		      0 )  ); 
      strcpy(  fileB,
	       appSw( absMisTape(OprogStatus.dataFiles[misNum]), 
		      1 )  ); 

      retValueA = saveMeasure(PN, fileA, misNum, 
			      "First measure file not openable on tape",
			      "First measure file not writeable on tape",
			      "First measure file not close, anyway datas probably not corrupted");
      retValueB = saveMeasure(PN, fileB, misNum,
			      "Second measure file not openable",
			      "Second measure file not writeable",
			      "Second measure file on tape not close");
      
      if ( (retValueA == -1) && (retValueB == -1) ) 
	{
	  sprintf(msgStrA, "Measure corrisponding to file %s lost on tape", 
		  OprogStatus.dataFiles[misNum]);
	  mdMsg(ALL, NOSYS, NULL, "ALERT", 
		msgStrA,
		NULL);
	}
      else if ( (retValueA == -2) && (retValueB == -2) ) 
	{
	  sprintf(msgStrA, "Both measure files %s and %s not regularly closed on tape", 
		  fileA, fileB);
	  mdMsg(ALL, NOSYS, NULL, "WARNING", NULL,
		msgStrA,
		"Datas probably not corrupted on tape",
		NULL);
	}
      popSteps(&OmeasHead[misNum].saveSteps); 
      /* restore OmeasHead[i].saveSteps */ 
    }
}

/* =========================== >>> saveOneXva <<< ========================== */
void saveOneXva(char* fileName)
{
  /* DESCRIPTION:
     save positions, velocities and accelerations at regular intervals in
     tape file (xva file).
     This file could be thinked as a measure file, where the single  measure
     contains a set of coordinates for all particles */ 
  
  int savedXva;
  int xvafd; 
  
  xvafd = mdOpen(fileName, NULL, "unble to open xva file (tape file)", 
		 EXIT | O_WRONLY);
  
  //printf("file aperto: %d\n", xvafd);

  /* calculate the number of xva savings, accordingly to the current 
     step number */
  
  /* 12/09/2000 */
  /* savedXva = Oparams.curStep / OxvaHead.saveSteps - 1;*/
  savedXva = OprogStatus.savedXva;
  /* -1 because the current measure is not yet saved */
  /* position just after the 'savedMeasure' measures saved */ 
  lseek(xvafd, sizeof(struct xvaHead) + 
	OxvaHead.size * savedXva, SEEK_SET);
  //printf("xva size: %d\n", OxvaHead.size);
  //printf("scritti %d bytes\n", OxvaHead.size);
  /* if an error occurs exit (see writeSegs proc) */
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[0];

  writeSegs(xvafd, NULL, "Error saving coordinates on tape file", EXIT,
	    SEGSIZE, XVA_LISTA,  
	    NULL); 
  SEGSIZE = sizeof(CORTYPE) * Oparams.parnum[1];
  
  writeSegs(xvafd, NULL, "Error saving coordinates on tape file", EXIT,
	    SEGSIZE, XVA_LISTB,  
	    NULL); 

  /* XVA_LIST is the reduced list of coordinates to save on tape file 
     (xva file) */
  
  mdClose(xvafd, NULL, "Unable to close xva file(tape file)", CONT);
  /* if only a close error occur don't exit */
  
  /* 10/5/99: sync() not necessary here, done when saving restore file */
  //sync(); /* <------------------------------------------------------- SYNC() */
}

/* ============================ >>> saveXva <<< =============================*/
void saveXva(int* times, int tapeTimes)
{
  ++(*times); /* increment the access counter */

  ++OprogStatus.savedXva;
  /* save on Hard Disk */
  
  /* 16/1/1998 SUBS: saveOneXva(absMisHD(OprogStatus.xvafile)); 
     WITH: */
  saveOneXva(absXvaHD(OprogStatus.xvafile));
  
  if ( (*times == OprogStatus.tapeTimes) && (OprogStatus.tapeTimes != 0) )
    {
      *times = 0; /* set the sccess counter to zerp */
      pushSteps(&OxvaHead.saveSteps);
      /* save on tape ( see doubleBufBak for details)*/
      
      /* 16/1/1998 SUBS: saveOneXva(absMisTape(OprogStatus.xvafile)); 
	 WITH: */   
      saveOneXva(absXvaTape(OprogStatus.xvafile));
      
      popSteps(&OxvaHead.saveSteps);
    }
}

/* ========================== >>> pushSteps <<< ==========================*/
int  stepsStack; 
/* "stack" to store the value of saveSteps used to save on HD or Tape*/

void pushSteps(int* saveSteps)
{ 
  /* DESCRIPTION:
     Save '*saveSteps' in 'stepsStack' 
     and store in 'saveSteps' the right value for savings on Tape 
     This procedure is necessary infact suppose for example that
     TAPE_TIMES is 2, on tape we save every 2 * Omeasure[misNum].saveSteps, 
     then If we want a correct header file on tape we must put this value 
     temporarly into Omeasure[misNum].saveSteps.
     After saving we restore the previuos value using popMeasHead()
     (see later) */
  
  stepsStack = *saveSteps;
  *saveSteps = stepsStack * OprogStatus.tapeTimes;

}

/* ========================== >>> popMeasHead <<< ===========================*/
void popSteps(int* saveSteps)
{
  /* DESCRIPTION: restore the value of 'saveSteps' */
  
  *saveSteps = stepsStack;

}

/* ============================= >>> calc<<< ===============================*/
void calc(void)
{
  /* DESCRIPTION:
     Call, if needed,  all "measuring functions" here.
     WARNING: actually calculation of measures are not parallelizeable
       because this proc is called only by the father */
  int i;
  int csteps;
  for (i=0; Omeasure[i].buf != NULL; ++i)
    {
      /* every 'OprogStatus.measCalc[i]' steps calculate the measure, generally
	 this number is equal to 'OprogStatus.measSteps[i]', which is the 
	 number of steps between two savings */ 
      csteps = abs(OprogStatus.measCalc[i]);
      if ( (csteps != 0) && 
	   ((Oparams.curStep % csteps) == 0) &&
	   (Oparams.curStep / csteps >= 
	    OprogStatus.initCalc[i]))
	/* NOTE: the initial step of calculation is given in unit
	   of OprogStatus.measCalc[i] */
	{
	  if (Omeasure[i].calcFunc != NULL) /* NULL = no function */
	    {
	      NUMCALCS = Oparams.curStep / csteps - 
		OprogStatus.initCalc[i] + 1;
	      (*Omeasure[i].calcFunc)(); /* call the "mesuring function" */
	    }
	}
    }
}
