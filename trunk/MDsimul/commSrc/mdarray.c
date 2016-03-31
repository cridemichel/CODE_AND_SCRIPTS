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
#ifdef EDHE_FLEX
extern int readBinCoord_heflex(int cfd);
extern void writeBinCoord_heflex(int cfd);
#endif

extern int SEGSIZE;        /* length in bytes of an array */

/* string used to send  messages */
extern char msgStrA[MSG_LEN],msgStrB[MSG_LEN],msgStrC[MSG_LEN]; 
								   
extern unsigned char BAK, STA, BAKT;

#ifdef EDHE_FLEX
extern void writeAllCor(FILE* , int);
#else
extern void writeAllCor(FILE* );
#endif
extern void readAllCor(FILE* );

const char sepStr[] = "@@@\n";
#ifdef ED_PARALL_DD
extern int dd_totBytes;
extern void *dd_coord_ptr;
#endif

#ifdef MD_ALLOC_POLY
/* ============================ >>> AllocCoord <<< ==========================*/
void AllocCoordPoly(int size, COORD_TYPE** pointer, ...)
{
  /* Allocate contigously memmory */
  va_list ap;
  COORD_TYPE** ptrs[MAXVARS];
  COORD_TYPE* sa;
  int offs[MAXVARS];
  int totBytes, totArr; 
  int i, j, num[MAXVARS];
  va_start(ap, pointer);
  
  ptrs[0] = pointer;
  num[0] = va_arg(ap, int); 
  offs[0] = 0;
  totBytes = size*num[0];
#ifdef ED_PARALL_DD
  dd_totBytes = totBytes;
#endif
  for (i=1; (ptrs[i] = va_arg(ap, COORD_TYPE**)) != NULL; ++i)
    { 
      num[i] = va_arg(ap, int);
      totBytes += size*num[i]; /* total bytes to allocate */
      offs[i] = offs[i-1] + size*num[i-1];
      //printf("offs[%d]:%d\n",i , offs[i]);
    }
  //printf("i=%d ptrs[i]=%p\n", i, ptrs[i]);
  totArr = i; /* total number of arrays to allocate */
  
  sa = (COORD_TYPE*)malloc(totBytes);
#ifdef ED_PARALL_DD
  dd_coord_ptr = sa; 
#endif
  //intf("totBytes = %d size: %d\n", totBytes, size);
  for (i = 0; i < totArr; ++i)
    {
      for (j = 0; j < num[i]; j++)
	{
	  ptrs[i][j] =(COORD_TYPE*)(((char*) sa) + offs[i] + size*j);  
	  //printf("pointer value ptrs[%d]:%p, ptrs[%d][%d]:%d\n", i, ptrs[i],i, j,(int)ptrs[i][j]);
	  //printf("*sptr relative %d\n\n",(int) *ptrs[i] - (int)*ptrs[0]);
	  //printf("totArr: %d num[%d]:%d offs[i]:%d\n", totArr, i, num[i], offs[i]);
	}
   }	
  va_end(ap);
}

#endif
/* ============================ >>> AllocCoord <<< ==========================*/
void AllocCoord(int size, COORD_TYPE** pointer, ...)
{
  /* Allocate contigously memmory */
  va_list ap;
  COORD_TYPE** ptrs[MAXVARS];
  int totBytes, totArr; 
  int i;
  va_start(ap, pointer);
  
  ptrs[0] = pointer;
  totBytes = size;
  for (i=1; (ptrs[i] = va_arg(ap, COORD_TYPE**)) != 0; ++i)
    { 
      totBytes += size; /* total bytes to allocate */
    }
#ifdef ED_PARALL_DD
  dd_totBytes = totBytes;
#endif
  totArr = i; /* total number of arrays to allocate */
  
  *ptrs[0] = malloc(totBytes);
#ifdef ED_PARALL_DD
  dd_coord_ptr = *ptrs[0]; 
#endif
 
  //intf("totBytes = %d size: %d\n", totBytes, size);
  for (i = 1; i < totArr; ++i)
    {
      *ptrs[i] =(COORD_TYPE*)((char*) *ptrs[i-1] + size);  
      //intf("pointer value ptrs[%d]:%d\n", i, (int)*ptrs[i]);
      //intf("*sptr relative %d\n\n",(int) *ptrs[i] - (int)*ptrs[0]);
      //intf("totArr: %d\n", totArr);
    }	
  va_end(ap);
}
/* ============================ >>> AllocCoord <<< ==========================*/
void AllocInteger(int size, int** pointer, ...)
{
  /* Allocate contigously memmory */
  va_list ap;
  int** ptrs[MAXVARS];
  int totBytes, totArr; 
  int i;
  va_start(ap, pointer);
  
  ptrs[0] = pointer;
  totBytes = size;
  for (i=1; (ptrs[i] = va_arg(ap, int**)) != 0; ++i)
    { 
      totBytes += size; /* total bytes to allocate */
    }
  totArr = i; /* total number of arrays to allocate */
  
  *ptrs[0] = malloc(totBytes);
  //intf("totBytes = %d size: %d\n", totBytes, size);
  for (i = 1; i < totArr; ++i)
    {
      *ptrs[i] =(int*)((char*) *ptrs[i-1] + size);  
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
  COORD_TYPE** sptr;
  va_start(ap, pointer);
  *pointer = malloc(size);
  while((sptr = va_arg(ap, COORD_TYPE**)) != NULL)
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
  v = (int**) malloc(size1 * sizeof(int*));
  v[0] = (int*) malloc(size1 * size2 * sizeof(int));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}
void FreeMatI(int** v)
{
  free(v[0]);
  free(v);
}
#ifdef MD_ALLOC_POLY
/* =========================== >>> setToZero <<< ========================== */
void setToZeroPoly(COORD_TYPE** ptr, ...)
{
  /* DESCRIPTION:
     This procedure set to zero all the coordinates, passed as arguments */
  va_list ap;
  COORD_TYPE** sptr;
  int i, num, n;
  
  va_start(ap, ptr);
  num = va_arg(ap, int);
 
  for (n = 0; n < num; n++)
    for(i=0; i<Oparams.parnum; ++i)
      ptr[n][i]=0.0;
  while ( (sptr = va_arg(ap, COORD_TYPE**)) != NULL)
    {
      printf("qui sptr: %p\n", sptr);
      num = va_arg(ap, int);
      for (n = 0; n < num; n++)
       for(i=0; i<Oparams.parnum; ++i) 
	 sptr[n][i]=0.0;
    }
  va_end(ap);
}

#endif
/* =========================== >>> setToZero <<< ========================== */
void setToZero(COORD_TYPE* ptr, ...)
{
  /* DESCRIPTION:
     This procedure set to zero all the coordinates, passed as arguments */
  va_list ap;
  COORD_TYPE* sptr;
  int i;
  
  va_start(ap, ptr);
  for(i=0; i<Oparams.parnum; ++i) ptr[i]=0.0;
  while ( (sptr = va_arg(ap, COORD_TYPE*)) != NULL)
    {
       for(i=0; i<Oparams.parnum; ++i) sptr[i]=0.0;
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
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
 /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
#ifdef EDHE_FLEX
  rerr |= -readBinCoord_heflex(cfd);
#endif
#ifdef MD_ALLOC_POLY
  AllocCoordPoly(SEGSIZE, ALLOC_LIST, NULL);
  /* loads all arrays from the file associated with the fdes descriptor */
  rerr |= -readSegsPoly(cfd, "Init", "Error reading coordinates", CONT,
		    SEGSIZE, SAVE_LIST,
		    NULL);    /* NULL means: 'no more pointers to load' */
#else
  AllocCoord(SEGSIZE, ALLOC_LIST,
	     NULL);
 /* loads all arrays from the file associated with the fdes descriptor */
  rerr |= -readSegs(cfd, "Init", "Error reading coordinates", CONT,
		    SEGSIZE, SAVE_LIST,
		    NULL);    /* NULL means: 'no more pointers to load' */
#endif
#ifdef EXT_SLST 
  rerr |= -readSegs(cfd, "Init", "Error reading extra coordinates", CONT,
		    sizeof(COORD_TYPE), EXT_SLST,
		    NULL);    /* NULL means: 'no more pointers to load' */
#endif
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
  SEGSIZE = sizeof(double) * Oparams.parnum;
  /*printf("Using SEGSIZE=%d\n", SEGSIZE);*/
  mdWrite(cfd, NULL, 
	  "Error writing the params struct.", EXIT,
	  sizeof(struct params), &Oparams);
#ifdef EDHE_FLEX
  writeBinCoord_heflex(cfd);
#endif

#ifdef MD_ALLOC_POLY
  /* writes all arrays to the disk physically making a sync() */
  writeSegsPoly(cfd, "End", "Error writing final coordinates", EXIT,
	    SEGSIZE, SAVE_LIST,
	    NULL);
#else
  /* writes all arrays to the disk physically making a sync() */
  writeSegs(cfd, "End", "Error writing final coordinates", EXIT,
	    SEGSIZE, SAVE_LIST,
	    NULL);
#endif
#ifdef EXT_SLST 
  /* NOTE:
     Actually the EXT_SLST should be a list of COORD_TYPE variables to save */
  writeSegs(cfd, "End", "Error writing file coordinates", EXIT, 
	    sizeof(COORD_TYPE), EXT_SLST,
	    NULL);
#endif
  mdClose(cfd, "End", 
	  "Final coordinates are in restore files on HD and Tape", EXIT);

  sync(); /* <-------------------------------------------SYNC() !!!!!!!!!*/
  
  mdMsg(ALL, NOSYS, "End", "NOTICE", NULL,
	"Final coordinates file successfully saved",
	NULL);
}
#ifdef MD_DYNAMIC_OPROG
extern int dyn_alloc_oprog(void);
extern void set_dyn_ascii(void);
#endif 
/* ======================== >>> ReadBak <<< ================================*/ 
int readBak(int bfd)
{	
  int br; /* bytes read by 'read' system call */
  unsigned char rerr = 0;
  /* Read structures params and filenames and coordinates from restore file
     whose file descriptor is 'bfd'. 
     NOTE: nel salvataggio rispettare l'ordine di caricamento <---------!!!!!
     Attualmente non viene effettuato nessun controllo sugli fread <----!!!!!*/
#ifdef MD_DYNAMIC_OPROG
  int size;
#endif 
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
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
  /* SEGSIZE is the size in bytes of an array of coordinates */

  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */

#ifdef MD_ALLOC_POLY
  AllocCoordPoly(SEGSIZE, ALLOC_LIST, NULL);
#else
  AllocCoord(SEGSIZE, ALLOC_LIST,
	     NULL);
#endif  
  /* reads filenames structure ( see TECH_INFO for details) */
  br = mdRead(bfd, "Restore",  "Error reading the program status.", CONT, 
	      sizeof(struct progStatus), &OprogStatus);
#ifdef MD_DYNAMIC_OPROG
  OprogStatus.dyn_alloc_oprog = dyn_alloc_oprog;
  OprogStatus.set_dyn_ascii = set_dyn_ascii;
  OprogStatus.ptr = NULL;
  size = OprogStatus.dyn_alloc_oprog(); 
  //printf("parnum=%d size=%d\n", Oparams.parnum, size);
  br = mdRead(bfd, "Restore",  "Error reading the program dynamic status.", CONT, 
	    size, OprogStatus.ptr);
#endif 
 if (br == -1)
    { 
      return 1; /* 1 = ERROR */
    }
#ifdef EDHE_FLEX
 rerr |= -readBinCoord_heflex(bfd);
#endif

 /* loads all arrays from the file associated with the fdes descriptor */
#ifdef MD_ALLOC_POLY
 rerr |= -readSegsPoly(bfd, "Restore", "Error reading restore file", CONT, 
		   SEGSIZE, SAVE_LIST,
		   NULL);   /* NULL means: 'no more pointers to load' */
#else
 rerr |= -readSegs(bfd, "Restore", "Error reading restore file", CONT, 
		   SEGSIZE, SAVE_LIST,
		   NULL);   /* NULL means: 'no more pointers to load' */
#endif
#ifdef EXT_SLST
 rerr |= -readSegs(bfd, "Restore", 
		   "Error reading extra coords from restore file", 
		   CONT, sizeof(COORD_TYPE), EXT_SLST,
		   NULL);   /* NULL means: 'no more pointers to load' */
#endif
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
  
  for  (i=0; strcmp(strutt[i].parName,""); i++) /* parname="" and ptr=NULL means END */
    {
      /* 16/05/2010: ptr puo' essere nullo anche in caso di allocazione dinamica quindi controllo anche parName che
	 sia diverso dalla stringa vuota "" */
      //if (!strcmp(stringA, "inifile"))
      //	printf("stringA: %s stringB: %s i: %d name: %s %s\n",
      //	     stringA, stringB, i, strutt[i].parName, strutt[i].type);
      //printf("Parsing...%s:%s\n", strutt[i].parName, stringA);
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

	      //mdMsg(ALL, NOSYS, "Parsing", "ERROR", NULL,
	      //	    "Not valid parameter type in singlePar array(in mdsimdep.h).",
	      //	    NULL);
	    }
	  //exit(-1);
	  return;
	}
    }

  /* no one paremeter defined in the singlePar array matches the read one =>
     ERROR */
  sprintf(msgStrA, "Parameter %s is not valid", stringA);
  mdMsg(ALL, NOSYS, "Parsing", "ERROR", NULL,
	msgStrA,
	NULL);
  //exit(-1);
}
void read_parnum(FILE *pfs)
{
  char *line;
  char str1[NAME_LENGTH], *str2;
  int ll, cpos;
  int atc=0;
#if defined(MC_SUS) && defined(MD_DYNAMIC_OPROG)
  int itemsfound;
#endif
 
  cpos = ftell(pfs);
#ifdef MAXPAR
  /* if we have an array of double MAXPAR long we assume here 20 bytes for each element (usually %.15G
     is a sensible choice for double numbers hence 20 bytes is a safe overestimate */
  ll = MAXPAR*20;
#else
  ll = MAXLINELEN;
#endif
  line = malloc((ll+NAME_LENGTH)*sizeof(char));
  str2 = malloc(ll*sizeof(char)); 
  /* scanning the file for the other params */
#if defined(MC_SUS) && defined(MD_DYNAMIC_OPROG)
  itemsfound = 0;
#endif
  while (!feof(pfs))
    {
      /* The syntax must be <parameter>:<value> 
         if <value>  is a '*' then it use the default value or 
         the value loaded from the coordinates file inifile, if
         specified. 
         SPACES ARE IGNORED <---------------------- !!!! */
      fscanf(pfs, "%[^\n] ", line);
       // fscanf(pfs, "%s ", line);
      //printf("%s\n", line);
      if (!strcmp(line, "@@@"))
	{
	  atc++;
	  if (atc==2)
	    break;
	}
#if !(defined(MD_DYNAMIC_OPROG) && defined(MC_SUS))
      if (atc==0)
	continue;
#endif
      //printf("line=%s\n", line);
      //printf("line: %s\n", line);
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n#] ", str1, str2) < 2)
	continue;
      if (!strcmp(str1,"parnum"))
	{
#ifdef BIMIX
	  sscanf(str2,"%d %d", &(Oparams.parnum[0]), &(Oparams.parnum[1]));
	  printf("[readBakAscii()->read_parnum()]: Oparams.parnum=%d %d\n", Oparams.parnum[0], Oparams.parnum[1]);
#else
	  Oparams.parnum = atoi(str2);
	  printf("[readBakAscii()->read_parnum()]: Oparams.parnum=%d\n", Oparams.parnum);
#endif
	  /* 17/05/2010: qui bisognava fare la free! a valgrind non sfugge nulla! */
#if defined (MD_DYNAMIC_OPROG) && defined(MC_SUS)
	  itemsfound++;
#else
	  fseek(pfs, cpos, SEEK_SET);
	  free(line);
	  free(str2);
	  return;
#endif
	}	
#if defined(MD_DYNAMIC_OPROG) && defined(MC_SUS)
      if (!strcmp(str1, "susnmin"))
	{
	  OprogStatus.susnmin = atoi(str2);	  
	  printf("[readBakAscii()->read_parnum()]: Oparams.susnmin=%d\n", OprogStatus.susnmin);
	  itemsfound++;
	}
      if (!strcmp(str1, "susnmax"))
	{
	  OprogStatus.susnmax = atoi(str2);	  
	  printf("[readBakAscii()->read_parnum()]: Oparams.susnmax=%d\n", OprogStatus.susnmax);
	  itemsfound++;
	}
      if (itemsfound==3)
	{
     	  fseek(pfs, cpos, SEEK_SET);
	  free(line);
	  free(str2);
	  return;
	}
#endif
    }
  fseek(pfs, cpos, SEEK_SET);
  free(line);
  free(str2);
#if defined(MC_SUS)
  printf("[WARNING] not all items found (parnum,susnmix,susnmax) not found in ascii file\n");
#else 
  printf("[WARNING] parnum not found in ascii file\n");
#endif
  exit(-1);
}
/* ======================== >>> readAsciiPars <<< ========================= */
void readAsciiPars(FILE* pfs, struct pascii strutt[])
{
  char *line;
  char str1[NAME_LENGTH], *str2;
  int ll;
#ifdef MAXPAR
  /* if we have an array of double MAXPAR long we assume here 20 bytes for each element (usually %.15G
     is a sensible choice for double numbers hence 20 bytes is a safe overestimate */
  ll = MAXPAR*20;
#else
  ll = MAXLINELEN;
#endif
  line = malloc((ll+NAME_LENGTH)*sizeof(char));
  str2 = malloc(ll*sizeof(char)); 
  /* scanning the file for the other params */
  while (!feof(pfs))
    {
      /* The syntax must be <parameter>:<value> 
         if <value>  is a '*' then it use the default value or 
         the value loaded from the coordinates file inifile, if
         specified. 
         SPACES ARE IGNORED <---------------------- !!!! */
      fscanf(pfs, "%[^\n] ", line);
       // fscanf(pfs, "%s ", line);
      //printf("%s\n", line);
      if (!strcmp(line, "@@@"))
	break;
      //printf("line=%s\n", line);
      //printf("line: %s\n", line);
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n#] ", str1, str2) < 2)
	continue;
      //fscanf(pfs, " %[^: ] : %[^'\n' ]", str1, str2);
      /* analyzes the parameter just read 
         This function depends strongly upon simulation */
      asciiParsing(strutt, str1, str2);
      /*      if (!strcmp(str1,"M"))
	{
	  printf("Oparams.M:%d\n", Oparams.M);
	  printf("str2:%s\n", str2);
	  }*/
    }
  free(line);
  free(str2);
}
#ifdef EDHE_FLEX
extern int is_valid_parname_progStatus(char *pn);
extern int is_valid_parname_params(char *pn);
extern int is_to_save(int i);
#endif
/* ======================== >>> readAsciiPars <<< ========================= */
void writeAsciiPars(FILE* fs, struct pascii strutt[])
{
  int i, n, e;
  double* bd;
  int *bi;
  char *bc;
  char *bs;
#ifdef MD_DYNAMIC_OPROG
  for  (i=0; strutt[i].qty!=0 && strlen(strutt[i].type)!=0; ++i) 
#else 
  for  (i=0; strutt[i].ptr != NULL; ++i) /* parname=NULL means END */
#endif
    {
      if (strutt[i].ptr == NULL)
	continue;
#ifdef EDHE_FLEX
      if (OprogStatus.stripStore && strutt == opro_ascii)
	{
	  if (!is_valid_parname_progStatus(strutt[i].parName))
	    continue;
	}
      if (OprogStatus.stripStore && strutt == opar_ascii)
	{
	  if (!is_valid_parname_params(strutt[i].parName))
	    continue;

	}
#endif
      fprintf(fs, "%s: ", strutt[i].parName);
      if (strchr(strutt[i].type,'f') || strchr(strutt[i].type, 'e')
	  || strchr(strutt[i].type, 'E') || strchr(strutt[i].type, 'g')
	  || strchr(strutt[i].type, 'G'))
	{
	  bd = (double *) strutt[i].ptr;
	  /* N.B. se nel .h si mette -MAXPAR come qty (vedi struct pascii)
	   * allora vengono salvati solo Oparams.parnum elementi risparmiando disco.*/
#ifdef MC_SUS
	  if (strutt[i].qty == -MAXSUSWINDOW)		
	    strutt[i].qty = OprogStatus.susnmax-OprogStatus.susnmin+1;

#endif
	  if (strutt[i].qty == -MAXPAR)		
	    strutt[i].qty = Oparams.parnum;
	  for (n = 0; n < strutt[i].qty; n++)
	    {
#ifdef EDHE_FLEX
	      if ((strutt[i].qty==Oparams.parnum || strutt[i].qty==MAXPAR) 
		  && !globSaveAll)
		{
		  if (!is_to_save(n))
		    continue;
		}
#endif
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
#ifdef EDHE_FLEX
extern int globSaveAll;
#endif
/* ========================= >>> SaveCorAscii <<< ========================== */
void saveCorAscii(void)
{
#ifdef EDHE_FLEX
  int globSaveAllBak, stripStoreBak;
#endif
  FILE *bf;
  char fn[NAME_LENGTH];
  sync();
  if (mgl_mode==2)
    { 
#ifdef POVRAY
      if (!strcmp(inifile_for_mgl,"*"))
	sprintf(fn ,"initconf.pov");
      else		
	sprintf(fn ,"%s.pov", inifile_for_mgl);
#else
      if (!strcmp(inifile_for_mgl,"*"))
	sprintf(fn ,"initconf.mgl");
      else		
	sprintf(fn ,"%s.mgl", inifile_for_mgl);
#endif
    }
  else
    {
#ifdef EDHE_FLEX
      sprintf(fn ,"CorFinal");
#else
#ifdef MDLLINT 
      sprintf(fn ,"CorT%.6G_%s_%lld", 
	      Oparams.T, 
	      OprogStatus.nRun, (long long int)Oparams.curStep);
#else
      sprintf(fn ,"CorT%.6G_%s_%d", 
	      Oparams.T, 
	      OprogStatus.nRun, Oparams.curStep);
#endif
#endif
    }
  if ( (bf = fopenMPI(absTmpAsciiHD(fn), "w")) == NULL )
    {
      sprintf(msgStrA, "Problem opening for writing cor file %s ", fn);
      mdMsg(ALL, NOSYS, "saveCorAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
    }
#ifdef EDHE_FLEX
  stripStoreBak = OprogStatus.stripStore;
  OprogStatus.stripStore = 0;
  globSaveAllBak = globSaveAll;
  globSaveAll = 1;
#endif

  if (mgl_mode==0)
    {
      writeAsciiPars(bf, opar_ascii);
      fprintf(bf, sepStr);
    }
#ifdef EDHE_FLEX
  writeAllCor(bf,1);
  OprogStatus.stripStore = stripStoreBak;
  globSaveAll = globSaveAllBak;
#else
  writeAllCor(bf);
#endif  
#ifdef EDHE_FLEX 
  OprogStatus.stripStore = stripStoreBak;
  globSaveAll = globSaveAllBak;
#endif

  fclose(bf);
  sync();
}

/* ========================= >>> readCorAscii <<< ========================== */
void readCorAscii(char *fn)
{
  FILE* fs; 

  if ((fs = fopenMPI(fn, "r")) == NULL)
    {
      sprintf(msgStrA, "Problem opening restart file %s ", fn);
      mdMsg(ALL, NOSYS, "ReadBakAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
    }

  readAsciiPars(fs, opar_ascii);
  
  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
  printf("SEGSIZE:%d\n", SEGSIZE);
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
#ifdef MD_ALLOC_POLY
  AllocCoordPoly(SEGSIZE, ALLOC_LIST, NULL);
#else
  AllocCoord(SEGSIZE, ALLOC_LIST, NULL);
#endif
  readAllCor(fs);
  
  fclose(fs);
}

/* ========================= >>> reaBakAscii <<< =========================== */
void readBakAscii(char* fn)
{
  FILE* fs; 
  if ((fs = fopenMPI(fn, "r")) == NULL)
    {
      sprintf(msgStrA, "Problem opening restart file %s ", fn);
      mdMsg(ALL, NOSYS, "ReadBakAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
      exit(-1);
    }
  /* N.B. 16/05/2010: prima parnum non era stato assegnato quando si chiamava readAsciiPars()!
     read_parnum() legge il numero di particelle dal file di restart ascii. */
  read_parnum(fs);
#ifdef MD_DYNAMIC_OPROG
  /* N.B. 16/05/2010: alloca la memoria dinamica per OprogStatus e inizializza i puntatori della struttura
     opro_ascii con i puntatori ottenuti */	
  OprogStatus.dyn_alloc_oprog();
#endif
  readAsciiPars(fs, opro_ascii);
  readAsciiPars(fs, opar_ascii);
 
  /* Entrambe queste macro sono definite nel file mono_DPT.h */
  /* read up to coordinates begin */

  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
#ifdef MD_ALLOC_POLY
  AllocCoordPoly(SEGSIZE, ALLOC_LIST, NULL);
#else
  AllocCoord(SEGSIZE, ALLOC_LIST, NULL);
#endif
  readAllCor(fs);

  fclose(fs);
}

/* ADDED 04/04/2001 */
/* =========================== >>> SaveBakAscii <<< ========================*/
void saveBakAscii(char *fn)
{
  /* Write restore file, it choose not the last file written (see TECH_INFO) */
  FILE* bf; /* restore file descriptor */
  char fileop[1024], fileop2[1024];

  //double T;
  //int i;
  /* <---------------------------------------------------- OPEN RESTORE FILE */
  
  /* 6/5/99 ADD
     Before saving the restore do a sync to update all the measure files!!!*/
  sync();
  
  if (mgl_mode==2)
    {
#ifdef POVRAY
      if (!strcmp(inifile_for_mgl,"*"))
	sprintf(fileop2 ,"initconf.pov");
      else		
	sprintf(fileop2 ,"%s.pov", inifile_for_mgl);

#else
      if (!strcmp(inifile_for_mgl,"*"))
	sprintf(fileop2 ,"initconf.mgl");
      else		
	sprintf(fileop2 ,"%s.mgl", inifile_for_mgl);
#endif
    }
  else
    {
      if (fn == NULL)
	{

#ifdef MDLLINT
	  sprintf(fileop2 ,"CnfT%.6G_%s_%lld", 
		  Oparams.T, 
		  OprogStatus.nRun, Oparams.curStep);
#else
	  sprintf(fileop2 ,"CnfT%.6G_%s_%d", 
		  Oparams.T, 
		  OprogStatus.nRun, Oparams.curStep);
#endif
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	}
      else
	{
	  strcpy(fileop, fn);
	}
    }
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  /* Questo e' l'header nel formato ASCII!!!!!! */
  /*
    fprintf(bf, "%d %d %.5f %.5f\n", Oparams.curStep, Oparams.parnum, 
    Vol, Oparams.T);*/
  if (mgl_mode == 0)
    {
      writeAsciiPars(bf, opro_ascii);
      fprintf(bf, sepStr);
      writeAsciiPars(bf, opar_ascii);
      fprintf(bf, sepStr);
    }

  /* Questa deve essere definita nei file dipendenti dalla simulazione */  
#ifdef EDHE_FLEX
  writeAllCor(bf, 1);
#else
  writeAllCor(bf);
#endif
  fclose(bf);
  sync();/* <--------------------------------------------------------  SYNC */
}

/* =========================== >>> SaveBak <<< =============================*/
void saveBak(char *fileName)
{
  /* Write restore file, it choose not the last file written (see TECH_INFO) */
  int bf; /* restore file descriptor */
#ifdef MD_DYNAMIC_OPROG
  int size;
#endif
  /* <---------------------------------------------------- OPEN RESTORE FILE */
  
  /* 6/5/99 ADD
     Before saving the restore do a sync to update all the measure files!!!*/
  sync();
  bf = mdCreat(fileName, NULL, 
	       "Error during opening restore file for writing.", EXIT);

  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
  /* <--------------------------------------------------------- WRITE DATAS */
  
  mdWrite(bf, NULL, "Error writing the params struct.", EXIT,
	  sizeof(struct params), &Oparams);
  /* write the progStatus struct (see TECH_INFO file) */
  mdWrite(bf, NULL, "Error writing the program status.", EXIT,
	  sizeof (struct progStatus), &OprogStatus);
  //printf("sizeof progStatus=%d\n", sizeof(struct progStatus));
#ifdef MD_DYNAMIC_OPROG
  size = OprogStatus.len; 
  mdWrite(bf, NULL,  "Error writing the program dynamic status.", EXIT, 
	  size, OprogStatus.ptr);
#endif 
 
#ifdef EDHE_FLEX
  writeBinCoord_heflex(bf);
#endif

  /* writes all arrays to the disk physically making a sync() */
#ifdef MD_ALLOC_POLY
  writeSegsPoly(bf, NULL, "Error writing coordinates array on restore file", EXIT,
	    SEGSIZE, SAVE_LIST,
	    NULL);             /* NULL means: 'no more pointers to load' */  
#else
  writeSegs(bf, NULL, "Error writing coordinates array on restore file", EXIT,
	    SEGSIZE, SAVE_LIST,
	    NULL);             /* NULL means: 'no more pointers to load' */  
#endif
#ifdef EXT_SLST 
  /* NOTE:
     Actually the EXT_SLST should be a list of COORD_TYPE variables to save */
  writeSegs(bf, NULL, "Error writing file coordinates on restore file", EXIT, 
	    sizeof(COORD_TYPE), EXT_SLST,
	    NULL);
#endif
  mdClose(bf, NULL, "Error closing the restore file", EXIT);

  /* Save on disk the restore file */
  sync();/* <--------------------------------------------------------  SYNC */
}
extern double *rx, *ry, *rz, ***R;
/* ========================= >>> doubleBufferBak <<< =======================*/
#ifdef EDHE_FLEX
extern void R2u(void);
#endif
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
#ifdef EDHE_FLEX
    R2u();
#endif
    saveBak(appSw( absTmpHD(BAK_FILE_NAME), *hdWhich )  );
#ifdef MD_COORDTMP_ASCII
    MD_COORDTMP_ASCII(*hdWhich);
#endif

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

  /* 23/04/2001 MODIFICATO */
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

  /* position just after the 'savedMeasure' measures saved */ 
  //lseek(mfd1, sizeof(struct measHead) + 
  //OmeasHead[misNum].size * savedMeasures, SEEK_SET);
  /* see TECH_INFO for measure file structure */

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
  
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
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
#ifdef MD_ALLOC_POLY
  writeSegsPoly(xvafd, NULL, "Error saving coordinates on tape file", EXIT,
	    SEGSIZE, XVA_LIST,  
	    NULL); 
#else
  writeSegs(xvafd, NULL, "Error saving coordinates on tape file", EXIT,
	    SEGSIZE, XVA_LIST,  
	    NULL); 
#endif
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
#ifdef MDLLINT
long long int stepsStack;
#else
int  stepsStack; 
#endif
/* "stack" to store the value of saveSteps used to save on HD or Tape*/
#ifdef MDLLINT
void pushSteps(long long int* saveSteps)
#else
void pushSteps(int* saveSteps)
#endif
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
#ifdef MDLLINT
void popSteps(long long int* saveSteps)
#else
void popSteps(int* saveSteps)
#endif
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
  for (i=0; Omeasure[i].buf != NULL; ++i)
    {
      /* every 'OprogStatus.measCalc[i]' steps calculate the measure, generally
	 this number is equal to 'OprogStatus.measSteps[i]', which is the 
	 number of steps between two savings */ 
      if ( (OprogStatus.measCalc[i] != 0) && 
	   ((Oparams.curStep % OprogStatus.measCalc[i]) == 0) &&
	   (Oparams.curStep / OprogStatus.measCalc[i] >= 
	    OprogStatus.initCalc[i]))
	/* NOTE: the initial step of calculation is given in unit
	   of OprogStatus.measCalc[i] */
	{
	  if (Omeasure[i].calcFunc != NULL) /* NULL = no function */
	    {
	      NUMCALCS = Oparams.curStep / OprogStatus.measCalc[i] - 
		OprogStatus.initCalc[i] + 1;
	      (*Omeasure[i].calcFunc)(); /* call the "measuring function" */
	    }
	}
    }
}




