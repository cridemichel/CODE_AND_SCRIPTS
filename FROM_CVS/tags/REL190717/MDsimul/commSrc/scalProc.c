#include<mdsimul.h>

/* GLOBAL VARIABLES */
extern char outFile[NAME_LENGTH];  /* output file (ascii file) name */
extern char inputFile[NAME_LENGTH];/* input file (measures file) */
extern char posFile[NAME_LENGTH];  /* name of the file with positions inside*/
extern int SEGSIZE, posBool, corBool, tBool, outBool, ihdr, ohdr, 
  Nm, infoBool, boxcmBool;

/* strings to store messages (use calling mdMsg()) */
extern COORD_TYPE T;
extern COORD_TYPE DECL_LIST;
extern COORD_TYPE EXT_DLST;

/* ============================ >>> help <<< ============================= */
void help(void)
{
  printf("Syntax: scaling [[-T: <temperature>] [-o <output file>]\n");
  printf(" [-Nm: <number of particles>] [-nooh] [-noih] [-cor] [-infos] [-h]]\n");
  printf(" <coordinate file>\n");
}

/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  help();
  printf("Invalid argument!\n");
  exit(-1);
}

/* =============================== >>> args <<< =============================*/
void args(int argc,char **argv)
{
  int i; /* fictitiuos counter */
  
  /* ==== set defaults =====*/
  
  T = 1.0; /* default temperature */

  strcpy(outFile, "scaled.cor"); 
  /* this is the defaukt output filename */

  /* =======================*/
  i = 0; 
  if (argc == 1) 
    invalArg();
  do /* scan all argument after the first one and before the last one */
    { 
      ++i;                            /* first i is 1, that is 
					 the second arg => OK */
      if (!strcmp(argv[i], "-h"))
	{
	  help();
	  exit(-1);
	}
      else if (!strcmp(argv[i], "-boxcm")) 
	{
	  boxcmBool = 1;
	}
      else if (!strcmp(argv[i], "-o")) 
	{
	  outBool = 1;
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  strcpy(outFile, argv[i]);
	} 
      else if (!strcmp(argv[i], "-pf:"))
	{
	  posBool = 1;
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  strcpy(posFile, argv[i]);
	}
      else if (!strcmp(argv[i], "-Nm:"))
	{
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  Nm = atoi(argv[i]);
	  
	}
      else if (!strcmp(argv[i], "-cor"))
	{
	  corBool = 1;
	}
      else if (!strcmp(argv[i], "-noih"))
	{
	  /* In this way the header is not read from the input file, that is 
	     the input file should be a raw coordinate file (raw = 
	     without header) */
	  ihdr = 0;
	}
      else if (!strcmp(argv[i], "-nooh"))
	{
	  /* omit the header in the output coordinate file */
	  ohdr = 0;
	}
      else if (!strcmp(argv[i], "-T:")) 
	{
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  tBool = 1;
	  T = atof(argv[i]); /* Temperature to switch to */
	} 
      else if (!strcmp(argv[i], "-infos"))
	{
	  /* Print all the infos contained in the header of the file */
	  infoBool = 1;
	}
      else                            /* option not existant */ 
	{
	  invalArg();
	}
    } while ( i < (argc - 2) );
  /* and now take the input file */
  
  /* there is another arg ?  if N => ERROR */
  if (++i == argc)  
    {
      invalArg();
    }
 
  /* last arg that is input file name */
  strcpy(inputFile, argv[i]);
}

/* ========================= >>> AllocCoord <<< ========================= */
void allocCor(int size, COORD_TYPE** pointer, ...)
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

/* ======================== >>> loadSegs  <<< ===============================*/
void rSegs(int fdes, int size, void* pointer, ...)
{
  void* sptr;
  va_list ap;
  va_start(ap, pointer);
  
  /* Load first pointerof the list */
  if (read(fdes, pointer, size) < size) 
    {
      printf("ERROR: Too few bytes read!");
      exit(-1);
    }
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      if ( read(fdes, sptr, size) < size )
	{
	  printf("ERROR: Too few bytes read!");
	  exit(-1);
	}
    }
  va_end(ap);
}

/* =========================== >>> loadCor <<< ============================= */
void loadCor(int cfd)
{
  /* DESCRIPTION:
     Read structure params and coordinates from file whose file descriptor is 
     cfd.
     NOTE: nel salvataggio rispettare l'ordine di caricamento <----!!!!!
     Attualmente non viene effettuato nessun controllo sugli fread <----!!!!
   */
  if (ihdr == 1)
    {
      if (read(cfd, &Oparams,  sizeof(struct params)) < sizeof(struct params))
	{
	  printf("ERROR: Too few bytes read in the header!");
	  exit(-1);
	}
      Nm = Oparams.parnum;
    }
  else if (Nm == 0)
    {
      printf("ERROR: You want to discard the header but you don't specify the particles number\n");
      scanf("Particles Number? %d", &Nm);
      printf("OK thanks\n");
    }
 
  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Nm;
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  allocCor(SEGSIZE, ALLOC_LIST,
	   NULL);

  /* loads all arrays from the file associated with the fdes descriptor */
  rSegs(cfd, SEGSIZE, SAVE_LIST,
	NULL);    
  /* NULL means: 'no more pointers to load' */
#ifdef EXT_SLST 
  /* loads extra variable such as the friction coefficents and its derivatives
     in the Hoover isotermal method */
  rSegs(cfd, sizeof(COORD_TYPE), EXT_SLST,
	NULL);    /* NULL means: 'no more pointers to load' */
#endif
}

/* ======================== >>> saveSegs  <<< ===============================*/
void wSegs(int fdes, int size, void* pointer, ...)
{
  void* sptr;
  va_list ap;
  va_start(ap, pointer);
  write(fdes, pointer, size); 
  /* at least one pointer should be present */
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      /* if mode == EXIT if an error occurs then exit */
      write(fdes, sptr, size);
    }
  va_end(ap);
}
/* =========================== >>> saveCor <<< ============================*/
void saveCor(char* fileName)
{
  /* DESCRIPTION:
     save last coordinates on file fileName */
 
  int cfd; /* descriptor of coordinate file named fileName */ 
 
  cfd = creat(fileName, 0666);
  if (ohdr == 1)
    {
      write(cfd, &Oparams, sizeof(struct params));
    }
  /* writes all arrays to the disk physically making a sync() */
  //printf("SEGSIZE: %d\n", SEGSIZE);
  wSegs(cfd, SEGSIZE, SAVE_LIST,
	NULL);
#ifdef EXT_SLST
  /* Write the extra variable (COORD_TYPE) */
  wSegs(cfd, sizeof(COORD_TYPE), EXT_SLST,
	NULL);
#endif
  close(cfd);

}
