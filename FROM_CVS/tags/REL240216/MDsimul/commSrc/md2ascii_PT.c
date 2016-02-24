/* ========================= >>> md2ascii <<< ================================
   this program convert a measure file containing COORD_TYPE numbers 
   into an ascii file of variuos forms */

#define MD2ASCII
#include<mdsimul.h>

int REPLICA, multi = 0, lambda = 0;
/* se multi e' diversa da zero allora vuol dire
   che ogni misura e' OmeasHaed.PTM misure
   del tipo type cioe' una per ogni 
   replica */
struct measHead OmeasureH; /* measure file header (see mdsimul.h)*/

COORD_TYPE Vol = 1000.0;
/* ========================== >>> CONVERTERS <<< =============================
 >>> ACTUALLY ( see TODO below for future plannings ): 
 - writeFp()
 - writeVec() 
 - writeInt()
*/

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

/* =================================== writeFp =============================*/
void writeFp(FILE* afs,  int step,  int size)
{
  /* DESCRIPTION:
     Convert  COORD_TYPE measure in ascii string of this form:
     <step> <COORD_TYPE value>
    NOTE: size not used in this case.
  */
  if (multi)
    {
      fprintf(afs, sum("%d %.", precision, "f\n", NULL),
	      step, ((COORD_TYPE*)mis)[REPLICA]);
    }
  else
    {
      fprintf(afs, sum("%d %.", precision, "f\n", NULL),
	      step, *((COORD_TYPE*)mis));
    }
}

/* =================================== writeInt =============================*/
void writeInt(FILE* afs,  int step,  int size)
{
  /* DESCRIPTION:
     Convert  'int'  measure in ascii string of this form:
     <step> <int value>
    NOTE: size not used in this case.
  */
  if (multi)
    {
      fprintf(afs, sum("%d %d\n", NULL),
		step, ((int*)mis)[REPLICA]);
      if (lambda)
	fprintf(afs, sum("%d %d\n", NULL),
		step + OmeasureH.saveSteps, ((int*)mis)[REPLICA]);
      
    }
  else
    {
      fprintf(afs, sum("%d %d\n", NULL),
	      step, *(int*)mis);
    }
}
/* =================================== writeVec =============================*/
void writeVec(FILE* afs,  int step,  int size)
{
  /* DESCRIPTION:
     convert an array of COORD_TYPE number to an ascii string os this form:
     <step> <value1> <value2> ... */
  int i;
  /* reformat the single measure opportunely */
  fprintf(afs, "%d ", step);
  for (i = 0; sizeof(COORD_TYPE) * i < size; ++i)
    {
      /* convert the void pointer mis to a COORD_TYPE* pointer, then 
	 dereferenciate it, taking the i-th element of this measure */  
      fprintf(afs, sum("%.", precision, "f ", NULL),
	      ((COORD_TYPE *)mis)[i]);
    }
  fprintf(afs,"\n");
}
/* TODO:
   probably could be useful to have 'converters' that build up an ascii string
   with errors, given a measure with errors */

/* ===================== >>> USER CONVERTERS <<< ========================= */
void writeS(FILE*afs, int step, int size)
{
  /* DESCRIPTION:
     convert an array of COORD_TYPE numbers, containing the structure factor 
     values for variuos n ( k = 2*pi*sqrt(n) ) to an ascii file of the form:
     <2*pi*sqrt(n)> <S[n]>
     .
     .
     & <---- This indicates to xmgr that a new set begin 
   */
  COORD_TYPE pi2, k, scalFact;
  int i;
  pi2 = 4.0*acos(0);
  //scalFact = (KEND - KBEG) / NUMK;
  scalFact = 0.5*pi2/cbrt(Vol);
  for (i=0; i < size / sizeof(COORD_TYPE); ++i)
    {
      k = pi2/cbrt(Vol) + scalFact *((COORD_TYPE)i+0.5);
      fprintf(afs, sum("%.", precision, "f ", "%.", precision, "f ", NULL),
	      k, ((COORD_TYPE *) mis)[i]);
      fprintf(afs,"\n");
    }
  fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
}

/* ============================ >>> writeg <<< ========================== */
void writeg(FILE* afs, int step, int size)
{
  COORD_TYPE r, DELR;
  int i;
  DELR = (REND - RBEG) / MAXBIN;
  for (i=0; i < size / sizeof(COORD_TYPE); ++i)
    {
      r = RBEG + ((COORD_TYPE) i) * DELR; 
      fprintf(afs, sum("%.", precision, "f ", "%.", precision, "f ", NULL),
	      r, ((COORD_TYPE *) mis)[i]);
      fprintf(afs,"\n");
    }
  fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
}

/* ============================ >>> writeMB <<< ======================== */
void writeMB(FILE* afs, int step, int size)
{
  COORD_TYPE v, DELV;
  int i;
  DELV = (VEND - VBEG) / NUMV;
  for (i=0; i < size / sizeof(int); ++i)
    {
      //printf("NUMV:%d, size:%d\n",NUMV, size);
       
      v = VBEG + ((COORD_TYPE) i) * DELV; 
      fprintf(afs, sum("%.", precision, "f ", "%.", precision, "d", NULL),
	      v, ((int *) mis)[i]);
      fprintf(afs,"\n");
    }
  fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
}

/* ============================= >>> writePE <<< =========================== */
void writePE(FILE* afs, int step, int size)
{
  int ss, iE;
  double EN;
  
  rewind(afs);
  for (ss = 0; ss < OmeasureH.PTM; ss++)
    {
      for (iE = 0; iE < OmeasureH.PEpoints; iE++)
	{ 
	  EN = OmeasureH.ENmin + 
	    (((double) iE) / ((double)OmeasureH.PEpoints)) * 
	    (OmeasureH.ENmax - OmeasureH.ENmin); 
	  // if ( ((double *) mis)[ss * PE_POINTS + iE] != 0)
	  //{
	  fprintf(afs, sum("%.", precision, "f ", "%.", precision, "f", NULL),
		  EN, ((double *) mis)[ss * OmeasureH.PEpoints + iE]);
	  fprintf(afs,"\n");
	  //}
	}
      fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
    }
}

/* ============================= >>> writePE <<< =========================== */
void writePEij(FILE* afs, int step, int size)
{
  int ss, iE, nn;
  double EN;
  //rewind(afs);
  for (ss = 0; ss < 2*OmeasureH.PTM; ss++)
    {
      nn = 0;
      for (iE = 0; iE < OmeasureH.PEpoints; iE++)
	{ 
	  EN = OmeasureH.ENmin + 
	    (((double) iE) / ((double)OmeasureH.PEpoints)) * 
	    (OmeasureH.ENmax - OmeasureH.ENmin); 
	  if( ((double *) mis)[ss * OmeasureH.PEpoints + iE] != 0)
	    {
	      nn++;
	      fprintf(afs, sum("%.", precision, "f ", "%E", NULL),
		      EN, ((double *) mis)[ss * OmeasureH.PEpoints + iE]);
	      fprintf(afs,"\n");
	    }
	}
      if (nn)
	fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
    }
}

/* ==========================================================================*/

/* GLOBAL VARIABLES */
char ofm[MAX_M][NAME_LENGTH];
char outFile[NAME_LENGTH];  /* output file (ascii file) name */
char inputFile[NAME_LENGTH];/* input file (measures file) */
char type[128];              /* type of input measure: this determines also 
			       the form of the output file */
/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  printf("Syntax: md2ascii [[-t <type>] [-o <output file>] [-p <precision>]] \n");
  printf(" <measure file>\n");
  printf("Invalid argument!\n");
  exit(-1);
}

/* ============================ >>> args <<< =============================== */
void args(int argc,char **argv)
{
  int i; /* fictitiuos counter */
  
  /* ==== set defaults =====*/
  
  strcpy(type, "float");  
  /* float means floating point, that is COORD_TYPE type set in mdsimedep.h */

  strcpy(outFile, "xmgr.out"); 
  /* this is the defaukt output filename */
  
  strcpy(precision, "10");
  /* this the default nnumber of digits after comma for floating 
     point values */

  /* =======================*/
  i = 0; 
  do /* scan all argument after the first one and before the last one */
    { 
      ++i;                            /* first i is 1, that is 
					 the second arg => OK */
      if (argc == 1) 
	invalArg();
      if (!strcmp(argv[i], "-multi"))
	{
	  /* NOTA BENE:
	     Per ora si applica solo a Double e Interi cioe'
	     a writeFp e writeInt */
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  
	  multi = 1;
	}
      if (!strcmp(argv[i], "-lambda"))
	{
	  /* NOTA BENE:
	     Per ora si applica solo a Double e Interi cioe'
	     a writeFp e writeInt */
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  /* Serve per plottare correttamente il random walk 
	     in temperatura, in tal modo infatti prima di ogni 
	     cambio di lambda plotta anche il lambda precedente */
	  lambda = 1;
	}

      if (!strcmp(argv[i], "-o")) 
	{
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  strcpy(outFile, argv[i]);
	} 
      else if (!strcmp(argv[i], "-t"))
	{
	  if ( ++i  == (argc - 1) )   /* see above */
	    {
	      invalArg(); 
	    }
	  strcpy(type, argv[i]);
	}
      else if (!strcmp(argv[i], "-p"))
	{
	  if ( ++i  == (argc - 1) )   /* see above */
	    {
	      invalArg(); 
	    }
	  strcpy(precision, argv[i]);
	}
      else if (!strcmp(argv[i], "-Vol:"))
	{
	  if ( ++i  == (argc - 1) )   /* see above */
	    {
	      invalArg(); 
	    }
	  Vol=strtod(argv[i], NULL);
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

/* ============================== >>> MAIN<<< ===============================*/
void main(int argc, char** argv)
{
  int mfd;                   /* measure file descriptor */
  int ss, i, ii;                 /* coounters */
  FILE* afs;                 /* ascii file descriptor */
  FILE *afsm[MAX_M];
  char suff[NAME_LENGTH], ofm[MAX_M][NAME_LENGTH];
  char riga[1024];
  args(argc,argv);
  /* open input file for reading and output file for writing */
  if ( (mfd = open(inputFile,O_RDONLY)) == -1 )
    {
      perror("open measures file");
      exit(-1);
    }
  
  /* read header file froma the measure file */
  if (read(mfd, &OmeasureH, sizeof(struct measHead)) ==  -1)
    {
      perror("read header");
      exit(-1);
    } 
 

  if (multi)
    {
      for (ss = 0; ss < OmeasureH.PTM; ss++)
	{
	  //printf("ss:%d\n", ss);
	  sprintf(suff, "_%d", ss);
	  strcpy(ofm[ss], outFile);
	  strcat(ofm[ss], suff);
	  if ( (afsm[ss] = fopen(ofm[ss],"w+")) == NULL)
	    {
	      perror("open ascii file");
	      exit(-1);
	    }
	}
    }
  else
    {
      if ( (afs = fopen(outFile,"w")) == NULL)
	{
	  perror("open ascii file");
	  exit(-1);
	}
    }
  
  /* And now, knowing the size of each measure, allocate memory for it*/
  mis = malloc(OmeasureH.size); /* pointer to a buffer to store on measure */
  
  for (i=1; read(mfd, mis, OmeasureH.size) != 0; ++i)
    {
      /* Each converters must corrispond to a command line arg for -t 
	 option */
      for (ii=0; OconvStruct[ii].converter != NULL; ++ii)/* NULL=end of list */
	{
	  if (!strcmp(type, OconvStruct[ii].type))
	    {
	      /* dereferentiate pointer to converter */
	      if (multi)
		{
		  for (ss = 0; ss < OmeasureH.PTM; ss++)
		    {
		      REPLICA = ss;
		      (*OconvStruct[ii].converter)(afsm[ss], 
						   i * OmeasureH.saveSteps, 
						   OmeasureH.size);
		    }
		}
	      else
		{
		  (*OconvStruct[ii].converter)(afs, 
					       i * OmeasureH.saveSteps, 
					       OmeasureH.size);
		}
	      /* the second number is the step number that corrispond to the 
		 measure, the third arg is the length of the measure 
		 in byes*/
	    break;
	    }
	}
      
      /* If type filed is NULL means that the we have reach 'End of List'
	 without any match => invalid type */
      if (OconvStruct[ii].converter == NULL)
        {
	  printf("Invalid type.\n");
	  exit(-1);
	}
    }
  
  free(mis);
  
  /* close input file and output file */
  close(mfd);

  if (multi)
    {
      if ( (afs = fopen(outFile,"w")) == NULL)
	{
	  perror("unable to open file");
	  exit(-1);
	}
      /* Copia tutti i file in un unico file  */
      for(ss = 0; ss < OmeasureH.PTM; ss++)
	{
	  
	  rewind(afsm[ss]);
	  while (!feof(afsm[ss]))
	    {
	      if (!(fscanf(afsm[ss], "%[^\n] ", riga) == EOF))
		{
		  //fgets(riga, 0, afsm[ss]);
		  //printf("%s\n",riga);
		  fprintf(afs, "%s\n", riga);
		}
	    }
	  if (ss < (OmeasureH.PTM - 1))
	    fprintf(afs, "&\n");
	}
      for(ss = 0; ss < OmeasureH.PTM; ss++)
	{
	  fclose(afsm[ss]);
	  /* cancella tutti i file dei coeff. di diffusione */
	  unlink(ofm[ss]);
	}
      fclose(afs);
    }
  else
    {
      fclose(afs);
    }
}
