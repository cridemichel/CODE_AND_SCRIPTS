#include <mdsimul.h>
#include <xva2ascii.h>
COORD_TYPE XVA_DLST;
/* 27/03/2001 TODO:
   oltre a questa routine bisognerebbe farne anche una che converte 
   un file di restart ascii in un file di coordinate o in un bak file */


struct xvaHead OxvaH; /* measure file header (see mdsimul.h)*/
int njob = 0;
char filename[10000][512];


/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  printf("Syntax: jump [-p <precision>] -f <xva_parameters_file>\n");
  printf("Invalid argument!\n");
  exit(-1);
}

/* ============================== >>> Parsing <<< ===========================*/
void xvaParsing(char stringA[NAME_LENGTH], char stringB[NAME_LENGTH])
{
  int i;
  char msgStrA[512];
  /* check if 'stra', read from the params file, is one of the strings in the 
     singlePar array, if N => ERROR */ 
  for (i=0; OxvaPar[i].ptr != NULL; ++i) /* parname=NULL menas END */
    {
      if (!strcmp(OxvaPar[i].parName, stringA))
	{
	  switch (OxvaPar[i].type)
	    {
	    case STR : 
	      strcpy((char *) OxvaPar[i].ptr, stringB);
	      break;
	    case INT :
	      *((int*) OxvaPar[i].ptr) = atoi(stringB);
	      break;
	    case CT :
	      ( *((COORD_TYPE*)OxvaPar[i].ptr) ) = 
		(COORD_TYPE) atof(stringB); 
	      break;
	    default: 
	      printf("Not valid parameter type in singlePar array(in mdsimdep.h).");
	      exit(-1);
	    }
	  /* successfully read prameter */
	  return;
	}
    }
   /* no one paremeter defined in the singlePar array matches the read one =>
     ERROR */
  sprintf(msgStrA, "Paramater %s is not valid", stringA);
  printf( "%s\n", msgStrA);
  exit(-1);
}

/* ========================== >>> readPars <<< ========================== */
void readPars(char* fileName)
{
  FILE* pfs;
  char str1[255], str2[255];
  char line[2048];

  if ( (pfs = fopen(fileName, "r")) == NULL)
    {
      printf("ERROR: Unable to open the parameters file %s\n", fileName);
      exit(-1);
    }

  strcpy(str1, "");
  strcpy(str2, "");

  while(!feof(pfs))
    {
      fscanf(pfs, "%[^\n] ", line);
      
      if (!strcmp(line, "")) /* If a void line */
	continue;
      
      if (sscanf(line, "%[^:# ] : %[^\n# ] ", str1, str2) < 2)
	continue;
      
      //printf("s1:*%s* s2:*%s*\n", str1, str2); /* DEBUG */
      if (!strcmp(str2, "*")) return; /* '*' = default value */
      xvaParsing(str1, str2); /* Check the new parameter */
    }
  fclose(pfs);
}

/* ========================== >>> defaults <<< ============================= */
void defaults(void)
{
  PTM = 0;
  strcpy(inputFile,"");
  strcpy(tempStr, "");
  nRun = 1;

}

/* ============================ >>> args <<< =============================== */
void args(int argc,char **argv)
{
  int i; /* fictitiuos counter */
  
  /* =======================*/
  i = 0; 
  do /* scan all argument after the first one and before the last one */
    { 
      if (i+1 == argc) 
	{
	  invalArg();
	}

      ++i;    
                        
      if (!strcmp(argv[i], "-p"))
	{
	  if ( i+1  == argc )   /* see above */
	    {
	      invalArg(); 
	    }
	  ++i;
	  strcpy(precision, argv[i]);
	}

      else if(!strcmp(argv[i], "-f")) 
	{
	  if ( (i+1) == argc )
	    {
	      invalArg();
	    }
	  ++i;
	  strcpy(xvaparsFile, argv[i]);
	}
      else                            /* option not existant */ 
	{
	  invalArg();
	}
    } while ( i < (argc - 1) );
  /* and now take the input file */
  return;

  /* there is another arg ?  if N => ERROR */
  if (++i == argc)  
    {
      invalArg();
    }
  
  /* last arg that is input file name */
  strcpy(inputFile, argv[i]);
  printf("read parameter file!\n");
}



/* ======================= >>> openMeasFile <<< ============================ */
int openMeasFile(char* FileName)
{
  int mfd;
  /* open input file for reading and output file for writing */

  if ( (mfd = open(FileName,O_RDONLY)) == -1 )
    {
      perror("open measures file");
      printf("ERROR: Unable to open input file %s\n", FileName);
      exit(-1);
    }

  /* read header file froma the measure file */
  if (read(mfd, &OxvaH, sizeof(struct xvaHead)) ==  -1)
    {
      perror("read header");
      exit(-1);
    } 

  return mfd;
}

/* =========================== >>> readMeas <<< ============================*/
int readMeas(int fd, int nmeas, COORD_TYPE* firstPtr, ...)
{
  /* DESCRIPTION:
     read the nmeas measure with nmeas >= 0 */
  int begByte, rs; 
  COORD_TYPE* sptr;
  va_list ap;
  va_start(ap, firstPtr);

  begByte = sizeof(struct xvaHead);
  lseek(fd, begByte + nmeas * OxvaH.size, SEEK_SET);

  rs = read(fd, firstPtr, OxvaH.size / XVA_NUM);

  if (rs == -1)   
    {
      printf("ERROR: Unable to read xva file\n");
      perror("read");
      exit(-1);
    }
  else if (rs < (OxvaH.size / XVA_NUM))
    {
      return 0; /* 0 = end of file */
    }

  while (  (sptr = va_arg(ap, COORD_TYPE*)) != NULL )
    {
      rs = read(fd, sptr, OxvaH.size / XVA_NUM);
      
      if (rs == -1)   
	{
	  printf("ERROR: Unable to read xva file\n");
	  perror("read");
	  exit(-1);
	}
      else if (rs < (OxvaH.size / XVA_NUM))
	{
	  return 0; /* 0 = end of file */
	}
    }

  va_end(ap);
  return 1; /* 1 == OK */
}

/* ========================== >>> Allocate <<< ============================ */ 
void Allocate(int totSize, COORD_TYPE** pointer, ...)
{
  va_list ap;
  COORD_TYPE** sptr;
  int size;

  size = totSize / XVA_NUM;
  //printf("totSize: %d, size: %d\n", totSize, size);
  /* All the variables related to the current measure are contigously 
     allocated */
  va_start(ap, pointer);
  *pointer = malloc(size);
  
  while((sptr = va_arg(ap, COORD_TYPE**)) != NULL)
    {
      *sptr = malloc(size);
    }
  va_end(ap);
}

/* =========================== >>> aopen <<< =============================== */
FILE* aopen(char *fileName)
{
  FILE *ofs;
  if ( (ofs = fopen(fileName, "w")) == NULL)
    {
      perror("open file");
      printf("ERROR: unable to open %s file\n", fileName);
      exit(-1);
    }
  return ofs;
}

/* ============================= >>> aclose <<< ===========================*/
void aclose(FILE* ofs)
{
  fclose(ofs);
}

/* =========================== >>> mkFormat <<< ========================= */ 
void mkFormat(char* fmtStr)
{
  strcpy(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f");
  strcat(fmtStr, " ");
  strcat(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f\n");
}


/* =========================== >>> genFiles <<< ========================== */
void genFiles(int fd, int Nm, int i)
{
  int t0, indice; /* indice individua univocamente uno snaoshot dati una 
		     una temperatura e nRun */
  char fn[1024];
  FILE* fout;
 
  indice = 0;
  for (t0 = 0; ;t0 = t0 + tgap)
    {
      njob++;
      indice++;
      sprintf(fn, "CnfT%.4f-%d-n.%d", Temp[i], nRun, indice);
      /* L'array filename contiene i nomi di tutti i files da processare */
      strcpy(filename[njob], fn);

      if (readMeas(fd, t0, XVA_LIST, NULL) == 0)
	{
	  return;// end of file
	}
      fout = aopen(fn);
      //for (i = 0; i < Nm)
      aclose(fout);
    }
}

/* =========================== >>> MD2Fra <<< ============================= */
void MD2Fra(void)
{
  int mfd, i;
  int FirstTime = 1;
  FILE* fout;
  char tmpStr[512];
  /*
    NOTA: Converte i files XVA in tanti files (uno per configurazione),
    questo perché la routine di Francesco richiede tanti diversi files
    da leggere.
  */
  /* Se PTM > 0 vuol dire che ho a che fare con una simulazione fatta con 
     PARALLEL TEMPERING ed un file xva in ingresso per ogni temperatura
     mentre se PTM <= 0  ai assume che in ingresso si abbia un solo
     file xva */ 
  njob = 0;
   
  if (PTM > 0)
    {
      for (i = 0; i < PTM; i++)
	{
	  sprintf(tmpStr, "%s_T%d", inputFile, i);
	  mfd = openMeasFile(inputFile); /* Open measures file (xva file) */

	  if (!FirstTime)
	    {
	      Nm = OxvaH.size / XVA_NUM / sizeof(COORD_TYPE);
	      printf("m: %f\n", m);
	      //xvadt = OxvaH.saveSteps * dt; 
	      Allocate(OxvaH.size, XVA_ALST, NULL);
	      FirstTime = 0;
	    }
	  
	  genFiles(mfd, Nm, i);
	  close(mfd);
	}

    }
  else
    {
      mfd = openMeasFile(inputFile); /* Open measures file (xva file) */
      Nm = OxvaH.size / XVA_NUM / sizeof(COORD_TYPE);
      Allocate(OxvaH.size, XVA_ALST, NULL);
      genFiles(mfd, Nm, 0);
      close(mfd);
    }

}

/* ============================= >>> get_temps <<< ======================= */
void get_temps(void)
{
  int i;
  char ptt[1024];
  char *strN;
  
  strcpy(ptt, tempStr);
  strN = strtok(ptt, ",");
  sscanf(strN, "%lf", &Temp[0]);
 
  for(i=1; (strN = strtok(NULL, ","))!=NULL; i++)
    {
      sscanf(strN, "%lf", &Temp[i]);
    }
  
  if ((strN == NULL) && (i != PTM))
    {
      printf("[ERROR] number of temperatures must be equal to PTM!!!\n");
      exit(-1);
    }
}

/* =========================== >>> main <<< =============== */
int main(int argc, char** argv)
{
  FILE* fi;
  int i;

  defaults();
  
  args(argc, argv);
  
  /* Legge i parametri dal parameter file in ingresso */
  readPars(xvaparsFile);

  get_temps();

  MD2Fra();

  fi = fopen(fileLista, "w");
  
  for (i = 0; i < njob; i++)
    {
      fprintf(fi, "%s\n", filename[i]);
    }
  
  fclose(fi);
  return 0;
}
