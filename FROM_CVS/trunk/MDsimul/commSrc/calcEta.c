/* ========================= >>> calcEta <<< ================================
   This take as input measurement of the off-diagonal element of pressure 
   tensor and of its time integral to calculate the viscosity as a function
   of time and the autocorrelation function of the pressure tensor.
   Each measure should consist of the three off diagonal components of the
   tensor, that is xy, yz, zx.
*/
#include<mdsimul.h>
#include<calcEta.h>

struct measHead OxvaH; /* measure file header (see mdsimul.h)*/

/* ============================== >>> Parsing <<< ===========================*/
void xvaParsing(char stringA[NAME_LENGTH], char stringB[NAME_LENGTH])
{
  int i;
  char msgStrA[512];
  /* check if 'stra', read from the params file, is one of the strings in the 
     singlePar array, if N => ERROR */ 
  for  (i=0; OxvaPar[i].ptr != NULL; ++i) /* parname=NULL means END */
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

/* =========================== >>> SetFlags <<< =============================*/
void setFlags(void)
{
  if (!strcmp(EtaFile, "!"))
    EtaFlag = 0;
  if (!strcmp(PFile, "!"))
    PFlag = 0;
  if (!strcmp(DQFile, "!"))
    DQFlag = 0;
  if (!strcmp(ddtDQFile, "!"))
    ddtDQFlag = 0;
}

/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  printf("Syntax: md2ascii [-p <precision>] -f <xva_parameters_file>\n");
  printf("Invalid argument!\n");
  exit(-1);
}

/* ========================== >>> defaults <<< ============================ */
void defaults(void)
{
  pi = acos(-1);
  strcpy(xvaparsFile, "xva.par");
  T = 1.0;
  tCor = 10;
  tgap = 1;
  tBeg = 0;     /* t iniziale per le  funzioni di correlazione */
  printEvery = 10;
  strcpy(Pinput, "Ptens-0");
  strcpy(DQinput, "DQtens-0");
  strcpy(DQFile, "DQ.a");
  strcpy(ddtDQFile, "ddtDQ.a");
  strcpy(PFile, "Pacf.a");
  strcpy(EtaFile, "eta.a");
  strcpy(precision,"6");
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
}

/* ======================= >>> openMeasFile <<< ============================ */
int openMeas(char* FileName)
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
  if (read(mfd, &OxvaH, sizeof(struct measHead)) ==  -1)
    {
      perror("read header");
      exit(-1);
    } 
  //printf("steps ----->%d\n", OxvaH.saveSteps);
  return mfd;
}

/* ========================== >>> closeMeasFile <<< ======================= */
void closeMeasFile(int fd)
{
  close(fd);
}

/* ======================= >>> readSMeas <<< ========================= */
int readSMeas(int fd, int nmeas, 
	      COORD_TYPE* Cxy, COORD_TYPE* Cyz, COORD_TYPE* Czx)
{
  /* DESCRIPTION:
     read the nmeas measure with nmeas >= 0,
     Cxy, ... are the off-diagonal element of a tensor at a certain time */
  int begByte, rs1, rs2, rs3; 
  
  begByte = sizeof(struct measHead);
  lseek(fd, begByte + nmeas * OxvaH.size, SEEK_SET);

  rs1 = read(fd, Cxy, sizeof(COORD_TYPE));
  rs2 = read(fd, Cyz, sizeof(COORD_TYPE));
  rs3 = read(fd, Czx, sizeof(COORD_TYPE));

  if ( (rs1 == -1) || (rs2 == -1) || (rs3 == -1))   
    {
      printf("ERROR: Unable to read file\n");
      perror("read");
      exit(-1);
    }
  else if (rs1+rs2+rs3 < OxvaH.size)
    {
      return 0; /* 0 = end of file */
    }

  return 1; /* 1 == OK */
}

/* =========================== readAllMeas <<< ============================ */
void readAllMeas(int fd, COORD_TYPE* Cxy, COORD_TYPE* Cyz, COORD_TYPE* Czx)
{
  int j;
  /* NOTE: Here Cxy, ... are array of tRun element */
  for (j = 0; j < tRun; ++j)
    {
      readSMeas(fd, j, &Cxy[j], &Cyz[j], &Czx[j]);
    }
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

/* ========================== >>> saveAcf <<< ============================== */
void saveAcf(char* fileName, COORD_TYPE* acf, int beg, int end)
{
  int t;
  char fmtStr[128];
  FILE* afs;

  /* '!' = don't calculate that autocorrelation function */ 
  if (!strcmp(fileName, "!") || !strcmp(fileName, ""))
    {
      return;
    }
  
  mkFormat(fmtStr);
  
  //printf("format: *%s*\n", fmtStr);

  afs = aopen(fileName);
  for(t = beg; t < end; ++t)
    {
      fprintf(afs, fmtStr, ((COORD_TYPE)(t+0.5))*xvadt, acf[t]);
      /* xvadt*t is the time in the simulation untis */
    }
  aclose(afs);
}

/* ============================= >>> MAIN <<< ===============================*/
int main(int argc, char** argv)
{
  int mfd, acfbytes;                /* measure file descriptor */
  int t, t0, tt0, tt0Max;     
  COORD_TYPE dummy1, dummy2, dummy3;

  defaults();

  /* Take arguments from input line */
  args(argc,argv);

  readPars(xvaparsFile);

  if (tgap <= 0) 
    {
      printf("ERROR: tgap should be > 0!\n");
      exit(-1);
    }
      
  printf("Read following values from %s:\n", xvaparsFile);
  printf("tCor:%d Vol:%f T:%f\n", tCor, Vol, T);
  printf("tBeg:%d\n", tBeg);
  printf("Using %s as pressure-tensor file\n", Pinput);
  printf("Using %s as DQ-tensor file\n", DQinput);
  
  /* Set Flags */
  setFlags();

  /* Normalizations */
  norm = malloc(sizeof(int)*tCor);
  acfbytes = sizeof(COORD_TYPE)*tCor; 
	
  /* Calculate the viscosity as a function of time */
  if (EtaFlag)
    {
      printf("Calculating the viscosity as a function of time\n");

      mfd = openMeas(DQinput);
      printf("dt: %f xvasteps: %d\n", dt, OxvaH.saveSteps);
      /* Usefule quantities */
      printf("size of each measure: %d\n", OxvaH.size);
      xvadt = OxvaH.saveSteps * dt; 

      /* Determine the maximmum t for which calculating the correlation 
	 functions */
      tRun = 0;
      while(readSMeas(mfd, tRun, &dummy1, &dummy2, &dummy3) > 0) ++tRun;

      /* Now dt is the time in the simulation units between two savings on
	 xva file */
      
      /* ========================= >>> ALLOCATION <<<< ===================== */
      /* tCor elements, each COORD_TYPE bytes long */

      Eta = malloc(acfbytes);
      DQ  = malloc(acfbytes);
      ddtDQ = malloc(acfbytes);
      
      DQxy = malloc(sizeof(COORD_TYPE)*tRun);
      DQyz = malloc(sizeof(COORD_TYPE)*tRun);
      DQzx = malloc(sizeof(COORD_TYPE)*tRun);
      readAllMeas(mfd, DQxy, DQyz, DQzx);
      
      /* ========== >>> Correlation functions initialization BEGIN <<< ===== */
      for(t=0; t < tCor; ++t)
	{
	  Eta[t] = 0.0;
	  norm[t] = 0;
	}

      printf("Begin...\n");
      printf("%d points for each autocorrelation function.\n", tCor);
      printf("t0 up to %d\n", tRun);

      /* ============ >>> Loop over time origins <<< ==================== */
      for (t0 = 0; t0 < tRun ; t0 = t0 + tgap) 
	{
	  
	  if ( ((t0 / tgap) % printEvery) == 0) 
	    printf("Current t0: %d (up to %d)\n", t0, tRun);
      
	  /* tt0Max = min(t0, tCor) */
	  if (tRun < t0 + tCor) 
	    tt0Max = tRun;
	  else 
	    tt0Max = t0 + tCor;
      
	  /* Calculate the value of correlation variables at the actual time
	     orgin */
	  DQxy0 = DQxy[t0];
	  DQyz0 = DQyz[t0];
	  DQzx0 = DQzx[t0];
	  //printf("t0=%d DQ=(%f,%f,%f)\n", t0, DQxy[t0], DQyz[t0], DQzx[t0]);
	  //fprintf(stderr, "%f %f %f %f\n", xvadt*t0, 
	  //  DQxy0/((double)t0)/xvadt, DQyz0/((double)t0)/xvadt, 
	  //  DQzx0/((double)t0)/xvadt);
	  
	  for( tt0 = t0; tt0 < tt0Max; ++tt0)
	    {
	      t = tt0 - t0;
	      if (t < tBeg) continue;
	      
	      DQxytt0 = DQxy[tt0];
	      DQyztt0 = DQyz[tt0];
	      DQzxtt0 = DQzx[tt0];
	      Eta[t] += (Sqr(DQxy0 - DQxytt0) + Sqr(DQyz0 - DQyztt0) +
			 Sqr(DQzx0 - DQzxtt0));
	      ++norm[t]; /* Normalization for viscosity */
	      /*
		if (t == 25000000)
	        {
		printf("accum Eta: %f norm[%d]: %d\n", Eta[t], t, norm[t]);
		printf("sum: %f\n", ( Sqr(DQxy0 - DQxytt0) + 
		Sqr(DQyz0 - DQyztt0) +
		Sqr(DQzx0 - DQzxtt0)) );
		}
	      */
	    }
	  /* ===================================== ============ */
	}
      /* ============ >>> End of loop over time origins <<< ============== */
  
      /* Normalization */
      for(t = 0; t < tCor; ++t)
	{
	  if (norm[t] == 0) 
	    {
	      norm[t] = 1;
	      //printf("norm[%d]: %d\n", t, norm[t]);
	    }
	  printf("norm[%d]: %d\n", t, norm[t]);
	  Eta[t] /= (COORD_TYPE) norm[t];
	  DQ[t] = Eta[t] / 6.0;
	  if (t > tBeg) ddtDQ[t] = (DQ[t] - DQ[t-1]) / xvadt;
	  if (ddtNormFlag) ddtDQ[t] *= Vol / T;
	  Eta[t] *= Vol / T;
	  Eta[t] /= ( 3.0 * xvadt * 2.0 * ((COORD_TYPE) t + 0.5) );
 
	}
    
      /* ==================== >>> SAVING <<< ============================= */
      saveAcf(EtaFile, Eta, tBeg, tCor);
      saveAcf(DQFile, DQ, tBeg + 1, tCor);
      saveAcf(ddtDQFile, ddtDQ, tBeg + 1, tCor);
      free(DQxy);
      free(DQyz);
      free(DQzx);
      free(Eta);
      close(mfd);
    }
 
  /* Making the Pacf autocorrelation function */
  if (PFlag)
    {
      printf("Making pressure-tensor autocorrelation function...\n");
     
      mfd = openMeas(Pinput);
      printf("dt: %f xvasteps: %d\n", dt, OxvaH.saveSteps);
      /* Usefule quantities */
      printf("size of each measure: %d\n", OxvaH.size);
      xvadt = OxvaH.saveSteps * dt; 

      /* Determine the maximmum t for which calculating the correlation 
	 functions */
      tRun = 0;
      while(readSMeas(mfd, tRun, &dummy1, &dummy2, &dummy3) > 0) ++tRun;

      /* Now dt is the time in the simulation units between two savings on
	 xva file */
      
      /* ========================= >>> ALLOCATION <<<< ===================== */
      Pacf = malloc(acfbytes);
      Pxy = malloc(sizeof(COORD_TYPE)*tRun);
      Pyz = malloc(sizeof(COORD_TYPE)*tRun);
      Pzx = malloc(sizeof(COORD_TYPE)*tRun);
      readAllMeas(mfd, Pxy, Pyz, Pzx);
      
      /* ========== >>> Correlation functions initialization BEGIN <<< ===== */
      for(t=0; t < tCor; ++t)
	{
	  Pacf[t] = 0.0;
	  norm[t] = 0;
	}

      printf("Begin...\n");
      printf("%d points for each autocorrelation function.\n", tCor);
      printf("t0 up to %d\n", tRun);
      
      /* ============ >>> Loop over time origins <<< ==================== */
      for (t0 = 0; t0 < tRun ; t0 = t0 + tgap) 
	{
	  
	  if ( ((t0 / tgap) % printEvery) == 0) 
	    printf("Current t0: %d (up to %d)\n", t0, tRun);
      
	  /* tt0Max = min(t0, tCor) */
	  if (tRun < t0 + tCor) 
	    tt0Max = tRun;
	  else 
	    tt0Max = t0 + tCor;
      
	  /* Calculate the value of correlation variables at the actual time
	     orgin */
	  Pxy0 = Pxy[t0];
	  Pyz0 = Pyz[t0];
	  Pzx0 = Pzx[t0];
	  
	  for( tt0 = t0; tt0 < tt0Max; ++tt0)
	    {
	      t = tt0 - t0;
	      if (t < tBeg) continue;
	      Pxytt0 = Pxy[tt0];
	      Pyztt0 = Pyz[tt0];
	      Pzxtt0 = Pzx[tt0];
	      Pacf[t] += Pxy0 * Pxytt0 + Pyz0 * Pyztt0 +
		Pzx0 * Pzxtt0;
	      ++norm[t];
	    }
	}

      /* Normalization */
      for(t = 0; t < tCor; ++t)
	{
	  if (norm[t] == 0) {
	    norm[t] = 1;
	    //printf("norm[%d]: %d\n", t, norm[t]);
	  }
	  Pacf[t] /= 3.0 * ((COORD_TYPE) norm[t]);
	  if (PNormFlag) Pacf[t] *= Vol/T; 
	}

      /* ======================= >>> SAVING <<< ========================== */
      saveAcf(PFile, Pacf, tBeg, tCor);
      
      free(Pxy);
      free(Pyz);
      free(Pzx);
      free(Pacf);
      /* close input file (xva file) */
      close(mfd);
    }
  return 0;
}


