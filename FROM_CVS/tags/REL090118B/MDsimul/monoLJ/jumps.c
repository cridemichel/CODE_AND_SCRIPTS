/* This program follow a paricle in its motion */
#include<mdsimul.h>
#include<jumputil.h>
COORD_TYPE XVA_DLST;

struct xvaHead OxvaH; /* measure file header (see mdsimul.h)*/

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
  for  (i=0; OxvaPar[i].ptr != NULL; ++i) /* parname=NULL menas END */
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
  strcpy(wp,"0,1");
  dt = 0.001;
  d = 0.5;
  m0 = 1.0;
  m1 = 1.0;
  tgap =1;
  tBeg = 0;
  tTot = -1; /* means tTot = tRun (see main)*/
  strcpy(precision,"6");
  pi = acos(-1);
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

/* ============================ >>> calcCM <<< ============================ */
void calcCMi(int ii, COORD_TYPE* rCMx, COORD_TYPE* rCMy, COORD_TYPE* rCMz)
{
  /* DESCRIPTION:
     Calculate the center of mass positions of molecule ii*/
  *rCMx = m0 * rx[0][ii] + m1 * rx[1][ii];
  *rCMy = m0 * ry[0][ii] + m1 * ry[1][ii];
  *rCMz = m0 * rz[0][ii] + m1 * rz[1][ii];
  *rCMx /= Mtot;
  *rCMy /= Mtot;
  *rCMz /= Mtot;
  //printf("rx[0][ii]: %.6f rCMx: %.6f\n", rx[0][ii], *rCMx);
  //printf("Mtot: %f m0: %f m1: %f\n", Mtot, m0, m1);
}
/* ========================= >>> calcOrientVect <<< ======================== */
void calcOrientVecti(int ii, 
		     COORD_TYPE* u01x, COORD_TYPE* u01y, COORD_TYPE *u01z)
{
  /* DESCRIPTION:
     Calculate the normalized orientational vector of ii-th molecule */ 
//  COORD_TYPE r01; /* Bond length */
  
  *u01x = rx[0][ii] - rx[1][ii];
  *u01y = ry[0][ii] - ry[1][ii];
  *u01z = rz[0][ii] - rz[1][ii]; 
  *u01x /= d;
  *u01y /= d;
  *u01z /= d;
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
void saveFunc(char* fileName, COORD_TYPE* acf, int mShift)
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
  for(t = tBeg + mShift; t < tTot-mShift; t=t+tgap)
    {
      fprintf(afs, fmtStr, ((COORD_TYPE)t)*xvadt, acf[t]);
      /* xvadt*t is the time in the simulation untis*/
    }
  aclose(afs);
}

/* ======================== >>> getMeanPosition <<< ========================*/
void getMeanPositions(int mfd, int Np, int t_cur,  
		      COORD_TYPE* Rx, COORD_TYPE* Ry, COORD_TYPE* Rz)
{
  /* Actually AVERAGE OVER THREE POSITIONS separated by Dt_mean chinks 
     Following Muranaka Hiwatari 'Spatial heterogeneity in supercooled liquids
     and glasses' */
  
  int i, nn;
  COORD_TYPE Rxtmp, Rytmp, Rztmp;
  int NMP = 2*DtTra_mean + 1;
  /* Set to zero CoM positions arrays */ 
  loop(i, 1, Np)
    {
      Rx[i] = 0.0;
      Ry[i] = 0.0;
      Rz[i] = 0.0;
    }

  for (nn = t_cur - DtTra_mean; nn <= t_cur + DtTra_mean; ++nn)
    {
      //printf("t_cur: %d nn: %d\n", t_cur, nn);
      /* if 0 => end of file */
      if (readMeas(mfd, nn, XVA_LIST, NULL) == 0)
	{
	  printf("ERROR: Unexpected end of file\n");
	  printf("nn:%d t_cur: %d\n", nn, t_cur);
	  exit(-1);
	}
      loop(i, 1, Np)
	{
	  calcCMi(np[i], &Rxtmp, &Rytmp, &Rztmp);
	  //if (nn == t_cur && i == 50)
	  //printf("[not avged](%f,%f,%f)\n", Rxtmp, Rytmp, Rztmp);
	  
	  Rx[i] += Rxtmp;
	  Ry[i] += Rytmp;
	  Rz[i] += Rztmp;
	}
    }
  
  loop(i, 1, Nm)
    {
      Rx[i] /= NMP;
      Ry[i] /= NMP;
      Rz[i] /= NMP;
    }
  //printf("[averaged] (%f,%f,%f)\n", Rx[50], Ry[50], Rz[50]);

}

/* ======================== >>> getMeanPosition <<< ========================*/
void getMeanAngPos(int mfd, int Np, int t_cur,  
		   COORD_TYPE* ux, COORD_TYPE* uy, COORD_TYPE* uz)
{
  /* Actually AVERAGE OVER THREE POSITIONS separated by Dt_mean chinks 
     Following Muranaka Hiwatari 'Spatial heterogeneity in supercooled liquids
     and glasses' 
     This is substantially a low-pass filter to reduce noise */
  
  int i, nn;
  COORD_TYPE uxtmp, uytmp, uztmp;
  int NMP = 2*DtTra_mean + 1;
  /* Set to zero CoM positions arrays */ 
  loop(i, 1, Np)
    {
      ux[i] = 0.0;
      uy[i] = 0.0;
      uz[i] = 0.0;
    }

  for (nn = t_cur - DtTra_mean; nn <= t_cur + DtTra_mean; ++nn)
    {
      //printf("t_cur: %d nn: %d\n", t_cur, nn);
      /* if 0 => end of file */
      if (readMeas(mfd, nn, XVA_LIST, NULL) == 0)
	{
	  printf("ERROR: Unexpected end of file\n");
	  printf("nn:%d t_cur: %d\n", nn, t_cur);
	  exit(-1);
	}
      loop(i, 1, Np)
	{
	  calcOrientVecti(np[i], &uxtmp, &uytmp, &uztmp);
	  //if (nn == t_cur && i == 50)
	  //printf("[not avged](%f,%f,%f)\n", Rxtmp, Rytmp, Rztmp);
	  
	  ux[i] += uxtmp;
	  uy[i] += uytmp;
	  uz[i] += uztmp;
	}
    }
  
  loop(i, 1, Nm)
    {
      ux[i] /= NMP;
      uy[i] /= NMP;
      uz[i] /= NMP;
    }
  //printf("[averaged] (%f,%f,%f)\n", Rx[50], Ry[50], Rz[50]);

}

/* ======================== >>> main <<< ===================================*/
void main(int argc, char* argv[])
{
  COORD_TYPE **sqrtdr2, xx0[MP], yy0[MP], 
    zz0[MP], xx[MP], yy[MP], zz[MP];
  COORD_TYPE ux[MP], uy[MP], uz[MP], ux0[MP], uy0[MP], uz0[MP];
  COORD_TYPE dotu, **teta;
  char wpstring[128];
  char fnam[1024];
  char* nps, npss[128];
  int tRun, t, nmbytes, tbytes, jj;
  int mfd;
  int n;/* MP is the maximum number of particles that is 
		possible to follow */
  FILE* ffs;

  defaults();
  
  args(argc, argv);
  
  readPars(xvaparsFile);
  mfd = openMeasFile(inputFile); /* Open measures file (xva file) */
  
  Nm = OxvaH.size / XVA_NUM / sizeof(COORD_TYPE);
  Mtot = m0 + m1;
  printf("Mtot: %f\n", Mtot);
  xvadt = OxvaH.saveSteps * dt; 


  Allocate(OxvaH.size, XVA_ALST, NULL);

  /* np is an int array with all particles to follow */
  nps = strtok(wp,",");
  np[0] = atoi(nps);
  for(jj=1; (nps = strtok(NULL,","))!= NULL; ++jj)
    {
      np[jj] = atoi(nps);
    }
 
  np[jj] = -1; /* -1 = END */
  numParts = jj;


  tRun = 0;
  while(readMeas(mfd, tRun, XVA_LIST, NULL) > 0)
    {
      ++tRun;
    }
  printf("Begin to follow molecules:");
  for(jj = 0; jj < numParts; ++jj) printf(" %d ", np[jj]);
  printf("\n");
  if (tTot == -1)
    tTot = tRun; 
  printf("tTot> %d\n", tTot);
  nmbytes = sizeof(COORD_TYPE)*Nm; 
  tbytes = sizeof(COORD_TYPE)*tRun; 
  
  sqrtdr2 = (COORD_TYPE**) malloc(sizeof(COORD_TYPE*)*numParts);
  teta = (COORD_TYPE**) malloc(sizeof(COORD_TYPE*)*numParts);


  for (jj=0; jj < numParts; ++jj)
    {
      sqrtdr2[jj] = malloc(tbytes);
      teta[jj] =  malloc(tbytes);
    }
  
  /* Find jumps, cluster and whatever you want.
     tBeg is the inital step to find from and tTot is the final one*/
  /* if 0 => end of file */
  
  /*================== $$$ Translational Jumps BEGIN $$$ ==================*/
  /* Filter implemented to eliminate high frequencies for 
     translational jumps */  
  getMeanPositions(mfd, numParts, tBeg+DtTra_mean, xx0, yy0, zz0);

  ffs = fopen("frames.cor", "w");
  fprintf(ffs, ".atomRad: 0.475, 0.5\n");
  for (t = tBeg + DtTra_mean; t < tTot - DtTra_mean; t=t+tgap)
    {
      getMeanPositions(mfd, numParts, t, xx, yy, zz);
      fprintf(ffs, ".newframe\n");
      
      for(jj=0; jj<numParts; ++jj)
	{
	  n = np[jj];
	  fprintf(ffs, "%f %f %f\n", rx[0][n], ry[0][n], rz[0][n]);
	  fprintf(ffs, "%f %f %f\n", rx[1][n], ry[1][n], rz[1][n]);
	  
	  /*sqrtdr2[jj][t] = sqrt( Sqr(xx-xx0[jj]) + Sqr(yy-yy0[jj]) 
	    + Sqr(zz-zz0[jj]));*/
	  sqrtdr2[jj][t] = sqrt(Sqr(xx[jj]-xx0[jj]) + Sqr(yy[jj]-yy0[jj]) 
	    + Sqr(zz[jj]-zz0[jj]));
	  
	}
    }

  /* Save all numParts trajectories */
  for(jj = 0; jj < numParts; ++jj)
    {
      sprintf(wpstring, "_%d", np[jj]);
     
      strcpy(fnam, sqrtdr2File);
      strcat(fnam, wpstring);
      saveFunc(fnam, sqrtdr2[jj], DtTra_mean); 
      /* modulus of displacement */
    }
  /*==================== $$$ Translational Jumps END $$$ ====================*/

  close(mfd);\
  fclose(ffs);\
  exit(-1);
  
 /* =================== $$$ Rotational Jumps BEGIN $$$ ==================== */

  getMeanAngPos(mfd, numParts, tBeg+DtAng_mean, ux0, uy0, uz0);
  
  for (t = tBeg + DtAng_mean; t < tTot - DtAng_mean; t=t+tgap)
    {
      
      getMeanAngPos(mfd, numParts, tBeg+DtAng_mean, ux, uy, uz);

      for(jj=0; jj<numParts; ++jj)
	{
	  dotu = ux0[jj] * ux[jj] + uy0[jj] * uy[jj] +
	    uz0[jj] * uz[jj];	  
	  teta[jj][t] = 180.0*acos(dotu)/pi;
	
	}
    }

  /* Save all numParts trajectories */
  for(jj = 0; jj < numParts; ++jj)
    {
      sprintf(wpstring, "_%d", np[jj]);
     
      strcpy(fnam, tetaFile);
      strcat(fnam, wpstring);
      saveFunc(fnam, teta[jj], DtAng_mean);      
      /* orientations as a functions of time */
    }
  /* =================== $$$ Rotational Jumps END $$$ ====================== */

  close(mfd);
  fclose(ffs);
  printf("Files with orientations and displacements saved\n");
}
