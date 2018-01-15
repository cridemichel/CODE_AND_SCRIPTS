#include<mdsimul.h>
#include<constrxva.h>

struct xvaHead OxvaH; /* measure file header (see mdsimul.h)*/

COORD_TYPE XVA_DLST;

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

/*=========================== >>> vectProd <<< =========================== */
void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z)
{
  /* DESCRIPTIOM:
     r3 = [ r1, r2 ] where [ , ] the vectorial product */
  *r3x = r1y * r2z - r1z * r2y; 
  *r3y = r1z * r2x - r1x * r2z;
  *r3z = r1x * r2y - r1y * r2x;
}

/* =========================== >>> SetFlags <<< =============================*/
void setFlags(void)
{
  /* '!' = don't calculate */
  if (!strcmp(velFile, "!"))
    velFlag = 0;
  if (!strcmp(psi1File, "!"))
    psi1Flag = 0;
  if (!strcmp(vhFile, "!"))
    vhFlag = 0;
  if (!strcmp(psi2File, "!"))
    psi2Flag = 0;
  if (!strcmp(C1File, "!"))
    C1Flag = 0;
  if (!strcmp(C2File, "!"))
    C2Flag = 0;
  if (!strcmp(C3File, "!"))
    C3Flag = 0;
  if (!strcmp(C4File, "!"))
    C4Flag = 0;
  if (!strcmp(dphiSqFile, "!"))
    dphiSqFlag = 0;
  if (!strcmp(phiFile, "!"))
    phiFlag = 0;
  if (!strcmp(drSqFile, "!"))
    drSqFlag = 0;
  if (!strcmp(DtFile, "!"))
    DtFlag = 0;
  if (!strcmp(DrFile, "!"))
    DrFlag = 0;
  if (!strcmp(AalphaFile, "!"))
    AalphaFlag = 0;
  if (!strcmp(GselfFile, "!"))
    GselfFlag = 0;
  if (!strcmp(wtdFile, "!"))
    wtdFlag = 0;
  if (!strcmp(GsGsgaussFile, "!"))
    GsGsgaussFlag = 0;
  if (!strcmp(FselfFile, "!"))
    FselfFlag = 0;
}

/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  printf("Syntax: xvautil [-p <precision>] -f <xva_parameters_file>\n");
  printf("Invalid argument!\n");
  exit(-1);
}

/* ========================== >>> defaults <<< ============================ */
void defaults(void)
{
  pi = acos(-1);
  strcpy(xvaparsFile, "xva.par");
  m0 = 1.0;
  m1 = 1.0;
  d = 0.5;
  T = 1.0;
  tCor = 10;
  tgap = 1;
  vhgap = 1;
  tBeg = 0;     /* t iniziale per le  funzioni di correlazione */
  nTeta = 30;
  dteta = 0.05; /* dteta in radiants */
  printEvery = 10;
  Gsnr = 50;
  GsrMax = 1.5;
  kMax = 5.0;
  /*
  strcpy(velFile, "vel.a");
  strcpy(psi1File, "psi1.a");
  strcpy(psi2File, "psi2.a");
  strcpy(vhFile, "vanHove.a");
  strcpy(C1File, "C1.a");
  strcpy(C2File, "C2.a");
  strcpy(C3File, "C3.a");
  strcpy(C4File, "C4.a");
  strcpy(dphiSqFile, "dphiSq.a");
  strcpy(drSqFile, "drSq.a");
  strcpy(DtFile, "Dt.a");
  strcpy(DrFile, "Dr.a");
  strcpy(AalphaFile, "Aa.a");
  strcpy(GselfFile, "Gs.a");
  strcpy(GsGsgaussFile, "GsGsg.a");
  strcpy(FselfFile, "Fs.a");*/
  /* By default don't do anything */
  strcpy(velFile, "!");
  strcpy(psi1File, "!");
  strcpy(psi2File, "!");
  strcpy(vhFile, "!");
  strcpy(C1File, "!");
  strcpy(C2File, "!");
  strcpy(C3File, "!");
  strcpy(C4File, "!");
  strcpy(dphiSqFile, "!");
  strcpy(drSqFile, "!");
  strcpy(ddtdrFile, "!");
  strcpy(ddtdphiFile, "!");
  strcpy(DtFile, "!");
  strcpy(DrFile, "!");
  strcpy(AalphaFile, "!");
  strcpy(GselfFile, "!");
  strcpy(wtdFile, "!");
  strcpy(GsGsgaussFile, "!");
  strcpy(FselfFile, "!");
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
  /* and now take the input file */
  return;

  /* there is another arg ?  if N => ERROR */
  if (++i == argc)  
    {
      invalArg();
    }
  
  /* last arg that is input file name */
  strcpy(inputFile, argv[i]);
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

/* ========================== >>> closeMeasFile <<< ======================= */
void closeMeasFile(int fd)
{
  close(fd);
}

/* ======================= >>> openMeasFileWR <<< ============================ */
int openMeasFileWR(char* FileName, int rank)
{
  int mfd;
  char s[MAX_LENGTH];
  /* open input file for reading and output file for writing */


  sprintf(s, "%s_R%d", FileName, rank);

  if ( (mfd = open(s,O_RDONLY)) == -1 )
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

/* ============================ >>> readMeas <<< =========================== */
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

/* ========================== >>> calcCMVel <<< ============================ */
void calcCMVel(COORD_TYPE *vxCM, COORD_TYPE *vyCM, COORD_TYPE *vzCM)
{
  int i;
  
  for(i = 0; i < Nm; ++i)
    {      
      vxCM[i] = (m0 * vx[0][i] + m1 * vx[1][i]) / Mtot;
      vyCM[i] = (m0 * vy[0][i] + m1 * vy[1][i]) / Mtot;
      vzCM[i] = (m0 * vz[0][i] + m1 * vz[1][i]) / Mtot;
    }
}

/* ========================== >>> calcAngVel <<< ========================== */
void calcAngVel(COORD_TYPE* ox, COORD_TYPE* oy, COORD_TYPE *oz)
{
  int i;
  COORD_TYPE  r01x, r01y, r01z, v01x, v01y, v01z, r01Sq, dSq;

  dSq = Sqr(d);
  
  for(i = 0; i < Nm; ++i)
    {
      
      r01x = rx[0][i] - rx[1][i];
      r01y = ry[0][i] - ry[1][i];
      r01z = rz[0][i] - rz[1][i];
      //printf("rx[%d]: %.10f ry: %.10f rz: %.20f\n", i, rx[0][i], ry[0][i], rz[1][i]);

      //r01x = r01x - L * rint(invL * r01x); /* reduced coordinates */
      //r01y = r01y - L * rint(invL * r01y);
      //r01z = r01z - L * rint(invL * r01z);
      /* Molecules are indivisible(???) */
      v01x = vx[0][i] - vx[1][i];
      v01y = vy[0][i] - vy[1][i];
      v01z = vz[0][i] - vz[1][i];
	
      vectProd(r01x, r01y, r01z, v01x, v01y, v01z, &ox[i], &oy[i], &oz[i]);
      /*r01Sq = Sqr(r01x) + Sqr(r01y) + Sqr(r01z);
      ox[i] = ox[i] / r01Sq;
      oy[i] = oy[i] / r01Sq;
      oz[i] = oz[i] / r01Sq;*/
      ox[i] = ox[i] / dSq;
      oy[i] = oy[i] / dSq;
      oz[i] = oz[i] / dSq;
    }
  //exit(-1);
}

/* ========================= >>> calcOrientVect <<< ======================== */
void calcOrientVect(COORD_TYPE* u01x, COORD_TYPE* u01y, COORD_TYPE *u01z)
{
  /* DESCRIPTION:
     Calculate the normalized orientational vector of all molecules */ 
  COORD_TYPE r01; /* Bond length */
  int i;
  
  for (i = 0; i < Nm; ++i)
    {
      u01x[i] = rx[0][i] - rx[1][i];
      u01y[i] = ry[0][i] - ry[1][i];
      u01z[i] = rz[0][i] - rz[1][i]; 
      /*r01 = Sqr(u01x[i]) + Sqr(u01y[i]) + 
	Sqr(u01z[i]);
      r01 = sqrt(r01);
      u01x[i] /= r01;
      u01y[i] /= r01;
      u01z[i] /= r01;*/
      //printf("r01: %f\n", r01);
      u01x[i] /= d;
      u01y[i] /= d;
      u01z[i] /= d;
    }
}
/* ============================ >>> calcCM <<< ============================ */
void calcCM(COORD_TYPE* rCMx, COORD_TYPE* rCMy, COORD_TYPE* rCMz)
{
  /* DESCRIPTION:
     Calculate the center of mass positions and store them in the array
     passed as arguments */
  int i;
  for(i = 0; i < Nm; ++i)
    {
      rCMx[i] = m0 * rx[0][i] + m1 * rx[1][i];
      rCMy[i] = m0 * ry[0][i] + m1 * ry[1][i];
      rCMz[i] = m0 * rz[0][i] + m1 * rz[1][i];
      rCMx[i] /= Mtot;
      rCMy[i] /= Mtot;
      rCMz[i] /= Mtot;
    }
}

/* ========================== >>> getPhi <<< =========================== */
void getPhi(COORD_TYPE* fix, COORD_TYPE* fiy, COORD_TYPE* fiz)
{
  int i;
  for(i = 0; i < Nm; ++i)
    {
      fix[i] = Dphix[i];
      fiy[i] = Dphiy[i];
      fiz[i] = Dphiz[i];
    }
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

/* ========================== >>> redTime <<< ============================ */
double redTime(int p, int NN, double base)
{


}

/* ========================== >>> saveAcf <<< ============================== */
void saveAcf(char* fileName, COORD_TYPE* acf)
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
  for(t = tBeg; t < tCor; ++t)
    {
      fprintf(afs, fmtStr, ((COORD_TYPE)t)*xvadt, acf[t]);
      /* xvadt*t is the time in the simulation untis*/
    }
  aclose(afs);
}
/* ========================== >>> saveVH <<< ============================== */
void saveVH(char *fileName, COORD_TYPE** vh)
{
  FILE* afs;
  int i, t;
  COORD_TYPE teta;
  char fmtStr[128];

  /* '!' = don't calculate that autocorrelation function */ 
  if (!strcmp(fileName, "!") || !strcmp(fileName, ""))
    {
      return;
    }

  mkFormat(fmtStr);

  afs = aopen(fileName);
  /* save every vhgap steps */
  for(t = tBeg; t < tCor; t += vhgap )
    {
      if (t == 0) continue;
			
      for(i = 0; i < nTeta; ++i) /* Loop over angles */
	{ 
	  teta = ((COORD_TYPE)i + 0.5) * (180.0 / nTeta);/* In degree !!! */ 
          //fprintf(afs, fmtStr, ((COORD_TYPE)(t+1))*xvadt, acf[t]);
	  /* xvadt*t is the time in the simulation untis */
	  fprintf(afs, fmtStr, teta, vh[i][t]);
	}
      fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
    }
  fclose(afs);
}

/* =========================== >>> saveGself <<< ========================= */
void saveGself(char* fileName, COORD_TYPE** gs)
{
  FILE* afs;
  int j, t;
  COORD_TYPE r;
  char fmtStr[128];

  /* '!' = don't calculate that autocorrelation function */ 
  if (!strcmp(fileName, "!") || !strcmp(fileName, ""))
    {
      return;
    }

  mkFormat(fmtStr);

  afs = aopen(fileName);
  /* save every vhgap steps */
  for(t = tBeg; t < tCor; t += vhgap )
    {
      if (t == 0) continue;
			
      for(j = 0; j < Gsnr; ++j) /* Loop over angles */
	{ 
	  r = ((COORD_TYPE)j + 0.5) * (GsrMax /  Gsnr);/* In degree !!! */ 
	  fprintf(afs, fmtStr, r, gs[j][t]);
	}
      fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
    }
  fclose(afs);

}

/* ======================= >>> saveFself <<< =============================== */
void saveFself(char* fileName, COORD_TYPE* fs)
{
  FILE* afs;
  int j, t;
  char fmtStr[128];

  /* '!' = don't calculate that autocorrelation function */ 
  if (!strcmp(fileName, "!") || !strcmp(fileName, ""))
    {
      return;
    }

  mkFormat(fmtStr);
  afs = aopen(fileName);
  for(t = tBeg; t < tCor; t += vhgap )
    {
      if (t == 0) continue;
      fprintf(afs, fmtStr, ((COORD_TYPE)t)*xvadt, fs[t]);
    }
  /* This indicates to xmgr that begins a new set */
  //}
  fclose(afs);
}

/* =========================== >>> get_tt0 <<< ========================= */
int get_tt0(int t0, int t)
{ 
  if (OxvaH.mode = 0)
    {
      /* Salvataggi lineari */
      return (t0 + t);
    }
  else if (OxvaH.mode = 1)
    {
      /* Salvataggi semilogaritmici */
      nplg = OxvaH.NN;
      //printf("tt0: %d\n", tt0);
      if (t > nplg)
	{
	  return (t0 + (t - nplg)*nplg + 1);
	  
	}
      else
	{
	  return (t0 + t);
	}
    }
}

/* ============================= >>> MAIN <<< ===============================*/
void main(int argc, char** argv)
{
  int mfd, acfbytes, nmbytes;                /* measure file descriptor */
  int i, j, t, t0, tt0, tt0Max, tRun, rInt;     
  int tetai, nplg;
  FILE *phifs;
  COORD_TYPE Normv, o0, ott0, doto, teta, r, deltar, Gsg, prefact, k;
  /* coounters */
  COORD_TYPE *vx0, *vy0, *vz0, *ox0, *oy0, *oz0,
    *vxtt0, *vytt0, *vztt0, *oxtt0, *oytt0, *oztt0;
  COORD_TYPE *ux0, *uy0, *uz0, *uxtt0, *uytt0, *uztt0, *rCM0x, *rCM0y,
    *rCM0z, *rCMtt0x, *rCMtt0y, *rCMtt0z;
  COORD_TYPE dotu, drx, dry, drz, dphix, dphiz, dphiy, dr2i;
  /* NOTE: This variables could be put in the Header of the file */

  defaults();

  /* Take arguments from input line */
  args(argc,argv);

  readPars(xvaparsFile);
  //exit(-1);
  printf("Read following values from %s:\n", xvaparsFile);
  printf("m0:%f m1:%f d:%f\n", m0, m1, d);
  printf("tCor:%d Vol:%f T:%f\n", tCor, Vol, T);
  printf("phiFile: %s\n", phiFile);
  printf("vhgap: %d tBeg:%d\n", vhgap, tBeg);
  printf("Using %s as xva file\n", inputFile);
  /* Set Flags */
  setFlags();
  
#ifdef MPI
  /* Usa il file del processo di rank 0 tanto tutti gli header sono uguali */
  mfd = openMeasFile(inputFile, 0);
#else
  mfd = openMeasFile(inputFile); /* Open measures file (xva file) */
#endif
  printf("dt: %f xvasteps: %d\n", dt, OxvaH.saveSteps);

  /* Usefule quantities */
  printf("size: %d\n", OxvaH.size);
  Nm = OxvaH.size / XVA_NUM / sizeof(COORD_TYPE);
  Mtot = m0 + m1;
  /* 14/09/2000: ora dt, xvaSaveMode, NN, base sonoo tutte state messe
     nell'Header del file xva */
  xvaSaveMode = OxvaH.mode;
  dt = OxvaH.dt;
  NN = OxvaH.NN;
  base = OxvaH.base;
  xvadt = OxvaH.saveSteps * dt; 

  printf("Nm: %d Mtot: %f\n", Nm, Mtot);
  /* Allocate the variables in each measure */
  Allocate(OxvaH.size, XVA_ALST, NULL);

  /* Determine the maximmum t for which calculating the correlation 
     functions */
  
  tRun = 0;
  /*
  while(readMeas(mfd, tRun, XVA_LIST, NULL) > 0)
  {
  ++tRun;
  }
  */
  //printf("dopo scan xva file\n");
 
  /* Now dt is the time in the simulation units between two savings on
     xva file */
				
  /* ========================== >>> ALLOCATION <<<< ======================== */
  
  /* Allocation of variable needed by various correletion funciotn 
     calculations */
  nmbytes = sizeof(COORD_TYPE)*Nm; 
  /*Nm elements, each COORD_TYPE bytes long*/  
 
  acfbytes = sizeof(COORD_TYPE)*tCor; 
  /* tCor elements, each COORD_TYPE bytes long */

  /* ========= >>> VELOCITY <<< ========== */
  vx0 = malloc(nmbytes);
  vy0 = malloc(nmbytes);
  vz0 = malloc(nmbytes);
  vxtt0 = malloc(nmbytes);
  vytt0 = malloc(nmbytes);
  vztt0 = malloc(nmbytes);

  velacf = malloc(acfbytes);

  /* ========== >>> ANGULAR VELOCITY <<< ============ */
  ox0 = malloc(nmbytes);
  oy0 = malloc(nmbytes);
  oz0 = malloc(nmbytes);
  oxtt0 = malloc(nmbytes);
  oytt0 = malloc(nmbytes);
  oztt0 = malloc(nmbytes);

  psi1acf = malloc(acfbytes);
  psi2acf = malloc(acfbytes);

  vanHove = (COORD_TYPE**) malloc(nTeta*sizeof(COORD_TYPE**));
  for (i = 0; i < nTeta; ++i)
    vanHove[i] = (COORD_TYPE*) malloc(acfbytes);

  wtd = (COORD_TYPE*) malloc(acfbytes); /* Waiting time distribution */
  Gself = (COORD_TYPE**) malloc(Gsnr*sizeof(COORD_TYPE*));
  GsGsgauss = (COORD_TYPE**) malloc(Gsnr*sizeof(COORD_TYPE*));
  Fself = (COORD_TYPE*) malloc(acfbytes);
  for (j = 0; j < Gsnr; ++j)
    {
      Gself[j] = (COORD_TYPE*) malloc(acfbytes);
      GsGsgauss[j] = (COORD_TYPE*) malloc(acfbytes);
  
    }

  /* The van Hove function depends upon teta and t, that is G = G(teta, t) */

  /* Normalized orientational vectors */
  ux0 = malloc(nmbytes);
  uy0 = malloc(nmbytes);
  uz0 = malloc(nmbytes);
  uxtt0 = malloc(nmbytes);
  uytt0 = malloc(nmbytes);
  uztt0 = malloc(nmbytes);

  /* Legendre polinomials of directions */
  C1acf = malloc(acfbytes);
  C2acf = malloc(acfbytes);
  C3acf = malloc(acfbytes);
  C4acf = malloc(acfbytes);

  /* Diffusion coefficents */
  rCM0x = malloc(nmbytes);
  rCM0y = malloc(nmbytes);
  rCM0z = malloc(nmbytes);
  rCMtt0x = malloc(nmbytes);
  rCMtt0y = malloc(nmbytes);
  rCMtt0z = malloc(nmbytes);
  
  /* Angular 'position' */
  phix = malloc(Nm*sizeof(COORD_TYPE));
  phiy = malloc(Nm*sizeof(COORD_TYPE));
  phiz = malloc(Nm*sizeof(COORD_TYPE));
  phi0x = malloc(Nm*sizeof(COORD_TYPE));
  phi0y = malloc(Nm*sizeof(COORD_TYPE));
  phi0z = malloc(Nm*sizeof(COORD_TYPE));
  phitt0x = malloc(Nm*sizeof(COORD_TYPE));
  phitt0y = malloc(Nm*sizeof(COORD_TYPE));
  phitt0z = malloc(Nm*sizeof(COORD_TYPE));

  /* Diffusion ceofficents */
  Dr = malloc(acfbytes);
  Dt = malloc(acfbytes);
  ddtdrSq = malloc(acfbytes);
  ddtdphiSq = malloc(acfbytes);

  /* Squared displacements */
  drSq = malloc(acfbytes);
  dphiSq = malloc(acfbytes);
  Aalpha = malloc(acfbytes);
  dr4 = malloc(acfbytes);

  /* Normalizations */
  norm = malloc(sizeof(int)*tCor);
  normv = malloc(sizeof(int)*tCor);

  /* ===============================================================  */

  /* ========== >>> Correlation functions initialization BEGIN <<< ========= */
  for(t=0; t < tCor; ++t)
    {
      velacf[t] = 0.0;   /* VELOCITY */
      psi1acf[t] = 0.0;/* ANGULAR VELOCITY */
      psi2acf[t] = 0.0;
  
      for (i = 0; i < nTeta; ++i) vanHove[i][t] = 0.0;
      for (j = 0; j < Gsnr; ++j) 
	{
	  Gself[j][t] = 0.0;
	}
      C1acf[t] = 0.0;
      C2acf[t] = 0.0;
      C3acf[t] = 0.0;
      C4acf[t] = 0.0;
      Dt[t]    = 0.0;
      Dr[t]    = 0.0;
      ddtdrSq[t] = 0.0;
      ddtdphiSq[t] = 0.0;
      drSq[t]  = 0.0;
      dr4[t] = 0.0;
      dphiSq[t]= 0.0;
      normv[t] = 0;
      norm[t] = 0;
    }
  
  printf("Begin...\n");
  printf("%d points for each autocorrelation function.\n", tCor);
  printf("t0 up to %d\n", tRun);
  
  //goto savings;
  //printf("tRun: %d\n", tRun);

  // tgap se xvaSaveMode = 1 sarebbe logBlock 
  if (OxvaH.mode = 1) 
    tgap = OxvH.NN + 1;

  /* 13/09/2000 NOTA: mode = 2 non ancora implementato!!!*/
  if (tgap <= 0) 
    {
      printf("ERROR: tgap should be > 0!\n");
      exit(-1);
    }
  
#ifdef MPI
  /* Legge tutti i files xva scritti dai vari processi MPI */
  closeMeasFile(mfd);
  for(rank = 0; rank < numOfProcs; rank++)
  {
    /*Apre l'xva file prodotto dal processo di rango "rank" */
    openMeasFileWR(inputFile, rank);
#endif 
  /* ============ >>> Loop over time origins <<< ==================== */
  for (t0 = 0;; t0 = t0 + tgap) 
    {
      // printf("inizio inizio  Eta: %.10f\n", Eta[t]);
      
      if ( ((t0 / tgap) % printEvery) == 0) 
	printf("Current t0: %d (up to %d)\n", t0, tRun);
      
      /* if 0 => end of file */
      if (readMeas(mfd, t0, XVA_LIST, NULL) == 0)
	{
	  break;
	  //printf("ERROR: Unexpected end of file\n");
	  //exit(-1);
	}
      
      /* Calculate the value of correlation variables at the actual time
	 orgin */
      /* =============== >>> VELOCITY <<< =============== */
      if (velFlag) 
	calcCMVel(vx0, vy0, vz0);
      
      /* ============= >>> ANGULAR VELOCITY <<< =============== */
      if (psi1Flag || psi2Flag) calcAngVel(ox0, oy0, oz0);
      
      /* ========== >>> NORM. ORIENTATIONAL VECTORS <<< ========== */
      if (C1Flag || C2Flag || C3Flag || C4Flag || vhFlag)
	calcOrientVect(ux0, uy0, uz0);

      /* ========== >>> Positions at time t0 <<< ================== */
      if (DtFlag || drSqFlag || AalphaFlag || GselfFlag || GsGsgaussFlag||
             FselfFlag) 
	calcCM(rCM0x, rCM0y, rCM0z); 

      if( DrFlag || dphiSqFlag)
	getPhi(phi0x, phi0y, phi0z);
	
      //printf("tt0 up to %d\n", tt0Max);
      for( t = 0; t < tCor; ++t)
	{
	
	  /* restituisce il blocco dati t e t0 */
	  tt0 = get_tt0(t0, t);

	  
	  //if (t < tBeg) continue;

	  //printf("-->%f\n", xvadt * t);
	  /* Reads the measures tt0 */
	  
	  /* All the values are put in the XVA_LIST list of variables */ 
	  if (readMeas(mfd, tt0, XVA_LIST, NULL) == 0) 
	    {
	      break;
	      // printf("ERROR: Unexpected end of file\n");
	      //exit(-1);
	    }
	  
	  if (velFlag) calcCMVel(vxtt0, vytt0, vztt0);
	  if (psi1Flag || psi2Flag) calcAngVel(oxtt0, oytt0, oztt0);
	  if (C1Flag || C2Flag || C3Flag || C4Flag || vhFlag)
	    calcOrientVect(uxtt0, uytt0, uztt0);
	  if (drSqFlag || DtFlag || AalphaFlag || GselfFlag || GsGsgaussFlag||
               FselfFlag || wtdFlag) 
	    calcCM(rCMtt0x, rCMtt0y, rCMtt0z); 
	  if( DrFlag || dphiSqFlag )
	    getPhi(phitt0x, phitt0y, phitt0z);
	    
	  
	  /* ============= >>> SUM OF LOOP OVER PARTICLES <<< ============ */
	  for (i = 0; i < Nm; ++i)
	    {
	      /* Velocity */
	      if (velFlag) 
		velacf[t] += vx0[i] * vxtt0[i] + vy0[i] * vytt0[i] +
		  vz0[i] * vztt0[i];
	      
	      /* Angular Velocity correlation functions */
	      if (psi1Flag || psi2Flag)
		{
		  /* Norms of the angular velocities at time t0 and tt0 */ 
		  o0 = sqrt(Sqr(ox0[i]) + Sqr(oy0[i]) + Sqr(oz0[i]));
		  ott0 = sqrt(Sqr(oxtt0[i]) + Sqr(oytt0[i]) + Sqr(oztt0[i]));
		  doto = ox0[i] * oxtt0[i] + oy0[i] * oytt0[i] +
		    oz0[i] * oztt0[i];
		  doto /= o0 * ott0; 
		}

	      if (psi1Flag)
		{
		  psi1acf[t] += doto;
		}
	      if (psi2Flag)
		{
		  psi2acf[t] += 0.5 * ( 3.0 * Sqr(doto) - 1.0);
		}
	
	      if (C1Flag || C2Flag || C3Flag || C4Flag || vhFlag)
		dotu = ux0[i] * uxtt0[i] + uy0[i] * uytt0[i] +
		  uz0[i] * uztt0[i];
	   
	      if (vhFlag)
		{
		  teta = acos(dotu);
		  //printf("teta: %f vhgap %d\n", teta, vhgap);
		  if ( (teta > dteta) && (teta < (pi - dteta)) )
		      {
			/* nTeta is the number of angles for which we 
 			   calculate the van Hove function G = G(t, teta) */
			tetai = (int) (nTeta * acos(dotu) / pi);
			vanHove[tetai][t] += 1.0; 
		      }
		}
	      if (wtdFlag || GselfFlag || GsGsgaussFlag||FselfFlag)
		{
		  deltar = sqrt(Sqr(rCMtt0x[i] - rCM0x[i]) + 
				Sqr(rCMtt0y[i] - rCM0y[i]) +
				Sqr(rCMtt0z[i] - rCM0z[i]));
		  j = (int) (Gsnr * deltar / GsrMax); 
		  if (j < Gsnr) Gself[j][t] += 1.0;
		  
		}
	      if (C1Flag) C1acf[t] += dotu;
	      if (C2Flag) C2acf[t] += 0.5 * (3.0 * Sqr(dotu) - 1.0);
	      if (C3Flag) C3acf[t] += 0.5 * (5.0 * Sqr(dotu) * dotu - 
					     3.0 * dotu);
	      if (C4Flag) C4acf[t] += 0.125 * (35.0 * Sqr(dotu) * Sqr(dotu) - 
					       30.0 * Sqr(dotu) + 3.0);
	      /* Squared displacements amd Self-diffusion coefficents 
		 (Dr & Dt)*/
	      if (drSqFlag || DtFlag || Aalpha || GsGsgaussFlag)
		{
		  drx = rCMtt0x[i] - rCM0x[i];
		  dry = rCMtt0y[i] - rCM0y[i];
		  drz = rCMtt0z[i] - rCM0z[i];
		  dr2i = Sqr(drx) + Sqr(dry) + Sqr(drz);
		  drSq[t] += dr2i;
		  dr4[t] += Sqr(dr2i);
		}
	     
	      if (dphiSqFlag || DrFlag)
		{
		  dphix = phi0x[i] - phitt0x[i];
		  dphiy = phi0y[i] - phitt0y[i];
		  dphiz = phi0z[i] - phitt0z[i];
		  dphiSq[t] += Sqr(dphix) + Sqr(dphiy) + Sqr(dphiz); 
		}
	      ++normv[t];
	    }
	  /* ============= >>> END OF LOOP OVER PARTICLES <<< ============ */
	  ++norm[t]; 
	  /* ===================================== ============ */
	}
    }

  /* Legge tutti i files xva scritti dai vari processi MPI */
  closeMeasFile(mfd);

#ifdef MPI
  }
#endif 
  /* ============ >>> End of loop over time origins <<< ==================== */

  /* Normalization */
  for(t = tBeg; t < tCor; ++t)
    {
      Normv = (COORD_TYPE) (normv[t]);
      if (Normv == 0.0) Normv += 0.000000001;
      velacf[t] /= Normv;
      psi1acf[t] /= Normv;
      psi2acf[t] /= Normv;
      for (i = 0; i < nTeta; ++i) 
	{
	  /* The interval between 0 and pi is divided in nTeta subintervals,
	     so the for each interval the corrisponding teta is the point
	     in the middle */
	  teta = ((COORD_TYPE)i + 0.5) * (pi / nTeta); 
	  /* 0 < teta < pi */  
	  vanHove[i][t] /= Normv * sin(teta);
	  
	  /* In this way when t -> +oo G(teta, t) = 1 for all teta */ 
	  vanHove[i][t] *= 2 * nTeta / pi;
	}
      for (j = 0; j < Gsnr; ++j)
	{
	  r = ((COORD_TYPE) j + 0.5) * (GsrMax / Gsnr);
	  Gself[j][t] *= Gsnr/ GsrMax;/* Gs is a density of probability */
	  Gself[j][t] /= Normv * Sqr(r) * 4.0 * pi;
	}
      
      C1acf[t] /= Normv;
      C2acf[t] /= Normv;
      C3acf[t] /= Normv;
      C4acf[t] /= Normv;
      drSq[t]  /= Normv;
      dr4[t] /= Normv;
      dphiSq[t]/= Normv;
      if (norm[t] == 0.0) norm[t] += 0.000000001;
    }

  /* ======================== >>> SAVINGS <<< ====================  
     Opens various streams for saving correlation functions */
 //savings:
  /* Velocity */ 
  saveAcf(velFile, velacf);
  
  /* Angular velocuty */
  saveAcf(psi1File, psi1acf);
  saveAcf(psi2File, psi2acf);

  /* Van Hove Function */
  saveVH(vhFile, vanHove);
  
  /* Legendre polinomials of orientations */
  saveAcf(C1File, C1acf);
  saveAcf(C2File, C2acf);
  saveAcf(C3File, C3acf);
  saveAcf(C4File, C4acf);
  
  /* Squared displacements */
  saveAcf(dphiSqFile, dphiSq);
  saveAcf(drSqFile, drSq);
  
  /* Calculate the waiting time distribution 'wtd' */
  for (t = tBeg; t < tCor; ++t)
    {
      wtd[t] = 0.0;
      for (i = 0; i < wtdRmax * Gsnr / GsrMax; ++i)/* <--------!!!!!!!!! */
	{
	  wtd[t] += Gself[i][t] * Sqr((i + 0.5) * GsrMax / Gsnr);
	}
      wtd[t] *= 4 * pi * GsrMax / Gsnr; /* = dr */
    }
  /* Calculation of diffusion coefficent */
  for (t = tBeg; t < tCor; ++t)
    {
      /* Normalized correctly yet */
      if (t > tBeg)
	{
	  if (ddtdphiFlag) ddtdphiSq[t] = (dphiSq[t] - dphiSq[t-1]) 
			     / xvadt / 4.0;
	  if (ddtdrFlag) ddtdrSq[t] = (drSq[t] - drSq[t-1]) / xvadt / 6.0; 
	}
      
      Dt[t] = drSq[t] / 
	(6.0 * xvadt * ((COORD_TYPE)t));

      if (DrFlag)
	{
	  Dr[t] = dphiSq[t] /
	    (4.0 * xvadt * ((COORD_TYPE)t));
	}
      if (AalphaFlag)
	{
	  Aalpha[t] = 3.0 * dr4[t] / 5.0 / Sqr(drSq[t]) - 1.0; 
	}
      
      if (GsGsgaussFlag)
	{
	  prefact = pow(3.0 / drSq[t] / pi / 2.0, 1.5);
	  for (i = 0; i < Gsnr; ++i)
	    {
	      r = ((COORD_TYPE) i + 0.5) * (GsrMax / Gsnr); 
	      Gsg = prefact * exp(-3.0*Sqr(r)/(2.0*drSq[t]));
	      /* Difference from a gaussian */
	      GsGsgauss[i][t] = (Gself[i][t] - Gsg) / Gsg;
	    }	
	}
      if (FselfFlag)
	{
	  Fself[t] = 0.0;
	  k = kMax;
	  for(rInt = 0; rInt < Gsnr; ++rInt)
	    {
	      r = ( ((COORD_TYPE) rInt) + 
		    0.5) * (GsrMax / ((COORD_TYPE)Gsnr));
	      Fself[t] += 4.0 * pi * sin(k*r) * Gself[rInt][t] 
		    * r / k;
	    }
	  Fself[t] *= GsrMax / Gsnr; /* multiply for dr */ 
	}
    }
  //printf("1:%s 2:%s\n", DtFile, DrFile);
  /* Diffusion coefficents savings */
  saveAcf(DtFile, Dt);
  saveAcf(DrFile, Dr);
  saveAcf(ddtdrFile, ddtdrSq);
  saveAcf(ddtdphiFile, ddtdphiSq);
  saveAcf(AalphaFile, Aalpha);
  saveAcf(wtdFile, wtd);
  saveGself(GselfFile, Gself);
  saveGself(GsGsgaussFile, GsGsgauss);
  saveFself(FselfFile, Fself);
 
}


