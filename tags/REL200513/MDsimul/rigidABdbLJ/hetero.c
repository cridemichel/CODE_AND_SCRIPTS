#include<mdsimul.h>
#include<heteroconstr.h>

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
  if (!strcmp(heteroTraFile, "!"))
    heteroTraFlag = 0;
  if (!strcmp(heteroAngFile, "!"))
    heteroAngFlag = 0;
}

/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  printf("Syntax: heteroutil [-mathematica] [-p <precision>] -f <hetero_parameters_file>\n");
  printf("Invalid argument!\n");
  exit(-1);
}

/* ========================== >>> defaults <<< ============================ */
void defaults(void)
{
  pi = acos(-1);
  strcpy(xvaparsFile, "xva.par");
  DtTra_mean = 0;
  DtAng_mean = 0;
  m0 = 1.0;
  m1 = 1.0;
  d = 0.5;
  T = 1.0;
  printEvery = 10;
  strcpy(heteroTraFile, "!");
  strcpy(heteroAngFile, "!");
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
      else if (!strcmp(argv[i], "-mathematica"))
	{
	  if ( i+1  == argc )   /* see above */
	    {
	      invalArg(); 
	    }
	  mathFlag = 1;
	  
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
  COORD_TYPE  r01x, r01y, r01z, v01x, v01y, v01z, dSq;

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


* ========================= >>> calcOrientVect <<< ======================== */
void calcOrientVect(COORD_TYPE* u01x, COORD_TYPE* u01y, COORD_TYPE *u01z)
{
  /* DESCRIPTION:
     Calculate the normalized orientational vector of all molecules */ 
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

void calcCMi(int i, COORD_TYPE *Rx, COORD_TYPE *Ry, COORD_TYPE *Rz)
{
  *Rx = m0 * rx[0][i] + m1 * rx[1][i];
  *Ry = m0 * ry[0][i] + m1 * ry[1][i];
  *Rz = m0 * rz[0][i] + m1 * rz[1][i];
  *Rx /= Mtot;
  *Ry /= Mtot;
  *Rz /= Mtot;
}

/* ========================== >>> getPhi <<< =========================== */
void getPhii(int i, COORD_TYPE* fix, COORD_TYPE* fiy, COORD_TYPE* fiz)
{
  *fix = Dphix[i];
  *fiy = Dphiy[i];
  *fiz = Dphiz[i];
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
  strcat(fmtStr, "f");
  strcat(fmtStr, " ");
  strcat(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f");
  strcat(fmtStr, " ");
  strcat(fmtStr, " @ ");
  strcat(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f");
  strcat(fmtStr, " $ ");
  strcat(fmtStr, "%d\n");

}

/* =========================== >>> mkFormatMath <<< ===================== */ 
void mkFormatMath(char* fmtStr)
{
  strcpy(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f");
  strcat(fmtStr, " ");
  strcat(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f");
  strcat(fmtStr, " ");
  strcat(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f");
  strcat(fmtStr, " ");
  strcat(fmtStr, "%.");
  strcat(fmtStr, precision);
  strcat(fmtStr, "f\n");

}

/* ========================== >>> saveAcf <<< ============================== */
void save_hetero(char* fileName, COORD_TYPE* sqrtDisp, int* mobileLst, 
		int mobiles, COORD_TYPE Dt)
{
  int i, mobp, greyLvl;
  COORD_TYPE  L, invL, meanSqr, meanSqrt, radius, Rmx, Rmy, Rmz;
  char fmtStr[128];
  FILE* afs;
  //int Acc;
  printf("fileName:%s\n", fileName);

  printf("t0 (chunks): %d\n", t0);
  printf("t0 (reduced units)= %f t1 = %f\n", ((COORD_TYPE) t0)*xvadt, 
	 ((COORD_TYPE) t0 + Dt)*xvadt);

  /* '!' = don't calculate that autocorrelation function */ 
  if (!strcmp(fileName, "!") || !strcmp(fileName, ""))
    {
      return;
    }
  if (mathFlag)
    mkFormatMath(fmtStr);
  else
    mkFormat(fmtStr);
  
  printf("fmtstr---> %s\n",fmtStr);
  meanSqr = 0.0;
  for (i = 0; i < Nm; ++i)
    {
      meanSqr += Sqr(sqrtDisp[i]);
    }
  
  meanSqr /= Nm;
  meanSqrt = sqrt(meanSqr);
  //printf("format: *%s*\n", fmtStr);
  //Acc = 0;

  afs = aopen(fileName);
  if (!mathFlag) 
    {
      fprintf(afs, ".numAt: 1\n");
      //fprintf(afs, ".atomCol: grey%d\n", greyLevel);
      fprintf(afs, ".Vol: %f\n", Vol);
    }
  L = cbrt(Vol);
  invL = 1.0 / L;

  for(i = 0; i < mobiles; ++i)
    {
      
      //Acc += wtdArr[t];
      //printf("wtdtra[%d]:%d\n", t, wtdArr[t]);
      mobp = mobileLst[i]; 
      Rmx = (rCMt0x[mobp] + rCMt1x[mobp]) / 2.0;
      Rmy = (rCMt0y[mobp] + rCMt1y[mobp]) / 2.0;;
      Rmz = (rCMt0z[mobp] + rCMt1z[mobp]) / 2.0; ;
      // Reduce to first box 
      Rmx = Rmx - L * rint(invL * Rmx);
      Rmy = Rmy - L * rint(invL * Rmy);
      Rmz = Rmz - L * rint(invL * Rmz);
      radius = sqrtDisp[mobp] / meanSqrt;
      greyLvl = baseGrey + rint((maxGrey-baseGrey) 
				* ( 1.0 - ((Rmz + L / 2.0) / L) ));
      if (mathFlag)
	fprintf(afs, fmtStr, Rmx, Rmy, Rmz, radius/RADFACT);
      else
	fprintf(afs, fmtStr, Rmx, Rmy, Rmz, radius/RADFACT, greyLvl);
      /* xvadt*t is the time in the simulation untis*/
    }
  aclose(afs);
  //printf("Acc: %d\n", Acc);
}

/* ======================= >>> getMeanPosition <<< ==========================*/
void getMeanPosition(int mfd, int Nm, int t_cur, 
		     COORD_TYPE* Rx, COORD_TYPE* Ry, COORD_TYPE* Rz)
{
  /* Actually AVERAGE OVER DtTra_mean*2+1 POSITIONS around t_cur
     Following Muranaka Hiwatari 'Spatial heterogeneity in supercooled liquids
     and glasses' */
  
  int i, nn;
  COORD_TYPE Rxtmp, Rytmp, Rztmp;
  int NMP = 2*DtTra_mean + 1;
  /* Set to zero CoM positions arrays */ 
  loop(i, 1, Nm)
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
      loop(i, 1, Nm)
	{
	  calcCMi(i, &Rxtmp, &Rytmp, &Rztmp);
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

/* ======================= >>> getMeanPosition <<< ==========================*/
void getMeanAngPos(int mfd, int Nm, int t_cur, 
		     COORD_TYPE* phix, COORD_TYPE* phiy, COORD_TYPE* phiz)
{
  /* Actually AVERAGE OVER DtTra_mean*2+1 POSITIONS around t_cur
     Following Muranaka Hiwatari 'Spatial heterogeneity in supercooled liquids
     and glasses' */
  
  int i, nn;
  COORD_TYPE phixtmp, phiytmp, phiztmp;
  int NMP = 2*DtAng_mean + 1;
  /* Set to zero CoM positions arrays */ 
  loop(i, 1, Nm)
    {
      phix[i] = 0.0;
      phiy[i] = 0.0;
      phiz[i] = 0.0;
    }

  for (nn = t_cur - DtAng_mean; nn <= t_cur + DtAng_mean; ++nn)
    {
      //printf("t_cur: %d nn: %d\n", t_cur, nn);
      /* if 0 => end of file */
      if (readMeas(mfd, nn, XVA_LIST, NULL) == 0)
	{
	  printf("ERROR: Unexpected end of file\n");
	  printf("nn:%d t_cur: %d\n", nn, t_cur);
	  exit(-1);
	}
      loop(i, 1, Nm)
	{
	  getPhii(i, &phixtmp, &phiytmp, &phiztmp);
	  //if (nn == t_cur && i == 50)
	  //printf("[not avged](%f,%f,%f)\n", Rxtmp, Rytmp, Rztmp);
	  
	  phix[i] += phixtmp;
	  phiy[i] += phiytmp;
	  phiz[i] += phiztmp;
	}
    }
  
  loop(i, 1, Nm)
    {
      phix[i] /= NMP;
      phiy[i] /= NMP;
      phiz[i] /= NMP;
    }
  //printf("[averaged] (%f,%f,%f)\n", Rx[50], Ry[50], Rz[50]);

}

/* ============================= >>> sumSquared <<< ====================== */
COORD_TYPE sumSquared(COORD_TYPE* sqrtDisp, int Nm)
{
  int i;
  COORD_TYPE sum = 0.0;
  for (i = 0; i < Nm; ++i)
    {
      sum += Sqr(sqrtDisp[i]);
    }
  return sum;
}
/* ============================= >>> MAIN <<< ===============================*/
void main(int argc, char** argv)
{
  int mfd, nmbytes;                /* measure file descriptor */
  int i, t1, k = 0;     
  int mobilesTra = 0, mobilesAng = 0;
  /* coounters */
  COORD_TYPE *uxt0, *uyt0, *uzt0, *uxt1, *uyt1, *uzt1, lc = 0.0;
  COORD_TYPE *phit0x, *phit0y, *phit0z, *phit1x, *phit1y, *phit1z;
  COORD_TYPE *DR, *DPhi, sumSqrMob = 0.0;
  COORD_TYPE mobileSqrDisp = 0.0;
  int *mobileTraLst, *mobileAngLst;
  /* NOTE: This variables could be put in the Header of the file */


  defaults();

  /* Take arguments from input line */
  args(argc,argv);
  
  readPars(xvaparsFile);
  //exit(-1);
  printf("Read following values from %s:\n", xvaparsFile);
  printf("m0:%f m1:%f d:%f\n", m0, m1, d);
  printf("Vol:%f T:%f\n", Vol, T);
  printf("Using %s as xva file\n", inputFile);
  /* Set Flags */
  setFlags();
  
  mfd = openMeasFile(inputFile); /* Open measures file (xva file) */

  printf("dt: %f xvasteps: %d\n", dt, OxvaH.saveSteps);

  /* Usefule quantities */
  printf("size: %d\n", OxvaH.size);
  Nm = OxvaH.size / XVA_NUM / sizeof(COORD_TYPE);
  Mtot = m0 + m1;
  xvadt = OxvaH.saveSteps * dt; 
  
  printf("Nm: %d Mtot: %f\n", Nm, Mtot);
  /* Allocate the variables in each measure */
  Allocate(OxvaH.size, XVA_ALST, NULL);

  /* Determine the maximmum t for which calculating the correlation 
     functions */
 
  /* tRun = 0;
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
 
  /* Normalized orientational vectors */
  uxt0 = malloc(nmbytes);
  uyt0 = malloc(nmbytes);
  uzt0 = malloc(nmbytes);
  uxt1 = malloc(nmbytes);
  uyt1 = malloc(nmbytes);
  uzt1 = malloc(nmbytes);

  /* Diffusion coefficents */
  rCMt0x = malloc(nmbytes);
  rCMt0y = malloc(nmbytes);
  rCMt0z = malloc(nmbytes);
  rCMt1x = malloc(nmbytes);
  rCMt1y = malloc(nmbytes);
  rCMt1z = malloc(nmbytes);
  
  /* Angular 'position' */
  phit0x = malloc(Nm*sizeof(COORD_TYPE));
  phit0y = malloc(Nm*sizeof(COORD_TYPE));
  phit0z = malloc(Nm*sizeof(COORD_TYPE));
  phit1x = malloc(Nm*sizeof(COORD_TYPE));
  phit1y = malloc(Nm*sizeof(COORD_TYPE));
  phit1z = malloc(Nm*sizeof(COORD_TYPE));


  /* DR_i */
  DR = malloc(nmbytes);
  DPhi = malloc(nmbytes);

  /* ===============================================================  */

   /* Arrays to store mobile particles */
  mobileTraLst = malloc(nmbytes);
  mobileAngLst = malloc(nmbytes);
  
  /* ========== >>> Correlation functions initialization BEGIN <<< ========= */
  for(i=0; i < Nm; ++i)
    {
      mobileTraLst[i] = 0;
      mobileAngLst[i] = 0;
    }
  
  printf("Begin...\n");
  //printf("t0 up to %d\n", tRun);
  
  //goto savings;
  //printf("tRun: %d\n", tRun);
  
  /* ============ >>> Loop over time <<< ==================== */
  if (heteroTraFlag)
    {
      // printf("inizio inizio  Eta: %.10f\n", Eta[t]);
      
      /* Calc quantities for t0 */
      
      //calcAngVel(ox0, oy0, oz0);
      //calcOrientVect(ux0, uy0, uz0);
      //printf("uno\n");
      getMeanPosition(mfd, Nm, t0, rCMt0x, rCMt0y, rCMt0z); 
      
      t1 = t0 + DtTra;
      //printf("t1:%d t0:%d\n", t1, t0);
     
      /*if (t1 + DtTra_mean >= tRun)
	{
	printf("ERROR: tRun excedeed!\n");
	exit(-1);
	}*/
      
      //printf("t1: %d tRun: %d due\n", t1, tRun);
      getMeanPosition(mfd, Nm, t1, rCMt1x, rCMt1y, rCMt1z);
      
      
      /* ============= >>> SUM OF LOOP OVER PARTICLES <<< ============ */
      for (i = 0; i < Nm; ++i)
	{
	  //printf("DR:%f\n",DR);
	  DR[i] = sqrt(Sqr(rCMt0x[i] - rCMt1x[i]) + 
		       Sqr(rCMt0y[i] - rCMt1y[i]) + 
		       Sqr(rCMt0z[i] - rCMt1z[i]));
	  
	  /* ============= >>> END OF LOOP OVER PARTICLES <<< ============ */
	}
      sumSqrMob = 0.0;
      mobileSqrDisp =  mobTraPerc * sumSquared(DR, Nm);
      for (k = 1; k <= Nm; ++k)
	{
	  lc = sqrt(mobileSqrDisp / ((double) k));
	  mobilesTra = 0;
	  sumSqrMob = 0.0;
	  for (i = 0; i < Nm; ++i)
	    {
	      if (DR[i] > lc)
		{
		  mobileTraLst[mobilesTra] = i;
		  sumSqrMob += Sqr(DR[i]);
		  ++mobilesTra;
		  if (sumSqrMob > mobileSqrDisp)
		    break;
		}
	    }
	  if (sumSqrMob > mobileSqrDisp)
	    break;
	}
      
      if (sumSqrMob < mobileSqrDisp)
	printf("FAILED: There aren't enough mobile particles\n");
      else
	{
	  printf("n. of mobile particles(k=%d): %d\n",k, mobilesTra);
	  printf("mobTraPerc: %f sumSqrMob: %f mobilrSqrDisp: %f\n", 
		 mobTraPerc, sumSqrMob, mobileSqrDisp);
	  printf("lc: %f\n", lc);
	  save_hetero(heteroTraFile, DR, mobileTraLst, mobilesTra, DtTra);
	}
    }
  
  //  close(mfd);
  //  exit(-1);
  /* ============ >>> Loop over time <<< ==================== */
  if (heteroAngFlag)
    {
      getMeanAngPos(mfd, Nm, t0, phit0x, phit0y, phit0z); 
      /* Calc quantities for t0 */
      
      t1 = t0 + DtAng;
      /*
	if (t1 + DtAng_mean >= tRun)
	{
	printf("ERROR: tRun excedeed!\n");
	exit(-1);
	}
      */
      
      getMeanAngPos(mfd, Nm, t1, phit1x, phit1y, phit1z); 
      /* ============= >>> SUM OF LOOP OVER PARTICLES <<< ============ */
      for (i = 0; i < Nm; ++i)
	{
	  /* dphi is in degs not in rads */
	  DPhi[i] = (180.0/pi)*sqrt(Sqr(phit0x[i] - phit1x[i]) + 
				    Sqr(phit0y[i] - phit1y[i]) + 
				    Sqr(phit0z[i] - phit1z[i]));
	}
      sumSqrMob = 0.0;
      mobileSqrDisp =  mobAngPerc * sumSquared(DPhi, Nm);
      for (k = 1; k <= Nm; ++k)
	{
	  lc = sqrt(mobileSqrDisp / ((double) k));
	  mobilesAng = 0;
	  sumSqrMob = 0.0;  
	  for (i = 0; i < Nm; ++i)
	    {
	      if (DPhi[i] > lc)
		{
		  mobileAngLst[mobilesAng] = i;
		  sumSqrMob += Sqr(DPhi[i]);
		  ++mobilesAng;
		  if (sumSqrMob > mobileSqrDisp)
		    break;
		}
	    }
	  if (sumSqrMob > mobileSqrDisp)
	    break;
	}
	  
      //printf("dphi: %f\n", dphi);
      
      /* ============= >>> END OF LOOP OVER PARTICLES <<< ============ */
    
      if (sumSqrMob < mobileSqrDisp)
	printf("FAILED: There aren't enough mobile particles\n");
      else
	{
	  printf("n. of angular mobile particles(k=%d): %d\n",k, mobilesAng);
	  printf("mobAngPerc: %f sumSqrMob: %f mobilrSqrDisp: %f\n", 
		 mobAngPerc, sumSqrMob, mobileSqrDisp);
	  printf("lc: %f\n", lc);
	  save_hetero(heteroAngFile, DPhi, mobileAngLst, mobilesAng, DtAng);
	}
    }

  
  /* ============ >>> End of loop over time origins <<< ==================== */

  /* Normalization 
     ACTUALLY NOT DONE!!!!! */
  
  close(mfd);
}
