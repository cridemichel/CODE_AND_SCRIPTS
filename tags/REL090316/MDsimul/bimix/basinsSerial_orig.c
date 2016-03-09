/* CONJGATE GRADIENT  
   VERSIONE PER BKS con PETER'S PATHCES
*/
#define BASINS
#define MAIN
#include<mdsimul.h>
/*#define MD_USE_RINT*/

#define ITER_ELAPSED 1
/*=== Alcune di queste sono variabili fittizie definite per poter linkare il 
  file mdarray.c === */
struct simStat OsimStat;
int ENDSIM;
char msgStrA[MSG_LEN];
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE Vc, V, W, K,
  WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, 
  T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz,
  Pxy, Pyz, Pzx, Wm, Wmxx, Wmyy,
  Wmzz, Wmxy, Wmyz, Wmzx, Pmxx, Pmyy,
  Pmzz, Pmxy, Pmyz, Pmzx, T1mxy, 
  Patxy, Patyz, Patzx, Patxx, Patyy, 
  Patzz, T1myz, T1mzx, T1mxx, T1myy,
  T1mzz;

int nperx, npery, nperz;
double anperx, anpery, anperz, invLH;

COORD_TYPE Mtot;
/* used by linked list routines */
int *head, *list, *map;  /* arrays of integer */
int NCell, mapSize, M;

/* neighbour list method variables */
COORD_TYPE dispHi;
int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
char TXTA[10][MSG_LEN];
char TXT[MSG_LEN];

// Potential anti-crystalline energy
double VLJ, WLJ;
int NNC;
int SEGSIZE;        /* length in bytes of an array */

/* string used to send  messages */
char msgStrA[MSG_LEN],msgStrB[MSG_LEN],msgStrC[MSG_LEN]; 
								   
unsigned char BAK, STA, BAKT;
/* ===========================================================================*/

extern void writeAllCor(FILE* );
extern void readAllCor(FILE* );
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);

struct params Oparams;
struct progStatus OprogStatus;

/* ======================== >>> minimumImageFB <<< ======================== */
double minimumImageFB(double c)
{
  if (abs(c) > (L / 2.0))
    {
      if (c > 0)
	c -= L;
      else
	c += L;
    }
  return c;
}

/* ======================= >>> DSIGN <<< ================================ */
double DSIGN(double a, double b)
{
  if (b > 0)
    return fabs(a);
  else
    return -fabs(a);
  /*return (fabs(b)/b) * fabs(a);*/
}

#define DABS fabs

int dietag;

char filewalter[512], filename[1024][256];

int rangsl, numwrk;

const int iwtagHdr = 1, iwtagCor = 2;
int njob; 
double Ftol, Epoten, EpotIni, Emin, L, invL, fnorm;

/*double Vc, V, W, VAC, WAC, VLJ;*/

double *xicomx[NA], *xicomy[NA], *xicomz[NA], *pcomx[NA], 
  *pcomy[NA], *pcomz[NA], *xix[NA], *xiy[NA], 
  *xiz[NA], *Gx[NA], *Gy[NA], *Gz[NA], *Hx[NA], *Hy[NA], *Hz[NA];

/*double  *rx, *ry, *rz, *Fx, *Fy, *Fz, Vol;*/
double pi;

char pwd[512], percorso[512];
char logFile[512];


/* ============ >>> neighbourlist <<< ===================== */
/*COORD_TYPE dispHi;*/
int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;

void initAll(int[NA]);
void Forces(int[NA]);
void cleanstring(char *);
void readRestart(char filewalter[132]);
void frprmn(int Nm[NA], double Ftol, double *Fret);
void writeconf(char filewalter2[132]);
void mnbrak(double *AX, double *BX, double *CX, double *FA, double *FB, double *FC,
	    double (*func)(double));
void writeEnergy(char filewalter[132], double T, int steps,
		 double VpotIni, double VpotEnd);
extern void setToZero(int size, COORD_TYPE *arrx, ...);

extern void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigma[NA][NA]);
extern void LJForce(COORD_TYPE epsab[NA][NA], COORD_TYPE sigab[NA][NA], COORD_TYPE rcut);
extern void readBakAscii(char* fn);
extern void saveBakAscii(char *fn);
extern int** AllocMatI(int size1, int size2);

/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintf(int mode, char *text, ...)
{
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;
  /*  FILE* of;*/
  va_start(ap, text);

  strcpy(tstri, text);
  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      strcat(tstri, sptr);
    }

 
  /*of = fopen(logFile, "a");
    fprintf(of, tstri); 
    fclose(of);*/
  /* e anche su stdout */
  printf(tstri);
 
  va_end(ap);    
    
}
/* ========================== >>> mdMsg <<< ===================================== */
void mdMsg(int mode, int errnum, char *when, char *errType, char* sysCall,
  char* text, ...)
{
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;
  /*FILE* of;*/

  va_start(ap, text);
  strcpy(tstri, text);
  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      strcat(tstri, sptr);
    }
  /*of = fopen(logFile, "a");
    fprintf(of, tstri); 
    fclose(of);
   */
  printf(tstri);
  va_end(ap);
}

char* appSw(char *fileName,const unsigned char which){return NULL;}
char* absTmpTape(const char* fileSrc){return NULL;}
char* absTmpHD(const char* fileSrc){return NULL;}
char* absTmpAsciiHD(const char* fileSrc){return NULL;}
char* absMisHD(const char* fileSrc){return NULL;}
char* absMisTape(const char* fileSrc){return NULL;}
char* absXvaTape(const char* fileSrc){return NULL;}
char* absXvaHD(const char* fileSrc){return NULL;}
int mdWrite(int fdes, char* when, char* errMsg, int mode,
	    int size, void *punta_arr){return -1;}
int  mdRead(int fdes, char* when, char* errMsg, int mode, 
	    int size, void *pointer){return -1;}
int readSegs(int fdes, char* when, char *errMsg, int mode,
	     int size, void* pointer, ...){return -1;}

FILE* fopenMPI(char* fn, char* how)
{
  return fopen(fn, how);
}

int mdCreat(char* fileName, char* when, char* errMsg, int mode){return -1;}
int mdClose(int fd, char* when, char* errMsg, int mode){return -1;}
int mdOpen( char* fileName, char* when, char* errMsg, int mode){return -1;}
void writeSegs(int fdes, char* when, char* errMsg, int mode, int size,
	       void* pointer, ...){}
struct measure Omeasure[]={{NULL,0,NULL}};
unsigned char sw(unsigned char onOff){return 0;}
void scaleVelocities(double fact){}
/* ---------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------*/

/* ============================= >>> main <<< ============================== */
int main(int argc, char** argv)
{
  int NUM_WORK_REQUEST;
  int iwork;
  /*int steps, stepsok;*/
  FILE* finp;
  FILE* pf;
  char listaFile[1024];
  double Fret;
 
  /* L'unico argomento deve essere la lista di file da minimizzare */
  if (argc == 2)
    {
      if ((pf = fopen(argv[1],"r"))==NULL)
	{
	  printf("Unable to open file: %s\n", argv[1]);
	  exit(-1);
	}

      /* printf("arg: %s\n",argv[1]);*/
      strcpy(listaFile, argv[1]);
    }
  else
    {
      printf("Devi fornire come argomenti il nome del file contenente il percorso\n");
      printf("e il file contenente la lista di files!\n");
      exit(-1);
    }
  
  if( (finp = fopen(listaFile, "r"))==NULL )
    {
      printf("Unable to open file: %s\n", listaFile);
      exit(-1);
    }
  
  njob = 1;
  
  /* NOTA:
	 la prima riga deve contenere il percorso dei file da quenchare mentre
	 le altre righe devono contenere i nomi dei file senza il percorso */
  
  fscanf(finp, "%s", percorso);
  printf("Leggo i file da quenchare dalla directory:%s\n", percorso);
  
  while(!feof(finp))
    {
      if (fscanf(finp, "%s", filename[njob]) == 1)
	{
	  njob++;
	}
    }
  
  njob--;
  
  printf(" Number of conf to quench= %d\n", njob);
  
  NUM_WORK_REQUEST = njob;
  
  for(iwork = 1; iwork <= njob; iwork++)
    {

      printf("filename[iwork]: %s\n", filename[iwork]);

      strcpy(filewalter, percorso);
      strcat(filewalter, filename[iwork]);
      printf("filenamefull: %s\n", filewalter);
      readRestart(filewalter); 
      printf("Vol: %f rNebrShell:%f\n", Vol, OprogStatus.rNebrShell);	
      initAll(Oparams.parnum);
      /*
	 sprintf(logFile,"%sQncT%.6G.%dIW%d.log", 
    	 percorso, 
	 Oparams.T, 
	 Oparams.curStep, iwork);
       */
      L = cbrt(Vol);
      invL = 1.0 / L;
      /*INITALCHARGES(); per i SILICA*/
      Forces(Oparams.parnum);

      printf(" EPOTEN INIZIALE =     %.5f | %.5f\n",
	      Epoten, Epoten/((double)Oparams.parnum[0]+Oparams.parnum[1]));

      EpotIni = Epoten;

      /*
	 minimization configuration
       */
      Ftol = 1.0E-8;
      printf("Using tol= %e\n", Ftol);
      frprmn(Oparams.parnum, Ftol, &Fret);
      Forces(Oparams.parnum);
      Emin = Epoten;
      printf("DONE  WORK N.%d\n", numwrk);

      writeconf(filewalter);
      writeEnergy(filewalter, Oparams.T, Oparams.curStep,
		  EpotIni, Emin); /*Questa ruotine va riscritta*/

    }

  printf(" Finito.\n");
  /*
     attende risultati da un processo generico
     until work requests have been exhausted.
   */
  
  return 0;
}
  
/* -----------------------------------------------------------------------------------
   NOTA BENE: pure questo dipende dalla simulazione */


/* ======================== >>> InitAl() <<< =========================== */
void initAll(int Nm[NA])
{
  int size;
  static int firstTime = 1;
  /*
     size = Nm[0] * sizeof(double);
     AllocCoord(size, ALLOC_LISTA, NULL);
     size = Nm[1] * sizeof(double);
     AllocCoord(size, ALLOC_LISTB, NULL);
   */
  
  NNC = 0;
  pi = 2.0*acos(0);
  
  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac * (Nm[0]+Nm[1]) * NA;
  printf("nebrTabMax: %d NA: %d\n", nebrTabMax, NA);
  nebrNow = 1;
  if (firstTime)
    {
      nebrTab = AllocMatI(2, nebrTabMax); 
      firstTime = 0; 
    }
}

/* ------------------------------------------------------------------------------- */
void readOne(char* fn)
{
  FILE* fs; 
  int firstTime = 1;
  
  if ((fs = fopen(fn, "r")) == NULL)
    {
      printf("Problem opening restart file %s ", fn);
    }

  readAsciiPars(fs, opro_ascii);

  readAsciiPars(fs, opar_ascii);

  /* Entrambe queste macro sono definite nel file mono_DPT.h */
  /* read up to coordinates begin */
  if (firstTime)
    {
      firstTime = 0;
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

    }

  readAllCor(fs);
  
  fclose(fs);
}


/* ========================== >>> readStart <<< ============================ */
void readRestart(char filewalter[132])
{
  /*equivalence (nrecord, recname(1));*/
  int i, a;
  
  readOne(filewalter);
	  
  /*
    check of periodic boundary conditions on the read data
  */
  L = cbrt(Vol);
  invL = 1.0/L;
  invLH = 1.0/(L/2.0); 
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  
	  /*  
	      if (abs(anperx)+abs(anpery)+abs(anperz).ne.0) then
	      print *, 'Fuori scatola in restart',i,nperx,npery,nperz
	      end if*/
	  /* Scala al first box */

#ifdef MD_USE_RINT
	  rx[a][i] -= L*rint(invL*rx[a][i]);
  	  ry[a][i] -= L*rint(invL*ry[a][i]);
  	  rz[a][i] -= L*rint(invL*rz[a][i]);
#else
	  nperx =(int) (fabs(rx[a][i])*invLH);                                         
	  anperx=DSIGN(((double)(nperx)*L),rx[a][i]);                                  
	  npery =(int) (fabs(ry[a][i])*invLH);                                          
	  anpery =DSIGN(((double)(npery)*L),ry[a][i]);                                  
	  nperz =(int) (fabs(rz[a][i])*invLH);           
	  anperz=DSIGN(((double)(nperz)*L),rz[a][i]);
	  rx[a][i] -= anperx;
	  ry[a][i] -= anpery;
	  rz[a][i] -= anperz;
#endif
	}
    }
}
/*
    *****************************************************************
    Substitute the first three characters in the file with Qnc
    ******************************************************************  
    */
void writeconf(char filewalter[132])
{
  char fina[132];
  char fina2[512];
  strcpy(fina, filewalter);
  fina[0]='Q';
  fina[1]='n';
  fina[2]='c'; 
  
  strcpy(fina2, percorso);
  strcat(fina2, fina);
  printf("Scritto il file quenchato:%s\n", fina2);
  saveBakAscii(fina2);
}

/* ========================= >>> writeEnergy <<< ========================= */
void writeEnergy(char filewalter[132], double T, int steps,
		 double VpotIni, double VpotEnd)
{
  char fina[132], fina2[512];
  FILE* ndat;

  strcpy(fina, filewalter);
  fina[0]='E';
  fina[1]='n';
  fina[2]='e'; 
  
  strcpy(fina2, percorso);
  strcat(fina2, fina);

  ndat = fopen(fina2, "a");
  fprintf(ndat,"%.5G %d %.15G %.15G\n", T, steps, VpotIni, VpotEnd);
  printf("Scritto il file Ene:%s\n", fina2);
  fclose(ndat);

}

/* ========================= >>> cleanstring <<< ========================== */
void cleanstring(char *filei)
{
  char fileo[132];
  char recname[132];
  int ex[132], kk, ns;

  //printf("qui\n");
  strcpy(fileo,filei);
  /*
    clean from empty spaces
  */
  
  for(kk = 0; kk < 132; kk++)
    {
      ex[kk] = 1; 
      if (recname[kk] == ' ') 
	ex[kk] = 0; 
    }  
  
  ns=0;
  for(kk = 0; kk < 132; kk++)
    {
      if (ex[kk]) 
	{
	  ns = ns+1;
	  recname[ns] = recname[kk];
	}
    }
  
  for(kk = ns+1; kk < 132; kk++)
    { 
      recname[kk]=' ';
    }

  strcpy(filei,fileo);
  printf("fine\n");
}

/* =========================== >>> min <<< =============================== */
double min(double a, double b)
{
  if (a >= b)
    {
      return b;
    }
  else
    {
      return a;
    }
}

/* =========================== >>> max <<< ================================= */
double max(double a, double b)
{
  if (a >= b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


/* ============================ >>> brent <<< ============================ */
double brent(double AX, double BX, double CX, double (*FFF)(double),
	     double tol, double *xmin)
{
  const int ITMAX = 100;
  const double CGOLD = 0.3819660, ZEPS = 1.0E-15;
  double FU, FV, FW, XM, X, U, D, ETEMP, FX, E, A, B, V, tol1, tol2, R, Q, P;
  int iter;

  A = min(AX, CX);
  B = max(AX, CX);

  D = 0.0; /*
	     verificare che questa inizializzazione sia giusta!!!!!!!!!!!!!!!!!!!!!
	     perché nel codice fortran non anche se viene fatta automaticamente 
	   */
  
  V = BX;
  W = V;
  X = V;
  E = 0.0;
  FX = FFF(X);
  FV = FX;
  FW = FX;

  for(iter = 0; iter < ITMAX; iter++)
    {
      XM = 0.5 * (A+B);
      tol1 = tol*fabs(X) + ZEPS;
      tol2 = 2.0*tol1;
        
      if ( fabs(X-XM) <= (tol2 - 0.5 * (B-A)) ) 
	goto tre;
      
      if(fabs(E) > tol1) 
	{
          R = (X-W) * (FX-FV);
          Q = (X-V) * (FX-FW);
          P = (X-V) * Q - (X-W) * R;
          Q = 2.0 * (Q-R);
          if(Q > 0.0) 
	    {
	      P = -P;
	    }
          Q = fabs(Q);
          ETEMP = E;
          E = D;
          if(fabs(P) >= fabs(0.5 * Q * ETEMP)|| 
	     P <= Q * (A-X) || P >= Q * (B-X)) 
	    goto uno;
	  D = P/Q;
	  U = X+D;
          if(U-A < tol2 || B-U < tol2) 
	    D = DSIGN(tol1,XM-X);
	 
	  goto due;
	}
    uno: 
      if(X >= XM)
	E = A-X;
      else
	E = B-X;
      D = CGOLD * E;
    due:
      if(fabs(D) >= tol1)
	U = X + D;
      else
	U = X + DSIGN(tol1,D);
      
      FU = FFF(U);
      if(FU <= FX)
	{
          if(U >= X)
            A = X;
	  else
            B = X;
	  
          V = W;
          FV = FW;
          W = X;
          FW = FX;
	  X = U;
          FX = FU;
	}
      else
	{
          if (U < X)
            A = U;
          else
            B = U;
          
          if (FU <= FW || W == X)
            {
	      V = W;
	      FV = FW;
	      W = U;
	      FW = FU;
	    }
          else if (FU <= FV || V == X || V == W)
            {
	      V = U;
	      FV = FU;
	    }
	  
	}
    }
  printf("Brent exceed maximum iterations.\n");
tre:     
  *xmin = X;
  return FX;
}

/* ============================ >>> F1dim <<< ========================== */
double F1dim(double X)
{
  /*const int Nmax=1000;*/
  int j, i, a;
  for (a = 0; a < NA; a++)
    {
      for( j = 0; j < Oparams.parnum[a]; j++)
	{
	  rx[a][j] = pcomx[a][j] + X*xicomx[a][j];
	  ry[a][j] = pcomy[a][j] + X*xicomy[a][j];
	  rz[a][j] = pcomz[a][j] + X*xicomz[a][j];
	  /*printf("X:%f xicomx:%f pcomx:%f\n", X, xicomx[j], pcomx[j]);*/
	}
    }
  /*
    check that all molecules are in the original box
    check of periodic boundary conditions on the read data
  */
  
  /* porta tutte le particelle nel primo box */
  for (a = 0;  a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
#ifdef MD_USE_RINT
	  rx[a][i] = rx[a][i] - L * rint(invL * rx[a][i]);
	  ry[a][i] = ry[a][i] - L * rint(invL * ry[a][i]);
	  rz[a][i] = rz[a][i] - L * rint(invL * rz[a][i]);
#else
	  nperx =(int) (fabs(rx[a][i])*invLH);                                         
	  anperx=DSIGN(((double)(nperx)*L),rx[a][i]);                                  
	  npery =(int) (fabs(ry[a][i])*invLH);                                          
	  anpery =DSIGN(((double)(npery)*L),ry[a][i]);                                  
	  nperz =(int) (fabs(rz[a][i])*invLH);           
	  anperz=DSIGN(((double)(nperz)*L),rz[a][i]);
	  rx[a][i] -= anperx;
	  ry[a][i] -= anpery;
	  rz[a][i] -= anperz;
#endif
	}   
    }
  /*ENERGY*/
  Forces(Oparams.parnum);
  
  return Epoten;
  /*  printf("x: %.6f f1dim: %.6f ncom:%d\n",x,f1dim,ncom)*/
}

/* ======================= >>> SDchkRebuild <<< ========================= */
void SDchkRebuild(void)
{
  int i, a;
  COORD_TYPE vv, vvMax = 0.0;

  for (a = 0; a < NA; a++)  
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  vv = Sqr(xix[a][i]) + Sqr(xiy[a][i]) + Sqr(xiz[a][i]);
	  if (vv > vvMax) 
	    vvMax = vv;
	}
    }
  dispHi = dispHi + sqrt(vvMax);
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;

}

/* ========================= >>> linmin <<< ======== ======================= */
void LinMin(int Nm[NA], double *Fret)
{
  int j, i, a;
  /*const nmax = 1000;*/
  const double tol = 1.0E-10;
  double AX, XX, BX, xmin, FX, FB, FA;
  /*int Ncom;*/
 
  /*Ncom = N;*/

  for (a = 0; a < NA; a++)
    {
      for (j = 0; j < Nm[a]; j++) 
	{
	  pcomx[a][j]  = rx[a][j];
	  pcomy[a][j]  = ry[a][j];
	  pcomz[a][j]  = rz[a][j];
	  xicomx[a][j] = xix[a][j];
	  xicomy[a][j] = xiy[a][j];
	  xicomz[a][j] = xiz[a][j];
	}
    }


  Forces(Nm);
  /*sprintf(TXT, "Prima di MNBRAK=%.6f\n", Epoten); 
    mdPrintf(ALL, TXT, NULL);
    printf("bx et al',box,boy,boz,boxhx,boxhy,boxhz*/
  AX = 0.0;
  XX = 0.001 / fnorm;
  BX = 0.002 / fnorm;
  
  mnbrak(&AX, &XX, &BX, &FA, &FX, &FB, F1dim);
  *Fret = brent(AX, XX, BX, F1dim, tol, &xmin);
  /*sprintf(TXT, "IN LINMIN: MIN: %.6f\n", xmin);
    mdPrintf(ALL, TXT, NULL);*/
  /*
    move the system
  */
  for (a = 0; a < NA; a++)
    {
      for(j = 0; j < Nm[a]; j++)
	{
	  xix[a][j] = xmin * xix[a][j];
	  xiy[a][j] = xmin * xiy[a][j];
	  xiz[a][j] = xmin * xiz[a][j];
	  /*if ((xmin == 0.0) && (xi[j] != 0.0)
	    ) printf *,j,k,xi(j,k)*/
	  rx[a][j] = pcomx[a][j] + xix[a][j];
	  ry[a][j] = pcomy[a][j] + xiy[a][j];
	  rz[a][j] = pcomz[a][j] + xiz[a][j];
	}
    }
  /* porta tutte le particelle nel primo box */
  for (a = 0; a < NA; a++)
    {
      for (i = 1; i < Nm[a]; i++)
	{
#ifdef MD_USE_RINT
	  rx[a][i]= rx[a][i] - L * rint(invL * rx[a][i]);
	  ry[a][i]= ry[a][i] - L * rint(invL * ry[a][i]);
	  rz[a][i]= rz[a][i] - L * rint(invL * rz[a][i]);
#else
	  nperx =(int) (fabs(rx[a][i])*invLH);                                         
	  anperx=DSIGN(((double)(nperx)*L),rx[a][i]);                                  
	  npery =(int) (fabs(ry[a][i])*invLH);                                          
	  anpery =DSIGN(((double)(npery)*L),ry[a][i]);                                  
	  nperz =(int) (fabs(rz[a][i])*invLH);           
	  anperz=DSIGN(((double)(nperz)*L),rz[a][i]);
	  rx[a][i] -= anperx;
	  ry[a][i] -= anpery;
	  rz[a][i] -= anperz;
#endif
	}   
    }
  /*printf(" XMIN= %.6f EN=%.6f\n", xmin, *Fret);
    */

  SDchkRebuild();

  Forces(Nm);
  /*printf("Dopo di MNBRAK=%.6f\n", Epoten); 
    */
  if (xmin == 0.0) 
    {
      printf("ATTENZIONE X MIN= 0\n");
    }

}

/* ======================= >>> frprmn <<< =========================== */
void frprmn(int Nm[NA], double Ftol, double *Fret)
{
  const int ITMAX = 20000;
  const double eps = 1E-14;
  double FP, GG, DGG, GAM;
  int a, j, iter;
  
  /*BuildNebrListNoLinked(Nm, Oparams.rcut, Oparams.sigma);*/
  Forces(Nm);
  FP = Epoten;
  /*
    calculate force modulus
  */
  fnorm = 0.0;

  for (a = 0; a < NA; a++)
    {
      for (j = 0; j < Nm[a]; j++)
	{
	  fnorm = fnorm + Fx[a][j]*Fx[a][j]+Fy[a][j]*Fy[a][j]+Fz[a][j]*Fz[a][j];
	}
    }
  fnorm = sqrt(fnorm);
  
  for (a = 0; a < NA; a++)
    {
      for(j = 0; j < Nm[a]; j++)
	{
	  Gx[a][j]  = Fx[a][j];
	  Hx[a][j]  = Gx[a][j];
	  xix[a][j] = Hx[a][j];
	  Gy[a][j]  = Fy[a][j];
	  Hy[a][j]  = Gy[a][j];
	  xiy[a][j] = Hy[a][j];
	  Gz[a][j]  = Fz[a][j];
	  Hz[a][j]  = Gz[a][j];
	  xiz[a][j] = Hz[a][j];
	}
    }
  
  for(iter = 0; iter < ITMAX; iter++)
    {
      LinMin(Nm, Fret);
      if (iter % ITER_ELAPSED == 0)
	{
	  printf("[#%d] ",iter);
	}
      Forces(Nm);
      
      if (iter % ITER_ELAPSED == 0)
	{
	  printf("E=%.15G\n", Epoten);
	}
      /*
	calculate force modulus
      */
      fnorm=0.0;
      for (a = 0; a < NA; a++)
	{
	  for(j = 0; j < Nm[a]; j++)
	    {
	      fnorm = fnorm + Fx[a][j]*Fx[a][j]+Fy[a][j]*Fy[a][j]+Fz[a][j]*Fz[a][j];
	    }
	}
      fnorm = sqrt(fnorm);
      /* SISTEMARE !!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10,*) iter,*Fret,fnorm;
       */
      if (2.0 * fabs(*Fret-FP) <= Ftol*(fabs(*Fret)+fabs(FP) + eps))
	return;
      
      FP = Epoten;
      
      GG = 0.0;
      DGG = 0.0;
      for (a = 0; a < NA; a++) 
	{
	  for(j = 0; j < Oparams.parnum[a]; j++)
	    {
	      GG = GG + Sqr(Gx[a][j]) + Sqr(Gy[a][j]) + Sqr(Gz[a][j]);
	      //DGG=DGG+F(J,K)**2
	      DGG = DGG + (-Fx[a][j] + Gx[a][j])*(-Fx[a][j]) + 
		(-Fy[a][j] + Gy[a][j])*(-Fy[a][j]) +
		(-Fz[a][j] + Gz[a][j])*(-Fz[a][j]);
	    }
	}
      if (GG == 0)
	return;
      
      GAM = DGG / GG;
      for (a = 0;  a < NA; a++)
	{
	  for(j = 0; j < Nm[a]; j++)
	    {
	      Gx[a][j]  = Fx[a][j];
	      Gy[a][j]  = Fy[a][j];
	      Gz[a][j]  = Fz[a][j];
	      Hx[a][j]  = Gx[a][j] + GAM * Hx[a][j];
	      Hy[a][j]  = Gy[a][j]+ GAM * Hy[a][j];
	      Hz[a][j]  = Gz[a][j]+ GAM * Hz[a][j];
	      xix[a][j] = Hx[a][j];
	      xiy[a][j] = Hy[a][j];
	      xiz[a][j] = Hz[a][j];    
	      //         per steepest discent
	    }
	}
    }
  printf("FRPR maximum iterations exceeded\n");
}


/* ==================== >>> mnbrak <<< =================================== */
void mnbrak(double *AX, double *BX, double *CX, double *FA, double *FB, double *FC,
	     double (*func)(double))
{
  const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0E-15;
  double Q, FU, U, ULIM, DUM, R;
  /*register double AX_, BX_, CX_, FA_, FB_, FC_;  */
  
  *FA = func(*AX);
  *FB = func(*BX);
  
  if( *FB > *FA)
    {
      DUM = *AX;
      *AX = *BX;
      *BX = DUM;
      DUM = *FB;
      *FB = *FA;
      *FA = DUM;
    }
  *CX = *BX + GOLD*(*BX-*AX);
  
  *FC = func(*CX);
 inizio:
  if (*FB >= *FC)
    {
      R = (*BX-*AX) * (*FB-*FC);
      Q = (*BX-*CX) * (*FB-*FA);
      U = *BX-((*BX-*CX)*Q - (*BX-*AX)*R)/(2.0*DSIGN(max(DABS(Q-R),TINY),Q-R));
      ULIM = *BX + GLIMIT*(*CX-*BX);
      if((*BX-U)*(U-*CX) > 0)
	{
	  FU = func(U);
          if (FU < *FC)
	    {
	      *AX = *BX;
	      *FA = *FB;
	      *BX = U;
	      *FB = FU;
	      goto inizio;
	    }
          else if(FU > *FB)
	    {
	      *CX = U;
	      *FC = FU;
	      goto inizio;
	    }
	  
	  U = *CX + GOLD*(*CX-*BX);
	  FU = func(U);
	}
      else if ((*CX-U)*(U-ULIM) > 0) 
	{
	  FU = func(U);
      
	  if (FU < *FC)
	    {
	      *BX = *CX;
	      *CX = U;
	      U = *CX + GOLD*(*CX-*BX);
	      *FB = *FC;
	      *FC = FU;
	      FU = func(U);
	    }
	}
      else if((U-ULIM)*(ULIM-*CX) >= 0)
	{
	  U = ULIM;
	  FU = func(U);
	}
      else
	{
	  U = *CX + GOLD*(*CX-*BX);
	  FU = func(U);
	}
      
      *AX = *BX;
      *BX = *CX;
      *CX = U;
      *FA = *FB;
      *FB = *FC;
      *FC = FU;
      goto inizio;
    }
}


/* =========================== >>> forces <<< ======================= */
void  Forces(int Nm[NA])
{
  int i, a;
  double L, invL, invLH;
  
  L = cbrt(Vol);
  invL = 1.0/L;
  invLH=1.0/(L/2.0);
  /* porta tutte le particelle nel primo box */
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Nm[a]; i++)
	{
#ifdef MD_USE_RINT
	  rx[a][i]=rx[a][i] - L*rint(invL*rx[a][i]);
	  ry[a][i]=ry[a][i] - L*rint(invL*ry[a][i]);
      	  rz[a][i]=rz[a][i] - L*rint(invL*rz[a][i]);
#else
	  nperx =(int) (fabs(rx[a][i])*invLH);                                         
	  anperx=DSIGN(((double)(nperx)*L),rx[a][i]);                                  
	  npery =(int) (fabs(ry[a][i])*invLH);                                          
	  anpery =DSIGN(((double)(npery)*L),ry[a][i]);                                  
	  nperz =(int) (fabs(rz[a][i])*invLH);           
	  anperz=DSIGN(((double)(nperz)*L),rz[a][i]);
	  rx[a][i] -= anperx;
	  ry[a][i] -= anpery;
	  rz[a][i] -= anperz;
#endif
	}   
    }

  
  /* azzera tutte le forze */
  /*setToZero(Oparams.parnum[0], Fx[0], Fy[0], Fz[0], 
    NULL);  
    
    setToZero(Oparams.parnum[1], Fx[0], Fy[0], Fz[0], 
    NULL); */

  /* Costruisce le neighbour list */
  
  if (nebrNow)
    {
      nebrNow = 0;
      dispHi = 0.0;
      BuildNebrListNoLinked(Oparams.rcut, Oparams.sigab);
    }

  LJForce(Oparams.epsab, Oparams.sigab, Oparams.rcut);
 
  /* Energia potenziale totale */
  Epoten = V; /* V = valore non shiftato Vc=valore non shiftato */

  /*printf(,"Epoten IN Forces: %.6f\n",Epoten);
    */
}

