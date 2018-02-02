#define BASINS
#define MAIN
#include<mdsimul.h>
extern void readAsciiPars(FILE* pfs, struct pascii strutt[]);
int SEGSIZE;
extern void readAllCor(FILE* fs);
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);
double L, invL, invLH;
struct params Oparams;
struct progStatus OprogStatus;
extern void readBakAscii(char* fn);
extern void saveBakAscii(char *fn);
struct simStat OsimStat;
int ENDSIM;
char msgStrA[MSG_LEN];
/* neighbour list method variables */
COORD_TYPE dispHi;
int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
char TXTA[10][MSG_LEN];
char TXT[MSG_LEN];
COORD_TYPE pi, s1t, Vol1t, L, invL, invLH, s1p, Elrc, Plrc;   
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

int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
extern void setToZero(int size, COORD_TYPE *arrx, ...);
int *head, *list, *map;  /* arrays of integer */
int NCell, mapSize, M;
unsigned char BAK, STA, BAKT;

extern void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigma[NA][NA]);
extern void LJForce(COORD_TYPE epsab[NA][NA], COORD_TYPE sigab[NA][NA], COORD_TYPE rcut);

void initAll(int[NA]);
void Forces(int[NA]);

/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintf(int mode, char *text, ...)
{return; }
/* ========================== >>> mdMsg <<< ===================================== */
void mdMsg(int mode, int errnum, char *when, char *errType, char* sysCall,
  char* text, ...)
{ return; }

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

/* ------------------------------------------------------------------------------- */
void readOne(char* fn)
{
  FILE* fs; 
  static int firstTime = 1;
  
  /*printf("reading: %s\n",fn);*/
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
#define njobmax 1000

char filename[njobmax][512],filewalter[512];
char fileslave[njobmax][512], percorso[512];
int steps[njobmax];

/* ========================== >>> readStart <<< ============================ */
void readRestart(char filewalter[132])
{
  /*readBakAscii(filewalter);*/
  readOne(filewalter);

  /*
    check of periodic boundary conditions on the read data
  */
  L = cbrt(Vol);
  invL = 1.0/L;
  invLH = 1.0/(L/2.0); 
}
double calcpress(void)
{
  BuildNebrListNoLinked(Oparams.rcut, Oparams.sigab);
  
  L = cbrt(Vol);
  invL = 1.0 / L;

  LJForce(Oparams.epsab, Oparams.sigab, Oparams.rcut);
  return (W / Vol);

}
int nacc[njobmax];
double acc[njobmax];
int tempi[njobmax];
double r1x[2000], r1y[2000], r1z[2000],
  r2x[2000], r2y[2000], r2z[2000];
int njob;
void init(void)
{
  static int firstTime = 1;
  nebrTabMax = OprogStatus.nebrTabFac * (Oparams.parnum[0]+Oparams.parnum[1])*NA;
  nebrNow = 1;
  if (firstTime)
    {
      /* Initialize variables for neighbour list method */
      fprintf(stderr,"INIT nebrTabMax: %d\n", nebrTabMax);
      nebrTab = AllocMatI(2, nebrTabMax); 
      firstTime = 0;
    }
}
int main(int argc, char **argv)
{
  FILE* pf;
  char listaFile[1024];
  double fpn;
  int PP=-1;
  int aggiungi, steps1, steps2;
  int ind, maxt, i, nn, in, t1, t2;
#if defined(SOFT_SPHERE)
#if defined(MD_STATIC_PP36)
  PP = 36;
  fprintf(stderr, "WARNING: MD_STATIC_PP36 defined!\n");
#elif defined(MD_STATIC_PP18)
  PP = 18;
  fprintf(stderr, "WARNING: MD_STATIC_PP18 defined!\n");
#elif defined(MD_STATIC_PP12)
  PP = 12;
  fprintf(stderr, "WARNING: MD_STATIC_PP12 defined!\n");
#elif defined(MD_STATIC_PP8)
  PP = 8;
  fprintf(stderr, "WARNING: MD_STATIC_PP8 defined!\n");
#elif defined(MD_STATIC_PP6)
  PP = 6;
  fprintf(stderr, "WARNING: MD_STATIC_PP6 defined!\n");
#else
  fprintf(stderr, "WARNING: DYNAMIC PP=%d\n", Oparams.PP);
#endif
#endif
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
  /*
     if( (finp = fopen(listaFile, "r"))==NULL )
     {
     fprintf(stderr,"Unable to open file: %s\n", listaFile);
     exit(-1);
    }
    */
  njob = 1;
  
  fscanf(pf, "%s\n", percorso);
  fprintf(stderr,"Leggo i file da quenchare dalla directory:%s\n", percorso);
  
  while(!feof(pf))
    {
      if (fscanf(pf, "%s", filename[njob]) == 1)
	{
	 /* printf("njob:%d\n", njob);*/
	  fprintf(stderr,"leggo file %s\n", filename[njob]);
	  njob++;
	}
    }
  fclose(pf);
  maxt = 0; 
  for (t1= 1; t1 < njob; t1++)
    {
      strcpy(filewalter, percorso);
      strcat(filewalter, filename[t1]);
      readRestart(filewalter);
      init();
#ifdef SOFT_SPHERE
      if (PP!=-1)
	Oparams.PP = PP;
#endif
      for (i=0; i < Oparams.parnum[0]; i++)
	{
	  r1x[i] = rx[0][i];
	  r1y[i] = ry[0][i];
	  r1z[i] = rz[0][i];
	}
      in  = Oparams.curStep;
      fpn = calcpress();
      /*printf("fpn: %f\n", fpn);*/
      aggiungi = 1;
      nn = maxt;
      for (i = 0; i < maxt; i++)
	{
	  if (in == tempi[i])
	    {
	      aggiungi = 0;
	      nn = i;
	      break;
	    }
	}
      
      if (aggiungi)
	{
	  tempi[maxt] = in;
	  fprintf (stderr, "aggiungo tempo %d [%d]\n", tempi[maxt], maxt);
	  acc[maxt] = fpn;
	  nacc[maxt] = 1;
	  maxt++;
	}
      else
	{
	  acc[nn] += fpn;
	  nacc[nn]++;
	}
    }
  for (ind = 0; ind < maxt; ind++)
    {
      fprintf(stdout, "%d %.15G %d\n", tempi[ind], 
	      acc[ind]/((double)nacc[ind]), 
	      nacc[ind]);
    }
  return 0;
}
