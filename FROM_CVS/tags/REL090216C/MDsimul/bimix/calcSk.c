#define BASINS
#define MAIN
#include<mdsimul.h>
#define NJOBMAX 1000
char filename[NJOBMAX][132];
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
char filewalter[512];
int mesh[][150][3]=
#include "./kmesh.dat"
int ntripl[]=
#include "./ntripl.dat"
int numk = 98;
int begk = 0;
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
double Sk[100];
double scalFact, L, invL;
void calcSk(void)
{
  int Na, Nb, ia, ib, n, mp;
  double sumRho, invNa, rCMk, reRho, imRho;
  double pi2, kbeg;

  Na = Oparams.parnum[0];
  Nb = Oparams.parnum[1];

  L = cbrt(Vol);
  invL = 1.0 / L;
  pi = acos(0)*2.0;
  pi2 = 2.0 * pi;
  kbeg = 0.0; /*pi2 * invL;*/
  scalFact = pi2 * invL;
  invNa = 1.0 / Na;
  for(n = begk; n < numk; n++)
    {
      sumRho = 0.0;
      for(mp = 0; mp < ntripl[n]; mp++)
	{
	  reRho = 0.0;
	  imRho = 0.0;
	  for(ia=0; ia < Na; ia++)
	    {
	      /* il passo della mesh e' 0.5*pi2/L */
	      if (mesh[n][mp][0]==0 && mesh[n][mp][1] == 0 && 
		  mesh[n][mp][2] == 0)
		{
		  fprintf(stderr,"ERRORE nella MESH!!!!!!!! n=%d mp=%d ntripl[n]:%d\n", n,
			 mp, ntripl[n]);
		  exit(-1);
		}
	      rCMk = kbeg + scalFact * 
		(rx[0][ia] * mesh[n][mp][0] + 
		 ry[0][ia] * mesh[n][mp][1] + 
		 rz[0][ia] * mesh[n][mp][2]);
	      reRho = reRho + cos(rCMk) ; 
	      imRho = imRho + sin(rCMk);
	      /* Imaginary part of exp(i*k*r) for the actual molecule*/
	    }
	  for(ib=0; ib < Nb; ib++)
	    {
	      /* il passo della mesh e' 0.5*pi2/L */
	      if (mesh[n][mp][0] == 0 && mesh[n][mp][1] == 0 && 
		  mesh[n][mp][2] == 0)
		{
		  fprintf(stderr,"ERRORE nella MESH!!!!!!!! n=%d mp=%d ntripl[n]:%d\n", n,
			 mp, ntripl[n]);
		  exit(-1);
		}
	      rCMk = kbeg + scalFact * 
		(rx[1][ib] * mesh[n][mp][0] + 
		 ry[1][ib] * mesh[n][mp][1] + 
		 rz[1][ib] * mesh[n][mp][2]);
	      reRho = reRho + cos(rCMk) ; 
	      imRho = imRho + sin(rCMk);
	      /* Imaginary part of exp(i*k*r) for the actual molecule*/
	    }
	  sumRho = sumRho + reRho*reRho + imRho*imRho;
	}
      Sk[n] += sumRho / ((double)(Na+Nb));  
    }
}
int main(int argc, char **argv)
{
  int ind, njob, iwork;
  FILE* pf, *finp;
  char listaFile[132];

  if (argc < 2)
    {
      printf ("Devi fornire il nome del file Cnf!\n");
      exit(-1);
    }
  fprintf(stderr,"Leggo: %s\n", filewalter);
  for (ind=0; ind < numk; ind++)
    {
      Sk[ind] = 0.0;
    }
  if (argc >= 2)
    {
      if ((pf = fopen(argv[1],"r"))==NULL)
	{
	  fprintf(stderr,TXT, "Unable to open file: %s\n", argv[1]);
	  exit(-1);
	}
      /* printf("arg: %s\n",argv[1]); */
      strcpy(listaFile, argv[1]);
    }
  else
    {
      fprintf(stderr,"Devi fornire come argomenti il nome del file contenente il percorso\n");
      fprintf(stderr,"e la lista di files!\n");
    }


  if ( (finp = fopen(listaFile, "r"))==NULL )
    {
      fprintf(stderr, "Unable to open file: %s\n", listaFile);
      exit(-1);
    }
 
  if (argc==3)
    {
      numk = atoi(argv[2]);
    }
  else if (argc==4)
    {
      begk = atoi(argv[2]);
      numk = atoi(argv[3]);
    } 
  njob = 0;
      
  /* NOTA:
     - la prima riga deve contenere il  percorso dei file da analizzare
     - le altre righe devono contenere i nomi dei file senza il percorso */
  while(!feof(finp))
    {
      if (fscanf(finp, "%s", filename[njob]) == 1)
	{
	  njob++;
	}
    }
  
  fprintf(stderr, " Number of restart files=%d\n", njob);
  /* read all configurations */
  for(iwork = 0; iwork < njob; iwork++)
    {
      /*sscanf(filename[iwork], "%[^_]_R%d", stri, &rango);*/
      fprintf(stderr, "filename[iwork]: %s\n", filename[iwork]);
      fprintf(stderr,"iwork:%d\n", iwork);
      readRestart(filename[iwork]);
      calcSk();
    }
  for (ind = 0; ind < numk; ind++)
    {
      fprintf(stdout, "%f %.15G\n",  scalFact * ((COORD_TYPE)ind*0.5+1.25), 
	      Sk[ind]/((double)ntripl[ind])/((double)njob));
    }
  return 0;
}
