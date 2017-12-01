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
char filewalter[512];
int mesh[][150][3]=
#include "./kmesh.dat"
int ntripl[]=
#include "./ntripl.dat"
const int numk = 98;

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
double calcSk(void)
{
  int Na, Nb, ia, ib, n, mp;
  double sumRho, invNaNb, rCMk, reRhoA, imRhoA, reRhoB, imRhoB;
  double pi2, kbeg;

  Na = Oparams.parnum[0];
  Nb = Oparams.parnum[1];

  L = cbrt(Vol);
  invL = 1.0 / L;
  pi = acos(0)*2.0;
  pi2 = 2.0 * pi;
  kbeg = 0.0; /*pi2 * invL;*/
  scalFact = pi2 * invL;
  invNaNb = 1.0 / pow(Na*Nb,0.5);
  for(n = 0; n < numk; n++)
    {
      sumRho = 0.0;
      for(mp = 0; mp < ntripl[n]; mp++)
	{
	  reRhoA = 0.0;
	  imRhoA = 0.0;
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
	      reRhoA = reRhoA + cos(rCMk) ; 
	      imRhoA = imRhoA + sin(rCMk);
	      /* Imaginary part of exp(i*k*r) for the actual molecule*/
	    }
	  reRhoB = 0.0;
	  imRhoB = 0.0;
	  for (ib=0; ib < Nb; ib++)
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
		(rx[1][ib] * mesh[n][mp][0] + 
		 ry[1][ib] * mesh[n][mp][1] + 
		 rz[1][ib] * mesh[n][mp][2]);
	      reRhoB = reRhoB + cos(rCMk) ; 
	      imRhoB = imRhoB + sin(rCMk);
	      /* Imaginary part of exp(i*k*r) for the actual molecule*/

	    }
	  sumRho = sumRho + reRhoA*reRhoB + imRhoA*imRhoB;

	}
      Sk[n] = sumRho  * invNaNb/ ((COORD_TYPE) ntripl[n]);  
    }
}
int main(int argc, char **argv)
{
  int ind;
 
  if (argc < 2)
    {
      printf ("Devi fornire il nome del file Cnf!\n");
      exit(-1);
    }
  strcat(filewalter, argv[1]);
  fprintf(stderr,"Leggo: %s\n", filewalter);
  readRestart(filewalter);
  calcSk();
  for (ind = 0; ind < 98; ind++)
    {
      fprintf(stdout, "%f %.15G\n",  scalFact *((COORD_TYPE)ind), Sk[ind]);
    }
  return 0;
}
