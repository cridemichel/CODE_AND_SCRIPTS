#define BASINS
#define MAIN
#include<mdsimul.h>

extern void readBakAscii(char* fn);
extern void saveBakAscii(char *fn);

int SEGSIZE;        /* length in bytes of an array */
int my_rank;

/* string used to send  messages */
char msgStrA[MSG_LEN],msgStrB[MSG_LEN],msgStrC[MSG_LEN]; 
								   
unsigned char BAK, STA, BAKT;
/* ===========================================================================*/

extern void writeAllCor(FILE* );
extern void readAllCor(FILE* );
extern void AllocCoord(int size, int** pointer, ...);

struct params Oparams;
struct progStatus OprogStatus;

/* ======================= >>> DSIGN <<< ================================ */
double DSIGN(double a, double b)
{
  if (b > 0)
    {
      return fabs(a);
    }
  else
    {
      return -fabs(a);
    }
}
double *kcx, *kcy, *kcz;

struct simStat OsimStat;
/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintf(int mode, char *text, ...){}
/* ========================== >>> mdMsg <<< ===================================== */
void mdMsg(int mode, int errnum, char *when, char *errType, char* sysCall,
  char* text, ...){}

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
int mdOpen_Tord( char* fileName, char* when, char* errMsg, int mode){return -1;}
void writeSegs(int fdes, char* when, char* errMsg, int mode, int size,
	       void* pointer, ...){}
struct measure Omeasure[]={{NULL,0,NULL}};
unsigned char sw(unsigned char onOff){return 0;}
void scaleVelocities(double fact){}

/* =========================== >>> saveRestart <<< ============================ */
void saveRestart(char *fn)
{
  double lato;
  FILE *of;
  of = fopen(fn, "w");
  /* QUI  VA L'HEADER DI FRANCESCO */
  lato = Oparams.lattice_a * Oparams.lattice_M;
  printf("lato: %.6f\n", lato);
  fprintf(of, "%d %d %d %.15f\n", Oparams.curStep, Oparams.curStep, Oparams.parnum[0]+Oparams.parnum[1], lato);
  fprintf(of, "%.15f %.15f\n", lato, lato);
  fprintf(of, "%.15E %.15f %.15E\n", 0.0, Oparams.T, 0.0);
  fprintf(of, "%.15f %.15E %d\n", 0.0, 0.0, Oparams.curStep); 
  fprintf(of, "%.15f %d\n", 0.0, Oparams.curStep);
  writeAllCor(of);
  fclose(of);

}

/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  /* pure queste informazioni devono seguire il formato di Francesco */
  int i,a;
  double la;

  la = Oparams.lattice_a;

  for (a = 0; a < NA; a++)
  {
    for (i = 0; i < Oparams.parnum[a]; i++)
      {
	fprintf(fs, "%.15G %.15G %.15G\n", la*((double)rx[a][i]), 
		la*((double)ry[a][i]), la*((double)rz[a][i]));
      }
  }
}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i, a;
  double xx, yy, zz, la;
  
  la = Oparams.lattice_a;

  for(a = 0;  a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &xx, &yy, &zz) < 3)
	    {
	      printf("ERROR[pos] reading ascii file\n");
	      exit(-1);
	    }
	  rx[a][i] = (int) rint(xx/la);
	  ry[a][i] = (int) rint(yy/la);
	  rz[a][i] = (int) rint(zz/la);
	}
    }
}

/* ======================== >>> main <<< ================================= */
void main(int argc, char** argv)
{
  /* il primo argomento � il nome del file da leggere */
  readBakAscii(argv[1]);
  printf("File letto nel formato mio\n");
  saveRestart(argv[2]);
  //printf("Vol: %f\n", Vol);
  printf("File scritto nel formato di Francesco.\n");
}
