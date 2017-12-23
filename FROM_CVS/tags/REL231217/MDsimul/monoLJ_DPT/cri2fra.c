#define BASINS
#define MAIN
#include <mdsimul.h>
#include <mpi.h>

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
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);

struct params Oparams;
struct progStatus OprogStatus;

/* ======================= >>> DSIGN <<< ================================ */
double DSIGN(double a, double b)
{
  return (fabs(b)/b) * fabs(a);
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
int mdOpen_Tord( char* fileName, char* when, char* errMsg, int mode){return -1;}
void writeSegs(int fdes, char* when, char* errMsg, int mode, int size,
	       void* pointer, ...){}
struct measure Omeasure[]={{NULL,0,NULL}};
unsigned char sw(unsigned char onOff){return 0;}
void scaleVelocities(double fact){}

/* =========================== >>> saveRestart <<< ============================ */
void saveRestart(char *fn)
{
  FILE *of;
  of = fopen(fn, "w");
  fprintf(of, "%d %.15G %.15G %.15G", Oparams.curStep, Oparams.T, Oparams.P, Vol);
  writeAllCor(of);
  fclose(of);

}

/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  /* pure queste informazioni devono seguire il formato di Francesco */
  int i;
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", rx[i], ry[i], rz[i]);
    }
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", vx[i], vy[i], vz[i]);
    }

  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", Fx[i], Fy[i], Fz[i]);
    }

#if 0
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", vxo1[i], vyo1[i], vzo1[i]);
    }

  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", vxo2[i], vyo2[i], vzo2[i]);
    }
  
  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G %.15G",  Vol, s, s1, s2, s1o1, s1o2);
#endif

}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i;

  for (i = 0; i < Oparams.parnum; i++)
    {
      fscanf(fs, "%lf %lf %lf\n", &rx[i], &ry[i], &rz[i]);
    }
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      fscanf(fs, "%lf %lf %lf\n", &vx[i], &vy[i], &vz[i]);
    }
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      fscanf(fs, "%lf %lf %lf\n", &Fx[i], &Fy[i], &Fz[i]);
    }

  for (i = 0; i < Oparams.parnum; i++)
    {
      fscanf(fs, "%lf %lf %lf\n", &vxo1[i], &vyo1[i], &vzo1[i]);
    }

  for (i = 0; i < Oparams.parnum; i++)
    {
      fscanf(fs, "%lf %lf %lf\n", &vxo2[i], &vyo2[i], &vzo2[i]);
    }

  fscanf(fs, "%lf %lf %lf %lf %lf %lf", &Vol, &s, &s1, &s2, &s1o1, &s1o2);
  
}

/* ======================== >>> main <<< ================================= */
void main(int argc, char** argv)
{
  int size;

  /* il primo argomento è il nome del file da leggere */
  readBakAscii(argv[1]);
  printf("File letto nel formato mio\n");
  saveRestart(argv[2]);
  printf("File scritto nel formato di Francesco.\n");
}
