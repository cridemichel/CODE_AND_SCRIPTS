/* ============================ >>> mdsimul.h <<< =============================
   in questo file vengono definite le strutture dati necessarie alla
   simulazione, che pero' non dipendono dalla simulazione specifica.
   Le variabili specifiche vanno nel file incluso alla fine (ved. fine file) */
#include<time.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdio.h>
#include<signal.h>
#include<string.h>
#include<sys/types.h>
#include<sys/wait.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<dirent.h>
#include<math.h>
#include<stdarg.h>
#include<errno.h>
#include<limits.h>
#if defined(MPI)
#if defined(ORIGIN) || defined(CRAY)
#undef alpha
#include <mpi.h>
#else
//#include <mpich/mpi.h>
#include <mpi.h>
#endif
#endif
/* ============================== >>> DEFINES <<< ===========================*/
#ifdef CRAY
#include<fp.h>
/*in fp.h  c'e' rint(...) altrimenti rint non funziona*/
#define cbrt(X) pow(X, (double)1.0/3.0)

#endif
#if defined(CRAY) || defined(ORIGIN)
#define MPI_INTEGER MPI_INT
#define MPI_DOUBLE_PRECISION MPI_DOUBLE
#endif

#define MAXPROCS 512
#if defined(MPI)
#define MPI_MSG_LEN 512
#define MD_MPI_PRINT 1
#define MD_MPI_END 0
#define MD_MPI_EXCHG 128
#define RXMD_MPI 2
#define MPI_GC 100
#define MPI_OUT_PROCESS 0
//#define exit(VAL) MPI_Abort(MPI_COMM_WORLD, VAL)
#endif
/* 03/05/2010: massima lunghezza della linea nella readAsciiPars. Se si simulano molte particelle
questo numero deve essere abbastanza grande altrimenti la ripartenza da file ascii
da segfault. Considerare tuttavia che se MAXPAR è definita allora non viene usata questa
define ma il valore MAXPAR*20 */
#define MAXLINELEN 4096000 
#define MEMBLOCK 32000000
#define EXT_DLST          
#define EXT_SLST NULL
/* The user should override these defines in the simulation dependent 
   header file */

#define MDINT int

#define MAXSIZE 1000000     /* max size in bytes of a log file (CHK_LOG 
			       and MDS_LOG), if its length bacomes  greater 
			       than this value the file is truncated to zero 
			       length */

#define COORD_TYPE double     /* All floating point variables int the program 
				 use this type, i recommend strongly to use 
				 the double type, that in Pentium and Pentium 
				 Pro processors with Linux is even faster than 
				 float type (but slower than long double) */
#define NAME_LENGTH 1024        /* maximum filename length */

#define MAX_LENGTH NAME_LENGTH

#define MAXVARS 1024            /* maximum number of variable you can put as 
				   arguments in AllocSh() and AllocCoord() */
#define MAXSEGS  128            /* maximum number of segments allocated by 
				   AllocSh() or AllocCoord() */


#define BAK_FILE_NAME "COORD_TMP"
#define STATUS_FILE_NAME "STATUS_FILE"

#define CHK_LOG "check.log"
#define MDS_LOG "mdsimul.log"
#define CF "CHK_FILE"

/* used by mdsimul OUTPUT ROUTINES */ 
#define STD 1                 /* see mdPrintf in mdinit.c file */  
#define LOG 3
#define ALL 255 
#define MSG_LEN 2048           /* maximum length of error messages */

#define NOSYS -65533           /* it is important that it is negative */ 

/* used by PARSING PROCEDURE in file mdinit.c and in mdsimdep.h */
#define STR  1
#define INT  2
#define CT   4
#define LLINT 5

/* used by mdCreat */
#define EXIT 0
#define CONT 1


/* ============================= >>> MACROS <<< ============================ */

/* loopShr(<var>, <begin>, <end>) { ... } poolShr; */ 
#ifdef MAIN
int* i_;  
/* pointer to the SC counter, in this way you haven't to specify 
   SC at the end*/
int newSim; /* if 0 => new continuing */
int NUMCALCS; /* Before calling a measuring function this variable is set
		 to the number of calculations performed for that measure, 
		 this is usefule for doind averaging */
int mgl_mode=0; /* 0=salva nel formato non-mgl 1=salva i file ascii nel formato mgl 
		   2=salva nel formato mgl ed esce appena salvato al tempo iniziale senza 
		   avanzare la simulazione (serve per convertire Store o file Cor in file 
		   mgl.*/
char inifile_for_mgl[NAME_LENGTH];
int mdseed=0;
#else 
extern int* i_;
extern int newSim; /* if 1 => new simulation */
extern int NUMCALCS;
extern int mdseed;
extern int mgl_mode;
extern char inifile_for_mgl[NAME_LENGTH];
#endif

/* ========================== >>> SHARED LOOP <<< ==========================*/
/* If END  = 500 for example the first process (process=0) loops until 
   sc_[0] = 249, while the second process loops until sc_[1] = 250.
   This solution is not efficent because the two processes loop not 
   concurrently, anyway if the time to perfom a loop is not much longer 
   than the time of a lock/unlock this solution is better than 'loopShr' */ 

#define ProcSync0()

/* show a variable value on the screen you can use this macro everywhere,
   anyway generally it is useful to put this in tha 'move()' function or in
   'measuring functions'.
   FMT is the 'printf' format string and MDVAR is a variable you want to
   print */ 
#define mdShow(FMT, MDVAR) sprintf(msgStrA, FMT, MDVAR);\
                           mdMsg(STD, NOSYS, NULL, "SHOW", NULL,\
                           msgStrA, NULL);
							  
/* 'loop(i, 1, 100)' is equivalent to 'for(i=0; i<100; ++i)' 
   (index Fortran convention), the Syntax is:
   loop(<var>, <begin>, <end>)
     {
       ...
     }
   pool;
*/
#define loop(MDCOUNT, BEGIN, END) for(MDCOUNT = BEGIN - 1; MDCOUNT < END;\
                                  ++MDCOUNT)
#define pool

/* loopFat is like loop but this kind of loop is performed only by the 
   father, the Syntax is:
   loopFat(<counter>, <begin>, <end>)
   { 
   ... 
   }
   poolFat;
*/

#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */

/* ============================== >>> TYPES <<< =============================*/

/* backup files (or restore files)
   used to restart are named 
   (BAK_FILE_NAME)0 e (BAK_FILE_NAME)1 */

#define RINT_(X_) (X_>0?((int)(X_+0.5)):((int)(X_-0.5)))
typedef unsigned char bool;

/* ========================= >>> struct sim_stat <<< ========================*/
struct simStat
{
  /*This structure contains infos about simulation status and it is saved into 
    the file STATUS_FILE, used by the program chk_status to control if 
    simulation is gooing on */ 
#ifdef MDLLINT
  long long int curStep;	/* current simulation step */
#else
  int curStep;
#endif
  bool running;		        /* if true means that simulation is  
				   going on(unused)*/
  time_t totTime;		/* total simulation time
				   (unused) */
};

/* ========================= >>> struct chk_str <<< =========================*/
struct chkStr
{
#ifdef MDLLINT
  int curStep;
#else
  long long int curStep;
#endif
  bool whichSf;
};

/* ========================= >>> struct singlePar <<< =======================*/
struct singlePar 
{
  /* In mdsimul.c is istanziated an array of this structure, which is 
     initialized by the SINGLE_PAR macro, see 'mdsimdep.h' header file */
  char parName[NAME_LENGTH];
  void* ptr;
  unsigned char type; /* see mdsimdep.h for details */
};

/* ============================ >>> struct pascii <<< =================== */
struct pascii
{
  char parName[NAME_LENGTH];
  void* ptr;
  int qty;
  int block_length;
  char type[32];
};

/* ========================= >>> struct Measure <<< =======================*/
struct measure 
{
  /* DESCRIPTION:
     This structure is used to specify your measure.
     It is quite generic intentionally, in fact the user can 
     take advantage of this, taking whatever kind of measure he wants */
  void* buf;             /* pointer to the buffer */
  int size;     /* size of the buffer    */
  void (*calcFunc)(void);/* function that performs calculation of measure */
}; 

/* ========================== >>> convStruct <<< ============================*/
struct convStruct 
{
  /* DESCRIPTION:
     'type' is a string to use in the command line with option '-t' 
     to choose the converter to use ( for ex. md2ascii -t <type> ...)
     convereter is a pointer to a converter, that is to the function that 
     perform the conversion */
  char* type;        
#ifdef MDLLINT
  void (*converter)(FILE* afs, long long int step, int size);  
#else
  void (*converter)(FILE* afs, int step, int size);  
#endif
};

/*============================= >>> PROCEDURES <<< ==========================*/


/* DISK I/O 
 NOTE: generally load...() function name refers to function, that open, 
 read from and close a file, and save...() refers to a function thet open,
 write to, close a file.
 On the contrary read...() for example refers to a function, that read only a 
 file yet opened, ...*/
void chooseMeasure(char* absFile, int (*readFunc)(int));
int chooseRestore(char* absFile, unsigned char* PbakSwitch, 
		   int readFunc(int));
void chooseStatus(char* absFile, int (*readFunc)(int));
int openMPI(char*, int);
int openNewMPI(char*, int, int);
int creatMPI(char*, int);
FILE* fopenMPI(char*, char*);
int mdOpen( char* fileName, char* when, char* errMsg, int mode);
int creatWithHead(char* fileName, char* when, char* errMsg, int mode,
		   int headSize, void* header);
int mdCreat(char* fileName, char* when, char* errMsg, int mode);
int readFile(int fdes);
int mdWrite(int fdes, char* when, char* errMsg, int mode, 
	     int size, void *punta_arr);
int mdRead(int fdes, char* when, char* errMsg, int mode,
		     int size, void *pointer);
int mdClose(int fd, char* when, char* errMsg, int mode);
void writeSegs(int fdes, char* when, char* errMsg, int mode, 
	       int size, void* pointer, ...);
int readSegs(int fdes, char* when, char* errMsg, int mode,
	     int size, void* pointer, ...);
int readBak(int fdes);
void saveCoord(char* fileName);
int readCoord(int cfd); 
void saveBak(char *fileName);
int saveMeasure(int PN, char* fileName, int misNum, char* msgOpen, char* msgWrite, 
		char* msgClose);
void saveOneXva(char* fileName);
//void doubleSaveMeasure(int PN, int misNum, int* times, 
//       int tapeTimes);
void doubleBufStatus(void);
void doubleBufBak(unsigned char* hdWhich, unsigned char* tapeWhich, 
		  int* times, int tapeTimes);
void saveXva(int* times, int tapeTimes);
void openLog(char* mode);
void existDir(char* dirName);
int fileExists(struct dirent*, const char*);
int  TryOlderFile(int bf1, int (*readFunc)(int), 
		  int which);

int MPIgetchar();
void  _printing(void);
/* FILE MANAGEMENT */
unsigned char sw(unsigned char onOff);
char* appSw(char* fileName,const unsigned char which);
char* absTmpTape(const char* fileSrc);
char* absTmpHD(const char* fileSrc);
char* absTmpAsciiHD(const char* fileSrc);
char* absMisHD(const char* fileSrc);
char* absMisTape(const char* fileSrc);
char* absXvaTape(const char* fileSrc);
char* absXvaHD(const char* fileSrc);
void delTmpHD(char *relName);
void delDoubleHD(char *fileName);
void delDoubleTape(char *fileName);
char* sum(char *str, ...);
int copy(char* src, char* dest);
#ifdef MDLLINT
void pushSteps(long long int *saveSteps);
void popSteps(long long int* saveSteps);
#else
void pushSteps(int *saveSteps);
void popSteps(int* saveSteps);
#endif
/* PROCESS CONTROL */
void SuspendJob(int sigid);
#if defined(ORIGIN) || defined(CRAY) 
void endJob(void);
#else
void endJob(int, void*);
#endif
/* MEMORY MANAGEMENT */
//void AllocCoord(int size, COORD_TYPE** pointer, ...);
//void DeleteCoord(COORD_TYPE** pointer, ...);
void initBefore(void);
void usrInitBef(void);
void usrInitAft(void);
#if defined(CRAY) || defined(ORIGIN)
void calc(void);
void move(void);
#else
void calc(void);
void move(void);
#endif
/* MEMORY MANAGEMENT */
/* See Rapaport pag. 364 */
COORD_TYPE* AllocVecR(int size);
void FreeVecR(COORD_TYPE* v);
int* AllocVecI(int size);
void FreeVecI(int* v);
COORD_TYPE** AllocMatR(int size1, int size2);
void FreeMatR(COORD_TYPE** v);
int** AllocMatI(int size1, int size2);
void FreeMatI(int** v);

/* INTIALIZATION */
void Newsimul(char *argv);
void Continue(void);
int progExist(char *name);
int argsMd(int argc, char** argv);
void AllocCT(int num, COORD_TYPE** ptr, ...);
void AllocInt(int num, int** ptr, ...);
void initCoord(void);
//void setToZero(COORD_TYPE* ptr, ...);

/* MESSAGE HANDLING */ 
void mdPrintfWr(int mode, char *text);
void mdPrintfR0(int mode, char *text, ...);
void mdPrintf(int mode, char* text, ...);
void mdPrintfSpc(int mode, char* text);
void mdPrintfSpcWr(int mode, char* text);
void mdMsg(int mode, int errnum, char *when, char *errType, char* sysCall,
	   char* text, ...);
void delDoubleTape(char* fileName);
void delDoubleHD(char* fileName);
int MDprintf(int mode, const char *format,  ...);
/* SIMULATION */
void calc(void);
void initStruct(void);
void usrInitBef(void);
void usrInitAft(void);
void initBefore(void);
void sumAllProc(int nProc, COORD_TYPE* addendum);
void sumAllProcMA(COORD_TYPE* addendum, ...);

int chkXvaSteps(void);
#define UPDATE_SYSTEM 
#define ADJUST_LASTCOL
#define MD_EXT_INIT
#define MD_EXT_END
/*===========================================================================*/

/* ============================ >>> mdsimdep.h <<< ==========================*/

#ifdef LMC
#include<mdsimdepLMC.h>
#else
#include<mdsimdep.h>            /* simulation dependent header file */
#endif
/*===========================================================================*/
