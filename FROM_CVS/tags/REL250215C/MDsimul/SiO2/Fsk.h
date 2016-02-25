/* In questo file si trovano tutte le variabili che devono essere definte 
   per poter linkare il programma buildPE ai files:
   mdarray_DPT.c forces_DTP.c mdini_DPT.c */
/*=== Alcune di queste sono variabili fittizie definite per poter linkare il 
  file mdarray.c === */
struct simStat OsimStat;
int ENDSIM;
char msgStrA[MSG_LEN];
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE EE[MAX_M], Vc, V, W, K,
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
C_T VAC, WAC, VLJ, WLJ;
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

//double  *rx, *ry, *rz, *Fx, *Fy, *Fz, Vol;
double pi;
double *kcx, *kcy, *kcz;

char pwd[512];

int my_rank, iniBakFormat, iniCorFormat;
/* ============ >>> neighbourlist <<< ===================== */
//COORD_TYPE dispHi;
int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;

void initAll(int);
void Forces(int);
void cleanstring(char *);
void readRestart(char filewalter[132]);
void frprmn(int Nm, double Ftol, double *Fret);
void writeconf(char filewalter2[132]);
void mnbrak(double *AX, double *BX, double *CX, double *FA, double *FB, double *FC,
	    double (*func)(double));
void writeEnergy(char filewalter[132], double T, int steps,
		 double VpotIni, double VpotEnd, double VACIni, double VACmin);
extern void zeroArrays(COORD_TYPE *arrx, COORD_TYPE *arry, COORD_TYPE *arrz, int N);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void ACForce(int Nm, COORD_TYPE kmax, COORD_TYPE Dk, C_T alpha, C_T S0, 
		    double Lambda);
extern void LJForce(int Nm, COORD_TYPE epsilon, COORD_TYPE sigma, COORD_TYPE rcut, 
		    double Lambda);
extern void readBakAscii(char* fn);
extern void saveBakAscii(char *fn);
extern int** AllocMatI(int size1, int size2);

/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintf(int mode, char *text, ...)
{
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;

  va_start(ap, text);

  strcpy(tstri, text);
  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      //printf("qui!\n");

      strcat(tstri, sptr);
    }
  
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

  va_start(ap, text);
  strcpy(tstri, text);
  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      //printf("qui!\n");

      strcat(tstri, sptr);
    }

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
int mdOpen_Tord( char* fileName, char* when, char* errMsg, int mode){return -1;}
void writeSegs(int fdes, char* when, char* errMsg, int mode, int size,
	       void* pointer, ...){return;}
struct measure Omeasure[]={{NULL,0,NULL}};
unsigned char sw(unsigned char onOff){return 0;}
void scaleVelocities(double fact){}
/* ---------------------------------------------------------------------------- */



