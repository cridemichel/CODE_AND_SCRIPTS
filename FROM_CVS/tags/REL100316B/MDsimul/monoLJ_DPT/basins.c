/* CONJGATE GRADIENT  
   VERSIONE PER BKS con PETER'S PATHCES
*/
#define BASINS
#define MAIN
#include<mdsimul.h>
#include<mpi.h>

#define ITER_ELAPSED 10
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

/* ======================= >>> DSIGN <<< ================================ */
double DSIGN(double a, double b)
{
  if (b > 0)
    return fabs(a);
  else
    return -fabs(a);
  
  //return (fabs(b)/b) * fabs(a);
}

#define DABS fabs

MPI_Status status;
int dietag;
#define njobmax 1000
//const int njobmax = 10000; 
/* numero massimo di files che si possono processare */ 

char filename[njobmax][132],filewalter[132];
char fileslave[njobmax][132], percorso[132];

int naccslave[njobmax];

const int iwtagHdr = 1, iwtagCor = 2;
int my_rank, islave;
int numOfProcs, njob; 
double tempe[njobmax], Ftol, Epoten, EpotIni, Emin, VACIni, VACmin, L, invL, fnorm;
int rangsl, numwrk, rango[njobmax], steps[njobmax];


//double Vc, V, W, VAC, WAC, VLJ;

double *xicomx, *xicomy, *xicomz, *pcomx, *pcomy, *pcomz, *xix, *xiy, 
  *xiz, *Gx, *Gy, *Gz, *Hx, *Hy, *Hz;

//double  *rx, *ry, *rz, *Fx, *Fy, *Fz, Vol;
double pi;
double *kcx, *kcy, *kcz;

char pwd[512];
char logFile[128][njobmax];

int iniBakFormat, iniCorFormat;
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
extern void readAsciiPars(FILE* pfs, struct pascii strutt[]);
extern void saveBakAscii(char *fn);
extern int** AllocMatI(int size1, int size2);

int chistampa = 0;
/* ============================ >>> mdPrintfNr <<< ========================*/
void mdPrintf(int mode, char *text, ...)
{
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;
  FILE* of;
  int quanti;

  va_start(ap, text);

  strcpy(tstri, text);
  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      //printf("qui!\n");

      strcat(tstri, sptr);
    }

 
  if (my_rank != 0)
    {
      quanti = strlen(tstri)+1;
      MPI_Send (&quanti, 1, MPI_INT, 0, 0,  MPI_COMM_WORLD );
      MPI_Send( &numwrk, 1, MPI_INT, 0, 0,  MPI_COMM_WORLD );
      MPI_Send( tstri, quanti, 
		MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);/*invia dati tipo=0!!!*/
      
    } 
  else
    {
      of = fopen(logFile[chistampa], "a");
      fprintf(of, tstri); 
      fclose(of);
      /* e anche su stdout */
      printf("[%d]%s",my_rank,tstri);
    }
 
  va_end(ap);    
    
}
/* ========================== >>> mdMsg <<< ===================================== */
void mdMsg(int mode, int errnum, char *when, char *errType, char* sysCall,
  char* text, ...)
{
  va_list ap;
  char tstri[MSG_LEN];
  char* sptr;
  FILE* of;
  int quanti;

  va_start(ap, text);
  strcpy(tstri, text);
  /* Concatenate the n char* args to obtain a single string */ 
  for (;(sptr = va_arg(ap, char*)) != NULL;)
    {
      //printf("qui!\n");

      strcat(tstri, sptr);
    }


  if (my_rank != 0)
    {
      quanti = strlen(tstri)+1;
      MPI_Send (&quanti, 1, MPI_INT, 0, 0,  MPI_COMM_WORLD );
      MPI_Send( &numwrk, 1, MPI_INT, 0, 0,  MPI_COMM_WORLD );
      MPI_Send( tstri, quanti, 
		MPI_CHAR, 0, 0, MPI_COMM_WORLD);/*invia dati tipo=0!!!*/
    }
  else
    {
      of = fopen(logFile[chistampa], "a");
      fprintf(of, tstri); 
      fclose(of);
    }
  
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
	       void* pointer, ...){}
struct measure Omeasure[]={{NULL,0,NULL}};
unsigned char sw(unsigned char onOff){return 0;}
void scaleVelocities(double fact){}
/* ---------------------------------------------------------------------------- */


/* ========================= >>> sendHeader <<< ============================= */
void sendHeader(int rank)
{

  MPI_Send((char*)&Oparams, sizeof(Oparams), MPI_CHAR, rank, iwtagHdr, 
	   MPI_COMM_WORLD);
  
  MPI_Send((char*)&OprogStatus, sizeof(OprogStatus), MPI_CHAR, rank, iwtagHdr, 
	   MPI_COMM_WORLD);

} 

/* ========================= >>> sendHeader <<< ============================= */
int recvHeader(void)
{
  /* l'header lo possono ricevere solo gli slave */
  int sorgente;
  MPI_Status status;
  
  MPI_Recv((char*)&Oparams, sizeof(Oparams), MPI_CHAR, 0,
	   iwtagHdr, MPI_COMM_WORLD, &status);
  sorgente = status.MPI_SOURCE;
  MPI_Recv((char*)&OprogStatus, sizeof(OprogStatus), MPI_CHAR, 0,
	   iwtagHdr, MPI_COMM_WORLD, &status);

  return sorgente;
}
 

/* -----------------------------------------------------------------------------------
   NOTA BENE: Queste due routine dipendono dalla simulazione!  */

/* =========================== >>> sendCoord <<< ============================ */
void sendCoord(int rank, int Nm)
{
  MPI_Send(&Nm, 1, MPI_INTEGER, rank, iwtagCor, MPI_COMM_WORLD);
  //printf("invio Nm: %d\n", Nm);
  MPI_Send(&Vol, 1, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD);
  MPI_Send(rx, Nm, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD);
  MPI_Send(ry, Nm, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD);
  MPI_Send(rz, Nm, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD);
}

/* =========================== >>> recvCoord <<< ============================= */
void recvCoord(int rank)
{
  /* le coordinate le possono ricevere sia gli slave che il master e mentre
     il master le puo' ricevere da qualsiasi slave, gli slave le possono ricevere
     solo dal master */
  int Nm;
  MPI_Status status;

  MPI_Recv(&Nm, 1, MPI_INTEGER, rank, iwtagCor, MPI_COMM_WORLD, &status);
  //printf("ricevo %d cord\n", Nm);
  MPI_Recv(&Vol, 1, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD, &status);
  MPI_Recv(rx, Nm, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD, &status);
  MPI_Recv(ry, Nm, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD, &status);
  MPI_Recv(rz, Nm, MPI_DOUBLE_PRECISION, rank, iwtagCor, MPI_COMM_WORLD, &status);
  //printf("[%d]ricevute\n", my_rank);
}

/* --------------------------------------------------------------------------------*/

/* =============================== >>> ricevi <<< =============================== */
void ricevi(int *iriceve, int *iwork_recv)
{
  int quanti;
  char testo[4096];
  do
    { 
      //printf("aspetto...\n");
      MPI_Recv(&quanti, 1, MPI_INT, MPI_ANY_SOURCE, 
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      *iriceve = status.MPI_SOURCE;
      //printf("ricevuto da %d\n", *iriceve);	      
      
      //printf("GIA' QUI?!?\n");
      if (quanti == -1)
	break;
      MPI_Recv(&chistampa, 1, MPI_INT, *iriceve, 
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      /* se quanti = 0 vuol dire che si tratta di dati altrimenti contiene il 
	 numero di caratteri della stringa da ricevere */
      MPI_Recv(testo, quanti,  MPI_UNSIGNED_CHAR, *iriceve, MPI_ANY_TAG, 
	       MPI_COMM_WORLD, &status);
      mdPrintf(ALL, testo, NULL);
      chistampa = 0;
    }
  while (1);
  
  MPI_Recv(&EpotIni, 1,  MPI_DOUBLE_PRECISION, 
	   *iriceve, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
  MPI_Recv(&Emin, 1,  MPI_DOUBLE_PRECISION, 
	   *iriceve, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&VACIni, 1,  MPI_DOUBLE_PRECISION, 
	   *iriceve, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&VACmin, 1,  MPI_DOUBLE_PRECISION, 
	   *iriceve, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(iwork_recv, 1,  MPI_INTEGER, 
	   *iriceve, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  
  recvCoord(*iriceve);
}

/* ======================= >>> invia_lavoro <<< ============================= */
void invia_lavoro()
{


}

/* ============================= >>> main <<< ============================== */
int main(int argc, char** argv)
{
  int NUM_WORK_REQUEST;
  int missing, ndone = 0, npr, iwork_recv, iwork, iriceve;
  //int steps, stepsok;
  FILE* finp;
  FILE* pf;
  char listaFile[1024], stri[512];
  int quanti = 0;
  double Fret;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs );
  
  islave = numOfProcs - 1; 
  if (my_rank == 0)  
    {
      sprintf(TXT, "MASTER%d.log", getpid());
      strcpy(logFile[0], TXT);
      
      sprintf(TXT, "Ci sono # slaves=%d\n", islave);  
      mdPrintf(ALL, TXT, NULL);
    }
  printf("my_rank: %d pid:%d\n", my_rank, getpid());
  
  if (my_rank == 0)
    {
      /* L'unico argomento deve essere la lista di file da minimizzare */
      if (argc == 2)
	{
	  if ((pf = fopen(argv[1],"r"))==NULL)
	    {
	      sprintf(TXT, "Unable to open file: %s\n", argv[1]);
	      mdPrintf(ALL, TXT, NULL);
	      MPI_Finalize();
	      exit(-1);
	    }
	  
	  // printf("arg: %s\n",argv[1]);
	  strcpy(listaFile, argv[1]);
	}
      else
	{
	  sprintf(TXTA[0],"Devi fornire come argomenti il nome del file contenente il percorso\n");
	  sprintf(TXTA[1],"e il file contenente la lista di files!\n");
	  mdPrintf(ALL, TXTA[0], TXTA[1], NULL);
	  MPI_Finalize();
	  exit(-1);
	}
    }
  
  if (my_rank == 0)
    {
      if( (finp = fopen(listaFile, "r"))==NULL )
	{
	  sprintf(TXT, "Unable to open file: %s\n", listaFile);
	  mdPrintf(ALL, TXT, NULL);
	  MPI_Finalize();
	  exit(-1);
	}
      
      njob = 1;
      
      /* NOTA:
	 la prima riga deve contenere il percorso dei file da quenchare mentre
	 le altre righe devono contenere i nomi dei file senza il percorso */
      
      fscanf(finp, "%s", percorso);
      sprintf(TXT, "Leggo i file da quenchare dalla directory:%s\n", percorso);
      mdPrintf(ALL, TXT, NULL);
      while(!feof(finp))
	{
	  if (fscanf(finp, "%s", filename[njob]) == 1)
	    {
	      //printf("file:%s\n", filename[njob]);
	      njob++;
	    }
	}

      njob--;

      sprintf(TXT, " Number of conf to quench= %d\n", njob);
      mdPrintf(ALL, TXT, NULL);
      NUM_WORK_REQUEST = njob;


      for(iwork = 1; iwork <= islave; iwork++)
	{
	  //printf("filename:%s\n", filename[njob]);
	  sscanf(filename[iwork], "%[^_]_R%d", stri, &rango[iwork]);
	  sprintf(TXT, "filename[iwork]: %s rango: %d\n", filename[iwork], rango[iwork]);
	  mdPrintf(ALL, TXT, NULL);
	  strcpy(filewalter, percorso);
	  strcat(filewalter, filename[iwork]);
	  strcpy(fileslave[iwork], filewalter);
	  //printf("iwork:%d filewalter: %s\n", iwork, filewalter);
	  
	  readRestart(filewalter); 
	  tempe[iwork] = Oparams.T / Oparams.lambda0[Oparams.lambdat[rango[iwork]]];
	  steps[iwork] = Oparams.curStep;

	  //printf("[%d]rango: %d perc:%s\n",my_rank, rangsl, percorso);
	  sprintf(logFile[iwork],"%sQncT%.6G.%d_R%dIW%d.log", 
		  percorso, 
		  Oparams.T/Oparams.lambda0[Oparams.lambdat[rango[iwork]]], 
		  Oparams.curStep, rango[iwork], iwork);
	  
	  sprintf(TXT, "logFile[%d]: %s\n", iwork, logFile[iwork]);
	  mdPrintf(ALL, TXT, NULL);

	  //printf("iwork:%d\n", iwork);
          sendHeader(iwork);
 	  sendCoord(iwork, Oparams.parnum);
	  MPI_Send(&rango[iwork], 1, MPI_INTEGER, iwork, iwtagCor,
		   MPI_COMM_WORLD);
	  MPI_Send(&iwork, 1, MPI_INTEGER, iwork, iwtagCor,
		   MPI_COMM_WORLD);
	  MPI_Send(percorso, 132, MPI_CHAR, iwork, iwtagCor,
		   MPI_COMM_WORLD);
	  

	}
      mdPrintf(ALL, " Finito di distribuire i processi\n", NULL);
      
      /*
	attende risultati da un processo generico
	until work requests have been exhausted.
      */
      
      for( iwork = islave + 1; iwork <= NUM_WORK_REQUEST; iwork++)
        {
	  //printf("ricevo\n");
	  ricevi(&iriceve, &iwork_recv);

	  ndone = ndone + 1;
	  
	  strcpy(filewalter, filename[iwork_recv]);
	  
	  writeconf(filewalter);
	  writeEnergy(filewalter, tempe[iwork_recv], steps[iwork_recv],
		      EpotIni, Emin, VACIni, VACmin); /*Questa ruotine va riscritta*/

	  /* now send new iwork to the process which had completed 
	     the job */

	  sscanf(filename[iwork], "%[^_]_R%d", stri, &rango[iwork]);
	  sprintf(TXT, "filename[iwork]: %s rango: %d\n", filename[iwork], rango[iwork]);
	  mdPrintf(ALL, TXT, NULL);
	  strcpy(filewalter, percorso);
	  strcat(filewalter, filename[iwork]);
	  strcpy(fileslave[iwork], filewalter);

	  readRestart(filewalter);            
	  
	  tempe[iwork] = Oparams.T / Oparams.lambda0[Oparams.lambdat[rango[iwork]]];
	  steps[iwork] = Oparams.curStep;
	  strcpy(fileslave[iwork], filewalter);
	  sprintf(logFile[iwork],"%sQncT%.6G.%d_R%dIW%d.log", 
		  percorso, 
		  Oparams.T/Oparams.lambda0[Oparams.lambdat[rango[iwork]]], 
		  Oparams.curStep, rango[iwork], iwork);
	  
	  sprintf(TXT, "logFile[%d]: %s\n", iwork, logFile[iwork]);
	  mdPrintf(ALL, TXT, NULL);
	  sendHeader(iriceve);
	  sendCoord(iriceve, Oparams.parnum);
	  
	  MPI_Send(&rango[iwork], 1, MPI_INTEGER, iriceve, iwtagCor,
		   MPI_COMM_WORLD);
	  MPI_Send(&iwork, 1, MPI_INTEGER, iriceve, iwtagCor,
		   MPI_COMM_WORLD);
	  MPI_Send(percorso, 132, MPI_CHAR, iriceve, iwtagCor,
		   MPI_COMM_WORLD);

	  sprintf(TXT, "Sent work n.%d(File:%s) to %d\n", iwork, filewalter, iriceve);
	  mdPrintf(ALL, TXT, NULL);
	}
      mdPrintf(ALL, "NOW ASPETTO I JOB MISSING\n", NULL);
      /*
	   wait for missing jobs
      */

      for(missing = ndone + 1; missing <= NUM_WORK_REQUEST; missing++)
	{

	  ricevi(&iriceve, &iwork_recv); 

	  strcpy(filewalter, filename[iwork_recv]);
	  writeconf(filewalter);

	  writeEnergy(filewalter, tempe[iwork_recv], steps[iwork_recv],
		      EpotIni, Emin, VACIni, VACmin); /*Questa ruotine va riscritta*/

	  sprintf(TXT, " Ricevo il missing %d da %d\n", missing, iriceve);
	  mdPrintf(ALL, TXT, NULL);
	}

      /*
	Tell all the slaves to exit.
      */
      sprintf(TXT, " DICO AGLI SLAVE DI USCIRE\n");
      mdPrintf(ALL, TXT, NULL);
      for(npr = 1; npr <= islave; npr++)
	{
	  sprintf(TXT, " Send 0 to %d\n",npr);
	  mdPrintf(ALL, TXT, NULL);
	  Oparams.curStep = 0;
	  sendHeader(npr);
	}
      
      mdPrintf(ALL, " MASTER HAS DONE EVERYTHING\n", NULL);
    }
  /*
    end of master
    
    -------------------------------------------------------------------------
    
    Now Slaves !
    Each slave process accepts work requests and returns
    results until a special termination request is received.
  */

  if (my_rank != 0) 
    {     
      for(njob = 1; njob <= njobmax; njob++)
	{

	  iriceve = recvHeader();

	  if (Oparams.curStep == 0)
	    break;
	  
	  initAll(Oparams.parnum);
	  
	  recvCoord(0);// riceve le coordinate dal processo master
	  buildMesh(OprogStatus.kmax, OprogStatus.Dk);
	  //printf("NNC:%d  kmax; %f Dk:%f\n", NNC, OprogStatus.kmax, OprogStatus.Dk);
	  
	  L = cbrt(Vol);
	  invL = 1.0 / L;

	  MPI_Recv(&rangsl, 1,  MPI_INTEGER, 
		   0, iwtagCor, MPI_COMM_WORLD, &status);
	  MPI_Recv(&numwrk, 1,  MPI_INTEGER, 
		   0, iwtagCor, MPI_COMM_WORLD, &status);
	  MPI_Recv(percorso, 132,  MPI_CHAR, 
		   0, iwtagCor, MPI_COMM_WORLD, &status);
	  
	  //printf("Ricevuti dati da: %d\n", iriceve);
	  
	  //INITALCHARGES(); //per i SILICA
	  Forces(Oparams.parnum);
	  
	  sprintf(TXT," EPOTEN INIZIALE =     %.5f | %.5f\n",
		  Epoten, Epoten/((double)Oparams.parnum));
	  mdPrintf(ALL, TXT, NULL);
	  
	  EpotIni = Epoten;
	  VACIni = VAC;

	  /*
	    minimization configuration
	  */
	  sprintf(TXT, "Using tol= 1E-12\n");
	  mdPrintf(ALL, TXT, NULL);
	  Ftol = 1.0E-12;
	  frprmn(Oparams.parnum, Ftol, &Fret);
	  Forces(Oparams.parnum);
	  Emin = Epoten;
	  VACmin = VAC;
	  sprintf(TXT, "DONE  WORK N.%d\n", numwrk);
	  mdPrintf(ALL, TXT, NULL);
	  
	  /* 0 = root process */

	    
          quanti = -1;
	  MPI_Send( &quanti, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	  MPI_Send( &EpotIni, 1, MPI_DOUBLE_PRECISION, 0, 0,
		    MPI_COMM_WORLD);
	  MPI_Send( &Emin, 1, MPI_DOUBLE_PRECISION, 0, 0,
		    MPI_COMM_WORLD);
	  MPI_Send( &VACIni, 1, MPI_DOUBLE_PRECISION, 0, 0,
		    MPI_COMM_WORLD);
	  MPI_Send( &VACmin, 1, MPI_DOUBLE_PRECISION, 0, 0,
		    MPI_COMM_WORLD);
	  MPI_Send( &numwrk, 1, MPI_INTEGER, 0, 0,
		    MPI_COMM_WORLD);
	  
	  sendCoord(0, Oparams.parnum);

	  // manda le coordinate al master
	  
	  /* NOTA: quando lo slava manda il messaggio di fine lavoro (quanti=-1)
	     non puo' piu' mandare messaggi di testo perchè a quel punto il master
	     manda un nuovo lavoro e non accetta del testo da loggare */
	  
	}
    }
  
  MPI_Finalize();
  return 0;
}
  
/* -----------------------------------------------------------------------------------
   NOTA BENE: pure questo dipende dalla simulazione */
  

/* ======================== >>> InitAl() <<< =========================== */
void initAll(int Nm)
{
  int size;
  char fina[512];

  static int ft = 1;
   
  nebrNow = 1;
  NNC = 0;
  pi = 2.0*acos(0.0);

  if (ft)
    ft = 0;
  else
    return;

  size = Nm * sizeof(double);
  AllocCoord(size, &xicomx, &xicomy, &xicomz,
	     &pcomx, &pcomy, &pcomz, 
	     &xix, &xiy, &xiz, &Gx, &Gy, &Gz, &Hx, &Hy, &Hz, NULL);
  AllocCoord(size, ALLOC_LIST, NULL);
  
  kcx = malloc(NK*sizeof(int));
  kcy = malloc(NK*sizeof(int));
  kcz = malloc(NK*sizeof(int));
  
  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac * Nm;
  nebrTab = AllocMatI(2, nebrTabMax); 
}

/* ------------------------------------------------------------------------------- */


/* ========================= >>> reaRestart <<< =========================== */
void readOne(char* fn)
{
  FILE* fs; 
  static int firstTime = 1;
  
  if ((fs = fopen(fn, "r")) == NULL)
    {
      sprintf(msgStrA, "Problem opening restart file %s ", fn);
      mdMsg(ALL, NOSYS, "ReadBakAscii", "ERROR", NULL,
	    msgStrA,
	    NULL);
    }
  
  //printf("fin qui ok: %s %p %p %p\n", fn, opro_ascii, opar_ascii, fs);
  readAsciiPars(fs, opro_ascii);
  readAsciiPars(fs, opar_ascii);

  /* Entrambe queste macro sono definite nel file mono_DPT.h */
  /* read up to coordinates begin */

  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum;
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  /* PATCH: solo la prima volta deve allocare le coordinate e non 
     ogni volta!!! */
  if (firstTime)
    {
      firstTime = 0;
      AllocCoord(SEGSIZE, ALLOC_LIST, NULL);
    }
  readAllCor(fs);

  fclose(fs);
}



/* ========================== >>> readStart <<< ============================ */
void readRestart(char filewalter[132])
{
  //equivalence (nrecord, recname(1));
  int i;
  
  /*readBakAscii(filewalter);*/
  readOne(filewalter); 
  /*
    check of periodic boundary conditions on the read data
  */
  L = cbrt(Vol);
  invL = 1.0/L;

  for (i=0; i < Oparams.parnum; i++)
    {
      
      /*  
	  if (abs(anperx)+abs(anpery)+abs(anperz).ne.0) then
	  print *, 'Fuori scatola in restart',i,nperx,npery,nperz
	  end if*/
      /* Scala al first box */
    
      rx[i] -= L*rint(invL*rx[i]);
      ry[i] -= L*rint(invL*ry[i]);
      rz[i] -= L*rint(invL*rz[i]);
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
  sprintf(TXT, "Scritto il file quenchato:%s\n", fina2);
  mdPrintf(ALL, TXT, NULL);
  saveBakAscii(fina2);
}

/* ========================= >>> writeEnergy <<< ========================= */
void writeEnergy(char filewalter[132], double T, int steps,
		 double VpotIni, double VpotEnd, double VACIni, double VACmin)
{
  char fina[132], fina2[512];
  FILE* ndat;

  strcpy(fina, filewalter);
  fina[0]='E';
  fina[1]='n';
  fina[2]='e'; 
  
  strcpy(fina2, percorso);
  strcat(fina2, fina);

  ndat = fopen(fina2, "w");
  fprintf(ndat,"%6.5G %d %.10G %.10G %.10G %.10G\n", 
	  T, steps, VpotIni, VpotEnd, VACIni, VACmin);
  sprintf(TXT, "Scritto il file Ene:%s\n", fina2);
  mdPrintf(ALL, TXT, NULL);
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
  //printf("fine\n");
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
  sprintf(TXT, "Brent exceed maximum iterations.\n");
  mdPrintf(ALL, TXT, NULL);
tre:     
  *xmin = X;
  return FX;
}

/* ============================ >>> F1dim <<< ========================== */
double F1dim(double X)
{
  //const int Nmax=1000;
  int j, i;
  for( j = 0; j < Oparams.parnum; j++)
    {
      rx[j] = pcomx[j] + X*xicomx[j];
      ry[j] = pcomy[j] + X*xicomy[j];
      rz[j] = pcomz[j] + X*xicomz[j];
      //printf("X:%f xicomx:%f pcomx:%f\n", X, xicomx[j], pcomx[j]);
    }
  /*
    check that all molecules are in the original box
    check of periodic boundary conditions on the read data
  */
  
  /* porta tutte le particelle nel primo box */
  for (i = 0; i < Oparams.parnum; i++)
    {
      rx[i] = rx[i] - L * rint(invL * rx[i]);
      ry[i] = ry[i] - L * rint(invL * ry[i]);
      rz[i] = rz[i] - L * rint(invL * rz[i]);
    }   

  //ENERGY
  Forces(Oparams.parnum);
  
  return Epoten;
  //  printf("x: %.6f f1dim: %.6f ncom:%d\n",x,f1dim,ncom)
}


/* ======================= >>> SDchkRebuild <<< ========================= */
int SDchkRebuild(void)
{
  int i;
  COORD_TYPE vv, vvMax = 0.0;

  for (i = 0; i < Oparams.parnum; i++)
    {
      vv = Sqr(xix[i]) + Sqr(xiy[i]) + Sqr(xiz[i]);
      if (vv > vvMax) 
	vvMax = vv;
    }
    
  dispHi = dispHi + sqrt(vvMax);
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;

}

/* ========================= >>> linmin <<< ======== ======================= */
void LinMin(int Nm, double *Fret)
{
  int j, i;
  //const nmax = 1000;
  const double tol = 1.0E-10;
  double AX, XX, BX, xmin, FX, FB, FA;
  //int Ncom;
 
  //Ncom = N;

  for (j = 0; j < Nm; j++) 
    {
      pcomx[j]  = rx[j];
      pcomy[j]  = ry[j];
      pcomz[j]  = rz[j];
      xicomx[j] = xix[j];
      xicomy[j] = xiy[j];
      xicomz[j] = xiz[j];
    }


  Forces(Nm);
  /*
    sprintf(TXT, "Prima di MNBRAK=%.6f\n", Epoten); 
    mdPrintf(ALL, TXT, NULL);*/

  //printf("bx et al',box,boy,boz,boxhx,boxhy,boxhz
  AX = 0.0;
  XX = 0.001 / fnorm;
  BX = 0.002 / fnorm;
  
  mnbrak(&AX, &XX, &BX, &FA, &FX, &FB, F1dim);
  *Fret = brent(AX, XX, BX, F1dim, tol, &xmin);
  /* sprintf(TXT, "IN LINMIN: MIN: %.6f\n", xmin);
     mdPrintf(ALL, TXT, NULL); */
  
  /*
    move the system
  */
  for(j = 0; j < Nm; j++)
    {
      xix[j] = xmin * xix[j];
      xiy[j] = xmin * xiy[j];
      xiz[j] = xmin * xiz[j];
      //if ((xmin == 0.0) && (xi[j] != 0.0)
      //  ) printf *,j,k,xi(j,k)
      rx[j] = pcomx[j] + xix[j];
      ry[j] = pcomy[j] + xiy[j];
      rz[j] = pcomz[j] + xiz[j];
    }
  /* porta tutte le particelle nel primo box */
  for (i = 1; i < Nm; i++)
    {
      rx[i]= rx[i] - L * rint(invL * rx[i]);
      ry[i]= ry[i] - L * rint(invL * ry[i]);
      rz[i]= rz[i] - L * rint(invL * rz[i]);
    }   

  /* sprintf(TXT, " XMIN= %.6f EN=%.6f\n", xmin, *Fret);
     mdPrintf(ALL, TXT, NULL);
  */
  SDchkRebuild();
  Forces(Nm);
  /*sprintf(TXT, "Dopo di MNBRAK=%.6f\n", Epoten); 
    mdPrintf(ALL, TXT, NULL);*/
  if (xmin == 0.0) 
    {
      sprintf(TXT, "ATTENZIONE X MIN= 0\n");
      mdPrintf(ALL, TXT, NULL);
    }

}

/* ======================= >>> frprmn <<< =========================== */
void frprmn(int Nm, double Ftol, double *Fret)
{
  const int ITMAX = 20000;
  const double eps = 1E-14;
  double FP, GG, DGG, GAM;
  int j, iter;
  
  //BuildNebrListNoLinked(Nm, Oparams.rcut, Oparams.sigma);
  Forces(Nm);
  FP = Epoten;
  /*
    calculate force modulus
  */
  fnorm = 0.0;
  for (j = 0; j < Nm; j++)
    {
      fnorm = fnorm + Fx[j]*Fx[j]+Fy[j]*Fy[j]+Fz[j]*Fz[j];
    }
  
  fnorm = sqrt(fnorm);
  
  //write(10,*) 0,fp,fnorm
  
  
  for(j = 0; j < Nm; j++)
    {
      Gx[j]  = Fx[j];
      Hx[j]  = Gx[j];
      xix[j] = Hx[j];
      Gy[j]  = Fy[j];
      Hy[j]  = Gy[j];
      xiy[j] = Hy[j];
      Gz[j]  = Fz[j];
      Hz[j]  = Gz[j];
      xiz[j] = Hz[j];
    }
  
  
  for(iter = 0; iter < ITMAX; iter++)
    {
      LinMin(Nm, Fret);
      
      if (iter%ITER_ELAPSED == 0)
	{
	  sprintf(TXTA[0],"[#%d] ",iter);
	}
      Forces(Nm);
      if (iter%ITER_ELAPSED == 0)
	{
	  sprintf(TXTA[1], "E=%.15G\n", Epoten);
	  mdPrintf(ALL, TXTA[0], TXTA[1], NULL);
	}
       /*
	 calculate force modulus
       */
      fnorm=0.0;
      for(j = 0; j < Nm; j++)
	{
	  fnorm = fnorm + Fx[j]*Fx[j]+Fy[j]*Fy[j]+Fz[j]*Fz[j];
	}
      fnorm = sqrt(fnorm);
      // SISTEMARE !!!!!!!!!!!!!!!!!!!!!!!!!!
      //write(10,*) iter,*Fret,fnorm;
    
      //print *,iter,*Fret,fnorm
      if (2.0 * fabs(*Fret-FP) <= Ftol*(fabs(*Fret)+fabs(FP) + eps))
	return;
      
      FP = Epoten;
      
      GG = 0.0;
      DGG = 0.0;
      for(j = 0; j < Oparams.parnum; j++)
	{
	  GG = GG + Sqr(Gx[j]) + Sqr(Gy[j]) + Sqr(Gz[j]);
	  //DGG=DGG+F(J,K)**2
	  DGG = DGG + (-Fx[j] + Gx[j])*(-Fx[j]) + (-Fy[j] + Gy[j])*(-Fy[j]) +
	    (-Fz[j] + Gz[j])*(-Fz[j]);
	}
      if (GG == 0)
	return;
      
      GAM = DGG / GG;
      for(j = 0; j < Nm; j++)
	{
	  Gx[j]  = Fx[j];
	  Gy[j]  = Fy[j];
	  Gz[j]  = Fz[j];
	  Hx[j]  = Gx[j] + GAM * Hx[j];
	  Hy[j]  = Gy[j]+ GAM * Hy[j];
	  Hz[j]  = Gz[j]+ GAM * Hz[j];
	  xix[j] = Hx[j];
	  xiy[j] = Hy[j];
	  xiz[j] = Hz[j];    
	  //         per steepest discent
	}
    }
  sprintf(TXT, "FRPR maximum iterations exceeded\n");
  mdPrintf(ALL, TXT, NULL);
}


/* ==================== >>> mnbrak <<< =================================== */
void mnbrak(double *AX, double *BX, double *CX, double *FA, double *FB, double *FC,
	     double (*func)(double))
{
  const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0E-15;
  double Q, FU, U, ULIM, DUM, R;

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
  
  //sprintf(TXT, "fine MNBRAK\n");
  //mdPrintf(ALL, TXT, NULL);
}

/* =========================== >>> forces <<< ======================= */
void  Forces(int Nm)
{
  int i;
  double L, invL;
  int old_rank;

  L = cbrt(Vol);
  invL = 1.0/L;
  
  /* porta tutte le particelle nel primo box */
  for (i=0; i < Nm; i++)
    {
      rx[i]=rx[i] - L*rint(invL*rx[i]);
      ry[i]=ry[i] - L*rint(invL*ry[i]);
      rz[i]=rz[i] - L*rint(invL*rz[i]);
    }   

  /* azzera tutte le forze */
  zeroArrays(Fx, Fy, Fz, Oparams.parnum);

  /* Costruisce le neighbour list */
  if (nebrNow)
  {
    nebrNow = 0;
    dispHi = 0.0;
    BuildNebrListNoLinked(Nm, Oparams.rcut, Oparams.sigma);
  }
  //printf("kmax: %f Dk:%f alpha:%f S0:%f", OprogStatus.kmax, OprogStatus.Dk, OprogStatus.alpha, OprogStatus.S0 );

  /* questo perchè nella simulazione di MD my_rank serve per identificare il processo
     ma qui my_rank è li rank del processo slave e quello che ci serve invece 
     è rangsl che li master invia agli slave */


  //printf("rangsl:%d\n", rangsl);
  ACForce(Nm, OprogStatus.kmax, OprogStatus.Dk, OprogStatus.alpha, OprogStatus.S0,
	  1.0);
  LJForce(Nm, Oparams.epsilon, Oparams.sigma, Oparams.rcut, 
	  1.0);
  /* Notare che le funzioni vengono chiamate ora con Lambda = 1.0 */

  /* Energia potenziale totale */
  Epoten = Vc + VAC;
  //printf("VAC:%f Vc:%f\n", VAC, Vc);
  //sprintf(TXT,"Epoten IN Forces: %.6f\n",Epoten);
  //mdPrintf(ALL, TXT, NULL);
}

