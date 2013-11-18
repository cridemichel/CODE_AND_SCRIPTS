#define BASINS
#define MAIN
#include <mdsimul.h>
#include "Fsk.h"
#define ITER_ELAPSED 10
#define NJOBMAX 1000
#define NMOLA 800
#define KMODMAX 100 /* modulo (intero) massimo per i vettori d'onda */
#define NKSHELL 150 /* numero di vettori d'onda con un modulo fissato */
double tempei, dt, twopi;
char filewalter[512], filename[NJOBMAX][132], pwdfilename[NJOBMAX][512];
const double k0 = 6.12;
int it, tempi[NJOBMAX];
double tempe;
double rrx[NJOBMAX][NMOLA], rry[NJOBMAX][NMOLA], rrz[NJOBMAX][NMOLA];
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
 int ntripl[]=
#include "./ntripl.dat"
/* tempo di riferimento per il set di configurazioni utilizzate */
const int t0 = 1;
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double sqRe[NJOBMAX][NJOBMAX], sqIm[NJOBMAX][NJOBMAX];

/*============================== >>> writePEs <<< =========================== */
void writeFsk(int nf)
{
  /* Scrive le PE su stdout */
  int it;
  FILE *ofi;
  /* scrive la Fs(k) su stdout */
  fprintf(stdout, "%.8f\n", dt);
  for (it = 0; it < nf; it++)
    {
      if ((sqRe[0][it]!=0.0 || sqIm[0][it]!=0.0) && (tempi[it]>=tempi[0]))
	fprintf(stdout, "%d %d %.8f %.8f\n", 0, tempi[it]-tempi[0], sqRe[0][it], 
		sqIm[0][it]);
    }
}

void calcFsk(int kmod, int *Nm, int nf)
{
  int nptime, nacc, i, iq, kk, jj, tt;
  double sumRe, sumIm, rxdummy;

  for (iq=0; iq < ntripl[kmod]; iq++)
    {
      qx[kmod][iq]=invL*twopi*mesh[kmod][iq][0];
      qy[kmod][iq]=invL*twopi*mesh[kmod][iq][1];
      qz[kmod][iq]=invL*twopi*mesh[kmod][iq][2];
    }
  
  /* Inizia correlazioni */
  
  nptime=nf;
  nacc=0;
  
  /* cleaning */
  for (kk=0; kk < nptime; kk++)
   for (jj = 0; jj < nptime; jj++)
     {
       sqRe[kk][jj]=0.0;
       sqIm[kk][jj]=0.0;
     }
         
  /* iq: media su tutti i k che hanno modulo kmod e che sono stati scelti 
   * come mesh */
  for(iq=1; iq < ntripl[kmod]; iq++)
    {
      for (jj=0; jj < nptime; jj++)
       	{
   	  sumRe=0.0;
	  sumIm=0.0;

	  for (i=0; i < Nm[0]; i++)
	    {
      	      rxdummy = (rrx[0][i]-rrx[jj][i])*mesh[kmod][iq][0]
		+(rry[0][i]-rry[jj][i])*mesh[kmod][iq][1]
		+(rrz[0][i]-rrz[jj][i])*mesh[kmod][iq][2];
	      sumRe += cos(rxdummy);
	      sumIm += sin(rxdummy);	
	    }
	  
   	  sqRe[0][jj] += sumRe/((double) Nm[0]);
	  sqIm[0][jj] += sumIm/((double) Nm[0]);
    	}
      nacc=nacc+1;
    }

 /*
    Divide per il numero di vettori d'onda
    */
 tt=nacc;
 for (kk=0; kk < nptime; kk++)
   for (jj=0; jj < nptime; jj++)
     {
       sqRe[kk][jj] = sqRe[kk][jj]/((double)tt);
       sqIm[kk][jj] = sqIm[kk][jj]/((double)tt);
     }
}

/* ============================= >>> main <<< ============================== */
int main(int argc, char** argv)
{
  FILE* finp;
  FILE* pf;
  char listaFile[1024], percorso[512];
  int iwork, i, njob;

  pi = 2.0*acos(0.0);
  twopi = 2.0*pi;

  if (argc == 2)
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
 
  
  njob = 0;
      
  /* NOTA:
     - la prima riga deve contenere il  percorso dei file da analizzare
     - le altre righe devono contenere i nomi dei file senza il percorso */
  fscanf(finp, "%s", percorso);
  fprintf(stderr, "Leggo i file di restart dalla directory:%s\n", percorso);
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
      strcpy(filewalter, percorso);
      strcat(filewalter, filename[iwork]);
      strcpy(pwdfilename[iwork], filewalter);
      fprintf(stderr,"iwork:%d filewalter: %s\n", iwork, filewalter);
      readRestart(filewalter); 
      tempi[iwork] = Oparams.curStep;
      tempe = Oparams.T;
      dt = Oparams.steplength;
      for (i = 0; i < Oparams.parnum[0]; i++)
	{
	  rrx[iwork][i] = rx[0][i];
	  rry[iwork][i] = ry[0][i];
	  rrz[iwork][i] = rz[0][i];
	}
      fprintf(stderr, "tempe: %f dt: %f\n", tempe, dt);
      /* =============== *** E qui si puo' fare cio' che si vuole *** ===============*/
  
    }
    
  /* il primo argomento è kmod e deve essere il modulo del vettore per cui
   * la S(k) è massima!!! */
  calcFsk(18, Oparams.parnum, njob);
  
  /* ============================================================================*/
  writeFsk(njob);
  fprintf(stderr, " Finito calcolo Fs(k)\n");
  return 0;
}

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

  /* PATCH: solo la prima volta deve allocare le coordinate e non 
     ogni volta!!! */
  if (firstTime)
    {
      firstTime = 0;
      /* allocate  coordinates needed by simulation (see. mdsimul_p) */
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[0];
      AllocCoord(SEGSIZE, ALLOC_LISTA, NULL);
      /* allocate  coordinates needed by simulation (see. mdsimul_p) */
      SEGSIZE = sizeof(COORD_TYPE) * Oparams.parnum[1];
      AllocCoord(SEGSIZE, ALLOC_LISTB, NULL);
    }
  readAllCor(fs);

  fclose(fs);
}


/* ========================== >>> readStart <<< ============================ */
void readRestart(char filewalter[132])
{
  int i, a;
  
  readOne(filewalter);
  
  /*
    check of periodic boundary conditions on the read data
  */
  L = cbrt(Vol);
  invL = 1.0/Vol;

  for (a = 0; a < NA; a++)
    {
      for (i=0; i < Oparams.parnum[a]; i++)
	{
	  /* Scala al first box */
	  rx[0][i] -= L*rint(invL*rx[0][i]);
	  ry[0][i] -= L*rint(invL*ry[0][i]);
	  rz[0][i] -= L*rint(invL*rz[0][i]);
	}
    }
}

