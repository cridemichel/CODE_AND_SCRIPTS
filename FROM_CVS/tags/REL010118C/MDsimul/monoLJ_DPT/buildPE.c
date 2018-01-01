#define BASINS
#define MAIN
#include <mdsimul.h>
#include <buildPE.h>
#define ITER_ELAPSED 10
#define MAX_TEMPS 14
#define NJOBMAX 1000
double tempe;
char filewalter[512], filename[NJOBMAX][132], pwdfilename[NJOBMAX][512];
double ENmin, ENmax;
int maxlt=0, lambdat, rango, nPE, iPEN[MAX_TEMPS][PE_POINTS];
double Epoten;


/*============================== >>> writePEs <<< =========================== */
void writePEs(void)
{
  /* Scrive le PE su stdout */
  int i, iE;
  FILE *ofi;
  double EN;

  maxlt++;
  printf("nPE: %d maxlt: %d\n", nPE, maxlt);
  ofi = fopen("energie.dat","w");
  for (i = 0; i < maxlt; i++)
    {
      for (iE = 0; iE < nPE; iE++)
	{
	   
	  EN = ENmin + (((double) iE) / ((double)nPE)) * (ENmax - ENmin); 
	  fprintf(ofi, "%.8f %d\n",EN, iPEN[i][iE] );
	}
      
      fprintf(ofi, "&");
    
    }
  fclose(ofi);

}


/* ============================= >>> main <<< ============================== */
int main(int argc, char** argv)
{
  FILE* finp;
  FILE* pf;
  char listaFile[1024], stri[512], percorso[512];
  int iwork, i, iE, iT, njob;

  
  if (argc == 2)
    {
      if ((pf = fopen(argv[1],"r"))==NULL)
	{
	  sprintf(TXT, "Unable to open file: %s\n", argv[1]);
	  mdPrintf(ALL, TXT, NULL);
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
    }


  if ( (finp = fopen(listaFile, "r"))==NULL )
    {
      sprintf(TXT, "Unable to open file: %s\n", listaFile);
      mdPrintf(ALL, TXT, NULL);
      exit(-1);
    }
  
  njob = 0;
      
  /* NOTA:
     - la prima riga deve contenere l'energia minima, quella massima e il numero 
       di punti in cui suddividere tale intervallo di enegie(divise per N)
     - la seconda riga il percorso dei file da analizzare
     - le altre righe devono contenere i nomi dei file senza il percorso */
  fscanf(finp, "%s", percorso);
  fscanf(finp, "%lf %lf %d\n", &ENmin, &ENmax, &nPE);
  sprintf(TXT, "Leggo i file di restart dalla directory:%s\n", percorso);
  mdPrintf(ALL, TXT, NULL);
  while(!feof(finp))
    {
      if (fscanf(finp, "%s", filename[njob]) == 1)
	{
	  //printf("file:%s\n", filename[njob]);
	  njob++;
	}
    }
  
 // njob--;

  sprintf(TXT, " Number of restart files=%d\n", njob);
  mdPrintf(ALL, TXT, NULL);

  buildMesh(OprogStatus.kmax, OprogStatus.Dk);
  
  for (iT = 0; iT < MAX_TEMPS; iT++)
    {
      for (iE = 0; iE < nPE; iE++ )
	{
	  iPEN[iT][iE] = 0; 
	}
    }

  for(iwork = 0; iwork < njob; iwork++)
    {
      //printf("filename:%s\n", filename[njob]);
      sscanf(filename[iwork], "%[^_]_R%d", stri, &rango);
      sprintf(TXT, "filename[iwork]: %s rango: %d\n", filename[iwork], rango);
      mdPrintf(ALL, TXT, NULL);
      strcpy(filewalter, percorso);
      strcat(filewalter, filename[iwork]);
      strcpy(pwdfilename[iwork], filewalter);
      printf("iwork:%d filewalter: %s\n", iwork, filewalter);
	  
      readRestart(filewalter); 
      initAll(Oparams.parnum);

      tempe = Oparams.T / Oparams.lambda0[Oparams.lambdat[rango]];
      lambdat = Oparams.lambdat[rango];
      if (lambdat > maxlt)
	maxlt = lambdat;
      printf("tempe: %f lambdat: %d\n", tempe, lambdat);
      /* =============== *** E qui si puo' fare cio' che si vuole *** ===============*/
      
      Forces(Oparams.parnum); 
      updatePE(lambdat, Oparams.parnum);
      /* ============================================================================*/
    }
  
  writePEs();
  
  mdPrintf(ALL, " Finito di distribuire i processi\n", NULL);
}

/* ======================== >>> InitAl() <<< =========================== */
void initAll(int Nm)
{
  int size;
  char fina[512];
  static int firstTime = 1;

  /* Questa routine va eseguita una sola volta!!!
     Si assume che tutti i file di restart abbiamo lo stesso
     numero di particelle e lo stesso nebrTabFac */

  if (firstTime)
    firstTime = 0;
  else
    return;

  size = Nm * sizeof(double);
  AllocCoord(size, ALLOC_LIST, NULL);
  
  kcx = malloc(NK*sizeof(int));
  kcy = malloc(NK*sizeof(int));
  kcz = malloc(NK*sizeof(int));
  NNC = 0;
  pi = 2.0*acos(0.0);
  
  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac * Nm;
  nebrNow = 1;
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
  int i;
  
  readOne(filewalter);
  
  /*
    check of periodic boundary conditions on the read data
  */
  L = cbrt(Vol);
  invL = 1.0/Vol;

  for (i=0; i < Oparams.parnum; i++)
    {
      /* Scala al first box */
      rx[i] -= L*rint(invL*rx[i]);
      ry[i] -= L*rint(invL*ry[i]);
      rz[i] -= L*rint(invL*rz[i]);
    }
}


/* =========================== >>> forces <<< ======================= */
void  Forces(int Nm)
{
  int i;
  double L, invL;
  int old_rank;

  L = cbrt(Vol);
  invL = 1.0/L;
  
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
