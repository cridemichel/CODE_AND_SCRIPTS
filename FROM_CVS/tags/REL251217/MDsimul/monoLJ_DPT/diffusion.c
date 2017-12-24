#define BASINS
#define MAIN
#include <mdsimul.h>
#include <diffusion.h>
#define ITER_ELAPSED 10
#define MAX_TEMPS 30
#define NJOBMAX 1000
#define PT_MSD
double tempe;
char filewalter[512], filename[NJOBMAX][132], pwdfilename[NJOBMAX][512];
double tempo_ur[NJOBMAX], MSD[NJOBMAX];
int firstTime[MAX_TEMPS], tempo_ini[MAX_TEMPS], 
  tempiacc[NJOBMAX],tempo[NJOBMAX], ntempi = 0;
#ifdef PT_MSD
double *rxI, *ryI, *rzI;
#else
double *rxI[MAX_TEMPS], *ryI[MAX_TEMPS], *rzI[MAX_TEMPS];
#endif
/*============================== >>> writePEs <<< =========================== */
void writeDiff(double Nm)
{
  /* Scrive MSD su stdout */
  int i;
  FILE *ofi;

  ofi = fopen("msd.dat","w");
  for (i = 0; i < ntempi; i++)
    {
      fprintf( ofi, "%.8f %.8f %.8f %d\n", tempo_ur[i], MSD[i]/(((double)tempiacc[i])*Nm),
	       MSD[i]/(6.0*tempo_ur[i]*Nm*((double)tempiacc[i])),
	    tempiacc[i] );
    }

  fclose(ofi);

}

void diffusione(int Nm, int curstp, double dt, int rango);

/* ============================= >>> main <<< ============================== */
void main(int argc, char** argv)
{
  FILE* finp;
  char listaFile[1024], stri[512], percorso[512];
  int iwork, i, iE, iT, njob;
  int rango;
  
  if (argc == 2)
    {
      // printf("arg: %s\n",argv[1]);
      strcpy(listaFile, argv[1]);
    }
  else
    {
      sprintf(TXTA[0],"Devi fornire come argomenti il nome del file contenente il percorso\n");
      sprintf(TXTA[1],"e la lista di files!\n");
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
  sprintf(TXT, "Leggo i file di restart dalla directory:%s\n", percorso);
  mdPrintf(ALL, TXT, NULL);

  for (i = 0; i < MAX_TEMPS; i++)
    firstTime[i] = 1;

  for (i=0; i < NJOBMAX; i++)
    {
      tempiacc[i] = 0;
      MSD[i] = 0.0;
    }
  ntempi = 0;
 
  while(!feof(finp))
    {
      if (fscanf(finp, "%s", filename[njob]) == 1)
	{
	  //printf("filename:%s\n", filename[njob]);
	  sscanf(filename[njob], "%[^_]_R%d", stri, &rango);
	  //sprintf(TXT, "filename[iwork]: %s rango: %d\n", filename[njob], rango);
	  //mdPrintf(ALL, TXT, NULL);
	  strcpy(filewalter, percorso);
	  strcat(filewalter, "/");
	  strcat(filewalter, filename[njob]);
	  strcpy(pwdfilename[njob], filewalter);
	  readRestart(pwdfilename[njob]);
	  if (Oparams.curStep == 1)
	    initAll(Oparams.parnum, rango, Oparams.curStep);
	  njob++;
	}
    }
  
  // njob--;

  sprintf(TXT, " Number of restart files=%d\n", njob);
  mdPrintf(ALL, TXT, NULL);

  ntempi = 0;
  printf("Work started:");
  for(iwork = 0; iwork < njob; iwork++)
    {
      //printf("iwork:%d filewalter: %s\n", iwork, filewalter);
      sscanf(filename[iwork], "%[^_]_R%d", stri, &rango);
      readRestart(pwdfilename[iwork]); 
      diffusione(Oparams.parnum,Oparams.curStep, Oparams.steplength, rango);
      tempe = Oparams.T / Oparams.lambda0[Oparams.lambdat[rango]];
      printf("#");
      fflush(NULL);
      /* ============================================================================*/
      
    }
  writeDiff( (double) Oparams.parnum);
  printf("\n Nm: %d Work finished ntempi: %d\n", Oparams.parnum, ntempi);
  
  mdPrintf(ALL, " Finito di distribuire i processi\n", NULL);
}

/* ======================== >>> InitAl() <<< =========================== */
void initAll(int Nm, int rango, int ti)
{
  int i;
  /* Questa routine va eseguita una sola volta!!!
     Si assume che tutti i file di restart abbiamo lo stesso
     numero di particelle e lo stesso nebrTabFac */
#ifdef PT_MSD
  rxI = malloc(Nm*sizeof(double));
  ryI = malloc(Nm*sizeof(double));
  rzI = malloc(Nm*sizeof(double));
  for (i = 0; i < Nm; i++)
    {
      rxI[i] = rx[i];
      ryI[i] = ry[i];
      rzI[i] = rz[i]; 
      //printf("(%f,%f,%f)\n", rxI[i], ryI[i], rzI[i]);
    }
#else
  if (firstTime[rango])
    {
      printf("rango: %d\n", rango);
      firstTime[rango] = 0;
    }
  else
    {
      printf("ERRORE: ci sono piu' file con lo stesso rango(%d) e step iniziale 1!!!\n",
	     rango);
      exit(1);
    }

  rxI[rango] = malloc(Nm*sizeof(double));
  ryI[rango] = malloc(Nm*sizeof(double));
  rzI[rango] = malloc(Nm*sizeof(double));
  tempo_ini[rango] = ti;

  for (i = 0; i < Nm; i++)
    {
      rxI[rango][i] = rx[i];
      ryI[rango][i] = ry[i];
      rzI[rango][i] = rz[i]; 
      //printf("(%f,%f,%f)\n", rxI[i], ryI[i], rzI[i]);
    }
#endif
  
}

/* ------------------------------------------------------------------------------- */

/* ========================= >>> diffusione <<< ========================== */
void diffusione(int Nm, int curstp, double dt, int rango)
{
  int i, it=0;
  
    
  for (i = 0; i < ntempi; i++)
    {
      if (curstp == tempo[i])
	{
	  /* il tempo esiste gia' quindi accumuliamo!*/
	  it = i;
	  break;
	}    
    }
  if (i == ntempi)
    {
      printf("nuovo tempo: %d\n", curstp);
      /* Si tratta di un nuovo tempo quindi memorizziamolo! */
      tempo[ntempi] = curstp;
      tempo_ur[ntempi] = ((double)curstp) * dt;
      it = ntempi;
      ntempi++;
    }
  ++tempiacc[it];
  MSD[it] = 0.0;
  for (i = 0; i < Nm; i++)
    {
#ifdef PT_MSD
      MSD[it] += Sqr(rx[i]-rxI[i]) + Sqr(ry[i]-ryI[i]) + 
	Sqr(rz[i]-rzI[i]);
#else
      MSD[it] += Sqr(rx[i]-rxI[rango][i]) + Sqr(ry[i]-ryI[rango][i]) + 
	Sqr(rz[i]-rzI[rango][i]);
#endif
    }
  printf("MSD[%d]:%.10f\n ", it, MSD[it]);
  //Diff[it] /= tempo_ur[it]*6.0*((double)Nm);

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
  
}
