#include<stdio.h>
#include<stdlib.h>
#define maxfiles 10000

char todof[maxfiles][132], qncarr[maxfiles][132], cnfarr[132];

void main(int argc, char** argv)
{
  FILE *tdf, *fcnf, *fqnc, *of;
  int aggiungi, i, ind, nqnc, ntodo;
  char tmpfn[132];

  if (argc < 3)
    {
      fprintf(stderr, "Devi fornire due file come argomenti(<Cnf> <Qnc>)!\n");
      exit(-1);
    }
  //fprintf(stderr, "opening %s\n", argv[2]);
  fqnc = fopen(argv[2], "r");
  
  ind = 0;
  
  /* legge il file contenente la lista (senza percorso) di tutti i files quanchati */
  while(!feof(fqnc))
    {
      if (fscanf(fqnc, "%s", qncarr[ind]) ==1)
	{
	  ind++;
	}
    }
  fprintf(stderr, "   File quenchati: %8d\n", ind);
  fclose(fqnc);
  nqnc = ind;
  /* confronta la lista dei file quenchati con quella dei file di configurazioni */
  ind = 0;
  ntodo = 0;
  fcnf = fopen(argv[1], "r");
  while (!feof(fcnf))
    {
      if (fscanf(fcnf, "%s", cnfarr) < 1)
	{
	  break;
	}
      strcpy(tmpfn, cnfarr);
      tmpfn[0] = 'Q';
      tmpfn[1] = 'n';
      tmpfn[2] = 'c';
      aggiungi = 1;
      for(i = 0; i < nqnc; i++ )
	{
	  if (!strcmp(tmpfn, qncarr[i]))
	      {
		aggiungi = 0; 
		break;
	      }
	}
      if (aggiungi)
	{
	  strcpy(todof[ntodo], cnfarr);
	  ntodo++;
	}
      ind++;
    }
  fclose(fcnf);
  fprintf(stderr, "File da quenchare: %8d\n", ntodo);
  /* scrive su stdout tutti i file ancora da quanchare */
  for ( i = 0; i < ntodo; i++ )
    {
      printf("%s\n", todof[i]);
    }
 }

