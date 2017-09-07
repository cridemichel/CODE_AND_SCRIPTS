#include<stdlib.h>
#include<stdio.h>
#define nmaxfiles 10000

char nomifiles[nmaxfiles][132];
int nacc[nmaxfiles];
double acc[nmaxfiles];
int tempi[nmaxfiles];
int main(int argc, char** argv)
{
  FILE *inpf, *of;
  int i, ind, aggiungi, nn, in, nfiles, maxt;
  double fpn2, fpn, T, VpotIni, VACIni, VACmin;

  inpf = fopen(argv[1],"r");
  if (inpf == NULL)
    {
      printf("errore nell'apertura del file %s\n", argv[1]);
      exit(-1);
    }

  ind = 0;
  while (1)
    {
      if (fscanf(inpf, "%s\n", nomifiles[ind]) < 1)
	break;
      ind++;
    }
  fclose(inpf); 
  nfiles = ind;
  maxt = 0;
  for (ind = 0; ind < nfiles; ind++)
    {
      inpf = fopen(nomifiles[ind], "r");
      if (inpf == NULL)
	{
	  printf("[%d] errore nell'apertura del file %s\n", ind, nomifiles[ind]);
	  exit(-1);
	}
      fprintf(stderr, "Leggo dal file %s\n", nomifiles[ind]);
      while (1)
	{
	  if (fscanf(inpf, "%lf %lf\n", &fpn2, &fpn) < 2)
	    break;
	  aggiungi = 1;
	  nn = maxt;
	  for (i = 0; i < maxt; i++)
	    {
	      if (in == tempi[i])
		{
		  aggiungi = 0;
		  nn = i;
		  break;
		}
	    }
	    
	  if (aggiungi)
	    {
	      tempi[maxt] = in;
	      fprintf (stderr, "aggiungo tempo %d [%d]\n", tempi[maxt], maxt);
	      acc[maxt] = fpn;
	      nacc[maxt] = 1;
	      maxt++;
	    }
	  else
	    {
	      acc[nn] += fpn;
	      nacc[nn]++;
	    }
	}
      fclose(inpf);
    }

  for (ind = 0; ind < maxt; ind++)
    {
      fprintf(stdout, "%d %.15G %d\n", tempi[ind], acc[ind]/((double)nacc[ind]), 
	      nacc[ind]);
    }
}

