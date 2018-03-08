#include<stdlib.h>
#include<stdio.h>
#define nmaxfiles 10000

char nomifiles[nmaxfiles][132];
int nacc[nmaxfiles];
double acc[nmaxfiles];
double tempi[nmaxfiles];
int main(int argc, char** argv)
{
  FILE *inpf, *of;
  int i, ind, aggiungi, nn, nfiles, maxt, nm;
  double fpn, T, VpotIni, VACIni, VACmini, in;

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
	  if (fscanf(inpf, "%lf %lf %d\n", &in, &fpn, &nm) < 2)
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
	      fprintf (stderr, "aggiungo tempo %f [%d]\n", tempi[maxt], maxt);
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
      fprintf(stdout, "%.15f %.15G %d\n", tempi[ind], acc[ind]/((double)nacc[ind]), 
	      nacc[ind]);
    }
}

