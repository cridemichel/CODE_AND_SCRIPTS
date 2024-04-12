#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define nmaxfiles 10000
#define ITMAX 30
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
	  if (fscanf(inpf, "%lf %lf\n", &in, &fpn) < 2)
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

  int it;
  double ltmin, logt, logdt=0.1, logtbin[ITMAX]={0};
  double avgacf[ITMAX]={0}, ccc[ITMAX]={0};

  ltmin=-2.0;
  logtbin[0] = ltmin;

  for (it = 1; it < ITMAX; it++)
    {
      logtbin[it] = logtbin[it-1] + logdt;
    }
  for (ind=0; ind < maxt; ind++)
    {
      logt = log10(tempi[ind]);
      for (it=1; it < ITMAX; it++)
        {
          if (logt > logtbin[it-1] && logt < logtbin[it])
            {
              avgacf[it] += acc[ind]/((double)nacc[ind]);
              ccc[it] += 1;  
            }
        }
    }
  for (it = 0; it < ITMAX; it++)
    {
      if (ccc[it] != 0.0)
        fprintf(stdout, "%.15f %.15G %d\n", logtbin[it]+logdt*0.5, avgacf[it]/ccc[it], nacc[ind]);
    }
}

