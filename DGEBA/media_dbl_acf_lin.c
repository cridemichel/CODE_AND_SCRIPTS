#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define nmaxfiles 100000
#define ITMAX 1000
char nomifiles[nmaxfiles][132];
int nacc[nmaxfiles];
double acc[nmaxfiles];
double tempi[nmaxfiles];
int main(int argc, char** argv)
{
  FILE *inpf, *of;
  int i, ind, aggiungi, nn, nfiles, maxt, nm;
  double fpn, T, VpotIni, VACIni, VACmini, in, dt;

  inpf = fopen(argv[1],"r");
  if (inpf == NULL)
    {
      printf("errore nell'apertura del file %s\n", argv[1]);
      exit(-1);
    }

  ind = 0;
  if (argc >= 3)
     dt = atof(argv[2]);	  
  else 
     dt = 1.0;
  fprintf(stderr, "Using dt=%.15G\n", dt);
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
  double tt, tmin, tbin[ITMAX]={0};
  double avgacf[ITMAX]={0}, ccc[ITMAX]={0};

  tmin=0.0;
  tbin[0] = tmin;

  for (it = 1; it < ITMAX; it++)
    {
      tbin[it] = tbin[it-1] + dt;
    }
  for (ind=0; ind < maxt; ind++)
    {
      tt = tempi[ind];
      for (it=1; it < ITMAX; it++)
        {
          if (tt >= tbin[it-1] && tt < tbin[it])
            {
              avgacf[it-1] += acc[ind]/((double)nacc[ind]);
              ccc[it-1] += 1;  
	      break;
            }
        }
    }
  for (it = 0; it < ITMAX; it++)
    {
      if (ccc[it] != 0.0)
        fprintf(stdout, "%.15f %.15G %f\n", tbin[it]+dt*0.5, avgacf[it]/ccc[it], ccc[it]);
    }
}

