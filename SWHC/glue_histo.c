#include<math.h>
#include<stdlib.h>
#include<stdio.h>
const long double NmaxAlloc=10000;
char dn[4096], cm[4096];
void main(int argc, char **argv)
{
  FILE *f, *f1;
  long double *histo, pn, fact, sum;
  int i, N=-1, oldN=-1, Nmin, Nmax;
  f1=fopen(argv[1],"r");

  //printf("fn=%s\n", argv[1]);  
  histo = malloc(sizeof(long double)*NmaxAlloc);
  for (i=0; i < Nmax; i++)
    histo[i]=0;
  while (!feof(f1))
	{

  fscanf(f1, "%[^\n]\n", dn);
  //sprintf(cm, "cd %s", dn);
  //system(cm);
  //printf("cm=%s\n", cm);
  chdir(dn);
  f = fopen("histo.dat", "r");
  while (!feof(f))
    {
      fscanf(f, "%d %LG\n", &N, &pn);
      if (oldN==-1)
	{
	  fact=1.0;
	  histo[N]=pn;
	  Nmin=N;
	}
      else if (N==oldN)
	{
	  fact = histo[N]/pn; 
       //printf("N=%d pn=%f fact=%f\n", N, pn, fact);
	}
      else 
	{
	  histo[N] = pn*fact;
	}
      oldN = N;
    }   
    fclose(f);
  //system("cd ..");
   chdir("..");
  }
 fclose(f1);
 Nmax=N;
  sum=0;
  for (i=Nmin; i < Nmax; i++)
    {
      sum += histo[i];
    }
  for (i=Nmin; i < Nmax; i++)
    {
      histo[i]/=sum;
    }
 f = fopen("histor.dat", "w+");
  for (i=Nmin; i < Nmax; i++)
    fprintf(f, "%d %.15LG\n", i, histo[i]);
  fclose(f);
}
