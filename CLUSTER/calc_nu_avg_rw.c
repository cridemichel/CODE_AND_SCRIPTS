#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define maxlen 1000
#define maxlenpn 100000
double avgcls[maxlen], count[maxlen];
double PN[maxlenpn];
int Nthr, N1, N2;
double norm, weight;
char fn[256], fake[256];
void main(int argc, char** argv)
{
  /* NOTA: calc_nu_avg_rw <listafiles> <file con PN rw> <threshold 1/2-1/2> */  
  FILE *f, *f2;
  int i, len, eml=0;
  double val, val0;
  f = fopen(argv[2], "r");
  
  Nthr = atof(argv[3]); 
  for (i=0; i < maxlenpn; i++)
    PN[i] = 0.0;
  while (!feof(f))
    {
      fscanf(f, "%lf %lf", &val0, &val);
      len=(int) val0;
      if (len < maxlenpn-1)
	PN[len] = val;
      //printf("len: %d val:%f PN=%f\n", len, val0, val);
    }
  fclose(f);
  for (i=Nthr; i < maxlenpn; i++)
    norm += PN[i];

  for (i=Nthr; i < maxlenpn; i++)
    PN[i] /= norm;
#if 0
  f =fopen("boh.dat","w");
  for (i=Nthr; i < maxlenpn; i++)
    fprintf(f, "%d %.15G\n", i, PN[i]);
  fclose(f);
#endif
  for (i=0; i < maxlen; i++)
    {	
      avgcls[i]=count[i]=0.0;
    }
  f = fopen(argv[1],"r");
  norm = 0.0;

  while (!feof(f)) 
   {
      fscanf(f,"%[^\n]\n", fn);
      sscanf(fn,"N_%d_%d/%[^\n]\n", &N1, &N2, fake);
      //printf("N1=%d N2=%d\n", N1, N2);
      if (N1 < Nthr) 
	continue;
      weight = (PN[N1]+PN[(N1+N2)/2]+PN[N2])/3.0;
      //weight = PN[(N1+N2)/2];
      f2=fopen (fn, "r");
      while (!feof(f2))
	{
          fscanf(f2,"%d %lf\n", &len, &val);
	  if (len < maxlen-1)
            {
              if (len > eml)
		eml=len;
              avgcls[len]+=val*weight;
              count[len]++; 
            }
         }
      fclose(f2);
   }


  fclose(f);
  for (i=0; i <= eml; i++)
    if (count[i]!=0)
       printf("%i %.15G\n", i, avgcls[i]);
}
