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
#define npmax maxlen

const int nlin=50;
int l1[npmax], l2[npmax];
double dlog[npmax], xlog[npmax];
int media_log=0;

void main(int argc, char** argv)
{
  /* NOTA: calc_nu_avg_rw <listafiles> <file con PN rw> <threshold 1/2-1/2> */  
  FILE *f, *f2;
  int i, len, eml=0;
  double val, val0, xmed, am;
  int kmax, kj, i3, kk;
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
  if (media_log)
    {
      for (kk=1; kk <= 51; kk++)
	l1[kk]=(int) nlin*pow(1.25,kk-1);

      for(kk=1; kk <= 50; kk++)
	{
	  l2[kk]=l1[kk+1]-1;
	  if (l2[kk] < npmax) 
	    kmax=kk;
	  //printf("l1=%d l2=%d kmax=%d\n", l1[kk], l2[kk], kmax);	
	}
      for(kk=1; kk <= kmax; kk++)
	{
	  dlog[kk]=0.0;
	  xlog[kk]=0.0;
	  for (kj=l1[kk]; kj <= l2[kk]; kj++)
	    {
	      if (avgcls[kj] !=0 && kj < npmax)
		{   
		  dlog[kk]=dlog[kk]+avgcls[kj];
		}		  
	      xlog[kk]=xlog[kk]+kj;
	    }
	}
      for (i3=0; i3 < nlin; i3++)
	{
	  if (avgcls[i3] != 0) printf("%d %.15G\n", i3,(double)avgcls[i3]);
	}
      for (kk=1; kk <= kmax; kk++)
	{
	  if (dlog[kk]!=0) 
	    {
	      am=l2[kk]-l1[kk]+1; 
	      xmed=xlog[kk]/am;
	      dlog[kk]=dlog[kk]/am; 
	      //printf("xmed=%f am=%f\n", xmed, am);
	      printf("%.15G %.15G\n", xmed,((double)dlog[kk]));
	    }
	}
    }
  else
    {
      for (i=0; i <= eml; i++)
	if (count[i]!=0)
	  printf("%i %.15G\n", i, avgcls[i]);
    }
}
