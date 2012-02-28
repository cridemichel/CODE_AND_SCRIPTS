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
  double accw, val, val0, xmed, am, ene, accene, cc;
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
  
  norm=0.0;
  for (i=Nthr; i < maxlenpn; i++)
    norm += PN[i];

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
  accene = accw= 0.0;
  cc=0;
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
      fscanf(f2, "%lf\n", &ene);
      //printf("ene=%.15G N1=%d , N2=%d\n", ene, N1, N2);
      accene += ene*weight;
      accw += weight;
      cc++;
      fclose(f2);
   }
  fclose(f);
  printf("Energia media=%.15G (cc=%f, accw=%.15G, Nthr=%d)\n", accene/accw, cc, accw, Nthr);
}
