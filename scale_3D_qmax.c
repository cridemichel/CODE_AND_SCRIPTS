#include <stdlib.h>
#include <stdio.h>
#include <string.h>
char X0[256], Phi[256];
char scal_arr_str_Phi[1000][256];
char scal_arr_str_X0[1000][256];
double scal_arr_dbl[1000];
double val;
int nf, type;
/* arg #1 file con gli scaling factors
 * arg #2 file da scalare (nel formato X0 Phi <valore>)*/
double cerca_X0(char *X0, char *Phi)
{
  int i;
  while (i < nf)
    {
      if (!strcmp(X0, scal_arr_str_X0[i]) && !strcmp(Phi, scal_arr_str_Phi[i]))
	return scal_arr_dbl[i];
      else
	{
	  fprintf(stderr,"ERROR: Scaling factor not found!\n");
	  exit(-1);
	}
      i++;
    }
}
int main(int argc, char **argv)
{
  FILE *f;
  int i;
  double sf;

  type = atoi(argv[1]);
  /* type = 0 moltiplica 
   * type = 1 dividi */
  f = fopen(argv[2], "r");
  i=0;
  while (!feof(f))
    {
      fscanf(f,"%s %s %lf\n", X0, Phi, &val);
      strcpy(X0, scal_arr_str_X0[i++]);
      strcpy(Phi, scal_arr_str_Phi[i]);
      scal_arr_dbl[i] = val;
    }
  nf = i;
  fclose(f);
  f = fopen(argv[3], "r");
  while (!feof(f))
   {
     fscanf(f,"%s %s %lf\n", X0, Phi, &val);
     sf = cerca_X0(X0, Phi);
     printf("sf=%f\n", sf);
     if (type==0)
       printf("%s %s %f\n", X0, Phi, val*sf); 
     else
       printf("%s %s %f\n", X0, Phi, val/sf); 
  }
  fclose(f);
}
