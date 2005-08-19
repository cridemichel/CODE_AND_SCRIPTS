/*
 * =====================================================================================
 * 
 *        Filename:  fitsq.c
 * 
 *     Description:  trova rozzamente il massimo con un fit quadratico. 
 * 
 *         Version:  1.0
 *         Created:  14/08/2005 12:54:34 CEST
 *        Revision:  none
 *        Compiler:  gcc
 * 
 *          Author:   (), 
 *         Company:  
 * 
 * =====================================================================================
 */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
char fname[1024];
double *func[2];
char line[4096], dummy[4096];
int main(int argc, char **argv)
{
  FILE *f;
  int maxxi, cc, pnts, first, ccini;
  double maxy, maxx;
  if (argc==1)
    {
      printf("devi fornire il nome del file come argomento\n");
      exit(-1);
    }
  if (argc==3)
    ccini = atoi(argv[2]);
  else
    ccini = 15;
  f = fopen(argv[1], "r");
  cc=0;
  while (!feof(f))
    {
      fscanf(f, "%[^\n]\n", line);
      cc++;
    }
  pnts = cc;
  func[0] = malloc(sizeof(double)*pnts);
  func[1] = malloc(sizeof(double)*pnts);
  rewind(f);
  for (cc = 0; cc < pnts; cc++)
    {
      fscanf(f, "%lf %lf%[^\n]\n", &func[0][cc], &func[1][cc], dummy);
      //printf("x: %.15G y: %.15G\n", func[0][cc], func[1][cc]);
    }
  fclose(f);
  first = 1;
  for (cc = ccini; cc < pnts; cc++)
    {
      //printf("maxx: %.15G maxy: %.15G x: %.15G y: %.15G\n", maxx, maxy, func[0][cc], func[1][cc]);
      if (first || func[1][cc] > maxy)
	{
	  first = 0;
	  maxy = func[1][cc];
	  maxx = func[0][cc];
	}
    }
  maxxi = rint(maxx);
  printf("%003d\n", maxxi);
  return 0;
}

