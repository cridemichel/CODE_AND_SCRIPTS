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
  int maxxi, cc, pnts, ccini=0;
  double first, soglia;
  double maxy, maxx, m;
  if (argc < 3)
    {
      printf("devi fornire il nome del file come argomento e la soglia (reale compreso tra 0 e 1)\n");
      exit(-1);
    }
  soglia = atof(argv[2]);
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
  first = 1.0;
  for (cc = ccini; cc < pnts; cc++)
    {
      //printf("maxx: %.15G maxy: %.15G x: %.15G y: %.15G\n", maxx, maxy, func[0][cc], func[1][cc]);
      //printf("value=%.15G thr=%.15G\n", func[1][cc], soglia*first);
      if (func[1][cc] < soglia*first)
	{
	  if (cc==ccini)
	    {
	      fprintf(stderr,"first point is already below threashold!\n");
	      fprintf(stderr,"Chage it!\n");
	      exit(-1);
	    }
	  m = (func[1][cc] - func[1][cc-1])/(func[0][cc]-func[0][cc-1]);
	  if (m==0.0)
	    {
	      fprintf(stderr,"Error: m=0!\n");
	      exit(-1);
	    }
	  /* perform linear interpolation */
	  printf("%.15G", func[0][cc] + (soglia-func[1][cc])/m);
	  exit(-1);
	}
    }
  return 0;
}

