#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
#define KMODMAX 599
#define NKSHELL 150
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi, qavg, Sq, scalFact,  L;
int main(int argc, char** argv)
{
  FILE *f;
  char *fn;
  int mp, qmod;
#if 1
  if (argc == 1)
    {
      printf("fixSq <sq.dat> <box length>\n");
      exit(-1);
    }
#endif
  twopi = acos(0.0)*4.0;
  fn=argv[1];

  f=fopen(fn, "r");

  L=atof(argv[2]);
  scalFact=twopi/L; 
  while (!feof(f))
    {
      fscanf(f, "%d %lf\n", &qmod, &Sq);
      qavg = 0;
      for (mp = 0; mp < ntripl[qmod]; mp++) 
	{
	  qavg += sqrt(Sqr(((double)mesh[qmod][mp][0]))+
  		       Sqr(((double)mesh[qmod][mp][1]))+
		       Sqr(((double)mesh[qmod][mp][2])));
	}
      qavg *= scalFact/((double)ntripl[qmod]);
      fprintf(stdout, "%.15G %.15G\n", qavg, Sq); 
    }
  return 0;
}
