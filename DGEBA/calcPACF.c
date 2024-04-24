#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXPTS 1000
char **fname; 
double time, *ACF, *tempi, *pointsArr, *cc, veltmp, *ti, L, refTime, *omACV, *velACV;
int points, assez, NP, NPA, npoints;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0, A1, B0, B1, C0, C1, ACFm, ACFSqmi, dummy1, dummy2;


int main(int argc, char **argv)
{
  FILE *f, *f2, *fA, *fB;
  int first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a;
  int skip, NN, fine, JJ, maxl, nfiles, nat, np, maxnp, npoints;
  if (argc <= 1)
    {
      printf("Usage: calcPACF <file_dati>\n");
      exit(-1);
    }
 
  c2 = 0;
  f2 = fopen(argv[1], "r");
  maxl = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", dummy); 
      c2++;
    }	
  npoints = c2;
  points=npoints;
  if (argc>=3)
    points = atoi(argv[2]);
  printf("npoints=%d\n", npoints);
  rewind(f2);
  tempi = malloc(sizeof(double)*npoints);
  pointsArr = malloc(sizeof(double)*npoints);

  skip=1;
  if (argc >=4)
     skip = atoi(argv[3]);

  for (ii=0; ii < npoints; ii++)
    {
      fscanf(f2, " %lf %lf %lf %lf ", &(tempi[ii]), &(pointsArr[ii]), &dummy1, &dummy2); 
      //printf("tempo=%.15G punto=%.15G\n", tempi[ii], pointsArr[ii]);
    }
  cc = malloc(sizeof(double)*npoints);
  ACF= malloc(sizeof(double)*npoints);
  ti = malloc(sizeof(double)*npoints);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;

  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    {
      ACF[ii] = 0.0;
      cc[ii] = 0.0;
    }
  c2 = 0;
  JJ = 0;
  NN=skip;
  for (nr1 = 0; nr1 < npoints; nr1=nr1+NN)
    {	
      fine = 0;
      for (nr2 = nr1; nr2 < npoints; nr2++)
	{
	  /* N.B. considera NN punti in maniera logaritmica e poi calcola i punti in maniera lineare 
	   * distanziati di NN punti. */
	  np = nr2-nr1;	      
	  if (nr2 >= npoints || np >= points)
	    {
	      fine = 1;
	      break;
	    }
	  if (np < points && ti[np] == -1.0)
	    {
	      ti[np] = tempi[np];
	      //printf("np=%d time=%.15G\n", np, ti[np]);
	    }
	  if (nr2 == nr1)
	    continue;
	  ACF[np] += pointsArr[nr2]*pointsArr[nr1];
	  cc[np]+=1.0;
	  //printf("points nr1=%f nr2=%f\n", pointsArr[nr1], pointsArr[nr2]);
	  //printf("cc[%d]=%f np=%d nr1=%d nr2=%d\n", np, cc[np], np, nr1, nr2);
	}
    }
  f = fopen("ACF.dat", "w+");
  for (ii=1; ii < points; ii++)
    {
      ACF[ii] = ACF[ii]/cc[ii];
      if (ti[ii] > -1.0)
	{
	  fprintf(f, "%.15G %.15G\n", ti[ii]-ti[0], ACF[ii]);
	  //printf("%.15G %.15G %.15G\n", ti[ii]-ti[0], ACF[ii], cc[ii]);
	}
    }
  fclose(f);
  return 0;
}
