#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXPTS 1000
char **fname; 
double time, *cc, veltmp, *ti, *vel0[3], *velt[3], 
       L[3], *velACV;
int points, assez, NP;
char parname[12800], parval[25600000], line[25600000];
char dummy[2048000];
double A1, B0, B1, C0, C1;
double dt=-1.0;
void readconf(char *fname, double *vel[3])
{
  FILE *f;
  double v[3];
  int i;
  f = fopen(fname, "r");
  //printf("fn=%s\n", fname);
  //cpos = ftell(f);
  //printf("cpos=%d\n", cpos);
  //fscanf(f, "%[^\n]\n",line);
  fscanf(f,"%d %lf %lf %lf\n", &NP, &L[0], &L[1], &L[2]);
  //printf("NP=%d\n", NP);
  for (i = 0; i < NP; i++) 
    fscanf(f, "%[^\n]\n",line);
  for (i = 0; i < NP; i++) 
    {
      fscanf(f, "%lf %lf %lf\n", &v[0], &v[1], &v[2]);
    }
  fclose(f);
}

int main(int argc, char **argv)
{
  FILE *f, *f2;
  int first=1, c2, i, ii, nr1, nr2, a;
  int NN=10, fine, JJ, maxl, nfiles, nat, np, maxnp;
  if (argc <= 1)
    {
      printf("Usage: calcACV <lista_file> [points] [incremento] [num. particelle]\n");
      printf("where points is the number of points of the correlation function\n");
      exit(-1);
    }
 
  c2 = 0;
  f2 = fopen(argv[1], "r");
  maxl = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", dummy); 
      if (strlen(dummy)+1 > maxl)
	maxl = strlen(dummy)+1;
      c2++;
    }	
  nfiles = c2;
  rewind(f2);
  fname = malloc(sizeof(char*)*nfiles);
  for (ii=0; ii < nfiles; ii++)
    {
      fname[ii] = malloc(sizeof(char)*maxl);
      fscanf(f2, "%[^\n]\n", fname[ii]); 
    }

  f = fopen(fname[0], "r");
  nat = 0;
  if (argc >= 3)
    points = atoi(argv[2]);
  else
    points = 1000;
  if (argc == 5)
    NP = atoi(argv[4]);
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  if (nfiles < NN)
    points = nfiles;
  for (a=0; a < 3; a++)
    {
      vel0[a] = malloc(sizeof(double)*NP);
      velt[a] = malloc(sizeof(double)*NP);
    }

  cc = malloc(sizeof(double)*points);
  velACV = malloc(sizeof(double)*points);

  ti = malloc(sizeof(double)*points);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;
  if (argc >= 4)
    NN=atoi(argv[3]);
  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    {
      velACV[ii] = 0.0;
      cc[ii] = 0.0;
    }
  c2 = 0;
  JJ = 0;
  dt=0.005;
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      printf("file=%s\n", fname[nr1]);
      readconf(fname[nr1], vel0);
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      /* N.B. considera NN punti in maniera logaritmica e poi calcola i punti in maniera lineare 
	       * distanziati di NN punti. */
              np = (JJ == 0)?nr2-nr1:NN-1+JJ;	      
	      if (nr2 >= nfiles || np >= points)
		{
		  fine = 1;
		  break;
		}
	      if (JJ > 0 && (nr2 - nr1) % NN != 0)
		continue;
	      readconf(fname[nr2], velt);
	      time = np*dt*10;
	      //printf("nr1=%d nr2=%d JJ=%d time=%.15G\n", nr1, nr2, JJ, time);
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time;
		  //printf("np=%d time=%.15G\n", np, ti[np]);
		}
	      if (nr2 == nr1)
		continue;
              //printf("nr2=%d nr1=%d NP=%d\n", nr1, nr2, NP);
	      for (i=0; i < NP; i++) 
		{
		  veltmp = 0.0;
		  for (a = 0; a < 3; a++)
		   {
	             veltmp += velt[a][i]*vel0[a][i];
		   }
		  velACV[np] += veltmp;
		  cc[np]+=1.0;
		}
	    }
	}
    }

  f = fopen("velACV.dat", "w+");
  for (ii=1; ii < points; ii++)
    {
      //printf("cc[%d]=%f\n", ii, cc[ii]);
      velACV[ii] = velACV[ii]/cc[ii];
      fprintf(f, "%.15G %.15G\n", ti[ii]-ti[0], velACV[ii]);
    }
  fclose(f);
  return 0;
}
