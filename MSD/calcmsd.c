#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXPTS 10000
#define MAXFILES 5000
char fname[MAXFILES][256]; 
double L, time, ti[MAXPTS], rotMSD[MAXPTS], MSD[MAXPTS], cc[MAXPTS];
double *r0[3], *w0[3], *rt[3], *wt[3], *rtold[3];
char parname[128], parval[256000], line[256000];
char dummy[1024];
int points;
void readconf(char *fname, double *ti, int NP, double *r[3], double *w[3])
{
  FILE *f;
  int nat=0, i;
  f = fopen(fname, "r");
  while (!feof(f) && nat < 2) 
    {
      fscanf(f, "%[^:]:", parname);
      if (!strcmp(parname,"sumox"))
	{
	  for (i=0; i < NP; i++)
	    {
	      fscanf(f, " %lf ", w[0][i]); 
	    }
	}
      else if (!strcmp(parname,"sumoy"))
	{
	  for (i=0; i < NP; i++)
	    {
	      fscanf(f, " %lf ", w[1][i]); 
	    }
	}
      else if (!strcmp(parname,"sumoz"))
	{
	  for (i=0; i < NP; i++)
	    {
	      fscanf(f, " %lf ", w[2][i]); 
	    }
	}
      else if (!strcmp(parname, "time"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  *ti = atof(parval);
	}	
      else
	fscanf(f, "%[^\n]\n", parval);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  continue;
	}
      if (nat<2)
	continue;
      for (i = 0; i < NP; i++) 
	fscanf(f, "%lf %lf %lf %[^\n]\n", r[0][i], r[1][i], r[2][i], dummy); 
      break;  
    }
  
  fclose(f);
}
int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  double Dr, Dw, A1, A2, A3, dr;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int NP, NN, fine, JJ, nat;
  if (argc <= 1)
    {
      printf("Usage: calcmsd <listafile>\n");
      //printf("where NN il the lenght of the logarithmic block\n");
      exit(-1);
    }
  for (ii=0; ii < MAXPTS; ii++)
    {
      ti[ii] = -1.0;
      rotMSD[ii] = 0.0;
      MSD[ii] = 0.0;
    }
  f2 = fopen(argv[1], "r");
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname[c2]);
      c2++;
    }	
  nfiles = c2;
  fclose(f2);

  f = fopen(fname[0], "r");
  nat = 0;
  while (!feof(f) && nat < 2) 
    {
      fscanf(f, "%[^\n]\n)", line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  if (nat==2)
	    {
	      for (i=0; i < 2*NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      fscanf(f, "%lf\n", &L);
	      break;
	    }
	  continue;
	}
      sscanf(line, "%[^:]:%[^\n]\n", parname, parval); 
      if (!strcmp(parname,"parnum"))
	NP = atoi(parval);
      if (!strcmp(parname,"NN"))
	NN = atoi(parval);
    }
  fclose(f);
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      w0[a] = malloc(sizeof(double)*NP);
      rt[a] = malloc(sizeof(double)*NP);
      rtold[a] = malloc(sizeof(double)*NP);
      wt[a] = malloc(sizeof(double)*NP);
      
    }
  printf("NP = %d\n", NP);
  if (c2 < MAXPTS && ti[c2] == -1.0)
    {
      ti[c2] = time;
      //printf("c2=%d time=%.15G\n", c2, ti[c2]);
    }
  //printf("%d fname: %s %.15G %.15G %.15G\n", c2, fname, P0[0], P0[1], P0[2]);
  c2++;
  fclose(f2);
  //c2 = 0;

  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      readconf(fname[nr1], &time, NP, r0, w0);
      if (nr1 < MAXPTS && ti[nr1] == -1.0)
	{
	  ti[nr1] = time;
	  //printf("c2=%d time=%.15G\n", c2, ti[c2]);
	}
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      if (nr2 == nlines || nr2 - nr1 >= points)
		{
		  fine = 1;
		  break;
		}
	      
	      if (nr2 > nr1)
		{
		  for (i=0; i < NP; i++)
		    for (a=0; a < 3; a++)
		      rtold[a][i] = rt[a][i];
		}
	      readconf(fname[nr2], &time, NP, rt, wt);
	      for (i = 0; i < NP; i++)
		for (a = 0; a < 3; a++)
		  {
		    Dw = wt[a][i] - w0[a][i];
		    Dr = rt[a][i] - r0[a][i];
		    dr = rt[a][i] - rtold[a][i];
		    if (nr2 > nr1 && fabs(dr) > L*0.5)
		      Dr -= L;
		    MSD[nr2-nr1] += Dr*Dr;
		    rotMSD[nr2-nr1] += Dw*Dw;
		  }
	      cc[nr2-nr1] += 1.0;
	    }
	}
    }
  f = fopen("MSD.dat", "w+");
  f2 = fopen("rotMSD.dat", "w+");

  for (ii=1; ii < points; ii++)
    {
      if (cc[ii] > 0 && ti[ii] > -1.0)
	{
	  fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], MSD[ii]/cc[ii], cc[ii]);
	  fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSD[ii]/cc[ii], cc[ii]);
	}
    }
  fclose(f);
  fclose(f2);
	
  return 0;
}
