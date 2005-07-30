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
  int nat=0, i, cpos;
  f = fopen(fname, "r");
  while (!feof(f) && nat < 2) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	}
      if (nat < 2)
	{
	  fseek(f, cpos, SEEK_SET);
	  fscanf(f, "%[^:]:", parname);
	  //printf("[%s] parname=%s\n", fname, parname);
	  if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoy"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[1][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoz"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[2][i]); 
		}
	    }
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else
	{
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
	    }
	  break; 
	}

    }
  fclose(f);
}
int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  double *adjDr[3], Dr, Dw, A1, A2, A3, dr;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int NP, NN, fine, JJ, nat;
  if (argc <= 1)
    {
      printf("Usage: calcmsd <listafile> [number of points]\n");
      //printf("where NN il the lenght of the logarithmic block\n");
      exit(-1);
    }
  for (ii=0; ii < MAXPTS; ii++)
    {
      ti[ii] = -1.0;
      rotMSD[ii] = 0.0;
      MSD[ii] = 0.0;
      cc[ii]=0.0;
    }
  f2 = fopen(argv[1], "r");
  c2 = 0;
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
  if (argc == 3)
    points=atoi(argv[2]);
  else
    points=NN;
 
  printf("points=%d files=%d NP = %d L=%.15G NN=%d\n", points, nfiles, NP, L, NN);
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      w0[a] = malloc(sizeof(double)*NP);
      rt[a] = malloc(sizeof(double)*NP);
      rtold[a] = malloc(sizeof(double)*NP);
      wt[a] = malloc(sizeof(double)*NP);
      adjDr[a] = malloc(sizeof(double)*NP); 
    }
 for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      for (i=0; i < NP; i++)
	for (a=0; a < 3; a++)
	  adjDr[a][i] = 0.0;

      readconf(fname[nr1], &time, NP, r0, w0);
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      if (nr2 >= nfiles || nr2 - nr1 >= points)
		{
		  fine = 1;
		  break;
		}
		
	      if (nr2==nr1)
		{
		  for (i=0; i < NP; i++)
		    for (a=0; a < 3; a++)
		      rtold[a][i] = r0[a][i];
		}
	      else
		{
		  for (i=0; i < NP; i++)
		    for (a=0; a < 3; a++)
		      rtold[a][i] = rt[a][i];
		}
	      readconf(fname[nr2], &time, NP, rt, wt);
	      if (nr2 < MAXPTS && ti[nr2] == -1.0)
		{
		  ti[nr2] = time;
		  printf("nr1=%d time=%.15G\n", nr2, ti[nr2]);
		}
  
	      if (nr2 == nr1)
		continue;
	      for (i = 0; i < NP; i++)
		for (a = 0; a < 3; a++)
		  {
		    Dw = wt[a][i] - w0[a][i];
		    Dr = rt[a][i] - r0[a][i];
		    dr = rt[a][i] - rtold[a][i];
		    if (nr2 > nr1 && fabs(dr) > L*0.5)
		      if (dr > 0.0)
			adjDr[a][i] -= L;
		      else
			adjDr[a][i] += L;
		    //printf("adjDr[%d][%d]:%f\n", a, i, adjDr[a][i]);
		    MSD[nr2-nr1] += (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
		    rotMSD[nr2-nr1] += Dw*Dw;
		  }
	      cc[nr2-nr1] += 1.0;
	      //printf("cc[%d]:%f\n", nr2-nr1, cc[nr2-nr1]);
	    }
	}
    }
  f = fopen("MSDcnf.dat", "w+");
  f2 = fopen("rotMSDcnf.dat", "w+");

  for (ii=1; ii < points; ii++)
    {
      printf("cc[%d]=%f ti=%f\n", ii, cc[ii], ti[ii]);
      if (cc[ii] > 0 && ti[ii] > -1.0)
	{
	  fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], MSD[ii]/cc[ii]/((double)NP), cc[ii]);
	  fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSD[ii]/cc[ii]/((double)NP), cc[ii]);
	}
    }
  fclose(f);
  fclose(f2);
	
  return 0;
}
