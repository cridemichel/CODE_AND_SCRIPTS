#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXPTS 1000
char **fname; 
double time, *cc,*C2, C1, *ti, *u0[3], *ut[3], L, refTime;
int points, assez, NP, NPA;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0, A1, B0, B1, C0, C1;

void readconf(char *fname, double *ti, double *refTime, int NP, double *u[3])
{
  FILE *f;
  double r0, r1, r2, R[3][3];
  int nat=0, i, cpos, a;
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
	  if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
  	  else if (!strcmp(parname, "refTime"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *refTime = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else
	{
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		     &r0, &r1, &r2, &R[0][0], &R[0][1], &R[0][2],
		     &R[1][0], &R[1][1], &R[1][2], &R[2][0], &R[2][1], &R[2][2]); 
	      for (a = 0; a < 3; a++)
		u[a][i] = R[assez][a];
	    }
  	  break; 
	}

    }
  fclose(f);
}

int main(int argc, char **argv)
{
  FILE *f, *f2;
  int first=1, firstp=1, c1, c2, c3, i, ii, nlines, nr1, nr2, a;
  int NN, fine, JJ, maxl, nfiles, nat;
  if (argc <= 1)
    {
      printf("Usage: calcCn <lista_file> [points] \n");
      printf("where points is the number of points of the correlation function\n");
      exit(-1);
    }
 
  c2 = 0;
  f2 = fopen(argv[1], "r");
  maxl = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", dummy); 
      if (strlen(dummy) > maxl)
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
      else if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      else if (!strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (!strcmp(parname, "a"))
       	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &A0, &A1);
	}
      else if (!strcmp(parname, "b"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &B0, &B1);
	}
      else if (!strcmp(parname, "c"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &C0, &C1);
	}
    }
  fclose(f);
  if (argc == 3)
    points = atoi(argv[2]);
  else
    points = NN;
  if (points > nfiles)
    points = nfiles;
  if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
    assez = 0;
  else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
    assez = 1;
  else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
    assez = 2;
  if (NPA == -1)
    NPA = NP;
  fprintf(stderr, "allocating %d items\n", c2);
  for (a=0; a < 3; a++)
    {
      u0[a] = malloc(sizeof(double)*NP);
      ut[a] = malloc(sizeof(double)*NP);
    }

  cc = malloc(sizeof(double)*points);
  C2 = malloc(sizeof(double)*points);
  ti = malloc(sizeof(double)*points);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;

  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    {
      C2[ii] = 0.0;
      cc[ii] = 0;
    }
  c2 = 0;

  for (nr1 = 0; nr1 < nlines; nr1=nr1+NN)
    {	
      readconf(fname[nr1], &time, &refTime, NP, u0);
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      if (nr2 >= nlines || nr2 - nr1 >= points)
		{
		  fine = 1;
		  break;
		}
	      readconf(fname[nr2], &time, &refTime, NP, ut);
	      if (nr2 < points && ti[nr2] == -1.0)
		{
		  ti[nr2] = time + refTime;
		  //printf("nr1=%d time=%.15G\n", nr2, ti[nr2]);
		}
  
	      if (nr2 == nr1)
		continue;

	      for (i=0; i < NP; i++) 
		{
		  for (a = 0; a < 3; a++)
		    C2[nr2-nr1] += ut[i][a]*u0[i][a];
		  cc[nr2-nr1] += 1.0;
		}
	    }
	}
    }

  f = fopen("Cn.dat", "w+");
  for (ii=1; ii < points; ii++)
    {
      C1 = C2[ii];
      C2[ii] = C2[ii]/cc[ii];
      C2[ii] = 1.5*C2[ii] - 0.5;
      if (ti[ii] > -1.0)
	{
	  fprintf(f, "%.15G %.15G %.15G %f\n", ti[ii]-ti[0], C1, C2[ii], cc[ii]);
	}
    }
  fclose(f);
  return 0;
}
