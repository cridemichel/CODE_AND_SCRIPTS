#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXPTS 1000
char **fname; 
double time, *ccA, *ccB, *C1A, *C2A, *C4A, *C6A, *C1B, *C2B, *C4B, *C6B, C1m, C2m, C4m, C6m,
       C1mA, C2mA, C4mA, C6mA, C1mB, C2mB, C4mB, C6mB,
       costh2, costh4, costh6, *ti, *u0[3], *ut[3], L, refTime;
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
  FILE *f, *f2, *fA, *fB;
  int first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a;
  int NN, fine, JJ, maxl, nfiles, nat, np, maxnp;
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
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  if (nfiles < NN)
    points = nfiles;
  if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
    assez = 0;
  else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
    assez = 1;
  else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
    assez = 2;
  if (NPA == -1)
    NPA = NP;
  fprintf(stderr, "allocating %d items NN=%d NP=%d num files=%d maxnp=%d assez=%d maxl=%d\n", points, NN, NP, nfiles, maxnp, assez, maxl);
  if (NPA < NP)
    {
      printf("NPA=%d\n", NPA);
    }
  for (a=0; a < 3; a++)
    {
      u0[a] = malloc(sizeof(double)*NP);
      ut[a] = malloc(sizeof(double)*NP);
    }

  ccA = malloc(sizeof(double)*points);
  C1A = malloc(sizeof(double)*points);
  C2A = malloc(sizeof(double)*points);
  C4A = malloc(sizeof(double)*points);
  C6A = malloc(sizeof(double)*points);
  if (NPA < NP)
    {
      C1B = malloc(sizeof(double)*points);
      C2B = malloc(sizeof(double)*points);
      C4B = malloc(sizeof(double)*points);
      C6B = malloc(sizeof(double)*points);
      ccB = malloc(sizeof(double)*points);
    }
  ti = malloc(sizeof(double)*points);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;

  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    {
      C1A[ii] = 0.0;
      C2A[ii] = 0.0;
      C4A[ii] = 0.0;
      C6A[ii] = 0.0;
      ccA[ii] = 0;
      if (NPA < NP)
	{
	  C1B[ii] = 0.0;
	  C2B[ii] = 0.0;
	  C4B[ii] = 0.0;
	  C6B[ii] = 0.0;
	  ccB[ii] = 0;
	}
    }
  c2 = 0;
  JJ = 0;
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      readconf(fname[nr1], &time, &refTime, NP, u0);
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
	      readconf(fname[nr2], &time, &refTime, NP, ut);
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("np=%d time=%.15G\n", np, ti[np]);
		}
	      if (nr2 == nr1)
		continue;
	      for (i=0; i < NP; i++) 
		{
		  costh2 = 0.0;
		  for (a = 0; a < 3; a++)
		   {
	             costh2 += ut[a][i]*u0[a][i];
		   }
		  if (i < NPA)
		    C1A[np] += costh2;
		  else
		    C1B[np] += costh2;

		  costh2 = costh2*costh2;
		  costh4 = costh2*costh2;
		  costh6 = costh4*costh2;
		  if (i < NPA)
		    {
		      C2A[np] += costh2;
		      C4A[np] += costh4;
		      C6A[np] += costh6;
		      ccA[np] += 1.0;
		    }
		  else
		    {
		      C2B[np] += costh2;
		      C4B[np] += costh4;
		      C6B[np] += costh6;
		      ccB[np] += 1.0;
		    }
		}
	    }
	}
    }
  if (NPA < NP)
    {
      fA = fopen("CnA.dat", "w+");
      fB = fopen("CnB.dat", "w+");
    }
  
  f = fopen("Cn.dat", "w+");
  for (ii=1; ii < points; ii++)
    {
      if (NPA < NP)
	{
	  C6m = (231.0*(C6A[ii]+C6B[ii])/(ccA[ii]+ccB[ii]) - 
	      	 315.0*(C4A[ii]+C4B[ii])/(ccA[ii]+ccB[ii]) + 105.0*(C2A[ii]+C2B[ii])/(ccA[ii]+ccB[ii]) -
		 5.0)/16.0;
	  C4m = (35.0*(C4A[ii]+C4B[ii])/(ccA[ii]+ccB[ii]) - 30.0*(C2A[ii]+C2B[ii])/(ccA[ii]+ccB[ii]) 
		 + 3.0) / 8.0;
	  C2m = (3.0*(C2A[ii]+C2B[ii])/(ccA[ii]+ccB[ii]) - 1.0)/2.0;
	  C1m = (C1A[ii] + C1B[ii]) / (ccA[ii]+ccB[ii]);
	  C6mA = (231.0*C6A[ii]/ccA[ii] - 315.0*C4A[ii]/ccA[ii] + 105.0*C2A[ii]/ccA[ii] - 5.0)/16.0;
	  C4mA = (35.0*C4A[ii]/ccA[ii] - 30.0*C2A[ii]/ccA[ii] + 3.0) / 8.0;
	  C2mA = (3.0*C2A[ii]/ccA[ii] - 1.0)/2.0;
	  C1mA = C1A[ii]/ccA[ii];
	  C6mB = (231.0*C6B[ii]/ccB[ii] - 315.0*C4B[ii]/ccB[ii] + 105.0*C2B[ii]/ccB[ii] - 5.0)/16.0;
	  C4mB = (35.0*C4B[ii]/ccB[ii] - 30.0*C2B[ii]/ccB[ii] + 3.0) / 8.0;
	  C2mB = (3.0*C2B[ii]/ccB[ii] - 1.0)/2.0;
	  C1mB = C1B[ii]/ccB[ii];
	}
      else
	{
	  C6m = (231.0*C6A[ii]/ccA[ii] - 
		 315.0*C4A[ii]/ccA[ii] + 105.0*C2A[ii]/ccA[ii] -
		 5.0)/16.0;
	  C4m = (35.0*C4A[ii]/ccA[ii] - 30.0*C2A[ii]/ccA[ii] + 3.0) / 8.0;
	  C2m = (3.0*C2A[ii]/ccA[ii] - 1.0)/2.0;
	  C1m = C1A[ii] / ccA[ii];
	}

      if (ti[ii] > -1.0)
	{
	  if (NPA < NP)
	    {
	      fprintf(fA, "%.15G %.15G %.15G %.15G %.15G\n", ti[ii]-ti[0], C1mA, C2mA, C4mA, C6mA);
	      fprintf(fB, "%.15G %.15G %.15G %.15G %.15G\n", ti[ii]-ti[0], C1mB, C2mB, C4mB, C6mB); 
	    }
	  fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", ti[ii]-ti[0], C1m, C2m, C4m, C6m);
	  //fprintf(f, "%.15G %.15G %.15G %.15G\n", ti[ii]-ti[0], C2m, C4m, C6m);
	}
    }
  fclose(f);
  if (NPA < NP)
    {
      fclose(fA);
      fclose(fB);
    }
  return 0;
}
