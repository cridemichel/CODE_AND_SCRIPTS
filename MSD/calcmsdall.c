#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 5000
char **fname; 
double L, time, ti, rotMSD, MSD;
double **DR, **DR0;
double *r0[3], *w0[3], *rt[3], *wt[3], *rtold[3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int points, foundDRs=0, foundrot=0;
char fname2[512];

void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double **DR)
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
	  if (!strcmp(parname,"DR"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf %lf %lf ", &DR[i][0], &DR[i][1], &DR[i][2]);
		}
	      foundDRs = 1;
	    }
	  else if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	      foundrot = 1;
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
	      fscanf(f, "%[^\n]\n", line); 
	      if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
		{
		  sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
		}
	    }
	  break; 
	}

    }
  fclose(f);
}
int main(int argc, char **argv)
{
  FILE *f, *f2, *f3, *fA, *fB, *f2A, *f2B;
  double *adjDr[3], Dr, Dw, A1, A2, A3, dr;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int NP, NPA=-1, NN, fine, JJ, nat, maxl, maxnp, np;
  double refTime=0.0;
  double ti0;
  if (argc <= 1)
    {
      printf("Usage: calcmsd <listafile> <particle number>\n");
      //printf("where NN il the lenght of the logarithmic block\n");
      exit(-1);
    }
  f2 = fopen(argv[1], "r");
  c2 = 0;
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
      if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      if (!strcmp(parname,"NN"))
	NN = atoi(parval);
    }
  fclose(f);
  if (NPA == -1)
    NPA = NP;
  np=atoi(argv[2]);
  
  DR = malloc(sizeof(double*)*NP);
  DR0= malloc(sizeof(double*)*NP);
  for (ii = 0; ii < NP; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*3);
      DR0[ii] = malloc(sizeof(double)*3);
    }
  for (ii=0; ii < points; ii++)
    {
      ti = -1.0;
      rotMSD = 0.0;
      MSD = 0.0;
    }
    
  if (NPA != NP)
    printf("[MIXTURE] files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", nfiles, NP, NPA, L, NN, maxl);
  else
    printf("[MONODISPERE] files=%d NP = %d L=%.15G NN=%d maxl=%d\n", nfiles, NP, L, NN, maxl);
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      w0[a] = malloc(sizeof(double)*NP);
      rt[a] = malloc(sizeof(double)*NP);
      rtold[a] = malloc(sizeof(double)*NP);
      wt[a] = malloc(sizeof(double)*NP);
      adjDr[a] = malloc(sizeof(double)*NP); 
    }
  nr1 = 0;
  sprintf(fname2, "MSD-%d.dat", np);
  f = fopen(fname2, "w+");
  if (foundrot)
    {
      sprintf(fname2, "rotMSD-%d.dat", np);
      f2 = fopen(fname2, "w+");
    }
  for (i=0; i < NP; i++)
    for (a=0; a < 3; a++)
      adjDr[a][i] = 0.0;

  readconf(fname[nr1], &time, &refTime, NP, r0, w0, DR0);
  fine = 0;
  ti0 = 0.0;
  i = np;
  for (a=0; a < 3; a++)
    rtold[a][i] = r0[a][i];
 
  for (nr2 = nr1+1; nr2 < nfiles; nr2++)
    {
      /* N.B. considera NN punti in maniera logaritmica e poi calcola i punti in maniera lineare 
       * distanziati di NN punti. */
      
      readconf(fname[nr2], &time, &refTime, NP, rt, wt, DR);
      ti = time + refTime;

      if (nr2 == nr1)
	continue;
      for (a = 0; a < 3; a++)
	{
	  Dw = wt[a][i] - w0[a][i];
	  Dr = rt[a][i] - r0[a][i];
	  dr = rt[a][i] - rtold[a][i];
	  if (foundDRs)
	    {
	      adjDr[a][i] = L*(DR[i][a]-DR0[i][a]); 
	    }
	  else
	    {
	      if (nr2 > nr1 && fabs(dr) > L*0.5)
		if (dr > 0.0)
		  adjDr[a][i] -= L;
		else
		  adjDr[a][i] += L;
	    }
	  //printf("adjDr[%d][%d]:%f\n", a, i, adjDr[a][i]);
	  
      	  MSD = (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
	  if (foundrot)
	    rotMSD = Dw*Dw;
	}
      //printf("cc[%d]=%f ti=%f\n", ii, cc[ii], ti[ii]);
      fprintf(f, "%.15G %.15G\n", ti-ti0, MSD);
      if (foundrot)
	fprintf(f2, "%.15G %.15G\n", ti-ti0, rotMSD);
      //printf("cc[%d]:%f\n", nr2-nr1, cc[nr2-nr1]);
      for (a=0; a < 3; a++)
	rtold[a][i] = rt[a][i];
       
    }
  fclose(f);
  if (foundrot)
    fclose(f2);
  return 0;
}
