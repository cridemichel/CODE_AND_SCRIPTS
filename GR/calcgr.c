#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define Sqr(x) ((x)*(x))
char line[100000], parname[124], parval[256000];
char dummy[2048];
int N;
double *x[3], L, ti, *w[3];
double **DR;
char fname[1024], inputfile[1024];
int points;
int foundDRs=0, foundrot=0;
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double **DR)
{
  FILE *f;
  int nat=0, i, cpos;
  double dt=-1;
  int curstp=-1;

  *ti = -1.0;
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
	  else if (!strcmp(parname, "curStep"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      curstp = atoi(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else if (!strcmp(parname, "steplength"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      dt = atof(parval);
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
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
  if (*ti == -1)
    *ti = ((double)curstp)*dt;
  fclose(f);
}

void print_usage(void)
{
  printf("calcgr <confs_file> [points]\n");
  exit(0);
}

void parse_param(int argc, char** argv)
{
  int cc=1;
  
  if (argc==1)
    {
      print_usage();
      exit(1);
    }
  strcpy(inputfile,argv[1]);
  if (argc == 3)
    points=atoi(argv[2]);
  else
    points=100;
}
double pi;
int main(int argc, char** argv)
{
  FILE *f, *f2;
  int nf, i, a, b, nat, NPA, NP, NN, j, ii, *g0, bin;
  double r, delr, tref=0.0, Dx[3], g0m, distSq, rlower, rupper, cost, nIdeal;
  double time, refTime;
#if 0
  double g2m, g4m, g6m;
  double *g2, *g4, *g6;
#endif
  parse_param(argc, argv);
  f2 = fopen(inputfile,"r");
  fscanf(f2, "%[^\n]\n", fname);
  pi = acos(0)*2.0;
  nat = 0;
  f = fopen(fname,"r");
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
  for (a = 0; a < 3; a++)
    {
      x[a] = malloc(sizeof(double)*NP);
      w[a] = malloc(sizeof(double)*NP);
    }
  g0 = malloc(sizeof(double)*points);
  DR = malloc(sizeof(double*)*NP);
  for (ii = 0; ii < NP; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*3);
    }
#if 0
  g2 = malloc(sizeof(double)*points);
  g4 = malloc(sizeof(double)*points);
  g6 = malloc(sizeof(double)*points);
#endif
  delr = (L*sqrt(3.0)) / ((double)points);
  rewind(f2);
  nf = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      nf++;
      readconf(fname, &time, &refTime, NP, x, w, DR);

      for (i = 0; i < NP; i++)
	for (j = 0; j < i; j++)
	  {
	    for (a = 0; a < 3; a++)
	      {
		Dx[a] = x[a][i] - x[a][j];
	      }
	    for (a = 0; a < 3; a++)
	      Dx[a] = Dx[a] - L * rint(Dx[a]/L);
	   distSq = 0.0;
	   for (a = 0; a < 3; a++)
	     distSq += Sqr(Dx[a]);
	   bin = ((int) (sqrt(distSq) / delr)); 
	   if (bin < points || bin >= 0)
	     {
	       g0[bin] += 2.0;
	     }
	   else 
	     {
	       printf("bin=%d\n", bin);
	       //exit(1);
	     }
	  }
    }
  fclose(f2); 
  f = fopen("gr.dat", "w+");
  r = 0.0;
  cost = 4.0 * pi * NP / 3.0 / (L*L*L);
  for (ii = 0; ii < points; ii++)
    {
      r += delr*(ii+0.5);
      rlower = ( (double) ii) * delr;
      rupper = rlower + delr;
      nIdeal = cost * (Sqr(rupper)*rupper - Sqr(rlower)*rlower);
      g0m = g0[ii]/((double)nf)/nIdeal;
#if 0
      g2m = (3.0*g2[ii]/cc[ii] - 1.0)/2.0;
      g4m = (35.0*g4[ii]/cc[ii] - 30.0*g2[ii]/cc[ii] + 3.0) / 8.0;
      g6m = (231.0*g6[ii]/cc[ii] - 315.0*g4[ii]/cc[ii] + 105.0*g2[ii]/cc[ii] - 5.0)/16.0;
      fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", r, g0m, g2m, g4m, g6m);
#endif
      fprintf(f, "%.15G %.15G\n", r, g0m);
    }
  fclose(f);
}
