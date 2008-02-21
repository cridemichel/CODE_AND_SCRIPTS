#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXPTS 1000
char **fname; 
double time, *cc, veltmp, omtmp, *ti, *vel0[3], *velt[3], *omega0[3], *omegat[3], 
       L, refTime=0.0, *omACV, *velACV;
int points, assez, NP, NPA=-1;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0=-1.0, A1, B0, B1, C0, C1;
double dt=-1.0;
void readconf(char *fname, double *ti, double *refTime, int NP, double *vel[3], 
	      double *omega[3])
{
  FILE *f;
  double r0, r1, r2, R[3][3], v[3], w[3];
  int foundtime=0, nat=0, i, cpos, a;
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
	  if (nat==1 && !foundtime && !strcmp(parname,"curStep"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atoi(parval)*dt;
	    }
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      foundtime=1;
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
	      if (A0 > 0.0)
		fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		       &r0, &r1, &r2, &R[0][0], &R[0][1], &R[0][2],
		       &R[1][0], &R[1][1], &R[1][2], &R[2][0], &R[2][1], &R[2][2]); 
	      else
		fscanf(f, "%lf %lf %lf\n", &r0, &r1, &r2); 

	      //for (a = 0; a < 3; a++)
		//u[a][i] = R[assez][a];
	    }
	  for (i = 0; i < NP; i++) 
	    {
	      if (A0 > 0.0)
		fscanf(f, "%lf %lf %lf %lf %lf %lf\n", 
		       &v[0], &v[1], &v[2], &w[0], &w[1], &w[2]); 
	      else
		fscanf(f, "%lf %lf %lf\n", 
		       &v[0], &v[1], &v[2]); 

	      for (a = 0; a < 3; a++)
		{
		  vel[a][i] = v[a];
		  if (A0 > 0.0)
		    omega[a][i] = w[a]; 
		}
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
      else if (!strcmp(parname,"NN") && nat==0)
	NN = atoi(parval);
      else if (!strcmp(parname,"steplength"))
	dt = atof(parval);
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
  if (argc == 5)
    NP = atoi(argv[4]);
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
  if (dt<0.0)
    fprintf(stderr, "[ED SIM] allocating %d items NN=%d NP=%d num files=%d maxnp=%d assez=%d maxl=%d\n", points, NN, NP, nfiles, maxnp, assez, maxl);
  else
    fprintf(stderr, "[MD SIM] allocating %d items NN=%d NP=%d num files=%d maxnp=%d assez=%d maxl=%d dt=%.15G\n", points, NN, NP, nfiles, maxnp, assez, maxl, dt);
  if (A0 > 0.0)
    printf("[ELLIPSOIDS]\n"); 
  if (NPA < NP)
    {
      printf("NPA=%d\n", NPA);
    }
  for (a=0; a < 3; a++)
    {
      omega0[a] = malloc(sizeof(double)*NP);
      omegat[a] = malloc(sizeof(double)*NP);
      vel0[a] = malloc(sizeof(double)*NP);
      velt[a] = malloc(sizeof(double)*NP);
    }

  cc = malloc(sizeof(double)*points);
  omACV = malloc(sizeof(double)*points);
  velACV = malloc(sizeof(double)*points);
  ti = malloc(sizeof(double)*points);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;
  if (argc==4)
    NN=atoi(argv[3]);
  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    {
      velACV[ii] = 0.0;
      omACV[ii] = 0.0;
    }
  c2 = 0;
  JJ = 0;
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      readconf(fname[nr1], &time, &refTime, NP, vel0, omega0);
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
	      readconf(fname[nr2], &time, &refTime, NP, velt, omegat);
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("np=%d time=%.15G\n", np, ti[np]);
		}
	      if (nr2 == nr1)
		continue;
	      for (i=0; i < NP; i++) 
		{
		  veltmp = 0.0;
		  omtmp = 0.0;
		  for (a = 0; a < 3; a++)
		   {
	             veltmp += velt[a][i]*vel0[a][i];
		     if (A0 > 0.0)
		       omtmp += omegat[a][i]*omega0[a][i];
		   }
		  velACV[np] += veltmp;
		  if (A0 > 0.0)
		    omACV[np] += omtmp;
		  cc[np]+=1.0;
		}
	    }
	}
    }
  f = fopen("velACV.dat", "w+");
  if (A0 > 0.0)
    f2= fopen("omACV.dat", "w+");
  for (ii=1; ii < points; ii++)
    {
      //printf("cc[%d]=%f\n", ii, cc[ii]);
      velACV[ii] = velACV[ii]/cc[ii];
      if (A0 > 0.0)
	omACV[ii] = omACV[ii]/cc[ii];
      if (ti[ii] > -1.0)
	{
	  fprintf(f, "%.15G %.15G\n", ti[ii]-ti[0], velACV[ii]);
	  if (A0 > 0.0)
	    fprintf(f2, "%.15G %.15G\n", ti[ii]-ti[0], omACV[ii]);
	}
    }
  fclose(f);
  if (A0 > 0.0)
    fclose(f2);
  return 0;
}
