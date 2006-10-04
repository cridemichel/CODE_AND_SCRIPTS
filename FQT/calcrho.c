#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<sys/types.h>
#include<dirent.h>
#include<errno.h>
#define MAXPTS 1000
char **fname; 
double time, *r0[3], L, refTime;
int points, assez, NP, NPA;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0, A1, B0, B1, C0, C1;

void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3])
{
  FILE *f;
  //double R[3][3];
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
#define KMODMAX 598
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double *cc[KMODMAX];
char fname2[512];
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
int main(int argc, char **argv)
{
  FILE *f, *f2, *fA, *fB;
  DIR* dir;
  char mode[10];
  int c2, i, ii, nr1, a;
  int iq, NN, maxl, nfiles, nat;
  int qmin = 0, qmax = KMODMAX, qmod; 
  double invL, rxdummy, sumImA, sumReA, sumImB, sumReB, scalFact;

  twopi = acos(0)*4.0;	  
  if (argc <= 1)
    {
      printf("Usage: calcrho <lista_file> [qmin] [qmax]\n");
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
      //printf("file=%s\n", fname[ii]);
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
		{
		  fscanf(f, "%[^\n]\n", line);
		}
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
  invL = 1.0/L;
  printf("qui NP=%d argv=%s %s %s %s\n", NP, argv[0], argv[1], argv[2], argv[3]);
  if (argc >= 2)
    qmin = atoi(argv[2]);
  if (argc == 4)
    qmax = atoi(argv[3]);
  
  scalFact = twopi * invL;
  //printf("maxnp=%d points=%d\n",maxnp, points);
  printf("qmin = %d qmax=%d\n", qmin, qmax);
  if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
    assez = 0;
  else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
    assez = 1;
  else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
    assez = 2;
  if (NPA == -1)
    NPA = NP;
  //fprintf(stderr, "allocating %d items NN=%d NP=%d num files=%d maxnp=%d\n", points, NN, NP, nfiles, maxnp);
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
    }
  for (qmod = qmin; qmod <= qmax; qmod++)
    {
      for (iq=0; iq < ntripl[qmod]; iq++)
	{
	  qx[qmod][iq]=invL*twopi*mesh[qmod][iq][0];
	  qy[qmod][iq]=invL*twopi*mesh[qmod][iq][1];
	  qz[qmod][iq]=invL*twopi*mesh[qmod][iq][2];
	}
    }
  if (NPA < NP)
    {
      dir=opendir("./RHOTMPA");
      if (errno==ENOENT)
	system("mkdir RHOTMPA/");
      else
	closedir(dir);
      dir=opendir("./RHOTMPB");
      if (errno==ENOENT)
	system("mkdir RHOTMPB/");
      else
	closedir(dir);
      f=fopen("./RHOTMPA/NN.dat","w");
      fprintf(f, "%d\n", NN);
      fclose(f);
      f=fopen("./RHOTMPB/NN.dat","w");
      fprintf(f, "%d\n", NN);
      fclose(f);
    }
  else
    {
      dir=opendir("./RHOTMP");
      if (errno==ENOENT)
	system("mkdir RHOTMP/");
      else
	closedir(dir);
      f=fopen("./RHOTMP/NN.dat","w");
      fprintf(f, "%d\n", NN);
      fclose(f);

    }

  for (nr1 = 0; nr1 < nfiles; nr1++)
    {	
      readconf(fname[nr1], &time, &refTime, NP, r0);
      for (qmod = qmin; qmod <= qmax; qmod++)
	{
	  if (qmod == qmin)
	    strcpy(mode,"w");
	  else
	    strcpy(mode,"a");
	  if (NPA < NP)
	    {
	      sprintf(fname2,"RHOTMPA/ro.00.k=%03d", qmod);
	      fA=fopen(fname2,mode);
	      sprintf(fname2,"RHOTMPB/ro.00.k=%03d", qmod);
	      fB=fopen(fname2,mode);
      	    }
	  else 
	    {
	      sprintf(fname2,"RHOTMP/ro.00.k=%03d", qmod);
	      f=fopen(fname2,mode);
	    }
	  if (NPA < NP)
	    {
	      //printf("refTime=%.15G nf=%d\n", refTime, nr1);
	      fprintf(fA, "%.15G %d\n", time+refTime, ntripl[qmod]);
	      fprintf(fB, "%.15G %d\n", time+refTime, ntripl[qmod]);
	    }
	  else
	    fprintf(f, "%.15G %d\n", time+refTime, ntripl[qmod]);
	  for(iq=0; iq < ntripl[qmod]; iq++)
	    {
	      sumReA=sumReB=0.0;
	      sumImA=sumImB=0.0;
	      for (i=0; i < NP; i++)
		{
		  rxdummy = scalFact*(r0[0][i]*mesh[qmod][iq][0]
	  			      +r0[1][i]*mesh[qmod][iq][1]
				      +r0[2][i]*mesh[qmod][iq][2]);
		  //printf("dummy:%.15G\n", rxdummy);
		  if (i < NPA)
		    {
		      sumReA += cos(rxdummy);
		      sumImA += sin(rxdummy);
		    }
		  else
		    {
		      sumReB += cos(rxdummy);
		      sumImB += sin(rxdummy);
		    }  
		}
	      if (NPA < NP)
		{
		  fprintf(fA, "(%.15G,%.15G) ", sumReA, sumImA);
		  fprintf(fB, "(%.15G,%.15G) ", sumReB, sumImB);
		}
	      else
		fprintf(f, "(%.15G,%.15G) ", sumReA, sumImA);
	    }

	  if (NPA < NP)
	    {
	      fprintf(fA, "\n");
	      fprintf(fB, "\n");
	      fclose(fA);
	      fclose(fB);
	    }
	  else
	    {
	      fprintf(f,"\n");
	      fclose(f);
	    }
	}
    }
  return 0;
}
