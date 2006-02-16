#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 1000
char **fname; 
double time, *ti, *r0[3], *r1[3], L, refTime;
int points, assez, NP, NPA;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0, A1, B0, B1, C0, C1;

void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3])
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
#define KMODMAX 99
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double *sqRe[KMODMAX], *sqIm[KMODMAX];
double *cc[KMODMAX];
char fname2[512];
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
int main(int argc, char **argv)
{
  FILE *f, *f2;
  int first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a;
  int iq, NN, fine, JJ, maxl, nfiles, nat, np, maxnp;
  int qmin = 5, qmax = 30, qmod; 
  double invL, rxdummy, sumIm, sumRe, scalFact;

  twopi = acos(0)*4.0;	  
  if (argc <= 1)
    {
      printf("Usage: calcfqself <lista_file> [points] [qmin] [qmax]\n");
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
  invL = 1.0/L;
  if (argc == 3)
    points = atoi(argv[2]);
  else
    points = NN;
  if (argc == 4)
    qmin = atoi(argv[3]);
  if (argc == 5)
    qmax = atoi(argv[4]);
  
  scalFact = twopi * invL;
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
    assez = 0;
  else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
    assez = 1;
  else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
    assez = 2;
  if (NPA == -1)
    NPA = NP;
  //fprintf(stderr, "allocating %d items NN=%d NP=%d num files=%d maxnp=%d\n", points, NN, NP, nfiles, maxnp);
  for (qmod = qmin; qmod <= qmax; qmod++)
    {
      sqRe[qmod] = malloc(sizeof(double)*points);
      sqIm[qmod] = malloc(sizeof(double)*points);
      cc[qmod] = malloc(sizeof(double)*points);
    }
  ti = malloc(sizeof(double)*points);
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      r1[a] = malloc(sizeof(double)*NP);
    }

  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;

  first = 0;
  fclose(f2);
  for (qmod = qmin; qmod <= qmax; qmod++)
    {
      for (ii=0; ii < points; ii++)
	{
	  sqRe[qmod][ii] = 0.0;
	  sqIm[qmod][ii] = 0.0;
	  cc[qmod][ii] = 0.0;
	}
      for (iq=0; iq < ntripl[qmod]; iq++)
	{
	  qx[qmod][iq]=invL*twopi*mesh[qmod][iq][0];
	  qy[qmod][iq]=invL*twopi*mesh[qmod][iq][1];
	  qz[qmod][iq]=invL*twopi*mesh[qmod][iq][2];
	}
    }
  c2 = 0;
  JJ = 0;
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      readconf(fname[nr1], &time, &refTime, NP, r0);
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
	      readconf(fname[nr2], &time, &refTime, NP, r1);
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("np=%d time=%.15G\n", np, ti[np]);
		}
  
	      //if (nr2 == nr1)
		//continue;

	      for (qmod = qmin; qmod <= qmax; qmod++)
		{
		  for(iq=0; iq < ntripl[qmod]; iq++)
		    {
		      sumRe=0.0;
		      sumIm=0.0;
		      for (i=0; i < NP; i++)
			{
			  rxdummy = scalFact*((r0[0][i]-r1[0][i])*mesh[qmod][iq][0]
			    +(r0[1][i]-r1[1][i])*mesh[qmod][iq][1]
			    +(r0[2][i]-r1[2][i])*mesh[qmod][iq][2]);
			  //printf("dummy:%.15G\n", rxdummy);
			  sumRe += cos(rxdummy);
			  sumIm += sin(rxdummy);	
			}
		      sqRe[qmod][np] += sumRe;
		      sqIm[qmod][np] += sumIm;
		    }
		  sqRe[qmod][np] /= ((double)NP);
		  sqIm[qmod][np] /= ((double)NP);
		  cc[qmod][np] += 1.0;
		  //printf("cc[%d][%d]=%.15G sqre=%.15G sqim=%.15G\n", qmod, np, cc[qmod][np],
		//	 sqRe[qmod][np], sqIm[qmod][np]);
		}
	    }
	}
    }
  for (qmod = qmin; qmod <= qmax; qmod++)
    {
      for (ii=0; ii < points; ii++)
	{
	  sqRe[qmod][ii] = sqRe[qmod][ii]/cc[qmod][ii];
	  sqIm[qmod][ii] = sqIm[qmod][ii]/cc[qmod][ii];
	 // printf("qmod=%d ii=%d sqre=%.15G sqIm=%.15G cc=%.15G\n",  qmod, ii, sqRe[qmod][ii], sqIm[qmod][ii],
	//	 cc[qmod][ii]);
	}
    }
  for (qmod = qmin; qmod < qmax; qmod++)
    {
      sprintf(fname2, "Fqs-%d",qmod);
      f = fopen (fname2, "w+");
      for (ii = 0; ii < points; ii++)
	{
	  if ((sqRe[qmod][ii]!=0.0 || sqIm[qmod][ii]!=0.0) && (ti[ii]> -1.0))
	    fprintf(f, "%15G %.15G %.15G\n", ti[ii]-ti[0], sqRe[qmod][ii]/sqRe[qmod][0], 
		    sqIm[qmod][ii]);
	}
    }
  fclose(f);
  return 0;
}
