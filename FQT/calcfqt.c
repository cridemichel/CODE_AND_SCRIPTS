/*
 * =====================================================================================
 * 
 *        Filename:  evalcorr.c
 * 
 *     Description:  time averageg correl. func.
 * 
 *         Version:  1.0
 *         Created:  19/07/2005 16:52:36 CEST
 *        Revision:  none
 *        Compiler:  gcc
 * 
 *          Author:   (), 
 *         Company:  
 * 
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXPTS 1000
#define MAXQ 50
#define NUMQ 100
char fname[1024]; 
double time, Cav0, Cav, rhoR0[MAXQ], rhoI0[MAXQ], *cc[MAXQ],*C[MAXQ], 
       rhoR1[MAXQ], rhoI1[MAXQ], *ti, *rhoRt[MAXQ], *rhoIt[MAXQ],
       *rhoRtp[MAXQ], *rhoItp[MAXQ], *tiall;
int NQarr[NUMQ];
int points;
int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  double A1, A2, A3 ;
  int first=1, firstp=1, NQ, nq, c1, c2, c3, i, ii, nlines, nr1, nr2, ll, mm, llp, mmp;
  int NN, fine, JJ, maxnp, np;
  if (argc <= 1)
    {
      printf("Usage: calcfqt <l> <m> <l'> <m'> <NN> [points] \n");
      printf("where NN is the number of configurations in a logarithmic block\n");
      exit(-1);
    }
  NN =  atoi(argv[5]);
  if (argc == 7)
    points = atoi(argv[6]);
  else
    points = NN;
  ll = atoi(argv[1]);
  mm = atoi(argv[2]);
  llp = atoi(argv[3]);
  mmp = atoi(argv[4]);
  c2 = 0;
  sprintf(fname, "RHOTMP/ro.%02d.k=%03d", 0, 2);
  f2 = fopen(fname, "r");
  while (!feof(f2))
    {
      fscanf(f2, "%lf %d ", &time, &NQ);
      //printf("time=%f nq=%d\n", time, NQ);
      for (i = 0; i < NQ; i++)
	fscanf(f2, "(%lf,%lf) ", &A1, &A2);
      c2++;
    }	
  maxnp = NN + (c2-NN)/NN;
  nlines = c2;
  if (points > maxnp)
    points = maxnp;
  fprintf(stderr, "allocating %d items\n", nlines);
  for (i=0; i < MAXQ; i++)
    {
      rhoRt[i] = malloc(sizeof(double)*nlines);
      rhoIt[i] = malloc(sizeof(double)*nlines);
      rhoRtp[i] = malloc(sizeof(double)*nlines);
      rhoItp[i] = malloc(sizeof(double)*nlines);
      cc[i] = malloc(sizeof(double)*points);
      C[i]  = malloc(sizeof(double)*points);
    }
  ti = malloc(sizeof(double)*points);
  tiall = malloc(sizeof(double)*nlines);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;
  for (ii=0; ii < nlines; ii++)
    tiall[ii] = -1.0;
  first = 0;
  fclose(f2);

  for (nq = 2; nq < NUMQ; nq++)
    {
      for (i=0; i < MAXQ; i++)
	for (ii=0; ii < points; ii++)
	  {
	    C[i][ii] = 0.0;
	    cc[i][ii] = 0;
	  }
      fprintf(stderr, "nf=%d\n", nq);
      sprintf(fname, "RHOTMP/ro.%02d.k=%03d", ll*10+mm,nq);
      printf("Opening fname=%s\n", fname);
      c2 = 0;
      f2 = fopen(fname, "r");
      while (!feof(f2))
	{
	  //printf("reading c2=%d\n", c2);
	  fscanf(f2, "%lf %d ", &time, &NQ);
	  
	  if (c2==0)
	    NQarr[nq] = NQ;
  	  for (i = 0; i < NQ; i++)
	    {
	      fscanf(f2, "(%lf,%lf) ", &(rhoRt[i][c2]), &(rhoIt[i][c2]));
	    }
	  if (tiall[c2] == -1.0)
	    {
	      tiall[c2] = time;
	      //printf("c2=%d time=%.15G\n", c2, ti[c2]);
	    }
	  //printf("%d fname: %s %.15G %.15G %.15G\n", c2, fname, P0[0], P0[1], P0[2]);
	  c2++;
	}
      fclose(f2);
      //c2 = 0;
      sprintf(fname, "RHOTMP/ro.%02d.k=%03d", llp*10+mmp,nq);
    
      c2 = 0;
      f2 = fopen(fname, "r");
      while (!feof(f2))
	{
	  //printf("reading c2=%d\n", c2);
	  fscanf(f2, "%lf %d ", &time, &NQ);
  	  for (i = 0; i < NQ; i++)
	    {
	      fscanf(f2, "(%lf,%lf) ", &(rhoRtp[i][c2]), &(rhoItp[i][c2]));
	     //printf("NQ=%d nq=%d QUI i=%d c2=%d rhoRtp: %.15G  rhoItp: %.15G \n",NQ, nq, i, c2,
	     //	     rhoRtp[i][c2], rhoItp[i][c2]);
	    }
  	  c2++;
	}
      fclose(f2);
      for (nr1 = 0; nr1 < nlines; nr1=nr1+NN)
	{	
	  for (i=0; i < NQ; i++)
	    {
	      rhoR0[i] = rhoRt[i][nr1];
	      rhoI0[i] = rhoIt[i][nr1]; 
	    }
	  fine = 0;
	  for (JJ = 0; fine == 0; JJ++)
	    {
	      for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
		{
		  np = (JJ == 0)?nr2-nr1:NN-1+JJ;	      
		  if (nr2 >= nlines || np >= points)
		    {
		      fine = 1;
		      break;
		    }
		  if (JJ > 0 && (nr2 - nr1) % NN != 0)
		    continue;
		  for (i=0; i < NQ; i++) 
		    {
		      //printf("i=%d NQ=%d rhoRtp=%f rhoItp=%f\n", i, NQ, rhoRtp[i][nr2], rhoItp[i][nr2]);
		      rhoR1[i] = rhoRtp[i][nr2];
		      rhoI1[i] = rhoItp[i][nr2];
		      C[i][np] += rhoR1[i]*rhoR0[i]+rhoI1[i]*rhoI0[i];
		      //printf("C[%d][%d]: %.15G cc:%f\n", i, nr2-nr1, C[i][nr2-nr1], cc[i][nr2-nr1]);
		      cc[i][np] += 1.0;
		      //printf("qui c3-c2=%d\n", c3-c2);
		    }
		  if (np < points && ti[np] == -1.0)
		    {
		      ti[np] = tiall[nr2];
		      //printf("np=%d time=%.15G\n", np, ti[np]);
		    }

		}
	    }
	}
      sprintf(fname, "sqt.%02d%02d.k=%03d",ll*10+mm,llp*10+mmp,nq);
      f = fopen(fname, "w+");
      sprintf(fname, "N-sqt.%02d%02d.k=%03d",ll*10+mm,llp*10+mmp,nq);
      f2 = fopen(fname, "w+");
      Cav0 = 0.0;
      for (i = 0; i < NQarr[nq]; i++)
	{
	  if (cc[i][0] > 0)
	    Cav0 += C[i][0]/cc[i][0];
	}
      Cav0 /= ((double)NQarr[nq]);
      
      for (ii=1; ii < points; ii++)
	{
	  Cav = 0.0;
	  for (i = 0; i < NQarr[nq]; i++)
	    {
	      if (cc[i][ii] > 0)
		Cav += C[i][ii]/cc[i][ii];
	    }
	  Cav /= ((double)NQarr[nq]);
	  if (ti[ii] > -1.0)
	    {
	      fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], Cav, cc[0][ii]*NQarr[nq]);
	      fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], Cav/Cav0, cc[0][ii]*NQarr[nq]);
	    }
	}
      fclose(f);
      fclose(f2);
    }
  for (i=0; i < MAXQ; i++)
    {
      free(rhoIt[i]);
      free(rhoRt[i]);
      free(rhoItp[i]);
      free(rhoRtp[i]);
    }
	
  return 0;
}
