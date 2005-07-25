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
#define MAXQ 100
#define NUMQ 100
char fname[1024]; 
double time, Cav, rhoR0[MAXQ], rhoI0[MAXQ], cc[3][MAXPTS],C[3][MAXPTS], 
       rhoR1[MAXQ], rhoI1[MAXQ], ti[MAXPTS], *rhoRt[MAXQ], *rhoIt[MAXQ];
int NQarr[NUMQ];
int points;
int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  double A1, A2, A3 ;
  int NQ, nq, c1, c2, c3, i, ii, nlines, nr1, nr2, ll, mm;
  if (argc <= 1)
    {
      printf("you must supply a filename!\n");
      exit(-1);
    }
  if (argc == 4)
    points = atoi(argv[2]);
  else
    points = MAXPTS;
  for (i=0; i < 3; i++)
    for (ii=0; ii < MAXPTS; ii++)
      {
	C[i][ii] = 0.0;
	cc[i][ii] = 0;
      }
  ll = atoi(argv[1]);
  mm = atoi(argv[2]);
  for (ii=0; ii < MAXPTS; ii++)
    ti[ii] = -1.0;
  for (nq = 2; nq < NUMQ; nq++)
    {
      fprintf(stderr, "nf=%d\n", nq);
      c2 = 0;
      sprintf(fname, "RHOTMP/ro.%2d.k=%3d", ll*10+mm,nq);
      f2 = fopen(fname, "r");
      while (!feof(f2))
	{
	  fscanf(f2, "%lf %d ", &time, &NQ);
	  for (i = 0; i < NQ; i++)
	    fscanf(f2, "(%lf,%lf) ", &A1, &A2);
	  c2++;
	}	 
      fprintf(stderr, "allocating %d items\n", c2);
      NQarr[nq] = NQ;
      for (i=0; i < NQ; i++)
	{
	  rhoRt[i] = malloc(sizeof(double)*c2);
	  rhoIt[i] = malloc(sizeof(double)*c2);
	}
      nlines = c2;
      fclose(f2);
    
      c2 = 0;
      f2 = fopen(fname, "r");
      while (!feof(f2))
	{
	  //printf("reading c2=%d\n", c2);
	  fscanf(f2, "%lf %d ", &time, &NQ);
  	  for (i = 0; i < NQ; i++)
	    fscanf(f2, "(%lf,%lf) ", &rhoRt[c2][i], &rhoIt[c2][i]);
	  //printf("dopo\n");
	  if (c2 < MAXPTS && ti[c2] == -1.0)
	    {
	      ti[c2] = time;
	      //printf("c2=%d time=%.15G\n", c2, ti[c2]);
	    }
	  //printf("%d fname: %s %.15G %.15G %.15G\n", c2, fname, P0[0], P0[1], P0[2]);
	  c2++;
	}
      fclose(f2);
      for (nr1 = 0; nr1 < nlines; nr1++)
	{	
	  for (i=0; i < NQ; i++)
	    {
	      rhoR0[i] = rhoRt[i][nr1];
	      rhoI0[i] = rhoIt[i][nr1]; 
	    }
	  for (nr2 = nr1; nr2-nr1 < points && nr2 < nlines; nr2++)
	    {
	      for (i=0; i < NQ; i++) 
		{
		  rhoR1[i] = rhoRt[i][nr2];
		  rhoI1[i] = rhoIt[i][nr2];
		  C[i][nr2-nr1] += rhoR1[i]*rhoR0[i]+rhoI1[i]*rhoI0[i];
		  cc[i][nr2-nr1] += 1.0;
		  //printf("qui c3-c2=%d\n", c3-c2);
		}
	    }
	}
      for (i=0; i < NQ; i++)
	{
	  free(rhoIt[i]);
	  free(rhoRt[i]);
	}
	  
    }
  for (nq = 0; nq < NUMQ; nq++)
    {
      sprintf(fname, "sqt.%2d.k=%3d",ll*10+mm,nq);
      f = fopen(fname, "w+");
      for (ii=1; ii < points; ii++)
	{
	  Cav = 0.0;
	  for (i = 0; i < NQarr[nq]; i++)
	    if (cc[i][ii] > 0)
	  Cav += C[i][ii]/cc[i][ii];
	  if (ti[ii] > -1.0)
	    printf("%.15G %.15G %f\n", ti[ii]-ti[0], Cav, cc[0][ii]);
	}
      fclose(f);
    }
  return 0;
}
