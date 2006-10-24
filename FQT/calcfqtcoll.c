#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 1000
#define MAXQ 500
#define NUMQ 100
#define MAXCOMPS 2
#define KMODMAX 600
#define NKSHELL 150
#define Sqr(x) ((x)*(x))
char fname[2][1024]; 
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"

double time, Cav, Cav0, CavAA0, CavBB0, CavAB0, 
       CavAA, CavAB, CavBB, CavBA0, CavBA, 
       rhoR0[MAXCOMPS][MAXQ], rhoI0[MAXCOMPS][MAXQ], 
       rhoR0[MAXCOMPS][MAXQ], rhoI0[MAXCOMPS][MAXQ],
       *cc[MAXQ],*CAA[MAXQ],
       *CBB[MAXQ], *CAB[MAXQ], *CBA[MAXQ],
       rhoR1[MAXCOMPS][MAXQ], rhoI1[MAXCOMPS][MAXQ], *ti, *rhoRt[MAXCOMPS][MAXQ], 
       *rhoIt[MAXCOMPS][MAXQ],
       //*rhoRtp[MAXCOMPS][MAXQ], *rhoItp[MAXCOMPS][MAXQ], 
       *rhoRt[MAXCOMPS][MAXQ], 
       *rhoIt[MAXCOMPS][MAXQ], *tiall;
int *NQarr;
int qmin, qmax, points=-1, comps=2;
char AB[2]={'A','B'};
void print_usage(void)
{
  printf("calcfqtcoll [--qmin/-qm <qmin> | --qmax/qM <qmax> |--help/-h] <lista_files> [qmin] [qmax] [points]\n");
  printf("where points is the number of points of the correlation function\n");
  exit(0);
}
double qavg[KMODMAX];

double qminpu = -1.0, qmaxpu = -1.0;

void parse_param(int argc, char** argv)
{
  int cc=1, extraparam=0;
  
  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      //printf("cc=%d extraparam=%d argc=%d\n", cc, extraparam, argc);
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"--ncomps") || !strcmp(argv[cc],"-nc"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  comps = atoi(argv[cc]);
	  if (comps < 1 || comps > 2)
	    print_usage();
	}
      else if (!strcmp(argv[cc],"--qmin") || !strcmp(argv[cc],"-qm"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  qmin = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--qmax") || !strcmp(argv[cc],"-qM"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  qmax = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--qpumin") || !strcmp(argv[cc],"-qpum"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  qminpu = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--qpumax") || !strcmp(argv[cc],"-qpuM"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  qmaxpu = atof(argv[cc]);
	}
      else if (cc == argc || extraparam == 4)
	print_usage();
      else if (extraparam == 2)
	{
	  extraparam++;
	  //printf("qui2 argv[%d]:%s\n",cc, argv[cc]);
	  points = atoi(argv[cc]);
	  //printf("points:%d\n", points);
	}
      else if (extraparam == 0)
	{
	  extraparam++;
	  qmin = atoi(argv[cc]);
	  //printf("qmin=%d\n", qmin);
	}
      else if (extraparam == 1)
	{
	  extraparam++;
	  qmax = atoi(argv[cc]);
	}
      else
	print_usage();
      cc++;
    }
}

void set_qmin_qmax_from_q_pu(double scalFact)
{
  int qmod, mp;
  for (qmod = 0; qmod < KMODMAX; qmod++)
    {
      qavg[qmod] = 0;
      for (mp = 0; mp < ntripl[qmod]; mp++) 
	{
	  qavg[qmod] += sqrt(Sqr(((double)mesh[qmod][mp][0]))+
			     Sqr(((double)mesh[qmod][mp][1]))+
			     Sqr(((double)mesh[qmod][mp][2])));
	}
      qavg[qmod] *= scalFact/((double)ntripl[qmod]);
      //printf("qavg[%d]:%.15G - %.15G\n", qmod, qavg[qmod], scalFact*(1.25+0.5*qmod));
    }
  if (qminpu != -1.0 && qminpu == qmaxpu)
    {
      qmin = rint(qminpu / (scalFact/2.0) - 2.0);
      qmax = qmin;
    }
  else 
    {
      if (qminpu != -1.0 && qminpu >= 0.0)
	qmin = ceil( qminpu / (scalFact/2.0) - 2.0);
      if (qmaxpu != -1.0 && qmaxpu > 0.0)
	qmax = floor( qmaxpu / (scalFact/2.0) - 2.0);
    }
}


int main(int argc, char **argv)
{
  FILE *f, *f2, *f3, *f4;
  double A1, A2, A3 ;
  int c, first=1, firstp=1, NQ, nq, c1, c2, c3, i, ii, nlines, nr1, nr2, ll, mm, llp, mmp;
  int NN, fine, JJ, maxnp, np;
#if 0
  if (argc <= 1)
    {
      printf("Usage: calcfqtcoll <qmin> <qmax> [points] \n");
      //printf("where NN is the number of configurations in a logarithmic block\n");
      exit(-1);
    }
#endif
  parse_param(argc, argv);
  //NN =  atoi(argv[1]);
#if 0
  qmin = atoi(argv[1]);
  qmax = atoi(argv[2]);
#endif
  NQarr = malloc(sizeof(int)*qmax);
  if (qmax >= KMODMAX)
    qmax = KMODMAX-1;
  if (qmax < 0)
    qmax = 0;
  if (qmin < 0)
    qmin = 0;
  if (qmin >= KMODMAX)
    qmin = KMODMAX-1;
  
  if (qmin > qmax)
    {
      printf("ERROR: qmin must be less than qmax\n");
      exit(-1);
    }

  if (points==-1)
    points = NN;
#if 0
  if (argc == 4)
    points = atoi(argv[3]);
  else
    points = NN;
#endif
  printf("points:%d\n", points);
  c2 = 0;
  if (comps==2)
    sprintf(fname[0], "RHOTMPA/ro.00.k=%03d", qmin);
  else
    sprintf(fname[0], "RHOTMP/ro.00.k=%03d", qmin);
  f2 = fopen(fname[0], "r");
  while (!feof(f2))
    {
      fscanf(f2, "%lf %d\n", &time, &NQ);
      //printf("time=%f nq=%d\n", time, NQ);
      for (i = 0; i < NQ; i++)
	fscanf(f2, "(%lf,%lf) ", &A1, &A2);
      c2++;
    }
  if (comps==2)
    sprintf(fname[0], "RHOTMPA/NN.dat");
  else
    sprintf(fname[0], "RHOTMP/NN.dat");
  if (!(f = fopen(fname[0],"r")))
    {
      printf("ERROR: file %s does not exist\n", fname[0]);
      exit(-1);
    }
  fscanf(f,"%d\n", &NN);
  fclose(f);
  //printf("NN=%d\n", NN);
  maxnp = NN + (c2-NN)/NN;
  nlines = c2;
  if (points > maxnp)
    points = maxnp;
  //printf("argv[3]:%s points=%d maxnp=%d c2=%d NN=%d\n", argv[3], points, maxnp, c2, NN);
  fprintf(stderr, "allocating %d items points=%d NN=%d\n", nlines, points, NN);
#if 0
  if ((f4 = fopen("RHOTMP/components.dat","r"))==NULL)
    comps = 1; /* monodisperse system (one component) */
  else
    {
      fscanf(f4, "%d\n", &comps);
      fclose(f4);
    }
#endif 
  for (i=0; i < MAXQ; i++)
    {
      for ( c = 0; c < 2; c++)
	{
	  rhoRt[c][i] = malloc(sizeof(double)*nlines);
	  rhoIt[c][i] = malloc(sizeof(double)*nlines);
#if 0
	  rhoRtp[c][i] = malloc(sizeof(double)*nlines);
	  rhoItp[c][i] = malloc(sizeof(double)*nlines);
#endif
	}
      cc[i] = malloc(sizeof(double)*points);
      CAA[i]  = malloc(sizeof(double)*points);
      if (comps==2)
	{
	  CBB[i]  = malloc(sizeof(double)*points);
	  CAB[i]  = malloc(sizeof(double)*points);
	  CBA[i]  = malloc(sizeof(double)*points);
	}
    }
  ti = malloc(sizeof(double)*points);
  tiall = malloc(sizeof(double)*nlines);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;
  for (ii=0; ii < nlines; ii++)
    tiall[ii] = -1.0;
  first = 0;
  fclose(f2);

  for (nq = qmin; nq <= qmax; nq++)
    {
      for (i=0; i < MAXQ; i++)
	for (ii=0; ii < points; ii++)
	  {
	    CAA[i][ii] = CBB[i][ii] = CAB[i][ii] = CBA[i][ii] = 0.0;
	    cc[i][ii] = 0;
	  }
      fprintf(stderr, "nq=%d\n", nq);
      for (c = 0; c < 2; c++)
	{
	  c2 = 0;
	  sprintf(fname[0], "RHOTMP%c/ro.00.k=%03d", AB[c],nq);
	  //printf("Opening fname=%s\n", fname[0]);
	  f2 = fopen(fname[0], "r");
	  while (!feof(f2))
	    {
	      //printf("reading c2=%d f2=%p feof(f2):%d\n", c2, f2, feof(f2));
	      fscanf(f2, "%lf %d\n", &time, &NQ);
	      //printf("time=%.15G c2=%d\n", time, c2);
	      if (c2==0)
		NQarr[nq] = NQ;
	      
	      for (i = 0; i < NQ; i++)
		{
		  fscanf(f2, "(%lf,%lf) ", &(rhoRt[c][i][c2]), &(rhoIt[c][i][c2]));
		}
	      //printf("qui1\n");
	      if (tiall[c2] == -1.0)
		{
		  tiall[c2] = time;
		  //printf("c2=%d time=%.15G\n", c2, ti[c2]);
		}
	      //printf("qui2\n");
	      //printf("%d fname: %s %.15G %.15G %.15G\n", c2, fname, P0[0], P0[1], P0[2]);
	      c2++;
	    }
	  fclose(f2);
	  //c2 = 0;
	}
      for (nr1 = 0; nr1 < nlines; nr1=nr1+NN)
	{	
	  for (i=0; i < NQ; i++)
	    {
	      for (c = 0; c < 2; c++)
		{
		  rhoR0[c][i] = rhoRt[c][i][nr1];
		  rhoI0[c][i] = rhoIt[c][i][nr1]; 
		}
	    }

	  fine = 0;
	  for (JJ = 0; fine == 0; JJ++)
	    {
	      for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
		{
		  //printf("qui nr1=%d nr2=%d NN=%d nlines=%d points=%d\n", nr1, nr2, NN, nlines, points);
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
		      for (c=0; c < 2; c++)
			{
			  rhoR1[c][i] = rhoRt[c][i][nr2];
			  rhoI1[c][i] = rhoIt[c][i][nr2];
			}
		      CAA[i][np] += rhoR1[0][i]*rhoR0[0][i]+rhoI1[0][i]*rhoI0[0][i];
		      //printf("CAA[%d][%d]=%.15G\n", i, np, CAA[i][np]);
		      if (comps==2)
			{
			  CBB[i][np] += rhoR1[1][i]*rhoR0[1][i]+rhoI1[1][i]*rhoI0[1][i];
			  CAB[i][np] += rhoR1[0][i]*rhoR0[1][i]+rhoI1[0][i]*rhoI0[1][i];
			  /* N.B.: check this!!! */
			  CBA[i][np] += rhoR1[1][i]*rhoR0[0][i]+rhoI1[1][i]*rhoI0[0][i];
			}
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
      sprintf(fname[0], "sqt.k=%03d",nq);
      f = fopen(fname[0], "w+");
      sprintf(fname[0], "N-sqt.k=%03d", nq);
      f2 = fopen(fname[0], "w+");
      Cav = CavAA0 = CavBB0 = CavAB0 = CavBA0 = 0.0;
      for (i = 0; i < NQarr[nq]; i++)
	{
	  if (cc[i][0] > 0)
	    {
	      CavAA0 += CAA[i][0]/cc[i][0];
	      if (comps == 2)
		{
		  CavBB0 += CBB[i][0]/cc[i][0];
		  CavAB0 += CAB[i][0]/cc[i][0];
		  CavBA0 += CBA[i][0]/cc[i][0];
		}
	    }
	}
      CavAA0 /= ((double)NQarr[nq]);
      //printf("CavAA0=%.15G NQarr[%d]=%d\n", CavAA0, nq, NQarr[nq]);
      if (comps==2)
	{
	  CavBB0 /= ((double)NQarr[nq]);
	  CavAB0 /= ((double)NQarr[nq]);
	  CavBA0 /= ((double)NQarr[nq]);
	  Cav0 = CavAA0 + CavAB0 + CavBB0 + CavBA0;
	}
      for (ii=1; ii < points; ii++)
	{
	  CavAA = 0.0;
	  if (comps==2)
	    CavBB = CavAB = 0.0;
	  for (i = 0; i < NQarr[nq]; i++)
	    {
	      if (cc[i][ii] > 0)
		{
		  CavAA += CAA[i][ii]/cc[i][ii];
		  if (comps==2)
		    {
		      CavBB += CBB[i][ii]/cc[i][ii];
		      CavAB += CAB[i][ii]/cc[i][ii];
		      CavBA += CBA[i][ii]/cc[i][ii];
		    }
		}
	    }
	  CavAA /= ((double)NQarr[nq]);
	  if (comps==2)
	    {
	      CavBB /= ((double)NQarr[nq]);
	      CavAB /= ((double)NQarr[nq]);
	      CavBA /= ((double)NQarr[nq]);
	      Cav = CavAA + CavAB + CavBB + CavBA;
	    }
	  if (ti[ii] > -1.0)
	    {
	    if (comps==2)
	      {
		fprintf(f, "%.15G %.15G %.15G %.15G %.15G %f\n", ti[ii]-ti[0], Cav, CavAA, 
		   	CavBB, CavAB+CavBA, cc[0][ii]*NQarr[nq]);
		fprintf(f2, "%.15G %.15G %.15G %.15G %.15G %f\n", ti[ii]-ti[0], Cav/Cav0, 
			CavAA/CavAA0, CavBB/CavBB0, (CavAB+CavBA)/(CavAB0+CavBA0),
		     	cc[0][ii]*NQarr[nq]);

	      }
	    else
	      {
		fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], CavAA, cc[0][ii]*NQarr[nq]);
		fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], CavAA/CavAA0, cc[0][ii]*NQarr[nq]);
	      }
	    }
	}
      fclose(f);
      fclose(f2);
    }
  return 0;
}
