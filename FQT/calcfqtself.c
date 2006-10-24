#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 1000
#define Sqr(x) ((x)*(x))
char **fname; 
double time, *ti, *r0[3], *r1[3], L, refTime;
int points=-1, assez, NP, NPA;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0, A1, B0, B1, C0, C1, storerate=-1.0;
int bakSaveMode = -1, eventDriven=0;
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3])
{
  FILE *f;
  double r0, r1, r2, R[3][3];
  int curstp, nat=0, i, cpos, a;
  double dt;
  *ti = -1.0;
  if (!(f = fopen(fname, "r")))
    {
      printf("ERROR: I can not open file %s\n", fname);
      exit(-1);
    }
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
  if (*ti == -1.0)
    *ti = ((double)curstp)*dt;
  fclose(f);
}
#define KMODMAX 599
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double *sqReA[KMODMAX], *sqImA[KMODMAX], *sqReB[KMODMAX], *sqImB[KMODMAX];
double *cc[KMODMAX];
char fname2[512];
char inputfile[1024];
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
void print_usage(void)
{
  printf("calcfqtself [ --qminpu/-qpum | --qmaxpu/-qpuM | --qmin/-qm <qmin> | --qmax/qM <qmax> |--help/-h] <lista_files> [points] [qmin] [qmax]\n");
  printf("where points is the number of points of the correlation function\n");
  exit(0);
}
double qavg[KMODMAX];

int qmin = 0, qmax = KMODMAX-1; 
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
      else if (extraparam == 0)
	{ 
	  extraparam++;
	  //printf("qui1 extraparam:%d\n", extraparam);
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam == 1)
	{
	  extraparam++;
	  //printf("qui2 argv[%d]:%s\n",cc, argv[cc]);
	  points = atoi(argv[cc]);
	  //printf("points:%d\n", points);
	}
      else if (extraparam == 2)
	{
	  extraparam++;
	  qmin = atoi(argv[cc]);
	  //printf("qmin=%d\n", qmin);
	}
      else if (extraparam == 3)
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
  FILE *f, *f2;
  int first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a;
  int iq, NN, fine, JJ, maxl, nfiles, nat, np, maxnp;
  int qmod; 
  double invL, rxdummy, sumImA, sumReA, sumImB, sumReB, scalFact;

  twopi = acos(0)*4.0;	  
#if 0
  if (argc <= 1)
    {
      printf("Usage: calcfqself <lista_file> [points] [qmin] [qmax]\n");
      printf("where points is the number of points of the correlation function\n");
      exit(-1);
    }
#endif 
  parse_param(argc, argv);
  //printf(">>qmin=%d\n", qmin);
  c2 = 0;
  if (!(f2 = fopen(inputfile, "r")))
    {
      printf("ERROR: I can not open file %s\n", inputfile);
      exit(-1);
    }
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


  if (!(f = fopen(fname[0], "r")))
    {
      printf("ERROR: I can not open file %s\n", fname[0]);
      exit(-1);
    }
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
      else if (!strcmp(parname,"storerate"))
	{
	  storerate = atof(parval);
	  eventDriven = 1;
	}
      else if (!strcmp(parname,"bakSaveMode"))
	bakSaveMode = atoi(parval);
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
#if 0
  if (argc >= 3)
    points = atoi(argv[2]);
  else
    points = NN;
  if (argc >= 4)
    qmin = atoi(argv[3]);
  if (argc == 5)
    qmax = atoi(argv[4]);
#endif
  if (points == -1)
    points = NN;
  scalFact = twopi * invL;
  set_qmin_qmax_from_q_pu(scalFact);
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
  if ((eventDriven==1 && storerate <= 0.0 && bakSaveMode <= 0)
      || (eventDriven==0 && bakSaveMode <= 0)) 
    NN = 1;
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;

  printf("qmin=%d qmax=%d invL=%.15G\n", qmin, qmax, invL);
  //printf("maxnp=%d points=%d\n",maxnp, points);
  if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
    assez = 0;
  else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
    assez = 1;
  else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
    assez = 2;
  if (NPA == -1)
    NPA = NP;
  if (eventDriven)
    printf("[ED] Event-Driven simulation\n");
  else
    printf("[MD] Time-Driven simulation\n");
 
  //fprintf(stderr, "allocating %d items NN=%d NP=%d num files=%d maxnp=%d\n", points, NN, NP, nfiles, maxnp);
  for (qmod = qmin; qmod <= qmax; qmod++)
    {
      sqReA[qmod] = malloc(sizeof(double)*points);
      sqImA[qmod] = malloc(sizeof(double)*points);
      cc[qmod] = malloc(sizeof(double)*points);
      if (NPA < NP)
	{
	  sqReB[qmod] = malloc(sizeof(double)*points);
	  sqImB[qmod] = malloc(sizeof(double)*points);
	}
      //ccB[qmod] = malloc(sizeof(double)*points);
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
	  sqReA[qmod][ii] = 0.0;
	  sqImA[qmod][ii] = 0.0;
	  cc[qmod][ii] = 0.0;
	  if (NPA < NP)
	    {
	      sqReB[qmod][ii] = 0.0;
	      sqImB[qmod][ii] = 0.0;
	    }
	  //ccB[qmod][ii] = 0.0;

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
		      sumReA=sumReB=0.0;
		      sumImA=sumImB=0.0;
		      for (i=0; i < NP; i++)
			{
			  rxdummy = scalFact*((r0[0][i]-r1[0][i])*mesh[qmod][iq][0]
			    +(r0[1][i]-r1[1][i])*mesh[qmod][iq][1]
			    +(r0[2][i]-r1[2][i])*mesh[qmod][iq][2]);
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
	    	      sqReA[qmod][np] += sumReA;
    		      sqImA[qmod][np] += sumImA;
		      sqReB[qmod][np] += sumReB;
		      sqImB[qmod][np] += sumImB;
		      cc[qmod][np] += 1.0;
		    }
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
	  //printf("cc[%d][%d]:%.15G\n", qmod, ii, cc[qmod][ii]);
	  sqReA[qmod][ii] = sqReA[qmod][ii]/cc[qmod][ii];
	  sqImA[qmod][ii] = sqImA[qmod][ii]/cc[qmod][ii];
	  if (NPA  < NP)
	    {
	      sqReB[qmod][ii] = sqReB[qmod][ii]/cc[qmod][ii];
	      sqImB[qmod][ii] = sqImB[qmod][ii]/cc[qmod][ii];
	    }
	 // printf("qmod=%d ii=%d sqre=%.15G sqIm=%.15G cc=%.15G\n",  qmod, ii, sqRe[qmod][ii], sqIm[qmod][ii],
	//	 cc[qmod][ii]);
	}
    }
  for (qmod = qmin; qmod <= qmax; qmod++)
    {
      sprintf(fname2, "Fqs-%d",qmod);
      if (!(f = fopen (fname2, "w+")))
	{
	  printf("ERROR: I can not open file %s\n", fname2);
	  exit(-1);
	}
      //printf("SqReA[%d][0]:%.15G SqReB[][0]: %.15G\n", qmod, sqReA[qmod][0], sqReB[qmod][0]);
      for (ii = 1; ii < points; ii++)
	{
	  if (NPA == NP)
	    {
	      if ((sqReA[qmod][ii]!=0.0 || sqImA[qmod][ii]!=0.0) && (ti[ii]> -1.0))
		fprintf(f, "%15G %.15G %.15G\n", ti[ii]-ti[0], sqReA[qmod][ii]/sqReA[qmod][0],
			sqImA[qmod][ii]);
	    }
	  else
	    {
	      if ((sqReA[qmod][ii]!=0.0 || sqImA[qmod][ii]!=0.0) && (ti[ii]> -1.0))
		fprintf(f, "%15G %.15G %.15G %.15G %.15G\n", ti[ii]-ti[0], sqReA[qmod][ii]/sqReA[qmod][0],
			sqReB[qmod][ii]/sqReB[qmod][0], 
			(sqReA[qmod][ii]+sqReB[qmod][ii])/(sqReA[qmod][0]+sqReB[qmod][0]), 
			sqImA[qmod][ii]+sqImB[qmod][ii]);
	    }
	}
      fclose(f);
    }
  return 0;
}
