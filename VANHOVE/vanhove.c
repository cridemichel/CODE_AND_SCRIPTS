#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define COORD_TYPE double
#define MAXPTS 1000
#define Sqr(x) ((x)*(x))
char **fname; 
int *isPercPart;
double time, *ti, *r0[3], *r1[3], L, refTime;
int points=-1, assez, NP, NPA, clusters=0;
char parname[128], parval[25600000], line[25600000];
char dummy[2048000], cluststr[2048000];
char *pnum;
double A0, A1, B0, B1, C0, C1, storerate=-1.0, *ccav;
int bakSaveMode = -1, eventDriven=0, skip=0;
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
	      //printf("fn=%s r[%d]=(%f,%f,%f)\n", fname, i, r[0][i], r[1][i], r[2][i]);
	    }
 	  break; 
	}
    }
  if (*ti == -1.0)
    *ti = ((double)curstp)*dt;
  fclose(f);
}
char fname2[512];
char inputfile[1024];
void print_usage(void)
{
  printf("vanhove [ --deltar/-dr | --skip/-s | --rmax/-rM | --tmin/-tm | --tmax/-tM ] <lista_files> [points]\n");
  printf("where points is the number of points of the correlation function\n");
  exit(0);
}

double delr=0.01, rmax = 10.0;
int tmax=100, tmin=0; 
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
      else if (!strcmp(argv[cc],"--deltar") || !strcmp(argv[cc],"-dr"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  delr = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--skip") || !strcmp(argv[cc],"-s"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  skip = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--tmax") || !strcmp(argv[cc],"-tM"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  tmax = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--tmin") || !strcmp(argv[cc],"-tm"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  tmin = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--rmax") || !strcmp(argv[cc],"-rM"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  rmax = atof(argv[cc]);
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
      else
	print_usage();
      cc++;
    }
}
int Gsnr;
double Deltar;
char fn[1024];
void saveGself(char* fileName, COORD_TYPE** gs)
{
  FILE* afs;
  int j, t;
  COORD_TYPE r;

  /* save every vhgap steps */
  for(t = 0; t < tmax; t++)
    {
      if (t == 0) continue;
      sprintf(fn, "t-%f-%s", ti[t], fileName);
      afs = fopen(fn, "w+");
      for(j = 0; j < Gsnr; j++) /* Loop over angles */
	{ 
	  r = ((COORD_TYPE)j + 0.5) * (rmax /  Gsnr);/* In degree !!! */ 
	  fprintf(afs, "%.15G %.15G\n", r, gs[t][j]);
	}
      //fprintf(afs, "&\n"); /* This indicates to xmgr that begins a new set */
      fclose(afs);
    }

}
char GselfFile[]="Gself.dat";
char GselfFileA[]="GselfA.dat";
char GselfFileB[]="GselfB.dat";
double twopi, pi;
double **Gself, **GselfA, **GselfB;
int main(int argc, char **argv)
{
  FILE *f, *f2;
  int t, j, first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a;
  int iq, NN, fine, JJ, maxl, nfiles, nat, np, maxnp;
  int NP1, NP2, kk, isperc; 
  double Normv, invL, rxdummy, sumImA, sumReA, sumImB, sumReB, scalFact;
  double r, costmp, sintmp;
  twopi = acos(0)*4.0;	  
  pi = twopi*0.5;
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
	{
	  if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
	    {
	      NP = atoi(parval);
	    }
	  else
	    {
	      NP = NP1+NP2;
	      NPA = NP1;
	    }
	}
      else if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      else if (nat==0 && !strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (!strcmp(parname,"storerate"))
	{
	  storerate = atof(parval);
	  eventDriven = 1;
	}
      else if (!strcmp(parname,"bakSavedMode"))
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
  if (!eventDriven)
    L = cbrt(L);
  invL = 1.0/L;
  if (points == -1)
    points = NN;
  
  if ((eventDriven==1 && storerate <= 0.0 && bakSaveMode <= 0)
      || (eventDriven==0 && bakSaveMode <= 0)) 
    NN = 1;
  if (NN!=1)
    skip = 0;
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  printf("rmax=%f invL=%.15G NN=%d points=%d\n", rmax, invL, NN, points);
  //printf("maxnp=%d points=%d\n",maxnp, points);
  if (NPA == -1)
    NPA = NP;
  if (eventDriven)
    printf("[ED] Event-Driven simulation\n");
  else
    printf("[MD] Time-Driven simulation\n");
  if (NP==NPA)
    {
      printf("[MONODISPERSE] NP=%d\n", NP);
    }
  else
    {
      printf("[MIXTURE] NP=%d NPA=%d\n", NP, NPA);
    }
  Gsnr = rmax / delr;
  printf("Gsnr=%d rmax=%f delr=%f\n", Gsnr, rmax, delr);
  Gself = (double**) malloc(tmax*sizeof(double*));
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      r1[a] = malloc(sizeof(double)*NP);
    }

  if (NPA!=NP)
    {  
      GselfA = (double**) malloc(tmax*sizeof(double*));
      GselfB = (double**) malloc(tmax*sizeof(double*));
      ccav = (double*) malloc(tmax*sizeof(double));
      ti = (double*) malloc(tmax*sizeof(double));
    }
  for (j = 0; j < tmax; j++)
    {
      Gself[j] = (COORD_TYPE*) malloc(Gsnr*sizeof(double));
      if (NPA!=NP)
	{  
	  GselfA[j] = (COORD_TYPE*) malloc(Gsnr*sizeof(double));
	  GselfB[j] = (COORD_TYPE*) malloc(Gsnr*sizeof(double));
	}
    } 
  for (ii=0; ii < tmax; ii++)
    ti[ii] = -1.0;
  first = 0;
  fclose(f2);
  c2 = 0;
  JJ = 0;
  for (t=0; t < tmax; ++t)
    {
      ccav[t] = 0.0;
      for (j = 0; j < Gsnr; j++) 
	{
	  Gself[t][j] = 0.0;
	  if (NPA!=NP)
	    {
	      GselfA[t][j] = 0.0;
	      GselfB[t][j] = 0.0;
	    }
	}
    }
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN+skip)
    {	
      //printf("nr1=%d fname=%s nfiles=%d\n", nr1, fname[nr1], nfiles);
      readconf(fname[nr1], &time, &refTime, NP, r0);
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      /* N.B. considera NN punti in maniera logaritmica e poi calcola i punti in maniera lineare 
	       * distanziati di NN punti. */
              np = (JJ == 0)?nr2-nr1:NN-1+JJ;	      
	      //printf("nr2=%d fname=%s nfiles=%d\n", nr2, fname[nr2], nfiles);
	      if (nr2 >= nfiles || np >= tmax)
		{
		  fine = 1;
		  //printf("qui3\n");
		  break;
		}
	      if (JJ > 0 && (nr2 - nr1) % NN != 0)
		continue;
	      if (np >= tmax) 
		continue;
	      //printf("qui1\n");
	      readconf(fname[nr2], &time, &refTime, NP, r1);
      	      //printf("qui2\n");
	      if (np < tmax && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("np=%d time=%.15G\n", np, ti[np]);
		}
  
	      //if (nr2 == nr1)
		//continue;
	      t=np;
	      ccav[t] += 1.0;
	      //printf("qui B\n");
	      for (i=0; i < NP; i++)
		{
		  Deltar = sqrt(Sqr(r1[0][i] - r0[0][i]) + 
				Sqr(r1[1][i] - r0[1][i]) +
				Sqr(r1[2][i] - r0[2][i]));
		  j = (int) ((((double)Gsnr) * Deltar) / rmax); 
		  //printf("qui C t=%d j=%d i=%d\n", t, j, i);

		  if (j < Gsnr)
		    {
		      Gself[t][j] += 2.0;
		      if (NPA != NP)
			{
			  if (i < NPA)
			    {
			      GselfA[t][j] +=2.0;
			    }
			  else
			    {
			      GselfB[t][j] += 2.0;
			    }
			}
		    }
		}
	    }
	}
    }
  /* Normalization */
  for(t = 0; t < tmax; t++)
    {
      for (j = 0; j < Gsnr; j++)
	{
	  r = ((COORD_TYPE) j + 0.5) * delr;
	  Normv = ccav[t];
	  printf("Normv=%f\n", Normv);
	  Gself[t][j] /= delr ;/* Gs is a density of probability */
	  Gself[t][j] /= Normv * NP;
	  //Gself[t][j] /= Normv * Sqr(r) * 4.0 * pi * NP;
	  if (NPA!=NP)
	    {
	      GselfA[t][j] /= delr;/* Gs is a density of probability */
	      GselfA[t][j] /= Normv * NPA;
	      //GselfA[t][j] /= Normv * Sqr(r) * 4.0 * pi * NPA;

	      GselfB[t][j] /= delr;/* Gs is a density of probability */
	      GselfB[t][j] /= Normv * (NP-NPA);
	      //GselfB[t][j] /= Normv * Sqr(r) * 4.0 * pi * (NP-NPA);
	    }
	}
    }

  //printf("eccoci\n");
  /* save Gself */
  saveGself(GselfFile, Gself);
  if (NPA!=NP)
    {
      saveGself(GselfFileA, GselfA);
      saveGself(GselfFileB, GselfB);
    }
  return 0;
}
