#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
char line[1000000], parname[124], parval[2560000];
int N, NA=-1;
double x[3];
char fname[1024], inputfile[1024];
int readCnf = 0, physunit=0;
int  eventDriven = 0;
double scalFact, invNm, invL, L, invNmAA, invNmBB, invNmAB;
void print_usage(void)
{
  printf("calcT [ --help/-h | --cnf/-c ] <confs_file>\n");
  exit(0);
}
void parse_param(int argc, char** argv)
{
  int cc=1;
  int extraparam = 0;  

  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"--cnf") || !strcmp(argv[cc],"-c" ))
	{
	  readCnf = 1;
	} 
      else if (cc == argc || extraparam == 1)
	print_usage();
      else if (extraparam == 0)
	{ 
	  extraparam = 1;
	  strcpy(inputfile,argv[cc]);
	}
      else
	print_usage();
      cc++;
    }
}
double twopi;
char dummy[2048];
int main(int argc, char** argv)
{
  FILE *f, *f2, *of;
  int nf, i, a, b, n, mp;
  double tiTD, dt, Tist, Tavg=0.0, ti=-1.0, tref=0.0, Vol, a1, a2, a3, a4;
  int first = 1, NP1, NP2;
#if 0
  if (argc == 1)
    {
      printf("file with all confs to read as input\n");
      exit(-1);
    }
#endif
  parse_param(argc, argv);
  if (!(f2 = fopen(inputfile,"r")))
    {
      printf("ERROR: I can not open file %s\n", inputfile);
      exit(-1);
    }
  nf = 0;
  of = fopen("temperature.dat", "w");
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      if (!(f = fopen(fname,"r")))
	{
	  printf("ERROR: I can not open file %s\n", fname);
	  exit(-1);
	}
      nf++;
      tref=0.0;
      if (first)
	{
	  /* legge L */
	  if (readCnf)
	    {
	      do 
		{
		  fscanf(f,"%[^\n]\n",line);
		}
	      while (strcmp(line,"@@@"));
	    }
	  do 
	    {
	      fscanf(f,"%[^\n]\n",line);
	      sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	      if (!strcmp(parname,"parnum"))
		{
	       	  if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
	    	    {
	    	      N = atoi(parval);
	    	    }
		  else
		    {
		      N = NP1+NP2;
		      NA = NP1;
		    }
		}
	      else if (!strcmp(parname, "time"))
		eventDriven = 1;
	      else if (!strcmp(parname, "steplength"))
		dt = atof(parval);
	    }
	  while (strcmp(line,"@@@"));
	  for (i=0; i < 2*N; i++)
	    {
	      fscanf(f, "%[^\n]\n", line); 
	    }	
 	  if (eventDriven == 1)
	    {
	      fscanf(f, "%lf\n", &L);
	      invL = 1.0/L;
	      //printf("qui L=%.15G\n", L);
	    }
	  else
	    {
	      fscanf(f, "%lf %lf %lf %lf %lf\n", &Vol, &a1, &a2, &a3, &a4);
	      L = cbrt(Vol);
	      invL = 1.0/L;
	    }
	}
      //printf("L=%.15G\n", L);
      /* -------- */
      rewind(f);
      if (readCnf)
	{
	  do 
	    {
	      fscanf(f,"%[^\n]\n",line);
	      sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	      if (!strcmp(parname,"refTime"))
		tref = atof(parval);
	    }
	  while (strcmp(line,"@@@"));
	}
      do 
	{
	  fscanf(f,"%[^\n]\n",line);
	  sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	  if (!strcmp(parname,"parnum"))
	    {
	      if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
		{
		  N = atoi(parval);
		}
	      else
		{
		  N = NP1+NP2;
		  NA = NP1;
		}
	    }
	  if (!strcmp(parname,"parnumA"))
	    NA = atoi(parval); 
	  if (!strcmp(parname,"time"))
	    ti = atof(parval);
	  if (!strcmp(parname,"curStep"))
	    tiTD = atof(parval)*dt;
  	}
      while (strcmp(line,"@@@"));
      //printf("fname=%s %d ellipsoids...\n", fname, N);

      for (i=0; i < N; i++)
	{
	   fscanf(f, "%[^\n]\n", line); 
	   //printf("line=%s\n", line);
	   if (!sscanf(line, "%lf %lf %lf\n", &a1, &a2, &a3)==3)
	     {
	       //printf("boh\n");
	       sscanf(line, "%lf %lf %lf %[^\n]\n", &a1, &a2, &a3, dummy); 
	     }
	   //printf("r=(%.15G,%.15G,%.15G)\n", r[0][i], r[1][i], r[2][i]);
	}
      Tist = 0.0;
      for (i=0; i < N; i++)
	{
	   fscanf(f, "%[^\n]\n", line); 
	   //printf("line=%s\n", line);
	   if (!sscanf(line, "%lf %lf %lf\n", &a1, &a2, &a3)==3)
	     {
	       //printf("boh\n");
	       sscanf(line, "%lf %lf %lf %[^\n]\n", &a1, &a2, &a3, dummy); 
	     }
	   Tist += Sqr(a1)+Sqr(a2)+Sqr(a3);
	}
      Tist /= 3.0*(N-3.0); 
      Tavg += Tist;
      fprintf(of, "%.15G %.15G\n", ((ti > 0.0)?ti:tiTD), Tist);
      if (NA == -1)
	NA = N;
      //NA = N;//force monodisperse
      fclose(f);
    }
  fclose(of);
  fclose(f2); 
  printf("Average Temperature=%.15G\n", Tavg/((double)nf));
  return 0;
}
