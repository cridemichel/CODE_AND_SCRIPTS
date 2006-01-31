#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
char line[100000], parname[124], parval[256000];
int N;
double x[3], *r[3];
char fname[1024], inputfile[1024];
int readCnf = 0;
#define KMODMAX 99
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
double Sq[KMODMAX], sumRho, reRho, imRho, rCMk, scalFact, invNm, invL, L;
void print_usage(void)
{
  printf("calcSq <confs_file>\n");
  exit(0);
}
void parse_param(int argc, char** argv)
{
  int cc=1;
  
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
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
}
double twopi;
char dummy[2048];
int main(int argc, char** argv)
{
  FILE *f, *f2, *of;
  int nf, i, a, b, n, mp;
  double ti, tref=0.0, kbeg=0.0;
  int qmod, first = 1;
#if 0
  if (argc == 1)
    {
      printf("file with all confs to read as input\n");
      exit(-1);
    }
#endif
  twopi = acos(0.0)*4.0;
  parse_param(argc, argv);
  f2 = fopen(inputfile,"r");
  nf = 0;
  for (qmod=0; qmod < KMODMAX; qmod++)
    Sq[qmod] = 0.0;      

  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      f = fopen(fname,"r");
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
		N = atoi(parval);
	    }
	  while (strcmp(line,"@@@"));
	  for (i=0; i < 2*N; i++)
	    {
	      fscanf(f, "%[^\n]\n", line); 
	    }	
	  fscanf(f, "%lf\n", &L);
	  invL = 1.0/L;
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
	    N = atoi(parval);
	  if (!strcmp(parname,"time"))
	    ti = atof(parval);
  	}
      while (strcmp(line,"@@@"));
      //printf("fname=%s %d ellipsoids...\n", fname, N);

      if (first)
	{
	  for (a = 0; a < 3; a++)
	    {
	      r[a] = malloc(sizeof(double)*N);
	    }
	  first = 0;
	}
      for (i=0; i < N; i++)
	{
	   fscanf(f, "%[^\n]\n", line); 
	   if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
	     {
	       sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
	     }
	   scalFact = twopi * invL;
	   invNm = 1.0 / ((double)N);
	   printf("invL=%.15G invNm:%.15G\n", invL, invNm);
	   for(n = 0; n < KMODMAX; n++)
	     {
	       sumRho = 0.0;
	       for(mp = 0; mp < ntripl[n]; mp++)
		 {
		   reRho = 0.0;
		   imRho = 0.0;
		   for(i=0; i < N; i++)
		     {
		       /* il passo della mesh e' 0.5*pi2/L */
		       if (mesh[n][mp][0]==0 && mesh[n][mp][1] == 0 && 
			   mesh[n][mp][2] == 0)
			 {
			   printf("ERRORE nella MESH!!!!!!!! n=%d mp=%d ntripl[n]:%d\n", n,
				  mp, ntripl[n]);
			   exit(-1);
			 }
		       rCMk = kbeg + scalFact * 
			 (r[0][i] * mesh[n][mp][0] + r[1][i] * mesh[n][mp][1] + 
			  r[2][i] * mesh[n][mp][2]);
		       reRho = reRho + cos(rCMk) ; 
		       imRho = imRho + sin(rCMk);
		       /* Imaginary part of exp(i*k*r) for the actual molecule*/
		     }
		   sumRho = sumRho + Sqr(reRho) + Sqr(imRho);
		 }
	       Sq[n] += sumRho;  
	     }
	}
      fclose(f);
    }
  fclose(f2); 
  of = fopen("Sq.dat", "w+");
  for (qmod = 0; qmod  < KMODMAX; qmod++)
    {
      Sq[qmod] = (Sq[qmod]  * invNm) / ((double) ntripl[qmod] / ((double)nf));  
      fprintf(of, "%d %.15G\n", qmod, Sq[qmod]); 
    }
  fclose(of);
}
