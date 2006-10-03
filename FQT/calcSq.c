#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
char line[100000], parname[124], parval[256000];
int N, NA=-1;
double x[3], *r[3];
char fname[1024], inputfile[1024];
int readCnf = 0;
#define KMODMAX 98
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
double Sq[KMODMAX], sumRho, reRho, imRho, rCMk, scalFact, invNm, invL, L, invNmAA, invNmBB, invNmAB;
double SqAA[KMODMAX], SqBB[KMODMAX], SqAB[KMODMAX], sumRhoAA, sumRhoAB, sumRhoBB, reRhoA, reRhoB, imRhoA, 
       imRhoB;
void print_usage(void)
{
  printf("calcSq [--help/-h] [--cnf/-c] <confs_file>\n");
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
	  if (!strcmp(parname,"parnumA"))
	    NA = atoi(parval); 
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
	}
      for (qmod=0; qmod < KMODMAX; qmod++)
	{
	  Sq[qmod] = 0.0;      
	  if (NA!=-1)
	    {
	      SqAA[qmod] = SqAB[qmod] = SqBB[qmod] = 0.0;
	    }
	}
  
      scalFact = twopi * invL;
      invNm = 1.0 / ((double)N);
      if (NA != -1)
	{
	  invNmAA = 1.0/((double)NA);
	  invNmBB = 1.0/((double)N-NA);
	  invNmAB = 1.0/sqrt(((double)N-NA)*((double)NA));
	}
      
      if (NA!=-1)
	printf("[MIXTURE N=%d NA=%d] ", N, NA);
      else 
	printf("[MONODISPERSE] ");
      printf("nf=%d twopi=%.15G N=%d invL=%.15G invNm:%.15G\n", nf, twopi, N, invL, invNm);
      if (NA == -1)
	{
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
		      reRho = reRho + cos(rCMk); 
		      imRho = imRho + sin(rCMk);
		      /* Imaginary part of exp(i*k*r) for the actual molecule*/
		    }
		  sumRho = sumRho + Sqr(reRho) + Sqr(imRho);
		}
	      Sq[n] += sumRho;  
	      //printf("sumRho=%.15G Sq[%d]=%.15G\n", sumRho, n, Sq[n]);
	    }
	  fclose(f);
	}
      else
	{
	  for(n = 0; n < KMODMAX; n++)
	    {
	      sumRhoAA = sumRhoBB = sumRhoAB = 0.0;
	      for(mp = 0; mp < ntripl[n]; mp++)
		{
		  reRhoA = reRhoB = 0.0;
		  imRhoA = imRhoB = 0.0;
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
		      reRhoA = reRhoA + cos(rCMk); 
		      imRhoA = imRhoA + sin(rCMk);
		      if (NA!=-1)
			{
			  reRhoA = reRhoA + cos(rCMk); 
			  imRhoA = imRhoA + sin(rCMk);
			}

		      /* Imaginary part of exp(i*k*r) for the actual molecule*/
		    }
		  sumRhoAA = sumRhoAA + Sqr(reRhoA) + Sqr(imRhoA);
		  sumRhoBB = sumRhoBB + Sqr(reRhoB) + Sqr(imRhoB);
		  sumRhoAB = sumRhoAB + reRhoA*reRhoB + imRhoA*imRhoB;
		}
	      Sq[n] += sumRhoAA + sumRhoBB + 2.0*sumRhoAB;  
	      SqAA[n] += sumRhoAA;
	      SqBB[n] += sumRhoBB;
	      SqAB[n] += sumRhoAB;

	      //printf("sumRho=%.15G Sq[%d]=%.15G\n", sumRho, n, Sq[n]);
	    }
	  fclose(f);

	}
    }
  fclose(f2); 
  of = fopen("Sq.dat", "w+");
  for (qmod = 0; qmod  < KMODMAX; qmod++)
    {
      Sq[qmod] = (Sq[qmod]  * invNm) / ((double) ntripl[qmod]) / ((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      fprintf(of, "%d %.15G\n", qmod, Sq[qmod]); 
    }
  fclose(of);
  if (NA != -1) 
    {
      of = fopen("SqAA.dat", "w+");
      for (qmod = 0; qmod  < KMODMAX; qmod++)
	{
	  SqAA[qmod] = (SqAA[qmod]  * invNmAA) / ((double) ntripl[qmod]) / ((double)nf);  
	  //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
	  fprintf(of, "%d %.15G\n", qmod, SqAA[qmod]); 
	}
      fclose(of);
      of = fopen("SqBB.dat", "w+");
      for (qmod = 0; qmod  < KMODMAX; qmod++)
	{
	  SqBB[qmod] = (SqBB[qmod]  * invNmBB) / ((double) ntripl[qmod]) / ((double)nf);  
	  //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
	  fprintf(of, "%d %.15G\n", qmod, SqBB[qmod]); 
	}
      fclose(of);
      of = fopen("SqAB.dat", "w+");
      for (qmod = 0; qmod  < KMODMAX; qmod++)
	{
	  SqAB[qmod] = (SqAB[qmod]  * invNmAB) / ((double) ntripl[qmod]) / ((double)nf);  
	  //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
	  fprintf(of, "%d %.15G\n", qmod, SqAB[qmod]); 
	}
      fclose(of);

    }
}
