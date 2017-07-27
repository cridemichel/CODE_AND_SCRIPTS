#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
const double patchdist=11.15;
char line[1000000], parname[124], parval[2560000];
int N, NA=-1;
double x[3], *r[3], *u[3], rx, ry, rz;
char fname[1024], inputfile[1024];
int readCnf = 0, physunit=0;
#define KMODMAX 599
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double qavg[KMODMAX];
int qmin=0, qmax=KMODMAX-1, eventDriven = 0, bigpatchsq=0;
double qminpu=-1.0, qmaxpu=-1.0;
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
double Sq[KMODMAX], sumRho, reRho, imRho, rCMk, scalFact[3], invNm, L[3], invNmAA, invNmBB, invNmAB;
double SqAA[KMODMAX], SqBB[KMODMAX], SqAB[KMODMAX], sumRhoAA, sumRhoAB, sumRhoBB, reRhoA, reRhoB, imRhoA, imRhoB;
void print_usage(void)
{
  printf("calcSq [ --bigpatchsq/-bp --qminpu/-qpum | --qmaxpu/-qpuM | --qmin/-qm <qmin> | --qmax/qM <qmax> |--help/-h | --cnf/-c | --phys-unit/-pu] <confs_file> [qmin] [qmax]\n");
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
      else if (!strcmp(argv[cc],"--phys-unit") || !strcmp(argv[cc],"-pu"))
	{
	  physunit = 1;
	}
      else if (!strcmp(argv[cc],"--bigpatchsq") || !strcmp(argv[cc],"-bp"))
	{
	  bigpatchsq = 1;
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
      else if (cc == argc || extraparam == 3)
	print_usage();
      else if (extraparam == 0)
	{ 
	  extraparam = 1;
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam == 1)
	{
	  extraparam = 2;
	  //printf("qui2 argv[%d]:%s\n",cc, argv[cc]);
	  qmin = atoi(argv[cc]);
	  //printf("qmin:%d\n", qmin);
	}
      else if (extraparam == 2)
	{
	  extraparam = 3;
	  qmax = atoi(argv[cc]);
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
  double qval, ti, tref=0.0, kbeg=0.0, Vol, a1, a2, a3, a4;
  int qmod, first = 1, NP1, NP2;
#if 0
  if (argc == 1)
    {
      printf("file with all confs to read as input\n");
      exit(-1);
    }
#endif
  twopi = acos(0.0)*4.0;
  parse_param(argc, argv);
  if (!(f2 = fopen(inputfile,"r")))
    {
      printf("ERROR: I can not open file %s\n", inputfile);
      exit(-1);
    }
  nf = 0;
  if (qmax >= KMODMAX)
    qmax = KMODMAX-1;
  if (qmin < 0)
    qmin = 0;
  if (qmin > qmax)
    {
      printf("ERROR: qmin must be less than qmax\n");
      exit(-1);
    }
  if (qminpu < 0)
    qminpu = 0.0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      printf("fname=%s argv[2]=%s\n",fname, argv[2]);
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
	    }
	  while (strcmp(line,"@@@"));
	  do 
    	    {
	      fscanf(f,"%[^\n]\n",line);
	    }
	  while (strcmp(line,"@@@"));
	  for (i=0; i < N; i++)
	    {
	      fscanf(f, "%[^\n]\n", line); 
	    }	
	  fscanf(f, "%lf %lf %lf\n", &L[0], &L[1], &L[2]);
      	  //invL = 1.0/L;
	  //printf("qui L=%.15G\n", L);
          if (qmaxpu < 0)
    	   qmaxpu = twopi/pow(L[0]*L[1]*L[2],1.0/3.0)*qmax;
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
  	}
      while (strcmp(line,"@@@"));
      //printf("fname=%s %d ellipsoids...\n", fname, N);

      if (first)
	{
	  for (a = 0; a < 3; a++)
	    {
	      r[a] = malloc(sizeof(double)*N);
	      u[a] = malloc(sizeof(double)*N);
	    }
	  //first = 0;
	}
      for (i=0; i < N; i++)
	{
	   fscanf(f, "%[^\n]\n", line); 
	   //printf("line=%s\n", line);
	   //if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
	     //{
	       //printf("boh\n");
	       sscanf(line, "%lf %lf %lf %lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], &u[0][i], &u[1][i], &u[2][i], dummy); 
	     //}
	   //printf("r=(%.15G,%.15G,%.15G)\n", r[0][i], r[1][i], r[2][i]);
	}
      if (NA == -1)
	NA = N;
      //NA = N;//force monodisperse
      for (a=0; a < 3; a++)
	scalFact[a] = twopi / L[a];
      if (first)
	{
	  if (physunit || qminpu != -1.0 || qmaxpu != -1.0)
	    {
	      for (qmod = 0; qmod < KMODMAX; qmod++)
		{
		  qavg[qmod] = 0.0;
		  for (mp = 0; mp < ntripl[qmod]; mp++) 
		    {
		      qval = sqrt(Sqr(((double)mesh[qmod][mp][0])*scalFact[0])+
		 		  Sqr(((double)mesh[qmod][mp][1])*scalFact[1])+
					 Sqr(((double)mesh[qmod][mp][2])*scalFact[2]));
		      qavg[qmod] += qval;
		      //sqrt(Sqr(((double)mesh[qmod][mp][0])*scalFact[0])+
		      //Sqr(((double)mesh[qmod][mp][1])*scalFact[1])+
		      //	Sqr(((double)mesh[qmod][mp][2])*scalFact[2]));
		    }
		  qavg[qmod] *= 1.0/((double)ntripl[qmod]);
		  //printf("qavg[%d]: %f\n", qmod, qavg[qmod]);
		  if (qmin==-1 && (qminpu < 0 || qavg[qmod] >= qminpu))
		    {
		      //printf("qui qmin=%d\n", qmod);
		      qmin = qmod;
		    }
		 
		  if (qmax==-1 && qmin != -1 && (qmaxpu > 0 && qavg[qmod] >= qmaxpu))
		    {
		      //printf("qui qmax=%d\n", qmod);
		      qmax = qmod-1;
		    }
		  //printf("qavg[%d]:%.15G - %.15G\n", qmod, qavg[qmod], scalFact*(1.25+0.5*qmod));
		}
	      if (qmax==-1)
		qmax = qmod-1;
#if 0
	      if (qminpu != -1.0 && qminpu == qmaxpu)
		{
		  //qmin = rint((qminpu-1.0) / (scalFact/2.0));
		  //qmax = qmin;
		}
	      else 
		{
		  //if (qminpu != -1.0 && qminpu >= 0.0)
		    //qmin = ceil( qminpu / (scalFact/2.0) - 2.0);
		  //if (qmaxpu != -1.0 && qmaxpu > 0.0)
		    //qmax = floor( qmaxpu / (scalFact/2.0) - 2.0);
		}
#endif
	      if (qmin < 0)
		qmin = 0;
	      if (qmin >= KMODMAX)
		qmax = KMODMAX-1; 
	      if (qmax >= KMODMAX)
		qmax = KMODMAX-1;
	      if (qmax < 0)
		qmax = 0;
	    }
	  printf("L: %f %f %f qmin: %d qmax: %d qminpu: %.15G qmaxpu: %.15G\n", L[0], L[1], L[2], 
		 qmin, qmax, qminpu, qmaxpu);
	  for (qmod=qmin; qmod <= qmax; qmod++)
	    {
	      Sq[qmod] = 0.0;      
	      if (NA < N)
		{
	    	  SqAA[qmod] = SqAB[qmod] = SqBB[qmod] = 0.0;
		}
	    }
	  //first = 0;
	}
      invNm = 1.0 / ((double)N);
      if (NA < N)
	{
	  invNmAA = 1.0/((double)NA);
	  invNmBB = 1.0/((double)N-NA);
	  invNmAB = 1.0/sqrt(((double)N-NA)*((double)NA));
	}
      if (first)
	{
	  if (NA < N)
	    printf("[MIXTURE N=%d NA=%d] ", N, NA);
	  else 
	    printf("[MONODISPERSE] ");
	  if (eventDriven)
	    printf("[ED] ");
	  else
	    printf("[MD]");
	  printf("twopi=%.15G N=%d invNm:%.15G qmin: %d qmax: %d\n", twopi, N, invNm,
		 qmin, qmax);
	  first = 0;
	}
      for(n = qmin; n <= qmax; n++)
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
		  if (bigpatchsq==1)
		    {
		      rx=r[0][i]+u[0][i]*patchdist;
		      ry=r[1][i]+u[1][i]*patchdist;
		      rz=r[2][i]+u[2][i]*patchdist;
		      rCMk = kbeg + 
			(scalFact[0]*rx * mesh[n][mp][0] + scalFact[1]*ry * mesh[n][mp][1] + 
			 scalFact[2]*rz * mesh[n][mp][2]);
		    }
		  else
		    {
		      rCMk = kbeg + 
			(scalFact[0]*r[0][i] * mesh[n][mp][0] + scalFact[1]*r[1][i] * mesh[n][mp][1] + 
			 scalFact[2]*r[2][i] * mesh[n][mp][2]);
		    }
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
  fclose(f2); 
  of = fopen("Sq.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      Sq[qmod] = (Sq[qmod]  * invNm) / ((double) ntripl[qmod]) / ((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	{
	  fprintf(of, "%.15G %.15G\n", qavg[qmod], Sq[qmod]); 
	}
      else
	{
	  fprintf(of, "%d %.15G\n", qmod, Sq[qmod]); 
	}
    }
  fclose(of);
  if (NA < N) 
    {
      of = fopen("SqAA.dat", "w+");
      for (qmod = qmin; qmod <= qmax; qmod++)
	{
	  SqAA[qmod] = (SqAA[qmod]  * invNmAA) / ((double) ntripl[qmod]) / ((double)nf);  
	  //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
	  if (physunit)
	    fprintf(of, "%.15G %.15G\n", qavg[qmod], SqAA[qmod]); 
	  else
	    fprintf(of, "%d %.15G\n", qmod, SqAA[qmod]); 
	}
      fclose(of);
      of = fopen("SqBB.dat", "w+");
      for (qmod = 0; qmod <= qmax; qmod++)
	{
	  SqBB[qmod] = (SqBB[qmod]  * invNmBB) / ((double) ntripl[qmod]) / ((double)nf);  
	  //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
	  if (physunit)
	    fprintf(of, "%.15G %.15G\n", qavg[qmod], SqBB[qmod]); 
	  else
	    fprintf(of, "%d %.15G\n", qmod, SqBB[qmod]); 
	}
      fclose(of);
      of = fopen("SqAB.dat", "w+");
      for (qmod = 0; qmod <= qmax; qmod++)
	{
	  SqAB[qmod] = (SqAB[qmod]  * invNmAB) / ((double) ntripl[qmod]) / ((double)nf);  
	  //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
	  if (physunit)
	    fprintf(of, "%.15G %.15G\n", qavg[qmod], SqAB[qmod]); 
	  else
	    fprintf(of, "%d %.15G\n", qmod, SqAB[qmod]); 
	}
      fclose(of);
    }
  return 0;
}
