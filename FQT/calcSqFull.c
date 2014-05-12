#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
char line[1000000], parname[124], parval[2560000];
int N, NA=-1, Npts=100, Nptseff;
double x[3], *r[3], sax[3], *R[3][3];
char fname[1024], inputfile[1024];
int readCnf = 0, physunit=0;
#define KMODMAX 599
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double qavg[KMODMAX];
int qmin=0, qmax=KMODMAX-1, eventDriven = 0;
double qminpu=-1.0, qmaxpu=-1.0;
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi;
double Sq[KMODMAX], sumRho, reRho, imRho, rCMk, scalFact, invNm, invL, L, invNmAA, invNmBB, invNmAB;
double SqAA[KMODMAX], SqBB[KMODMAX], SqAB[KMODMAX], sumRhoAA, sumRhoAB, sumRhoBB, reRhoA, reRhoB, imRhoA, imRhoB;
void print_usage(void)
{
  printf("calcSq [ --qminpu/-qpum | --qmaxpu/-qpuM | --qmin/-qm <qmin> | --qmax/qM <qmax> |--help/-h | --cnf/-c | --phys-unit/-pu] <confs_file> [qmin] [qmax]\n");
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
      else if (!strcmp(argv[cc],"--Npts") || !strcmp(argv[cc],"-Np"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  Npts = atoi(argv[cc]);
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
double *rmesh[3];

int generate_mesh(int N, double sax[3])
{
  int cc, k1, k2, k3, fine=0;
  double dx, dy, dz, rx, ry, rz;
  double dens, pi;

  pi = 2.0*acos(0.0);
  for (k1=0; k1 < 3; k1++)
    rmesh[k1]=malloc(sizeof(double)*N); 
  dens = N / (sax[0]*sax[1]*sax[2]*4.0*pi/3.0);
  dx=dy=dz=pow((sax[0]*sax[1]*sax[2]*4.0*pi/3.0)/N,1.0/3.0);
  cc=0;
  printf("vol=%f dx=%f sax=%f %f %f N=%d\n",sax[0]*sax[1]*sax[2]*4.0*pi/3.0, dx, sax[0], sax[1], sax[2], N);
  rx = -sax[0];
  ry = -sax[1];
  rz = -sax[2];
  while (!fine)
    {
      if (Sqr(rx/sax[0]) + Sqr(ry/sax[1]) + Sqr(rz/sax[2]) <= 1.0)
	{
	  printf("assigned cc=%d r=%f %f %f\n", cc, rx, ry, rz);
	  rmesh[0][cc] = rx;
	  rmesh[1][cc] = ry;
	  rmesh[2][cc] = rz;
	  cc++;
	}
      if (cc >= N) 
	{
	  printf("OK all points set!\n");
	  fine=1;
	  break;
	}
      rx += dx;
      if (rx > sax[0])
	{
	  rx = -sax[0];
	  ry += dy;
	  if (ry > sax[1])
	    {
	      ry = -sax[1];
	      rz += dz;
	      if (rz > sax[2])
		fine=1;
	    }
	}
    }
  printf("Assigned # %d points out of %d\n", cc, N);
  return cc; 
}
double **rmeshLab[3];
void body2lab(int N, int Npts)
{
  int k1, k2, a;
  int i;

  for (i=0; i < N; i++)
    {
      for (a = 0; a < Npts; a++)
	{
	  for (k1=0; k1 < 3; k1++)
	    {
	      rmeshLab[k1][i][a] = 0;
	      for (k2=0; k2 < 3; k2++)
	     	{
		  rmeshLab[k1][i][a] += R[k2][k1][i]*rmesh[k2][a];
		} 
	      rmeshLab[k1][i][a] += r[k1][i];
	    }
	}
    }
}
int main(int argc, char** argv)
{
  FILE *f, *f2, *of;
  int nf, i, a, b, n, mp, cc=0, type, k1, k2, ccq;
  double ti, tref=0.0, kbeg=0.0, Vol, a1, a2, a3, a4;
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

	  while (strcmp(line,"@@@"));
	  cc=0;
	  do 
	    {
	      if (cc==1)
		fscanf(f, "%lf %lf %lf\n", &sax[0], &sax[1], &sax[2]);
	      else
		fscanf(f,"%[^\n]\n",line);
	      cc++;
	    }
	  while (strcmp(line,"@@@"));
	  Nptseff=generate_mesh(Npts, sax);
	  printf("===>N=%d Npts=%d\n", N, Npts);
	  for (k1=0; k1 < 3; k1++)
	    { 
	      rmeshLab[k1]=malloc(sizeof(double*)*N); 
	      for (k2=0; k2 < 3; k2++)
		R[k1][k2] = malloc(sizeof(double)*N);
	      for (k2 = 0; k2 < N; k2++)
		rmeshLab[k1][k2]=malloc(sizeof(double)*Npts); 
	    }
	  for (i=0; i < N; i++)
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
  	}
      while (strcmp(line,"@@@"));
      //printf("fname=%s %d ellipsoids...\n", fname, N);

      if (first)
	{
	  for (a = 0; a < 3; a++)
	    {
	      r[a] = malloc(sizeof(double)*N);
	    }
	  //first = 0;
	}
      for (i=0; i < N; i++)
	{
	   fscanf(f, "%[^\n]\n", line); 
	   //printf("line=%s\n", line);
       	   sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &r[0][i], &r[1][i], &r[2][i], 
		  &R[0][0][i], &R[0][1][i], &R[0][2][i], 
		  &R[1][0][i], &R[1][1][i], &R[1][2][i], 
		  &R[2][0][i], &R[2][1][i], &R[2][2][i], &type); 
	   //printf("r=(%.15G,%.15G,%.15G)\n", r[0][i], r[1][i], r[2][i]);
	}
      body2lab(N, Npts);
      if (NA == -1)
	NA = N;
      //NA = N;//force monodisperse
      scalFact = twopi * invL;
      if (first)
	{
	  if (physunit || qminpu != -1.0 || qmaxpu != -1.0)
	    {
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
		  qmin = rint((qminpu-1.0) / (scalFact/2.0));
		  qmax = qmin;
		}
	      else 
		{
		  if (qminpu != -1.0 && qminpu >= 0.0)
		    qmin = ceil( qminpu / (scalFact/2.0) - 2.0);
		  if (qmaxpu != -1.0 && qmaxpu > 0.0)
		    qmax = floor( qmaxpu / (scalFact/2.0) - 2.0);
		}
	      if (qmin < 0)
		qmin = 0;
	      if (qmin >= KMODMAX)
		qmax = KMODMAX-1; 
	      if (qmax >= KMODMAX)
		qmax = KMODMAX-1;
	      if (qmax < 0)
		qmax = 0;
	    }
	  printf("sax= %f %f %f scalFact: %.15G qmin: %d qmax: %d qminpu: %.15G qmaxpu: %.15G\n", sax[0], sax[1], sax[2], scalFact, qmin, qmax, qminpu, qmaxpu);
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
	  printf("twopi=%.15G N=%d invL=%.15G invNm:%.15G qmin: %d qmax: %d\n", twopi, N, invL, invNm,
		 qmin, qmax);
	  first = 0;
	}
      if (NA == N)
	{
	  for(n = qmin; n <= qmax; n++)
	    {
	      sumRho = 0.0;
	      for(mp = 0; mp < ntripl[n]; mp++)
		{
		  reRho = 0.0;
		  imRho = 0.0;
		  for(i=0; i < N; i++)
		    {
		      for (a = 0; a < Npts; a++)
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
			    (rmeshLab[0][i][a] * mesh[n][mp][0] + rmeshLab[1][i][a] * mesh[n][mp][1] + 
			     rmeshLab[2][i][a] * mesh[n][mp][2]);
			  reRho = reRho + cos(rCMk); 
			  imRho = imRho + sin(rCMk);
			  /* Imaginary part of exp(i*k*r) for the actual molecule*/
			}
		      sumRho = sumRho + Sqr(reRho) + Sqr(imRho);
		    }
		}
	      Sq[n] += sumRho;  
	      //printf("sumRho=%.15G Sq[%d]=%.15G\n", sumRho, n, Sq[n]);
	    }
	  fclose(f);
	}
      else
	{
	  for(n = qmin; n <= qmax; n++)
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
		      if (i < NA)
			{
			  reRhoA = reRhoA + cos(rCMk); 
			  imRhoA = imRhoA + sin(rCMk);
		       	}
		      else 
			{
			  reRhoB = reRhoB + cos(rCMk); 
			  imRhoB = imRhoB + sin(rCMk);
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
  of = fopen("SqFull.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      Sq[qmod] = Sq[qmod] * (1.0 / (((double)Nptseff)*((double)N))) / ((double) ntripl[qmod]) / ((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], Sq[qmod]); 
      else
	fprintf(of, "%d %.15G\n", qmod, Sq[qmod]); 
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
