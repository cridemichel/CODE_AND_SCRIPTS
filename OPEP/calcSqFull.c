#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
char line[1000000], parname[124], parval[2560000];
int N, NA=-1;
double x[3], *r[3];
char fname[1024], inputfile[1024];
int readCnf = 0, physunit=0;
const int numat=1038, numprot=64;
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
double Sq[KMODMAX], SqCM[KMODMAX], sumRho, reRho, imRho, rCM[3], rCMk, scalFact, invNm, invL, L, invNmAA, invNmBB, invNmAB;
double P[KMODMAX], reF[KMODMAX], imF[KMODMAX];
double SqAA[KMODMAX], SqBB[KMODMAX], SqAB[KMODMAX], sumRhoAA, sumRhoAB, sumRhoBB, reRhoA, reRhoB, imRhoA, imRhoB;
double rr[3], rpk, rp[3], reRhoCM, imRhoCM, sumRhoCM, reRhoFns, imRhoFns, sumRhoFns, sumreRhoFns, sumimRhoFns, rk;
double SqAAovFF, Fq;
void print_usage(void)
{
  printf("calcSqFull [ --qminpu/-qpum | --qmaxpu/-qpuM | --qmin/-qm <qmin> | --qmax/qM <qmax> |--help/-h | --cnf/-c | --phys-unit/-pu] <confs_file> [qmin] [qmax]\n");
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
  double ti, tref=0.0, kbeg=0.0, Vol, a1, a2, a3, a4, shift[3];
  int qmod, first = 1, NP1, NP2, kk, np;
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
	      L = Vol;
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
	   if (!(sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3))
	     {
	       //printf("boh\n");
	       sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
	     }
	   //printf("r=(%.15G,%.15G,%.15G)\n", r[0][i], r[1][i], r[2][i]);
	}
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
	  printf("N:%d L:%f scalFact: %.15G qmin: %d qmax: %d qminpu: %.15G qmaxpu: %.15G\n", N, L, scalFact, qmin, qmax, qminpu, qmaxpu);
	  for (qmod=qmin; qmod <= qmax; qmod++)
	    {
	      SqAA[qmod] = SqCM[qmod] = P[qmod] = reF[qmod] = imF[qmod] = 0.0;      
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
	  printf("twopi=%.15G N=%d L=%f invL=%.15G invNm:%.15G qmin: %d qmax: %d\n", twopi, N, L, invL, invNm,
		 qmin, qmax);
	  first = 0;
	}
      if (NA == N)
	{
	  for(n = qmin; n <= qmax; n++)
	    {
	      sumRho = 0.0;
	      sumRhoFns = 0.0;
	      sumreRhoFns = 0.0;
	      sumimRhoFns = 0.0;
	      sumRhoCM = 0.0;
	      for(mp = 0; mp < ntripl[n]; mp++)
		{
		  reRho = 0.0;
		  imRho = 0.0;
		  reRhoCM=0.0;
		  imRhoCM=0.0;
		  for(np=0; np < numprot; np++)
		    {
		      rCM[0]=rCM[1]=rCM[2]=0;
		      for (a=0; a < numat; a++)
			{
			  i = np*numat + a;	  
			  for (kk=0; kk < 3; kk++)
			    rr[kk] = r[kk][i];
			  if (a > 0)
			    { 
			      shift[0] = L*rint((r[0][i]-r[0][np*numat])/L);
			      shift[1] = L*rint((r[1][i]-r[1][np*numat])/L);
			      shift[2] = L*rint((r[2][i]-r[2][np*numat])/L);
			      rr[0] -= shift[0];
			      rr[1] -= shift[1];
			      rr[2] -= shift[2];
			    }
			  for (kk=0; kk < 3; kk++)  
			    rCM[kk] += rr[kk];
			}	
		      for (kk=0; kk < 3; kk++)  
	    		rCM[kk] /= ((double)numat);
		     // printf("2)CM[%d] %f %f %f\n", np, rCM[0], rCM[1], rCM[2]);
		      reRhoFns=0.0;
		      imRhoFns=0.0;
		      for (a=0; a < numat; a++)
			{
			  i = np*numat + a;	  
			  /* il passo della mesh e' 0.5*pi2/L */
			  if (mesh[n][mp][0]==0 && mesh[n][mp][1] == 0 && 
			      mesh[n][mp][2] == 0)
			    {
			      printf("ERRORE nella MESH!!!!!!!! n=%d mp=%d ntripl[n]:%d\n", n,
				     mp, ntripl[n]);
			      exit(-1);
			    }
		  	  for (kk=0; kk < 3; kk++)
			    rr[kk] = r[kk][i];
			  if (a > 0)
			    { 
			      shift[0] = L*rint((r[0][i]-r[0][np*numat])/L);
			      shift[1] = L*rint((r[1][i]-r[1][np*numat])/L);
			      shift[2] = L*rint((r[2][i]-r[2][np*numat])/L);
			      rr[0] -= shift[0];
			      rr[1] -= shift[1];
			      rr[2] -= shift[2];
			    }
			  for (kk=0; kk < 3; kk++)
			    rp[kk] = rr[kk] - rCM[kk];
			  rpk = kbeg + scalFact * 
			    (rp[0] * mesh[n][mp][0] + rp[1] * mesh[n][mp][1] + 
			     rp[2] * mesh[n][mp][2]);
			  rk = kbeg + scalFact * 
			    (r[0][i] * mesh[n][mp][0] + r[1][i] * mesh[n][mp][1] + 
			     r[2][i] * mesh[n][mp][2]);
			  reRhoFns += cos(rpk);
			  imRhoFns += sin(rpk);
			  reRho = reRho + cos(rk); 
			  imRho = imRho + sin(rk);
			  /* Imaginary part of exp(i*k*r) for the actual molecule*/
			}
		      rCMk = scalFact * 
			    (rCM[0] * mesh[n][mp][0] + rCM[1] * mesh[n][mp][1] + 
			     rCM[2] * mesh[n][mp][2]);
		      reRhoCM += cos(rCMk);
		      imRhoCM += sin(rCMk);
		      sumRhoFns += Sqr(reRhoFns) + Sqr(imRhoFns);
		      sumreRhoFns += reRhoFns;
		      sumimRhoFns += imRhoFns;
		    }
		  sumRho = sumRho + Sqr(reRho) + Sqr(imRho);
		  sumRhoCM = sumRhoCM + Sqr(reRhoCM) + Sqr(imRhoCM);
		}
	      SqAA[n] += sumRho;  
	      SqCM[n] += sumRhoCM;
	      P[n] += sumRhoFns;
	      reF[n] += sumreRhoFns;
	      imF[n] += sumimRhoFns;	
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
  of = fopen("SqAA.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      SqAA[qmod] = SqAA[qmod]/(((double)numprot)*numat) / ((double) ntripl[qmod]) / ((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], SqAA[qmod]); 
      else
	fprintf(of, "%d %.15G\n", qmod, SqAA[qmod]); 
    }
  fclose(of);
 
  of = fopen("SqCM.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      SqCM[qmod] = (SqCM[qmod]/((double)numprot)) / ((double) ntripl[qmod]) / ((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], SqCM[qmod]); 
      else
	fprintf(of, "%d %.15G\n", qmod, SqCM[qmod]); 
    }
  fclose(of);
  
  of = fopen("SqAAovFF.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      P[qmod] = P[qmod] / ((double) ntripl[qmod]) / ((double)nf) / ((double) numprot);  
      SqAAovFF = SqAA[qmod] * ((double)numat) / P[qmod];
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], SqAAovFF); 
      else
	fprintf(of, "%d %.15G\n", qmod, SqAAovFF); 
    }
  fclose(of);
  
  of = fopen("SqAAovFFbetter.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      reF[qmod] = reF[qmod] / ((double) ntripl[qmod]) / ((double)nf) /((double) numprot);
      imF[qmod] = imF[qmod] / ((double) ntripl[qmod]) / ((double)nf) /((double) numprot); 
      Fq = Sqr(reF[qmod])+Sqr(imF[qmod]);
      SqAAovFF = (SqAA[qmod]*numat - P[qmod])/Fq + 1.0;  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G %.15G\n", qavg[qmod], SqAAovFF, Fq); 
      else
	fprintf(of, "%d %.15G %.15G\n", qmod, SqAAovFF, Fq); 
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
