#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
#define Sqr(x) ((x)*(x))
char line[1000000], parname[124], parval[2560000];
int N, NA=-1;
double L[3], invL[3];
double x[3], *r[3], nemthr=0.2;
char fname[1024], inputfile[1024];
int readCnf = 0, physunit=0;
#define KMODMAX 599
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
int nonem=0;
double qavg[KMODMAX], nv[3], scalprodnvq, meshvec[3];
int qmin=0, qmax=KMODMAX-1, eventDriven = 0;
double qminpu=-1.0, qmaxpu=-1.0;
int ntripl[]=
#include "./ntripl.dat"
int mesh[][NKSHELL][3]= 
#include "./kmesh.dat"
double twopi, sinrCMk, cosrCMk;
double Sq[KMODMAX], SqPara[KMODMAX], SqPerp[KMODMAX], ccperp[KMODMAX], ccpara[KMODMAX],
       sumRho, reRho, sumRhoPara, sumRhoPerp, imRho, rCMk, scalFact, invNm, invNmAA, invNmBB, invNmAB;
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
      else if (!strcmp(argv[cc],"--nemvector")||!strcmp(argv[cc],"-nv"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  sscanf(argv[cc],"(%lf,%lf,%lf)", &(nv[0]), &(nv[1]), &(nv[2]));
	}
      else if (!strcmp(argv[cc],"--nemthr") || !strcmp(argv[cc],"-nt"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  nemthr = atof(argv[cc]);
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
double scalProd(double A[], double B[])
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
double calc_norm(double vec[3])
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

int main(int argc, char** argv)
{
  FILE *f, *f2, *of;
  int nf, i, a, b, n, mp, k1;
  double ti, tref=0.0, kbeg=0.0, Vol, a1, a2, a3, a4;
  int qmod, first = 1, NP1, NP2;
  nv[0] = 1.0;
  nv[1] = 0.0;
  nv[2] = 0.0;
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
      fscanf(f, "%d %lf %lf %lf\n", &N, &(L[0]), &(L[1]), &(L[2]));
      invL[0] = 1.0/L[0];
      invL[1] = 1.0/L[1];
      invL[2] = 1.0/L[2];

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
      scalFact = twopi;
      if (physunit || qminpu != -1.0 || qmaxpu != -1.0)
        {
          for (qmod = 0; qmod < KMODMAX; qmod++)
            {
              qavg[qmod] = 0;
              for (mp = 0; mp < ntripl[qmod]; mp++) 
                {
                  qavg[qmod] += sqrt(Sqr(invL[0]*((double)mesh[qmod][mp][0]))+
                                     Sqr(invL[1]*((double)mesh[qmod][mp][1]))+
                                     Sqr(invL[2]*((double)mesh[qmod][mp][2])));
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
	  printf("scalFact: %.15G qmin: %d qmax: %d qminpu: %.15G qmaxpu: %.15G\n", scalFact, qmin, qmax, qminpu, qmaxpu);
	  for (qmod=qmin; qmod <= qmax; qmod++)
	    {
	      Sq[qmod] = SqPerp[qmod] = SqPara[qmod] = 0.0;      

	      if (NA < N)
		{
	    	  SqAA[qmod] = SqAB[qmod] = SqBB[qmod] = 0.0;
		}
	    }
	  //first = 0;
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
	  printf("twopi=%.15G N=%d invL=%.15G %.15G %.15G invNm:%.15G qmin: %d qmax: %d\n", twopi, N, 
                 invL[0], invL[1], invL[2], invNm,
		 qmin, qmax);
	  first = 0;
	}
      if (NA == N)
	{
	  for(n = qmin; n <= qmax; n++)
	    {
	      sumRho = sumRhoPerp = sumRhoPara = 0.0;
	      ccperp[n] = ccpara[n] = 0.0;
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
			(invL[0]*r[0][i] * mesh[n][mp][0] + invL[1]*r[1][i] * mesh[n][mp][1] + 
			 invL[2]*r[2][i] * mesh[n][mp][2]);
		      //printf("i=%d r=%f %f %f\n", i, r[0][i], r[1][i], r[2][i]);
		      cosrCMk = cos(rCMk);
		      sinrCMk = sin(rCMk);
		      reRho = reRho + cosrCMk; 
		      imRho = imRho + sinrCMk;
		      /* Imaginary part of exp(i*k*r) for the actual molecule*/
		    }
		  sumRho = sumRho + Sqr(reRho) + Sqr(imRho);
		  for (k1=0; k1 < 3; k1++)
		    meshvec[k1] = ((double)mesh[n][mp][k1]);
	    	  scalprodnvq=fabs(scalProd(meshvec, nv)/calc_norm(nv)/calc_norm(meshvec));
    		  if (scalprodnvq < nemthr) 
		    {
		      sumRhoPerp += Sqr(reRho) + Sqr(imRho);
		      ccperp[n]++;
		    }
		  else if (fabs(1.0-scalprodnvq) < nemthr)
		    {
		      sumRhoPara += Sqr(reRho) + Sqr(imRho);
		      ccpara[n]++;
	    	    }
		}
	      Sq[n] += sumRho;  
	      SqPerp[n] += sumRhoPerp;
	      SqPara[n] += sumRhoPara; 
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
  of = fopen("Sq.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      Sq[qmod] = (Sq[qmod]  * invNm) / ((double) ntripl[qmod]) / ((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], Sq[qmod]); 
      else
	fprintf(of, "%d %.15G\n", qmod, Sq[qmod]); 
    }
  fclose(of);
  of = fopen("SqPerp.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      if (ccperp[qmod]==0)
        continue;
      SqPerp[qmod] = (SqPerp[qmod]  * invNm) / ccperp[qmod] /((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], SqPerp[qmod]); 
      else
	fprintf(of, "%d %.15G\n", qmod, SqPerp[qmod]); 
    }
  fclose(of);
  of = fopen("SqPara.dat", "w+");
  for (qmod = qmin; qmod  <= qmax; qmod++)
    {
      if (ccpara[qmod]==0)
        continue;
       SqPara[qmod] = (SqPara[qmod]  * invNm) / ccpara[qmod] /((double)nf);  
      //printf("nf=%d ntripl[%d]=%d\n", nf, qmod, ntripl[qmod]);
      if (physunit)
	fprintf(of, "%.15G %.15G\n", qavg[qmod], SqPara[qmod]); 
      else
	fprintf(of, "%d %.15G\n", qmod, SqPara[qmod]); 
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
