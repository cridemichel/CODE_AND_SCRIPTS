#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define Sqr(x) ((x)*(x))
//#define SQ_BINNING
#ifdef ALL_N_INT
#undef ALL_N_INT
#endif
double lenHC=4.0, surfRad=10.0, twopiL[3], dk[3], nemv[3], k13max[3], k13sqmax[3], dksq[3];
char line[1000000], parname[124], parval[1000000];
char dummy[2048];
int mcsim=0, cubic_box=1, nonem=0, calcnv=0, npoints[3];
int N, particles_type=1, k1, k2, nmax[3]={-1,-1,-1}, n13x, n13z, nv[3];
double *x[3], L, ti, *w[3], storerate, k13x, k13z, k13sqx, k13sqz, kmaxmod, kv[3]; 
double Lx, Ly, Lz, V;
double vecx[3], vecy[3], vecz[3], *R[3][3];
double *DR[3], deltaAA=-1.0, deltaBB=-1.0, deltaAB=-1.0, sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0, 
       Dr, theta, sigmaSticky=-1.0, sa[2], sb[2], sc[2], maxax0, maxax1, maxax, maxsax, maxsaxAA, maxsaxAB, maxsaxBB;
char fname[1024], inputfile[1024];
int eventDriven=0;
int points1=50, points3=-1;
int foundDRs=0, foundrot=0, minbins=0;

void readconf(char *fname, double *ti, double *refTime, int *NP, int *NPA, double *r[3], double *w[3], double *DR[3])
{
  FILE *f;
  int nat=0, i, cpos, dummyint;
  double dt=-1;
  int curstp=-1, NP1, NP2;
  *ti = -1.0;
  f = fopen(fname, "r");
  while (!feof(f)) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n] ",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  //printf("qui nat=%d\n", nat);
	  continue;
	}
	//printf("line=%s\n", line);
      if (nat < 2)
	{
	  fseek(f, cpos, SEEK_SET);
	  fscanf(f, "%[^:]:", parname);
	  //printf("[%s] parname=%s\n", fname, parname);
	  if (!strcmp(parname,"parnum"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
		{
		  *NP = atoi(parval);
		}
	      else
		{
		  *NP = NP1+NP2;
		  *NPA = NP1;
		}
	    }
	  else if (!strcmp(parname,"parnumA"))
	    *NPA = atoi(parval);
	  else if (!strcmp(parname,"DR"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf %lf %lf ", &DR[0][i], &DR[1][i], &DR[2][i]);
		}
	      foundDRs = 1;
	    }
	  else if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	      foundrot = 1;
	    }
	  else if (!strcmp(parname,"sumoy"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf ", &w[1][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoz"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf ", &w[2][i]); 
		}
	    }
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
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
	  else if (!strcmp(parname, "refTime"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *refTime = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else if (nat==3)
	{
	  if (*NPA==-1)
	    *NPA = *NP;
	  fseek(f, cpos, SEEK_SET);
	  for (i = 0; i < *NP; i++) 
	    {
	      fscanf(f, "%[^\n]\n", line); 
#if 0	      
	      if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
		{
		  sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
		}
#endif
	      sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
	    	      &(r[0][i]), &(r[1][i]), &(r[2][i]), 
	    	      &(R[0][0][i]), &(R[0][1][i]), &(R[0][2][i]), &(R[1][0][i]), &(R[1][1][i]), 
	    	      &(R[1][2][i]), &(R[2][0][i]), &(R[2][1][i]), &(R[2][2][i]), &dummyint);

	      //printf("r[%d]=%f %f %f\n", i, r[0][i], r[1][i], r[2][i]);
	    }
	  if (mcsim==1)
	    {
	      if (fscanf(f, "%lf %lf %lf\n", &Lx, &Ly, &Lz)==1)
		{
		  L=Lx;
		}
	      break;
	    }
	  else
	    {
	      /* if MD read vels */	
	      for (i=0; i < *NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      // fscanf(f, "%lf\n", &L);
	      if (fscanf(f, "%lf %lf %lf\n", &Lx, &Ly, &Lz)==1)
		{
		  L=Lx;
		}
	      break;
	    }

	  break; 
	}
    }
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
  if (*ti == -1)
    *ti = ((double)curstp)*dt;
  fclose(f);
}

void print_usage(void)
{
#ifdef SQ_BINNING
  printf("calc_elconst_W [--lenhc/-lhc <val>] [--minbins/-mb <val>] [-k1M/--k1sqmax <val>]  [-k3M/--k3sqmax <val>] [--nxmax/nxM <val>] [--nymax/nyM <val>] [--nzmax/nzM <val>] [--mcsim/-mc] <confs_file> [points1] [points3]\n");
#else
  printf("calc_elconst_W [--lenhc/-lhc <val>] [--minbins/-mb <val>] [-k1M/--k1max <val>]  [-k3M/--k3max <val>] [--nxmax/nxM <val>] [--nymax/nyM <val>] [--nzmax/nzM <val>] [--mcsim/-mc] <confs_file> [points1] [points3]\n");
#endif
  exit(0);
}
double threshold=0.05;
int gnuplot=0, isoav=0;
int parse_param(int argc, char** argv)
{
  int extraparam=0;  
  int cc=1;
  points1=100;
  points3=-1;
#ifdef SQ_BINNING
  k13sqmax[0]=k13sqmax[1]=k13sqmax[2]=-1.0;
#else
  k13max[0]=k13max[1]=k13max[2]=-1.0;
#endif
  if (argc==1)
    {
      print_usage();
      exit(1);
    }
  while (cc < argc)
    {

      printf("argc=%d argv[argc]=%s cc=%d\n", argc, argv[cc], cc);
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
       else if (!strcmp(argv[cc],"--mcsim")||!strcmp(argv[cc],"-mc"))
	{
	  mcsim=1;
	}
      else if (!strcmp(argv[cc],"--gnuplot")||!strcmp(argv[cc],"-gp"))
	{
	  gnuplot=1;
	}
      else if (!strcmp(argv[cc],"--minbins")||!strcmp(argv[cc],"-mb"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  minbins = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--lenhc")||!strcmp(argv[cc],"-lhc"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  lenHC = atof(argv[cc]);
	}
#ifdef SQ_BINNING
      else if (!strcmp(argv[cc],"--k1sqmax")||!strcmp(argv[cc],"-k1M"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  k13sqmax[0] = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--k3sqmax")||!strcmp(argv[cc],"-k3M"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  k13sqmax[2] = atof(argv[cc]);
	}

#else
      else if (!strcmp(argv[cc],"--k1max")||!strcmp(argv[cc],"-k1M"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  k13max[0] = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--k3max")||!strcmp(argv[cc],"-k3M"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  k13max[2] = atof(argv[cc]);
	}
#endif
      else if (!strcmp(argv[cc],"--nxmax")||!strcmp(argv[cc],"-nxM"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  nmax[0] = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--nymax")||!strcmp(argv[cc],"-nyM"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  nmax[1] = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--nzmax")||!strcmp(argv[cc],"-nzM"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  nmax[2] = atoi(argv[cc]);
	}
      else if (cc==argc|| extraparam==3)
	print_usage();
      else if (extraparam == 0)
	{
	  extraparam++;
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam==1)
	{
	  extraparam++;
	  points1=atoi(argv[cc]);
	}
      else if (extraparam==2)
	{
	  extraparam++;
	  points3=atoi(argv[cc]);
	}
      else 
	print_usage();
      cc++;
    }
  return 0;
}
double pi;
int *inCell[3]={NULL,NULL,NULL}, *cellList=NULL, cellsx, cellsy, cellsz;
int NP, NPA;
#if 0
void build_linked_list(void)
{
  double L2;
  int j, n;
  L2 = 0.5 * L;

  for (j = 0; j < cellsx*cellsy*cellsz + NP; j++)
    cellList[j] = -1;

  for (n = 0; n < NP; n++)
    {
      inCell[0][n] =  (x[0][n] + L2) * cellsx / L;
      inCell[1][n] =  (x[1][n] + L2) * cellsy / L;
      inCell[2][n] =  (x[2][n] + L2) * cellsz / L;
      j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + NP;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}
#endif
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
void lab2nem(double v[3], double vp[3])
{
  int k1, k2;
  double R[3][3];

  for (k1=0; k1 < 3; k1++)
    {
      R[0][k1] = vecx[k1];
      R[1][k1] = vecy[k1];
      R[2][k1] = vecz[k1];
      //if (isnan(cos(R[0][k1]))||isnan(R[1][k1])||isnan(R[2][k1]))
	//printf("mmmah\n");
    } 
  for (k1=0; k1 < 3; k1++)
    {
      vp[k1] = 0.0;
      for (k2=0; k2 < 3; k2++)
	vp[k1] += R[k1][k2]*v[k2];
    }
  //if(isinf(vp[0])|| isinf(vp[1]) || isinf(vp[2]))
    //printf("vp=%f %f %f\n", vp[0], vp[1], vp[2]);
  //printf("v=%f %f %f\n", v[0], v[1], v[2]);
}
double min3(double a, double b, double c)
{
  double m;
  m = a;
  if (b < m)
    m = b;
  if (c < m)
    m = c;
  return m;
}
double max3(double a, double b, double c)
{
  double m;
  m = a;
  if (b > m)
    m = b;
  if (c > m)
    m = c;
  return m;
}
double eigvec[3][3], *eigvec_n[3][3], eigvec_t[3][3];
int eigenvectors=1;
double S=0, Q[3][3], reQ[3][3], imQ[3][3], reQ13[3][3], imQ13[3][3], **accQsq13, **accQsq23;
double **nbins;
void diagonalize(double M[3][3], double ev[3])
{
  double a[9], work[45];
  char jobz, uplo;
  int info, i, j, lda, lwork;
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) a[j+3*i]=M[j][i];		
      //for(j=0; j<3; j++) a[j][i]=M[j][i];		
    }	
  lda = 3;
  if (eigenvectors)
    jobz='V';
  else
    jobz='N';
  uplo='U';
  lwork = 45;
  dsyev_(&jobz, &uplo, &lda, a, &lda, ev, work, &lwork,  &info);  
  if (!eigenvectors)
    return;
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) eigvec[i][j]=a[j+3*i];		
    }	
}
void build_ref_system(double kv[3], double *k13x, double *k13z)
{
  int k;
  double sp, norm;
  /* build the nematic reference system */
  if (nemv[0] > 0.0)
    {
      vecz[0] = nemv[0]; /* z axis is the nematic axis */
      vecz[1] = nemv[1];
      vecz[2] = nemv[2];
    }
  else
    {
      vecz[0] = -nemv[0]; /* z axis is the nematic axis */
      vecz[1] = -nemv[1];
      vecz[2] = -nemv[2];
    }
  //printf("vecz=%f %f %f\n", vecz[0],vecz[1],vecz[2]);
  sp = 0.0;
  for (k=0; k < 3 ; k++)
    sp+=kv[k]*vecz[k];
  *k13z = sp;
  for (k=0; k < 3 ; k++)
    vecx[k] = kv[k];

#ifdef FIXED_Z
  if (kv[1]==0.0 && kv[2]==0.0)
    {
      vecx[0]=0;
      vecx[1]=0;
      vecx[2]=1;
    }
  else
    {
      for (k=0; k < 3 ; k++)
	vecx[k] -= sp*vecz[k];
      norm = calc_norm(vecx);
      *k13x = norm;
      //printf("norm kv = %f norm k13=%f\n", calc_norm(kv), sqrt(Sqr(*k13x)+Sqr(*k13z)));
      for (k=0; k < 3 ; k++)
	vecx[k] /= norm;
      //printf("k13=%f %f vecz=%f %f %f kv=%f %f %f\n", *k13x, *k13z, vecz[0], vecz[1], vecz[2], kv[0], kv[1], kv[2]);
  //printf("delr=%.15G nematic=%f %f %f\n", delr, nv[0], nv[1], nv[2]);
    }
#else
  for (k=0; k < 3 ; k++)
    vecx[k] -= sp*vecz[k];
  norm = calc_norm(vecx);
  *k13x = norm;
  //printf("k.assnem=%f k13z=%f\n", kv[0]*vecz[0]+kv[1]*vecz[1]+kv[2]*vecz[2], *k13z);
  //printf("norm kv = %f norm k13=%f\n", calc_norm(kv), sqrt(Sqr(*k13x)+Sqr(*k13z)));
  for (k=0; k < 3 ; k++)
    vecx[k] /= norm;
  //printf("k.assex=%f  k13x=%f\n",          kv[0]*vecx[0]+kv[1]*vecx[1]+kv[2]*vecx[2], *k13x); 
#endif
#if 0
  vecx[0] = 1.0;
  vecx[1] = 1.0;
  vecx[2] = 1.0;
  if (vecx[0] == vecz[0] && vecx[1]==vecz[1] && vecx[2]==vecz[2])
    {
      vecx[0] = 1.0;
      vecx[1] = -1.0;
      vecx[2] = 1.0;
    }
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=vecx[k]*vecz[k];
  for (k=0; k < 3 ; k++)
    vecx[k] -= sp*vecz[k];
  norm = calc_norm(vecx);
  for (k=0; k < 3 ; k++)
    vecx[k] = vecx[k]/norm;
#endif

  vectProdVec(vecz, vecx, vecy);
#if 0
  if (isnan(vecx[0])|| isnan(vecy[0]) || isnan(vecz[0])||
      isnan(vecx[1])|| isnan(vecy[1]) || isnan(vecz[1])
      isnan(vecx[2])|| isnan(vecy[2]) || isnan(vecz[2]) )
      printf("boh\n");
#endif
}
void calc_nem_vec(void)
{
  int a, b, numev, i;
  double St, ev[3];
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Q[a][b] = 0.0;      
      }
  for (i=0; i < NP; i++)
    {
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  {
	    Q[a][b] += 1.5 * R[0][a][i]*R[0][b][i];
	    if (a==b)
	      Q[a][a] -= 0.5;
	  }
    }
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Q[a][b] /= ((double)NP);
      } 
  diagonalize(Q, ev);
  if (fabs(ev[0]) > fabs(ev[1]))
   { 
     St = ev[0];
     numev=0;
   }
  else
    {
      St = ev[1];
      numev=1;
    }  
  if (fabs(ev[2]) > St)
    {
      St = ev[2];
      numev=2;
    }
  S+=St;
  for (a=0; a < 3; a++)
    nemv[a] = eigvec[numev][a];
  //printf("St=%f nemv=%f %f %f\n", St, nemv[0], nemv[1], nemv[2]);
}

void calc_fourier_Q(double kv[3])
{
  int a, b, numev, i, k;
  double u[3], up[3];
  double St, ev[3], kr, coskr, sinkr;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	reQ[a][b] = imQ[a][b] = 0.0;      
      }
  for (i=0; i < NP; i++)
    {
      for (k=0; k < 3; k++)
	u[k] = R[0][k][i];
#ifdef TRANSFORM_Q
      for (k=0; k < 3; k++)
	up[k] = u[k];
#else
      lab2nem(u, up);
#endif
      kr = kv[0]*x[0][i]+kv[1]*x[1][i]+kv[2]*x[2][i];
#ifdef TRANSFORM_Q
      for (a = 0; a < 3; a++)
	for (b = a; b < 3; b++)
	  {
	    Q[a][b] = 1.5 * up[a]*up[b];
	    if (a==b)
	      Q[a][b] += -0.5 * up[a]*up[b];
	  }
      Q[1][0] = Q[0][1];
      Q[2][0] = Q[0][2];
      Q[2][1] = Q[1][2];
#else
      Q[0][2] = 1.5 * up[0]*up[2];
      Q[1][2] = 1.5 * up[1]*up[2];
#endif
      coskr = cos(kr);
      sinkr = sin(kr);
#ifdef TRANSFORM_Q
      for (a = 0; a < 3; a++)
	for (b = a; b < 3; b++)
	  {
	    reQ[a][b] += Q[a][b]*coskr;
	    imQ[a][b] += Q[a][b]*sinkr;
	  }
#else
      reQ[0][2] += Q[0][2]*coskr;
      imQ[0][2] += Q[0][2]*sinkr;
      reQ[1][2] += Q[1][2]*coskr;
      imQ[1][2] += Q[1][2]*sinkr;
#endif
      //if (isinf(cos(kr))||isinf(sin(kr))||isinf(Q[0][2])||isinf(Q[1][2]))
	//printf("up=%f %f %f kr=%f cos(kr)=%f sin(kr)=%f NP=%d\n", up[0], up[1], up[2], kr, cos(kr), sin(kr), NP);
    }
#ifdef TRANSFORM_Q
  for (a = 0; a < 3; a++)
    for (b = a; b < 3; b++)
      {
	reQ[a][b] *= V/((double)NP);
	imQ[a][b] *= V/((double)NP);
      }	
  reQ[1][0] = reQ[0][1];
  reQ[2][0] = reQ[0][2];
  reQ[2][1] = reQ[1][2];
  imQ[1][0] = imQ[0][1];
  imQ[2][0] = imQ[0][2];
  imQ[2][1] = imQ[1][2];
#else
  reQ[0][2] *= V/((double)NP);
  imQ[0][2] *= V/((double)NP);
  reQ[1][2] *= V/((double)NP);
  imQ[1][2] *= V/((double)NP); 
#endif
  //if (isinf(reQ[0][2])||isinf(imQ[1][2]))
    //printf("V=%f reQ=%f imQ=%f\n", V, reQ[0][2], imQ[1][2]);
}

int main(int argc, char **argv)
{
  FILE *f, *f2, *f1;
  int binx, biny, binz;
  int k, nf=0, i, a, b, nat, NN, j, ii, bin, binParaAv, binPerpAv;
  double rx, ry, rz, norm, g0para, g0perp, minpara, maxpara, minperp, maxperp;
  double r, delr, tref=0.0, Dx[3], DxNem[3], *g0, **g0Perp, **g0Parall, *g0ParaAv, *g0PerpAv,
	 g0m, distSq, rlower, rupper, cost, nIdeal, costParaAv, costPerpAv, kmax;
  double sp, time, refTime, RCUT;
  int iZ, jZ, iX, jX, iY, jY, NP1, NP2;
  double shift[3];
  double ene=0.0;
#ifdef TRANSFORM_Q
  int k1, k2, k3;
  double Rl[3][3];
#endif
  parse_param(argc, argv);
  f2 = fopen(inputfile,"r");
  fscanf(f2, "%[^\n]\n", fname);
  pi = acos(0)*2.0;
  nat = 0;
  f = fopen(fname,"r");
  while (!feof(f) && nat < 3) 
    {
      fscanf(f, "%[^\n]\n)", line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  if (nat==3)
	    {
	      if (mcsim==1)
	       {
		 for (i=0; i < NP; i++)
		   {
		     fscanf(f, "%[^\n]\n", line);
		   }
		 if (fscanf(f, "%lf %lf %lf\n", &Lx, &Ly, &Lz)==1)
		   {
		     cubic_box=1;
		     L=Lx;
		   }
		 else
		   {
		     cubic_box=0;
		   }
	         break;
               }
              else
               {
	         for (i=0; i < 2*NP; i++)
		   fscanf(f, "%[^\n]\n", line);
	         fscanf(f, "%lf\n", &L);
	         break;
                }
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
      else if (nat == 0 && !strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (!strcmp(parname,"storerate"))
	{
	  storerate = atof(parval);
	  eventDriven = 1;
	}
      else if (nat==1 && !strcmp(parname,"a"))
	sscanf(parval, "%lf %lf\n", &sa[0], &sa[1]);	
      else if (nat==1 && !strcmp(parname,"b"))
	sscanf(parval, "%lf %lf\n", &sb[0], &sb[1]);	
      else if (nat==1 && !strcmp(parname,"c"))
	sscanf(parval, "%lf %lf\n", &sc[0], &sc[1]);	
      else if (nat==1 && !strcmp(parname,"sigma"))
	sscanf(parval, "%lf %lf %lf %lf\n", &sigmaAA, &sigmaAB, &sigmaAB, &sigmaBB);	
      else if (nat==1 && !strcmp(parname,"delta"))
	sscanf(parval, "%lf %lf %lf %lf\n", &deltaAA, &deltaAB, &deltaAB, &deltaBB);	
      else if (nat==1 && !strcmp(parname,"sigmaSticky"))
	sigmaSticky = atof(parval);
      else if (nat==1 && !strcmp(parname,"theta"))
	theta = atof(parval);
      else if (nat==1 && !strcmp(parname,"Dr"))
	Dr = atof(parval);
    }
  fclose(f);
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
  for (a = 0; a < 3; a++)
    {
      x[a] = malloc(sizeof(double)*NP);
      w[a] = malloc(sizeof(double)*NP);
      /* allocate orientations */
      for (b=0; b < 3; b++)      
	R[a][b] = malloc(sizeof(double)*NP);
    }

  if (!eventDriven)
    L = cbrt(L);
  if (deltaAA!=-1)
    particles_type = 2; /* 2 means square well system*/
  else if (sigmaAA != -1.0)
    particles_type = 0;

  if (particles_type == 1)
    {
      maxax0 = sa[0];
      if (sb[0] > maxax0)
	maxax0 = sb[0];
      if (sc[0] > maxax0)
	maxax0 = sc[0];
      maxax1 = sa[1];
      if (sb[1] > maxax1)
	maxax1 = sb[1];
      if (sc[0] > maxax1)
	maxax1 = sc[1];
      maxsax = fabs(maxax1)+fabs(maxax0)+2.0*sigmaSticky;
      //printf("maxsax=%.15G\n", maxsax);
    }
  else if (particles_type == 0)
    {
      maxsaxAA = fabs(sigmaAA)+2.0*sigmaSticky;
      maxsaxAB = fabs(sigmaAB)+2.0*sigmaSticky;
      maxsaxBB = fabs(sigmaBB)+2.0*sigmaSticky;
    }
  else if (particles_type == 2)
    {
      maxsaxAA = sigmaAA + deltaAA;
      maxsaxAB = sigmaAB + deltaAB;
      maxsaxBB = sigmaBB + deltaBB;
      maxsax = maxsaxAA;
      if (maxsaxAB > maxsax)
	maxsax = maxsaxAB;
      if (maxsaxBB > maxsax)
	maxsax = maxsaxBB;
    }
  /* WARNING: se i diametri sono diversi va cambiato qua!! */ 
  if (particles_type == 1)
    RCUT = maxsax;
  else if (particles_type == 0)
    RCUT = maxsaxAA*1.01;
  else if (particles_type == 2)
    RCUT = maxsax*1.01;
  if (particles_type==1)
    {
      printf("SYSTEM: ELLISPOIDS - DGEBA\n");
    }
  else if (particles_type==0)
    {
      printf("SYSTEM: SPHERES 3-2\n");
    }
  else if (particles_type == 2)
    {
      printf("SYSTEM: SQUARE WELL\n");
    }
  for (ii = 0; ii < 3; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*NP);
    }

#if 0
  cellsx = L / RCUT;
  cellsy = L / RCUT;
  cellsz = L / RCUT;
  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+NP));
  inCell[0] = malloc(sizeof(int)*NP);
  inCell[1] = malloc(sizeof(int)*NP);
  inCell[2] = malloc(sizeof(int)*NP);
#endif
  if (cubic_box)
    {
      printf("L=%.15G\n", L);
      V = L*L*L;
      twopiL[0] = twopiL[1] = twopiL[2] = 2.0*M_PI/L;
    }
  else
    {
      printf("Lx=%.15G Ly=%.15G Lz=%.15G\n", Lx, Ly, Lz);
      V = Lx*Ly*Lz;
      twopiL[0] = 2.0*M_PI/Lx;
      twopiL[1] = 2.0*M_PI/Ly;
      twopiL[2] = 2.0*M_PI/Lz;
    }
  if (nmax[0]==-1)
    nmax[0] = 40;
  if (nmax[1]==-1)
    nmax[1] = 20;
  if (nmax[2]==-1)
    nmax[2] = 20;
  if (points3==-1)
    points3=points1;
  npoints[0] = points1;
  npoints[2] = points3;
  //kmaxmod = (sqrt(Sqr(twopiL[0]*nmax[0])+Sqr(twopiL[1]*nmax[1])+Sqr(twopiL[2]*nmax[2])));
  //kmaxmod = max3(twopiL[0]*nmax[0],twopiL[1]*nmax[1],twopiL[2]*nmax[2]);
#ifdef SQ_BINNING
  if (k13sqmax[2] == -1.0)
    k13sqmax[2] = Sqr(twopiL[0]*nmax[0]);
  if (k13sqmax[0] == -1.0)
    k13sqmax[0] = Sqr(twopiL[1]*nmax[1]);
  dksq[0] = k13sqmax[0]/((double)npoints[0]);
  dksq[2] = k13sqmax[2]/((double)npoints[2]);
  printf("k13sqmax=%f %f npoints13= %d %d dksq13=%f %f\n", k13sqmax[0], k13sqmax[2], npoints[0], npoints[2], dksq[0], dksq[2]);
#else
  if (k13max[2]==-1.0)
    k13max[2] = twopiL[0]*nmax[0];
  if (k13max[0]==-1.0)
    k13max[0] = twopiL[1]*nmax[1];
  dk[0] = k13max[0]/((double)npoints[0]);
  dk[2] = k13max[2]/((double)npoints[2]);
  printf("kmax13=%f %f npoints13= %d %d dk13=%f %f\n", k13max[0], k13max[2], npoints[0], npoints[2], dk[0], dk[2]);
#endif
  accQsq13 = malloc(sizeof(double*)*(npoints[0]+1));
  accQsq23 = malloc(sizeof(double*)*(npoints[0]+1));
  nbins =    malloc(sizeof(double*)*(npoints[0]+1));
  for (k=0; k < npoints[0]; k++) 
    {
      accQsq13[k] = malloc(sizeof(double)*(npoints[2]+1));
      accQsq23[k] = malloc(sizeof(double)*(npoints[2]+1)); 
      nbins[k]    = malloc(sizeof(double)*(npoints[2]+1));
    }
  for (nv[0] = 0; nv[0] < npoints[0]; (nv[0])++)
    for (nv[2] = 0; nv[2] < npoints[2]; (nv[2])++)
      {
	accQsq13[nv[0]][nv[2]] = 0.0;
	accQsq23[nv[0]][nv[2]] = 0.0;
	nbins[nv[0]][nv[2]] = 0.0;
      }

  nf = 0;
  rewind(f2);
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      nf++;
      printf("FILE #%d\n", nf);
      NPA = -1;
      readconf(fname, &time, &refTime, &NP, &NPA, x, w, DR);
      calc_nem_vec();
#ifdef FIXED_Z
      nemv[0] = 1.0;
      nemv[1] = 0.0;
      nemv[2] = 0.0;
#endif
#ifndef ALL_N_INT
      for (nv[0] = 0; nv[0] <= nmax[0]; (nv[0])++)
	for (nv[1] = 0; nv[1] <= nmax[1]; (nv[1])++)
	  for (nv[2] = 0; nv[2] <= nmax[2]; (nv[2])++)
#else
      for (nv[0] = -nmax[0]; nv[0] <= nmax[0]; (nv[0])++)
	for (nv[1] = -nmax[1]; nv[1] <= nmax[1]; (nv[1])++)
	  for (nv[2] = 0; nv[2] <= nmax[2]; (nv[2])++)
#endif
	    {
	      if (nv[0]==0 && nv[1]==0 && nv[2]==0)
		continue;
	      for (k=0; k < 3; k++)
		kv[k] = twopiL[k]*nv[k]; 
	      /* da eliminare !!! */	
	      //kv[1] = -kv[1];
	      /* ********* */
	      build_ref_system(kv, &k13x, &k13z);
//printf("qui\n");

#ifdef SQ_BINNING
#ifdef USE_RINT
	      n13x = rint(Sqr(k13x)/dksq[0]);
	      n13z = rint(Sqr(k13z)/dksq[2]);
#else
	      n13x = ((int) (Sqr(k13x)/dksq[0]));
	      n13z = ((int) (Sqr(k13z)/dksq[2]));
#endif
#else
#ifdef USE_RINT
	      n13x = rint(fabs(k13x)/dk[0]);
	      n13z = rint(fabs(k13z)/dk[2]);
#else
	      n13x = ((int) (fabs(k13x)/dk[0]));
	      n13z = ((int) (fabs(k13z)/dk[2]));
#endif
#endif
#if 0
      	      if (nv[0]==15 && nv[1]==0 && nv[2]==0)
      		{
      		  printf("*** n13x=%d n13z=%d\n", n13x, n13z);
      		  printf("*** k13x=%f k13z=%f\n", Sqr(k13x), Sqr(k13z));
      		}
      	      if ( rint(Sqr(k13x)) < 5 && rint(Sqr(k13z)) > 255 && ((int)Sqr(k13z))<260)
      		{
		  printf("boh n13x=%d n13z=%d\n", n13x, n13z);
		  printf("k1=%f k3=%f\n", (n13x)*dksq[0], (n13z)*dksq[2]);
		  printf("k13x=%f k13z=%f\n", Sqr(k13x), Sqr(k13z));
		}
#endif
	      if (n13x < 0 || n13z < 0)
		{
		  //printf("BOH\n");
		}
	      if (n13x < npoints[0] && n13z < npoints[2] && n13x >= 0 && n13z >=0 )
		{
		  calc_fourier_Q(kv);
//printf("qui2 n13x=%d n13z=%d reQ[0][2]=%f imQ[0][2]=%f\n", n13x, n13z, Sqr(reQ[0][2])+Sqr(imQ[0][2]),Sqr(reQ[1][2])+Sqr(imQ[1][2]));
#ifdef TRANSFORM_Q
		  for (k1=0; k1 < 3; k1++)
		    {
		      Rl[0][k1] = vecx[k1];
		      Rl[1][k1] = vecy[k1];
		      Rl[2][k1] = vecz[k1];
		      //if (isnan(cos(R[0][k1]))||isnan(R[1][k1])||isnan(R[2][k1]))
		      //printf("mmmah\n");
		    }   
		  for (k1 = 0; k1 < 3; k1++)
		    for (k2 = 0; k2 < 3; k2++)
		      {
			reQ13[k1][k2] = 0;
			imQ13[k1][k2] = 0;
			for (k3 = 0; k3 < 3; k3++)
			  //Ro[k1][k2] += M[k1][k3]*R[i][k3][k2];
			  {
			    reQ13[k1][k2] += reQ[k1][k3]*Rl[k2][k3];
			    imQ13[k1][k2] += imQ[k1][k3]*Rl[k2][k3];
			  }
		      }
		  for (k1 = 0; k1 < 3; k1++)
		    for (k2 = 0; k2 < 3; k2++)
		      {
			reQ[k1][k2] = reQ13[k1][k2];
			imQ[k1][k2] = imQ13[k1][k2];
		      } 
		  for (k1 = 0; k1 < 3; k1++)
		    for (k2 = 0; k2 < 3; k2++)
		      {
			reQ13[k1][k2] = 0.0;
			imQ13[k1][k2] = 0.0;
			for (k3 = 0; k3 < 3; k3++)
			  //Ro[k1][k2] += M[k1][k3]*R[i][k3][k2];
			  {
			    reQ13[k1][k2] += Rl[k1][k3]*reQ[k3][k2];
			    imQ13[k1][k2] += Rl[k1][k3]*imQ[k3][k2];
			  }
		      }

		  accQsq13[n13x][n13z] += Sqr(reQ13[0][2])+Sqr(imQ13[0][2]);
		  accQsq23[n13x][n13z] += Sqr(reQ13[1][2])+Sqr(imQ13[1][2]);
		  nbins[n13x][n13z] += 1.0;
#else
		  accQsq13[n13x][n13z] += Sqr(reQ[0][2])+Sqr(imQ[0][2]);
		  accQsq23[n13x][n13z] += Sqr(reQ[1][2])+Sqr(imQ[1][2]);
		  nbins[n13x][n13z] += 1.0;
#endif
		  //printf("n13x=%d n13z=%d acc13=%f acc23=%f\n", n13x, n13z, accQsq13[n13x][n13z], accQsq23[n13x][n13z]);
#if 0
		  if (isnan(accQsq13[n13x][n13z]))
		    {
		      printf("reQ=%f imQ=%f\n",reQ[0][2], imQ[0][2]);
		    }
#endif
//printf("qui3\n");
		}
	      //printf("acc13=%f acc23=%f\n", accQsq13[119][49], accQsq23[119][49]);
	    }
    }
  fclose(f2);	
  //printf("eccolo acc13=%f acc23=%f\n", accQsq13[119][49], accQsq23[119][49]);
  for (nv[0] = 0; nv[0] < npoints[0]; (nv[0])++)
    for (nv[2] = 0; nv[2] < npoints[2]; (nv[2])++)
      {
	if (nbins[nv[0]][nv[2]] > 0)
	  {
	    accQsq13[nv[0]][nv[2]] /= nbins[nv[0]][nv[2]];
	    accQsq23[nv[0]][nv[2]] /= nbins[nv[0]][nv[2]];
	  }
       	//printf("eccolo acc13=%f acc23=%f\n", accQsq13[nv[0]][nv[2]], accQsq23[nv[0]][nv[2]]);
	//if (nv[0]==119&&nv[2]==49)
	//printf("eccolo2 acc13=%f acc23=%f\n", accQsq13[119][49], accQsq23[119][49]);
      }

  //printf("eccolo3 acc13=%f acc23=%f\n", accQsq13[119][49], accQsq23[119][49]);
  f = fopen("Welconst.dat", "w+");
  S/=((double)nf);
  printf("S=%f\n", S);
  for (nv[0] = 0; nv[0] < npoints[0]; (nv[0])++)
    for (nv[2] = 0; nv[2] < npoints[2]; (nv[2])++)
      {
#ifdef SQ_BINNING
#ifdef USE_RINT
	k13sqx = nv[0]*dksq[0];
	k13sqz = nv[2]*dksq[2];
#else
	k13sqx = (nv[0]+0.5)*dksq[0];
	k13sqz = (nv[2]+0.5)*dksq[2];
#endif
	if (nbins[nv[0]][nv[2]] > minbins)
	  fprintf(f, "%f %f %f %f %.0f\n", k13sqx, k13sqz, Sqr(S)*V*9.0/accQsq13[nv[0]][nv[2]]/4.0, Sqr(S)*V*9.0/accQsq23[nv[0]][nv[2]]/4.0, nbins[nv[0]][nv[2]]);
#else
#ifdef USE_RINT
	k13x = nv[0]*dk[0];
	k13z = nv[2]*dk[2];
#else
	k13x = (nv[0]+0.5)*dk[0];
	k13z = (nv[2]+0.5)*dk[2];
#endif
	if (nbins[nv[0]][nv[2]] > minbins)
	  fprintf(f, "%f %f %f %f %.0f\n", Sqr(k13x), Sqr(k13z), Sqr(S)*V*9.0/accQsq13[nv[0]][nv[2]]/4.0, Sqr(S)*V*9.0/accQsq23[nv[0]][nv[2]]/4.0, nbins[nv[0]][nv[2]]);
#endif
	//else 
	  //printf("nv[0]=%d nv[2]=%d k13x=%f k13z=%f is zero!\n", nv[0], nv[2], k13x, k13z);
      }  
  fclose(f);
  return 1;
}
