#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define Sqr(x) ((x)*(x))
double lenHC=4.0, surfRad=10.0;
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );
char line[1000000], parname[124], parval[1000000];
char dummy[2048];
int mcsim=0, cubic_box=1, nonem=0, calcnv=0, no2d=0;
int N, particles_type=1, k1, k2;
double *x[3], L, ti, *w[3], storerate; 
double Lx, Ly, Lz;
double nv[3], vecx[3], vecy[3], vecz[3], *R[3][3];
double *DR[3], deltaAA=-1.0, deltaBB=-1.0, deltaAB=-1.0, sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0, 
       Dr, theta, sigmaSticky=-1.0, sa[2], sb[2], sc[2], maxax0, maxax1, maxax, maxsax, maxsaxAA, maxsaxAB, maxsaxBB;
char fname[1024], inputfile[1024];
int eventDriven=0;
int points;
int foundDRs=0, foundrot=0;

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

	     // printf("r[%d]=%f %f %f\n", i, r[0][i], r[1][i], r[2][i]);
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
  printf("calcgr [ --no2d/-n2 ] [--lenhc/-lhc] [--nonem/-nn] [--calcnemvec/-cv] [--mcsim/-mc] [--isoav/-ia] [--nemvector/-nv (x,y,z) ] [-gp/-gnuplot] <confs_file> [points]\n");
  exit(0);
}
double threshold=0.05;
int gnuplot=0, isoav=0;
void parse_param(int argc, char** argv)
{
  int extraparam=0;  
  int cc=1;
  points=100;
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
      else if (!strcmp(argv[cc],"--no2d")||!strcmp(argv[cc],"-n2"))
	{
	  no2d=1;
	}
      else if (!strcmp(argv[cc],"--gnuplot")||!strcmp(argv[cc],"-gp"))
	{
	  gnuplot=1;
	}
      else if (!strcmp(argv[cc],"--isoav")||!strcmp(argv[cc],"-ia"))
	{
	  isoav=1;
	}
      else if (!strcmp(argv[cc],"--nonem")||!strcmp(argv[cc],"-nn"))
	{
	  nonem=1;
	}
      else if (!strcmp(argv[cc],"--calcnemvec")||!strcmp(argv[cc],"-cv"))
	{
	  calcnv=1;
	}
      else if (!strcmp(argv[cc],"--lenhc")||!strcmp(argv[cc],"-lhc"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  lenHC = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--nemvector")||!strcmp(argv[cc],"-nv"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  sscanf(argv[cc],"(%lf,%lf,%lf)", &(nv[0]), &(nv[1]), &(nv[2]));
	}
      else if (cc==argc|| extraparam==2)
	print_usage();
      else if (extraparam == 0)
	{
	  extraparam++;
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam==1)
	{
	  extraparam++;
	  points=atoi(argv[cc]);
	}
      else 
	print_usage();
      cc++;
    }
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
    } 
  for (k1=0; k1 < 3; k1++)
    {
      vp[k1] = 0.0;
      for (k2=0; k2 < 3; k2++)
	vp[k1] += R[k1][k2]*v[k2];
    }
//  printf("vp=%f %f %f\n", vp[0], vp[1], vp[2]);
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
double S=0, Q[3][3];
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
void build_ref_system(void)
{
  int k;
  double sp, norm;
  /* build the nematic reference system */
  vecz[0] = nv[0];
  vecz[1] = nv[1];
  vecz[2] = nv[2];
  //printf("delr=%.15G nematic=%f %f %f\n", delr, nv[0], nv[1], nv[2]);
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
  vectProdVec(vecz, vecx, vecy);
#if 1
  printf("normx=%f\n", sqrt(Sqr(vecx[0])+Sqr(vecx[1])+Sqr(vecx[2])));
  printf("{{%f,%f,%f},\n", vecx[0], vecx[1], vecx[2]);
  printf("{%f,%f,%f},\n", vecy[0], vecy[1], vecy[2]);
  printf("{%f,%f,%f}}\n", vecz[0], vecz[1], vecz[2]);
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
    nv[a] = eigvec[numev][a];
}
double pos1[3], pos2[3], uhc[3], pos1Nem[3], pos2Nem[3], numsd=0.0;
int main(int argc, char** argv)
{
  FILE *f, *f2, *f1;
  int binx, biny, binz;
  int k, nf, i, a, b, nat, NN, j, ii, bin, binParaAv, binPerpAv;
  double rx, ry, rz, norm, g0para, g0perp, minpara, maxpara, minperp, maxperp;
  double r, delr, tref=0.0, Dx[3], DxNem[3], *g0, **g0Perp, **g0Parall, *g0ParaAv, *g0PerpAv,
	 g0m, distSq, rlower, rupper, cost, nIdeal, costParaAv, costPerpAv;
  double sp, time, refTime, RCUT;
  int iZ, jZ, iX, jX, iY, jY, NP1, NP2;
  double shift[3];
  double ene=0.0;

#if 0
  double g2m, g4m, g6m;
  double *g2, *g4, *g6;
#endif
  parse_param(argc, argv);
  f2 = fopen(inputfile,"r");
  fscanf(f2, "%[^\n]\n", fname);
  pi = acos(0)*2.0;
  printf("pi=%f\n", pi);
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

  /* max sus window */
  NP+=100;
  NPA+=100;
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
    printf("L=%.15G\n", L);
  else
    printf("Lx=%.15G Ly=%.15G Lz=%.15G\n", Lx, Ly, Lz);
  if (!nonem)
    {
      if (!no2d)
	{
	  g0Perp = malloc(sizeof(double*)*4*points);
	  g0Parall= malloc(sizeof(double*)*4*points);
	}

      for (k1=0; k1 < 4*points; k1++)
	{
	  g0Perp[k1] = malloc(sizeof(double)*points*4);
	  g0Parall[k1] = malloc(sizeof(double)*points*4);
	}
    }
  g0 = malloc(sizeof(double)*points*4);
  if (isoav)
    {
      g0ParaAv = malloc(sizeof(double)*points*4);
      g0PerpAv = malloc(sizeof(double)*points*4);
      //normPerpAv = malloc(sizeof(double)*points*4);
      //normParaAv = malloc(sizeof(double)*points*4);
    }
  for (k1=0; k1 < 4*points; k1++)
    {
      g0[k1] = 0.0;
      if (isoav)
	{
	  g0ParaAv[k1] = 0;
	  g0PerpAv[k1] = 0;
	  //normPerpAv[k1] = normParaAv[k1] = 0.0;
 
	}
      if (!nonem && !no2d)
	{
	  for (k2=0; k2 < 4*points; k2++)
	    {
	      g0Perp[k1][k2] = 0.0;
	      g0Parall[k1][k2] = 0.0;
	    }
	}
    }
  for (ii = 0; ii < 3; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*NP);
    }
#if 0
  g2 = malloc(sizeof(double)*points);
  g4 = malloc(sizeof(double)*points);
  g6 = malloc(sizeof(double)*points);
#endif
  if (cubic_box)
    {
      delr = L / 2.0 / ((double)points);
    }
  else
    {
      delr = max3(Lx,Ly,Lz)/2.0/((double)points);
    }
  rewind(f2);
  nf = 0;

  if (!nonem && calcnv==0)
    {
      /* build the nematic reference system */
      vecz[0] = nv[0];
      vecz[1] = nv[1];
      vecz[2] = nv[2];
      printf("delr=%.15G nematic=%f %f %f\n", delr, nv[0], nv[1], nv[2]);
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
      vectProdVec(vecz, vecx, vecy);
      printf("{{%f,%f,%f},\n", vecx[0], vecx[1], vecx[2]);
      printf("{%f,%f,%f},\n", vecy[0], vecy[1], vecy[2]);
      printf("{%f,%f,%f}}\n", vecz[0], vecz[1], vecz[2]);
    }
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      nf++;
      NPA = -1;
      readconf(fname, &time, &refTime, &NP, &NPA, x, w, DR);
      if (calcnv)
	{
	  calc_nem_vec();
	  build_ref_system();
	}
      //printf("NP=%d NPA=%d\n", NP, NPA);
#if 0
      for (i=0; i < NP; i++)
	{
	  for (a=0; a < 3; a++)
	    {  uhc[a] = R[0][a][i];
	      pos1[a]=x[a][i] + uhc[a]*lenHC/2.0; 
	      pos2[a]=x[a][i] - uhc[a]*lenHC/2.0;
	    }
  	  lab2nem(pos1,pos1Nem);
	  lab2nem(pos2,pos2Nem);
	  if (pos1Nem[2]*pos2Nem[2] < 0.0 && (Sqr(pos1Nem[0])+Sqr(pos2Nem[1]) < Sqr(surfRad)))
	    numsd++;
	}
#endif
      for (i=0; i < NP-1; i++)
	for (j = i+1; j < NP; j++)
	  {
	    for (a = 0; a < 3; a++)
	      {
		Dx[a] = x[a][i] - x[a][j];
	      }

	    if (cubic_box)
	      {
		for (a = 0; a < 3; a++)
		  Dx[a] = Dx[a] - L * rint(Dx[a]/L);
	      }
	    else
	      {
		Dx[0] = Dx[0] - Lx * rint(Dx[0]/Lx);
	      	Dx[1] = Dx[1] - Ly * rint(Dx[1]/Ly);
      		Dx[2] = Dx[2] - Lz * rint(Dx[2]/Lz);
	      }
	    if (!nonem)
	      {
		lab2nem(Dx,DxNem);
		//printf("Dxnem[%d]=%f %f %f\n",i, DxNem[0], DxNem[1], DxNem[2] );
		//printf("Dx   [%d]=%f %f %f\n",i, Dx[0], Dx[1], Dx[2] );
	      }	      
	   //printf("DxNem=%f %f %f\n", DxNem[0], DxNem[1], DxNem[2]);
	    distSq = 0.0;
	    for (a = 0; a < 3; a++)
	      distSq += Sqr(Dx[a]);
	    if (!nonem)
	      {
		if (cubic_box)
		  {
		    binx = (int) floor(DxNem[0] / delr);
		    biny = (int) floor(DxNem[1] / delr);
		    binz = (int) floor(DxNem[2] / delr);
		  }
		else
		  {
		    binx = (int) floor(DxNem[0] / delr);
		    biny = (int) floor(DxNem[1] / delr);
		    binz = (int) floor(DxNem[2] / delr);
		  }
	      }
	    bin = (int) (sqrt(distSq)/delr);
	    //printf("(%d-%d) bin=%d\n", i, j, bin);
	    //printf("(%d-%d) %d %d %d\n", i, j, binx, biny, binz);
	    //printf("DrNemSq=%f Dr=%f\n", calc_norm(DxNem), calc_norm(Dx));
#if 0
	    if (binx < -2*points )
	      printf("boh=%d\n", binx+2*points);
#endif
	    if (!nonem)
	      {
		if (binx < 2*points && biny < 2*points && binz < 2*points && 
		    binx >= -2*points && biny >= -2*points  && binz >= -2*points)
		  {
		    if (isoav && (binx==0 || binx==-1) && (biny==0 || biny==-1))
		      {
			binParaAv = (int) (fabs(DxNem[2])/delr);
    			g0ParaAv[binParaAv]+=2.0;
		      } 
		    if (binx==0 || binx==-1)
		      {
			if (!no2d)
			  g0Parall[biny+2*points][binz+2*points] += 2.0;
		      }
		    if (binz==0 || binz==-1)
		      {
			if (!no2d)
			  g0Perp[binx+2*points][biny+2*points] += 2.0;
			//if (binx < 5)	
			//  printf("binx=%d biny=%d\n", binx, biny);
			if (isoav)
			  {
			    binPerpAv = (int) (sqrt(Sqr(DxNem[0])+Sqr(DxNem[1]))/delr);
			    g0PerpAv[binPerpAv]+=2.0;
 			  }
		      }
		    //printf("g0[%d]=%.15G\n", bin, g0[bin]);
		  }
	      }
    	    if (bin >=0 && bin < 4*points)
	      g0[bin] += 2.0;
		//printf("bin=%d\n", bin);
	    //exit(1);
	      
	  }
    }
  fclose(f2); 
  f = fopen("gr.dat", "w+");
  if (isoav)
    {
      f1 = fopen("grParaAv.dat", "w+");
      f2 = fopen("grPerpAv.dat", "w+");
    }
  r = delr*0.5;
  if (cubic_box)
    {
      costPerpAv = (L*L*L)/((double)NP)/((double)NP);// /(delr*delr*delr);
      costParaAv = (L*L*L)/((double)NP)/((double)NP)/(delr*delr*delr);
    }
  else
    {
      costPerpAv = (Lx*Ly*Lz)/((double)NP)/((double)NP);// /(delr*delr*delr);
      costParaAv = (Lx*Ly*Lz)/((double)NP)/((double)NP)/(delr*delr*delr);
    }	
  
  costParaAv /= 8.0; /* 4 poichè ho considerato un "tubo" 2x2 ossia binx=-1,0 e biny=-1,0 inoltre
			sopra e sotto sono equivalenti (vedi il fabs() nel loop sopra) quindi un altro due */
  costPerpAv /= 2.0;
  printf("Lx=%f NP=%d delr=%f\n", Lx, NP, delr);
  if (isoav)
    {
      for (k1 = 0; k1 < 2*points; k1++)
	{
	  rlower = ( (double) k1) * delr;
	  rupper = rlower + delr;
	  //nIdeal =  pi * delr * (Sqr(rupper) - Sqr(rlower));
	  nIdeal =  pi * delr * (Sqr(rupper) - Sqr(rlower));
	  g0PerpAv[k1] *= costPerpAv / nIdeal / ((double)nf);
	  g0ParaAv[k1] *= costParaAv / ((double) nf);
	  //printf("normPerpAv[%d]=%f\n", k1, normPerpAv[k1]);
	}   
    }

  if (cubic_box)
    {
      cost = 4.0 * pi * NP / 3.0 / (L*L*L);
    }
  else
    {
      cost = 4.0 * pi * NP / 3.0 / (Lx*Ly*Lz);
    }
  for (ii = 0; ii < 2*points; ii++)
    {
      rlower = ( (double) ii) * delr;
      rupper = rlower + delr;
      nIdeal = cost * (Sqr(rupper)*rupper - Sqr(rlower)*rlower);
      //printf("nf=%d nIdeal=%.15G g0[%d]=%.15G\n", nf, nIdeal, ii, g0[ii]);
      g0m = g0[ii]/((double)nf)/((double)NP)/nIdeal;
#if 0
      g2m = (3.0*g2[ii]/cc[ii] - 1.0)/2.0;
      g4m = (35.0*g4[ii]/cc[ii] - 30.0*g2[ii]/cc[ii] + 3.0) / 8.0;
      g6m = (231.0*g6[ii]/cc[ii] - 315.0*g4[ii]/cc[ii] + 105.0*g2[ii]/cc[ii] - 5.0)/16.0;
      fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", r, g0m, g2m, g4m, g6m);
#endif
      fprintf(f, "%.15G %.15G %.15G\n", r, g0m, g0[ii]);
      if (isoav)
	{
	  fprintf(f1, "%.15G %.15G\n", r, g0ParaAv[ii]);
	  fprintf(f2, "%.15G %.15G\n", r, g0PerpAv[ii]);
	}
      r += delr;
      if (cubic_box)
	{
	  if (r > L*0.5)
	    break;
	}
      else
	{
	  if (r > max3(Lx,Ly,Lz)*0.5) 
	    break;
	}
    }
  fclose(f);
  if (isoav)
    {
      fclose(f1);
      fclose(f2);
    }
  if (nonem)
    return 0;
  if (!no2d)
    {
      f1 = fopen("grpara.dat", "w+");
      f2 = fopen("grperp.dat", "w+");
      r = delr*0.5;
      if (cubic_box)
	cost = (L*L*L)/((double)NP)/((double)NP)/(delr*delr*delr);
      else
	cost = (Lx*Ly*Lz)/((double)NP)/((double)NP)/(delr*delr*delr);
      cost /= 2; /* N.B. divido per due poiché sto considerando in realtà due piani: uno appena
		    sopra e uno appena sotto (z=0 o x=0) */

      for (k1 = 0; k1 < 4*points; k1++)
	{
	  for (k2 = 0; k2 < 4*points; k2++)
	    {
	      //printf("nf=%d nIdeal=%.15G g0[%d]=%.15G\n", nf, nIdeal, ii, g0[ii]);
	      g0perp = cost*g0Perp[k1][k2]/((double)nf);
	      g0para = cost*g0Parall[k1][k2]/((double)nf);
#if 0
	      g2m = (3.0*g2[ii]/cc[ii] - 1.0)/2.0;
	      g4m = (35.0*g4[ii]/cc[ii] - 30.0*g2[ii]/cc[ii] + 3.0) / 8.0;
	      g6m = (231.0*g6[ii]/cc[ii] - 315.0*g4[ii]/cc[ii] + 105.0*g2[ii]/cc[ii] - 5.0)/16.0;
	      fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", r, g0m, g2m, g4m, g6m);
#endif
	      if (cubic_box)
		{
		  rx = (((double)k1)-2*points)*delr;
		  ry = (((double)k2)-2*points)*delr;
		}
	      else
		{
		  rx = (((double)k1)-2*points)*delr;
		  ry = (((double)k2)-2*points)*delr;
		}
	      fprintf(f1, "%.15G %.15G %.15G\n", rx, ry, g0para);
	      fprintf(f2, "%.15G %.15G %.15G\n", rx, ry, g0perp);
	      if (k2==0 && k1==0)
		{
		  minpara = maxpara = g0para;
		  minperp = maxperp =g0perp;
		}
	      else
		{
		  if (g0para < minpara)
		    minpara = g0para;
		  if (g0para > maxpara)
		    maxpara = g0para;
		  if (g0perp < minperp)
		    minperp = g0perp;
		  if (g0perp > maxperp)
		    maxperp = g0perp;
		}

	    }
	  if (gnuplot && k1 < 4*points-1)
	    {
	      fprintf(f1, "\n");
	      fprintf(f2, "\n");
	    }
	}
      //printf("gnuplot=%d\n", gnuplot);
      fclose(f1);
      fclose(f2);
    }
  printf("Parallel Range [%.15G:%.15G]\n", minpara, maxpara);
  printf("Perpendicular Range [%.15G:%.15G]\n", minperp, maxperp);
  printf("Average S=%f\n", S/((double)nf));
  //printf("surface number density=%f\n", numsd/3.1415926535897932385/Sqr(surfRad)/((double)nf));
  return 0;
}

