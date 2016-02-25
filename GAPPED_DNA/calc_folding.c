#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define Sqr(x) ((x)*(x))
double lenHC=4.0, surfRad=10.0;
char line[1000000], parname[124], parval[1000000];
char dummy[2048];
int mcsim=0, cubic_box=1, nonem=0, calcnv=0, no2d=0;
int N, particles_type=1, k1, k2;
double *x[3], L, ti, *w[3], storerate, *angdistr; 
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
  printf("calc_folding [--mcsim/-mc] [--threshold/-th <threshold>] <confs_file> [points]\n");
  exit(0);
}
double threshold=30.0;
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
      else if (!strcmp(argv[cc],"--threshold")||!strcmp(argv[cc],"-tr"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  threshold = atof(argv[cc]);
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
#if 0
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
#endif
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
#if 0
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
#endif
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
double pos1[3], pos2[3], uhc[3], pos1Nem[3], pos2Nem[3], numsd=0.0;
int main(int argc, char** argv)
{
  FILE *f, *f2, *f1;
  int binx, biny, binz;
  int k, nf, i, a, b, nat, NN, j, ii, bin, n, kk;
  double angle, rx, ry, rz, norm, dang, sum;
  double ni[3], nj[3];
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

  for (a = 0; a < 3; a++)
    {
      x[a] = malloc(sizeof(double)*NP);
      w[a] = malloc(sizeof(double)*NP);
      /* allocate orientations */
      for (b=0; b < 3; b++)      
	R[a][b] = malloc(sizeof(double)*NP);
    }
  for (ii = 0; ii < 3; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*NP);
    }
  angdistr = malloc(sizeof(double)*points);

  if (cubic_box)
    printf("L=%.15G\n", L);
  else
    printf("Lx=%.15G Ly=%.15G Lz=%.15G\n", Lx, Ly, Lz);
  dang = 180.0/((double)points);
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      nf++;
      NPA = -1;
      readconf(fname, &time, &refTime, &NP, &NPA, x, w, DR);
      for (i=0; i < NP; i+=2)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      ni[kk] = R[0][kk][i];
	      nj[kk] = R[0][kk][i+1];
	    }
	  angle=180.0*acos(scalProd(ni,nj))/M_PI;
	  bin = (int) (angle/dang);
	  if (bin>=0 && bin < points)
	    (angdistr[bin])++;
	  //printf("angle=%f bin=%d M_PI=%f\n", angle, bin, M_PI);
	}
    }
  fclose(f2); 
  f1 = fopen("angdistr.dat", "w+");
  for (n=0; n < points; n++)
    {
      angdistr[n] = angdistr[n]/((double)NP)/2.0/((double)nf);
    }
  sum = 0.0; 
  for (n=0; n < points; n++)
    {
      sum += angdistr[n]*dang;
    } 
  for (n=0; n < points; n++)
    {
      angdistr[n] /= sum;
    }
  for (n=0; n < points; n++)
    {
      //printf("%f %.15G f=%p\n", n*dang-90.0, angdistr[n], f);
      fprintf(f1, "%f %.15G\n", n*dang, angdistr[n]);
    }
  fclose(f1);  
  sum = 0.0;
  for (n=0; n < points; n++)
    {
      if (n*dang < threshold)
	sum += angdistr[n]*dang;
    }
  f1 = fopen("foldfract.dat", "w+");
  fprintf(f1, "%.15G\n", sum);
  fclose(f1);
  return 0;
}

