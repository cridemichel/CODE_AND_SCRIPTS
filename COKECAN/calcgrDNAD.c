#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define Sqr(x) ((x)*(x))
char line[1000000], parname[124], parval[1000000];
char dummy[2048];
int N, particles_type=1, k1, k2;
double *x[3], L, ti, *w[3], storerate;
double nv[3], vecx[3], vecy[3], vecz[3];
double *DR[3], deltaAA=-1.0, deltaBB=-1.0, deltaAB=-1.0, sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0, 
       Dr, theta, sigmaSticky=-1.0, sa[2], sb[2], sc[2], maxax0, maxax1, maxax, maxsax, maxsaxAA, maxsaxAB, maxsaxBB;
char fname[1024], inputfile[1024];
int eventDriven=0;
int points;
int foundDRs=0, foundrot=0;
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double *DR[3])
{
  FILE *f;
  int nat=0, i, cpos;
  double dt=-1;
  int curstp=-1;
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
	  if (!strcmp(parname,"DR"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf %lf %lf ", &DR[0][i], &DR[1][i], &DR[2][i]);
		}
	      foundDRs = 1;
	    }
	  else if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	      foundrot = 1;
	    }
	  else if (!strcmp(parname,"sumoy"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[1][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoz"))
	    {
	      for (i=0; i < NP; i++)
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
	  fseek(f, cpos, SEEK_SET);
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%[^\n]\n", line); 
	      if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
		{
		  sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
		}
	      //printf("r[%d]=%f %f %f\n", i, r[0][i], r[1][i], r[2][i]);
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
  printf("calcgr [--nemvector/-nv (x,y,z) ] <confs_file> [points]\n");
  exit(0);
}
double threshold=0.05;
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
      else if (!strcmp(argv[cc],"--threshold")||!strcmp(argv[cc],"-t"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  threshold=atof(argv[cc]);
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
      R[k1][0] = vecx[k1];
      R[k1][1] = vecy[k1];
      R[k1][2] = vecz[k1];
    } 
  for (k1=0; k1 < 3; k1++)
    {
      vp[k1] = 0.0;
      for (k2=0; k2 < 3; k2++)
	vp[k1] += R[k2][k1]*v[k2];
    }
}
int main(int argc, char** argv)
{
  FILE *f, *f2, *f1;
  int binx, biny, binz;
  int k, nf, i, a, b, nat, NN, j, ii, bin;
  double rx, ry, rz, norm, g0para, g0perp;
  double r, delr, tref=0.0, Dx[3], DxNem[3], *g0, **g0Perp, **g0Parall, g0m, distSq, rlower, rupper, cost, nIdeal;
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
	      for (i=0; i < 2*NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      fscanf(f, "%lf\n", &L);
	      break;
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
  printf("L=%.15G\n", L);
  g0Perp = malloc(sizeof(double*)*points);
  g0Parall= malloc(sizeof(double*)*points);

  for (k1=0; k1 < points; k1++)
    {
      g0Perp[k1] = malloc(sizeof(double)*points);
      g0Parall[k1] = malloc(sizeof(double)*points);
    }
  g0 = malloc(sizeof(double)*points);
  for (k1=0; k1 < points; k1++)
    {
      g0[k1] = 0.0;
      for (k2=0; k2 < points; k2++)
	{
	  g0Perp[k1][k2] = 0.0;
	  g0Parall[k1][k2] = 0.0;
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
  delr = L / 2.0 / ((double)points);
  rewind(f2);
  nf = 0;
  /* build the nematic reference system */
  vecz[0] = nv[0];
  vecz[1] = nv[1];
  vecz[2] = nv[2];
  printf("nematic=%f %f %f\n", nv[0], nv[1], nv[2]);
  vecx[0] = 1.0;
  vecx[1] = 0;
  vecx[2] = 0;
  if (vecx[0] == vecz[0] && vecx[1]==vecz[1] && vecx[2]==vecz[2])
    {
      vecx[0] = 0.0;
      vecx[1] = 1.0;
      vecx[2] = 0.0;
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
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      nf++;
      readconf(fname, &time, &refTime, NP, x, w, DR);
      for (i=0; i < NP-1; i++)
	for (j = i+1; j < NP; j++)
	  {
	    for (a = 0; a < 3; a++)
	      {
		Dx[a] = x[a][i] - x[a][j];
	      }
	    for (a = 0; a < 3; a++)
	      Dx[a] = Dx[a] - L * rint(Dx[a]/L);
	    lab2nem(Dx,DxNem);
	    distSq = 0.0;
	    for (a = 0; a < 3; a++)
	      distSq += Sqr(Dx[a]);
	    binx = (int) (fabs(DxNem[0]+L*0.5) / delr);
	    biny = (int) (fabs(DxNem[1]+L*0.5) / delr);
	    binz = (int) (fabs(DxNem[2]+L*0.5) / delr);
	    bin = (int) (sqrt(distSq)/delr);
	    //printf("(%d-%d) bin=%d\n", i, j, bin);
	    if (binx < points && biny < points && binz < points && 
		binx >= 0 && biny >= 0  && binz >=0)
	      {
		if (binx==0)
		  g0Parall[biny][binz] += 2.0;
		else if (binz==0)
		  g0Perp[binx][biny] += 2.0;
		//printf("g0[%d]=%.15G\n", bin, g0[bin]);
	      }
	    if (bin >=0 && bin < points)
	      g0[bin] += 2.0;
    	    //printf("bin=%d\n", bin);
	    //exit(1);
	      
	  }
    }
  fclose(f2); 
  f = fopen("gr.dat", "w+");
  r = delr*0.5;
  cost = 4.0 * pi * NP / 3.0 / (L*L*L);
  for (ii = 0; ii < points; ii++)
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
      r += delr;
    }
  fclose(f);
  f1 = fopen("grpara.dat", "w+");
  f2 = fopen("grperp.dat", "w+");
  r = delr*0.5;
  cost = (L*L*L)/((double)NP)/((double)NP);
  for (k1 = 0; k1 < points; k1++)
    {
      for (k2 = 0; k2 < points; k2++)
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
	  rx = (((double)k1)+0.5)*delr;
	  ry = (((double)k2)+0.5)*delr;
	  fprintf(f1, "%.15G %.15G %.15G\n", rx, ry, g0para);
	  fprintf(f2, "%.15G %.15G %.15G\n", rx, ry, g0perp);
	}
    }
  fclose(f1);
  fclose(f2);
  return 0;
}

