#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define Sqr(x) ((x)*(x))
char line[1000000], parname[124], parval[1000000];
char dummy[2048];
int mcsim=0, calciso=0;
int N, particles_type=1, k1, k2;
double *x[3], L, ti, *w[3], storerate;
double rv[3], vecx[3], vecy[3], vecz[3], *u[3];
double *DR[3], deltaAA=-1.0, deltaBB=-1.0, deltaAB=-1.0, sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0, 
       Dr, theta, sigmaSticky=-1.0, sa[2], sb[2], sc[2], maxax0, maxax1, maxax, maxsax, maxsaxAA, maxsaxAB, maxsaxBB;
char fname[1024], inputfile[1024];
int eventDriven=0;
int points, angpoints=-1, perppoints=10;
int foundDRs=0, foundrot=0;
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double *DR[3], double *u[3])
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
	      if (!sscanf(line, "%lf %lf %lf %lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i], &u[0][i], &u[1][i],
			 &u[2][i])==6)
		{
		  sscanf(line, "%lf %lf %lf %lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], &u[0][i], &u[1][i],
			 &u[2][i], dummy); 
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
double calc_norm(double *vec);
void vectProdVec(double *A, double *B, double *C);
void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;

  /* first row vector */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 1; u[1] = 1; u[2] = 1;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = -1; u[1] = -1; u[2] = 1;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[0][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[0][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
#if 0
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    Rt[k1][k2]=R[k2][k1];
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    R[k1][k2]=Rt[k1][k2];
#endif

  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
}

void print_usage(void)
{
  printf("calcgr [ --calciso/-ci <0|1> ] [ --perppoints/-pp <number> ] [ --angpoints/ap <number>] [--mcsim/-mc] [--rversor/-rv (x,y,z) ] <confs_file> [points]\n");
  exit(0); 
}
double threshold=0.05;
int gnuplot=0;
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
      //printf("argc=%d argv[argc]=%s cc=%d\n", argc, argv[cc], cc);
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
       else if (!strcmp(argv[cc],"--mcsim")||!strcmp(argv[cc],"-mc"))
	{
	  mcsim=1;
	}
      //else if (!strcmp(argv[cc],"--gnuplot")||!strcmp(argv[cc],"-gp"))
	//{
	  //gnuplot=1;
	//}
      else if (!strcmp(argv[cc], "--calciso")||!strcmp(argv[cc],"-ci"))
	{
	  calciso = 1;
	} 
      else if (!strcmp(argv[cc],"--rversor")||!strcmp(argv[cc],"-rv"))
	{
	  /* rv è la direzione lungo la quale calcolare la g(r,\Omega_12) vettoriale nel 
	     riferimento del corpo rigido (cioè quello individuato dai suoi assi principali) */
	  cc++;
	  if (cc==argc)
	    print_usage();
	  sscanf(argv[cc],"(%lf,%lf,%lf)", &(rv[0]), &(rv[1]), &(rv[2]));
	}
      else if (!strcmp(argv[cc],"--angpoints")||!strcmp(argv[cc],"-ap"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  angpoints=atoi(argv[cc]);
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
void lab2body(int i, double v[3], double vp[3])
{
  int k1, k2;
  double R[3][3], uL[3];

  for (k1=0; k1 < 3; k1++)
    uL[k1] = u[k1][i];
  versor_to_R(uL[0], uL[1], uL[2], R);
  for (k1=0; k1 < 3; k1++)
    {
      vp[k1] = 0.0;
      for (k2=0; k2 < 3; k2++)
	vp[k1] += R[k1][k2]*v[k2];
    }
//  printf("vp=%f %f %f\n", vp[0], vp[1], vp[2]);
  //printf("v=%f %f %f\n", v[0], v[1], v[2]);
}
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}

void projectVec(double DxP[3], double vec[3], int *binx, int *biny, int *binz, double delr, double delrp)
{
  double Rl[3][3], vec2[3];
  int kk;
  /* il vettore vec viene assunto come asse lungo cui calcolare la g(r) vettoriale */
  *binx = (int) floor(scalProd(DxP, vec)/delr);
  versor_to_R(vec[0], vec[1], vec[2], Rl);
  vectProdVec(DxP, vec, vec2);
  *biny = (int) floor(scalProd(vec2, Rl[1])/delrp);
  *binz = (int) floor(scalProd(vec2, Rl[2])/delrp);  
}
int main(int argc, char** argv)
{
  FILE *f, *f2, *f1;
  int jj, bina, binx, biny, binz;
  int k, nf, i, a, b, nat, NN, j, ii, bin;
  double delrp, pi, rx, ry, rz, norm, g0para, g0perp, minpara, maxpara, minperp, maxperp;
  double dtheta, angle, sp, r, delr, tref=0.0, DxP[3], Dx[3], DxNem[3], **g0, **gav, **g0Perp, **g0Parall, g0m, distSq, rlower, rupper, cost, nIdeal, cost2;
  double time, refTime, RCUT;
  int iZ, jZ, iX, jX, iY, jY, NP1, NP2;
  double shift[3];
  double ene=0.0;

  rv[0]=1;
  rv[1]=rv[2]=0;
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
	      if (mcsim==1)
	       {
	         for (i=0; i < NP; i++)
		   fscanf(f, "%[^\n]\n", line);
	         fscanf(f, "%lf\n", &L);
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
      u[a] = malloc(sizeof(double)*NP);
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
  delrp= L / 2.0 / ((double)perppoints);
  norm = calc_norm(rv);
  for (k=0; k < 3; k++)
    {
      rv[k] /= norm;
    }
  rewind(f2);
  nf = 0;
  if (angpoints==-1)
    angpoints = points;
 
  g0 = malloc(sizeof(double*)*points*4);
  if (calciso)
    gav = malloc(sizeof(double*)*points*4);
  for (k1=0; k1 < 4*points; k1++)
    {
      g0[k1] = malloc(sizeof(double)*angpoints);
      if (calciso)
	gav[k1] = malloc(sizeof(double)*angpoints);
	for (k2=0; k2 < angpoints; k2++)
	  {
	    g0[k1][k2] = 0.0;
	    if (calciso)
	      gav[k1][k1] = 0.0;
	  }
    }
  pi = acos(0.0)*2.0;
  dtheta = pi/angpoints; 
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      nf++;
      readconf(fname, &time, &refTime, NP, x, w, DR, u);
      for (i=0; i < NP-1; i++)
	for (j = i+1; j < NP; j++)
	  {
	    for (a = 0; a < 3; a++)
	      {
		Dx[a] = x[a][i] - x[a][j];
	      }

	    for (a = 0; a < 3; a++)
	      Dx[a] = Dx[a] - L * rint(Dx[a]/L);
	    if (calciso)
	      {
		distSq = 0.0;
		for (a = 0; a < 3; a++)
		  distSq += Sqr(Dx[a]);
		bin = (int) (sqrt(distSq)/delr);
	      }
	    lab2body(i, Dx, DxP);
	    projectVec(DxP, rv, &binx, &biny, &binz, delr, delrp);
	    //printf("P bin=%d %d %d\n", binx, biny, binz);
	    //binx = (int) floor(DxP[0] / delr);
	    //biny = (int) floor(DxP[1] / delrp);
	    //binz = (int) floor(DxP[2] / delrp);
	    //printf("D bin=%d %d %d\n", binx, biny, binz);
	    //printf("bin=%d %d %d\n", binx, biny, binz);
	    sp = 0.0;
	    for (a=0; a < 3; a++)
	      sp += u[a][i]*u[a][j];
	    angle = acos(sp);
	    bina = (int) (angle/dtheta); 
	    if (bina < angpoints)
	      {
		if ((biny==0 || biny==-1) && (binz==0 || binz==-1))
		  g0[binx+2*points][bina] += 2.0;
		if (calciso && bin >= 0 && bin < 4*points)
		  gav[bin][bina] += 2.0;
	      }
    	    //printf("bin=%d\n", bin);
	    //exit(1);
	      
	  }
    }
  fclose(f2); 
  //cost = 4.0 * pi * NP / 3.0 / (L*L*L);
  angle = dtheta*0.5;
  cost = (L*L*L)/((double)NP)/((double)NP)/(delrp*delrp*delr);
  cost /= 2;
  cost2 = 4.0 * pi * NP / 3.0 / (L*L*L);

  for (jj=0; jj < angpoints; jj++)
    {
      r = delr*0.5;
      sprintf(fname, "gr_theta%f.dat", angle);
      f = fopen(fname, "w+");
      for (ii = 0; ii < 4*points; ii++)
	{
	  g0m = cost*g0[ii][jj]/((double)nf)/sin(angle)/dtheta;
	  r = (((double)ii)-2*points)*delr;
	 
	  fprintf(f, "%.15G %.15G\n", r, g0m);
	  //r += delr;
    	}
      fclose(f);
      if (calciso)
	{
	  sprintf(fname, "griso_theta%f.dat", angle);
	  f = fopen(fname, "w+");
	  r = 0.5*delr;
	  for (ii = 0; ii < 4*points; ii++)
	    {
	      rlower = ( (double) ii) * delr;
	      rupper = rlower + delr;
	      nIdeal = cost2 * (Sqr(rupper)*rupper - Sqr(rlower)*rlower);
	      g0m = gav[ii][jj]/((double)nf)/((double)NP)/sin(angle)/dtheta/nIdeal;
	      //g0m = g0[ii][jj]/((double)nf)/((double)NP)/nIdeal;
#if 0
	      g2m = (3.0*g2[ii]/cc[ii] - 1.0)/2.0;
	      g4m = (35.0*g4[ii]/cc[ii] - 30.0*g2[ii]/cc[ii] + 3.0) / 8.0;
	      g6m = (231.0*g6[ii]/cc[ii] - 315.0*g4[ii]/cc[ii] + 105.0*g2[ii]/cc[ii] - 5.0)/16.0;
	      fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", r, g0m, g2m, g4m, g6m);
#endif
	      fprintf(f, "%.15G %.15G\n", r, g0m);
	      r += delr;
	    }
	}
      fclose(f);
      angle += dtheta;
    }
  return 0;
}

