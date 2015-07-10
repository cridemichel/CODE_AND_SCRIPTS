#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
int allocated_cyls=100;
int numcyls=0, ncyltot, ibin;
int npts = 1000;
double D, L, Lbox, r[3], u[3], delvol;
#define Sqr(x) ((x)*(x))
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
struct cylstr 
{
  double u[3];
  double r[3];
  double L;
  double D;
  int i1;
  int i2;
} *cylinders;

void add_cylinder(double r[], double u[], double L, double D, int i, int j)
{
  int k;
  numcyls++; 
  if (numcyls >= allocated_cyls)
    {
      allocated_cyls += 100;
      cylinders = realloc(cylinders, sizeof(struct cylstr)*allocated_cyls);
    }
  for (k=0; k < 3; k++)
    {
      cylinders[numcyls-1].r[k] = r[k];
      cylinders[numcyls-1].u[k] = u[k];  
    }
  cylinders[numcyls-1].i1 = i;
  cylinders[numcyls-1].i2 = j;
  cylinders[numcyls-1].D = D;
  cylinders[numcyls-1].L = L;
  //printf("allocated_cyls=%d numcyls-1=%d D=%f\n", allocated_cyls, numcyls-1, cylinders[numcyls-1].D);  
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

int point_is_inside(double p[3], int icyl)
{
  double norm, sp, distsq, C[2][3], delp[3];
  int kk, k1, k2;
  
  for (kk = 0; kk < 3; kk++)
    {
      C[0][kk] = cylinders[icyl].r[kk] + cylinders[icyl].u[kk]*cylinders[icyl].L*0.5;
      C[1][kk] = cylinders[icyl].r[kk] - cylinders[icyl].u[kk]*cylinders[icyl].L*0.5;
    }
  /* punto dentro le sfere */
  for (k1=0; k1 < 2; k1++)
    {
      distsq = Sqr(p[0]-C[k1][0])+Sqr(p[1]-C[k1][1])+Sqr(p[2]-C[k1][2]); 
      if (distsq < Sqr(cylinders[icyl].D*0.5))
	return 1;
    }

  /* punto dentro al cilindro */	
  for (kk = 0; kk < 3; kk++)
    delp[kk] = p[kk] - cylinders[icyl].r[kk];

  sp = scalProd(delp, cylinders[icyl].u);

  for (kk = 0; kk < 3; kk++)
    delp[kk] -= sp*cylinders[icyl].u[kk];

  norm = calc_norm(delp);

  if (norm  < cylinders[icyl].D*0.5 && fabs(sp) < cylinders[icyl].L*0.5)
    return 1;

  return 0;
}

void main(int argc, char **argv)
{
  FILE *cf, *f;
  double pp[3], vol, maxvol=0.0, PVtot;
  int ncyl, nc, i, j;
  long long int max_MC_trials=100000;
  double ov, *PV;
  long long int tt;

  if (argc==1)
    {
      printf("calc_cls_vol <file> <MC trials> [histogram points]\n");
      exit(-1);
    }
  cf = fopen(argv[1], "r");
  max_MC_trials = atoll(argv[2]);
  if (argc==4)
    npts = atoi(argv[3]);
  PV = malloc(sizeof(double)*npts);
  srand((int)time(NULL));
  cylinders = malloc(sizeof(struct cylstr)*allocated_cyls);
  maxvol = 0.0;
  while(!feof(cf))
    {
      fscanf(cf, "%d  %lf\n", &ncyl, &Lbox);
      ncyltot += ncyl;
      vol = 0.0;
      for (nc = 0; nc < ncyl; nc++)
	{
	  fscanf(cf, "%lf %lf %lf %lf %lf %lf %lf %lf ", &r[0], &r[1], &r[2],
	    	 &u[0], &u[1], &u[2], &L, &D);
	  vol += L*(D/2.)*(D/2.)*M_PI;
	}
      if (vol > maxvol)
	maxvol = vol;
    }
  delvol = maxvol / npts; 
  rewind(cf);
  f = fopen("all_clusters_volumes.txt", "w+");
  while (!feof(cf))
    {
      ov = 0.0;
      fscanf(cf, "%d  %lf\n", &ncyl, &Lbox);
      for (nc = 0; nc < ncyl; nc++)
	{
	  fscanf(cf, "%lf %lf %lf %lf %lf %lf %lf %lf ", &r[0], &r[1], &r[2],
		      &u[0], &u[1], &u[2], &L, &D);
	  //printf("r=%f %f %f u=%f %f %f L=%f D=%f\n", r[0], r[1], r[2], u[0], u[1], u[2], L, D);
	  add_cylinder(r, u, L, D, nc, nc);
	}	
      for (tt=0; tt < max_MC_trials; tt++)
	{
	  pp[0] = Lbox*(drand48()-0.5); 
	  pp[1] = Lbox*(drand48()-0.5); 
	  pp[2] = Lbox*(drand48()-0.5);

	  for (j = 0; j < ncyl; j++)
	    {
	      if (point_is_inside(pp, j))
		{
		  //printf("qui\n");
		  ov=ov+1.0;
		}
	    }
	}
      printf("maxvol=%f delvol=%.15G\n", maxvol, delvol);
      vol = Lbox*Lbox*Lbox*ov/((double)tt); 
      ibin = (int) (vol/delvol);
      PV[ibin]++;
      fprintf(f, "%d %.15G\n", nc, vol);
    }

  fclose(f);
  fclose(cf);
  f = fopen("avg_vol_distr.dat", "w+");

  PVtot = 0.0;
  for (i=0; i < npts; i++)
    PVtot += PV[i];

  for (i=0; i < npts; i++)
    {
      fprintf(f, "%.15G %.15G\n", i*delvol, PV[i]/PVtot);
    }
  fclose(f);

}
