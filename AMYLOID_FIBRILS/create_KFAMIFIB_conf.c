#include<stdlib.h>
#include<stdio.h>
#include<math.h>
double nx, ny, nz, L, *rx, *ry, *rz;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
int full, ibeg, numpoly;
int main(int argc, char **argv)
{
  FILE *f;
  double phi, Diam, del0, maxL, pi;
  int numpoly, parnum=2800, i, j, polylen=20;

  if (argc == 1)
   {

     printf("create_KFpoly_conf <conf_file_name> <num. polymers> <volume fraction> <polymer length>\n"); 
     exit(-1);
   }
  f = fopen(argv[1], "w+");
  Diam=2.0;
  pi = acos(0.0)*2.0;

  if (argc > 4) 
    polylen = atoi(argv[4]);

  numpoly=500;
  if (argc > 2)
    {
      numpoly = atoi(argv[2]);
      parnum=numpoly*polylen;
    }
  else
    parnum =500*polylen; 
  if (argc > 3)
    {
      phi = atof(argv[3]);
      L= pow(((double)2.0),2.0/3.0)*pow(((double)numpoly),1.0/3.0)*pow(pi/3.0,1.0/3.0)*pow(((double)polylen),1.0/3.0)/pow(phi,1.0/3.0);
    }
  else
    L=50.0;
  printf("L=%f phi=%f argc=%d polylen=%d\n", L, phi, argc, polylen);
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxc =malloc(sizeof(double)*polylen);
  ryc =malloc(sizeof(double)*polylen);
  rzc =malloc(sizeof(double)*polylen);
  del0 = Diam*0.5+0.00001;
  for (i=0; i < polylen; i++)
    {
      rxc[i] = i*Diam;
      ryc[i] = 0.0;
      rzc[i] = 0.0;
    }
  fprintf(f, "parnum: %d\n", parnum);
  fprintf(f,"ninters: 11\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"1 1 1\n");
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"6 0\n");
  fprintf(f,"1 0 0 0.119\n");/* 0: along x axis (permanent) */
  fprintf(f,"-1 0 0 0.119\n");/* 1: along x axis (permanent) */
  fprintf(f,"0  1  0 0.5\n");/* 2: along y axis */ 
  fprintf(f,"0 -1  0 0.5\n");/* 3: along y axis */
  fprintf(f,"0  0  1 0.5\n");/* 4: along z axis */
  fprintf(f,"0  0 -1 0.5\n");/* 5: along z axis */
  fprintf(f,"0 0 0 0 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 0 0 1 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 1 0 1 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 2 0 2 0.323 0 0 100000\n"); /* y-axis kern-frenkel patches */
  fprintf(f,"0 2 0 3 0.323 0 0 100000\n"); /* y-axis kern-frenkel patches */
  fprintf(f,"0 3 0 3 0.323 0 0 100000\n"); /* y-axis kern-frenkel patches */
  fprintf(f,"0 4 0 4 1.67 0 0 100000\n"); /* z-axiz KF patch */ 
  fprintf(f,"0 5 0 5 1.67 0 0 100000\n");  /* z-axiz KF patch */
  fprintf(f,"0 4 0 5 4.65 0 0 100000\n");  /* z-axiz KF patch */
  fprintf(f,"0 3 0 4 0.8 0 0 100000\n");  /* perp conf */
  fprintf(f,"0 2 0 5 0.8 0 0 100000\n");  /* perp conf */
  fprintf(f,"@@@\n");  
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  maxL = L;
  numpoly=parnum/polylen;
  printf("numpoly=%d numpoly*polylen=%d\n", numpoly, numpoly*polylen);
  drx = Diam*((double)polylen); //0.500000000001;
  dry = Diam;
  drz = Diam;
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx+del0;
      if (rxl+Diam*0.5 > L)
	{
	  ny++;
	  nx=0;
	  rxl = rxc[polylen-1]+nx*drx+del0;
	}
      ryl = ryc[polylen-1]+ny*dry+del0;
      if (ryl+Diam*0.5 > L)
	{
	  nz++;
	  ny=0;
	  ryl = ryc[polylen-1]+ny*dry+del0;
      	}
      rzl = rzc[polylen-1]+nz*drz+del0;
      if (rzl+Diam*0.5 > L)
	{
	  full=1;
	}

      for (j = 0; j < polylen; j++)
	{
	  rx[i*polylen+j] = rxc[j]+nx*drx+del0;
	  //printf("nx=%f i=%d j=%d x=%f\n", nx, i, j,  rx[i*polylen+j]);
	  ry[i*polylen+j] = ryc[j]+ny*dry+del0;
	  rz[i*polylen+j] = rzc[j]+nz*drz+del0;
	  //printf("np=%d\n", i*polylen+j);
	}
      if (full==1)
	{
	  printf("I could only place %d particles!\n", i);
	  break;
	}
      nx++;
    }
  //  if (full)
  //   exit(-1);
  for (i=0; i < parnum; i++)
    {
      rx[i] -= L*0.5;
      ry[i] -= L*0.5;
      rz[i] -= L*0.5;
      fprintf(f, "%f %f %f 1 0 0 0 1 0 0 0 1 0\n", rx[i], ry[i], rz[i]);
    }
  fprintf(f, "%.15G %.15G %.15G\n", L, L, L);
  return 1;
} 
