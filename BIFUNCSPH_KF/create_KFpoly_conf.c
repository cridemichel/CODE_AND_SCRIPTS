#include<stdlib.h>
#include<stdio.h>
double nx, ny, nz, L[3], *rx, *ry, *rz;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
int full, ibeg, numpoly;
int main(int argc, char **argv)
{
  FILE *f;
  double Diam, del0, distKF=0.546;
  int parnum=5848; //phi=0.299 i.e. coexistence according to Thuy's results using SUS
  int i, j, polylen=1;
  f = fopen(argv[1], "w+");
  Diam=1.0;
  /* elongated box */
  L[0] = 40.0*Diam;
  L[1] = 16.0*Diam;
  L[2] = 16.0*Diam;
  printf("argc=%d\n", argc);
  if (argc==3)
    parnum=atoi(argv[2]);
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
  fprintf(f,"ninters: 3\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"%f %f %f\n", Diam*0.5, Diam*0.5, Diam*0.5);
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n",Diam*0.5, Diam*distKF);
  fprintf(f,"0 0 %f %f\n",-Diam*0.5, Diam*distKF);
  fprintf(f,"0 0 0 0 1 0 0 100000\n");
  fprintf(f,"0 0 0 1 1 0 0 100000\n");
  fprintf(f,"0 1 0 1 1 0 0 100000\n");
  fprintf(f,"@@@\n");  
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  numpoly=parnum/polylen;
  printf("numpoly=%d numpoly*polylen=%d\n", numpoly, numpoly*polylen);
  drx = Diam*polylen*1.00001; //0.500000000001;
  dry = Diam*1.00001;
  drz = Diam*1.00001;
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx+del0;
      if (rxl+Diam*0.5 > L[0])
	{
	  ny++;
	  nx=0;
	  rxl = rxc[polylen-1]+nx*drx+del0;
	}
      ryl = ryc[polylen-1]+ny*dry+del0;
      if (ryl+Diam*0.5 > L[1])
	{
	  nz++;
	  ny=0;
	  ryl = ryc[polylen-1]+ny*dry+del0;
      	}
      rzl = rzc[polylen-1]+nz*drz+del0;
      if (rzl+Diam*0.5 > L[2])
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
      rx[i] -= L[0]*0.5;
      ry[i] -= L[1]*0.5;
      rz[i] -= L[2]*0.5;
      fprintf(f, "%f %f %f 1 0 0 0 1 0 0 0 1 0\n", rx[i], ry[i], rz[i]);
    }
  fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  return 1;
} 
