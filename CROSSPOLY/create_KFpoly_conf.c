#include<stdlib.h>
#include<stdio.h>
double nx, ny, nz, L, *rx, *ry, *rz;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
int full, ibeg, numpoly;
int main(int argc, char **argv)
{
  FILE *f;
  double del, del0, maxL;
  int parnum=2800, i, j, polylen=4;
  f = fopen(argv[1], "w+");
  L = 30.75;
  printf("argc=%d\n", argc);
  if (argc==3)
    parnum=atoi(argv[2]);
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxc =malloc(sizeof(double)*polylen);
  ryc =malloc(sizeof(double)*polylen);
  rzc =malloc(sizeof(double)*polylen);
  del=1.51;
  del0 = 0.01;
  for (i=0; i < polylen; i++)
    {
      rxc[i] = i;
      ryc[i] = 0.0;
      rzc[i] = 0.0;
    }
  fprintf(f, "parnum: %d\n", parnum);
  fprintf(f,"ninters: 4\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"0.5 0.5 0.5\n");
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"3 0\n");
  fprintf(f,"0.5 0 0 0.5\n");
  fprintf(f,"-0.5 0 0 0.5\n");
  fprintf(f,"0.0 0.5 0 0.25\n");
  fprintf(f,"0 0 0 0 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 0 0 1 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 1 0 1 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 2 0 2 1 0 0 100000\n"); /* kern-frenkel patches */
  fprintf(f,"@@@\n");  
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  maxL = L;
  numpoly=parnum/polylen;
  printf("numpoly=%d numpoly*polylen=%d\n", numpoly, numpoly*polylen);
  drx = 5.6;
  dry = 1.01;
  drz = 1.01;
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx+del0;
      if (rxl >= maxL)
	{
	  ny++;
	  nx=0;
	  rxl = rxc[polylen-1]+nx*drx+del0;
	}
      ryl = ryc[polylen-1]+ny*dry+del0;
      if (ryl >= L-1.0)
	{
	  nz++;
	  ny=0;
	  ryl = ryc[polylen-1]+ny*dry+del0;
      	}
      rzl = rzc[polylen-1]+nz*drz+del0;
      if (rzl >= L-1.0)
	{
	  full=1;
	  break;
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
