#include<stdlib.h>
#include<stdio.h>
double L, *rx, *ry, *rz;
double *rxc, *ryc, *rzc, rxl , drx, dry, drz, ryl, rzl;
int full, ibeg, numpoly;
int main(int argc, char **argv)
{
  FILE *f;
  double nx, ny, nz, del, del0;
  int parnum=2800, i, j, polylen=4;
  f = fopen(argv[1], "w+");
  L = 30.75;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxc =malloc(sizeof(double)*polylen);
  ryc =malloc(sizeof(double)*polylen);
  rzc =malloc(sizeof(double)*polylen);
  del=1.51;
  del0 = 0.0001;
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
  fprintf(f,"maxbondsSaved: -1\n");
  fprintf(f,"@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"0.5 0.5 0.5\n");
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"3 0\n");
  fprintf(f,"0.5 0 0 0.5\n");
  fprintf(f,"-0.5 0 0 0.5\n");
  fprintf(f,"0.0 0.5 0 0.25\n");
  fprintf(f,"0 0 0 0 1 0 0 100000\n");
  fprintf(f,"0 0 0 1 1 0 0 100000\n");
  fprintf(f,"0 1 0 1 1 0 0 100000\n");
  fprintf(f,"0 2 0 2 1 0 0 100000\n"); /* kern-frenkel patches */
  fprintf(f,"@@@\n");  
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  numpoly=parnum/polylen;
  drx = 4.6;
  dry = 1.01;
  drz = 1.01;
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx+del;
      nx++;
      if (rxl > L*0.5)
	{
	  ny++;
	  nx=0;
	}
      ryl = ryc[polylen-1]+ny*dry;
      if (ryl > L*0.5)
	{
	  nz++;
	  ny=0;
      	}
      rzl = rxc[polylen-1]+nz*drz;
      if (rzl  > L*0.5)
	{
	  full=1;
	}

      for (j = 0; j < polylen; j++)
	{
	  rx[i*polylen+j] =  rxc[j]+nx*drx+del;
	  ry[i*polylen+j] =  ryc[j]+ny*dry;
	  rz[i*polylen+j] =  rzc[j]+nz*drz;
	}
      if (full==1)
	{
	  printf("I could only place %d particles!\n", i);
	  break;
	}
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
