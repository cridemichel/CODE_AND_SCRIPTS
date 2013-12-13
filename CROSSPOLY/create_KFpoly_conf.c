#include<stdlib.h>
#include<stdio.h>
double L, *rx, *ry, *rz;

int main(int argc, char **argv)
{
  FILE *f;
  double nx, ny, nz, del, del0;
  int parnum=2800, i, polylen=4;
  f = fopen(argv[1], "w+");
  L = 30.75;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);

  fprintf(f, "parnum: %d\n", parnum);
  fprintf(f,"ninters: 4\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f,"maxbondsSaved: -1\n");
  fprintf(f,"@@@\n");
  fprintF(f, "%d\n", parnum);
  fprintf(f,"0.5 0.5 0.5\n");
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"3 0\n");
  fprintf(f,"0.5 0 0 0.5\n");
  fprintf(f,"-0.5 0 0 0.5\n");
  fprintf(f,"0 0 0 0 1 0 0 100000\n");
  fprintf(f,"0 0 0 1 1 0 0 100000\n");
  fprintf(f,"0 1 0 1 1 0 0 100000\n");
  fprintf(f,"0 2 0 2 1 0 0 100000\n"); /* kern-frenkel patches */
  fprintf(f,"@@@\n");  
  nx=ny=nz=0;
  full=0;
  del=0.0;
  del0 = 0.0001;
  ibeg = 0;
  for (i=0; i < parnum; i++)
    {
      if (parnum % polylen == 0)
	{
	  del = 1.51;
	  ibeg = i+1;
	}
      else
	{ 
	  del = 0.0;
	}
      rx[i] = nx+del0+del;
      nx++;
      if (rx[i] > L*0.5)
	{ 
	  if (del==0.0)
	    {
	      for (j=ibeg; j < i; j++)
		{
		  ry[j] += 1.0+del0;
 		}
	    }
	  nx--;
	  rx[i] = nx+del0+del;

	  ny++;
	  ry[i] = ny+del0;
	  if (ry[i] > L*0.5)
	    {
	      ny--;
	      ry[i] = ny+del0;
	      nz++; 
	      if (del==0.0)
		{
		  for (j=ibeg; j < i; j++)
		    {
		      rz[j] += 1.0+del0;
		    }
		}
	      rz[i] = nz+del0;
	      if (rz[i] > L*0.5)
		full=1;
	    }
	}
      if (full==1)
	{
	  printf("I could only place %d particles!\n", i);
	  break;
	}
      rx -= L*0.5;
      ry -= L*0.5;
      rz -= L*0.5;
      fprintf("%f %f %f 1 0 0 0 1 0 0 0 1 0\n", rx, ry, rz);
    }
  fprintf(f, "%.15G %.15G %15G\n", L, L, L);
  return 1;
} 
