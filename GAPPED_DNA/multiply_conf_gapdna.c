#include <stdlib.h>
#include <stdio.h>
#include <math.h>
char line[16384];
char par[16384], val[16384];
int factx, facty, factz;
int main(int argc, char **argv)
{
  int numat, parnum;
  double *rx, *ry, *rz, rCMx, rCMy, rCMz, *R[3][3], L, Lnew, Lxnew, Lynew, Lznew;
  double Lx, Ly, Lz, vx, vy, vz, wx, wy, wz, dx, dy, dz, rxn, ryn, rzn;
  FILE *f;
  long pos;
  int a, b, type, jx, jy, jz, i, k1, k2;
  int fact, fact3;
  f = fopen(argv[1],"r");
  while (!feof(f))
    fscanf(f, "%[^\n] ", line);
  sscanf(line, "%lf %lf %lf\n", &Lx, &Ly, &Lz);
  fclose(f);
  f = fopen(argv[1],"r");
  numat = 0;
  if (argc==3)
    {
      factx = facty = factz =atoi(argv[2]);

    }
  else if (argc==5)
    {
      factx = atoi(argv[2]);
      facty = atoi(argv[3]);
      factz = atoi(argv[4]); 
    }
  else 
    {
      printf("multiply_conf_mc <conf file> <fact_x> <fact_y> <fact_z>\n");
      exit-(1);
    }
  
  fact3 = factx*facty*factz;
	
  for (;numat < 2;)
    {
      fscanf(f,"%[^\n] ", line);
      sscanf(line, "%[^:]:%[^\n] ", par, val);
      if (!strcmp(par,"parnum"))
	{	
	  parnum = atoi(val);
	  printf("parnum: %d\n", atoi(val)*fact3);
	}
      else if (!strcmp(par,"@@@"))
	{
	  numat++;
	  printf("@@@\n");
	  if (numat==1)
	    {
	      fscanf(f, "%[^\n] ", line);
	      printf("%d\n", parnum*fact3);
	    }
	}
      else
	printf("%s\n", line);
    }
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  for (a=0; a < 3; a++)	
    for (b=0; b < 3; b++)	
 	R[a][b] = malloc(sizeof(double)*parnum);
  Lxnew = Lx*factx;
  Lynew = Ly*facty;
  Lznew = Lz*factz;
  pos = ftell(f);
   
  rCMx = rCMy= rCMz = 0;
  for (jx = 0; jx < factx; jx++)
    {
      for (jy = 0; jy < facty; jy++)
      	{
	  for (jz = 0; jz < factz; jz++)
	    {
	      rCMx += Lx*jx;
	      rCMy += Ly*jy;
	      rCMz += Lz*jz;
	    }
	}
    }
  rCMx /= fact3;
  rCMy /= fact3;
  rCMz /= fact3;
  /* read all coords and orient */
  for (i=0; i < parnum; i++)
    {
      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
	    &rx[i], &ry[i], &rz[i], &R[0][0][i], &R[0][1][i], &R[0][2][i], &R[1][0][i], &R[1][1][i], &R[1][2][i], 
	    &R[2][0][i], &R[2][1][i], &R[2][2][i], &type); 
    } 
  /* join molecules first */
  for (i=0; i < parnum; i+=2)
    {
       dx = rx[i+1] - rx[i];
       dy = ry[i+1] - ry[i];
       dz = rz[i+1] - rz[i];
       rx[i+1] -= Lx*rint(dx/Lx);
       ry[i+1] -= Ly*rint(dy/Ly);
       rz[i+1] -= Lz*rint(dz/Lz);	 
    }
  for (jx = 0; jx < factx; jx++)
    {
      for (jy = 0; jy < facty; jy++)
	{
	  for (jz = 0; jz < factz; jz++)
	    {
              for (i=0; i < parnum; i++)
                {
                 rxn = rx[i] + Lx*jx - rCMx;
	         ryn = ry[i] + Ly*jy - rCMy;
		 rzn = rz[i] + Lz*jz - rCMz;
  		 /* apply periodic boundary conditions in the new box */
		 rxn = rxn - Lxnew*rint(rxn/Lxnew); 
		 ryn = ryn - Lynew*rint(ryn/Lynew);
		 rzn = rzn - Lznew*rint(rzn/Lznew);
                 printf("%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n",
			 rxn, ryn, rzn,
			 R[0][0][i], R[0][1][i], R[0][2][i], R[1][0][i], R[1][1][i], R[1][2][i], 
			 R[2][0][i], R[2][1][i], R[2][2][i], type); 
		}
	    }
	}
    }


#if 0
  for (i=0; i < parnum; i++)
    {
      fscanf(f, "%lf %lf %lf %lf %lf %lf\n",
	    &vx, &vy, &wz, &wx, &wy, &wz); 
      for (jx = 0; jx < fact; jx++)
	{
	  for (jy = 0; jy < fact; jy++)
	    {
	      for (jz = 0; jz < fact; jz++)
		{
		  printf("%.15G %.15G %.15G %.15G %.15G %.15G\n",
			 vx, vy, vz, wx, wy, wz); 
		}
	    }
	}

    }
#endif 
  printf("%.15G %.15G %.15G\n", Lxnew, Lynew, Lznew);
  fclose(f);
}
