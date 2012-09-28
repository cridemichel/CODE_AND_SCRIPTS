#include <stdlib.h>
#include <stdio.h>
char line[16384];
char par[16384], val[16384];
int factx, facty, factz;
int main(int argc, char **argv)
{
  int numat, parnum;
  double rx, ry, rz, rCMx, rCMy, rCMz, R[3][3], L, Lnew, Lxnew, Lynew, Lznew;
  double Lx, Ly, Lz, vx, vy, vz, wx, wy, wz;
  FILE *f;
  long pos;
  int type, jx, jy, jz, i, k1, k2;
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
  for (i=0; i < parnum; i++)
    {
      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
	    &rx, &ry, &rz, &R[0][0], &R[0][1], &R[0][2], &R[1][0], &R[1][1], &R[1][2], 
	    &R[2][0], &R[2][1], &R[2][2], &type); 
      for (jx = 0; jx < factx; jx++)
	{
	  for (jy = 0; jy < facty; jy++)
	    {
	      for (jz = 0; jz < factz; jz++)
		{
		  printf("%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n",
			 rx+Lx*jx-rCMx, ry+Ly*jy-rCMy, rz+Lz*jz-rCMz,
			 R[0][0], R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], 
			 R[2][0], R[2][1], R[2][2], type); 
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
