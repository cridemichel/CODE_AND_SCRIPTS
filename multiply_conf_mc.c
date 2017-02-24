#include <stdlib.h>
#include <stdio.h>
#include <string.h>
char line[16384];
char par[16384], val[16384];
int factx, facty, factz;
double ex_x=1.0, ex_y=1.0, ex_z=1.0;
char inputfile[1024];
void print_usage(void)
{
  printf("multiply_conf_mc [--expand|-ex <val>][--exp-x|-ex <val>] [--exp-y|-ey <val>] [--exp-z|-ez <val>] <conf file> <fact_x> <fact_y> <fact_z>\n");
  exit(0);
}
int parse_param(int argc, char** argv)
{
  int extraparam=0;  
  int cc=1;
  //printf("argc==%d\n", argc);	
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
      else if (!strcmp(argv[cc],"--expand")||!strcmp(argv[cc],"-ea"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  ex_x = ex_y = ex_z = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--expandx")||!strcmp(argv[cc],"-ex"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  ex_x = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--expandx")||!strcmp(argv[cc],"-ey"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  ex_y = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--expandz")||!strcmp(argv[cc],"-ez"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  ex_z = atof(argv[cc]);
	}
      else if (cc==argc|| extraparam==4)
	print_usage();
      else if (extraparam == 0)
	{
	  extraparam++;
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam==1)
	{
	  extraparam++;
	  factx=facty=factz=atof(argv[cc]);
	}
      else if (extraparam==2)
	{
	  extraparam++;
	  facty=atof(argv[cc]);
	}
      else if (extraparam==3)
	{
	  extraparam++;
	  factz=atof(argv[cc]);
	}
      else 
	print_usage();
      cc++;
    }
  return 0;
}

int main(int argc, char **argv)
{
  int numat, parnum;
  double rx, ry, rz, rCMx, rCMy, rCMz, R[3][3], L, Lnew, Lxnew, Lynew, Lznew;
  double Lx, Ly, Lz, vx, vy, vz, wx, wy, wz;
  FILE *f;
  long pos;
  int type, jx, jy, jz, i, k1, k2;
  int fact, fact3;
  parse_param(argc, argv);
  f = fopen(inputfile,"r");
  while (!feof(f))
    fscanf(f, "%[^\n] ", line);
  sscanf(line, "%lf %lf %lf\n", &Lx, &Ly, &Lz);
  fclose(f);
  f = fopen(inputfile,"r");
  numat = 0;
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

  if (ex_x > 1.0)
    Lx*= ex_x; 
  if (ex_y > 1.0) 
    Ly*= ex_y;
  if (ex_z > 1.0)
    Lz*= ex_z;
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
		  if (ex_x > 1.0 || ex_y > 1.0 || ex_z > 1.0)
		    printf("%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n",
			   ex_x*rx+Lx*jx-rCMx, ex_y*ry+Ly*jy-rCMy, ex_z*rz+Lz*jz-rCMz,
			   R[0][0], R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], 
			   R[2][0], R[2][1], R[2][2], type); 
		  else
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
