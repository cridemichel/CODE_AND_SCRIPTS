#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
char line[100000], parname[124], parval[256000];
int N, mcsim=0;
double x[3], R[3][3], Q[3][3];
double r1[3], r2[3], r3[3], u1[3], u2[3], u3[3], dt;
char fname[1024], inputfile[1024];
int readCnf = 0, timeEvol = 0, ordmatrix=0, eigenvectors=0, curStep, readLapo=0, heflex=0;
#define Sqr(x) ((x)*(x))
void vectProd(double r1x, double r1y, double r1z, 
	 double r2x, double r2y, double r2z, 
	 double* r3x, double* r3y, double* r3z)
{
  /* DESCRIPTIOM:
     r3 = [ r1, r2 ] where [ , ] the vectorial product */
  *r3x = r1y * r2z - r1z * r2y; 
  *r3y = r1z * r2x - r1x * r2z;
  *r3z = r1x * r2y - r1y * r2x;
}
double eigvec[3][3];
void diagonalize(double M[3][3], double ev[3])
{
  double a[9], work[45];
  char jobz, uplo;
  int info, i, j, lda, lwork;
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) a[j+3*i]=M[j][i];		
      //for(j=0; j<3; j++) a[j][i]=M[j][i];		
    }	
  lda = 3;
  if (eigenvectors)
    jobz='V';
  else
    jobz='N';
  uplo='U';
  lwork = 45;
  dsyev_(&jobz, &uplo, &lda, a, &lda, ev, work, &lwork,  &info);  
  if (!eigenvectors)
    return;
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) eigvec[i][j]=a[j+3*i];		
    }	
}
void print_usage(void)
{
  printf("order_param [ -mc | --lapo/-l | --heflex/-fl |--cnf/-c | --time/-t | --ordmatrix/-Q ] [--eigenvec/-ev] <confs_file>\n");
  exit(0);
}
void parse_param(int argc, char** argv)
{
  int cc=1;
  
  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"-mc"))
	{
	  mcsim = 1;
	} 
      else if (!strcmp(argv[cc],"--cnf") || !strcmp(argv[cc],"-c" ))
	{
	  readCnf = 1;
	} 
      else if (!strcmp(argv[cc],"--heflex") || !strcmp(argv[cc],"-fl" ))
	{
	  heflex = 1;
	} 
      else if (!strcmp(argv[cc],"--lapo") || !strcmp(argv[cc],"-l" ))
	{
	  readLapo = 1;
	} 
      else if (!strcmp(argv[cc],"--time") || !strcmp(argv[cc],"-t"))
	{
	  timeEvol = 1;
	}
      else if (!strcmp(argv[cc],"--ordmatrix") || !strcmp(argv[cc],"-Q"))
	{
	  ordmatrix = 1;
	}
      else if (!strcmp(argv[cc],"--eigenvec") || !strcmp(argv[cc],"-ev"))
	{
	  eigenvectors = 1;
	}
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
}
double calc_norm(double *v)
{
  int a;
  double norm=0;
  for (a = 0; a < 3; a++)
   norm += Sqr(v[a]);
  return sqrt(norm); 
}
void build_ref_axes2(double u1[3], double u2[3], double u3[3])
{
  int a;
  double norm1, u1u2, norm2;
  u1[0] = r2[0] - r1[0];
  u1[1] = r2[1] - r1[1];
  u1[2] = r2[2] - r1[2];
  u2[0] = r3[0] - r1[0];
  u2[1] = r3[1] - r1[1];
  u2[2] = r3[2] - r1[2];
  norm2 = calc_norm(u2);
  for (a=0;  a < 3; a++)
    u2[a] /= norm2;
  u1u2 = u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
  for (a=0;  a < 3; a++)
    {
      u1[a] = u1[a] - u1u2*u2[a];
    }
  norm1 = calc_norm(u1);
  for (a=0;  a < 3; a++)
    {
      u1[a] /= norm1;
    }
  vectProd(u1[0],u1[1],u1[2],u2[0],u2[1], u2[2],&u3[0], &u3[1], &u3[2]); 
}

void build_ref_axes(double u1[3], double u2[3], double u3[3])
{
  int a;
  double n;
  u1[0] = r2[0] - r1[0];
  u1[1] = r2[1] - r1[1];
  u1[2] = r2[2] - r1[2];
  n = calc_norm(u1);
  for (a = 0; a < 3; a++)
    u1[a] /= n;
    
  u2[0] = 2.0*r3[0]-(r1[0] + r2[0]);
  u2[1] = 2.0*r3[1]-(r1[1] + r2[1]);
  u2[2] = 2.0*r3[2]-(r1[2] + r2[2]);
  
  n = calc_norm(u2);
  for (a = 0; a < 3; a++)
    u2[a] /= n;
   
  //fprintf(stderr,"u1.u2=%.15G\n", u1[0]*u2[0]+u1[1]*u2[1]+u1[2]*u2[2]);
  vectProd(u1[0], u1[1], u1[2], u2[0], u2[1], u2[2], &u3[0], &u3[1], &u3[2]);  

}

int main(int argc, char** argv)
{
  FILE *f, *f2;
  int nf, i, a, b, dummyint;
  double ev[3], ti=-1, S, tref=0.0;
#if 0
  if (argc == 1)
    {
      printf("file with all confs to read as input\n");
      exit(-1);
    }
#endif
  parse_param(argc, argv);
  f2 = fopen(inputfile,"r");
  nf = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      f = fopen(fname,"r");
      nf++;
      tref=0.0;
      if (readCnf)
	{
	  do 
	    {
	      fscanf(f,"%[^\n]\n",line);
	      sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	      if (!strcmp(parname,"refTime"))
		tref = atof(parval);
	    }
	  while (strcmp(line,"@@@"));
	}
      do 
	{
	  fscanf(f,"%[^\n]\n",line);
	  sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	  if (!strcmp(parname,"parnum"))
	    N = atoi(parval);
	  if (!strcmp(parname,"time"))
	    ti = atof(parval);
	  if ((readLapo||mcsim) && !strcmp(parname,"curStep"))
	    curStep = atoi(parval);
	  if (readLapo && !strcmp(parname,"steplength"))
	    dt = atof(parval);
  
  	}
      while (strcmp(line,"@@@"));
      if (heflex)
	{
	  do 
	    {
	      fscanf(f,"%[^\n]\n",line);
	    }
	  while (strcmp(line,"@@@"));
	}
      //printf("fname=%s %d ellipsoids...\n", fname, N);
      if (timeEvol)
	{
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      Q[a][b] = 0.0;      
	}
      for (i=0; i < N; i++)
	{
	  if (readLapo)
	    {
	      fscanf (f, "%lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n",
	  	    &(r1[0]), &(r1[1]), &(r1[2]), 
	  	    &(r2[0]), &(r2[1]), &(r2[2]), 
		    &(r3[0]), &(r3[1]), &(r3[2]));
	      build_ref_axes(u1, u2, u3);
	      for (a=0; a < 3; a++)
		{
		  R[0][a] = u1[a];
		  R[1][a] = u2[a];
		  R[2][a] = u3[a];
		}
	    }
	  else if (heflex)
	    fscanf (f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
	  	    &(x[0]), &(x[1]), &(x[2]), 
	  	    &(R[0][0]), &(R[0][1]), &(R[0][2]), &(R[1][0]), &(R[1][1]), 
	  	    &(R[1][2]), &(R[2][0]), &(R[2][1]), &(R[2][2]), &dummyint);
	  else
	    fscanf (f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
	  	    &(x[0]), &(x[1]), &(x[2]), 
	  	    &(R[0][0]), &(R[0][1]), &(R[0][2]), &(R[1][0]), &(R[1][1]), 
	  	    &(R[1][2]), &(R[2][0]), &(R[2][1]), &(R[2][2]));

	  //printf("R=%.15G %.15G %.15G\n", R[0][1], R[0][1], R[0][2]);
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
		if (readLapo)
		  Q[a][b] += 1.5 * R[2][a]*R[2][b];
		else
		  Q[a][b] += 1.5 * R[0][a]*R[0][b];
		if (a==b)
		  Q[a][a] -= 0.5;
	      }
	}
      if (timeEvol)
	{
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      Q[a][b] /= ((double)N); 
	  diagonalize(Q, ev);
	  if (mcsim)
	    printf("%d %.15G %.15G %.15G\n", curStep, ev[0], ev[1], ev[2]); 
	  else if (readLapo)
	    printf("%.15G %.15G %.15G %.15G\n", curStep*dt, ev[0], ev[1], ev[2]); 
	  else
	    printf("%.15G %.15G %.15G %.15G\n", ti+tref, ev[0], ev[1], ev[2]); 
	}

      fclose(f);
    }
  fclose(f2); 
  if (timeEvol)
    return 0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      Q[a][b] /= ((double)N)*((double)nf); 
  diagonalize(Q, ev);
  if (fabs(ev[0]) > fabs(ev[1]))
    S = ev[0];
  else
    S = ev[1];  
  if (fabs(ev[2]) > S)
    S = ev[2];
  printf("Order Parameter: %.8G (of %.8G %.8G %.8G)\n", S, 
	 ev[0], ev[1], ev[2]);
#if 1
  if (!timeEvol && ordmatrix)
    {
      printf("Ordering matrix Q:\n");
      printf("{");
      for (a=0; a < 3; a++)
	{
	  printf("{");
	  for (b=0; b < 3; b++)
	    {
	      if (b < 2)
		printf("%.15f,", Q[a][b]);	  
	      else
		printf("%.15f",Q[a][b]);
	    }
	  if (a < 2)
	    printf("},\n");
	  else if (b < 2)
	    printf("}\n");
	  else 
	    printf("}");
	}
      printf("}\n");
    }
  if (!timeEvol && eigenvectors)
    {
      printf("Eigenvectors matrix:\n");
      for (a=0; a < 3; a++)
	{
	  printf("lambda=%.15G {%.15f,%.15G,%.15G}\n",ev[a], eigvec[a][0], eigvec[a][1], eigvec[a][2]);	  
       	}
    }
#endif
}
