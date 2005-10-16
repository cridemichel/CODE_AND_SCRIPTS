#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>
char line[100000], parname[124], parval[1024];
int N;
double x[3], R[3][3], Q[3][3];
char fname[1024], inputfile[1024];
int readCnf = 0, timeEvol = 0, ordmatrix=0;
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
  jobz='N';
  uplo='U';
  lwork = 45;
  dsyev_(&jobz, &uplo, &lda, a, &lda, ev, work, &lwork,  &info);  
}
void print_usage(void)
{
  printf("order_param [--cnf/-c | --time/-t | --ordmatrix/-Q ] <confs_file>\n");
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
      else if (!strcmp(argv[cc],"--cnf") || !strcmp(argv[cc],"-c" ))
	{
	  readCnf = 1;
	} 
      else if (!strcmp(argv[cc],"--time") || !strcmp(argv[cc],"-t"))
	{
	  timeEvol = 1;
	}
      else if (!strcmp(argv[cc],"--ordmatrix") || !strcmp(argv[cc],"-Q"))
	{
	  ordmatrix = 1;
	}
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
}
int main(int argc, char** argv)
{
  FILE *f, *f2;
  int nf, i, a, b;
  double ev[3], ti, S, tref=0.0;
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
	      //sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
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
	  if (!strcmp(parname,"refTime"))
	    tref = atof(parval);
	}
      while (strcmp(line,"@@@"));

      //printf("fname=%s %d ellipsoids...\n", fname, N);
      if (timeEvol)
	{
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      Q[a][b] = 0.0;      
	}
      for (i=0; i < N; i++)
	{
	  fscanf (f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
		  &(x[0]), &(x[1]), &(x[2]), 
		  &(R[0][0]), &(R[0][1]), &(R[0][2]), &(R[1][0]), &(R[1][1]), 
		  &(R[1][2]), &(R[2][0]), &(R[2][1]), &(R[2][2]));

	  //printf("R=%.15G %.15G %.15G\n", R[0][1], R[0][1], R[0][2]);
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
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
#endif
}
