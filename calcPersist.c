#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 1000
char **fname; 
double time, *ccA, *ccB, *GA, *GB, Gm,
       GmA, GmB, costh2, *ti, *u0[3], *ut[3], L, refTime, *rr[3];
int points=-1, assez, NP, NPA=-1, type=0;
char inputfile[2048], parname[128], parval[256000], line[256000];
char dummy[2048];
double A0=-1, A1=-1, B0=-1, B1=-1, C0=-1, C1=-1;
inline void normalize(double **u, int i)
{
  int a;
  double n=0.0;
  for (a=0; a < 3; a++)
    {
      n+=u[a][i]*u[a][i]; 
    } 
  n=sqrt(n);
  for (a=0; a < 3; a++)
    {
      u[a][i]/=n; 
    } 
}
void readconf(char *fname, double *ti, double *refTime, int NP, double *u[3])
{
  FILE *f;
  double r0, r1, r2, R[3][3];
  int nat=0, i, cpos, a, dummyint;
  f = fopen(fname, "r");
  while (!feof(f) && nat < 3) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	}
      if (nat < 2)
	{
	  fseek(f, cpos, SEEK_SET);
	  fscanf(f, "%[^:]:", parname);
	  //printf("[%s] parname=%s\n", fname, parname);
	  if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
  	  else if (!strcmp(parname, "refTime"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *refTime = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else if (nat == 3)
	{
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
		     &r0, &r1, &r2, &R[0][0], &R[0][1], &R[0][2],
		     &R[1][0], &R[1][1], &R[1][2], &R[2][0], &R[2][1], &R[2][2], &dummyint); 

	      //printf("%f %f %f\n",r0, r1, r2);
	      for (a = 0; a < 3; a++)
		{
		  u[a][i] = R[assez][a];
		}
	      rr[0][i] = r0;
	      rr[1][i] = r1;
	      rr[2][i] = r2; 
	      normalize(u, i);
	    }
  	  break; 
	}

    }
  fclose(f);
}
void print_usage(void)
{
  printf("Usage: calcPersist [--type/-t <0|1>] <lista_file> [points] \n");
  printf("where points is the number of points of the correlation function\n");
  exit(-1);
}
void parse_params(int argc, char** argv)
{
  int extraparam=0;  
  int cc=1;
  if (argc==1)
    {
      print_usage();
      exit(1);
    }
  while (cc < argc)
    {

      printf("argc=%d argv[argc]=%s cc=%d\n", argc, argv[cc], cc);
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"--type")||!strcmp(argv[cc],"-t"))
	{
	  cc++;
	  if (cc==argc)
	    print_usage();
	  type=atoi(argv[cc]);
	}
      else if (cc==argc|| extraparam==2)
	print_usage();
      else if (extraparam == 0)
	{
	  extraparam++;
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam==1)
	{
	  extraparam++;
	  points=atoi(argv[cc]);
	}
      else 
	print_usage();
      cc++;
    }
}
#define Sqr(x) ((x)*(x))
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
int main(int argc, char **argv)
{
  FILE *f, *f2, *fA, *fB;
  int del, first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a, j;
  int NN, fine, JJ, maxl, nfiles, nat, np, maxnp;
  double Ri[3], Rj[3], normi, normj;
  parse_params(argc, argv); 
  c2 = 0;
  f2 = fopen(inputfile, "r");
  maxl = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", dummy); 
      if (strlen(dummy)+1 > maxl)
	maxl = strlen(dummy)+1;
      c2++;
    }	
  nfiles = c2;
  rewind(f2);
  fname = malloc(sizeof(char*)*nfiles);
  for (ii=0; ii < nfiles; ii++)
    {
      fname[ii] = malloc(sizeof(char)*maxl);
      fscanf(f2, "%[^\n]\n", fname[ii]); 
    }

  f = fopen(fname[0], "r");
  nat = 0;
  while (!feof(f) && nat < 3) 
    {
      fscanf(f, "%[^\n]\n)", line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  if (nat==3)
	    {
	      for (i=0; i < 2*NP; i++)
		{
		  fscanf(f, "%[^\n]\n", line);
		  //printf("%s\n",line);
		}
	      fscanf(f, "%lf\n", &L);
	      break;
	    }
	  continue;
	}
      sscanf(line, "%[^:]:%[^\n]\n", parname, parval); 
      if (!strcmp(parname,"parnum"))
	NP = atoi(parval);
      else if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      else if (!strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (!strcmp(parname, "a"))
       	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &A0, &A1);
	}
      else if (!strcmp(parname, "b"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &B0, &B1);
	}
      else if (!strcmp(parname, "c"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &C0, &C1);
	}
    }
  fclose(f);
  if (points==-1)
    points = NP;
  if (A0==-1)
    {
      assez=0;
    }
  else
    {
      if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
	assez = 0;
      else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
	assez = 1;
      else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
	assez = 2;
    }
  assez = 0; /* we take x along which two spots lie*/
  if (NPA == -1)
    NPA = NP;
  fprintf(stderr, "allocating %d items NN=%d NP=%d num files=%d maxnp=%d assez=%d maxl=%d\n", points, NN, NP, nfiles, maxnp, assez, maxl);
  if (NPA < NP)
    {
      printf("NPA=%d\n", NPA);
    }
  for (a=0; a < 3; a++)
    {
      u0[a] = malloc(sizeof(double)*NP);
      ut[a] = malloc(sizeof(double)*NP);
      rr[a] = malloc(sizeof(double)*NP);
    }

  ccA = malloc(sizeof(double)*points);
  GA = malloc(sizeof(double)*points);
  if (NPA < NP)
    {
      GB = malloc(sizeof(double)*points);
      ccB = malloc(sizeof(double)*points);
    }
  ti = malloc(sizeof(double)*points);
  for (ii=0; ii < points; ii++)
    ti[ii] = -1.0;

  ti[0]=0.0;
  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    {
      GA[ii] = 0.0;
      ccA[ii] = 0.0;
      if (NPA < NP)
	{
	  GB[ii] = 0.0;
	  ccB[ii] = 0;
	}
    }
  c2 = 0;
  JJ = 0;
  if (type==1)
    del = 1;
  else
    del = 0;
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+1)
    {	
      readconf(fname[nr1], &time, &refTime, NP, u0);
      for (i=0; i < NP-del-1; i++) 
	{
	  for (j=i+1; j < NP-del; j++) 
	    {
	      costh2 = 0.0;
	      np = abs(j-i);
	      if (np > points) 
		continue;
	      ti[np] = np;
	      if (type==1)
		{
		  for (a=0; a < 3; a++)
		    {
		      Ri[a] = rr[a][i+1]-rr[a][i];
		      Rj[a] = rr[a][j+1]-rr[a][j]; 
		    }
		  normi=calc_norm(Ri);
		  normj=calc_norm(Rj);
		  for (a=0; a < 3; a++)
		    {
		      Ri[a] /= normi;
		      Rj[a] /= normj;
		    } 
		  for (a=0; a < 3; a++)
		    {
		      costh2 += Ri[a]*Rj[a]; 
		    }	
		}
	      else
		{  
		  for (a = 0; a < 3; a++)
		    {
		      costh2 += u0[a][i]*u0[a][j];
		    }
		}
	      //printf("costh2=%.15G\n", costh2);
	      if (i < NPA)
		GA[np] += costh2;
	      else
		GB[np] += costh2;

	      if (i < NPA)
		{
		  ccA[np] += 1.0;
		}
	      else
		{
		  ccB[np] += 1.0;
		}
	    }
	}
    }
  if (NPA < NP)
    {
      fA = fopen("persistA.dat", "w+");
      fB = fopen("persistB.dat", "w+");
    }
  
  f = fopen("persist.dat", "w+");
  for (ii=1; ii < points; ii++)
    {
      if (NPA < NP)
	{
	  Gm = (GA[ii] + GB[ii]) / (ccA[ii]+ccB[ii]);
	  GmA = GA[ii]/ccA[ii];
	  GmB = GB[ii]/ccB[ii];
	}
      else
	{
	  Gm = GA[ii] / ccA[ii];
	  //printf("ccA[%d]=%.15G\n", ii, ccA[ii]);
	}

      if (ti[ii] > -1.0)
	{
	  if (NPA < NP)
	    {
	      fprintf(fA, "%.15G %.15G\n", ti[ii]-ti[0], GmA);
	      fprintf(fB, "%.15G %.15G\n", ti[ii]-ti[0], GmB); 
	    }
	  fprintf(f, "%.15G %.15G\n", ti[ii]-ti[0], Gm);
	  //fprintf(f, "%.15G %.15G %.15G %.15G\n", ti[ii]-ti[0], C2m, C4m, C6m);
	}
    }
  fclose(f);
  if (NPA < NP)
    {
      fclose(fA);
      fclose(fB);
    }
  return 0;
}
