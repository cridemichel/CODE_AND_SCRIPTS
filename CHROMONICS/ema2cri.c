#include<stdlib.h>
#include<stdio.h>
#include<string.h>
const int Npart=1408;
const double Lbox=8.252567086;
char line[2048];
void writehdr(FILE *f)
{
  fprintf(f, "@@@\n");
  fprintf(f, "parnum: %d\n", Npart);
  fprintf(f, "ntypes: 1\n");
  fprintf(f, "ninters: 3\n");
  fprintf(f, "nintersIJ: 0\n");
  fprintf(f, "saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", Npart);
  fprintf(f, "0.17 0.55 0.55\n");
  fprintf(f, "2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f, "2 0\n");
  fprintf(f, "0.18 0 0 0.3\n");
  fprintf(f, "-0.18 0 0 0.3\n");
  fprintf(f, "0 0 0 0 1 0 0 100000\n");
  fprintf(f, "0 0 0 1 1 0 0 100000\n");
  fprintf(f, "0 1 0 1 1 0 0 100000\n");

  fprintf(f, "@@@\n");
}
int main(int argc, char**argv)
{
  char fname[1024], fon[1024];
  double rx, ry, rz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz;
  int cc=0, t;
  FILE *f, *fo;
  if (argc<2)
    {
      printf("You must supply file name!\n");
      exit(-1);
    }
  strcpy(fname,argv[1]);
  f=fopen(fname,"r");
  sprintf(fon,"Cnf-%d.cnf", cc);
  fo = fopen(fon, "w+");
  writehdr(fo);
  while (!feof(f))
    {
      fscanf(f,"%[^\n]\n",line); 
      if (!strcmp(line,"0"))
	{
	  cc++;
	  fprintf(fo,"10 10 10\n");
	  fclose(fo);
	  sprintf(fon,"Cnf-%d.cnf", cc);
	  fo = fopen(fon, "w+");
	  writehdr(fo);
	}
      else
	{
	  sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &rx, &ry, &rz, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz,
		 &uzx, &uzy, &uzz, &t);
	  fprintf(fo, "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n", rx, ry, rz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, t);
	}
    }
  fprintf(fo, "%.15G %.15G %.15G\n", Lbox, Lbox, Lbox);
  fclose(fo);
  fclose(f);
  exit(0);
}
