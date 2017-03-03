#include<stdlib.h>
#include<stdio.h>
#include<string.h>
const int Npart=8000;
char line[2048];
void writehdr(FILE *f)
{
  fprintf(f, "@@@\n");
  fprintf(f, "parnum: 8000\n");
  fprintf(f, "@@@\n");
  fprintf(f, "@@@\n");
}
int main(int argc, char**argv)
{
  char fname[1024], fon[1024];
  double rx, ry, rz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz;
  int cc=0;
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
	  sscanf(line,"%lf %lf %lf %lf %lf %lf\n", &rx, &ry, &rz, &uxx, &uxy, &uxz, &ayx, &uyy, &uyz,
		 &uzx, &uzy, &uzz);
	  fprintf(fo, "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", rx, ry, rz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz);
	}
    }
  fclose(fo);
  fclose(f);
  exit(0);
}
