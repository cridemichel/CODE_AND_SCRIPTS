#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
const double L[3]={120.0,120.0,120.0};
const int numprot=8;
struct anatom {
  char resname[32];
  double x;
  double y;
  double z;
  int idx;
} anatomStruct;
int cmpfunc(const void *a, const void *b)
{
   return ((struct anatom*)a)->idx - ((struct anatom*)b)->idx;
}
char fname[256], strout[1024], line[1024], fnout[1024];
struct anatom p;
int numatoms, natprot;

int main(int argc, char **argv)
{
  FILE *fin, *fout;
  int nr=0, i, ifirst, jj, j, nframe;
  double dx, dy, dz;
  strcpy(fname,argv[1]);
  fin=fopen(fname,"r");
  nframe=0;

  fscanf(fin, "%d ", &numatoms);
  fscanf(fin,"%[^\n] ",line);
  rewind(fin);
  while (!feof(fin))
    {
      if (nr % numatoms == 0)
	{
	  sprintf(fnout,"Store-%d", nframe);
	  fout = fopen(fnout, "w+");
	  fprintf(fout, "parnum:%d\n", numatoms);
	  fprintf(fout, "@@@\n");
	}
      fscanf(fin,"%[^\n] ",line);
      
      if (nr % numatoms > 1)
	{
	  sscanf(line, "%s %lf %lf %lf %d ", p.resname, &(p.x), &(p.y), &(p.z),
		 &(p.idx));
  	  fprintf(fout, "%lf %lf %lf\n", p.x, p.y, p.z);
	}
      //fprintf(fout, "%s\n", line);
      nr++;
      if (nr % numatoms ==0)
	{
	  fclose(fout);
	  nframe++;
	}	
    }
  /* sorting particles */
  fclose(fin);
}
