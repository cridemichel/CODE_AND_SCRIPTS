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
struct anatom *p;
int numatoms, natprot;

int main(int argc, char **argv)
{
  FILE *fin, *fout;
  int i, ifirst, jj, j, nframe, nhdr;
  double dx, dy, dz;
  strcpy(fname,argv[1]);
  fin=fopen(fname,"r");
  nframe=0;

  nhdr=2; /* numero di linee che costituiscono l'header del file OPEP */
  fscanf(fin, "%d ", &numatoms);
  fscanf(fin,"%[^\n] ",line);
  rewind(fin);
  p=(struct anatom*)malloc(sizeof(struct anatom)*numatoms);

  while (!feof(fin))
    {
      sprintf(fnout,"Store-%d", nframe);
      fout = fopen(fnout, "w+");
      fprintf(fout, "parnum:%d\n", numatoms);
      fprintf(fout, "@@@\n");
	
      fscanf(fin,"%[^\n] ",line);
      fscanf(fin,"%[^\n] ",line);

      for (i=0; i < numatoms; i++) 
	sscanf(line, "%s %lf %lf %lf %d ", p[i].resname, &(p[i].x), &(p[i].y), &(p[i].z),
	       &(p[i].idx));
      qsort(p, numatoms,sizeof(struct anatom), cmpfunc);
      for (i=0; i < numatoms; i++) 
	fprintf(fout, "%12s %12lf %12lf %12lf %12d\n", p[i].resname, p[i].x, p[i].y, p[i].z, p[i].idx);
      //fprintf(fout, "%s\n", line);
      fclose(fout);
      nframe++;
    }
  /* sorting particles */
  fclose(fin);
}
