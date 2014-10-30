#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
const double L[3]={260.0,260.0,260.0};
const int numprot=64;
const int stepinterval=20000;
const double dt=0.5;
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
double comx, comy, comz;
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

  natprot=numatoms/numprot;

  while (!feof(fin))
    {
      sprintf(fnout,"Store-%d", nframe);
      fout = fopen(fnout, "w+");
      fprintf(fout, "initFormat:0\n");
      fprintf(fout, "@@@\n");
#ifdef COM
      fprintf(fout, "parnum:%d\n", numprot);
#else
      fprintf(fout, "parnum:%d\n", numatoms);
#endif
      fprintf(fout, "steplength:%f\n", dt);
      fprintf(fout, "curStep: %d\n", stepinterval*nframe); 
      fprintf(fout, "@@@\n");
	
      fscanf(fin,"%[^\n] ",line);
      fscanf(fin,"%[^\n] ",line);

      for (i=0; i < numatoms; i++) 
	{
	  fscanf(fin, " %s %lf %lf %lf %d ", p[i].resname, &(p[i].x), &(p[i].y), &(p[i].z),
	     	 &(p[i].idx));
	}
      qsort(p, numatoms, sizeof(struct anatom), cmpfunc);
#ifdef COM
      for (i=0; i < numprot; i++)
	{
	  comx=comy=comz=0.0;
	  for (j=0; j < natprot; j++)
	    {
	      comx+=p[i*natprot+j].x;
	      comy+=p[i*natprot+j].y;
	      comz+=p[i*natprot+j].z;
	    }
	  comx /= (double) natprot;
	  comy /= (double) natprot;
	  comz /= (double) natprot;
	  fprintf(fout, "%.12G %.12G %.12G\n", comx, comy, comz);
	}
      for (i=0; i < numprot; i++) 
	fprintf(fout, "0 0 0\n");
#else
      for (i=0; i < numatoms; i++) 
	fprintf(fout, "%.12G %.12G %.12G\n", p[i].x, p[i].y, p[i].z);
      for (i=0; i < numatoms; i++) 
	fprintf(fout, "0 0 0\n");
#endif
      fprintf(fout, "%f\n", L[0]);
      //fprintf(fout, "%s\n", line);
      fclose(fout);
      nframe++;
    }
  /* sorting particles */
  fclose(fin);
}
