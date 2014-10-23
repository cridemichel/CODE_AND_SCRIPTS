#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
const double L[3]={120.0,120.0,120.0};
const int numprot=8;
int numatprot=0;
const int stepinterval=20000;
const double dt=0.5;
const double expandFact=10.0;
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
double *comx, *comy, *comz;
int main(int argc, char **argv)
{
  FILE *fin, *fout;
  int k, i, ifirst, jj, j, nframe, nhdr;
  double dx, dy, dz;
  strcpy(fname,argv[1]);
  fin=fopen(fname,"r");
  nframe=0;

  comx = malloc(sizeof(double)*numprot);
  comy = malloc(sizeof(double)*numprot);
  comz = malloc(sizeof(double)*numprot);

  nhdr=2; /* numero di linee che costituiscono l'header del file OPEP */
  fscanf(fin, "%d ", &numatoms);
  fscanf(fin,"%[^\n] ",line);
  rewind(fin);
  p=(struct anatom*)malloc(sizeof(struct anatom)*numatoms);
  numatprot=numatoms/numprot;

  while (!feof(fin))
    {
      sprintf(fnout,"StoreFF-%d", nframe);
      fout = fopen(fnout, "w+");
      fprintf(fout, "initFormat:0\n");
      fprintf(fout, "@@@\n");
      fprintf(fout, "parnum:%d\n", numatoms);
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
      for (k=0; k < numprot; k++)
	{
	  comx[k]=0.0;
	  comy[k]=0.0;
	  comz[k]=0.0;
	}	  
      for (i=0; i < numatoms; i++) 
	{
	  k = i / numatprot;
	  comx[k]+=p[i].x;
	  comy[k]+=p[i].y;
	  comz[k]+=p[i].z;
	}
      for (k=0; k < numprot; k++)
	{
	  comx[k] /= ((double)numatprot);
	  comy[k] /= ((double)numatprot);
	  comz[k] /= ((double)numatprot);
	}
      qsort(p, numatoms, sizeof(struct anatom), cmpfunc);
      for (i=0; i < numatoms; i++) 
	{
	  k = i / numatprot;
	  fprintf(fout, "%.12G %.12G %.12G\n", (p[i].x-comx[k])+comx[k]*expandFact, (p[i].y-comy[k])+comy[k]*expandFact, 
		  (p[i].z-comz[k])+comz[k]*expandFact);
	}
      for (i=0; i < numatoms; i++) 
	fprintf(fout, "0 0 0\n");
      fprintf(fout, "%f\n", L[0]*expandFact);
      //fprintf(fout, "%s\n", line);
      fclose(fout);
      nframe++;
    }
  /* sorting particles */
  fclose(fin);
}
