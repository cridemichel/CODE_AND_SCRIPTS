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
char fname[256], dummy[1024];
struct anatom *p=NULL;
int numatoms, natprot;

int main(int argc, char **argv)
{
  FILE *fin;
  int nr=0, i, ifirst, jj, j;
  double dx, dy, dz;
  strcpy(fname,argv[1]);
  fin=fopen(fname,"r");
  while (!feof(fin))
    {
      if (nr==1) 
	fscanf(fin,"%[^\n] ",dummy);
      else if (nr==0)
	{
	  fscanf(fin, "%d ", &numatoms);
	  p=(struct anatom*)malloc(sizeof(struct anatom)*numatoms);
	  fprintf(stderr,"numatoms=%d\n", numatoms);
	}
      else 
	{
	  //printf("nr-2=%d\n", nr-2);
	  fscanf(fin, "%s %lf %lf %lf %d ", p[nr-2].resname, &(p[nr-2].x), &(p[nr-2].y), &(p[nr-2].z),
		 &(p[nr-2].idx));
	}
      nr++;
    }
  /* sorting particles */
  qsort(p, numatoms, sizeof(struct anatom), cmpfunc);
  fclose(fin);
  natprot=numatoms/numprot;
  fprintf(stderr,"natprot=%d\n", natprot);
  for (i=0; i < numprot; i++)
    {
      ifirst= i*natprot;
      fprintf(stderr,"ifirst=%d\n", ifirst);
      for (j=1; j < natprot; j++)
	{
	  jj = ifirst + j;
	  dx = p[jj].x - p[ifirst].x;
	  dy = p[jj].y - p[ifirst].y;
	  dz = p[jj].z - p[ifirst].z;
	  p[jj].x = p[jj].x - L[0]*rint(dx/L[0]);
	  p[jj].y = p[jj].y - L[1]*rint(dy/L[1]);
	  p[jj].z = p[jj].z - L[2]*rint(dz/L[2]);
	}
    }

  printf("%12d\n", numatoms);
  printf("%12s\n", dummy);
  for (i=0; i < numatoms; i++)
    {
      printf("%12s %12lf %12lf %12lf %12d\n", p[i].resname, p[i].x, p[i].y, p[i].z, p[i].idx);
    }
}
