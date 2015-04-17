#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#define Sqr(x) ((x)*(x))
const double L[3]={168.0,152.0,232.0};
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

char fname[256], strout[1024], line[1024], fnout[1024];
struct anatom *p;
int numatoms, natprot;
double comx, comy, comz;

int cmpfunc(const void *a, const void *b)
{
   return ((struct anatom*)a)->idx - ((struct anatom*)b)->idx;
}
#ifdef ORIENT
#define COM
#include <float.h>
double eigvec[3][3], *eigvec_n[3][3], eigvec_t[3][3];
int eigenvectors=1;
double S=0, Q[3][3];
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

void calcItens(double I[3][3], double RCM[3], int iINI, int iEND)
{
  int i, j, k;
  double distSq, ri[3], rj[3];
  double Icom[3][3];
  for (j=0; j < 3; j++)
    for (k=0; k < 3; k++)
      I[j][k] = 0.0;
  /* moment of inertia of centers of mass */
  for (j=0; j < 3; j++)
    for (k=0; k < 3; k++)
      {
	I[j][k] = 0.0;
	for (i=iINI; i < iEND; i++)
	  {
	    ri[0] = p[i].x-RCM[0];
	    ri[1] = p[i].y-RCM[1];
	    ri[2] = p[i].z-RCM[2];
	    distSq = Sqr(ri[0])+Sqr(ri[1])+Sqr(ri[2]);
	    I[j][k] += 1.0*(((j==k)?distSq:0.0) - ri[j]*ri[k]);
	  }
      }
#if 0
  /* moment of inertia with respect to centes of mass */
  for (i=0; i < Namino; i++)
    {
      tRDiagR(Icom, Iamino[0], Iamino[1], Iamino[2], Ri[i]);
      for (j=0; j < 3; j++)
    	for (k=0; k < 3; k++)
	  {
	    I[j][k] += Icom[j][k]; 
	  }
      }
#endif	
}
#endif
double max3(double a, double b, double c)
{
  double m;
  m = a;
  if (b > m)
    m = b;
  if (c > m)
    m = c;
  return m;
}
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
struct evStruct {
  double eigvec[3];
  double ev;
  int idx;
} evstruct[3];
int cmpfuncev (const void *p1, const void *p2)
{
  if (((struct evStruct*)p1)->ev > ((struct evStruct*)p2)->ev) 
    return 1;
  else
    return -1;
}
int main(int argc, char **argv)
{
  FILE *fin, *fout;
  int i, ifirst, jj, j, nframe, nhdr, a, b;
  double distx, disty, distz;
  double norm, vp[3], pp[3], dx, dy, dz, CM[3], Itens[3][3], ev[3], maxpx, maxpy, maxpz;
  double shift[3];
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
	      if (j > 0)
		{ 
		  shift[0] = L[0]*rint((p[i*natprot+j].x-p[i*natprot].x)/L[0]);
		  shift[1] = L[1]*rint((p[i*natprot+j].y-p[i*natprot].y)/L[1]);
		  shift[2] = L[2]*rint((p[i*natprot+j].z-p[i*natprot].z)/L[2]);
		  p[i*natprot+j].x -= shift[0];
		  p[i*natprot+j].y -= shift[1];
		  p[i*natprot+j].z -= shift[2];
		}
	     
	      comx+=p[i*natprot+j].x;
	      comy+=p[i*natprot+j].y;
	      comz+=p[i*natprot+j].z;
	    }
	  comx /= (double) natprot;
	  comy /= (double) natprot;
	  comz /= (double) natprot;
	  CM[0] = comx;
	  CM[1] = comy;
	  CM[2] = comz;
  	  calcItens(Itens, CM, i*natprot, i*natprot+natprot);
	  diagonalize(Itens, ev);
	  evstruct[0].ev=ev[0];
	  evstruct[1].ev=ev[1];
	  evstruct[2].ev=ev[2];
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      evstruct[a].eigvec[b] = eigvec[a][b];
	  evstruct[a].idx = a;
	  qsort(&evstruct, 3, sizeof(struct evStruct), cmpfuncev);
	  for (j=0; j < natprot; j++)
	    {
	      pp[0] = p[i*natprot+j].x-CM[0];
	      pp[1] = p[i*natprot+j].y-CM[1];
	      pp[2] = p[i*natprot+j].z-CM[2];
	      //vectProdVec(pp,eigvec[0],vp);
	      distx=fabs(scalProd(pp,evstruct[0].eigvec));
	      disty=fabs(scalProd(pp,evstruct[1].eigvec));
	      distz=fabs(scalProd(pp,evstruct[2].eigvec));
	      if (distx > maxpx || j==0)
		maxpx = distx;
	      if (disty > maxpy || j==0)
		maxpy = disty;
	      if (distz > maxpz || j==0)
		maxpz = distz;
	    }
	 
	  fprintf(fout, "%.12G %.12G %.12G %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", comx, comy, comz,
		  evstruct[0].eigvec[0], evstruct[0].eigvec[1], evstruct[0].eigvec[2], 
		  evstruct[1].eigvec[0], evstruct[1].eigvec[1], evstruct[1].eigvec[2], 
		  evstruct[2].eigvec[0], evstruct[2].eigvec[1], evstruct[2].eigvec[2], 
		  sqrt(evstruct[0].ev/((double)natprot)), sqrt(evstruct[1].ev/((double)natprot)), sqrt(evstruct[2].ev/((double)natprot)),
		  maxpx, maxpy, maxpz);
	}
      for (i=0; i < numprot; i++) 
	fprintf(fout, "0 0 0\n");
#else
      for (i=0; i < numatoms; i++) 
	fprintf(fout, "%.12G %.12G %.12G\n", p[i].x, p[i].y, p[i].z);
      for (i=0; i < numatoms; i++) 
	fprintf(fout, "0 0 0\n");
#endif
      fprintf(fout, "%f %f %f\n", L[0], L[1], L[2]);
      //fprintf(fout, "%s\n", line);
      fclose(fout);
      nframe++;
    }
  /* sorting particles */
  fclose(fin);
}
