#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define Sqr(x) ((x)*(x))
long long int maxtrials=1000000000, tt, overlaps=0;
int outits=500000;
/* ============================ >>> ranf <<< =============================== */
double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

double Lhc, Dhc=2.0;
double PI;
double Lx=5, Ly=5, Lz=5;
int main(int argc, char** argv)
{
  FILE *f;
  int dontcheck;
  char fn[255];
  double angle, norm, pos1[3], pos2[3];
  double n1[3], n2[3], vp[3], fact, delx, rp[3], drp[3];

  int i;
  if (argc==1)
    {
      printf("You have to supplly the angle!\n");
      exit(-1);
    }
  if (argc>1)
    maxtrials=atoi(argv[2]);
  if (argc>2)
    outits=atoi(argv[3]);
  printf("maxtrials=%lld\n", maxtrials);
  PI = 2.0*acos(0.0);
  angle=PI*(180.0 - atof(argv[1]))/180.0;
  n1[0] = cos(angle);
  n1[1] = sin(angle);
  n1[2] = 0.0;
  n2[0] = -1.0;
  n2[1] = 0.0;
  n2[2] = 0.0;

  delx = tan(angle/2.0)*Dhc/2.0; 
  Lhc = (Dhc + delx);
  printf("delx=%f\n", delx);
  /* com1 = uvec[\[Theta]] (X0*D*0.5 - delx*0.5); */
  fact= Dhc*0.5 - delx*0.5;
  pos1[0]=n1[0]*fact;
  pos1[1]=n1[1]*fact;
  pos1[2]=n1[2]*fact;
  printf("n=%f %f %f pos1=%f %f %f (norm=%f)\n", n1[0], n1[1], n1[2], pos1[0], pos1[1], pos1[2], calc_norm(pos1));
  /* com2 = {-(X0 D /2 - delx*0.5), 0, 0};*/
  pos2[0]=-(Dhc/2.0 - delx*0.5);
  pos2[1]=0.0;
  pos2[2]=0.0;
  
  printf("n=%f %f %f pos2=%f %f %f (norm=%f)\n", n2[0], n2[1], n2[2], pos2[0], pos2[1], pos2[2], calc_norm(pos2));
  overlaps=0;
  while (tt < maxtrials)
    {
      if (tt % outits == 0 && overlaps!=0)
	printf("bhc vol= %.15G ( %lld/%lld) \n", (((double)overlaps)/((double)tt))*Lx*Ly*Lz, overlaps, tt);

      /* pick a random point */
      rp[0] = 0.5*Lx*(2.0*ranf()-1.0);
      rp[1] = 0.5*Ly*(2.0*ranf()-1.0);
      rp[2] = 0.5*Lz*(2.0*ranf()-1.0);
      /* check overlap here */
      for (i=0; i < 3; i++)
	drp[i] = rp[i] - pos1[i];
      dontcheck=0;
      if (fabs(scalProd(drp, n1)) < Lhc * 0.5)
	{
	  vectProdVec(drp, n1, vp);
	  norm=calc_norm(vp);
	  if (norm < Dhc*0.5)
	    {
	      overlaps+=1;
	      dontcheck=1;
	    }
	}
      for (i=0; i < 3; i++)
	drp[i] = rp[i] - pos2[i];
      if (dontcheck==0)
	{	
	  if (fabs(scalProd(drp, n2)) < Lhc * 0.5)
	    {
	      vectProdVec(drp, n2, vp);
	      norm=calc_norm(vp);
	      if (norm < Dhc*0.5)
		overlaps+=1;
	    }
	}
      tt++;
    }
  sprintf(fn, "volume-deg%3.0f.dat",angle);
  f = fopen(fn, "w+");
  fprintf(f, "%.15G\n", (((double)overlaps)/((double)tt))*Lx*Ly*Lz);
  fclose(f);
  return 0;
}
