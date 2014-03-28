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
#define nmaxNxi 1000
const int Nxi=15;
const double sighelix=1.0, pitch=4.0, radhelix=0.2, lenhelix=10.0;
double xihel[nmaxNxi], xhel[nmaxNxi][3], xhelA[nmaxNxi][3], xhelB[nmaxNxi][3];
void build_helix(void)
{
  double temp, pi, npitch, deltaxi, length_eucl, radius;
  int jj, kk;
  double xcm[3];

  radius = radhelix; /* x,y perpendicular to helix axis in body reference frame */
  length_eucl = lenhelix; /* helix axis along z in body reference frame */
  pi = acos(0.0)*2.0;
  npitch = length_eucl/pitch;
  temp=npitch*sqrt(Sqr(radius)+Sqr(pitch/(2.0*pi)));
  deltaxi=2.0*pi*npitch/((double)(Nxi-1));
 // printf("pi=%f npitch=%f Nxi-1=%d\n", pi, npitch, Nxi-1);
 // printf("deltaxi=%f\n", deltaxi);
  // length_eucl=OprogStatus.npitch*OprogStatus.pitch;
  for (jj=0; jj < Nxi; jj++)
    {
      xihel[jj]=((double)jj)*deltaxi;
      xhel[jj][0]=radius*cos(xihel[jj]);
      xhel[jj][1]=radius*sin(xihel[jj]);
      xhel[jj][2]=pitch*xihel[jj]/(2.0*pi);
    }
  xcm[0]=xcm[1]=xcm[2]=0.0;
  for (jj=0; jj < Nxi; jj++)
    {
      for (kk=0; kk < 3; kk++)
	xcm[kk] += xhel[jj][kk]; 
    } 

  for (kk=0; kk < 3; kk++)
    xcm[kk] /= Nxi;
  for (jj=0; jj < Nxi; jj++)
    {
      for (kk=0; kk < 3; kk++)
	xhel[jj][kk] -= xcm[kk];
    }   
}
#define Sqr(x) ((x)*(x))
double Lhc, Dhc=2.0;
double PI;
double Lx=2, Ly=2, Lz=11;
int main(int argc, char** argv)
{
  FILE *f;
  int jj;
  char fn[255];
  double angle, norm, pos1[3], pos2[3];
  double n1[3], n2[3], vp[3], fact, delx, rp[3], drp[3];

  int i;
  if (argc==1)
    {
      printf("You have to supply the angle!\n");
      exit(-1);
    }
  if (argc>1)
    maxtrials=atoi(argv[1]);
  if (argc>2)
    outits=atoi(argv[2]);
  printf("maxtrials=%lld\n", maxtrials);
  /* com1 = uvec[\[Theta]] (X0*D*0.5 - delx*0.5); */
 // printf("n=%f %f %f pos2=%f %f %f (norm=%f)\n", n2[0], n2[1], n2[2], pos2[0], pos2[1], pos2[2], calc_norm(pos2));

  overlaps=0;
  build_helix();
  while (tt < maxtrials)
    {
      if (tt % outits == 0 && overlaps!=0)
	printf("bhc vol= %.15G ( %lld/%lld) \n", (((double)overlaps)/((double)tt))*Lx*Ly*Lz, overlaps, tt);

      /* pick a random point */
      rp[0] = 0.5*Lx*(2.0*ranf()-1.0);
      rp[1] = 0.5*Ly*(2.0*ranf()-1.0);
      rp[2] = 0.5*Lz*(2.0*ranf()-1.0);
      /* check overlap here */
      for (jj=0; jj < Nxi; jj++)
	{
	  if (Sqr(rp[0]-xhel[jj][0])+Sqr(rp[1]-xhel[jj][1])+Sqr(rp[2]-xhel[jj][2]) <= Sqr(0.5*sighelix))
	    {
	      overlaps+=1.0;
	      break;
	    }
	}
      tt++;
    }
  sprintf(fn, "volume-deg%s.dat",argv[1]);
  f = fopen(fn, "w+");
  fprintf(f, "%.15G\n", (((double)overlaps)/((double)tt))*Lx*Ly*Lz);
  fclose(f);
  return 0;
}
