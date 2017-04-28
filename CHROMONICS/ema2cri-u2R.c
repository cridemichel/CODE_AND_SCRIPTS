#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
const int Npart=2;
const double Lbox=20.0;
double R[3][3];
double rx, ry, rz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, ox, oy, oz;
char line[2048];
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */

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
void R2u(void)
{
  uxx = R[0][0];
  uxy = R[0][1];
  uxz = R[0][2];
  uyx = R[1][0];
  uyy = R[1][1];
  uyz = R[1][2];
  uzx = R[2][0];
  uzy = R[2][1];
  uzz = R[2][2];
}

void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
#ifdef MC_BENT_DBLCYL
  double Rout[3][3];
  int k1, k2;
#endif
  /* first row vector */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[0][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[0][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
#if 0
  if (typesArr[0].nspots==3 && type==0)
    {
      for (k=0; k < 3 ; k++)
	u[k] = R[1][k];
      vectProdVec(R[0], u, up);
      /* rotate randomly second axis */
      angle=4.0*acos(0.0)*ranf_vb();
      xx = cos(angle);
      yy = sin(angle);
      for (k=0; k < 3 ; k++)
	R[1][k] = u[k]*xx + up[k]*yy;
      //printf("calc_norm(R[1])=%.15G\n", calc_norm(R[1]));
    }
#endif
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
#ifdef MC_BENT_DBLCYL
  /* add a random rotation around the axis (ox, oy, oz) */
  add_rotation_around_axis(ox, oy, oz, R, Rout);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R[k1][k2] = Rout[k1][k2];
#endif
#if 0
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    Rt[k1][k2]=R[k2][k1];
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    R[k1][k2]=Rt[k1][k2];
#endif
  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
}
void writehdr(FILE *f)
{
  fprintf(f, "@@@\n");
  fprintf(f, "parnum: %d\n", Npart);
  fprintf(f, "ntypes: 1\n");
  fprintf(f, "ninters: 0\n");
  fprintf(f, "nintersIJ: 0\n");
  fprintf(f, "saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", Npart);
  fprintf(f, "0.17 0.55 0.55\n");
  fprintf(f, "2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f, "0 0\n");

  fprintf(f, "@@@\n");
}
int main(int argc, char**argv)
{
  char fname[1024], fon[1024];
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
	  sscanf(line,"%lf %lf %lf %lf %lf %lf\n", &rx, &ry, &rz,  &ox, &oy, &oz);
	  versor_to_R(ox,oy,oz, R);
	  R2u();
	  fprintf(fo, "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n", rx, ry, rz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, t);
	}
    }
  fprintf(fo, "%.15G %.15G %.15G\n", Lbox, Lbox, Lbox);
  fclose(fo);
  fclose(f);
  exit(0);
}
