#define MD_SP_DELR 0.0
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2

double spApos[MD_STSPOTS_A][3] = {{MD_SP_DELR, 0.54, 0.0},{MD_SP_DELR, 0.54, 3.14159},{MD_SP_DELR, 2.60159,0.0},
    {MD_SP_DELR, 2.60159, 3.14159},{MD_SP_DELR, 1.5708, 0.0}};
double spBpos[MD_STSPOTS_B][3] = {{MD_SP_DELR, 0.0, 0.0},{MD_SP_DELR, 3.14159, 0.0}};
double spXYZ_A[5][3];
double spXYZ_B[2][3];
double theta, Dr, sigmaSticky;
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double a[2]={1.0,2.0}, b[2]={1.0,2.0}, c[2]={5.0,10.0};
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}


void build_atom_positions(void)
{
 /* N.B. le coordinate spXpos sono del tipo (Dr, theta, phi),
  * dove se Dr=0 la sfera sticky viene posizionata esattamente in 
  * maniera tangente e theta (0 <= theta <= Pi) e phi (0 <= phi < 2Pi)
  * sono gli angoli in coordinate sferiche che individuano il punto di contatto
  * tra sticky sphere ed ellissoide.
  * Tale routine converte le coordinate spXpos in coordinate cartesiane 
  * riferite al riferimento del corpo rigido. */  
  int kk, k1, aa;
  double pi, x,y,z, grad[3], ng, dd[3];

  sigmaSticky = 0.4;
  Dr = -0.2;
  theta = 0.54;
  pi = acos(0.0)*2.0;
  spApos[0][1] = theta;
  spApos[1][1] = theta;
  spApos[2][1] = pi - theta;
  spApos[3][1] = pi - theta;
  spApos[0][0] = Dr;
  spApos[1][0] = Dr;
  spApos[2][0] = Dr;
  spApos[3][0] = Dr;
  spApos[4][0] = Dr;
  spBpos[0][0] = Dr;
  spBpos[1][0] = Dr;
  for (k1 = 0; k1 < MD_STSPOTS_A; k1++)
    {
      x = a[0]*cos(spApos[k1][2])*sin(spApos[k1][1]);
      y = b[0]*sin(spApos[k1][2])*sin(spApos[k1][1]);
      z = c[0]*cos(spApos[k1][1]);
      //printf("xyz=%f %f %f\n", x, y, z);
      grad[0] = 2.0 * x / Sqr(a[0]);
      grad[1] = 2.0 * y / Sqr(b[0]);
      grad[2] = 2.0 * z / Sqr(c[0]);
      ng = calc_norm(grad);
      for (aa = 0; aa < 3; aa++)
	grad[aa] /= ng;
      spXYZ_A[k1][0] = x + grad[0]*(sigmaSticky*0.5 + spApos[k1][0]);
      spXYZ_A[k1][1] = y + grad[1]*(sigmaSticky*0.5 + spApos[k1][0]);
      spXYZ_A[k1][2] = z + grad[2]*(sigmaSticky*0.5 + spApos[k1][0]);
	
	      //printf("k1=%d %f %f %f \n", k1,  spXYZ_A[k1][0] ,    spXYZ_A[k1][1] ,  spXYZ_A[k1][2]  );
    }
  for (kk=0; kk < 3; kk++)
    dd[kk] = spXYZ_A[0][kk] - spXYZ_A[1][kk];
  printf("Molecule A distance between Atoms 0 and 1: %.15G\n", calc_norm(dd));;
  for (kk=0; kk < 3; kk++)
    dd[kk] = spXYZ_A[2][kk] - spXYZ_A[3][kk];
  printf("Molecule A distance between Atoms 2 and 3: %.15G\n", calc_norm(dd));;

  for (k1 = 0; k1 < MD_STSPOTS_B; k1++)
    {
      x = a[1]*cos(spBpos[k1][2])*sin(spBpos[k1][1]);
      y = b[1]*sin(spBpos[k1][2])*sin(spBpos[k1][1]);
      z = c[1]*cos(spBpos[k1][1]);
      grad[0] = 2.0 * x / Sqr(a[1]);
      grad[1] = 2.0 * y / Sqr(b[1]);
      grad[2] = 2.0 * z / Sqr(c[1]);
      ng = calc_norm(grad);
      for (aa = 0; aa < 3; aa++)
	grad[aa] /= ng;
      spXYZ_B[k1][0] = x + grad[0]*(sigmaSticky*0.5 + spBpos[k1][0]);
      spXYZ_B[k1][1] = y + grad[1]*(sigmaSticky*0.5 + spBpos[k1][0]);
      spXYZ_B[k1][2] = z + grad[2]*(sigmaSticky*0.5 + spBpos[k1][0]) ;
    }
  printf("PARTICLES A\n\n");
  for (k1 = 0; k1 < MD_STSPOTS_A; k1++)
    {
      printf("%.15G %.15G %.15G %f\n", spXYZ_A[k1][0], spXYZ_A[k1][1], spXYZ_A[k1][2], sigmaSticky);
    }
  printf("\n PARTICLES B\n\n");
  for (k1 = 0; k1 < MD_STSPOTS_B; k1++)
    {
      printf("%.15G %.15G %.15G %f\n", spXYZ_B[k1][0], spXYZ_B[k1][1], spXYZ_B[k1][2], sigmaSticky);
    }
}
main ()
{
build_atom_positions();
}
