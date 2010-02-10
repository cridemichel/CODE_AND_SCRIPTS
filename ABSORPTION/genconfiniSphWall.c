#include<stdlib.h>
#include<stdio.h>
#include<math.h>
int Nrecept;
int Nprot, brownian, NsphWall, NsphWallOuter;
double bufHeight, Lx, Ly, Lz, sigRecept, massRecept, sigProt, massProt, T, sigSphWall, sigSphWallOuter;
double rx, ry, rz, ox, oy, oz, modr;
double ranf(void)
{
 return rand() / ( (double) RAND_MAX );
}
void solidAngUniform(double *oox, double *ooy, double* ooz)
{
  int a;
  double inert;                 /* momentum of inertia of the molecule */
  double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, Mtot;
  double r21[3], symax[3], wsz, ww[3]; 

  xisq=0.5;
  while (xisq >= 0.5) /* con questa condizione oz > 0 */
    {
      xi1  = ranf() * 2.0 - 1.0;
      xi2  = ranf() * 2.0 - 1.0;
      xisq = xi1 * xi1 + xi2 * xi2;
    }

  xi = sqrt (fabs(1.0 - xisq));
  ox = 2.0 * xi1 * xi;
  oy = 2.0 * xi2 * xi;
  oz = 1.0 - 2.0 * xisq;
  /* Renormalize */
  osq   = ox * ox + oy * oy + oz * oz;
  norm  = sqrt(fabs(osq));
  ox    = ox / norm;
  oy    = oy / norm;
  oz    = oz / norm;
  *oox = ox;
  *ooy = oy;
  *ooz = oz;
  //printf("oo=%.15G %.15G %.15G\n", *oox, *ooy, *ooz);
}

double gauss(void)
{
  double  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0; i < 12; i++)
    {
      sum = sum + ranf();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;

}
int main(int argc, char** argv)
{
 FILE *f;
 double K, DzTop, DzBot, DzTot;
 int i;
 T = 1.0;
 Nrecept = 1;
 Nprot=4000;
 NsphWall = 1;
 NsphWallOuter = 1;
 sigRecept = 1.0;
 sigProt = 0.005;
 massProt = 1.0;
 massRecept = 1.0E200;
 brownian = 1;
 bufHeight = 0.05;
 Lx=Ly=40.0;
 sigSphWall = Lx - bufHeight;
 sigSphWallOuter = sigSphWall + bufHeight;
 Lz=Lx/2.0+bufHeight;
 f=fopen("absorb.cnf","w");
 fprintf(f,"parnum:%d\n", Nrecept+Nprot+NsphWall+NsphWallOuter); 
 fprintf(f,"totStep:%d\n", 20000);
 fprintf(f,"time:%.15G\n",0.0); 
 fprintf(f,"curStep:%d\n", 1); 
 fprintf(f,"P:%.15G\n", 1.0);
 fprintf(f,"T:%.15G\n", 1.0);
 fprintf(f,"rcut:%.15G\n", 1.05);
 fprintf(f,"equilibrat: %d\n", 0); 
 fprintf(f,"Dt:%.15G\n", 0.05); 
 fprintf(f,"ninters:%d\n", 4);
 fprintf(f,"ntypes:%d\n", 5); 
 fprintf(f,"@@@\n");
 fprintf(f,"%d %d %d %d %d\n", Nrecept, Nprot,0,NsphWall,NsphWallOuter);
 /* Receptor type */
 fprintf(f, "%.15G %.15G %.15G\n", sigRecept/2.0, sigRecept/2.0, sigRecept/2.0);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massRecept, 1.0, 1.0, 1.0, 0, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigRecept);
 /* protein type */ 
 fprintf(f, "%.15G %.15G %.15G\n", sigProt/2.0, sigProt/2.0, sigProt/2.0);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massProt, 1.0, 1.0, 1.0, brownian, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigProt);
 /* buffer type */ 
 fprintf(f, "%.15G %.15G %.15G\n", sigProt/2.0, sigProt/2.0, sigProt/2.0);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massProt, 1.0, 1.0, 1.0, brownian, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigProt);
 /* spherical wall type */ 
 fprintf(f, "%.15G %.15G %.15G\n", 0.1, 0.1, 0.1);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massRecept, 1.0, 1.0, 1.0, 0, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigSphWall);
 /* outer spherical wall type */ 
 fprintf(f, "%.15G %.15G %.15G\n", 0.1, 0.1, 0.1);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massRecept, 1.0, 1.0, 1.0, 0, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigSphWallOuter);
 /* interactions */ 
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 0, 1, 0, 1.0, 0.0, 0.0, 1000000);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 1, 0, 3, 0, 0.0001, 0.0, 100000.0, 1000000);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 2, 0, 3, 0, 0.000000003, 0.0, 100000.0, 1000000);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 2, 0, 4, 0, 0.000000003, 0.0, 100000.0, 1000000);
 fprintf(f, "@@@\n");
 /* positions */
 /* receptor */
 fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 0\n", 0.0, 0.0, -Lz*0.5);
 DzTop = 1.2+sigProt*0.5;		
 DzBot = sigProt*0.5+sigRecept*0.5; 
 DzTot = DzTop+DzBot;
 for (i=0; i < Nprot; i++)
   {
#if 0
     fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 1\n", 
	Lx*(ranf()-0.5), Ly*(ranf()-0.5),
	     (Lz-DzTot)*ranf()-Lz*0.5+DzBot); 
#else
     solidAngUniform(&ox, &oy, &oz);
     //printf("o=%f %f %f modo=%.15G\n", ox, oy, oz, sqrt(ox*ox+oy*oy+oz*oz));
     modr = 0.5*(sigSphWall- sigRecept - DzBot)*ranf()+0.5*sigRecept;
     //printf("modr=%.15G\n", modr);
     rx = ox * modr;
     ry = oy * modr;
     rz = oz * modr;
     //fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 1\n", 
     //	     0.0, 0.0, 0.0); 

     fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 1\n", 
	rx, ry, rz - Lz*0.5 + DzBot); 
#endif
   }
 /* spherical walls */
 fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 3\n", 0.0, 0.0, -Lz*0.5);
 fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 4\n", 0.0, 0.0, -Lz*0.5);
 /* velocities */
 /* receptor */
 fprintf(f, "0.0 0.0 0.0 0.0 0.0 0.0\n");
 K = sqrt(T/massProt);
 for (i=0; i < Nprot; i++)
   {
     fprintf(f, "%.15G %.15G %.15G 0.0 0.0 0.0\n", K*gauss(), K*gauss(), K*gauss());
   } 
 /* walls */
 fprintf(f, "0.0 0.0 0.0 0.0 0.0 0.0\n");
 fprintf(f, "0.0 0.0 0.0 0.0 0.0 0.0\n");
 fprintf(f, "%.15G %.15G %.15G\n", Lx, Ly, Lz);
 fclose(f);
}
