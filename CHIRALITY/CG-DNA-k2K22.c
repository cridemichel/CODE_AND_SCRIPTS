//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}

char fn[1024];
struct DNA {
  double x;
  double y;
  double z;
  double rad;
} *DNAchain;

struct DNA *DNADs[2];
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

#define MC_BENT_DBLCYL

char dummy1[32], dummy2[32], atname[32], nbname[8];
int nat, atnum, nbnum, len, tot_trials, tt=0;
double L, rx, ry, rz, alpha;
/*
                pdb        radius (angstrom)
    sugar       Xe          3.5        
    phosphate   B           3.0         
    base        Se          4.0    
*/
void body2lab(double xp[3], double x[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
      x[k1] += rO[k1];
    }
}
#ifdef MC_BENT_DBLCYL
/* apply a random rotation around the supplied axis because 
   bent cylinders do not have azimuthal symmetry */
double thetaGlobalBondangle;
void add_rotation_around_axis(double ox, double oy, double oz, double Rin[3][3], double Rout[3][3])
{
  double theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random rotation angle between 0 and 2*pi*/
  theta = 4.0*acos(0.0)*drand48();
  /* set to be used in az. angle distro calculation */
  thetaGlobalBondangle = theta;

  thetaSq=Sqr(theta);
  sinw = sin(theta);
  cosw = (1.0 - cos(theta));
  Omega[0][0] = 0;
  Omega[0][1] = -oz;
  Omega[0][2] = oy;
  Omega[1][0] = oz;
  Omega[1][1] = 0;
  Omega[1][2] = -ox;
  Omega[2][0] = -oy;
  Omega[2][1] = ox;
  Omega[2][2] = 0;
  OmegaSq[0][0] = -Sqr(oy) - Sqr(oz);
  OmegaSq[0][1] = ox*oy;
  OmegaSq[0][2] = ox*oz;
  OmegaSq[1][0] = ox*oy;
  OmegaSq[1][1] = -Sqr(ox) - Sqr(oz);
  OmegaSq[1][2] = oy*oz;
  OmegaSq[2][0] = ox*oz;
  OmegaSq[2][1] = oy*oz;
  OmegaSq[2][2] = -Sqr(ox) - Sqr(oy);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = Rin[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += Rin[k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     Rout[k1][k2] = Ro[k1][k2]; 
}
#endif

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

void place_DNAD(double x, double y, double z, double ux, double uy, double uz, int which)
{
  double xp[3], rO[3], xl[3];
  double R[3][3];
  int i; 
  rO[0] = x;
  rO[1] = y;
  rO[2] = z;
  /* build R here from the orientation (ux,uy,uz) */
  versor_to_R(ux, uy, uz, R);
  /* ============ */
  for (i=0; i < nat; i++)
    {
      xp[0] = DNAchain[i].x;
      xp[1] = DNAchain[i].y;
      xp[2] = DNAchain[i].z;
      body2lab(xp, xl, rO, R);
      DNADs[which][i].x = xl[0];
      DNADs[which][i].y = xl[1];
      DNADs[which][i].z = xl[2];
      DNADs[which][i].rad = DNAchain[i].rad;
    }
}
/* ============================ >>> ranf <<< =============================== */
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}

double fons(double theta, double alpha)
{
  double pi;
  pi = acos(0.0)*2.0;
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return cosh(alpha*cos(theta))*alpha/(4.0*pi*sinh(alpha));
}
/* return an angle theta sampled from an Onsager angular distribution */
double theta_onsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta, dtheta;
  //printf("alpha=%f\n", alpha);
  pi = acos(0.0)*2.0;
  if (first == 1)
    {
      first=0;
      f0 = 1.01*fons(0.0,alpha);
    }

  do 
    {
      /* uniform theta between 0 and pi */
      theta = pi*ranf_vb();
      /* uniform y between 0 and 1 (note that sin(theta) <= 1 for 0 < theta < pi)*/
      y = f0*ranf_vb();
      f = sin(theta)*fons(theta,alpha);
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}
double distro[10000];
const int nfons=100;
void orient_onsager(double *omx, double *omy, double* omz, double alpha)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_onsager(alpha);
  //printf("thos=%f\n", thons);
  distro[(int) (thons/(pi/((double)nfons)))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1;
#if 1 /* along z */
  *omx = verso*sin(thons)*cos(phi);
  *omy = verso*sin(thons)*sin(phi);
  *omz = verso*cos(thons); 
#else /* or along x (but it has to be same of course!) */
  *omy = verso*sin(thons)*cos(phi);
  *omz = verso*sin(thons)*sin(phi);
  *omx = verso*cos(thons); 
#endif
  //printf("norma=%f\n", sqrt(Sqr(*omx)+Sqr(*omy)+Sqr(*omz)));
}

/* first derivative of Onsager distribution */
double dfons(double theta, double alpha)
{
  double pi;
  pi = acos(0.0)*2.0;
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return sin(theta)*cosh(alpha*cos(theta))*alpha*alpha/(4.0*pi*sinh(alpha));
}

/* return an angle theta sampled from an Onsager angular distribution */
double theta_donsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta, dtheta;
  //printf("alpha=%f\n", alpha);
  pi = acos(0.0)*2.0;
  if (first == 1)
    {
      first=0;
      /* qui va fornito il massimo della funzione nell'intervallo
	 [0,Pi] */
      f0 = 1.01*alpha*fons(0.0,alpha);
    }

  do 
    {
      /* uniform theta between 0 and pi */
      theta = pi*ranf_vb();
      /* uniform y between 0 and 1 (note that sin(theta) <= 1 for 0 < theta < pi)*/
      y = f0*ranf_vb();
      f = sin(theta)*dfons(theta,alpha);
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}

extern const int nfons;
void orient_donsager(double *omx, double *omy, double* omz, double alpha)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_donsager(alpha);
  //printf("thos=%f\n", thons);
  //distro[(int) (thons/(pi/((double)nfons)))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1;
#if 1 /* along z */
  *omx = verso*sin(thons)*cos(phi);
  *omy = verso*sin(thons)*sin(phi);
  *omz = verso*cos(thons); 
#else /* or along x (but it has to be same of course!) */
  *omy = verso*sin(thons)*cos(phi);
  *omz = verso*sin(thons)*sin(phi);
  *omx = verso*cos(thons); 
#endif
  //printf("norma=%f\n", sqrt(Sqr(*omx)+Sqr(*omy)+Sqr(*omz)));
}

int main(int argc, char**argv)
{
  FILE *fin, *fout;
  int k, i, j, overlap, type, outits, thetapts, fileoutits;
  char fnin[1024],fnout[256];
  double u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, factor, dth, th;
  /* syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */
  if (argc < 7)
    {
      printf("syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <fileoutits> [outits]\n");
      exit(1);
    }
  strcpy(fnin,argv[1]);
  fin=fopen(fnin,"r");
  len=atoi(argv[2]);
  alpha = atof(argv[3]);
  tot_trials=atoi(argv[4]);
  type = atoi(argv[5]);
  fileoutits = atoi(argv[6]);
  
  if (argc == 7)
    outits=100*fileoutits;
  else
    outits = atoi(argv[7]);
   /* ATOM    39   Xe   G A   14      -5.687  -8.995  37.824 */
  nat = 70*len;
  DNAchain = (struct DNA*) malloc(sizeof(struct DNA)*nat); 
  for (k=0; k < 2; k++)
    DNADs[k] = (struct DNA*) malloc(sizeof(struct DNA)*nat);
  L = 1.05*2.0*40*len; /* 4 nm is approximately the length of a 12 bp DNAD */ 
  /* read the CG structure */
  while (!feof(fin))
    {
      fscanf(fin, "%s %d %s %s %s %d %lf %lf %lf ", dummy1, &atnum, atname, nbname, dummy2, &nbnum, &rx, &ry, &rz);
      DNAchain[atnum].x = rx;
      DNAchain[atnum].y = ry;
      DNAchain[atnum].z = rz;
      if (!strcmp(atname, "Xe"))
	{
	  DNAchain[atnum].rad = 3.5;
	}
      else if (!strcmp(atname, "B"))
	{
	  DNAchain[atnum].rad = 3.0;
	}
      else if (!strcmp(atname, "Se"))
	{
	  DNAchain[atnum].rad = 4.0;
	}
      else
	{
	  printf("Unrecognized atom name, exiting...\n");
	  exit(1);
	}
    };
  printf("nat=%d L=%f alpha=%f I am going to calculate v%d and I will do %d trials\n", nat, L, alpha, type, tot_trials);
  rcmx=rcmy=rcmz=0.0;
  for (i=0; i < nat; i++)
    {
      rcmx += DNAchain[i].x;
      rcmy += DNAchain[i].y;
      rcmz += DNAchain[i].z;
    }
  rcmx /= (double) nat;
  rcmy /= (double) nat;
  rcmz /= (double) nat;
  for (i=0; i < nat; i++)
    {
      DNAchain[i].x -= rcmx;
      DNAchain[i].y -= rcmy;
      DNAchain[i].z -= rcmz;
    }
  fclose(fin);
  srand48((int)time(NULL));
  sprintf(fnout, "v%d.dat", type);
  factor=0.0;

  thetapts = 100000;
  dth=2.0*(acos(0.0))/((double)thetapts);
  th=0.0;
  for (i=0; i < thetapts; i++)
    {
      factor += 0.5*dth*sin(th)*(dfons(th, alpha) + dfons(th+dth,alpha));
      th += dth;
      //printf("%f %.15G\n", th, dfons(th, alpha));
    }
  factor *= 4.0*acos(0.0);
  printf("dth=%f factor=%.15G\n", dth, factor);
  fout = fopen(fnout, "w+");
  fclose(fout);  
  for (tt=0; tt < tot_trials; tt++)
    {
      /* place first DNAD in the origin oriented according to the proper distribution */

      if (type==0||type==1)
	orient_onsager(&u1x, &u1y, &u1z, alpha);
      else
	orient_donsager(&u1x, &u1y, &u1z, alpha);
      place_DNAD(0.0, 0.0, 0.0, u1x, u1y, u1z, 0);      
      /* place second DNAD randomly */
      rcmx = L*(drand48()-0.5);
      rcmy = L*(drand48()-0.5);
      rcmz = L*(drand48()-0.5);
      if (type==0)
	orient_onsager(&u2x, &u2y, &u2z, alpha);
      else
	orient_donsager(&u2x, &u2y, &u2z, alpha);
      place_DNAD(rcmx, rcmy, rcmz, u2x, u2y, u2z, 1);
      /* check overlaps */
      overlap=0;
      for (i=0; i < nat; i++)
	{
	  for (j=0; j < nat; j++)
	    {
	      distsq = Sqr(DNADs[0][i].x-DNADs[1][j].x) + Sqr(DNADs[0][i].y-DNADs[1][j].y) + Sqr(DNADs[0][i].z-DNADs[1][j].z) ;
	      sigijsq = Sqr(DNADs[0][i].rad + DNADs[1][j].rad);
	      if (distsq < sigijsq)
		{
		  overlap=1;
		  break;
		}
	    }
	  if (overlap)
	    break;
	}
      if (overlap)
	continue;
      /* otherwise calculate the integrand */
      if (type==0)
	vexcl += 1.0;
      else if (type==1)
	vexcl += -u2x*rcmy; /* questo '-' rende negativa la k2 e viene dalla derivata della funzione di Onsager! */
      else 
	vexcl += u1x*u2x*rcmy*rcmy;
      if (tt % fileoutits == 0)
	{
	 fout = fopen(fnout, "a+");
	 if (type==0)
	   fprintf(fout,"%d %.15G\n", tt, L*L*L*vexcl/((double)tt));
	 else if (type==1)
	   fprintf(fout,"%d %.15G\n", tt, L*L*L*vexcl/((double)tt)/factor);
	 else
	   fprintf(fout,"%d %.15G\n", tt, L*L*L*vexcl/((double)tt)/Sqr(factor));
	 fclose(fout);
	}
      if (tt % outits==0)
	printf("trials: %d/%d\n", tt, tot_trials);
    }
}
