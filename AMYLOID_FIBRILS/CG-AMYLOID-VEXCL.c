//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#define USEGSL
//#define ALBERTA
#ifdef USEGSL
#include <gsl/gsl_qrng.h>
#endif
void print_matrix(double M[3][3])
{
  int k1, k2;
  int n=3;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}

static double maxarg1,maxarg2;
double CHROMheight=0.34;
#ifdef AMYLOID_ELEC
struct CHROM {
  double x;
  double y;
  double z;
  double rad;
  int atype;
} *CHROMchain;
struct CHROM *CHROMs[2];
#endif
#ifdef AMYLOID_ELEC
double kD, yukcut, yukcutkD, yukcutkDsq;
#ifdef PARALLEL
struct kDsortS {
int k1;
int k2;
double invkD;
} *kD_sorted;
double *cchrom_arr, *beta_arr;
int numtemps, numconcs;
double **yukcutkD_arr, **kD_arr, **yuk_corr_fact_arr, **yukcutkDsq_arr, **uel_arr, **vexclel_arr;
double num_kD=0;
double maxyukcutkDsq, maxyukcutkD;
#endif
double delta_rab0=0.2/*in nm!*/, epsr_prime=1.8, yuk_corr_fact;
double esq_eps, esq_eps_prime; /* = e^2 / (4*pi*epsilon0*epsilon*kB) in J*m */
double esq_eps10, esq_eps_prime10;
const double bmann = 1E-9*0.34/2.0; /* spacing between charged phosphate groups for manning theory */ 
const double Dalton = 1.660538921E-27;
const double kB = 1.3806503E-23, eps0=8.85E-12; /* boltzmann constant */
const double qel = 1.602176565E-19, Nav=6.02214129E23;
const double qchrom = 1.0, qsalt = 1.0; /* qsalt è la valenza del sale aggiunto (tipicamente 1 poiché si tratta di NaCl */
const double CHROMmolmass = 452.0; /* espressa in Dalton, la massa molare fornita da Lavrentovich è 0.452 Kg/mol */
double cchrom, csalt = 0.0; /* concentrazione del sale aggiunto molare */
double ximanning, deltamann; /* Debye screening length */
/* charge on phosphate groups */
double zeta_a, zeta_b;

double Ucoul(double rab)
{
  //return esq_eps_prime10/rab;
  return esq_eps_prime10*zeta_a*zeta_b/rab;

}
double Uyuk(double rab)
{
#if 0
  printf("qui esq_eps10=%.15G zeta_a=%f exp(-kD*rab)=%.15G rab=%.15G\n", esq_eps10, zeta_a, exp(-kD*rab), rab);
  printf("Uyuk=%.15G\n",esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab );
#endif
  
  return yuk_corr_fact*esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab; 
} 
#ifdef PARALLEL
double Uyuk_arr(double rab)
{
#if 0
  printf("qui esq_eps10=%.15G zeta_a=%f exp(-kD*rab)=%.15G rab=%.15G\n", esq_eps10, zeta_a, exp(-kD*rab), rab);
  printf("Uyuk=%.15G\n",esq_eps10*zeta_a*zeta_b*exp(-kD*rab)/rab );
#endif
  return yuk_corr_fact*esq_eps10*zeta_a*zeta_b/rab; 
}
double calc_yukawa_arr(int i, int j, double distsq, int *kks)
{
  double ret, rab0, rab, sigab;
  int kk, k1, k2;
  rab = sqrt(distsq);
  sigab = CHROMs[0][i].rad + CHROMs[1][j].rad;
  rab0 = sigab + delta_rab0; /* we are using Angstrom units here (see pag. S2 in the SI of Frezza Soft Matter 2011) */ 
  *kks=-2;
  if (distsq > maxyukcutkDsq)
    return 0.0;
  for (kk=0; kk < num_kD; kk++)
    {    
      if (rab < rab0)
	{
	  *kks=-1;
	  //printf("interp=%.15G\n",  Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab));
#ifdef NO_INTERP
	  return Ucoul(rab);
#else
	  return Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab);
#endif
#if 0
	  if (isnan(ret))
	    {
	      printf("a=%f rab=%f boh ret=%f Ucoul(rab)=%f Uyuk(rab)=%f\n", zeta_a, rab, ret, Ucoul(rab), Uyuk(rab));
	      exit(-1);
	    }
	  return ret;
#endif    
	}
      /* we set a cutoff for electrostatic interactions */
      else if (rab < yukcut*kD_sorted[kk].invkD)
	{
	  //printf("Yuk=%.15G\n", Uyuk(rab));
	  *kks = kk; 
	  return Uyuk_arr(rab);
	} 
    }
}
#endif
double calc_yukawa(int i, int j, double distsq)
{
  double ret, rab0, rab, sigab;
  rab = sqrt(distsq);
  sigab = CHROMs[0][i].rad + CHROMs[1][j].rad;
  if (sigab==0.0)
    delta_rab0 = 0.0;

  rab0 = sigab + delta_rab0; /* we are using Angstrom units here (see pag. S2 in the SI of Frezza Soft Matter 2011) */ 
  
  if (rab < rab0)
    {
      //printf("interp=%.15G\n",  Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab));
#ifdef NO_INTERP
      return Ucoul(rab);
#else
      return Ucoul(sigab) + (rab-sigab)*(Uyuk(rab0) - Ucoul(sigab))/(rab0-sigab);
#endif
#if 0
      if (isnan(ret))
	{
	  printf("a=%f rab=%f boh ret=%f Ucoul(rab)=%f Uyuk(rab)=%f\n", zeta_a, rab, ret, Ucoul(rab), Uyuk(rab));
	  exit(-1);
	}
      return ret;
#endif    
    }
  /* we set a cutoff for electrostatic interactions */
  else if (rab < yukcutkD)
    {
      //printf("Yuk=%.15G\n", Uyuk(rab));
      return Uyuk(rab);
    } 
  else return 0.0;
}
#endif
struct boxS {
  double R[3][3];
  double x[3];
  double sax[3];
};

typedef struct amyloid {
  double nD; /* spessore della fibra in protofilamenti*/
  double nL; /* numero di box che costituiscono la fibra */
  double nP; /* persistence length in numero di boxes */
  double ribthick; /* profondità box (lungo direzione ortogonale al piano dell'amyloid ribbon) */
  double npitch;
  double Lbox;
  double Dproto; /* diametro del protofilamento  Dproto*nD è la larghezza del box */
  struct boxS *boxes;
  struct boxS *boxesLab;
  double R[3][3];
  double rcm[3];
  double boxsax[3];
} amyloidS;

amyloidS amyloids[2];

 
void build_amyloid(int nL)
{
  double npitch, dth;
  int jj, kk, i;
  double xcm[3];

  amyloids[0].nL=nL;
  amyloids[1].nL=nL;
  amyloids[0].boxes = malloc(sizeof(struct boxS)*amyloids[0].nL);
  amyloids[1].boxes = malloc(sizeof(struct boxS)*amyloids[1].nL);
  amyloids[0].boxesLab = malloc(sizeof(struct boxS)*amyloids[0].nL);
  amyloids[1].boxesLab = malloc(sizeof(struct boxS)*amyloids[1].nL);

  for (i=0; i < 2; i++)
    {
      amyloids[i].nD=6;
      amyloids[i].nP=652;
      amyloids[i].npitch=35;
      amyloids[i].Lbox=2.3;/* altezza lungo l'asse dell'elica del box in nm (qui considero 4 ILQINS stacked 
			      per ridurre il numero di box totali */ 
      amyloids[i].Dproto=2.0;
      amyloids[i].ribthick = 2.0; /* per ora lo spessore del foglietto è di un solo ILQINS quindi circa 2 nm */
      amyloids[i].boxsax[0] = amyloids[i].Lbox*(amyloids[i].nL+1)/2.0;
      amyloids[i].boxsax[1] = (amyloids[i].nD+1)*amyloids[i].Dproto/2.0;
      amyloids[i].boxsax[2] = (amyloids[i].nD+1)*amyloids[i].Dproto/2.0;
    
      dth=2.0*M_PI/amyloids[i].npitch;
      // length_eucl=OprogStatus.npitch*OprogStatus.pitch;
      for (jj=0; jj < amyloids[i].nL; jj++)
	{
	  amyloids[i].boxes[jj].x[0]=((double)jj)*amyloids[i].Lbox;
	  //printf("x=%f\n", amyloids[i].boxes[jj].x[0]);
	  amyloids[i].boxes[jj].x[1]=0.0;
	  amyloids[i].boxes[jj].x[2]=0.0;
	  printf("rcm %f %f %f\n", amyloids[i].boxes[jj].x[0],amyloids[i].boxes[jj].x[1],amyloids[i].boxes[jj].x[2]);
	  amyloids[i].boxes[jj].sax[0] = amyloids[i].Lbox/2.0;
	  amyloids[i].boxes[jj].sax[1] = amyloids[i].nD*amyloids[i].Dproto/2.0;
	  amyloids[i].boxes[jj].sax[2] = amyloids[i].ribthick/2.0;
	  /* ...and now set orientation */
	  amyloids[i].boxes[jj].R[0][0]=1.0;
	  amyloids[i].boxes[jj].R[0][1]=0.0;
	  amyloids[i].boxes[jj].R[0][2]=0.0;
	  amyloids[i].boxes[jj].R[1][0]=0.0;
	  amyloids[i].boxes[jj].R[1][1]=cos(jj*dth);
	  amyloids[i].boxes[jj].R[1][2]=-sin(jj*dth);
	  amyloids[i].boxes[jj].R[2][0]=0.0;
	  amyloids[i].boxes[jj].R[2][1]=sin(jj*dth);
	  amyloids[i].boxes[jj].R[2][2]=cos(jj*dth);
	}
      xcm[0]=xcm[1]=xcm[2]=0.0;
      for (jj=0; jj < amyloids[i].nL; jj++)
	{
	  for (kk=0; kk < 3; kk++)
	    xcm[kk] += amyloids[i].boxes[jj].x[kk]; 
	} 

      for (kk=0; kk < 3; kk++)
	xcm[kk] /= amyloids[i].nL;
      for (jj=0; jj < amyloids[i].nL; jj++)
	{
	  for (kk=0; kk < 3; kk++)
	    amyloids[i].boxes[jj].x[kk] -= xcm[kk];
    	}   
    }
}


#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
#define SYMMETRY
#ifdef QUASIMC
  double sv[10];
  int nsv;
#endif
struct CHROMallStr {
  double R[3][3];
  double rcm[3];
  double sax[3];
  double boxsax[3];
} CHROMall[2];
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
int check_convergence(double Told[3], double Tnew[3])
{
  double test=0.0, temp;
  int i;
  for (i=0;i<3;i++) 
    {
      temp=(fabs(Tnew[i]-Told[i]))/FMAX(fabs(Tnew[i]),1.0); 
      //temp=(fabs(x[i]-xold[i]))/fabs(x[i]); 
      if (temp > test) 
	test=temp; 
    }
  if (test < 1.0E-13)
    {
      //printf("convergence reached! test=%.15G\n", test);
      return 1;
    }
  else 
    return 0;
}

void versor_to_R(double ox, double oy, double oz, double R[3][3], double gamma);
void versor_to_R_sym(double ox, double oy, double oz, double R[3][3]);


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
double calcDistBox(double rcmA[3], double saxA[3], double RA[3][3],double rcmB[3], double saxB[3], double RB[3][3])
{
  double RR, R0, R1, cij[3][3], fabscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k, k1, k2, existsParallelPair = 0;
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  //return -1;
  for (k=0; k < 3; k++)
    {
      rA[k] = rcmA[k];
      rB[k] = rcmB[k];
      EA[k] = saxA[k];
      EB[k] = saxB[k];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  AA[k1][k2] = RA[k1][k2];
	  BB[k1][k2] = RB[k1][k2];
	}
    	DD[k1] = rA[k1] - rB[k1];
    }
  /* axis C0+s*A0 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[0][k1] =  scalProd(AA[0], BB[k1]);
      fabscij[0][k1] = fabs(cij[0][k1]);
      if ( fabscij[0][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[0] = scalProd(AA[0],DD);
  RR = fabs(AD[0]);
  R1 = EB[0]*fabscij[0][0]+EB[1]*fabscij[0][1]+EB[2]*fabscij[0][2];
  R01 = EA[0] + R1;
  if ( RR > R01 )
    return 1.0; /* non si intersecano */
  /* axis C0+s*A1 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[1][k1] = scalProd(AA[1],BB[k1]);
      fabscij[1][k1] = fabs(cij[1][k1]);
      if ( fabscij[1][k1] == 1.0  )
	existsParallelPair = 1;
    }
  AD[1] = scalProd(AA[1],DD);
  RR = fabs(AD[1]);
  R1 = EB[0]*fabscij[1][0]+EB[1]*fabscij[1][1]+EB[2]*fabscij[1][2];
  R01 = EA[1] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*A2 */
  for (k1= 0; k1 < 3; k1++)
    {
      cij[2][k1] = scalProd(AA[2], BB[k1]);
      fabscij[2][k1] = fabs(cij[2][k1]);
      if ( fabscij[2][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[2] = scalProd(AA[2],DD);
  RR = fabs(AD[2]);
  R1 = EB[0]*fabscij[2][0]+EB[1]*fabscij[2][1]+EB[2]*fabscij[2][2];
  R01 = EA[2] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*B0 */
  RR = fabs(scalProd(BB[0],DD));
  R0 = EA[0]*fabscij[0][0]+EA[1]*fabscij[1][0]+EA[2]*fabscij[2][0];
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B1 */
  RR = fabs(scalProd(BB[1],DD));
  R0 = EA[0]*fabscij[0][1]+EA[1]*fabscij[1][1]+EA[2]*fabscij[2][1];
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;
  
  /* axis C0+s*B2 */
  RR = fabs(scalProd(BB[2],DD));
  R0 = EA[0]*fabscij[0][2]+EA[1]*fabscij[1][2]+EA[2]*fabscij[2][2];
  R01 = R0 + EB[2];
  if ( RR > R01 )
    return 1.0;

  /* At least one pair of box axes was parallel, therefore the separation is
   * effectively in 2D, i.e. checking the "edge" normals is sufficient for
   * the separation of the boxes. 
   */
  if ( existsParallelPair )
    return -1.0;

  /* axis C0+s*A0xB0 */
  RR = fabs(AD[2]*cij[1][0]-AD[1]*cij[2][0]);
  R0 = EA[1]*fabscij[2][0] + EA[2]*fabscij[1][0];
  R1 = EB[1]*fabscij[0][2] + EB[2]*fabscij[0][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB1 */
  RR = fabs(AD[2]*cij[1][1]-AD[1]*cij[2][1]);
  R0 = EA[1]*fabscij[2][1] + EA[2]*fabscij[1][1];
  R1 = EB[0]*fabscij[0][2] + EB[2]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB2 */
  RR = fabs(AD[2]*cij[1][2]-AD[1]*cij[2][2]);
  R0 = EA[1]*fabscij[2][2] + EA[2]*fabscij[1][2];
  R1 = EB[0]*fabscij[0][1] + EB[1]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB0 */
  RR = fabs(AD[0]*cij[2][0]-AD[2]*cij[0][0]);
  R0 = EA[0]*fabscij[2][0] + EA[2]*fabscij[0][0];
  R1 = EB[1]*fabscij[1][2] + EB[2]*fabscij[1][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB1 */
  RR = fabs(AD[0]*cij[2][1]-AD[2]*cij[0][1]);
  R0 = EA[0]*fabscij[2][1] + EA[2]*fabscij[0][1];
  R1 = EB[0]*fabscij[1][2] + EB[2]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB2 */
  RR = fabs(AD[0]*cij[2][2]-AD[2]*cij[0][2]);
  R0 = EA[0]*fabscij[2][2] + EA[2]*fabscij[0][2];
  R1 = EB[0]*fabscij[1][1] + EB[1]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB0 */
  RR = fabs(AD[1]*cij[0][0]-AD[0]*cij[1][0]);
  R0 = EA[0]*fabscij[1][0] + EA[1]*fabscij[0][0];
  R1 = EB[1]*fabscij[2][2] + EB[2]*fabscij[2][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB1 */
  RR = fabs(AD[1]*cij[0][1]-AD[0]*cij[1][1]);
  R0 = EA[0]*fabscij[1][1] + EA[1]*fabscij[0][1];
  R1 = EB[0]*fabscij[2][2] + EB[2]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB2 */
  RR = fabs(AD[1]*cij[0][2]-AD[0]*cij[1][2]);
  R0 = EA[0]*fabscij[1][2] + EA[1]*fabscij[0][2];
  R1 = EB[0]*fabscij[2][1] + EB[1]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  return -1.0;
}

/* cylinder overlap routines here */
double diamHC=2.0, lengthHC=2.0;
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
void saveAmyloidPovray(amyloidS amy, char *fname)
{
  double L=10;
  FILE *f, *fs;
  int i, jj, kk, k1, k2, k3;
  double R[3][3], ux[3], uy[3], uz[3], x1[3], x2[3], xbox[3];
  f = fopen(fname, "w+");
  fs=f;  
  /* povray preamble */
  fprintf(fs, "#declare RAD = off;\n");
  fprintf(fs, "#declare DIFF=0.475;\n");
  fprintf(fs, "#declare AMB=0.475;\n");
  fprintf(fs, "#declare REFL=0.025;\n");
  fprintf(fs, "#declare TRAN=0.0;\n");
  fprintf(fs, "#declare ROUGH=0.001;\n");
  fprintf(fs,"#declare SPEC=0.0;\n");
  fprintf(fs,"#declare PHONG=1;\n");
  fprintf(fs,"#declare PHONG_SIZE=10;\n");
  fprintf(fs,"#include \"colors.inc\"\n");
  fprintf(fs,"#include \"textures.inc\"\n");
  fprintf(fs,"#include \"stones.inc\"\n");
  fprintf(fs,"#include \"golds.inc\"\n");
  fprintf(fs,"global_settings {\n");
  fprintf(fs,"#if (RAD)\n");
  fprintf(fs,"radiosity {\n");
  fprintf(fs,"pretrace_start 0.08\n");
  fprintf(fs,"pretrace_end 0.01\n");
  fprintf(fs,"count 500\n");
  fprintf(fs,"nearest_count 10\n");
  fprintf(fs,"error_bound 0.02\n");
  fprintf(fs,"recursion_limit 1\n");
  fprintf(fs,"low_error_factor 0.2\n");
  fprintf(fs,"gray_threshold 0.0\n");	
  fprintf(fs,"minimum_reuse 0.015\n");
  fprintf(fs, "brightness 0.6\n");
  fprintf(fs, "adc_bailout 0.01/2\n");
  fprintf(fs, "}\n");
  fprintf(fs, "#end\n");
  fprintf(fs, "}\n");
  fprintf(fs, "background{White}\n");
  fprintf(fs, "camera {\n");	
  fprintf(fs, "angle %f\n", 30.0);
  fprintf(fs, "location <%f,%f,%f>\n", 2.5*L, 2.0*L, -5.2*L);
  fprintf(fs, "look_at <0,0,0>\n");
  fprintf(fs, "//focal_point < 1, 1, -6> // pink sphere in focus\n");
  fprintf(fs, "//aperture 0.4\n");
  fprintf(fs, "//blur_samples 20\n");
  fprintf(fs, "}\n");
  fprintf(fs, "#if(1)\n");
  fprintf(fs, "// spotlight light source\n");
  fprintf(fs, "light_source {\n");
  fprintf(fs, "//<0, 10, -3>\n");
  fprintf(fs, "<%f, %f, %f>\n", 10.0*L, 10.0*L, -10.0*L);
  fprintf(fs, "color White\n");
  fprintf(fs, "spotlight\n");
  fprintf(fs, "radius 15\n");
  fprintf(fs, "falloff 20\n");
  fprintf(fs, "tightness 10\n");
  fprintf(fs, "point_at <0, 0, 0>\n");
  fprintf(fs, "}\n");
  fprintf(fs, "#else\n");
  fprintf(fs, "//point light source\n");
  fprintf(fs, "light_source { <10, 20, -10> color White }\n");
  fprintf(fs, "#end\n");	
#if 0
  fprintf(fs,"#include \"shapes.inc\"\n");
  fprintf(fs,"object {\n");
  // Wire_Box(A, B, WireRadius, Merge)
#ifdef MD_LXYZ
  fprintf(fs,"Wire_Box(<%f,%f,%f>,<%f,%f,%f>, %f, 0",-L[0]*0.5,-L[1]*0.5,-L[2]*0.5,L[0]*0.5,L[1]*0.5,L[2]*0.5,
	  OprogStatus.povboxthickness);
#else
  fprintf(fs,"Wire_Box(<%f,%f,%f>,<%f,%f,%f>, %f, 0",-L*0.5,-L*0.5,-L*0.5,L*0.5,L*0.5,L*0.5,
	  OprogStatus.povboxthickness);
#endif
  fprintf(fs,"texture{");
  fprintf(fs,"pigment{ color rgb<0,0,0>}");
  // fprintf(fs,"finish {  phong 1}\n");
  fprintf(fs,"}\n"); // end of texture
  fprintf(fs,"scale<1,1,1>\n");
  fprintf(fs,"  rotate<0,0,0>\n");
  fprintf(fs,"  translate<0,0,0>");
  fprintf(fs," } //---------------------------------\n");
#endif
  for (jj=0; jj < amy.nL; jj++)
    {
      body2lab(amy.boxes[jj].x,xbox,amy.rcm,amy.R);

      printf(">>>rcm %f %f %f\n", amy.boxes[jj].x[0],amy.boxes[jj].x[1],amy.boxes[jj].x[2]);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      R[k1][k2] = 0.0;//R[i][k1][k2];
	      for (k3 = 0; k3 < 3; k3++)
		{
		  /* matrix multiplication: riga * colonna */
		  R[k1][k2] += amy.boxes[jj].R[k1][k3]*amy.R[k3][k2];
		}
	    }
	}

      //print_matrix(R);
      fprintf(f,"box\n");	
      fprintf(f,"{\n");
      for (kk=0; kk < 3; kk++)
	{
  	  ux[kk]=R[kk][0];
	  uy[kk]=R[kk][1];
	  uz[kk]=R[kk][2];
	}
      printf("ux=%f %f %f\n", ux[0], ux[1], ux[2]);
      printf("uy=%f %f %f\n", uy[0], uy[1], uy[2]);
      printf("uz=%f %f %f\n", uz[0], uz[1], uz[2]);
      printf("xbox %f %f %f\n", xbox[0], xbox[1], xbox[2]);
      for (kk=0; kk < 3; kk++)
	{
	  x1[kk] = xbox[kk]+ux[kk]*amy.boxes[jj].sax[0]+uy[kk]*amy.boxes[jj].sax[1]+uz[kk]*amy.boxes[jj].sax[2];
	  x2[kk] = xbox[kk]-ux[kk]*amy.boxes[jj].sax[0]-uy[kk]*amy.boxes[jj].sax[1]-uz[kk]*amy.boxes[jj].sax[2];
	}
      fprintf(f,"<%.15G,%.15G,%.15G>, <%.15G,%.15G,%.15G>\n",x1[0],x1[1],x1[2], x2[0],x2[1],x2[2]);
      fprintf(f," pigment { color rgb<0.0,0.9,0.2> transmit 0.0 }\n");
      //pigment { Green transmit 0.7 }
      fprintf(f,"rotate <0,0,90>\n");
      //fprintf(fs, "matrix <%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,0,0,0>\n",
	//      R[0][0],R[0][1],R[0][2],R[1][0],R[1][1],R[1][2],R[2][0],R[2][1],R[2][2]);
      fprintf(f,"translate <%.15G, %.15G, %.15G>", xbox[0], xbox[1], xbox[2]);
      fprintf(f,"finish { phong PHONG phong_size PHONG_SIZE reflection REFL\n"); 
      fprintf(f,"ambient AMB diffuse DIFF\n");
      fprintf(f,"specular SPEC roughness ROUGH}\n");
      fprintf(f, "}\n");
    }
  fclose(f);
}
double calcDistAmyloid()
{
  int iA, iB, k1, k2, k3, jj;
  
  /* calc orientation matrix and position of boxes in the laboratory frame */

  for (jj=0; jj < amyloids[0].nL; jj++)
    {
      body2lab(amyloids[0].boxes[jj].x,amyloids[0].boxesLab[jj].x,amyloids[0].rcm,amyloids[0].R);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      amyloids[0].boxesLab[jj].R[k1][k2] = 0.0;//R[i][k1][k2];
	      for (k3 = 0; k3 < 3; k3++)
		{
		  /* matrix multiplication: riga * colonna */
		  amyloids[0].boxesLab[jj].R[k1][k2] += amyloids[0].boxes[jj].R[k1][k3]*amyloids[0].R[k3][k2];
		}  
	    }
	}
    }
  for (jj=0; jj < amyloids[1].nL; jj++)
    {
      body2lab(amyloids[1].boxes[jj].x,amyloids[1].boxesLab[jj].x,amyloids[1].rcm,amyloids[1].R);
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      amyloids[1].boxesLab[jj].R[k1][k2] = 0.0;//R[i][k1][k2];
	      for (k3 = 0; k3 < 3; k3++)
		{
		  /* matrix multiplication: riga * colonna */
		  amyloids[1].boxesLab[jj].R[k1][k2] += amyloids[1].boxes[jj].R[k1][k3]*amyloids[1].R[k3][k2];
		}  
	    }
	}
    }

  for (iA=0; iA < amyloids[0].nL; iA++)
    for (iB=0; iB < amyloids[0].nL; iB++)
      {

	if (calcDistBox(amyloids[0].boxesLab[iA].x,amyloids[0].boxes[iA].sax,amyloids[0].boxesLab[iA].R,
		        amyloids[1].boxesLab[iB].x,amyloids[1].boxes[iB].sax,amyloids[1].boxesLab[iB].R) < 0.0)
	  return -1;
      }
  return 1;
}


/* ------------------------------- */
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
#define MAXBIT 30 
#define MAXDIM 6
void sobseq(int *n, double x[])
/*When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[1..n] the next values from n of these sequences. (n must not be changed between initializations.)*/
{
  int j,k,l;
  unsigned long i,im,ipp;
  static double fac;
  static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
  static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
  static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4}; 
  static unsigned long iv[MAXDIM*MAXBIT+1]={0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
  if (*n < 0) 
    { 
      /*Initialize, don’t return a vector. */
      for (k=1;k<=MAXDIM;k++) ix[k]=0;
      in=0;
      if (iv[1] != 1) return;
      fac=1.0/(1L << MAXBIT);
      for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) 
	iu[j] = &iv[k];/* To allow both 1D and 2D addressing.*/
      for (k=1;k<=MAXDIM;k++) 
	{
	  for (j=1;j<=mdeg[k];j++)
	    iu[j][k] <<= (MAXBIT-j); /*Stored values only require normalization.*/
	  for (j=mdeg[k]+1;j<=MAXBIT;j++) 
	    {
	      ipp=ip[k]; i=iu[j-mdeg[k]][k];
	      i ^= (i >> mdeg[k]);
	      for (l=mdeg[k]-1;l>=1;l--) 
		{
		  if (ipp & 1) i ^= iu[j-l][k];
		  ipp >>= 1; 
		}
	      iu[j][k]=i;
	    }
	}
    } 
  else 
    {
      im=in++;
      for (j=1;j<=MAXBIT;j++) {
	if (!(im & 1)) break;
	im >>= 1; }
      if (j > MAXBIT) {
	printf("MAXBIT too small in sobseq");
	exit(-1);
      } 
      im=(j-1)*MAXDIM;
      for (k=1;k<=IMIN(*n,MAXDIM);k++)
	{
	  ix[k] ^= iv[im+k]; 
	  x[k]=ix[k]*fac;
	}
      /*XOR the appropriate direction num- ber into each component of the vector and convert to a floating number.
       */
    }
}

//#define NO_INTERP

char dummy1[32], dummy2[32], atname[32], nbname[8];
int nat, atnum, nbnum, len;
long long int tot_trials, tt=0, ttini=0;
double gamma1, gamma2, L, rx, ry, rz, alpha, dfons_sinth_max, fons_sinth_max;
const double thetapts=100000;
int type;


//#define ALBERTA
char fn[1024];
#define MC_BENT_DBLCYL

/*
                pdb        radius (angstrom)
    sugar       Xe          3.5        
    phosphate   B           3.0         
    base        Se          4.0    
*/

#ifdef AMYLOID_ELEC
/* apply a random rotation around the supplied axis because 
   bent cylinders do not have azimuthal symmetry */
double thetaGlobalBondangle;
void add_rotation_around_axis(double ox, double oy, double oz, double Rin[3][3], double Rout[3][3], double gamma)
{
  double theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random rotation angle between 0 and 2*pi*/
  //theta = 4.0*acos(0.0)*drand48();
  theta = gamma;
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
//	  Ro[k1][k2] += M[k1][k3]*Rin[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     Rout[k1][k2] = Ro[k1][k2]; 
}
#endif
void versor_to_R_sym(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
  /* first row vector (note that cylinder symmetry axis which is x) */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = 0.0; u[1] = 0.0; u[2] = 1.0;
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

void versor_to_R(double ox, double oy, double oz, double R[3][3], double gamma)
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
#ifdef AMYLOID_ELEC
  double Rout[3][3];
  int k1, k2;
#endif
  /* first row vector (note that cylinder symmetry axis which is x) */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = 0.0; u[1] = 0.0; u[2] = 1.0;
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
#ifdef AMYLOID_ELEC
  /* add a random rotation around the axis (ox, oy, oz) */
  add_rotation_around_axis(ox, oy, oz, R, Rout, gamma);
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
double RMDNA[2][3][3];
void place_AMYLOID(double x, double y, double z, double ux, double uy, double uz, int which, double gamma)
{
  double xp[3], rO[3], xl[3];
  double R[3][3];
  int k1, k2;
  FILE *f;
  int i; 
  rO[0] = x;
  rO[1] = y;
  rO[2] = z;
  /* build R here from the orientation (ux,uy,uz) */
  versor_to_R(ux, uy, uz, R, gamma);

  /* ============ */
  for (k1=0; k1 < 3; k1++)
    {
      amyloids[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	amyloids[which].R[k1][k2] = R[k1][k2];
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
#if 0
      f0 = 1.01*fons(0.0,alpha);
#else      
      f0 = 1.01*fons_sinth_max;
#endif
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
void orient(double *omx, double *omy, double* omz)
{
  int i;
  //double inert;                 /* momentum of inertia of the molecule */
  //double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, osq, norm;
  
  //Mtot = m; /* total mass of molecule */

  //inert = I; /* momentum of inertia */
 
  //mean = 3.0*temp / inert;

  xisq = 1.0;

  while (xisq >= 1.0)
    {
      xi1  = ranf_vb() * 2.0 - 1.0;
      xi2  = ranf_vb() * 2.0 - 1.0;
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

  *omx = ox;
  *omy = oy;
  *omz = oz; 
  //distro[(int) (acos(oz)/(pi/1000.0))] += 1.0;

#if 0
  /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/

  osq   = - mean * log(ranf());
  o     = sqrt(fabs(osq));
  ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      *wx = ox;
      *wy = oy;
      *wz = oz;
#endif 
}

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
  return sinh(alpha*cos(theta))*alpha*alpha/(4.0*pi*sinh(alpha));
}

/* return an angle theta sampled from an Onsager angular distribution */
double theta_donsager(double alpha, int domain)
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
      //f0 = 1.01*alpha*fons(0.0,alpha);
      f0 = 1.01*dfons_sinth_max;
    }
  do 
    {
      /* uniform theta between 0 and pi/2 (domain=0) or pi/2 and pi (domain=1) */
      if (domain==0)
	theta = ranf_vb()*pi/2.;
      else
	theta = (pi/2.)*(1.+ranf_vb());
      /* uniform y between 0 and 1 (note that sin(theta) <= 1 for 0 < theta < pi)*/
      y = f0*ranf_vb();
      f = fabs(sin(theta)*dfons(theta,alpha));
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}

//extern const int nfons;
void orient_donsager(double *omx, double *omy, double* omz, double alpha, int domain)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_donsager(alpha, domain);
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
double estimate_maximum_dfons(double alpha)
{
  double th, dth, maxval, m;
  int i;
  dth=2.0*(acos(0.0))/((double)thetapts);
  th=0.0;
  for (i=0; i < thetapts; i++)
    {
      m=sin(th)*dfons(th,alpha);
      if (i==0 || maxval < m)
	maxval = m;
      th += dth;
      //printf("%f %.15G\n", th, sin(th)*dfons(th, alpha));
    }
  // printf("maxval=%f\n", maxval);
  return maxval;
}

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
double max2(double a, double b)
{
  if (a > b)
    return a;
  else 
    return b;
}
#if 0
double epsrtbl[21][2]={{0,87.740},{5,85.763},{10,83.832},{15,81.946},{20,80.103},{25,78.304},{30,76.546},{35,74.828},{40,73.151},{45,71.512},{50,69.910},{55,68.345},{60,66.815},{65,65.319},{70,63.857},{75,62.427},{80,61.027},{85,59.659},{90,58.319},{95,57.007},{100,55.720}};
double epsr(double T)
{
  int i;
  i = 0;
  T = T - 273.15;
  if (T <= 0.0 || T >= 100.0)
    {
      printf("Temperature must be between 0 and 100!\n");
      exit(-1);
    }
  for (i=0; i < 21; i++)
    {
      if (T > epsrtbl[i][0])
	{
	  /* linear interpolation */
	  return epsrtbl[i][1]+(T-epsrtbl[i][0])*(epsrtbl[i+1][1]-epsrtbl[i][1])/(epsrtbl[i+1][0]-epsrtbl[i][0]);
	}
    }
}
#endif
int main(int argc, char**argv)
{
  double aaa;
  int nL;
#ifdef QUASIMC
#ifdef USEGSL
  gsl_qrng *qsob;
#endif
#endif
#ifdef AMYLOID_ELEC
  double uel, beta;
  int interact;
#endif
#ifdef PARALLEL
  FILE *fp;
  double sigab, rab0, rab0sq, uelcontrib, tempfact;
  int k1, k2, kk;
#endif
  double Lx, Ly, Lz;
  FILE *fin, *fout, *f, *fread;
  int ncontrib, cc, k, i, j, overlap, contrib, cont=0, nfrarg, kk;
  long long int fileoutits, outits;
  char fnin[1024],fnout[256];
  double dummydbl, segno, u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, vexclel=0.0, factor, dth, th;
  /* syntax:  CG-DNA-k2K22 <pdb file> <CHROM length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */
  if (argc < 7)
    {
#ifdef AMYLOID_ELEC
      printf("syntax:  AMYLOID-K11K22K33 <AMYLOID diam> <AMYLOID len> <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2> <fileoutits> [outits] [Temperature (in K)] [CHROM concentration in mg/ml] [yukawa cutoff in units of 1/kD] [epsr_prime (1.0-3.0, default=2 ] [delta_rab0 (default=2) ]\n");
#else
      printf("syntax:  AMYLOID-K11K22K33 <AMYLOID diam> <AMYLOID len> <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2> <fileoutits> [outits]\n");
#endif
      exit(1);
    }
  diamHC=atof(argv[1]);
  nL=atoi(argv[2]);
  alpha = atof(argv[3]);
  tot_trials=atoll(argv[4]);
  type = atoi(argv[5]);
  fileoutits = atoll(argv[6]);
 
  printf("length=%d\n", nL);
  if (argc == 7)
    outits=100*fileoutits;
  else
    outits = atoll(argv[7]);
#ifdef AMYLOID_ELEC
#ifdef PARALLEL
  if (argc <= 8)
    {
      numtemps=1;
      beta = 1.0;
    }
  else
    {
      if (sscanf(argv[8], "%lf", &dummydbl) < 1)
	{
	  fp=fopen(argv[8],"r");
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      cc++;
	    }
	  beta_arr = malloc(sizeof(double)*cc);
	  rewind(fp);
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      beta_arr[cc] = 1.0/dummydbl;
	      cc++;
	    }
	  fclose(fp);
	  numtemps=cc;
	}
      else
	{
	  numtemps = 1;
	  beta = 1.0/atof(argv[8]);
	} 
    }
  if (argc <= 9)
    {
      numconcs = 1;
      cchrom = 600;
    }
  else
    {
      if (!sscanf(argv[9], "%lf", &dummydbl))
	{
	  fp=fopen(argv[9],"r");
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      cc++;
	    }
	  cchrom_arr = malloc(sizeof(double)*cc);
	  rewind(fp);
	  numconcs=cc;
	  cc=0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf ", &dummydbl);
	      cchrom_arr[cc] = dummydbl;
	      cc++;
	    }
	  fclose(fp);
	}
      else
	{
	  numconcs = 1;
	  cchrom = atof(argv[9]);
	}
    }
#else
  if (argc <= 8)
    beta = 1.0;
  else  
    beta = 1.0/atof(argv[8]);

  if (argc <= 9)
    cchrom =  600; /* mg/ml */
  else
    cchrom = atof(argv[9]);
#endif
  if (argc <= 10)
    yukcut = 2.0;
  else 
    yukcut = atof(argv[10]);

  if (argc <= 11)
    epsr_prime = 2.0;
  else
    epsr_prime = atof(argv[11]);
  if (argc <= 12)
    delta_rab0 = 0.2;
  else
    delta_rab0 = atof(argv[12]);

#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      esq_eps = Sqr(qel)/(4.0*M_PI*eps0)/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
    }
  else
    {
      esq_eps = Sqr(qel)/(4.0*M_PI*eps0*epsr(1.0/beta))/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
    }
#else
  esq_eps = Sqr(qel)/(4.0*M_PI*eps0*epsr(1.0/beta))/kB; /* epsilon_r per l'acqua a 20°C vale 80.1 */
#endif
  esq_eps10 = esq_eps*1E9;
  esq_eps_prime = Sqr(qel)/(4.0*M_PI*eps0*epsr_prime)/kB;
  esq_eps_prime10 = esq_eps_prime*1E9;
#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      ximanning = esq_eps/bmann;
      deltamann = 1.0/ximanning;
      zeta_a = deltamann;
      zeta_b = deltamann;
      //printf("zeta_a=%f zeta_b:%f\n", zeta_a, zeta_b);
    }
  else
    {
      ximanning = esq_eps*beta/bmann;
      deltamann = 1.0/ximanning;
      zeta_a = deltamann;
      zeta_b = deltamann;
    }
#else
  ximanning = esq_eps*beta/bmann;
  deltamann = 1.0/ximanning;
  zeta_a = deltamann;
  zeta_b = deltamann;
#endif
  /*
     rho_salt =2 csalt Nav 1000;
     rho_counter[cdna_]:=(2 cdna)/(660*Dalton);
     (* cdna in mg/ml e csalt Molare (=moli/litro), nu=numero di cariche per unità di lunghezza *)
     InvDebyeScrLen[T_,qdna_,cdna_,qsalt_,csalt_,\[Epsilon]rel_]:= Sqrt[qdna^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) ( 2cdna)/(660*Dalton)+qsalt^2 qel^2/(kB T \[Epsilon]0 \[Epsilon]rel ) 2 csalt Nav 1000 ];
     InvDebyeScrLen[300, 2, 200, 2, 1, 20]^-1*10^9

   */
  /* qdna è la carica rilasciata da ogni gruppo fosfato in soluzione (tipicamente=1) */
#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      kD_arr = malloc(sizeof(double*)*numtemps);
      yukcutkD_arr = malloc(sizeof(double*)*numtemps);
      yukcutkDsq_arr = malloc(sizeof(double*)*numtemps); 
      vexclel_arr = malloc(sizeof(double*)*numtemps); 
      uel_arr = malloc(sizeof(double*)*numtemps);
      yuk_corr_fact_arr = malloc(sizeof(double*)*numtemps); 
      for (k1=0; k1 < numtemps; k1++)
	{
	  kD_arr[k1] = malloc(sizeof(double)*numconcs);
	  yukcutkD_arr[k1] = malloc(sizeof(double)*numconcs);
	  yukcutkDsq_arr[k1] = malloc(sizeof(double)*numconcs);
	  vexclel_arr[k1] = malloc(sizeof(double)*numconcs);
	  uel_arr[k1] = malloc(sizeof(double)*numconcs);
	  yuk_corr_fact_arr[k1] = malloc(sizeof(double)*numconcs); 
	}
      for (k1 = 0; k1 < numtemps; k1++)
	for (k2 = 0; k2 < numconcs; k2++)
	  {
	    kD_arr[k1][k2] = sqrt((4.0*M_PI*esq_eps*(1.0/epsr(1.0/beta_arr[k1])))*beta_arr[k1]*(Sqr(qchrom)*2.0*epsr(1.0/beta_arr[k1])*(deltamann/beta_arr[k1])*cchrom_arr[k2]/CHROMmolmass/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E9;
	    printf("numtemps=%d numconcs=%d kD:%f beta_arr:%f cchrm_arr: %f\n", numtemps, numconcs, kD_arr[k1][k2], beta_arr[k1], cchrom_arr[k2]);
	    printf("esq_eps: %f qdna=%f deltamann=%f qsalt=%f csalt=%f\n", esq_eps, qchrom, deltamann, qsalt, csalt);
	    /* 6.0 Angstrom is the closest distance between phosphate charges */
	    yukcutkD_arr[k1][k2] = yukcut/kD_arr[k1][k2];
	    yukcutkDsq_arr[k1][k2] = Sqr(yukcutkD_arr[k1][k2]);	
	    yuk_corr_fact_arr[k1][k2] = exp(kD_arr[k1][k2]*6.0)/(1.0+kD_arr[k1][k2]*6.0);
	  }
      num_kD = numtemps*numconcs;
      kD_sorted = malloc(sizeof(struct kDsortS)*num_kD);
      cc=0;
      for (k1 = 0; k1 < numtemps; k1++)
	for (k2 = 0; k2 < numconcs; k2++)
	  {
	    kD_sorted[cc].invkD = 1.0/kD_arr[k1][k2];
	    kD_sorted[cc].k1 = k1;
	    kD_sorted[cc].k2 = k2;
	    cc++;
	  }
      qsort(kD_sorted, cc, sizeof(struct kDsortS), compare_func);
      yuk_corr_fact = 1.0;//exp((1.0/kD_sorted[cc-1].invkD)*6.0)/(1.0+(1.0/kD_sorted[cc-1].invkD)*6.0);
      maxyukcutkD = yukcut*kD_sorted[cc-1].invkD;
      maxyukcutkDsq = Sqr(yukcut*kD_sorted[cc-1].invkD);
      printf("min: %f max: %f maxyukcutkD=%f\n",kD_sorted[0].invkD,kD_sorted[cc-1].invkD, sqrt(maxyukcutkDsq));
    }
  else
    {
      kD = sqrt((4.0*M_PI*esq_eps)*beta*(Sqr(qchrom)*2.0*deltamann*cchrom/CHROMmolmass/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E9;
      /* 6.0 Angstrom is the closest distance between phosphate charges */
      yuk_corr_fact = 1.0;//exp(kD*6.0)/(1.0+kD*6.0);
      yukcutkD = yukcut/kD;
      yukcutkDsq = Sqr(yukcutkD);
    }
#else
  kD = sqrt((4.0*M_PI*esq_eps)*beta*(Sqr(qchrom)*2.0*deltamann*cchrom/CHROMmolmass/Dalton + Sqr(qsalt)*2.0*csalt*Nav*1000.))/1E9;
  /* 6.0 Angstrom is the closest distance between phosphate charges */
  yuk_corr_fact = 1.0;//exp(kD*6.0)/(1.0+kD*6.0);
  yukcutkD = yukcut/kD;
  yukcutkDsq = Sqr(yukcutkD);
#endif
#ifdef PARALLEL
  printf("epsr_prime=%f beta=%f deltamanning=%.15G kB=%.15G kD=%.15G (in nm^-1) esq_eps=%.15G esq_eps_prime=%.15G yukcut=%f\n", epsr_prime, beta, deltamann, kB, kD, esq_eps, esq_eps_prime, yukcut);
  printf("yukawa cutoff=%.15G yuk_corr_fact=%.15G\n", yukcutkD, yuk_corr_fact);
#else
  printf("epsr_prime=%f epsr=%f beta=%f deltamanning=%.15G kB=%.15G kD=%.15G (in nm^-1) esq_eps=%.15G esq_eps_prime=%.15G yukcut=%f\n", epsr_prime, epsr(1.0/beta), beta, deltamann, kB, kD, esq_eps, esq_eps_prime, yukcut);
  printf("yukawa cutoff=%.15G yuk_corr_fact=%.15G\n", yukcutkD, yuk_corr_fact);
#endif
#endif
  cont=0;
#ifdef AMYLOID_ELEC
  nfrarg = 14;
#else
  nfrarg = 9;
#endif
  if (argc == nfrarg)
    {
      cont=1;
      fread = fopen(argv[nfrarg-1], "r");
      printf("reading file = %s\n", argv[nfrarg-1]);
      while (!feof(fread))
	{
#ifdef AMYLOID_ELEC
	  fscanf(fread, "%lld %lf %lf %lf\n", &ttini, &dummydbl, &vexcl, &vexclel);
#else
	  fscanf(fread, "%lld %lf\n", &ttini, &vexcl);
#endif
	}
      fclose(fread);
#ifdef AMYLOID_ELEC
      printf("restarting tt=%lld vexcltot=%.15G vexcl=%.15G vexclel=%.15G\n", ttini, vexcl+vexclel, vexcl, vexclel);
#else
      printf("restarting tt=%lld vexcl=%.15G\n", ttini, vexcl);
#endif
    }
  else
    {
      vexcl = 0.0;
#ifdef AMYLOID_ELEC
      vexclel = 0.0;
#endif
      ttini = 0;
    }

  /* ELISA: ATOM    39   Xe   G A   14      -5.687  -8.995  37.824 */
  /* ALBERTA: HETATM    1  B            1     -1.067  10.243 -35.117 */
  /* len here is the number of dodecamers, where 70 is the number of atoms per dodecamers
     in our CG model */
  //nat = 1;//70*len;
  //L = 1.05*3.0*40*len; /* 4 nm is approximately the length of a 12 bp DNAD */ 
  /* read the CG structure */
  cc=0;
  
  build_amyloid(nL);

  L=1.05*2.0*sqrt(Sqr(amyloids[0].boxsax[0])+Sqr(amyloids[0].boxsax[1])+Sqr(amyloids[0].boxsax[2]))*3.0;
  for (kk=0; kk < 2; kk++)
    {
      amyloids[kk].R[0][0]=1.0;
      amyloids[kk].R[0][1]=0.0;
      amyloids[kk].R[0][2]=0.0;
      amyloids[kk].R[1][0]=0.0;
      amyloids[kk].R[1][1]=1.0;
      amyloids[kk].R[1][2]=0.0;
      amyloids[kk].R[2][0]=0.0;
      amyloids[kk].R[2][1]=0.0;
      amyloids[kk].R[2][2]=1.0;
      amyloids[kk].rcm[0] = 0.0;
      amyloids[kk].rcm[1] = 0.0;
      amyloids[kk].rcm[2] = 0.0;
    }
  saveAmyloidPovray(amyloids[0],"amyfibr.pov");
  exit(-1);
  printf("nat=%d L=%f alpha=%f I am going to calculate v%d and I will do %lld trials\n", nat, L, alpha, type, tot_trials);
  printf("box semiaxes=%f %f %f\n", amyloids[0].boxsax[0], amyloids[0].boxsax[1], amyloids[0].boxsax[2]);
#if 1
  srand48((int)time(NULL));
#else
  srand48(0);
#endif
  sprintf(fnout, "v%d.dat", type);
  factor=0.0;
  dfons_sinth_max=estimate_maximum_dfons(alpha);
  fons_sinth_max=dfons_sinth_max/alpha;
  printf("Estimated maximum of dfons is %f\n", dfons_sinth_max);
  //exit(-1);
  /* avendo diviso l'integrazione in theta negli intervalli [0,pi/2] e [pi/2,pi]
     il fattore si deve ottenere integrando fra 0 e pi/2 */
#if 0
  dth=acos(0.0)/((double)thetapts);
  th=0.0;
  //f = fopen("dfons.dat", "w+");
  for (i=0; i < thetapts; i++)
    {
      factor += 0.5*dth*sin(th)*(dfons(th, alpha) + dfons(th+dth,alpha));
      th += dth;
      //fprintf(f,"%f %.15G\n", th, dfons(th, alpha));
    };
  //fclose(f);
  factor= fabs(factor);
  factor *= 4.0*acos(0.0);
#else
  factor = alpha/2.0;
#endif
  printf("factor=%.15G\n", factor);
  if (cont)
    { 
#ifdef CHORM_ELEC
      if (type == 0)
	vexclel *= 1E3*((double)ttini)/(L*L*L);
      else if (type==1)
	vexclel *= 1E4*((double)ttini)/(L*L*L)/factor;
      else
	vexclel *= 1E5*((double)ttini)/(L*L*L)/Sqr(factor);
#endif    
      if (type == 0)
	vexcl *= 1E3*((double)ttini)/(L*L*L);
      else if (type==1)
	vexcl *= 1E4*((double)ttini)/(L*L*L)/factor;
      else
	vexcl *= 1E5*((double)ttini)/(L*L*L)/Sqr(factor);
    }
 // printf("VEXCL=%f VEXCLEL=%f\n", vexcl, vexclel);
  if (!cont)
    {
#ifdef PARALLEL
      if (numtemps > 1 || numconcs > 1)
	{
	  for (k1=0; k1 < numtemps; k1++)
    	    for (k2=0; k2 < numconcs; k2++)
	      {
		sprintf(fnout, "v%d_c%.0f_T%.0f.dat", type, cchrom_arr[k2], 1.0/beta_arr[k1]);
		fout = fopen(fnout, "w+");
      		fclose(fout);
	      }
	}
      else
	{    
	  fout = fopen(fnout, "w+");
	  fclose(fout);
	}
#else
      fout = fopen(fnout, "w+");
      fclose(fout);
#endif
    }  
  if (type==0)
    ncontrib=1;
  else if (type==1)
    ncontrib=2;
  else
#ifdef SYMMETRY
    ncontrib=2;
#else
  ncontrib=4;
#endif  
#if 0
  if (type==1)
    {
      Lx=Ly=max3(DNADall[0].sax[0],DNADall[0].sax[1],3.0*2.0*DNADall[0].sax[2]*sin(2.0*sqrt(2.0/alpha)));
      Lz=L;
      printf("Lx=%f Ly=%f Lz=%f\n", Lx, Ly, Lz);
    }
  else
    {}
#endif
  Lx=Ly=Lz=L;
#ifdef PARALLEL
  if (numtemps > 1 || numconcs > 1)
    {
      sigab = CHROMs[0][0].rad + CHROMs[1][0].rad;
      rab0 = sigab + delta_rab0; 
      rab0sq = Sqr(rab0);
    }
#endif
  printf("Lx=%f Ly=%f Lz=%f\n", Lx, Ly, Lz);
  printf("type=%d ncontrib=%d\n", type, ncontrib);
#ifdef QUASIMC
#ifdef USEGSL
 nsv = 5; 
 qsob = gsl_qrng_alloc (gsl_qrng_sobol, nsv);
#else
  /* initialization */
  nsv = -1;  
  sobseq(&nsv, sv);
  nsv = 5;
#endif
#endif
  for (tt=ttini+1; tt < tot_trials; tt++)
    {
      /* place first DNAD in the origin oriented according to the proper distribution */
      /* per la v2: contrib = 0 e 3 sono contributi con lo stesso segno ossia tra [0,Pi/2] e [0,Pi2] e tra [Pi/2,Pi] e [Pi/2,P]
	 invece 1 e 2 sono quelli misti che danno un segno meno
	 per la v1: contrib = 0 è il contributo positivo tra [0,Pi/2] e 1 quello negativo tra [Pi/2,Pi] */
      for (contrib=0; contrib < ncontrib; contrib++)
	{
	  if (type==0||type==1)
	    {
	      if (alpha > 0.0)
		orient_onsager(&u1x, &u1y, &u1z, alpha);
	      else
		{
		  u1x = 0.0;
		  u1y = 0.0;
		  u1z = 1.0;
		}
	    }
	  else
	    {
	      if (contrib==0||contrib==1)
		orient_donsager(&u1x, &u1y, &u1z, alpha, 0);
	      else 
		orient_donsager(&u1x, &u1y, &u1z, alpha, 1);
	    }
#ifdef QUASIMC
	  /* quasi-MC per rcmx, rcmy, rcmz e gamma2, gamma1 rimane random per 
	     poter fare run indipendenti */
#ifdef USEGSL
	  gsl_qrng_get (qsob, sv);
	  //printf("sv=%f %f %f %f %f\n",sv[1], sv[2], sv[3], sv[4], sv[5]);
	  rcmx = Lx*(sv[0]-0.5);
	  rcmy = Ly*(sv[1]-0.5);
	  rcmz = Lz*(sv[2]-0.5);
	  gamma1 = 2.0*M_PI*sv[3];
	  gamma2 = 2.0*M_PI*sv[4];
#else
	  sobseq(&nsv, sv);
	  //printf("sv=%f %f %f %f %f\n",sv[1], sv[2], sv[3], sv[4], sv[5]);
	  rcmx = Lx*(sv[1]-0.5);
	  rcmy = Ly*(sv[2]-0.5);
	  rcmz = Lz*(sv[3]-0.5);
	  gamma1 = 2.0*M_PI*sv[4];
	  gamma2 = 2.0*M_PI*sv[5];
#endif
#else
	  /* place second DNAD randomly */
	  rcmx = Lx*(drand48()-0.5);
	  rcmy = Ly*(drand48()-0.5);
	  rcmz = Lz*(drand48()-0.5);
	  gamma1 = 2.0*M_PI*drand48();
	  gamma2 = 2.0*M_PI*drand48();
#endif
	  //if (type==1 && rcmy > 2.0*max2(DNADall[k].sax[0],DNADall[k].sax[1]))
	  //break;
	  place_AMYLOID(0.0, 0.0, 0.0, u1x, u1y, u1z, 0, gamma1);      
	  if (type==0)
	    {
	      if (alpha == 0.0)
		orient(&u2x, &u2y, &u2z);
	      else
		orient_onsager(&u2x, &u2y, &u2z, alpha);
	      //printf("u=%f %f %f\n",u2x, u2y, u2z);
	    }
	  else
	    {
	      if (type==1)
		{
		  if (contrib==0)
		    orient_donsager(&u2x, &u2y, &u2z, alpha,0);
		  else
		    orient_donsager(&u2x, &u2y, &u2z, alpha,1);
		}
	      else
		{
		  if (contrib==0||contrib==2)
		    orient_donsager(&u2x, &u2y, &u2z, alpha,0);
		  else 
		    orient_donsager(&u2x, &u2y, &u2z, alpha,1);
		}
	    }
	  place_AMYLOID(rcmx, rcmy, rcmz, u2x, u2y, u2z, 1, gamma2);
	  //place_AMYLOID(0, 0, 0, u1x, u1y, u1z, 1, gamma1);
#ifdef DEBUG
	  exit(-1);
#endif
	  /* check overlaps */
	  overlap=0;
#ifdef AMYLOID_ELEC
#ifdef PARALLEL
	  if (numtemps > 1 || numconcs > 1)
	    {
	      for (k1=0; k1 < numtemps; k1++)
		for (k2=0; k2 < numconcs; k2++)
		  {
		    uel_arr[k1][k2]=0.0;
		  }		  
	    }
	  else
	    {
	      uel = 0.0;
	    }
	  interact = 0;
#else
	  uel = 0.0;
	  interact = 0;
#endif
#endif

#if 0
	  printf("res=%f\n", calcDistBox());

	  printf("{%f,%f,%f},", DNADall[0].rcm[0], DNADall[0].rcm[1], DNADall[0].rcm[2]);
	  print_matrix(DNADall[0].R,3);
	  printf(",{%f,%f,%f},", DNADall[1].rcm[0], DNADall[1].rcm[1], DNADall[1].rcm[2]);
	  print_matrix(DNADall[1].R,3);
	  printf(",{%f,%f,%f}\n", DNADall[0].sax[0],DNADall[0].sax[1],DNADall[0].sax[2]);
	  exit(-1);
#endif
	  if (calcDistBox(amyloids[0].rcm,amyloids[0].boxsax,amyloids[0].R,
			  amyloids[1].rcm,amyloids[1].boxsax,amyloids[1].R) < 0.0)
	    {
	      if (calcDistAmyloid() < 0.0)
		{
		  //printf("qui\n");
		  overlap = 1;
		}
	    }
#ifdef AMYLOID_ELEC
	  if (!overlap)
	    {
	      for (i=0; i < nat; i++)
		{
		  for (j=i; j < nat; j++)
		    {
		      distsq = Sqr(CHROMs[0][i].x-CHROMs[1][j].x)+Sqr(CHROMs[0][i].y-CHROMs[1][j].y)+
			Sqr(CHROMs[0][i].z-CHROMs[1][j].z);
#if 0
		      if (distsq<0.00001)
			{
			  printf("BOH?!?\n");
			  exit(-1);
			}
#endif
		      if (distsq <= Sqr(CHROMs[0][i].rad+CHROMs[1][j].rad))
			{	
			  overlap = 1;
			  interact = 0;
			  break;
			}
#ifdef PARALLEL
			
		      if (numtemps > 1 || numconcs > 1)
			{
			  yukcutkDsq = maxyukcutkDsq;
			}
#endif
		      if (distsq < yukcutkDsq)
			{
			  interact = 1;
			  //printf("dist=%f yukcutkDsq=%f uel=%f\n", sqrt(distsq), sqrt(yukcutkDsq),calc_yukawa(i, j, distsq));
#if 0
			  if (distsq < Sqr(yukcut/kD))
			    printf("tt=%lld boh... dist=%f sigij=%f yukcut/kD=%f\n", tt, sqrt(distsq), sqrt(sigijsq), yukcut/kD);
#endif
#if 0
			  if (calc_yukawa(i,j,distsq) < 0.0)
			    printf("tt=%lld boh... dist=%f sigij=%f yukcut/kD=%f yuk=%f\n", tt, sqrt(distsq), sqrt(sigijsq), yukcut/kD, calc_yukawa(i,j,distsq));
#endif
#ifdef PARALLEL
			  if (numtemps > 1 || numconcs > 1)
			    {
			      uelcontrib=calc_yukawa_arr(i, j, distsq, &kk);
			      //printf("uelcontrib:%f\n", uelcontrib);
			      if (uelcontrib != 0.0)
				{
				  for (k1=0; k1 < numtemps; k1++)
				    for (k2=0; k2 < numconcs; k2++)
				      {
					if (kk==-1 || 1.0/kD_arr[k1][k2] >= kD_sorted[kk].invkD)  
					  {
					    if (kk==-1)
					      uel_arr[k1][k2] += uelcontrib;
					    else
					      {
#ifdef YUK_CORR
						uel_arr[k1][k2] += yuk_corr_fact_arr[k1][k2]*exp(-kD_arr[k1][k2]*sqrt(distsq))*uelcontrib/epsr(1.0/beta_arr[k1]);
#else
						uel_arr[k1][k2] += exp(-kD_arr[k1][k2]*sqrt(distsq))*uelcontrib/epsr(1.0/beta_arr[k1]);
#endif
					      }
					  }
				      }
				}
			    }
			  else
			    uel += calc_yukawa(i, j, distsq); 
#else
			  uel += calc_yukawa(i, j, distsq); 
#endif
			}

		      /* if no overlap calculate electrostatics contributions */
		    }
		  if (overlap)
		    break;
		}
	    }
#endif
#if 0
	  if (overlap && interact)
	    {
  	      printf("BOH\n");
	      exit(-1);
	    }
#endif
	  if (overlap) 
	    {
	      if (type==1) 
		{
		  if (contrib==0)
		    segno = 1.0;
		  else
		    segno = -1.0;
		}
	      if (type>=2)
		{
#ifdef SYMMETRY
		  if (contrib==0||contrib==3)
		    segno = 2.0; 
		  else
		    segno = -2.0;

#else
		  if (contrib==0||contrib==3)
		    segno = 1.0; 
		  else
		    segno = -1.0;
#endif
		}
	      /* otherwise calculate the integrand */
	      if (type==0)
		vexcl += 1.0;
	      else if (type==1)
		vexcl += segno*u2x*rcmy; /* questo '-' rende negativa la k2 e viene dalla derivata della funzione di Onsager! */
	      else if (type==2) /* K22 */
		{
		  /* il secondo contributo è equivalente al primo poiché posso scambiare gli assi x e y, quindi lo uso 
		   * per mediare */
		  aaa = - segno*u1x*u2x*rcmy*rcmy - segno*u1y*u2y*rcmx*rcmx;// - segno*(u1x*u2y+u1y*u2x)*rcmy*rcmx;
		  //aaa = - segno*u1x*u2x*rcmy*rcmy;
		  vexcl += aaa/2.0;
		}
	      /* NOTA:per ottenere le seguenti espressioni per K11 e K22 basta considerare l'eq. (10)
		 del Phys. Rev. A di Straley del 1976, notando che il versore y nel nostro caso diventa il versore
		 x e che i termini misti (rcmx*rcmy, rcmx*rcmz e rcmy*rcmz) fanno zero per simmetria.
		 Tale equazione di Straley si ottiene considerando che grad(n.n) = 0 dove n è il versore parallelo al direttore
		 nematico e che si puo' sempre scegliere n parallelo all'asse z nel centro di massa della prima particella (R1)
		 e si puo' sempre fare una rotazione intorno all'asse z per annullare i contributi proporzionali a u1y e u2y.
		 Le espressioni che ottenute da Straley vanno bene per colesterici, se si tratta di sistemi nematici 
		 si devono anche considerare i termini nell'espressione di Frank che si ottengono derivando ny. 
		 */
	      else if (type == 3) /* K11 */
		{
		  /* il secondo contributo è equivalente al primo poiché posso scambiare gli assi x e y, quindi lo uso 
		   * per mediare */
#if 1
		  //aaa = - segno*u1x*u2x*rcmx*rcmx;//- segno*(u1x*u2y+u1y*u2x)*rcmy*rcmx;
		  aaa = -segno*u1x*u2x*rcmx*rcmx - segno*u1y*u2y*rcmy*rcmy;//- segno*(u1x*u2y+u1y*u2x)*rcmy*rcmx;
		  vexcl += aaa/2.0;
#else
		  aaa = -segno*(u1x*u2y+u1y*u2x)*rcmy*rcmx;
		  vexcl += aaa;

#endif

		}
	      else /* K33 */
		{
		  vexcl += (-segno*u1x*u2x*rcmz*rcmz - segno*u1y*u2y*rcmz*rcmz)/2.0;
		}
	    }
#ifdef AMYLOID_ELEC
	  else if (interact)
	    {
	      // printf("boh?!? tt=%lld uel=%f\n", tt, uel);	
	      if (type==1) 
		{
		  if (contrib==0)
		    segno = 1.0;
		  else
		    segno = -1.0;
		}
	      if (type>=2)
		{
#ifdef SYMMETRY
		  if (contrib==0||contrib==3)
		    segno = 2.0; 
		  else
		    segno = -2.0;
#else
		  if (contrib==0||contrib==3)
		    segno = 1.0; 
		  else
		    segno = -1.0;
#endif
		}
#ifdef PARALLEL
	      if (numtemps > 1 || numconcs > 1)
		{
		  for (k1=0; k1 < numtemps; k1++)
		    {
		      tempfact = Sqr(epsr(1.0/beta_arr[k1])/beta_arr[k1]); 
		      for (k2=0; k2 < numconcs; k2++)
			{
			  /* otherwise calculate the integrand */
			  if (type==0)
			    vexclel_arr[k1][k2] += (1.0-exp(-beta_arr[k1]*tempfact*uel_arr[k1][k2]));
			  else if (type==1)
			    vexclel_arr[k1][k2] += segno*u2x*rcmy*(1.0-exp(-beta_arr[k1]*tempfact*uel_arr[k1][k2])); /* questo '-' rende negativa la k2 e viene dalla derivata della funzione di Onsager! */
			  else if (type==2)
			    vexclel_arr[k1][k2] += -segno*u1x*u2x*rcmy*rcmy*(1.0-exp(-beta_arr[k1]*tempfact*uel_arr[k1][k2]));
			  else if (type==3)
			    vexclel_arr[k1][k2] += -segno*u1x*u2x*rcmx*rcmx*(1.0-exp(-beta_arr[k1]*tempfact*uel_arr[k1][k2]));
			  else if (type==4)
			    vexclel_arr[k1][k2] += -segno*u1x*u2x*rcmz*rcmz*(1.0-exp(-beta_arr[k1]*tempfact*uel_arr[k1][k2]));
			}
		    }
		}
	      else
		{
		  /* otherwise calculate the integrand */
		  if (type==0)
		    vexclel += (1.0-exp(-beta*uel));
		  else if (type==1)
		    vexclel += segno*u2x*rcmy*(1.0-exp(-beta*uel)); /* questo '-' rende negativa la k2 e viene dalla derivata della funzione di Onsager! */
		  else  if (type==2)
		    vexclel += -segno*u1x*u2x*rcmy*rcmy*(1.0-exp(-beta*uel));
		  else if (type==3)
		    vexclel += -segno*u1x*u2x*rcmx*rcmx*(1.0-exp(-beta*uel));
		  else
		    vexclel += -segno*u1x*u2x*rcmz*rcmz*(1.0-exp(-beta*uel));
		}		
#else
	      /* otherwise calculate the integrand */
	      if (type==0)
		vexclel += (1.0-exp(-beta*uel));
	      else if (type==1)
		vexclel += segno*u2x*rcmy*(1.0-exp(-beta*uel)); /* questo '-' rende negativa la k2 e viene dalla derivata della funzione di Onsager! */
	      else if (type==2)
		vexclel += -segno*u1x*u2x*rcmy*rcmy*(1.0-exp(-beta*uel));
	      else if (type==3)
		vexclel += -segno*u1x*u2x*rcmx*rcmx*(1.0-exp(-beta*uel));
	      else
		vexclel += -segno*u1x*u2x*rcmz*rcmz*(1.0-exp(-beta*uel));
#endif
	      //printf("vexcl:%f vexclel:%f uel=%f\n", vexcl, vexclel, uel);
	    }
#endif
	}
      if (tt > 0 && tt % fileoutits == 0)
	{
#ifdef AMYLOID_ELEC
#ifdef PARALLEL
	  if (numtemps > 1 || numconcs > 1)
	    {
	      for (k1=0; k1 < numtemps; k1++)
		for (k2=0; k2 < numconcs; k2++)
		  {
		    sprintf(fnout, "v%d_c%.0f_T%.0f.dat", type, cchrom_arr[k2], 1.0/beta_arr[k1]);
		    fout = fopen(fnout, "a+");
		    if (type==0)
		      //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
		      fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, Lx*Ly*Lz*(vexcl+vexclel_arr[k1][k2])/((double)tt), Lx*Ly*Lz*vexcl/((double)tt)/1E3,
			      Lx*Ly*Lz*vexclel_arr[k1][k2]/((double)tt));
		    else if (type==1)
		      fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel_arr[k1][k2])/((double)tt))*factor,
			      (Lx*Ly*Lz*vexcl/((double)tt))*factor,
			      (Lx*Ly*Lz*vexclel_arr[k1][k2]/((double)tt))*factor); /* divido per 10^4 per convertire in nm */
		    else
		      fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel_arr[k1][k2])/((double)tt))*Sqr(factor),
			      (Lx*Ly*Lz*vexcl/((double)tt))*Sqr(factor),
			      (Lx*Ly*Lz*vexclel_arr[k1][k2]/((double)tt))*Sqr(factor)); /* divido per 10^5 per convertire in nm */
		    fclose(fout);
		  }
	    }
	  else
	    {
	      fout = fopen(fnout, "a+");
      	      if (type==0)
	    	//fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
		fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, Lx*Ly*Lz*(vexcl+vexclel)/((double)tt), Lx*Ly*Lz*vexcl/((double)tt),
		      	Lx*Ly*Lz*vexclel/((double)tt));
	      else if (type==1)
		fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel)/((double)tt))*factor,
			(Lx*Ly*Lz*vexcl/((double)tt))*factor,
			(Lx*Ly*Lz*vexclel/((double)tt))*factor); /* divido per 10^4 per convertire in nm */
	      else
		fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel)/((double)tt))*Sqr(factor),
			(Lx*Ly*Lz*vexcl/((double)tt))*Sqr(factor),
			(Lx*Ly*Lz*vexclel/((double)tt))*Sqr(factor)); /* divido per 10^5 per convertire in nm */
	      fclose(fout);
      	    }
#else 
	  fout = fopen(fnout, "a+");
	  if (type==0)
	    //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	    fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, Lx*Ly*Lz*(vexcl+vexclel)/((double)tt), Lx*Ly*Lz*vexcl/((double)tt),
		    Lx*Ly*Lz*vexclel/((double)tt));
	  else if (type==1)
	    fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel)/((double)tt))*factor,
		    (Lx*Ly*Lz*vexcl/((double)tt))*factor,
		    (Lx*Ly*Lz*vexclel/((double)tt))*factor); /* divido per 10^4 per convertire in nm */
	  else
	    fprintf(fout,"%lld %.15G %.15G %.15G\n", tt, (Lx*Ly*Lz*(vexcl+vexclel)/((double)tt))*Sqr(factor),
		    (Lx*Ly*Lz*vexcl/((double)tt))*Sqr(factor),
		    (Lx*Ly*Lz*vexclel/((double)tt))*Sqr(factor)); /* divido per 10^5 per convertire in nm */
	  fclose(fout);
#endif
#else
	  fout = fopen(fnout, "a+");
	  if (type==0)
	    //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	    fprintf(fout,"%lld %.15G\n", tt, Lx*Ly*Lz*vexcl/((double)tt));
	  else if (type==1)
	    fprintf(fout,"%lld %.15G\n", tt, (Lx*Ly*Lz*vexcl/((double)tt))*factor); /* divido per 10^4 per convertire in nm */
	  else
	    fprintf(fout,"%lld %.15G\n", tt, (Lx*Ly*Lz*vexcl/((double)tt))*Sqr(factor)); /* divido per 10^5 per convertire in nm */
	  fclose(fout);
#endif
	}
      if (tt % outits==0)
	printf("trials: %lld/%lld\n", tt, tot_trials);
    }
#if defined(QUASIMC) && defined(USEGSL)
  gsl_qrng_free (qsob);   
#endif
}
