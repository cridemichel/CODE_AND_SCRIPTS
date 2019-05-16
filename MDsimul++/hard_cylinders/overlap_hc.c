/* ====================================== >>> PARAMETRI DEFAULT HC <<< ===================================== */
//#define POVRAY
double diamHC=2.0, lengthHC=2.0;
/* ========================================================================================================== */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <float.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
//#define USEGSL
#define GAUSS
#define HC_FAST
#define MC_QUART_HYBRID 
#define EULER_ROT
#define MCGAMMA
#define PRINC_AXES
#define QUASIMC
#define SOBOL_LL /* use long long to have MAXBIT=60! */
#ifndef CMPLX
#define CMPLX(x,y) ((x)+I*(y))
#endif
#ifdef USEGSL
#include <gsl/gsl_qrng.h>
#endif
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#if 0
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#endif
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))
#if 0
static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

#endif
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */

//#define ELEC
//#define NO_INTERP
double **XI1, **XI2, **XI3, **XI4, **XI5, **XI6;
/* ============ MISER ============== */
#define PFAC 0.1
#define MNPT 15
#define MNBS 60
#define TINY 1.0e-30
#define BIG 1.0e30
struct boxS {
  double R[3][3];
  double x[3];
  double sax[3];
};
struct HardCylStr {
  double R[3][3];
  double rcm[3];
  double sax[3];
  double boxsax[3];
  double n[3];
} HardCyl[2];

double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
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

double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
void init_distbox(void)
{
  int k;
  double max_x, max_y, max_z;
  max_x = 0.5*diamHC;
  //printf("max_x=%f len=%f\n", max_x, lengthHC);
  max_y = 0.5*diamHC;
  max_z = 0.5*lengthHC;
 
  //printf("maxax=%f %f %f\n", max_x, max_y, max_z);
  for (k=0; k < 2; k++)
    {
      HardCyl[k].sax[0] = diamHC*0.5;
      HardCyl[k].sax[1] = diamHC*0.5;
      HardCyl[k].sax[2] = lengthHC*0.5;
      HardCyl[k].boxsax[0] = max_x;
      HardCyl[k].boxsax[1] = max_y;
      HardCyl[k].boxsax[2] = max_z;
    }
  //printf("maxx=%f %f\n",CHROMall[0].sax[0],CHROMall[0].sax[1]);
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
void print_pov_header(FILE *fs, double L)
{
  /* povray preamble */
#ifdef POVRAY
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
  fprintf(fs, "location <%f,%f,%f>\n", -L, L, 0.0);
  fprintf(fs, "look_at <0,0,0>\n");
  fprintf(fs, "//focal_point < 1, 1, -6> // pink sphere in focus\n");
  fprintf(fs, "//aperture 0.4\n");
  fprintf(fs, "//blur_samples 20\n");
  fprintf(fs, "}\n");
  fprintf(fs, "#if(1)\n");
  fprintf(fs, "// spotlight light source\n");
  fprintf(fs, "light_source {\n");
  fprintf(fs, "//<0, 10, -3>\n");
  fprintf(fs, "<%f, %f, %f>\n", 0*L, 600.0*L, -600.0*L);
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
#else
  fprintf(fs,".Vol: %G\n", L*L*L);
#endif
}
#if 0
void print_one_amyloid(amyloidS amy, FILE *f, double dr[3])
{
  int i, jj, kk, k1, k2, k3;
  double R[3][3], ux[3], uy[3], uz[3], x1[3], x2[3], x1p[3], x2p[3], x1pp[3], x2pp[3], xbox[3], rb[3];

  
  for (jj=0; jj < amy.nL; jj++)
    {
      for (kk=0; kk < 3; kk++)
	rb[kk] = amy.rcm[kk]-dr[kk];

      body2lab(amy.boxes[jj].x,xbox,rb,amy.R);

      //printf(">>>rcm %f %f %f\n", amy.boxes[jj].x[0],amy.boxes[jj].x[1],amy.boxes[jj].x[2]);
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
#ifdef POVRAY
      fprintf(f,"box\n");	
      fprintf(f,"{\n");
      for (kk=0; kk < 3; kk++)
	{
	  x1p[kk] = amy.boxes[jj].sax[kk];
	  x2p[kk] = -amy.boxes[jj].sax[kk];
	}	  
      fprintf(f,"<%.15G,%.15G,%.15G>, <%.15G,%.15G,%.15G>\n",x1p[0],x1p[1],x1p[2], x2p[0],x2p[1],x2p[2]);
      //fprintf(f," pigment { color rgb<0.0,0.9,0.2> transmit 0.0 }\n");
      fprintf(f,"pigment { Red transmit 0.0 }\n");
      //fprintf(f,"rotate <0,0,90>\n");
      fprintf(f, "matrix <%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,0,0,0>\n",
	      R[0][0],R[0][1],R[0][2],R[1][0],R[1][1],R[1][2],R[2][0],R[2][1],R[2][2]);
      fprintf(f,"translate <%.15G, %.15G, %.15G>\n", xbox[0], xbox[1], xbox[2]);
      /* left-handed to right-handed, since POVRAY coordinate system is left-handed */
      fprintf(f, "matrix <1,0,0,0,1,0,0,0,-1,0,0,0>\n");
      fprintf(f,"finish { phong PHONG phong_size PHONG_SIZE reflection REFL\n"); 
      fprintf(f,"ambient AMB diffuse DIFF\n");
      fprintf(f,"specular SPEC roughness ROUGH}\n");
      fprintf(f, "}\n");
#else
      fprintf(f, "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G B %f %f %f C[%s]\n", 
	   xbox[0], xbox[1], xbox[2], 
	   R[0][0],R[0][1],R[0][2],R[1][0],R[1][1],R[1][2],R[2][0],R[2][1],R[2][2],
	   2.0*amy.boxes[jj].sax[0], 2.0*amy.boxes[jj].sax[1], 2.0*amy.boxes[jj].sax[2], "red");
#endif
    }
}


void saveAmyloidPovray(char *fname, double L)
{
  FILE *f;
  int kk;
  double dr[3]={0.0,0.0,0.0};
  f = fopen(fname, "w+");
  print_pov_header(f, L);

  for (kk=0; kk < 3; kk++)
    dr[kk] = (amyloids[0].rcm[kk]+amyloids[1].rcm[kk])*0.5;

  print_one_amyloid(amyloids[0], f, dr);
  //print_one_amyloid(amyloids[1], f, dr);
  fclose(f);
}
#endif
//void versor_to_R(double ox, double oy, double oz, double R[3][3], double gamma);
void versor_to_R_sym(double ox, double oy, double oz, double R[3][3]);


void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
void print_vec(char *lbl, double *V)
{
  printf("%s(%.15G,%.15G,%.15G)\n", lbl, V[0], V[1], V[2]);
}
double perpcomp(double *V, double *C, double *n)
{
  int kk2;
  double dsc[3], sp, dscperp[3];
  for (kk2=0; kk2 < 3; kk2++)
    dsc[kk2] = V[kk2] - C[kk2]; 
  sp = scalProd(dsc, n);
  for (kk2=0; kk2 < 3; kk2++)
    dscperp[kk2] = dsc[kk2]-sp*n[kk2];
  return calc_norm(dscperp);
}
void body2labHC(double xp[3], double x[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      /* NOTE: k2 starts from 1 because xp[0] = 0.0 see function find_initial_guess() below */
      for (k2=1; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
      x[k1] += rO[k1];
    }
}
void versor_to_R_sym(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double u[3], sp, norm;
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
#ifdef DEBUG
  printf("==============\n");
  print_matrix(R, 3);
  printf("==============\n");
#endif
}
void find_initial_guess(double *Ai, double Ci[3], double ni[3], double Dj[3], double nj[3], double D)
{
  const int meshpts = 8;
  double Pj[3], Rj[3][3];
  int kk, nn;
  static int firstcall=1;
  double th, dth, xp[3], Ui[3], UiPj[3];
  static double **mesh; /* {{1,0},{0.707106781186547, 0.707106781186547},{0,1},
      {-0.707106781186547,0.707106781186547},{-1,0},{-0.707106781186547,-0.707106781186547},
      {0,-1},{0.707106781186547,-0.707106781186547}};*/
  double PjCini, PjCi[3], d, mindist=-1.0; 
  versor_to_R_sym(nj[0],nj[1],nj[2], Rj); 
#if 1
  if (firstcall)
    {
      mesh = (double **)malloc(sizeof(double*)*meshpts);
      for (nn=0; nn < meshpts; nn++)
	mesh[nn] = (double*)malloc(sizeof(double)*3);
      firstcall=0;
      dth = acos(0)*4.0/((double)meshpts);

      th=0.0;
      for (nn=0; nn < meshpts; nn++)
	{
	  mesh[nn][0] = cos(th);
	  mesh[nn][1] = sin(th);
	  th += dth;
	}
    }
#endif
  for (nn=0; nn < meshpts; nn++)
    {
      //xp[0] = 0.0;
      xp[1] = D*0.5*mesh[nn][0];
      xp[2] = D*0.5*mesh[nn][1];
      body2labHC(xp, Pj, Dj, Rj);    
      for (kk=0; kk < 3; kk++)
	PjCi[kk] = Pj[kk] - Ci[kk];
      //normPjCi = calc_norm(PjCi);
      PjCini = scalProd(PjCi,ni);
      for (kk=0; kk < 3; kk++)
	{
	  Ui[kk] = Ci[kk] + PjCini*ni[kk];
	  UiPj[kk] = Ui[kk]-Pj[kk];
	}
      if ((d=calc_norm(UiPj)) < mindist || nn==0)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      Ai[kk] = Ui[kk];
    	    }
	  mindist=d;
	  //printf("nn=%d mindist=%.15G d=%.15G\n", nn, mindist, d);
	  //printf("Ui=%f %f %f Pi=%f %f %f\n", Ui[0],Ui[1], Ui[2], Pj[0], Pj[1], Pj[2]);
	}
    }
  //printf("done\n");
#if 0
   for (kk=0; kk < 3; kk++)
     AiCi[kk]  = Ai[kk] - Ci[kk]; 
  printf("norm AiCi=%.15G sp=%.15G\n", calc_norm(AiCi), scalProd(AiCi,ni)/calc_norm(AiCi));
  for (kk=0; kk < 3; kk++)
    AiCi[kk]  = Pj[kk] - Dj[kk]; 

  printf("norm AiCi=%.15G sp=%.15G\n", calc_norm(AiCi), scalProd(AiCi,nj));
#endif 
}
double diskdisk(double D, double L, double Di[2][3], double Ci[3], double ni[3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, j2, kk; 
  double sp, normCiCj, CiCj[3], VV[3], Q1;
  double DiN, DjN, niN[3], njN[3], Djnj, Dini, Pi[3], N[3], Pj[3], normN;
  double normNSq, PiPj[3], normPiPj, Q2, PjDj[3], normPiDi, normPjDj, PiDi[3];
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= L)
	    return -1;
	  else
	    return 1;
	}

      /* parallel disks */
      for (j1=0; j1 < 2; j1++)
	for (j2=j1; j2 < 2; j2++)
	  {
	    sp=0.0;
	    for (kk=0; kk < 3; kk++)
	      {
		VV[kk] = Di[j1][kk]-Dj[j2][kk];
		sp += ni[kk]*VV[kk];
	      }
	    if (sp == 0 && calc_norm(VV) < D)
	      {
		return -1;
	      }
	  }
    }
  else 
    {
   
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
#if 0
	    sp=0;
	    for (kk=0; kk < 3; kk++)
	      {
		sp  += Sqr(Di[j1][kk]-Dj[j2][kk]);
	      }	  
	    if (sp > D*D)
	      continue;
#endif
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
	    if (normPiDi <= 0.5*D && normPjDj <= 0.5*D)
	      {
		Q1 = sqrt(Sqr(D)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  return 0;
}
double rimrim(double D, double L, double Ci[3],double ni[3], double Cj[3], double nj[3])
{
  int kk;
  double ViVj[3], lambdai, lambdaj, ninj;
  double CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3]; 
  /* case A.3 rim-rim overlap */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if ( //Sqr(ViVj[0])+Sqr(ViVj[1])+Sqr(ViVj[2])<Sqr(D) 
       calc_norm(ViVj) < D 
       && fabs(lambdai) < 0.5*L && fabs(lambdaj) < 0.5*L)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 0;
}
void solve_cubic_analytic(double *coeff, complex double* sol)
{
  /* solve the cubic coeff[3]*x^3 + coeff[2]*x^2 +  coeff[1]*x + coeff[0] = 0
   * according to the method described in Numerical Recipe book */  
  double a, b, c, Q, R, theta, Q3, R2, A, B;
  const double sqrt32=sqrt((double)3.0)/2.0;
  a = coeff[2]/coeff[3];
  b = coeff[1]/coeff[3];
  c = coeff[0]/coeff[3];
  Q = (Sqr(a) - 3.0*b)/9.0;
  R = (2.0*Sqr(a)*a - 9.0*a*b + 27.0*c)/54.0;
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sol[0] = -2.0*sqrt(Q)*cos(theta/3.0)- a/3.0;
      sol[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a/3.0;
      sol[2] = -2.0*sqrt(Q)*cos((theta-2.0*M_PI)/3.0) - a/3.0;
    }
  else
    {
      A = -copysign((double)1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = (A+B) - a/3.0;
      sol[1] = CMPLX(-0.5*(A+B)-a/3.0,sqrt32*(A-B));
      sol[2] = conj(sol[1]);
      //sol[1] = -0.5*(A+B)-a/3.0+cmplx(0,1)*sqrt32*(A-B);
      //sol[2] = -0.5*(A+B)-a/3.0-cmplx(0,1)*sqrt32*(A-B);
    }
}
void oqs_solve_quadratic(double a, double b, complex double roots[2]);

void solve_cubic(double *coeff, int *numsol, double sol[3])
{
  int k;
  /* solution from Abramovitz */
  double a, b;
  double complex csol[3];
  if (coeff[3] == 0.0)
    {
      printf("[WARNING] fallback to quadratic from cubic\n");
      //solve_quadratic(coeff, numsol, sol);
      a = coeff[1]/coeff[2];
      b = coeff[0]/coeff[2];
      oqs_solve_quadratic(a, b, csol);
      for (k=0; k < 2; k++)
        {
          if (cimag(csol[k])==0)
            {
              sol[*numsol] = creal(csol[k]);
              (*numsol)++;
            }
        }     
      return;
    }
  /* la soluzione numnerica è comunque più accurata di quella analitica */
  solve_cubic_analytic(coeff, csol); 
  for (k=0; k < 3; k++)
    {
      if (cimag(csol[k])==0)
	{
	  sol[*numsol] = creal(csol[k]);
	  (*numsol)++;
	}
    }
}
const double cubic_rescal_fact = 3.488062113727083E+102; //= pow(DBL_MAX,1.0/3.0)/1.618034;
const double quart_rescal_fact = 7.156344627944542E+76; // = pow(DBL_MAX,1.0/4.0)/1.618034;
const double macheps = 2.2204460492503131E-16; // DBL_EPSILON
double oqs_max2(double a, double b)
{
  if (a >= b)
    return a;
  else
    return b;
}
double oqs_max3(double a, double b, double c)
{
  double t;
  t = oqs_max2(a,b);
  return oqs_max2(t,c);
}

void oqs_solve_cubic_analytic_depressed_handle_inf(double b, double c, double *sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  double Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const double PI2=M_PI/2.0, TWOPI=2.0*M_PI;
  Q = -b/3.0;
  R = 0.5*c;
  if (R==0)
    {
      if (b <= 0)
        {
          *sol=sqrt(-b);
        }
      else
        {
          *sol=0;
        }
      return;
    }

  if (fabs(Q) < fabs(R))
    {
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      RQ = R/Q;
      KK = copysign(1.0,Q)*(RQ*RQ/Q-1.0);
    }

  if (KK < 0.0)
    {
      sqrtQ=sqrt(Q);
      theta = acos((R/fabs(Q))/sqrtQ);
      if (theta < PI2) 
        *sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
        *sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (fabs(Q) < fabs(R))
        A = -copysign(1.0,R)*cbrt(fabs(R)*(1.0+sqrt(KK)));
      else
        {
          A = -copysign(1.0,R)*cbrt(fabs(R)+sqrt(fabs(Q))*fabs(Q)*sqrt(KK));
        }
      if (A==0.0)
        B=0.0;
      else
        B = Q/A;
      *sol = A+B;
    }
}
void oqs_solve_cubic_analytic_depressed(double b, double c, double *sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  double Q, R, theta, Q3, R2, A, B, sqrtQ;
  Q = -b/3.0;
  R = 0.5*c;
  if (fabs(Q) > 1E102 || fabs(R) > 1E154)
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sqrtQ=-2.0*sqrt(Q);
      if (theta < M_PI/2) 
        *sol = sqrtQ*cos(theta/3.0);
      else 
        *sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
        B=0.0;
      else
        B = Q/A;
      *sol = A+B; /* this is always largest root even if A=B */
    }
}
void oqs_calc_phi0(double a, double b, double c, double d, double *phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (79) (see also the discussion in sec. 2.2 of the manuscript) */
  double rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  double maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  double ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
  int iter;
  diskr=9*a*a-24*b;                    
  /* eq. (87) */
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
        s=-2*b/(3*a+diskr);                     
      else
        s=-2*b/(3*a-diskr);                      
    }
  else
    {      
      s=-a/4;                                    
    }
  /* eqs. (83) */
  aq=a+4*s;                                      
  bq=b+3*s*(a+2*s);                              
  cq=c+s*(2*b+s*(3*a+4*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9;
  hh=aq*cq;     

  g=hh-4*dq-3*gg;                       /* eq. (85) */  
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq; /* eq. (86) */          
  oqs_solve_cubic_analytic_depressed(g, h, &rmax);
  if (isnan(rmax) || isinf(rmax))
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
      if ((isnan(rmax) || isinf(rmax)) && scaled)
        {
          // try harder: rescale also the depressed cubic if quartic has been already rescaled
          rfact = cubic_rescal_fact; 
          rfactsq = rfact*rfact;
          ggss = gg/rfactsq;
          hhss = hh/rfactsq;
          dqss = dq/rfactsq;
          aqs = aq/rfact;
          bqs = bq/rfact;
          cqs = cq/rfact;
          ggss=bqs*bqs/9.0;
          hhss=aqs*cqs;   
          g=hhss-4.0*dqss-3.0*ggss;                       
          h=(8.0*dqss+hhss-2.0*ggss)*bqs/3-cqs*(cqs/rfact)-(dq/rfact)*aqs*aqs; 
          oqs_solve_cubic_analytic_depressed(g, h, &rmax);
          if (isnan(rmax) || isinf(rmax))
            {
              oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
            }
          rmax *= rfact;
        }
    }
  /* Newton-Raphson used to refine phi0 (see end of sec. 2.2 in the manuscript) */
  x = rmax;
  xsq=x*x;
  xxx=x*xsq;
  gx=g*x;
  f = x*(xsq + g) + h;
  if (fabs(xxx) > fabs(gx))
    maxtt = fabs(xxx);
  else
    maxtt = fabs(gx);
  if (fabs(h) > maxtt)
    maxtt = fabs(h);

  if (fabs(f) > macheps*maxtt)
    {
      for (iter=0; iter < 8; iter++)
        {   
          df =  3.0*xsq + g;
          if (df==0)
            {
              break;
            }
          xold = x;
          x += -f/df;
          fold = f;
          xsq = x*x;
          f = x*(xsq + g) + h;
          if (f==0)
            {
              break;
            } 

          if (fabs(f) >= fabs(fold))
            {
              x = xold;
              break;
            }
        }
    }
  *phi0 = x;
}
double oqs_calc_err_ldlt(double b, double c, double d, double d2, double l1, double l2, double l3)
{
  /* Eqs. (29) and (30) in the manuscript */
  double sum;
  sum =  (b==0)?fabs(d2 + l1*l1 + 2.0*l3):fabs(((d2 + l1*l1 + 2.0*l3)-b)/b);
  sum += (c==0)?fabs(2.0*d2*l2 + 2.0*l1*l3):fabs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  sum += (d==0)?fabs(d2*l2*l2 + l3*l3):fabs(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
}
double oqs_calc_err_abcd_cmplx(double a, double b, double c, double d, 
                               complex double aq, complex double bq, complex double cq, complex double dq)
{
  /* Eqs. (68) and (69) in the manuscript for complex alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq) */
  double sum;
  sum = (d==0)?cabs(bq*dq):cabs((bq*dq-d)/d);
  sum += (c==0)?cabs(bq*cq + aq*dq):cabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?cabs(bq + aq*cq + dq):cabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?cabs(aq + cq):cabs(((aq + cq) - a)/a);
  return sum;
}
double oqs_calc_err_abcd(double a, double b, double c, double d, double aq, double bq, double cq, double dq)
{
  /* Eqs. (68) and (69) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  double sum;
  sum = (d==0)?fabs(bq*dq):fabs((bq*dq-d)/d);
  sum += (c==0)?fabs(bq*cq + aq*dq):fabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabs(bq + aq*cq + dq):fabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabs(aq + cq):fabs(((aq + cq) - a)/a);
  return sum;
}
double oqs_calc_err_abc(double a, double b, double c, double aq, double bq, double cq, double dq)
{
  /* Eqs. (48)-(51) in the manuscript */
  double sum;
  sum = (c==0)?fabs(bq*cq + aq*dq):fabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabs(bq + aq*cq + dq):fabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabs(aq + cq):fabs(((aq + cq) - a)/a);
  return sum;
}
void oqs_NRabcd(double a, double b, double c, double d, double *AQ, double *BQ, double *CQ, double *DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  double x02, errf, errfold, xold[4], x[4], dx[4], det, Jinv[4][4], fvec[4], vr[4];
  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
  vr[0] = d;
  vr[1] = c;
  vr[2] = b;
  vr[3] = a;
  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  for (k1=0; k1 < 4; k1++)
    {
      errf += (vr[k1]==0)?fabs(fvec[k1]):fabs(fvec[k1]/vr[k1]);
    }
  for (iter = 0; iter < 8; iter++)
    {
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - 2.0*x[3]) + x[3]*(x[0]*x02 + x[3]);
      if (det==0.0)
        break;
      Jinv[0][0] = x02;
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = -x[1]*Jinv[0][1] - x[0]*Jinv[0][2]; 
      Jinv[1][0] = x[0]*Jinv[0][0] + Jinv[0][1];
      Jinv[1][1] = -x[1]*Jinv[0][0];
      Jinv[1][2] = -x[1]*Jinv[0][1];   
      Jinv[1][3] = -x[1]*Jinv[0][2];
      Jinv[2][0] = -Jinv[0][0];
      Jinv[2][1] = -Jinv[0][1];
      Jinv[2][2] = -Jinv[0][2];
      Jinv[2][3] = Jinv[0][2]*x[2] + Jinv[0][1]*x[3];
      Jinv[3][0] = -x[2]*Jinv[0][0] - Jinv[0][1];
      Jinv[3][1] = Jinv[0][0]*x[3];
      Jinv[3][2] = x[3]*Jinv[0][1];
      Jinv[3][3] = x[3]*Jinv[0][2];
      for (k1=0; k1 < 4; k1++)
        {
          dx[k1] = 0;
          for (k2=0; k2 < 4; k2++)
            dx[k1] += Jinv[k1][k2]*fvec[k2];
        }
      for (k1=0; k1 < 4; k1++)
        xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
        {
          x[k1] += -dx[k1]/det;
        }
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 
      errfold = errf;
      errf=0;
      for (k1=0; k1 < 4; k1++)
        {
          errf += (vr[k1]==0)?fabs(fvec[k1]):fabs(fvec[k1]/vr[k1]);
        }
      if (errf==0)
        break;
      if (errf >= errfold)
        {
          for (k1=0; k1 < 4; k1++)
            x[k1] = xold[k1];
          break;
        }
    }
  *AQ=x[0];
  *BQ=x[1];
  *CQ=x[2];
  *DQ=x[3];
}
void oqs_solve_quadratic(double a, double b, complex double roots[2])
{ 
  double div,sqrtd,diskr,zmax,zmin;
  diskr=a*a-4*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
        div=-a-sqrt(diskr);
      else
        div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
        zmin=0.0;
      else
        zmin=b/zmax;

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=CMPLX(-a/2,sqrtd/2);
      roots[1]=CMPLX(-a/2,-sqrtd/2);      
    }   
}
void oqs_quartic_solver(double coeff[5], complex double roots[4])      
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * 
   * if coeff[4] != 0 
   *
   * the four roots will be stored in the complex array roots[] 
   *
   * */
  complex double acx1, bcx1, ccx1, dcx1,acx,bcx,ccx,dcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
  double l2m[12], d2m[12], res[12], resmin, bl311, dml3l3, err0=0, err1=0, aq1, bq1, cq1, dq1; 
  double a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3],gamma,del2;
  int realcase[2], whichcase, k1, k, kmin, nsol;
  double rfactsq, rfact=1.0;

  if (coeff[4]==0.0)
    {
      printf("That's not a quartic!\n");
      return;
    }
  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  oqs_calc_phi0(a,b,c,d,&phi0, 0);

  // simple polynomial rescaling
  if (isnan(phi0)||isinf(phi0))
    {
      rfact = quart_rescal_fact; 
      a /= rfact;
      rfactsq = rfact*rfact;
      b /= rfactsq;
      c /= rfactsq*rfact;
      d /= rfactsq*rfactsq;
      oqs_calc_phi0(a,b,c,d,&phi0, 1);
    }
  l1=a/2;          /* eq. (16) */                                        
  l3=b/6+phi0/2;   /* eq. (18) */                                
  del2=c-a*l3;     /* defined just after eq. (27) */                             
  nsol=0;
  bl311 =2.*b/3.-phi0-l1*l1;   /* This is d2 as defined in eq. (20)*/ 
  dml3l3 = d-l3*l3;            /* dml3l3 is d3 as defined in eq. (15) with d2=0 */ 

  /* Three possible solutions for d2 and l2 (see eqs. (28) and discussion which follows) */
  if (bl311!=0.0)
    {
      d2m[nsol] = bl311;  
      l2m[nsol] = del2/(2.0*d2m[nsol]);   
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }
  if (del2!=0)
    {
      l2m[nsol]=2*dml3l3/del2;
      if (l2m[nsol]!=0)
        {
          d2m[nsol]=del2/(2*l2m[nsol]);
          res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
          nsol++;
        }

      d2m[nsol] = bl311;
      l2m[nsol] = 2.0*dml3l3/del2;
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }

  if (nsol==0)
    {
      l2=d2=0.0;
    }
  else
    {
      /* we select the (d2,l2) pair which minimizes errors */
      for (k1=0; k1 < nsol; k1++)
        {
          if (k1==0 || res[k1] < resmin)
            {
              resmin = res[k1];
              kmin = k1;      
            }
        }
      d2 = d2m[kmin];
      l2 = l2m[kmin];
    }
  whichcase = 0; 
  if (d2 < 0.0) 
    {
      /* Case I eqs. (37)-(40) */
      gamma=sqrt(-d2);                               
      aq=l1+gamma;                                  
      bq=l3+gamma*l2;                              

      cq=l1-gamma;                                
      dq=l3-gamma*l2;                            
      if(fabs(dq) < fabs(bq))
        dq=d/bq;                                
      else if(fabs(dq) > fabs(bq))
        bq=d/dq;                               
      if (fabs(aq) < fabs(cq))
        {
          nsol=0;
          if (dq !=0)
            {
              aqv[nsol] = (c - bq*cq)/dq;    /* see eqs. (47) */
              errv[nsol]=oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
              nsol++;
            }
          if (cq != 0) 
            {
              aqv[nsol] = (b - dq - bq)/cq;  /* see eqs. (47) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
              nsol++;
            }
          aqv[nsol] = a - cq;                /* see eqs. (47) */
          errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
          nsol++;
          /* we select the value of aq (i.e. alpha1 in the manuscript) which minimizes errors */
          for (k=0; k < nsol; k++)
            {
              if (k==0 || errv[k] < errmin)
                {
                  kmin = k;
                  errmin = errv[k];
                }
            }
          aq = aqv[kmin];
        }
      else 
        {
          nsol = 0;
          if (bq != 0)
            { 
              cqv[nsol] = (c - aq*dq)/bq;              /* see eqs. (53) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
              nsol++;
            }
          if (aq != 0)
            {
              cqv[nsol] = (b - bq - dq)/aq;            /* see eqs. (53) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
              nsol++;
            }
          cqv[nsol] = a - aq;                          /* see eqs. (53) */
          errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
          nsol++;   
          /* we select the value of cq (i.e. alpha2 in the manuscript) which minimizes errors */
          for (k=0; k < nsol; k++)
            {
              if (k==0 || errv[k] < errmin)
                {
                  kmin = k;
                  errmin = errv[k];
                }
            }
          cq = cqv[kmin];
        }
      realcase[0]=1;
    }
  else if (d2 > 0)   
    {
      /* Case II eqs. (53)-(56) */
      gamma=sqrt(d2); 
      acx=CMPLX(l1,gamma);  
      bcx=CMPLX(l3,gamma*l2);
      ccx = conj(acx);
      dcx = conj(bcx);
      realcase[0] = 0; 
    }
  else 
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is better) */
  if (realcase[0]==-1 || (fabs(d2) <= macheps*oqs_max3(fabs(2.*b/3.), fabs(phi0), l1*l1))) 
    {
      d3 = d - l3*l3;
      if (realcase[0]==1)
        err0 = oqs_calc_err_abcd(a, b, c, d, aq, bq, cq, dq);
      else if (realcase[0]==0)
        err0 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);
      if (d3 <= 0)
        {
          realcase[1] = 1;
          aq1 = l1;   
          bq1 = l3 + sqrt(-d3);
          cq1 = l1;
          dq1 = l3 - sqrt(-d3);
          if(fabs(dq1) < fabs(bq1))  
            dq1=d/bq1;                                        
          else if(fabs(dq1) > fabs(bq1))
            bq1=d/dq1;                                       
          err1 = oqs_calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (68) */
        }
      else /* complex */
        {
          realcase[1] = 0;
          acx1 = l1;
          bcx1 = l3 + I*sqrt(d3);
          ccx1 = l1;
          dcx1 = conj(bcx1);
          err1 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1); 
        }
      if (realcase[0]==-1 || err1 < err0)
        {
          whichcase=1; // d2 = 0
          if (realcase[1]==1)
            {
              aq = aq1;
              bq = bq1;
              cq = cq1;
              dq = dq1;
            }
          else
            {
              acx = acx1;
              bcx = bcx1;
              ccx = ccx1;
              dcx = dcx1;
            }
        }
    }
  if (realcase[whichcase]==1)
    {
      /* if alpha1, beta1, alpha2 and beta2 are real first refine 
       * the coefficient through a Newton-Raphson */
      oqs_NRabcd(a,b,c,d,&aq,&bq,&cq,&dq);      
      /* finally calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
      oqs_solve_quadratic(aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      oqs_solve_quadratic(cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];
    }
  else
    {
      /* complex coefficients of p1 and p2 */
      if (whichcase==0) // d2!=0
        {
          cdiskr=acx*acx/4-bcx;               
          /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
          zx1=-acx/2+csqrt(cdiskr);
          zx2=-acx/2-csqrt(cdiskr);
          if(cabs(zx1) > cabs(zx2))
            zxmax=zx1;
          else
            zxmax=zx2;
          zxmin=bcx/zxmax;        
          roots[0]=zxmin;
          roots[1]=conj(zxmin);
          roots[2]=zxmax;
          roots[3]=conj(zxmax);
        }
      else // d2 ~ 0
        {
          /* never gets here! */
          cdiskr=csqrt(acx*acx-4.0*bcx);
          zx1 = -0.5*(acx+cdiskr);
          zx2 = -0.5*(acx-cdiskr);
          if (cabs(zx1) > cabs(zx2))
            zxmax = zx1;
          else
            zxmax = zx2;
          zxmin = bcx/zxmax;
          roots[0] = zxmax;
          roots[1] = zxmin;
          cdiskr=csqrt(ccx*ccx-4.0*dcx);
          zx1 = -0.5*(ccx+cdiskr);
          zx2 = -0.5*(ccx-cdiskr);
          if (cabs(zx1) > cabs(zx2))
            zxmax = zx1;
          else
            zxmax = zx2;
          zxmin = dcx/zxmax;
          roots[2]= zxmax;
          roots[3]= zxmin;
        }
    }
  if (rfact!=1.0)
    {
      for (k=0; k < 4; k++)
        roots[k] *= rfact;
    }
}
void wrap_oqs_quartic_solver(double coeff[5], int *numsol, double solqua[4])
{
  //double ERRTHR=1E-10;
  double complex csol[4];
  int k; //, printsol;
  oqs_quartic_solver(coeff, csol);
  //quarticRoots(coeff, numsol, csol);
  *numsol=0;
#if 0
  for (k=0; k < 4; k++)
    {
      if (isnan(creal(csol[k]))||isnan(cimag(csol[k]))|| 
	  isinf(creal(csol[k]))||isinf(cimag(csol[k])))
	planB=1;
    }
 // if error is too large use quadruple precision HQR routine as a safe (although much slower) backup */ 
  if (planB)
    {
      //printf("[USING HQRL] BEFORE sum=%.15G\n", sum);
      for (k=0; k < 4; k++)
	{
	  printf("BEFORE csol=%.16G + I*(%.16G)\n", creal(csol[k]), cimag(csol[k]));
	}
      for (k=0; k < 5; k++)
	coeffl[k]=coeff[k];
      solve_numrecl_cmplx(coeffl, csoll, 4);
      for (k=0; k < 4; k++)
	{
	  csol[k] = csoll[k];
	  printf("AFTER csol=%.16G + I*(%.16G)\n", creal(csol[k]), cimag(csol[k]));
	}
      sum2=0;
      for (k=0; k < 4; k++)
	{
	  sum2 += cabs(quarticfunc(coeff, csol[k]));
	}
      printf("[USING HQRL] AFTER sum=%.20G\n", sum2);
    }
#endif
  for (k=0; k < 4; k++)
    {
      if (cimag(csol[k])==0)
	{
	  solqua[*numsol] = creal(csol[k]);
	  (*numsol)++;
	}
    }
}
void solve_quartic(double coeff[5], int *numsol, double solqua[4])
{
  //double EPS=2.2204460492503131E-16;
  if (coeff[4]==0)
    {
      solve_cubic(coeff, numsol, solqua);
      return;
    }
  wrap_oqs_quartic_solver(coeff, numsol, solqua);   
}
void discard_spurious(double *solqua, int *numsol)
{
  /* each solution x must be such that |x| <= 1 */
  int k, nsol;
  const double EPS=5E-8; // it should be at least not greater than DIST_THR below
  double solL[4];
  nsol=0;
  for (k=0; k < *numsol; k++)
    solL[k] = solqua[k];
  for (k=0; k < *numsol; k++)
    {
      if (fabs(solL[k]) < 1.0+EPS)
	{
	  solqua[nsol] = solL[k];
	  nsol++;
	}
    }
#if 0
  if (numsol==4)
    printf("numsol=%d nsol=%d\n", *numsol, nsol);
#endif
  *numsol = nsol;
}
int test_for_fallback(double *P, double *Cip, double *nip, double D2, double *diff)
{
#ifdef MC_QUART_USE_ANALYTIC
  const double DIST_THR=1E-4;
#else
#ifdef MC_QUART_HYBRID
#ifdef FAST_QUARTIC_SOLVER
  const double DIST_THR=5E-13;// old value 5.0E-13 -- NOTA 22/02/18 con il nuovo quartic solver ci si può spingere a 5E-14 volendo!
#else
  const double DIST_THR=5E-12;
#endif
#else
  const double DIST_THR=5E-12;
#endif
#endif
  double diff1, diff2;
  diff1=fabs(perpcomp(P, Cip, nip)-D2); // qui D2 è il diametro del rim
  diff2=fabs(sqrt(Sqr(P[1])+Sqr(P[2]))-D2);// qui D2 è il diametro del disco

  *diff=diff1+diff2;
  if (diff1 > DIST_THR*D2 || diff2 > DIST_THR*D2)
    return 1;
  else
    return 0;
#if 0
  const double DIST_THR=5E-8;
  if ((*diff=fabs(perpcomp(P, Cip, nip)-D2)) > DIST_THR*D2)
    return 1;
  else 
    return 0;
#endif
}
double ranf(void)
{
  return drand48();
}
void rotate_axes_on_plane(double RR[3][3])
{
  double Rin[3][3];
  int k1, k2, k3;
  double ox, oy, oz, theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];

  ox = RR[0][0];
  oy = RR[0][1];
  oz = RR[0][2];
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Rin[k1][k2]=RR[k1][k2];
      }
  //NOTA 15/01/2018
  //se uso ranf altero la sequenza casuale e perdo il confronto la simulazione da 50x10^6 già fatta con il metodo di Alberto 
  //theta = (ranf()>0.5?1.:-1.)*M_PI/4.0;
  theta = M_PI/4.0;
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
     RR[k1][k2] = Ro[k1][k2];
}
void versor_to_R_alt(double *Ci, double *ni, double *Dj, double *nj, double R[3][3], double D)
{
  int k, kk1, kk;
  double u[3]; 
  double normDjCi, DjCi[3], DjCini, Ai[3], AiDjnj, AiDj[3], VV[3], ragg;
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

  /* N.B. qui faccio uno step dell'algorito di Ibarra semplificato per
   * determinare l'asse y del riferimenti del disco (l'asse x è l'asse perpendicolare
   * al disco e l'asse z si ottiene con il prodotto vettore dell'asse x e y) */
#if 1
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  normDjCi = calc_norm(DjCi);
  DjCini = scalProd(DjCi,ni);

  for (kk1 = 0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
  for (kk1=0; kk1 < 3; kk1++)
    AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
  AiDjnj = scalProd(AiDj, nj);
  for (kk1=0; kk1 < 3; kk1++)
    {
      VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
    }
  //for (kk1=0; kk1 < 3; kk1++)
    //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
  ragg = calc_norm(VV);

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k]/ragg;
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
#if 0
  ragg = scalProd(TnCi,ni);
  for (k=0;k<3;k++)
    Ai[k] = Ci[k] + ragg*ni[k];
#endif
#else

 for (k=0; k < 3; k++)
    dsc[k] = Ci[k] - Dj[k]; 
  sp = scalProd(dsc, nj);
  for (k=0; k < 3; k++)
    R[1][k] = dsc[k] - sp*nj[k];
#endif  
  //printf("scalProd=%.15G\n", scalProd(R[1],R[0]));
  vectProdVec(R[0], R[1], u);
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
void versor_to_R_alt_fb(double *Ci, double *ni, double *Dj, double *nj, double R[3][3], double D, double *Tj, int MAXITS)
{
  int k, kk1, kk, its;
  double u[3]; 
  double normDjCi, DjCi[3], DjCini, Ai[3], AiDjnj, AiDj[3], VV[3], ragg;
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

  /* N.B. Questa viene utilizzata in fallback e faccio un ulteriore step dell'algorito di Ibarra semplificato per
   * determinare l'asse y del riferimenti del disco (l'asse x è l'asse perpendicolare
   * al disco e l'asse z si ottiene con il prodotto vettore dell'asse x e y) */
#if 1
  for (its = 0; its < MAXITS; its++)
    {
      for (kk=0; kk < 3; kk++)
	DjCi[kk] = Dj[kk] - Ci[kk];
      normDjCi = calc_norm(DjCi);
      DjCini = scalProd(DjCi,ni);

      for (kk1 = 0; kk1 < 3; kk1++)
	Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
      for (kk1=0; kk1 < 3; kk1++)
	AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
      AiDjnj = scalProd(AiDj, nj);
      for (kk1=0; kk1 < 3; kk1++)
	{
	  VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
	}
      //for (kk1=0; kk1 < 3; kk1++)
      //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
      ragg = calc_norm(VV);
      for(k=0;k<3;k++)
	{
	  VV[k] /= ragg;
	  Tj[k] = D*0.5*VV[k];

	}
    }

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k];
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
#if 0
      ragg = scalProd(TnCi,ni);
      for (k=0;k<3;k++)
	Ai[k] = Ci[k] + ragg*ni[k];
#endif
#else

      for (k=0; k < 3; k++)
	dsc[k] = Ci[k] - Dj[k]; 
      sp = scalProd(dsc, nj);
      for (k=0; k < 3; k++)
	R[1][k] = dsc[k] - sp*nj[k];
#endif  
  //printf("scalProd=%.15G\n", scalProd(R[1],R[0]));
  vectProdVec(R[0], R[1], u);
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
double test_overlap_parall_cyl(double *Ci, double *ni, double *Dj, double *nj, double Li, double Diami, double Diamj)
{
  int kk;
  double DjCi[3], DjCini, Ui[3], DjUi[3], normDjUi;
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  //normDjCi = calc_norm(DjCi);
  DjCini = scalProd(DjCi,ni);

  for (kk=0; kk < 3; kk++)
    {
      Ui[kk] = Ci[kk] + DjCini*ni[kk];
      DjUi[kk] = Dj[kk] - Ui[kk];
    }
  normDjUi = calc_norm(DjUi);

  if (normDjUi <= 0.5*(Diamj+Diami) && fabs(DjCini) <= Li*0.5)
    return -1.0;
  else
    return 1.0;
}

double rimdiskone_hybrid(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], fallback, solset;
#ifdef MC_QUART_VERBOSE
  static long int numfb=0;
#endif
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4], uy[3];
  double dsc[3], c0, c1, c2, c3, c02, c12, c22, coeffEr[6], c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2];
  double Cip[2][3], nip[2][3], norm, Rl[3][3];
  double nip02,nip12,nip22,Cip02,Cip12,Cip22, temp;
  double omnip02, omnip12, omnip22;
  //double nip04,nip14,nip24,nip03,nip13,nip23;

  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2, nip1nip2, nip0nip2, nip0nip1; 
  /* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  D2 = D*0.5; 
  D2sq = Sqr(D2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  for (kk1=0; kk1 < 3; kk1++)
    uy[kk1]=Dj[kk1];
  versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 1); 

  //versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[0][kk1] = 0;
      //Aip[kk1] = 0;
      Cip[0][kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[0][kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[0][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip[0]);
  nip0 = nip[0][0]/norm;
  nip1 = nip[0][1]/norm;
  nip2 = nip[0][2]/norm;
  Cip0 = Cip[0][0];
  Cip1 = Cip[0][1];
  Cip2 = Cip[0][2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
#if 0
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
#endif
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);
#if 1
  /* with some simplifications we save a bunch of FLOPS... */
  omnip02 = 1.0 - nip02;
  omnip12 = 1.0 - nip12;
  omnip22 = 1.0 - nip22;
  nip1nip2 = nip1*nip2;
  nip0nip2 = nip0*nip2;
  nip0nip1 = nip0*nip1;

  coeffEr[0] = omnip12;
  coeffEr[1] = omnip22;
  coeffEr[2] = -2.0*nip1nip2;  
  coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
								    Cip1*Cip2*nip1nip2) - D2sq;
  coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
  coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  
#elif 0
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
				    2*nip1*nip23);

  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;

  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
  /* check ellipse */
#if 0
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  double cq[3], x, lam, solq[2], p[3];
  int numsol;
  lam = -Cip0/nip0;
  x = Cip1 + lam*nip1;
  cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
  cq[1] = coeffEr[5]+coeffEr[2]*x;
  cq[2] = coeffEr[1];
  solve_quadratic(cq, &numsol, solq);
  p[0] = 0.0;
  p[1] = x;
  p[2] = solq[0];
  if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
    {
      printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
      printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
    }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2sq;
  coeffEr[1] *= D2sq; 
  coeffEr[2] *= D2sq;
  coeffEr[4] *= D2;
  coeffEr[5] *= D2;
  //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
#ifndef MC_EXCHG_QUART_SOL
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
  coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
  coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
  coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
  printf("coeff[C]=");
  for (kk1=0; kk1 < 5; kk1++)
    printf(" %.15G ", coeff[kk1]);
  printf("\n");
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
      if (temp==0.0) 
	{
	  fallback=1;
	}
    }
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2*solec[kk1][0];
      solarr[0][kk1][2] = D2*solec[kk1][1];
      if (test_for_fallback(solarr[0][kk1], Cip[0], nip[0], D2, &(diff[0][kk1])))
	{
	  fallback=1;
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
    }
#if 0
  if (tinyimagGBL)
    {
      //printf("qui\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
#endif
  solset=0;
  if (fallback)
    {
#if defined(MC_QUART_VERBOSE) && 0
      printf("numsol=%d,", numsol[0]);
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  printf("sol=%.16G ", solqua[kk1]);
	}
      printf("\n");
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  temp = c4 + c2*solqua[kk1];
	  printf("temp[%d]=%.16G\n", kk1, temp);
	}
      printf("c2=%.16G c4=%.16G c5=%.16G\n", c2, c4, c5);
      store_bump(iGbl,jGbl);
#endif

#if 0
      for (kk1=0; kk1 < 3; kk1++)
	uy[kk1] = D2*Rl[1][kk1];
      versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 2); 
#else
      rotate_axes_on_plane(Rl);
#endif
      for (kk1=0; kk1 < 3; kk1++)
	{
	  nip[1][kk1] = 0;
	  //Aip[kk1] = 0;
	  Cip[1][kk1] = 0;
	  for (kk2=0; kk2 < 3; kk2++)
	    {
	      nip[1][kk1] += Rl[kk1][kk2]*ni[kk2];
	      Cip[1][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	      //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	    } 
	}
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
      norm = calc_norm(nip[1]);
      nip0 = nip[1][0]/norm;
      nip1 = nip[1][1]/norm;
      nip2 = nip[1][2]/norm;
      Cip0 = Cip[1][0];
      Cip1 = Cip[1][1];
      Cip2 = Cip[1][2];
      nip02=Sqr(nip0);
      nip12=Sqr(nip1);
      nip22=Sqr(nip2);
      //nip04=Sqr(nip02);
      //nip14=Sqr(nip12);
      //nip24=Sqr(nip22);
      //nip03=nip02*nip0;
      //nip13=nip12*nip1;
      //nip23=nip22*nip2;
      Cip02=Sqr(Cip0);
      Cip12=Sqr(Cip1);
      Cip22=Sqr(Cip2);   
#if 1
      omnip02 = 1.0 - nip02;
      omnip12 = 1.0 - nip12;
      omnip22 = 1.0 - nip22;
      nip1nip2 = nip1*nip2;
      nip0nip2 = nip0*nip2;
      nip0nip1 = nip0*nip1;
      coeffEr[0] = omnip12;
      coeffEr[1] = omnip22;
      coeffEr[2] = -2.0*nip1nip2;  
      coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
									Cip1*Cip2*nip1nip2) - D2sq;
      coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
      coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  

#elif 0
      coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
	nip12*nip22;
      coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
	nip12*nip22 + nip24;
      coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
	2*nip1*nip23;
      coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
	2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
	2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
	4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
	2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
	Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
	2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
      coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
	4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
	2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;
      coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
	4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
	- 2*Cip2*nip24;
#else
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

      coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
      coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
      coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
					2*nip1*nip23);

      coeffEr[3] = 
	(- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
	 Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
	 + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
	 Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
	+ Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
	2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
	+ Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

      coeffEr[4] =
	(4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
	 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
	 - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
	- 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;

      coeffEr[5] = 
	(4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
	 - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
      /* check ellipse */

      /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
      coeffEr[0] *= D2sq;
      coeffEr[1] *= D2sq; 
      coeffEr[2] *= D2sq;
      coeffEr[4] *= D2;
      coeffEr[5] *= D2;
      //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
      c0 = coeffEr[0];
      c1 = coeffEr[1];
      c2 = coeffEr[2];
      c3 = coeffEr[3];
      c4 = coeffEr[4];
      c5 = coeffEr[5];
      c02 = Sqr(c0);
      c12 = Sqr(c1);
      c22 = Sqr(c2);
      c32 = Sqr(c3);
      c42 = Sqr(c4);
      c52 = Sqr(c5);

#if 0
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));
#ifdef MC_QUART_VERBOSE
      printf("falling back [#%ld] type=%d numsol=%d %d\n", numfb++,fallback, numsol[0], numsol[1]);
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{

#if 0 
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif
	  temp = c4 + c2*solqua[kk1];
	  //printf("tempnew=%.16G\n", temp);
	  if (temp==0)
	    {
	      printf("[WARNING] temp is 0 in fallback hybrid numsol=%d %d\n", numsol[0], numsol[1]);
	    }
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2*solec[kk1][0];
	  solarr[1][kk1][2] = D2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
	  test_for_fallback(solarr[1][kk1], Cip[1], nip[1], D2, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  
	}
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	    //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
	}
    }
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[solset][kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallback(solarr[solset][kk1], Cip[solset], nip[solset], D2, &tmp)) 
	{
	  printf("# %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip[solset], nip[solset]));

	  if (kk1 < numsol[1-solset])
	    {
	      printf("DISCARDED SOLSET [%d]\n", 1-solset);
	      printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[1-solset][kk1]));
	      printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip[1-solset], nip[1-solset]));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[solset][1], Cip[solset][2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
#endif
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProd(dsc, nip[solset]);
      if (fabs(sp) < L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
double rimdiskone_solvxy(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], fallback, solset;
#ifdef MC_QUART_VERBOSE
  static long int numfb=0;
#endif
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4];
  double dsc[3], c0, c1, c2, c3, c02, c12, c22, coeffEr[6], c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2];
  double Cip[3], nip[3], norm, Rl[3][3];
  double nip02,nip12,nip22,nip03,nip13,nip23,nip04,nip14,nip24,Cip02,Cip12,Cip22, temp;
  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2; 
/* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  D2 = D*0.5; 
  D2sq = Sqr(D2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[kk1] = 0;
      //Aip[kk1] = 0;
      Cip[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip);
  nip0 = nip[0]/norm;
  nip1 = nip[1]/norm;
  nip2 = nip[2]/norm;
  Cip0 = Cip[0];
  Cip1 = Cip[1];
  Cip2 = Cip[2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);   
#if 1
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
    2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
    2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
    4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
    2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
    Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
    2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
  coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
    4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
    2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
  coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
    4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
    - 2*Cip2*nip24;
#else
 /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

  coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
  coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
  coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23);
 
  coeffEr[3] = 
    (- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
     Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
     + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
     Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
    + Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
    2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
    + Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
 
  coeffEr[4] =
    (4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
     - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
     - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
    - 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
    2*Cip2*nip1*nip23;
 
  coeffEr[5] = 
    (4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
     - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
    4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;
 
#endif
  /* check ellipse */
#if 0
  {
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
    double cq[3], x, lam, solq[2], p[3];
    int numsol;
    lam = -Cip0/nip0;
    x = Cip1 + lam*nip1;
    cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
    cq[1] = coeffEr[5]+coeffEr[2]*x;
    cq[2] = coeffEr[1];
    solve_quadratic(cq, &numsol, solq);
    p[0] = 0.0;
    p[1] = x;
    p[2] = solq[0];
    if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
      {
	printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      	printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
	printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
      }
  }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2sq;
  coeffEr[1] *= D2sq; 
  coeffEr[2] *= D2sq;
  coeffEr[4] *= D2;
  coeffEr[5] *= D2;
  //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
  //xC=yC=0;
#ifndef MC_EXCHG_QUART_SOL
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
  coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
  coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
  coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
#if 0
  if (numsol==1)
    {
      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
      printf("sol=%.15G\n", solqua[0]);
      printf("BOH\n");
    }
#endif
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
#if 0
      c0l = c0;
      c1l = c1;
      c2l = c2;
      c3l = c3;
      c4l = c4;
      c5l = c5;
      solqual=solqua[kk1];
      templ = c4l + c2l*solqual;
      solec[kk1][0]=((double)((-c0l - c3l - c5l*solqual + (c0l - c1l)*Sqr(solqual))/templ));
      solec[kk1][1] = solqua[kk1];
#else
#ifndef MC_EXCHG_QUART_SOL
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
#else
      temp = c5 + c2*solqua[kk1];
      solec[kk1][0] = solqua[kk1];
      solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif     
 //printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
#if 0
      if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	{
  	  printf("solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	  printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
	  printf("numeratore=%.16G\n", -c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]));
	}
#endif
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
#if 0
      if (test_solution_xy(solec[kk1], &(diffxy[0][kk1])))
	  fallback=1;
#endif
      if (temp==0.0) 
	{
	  printf("[WARNING] temp is 0 in hybrid numsol=%d\n", numsol[0]);
	  fallback=1;
	}
#if 0
      if (fabs(temp) < FALLBACK_THR)
	{
	  fallback = 1;
	  break;
	}
#endif
    }
#if 0
      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
	     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
	     coeff[1]*solqua[kk1]+coeff[0]);
      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
#if 0
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      solve_quartic(coeff, &numsol, solqua);
      for (kk1=0; kk1 < numsol; kk1++)
	{
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	}
    }
#endif
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2*solec[kk1][0];
      solarr[0][kk1][2] = D2*solec[kk1][1];
#if 1
      if (test_for_fallback(solarr[0][kk1], Cip, nip, D2, &(diff[0][kk1])))
	{
	  fallback=1;
#if 0
	  if (numsol==4)
	    {
	      printf("%d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		     coeff[1]*solqua[kk1]+coeff[0]);
	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[0][kk1]);
	      printf(">>>> =====================\n");
	    }
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
#endif
    }
#if 0
  if (tinyimagGBL)
    {
      //printf("BOHHHH\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
#endif

  solset=0;
#if 1
  if (fallback)
    {
#if defined(MC_QUART_VERBOSE) && 0
      printf("numsol=%d,", numsol[0]);
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  printf("sol=%.16G ", solqua[kk1]);
	}
      printf("\n");
      printf("temp=%.16G\n", temp);
      store_bump(iGbl,jGbl);
#endif

#ifndef MC_EXCHG_QUART_SOL
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#else
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#endif
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));
#ifdef MC_QUART_VERBOSE
      printf("falling back [#%ld] type=%d numsol=%d %d\n", numfb++,fallback, numsol[0], numsol[1]);
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0
    	  c0l = c0;
	  c1l = c1;
	  c2l = c2;
	  c3l = c3;
	  c4l = c4;
	  c5l = c5;
	  solqual=solqua[kk1];
	  templ = c5l + c2l*solqual;
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = ((double)((-c1l - c3l - c4l*solqual + (c1l - c0l)*Sqr(solqual))/templ)); 
#else
#ifndef MC_EXCHG_QUART_SOL
	  temp = c5 + c2*solqua[kk1];
	  //printf("tempnew=%.16G\n", temp);
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#else
    	  temp = c4 + c2*solqua[kk1];
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
#endif
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2*solec[kk1][0];
	  solarr[1][kk1][2] = D2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
#if 1
	  test_for_fallback(solarr[1][kk1], Cip, nip, D2, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  

#endif
#if 0
  	  if (numsol==4)
  	    {
  	      printf("FALLBACK %d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
  	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
  	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
  	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
  	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
  	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
  		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
  		     coeff[1]*solqua[kk1]+coeff[0]);
  	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[1][kk1]);
	      printf(">>>> =====================\n");
	  }
#endif
	
	}
#if 1
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	  //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
#if 0
	  if (fallback==3 && solset != 1)
	    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
	}
#endif
    }
#endif
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
#if 0
  construct_inner_points(solarr, Ci, ni, Dj, nj, D);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
#if 0
      printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
      printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 0
      if (calc_normsq(solarr[kk1])-Sqr(D2) > NEWT_THR)
	{
	  newt2Dquartic(coeffEr, solarr[kk1], D2);
	}
#endif
#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallback(solarr[solset][kk1], Cip, nip, D2, &tmp)) 
	{
	  printf("# %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip, nip));
	  
	  if (kk1 < numsol[1-solset])
	    {
	      printf("DISCARDED SOLSET [%d]\n", 1-solset);
    	      printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[1-solset][kk1]));
    	      printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip, nip));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[1], Cip[2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.14G\n", scalProd(ni,nj));
	  printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
#endif
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = scalProd(dsc, nip);
      if (fabs(sp) < L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
double rimdiskone(double D, double L, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
#ifdef MC_QUART_HYBRID
  return rimdiskone_hybrid(D, L, Ci, ni, Dj, nj, DjCini);
#else
  return rimdiskone_solvxy(D, L, Ci, ni, Dj, nj, DjCini);
#endif

}
double rimdisk(double D, double L, double Ci[3], double ni[3], double Di[2][3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, kk, j2, k2, ignore[2];
  double DjUini, DjTmp[2][3], DjCi[3], DjUi[3], niTmp[3], njTmp[3], perpdist[2];
  double normDjUi, normDjCi, DjCini, Ui[3], CiTmp[3];


  for (j1=0; j1 < 2; j1++)
    {
      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	    }
	}
      
      for (k2=0; k2 < 2; k2++)
	{
  	  perpdist[k2]=perpcomp(Dj[k2], Ci, ni);
	  ignore[k2] = 0;
	}

      if (perpdist[0] < perpdist[1])
	{ 
    	  ignore[1] = 1;
	}
      else
	{
	  ignore[0] = 1;
	}

      for (j2=0; j2 < 2; j2++)
	{
	  if (ignore[j2])
	    continue;
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);
#if 0
	  if (dostorebump)
	    {
	      printf("normDjUi=%.15G DjUini=%.15G\n", normDjUi, DjUini);
	      printf("Ci=%f %f %f Dj=%f %f %f\n", Ci[0], Ci[1], Ci[2], Dj[0], Dj[1], Dj[2]);
	      printf("DjUi=%.15G %.15G %.15G\n", DjUi[0], DjUi[1], DjUi[2]); 
	      printf("Uj=%.15G %.15G %.15G\n", Ui[0], Ui[1], Ui[2]); 
	      printf("nj=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
	      printf("DjCini= %.15G\n", DjCini);
	    }
#endif 
	  if (normDjUi > D)
            {
              printf("[C] qui\n");
              continue;
            }
	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Ui| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

*/
#ifndef MC_IBARRA_SIMPLER
	  /* se sono paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) == 0)
	  //if (1.0-fabs(scalProd(ni,nj)) < 1.0E-8)
	    {
	      if (normDjUi <= D && fabs(DjCini) <= L*0.5)
		return -1;
	      else
		continue;
	    }
#else
	  /* se sono paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) == 0)
	  //if (1.0-fabs(scalProd(ni,nj)) < 3E-16)
	    {
	      if (normDjUi <= D && fabs(DjCini) <= L*0.5)
		return -1;
	      else
		continue;
	    }
#endif
	  if (normDjUi < D*0.5 && fabs(DjCini) > L*0.5)
	    continue;

	  if (normDjUi < D*0.5 && fabs(DjCini) <= L*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }

	  if (rimdiskone(D, L, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
            return -1;
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	    }
	}

    }
  return 0;
}

double calcDistNegHCopt(void)
{
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int kk;
  double Ci[3], Cj[3], L, D, Di[2][3], Dj[2][3], ni[3], nj[3], ret;

#if 0
  if (typesArr[typeOfPart[i]].sax[0]!=typesArr[typeOfPart[j]].sax[0]
      || typesArr[typeOfPart[i]].sax[1] != typesArr[typeOfPart[j]].sax[1])
    return calcDistNegHCoptdiff(i, j, shift, retchk);
#endif

  for (kk=0; kk < 3; kk++)
    {
      ni[kk] = HardCyl[0].n[kk];
      nj[kk] = HardCyl[1].n[kk];
    }
  Ci[0] = HardCyl[0].rcm[0];
  Ci[1] = HardCyl[0].rcm[1];
  Ci[2] = HardCyl[0].rcm[2];
  Cj[0] = HardCyl[1].rcm[0];
  Cj[1] = HardCyl[1].rcm[1];
  Cj[2] = HardCyl[1].rcm[2]; 
  L = 2.0*HardCyl[0].sax[2];// elongated along z axis in the body reference frame
  D = 2.0*HardCyl[0].sax[0];
  
  //meshptsGbl = MESH_PTS;
  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*L*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*L*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*L*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*L*nj[kk];
    }
#ifdef MC_HC_SPHERO_OPT
  if ((sphov=check_spherocyl(CiCj, D, L, Di, Ci, ni, Dj, Cj, nj, &rim)) > 0.0)
    return 1;
#endif
#ifdef MC_DEBUG_HCALGO
  iGbl=i;
  jGbl=j;
#endif
  if (L >= D) // prolate
    {
      if ((ret=diskdisk(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
        {
          //printf("[C] diskdisk\n");
          return ret;
        }
      if ((ret=rimrim(D, L, Ci, ni, Cj, nj)) != 0.0)
        {
          //printf("[C] rimrim\n");
          return ret;
        }
      if ((ret=rimdisk(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
        {
          //printf("[C] rimdisk\n");
          return ret;
        }
      //printf("[C] no overlap\n");
     }
  else // oblate
    {
      if ((ret=diskdisk(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimrim(D, L, Ci, ni, Cj, nj)) != 0.0)
	return ret;

      if ((ret=rimdisk(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
	return ret;
    }
  /* case A.2 overlap of rim and disk */
  /* =================================== >>> Part A <<< ========================= */
 /* =================================== >>> Part B <<< ========================= */
  return 1;
}
/* cylinder overlap routines here */
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
double calcDistNegHCdiff(void)
{
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int it, k2;
  double normDjCi, normNSq, ViVj[3], lambdai, lambdaj, Li, Diami, Lj, Diamj; 
  double LiTmp, LjTmp, DiamiTmp, DiamjTmp;
  double sp, Q1, Q2, normPiDi, normPjDj, normN, DiN, DjN, niN[3], njN[3], Djnj;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DjCi[3];
  double PiDi[3], PjDj[3], Ai[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3];
  double Dini, Tjm_perp[3];
  double DjCini;
  double Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp;
  double normCiCj;	
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3];
  int kk, j1, j2;
  // *retchk = 0; 

  // return calcDistNegHCsame(i, j, shift, retchk);
  for (kk=0; kk < 3; kk++)
    {
      ni[kk] = HardCyl[0].n[kk];
      nj[kk] = HardCyl[1].n[kk];
    }
  
  /* qui va cambiato */
  Ci[0] = HardCyl[0].rcm[0];
  Ci[1] = HardCyl[0].rcm[1];
  Ci[2] = HardCyl[0].rcm[2];
  Cj[0] = HardCyl[1].rcm[0];
  Cj[1] = HardCyl[1].rcm[1];
  Cj[2] = HardCyl[1].rcm[2]; 
  Li = 2.0*HardCyl[0].sax[2];
  Diami = 2.0*HardCyl[0].sax[1];
  Lj = 2.0*HardCyl[1].sax[2];
  Diamj = 2.0*HardCyl[1].sax[1];
  //printf("r1=%f %f %f r2=%f %f %f\n",CHROMall[0].rcm[0],CHROMall[0].rcm[1],CHROMall[0].rcm[2],
  //  CHROMall[1].rcm[0],CHROMall[1].rcm[1],CHROMall[1].rcm[2]);
  //printf("Li=%f Di=%f Lj=%f Dj=%f\n", Li, Diami, Li, Diamj);
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }

  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*Li*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*Li*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*Lj*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*Lj*nj[kk];
    }
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= 0.5*(Li+Lj))
	    return -1;
	  else
	    return 1;
	}

      /* parallel disks */
      for (j1=0; j1 < 2; j1++)
	for (j2=j1; j2 < 2; j2++)
	  {
	    sp=0.0;
	    for (kk=0; kk < 3; kk++)
	      {
		VV[kk] = Di[j1][kk]-Dj[j2][kk];
		sp += ni[kk]*VV[kk];
	      }
	    if (sp == 0 && calc_norm(VV) < 0.5*(Diami+Diamj))
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
#ifdef DEBUG_HCMC
	    printf("Di=%f %f %f\n", Di[j1][0], Di[j1][1], Di[j1][2]);
	    printf("Dj=%f %f %f\n", Dj[j2][0], Dj[j2][1], Dj[j2][2]);
	    printf("normPiDi: %f normPjDj=%f\n", normPiDi, normPjDj);
	    printf("0.5*Diami=%f 0.5*Diamj=%f\n", 0.5*Diami, 0.5*Diamj);
#endif
	    if (normPiDi <= 0.5*Diami && normPjDj <= 0.5*Diamj)
	      {
		Q1 = sqrt(Sqr(Diami)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(Diamj)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  /* =================================== >>> Part A <<< ========================= */
  for (j1=0; j1 < 2; j1++)
    {

      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      DiamiTmp = Diami;
	      DiamjTmp = Diamj;
	      LiTmp = Li;
	      LjTmp = Lj;
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	      Diami = Diamj;
	      Diamj = DiamiTmp;
	      Li = Lj;
	      Lj = LiTmp;
	    }
	}
      for (j2=0; j2 < 2; j2++)
	{
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);

	  if (normDjUi > 0.5*(Diami+Diamj))
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Uj| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

	   */
	  if (normDjUi < Diami*0.5 && fabs(DjCini) > Li*0.5)
	    continue;

	  if (normDjUi < Diami*0.5 && fabs(DjCini) <= Li*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }
#if 1
	  find_initial_guess(Ai, Ci, ni, Dj[j2], nj, Diamj);
#else
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
	  for (it = 0; it < MAX_ITERATIONS; it++)
	    {
	      for (kk=0; kk < 3; kk++)
		{
		  AiDj[kk] = Ai[kk] - Dj[j2][kk];
		}
	      AiDjnj = scalProd(AiDj,nj);
	      vectProdVec(AiDj,nj,AiDjnjvec);
	      for (kk=0; kk < 3; kk++)
		VV[kk] =  0.5*Diamj*(AiDj[kk]-AiDjnj*nj[kk])/calc_norm(AiDjnjvec);
	      for (kk=0; kk < 3; kk++)
		{
		  Tjp[kk] = Dj[j2][kk] + VV[kk];
		  Tjm[kk] = Dj[j2][kk] - VV[kk];
		  TjpCi[kk] = Tjp[kk] - Ci[kk];
		  TjmCi[kk] = Tjm[kk] - Ci[kk];
		}
	      TjpCini = scalProd(TjpCi,ni);  
	      TjmCini = scalProd(TjmCi,ni);
	      for (kk=0; kk < 3; kk++)
		{
		  Tjp_perp[kk] = TjpCi[kk]-TjpCini*ni[kk];
		  Tjp_para[kk] = TjpCini*ni[kk];
		  Tjm_perp[kk] = TjmCi[kk]-TjmCini*ni[kk];
		  Tjm_para[kk] = TjmCini*ni[kk];
		} 
	      normTjm_perp = calc_norm(Tjp_perp);
	      for (kk=0; kk < 3; kk++)
		TjOld[kk] = TjNew[kk];
	      if (calc_norm(Tjm_perp) < calc_norm(Tjp_perp))
		{
		  for (kk=0; kk < 3; kk++)
		    TjNew[kk] = Tjm[kk];
		}	  
	      else
		{
		  for (kk=0; kk < 3; kk++)
		    TjNew[kk] = Tjp[kk];
		}

	      for (kk=0; kk < 3; kk++)
		TjNewCi[kk] = TjNew[kk] - Ci[kk];
	      TjNewCini = scalProd(TjNewCi,ni);

#ifdef DEBUG_HCMC
	      printf("j1=%d A it=%d Aiold=%.15G %.15G %.15G\n", j1, it, Ai[0], Ai[1], Ai[2]);
#endif
	      for (kk=0; kk < 3; kk++)
		Ai[kk] = TjNewCini*ni[kk] + Ci[kk]; 
#ifdef DEBUG_HCMC
	      printf("A it=%d Ainew=%.15G %.15G %.15G TjNewCini=%.15G\n", it, Ai[0], Ai[1], Ai[2], TjNewCini);
	      printf("A Ci=%.15G %.15G %.15G\n", Ci[0], Ci[1], Ci[2]);
	      printf("A ni=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
#endif
	      if ( it > 0 && check_convergence(TjOld,TjNew) ) 
		break;
	    }
	  //totitsHC += it;
#ifdef DEBUG_HCMC
	  printf("A #1 number of iterations=%d Tjold=%.15G %.15G %.15G Tjnew=%.15G %.15G %.15G\n",it, 
		 TjOld[0], TjOld[1], TjOld[2], TjNew[0], TjNew[1], TjNew[2]);
#endif
	  if (it >= MAX_ITERATIONS)
	    {
	      printf("MAX ITERATIONS REACHED in A!\n");
	      //*retchk=1;
	      return -1;
	    }
	  if ( (calc_norm(Tjp_para) <= Li*0.5 && calc_norm(Tjp_perp) <= Diami*0.5)||
	       (calc_norm(Tjm_para) <= Li*0.5 && calc_norm(Tjm_perp) <= Diami*0.5) )
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #2 disk-rim\n");
#endif	   
	      return -1;
	    }
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	      Diami = DiamiTmp;
	      Diamj = DiamjTmp;
	      Li = LiTmp;
	      Lj = LjTmp;
	    }
	}

    }
  /* =================================== >>> Part B <<< ========================= */
  //numcallsHC += 4.0; 

  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < 0.5*(Diami+Diamj) && fabs(lambdai) < 0.5*Li && fabs(lambdaj) < 0.5*Lj)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 1;
}

static long iran=0;
double *vector(int n1, int n2)
{
  return (double*)malloc(sizeof(double)*(n2+1));
}
void free_vector(double* p, int n1, int n2)
{
  free(p);
}
void sobseq(int *n, double x[]);
#define SQR(x) ((x)*(x))
double ran1(void)
{
  return drand48();
}
#if defined(QUASIMC) || defined(QFGAUSS)
double sv[10];
int nsv;
#endif
long long int fileoutits, outits;
long long int tot_trials, tt=0, ttini=0;
double rSugar=3.5;
void ranpt(double pt[], double regn[], int n)
{
  //	float ran1(long *idum);
  int j;
#ifdef QUASIMC
  double sv[10];

  if (nsv!=n)
    {
      printf("nsv=%d is different from n=%d\n", nsv, n);
      exit(-1);
    }
  sobseq(&nsv, sv);
  for (j=1;j<=n;j++)
    pt[j]=regn[j]+(regn[n+j]-regn[j])*sv[j];
#else
  for (j=1;j<=n;j++)
    pt[j]=regn[j]+(regn[n+j]-regn[j])*ran1();
#endif
}
void miser(double (*func)(double []), double regn[], int ndim, unsigned long npts,
	double dith, double *ave, double *var)
{
	double *regn_temp;
	unsigned long n,npre,nptl,nptr;
	int j,jb;
	double avel,varl;
	double fracl,fval;
	double rgl,rgm,rgr,s,sigl,siglb,sigr,sigrb;
	double sum,sumb,summ,summ2;
	double *fmaxl,*fmaxr,*fminl,*fminr;
	double *pt,*rmid;

	pt=vector(1,ndim);
	//pt = malloc(sizeof(double)*(ndim+1));
	if (npts < MNBS) {
		summ=summ2=0.0;
		for (n=1;n<=npts;n++) {
			ranpt(pt,regn,ndim);
			fval=(*func)(pt);
			summ += fval;
			summ2 += fval * fval;
		}
		*ave=summ/npts;
		*var=FMAX(TINY,(summ2-summ*summ/npts)/(npts*npts));
		//printf("[npts < MNBS npts=%d] ave=%.15G  var=%.15G\n", npts, *ave, *var);
	}
	else {
		rmid=vector(1,ndim);
		npre=LMAX((unsigned long)(npts*PFAC),MNPT);
		fmaxl=vector(1,ndim);
		fmaxr=vector(1,ndim);
		fminl=vector(1,ndim);
		fminr=vector(1,ndim);
		for (j=1;j<=ndim;j++) {
			iran=(iran*2661+36979) % 175000;
			s=SIGN(dith,(float)(iran-87500));
			rmid[j]=(0.5+s)*regn[j]+(0.5-s)*regn[ndim+j];
			fminl[j]=fminr[j]=BIG;
			fmaxl[j]=fmaxr[j] = -BIG;
		}
		for (n=1;n<=npre;n++) 
		  {
		    //if (n % outits == 0)
		      //printf("step %d/%lld\n", n, tot_trials);
		    ranpt(pt,regn,ndim);
		    fval=(*func)(pt);
	    	    for (j=1;j<=ndim;j++) {
		      if (pt[j]<=rmid[j]) {
					fminl[j]=FMIN(fminl[j],fval);
					fmaxl[j]=FMAX(fmaxl[j],fval);
				}
				else {
					fminr[j]=FMIN(fminr[j],fval);
					fmaxr[j]=FMAX(fmaxr[j],fval);
				}
			}
		}
		sumb=BIG;
		jb=0;
		siglb=sigrb=1.0;
		for (j=1;j<=ndim;j++) {
			if (fmaxl[j] > fminl[j] && fmaxr[j] > fminr[j]) {
				sigl=FMAX(TINY,pow(fmaxl[j]-fminl[j],2.0/3.0));
				sigr=FMAX(TINY,pow(fmaxr[j]-fminr[j],2.0/3.0));
				sum=sigl+sigr;
				if (sum<=sumb) {
					sumb=sum;
					jb=j;
					siglb=sigl;
					sigrb=sigr;
				}
			}
		}
		free_vector(fminr,1,ndim);
		free_vector(fminl,1,ndim);
		free_vector(fmaxr,1,ndim);
		free_vector(fmaxl,1,ndim);
		if (!jb) jb=1+(ndim*iran)/175000;
		rgl=regn[jb];
		rgm=rmid[jb];
		rgr=regn[ndim+jb];
		fracl=fabs((rgm-rgl)/(rgr-rgl));
		nptl=(unsigned long)(MNPT+(npts-npre-2*MNPT)*fracl*siglb
			/(fracl*siglb+(1.0-fracl)*sigrb));
		nptr=npts-npre-nptl;
		regn_temp=vector(1,2*ndim);
		for (j=1;j<=ndim;j++) {
			regn_temp[j]=regn[j];
			regn_temp[ndim+j]=regn[ndim+j];
		}
		regn_temp[ndim+jb]=rmid[jb];
		miser(func,regn_temp,ndim,nptl,dith,&avel,&varl);
		regn_temp[jb]=rmid[jb];
		regn_temp[ndim+jb]=regn[ndim+jb];
		miser(func,regn_temp,ndim,nptr,dith,ave,var);
		free_vector(regn_temp,1,2*ndim);
		*ave=fracl*avel+(1-fracl)*(*ave);
		*var=fracl*fracl*varl+(1-fracl)*(1-fracl)*(*var);
		//printf("[npts > MNBS npts=%d] ave=%.15G  var=%.15G\n", npts, *ave, *var);
		free_vector(rmid,1,ndim);
	}
	free_vector(pt,1,ndim);
}
#undef MNPT
#undef MNBS
#undef TINY
#undef BIG
#undef PFAC
#undef NRANSI
/* =================== END MISER ================= */
#if 0
/* NOTA: la routine vegas fa un importance sampling e questo non 
   dovrebbe essere molto utile nel caso presente poiché l'integrazione
   su \phi_{12 }e \theta_{12} viene fatta con la quadratura di gauss */
/* ================== BEGIN VEGAS ================ */
#define ALPH 1.5
#define NDMX 50
#define MXDIM 10
#define TINY 1.0e-30

long idum;
double ran2(long *idum)
{
  return drand48();
}

void rebin(double rc, int nd, double r[], double xin[], double xi[]);
void vegas(double regn[], int ndim, double (*fxn)(double [], double), int init,
	unsigned long ncall, int itmx, int nprn, double *tgral, double *sd,
	double *chi2a)
{
	//double ran2(long *idum);
	static int i,it,j,k,mds,nd,ndo,ng,npg,ia[MXDIM+1],kg[MXDIM+1];
	static double calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo;
	static double d[NDMX+1][MXDIM+1],di[NDMX+1][MXDIM+1],dt[MXDIM+1],
		dx[MXDIM+1], r[NDMX+1],x[MXDIM+1],xi[MXDIM+1][NDMX+1],xin[NDMX+1];
	static double schi,si,swgt;

	if (init <= 0) {
		mds=ndo=1;
		for (j=1;j<=ndim;j++) xi[j][1]=1.0;
	}
	if (init <= 1) si=swgt=schi=0.0;
	if (init <= 2) {
		nd=NDMX;
		ng=1;
		if (mds) {
			ng=(int)pow(ncall/2.0+0.25,1.0/ndim);
			mds=1;
			if ((2*ng-NDMX) >= 0) {
				mds = -1;
				npg=ng/NDMX+1;
				nd=ng/npg;
				ng=npg*nd;
			}
		}
		for (k=1,i=1;i<=ndim;i++) k *= ng;
		npg=IMAX(ncall/k,2);
		calls=(float)npg * (float)k;
		dxg=1.0/ng;
		for (dv2g=1,i=1;i<=ndim;i++) dv2g *= dxg;
		dv2g=SQR(calls*dv2g)/npg/npg/(npg-1.0);
		xnd=nd;
		dxg *= xnd;
		xjac=1.0/calls;
		for (j=1;j<=ndim;j++) {
			dx[j]=regn[j+ndim]-regn[j];
			xjac *= dx[j];
		}
		if (nd != ndo) {
			for (i=1;i<=IMAX(nd,ndo);i++) r[i]=1.0;
			for (j=1;j<=ndim;j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
			ndo=nd;
		}
		if (nprn >= 0) {
			printf("%s:  ndim= %3d  ncall= %8.0f\n",
				" Input parameters for vegas",ndim,calls);
			printf("%28s  it=%5d  itmx=%5d\n"," ",it,itmx);
			printf("%28s  nprn=%3d  ALPH=%5.2f\n"," ",nprn,ALPH);
			printf("%28s  mds=%3d  nd=%4d\n"," ",mds,nd);
			for (j=1;j<=ndim;j++) {
				printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
					" ",j,regn[j],j,regn[j+ndim]);
			}
		}
	}
	for (it=1;it<=itmx;it++) {
		ti=tsi=0.0;
		for (j=1;j<=ndim;j++) {
			kg[j]=1;
			for (i=1;i<=nd;i++) d[i][j]=di[i][j]=0.0;
		}
		for (;;) {
			fb=f2b=0.0;
			for (k=1;k<=npg;k++) {
				wgt=xjac;
				for (j=1;j<=ndim;j++) {
					xn=(kg[j]-ran2(&idum))*dxg+1.0;
					ia[j]=IMAX(IMIN((int)(xn),NDMX),1);
					if (ia[j] > 1) {
						xo=xi[j][ia[j]]-xi[j][ia[j]-1];
						rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
					} else {
						xo=xi[j][ia[j]];
						rc=(xn-ia[j])*xo;
					}
					x[j]=regn[j]+rc*dx[j];
					wgt *= xo*xnd;
				}
				f=wgt*(*fxn)(x,wgt);
				f2=f*f;
				fb += f;
				f2b += f2;
				for (j=1;j<=ndim;j++) {
					di[ia[j]][j] += f;
					if (mds >= 0) d[ia[j]][j] += f2;
				}
			}
			f2b=sqrt(f2b*npg);
			f2b=(f2b-fb)*(f2b+fb);
			if (f2b <= 0.0) f2b=TINY;
			ti += fb;
			tsi += f2b;
			if (mds < 0) {
				for (j=1;j<=ndim;j++) d[ia[j]][j] += f2b;
			}
			for (k=ndim;k>=1;k--) {
				kg[k] %= ng;
				if (++kg[k] != 1) break;
			}
			if (k < 1) break;
		}
		tsi *= dv2g;
		wgt=1.0/tsi;
		si += wgt*ti;
		schi += wgt*ti*ti;
		swgt += wgt;
		*tgral=si/swgt;
		*chi2a=(schi-si*(*tgral))/(it-0.9999);
		if (*chi2a < 0.0) *chi2a = 0.0;
		*sd=sqrt(1.0/swgt);
		tsi=sqrt(tsi);
		if (nprn >= 0) {
			printf("%s %3d : integral = %14.7g +/-  %9.2g\n",
				" iteration no.",it,ti,tsi);
			printf("%s integral =%14.7g+/-%9.2g chi**2/IT n = %9.2g\n",
				" all iterations:  ",*tgral,*sd,*chi2a);
			if (nprn) {
				for (j=1;j<=ndim;j++) {
					printf(" DATA FOR axis  %2d\n",j);
					printf("%6s%13s%11s%13s%11s%13s\n",
						"X","delta i","X","delta i","X","delta i");
					for (i=1+nprn/2;i<=nd;i += nprn+2) {
						printf("%8.5f%12.4g%12.5f%12.4g%12.5f%12.4g\n",
							xi[j][i],di[i][j],xi[j][i+1],
							di[i+1][j],xi[j][i+2],di[i+2][j]);
					}
				}
			}
		}
		for (j=1;j<=ndim;j++) {
			xo=d[1][j];
			xn=d[2][j];
			d[1][j]=(xo+xn)/2.0;
			dt[j]=d[1][j];
			for (i=2;i<nd;i++) {
				rc=xo+xn;
				xo=xn;
				xn=d[i+1][j];
				d[i][j] = (rc+xn)/3.0;
				dt[j] += d[i][j];
			}
			d[nd][j]=(xo+xn)/2.0;
			dt[j] += d[nd][j];
		}
		for (j=1;j<=ndim;j++) {
			rc=0.0;
			for (i=1;i<=nd;i++) {
				if (d[i][j] < TINY) d[i][j]=TINY;
				r[i]=pow((1.0-d[i][j]/dt[j])/
					(log(dt[j])-log(d[i][j])),ALPH);
				rc += r[i];
			}
			rebin(rc/xnd,nd,r,xin,xi[j]);
		}
	}
}

void rebin(double rc, int nd, double r[], double xin[], double xi[])
{
	int i,k=0;
	double dr=0.0,xn=0.0,xo=0.0;

	for (i=1;i<nd;i++) {
		while (rc > dr)
			dr += r[++k];
		if (k > 1) xo=xi[k-1];
		xn=xi[k];
		dr -= rc;
		xin[i]=xn-(xn-xo)*dr/r[k];
	}
	for (i=1;i<nd;i++) xi[i]=xin[i];
	xi[nd]=1.0;
}
#undef ALPH
#undef NDMX
#undef MXDIM
#undef TINY
#undef NRANSI
/* ================== END VEGAS ================== */ 
#endif
char dummy1[32], dummy2[32], atname[32], nbname[8];
int type, nat, atnum, nbnum, len;
double L, rx, ry, rz, alpha, dfons_sinth_max, fons_sinth_max, ROMBTOL, Lx, Ly, Lz;
const double thetapts=100000;
#ifdef GAUSS
double gammln(double xx)
  /*Returns the value ln[Γ(xx)] for xx > 0.*/
{
  /*Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.*/
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5}; 
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp); ser=1.000000000190015;
  for (j=0;j<=5;j++) 
    ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x);
}

#if 0
#define EPSJAC 3.0e-14 /*Increase EPS if you don’t have this precision */ 
#define MAXITJAC 10 
void gaujac(double x[], double w[], int n, double alf, double bet)
/*Given alf and bet, the parameters α and β of the Jacobi polynomials, this routine returns arrays x[1..n] and w[1..n] containing the abscissas and weights of the n-point Gauss-Jacobi quadrature formula. The largest abscissa is returned in x[1], the smallest in x[n].*/
{
  double gammln(double xx);
  /*void nrerror(char error_text[]); */
  int i,its,j;
  double alfbet,an,bn,r1,r2,r3;
  double a,b,c,p1,p2,p3,pp,temp,z,z1;
  for (i=1;i<=n;i++) { 
    if (i == 1) { 
      an=alf/n;
      /* High precision is a good idea for this rou- tine.
	 Loop over the desired roots. Initial guess for the largest root.*/
      bn=bet/n; 
      r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n); 
      r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn; 
      z=1.0-r1/r2;
    } 
    else if (i == 2) 
      { /* Initial guess for the second largest root.*/
	r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf)); 
	r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n; 
	r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
	z -= (1.0-z)*r1*r2*r3;
      } 
    else if (i == 3) 
      { /* Initial guess for the third largest root.*/
	r1=(1.67+0.28*alf)/(1.0+0.37*alf); 
	r2=1.0+0.22*(n-8.0)/n; 
	r3=1.0+8.0*bet/((6.28+bet)*n*n);
	z -= (x[1]-z)*r1*r2*r3;
      } 
    else if (i == n-1) 
      { /* Initial guess for the second smallest root. */ 
	r1=(1.0+0.235*bet)/(0.766+0.119*bet); 
	r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0))); 
	r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
	z += (z-x[n-3])*r1*r2*r3;
      } 
    else if (i == n) 
      { /* Initial guess for the smallest root.*/
	r1=(1.0+0.37*bet)/(1.67+0.28*bet); 
	r2=1.0/(1.0+0.22*(n-8.0)/n); 
	r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n)); 
	z += (z-x[n-2])*r1*r2*r3;
      } 
    else 
      { /* Initial guess for the other roots.*/
       	z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
      }
    alfbet=alf+bet;
    for (its=1;its<=MAXITJAC;its++) 
      {
	temp=2.0+alfbet; p1=(alf-bet+temp*z)/2.0; p2=1.0;
	for (j=2;j<=n;j++) {
	  p3=p2;
	  p2=p1;
	  temp=2*j+alfbet;
	  a=2*j*(j+alfbet)*(temp-2.0);
	  b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z); 
	  c=2.0*(j-1+alf)*(j-1+bet)*temp;
	  p1=(b*p2-c*p3)/a;
	} 
	pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z)); 
	/* p1 is now the desired Jacobi polynomial. We next compute pp, its derivative, 
	   by a standard relation involving also p2, the polynomial of one lower order. */
	z1=z;
	z=z1-p1/pp; /*Newton’s formula.*/
	if (fabs(z-z1) <= EPSJAC) break;
      }
    if (its > MAXITJAC)
      {
	printf("too many iterations in gaujac"); 
	exit(-1);
      }
    x[i]=z; /* Store the root and the weight.*/
    w[i]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
	     gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
  }
}
#endif
int ntheta, nphi, ngamma;
double *xtheta, *xphi, *wtheta, *wphi, *xgamma, *xphi, *wgamma;
#define EPSGAULEG 3.0e-11 
/*EPS is the relative precision.*/
void gauleg(double x1, double x2, double x[], double w[], int n)
  /*Given the lower and upper limits of integration x1 and x2, and given n, this routine returns arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss- Legendre n-point quadrature formula.*/
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  /*
     High precision is a good idea for this rou- tine.
     The roots are symmetric in the interval, so we only have to find half of them.
     Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
     Loop over the desired roots. */
  for (i=1;i<=m;i++) 
    {
      z=cos(3.141592654*(i-0.25)/(n+0.5));
      /*Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.*/
      do {
	p1=1.0;
	p2=0.0;
	for (j=1;j<=n;j++) {
	  p3=p2;
	  p2=p1; p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	}
	/*
	   p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a 
	   standard relation involving also p2, the polynomial of one lower order. 
	 */	   
	pp=n*(z*p1-p2)/(z*z-1.0);
	z1=z;
	z=z1-p1/pp;
      } while (fabs(z-z1) > EPSGAULEG); 
      x[i]=xm-xl*z;
      x[n+1-i]=xm+xl*z; 
      w[i]=2.0*xl/((1.0-z*z)*pp*pp); 
      w[n+1-i]=w[i];
      /*Scale the root to the desired interval, and put in its symmetric counterpart. 
	Compute the weight and its symmetric counterpart.*/
    }
}
#if 0
int j;
float xr,xm,dx,s;
static float x[]={0.0,0.1488743389,0.4333953941,
0.6794095682,0.8650633666,0.9739065285}; static float w[]={0.0,0.2955242247,0.2692667193, 0.2190863625,0.1494513491,0.0666713443};
The abscissas and weights. First value of each array not used.
xm=0.5*(b+a); xr=0.5*(b-a);
s=0;
for (j=1;j<=5;j++) {
dx=xr*x[j];
Will be twice the average value of the function, since the ten weights (five numbers above each used twice) sum to 2.
s += w[j]*((*func)(xm+dx)+(*func)(xm-dx)); }
return s *= xr;
#endif
double qgaus(double (*func)(double, int), double a, double b, double *x, double *w, int np)
{
#if 0
  static const double x[]={0.1488743389816312,0.4333953941292472,
    0.6794095682990244,0.8650633666889845,0.9739065285171717};
  static const double w[]={0.2955242247147529,0.2692667193099963,
    0.2190863625159821,0.1494513491505806,0.0666713443086881};
#endif
  int j;
  double s;
  s=0.0;
  for (j=1;j<=np;j++) 
    {
      s += w[j]*func(x[j],j);
    }
  return s;
}
#endif
#define JMAX 15
#define JMAXP (JMAX+1) 
#define KROMB 5
int polinterr=0;
/*Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate; JMAX limits the total number of steps; K is the number of points used in the extrapolation.*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y,
 * and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xai) = yai, 
 * i = 1, . . . , n, then the returned value y = P(x).*/
{ 
  int i,m,ns=1; 
  double den,dif,dift,ho,hp,w;
  double c[KROMB+1], d[KROMB+1];
#if 0
  for (i=0; i < n; i++)
    {
      xa[i+1] = xain[i];
      ya[i+1] = yain[i];
    }
#endif
  dif=fabs(x-xa[1]); 
  //c=vector(n); 
  //d=vector(n); 
  for (i=1;i<=n;i++) 
    { 
      /* Here we find the index ns of the closest table entry,*/
      if ( (dift=fabs(x-xa[i])) < dif) 
	{ 
	  ns=i; 
	  dif=dift;
	} 
      c[i]=ya[i];
      /* and initialize the tableau of c s and d s.*/
      d[i]=ya[i]; 
    } 
  *y=ya[ns--];
  /* This is the initial approximation to y.*/
  for (m=1;m<n;m++) 
    { 
      /* For each column of the tableau,*/
      for (i=1;i<=n-m;i++)
	{
	  /* we loop over the current c s and d s and update them.*/
	  ho=xa[i]-x; 
	  hp=xa[i+m]-x; 
	  w=c[i+1]-d[i]; 
	  if ( (den=ho-hp) == 0.0) 
	    {
	      polinterr=1;
	      printf("error in routinr polint\n");
	      return ;
	      //nrerror("Error in routine polint"); 
	    }
	  /* This error can occur only if two input xa s are (to within roundoff)*/
	  den=w/den; d[i]=hp*den; 
	  /*Here the c s and d s are updated. */
	  c[i]=ho*den; 
	} 
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--])); 
      /* After each column in the tableau is completed, we decide which correction, 
       * c or d, we want to add to our accumulating value of y, i.e., which path to take through the tableau 
       * forking up or down. We do this in such a way as to take the most  straight line  route through the 
       * tableau to its apex, updating ns accordingly to keep track of where we are. 
       * This route keeps the partial approximations centered (insofar as possible) on the target x. 
       * The last dy added is thus the error indication. */
    } 
  //free_vector(d); 
  //free_vector(c); 
}
#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
/*This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
as a pointer to the function to be integrated between limits a and b, also input. When called with
n=1, the routine returns the crudest estimate of b f (x)dx. Subsequent calls with n=2,3,... a
(in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.*/
{
  double x,tnm,sum,del; 
  static double s;
  int it,j;
  if (n == 1) 
    {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } 
  else
    {
      for (it=1,j=1;j<n-1;j++) 
	it <<= 1;
      tnm=it;
      del=(b-a)/tnm; /*This is the spacing of the points to be added.*/ 
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) 
	sum += FUNC(x); 
      s=0.5*(s+(b-a)*sum/tnm); /* This replaces s by its refined value. */
      return s;
    } 
}

double qromb(double (*func)(double), double a, double b)
/*Returns the integral of the function func from a to b. Integration is performed by Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.*/
{
  //void nrerror(char error_text[]);
  double ss,dss;
  double s[JMAXP],h[JMAXP+1]; 
  int j;
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd(func,a,b,j); 
    if (j >= KROMB) {
      /* These store the successive trapezoidal approxi- mations and their relative stepsizes.*/
      polint(&h[j-KROMB],&s[j-KROMB],KROMB,0.0,&ss,&dss);
      if (fabs(dss) <= ROMBTOL*fabs(ss)) 
	return ss; 
    }
    h[j+1]=0.25*h[j];
    /*This is a key step: The factor is 0.25 even though the stepsize is decreased by only 0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1), not just a polynomial in h.*/
  }
  printf("Too many steps in routine qromb\n"); 
  exit(-1);
  return 0.0; /*Never get here.*/
}

char fn[1024];


#define MC_BENT_DBLCYL

#ifdef MC_BENT_DBLCYL
/* apply a random rotation around the supplied axis because 
   bent cylinders do not have azimuthal symmetry */
double thetaGlobalBondangle;
void add_rotation_around_axis(double ox, double oy, double oz, double Rin[3][3], double Rout[3][3], double theta)
{
  double thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random rotation angle between 0 and 2*pi*/
 
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
void print_matrix(double M[3][3], int n)
{
  int k1, k2;
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

void versor_to_R(double ox, double oy, double oz, double gamma, double R[3][3])
{
  int k;
  double u[3], sp, norm;
#ifdef MC_BENT_DBLCYL
  double Rout[3][3];
  int k1, k2;
#endif
  /* first row vector */
  R[2][0] = ox;
  R[2][1] = oy;
  R[2][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[2][0] && u[1]==R[2][1] && u[2]==R[2][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[2][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[2][k];
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
  vectProdVec(R[1], R[2], u);
 
  for (k=0; k < 3 ; k++)
    R[0][k] = u[k];
#ifdef MC_BENT_DBLCYL
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
#ifdef DEBUG
  printf("==============\n");
  print_matrix(R, 3);
  printf("==============\n");
#endif
}
#if defined(EULER_ROT)
void place_HardCyl(double x, double y, double z, double cosphi12, double sinphi12, double costheta12, 
		double sintheta12, double cosgamma12, double singamma12,
		int which)
{
  double rO[3];
  double R[3][3];
  int k1, k2;
  rO[0] = x;
  rO[1] = y;
  rO[2] = z;
  R[0][0] = cosgamma12*costheta12*cosphi12 - singamma12*sinphi12;
  R[0][1] = cosphi12*singamma12 + cosgamma12*costheta12*sinphi12;
  R[0][2] = -cosgamma12*sintheta12;
  R[1][0] = -costheta12*cosphi12*singamma12-cosgamma12*sinphi12;
  R[1][1] = cosgamma12*cosphi12-costheta12*singamma12*sinphi12;
  R[1][2] = singamma12*sintheta12;
  R[2][0] = cosphi12*sintheta12;
  R[2][1] = sintheta12*sinphi12;
  R[2][2] = costheta12; 
  /* ============ */
  for (k1=0; k1 < 3; k1++)
    {
      HardCyl[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	HardCyl[which].R[k1][k2] = R[k1][k2];
    }
} 
#else
void place_HardCyl(double x, double y, double z, double ux, double uy, double uz, double gamma, int which)
{
  double xp[3], rO[3], xl[3];
  double R[3][3];
  int k1, k2;
#ifdef DEBUG
  FILE *fd;
  char fn[128];
#endif
  FILE *f;
  int i; 
  rO[0] = x;
  rO[1] = y;
  rO[2] = z;
  /* build R here from the orientation (ux,uy,uz) */
  versor_to_R(ux, uy, uz, gamma, R);
  /* ============ */
  for (k1=0; k1 < 3; k1++)
    {
      HardCyl[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	HardCyl[which].R[k1][k2] = R[k1][k2];
    }
}
#endif
/* ============================ >>> ranf <<< =============================== */
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
double fonsfact, dfonsfact; 

double fons(double costheta, double alpha)
{
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return cosh(alpha*costheta);
}
/* return an angle theta sampled from an Onsager angular distribution */
double theta_onsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta;
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
void angles_to_R(double *omx, double *omy, double* omz, double alpha)
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
double dfons(double costheta, double alpha)
{
  /* ho aggiunto un sin(theta) come giustamente fatto notare da Thuy, infatti la distribuzione 
     di Onsager si riduce a 1/(4*pi) e se non c'è il sin(theta) non è uniforma sull'angolo solido */
  return sinh(alpha*costheta);
}

/* return an angle theta sampled from an Onsager angular distribution */
double theta_donsager(double alpha, int domain)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  static int first = 1;
  static double f0;
  double pi, y, f, theta;
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
double calcDistNegHC(void)
{
#ifdef HC_FAST
  return calcDistNegHCopt();
#else
  return calcDistNegHCdiff();
#endif
}
double integrandv1(double rcmx, double rcmy, double rcmz, 
		    double phi12, int nphi12, double theta12, int ntheta12, double gamma12, int ngamma12,
		    double alpha)
{
  double sintheta12, costheta12, sinphi12, cosphi12, cosgamma12, singamma12;

  costheta12 = cos(theta12);
  sintheta12 = sin(theta12);
  cosphi12 = cos(phi12);
  sinphi12 = sin(phi12);
  cosgamma12 = cos(gamma12);
  singamma12 = sin(gamma12);
  //versor_to_R(u1x, u1y, u1z, gamma1, DNADall[0].R);
  //versor_to_R(u2x, u2y, u2z, gamma2, DNADall[1].R);
#ifdef EULER_ROT
  place_HardCyl(rcmx, rcmy, rcmz, cosphi12, sinphi12, costheta12, sintheta12, cosgamma12, singamma12, 1);
#else
  u2x = sintheta12*cosphi12;
  u2y = sintheta12*sinphi12;
  u2z = costheta12;  
  place_HardCyl(rcmx, rcmy, rcmz, u2x, u2y, u2z, gamma12, 1);
#endif
  if (calcDistBox(HardCyl[0].rcm,HardCyl[0].boxsax,HardCyl[0].R,
		  HardCyl[1].rcm,HardCyl[1].boxsax,HardCyl[1].R) < 0.0)
    {
      if (calcDistNegHC() < 0.0)
	{
	  switch (type)
	    {
	    case 0:
	      return XI1[nphi12][ntheta12];
	      break;
	    case 1:
	      return rcmx*XI1[nphi12][ntheta12]+
		rcmy*XI2[nphi12][ntheta12]+
		rcmz*XI3[nphi12][ntheta12];
	      break;
	    case 2: /* K22 */
	    case 3: /* K11 */
	    case 4: /* K33 */
	      return -(Sqr(rcmx)*XI1[nphi12][ntheta12]+
		       Sqr(rcmy)*XI2[nphi12][ntheta12]+
		       Sqr(rcmz)*XI3[nphi12][ntheta12]+rcmx*rcmy*XI4[nphi12][ntheta12]+
		       rcmx*rcmz*XI5[nphi12][ntheta12]+rcmy*rcmz*XI6[nphi12][ntheta12]);
	      break;
	    }
	}
    }
  return 0.0;
}

double phi12sav, theta12sav, gamma12sav;
#ifdef QFGAUSS
double rcmxsav, nrcmxsav, rcmysav, nrcmysav, rcmzsav, nrcmzsav;
int nrcmx, nrcmy, nrcmz;
double *xrcmx, *wrcmx, *xrcmy, *wrcmy, *xrcmz, *wrcmz;
#endif
int nphi12sav, ntheta12sav, ngamma12sav;
double ftheta12(double theta12, int ntheta12);
double fphi12(double phi12, int nphi12);
double fgamma12(double gamma12, int ngamma12);
#ifdef QFGAUSS
double frcmx(double rcmx, int nrcmx);
double frcmy(double rcmy, int nrcmy);
double frcmz(double rcmz, int nrcmz);
double (*nrfunc)(double,int,double,int,double,int,double,int,double,int);
double quad3d(double (*func)(double,int,double,int,double,int,double,int,double,int), 
	      double phi12_1, double phi12_2)
{
  nrfunc=func;
#ifdef GAUSS
  return qgaus(fphi12,phi12_1,phi12_2,xphi,wphi,nphi);
#else
  return qromb(fphi12,phi12_1,phi12_2);
#endif
}
double fphi12(double phi12, int nphi12) 
{
  phi12sav=phi12;
  nphi12sav=nphi12;
#ifdef GAUSS
  return qgaus(ftheta12,0.0,M_PI, xtheta, wtheta, ntheta); 
#else
  return qromb(ftheta12,0.0,M_PI); 
#endif
}
double ftheta12(double theta12, int ntheta12) 
{
  theta12sav=theta12;
  ntheta12sav=ntheta12;
#ifdef GAUSS
  /* notare che le ascisse e ordinate di phi vanno bene anche per theta poiché 
     gamma varia tra 0 e 2*pi come phi */
  return qgaus(frcmx,-Lx/2.,Lx/2., xrcmx, wrcmx, nrcmx); 
#else
  return qromb(frcmx,-Lx/2.,Lx/2.); 
#endif
}
double frcmx(double rcmx, int nrcmx) 
{
  rcmxsav=rcmx;
  nrcmxsav=nrcmx;
#ifdef GAUSS
  /* notare che le ascisse e ordinate di phi vanno bene anche per theta poiché 
     gamma varia tra 0 e 2*pi come phi */
  return qgaus(frcmy,-Ly/2.,Ly/2., xrcmy, wrcmy, nrcmy); 
#else
  return qromb(frcmy,-Ly/2.,Ly/2.); 
#endif
}
double frcmy(double rcmy, int nrcmy) 
{
  rcmysav=rcmy;
  nrcmysav=nrcmy;
#ifdef GAUSS
  /* notare che le ascisse e ordinate di phi vanno bene anche per theta poiché 
     gamma varia tra 0 e 2*pi come phi */
  return qgaus(frcmz,-Lz/2.,Lz/2., xrcmz, wrcmz, nrcmz); 
#else
  return qromb(frcmz,-Lz/2.,Lz/2.); 
#endif
}
double frcmz(double rcmz, int nrcmz) 
{
  return (*nrfunc)(phi12sav,nphi12sav,theta12sav,ntheta12sav,rcmxsav,nrcmxsav,rcmysav,nrcmysav,rcmz,nrcmz);
}
double rcmxsav, rcmysav, rcmzsav, alphasav;
double intfunc(double phi12, int nphi12, double theta12, int ntheta12, 
	       double rcmx, int nrcmx, double rcmy, int nrcmy, double rcmz, int nrcmz)
{
  return integrandv1(rcmx, rcmy, rcmz, phi12, nphi12, theta12, ntheta12, gamma12sav, 0, alphasav);
}
#elif defined(MCGAMMA)
double (*nrfunc)(double,int,double,int);
double quad3d(double (*func)(double,int,double,int), 
	      double phi12_1, double phi12_2)
{
  nrfunc=func;
#ifdef GAUSS
  return qgaus(fphi12,phi12_1,phi12_2,xphi,wphi,nphi);
#else
  return qromb(fphi12,phi12_1,phi12_2);
#endif
}
double fphi12(double phi12, int nphi12) 
{
  phi12sav=phi12;
  nphi12sav=nphi12;
#ifdef GAUSS
  return qgaus(ftheta12,0.0,M_PI, xtheta, wtheta, ntheta); 
#else
  return qromb(ftheta12,0.0,M_PI); 
#endif
}
double ftheta12(double theta12, int ntheta12) 
{
  return (*nrfunc)(phi12sav,nphi12sav,theta12, ntheta12);
}
double rcmxsav, rcmysav, rcmzsav, alphasav;
double intfunc(double phi12, int nphi12, double theta12, int ntheta12)
{
  return integrandv1(rcmxsav, rcmysav, rcmzsav, phi12, nphi12, theta12, ntheta12, gamma12sav, 0, alphasav);
}
#else
double (*nrfunc)(double,int,double,int,double,int);
double quad3d(double (*func)(double,int,double,int,double,int), 
	      double phi12_1, double phi12_2)
{
  nrfunc=func;
#ifdef GAUSS
  return qgaus(fphi12,phi12_1,phi12_2,xphi,wphi,nphi);
#else
  return qromb(fphi12,phi12_1,phi12_2);
#endif
}
double fphi12(double phi12, int nphi12) 
{
  phi12sav=phi12;
  nphi12sav=nphi12;
#ifdef GAUSS
  return qgaus(ftheta12,0.0,M_PI, xtheta, wtheta, ntheta); 
#else
  return qromb(ftheta12,0.0,M_PI); 
#endif
}
double ftheta12(double theta12, int ntheta12) 
{
  theta12sav=theta12;
  ntheta12sav=ntheta12;
#ifdef GAUSS
  /* notare che le ascisse e ordinate di phi vanno bene anche per theta poiché 
     gamma varia tra 0 e 2*pi come phi */
  return qgaus(fgamma12,0.0,M_PI, xgamma, wgamma, ngamma); 
#else
  return qromb(fgamma12,0.0,M_PI); 
#endif
}
double fgamma12(double gamma12, int ngamma12) 
{
  return (*nrfunc)(phi12sav,nphi12sav,theta12sav, ntheta12sav, gamma12, ngamma12);
}
double rcmxsav, rcmysav, rcmzsav, alphasav;
double intfunc(double phi12, int nphi12, double theta12, int ntheta12, double gamma12, int ngamma12)
{
  return integrandv1(rcmxsav, rcmysav, rcmzsav, phi12, nphi12, theta12, ntheta12, gamma12, ngamma12, alphasav);
}
#endif
#ifdef SOBOL_LL
static int iminarg1,iminarg2;
/*#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))*/
#define MAXBIT 60 
#define MAXDIM 6
void sobseq(int *n, double x[])
/*When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[1..n] the next values from n of these sequences. (n must not be changed between initializations.)*/
{
  int j,k,l;
  unsigned long long i,im,ipp;
  static double fac;
  static unsigned long long in,ix[MAXDIM+1],*iu[MAXBIT+1];
  static unsigned long long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
  static unsigned long long ip[MAXDIM+1]={0,0,1,1,2,1,4}; 
  static unsigned long long iv[MAXDIM*MAXBIT+1]={0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
  if (*n < 0) 
    { 
      /*Initialize, don’t return a vector. */
      for (k=1;k<=MAXDIM;k++) ix[k]=0;
      in=0;
      if (iv[1] != 1) return;
      fac=1.0/(1LL << MAXBIT);
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
#else
static int iminarg1,iminarg2;
/*#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))*/
#define MAXBIT 30 
#define MAXDIM 6
void sobseq(int *n, double x[])
/*When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[1..n] the next values from n of these sequences. (n must not be changed between initializations.)*/
{
  int j,k,l;
  unsigned long i,im,ipp;
  static double fac;
  static unsigned long  in,ix[MAXDIM+1],*iu[MAXBIT+1];
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
#endif
#ifdef SOBOLBF
void i8_sobol ( int dim_num, long long int *seed, double quasi[ ] );
#endif
#ifdef USE_MISER
double miser_func(double x[])
{
  rcmxsav = x[1];
  rcmysav = x[2];
  rcmzsav = x[3];
  gamma12sav = x[4];
  return quad3d(intfunc, 0.0, 2.0*M_PI);
}
#endif
double overlap_hc_c(double rA[3], double nA[3], double LA, double DA, 
                    double rB[3], double nB[3], double LB, double DB)
{
  HardCyl[0].sax[0] = DA*0.5;
  HardCyl[0].sax[1] = DA*0.5;
  HardCyl[0].sax[2] = LA*0.5;
  HardCyl[0].rcm[0] = rA[0];
  HardCyl[0].rcm[1] = rA[1];
  HardCyl[0].rcm[2] = rA[2];
  HardCyl[0].n[0] = nA[0];
  HardCyl[0].n[1] = nA[1];
  HardCyl[0].n[2] = nA[2];
  
  HardCyl[1].sax[0] = DB*0.5;
  HardCyl[1].sax[1] = DB*0.5;
  HardCyl[1].sax[2] = LB*0.5;
  HardCyl[1].rcm[0] = rB[0];
  HardCyl[1].rcm[1] = rB[1];
  HardCyl[1].rcm[2] = rB[2];
  HardCyl[1].n[0] = nB[0];
  HardCyl[1].n[1] = nB[1];
  HardCyl[1].n[2] = nB[2];
  return calcDistNegHC();
}
#if 0
int main(int argc, char**argv)
{
#ifdef USE_MISER
  double regn[4];
  double mis_ave, mis_var;
#endif
#ifdef QUASIMC
#ifdef SOBOLBF
  long long bfseed=0;
#endif
#ifdef USEGSL
  gsl_qrng *qsob;
#endif
#endif
#ifdef MPI
  MPI_Status status;
#endif
  char fn[256];
  int aa, bb;
  double ccc, totfact;
  FILE *fout, *fxi1, *fxi2, *fxi3, *fxi4, *fxi5, *fxi6;
  int i, j, cont=0, nfrarg;
  char fnout[256];
  double rcmx, rcmy, rcmz;
  double vexcl=0.0, factor;
  /* syntax:  CG-DNA-k2K22 <pdb file> <DNAD length> <tot_trials> <alpha> <type:0=v0, 1=v1, 2=v2> <outits> */
#if defined(MPI) 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
  //sprintf(TXT, "rank:%d\n", my_rank);
#endif

  if (argc < 7)
    {
#ifdef GAUSS
#ifdef QFGAUSS
      printf("syntax:  ELCONST_HC <diameter> <length>  <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [nphi] [ntheta] [rcmx] [rcmy] [rcmz] \n");
#else
#ifdef MCGAMMA
      printf("syntax:  ELCONST_HC <diameter> <length>  <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [nphi] [ntheta]\n");
#else
      printf("syntax:  ELCONST_HC <diameter> <length> <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [nphi] [ntheta] [ngamma]\n");
#endif
#endif
#else
      printf("syntax:  ELCONST_HC <diameter> <length> <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [romb-tol]\n");
#endif
      exit(1);
    }
  diamHC=atof(argv[1]);
  lengthHC=atof(argv[2]);
  alpha = atof(argv[3]);
  tot_trials=atoll(argv[4]);
  type = atoi(argv[5]);
  fileoutits = atoll(argv[6]);
  
  if (argc <= 7)
    outits=100*fileoutits;
  else
    outits = atoll(argv[7]);
#ifdef GAUSS 
  if (argc <= 8)
    nphi = 10;
  else
    nphi = atoi(argv[8]);
  if (argc <= 9)
    ntheta = 10;
  else
    {
      ntheta = atoi(argv[9]);
    }
  printf("qui argc=%d\n", argc);
#if !defined(MCGAMMA)  && !defined(QFGAUSS)
  if (argc <= 10)
    ngamma = 10;
  else
    ngamma = atoi(argv[10]);
#endif

#ifdef QFGAUSS
  if (argc <= 9)
    nrcmx = 10;
  else
    nrcmx = atoi(argv[9]);

 if (argc <= 10)
    nrcmy = 10;
  else
    nrcmy = atoi(argv[10]);

  if (argc <= 11)
    nrcmz = 10;
  else
    nrcmz = atoi(argv[11]);
#endif
#else
  if (argc <= 9)
    ROMBTOL = 1.0E-2;
  else
    ROMBTOL = atof(argv[9]);
#endif
  cont=0;
#ifdef GAUSS
  nfrarg = 11;
#else
  nfrarg = 10;
#endif
  
  vexcl = 0.0;
  ttini = 0;

  init_distbox();
  L=Lx=Ly=Lz=1.05*2.0*sqrt(Sqr(HardCyl[0].boxsax[0])+Sqr(HardCyl[0].boxsax[1])+Sqr(HardCyl[0].boxsax[2]))*2.0;
#if 0 
  Lx=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*2.0+2.0*DNADall[0].sax[0];
  Ly=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*2.0+2.0*DNADall[0].sax[1];
  Lz=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*2.0+2.0*DNADall[0].sax[2];
#endif
  printf("L=%f alpha=%f I am going to calculate v%d and I will do %lld trials\n", L, alpha, type, tot_trials);
  printf("box semiaxes=%f %f %f alpha=%f ntheta=%d nphi=%d\n", HardCyl[0].boxsax[0], HardCyl[0].boxsax[1], 
	 HardCyl[0].boxsax[2], alpha, ntheta, nphi);
#ifdef MPI
  srand48(((int)time(NULL))+my_rank);
#else
  srand48((int)time(NULL));
#endif
  sprintf(fnout, "v%d.dat", type);
  factor=0.0;
#if 0
  dfons_sinth_max=estimate_maximum_dfons(alpha);
  fons_sinth_max=dfons_sinth_max/alpha;
  printf("Estimated maximum of dfons is %f\n", dfons_sinth_max);
#endif
#ifdef QFGAUSS
  printf("Quasi Full Gauss quadrature with nrcmx=%d nrcmy=%d nrcmz=%d points\n", nrcmx, nrcmy, nrcmz);
#endif
#ifdef GAUSS
#ifdef MCGAMMA
  printf("Gauss quadrature with nphi=%d ntheta=%d points\n", nphi, ntheta);
#else
  printf("Gauss quadrature with nphi=%d ntheta=%d ngamma=%d points\n", nphi, ntheta, ngamma);
#endif
#else
  printf("Romberg method with %.15G tolerance\n", ROMBTOL);
#endif
#ifdef QUASIMC
  printf("I will generate a Quasi Monte Carlo sequence\n");
#endif
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
  fout = fopen(fnout, "w+");
  fclose(fout);
  //Lx=Ly=Lz=L;
  printf("Lx=%f Ly=%f Lz=%f\n", Lx, Ly, Lz);
  printf("type=%d\n", type);
  alphasav=alpha;
#if 0
  rcmysav=0.1;
  rcmxsav=rcmzsav=gamma1sav=gamma2sav=0.0;
  printf("val=%.15G\n", intfunc(0.2, 0.1, 0.2, 0.2));
  exit(-1);
#endif
  //nrfunc = intfunc;
#ifdef GAUSS
  xtheta = malloc(sizeof(double)*(ntheta+1));
  xphi = malloc(sizeof(double)*(nphi+1));
  xgamma = malloc(sizeof(double)*(ngamma+1));
  wtheta = malloc(sizeof(double)*(ntheta+1));
  wphi = malloc(sizeof(double)*(nphi+1));
  wgamma = malloc(sizeof(double)*(ngamma+1));
#ifdef QFGAUSS
  xrcmx = malloc(sizeof(double)*(nrcmx+1));
  wrcmx = malloc(sizeof(double)*(nrcmx+1));
  xrcmy = malloc(sizeof(double)*(nrcmy+1));
  wrcmy = malloc(sizeof(double)*(nrcmy+1));
  xrcmz = malloc(sizeof(double)*(nrcmz+1));
  wrcmz = malloc(sizeof(double)*(nrcmz+1));
#endif
  gauleg(0.0, M_PI, xtheta, wtheta, ntheta);
#if 0
  printf("x=%.15G %.15G %.15G %.15G %.15G\n w=%.15G %.15G %.15G %.15G %.15G\n",
       (M_PI/2.0-xtheta[1])/(M_PI/2.), xtheta[2]/M_PI, xtheta[3]/M_PI, xtheta[4]/M_PI, ((M_PI/2.)-xtheta[5])/(M_PI/2.0),
       wtheta[1]/(M_PI/2.), wtheta[2]/(M_PI/2.), wtheta[3]/(M_PI/2.), wtheta[4]/(M_PI/2.), wtheta[5]/(M_PI/2.));
  exit(-1);
#endif
  gauleg(0.0, 2.0*M_PI, xphi, wphi, nphi);
  gauleg(0.0, 2.0*M_PI, xgamma, wgamma, ngamma);
#endif
#ifdef QFGAUSS
  gauleg(-Lx/2., Lx/2., xrcmx, wrcmx, nrcmx);
  gauleg(-Ly/2., Ly/2., xrcmy, wrcmy, nrcmy);
  gauleg(-Lz/2., Lz/2., xrcmz, wrcmz, nrcmz);
  nsv = -1;
  sobseq(&nsv, sv);
  nsv = 1;
#endif
#ifdef QUASIMC
#ifdef USEGSL
#ifdef MCGAMMA
  nsv = 4; 
#else
  nsv = 3; 
#endif
  qsob = gsl_qrng_alloc (gsl_qrng_sobol, nsv);
#elif defined(SOBOLBF)
  // DO NOTHING
#else
  /* initialization */
  nsv = -1;  
  sobseq(&nsv, sv);
#ifdef MCGAMMA
  nsv = 4;
#else
  nsv = 3;
#endif
#endif
#endif
  
  XI1=malloc(sizeof(double)*(nphi+1));
  XI2=malloc(sizeof(double)*(nphi+1));
  XI3=malloc(sizeof(double)*(nphi+1));
  XI4=malloc(sizeof(double)*(nphi+1));
  XI5=malloc(sizeof(double)*(nphi+1));
  XI6=malloc(sizeof(double)*(nphi+1));


  for (i=1; i <= ntheta; i++)
    {
      XI1[i] = malloc(sizeof(double)*(ntheta+1));
      XI2[i] = malloc(sizeof(double)*(ntheta+1));
      XI3[i] = malloc(sizeof(double)*(ntheta+1));
      XI4[i] = malloc(sizeof(double)*(ntheta+1));
      XI5[i] = malloc(sizeof(double)*(ntheta+1));
      XI6[i] = malloc(sizeof(double)*(ntheta+1));
	
    }
  /* read XI1, X2 and X3 */
  
  sprintf(fn, "XI1_v%d.dat", type);
  if ((fxi1=fopen(fn, "r"))==NULL)
    {
      printf("You have to supply %s file\n", fn);
      exit(-1);
    }
  if (type >= 1)
    { 
      sprintf(fn, "XI2_v%d.dat", type);
      if ((fxi2=fopen(fn, "r"))==NULL)
	{
	  printf("You have to supply %s file\n", fn);
	  exit(-1);
	}
      sprintf(fn, "XI3_v%d.dat", type);
      if ((fxi3=fopen(fn, "r"))==NULL)
	{
	  printf("You have to supply %s file\n", fn);
	  exit(-1);
	}
    }
  if (type >= 2)
    {
      sprintf(fn, "XI4_v%d.dat", type);
      if ((fxi4=fopen(fn, "r"))==NULL)
	{
	  printf("You have to supply %s file\n", fn);
	  exit(-1);
	} 
      sprintf(fn, "XI5_v%d.dat", type);
      if ((fxi5=fopen(fn, "r"))==NULL)
	{
	  printf("You have to supply %s file\n", fn);
	  exit(-1);
	}
      sprintf(fn, "XI6_v%d.dat", type);
      if ((fxi6=fopen(fn, "r"))==NULL)
	{
	  printf("You have to supply %s file\n", fn);
	  exit(-1);
	}

    }
  fscanf(fxi1,"%lf %d %d\n", &ccc, &aa, &bb);
  if (aa!=nphi || bb!=ntheta|| ccc!= alpha)
    {
      printf("Wrong numbers of abscissas or wrong alpha!\n");
      printf("nphi=%d ntheta=%d aa=%d bb=%d alpha=%f/%f", nphi, ntheta, aa, bb, alpha, ccc);
      exit(-1);
    };
  if (type >= 1)
    {
      fscanf(fxi2,"%lf %d %d\n", &ccc, &aa, &bb);
      if (aa!=nphi || bb!=ntheta|| ccc!= alpha)
	{
	  printf("Wrong numbers of abscissas or wrong alpha!\n");
	  printf("nphi=%d ntheta=%d aa=%d bb=%d alpha=%f/%f", nphi, ntheta, aa, bb, alpha, ccc);
	  exit(-1);
	};
      fscanf(fxi3,"%lf %d %d\n", &ccc, &aa, &bb);
      if (aa!=nphi || bb!=ntheta || ccc!= alpha)
	{
	  printf("Wrong numbers of abscissas or wrong alpha!\n");
	  printf("nphi=%d ntheta=%d aa=%d bb=%d alpha=%f/%f", nphi, ntheta, aa, bb, alpha, ccc);
	  exit(-1);
	};
    }
  if (type >= 2)
    {
      fscanf(fxi4,"%lf %d %d\n", &ccc, &aa, &bb);
      if (aa!=nphi || bb!=ntheta|| ccc!= alpha)
	{
	  printf("Wrong numbers of abscissas or wrong alpha!\n");
	  printf("nphi=%d ntheta=%d aa=%d bb=%d alpha=%f/%f", nphi, ntheta, aa, bb, alpha, ccc);
	  exit(-1);
	};
      fscanf(fxi5,"%lf %d %d\n", &ccc, &aa, &bb);
      if (aa!=nphi || bb!=ntheta|| ccc!= alpha)
	{
	  printf("Wrong numbers of abscissas or wrong alpha!\n");
	  printf("nphi=%d ntheta=%d aa=%d bb=%d alpha=%f/%f", nphi, ntheta, aa, bb, alpha, ccc);
	  exit(-1);
	};
      fscanf(fxi6,"%lf %d %d\n", &ccc, &aa, &bb);
      if (aa!=nphi || bb!=ntheta || ccc!= alpha)
	{
	  printf("Wrong numbers of abscissas or wrong alpha!\n");
	  printf("nphi=%d ntheta=%d aa=%d bb=%d alpha=%f/%f", nphi, ntheta, aa, bb, alpha, ccc);
	  exit(-1);
	};

    }
  for (i=0; i < nphi; i++)
    for (j=0; j < ntheta; j++)
      {
	fscanf(fxi1, "%lf ", &(XI1[i+1][j+1]));
	if (type >= 1)
	  {
	    fscanf(fxi2, "%lf ", &(XI2[i+1][j+1]));
	    fscanf(fxi3, "%lf ", &(XI3[i+1][j+1]));
	  }
	if (type>=2)
	  {
	    fscanf(fxi4, "%lf ", &(XI4[i+1][j+1]));
	    fscanf(fxi5, "%lf ", &(XI5[i+1][j+1]));
	    fscanf(fxi6, "%lf ", &(XI6[i+1][j+1]));
	  }
      }
  
  fclose(fxi1);
  if (type >= 1)
    {
      fclose(fxi2);
      fclose(fxi3);
    }
  if (type >= 2)
    {
      fclose(fxi4);
      fclose(fxi5);
      fclose(fxi6);
    }
  //printf("XI1[7][8]:%.15G \n", XI1[7][8]);
  /* we use as the reference system the body reference system of first particle */
#ifdef EULER_ROT
  place_HardCyl(0.0, 0.0, 0.0, 1., 0., 1., 0., 1., 0., 0);      
#else
  place_HardCyl(0.0, 0.0, 0.0, 0., 0., 1., 0., 0);      
#endif
#if 0
  fonsfact= alpha/(4.0*M_PI*sinh(alpha));
  dfonsfact = alpha*alpha/(4.0*M_PI*sinh(alpha));
#endif
  totfact = 1.0/(2.0*M_PI);

#ifdef USE_MISER
  regn[1] = -Lx/2.;
  regn[2] = -Ly/2.;
  regn[3] = -Lz/2.;
  regn[4] = 0.;
  regn[5] = Lx/2.;
  regn[6] = Ly/2.;
  regn[7] = Lz/2.;
  regn[8] = 2.0*M_PI;  
  printf("I will use MISER algorithm\n");
  miser(miser_func, regn, 4, tot_trials, 0.1, &mis_ave, &mis_var);
  vexcl = mis_ave;
  printf("mis_ave=%.15G\n", mis_ave);
  fout = fopen(fnout, "a+");
  if (type==0)
    //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
    fprintf(fout,"%.15 %.15G\n", Lx*Ly*Lz*vexcl, mis_var);
  else if (type==1)
    fprintf(fout,"%.15G %.15G\n", Lx*Ly*Lz*vexcl, mis_var); /* divido per 10^4 per convertire in nm */
  else
    fprintf(fout,"%.15G %.15G\n", Lx*Ly*Lz*vexcl, mis_var); /* divido per 10^5 per convertire in nm */

  fclose(fout);
  return;
#endif
  for (tt=ttini+1; tt < tot_trials; tt++)
    {
      /* place second DNAD randomly */
#ifdef QFGAUSS
      sobseq(&nsv, sv);
      gamma12sav = 2.0*M_PI*sv[1];
#else
#ifdef QUASIMC
#ifdef USEGSL
      gsl_qrng_get (qsob, sv);
      //printf("sv=%f %f %f %f %f\n",sv[1], sv[2], sv[3], sv[4], sv[5]);
      rcmx = Lx*(sv[0]-0.5);
      rcmy = Ly*(sv[1]-0.5);
      rcmz = Lz*(sv[2]-0.5);
#ifdef MCGAMMA
      gamma12sav = 2.0*M_PI*sv[3];
#endif
#elif defined(SOBOLBF)
#ifdef MCGAMMA
      nsv = 4; 
#else
      nsv = 3; 
#endif
      i8_sobol(nsv, &bfseed, sv);
      //printf("sv=%f %f %f %f %f\n",sv[1], sv[2], sv[3], sv[4], sv[5]);
      rcmx = Lx*(sv[0]-0.5);
      rcmy = Ly*(sv[1]-0.5);
      rcmz = Lz*(sv[2]-0.5);
#ifdef MCGAMMA
      gamma12sav = 2.0*M_PI*sv[3];
#endif
#else
      /* quasi-MC */
      sobseq(&nsv, sv);
      //printf("sv=%f %f %f %f %f\n",sv[1], sv[2], sv[3], sv[4], sv[5]);
      rcmx = Lx*(sv[1]-0.5);
      rcmy = Ly*(sv[2]-0.5);
      rcmz = Lz*(sv[3]-0.5);
#ifdef MCGAMMA
      gamma12sav = 2.0*M_PI*sv[4];
#endif
#endif
#else
      rcmx = Lx*(drand48()-0.5);
      rcmy = Ly*(drand48()-0.5);
      rcmz = Lz*(drand48()-0.5);
      //gamma1 = 2.0*M_PI*drand48();
      //gamma2 = 2.0*M_PI*drand48();
#ifdef MCGAMMA
      gamma12sav = 2.0*M_PI*drand48();
#endif
#endif
      rcmxsav = rcmx;
      rcmysav = rcmy;
      rcmzsav = rcmz;
#endif

#if defined(MCGAMMA) || defined(QFGAUSS)
      vexcl += quad3d(intfunc, 0.0, 2.0*M_PI);
#else
      vexcl += quad3d(intfunc, 0.0, 2.0*M_PI)*totfact;
#endif
#if 1
      if (tt <= 1)
	{
#if 0
#ifdef POVRAY
	  sprintf(fname, "amyfibrils-%lld.pov", tt);
#else
	  sprintf(fname, "amyfibrils-%lld.mgl", tt);
#endif
	  //saveAmyloidPovray(fname, L);
#endif
	  //exit(-1);
	}
#endif

      if (tt > 0 && tt % fileoutits == 0)
	{
	  fout = fopen(fnout, "a+");
#ifdef QFGAUSS
  	  if (type==0)
	    //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	    fprintf(fout,"%lld %.15G\n", tt, vexcl/((double)tt));
	  else if (type==1)
	    fprintf(fout,"%lld %.15G\n", tt, (vexcl/((double)tt))); /* divido per 10^4 per convertire in nm */
	  else
	    fprintf(fout,"%lld %.15G\n", tt, (vexcl/((double)tt))); /* divido per 10^5 per convertire in nm */

#else
	  if (type==0)
	    //fprintf(fout,"%d %.15G %f %d\n", tt, L*L*L*vexcl/((double)tt)/1E3, vexcl, tt);
	    fprintf(fout,"%lld %.15G\n", tt, Lx*Ly*Lz*vexcl/((double)tt));
	  else if (type==1)
	    fprintf(fout,"%lld %.15G\n", tt, (Lx*Ly*Lz*vexcl/((double)tt))); /* divido per 10^4 per convertire in nm */
	  else
	    fprintf(fout,"%lld %.15G\n", tt, (Lx*Ly*Lz*vexcl/((double)tt))); /* divido per 10^5 per convertire in nm */
#endif
	  fclose(fout);
	}
      if (tt % outits==0)
	printf("trials: %lld/%lld\n", tt, tot_trials);
    } 
#if defined(USEGSL) && defined(QUASIMC)
  gsl_qrng_free(qsob);
#endif
}
#endif
