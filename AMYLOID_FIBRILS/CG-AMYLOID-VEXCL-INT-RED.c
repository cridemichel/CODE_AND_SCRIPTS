//#include "./G-DNA-k2K22.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
//#define USEGSL
#define GAUSS
#define EULER_ROT
#define MCGAMMA
#define PRINC_AXES
#define QUASIMC
#define SOBOL_LL /* use long long to have MAXBIT=60! */
#ifdef USEGSL
#include <gsl/gsl_qrng.h>
#endif
#if defined(MPI)
int MPIpid;
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
#endif 
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

double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
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
      amyloids[i].nD=2; /* la larghezza del foglietto (perpendicolare all'elica) deve essere di 4 nm */
      amyloids[i].Lbox=2.0;/* altezza lungo l'asse dell'elica del box in nm, l'unico vincolo è che sia abbastanza più piccola
			      del pitch di 70 nm delle fibre */
      amyloids[i].npitch=(int)(70.0/amyloids[i].Lbox); /* il pitch di ogni fibra deve essere di 70 nm circa */
      amyloids[i].nP=(int)(1980.0/amyloids[i].Lbox); /* persistence length in unità di Lbox */
      amyloids[i].Dproto=2.0; /* ogni protofilamente ha un diametro di circa 2 nm */
      amyloids[i].ribthick = 2.0; /* lo spessore del foglietto è di circa 2 nm */
      amyloids[i].boxsax[0] = (amyloids[i].nD+1)*amyloids[i].Dproto/2.0;
      amyloids[i].boxsax[1] = (amyloids[i].nD+1)*amyloids[i].Dproto/2.0;
      amyloids[i].boxsax[2] = amyloids[i].Lbox*(amyloids[i].nL+1)/2.0;
      dth=2.0*M_PI/amyloids[i].npitch;
      // length_eucl=OprogStatus.npitch*OprogStatus.pitch;
      for (jj=0; jj < amyloids[i].nL; jj++)
	{
	  //printf("x=%f\n", amyloids[i].boxes[jj].x[0]);
	  amyloids[i].boxes[jj].x[0]=0.0;
	  amyloids[i].boxes[jj].x[1]=0.0;
	  amyloids[i].boxes[jj].x[2]=((double)jj)*amyloids[i].Lbox;
	  //printf("rcm %f %f %f\n", amyloids[i].boxes[jj].x[0],amyloids[i].boxes[jj].x[1],amyloids[i].boxes[jj].x[2]);
	  amyloids[i].boxes[jj].sax[0] = amyloids[i].ribthick/2.0;
	  amyloids[i].boxes[jj].sax[1] = amyloids[i].nD*amyloids[i].Dproto/2.0;
	  amyloids[i].boxes[jj].sax[2] = amyloids[i].Lbox/2.0;
	  /* ...and now set orientation */
	  
	  amyloids[i].boxes[jj].R[0][0]=cos(jj*dth);
	  amyloids[i].boxes[jj].R[0][1]=-sin(jj*dth);
	  amyloids[i].boxes[jj].R[0][2]=0.0;
	  amyloids[i].boxes[jj].R[1][0]=sin(jj*dth);
	  amyloids[i].boxes[jj].R[1][1]=cos(jj*dth);
	  amyloids[i].boxes[jj].R[1][2]=0.0;
	  amyloids[i].boxes[jj].R[2][0]=0.0;
	  amyloids[i].boxes[jj].R[2][1]=0.0;
	  amyloids[i].boxes[jj].R[2][2]=1.0;}
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

double calcDistAmyloid(void)
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
    for (iB=0; iB < amyloids[1].nL; iB++)
      {

	if (calcDistBox(amyloids[0].boxesLab[iA].x,amyloids[0].boxes[iA].sax,amyloids[0].boxesLab[iA].R,
		        amyloids[1].boxesLab[iB].x,amyloids[1].boxes[iB].sax,amyloids[1].boxesLab[iB].R) < 0.0)
	  {
	    return -1;
	  }
      }
  return 1;
}

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

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


void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
//#define ALBERTA
char fn[1024];
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

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
  double angle, u[3], sp, norm, up[3], xx, yy;
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
void place_AMYLOID(double x, double y, double z, double cosphi12, double sinphi12, double costheta12, 
		double sintheta12, double cosgamma12, double singamma12,
		int which)
{
  FILE *fd;
  char fn[256];
  double xp[3], rO[3], xl[3];
  double R[3][3];
  int i, k1, k2;
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
      amyloids[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	amyloids[which].R[k1][k2] = R[k1][k2];
    }
} 
#else
void place_AMYLOID(double x, double y, double z, double ux, double uy, double uz, double gamma, int which)
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
      amyloids[which].rcm[k1] = rO[k1];
      for (k2=0; k2 < 3; k2++)
	amyloids[which].R[k1][k2] = R[k1][k2];
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

double integrandv1(double rcmx, double rcmy, double rcmz, 
		    double phi12, int nphi12, double theta12, int ntheta12, double gamma12, int ngamma12,
		    double alpha)
{
  int i, j;
  double sigsq, distsq, sigijsq, u1z, u2x, u2y, u2z;
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
  place_AMYLOID(rcmx, rcmy, rcmz, cosphi12, sinphi12, costheta12, sintheta12, cosgamma12, singamma12, 1);
#else
  u2x = sintheta12*cosphi12;
  u2y = sintheta12*sinphi12;
  u2z = costheta12;  
  place_AMYLOID(rcmx, rcmy, rcmz, u2x, u2y, u2z, gamma12, 1);
#endif
  if (calcDistBox(amyloids[0].rcm,amyloids[0].boxsax,amyloids[0].R,
		  amyloids[1].rcm,amyloids[1].boxsax,amyloids[1].R) < 0.0)
    {
      if (calcDistAmyloid() < 0.0)
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
  unsigned long long i,im,ipp;
  static double fac;
  static unsigned long long  in,ix[MAXDIM+1],*iu[MAXBIT+1];
  static unsigned long long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
  static unsigned long long ip[MAXDIM+1]={0,0,1,1,2,1,4}; 
  static unsigned long long iv[MAXDIM*MAXBIT+1]={0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
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
  int aa, bb, nL;
  double ccc, totfact;
  int cc;
  double gamma1, gamma2;
  FILE *fin, *fout, *f, *fread, *fxi1, *fxi2, *fxi3, *fxi4, *fxi5, *fxi6;
  int ncontrib, k, i, j, overlap, contrib, cont=0, nfrarg;
  char fnin[1024],fnout[256];
  double dummydbl, segno, u1x, u1y, u1z, u2x, u2y, u2z, rcmx, rcmy, rcmz;
  double sigijsq, distsq, vexcl=0.0, vexclel=0.0, factor, dth, th;
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
      printf("syntax:  CG-AMYLOID-k2K22 <AMYLOID length>  <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [nphi] [ntheta] [rcmx] [rcmy] [rcmz] \n");
#else
#ifdef MCGAMMA
      printf("syntax:  CG-AMYLOID-k2K22 <AMYLOID length>  <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [nphi] [ntheta]\n");
#else
      printf("syntax:  CG-AMYLOID-k2K22 <AMYLOID length>  <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [nphi] [ntheta] [ngamma]\n");
#endif
#endif
#else
      printf("syntax:  CG-DNA-k2K22 <AMYLOID length> <alpha> <tot_trials> <type:0=v0, 1=v1, 2=v2(K22) 3=v2(K11) 4=v2(K33)> <fileoutits> [outits] [romb-tol]\n");
#endif
      exit(1);
    }
  nL=len=atoi(argv[1]);
  alpha = atof(argv[2]);
  tot_trials=atoll(argv[3]);
  type = atoi(argv[4]);
  fileoutits = atoll(argv[5]);
  
  if (argc <= 6)
    outits=100*fileoutits;
  else
    outits = atoll(argv[6]);
#ifdef GAUSS 
  if (argc <= 7)
    nphi = 10;
  else
    nphi = atoi(argv[7]);
  if (argc <= 8)
    ntheta = 10;
  else
    {
      ntheta = atoi(argv[8]);
    }
  printf("qui argc=%d\n", argc);
#if !defined(MCGAMMA)  && !defined(QFGAUSS)
  if (argc <= 9)
    ngamma = 10;
  else
    ngamma = atoi(argv[9]);
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

  build_amyloid(nL);

  L=Lx=Ly=Lz=1.05*2.0*sqrt(Sqr(amyloids[0].boxsax[0])+Sqr(amyloids[0].boxsax[1])+Sqr(amyloids[0].boxsax[2]))*2.0;
#if 0 
  Lx=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*2.0+2.0*DNADall[0].sax[0];
  Ly=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*2.0+2.0*DNADall[0].sax[1];
  Lz=1.05*2.0*sqrt(Sqr(DNADall[0].sax[0])+Sqr(DNADall[0].sax[1])+Sqr(DNADall[0].sax[2]))*2.0+2.0*DNADall[0].sax[2];
#endif
  printf("nL=%d L=%f alpha=%f I am going to calculate v%d and I will do %lld trials\n", nL, L, alpha, type, tot_trials);
  printf("box semiaxes=%f %f %f alpha=%f ntheta=%d nphi=%d\n", amyloids[0].boxsax[0], amyloids[0].boxsax[1], 
	 amyloids[0].boxsax[2], alpha, ntheta, nphi);
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
  place_AMYLOID(0.0, 0.0, 0.0, 1., 0., 1., 0., 1., 0., 0);      
#else
  place_AMYLOID(0.0, 0.0, 0.0, 0., 0., 1., 0., 0);      
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
