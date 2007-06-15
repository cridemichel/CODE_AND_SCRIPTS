#include <stdio.h> /* standart input,output */
#include <ctype.h> /* for charachter recognition */
#include <stdlib.h> /* for conversion from char to dec */
#include <strings.h>
#include <time.h>
#include <math.h>
#include <float.h>
#undef DEBUG
#undef DEBUGP
#define NAMINO 60
double *rx, *ry, *rz, *vx, *vy, *vz, *wx, *wy, *wz, ***Ri, omega[3];
double eigenVect[3][3];
int Namino;
int brownian;
double temp, L, massAmino, T, Iamino[3], x[4][3], tetraEdge, dL, PBw, Dh, Dh2;
double sigCA, sigBC, sigNC, sigCC, sigAC, sigAN, sigNN, sigBA, sigBN, sigAA, sigBB,
       sigCN1, sigCN2, sigCCA1, sigCCA2, sigCACA1, sigCACA2, sigCAN1, sigCAN2,
       sigNCA1, sigNCA2, sigNC1, sigNC2, sigCAC2, sigCAC2, sigCAC1, sigCAC2, 
       A[3], B[3], C[3]; 
double  xT, rT, dT;
#ifdef PEPTIDE_PLATE
double massPlate, sigPepCA, sigPepC, sigPepN, sigCC_Ami, sigCACA_Ami, sigNN_Ami;
double Iplate[3], xpep[4][3];
double *rxP, *ryP, *rzP, *vxP, *vyP, *vzP, *wxP, *wyP, *wzP, ***RiP;
int Nplate;
#endif
#include <float.h>

void diag_sym(int n, double ** A, double *lambda, double ** B)
{
  unsigned int mult= 314159261;
  unsigned int add =907633385;
  double big=4294967296.0;
  unsigned int rn;
  unsigned int * rn1;
  int i,j,k;
  double *a,*b,*c, s,ss,sb,sss,lnew,lold,eps=1.0e-12;
  double fact;
  a=(double *)malloc(n*sizeof(double));
  b=(double *)malloc(n*sizeof(double));
  fact=1/big;
  rn1=(unsigned int *)A[0];
  rn=rn1[0];

  for(k=0;k<n;k++)
    {
      for(i=0;i<n;i++)
	{
	  rn=rn*mult+add;
	  a[i]=fact*rn;
	}
      lnew=DBL_MAX;
      do{
	if(k>0)
	  for(j=0;j<k;j++)
	    {
	      c=B[j];
	      s=c[0]*a[0];
	      for(i=1;i<n;i++)
		s+=c[i]*a[i];
	      for(i=0;i<n;i++)
		a[i]-=s*c[i];
	    }

	ss=0;
	for(i=0;i<n;i++)
	  ss+=a[i]*a[i];
	ss=1/sqrt(ss);


	for(i=0;i<n;i++)
	  a[i]*=ss;
	ss=0;
	sss=0;
	for(i=0;i<n;i++)
	  {
	    c=A[i];
	    s=a[0]*c[0];
	    for(j=1;j<n;j++)
	      s+=a[j]*c[j];
	    ss+=s*s;
	    sss+=s*a[i];
	    b[i]=s;      
	  }
	ss=sqrt(ss);
	if(sss<0)ss=-ss;
	lold=lnew;
	lnew=ss;
	ss=1/ss;
	sss=0;
	for(i=0;i<n;i++)
	  {
	    sb=ss*b[i];
	    s=a[i]-sb;
	    sss+=s*s;
	    a[i]=sb;

	  }
	sss=sqrt(sss);
      }while(fabs(lnew-lold)>eps);
      c=B[k];
      lambda[k]=lnew;
      for(i=0;i<n;i++)
	{
	  c[i]=a[i];    
	  /*    printf("%lf\n", c[i]);*/
	}
    }
  free(a);
  free(b);
  return;
}

/* solve a linear system of equations */
void wrap_dgesv(double a[3][3], double x[3], int *ok)
{
  double AT[9];
  int i, j, c1, c2, pivot[3];
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) AT[j+3*i]=a[j][i];		
    }						
  c1 = 3;
  c2 = 1;
  dgesv_(&c1, &c2, AT, &c1, pivot, x, &c1, ok);      
}
int SolveLineq (double a[3][3], double x[3]) 
{
  int indx[3], ok;
  double dd;
  wrap_dgesv(a, x, &ok);
  return 0;
}

/* find eigenvectors and eigenvalues */
void wrap_dsyev(double a[3][3], double b[3][3], double x[3], int *ok)
{
  char JOBZ, UPLO;
  double AT[9], work[10];
  int i, j, c1, c2, c3;
  JOBZ='V';
  UPLO='U';
  extern void dsyev_(char* , char*, int*, double* , int*, double*, double*, int*, int*);
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) AT[j+3*i]=a[j][i];		
    }						
  c1 = 3;
  c2 = 1;
  c3 = 8;
  dsyev_(&JOBZ, &UPLO, &c1, AT, &c1, x, work, &c3, ok);      
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) b[i][j]=AT[j+3*i];		
    }	
  if (*ok != 0)
    printf("not ok (%d)\n", *ok);
}

void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
void lab2body(double x[3], double xp[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      xp[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  xp[k1] += R[k1][k2]*(x[k2]-rO[k2]);
       	} 
    }
}

void body2labR(double xp[3], double x[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
    }
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

#define Sqr(x) ((x)*(x))
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
double ranf(void)
{
 return rand() / ( (double) RAND_MAX );
}
void evalfourth(double *r1, double *r2, double *r3, double *r4, double radius)
{
  double nr, r21[3], r31[3];
  int kk;
  for (kk = 0; kk < 3; kk++)
    {
      r21[kk] = r2[kk]-r1[kk];
      r31[kk] = r3[kk]-r1[kk];
    }
#if 0
  printf("%f %f %f @ 0.5 C[green]\n", r1[0], r1[1], r1[2]);
  printf("%f %f %f @ 0.5 C[red]\n", r2[0], r2[1], r2[2]);
  printf("%f %f %f @ 0.5 C[orange]\n", r3[0], r3[1], r3[2]);
#endif
  vectProdVec(r21, r31, r4);
  nr = calc_norm(r4);
  for (kk = 0; kk < 3; kk++)
    r4[kk] *= radius/nr;
#if 0
  printf("%f %f %f @ 0.5 C[blue]\n", r4[0], r4[1], r4[2]);
#endif
  //printf("norm_r4=%.15G\n", calc_norm(r4));


}
#if 0
double calc_energy(int nummol)
{
  int i, k1;
  double K,wt[3],m,I;
  K = 0;
  for (i=0; i < nummol; i++)
    {
	 /* calcola tensore d'inerzia e le matrici delle due quadriche */
	  m = massAmino;
	  I = Iamino;
	  K += m*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    K += Sqr(wt[k1])*I;
    }
  K *= 0.5;
  return K;
}
#endif

double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
void init_parameters(void)
{
  double bondDistAvg;
  L=230.0;
  xT = sqrt(3.0)/3.0;
  rT = sqrt(6.0)/12.0;
  dT=sqrt(3.0)/6.0;
  T = 1.0;
  tetraEdge = 3.25;
  dL = 1.5;//1.5*tetraEdge;
  Namino=60;
  PBw = 0.1; /* larghezza della buca peptidica */
  sigNN= 3.3800000;
  sigCC= 3.5200000;
  sigAN= 3.3600000;
  sigAC= 3.4500000;
  sigAA= 3.4100000;
  sigBN= 3.1700000;
  sigBC= 3.2600000;
  sigBA= 3.2400000;
  sigBB= 3.0700000;
  sigNC= 3.1900000;
  /* values given by Sergey 06/06/07 */
#ifdef PEPTIDE_PLATE
  sigCC_Ami = 0.1;
  sigCACA_Ami = 0.1;
  sigNN_Ami = 0.1;
#else
  /* diametri spot per legame peptidico */
  sigCN1 = 1.2985;
  sigCN2 = 1.3515;
  sigCAN1 = 2.35788;
  sigCAN2 = 2.45412;
  sigCACA1 = 3.70832;
  sigCACA2 = 3.85968;
  sigCCA1 = 2.38336;
  sigCCA2 = 2.48064;
  sigNCA1 = 2.35788;
  sigNCA2 = 2.45412;
  sigNC1 = 1.2985;
  sigNC2 = 1.3515;
  sigCAC1 = 2.38336;
  sigCAC2 = 2.48064;
#endif
#if 0
  sigCAN1 = dL - tetraEdge*0.5-PBw;
  sigCAN2 = dL - tetraEdge*0.5+PBw;
  sigNCA1 = sigCAN1;
  sigNCA2 = sigCAN2;
  sigCCA1 = dL - tetraEdge*0.5-PBw;
  sigCCA2 = dL - tetraEdge*0.5+PBw;
  sigCAC1 = sigCCA1;
  sigCAC2 = sigCCA2;
  bondDistAvg = sqrt(Sqr((xT+dT)*tetraEdge)+Sqr(dL-tetraEdge));
  sigCN1 = bondDistAvg - PBw;
  sigCN2 = bondDistAvg + PBw;
  sigNC1 = sigCN1;
  sigNC2 = sigCN2;
  bondDistAvg = sqrt(Sqr(dL)+Sqr((xT+dT)*tetraEdge));
  sigCACA1 = bondDistAvg - PBw;
  sigCACA2 = bondDistAvg + PBw;
#endif
  massAmino = 1.0;
#ifdef PEPTIDE_PLATE
  massPlate = 1.0;
  Iplate[0] = Iplate[1] = Iplate[2] = 1.0;
  sigPepCA = 0.1;
  sigPepC = 0.1;
  sigPepN = 0.1;
#endif
  Iamino[0] = Iamino[1] = Iamino[2] = 1.0;
  brownian = 0;
}
void angvel(int i, double *wx, double *wy, double* wz, double temp)
{
  int a;
  double inert;                 /* momentum of inertia of the molecule */
  double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, Mtot;
  double r21[3], symax[3], wsz, ww[3]; 
  //Mtot = m; /* total mass of molecule */

  //inert = I; /* momentum of inertia */

  Mtot = massAmino;
  inert = Iamino[0];
  mean = 2.0*temp / inert;
  xisq = 1.0;

  while (xisq >= 1.0)
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

  /* Choose the magnitude of the angular velocity
   * NOTE: consider that it is an exponential distribution 
   (i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/

  osq   = - mean * log(ranf());
  o     = sqrt(fabs(osq));
  ox    = o * ox;
  oy    = o * oy;
  oz    = o * oz;
  *wx = ox;
  *wy = oy;
  *wz = oz;
}
double polyAlaCoord[60][4][3];
/* ordine sergey N CA C CB */
double aminoXYZLab[4][3]={
      {1.211991, 0.805049, 0.000000}, /* CA */
      {1.218941, 1.723877, 1.227108},/* CB */
      {0.000000, 0.000000, 0.000000}, /* N  */
      {2.443014, -0.069412, 0.000000}};/* C  */
void calcCOM(double xin[4][3], double com[3])
{
  int a, kk;

  for (kk=0; kk < 3; kk++)
    com[kk] = 0.0;

  for (a=0; a < 4; a++)
    for (kk=0; kk < 3; kk++)
      com[kk] += xin[a][kk]; 
  
  for (kk=0; kk < 3; kk++)
    com[kk] /= 4;
}

/* N.B. mass of atoms CA, CB, N, C, now for simplicity are all assumed to be 1 */
double mass[4]={1.0,1.0,1.0,1.0};
void calcItensSergey(double pos[4][3], double m[4], double I[3][3])
{
  double s;
  int j, k, l, n_dim=3;
  s=0;
  for(k=0;k<n_dim;k++)
    {
      for(l=0;l<n_dim;l++)
	{
	  I[k][l]=0; 
	  for(j=0;j<4;j++)
	    I[k][l]-=pos[j][k]*pos[j][l];
	}
      s-=I[k][k];
    }
  for(k=0;k<n_dim;k++)
    {
      I[k][k]+=s;
    }
}
void calcItens(double pos[4][3], double m[4], double I[3][3])
{
  int i, j, k;
  double distSq;
  for (j=0; j < 3; j++)
    for (k=0; k < 3; k++)
      {
	I[j][k] = 0.0;
	for (i=0; i < 4; i++)
	  {
	    distSq = Sqr(pos[i][0])+Sqr(pos[i][1])+Sqr(pos[i][2]);
	    I[j][k] += m[i]*(((j==k)?distSq:0.0) - pos[i][j]*pos[i][k]);
	  }
      }
}
void sort_eigenvect(double R[3][3], double EV[3])
{
  /* sort from smallest eigenvalues to bigges */
  int i, j, k;
  double u[3], ev;
  for (i=0; i < 2; i++)
    for (j=i+1; j < 3; j++)
      {
	if (EV[i] > EV[j])
	  {
	    /* swap eigenvectors */
	    for (k=0; k < 3; k++)
	      u[k] = R[i][k];
	    for (k=0; k < 3; k++)
	      R[i][k] = R[j][k];
	    for (k=0; k < 3; k++)
	      R[j][k] = u[k];
	    /* swap eigenvalues */
	    ev = EV[i];
	    EV[i] = EV[j];
	    EV[j] = ev;
	  }
      }
}
void calcEigenValVect(double I[3][3], double R[3][3], double EV[3])
{
  int ok, kk;
  double u1[3], u2[3];
  wrap_dsyev(I, R, EV, &ok);
  sort_eigenvect(R, EV);

  /* we want a right-handed reference system */
  for (kk=0; kk < 3; kk++)
    u1[kk] = R[0][kk];
  for (kk=0; kk < 3; kk++)
    u2[kk] = R[1][kk];
  vectProdVec(u1, u2, R[2]);
}
double** AllocMatR(int size1, int size2);
const double minDbl = 1E-10;
double deterM(double M[3][3])
{
  return M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[2][0]*M[1][1]*M[0][2] 
    - M[2][1]*M[1][2]*M[0][0] - M[2][2]*M[1][0]*M[0][1]; 
}
void check_eigenval(int type, double xlab[4][3], double R[3][3], double EV[3], double com[3], double x[4][3])
{
  double xbody[4][3], Rini[3][3];
  int jj, kk, ii, l;
  /* computed eigenvectors are unitary but may have opposite 
   * versus with respect to the desired body fixed frame */
	
  for (jj=0; jj < 4; jj++)
    {
      lab2body(xlab[jj], xbody[jj], com, R);
    }
  //printf("prima detR=%.15G\n", deterM(R));
  for (l=0; l < 3; l++)
    for (jj=0; jj < 4; jj++)
      {
	//printf("xbody[%d][%d] %.15G x=%.15G\n", jj, l, xbody[jj][l], x[jj][l]);
	if (fabs(xbody[jj][l]-x[jj][l]) > minDbl && xbody[jj][l]*x[jj][l] < 0.0)
	  {
	    for (kk=0; kk < 3; kk++)
	      R[l][kk] = -R[l][kk];
	    break;
	  }
      }
  /* dirty fix for plates: they have z-coords=0 in the rigid body ref. system */
  if (deterM(R) < 0.0)
    {
      for (l=0; l < 3; l++)
	{
	  if (fabs(xbody[0][l]) < minDbl)
	    {
	      for (kk=0; kk < 3; kk++)
		R[l][kk] = -R[l][kk];
	    }
	}
    }
  //printf("dopo detR=%.15G\n", deterM(R));
}
void calcEigenValVectSergey(double I[3][3], double R[3][3], double EV[3])
{
  int ok, kk, ii;
  double **Il, **Rl;
  Il=AllocMatR(3,3);
  Rl=AllocMatR(3,3);
  for (ii=0; ii < 3; ii++)
    for (kk=0; kk < 3; kk++)
      Il[ii][kk] = I[ii][kk];
  diag_sym(3, Il, EV, Rl);
  for (ii=0; ii < 3; ii++)
    for (kk=0; kk < 3; kk++)
	R[ii][kk] = Rl[ii][kk];

  free(Rl);
  free(Il);
}
double aminoPos[60][4][3];
#ifdef PEPTIDE_PLATE
double platePos[60][4][3];
#endif
void print_distances(double xx[4][3])
{
  int jj, kk;
  for (jj = 0; jj < 4; jj++)
    for (kk = 0; kk < 4; kk++)
      printf("distance %d-%d: %.15G\n", jj, kk, sqrt(Sqr(xx[jj][0]-xx[kk][0])+Sqr(xx[jj][1]-xx[kk][1])+
						     Sqr(xx[jj][2]-xx[kk][2])));
}

double spotsPos[4][3]={
    {-0.119826755985092,0.689377462631642,-1.25388776798239},/* CA */
    {-0.115327445521127,0.725512097842898,1.2398463722442},/* CB*/ 
    {-0.125004064249858,-1.39306360336586,0.0253895049617639}, /* N */
    {0.360158265756074,-0.0218259571086818,-0.0113481092235923}};/* C */
void buildAminoSpotsSergey(double xout[4][3])
{
  int ii, jj;
  for (ii=0; ii < 4; ii++)
    for (jj=0; jj < 3; jj++)
      xout[ii][jj] = spotsPos[ii][jj];
}
#ifdef PEPTIDE_PLATE
void print_matrix(char *txt, double M[3][3]);
double massPlateSpots[4]={1.0,1.0,1.0,1.0};
void buildPeptidePlate(double xout[4][3])
{
  int a, kk, jj, jj2;
  double  com[3], R[3][3], xtmp[4][3], Itens[3][3], Iev[3], Rt[3][3];
  double xlab[4][3];
  int i;
#ifdef DEBUGP
  for (kk=0; kk < 4; kk++)
    printf("[BUILDPLATEPOS]platePos[0][%d]:%f %f %f\n", kk, platePos[0][kk][0],platePos[0][kk][1],platePos[0][kk][2]);
#endif  
  calcCOM(platePos[0],com); 
  for (a=0; a < 4; a++)
    for (kk=0; kk < 3; kk++)
      xtmp[a][kk] = platePos[0][a][kk] - com[kk];
  
  calcItens(xtmp, massPlateSpots, Itens);
  /* diagonalizza il tensore d'inerzia e ottiene gli autovettori */
  calcEigenValVect(Itens, R, Iev);
  //printf("Using %d Eigenvalues=%.15G %.15G %.15G\n", use, Iev[0], Iev[1], Iev[2]);
  if (deterM(R) < 0.0)
    {
      for (jj=0; jj < 3; jj++)
	{
	  R[2][jj] = -R[2][jj];
	}
    }
  for (jj=0; jj < 3; jj++)
    {
      for (kk=0; kk < 3; kk++)
	{
  	  eigenVect[jj][kk] = R[jj][kk];
	}
    }
  for (jj=0; jj < 4; jj++)
    lab2body(platePos[0][jj], xout[jj], com, R);

#ifdef DEBUGP
  for (jj=0; jj < 4; jj++)
   {
     printf("plate spots[%d]=(%.15G,%.15G,%.15G)\n", jj, xout[jj][0], xout[jj][1], xout[jj][2]);
   }
  print_distances(xout);
#endif
#if 1
#ifdef DEBUGP
     printf("[BUILPLATEPOS] i=0\n");
     i=0;
     for (jj=0; jj < 4; jj++)
       {
	 print_matrix("i=0 Orientation Matrix",R);
	 printf("i=0 COM=%f %f %f\n", com[0], com[1], com[2]);
	 body2lab(xout[jj], xlab[jj], com, R);
	 printf("xout[%d]=%f %f %f\n", jj, xout[jj][0], xout[jj][1], xout[jj][2]);
	 printf("xout=%.15G %.15G %.15G\n", xlab[jj][0], xlab[jj][1], xlab[jj][2]);
	 printf("sergey=%.15G %.15G %.15G\n", platePos[i][jj][0],platePos[i][jj][1],platePos[i][jj][2]);
	 lab2body(platePos[0][jj], xlab[jj], com, R);
	 printf("plate spots=%.15G %.15G %.15G\n", xlab[jj][0], xlab[jj][1], xlab[jj][2]);
       }
#endif
#ifdef DEBUGP
      if (i==0 || i==1)
       {
	 printf("BEG PLATE[BUILDPLATEPOS] ====> eigenvalues: %f %f %f\n", Iev[0], Iev[1], Iev[2]);
	 printf("COM=%.15G %.15G %.15G\n", com[0], com[1], com[2]);
	 printf("i=%d ------------ \n", i);
	 print_matrix("Itens",Itens);
	 printf("===================== END[BUILDPLATEPOS] ================ \n");
       }
#endif
 #endif
 }
#endif
void buildAminoSpots(int use, double xout[4][3])
{
  int a, kk, jj, jj2;
  double  com[3], R[3][3], xtmp[4][3], Itens[3][3], Iev[3], Rt[3][3];
      
  calcCOM(aminoPos[use],com); 
  /* calcola il tensore d'inerzia */
  /* i vettori riga di R sono i versori che individuano il 
     sistema di riferimento solidale con l'amminoacido cioè
     con il corpo rigido */
  for (a=0; a < 4; a++)
    for (kk=0; kk < 3; kk++)
      xtmp[a][kk] = aminoPos[use][a][kk] - com[kk];
#if 1
  calcItens(xtmp, mass, Itens);
#else
  calcItensSergey(xtmp, mass, Itens); 
#endif
  /* diagonalizza il tensore d'inerzia e ottiene gli autovettori */
  calcEigenValVect(Itens, R, Iev);
  //printf("Using %d Eigenvalues=%.15G %.15G %.15G\n", use, Iev[0], Iev[1], Iev[2]);
  for (jj=0; jj < 3; jj++)
    {
      for (kk=0; kk < 3; kk++)
	eigenVect[jj][kk] = R[jj][kk];
    }

#if 0
  for (jj=0; jj < 3; jj++)
    {
      for (kk=0; kk < 3; kk++)
	printf("%.15G ", R[jj][kk]);
      printf("\n");
    }
  for (jj=0; jj < 3; jj++)
    printf("EV[%d]=%.15G\n", jj, Iev[jj]);
  //exit(-1);
#endif
#if 0
  for (kk = 0; kk < 3; kk++)
    for (jj = 0; jj < 3; jj++)
      Rt[kk][jj] = R[jj][kk];
#endif
  for (jj=0; jj < 4; jj++)
    lab2body(aminoPos[use][jj], xout[jj], com, R);
#ifdef DEBUG
  for (jj=0; jj < 4; jj++)
   {
     printf("spots[%d]=(%.15G,%.15G,%.15G)\n", jj, xout[jj][0], xout[jj][1], xout[jj][2]);
   }
  print_distances(xout);
#endif
}
void buildTetrahedra(double xout[4][3], double alpha)
{
  int i;
  double x, r, d;
  const double Kl = sqrt(8.0/3.0), Kdh = 1.0/3.0, Ktr = sqrt(8.0)/6.0;
  double Oangle;
  double radius; 
  //Oangle = acos(0) * 2.0 * 145.8 / 180.0;
  Oangle = acos(-1.0/3.0);
  radius = alpha*sqrt(6.0)/4.0; /* typical dimension of aminoacid is 3.5 Angstrom */
  x = xT * alpha;
  d = dT * alpha;
  r = rT* alpha;
#if 1
  /* il raggio è quello dell'interazione Si-O */
  /* CA */
  xout[0][0] = x;
  xout[0][1] = 0;
  xout[0][2] = -r;
  /* CB */
  xout[1][0] = 0.0;
  xout[1][1] = 0.0;
  xout[1][2] = radius;
  //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
  /* N */
  xout[2][0] = -d;
  xout[2][1] = alpha*0.5;
  xout[2][2] = -r;
  //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
  /* C */
  xout[3][0] = -d;
  xout[3][1] = -alpha*0.5;
  xout[3][2] = -r;
    
 //printf("W xout %f %f %f\n", xout[3][0], xout[3][1], xout[3][2]);
  //evalfourth(xout[0], xout[1], xout[2], xout[3], radius);
  //printf("M xout %f %f %f\n", xout[3][0], xout[3][1], xout[3][2]);
  //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
#else
  /* il raggio è quello dell'interazione Si-O */
  /* CA */
  xout[0][0] = Ktr * 2.0 * radius; /* =  x */
  xout[0][1] = 0.0;
  xout[0][2] = -Kdh * radius;      /* = -r */
  /* N */
  xout[2][0] = -Ktr * radius;     /* = -d */
  xout[2][1] = Kl * radius / 2.0; /* = +alpha/2 */
  xout[2][2] = -Kdh * radius;     /* = -r */
  //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
  /* C */
  xout[3][0] = -Ktr * radius;      /* = -d  */
  xout[3][1] = -Kl * radius / 2.0; /* = -alpha/2 */
  xout[3][2] = -Kdh * radius;      /* = -r */
  //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
  evalfourth(xout[0], xout[2], xout[3], xout[1], radius);
  //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
#endif
#if 0
{double dr[3];
      int a, b, kk;
    for (a=0; a < 4; a++)
   for (b=0; b  < 4; b++)
     {
       for (kk=0; kk < 3; kk++)
	 dr[kk] = xout[a][kk] - xout[b][kk];
       printf("dist(%d,%d)=%.15G\n", a, b, calc_norm(dr));
     }
    }
#endif
}
void fitTetrahedra (double xin[4][3], double tetrahedraBF[4][3], double rcm[3], double R[3][3])
{
/* given the coordinates of 4 atoms and the coordinates of the 4 vertices of a regular tetraheadra 
 * this routine evaluate the best orientation and center of mass position of the tetrahedra in order
 * to minimize the least squared distances between atoms and tetrahedra vertices. */

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
void readSergeyPos(char *fn)
{
  FILE *f;
  int i, kk;
  if (!(f = fopen(fn,"r")))
    {
      printf("problem opening file %s\n", fn);
      exit(-1);	
    }
  i=0;
  while (!feof(f))
    {
      /* N.B. nella mail diceva che la sequenza era N, CA, C, CB! */
      /* ordine sergey N, C, CA, CB */
      fscanf(f,"%lf %lf %lf\n", &aminoPos[i][2][0], &aminoPos[i][2][1], &aminoPos[i][2][2]);  
      fscanf(f,"%lf %lf %lf\n", &aminoPos[i][3][0], &aminoPos[i][3][1], &aminoPos[i][3][2]);
      fscanf(f,"%lf %lf %lf\n", &aminoPos[i][0][0], &aminoPos[i][0][1], &aminoPos[i][0][2]);
      fscanf(f,"%lf %lf %lf\n", &aminoPos[i][1][0], &aminoPos[i][1][1], &aminoPos[i][1][2]);  
      i++;

    }
#ifdef PEPTIDE_PLATE
  for (i = 0; i < NAMINO-1; i++)
    {
      /* CA of aminoacid i */
      for (kk=0; kk < 3; kk++)
	platePos[i][0][kk] = aminoPos[i][0][kk];
      /* C of aminoacid i */
      for (kk=0; kk < 3; kk++)
	platePos[i][1][kk] = aminoPos[i][3][kk];
      /* CA of aminoacid i+1*/
      for (kk=0; kk < 3; kk++)
	platePos[i][2][kk] = aminoPos[i+1][0][kk];
      /* N of aminaocid i+1 */
      for (kk=0; kk < 3; kk++)
	platePos[i][3][kk] = aminoPos[i+1][2][kk];
    }
#endif
}

/* Allocate memory for a matrix of COORD_TYPE */
double** AllocMatR(int size1, int size2)
{
  double ** v;
  int k;
  v = (double**) malloc(size1 * sizeof(double*));
  v[0] = (double*) malloc(size1 * size2 * sizeof(double));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}
void print_matrix(char *txt, double M[3][3])
{
  int ii, kk;
  printf("%s\n", txt);
  printf("{");
  for (ii=0; ii < 3; ii++)
    {
      printf("{");
      for (kk=0; kk < 3; kk++)
	{
	  printf("%.15G", M[ii][kk]);
	  if (kk < 2)
	    printf(",");
	}
      printf("}");
      if (ii < 2)
	printf(",\n");
    }
  printf("}\n");
}
void tRDiagR(double M[3][3], double a, double b, double c, double **Ri)
{
  int na;
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = a;
  Di[1][1] = b;
  Di[2][2] = c;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}

void calcItensTot(double I[3][3], double RCM[3])
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
	for (i=0; i < Namino; i++)
	  {
	    ri[0] = rx[i]-RCM[0];
	    ri[1] = ry[i]-RCM[1];
	    ri[2] = rz[i]-RCM[2];
	    distSq = Sqr(ri[0])+Sqr(ri[1])+Sqr(ri[2]);
	    I[j][k] += massAmino*(((j==k)?distSq:0.0) - ri[j]*ri[k]);
	  }
      }
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
#ifdef PEPTIDE_PLATE
  /* moment of inertia of centers of mass */
  for (j=0; j < 3; j++)
    for (k=0; k < 3; k++)
      {
	for (i=0; i < Nplate; i++)
	  {
	    ri[0] = rxP[i]-RCM[0];
	    ri[1] = ryP[i]-RCM[1];
	    ri[2] = rzP[i]-RCM[2];
	    distSq = Sqr(ri[0])+Sqr(ri[1])+Sqr(ri[2]);
	    I[j][k] += massPlate*(((j==k)?distSq:0.0) - ri[j]*ri[k]);
	  }
      }
  /* moment of inertia with respect to centes of mass */
  for (i=0; i < Nplate; i++)
    {
      tRDiagR(Icom, Iplate[0], Iplate[1], Iplate[2], RiP[i]);
      for (j=0; j < 3; j++)
    	for (k=0; k < 3; k++)
	  {
	    I[j][k] += Icom[j][k]; 
	  }
      }	
#endif
}
void calcTotAngMom(double Mtot[3], double RCM[3], double VCM[3])
{
  int i, kk, ll;
  double A[3], B[3], C[3], Ia[3][3], om[3], Lrb[3];
  Mtot[0]=Mtot[1]=Mtot[2]=0.0;
  for (i=0; i < Namino; i++)
   {
     A[0] = vx[i]-VCM[0];
     A[1] = vy[i]-VCM[1];
     A[2] = vz[i]-VCM[2];
     B[0] = rx[i]-RCM[0];
     B[1] = ry[i]-RCM[1];
     B[2] = rz[i]-RCM[2];
     vectProdVec(B, A, C);
     om[0]=wx[i];
     om[1]=wy[i];
     om[2]=wz[i];
#if 1
     if (Iamino[0]==Iamino[1] && Iamino[1]==Iamino[2])
       {
	 for (kk=0; kk < 3; kk++)
	   {
	     Lrb[kk] = Iamino[kk]*om[kk];
	   }
       }
     else
       {
	 tRDiagR(Ia,Iamino[0],Iamino[1],Iamino[2],Ri[i]);
	 for (kk=0; kk < 3; kk++)
	   {
	     Lrb[kk] = 0;
	     for (ll=0; ll < 3; ll++)
	       {
		 Lrb[kk]+=Ia[kk][ll]*om[ll];
	       }
	   }
       }
     	 
#endif
     Mtot[0] += Lrb[0]+massAmino*C[0];
     Mtot[1] += Lrb[1]+massAmino*C[1];
     Mtot[2] += Lrb[2]+massAmino*C[2];
   }
#ifdef PEPTIDE_PLATE
  for (i=0; i < Nplate; i++)
   {
     A[0] = vxP[i]-VCM[0];
     A[1] = vyP[i]-VCM[1];
     A[2] = vzP[i]-VCM[2];
     B[0] = rxP[i]-RCM[0];
     B[1] = ryP[i]-RCM[1];
     B[2] = rzP[i]-RCM[2];
     vectProdVec(B, A, C);
     om[0]=wxP[i];
     om[1]=wyP[i];
     om[2]=wzP[i];
#if 1
     if (Iplate[0]==Iplate[1] && Iplate[1]==Iplate[2])
       {
	 for (kk=0; kk < 3; kk++)
	   {
	     Lrb[kk] = Iplate[kk]*om[kk];
	   }
       }
     else
       {
	 tRDiagR(Ia,Iplate[0],Iplate[1],Iplate[2],RiP[i]);
	 for (kk=0; kk < 3; kk++)
	   {
	     Lrb[kk] = 0;
	     for (ll=0; ll < 3; ll++)
	       {
		 Lrb[kk]+=Ia[kk][ll]*om[ll];
	       }
	   }
       }
     	 
#endif
     Mtot[0] += Lrb[0]+massPlate*C[0];
     Mtot[1] += Lrb[1]+massPlate*C[1];
     Mtot[2] += Lrb[2]+massPlate*C[2];
   }

#endif
}
int main(int argc, char** argv)
{
 FILE *f;
 double masstot;
 double RCM[3], VCM[3];
 double K, rr, PBdepth, mod, dirx, diry, com[3], R[3][3];
 double Vx, Vy, Vz, Mtot[3], ItensTot[3][3],
	Rx, Ry, Rz, xtmp[4][3], Itens[3][3], Iev[3],
	xlab[4][3];
 int i, jj, kk, a, ii;
 if (argc==1)
   {
     printf("You must supply the input file!\n");
     exit(-1);
   }
 srand(2);
 readSergeyPos(argv[1]);
 f=fopen("polyAla.cnf","w");
 init_parameters();
 /* i bond peptidici sono ognuno una particella */
#ifdef PEPTIDE_PLATE
 fprintf(f,"parnum:%d\n", 2*Namino-1); 
#else
 fprintf(f,"parnum:%d\n", Namino); 
#endif 
 fprintf(f,"totStep:%d\n", 20000);
 fprintf(f,"time:%.15G\n",0.0); 
 fprintf(f,"curStep:%d\n", 0); 
 fprintf(f,"P:%.15G\n", 1.0);
 fprintf(f,"T:%.15G\n", 1.0);
 fprintf(f,"rcut:%.15G\n", -1.0);
 fprintf(f,"equilibrat: %d\n", 0); 
 fprintf(f,"Dt:%.15G\n", 0.05); 
#ifdef PEPTIDE_PLATE
 fprintf(f,"ninters:%d\n", 14);
 fprintf(f,"ntypes:%d\n", 2); 
#else
 fprintf(f,"ninters:%d\n", 18);
 fprintf(f,"ntypes:%d\n", 1); 
#endif 
 fprintf(f,"@@@\n");
#ifdef PEPTIDE_PLATE
 fprintf(f,"%d %d\n", Namino, Namino-1);
#else
 fprintf(f,"%d\n", Namino);
#endif
 /* amminoacid types (Alanyn in this case) */ 
 fprintf(f, "%.15G %.15G %.15G\n", 0.01, 0.01, 0.01);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massAmino, Iamino[0], Iamino[1], Iamino[2], brownian, 1);
#ifdef PEPTIDE_PLATE 
 fprintf(f, "%d %d\n", 19, 0);
#else
 fprintf(f, "%d %d\n", 30, 0);
#endif
#ifdef TETRAHEDRON
 buildTetrahedra(x, tetraEdge);
 Dh = sqrt(6.0)/12.0*tetraEdge;/* =r; */
 Dh2 = (sqrt(3.0)/3.0 - sqrt(3.0)/6.0)*tetraEdge;/* = x-d */
#endif
 buildAminoSpots(0, x);
/* spots for steric hindrance */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigAA);/* CA (0) - 0 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigAC);/* CA (0) - 1 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigBA);/* CA (0) - 2 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigAN);/* CA (0) - 3 */

 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[1][0], x[1][1], x[1][2], sigBB);/* CB (1) - 4 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[1][0], x[1][1], x[1][2], sigBC);/* CB (1) - 5 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[1][0], x[1][1], x[1][2], sigBA);/* CB (1) - 6 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[1][0], x[1][1], x[1][2], sigBN);/* CB (1) - 7 */

 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigNN);/* N  (2) - 8 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigNC);/* N  (2) - 9 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigAN);/* N  (2) - 10 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigBN);/* N  (2) - 11*/

 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigCC);/* C  (3) - 12 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigNC);/* C  (3) - 13 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigAC);/* C  (3) - 14 */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigBC);/* C  (3) - 15 */
#ifdef PEPTIDE_PLATE
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigPepC);  /* 16 C */ 
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigPepCA);/* 17 CA */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigPepN); /* 18 N */
#else
 /* spots for peptide bond (centers coincide with previous spots) */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigCN1);  /* 16 C */ 
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigCN2);  /* 17 C */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigCCA1); /* 18 C */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[3][0], x[3][1], x[3][2], sigCCA2); /* 19 C */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigCACA1);/* 20 CA */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigCACA2);/* 21 CA */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigCAN1); /* 22 CA */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigCAN2); /* 23 CA */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigNCA1); /* 24 N */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigNCA2); /* 25 N */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigNC1);  /* 26 N */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[2][0], x[2][1], x[2][2], sigNC2);  /* 27 N */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigCAC1); /* 28 CA */
 fprintf(f, "%.15G %.15G %.15G %.15G\n", x[0][0], x[0][1], x[0][2], sigCAC2); /* 29 CA */
#endif
#ifdef PEPTIDE_PLATE
 buildPeptidePlate(xpep);
 fprintf(f, "%.15G %.15G %.15G\n", 0.01, 0.01, 0.01);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massPlate, Iplate[0], Iplate[1], Iplate[2], brownian, 1);
 fprintf(f, "%d %d\n", 4, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", xpep[0][0], xpep[0][1], xpep[0][2], sigPepCA);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", xpep[1][0], xpep[1][1], xpep[1][2], sigPepC);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", xpep[2][0], xpep[2][1], xpep[2][2], sigPepCA);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", xpep[3][0], xpep[3][1], xpep[3][2], sigPepN);
#endif
 
 /* all interactions */ 
 /* N.B. if barrier is higher than a certain threshold optimize bump routine! */ 
 /* hard core interactions */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 0, 0, 0,  0.0, 1E10, 0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 1, 0, 14, 0.0, 1E10, 0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 2, 0, 6,  0.0, 1E10, 0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 3, 0, 10, 0.0, 1E10, 0.0, 10);
 
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 4, 0, 4, 0.0, 1E10, 0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 5, 0, 15, 0.0, 1E10, 0.0, 10);
 //fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 6, 0, 2, 0.0, 1E10, 0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 7, 0, 11, 0.0, 1E10, 0.0, 10);

 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 8, 0, 8, 0.0, 1E10, 0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 9, 0, 13, 0.0, 1E10, 0.0, 10);
 //fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,10, 0, 3, 0.0, 1E10, 0.0, 10);
 //fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,11, 0, 7, 0.0, 1E10, 0.0, 10);

 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,12, 0, 12, 0.0, 1E10, 0.0, 10);
 //fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,13, 0, 9, 0.0, 1E10, 0.0, 10);
 //fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,14, 0, 1, 0.0, 1E10, 0.0, 10);
 //fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,15, 0, 5, 0.0, 1E10, 0.0, 10);

 PBdepth = 0.0001;
#ifdef PEPTIDE_PLATE
 /* C(Amino)-C(Plate)*/
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 16, 1, 1, PBdepth, 0.0,  1E10, 1);
 /* CA(Amino)-CA(Plate)*/
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 17, 1, 0, PBdepth,  0.0, 1E10, 1);
 /* CA(Amino)-CA(Plate) */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  17, 1, 2, PBdepth, 0.0,  1E10, 1);
 /* N(Amino)-N(Plate) */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  18, 1, 3, PBdepth,  0.0, 1E10, 1);
#else
 /* peptide covalent interactions (permanent) */
 /* C-N */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  16, 0, 26, 0.0, 1E10,  0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  17, 0, 27, PBdepth,  0.0, 1E10, 1);
 /* C-CA */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  18, 0, 28, 0.0, 1E10,  0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  19, 0, 29, PBdepth,  0.0, 1E10, 1);
 /* CA-CA */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  20, 0,  20, 0.0, 1E10,  0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0,  21, 0,  21, PBdepth,  0.0, 1E10, 1);
 /* CA-N */
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 22, 0, 24, 0.0, 1E10,  0.0, 10);
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 23, 0, 25, PBdepth,  0.0, 1E10, 1);
#endif
 fprintf(f, "@@@\n");
 /* positions of amminoacids */
 dirx = 1.0;
 diry = 0.0;
#ifdef PEPTIDE_PLATE
 Nplate = Namino-1;
 RiP = malloc(sizeof(double**)*Nplate);
 for (i=0; i < Nplate; i++)
   {
     RiP[i] = AllocMatR(3,3);
   }
 rxP = malloc(sizeof(double)*Nplate);
 ryP = malloc(sizeof(double)*Nplate);
 rzP = malloc(sizeof(double)*Nplate);
 vxP = malloc(sizeof(double)*Nplate);
 vyP = malloc(sizeof(double)*Nplate);
 vzP = malloc(sizeof(double)*Nplate);
 wxP = malloc(sizeof(double)*Nplate);
 wyP = malloc(sizeof(double)*Nplate);
 wzP = malloc(sizeof(double)*Nplate);
#endif
 Rx = Ry = Rz = 0.0;
 Ri = malloc(sizeof(double**)*Namino);
 for (i=0; i < Namino; i++)
   {
     Ri[i] = AllocMatR(3,3);
   }
 rx = malloc(sizeof(double)*Namino);
 ry = malloc(sizeof(double)*Namino);
 rz = malloc(sizeof(double)*Namino);
 vx = malloc(sizeof(double)*Namino);
 vy = malloc(sizeof(double)*Namino);
 vz = malloc(sizeof(double)*Namino);
 wx = malloc(sizeof(double)*Namino);
 wy = malloc(sizeof(double)*Namino);
 wz = malloc(sizeof(double)*Namino);
 masstot = 0;
 for (i=0; i < Namino; i++)
   {
#if 0
     mod =-L*0.5 + dL*i;
      if (i%2==0)
       //fprintf(f, "%.15G %.15G %.15G 1  0  0  0  1  0  0  0  1  0\n", mod*dirx, mod*diry, Dh); 
       fprintf(f, "%.15G %.15G %.15G 0  1  0  -1  0  0  0  0  1  0\n", mod*dirx, mod*diry, Dh); 
     else
       //fprintf(f, "%.15G %.15G %.15G 1  0  0  0  -1  0  0  0 -1  0\n",mod*dirx, mod*diry, -Dh); 
       fprintf(f, "%.15G %.15G %.15G 0  -1  0  -1  0  0  0  0 -1  0\n",mod*dirx, mod*diry+Dh2, -Dh); 
#else
     calcCOM(aminoPos[i],com); 
     /* calcola il tensore d'inerzia */
     /* i vettori riga di R sono i versori che individuano il 
	sistema di riferimento solidale con l'amminoacido cioè
	con il corpo rigido */
     for (a=0; a < 4; a++)
       for (kk=0; kk < 3; kk++)
	 xtmp[a][kk] = aminoPos[i][a][kk] - com[kk];
#if 1
     calcItens(xtmp, mass, Itens);
#else
     calcItens(xtmp, mass, Itens);
     if (i==0|| i==1)
       print_matrix("MYItens",Itens);
     calcItensSergey(xtmp, mass, Itens); 
#endif
     /* diagonalizza il tensore d'inerzia e ottiene gli autovettori */
     calcEigenValVect(Itens, R, Iev);
     check_eigenval(0, aminoPos[i], R, Iev, com, x);
     if (deterM(R) < 0.0)
       {
	 printf("[PLATE] i=%d detR=%.15G is less than 0!\n", i, deterM(R));
	 exit(-1);
       }

     for (ii=0; ii < 3; ii++)
       for (kk=0; kk < 3; kk++)
	 Ri[i][ii][kk] = R[ii][kk];
#ifdef DEBUG
      if (i==0 || i==1)
       {
	 printf("====> eigenvalues: %f %f %f\n", Iev[0], Iev[1], Iev[2]);
	 printf("COM=%.15G %.15G %.15G\n", com[0], com[1], com[2]);
	 printf("i=%d ------------ \n", i);
	 print_matrix("Itens",Itens);
	 print_matrix("Orientation Matrix",R);
       }
#endif
#if 0
     calcCOM(aminoPos[i], com);
     calcR(aminoPos[i], R); 
#endif
     rx[i] = com[0];
     ry[i] = com[1];
     rz[i] = com[2];
#ifdef DEBUG
     printf("i=%d ------------\n", i);
     for (jj=0; jj < 4; jj++)
       {
	 body2lab(x[jj], xlab[jj], com, R);
	 printf("xout=%.15G %.15G %.15G\n", xlab[jj][0], xlab[jj][1], xlab[jj][2]);
	 printf("sergey=%.15G %.15G %.15G\n", aminoPos[i][jj][0],aminoPos[i][jj][1],aminoPos[i][jj][2]);
	 lab2body(aminoPos[i][jj], xlab[jj], com, R);
	 printf("spots=%.15G %.15G %.15G\n", xlab[jj][0], xlab[jj][1], xlab[jj][2]);
       }
#endif
     masstot += massAmino;
     Rx += massAmino*com[0];
     Ry += massAmino*com[1];
     Rz += massAmino*com[2];
#if 0
     printf("scalProd R[2]*R[0]=%.15G\n", scalProd(R[2],R[0]));
     printf("scalProd R[1]*R[0]=%.15G\n", scalProd(R[1],R[0]));
#endif
#endif
   }
#ifdef PEPTIDE_PLATE
 /* --------------- PLATES POS BEG ------------ */
 for (i=0; i < Nplate; i++)
   {
     calcCOM(platePos[i],com); 
#ifdef DEBUGP
     for (kk=0; kk < 4; kk++)
       printf("[MAIN]platePos[%d][%d]:%f %f %f\n", i, kk, platePos[i][kk][0],platePos[i][kk][1],platePos[i][kk][2]);
#endif
     /* calcola il tensore d'inerzia */
     /* i vettori riga di R sono i versori che individuano il 
	sistema di riferimento solidale con l'amminoacido cioè
	con il corpo rigido */
     for (a=0; a < 4; a++)
       for (kk=0; kk < 3; kk++)
	 xtmp[a][kk] = platePos[i][a][kk] - com[kk];
#if 1
     calcItens(xtmp, massPlateSpots, Itens);
#else
     calcItensSergey(xtmp, massPlate, Itens); 
#endif
     /* diagonalizza il tensore d'inerzia e ottiene gli autovettori */
     calcEigenValVect(Itens, R, Iev);
   
     check_eigenval(1, platePos[i], R, Iev, com, xpep);
     if (deterM(R) < 0.0)
       {
	 printf("[PLATE] i=%d detR=%.15G is less than 0!\n", i, deterM(R));
	 exit(-1);
       }

#ifdef DEBUGP
      printf("PLATE i=%d ------------\n", i);
      for (jj=0; jj < 4; jj++)
	{
	  printf("i=%d COM=%f %f %f\n", i, com[0], com[1], com[2]);
	  print_matrix("[MAIN]Orientation Matrix",R);
	  body2lab(xpep[jj], xlab[jj], com, R);
	  printf("xpep[%d]=%f %f %f\n", jj, xpep[jj][0], xpep[jj][1], xpep[jj][2]);
	  printf("xout=%.15G %.15G %.15G\n", xlab[jj][0], xlab[jj][1], xlab[jj][2]);
	  printf("plate from sergey[%d]=%.15G %.15G %.15G\n", i,platePos[i][jj][0],platePos[i][jj][1],platePos[i][jj][2]);
	  printf("sergey[%d]=%.15G %.15G %.15G\n", i,aminoPos[i][jj][0],aminoPos[i][jj][1],aminoPos[i][jj][2]);
	  if (i < NAMINO)
	    printf("sergey[%d]=%.15G %.15G %.15G\n",i+1,aminoPos[i+1][jj][0],aminoPos[i+1][jj][1],aminoPos[i+1][jj][2]);
	  //lab2body(aminoPos[i][jj], xlab[jj], com, R);
	 //printf("spots=%.15G %.15G %.15G\n", xlab[jj][0], xlab[jj][1], xlab[jj][2]);
       }
#endif
     for (ii=0; ii < 3; ii++)
       for (kk=0; kk < 3; kk++)
	 RiP[i][ii][kk] = R[ii][kk];
#ifdef DEBUGP
      if (1 || i==0 || i==1)
       {
	 printf("i=%d BEG PLATE ====> eigenvalues: %f %f %f\n", i, Iev[0], Iev[1], Iev[2]);
	 printf("COM=%.15G %.15G %.15G\n", com[0], com[1], com[2]);
	 print_matrix("Itens",Itens);
	 //print_matrix("Orientation Matrix",R);
	 printf("===================== END ================ \n");
       }
#endif
      rxP[i] = com[0];
      ryP[i] = com[1];
      rzP[i] = com[2];
      masstot += massPlate;
      Rx += massPlate*com[0];
      Ry += massPlate*com[1];
      Rz += massPlate*com[2];
   }
 Rx /= masstot;
 Ry /= masstot;
 Rz /= masstot;

  /* -------------- PLATES ---------------- */
#else
 Rx /= masstot;
 Ry /= masstot;
 Rz /= masstot;
#endif
 for (i=0; i < Namino; i++)
   {
     rx[i] -= Rx;
     ry[i] -= Ry;
     rz[i] -= Rz;
     fprintf(f, "%.15G %.15G %.15G %.15G %.15G  %.15G %.15G %.15G %.15G %.15G %.15G %.15G 0\n", rx[i], ry[i], rz[i],
	     Ri[i][0][0], Ri[i][0][1], Ri[i][0][2], Ri[i][1][0], Ri[i][1][1], 
	     Ri[i][1][2], Ri[i][2][0],Ri[i][2][1],Ri[i][2][2]); 
   } 
#ifdef PEPTIDE_PLATE
 for (i=0; i < Nplate; i++)
   {
     rxP[i] -= Rx;
     ryP[i] -= Ry;
     rzP[i] -= Rz;
     /* fix quick and dirty: we now that orientations of plate are all equal */
     fprintf(f, "%.15G %.15G %.15G %.15G %.15G  %.15G %.15G %.15G %.15G %.15G %.15G %.15G 1\n", rxP[i], ryP[i], rzP[i],
	     RiP[i][0][0], RiP[i][0][1], RiP[i][0][2], RiP[i][1][0], RiP[i][1][1], 
	     RiP[i][1][2], RiP[i][2][0], RiP[i][2][1], RiP[i][2][2]); 
   }
#endif
  /* velocities and angular velocities of aminoacids */
 K = sqrt(T/massAmino);
 Mtot[0] = Mtot[1] = Mtot[2] = Vx = Vy = Vz = 0;
 for (i=0; i < Namino; i++)
   {
     vx[i] =  K*gauss();
     vy[i] =  K*gauss();
     vz[i] =  K*gauss();
     Vx += massAmino*vx[i];
     Vy += massAmino*vy[i];
     Vz += massAmino*vz[i];
     angvel(i, &wx[i], &wy[i], &wz[i], T);
   }
#ifdef PEPTIDE_PLATE
 for (i=0; i < Nplate; i++)
   {
     vxP[i] =  K*gauss();
     vyP[i] =  K*gauss();
     vzP[i] =  K*gauss();
     Vx += massPlate*vxP[i];
     Vy += massPlate*vyP[i];
     Vz += massPlate*vzP[i];
     angvel(i, &wxP[i], &wyP[i], &wzP[i], T);
   }
 Vx /= masstot;
 Vy /= masstot;
 Vz /= masstot;
#else 
 Vx /= masstot;
 Vy /= masstot;
 Vz /= masstot;
#endif 
 RCM[0] = 0.0;
 RCM[1] = 0.0;
 RCM[2] = 0.0;
 VCM[0] = Vx;
 VCM[1] = Vy;
 VCM[2] = Vz;
 for (i=0; i < Namino; i++)
   {
     vx[i] -= VCM[0];
     vy[i] -= VCM[1];
     vz[i] -= VCM[2];
   }
#ifdef PEPTIDE_PLATE
 for (i=0; i < Nplate; i++)
   {
     vxP[i] -= VCM[0];
     vyP[i] -= VCM[1];
     vzP[i] -= VCM[2];
   }

#endif
 VCM[0]=VCM[1]=VCM[2]=0.0;
 calcTotAngMom(Mtot, RCM, VCM);
 printf("[Prima]Total Angular Momentum= %.15G %.15G %.15G\n", Mtot[0], Mtot[1], Mtot[2]);
 calcItensTot(ItensTot, RCM); 
 for (kk = 0; kk < 3; kk++)
   omega[kk] = Mtot[kk];
 SolveLineq(ItensTot, omega);
   
 printf("wxr=%.15G %.15G %.15G  omega= %.15G %.15G %.15G\n", C[0], C[1], C[2], omega[0], omega[1], omega[1]);
 for (i=0; i < Namino; i++)
   {
     B[0] = rx[i];
     B[1] = ry[i];
     B[2] = rz[i];
     vectProdVec(omega, B, C);
     vx[i] -= C[0];
     vy[i] -= C[1];
     vz[i] -= C[2];
     wx[i] -= omega[0];
     wy[i] -= omega[1];
     wz[i] -= omega[2];
     fprintf(f, "%.15G %.15G %.15G %.15G %.15G %.15G\n", vx[i], vy[i], vz[i], 
	     wx[i], wy[i], wz[i]);
   } 
#ifdef PEPTIDE_PLATE
 for (i=0; i < Nplate; i++)
   {
     B[0] = rxP[i];
     B[1] = ryP[i];
     B[2] = rzP[i];
     vectProdVec(omega, B, C);
     vxP[i] -= C[0];
     vyP[i] -= C[1];
     vzP[i] -= C[2];
     wxP[i] -= omega[0];
     wyP[i] -= omega[1];
     wzP[i] -= omega[2];
     fprintf(f, "%.15G %.15G %.15G %.15G %.15G %.15G\n", vxP[i], vyP[i], vzP[i], 
	     wxP[i], wyP[i], wzP[i]);
   } 
#endif
 calcTotAngMom(Mtot, RCM, VCM);
 printf("[dopo]Total Angular Momentum= %.15G %.15G %.15G\n", Mtot[0], Mtot[1], Mtot[2]);
 fprintf(f, "%.15G\n", L);
 fclose(f);
}

