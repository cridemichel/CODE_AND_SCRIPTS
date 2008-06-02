#include <stdio.h> /* standart inputoutput */
#include <ctype.h> /* for charachter recognition */
#include <stdlib.h> /* for conversion from char to dec */
#include <strings.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <sys/types.h>
#include <inttypes.h>
#include <errno.h>
#undef DEBUG
#undef DEBUGP
#define NAMINO 60
/* DIAMETRI DELLE PATCH (Unita' di misura in Angstrom) */
#define DP 2.0*0.11965683746373795116
#define DPIV 4.0*0.11965683746373795116
#define DREP 9.0
#define DIAM_AMINO 4.0
#define NUM_PATCH 7
double sigPatch[NUM_PATCH]={DP,DP,DP,DP,DPIV,DPIV,DREP};
typedef struct goPair {
int i;
int j;
double u0;
} goPairType;
goPairType *goArr;
int ntypesOfAmino[20];
int *typeOfAmino;
int ngointer=0;
const double infbarr = 1E10;
double *rx, *ry, *rz, *vx, *vy, *vz, *wx, *wy, *wz, ***Ri, omega[3];
double eigenVect[3][3];
int Namino;
int brownian, ntypes, ninters;
double fact, temp, L, massAmino, T, Iamino[3], x[4][3], tetraEdge, dL, PBw, Dh, Dh2;
double sigCA, sigBC, sigNC, sigCC, sigAC, sigAN, sigNN, sigBA, sigBN, sigAA, sigBB,
       sigCN1, sigCN2, sigCCA1, sigCCA2, sigCACA1, sigCACA2, sigCAN1, sigCAN2,
       sigNCA1, sigNCA2, sigNC1, sigNC2, sigCAC2, sigCAC2, sigCAC1, sigCAC2, sigCB2GM,
       A[3], B[3], C[3], delCANbond, delCCAbond; 
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
  L=250.0;
  massAmino = 1.0;
  Iamino[0] = Iamino[1] = 1.0;
  Iamino[2] = 1.0;
  brownian = 0;
  T=0.1;	
  ntypes=20;
}
double gauss(void);

void angvel(int i, double *wx, double *wy, double* wz, double temp)
{
  double mean, inert;                 /* momentum of inertia of the molecule */
  inert = Iamino[0];
  mean = sqrt(temp / inert);
  *wx = mean*gauss();
  *wy = mean*gauss();
  *wz = mean*gauss();
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
/* GENERA LE POSIZIONI DELLE PATCH */
double pos_spot[NUM_PATCH][3]={{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}};
void buildAminoSpots(double xout[NUM_PATCH][3])
{
  int a, kk, jj, jj2;
      
  for (a=0; a < NUM_PATCH; a++)
    for (kk=0; kk < 3; kk++)
      xout[a][kk] = pos_spot[a][kk]*DIAM_AMINO*0.5;
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
     Mtot[0] += Lrb[0]+massAmino*C[0];
     Mtot[1] += Lrb[1]+massAmino*C[1];
     Mtot[2] += Lrb[2]+massAmino*C[2];
   }
}
char line[1024];
int main(int argc, char** argv)
{
 FILE *f, *ff;
 int seqlen;
 double masstot;
 double RCM[3], VCM[3];
 double K, rr, PBdepth, mod, dirx, diry, com[3], R[3][3];
 double Vx, Vy, Vz, Mtot[3], ItensTot[3][3],
	Rx, Ry, Rz, xtmp[4][3], Itens[3][3], Iev[3],
	xlab[4][3];
 int i, jj, kk, a, ii, GMspot;
#if 0
 if (argc < 3)
   {
     printf("You must supply the input file and the go matrix file!\n");
     exit(-1);
   }
#endif
 srand(2);
 f=fopen("ivan.cnf","w");
 init_parameters();
 printf("Reading sequence file...");
 ff = fopen("sequence.dat","r");
 fscanf(ff, "%d ", &seqlen);
 Namino = seqlen;
 typeOfAmino = (int *) malloc(sizeof(int)*seqlen);
 for (i = 0; i < seqlen; i++)
   {
     fscanf(ff, "%d ", &(typeOfAmino[i]));
   }
 fclose(ff);
 for (i = 0; i < 20; i++)
   ntypesOfAmino[i]=0;
 for (i = 0; i < seqlen; i++)
   { 
     typeOfAmino[i]-=1;
     ++(ntypesOfAmino[typeOfAmino[i]]);
   }
 //printf("===>T=%.15G seqlen=%d",T, seqlen);
 printf("...done\n");
 ff=fopen("Interactions.dat","r");
 ninters = 0;
 while (!feof(ff))
   {
     fscanf(ff, "%[^\n] ", line);
     ninters++;
   }
 fclose(ff);/* i bond peptidici sono ognuno una particella */
 fprintf(f,"parnum:%d\n", Namino); 
 fprintf(f,"totStep:%d\n", 20000);
 fprintf(f,"time:%.15G\n",0.0); 
 fprintf(f,"curStep:%d\n", 0); 
 fprintf(f,"P:%.15G\n", 1.0);
 fprintf(f,"T:%.15G\n", 1.0);
 fprintf(f,"rcut:%.15G\n", -1.0);
 fprintf(f,"equilibrat: %d\n", 0); 
 fprintf(f,"Dt:%.15G\n", 0.05); 
 fprintf(f,"ninters:%d\n", ninters);
 fprintf(f,"nintersIJ:%d\n", ngointer);
 fprintf(f,"ntypes:%d\n", ntypes); 
 fprintf(f,"@@@\n");

 /* numero di amminoacidi di ogni tipo */
 for (i=0; i < 20; i++)
   fprintf(f,"%d ", ntypesOfAmino[i]);
 fprintf(f, "\n");
 /* geometria dei 20 amminoacidi */
 for (i=0; i < 20; i++)
   {
     fprintf(f, "%.15G %.15G %.15G\n", DIAM_AMINO*0.5, DIAM_AMINO*0.5, DIAM_AMINO*0.5);
     fprintf(f, "%d %d %d\n", 1, 1, 1);/* ignore*/
     fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massAmino, Iamino[0], Iamino[1], Iamino[2], brownian, 1);
     /* numero di patch del modello per ogni amminoacido */
     fprintf(f, "%d %d\n", NUM_PATCH, 0);
     buildAminoSpots(x);
     /* spots for steric hindrance */
     for (a = 0; a < NUM_PATCH; a++)
       fprintf(f, "%.15G %.15G %.15G %.15G\n", x[a][0], x[a][1], x[a][2], sigPatch[a]);
   }
  /* PER IVAN: metti qui la tua routine che genera le interazioni */
 printf("Reading interaction file...");
 ff=fopen("Interactions.dat","r");
 ninters = 0;
 while (!feof(ff))
   {
     fscanf(ff, "%[^\n] ", line);
     fprintf(f,"%s\n", line);
     ninters++;
   }
 fclose(ff);
 printf("...done\n");
 /* qui finisce la geometria */
 fprintf(f, "@@@\n");
 /* positions of amminoacids */
 dirx = 1.0;
 diry = 0.0;
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
     rx[i] = -L*0.5+DIAM_AMINO*0.5+5.0+i*(DIAM_AMINO*(1+0.01));
     ry[i] = 0.0;
     rz[i] = 0.0;
     fprintf(f, "%.15G %.15G %.15G ", rx[i], ry[i], rz[i]);
     fprintf(f, "1 0 0 0 1 0 0 0 1 %d\n", typeOfAmino[i]);
     masstot += massAmino;
     Ri[i][0][0] = 1.0;
     Ri[i][0][1] = 0.0;
     Ri[i][0][2] = 0.0;
     Ri[i][1][0] = 0.0;
     Ri[i][1][1] = 1.0;
     Ri[i][1][2] = 0.0;
     Ri[i][2][0] = 0.0;
     Ri[i][2][1] = 0.0;
     Ri[i][2][2] = 1.0;
     Rx += massAmino*rx[i];
     Ry += massAmino*ry[i];
     Rz += massAmino*rz[i];
   }
 Rx /= masstot;
 Ry /= masstot;
 Rz /= masstot;
 for (i=0; i < Namino; i++)
   {
     rx[i] -= Rx;
     ry[i] -= Ry;
     rz[i] -= Rz;
   } 
 /* velocities and angular velocities of aminoacids */
 K = sqrt(T/massAmino);
 printf("T=%.15G massAmin: %.15G\n", T, massAmino);
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
 Vx /= masstot;
 Vy /= masstot;
 Vz /= masstot;

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
 calcTotAngMom(Mtot, RCM, VCM);
 printf("[dopo]Total Angular Momentum= %.15G %.15G %.15G\n", Mtot[0], Mtot[1], Mtot[2]);
 fprintf(f, "%.15G\n", L);
 fclose(f);
 return 0;
}

