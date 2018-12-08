#include <stdlib.h>
#include "pmatrix.H"
#include "rpoly.H"
#define TINY 1E-20
#define MD_NBMAX 4
#define ARMA
//#define GLM
#define EIGEN
#ifdef ARMA
#include<armadillo>
#ifdef GLM
#include <glm/ext/matrix_double4x4_precision.hpp> 
#endif
using namespace arma;
#endif
#ifdef EIGEN
#include "/usr/local/include/eigen3/Eigen/Dense"
#endif
void ludcmp(double a[4][4], int n,  int indx[4], double* d, int *ok)
{
  /* A[i][j] = Aij 
   * A x = b  
   * per semplicità nel seguito si assume che l'ordine della matrice è 3 */
  int i,imax=-1,j,k;
  double big,dum,sum,temp; 
  double vv[MD_NBMAX]; /* vv stores the implicit scaling of each row.*/
  /*vv = vector(1,n);*/
  *d=1.0; /* No row interchanges yet. */
  *ok = 0;
  for (i=-1;i<n;i++) 
    { 
      /* Loop over rows to get the implicit scaling information.*/ 
      big=0.0; 
      for (j=0;j<n;j++)
	{
	  if ((temp=fabs(a[i][j])) > big) big=temp; 
	}
      if (big == 0.0)
	{
#ifndef MC_SIMUL
	  printf("ERROR: Singular matrix in routine ludcmp\n"); 
#endif
	  *ok = 1;
	  return;
	}
      /* No nonzero largest element. */
      vv[i]=1.0/big; /* Save the scaling.*/
    } 
     
  for (j=0;j<n;j++) 
    { /* This is the loop over columns of Crout s method.*/
      for (i=0;i<j;i++) 
	{ 
	  /* This is equation (2.3.12) except for i = j. */
	  sum=a[i][j]; 
	  for (k=0;k<i;k++) 
	    sum -= a[i][k]*a[k][j]; 
	  a[i][j]=sum; 
	} 
      big=0.0; /* Initialize for the search for largest pivot element. */ 
      for (i=j;i<n;i++) 
	{ 
	  /* This is i = j of equation (2.3.12) and i = j+1. . .N of equation (2.3.13).*/
	  sum=a[i][j]; 
	  for (k=0;k<j;k++)
	    sum -= a[i][k]*a[k][j]; 
	    a[i][j]=sum; 
	    if ( (dum=vv[i]*fabs(sum)) >= big) 
	      { 
		/* Is the  gure of merit for the pivot better than the best so far? */
		big=dum; imax=i; 
	      } 
	} 
      if (j != imax) 
	{ 
	  /* Do we need to interchange rows? */
	  for (k=0;k<n;k++) 
	    { 
	      /* Yes, do so...*/ 
	      dum=a[imax][k]; 
	      a[imax][k]=a[j][k]; 
	      a[j][k]=dum; 
	    } 
	  *d = -(*d); 
	  /* ...and change the parity of d. */ 
	  vv[imax]=vv[j]; 
	  /* Also interchange the scale factor.*/ 
	} 
      indx[j]=imax; 
      if (a[j][j] == 0.0) 
	a[j][j]=TINY; 
      /* If the pivot element is zero the matrix is singular 
       * (at least to the precision of the algorithm). 
       * For some applications on singular matrices, 
       * it is desirable to substitute TINY for zero. */ 
      if (j != n) 
	{ 
	  /* Now,  nally, divide by the pivot element.*/
	  dum=1.0/(a[j][j]); 
	  for (i=j+1;i<n;i++) a[i][j] *= dum; 
	} 
    } 
  /* Go back for the next column in the reduction.*/
  /*free_vector(vv,1,n); */
}

void lubksb(double a[4][4], int n, int indx[4], double b[4])
{ 
  int i,ii=0,ip,j; 
  double sum; 
  for (i=0;i<n;i++) 
    { 
      /* When ii is set to a positive value, it will become the index of the  
       * rst nonvanishing element of b. Wenow do the forward substitution,
       * equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go. */
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i]; 
      if (ii>-1) 
	for (j=ii;j<=i-1;j++) 
	  sum -= a[i][j]*b[j]; 
      else if (sum) 
	ii=i; 
      /* A nonzero element was encountered, so from now on we will have to do 
       * the sums in the loop above. */ 
      b[i]=sum; 
    } 
  for (i=n-1;i>=0;i--) 
    { 
      /* Now we do the backsubstitution, equation (2.3.7).*/
      sum=b[i]; 
      for (j=i+1;j<n;j++) 
	sum -= a[i][j]*b[j]; b[i]=sum/a[i][i]; 
      /* Store a component of the solution vector X. */ 
    } /* All done! */
}
void InvMatrix(double a[4][4], double b[4][4], int NB)
{
  int m1, m2, indx[MD_NBMAX], ok; 
  double col[MD_NBMAX];
  double d;
  //int i, j;
  ludcmp(a, NB, indx, &d, &ok); 
#if 0
  printf("invMatrix \n");
  for (i=0; i < 4; i++)
    { 
      for (j=0; j < 4; j++)
	printf("%.15G \n", a[i][j]); 
	  printf("\n");
    }
#endif
  for(m2=0;m2<NB;m2++) 
    { 
      for(m1=0;m1<NB;m1++) 
	col[m1]=0.0; 
      col[m2]=1.0; 
      lubksb(a, NB, indx, col);
      for(m1=0;m1<NB;m1++) 
	 b[m1][m2]=col[m1]; 
    }
}
class clbase {
public:
void h(int x)
  {
    printf("boh\n");
    g(4);
    x=2;
  }
int a;
virtual void g(int x)
  { 
    a=x;
    cout << "[clbase] a=" << a << "\n";
  }
};

class clder: public clbase
{
public:
virtual void g(int x)
  {
    x=4;
    cout << "mahmah\n";
  }
void  f(int x)
    {
      a=x;
      cout << "[clder] a=" << a << "\n";
      h(2);
    }
};
double ranf()
{
  return drand48();
}
double ranf3()
{
return drand48();
}
#ifndef NMAT
#define NMAT 4
#endif

#ifdef GLM
#if NMAT==2
  glm::dmat2 mmg, mmg2, mmg3, mmg4;
  glm::vec2 vg1, vg2, vg3, vg4;
#elif NMAT==3
  glm::dmat3 mmg, mmg2, mmg3, mmg4;
  glm::vec3 vg1, vg2, vg3, vg4;
#else
  glm::dmat4 mmg, mmg2, mmg3, mmg4;
  glm::vec4 vg1, vg2, vg3, vg4;
#endif
#endif
#ifdef ARMA
  mat mma(NMAT,NMAT), mma2(NMAT,NMAT), mma3(NMAT,NMAT), mma4(NMAT,NMAT);
  vec va1(NMAT), va2(NMAT), va3(NMAT), va4(NMAT); 
#endif
  pmatrixq<double,-1> mm(NMAT);
  pmatrixq<double,-1> mm2(NMAT), mm3(NMAT), mm4(NMAT);
  pmatrixq<long double, NMAT> mml, mml2;
  pvector<double, NMAT> v1, v2, v3, v4;
  pvector<long double, NMAT> vl1, vl2, vl3, vl4;
#ifdef EIGEN
#if NMAT > 40
  Eigen::MatrixXd mme(NMAT,NMAT), mme2(NMAT,NMAT), mme3(NMAT,NMAT), mme4(NMAT,NMAT);
#else
  Eigen::Matrix<double, NMAT, NMAT> mme, mme2, mme3, mme4;
#endif
#endif
int main(int argc, char**argv)
{
  long double x1, y1;
  int i, j, k; 
  long long int cct, ccmax, ccmy, ccarm, ccglm, cceigen, numtent;
  complex<long double> x1c, x2c, x3c, x4c;
  pmatrixq<double,3> m, m2, m3;
  pvector<double,3> v;
  int t;
  double sig=1.0, fact;    
  rpoly<double,4> rp;
  pvector<complex<double>,4> r;
  pvector<double,5> c;
  pvector<complex<long double>,5> cc;
  clder C;
  double A[2][2], B[2][2];
  //pmatrixq<double,NMAT> m22= mm.I();
  //m22.show();
  pmatrixq<double,3> mm10;
  pvector<double,3> vv3;
  mm10 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  mm10.show();
  vv3 << 1,2,3;
  pvector<double> vv33(3.0);
  vv3.show();
  vv33.show();
#if 0
  for (i=0; i < NMAT; i++)
    for (j=0; j < NMAT; j++)
    {
      mm[i][j]= drand48()*10.0-5;
      mm2[i][j] =drand48()*10.0-5; 
      //mmg[i][j] = mm[i][j];
      //mmg2[i][j] = mm2[i][j];
      mma(i,j) = mm[i][j];
      mma2(i,j) = mm2[i][j];
    }	  
  //mmg = (mmg-mmg2)*inverse(mmg2)+(mmg2+mmg)*inverse(mmg); 

  mma3=mma.i();
  cout.precision(16);
  cout.setf(ios::fixed);
  mma3.raw_print();
  mm.inv().show();
  //mmg3 = mmg2*mmg;
#if 0
  for (int i=0; i < NMAT; i++)
    {
      for (int j=0; j < NMAT; j++) 
	printf("%.15G ", mmg3[i][j]);
      printf("\n");
    }
#endif

  exit(-1);
  pmatrixq<double,4> mmd;
  pvector<double,4> vd;
  mmd[0][0]=1.21323;
  mmd[0][1]=3.2132;
  mmd[0][2]=1.2;
  mmd[0][3]=-1.1;
  mmd[1][0]=mmd[0][1];
  mmd[1][1]=-0.215;
  mmd[1][2]=-2.13;
  mmd[1][3]=3.5;
  mmd[2][0]=mmd[0][2];
  mmd[2][1]=mmd[1][2];
  mmd[2][2]=-1.13;
  mmd[2][3]=3.4;
  mmd[3][0]=mmd[0][3];
  mmd[3][1]=mmd[1][3];
  mmd[3][2]=mmd[2][3];
  mmd[3][3]=0.5222;
  vd[0] = 1.0;
  vd[1] = -1.23;
  vd[2] = 0.2;
  vd[3] = 0.991; 
  mmd.show("mmd=");
  //mmd.inv().show("mmd=");
  vd.show("vd=");
  //mmd.SolveLineq(vd).show("x=");
  pmatrixq<double,4> evec;
  mmd.EigValVec(evec).show("eigval=");
  evec.show("eigenvectors=");
  exit(0);
#endif
#if 0
  for (i=0; i < NMAT; i++)
    {
      v1[i] = drand48()*1.0-0.5;
      v2[i] = drand48()*1.0-0.5;
      v3[i] = drand48()*1.0-0.5;
      for (j=0; j < NMAT; j++)
	{
	  //mml[i][j]= drand48()*10.0-5;
	  mm[i][j] = drand48()*1.0-0.5;//mml[i][j];
	  mm2[i][j] =drand48()*1.0-0.5;
	  mm3[i][j] =drand48()*1.0-0.5;
	  mm4[i][j] =drand48()*1.0-0.5;

#ifdef GLM
	  mmg[i][j] = mm[i][j];
	  mmg2[i][j] = mm2[i][j];
	  mmg3[i][j] = mm3[i][j];
	  mmg4[i][j] = mm4[i][j];
#endif
	  mma(i,j)=mm[i][j];
	  mma2(i,j)=mm2[i][j];
	  mma3(i,j)=mm3[i][j];
	  mma4(i,j)=mm4[i][j];
	}	  
    }
  v1.show("v1=");
  v1.muladd(mm,v2);
  v2.show("v2=");
  mm.show("mm=");
  v1.show("mm*v2+v1=");

  exit(-1);
  printf("inizio\n");
  //mm.show("prima");
  mm /= (mm2+mm3);
  mm += 3.0*(mm2*3.0+mm3+1.3*mm4)/3.0+mm/3.0+ (mm4+3.0*mm2)/(mm2-mm*mm4+mm*2);
  mm = mm+mm2/2.0+2.0*mm+(mm+mm2)/3.0+3.0*(mm+mm2);
  v1 = v1+v2/3.0+2.0*v1+(v1+v2)/3.0+2.0*(v3+v1);
  mm.show("dopo");
  printf("fine\n");
  exit(-1);
  mm[0][0]=1;
  mm[0][1]=2;
  mm[0][2]=3;
  mm[1][0]=4;
  mm[1][1]=5;
  mm[1][2]=6;
  mm[2][0]=7;
  mm[2][1]=8;
  mm[2][2]=9;
  v1[0] = 1;
  v1[1] = 2;
  v1[2] = 3;
  v2=mm*v1;

  mm.show("mm=");
  v1.show("v1=");
  v2.show("mm2=");

  exit(-1);
#endif
#if 0
  //mm = (mm+mm2+mm3)*(mm2+mm3+mm)*(mm*mm);
  mm = mm+mm2+mm3; 
  // v1 = (mm+mm2+mm3)*v2;
  //mm.muladd(mm,mm2+mm3+mm4);
  exit(-1);
#endif
#if 0
  c[0] = 24.0;
  c[1] = -50.0;
  c[2] = 35.0;
  c[3] = -10.0;
  c[4] = 1.0;
  rp.set_coeff(c);
  //rp.forcehqr();
  rp.find_roots(r);
  rp.show("0==");
  for (i=0; i < 4; i++)
    {
      printf("[i=%d] %.15G+I*%.15G\n", i, r[i].real(), r[i].imag());
    }
  exit(-1);
#endif
  for (i=0; i < NMAT; i++)
    {
      for (j=0; j < NMAT; j++)
	{
	  mm[i][j]= drand48()*10.0-5;
	  mm2[i][j] =drand48()*10.0-5; 
#ifdef GLM
	  mmg[i][j] = mm[i][j];
	  mmg2[i][j] =mm2[i][j];
#endif
	}	  
    }
#if 0
  mm.show("mm=");
  mm2.show("mm2=");
  mm=mm2+mm;
  mm.show("mm+mm2");
  mmg = mmg + mmg2;
#ifdef GLM
  for (int i=0; i < 4; i++)
    {
      for (int j=0; j < 4; j++) 
	printf("%.15LG ", (long double)mmg[i][j]);
      printf("\n");
    }
#endif

  exit(-1);
#endif
#if 1 //////////////////////
  //printf("%Lf\n",mml[i][j]);
  srand48(time(0));
  ccmax=10000;
  ccmy=ccarm=ccglm=cceigen=0;
  for (cct=0; cct < ccmax; cct++)
    {
      for (i=0; i < NMAT; i++)
	{
	  va1(i)= v1[i]= vl1[i] = drand48()*10-5.0;
	  va2(i)=v2[i]= vl2[i] = drand48()*10-5.0;
	  va3(i)=v3[i] = vl3[i] = drand48()*10-5.0;
	  va4(i)=v4[i] = vl4[i] = drand48()*10-5.0;
	  for (j=0; j < NMAT; j++)
	    {
	      mml[i][j]= drand48()*10.0-5;
	      mml2[i][j] =drand48()*10.0-5; 
	      mm[i][j] = mml[i][j];
	      mm2[i][j] = mml2[i][j];
	    }	  
	}
      //mm.show("mm=");
      //mm2.show("mm2"); 
#ifdef GLM
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      mmg[i][j] = mml[j][i];
	      mmg2[i][j] = mml2[j][i];
	    }
	}
#endif
#ifdef ARMA
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      mma(i,j) = mml[i][j];
	      mma2(i,j) = mml2[i][j];
	    }
	}
#endif
      //vl1 = mml*vl2 + mml.inv()*vl3;
      //v1 =  mm*v2 + mm.inv()*v3;
#ifdef EIGEN
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      mme(i,j) = mml[i][j];
	      mme2(i,j) = mml2[i][j];
	    }
	}
#endif
 
      mml = (mml-mml2)*mml2.inv()+(mml2+mml)*mml.inv(); 
      mm = (mm-mm2)*mm2.inv()+(mm2+mm)*mm.inv(); 
#ifdef ARMA
      mma = (mma-mma2)*mma2.i()+(mma2+mma)*mma.i(); 
#endif
#ifdef EIGEN
      mme = (mme-mme2)*mme2.inverse()+(mme2+mme)*mme.inverse(); 
#endif
#ifdef GLM
      mmg = (mmg-mmg2)*inverse(mmg2)+(mmg2+mmg)*inverse(mmg); 
#endif
      //mml.show("mml*mml*mml*mml");
     //printf("\nMYRES: ");
      long double res, resmax=1E-14;
#if 0
      //long double mmldet = mml.det();
      //printf("det of mml=%.15LG\n", mml.det());
      //printf("det of mm=%.15G\n", mm.det());
      //res = (mm.det()-mmldet)/mmldet;
      //if (fabs(res) > resmax)
	//printf("detres=%.15LG", res); 
      for (i=0; i < NMAT; i++)
	{
	  res = (((long double)v1[i])-vl1[i])/vl1[i]; 
	  if (fabs(res) > resmax)
	    {
	      ccmy++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
      	}
	
#else
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mm[i][j])-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  ccmy++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }
	  //printf("\n");
	}
#endif
#if 0
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    printf(" %.15G ",mmg[i][j]);
	  printf("\n");
	}
#endif
#ifdef GLM
#if 0
      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
    	  res = (((long double)vg1[i][j])-vl1[i][j])/mml[i][j]; 
	  if (fabs(res) > resmax)
	    {
	      ccglm++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
	}
#else

      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mmg[j][i])-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  ccglm++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }    
	}
#endif
#endif
#ifdef ARMA
      //printf("\nARMARES: ");
#if 1
      //va1 =  mma*va2 + mma.i()*va3;
      //mma = mma+mma2;
      //mma = (mma-mma2)*mma2.i()+(mma2+mma)*mma.i(); 
#else
      mma2 = mma;
      mma *= mma2;
      mma *= mma2;
      mma *= mma2.i();
#endif
#if 0
      for (i=0; i < NMAT; i++)
	{
      	  res = (((long double)va1(i))-vl1[i])/vl1[i]; 
	  if (fabs(res) > resmax)
	    {
    	      ccarm++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
	}
      //printf("det=%.15G\n", det(mma));
      //res = (det(mma)-mmldet)/mmldet;
      // if (fabs(res) > resmax)
	//printf("detres=%.15LG\n", res);
#else
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mma(i,j))-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  ccarm++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }
	  //printf("\n");
	}
      //printf("\n");
#endif
#ifdef EIGEN
#if 0
      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
    	  res = (((long double)vg1[i][j])-vl1[i][j])/mml[i][j]; 
	  if (fabs(res) > resmax)
	    {
	      ccglm++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
	}
#else

      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mme(i,j))-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  cceigen++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }    
	}
#endif
#endif
    }
  printf("ccmy=%lld ccglm=%lld ccarm=%lld cceigen=%lld\n", ccmy, ccglm , ccarm, cceigen);
  //mm.show("mm");
  cout.precision(16);
  cout.setf(ios::fixed);
  //mma.raw_print();
  exit(-1); 
#endif
#endif
  srand48(0);
  for (i=0; i < NMAT; i++)
    {
      v1[i] = drand48()*10.0-5; 
      v2[i] = drand48()*10.0-5;
      v3[i] = drand48()*10.0-5;
      v4[i] = drand48()*10.0-5;
      va1(i) = va1[i];
      va2(i) = va2[i];
      va3(i) = va3[i];
      va4(i) = va4[i];
#ifdef GLM
      vg1[i] = v1[i];
      vg2[i] =v2[i];
      vg3[i] =v3[i];
      vg4[i] = v4[i];
#endif
      for (j=0; j < NMAT; j++)
	{
	  //mml[i][j]= drand48()*10.0-5;
	  mm[i][j] = drand48()*10.0-5;//mml[i][j];
#ifdef GLM
	  mmg[i][j] = mm[i][j];
#endif
	  mma(i,j)=mm[i][j];
	}	  
    }
  //mm.show("mm=");
  srand48(0);
  if (argc==2)
    numtent=atoi(argv[1]);
  else
    numtent = 200;
  for (t=0; t < numtent; t++)
    {
#if 0
      x1 = sig*(ranf()-0.5);
      y1 = sig*(ranf()-0.5);
      x1c = complex<long double>(x1,0);
      x2c = complex<long double>(y1,0);
      x1 = sig*(ranf()-0.5);
      y1 = sig*(ranf()-0.5);
      x3c = complex<long double>(x1,y1);
      x4c = complex<long double>(x1,-y1);

      cc[4] = 1.0;
      cc[3] = -(x1c+x2c+x3c+x4c);
      cc[2] = (x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
      cc[1] = (-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
      cc[0] = (x1c*x2c*x3c*x4c);
      for (j=0; j < 5; j++)
	c[j] = cc[j].real();
      
      rp.set_coeff(c);
      //rp.forcehqr();
      rp.find_roots(r);
#endif 
  for (i=0; i < NMAT; i++)
    {
      for (j=0; j < NMAT; j++)
	{
	  //mml[i][j]= drand48()*10.0-5;
	  mm[i][j] = drand48()*1.0-0.5;//mml[i][j];
	  mm2[i][j] =drand48()*1.0-0.5;
	  mm3[i][j] =drand48()*1.0-0.5;
	  mm4[i][j] =drand48()*1.0-0.5;

#ifdef GLM
	  mmg[i][j] = mm[i][j];
	  mmg2[i][j] = mm2[i][j];
	  mmg3[i][j] = mm3[i][j];
	  mmg4[i][j] = mm4[i][j];
#endif
	  mma(i,j)=mm[i][j];
	  mma2(i,j)=mm2[i][j];
	  mma3(i,j)=mm3[i][j];
	  mma4(i,j)=mm4[i][j];
#ifdef EIGEN
	  mme(i,j)=mm[i][j];
	  mme2(i,j)=mm2[i][j];
	  mme3(i,j)=mm3[i][j];
	  mme4(i,j)=mm4[i][j];
#endif
	}	  
    }
#ifndef SEL
#define SEL 1
#endif
#if SEL==1
      //mm[0][0] += sin(t)/10000.0;
      // mm2 = mm;
      // mm2 *= mm;
      // mm2 *= mm;
      // mm2 *= mm.inv();
      //mm2 = (mm*mm2)*(mm*mm2);
      for (k=0; k < 100; k++)
	{
	  //v3=(v1+v2)+(v3+v4);
	  //mm += mm+mm2+mm3;
	  mm += mm2*mm3*mm4.inv()+(mm+mm4+mm3+mm4);
#if 1
	  //mm += mm+mm2+mm3;
	  //mm = (mm+mm2+mm3)*(mm2+mm3+mm)*(mm*mm);
	  //mm = (mm+mm2)+(mm3+mm4);
	  //mm.muladd(mm,mm2+mm3+mm4);
	  //mm = mm+mm4+mm2+mm3;
	  //mm.muladd(mm2, mm3);
	  //mm += mm2;
	  //mm += mm3;
	  //mm.fastadd(mm2,mm3);
#else
	  for (i=0; i < NMAT; i++)
	    for (j=0; j < NMAT; j++)
	      mm[i][j] += (mm[i][j] + mm3[i][j]+mm2[i][j]);//(mm[i][j]+mm4[i][j]);
#endif
	}
#elif SEL==2
      for (k=0; k < 100; k++)
	{
	  mmg += mmg2*mmg3*mmg4.inv()+(mmg+mm4g+mmg3+mmg4);
	  //mmg += mmg+mmg2+mmg3;
	  //mmg += mmg2*mmg3+inverse(mmg4);
	 //mmg =mmg2*mmg3*mmg4;
	  //mmg += inverse(mmg3)+mmg*mmg3;
	  //mmg += inverse(mmg4)+mmg*mmg4;
	  //mmg = (mmg+mmg4)+(mmg2+mmg3);
	}
#if 1
      //mmg2 += ((mmg+mmg)+(mmg+mmg))/mmg;
      //mmg2 += ((mmg*inverse(mmg))*(mmg*inverse(mmg)));
#else
      mmg2 = mmg;
      mmg2 *= mmg;
      mmg2 *= mmg;
      mmg2 *= inverse(mmg);
#endif
#elif SEL==3
      //mma(0,0) += sin(t)/(10000.0);
      //mma2 = ((mma*mma2)*(mma*mma2));
      for (k=0; k < 100; k++)
	{
	  //va3 = cross(cross(va1,va2),cross(va3,va4));
	  //va1 += mma2*va2+mma3*va3+mma4*va4;
	  //mma += mma*mma2*mma3.i()+mma4.i();
	  //mma = (mma3+mma2)+(mma+mma4);
	  //va3=(va1+va2)*3.0+(va3+va4)*5.0;
	  mma += mma2*mma3*mma4.i()+(mma+mma4+mma3+mma4);
	  //mma += mma2*mma3+mma4.i();  
	}
      //mma2 = ((mma*mma.i())*(mma*mma.i()));
#elif SEL==4
      for (k=0; k < 100; k++)
	{
	  //va3 = cross(cross(va1,va2),cross(va3,va4));
	  //va1 += mma2*va2+mma3*va3+mma4*va4;
	  //mma += mma*mma2*mma3.i()+mma4.i();
	  //mma = (mma3+mma2)+(mma+mma4);
	  //va3=(va1+va2)*3.0+(va3+va4)*5.0;
	  mme += mme2*mme3*mme4.inverse()+(mme+mme4+mme3+mme4);
	  //mme += mme2*mme3+mme4.inverse();  
	}

#endif
    }
#ifdef GLM
  printf("GLM last\n");
  for (int i=0; i < 3; i++)
    {
      for (int j=0; j < 3; j++) 
	printf("%.15G ", mmg2[i][j]);
      printf("\n");
    }
#endif
#if 0
#if SEL==1
  mm.show("mm=");
#elif SEL==3
#ifdef ARMA
  printf("lastarma=\n");
  mma.raw_print();
#endif
#endif
  mm=mm2+mm3;
  mm.show("mmfinal=");
  mma = mma2+mma3;
  printf("lastarmafinal=\n");
  mma.raw_print();
  //mm2 = ((mm*mm2)*(mm*mm2));
  //mm2.show("mm2last=");
#endif
  return 0;
}
