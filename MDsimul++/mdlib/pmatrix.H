#ifndef _PMATRIX_
/* N.B. le prestazione sono quelle di armadillo per N > 5 ma molto migliori per N < 5 */
#define _PMATRIX_
//#define USE_STRASSEN
// matrici NxM
#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include "pvector.H"
#define OPT_INV
#define USE_LAPACK
#ifdef USE_LAPACK
#undef __APPLE__ // undefining __APPLE__ one can prevent the usage of macosx accelerated framework switching to standard lapack blas
#ifdef __APPLE__
//  /System/Library/Frameworks/vecLib.framework/Headers/cblas.h
//  g++ -framework accelerate to compile
//#include <cblas.h>
//#include <lapacke.h>
#include <Accelerate/Accelerate.h>
#endif
#define MOVE_SEMANTIC
#define MAT_LAZY_EVAL
#ifndef __APPLE__
extern "C" {
extern void dgemm_(char* TRANSA, char* TRANSB, int* M, int *N, int*K, double* ALPHA,
	  const double *A, int *LDA, const double *B, int *LDB, double *BETA, double *C, int *LDC);      
extern void sgemm_(char* TRANSA, char* TRANSB, int* M, int *N, int*K, float* ALPHA,
	  const float *A, int *LDA, const float *B, int *LDB, float *BETA, float *C, int *LDC);      
extern int dgetrf_(int *__m, int *__n, double *__a, int *__lda, int *__ipiv, int *__info);
extern int sgetrf_(int *__m, int *__n, float *__a,  int *__lda, int *__ipiv, int *__info);
extern int dgetri_(int *__n, double *__a, int *__lda,
        int *__ipiv, double *__work, int*__lwork, int *__info);
extern int sgetri_(int *__n, float *__a, int *__lda,
        int *__ipiv, float *__work, int*__lwork, int *__info);
extern void dgemv_(char* TRANSA, int* M, int *N, double* ALPHA,
	  const double *A, int *LDA, const double *X, int *incx, double *BETA, double *Y, int *incy);  
extern void sgemv_(char* TRANSA, int* M, int *N, float* ALPHA,
	  const float *A, int *LDA, const float *X, int *incx, float *BETA, float *Y, int *incy);      
}
#endif
//Matrix times vector
void wrap_dgemv(char ta, double *Y, const double *A, const double *X, int n, double alpha=1.0, double beta=0.0)
{
  double ALPHA, BETA;
  int LDA, M, N, incx=1, incy=1;
  M=N=n;
  ALPHA=alpha;
  BETA=beta; 
  LDA=n;
#ifdef __APPLE__
  cblas_dgemv(CblasRowMajor, (ta=='t')?CblasNoTrans:CblasTrans, n, n, alpha, A, LDA, X, 1, beta, Y, 1);
#else
  dgemv_(&ta, &M, &N, &ALPHA, A, &LDA, X, &incx, &BETA, Y, &incy);
#endif
}

void wrap_sgemv(char ta, float *Y, const float *A, const float *X, int n, float alpha=1.0, float beta=0.0)
{
  float ALPHA, BETA;
  int LDA, M, N, incx=1, incy=1;
  M=N=n;
  ALPHA=alpha;
  BETA=beta; 
  LDA=n;
#ifdef __APPLE__
  cblas_sgemv(CblasRowMajor, (ta=='t')?CblasNoTrans:CblasTrans, n, n, alpha, A, LDA, X, 1, beta, Y, 1);
#else
  sgemv_(&ta, &M, &N, &ALPHA, A, &LDA, X, &incx, &BETA, Y, &incy);
#endif
}
//LU decomposition lapack routines
void wrap_dgetrf(double *A, int *IPIV, int n)
{
  //int *IPIV = new int[N+1];
  int N=n;
  int INFO;
  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
}
void wrap_sgetrf(float *A, int *IPIV, int n)
{
  int N=n;
  int INFO;
  sgetrf_(&N,&N,A,&N,IPIV,&INFO);
}
void wrap_dgetri(double *A, int *IPIV, int n)
{
  int N=n;
  int LWORK = n*n;
  int INFO;
  double *WORK = new double[LWORK];
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
  delete [] WORK;
}
void wrap_sgetri(float *A, int *IPIV, int n)
{
  int LWORK = n*n;
  int N=n;
  int INFO;
  float *WORK = new float[LWORK];
  sgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
  delete [] WORK;
}
void wrap_dgemm(double *C, const double *A, const double *B, const int n, double alpha=1.0, double beta=0.0)
{
  //C=A*B + beta*C if beta=0 C is not used as input
  double ALPHA, BETA;
  int LDA, LDB, LDC, M, N, K;
  char ta='n', tb='n';
  M=N=K=n;
  ALPHA=alpha;
  BETA=beta; 
  LDA=LDB=LDC=n;
#ifdef __APPLE__
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);      
#else
  // LAPACK and BLAS use column major convention hence we have to swap A and B
  dgemm_(&ta, &tb, &M, &N, &K, &ALPHA, B, &LDA, A, &LDB, &BETA, C, &LDC);      
#endif
}
void wrap_sgemm(float *C, const float *A, const float *B, int n, float alpha=1.0, float beta=0.0)
{
  //C=A*B + beta*C if beta=0 C is not used as input
  int LDA, LDB, LDC, M, N, K;
  char ta='n', tb='n';
  float ALPHA, BETA;
  M=N=K=n;
  ALPHA=alpha;
  BETA=beta; 
  LDA=LDB=LDC=n;
#ifdef __APPLE__
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);      
#else
  // LAPACK and BLAS use column major convention hence we have to swap A and B
  sgemm_(&ta, &tb, &M, &N, &K, &ALPHA, B, &LDA, A, &LDB, &BETA, C, &LDC);      
#endif
}
#endif
#ifdef MAT_LAZY_EVAL
//REMARK on Lazyness Mechanics:
//each class it is an operation of Lhs and Rhs, this class stores the values of the operands
//and only when the method get_m is called the operation is performed.
//Evaluation happens only as assign operator of pmatrixq class gets called (see below)
//
//Matrix operation types
struct MOpTypes
{  
  static const int MatPlusMat=0, MatMinusMat=1, MatTimesScal=2, ScalTimesMat=3;
};
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>
struct AnOpM 
{
  Lhs const& lhs;
  Rhs const& rhs;
  AnOpM(Lhs const& lhs, Rhs const& rhs):lhs(lhs), rhs(rhs) {
    // empty body
  }

  inline ntype get_m(int i, int j) const {
    if constexpr (tipo == MOpTypes::MatPlusMat)// mat + mat
      {
	return lhs.get_m(i,j) + rhs.get_m(i,j);
      }
    else if constexpr (tipo == MOpTypes::MatMinusMat) // mat - mat
      {
	return lhs.get_m(i,j) - rhs.get_m(i,j);
      }
    else if constexpr (tipo == MOpTypes::MatTimesScal) // mat * scalar [NOT IMPLEMENTED YET]
      {
	return lhs.get_m(i,j)*rhs;
      }
    else if constexpr (tipo == MOpTypes::ScalTimesMat) // scalar * mat [NOT IMPLEMENTED YET]
      {
	return lhs*rhs.get_m(i,j);
      }
    return 0;
  }
  Lhs const& get_lhs() const { return lhs; }//unused
  Rhs const& get_rhs() const { return rhs; }
};
#endif
template <class ntype, int N, int M> class pmatrix {
  //int nr, nc;
public:
  ntype m[N][M];
  constexpr static int dynamic = false;
#if 0
  pmatrix<ntype,N,M>& operator=(const pmatrix<ntype,N,M>& m1)
    {
      int i,j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] = m1.m[i][j];
      return (*this);
    }
#endif
};
template <class ntype, int N, int M> class pmatrixd {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  ntype (*m)[M]; // N riga M colonna
  // copy assignment
 pmatrixd<ntype,N,M>& operator=(const pmatrixd<ntype,N,M>& m1)
    {
      int i,j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] = m1.m[i][j];
      return (*this);
    }
#ifdef MOVE_SEMANTIC
  // move assignment
  pmatrixd<ntype,N,M>& operator=(pmatrixd<ntype,N,M>&& m1)
    {
      swap(m, m1.m);
      return (*this);
    }
#endif
  // copy constructor 
  pmatrixd(const pmatrixd<ntype,N,M>& m1)
    {
      //nr=N;
      //nc=M;
      m = new ntype[N][M];
      (*this) = m1;
    }
#ifdef MOVE_SEMANTIC
  // move constructor
  pmatrixd(pmatrixd<ntype,N,M>&& m1)
    {
      //nr=N;
      //nc=M;
      //cout << "qui\n";
      (*this).m = m1.m;
      m1.m=nullptr;
    }
#endif
  // default constructor
  pmatrixd()
    {
      //nr=N;
      //nc=M;
      //cout << " //nel costruttore//\n ";
      m = new ntype[N][M];
    }
  // destructor
  ~pmatrixd()
    {
      delete[] m;
    }
};	
// redefine also *= += and
// define inv() det() i.e. inverse and determinant
// matrici quadrate NxN
using namespace std;
#if 0
namespace matrixns
{
  template <class ntype> ntype abs(ntype a)
    {
      return abs(a);
    }
}
using namespace matrixns;
#endif
#if 0
template<int N>
typename enable_if<N == 4>::type
inverse4(void) {
    return ;
}
#endif
//NMAX_=45, int MAXSTRA_=40, int NMAXINV_=40, int NMAXMUL_=10, int NSTA_=8, int NLAZY_=5
struct matpars
{
  // params for my class
  static const int NMAX = 45; // if N>NMAX use dynamic allocation of matrices
  static const int NMAXSTRA = 40;
  static const int NMAXINV = 40; // if N>NAMXINV use LAPACK LU DECOMP
  static const int NMAXMUL = 10; // if N>NMAXMUL use BLAS multiplication routines
  static const int NSTA = 8;// value below which declares static object in operator funcs (only if m is not dynamically allocated)
  static const int NLAZY = 5;// value below which does not use lazy evaluation for addition
};
// identity matrix recursively defined
//


// for large matrices use dynamic allocation to avoid stack overflow!
template <class ntype, int N> using pmatrixb = typename std::conditional<(N>matpars::NMAX),pmatrixd <ntype, N, N>,pmatrix <ntype, N, N>>::type;

template <class ntype, int N> class pmatrixq: public pmatrixb<ntype,N>, public matpars
  //public std::conditional<(N>100),pmatrixd <ntype, N, N>,pmatrix <ntype, N, N>>::type 
{
  using pmatrixb<ntype,N>::m;
  //int n, maxdigits;
  // il const vuol dire che il puntore this, passato ai membri di una classe implicitamente, viene passato come const
  //using type = pmatrixq;
  inline pmatrixq<ntype,N> ludcmp(pvector<int,N>& indx, double& d) const
    {
      /* return the LU decomposition of calling matrix (i.e. *this) */
      pmatrixq <ntype,N> m1;
      const double TINY=1E-40;
      int i,imax=-1,j,k,n;
      ntype big,dum,sum,temp; 
      ntype vv[N]; /* vv stores the implicit scaling of each row.*/
      d=1.0; /* No row interchanges yet. */
      n=N;
      m1 = (*this);
      for (i=0;i<n;i++) 
	{ 
	  /* Loop over rows to get the implicit scaling information.*/ 
	  big=0.0; 
	  for (j=0;j<n;j++)
	    {
	      if ((temp=abs(m1[i][j])) > big) big=temp; 
	    }
	  if (big == 0.0)
	    {
	      return m1;
	    }
	  /* No nonzero largest element. */
	  vv[i]=1.0/big; /* Save the scaling.*/
	} 
      for (j=0;j<n;j++) 
	{ /* This is the loop over columns of Crout s method.*/
	  for (i=0;i<j;i++) 
	    { 
	      /* This is equation (2.3.12) except for i = j. */
	      sum=m1[i][j]; 
	      for (k=0;k<i;k++) 
		sum -= m1[i][k]*m1[k][j]; 
	      m1[i][j]=sum; 
	    } 
	  big=0.0; /* Initialize for the search for largest pivot element. */ 
	  for (i=j;i<n;i++) 
	    { 
	      /* This is i = j of equation (2.3.12) and i = j+1. . .N of equation (2.3.13).*/
	      sum=m1[i][j]; 
	      for (k=0;k<j;k++)
		sum -= m1[i][k]*m1[k][j]; 
	      m1[i][j]=sum; 
	      if ( (dum=vv[i]*abs(sum)) >= big) 
		{ 
		  /* Is the  gure of merit for the pivot better than the best so far? */
		  big=dum; 
		  imax=i; 
		} 
	    } 
	  if (j != imax) 
	    { 
	      /* Do we need to interchange rows? */
	      for (k=0;k<n;k++) 
		{ 
		  /* Yes, do so...*/ 
		  dum=m1[imax][k]; 
		  m1[imax][k]=m1[j][k]; 
		  m1[j][k]=dum; 
		} 
	      d = -d; 
	      /* ...and change the parity of d. */ 
	      vv[imax]=vv[j]; 
	      /* Also interchange the scale factor.*/ 
	    } 
	  indx[j]=imax; 
	  if (m1[j][j] == 0.0) 
	    m1[j][j]=TINY; 
	  /* If the pivot element is zero the matrix is singular 
	   * (at least to the precision of the algorithm). 
	   * For some applications on singular matrices, 
	   * it is desirable to substitute TINY for zero. */ 
	  if (j != n) 
	    { 
	      /* Now,  nally, divide by the pivot element.*/
	      dum=1.0/(m1[j][j]); 
	      for (i=j+1;i<n;i++) m1[i][j] *= dum; 
	    } 
	} 
      return m1;
    }

  inline pvector<ntype,N> lubksb(pvector <int,N> indx, pvector<ntype,N> b) const
    { 
      /* the calling matrix (*this) must be the LU decomposition of a given matrix M
       * with indx the raw permutations done to obtain LU from M, 
       * the return vector is the solution of A*x + b = 0 */
      int i,ii=0,ip,j; 
      double sum; 
      for (i=0;i<N;i++) 
	{ 
	  /* When ii is set to a positive value, it will become the index of the  
	   * rst nonvanishing element of b. Wenow do the forward substitution,
	   * equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go. */
	  ip=indx[i];
	  sum=b[ip];
	  b[ip]=b[i]; 
	  if (ii>-1) 
	    for (j=ii;j<=i-1;j++) 
	      sum -= m[i][j]*b[j]; 
	  else if (sum) 
	    ii=i; 
	  /* A nonzero element was encountered, so from now on we will have to do 
	   * the sums in the loop above. */ 
	  b[i]=sum; 
	} 
      for (i=N-1;i>=0;i--) 
	{ 
	  /* Now we do the backsubstitution, equation (2.3.7).*/
	  sum=b[i]; 
	  for (j=i+1;j<N;j++) 
	    sum -= m[i][j]*b[j]; b[i]=sum/m[i][i]; 
	  /* Store a component of the solution vector X. */ 
	} /* All done! */
      return b;
    }
public:
#if 1
  // TODO: make this static
  inline pmatrixq<ntype,N> I() const
    {
      pmatrixq<ntype,N> MI;
      for (auto i=0; i < N; i++)
	for (auto j=0; j < N; j++)
	  MI[i][j] = (i==j)?1.0:0.0;
     return MI;
    }
#endif
#ifdef MAT_LAZY_EVAL
  template<typename Lhs, typename Rhs, int tipo>
    pmatrixq<ntype,N>& operator=(AnOpM<ntype,N, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	for (int j=0; j < N; j++)
	  m[i][j] = op.get_m(i,j);
      return (*this);
    }
  template<typename Lhs, typename Rhs, int tipo>
    pmatrixq<ntype,N>(AnOpM<ntype,N, tipo, Lhs, Rhs> const& op) {
      //cout << "qui\n";
       for (int i=0; i < N; i++)
	for (int j=0; j < N; j++)
	  m[i][j] = op.get_m(i,j);
    }
  pmatrixq<ntype,N>() = default;
  inline ntype get_m(int i, int j) const {return m[i][j];}
  inline void* get_mptr() const {return m;}
#endif
#ifndef MAT_LAZY_EVAL
  inline pmatrixq<ntype,N> operator+(const pmatrixq<ntype,N>& m1) const
    {
      //cout << "sono qui<<<<<<<<<<<<\n";
      return addition(m1);
    }
  inline pmatrixq<ntype,N> operator-(const pmatrixq<ntype,N>& m1) const
    {
      return subtraction(m1);
    }
#endif
  // addition
  inline pmatrixq<ntype,N> addition(const pmatrixq<ntype,N>& m1) const
    {
      if constexpr (N <= NSTA && pmatrixq<ntype,N>::dynamic==false)
	{
	  static pmatrixq<ntype, N> m2;
	  int i, j;
	  for (i=0; i < N; i++)
	    for (j=0; j < N; j++)
	      {
#if 1
      		m2.m[i][j] = m[i][j] + m1.m[i][j];
#else
    		// NOTA: const_cast rimuove l'attributo const ad m1 e rende possibile l'applicazione dell'operatore [] 
    		// che restituisce ntype*, che non è const (altimenti si ottiene un errore nella compilazione).
    		(const_cast<pmatrixq<ntype,N>&>(m2))[i][j] = m[i][j] + (const_cast<pmatrixq<ntype,N>&>(m1))[i][j];
#endif
    	      }
	  return m2;
	}
      else
	{
	  pmatrixq<ntype, N> m2;
	  int i, j;
	  for (i=0; i < N; i++)
	    for (j=0; j < N; j++)
    	      {
      		m2.m[i][j] = m[i][j] + m1.m[i][j];
    	      }
	  if constexpr (N > NMAX)
	    return std::move(m2);
	  else
	    return m2;
	}
    }
  inline pmatrixq<ntype,N> subtraction(const pmatrixq<ntype,N>& m1) const
    {
      int i, j;
      if constexpr (N <= NSTA && pmatrixq<ntype,N>::dynamic==false)
	{
	  static pmatrixq<ntype, N> m2;
	  for (i=0; i < N; i++)
	    for (j=0; j < N; j++)
	      {
#if 1
      		m2.m[i][j] = m[i][j] - m1.m[i][j];
#else
    		// NOTA: const_cast rimuove l'attributo const ad m1 e rende possibile l'applicazione dell'operatore [] 
    		// che restituisce ntype*, che non è const (altimenti si ottiene un errore nella compilazione).
    		(const_cast<pmatrixq<ntype,N>&>(m2))[i][j] = m[i][j] - (const_cast<pmatrixq<ntype,N>&>(m1))[i][j];
#endif
	  }
	  return m2;
	}
      else
	{
	  pmatrixq<ntype, N> m2;
	  for (i=0; i < N; i++)
	    for (j=0; j < N; j++)
	      {
#if 1
      		m2.m[i][j] = m[i][j] - m1.m[i][j];
#else
    		// NOTA: const_cast rimuove l'attributo const ad m1 e rende possibile l'applicazione dell'operatore [] 
    		// che restituisce ntype*, che non è const (altimenti si ottiene un errore nella compilazione).
    		(const_cast<pmatrixq<ntype,N>&>(m2))[i][j] = m[i][j] - (const_cast<pmatrixq<ntype,N>&>(m1))[i][j];
#endif
	  }
	  if constexpr (N > NMAX)
	    return std::move(m2);
	  else 
	    return m2;
	}
    }

  inline double detgenNR() const
    {
      double d;
      static pmatrixq<ntype,N> m1;
      int j;
      static pvector<int,N> indx;
      // ludcmp return the LU decomposition of actual matrix (*this)
      m1=ludcmp(indx,d);
      for(j=0;j<N;j++) 
	d *= m1[j][j];
      return d;
    }
  inline ntype detgen() const
    {
#ifdef USE_LAPACK
      if constexpr (is_same<ntype, double>::value)
	{
	  double d=1.0;
	  int j;	
	  static pmatrixq<ntype,N> om;
	  int *IPIV = new int[N];
	  om = *this;
	  wrap_dgetrf(&(om.m[0][0]), IPIV, N);
	  for(j=0;j<N;j++)
	    { 
	      d = d*om[j][j];
	      if(IPIV[j]!=j+1) 
		d*=-1.0;
	    }
	  delete [] IPIV;
	  return d;
	}
      else if constexpr (is_same<ntype, float>::value)
	{
	  double d=1.0;
	  int j;	
	  static pmatrixq<ntype,N> om;
	  int *IPIV = new int[N];
	  om = *this;
	  wrap_fgetrf(&(om.m[0][0]), IPIV, N);
	  for(j=0;j<N;j++)
	    { 
	      d = d*om[j][j];
	      if(IPIV[j]!=j+1) 
		d*=-1.0;
	    }
	  delete [] IPIV;
	  return d;
	}
      else
	{
	  return detgenNR();
	}
#else
      return detgenNR();
#endif 
    }
  inline void invLUNR(pmatrixq<ntype,N>& om) const
    {
      static pvector <ntype,N> col;
      static pmatrixq<ntype,N> LU;
      int m2, m1;
      double d;
      static pvector<int,N> indx;
      LU=ludcmp(indx, d); 
      for(m2=0;m2<N;m2++) 
	{ 
	  for(m1=0;m1<N;m1++) 
	    col[m1]=0.0; 
	  col[m2]=1.0; 
	  col=LU.lubksb(indx,col);
	  for(m1=0;m1<N;m1++) 
	    om[m1][m2]=col[m1]; 
	}
    }
  inline ntype det() const
    {
      if constexpr (N==1)
	{
	  return m[0][0];
	}
      else if constexpr (N==2)
	{
	  return m[0][0]*m[1][1]-m[0][1]*m[1][0];
	}
      else if constexpr (N==3)
	{
	  return -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + 
	    m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - 
	    m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2];
	} 
      else if constexpr (N == 4)
	{	  
	  return m[0][1]*m[1][3]*m[2][2]*m[3][0] - m[0][1]*m[1][2]*m[2][3]*m[3][0] - 
	    m[0][0]*m[1][3]*m[2][2]*m[3][1] + m[0][0]*m[1][2]*m[2][3]*m[3][1] - 
	    m[0][1]*m[1][3]*m[2][0]*m[3][2] + m[0][0]*m[1][3]*m[2][1]*m[3][2] + 
	    m[0][1]*m[1][0]*m[2][3]*m[3][2] - m[0][0]*m[1][1]*m[2][3]*m[3][2] + 
	    m[0][3]*(m[1][2]*m[2][1]*m[3][0] - m[1][1]*m[2][2]*m[3][0] - m[1][2]*m[2][0]*m[3][1] + 
		     m[1][0]*m[2][2]*m[3][1] + m[1][1]*m[2][0]*m[3][2] - m[1][0]*m[2][1]*m[3][2]) + 
	    (m[0][1]*m[1][2]*m[2][0] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] + 
	     m[0][0]*m[1][1]*m[2][2])*m[3][3] + 
	    m[0][2]*(-(m[1][3]*m[2][1]*m[3][0]) + m[1][1]*m[2][3]*m[3][0] + m[1][3]*m[2][0]*m[3][1] - 
		     m[1][0]*m[2][3]*m[3][1] - m[1][1]*m[2][0]*m[3][3] + m[1][0]*m[2][1]*m[3][3]);
	}
      else 
	return detgen();
    }

  inline pmatrixq<ntype,N> invgen() const
    {
      pmatrixq<ntype,N> om;
      // ludcmp returns the LU decomposition of actual matrix (*this)
      // where indx will store the raw permutiation effected by partial pivoting (see Numerical Recipe
      // for more details c2.3)
#ifdef USE_LAPACK
      om = (*this);
      if constexpr (N > NMAXINV)
	{
	  if constexpr (is_same<ntype, double>::value)
	    {    
	      int *IPIV = new int[N];
	      wrap_dgetrf(&(om.m[0][0]), IPIV, N);
	      wrap_dgetri(&(om.m[0][0]), IPIV, N);
	      delete[] IPIV;
	    }
	  else if constexpr (is_same<ntype, float>::value)
	    {
	      int *IPIV = new int[N];
	      wrap_fgetrf(&(om.m[0][0]), IPIV, N);
	      wrap_fgetri(&(om.m[0][0]), IPIV, N);
	      delete[] IPIV;
	    }
	  else 
	    {
	      invLUNR(om);
	    }
	}
      else
	invLUNR(om);
#else
      invLUNR(om);
#endif
      return om;
    }
  pmatrixq<ntype,N> inv() const
    {
      double invd;
      pmatrixq <ntype,N> m1;
      if constexpr (N==1)
	{
	  m1[0][0]=1.0/m[0][0];
	}
      else if constexpr (N==2)
	{
	  invd = 1.0/det();
	  m1[0][0] = m[1][1]*invd;
	  m1[0][1] = m[1][0]*invd;
	  m1[1][0] = m[0][1]*invd;
	  m1[1][1] = m[0][0]*invd;
	}
      else if constexpr (N==3)
	{
#ifdef OPT_INV
	  //return -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + 
	    //m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - 
	   // m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2];

	  ntype m11m20= m[1][1]*m[2][0];
	  ntype m12m20= m[1][2]*m[2][0];
	  ntype m10m21= m[1][0]*m[2][1];
	  ntype m12m21= m[1][2]*m[2][1];
	  ntype m10m22= m[1][0]*m[2][2];
	  ntype m11m22= m[1][1]*m[2][2];
	  invd = 1.0/(-m[0][2]*m11m20 + m[0][1]*m12m20 + m[0][2]*m10m21
		      - m[0][0]*m12m21 - m[0][1]*m10m22 + m[0][0]*m11m22);

	  m1[0][0] = (-m12m21 + m11m22)*invd;
	  m1[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2])*invd;
	  m1[0][2] = (-m[0][2]*m[1][1] + m[0][1]*m[1][2])*invd;
	  m1[1][0] = (m12m20 - m10m22)*invd;
	  m1[1][1] = (-m[0][2]*m[2][0] + m[0][0]*m[2][2])*invd; 
	  m1[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2])*invd; 
	  m1[2][0] = (-m11m20 + m10m21)*invd; 
	  m1[2][1] = (m[0][1]*m[2][0] - m[0][0]*m[2][1])*invd;
	  m1[2][2] = (-m[0][1]*m[1][0] + m[0][0]*m[1][1])*invd;
#else
	  invd = 1.0/det();
	  m1[0][0] = -m[1][2]*m[2][1] + m[1][1]*m[2][2];
	  m1[0][1] =  m[0][2]*m[2][1] - m[0][1]*m[2][2];
	  m1[0][2] = -m[0][2]*m[1][1] + m[0][1]*m[1][2];
	  m1[1][0] =  m[1][2]*m[2][0] - m[1][0]*m[2][2];
	  m1[1][1] = -m[0][2]*m[2][0] + m[0][0]*m[2][2]; 
	  m1[1][2] =  m[0][2]*m[1][0] - m[0][0]*m[1][2]; 
	  m1[2][0] = -m[1][1]*m[2][0] + m[1][0]*m[2][1]; 
	  m1[2][1] =  m[0][1]*m[2][0] - m[0][0]*m[2][1];
	  m1[2][2] = -m[0][1]*m[1][0] + m[0][0]*m[1][1];
	  m1[0][0] *= invd;
	  m1[0][1] *= invd;
	  m1[0][2] *= invd;
	  m1[1][0] *= invd; 
	  m1[1][1] *= invd;
	  m1[1][2] *= invd;
	  m1[2][0] *= invd;
	  m1[2][1] *= invd; 
	  m1[2][2] *= invd; 	    
#endif
	}
      else if constexpr (N==4)
	{
#ifdef OPT_INV
	  ntype s0 = m[0][0] * m[1][1] - m[1][0] * m[0][1];
      	  ntype s1 = m[0][0] * m[1][2] - m[1][0] * m[0][2];
	  ntype s2 = m[0][0] * m[1][3] - m[1][0] * m[0][3];
	  ntype s3 = m[0][1] * m[1][2] - m[1][1] * m[0][2];
	  ntype s4 = m[0][1] * m[1][3] - m[1][1] * m[0][3];
	  ntype s5 = m[0][2] * m[1][3] - m[1][2] * m[0][3];

	  ntype c5 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
	  ntype c4 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
	  ntype c3 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
	  ntype c2 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
	  ntype c1 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
	  ntype c0 = m[2][0] * m[3][1] - m[3][0] * m[2][1];

	  // Should check for 0 determinant
	  invd = 1.0 / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);


	  m1[0][0] = ( m[1][1] * c5 - m[1][2] * c4 + m[1][3] * c3) * invd;
	  m1[0][1] = (-m[0][1] * c5 + m[0][2] * c4 - m[0][3] * c3) * invd;
	  m1[0][2] = ( m[3][1] * s5 - m[3][2] * s4 + m[3][3] * s3) * invd;
	  m1[0][3] = (-m[2][1] * s5 + m[2][2] * s4 - m[2][3] * s3) * invd;

	  m1[1][0] = (-m[1][0] * c5 + m[1][2] * c2 - m[1][3] * c1) * invd;
	  m1[1][1] = ( m[0][0] * c5 - m[0][2] * c2 + m[0][3] * c1) * invd;
	  m1[1][2] = (-m[3][0] * s5 + m[3][2] * s2 - m[3][3] * s1) * invd;
	  m1[1][3] = ( m[2][0] * s5 - m[2][2] * s2 + m[2][3] * s1) * invd;

	  m1[2][0] = ( m[1][0] * c4 - m[1][1] * c2 + m[1][3] * c0) * invd;
	  m1[2][1] = (-m[0][0] * c4 + m[0][1] * c2 - m[0][3] * c0) * invd;
	  m1[2][2] = ( m[3][0] * s4 - m[3][1] * s2 + m[3][3] * s0) * invd;
      	  m1[2][3] = (-m[2][0] * s4 + m[2][1] * s2 - m[2][3] * s0) * invd;

	  m1[3][0] = (-m[1][0] * c3 + m[1][1] * c1 - m[1][2] * c0) * invd;
	  m1[3][1] = ( m[0][0] * c3 - m[0][1] * c1 + m[0][2] * c0) * invd;
	  m1[3][2] = (-m[3][0] * s3 + m[3][1] * s1 - m[3][2] * s0) * invd;
	  m1[3][3] = ( m[2][0] * s3 - m[2][1] * s1 + m[2][2] * s0) * invd;

#else
	  invd = 1.0/det();
	  
	  //printf("invD=%.15G\n", invd);
	  m1[0][0] = -(m[1][3]*m[2][2]*m[3][1]) + m[1][2]*m[2][3]*m[3][1] + m[1][3]*m[2][1]*m[3][2] - 
	    m[1][1]*m[2][3]*m[3][2] - m[1][2]*m[2][1]*m[3][3] + m[1][1]*m[2][2]*m[3][3];
	  m1[0][1] =m[0][3]*m[2][2]*m[3][1] - m[0][2]*m[2][3]*m[3][1] - m[0][3]*m[2][1]*m[3][2] + 
	    m[0][1]*m[2][3]*m[3][2] + m[0][2]*m[2][1]*m[3][3] - m[0][1]*m[2][2]*m[3][3];
	  m1[0][2] =-(m[0][3]*m[1][2]*m[3][1]) + m[0][2]*m[1][3]*m[3][1] + m[0][3]*m[1][1]*m[3][2] - 
	    m[0][1]*m[1][3]*m[3][2] - m[0][2]*m[1][1]*m[3][3] + m[0][1]*m[1][2]*m[3][3];
	  m1[0][3] = m[0][3]*m[1][2]*m[2][1] - m[0][2]*m[1][3]*m[2][1] - m[0][3]*m[1][1]*m[2][2] +
	    m[0][1]*m[1][3]*m[2][2] + m[0][2]*m[1][1]*m[2][3] - m[0][1]*m[1][2]*m[2][3];
	  m1[1][0] = m[1][3]*m[2][2]*m[3][0] - m[1][2]*m[2][3]*m[3][0] - m[1][3]*m[2][0]*m[3][2] + 
	    m[1][0]*m[2][3]*m[3][2] + m[1][2]*m[2][0]*m[3][3] - m[1][0]*m[2][2]*m[3][3];
	  m1[1][1] = -(m[0][3]*m[2][2]*m[3][0]) + m[0][2]*m[2][3]*m[3][0] + m[0][3]*m[2][0]*m[3][2] - 
	    m[0][0]*m[2][3]*m[3][2] - m[0][2]*m[2][0]*m[3][3] + m[0][0]*m[2][2]*m[3][3];
	  m1[1][2] = m[0][3]*m[1][2]*m[3][0] - m[0][2]*m[1][3]*m[3][0] - m[0][3]*m[1][0]*m[3][2] + 
	    m[0][0]*m[1][3]*m[3][2] + m[0][2]*m[1][0]*m[3][3] - m[0][0]*m[1][2]*m[3][3] ;
	  m1[1][3] = -(m[0][3]*m[1][2]*m[2][0]) + m[0][2]*m[1][3]*m[2][0] + m[0][3]*m[1][0]*m[2][2] - 
	    m[0][0]*m[1][3]*m[2][2] - m[0][2]*m[1][0]*m[2][3] + m[0][0]*m[1][2]*m[2][3];
	  m1[2][0] = -(m[1][3]*m[2][1]*m[3][0]) + m[1][1]*m[2][3]*m[3][0] + m[1][3]*m[2][0]*m[3][1] - 
	    m[1][0]*m[2][3]*m[3][1] - m[1][1]*m[2][0]*m[3][3] + m[1][0]*m[2][1]*m[3][3];
	  m1[2][1] = m[0][3]*m[2][1]*m[3][0] - m[0][1]*m[2][3]*m[3][0] - m[0][3]*m[2][0]*m[3][1] + 
	    m[0][0]*m[2][3]*m[3][1] + m[0][1]*m[2][0]*m[3][3] - m[0][0]*m[2][1]*m[3][3];
	  m1[2][2] = -(m[0][3]*m[1][1]*m[3][0]) + m[0][1]*m[1][3]*m[3][0] + m[0][3]*m[1][0]*m[3][1] - 
	    m[0][0]*m[1][3]*m[3][1] - m[0][1]*m[1][0]*m[3][3] + m[0][0]*m[1][1]*m[3][3];
	  m1[2][3] =m[0][3]*m[1][1]*m[2][0] - m[0][1]*m[1][3]*m[2][0] - m[0][3]*m[1][0]*m[2][1] + 
	    m[0][0]*m[1][3]*m[2][1] + m[0][1]*m[1][0]*m[2][3] - m[0][0]*m[1][1]*m[2][3];
	  m1[3][0] = m[1][2]*m[2][1]*m[3][0] - m[1][1]*m[2][2]*m[3][0] - m[1][2]*m[2][0]*m[3][1] + 
	    m[1][0]*m[2][2]*m[3][1] + m[1][1]*m[2][0]*m[3][2] - m[1][0]*m[2][1]*m[3][2];
	  m1[3][1] =-(m[0][2]*m[2][1]*m[3][0]) + m[0][1]*m[2][2]*m[3][0] + m[0][2]*m[2][0]*m[3][1] - 
	    m[0][0]*m[2][2]*m[3][1] - m[0][1]*m[2][0]*m[3][2] + m[0][0]*m[2][1]*m[3][2];
	  m1[3][2] = m[0][2]*m[1][1]*m[3][0] - m[0][1]*m[1][2]*m[3][0] - m[0][2]*m[1][0]*m[3][1] + 
	    m[0][0]*m[1][2]*m[3][1] + m[0][1]*m[1][0]*m[3][2] - m[0][0]*m[1][1]*m[3][2];
	  m1[3][3] = -(m[0][2]*m[1][1]*m[2][0]) + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - 
	    m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2];
	  m1[0][0] *= invd;
	  m1[0][1] *= invd;
	  m1[0][2] *= invd;
	  m1[0][3] *= invd;
	  m1[1][0] *= invd; 
	  m1[1][1] *= invd;
	  m1[1][2] *= invd;
	  m1[1][3] *= invd;
	  m1[2][0] *= invd;
	  m1[2][1] *= invd; 
	  m1[2][2] *= invd;
	  m1[2][3] *= invd;
	  m1[3][0] *= invd;
	  m1[3][1] *= invd; 
	  m1[3][2] *= invd;
	  m1[3][3] *= invd;
#endif
	}
      else
	m1=invgen(); 

      if constexpr (N > NMAX)
	return std::move(m1);
      else
	return m1;
    }
  inline pmatrixq<ntype, N>& operator += (const pmatrixq<ntype, N>& m1)
    {
      int i, j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] += m1.m[i][j];
      return (*this);
    }
  
  inline pmatrixq<ntype, N>& operator -= (const pmatrixq<ntype, N>& m1)
    {
      int i, j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] -= m1.m[i][j];
      return (*this);
    }
  // multiplication with scalars
  inline pmatrixq<ntype, N>& operator *=(const double& param) 
    {
      int i,j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] *= param;
      return (*this);
    }
  inline pmatrixq<ntype, N> operator *(const double& param) 
    {
      int i,j;
      pmatrixq<ntype,N> m1;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m1[i][j] = m[i][j]*param;
      return m1;
    }
  inline pmatrixq<ntype, N>& operator /=(const double& param) 
    {
      int i,j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] /= param;
      return (*this);
    }
  inline pmatrixq<ntype, N> operator /(const double& param) 
    {
      int i,j;
      pmatrixq<ntype,N> m1;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m1[i][j] = m[i][j]/param;
      return m1;
    }
  friend pmatrixq<ntype,N> operator *(const double& param, const pmatrixq<ntype,N>& m1) 
    {
      int i, j;
      pmatrixq<ntype,N> m2;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m2.m[i][j] = m1.m[i][j]*param;
      return m2;
    }
  // multiply two NxN matrices
  inline pmatrixq<ntype,N>& operator *=(const pmatrixq<ntype,N>& m1) 
    {
      *this = (*this).mul(m1);
      return *this;
    }
  inline pmatrixq<ntype,N> operator *(const pmatrixq<ntype,N>& m1) 
    {
      if constexpr (N > NMAX)
   	return std::move((*this).mul(m1));// std:move is needed otherwise move constructor will not be used!!
      else
   	return (*this).mul(m1);// std:move is needed otherwise move constructor will not be used!!
    }
#if 0
  inline pmatrixq<ntype,N> operator *(const pmatrixq<ntype,N>&& m1) 
    {
      //cout << "qui\n";
      return (*this).mul(m1);	
    }
#endif
   inline pmatrixq<ntype,N>& operator /=(const pmatrixq<ntype,N>& m1) 
    {
      *this = (*this).mul(m1.inv());	
      return *this;
    }
  
  inline pmatrixq<ntype,N> operator /(const pmatrixq<ntype,N>& m1) 
    {
      if constexpr (N > NMAX)
   	return std::move((*this).mul(m1.inv()));// std:move is needed otherwise move constructor will not be used!!	
      else 
	return (*this).mul(m1.inv());
    }
  
  // matrix times vector
  inline pvector<ntype,N> operator *(const pvector<ntype,N>& v1) 
    {
      int i, j;
      pvector <ntype,N> v2;
#ifdef USE_LAPACK
      if constexpr (N > NMAXMUL) 
	{
	  if constexpr (is_same<ntype, double>::value)
	    wrap_dgemv('t', &(v2.v[0]), &(m[0][0]), &(v1.v[0]), N);
	  else if constexpr (is_same<ntype, float>::value)
	    wrap_sgemv('t', &(v2.v[0]), &(m[0][0]), &(v1.v[0]), N);
	  else
	    {
	      for (i=0; i < N; i++)
		{
		  v2.v[i]=0;
		  for (j=0; j < N; j++)
		    v2.v[i] += m[i][j]*v1.v[j];
		}
	    }
	}
      else
	{
    	  for (i=0; i < N; i++)
	    {
    	      v2.v[i]=0;
	      for (j=0; j < N; j++)
		v2.v[i] += m[i][j]*v1.v[j];
	    }
	}
#else	
      for (i=0; i < N; i++)
	{
	  v2.v[i]=0;
	  for (j=0; j < N; j++)
	    v2.v[i] += m[i][j]*v1.v[j];
	}
#endif
      return v2;
    }
  // transpose of v1 times m1
  friend pvector<ntype,N> operator *(const pvector<ntype,N> v1, const pmatrixq<ntype,N>& m1) 
    {
      int i, j;
      pvector <ntype,N> v2;
#ifdef USE_LAPACK
      if constexpr (N > NMAXMUL) 
	{
	  if constexpr (is_same<ntype, double>::value)
	    wrap_dgemv('n', &(v2.v[0]), &(m[0][0]), &(v1.v[0]), N);
	  else if constexpr (is_same<ntype, float>::value)
	    wrap_sgemv('n', &(v2.m[0][0]), &(m[0][0]), &(v1.v[0]), N);
	  else
	    {
    	      v2.v[i]=0;
    	      for (j=0; j < N; j++)
    		v2.v[i] += v1.v[j]*m1.m[j][i];
	    }
	}
      else
	{
    	  v2.v[i]=0;
	  for (j=0; j < N; j++)
	    v2.v[i] += v1.v[j]*m1.m[j][i];
	}
#else	
      for (i=0; i < N; i++)
	{
	  v2.v[i]=0;
	  for (j=0; j < N; j++)
	    v2.v[i] += v1.v[j]*m1.m[j][i];
	}
#endif
      return v2;
    }
  inline ntype* operator[](const int& i)
    {
      return &(m[i][0]);
    }

  inline  void mulblock(pmatrixq<ntype,N>& m2, const pmatrixq<ntype,N>& m1)
    {
      /* tiled matrix multiplication to optimize cache usage (from wikipedia!) 
      */
      int i, j, k;
      int T=(int)sqrt(N);
      int I, J, K;
      ntype sum;
      for (I=0; I < N; I+=T)
	{
	  for (J=0; J < N; J+=T)
	    {
	      for (K=0; K < N; K+=T)
		{
		  for (i=I; i < min(I+T,N); i++)
		    {
		      for (j=J; j < min(J+T,N); j++)
			{
			  sum = 0.0;
			  for (k=K; k < min(K+T,N); k++)
			    sum += pmatrixq::m[i][k]*m1.m[k][j];
			  if (K==0)
			    m2.m[i][j]=0;
			  m2.m[i][j] += sum;
			}
		    }
		}
	    }
	}
    }
  inline void mulstra(pmatrixq<ntype,N>& m2, const pmatrixq<ntype,N>& m1)
    {
      int i, j;
      using selftype  = typename std::conditional<N%2==1,pmatrixq<ntype,(N+1)/2>,pmatrixq<ntype,N/2>>::type; 
      static constexpr int dN = N % 2;
      static constexpr int Np = N+dN;
      selftype A11, A12, A21, A22, B11, B12, B22, B21,
	       P1, P2, P3, P4, P5, P6, P7;
      int ip, jp;
      for (i=0; i < Np/2; i++) 
	for (j=0; j <  Np/2; j++) 
	  {
	    A11[i][j] = m[i][j];
	    jp = j+Np/2;
	    ip = i+Np/2; 
	    A12[i][j] = (dN==1 && (jp==N))?0.0:m[i][jp];
	    A21[i][j] = (dN==1 && (ip==N))?0.0:m[ip][j]; 
	    A22[i][j] = (dN==1 && (ip==N||jp==N))?0.0:m[ip][jp];
	    B11[i][j] = m1.m[i][j]; 
	    B12[i][j] = (dN==1 && (jp==N))?0.0:m1.m[i][jp];
	    B21[i][j] = (dN==1 && (ip==N))?0.0:m1.m[ip][j]; 
	    B22[i][j] = (dN==1 && (ip==N||jp==N))?0.0:m1.m[ip][jp];
	  }
      P1 = (A11+A22)*(B11+B22);
      P2 = (A21+A22)*B11;
      P3 = A11*(B12-B22);
      P4 = A22*(B21-B11);
      P5 = (A11+A12)*B22;
      P6 = (A21-A11)*(B11+B12);
      P7 = (A12-A22)*(B21+B22);
      for (i=0; i < Np/2; i++) 
	for (j=0; j < Np/2; j++) 
	  {
	    m2[i][j] = P1[i][j] + P4[i][j] - P5[i][j] + P7[i][j];
	    jp = j+Np/2;
	    ip = i+Np/2; 
	    if (jp < N)
	      m2[i][jp] = P3[i][j] + P5[i][j];
	    if (ip < N)
	      {
		m2[ip][j] = P2[i][j] + P4[i][j];
		if (jp < N)
		  m2[ip][jp] = P1[i][j] + P3[i][j] - P2[i][j] + P6[i][j];
	      }
	  }
    }
  void muladd(const pmatrixq<ntype,N>& m1, const pmatrixq<ntype,N>& m2, ntype alpha=1.0, ntype beta=0.0) 
    {
      // m=alpha*m*m1 + beta*m2
#ifdef USE_LAPACK
      if constexpr (N > NMAXMUL)
	{
	  if constexpr (is_same<ntype, double>::value)
	    {
	      wrap_dgemm(&(m[0][0]), &(m1.m[0][0]), &(m2.m[0][0]), N, alpha, beta);
	    }
	  else if constexpr (is_same<ntype, float>::value)
	    {
	      wrap_sgemm(&(m[0][0]), &(m1.m[0][0]), &(m2.m[0][0]), N, alpha, beta);
	    }
	  else
	    {
	      if (alpha==1.0 && beta==0)
		(*this) = m1*m2; 
	      else if (alpha==1.0 && beta==1.0)
		(*this) = m1*m2 + (*this);
	      else
		(*this) = alpha*m1*m2 + beta*(*this);
	    }
	}
      else
	{
       	  if (alpha==1.0 && beta==0)
	    (*this) = m1*m2; 
	  else if (alpha==1.0 && beta==1.0)
	    (*this) = m1*m2 + (*this);
	  else
	    (*this) = alpha*m1*m2 + beta*(*this);
	}
#endif
    }
#if 0
  inline void fastadd(const pmatrixq<ntype,N>& m1,const pmatrixq<ntype,N>& m2)
    {
      int i, j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] += m1.m[i][j] + m2.m[i][j];
    }
#endif
  inline pmatrixq<ntype, N> mul(const pmatrixq<ntype,N>& m1)
    {
      int i, j, k;
      pmatrixq<ntype,N> m2;
      if constexpr (N==1)
	{
	  m2.m[0][0] = pmatrixq::m[0][0]*m1.m[0][0];
	}
      else if constexpr (N==2)
	{
	  m2.m[0][0] = pmatrixq::m[0][0]*m1.m[0][0] + pmatrixq::m[0][1]*m1.m[1][0];
	  m2.m[0][1] = pmatrixq::m[0][0]*m1.m[0][1] + pmatrixq::m[0][1]*m1.m[1][1];
	  m2.m[1][0] = pmatrixq::m[1][0]*m1.m[0][0] + pmatrixq::m[1][1]*m1.m[1][0];
	  m2.m[1][1] = pmatrixq::m[1][0]*m1.m[0][1] + pmatrixq::m[1][1]*m1.m[1][1]; 
	}
      else if constexpr (N==3)
	{
	  m2.m[0][0] = pmatrixq::m[0][0]*m1.m[0][0] + pmatrixq::m[0][1]*m1.m[1][0] + pmatrixq::m[0][2]*m1.m[2][0];
	  m2.m[0][1] = pmatrixq::m[0][0]*m1.m[0][1] + pmatrixq::m[0][1]*m1.m[1][1] + pmatrixq::m[0][2]*m1.m[2][1];
	  m2.m[0][2] = pmatrixq::m[0][0]*m1.m[0][2] + pmatrixq::m[0][1]*m1.m[1][2] + pmatrixq::m[0][2]*m1.m[2][2];
	  m2.m[1][0] = pmatrixq::m[1][0]*m1.m[0][0] + pmatrixq::m[1][1]*m1.m[1][0] + pmatrixq::m[1][2]*m1.m[2][0];
	  m2.m[1][1] = pmatrixq::m[1][0]*m1.m[0][1] + pmatrixq::m[1][1]*m1.m[1][1] + pmatrixq::m[1][2]*m1.m[2][1];
	  m2.m[1][2] = pmatrixq::m[1][0]*m1.m[0][2] + pmatrixq::m[1][1]*m1.m[1][2] + pmatrixq::m[1][2]*m1.m[2][2];
	  m2.m[2][0] = pmatrixq::m[2][0]*m1.m[0][0] + pmatrixq::m[2][1]*m1.m[1][0] + pmatrixq::m[2][2]*m1.m[2][0];
	  m2.m[2][1] = pmatrixq::m[2][0]*m1.m[0][1] + pmatrixq::m[2][1]*m1.m[1][1] + pmatrixq::m[2][2]*m1.m[2][1];
	  m2.m[2][2] = pmatrixq::m[2][0]*m1.m[0][2] + pmatrixq::m[2][1]*m1.m[1][2] + pmatrixq::m[2][2]*m1.m[2][2];  
	} 
      else if constexpr (N==4)
	{
	  //cout << pmatrix<ntype,N,N>::dynamic;
	  m2.m[0][0] = pmatrixq::m[0][0]*m1.m[0][0] + pmatrixq::m[0][1]*m1.m[1][0] + pmatrixq::m[0][2]*m1.m[2][0] + pmatrixq::m[0][3]*m1.m[3][0];
	  m2.m[0][1] = pmatrixq::m[0][0]*m1.m[0][1] + pmatrixq::m[0][1]*m1.m[1][1] + pmatrixq::m[0][2]*m1.m[2][1] + pmatrixq::m[0][3]*m1.m[3][1];
	  m2.m[0][2] = pmatrixq::m[0][0]*m1.m[0][2] + pmatrixq::m[0][1]*m1.m[1][2] + pmatrixq::m[0][2]*m1.m[2][2] + pmatrixq::m[0][3]*m1.m[3][2];
	  m2.m[0][3] = pmatrixq::m[0][0]*m1.m[0][3] + pmatrixq::m[0][1]*m1.m[1][3] + pmatrixq::m[0][2]*m1.m[2][3] + pmatrixq::m[0][3]*m1.m[3][3];
	  m2.m[1][0] = pmatrixq::m[1][0]*m1.m[0][0] + pmatrixq::m[1][1]*m1.m[1][0] + pmatrixq::m[1][2]*m1.m[2][0] + pmatrixq::m[1][3]*m1.m[3][0];
	  m2.m[1][1] = pmatrixq::m[1][0]*m1.m[0][1] + pmatrixq::m[1][1]*m1.m[1][1] + pmatrixq::m[1][2]*m1.m[2][1] + pmatrixq::m[1][3]*m1.m[3][1];
	  m2.m[1][2] = pmatrixq::m[1][0]*m1.m[0][2] + pmatrixq::m[1][1]*m1.m[1][2] + pmatrixq::m[1][2]*m1.m[2][2] + pmatrixq::m[1][3]*m1.m[3][2];
	  m2.m[1][3] = pmatrixq::m[1][0]*m1.m[0][3] + pmatrixq::m[1][1]*m1.m[1][3] + pmatrixq::m[1][2]*m1.m[2][3] + pmatrixq::m[1][3]*m1.m[3][3];
	  m2.m[2][0] = pmatrixq::m[2][0]*m1.m[0][0] + pmatrixq::m[2][1]*m1.m[1][0] + pmatrixq::m[2][2]*m1.m[2][0] + pmatrixq::m[2][3]*m1.m[3][0];
	  m2.m[2][1] = pmatrixq::m[2][0]*m1.m[0][1] + pmatrixq::m[2][1]*m1.m[1][1] + pmatrixq::m[2][2]*m1.m[2][1] + pmatrixq::m[2][3]*m1.m[3][1];
	  m2.m[2][2] = pmatrixq::m[2][0]*m1.m[0][2] + pmatrixq::m[2][1]*m1.m[1][2] + pmatrixq::m[2][2]*m1.m[2][2] + pmatrixq::m[2][3]*m1.m[3][2];
	  m2.m[2][3] = pmatrixq::m[2][0]*m1.m[0][3] + pmatrixq::m[2][1]*m1.m[1][3] + pmatrixq::m[2][2]*m1.m[2][3] + pmatrixq::m[2][3]*m1.m[3][3];
	  m2.m[3][0] = pmatrixq::m[3][0]*m1.m[0][0] + pmatrixq::m[3][1]*m1.m[1][0] + pmatrixq::m[3][2]*m1.m[2][0] + pmatrixq::m[3][3]*m1.m[3][0];
	  m2.m[3][1] = pmatrixq::m[3][0]*m1.m[0][1] + pmatrixq::m[3][1]*m1.m[1][1] + pmatrixq::m[3][2]*m1.m[2][1] + pmatrixq::m[3][3]*m1.m[3][1];
	  m2.m[3][2] = pmatrixq::m[3][0]*m1.m[0][2] + pmatrixq::m[3][1]*m1.m[1][2] + pmatrixq::m[3][2]*m1.m[2][2] + pmatrixq::m[3][3]*m1.m[3][2];
	  m2.m[3][3] = pmatrixq::m[3][0]*m1.m[0][3] + pmatrixq::m[3][1]*m1.m[1][3] + pmatrixq::m[3][2]*m1.m[2][3] + pmatrixq::m[3][3]*m1.m[3][3];
	}
#ifdef USE_LAPACK
      else if constexpr (N > NMAXMUL) 
	{
	  if constexpr (is_same<ntype, double>::value)
	    wrap_dgemm(&(m2.m[0][0]), &(m[0][0]), &(m1.m[0][0]), N);
	  else if constexpr (is_same<ntype, float>::value)
	    wrap_sgemm(&(m2.m[0][0]), &(m[0][0]), &(m1.m[0][0]), N);
	  else
	    {
	      if constexpr (N > NMAXSTRA)
		{
#ifdef USE_STRASSEN
		  mulstra(m2, m1);
#else
    		  mulblock(m2, m1);
#endif
		}
	      else
		{
		  for (i=0; i < N; i++)
		    {
		      for (j=0; j < N; j++)
			{
			  m2.m[i][j]=0.0;
			  for (k=0; k < N; k++)
			    m2.m[i][j] += pmatrixq::m[i][k]*m1.m[k][j];
			}
		    }
		}
	    }
	}
#else
      else if constexpr (N > NMAXSTRA) // 40 has been chose optimizing performance on my macbook pro 13" but can be architecture dependent
	{	
#ifdef USE_STRASSEN
	  mulstra(m2, m1);
#else
	  mulblock(m2, m1);
#endif
	}
#endif
      else
	{
	  for (i=0; i < N; i++)
	    {
	      for (j=0; j < N; j++)
		{
		  m2.m[i][j]=0.0;
		  for (k=0; k < N; k++)
		    m2.m[i][j] += pmatrixq::m[i][k]*m1.m[k][j];
		}
	    }
	}
      if constexpr (N > NMAX)
	return std::move(m2);
      else
	return m2;
    }
   // transpose
  pmatrixq<ntype,N> transp()
    {
      int i, j;
      pmatrixq<ntype,N> m1;
      for (i=0; i < N; i++)
	{
	  for (j=0; j < N; j++)
	    {
	      m1[i][j] = m[j][i];
	    }
	}
      return m1;
    }
  // get i-th column vector 
  pvector<ntype,N> col(int i)
    {
      pvector<ntype,N> v;
      int j;
      for (j=0; j < N; j++)
	{
	  v[j] = m[j][i];
	}
      return v;
    }
  // get i-th row vector 
  pvector<ntype,N> row(int i)
    {
      pvector<ntype,N> v;
      int j;
      for (j=0; j < N; j++)
	{
	  v[j] = m[i][j];
	}
      return v;
    }
  
  void show(const char* str)
    {
      int i, j;
      if (str!=NULL)
	cout << str;
      cout << "{";
      for (i=0; i < N; i++)
	{
	  cout << "{";
	  for (j=0; j < N; j++)
	    { 
	      cout << setprecision(std::numeric_limits<ntype>::digits10) << m[i][j];
		if (j < N-1)
		  cout << ",";
	      }
	    cout << "}";
	    if (i < N-1)
	      cout << ",\n";
	  }
	cout << "}\n";
      }
    void show(void)
      {
	show(NULL);
      }

    int size()
      {
	return N;
      } 
#if 0
    pmatrixq()
      {
	//n=N;
	//maxdigits=std::numeric_limits<ntype>::digits10;
      }
    ~pmatrixq()
      {
      }
#endif
};
#ifdef MAT_LAZY_EVAL
// ADDITION
// add expression template with point at the right
  template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(AnOpM<ntype, N, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,N> const& p)
{
  return AnOpM<ntype, N, MOpTypes::MatPlusMat, AnOpM<ntype, N, tipo, Lhs, Rhs>, pmatrixq<ntype,N>>(lhs, p);
} 
// add expression template with point at the left
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(pmatrixq<ntype,N> const& p, AnOpM<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpM< ntype, N, MOpTypes::MatPlusMat, pmatrixq<ntype,N>, AnOpM<ntype, N, tipo, Lhs, Rhs> >(p, rhs);
}

// add two points, yield a expression template    
// se N < NSTA restituisce pmatrixq e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpM e usa la lazy evaluation
  template <typename ntype, int N> 
  // conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  typename std::conditional<(N<=matpars::NLAZY),pmatrixq<ntype,N>,AnOpM<ntype, N, MOpTypes::MatPlusMat, pmatrixq<ntype,N>, pmatrixq<ntype,N>>>::type
inline operator+(pmatrixq<ntype,N> const& lhs, pmatrixq<ntype,N> const& rhs) 
{
  if constexpr (N <= pmatrixq<ntype,N>::NLAZY)
    {
      return lhs.addition(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << NSTA << "\n";
      return AnOpM<ntype, N, MOpTypes::MatPlusMat, pmatrixq<ntype,N>, pmatrixq<ntype,N>>(lhs, rhs); // lazy evaluation
    }
}
// caso aggiunto secondo un post in https://stackoverflow.com/questions/414243/lazy-evaluation-in-c
template <typename ntype, int N, int tipo, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator+(const AnOpM<ntype, N, tipo, LLhs, LRhs> & leftOperandconst, const AnOpM<ntype, N, tipo, RLhs, RRhs> & rightOperand)
{
  return  AnOpM<ntype, N, MOpTypes::MatPlusMat, AnOpM<ntype, N, tipo, LLhs, LRhs>, AnOpM<ntype, N, tipo, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
// SUBTRACTION
// add expression template with point at the right
  template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator-(AnOpM<ntype,N, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,N> const& p)
{
  return AnOpM<ntype, N, MOpTypes::MatMinusMat, AnOpM<ntype, N, tipo, Lhs, Rhs>, pmatrixq<ntype,N>>(lhs, p);
} 
// add expression template with point at the left
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator-(pmatrixq<ntype,N> const& p, AnOpM<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpM< ntype, N, MOpTypes::MatMinusMat, pmatrixq<ntype,N>, AnOpM<ntype, N, tipo, Lhs, Rhs> >(p, rhs);
}

// add two points, yield a expression template    
// se N < NSTA restituisce pmatrixq e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpM e usa la lazy evaluation
  template <typename ntype, int N> 
  // conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  typename std::conditional<(N<=matpars::NLAZY),pmatrixq<ntype,N>,AnOpM<ntype, N, MOpTypes::MatMinusMat, pmatrixq<ntype,N>, pmatrixq<ntype,N>>>::type
inline operator-(pmatrixq<ntype,N> const& lhs, pmatrixq<ntype,N> const& rhs) 
{
  if constexpr (N <= matpars::NLAZY)
    {
      return lhs.subtraction(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << NSTA << "\n";
      return AnOpM<ntype, N, MOpTypes::MatMinusMat, pmatrixq<ntype,N>, pmatrixq<ntype,N>>(lhs, rhs); // lazy evaluation
    }
}
// caso aggiunto secondo un post in https://stackoverflow.com/questions/414243/lazy-evaluation-in-c
  template <typename ntype, int N, int tipo, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator-(const AnOpM<ntype, N, tipo, LLhs, LRhs> & leftOperandconst, const AnOpM<ntype, N, tipo, RLhs, RRhs> & rightOperand)
{
  return  AnOpM<ntype, N, MOpTypes::MatMinusMat, AnOpM<ntype, N, tipo, LLhs, LRhs>, AnOpM<ntype, N, tipo, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
   
///
// MULTIPLICATION OF TWO AnOpM
//
  template <typename ntype, int N, int tipo, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator*(const AnOpM<ntype, N, tipo, LLhs, LRhs> & leftOperand, const AnOpM<ntype, N, tipo, RLhs, RRhs> & rightOperand)
{
  pmatrixq<ntype, N> mL, mR;
  for (int i=0; i < N; i++)
    for (int j=0; j < N; j++)
      {
    	mR[i][j] = leftOperand.get_m(i,j);
	mL[i][j] = rightOperand.get_m(i,j);
      }
  return mR*mL;
}
// AnOpM times pmatrixq
//
// left AnOpM
   template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(AnOpM<ntype,N, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,N> const& mR)
{
  pmatrixq<ntype, N> mL;
  for (int i=0; i < N; i++)
    for (int j=0; j < N; j++)
      {
    	mL[i][j] = lhs.get_m(i,j);
      }
  return mL*mR;
} 
// right AnOpM
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>  
inline auto operator*(pmatrixq<ntype,N> const& mL, AnOpM<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  pmatrixq<ntype, N> mR;
  for (int i=0; i < N; i++)
    for (int j=0; j < N; j++)
      {
    	mR[i][j] = rhs.get_m(i,j);
      }
  return mL*mR;
}
// vector times AnOpM
template <typename ntype, int N, int tipo, typename Lhs, typename Rhs>
inline auto operator *(const pvector<ntype,N> vL, AnOpM<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  pmatrixq<ntype, N> mR;
  for (int i=0; i < N; i++)
    for (int j=0; j < N; j++)
      {
    	mR[i][j] = rhs.get_m(i,j);
      }
  return vL*mR;
}
// AnOpM times vector
template <typename ntype, int N, int tipo, typename Lhs, typename Rhs>
inline auto operator *(AnOpM<ntype, N, tipo, Lhs, Rhs> const& lhs, const pvector<ntype,N> vR) 
{
  pmatrixq<ntype, N> mL;
  for (int i=0; i < N; i++)
    for (int j=0; j < N; j++)
      {
    	mL[i][j] = lhs.get_m(i,j);
      }
  return mL*vR;
}
// AnOpM times scalar
template <typename ntype, int N, int tipo, typename Lhs, typename Rhs>
 inline auto operator *(AnOpM<ntype, N, tipo, Lhs, Rhs> const& lhs, const double& param) 
    {
      int i,j;
      pmatrixq<ntype,N> mL;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  mL[i][j] = lhs.get_m(i,j)*param;
      return mL;
    }
  
// AnOpM divided by scalar
template <typename ntype, int N, int tipo, typename Lhs, typename Rhs>
 inline auto operator /(AnOpM<ntype, N, tipo, Lhs, Rhs> const& lhs, const double& param) 
    {
      int i,j;
      pmatrixq<ntype,N> mL;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  mL[i][j] = lhs.get_m(i,j)/param;
      return mL;
    }
#endif
typedef pmatrixq<double,3> pmatrix3d;
typedef pmatrixq<int,3>    pmatrix3i;
#endif
