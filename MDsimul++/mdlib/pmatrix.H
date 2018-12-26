#ifndef _PMATRIX_
/* N.B. current status: better than armadillo and glm */
#define _PMATRIX_
#define USE_STRASSEN
// matrici NxM
#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include "pvector_mp.H"
#define OPT_INV
#define MAT_USE_LAPACK
#define MAT_COMMA_INIT
//#define LAZY_MUL
//N.B. Lazy evaluation of matrix multiplication has not been fully implemented yet
#ifdef MAT_USE_LAPACK
#undef __APPLE__ // undefining __APPLE__ one can prevent the usage of macosx accelerated framework switching to standard lapack blas
#ifdef __APPLE__
//  /System/Library/Frameworks/vecLib.framework/Headers/cblas.h
//  g++ -framework accelerate to compile
//#include <cblas.h>
//#include <lapacke.h>
#include <Accelerate/Accelerate.h>
#endif
#define MAT_MOVE_SEMANTIC
#define MAT_LAZY_EVAL
#define USE_GAUSS_INV
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
extern void dsyev_(char *JOBZ, char *UPLO, int *N, double *AT, int *LDA, double *B, double *WORK, int *LWORK, int *INFO);      
extern void dgesv_(int *N, int *NRHS, double* AT, int *LDA, int *pivot, double *B, int *LDB, int *INFO);      
extern void ssyev_(char *JOBZ, char *UPLO, int *N, float *AT, int *LDA, float *B, float *WORK, int *LWORK, int *INFO);      
extern void sgesv_(int *N, int *NRHS, float* AT, int *LDA, int *pivot, float *B, int *LDB, int *INFO);      
}
#endif
//Matrix times vector
void wrap_dgemv(char ta, double *Y, const double *A, const double *X, int n, double alpha=1.0, double beta=0.0)
{
  // Y = A*X + Y
  double ALPHA, BETA;
  int LDA, M,N, incx=1, incy=1;
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
  int LDA, M,N, incx=1, incy=1;
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
  int LDA, LDB, LDC, M,N, K;
  char ta='n', tb='n';
  M=N=K=n;
  ALPHA=alpha;
  BETA=beta; 
  LDA=LDB=LDC=n;
#ifdef __APPLE__
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M,N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);      
#else
  // LAPACK and BLAS use column major (matrix is a sequence of columns, i.e.
  // ABC|DEF|GHI is a matrix with columns ABC, DEF and GHI) convention hence we have to swap A and B
  dgemm_(&ta, &tb, &M, &N, &K, &ALPHA, B, &LDA, A, &LDB, &BETA, C, &LDC);      
#endif
}
void wrap_sgemm(float *C, const float *A, const float *B, int n, float alpha=1.0, float beta=0.0)
{
  //C=A*B + beta*C if beta=0 C is not used as input
  int LDA, LDB, LDC, M,N, K;
  char ta='n', tb='n';
  float ALPHA, BETA;
  M=N=K=n;
  ALPHA=alpha;
  BETA=beta; 
  LDA=LDB=LDC=n;
#ifdef __APPLE__
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M,N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);      
#else
  // LAPACK and BLAS use column major convention hence we have to swap A and B
  sgemm_(&ta, &tb, &M, &N, &K, &ALPHA, B, &LDA, A, &LDB, &BETA, C, &LDC);      
#endif
}
// eigenvalues eigenvectors
/* solve a linear system of equations */
void wrap_dgesv(double *a, double* x, int N, int *ok)
{
  int i, j, c1;
  auto pivot = new int[N];
  auto AT = new double[N*N];
  for (i=0; i<N; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix from column major to row major*/
      for(j=0; j<N; j++) AT[j+N*i]=a[j*N+i];		
    }						
  c1 = N;
  // this routine uses LU decomposition hence if MAT_USE_LAPACK is not defined
  // the LU decomposition of Numerical Recipe can be used
  dgesv_(&c1, &c1, AT, &c1, pivot, x, &c1, ok);      
  delete [] pivot;
  delete [] AT;
}
/* solve a linear system of equations */
void wrap_sgesv(float *a, float* x, int N, int *ok)
{
  int i, j, c1;
  auto pivot = new int[N];
  auto AT = new float[N*N];
  for (i=0; i<N; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix from column major to row major*/
      for(j=0; j<N; j++) AT[j+N*i]=a[j*N+i];		
    }						
  c1 = N;
  // this routine uses LU decomposition hence if MAT_USE_LAPACK is not defined
  // the LU decomposition of Numerical Recipe can be used
  sgesv_(&c1, &c1, AT, &c1, pivot, x, &c1, ok);      
  delete [] pivot;
  delete [] AT;
}
/* find eigenvectors (in b) and eigenvalues (in 1d array x) */
void wrap_dsyev(double *a, double *b, double *x, int N, int *ok, bool ev=false)// if ev=1 find also eigenvectors
{
  char JOBZ, UPLO;
  double *AT, *work;
  int i, j, c1, c3;

  if (ev==true)
    JOBZ='V';
  else
    JOBZ='N';
  UPLO='U';
  AT = new double[N*N];
  c3 = 3*N-1;
  work = new double[c3];
  for (i=0; i<N; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<N; j++)
       	AT[j+N*i]=a[j*N+i];		
    }						
  c1 = N;
  dsyev_(&JOBZ, &UPLO, &c1, AT, &c1, x, work, &c3, ok);      
  if (ev==true)
    {
      for (i=0; i<N; i++)		/* to call a Fortran routine from C we */
	{				/* have to transform the matrix */
	  for(j=0; j<N; j++) 
	    b[i*N+j]=AT[j+N*i];		
	}	
    }
  //if (*ok != 0)
    //printf("not ok (%d)\n", *ok);
  delete [] AT;
  delete [] work;
}
/* find eigenvectors (in b) and eigenvalues (in 1d array x) */
void wrap_ssyev(float *a, float *b, float *x, int N, int *ok, bool ev=false)// if ev=1 find also eigenvectors
{
  char JOBZ, UPLO;
  float *AT, *work;
  int i, j, c1, c3;

  if (ev==true)
    JOBZ='V';
  else
    JOBZ='N';
  UPLO='U';
  AT = new float[N*N];
  c3 = 3*N-1;
  work = new float[c3];
  for (i=0; i<N; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<N; j++)
       	AT[j+N*i]=a[j*N+i];		
    }						
  c1 = N;
  ssyev_(&JOBZ, &UPLO, &c1, AT, &c1, x, work, &c3, ok);      
  if (ev==true)
    {
      for (i=0; i<N; i++)		/* to call a Fortran routine from C we */
	{				/* have to transform the matrix */
	  for(j=0; j<N; j++) 
	    b[i*N+j]=AT[j+N*i];		
	}	
    }
  //if (*ok != 0)
    //printf("not ok (%d)\n", *ok);
  delete [] AT;
  delete [] work;
}
#endif
#ifdef MAT_LAZY_EVAL
//REMARK on Lazyness Logic:
//each class is an operation with operands Lhs and Rhs, the instance of the class (i.e. the object created
//when an operator is parsed by the compiler)
//stores the values of the operands and only when the method get_m is called the operation is performed.
//Evaluation is triggered as assign operator of pmatrixq class gets called (see below) or constructor is called 
//with an AnOpM argument (i.e. a casting from AnOpM to pmatrixq type).
//
//Matrix operation types
struct MOpTypes
{  
  static const int MatPlusMat=0, MatMinusMat=1, MatTimesScal=2, ScalTimesMat=3,
	       MatDivScal=4, MatTimesMat=5;
};
template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
struct AnOpM 
{
  Lhs const& lhs;
  Rhs const& rhs;
  AnOpM(Lhs const& lhs, Rhs const& rhs):lhs(lhs), rhs(rhs) {
    // empty body only initialize operands 
  }
  // LAZY operations
  inline ntype get_m(int i, int j) const {
    if constexpr (tipo == MOpTypes::MatPlusMat)// mat + mat
      {
	return lhs.get_m(i,j) + rhs.get_m(i,j);
      }
    else if constexpr (tipo == MOpTypes::MatMinusMat) // mat - mat
      {
	return lhs.get_m(i,j) - rhs.get_m(i,j);
      }
    else if constexpr (tipo == MOpTypes::MatTimesScal) // mat * scalar 
      {
	return lhs.get_m(i,j)*rhs;
      }
    else if constexpr (tipo == MOpTypes::ScalTimesMat) // scalar * mat 
      {
	return lhs*rhs.get_m(i,j);
      }
    else if constexpr (tipo == MOpTypes::MatDivScal) // mat / scalar 
      {
	return lhs.get_m(i,j)/rhs;
      }
#ifdef LAZY_MUL
    // N.B. lazy multiplication does not improve performance and it can provide wrong results
    // if the assigned matrix appears also on right hand side of assignment
    else if constexpr (tipo == MOpTypes::MatTimesMat && NT==2) 
      {
	// i,j element of resulting matrix
	ntype sum=0;
    	for (auto k=0; k < lhs.get_N(); k++) 
	  sum +=lhs.get_m(i,k)*rhs.get_m(k,j);
	return sum;
      }
#endif
    return 0;
  }
  int get_N() const { return lhs.get_N(); }
  Lhs const& get_lhs() const { return lhs; }//unused
  Rhs const& get_rhs() const { return rhs; }
};
#endif
// static base class
template <class ntype,int NT, int MT> class pmatrix {
  //int nr, nc;
public:
  ntype m[NT][MT];
  static const int N=NT, M=MT;
  constexpr static int dynamic = false;
  pmatrix() = default;
  pmatrix(int NN, int MM): pmatrix()
    {
      NN=NT;// just to avoid warnings
      MM=MT;
      //empty body
    }
};
// dynamic base class with static dimension (with move assignment and move constructor properly defined)
template <class ntype,int NT, int MT> class pmatrixdn {
  //int nr, nc;
public:
  ntype (*m)[MT]; // N riga M colonna
  static const int N=NT, M=MT;
  // copy assignment
 pmatrixdn<ntype,NT,MT>& operator=(const pmatrixdn<ntype,NT,MT>& m1)
    {
      int i,j;
      for (i=0; i < NT; i++)
	for (j=0; j < NT; j++)
	  m[i][j] = m1.m[i][j];
      return (*this);
    }
#ifdef MAT_MOVE_SEMANTIC
  // move assignment
  pmatrixdn<ntype,NT,MT>& operator=(pmatrixdn<ntype,NT,MT>&& m1)
    {
      swap(m, m1.m);
      return (*this);
    }
#endif
  // copy constructor 
  pmatrixdn(const pmatrixdn<ntype,NT,MT>& m1)
    {
      m = new ntype[NT][MT];
      (*this) = m1;
    }
#ifdef MAT_MOVE_SEMANTIC
  // move constructor
  pmatrixdn(pmatrixdn<ntype,NT,MT>&& m1)
    {
      (*this).m = m1.m;
      m1.m=nullptr;
    }
#endif
  // default constructor
  pmatrixdn(int NN, int MM)
    {
      NN=NT;
      MM=MT;
      m = new ntype[NT][MT];
    }

  // default constructor
  pmatrixdn()
    {
      m = new ntype[NT][MT];
    }
  // destructor
  ~pmatrixdn()
    {
      delete[] m;
    }
};	
// dynamic base class using ptr to ptr(with move assignment and move constructor properly defined)
template <class ntype,int NT, int MT> class pmatrixdp {
  //int nr, nc;
public:
  ntype **m; // N riga M colonna
  ntype *pm;
  int N, M;
  // copy assignment
  pmatrixdp<ntype,NT,MT>& operator=(const pmatrixdp<ntype,NT,MT>& m1)
    {
      int i,j;
      for (i=0; i < N; i++)
	for (j=0; j < M; j++)
	  m[i][j] = m1.m[i][j];
      return (*this);
    }
#ifdef MAT_MOVE_SEMANTIC
  // move assignment
  pmatrixdp<ntype,NT,MT>& operator=(pmatrixdp<ntype,NT,MT>&& m1)
    {
      //cout << "move assignment\n";
      swap(m, m1.m);
      swap(pm, m1.pm);
      return (*this);
    }
#endif
  // copy constructor 
  pmatrixdp(const pmatrixdp<ntype,NT,MT>& m1)
    {
      N=m1.N;
      M=m1.M;
      //ntype *pt = (ntype *)malloc(sizeof(ntype)*N*M);
      //m = (ntype**)malloc(sizeof(ntype*)*N);
      pm = new ntype[N*M];//(ntype *)malloc(sizeof(ntype)*N*M);
      m = new ntype*[N];//(ntype**)malloc(sizeof(ntype*)*N);
 
      for (auto i=0; i < N; i++)
	{
	  m[i] = pm + i*M;
	  //cout << "m[" << i << "]=" << m[i] << "\n";
	}

      for (auto i=0; i < N; i++)
	for (auto j=0; j < M; j++)
	  m[i][j] = m1.m[i][j];
    }
#ifdef MAT_MOVE_SEMANTIC
  // move constructor
  pmatrixdp(pmatrixdp<ntype,NT,MT>&& m1)
    {
      N=m1.N;
      M=m1.M;
      (*this).m = m1.m;
      (*this).pm = m1.pm;
      m1.pm = nullptr;
      m1.m=nullptr;
    }
#endif
  // default constructor deleted: one has to supply size!
  pmatrixdp() = default;
#if 0
    {
      ntype *pt = malloc(sizeof(ntype)*N*M);
      m = malloc(sizeof(ntype*)*N);
      for (auto i=0; i < N; i++)
	m[i] = pt + i*sizeof(ntype)*M;
    }
#endif
  pmatrixdp(int NN, int MM)
    {
      N=NN;
      M=MM;
      pm = new ntype[N*M];//(ntype *)malloc(sizeof(ntype)*N*M);
      m = new ntype*[N];//(ntype**)malloc(sizeof(ntype*)*N);
      for (auto i=0; i < N; i++)
	{
	  m[i] = pm + i*M;
	  //cout << "m[" << i << "]=" << m[i] << "\n";
	}
    }
  void deallocate(void)
    {
      delete [] pm;
      delete [] m;
    }
  void allocate(int NN, int MM)
    {
      N=NN;
      M=MM;
      pm = new ntype[N*M];
      m = new ntype*[N];
      for (auto i=0; i < N; i++)
	{
	  m[i] = pm + i*M;
	  //cout << "m[" << i << "]=" << m[i] << "\n";
	}
    }
  void resize(int NN, int MM)
    {
      ntype* pmt = new ntype[NN*MM];
      ntype** mt = new ntype*[NN];
      for (auto i=0; i < NN; i++)
	{
	  mt[i] = pmt + i*MM;
	  //cout << "m[" << i << "]=" << m[i] << "\n";
	}
      // initialize with old values
      for (auto i=0; i < NN; i++)
	for (auto j=0; j < MM; j++)
	  {
	    mt[i][j] = (i < N && j < M)?m[i][j]:0.0;
	  }
      delete [] m;
      delete [] pm;
      N=NN;
      M=MM;
      m = mt;
      pm = pmt;
    }
  // destructor
  ~pmatrixdp()
    {
      delete [] pm;
      delete [] m;
    }
};	
template <class ntype, int NT> using pmatbasedyn = 
typename std::conditional<(NT>0), pmatrixdn <ntype, NT, NT>,
	 pmatrixdp <ntype, NT, NT>>::type;

// dynamic base class
template <class ntype,int NT> class pmatrixd: public pmatbasedyn<ntype,NT> {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  pmatrixd<ntype,NT>(int NN): pmatbasedyn<ntype,NT>(NN,NN)
  {
    // empty body
  }
  pmatrixd<ntype,NT>(): pmatbasedyn<ntype,NT>()
  {
    // empty body
  }
  void allocate(int NN)
    {
      if constexpr (NT < 0)
	pmatbasedyn<ntype,NT>::allocate(NN,NN);
    }
  void deallocate(int NN)
    {
      if constexpr (NT < 0)
	pmatbasedyn<ntype,NT>::deallocate(NN,NN);
    }
  void resize(int NN)
    {
      if constexpr (NT < 0)
	pmatbasedyn<ntype,NT>::resize(NN,NN);
    }


};
// matrici quadrate NxN
using namespace std;
struct matpars
{
  // params for my class
  static const int NMAX = 45;     // if N>NMAX use dynamic allocation of matrices
  static const int NMAXSTRA = 40; // if N > NMAXSTRA use Strassen multiplication algorithm
  static const int NMAXINV = 25;  // if N>NAMXINV use LAPACK LU DECOMP
  static const int NMAXMUL = 25;  // if N>NMAXMUL use BLAS multiplication routines
  static const int NLAZY = 10;     // value below which does not use lazy evaluation for addition
  static const int Dynamic = -1;
  static const int DynNoLazy = -2;
};
// for large matrices use dynamic allocation to avoid stack overflow!
// std::conditional metafunction returns either one of the two types depending on the condition (N>matpars::NMAX in this 
// case)
template <class ntype, int NT=-1> using pmatrixb = 
typename std::conditional<(NT>matpars::NMAX||NT <=0), pmatrixd <ntype, NT>,
	 pmatrix <ntype, NT, NT>>::type;

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

template <class ntype, int NT> class pmatrixq: public pmatrixb<ntype,NT>, public matpars
{
#ifdef MAT_COMMA_INIT
  unsigned int curidx;
#endif
public:
  // solve the linear system (*this)*x = b (b is replaced by x)
  // if b==NULL then returns just the inverse of (*this)

  using pvec = pvector<ntype,NT>;
  template <class T=pvec>
    pmatrixq<ntype,NT> gaussj(typename  std::conditional< std::is_same<T,pvec>::value, T&, int>::type b) const
      {
        // a matrix with m columns (instead of a vector) could be supplied to solve many linear systems at once  
        int i,icol,irow,j,k,l,ll;
        ntype big,dum,pivinv;
        static const int n=N;
        //static const int m=N;
        pvector<int,NT> indxc,indxr,ipiv;
        pmatrixq<ntype,NT> a;
        if constexpr (NT < 0)
          {
            a.allocate(N);
            indxc.allocate(N);
            indxr.allocate(N);
            ipiv.allocate(N);
          }
        a = *this;
        for (j=0;j<n;j++) ipiv[j]=0;
        for (i=0;i<n;i++) {
          big=0.0;
          for (j=0;j<n;j++)
            if (ipiv[j] != 1)
              for (k=0;k<n;k++) {
                if (ipiv[k] == 0) {
                  if (fabs(a[j][k]) >= big) {
                    big=fabs(a[j][k]);
                    irow=j;
                    icol=k;
                  }
                }
              }
          ++(ipiv[icol]);
          if (irow != icol) {
            for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
            if constexpr ( std::is_same<pvec,T>::value)
              SWAP(b[irow],b[icol]);
            // a matrix with m columns could be supplied to solve many linear systems at once  
            //for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
          }
          indxr[i]=irow;
          indxc[i]=icol;
          if (a[icol][icol] == 0.0)
            {
              //nrerror("gaussj: Singular Matrix");
              // throw an exception here instead  
              cout << "gaussj: Singular Matrix\n";
              return a;
            }
          pivinv=1.0/a[icol][icol];
          a[icol][icol]=1.0;
          for (l=0;l<n;l++) a[icol][l] *= pivinv;
          if constexpr ( std::is_same<pvec,T>::value)
            b[icol] *= pivinv;
          // a matrix with m columns could be supplied to solve many linear systems at once  
          //for (l=0;l<m;l++) b[icol][l] *= pivinv;
          for (ll=0;ll<n;ll++)
            if (ll != icol) {
              dum=a[ll][icol];
              a[ll][icol]=0.0;
              for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
              if constexpr ( std::is_same<pvec,T>::value)
                b[ll][l] -= b[icol]*dum;
              // a matrix with m columns could be supplied to solve many linear systems at once  
              //   for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
            }
        }
        for (l=n-1;l>=0;l--) {
          if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
              SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
        return a;
      }
  pmatrixq<ntype,NT> gaussj() const
    {
      //Overloaded version with no right-hand sides. Replaces a by its inverse. 
      //MatDoub b(a.nrows(),0); //Dummy vector with zero columns.
      return gaussj<int>(0);
    }
  inline pmatrixq<ntype,NT> ludcmp(pvector<int,NT>& indx, ntype& d) const
    {
      /* LU decomposition routine  borrowed from Numerical Recipe */
      /* return the LU decomposition of calling matrix (i.e. *this) */
      pmatrixq <ntype,NT> m1;
      if constexpr (NT < 0)
        m1.allocate(N);
      static const double TINY=1E-40;
      int i,imax=-1,j,k;
      ntype big,dum,sum,temp; 
      ntype vv[N]; /* vv stores the implicit scaling of each row.*/
      d=1.0; /* No row interchanges yet. */
      m1 = (*this);
      for (i=0;i<N;i++) 
        { 
          /* Loop over rows to get the implicit scaling information.*/ 
          big=0.0; 
          for (j=0;j<N;j++)
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
      for (j=0;j<N;j++) 
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
          for (i=j;i<N;i++) 
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
              for (k=0;k<N;k++) 
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
          if (j != N) 
            { 
              /* Now,  nally, divide by the pivot element.*/
              dum=1.0/(m1[j][j]); 
              for (i=j+1;i<N;i++) 
                m1[i][j] *= dum; 
            } 
        } 
      return m1;
    }

  inline pvector<ntype,NT> lubksb(pvector <int,NT> indx, pvector<ntype,NT> b) const
    { 
      /* the calling matrix (*this) must be the LU decomposition of a given matrix M
       * with indx the raw permutations done to obtain LU from M, 
       * the return vector is the solution of A*x + b = 0 */
      int i,ii=0,ip,j; 
      ntype sum; 
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
            sum -= m[i][j]*b[j]; 
          b[i]=sum/m[i][i]; 
          /* Store a component of the solution vector X. */ 
        } /* All done! */
      return b;
    }
  using pmatrixb<ntype,NT>::m;
  using pmatrixb<ntype,NT>::N;
  // Build identity matrix, simple as that
#if 1
  pmatrixq<ntype,NT>(int NN): pmatrixb<ntype,NT>(NN)
  {
    // empty body
  }
#endif
  pmatrixq<ntype,NT>()= default; //pmatrixb<ntype,NT>()
  inline pmatrixq<ntype,NT> I() const
    {
      pmatrixq<ntype,NT> MI;
      if constexpr (NT < 0)
        MI.allocate(N);
      for (auto i=0; i < N; i++)
        for (auto j=0; j < N; j++)
          MI[i][j] = (i==j)?1.0:0.0;
      return MI;
    }

#ifdef MAT_LAZY_EVAL
  // -= and += are overloaded in such a way that if rvalue is an AnOpM (lazy class)
  // they trigger lazy expression evalution 
  template<typename Lhs, typename Rhs, int tipo>
    inline pmatrixq<ntype,NT>& operator -= (AnOpM<ntype,NT, tipo, Lhs, Rhs> const& op)
      {
        int i, j;
        for (i=0; i < N; i++)
          for (j=0; j < N; j++)
            m[i][j] -= op.get_m(i,j);
        return (*this);
      }

  template<typename Lhs, typename Rhs, int tipo>
    inline pmatrixq<ntype,NT>& operator += (AnOpM<ntype,NT, tipo, Lhs, Rhs> const& op)
      {
        int i, j;
        for (i=0; i < N; i++)
          for (j=0; j < N; j++)
            m[i][j] += op.get_m(i,j);
        return (*this);
      }
  // assignment operator triggers evaluation of lazy expressions
  template<typename Lhs, typename Rhs, int tipo>
    pmatrixq<ntype,NT>& operator=(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
        for (int j=0; j < N; j++)
          m[i][j] = op.get_m(i,j);
      return (*this);
    }
  template<typename Lhs, typename Rhs, int tipo>
    pmatrixq<ntype,NT>(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
        for (int j=0; j < N; j++)
          m[i][j] = op.get_m(i,j);
    }

  //pmatrixq<ntype,NT>() = default;
  inline ntype get_m(int i, int j) const {return m[i][j];}
  inline void* get_mptr() const {return m;}
  inline int get_N() const {return N;}
#endif

#ifndef MAT_LAZY_EVAL
  // if lazyness is disabled (e.g. matrices are small) these functions get called
  inline pmatrixq<ntype,NT> operator+(const pmatrixq<ntype,NT>& m1) const
    {
      return addition(m1);
    }
  inline pmatrixq<ntype,NT> operator-(const pmatrixq<ntype,NT>& m1) const
    {
      return subtraction(m1);
    }
#endif
  // addition
  inline pmatrixq<ntype,NT> addition(const pmatrixq<ntype,NT>& m1) const
    {
      if constexpr (pmatrixq<ntype,NT>::dynamic==false)
        {
          pmatrixq<ntype,NT> m2;
          int i, j;
          for (i=0; i < N; i++)
            for (j=0; j < N; j++)
              {
#if 1
                m2.m[i][j] = m[i][j] + m1.m[i][j];
#else
                // NOTA: const_cast rimuove l'attributo const ad m1 e rende possibile l'applicazione dell'operatore [] 
                // che restituisce ntype*, che non è const (altimenti si ottiene un errore nella compilazione).
                (const_cast<pmatrixq<ntype,NT>&>(m2))[i][j] = m[i][j] + (const_cast<pmatrixq<ntype,NT>&>(m1))[i][j];
#endif
              }
          return m2;
        }
      else
        {
          pmatrixq<ntype,NT> m2;
          if constexpr (NT < 0)
            m2.allocate(N);
          int i, j;
          for (i=0; i < N; i++)
            for (j=0; j < N; j++)
              {
                m2.m[i][j] = m[i][j] + m1.m[i][j];
              }
          if (N > NMAX)
            return std::move(m2);
          else
            return m2;
        }
    }
  inline pmatrixq<ntype,NT> subtraction(const pmatrixq<ntype,NT>& m1) const
    {
      int i, j;
      if constexpr (pmatrixq<ntype,NT>::dynamic==false)
        {
          pmatrixq<ntype,NT> m2;
          for (i=0; i < N; i++)
            for (j=0; j < N; j++)
              {
#if 1
                m2.m[i][j] = m[i][j] - m1.m[i][j];
#else
                // NOTA: const_cast rimuove l'attributo const ad m1 e rende possibile l'applicazione dell'operatore [] 
                // che restituisce ntype*, che non è const (altimenti si ottiene un errore nella compilazione).
                (const_cast<pmatrixq<ntype,NT>&>(m2))[i][j] = m[i][j] - (const_cast<pmatrixq<ntype,NT>&>(m1))[i][j];
#endif
              }
          return m2;
        }
      else
        {
          pmatrixq<ntype,NT> m2;
          if constexpr (NT < 0)
            m2.allocate(N);

          for (i=0; i < N; i++)
            for (j=0; j < N; j++)
              {
#if 1
                m2.m[i][j] = m[i][j] - m1.m[i][j];
#else
                // NOTA: const_cast rimuove l'attributo const ad m1 e rende possibile l'applicazione dell'operatore [] 
                // che restituisce ntype*, che non è const (altimenti si ottiene un errore nella compilazione).
                (const_cast<pmatrixq<ntype,NT>&>(m2))[i][j] = m[i][j] - (const_cast<pmatrixq<ntype,NT>&>(m1))[i][j];
#endif
              }
          if (N > NMAX)
            return std::move(m2);
          else 
            return m2;
        }
    }

  inline ntype detgenNR() const
    {
      ntype d;
      pmatrixq<ntype,NT> m1;
      int j;
      pvector<int,NT> indx;
      if constexpr (NT < 0)
        {
          m1.allocate(N);
          indx.allocate(N);
        }
      // ludcmp return the LU decomposition of actual matrix (*this)
      m1=ludcmp(indx,d);
      for(j=0;j<N;j++) 
        d *= m1[j][j];
      return d;
    }
  inline ntype detgen() const
    {
#ifdef MAT_USE_LAPACK
      // LAPACK and BLAS are available only for float and double
      if constexpr ( std::is_same<ntype, double>::value)
        {
          double d=1.0;
          int j;	
          pmatrixq<ntype,NT> om;
          int *IPIV = new int[N];
          if constexpr (NT < 0)
            {
              om.allocate(N);
            }
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
      else if constexpr ( std::is_same<ntype, float>::value)
        {
          double d=1.0;
          int j;	
          pmatrixq<ntype,NT> om;
          int *IPIV = new int[N];
          if constexpr (NT < 0)
            {
              om.allocate(N);
            }
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
  inline pmatrixq<ntype,NT> invLUNR(void) const
    {
      pvector <ntype,NT> col;
      pmatrixq<ntype,NT> LU;
      pmatrixq<ntype,NT> om;
      int m2, m1;
      double d;
      pvector<int,NT> indx;
      if constexpr (NT < 0)
        {
          LU.allocate(N);
          col.allocate(N);
          indx.allocate(N);
          om.allocate(N);
        }	
      LU=ludcmp(indx, d); 
      // inversion of a matrix through numerical recipes routines
      for(m2=0;m2<N;m2++) 
        { 
          for(m1=0;m1<N;m1++) 
            col[m1]=0.0; 
          col[m2]=1.0; 
          col=LU.lubksb(indx,col);
          for(m1=0;m1<N;m1++) 
            om[m1][m2]=col[m1]; 
        }
      return om;
    }
  inline ntype det() const
    {
      // for small matrices determinant is evaluated explicitly 
      // N.B. this first condition apparently is not needed
      // but in practice compiler optimizations can get fooled
      // if removed!
      ntype det;
      if (N <= 4)
        {
          if (N==1)
            {
              return m[0][0];
            }
          else if (N==2)
            {
              return m[0][0]*m[1][1]-m[0][1]*m[1][0];
            }
          else if (N==3)
            {
              return -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + 
                m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - 
                m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2];
            } 
          else if (N == 4)
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
        }
      else 
        // use numerical recipe or LAPACK LU decomposition (see detgen)
        det=detgen();
      return det;
    }

  inline pmatrixq<ntype,NT> invgen() const
    {
      // ludcmp returns the LU decomposition of actual matrix (*this)
      // where indx will store the raw permutiation effected by partial pivoting (see Numerical Recipe
      // for more details c2.3)
#ifdef MAT_USE_LAPACK
      if (N > NMAXINV)
        {
          if constexpr ( std::is_same<ntype, double>::value)
            {    
              // Transpose(Inverse(M))=Inverse(Transpose(M))
              // hence we do not have to do anything to account
              // for column major order of lapack/blas routines
              pmatrixq<ntype,NT> om;
              if constexpr (NT < 0)
                {
                  om.allocate(N);
                }	  
              om = (*this);
              int *IPIV = new int[N];
              wrap_dgetrf(&(om.m[0][0]), IPIV, N);
              wrap_dgetri(&(om.m[0][0]), IPIV, N);
              delete[] IPIV;
              return om;
            }
          else if constexpr ( std::is_same<ntype, float>::value)
            {
              pmatrixq<ntype,NT> om;
              if constexpr (NT < 0)
                {
                  om.allocate(N);
                }	  
              om = (*this);
              int *IPIV = new int[N];
              wrap_fgetrf(&(om.m[0][0]), IPIV, N);
              wrap_fgetri(&(om.m[0][0]), IPIV, N);
              delete[] IPIV;
              return om;
            }
          else 
            {
              // gauss eliminatio is much faster!
#ifdef USE_GAUSS_INV
              return gaussj();
#else
              return invLUNR();
#endif
            }
        }
      else
        {
#ifdef USE_GAUSS_INV
          return gaussj();
#else
          return invLUNR();
#endif
        }
#else
#ifdef USE_GAUSS_INV
      return gaussj();
#else
      return invLUNR();
#endif
#endif
    }

  pmatrixq<ntype,NT> inv(void) const
    {
      // N.B. this first condition apparently is not needed
      // but in practice compiler optimizations get fooled nuts
      // if removed!
      if (N <= 4)
        {
          ntype invd;
          pmatrixq <ntype,NT> m1;
          //initialize dynamic class
          if constexpr (NT < 0)
            m1.allocate(N); 
          if (N==1)
            {
              m1[0][0]=1.0/m[0][0];
            }
          else if (N==2)
            {
#ifdef OPT_INV
              auto src = &(m[0][0]);
              auto dst = &(m1[0][0]);
              // 0 0 -> 0;  0 1 -> 1 
              // 1 0 -> 2;  1 1 -> 3 
              dst[0] = src[3];
              dst[1] = -src[1];
              dst[2] = -src[2];
              dst[3] = src[0];
              invd = src[0]*src[3] - src[1]*src[2];
              invd = 1.0/invd;
              dst[0] *= invd;
              dst[1] *= invd;
              dst[2] *= invd;
              dst[3] *= invd;
#else
              m1[0][0] = m[1][1];
              m1[0][1] = -m[0][1];
              m1[1][0] = -m[1][0];
              m1[1][1] = m[0][0];

              invd = m[0][0]*m[1][1];
              invd -= m[0][1]*m[1][0];
              invd = 1.0/invd;

              m1[0][0] *= invd;
              m1[0][1] *= invd;
              m1[1][0] *= invd;
              m1[1][1] *= invd;
#endif
            }
          else if (N==3)
            {
#ifdef OPT_INV
#if 0
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
              ntype m11m20= m[1][1]*m[2][0];
              ntype m12m20= m[1][2]*m[2][0];
              ntype m10m21= m[1][0]*m[2][1];
              ntype m12m21= m[1][2]*m[2][1];
              ntype m10m22= m[1][0]*m[2][2];
              ntype m11m22= m[1][1]*m[2][2];
              ntype A1 = -m[0][2]*m11m20 + m[0][1]*m12m20;
              ntype A2 = m[0][2]*m10m21 - m[0][0]*m12m21;
              ntype A3 = - m[0][1]*m10m22 + m[0][0]*m11m22;

              invd = A1 + A2 + A3;
              invd = 1.0/invd;
              //invd = 1.0/(-m[0][2]*m11m20 + m[0][1]*m12m20 + m[0][2]*m10m21
              //	  - m[0][0]*m12m21 - m[0][1]*m10m22 + m[0][0]*m11m22);

              m1[0][0] = -m12m21 + m11m22;
              m1[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
              m1[0][2] = -m[0][2]*m[1][1] + m[0][1]*m[1][2];
              m1[1][0] = m12m20 - m10m22;
              m1[1][1] = -m[0][2]*m[2][0] + m[0][0]*m[2][2]; 
              m1[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2]; 
              m1[2][0] = -m11m20 + m10m21; 
              m1[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
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
          else if (N==4)
            {
#ifdef OPT_INV
#if 0
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

              ntype A1 =   s0 * c5 - s1 * c4;
              ntype A2 =   s2 * c3 + s3 * c2; 
              ntype A3 = - s4 * c1 + s5 * c0; 

              invd = A1 + A2 + A3;
              //invd = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
              // Should check for 0 determinant
              invd = 1.0 / invd;
              ntype A4;

              A1 = m[1][1] * c5 - m[1][2] * c4;
              A2 = + m[1][3] * c3;
              A3 = -m[0][1] * c5 + m[0][2] * c4;
              A4 = - m[0][3] * c3;
              m1[0][0] = A1 + A2;
              m1[0][1] = A3 + A4;

              A1 = m[3][1] * s5 - m[3][2] * s4;
              A2 = m[3][3] * s3;
              m1[0][2] = A1 + A2;
              A3 = -m[2][1] * s5 + m[2][2] * s4;
              A4 = - m[2][3] * s3;
              m1[0][3] = A3 + A4;
              //m1[0][0] =  m[1][1] * c5 - m[1][2] * c4 + m[1][3] * c3;
              //m1[0][1] = -m[0][1] * c5 + m[0][2] * c4 - m[0][3] * c3;
              //m1[0][2] =  m[3][1] * s5 - m[3][2] * s4 + m[3][3] * s3;
              //m1[0][3] = -m[2][1] * s5 + m[2][2] * s4 - m[2][3] * s3;

              A1 =-m[1][0] * c5 + m[1][2] * c2;
              A2 = - m[1][3] * c1;
              A3 =m[0][0] * c5 - m[0][2] * c2;
              A4 = + m[0][3] * c1;
              m1[1][0] = A1 + A2;
              m1[1][1] = A3 + A4;

              A1 = -m[3][0] * s5 + m[3][2] * s2; 
              A2 = - m[3][3] * s1;
              m1[1][2] = A1 + A2;
              A3 =  m[2][0] * s5 - m[2][2] * s2;
              A4 = m[2][3] * s1;
              m1[1][3] = A3 + A4;

              //m1[1][0] = -m[1][0] * c5 + m[1][2] * c2 - m[1][3] * c1;
              //m1[1][1] =  m[0][0] * c5 - m[0][2] * c2 + m[0][3] * c1;
              //m1[1][2] = -m[3][0] * s5 + m[3][2] * s2 - m[3][3] * s1;
              //m1[1][3] =  m[2][0] * s5 - m[2][2] * s2 + m[2][3] * s1;

              A1 = m[1][0] * c4 - m[1][1] * c2;
              A2 = m[1][3] * c0;
              A3 =-m[0][0] * c4 + m[0][1] * c2;
              A4 = - m[0][3] * c0;
              m1[2][0] = A1 + A2;
              m1[2][1] = A3 + A4;

              A1 = m[3][0] * s4 - m[3][1] * s2;
              A2 = + m[3][3] * s0;
              m1[2][2] = A1 + A2;
              A3 = -m[2][0] * s4 + m[2][1] * s2;
              A4 = - m[2][3] * s0;
              m1[2][3] = A3 + A4;

              //m1[2][0] =  m[1][0] * c4 - m[1][1] * c2 + m[1][3] * c0;
              //m1[2][1] = -m[0][0] * c4 + m[0][1] * c2 - m[0][3] * c0;
              //m1[2][2] =  m[3][0] * s4 - m[3][1] * s2 + m[3][3] * s0;
              //m1[2][3] = -m[2][0] * s4 + m[2][1] * s2 - m[2][3] * s0;

              A1 = -m[1][0] * c3 + m[1][1] * c1;
              A2 = - m[1][2] * c0;
              A3 = m[0][0] * c3 - m[0][1] * c1;
              A4 = + m[0][2] * c0;
              m1[3][0] = A1 + A2;
              m1[3][1] = A3 + A4;

              A1 = -m[3][0] * s3 + m[3][1] * s1;
              A2 = -m[3][2] * s0;
              m1[3][2] = A1 + A2;
              A3 = m[2][0] * s3 - m[2][1] * s1;
              A4 = m[2][2] * s0;
              m1[3][3] = A3 + A4;

              //m1[3][0] = -m[1][0] * c3 + m[1][1] * c1 - m[1][2] * c0;
              //m1[3][1] =  m[0][0] * c3 - m[0][1] * c1 + m[0][2] * c0;
              //m1[3][2] = -m[3][0] * s3 + m[3][1] * s1 - m[3][2] * s0;
              //m1[3][3] =  m[2][0] * s3 - m[2][1] * s1 + m[2][2] * s0;

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

#else
              invd = 1.0/det();
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
          return m1;
        }
      else
        {
          auto m1=invgen(); 

          if (N <= NMAX)
            return m1;
          else
            return move(m1);
        }
    }
  inline pmatrixq<ntype,NT>& operator += (const pmatrixq<ntype,NT>& m1)
    {
      int i, j;
      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m[i][j] += m1.m[i][j];
      return (*this);
    }

  inline pmatrixq<ntype,NT>& operator -= (const pmatrixq<ntype,NT>& m1)
    {
      int i, j;
      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m[i][j] -= m1.m[i][j];
      return (*this);
    }
  // multiplication with scalars
  inline pmatrixq<ntype,NT>& operator *=(const double& param) 
    {
      int i,j;
      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m[i][j] *= param;
      return (*this);
    }
#ifdef MAT_USE_LAPACK
  pvector<ntype,NT> EigVal (void) 
    {
      pvector<ntype, NT> eval;
      int ok;
      if constexpr ( std::is_same<ntype, double>::value)
        wrap_dsyev(&(m[0][0]), NULL, &(eval[0]), N, &ok);
      else if constexpr ( std::is_same<ntype, double>::value)
        wrap_ssyev(&(m[0][0]), NULL, &(eval[0]), N, &ok);
      return eval;
    }
  pvector<ntype,NT> EigValVec(pmatrixq<ntype, NT>& m1) 
    {
      pvector<ntype, NT> eval;
      int ok;
      if constexpr ( std::is_same<ntype, double>::value)
        wrap_dsyev(&(m[0][0]), &(m1[0][0]), &(eval[0]), N, &ok, true);
      else if constexpr ( std::is_same<ntype, double>::value)
        wrap_ssyev(&(m[0][0]), &(m1[0][0]), &(eval[0]), N, &ok,true);
      return eval;
    }
#endif
  // il vettore restituito è quello con le soluzioni
  pvector<ntype,NT> SolveLineq (pvector<ntype,NT> x) 
    {
#ifdef MAT_USE_LAPACK
      int ok;
      if constexpr ( std::is_same<ntype, double>::value)
        {
          pvector<ntype,NT> sol;
          if constexpr (NT < 0)
            sol.allocate(N); 
          sol = x;
          wrap_dgesv(&(m[0][0]), &(sol[0]), N, &ok);
          return sol;
        }
      else if constexpr ( std::is_same<ntype, float>::value)
        {
          pvector<ntype,NT> sol;
          if constexpr (NT < 0)
            sol.allocate(N); 
          sol = x;
          wrap_sgesv(&(m[0][0]), &(sol[0]), N, &ok);
          return sol;
        }
      else
        {
          pmatrixq<ntype, NT> mt;
          ntype dd;
          pvector<int,NT> indx;
          if constexpr (NT < 0)
            {
              mt.allocate(N);
              indx.allocate(N);
            }
          mt = (*this).ludcmp(indx, &dd);
          return mt.lubksb(indx, x);
        }
#else
      pmatrixq<ntype, NT> mt;
      ntype dd;
      pvector<int, NT> indx;
      if constexpr (NT < 0)
        {
          mt.allocate(N);
          indx.allocate(N);
        }
      mt = (*this).ludcmp(indx, dd);
      //if (ok==1) // for now do not handle errors
      //return 1;
      return mt.lubksb(indx, x);
#endif
    }
#ifdef MAT_LAZY_EVAL
  // used if lazyness is disable (e.g. for small matrices)
  pmatrixq<ntype,NT> matscal(const ntype& param)
    {
      int i,j;
      pmatrixq<ntype,NT> m1;
      //initialize dynamic class
      if constexpr (NT < 0)
        m1.allocate(N); 

      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m1[i][j] = m[i][j]*param;
      if (N > NMAX)
        return std::move(m1);
      else 
        return m1;
    }
#endif
  inline pmatrixq<ntype,NT> operator *(const double& param) 
    {
      int i,j;
      pmatrixq<ntype,NT> m1;
      //initialize dynamic class
      if constexpr (NT < 0)
        m1.allocate(N); 

      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m1[i][j] = m[i][j]*param;
      return m1;
    }
  inline pmatrixq<ntype,NT>& operator /=(const double& param) 
    {
      int i,j;
      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m[i][j] /= param;
      return (*this);
    }
  inline pmatrixq<ntype,NT> operator /(const double& param) 
    {
      int i,j;
      pmatrixq<ntype,NT> m1;
      //initialize dynamic class
      if constexpr (NT < 0)
        m1.allocate(N); 

      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m1[i][j] = m[i][j]/param;
      return m1;
    }
#ifdef MAT_COMMA_INIT
  /* NOTE: comma operator has the lowest priority among operators 
   * hence it has less operators than <<.
   * In this case << operator initialize first element of matrix
   * and return matrix itself, therefore we need to overload comma operator
   * for pmatrixq object returned by << operator et voilà we have the comma initializer!  
   * */
  inline pmatrixq<ntype,NT>& operator,(const ntype& mr)
    {
      curidx++;
      if (curidx >= N*N)
        {
          cout << "Too many elements in comma initialization!\n";
          exit(-1);
        }
      m[curidx/N][curidx%N]=mr;
      return (*this);
    }
  inline pmatrixq<ntype,NT>& operator<<(const ntype& mr)
    {
      m[0][0] = mr;
      curidx=0;
      return (*this);
    } 
#endif
  friend pmatrixq<ntype,NT> operator *(const double& param, const pmatrixq<ntype,NT>& m1) 
    {
      int i, j;
      pmatrixq<ntype,NT> m2;
      if constexpr (NT < 0)
        {
          m2.allocate(N);
        }
      for (i=0; i < N; i++)
        for (j=0; j < N; j++)
          m2.m[i][j] = m1.m[i][j]*param;
      return m2;
    }
  // multiply two NxN matrices
  inline pmatrixq<ntype,NT>& operator *=(const pmatrixq<ntype,NT>& m1) 
    {
      *this = (*this).mul(m1);
      return *this;
    }
#ifndef LAZY_MUL
  inline pmatrixq<ntype,NT> operator *(const pmatrixq<ntype,NT>& m1) 
    {
      if (N > NMAX)
        return std::move((*this).mul(m1));// std:move is needed otherwise move constructor will not be used!!
      else
        return (*this).mul(m1);
    }
#endif
  inline pmatrixq<ntype,NT>& operator /=(const pmatrixq<ntype,NT>& m1) 
    {
      *this = (*this).mul(m1.inv());	
      return *this;
    }

  inline pmatrixq<ntype,NT> operator /(const pmatrixq<ntype,NT>& m1) 
    {
      if (N > NMAX)
        return std::move((*this).mul(m1.inv()));// std:move is needed otherwise move constructor will not be used!!	
      else 
        return (*this).mul(m1.inv());
    }

  // matrix times vector
  inline pvector<ntype,NT> operator *(const pvector<ntype,NT>& v1) 
    {
      int i, j;
      pvector <ntype,NT> v2;
      //initialize dynamic class
      if constexpr (NT < 0)
        v2.allocate(N); 

#ifdef MAT_USE_LAPACK
      if (N > NMAXMUL) 
        {
          if constexpr ( std::is_same<ntype, double>::value)
            wrap_dgemv('t', &(v2.v[0]), &(m[0][0]), &(v1.v[0]), N);
          else if constexpr ( std::is_same<ntype, float>::value)
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
      if (N > vecpars::NMAX)
        return std::move(v2);
      else
        return v2;
    }
  // transpose of v1 times m1
  friend pvector<ntype,NT> operator *(const pvector<ntype,NT> v1, const pmatrixq<ntype,NT>& m1) 
    {
      int i, j;
      pvector <ntype,NT> v2;
      if constexpr (NT < 0)
        v2.allocate(N); 

#ifdef MAT_USE_LAPACK
      if (N > NMAXMUL) 
        {
          if constexpr ( std::is_same<ntype, double>::value)
            wrap_dgemv('n', &(v2.v[0]), &(m[0][0]), &(v1.v[0]), N);
          else if constexpr ( std::is_same<ntype, float>::value)
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
      if (N > vecpars::NMAX)
        return std::move(v2);
      else
        return v2;
    }
  inline ntype* operator[](const int& i)
    {
      return &(m[i][0]);
    }

  inline void mulblock(pmatrixq<ntype,NT>& m2, const pmatrixq<ntype,NT>& m1)
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
  inline void mulstra(pmatrixq<ntype,NT>& m2, const pmatrixq<ntype,NT>& m1)
    {
      // strassen algorithm for matrix multiplication recursively defined (see Matrix Computations book)
      int i, j;
      using selftype  = typename std::conditional<N%2==1,pmatrixq<ntype,(N+1)/2>,pmatrixq<ntype,N/2>>::type; 
      static constexpr int dN = N % 2;
      static constexpr int Np = N+dN;
      selftype A11, A12, A21, A22, B11, B12, B22, B21,
               P1, P2, P3, P4, P5, P6, P7;
      if constexpr (NT < 0)
        {
          A11.allocate(Np);
          A12.allocate(Np);
          A21.allocate(Np);
          A22.allocate(Np);
          B11.allocate(Np);
          B12.allocate(Np);
          B22.allocate(Np);
          B21.allocate(Np);
          P1.allocate(Np);
          P2.allocate(Np);
          P3.allocate(Np);
          P4.allocate(Np);
          P5.allocate(Np);
          P6.allocate(Np);
          P7.allocate(Np);
        }
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

  void muladd(const pmatrixq<ntype,NT>& m1, const pmatrixq<ntype,NT>& m2, ntype alpha=1.0, ntype beta=1.0) 
    {
      // m=alpha*m*m1 + beta*m2
#ifdef MAT_USE_LAPACK
      if (N > NMAXMUL)
        {
          if constexpr ( std::is_same<ntype, double>::value)
            {
              wrap_dgemm(&(m[0][0]), &(m1.m[0][0]), &(m2.m[0][0]),NT, alpha, beta);
            }
          else if constexpr ( std::is_same<ntype, float>::value)
            {
              wrap_sgemm(&(m[0][0]), &(m1.m[0][0]), &(m2.m[0][0]),NT, alpha, beta);
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
  inline pmatrixq<ntype,NT> mul(const pmatrixq<ntype,NT>& m1)
    {
      int i, j, k;
      //this first condition apparently is not needed
      // but in practice compiler optimizations can get fooled 
      // if removed!
      if (N <= 4)
        {
          pmatrixq<ntype,NT> m2;
          //initialize dynamic class
          if constexpr (NT < 0)
            m2.allocate(N); 

          // for small matrices we use explicit expressions for multiplication
          if (N==1)
            {
              m2.m[0][0] = pmatrixq::m[0][0]*m1.m[0][0];
            }
          else if (N==2)
            {
#if 0
              auto src0 = &(pmatrixq::m[0][0]);
              auto src1 = &(m1.m[0][0]);
              auto dst = &(m2.m[0][0]);
              // 0 0 -> 0;  0 1 -> 1 
              // 1 0 -> 2;  1 1 -> 3 
              dst[0] = src0[0]*src1[0] + src0[1]*src1[2];
              dst[1] = src0[0]*src1[1] + src0[1]*src1[3];
              dst[2] = src0[2]*src1[0] + src0[3]*src1[2];
              dst[3] = src0[2]*src1[1] + src0[3]*src1[3]; 
#else
              m2.m[0][0] = pmatrixq::m[0][0]*m1.m[0][0] + pmatrixq::m[0][1]*m1.m[1][0];
              m2.m[0][1] = pmatrixq::m[0][0]*m1.m[0][1] + pmatrixq::m[0][1]*m1.m[1][1];
              m2.m[1][0] = pmatrixq::m[1][0]*m1.m[0][0] + pmatrixq::m[1][1]*m1.m[1][0];
              m2.m[1][1] = pmatrixq::m[1][0]*m1.m[0][1] + pmatrixq::m[1][1]*m1.m[1][1]; 
#endif
            }
          else if (N==3)
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
          else if (N==4)
            {
              //cout << pmatrix<ntype,NT,NT>::dynamic;
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
          return m2;
        }
      else
        {
          pmatrixq<ntype,NT> m2;
          //initialize dynamic class
          if constexpr (NT < 0)
            m2.allocate(N); 
#ifdef MAT_USE_LAPACK
          if (N > NMAXMUL)  
            {
              if constexpr ( std::is_same<ntype, double>::value)
                wrap_dgemm(&(m2.m[0][0]), &(m[0][0]), &(m1.m[0][0]), N);
              else if constexpr ( std::is_same<ntype, float>::value)
                wrap_sgemm(&(m2.m[0][0]), &(m[0][0]), &(m1.m[0][0]), N);
              else
                {
                  if (N > NMAXSTRA)
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
          if (N > NMAXSTRA) // 40 has been chose optimizing performance on my macbook pro 13" but can be architecture dependent
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
          if (N > NMAX)
            return std::move(m2);
          else
            return m2;
        }
    }
  // transpose
  pmatrixq<ntype,NT> transp()
    {
      int i, j;
      pmatrixq<ntype,NT> m1;
      //initialize dynamic class
      if constexpr (NT < 0)
        m1.allocate(N); 

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
  pvector<ntype,NT> col(int i)
    {
      pvector<ntype,NT> v;
      if constexpr (NT < 0)
        v.allocate(N); 

      int j;
      for (j=0; j < N; j++)
        {
          v[j] = m[j][i];
        }
      return v;
    }
  // get i-th row vector 
  pvector<ntype,NT> row(int i)
    {
      pvector<ntype,NT> v;
      int j;
      if constexpr (NT < 0)
        v.allocate(N); 
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
              cout << setprecision(std::numeric_limits<ntype>::digits10-1) << m[i][j];
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
};
#ifdef MAT_LAZY_EVAL
#ifdef LAZY_MUL
// MATRIX TIMES MATRIX
//
//
// AnOpM times matrix
  template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,NT> const& p)
{
    return AnOpM<ntype,NT, MOpTypes::MatTimesMat, AnOpM<ntype,NT, tipo, Lhs, Rhs>, pmatrixq<ntype,NT>>(lhs, p);
} 
// matrix times AnOpM 
template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(pmatrixq<ntype,NT> const& p, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpM< ntype,NT, MOpTypes::MatTimesMat, pmatrixq<ntype,NT>, AnOpM<ntype,NT, tipo, Lhs, Rhs> >(p, rhs);
}

// matrix times matrix
// se N < NLAZY restituisce pmatrixq e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpM e usa la lazy evaluation
  template <typename ntype, int NT> 
inline auto operator*(pmatrixq<ntype,NT> const& lhs, pmatrixq<ntype,NT> const& rhs) 
{
  if constexpr (NT != 2)
    {
      return lhs.mul(rhs); //normal evaluation
    }
  else
    {
      return AnOpM<ntype,NT, MOpTypes::MatTimesMat, pmatrixq<ntype,NT>, pmatrixq<ntype,NT>>(lhs, rhs); // lazy evaluation
    }
}
// AnOpM times AnOpM
template <typename ntype,int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator*(const AnOpM<ntype,NT, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpM<ntype,NT, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpM<ntype,NT, MOpTypes::MatTimesMat, AnOpM<ntype,NT, tipoL, LLhs, LRhs>, AnOpM<ntype,NT, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}

#endif
// ADDITION
// AnOpM plus matrix
  template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,NT> const& p)
{
  return AnOpM<ntype,NT, MOpTypes::MatPlusMat, AnOpM<ntype,NT, tipo, Lhs, Rhs>, pmatrixq<ntype,NT>>(lhs, p);
} 
// matrix plus AnOpM 
template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(pmatrixq<ntype,NT> const& p, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpM< ntype,NT, MOpTypes::MatPlusMat, pmatrixq<ntype,NT>, AnOpM<ntype,NT, tipo, Lhs, Rhs> >(p, rhs);
}

// matrix plus matrix
// se N < NLAZY restituisce pmatrixq e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpM e usa la lazy evaluation
  template <typename ntype, int NT> 
  // std::conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
 // typename std::conditional<((NT<=matpars::NLAZY&&NT >=0)||NT==matpars::DynNoLazy),pmatrixq<ntype,NT>,AnOpM<ntype,NT, MOpTypes::MatPlusMat, pmatrixq<ntype,NT>, pmatrixq<ntype,NT>>>::type
inline auto operator+(pmatrixq<ntype,NT> const& lhs, pmatrixq<ntype,NT> const& rhs) 
{
  if constexpr ((NT >= 0 && NT <= pmatrixq<ntype,NT>::NLAZY)|| NT==matpars::DynNoLazy)
    {
      // mettere qui direttamente l'addizione e provare se aumentano le performance
#if 0
      pmatrixq<ntype,NT> m2;
	  int i, j;
	  for (i=0; i < lhs.N; i++)
	    for (j=0; j < lhs.N; j++)
	      {
		m2.m[i][j] = lhs.m[i][j] + rhs.m[i][j];
	      }
	  return m2;
#else
	  return lhs.addition(rhs); //normal evaluation
#endif
    }
  else
    {
      return AnOpM<ntype,NT, MOpTypes::MatPlusMat, pmatrixq<ntype,NT>, pmatrixq<ntype,NT>>(lhs, rhs); // lazy evaluation
    }
}
// AnOpM plus AnOpM
template <typename ntype,int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator+(const AnOpM<ntype,NT, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpM<ntype,NT, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpM<ntype,NT, MOpTypes::MatPlusMat, AnOpM<ntype,NT, tipoL, LLhs, LRhs>, AnOpM<ntype,NT, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
// SUBTRACTION
// AnOpM minus matrix
  template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator-(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,NT> const& p)
{
  return AnOpM<ntype,NT, MOpTypes::MatMinusMat, AnOpM<ntype,NT, tipo, Lhs, Rhs>, pmatrixq<ntype,NT>>(lhs, p);
} 
// matrix minus AnOpM
template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator-(pmatrixq<ntype,NT> const& p, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpM< ntype,NT, MOpTypes::MatMinusMat, pmatrixq<ntype,NT>, AnOpM<ntype,NT, tipo, Lhs, Rhs> >(p, rhs);
}

// matrix plus matrix 
// se N < NLAZY restituisce pmatrixq e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpM e usa la lazy evaluation
  template <typename ntype, int NT> 
  // std::conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  //typename std::conditional <((NT<=matpars::NLAZY&&NT >=0)||NT==matpars::DynNoLazy), pmatrixq<ntype,NT>,AnOpM<ntype,NT, MOpTypes::MatMinusMat, pmatrixq<ntype,NT>, pmatrixq<ntype,NT>>>::type
inline auto operator-(pmatrixq<ntype,NT> const& lhs, pmatrixq<ntype,NT> const& rhs) 
{
  if constexpr ((NT <= matpars::NLAZY && NT >=0)|| NT==matpars::DynNoLazy)
    {
      return lhs.subtraction(rhs); //normal evaluation
    }
  else
    {
      return AnOpM<ntype,NT, MOpTypes::MatMinusMat, pmatrixq<ntype,NT>, pmatrixq<ntype,NT>>(lhs, rhs); // lazy evaluation
    }
}
// AnOpM plus AnOpM
  template <typename ntype,int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator-(const AnOpM<ntype,NT, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpM<ntype,NT, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpM<ntype,NT, MOpTypes::MatMinusMat, AnOpM<ntype,NT, tipoL, LLhs, LRhs>, AnOpM<ntype,NT, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
// DIVISION
// AnOpM divided by AnOpM
 template <typename ntype,int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator/(AnOpM<ntype,NT, tipoL, LLhs, LRhs> const& lhs, AnOpM<ntype,NT, tipoR, RLhs, RRhs> const& rhs)
{
  pmatrixq<ntype,NT> mL, mR;
  int i, j;
  if constexpr (NT < 0)
    {
      int N=rhs.get_N();
      mL.allocate(N);
      mR.allocate(N);
    }
 
  for (i=0; i < mL.N; i++)
    for (j=0; j < mL.N; j++)
      {
	mL[i][j] = lhs.get_m(i,j);
	mR[i][j] = rhs.get_m(i,j);
      }
  return mL/mR;
}
// AnOpM divided by matrix
 template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator/(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,NT> const& mR)
{
  pmatrixq<ntype,NT> mL;
  int i, j;
  if constexpr (NT < 0)
    {
      mL.allocate(mR.N);
    }

  for (i=0; i < mR.N; i++)
    for (j=0; j < mR.N; j++)
      {
	mL[i][j] = lhs.get_m(i,j);
      }
  return mL*mR;
}
// matrix divided by AnOpM 
 template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator/(pmatrixq<ntype,NT> const& mL, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs)
{
  pmatrixq<ntype,NT> mR;
  int i, j;
  if constexpr (NT < 0)
    {
      mR.allocate(mL.N);
    }
  for (i=0; i < mL.N; i++)
    for (j=0; j < mL.N; j++)
      {
	mL[i][j] = rhs.get_m(i,j);
      }
  return mL/mR;
}
///
// MULTIPLICATION 
//
// AnOpM times AnOpM
#ifndef LAZY_MUL
  template <typename ntype,int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator*(const AnOpM<ntype,NT, tipoL, LLhs, LRhs> &lhs, const AnOpM<ntype,NT, tipoR, RLhs, RRhs> &rhs)
{
  pmatrixq<ntype,NT> mL, mR;
  if constexpr (NT < 0)
    {
      int N=lhs.get_N();
      mL.allocate(N);
      mR.allocate(N);
    }
 
  for (int i=0; i < mL.N; i++)
    for (int j=0; j < mL.N; j++)
      {
    	mR[i][j] = lhs.get_m(i,j);
	mL[i][j] = rhs.get_m(i,j);
      }
  return mR*mL;
}
//AnOpM times matrix
   template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, pmatrixq<ntype,NT> const& mR)
{
  pmatrixq<ntype,NT> mL;
  if constexpr (NT < 0)
    {
      mL.allocate(mR.N);
    }
 
  for (int i=0; i < mR.N; i++)
    for (int j=0; j < mR.N; j++)
      {
    	mL[i][j] = lhs.get_m(i,j);
      }
  return mL*mR;
} 
// matrix times AnOpM
template<typename ntype,int NT, int tipo, typename Lhs, typename Rhs>  
inline auto operator*(pmatrixq<ntype,NT> const& mL, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs) 
{
  pmatrixq<ntype,NT> mR;
  if constexpr (NT < 0)
    {
      mR.allocate(mL.N);
    }
  
  for (int i=0; i < mL.N; i++)
    for (int j=0; j < mL.N; j++)
      {
    	mR[i][j] = rhs.get_m(i,j);
      }
  return mL*mR;
}
#endif
// vector times AnOpM
template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator *(const pvector<ntype,NT> vL, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs) 
{
  pmatrixq<ntype,NT> mR;
  if constexpr (NT < 0)
    {
      mR.allocate(rhs.get_N());
    }
  
  for (int i=0; i < mR.N; i++)
    for (int j=0; j < mR.N; j++)
      {
    	mR[i][j] = rhs.get_m(i,j);
      }
  return vL*mR;
}
// AnOpM times vector
template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator *(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, const pvector<ntype,NT> vR) 
{
  pmatrixq<ntype,NT> mL;
  if constexpr (NT < 0)
    {
      mL.allocate(lhs.get_N());
    }
  for (int i=0; i < mL.N; i++)
    for (int j=0; j < mL.N; j++)
      {
    	mL[i][j] = lhs.get_m(i,j);
      }
  return mL*vR;
}
// scalar times matrix
  template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator *(const ntype& lhs, pmatrixq<ntype,NT> const& rhs) 
{
  if constexpr (NT <= matpars::NLAZY&&NT>=0)
    return rhs.matscal(lhs);
  else
    return AnOpM<ntype,NT, MOpTypes::ScalTimesMat, ntype, pmatrixq<ntype,NT>>(lhs,rhs);
}

// scalar times AnOpM 
template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
 inline auto operator *(const ntype& lhs, AnOpM<ntype,NT, tipo, Lhs, Rhs> const& rhs) 
    {
      return AnOpM<ntype,NT,MOpTypes::ScalTimesMat,ntype,AnOpM<ntype,NT,tipo,Lhs,Rhs>>(lhs,rhs);
    }
// matrix times scalar
  template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator *(pmatrixq<ntype,NT> const& lhs, const ntype& rhs) 
{
  if constexpr ((NT <= matpars::NLAZY&&NT>=0)||NT==matpars::DynNoLazy)
    return lhs.matscal(rhs);
  else
    return AnOpM<ntype,NT, MOpTypes::MatTimesScal, pmatrixq<ntype,NT>, ntype>(lhs,rhs);
}
// AnOpM times scalar
  template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator *(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, const ntype& rhs) 
{
  return AnOpM<ntype,NT,MOpTypes::MatTimesScal,AnOpM<ntype,NT,tipo,Lhs,Rhs>, ntype>(lhs,rhs);
}

// matrix divided by scalar
  template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator /(pmatrixq<ntype,NT> const& lhs, const ntype& rhs) 
{
  if constexpr ((NT <= matpars::NLAZY && NT >=0)||NT==matpars::DynNoLazy)
    return lhs.matscal(1.0/rhs);
  else
    return AnOpM<ntype,NT, MOpTypes::MatDivScal, pmatrixq<ntype,NT>, ntype>(lhs,rhs);
}
// AnOpM divided by scalar
  template <typename ntype,int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator /(AnOpM<ntype,NT, tipo, Lhs, Rhs> const& lhs, const ntype& rhs) 
{
  return AnOpM<ntype,NT,MOpTypes::MatDivScal,AnOpM<ntype,NT,tipo,Lhs,Rhs>, ntype>(lhs,rhs);
}
#endif

typedef pmatrixq<double,3> pmatrix3d;
typedef pmatrixq<int,3>    pmatrix3i;
#endif
