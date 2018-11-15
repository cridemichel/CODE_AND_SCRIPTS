#ifndef _PMATRIX_
/* N.B. fino a N=10 vanno meglio di GLM e armadillo in tutte le operazioni (somma, moltiplicazione e inversione). */
#define _PMATRIX_
// matrici NxM
#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include<cmath>
#include "pvector.H"
#define OPT_INV
//#define MOVE_SEMANTIC
template <class ntype, int N, int M> class pmatrix {
  //int nr, nc;
public:
#ifndef MOVE_SEMANTIC 
  ntype m[N][M];
#else
  ntype (*m)[M]; // N riga M colonna
  // example of operator = overload <-- not needed anyway!!! 
  pmatrix<ntype,N,M>& operator=(const pmatrix<ntype,N,M>& m1)
    {
      int i,j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] = m1.m[i][j];
      return (*this);
    }
  pmatrix<ntype,N,M>& operator=(pmatrix<ntype,N,M>&& m1)
    {
      auto mt = (*this).m;
      (*this).m=m1.m;
      m1.m = mt;
      return (*this);
    }
  pmatrix(pmatrix<ntype,N,M>&& m1)
    {
      //nr=N;
      //nc=M;
      (*this).m = m1.m;
      m1.m=nullptr;
    }
  pmatrix()
    {
      //nr=N;
      //nc=M;
      m = new ntype[N][M];
    }
  ~pmatrix()
    {
      delete[] m;
    }
#endif
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
template <class ntype, int N> class pmatrixq: public pmatrix <ntype, N, N> 
{
  using pmatrix<ntype,N,N>::m;
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
  inline pmatrixq<ntype,N> operator +(const pmatrixq<ntype,N>& m1)
    {
      pmatrixq<ntype, N> m2;
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
  inline pmatrixq<ntype,N> operator -(const pmatrixq<ntype,N>& m1)
    {
      pmatrixq<ntype, N> m2;
      int i, j;
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


  inline double detgen() const
    {
      double d;
      pmatrixq<ntype,N> m1;
      int j;
      pvector<int,N> indx;
      // ludcmp return the LU decomposition of actual matrix (*this)
      m1=ludcmp(indx,d);
      for(j=0;j<N;j++) 
	d *= m1[j][j];
      return d;
    }
  inline double det() const
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
      int m1, m2; 
      double d;
      pmatrixq<ntype,N> om, LU;
      pvector <ntype,N> col, col2;
      pvector<int,N> indx;
      // ludcmp returns the LU decomposition of actual matrix (*this)
      // where indx will store the raw permutiation effected by partial pivoting (see Numerical Recipe
      // for more details c2.3)
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
      return (*this).mul(m1);	
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
      return (*this).mul(m1.inv());	
    }
  
  // matrix times vector
  inline pvector<ntype,N> operator *(const pvector<ntype,N>& v1) 
    {
      int i, j;
      pvector <ntype,N> v2;
      for (i=0; i < N; i++)
	{
	  v2.v[i]=0;
	  for (j=0; j < N; j++)
	    v2.v[i] += m[i][j]*v1.v[j];
	}
      return v2;
    }
  // transpose of v1 times m1
  friend pvector<ntype,N> operator *(const pvector<ntype,N> v1, const pmatrixq<ntype,N>& m1) 
    {
      int i, j;
      pvector <ntype,N> v2;
      for (i=0; i < N; i++)
	{
	  v2.v[i]=0;
	  for (j=0; j < N; j++)
	    v2.v[i] += v1.v[j]*m1.m[j][i];
	}
      return v2;
    }
  inline ntype* operator[](const int& i)
    {
      return &(m[i][0]);
    }
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
      else
	{
	  if constexpr (N > 10)
	    {
#if 1
	      using selftype  = typename std::conditional<N%2==1,pmatrixq<ntype,(N+1)/2>,pmatrixq<ntype,N/2>>::type; 
	      selftype A11, A12, A21, A22, B11, B12, B22, B21,
		  P1, P2, P3, P4, P5, P6, P7;
	      static constexpr int dN = N % 2;
	      static constexpr int Np = N+dN;
	      //if constexpr (dN==1)
		//cout << "dN=" << dN << "\n";
	      // Strassen method applied by recursion
	      // we pad with 0s if matrix is odd
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
#else
	      /* tiled matrix multiplication to optimize cache usage (from wikipedia!) 
	       * [ Strassen is much better though! ]*/
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

typedef pmatrixq<double,3> pmatrix3d;
typedef pmatrixq<int,3>    pmatrix3i;
#endif
