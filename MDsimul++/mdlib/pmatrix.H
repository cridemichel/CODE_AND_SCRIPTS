// matrici NxM
#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include<cmath>
#include "pvector.H"
 template <class ntype, int N, int M> class pmatrix {
    int nr, nc;
public:
    ntype m[N][M]; // N riga M colonna
#if 0 // example of operator = overload <-- not needed anyway!!!
    pmatrix<ntype,N,M>& operator=(const pmatrix<ntype,N,M>& m1)
      {
	int i,j;
	for (i=0; i < N; i++)
	  for (j=0; j < N; j++)
	    m[i][j] = m1.m[i][j];
	return (*this);
      }
#endif
    pmatrix()
      {
	nr=N;
	nc=M;
      }
    ~pmatrix()
      {

      }
};	
// redefine also *= += and
// define inv() det() i.e. inverse and determinant
// matrici quadrate NxN

//using namespace std;
template <class ntype> double abs(ntype a)
{
  return abs(a);
}

template <class ntype, int N> class pmatrixq: public pmatrix <ntype, N, N> 
{
  using pmatrix<ntype,N,N>::m;
  int n;
  pmatrixq<ntype,N> ludcmp(pvector<int,N>& indx, double& d)
    {
      /* return the LU decomposition of calling matrix (i.e. *this) */
      pmatrixq <ntype,N> m1;
      const double TINY=1E-20;
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

  pvector<ntype,N> lubksb(pvector <int,N> indx, pvector<ntype,N> b)
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
      for (i=n-1;i>=0;i--) 
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
      static pmatrixq<ntype, N> m2;
      int i, j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m2[i][j] = m[i][j] +  m1[i][j];
      return m2;
    }
  double detgen()
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
  double det()
    {
      switch ((*this).size())
	{
	case 1:
	  return m[0][0];
	  break;
	case 2:
	  return m[0][0]*m[1][1]-m[0][1]*m[1][0];
	  break;
	case 3:
	  return -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + 
	    m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - 
	    m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2]; 
	  break;
	default:
	  return detgen();
	  break;
	}
    }

  pmatrixq<ntype,N> invgen()
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
  pmatrixq<ntype,N> inv()
    {
      double invd;
      pmatrixq <ntype,N> m1;
      switch ((*this).size())
	{
	case 1: 
	  m1[0][0]=1.0/m[0][0];
	  break; 
	case 2:
	  invd = 1.0/det();
	  m1[0][0] = m[1][1]*invd;
	  m1[0][1] = m[1][0]*invd;
	  m1[1][0] = m[0][1]*invd;
	  m1[1][1] = m[0][0]*invd;
	  break;
	case 3:
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
	  break;
	default:
	  m1=invgen(); 
	  break;
	}
      return m1;
    }
  inline pmatrixq<ntype, N> operator += (const pmatrixq<ntype, N>& m1)
    {
      int i, j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m[i][j] += m1.m[i][j];
      return (*this);
    }

  inline pmatrixq<ntype, N> operator *(const double& param) 
    {
      int i,j;
      static pmatrixq<ntype,N> m1;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m1[i][j] = m[i][j]*param;
      return m1;
    }

  friend pmatrixq<ntype,N> operator *(const double& param, const pmatrixq<ntype,N>& m1) 
    {
      int i, j;
      static pmatrixq<ntype,N> m2;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m2.m[i][j] = m1.m[i][j]*param;
      return m2;
    }
  // multiply two NxN matrices
  inline pmatrixq<ntype,N> operator *(const pmatrixq<ntype,N>& m1) 
    {
      return (*this).mul(m1);	
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
  inline double* operator[](const int& i)
    {
      return &(m[i][0]);
    }
  inline pmatrixq<ntype, N> mul(const pmatrixq<ntype,N>& m1)
    {
      int i, j, j2;
      pmatrixq<ntype,N> m2;
      for (i=0; i < N; i++)
	{
	  for (j=0; j < N; j++)
	    {
	      m2[i][j]=0.0;
	      for (j2=0; j2 < N; j2++)
		m2[i][j] += pmatrixq::m[i][j2]*m1.m[j2][j];
	    }
	}
      return m2;
    }
  void show(const char* str)
    {
      int i, j;
      if (str!=NULL)
	printf("%s", str);
      printf("{");
      for (i=0; i < N; i++)
	{
	  printf("{");
	  for (j=0; j < N; j++)
	    { 
	      printf("%.15G", m[i][j]);
		if (j < N-1)
		  printf(",");
	      }
	    printf("}");
	    if (i < N-1)
	      printf(",\n");
	  }
	printf("}\n");
      }
    void show(void)
      {
	show(NULL);
      }

    int size()
      {
	return N;
      } 
    pmatrixq()
      {
	n=N;
      }
   
    ~pmatrixq()
      {
      }
};

typedef pmatrixq<double,3> pmatrix3d;
typedef pmatrixq<int,3>    pmatrix3i;
