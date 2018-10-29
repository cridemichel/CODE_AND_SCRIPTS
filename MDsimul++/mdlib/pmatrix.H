// matrici NxM
#include<stdio.h>
#include<stdlib.h>
#include "pvector.H"
template <class ntype, int N, int M> class pmatrix {
    int nr, nc;
public:
    ntype m[N][M]; // N riga M colonna
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
template <class ntype, int N> class pmatrixq: public pmatrix <ntype, N, N> 
      {
	using pmatrix<ntype,N,N>::m;
	int n;
public:
    inline pmatrixq<ntype,N> operator +(const pmatrixq<ntype,N>& m1)
    {
      static pmatrixq<ntype, N> m2;
      int i, j;
      for (i=0; i < N; i++)
	for (j=0; j < N; j++)
	  m2.m[i][j] = m[i][j] +  m1.m[i][j];
      return m2;
    }
    double det()
      {
	switch (N)
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
	  }
      }
    pmatrixq inv()
      {
	double invd;
	switch (N)
	  {
	  case 1: 
	    return 1.0/pmatrixq::m[0][0];
	    break; 
	  case 2:
	    invd = 1.0/det();
	    m[0][0] = m[1][1]*invd;
	    m[0][1] = m[1][0]*invd;
	    m[1][0] = m[0][1]*invd;
	    m[1][1] = m[0][0]*invd;
	    break;
	  case 3:
	    invd = 1.0/det();
	    m[0][0] = -m[1][2]*m[2][1] + m[1][1]*m[2][2];
	    m[0][1] =  m[0][2]*m[2][1] - m[0][1]*m[2][2];
	    m[0][2] = -m[0][2]*m[1][1] + m[0][1]*m[1][2];
	    m[1][0] =  m[1][2]*m[2][0] - m[1][0]*m[2][2];
	    m[1][1] = -m[0][2]*m[2][0] + m[0][0]*m[2][2]; 
	    m[1][2] =  m[0][2]*m[1][0] - m[0][0]*m[1][2]; 
	    m[2][0] = -m[1][1]*m[2][0] + m[1][0]*m[2][1]; 
	    m[2][1] =  m[0][1]*m[2][0] - m[0][0]*m[2][1];
	    m[2][2] = -m[0][1]*m[1][0] + m[0][0]*m[1][1];
	    m[0][0] *= invd;
	    m[0][1] *= invd;
	    m[0][2] *= invd;
	    m[1][0] *= invd; 
	    m[1][1] *= invd;
	    m[1][2] *= invd;
	    m[2][0] *= invd;
	    m[2][1] *= invd; 
	    m[2][2] *= invd; 	    
	    break;
	  }
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
