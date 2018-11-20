#ifndef _PVECTOR_
#define _PVECTOR_
#include <limits>
#include<cstdlib>
#include<iostream>
#include<iomanip>
using namespace std;
#define MOVE_SEMANTIC
#define LAZY_EVAL
#ifdef LAZY_EVAL
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>
struct PmOpV {
    Lhs const& lhs;
    Rhs const& rhs;
    PmOpV(Lhs const& lhs, Rhs const& rhs):lhs(lhs), rhs(rhs) {
      // empty body
      //cout << "qui\n";
    }
    inline ntype get_v(int i) const {
      if constexpr (tipo == 0)// vec + vec
	return lhs.get_v(i) + rhs.get_v(i);
      else if constexpr (tipo == 1) // vec - vec
	return lhs.get_v(i) - rhs.get_v(i);
      else if constexpr (tipo == 2) // vec * scalar
	return lhs.get_v(i)*rhs;
      else if constexpr (tipo == 3) // scalar * vec NOT IMPLEMENTED YET
	return rhs.get_v(i)*lhs;
      else if constexpr (tipo == 4) // vec * vec NOT IMPLEMENTED YET
	return rhs.get_v(i)*lhs;
      return 0;
    }
    Lhs const& get_lhs() const { return lhs; }
    Rhs const& get_rhs() const { return rhs; }
};
#endif
//NMAX_=45, int MAXSTRA_=40, int NMAXINV_=40, int NMAXMUL_=10, int NSTA_=8, int NADD_=5
struct vecpars
{
  // params for my class
  static const int NMAX = 45; // if N>NMAX use dynamic allocation of matrices
  static const int NMAXSTRA = 40;
  static const int NMAXINV = 40; // if N>NAMXINV use LAPACK LU DECOMP
  static const int NMAXMUL = 10; // if N>NMAXMUL use BLAS multiplication routines
  static const int NSTA = 8;// value below which declares static object in operator funcs (only if m is not dynamically allocated)
  static const int NADD = 5;// value below which does not use lazy evaluation for addition
};

// template for a static vector
template <class ntype, int N> class pvecsta {
public:
  ntype v[N];
  constexpr static int dynamic = false;
  pvecsta<ntype,N>& operator=(const pvecsta<ntype,N>& v1)
    {
      int i;
      for (i=0; i < N; i++)
	  v[i] = v1.v[i];
      return (*this);
    }};
// template for a dynamic vector (i.e. allocated with new)
template <class ntype, int N> class pvecdyn {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  ntype (*v); // N elements
  // copy assignment
  pvecdyn<ntype,N>& operator=(const pvecdyn<ntype,N>& v1)
    {
      int i;
      for (i=0; i < N; i++)
	  v[i] = v1.v[i];
      return (*this);
    }
#ifdef MOVE_SEMANTIC
  // move assignment
  pvecdyn<ntype,N>& operator=(pvecdyn<ntype,N>&& v1)
    {
      swap(v, v1.v);
      return (*this);
    }
#endif
  // copy constructor 
  pvecdyn(const pvecdyn<ntype,N>& v1)
    {
      //nr=N;
      //nc=M;
      v = new ntype[N];
      (*this) = v1;
    }
#ifdef MOVE_SEMANTIC
  // move constructor
  pvecdyn(pvecdyn<ntype,N>&& v1)
    {
      (*this).v = v1.v;
      v1.v=nullptr;
    }
#endif
  // default constructor
  pvecdyn()
    {
      v = new ntype[N];
    }
  // destructor
  ~pvecdyn()
    {
      delete[] v;
    }
};	
template <class ntype, int N> using pvecbase = typename std::conditional<(N>vecpars::NMAX),pvecdyn <ntype, N>,pvecsta <ntype, N>>::type;
template <class ntype, int N> class pvector : public pvecbase<ntype,N>, vecpars {
  using pvecbase<ntype,N>::v;
  //int vsize;
  static const int maxdigits=numeric_limits<ntype>::digits10;
public:
#ifdef LAZY_EVAL
  template<typename Lhs, typename Rhs, int tipo>
    pvector<ntype,N>& operator=(PmOpV<ntype,N, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	v[i] = op.get_v(i);
      return (*this);
    }
  template<typename Lhs, typename Rhs, int tipo>
    pvector<ntype,N>(PmOpV<ntype,N, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	v[i] = op.get_v(i);
    }
  pvector<ntype,N>() = default;
  inline ntype get_v(int i) const {return v[i];}
#endif


#ifndef LAZY_EVAL
  //ntype v[N];
  inline pvector<ntype,N> operator +(const pvector<ntype,N>& param)
    {
      return addition(param);
    }
  inline pvector<ntype,N> operator -(const pvector<ntype,N>& param)
    {
      return subtraction(param);
    }
#endif
  // addition
  inline pvector<ntype,N> addition(const pvector<ntype,N>& v1) const
    {
      if constexpr (N <= NSTA && pvector<ntype,N>::dynamic==false)
	{
	  static pvector<ntype, N> v2;
	  int i;
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] + v1.v[i];
	  return v2;
	}
      else
	{
	  pvector<ntype, N> v2;
	  int i;
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] + v1.v[i];
	  return std::move(v2);
	}
    }
  inline pvector<ntype,N> subtraction(const pvector<ntype,N>& v1) const
    {
      int i;
      if constexpr (N <= NSTA && pvector<ntype,N>::dynamic==false)
	{
	  static pvector<ntype, N> v2;
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] - v1.v[i];
	  return v2;
	}
      else
	{
	  pvector<ntype, N> v2;
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] - v1.v[i];
	  return std::move(v2);
	}
    }

  inline pvector<ntype,N> operator += (const pvector<ntype,N>& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] += param.v[i];
      return *this;
    }

  // vec times scalar
#ifndef LAZY_EVAL
  inline pvector<ntype,N> operator *(const double& param) 
    {
      int i;
      pvector<ntype,N> lv;
      for (i=0; i < N; i++)
	lv.v[i] = v[i]*param;
      return std::move(lv);
    }
#endif
  inline pvector<ntype,N> vecscal(const double& param) const
    {
      int i;
      pvector<ntype,N> lv;
      for (i=0; i < N; i++)
	lv.v[i] = v[i]*param;
      return std::move(lv);
    }
  friend pvector<ntype,N> operator *(const ntype& param, const pvector<ntype,N>& v) 
    {
      int i;
      pvector<ntype,N> lv;
      for (i=0; i < N; i++)
	lv.v[i] = v.v[i]*param;
      return std::move(lv);
    }

  inline ntype& operator[](const int& i)
    {
      return v[i];
    }
  inline ntype dot(const pvector<ntype,N>& v1)
    {
      int i;
      ntype sum=0;
      for (i=0; i < N; i++)
	{
	  sum+=v[i]*v1.v[i]; 
	}
      return sum;
    }
  inline double norm()
    {
      return (*this).dot(*this);
    }
  // cross product
  inline pvector<ntype,N> cross(const pvector<ntype,N>& v1)
    {
      static pvector<ntype,N> v2;
      v2.v[0] = v[1]*v1.v[2] - v[2]*v1.v[1];
      v2.v[1] = v[2]*v1.v[0] - v[0]*v1.v[2];
      v2.v[2] = v[0]*v1.v[1] - v[1]*v1.v[0];
      return v2;
    }
  //^ = cross product 
  inline pvector<ntype,N> operator ^(const pvector<ntype,N> &vv)
    {
      return (*this).cross(vv);
    }
  //* dot product
  inline double operator *(const pvector<ntype,N> &vv)
    {
      return (*this).dot(vv);
    }
  //friend double scalprod(const pvector<ntype,N>& v1, const pvector<ntype,N>& v2);
  void show(void)
    {
      show(NULL);
    }
  void show(const char* str)
    {
      int i;
      if (str!=NULL)
	cout << str;
      cout << "{";
      for (i=0; i < N; i++)
	{
	  cout << setprecision(maxdigits) << v[i];
	  if (i < N-1)
	    cout << ",";
	}
      cout << "}\n";
    }
  int size()
    {
      return N;
    } 
};
template <class ntype, int N>
double scalprod(const pvector<ntype,N>& v1, const pvector<ntype,N>& v2)
{
  int i;
  double sum=0;
  for (i=0; i < N; i++)
    {
      sum+=v1.v[i]*v2.v[i]; 
    }
  return sum;
}
#ifdef LAZY_EVAL
// ADDITION
// add expression template with point at the right
  template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> PmOpV<ntype,N, 0, PmOpV<ntype,N, tipo, Lhs, Rhs>, pvector<ntype,N>> 
inline operator+(PmOpV<ntype, N, tipo, Lhs, Rhs> const& lhs, pvector<ntype,N> const& p)
{
  return PmOpV<ntype, N, 0, PmOpV<ntype, N, tipo, Lhs, Rhs>, pvector<ntype,N>>(lhs, p);
} 
// add expression template with point at the left
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> PmOpV< ntype, N, 0, pvector<ntype,N>, 
  PmOpV<ntype, N, tipo, Lhs, Rhs> > 
inline operator+(pvector<ntype,N> const& p, PmOpV<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  return PmOpV< ntype, N, 0, pvector<ntype,N>, PmOpV<ntype, N, tipo, Lhs, Rhs> >(p, rhs);
}

// add two points, yield a expression template    
// se N < NSTA restituisce pvector e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce PmOpV e usa la lazy evaluation

  template <typename ntype, int N> 
  // conditional restituisce il giusto tipo in base alla condizione (N <= NADD)
  typename std::conditional<(N<=vecpars::NADD),pvector<ntype,N>,PmOpV<ntype, N, 0, pvector<ntype,N>, pvector<ntype,N>>>::type
inline operator+(pvector<ntype,N> const& lhs, pvector<ntype,N> const& rhs) 
{
  if constexpr (N <= vecpars::NADD)
    {
      return lhs.addition(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return PmOpV<ntype, N, 0, pvector<ntype,N>, pvector<ntype,N>>(lhs, rhs); // lazy evaluation
    }
}
// caso aggiunto secondo un post in https://stackoverflow.com/questions/414243/lazy-evaluation-in-c
template <typename ntype, int N, int tipo, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline PmOpV<ntype, N, 0, PmOpV<ntype, N, tipo, LLhs, LRhs>, PmOpV<ntype, N, tipo, RLhs, RRhs>> operator+(const PmOpV<ntype, N, tipo, LLhs, LRhs> & leftOperandconst, const PmOpV<ntype, N, tipo, RLhs, RRhs> & rightOperand)
{
  return  PmOpV<ntype, N, 0, PmOpV<ntype, N, tipo, LLhs, LRhs>, PmOpV<ntype, N, tipo, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
// SUBTRACTION
// add expression template with point at the right
  template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> PmOpV<ntype,N, 1, PmOpV<ntype,N, tipo, Lhs, Rhs>, pvector<ntype,N>> 
inline operator-(PmOpV<ntype,N, tipo, Lhs, Rhs> const& lhs, pvector<ntype,N> const& p)
{
  return PmOpV<ntype, N, 1, PmOpV<ntype, N, tipo, Lhs, Rhs>, pvector<ntype,N>>(lhs, p);
} 
// add expression template with point at the left
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> PmOpV< ntype, N, 1, pvector<ntype,N>, 
  PmOpV<ntype, N, tipo, Lhs, Rhs> > 
inline operator-(pvector<ntype,N> const& p, PmOpV<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  return PmOpV< ntype, N, 1, pvector<ntype,N>, PmOpV<ntype, N, tipo, Lhs, Rhs> >(p, rhs);
}

// add two points, yield a expression template    
// se N < NSTA restituisce pvector e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce PmOpV e usa la lazy evaluation
template <typename ntype, int N> 
// conditional restituisce il giusto tipo in base alla condizione (N <= NADD)
  typename std::conditional<(N<=vecpars::NADD),pvector<ntype,N>,PmOpV<ntype, N, 1, pvector<ntype,N>, pvector<ntype,N>>>::type
inline operator-(pvector<ntype,N> const& lhs, pvector<ntype,N> const& rhs) 
{
  if constexpr (N <= vecpars::NADD)
    {
      return lhs.subtraction(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << NSTA << "\n";
      return PmOpV<ntype, N, 1, pvector<ntype,N>, pvector<ntype,N>>(lhs, rhs); // lazy evaluation
    }
}
// caso aggiunto secondo un post in https://stackoverflow.com/questions/414243/lazy-evaluation-in-c
  template <typename ntype, int N, int tipo, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline PmOpV<ntype, N, 1, PmOpV<ntype, N, tipo, LLhs, LRhs>, PmOpV<ntype, N, tipo, RLhs, RRhs>> operator-(const PmOpV<ntype, N, tipo, LLhs, LRhs> & leftOperandconst, const PmOpV<ntype, N, tipo, RLhs, RRhs> & rightOperand)
{
  return  PmOpV<ntype, N, 1, PmOpV<ntype, N, tipo, LLhs, LRhs>, PmOpV<ntype, N, tipo, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
   
///
// Scalar product of two PmOpV
//
  template <typename ntype, int N, int tipo, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline ntype operator*(const PmOpV<ntype, N, tipo, LLhs, LRhs> & leftOperand, const PmOpV<ntype, N, tipo, RLhs, RRhs> & rightOperand)
{
  pvector<ntype, N> vR, vL;
  for (int i=0; i < N; i++)
    {
      vL[i] = leftOperand.get_v(i);
      vR[i] = rightOperand.get_v(i);
    }
  return vL*vR;
}
// PmOpV times pvector
//
// left PmOpV
   template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> ntype 
inline operator*(PmOpV<ntype,N, tipo, Lhs, Rhs> const& lhs, pvector<ntype,N> const& vR)
{
  pvector<ntype, N> vL;
  for (int i=0; i < N; i++)
    vL[i] = lhs.get_v(i);
  return vL*vR;
} 
// right PmOpV
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>  
inline ntype operator*(pvector<ntype,N> const& vL, PmOpV<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  pvector<ntype, N> vR;
  for (int i=0; i < N; i++)
    {
      vR[i] = rhs.get_v(i);
    }
  return vL*vR;
}
// PmOvV times scalar
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> PmOpV<ntype,N, 2, PmOpV<ntype,N, tipo, Lhs, Rhs>, ntype> 
inline operator*(PmOpV<ntype, N, tipo, Lhs, Rhs> const& lhs, const ntype& p)
{
  return PmOpV<ntype, N, 2, PmOpV<ntype, N, tipo, Lhs, Rhs>, ntype>(lhs, p);
}

template<typename ntype, int N>
  typename std::conditional<(N<=vecpars::NADD),pvector<ntype,N>,PmOpV<ntype,N, 2, pvector<ntype,N>, ntype>>::type
inline operator*(pvector<ntype,N> const& lhs, const ntype& rhs)
{
  if constexpr (N <= vecpars::NADD)
    {
      return lhs.vecscal(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return PmOpV<ntype, N, 2, pvector<ntype,N>, ntype>(lhs, rhs);
    }
}
// PmOpV divided by scalar
template <typename ntype, int N, int tipo, typename Lhs, typename Rhs>
 inline pvector<ntype, N> operator /(PmOpV<ntype, N, tipo, Lhs, Rhs> const& lhs, const ntype& param) 
    {
      int i;
      pvector<ntype,N> vL;
      for (i=0; i < N; i++)
	  vL[i] = lhs.get_v(i)/param;
      return vL;
    }
#endif

/* some predefined vectors useful for simulations */
typedef pvector<double,2> vecd2;
typedef pvector<double,3> vecd3;
typedef pvector<double,4> vecd4;
#endif
