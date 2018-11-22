#ifndef _PVECTOR_
#define _PVECTOR_
#include <limits>
#include<cstdlib>
#include<iostream>
#include<iomanip>
using namespace std;
#define VEC_MOVE_SEMANTIC
#define VEC_LAZY_EVAL
#ifdef VEC_LAZY_EVAL
// possible types of operation that can be made lazy
struct VOpTypes
{  
  static const int VecPlusVec=0, VecMinusVec=1, VecTimesScal=2, ScalTimesVec=3, VecDivScal=4;
};
//REMARK on Lazyness Logic:
//each class it is an operation of Lhs and Rhs which store the values of the operands
//only when the method get_v is called the operation is performed 
//and this happens only when assign operator of pvector class is called (see below)
// or constructor with AnOpV argument is called
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>
struct AnOpV {
    Lhs const& lhs;
    Rhs const& rhs;
    AnOpV(Lhs const& lhs, Rhs const& rhs):lhs(lhs), rhs(rhs) {
      // empty body
      //cout << "qui\n";
    }
    // lazy operations
    inline ntype get_v(int i) const {
      if constexpr (tipo == VOpTypes::VecPlusVec)// vec + vec
	return lhs.get_v(i) + rhs.get_v(i);
      else if constexpr (tipo == VOpTypes::VecMinusVec) // vec - vec
	return lhs.get_v(i) - rhs.get_v(i);
      else if constexpr (tipo == VOpTypes::VecTimesScal) // vec * scalar
	{
	  return lhs.get_v(i)*rhs;
	}
      else if constexpr (tipo == VOpTypes::ScalTimesVec) // scalar * vec 
	{
	  return rhs.get_v(i)*lhs;
	}
      else if constexpr (tipo == VOpTypes::VecDivScal) // scalar div vec 
	{
	  return lhs.get_v(i)/rhs;
	}
      return 0;
#ifdef _I_DO_NOT_WANT_TO_IMPLEMENT_THIS_
      else if constexpr (tipo == 5) // mat * vec 
	{
	  // get two matrices and the multiply them
	  Lhs m1((double (*)[N])lhs.get_mptr()),
	  Rhs m2((double*)lhs.get_vptr());
	  return matvec(Lhs,Rhs)[i];
	}
      else if constexpr (tipo == 6) // vec * mac 
	return 0;
#endif
    }
    Lhs const& get_lhs() const { return lhs; }
    Rhs const& get_rhs() const { return rhs; }
};
#endif
//NMAX_=45, int MAXSTRA_=40, int NMAXINV_=40, int NMAXMUL_=10, int NSTA_=8, int NLAZY_=5
struct vecpars
{
  // params for my class
  static const int NMAX = 1000; // if N>NMAX use dynamic allocation of matrices
  static const int NSTA = 50;// value below which declares static object in operator funcs (only if m is not dynamically allocated)
  static const int NLAZY = 25;// value below which does not use lazy evaluation for addition
};

// template for a static vector
template <class ntype, int N> class pvecsta {
public:
  ntype v[N];
  constexpr static int dynamic = false;
};
// template for a dynamic vector (i.e. allocated with new)
template <class ntype, int N> class pvecdyn {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  ntype (*v); 
  // copy assignment
  pvecdyn<ntype,N>& operator=(const pvecdyn<ntype,N>& v1)
    {
      int i;
      for (i=0; i < N; i++)
	  v[i] = v1.v[i];
      return (*this);
    }
#ifdef VEC_MOVE_SEMANTIC
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
      v = new ntype[N];
      (*this) = v1;
    }
#ifdef VEC_MOVE_SEMANTIC
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
//choose dynamnic or static base class depending on N
template <class ntype, int N> using pvecbase = typename std::conditional<(N>vecpars::NMAX),pvecdyn <ntype, N>,pvecsta <ntype, N>>::type;
template <class ntype, int N> class pvector : public pvecbase<ntype,N>, vecpars {
  //int vsize;
public:
  using pvecbase<ntype,N>::v;
  static const int maxdigits=numeric_limits<ntype>::digits10;
#ifdef VEC_LAZY_EVAL
  template<typename Lhs, typename Rhs, int tipo>
  // += operator triggers evaluation of lazy expression
  inline pvector<ntype,N> operator += (AnOpV<ntype,N, tipo, Lhs, Rhs> const& op)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] += op.get_v[i];
      return *this;
    }

  // -= operator triggers evaluation of lazy expression
  template<typename Lhs, typename Rhs, int tipo>
  inline pvector<ntype,N> operator -= (AnOpV<ntype,N, tipo, Lhs, Rhs> const& op)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] -= op.get_v[i];
      return *this;
    }
  // assignment operator triggers evaluation of lazy expression
  template<typename Lhs, typename Rhs, int tipo>
    pvector<ntype,N>& operator=(AnOpV<ntype,N, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	v[i] = op.get_v(i);
      return (*this);
    }
  // constructor with AnOpV argument triggers evaluation of lazy expression
  template<typename Lhs, typename Rhs, int tipo>
    pvector<ntype,N>(AnOpV<ntype,N, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	v[i] = op.get_v(i);
    }
  pvector<ntype,N>() = default;
  inline ntype get_v(int i) const {return v[i];}
  inline ntype get_vptr() const {return v;}
#endif

#ifndef VEC_LAZY_EVAL
  // these methods are used if lazyness is not used at all
  inline pvector<ntype,N> operator +(const pvector<ntype,N>& param)
    {
      return addition(param);
    }
  inline pvector<ntype,N> operator -(const pvector<ntype,N>& param)
    {
      return subtraction(param);
    }
#endif
  // methods used if lazyness is disabled (e.g. for vectors with few elemenits)
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

	  if constexpr (N > NMAX)
	    return std::move(v2);
	  else
	    return v2;
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
	  if constexpr (N > NMAX)
	    return std::move(v2);
	  else
	    return v2;
	}
    }

  // += for pvectors 
  inline pvector<ntype,N> operator += (const pvector<ntype,N>& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] += param.v[i];
      return *this;
    }
  // -= for pvectors
  inline pvector<ntype,N> operator -= (const pvector<ntype,N>& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] -= param.v[i];
      return *this;
    }

#ifndef VEC_LAZY_EVAL
  // used if lazyned is completely disabled
  // vec times scalar
  inline pvector<ntype,N> operator *(const double& param) 
    {
      int i;
      pvector<ntype,N> lv;
      for (i=0; i < N; i++)
	lv.v[i] = v[i]*param;

      if constexpr (N > NMAX)
	return std::move(lv);
      else
	return lv;
    }
#endif
  inline pvector<ntype,N> vecscal(const double& param) const
    {
      int i;
      pvector<ntype,N> lv;
      for (i=0; i < N; i++)
	lv.v[i] = v[i]*param;
      if constexpr (N > NMAX)
	return std::move(lv);
      else
	return lv;
    }
#ifndef VEC_LAZY_EVAL
  friend pvector<ntype,N> operator *(const ntype& param, const pvector<ntype,N>& v) 
    {
      int i;
      pvector<ntype,N> lv;
      for (i=0; i < N; i++)
	lv.v[i] = v.v[i]*param;

      if constexpr (N > NMAX)
	return std::move(lv);
      else
	return lv;
    }
#endif
  inline ntype& operator[](const int& i)
    {
      return v[i];
    }
  // scalar product
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
  // norm 
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
#ifdef VEC_LAZY_EVAL
// ADDITION
// AnOpV plus vector 
  template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(AnOpV<ntype, N, tipo, Lhs, Rhs> const& lhs, pvector<ntype,N> const& p)
{
  return AnOpV<ntype, N, VOpTypes::VecPlusVec, AnOpV<ntype, N, tipo, Lhs, Rhs>, pvector<ntype,N>>(lhs, p);
} 
// add expression template with point at the left
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(pvector<ntype,N> const& p, AnOpV<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpV< ntype, N, VOpTypes::VecPlusVec, pvector<ntype,N>, AnOpV<ntype, N, tipo, Lhs, Rhs> >(p, rhs);
}

// vector plus vector 
// se N < NSTA restituisce pvector e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpV e usa la lazy evaluation
template <typename ntype, int N> 
// conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  typename std::conditional<(N<=vecpars::NLAZY),pvector<ntype,N>,AnOpV<ntype, N, VOpTypes::VecPlusVec, pvector<ntype,N>, pvector<ntype,N>>>::type
inline operator+(pvector<ntype,N> const& lhs, pvector<ntype,N> const& rhs) 
{
  if constexpr (N <= vecpars::NLAZY)
    {
      return lhs.addition(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, N, VOpTypes::VecPlusVec, pvector<ntype,N>, pvector<ntype,N>>(lhs, rhs); // lazy evaluation
    }
}
//AnOpV plus AnOpV
template <typename ntype, int N, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator+(const AnOpV<ntype, N, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpV<ntype, N, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpV<ntype, N, VOpTypes::VecPlusVec, AnOpV<ntype, N, tipoL, LLhs, LRhs>, AnOpV<ntype, N, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
// SUBTRACTION
// AnOpV minus vector
  template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator-(AnOpV<ntype,N, tipo, Lhs, Rhs> const& lhs, pvector<ntype,N> const& p)
{
  return AnOpV<ntype, N, VOpTypes::VecMinusVec, AnOpV<ntype, N, tipo, Lhs, Rhs>, pvector<ntype,N>>(lhs, p);
} 
// vector minus AnOpV 
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>
inline auto operator-(pvector<ntype,N> const& p, AnOpV<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpV< ntype, N, VOpTypes::VecMinusVec, pvector<ntype,N>, AnOpV<ntype, N, tipo, Lhs, Rhs> >(p, rhs);
}

// vector minus vector
// se N < NSTA restituisce pvector e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpV e usa la lazy evaluation
template <typename ntype, int N> 
// conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  typename std::conditional<(N<=vecpars::NLAZY),pvector<ntype,N>,AnOpV<ntype, N, VOpTypes::VecMinusVec, pvector<ntype,N>, pvector<ntype,N>>>::type
inline operator-(pvector<ntype,N> const& lhs, pvector<ntype,N> const& rhs) 
{
  if constexpr (N <= vecpars::NLAZY)
    {
      return lhs.subtraction(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << NSTA << "\n";
      return AnOpV<ntype, N, VOpTypes::VecMinusVec, pvector<ntype,N>, pvector<ntype,N>>(lhs, rhs); // lazy evaluation
    }
}
// AnOpV time AnOpV
  template <typename ntype, int N, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator-(const AnOpV<ntype, N, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpV<ntype, N, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpV<ntype, N, VOpTypes::VecMinusVec, AnOpV<ntype, N, tipoL, LLhs, LRhs>, AnOpV<ntype, N, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
   
///
// Scalar product of two AnOpV
//
  template <typename ntype, int N, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline ntype operator*(const AnOpV<ntype, N, tipoL, LLhs, LRhs> & leftOperand, const AnOpV<ntype, N, tipoR, RLhs, RRhs> & rightOperand)
{
  pvector<ntype, N> vR, vL;
  for (int i=0; i < N; i++)
    {
      vL[i] = leftOperand.get_v(i);
      vR[i] = rightOperand.get_v(i);
    }
  return vL*vR;
}
// AnOpV times pvector
//
   template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> ntype 
inline operator*(AnOpV<ntype,N, tipo, Lhs, Rhs> const& lhs, pvector<ntype,N> const& vR)
{
  pvector<ntype, N> vL;
  for (int i=0; i < N; i++)
    vL[i] = lhs.get_v(i);
  return vL*vR;
} 
// vector times right AnOpV
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs>  
inline ntype operator*(pvector<ntype,N> const& vL, AnOpV<ntype, N, tipo, Lhs, Rhs> const& rhs) 
{
  pvector<ntype, N> vR;
  for (int i=0; i < N; i++)
    {
      vR[i] = rhs.get_v(i);
    }
  return vL*vR;
}
// scalar times AnOvV 
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(const ntype& lhs, AnOpV<ntype, N, tipo, Lhs, Rhs> const& rhs)
{
  return AnOpV<ntype, N, VOpTypes::ScalTimesVec, ntype, AnOpV<ntype, N, tipo, Lhs, Rhs>>(lhs, rhs);
}
// scalar times vector
template<typename ntype, int N>
  typename std::conditional<(N<=vecpars::NLAZY),pvector<ntype,N>,AnOpV<ntype,N, VOpTypes::ScalTimesVec, ntype, pvector<ntype,N>>>::type
inline operator*(const ntype& lhs, pvector<ntype,N> const& rhs)
{
  if constexpr (N <= vecpars::NLAZY)
    {
      return rhs.vecscal(lhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, N, VOpTypes::ScalTimesVec, ntype, pvector<ntype,N>>(lhs, rhs);
    }
}

// AnOpV times scalar
template<typename ntype, int N, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(AnOpV<ntype, N, tipo, Lhs, Rhs> const& lhs, const ntype& p)
{
  return AnOpV<ntype, N, VOpTypes::VecTimesScal, AnOpV<ntype, N, tipo, Lhs, Rhs>, ntype>(lhs, p);
}
// PmOvV times scalar
template<typename ntype, int N>
  typename std::conditional<(N<=vecpars::NLAZY),pvector<ntype,N>,AnOpV<ntype,N, VOpTypes::VecTimesScal, pvector<ntype,N>, ntype>>::type
inline operator*(pvector<ntype,N> const& lhs, const ntype& rhs)
{
  if constexpr (N <= vecpars::NLAZY)
    {
      return lhs.vecscal(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, N, VOpTypes::VecTimesScal, pvector<ntype,N>, ntype>(lhs, rhs);
    }
}
// vector divided by scalar 
  template <typename ntype, int N >
typename std::conditional<(N<=vecpars::NLAZY),pvector<ntype,N>,AnOpV<ntype,N, VOpTypes::VecDivScal, pvector<ntype,N>, ntype>>::type
inline operator /(pvector<ntype,N> const& lhs, const ntype& rhs) 
{
  if constexpr (N <= vecpars::NLAZY)
    {
      return lhs.vecscal(1.0/rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, N, VOpTypes::VecDivScal, pvector<ntype,N>, ntype>(lhs, rhs);
    }
}

// AnOpV divided by scalar
  template <typename ntype, int N, int tipo, typename Lhs, typename Rhs>
inline auto operator /(AnOpV<ntype, N, tipo, Lhs, Rhs> const& lhs, const ntype& rhs) 
{
  return AnOpV<ntype, N, VOpTypes::VecDivScal, AnOpV<ntype, N, tipo, Lhs, Rhs>, ntype>(lhs, rhs);
}
#endif

/* some predefined vectors useful for simulations */
typedef pvector<double,2> vecd2;
typedef pvector<double,3> vecd3;
typedef pvector<double,4> vecd4;
#endif
