#include "./pvector.H"
#include "pmatrix.H"
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<cmath>
//#include "pvector.H"
#include <algorithm> 
#include <limits>
#include <cstdlib>
#include <vector>
#include <array>
//#include <boost/multiprecision/mpfr.hpp>

//#include <complex>
//fix to avoid crash of clang
#define Sqr(x) ((x)*(x))
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
//#define SIGN(a,b) ((b) >= 0.0 ? abs(a) : -abs(a))
template <class ntype>
inline ntype SIGN(ntype a, ntype b)
{
  if ((b) > 0)
    return abs(a);
  else
    return -abs(a);
}
//#define USE_DEFL_NR // metterlo o no non cambia nulla
//#define USE_FULL_NR // peggiora l'accuratezza attivandolo
//#define LAGUERRE_POLISH // migliora drasticamente con maxiter=8 da tenere!
//#define ENABLE_FORW_DEFL // non usando la forward deflation è più accurato!
#define USE_BRENT
#define USE_LAGUERRE_ONLY 
#define POLISH_NR_REAL
#define BINI_CONV_CRIT
//#define FAST_MATH
using namespace std;
template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

template <class ntype, int N> 
class rpoly_base_static 
{
public:
  int n;
  constexpr static int dynamic = false;
  pvector<ntype, N+1> coeff;
  pvector<ntype, N+1> cmon;
 // pvector<ntype, N> cmonder;
 // pvector<ntype, N-1> cmondder;
  pvector<ntype, N-1> bd;
  pvector<ntype, N-1> bd0;
  pvector<ntype, N-1> bd1;
  pvector<ntype, N> deflcoeff0;
  pvector<ntype, N> deflcoeff1;
  pvector<ntype, N> deflcoeff;
  void set_coeff(pvector<ntype,N+1> v)
    {
      coeff = v;
      //coeff.show("coefficienti=");
    }
  rpoly_base_static()
    {
      n=N;
    }
  rpoly_base_static(int nc): rpoly_base_static()
  {
    n=nc;
    // empty body
  }
};
template <class ntype, int N> 
class rpoly_base_dynamic 
{
public:
  int n;
  constexpr static int dynamic = true;
  pvector<ntype> coeff;
  pvector<ntype> cmon;
  //pvector<ntype> cmonder;
  //pvector<ntype> cmondder;
  pvector<ntype> bd0;
  pvector<ntype> bd1;
  pvector<ntype> bd;
  pvector<ntype> deflcoeff0;
  pvector<ntype> deflcoeff1;
  pvector<ntype> deflcoeff;
  rpoly_base_dynamic() = default;
  void set_coeff(pvector<ntype,-1> v)
    {
      coeff = v;
#if 0
      coeff.show();
      cout << "n=" << n << "\n";
      cout << "v=" << v[n]<< "\n";
        v.show();
#endif
    }

  void use_vec(int nc, ntype* coeffv, ntype* cmonv,
               ntype* bd0v, ntype *bd1v, ntype* bdv, ntype *deflcoeff0v, ntype* deflcoeff1v, ntype* deflcoeffv)
    {
      coeff.use_vec(nc+1,coeffv);
      cmon.use_vec(nc+1,cmonv);
      //cmonder.use_vec(nc,cmonderv);
      //cmondder.use_vec(nc-1,cmondderv);
      bd0.use_vec(nc-1,bd0v);
      bd1.use_vec(nc-1,bd1v);
      bd.use_vec(nc-1,bdv);
      deflcoeff0.use_vec(nc,deflcoeff0v);
      deflcoeff1.use_vec(nc,deflcoeff1v);
      deflcoeff.use_vec(nc,deflcoeffv);
      n=nc;
    }
 
  rpoly_base_dynamic(int nc): coeff(nc+1), cmon(nc+1), bd0(nc-1), bd1(nc-1), bd(nc-1),
  deflcoeff0(nc), deflcoeff1(nc), deflcoeff(nc)
    {
      n=nc;
    }
  ~rpoly_base_dynamic() = default;
  void deallocate(void)
    {
      coeff.deallocate();
      cmon.deallocate();
      //cmonder.deallocate();
      //cmondder.deallocate();
      bd0.deallocate();
      bd1.deallocate();
      bd.deallocate();
      deflcoeff0.deallocate();
      deflcoeff1.deallocate();
      deflcoeff.deallocate();
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      cmon.allocate(n+1);
      //cmonder.allocate(n);
      //cmondder.allocate(n-1);
      bd0.allocate(n-1);
      bd1.allocate(n-1);
      bd.allocate(n-1);
      deflcoeff0.allocate(n);
      deflcoeff1.allocate(n);
      deflcoeff.allocate(n);
    }
};

template <class ntype, int N> using rpolybase = 
typename std::conditional<(N>0), rpoly_base_static <ntype, N>,
	 rpoly_base_dynamic <ntype, N>>::type;
 
template <class ntype, int N=-1, bool force_hqr=false, class cmplx=complex<ntype>> 
class rpoly: public numeric_limits<ntype>, public rpolybase<ntype,N> {
  using rpolybase<ntype,N>::n;
  using rpolybase<ntype,N>::coeff;
  using rpolybase<ntype,N>::cmon;
  //using rpolybase<ntype,N>::cmonder;
  //using rpolybase<ntype,N>::cmondder;
  using rpolybase<ntype,N>::bd0;
  using rpolybase<ntype,N>::bd1;
  using rpolybase<ntype,N>::bd;
  using rpolybase<ntype,N>::deflcoeff0;
  using rpolybase<ntype,N>::deflcoeff1;
  using rpolybase<ntype,N>::deflcoeff;
#if 0
#ifdef USE_MULTIPREC
  static const int numlin=std::numeric_limits<ntype>::digits10;
#ifdef MP1
  using cmplx=typename std::conditional<boost::is_same<ntype,number<mpfr_float_backend<numlin>>>::value,
    number<mpc_complex_backend<std::numeric_limits<ntype>::digits10>>,complex<ntype>>::type;
#else
  using cmplx=typename std::conditional<boost::is_same<ntype,number<cpp_bin_float<numlin>>>::value,
    cpp_complex<std::numeric_limits<ntype>::digits10>,complex<ntype>>::type;
#endif
#endif
#endif
  //using cmplx=cmplx;
  template <class vtype>
  using pvecNm1 = typename std::conditional<(N>0),  pvector<vtype, N-1>,
	 pvector<vtype, -1>>::type;
  template <class vtype>
  using pvecNp1 = typename std::conditional<(N>0),  pvector<vtype, N+1>,
	 pvector<vtype, -1>>::type;
  template <class vtype>
  using pvecNm2 = typename std::conditional<(N>0),  pvector<vtype, N-2>,
	 pvector<vtype, -1>>::type;
  using rpolyNm1 = typename std::conditional<(N>0),  rpoly<ntype, N-1,force_hqr,cmplx>,
	 rpoly<ntype, -1,force_hqr,cmplx>>::type;
  using rpolyNm2 = typename std::conditional<(N>0),  rpoly<ntype, N-2,force_hqr,cmplx>,
	 rpoly<ntype, -1,force_hqr,cmplx>>::type;
  template <class vtype>
  using pmatNm1 = typename std::conditional<(N>0),  pmatrixq<vtype, N-1>,
	 pmatrixq<vtype, -1>>::type;

  template <class vtype>
  using pmatNm2 = typename std::conditional<(N>0),  pmatrixq<vtype, N-2>,
	 pmatrixq<vtype, -1>>::type;
 template <class vtype, int NT>
  using stlarr = typename std::conditional<(NT>0), std::array<vtype,NT-1>, std::vector<vtype>>::type;

  ntype px, dpx; 

  const int maxiter_polish=8;
  //force_hqr;
  int imaxarg1,imaxarg2;
  ntype eps05, meps, maxf, maxf2, maxf3, scalfact, cubic_rescal_fact;
  int maxdigits;
  ntype goaleps;
  bool deflated;
#ifdef BINI_INIGUESS
  cmplx *iniguess;
#endif
  ntype rtsafe(const ntype x1, const ntype x2, const ntype xacc)
   {
      const int MAXIT=1000000;
      int j;
      ntype df,dx,dxold,f,fh,fl,temp,xh,xl,rts;

      fl=evalpoly(x1);
      fh=evalpoly(x2);
#if 0
      this->show("p[x_]:=");
      cout << "p[ "<< x2 << "]=" << evalpoly(x2) << "\n";
      cout << "dp[" << x2 << "]=" << evaldpoly(x2) << "\n";
      cmonder.show();
      exit(-1);
#endif
      if ((fl > 0.0  && fh > 0.0) || (fl < 0.0 && fh < 0.0))
        {
          cout << "Root must be bracketed in rtsafe\n";
          exit(-1);
        }
      if (fl == 0.0) 
        return x1;
      if (fh == 0.0) 
        return x2;
      if (fl < 0.0) {
        xl=x1;
        xh=x2;
      } else {
        xh=x1;
        xl=x2;
      }
      rts=0.5*(x1+x2);
      dxold=fabs(x2-x1);
      dx=dxold;
      //funcd(rts,f,df);
      f = evalpoly(rts);
      df = evaldpoly(rts);
      for (j=0;j<MAXIT;j++) {
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
            || (fabs(2.0*f) > fabs(dxold*df))) {
          dxold=dx;
          dx=0.5*(xh-xl);
          rts=xl+dx;
          if (xl == rts) return rts;
        } else {
          dxold=dx;
          dx=f/df;
          temp=rts;
          rts -= dx;
          if (temp == rts) return rts;
        }
        if (fabs(dx) < xacc) return rts;
        //funcd(rts,f,df);
        f = evalpoly(rts);
        df = evaldpoly(rts);

        if (f < 0.0)
          xl=rts;
        else
          xh=rts;
      }
      cout << "Maximum number of iterations exceeded in rtsafe\n";
      return 0.0;
    }
  //ntype (rpoly<ntype,N>::*func) (ntype x);
  ntype zbrent(const ntype x1, const ntype x2, const ntype tol)
  {
    const int ITMAX=1000000;
    const ntype EPS=numeric_limits<ntype>::epsilon();
    int iter;
    ntype a=x1,b=x2,c=x2,d,e,min1,min2;
    ntype fa=evalpoly(a),fb=evalpoly(b),fc,p,q,r,s,tol1,xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
      {	
        cout << "Root must be bracketed in zbrent";
        cout << "fa=" << fa << " fb=" << fb << "\n";
        exit(-1);
      }
    fc=fb;
    for (iter=0;iter<ITMAX;iter++)
      {
       if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
          c=a;
          fc=fa;
          e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) 
          {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
          }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
          s=fb/fa;
          if (a == c) 
            {
              p=2.0*xm*s;
              q=1.0-s;
            } 
          else 
            {
              q=fa/fc;
              r=fb/fc;
              p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
              q=(q-1.0)*(r-1.0)*(s-1.0);
            }
          if (p > 0.0) q = -q;
          p=fabs(p);
          min1=3.0*xm*q-fabs(tol1*q);
          min2=fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2)) 
            {
              e=d;
              d=p/q;
            } 
          else 
            {
              d=xm;
              e=d;
            }
        } 
        else
          {
            d=xm;
            e=d;
          }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
          b += d;
        else
          b += SIGN<ntype>(tol1,xm);
        fb=evalpoly(b);
      }
    cout << "Maximum number of iterations exceeded in zbrent";
    //exit(-1);
    return 0.0;
  }
  ntype oqs_max2(ntype a, ntype b)
    {
      if (a >= b)
	return a;
      else
	return b;
    }
  ntype oqs_max3(ntype a, ntype b, ntype c)
    {
      ntype t;
      t = oqs_max2(a,b);
      return oqs_max2(t,c);
    }
  //void oqs_quartic_solver(pvector<cmplx,N>& roots);
  inline void solve_quadratic(pvector<cmplx,N>& sol);
  inline void solve_cubic_analytic(pvector<cmplx,N>& sol);
  inline void oqs_quartic_solver(pvector<cmplx,N>& roots);
  inline void oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol);
  inline void oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol);
  inline void oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled);
  inline ntype oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3);
  inline ntype oqs_calc_err_abcd_cmplx(ntype a, ntype b, ntype c, ntype d, 
	    			       cmplx aq, cmplx bq, cmplx cq, cmplx dq);
  inline ntype oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq);
  inline ntype oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq); 
  inline void oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ);
  inline void oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2]);
  inline void balance(pmatrixq<ntype,N>& a);
  inline void hqr(pmatrixq<ntype,N>& a, pvector<cmplx,N>& wri, int& ok);
  inline void QRfactorization(pmatrixq<ntype,N>& hess, pvector<cmplx,N>& sol, int& ok);
  inline void solve_numrec(pvector<cmplx,N>& csol, int& ok);
public:
  void show(void)
    {
      show(NULL);
    }
  void show(const char* str)
    {
      int i;
      if (str!=NULL)
	cout <<  str;
      for (i=n; i >= 0; i--)
	{
	  if (coeff[i] > 0)
    	    {
	      if (i < n)
		cout << "+";
	    }
	  else
	    { 
	      cout << "-";
	    }
	  if (i==0)
	    cout << setprecision(maxdigits) << abs(coeff[i]);
	  else if (i > 0 && abs(coeff[i]) != 1.0)
	    cout << setprecision(maxdigits) << abs(coeff[i])<< "*";
	 
	  if ( i > 1)
	    {
	      cout << "x^" << i;
	    }
	  else if (i==1)
	    cout << "x";
	}
      cout << "\n";
    }
   cmplx evalpoly(cmplx x)
    {
      // evaluate polynomail via Horner's formula 
      cmplx bn=cmplx(0.0);
      for (int i=n; i >= 0; i--)
        {
          bn = cmon[i] + bn*x;
        }
      return bn;
    }
   cmplx evaldpoly(cmplx x)
    {
      // evaluate polynomail via Horner's formula 
      cmplx bn=0.0;
      for (int i=n-1; i >= 0; i--)
        {
          bn = (i+1)*cmon[i+1] + bn*x;
        }
      return bn;
    }


   ntype evalpoly(ntype x)
    {
      // evaluate polynomail via Horner's formula 
      ntype bn=0.0;
      for (int i=n; i >= 0; i--)
        {
          bn = cmon[i] + bn*x;
        }
      return bn;
    }

  ntype evaldpoly(ntype x)
    {
      // evaluate first derivative of polynomail via Horner's formula 
      ntype bn=0.0;

      for (int i=n-1; i >= 0; i--)
        {
          bn = (i+1)*cmon[i+1] + bn*x;
        }
      return bn;
    }
  ntype evalddpoly(ntype x)
    {
      // evaluate second derivative of polynomail via Horner's formula 
      ntype bn=0.0;
      if (n == 1)
        return 0;
      for (int i=n-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmon[i+2] + bn*x;
        }
      return bn;
    }
  ntype calc_err_q(pvector<ntype,N> c, ntype r0)
    {
      int i;
      ntype sum=0;
      for (i=0; i < n; i++)
        {
          if (i == 0)
            {
              if (cmon[i]==0)
                {
                  sum+=abs(r0*c[i]+cmon[i]);
                } 
              else 
                {
                  sum+=abs((r0*c[i]+cmon[i])/cmon[i]);
                }
            }
          else 
            {
              if (cmon[i]==0)
                sum +=abs(-c[i-1] + r0*c[i] + cmon[i]);
              else
                sum +=abs((-c[i-1] + r0*c[i] + cmon[i])/cmon[i]);
            }
        }
      //sum += (cmon[0]==0.0)?abs(r0*c[0]+cmon[0]):abs((r0*c[0]+cmon[0])/cmon[0]);
      return sum;
    }
#ifdef USE_FULL_NR
  void use_NR_defl_lin(pvector<ntype,N>& defc, cmplx& r0)
    {
      ntype errf, errfold, det;
      pvector<ntype,N> c, cold;
      pvector<ntype,N> fvec;
      pmatrixq<ntype,N> Jinv, J;
      int i, j, iter;

      if constexpr (N < 0)
        {
          c.allocate(n);
          cold.allocate(n);
          fvec.allocate(n);
          Jinv.allocate(n);
          J.allocate(n);
        }
      for (i=0; i < n-1; i++)
        c[i]=defc[i];
      c[N-1] = r0.real();
      for (i=0; i < n; i++)
        {
          fvec[i] = ((i>=1)?(c[i-1]):0.0)-r0.real()*((i<n-1)?c[i]:1.0)-cmon[i];
        } 
      //cmon.show("monic coeff=");
      //defc.show("deflcoeff");
      //fvec.show("errori=");
      errf = 0;
      for (i=0; i < n; i++)
        {
          errf += abs((cmon[i]==0)?fvec[i]:(fvec[i]/cmon[i]));
        }
      //cout << "errf inizio=" << errf << "\n";
     for (iter = 0; iter < 8; iter++)
        {
          cold = c;
          for (i=0; i < n; i++)
            for (j=0; j < n; j++)
              {
                if (i==j)
                  {
                    J[i][j] = (i<n-1)?(-c[n-1]):((ntype)-1.0);// c[N-1]=r0
                  }                
                else if (j==i+1)
                  J[i][j] = 1.0;
                else 
                  J[i][j] = 0;
                if (j==n-1&& i!=j)
                  J[i][j] = -((i<n-1)?c[i]:((ntype)1.0));
              }
          //J.show();
          if (J.det()==0)
            break;
          Jinv=J.inv();

          c = c - Jinv*fvec;
          for (i=0; i < n; i++)
            {
              fvec[i] = ((i>=1)?(c[i-1]):0.0)-c[n-1]*((i<n-1)?c[i]:1.0)-cmon[i];
            } 
          errfold = errf;
          errf = 0;
          for (i=0; i < n; i++)
            {
              errf += abs((cmon[i]==0)?fvec[i]:(fvec[i]/cmon[i]));
            }
          //cout << "errf iter=" <<iter << " :"<< errf << " errfold=" << errfold << "\n";
          if (errf==0)
            break;
          if (errf > errfold)
            {
              //cout << "iter=" << iter << "errf=" << errf << " errfold=" << errfold << "\n";
              c = cold;
              break;
            }
          
        } 
      for (i=0; i < n-1; i++)
        defc[i] = c[i];
      r0=cmplx(c[n-1],r0.imag());
    }

#else
void use_NR_defl_lin(pvector<ntype,N>& defc, cmplx r0)
  {
    ntype errf, errfold, det;
    pvecNm1<ntype> c, cold;
    pvecNm1<ntype> fvec;
    pmatNm1<ntype> Jinv;
    int i, j, iter;
    if constexpr (N < 0)
      {
        c.allocate(n-1);
        cold.allocate(n-1);
        fvec.allocate(n-1);
        Jinv.allocate(n-1);
      }
    for (i=n-2; i >= 0; i--)
      {
        c[i]=defc[i];
        fvec[i] = defc[i]-r0.real()*((i==n-2)?1:defc[i+1])-cmon[i+1];
      } 
    //cmon.show("monic coeff=");
    //defc.show("deflcoeff");
    //fvec.show("errori=");
    errf = 0;
    for (i=0; i < n-1; i++)
      {
        errf += abs((cmon[i+1]==0)?fvec[i]:(fvec[i]/cmon[i+1]));
      }
    //cout << "errf inizio=" << errf << "\n";
    for (i=0; i < n-1; i++)
      for (j=0; j < n-1; j++)
        Jinv[i][j] = (i==j)?1.0:0.0;
    for (j=0; j < n-2; j++)
      for (i=j+1; i < n-1; i++)
        {
          Jinv[i][j]=r0.real()*Jinv[i-1][j];
        }
    // Jinv.show();
    for (iter = 0; iter < 8; iter++)
      {
        cold = c;
        c = c - Jinv*fvec;
        for (i=n-2; i >= 0; i--)
          {
            fvec[i] = c[i]-r0.real()*((i==n-2)?1:c[i+1])-cmon[i+1];
          } 
        errfold = errf;
        errf = 0;
        for (i=0; i < n-1; i++)
          {
            errf += abs((cmon[i+1]==0)?fvec[i]:(fvec[i]/cmon[i+1]));
          }
        //cout << "errf iter=" <<iter << " :"<< errf << " errfold=" << errfold << "\n";
        if (errf==0)
          break;
        if (errf > errfold)
          {
            //cout << "iter=" << iter << "errf=" << errf << " errfold=" << errfold << "\n";
            c = cold;
            break;
          }

      } 
    for (i=0; i < n-1; i++)
      defc[i] = c[i];
  }

#endif
#ifdef USE_FULL_NR
void use_NR_defl_quad(pvecNm1<ntype>& bd, pvector<ntype,3>& ad)
  {
    // bd[N-2] = 1 since all polynomials are in their minic form
    // hence we do not need to refine this coefficient here
    ntype term, errf, errfold, det;
    pvector<ntype,N> cold, c;
    pvector<ntype,N> fvec;
    pmatrixq<ntype,N> Jinv, J;
    int i, j, iter;
    ntype sum=0;
    if constexpr (N < 0)
      {
        cold.allocate(n);
        c.allocate(n);
        fvec.allocate(n);
        Jinv.allocate(n);
        J.allocate(n);
      }
    for (i=0; i < n-2; i++)
      c[i]=bd[i];
    c[n-2] = ad[0];
    c[n-1] = ad[1];
    for (i=0; i < n; i++)
      {
        fvec[i] = -cmon[i];
        if (i <= n-2)
          {
            fvec[i] += ((i==n-2)?1.0:c[i])*ad[0];
          }
        if (i >= 1) 
          fvec[i] += ((i==n-1)?1.0:c[i-1])*ad[1];
        if (i >=2)
          fvec[i] += c[i-2];
      }
    //cmon.show("monic coeff=");
    //defc.show("deflcoeff");
    //fvec.show("errori=");
    errf = 0;
    for (i=0; i < n; i++)
      {
        errf += abs((cmon[i]==0)?fvec[i]:(fvec[i]/cmon[i]));
      }
    //cout << "errf inizio=" << errf << "\n";
   #if 0
    pvector<ntype, N-2> b;
    for (i=0; i < N-2; i++)
      b[i] = -cmon[i+2];
    c = Jinv*b;
    for (i=0; i < N-2; i++)
      bd[i]=c[i];
    return;
#endif
    // Jinv.show();
    for (iter = 0; iter < 8; iter++)
      {
        cold = c;
        for (i=0; i < n; i++)
          for (j=0; j < n; j++)
            {
              if (i==j)
                { 
                  if (i < n-2)
                    J[i][j] = c[n-2]; // c[N-2]->ad[0]
                  else 
                    J[i][j] = 1.0;
                }
              else if (j==i+1)    // c[N-1]->ad[1]
                {
                  J[i][j] = c[n-1];
                }
              else if (j==i+2)
                J[i][j] = 1.0;
              else
                J[i][j] = 0.0;
              if (j==n-2 && j!=i)
                J[i][j] = (i < n-2)?c[i]:1.0;
              if (j==n-1 && j!=i && i >= 1)
                J[i][j] = (i < n-1)?c[i-1]:1.0; 
            }
        if (J.det()==0)
          break;
        Jinv = J.inv();
       
        c = c - Jinv*fvec;
        for (i=0; i < n; i++)
          {
            fvec[i] = -cmon[i];
            if (i <= n-2)
              fvec[i] += ((i==n-2)?1.0:c[i])*c[n-2];
            if (i >= 1) 
              fvec[i] += ((i==n-1)?1.0:c[i-1])*c[n-1];
            if (i >=2)
              fvec[i] += c[i-2];
          }
        errfold = errf;
        for (i=0; i < n; i++)
          {
            errf += abs((cmon[i]==0)?fvec[i]:(fvec[i]/cmon[i]));
          }
        //cout << "errf iter=" <<iter << " :"<< errf << " errfold=" << errfold << "\n";
        if (errf==0)
          break;
        if (errf > errfold)
          {
            //cout << "iter=" << iter << "errf=" << errf << " errfold=" << errfold << "\n";
            c = cold;
            break;
          }
      } 
    for (i=0; i < n-2; i++)
      bd[i]=c[i];
    ad[0]=c[n-2];
    ad[1]=c[n-1];
  }


#else
  void use_NR_defl_quad(pvecNm1<ntype>& bd, pvector<ntype,3> ad)
    {
      // bd[N-2] = 1 since all polynomials are in their minic form
      // hence we do not need to refine this coefficient here
      ntype term, errf, errfold, det;
      pvecNm2<ntype> cold, c;
      pvecNm2<ntype> fvec;
      pmatNm2<ntype> Jinv, J;
      int i, j, iter;
      ntype sum=0;
      if constexpr (N < 0)
        {
          cold.allocate(n-2);
          c.allocate(n-2);
          fvec.allocate(n-2);
          Jinv.allocate(n-2);
          J.allocate(n-2);
        }
      for (i=0; i < n-2; i++)
        c[i]=bd[i];
      for (i=2; i < n; i++)
        {
          fvec[i-2] = cmon[i];
          if (i < n-1)
            fvec[i-2] -= ((i==n-2)?1.0:c[i])*ad[0];
          fvec[i-2] -= ((i==n-1)?1.0:c[i-1])*ad[1];
          fvec[i-2] -= c[i-2];
        }
      //cmon.show("monic coeff=");
      //defc.show("deflcoeff");
      //fvec.show("errori=");
      errf = 0;
      for (i=2; i <n; i++)
        {
          errf += abs((cmon[i]==0)?fvec[i-2]:(fvec[i-2]/cmon[i]));
        }
      //cout << "errf inizio=" << errf << "\n";
      for (i=0; i < n-2; i++)
        for (j=0; j < n-2; j++)
          {
            if (i==j)
              J[i][j] = -1.0;
            else if (i==j+1)
              J[i][j] = -ad[1];
            else if (i==j+2)
              J[i][j] = -ad[0];
            else 
              J[i][j] = 0.0;
          }
      Jinv = J.inv();
#if 0
      pvector<ntype, N-2> b;
      for (i=0; i < N-2; i++)
        b[i] = -cmon[i+2];
      c = Jinv*b;
      for (i=0; i < N-2; i++)
        bd[i]=c[i];
      return;
#endif
      // Jinv.show();
      for (iter = 0; iter < 8; iter++)
        {
          cold = c;
          c = c - Jinv*fvec;
          for (i=2; i < n; i++)
            {
              fvec[i-2] = cmon[i];
              if (i < n-1)
                fvec[i-2] -= ((i==n-2)?1.0:c[i])*ad[0];
              fvec[i-2] -= ((i==n-1)?1.0:c[i-1])*ad[1];
              fvec[i-2] -= c[i-2];
            }

          errfold = errf;
          for (i=2; i < n; i++)
            {
              errf += abs((cmon[i]==0)?fvec[i-2]:(fvec[i-2]/cmon[i]));
            }
          //cout << "errf iter=" <<iter << " :"<< errf << " errfold=" << errfold << "\n";
          if (errf==0)
            break;
          if (errf > errfold)
            {
              //cout << "iter=" << iter << "errf=" << errf << " errfold=" << errfold << "\n";
              c = cold;
              break;
            }
        } 
      for (i=0; i < n-2; i++)
        bd[i]=c[i];
    }
#endif
void cmplx_div(ntype c[2], ntype a[2], ntype b[2])
  {
    cmplx ac,bc, acbc;
    ac=cmplx(a[0],a[1]);
    bc=cmplx(b[0],b[1]);
    acbc=ac/bc;
    c[0]=(acbc).real();
    c[1]=(acbc).imag();
  } 
void cmplx_mul(ntype c[2], ntype a[2], ntype b[2])
  {
    c[0] = a[0]*b[0]-a[1]*b[1];
    c[1] = a[0]*b[1]+a[1]*b[0];
  } 
#if 0
bool laguer_real(cmplx &xc, int& flag)
  {
    const int MR=8,MT=10,MAXIT=MT*MR;
    const ntype EPS=goaleps;
    static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    int iter,j, numimp,jrev;
    ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
    ntype abx,abp,abm, phi, rho, s, absb, absbold, xbak[2];
    cmplx rhoc;
    bool rev;
#ifndef BINI_CONV_CRIT
    ntype err;
#endif
    //cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    ntype dx[2],x1[2],b[2],d[2],f[2],g[2],h[2],sq[2],gp[2],gm[2],g2[2], fb[2], pol[2];
    ntype baux[2], x[2], f0, d0, b0, invn, g0, invgp[2], sarg[2], xg[2], xsq[2], xsqh[2];
    int m=n;
    //cout << "r0 ini = "<<xc << "\n"; 
    x[0] = xc.real();
    x[1] = xc.imag();
    flag=0;
    numimp=0;
    for (iter=1;iter<=MAXIT;iter++) 
      {
        //its=iter;
        abx=sqrt(x[0]*x[0]+x[1]*x[1]);
        xbak[0]=x[0];
        xbak[1]=x[1];
        if (abx > 1)
          {
            rev=true;
            rhoc=1.0/cmplx(x[0],x[1]);
            x[0]=real(rhoc);
            x[1]=imag(rhoc);
          }
        else
          {
            rev=false;
          }
        if (rev)
          b[0]=cmon[0];  
        else
          b[0]=cmon[m];
        b[1]=0.0;
#ifndef BINI_CONV_CRIT
        err=sqrt(b[0]*b[0]+b[1]*b[1]);
#endif
        d[0]=f[0]=0.0;
        d[1]=f[1]=0.0;
        for (j=m-1;j>=0;j--) {
          f0 = f[0];
          f[0] = x[0]*f[0] - x[1]*f[1] + d[0];
          f[1] = x[0]*f[1] + x[1]*f0 + d[1];
          d0 = d[0];
          d[0] = x[0]*d[0] - x[1]*d[1] + b[0];
          d[1] = x[0]*d[1] + x[1]*d0 + b[1];
          b0 = b[0];
          if (rev)
            jrev=m-j;
          else
            jrev=j;
          b[0] = x[0]*b[0] - x[1]*b[1] + cmon[jrev];
          b[1] = x[0]*b[1] + x[1]*b0;
          // note that f is P''/2
#ifndef BINI_CONV_CRIT
          err=sqrt(b[0]*b[0]+b[1]*b[1])+abx*err;
#endif
        }
#ifndef BINI_CONV_CRIT
        err *= EPS;
#endif
       if (!rev && (isinf(b[0]) || isnan(b[0]) || 
            isinf(b[1]) || isnan(b[1])))
          {
            //cout << "f= " << f[0] << " " << f[1]  << "d= " << d[0] << " " << d[1] << "b=" << b[0] << " " << b[1] << "\n";
            //cout << "x=" << x[0] << " " << x[1] << "\n";
            //cout << "p(x)=" << evalpoly(x) << "\n";
            flag=2;
            return false;
          }
#ifdef BINI_CONV_CRIT
        //if (iter > 30)
          //cout << "#" << iter << " p(x)=" << abs(cmplx(b[0],b[1])) << "\n";
        s=coeff[n];
        for (j=n-1; j >=0; j--) 
          {
            s=abx*s+coeff[j];
          }
        if (iter > 1)
          absbold=absb;
        if (rev)
          {
            baux[0]=cmon[m];
            baux[1]=0.0;
            for (j=m-1;j>=0;j--)
              {
                b0 = baux[0];
                baux[0] = xbak[0]*baux[0] - xbak[1]*baux[1] + cmon[j];
                baux[1] = xbak[0]*baux[1] + xbak[1]*b0;
              }
            absb=abs(cmplx(baux[0],baux[1]));
          }
        else
          absb=abs(cmplx(b[0],b[1]));
        if (rev && (isinf(b[0]) || isnan(b[0]) || 
            isinf(b[1]) || isnan(b[1])))
          {
            //cout << "f= " << f[0] << " " << f[1]  << "d= " << d[0] << " " << d[1] << "b=" << b[0] << " " << b[1] << "\n";
            //cout << "x=" << x[0] << " " << x[1] << "\n";
            //cout << "p(x)=" << evalpoly(x) << "\n";
            flag=2;
            return false;
          }
        if (absb <= K*s) // stopping criterion of bini 
          {
            xc = cmplx(x[0],x[1]);
            //if (iter > 20)
              //cout << "lag iter=" << iter << " numimp=" << numimp << "\n";
            return true;
          }
#if 0
       if (iter > 1 && absb > absbold)
          numimp=0;
        else 
          numimp++;
        if (iter > 30)
          {
            //cout << "iter=" << iter << " numimp = " << numimp << " absb=" << absb << " absbold=" << absbold << "\n";
          }
#endif
#if 0
        if (iter > 30 && numimp < 2)
         {
           flag=1;
           return false;
         }
#endif
#else       
        if (abs(cmplx(b[0],b[1])) <= err)
          {
            //if (iter > 15)
            //cout << "iter=" << iter << "\n";
            xc = cmplx(x[0],x[1]);
            return true;
          }
#endif
        if (rev)
          {
            cmplx_div(g,d,b);
            cmplx_mul(xg,x,g); 
            xg[0] = ntype(m)-xg[0];
            xg[1] = -xg[1];
            cmplx_mul(g,x,xg);
            cmplx_div(fb,f,b);
            g2[0] = g[0]*g[0] - g[1]*g[1];
            g2[1] = 2.0*g[0]*g[1];
            h[0] = g2[0] - 2.0*fb[0];
            h[1] = g2[1] - 2.0*fb[1];
            cmplx_mul(xsq,x,x);
            cmplx_mul(xsqh,xsq,h);
            xsqh[0]+=2.0*xg[0];
            xsqh[1]+=2.0*xg[1];
            xsqh[0]+=ntype(m);
            cmplx_mul(h,xsq,xsqh);
          }
        else
          {
            cmplx_div(g,d,b);
            cmplx_div(fb,f,b);
            g2[0] = g[0]*g[0] - g[1]*g[1];
            g2[1] = 2.0*g[0]*g[1];
            h[0] = g2[0] - 2.0*fb[0];
            h[1] = g2[1] - 2.0*fb[1];
          }
        x[0] = xbak[0];
        x[1] = xbak[1];
       
        sarg[0] = ((ntype)(m-1))*( ((ntype)m) *h[0] - g2[0]);
        sarg[1] = ((ntype)(m-1))*( ((ntype)m) *h[1] - g2[1]);
        rho = sqrt(sarg[0]*sarg[0] + sarg[1]*sarg[1]);
        if (rho==0)
          {
            sq[0] = sq[1] = 0.0;
          }
        else
          {
            phi = acos(sarg[0]/rho);
            if (sarg[1] < 0.0)
              {
                phi = -phi;
              }
            rho = sqrt(rho);
            sq[0] = rho*cos(phi/2.0);
            sq[1] = rho*sin(phi/2.0);
          }
        gp[0] = g[0] + sq[0];
        gp[1] = g[1] + sq[1];
        gm[0] = g[0] - sq[0];
        gm[1] = g[1] - sq[1]; 
        abp = sqrt(gp[0]*gp[0] + gp[1]*gp[1]);
        abm = sqrt(gm[0]*gm[0] + gm[1]*gm[1]);
        if (abp < abm) 
          {
            gp[0]=gm[0];
            gp[1]=gm[1];
          }
        //cout << "gp = " << gp[0] << " " << gp[1]  << " gm=" << gm[0] << " " << gm[1] << "\n";
        invn = 1.0/(gp[0]*gp[0]+gp[1]*gp[1]);
        invgp[0] = gp[0]*invn;
        invgp[1] = -gp[1]*invn;
        //cout << "abp=" << abp << " abm=" << abm << "\n";
        //cout << "1/gp = " << invgp[0] << " " << invgp[1] << "\n";
        if (MAX(abp,abm) > 0.0)
          {
            dx[0]=((ntype)m)*invgp[0];
            dx[1]=((ntype)m)*invgp[1];
          }
        else
          {
            pol[0] = (1.0+abx)*cos((ntype)iter);
            pol[1] = (1.0+abx)*sin((ntype)iter); 
            dx[0]=pol[0];
            dx[1]=pol[1];
          }
        //cout << "dx=" << dx[0] << " " << dx[1] << "\n";
        //x1=x-dx;
        x1[0] = x[0] - dx[0];
        x1[1] = x[1] - dx[1];
        if (x[0] == x1[0] && x[1]==x1[1]) 
          {
            xc = cmplx(x[0],x[1]);
            return true;
          }
        if (iter % MT != 0) 
          {
            x[0]=x1[0];
            x[1]=x1[1];
          }
        else 
          {
            x[0] -= frac[iter/MT]*dx[0];
            x[1] -= frac[iter/MT]*dx[1];
          }
        //cout << "x=" << x[0] << "I*" << x[1] << " iter=" << iter << "\n";
      }
    //cout << "too many iterations in laguer\n";
    xc = cmplx(x[0],x[1]);
    flag=1;
    return false;
  }
#else
bool laguer_real(cmplx &xc, int& flag)
  {
    const int MR=8,MT=10,MAXIT=MT*MR;
    const ntype EPS=goaleps;
    static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    int iter,j, numimp;
    ntype K=EPS*(4.0*ntype(n)+1.0);
    ntype abx,abp,abm, phi, rho, s, absb, absbold;
#ifndef BINI_CONV_CRIT
    ntype err;
#endif
    //cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    ntype dx[2],x1[2],b[2],d[2],f[2],g[2],h[2],sq[2],gp[2],gm[2],g2[2], fb[2], pol[2];
    ntype x[2], f0, d0, b0, invn, g0, invgp[2], sarg[2];
    int m=n;
    //cout << "r0 ini = "<<xc << "\n"; 
    x[0] = xc.real();
    x[1] = xc.imag();
    flag=0;
    numimp=0;
    for (iter=1;iter<=MAXIT;iter++) 
      {
        //its=iter;
        b[0]=cmon[m];
        b[1]=0.0;
#ifndef BINI_CONV_CRIT
        err=sqrt(b[0]*b[0]+b[1]*b[1]);
#endif
        d[0]=f[0]=0.0;
        d[1]=f[1]=0.0;
        abx=sqrt(x[0]*x[0]+x[1]*x[1]);
        for (j=m-1;j>=0;j--) {
          f0 = f[0];
          f[0] = x[0]*f[0] - x[1]*f[1] + d[0];
          f[1] = x[0]*f[1] + x[1]*f0 + d[1];
          d0 = d[0];
          d[0] = x[0]*d[0] - x[1]*d[1] + b[0];
          d[1] = x[0]*d[1] + x[1]*d0 + b[1];
          b0 = b[0];
          b[0] = x[0]*b[0] - x[1]*b[1] + cmon[j];
          b[1] = x[0]*b[1] + x[1]*b0;
          // note that f is P''/2
#ifndef BINI_CONV_CRIT
          err=sqrt(b[0]*b[0]+b[1]*b[1])+abx*err;
#endif
        }
#ifndef BINI_CONV_CRIT
        err *= EPS;
#endif
        if (isinf(b[0]) || isnan(b[0]) || 
            isinf(b[1]) || isnan(b[1]))
          {
            //cout << "f= " << f[0] << " " << f[1]  << "d= " << d[0] << " " << d[1] << "b=" << b[0] << " " << b[1] << "\n";
            //cout << "x=" << x[0] << " " << x[1] << "\n";
            //cout << "p(x)=" << evalpoly(x) << "\n";
            flag=2;
            return false;
          }
#ifdef BINI_CONV_CRIT
        //if (iter > 30)
          //cout << "#" << iter << " p(x)=" << abs(cmplx(b[0],b[1])) << "\n";
        s=coeff[n];
        for (j=n-1; j >=0; j--) 
          {
            s=abx*s+coeff[j];
          }
        if (iter > 1)
          absbold=absb;
        absb=abs(cmplx(b[0],b[1]));
        if (absb <= K*s) // stopping criterion of bini 
          {
            xc = cmplx(x[0],x[1]);
            //if (iter > 20)
              //cout << "lag iter=" << iter << " numimp=" << numimp << "\n";
            return true;
          }
#if 0
       if (iter > 1 && absb > absbold)
          numimp=0;
        else 
          numimp++;
        if (iter > 30)
          {
            //cout << "iter=" << iter << " numimp = " << numimp << " absb=" << absb << " absbold=" << absbold << "\n";
          }
#endif
#if 0
        if (iter > 30 && numimp < 2)
         {
           flag=1;
           return false;
         }
#endif
#else       
        if (abs(cmplx(b[0],b[1])) <= err)
          {
            //if (iter > 15)
            //cout << "iter=" << iter << "\n";
            xc = cmplx(x[0],x[1]);
            return true;
          }
#endif
        cmplx_div(g,d,b);
        cmplx_div(fb,f,b);
        g2[0] = g[0]*g[0] - g[1]*g[1];
        g2[1] = 2.0*g[0]*g[1];
        h[0] = g2[0] - 2.0*fb[0];
        h[1] = g2[1] - 2.0*fb[1];
        sarg[0] = ((ntype)(m-1))*( ((ntype)m) *h[0] - g2[0]);
        sarg[1] = ((ntype)(m-1))*( ((ntype)m) *h[1] - g2[1]);
        rho = sqrt(sarg[0]*sarg[0] + sarg[1]*sarg[1]);
        if (rho==0)
          {
            sq[0] = sq[1] = 0.0;
          }
        else
          {
            phi = acos(sarg[0]/rho);
            if (sarg[1] < 0.0)
              {
                phi = -phi;
              }
            rho = sqrt(rho);
            sq[0] = rho*cos(phi/2.0);
            sq[1] = rho*sin(phi/2.0);
          }
        gp[0] = g[0] + sq[0];
        gp[1] = g[1] + sq[1];
        gm[0] = g[0] - sq[0];
        gm[1] = g[1] - sq[1]; 
        abp = sqrt(gp[0]*gp[0] + gp[1]*gp[1]);
        abm = sqrt(gm[0]*gm[0] + gm[1]*gm[1]);
        if (abp < abm) 
          {
            gp[0]=gm[0];
            gp[1]=gm[1];
          }
        //cout << "gp = " << gp[0] << " " << gp[1]  << " gm=" << gm[0] << " " << gm[1] << "\n";
        invn = 1.0/(gp[0]*gp[0]+gp[1]*gp[1]);
        invgp[0] = gp[0]*invn;
        invgp[1] = -gp[1]*invn;
        //cout << "abp=" << abp << " abm=" << abm << "\n";
        //cout << "1/gp = " << invgp[0] << " " << invgp[1] << "\n";
        if (MAX(abp,abm) > 0.0)
          {
            dx[0]=((ntype)m)*invgp[0];
            dx[1]=((ntype)m)*invgp[1];
          }
        else
          {
            pol[0] = (1.0+abx)*cos((ntype)iter);
            pol[1] = (1.0+abx)*sin((ntype)iter); 
            dx[0]=pol[0];
            dx[1]=pol[1];
          }
        //cout << "dx=" << dx[0] << " " << dx[1] << "\n";
        //x1=x-dx;
        x1[0] = x[0] - dx[0];
        x1[1] = x[1] - dx[1];
        if (x[0] == x1[0] && x[1]==x1[1]) 
          {
            xc = cmplx(x[0],x[1]);
            return true;
          }
        if (iter % MT != 0) 
          {
            x[0]=x1[0];
            x[1]=x1[1];
          }
        else 
          {
            x[0] -= frac[iter/MT]*dx[0];
            x[1] -= frac[iter/MT]*dx[1];
          }
        //cout << "x=" << x[0] << "I*" << x[1] << " iter=" << iter << "\n";
      }
    //cout << "too many iterations in laguer\n";
    xc = cmplx(x[0],x[1]);
    flag=1;
    return false;
  }
#endif
bool laguer_cmplx_mod(cmplx &x, pvector<cmplx,N> roots, int iac)
  {
    const int MR=8,MT=10,MAXIT=MT*MR;
    const ntype EPS=goaleps;
    static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    int iter,j, flag;
    ntype abx,abp,abm, M, s, pr, absb;
#ifndef BINI_CONV_CRIT
    ntype err;
#else
    ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
#endif
    cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    int m;
    m=n;
    //cout << "r0 ini = "<< x << "\n"; 
    flag=0;
    for (iter=1;iter<=MAXIT;iter++) 
      {
        //its=iter;
        b=cmon[m];
#ifndef BINI_CONV_CRIT
        err=abs(b);
#endif
        d=f=0.0;
        abx=abs(x);
        for (j=m-1;j>=0;j--) {
          f=x*f+d;
          d=x*d+b;
          b=x*b+cmon[j];
#ifndef BINI_CONV_CRIT
          err=abs(b)+abx*err;
#endif
        }
#ifndef BINI_CONV_CRIT
        err *= EPS;
#endif
        //cout << "f= " << f << "d= " << d << "b=" << b << "\n";
        absb=abs(b);
        if (isinf(absb) || isnan(absb))
          {
            //cout << "f= " << f << "d= " << d << "b=" << b << "\n";
            //cout << "x=" << x << "\n";
            //cout << "p(x)=" << evalpoly(x) << "\n";
            flag=2;
            return false;
          }
#ifdef BINI_CONV_CRIT
        s=0.0;
        pr=1.0; 
        for (j=0; j <= n; j ++) 
          {
            s+=coeff[j]*pr;
            pr *= abx;
          }
        if (absb <= K*s) // stopping criterion of bini 
          return true;
#else
        if (absb <= err) 
          {
            return true;
          }
#endif
        cmplx gc=0.0, hc=0.0, term;
        for (j=0; j < n; j++)
          {
            if (j==iac)
              continue;
            term=1.0/(x-roots[j]);
            gc += term;
            hc += term*term; 
          }
        g=d/b-gc;
        g2=g*g;
        h=g2-((ntype)2.0)*f/b-hc;
        sq=sqrt( ((ntype)(m-1))*( ((ntype)m)*h-g2));
        gp=g+sq;
        gm=g-sq;
        abp=abs(gp);
        abm=abs(gm);

        if (abp < abm) gp=gm;
        //cout << "abp=" << abp << " abm=" << abm << "\n";
        //cout << "1/gp=" << 1.0/gp << "\n";
        if (MAX(abp,abm) > 0.0)
          {
            dx =  ((cmplx)m)/((cmplx)gp);
          }
        else
          {
            dx = polar(((ntype)(1+abx)),((ntype)iter));
          }
        //dx=(MAX(abp,abm) > 0.0) ? ((cmplx)m)/((cmplx)gp) : polar(((ntype)(1+abx)),((ntype)iter));
        //cout << "dx=" << dx << "\n";
        x1=x-dx;
        if (x == x1) 
          {
            return true;
          }

        if (iter % MT != 0) 
          x=x1;
        else 
          x -= frac[iter/MT]*dx;
        //cout << "iter=" << iter << " x=" << x << "\n";        
      }
    //refine_root(x);
    //cout << "too many iterations in laguer\n";
    flag=1;
    return false;
  }
 bool laguer_cmplx(cmplx &x, int& flag)
    {
      const int MR=8,MT=10,MAXIT=MT*MR;
      const ntype EPS=goaleps;
      static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
      int iter,j;
      ntype abx,abp,abm, M, s, pr, absb;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
#endif
      cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
      int m;
      m=n;
      //cout << "r0 ini = "<< x << "\n"; 
      flag=0;
      for (iter=1;iter<=MAXIT;iter++) 
        {
          //its=iter;
          b=cmon[m];
#ifndef BINI_CONV_CRIT
          err=abs(b);
#endif
          d=f=0.0;
          abx=abs(x);
          for (j=m-1;j>=0;j--) {
            f=x*f+d;
            d=x*d+b;
            b=x*b+cmon[j];
#ifndef BINI_CONV_CRIT
            err=abs(b)+abx*err;
#endif
          }
#ifndef BINI_CONV_CRIT
          err *= EPS;
#endif
          //cout << "f= " << f << "d= " << d << "b=" << b << "\n";
          absb=abs(b);
          if (isinf(absb) || isnan(absb))
            {
              //cout << "f= " << f << "d= " << d << "b=" << b << "\n";
              //cout << "x=" << x << "\n";
              //cout << "p(x)=" << evalpoly(x) << "\n";
              flag=2;
              return false;
            }
#ifdef BINI_CONV_CRIT
          s=0.0;
          pr=1.0; 
          for (j=0; j <= n; j ++) 
            {
              s+=coeff[j]*pr;
              pr *= abx;
            }
          if (absb <= K*s) // stopping criterion of bini 
            return true;
#else
          if (absb <= err) 
            {
              return true;
            }
#endif
          g=d/b;
          g2=g*g;
          h=g2-((ntype)2.0)*f/b;
          sq=sqrt( ((ntype)(m-1))*( ((ntype)m)*h-g2));
          gp=g+sq;
          gm=g-sq;
          abp=abs(gp);
          abm=abs(gm);

          if (abp < abm) gp=gm;
          //cout << "abp=" << abp << " abm=" << abm << "\n";
          //cout << "1/gp=" << 1.0/gp << "\n";
          if (MAX(abp,abm) > 0.0)
            {
              dx =  ((cmplx)m)/((cmplx)gp);
            }
          else
            {
              dx = polar(((ntype)(1+abx)),((ntype)iter));
            }
          //dx=(MAX(abp,abm) > 0.0) ? ((cmplx)m)/((cmplx)gp) : polar(((ntype)(1+abx)),((ntype)iter));
          //cout << "dx=" << dx << "\n";
          x1=x-dx;
          if (x == x1) 
            {
              return true;
            }

          if (iter % MT != 0) 
            x=x1;
          else 
            x -= frac[iter/MT]*dx;
          //cout << "iter=" << iter << " x=" << x << "\n";        
        }
      //refine_root(x);
      //cout << "too many iterations in laguer\n";
      flag=1;
      return false;
    }
#if 0
 bool laguer_wrap(cmplx &x, int &flag)
   {
     if (is_same<ntype,double>::value || is_same<ntype,long double>::value)
       return laguer_cmplx(x,flag);
     else
       return laguer_real(x,flag);
   }
#endif
 void laguer(cmplx *a, cmplx &x, int m)
  {
    const int MR=8,MT=10,MAXIT=MT*MR;
    const ntype EPS=goaleps;//numeric_limits<ntype>::epsilon();
    static const ntype frac[MR+1]=
      {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    int iter,j;
    ntype abx,abp,abm,err;
    cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;

    for (iter=1;iter<=MAXIT;iter++) {
      b=a[m];
      err=abs(b);
      d=f=0.0;
      abx=abs(x);
      // note that f is P''/2 b is p(x) and d is p'(x)
      for (j=m-1;j>=0;j--) {
        f=x*f+d;
        d=x*d+b;
        b=x*b+a[j];
        err=abs(b)+abx*err;
      }
      err *= EPS;
#if 0
      if (isnan(abs(b)) || isinf(abs(b)))
        {
          cout << "ecco qua\n";
          exit(1);
        }
#endif
      if (abs(b) <= err) return;
      g=d/b;
      g2=g*g;
      h=g2-2.0*f/b;
      sq=sqrt(((ntype)(m-1))*(((ntype)m)*h-g2));
      gp=g+sq;
      gm=g-sq;
      abp=abs(gp);
      abm=abs(gm);
      if (abp < abm) gp=gm;
      dx=MAX(abp,abm) > 0.0 ? ((ntype)m)/gp : polar(1+abx,((ntype)iter));
      x1=x-dx;
      if (x == x1) return;
      if (iter % MT != 0) x=x1;
      else x -= frac[iter/MT]*dx;
    }
    //nrerror("too many iterations in laguer");
    return;
  }
#if 1
 void zroots(pvector<cmplx,N> &rootsa, const bool polish)
   {
     const ntype EPS=meps;
     int i,j,jj;
     cmplx x,b,c;
     cmplx *a;
     int m;
     cmplx *ad;
     m=n;
     ad = new cmplx[n+1];
     a = new cmplx[n+1];
     if (coeff[n]==0.0)
       {
         cout << "That's not an " << n << "degree polynomial!\n";
         delete [] a;
         delete [] ad;
         return;
       }
     cmon[n]=1.0;
     for (i=n-1; i >=0; i--)
       cmon[i]=coeff[i]/coeff[n];
     //cmon.show("bboh");
     for (j=0;j<=m;j++) 
       {
         a[j]=cmon[j];
         ad[j]=a[j];
       }
     for (j=m-1;j>=0;j--) {
       x=0.0;
       cmplx *ad_v = new cmplx[j+2];
       for (jj=0;jj<j+2;jj++) ad_v[jj]=ad[jj];
       laguer(ad_v,x,j+1);
       delete [] ad_v;
       if (fabs(imag(x)) <= 2.0*EPS*fabs(real(x)))
         x=cmplx(real(x),0.0);
       rootsa[j]=x;
       b=ad[j+1];
       for (jj=j;jj>=0;jj--) {
         c=ad[jj];
         ad[jj]=b;
         b=x*b+c;
       }
     }
     if (polish)
       for (j=0;j<m;j++)
         {
           laguer(a,rootsa[j],m);
         }
     for (j=1;j<m;j++) {
       x=rootsa[j];
       for (i=j-1;i>=0;i--) {
         if (real(rootsa[i]) <= real(x)) break;
         rootsa[i+1]=rootsa[i];
       }
       rootsa[i+1]=x;
#if 0
       if (isnan(x.real()) || isinf(x.real()))
         {
           cout << "x=" << x << "\n";
           exit(1);
         }
#endif
     }
     delete [] a;
     delete [] ad;
   }
#endif
 ntype calc_err_defl_quad(pvecNm1<ntype>& bd, pvector<ntype,3>& ad)
    {
      int i;
      ntype sum=0, term;
      for (i=0; i < n; i++)
        {
          term=cmon[i];
          if (i <= n-2) 
            term -= bd[i]*ad[0];
          if (i >= 1 && i <= n-1)
            term -= bd[i-1]*ad[1];
          if (i >= 2)
            term -= bd[i-2];
          if (cmon[i]==0)
            sum+=abs(term);
          else
            sum+=abs(term/cmon[i]);
        }
      return sum;
    }
 void refine_root_maehly(cmplx &r0, pvector<cmplx,N> &roots, int iac)
   {
     int iter,i;
     ntype err;
     cmplx povp1, denom;
     // root polishing by NR
     cmplx r0old, p, sum;
     cmplx p1, p1c;
     ntype errold;
     for (iter=0; ; iter++)
       {
         p=cmon[n]*r0+cmon[n-1];
         p1=cmon[n];
         for(i=n-2;i>=0;i--) {
           p1=p+p1*r0;
           p=cmon[i]+p*r0;
         }

         if (iter > 0)
           errold=err;

         err = abs(p);//abs(p.real())+abs(p.imag());
         if (err==0)
           break;

         if (iter > 0 && err >= errold)
           {
             r0=r0old;
             break;
           }

         if (p1==cmplx(0,0))
           {
             break;
           }
         if (iter==maxiter_polish)
           break;
         r0old=r0;
         sum=0;
         for (i=0; i < iac; i++)
           sum += 1.0/(r0-roots[i]);
         sum *= p;
         p1c = p1 - sum;
         r0 -= p/p1c;
         if (isnan(abs(r0)) || isinf(abs(r0)))
           {
             r0=r0old;
             break;
           }
       }

     //if (iter > 0)
       //cout << "NR iter=" << iter << " err=" << err << "\n";
   }
 void refine_root(cmplx &r0)
   {
     int iter,i;
     ntype err;
     cmplx povp1;
     if (r0.imag()!=0)
       {
#ifdef POLISH_NR_REAL
         // root polishing by NR
         ntype r0old[2], p[2], r0new[2];
         ntype p1[2], p0, p10;
         ntype errold, invnp1;
         cmplx pc, p1c;
         r0new[0]=r0.real();
         r0new[1]=r0.imag();
         //cout << "r0new=" << r0new[0] << " " << r0new[1] << "\n";
         for (iter=0; ;iter++)
           {
             p[0]=cmon[n]*r0new[0]+cmon[n-1];
             p[1]=cmon[n]*r0new[1];    
             p1[0]=cmon[n];
             p1[1]=0.0;
             for(i=n-2;i>=0;i--) {
               p10=p1[0];
               p1[0]=p[0]+p1[0]*r0new[0]-p1[1]*r0new[1];
               p1[1]=p[1]+p10*r0new[1]+p1[1]*r0new[0];
               p0=p[0];
               p[0]=cmon[i]+p[0]*r0new[0]-p[1]*r0new[1];
               p[1]=p0*r0new[1]+p[1]*r0new[0];
             }
             if (iter > 0)
               errold=err;
             err = abs(cmplx(p[0],p[1]));

             //cout << "iter=" << iter << " err=" << err << "\n";
             if (err==0)
               {
                 //cout << "QUI1 iter=" << iter << " err=" << err << "\n";
                 break;
               }
             if (iter > 0 && err >= errold)
               {
                 //cout << "QUI2 iter=" << iter << " err=" << err << "\n";
                 r0new[0]=r0old[0];
                 r0new[1]=r0old[1];
                 break;
               }
             if (p1[0]==0 && p1[1]==0)
               {
                 // cout << "QUI3 iter=" << iter << " err=" << err << "\n";
                 break;
               }
             if (iter == maxiter_polish)
               break;
             r0old[0]=r0new[0];
             r0old[1]=r0new[1];
#if 0
             invnp1 = 1.0/(p1[0]*p1[0]+p1[1]*p1[1]); 
             r0new[0] -= (p[1]*p1[1] + p[0]*p1[0])*invnp1;
             r0new[1] -= (p[1]*p1[0] - p[0]*p1[1])*invnp1;
#else
             povp1 = cmplx(p[0],p[1])/cmplx(p1[0],p1[1]);
             r0new[0] -= povp1.real();
             r0new[1] -= povp1.imag();
#endif
             if (isnan(r0new[0]) || isnan(r0new[1]) ||isinf(r0new[0]) || isinf(r0new[1]))
               {
                 r0new[0]=r0old[0];
                 r0new[1]=r0old[1];
                 break; 
               }
           }
         r0=cmplx(r0new[0], r0new[1]);
#else
         // root polishing by NR
         cmplx r0old, p;
         cmplx p1;
         ntype errold;
         for (iter=0; ; iter++)
           {
             p=cmon[n]*r0+cmon[n-1];
             p1=cmon[n];
             for(i=n-2;i>=0;i--) {
               p1=p+p1*r0;
               p=cmon[i]+p*r0;
             }
             if (iter > 0)
               errold=err;

             err = abs(p);//abs(p.real())+abs(p.imag());
             if (err==0)
               break;

             if (iter > 0 && err >= errold)
               {
                 r0=r0old;
                 break;
               }

             if (p1==cmplx(0,0))
               {
                 break;
               }
             if (iter==maxiter_polish)
               break;
             r0old=r0;
             r0 -= p/p1;
             if (isnan(abs(r0)) || isinf(abs(r0)))
               {
                 r0=r0old;
                 break;
               }
           }
#endif
       }
     else
       {
         // root polishing by NR
         ntype p;
         ntype p1;
         ntype r0r = r0.real();
         ntype errold, r0old;
         for (iter=0; ;iter++)
           {
             p=cmon[n]*r0.real()+cmon[n-1];
             p1 = cmon[n];
             for(i=n-2;i>=0;i--) {
               p1=p+p1*r0r;
               p=cmon[i]+p*r0r;
             }
             if (iter > 0)
               errold=err;
             err = abs(p);
             if (err==0)
               break;

             if (iter > 0 && err >= errold)
               {
                 r0r=r0old;
                 break;
               }
             if (p1==cmplx(0,0))
               break;
             if (iter==maxiter_polish)
               break;
             r0old=r0r;
             r0r -= p/p1;
             if (isnan(r0r)|| isinf(r0r))
               {
                 r0r=r0old;
                 break;
               }
           }

         r0 = cmplx(r0r,0);
       }

     //if (iter > 0)
       //cout << "NR iter=" << iter << " err=" << err << "\n";
   }
 void rescale_coeff(ntype xi)
   {
     ntype pr=xi;
     int i;
     for (i=n-1; i >=0; i--)
       {
         cmon[i] /= pr;
         pr *= xi;
       }
   }
 ntype calc_upper_bound_kal(void)
   {
     // Kalantari's formula as found in McNamee Pan Vol. 1
     int k;
     static const ntype K= 1.0/0.682338;
     ntype ximax, xi;
     ntype cnsq, term;

     cnsq=cmon[n-1]*cmon[n-1];
     for (k=4; k <= n+3; k++)
       {
         term=0;
         if (n-k+3 >= 0)
           term += cnsq*cmon[n-k+3]-cmon[n-2]*cmon[n-k+3];
         if (n-k+2 >= 0)
           term +=-cmon[n-1]*cmon[n-k+2];
         if (n-k+1 >= 0)
           term += cmon[n-k+1];
         xi = pow(abs(term),1.0/(k-1));
         if (k==4 || xi > ximax)
           ximax=xi;
       }
     ximax=K*ximax;
     if (isinf(ximax)||isnan(ximax))
       {
         //cout << "I cannot calculate the upper bound...\n";
         //cmon.show("coeff");     
         //exit(-1);
         return pow(maxf/1.618034,1.0/n);
       }
           
     return ximax;
   }
 ntype calc_upper_bound_mig(void)
   {
     // Mignotte et al 2017
     int i, imax[2], h;
     ntype CA, CB, ai, maxai[2], maxai2, R, rho; 
     //maxai[0] is the largest and maxai[1] is the second largest
     maxai[0]=pow(abs(cmon[0]),1.0/n); 
     imax[0] = n;
     maxai[1]=pow(abs(cmon[1]),1.0/(n-1));
     imax[1] = n-1;
     if (maxai[1] > maxai[0])
       {
         swap(maxai[0],maxai[1]);
         swap(imax[0], imax[1]);
       }
     for (i=n-3; i >= 0; i--)
       {
         ai=pow(abs(cmon[n-i-1]),1.0/(i+1));
         if (ai > maxai[0])
           {
             maxai[0]=ai;
             imax[0]=i+1;
           }
        else if (ai > maxai[1])
           {
             maxai[1]=ai;
             imax[1]=i+1;
           }
       }
     if (imax[0]==imax[1] || imax[0]==0 || imax[1]==0)
       {
         cout << "non puo' essere!\n!";
         exit(-1);
       }
     R=maxai[0];
     rho=maxai[1];
     //cout << "R= " << R << " rho=" << rho << "\n";
     //cout << "imax0=" << imax[0] << "imax1=" << imax[1] << "\n";
     ntype arg=R*R - 2.0*rho*R + 5.0*rho*rho;
     CB=R+rho;
     if (arg>=0)
       {
         CB += sqrt(arg);
       }
     CB /= 2.0;
     if (imax[0] == 1)
      return CB;
     CA=0.0;
     h=imax[0];
     // calcolare in maniera più efficiente quanto segue
#if 0
     ntype *rhoarr = new ntype[h];
     ntype pr;
     pr=1.0;
     for (i=h-1; i >= 0; i--)
       {
         rhoarr[i] = pr;
         pr *= rho;
       }
      CA=0;
      for (i=h-1; i >= 0; i--)
        {
          CA = rhoarr[i] + CA*R;
        }
 
    delete [] rhoarr; 
#else
     for (i=h-1; i >= 0; i--)
       { 
         CA +=pow(R,i)*pow(rho,h-1-i);
       }
#endif
     if (isinf(CA)||isnan(CA))
       {
         // close to our actual infinite...
         CA=maxf/1.618034;
       }
     //cout << "Ph=" << h << " CA=" << CA << "\n";
     CA=pow(CA,1.0/(h-1));
     //cout << "h=" << h << " CA=" << CA << " R=" << R << " rho=" << rho << "\n";
    if (CA > CB)
       return CA;
     else 
       return CB;
   }
#ifdef BINI_INIGUESS
#include<vector>
#include<set>
#include<utility>

 // A divide and conquer program to find convex 
 // hull of a given set of points. 
 // stores the centre of polygon (It is made 
 // global because it is used in compare function) 
 using point=pair<int,ntype>; 
 point<ntype> mid; 
 // determines the quadrant of a point 
 // (used in compare()) 
 int quad(point p) 
   { 
     if (p.first >= 0 && p.second >= 0) 
       return 1; 
     if (p.first <= 0 && p.second >= 0) 
       return 2; 
     if (p.first <= 0 && p.second <= 0) 
       return 3; 
     return 4; 
   } 

 // Checks whether the line is crossing the polygon 
 int orientation(point a, point b, 
                 point c) 
   { 
     int res = (b.second-a.second)*(c.first-b.first) - 
       (c.second-b.second)*(b.first-a.first); 

     if (res == 0) 
       return 0; 
     if (res > 0) 
       return 1; 
     return -1; 
   } 

 // compare function for sorting 
 bool compare(point p1, point q1) 
   { 
     point p = make_pair(p1.first - mid.first, 
                         p1.second - mid.second); 
     point q = make_pair(q1.first - mid.first, 
                         q1.second - mid.second); 

     int one = quad(p); 
     int two = quad(q); 

     if (one != two) 
       return (one < two); 
     return (p.second*q.first < q.second*p.first); 
   } 

 // Finds upper tangent of two polygons 'a' and 'b' 
 // represented as two vectors. 
 vector<point> merger(vector<point > a, vector<point > b) 
   { 
     // n1 -> number of points in polygon a 
     // n2 -> number of points in polygon b 
     int n1 = a.size(), n2 = b.size(); 

     int ia = 0, ib = 0; 
     for (int i=1; i<n1; i++) 
       if (a[i].first > a[ia].first) 
         ia = i; 

     // ib -> leftmost point of b 
     for (int i=1; i<n2; i++) 
       if (b[i].first < b[ib].first) 
         ib=i; 

     // finding the upper tangent 
     int inda = ia, indb = ib; 
     bool done = 0; 
     while (!done) 
       { 
         done = 1; 
         while (orientation(b[indb], a[inda], a[(inda+1)%n1]) >=0) 
           inda = (inda + 1) % n1; 

         while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) <=0) 
           { 
             indb = (n2+indb-1)%n2; 
             done = 0; 
           } 
       } 

     int uppera = inda, upperb = indb; 
     inda = ia, indb=ib; 
     done = 0; 
     int g = 0; 
     while (!done)//finding the lower tangent 
       { 
         done = 1; 
         while (orientation(a[inda], b[indb], b[(indb+1)%n2])>=0) 
           indb=(indb+1)%n2; 

         while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1])<=0) 
           { 
             inda=(n1+inda-1)%n1; 
             done=0; 
           } 
       } 

     int lowera = inda, lowerb = indb; 
     vector<point> ret; 

     //ret contains the convex hull after merging the two convex hulls 
     //with the points sorted in anti-clockwise order 
     int ind = uppera; 
     ret.push_back(a[uppera]); 
     while (ind != lowera) 
       { 
         ind = (ind+1)%n1; 
         ret.push_back(a[ind]); 
       } 

     ind = lowerb; 
     ret.push_back(b[lowerb]); 
     while (ind != upperb) 
       { 
         ind = (ind+1)%n2; 
         ret.push_back(b[ind]); 
       } 
     return ret; 

   } 

 // Brute force algorithm to find convex hull for a set 
 // of less than 6 points 
 vector<point> bruteHull(vector<point> a) 
   { 
     // Take any pair of points from the set and check 
     // whether it is the edge of the convex hull or not. 
     // if all the remaining points are on the same side 
     // of the line then the line is the edge of convex 
     // hull otherwise not 
     set <point >s; 

     for (int i=0; i<a.size(); i++) 
       { 
         for (int j=i+1; j<a.size(); j++) 
           { 
             int x1 = a[i].first, x2 = a[j].first; 
             int y1 = a[i].second, y2 = a[j].second; 

             int a1 = y1-y2; 
             int b1 = x2-x1; 
             int c1 = x1*y2-y1*x2; 
             int pos = 0, neg = 0; 
             for (int k=0; k<a.size(); k++) 
               { 
                 if (a1*a[k].first+b1*a[k].second+c1 <= 0) 
                   neg++; 
                 if (a1*a[k].first+b1*a[k].second+c1 >= 0) 
                   pos++; 
               } 
             if (pos == a.size() || neg == a.size()) 
               { 
                 s.insert(a[i]); 
                 s.insert(a[j]); 
               } 
           } 
       } 

     vector<point>ret; 
     for (auto e:s) 
       ret.push_back(e); 

     // Sorting the points in the anti-clockwise order 
     mid = {0, 0}; 
     int n = ret.size(); 
     for (int i=0; i<n; i++) 
       { 
         mid.first += ret[i].first; 
         mid.second += ret[i].second; 
         ret[i].first *= n; 
         ret[i].second *= n; 
       } 
     sort(ret.begin(), ret.end(), compare); 
     for (int i=0; i<n; i++) 
       ret[i] = make_pair(ret[i].first/n, ret[i].second/n); 

     return ret; 
   } 

 // Returns the convex hull for the given set of points 
 vector<point> divide(vector<point> a) 
   { 
     // If the number of points is less than 6 then the 
     // function uses the brute algorithm to find the 
     // convex hull 
     if (a.size() <= 5) 
       return bruteHull(a); 

     // left contains the left half points 
     // right contains the right half points 
     vector<point>left, right; 
     for (int i=0; i<a.size()/2; i++) 
       left.push_back(a[i]); 
     for (int i=a.size()/2; i<a.size(); i++) 
       right.push_back(a[i]); 

     // convex hull for the left and right sets 
     vector<point>left_hull = divide(left); 
     vector<point>right_hull = divide(right); 

     // merging the convex hulls 
     return merger(left_hull, right_hull); 
   } 
 void calc_upper_convex_hull(pvector <ntype,-1> pts, int& q, pvector<int,-1>& k)
   {
     vector<point > a; 
     int i;
     for (i=0; i < n; i++) 
       a.push_back(make_pair(int(pts[0]), pts[1])); 
     q = a.size(); 
     k.allocate(q);
     // sorting the set of points according 
     // to the x-coordinate 
     sort(a.begin(), a.end()); 
     vector<point >ans = divide(a); 

#if 1
    // cout << "convex hull:\n"; 
     i=0;
     for (auto e:ans) 
       {
         k[i]=e.first;
         i++;
       }
       //cout << e.first << " "
         //<< e.second << endl; 
#endif
   }
 void initial_guess(vector<cmplx>& rg)
   {
     pvector<ntype,-1> pts;
     pvector<int,-1> k;
     pts.allocate(n);
     k.allocate(n+1);//#verts will be < n
     int q, i, j;
     for (i = 0; i < n; i++)
       {
         pts[i][0] = i;
         pts[i][1] = log(abs(cmon[i]));
       }
     calc_upper_convex_hull(pts, q, k);
     pvector<ntype,-1> uk(q-1);
     for (i=1; i <= q-1; i++) 
       uk[i] = pow(cmon[k[i]]/cmon[k[i-1]],1.0/(k[i+1]-k[i]));

     for (i=1; i <= q-1; i++)
       {
         for (j=0; j < k[i+1]-k[i]; j++)
           {
             rg[k[i]+j] = polar(uk[i+1],j*2.0*M_PI/(k[i+1]-k[i])+2.0*M_PI*ntype(i)/ntype(n)+2.0*M_PI*(drand48()-0.5));
           }
       }
     // sort in descending order
     std::sort(rg.begin(),rg.end(), [&] (ntype a, ntype b)-> bool {return abs(a) > abs(b);});
   }
#endif
 void laguerre_defl(pvector<cmplx,N>& roots,bool polish=false, bool backup=false)
   {
     // find initial guess for x here
     cmplx r0;
     const ntype goldfact=1.618034;
     ntype sig, sigold, xi;
     //bool mulsig=true;
     ntype fact, err0,err1, errA, errB;
     pvector<ntype,3> ad;
     typename std::conditional<(N < 0),pvector<ntype>,bool>::type cbak;
     int iter, i, flag;
#if 0
      if constexpr (N < 0)
        {
          bd.allocate(n-1);
          bd0.allocate(n-1);
          bd1.allocate(n-1);
          // ro.allocate(n-1);
          deflcoeff0.allocate(n);
          deflcoeff1.allocate(n);
          deflcoeff.allocate(n);
        }
#endif

      //cout << "coeff[" << n << "]=" << coeff[n] << "\n"; 
      if (coeff[n]==0.0)
        {
          cout << "That's not an " << n << "degree polynomial!\n";
          return;
        }

      cmon[n]=1.0;
      for (i=n-1; i >=0; i--)
        cmon[i]=coeff[i]/coeff[n];
      if (cmon[0]==0)
        {
          if constexpr (N < 0)
            {
              cbak.allocate(n+1);
              if (backup||polish)
                {
                  cbak = coeff;
                }
            }

          roots[0]=0.0;
          for (i=n; i >=1; i--)
            {
              deflcoeff[i-1] = cmon[i];
            }
          rpolyNm1 qp;
          // pvector<cmplx,N-1> qr;
          pvecNm1<cmplx> qr;

          if constexpr (N < 0)
            {
              //qp.allocate(n-1);
              //qr.allocate(n-1);
              qr.use_vec(n-1,&(roots[1])); 
              qp.use_vec(n-1,coeff.v, cmon.v, bd0.v, bd1.v, bd.v, deflcoeff0.v, deflcoeff1.v, deflcoeff.v);
            }

          //refine deflated quartic coefficients thourh Newton-Raphson method
          qp.set_output_prec(this->get_output_prec());
          qp.set_coeff(deflcoeff);

          qp.find_roots(qr, false, false);
          if constexpr (N > 0)
            {
              for (i=1; i < n; i++)
                roots[i] = qr[i-1];
            }
        }
      else 
        {
          iter = 0;
          // find first zero using laguerre
          // if it does not converge generate a random initial guess
          // and retry (a limited number of times though!)
          //cout << "mig=" << calc_upper_bound() << " kal=" << calc_upper_bound_kal() << "\n";
#if 1
          if constexpr (N < 0)
            {
              cbak.allocate(n+1);
              if (backup||polish)
                {
                  cbak = coeff;
                }
            }
#endif
#ifdef BINI_INIGUESS
          if (!deflated)
            {
              vector<cmplx> rg;
              rg.resize(n);
              initial_guess(rg);
              for (i=0; i < n; i++)
                iniguess=rg[i];
            }
          r0 = rg[n-1];//guess with smallest magnitude 
#else
          r0=0.0;
#endif
          //absolute value of coefficients are used in Bini's stopping criterion
          //we use coeff vector to store them since it won't be used anymore from here on
          for (i=0; i <= n; i++) 
            coeff[i] = abs(cmon[i]);   

          int MAXSTPS=1000;// MT=100;
          sigold=0.0;
          sig=1.0;

          //cout << "Upper bound xi=" << xi << "\n";
          bool rescaled=false;
          //cout << "scfact=" << scfact << "0.1/maxf=" << 0.1/maxf << " n=" << n << "\n";
          while (!laguer_real(r0,flag))
            {
              // if we got an overflow rescale coefficients
              if (iter==0 && flag==2) // flag==2 means inf or nan
                {
                  xi=calc_upper_bound_kal();
                  //cout << "xi^n=" << pow(xi,n) << "\n";
                  if (xi > 1.0)
                    {
                      //cout << "qui rescaling" << "n=" << n << "\n";
                      rescaled=true;
                      rescale_coeff(xi);
                      fact=1.0;
                    }
                  else 
                    fact=xi;
                }
              else
                {
                  fact=calc_upper_bound_kal();
                }
              if (iter==0)
                {
                  if (rescaled)
                    r0=0.0;
                  else
                    r0=cmplx(0.0,1.0);
                }
              else
                {
                  r0 = polar(fact*drand48(),ntype(2.0*M_PI*(drand48()-0.5)));
                }
              //r0 = polar(ntype(drand48()),ntype(2.0*M_PI*(drand48()-0.5)));
              //r0 = polar(sigold+(sig-sigold)*drand48(),ntype(2.0*M_PI*(drand48()-0.5)));
              iter++;
              //cout << "r0=" << r0 << "xi=" << xi << "\n";
#if 0
              if (iter%MT==0)
                {
                  sigold=sig;
                  sig*=goldfact;
                }
#endif
              if (iter > MAXSTPS)
                {
                  cout << "[ERROR] Laguerre method did not converge:\n";
                  cout << "this problem can be fixed by using floating point number\n";
                  cout << "with higher precision (e.g. long double and complex<long double>\n";
                  cout << "or reducing the convergence precision\n";
                  cout << "xi=" << xi << "\n";
                  exit(-1);
                }
            };
#if 0
        if (iter>0)
          {
            cout << "flag=" << flag << " iter=" << iter << "\n";
          }
#endif
#ifdef LAGUERRE_POLISH
          if (is_same<ntype,double>::value || is_same<ntype,long double>::value)
            {
              refine_root(r0);
            }
#endif
          if (//abs(evalpoly(r0.real())) <= abs(evalpoly(r0)) //&& abs(evalpoly(r0.real())) <= abs(evalpoly(conj(r0)))) 
            abs(r0.imag()) <= meps*abs(r0))
              {
                rpolyNm1 qp;
                // pvector<cmplx,N-1> qr;
                pvecNm1<cmplx> qr;
#if 0

                if constexpr (N < 0)
                  {
                    //qp.allocate(n-1);
                    //qr.allocate(n-1);
                    qr.use_vec(n-1,&(roots[1])); 
                    qp.use_vec(n-1,coeff.v, cmon.v, cmonder.v, cmondder.v);
                  }
#endif
                //cout << "minus one r0=" << r0 << "\n";
                roots[0] = r0.real();
                for (i=n; i >=1; i--)
                  {
                    deflcoeff0[i-1] = cmon[i];
                    if (i <= n-1)
                      deflcoeff0[i-1] += r0.real()*deflcoeff0[i];
                  }
#ifdef ENABLE_FORW_DEFL
                if (r0.real()!=0)
                  {
                    //forward deflation
#if 0
                    deflcoeff1[n-1]=1.0;
                    for (i=0; i < n-1; i++)
                      {
                        deflcoeff1[i] = cmon[i]; 
                        if (i >= 1)
                          deflcoeff1[i] -= deflcoeff1[i-1];
                        deflcoeff1[i] /= -r0.real();
                        //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff1[" << i << "]=" << deflcoeff1[i]<< "\n";
                      }
#else
                    deflcoeff1[0]=-cmon[0]/r0.real();
                    deflcoeff1[n-1] = 1.0;
                    for (i=1; i < n-1; i++)
                      {
                        deflcoeff1[i] = (deflcoeff1[i-1] - cmon[i])/r0.real();
                      }
#endif
                    // choose the set which yields smaller error 
                    err0=calc_err_q(deflcoeff0,r0.real());
                    err1=calc_err_q(deflcoeff1,r0.real());

                    if (isnan(err1) || isinf(err1))
                      deflcoeff=deflcoeff0;
                    else if (isnan(err0) || isinf(err0))
                      deflcoeff=deflcoeff1;
                    else if (err0 < err1)
                      {
                        //errA = err0;
                        deflcoeff = deflcoeff0;
                      }
                    else
                      {
                        //errA = err1;
                        deflcoeff = deflcoeff1;
                      }
                  }
                else
                  {
                    //errA=calc_err_q(deflcoeff0,r0.real());
                    deflcoeff = deflcoeff0;
                  }
#else
                deflcoeff = deflcoeff0;
#endif
#ifdef USE_DEFL_NR
                //ntype errP, errD;
                //errP = calc_err_q(deflcoeff, r0.real());
                use_NR_defl_lin(deflcoeff,r0);
                //errD = calc_err_q(deflcoeff, r0.real());
#endif

                if constexpr (N < 0)
                  {
                    //qp.allocate(n-1);
                    //qr.allocate(n-1);
                    qr.use_vec(n-1,&(roots[1])); 
                    qp.use_vec(n-1,coeff.v, cmon.v, bd0.v, bd1.v, bd.v, deflcoeff0.v, deflcoeff1.v, deflcoeff.v);
                  }
#ifdef BINI_INIGUESS
                qp.set_iniguess(iniguess);
                qp.set_deflated();
#endif
                //refine deflated quartic coefficients thourh Newton-Raphson method
                qp.set_output_prec(this->get_output_prec());
                qp.set_coeff(deflcoeff);

                qp.find_roots(qr, false, false);
                if constexpr (N > 0)
                  {
                    for (i=1; i < n; i++)
                      roots[i] = qr[i-1];
                  }
              }
          else
            {
              rpolyNm2 rp;
              //pvector<cmplx,N-2> ro_quad;
              pvecNm2<cmplx> ro_quad;
#if 0
              if constexpr (N < 0)
                {
                  ro_quad.allocate(n-2);
                  rp.allocate(n-2);
                }
#endif
              //cout << "2) r0=" << r0 << "\n";
              roots[0] = r0;
              //cout << "r0=" << r0 << " p(real(r0))=" << evalpoly(r0.real()) << " p(r0)=" << evalpoly(r0)<< "\n";

              roots[1] = conj(r0);
              // find coefficients of quadratic polynomial for deflation

              ad[2] = 1.0;
              ad[1] = -2*r0.real();
              ad[0] = abs(r0)*abs(r0);
              //do forward and backward deflation and then choose the one which provides smaller error
              //backward deflation 
              //cout << "minus two r0=" << r0 << "\n";
              for (i=n; i >=2;  i--)
                {
                  bd0[i-2] = cmon[i];

                  if (i <= n-1)
                    bd0[i-2] -= ad[1]*bd0[i-1];
                  if (i <= n-2)
                    bd0[i-2] -= ad[0]*bd0[i];

                }
#ifdef ENABLE_FORW_DEFL
              if (ad[0]!=0)
                {
                  //forward deflation
                  bd1[n-2]=1.0;
                  for (i=0; i < n-2; i++)
                    {
                      bd1[i] = cmon[i];
                      if (i >= 1)
                        bd1[i] -= bd1[i-1]*ad[1];
                      if (i >= 2)
                        bd1[i] -= bd1[i-2];
                      bd1[i] /= ad[0];
                    }
                  err0=calc_err_defl_quad(bd0, ad);
                  err1=calc_err_defl_quad(bd1, ad);

                  //cout << "AD err0=" << err0 << " err1=" << err1 << "\n";
                  if (isnan(err1) || isinf(err1))
                    bd=bd0;
                  else if (isnan(err0) || isinf(err0))
                    bd=bd1;
                  else if (err0 < err1)
                    {
                      bd = bd0;
                      //errB = err0;
                    }
                  else
                    {
                      bd = bd1;
                      //errB = err1;
                    }
                }
              else
                {
                  //errB=calc_err_defl_quad(bd0, ad);
                  bd = bd0;
                }
#else
              bd=bd0;
#endif
              //cout << "N=" << n << "\n";
#ifdef USE_DEFL_NR
              //ntype errP, errD;
              //errP =calc_err_defl_quad(bd, ad);
              use_NR_defl_quad(bd,ad);
#endif
              if constexpr (N < 0)
                {
                  // reuse already allocated vectors
                  rp.use_vec(n-2,coeff.v,cmon.v, bd0.v, bd1.v, bd.v, deflcoeff0.v, deflcoeff1.v, deflcoeff.v);
                  ro_quad.use_vec(n-2,&(roots[2]));
                }
              // find zeroes of N-2 degree polynomial
#ifdef BINI_INIGUESS
              rp.set_iniguess(iniguess);
              rp.set_deflated();
#endif
              rp.set_output_prec(this->get_output_prec());
              rp.set_coeff(bd);
              rp.find_roots(ro_quad, false, false);
              //ro_quad.show("quartic roots=");
              if constexpr (N > 0)
                {
                  for (i=2; i < n; i++)
                    roots[i] = ro_quad[i-2];
                }
            }
          if (rescaled)
            {
              for (i=0; i < n; i++)
                roots[i] *= xi;
#if 0
              for (i=0; i < n; i++)
                {
                  cout << "ratio=" <<  abs(ub/roots[i]) << "\n"; 
                } 
#endif
            }
        }
      if constexpr (N < 0)
        {
          if (backup||polish)
            {
              coeff = cbak;
            }
        }
      bool copyroot=false;
      if (polish==true)
        {
          if constexpr (N < 0)
            {
              cmon[n]=1.0;
              for (i=n-1; i >=0; i--)
                cmon[i]=coeff[i]/coeff[n];
            }
          for (i=0; i < n; i++)
            {
              // refine conjuigate pairs only once!
              if (i < n-1 && roots[i] == conj(roots[i+1]))
                copyroot=true;
              else 
                copyroot=false;
              refine_root_maehly(roots[i], roots, i);
              if (copyroot==true)
                {
                  roots[i+1] = conj(roots[i]);
                  i++;
                }
            }
        }
   } 
#if 0
 void laguerre_defl_lag_parall()
   {
     bool *found = new bool[n];
     for (i=0;i < n; i++)
       found[i]=false;
     int nf=0;
     for (iter=0; iter < 1000; iter++)
       {
         for (i=0; i < n; i++)
           {
             // refine conjuigate pairs only once!
             //if (i < n-1 && roots[i] == conj(roots[i+1]))
             //copyroot=true;
             //else 
             //copyroot=false;
             //refine_root_maehly(roots[i], roots, i);
             if (found[i]) 
               continue;
             if (laguer_cmplx_mod(roots[i], roots, i))
               {
                 found[i]=true;
                 nf++;
                 if (nf==n)
                   break;
               }
             if (copyroot==true)
               {
                 //roots[i+1] = conj(roots[i]);
                 //i++;
               }
           }
       }
     cout << "# found=" << nf << "\n";
     delete []found;

   }
#endif
  void odd_degree(pvector<cmplx,N>& roots)
    {
      ntype *coeffder;
      rpolyNm1 quar;
      //pvector<cmplx,N-1> qroots;
      pvecNm1<cmplx> qroots;
      //pvector<ntype,N> deflcoeff0, deflcoeff1, deflcoeff;
      cmplx quadr[2];
      stlarr<ntype,N> rqroots;
      ntype goldfact = 1.618034;
      ntype delx, c;
      pvector<ntype, N> cmonder(N);
      pvecNm1<ntype> cmondder(N-1);
 
      int i, cc;
      if constexpr (N < 0)
        {
          qroots.allocate(n-1);
          //deflcoeff0.allocate(n);
          //deflcoeff1.allocate(n);
          //deflcoeff.allocate(n);
          quar.allocate(n-1);
          rqroots.resize(n-1);
        }

      if (coeff[n]==0.0)
        {
          cout << "That's not an " << n  << " degree polynomial!\n";
          return;
        }
      cmon[n]=1.0;
      for (i=n-1; i >=0; i--)
        cmon[i]=coeff[i]/coeff[n];

      cmonder[n-1] = n;
      for (i = n-2; i >= 0; i--)
        {
          cmonder[i] = (i+1)*cmon[i+1]; 
        }
      quar.set_coeff(cmonder);

      quar.find_roots(qroots);
      //quar.show("quartic=");
      //qroots.show("radici quartica");
      cc=0;
      for (i=0; i < n-1; i++)
        if (qroots[i].imag()==0.0 && evalddpoly(qroots[i].real())!=0)
          rqroots[cc++] = qroots[i].real();
      int numr = cc;
      //for (i=0; i < numr; i++)
        //cout << "PRIMA rqr[" << i <<  "]" << rqroots[i]<< "\n";
      std::sort(rqroots.begin(),rqroots.end(), [&] (ntype a, ntype b)-> bool {return a < b;});
      //for (i=0; i < numr; i++)
        //cout << "DOPO rqr[" << i <<  "]" << rqroots[i] << "\n";
    
      //cout << "----> "<< rqroots[0] << " "<< rqroots[1] << "\n";
      bool foundleft = false;
      bool foundright = false;
      ntype xleft, xright;
      if (numr == 0)
        {
          rpolyNm2 cub;
          ntype r0;
          //pvector<cmplx,N-2> cubroots;
          pvecNm2<cmplx> cubroots;
          if constexpr (N < 0)
            cubroots.allocate(n-2);

          cmondder[n-2] = n*(n-1);
          //cout << "qui\n";
          for (i = n-3; i >= 0; i--)
            {
              cmondder[i] = (i+1)*(i+2)*cmon[i+2]; 
            }
          cub.set_coeff(cmondder);
          //cmondder.show();
          cub.find_roots(cubroots); 
          for (i=0; i < n-2; i++)
            {
              if (cubroots[i].imag()==0)
                {
                  r0=cubroots[i].real();
                  break;
                }
            } 
          //cout << "real cubic root=" << r0 << "val=" << evalpoly(r0) << "\n";
          if ((c=evalpoly(r0)) > 0.0)
            {
              xright = r0;
              // quadratic extrapolation
              xleft = xright-pow(c,1.0/n);
              while (evalpoly(xleft) >=0)
                {
                  delx = xleft - xright;
                  delx *= goldfact;
                  xleft = xright+delx;
                };
              //cout << "QUIQUI-1xleft=" << xleft << " xright=" << xright << "\n";
            }
          else
            {
              xleft = r0;
              c = evalpoly(xleft);///evalddpoly(xleft);
              //cout << evalpoly(xleft) << "dd=" << evalddpoly(xleft)<<"\n"; 
              xright = xleft+pow(-c,1.0/n);
              while (evalpoly(xright) <=0)
                {
                  delx = xright - r0;
                  delx *= goldfact;
                  xright = r0+delx;
                };
              //cout << "QUIQUI0 xleft=" << xleft << " xright=" << xright << "\n";
            }
        }
      else
        {
          if ((c=evalpoly(rqroots[0])) > 0.0)
            {
              foundright = true;
              xright = rqroots[0];
              // find xleft here
              // quadratic extrapolation
              xleft = xright-pow(c,1.0/n);
              while (evalpoly(xleft) >=0)
                {
                  delx = xleft - xright;
                  delx *= goldfact;
                  xleft = xright+delx;
                };
              foundleft = true;
              //cout << "QUIQUI2 xleft=" << xleft << " xright=" << xright << "\n";
            }
          else if ((c=evalpoly(rqroots[numr-1])) < 0.0)  
            {
              foundleft = true;
              xleft = rqroots[numr-1];
              //cout << "xleft=" << xleft << "\n";
              // find xright here
              //cout << "c=" << c << "\n";
              xright = xleft+pow(-c,1.0/n);
              //cout << "xright temp=" << xright << "\n";
              cc=0;
              while (evalpoly(xright) <=0)
                {
                  delx = xright - xleft;
                  delx *= goldfact;
                  xright = xleft+delx;
                  cc++;
                  //if (cc > 10)
                  // cout << "actung actung cc="<< cc <<"\n"; 
                };

              //cout << "QUIQUI3 xleft=" << xleft << " xright=" << xright << "\n";
              //cout << "QUIQUI3 P(xl)=" << evalpoly(xleft) << "P(xr)=" << evalpoly(xright)<< "\n";
              //cout << "xright=" << xright << "\n";
              foundright = true;
            }
          if (foundleft == false && foundright == false)
            {
              for (i=0; i < numr; i++)
                {
                  if (evalpoly(rqroots[i])*evalpoly(rqroots[i+1]) < 0.0)
                    {
                      xleft = rqroots[i];
                      xright = rqroots[i+1];
                      foundleft = true;
                      foundright = true;
                    }
                }
              if (foundright==false && foundleft==false)
                {
                  cout << "num=" << numr << "\n";
                  cout << "QUIQUI4 P(xl)=" << evalpoly(xleft) << "P(xr)=" << evalpoly(xright)<< "\n";
                  this->show("Prob quintic=");
                  exit(-1);
                }
            }
        } 
      // first root
      //cout << "BOH xleft=" << xleft << " xright=" << xright << "\n";
#ifdef USE_BRENT
      ntype r0 = zbrent(xleft, xright, numeric_limits<ntype>::epsilon());
#else
      ntype r0 = rtsafe(xleft, xright, 1E-8);//1E-8 100*numeric_limits<ntype>::epsilon());
#endif
      roots[0] = r0;
      //cout << "quintic r0=" << r0 << "\n";
      // deflate quintic and use OQS
      // first set of coefficients
      deflcoeff0[n-1]=1.0;
      //cout << "cmon=" << cmon[N-1] << "\n";
      deflcoeff0[n-2] = cmon[n-1]+r0;
      //cout << "deflcoeff[" << N-1 << "]=" << deflcoeff[N-1]<< "\n";
      for (i=n-3; i >=0; i--)
        {
          //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff[" << i+1 << "]=" << deflcoeff[i+1]<< "\n";
          deflcoeff0[i] = cmon[i+1] + r0*deflcoeff0[i+1];
        }
      if (r0!=0)
        {
          deflcoeff1[0]=-cmon[0]/r0;
          deflcoeff1[n-1] = 1.0;
          //cout << "cmon[0]=" << cmon[0] << " deflcoeff1[0]=" << deflcoeff1[0]<< "\n";
          for (i=1; i < n-1; i++)
            {
              deflcoeff1[i] = (deflcoeff1[i-1] - cmon[i])/r0;
              //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff1[" << i << "]=" << deflcoeff1[i]<< "\n";
            }
          ntype err0,err1;
          // choose the set which yields smaller error 
          err0=calc_err_q(deflcoeff0,r0);
          err1=calc_err_q(deflcoeff1,r0);
          if (err0 < err1)
            deflcoeff = deflcoeff0;
          else
            deflcoeff = deflcoeff1;
          //cout << "err0= " << err0 << " err1=" << err1 << "\n";
        }
      else
        {
          //cout << "err0=" << calc_err_q(deflcoeff0, r0) << "\n";
          deflcoeff = deflcoeff0;
        }
      // refine deflated quartic coefficients here
      //use_NR_defl_lin(deflcoeff,r0);

      // calculate roots of quartic
      quar.set_coeff(deflcoeff);
      quar.find_roots(qroots);
      //quar.show("deflq=");
      //qroots.show("quartic roots=");
      for (i=1; i < n; i++)
        {
          roots[i] = qroots[i-1];
        }
    }
  pvecNp1<ntype> get_coeff()
    {
      return coeff;
    }
  int degree()
    {
      return n; 
    }
  inline void find_roots(pvector<cmplx,N>& roots, bool polish=false, bool backup=true)
    {
      // if backup is true the dynamic version preserve the coefficients after calling this method

      if constexpr (N < 0)
        {
          if (n==1)
            {
              cout << "What?!? You are not able to solve a linear equation come-on!";
              exit(-1);
            }
          else if (n==2)
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  solve_quadratic(roots);
                }
            }
          else if (n==3) 
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  solve_cubic_analytic(roots);
                }
            }
          else if (n==4)
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  oqs_quartic_solver(roots);
                }
            }
          else if (n==5)
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
#ifdef USE_LAGUERRE_ONLY
                  laguerre_defl(roots, polish, backup);
#else
                  odd_degree(roots);
#endif
                }
            }
          else 
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  laguerre_defl(roots, polish, backup);
                }  
            }
        }
      else
        {
          if constexpr (N==1)
            {
              cout << "What?!? You are not able to solve a linear equation come-on!";
              exit(-1);
            }
          else if constexpr (N==2)
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  solve_quadratic(roots);
                }
            }
          else if constexpr (N==3) 
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  solve_cubic_analytic(roots);
                }
            }
          else if constexpr (N==4)
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
                  oqs_quartic_solver(roots);
                }
            }
          else if constexpr (N==5)
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
#ifdef USE_LAGUERRE_ONLY
                  laguerre_defl(roots, polish);
#else
                  odd_degree(roots);
#endif
                }
            }
          else 
            {
              if constexpr (force_hqr)
                {
                  int ok=1;
                  solve_numrec(roots, ok);
                }
              else
                {
#ifdef USE_LAGUERRE_ONLY
                  laguerre_defl(roots, polish);
#else
                  if constexpr (N % 2 == 0)
                    laguerre_defl(roots, polish);
                  else 
                    odd_degree(roots);
                  //quintic(roots);     
#endif 
                }  
            }
        }
    }
#if 0
  void forcehqr()
    {
      force_hqr=1;
    }
#endif
  // get machine precision for "ntype" type (ntype can float, double, long double)
  ntype epsilon()
    {
      return numeric_limits<ntype>::epsilon(); 
    }
  ntype getmax()
    {
      return numeric_limits<ntype>::max();
    }
   void init_const(void)
    {
      meps = epsilon();
      eps05 = pow(numeric_limits<ntype>::epsilon(),0.5);
      maxf= getmax();
      maxdigits = numeric_limits<ntype>::digits10-1;
      //force_hqr=0;
      maxf2 = pow(maxf,0.5)/10.0;
      maxf3 = pow(maxf,1.0/3.0)/10.0;
      scalfact = pow(maxf,1.0/4.0)/1.618034;
      cubic_rescal_fact = pow(maxf, 1.0/3.0)/1.618034;
      goaleps=numeric_limits<ntype>::epsilon();   
      deflated=false;
#ifdef BINI_INIGUESS
      iniguess = new cmplx[n];
#endif
    }
   void set_deflated(void)
     {
       deflated=true;
     }
#ifdef BINI_INIGUESS
   void set_iniguess(cmplx *ig)
     {
       iniguess=ig;
     }
#endif

  void set_output_prec(ntype e)
    {
      goaleps=e;
    } 
  ntype get_output_prec(void)
    {
      return goaleps;
    } 
  rpoly(): rpoly_base<type,N>()
    {
      init_const();
     //std::cout << "max2= " << maxf2<< " max3=" << maxf3 << "\n";
      //printf("%.15G\n", maxf);
      //std::cout << "macheps= " << std::setprecision(35) << meps << "\n";
    }

  rpoly(int nc): rpolybase<ntype,N>(nc)
    {
      init_const(); 
    }
  //rpoly() = default;
#if 0
    {
      //coeff=v;
      //rpoly();
      //cout << "macheps= " << setprecision(35) << eps << "\n";
    }
#endif
#ifdef BINI_INIGUESS
  ~rpoly()
    {
      if (!deflated)
        delete [] iniguess;
    }
#endif
};
// quadratic equation
template<class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N,force_hqr,cmplx>::solve_quadratic(pvector<cmplx,N>&sol)
{
  cmplx r[2];
  ntype a,b;
  a = coeff[1]/coeff[2];
  b = coeff[0]/coeff[2];
  oqs_solve_quadratic(a, b, r);
  sol[0] = r[0];
  sol[1] = r[1];
}
// cubic polynomial
template <class ntype, int N, bool force_hqr,class cmplx> void rpoly<ntype,N,force_hqr, cmplx>::solve_cubic_analytic(pvector<cmplx,N>& sol)
{
  /* solve the cubic coeff[3]*x^3 + coeff[2]*x^2 +  coeff[1]*x + coeff[0] = 0
   * according to the method described in Numerical Recipe book */  
  ntype a, b, c, Q, R, theta, Q3, R2, A, B;
  const ntype sqrt32=sqrt((ntype)3.0)/2.0;
  a = coeff[2]/coeff[3];
  b = coeff[1]/coeff[3];
  c = coeff[0]/coeff[3];
  Q = (Sqr(a) - 3.0*b)/9.0;
  R = (2.0*Sqr(a)*a - 9.0*a*b + 27.0*c)/54.0;
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sol[0] = -2.0*sqrt(Q)*cos(theta/3.0)- a/3.0;
      sol[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a/3.0;
      sol[2] = -2.0*sqrt(Q)*cos((theta-2.0*M_PI)/3.0) - a/3.0;
    }
  else
    {
      A = -copysign((ntype)1.0,R)*pow(abs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = (A+B) - a/3.0;
      sol[1] = cmplx(-0.5*(A+B)-a/3.0,sqrt32*(A-B));
      sol[2] = conj(sol[1]);
      //sol[1] = -0.5*(A+B)-a/3.0+cmplx(0,1)*sqrt32*(A-B);
      //sol[2] = -0.5*(A+B)-a/3.0-cmplx(0,1)*sqrt32*(A-B);
    }
}

// quartics with OQS
template <class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N,force_hqr,cmplx>::oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol)
{
 /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
  * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const ntype PI2=M_PI/2.0, TWOPI=2.0*M_PI;
#ifdef FAST_MATH
  ntype rq3;
#endif
 
  Q = -b/3.0;
  R = 0.5*c;
  if (R==0)
    {
      if (b <= 0)
	{
	  sol=sqrt(-b);
	}
      else
	{
	  sol=0;
	}
      return;
    }
  
  if (abs(Q) < abs(R))
    {
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      RQ = R/Q;
      KK = copysign((ntype)1.0,Q)*(RQ*RQ/Q-1.0);
    }

  if (KK < 0.0)
    {
      sqrtQ=sqrt(Q);
#ifdef FAST_MATH
      // se si use -Ofast rq3 può essere >1 o < -1 di quantità minori di machine epsilon 
      // causando dei NaN
      rq3 = (R/abs(Q))/sqrtQ;
      if (rq3 > 1.0)
        theta = 1.0;
      else if (rq3 < -1.0)
        theta  = M_PI;
      else
        theta = acos(rq3);
#else
      theta = acos((R/abs(Q))/sqrtQ);
#endif
      if (theta < PI2) 
	sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
	sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (abs(Q) < abs(R))
	A = -copysign((ntype)1.0,R)*cbrt(abs(R)*(1.0+sqrt(KK)));
      else
	{
	  A = -copysign((ntype)1.0,R)*cbrt(abs(R)+sqrt(abs(Q))*abs(Q)*sqrt(KK));
	}
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol = A+B;
    }
}
template <class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N,force_hqr,cmplx>::oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, Q3, R2, A, B, sqrtQ;
#ifdef FAST_MATH
  ntype rq3;
#endif
  Q = -b/3.0;
  R = 0.5*c;
  // these number could be made larger for long double */
  //if (abs(Q) > 1E102 || abs(R) > 1E154)
  if (abs(Q) > maxf3 || abs(R) > maxf2)
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
#ifdef FAST_MATH
      // se si use -Ofast rq3 può essere >1 o < -1 di quantità minori di machine epsilon 
      // causando dei NaN
      rq3 = R/sqrt(Q3);
      if (rq3 > 1.0)
        theta = 1.0;
      else if (rq3 < -1.0)
        theta  = M_PI;
      else
        theta = acos(rq3);
#else
      theta = acos(R/sqrt(Q3));
#endif
      sqrtQ=-2.0*sqrt(Q);
      if (theta < M_PI/2) 
	sol = sqrtQ*cos(theta/3.0);
      else 
	sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign((ntype)1.0,R)*pow(abs(R) + sqrt(R2 - Q3),1.0/3.0);
       if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol = A+B; /* this is always largest root even if A=B */
    }
}
template <class ntype, int N, bool force_hqr, class cmplx> void  rpoly<ntype,N,force_hqr, cmplx>::oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (64) (see also the discussion in sec. 2.2 of the manuscript) */
  ntype rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  ntype maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  ntype ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
  int iter;

  diskr=9*a*a-24*b;                    
  /* eq. (67) */
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
	s=-2*b/(3*a+diskr);                     
      else
	s=-2*b/(3*a-diskr);                      
    }
  else
    {      
      s=-a/4;                                    
    }
  /* eqs. (63) */
  aq=a+4*s;                                      
  bq=b+3*s*(a+2*s);                              
  cq=c+s*(2*b+s*(3*a+4*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9;
  hh=aq*cq;     
  
  g=hh-4*dq-3*gg;                       /* eq. (60) */  
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq; /* eq. (61) */          
  oqs_solve_cubic_analytic_depressed(g, h, rmax);
  if (isnan(rmax) || isinf(rmax))
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, rmax);
      if ((isnan(rmax) || isinf(rmax)) && scaled)
	{
	  // try harder: rescale also the depressed cubic if quartic has been already rescaled
	  rfact = cubic_rescal_fact; 
	  rfactsq = rfact*rfact;
	  ggss = gg/rfactsq;
	  hhss = hh/rfactsq;
	  dqss = dq/rfactsq;
	  aqs = aq/rfact;
	  bqs = bq/rfact;
	  cqs = cq/rfact;
	  ggss=bqs*bqs/9.0;
	  hhss=aqs*cqs;   
	  g=hhss-4.0*dqss-3.0*ggss;                       
	  h=(8.0*dqss+hhss-2.0*ggss)*bqs/3-cqs*(cqs/rfact)-(dq/rfact)*aqs*aqs; 
	  oqs_solve_cubic_analytic_depressed(g, h, rmax);
	  if (isnan(rmax) || isinf(rmax))
	    {
	      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, rmax);
	    }
	  rmax *= rfact;
	}
    }
  /* Newton-Raphson used to refine phi0 (see end of sec. 2.2 in the manuscript) */
  x = rmax;
  xsq=x*x;
  xxx=x*xsq;
  gx=g*x;
  f = x*(xsq + g) + h;
  if (abs(xxx) > abs(gx))
    maxtt = abs(xxx);
  else
    maxtt = abs(gx);
  if (abs(h) > maxtt)
    maxtt = abs(h);

  if (abs(f) > maxtt)
    {
      for (iter=0; iter < 8; iter++)
	{   
	  df =  3.0*xsq + g;
	  if (df==0)
	    {
	      break;
	    }
	  xold = x;
	  x += -f/df;
	  fold = f;
	  xsq = x*x;
	  f = x*(xsq + g) + h;
	  if (f==0)
	    {
	      break;
	    } 

	  if (abs(f) >= abs(fold))
	    {
	      x = xold;
	      break;
	    }
    	}
    }
  phi0 = x;
}
template <class ntype, int N, bool force_hqr, class cmplx> ntype  rpoly<ntype,N,force_hqr,cmplx>::oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3)
{
  /* Eq. (21) in the manuscript */
  ntype sum;
  if (b==0.0)
    sum =  abs(d2 + l1*l1 + 2.0*l3);
  else 
    sum =  abs(((d2 + l1*l1 + 2.0*l3)-b)/b);
  if (c==0.0)
    sum += abs(2.0*d2*l2 + 2.0*l1*l3);
  else 
    sum += abs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  if (d==0.0)
    sum += abs(d2*l2*l2 + l3*l3);
  else 
    sum += abs(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
}
template <class ntype, int N, bool force_hqr, class cmplx> 
ntype rpoly<ntype,N, force_hqr,cmplx>::oqs_calc_err_abcd_cmplx(ntype a, ntype b,  ntype c, ntype d, cmplx aq, 
                                                         cmplx bq, cmplx cq, cmplx dq)
{
  /* Eq. (53) in the manuscript for complex alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq) */
  ntype sum;
  sum = (d==0)?abs(bq*dq):abs((bq*dq-d)/d);
  sum += (c==0)?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?abs(aq + cq):abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, bool force_hqr, class cmplx> ntype rpoly<ntype,N,force_hqr,cmplx>::oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq)
{
  /* Eq. (53) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  ntype sum;

  if (d==0.0)
    sum = abs(bq*dq);
  else
    sum = abs((bq*dq-d)/d);

  if (c==0.0)
    sum += abs(bq*cq + aq*dq);
  else 
    sum += abs(((bq*cq + aq*dq) - c)/c);

  if (b==0.0)
    sum +=abs(bq + aq*cq + dq);
  else 
    sum +=abs(((bq + aq*cq + dq) - b)/b);
  if (a==0.0)
    sum +=abs(aq + cq);
  else 
    sum +=abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, bool force_hqr, class cmplx> ntype  rpoly<ntype,N,force_hqr,cmplx>::oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq)
{
  /* Eq. (40) in the manuscript */
  ntype sum;
  if (c==0.0)
    sum = abs(bq*cq + aq*dq);
  else
    sum = abs(((bq*cq + aq*dq) - c)/c);
  if (b==0.0)
  sum +=abs(bq + aq*cq + dq);
  else 
  sum +=abs(((bq + aq*cq + dq) - b)/b);
  if (a==0.0)
    sum +=abs(aq + cq);
  else 
    sum +=abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N,force_hqr, cmplx>::oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  ntype x02, errf, errfold, xold[4], x[4], dx[4], det, Jinv[4][4], fvec[4], vr[4];
  x[0] = AQ;
  x[1] = BQ;
  x[2] = CQ;
  x[3] = DQ;
  vr[0] = d;
  vr[1] = c;
  vr[2] = b;
  vr[3] = a;
  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  for (k1=0; k1 < 4; k1++)
    {
      if (vr[k1]==0)
        errf += abs(fvec[k1]);
      else
        errf +=abs(fvec[k1]/vr[k1]);
    }
  for (iter = 0; iter < 8; iter++)
    {
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - 2.0*x[3]) + x[3]*(x[0]*x02 + x[3]);
      if (det==0.0)
	break;
      Jinv[0][0] = x02;
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = -x[1]*Jinv[0][1] - x[0]*Jinv[0][2]; 
      Jinv[1][0] = x[0]*Jinv[0][0] + Jinv[0][1];
      Jinv[1][1] = -x[1]*Jinv[0][0];
      Jinv[1][2] = -x[1]*Jinv[0][1];   
      Jinv[1][3] = -x[1]*Jinv[0][2];
      Jinv[2][0] = -Jinv[0][0];
      Jinv[2][1] = -Jinv[0][1];
      Jinv[2][2] = -Jinv[0][2];
      Jinv[2][3] = Jinv[0][2]*x[2] + Jinv[0][1]*x[3];
      Jinv[3][0] = -x[2]*Jinv[0][0] - Jinv[0][1];
      Jinv[3][1] = Jinv[0][0]*x[3];
      Jinv[3][2] = x[3]*Jinv[0][1];
      Jinv[3][3] = x[3]*Jinv[0][2];
      for (k1=0; k1 < 4; k1++)
	{
	  dx[k1] = 0;
	  for (k2=0; k2 < 4; k2++)
	    dx[k1] += Jinv[k1][k2]*fvec[k2];
	}
      for (k1=0; k1 < 4; k1++)
      	xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1]/det;
	}
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 
      errfold = errf;
      errf=0;
      for (k1=0; k1 < 4; k1++)
	{
        if (vr[k1]==0)
	  errf += abs(fvec[k1]);
        else
	  errf += abs(fvec[k1]/vr[k1]);
	}
      if (errf==0)
	break;
      if (errf >= errfold)
	{
	  for (k1=0; k1 < 4; k1++)
	    x[k1] = xold[k1];
	  break;
	}
    }
  AQ=x[0];
  BQ=x[1];
  CQ=x[2];
  DQ=x[3];
}
template <class ntype, int N, bool force_hqr, class cmplx> void  rpoly<ntype,N, force_hqr,cmplx>::oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2])
{ 
  ntype div,sqrtd,diskr,zmax,zmin;
  diskr=a*a-4*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
	zmin=0.0;
      else
	zmin=b/zmax;
      roots[0]=cmplx(zmax,0.0);
      roots[1]=cmplx(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=cmplx(-a/2,sqrtd/2);
      roots[1]=cmplx(-a/2,-sqrtd/2);      
    }   
}

template <class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N,force_hqr,cmplx>::oqs_quartic_solver(pvector<cmplx,N>& roots)
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * 
   * if coeff[4] != 0 
   *
   * the four roots will be stored in the complex array roots roots[] 
   *
   * */
  cmplx acx1, bcx1, ccx1, dcx1,acx,bcx,ccx,dcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
  ntype l2m[12], d2m[12], res[12], resmin, bl311, dml3l3, err0=0, err1=0, aq1, bq1, cq1, dq1; 
  ntype a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3],gamma,del2;
  int realcase[2], whichcase, k1, k, kmin, nsol;
  ntype rfactsq, rfact=1.0;

  if (coeff[4]==0.0)
    {
      cout << "That's not a quartic!\n";
      return;
    }

  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  oqs_calc_phi0(a,b,c,d,phi0,0);
  //cout << "phi0=" << phi0 << "\n";
  // simple polynomial rescaling
  if (isnan(phi0)||isinf(phi0))
    {
      rfact = scalfact;
      a /= rfact;
      rfactsq = rfact*rfact;
      b /= rfactsq;
      c /= rfactsq*rfact;
      d /= rfactsq*rfactsq;
      oqs_calc_phi0(a,b,c,d,phi0,1);
    }

  l1=a/2;          /* eq. (4) */                                        
  l3=b/6+phi0/2;   /* eq. (6) */                                
  del2=c-a*l3;     /* defined just after eq. (20) */                             
  nsol=0;
  bl311 =2.*b/3.-phi0-l1*l1;   /* This is d2 as defined in eq. (18)*/ 
  dml3l3 = d-l3*l3;            /* dml3l3 is d3 as defined in eq. (9) with d2=0 */ 
  
  /* Three possible solutions for d2 and l2 (see eqs. (18)-(20) and discussion which follows) */
  if (bl311!=0.0)
    {
      d2m[nsol] = bl311;  
      l2m[nsol] = del2/(2.0*d2m[nsol]);   
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }
  if (del2!=0)
    {
      l2m[nsol]=2*dml3l3/del2;
      if (l2m[nsol]!=0)
	{
  	  d2m[nsol]=del2/(2*l2m[nsol]);
	  res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
	  nsol++;
	}

      d2m[nsol] = bl311;
      l2m[nsol] = 2.0*dml3l3/del2;
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }

  if (nsol==0)
    {
      l2=d2=0.0;
    }
  else
    {
      /* we select the (d2,l2) pair which minimizes errors */
      for (k1=0; k1 < nsol; k1++)
	{
	  if (k1==0 || res[k1] < resmin)
	    {
	      resmin = res[k1];
	      kmin = k1;	
	    }
	}
      d2 = d2m[kmin];
      l2 = l2m[kmin];
    }
  whichcase = 0; 
  if (d2 < 0.0) 
    {
      /* Case I eqs. (27)-(30) */
      gamma=sqrt(-d2);                               
      aq=l1+gamma;                                  
      bq=l3+gamma*l2;                              

      cq=l1-gamma;                                
      dq=l3-gamma*l2;                            
      if(abs(dq) < abs(bq))
	dq=d/bq;                                
      else if(abs(dq) > abs(bq))
	bq=d/dq;                               
      if (abs(aq) < abs(cq))
	{
	  nsol=0;
	  if (dq !=0)
	    {
	      aqv[nsol] = (c - bq*cq)/dq;    /* eq. (37) */
	      errv[nsol]=oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	      nsol++;
	    }
	  if (cq != 0) 
	    {
	      aqv[nsol] = (b - dq - bq)/cq;  /* eq. (38) */
	      errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	      nsol++;
	    }
	  aqv[nsol] = a - cq;                /* eq. (39) */
	  errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	  nsol++;
	  /* we select the value of aq (i.e. alpha1 in the manuscript) which minimizes errors */
	  for (k=0; k < nsol; k++)
	    {
	      if (k==0 || errv[k] < errmin)
		{
		  kmin = k;
		  errmin = errv[k];
		}
	    }
	  aq = aqv[kmin];
	}
      else 
	{
	  nsol = 0;
	  if (bq != 0)
	    { 
	      cqv[nsol] = (c - aq*dq)/bq;              /* eq. (44) */
	      errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	      nsol++;
	    }
	  if (aq != 0)
	    {
	      cqv[nsol] = (b - bq - dq)/aq;            /* eq. (45) */
	      errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	      nsol++;
	    }
	  cqv[nsol] = a - aq;                          /*  eq. (46) */
	  errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	  nsol++;	  
	  /* we select the value of cq (i.e. alpha2 in the manuscript) which minimizes errors */
	  for (k=0; k < nsol; k++)
	    {
	      if (k==0 || errv[k] < errmin)
		{
		  kmin = k;
		  errmin = errv[k];
		}
	    }
	  cq = cqv[kmin];
	}
      realcase[0]=1;
    }
   else if (d2 > 0)   
    {
      /* Case II eqs. (47)-(50) */
      gamma=sqrt(d2); 
      acx=cmplx(l1,gamma);  
      bcx=cmplx(l3,gamma*l2);
      ccx = conj(acx);
      dcx = conj(bcx);
      realcase[0] = 0; 
    }
  else 
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is better) */
  if (realcase[0]==-1 || (abs(d2) <= meps*oqs_max3(abs(2.*b/3.), abs(phi0), l1*l1))) 
    {
      d3 = d - l3*l3;
      if (realcase[0]==1)
	err0 = oqs_calc_err_abcd(a, b, c, d, aq, bq, cq, dq);
      else if (realcase[0]==0)
	err0 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);
      if (d3 <= 0)
	{
	  realcase[1] = 1;
	  aq1 = l1;   
	  bq1 = l3 + sqrt(-d3);
	  cq1 = l1;
	  dq1 = l3 - sqrt(-d3);
	  if(abs(dq1) < abs(bq1))  
	    dq1=d/bq1;                                        
	  else if(abs(dq1) > abs(bq1))
	    bq1=d/dq1;                                       
	  err1 = oqs_calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (53) */
	}
      else /* complex */
	{
	  realcase[1] = 0;
	  acx1 = l1;
	  bcx1 = cmplx(l3,sqrt(d3));
	  ccx1 = l1;
	  dcx1 = conj(bcx1);
	  err1 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1); 
	}
      if (realcase[0]==-1 || err1 < err0)
	{
          whichcase=1; // d2 = 0
	  if (realcase[1]==1)
	    {
	      aq = aq1;
	      bq = bq1;
	      cq = cq1;
	      dq = dq1;
	    }
	  else
	    {
	      acx = acx1;
	      bcx = bcx1;
	      ccx = ccx1;
	      dcx = dcx1;
	    }
	}
    }
  if (realcase[whichcase]==1)
    {
      /* if alpha1, beta1, alpha2 and beta2 are real first refine 
       * the coefficient through a Newton-Raphson */
      oqs_NRabcd(a,b,c,d,aq,bq,cq,dq);      
      /* finally calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
      oqs_solve_quadratic(aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      oqs_solve_quadratic(cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];
    }
  else
    {
      /* complex coefficients of p1 and p2 */
      if (whichcase==0) // d2!=0
	{
	  cdiskr=acx*acx/((cmplx)4.0)-bcx;               
	  /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
	  zx1=-acx/((cmplx)2.0)+sqrt(cdiskr);
	  zx2=-acx/((cmplx)2.0)-sqrt(cdiskr);
	  if(abs(zx1) > abs(zx2))
	    zxmax=zx1;
	  else
	    zxmax=zx2;
	  zxmin=bcx/zxmax;        
	  roots[0]=zxmin;
	  roots[1]=conj(zxmin);
	  roots[2]=zxmax;
	  roots[3]=conj(zxmax);
	}
      else // d2 ~ 0
	{
	  /* never gets here! */
	  cdiskr=sqrt(acx*acx-cmplx(4.0)*bcx);
	  zx1 = -cmplx(0.5)*(acx+cdiskr);
	  zx2 = -cmplx(0.5)*(acx-cdiskr);
	  if (abs(zx1) > abs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = bcx/zxmax;
	  roots[0] = zxmax;
	  roots[1] = zxmin;
	  cdiskr=sqrt(ccx*ccx-cmplx(4.0)*dcx);
	  zx1 = -cmplx(0.5)*(ccx+cdiskr);
	  zx2 = -cmplx(0.5)*(ccx-cdiskr);
	  if (abs(zx1) > abs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = dcx/zxmax;
	  roots[2]= zxmax;
	  roots[3]= zxmin;
	}
    }
  if (rfact!=1.0)
    {
      for (k=0; k < 4; k++)
        {
          roots[k] *= rfact;
        }
    }
}
// NR method based on HQR decomposition
template<class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N,force_hqr,cmplx>::balance(pmatrixq<ntype,N>& a)
{
  const ntype RADIX=numeric_limits<ntype>::radix;
  int i, j;
  ntype *scale;
  int done=0;
  ntype r, c, g, f, s, sqrdx=RADIX*RADIX;
  scale = new ntype[n];
  for (i=0; i < n; i++)
    scale[i] = 1.0;
  while (!done) 
    {
      done=1;
      for (i=0;i<n;i++) 
	{
	  //Calculate row and column norms.
	  //If both are nonzero,
	  //find the integer power of the machine radix that comes closest to balancing the matrix.
	  r=0.0;
	  c=0.0;
	  for (j=0;j<n;j++)
	    if (j != i) 
	      {
		c += abs(a[j][i]);
		r += abs(a[i][j]);
	      }
	  if (c != 0.0 && r != 0.0) 
	    {
	      g=r/RADIX;
	      f=1.0;
	      s=c+r;
	      while (c<g) {
		f *= RADIX;
		c *= sqrdx;
	      }
	      g=r*RADIX;
	      while (c>g) 
		{
		  f /= RADIX;
		  c /= sqrdx; 
		}
	      if ((c+r)/f < 0.95*s) 
		{
		  done=0;
		  g=1.0/f;
		  scale[i] *= f;
		  for (j=0;j<n;j++) a[i][j] *= g; //Apply similarity transformation
		  for (j=0;j<n;j++) a[j][i] *= f;
		}
	    }
	}
    }
  delete [] scale;
}

template<class ntype, int N, bool force_hqr,class cmplx> void rpoly<ntype,N, force_hqr,cmplx>::hqr(pmatrixq<ntype,N>& a, pvector<cmplx,N>& wri, int& ok)
{
  int nn,m,l,k,j,its,i,mmin;
  ntype z,y,x,w,v,u,t,s,r=0.0,q=0.0,p=0.0, anorm=0.0;
  //const double EPS=2.2204460492503131E-16;
  for (i=0;i<n;i++)
    //Compute matrix no rm for possible use in lo- cating single small sub diagonal element.
    for (j=IMAX(i-1,0);j<n;j++)
      anorm += abs(a[i][j]);
  nn=n-1;
  t=0.0;
  ok = 1;
  //Gets changed only by an exceptional shift.
  while (nn >= 0) 
    {
      //Begin search for next eigenvalue.
      its=0;
      do 
	{
	  for (l=nn;l>0;l--)
	    {
	      //Begin iteration: look for single small sub di- agonal element.
	      s=abs(a[l-1][l-1])+abs(a[l][l]);
	      if (s == 0.0)
		s=anorm;
	      if (abs(a[l][l-1]) + s == s)
	       	{
	  	  a[l][l-1]=0.0;
	  	  break; 
	    	}
	    }
	  x=a[nn][nn];
	  if (l == nn)
	    {
	      //One root found.  
	      wri[nn--]=x+t;
	    } 
	  else
	    {
	      y=a[nn-1][nn-1];
	      w=a[nn][nn-1]*a[nn-1][nn];
	      if (l == nn-1)
		{
		  //Two roots found...
		  p=0.5*(y-x);
		  q=p*p+w;
		  z=sqrt(abs(q));
		  x += t;
		  if (q >= 0.0)
		    {
		      //...a real pair.
		      z=p+SIGN<ntype>(z,p);
		      wri[nn-1]=wri[nn]=x+z;
		      if (z != 0.0)
			wri[nn]=x-w/z;
		    } 
		  else
		    {
		      //...a complex pair.
		      wri[nn]=cmplx(x+p,-z);
		      wri[nn-1]=conj(wri[nn]);
		    }
		  nn -= 2;
		} 
	      else
		{
		  //No roots found.  Continue iteration.
		  if (its == 480)
		    {
		      cout << "Too many iterations in hqr";
		      ok=0;
		      return;
		      //exit(-1);
		    }
		  if (its % 10 == 0 && its > 0)
		    {
		      //Form exceptional shift.
		      t += x;
		      for (i=0;i<nn+1;i++)
			a[i][i] -= x;
		      s=abs(a[nn][nn-1])+abs(a[nn-1][nn-2]);
		      y=x=0.75*s;
		      w = -0.4375*s*s;
		    }
		  ++its;
		  for (m=nn-2;m>=l;m--)
		    {
		      //Form shift and then look for 2 consecutive small sub- diagonal elements.
		      z=a[m][m];
		      r=x-z;
		      s=y-z;
		      p=(r*s-w)/a[m+1][m]+a[m][m+1];
		      //Equation (W ebnote 16.21).
		      q=a[m+1][m+1]-z-r-s;
		      r=a[m+2][m+1];
		      s=abs(p)+abs(q)+abs(r);
		      //Scale to prevent over flow or under flow.
		      p /= s;
		      q /= s;
		      r /= s;
		      if (m == l) 
			break;
		      u=abs(a[m][m-1])*(abs(q)+ abs(r));
		      v=abs(p)*(abs(a[m-1][m-1])+abs(z)+abs(a[m+1][m+1]));
		      if (u <= meps*v)
			break;
		      //Equation (W ebnote 16.24).
		    }
		  for (i=m;i<nn-1;i++)
		    {
		      a[i+2][i]=0.0;
		      if (i != m) a[i+2][i-1]=0.0;
		    }
		  for (k=m;k<nn;k++)
		    {
		      //Double QR step on rows l to nn and columns m to nn .
		      if (k != m) 
			{
			  p=a[k][k-1];
			  //Begin setup of Householder vector.
			  q=a[k+1][k-1];
			  r=0.0;
			  if (k+1 != nn) 
			    r=a[k+2][k-1];
			  if ((x=abs(p)+abs(q)+abs(r)) != ((ntype)0.0))
			    {
			      p /= x;
			      //Scale to prevent over flow or under flow.
			      q /= x;
			      r /= x;
			    }
			}
		      if ((s=SIGN<ntype>(sqrt(p*p+q*q+ r*r),p)) != ((ntype)0.0))
			{
			  if (k == m) 
			    {
			      if (l != m)
				a[k][k-1] = -a[k][k-1];
			    } 
			  else
			    a[k][k-1] = -s*x;
			  p += s;
			  //Equations (Webnote 16.22).
			  x=p/s;
			  y=q/s;
			  z=r/s;
			  q /= p;
			  r /= p;
			  for (j=k;j<nn+1;j++)
			    {
			      //Row modification.
			      p=a[k][j]+q*a[k+1][j];
			      if (k+1 != nn)
				{
				  p += r*a[k+2][j];
				  a[k+2][j] -= p*z;
				}
			      a[k+1][j] -= p*y;
			      a[k][j] -= p*x;
                            }
			  mmin = nn < k+3 ? nn : k+3;
			  for (i=l;i<mmin+1;i++)
			    {
			      //Column modification.
			      p=x*a[i][k]+y*a[i][k+1 ];
			      if (k+1 != nn) {
				p += z*a[i][k+2];
				a[i][k+2] -= p*r;
			      }
			      a[i][k+1] -= p*q;
			      a[i][k] -= p;
			    }
			}
		    }
		}
	    }
	} 
      while (l+1 < nn);
    }
}
template<class ntype, int N, bool force_hqr, class cmplx> void rpoly<ntype,N, force_hqr,cmplx>::QRfactorization(pmatrixq<ntype,N>& hess, pvector<cmplx,N>& sol, int& ok)
{
  /* pag. 615 Num. Rec. */  
  balance(hess);
  hqr(hess, sol, ok);
}
template<class ntype, int N, bool force_hqr,class cmplx> void rpoly<ntype,N, force_hqr,cmplx>::solve_numrec(pvector<cmplx,N>& csol, int& ok)
{
  /* Find all the roots of a polynomial with real coefficients, 
   * coeff[4]*x^4+coeff[3]*x^3+coeff[2]*x^2+coeff[1]*x+coeff[0], 
   * The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pag. 497 Num. Rec. */
  int m=n;
  pmatrixq<ntype,N> hess;
  int j, k;

  if constexpr (N < 0)
    {
      hess.allocate(n);
    }
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorization(hess, csol, ok);
}
