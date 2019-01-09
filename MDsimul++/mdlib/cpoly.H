#include "pvector.H"
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
//#define LAGUERRE_POLISH // migliora drasticamente con maxiter=8 ma per polinomi di alto grado (500) crea problemi
#define USE_BRENT
#define USE_LAGUERRE_ONLY 
//#define POLISH_NR_REAL 
//#define USE_LAGUER_REAL
//#define FAST_MATH
using namespace std;
template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

template <class cmplx, int N> 
class cpoly_base_static 
{
public:
  int n;
  constexpr static int dynamic = false;
  pvector<cmplx, N+1> coeff;
  pvector<cmplx, N+1> cmon;
 // pvector<ntype, N> cmonder;
 // pvector<ntype, N-1> cmondder;
  pvector<cmplx, N-1> bd;
  pvector<cmplx, N-1> bd0;
  pvector<cmplx, N-1> bd1;
  pvector<cmplx, N> deflcoeff0;
  pvector<cmplx, N> deflcoeff1;
  pvector<cmplx, N> deflcoeff;
  void set_coeff(pvector<cmplx,N+1> v)
    {
      coeff = v;
      //coeff.show("coefficienti=");
    }
  cpoly_base_static()
    {
      n=N;
    }
  cpoly_base_static(int nc): cpoly_base_static()
  {
    n=nc;
    // empty body
  }
};
template <class cmplx, int N> 
class cpoly_base_dynamic 
{
public:
  int n;
  constexpr static int dynamic = true;
  pvector<cmplx> coeff;
  pvector<cmplx> cmon;
  //pvector<ntype> cmonder;
  //pvector<ntype> cmondder;
  pvector<cmplx> bd0;
  pvector<cmplx> bd1;
  pvector<cmplx> bd;
  pvector<cmplx> deflcoeff0;
  pvector<cmplx> deflcoeff1;
  pvector<cmplx> deflcoeff;
  cpoly_base_dynamic() = default;
  void set_coeff(pvector<cmplx,-1> v)
    {
      coeff = v;
#if 0
      coeff.show();
      cout << "n=" << n << "\n";
      cout << "v=" << v[n]<< "\n";
        v.show();
#endif
    }
  void use_vec(int nc, cmplx* coeffv, cmplx* cmonv,
               cmplx* bd0v, cmplx *bd1v, cmplx* bdv, cmplx *deflcoeff0v, cmplx* deflcoeff1v, cmplx* deflcoeffv)
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
 
  cpoly_base_dynamic(int nc): coeff(nc+1), cmon(nc+1), bd0(nc-1), bd1(nc-1), bd(nc-1),
  deflcoeff0(nc), deflcoeff1(nc), deflcoeff(nc)
    {
      n=nc;
    }
  ~cpoly_base_dynamic() = default;
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

template <class cmplx, int N> using cpolybase = 
typename std::conditional<(N>0), cpoly_base_static <cmplx, N>,
	 cpoly_base_dynamic <cmplx, N>>::type;
 
template <class cmplx, int N=-1, class ntype=double> 
class cpoly: public numeric_limits<ntype>, public cpolybase<cmplx,N> {
  using cpolybase<cmplx,N>::n;
  using cpolybase<cmplx,N>::coeff;
  using cpolybase<cmplx,N>::cmon;
  //using cpolybase<ntype,N>::cmonder;
  //using cpolybase<ntype,N>::cmondder;
  using cpolybase<cmplx,N>::bd0;
  using cpolybase<cmplx,N>::bd1;
  using cpolybase<cmplx,N>::bd;
  using cpolybase<cmplx,N>::deflcoeff0;
  using cpolybase<cmplx,N>::deflcoeff1;
  using cpolybase<cmplx,N>::deflcoeff;
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
  template <class vtype>
  using pvecNm1 = typename std::conditional<(N>0),  pvector<vtype, N-1>,
	 pvector<vtype, -1>>::type;
  template <class vtype>
  using pvecNp1 = typename std::conditional<(N>0),  pvector<vtype, N+1>,
	 pvector<vtype, -1>>::type;
  template <class vtype>
  using pvecNm2 = typename std::conditional<(N>0),  pvector<vtype, N-2>,
	 pvector<vtype, -1>>::type;
  using cpolyNm1 = typename std::conditional<(N>0),  cpoly<cmplx, N-1,ntype>,
	 cpoly<cmplx, -1,ntype>>::type;
  using cpolyNm2 = typename std::conditional<(N>0),  cpoly<cmplx, N-2,ntype>,
	 cpoly<cmplx, -1,ntype>>::type;
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

  inline void solve_quadratic(pvector<cmplx,N>& sol);
#if 0 
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
#endif
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
      for (i=N; i >= 0; i--)
	{
	  if (abs(coeff[i]) > 0.0)
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

  cmplx evalddpoly(cmplx x)
    {
      // evaluate second derivative of polynomail via Horner's formula 
      cmplx bn=0.0;
      if (n == 1)
        return 0;
      for (int i=n-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmon[i+2] + bn*x;
        }
      return bn;
    }
  ntype calc_err_q(pvector<cmplx,N> c, cmplx r0)
    {
      int i;
      ntype sum=0;
      for (i=0; i < n; i++)
        {
          if (i == 0){
            if (cmon[i]==cmplx(0.0,0.0))
              {
                sum+=abs(r0*c[i]+cmon[i]);} 
            else 
              {
                sum=+abs((r0*c[i]+cmon[i])/cmon[i]);}
          }
          else 
            if (cmon[i]==cmplx(0.0,0.0))
              sum +=abs(-c[i-1] + r0*c[i] + cmon[i]);
            else
              sum +=abs((-c[i-1] + r0*c[i]+cmon[i])/cmon[i]);
        }
      //sum += (cmon[0]==0.0)?abs(r0*c[0]+cmon[0]):abs((r0*c[0]+cmon[0])/cmon[0]);
      return sum;
    }
#if 0
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
bool laguer_real(cmplx &xc)
  {
    const int MR=8,MT=10,MAXIT=MT*MR;
    const ntype EPS=numeric_limits<ntype>::epsilon();
    static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    int iter,j;
    ntype abx,abp,abm,err, phi, rho,s, pr;
    //cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    ntype dx[2],x1[2],b[2],d[2],f[2],g[2],h[2],sq[2],gp[2],gm[2],g2[2], fb[2], pol[2];
    ntype x[2], f0, d0, b0, invn, g0, invgp[2], sarg[2];
    cmplx bc,dc,fc,sqrtcmplx, fovb, dovb;
    int m=n;
    //cout << "r0 ini = "<<xc << "\n"; 
    x[0] = xc.real();
    x[1] = xc.imag();
    for (iter=1;iter<=MAXIT;iter++) 
      {
        //its=iter;
        b[0]=real(cmon[m]);
        b[1]=imag(cmon[m]);
        err=sqrt(b[0]*b[0]+b[1]*b[1]);
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
          b[0] = x[0]*b[0] - x[1]*b[1] + real(cmon[j]);
          b[1] = x[0]*b[1] + x[1]*b0 + imag(cmon[j]);
          err=sqrt(b[0]*b[0]+b[1]*b[1])+abx*err;
        }
        err *= EPS;
        if (isinf(b[0]) || isnan(b[0]) || 
            isinf(b[1]) || isnan(b[1]))
          {
            //cout << "f= " << f << "d= " << d << "b=" << b << "\n";
            //cout << "x=" << x << "\n";
            //cout << "p(x)=" << evalpoly(x) << "\n";
            return false;
          }

#if 1
        s=0;
        pr=1.0;
        for (j=0; j <= n; j ++) 
          {
            s+=coeff[j].real()*pr;
            pr *= abx;
          }
        if (abs(cmplx(b[0],b[1])) <= 2.0*EPS*(4.0*ntype(n)+1)*s) // stopping criterion of bini 
          {
            xc = cmplx(x[0],x[1]);
            return true;
          }
#else
        if (abs(cmplx(b[0],b[1])) <= err)
          {
            xc = cmplx(x[0],x[1]);
            return true;
          }
#endif          
        bc=cmplx(b[0],b[1]);
        dc=cmplx(d[0],d[1]);
        fc=cmplx(f[0],f[1]);
        dovb=dc/bc;
        fovb=fc/bc;
        g[0]=dovb.real();
        g[1]=dovb.imag();
        g2[0] = g[0]*g[0] - g[1]*g[1];
        g2[1] = 2.0*g[0]*g[1];
        fb[0] = fovb.real();
        fb[1] = fovb.imag();
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
    return false;
  }
bool laguer_cmplx(cmplx &x)
  {
    const int MR=8,MT=10,MAXIT=MT*MR;// MT changed from 10 to 30 
    const ntype EPS=numeric_limits<ntype>::epsilon();
    static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    int iter,j;
    ntype abx,abp,abm,err, M, s, pr;
    cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    int m;
    m=n;
    //cout << "r0 ini = "<< x << "\n"; 
    for (iter=1;iter<=MAXIT;iter++) 
      {
        //its=iter;
        b=cmon[m];
        err=abs(b);
        d=f=0.0;
        abx=abs(x);
        for (j=m-1;j>=0;j--) {
          f=x*f+d;
          d=x*d+b;
          b=x*b+cmon[j];
          err=abs(b)+abx*err;
        }
        err *= EPS;
        //cout << "f= " << f << "d= " << d << "b=" << b << "\n";
        if (isinf(b.real()) || isinf(b.imag()) || isnan(b.real()) ||
            isnan(b.imag()))
          return false;
#if 1
        s=0;
        pr=1.0;
        for (j=0; j <= n; j ++) 
          {
            s+=coeff[j].real()*pr;
            pr *= abx;
          }
        if (abs(b) <= 2.0*EPS*(4.0*ntype(n)+1)*s) // stopping criterion of bini 
          return true;
#else
 
        if (abs(b) <= err) 
          {
            //if (iter > 15)
            //cout << "iter=" << iter << "\n";
            //cout << "x=" << x << " iter=" << iter << " abs(M-1)=" << abs(M-1) << "\n";
            return true;
          }
#endif
        g=d/b;
        g2=g*g;
        h=g2-((ntype)2.0)*f/b;
        // SONO ARRIVATO QUI CON IL CASTING
        sq=sqrt( ((ntype)(m-1))*( ((ntype)m)*h-g2));
        gp=g+sq;
        gm=g-sq;
        abp=abs(gp);
        abm=abs(gm);

        if (abp < abm) gp=gm;
        //cout << "abp=" << abp << " abm=" << abm << "\n";
        //cout << "1/gp=" << 1.0/gp << "\n";
        dx=(MAX(abp,abm) > 0.0) ? ((cmplx)m)/((cmplx)gp) : polar(((ntype)(1+abx)),((ntype)iter));
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
    return false;
  }
#if 0
bool laguer_wrap(cmplx &x)
  {
    if (is_same<ntype,double>::value || is_same<ntype,long double>::value)
      return laguer_real(x);
    else
      return laguer_real(x);
  }
#endif
void refine_root(cmplx &r0)
  {
     int iter,i;
     ntype err;
     // root polishing by NR
#ifdef POLISH_NR_REAL
     // root polishing by NR
     ntype r0old[2], p[2], r0new[2];
     ntype p1[2], p0, p10;
     ntype errold, invnp1;
     cmplx povp1;
     r0new[0]=r0.real();
     r0new[1]=r0.imag();
     //cout << "r0new=" << r0new[0] << " " << r0new[1] << "\n";
     for (iter=0; ;iter++)
       {
         p[0]=real(cmon[n])*r0new[0]-imag(cmon[n])*r0new[1] + real(cmon[n-1]);
         p[1]=real(cmon[n])*r0new[1]+imag(cmon[n])*r0new[0] + imag(cmon[n-1]);    
         p1[0]=real(cmon[n]);
         p1[1]=imag(cmon[n]);
         for(i=n-2;i>=0;i--) {
           p10=p1[0];
           p1[0]=p[0]+p1[0]*r0new[0]-p1[1]*r0new[1];
           p1[1]=p[1]+p10*r0new[1]+p1[1]*r0new[0];
           p0=p[0];
           p[0]=real(cmon[i])+p[0]*r0new[0]-p[1]*r0new[1];
           p[1]=imag(cmon[i])+p0*r0new[1]+p[1]*r0new[0];
         }
         if (iter > 0)
           errold=err;
         err = abs(cmplx(p[0],p[1]));// abs(p[0])+abs(p[1]);

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
//if (iter > 2)
  //   cout << "NR iter=" << iter << " err=" << err << "\n";

#else
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
     cout << "NR iter=" << iter << " err=" << err << "\n";
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
     //if (iter > 2)
     //cout << "NR iter=" << iter << " err=" << err << "\n";
#endif
   }

 void laguerre_defl(pvector<cmplx,N>& roots,bool polish=false, bool backup=false)
   {
     // find initial guess for x here
     cmplx r0;
     const ntype goldfact=1.618034;
     ntype sig;// sigold;
     //bool mulsig=true;
     ntype err0,err1, errA, errB;
     pvector<cmplx,3> ad;
     typename std::conditional<(N < 0),pvector<cmplx>,bool>::type cbak;
     int iter, i;
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
      if (coeff[n]==cmplx(0.0,0.0))
        {
          cout << "That's not an " << n << "degree polynomial!\n";
          return;
        }
      cmon[n]=1.0;
      for (i=n-1; i >=0; i--)
        cmon[i]=coeff[i]/coeff[n];
      //absolute value of coefficients are used in Bini's stopping criterion
      //we use coeff vector to store them since it won't be used anymore from here on
      for (i=0; i <= n; i++) 
        coeff[i] = cmplx(abs(cmon[i]),0); 

      sig=1.0;//pow(std::numeric_limits<ntype>::max(),0.1)/goldfact;
      iter = 0;
      // find first zero using laguerre
      // if it does not converge generate a random initial guess
      // and retry (a limited number of times though!)
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
      //cout << "r0=" << r0 << "\n";
      r0=cmplx(0.0,0.0);
      int MAXSTPS=1000;
      while (!laguer_real(r0))
        {
          if (iter==0)
            r0=cmplx(0.0,1.0);
          else
            r0=polar(ntype(drand48()),ntype(2.0*M_PI*(drand48()-0.5)));
          iter++;
          //if (iter%MT==0)
            //sig*=goldfact;
          if (iter > MAXSTPS)
            {
              cout << "Laguerre can not converge, I am sorry...\n";
              cmon.show("coeff=");
              cout << "size=" << cmon.size()<<"\n";
              exit(-1); 
            }
        };
      //cout << "p(" << r0 << ")=" << evalpoly(r0) << "\n";
#if 0
      if (iter > 0)
        cout << "N=" << n << "iter=" << iter << "sig=" << sig << "\n";
      if (isnan(r0.real())|| isinf(r0.real()))
        {
          cout << "bah\n";
          exit(-1);
        }
#endif
#ifdef LAGUERRE_POLISH
      if (is_same<ntype,double>::value || is_same<ntype,long double>::value)
        {
          refine_root(r0);
        }
#endif
#if 0
      if (isnan(r0.real()||isinf(r0.real())))
        {
          cout << "mah...\n"; 
          exit(1);
        }
#endif
          cpolyNm1 qp;
          // pvector<cmplx,N-1> qr;
          pvecNm1<cmplx> qr;
          //cout << "minus one r0=" << r0 << "\n";
#if 0
          if (abs(imag(r0)) <= 2.0*meps*abs(real(r0)))
            {
              r0=cmplx(real(r0),0.0);
            }
#endif
          //cout << "n=" << n << " r0=" << r0 << "\n";
          roots[0] = r0;
#if 0
          cmplx a,b,c;
          b=cmon[n];
          for (i=n-1;i>=0;i--) {
            c=cmon[i];
            deflcoeff0[i]=b;
            b=r0*b+c;
          }
#else
          for (i=n; i >=1; i--)
            {
              deflcoeff0[i-1] = cmon[i];
              if (i <= n-1)
                deflcoeff0[i-1] += r0*deflcoeff0[i];
            }
#endif
#if 0
          if (r0!=cmplx(0.0,0.0))
            {
              //forward deflation
#if 0
              deflcoeff1[N-1]=1.0;
              for (i=0; i < N-1; i++)
                {
                  deflcoeff1[i] = cmon[i]; 
                  if (i >= 1)
                    deflcoeff1[i] -= deflcoeff1[i-1];
                  deflcoeff1[i] /= -r0.real();
                  //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff1[" << i << "]=" << deflcoeff1[i]<< "\n";
                }
#else
              deflcoeff1[0]=-cmon[0]/r0;
              deflcoeff1[n-1] = 1.0;
              for (i=1; i < n-1; i++)
                {
                  deflcoeff1[i] = (deflcoeff1[i-1] - cmon[i])/r0;
                }
#endif
              // choose the set which yields smaller error 
              err0=calc_err_q(deflcoeff0,r0);
              err1=calc_err_q(deflcoeff1,r0);
              
              if (err0 < err1)
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
              deflcoeff = deflcoeff0;
            }
#else
          deflcoeff = deflcoeff0;
#endif
          //cout << "N=" << n << " r0=" << r0<< "\n";
          //deflcoeff0.show("boh2A");
          //deflcoeff1.show("boh2B");
#ifdef USE_DEFL_NR
          //ntype errP, errD;
          //errP = calc_err_q(deflcoeff, r0.real());
          use_NR_defl_lin(deflcoeff,r0);
          //errD = calc_err_q(deflcoeff, r0.real());
#if 0
          if (errD < errP)
            {
              cout << "errP=" << errP << " errD=" << errD << "\n";

            }
#endif
#endif

          if constexpr (N < 0)
            {
              //qp.allocate(n-1);
              //qr.allocate(n-1);
              qr.use_vec(n-1,&(roots[1])); 
              qp.use_vec(n-1,coeff.v, cmon.v, bd0.v, bd1.v, bd.v, deflcoeff0.v, deflcoeff1.v, deflcoeff.v);
            }


          //refine deflated quartic coefficients thourh Newton-Raphson method
          qp.set_coeff(deflcoeff);
          //qp.show("mah");
          //cout << "size=" << n << "deflsize=" << deflcoeff.size()<< "\n";
          //deflcoeff.show("arimah");
          
          qp.find_roots(qr, false, false);
          if constexpr (N > 0)
            {
              for (i=1; i < n; i++)
                roots[i] = qr[i-1];
            }
#if 1
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
              //refine_root(roots[i]);
              laguer_real(roots[i]);
              if (copyroot==true)
                {
                  roots[i+1] = conj(roots[i]);
                }
            }
        }
#endif
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
              solve_quadratic(roots);
            }
#if 0
          else if (n==3) 
            {
              solve_cubic_analytic(roots);
            }
          else if (n==4)
            {
              oqs_quartic_solver(roots);
            }
#endif
          else
            {
              laguerre_defl(roots, polish, backup);
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
              solve_quadratic(roots);
            }
#if 0
          else if constexpr (N==3) 
            {
              solve_cubic_analytic(roots);
            }
          else if constexpr (N==4)
            {
              oqs_quartic_solver(roots);
            }
#endif
          else 
            {
              laguerre_defl(roots, polish);
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
  cpoly()
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
      //std::cout << "max2= " << maxf2<< " max3=" << maxf3 << "\n";
      //printf("%.15G\n", maxf);
      //std::cout << "macheps= " << std::setprecision(35) << meps << "\n";
    }

  cpoly(int nc): cpolybase<cmplx,N>(nc)
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
    }
  //cpoly() = default;
#if 0
    {
      //coeff=v;
      //cpoly();
      //cout << "macheps= " << setprecision(35) << eps << "\n";
    }
  ~cpoly()
    {

    }
#endif
};
// quadratic equation
template<class cmplx, int N, class ntype> void cpoly<cmplx,N,ntype>::solve_quadratic(pvector<cmplx,N>&sol)
{
  cmplx acx,bcx,zx1,zx2,cdiskr,zxmax,zxmin;
  acx = coeff[1]/coeff[2];
  bcx = coeff[0]/coeff[2];
  cdiskr=sqrt(acx*acx-4.0*bcx);
  zx1 = -0.5*(acx+cdiskr);
  zx2 = -0.5*(acx-cdiskr);
  if (abs(zx1) > abs(zx2))
    zxmax = zx1;
  else
    zxmax = zx2;
  if (zxmax==cmplx(0.0,0.0))
    zxmin=0;
  else
    zxmin = bcx/zxmax;

  sol[0] = zxmin;
  sol[1] = zxmax;
}
#if 0
// cubic polynomial
template <class ntype, int N> void cpoly<ntype,N>::solve_cubic_analytic(pvector<cmplx,N>& sol)
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
template <class ntype, int N> void cpoly<ntype,N>::oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol)
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
template <class ntype, int N> void cpoly<ntype,N>::oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol)
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
template <class ntype, int N> void  cpoly<ntype,N>::oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled)
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
template <class ntype, int N> ntype  cpoly<ntype,N>::oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3)
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
template <class ntype, int N> 
ntype cpoly<ntype,N>::oqs_calc_err_abcd_cmplx(ntype a, ntype b,  ntype c, ntype d, cmplx aq, 
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
template <class ntype, int N> ntype cpoly<ntype,N>::oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq)
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
template <class ntype, int N> ntype  cpoly<ntype,N>::oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq)
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
template <class ntype, int N> void cpoly<ntype,N>::oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ)
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
template <class ntype, int N> void  cpoly<ntype,N>::oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2])
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

template <class ntype, int N> void cpoly<ntype,N>::oqs_quartic_solver(pvector<cmplx,N>& roots)
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
	  cdiskr=acx*acx/cmplx(4.0)-bcx;               
	  /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
	  zx1=-acx/cmplx(2.0)+sqrt(cdiskr);
	  zx2=-acx/cmplx(2.0)-sqrt(cdiskr);
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
#endif
