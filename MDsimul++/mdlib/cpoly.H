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

#define Sqr(x) ((x)*(x))
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
                   (imaxarg1) : (imaxarg2))
  template <class ntype>
inline ntype SIGN(ntype a, ntype b)
{
  if ((b) > 0)
    return abs(a);
  else
    return -abs(a);
}
//#define LAGUERRE_POLISH // migliora drasticamente con maxiter=8 ma per polinomi di alto grado (500) crea problemi
//#define BINI_CONV_CRIT
//#define USE_CONVEX_HULL 
#define USE_LAGUER_REAL
//#define PETKOVIC_6TH
using namespace std;
  template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template <class cmplx, class ntype, int N> 
class cpoly_base_static 
{
public:
  int n;
  constexpr static int dynamic = false;
  pvector<cmplx, N+1> coeff;
  pvector<cmplx, N+1> cmon;
  pvector<cmplx, N-1> abscmon;
  pvector<cmplx, N-1> alpha;
  pvector<cmplx, N> deflcoeff0;
  pvector<cmplx, N> deflcoeff1;
  pvector<cmplx, N> deflcoeff;
  void set_coeff(pvector<cmplx,N+1> v)
    {
      coeff = v;
    }
  cpoly_base_static()
    {
      n=N;
    }
  cpoly_base_static(int nc): cpoly_base_static()
  {
    n=nc;
  }
};
template <class cmplx, class ntype, int N> 
class cpoly_base_dynamic 
{
public:
  int n;
  constexpr static int dynamic = true;
  pvector<cmplx> coeff;
  pvector<cmplx> cmon;
  pvector<ntype> acmon;
  pvector<ntype> alpha;
  pvector<cmplx> deflcoeff0;
  pvector<cmplx> deflcoeff1;
  pvector<cmplx> deflcoeff;
  cpoly_base_dynamic() = default;
  void set_coeff(pvector<cmplx,-1> v)
    {
      coeff = v;
    }
  void use_vec(int nc, cmplx* coeffv, cmplx* cmonv,
               ntype* acmonv, ntype *alphav, cmplx *deflcoeff0v, cmplx* deflcoeff1v, cmplx* deflcoeffv)
    {
      coeff.use_vec(nc+1,coeffv);
      cmon.use_vec(nc+1,cmonv);
      acmon.use_vec(nc+1,acmonv);
      alpha.use_vec(nc+1,alphav);
      deflcoeff0.use_vec(nc,deflcoeff0v);
      deflcoeff1.use_vec(nc,deflcoeff1v);
      deflcoeff.use_vec(nc,deflcoeffv);
      n=nc;
    }

  cpoly_base_dynamic(int nc): coeff(nc+1), cmon(nc+1), acmon(nc+1), alpha(nc+1),
  deflcoeff0(nc), deflcoeff1(nc), deflcoeff(nc)
  {
    n=nc;
  }
  ~cpoly_base_dynamic() = default;
  void deallocate(void)
    {
      coeff.deallocate();
      cmon.deallocate();
      acmon.deallocate();
      alpha.deallocate();
      deflcoeff0.deallocate();
      deflcoeff1.deallocate();
      deflcoeff.deallocate();
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      cmon.allocate(n+1);
      acmon.allocate(n+1);
      alpha.allocate(n+1);
      deflcoeff0.allocate(n);
      deflcoeff1.allocate(n);
      deflcoeff.allocate(n);
    }
};

template <class cmplx, class ntype, int N> using cpolybase = 
typename std::conditional<(N>0), cpoly_base_static <cmplx, ntype, N>,
         cpoly_base_dynamic <cmplx, ntype, N>>::type;

template <class cmplx, int N=-1, class ntype=double> 
class cpoly: public numeric_limits<ntype>, public cpolybase<cmplx,ntype,N> 
{
  using cpolybase<cmplx,ntype,N>::n;
  using cpolybase<cmplx,ntype,N>::coeff;
  using cpolybase<cmplx,ntype,N>::cmon;
  using cpolybase<cmplx,ntype,N>::acmon;
  using cpolybase<cmplx,ntype,N>::alpha;
  using cpolybase<cmplx,ntype,N>::deflcoeff0;
  using cpolybase<cmplx,ntype,N>::deflcoeff1;
  using cpolybase<cmplx,ntype,N>::deflcoeff;
  template <class vtype>
    using pvecNm1 = typename std::conditional<(N>0),  pvector<vtype, N-1>,
          pvector<vtype, -1>>::type;
  template <class vtype>
    using pvecNp1 = typename std::conditional<(N>0),  pvector<vtype, N+1>,
          pvector<vtype, -1>>::type;
  using cpolyNm1 = typename std::conditional<(N>0),  cpoly<cmplx, N-1,ntype>,
        cpoly<cmplx, -1,ntype>>::type;

  const int maxiter_polish=8;
  int imaxarg1,imaxarg2;
  ntype eps05, meps, maxf, maxf2, maxf3, scalfact, cubic_rescal_fact;
  int maxdigits;
  ntype goaleps;
  ntype Kconv;
  inline void solve_quadratic(pvector<cmplx,N>& sol);
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

#if 1
  bool laguer_real(cmplx &xc)
    {
      const int MR=8,MT=10,MAXIT=MT*MR;
      const ntype EPS=goaleps;
      static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
      int iter,j;
      ntype abx,abp,abm,err, phi, rho,s, pr;
      ntype dx[2],x1[2],b[2],d[2],f[2],g[2],h[2],sq[2],gp[2],gm[2],g2[2], fb[2], pol[2];
      ntype x[2], f0, d0, b0, invn, g0, invgp[2], sarg[2];
      cmplx bc,dc,fc,sqrtcmplx, fovb, dovb;
      int m=n;
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
              return false;
            }

#ifdef BINI_CONV_CRIT
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
          invn = 1.0/(gp[0]*gp[0]+gp[1]*gp[1]);
          invgp[0] = gp[0]*invn;
          invgp[1] = -gp[1]*invn;

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
        }
      xc = cmplx(x[0],x[1]);
      return false;
    }
#endif
  bool laguer_cmplx(cmplx &x)
    {
      const int MR=8,MT=10,MAXIT=MT*MR;// MT changed from 10 to 30 
      const ntype EPS=goaleps;
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
          if (isinf(b.real()) || isinf(b.imag()) || isnan(b.real()) ||
              isnan(b.imag()))
            return false;
#ifdef BINI_CONV_CRIT
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
          dx=(MAX(abp,abm) > 0.0) ? ((cmplx)m)/((cmplx)gp) : polar(((ntype)(1+abx)),((ntype)iter));
          x1=x-dx;
          if (x == x1) 
            {
              return true;
            }

          if (iter % MT != 0) 
            x=x1;
          else 
            x -= frac[iter/MT]*dx;
        }
      return false;
    }
  void refine_root_maehly(cmplx &r0, pvector<cmplx,N> &roots, int iac)
    {
      int iter,i, j;
      ntype err;
      cmplx povp1, denom;
      // root polishing by NR
      cmplx r0old, p, sum;
      cmplx p1, p1c;
      ntype b[2], d[2], d0, b0, x[2];
      ntype errold;
      for (iter=0; ; iter++)
        {
#if 0
          b[0]=real(cmon[n]);
          b[1]=imag(cmon[n]);
          d[0]=0.0;
          d[1]=0.0;
          x[0] = real(r0);
          x[1] = imag(r0);
          for (j=n-1;j>=0;j--) {
            d0 = d[0];
            d[0] = x[0]*d[0] - x[1]*d[1] + b[0];
            d[1] = x[0]*d[1] + x[1]*d0 + b[1];
            b0 = b[0];
            b[0] = x[0]*b[0] - x[1]*b[1] + real(cmon[j]);
            b[1] = x[0]*b[1] + x[1]*b0 + imag(cmon[j]);
          }
          p=cmplx(b[0],b[1]);
          p1=cmplx(d[0],d[1]);
#else
          p=cmon[n]*r0+cmon[n-1];
          p1=cmon[n];
          for(i=n-2;i>=0;i--) {
            p1=p+p1*r0;
            p=cmon[i]+p*r0;
          }
#endif
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
            sum += ntype(1.0)/(r0-roots[i]);
          sum *= p;
          p1c = p1 - sum;
          r0 -= p/p1c;
          if (isnan(abs(r0)) || isinf(abs(r0)))
            {
              r0=r0old;
              break;
            }
        }
    }
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

          if (err==0)
            {
              break;
            }
          if (iter > 0 && err >= errold)
            {
              r0new[0]=r0old[0];
              r0new[1]=r0old[1];
              break;
            }
          if (p1[0]==0 && p1[1]==0)
            {
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
          //cout << "NR iter=" << iter << " err=" << err << "\n";
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
  bool laguer_cmplx_mod_rev(cmplx &xini, pvector<cmplx,N> roots, int iac)
    {
      const ntype EPS=goaleps;
      int j;
      ntype abx,abp,abm, M, s, pr, absb;
      cmplx term;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
#endif
      cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2, ba, xsq, x, ddb, db;
      int m;
      m=n;
      b=cmon[0];
      d=f=0.0;
      x=ntype(1.0)/xini;
      for (j=1;j<=m;j++) {
        f=x*f+d;
        d=x*d+b;
        b=x*b+cmon[j];
      }

#ifdef BINI_CONV_CRIT
      s=acmon[m];
      ba=cmon[m];
      abx=abs(xini);
      for (j=m-1; j >=0; j--) 
        {
          s=abx*s+acmon[j];
          ba=xini*ba+cmon[j];
        }
      absb=abs(ba);
      if (absb <= K*s) // stopping criterion of bini 
        {
          return true;
        }
#else
      err=alpha[0];//abs(cmon[0])*(Kconv*m+1);
      abx=1.0/abs(xini);
      for (j=1;j<=m;j++) {
        err=abx*err+alpha[j];//abs(cmon[j])*(Kconv*j+1);
      }
      if (abs(b) <= EPS*err) 
        {
          return true;
        }
#endif
      db=d/b;
      ddb=db*db-((ntype)2.0)*f/b;
      g=x*(ntype(m)-x*db);
      xsq=x*x;
      h=xsq*(ntype(m)-ntype(2.0)*x*db+xsq*ddb);

      for (j=0; j < n; j++)
        {
          if (j==iac)
            continue;
          term=ntype(1.0)/(xini-roots[j]);
          g -= term;
          h -= term*term; 
        }

      g2=g*g;
      sq=sqrt( ((ntype)(m-1))*(((ntype)m)*h-g2));
      gp=g+sq;
      gm=g-sq;
      abp=abs(gp);
      abm=abs(gm);

      if (abp < abm) gp=gm;
      dx =  ((cmplx)m)/((cmplx)gp);
      xini=xini-dx;
#if 0
      if (x1 == xini) 
        {
          return true;
        }
#endif
      return false;
    }
  bool laguer_real_mod_rev(cmplx &xc, pvector<cmplx,N> roots, int iac)
   {
     const ntype EPS=goaleps;
     int j;
     ntype abx,abp,abm,err, phi, rho,s, pr;
     ntype dx[2],x1[2],b[2],d[2],f[2],g[2],h[2],sq[2],gp[2],gm[2],g2[2], fb[2], pol[2];
     ntype xcR, xcI, x[2], f0, d0, b0, invn, g0, invgp[2], sarg[2], aR, aI, invden, term[2];
     ntype dovb[2], fovb[2], tt[2], ddb[2], tt2[2], ba[2];
     cmplx bc,dc,fc,sqrtcmplx;
     int m=n;
     xcR = xc.real();
     xcI = xc.imag();
     //cout << "xc=" << xc << "\n";
     invden = 1.0/(xcR*xcR+xcI*xcI);
     x[0] = xcR*invden;// 1/x
     x[1] = -xcI*invden;
     //its=iter;
     b[0]=real(cmon[0]);
     b[1]=imag(cmon[0]);
     d[0]=f[0]=0.0;
     d[1]=f[1]=0.0;
     for (j=1;j<=m;j++) {
       f0 = f[0];
       f[0] = x[0]*f[0] - x[1]*f[1] + d[0];
       f[1] = x[0]*f[1] + x[1]*f0 + d[1];
       d0 = d[0];
       d[0] = x[0]*d[0] - x[1]*d[1] + b[0];
       d[1] = x[0]*d[1] + x[1]*d0 + b[1];
       b0 = b[0];
       b[0] = x[0]*b[0] - x[1]*b[1] + real(cmon[j]);
       b[1] = x[0]*b[1] + x[1]*b0 + imag(cmon[j]);
     }

#ifdef BINI_CONV_CRIT
     s=acmon[m];
     abx = sqrt(xcR*xcR+xcI*xcI);
     ba[0]=real(cmon[m]);
     ba[1]=imag(cmon[m]);
     for (j=m-1; j >=0; j--) 
       {
         s=abx*s+acmon[j];
         b0 = ba[0];
         ba[0] = xcR*ba[0] - xcI*ba[1] + real(cmon[j]);
         ba[1] = xcR*ba[1] + xcI*b0 + imag(cmon[j]);
       }
     if (abs(cmplx(ba[0],ba[1])) <= 2.0*EPS*(4.0*ntype(n)+1)*s) // stopping criterion of bini 
       {
         //xc = cmplx(xcR, xcI);
         return true;
       }
#else
     //if (sqrt(b[0]*b[0]+b[1]*b[1]) <= EPS*err)
     err=alpha[0];//abs(cmon[0])*(Kconv*m+1);
     abx=1.0/sqrt(xcR*xcR+xcI*xcI);
     for (j=1;j<=m;j++) {
       err=abx*err+alpha[j];//abs(cmon[j])*(Kconv*j+1);
     }
     //cout << "abs=" << abs(cmplx(b[0],b[1])) << " EPS*err=" << EPS*err << "\n";
     if (abs(cmplx(b[0],b[1])) <= EPS*err)
       {
         //xc = cmplx(x[0],x[1]);
         return true;
       }
#endif          
     bc=cmplx(b[0],b[1]);
     dc=cmplx(d[0],d[1]);
     fc=cmplx(f[0],f[1]);
     //dcbc=dc/bc;
     //fcbc=fc/bc
     dovb[0]=real(dc/bc);
     dovb[1]=imag(dc/bc);
     fb[0]=real(fc/bc);
     fb[1]=imag(fc/bc);

     tt[0] = ntype(m)-x[0]*dovb[0]+x[1]*dovb[1];
     tt[1] = -x[0]*dovb[1]-x[1]*dovb[0];
     g[0] = x[0]*tt[0] - x[1]*tt[1];
     g[1] = x[1]*tt[0] + x[0]*tt[1];
     tt[0] = x[0]*x[0] - x[1]*x[1];
     tt[1] = 2.0*x[0]*x[1];
     ddb[0] = dovb[0]*dovb[0] - dovb[1]*dovb[1] - ntype(2.0)*fb[0];
     ddb[1] = 2.0*dovb[0]*dovb[1] - ntype(2.0)*fb[1];

     h[0] = tt[0]*ddb[0] - tt[1]*ddb[1];
     h[1] = tt[0]*ddb[1] + tt[1]*ddb[0]; 
     h[0] -= 2.0*(x[0]*dovb[0]-x[1]*dovb[1]);
     h[1] -= 2.0*(x[0]*dovb[1]+x[1]*dovb[0]);  
     h[0] += ntype(m);
     tt2[0] = h[0];
     tt2[1] = h[1];
     cmplx_mul(h,tt,tt2);
     for (j=0; j < n; j++)
       {
         if (j==iac)
           continue;
         aR = xcR-real(roots[j]);
         aI = xcI-imag(roots[j]);
         invden = 1.0/(aR*aR+aI*aI);
         term[0]=aR*invden;
         term[1]=-aI*invden;
         g[0] -= term[0];
         g[1] -= term[1];
         h[0] -= term[0]*term[0]-term[1]*term[1];
         h[1] -= 2.0*term[0]*term[1]; 
       }
     g2[0] = g[0]*g[0] - g[1]*g[1];
     g2[1] = 2.0*g[0]*g[1];
     sarg[0] = ((ntype)(m-1))*( ((ntype)m)*h[0] - g2[0]);
     sarg[1] = ((ntype)(m-1))*( ((ntype)m)*h[1] - g2[1]);
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
     invn = 1.0/(gp[0]*gp[0]+gp[1]*gp[1]);
     invgp[0] = gp[0]*invn;
     invgp[1] = -gp[1]*invn;

     dx[0]=((ntype)m)*invgp[0];
     dx[1]=((ntype)m)*invgp[1];
     //x1=x-dx;
     xc = cmplx(xcR - dx[0], xcI - dx[1]);

     //cout << "F xc=" << xc << "\n";
#if 0
     if (x[0] == x1[0] && x[1]==x1[1]) 
       {
         xc = cmplx(x[0],x[1]);
         return true;
       }
#endif
     return false;
   }

 bool laguer_real_mod(cmplx &xc, pvector<cmplx,N> roots, int iac)
   {
     const ntype EPS=goaleps;
     int j;
     ntype abx,abp,abm,err, phi, rho,s, pr;
     ntype dx[2],x1[2],b[2],d[2],f[2],g[2],h[2],sq[2],gp[2],gm[2],g2[2], fb[2], pol[2];
     ntype x[2], f0, d0, b0, invn, g0, invgp[2], sarg[2], aR, aI, invden, term[2];
     cmplx bc,dc,fc,sqrtcmplx, fovb, dovb;
     int m=n;
     x[0] = xc.real();
     x[1] = xc.imag();
     //its=iter;
     b[0]=real(cmon[m]);
     b[1]=imag(cmon[m]);
#ifndef BINI_CONV_CRIT
     err=alpha[m];
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
       b[0] = x[0]*b[0] - x[1]*b[1] + real(cmon[j]);
       b[1] = x[0]*b[1] + x[1]*b0 + imag(cmon[j]);
#ifndef BINI_CONV_CRIT
       err=abx*err+alpha[j];
#endif
     }

#ifdef BINI_CONV_CRIT
     s=acmon[n];
     for (j=n-1; j >=0; j--) 
       {
         s=abx*s+acmon[j];
       }

     if (abs(cmplx(b[0],b[1])) <= 2.0*EPS*(4.0*ntype(n)+1)*s) // stopping criterion of bini 
       {
         return true;
       }
#else
     //if (sqrt(b[0]*b[0]+b[1]*b[1]) <= EPS*err)
     if (abs(cmplx(b[0],b[1])) <= EPS*err)
       {
         return true;
       }
#endif          
#if 0
     cmplx_div(g,d,b);
     cmplx_div(fb,f,b);
     g2[0] = g[0]*g[0] - g[1]*g[1];
     g2[1] = 2.0*g[0]*g[1];
#else
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
#endif
     h[0] = g2[0] - 2.0*fb[0];
     h[1] = g2[1] - 2.0*fb[1];
     for (j=0; j < n; j++)
       {
         if (j==iac)
           continue;
         aR = x[0]-real(roots[j]);
         aI = x[1]-imag(roots[j]);
         invden = 1.0/(aR*aR+aI*aI);
         term[0]=aR*invden;
         term[1]=-aI*invden;
         g[0] -= term[0];
         g[1] -= term[1];
         h[0] -= term[0]*term[0]-term[1]*term[1];
         h[1] -= 2.0*term[0]*term[1]; 
       }
     g2[0] = g[0]*g[0] - g[1]*g[1];
     g2[1] = 2.0*g[0]*g[1];
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
     invn = 1.0/(gp[0]*gp[0]+gp[1]*gp[1]);
     invgp[0] = gp[0]*invn;
     invgp[1] = -gp[1]*invn;

     dx[0]=((ntype)m)*invgp[0];
     dx[1]=((ntype)m)*invgp[1];
     //x1=x-dx;
     x[0] = x[0] - dx[0];
     x[1] = x[1] - dx[1];
#if 0
     if (x[0] == x1[0] && x[1]==x1[1]) 
       {
         xc = cmplx(x[0],x[1]);
         return true;
       }
#endif
     xc = cmplx(x[0],x[1]);
     return false;
   }
  bool laguer_cmplx_mod(cmplx &x, pvector<cmplx,N> roots, int iac)
    {
      const ntype EPS=goaleps;
      cmplx term;
      //static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
      int j;
      ntype be,abx,abp,abm, M, s, pr, absb;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
#endif
      cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
      int m;
      m=n;
      b=cmon[m];
#ifndef BINI_CONV_CRIT
      err=alpha[m];//abs(cmon[m])*(Kconv*m+1);
#endif
      d=f=0.0;
      abx=abs(x);
      for (j=m-1;j>=0;j--) {
        f=x*f+d;
        d=x*d+b;
        b=x*b+cmon[j];
#ifndef BINI_CONV_CRIT
        err=abx*err + alpha[j];//abs(cmon[j])*(Kconv*ntype(j)+1.0);
#endif
      }
      absb=abs(b);
#ifdef BINI_CONV_CRIT
      s=acmon[n];
      for (j=n-1; j >=0; j--) 
        {
          s=abx*s+acmon[j];
        }
      if (absb <= K*s) // stopping criterion of bini 
        {
          return true;
        }
#else
      if (absb <= err*EPS) 
        {
          return true;
        }
#endif
      g=d/b;
      h=(g*g-((ntype)2.0)*f/b);
      for (j=0; j < n; j++)
        {
          if (j==iac)
            continue;
          term=ntype(1.0)/(x-roots[j]);
          g -= term;
          h -= term*term; 
        }

      g2=g*g;
      sq=sqrt( ((ntype)(m-1))*( ((ntype)m)*h-g2));
      gp=g+sq;
      gm=g-sq;
      abp=abs(gp);
      abm=abs(gm);

      if (abp < abm) gp=gm;
      dx =  ((cmplx)m)/((cmplx)gp);
      x=x-dx;
#if 0
      if (x == x1) 
        {
          return true;
        }
#endif
      return false;
    }

  ntype calc_upper_bound_kal(void)
    {
      // Kalantari's formula as found in McNamee Pan Vol. 1
      int k;
      static const ntype K= 1.0/0.682338;
      ntype ximax, xi;
      cmplx cnsq, term;

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
  using point=pair<ntype,ntype>; 

  typedef ntype coord_t;         // coordinate type
  typedef ntype coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

  struct Point {
    coord_t x;
    coord_t y;
    bool operator <(const Point &p) const {
      return x < p.x || (x == p.x && y < p.y);
    }
  };
  // 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
  // Returns a positive value, if OAB makes a clockwise turn,
  // negative for counter-clockwise turn, and zero if the points are collinear.
  coord2_t cross(const Point &O, const Point &A, const Point &B)
    {
      return ((A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x));
    }

  // Returns a list of points on the convex hull in counter-clockwise order.
  // Note: the last point in the returned list is the same as the first one.
  // Method: Andrew’s monotone chain algorithm O(n) in this case since we do not need sorting
  vector<Point> upper_convex_hull(vector<Point> P)
    {
      int n = P.size(), k = 0;
      if (n <= 3) 
        return P;
      vector<Point> H(n+1);
      // points are already ordered in the present case!
      //sort(P.begin(), P.end());
      // we need just upper hull for Bini's algorithm
      // Build upper hull
      for (int i = n-1; i >= 0; i--) {
        while (k >= 2 && cross(H[k-1], H[k], P[i]) < 0) 
          {
            k--;
          }
        k++;
        H[k] = P[i];
      }
      H.resize(k+1);
      return H;
    }
  void initial_guess(pvector<cmplx,N>& roots)
    {
#ifdef USE_CONVEX_HULL
      vector<Point> pts, convh;
#endif
      vector<cmplx> rg;
      vector<int> k;
      ntype sigma, ukc, vkc, maxt, t;
      Point p;
      int q, i, j;
#ifdef USE_CONVEX_HULL
      pts.resize(n+1);
      for (i = 0; i <= n; i++)
        {
          p.x=ntype(i);
          if (abs(cmon[i])!=0)
            {
              p.y=log(abs(cmon[i]));
            }
          else
            {
              p.y=-maxf/1.618034;
            }
          pts[i] = p;
        }
      convh=upper_convex_hull(pts);
      q=convh.size()-1;
      if (q<=2)
        {
          q=2;
          k.resize(q+1);
          k[1]=n;
          k[2]=0;
        }
      else
        {
          k.resize(q+1);
        }
      i=1; 
      for (i=1; i <= q; i++)
        {
          k[i] = ((int)convh[i].x);
        }
#endif
      ntype ukip1;
#ifndef USE_CONVEX_HULL
      vector<ntype> vk, uk;
      int kk;
      vk.resize(n+1);
      uk.resize(n+1);
      for (i=1; i <= n; i++) 
        {
          t=abs(cmon[i]/cmon[0]);
          if (i==1 || t > maxt)
            maxt=t;
        }
      uk[0] = 1.0/(1.0+maxt);   

      for (i=0; i < n; i++) 
        {
          t=abs(cmon[i]/cmon[n]);
          if (i==0 || t > maxt)
            maxt=t;
        }
      vk[n] = 1.0 + maxt;

      for (kk=1; kk <= n; kk++)
        {
          for (i=0; i < kk; i++) 
            {
              //if (cmon[kk]==cmplx(0.0,0.0))
                //continue;
              ukc=pow(abs(cmon[i]/cmon[kk]),1.0/(kk-i));
              if (i==0 || ukc > uk[kk])
                uk[kk]=ukc;
            }
        }
      for (kk=0; kk < n; kk++)
        {
          bool first=true;
          for (i=kk+1; i <= n; i++) 
            {
              if (cmon[i]==cmplx(0.0,0.0)) 
                continue;
              //if (cmon[kk]==cmplx(0.0,0.0))
                //continue;
              vkc=pow(abs(cmon[i]/cmon[kk]),1.0/(kk-i));
              if (first || vkc < vk[kk])
                {
                  first=false;
                  vk[kk]=vkc;
                }
            }
        } 
      k.resize(n+2); 
      q=1;
      for (kk=0; kk <= n; kk++)
        {
          if (uk[kk] <= vk[kk])
            {
              k[q] = kk; 
              q++;
            }
        }  
      q--;
      //sigma=M_PI+2.0*M_PI*(drand48()-0.5);
      sigma=2.0*M_PI*(drand48()-0.5);
      rg.resize(n); 
      int cc=0, dk;
      for (i=1; i < q; i++)
        {
          dk = k[i+1]-k[i];
          for (j=1; j <= dk; j++)
            {
              rg[(k[i]+j-1)] = polar(uk[k[i+1]],j*2.0*M_PI/dk+2.0*M_PI*ntype(i)/ntype(n)+sigma);
              cc++;
            }
        }
#else
      //sigma=2.0*M_PI*(drand48()-0.5)-M_PI;
      sigma=0.6;//2.0*M_PI*(drand48()-0.5)-M_PI;
      rg.resize(n); 
      int cc=0;

      for (i=q-1; i >= 1; i--)
        {
          ukip1 = pow(abs(cmon[k[i+1]]/cmon[k[i]]),1.0/(k[i]-k[i+1]));
          for (j=1; j <= k[i]-k[i+1]; j++)
            {
              rg[cc] = polar(ukip1,j*2.0*M_PI/(k[i]-k[i+1])+2.0*M_PI*ntype(i)/ntype(n)+sigma);
              cc++;
            }
        }
#endif
      //std::sort(rg.begin(),rg.end(), [&] (cmplx a, cmplx b)-> bool {return abs(a) > abs(b);});
      for (i=0; i < n; i++)
        {
          roots[i]=rg[i]; //dal più piccolo al più grande
          //cout << setprecision(20) << "(" << real(roots[i]) << "," << imag(roots[i]) << ")" << "r=" << abs(roots[i]) << "\n";
        }
    }
#ifdef PETKOVIC_6TH
  bool halley_mod_rev(cmplx &x, pvector<cmplx,N> roots, int iac)
    {
      const ntype EPS=goaleps;
      cmplx term;
      //static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
      int j;
      ntype be,abx,abp,abm, M, s, pr, absb;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
#endif
      cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
      int m;
      m=n;
      b=cmon[m];
#ifndef BINI_CONV_CRIT
      err=alpha[m];//abs(cmon[m])*(Kconv*m+1);
#endif
      d=f=0.0;
      abx=abs(x);
      for (j=m-1;j>=0;j--) {
        f=x*f+d;
        d=x*d+b;
        b=x*b+cmon[j];
#ifndef BINI_CONV_CRIT
        err=abx*err + alpha[j];//abs(cmon[j])*(Kconv*ntype(j)+1.0);
#endif
      }
      absb=abs(b);
#ifdef BINI_CONV_CRIT
      s=0.0;
      pr=1.0; 
      for (j=0; j <= n; j ++) 
        {
          s+=acmon[j]*pr;
          pr *= abx;
        }

      if (absb <= K*s) // stopping criterion of bini 
        {
          return true;
        }
#else
      if (absb <= err*EPS) 
        {
          return true;
        }
#endif
      g=d/b;
      h=(g*g-((ntype)2.0)*f/b);
      for (j=0; j < n; j++)
        {
          if (j==iac)
            continue;
          term=ntype(1.0)/(x-roots[j]);
          g -= term;
          h -= term*term; 
        }

      g2=g*g;
      sq=sqrt( ((ntype)(m-1))*( ((ntype)m)*h-g2));
      gp=g+sq;
      gm=g-sq;
      abp=abs(gp);
      abm=abs(gm);

      if (abp < abm) gp=gm;
      dx =  ((cmplx)m)/((cmplx)gp);
      x=x-dx;
#if 0
      if (x == x1) 
        {
          return true;
        }
#endif
      return false;
    }
  bool halley_mod(cmplx &x, pvector<cmplx,N> roots, int iac)
    {
      const ntype EPS=goaleps;
      cmplx term;
      //static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
      int j;
      ntype be,abx,abp,abm, M, s, pr, absb;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
#endif
      cmplx dx,x1,b,d,f,g,h,sq,gp,gm,g2;
      int m;
      m=n;
      b=cmon[m];
#ifndef BINI_CONV_CRIT
      err=alpha[m];//abs(cmon[m])*(Kconv*m+1);
#endif
      d=f=0.0;
      abx=abs(x);
      for (j=m-1;j>=0;j--) {
        f=x*f+d;
        d=x*d+b;
        b=x*b+cmon[j];
#ifndef BINI_CONV_CRIT
        err=abx*err + alpha[j];//abs(cmon[j])*(Kconv*ntype(j)+1.0);
#endif
      }
      absb=abs(b);
#ifdef BINI_CONV_CRIT
      s=0.0;
      pr=1.0; 
      for (j=0; j <= n; j ++) 
        {
          s+=acmon[j]*pr;
          pr *= abx;
        }

      if (absb <= K*s) // stopping criterion of bini 
        {
          return true;
        }
#else
      if (absb <= err*EPS) 
        {
          return true;
        }
#endif
      g=d/b;
      h=(g*g-((ntype)2.0)*f/b);
      for (j=0; j < n; j++)
        {
          if (j==iac)
            continue;
          term=ntype(1.0)/(x-roots[j]);
          g -= term;
          h -= term*term; 
        }

      g2=g*g;
      sq=sqrt( ((ntype)(m-1))*( ((ntype)m)*h-g2));
      gp=g+sq;
      gm=g-sq;
      abp=abs(gp);
      abm=abs(gm);

      if (abp < abm) gp=gm;
      dx =  ((cmplx)m)/((cmplx)gp);
      x=x-dx;
#if 0
      if (x == x1) 
        {
          return true;
        }
#endif
      return false;
    }
  void petkovic6th(pvector<cmplx,N>& roots)
    {
      bool *found = new bool[n], fine=false, ret;
      int i,  nf=0, iter;
      int itmax=100;
      pvector<cmplx> rold;
      ntype r;
      rold.allocate(n);
      if (coeff[n]==cmplx(0.0,0.0))
        {
          cout << "That's not an " << n << "degree polynomial!\n";
          return;
        }
      cmon[n]=1.0;
      for (i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }
      for (i=0; i <=n; i++)
        {
          acmon[i] = abs(cmon[i]);
          alpha[i] = acmon[i]*Kconv*(ntype(i)+ntype(1.0));
        }
      //absolute value of coefficients are used in Bini's stopping criterion
      //we use coeff vector to store them since it won't be used anymore from here on

      for (i=0;i < n; i++)
        found[i]=false;
      initial_guess(roots);
      for (iter=0; iter < itmax && !fine; iter++)
        {
          rold=roots;
          for (i=0; i < n; i++)
            {
              // refine conjuigate pairs only once!
              if (found[i]) 
                continue;
              r=abs(roots[i]);
              if (r > 1)
                ret = halley_mod_rev(roots[i], rold, i);
              else
                ret = halley_mod(roots[i], rold, i);
              if (ret)
                {
                  found[i]=true;
                  nf++;
                  if (nf==n)
                    {
                      fine=true;
                      break;
                    }
                }
            }
        }
      if (nf < n)
        {
          cout << "[ERROR] We had a problem, I was not able to converge to all roots\n";
          cout << "Only " << nf << " out of " << n << " were found\n";
          cout << "Try to increase floating-point precision\n";
          for (i=0; i <= n; i++)
            {
              cout << setprecision(80) << real(cmon[i]) << "\n";
            }
          exit(1);
        }
      //cout << "iter=" << iter << " # found=" << nf << "\n";
      delete []found;
    }
#endif
  void laguerre_aberth(pvector<cmplx,N>& roots, bool polish=false)
    {
      bool *found = new bool[n], fine=false, ret;
      int i,  nf=0, iter, is;
      int itmax=1000, nold;
      //pvector<cmplx> rold;
      ntype r;
      cmplx cn;
      //rold.allocate(n);
      if (coeff[n]==cmplx(0.0,0.0))
        {
          cout << "That's not an " << n << " degree polynomial!\n";
          return;
        }
      i=0;
      cn = coeff[n];
      // if first is coefficients are zero we x=0 root with multiplicity is
      nold=n;
      while (coeff[i]==cmplx(0.0,0.0))
        {
          roots[n-1]=0;
          n--;
          i++;
        }
      is=i;
      //cout << "n=" << n << " is=" << is << "\n";
      cmon[n]=1.0;
      for (i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i+is]/cn;
          //cout << "cmon[" << i << "]=" << cmon[i] << "\n";
        }
      for (i=0; i <=n; i++)
        {
          acmon[i] = abs(cmon[i]);
          alpha[i] = acmon[i]*Kconv*(ntype(i)+ntype(1.0));
        }
      //absolute value of coefficients are used in Bini's stopping criterion
      //we use coeff vector to store them since it won't be used anymore from here on

      for (i=0;i < n; i++)
        found[i]=false;
      initial_guess(roots);
      for (iter=0; iter < itmax && !fine; iter++)
        {
          //rold=roots;
          for (i=0; i < n; i++)
            {
              // refine conjuigate pairs only once!
              if (found[i]) 
                continue;
              r=abs(roots[i]);
#ifdef USE_LAGUER_REAL
              if (r > 1)
                ret = laguer_real_mod_rev(roots[i], roots, i);
              else
                ret = laguer_real_mod(roots[i], roots, i);
#else
              if (r > 1)
                ret = laguer_cmplx_mod_rev(roots[i], roots, i);
              else
                ret = laguer_cmplx_mod(roots[i], roots, i);
#endif
              if (ret)
                {
                  found[i]=true;
                  nf++;
                  if (nf==n)
                    {
                      fine=true;
                      break;
                    }
                }
            }
        }
      //cout << "iter=" << iter << "\n";
      // refine 
      if (polish==true)
        {
          for (i=0; i < n; i++)
            refine_root_maehly(roots[i], roots, i);
        }
      if (nf < n)
        {
#if 1
          cout << "[LAGABERTH-ERROR] We had a problem, I was not able to converge to all roots\n";
          cout << "Only " << nf << " out of " << n << " were found\n";
          cout << "Try to increase floating-point precision\n";
          for (i=0; i <= n; i++)
            {
              cout << setprecision(80) << real(cmon[i]) << "\n";
            }
          exit(1);
#endif
        }
      n=nold;
      //cout << "iter=" << iter << " # found=" << nf << "\n";
      delete []found;
    }
  void laguerre_defl(pvector<cmplx,N>& roots,bool polish=false, bool backup=false)
    {
      // find initial guess for x here
      cmplx r0;
      const ntype goldfact=1.618034;
      ntype sig, sigold;
      //bool mulsig=true;
      ntype err0,err1, errA, errB;
      pvector<cmplx,3> ad;
      typename std::conditional<(N < 0),pvector<cmplx>,bool>::type cbak;
      int iter, i;
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
      if constexpr (N < 0)
        {
          cbak.allocate(n+1);
          if (backup||polish)
            {
              cbak = coeff;
            }
        }
      r0=cmplx(0.0,0.0);
      sigold=0.0;
      sig=1.0;

      int MAXSTPS=1000,MT=100;
      while (!laguer_real(r0))
        {
          if (iter==0)
            r0=cmplx(0.0,1.0);
          else
            //r0 = polar(1.0,ntype(2.0*M_PI*(drand48()-0.5)));
            //r0 = polar(ntype(drand48()),ntype(2.0*M_PI*(drand48()-0.5)));
            r0 = polar(sigold+(sig-sigold)*drand48(),ntype(2.0*M_PI*(drand48()-0.5)));
          iter++;

          if (iter%MT==0)
            {
              sigold=sig;
              sig*=goldfact;
            }
          if (iter > MAXSTPS)
            {
              cout << "[ERROR] Laguerre method did not converge:\n";
              cout << "this problem can be fixed by using floating point number\n";
              cout << "with higher precision (e.g. long double and complex<long double>\n";
              exit(-1); 
            }
        };

#ifdef LAGUERRE_POLISH
      if (is_same<ntype,double>::value || is_same<ntype,long double>::value)
        {
          refine_root(r0);
        }
#endif
      cpolyNm1 qp;
      // pvector<cmplx,N-1> qr;
      pvecNm1<cmplx> qr;
      //cout << "minus one r0=" << r0 << "\n";
      roots[0] = r0;
      for (i=n; i >=1; i--)
        {
          deflcoeff0[i-1] = cmon[i];
          if (i <= n-1)
            deflcoeff0[i-1] += r0*deflcoeff0[i];
        }
      deflcoeff = deflcoeff0;
      if constexpr (N < 0)
        {
          //qp.allocate(n-1);
          //qr.allocate(n-1);
          qr.use_vec(n-1,&(roots[1])); 
          qp.use_vec(n-1,coeff.v, cmon.v, acmon.v, alpha.v, deflcoeff0.v, deflcoeff1.v, deflcoeff.v);
        }

      qp.set_output_prec(this->get_output_prec());
      qp.set_coeff(deflcoeff);

      qp.find_roots_defl(qr, false, false);
      if constexpr (N > 0)
        {
          for (i=1; i < n; i++)
            roots[i] = qr[i-1];
        }
      if constexpr (N < 0)
        {
          if (backup||polish)
            {
              coeff = cbak;
            }
        }
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
              refine_root_maehly(roots[i], roots,  (i==0)?0:i-1);
            }
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
  inline void find_roots(pvector<cmplx,N>& roots)
    {
      laguerre_aberth(roots);
    }
  inline void find_roots_defl(pvector<cmplx,N>& roots, bool polish=false, bool backup=true)
    {
      // if backup is true the dynamic version preserve the coefficients after calling this method
      // solve by deflation as zroots routine of Numerical Recipe
      // This algorithm is not backward stable hence not usable for n > 50-60
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
          else 
            {
              laguerre_defl(roots, polish);
            }
        }
    }
  // get machine precision for "ntype" type (ntype can float, double, long double)
  ntype epsilon()
    {
      return numeric_limits<ntype>::epsilon(); 
    }
  ntype getmax()
    {
      return numeric_limits<ntype>::max();
    }
  void set_output_prec(ntype e)
    {
      goaleps=e;
    } 
  ntype get_output_prec(void)
    {
      return goaleps;
    } 

  void init_const(void)
    {
      meps = epsilon();
      eps05 = pow(numeric_limits<ntype>::epsilon(),0.5);
      maxf= getmax();
      maxdigits = numeric_limits<ntype>::digits10-1;
      maxf2 = pow(maxf,0.5)/10.0;
      maxf3 = pow(maxf,1.0/3.0)/10.0;
      scalfact = pow(maxf,1.0/4.0)/1.618034;
      cubic_rescal_fact = pow(maxf, 1.0/3.0)/1.618034;
      goaleps=numeric_limits<ntype>::epsilon(); 
      //Kconv=sqrt(ntype(2.0))*ntype(2.0)+ntype(1.0);  
      Kconv=3.8;
    }
  cpoly()
    {
      init_const();  
    }

  cpoly(int nc): cpolybase<cmplx,ntype,N>(nc)
  {
    init_const();
  }
};
// quadratic equation
template<class cmplx, int N, class ntype> void cpoly<cmplx,N,ntype>::solve_quadratic(pvector<cmplx,N>&sol)
{
  cmplx acx,bcx,zx1,zx2,cdiskr,zxmax,zxmin;
  acx = coeff[1]/coeff[2];
  bcx = coeff[0]/coeff[2];
  cdiskr=sqrt(acx*acx-ntype(4.0)*bcx);
  zx1 = -ntype(0.5)*(acx+cdiskr);
  zx2 = -ntype(0.5)*(acx-cdiskr);
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
