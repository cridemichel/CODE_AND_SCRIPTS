// matrici NxM
//#include<stdio.h>
//#include<stdlib.h>
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<complex>
//#include "pvector.H"
#include "pmatrix.H"
#include <limits>
#include <cstdlib>
#include <array>
#define Sqr(x) ((x)*(x))
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? abs(a) : -abs(a))
using namespace std;
template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}
template <class ntype, int N, bool force_hqr=false> 
class rpoly: public numeric_limits<ntype> {
  pvector<ntype, N+1> coeff;
  pvector<ntype, N> cmon;
  pvector<ntype, N> cmonder;
  pvector<ntype, N-1> cmondder;
  int n;
  //force_hqr;
  int imaxarg1,imaxarg2;
  ntype meps, maxf, maxf2, maxf3, scalfact, cubic_rescal_fact;
  int maxdigits;

  //ntype (rpoly<ntype,N>::*func) (ntype x);
  ntype zbrent(const ntype x1, const ntype x2, const ntype tol)
  {
    const int ITMAX=100;
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
          b += SIGN(tol1,xm);
        fb=evalpoly(b);
      }
    cout << "Maximum number of iterations exceeded in zbrent";
    exit(-1);
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
  //void oqs_quartic_solver(pvector<complex<ntype>,N>& roots);
  inline void solve_quadratic(pvector<complex<ntype>,N>& sol);
  inline void solve_cubic_analytic(pvector<complex<ntype>,N>& sol);
  inline void oqs_quartic_solver(pvector<complex<ntype>,N>& roots);
  inline void oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol);
  inline void oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol);
  inline void oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled);
  inline ntype oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3);
  inline ntype oqs_calc_err_abcd_cmplx(ntype a, ntype b, ntype c, ntype d, 
	    			       complex<ntype> aq, complex<ntype> bq, complex<ntype> cq, complex<ntype> dq);
  inline ntype oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq);
  inline ntype oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq); 
  inline void oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ);
  inline void oqs_solve_quadratic(ntype a, ntype b, complex<ntype> roots[2]);
  inline void balance(ntype a[N][N]);
  inline void hqr(ntype a[N][N], pvector<complex<ntype>,N>& wri, int& ok);
  inline void QRfactorization(ntype hess[N][N], pvector<complex<ntype>,N>& sol, int& ok);
  inline void solve_numrec(pvector<complex<ntype>,N>& csol, int& ok);
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
	  if (coeff[i] > 0)
    	    {
	      if (i < N)
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
  void set_coeff(pvector<ntype,N+1> v)
    {
      coeff = v;
    }

  ntype evalpoly(ntype x)
    {
      // evaluate polynomail via Horner's formula 
      ntype bn=0.0;
      for (int i=N; i >= 0; i--)
        {
          bn = cmon[i] + bn*x;
        }
      return bn;
    }
  ntype evaldpoly(ntype x)
    {
      // evaluate first derivative of polynomail via Horner's formula 
      ntype bn=0.0;

      for (int i=N-1; i >= 0; i--)
        {
          bn = (i+1)*cmonder[i+1] + bn*x;
        }
      return bn;
    }
  ntype evalddpoly(ntype x)
    {
      // evaluate second derivative of polynomail via Horner's formula 
      ntype bn=0.0;
      if (N == 1)
        return 0;
      for (int i=N-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmondder[i+2] + bn*x;
        }
      return bn;
    }
  ntype calc_err_q(pvector<ntype,N> c, ntype r0)
    {
      int i;
      ntype sum;
      sum = (cmon[N-1]==0.0)?abs((c[N-2] - r0)-cmon[N-1]):abs(((c[N-2] - r0)-cmon[N-1])/cmon[N-1]);
      for (i=N-2; i >=0; i--)
        {
          sum += (cmon[i]==0)?abs((c[i-1] - r0)-cmon[i]):abs(((c[i-1] - r0)-cmon[i])/cmon[i]);
        }
      sum += (cmon[0]==0.0)?abs(r0*c[0]+cmon[0]):abs((r0*c[0]+cmon[0])/cmon[0]);
      return sum;
    }
#if 0
  void use_NR(pvector<ntype,N>& defc, ntype r0)
    {
      ntype errf, errfold, det;
      pvector<ntype,N-1> c, cold;
      pvector<ntype,N-1> fvec;
      pmatrixq<ntype,N-1> Jinv;
      int i, j, iter;
      for (i=N-2; i >= 0; i--)
        {
          c[i]=defc[i];
          fvec[i] = defc[i]-r0*((i==N-2)?1:defc[i+1])-cmon[i+1];
        } 
      //cmon.show("monic coeff=");
      //defc.show("deflcoeff");
      //fvec.show("errori=");
      errf = 0;
      for (i=0; i < N-1; i++)
        {
          errf += abs((cmon[i+1]==0)?fvec[i]:(fvec[i]/cmon[i+1]));
        }
      //cout << "errf inizio=" << errf << "\n";
      for (i=0; i < N-1; i++)
        for (j=0; j < N-1; j++)
          Jinv[i][j] = (i==j)?1.0:0.0;
      for (j=0; j < N-2; j++)
        for (i=j+1; i < N-1; i++)
          {
            Jinv[i][j]=r0*Jinv[i-1][j];
          }
      // Jinv.show();
      for (iter = 0; iter < 8; iter++)
        {
          cold = c;
          c = c - Jinv*fvec;
          for (i=N-2; i >= 0; i--)
            {
              fvec[i] = c[i]-r0*((i==N-2)?1:defc[i+1])-cmon[i+1];
            } 
          errfold = errf;
          errf = 0;
          for (i=0; i < N-1; i++)
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
      for (i=0; i < N-1; i++)
        defc[i] = c[i];
    }
#endif
  void use_NR_defl_lin(pvector<ntype,N>& bd, ntype r0)
    {
      ntype term, errf, errfold;
      pvector<ntype,N-1> c, cold;
      pvector<ntype,N-1> fvec;
      pmatrixq<ntype,N-1> Jinv, J;
      int i, j, iter;
      ntype sum=0;

      for (i=0; i < N-1; i++)
        c[i] = bd[i];
      for (i=1; i < N-1; i++)
        {
          if (i < N-1)
            fvec[i-1] = -c[i]*r0;
          fvec[i-1] += c[i-1];
          fvec[i-1] -= cmon[i];
        }
      //cmon.show("monic coeff=");
      //defc.show("deflcoeff");
      //fvec.show("errori=");
      errf = 0;
      for (i=1; i < N-1; i++)
        {
          errf += abs((cmon[i]==0)?fvec[i-1]:(fvec[i-1]/cmon[i]));
        }
      for (i=0; i < N-1; i++)
        for (j=0; j < N-1; j++)
          Jinv[i][j] = (i==j)?1.0:0.0;
      for (i=N-3; i >= 0; i--)
        for (j=i+1; j < N-1; j++)
          {
            Jinv[i][j]=r0*Jinv[i+1][j];
          }
      //cout << "errf inizio=" << errf << "\n";
      // Jinv.show();
      for (iter = 0; iter < 8; iter++)
        {
          cold = c;
          c = c - Jinv*fvec;
          for (i=1; i < N-1; i++)
            {
              if (i < N-1)
                fvec[i-1] = -c[i]*r0;
              fvec[i-1] += c[i-1];
              fvec[i-1] -= cmon[i];
            }
          errfold = errf;
          for (i=1; i < N-1; i++)
            {
              errf += abs((cmon[i]==0)?fvec[i-1]:(fvec[i-1]/cmon[i]));
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
      for (i=0; i < N-1; i++)
        bd[i] = c[i];
      
    }
  void use_NR_defl_quad(pvector<ntype,N-1>& bd, pvector<ntype,3> ad)
    {
      // bd[N-2] = 1 since all polynomials are in their minic form
      // hence we do not need to refine this coefficient here
      ntype term, errf, errfold, det;
      pvector<ntype,N-2> cold, c;
      pvector<ntype,N-2> fvec;
      pmatrixq<ntype,N-2> Jinv, J;
      int i, j, iter;
      ntype sum=0;
      for (i=0; i < N-2; i++)
        c[i]=bd[i];
      for (i=2; i < N; i++)
        {
          if (i < N-2)
            fvec[i-2] = -c[i]*ad[0];
          if (i < N-1)
            fvec[i-2] -= c[i-1]*ad[1];
          fvec[i-2] -= c[i-2];
          fvec[i-2] += cmon[i];
        }
      //cmon.show("monic coeff=");
      //defc.show("deflcoeff");
      //fvec.show("errori=");
      errf = 0;
      for (i=2; i < N; i++)
        {
          errf += abs((cmon[i]==0)?fvec[i-2]:(fvec[i-2]/cmon[i]));
        }
      //cout << "errf inizio=" << errf << "\n";
      for (i=0; i < N-2; i++)
        for (j=0; j < N-2; j++)
          {
            if (i==j)
              J[i][j] = 1.0;
            else if (i+1==j)
              J[i][j] = ad[1];
            else if (i+2==j)
              J[i][j] = ad[0];
            else 
              J[i][j] = 0.0;
          }
      Jinv = J.inv();
      // Jinv.show();
      for (iter = 0; iter < 8; iter++)
        {
          cold = c;
          c = c - Jinv*fvec;
          for (i=2; i < N; i++)
            {
              if (i < N-2)
                fvec[i-2] = -c[i]*ad[0];
              if (i < N-1)
                fvec[i-2] -= c[i-1]*ad[1];
              fvec[i-2] -= c[i-2];
              fvec[i-2] += cmon[i];
            }

          errfold = errf;
          for (i=0; i < N-2; i++)
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
      for (i=0; i < N-2; i++)
        bd[i]=c[i];
    }
 bool laguer(complex<ntype> &x)
    {
      const int MR=8,MT=10,MAXIT=MT*MR;
      const ntype EPS=numeric_limits<ntype>::epsilon();
      static const ntype frac[MR+1]= {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
      int iter,j;
      ntype abx,abp,abm,err;
      complex<ntype> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
      int m=N;
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
          if (abs(b) <= err) return true;
          g=d/b;
          g2=g*g;
          h=g2-2.0*f/b;
          sq=sqrt(ntype(m-1)*(ntype(m)*h-g2));
          gp=g+sq;
          gm=g-sq;
          abp=abs(gp);
          abm=abs(gm);
          if (abp < abm) gp=gm;
          dx=MAX(abp,abm) > 0.0 ? ntype(m)/gp : polar(1+abx,ntype(iter));
          x1=x-dx;
          if (x == x1) 
            return true;
          if (iter % MT != 0) 
            x=x1;
          else 
            x -= frac[iter/MT]*dx;
        }
      //cout << "too many iterations in laguer\n";
      return false;
    }
  ntype calc_err_defl_quad(pvector<ntype,N-1>& bd, pvector<ntype,3>& ad)
    {
      int i;
      ntype sum=0, term;
      for (i=0; i < N; i++)
        {
          term=cmon[i];
          if (i <= N-2) 
            term -= bd[i]*ad[0];
          if (i >= 1 && i <= N-1)
            term -= bd[i-1]*ad[1];
          if (i >= 2)
            term -= bd[i-2]*ad[1];
          if (cmon[i]==0)
            sum+=abs(term);
          else
            sum+=abs(term/cmon[i]);
        }
      return sum;
    }
  void deg6roots(pvector<complex <ntype>,N>& roots)
    {
      // find initial guess for x here
      complex<double> r0;
      const ntype goldfact=1.618034;
      ntype sig;
      pvector<ntype,N> deflcoeff0, deflcoeff1, deflcoeff;
      ntype err0,err1;
      int iter; 
      int i;
      if (coeff[N]==0.0)
        {
          cout << "That's not a quintic!\n";
          return;
        }
      cmon[N]=1.0;
      for (i=N-1; i >=0; i--)
        cmon[i]=coeff[i]/coeff[N];

      r0 = 0.0;
      iter = 0;
      sig=1.0;
      while (!laguer(r0))
        {
          r0 = sig*complex<ntype>(drand48()-0.5,drand48()-0.5);
          sig *= goldfact; 
          iter++;
          if (iter > 1000)
            {
              cout << "Laguerre can not converge, I am sorry...\n";
              exit(-1); 
            }
        };
      roots[0] = r0;
      // la forward deflation è sbagliata nel secondo caso!! Fix it!!!
      if (r0.imag()==0)
        {
          // deflate by (x-r0) and call quintic solver
          // backward deflation
          for (i=N; i >=1; i--)
            {
              //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff[" << i+1 << "]=" << deflcoeff[i+1]<< "\n";
              deflcoeff0[i-1] = cmon[i];
              if (i <= N-1)
                deflcoeff0[i-1] += r0.real()*deflcoeff0[i];
             }
          if (r0.real()!=0.0)
            {
              //forward deflation
              //cout << "cmon[0]=" << cmon[0] << " deflcoeff1[0]=" << deflcoeff1[0]<< "\n";
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
              deflcoeff1[0]=-cmon[0]/r0.real();
              deflcoeff1[N-1] = 1.0;
              //cout << "cmon[0]=" << cmon[0] << " deflcoeff1[0]=" << deflcoeff1[0]<< "\n";
              for (i=1; i < N-1; i++)
                {
                  deflcoeff1[i] = (deflcoeff1[i-1] - cmon[i])/r0.real();
                  //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff1[" << i << "]=" << deflcoeff1[i]<< "\n";
                }
#endif
              //deflcoeff1.show();
              // choose the set which yields smaller error 
              err0=calc_err_q(deflcoeff0,r0.real());
              err1=calc_err_q(deflcoeff1,r0.real());
              if (err0 < err1)
                deflcoeff = deflcoeff0;
              else
                deflcoeff = deflcoeff1;
            }
          else
            deflcoeff = deflcoeff0;
          // refine deflated quartic coefficients here
          //deflcoeff.show("prima");
          //pvector<ntype,N-1> bd;
          //deflcoeff=deflcoeff0;
          use_NR_defl_lin(deflcoeff,r0.real());
          rpoly<ntype,N-1> qp;
          pvector<complex<ntype>,N-1> qr;
          qp.set_coeff(deflcoeff);
          qp.find_roots(qr);
          for (i=1; i < N; i++)
            roots[i] = qr[i-1];
          //roots.show();
        }
      else
        {
          // qui ci sono ancora casini e non funziona
          // deflate by (x-r0)*(x-conj(r0)) and call quartic solver
          roots[1] = conj(r0);
          cout << "lag roots=" << roots[1] << " " << roots[0] << "\n";
          rpoly<ntype,N-2> rp;
          pvector<ntype,3> ad;
          pvector<ntype,N-1> bd, bd0, bd1;
          pvector<ntype,N-1> ro;
          // find coefficients for deflation
          ad[2] = 1.0;
          ad[1] = -2*roots[0].real();
          ad[0] = abs(roots[0])*abs(roots[0]);
          //backward deflation 
          for (i=N; i >=2;  i--)
            {
              bd0[i-2] = cmon[i];
              if (i <= N-1)
                bd0[i-2] -= ad[1]*bd0[i-1];
              if (i <= N-2)
                bd0[i-2] -= ad[0]*bd0[i];
            }
          if (ad[0]!=0)
            {
              //forward deflation
              bd1[N-2]=1.0;
              for (i=0; i < N-2; i++)
                {
                  bd1[i] = cmon[i];
                  if (i >= 1)
                    bd1[i] -= bd[i-1]*ad[1];
                  if (i >= 2)
                    bd1[i] -= bd[i-2];
                  bd1[i] /= ad[0];
                }
            }
          err0=calc_err_defl_quad(bd0, ad);
          err1=calc_err_defl_quad(bd1, ad);
          if (err0 < err1)
            {
              bd = bd0;
            }
          else
            {
              bd = bd1;
            }
          bd=bd0;
          pvector<complex<ntype>,4> ro_quad;
          use_NR_defl_quad(bd,ad);
          if (bd[4]==0.0)
            {
              bd.show("bd=");
              bd0.show("bd0=");
              bd1.show("bd1=");
              exit(-1);
            }
          rp.set_coeff(bd);
          rp.find_roots(ro_quad);
          ro_quad.show("quartic roots=");
          for (i=2; i < N; i++)
            roots[i] = ro_quad[i-2];
        }
      //cout << "x=" << r0 << "\n";
    } 
  void quintic(pvector<complex <ntype>,N>& roots)
    {
      ntype coeffder[N-1];
      rpoly<ntype, N-1> quar;
      pvector<complex<ntype>,N-1> qroots;
      pvector<ntype,N> deflcoeff0, deflcoeff1, deflcoeff;
      complex<ntype> quadr[2];
      std::array<ntype,N-1> rqroots;
      ntype goldfact = 1.618034;
      ntype delx, c;
      int i, cc;
      if (coeff[N]==0.0)
        {
          cout << "That's not a quintic!\n";
          return;
        }
      cmon[N]=1.0;
      for (i=N-1; i >=0; i--)
        cmon[i]=coeff[i]/coeff[N];

      cmonder[N-1] = N;
      for (i = N-2; i >= 0; i--)
        {
          cmonder[i] = (i+1)*cmon[i+1]; 
        }
      quar.set_coeff(cmonder);

      quar.find_roots(qroots);
      //quar.show("quartic=");
      //qroots.show("radici quartica");
      cc=0;
      for (i=0; i < N-1; i++)
        if (qroots[i].imag()==0.0 && evalddpoly(qroots[i].real())!=0)
          rqroots[cc++] = qroots[i].real();
      int numr = cc;
      //qroots.show();
      //qsort(rqroots, numr, sizeof(ntype),
      //    [&] (ntype *x, ntype *y){if (*x < *y) return -1; else return 1;});
      // sort(rqroots.begin(),rqroots.end());
      if (rqroots[1] < rqroots[0])
        {
          swap(rqroots[0],rqroots[1]);
        }
      //cout << "----> "<< rqroots[0] << " "<< rqroots[1] << "\n";
      bool foundleft = false;
      bool foundright = false;
      ntype xleft, xright;
      if (numr == 0)
        {
          rpoly<ntype,N-2> cub;
          ntype r0;
          pvector<complex<ntype>,N-2> cubroots;
          cmondder[N-2] = N*(N-1);
          //cout << "qui\n";
          for (i = N-3; i >= 0; i--)
            {
              cmondder[i] = (i+1)*(i+2)*cmon[i+2]; 
            }
          cub.set_coeff(cmondder);
          //cmondder.show();
          cub.find_roots(cubroots); 
          for (i=0; i < N-2; i++)
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
              xleft = xright-pow(c,1.0/5.0);
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
              xright = xleft+pow(-c,1.0/5.0);
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
              xleft = xright-pow(c,1.0/5.0);
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
              xright = xleft+pow(-c,1.0/5.0);
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
                if (evalpoly(rqroots[i])*evalpoly(rqroots[i+1]) < 0.0)
                  {
                    xleft = rqroots[i];
                    xright = rqroots[i+1];
                    foundleft = true;
                    foundright = true;
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
      ntype r0 = zbrent(xleft, xright, numeric_limits<ntype>::epsilon());
      roots[0] = r0;
      //cout << "r0=" << r0 << "\n";
      // deflate quintic and use OQS
      // first set of coefficients
      deflcoeff0[N-1]=1.0;
      //cout << "cmon=" << cmon[N-1] << "\n";
      deflcoeff0[N-2] = cmon[N-1]+r0;
      //cout << "deflcoeff[" << N-1 << "]=" << deflcoeff[N-1]<< "\n";
      for (i=N-3; i >=0; i--)
        {
          //cout << "cmon[" << i << "]=" << cmon[i] << " deflcoeff[" << i+1 << "]=" << deflcoeff[i+1]<< "\n";
          deflcoeff0[i] = cmon[i+1] + r0*deflcoeff0[i+1];
        }
      if (r0!=0)
        {
          deflcoeff1[0]=-cmon[0]/r0;
          deflcoeff1[N-1] = 1.0;
          //cout << "cmon[0]=" << cmon[0] << " deflcoeff1[0]=" << deflcoeff1[0]<< "\n";
          for (i=1; i < N-1; i++)
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
        }
      else
        deflcoeff = deflcoeff0;
      // refine deflated quartic coefficients here
      //deflcoeff.show("prima");
#if 1
       use_NR_defl_lin(deflcoeff,r0);
#else
       use_NR(deflcoeff,r0);
#endif
        //deflcoeff.show("dopo");
      // calculate roots of quartic
      quar.set_coeff(deflcoeff);
      quar.find_roots(qroots);
      //quar.show("deflq=");
      for (i=1; i < N; i++)
        {
          roots[i] = qroots[i-1];
        }
    }
  pvector<ntype,N+1> get_coeff()
    {
      return coeff;
    }
  int degree()
    {
      return n; 
    }
  inline void find_roots(pvector<complex<ntype>,N>& roots)
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
              quintic(roots);     
            }   
        }
      else if constexpr (N==6)
        {
          if constexpr (force_hqr)
            {
              int ok=1;
	      solve_numrec(roots, ok);
            }
          else
            {
              deg6roots(roots);
            }
        }
      else
        {
          // use HQR from NUMERICAL RECIPE
          int ok=1;
	  solve_numrec(roots, ok);
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
  rpoly()
    {
      n=N;
      meps = epsilon();
      maxf= getmax();
      maxdigits = numeric_limits<ntype>::digits10;
      //force_hqr=0;
      maxf2 = pow(maxf,0.5)/10.0;
      maxf3 = pow(maxf,1.0/3.0)/10.0;
      scalfact = pow(maxf,1.0/4.0)/1.618034;
      cubic_rescal_fact = pow(maxf, 1.0/3.0)/1.618034;
      //std::cout << "max2= " << maxf2<< " max3=" << maxf3 << "\n";
      //printf("%.15G\n", maxf);
      //std::cout << "macheps= " << std::setprecision(35) << meps << "\n";
    }

  rpoly(pvector<ntype,N> v)
    {
      coeff=v;
      rpoly();
      //cout << "macheps= " << setprecision(35) << eps << "\n";
    }

  ~rpoly()
    {

    }
};
// quadratic equation
template<class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::solve_quadratic(pvector<complex<ntype>,N>&sol)
{
  complex<ntype> r[2];
  ntype a,b;
  a = coeff[1]/coeff[2];
  b = coeff[0]/coeff[2];
  oqs_solve_quadratic(a, b, r);
  sol[0] = r[0];
  sol[1] = r[1];
}
// cubic polynomial
template <class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::solve_cubic_analytic(pvector<complex<ntype>,N>& sol)
{
  /* solve the cubic coeff[3]*x^3 + coeff[2]*x^2 +  coeff[1]*x + coeff[0] = 0
   * according to the method described in Numerical Recipe book */  
  ntype a, b, c, Q, R, theta, Q3, R2, A, B;
  const ntype sqrt32=sqrt(3)/2.0;
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
      A = -copysign(1.0,R)*pow(abs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = (A+B) - a/3.0;
      sol[1] = -0.5*(A+B)-a/3.0+complex<ntype>(0,1)*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-a/3.0-complex<ntype>(0,1)*sqrt32*(A-B);
    }
}

// quartics with OQS
template <class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol)
{
 /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
  * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const ntype PI2=M_PI/2.0, TWOPI=2.0*M_PI;
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
      KK = copysign(1.0,Q)*(RQ*RQ/Q-1.0);
    }

  if (KK < 0.0)
    {
      sqrtQ=sqrt(Q);
      theta = acos((R/abs(Q))/sqrtQ);
      if (theta < PI2) 
	sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
	sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (abs(Q) < abs(R))
	A = -copysign(1.0,R)*cbrt(abs(R)*(1.0+sqrt(KK)));
      else
	{
	  A = -copysign(1.0,R)*cbrt(abs(R)+sqrt(abs(Q))*abs(Q)*sqrt(KK));
	}
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol = A+B;
    }
}
template <class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, Q3, R2, A, B, sqrtQ;
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
      theta = acos(R/sqrt(Q3));
      sqrtQ=-2.0*sqrt(Q);
      if (theta < M_PI/2) 
	sol = sqrtQ*cos(theta/3.0);
      else 
	sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign(1.0,R)*pow(abs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol = A+B; /* this is always largest root even if A=B */
    }
}
template <class ntype, int N, bool force_hqr> void  rpoly<ntype,N,force_hqr>::oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (64) (see also the discussion in sec. 2.2 of the manuscript) */
  ntype rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  ntype maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  double ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
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
template <class ntype, int N, bool force_hqr> ntype  rpoly<ntype,N,force_hqr>::oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3)
{
  /* Eq. (21) in the manuscript */
  ntype sum;
  sum =  (b==0)?abs(d2 + l1*l1 + 2.0*l3):abs(((d2 + l1*l1 + 2.0*l3)-b)/b);
  sum += (c==0)?abs(2.0*d2*l2 + 2.0*l1*l3):abs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  sum += (d==0)?abs(d2*l2*l2 + l3*l3):abs(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
}
template <class ntype, int N, bool force_hqr> 
ntype rpoly<ntype,N, force_hqr>::oqs_calc_err_abcd_cmplx(ntype a, ntype b,  ntype c, ntype d, complex<ntype> aq, 
                                                         complex<ntype> bq, complex<ntype> cq, complex<ntype> dq)
{
  /* Eq. (53) in the manuscript for complex alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq) */
  ntype sum;
  sum = (d==0)?abs(bq*dq):abs((bq*dq-d)/d);
  sum += (c==0)?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?abs(aq + cq):abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, bool force_hqr> ntype rpoly<ntype,N,force_hqr>::oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq)
{
  /* Eq. (53) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  ntype sum;
  sum = (d==0)?abs(bq*dq):abs((bq*dq-d)/d);
  sum += (c==0)?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?abs(aq + cq):abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, bool force_hqr> ntype  rpoly<ntype,N,force_hqr>::oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq)
{
  /* Eq. (40) in the manuscript */
  ntype sum;
  sum = (c==0)?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?abs(aq + cq):abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ)
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
      errf += (vr[k1]==0)?abs(fvec[k1]):abs(fvec[k1]/vr[k1]);
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
	  errf += (vr[k1]==0)?abs(fvec[k1]):abs(fvec[k1]/vr[k1]);
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
template <class ntype, int N, bool force_hqr> void  rpoly<ntype,N, force_hqr>::oqs_solve_quadratic(ntype a, ntype b, complex<ntype> roots[2])
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
      roots[0]=complex<ntype>(zmax,0.0);
      roots[1]=complex<ntype>(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=complex<ntype>(-a/2,sqrtd/2);
      roots[1]=complex<ntype>(-a/2,-sqrtd/2);      
    }   
}

template <class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::oqs_quartic_solver(pvector<complex<ntype>,N>& roots)
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
  complex<ntype> acx1, bcx1, ccx1, dcx1,acx,bcx,ccx,dcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
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
      acx=complex<ntype>(l1,gamma);  
      bcx=complex<ntype>(l3,gamma*l2);
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
	  bcx1 = complex<ntype>(l3,sqrt(d3));
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
	  cdiskr=acx*acx/complex<ntype>(4.0)-bcx;               
	  /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
	  zx1=-acx/complex<ntype>(2.0)+sqrt(cdiskr);
	  zx2=-acx/complex<ntype>(2.0)-sqrt(cdiskr);
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
	  cdiskr=sqrt(acx*acx-complex<ntype>(4.0)*bcx);
	  zx1 = -complex<ntype>(0.5)*(acx+cdiskr);
	  zx2 = -complex<ntype>(0.5)*(acx-cdiskr);
	  if (abs(zx1) > abs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = bcx/zxmax;
	  roots[0] = zxmax;
	  roots[1] = zxmin;
	  cdiskr=sqrt(ccx*ccx-complex<ntype>(4.0)*dcx);
	  zx1 = -complex<ntype>(0.5)*(ccx+cdiskr);
	  zx2 = -complex<ntype>(0.5)*(ccx-cdiskr);
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
template<class ntype, int N, bool force_hqr> void rpoly<ntype,N,force_hqr>::balance(ntype a[N][N])
{
  const ntype RADIX=numeric_limits<ntype>::radix;
  int i, j;
  ntype scale[N];
  int done=0;
  ntype r, c, g, f, s, sqrdx=RADIX*RADIX;
  const int n=N;
  for (i=0; i < N; i++)
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
}

template<class ntype, int N, bool force_hqr> void rpoly<ntype,N, force_hqr>::hqr(ntype a[N][N], pvector<complex<ntype>,N>& wri, int& ok)
{
  int nn,m,l,k,j,its,i,mmin;
  ntype z,y,x,w,v,u,t,s,r=0.0,q=0.0,p=0.0, anorm=0.0;
  //const double EPS=2.2204460492503131E-16;
  const int n=N;
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
		      z=p+SIGN(z,p);
		      wri[nn-1]=wri[nn]=x+z;
		      if (z != 0.0)
			wri[nn]=x-w/z;
		    } 
		  else
		    {
		      //...a complex pair.
		      wri[nn]=complex<ntype>(x+p,-z);
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
			  if ((x=abs(p)+abs(q)+abs(r)) != 0.0)
			    {
			      p /= x;
			      //Scale to prevent over flow or under flow.
			      q /= x;
			      r /= x;
			    }
			}
		      if ((s=SIGN(sqrt(p*p+q*q+ r*r),p)) != 0.0)
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
template<class ntype, int N, bool force_hqr> void rpoly<ntype,N, force_hqr>::QRfactorization(ntype hess[N][N], pvector<complex<ntype>,N>& sol, int& ok)
{
  /* pag. 615 Num. Rec. */  
  balance(hess);
  hqr(hess, sol, ok);
}
template<class ntype, int N, bool force_hqr> void rpoly<ntype,N, force_hqr>::solve_numrec(pvector<complex<ntype>,N>& csol, int& ok)
{
  /* Find all the roots of a polynomial with real coefficients, 
   * coeff[4]*x^4+coeff[3]*x^3+coeff[2]*x^2+coeff[1]*x+coeff[0], 
   * The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pag. 497 Num. Rec. */
  const int m=N;
  ntype hess[N][N];
  int j, k;
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorization(hess, csol, ok);
}

