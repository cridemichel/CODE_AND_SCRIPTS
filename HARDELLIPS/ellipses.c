#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#define Sqr(x) ((x)*(x))
#ifndef CMPLX
#define CMPLX(x,y) (x)+(y)*I
#endif
double R[2][2][2];
double rx[2], ry[2], rz[2], sax[2][2];
/* perram wertheim overlap ellissoidi */
void tRDiagRpw2d(int i, double M[2][2], double D[2], double Ri[2][2])
{
  int k1, k2, k3;
  double Di[2][2];
  double Rtmp[2][2];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = D[0];
  Di[1][1] = D[1];
  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 2; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 2; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
void xlambda2d(double lambda, double rA[2], double A[2][2], double rB[2], double B[2][2], double x[2])
{
  double lamA[2][2], onemlamB[2][2], ABL[2][2], invABL[2][2];
  double x1[2], x2[2], x3[2], detinvABL;
  int k1, k2;
  /* calcola xlambda, vedi L. Paramonov and S. N. Yaliraki J. Chem. Phys. 123, 194111 (2005) */
  for (k1=0; k1 < 2; k1++)
    {
      for (k2=0; k2 < 2; k2++)
	{
	  lamA[k1][k2] = lambda*A[k1][k2];
	  onemlamB[k1][k2] = (1.0-lambda)*B[k1][k2];
 	  ABL[k1][k2] = lamA[k1][k2] + onemlamB[k1][k2];
	}
    }
  for (k1=0; k1 < 2; k1++)
    {
      x1[k1]=0;
      x2[k1]=0;
      for (k2 = 0; k2 < 2; k2++)
	{
	  x1[k1] += lamA[k1][k2]*rA[k2];
	  x2[k1] += onemlamB[k1][k2]*rB[k2];
	}
      x3[k1] = x1[k1] + x2[k1];
    }
  detinvABL=-Sqr(ABL[0][1]) + ABL[0][0]*ABL[1][1];
  //{{m11, -m01}, {-m01, m00}}
  invABL[0][0] = ABL[1][1];
  invABL[0][1] = -ABL[0][1];
  invABL[1][0] = -ABL[0][1];
  invABL[1][1] = ABL[0][0];

  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      invABL[k1][k2] /= detinvABL;

  for (k1 = 0; k1 < 2; k1++)
    {
      x[k1] = 0.0;
      for (k2 = 0; k2 < 2; k2++)
	{
	  x[k1] += invABL[k1][k2]*x3[k2];
	}
    }
}

double Slam2d(double lambda, double rA[3], double A[2][2], double rB[2], double B[2][2])
{
  int k1, k2;
  double xlam[2], fA[2], fB[2], SA, SB;

  xlambda2d(lambda, rA, A, rB, B, xlam);

  for (k1=0; k1 < 2; k1++)
    {
      fA[k1] = 0;
      fB[k1] = 0;
      for (k2=0; k2 < 2; k2++)
	{
	  fA[k1] += A[k1][k2]*(xlam[k2]-rA[k2]);
	  fB[k1] += B[k1][k2]*(xlam[k2]-rB[k2]);
	}
    }

  SA = SB = 0.0;
  for (k1=0; k1 < 2; k1++)
    {
      SA += lambda*(xlam[k1]-rA[k1])*fA[k1];
      SB += (1.0-lambda)*(xlam[k1]-rB[k1])*fB[k1];
    }
  /* ho messo un - così la funzione ha un minimo invece
     che un massimo e questo minimo viene trovato dalla funzione brentPW */
  return -(SA+SB);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
int brentPWTooManyIter=0;
double brentPW2d(double ax, double bx, double cx, double tol, double *xmin, double rA[2], double A[2][2], double rB[2], double B[2][2])
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
 * (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)),
 * this routine isolates the minimum to a fractional precision of about tol using Brent's
 * method. The abscissa of the minimum is returned as xmin, and the minimum function value 
 * is returned as brent, the returned function value. */
{ 
  int iter, ITMAXBR=100;
  const double CGOLD=0.3819660;
  const double ZEPSBR=1E-20;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0, fuold;
  brentPWTooManyIter=0;
  /* This will be the distance moved on the step before last.*/
  a=(ax < cx ? ax : cx); /*a and b must be in ascending order, 
			   but input abscissas need not be.*/
  b=(ax > cx ? ax : cx);
  x=w=v=bx; /*Initializations...*/
  fw=fv=fx=Slam2d(x, rA, A, rB, B); 
  if (fw < -1.0)
    {
      /* non-overlap! */
      *xmin=x;
      return -100.0;
    }
  fuold = fv;
  for (iter=1;iter<=ITMAXBR;iter++)
    { 
      /*Main program loop.*/
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPSBR); 
      if (fabs(x-xm) <= (tol2-0.5*(b-a)))
	{ /*Test for done here.*/
	  *xmin=x;
	  return fx;
	} 
      if (fabs(e) > tol1) 
	{ /*Construct a trial parabolic fit.*/
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if (q > 0.0)
	    p = -p; 
	  q=fabs(q);
	  etemp=e; 
	  e=d; 
	  if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=CGOLD*(e=(x >= xm ? a-x : b-x)); 
	    /*The above conditions determine the acceptability of the parabolic fit.
	     * Here we take the golden section step into the larger of the two segments.*/
	  else
	    {
	      d=p/q; /* Take the parabolic step.*/
	      u=x+d; 
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x); 
	    }
	}
      else
	{
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	} 
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=Slam2d(u, rA, A, rB, B); /*This is the one function evaluation per iteration.*/
      if (fu < -1.0)
	{
	  /* non overlap! */
	  *xmin=x;
	  return -100.0;
	}
#if 0
      if (2.0*fabs(fuold-fu) <= tol*(fabs(fuold)+fabs(fu)+ZEPSBR)) 
	{ 
	  *xmin=u;
	  return fu;
	}
#endif
      fuold = fu;//
      if (fu <= fx)
	{ /*Now decide what to do with our function evaluation.*/
	  if (u >= x) 
	    a=x;
	  else
	    b=x;
	  SHFT(v,w,x,u); /* Housekeeping follows:*/
	  SHFT(fv,fw,fx,fu); 
	} 
      else
	{ 
	  if (u < x) 
	    a=u; 
	  else 
	    b=u; 
	  if (fu <= fw || w == x)
	    {
	      v=w; w=u; fv=fw; fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w)
	    { 
	      v=u; fv=fu;
	    }
	} /* Done with housekeeping. Back for another iteration.*/
    }
  printf("Too many iterations in brent!\n");
  brentPWTooManyIter=1;
  //nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}
double calcfel2d(double M[2][2], double r0[2], double x[2]);

double check_overlap_pw2d(int i, int j, double shift[2])
{
  const double tolPW=DBL_EPSILON;
  double res, A[2][2], B[2][2], xmin; 
  int k1;
  double  DA[2], DB[2], rA[2], rB[2];
  //int typei, typej;

  //typei = typeOfPart[i];
  //typej = typeOfPart[j];

  rA[0] = rx[i];
  rA[1] = ry[i];

  rB[0] = rx[j]+shift[0];
  rB[1] = ry[j]+shift[1];

  for (k1=0; k1 < 2; k1++)
    {
      DA[k1]= 1.0/Sqr(sax[i][k1]);
      DB[k1]= 1.0/Sqr(sax[j][k1]);
    }
  tRDiagRpw2d(i, A, DA, R[i]);
  tRDiagRpw2d(j, B, DB, R[j]);
#if 0
  /* verifico che il centro di i non appartenga a j e viceversa come check preliminare */
  if (calcfel2d(A,rA,rB) < 0.0)
    return -1.0;
  if (calcfel2d(B,rB,rA) < 0.0)
    return -1.0;
#endif

  res = -brentPW2d(0, 0.5, 1.0, tolPW, &xmin, rA, A, rB, B);
  if (brentPWTooManyIter)
    {
      printf("res=%f xmin=%f\n", res, xmin);
      exit(-1);
    }
  //printf("res=%f\n", res);
  return res - 1.0;
}

const double cubic_rescal_fact = 3.488062113727083E+102; //= pow(DBL_MAX,1.0/3.0)/1.618034;
const double quart_rescal_fact = 7.156344627944542E+76; // = pow(DBL_MAX,1.0/4.0)/1.618034;
const double macheps = 2.2204460492503131E-16; // DBL_EPSILON
double oqs_max2(double a, double b)
{
  if (a >= b)
    return a;
  else
    return b;
}
double oqs_max3(double a, double b, double c)
{
  double t;
  t = oqs_max2(a,b);
  return oqs_max2(t,c);
}
#if 0
void solve_quadratic(double coeff[3], int *numsol, double *sol)
{
  double delta, a2inv, sqrtd;
  delta = Sqr(coeff[1]) - 4.0*coeff[2]*coeff[0];
  if (delta > 0.0)
    {
      sqrtd = sqrt(delta);
      a2inv = 1.0/(2.0*coeff[2]);
      sol[0] = (-coeff[1]+sqrtd)*a2inv;
      sol[1] = (-coeff[1]-sqrtd)*a2inv; 
      *numsol = 2;
    } 
  else if (delta == 0)
    {
      sol[0] = -coeff[1]/(2.0*coeff[2]);
      *numsol = 1;
    }
  else
    {
      *numsol = 0;
    }
}
#else
void solve_quadratic(double coeff[3], int *numsol, double *sol)
{
  /* numeric error safe version of solve_quadratic from Numerical Recipe */
  double delta, a, b, c, q;
  a = coeff[2];
  b = coeff[1];
  c = coeff[0];
  delta = Sqr(b) - 4.0*a*c;
  if (delta > 0.0)
    {
      q = -0.5*(b+copysign(1.0,b)*sqrt(delta));
      sol[0] = q/a;
      sol[1] = c/q;
      *numsol = 2;
    } 
  else if (delta == 0)
    {
      sol[0] = -b/(2.0*a);
      *numsol = 1;
    }
  else
    {
      *numsol = 0;
    }
}
#endif

void oqs_solve_cubic_analytic_depressed_handle_inf(double b, double c, double *sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  double Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const double PI2=M_PI/2.0, TWOPI=2.0*M_PI;
  Q = -b/3.0;
  R = 0.5*c;
  if (R==0)
    {
      if (b <= 0)
        {
          *sol=sqrt(-b);
        }
      else
        {
          *sol=0;
        }
      return;
    }

  if (fabs(Q) < fabs(R))
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
      theta = acos((R/fabs(Q))/sqrtQ);
      if (theta < PI2) 
        *sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
        *sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (fabs(Q) < fabs(R))
        A = -copysign(1.0,R)*cbrt(fabs(R)*(1.0+sqrt(KK)));
      else
        {
          A = -copysign(1.0,R)*cbrt(fabs(R)+sqrt(fabs(Q))*fabs(Q)*sqrt(KK));
        }
      if (A==0.0)
        B=0.0;
      else
        B = Q/A;
      *sol = A+B;
    }
}
void oqs_solve_cubic_analytic_depressed(double b, double c, double *sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  double Q, R, theta, Q3, R2, A, B, sqrtQ;
  Q = -b/3.0;
  R = 0.5*c;
  if (fabs(Q) > 1E102 || fabs(R) > 1E154)
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
        *sol = sqrtQ*cos(theta/3.0);
      else 
        *sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
        B=0.0;
      else
        B = Q/A;
      *sol = A+B; /* this is always largest root even if A=B */
    }
}
void oqs_calc_phi0(double a, double b, double c, double d, double *phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (79) (see also the discussion in sec. 2.2 of the manuscript) */
  double rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  double maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  double ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
  int iter;
  diskr=9*a*a-24*b;                    
  /* eq. (87) */
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
  /* eqs. (83) */
  aq=a+4*s;                                      
  bq=b+3*s*(a+2*s);                              
  cq=c+s*(2*b+s*(3*a+4*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9;
  hh=aq*cq;     

  g=hh-4*dq-3*gg;                       /* eq. (85) */  
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq; /* eq. (86) */          
  oqs_solve_cubic_analytic_depressed(g, h, &rmax);
  if (isnan(rmax) || isinf(rmax))
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
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
          oqs_solve_cubic_analytic_depressed(g, h, &rmax);
          if (isnan(rmax) || isinf(rmax))
            {
              oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
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
  if (fabs(xxx) > fabs(gx))
    maxtt = fabs(xxx);
  else
    maxtt = fabs(gx);
  if (fabs(h) > maxtt)
    maxtt = fabs(h);

  if (fabs(f) > macheps*maxtt)
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

          if (fabs(f) >= fabs(fold))
            {
              x = xold;
              break;
            }
        }
    }
  *phi0 = x;
}
double oqs_calc_err_ldlt(double b, double c, double d, double d2, double l1, double l2, double l3)
{
  /* Eqs. (29) and (30) in the manuscript */
  double sum;
  sum =  (b==0)?fabs(d2 + l1*l1 + 2.0*l3):fabs(((d2 + l1*l1 + 2.0*l3)-b)/b);
  sum += (c==0)?fabs(2.0*d2*l2 + 2.0*l1*l3):fabs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  sum += (d==0)?fabs(d2*l2*l2 + l3*l3):fabs(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
}
double oqs_calc_err_abcd_cmplx(double a, double b, double c, double d, 
                               complex double aq, complex double bq, complex double cq, complex double dq)
{
  /* Eqs. (68) and (69) in the manuscript for complex alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq) */
  double sum;
  sum = (d==0)?cabs(bq*dq):cabs((bq*dq-d)/d);
  sum += (c==0)?cabs(bq*cq + aq*dq):cabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?cabs(bq + aq*cq + dq):cabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?cabs(aq + cq):cabs(((aq + cq) - a)/a);
  return sum;
}
double oqs_calc_err_abcd(double a, double b, double c, double d, double aq, double bq, double cq, double dq)
{
  /* Eqs. (68) and (69) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  double sum;
  sum = (d==0)?fabs(bq*dq):fabs((bq*dq-d)/d);
  sum += (c==0)?fabs(bq*cq + aq*dq):fabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabs(bq + aq*cq + dq):fabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabs(aq + cq):fabs(((aq + cq) - a)/a);
  return sum;
}
double oqs_calc_err_abc(double a, double b, double c, double aq, double bq, double cq, double dq)
{
  /* Eqs. (48)-(51) in the manuscript */
  double sum;
  sum = (c==0)?fabs(bq*cq + aq*dq):fabs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabs(bq + aq*cq + dq):fabs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabs(aq + cq):fabs(((aq + cq) - a)/a);
  return sum;
}
void oqs_NRabcd(double a, double b, double c, double d, double *AQ, double *BQ, double *CQ, double *DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  double x02, errf, errfold, xold[4], x[4], dx[4], det, Jinv[4][4], fvec[4], vr[4];
  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
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
      errf += (vr[k1]==0)?fabs(fvec[k1]):fabs(fvec[k1]/vr[k1]);
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
          errf += (vr[k1]==0)?fabs(fvec[k1]):fabs(fvec[k1]/vr[k1]);
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
  *AQ=x[0];
  *BQ=x[1];
  *CQ=x[2];
  *DQ=x[3];
}
void oqs_solve_quadratic(double a, double b, complex double roots[2])
{ 
  double div,sqrtd,diskr,zmax,zmin;
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

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=CMPLX(-a/2,sqrtd/2);
      roots[1]=CMPLX(-a/2,-sqrtd/2);      
    }   
}
void oqs_quartic_solver(double coeff[5], complex double roots[4])      
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * 
   * if coeff[4] != 0 
   *
   * the four roots will be stored in the complex array roots[] 
   *
   * */
  complex double acx1, bcx1, ccx1, dcx1,acx,bcx,ccx,dcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
  double l2m[12], d2m[12], res[12], resmin, bl311, dml3l3, err0=0, err1=0, aq1, bq1, cq1, dq1; 
  double a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3],gamma,del2;
  int realcase[2], whichcase, k1, k, kmin, nsol;
  double rfactsq, rfact=1.0;

  if (coeff[4]==0.0)
    {
      printf("That's not a quartic!\n");
      return;
    }
  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  oqs_calc_phi0(a,b,c,d,&phi0, 0);

  // simple polynomial rescaling
  if (isnan(phi0)||isinf(phi0))
    {
      rfact = quart_rescal_fact; 
      a /= rfact;
      rfactsq = rfact*rfact;
      b /= rfactsq;
      c /= rfactsq*rfact;
      d /= rfactsq*rfactsq;
      oqs_calc_phi0(a,b,c,d,&phi0, 1);
    }
  l1=a/2;          /* eq. (16) */                                        
  l3=b/6+phi0/2;   /* eq. (18) */                                
  del2=c-a*l3;     /* defined just after eq. (27) */                             
  nsol=0;
  bl311 =2.*b/3.-phi0-l1*l1;   /* This is d2 as defined in eq. (20)*/ 
  dml3l3 = d-l3*l3;            /* dml3l3 is d3 as defined in eq. (15) with d2=0 */ 

  /* Three possible solutions for d2 and l2 (see eqs. (28) and discussion which follows) */
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
      /* Case I eqs. (37)-(40) */
      gamma=sqrt(-d2);                               
      aq=l1+gamma;                                  
      bq=l3+gamma*l2;                              

      cq=l1-gamma;                                
      dq=l3-gamma*l2;                            
      if(fabs(dq) < fabs(bq))
        dq=d/bq;                                
      else if(fabs(dq) > fabs(bq))
        bq=d/dq;                               
      if (fabs(aq) < fabs(cq))
        {
          nsol=0;
          if (dq !=0)
            {
              aqv[nsol] = (c - bq*cq)/dq;    /* see eqs. (47) */
              errv[nsol]=oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
              nsol++;
            }
          if (cq != 0) 
            {
              aqv[nsol] = (b - dq - bq)/cq;  /* see eqs. (47) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
              nsol++;
            }
          aqv[nsol] = a - cq;                /* see eqs. (47) */
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
              cqv[nsol] = (c - aq*dq)/bq;              /* see eqs. (53) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
              nsol++;
            }
          if (aq != 0)
            {
              cqv[nsol] = (b - bq - dq)/aq;            /* see eqs. (53) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
              nsol++;
            }
          cqv[nsol] = a - aq;                          /* see eqs. (53) */
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
      /* Case II eqs. (53)-(56) */
      gamma=sqrt(d2); 
      acx=CMPLX(l1,gamma);  
      bcx=CMPLX(l3,gamma*l2);
      ccx = conj(acx);
      dcx = conj(bcx);
      realcase[0] = 0; 
    }
  else 
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is better) */
  if (realcase[0]==-1 || (fabs(d2) <= macheps*oqs_max3(fabs(2.*b/3.), fabs(phi0), l1*l1))) 
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
          if(fabs(dq1) < fabs(bq1))  
            dq1=d/bq1;                                        
          else if(fabs(dq1) > fabs(bq1))
            bq1=d/dq1;                                       
          err1 = oqs_calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (68) */
        }
      else /* complex */
        {
          realcase[1] = 0;
          acx1 = l1;
          bcx1 = l3 + I*sqrt(d3);
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
      oqs_NRabcd(a,b,c,d,&aq,&bq,&cq,&dq);      
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
          cdiskr=acx*acx/4-bcx;               
          /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
          zx1=-acx/2+csqrt(cdiskr);
          zx2=-acx/2-csqrt(cdiskr);
          if(cabs(zx1) > cabs(zx2))
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
          cdiskr=csqrt(acx*acx-4.0*bcx);
          zx1 = -0.5*(acx+cdiskr);
          zx2 = -0.5*(acx-cdiskr);
          if (cabs(zx1) > cabs(zx2))
            zxmax = zx1;
          else
            zxmax = zx2;
          zxmin = bcx/zxmax;
          roots[0] = zxmax;
          roots[1] = zxmin;
          cdiskr=csqrt(ccx*ccx-4.0*dcx);
          zx1 = -0.5*(ccx+cdiskr);
          zx2 = -0.5*(ccx-cdiskr);
          if (cabs(zx1) > cabs(zx2))
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
        roots[k] *= rfact;
    }
}

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void balance(double a[6][6])
{
  const double RADIX=FLT_RADIX;
  int i, j;
  double scale[6]={1.0,1.0,1.0,1.0,1.0,1.0};
  int done=0;
  double r, c, g, f, s, sqrdx=RADIX*RADIX;
  const int n=6;
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
		c += fabs(a[j][i]);
		r += fabs(a[i][j]);
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

void hqr(double a[6][6], complex double wri[6], int *ok)
{
  int nn,m,l,k,j,its,i,mmin;
  double z,y,x,w,v,u,t,s,r=0.0,q=0.0,p=0.0, anorm=0.0;
  const double EPS=2.2204460492503131E-16;
  const int n=6;
  for (i=0;i<n;i++)
    //Compute matrix no rm for possible use in lo- cating single small sub diagonal element.
    for (j=IMAX(i-1,0);j<n;j++)
      anorm += fabs(a[i][j]);
  nn=n-1;
  t=0.0;
  *ok = 1;
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
	      s=fabs(a[l-1][l-1])+fabs(a[l][l]);
	      if (s == 0.0)
		s=anorm;
	      if (fabs(a[l][l-1]) + s == s)
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
		  z=sqrt(fabs(q));
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
		      wri[nn]=CMPLX(x+p,-z);
		      wri[nn-1]=conj(wri[nn]);
		    }
		  nn -= 2;
		} 
	      else
		{
		  //No roots found.  Continue iteration.
		  if (its == 480)
		    {
		      printf("Too many iterations in hqr");
		      *ok=0;
		      return;
		      //exit(-1);
		    }
		  if (its % 10 == 0 && its > 0)
		    {
		      //Form exceptional shift.
		      t += x;
		      for (i=0;i<nn+1;i++)
			a[i][i] -= x;
		      s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
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
		      s=fabs(p)+fabs(q)+fabs(r);
		      //Scale to prevent over flow or under flow.
		      p /= s;
		      q /= s;
		      r /= s;
		      if (m == l) 
			break;
		      u=fabs(a[m][m-1])*(fabs(q)+ fabs(r));
		      v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
		      if (u <= EPS*v)
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
			  if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0)
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
void QRfactorization( double hess[6][6], complex double sol[6], int *ok)
{
  /* pag. 615 Num. Rec. */  
  balance(hess);
  hqr(hess, sol, ok);
}
void solve_numrec (double coeff[7], complex double csol[6], int *ok)
{
  /* Find all the roots of a polynomial with real coefficients, 
   * coeff[4]*x^4+coeff[3]*x^3+coeff[2]*x^2+coeff[1]*x+coeff[0], 
   * The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pag. 497 Num. Rec. */
  const int m=6;
  double hess[6][6];
  int j, k;
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorization(hess, csol, ok);
}
void tRDiagRqe2d(double M[2][2], double D[2], double Ri[2][2])
{
  int k1, k2, k3;
  double Di[2][2];
  double Rtmp[2][2];
  Di[0][0] = D[0];
  Di[1][1] = D[1];
  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 2; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 2; k1++)
    for (k2 = 0; k2 < 2; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 2; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
double distSq2origM2d(double Alpha, double m00, double m01, double m11, double x0, double y0)
{
  double x[2], detMa;
  int i;
  detMa =-m01*m01 + (m00 + Alpha)*(m11 + Alpha);
  x[0] = x0*(-m01*m01 + m00*(m11 + Alpha)) + y0*(-m01*m11 + m01*(m11 + Alpha));
  x[1] = x0*(-m00*m01 + m01*(m00 + Alpha)) + y0*(-m01*m01 + m11*(m00 + Alpha));
  for (i=0; i< 2; i++)
    x[i] /= detMa;
 return x[0]*x[0]+x[1]*x[1]; 
}
double calcfel2d(double M[2][2], double r0[2], double x[2])
{
  int i, j;
  double res, v[2], xr0[2];
  
  for (i=0; i < 2; i++)
    xr0[i] = x[i] - r0[i];
  for (i=0; i < 2; i++)
    {
      v[i]=0;
      for (j=0; j < 2; j++)
        {
          v[i] += M[i][j]*xr0[j];
        }
    }
  res=0.0;
  for (i=0; i < 2; i++)
    {
      res += xr0[i]*v[i];
    }
  return res-1.0;
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}

double check_overlap_polyell_2D(int i, int j, double shift[3])
{
  const double GOLD=1.618034;
  double Rjp[2][2], r0jp[2], ri[2], rj[2], Di[2], Dj[2], Mjp[2][2], gradj[2], xppg[2], M2I[2][2], invM2I[2][2], oax[2], theta, norm;
  double RM[2][2], Mtmp[2][2], Mjpp[2][2], r0jpp[2], Mjp3[2][2], r0jp3[2], xp3g[2];
  double Mip[2][2], Mipp[2][2];
  double evec[2][2], eval[2], coeffpa[5], sx, sy, x0, y0, docc[2], doccn;
  double sx2, sy2, sx4, sy4, sz4, x02, y02, sai[2];
  double vt, at[2], Mi[2][2], Mj[2][2], Ri[2][2], Rj[2][2], dist, m00, m01, m11, b2d, m002, m012, m112;
  double alpha;
  int np, typei, typej;
  int kk1, kk2, kk3, k1, k2, k3, ok;
  complex double roots[4];
 
  //typei = typeOfPart[i];
  //typej = typeOfPart[j];
  
  /* apply affinity to reduce first ellipsoid to a sphere */
  for (k1=0; k1 < 2; k1++)
    {
      Di[k1]= 1.0/Sqr(sax[0][k1]);
      Dj[k1]= 1.0/Sqr(sax[1][k1]);
      for (k2=0; k2 < 2; k2++)
        {
          Ri[k1][k2] = R[i][k1][k2];
          Rj[k1][k2] = R[j][k1][k2];
        }
    }
  /* sai[0] e sai[1] sono i semiassi della particella i */
  for (k1=0; k1 < 2; k1++)
    {
      sai[k1] = sax[0][k1];
    } 
  tRDiagRqe2d(Mi, Di, Ri);
  tRDiagRqe2d(Mj, Dj, Rj);
  ri[0] = rx[i];
  ri[1] = ry[i];
  rj[0] = rx[j]+shift[0];
  rj[1] = ry[j]+shift[1];

  /* verifico che il centro di i non appartenga a j e viceversa come check preliminare */
  if (calcfel2d(Mi,ri,rj) < 0.0)
    return -1.0;
  if (calcfel2d(Mj,rj,ri) < 0.0)
    return -1.0;

  //printf("coords i: %f %f %f j: %f %f %f\n", ri[0], ri[1], ri[2], rj[0], rj[1], rj[2]);
  /* switch to ellipsoid i reference system */
  for (kk1=0; kk1 < 2; kk1++)
    {
      r0jp[kk1] = 0;
      for (kk2=0; kk2 < 2; kk2++)
        {
          r0jp[kk1] += R[i][kk1][kk2]*(rj[kk2]-ri[kk2]);
          Rjp[kk1][kk2] = 0;
          //Aip[kk1] = 0;
          for (kk3=0; kk3 < 2; kk3++)
            {
              Rjp[kk1][kk2] += R[j][kk1][kk3]*R[i][kk2][kk3];
              //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
            } 
        }
    }
  //tRDiagRpw(i, Mi, DA, R[i]);
  tRDiagRqe2d(Mjp, Dj, Rjp);

  /* calculate matrix and position of ellipsoid j after application of affinity
   * which reduces ellipsoid i to a sphere */   
  Mjpp[0][0] = Mjp[0][0]*Sqr(sai[0]);
  Mjpp[0][1] = Mjp[0][1]*sai[0]*sai[1];
  //Mjpp[1][0] = Mjp[1][0]*sai[0]*sai[1];
  Mjpp[1][1] = Mjp[1][1]*Sqr(sai[1]);
  r0jpp[0] = r0jp[0]/sai[0];
  r0jpp[1] = r0jp[1]/sai[1];
  x0 = r0jpp[0];
  y0 = r0jpp[1];
  x02 = Sqr(x0);
  y02 = Sqr(y0);
  m00 = Mjpp[0][0];
  m01 = Mjpp[0][1];
  m11 = Mjpp[1][1];
  m002 = Sqr(m00);
  m012 = Sqr(m01);
  m112 = Sqr(m11);
  //printf("M2d={{%.15G,%.15G},{%.15G,%.15G}}\n", m00, m01, m01, m11);
  //printf("xv02d={%.15G,%.15G}", x0, y0);
  b2d = -1.0 + m00*x02 + 2.0*m01*x0*y0 + m11*y02;

  //printf("b2d=%.15G\n", b2d);
  coeffpa[0] = -Sqr(m012 - m00*m11)*(-b2d + m00*x02 + y0*(2*m01*x0 + m11*y0)); 
  coeffpa[1] = -2*(m00 + m11)*(-m01*m01 + m00*m11)*(-b2d + m00*x02 + 
   y0*(2*m01*x0 + m11*y0));
  coeffpa[2] =b2d*(m002 - 2*m012 + 4*m00*m11 + m112) - (m002*m00 - 2*m00*m012 + 
    4*m002*m11 + m012*m11)*x02 - 
    2*m01*(m002 - 3*m012 + 5*m00*m11 + m112)*x0*y0 - (-2*m012*m11 + 
    m11*m112 + m00*(m012 + 4*m112))*y02;
  coeffpa[3] = 2*b2d*(m00 + m11) - 2*(m002 + m012)*x02 - 
    4*m01*(m00 + m11)*x0*y0 - 2*(m012 + m112)*y02;
  coeffpa[4] = b2d;
  oqs_quartic_solver(coeffpa, roots);
  for (kk1=0; kk1 < 4; kk1++)
    {
      //printf("root #%d=%.15G + I*(%.15G)\n", kk1, creal(roots[kk1]), cimag(roots[kk1]));
      if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
        {
          alpha=creal(roots[kk1]);
          break;
        }
    }
  dist=distSq2origM2d(alpha, m00, m01, m11, x0, y0) - 1.0;
  //printf("alpha=%.15G dist=%.15G\n", alpha, dist);
  /* trasformando tramite l'affinità inversa i punti che individuano la distanza tra sfera ed ellissoide
   * si avrà la distanza tra i due ellissoidi che si può usare nella dinamica event-driven */
  if (dist < 0.0)
    return -1.0;
  else
    return 1.0;
}

int main(int argc, char** argv)
{
  double L, theta, shift[2], over, over2;
  int tt, maxtrials, algo=1; // 0 = dryrun 1=polyell 2=PW 3=both of them to check consistency
  /* 
a1: 0.4823216982790199 b1:0.4306991794870728
a2: 0.4730194070139193 b2:0.3599600537746439
k1ax:-0.7586634226338149 k1ay:0.6514827788649102
k1bx:-0.6514827788649101 k1by:-0.7586634226338151
k2ax:-0.9931796094246184, k2ay:0.1165944399324530
k2bx:-0.1165944399324531 k2by:-0.9931796094246185
r1x:-16.8820392119991105 r1y:-19.3780665482635968
r2x:-16.2362572094466024 r2y:-18.8707434949179884
shiftx:0.0000000000000000 shifty:0.0000000000000000
a1,b1,a1,b2 = semisaani ellissi 1,2
k1ax, k1ay,k1bx,k1by= componenti versori dei semiassi ellisse 1
k2ax, k2ay,k2bx,k2by= componenti versori dei semiassi ellisse 2
r1x,r1y,r2x,r2y= componenti x,y centri ellissi 1,2
shiftx,shifty = shift centro ellisse 2

*/
  if (argc > 1)
    maxtrials=atoi(argv[1]);
  else
    maxtrials=10000;
  if (argc > 2)
    L = atof(argv[2]);
  else
    L = 6;
  if (argc > 2)
    algo=atoi(argv[3])>0?1:0;
  printf("L=%f\n", L);
#if 1 
  sax[0][0] = 1.0;
  sax[0][1] = 0.5;
  sax[1][0] = 1.0;
  sax[1][1] = 0.5;
  for (tt=0; tt < maxtrials; tt++)
    {
      rx[0] = (drand48()-0.5)*L;
      ry[0] = (drand48()-0.5)*L;
      rx[1] = (drand48()-0.5)*L;
      ry[1] = (drand48()-0.5)*L;
      theta = drand48()*2.0*M_PI;
      R[0][0][0] = cos(theta); // k1ax
      R[0][0][1] = sin(theta); // k1ay
      R[0][1][0] = -R[0][0][1]; // k1bx
      R[0][1][1] =  R[0][0][0]; // k1by
      theta = drand48()*2.0*M_PI;
      R[1][0][0] = cos(theta);// k2ax
      R[1][0][1] = sin(theta); // k2ay
      R[1][1][0] = -R[1][0][1]; // k2bx
      R[1][1][1] =  R[1][0][0]; // k2by
      shift[0] = 0.0;
      shift[1] = 0.0;
//#define POLYELL
     switch (algo) 
        {
        case 1:
          over= check_overlap_polyell_2D(0, 1, shift);
          break;
#if 0
          over2= check_overlap_pw2d(0, 1, shift);
          
#endif
        case 2:
          over = check_overlap_pw2d(0, 1, shift);
          break;
        case 3:
          over= check_overlap_polyell_2D(0, 1, shift);
          over2 = check_overlap_pw2d(0, 1, shift);
          if (over*over2 < 0.0)
            {
              printf("problem\n");
              exit(1);
            }
          break;
        default:
          break;
          // do nothin 
        }
    }
#else
  sax[0][0] = 0.4823216982790199;
  sax[0][1] = 0.4306991794870728;
  sax[1][0] = 0.4730194070139193;
  sax[1][1] = 0.3599600537746439;
  rx[0] = -16.8820392119991105;
  ry[0] = -19.3780665482635968;
  rx[1] = -16.2362572094466024;
  ry[1] = -18.8707434949179884;
  R[0][0][0] = -0.7586634226338149; // k1ax
  R[0][0][1] = 0.6514827788649102; // k1ay
  R[0][1][0] = -0.6514827788649101; // k1bx
  R[0][1][1] = -0.7586634226338151; // k1by
  R[1][0][0] = -0.9931796094246184; // k2ax
  R[1][0][1] = 0.1165944399324530; // k2ay
  R[1][1][0] = -0.1165944399324531; // k2bx
  R[1][1][1] = -0.9931796094246185; // k2by

  shift[0] = 0.0;
  shift[1] = 0.0;
  printf("overlap=%f\n", over);
#endif
  return 0;

}
