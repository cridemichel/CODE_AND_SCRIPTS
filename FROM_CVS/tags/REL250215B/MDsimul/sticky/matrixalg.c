#define TINY 1E-20
#define MD_NBMAX 8 
#include<mdsimul.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define MD_DEBUG(X) 
#define MD_DEBUG10(X) 
#define MD_DEBUG18(X) 
int ncom;
double (*nrfunc)(double []); 
int *ivector(int n)
{
  return calloc(n, sizeof(int)); 
}
double *vector(int n)
{
  return calloc(n, sizeof(double)); 
}
void polintRyck(double xain[], double yain[], int n, double x, double *y, double *dy);
void projectgrad(double *p, double *xi, double *gradf, double *gradg);
void projonto(double* ri, double *dr, double* rA, double **Xa, double *gradf, double *sfA, double dist);
void polint(double xain[], double yain[], int n, double x, double *y, double *dy);

double **matrix(int n, int m)
{
  double **M;
  int i;
  M = malloc(sizeof(double*)*n);
  for (i=0; i < n; i++)
    M[i] = calloc(m, sizeof(double));
  return M;
}
void free_vector(double *vec)
{
  free(vec);
}
void free_ivector(int *vec)
{
  free(vec);
}
void free_matrix(double **M, int n)
{
  int i;
  for (i=0; i < n; i++)
    {
      free(M[i]);
    }
  free(M);
}
extern void funcs2beZeroed(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroedDist(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroedDistNeg(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroedDistNeg5(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void fdjacDistNeg5(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3]);
extern void fdjacDistNeg(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3]);
extern void fdjacDist(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3]);

void calc_grad(double *rC, double *rA, double **Xa, double *grad);
double calc_norm(double *vec);
void (*nrdfun)(double [], double []);
/* =========================== >>> min <<< =============================== */
double min(double a, double b)
{
  if (a >= b)
    {
      return b;
    }
  else
    {
      return a;
    }
}
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
double sfA, sfB;
/* =========================== >>> max <<< ================================= */
double max(double a, double b)
{
  if (a >= b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


void nrerror(char *msg)
{
  printf(msg);
  exit(-1);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define DABS fabs
long long int itsfrprmn=0, callsfrprmn=0,callsok=0, callsprojonto=0, itsprojonto=0;
double accngA=0, accngB=0;
double xicom[8], pcomI[8], pcom[8], pcom2[8], xi[8], G[8], H[8], grad[8];//, vec[6];
double Ftol, Epoten, Emin, fnorm;
int cghalfspring, icg, jcg, doneryck;
double shiftcg[3], lambdacg, minaxicg, minaxjcg;
double gradfG[3], gradgG[3], dxG[6];
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, rA[3], rB[3], **RtA, **RtB;
extern int polinterr, polinterrRyck;
/* ============================ >>> brent <<< ============================ */
void  conjgradfunc(void);

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
/*Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the maximum
 * magnification allowed for a parabolic-fit step.*/
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
  /*Given a function func, and given distinct initial points ax and bx, this routine searches in the downhill
   * direction (defined by the function as evaluated at the initial points) and returns new points ax, bx, cx
   * that bracket a minimum of the function. Also returned are the function values at the three points, fa, fb,
   * and fc.*/
{
  const double GOLD= 1.618034, GLIMIT=100.0;
  double ulim,u,r,q,fu,dum; 
  *fa=(*func)(*ax); 
  *fb=(*func)(*bx);
  if (*fb > *fa) 
    { /*Switch roles of a and b so that we can go downhill in the direction from a to b. */
      SHFT(dum,*ax,*bx,dum); 
      SHFT(dum,*fb,*fa,dum); 
    } 
  *cx=(*bx)+GOLD*(*bx-*ax); /*First guess for c.*/
  *fc=(*func)(*cx);
  while (*fb > *fc) 
    { /*Keep returning here until we bracket.*/
      r=(*bx-*ax)*(*fb-*fc); /*Compute u by parabolic extrapolation from a, b, c. TINY is used to
			       prevent any possible division by zero.*/
      q=(*bx-*cx)*(*fb-*fa); 
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
	(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r)); 
      ulim=(*bx)+GLIMIT*(*cx-*bx); /*We won t go farther than this. Test various possibilities: */
      if ((*bx-u)*(u-*cx) > 0.0) 
	{ /*Parabolic u is between b and c: try it.*/
	  fu=(*func)(u);
	  if (fu < *fc) 
	    { /*Got a minimum between b and c. */
	      *ax=(*bx); *bx=u; *fa=(*fb); *fb=fu;
	      //printf("fa: %.15G fb: %.15G fc: %.15G\n", *fa, *fb, *fc);
	      return;
	    } 
	  else if (fu > *fb) 
	    { /*Got a minimum between between a and u.*/
	      *cx=u; *fc=fu;
	      //printf("fa: %.15G fb: %.15G fc: %.15G\n", *fa, *fb, *fc);
	      return;
	    } 
	  u=(*cx)+GOLD*(*cx-*bx);/* Parabolic fit was no use. Use default magnification.*/
	  fu=(*func)(u); 
	} 
      else if ((*cx-u)*(u-ulim) > 0.0) 
	{ /*Parabolic fit is between c and its allowed limit.*/
	  fu=(*func)(u); 
	  if (fu < *fc) 
	    {
	      SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
	      SHFT(*fb,*fc,fu,(*func)(u));
	    } 
	}
      else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	{ /*Limit parabolic u to maximum allowed value.*/
	  u=ulim; fu=(*func)(u); 
	} 
      else
	{ /*Reject parabolic u, use default magnification.*/
	  u=(*cx)+GOLD*(*cx-*bx); fu=(*func)(u); 
	} 
      SHFT(*ax,*bx,*cx,u);  /*Eliminate oldest point and continue.*/
      SHFT(*fa,*fb,*fc,fu);
    }
}
int powmeth;
double brentRyck(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
 * (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)),
 * this routine isolates the minimum to a fractional precision of about tol using Brent's
 * method. The abscissa of the minimum is returned as xmin, and the minimum function value 
 * is returned as brent, the returned function value. */
{ 
  int iter, ITMAXBR=100;
  const double CGOLD=0.3819660;
  const double ZEPSBR=1E-10;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0, fuold;
  /* This will be the distance moved on the step before last.*/
  a=(ax < cx ? ax : cx); /*a and b must be in ascending order, 
			   but input abscissas need not be.*/
  b=(ax > cx ? ax : cx);
  x=w=v=bx; /*Initializations...*/
  fw=fv=fx=(*f)(x); 
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
      fu=(*f)(u); /*This is the one function evaluation per iteration.*/
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
  polinterrRyck=1;
  nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
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
  /* This will be the distance moved on the step before last.*/
  a=(ax < cx ? ax : cx); /*a and b must be in ascending order, 
			   but input abscissas need not be.*/
  b=(ax > cx ? ax : cx);
  x=w=v=bx; /*Initializations...*/
  fw=fv=fx=(*f)(x); 
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
      fu=(*f)(u); /*This is the one function evaluation per iteration.*/
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
  nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}
#if 0
void linminConstr(double p[], double xi[], int n, double *fret, double (*func)(double []))
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
 * resets p to where the function func(p) takes on a minimum along the direction xi from p,
 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
 * the value of func at the returned location p. This is actually all accomplished by calling
 * the routines mnbrak and brent. */
{ 
  const double TOLLM=OprogStatus.tolSD/10;
  double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
  double f1dimConstr(double x); 
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, 
	      double (*func)(double));
  int j; 
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n; /*Define the global variables.*/
  //pcom=vector(1,n);
  //xicom=vector(1,n); 
  nrfunc=func; 
  for (j=0;j<n;j++)
    { 
      pcom[j]=p[j];
    }
  for (j=0;j<n;j++)
    { 
      xicom[j]=xi[j];
    } 
  //projectgrad(p,xicom,gradfG,gradgG);
  ax=0.0; /*Initial guess for brackets.*/
  xx=1.0; 
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dimConstr); 
  *fret=brent(ax,xx,bx,f1dimConstr,TOLLM,&xmin);
  for (j=0;j<n;j++)
    { /*Construct the vector results to return. */
      xi[j] *= xmin;
    } 
  projectgrad(p,xi,gradfG,gradgG);
  for (j=0;j<n;j++)
    {  
      p[j] += xi[j]; 
    }
  //free_vector(xicom,1,n); free_vector(pcom,1,n);
}
#endif
#if 0
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
 * resets p to where the function func(p) takes on a minimum along the direction xi from p,
 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
 * the value of func at the returned location p. This is actually all accomplished by calling
 * the routines mnbrak and brent. */
{ 
  const double TOLLM=1.0E-10;
  double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
  double f1dim(double x); 
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, 
	      double (*func)(double));
  int j; 
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n; /*Define the global variables.*/
  //pcom=vector(1,n);
  //xicom=vector(1,n); 
  nrfunc=func; 
  for (j=0;j<n;j++)
    { 
      pcom[j]=p[j];
      xicom[j]=xi[j];
    } 
  ax=0.0; /*Initial guess for brackets.*/
  xx=1.0; 
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim); 
  //printf("ax=%.15G xx=%.15G bx=%.15G\n", ax, xx, bx);
  *fret=brent(ax,xx,bx,f1dim,TOLLM,&xmin);
  //printf("xmin: %.15G\n", xmin);
  //printf("xi=%.15G %.15G %.15G %.15G %.15G %.15G \n",
	// xi[0], xi[1], xi[2], xi[3], xi[4], xi[5]);
  //printf("p=%.15G %.15G %.15G %.15G %.15G %.15G \n",
	// p[0], p[1], p[2], p[3], p[4], p[5]);
  for (j=0;j<n;j++)
    { /*Construct the vector results to return. */
      xi[j] *= xmin;
      p[j] += xi[j]; 
    } 
  //free_vector(xicom,1,n); free_vector(pcom,1,n);
}
#endif
#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f);
double dbrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double), double tol, double *xmin) 
  /* Given a function f and its derivative function df, and given a bracketing triplet of abscissas ax, bx, cx
   * [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)], this routine isolates the
   * minimum to a fractional precision of about tol using a modi cation of Brent s method that uses derivatives.
   * The abscissa of the minimum is returned as xmin, and 
   the minimum function value is returned as dbrent, the returned function value. */
{ 
  const double ZEPSDBR = 1E-10;
  const int ITMAXDBR=100; 
  int iter,ok1,ok2; /*Will be used as  ags for whether proposed steps are acceptable or not.*/
  double a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0; double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm; 
  a=(ax < cx ? ax : cx); b=(ax > cx ? ax : cx); x=w=v=bx; fw=fv=fx=(*f)(x); dw=dv=dx=(*df)(x); 
  /*All our housekeeping chores are doubled by the necessity of moving derivative values around as well 
   * as function values. */
  for (iter=1;iter<=ITMAXDBR;iter++)
    { 
      printf("x=%.15G iter=%d\n", x, iter );
      xm=0.5*(a+b); 
      tol1=tol*fabs(x)+ZEPSDBR; tol2=2.0*tol1; 
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) 
	{ 
	  *xmin=x; return fx; 
	} 
      if (fabs(e) > tol1)
	{ 
	  d1=2.0*(b-a); /*Initialize these d s to an out-of-bracket value.*/
	  d2=d1; 
	  if (dw != dx)
	    d1=(w-x)*dx/(dx-dw); /*Secant method with one point.*/
	  if (dv != dx)
	    d2=(v-x)*dx/(dx-dv); /*And the other. Which of these two estimates of d shall we take?
				   We will insist that they be within the bracket, and on the side
				   pointed to by the derivative at x: */
	  u1=x+d1; u2=x+d2; 
	  ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0; 
	  ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0; 
	  olde=e; 
	  /*Movement on the step before last.*/
	  e=d; 
	  if (ok1 || ok2) 
	    {
	      /* both are acceptable, then take the smallest one. */
	      if (ok1 && ok2) 
		d=(fabs(d1) < fabs(d2) ? d1 : d2); 
	      else if (ok1) 
		d=d1; 
	      else d=d2; 
	      if (fabs(d) <= fabs(0.5*olde)) 
		{ 
		  u=x+d; 
		  if (u-a < tol2 || b-u < tol2)
		    d=SIGN(tol1,xm-x); 
		} else 
		{ /* Bisect, not golden section.*/
		  d=0.5*(e=(dx >= 0.0 ? a-x : b-x)); /*Decide which segment by the sign of the derivative.*/
		} 
	    } 
	  else
	    {
	      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	    } 
	} 
      else
	{ 
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	} 
      if (fabs(d) >= tol1) 
	{ 
	  u=x+d; fu=(*f)(u);
	} 
      else
	{ 
	  u=x+SIGN(tol1,d); 
	  fu=(*f)(u);
	  if (fu > fx) 
	    { /*If the minimum step in the downhill direction takes us uphill, then we are done.*/
	      *xmin=x;
	      return fx;
	    }
	}
      du=(*df)(u); /*Now all the housekeeping, sigh.*/
      if (fu <= fx) 
	{ 
	  if (u >= x) a=x;
	  else
	    b=x;
	  MOV3(v,fv,dv, w,fw,dw); 
	  MOV3(w,fw,dw, x,fx,dx); 
	  MOV3(x,fx,dx, u,fu,du);
	} 
      else
	{
	  if (u < x)
	    a=u; 
	  else
	    b=u; 
	  if (fu <= fw || w == x)
	    { 
	      MOV3(v,fv,dv, w,fw,dw); 
	      MOV3(w,fw,dw, u,fu,du); }
	  else if (fu < fv || v == x || v == w)
	    {
	      MOV3(v,fv,dv, u,fu,du); 
	    }
	}
    } 
  nrerror("Too many iterations in routine dbrent"); 
  return 0.0; /*Never get here.*/
}

#if 1
extern double *axa, *axb, *axc;

double f1dimPow(double x) 
  /*Must accompany linmin.*/
{
  int j; 
  double f, xt[6];
  for (j=0;j<ncom;j++) 
    xt[j]=pcom[j]+x*xicom[j]; 
  f=(*nrfunc)(xt); 
  return f;
}
#if 0
void linminPow(double p[], double xi[], int n, double *fret, double (*func)(double []))
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
 * resets p to where the function func(p) takes on a minimum along the direction xi from p,
 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
 * the value of func at the returned location p. This is actually all accomplished by calling
 * the routines mnbrak and brent. */
{ 
  const double TOLLM=OprogStatus.tolSD;
  double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
  double f1dimPow(double x); 
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, 
	      double (*func)(double));
  int j; 
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n; /*Define the global variables.*/
  nrfunc=func; 
  for (j=0;j<n;j++)
    { 
      pcom[j]=p[j];
      xicom[j]=xi[j];
    } 
  ax=0; /*Initial guess for brackets.*/
  xx=1.0; 
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dimPow); 
  *fret=brent(ax,xx,bx,f1dimPow,TOLLM,&xmin);
  for (j=0;j<n;j++)
    { /*Construct the vector results to return. */
      xi[j] *= xmin;
      p[j] += xi[j]; 
    } 
}
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, 
	    double (*func)(double []))
  /* Minimization of a function func of n variables. Input consists of an initial starting
   * point p[1..n]; an initial matrix xi[1..n][1..n], whose columns contain the initial set
   * of directions (usually the n unit vectors); and ftol, the fractional tolerance in the
   * function value such that failure to decrease by more than this amount on one iteration
   * signals doneness. On output, p is set to the best point found, xi is the then-current 
   * direction set, fret is the returned function value at p, and iter is the number of 
   * iterations taken. The routine linmin is used.*/
{
  int i,ibig,j; 
  double del,fp,fptt,t,pt[6],ptt[6],xit[6]; 
  const int ITMAXPOW=OprogStatus.maxitsSD;
  //pt=vector(1,n);
  //ptt=vector(1,n); xit=vector(1,n);
  *fret=(*func)(p);
  for (j=0;j<n;j++) 
    pt[j]=p[j]; /*Save the initial point.*/
  for (*iter=1;;++(*iter))
    { 
      fp=(*fret);
      ibig=0;
      del=0.0;
      /*Will be the biggest function decrease.*/
      for (i=0;i<n;i++)
	{ /*In each iteration, loop over all directions in the set.*/
	  for (j=0;j<n;j++)
	    xit[j]=xi[j][i]; /*Copy the direction,*/
	  fptt=(*fret); 
	  linminPow(p,xit,n,fret,func); /*minimize along it,*/
	  if (fptt-(*fret) > del)
	    { /*and record it if it is the largest decrease so far.*/
	      del=fptt-(*fret); 
	      ibig=i;
	    }
	} 
      if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY)
	{
	  //free_vector(xit,1,n); /*Termination criterion.*/
	  //free_vector(ptt,1,n); free_vector(pt,1,n);
	  //printf("powell iter=%d\n", *iter);
	  return;
	}

      if (*iter == ITMAXPOW) 
	{
	  printf("powell iter=%d\n", *iter);
	  return;
	  nrerror("powell exceeding maximum iterations."); 
	}
      for (j=0;j<n;j++) 
	{ /*Construct the extrapolated point and the average direction moved.
	    Save the old starting point.*/
	  ptt[j]=2.0*p[j]-pt[j];
	  xit[j]=p[j]-pt[j];
	  pt[j]=p[j];
	} 
      fptt=(*func)(ptt); /*Function value at extrapolated point.*/
      if (fptt < fp)
	{ 
	  t=2.0*(fp-2.0*(*fret)+fptt)*Sqr(fp-(*fret)-del)-del*Sqr(fp-fptt); 
	  if (t < 0.0)
	    {
	      linminPow(p,xit,n,fret,func);/* Move to the minimum of the new direction, 
					   and save the new direction.*/
	      for (j=0;j<n;j++)
		{
		  xi[j][ibig]=xi[j][n-1]; 
		  xi[j][n-1]=xit[j];
		}
	    }
	}
    }
  /*Back for another iteration.*/
}
#endif
#endif
void lab2body(int i, double x[], double xp[], double *rO, double **R)
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      xp[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  xp[k1] += R[k1][k2]*(x[k2]-rO[k2]);
       	} 
    }
}
void body2labR(int i, double xp[], double x[], double *rO, double **R)
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
    }
}

void body2lab(int i, double xp[], double x[], double *rO, double **R)
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
      x[k1] += rO[k1];
    }
}

void angs2coord(double angs[], double p[])
{
  double sin0, sin2;
  sin0 = sin(angs[0]);
  p[0] = axa[icg]*cos(angs[1])*sin0;
  p[1] = axb[icg]*sin(angs[1])*sin0;
  p[2] = axc[icg]*cos(angs[0]);
  
  sin2 = sin(angs[2]);
  p[3] = axa[jcg]*cos(angs[3])*sin2;
  p[4] = axb[jcg]*sin(angs[3])*sin2;
  p[5] = axc[jcg]*cos(angs[2]);
}
double funcPowell(double angs[])
{
  double vec[6], vecP[6], gx2[3], fx2[3];
  double A, B, S=1.0, F;
  int k1, k2;
  angs2coord(angs, vecP);
  body2lab(icg, vecP, vec, rA, RtA);
  body2lab(jcg, &vecP[3], &vec[3], rB, RtB);
  if (OprogStatus.forceguess)
    {
      for (k1 = 0; k1 < 3; k1++)
	{
	  gx2[k1] = 0;
	  for (k2 = 0; k2 < 3; k2++)
	    gx2[k1] += 2.0*Xb[k1][k2]*(vec[k2] - rB[k2]);
       	  fx2[k1] = 0;
	  for (k2 = 0; k2 < 3; k2++)
	    fx2[k1] += 2.0*Xa[k1][k2]*(vec[k2+3] - rA[k2]);
	}
      A = B = 0.0;
      for (k1 = 0; k1 < 3; k1++)
	{
	  A += (vec[k1+3]-rA[k1])*fx2[k1];
	  B += (vec[k1]-rB[k1])*gx2[k1];
	}
      A = 0.5*A - 1.0;
      B = 0.5*B - 1.0;
      if (A<0 && B<0)
	S = -1.0;
      else
	S = 1.0;
    }
  F = 0;
  for (k1 = 0; k1 < 3; k1++)
    F += S*Sqr(vec[k1+3]-vec[k1]);
  return F;
}
extern double pi, **powdirs;
double powdirsI[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
#if 0
void powellmethod(double *vec)
{
  double Fret, r1p[3], r2p[3], vecP[6], angs[4], sinth;
  int iter, k1, k2;
  //powdirs = matrix(4,4);
  powmeth = 1;
  for (k1=0; k1 < 4; k1++)
    for (k2=0; k2 < 4; k2++)
      powdirs[k1][k2] = powdirsI[k1][k2];
  /* trovo le coordinate nel riferimento del corpo rigido (assi principali) */
  lab2body(icg, vec, vecP, rA, RtA);
  lab2body(jcg, &vec[3], &vecP[3], rB, RtB);
  /* determino theta e phi per entrambi gli ellissoidi 
   * angs[] = (thetaA,phiA,thetaB,phiB)
   * 0 < theta < 2PI
   * 0 < phi < PI */
  angs[0] = acos(vecP[2]/axc[icg]);
  angs[1] = acos(vecP[0]/axa[icg]/sin(angs[0]));
  sinth = vecP[1]/axb[icg]/sin(angs[0]);
  if (sinth < 0)
    angs[1] = 2.0*pi-angs[1];
  angs[2] = acos(vecP[5]/axc[jcg]);
  angs[3] = acos(vecP[3]/axa[jcg]/sin(angs[2]));
  sinth = vecP[4]/axb[jcg]/sin(angs[2]);
  if (sinth < 0)
    angs[3] = 2.0*pi-angs[3];
  powell(angs, powdirs, 4, OprogStatus.tolSD, &iter, &Fret, funcPowell);

  /* trovo le coordinate cartesiane dei punti trovati */
  angs2coord(angs, vecP);

  /* torno nel riferimento del laboratorio */
  body2lab(icg, vecP, vec, rA, RtA);
  body2lab(jcg, &vecP[3], &vec[3], rB, RtB);
  //free_matrix(powdirs, 4);
  powmeth = 0;
}
#endif
double powdirsP[6][6]={{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},
    {0,0,0,0,1,0},{0,0,0,0,0,1}};
double  cgfunc(double *vec);
#if 0
void powellmethodPenalty(double *vec)
{
  int k1, k2, iter; 
  double Fret;
  for (k1=0; k1 < 6; k1++)
    for (k2=0; k2 < 6; k2++)
      powdirs[k1][k2] = powdirsP[k1][k2];
  powell(vec, powdirs, 6, OprogStatus.tolSD, &iter, &Fret, cgfunc);
}
#endif
int check_point(char* msg, double *p, double *rc, double **XX)
{
  int k1, k2;
  double Q;
  Q=0;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  Q += (p[k1]-rc[k1])*XX[k1][k2]*(p[k2]-rc[k2]);  
	}
    }
  Q -= 1.0;
  if (msg)
    printf("[%s] should be zero: %.15G\n",msg, Q);
  else
    printf("[check_point] should be zero: %.15G\n", Q);
  if (fabs(Q)>1E-8)
    return 0;
  else
    return 1;
}
#if 0
void projonto(double* ri, double *dr, double* rA, double **Xa, double *gradf, double *sfA, double dist)
{
  int kk, its, done=0, k1, k2, MAXITS=50;
  const double GOLD=1.618034;
  double r1[3], r1A[3], sf, s1, s2, sqrtDelta, A2;
  double A, B, C, Delta, sol=0.0, ng;
  sf = *sfA;
  its = 0;
 
  callsprojonto++;
  while (!done && its <= MAXITS)
    {
      itsprojonto++;
      //printf("sf*dr=%.15G %.15G %.15G dr=%.15G\n", dr[0], dr[1], dr[2], calc_norm(dr));
#if 0
      if (!check_point("inside loop", ri, rA, Xa))
	{
	  exit(-1);
	}
#endif
      for (kk=0; kk < 3; kk++)
	{
	  r1[kk] = ri[kk] + dr[kk]*sf; 
	  r1A[kk] = r1[kk] - rA[kk];
	}
      A=0;
      B=0;
      C=0;
      for (k1=0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      A += gradf[k1]*Xa[k1][k2]*gradf[k2];
	      B += r1A[k1]*Xa[k1][k2]*gradf[k2];
	      //  printf("riA[%d]=%f Xa[%d][%d]=%f\n", r1A[k1], Xa[k1][k2]);
	      C += r1A[k1]*Xa[k1][k2]*r1A[k2];
	    }
	}
      B *= 2.0;
      C -= 1.0;
      Delta = Sqr(B) - 4.0*A*C;
      if (Delta < 0 || A==0.0)
	{
	  sf /= GOLD;
	  its++;
	  continue;
	}
      sqrtDelta = sqrt(Delta);
      A2 = 2.0*A;
#if 1
      s1 = (-B - sqrtDelta)/A2;
      s2 = (-B + sqrtDelta)/A2;
      if (fabs(s1) < fabs(s2))
	sol = s1;
      else
	sol = s2;
#else
      if (B > 0)
	sol = (-B + sqrtDelta)/A2; 
      else
	sol = (-B - sqrtDelta)/A2;
#endif
      ng = calc_norm(gradf);
      //if (dist > OprogStatus.epsd && fabs(sol)*ng > OprogStatus.tolSDconstr*sf*calc_norm(dr))
      if (dist > OprogStatus.epsd && fabs(sol)*ng > OprogStatus.tolSDconstr*dist/2.0)
	{
	  sf /= GOLD;
	  its++;
	  continue;
	}
    
      done = 1;
    
    }
 
  if (!done)
    {
      printf("maximum number of iterations reached in projont! Aborting...\n");
      printf("sol=%.15G norm(dr)=%.15G sf=%.15G\n", sol, calc_norm(dr), sf);
      exit(-1);
    }
#if 0
    {
      double Q= 0;
      int k1, k2;
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      Q += (ri[k1]+sf*dr[k1]+OprogStatus.stepSD*sol*gradf[k1]-rA[k1])*Xa[k1][k2]*
		(ri[k2]+ sf*dr[k2]+OprogStatus.stepSD*sol*gradf[k2]-rA[k2]);  
	    }
	}
      Q -= 1.0;
      printf("should be zero: %.15G sf:%.15G\n", Q, sf);
    }
  printf("done! sol=%.15G\n", sol);
#endif
  for (kk = 0; kk < 3; kk++)
    {
      dr[kk] = sol*gradf[kk] + sf*dr[kk]; 
    }
  *sfA = sf;
}
#endif
#if 0
void projectgrad(double *p, double *xi, double *gradf, double *gradg)
{
  int kk, k1, k2, k3;
  double r1[3], r2[3], rIf[3], rIg[3], pp[3], dist;
  double r1A[3], r2B[3], A1, B1, C1, Delta1, A2, B2, C2, Delta2; 
  //double gradf[3], gradg[3];
#if 0
  for (kk=0; kk < 3; kk++)
    {
      r1[kk] = p[kk]+xi[kk];
      r2[kk] = p[kk+3]+xi[kk+3];
    }
  calc_intersec(r1, rA, Xa, rIf);
  calc_intersec(r2, rB, Xb, rIg);
  for (kk=0; kk < 3; kk++)
    {
      r1[kk] += rIf[kk] - r1[kk];
      r2[kk] += rIg[kk] - r2[kk];
      xi[kk] = r1[kk] - p[kk];
      xi[kk+3] = r2[kk] - p[kk+3];
    }
  calc_grad(p, rA, Xa, gradf);
  calc_grad(&p[3], rB, Xb, gradg);
#endif
  dist = 0;
  for (kk=0; kk < 3; kk++)
    {
      dist+=Sqr(p[kk+3]-p[kk]);
    }
  dist = sqrt(dist);
  projonto(p, xi, rA, Xa, gradf, &sfA, dist);
  projonto(&p[3], &xi[3], rB, Xb, gradg, &sfB, dist);
}
#endif
int maxitsRyck;

double xaRyck[3], yaRyck[3];
double polintfuncRyck(double x)
{
  double dy=0.0, y;
  polintRyck(xaRyck, yaRyck, 3, x, &y, &dy);
  if (polinterrRyck==1)
    return 0.0;
  if (dy > OprogStatus.epsd)
    {
      polinterrRyck = 1;
    }
  else 
    polinterrRyck = 0;
  return y;
}
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
double Asd, Bsd, OmegaSqAsd[3][3], OmegaAsd[3][3], OmegaBsd[3][3], OmegaSqBsd[3][3];
extern void calc_intersec(double *rB, double *rA, double **Xa, double* rI);
double calc_dist(double *p)
{
  int kk;
  double D;
  D = 0;
  for (kk=0; kk < 3; kk++)
    D += Sqr(p[kk+3]-p[kk]);
  D = sqrt(D);
  return D;
}
#if 0
void updateByRot(double p[], double xi[])
{
  int kk, k1, k2;
  double sinwA, sinwB, coswA, coswB, phiA, phiB;
  double MA[3][3], pn[6], MB[3][3], A, B, distini, distfine, scp;
  double omA[3], omB[3], vA[3], vB[3], r1[3], r2[3], nA, nB, pi[6];
  double F, S;
  for (kk =0; kk < 3; kk++)
    {
      r1[kk] = p[kk] - rA[kk];
      r2[kk] = p[kk+3] - rB[kk];
    }

  vectProdVec(r1, xi, omA);
  vectProdVec(r2, &xi[3], omB);
  nA = calc_norm(omA);
  nB = calc_norm(omB);
  for (kk=0; kk < 3; kk++)
    { 
      omA[kk] /= nA;
      omB[kk] /= nB;
    }
  OmegaAsd[0][0] = 0;
  OmegaAsd[0][1] = -omA[2];
  OmegaAsd[0][2] = omA[1];
  OmegaAsd[1][0] = omA[2];
  OmegaAsd[1][1] = 0;
  OmegaAsd[1][2] = -omA[0];
  OmegaAsd[2][0] = -omA[1];
  OmegaAsd[2][1] = omA[0];
  OmegaAsd[2][2] = 0;
  OmegaSqAsd[0][0] = -Sqr(omA[1]) - Sqr(omA[2]);
  OmegaSqAsd[0][1] = omA[0]*omA[1];
  OmegaSqAsd[0][2] = omA[0]*omA[2];
  OmegaSqAsd[1][0] = omA[0]*omA[1];
  OmegaSqAsd[1][1] = -Sqr(omA[0]) - Sqr(omA[2]);
  OmegaSqAsd[1][2] = omA[1]*omA[2];
  OmegaSqAsd[2][0] = omA[0]*omA[2];
  OmegaSqAsd[2][1] = omA[1]*omA[2];
  OmegaSqAsd[2][2] = -Sqr(omA[0]) - Sqr(omA[1]);
  OmegaBsd[0][0] = 0;
  OmegaBsd[0][1] = -omB[2];
  OmegaBsd[0][2] = omB[1];
  OmegaBsd[1][0] = omB[2];
  OmegaBsd[1][1] = 0;
  OmegaBsd[1][2] = -omB[0];
  OmegaBsd[2][0] = -omB[1];
  OmegaBsd[2][1] = omB[0];
  OmegaBsd[2][2] = 0;
  OmegaSqBsd[0][0] = -Sqr(omB[1]) - Sqr(omB[2]);
  OmegaSqBsd[0][1] = omB[0]*omB[1];
  OmegaSqBsd[0][2] = omB[0]*omB[2];
  OmegaSqBsd[1][0] = omB[0]*omB[1];
  OmegaSqBsd[1][1] = -Sqr(omB[0]) - Sqr(omB[2]);
  OmegaSqBsd[1][2] = omB[1]*omB[2];
  OmegaSqBsd[2][0] = omB[0]*omB[2];
  OmegaSqBsd[2][1] = omB[1]*omB[2];
  OmegaSqBsd[2][2] = -Sqr(omB[0]) - Sqr(omB[1]);
  phiA = sfA;
  phiB = sfB;
  sinwA = sin(phiA);
  coswA = (1.0 - cos(phiA));
  sinwB = sin(phiB);
  coswB = (1.0 - cos(phiB));

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  MA[k1][k2] = sinwA*OmegaAsd[k1][k2]+coswA*OmegaSqAsd[k1][k2];
	  MB[k1][k2] = sinwB*OmegaBsd[k1][k2]+coswB*OmegaSqBsd[k1][k2];
	}
    }

  for (k1 = 0; k1 < 3; k1++)
    {
      vA[k1] = p[k1]-rA[k1];
      vB[k1] = p[k1+3]-rB[k1];
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      pn[k1] = vA[k1];
      pn[k1+3] = vB[k1];

      for (k2 = 0; k2 < 3; k2++)
	{
	  pn[k1] += MA[k1][k2]*vA[k2];
	  pn[k1+3] += MB[k1][k2]*vB[k2]; 
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      pn[k1] += rA[k1];
      pn[k1+3] += rB[k1];
    }	
  calc_intersec(pn, rA, Xa, p);
  calc_intersec(&pn[3], rB, Xb, &p[3]);

  //printf("scal prod1=%.15G scalprod2=%.15G\n", scalProd(pn,omA), scalProd(&pn[3],omB));
  //printf("pcomI = (%.15G,%.15G,%.15G)\n", pcomI[0], pcomI[1], pcomI[2]);
  //distfine= calc_dist(p);
  //for (kk=0; kk < 6; kk++)
    //pi[kk] = p[kk]-pi[kk];
  //scp = scalProd(pi,xi);
  //printf("scp=%.15G\n", scp);
  //if (scp < 0)
    //printf("scp=%.15G distini=%.15G distfine=%.15G phiA=%.15G phiB=%.15G\n", scp, distini, distfine, phiA,phiB);
}
#endif
double  cgfuncRyck(double *vec);
double zbrentRyck(double (*func)(double), double x1, double x2, double tol);
double get_sign(double *vec);
#if 0
int check_done(double fp, double fpold, double minax)
{
  const double EPSFR=1E-10;
  if (OprogStatus.tolSDgrad > 0)
    {
      if (fp > Sqr(OprogStatus.epsd)) 
	{
	  if (doneryck == 1 || 
	      2.0*fabs(fpold-fp) < OprogStatus.tolSDlong*(fabs(fpold)+fabs(fp)+EPSFR))
	  //fabs(fpold-fp) < OprogStatus.tolSDlong*Sqr(minax*2))
	    return 1;
	}
      else 
	{
	  if (2.0*fabs(fpold-fp) < OprogStatus.tolSD*(fabs(fpold)+fabs(fp)+EPSFR))
	    return 1;
	}
    }
  else
    {
      if (fp < Sqr(OprogStatus.epsd))
	{
	  if (2.0*fabs(fpold-fp) < OprogStatus.tolSD*(fabs(fpold)+fabs(fp)+EPSFR))
	    return 1;
	}
      else
	{
	  if (2.0*fabs(fpold-fp) < OprogStatus.tolSDlong*(fabs(fpold)+fabs(fp)+EPSFR))
	    return 1;
	  //if (fabs(fpold-fp) < OprogStatus.tolSDlong*Sqr(minax*2))
	    //return 1;
	}
    }
  return 0;
}
#endif
#if 0
void frprmnRyck(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), double (*dfunc)(double [], double [], double [], double [], double*, double*))
  /*Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func,
   * using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is
   * input as ftol. Returned quantities are p (the location of the minimum), iter
   * (the number of iterations that were performed), and fret (the minimum value of the function).
   * The routine linmin is called to perform line minimizations. */
{ 
  int j,its,kk;
  const int ITMAXFR = OprogStatus.maxitsSD;
  const double EPSFR=1E-10, GOLD=1.618034;
  double normxi,gg,gam,fp,dgg,norm1,norm2, sp, fpold, gradf[3], gradg[3], signA, signB;
  double minax, distini, distfin, dist, g[6],h[6],xi[6], dx[3], fx[3], gx[3], dd[3], xiold[6];
  double pm[6], fpm, signAold, signBold, pold[6];
  double xmin, xim[6];
  //printf("primaprima p= %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2], p[3], p[4], p[5]);
 
  minax = min(minaxicg,minaxjcg);
  sfA = icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
  sfB = jcg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
  callsfrprmn++;
  /*Initializations.*/
  fp = (*dfunc)(p,xi,gradfG,gradgG, &signA, &signB); 
  if (doneryck==2)
    {
      callsok++;
      return;
    }
  projectgrad(p,xi,gradfG,gradgG);  
  for (its=1;its<=ITMAXFR;its++)
    { 
      itsfrprmn++;      
      *iter=its;
#if 0
      linminRyck(p, xi, &fp);  
      //updateByRot(p, xi);
      //linminConstr(p, xi, 6, &fret, cgfuncRyck);
#else
      for (j=0; j < n; j++)
	{
	  pold[j] = p[j];
	  xiold[j] = xi[j];
	  p[j] += xi[j];
	}
#endif
      signAold = signA;
      signBold = signB;
      fpold = fp; 
      fp = (*dfunc)(p,xi,gradfG, gradgG, &signA, &signB);
      if (fp > fpold)
	{
	  sfA /= GOLD;
	  sfB /= GOLD;
	}      
      projectgrad(p, xi, gradfG, gradgG);
      if (doneryck==2)
	{
	  callsok++;
	  return;
	 }
      if (check_done(fp, fpold, minax))
	{
	  callsok++;
	  return;
	}
#if 0
	if ((OprogStatus.tolSDgrad <=0|| 
	    doneryck==1) 
	   && (OprogStatus.tolSD <=0 || 
	       2.0*fabs(fpold-fp) < ftol*(fabs(fpold)+fabs(fp)+EPSFR)))
	 {
#if 0
	   double ngA, ngB;
	   (*dfunc)(p,xi,gradfG, gradgG, &signA, &signB);
	   ngA = ngB = 0;  
	   for (kk=0; kk < 3; kk++)
	     {
	       ngA += Sqr(xi[kk]);
	       ngB += Sqr(xi[kk+3]); 
	     }
	   ngA = sqrt(ngA);
	   ngB = sqrt(ngB);
	   accngA += ngA/(icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB);
	   accngB += ngB/(jcg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB);
#endif
	   callsok++;
	   return;
	 }
#endif
    } 
  return; 
  nrerror("Too many iterations in frprmn");
  
}
#endif
double min3(double a, double b, double c)
{
  double m;
  m = a;
  if (b < m)
    m = b;
  if (c < m)
    m = c;
  return m;
}
void estimate_lambda(double *vec)
{
  int kk, k1, k2;
  double gx[3], fx[3], grad[6];
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(vec[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(vec[k2+3] - rB[k2]);
    }
  for (kk=0; kk < 3; kk++)
    {
      grad[kk]= -2.0*(vec[kk+3]-vec[kk]);
      grad[kk+3]= 2.0*(vec[kk+3]-vec[kk]);
    }
  vec[6] = 0;
  for (k1=0 ; k1 < 3; k1++)
    vec[6] += fx[k1]*grad[k1];
  vec[7] = 0;
  for (k1=0 ; k1 < 3; k1++)
    vec[7] += gx[k1]*grad[k1+3];
  vec[6] = fabs(vec[6]);
  vec[7] = fabs(vec[7]);
}
double get_sign(double *vec)
{
  int k1, k2;
  double A, B, S, fx2[3], gx2[3];
  for (k1 = 0; k1 < 3; k1++)
    {
      gx2[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx2[k1] += 2.0*Xb[k1][k2]*(vec[k2] - rB[k2]);
      fx2[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx2[k1] += 2.0*Xa[k1][k2]*(vec[k2+3] - rA[k2]);
	}
  A = B = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      A += (vec[k1+3]-rA[k1])*fx2[k1];
      B += (vec[k1]-rB[k1])*gx2[k1];
    }
  A = 0.5*A - 1.0;
  B = 0.5*B - 1.0;
  
  if (A<0 && B<0)
    S = -1.0;
  else
    S = 1.0;
  
  return S;
}
typedef struct {
	double point[3];
//	double grad[3];
	struct {
	  int i;	  
	  int j;
	} neigh[4];
} MESHXYZ;

MESHXYZ **ellips_mesh[2];

int choose_neighbour2(double *grad, int *th1, int *phi1, int *th2, int *phi2,
		     double *distSq, double *S, double maxstA, double maxstB, double *vec, 
		     int calc_sign)
{
  int k1, k2, kk, cth1, cphi1, cth2, cphi2, mA, mB;
  int nphi1, nth1, nphi2, nth2, firstdist=1;
  double distSqold, distSqT, distSqmin=0.0;
  double vecini[6], cxmesh[6],xmesh[6], xmeshP1[3], xmeshP2[3];

  cth1 = *th1;
  cphi1 = *phi1;
  cth2 = *th2;
  cphi2 = *phi2;
  mA = (icg<Oparams.parnumA)?0:1;
  mB = (jcg<Oparams.parnumA)?0:1;
  //printf("maxstA=%.15G maxstB=%.15G\n", maxstA, maxstB);
  for (kk=0; kk < 6; kk++)
    vecini[kk] = vec[kk];
#if 0
  printf("$$$PR VECA=(%.15G,%.15G, %.15G)\n", vec[0], vec[1], vec[2]);
  printf("$$$PR VECB=(%.15G,%.15G, %.15G)\n", vec[3], vec[4], vec[5]);
#endif
  for (k1 = 0; k1 < 4; k1++)
    {
      for (k2=0; k2 < 4; k2++)
	{
	  nth1  = ellips_mesh[mA][*th1][*phi1].neigh[k1].i;
	  nphi1 = ellips_mesh[mA][*th1][*phi1].neigh[k1].j;
	  nth2  = ellips_mesh[mB][*th2][*phi2].neigh[k2].i;
	  nphi2 = ellips_mesh[mB][*th2][*phi2].neigh[k2].j;
      	  if (nth1==-1 || nphi1==-1|| nth2 ==-1 || nphi2 ==-1)
	    continue;
	  for (kk=0; kk < 3; kk++) 
	    {
	      xmeshP1[kk] = ellips_mesh[mA][nth1][nphi1].point[kk];
	      //dxP[kk] = xmeshP1[kk] - ellips_mesh[mA][*th1][*phi1].point[kk]; 
	      xmeshP2[kk] = ellips_mesh[mB][nth2][nphi2].point[kk];
	      //dxP[kk] = xmeshP[kk] - ellips_mesh[mB][*th2][*phi2].point[kk]; 
	    }
	  body2lab(icg, xmeshP1, xmesh, rA, RtA);
	  body2lab(jcg, xmeshP2, &xmesh[3], rB, RtB);
	  //normdx = calc_norm(dx);
	  distSqT = 0;
	  for (kk=0; kk < 3; kk++) 
	    distSqT += Sqr(xmesh[kk+3]-xmesh[kk]);
	  if (calc_sign)
	    {
	      *S = get_sign(xmesh);
	      distSqT *= *S;
	    }
	  distSqT *= *S;
	  if (firstdist || distSqT < distSqmin )
	    {
	      firstdist=0;
	      distSqmin = distSqT;
	      cth1 = nth1;
	      cphi1 = nphi1;
	      cth2= nth2;
	      cphi2 = nphi2;
	      for (kk=0; kk < 6; kk++)
		{
		  cxmesh[kk] = xmesh[kk];
		}
	    }
	}
    }
  /* calcola la nuova distanza, il nuovo gradiente e le nuove coordinate del punto */
  //printf("nth2 nphi2 = %d %d *th2 *nphi2 = %d %d\n", nth2, nphi2, *th2, *phi2);
  //printf("nth1 nphi1 = %d %d *th1 *nphi1 = %d %d\n", nth1, nphi1, *th1, *phi1);
  for (kk = 0; kk < 6; kk++)
    {
      vec[kk] =  cxmesh[kk];
    }
#if 0
    {
      double VP[6], VL[6];
      for (kk=0; kk < 3; kk++)
	{
	  VP[kk] = ellips_mesh[mA][*th1][*phi1].point[kk]; 
	  VP[kk+3]=ellips_mesh[mB][*th2][*phi2].point[kk];  
	}
      body2lab(icg, VP, VL, rA, RtA);
      body2lab(jcg, &VP[3], &VL[3], rB, RtB);
      printf("$$$ VA= (%.15G,%.15G, %.15G)\n", VL[0], VL[1], VL[2]);
      printf("$$$ VB= (%.15G,%.15G, %.15G)\n", VL[3], VL[4], VL[5]);
      printf("$$$ VECA=(%.15G,%.15G, %.15G)\n", vec[0], vec[1], vec[2]);
      printf("$$$ VECB=(%.15G,%.15G, %.15G)\n", vec[3], vec[4], vec[5]);
    }
#endif
  distSqold = *distSq;
  *distSq = distSqmin;
  //printf("distSqold=%.15G distSq=%.15G\n",sqrt(fabs(distSqold)), sqrt(fabs(*distSq)));
  if (*distSq > distSqold)
    {
      for (kk=0; kk < 6; kk++)
	vec[kk] = vecini[kk];
      *distSq = distSqold;
      return 1;
    }
  *th1 = cth1;
  *phi1 = cphi1;
  *th2 = cth2;
  *phi2 = cphi2;
  //printf("USCENDO *th1 *phi1 = %d %d *th2 *phi2 = %d %d\n", *th1, *phi1, *th2, *phi2);
  return 0;	
}

int choose_neighbour(double *grad, int *th1, int *phi1, int *th2, int *phi2,
		     double *distSq, double *S, double maxstA, double maxstB, double *vec, 
		     int calc_sign)
{
  int k1, kk, cth1, cphi1, cth2, cphi2, mA, mB;
  int nphi1, nth1, nphi2, nth2, firstspA=1, firstspB=1;
  double sp, spmaxA=0.0, spmaxB=0.0, dx[3], dxP[3], cdxA[3], cdxB[3], distSqold;
  double vecini[6], xmesh[3], xmeshP[3], xmeshA[3], xmeshB[3];

  cth1 = *th1;
  cphi1 = *phi1;
  cth2 = *th2;
  cphi2 = *phi2;
  mA = (icg<Oparams.parnumA)?0:1;
  mB = (jcg<Oparams.parnumA)?0:1;
  //printf("maxstA=%.15G maxstB=%.15G\n", maxstA, maxstB);
  for (kk=0; kk < 6; kk++)
    vecini[kk] = vec[kk];
#if 0
  printf("$$$PR VECA=(%.15G,%.15G, %.15G)\n", vec[0], vec[1], vec[2]);
  printf("$$$PR VECB=(%.15G,%.15G, %.15G)\n", vec[3], vec[4], vec[5]);
#endif

  for (k1 = 0; k1 < 4; k1++)
    {
      nth1  = ellips_mesh[mA][*th1][*phi1].neigh[k1].i;
      nphi1 = ellips_mesh[mA][*th1][*phi1].neigh[k1].j;
      if (nth1==-1 || nphi1==-1)
	continue;
      for (kk=0; kk < 3; kk++) 
	{
	  xmeshP[kk] = ellips_mesh[mA][nth1][nphi1].point[kk];
	  dxP[kk] = xmeshP[kk] - ellips_mesh[mA][*th1][*phi1].point[kk]; 
	}
      body2labR(icg, dxP, dx, rA, RtA);
#if 0
      body2lab(jcg, xmeshP, xmesh, rA, RtA);
      for (kk=0; kk < 3; kk++)
	{
	  dx[kk] = xmesh[kk] - vec[kk];
	}
#endif
      if (calc_norm(dxP) > 0.7)
      	{
	  printf("nth1 nphi1 = %d %d *th1 *phi1 = %d %d\n", nth1, nphi1, *th1, *phi1);
	  printf("norm(dxP) %.15G\n",calc_norm(dxP));   
	  printf("norm(dx) %.15G\n",calc_norm(dx));   
	}
#if 0
      for (kk=0; kk < 3; kk++) 
	{
	  dxP[kk] = ellips_mesh[mA][nth1][nphi1].point[kk] -
	    ellips_mesh[mA][*th1][*phi1].point[kk];
	  if (fabs(dxP[kk]) > 1.0)
	    printf("dxP[%d]: %.15G\n", kk, dxP[kk]);   
	}
      body2labR(icg, dxP, dx, rA, RtA);
#endif
      //normdx = calc_norm(dx);
      //for (kk=0; kk < 3; kk++) 
      //dxN[kk] = dx[kk] / normdx;

      sp = scalProd(grad, dx);
      if (firstspA || sp > spmaxA)
	{
	  firstspA=0;
	  cth1 = nth1;
	  cphi1 = nphi1;
	  spmaxA = sp;
	  for (kk=0; kk < 3; kk++)
	   {
	     xmeshA[kk] = xmesh[kk];
	     cdxA[kk] = dx[kk];
	   }
	}
    }
  for (k1 = 0; k1 < 4; k1++)
    {
      nth2  = ellips_mesh[mB][*th2][*phi2].neigh[k1].i;
      nphi2 = ellips_mesh[mB][*th2][*phi2].neigh[k1].j;
      if (nth2==-1 || nphi2==-1)
	continue;
      for (kk=0; kk < 3; kk++) 
	{
	  xmeshP[kk] = ellips_mesh[mB][nth2][nphi2].point[kk];
	  dxP[kk] = xmeshP[kk] - ellips_mesh[mB][*th2][*phi2].point[kk]; 
	}
      body2labR(jcg, dxP, dx, rB, RtB);
#if 0
      body2lab(jcg, xmeshP, xmesh, rB, RtB);
      for (kk=0; kk < 3; kk++)
	{
	  dx[kk] = xmesh[kk] - vec[kk+3];
	}
#endif
      if (calc_norm(dxP) > 0.7)
	{
	  printf("nth2 nphi2 = %d %d *th2 *phi2 = %d %d\n", nth2, nphi2, *th2, *phi2);
	  printf("norm(dxP): %.15G\n", calc_norm(dxP));   
	  printf("norm(dx): %.15G\n", calc_norm(dx));   
	}
      //normdx = calc_norm(dx);
      //for (kk=0; kk < 3; kk++) 
      //dxN[kk] = dx[kk] / normdx;
      sp = scalProd(&grad[3], dx);
      //printf("B sp=%.15G\n", sp);
      if (firstspB || sp > spmaxB)
	{
	  firstspB = 0;
	  cth2 = nth2;
	  cphi2 = nphi2;
	  spmaxB = sp;
	  for (kk=0; kk < 3; kk++)
	    {
	      cdxB[kk] = dx[kk];
	      xmeshB[kk]= xmesh[kk];
	    }
	}
    }
  /* calcola la nuova distanza, il nuovo gradiente e le nuove coordinate del punto */
  //printf("nth2 nphi2 = %d %d *th2 *nphi2 = %d %d\n", nth2, nphi2, *th2, *phi2);
  //printf("nth1 nphi1 = %d %d *th1 *nphi1 = %d %d\n", nth1, nphi1, *th1, *phi1);
  for (kk = 0; kk < 3; kk++)
    {
      //printf("dxA[%d] = %.15G\n", kk, cdxA[kk]);
      vec[kk] += cdxA[kk];//xmeshA[kk];//cdxA[kk];
      vec[kk+3]+= cdxB[kk];//xmeshB[kk];//cdxB[kk];
    }
#if 0
    {
      double VP[6], VL[6];
      for (kk=0; kk < 3; kk++)
	{
	  VP[kk] = ellips_mesh[mA][*th1][*phi1].point[kk]; 
	  VP[kk+3]=ellips_mesh[mB][*th2][*phi2].point[kk];  
	}
      body2lab(icg, VP, VL, rA, RtA);
      body2lab(jcg, &VP[3], &VL[3], rB, RtB);
      printf("$$$ VA= (%.15G,%.15G, %.15G)\n", VL[0], VL[1], VL[2]);
      printf("$$$ VB= (%.15G,%.15G, %.15G)\n", VL[3], VL[4], VL[5]);
      printf("$$$ VECA=(%.15G,%.15G, %.15G)\n", vec[0], vec[1], vec[2]);
      printf("$$$ VECB=(%.15G,%.15G, %.15G)\n", vec[3], vec[4], vec[5]);
    }
#endif
  distSqold = *distSq;
  *distSq = 0;
  for (kk = 0; kk < 3; kk++)
    {
      *distSq += Sqr(vec[kk+3]-vec[kk]);
    }
  //printf("distSqold=%.15G distSq=%.15G\n",sqrt(fabs(distSqold)), sqrt(fabs(*distSq)));
  if (*distSq > distSqold)
    {
      for (kk=0; kk < 6; kk++)
	vec[kk] = vecini[kk];
      *distSq = distSqold;
      return 1;
    }
  if (calc_sign)
    {
      *S = get_sign(vec);
      *distSq *= *S;
    }
  for (kk = 0; kk < 3; kk++)
    {
      grad[kk] = *S*(vec[kk+3]-vec[kk]);
      grad[kk+3] = -grad[kk];
    }
  *th1 = cth1;
  *phi1 = cphi1;
  *th2 = cth2;
  *phi2 = cphi2;
  //printf("USCENDO *th1 *phi1 = %d %d *th2 *phi2 = %d %d\n", *th1, *phi1, *th2, *phi2);
  return 0;	
}
extern double *maxax;

void findminMesh(double *vec)
{
  int mA, mB, kk, th1, th2, phi1, phi2, calc_sign=0, its;
  double S, angs[4], vecP[6], sinth, grad[6], maxstA, maxstB, distSq;
  const double TWOPI = 2.0*pi;
  lab2body(icg, vec, vecP, rA, RtA);
  lab2body(jcg, &vec[3], &vecP[3], rB, RtB);
  /* angs = (theta1,phi1,theta2,phi2) 
   * 0 < theta < PI
   * 0< phi < 2PI */ 
  angs[0] = acos(vecP[2]/axc[icg]);
  angs[1] = acos(vecP[0]/axa[icg]/sin(angs[0]));
  sinth = vecP[1]/axb[icg]/sin(angs[0]);
  callsfrprmn++;
  if (sinth < 0)
    angs[1] = 2.0*pi-angs[1];
  angs[2] = acos(vecP[5]/axc[jcg]);
  angs[3] = acos(vecP[3]/axa[jcg]/sin(angs[2]));
  sinth = vecP[4]/axb[jcg]/sin(angs[2]);
  if (sinth < 0)
    angs[3] = 2.0*pi-angs[3];
#if 0 
    {
      double vv[6], vl[6];
      printf("INIZIO VECA=(%.15G,%.15G, %.15G)\n", vec[0], vec[1], vec[2]);
      printf("INIZIO VECB=(%.15G,%.15G, %.15G)\n", vec[3], vec[4], vec[5]);
      angs2coord(angs, vv);
      body2lab(icg, vv, vl, rA, RtA);
      body2lab(jcg, &vv[3], &vl[3], rB, RtB);
      printf("INIZIO2 VECA=(%.15G,%.15G, %.15G)\n", vl[0], vl[1], vl[2]);
      printf("INIZIO2 VECB=(%.15G,%.15G, %.15G)\n", vl[3], vl[4], vl[5]);
   }
#endif
  /* determino lo starting point sulla mesh */	
  th1 = OprogStatus.n1*angs[0]/TWOPI; 
  phi1 = OprogStatus.n2*angs[1]/TWOPI;
  th2 = OprogStatus.n1*angs[2]/TWOPI; 
  phi2 = OprogStatus.n2*angs[3]/TWOPI;
  if (th1 == 0) 
    th1 = 1;
  if (th2 == 0)
    th2 = 1;
  if (th1 == OprogStatus.n1/2)
    th1 = OprogStatus.n1/2-1;
  if (th2 == OprogStatus.n1/2)
    th2 = OprogStatus.n1/2-1;
#if 0
  angs[0] = th1 * TWOPI /OprogStatus.n1;
  angs[1] = phi1* TWOPI /OprogStatus.n2;
  angs[2] = th2 * TWOPI /OprogStatus.n1;
  angs[3] = phi2* TWOPI /OprogStatus.n2;
#endif
  mA = (icg<Oparams.parnumA)?0:1;
  mB = (jcg<Oparams.parnumA)?0:1;
  for (kk = 0; kk < 3; kk++)
    {
      vecP[kk] = ellips_mesh[mA][th1][phi1].point[kk];
      vecP[kk+3] = ellips_mesh[mB][th2][phi2].point[kk];
    }
  body2lab(icg, vecP, vec, rA, RtA);
  body2lab(jcg, &vecP[3], &vec[3], rB, RtB);
  //printf("INIZIO: (%d,%d,%d,%d)\n", th1, phi1, th2, phi2);
  /* maggiorazione degli step sui due ellissoidi */
  if (OprogStatus.n1 > OprogStatus.n2)
    maxstA = maxax[icg]*TWOPI/OprogStatus.n2;
  else
    maxstA = maxax[icg]*TWOPI/OprogStatus.n1;

  if (OprogStatus.n1 > OprogStatus.n2)
    maxstB = maxax[jcg]*TWOPI/OprogStatus.n2;
  else
    maxstB = maxax[jcg]*TWOPI/OprogStatus.n2;
  /* search begin */
  S = get_sign(vec); 
  distSq = 0;
  for (kk = 0; kk < 3; kk++)
    {
      grad[kk] = S*(vec[kk+3]-vec[kk]);
      grad[kk+3] = -S*(vec[kk+3]-vec[kk]); 
      distSq += S*Sqr(vec[kk+3]-vec[kk]);
    }
  if (S < 0)
    calc_sign=1;
#if 0
  printf("INIZIO VECA=(%.15G,%.15G, %.15G)\n", vec[0], vec[1], vec[2]);
  printf("INIZIO VECB=(%.15G,%.15G, %.15G)\n", vec[3], vec[4], vec[5]);
#endif
  //printf(">>> distSq: %.15G\n", sqrt(distSq));
  its = 0;
  while (1)
    {
      its++;
      itsfrprmn++;
      /* calcola il segno se necessario */
      if (!calc_sign && distSq < Sqr(maxstA+maxstB))
	{
	  S = get_sign(vec);
	  calc_sign = 1;
	}
      if (choose_neighbour2(grad, &th1, &phi1, &th2, &phi2, &distSq, &S, 
      			   maxstA, maxstB, vec, calc_sign))
	{
#if 0
	    /* we're done here! */
	    for (kk=0; kk < 3; kk++)
	      {
		vecP[kk] = ellips_mesh[mA][th1][phi1].point[kk];
		vecP[kk+3] = ellips_mesh[mB][th2][phi2].point[kk];
	      }
	    /* torno nel riferimento del laboratorio */
	    body2lab(icg, vecP, vec, rA, RtA);
	    body2lab(jcg, &vecP[3], &vec[3], rB, RtB);
	    //printf("BOH distSq=%.15G\n", sqrt(distSq));
	    distSq = 0;
	    for (kk = 0; kk < 3; kk++)
	      {
		distSq += S*Sqr(vec[kk+3]-vec[kk]);
	      }
	    //printf("FINE S=%f its= %d dist=%.15G\n", S,its, sqrt(distSq));
#endif
	    break;
	  } 
    }
}
#define ZBFACTOR 1.6
#define ZBNTRY 50 
int zbrac(double (*func)(double), double *x1, double *x2) 
/* Given a function func and an initial guessed range x1 to x2, the routine expands the range 
 * geometrically until a root is bracketed by the returned values x1 and x2 (in which case 
 * zbrac returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0).
 */
{ 
 //void nrerror(char error_text[]); 
 int j; 
 double f1,f2; 
 if (*x1 == *x2) 
   {
     return 0;
   }
 //nrerror("Bad initial range in zbrac"); 
 f1=(*func)(*x1); 
 f2=(*func)(*x2); 
 for (j=1;j<=ZBNTRY;j++)
   { 
     if (f1*f2 < 0.0) 
       return 1; 
     if (fabs(f1) < fabs(f2)) 
       f1=(*func)(*x1 += ZBFACTOR*(*x1-*x2)); 
     else
       f2=(*func)(*x2 += ZBFACTOR*(*x2-*x1)); 
   } 
 return 0;
}
#define ITMAXZB 1000 
/* Maximum allowed number of iterations.*/
#define EPSP 3.0E-16/* Machine floating-point precision.*/
void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[], double xb2[], 
	   int *nb)
/* Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
 * spaced segments, and search for zero crossings of the function. nb is input as the maximum 
 * number of roots sought, and is reset to the number of bracketing pairs xb
 * 1[1..nb], xb2[1..nb] that are found. */
{
  int nbb,i; 
  double x,fp,fc,dx; 
  nbb=0; dx=(x2-x1)/n;
  /* Determine the spacing appropriate to the mesh.*/
  fp=(*fx)(x=x1); 
  for (i=0;i<n;i++) { /* Loop over all intervals*/
    x += dx;
    fc=(*fx)(x); 
    if (fp >= 0 && fc <= 0.0) 
      { /* If a sign change occurs then record values for the bounds.*/
	xb1[nbb]=x1; xb2[nbb]=x;
	nbb++;
	if(*nb == nbb) 
	  return;
      } 
    //fp=fc;
  } 
  *nb = nbb;
}
double zbrentRyck(double (*func)(double), double x1, double x2, double tol)
/* Using Brent s method, find the root of a function func known to lie between x1 and x2. 
 * The root, returned as zbrent, will be refined until its accuracy is tol.*/
{
  int iter; 
  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2; 
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm; 
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    {
      polinterrRyck = 1;
      return 0.0;
      //nrerror("Root must be bracketed in zbrent");
    }
  fc=fb;
  for (iter=0;iter<ITMAXZB;iter++) 
    { 
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
	{ 
	  c=a; /* Rename a, b, c and adjust bounding interval d.*/
	  fc=fa; e=d=b-a;
	} 
      if (fabs(fc) < fabs(fb)) 
	{
	  a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;
	}
      tol1=2.0*EPSP*fabs(b)+0.5*tol;
      /* Convergence check. */
      xm=0.5*(c-b); 
      if (fabs(xm) <= tol1 || fb == 0.0) 
	return b;
      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
	{
	  s=fb/fa;/* Attempt inverse quadratic interpolation.*/
	  if (a == c) 
    	    { 
	      p=2.0*xm*s; q=1.0-s;
	    } 
	  else 
	    { 
	      q=fa/fc; r=fb/fc; 
	      p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	      q=(q-1.0)*(r-1.0)*(s-1.0);
	    }
	  if (p > 0.0)
	    q = -q;  /* Check whether in bounds. */
	  p=fabs(p); 
	  min1=3.0*xm*q-fabs(tol1*q); 
	  min2=fabs(e*q);
	  if (2.0*p < (min1 < min2 ? min1 : min2)) 
	    {
	      e=d; /*Accept interpolation. */ 
	      d=p/q; 
	    } 
	  else 
	    { 
	      d=xm; /*Interpolation failed, use bisection.*/
	      e=d; 
	    } 
	} 
      else 
	{ 
	  /* Bounds decreasing too slowly, use bisection.*/
	  d=xm; e=d;
	} 
      a=b; /* Move last best guess to a. */
      fa=fb;
      if (fabs(d) > tol1) /* Evaluate new trial root.*/
	b += d; 
      else 
	b += SIGN(tol1,xm); 
      fb=(*func)(b); 
      if (polinterrRyck)
	return 0.0;
    } 

  polinterrRyck = 1;
  return 0.0;
  //nrerror("Maximum number of iterations exceeded in zbrent"); 
  return 0.0; /* Never get here.*/ 
}

double zbrent(double (*func)(double), double x1, double x2, double tol)
/* Using Brent s method, find the root of a function func known to lie between x1 and x2. 
 * The root, returned as zbrent, will be refined until its accuracy is tol.*/
{
  int iter; 
  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2; 
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm; 
  MD_DEBUG(printf("==>>>a=%f b=%f fa: %f fb: %f\n",a, b, fa, fb));
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    {
      polinterr = 1;
      //printf("bohbohboh\n");
      return 0.0;
      //nrerror("Root must be bracketed in zbrent");
    }
  fc=fb;
  for (iter=0;iter<ITMAXZB;iter++) 
    { 
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
	{ 
	  c=a; /* Rename a, b, c and adjust bounding interval d.*/
	  fc=fa; e=d=b-a;
	} 
      if (fabs(fc) < fabs(fb)) 
	{
	  a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;
	}
      tol1=2.0*EPSP*fabs(b)+0.5*tol;
      /* Convergence check. */
      xm=0.5*(c-b); 
      if (fabs(xm) <= tol1 || fb == 0.0) 
	{
	  return b;
	}
      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
	{
	  s=fb/fa;/* Attempt inverse quadratic interpolation.*/
	  if (a == c) 
    	    { 
	      p=2.0*xm*s; q=1.0-s;
	    } 
	  else 
	    { 
	      q=fa/fc; r=fb/fc; 
	      p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	      q=(q-1.0)*(r-1.0)*(s-1.0);
	    }
	  if (p > 0.0)
	    q = -q;  /* Check whether in bounds. */
	  p=fabs(p); 
	  min1=3.0*xm*q-fabs(tol1*q); 
	  min2=fabs(e*q);
	  if (2.0*p < (min1 < min2 ? min1 : min2)) 
	    {
	      e=d; /*Accept interpolation. */ 
	      d=p/q; 
	    } 
	  else 
	    { 
	      d=xm; /*Interpolation failed, use bisection.*/
	      e=d; 
	    } 
	} 
      else 
	{ 
	  /* Bounds decreasing too slowly, use bisection.*/
	  d=xm; e=d;
	} 
      a=b; /* Move last best guess to a. */
      fa=fb;
      if (fabs(d) > tol1) /* Evaluate new trial root.*/
	b += d; 
      else 
	b += SIGN(tol1,xm); 
      fb=(*func)(b); 
      if (polinterr)
	return 0.0;
    } 

  polinterr = 1;
  return 0.0;
  //nrerror("Maximum number of iterations exceeded in zbrent"); 
  return 0.0; /* Never get here.*/ 
}
void polintRyck(double xain[], double yain[], int n, double x, double *y, double *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y,
 * and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xai) = yai, 
 * i = 1, . . . , n, then the returned value y = P(x).*/
{ 
  int i,m,ns=1; 
  double den,dif,dift,ho,hp,w, xa[4], ya[4];
  double c[4], d[4];
  for (i=0; i < n; i++)
    {
      xa[i+1] = xain[i];
      ya[i+1] = yain[i];
    }
  dif=fabs(x-xa[1]); 
  //c=vector(n); 
  //d=vector(n); 
  for (i=1;i<=n;i++) 
    { 
      /* Here we find the index ns of the closest table entry,*/
      if ( (dift=fabs(x-xa[i])) < dif) 
	{ 
	  ns=i; 
	  dif=dift;
	} 
      c[i]=ya[i];
      /* and initialize the tableau of c s and d s.*/
      d[i]=ya[i]; 
    } 
  *y=ya[ns--];
  /* This is the initial approximation to y.*/
  for (m=1;m<n;m++) 
    { 
      /* For each column of the tableau,*/
      for (i=1;i<=n-m;i++)
	{
	  /* we loop over the current c s and d s and update them.*/
	  ho=xa[i]-x; 
	  hp=xa[i+m]-x; 
	  w=c[i+1]-d[i]; 
	  if ( (den=ho-hp) == 0.0) 
	    {
	      polinterrRyck=1;
	      return ;
	      //nrerror("Error in routine polint"); 
	    }
	  /* This error can occur only if two input xa s are (to within roundoff)*/
	  den=w/den; d[i]=hp*den; 
	  /*Here the c s and d s are updated. */
	  c[i]=ho*den; 
	} 
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--])); 
      /* After each column in the tableau is completed, we decide which correction, 
       * c or d, we want to add to our accumulating value of y, i.e., which path to take through the tableau 
       * forking up or down. We do this in such a way as to take the most  straight line  route through the 
       * tableau to its apex, updating ns accordingly to keep track of where we are. 
       * This route keeps the partial approximations centered (insofar as possible) on the target x. 
       * The last dy added is thus the error indication. */
    } 
  //free_vector(d); 
  //free_vector(c); 
}

void polint(double xain[], double yain[], int n, double x, double *y, double *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y,
 * and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xai) = yai, 
 * i = 1, . . . , n, then the returned value y = P(x).*/
{ 
  int i,m,ns=1; 
  double den,dif,dift,ho,hp,w, xa[4], ya[4];
  double c[4], d[4];
  for (i=0; i < n; i++)
    {
      xa[i+1] = xain[i];
      ya[i+1] = yain[i];
    }
  dif=fabs(x-xa[1]); 
  //c=vector(n); 
  //d=vector(n); 
  for (i=1;i<=n;i++) 
    { 
      /* Here we find the index ns of the closest table entry,*/
      if ( (dift=fabs(x-xa[i])) < dif) 
	{ 
	  ns=i; 
	  dif=dift;
	} 
      c[i]=ya[i];
      /* and initialize the tableau of c s and d s.*/
      d[i]=ya[i]; 
    } 
  *y=ya[ns--];
  /* This is the initial approximation to y.*/
  for (m=1;m<n;m++) 
    { 
      /* For each column of the tableau,*/
      for (i=1;i<=n-m;i++)
	{
	  /* we loop over the current c s and d s and update them.*/
	  ho=xa[i]-x; 
	  hp=xa[i+m]-x; 
	  w=c[i+1]-d[i]; 
	  if ( (den=ho-hp) == 0.0) 
	    {
	      polinterr=1;
	      return ;
	      //nrerror("Error in routine polint"); 
	    }
	  /* This error can occur only if two input xa s are (to within roundoff)*/
	  den=w/den; d[i]=hp*den; 
	  /*Here the c s and d s are updated. */
	  c[i]=ho*den; 
	} 
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--])); 
      /* After each column in the tableau is completed, we decide which correction, 
       * c or d, we want to add to our accumulating value of y, i.e., which path to take through the tableau 
       * forking up or down. We do this in such a way as to take the most  straight line  route through the 
       * tableau to its apex, updating ns accordingly to keep track of where we are. 
       * This route keeps the partial approximations centered (insofar as possible) on the target x. 
       * The last dy added is thus the error indication. */
    } 
  //free_vector(d); 
  //free_vector(c); 
}
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void gaussj(double **a, int n, double *bb) 
/* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n] 
 * is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. 
 * On output, a is replaced by its matrix inverse,
 * and b is replaced by the corresponding set of solution vectors. */
{ 
  const int m=1;
  int *indxc,*indxr,*ipiv;
  int i,icol=0,irow=0,j,k,l,ll; 
  double **b,big,dum,pivinv,temp; 
  b = matrix(n,1);
  for (i=0; i < n; i++)
    b[i][0] = bb[i];
  indxc=ivector(n); 
  /* The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting. */
  indxr=ivector(n); 
  ipiv=ivector(n); 
  for (j=0;j<n;j++) 
    ipiv[j]=0; 
  for (i=0;i<n;i++) 
    { /* This is the main loop over the columns to be reduced. */
      big=0.0; 
      for (j=0;j<n;j++) /* This is the outer loop of the search for a pivot element.*/
      if (ipiv[j] != 1) 
	for (k=0;k<n;k++)
	  { 
	    if (ipiv[k] == 0)
	      { 
		if (fabs(a[j][k]) >= big)
		  { 
		    big=fabs(a[j][k]); 
		    irow=j;
		    icol=k; 
		  }
	      }
	  } 
      ++(ipiv[icol]);
      /* We now have the pivot element, so we interchange rows, if needed, to put the pivot 
       * element on the diagonal. The columns are not physically interchanged, 
       * only relabeled: indxc[i], the column of the ith pivot element, is the ith column 
       * that is reduced, while indxr[i] is the row in which that pivot element was originally
       * located. If indxr[i]  = indxc[i] there is an implied column interchange. 
       * With this form of bookkeeping, the solution b s will end up in the correct order, 
       * and the inverse matrix will be scrambled by columns. */
      if (irow != icol) 
	{ 
	  for (l=0;l<n;l++) 
	    SWAP(a[irow][l],a[icol][l]); 
	  for (l=0;l<m;l++) 
	    SWAP(b[irow][l],b[icol][l]);
	}
      indxr[i]=irow;
      /* We are now ready to divide the pivot row by the pivot element,
       * located at irow and icol. */
      indxc[i]=icol; 
      if (a[icol][icol] == 0.0)
	nrerror("gaussj: Singular Matrix"); 
      pivinv=1.0/a[icol][icol]; 
      a[icol][icol]=1.0; 
      for (l=0;l<n;l++) 
	a[icol][l] *= pivinv; 
      for (l=0;l<m;l++) 
	b[icol][l] *= pivinv;
      for (ll=0;ll<n;ll++) 
	/* Next, we reduce the rows...*/
	if (ll != icol) 
	  { 
	    /* ...except for the pivot one, of course.*/
	    dum=a[ll][icol];
	    a[ll][icol]=0.0;
	    for (l=0;l<n;l++) 
	      a[ll][l] -= a[icol][l]*dum; 
	    for (l=0;l<m;l++) 
	      b[ll][l] -= b[icol][l]*dum; 
	  }
    } 
  /* This is the end of the main loop over columns of the reduction. It only remains to 
   * unscramble the solution in view of the column interchanges. We do this by 
   * interchanging pairs of columns in the reverse order that the permutation was built up.
   * */
  for (l=n-1;l>=0;l--) 
    { 
      if (indxr[l] != indxc[l]) 
	for (k=0;k<n;k++) 
	  SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
  /* And we are done.*/
  for (i=0; i < n; i++)
    bb[i] = b[i][0];
  free_matrix(b,n);
  free_ivector(ipiv); 
  free_ivector(indxr); 
  free_ivector(indxc);
} 
void fdjacFD(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3]);
void ludcmp(double **a, int n,  int* indx, double* d, int *ok)
{
  /* A[i][j] = Aij 
   * A x = b  
   * per semplicità nel seguito si assume che l'ordine della matrice è 3 */
  int i,imax=-1,j,k;
  double big,dum,sum,temp; 
  double vv[MD_NBMAX]; /* vv stores the implicit scaling of each row.*/
  /*vv = vector(1,n);*/
  *d=1.0; /* No row interchanges yet. */
  *ok = 0;
  for (i=0;i<n;i++) 
    { 
      /* Loop over rows to get the implicit scaling information.*/ 
      big=0.0; 
      for (j=0;j<n;j++)
	{
	  if ((temp=fabs(a[i][j])) > big) big=temp; 
	}
      if (big == 0.0)
	{
	  printf("ERROR: Singular matrix in routine ludcmp\n"); 
	  *ok = 1;
	  return;
	}
      /* No nonzero largest element. */
      vv[i]=1.0/big; /* Save the scaling.*/
    } 
  for (j=0;j<n;j++) 
    { /* This is the loop over columns of Crout s method.*/
      for (i=0;i<j;i++) 
	{ 
	  /* This is equation (2.3.12) except for i = j. */
	  sum=a[i][j]; 
	  for (k=0;k<i;k++) 
	    sum -= a[i][k]*a[k][j]; 
	  a[i][j]=sum; 
	} 
      big=0.0; /* Initialize for the search for largest pivot element. */ 
      for (i=j;i<n;i++) 
	{ 
	  /* This is i = j of equation (2.3.12) and i = j+1. . .N of equation (2.3.13).*/
	  sum=a[i][j]; 
	  for (k=0;k<j;k++)
	    sum -= a[i][k]*a[k][j]; 
	    a[i][j]=sum; 
	    if ( (dum=vv[i]*fabs(sum)) >= big) 
	      { 
		/* Is the  gure of merit for the pivot better than the best so far? */
		big=dum; imax=i; 
	      } 
	} 
      if (j != imax) 
	{ 
	  /* Do we need to interchange rows? */
	  for (k=0;k<n;k++) 
	    { 
	      /* Yes, do so...*/ 
	      dum=a[imax][k]; 
	      a[imax][k]=a[j][k]; 
	      a[j][k]=dum; 
	    } 
	  *d = -(*d); 
	  /* ...and change the parity of d. */ 
	  vv[imax]=vv[j]; 
	  /* Also interchange the scale factor.*/ 
	} 
      indx[j]=imax; 
      if (a[j][j] == 0.0) 
	a[j][j]=TINY; 
      /* If the pivot element is zero the matrix is singular 
       * (at least to the precision of the algorithm). 
       * For some applications on singular matrices, 
       * it is desirable to substitute TINY for zero. */ 
      if (j != n) 
	{ 
	  /* Now,  nally, divide by the pivot element.*/
	  dum=1.0/(a[j][j]); 
	  for (i=j+1;i<n;i++) a[i][j] *= dum; 
	} 
    } 
  /* Go back for the next column in the reduction.*/
  /*free_vector(vv,1,n); */
}

void lubksb(double **a, int n, int* indx, double *b)
{ 
  int i,ii=0,ip,j; 
  double sum; 
  for (i=0;i<n;i++) 
    { 
      /* When ii is set to a positive value, it will become the index of the  
       * rst nonvanishing element of b. Wenow do the forward substitution,
       * equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go. */
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i]; 
      if (ii>-1) 
	for (j=ii;j<=i-1;j++) 
	  sum -= a[i][j]*b[j]; 
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
      for (j=i+1;j<n;j++) 
	sum -= a[i][j]*b[j]; b[i]=sum/a[i][i]; 
      /* Store a component of the solution vector X. */ 
    } /* All done! */
}
#ifdef MD_USE_LAPACK
void wrap_dgesv(double **a, double *x, int n, int *ok)
{
  double AT[MD_NBMAX*MD_NBMAX];
  int i, j, c1, c2, ok, pivot[MD_NBMAX];
  for (i=0; i<n; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<n; j++) AT[j+n*i]=a[j][i];		
    }						
  c1 = n;
  c2 = 1;
  dgesv_(&c1, &c2, AT, &c1, pivot, x, &c1, ok);      
}
#endif
int SolveLineq (double **a, double *x, int n) 
{
  int indx[MD_NBMAX], ok;
  double dd;
#ifdef MD_USE_LAPACK
  wrap_dgesv(a, x, n, &ok);
#else
  ludcmp(a, n, indx, &dd, &ok);
  if (ok==1)
    return 1;
  lubksb(a, n, indx, x);
#endif
  return 0;
}

void InvMatrix(double **a, double **b, int NB)
{
  int m1, m2, indx[MD_NBMAX], ok; 
  double col[MD_NBMAX];
  double d;
  ludcmp(a, NB, indx, &d, &ok); 
  for(m2=0;m2<NB;m2++) 
    { 
      for(m1=0;m1<NB;m1++) 
	col[m1]=0.0; 
      col[m2]=1.0; 
      lubksb(a, NB, indx, col);
      for(m1=0;m1<NB;m1++) 
	 b[m1][m2]=col[m1]; 
    }
}

#define ALF 1.0e-4 /* Ensures sufficient decrease in function value.*/
#define TOLX 1.0E-12//1.0e-7 /* Convergence criterion on  x.*/ 
#define TOLX2 1.E-6
#define TOLXD 1.0E-7
#define MAXITS 50 // se le particelle non si urtano il newton-raphson farà MAXITS iterazioni
#define MAXITS2 20
#define TOLF 1.0e-9// 1.0e-4
#define TOLF2 1.0E-3
#define TOLFD 1.0E-6
#define TOLMIN 1.0E-7//1.0e-6 
#define STPMX 100.0
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], 
	    double *f, double stpmax, int *check, 
	    double (*func)(double [], int, int, double[]), int iA, int iB, double shift[3],
	    double tolx)
/*
   Given an n-dimensional point xold[1..n], the value of the function and gradient there, 
   fold and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p
   from xold where the function func has decreased  "sufficiently".  The new function value is 
   returned in f. stpmax is an input quantity that limits the length of the steps so that 
   you do not try to evaluate the function in regions where it is unde ned or subject 
   to overflow.
   p is usually the Newton direction. The output quantity check is false (0) on a normal exit. 
   It is true (1) when x is too close to xold. In a minimization algorithm, this usually 
   signals convergence and can be ignored. 
   However, in a zero-finding algorithm the calling program 
   should check whether the convergence is spurious. Some  difficult  problems may require 
   double precision in this routine.*/
{
  int i; 
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0,rhs1,rhs2,slope,sum,temp, test,tmplam; 
  *check=0; 
  for (sum=0.0,i=0;i<n;i++) 
    sum += p[i]*p[i]; 
  sum=sqrt(sum); 
  if (sum > stpmax) 
    for (i=0;i<n;i++) 
      p[i] *= stpmax/sum; /*Scale if attempted step is too big.*/ 
  for (slope=0.0,i=0;i<n;i++) 
    slope += g[i]*p[i]; 
  if (slope >= 0.0) 
    nrerror("Roundoff problem in lnsrch."); 
  test=0.0; /*Compute lambda_min.*/
  for (i=0;i<n;i++) 
    {
      temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0); 
      if (temp > test) 
	test=temp; 
    } 
  alamin=tolx/test; alam=1.0;
  for (;;) 
    { 
      for (i=0;i<n;i++) 
	x[i]=xold[i]+alam*p[i]; 
      *f=(*func)(x,iA,iB,shift); 
      if (alam < alamin) 
	{ /* Convergence on  x. For zero  nding, the calling program 
	     should verify the convergence.*/ 
	  for (i=0;i<n;i++) 
	    x[i]=xold[i]; 
	  *check=1; 
	  return;
	}
      else if (*f <= fold+ALF*alam*slope) 
	return; 
	/* Su cient function decrease.*/
      else 
	{ /* Backtrack. */
	  if (alam == 1.0) 
	    tmplam = -slope/(2.0*(*f-fold-slope));/* First time.*/
	  else
	    { /* Subsequent backtracks.*/
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2); 
	      if (a == 0.0) 
		tmplam = -slope/(2.0*b); 
	      else 
		{
		  disc=b*b-3.0*a*slope; 
		  if (disc < 0.0) 
		    tmplam=0.5*alam; 
		  else if (b <= 0.0) 
		    tmplam=(-b+sqrt(disc))/(3.0*a); 
		  else 
		    tmplam=-slope/(b+sqrt(disc)); 
		} 
	      if (tmplam > 0.5*alam) 
		tmplam=0.5*alam; /* lambda <= 0.5 lambda_1.*/
	    } 
	}
      alam2=alam;
      f2 = *f; 
      alam=FMAX(tmplam,0.1*alam); /* lambda >= 0.1 lambda_1.*/
    }/* Try again.*/
}
int nn, nn2, nnD; /* Global variables to communicate with fmin.*/
double *fvec, *fvecG, *fvecD; 
#ifdef MD_GLOBALNR2
#define FREERETURN {MD_DEBUG(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d\n", x[0], x[1], x[2], x[3], x[4], test, its, *check));\
free_vector(fvec);free_vector(xold);free_vector(g2); free_vector(xold2); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);free_vector(fvecG);return;}
#else
#define FREERETURN {MD_DEBUG10(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d fvec=(%.15G,%.15G,%.15G,%.15G,%.15G)\n", x[0], x[1], x[2], x[3], x[4], test, its, *check, fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));\
free_vector(fvec);free_vector(xold); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);free_vector(fvecG);return;}
#endif
#define FREERETURND {MD_DEBUG(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d\n", x[0], x[1], x[2], x[3], x[4], test, its, *check));\
free_vector(fvecD);free_vector(xold); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);return;}

void (*nrfuncv)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncv2)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncvD)(int n, double v[], double fvec[], int i, int j, double shift[3]);


extern void fdjac(int n, double x[], double fvec[], double **fjac, 
		  void (*vecfunc)(int n, double v[], double fvec[], int i, int j, double shift[3]), int iA, int iB, double shift[3]); 
double fminD(double x[], int iA, int iB, double shift[3]);
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], double *f, 
	    double stpmax, int *check, double (*func)(double [], int, int, double []),
	    int iA, int iB, double shift[3], double tolx);
void lubksb(double **a, int n, int *indx, double b[]); 
void ludcmp(double **a, int n, int *indx, double *d, int *ok); 
extern void funcs2beZeroedGuess(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroed(int n, double x[], double fvec[], int i, int j, double shift[3]);

extern void upd2tGuess(int i, int j, double shift[3], double tGuess);

#define EPS 3.0E-8 /* Approximate square root of the machine precision.*/

extern int nn, nn2; 
extern void (*nrfuncv)(int n, double v[], double f[], int iA, int iB, double shift[3]); 
extern void (*nrfuncv2)(int n, double v[], double f[], int iA, int iB, double shift[3]); 
