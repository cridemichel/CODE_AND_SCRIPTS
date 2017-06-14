#define TINY 1E-20
#define MD_NBMAX 8 
#define MD_NNLPLANES
#include<mdsimul.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define MD_DEBUG10(x) 
#define MD_DEBUG18(x)
#define MD_DEBUG20(x) 
#define MD_NEW_NR_CHECKS
#ifdef EDHE_FLEX 
extern int *typeOfPart;
#endif
int ncom;
double (*nrfunc)(double []); 
double f1dim(double x);
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));
double fminMD(double x[], int iA, int iB, double shift[3]);
void polintRyck(double xain[], double yain[], int n, double x, double *y, double *dy);
extern void fdjacDistNegNeighPlane5(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int), int iA);
extern void fdjacDistNegNeighPlane(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA);
extern void funcs2beZeroedNeighPlane(int n, double x[], double fvec[], int i);
extern void funcs2beZeroedNeigh(int n, double x[], double fvec[], int i);
extern void fdjacNeighPlane(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int), int iA);

int *ivector(int n)
{
  return calloc(n, sizeof(int)); 
}
double *vector(int n)
{
  return calloc(n, sizeof(double)); 
}
void projectgrad(double *p, double *xi, double *gradf, double *gradg);
void projonto(double* ri, double *dr, double* rA, double **Xa, double *gradf, double *sfA, double dist);
void polint(double xain[], double yain[], int n, double x, double *y, double *dy);
extern void funcs2beZeroedDistNegNeighPlane5(int n, double x[], double fvec[], int i);
extern void funcs2beZeroedDistNegNeighPlane(int n, double x[], double fvec[], int i);
extern void funcs2beZeroedDistNegNeigh(int n, double x[], double fvec[], int i);
void fdjacDistNegNeigh(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int), int iA);
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
extern void funcs2beZeroedDistNegNew(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroedDistNeg5(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void fdjacDistNeg5(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3], double *fx, double *gx);
extern void fdjacDistNeg(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3], double *fx, double *gx);
extern void fdjacDistNegNew(int n, double x[], double fvec[], double **df, 
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
long long accngA=0, accngB=0;
double xicom[8], pcomI[8], pcom[8], pcom2[8], xi[8], G[8], H[8], grad[8];//, vec[6];
double Ftol, Epoten, Emin, fnorm;
int cghalfspring, icg, jcg, doneryck;
double shiftcg[3], lambdacg, minaxicg, minaxjcg;
double gradfG[3], gradgG[3], dxG[6];
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, rA[3], rB[3], **RtA, **RtB;
extern double minaxA, minaxB, minaxAB;
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
int brentTooManyIter;
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
  brentTooManyIter=0;
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
  printf("Too many iterations in brent!\n");
  brentTooManyIter=1;
  //nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}
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
double f1dimConstr(double x) 
  /*Must accompany linmin.*/
{
  int j; 
  double f, xt[6];
  // xt=vector(1,ncom);
  for (j=0;j<ncom;j++) 
    {
      dxG[j]=x*xicom[j];
    }
  projonto(pcom, dxG, rA, Xa, gradfG, &sfA, 0);
  projonto(&pcom[3], &dxG[3], rB, Xb, gradgG, &sfB, 0);
  for (j=0;j<ncom;j++) 
    {
      xt[j]=pcom[j]+dxG[j]; 
    }
  f=(*nrfunc)(xt); 
  // free_vector(xt,1,ncom); 
  return f;
}

double df1dim(double x) 
{ 
  int j;
  double df1=0.0; 
  double xt[8], df[8];
  //xt=vector(1,ncom); df=vector(1,ncom);
  for (j=0;j<ncom;j++) 
    xt[j]=pcom[j]+x*xicom[j]; 
  (*nrdfun)(xt,df); 
  for (j=0;j<ncom;j++) 
  df1 += df[j]*xicom[j]; 
  //free_vector(df,1,ncom); free_vector(xt,1,ncom); 
  return df1;
}
void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),
	     void (*dfunc)(double [], double [])) 
  /*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and resets p to where 
   * the function func(p) takes on a minimum along the direction xi from p, and replaces xi by the actual vector
   * displacement that p was moved. Also returns as fret the value of func at the returned location p.
   * This is actually all accomplished by calling the routines mnbrak and dbrent. */
{ 
  double dbrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double), double tol,
		double *xmin);
  const double TOLLINMIN = 1E-10;
  double f1dim(double x); 
  double df1dim(double x); 
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
  int j; 
  double xx,xmin,fx,fb,fa,bx,ax; 
  ncom=n; /*Define the global variables.*/
 // pcom=vector(1,n);
 // xicom=vector(1,n); 
  nrfunc=func; nrdfun=dfunc;
  for (j=0;j<n;j++) 
    { 
      pcom[j]=p[j];
      xicom[j]=xi[j]; 
    } 
  ax=0.0; /*Initial guess for brackets.*/
  xx=1.0; 
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=dbrent(ax,xx,bx,f1dim,df1dim,TOLLINMIN,&xmin);
  for (j=0;j<n;j++) 
    { /* Construct the vector results to return. */
      xi[j] *= xmin; p[j] += xi[j]; 
    } 
  //free_vector(xicom,1,n); free_vector(pcom,1,n); 
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
  double sin1, sin3;
  sin1 = sin(angs[1]);
  p[0] = axa[icg]*cos(angs[0])*sin1;
  p[1] = axb[icg]*sin(angs[0])*sin1;
  p[2] = axc[icg]*cos(angs[1]);
  
  sin3 = sin(angs[3]);
  p[3] = axa[jcg]*cos(angs[2])*sin3;
  p[4] = axb[jcg]*sin(angs[2])*sin3;
  p[5] = axc[jcg]*cos(angs[3]);
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

void powellmethod(double *vec)
{
  double Fret, vecP[6], angs[4], sinth;
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
double powdirsP[6][6]={{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},
    {0,0,0,0,1,0},{0,0,0,0,0,1}};
double  cgfunc(double *vec);

void powellmethodPenalty(double *vec)
{
  int k1, k2, iter; 
  double Fret;
  for (k1=0; k1 < 6; k1++)
    for (k2=0; k2 < 6; k2++)
      powdirs[k1][k2] = powdirsP[k1][k2];
  powell(vec, powdirs, 6, OprogStatus.tolSD, &iter, &Fret, cgfunc);
}
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
extern double costolAngSD;
void projonto(double* ri, double *dr, double* rA, double **Xa, double *gradf, double *sfA, double dist)
{
  int kk, its, done=0, k1, k2, MAXITS=50;
  const double GOLD=1.618034;
  double Xag[3], r1AXa[3], r1[3], r1A[3], sf, sqrtDelta, A2;
  double curv2, A, B, C, Delta, sol=0.0, ng, curv, lambda;
  sf = *sfA;
  its = 0;
 
  callsprojonto++;
  while (!done && its <= MAXITS)
    {
      itsprojonto++;
      ng = calc_norm(gradf);
      for (k1 = 0; k1 < 3; k1++)
	{
	  Xag[k1] = 0.0;
	  for (k2=0; k2 < 3; k2++)
	    {
	      Xag[k1] += Xa[k1][k2]*gradf[k2];
	    }	      
	}
      if ((OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 4) && OprogStatus.tolAngSD > 0.0)
	{
	  curv = 0.0;
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      curv += dr[k1]*Xag[k1];
	    }
	  lambda = min(sf, Sqr(ng)*fabs(costolAngSD/curv));
	  //printf("curv=%.15G lambda= %.15G\n", curv, lambda);
	  sf = lambda;
	}
      if ((OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 4) && OprogStatus.tolSDconstr > 0.0)
	{
	  //dh = 0.0;
	  curv2 = 0.0;
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      //dh += dr[k1]*gradf[k1];
	      for (k2 = 0; k2 < 3; k2++)
		{
		  curv2 += dr[k1]*Xa[k1][k2]*dr[k2];
		}	
	    }
	  curv2 /= 2.0;
	  //dh /= 2.0;
	  lambda = min(sf,sqrt(OprogStatus.tolSDconstr/fabs(curv2)));
	  //printf("boh=%.15G\n",sqrt(OprogStatus.tolSDconstr/fabs(curv2)) );
	  sf = lambda;
	}
      for (kk=0; kk < 3; kk++)
	{
	  r1[kk] = ri[kk] + dr[kk]*sf; 
	  r1A[kk] = r1[kk] - rA[kk];
	}
      //dr1par = sf*calc_norm(dr);
      A=0;
      B=0;
      C=0;
#if 1
      for (k1 = 0; k1 < 3; k1++)
	{
	  //Xag[k1] = 0.0;
	  r1AXa[k1] = 0.0;
	  for (k2=0; k2 < 3; k2++)
	    {
	      //Xag[k1] += Xa[k1][k2]*gradf[k2];
	      r1AXa[k1] += r1A[k2]*Xa[k2][k1]; 
	    }
	} 
#endif
      for (k1=0; k1 < 3; k1++)
	{

#if 0
	  for (k2=0; k2 < 3; k2++)
	    {
	      A += gradf[k1]*Xa[k1][k2]*gradf[k2];
	      B += r1A[k1]*Xa[k1][k2]*gradf[k2];
	      //  printf("riA[%d]=%f Xa[%d][%d]=%f\n", r1A[k1], Xa[k1][k2]);
	      C += r1A[k1]*Xa[k1][k2]*r1A[k2];
	    }
#else
	  A += gradf[k1]*Xag[k1];
	  B += r1AXa[k1]*gradf[k1];
	  C += r1AXa[k1]*r1A[k1];
#endif
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
#if 0
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
      //if (dist > OprogStatus.epsd && fabs(sol)*ng > OprogStatus.tolSDconstr*sf*calc_norm(dr))
      if (OprogStatus.SDmethod==1 || OprogStatus.SDmethod==3)
	{
	  if (dist > OprogStatus.epsd && fabs(sol)*ng > OprogStatus.tolSDconstr*dist/2.0)
	    {
	      sf /= GOLD;
	      its++;
	      continue;
	    }
	}
#if 0
      /* WARNING: APPARENTEMENTE NON CONVIENE RISCALARE sol e sf se stepSDA è abbastanza piccolo,
       * ma verificare che questo è vero */
      else
	{
#if 0
	  factor = min(1.0, OprogStatus.tolSDconstr*dr1par/ng);
	  sol *= factor;
	  sf *= factor;
#endif
	  if (fabs(sol)*ng > OprogStatus.tolSDconstr*dr1par)
	    {
	      sf /= GOLD;
	      its++;
	      continue;
	    }
	}
#endif
     done = 1;
    }
 
  if (!done)
    {
      printf("maximum number of iterations reached in projont! Aborting...\n");
      printf("sol=%.15G norm(dr)=%.15G sf=%.15G\n", sol, calc_norm(dr), sf);
      exit(-1);
    }
  for (kk = 0; kk < 3; kk++)
    {
      dr[kk] = sol*gradf[kk] + sf*dr[kk]; 
    }
  /* commentando questa riga il valore di sf usato per rimanere "aderenti" alla superficie
   * non viene mantenuto.
   * In tal modo il passo non puo' decrescere in maniera irreversibile se non intorno al minimo. */
  if (OprogStatus.SDmethod == 1 || OprogStatus.SDmethod==3)
    *sfA = sf;
}
void projectgrad(double *p, double *xi, double *gradf, double *gradg)
{
  int kk;
  double dist;
  dist = 0;
  for (kk=0; kk < 3; kk++)
    {
      dist+=Sqr(p[kk+3]-p[kk]);
    }
  dist = sqrt(dist);
  projonto(p, xi, rA, Xa, gradf, &sfA, dist);
  projonto(&p[3], &xi[3], rB, Xb, gradg, &sfB, dist);
}
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
double funcPowellRyck(double phi[])
{
  int kk, k1, k2;
  double sinwA, sinwB, coswA, coswB;
  double MA[3][3], pn[6], MB[3][3], A, B;
  double vA[3], vB[3], gx2[3], fx2[3];
  double F, S;
  sinwA = sin(phi[0]);
  coswA = (1.0 - cos(phi[0]));
  sinwB = sin(phi[1]);
  coswB = (1.0 - cos(phi[1]));

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
       vA[k1] = pcomI[k1]-rA[k1];
       vB[k1] = pcomI[k1+3]-rB[k1];
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
  //printf("scal prod1=%.15G scalprod2=%.15G\n", scalProd(pn,omA), scalProd(&pn[3],omB));
  //printf("pcomI = (%.15G,%.15G,%.15G)\n", pcomI[0], pcomI[1], pcomI[2]);
  for (k1 = 0; k1 < 3; k1++)
    {
      pn[k1] += rA[k1];
      pn[k1+3] += rB[k1];
    }	
   calc_intersec(pn, rA, Xa, pcom2);
   calc_intersec(&pn[3], rB, Xb, &pcom2[3]);
   if (OprogStatus.forceguess)
     {
       for (k1 = 0; k1 < 3; k1++)
 	 {
 	   gx2[k1] = 0;
 	   for (k2 = 0; k2 < 3; k2++)
 	     gx2[k1] += 2.0*Xb[k1][k2]*(pcom2[k2] - rB[k2]);
 	 }
       B = 0.0;
       for (k1 = 0; k1 < 3; k1++)
 	 {
 	   B += (pcom2[k1]-rB[k1])*gx2[k1];
 	 }
       B = 0.5*B - 1.0;
       for (k1 = 0; k1 < 3; k1++)
 	 {
  	   fx2[k1] = 0;
 	   for (k2 = 0; k2 < 3; k2++)
 	     fx2[k1] += 2.0*Xa[k1][k2]*(pcom2[k2+3] - rA[k2]);
 	 }
       A = 0.0;
       for (k1 = 0; k1 < 3; k1++)
 	 {
 	   A += (pcom2[k1+3]-rA[k1])*fx2[k1];
 	 }
       A = 0.5*A - 1.0;
       
       if (A<0 && B<0)
	 S = - OprogStatus.springkSD;
       else
	 S = OprogStatus.springkSD;
     }
  else
    S = OprogStatus.springkSD;
 
  F = 0.0;
  for (kk=0; kk < 3; kk++)
    F += S*Sqr(pcom2[kk]-pcom2[kk+3]);
  //printf("S=%f vec: %f %f dist: %.15G\n",S, phi[0], phi[1], sqrt(fabs(F)));
  return F;
}
double powdirsIRyck[2][2]={{1,0},{0,1}};

void powellmethodRyck(double *fret)
{
  int k1, k2, iter;
  double angs[2]={0,0};
  for (k1=0; k1 < 2; k1++)
    for (k2=0; k2 < 2; k2++)
      powdirs[k1][k2] = powdirsIRyck[k1][k2];
  powell(angs, powdirs, 2, OprogStatus.tolSD, &iter, fret, funcPowellRyck);
}
void linminRyck(double p[], double xi[], double *fret)
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
 * resets p to where the function func(p) takes on a minimum along the direction xi from p,
 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
 * the value of func at the returned location p. This is actually all accomplished by calling
 * the routines mnbrak and brent. */
{ 
  double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
  double f1dimPhi(double x); 
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, 
	      double (*func)(double));
  int j, kk; 
  double r1[3], r2[3];
  double nA, nB;
  double omA[3], omB[3]; 
#if 1
  for (j=0;j<6;j++)
    { 
      pcomI[j]=pcom[j]=p[j];
      xicom[j]=xi[j];
    }
#endif
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
  //printf("pcom ini (%f %f %f)\n", pcom[0], pcom[1], pcom[2]);
  powellmethodRyck(fret);
  //printf("pcom2 fine (%f %f %f)\n", pcom2[0], pcom2[1], pcom2[2]);
  for (j=0;j<6;j++)
    p[j] = pcom2[j]; 
  //printf("xminA=%f xminB=%f dist=%.15G\n", xminA, xminB, *fret);
}

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
void updateByRot(double p[], double xi[])
{
  int kk, k1, k2;
  double sinwA, sinwB, coswA, coswB, phiA, phiB;
  double MA[3][3], pn[6], MB[3][3];
  double omA[3], omB[3], vA[3], vB[3], r1[3], r2[3], nA, nB;
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
double  cgfuncRyck(double *vec);
double zbrentRyck(double (*func)(double), double x1, double x2, double tol);
double get_sign(double *vec);
int check_done(double fp, double fpold, double minax)
{
  const double EPSFR=1E-10;
  double dold=0.0, d=0.0;
  if (OprogStatus.SDmethod != 2 && OprogStatus.SDmethod != 4)
    {
      dold = sqrt(fpold/OprogStatus.springkSD);
      d = sqrt(fp/OprogStatus.springkSD);
    }
  if (OprogStatus.tolSDgrad > 0)
    {
      if (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod==4)
	{
	  if (OprogStatus.epsdSD > 0.0 && fp < Sqr(OprogStatus.epsdSD))
	    return 1;
	  if (doneryck == 1) 
	      //|| 2.0*fabs(dold-d) < OprogStatus.tolSDlong*(fabs(dold)+fabs(d)+EPSFR))
	    {
	      accngA++;
	      accngB++;
	      return 1;
	    }
	}
      else
	{
	  if (fp > Sqr(OprogStatus.epsdSD))//Sqr(OprogStatus.epsd)) 
	    {
	      if (doneryck == 1 || 
		  2.0*fabs(dold-d) < OprogStatus.tolSDlong*(fabs(dold)+fabs(d)+EPSFR))
		{
		  accngA++;
		  accngB++;
		  return 1;
		}
	    }
	  else 
	    {
	      if (2.0*fabs(dold-d) < OprogStatus.tolSD*(fabs(dold)+fabs(d)+EPSFR))
		return 1;
	    }
	}
    }
  else
    {
      if (fp < Sqr(OprogStatus.epsdSD))//Sqr(OprogStatus.epsd))
	{
	  if (2.0*fabs(dold-d) < OprogStatus.tolSD*(fabs(dold)+fabs(d)+EPSFR))
	    return 1;
	}
      else
	{
	  if (2.0*fabs(dold-d) < OprogStatus.tolSDlong*(fabs(dold)+fabs(d)+EPSFR))
	    return 1;
	}
    }
  return 0;
}

void frprmnRyck(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), double (*dfunc)(double [], double [], double [], double [], double*, double*))
  /*Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func,
   * using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is
   * input as ftol. Returned quantities are p (the location of the minimum), iter
   * (the number of iterations that were performed), and fret (the minimum value of the function).
   * The routine linmin is called to perform line minimizations. */
{ 
  int j,its;
  const int ITMAXFR = OprogStatus.maxitsSD;
  const double GOLD=1.618034;
  double fp, fpold=0.0, signA, signB;
  double minax, xi[6], xiold[6];
  double signAold, signBold, pold[6];
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
#if 1
  if ((OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 4) &&
      check_done(fp, fpold, minax))
    {
      callsok++;
      return;
    }
#endif
  projectgrad(p,xi,gradfG,gradgG);  
  for (its=1;its<=ITMAXFR;its++)
    { 
      itsfrprmn++;      
      *iter=its;
      for (j=0; j < n; j++)
	{
	  pold[j] = p[j];
	  xiold[j] = xi[j];
	  p[j] += xi[j];
	}
      signAold = signA;
      signBold = signB;
      fpold = fp; 
      fp = (*dfunc)(p,xi,gradfG, gradgG, &signA, &signB);

      if ((OprogStatus.SDmethod == 1 || OprogStatus.SDmethod == 3) && fp > fpold)
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
    } 
  return; 
  nrerror("Too many iterations in frprmn");
  
}

void frprmn(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double []))
  /*Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func,
   * using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is
   * input as ftol. Returned quantities are p (the location of the minimum), iter
   * (the number of iterations that were performed), and fret (the minimum value of the function).
   * The routine linmin is called to perform line minimizations. */
{ 
  void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),
	       void (*dfunc)(double [], double [])); 
  int j,its;
  const int ITMAXFR = OprogStatus.maxitsSD;
  const double EPSFR=1E-10;
  double gg,gam,fp,dgg;
  double g[8],h[8],xi[8];
  fp=(*func)(p); /*Initializations.*/
  (*dfunc)(p,xi); 
  callsfrprmn++;
  // printf("P p=%.15G %.15G %.15G %.15G %.15G %.15G \n",
  //	     p[0], p[1], p[2], p[3], p[4], p[5]);
  // printf("P xi=%.15G %.15G %.15G %.15G %.15G %.15G \n",
  //	 xi[0], xi[1], xi[2], xi[3], xi[4], xi[5]);

  for (j=0;j<n;j++)
    { 
      g[j] = -xi[j]; 
      xi[j]=h[j]=g[j];
    }
  for (its=1;its<=ITMAXFR;its++)
    { /* Loop over iterations.*/
       itsfrprmn++;      
      *iter=its;
      linmin(p,xi,n,fret,func); /* Next statement is the normal return: */
      //printf("its=%d 2.0*fabs(*fret-fp):%.15G rs: %.15G fp=%.15G fret: %.15G\n",its, 2.0*fabs(*fret-fp),ftol*(fabs(*fret)+fabs(fp)+EPSFR),fp,*fret );
      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPSFR)) 
	{ 
	  callsok++;
	  return;
	} 
      fp= *fret; 
      (*dfunc)(p,xi);
      //printf("p=%.15G %.15G %.15G %.15G %.15G %.15G \n",
      //     p[0], p[1], p[2], p[3], p[4], p[5]);
      //printf("xi=%.15G %.15G %.15G %.15G %.15G %.15G \n",
      // xi[0], xi[1], xi[2], xi[3], xi[4], xi[5]);
      dgg=gg=0.0;
      for (j=0;j<n;j++) 
	{
	  gg += g[j]*g[j]; 
	  //dgg += xi[j]*xi[j];  /* This statement for Fletcher-Reeves.*/
	  dgg += (xi[j]+g[j])*xi[j]; /*This statement for Polak-Ribiere.*/
	} 
      if (gg == 0.0) 
	{ /* Unlikely. If gradient is exactly zero then we are already done.*/
	  callsok++;
	  return; 
	} 
      gam=dgg/gg; 
      for (j=0;j<n;j++) 
	{ 
	  g[j] = -xi[j]; xi[j]=h[j]=g[j]+gam*h[j]; 
	} 
    } 
  return;
  nrerror("[frprmn]Too many iterations in frprmn");
}

/* =========================== >>> forces <<< ======================= */
double  cgfunc2(double *vec)
{
  int k1, k2;
  double fx[3], gx[3], fx2[3], gx2[3];
  double Q1, Q2, A, B, F;
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
  Q1 = 0.0;
  Q2 = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      Q1 += (vec[k1]-rA[k1])*fx[k1];
      Q2 += (vec[k1+3]-rB[k1])*gx[k1];
    }
  Q1 = 0.5*Q1 - 1.0;
  Q2 = 0.5*Q2 - 1.0;
  F = A+B;
  F += OprogStatus.springkSD*Sqr(Q1);
  F += OprogStatus.springkSD*Sqr(Q2);
  return F;
}

void gradcgfunc2(double *vec, double *grad)
{
  int kk, k1, k2; 
  double fx[3], gx[3], fx2[3], gx2[3];
  double Q1, Q2;
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
  Q1 = 0.0;
  Q2 = 0.0;
 
  for (k1 = 0; k1 < 3; k1++)
    {
      Q1 += (vec[k1]-rA[k1])*fx[k1];
      Q2 += (vec[k1+3]-rB[k1])*gx[k1];
    }
  Q1 = 0.5*Q1 - 1.0;
  Q2 = 0.5*Q2 - 1.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      gx2[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx2[k1] += 2.0*Xb[k1][k2]*(vec[k2] - rB[k2]);
      fx2[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx2[k1] += 2.0*Xa[k1][k2]*(vec[k2+3] - rA[k2]);
    }
#if 0
  A = B = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      A += (vec[k1+3]-rA[k1])*fx2[k1];
      B += (vec[k1]-rB[k1])*gx2[k1];
    }
  A = 0.5*A - 1.0;
  B = 0.5*B - 1.0;
#endif
  for (kk=0; kk < 3; kk++)
    {
      grad[kk] = gx2[kk];
      grad[kk] += 2.0*OprogStatus.springkSD*fx[kk]*Q1;
      grad[kk+3] = fx2[kk];
      grad[kk+3] += 2.0*OprogStatus.springkSD*gx[kk]*Q2;
    }
}
void gradcgfunc(double *vec, double *grad)
{
  int kk, k1, k2; 
  double fx[3], gx[3], fx2[3], gx2[3];
  double S=1.0, Q1, Q2, A, B;
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
  Q1 = 0.0;
  Q2 = 0.0;
  A = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      Q1 += (vec[k1]-rA[k1])*fx[k1];
      Q2 += (vec[k1+3]-rB[k1])*gx[k1];
    }
  Q1 = 0.5*Q1 - 1.0;
  Q2 = 0.5*Q2 - 1.0;
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
	S = - OprogStatus.springkSD;
      else
	S = OprogStatus.springkSD;
    }
  else
    S= OprogStatus.springkSD;
  for (kk=0; kk < 3; kk++)
    {
      //grad[kk]=-2.0*(vec[kk+3]-vec[kk])*A + vec[6]*fx[kk];
      //grad[kk+3]=2.0*(vec[kk+3]-vec[kk])*A + vec[7]*gx[kk];
      grad[kk]= -2.0*(vec[kk+3]-vec[kk])*S;
      
      if (cghalfspring && Q1 > 0)
	grad[kk] +=  2.0*lambdacg*fx[kk]*Q1;
      //grad[kk] += vec[6]*fx[kk];

      grad[kk+3]= 2.0*(vec[kk+3]-vec[kk])*S;
      if (cghalfspring && Q2 > 0)
	grad[kk+3] +=  2.0*lambdacg*gx[kk]*Q2;
      //grad[kk+3] += vec[7]*gx[kk];
    }
  //grad[6] = Sqr(Q1);
  //grad[7] = Sqr(Q2);
  //grad[6] = -Q1;
 // grad[7] = -Q2;
}
/* =========================== >>> forces <<< ======================= */
double  cgfunc(double *vec)
{
  int kk, k1, k2;
  double fx[3], gx[3], fx2[3], gx2[3];
  double Q1, Q2, A, B, F;
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
	A = - OprogStatus.springkSD;
      else
	A = OprogStatus.springkSD;
    }
  else
    A= OprogStatus.springkSD;
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
  Q1 = 0.0;
  Q2 = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      Q1 += (vec[k1]-rA[k1])*fx[k1];
      Q2 += (vec[k1+3]-rB[k1])*gx[k1];
    }
  Q1 = 0.5*Q1 - 1.0;
  Q2 = 0.5*Q2 - 1.0;
  F = 0.0;
  for (kk=0; kk < 3; kk++)
    F += A*Sqr(vec[kk]-vec[kk+3]);
#if 1
  if (cghalfspring && Q1 > 0)
    F += lambdacg*Sqr(Q1);
  if (cghalfspring && Q2 > 0)
    F += lambdacg*Sqr(Q2);
#endif
  //F += vec[6]*Q1;
  //F += vec[7]*Q2;
  //printf("A=%f vec: %f %f %f, %f %f %f Epoten: %.15G\n", A,vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], F);
  //printf("Q1=%.15G Q2=%.15G\n",  Q1, Q2);
  return F;
}
extern double costolSDgrad;
double gradcgfuncRyck(double *vec, double *grad, double *fx, double *gx, double *signA, double *signB)
{
  int kk, k1, k2; 
  double K1, K2, F, nf, ng, gx2[3], fx2[3], dd[3], normdd, ngA, ngB;
  double S=1.0, A=1.0, B, gradfx, gradgx;
  doneryck = 0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += Xa[k1][k2]*(vec[k2] - rA[k2]);
      fx[k1] *= 2.0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += Xb[k1][k2]*(vec[k2+3] - rB[k2]);
      gx[k1] *= 2.0;
      dd[k1] = vec[k1+3]-vec[k1];
    }

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
      *signA = A = 0.5*A - 1.0;
      *signB = B = 0.5*B - 1.0;
      
      if (A<0 && B<0)
	S = -1.0;
      else
	S = 1.0;
    }
  normdd = calc_norm(dd);
  if (normdd==0)
    {
      doneryck = 2;
      return 0;
    }
  /* la norma dei gradienti è sempre stepSDA e stepSDB*/ 
  if (OprogStatus.SDmethod==1 || OprogStatus.SDmethod==3)
    {
      K1= icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
      K2= jcg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
    }
  else
    {
      K1 = OprogStatus.springkSD;
      K2 = OprogStatus.springkSD;
    }
  for (kk=0; kk < 3; kk++)
    {
      if (OprogStatus.SDmethod == 1 || OprogStatus.SDmethod==3)
	grad[kk] = S*dd[kk]/normdd;
      else
	grad[kk] = S*dd[kk];
      grad[kk+3]= -K1*grad[kk];
      grad[kk] *= K2;
    }
  nf = calc_norm(fx);
  ng = calc_norm(gx);
  for (k1=0; k1 < 3; k1++)
    {
      fx[k1] /= nf;
      gx[k1] /= ng;
    }
  gradfx = 0;
  gradgx = 0;
  for (k1=0; k1 < 3; k1++)
    {
      gradfx += grad[k1]*fx[k1]; 
      gradgx += grad[k1+3]*gx[k1];
    }
  for (kk=0; kk < 3; kk++)
    {
      grad[kk] -= gradfx*fx[kk];
      grad[kk+3] -= gradgx*gx[kk];
    }
  if (OprogStatus.tolSDgrad > 0.0)
    {
      if (OprogStatus.SDmethod==1 || OprogStatus.SDmethod==3)
	{
	  ngA = ngB = 0;  
	  for (kk=0; kk < 3; kk++)
	    {
	      ngA += Sqr(grad[kk]);
	      ngB += Sqr(grad[kk+3]); 
	      
	    }
	  if (sqrt(ngA) < OprogStatus.tolSDgrad*(icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB)
	      && sqrt(ngB) < OprogStatus.tolSDgrad*(jcg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB))
	    {
	      //accngA++;
	      //accngB++;
	      doneryck = 1;
	    }
	}
      else
	{
	  if (fabs(scalProd(fx,dd) > normdd*costolSDgrad && fabs(scalProd(gx,dd)) > normdd*costolSDgrad))
	    {
	      //accngA++;
	      //accngB++;
	      doneryck = 1;
	    }
	}
    }
  S *= OprogStatus.springkSD;
  F = S*Sqr(normdd);
  return F; 
}
/* =========================== >>> forces <<< ======================= */
double  cgfuncRyck(double *vec)
{
  int kk, k1, k2;
  double fx2[3], gx2[3];
  double A, B, F;
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
	A = - OprogStatus.springkSD;
      else
	A = OprogStatus.springkSD;
    }
  else
    A = OprogStatus.springkSD;
 
  F = 0.0;
  for (kk=0; kk < 3; kk++)
    F += A*Sqr(vec[kk]-vec[kk+3]);
  return F;
}
double max3(double a, double b, double c)
{
  double m;
  m = a;
  if (b > m)
    m = b;
  if (c > m)
    m = c;
  return m;
}

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
  double distSqold, distSqT, distSqmin=1E200;
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
void guessdistByMesh(int i, int j, double shift[3], double *vecg)
{
  int kk;
  double vec[8];
  icg = i;
  jcg = j;
#ifdef MD_POLYDISP
  minaxicg = min3(axaP[i],axbP[i],axcP[i]);
  minaxjcg = min3(axaP[j],axbP[j],axcP[j]);
#else
  if (i < Oparams.parnumA)
    minaxicg = min3(Oparams.a[0],Oparams.b[0],Oparams.c[0]);
  else 
    minaxicg = min3(Oparams.a[1],Oparams.b[1],Oparams.c[1]);
  if (j < Oparams.parnumA)
    minaxjcg = min3(Oparams.a[0],Oparams.b[0],Oparams.c[0]);
  else 
    minaxjcg = min3(Oparams.a[1],Oparams.b[1],Oparams.c[1]);
#endif
  for (kk=0; kk < 3; kk++)
    {
      shiftcg[kk] = shift[kk];
    }
  for (kk=0; kk < 6; kk++)
    {
      vec[kk] = vecg[kk];
    }
 
  findminMesh(vec);
  //frprmnRyck(vec, 6, OprogStatus.tolSD, &iter, &Fret, cgfuncRyck, gradcgfuncRyck);
  for (kk=0; kk < 6; kk++)
    {
      vecg[kk] = vec[kk];
    }
}
double conjgrad_func(double angs[4]);

double f1dim(double x) 
  /*Must accompany linmin.*/
{
  int j; 
  double f, xt[8];
  for (j=0;j<ncom;j++) 
    xt[j]=pcom[j]+x*xicom[j]; 
  f=(*nrfunc)(xt); 
  return f;
}
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
/* Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
 * resets p to where the function func(p) takes on a minimum along the direction xi from p,
 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
 * the value of func at the returned location p. This is actually all accomplished by calling
 * the routines mnbrak and brent. */
{ 
  const double TOLLM=1.0E-10;
  int j; 
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n; /*Define the global variables.*/
  nrfunc=func; 
  for (j=0;j<n;j++)
    { 
      pcom[j]=p[j];
      xicom[j]=xi[j];
    } 
  ax=0.0; /*Initial guess for brackets.*/
  xx=1.0; 
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim); 
  *fret=brent(ax,xx,bx,f1dim,TOLLM,&xmin);
  for (j=0;j<n;j++)
    { 
      /*Construct the vector results to return. */
      xi[j] *= xmin;
      p[j] += xi[j]; 
    } 
}
#ifdef MD_CONJGRAD
double conjgrad_func(double angs[4])
{
  double dist, r1[3], r2[3], xp[6];
  angs2coord(angs, xp);
  body2lab(icg, xp, r1, rA, RtA);
  body2lab(jcg, &xp[3], r2, rB, RtB);
  printf("angs=(%.15G,%.15G,%.15G,%.15G)\n", angs[0], angs[1], angs[2], angs[3]);
  dist = Sqr(r1[0]-r2[0])+Sqr(r1[1]-r2[1])+Sqr(r1[2]-r2[2]); 
  return dist;
}
void conjgrad_grad(double *angs, double *grad)
{
  double r1[3], r2[3], xp[6], dd[3], jac[6][4];
  double sin0, sin1, cos0, cos1, sin2, sin3, cos2, cos3;
  int kk, k1, k2;

  sin0 = sin(angs[0]);
  sin1 = sin(angs[1]);
  cos0 = cos(angs[0]);
  cos1 = cos(angs[1]);
  sin2 = sin(angs[2]);
  sin3 = sin(angs[3]);
  cos2 = cos(angs[2]);
  cos3 = cos(angs[3]);
  xp[0] = axa[icg]*cos0*sin1;
  xp[1] = axb[icg]*sin0*sin1;
  xp[2] = axc[icg]*cos1;
  xp[3] = axa[jcg]*cos2*sin3;
  xp[4] = axb[jcg]*sin2*sin3;
  xp[5] = axc[jcg]*cos3;

  body2lab(icg, xp, r1, rA, RtA);
  body2lab(jcg, &xp[3], r2, rB, RtB);

  for (kk = 0; kk < 3; kk++)
    dd[kk] = -2.0*(r2[kk] - r1[kk]);   
  printf("[CG grad] i=%d j=%d dist=%.15G\n", icg, jcg, calc_norm(dd)/2.0);
  /* NOTA: jac[][] è lo jacobiano del cambio di coordinate 
   * (theta1,phi1,theta2,phi2)->(x1,y1,z1,x2,y2,z2) cioè 
   * \frac{\delta(x1,y1,z1,x2,y2,z2)}{\delta(theta1,phi1,theta2,phi2)}*/
  jac[0][0] = -axa[icg]*sin0*sin1;
  jac[0][1] = axa[icg]*cos0*cos1;
  jac[0][2] = 0.0;
  jac[0][3] = 0.0;
  jac[1][0] = axb[icg]*cos0*sin1;
  jac[1][1] = axb[icg]*sin0*cos1; 
  jac[1][2] = 0.0;
  jac[1][3] = 0.0;
  jac[2][0] = 0.0;
  jac[2][1] = -axc[icg]*sin1; 
  jac[2][2] = 0.0;
  jac[2][3] = 0.0;

  jac[3][0] = 0.0;
  jac[3][1] = 0.0;
  jac[3][2] = -axa[jcg]*sin2*sin3;
  jac[3][3] = axa[jcg]*cos2*cos3;
  jac[4][0] = 0.0;
  jac[4][1] = 0.0;
  jac[4][2] = axb[jcg]*cos2*sin3;
  jac[4][3] = axb[jcg]*sin2*cos3;
  jac[5][0] = 0.0;
  jac[5][1] = 0.0;
  jac[5][2] = 0.0;
  jac[5][3] = -axc[jcg]*sin3;

  for (kk = 0; kk < 2; kk++)
    {
      grad[kk] = 0.0;
      grad[kk+2] = 0.0;
      for (k1 = 0; k1 < 3; k1++)
	{
	  grad[kk] += -dd[k1]*jac[k1][kk];
	  grad[kk+2] += dd[k1]*jac[k1+3][kk+2];
	}
    }
  if (angs[1]+grad[1] > pi)
    grad[1] = 2.0*pi - grad[1];  
  if (angs[1]+grad[1] < 0.0)
    grad[1] = -grad[1];
  if (angs[3]+grad[3] > pi)
    grad[3] = 2.0*pi - grad[3];  
  if (angs[3]+grad[3] < 0.0)
    grad[3] = -grad[3];
   //printf("grad=(%.15G,%.15G,%.15G,%.15G\n", grad[0], grad[1], grad[2], grad[3]);
}
double rDcg[3];
int conjgrad(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double []))
{
  const int ITMAX=10;
  const double EPS = 1.0E-12;
  int j,its;
  double gg,gam,fp,dgg, g[8],h[8],xi[8], dist, distold;
  /* Initializations.*/
  fp=(*func)(p); 
  (*dfunc)(p,xi);
  //distold = Sqr(rDcg[0]-p[0])+Sqr(rDcg[1]-p[1])+Sqr(rDcg[2]-p[2]);

  for (j=0; j < n; j++) 
    {
      g[j] = -xi[j];
      xi[j] = h[j] = g[j];
    }
  for (its = 0; its < ITMAX; its++) 
    { 
      *iter=its;
      linmin(p, xi, n, fret, func); /* Next statement is the normal return: */
      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) 
	{
	  return 0;
	}
      fp= *fret;
      (*dfunc)(p,xi);
      dgg = gg = 0.0;
      for (j = 0; j < n; j++)
	{
	  gg += g[j]*g[j];
	  /* dgg += xi[j]*xi[j]; */ /* This statement for Fletcher-Reeves.*/
	  dgg += (xi[j]+g[j])*xi[j]; /* This statement for Polak-Ribiere.*/
	}
      if (gg == 0.0) 
	{ 
	  /* Unlikely. If gradient is exactly zero then
	     we are already done. */
	  return 0;
	}
      gam=dgg/gg;
      for (j = 0; j < n; j++) 
	{
	  g[j] = -xi[j];
	  xi[j]=h[j]=g[j]+gam*h[j];
	}
    }
  //printf("ERROR: Too many iterations in frprmn *fret=%.15G\n", *fret);
  //exit(-1);
  return 0;
}
extern double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
double conjgrad_dist5(double x[])
{
  int k1, k2; 
  double fx[3], gx[3], rD[3];
  /* x = (r, alpha, t) */ 
  double val, dtmp1, dtmp2;
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
	}
      rDcg[k1] = rD[k1] = x[k1] + fx[k1]*x[4];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(rD[k2] - rB[k2]);
    }
  val = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      dtmp1 = fx[k1] + Sqr(x[3])*gx[k1];
      val += Sqr(dtmp1);
    }
  dtmp1 = 0.0;
  dtmp2 = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      dtmp1 += (x[k1]-rA[k1])*fx[k1];
      dtmp2 += (rD[k1]-rB[k1])*gx[k1];
    }
  dtmp1 = 0.5*dtmp1-1.0;
  dtmp2 = 0.5*dtmp2-1.0;
  val += Sqr(dtmp1)+Sqr(dtmp2);
  //printf("val=%.15G dtmp1:%.15G dtmp2=%.15G x=%.15G %.15G %.15G %.15G %.15G\n", val,
//	 dtmp1, dtmp2,
//	 x[0],x[1],x[2],x[3],x[4]);
  //printf("val=%.15G\n", val);
  //printf("BAHBOH dist=%.15G\n", sqrt(Sqr(x[0]-rD[0])+Sqr(x[1]-rD[1])+Sqr(x[2]-rD[2])));
  return val;
}
void conjgrad_graddist5(double *x, double *xi)
{
  double fx[3], xit[5], gx[3], rD[3], A[3][3], b[3], c[3], f, g;
  int k1, k2, k3;
  double df[5][5], fgx[3];
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  A[k1][k2] = XbXa[k1][k2];
	  A[k1][k2] *= 4.0*Sqr(x[3])*x[4];
	  A[k1][k2] += 2.0*Xb[k1][k2]*Sqr(x[3]);
	}
    }	
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2] + A[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	}
      rD[k1] = x[k1] + fx[k1]*x[4];
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  gx[k1] += 2.0*Xb[k1][k2]*(rD[k2]-rB[k2]);
	}
    } 
  f = 0.0;
  g = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      f += (x[k1]-rA[k1])*fx[k1];
      g += (rD[k1]-rB[k1])*gx[k1];
    }
  f = 0.5*f - 1.0;
  g = 0.5*g - 1.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      b[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  b[k1] += Xb[k1][k2]*fx[k2];
	}
      b[k1] *= 2.0*Sqr(x[3]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      c[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  c[k1] += gx[k2]*Xa[k2][k1];
	}
      c[k1] *= 2.0*x[4];
      c[k1] += gx[k1];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  df[3][3] = 0.0;
  df[3][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = c[k1];
    } 
  df[4][3] = 0.0;
  df[4][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    df[4][4] += gx[k1]*fx[k1];
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = 2.0*x[3]*gx[k1];
      df[k1][4] = b[k1];
    }
  
  for (k1 = 0; k1 < 3; k1++)
    {
      xit[k1] = fx[k1]+Sqr(x[3])*gx[k1];
      xit[k1] *= 2.0;
    }
  xit[3] = 2.0*f;
  xit[4] = 2.0*g;

  for (k1 = 0; k1 < 5; k1++)
    {
      xi[k1] = 0.0;
      for (k2 = 0; k2 < 5; k2++)
	xi[k1] += xit[k2]*df[k2][k1];
    }
}
double conjgrad_dist(double *vec)
{
  return 0.0;
}
void conjgrad_graddist(double *vec, double *xi)
{

}
void preNR_conjgrad(int i, int j, double *vec)
{
  int iter;
  double Fret;
  icg = i;
  jcg = j;
  //printf("========= >>> INIZIO CONJ GRAD <<< =======\n");
  if (OprogStatus.dist5)
    conjgrad(vec, 5, 1E-12, &iter, &Fret, conjgrad_dist5, conjgrad_graddist5);
  else
    conjgrad(vec, 8, 1E-14, &iter, &Fret, conjgrad_dist, conjgrad_graddist);
  //printf("[preNR CG] iter=%d\n", iter);
}


void distconjgrad(int i, int j, double shift[3], double *vec)
{
  int kk, iter;
  double Fret, vecP[6], angs[4], sinth;
  icg = i;
  jcg = j;
  for (kk=0; kk < 3; kk++)
    {
      shiftcg[kk] = shift[kk];
    }
  /* trovo le coordinate nel riferimento del corpo rigido (assi principali) */
  lab2body(icg, vec, vecP, rA, RtA);
  lab2body(jcg, &vec[3], &vecP[3], rB, RtB);
  /* determino theta e phi per entrambi gli ellissoidi 
   * angs[] = (thetaA,phiA,thetaB,phiB)
   * 0 < theta < 2PI
   * 0 < phi < PI */
  angs[1] = acos(vecP[2]/axc[icg]);
  angs[0] = acos(vecP[0]/axa[icg]/sin(angs[1]));
  sinth = vecP[1]/axb[icg]/sin(angs[1]);
  if (sinth < 0.0)
    angs[0] = 2.0*pi-angs[0];
  angs[3] = acos(vecP[5]/axc[jcg]);
  angs[2] = acos(vecP[3]/axa[jcg]/sin(angs[3]));
  sinth = vecP[4]/axb[jcg]/sin(angs[3]);
  if (sinth < 0)
    angs[2] = 2.0*pi-angs[2];
  printf("[CG] PRIMA dist: %.15G\n", sqrt(conjgrad_func(angs)));
  printf("[CG] PRIMA angs=%f %f %f %f\n", angs[0], angs[1], angs[2], angs[3]);
  conjgrad(angs, 4, 1E-14, &iter, &Fret, conjgrad_func, conjgrad_grad);
  angs2coord(angs, vecP);
  printf("[CG] DOPO dist: %.15G iter=%d\n", sqrt(conjgrad_func(angs)), iter);
  /* torno nel riferimento del laboratorio */
  body2lab(icg, vecP, vec, rA, RtA);
  body2lab(jcg, &vecP[3], &vec[3], rB, RtB);
}
#endif
void distSD(int i, int j, double shift[3], double *vecg, double lambda, int halfspring)
{
  int kk;
  double Fret;
  int iter;
  double vec[8];
#ifdef EDHE_FLEX
  int typei, typej;
  double axaiF, axbiF, axciF, axajF, axbjF, axcjF;
#endif
  icg = i;
  jcg = j;
#ifdef EDHE_FLEX
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  axaiF = typesArr[typei].sax[0];
  axbiF = typesArr[typei].sax[1];
  axciF = typesArr[typei].sax[2];
  minaxicg = min3(axaiF, axbiF, axciF);
  axajF = typesArr[typej].sax[0];
  axbjF = typesArr[typej].sax[1];
  axcjF = typesArr[typej].sax[2];
  minaxjcg = min3(axajF, axbjF, axcjF);
#elif defined(MD_POLYDISP)
  minaxicg = min3(axaP[i], axbP[i], axcP[i]);
  minaxjcg = min3(axaP[j], axbP[j], axcP[j]);
#else
  if (i < Oparams.parnumA)
    minaxicg = min3(Oparams.a[0],Oparams.b[0],Oparams.c[0]);
  else 
    minaxicg = min3(Oparams.a[1],Oparams.b[1],Oparams.c[1]);
  if (j < Oparams.parnumA)
    minaxjcg = min3(Oparams.a[0],Oparams.b[0],Oparams.c[0]);
  else 
    minaxjcg = min3(Oparams.a[1],Oparams.b[1],Oparams.c[1]);
#endif
  lambdacg = lambda;
  cghalfspring = halfspring;
  for (kk=0; kk < 3; kk++)
    {
      shiftcg[kk] = shift[kk];
    }
  for (kk=0; kk < 6; kk++)
    {
      vec[kk] = vecg[kk];
    }
  
  //estimate_lambda(vec);
#if 0
  printf(">>> vec: %.15G %.15G %.15G %.15G %.15G %.15G\n", vec[0], vec[1], vec[2],
  	 vec[3], vec[4], vec[5]);
  printf(">>> vec[6]:%.15G vec[7]: %.15G\n", vec[6], vec[7]);
#endif
  //frprmn(vec, 6, OprogStatus.tolSD, &iter, &Fret, cgfunc2, gradcgfunc2);
  //powellmethodPenalty(vec);
  frprmnRyck(vec, 6, OprogStatus.tolSD, &iter, &Fret, cgfuncRyck, gradcgfuncRyck);
  //powellmethod(vec);
  //findminMesh(vec);
  for (kk=0; kk < 6; kk++)
    {
      vecg[kk] = vec[kk];
    }
}

#include <math.h>
#define ZBRAC_FACTOR 1.6
#define ZBRAC_NTRY 50
int zbrac(double (*func)(double), double *x1, double *x2)
{
  int j;
  double f1,f2;
  if (*x1 == *x2) 
    {
      printf("Bad initial range in zbrac");
      return 1;
    }
  f1=(*func)(*x1);
  f2=(*func)(*x2);
  for (j=1;j<=ZBRAC_NTRY;j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2))
      f1=(*func)(*x1 += ZBRAC_FACTOR*(*x1-*x2));
    else
      f2=(*func)(*x2 += ZBRAC_FACTOR*(*x2-*x1));
  }
  return 0;
}

#define ITMAXZB 100 
/* Maximum allowed number of iterations.*/
#define EPSP 1.11e-16 /* Machine floating-point precision (see http://en.wikipedia.org/wiki/Machine_epsilon).*/
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
extern int ibr;
double zbrent(double (*func)(double), double x1, double x2, double tol)
/* Using Brent s method, find the root of a function func known to lie between x1 and x2. 
 * The root, returned as zbrent, will be refined until its accuracy is tol.*/
{
  int iter; 
  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2; 
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm; 
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    {
      MD_DEBUG(printf("BRENT BAD BRACKETING fa(%.15G)=%.15G fb(%.15G)=%.15G\n", a, fa, b, fb));
      polinterr = 1;
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
      if (polinterr)
	return 0.0;
    } 

  printf("[zbrent] BRENT TOO MANY ITERATIONS\n");
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
#ifndef MC_SIMUL
	  printf("ERROR: Singular matrix in routine ludcmp\n"); 
#endif
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
  int i, j, c1, c2, pivot[MD_NBMAX];
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
#if 0
#define TOLX2 1.E-6
#define TOLF2 1.0E-3
#endif
#define ALF 1.0e-4 /* Ensures sufficient decrease in function value.*/
#ifdef MD_NEW_NR_CHECKS
#ifdef MD_SUPERELLIPSOID
#define NR_DAMP_FACT 1E-2
#define TOLX 1.0E-14//1.0e-7 /* Convergence criterion on  x.*/ 
#define TOLXD 1.0E-14
/* NOTA 13/04/2010: TOLXDNL sembra essere l'unico parametro critico
   per la convergenza del NR nel caso d'urto fra una SQ ed un piano 
   infatti accade spesso che la convergenza è su x cioè il punto non cambia
   in maniera significativa (questo dovrebbe accadere quando uno spigolo della SQ 
   urta un piano).
 */
#define TOLXDNL 1.0E-14

/* NOTA 13/04/2010: per ora TOLXNL non lo uso */
#define TOLXNL 1.0E-14
#else
#define TOLX 1.0E-14//1.0e-7 /* Convergence criterion on  x.*/ 
#define TOLXD 1.0E-14
#define TOLXDNL TOLXD
#define TOLXNL TOLX
#endif
#else
#define TOLX 1.0E-14//1.0e-7 /* Convergence criterion on  x.*/ 
#define TOLXD 1.0E-14
#define TOLXNL TOLX
#define TOLXDNL TOLXD
#endif
#define MAXITS 100 // se le particelle non si urtano il newton-raphson farà MAXITS iterazioni
#define MAXITS2 100
/* N.B. se dovessero esserci problemi aumentare TOLF e TOLFD */
#ifdef MD_NEW_NR_CHECKS
#ifdef MD_SUPERELLIPSOID
/* NOTA 13/04/2010: per ora TOLFNL e TOLXNL non li uso poiché in newtNeigh 
   sembra che non ci siano problemi ad usare TOLF e TOLX come nell'urto tra SQ. */
#define TOLF 1.0E-11// 1.0e-4
#define TOLFD 1.0E-11
#define TOLFNL 1.0E-11
/* NOTA 13/04/2010: usando TOLFDNL come per l'urto tra SQ non ci sono problemi */
#define TOLFDNL 1.0E-11
#else
#define TOLF 1.0E-11// 1.0e-4
#define TOLFD 1.0E-11
#define TOLFNL TOLF
#define TOLFDNL TOLFD
#endif
#else
#define TOLF 1.0E-10// 1.0e-4
#define TOLFD 1.0E-10
#define TOLFNL TOLF
#define TOLFDNL TOLFD
#endif
#ifdef MD_SUPERELLIPSOID
#define TOLMIN 1.0E-12
#define TOLMINNL 1.0E-14
#define TOLMINDNL 1.0E-14
#else
#define TOLMIN 1.0E-12//1.0e-6 
#define TOLMINDNL TOLMIN
#define TOLMINNL TOLMIN
#endif
#define STPMX 100.0
void lnsrchNeigh(int n, double xold[], double fold, double g[], double p[], double x[], 
	    double *f, double stpmax, int *check, 
	    double (*func)(double [], int), int iA, 
	    double tolx)
/*
   Given an n-dimensional point xold[1..n], the value of the function and gradient there, 
   fold and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p
   from xold where the function func has decreased  "sufficiently".  The new function value is 
   returned in f. stpmax is an input quantity that limits the length of the steps so that 
   you do not try to evaluate the function in regions where it is undefined or subject 
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
    {
      *check=2;
      return;
     // nrerror("Roundoff problem in lnsrch."); 
    }
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
      *f=(*func)(x,iA); 
      if (alam < alamin) 
	{ /* Convergence on  x. For zero finding, the calling program 
	     should verify the convergence.*/ 
	  for (i=0;i<n;i++) 
	    x[i]=xold[i]; 
	  *check=1; 
	  return;
	}
      else if (*f <= fold+ALF*alam*slope) 
	return; 
	/* Sufficient function decrease.*/
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
    {
      *check=2;
      return;
      //nrerror("Roundoff problem in lnsrch."); 
    }
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
	{ /* Convergence on  x. For zero finding, the calling program 
	     should verify the convergence.*/ 
	  for (i=0;i<n;i++) 
	    x[i]=xold[i]; 
	  *check=1; 
	  return;
	}
      else if (*f <= fold+ALF*alam*slope) 
	return; 
	/* Sufficient function decrease.*/
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
double **fjac,*g,*p,*xold;
int *indx;
#if 0
#define FREERETURN {MD_DEBUG10(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d fvec=(%.15G,%.15G,%.15G,%.15G,%.15G)\n", x[0], x[1], x[2], x[3], x[4], test, its, *check, fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));\
free_vector(fvec);free_vector(xold); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);free_vector(fvecG);return;}
#define FREERETURND {MD_DEBUG(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d\n", x[0], x[1], x[2], x[3], x[4], test, its, *check));\
free_vector(fvecD);free_vector(xold); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);return;}
#else
#define FREERETURND return;
#define FREERETURN  return;
#endif
void (*nrfuncv)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncv2)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncvD)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncvNeigh)(int n, double v[], double f[], int iA); 
void (*nrfuncvDNeigh)(int n, double v[], double f[], int iA);

extern void fdjac(int n, double x[], double fvec[], double **fjac, 
		  void (*vecfunc)(int n, double v[], double fvec[], int i, int j, double shift[3]), int iA, int iB, double shift[3]); 
double fminNeigh(double x[], int iA);
double fminDNeigh(double x[], int iA);
double fminD(double x[], int iA, int iB, double shift[3]);
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], double *f, 
	    double stpmax, int *check, double (*func)(double [], int, int, double []),
	    int iA, int iB, double shift[3], double tolx);
void lubksb(double **a, int n, int *indx, double b[]); 
void ludcmp(double **a, int n, int *indx, double *d, int *ok); 
extern void funcs2beZeroedGuess(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroed(int n, double x[], double fvec[], int i, int j, double shift[3]);

extern void upd2tGuess(int i, int j, double shift[3], double tGuess);
//#define MD_GLOBALNR
#ifdef MD_NEW_NR_CHECKS
double test_func_values(double *fvec, int n)
{
  int i;
  double test=0.0;
  for (i=0;i<n;i++) 
    if (fabs(fvec[i]) > test) 
      test=fabs(fvec[i]);
  return test;
}
double test_xvalues(double *xold, double *x, int n)
{
  double test=0.0;
  int i;
  for (i=0;i<n;i++) 
    {
      temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
      //temp=(fabs(x[i]-xold[i]))/fabs(x[i]); 
      if (temp > test) 
	test=temp; 
    }
  return test;
}
#endif
void newtNeigh(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int),
	  int iA)
{
  int i,its,ok;
  double d,f,stpmax,sum,test, den; 
#ifdef MD_GLOBALNRNL
  int j;
  double fold;
#endif
#if 0
  int *indx;
  double **fjac,*g,*p,*xold;
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvec=vector(n); 
  fvecG=vector(n);
#endif
  /*Define global variables.*/
  nn=n; 
  nn2=n-1;
  nrfuncvNeigh=vecfunc; 
  f=fminNeigh(x,iA); /*fvec is also computed by this call.*/
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
  for (i=0;i<n;i++) 
    if (fabs(fvec[i]) > test)
      test=fabs(fvec[i]); 
  if (test < 0.01*TOLFNL)
    {
      *check=0; 
      FREERETURN;
    }
  for (sum=0.0,i=0;i<n;i++) 
    sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
#ifdef MD_NNLPLANES
  funcs2beZeroedNeighPlane(n,x,fvec,iA);
#else
  funcs2beZeroedNeigh(n,x,fvec,iA);
#endif
  for (its=0;its<MAXITS;its++)
    { /* Start of iteration loop. */
       //funcs2beZeroed(n,x,fvec,iA,iB,shift);
#ifdef MD_NNLPLANES
       fdjacNeighPlane(n,x,fvec,fjac,vecfunc, iA); 
#else
       fdjacNeigh(n,x,fvec,fjac,vecfunc, iA); 
#endif
       /* If analytic Jacobian is available, you can 
	  replace the routine fdjac below with your own routine.*/
#ifdef MD_GLOBALNRNL
       for (i=0;i<n;i++) { /* Compute  f for the line search.*/
	 for (sum=0.0,j=0;j<n;j++)
	  sum += fjac[j][i]*fvec[j]; 
	g[i]=sum; 
      } 
      for (i=0;i<n;i++) 
	xold[i]=x[i]; /* Store x,*/ 
      fold=f; /* and f. */
#else
#ifdef MD_NEW_NR_CHECKS
      test = test_func_values(fvec, n);
#else
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvec[i]); 
#endif
      if (test < TOLFNL)
	{
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLF\n"));
	  FREERETURN;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvec[i]; /* Right-hand side for linear equations.*/
#ifdef MD_USE_LAPACK
      SolveLineq(fjac,p,n);
#else
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#endif 
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNRNL
      lnsrchNeigh(n,xold,fold,g,p,x,&f,stpmax,check,fminNeigh,iA, TOLXNL); 
      MD_DEBUG(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvec[i]) > test) 
	  test=fabs(fvec[i]); 
      if (test < TOLFNL) 
	{ 
	  *check=0; 
	  MD_DEBUG(printf("test < TOLF\n"));
	  FREERETURN
	}
      if (*check) 
	{ /* Check for gradient of f zero, i.e., spurious convergence.*/
#ifdef MD_SUPERELLIPSOID
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMINNL ? 2 : 0);
#endif
#ifdef MD_SUPERELLIPSOID
	  /* NOTA 13/04/2010: nel caso dei superellissoidi quando la distanza con il piano
	     diventa negativa è molto probabile una convergenza non spuria in x
	     quindi in tal caso ho ripristinato il check sulla convergenza spuria,
	     inoltre ho alzato le tolleranze per permettere una più facile convergenza
	     del NR. */
	  if (*check==0) /* cioè se non si tratta di convergenza 
			    spuria secondo i criteri precedenti */
	    FREERETURND
	  else
	    {
	      /* N.B. se si tratta di convergenza spuria 
		 prova a fare un normale NR incrociando le dita :-) */
	      //printf("CONVERGENZA SPURIA\n");
	      for (i=0; i < n; i++)
		{
		  x[i] = xold[i];
		  x[i] += NR_DAMP_FACT*p[i]; 
		}
	      *check = 0;
	    }

#else

	  /* se c'è anche il sospetto di un minimo locale allora fai
	   * un newton-raphson semplice */
	  for (i=0; i < n; i++)
	    {
	      x[i] = xold[i];
	      x[i] += p[i]; 
	    }
	  *check = 0;
#endif
	  MD_DEBUG(printf("*check:%d test=%f\n", *check, test));
  	  //FREERETURN 
	} 
      test=0.0; /* Test for convergence on x. */
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	} 
      if (test < TOLXNL) 
	{
	  MD_DEBUG(printf("test<TOLX test=%.15f\n", test));
	  FREERETURN;
	}
#if 1
      if (*check==2)
	{
	  MD_DEBUG(printf("spurious convergence\n"));
	  FREERETURN;
	}
#endif
#else
#ifdef MD_NEW_NR_CHECKS
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      test = test_xvalues(xold, x, n);
#else
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
#endif
      if (test < TOLXNL) 
	{ 
	  *check = 0;
	  MD_DEBUG(printf("test < TOLX\n"));
	  FREERETURN; 
	}
#endif
    } 
  MD_DEBUG10(printf("maxits!!!\n"));
  *check = 2;
  FREERETURN; 
  return;
  nrerror("MAXITS exceeded in newt"); 
}
void newt(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3])
{
  int i,j, its,ok;
  double d,f,stpmax,sum,test, fold; 
#if 0
  int *indx;
  double **fjac,*g,*p,*xold;
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvec=vector(n); 
  fvecG=vector(n);
#endif
  /*Define global variables.*/
  nn=n; 
  nn2=n-1;
  nrfuncv=vecfunc; 
  f=fminMD(x,iA,iB,shift); /*fvec is also computed by this call.*/
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
  for (i=0;i<n;i++) 
    if (fabs(fvec[i]) > test)
      test=fabs(fvec[i]); 
  if (test < 0.01*TOLF)
    {
      *check=0; 
      FREERETURN;
    }
  for (sum=0.0,i=0;i<n;i++) 
    sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  funcs2beZeroed(n,x,fvec,iA,iB,shift);
  for (its=0;its<MAXITS;its++)
    { 
      /* Start of iteration loop. */
      /* Stabilization */
       /* ============ */
       //funcs2beZeroed(n,x,fvec,iA,iB,shift);
       fdjac(n,x,fvec,fjac,vecfunc, iA, iB, shift); 
       /* If analytic Jacobian is available, you can 
	  replace the routine fdjac below with your own routine.*/
#ifdef MD_GLOBALNR
       for (i=0;i<n;i++) { /* Compute  f for the line search.*/
	 for (sum=0.0,j=0;j<n;j++)
	  sum += fjac[j][i]*fvec[j]; 
	g[i]=sum; 
      } 
      for (i=0;i<n;i++) 
	xold[i]=x[i]; /* Store x,*/ 
      fold=f; /* and f. */
#else
#ifdef MD_NEW_NR_CHECKS
      test = test_func_values(fvec, n);
#else
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvec[i]); 
#endif
      if (test < TOLF)
	{
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLF\n"));
	  FREERETURN;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvec[i]; /* Right-hand side for linear equations.*/
#ifdef MD_USE_LAPACK
      SolveLineq(fjac,p,n);
#else
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#endif 
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNR
      lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fminMD,iA,iB,shift, TOLX); 
      MD_DEBUG20(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvec[i]) > test) 
	  test=fabs(fvec[i]); 
      if (test < TOLF) 
	{ 
	  *check=0; 
	  MD_DEBUG(printf("test < TOLF\n"));
	  FREERETURN
	}
      if (*check) 
	{ 
	  /* Check for gradient of f zero, i.e., spurious convergence.*/
	  MD_DEBUG(printf("*check:%d test=%f\n", *check, test));
	  /* se c'è anche il sospetto di un minimo locale allora fai
	   * un newton-raphson semplice */
	  for (i=0; i < n; i++)
	    {
	      x[i] = xold[i];
	      x[i] += p[i]; 
	    }
	  //printf("QUIIIII <========\n");
	  *check = 0;
	  /* WARNING */
	  //continue;
	  //FREERETURN 
	} 
      test=0.0; /* Test for convergence on x. */
#if 0
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	} 
#else
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(p[i]))/FMAX(fabs(xold[i]+p[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	} 
#endif
      if (test < TOLX) 
	{
	  MD_DEBUG(printf("test<TOLX test=%.15f\n", test));
	  FREERETURN;
	}
#if 1
      if (*check==2)
	{
	  MD_DEBUG(printf("spurious convergence\n"));
	  FREERETURN;
	}
#endif
#else
#ifdef MD_NEW_NR_CHECKS
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      test = test_xvalues(xold, x, n);
#else
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
#endif
      MD_DEBUG20(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      MD_DEBUG20(printf("fvec = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLX) 
	{ 
	  *check = 0;
	  MD_DEBUG(printf("test < TOLX\n"));
	  FREERETURN; 
	}
#endif
    } 
  MD_DEBUG20(printf("maxits!!!\n"));
  *check = 2;
  FREERETURN; 
  return;
  nrerror("MAXITS exceeded in newt"); 
  
}
//#define MD_GLOBALNRD
#define MAXITS3 200
#if 0
extern double gradplane[];
extern double calc_norm(double *vec);
extern double scalProd(double *A, double *B);
#endif
extern void print_matrix(double **M, int n);

void newtDistNegNeighPlane(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int),
	  int iA)
{
  int i,its=0,ok;
  double d,stpmax,sum,test;
#ifdef MD_GLOBALNRDNL
  double fold,f, den;
  int j;
#endif
#if 0
  int * indx;
  double **fjac,*g,*p,*xold; 
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvecD=vector(n); 
#endif
  /*Define global variables.*/
  nnD=n; 
  nrfuncvDNeigh=vecfunc; 
#ifdef MD_GLOBALNRDNL
  f=fminDNeigh(x,iA); /*fvec is also computed by this call.*/
#else
  if (OprogStatus.dist5NL)
    funcs2beZeroedDistNegNeighPlane5(n,x,fvecD,iA);
  else
    funcs2beZeroedDistNegNeighPlane(n,x,fvecD,iA);
#endif
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
  for (i=0;i<n;i++) 
    if (fabs(fvecD[i]) > test)
      test=fabs(fvecD[i]); 
  //printf("INI fvec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", fvecD[0], fvecD[1], fvecD[2], fvecD[3], fvecD[4], fvecD[5], fvecD[6], fvecD[7]);
  if (test < 0.01*TOLFDNL)
    {
      *check=0; 
      FREERETURND;
    }
  for (sum=0.0,i=0;i<n;i++) 
    sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  //printf("INIZIO newt\n");
  for (its=0;its<MAXITS3;its++)
    { /* Start of iteration loop. */
       /* ============ */
      //if (x[0] > 1E4)
#if 0
       	printf("A its=%d check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",its, *check, test, x[0], x[1], x[2], x[3],x[4]);
	printf("LOOP fvec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", fvecD[0], fvecD[1], fvecD[2], fvecD[3], fvecD[4], fvecD[5], fvecD[6], fvecD[7]);
#endif
     if (OprogStatus.dist5NL)
	fdjacDistNegNeighPlane5(n,x,fvecD,fjac,vecfunc, iA);
      else
	fdjacDistNegNeighPlane(n,x,fvecD,fjac,vecfunc, iA);
      	/* If analytic Jacobian is available, you can 
	   replace the routine fdjac below with your own routine.*/
#ifdef MD_GLOBALNRDNL
       for (i=0;i<n;i++) { /* Compute  f for the line search.*/
	 for (sum=0.0,j=0;j<n;j++)
	  sum += fjac[j][i]*fvecD[j]; 
	g[i]=sum; 
      } 
      for (i=0;i<n;i++) 
	xold[i]=x[i]; /* Store x,*/ 
       //printf("BOH PRIMA lnsrchNeigh: xold=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", x[0], x[1], x[2],
	//      	 x[3], x[4], x[5], x[6], x[7]);
      fold=f; /* and f. */
#else
#ifdef MD_NEW_NR_CHECKS
      test = test_func_values(fvecD, n);
#else
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
#endif
      if (test < TOLFDNL)
	{
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLF\n"));
	  FREERETURND;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvecD[i]; /* Right-hand side for linear equations.*/
#if 0
      printf("BOH PRIMA p=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2],
	     p[3], p[4], p[5], p[6], p[7]);
      printf("xold=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", x[0], x[1], x[2],
	     x[3], x[4], x[5], x[6], x[7]);
#endif
#ifdef MD_USE_LAPACK
      SolveLineq(fjac,p,n);
#else
#if 1
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#else
      gaussj(fjac,n,p);
#endif
#endif 
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNRDNL
#if 0
      printf("PRIMA fvec.fvec=%.15G\n", scalProd(fvecD, fvecD));
      printf("prima lnsrchNeigh: xold=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", x[0], x[1], x[2],
	      	 x[3], x[4], x[5], x[6], x[7]);
#endif 
      lnsrchNeigh(n,xold,fold,g,p,x,&f,stpmax,check,fminDNeigh,iA, TOLXD); 
#if 0
      //printf("BOH (CHI CALCOLA P?) p=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2],
	//     p[3], p[4], p[5], p[6], p[7]);
      printf("its=%d check=%d xnew=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", its, *check, x[0], x[1], x[2],
            	 x[3], x[4], x[5], x[6], x[7]);
      printf("its=%d fvec: (%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G)\n",
	     its, fvecD[0], fvecD[1], fvecD[2], fvecD[3], 
	     fvecD[4], fvecD[5], fvecD[6], fvecD[7]);

#endif     
      MD_DEBUG(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));
#if 0
      printf("its=%d check=%d test = %.15f x = (%.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G )\n",its, *check, test, x[0], x[1], x[2], x[3],x[4], x[5], x[6], x[7]);
      printf("fvec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", fvecD[0], fvecD[1], fvecD[2], fvecD[3], fvecD[4],
	      fvecD[5], fvecD[6], fvecD[7]);
      //printf("gradplane=%.15G %.15G %.15G\n", gradplane[0], gradplane[1], gradplane[2]);      
#endif
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvecD[i]) > test) 
	  test=fabs(fvecD[i]); 
#if 0
      if (*check==2)
	{
	  double xpA[3], nf, fx[3], fxp[3];
	  printf("BAHBAH\n");
	  printf("g=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]);
	  printf("p=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
	  printf("fvec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", fvecD[0], fvecD[1], fvecD[2], fvecD[3], fvecD[4],
	      fvecD[5], fvecD[6], fvecD[7]);
	  //lab2body(iA, &(x[0]), xpA, rA, RtA);
	  //calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
	  //body2lab_fx(iA, fxp, fx, RtA);
	  //nf=calc_norm(fx);
          //printf("gradplane=%.15G %.15G %.15G fx=%.15G %.15G %.15G\n", gradplane[0], gradplane[1], gradplane[2],
	//	 fx[0]/nf, fx[1]/nf, fx[2]/nf);

	  printf("DOPO fvec.fvec=%.15G\n", scalProd(fvecD, fvecD));
	}
#endif
      if (test < TOLFDNL) 
	{ 
	  *check=0; 
	  MD_DEBUG(printf("test < TOLF\n"));
#if 0
	  printf("its=%d TOLFDNL fvec: (%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G)\n",
			  its, fvecD[0], fvecD[1], fvecD[2], fvecD[3], 
			  fvecD[4], fvecD[5], fvecD[6], fvecD[7]);
#endif
	  FREERETURND
	}
      if (*check==1) 
	{ /* Check for gradient of f zero, i.e., spurious convergence.*/
#ifdef MD_SUPERELLIPSOID
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    }
#if 0
	  //printf("qui test=%.15G g=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n",
	//	 test, g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]);
	  printf("its=%d qui test=%.15G x=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n",
		 its, test, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
	  printf("its=%d qui test=%.15G xold=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n",
		 its, test, xold[0], xold[1], xold[2], xold[3], xold[4], xold[5], xold[6], xold[7]);
#endif

  	  *check=(test < TOLMINDNL ? 2 : 0);
	  //printf("check=%d\n", *check);
	  //printf("MAAAH test=%.15G\n", test);
#endif
	  /* se c'è anche il sospetto di un minimo locale allora fai
	   * un newton-raphson semplice */
#ifdef MD_SUPERELLIPSOID
	  /* NOTA 13/04/2010: nel caso dei superellissoidi quando la distanza con il piano
	     diventa negativa è molto probabile una convergenza non spuria in x
	     quindi in tal caso ho ripristinato il check sulla convergenza spuria,
	     inoltre ho alzato le tolleranze per permettere una più facile convergenza
	     del NR. */
	  if (*check==0) /* cioè se non si tratta di convergenza 
			    spuria secondo i criteri precedenti */
	    {
	      //printf("HIPP its=%d qui test=%.15G fvec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n",
		// its, test, fvecD[0], fvecD[1], fvecD[2], fvecD[3], fvecD[4], fvecD[5], fvecD[6], fvecD[7]);


	      FREERETURND
	    }
	  else
	    {
	      /* N.B. se si tratta di convergenza spuria 
		 prova a fare un normale NR però con passo "damped" 
                 (incrociando le dita) :-) */
	      //printf("CONVERGENZA SPURIA\n");
	      for (i=0; i < n; i++)
		{
		  x[i] = xold[i];
		  x[i] += NR_DAMP_FACT*p[i]; 
		}
	      *check = 0;
	    }
#else
     	  for (i=0; i < n; i++)
	    {
	      x[i] = xold[i];
	      x[i] += p[i]; 
	    }
	  *check = 0;
 	  MD_DEBUG(printf("*check:%d test=%f\n", *check, test));
  	  //FREERETURND 
#endif
	} 
      test=0.0; /* Test for convergence on x. */
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	} 
      if (test < TOLXDNL) 
	{
	  //printf("test====== %.15G check=%d\n", test, *check);
	  MD_DEBUG(printf("test<TOLXD test=%.15f\n", test));
	  MD_DEBUG(printf("fvec: (%f,%f,%f,%f,%f,%f,%f,%f)\n",
			  fvecD[0], fvecD[1], fvecD[2], fvecD[2], fvecD[3], 
			  fvecD[4], fvecD[5], fvecD[6], fvecD[7]));
#if 0
	  printf("its=%d TOLXDNL fvec: (%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G,%.15G)\n",
			  its, fvecD[0], fvecD[1], fvecD[2], fvecD[3], 
			  fvecD[4], fvecD[5], fvecD[6], fvecD[7]);
#endif
	  *check=0;
	  FREERETURND;
	}
#if 1
      if (*check==2)
	{
	  MD_DEBUG(printf("spurious convergence\n"));

	  FREERETURND;
	}
#endif
#else
#ifdef MD_NEW_NR_CHECKS
#if 0
      if (iA==121)
	printf("its=%d p=%f %f %f %f %f %f %f %f\n", its, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
#endif
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      
      test = test_xvalues(xold, x, n);
#else
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
#endif
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLXDNL) 
	{ 
	  //printf("SSSSDDD test=%.15G\n", test);
	  *check = 0;
	  MD_DEBUG(printf("test < TOLX\n"));
	  FREERETURND; 
	}
#endif
    } 
  MD_DEBUG18(printf("maxits!!!\n"));
  *check = 2;
  FREERETURND;
  return;
  nrerror("MAXITS exceeded in newt"); 
  
}

void newtDistNegNeigh(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int),
	  int iA)
{
  int i,its=0,ok;
  double d,stpmax,sum,test; 
#ifdef MD_GLOBALNRDNL
  double f, fold;
  int j;
#endif
#if 0
  int *indx;
  double **fjac,*g,*p,*xold;
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvecD=vector(n); 
#endif
  /*Define global variables.*/
  nnD=n; 
  nrfuncvDNeigh=vecfunc; 
#ifdef MD_GLOBALNRDNL
  f=fminDNeigh(x,iA); /*fvec is also computed by this call.*/
#else
  funcs2beZeroedDistNegNeigh(n,x,fvecD,iA);
#endif
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
  for (i=0;i<n;i++) 
    if (fabs(fvecD[i]) > test)
      test=fabs(fvecD[i]); 
  if (test < 0.01*TOLFD)
    {
      *check=0; 
      FREERETURND;
    }
  for (sum=0.0,i=0;i<n;i++) 
    sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=0;its<MAXITS3;its++)
    { /* Start of iteration loop. */
       /* ============ */
      fdjacDistNegNeigh(n,x,fvecD,fjac,vecfunc, iA);
      /* If analytic Jacobian is available, you can 
	  replace the routine fdjac below with your own routine.*/
#ifdef MD_GLOBALNRDNL
       for (i=0;i<n;i++) { /* Compute  f for the line search.*/
	 for (sum=0.0,j=0;j<n;j++)
	  sum += fjac[j][i]*fvecD[j]; 
	g[i]=sum; 
      } 
      for (i=0;i<n;i++) 
	xold[i]=x[i]; /* Store x,*/ 
      fold=f; /* and f. */
#else
#ifdef MD_NEW_NR_CHECKS
      test = test_func_values(fvecD, n);
#else
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
#endif
      if (test < TOLFD)
	{
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLF\n"));
	  FREERETURND;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvecD[i]; /* Right-hand side for linear equations.*/
#ifdef MD_USE_LAPACK
      SolveLineq(fjac,p,n);
#else
#if 1
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#else
      gaussj(fjac,n,p);
#endif
#endif 
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNRDNL
      lnsrchNeigh(n,xold,fold,g,p,x,&f,stpmax,check,fminDNeigh,iA, TOLXD); 
      MD_DEBUG(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvecD[i]) > test) 
	  test=fabs(fvecD[i]); 
      if (test < TOLFD) 
	{ 
	  *check=0; 
	  MD_DEBUG(printf("test < TOLF\n"));
	  FREERETURND
	}
      if (*check) 
	{ /* Check for gradient of f zero, i.e., spurious convergence.*/
#if 0
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMIN ? 2 : 0);
#endif
	  /* se c'è anche il sospetto di un minimo locale allora fai
	   * un newton-raphson semplice */
	  for (i=0; i < n; i++)
	    {
	      x[i] = xold[i];
	      x[i] += p[i]; 
	    }
	  *check = 0;
 	  MD_DEBUG(printf("*check:%d test=%f\n", *check, test));
  	  //FREERETURND 
	} 
      test=0.0; /* Test for convergence on x. */
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	} 
      if (test < TOLXD) 
	{
	  MD_DEBUG(printf("test<TOLXD test=%.15f\n", test));
	  MD_DEBUG(printf("fvec: (%f,%f,%f,%f,%f,%f,%f,%f)\n",
			  fvecD[0], fvecD[1], fvecD[2], fvecD[2], fvecD[3], 
			  fvecD[4], fvecD[5], fvecD[6], fvecD[7]));
	  FREERETURND;
	}
#if 1
      if (*check==2)
	{
	  MD_DEBUG(printf("spurious convergence\n"));
	  FREERETURND;
	}
#endif
#else
#ifdef MD_NEW_NR_CHECKS
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      test = test_xvalues(xold, x, n);
#else
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
#endif
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLXD) 
	{ 
	  *check = 0;
	  MD_DEBUG(printf("test < TOLX\n"));
	  FREERETURND; 
	}
#endif
    } 
  MD_DEBUG18(printf("maxits!!!\n"));
  *check = 2;
  FREERETURND;
  return;
  nrerror("MAXITS exceeded in newt"); 
  
}
double costhrNR;
void adjust_step_dist5(double *x, double *dx, double *fx, double *gx)
{
  int k1;
  double fxfxold, gxgxold, minst, ngx, nfx, minst1, minst2;
  static int first = 1;
  static double fxold[3], gxold[3];
  nfx = calc_norm(fx);
  ngx = calc_norm(gx);
  //norm = calc_norm(dx);
  //minst = OprogStatus.toldxNR*min3(minaxA/norm,minaxAB/fabs(dx[4])/calc_norm(fx), minaxAB/fabs(x[4])/calc_norm(fx));
  //minst = OprogStatus.toldxNR*minaxAB/fabs(dx[4])/calc_norm(fx);
#if 1
  minst1 = min3(sqrt(OprogStatus.toldxNR*nfx/ngx/Sqr(dx[3])), OprogStatus.toldxNR*minaxAB/fabs(dx[4])/nfx, OprogStatus.toldxNR*nfx/ngx/2.0/fabs(dx[3]));
#else
  minst1 = min(minaxAB/fabs(dx[4])/nfx, nfx/ngx/2.0/fabs(dx[3]));
  minst1 *= OprogStatus.toldxNR;
#endif
  if (!first && OprogStatus.tolAngNR > 0.0)
    {
      fxfxold = scalProd(fx,fxold)/nfx;
      gxgxold = scalProd(gx,gxold)/ngx;
      minst2 = min(fabs((costhrNR-1.0)/(fxfxold-1.0)),
		   fabs((costhrNR-1.0)/(gxgxold-1.0)));
      minst = min(minst1,minst2);
    }
  else
    minst = minst1;
  for (k1 = 0; k1 < 5; k1++)
    dx[k1] *= min(1.0, minst);
  if (OprogStatus.tolAngNR > 0.0)
    {
      for (k1 = 0; k1 < 3; k1++)
	{
	  fxold[k1] = fx[k1]/nfx;
	  gxold[k1] = gx[k1]/ngx;
	}
    }
  if (first)
    first = 0;
}
void adjust_step_dist8(double *x, double *dx, double *fx, double *gx)
{
  int k1;
  double ngx, minst1, minst2, fxfxold, gxgxold, minst, nfx;
  static int first = 0;
  static double fxold[3], gxold[3];
  //normA = calc_norm(dx);
  //normB = calc_norm(&dx[3]);
  //minst = OprogStatus.toldxNR*min3(minaxA/normA, minaxA/normB, minaxAB/fabs(dx[7])/calc_norm(fx));
  nfx = calc_norm(fx);
  ngx = calc_norm(gx);
#if 1
  minst1 = min3(OprogStatus.toldxNR*minaxAB/fabs(dx[7])/nfx,
	       OprogStatus.toldxNR*nfx/ngx/2.0/fabs(dx[6]),
	       sqrt(OprogStatus.toldxNR*nfx/ngx/Sqr(dx[6])));
#else
  minst1 = min(minaxAB/fabs(dx[7])/nfx,nfx/ngx/2.0/fabs(dx[6]));
  minst1 *= OprogStatus.toldxNR;
#endif
  if (!first && OprogStatus.tolAngNR > 0.0)
    {
      fxfxold = scalProd(fx,fxold)/nfx;
      gxgxold = scalProd(gx,gxold)/ngx;
      minst2 = min(fabs((costhrNR-1.0)/(fxfxold-1.0)),
		   fabs((costhrNR-1.0)/(gxgxold-1.0)));
      minst = min(minst1,minst2);
    }
  else
    minst = minst1;
  for (k1 = 0; k1 < 8; k1++)
    dx[k1]*= min(minst,1.0);
  if (OprogStatus.tolAngNR > 0.0)
    {
      for (k1 = 0; k1 < 3; k1++)
	{
	  fxold[k1] = fx[k1]/nfx;
	  gxold[k1] = gx[k1]/ngx;
	}
    }
  if (first)
    first = 0;
}
extern long long itsNRdist, callsdistNR;
extern int fdjac_disterr;
#if 0
extern int check_overlp_in_calcdist(double *x, double *fx, double *gx, int iA, int iB);
#endif
void newtDistNeg(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3], int tryagain)
{
  int i,its=0,ok;
  double fx[3], gx[3];
  double d,stpmax,sum,test; 
#ifdef MD_GLOBALNRD
  double f, fold;
  int j;
#endif
#if 0
  int *indx;
  double **fjac,*g,*p,*xold;
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvecD=vector(n);
#endif
  /*Define global variables.*/
  nnD=n; 
  nrfuncvD=vecfunc; 
#ifdef MD_GLOBALNRD
  f=fminD(x,iA,iB,shift); /*fvec is also computed by this call.*/
#else
  if (OprogStatus.dist5)
    funcs2beZeroedDistNeg5(n,x,fvecD,iA,iB,shift);
  else
    funcs2beZeroedDistNeg(n,x,fvecD,iA,iB,shift);
#endif
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
  for (i=0;i<n;i++) 
    if (fabs(fvecD[i]) > test)
      test=fabs(fvecD[i]); 
  //printf("newtDistNeg BEGIN: fvec= %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", fvecD[0], fvecD[1], fvecD[2], fvecD[3], fvecD[4], fvecD[5], fvecD[6], fvecD[7]);
  if (test < 0.01*TOLFD)
    {
      *check=0; 
      FREERETURND;
    }
  for (sum=0.0,i=0;i<n;i++) 
    sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  callsdistNR++;
  
  for (its=0;its<MAXITS3;its++)
    { /* Start of iteration loop. */
       /* ============ */
      //fdjacFD(n,x,fvecD,fjac,vecfunc, iA, iB, shift); 
      if (n==8 && OprogStatus.dist8stps != 0 && its >= abs(OprogStatus.dist8stps))
	{
	  *check = 0;
	  FREERETURND;
	}
      fdjac_disterr = 0;
      if (OprogStatus.dist5)
	fdjacDistNeg5(n,x,fvecD,fjac,vecfunc, iA, iB, shift, fx, gx);
      else
	fdjacDistNeg(n,x,fvecD,fjac,vecfunc, iA, iB, shift, fx, gx);
      if (fdjac_disterr && !tryagain)
	{
	  *check = 2;
	  FREERETURND;
	}
#if 0
      /* se 1 vuol dire che almeno un punto di una SQ è dentro l'altra */ 
      if (check_overlp_in_calcdist(x, fx, gx, iA, iB))
	{
	  //printf("CHECK OVER its=%d\n", its);
	  *check=0;
	  FREERETURND;
	}
#endif
      /* If analytic Jacobian is available, you can 
	 replace the routine fdjac below with your own routine.*/
#ifdef MD_GLOBALNRD
       for (i=0;i<n;i++) { /* Compute  f for the line search.*/
	 for (sum=0.0,j=0;j<n;j++)
	  sum += fjac[j][i]*fvecD[j]; 
	g[i]=sum; 
      } 
      for (i=0;i<n;i++) 
	xold[i]=x[i]; /* Store x,*/ 
      fold=f; /* and f. */
#else
#ifdef MD_NEW_NR_CHECKS
      test = test_func_values(fvecD, n);
#else
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
      //test /= ((double) n);
#endif
      if (test < TOLFD)
	{
	  //printf("[newtDistNeg] test < TOLFD\n");
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLFD\n"));
	  FREERETURND;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvecD[i]; /* Right-hand side for linear equations.*/
#ifdef MD_USE_LAPACK
      SolveLineq(fjac,p,n);
#else
#if 1
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#ifdef MC_SIMUL
      if (ok)
	{
	  *check=2;
	  FREERETURN;
	}
#endif

#else
      gaussj(fjac,n,p);
#endif
#endif 
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNRD
      if (OprogStatus.toldxNR > 0.0)
	{
	  if (OprogStatus.dist5)
	    adjust_step_dist5(x, p, fx, gx);
	  else
	    adjust_step_dist8(x, p, fx, gx);
	} 
      lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fminD,iA,iB,shift, TOLXD); 
      MD_DEBUG(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));

      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvecD[i]) > test) 
	  test=fabs(fvecD[i]); 
      if (test < TOLFD) 
	{ 
	  *check=0; 
	  MD_DEBUG(printf("test < TOLF\n"));
	  FREERETURND
	}
      if (*check) 
	{ /* Check for gradient of f zero, i.e., spurious convergence.*/
#if 0
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMIN ? 2 : 0);
#endif
	  /* se c'è anche il sospetto di un minimo locale allora fai
	   * un newton-raphson semplice */
	  for (i=0; i < n; i++)
	    {
	      x[i] = xold[i];
	      x[i] += p[i]; 
	    }
	  //printf("=========> QUI *check=%d\n", *check);
	  *check = 0;
 	  MD_DEBUG(printf("*check:%d test=%f\n", *check, test));
	  /* WARNING */
	  //continue;
  	  //FREERETURND 
	} 
#if 0
      test=0.0; /* Test for convergence on x. */
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	}
#else
      test=0.0; /* Test for convergence on x using standard NR step and not
		   those obtained by lnsrch that can be a local minimum. */
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(p[i]))/FMAX(fabs(xold[i]+p[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	}

#endif 
#if 1
      if (test < TOLXD) 
	{
#if 0
    	  test=0.0; /* Test for convergence on function values.*/
	  for (i=0;i<n;i++) 
	    if (fabs(fvecD[i]) > test) 
	      test=fabs(fvecD[i]); 
 
	  printf("QUI testTOLF=%.15G?!?\n", test);
#endif
	  MD_DEBUG(printf("test<TOLXD test=%.15f\n", test));
	  MD_DEBUG(printf("fvec: (%f,%f,%f,%f,%f,%f,%f,%f)\n",
			  fvecD[0], fvecD[1], fvecD[2], fvecD[2], fvecD[3], 
			  fvecD[4], fvecD[5], fvecD[6], fvecD[7]));
#ifdef MD_SUPERELLIPSOID
	  *check = 0;
#endif
	  FREERETURND;
	}
#endif
#if 1
      if (*check==2)
	{
	  MD_DEBUG(printf("spurious convergence\n"));
	  FREERETURND;
	}
#endif
#else
      if (OprogStatus.toldxNR > 0.0)
	{
	  MD_DEBUG(printf("qui?!?\n"));
	  if (OprogStatus.dist5)
	    adjust_step_dist5(x, p, fx, gx);
	  else
	    adjust_step_dist8(x, p, fx, gx);
	}
#ifdef MD_NEW_NR_CHECKS
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      test = test_xvalues(xold, x, n);
#else
      test = 0;
      for (i=0;i<n;i++) 
	{ 
	  test += fabs(p[i]);
	  x[i] += p[i];
	}
      //test /= ((double) n);
#endif
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLXD) 
	{ 
	  //printf("test<TOLXD\n");
	  *check = 0;
	  MD_DEBUG(printf("test < TOLXD\n"));
	  FREERETURND; 
	}
#endif
      itsNRdist++;
    } 
  MD_DEBUG18(printf("maxits!!!\n"));
  *check = 2;
  FREERETURND;
  return;
  nrerror("MAXITS exceeded in newt"); 
  
}

void newtDist(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3])
{
  int i, its,ok;
  double d,stpmax,sum,test; 
#ifdef MD_GLOBALNRD
  double f, fold;
  int j;
#endif
#if 0
  int *indx;
  double **fjac,*g,*p,*xold;
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvecD=vector(n); 
#endif
  /*Define global variables.*/
  nnD=n; 
  nrfuncvD=vecfunc; 
#ifdef MD_GLOBALNRD
  f=fminD(x,iA,iB,shift); /*fvec is also computed by this call.*/
#else
  funcs2beZeroedDist(n,x,fvecD,iA,iB,shift);
#endif
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
  for (i=0;i<n;i++) 
    if (fabs(fvecD[i]) > test)
      test=fabs(fvecD[i]); 
  if (test < 0.01*TOLFD)
    {
      *check=0; 
      FREERETURND;
    }
  for (sum=0.0,i=0;i<n;i++) 
    sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=0;its<MAXITS3;its++)
    { /* Start of iteration loop. */
       /* ============ */
      //fdjacFD(n,x,fvecD,fjac,vecfunc, iA, iB, shift); 
      fdjacDist(n,x,fvecD,fjac,vecfunc, iA, iB, shift);
       /* If analytic Jacobian is available, you can 
	  replace the routine fdjac below with your own routine.*/
#ifdef MD_GLOBALNRD
       for (i=0;i<n;i++) { /* Compute  f for the line search.*/
	 for (sum=0.0,j=0;j<n;j++)
	  sum += fjac[j][i]*fvecD[j]; 
	g[i]=sum; 
      } 
      for (i=0;i<n;i++) 
	xold[i]=x[i]; /* Store x,*/ 
      fold=f; /* and f. */
#else
#ifdef MD_NEW_NR_CHECKS
      test = test_func_values(fvecD, n);
#else
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
#endif
      if (test < TOLFD)
	{
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLF\n"));
	  FREERETURND;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvecD[i]; /* Right-hand side for linear equations.*/
#ifdef MD_USE_LAPACK
      SolveLineq(fjac,p,n);
#else
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#endif
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNRD
      lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fminD,iA,iB,shift, TOLXD); 
      MD_DEBUG(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvecD[i]) > test) 
	  test=fabs(fvecD[i]); 
      if (test < TOLFD) 
	{ 
	  *check=0; 
	  MD_DEBUG(printf("test < TOLF\n"));
	  FREERETURND
	}
      if (*check) 
	{ /* Check for gradient of f zero, i.e., spurious convergence.*/
#if 0
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMIN ? 2 : 0);
#endif
	  /* se c'è anche il sospetto di un minimo locale allora fai
	   * un newton-raphson semplice */
	  for (i=0; i < n; i++)
	    {
	      x[i] = xold[i];
	      x[i] += p[i]; 
	    }
	  *check = 0;
	  MD_DEBUG(printf("*check:%d test=%f\n", *check, test));
	  /* WARNING */
	  //continue;
  	  //FREERETURND 
	} 
      test=0.0; /* Test for convergence on x. */
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0); 
	  if (temp > test) 
	    test=temp; 
	} 
      if (test < TOLXD) 
	{
	  MD_DEBUG(printf("test<TOLXD test=%.15f\n", test));
	  MD_DEBUG(printf("fvec: (%f,%f,%f,%f,%f,%f,%f,%f)\n",
			  fvecD[0], fvecD[1], fvecD[2], fvecD[2], fvecD[3], 
			  fvecD[4], fvecD[5], fvecD[6], fvecD[7]));
	  FREERETURND;
	}
#if 1
      if (*check==2)
	{
	  MD_DEBUG(printf("spurious convergence\n"));
	  FREERETURND;
	}
#endif
#else
#ifdef MD_NEW_NR_CHECKS
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      test = test_xvalues(xold, x, n);
#else
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
#endif
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLXD) 
	{ 
	  *check = 0;
	  MD_DEBUG(printf("test < TOLX\n"));
	  FREERETURND; 
	}
#endif
    } 
  MD_DEBUG(printf("maxits!!!\n"));
  *check = 2;
  FREERETURND; 
  return;
  nrerror("MAXITS exceeded in newt"); 
  
}


#define EPS 3.0E-8 /* Approximate square root of the machine precision.*/
void fdjacFD(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int, int, double []), int iA, int iB, double shift[3])
{ int i,j; 
  double h,temp,*f; 
  f=vector(n); 
  for (j=0;j<n;j++) 
    {
      temp=x[j]; 
      h=EPS*fabs(temp);
      if (h == 0.0)
	h=EPS; 
      x[j]=temp+h; 
      /* Trick to reduce  nite precision error.*/
      h=x[j]-temp; 
      (*vecfunc)(n,x,f, iA, iB, shift); 
      x[j]=temp;
      for (i=0;i<n;i++)
	df[i][j]=(f[i]-fvec[i])/h; /* Forward difference*/
    }
  free_vector(f); 
}
extern int nn, nn2; 
extern void (*nrfuncv)(int n, double v[], double f[], int iA, int iB, double shift[3]); 
extern void (*nrfuncv2)(int n, double v[], double f[], int iA, int iB, double shift[3]); 
//extern void (*nrfuncvNeigh)(int n, double v[], double f[], int iA); 
//extern void (*nrfuncvDNeigh)(int n, double v[], double f[], int iA);
double fminNeigh(double x[], int iA)
{
  int i;
  double sum;
  (*nrfuncvNeigh)(nn,x,fvec,iA);
  for (sum=0.0,i=0;i<nn;i++)
    sum += Sqr(fvec[i]); 
    return 0.5*sum; 
}
double fminDNeigh(double x[], int iA) 
/* Returns f = 1 2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the 
calling program. Global variables also communicate the function values back to 
the calling program.*/
{
  int i;
  double sum;
  (*nrfuncvDNeigh)(nnD,x,fvecD,iA);
  for (sum=0.0,i=0;i<nnD;i++)
    sum += Sqr(fvecD[i]); 
  return 0.5*sum; 
}
double fminMD(double x[], int iA, int iB, double shift[3]) 
/* Returns f = 1/2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the 
calling program. Global variables also communicate the function values back to 
the calling program.*/
{
  int i;
  double sum;
  (*nrfuncv)(nn,x,fvec,iA,iB,shift);
  for (sum=0.0,i=0;i<nn;i++)
    sum += Sqr(fvec[i]); 
    return 0.5*sum; 
}
double fminD(double x[], int iA, int iB, double shift[]) 
/* Returns f = 1 2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the 
calling program. Global variables also communicate the function values back to 
the calling program.*/
{
  int i;
  double sum;
  (*nrfuncvD)(nnD,x,fvecD,iA,iB,shift);
  for (sum=0.0,i=0;i<nnD;i++)
    sum += Sqr(fvecD[i]); 
  return 0.5*sum; 
}
