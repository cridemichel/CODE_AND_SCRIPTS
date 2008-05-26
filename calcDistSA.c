#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define Sqr(x) ((x)*(x))
#define MD_ASYM_ITENS
#define TINY 1E-20
#define MD_NBMAX 8 
#define FREERETURND return;
#define FREERETURN  return;
#define ALF 1.0e-4 /* Ensures sufficient decrease in function value.*/
#define TOLX 1.0E-14//1.0e-7 /* Convergence criterion on  x.*/ 
#define TOLXD 1.0E-14
#define MAXITS 100 // se le particelle non si urtano il newton-raphson farà MAXITS iterazioni
#define MAXITS2 100
#define TOLF 1.0E-10// 1.0e-4
#define TOLFD 1.0E-10
#define TOLMIN 1.0E-12//1.0e-6 
#define STPMX 100.0
#define MAXITS3 200
double maxarg1,maxarg2;
double costhrNR;
long long itsNRdist=0, callsdistNR=0;
double minaxA, minaxB, minaxAB;
void calc_intersec(double *rB, double *rA, double **Xa, double* rI);

#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

long long int itsF=0, timesF=0, itsS=0, timesS=0, numcoll=0, itsFNL=0, timesFNL=0, 
     timesSNL=0, itsSNL=0, numcalldist=0, numdisttryagain=0;
const double springkSD=1.0, stepSD=0.05, tolSDgrad=0.2, tolAngNR=0.0, epsd=1E-5;
const int SDmethod  = 2, forceguess=1, guessDistOpt=1, dist8stps=0, maxitsSD=1000;
const double epsdGDO=0.0, toldxNR=0.0, tolSD=0.01, tolAngSD=0.0, tolSDconstr=0.002, tolSDlong=0.01;
const double epsdSD = 0.01;
int fdjac_disterr;
long long int itsfrprmn=0, callsfrprmn=0,callsok=0, callsprojonto=0, itsprojonto=0;
long long accngA=0, accngB=0;
double costolSDgrad, saxA[3]={2,1,1}, saxB[3]={2,1,1}, posA[3]={0.0,0.0,-3.0}, posB[3]={0.0,0.0,+3.0};
double RinpA[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double RinpB[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double gradfG[3], gradgG[3], dxG[6];
double sfA, sfB;
double invaSq, invbSq, invcSq, costolAngSD;
double **XbXa, **Xa, **Xb, **RA, **RB, **Rt, **RtA, **RtB, r1[3], r2[3];
double rA[3], rB[3];
void (*nrfuncv)(int n, double v[], double fvec[], double shift[3]);
void (*nrfuncvD)(int n, double v[], double fvec[], double shift[3]);
int nn, nn2, nnD; /* Global variables to communicate with fmin.*/
double *fvecD;
double **fjac,*g,*p,*xold;
int *indx;
double shiftcg[3], lambdacg, minaxicg, minaxjcg;

int cghalfspring, icg, jcg, doneryck;
double *vector(int n)
{
  return calloc(n, sizeof(double)); 
}
double **matrix(int n, int m)
{
  double **M;
  int i;
  M = malloc(sizeof(double*)*n);
  for (i=0; i < n; i++)
    M[i] = calloc(m, sizeof(double));
  return M;
}

/* Returns f = 1 2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the 
calling program. Global variables also communicate the function values back to 
the calling program.*/
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
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
double fminD(double x[], double shift[]) 
{
  int i;
  double sum;
  (*nrfuncvD)(nnD,x,fvecD,shift);
  for (sum=0.0,i=0;i<nnD;i++)
    sum += Sqr(fvecD[i]); 
  return 0.5*sum; 
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
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
  minst1 = min3(toldxNR*minaxAB/fabs(dx[7])/nfx,
	       toldxNR*nfx/ngx/2.0/fabs(dx[6]),
	       sqrt(toldxNR*nfx/ngx/Sqr(dx[6])));
#else
  minst1 = min(minaxAB/fabs(dx[7])/nfx,nfx/ngx/2.0/fabs(dx[6]));
  minst1 *= toldxNR;
#endif
  if (!first && tolAngNR > 0.0)
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
  if (tolAngNR > 0.0)
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

void calc_grad(double *rC, double *rA, double **Xa, double *grad)
{
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      grad[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	grad[k1] += 2.0*Xa[k1][k2]*(rC[k2]-rA[k2]); 
    }
}
void projonto(double* ri, double *dr, double* rA, double **Xa, double *gradf, double *sfA, double dist)
{
  int kk, its, done=0, k1, k2, LMAXITS=50;
  const double GOLD=1.618034;
  double Xag[3], r1AXa[3], r1[3], r1A[3], sf, sqrtDelta, A2;
  double curv2, A, B, C, Delta, sol=0.0, ng, curv, lambda;
  sf = *sfA;
  its = 0;
 
  while (!done && its <= LMAXITS)
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
      if ((SDmethod == 2 || SDmethod == 4) && tolAngSD > 0.0)
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
      if ((SDmethod == 2 || SDmethod == 4) && tolSDconstr > 0.0)
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
	  lambda = min(sf,sqrt(tolSDconstr/fabs(curv2)));
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
      if (SDmethod==1 || SDmethod==3)
	{
	  if (dist > epsd && fabs(sol)*ng > tolSDconstr*dist/2.0)
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
  if (SDmethod == 1 || SDmethod==3)
    *sfA = sf;
}

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

  if (forceguess)
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
  if (SDmethod==1 || SDmethod==3)
    {
      K1= stepSD;
      K2= stepSD;
    }
  else
    {
      K1 = springkSD;
      K2 = springkSD;
    }
  for (kk=0; kk < 3; kk++)
    {
      if (SDmethod == 1 || SDmethod==3)
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
  if (tolSDgrad > 0.0)
    {
      if (SDmethod==1 || SDmethod==3)
	{
	  ngA = ngB = 0;  
	  for (kk=0; kk < 3; kk++)
	    {
	      ngA += Sqr(grad[kk]);
	      ngB += Sqr(grad[kk+3]); 
	      
	    }
	  if (sqrt(ngA) < tolSDgrad*stepSD
	      && sqrt(ngB) < tolSDgrad*stepSD)
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
  S *= springkSD;
  F = S*Sqr(normdd);
  return F; 
}
/* =========================== >>> forces <<< ======================= */
double  cgfuncRyck(double *vec)
{
  int kk, k1, k2;
  double fx2[3], gx2[3];
  double A, B, F;
  if (forceguess)
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
	A = - springkSD;
      else
	A = springkSD;
    }
  else
    A = springkSD;
 
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

int check_done(double fp, double fpold, double minax)
{
  const double EPSFR=1E-10;
  double dold=0.0, d=0.0;
  if (SDmethod != 2 && SDmethod != 4)
    {
      dold = sqrt(fpold/springkSD);
      d = sqrt(fp/springkSD);
    }
  if (tolSDgrad > 0)
    {
      if (SDmethod == 2 || SDmethod==4)
	{
	  if (epsdSD > 0.0 && fp < Sqr(epsdSD))
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
	  if (fp > Sqr(epsdSD))//Sqr(OprogStatus.epsd)) 
	    {
	      if (doneryck == 1 || 
		  2.0*fabs(dold-d) < tolSDlong*(fabs(dold)+fabs(d)+EPSFR))
		{
		  accngA++;
		  accngB++;
		  return 1;
		}
	    }
	  else 
	    {
	      if (2.0*fabs(dold-d) < tolSD*(fabs(dold)+fabs(d)+EPSFR))
		return 1;
	    }
	}
    }
  else
    {
      if (fp < Sqr(epsdSD))//Sqr(OprogStatus.epsd))
	{
	  if (2.0*fabs(dold-d) < tolSD*(fabs(dold)+fabs(d)+EPSFR))
	    return 1;
	}
      else
	{
	  if (2.0*fabs(dold-d) < tolSDlong*(fabs(dold)+fabs(d)+EPSFR))
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
  const int ITMAXFR = maxitsSD;
  const double GOLD=1.618034;
  double fp, fpold=0.0, signA, signB;
  double minax, xi[6], xiold[6];
  double signAold, signBold, pold[6];
  //printf("primaprima p= %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2], p[3], p[4], p[5]);
 
  minax = min(minaxicg,minaxjcg);
  sfA = stepSD;
  sfB = stepSD;
  /*Initializations.*/
  fp = (*dfunc)(p,xi,gradfG,gradgG, &signA, &signB); 
  
  if (doneryck==2)
    {
      callsok++;
      return;
    }
#if 1
  if ((SDmethod == 2 || SDmethod == 4) &&
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

      if ((SDmethod == 1 || SDmethod == 3) && fp > fpold)
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
  
}


void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], 
	    double *f, double stpmax, int *check, 
	    double (*func)(double [], double[]), double shift[3],
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
      printf("Roundoff problem in lnsrch."); 
      exit(-1);
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
      *f=(*func)(x,shift); 
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
void fdjacFD(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double [], int, int, double []), double shift[3]);
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

void funcs2beZeroedDistNeg(int n, double x[], double fvec[], double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3];
  /* x = (r, alpha, t) */ 
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
}
void fdjacDistNeg(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], double []), double shift[3], double *fx, double *gx)
{
  int kk;
  double axi[3], axj[3];
  double rDC[3];
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = 2.0*Sqr(x[6])*Xb[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3]-rB[k2]);
	}
    }
#if 1
  if (SDmethod == 2 || SDmethod == 3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = x[k1+3] - x[k1];
      for (kk=0; kk < 3; kk++)
	{
	  axi[kk] = saxA[kk];
	  axj[kk] = saxB[kk];
	}
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axi[0],axi[1],axi[2])+max3(axj[0],axj[1],axj[2])))
	{
	  fdjac_disterr = 1;	
	}
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  for (k1 = 0; k1 < 5; k1++)
    {
      df[3][k1+3] = 0;
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = 0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1+3] = gx[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = 2.0*x[6]*gx[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1 + 2.0*x[7]*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 2.0*x[7]*Xa[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = fx[k1];
#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
}


void newtDistNeg(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], double []),
	  double shift[3], int tryagain)
{
  int i,its=0,ok;
  double fx[3], gx[3];
  double d,stpmax,sum,test; 
  /*Define global variables.*/
  nnD=n; 
  nrfuncvD=vecfunc; 
#ifdef MD_GLOBALNRD
  f=fminD(x,shift); /*fvec is also computed by this call.*/
#else
  funcs2beZeroedDistNeg(n,x,fvecD,shift);
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
  //callsdistNR++;
  
  for (its=0;its<MAXITS3;its++)
    { /* Start of iteration loop. */
       /* ============ */
      //fdjacFD(n,x,fvecD,fjac,vecfunc, iA, iB, shift); 
      if (n==8 && dist8stps != 0 && its >= abs(dist8stps))
	{
	  *check = 0;
	  FREERETURND;
	}
      fdjac_disterr = 0;
      fdjacDistNeg(n,x,fvecD,fjac,vecfunc, shift, fx, gx);
      if (fdjac_disterr && !tryagain)
	{
	  *check = 2;
	  FREERETURND;
	}
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
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
      if (test < TOLFD)
	{
	  *check = 0;
	  FREERETURND;
	}
#endif 
      for (i=0;i<n;i++) 
	p[i] = -fvecD[i]; /* Right-hand side for linear equations.*/
#if 1
      ludcmp(fjac,n,indx,&d, &ok); /* Solve linear equations by LU decomposition.*/
      lubksb(fjac,n,indx,p);
#else
      gaussj(fjac,n,p);
#endif
      /* lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
#ifdef MD_GLOBALNRD
      if (toldxNR > 0.0)
	{
	  adjust_step_dist8(x, p, fx, gx);
	} 
      lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fminD,shift, TOLXD); 

      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	if (fabs(fvecD[i]) > test) 
	  test=fabs(fvecD[i]); 
      if (test < TOLFD) 
	{ 
	  *check=0; 
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
	  FREERETURND;
	}
#if 1
      if (*check==2)
	{
	  FREERETURND;
	}
#endif
#else
      if (toldxNR > 0.0)
	{
	  adjust_step_dist8(x, p, fx, gx);
	}
      test = 0;
      for (i=0;i<n;i++) 
	{ 
	  test += fabs(p[i]);
	  x[i] += p[i];
	}
      if (test < TOLXD) 
	{ 
	  *check = 0;
	  FREERETURND; 
	}
#endif
      itsNRdist++;
    } 
  *check = 2;
  FREERETURND;
  return;
}

void tRDiagR(double **M, double a, double b, double c, double **Ri)
{
  int na;
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = a;
  Di[1][1] = b;
  Di[2][2] = c;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
void print_matrix(double **M, int n)
{
  int k1, k2;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}


void guess_dist(double *rA, double *rB, double **Xa, double **Xb, double *rC, double *rD,
		double **RA, double **RB)
{
  double gradA[3], gradB[3], gradaxA[3], gradaxB[3], dA[3], dB[3];
  int k1, n;
  double saA[3], saB[3];
  int typei, typej;
  saA[0] = saxA[0];
  saA[1] = saxA[1];
  saA[2] = saxA[2];
  saB[0] = saxB[0];
  saB[1] = saxB[1];
  saB[2] = saxB[2];
  //printf("axes[%d]=%f,%f,%f axes[%d]=%f,%f,%f\n", i, axa[i], axb[i], axc[i], j, axa[j], axb[j], axc[j]);
  for (k1 = 0; k1 < 3; k1++)
    {
      gradA[k1] =  (rB[k1]-rA[k1]);
      gradB[k1] = -(rB[k1]-rA[k1]);
    }
  for (n = 0; n < 3; n++)
    {
      gradaxA[n] = 0;
      gradaxB[n] = 0;
      for (k1 = 0; k1 < 3; k1++) 
	{
	  gradaxA[n] += gradA[k1]*RA[n][k1];
	  gradaxB[n] += gradB[k1]*RB[n][k1];
	}
    }
  for (k1=0; k1 < 3; k1++)
    {
      dA[k1] = rA[k1];
      dB[k1] = rB[k1];
      for (n=0; n < 3;n++)
	{
	  dA[k1] += gradaxA[n]*RA[n][k1]*saA[n]/2.0; 
	  dB[k1] += gradaxB[n]*RB[n][k1]*saB[n]/2.0;
	}
    }
  calc_intersec(dA, rA, Xa, rC);
  calc_intersec(dB, rB, Xb, rD);
}

void calc_intersec(double *rB, double *rA, double **Xa, double* rI)
{
  double A, B=0.0, C=0.0, D=0.0, tt=0.0;
  double rBA[3];
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = rB[k1] - rA[k1];
  A = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
    	A += rBA[k1]*Xa[k1][k2]*rBA[k2];
      }
  if (A <= 0)
    {
      printf("[calc_intersec] Serious problem guessing distance, aborting...\n");
      printf("tt = %f D=%f A=%f B=%f C=%f\n", tt, D, A, B, C);
      printf("distance: %f\n", sqrt(Sqr(rBA[0])+Sqr(rBA[1])+Sqr(rBA[2])));
      print_matrix(Xa,3);
      print_matrix(Xb,3);
      exit(-1);
    }
  tt = sqrt(1 / A); 
  for (k1 = 0; k1 < 3; k1++)
    {
      rI[k1] = rA[k1] + tt*rBA[k1];  
    }
}
void distSD(double shift[3], double *vecg, double lambda, int halfspring)
{
  int kk;
  double Fret;
  int iter;
  double vec[8];
  int typei, typej;
  double axaiF, axbiF, axciF, axajF, axbjF, axcjF;
  axaiF = saxA[0];
  axbiF = saxA[1];
  axciF = saxA[2];
  minaxicg = min3(axaiF, axbiF, axciF);
  axajF = saxB[0];
  axbjF = saxB[1];
  axcjF = saxB[2];
  minaxjcg = min3(axajF, axbjF, axcjF);
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
  frprmnRyck(vec, 6, tolSD, &iter, &Fret, cgfuncRyck, gradcgfuncRyck);
  for (kk=0; kk < 6; kk++)
    {
      vecg[kk] = vec[kk];
    }
}


double calcDistNeg(double rrA[3], double RRA[3][3], double ssaxA[3], double rrB[3], double RRB[3][3], 
		   double ssaxB[3])
{
  /* SDmethod=1 usa la riduzione del passo nello Steepest Descent (SD) e applica lo SD sempre
     SDmethod=2 non usa la riduzione del passo e applica lo SD solo se il calcolo della distanza fallisce 
     SDmethod=3 usa la riduzione del passo e applica lo SD solo se il calcolo della distanza fallisce 
     SDmethod=4 non usa la riduzione del passo e applica lo SD sempre */
  double shift[3]={0.0,0.0,0.0}; 
  double vecgsup[8]={0,0,0,0,0,0,0,0};
  double alpha;
  int aa, bb;
  int calcguess=1;
  double vecg[8], rC[3], rD[3], rDC[3], r12[3], vecgcg[6], fx[3];
  double ti, segno, segno2;
  double g1=0.0, g2=0.0, SP, nrDC, vecnf[3], nvecnf;
  int retcheck, tryagain = 0;
  double axaiF, axbiF, axciF;
  double axajF, axbjF, axcjF;
#ifndef MD_ASYM_ITENS
  double Omega[3][3];
#endif
  double nf, ng, gradf[3], gradg[3];
  int k1, k2, na, k3;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  costolSDgrad = cos(tolSDgrad);
  costhrNR = cos(tolAngNR);
  costolAngSD = cos(tolAngSD);
  fvecD=vector(8);
  XbXa = matrix(3, 3);
  Xa = matrix(3, 3);
  Xb = matrix(3, 3);
  RA = matrix(3, 3);
  RB = matrix(3, 3);
  Rt = matrix(3, 3);
  RtA = matrix(3, 3);
  RtB = matrix(3, 3);
  for (aa=0; aa < 3; aa++)
    for (bb=0; bb < 3; bb++)
      {
	RtA[aa][bb] = RRA[aa][bb];
	RtB[aa][bb] = RRB[aa][bb];	
      }
  axaiF = ssaxA[0];
  axbiF = ssaxA[1];
  axciF = ssaxA[2];
  minaxA = min3(axaiF,axbiF,axciF);
  axajF = ssaxB[0];
  axbjF = ssaxB[1];
  axcjF = ssaxB[2];
  minaxB = min3(axajF,axbjF,axcjF);
  minaxAB = min(minaxA,minaxB);
  rA[0] = rrA[0];
  rA[1] = rrA[1];
  rA[2] = rrA[2];
  /* ...and now orientations */
  invaSq = 1/Sqr(ssaxA[0]);
  invbSq = 1/Sqr(ssaxA[1]);
  invcSq = 1/Sqr(ssaxA[2]);
  tRDiagR(Xa, invaSq, invbSq, invcSq, RtA);

  rB[0] = rrB[0];
  rB[1] = rrB[1];
  rB[2] = rrB[2];
  invaSq = 1.0/Sqr(ssaxB[0]);
  invbSq = 1.0/Sqr(ssaxB[1]);
  invcSq = 1.0/Sqr(ssaxB[2]);
  tRDiagR(Xb, invaSq, invbSq, invcSq, RtB);

retry:
  if (forceguess)
    calcguess = 1;
  if (calcguess || tryagain)
    {
      if (guessDistOpt==1)
	{
	  guess_dist(rA, rB, Xa, Xb, rC, rD, RtA, RtB);
	}
      else
	{
	  calc_intersec(rB, rA, Xa, rC);
	  calc_intersec(rA, rB, Xb, rD);
	}
#if 1
      for(k1=0; k1 < 3; k1++)
	r12[k1] = rC[k1]-rD[k1]; 
      if ((SDmethod==1 || SDmethod==4) || tryagain)
	{
	  for (k1=0; k1 < 3; k1++)
	    {
	      vecgcg[k1] = rC[k1];
	      vecgcg[k1+3] = rD[k1];
	    }
	  distSD(shift, vecgcg, springkSD, 1);
	  for (k1=0; k1 < 3; k1++)
	    {
	      rC[k1] = vecgcg[k1];
	      rD[k1] = vecgcg[k1+3];
	    }	 
#endif
	}
      calc_grad(rC, rA, Xa, gradf);
      calc_grad(rD, rB, Xb, gradg);
      nf = calc_norm(gradf);
      ng = calc_norm(gradg);
      vecg[6] = sqrt(nf/ng);
      for (k1=0; k1 < 3; k1++)
	{
	  vecg[k1] = rC[k1];
	  vecg[k1+3] = rD[k1];
	  rDC[k1] = rD[k1] - rC[k1];
	}
      if (epsdGDO > 0.0)
	{
	  g1 = calc_norm(rDC)/nf;
	  nrDC = calc_norm(rDC);
	  SP = scalProd(rDC,gradf)/nf;
	  for (k1=0; k1 < 3; k1++)
	    {
	      vecnf[k1] = rDC[k1] - SP*gradf[k1]; 
	    }
	  nvecnf = calc_norm(vecnf);
	  if ( nvecnf > 0.0)
	    g2 = epsdGDO*min3(axajF,axbjF,axcjF)/calc_norm(vecnf); 
	  else 
	    g2 = g1;
	}	  
      
#if 1
      if ((SDmethod==1 || SDmethod==4) || tryagain)
	{
	  if (scalProd(gradf, rDC) < 0.0)
	    vecg[7] = 0.0;
	  else
	    vecg[7] = calc_norm(rDC)/nf;  
	}
      else
	{
	  if (epsdGDO > 0.0)
	    {
	      if (scalProd(gradf, rDC) < 0.0)
		vecg[7] = 0.0;
	      else
		vecg[7] = min(g1,g2);
	    }
	  else
	    vecg[7] = 0.0;
	}
#endif
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  
  newtDistNeg(vecg, 8, &retcheck, funcs2beZeroedDistNeg, shift, tryagain); 
  //numcalldist++;
  if (retcheck != 0)
    {
      if (tryagain)
	{
	  printf("[ERROR] I'm sorry but I can't really calculate distance\n");
	  exit(-1);
     	}
	 
      if (!tryagain && (SDmethod == 2 || SDmethod==3))
	{
	  numdisttryagain++;
	  tryagain = 1; 
	  goto retry;
	}
      if (calcguess==0)
	{
	  calcguess=2;
	  goto retry;
	} 
      printf("[calcDistNeg] I couldn't calculate distance between HEs exiting....\n");
      exit(-1);
    }
  for (k1 = 0; k1 < 8; k1++)
    {
      vecgsup[k1] = vecg[k1]; 
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      r2[k1] = vecg[k1+3];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      r12[k1] = r1[k1] - r2[k1];
    } 
  
  alpha = vecg[3];
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (r2[k1]-rA[k1])*Xa[k1][k2]*(r2[k2]-rA[k2]); 

  if (segno*vecg[7]<0 && fabs(segno*vecg[7])>3E-8)
    {

      if (tryagain)
	{
	  printf("[ERROR] I'm sorry but I can't really calculate distance\n");
	  exit(-1);
	} 
      if (!tryagain && ( SDmethod==2 || SDmethod==3 ))
	{
	  numdisttryagain++;
	  tryagain = 1; 
	  goto retry;
	}
      printf("1) segno: %.8G vecg[7]: %.8G dist=%.15G\n", segno, vecg[7], calc_norm(r12));
      return -1.0E10;//calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
    }
  segno2 = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno2 += (r1[k1]-rB[k1])*Xb[k1][k2]*(r1[k2]-rB[k2]); 
  if (segno2*segno < 0.0 && fabs(segno*segno2) > 3E-8)
    {
      if (tryagain)
	{
	  printf("[ERROR segno*segno2] I'm sorry but I can't really calculate distance\n");
	  exit(-1);
	} 
      if (!tryagain && ( SDmethod==2 || SDmethod==3 ))
	{
	  numdisttryagain++;
	  tryagain = 1; 
	  goto retry;
	}
      printf("2) segno: %.8G segno2: %.15G dist=%.15G\n", segno, segno2, calc_norm(r12));
      return 1.0E-10;//calcDist(t, t1, i, j, shift, r1, r2, alpha, vecgsup, 1);
    }
  if (segno > 0)
    return calc_norm(r12);
  else
    return -calc_norm(r12);
}
int main(int argc, char **argv)
{
  printf("Dist=%.15G\n", calcDistNeg(posA, RinpA, saxA, posB, RinpB, saxB));
  return 0;
}
