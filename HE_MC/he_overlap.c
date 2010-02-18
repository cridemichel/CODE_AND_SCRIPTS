#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define TINY 1E-20
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define DABS fabs
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
void nrerror(char *msg)
{
  printf(msg);
  exit(-1);
}
int brentTooManyIter;
const TOL=1.0E-15, LAM0=0.0, LAM1=0.5, LAM2=1.0;
double saA[3], saB[3], RA[3][3], RB[3][3];
double COMA[3], COMB[3];
void ludcmp(double a[3][3], int n,  int* indx, double* d, int *ok)
{
  /* A[i][j] = Aij 
   * A x = b  
   * per semplicità nel seguito si assume che l'ordine della matrice è 3 */
  int i,imax=-1,j,k;
  double big,dum,sum,temp; 
  double vv[3]; /* vv stores the implicit scaling of each row.*/
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

void lubksb(double a[3][3], int n, int* indx, double *b)
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

#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f);
double dbrent(double ax, double bx, double cx, void (*f)(double,double*,double*), double tol, double *xmin) 
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
  a=(ax < cx ? ax : cx); b=(ax > cx ? ax : cx); x=w=v=bx; 
  
  (*f)(x,&fx,&dx);
  fw=fv=fx; 
  dw=dv=dx; 
  /*All our housekeeping chores are doubled by the necessity of moving derivative values around as well 
   * as function values. */
  for (iter=1;iter<=ITMAXDBR;iter++)
    { 
      //printf("x=%.15G iter=%d\n", x, iter );
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
	  u=x+d; 
	  (*f)(u,&fu,&du);
	} 
      else
	{ 
	  u=x+SIGN(tol1,d); 
	  (*f)(u, &fu, &du);
	  if (fu > fx) 
	    { /*If the minimum step in the downhill direction takes us uphill, then we are done.*/
	      *xmin=x;
	      return fx;
	    }
	}
      //du=(*df)(u); /*Now all the housekeeping, sigh.*/
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
const int USEDBRENT=1, READFROMFILE=1,USEDIRINV=1;
void MATINV(double A[3][3],double Ainv[3][3])
{
  double D, col[3];
  int i,j,INDX[3], ok;
  for (i=0; i < 3; i++)
    {
      for (j=0; j < 3; j++)
	Ainv[i][j]=0.0;
      Ainv[i][i]=1.0;
    }
  ludcmp(A,3,INDX,&D, &ok);
  for (j=0; j < 3; j++)
    {
      for (i=0; i < 3; i++)
	col[i]=0.0;
      col[j]=1.0;
      lubksb(A,3,INDX,col);
      for (i=0; i < 3; i++)
	Ainv[i][j] = col[i];
    }
}
void MATINV33(double A[3][3],double Ainv[3][3])
{
  double DET;
  int i, j;
  DET = A[0][0]*A[1][1]*A[2][2]-A[0][0]*A[1][2]*A[2][1]-A[0][1]*A[1][0]
    *A[2][2]+A[0][1]*A[1][2]*A[2][0]+A[0][2]*A[1][0]*A[2][1]-A[0][2]*A[1][1]
    *A[2][0];
  Ainv[0][0]=A[1][1]*A[2][2]-A[1][2]*A[2][1];
  Ainv[0][1]=A[0][2]*A[2][1]-A[0][1]*A[2][2];
  Ainv[0][1]=A[0][1]*A[1][2]-A[1][1]*A[0][2];
  Ainv[1][0]=A[1][2]*A[2][0]-A[2][2]*A[1][0];
  Ainv[1][1]=A[0][0]*A[2][2]-A[2][0]*A[0][2];
  Ainv[1][2]=A[0][2]*A[1][0]-A[0][0]*A[1][2];
  Ainv[2][0]=A[1][0]*A[2][1]-A[1][1]*A[2][0];
  Ainv[2][1]=A[0][1]*A[2][0]-A[0][0]*A[2][1];
  Ainv[2][2]=A[0][0]*A[1][1]-A[0][1]*A[1][0];
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      Ainv[i][j]=Ainv[i][j]/DET;

}
void SpwDer(double lambda,double *Sl, double *SlP);
void SpwDer_wrap(double x, double *Sl, double *SlP)
{
  SpwDer(x,Sl,SlP);
}
double Spw(double x);
double Spw_wrap(double x)
{
  return Spw(x);
}
void SpwDer(double lambda, double *Sl, double *SlP)
{
  int a,i,j,n;
  double  GinvR[3], B[3],H[3][3];
  double Ginv[3][3],G[3][3],Ainv[3][3],Binv[3][3];
  double R[3];
  for (a=0; a < 3; a++)
    R[a] = COMB[a]-COMA[a];

  for (i=0; i < 3; i++)
    {
      for (j=0; j < 3; j++)
      	{
 	  Ainv[i][j]=0.;
 	  Binv[i][j]=0.;
 	  for (n=0; n < 3; n++)
 	    {
 	      Ainv[i][j]=Ainv[i][j]+saA[n]*saA[n]*RA[n][i]*RA[n][j];
 	      Binv[i][j]=Binv[i][j]+saB[n]*saB[n]*RB[n][i]*RB[n][j];
 	    }
 	  G[i][j]=(1.0 - lambda)*Ainv[i][j]+lambda*Binv[i][j];
	  H[i][j]=(1.0 - lambda)*(1.0-lambda)*Ainv[i][j]-(lambda*lambda)*Binv[i][j];
  	}
    }
  if (USEDIRINV)
    MATINV33(G,Ginv);
  else
    MATINV(G,Ginv);

  for (i=0; i < 3; i++)
    {
      GinvR[i]=0.0;
      for (j=0; j < 3; j++) 
	GinvR[i]=GinvR[i]+Ginv[i][j]*R[j];
    }
  /*     se il valore di Slp passato è maggiore di 0 calcola la
	 derivata di S(lambda) */
  *SlP = 0.0;
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      *SlP += GinvR[i]*H[i][j]*GinvR[j]; 
  *SlP = -*SlP;

  /* calculate also S(lambda) if Sl > 0.0*/
  *Sl = 0.0;
  for (i=0; i < 3; i++)
    *Sl += R[i]*GinvR[i];  
  *Sl = lambda*(1.0-lambda)*(*Sl);
  *Sl = -(*Sl);
}
double Spw(double lambda)
{
  int a,i,j,n;
  double Spw;
  double AB[3][3],ABinv[3][3],Ainv[3][3],Binv[3][3];
  double R[3];
  for (a=0; a < 3; a++)
    R[a] = COMB[a]-COMA[a];
  for (i=0; i < 3; i++)
    {
      for (j=0; j < 3; j++)
	{
	  Ainv[i][j]=0.;
	  Binv[i][j]=0.;
	  for (n=0; n < 3; n++)
	    {
	      Ainv[i][j]=Ainv[i][j]+saA[n]*saA[n]*RA[n][i]*RA[n][j]; 
	      Binv[i][j]=Binv[i][j]+saB[n]*saB[n]*RB[n][i]*RB[n][j];
	    }
	  ABinv[i][j]=(1.0 - lambda)*Ainv[i][j]+lambda*Binv[i][j];
	}
    }
  //      INVERT ABinv(3,3) matrix
    if (USEDIRINV) 
      MATINV33(ABinv,AB);
    else
      MATINV(ABinv,AB);
  Spw=0.0;
  for (i=0; i < 3; i++)
    {
      for (j=0; j < 3; j++)
	{
	  Spw+=R[i]*AB[i][j]*R[j];
	}
    }
  Spw = Spw*lambda*(1.0-lambda);
  return -Spw;
}
int main (int argc, char** argv)
{
  FILE *f;
  int ncall;
  double bret, xmin;
  if (USEDBRENT) 
    printf("Using DBRENT");
  else
    printf("Using BRENT"); 
  if (READFROMFILE) 
    {
      f = fopen("ellips.pos","r"); 
      fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
	     &COMA[0], &COMA[1], &COMA[2], 
	     &RA[0][0], &RA[0][1], &RA[0][2], &RA[1][0], &RA[1][1], &RA[1][2],
	     &RA[2][0], &RA[2][1], &RA[2][2], &saA[0], &saA[1], &saA[2]);
      fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
	     &COMB[0], &COMB[1], &COMB[2], 
	     &RB[0][0], &RB[0][1], &RB[0][2], &RB[1][0], &RB[1][1], &RB[1][2],
	     &RB[2][0], &RB[2][1], &RB[2][2], &saB[0], &saB[1], &saB[2]);
      fclose(f);
      //      it should be the 6th element read 
      printf("RA(1,...)=%.15G %.15G %.15G\n", RA[0][0], RA[0][1], RA[0][2]);
      printf("RA(2,...)=%.15G %.15G %.15G\n", RA[1][0], RA[1][1], RA[1][2]);
      printf("RA(3,...)=%.15G %.15G %.15G\n", RA[2][0], RA[2][1], RA[2][2]); 
      printf("semiaxes=%.15G %.15G %.15G\n",saA[0], saA[1], saA[2]);
    }
  else
    {
      saA[0]=2.;
      saA[1]=2.;
      saA[2]=1.;
      saB[0]=2.;
      saB[1]=2.;
      saB[2]=1.;
      COMA[0]=-1.99;
      COMA[1]=0.;
      COMA[2]=0.;
      COMB[0]=1.99;
      COMB[1]=0.;
      COMB[2]=0.;
      RA[0][0]=1.;
      RA[0][1]=0.;
      RA[0][2]=0.;
      RA[1][0]=0.;
      RA[1][1]=1.;
      RA[1][2]=0.;
      RA[2][0]=0.;
      RA[2][1]=0.;
      RA[2][2]=1.;
      RB[0][0]=1.;
      RB[0][1]=0.;
      RB[0][2]=0.;
      RB[1][0]=0.;
      RB[1][1]=1.;
      RB[1][2]=0.;
      RB[2][0]=0.;
      RB[2][1]=0.;
      RB[2][2]=1.;
    } 
  for (ncall=0; ncall < 10; ncall++)
    {
      if (USEDBRENT) 
	bret=dbrent(LAM0,LAM1,LAM2,SpwDer_wrap,TOL,&xmin);
      else
	bret=brent(LAM0,LAM1,LAM2,Spw_wrap,TOL,&xmin);
    }
  printf("F(A,B)=%.15G\n", -bret); 
}
