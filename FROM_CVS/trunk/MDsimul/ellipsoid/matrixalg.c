#define TINY 1E-20
#define MD_NBMAX 8 
#include<mdsimul.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define MD_DEBUG(X) 
#define MD_DEBUG10(X) 
void nrerror(char *msg)
{
  printf(msg);
  exit(-1);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define ITMAXZB 100 
/* Maximum allowed number of iterations.*/
#define EPSP 3.0e-8 /* Machine floating-point precision.*/
extern int polinterr;
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

double zbrent(double (*func)(double), double x1, double x2, double tol)
/* Using Brent s method, find the root of a function func known to lie between x1 and x2. 
 * The root, returned as zbrent, will be refined until its accuracy is tol.*/
{
  int iter; 
  double a=x1,b=x2,c=x2,d,e,min1,min2; 
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm; 
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    {
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

  polinterr = 1;
  return 0.0;
  //nrerror("Maximum number of iterations exceeded in zbrent"); 
  return 0.0; /* Never get here.*/ 
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
  int i,icol,irow,j,k,l,ll; 
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
#define FMAX(A,B) ((A)>(B)?(A):(B))
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
int *ivector(int n)
{
  return calloc(n, sizeof(int)); 
}
double *vector(int n)
{
  return calloc(n, sizeof(double)); 
}
double **matrix(int n, int m)
{
  double **M;
  int i, j;
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
int nn, nn2, nnD; /* Global variables to communicate with fmin.*/
double *fvec, *fvecG, *fvecD; 
#ifdef MD_GLOBALNR2
#define FREERETURN {MD_DEBUG(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d\n", x[0], x[1], x[2], x[3], x[4], test, its, *check));\
free_vector(fvec);free_vector(xold);free_vector(g2); free_vector(xold2); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);free_vector(fvecG);return;}
#else
#define FREERETURN {MD_DEBUG10(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d fvec=(%.15G,%.15G,%.15G,%.15G,%.15G)\n", x[0], x[1], x[2], x[3], x[4], test, its, *check, fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));\
free_vector(fvec);free_vector(xold); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);free_vector(fvecG);return;}
#endif
#define FREERETURND {if(its>10)printf("calcDistNeg its=%d\n", its); MD_DEBUG(printf("x=(%f,%f,%f,%f,%f) test: %f its: %d check:%d\n", x[0], x[1], x[2], x[3], x[4], test, its, *check));\
free_vector(fvecD);free_vector(xold); free_vector(p); free_vector(g);free_matrix(fjac,n);free_ivector(indx);return;}

void (*nrfuncv)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncv2)(int n, double v[], double fvec[], int i, int j, double shift[3]);
void (*nrfuncvD)(int n, double v[], double fvec[], int i, int j, double shift[3]);


extern void fdjac(int n, double x[], double fvec[], double **fjac, 
		  void (*vecfunc)(int n, double v[], double fvec[], int i, int j, double shift[3]), int iA, int iB, double shift[3]); 
extern void fdjacDist(int n, double x[], double fvec[], double **fjac, 
		  void (*vecfunc)(int n, double v[], double fvec[], int i, int j, double shift[3]), int iA, int iB, double shift[3]); 
double fmin(double x[], int iA, int iB, double shift[3]);
double fminD(double x[], int iA, int iB, double shift[3]);
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], double *f, 
	    double stpmax, int *check, double (*func)(double [], int, int, double []),
	    int iA, int iB, double shift[3], double tolx);
void lubksb(double **a, int n, int *indx, double b[]); 
void ludcmp(double **a, int n, int *indx, double *d, int *ok); 
extern void funcs2beZeroedGuess(int n, double x[], double fvec[], int i, int j, double shift[3]);
extern void funcs2beZeroed(int n, double x[], double fvec[], int i, int j, double shift[3]);

extern void upd2tGuess(int i, int j, double shift[3], double tGuess);
#define MD_GLOBALNR
#undef MD_GLOBALNR2
#ifdef MD_GLOBALNR2
double fmin2(double x[], int iA, int iB, double shift[3]);
#endif
void newt(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3])
{
  int ii, i,its, its2,j,*indx, ok;
  double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold, alphaold; 
  double r1[3], r2[3], alpha;
#ifdef MD_GLOBALNR2
  int check2;
  double f2, stpmax2, fold2, *xold2, *g2;
  xold2 = vector(n);
  g2 = vector(n);
#endif
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvec=vector(n); 
  fvecG=vector(n);
  /*Define global variables.*/
  nn=n; 
  nn2=n-1;
  nrfuncv=vecfunc; 
  f=fmin(x,iA,iB,shift); /*fvec is also computed by this call.*/
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
#ifndef MD_GLOBALNR2
  funcs2beZeroed(n,x,fvec,iA,iB,shift);
#endif
  for (its=0;its<MAXITS;its++)
    { /* Start of iteration loop. */
      /* Stabilization */
#if 0
      MD_DEBUG(printf("its=%d time = %.15f\n",its, x[4]));
      upd2tGuess(iA, iB, shift, x[4]);
#ifdef MD_GLOBALNR2
      nn2=n-1; 
      nrfuncv2=funcs2beZeroedGuess; 
      f2=fmin2(x,iA,iB,shift); /*fvec is also computed by this call.*/
      test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/
      for (i=0;i<n-1;i++) 
	if (fabs(fvec[i]) > test)
	  test=fabs(fvec[i]); 
      for (sum=0.0,i=0;i<n-1;i++) 
	sum += Sqr(x[i]); /* Calculate stpmax for line searches.*/
      stpmax2=STPMX*FMAX(sqrt(sum),(double)(n-1));
#else
      funcs2beZeroedGuess(n-1,x,fvecG,iA,iB,shift);
#endif
      for (its2=0; its2 < MAXITS2 ; its2++)
	{
#ifdef MD_GLOBALNR2
	  if (its2=0 && test < 0.01*TOLF)
	    {
	      check2=0; 
	      break;
	    }
#endif
	  MD_DEBUG(printf("Guessing its2=%d x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", its2, x[0], x[1], x[2], x[3],x[4]));
	  fdjacGuess(n-1,x,fvecG,fjac,funcs2beZeroedGuess, iA, iB, shift); 
#ifdef MD_GLOBALNR2
	  for (i=0;i<n-1;i++) 
	    { /* Compute  f for the line search.*/
	      for (sum=0.0,j=0;j<n-1;j++)
		sum += fjac[j][i]*fvecG[j]; 
	      g2[i]=sum; 
	    } 
	  for (i=0;i<n-1;i++) 
	    xold2[i]=x[i]; /* Store x,*/ 
	  fold2=f2; /* and f. */
#else
	  MD_DEBUG(printf("Guessing fvecG = (%.15f, %.15f, %.15f, %.15f)\n", fvecG[0], fvecG[1], fvecG[2], fvecG[3]));
	  test=0.0; /* Test for convergence on function values.*/
	  for (i=0;i<(n-1);i++) 
	    {
	      test +=fabs(fvecG[i]); 
	    }
	  if (test < TOLF2)
	    {
	      MD_DEBUG(printf("test=%.15f GUESS test < TOLF\n", test));
	      break;
	    }
#endif
	  for (i=0;i<n-1;i++) 
	    p[i] = -fvecG[i]; /* Right-hand side for linear equations.*/
	  ludcmp(fjac,n-1,indx,&d); /* Solve linear equations by LU decomposition.*/
	  lubksb(fjac,n-1,indx,p);
#ifdef MD_GLOBALNR2
	  lnsrch(n-1,xold2,fold2,g2,p,x,&f2,stpmax2,&check2,fmin2,iA,iB,shift); 
	  test=0.0; /* Test for convergence on function values.*/
	  for (i=0;i<n;i++) 
	    if (fabs(fvecG[i]) > test) 
	      test=fabs(fvecG[i]); 
	  if (test < TOLF2) 
	    { 
	      check2=0; 
	      MD_DEBUG(printf("test < TOLF\n"));
	      break;
	    }
	  if (check2) 
	    { /* Check for gradient of f zero, i.e., spurious convergence.*/
	      test=0.0; 
	      den=FMAX(f,0.5*(n-1));
	      for (i=0;i<n;i++)
		{
		  temp=fabs(g2[i])*FMAX(fabs(x[i]),1.0)/den;
		  if (temp > test) 
		    test=temp; 
		} 
	      check2=(test < TOLMIN ? 2 : 0);
	      MD_DEBUG(printf("*check2:%d test=%f\n", check2, test));
	      break;
	    } 
	  test=0.0; /* Test for convergence on x. */
	  for (i=0;i<n-1;i++) 
	    {
	      temp=(fabs(x[i]-xold2[i]))/FMAX(fabs(x[i]),1.0); 
	      if (temp > test) 
		test=temp; 
	    } 
	  if (test < TOLX2) 
	    {
	      check2 = 0;
	      MD_DEBUG(printf("test<TOLX test=%.15f\n", test));
	      break;
	    }
	  if (check2==2)
	    {
	      MD_DEBUG(printf("spurious convergence\n"));
	      break;
	    }
#else
	  test=0.0; /* Test for convergence on function values.*/
	  for (i=0;i<n-1;i++) 
	    { 
	      test += fabs(p[i]);
	      x[i] += p[i];
	    }
	  //MD_DEBUG(printf("iA: %d iB: %d test: %.15f\n",iA, iB,  test));
	  MD_DEBUG(printf("GUESSED x = (%.15f, %.15f, %.15f, %.15f\n", x[0], x[1], x[2], x[3]));
	  if (test < TOLX2) 
	    { 
	      MD_DEBUG(printf("GUESS test < TOLX\n"));
	      break;
	    }
#endif
	}
#ifdef MD_GLOBALNR2
      if (its2 == MAXITS2 || check2==2)
	{
	  *check = 2;
	  MD_DEBUG(printf("MAXITS2!!\n"));
	  FREERETURN;
	}
#else
       if (its2 == MAXITS2)
	{
	  *check = 2;
	  MD_DEBUG(printf("MAXITS2!!\n"));
	  FREERETURN;
	}

#endif
#endif
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
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvec[i]); 
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
      lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin,iA,iB,shift, TOLX); 
      MD_DEBUG(printf("check=%d test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n",*check, test, x[0], x[1], x[2], x[3],x[4]));
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
	{ /* Check for gradient of f zero, i.e., spurious convergence.*/
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMIN ? 2 : 0);
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
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLX) 
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
#undef MD_GLOBALNRD
#define MAXITS3 200
void newtDistNeg(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3])
{
  int ii, i,its, its2,j,*indx, ok;
  double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold, alphaold; 
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvecD=vector(n); 
  /*Define global variables.*/
  nnD=n; 
  nrfuncvD=vecfunc; 
#ifdef MD_GLOBALNRD
  f=fminD(x,iA,iB,shift); /*fvec is also computed by this call.*/
#else
#ifdef MD_DIST5
  funcs2beZeroedDistNeg5(n,x,fvecD,iA,iB,shift);
#else
  funcs2beZeroedDistNeg(n,x,fvecD,iA,iB,shift);
#endif
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
#ifdef MD_DIST5
      fdjacDistNeg5(n,x,fvecD,fjac,vecfunc, iA, iB, shift);
#else
      fdjacDistNeg(n,x,fvecD,fjac,vecfunc, iA, iB, shift);
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
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
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
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMIN ? 2 : 0);
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
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
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

void newtDist(double x[], int n, int *check, 
	  void (*vecfunc)(int, double [], double [], int, int, double []),
	  int iA, int iB, double shift[3])
{
  int ii, i,its, its2,j,*indx, ok;
  double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold, alphaold; 
  indx=ivector(n); 
  fjac=matrix(n, n);
  g=vector(n);
  p=vector(n); 
  xold=vector(n); 
  fvecD=vector(n); 
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
      test=0.0; /* Test for convergence on function values.*/
      for (i=0;i<n;i++) 
	test +=fabs(fvecD[i]); 
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
	  test=0.0; 
	  den=FMAX(f,0.5*n);
	  for (i=0;i<n;i++)
	    {
	      temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		test=temp; 
	    } 
	  *check=(test < TOLMIN ? 2 : 0);
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
      test = 0;
      for (i=0;i<n;i++) 
	{ 
      	  test += fabs(p[i]);
	  x[i] += p[i];
	}
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
double fmin(double x[], int iA, int iB, double shift[3]) 
/* Returns f = 1 2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
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
#ifdef MD_GLOBALNR2
double fmin2(double x[], int iA, int iB, double shift[3]) 
/* Returns f = 1 2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the 
calling program. Global variables also communicate the function values back to 
the calling program.*/
{
  int i;
  double sum;
  (*nrfuncv2)(nn2,x,fvecG,iA,iB,shift);
  for (sum=0.0,i=0;i<nn2;i++)
    sum += Sqr(fvecG[i]); 
    return 0.5*sum; 
}
#endif
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
