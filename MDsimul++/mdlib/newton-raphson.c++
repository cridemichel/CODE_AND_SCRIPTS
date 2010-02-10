void NR(class NR-sys, Vector iniGuess)
{
  x = iniGuess;
}

void ~NR(void)
{

}

Vector<int n> FindZero(class NR-sys)
{
  int i,its=0,ok;
  double d,stpmax,sum,test; 
  NR::fvec = NR-sys.funcs2beZeroed(NR::x);
  test=0.0; /* Test for initial guess being a root. Use more stringent test than simply TOLF.*/

  if (testConvFunc())
    return;

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
  callsdistNR++;
  
  for (its=0;its<MAXITS;its++)
    { 
      /* Start of iteration loop. */
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
      /* If analytic Jacobian is available, you can 
	 replace the routine fdjac below with your own routine.*/
      test = test_func_values(fvecD, n);
      if (test < TOLFD)
	{
	  //printf("[newtDistNeg] test < TOLFD\n");
	  *check = 0;
	  MD_DEBUG(printf(" test < TOLFD\n"));
	  FREERETURND;
	}
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

      if (OprogStatus.toldxNR > 0.0)
	{
	  MD_DEBUG(printf("qui?!?\n"));
	  if (OprogStatus.dist5)
	    adjust_step_dist5(x, p, fx, gx);
	  else
	    adjust_step_dist8(x, p, fx, gx);
	}
      for (i=0;i<n;i++) 
	{ 
	  xold[i] = x[i];
	  x[i] += p[i];
	}
      test = test_xvalues(xold, x, n);
      MD_DEBUG(printf("test = %.15f x = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", test, x[0], x[1], x[2], x[3],x[4]));
      //MD_DEBUG(printf("iA: %d iB: %d test: %f\n",iA, iB,  test));
      if (test < TOLXD) 
	{ 
	  //printf("test<TOLXD\n");
	  *check = 0;
	  MD_DEBUG(printf("test < TOLXD\n"));
	  FREERETURND; 
	}
      itsNRdist++;
    } 
  MD_DEBUG18(printf("maxits!!!\n"));
  *check = 2;
  FREERETURND;
  return;
  nrerror("MAXITS exceeded in newt"); 
}
