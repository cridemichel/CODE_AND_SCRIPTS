#if defined(POLYELLIPS) || defined(HE_FAST_DIST)
#include<mdsimul.h>
extern double scalProd(double *a, double *b);
// rotazione intorno ad asse (ox,oy,oz) di un angolo theta
// (serve per l'algoritmo per trovare alpha risolvendo una quartica)
/* perram wertheim overlap ellissoidi */
#include<complex.h>
#include<float.h>
#define HE_STOP_CAMERON// or HE_STOP_BINI
//#define HE_USE_HALLEY
extern double ***R;
#define USE_LAPACK
extern void solve_quadratic(double coeff[], int *ns, double sol[2]);
void tRDiagRqe(int i, double M[3][3], double D[3], double Ri[3][3])
{
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = D[0];
  Di[1][1] = D[1];
  Di[2][2] = D[2];
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
/* calculate matrix RM associated to rotation around axis oax by an angle theta */
void calc_rotation_matrix(double oax[3], double theta, double RM[3][3])
{
  int k1, k2, k3;
  double ox, oy, oz, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];

  ox = oax[0];
  oy = oax[1];
  oz = oax[2];
  //NOTA 15/01/2018
  //se uso ranf altero la sequenza casuale e perdo il confronto la simulazione da 50x10^6 già fatta con il metodo di Alberto 
  //theta = (ranf()>0.5?1.:-1.)*M_PI/4.0;
  theta = M_PI/4.0;
  thetaSq=Sqr(theta);
  sinw = sin(theta);
  cosw = (1.0 - cos(theta));
  Omega[0][0] = 0;
  Omega[0][1] = -oz;
  Omega[0][2] = oy;
  Omega[1][0] = oz;
  Omega[1][1] = 0;
  Omega[1][2] = -ox;
  Omega[2][0] = -oy;
  Omega[2][1] = ox;
  Omega[2][2] = 0;
  OmegaSq[0][0] = -Sqr(oy) - Sqr(oz);
  OmegaSq[0][1] = ox*oy;
  OmegaSq[0][2] = ox*oz;
  OmegaSq[1][0] = ox*oy;
  OmegaSq[1][1] = -Sqr(ox) - Sqr(oz);
  OmegaSq[1][2] = oy*oz;
  OmegaSq[2][0] = ox*oz;
  OmegaSq[2][1] = oy*oz;
  OmegaSq[2][2] = -Sqr(ox) - Sqr(oy);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
          /* RM = I - sinw*Omega + cosw*Omega*Omega */ 
	  RM[k1][k2] = (k1==k2)?1.0:0.0 - sinw*Omega[k1][k2] + cosw*OmegaSq[k1][k2];
	}
    }
#if 0
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = Rin[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += Rin[k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     RR[k1][k2] = Ro[k1][k2];
#endif
}
void invM(double m[3][3], double invM[3][3])
{
  double m11m20, m12m20, m12m21, m10m22, m10m21, m11m22, A1, A2, A3;
  double invd;

  m11m20= m[1][1]*m[2][0];
  m12m20= m[1][2]*m[2][0];
  m10m21= m[1][0]*m[2][1];
  m12m21= m[1][2]*m[2][1];
  m10m22= m[1][0]*m[2][2];
  m11m22= m[1][1]*m[2][2];
  A1 = -m[0][2]*m11m20 + m[0][1]*m12m20;
  A2 = m[0][2]*m10m21 - m[0][0]*m12m21;
  A3 = - m[0][1]*m10m22 + m[0][0]*m11m22;

  invd = A1 + A2 + A3;
  invd = 1.0/invd;
  //invd = 1.0/(-m[0][2]*m11m20 + m[0][1]*m12m20 + m[0][2]*m10m21
  //	  - m[0][0]*m12m21 - m[0][1]*m10m22 + m[0][0]*m11m22);

  invM[0][0] = -m12m21 + m11m22;
  invM[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
  invM[0][2] = -m[0][2]*m[1][1] + m[0][1]*m[1][2];
  invM[1][0] = m12m20 - m10m22;
  invM[1][1] = -m[0][2]*m[2][0] + m[0][0]*m[2][2]; 
  invM[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2]; 
  invM[2][0] = -m11m20 + m10m21; 
  invM[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
  invM[2][2] = -m[0][1]*m[1][0] + m[0][0]*m[1][1];

  invM[0][0] *= invd;
  invM[0][1] *= invd;
  invM[0][2] *= invd;
  invM[1][0] *= invd;
  invM[1][1] *= invd;
  invM[1][2] *= invd;
  invM[2][0] *= invd;
  invM[2][1] *= invd;
  invM[2][2] *= invd;
}
#ifdef USE_LAPACK
/* find eigenvectors and eigenvalues */
void sort_eigenvect(double R[3][3], double EV[3])
{
  /* sort from smallest eigenvalues to bigges */
  int i, j, k;
  double u[3], ev;
  for (i=0; i < 2; i++)
    for (j=i+1; j < 3; j++)
      {
	if (EV[i] > EV[j])
	  {
	    /* swap eigenvectors */
	    for (k=0; k < 3; k++)
	      u[k] = R[i][k];
	    for (k=0; k < 3; k++)
	      R[i][k] = R[j][k];
	    for (k=0; k < 3; k++)
	      R[j][k] = u[k];
	    /* swap eigenvalues */
	    ev = EV[i];
	    EV[i] = EV[j];
	    EV[j] = ev;
	  }
      }
}
#if 0
void calcEigenValVect(double I[3][3], double R[3][3], double EV[3])
{
  int ok, kk;
  double u1[3], u2[3];
  wrap_dsyev_he(I, R, EV, &ok);
  sort_eigenvect(R, EV);

  /* we want a right-handed reference system */
  for (kk=0; kk < 3; kk++)
    u1[kk] = R[0][kk];
  for (kk=0; kk < 3; kk++)
    u2[kk] = R[1][kk];
  vectProdVec(u1, u2, R[2]);
}
#endif
void wrap_dsyev_he(double a[3][3], double b[3][3], double x[3], int *ok)
{
  char JOBZ, UPLO;
  double AT[9], work[10];
  int i, j, c1, c2, c3;
  JOBZ='V';
  UPLO='U';
  extern void dsyev_(char* , char*, int*, double* , int*, double*, double*, int*, int*);
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) AT[j+3*i]=a[j][i];		
    }						
  c1 = 3;
  c2 = 1;
  c3 = 8;
  dsyev_(&JOBZ, &UPLO, &c1, AT, &c1, x, work, &c3, ok);      
  /* il primo autovettore è b[0][], il secondo b[1][] ed il terzo b[2][] */
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) b[i][j]=AT[j+3*i];		
    }	
  if (*ok != 0)
    printf("not ok (%d)\n", *ok);
}
#endif

/* Mout=MA*MB */
void mul_matrix(double MA[3][3], double MB[3][3], double Mout[3][3])
{
  int k1, k2, k3;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Mout[k1][k2] = 0;
	for (k3 = 0; k3 < 3; k3++)
	  Mout[k1][k2] += MA[k1][k3]*MB[k3][k2];
      }
}
// cubic polynomial
// we know that 3 zeros have to be real!
void solve_cubic_analytic(double coeff[4], double sol[3])
{
  /* solve the cubic coeff[3]*x^3 + coeff[2]*x^2 +  coeff[1]*x + coeff[0] = 0
   * according to the method described in Numerical Recipe book */  
  double a, b, c, Q, R, theta, Q3, R2, A, B;
  const double sqrt32=sqrt((double)3.0)/2.0;
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
      printf("here?!?\n");
      exit(-1);
      A = -copysign((double)1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
#if 0
      sol[0] = (A+B) - a/3.0;
      sol[1] = cmplx(-0.5*(A+B)-a/3.0,sqrt32*(A-B));
      sol[2] = conj(sol[1]);
#endif
      //sol[1] = -0.5*(A+B)-a/3.0+cmplx(0,1)*sqrt32*(A-B);
      //sol[2] = -0.5*(A+B)-a/3.0-cmplx(0,1)*sqrt32*(A-B);
    }
}
#if defined(HE_STOP_CAMERON)||defined(HE_STOP_BINI)
double alpha_HE[7];
#endif
extern void oqs_quartic_solver(double coeff[], complex double roots[4]);
double polyalpha(double cmon[7], double x)
{
  // evaluate polynomail via Horner's formula 
  double bn=0.0;
  int i;
  for (i=6; i >= 0; i--)
    {
      bn = cmon[i] + bn*x;
    }
  return bn;
}
double polyalphad(double x, double c[7], double *p, double *p1)
{
  // evaluate polynomail via Horner's formula 
  int i;
  int n=6;
#if defined(HE_STOP_CAMERON)
  double err, abx;
#endif
  *p=c[n]*x+c[n-1];
  *p1=c[n];
#if defined(HE_STOP_CAMERON)
  abx = fabs(x);
  err = alpha_HE[n-2]+abx*(alpha_HE[n-1]+alpha_HE[n]*abx);
#endif
  for(i=n-2;i>=0;i--) {
    *p1=*p+*p1*x;
    *p=c[i]+*p*x;
#if defined(HE_STOP_CAMERON)
    err=abx*err+alpha_HE[i];
#endif
  }
#if defined(HE_STOP_CAMERON)
  return err;
#else
  return 0;
#endif
}
double polyalphadd(double x, double c[7], double *p, double *p1, double *p2)
{
  // evaluate polynomial via Horner's formula 
  int i;
  int n=6;
#if defined(HE_STOP_CAMERON)
  double err, abx;
#endif
  /*c0 + c1 x + x^2 (c2 + x (c3 + x (c4 + x (c5 + c6 x)))) */
  *p=c[n-2]+x*(c[n-1]+c[n]*x);
#if defined(HE_STOP_CAMERON)
  abx = fabs(x);
  err = alpha_HE[n-2]+abx*(alpha_HE[n-1]+alpha_HE[n]*abx);
#endif
  *p1=c[n-1] + 2.0*c[n]*x;
  *p2=2.0*c[n];
  for(i=n-3;i>=0;i--) {
    (*p2) = 2.0*(*p1) + (*p2)*x;
    (*p1) = (*p) + (*p1)*x;
    (*p) = c[i] + (*p)*x;
#if defined(HE_STOP_CAMERON)
    err=abx*err+alpha_HE[i];
#endif
  }
#if defined(HE_STOP_CAMERON)
  return err;
#else 
  return 0;
#endif
}

inline double PowerI(double x, int n)
{
  int i;
  double ret=1.0;
  switch (n)
    {
    case 1:
      ret = x;
      break;
    case 2:
      ret =  x*x;
      break;
    case 3:
      ret = x*x*x;
      break;
    case 4:
      ret = x*x*x*x;
      break;
    default:
      for (i=0; i < n; i++)
        ret *= x;
    }
  //printf("x=%f n=%d ret=%f\n", x, n, ret);
  return ret;

}

double distSq2orig(double alpha, double sx, double sy, double sz, double x0, double y0, double z0)
{
  return x0*x0/Sqr(1.0+sx*sx*alpha) + y0*y0/Sqr(1.0+sy*sy*alpha) + z0*z0/Sqr(1.0+sz*sz*alpha);
}
extern double calc_norm(double *v);
void intersectPoint(double Alpha, double m00, double m01, double m02, double m11, double m12, double m22, 
                    double m012, double m022, double m122, double m222, double x0, double y0, double z0, double sai[3],
                    double pi[3], double pj[3])
{
  double d[3], xi[3], norm;
  double x[3], detMa, Al2, Al3, Al4;
  int i;
  Al2 = Alpha*Alpha;
  Al3 = Al2*Alpha;
  Al4 = Al2*Al2;
  detMa = -m022*m11 + 2*m01*m02*m12 - m00*m12*m12 - m012*m22 + m00*m11*m22 - 
    m012*Alpha - m022*Alpha + m00*m11*Alpha - m122*Alpha +
    m00*m22*Alpha + m11*m22*Alpha + m00*Al2 + 
    m11*Al2 + m22*Al2 + Al3;
  x[0] = -(m022*x0*(m11 + Alpha)) + 
    m02*(2*m01*m12*x0 - m12*y0*Alpha + z0*Alpha*(m11 + Alpha)) - 
    m01*(m12*z0*Alpha + m01*x0*(m22 + Alpha) - y0*Alpha*(m22 + Alpha)) + 
    m00*x0*(-m122 + (m11 + Alpha)*(m22 + Alpha));
  x[1] = -(m022*m11*y0) + 2*m01*m02*m12*y0 - m012*m22*y0 + 
    m01*m22*x0*Alpha - (m012 + m122 - m11*m22)*y0*Alpha - 
    m02*(m12*x0 + m01*z0)*Alpha + (m01*x0 + m11*y0 + m12*z0)*Al2 + 
    m00*(-(m122*y0) + m12*z0*Alpha + m11*y0*(m22 + Alpha));
  x[2] = -(m012*m22*z0) - (m01*m12*x0 + (m122 - m11*m22)*z0)*Alpha + 
   (m12*y0 + m22*z0)*Al2 - m022*z0*(m11 + Alpha) + 
   m00*(-(m122*z0) + m12*y0*Alpha + m22*z0*(m11 + Alpha)) + 
   m02*(x0*Alpha*(m11 + Alpha) + m01*(2*m12*z0 - y0*Alpha));
  for (i=0; i< 3; i++)
    x[i] /= detMa;
  norm = calc_norm(x); 
  for (i=0; i < 3; i++)
    xi[i] = x[i]/norm;
#if 1
  for (i=0; i < 3; i++)
    d[i] = x[i] - xi[i];
  for (i=0; i < 3; i++)
    d[i] *= sai[i];
#endif
  //printf("sai=%f %f %f norm=%f norm(d)=%f\n", sai[0], sai[1], sai[2], norm, calc_norm(d));
  //printf("d=%f %f %f x=%f %f %f xi=%f %f %f\n", d[0], d[1], d[2], x[0], x[1], x[2], xi[0], xi[1], xi[2]);

  for (i=0; i < 3; i++)
    {
      pi[i] = xi[i]*sai[i]; /* questo è il punto di tangenza dei due ellissoidi se traslo il secondo di un vettore pari 
                               alla distanza */
      pj[i] = x[i]*sai[i];
    }
  // ora devo calcolare il punto sull'altro ellissoide
  //return (norm > 1.0)?(calc_norm(d)):(-calc_norm(d));
}
double distSq2origMopt(double Alpha, double m00, double m01, double m02, double m11, double m12, double m22, 
                    double m012, double m022, double m122, double m222, double x0, double y0, double z0)
{
  double x[3], detMa, Al2, Al3, Al4;
  int i;
  Al2 = Alpha*Alpha;
  Al3 = Al2*Alpha;
  Al4 = Al2*Al2;
  detMa = -m022*m11 + 2*m01*m02*m12 - m00*m12*m12 - m012*m22 + m00*m11*m22 - 
    m012*Alpha - m022*Alpha + m00*m11*Alpha - m122*Alpha +
    m00*m22*Alpha + m11*m22*Alpha + m00*Al2 + 
    m11*Al2 + m22*Al2 + Al3;
  x[0] = -(m022*x0*(m11 + Alpha)) + 
    m02*(2*m01*m12*x0 - m12*y0*Alpha + z0*Alpha*(m11 + Alpha)) - 
    m01*(m12*z0*Alpha + m01*x0*(m22 + Alpha) - y0*Alpha*(m22 + Alpha)) + 
    m00*x0*(-m122 + (m11 + Alpha)*(m22 + Alpha));
  x[1] = -(m022*m11*y0) + 2*m01*m02*m12*y0 - m012*m22*y0 + 
    m01*m22*x0*Alpha - (m012 + m122 - m11*m22)*y0*Alpha - 
    m02*(m12*x0 + m01*z0)*Alpha + (m01*x0 + m11*y0 + m12*z0)*Al2 + 
    m00*(-(m122*y0) + m12*z0*Alpha + m11*y0*(m22 + Alpha));
  x[2] = -(m012*m22*z0) - (m01*m12*x0 + (m122 - m11*m22)*z0)*Alpha + 
   (m12*y0 + m22*z0)*Al2 - m022*z0*(m11 + Alpha) + 
   m00*(-(m122*z0) + m12*y0*Alpha + m22*z0*(m11 + Alpha)) + 
   m02*(x0*Alpha*(m11 + Alpha) + m01*(2*m12*z0 - y0*Alpha));
  for (i=0; i< 3; i++)
    x[i] /= detMa;
 return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; 
}
double distSq2origM(double Alpha, double m00, double m01, double m02, double m11, double m12, double m22, 
                    double x0, double y0, double z0)
{
  double x[3], detMa;
  int i;
  detMa = -m02*m02*m11 + 2*m01*m02*m12 - m00*m12*m12 - m01*m01*m22 + m00*m11*m22 - 
    m01*m01*Alpha - m02*m02*Alpha + m00*m11*Alpha - m12*m12*Alpha +
    m00*m22*Alpha + m11*m22*Alpha + m00*Alpha*Alpha + 
    m11*Alpha*Alpha + m22*Alpha*Alpha + Alpha*Alpha*Alpha;
  x[0] = -(PowerI(m02,2)*x0*(m11 + Alpha)) + 
    m02*(2*m01*m12*x0 - m12*y0*Alpha + z0*Alpha*(m11 + Alpha)) - 
    m01*(m12*z0*Alpha + m01*x0*(m22 + Alpha) - y0*Alpha*(m22 + Alpha)) + 
    m00*x0*(-PowerI(m12,2) + (m11 + Alpha)*(m22 + Alpha));
  x[1] = -(PowerI(m02,2)*m11*y0) + 2*m01*m02*m12*y0 - PowerI(m01,2)*m22*y0 + 
    m01*m22*x0*Alpha - (PowerI(m01,2) + PowerI(m12,2) - m11*m22)*y0*Alpha - 
    m02*(m12*x0 + m01*z0)*Alpha + (m01*x0 + m11*y0 + m12*z0)*PowerI(Alpha,2) + 
    m00*(-(PowerI(m12,2)*y0) + m12*z0*Alpha + m11*y0*(m22 + Alpha));
  x[2] = -(PowerI(m01,2)*m22*z0) - (m01*m12*x0 + (PowerI(m12,2) - m11*m22)*z0)*Alpha + 
   (m12*y0 + m22*z0)*PowerI(Alpha,2) - PowerI(m02,2)*z0*(m11 + Alpha) + 
   m00*(-(PowerI(m12,2)*z0) + m12*y0*Alpha + m22*z0*(m11 + Alpha)) + 
   m02*(x0*Alpha*(m11 + Alpha) + m01*(2*m12*z0 - y0*Alpha));
  for (i=0; i< 3; i++)
    x[i] /= detMa;
 return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; 
}

double calcfel(double M[3][3], double r0[3], double x[3])
{
  int i, j;
  double res, v[3], xr0[3];
  
  for (i=0; i < 3; i++)
    xr0[i] = x[i] - r0[i];
  for (i=0; i < 3; i++)
    {
      v[i]=0;
      for (j=0; j < 3; j++)
        {
          v[i] += M[i][j]*xr0[j];
        }
    }
  res=0.0;
  for (i=0; i < 3; i++)
    {
      res += xr0[i]*v[i];
    }
  return res-1.0;
}
extern double zbrent(double (*func)(double), double x1, double x2, double tol);
double coeffpa_HE[7];
double rtsafe(double c[7], double xg, double x1, double x2, double  xacc, int guess, int *ok)
{
  /* xg is the initial guess and x1, x2 must bracket the solution */
  const int MAXIT=100; //Maximum allowed number of iterations. Doub xh,xl;
  const double Kconv=3.8;
  double fl, fh, xl, xh;
  double err, rts, dx, dxold, df, f, temp, absp;
#ifdef HE_USE_HALLEY
  double ddf, corr;
#endif
#if defined(HE_STOP_BINI)
  double abx, s;
#endif
  int j, j2;
  fl=polyalpha(c,x1);
  fh=polyalpha(c,x2);
  *ok=1;
#if defined(HE_STOP_CAMERON)
  for (j=0; j <= 6; j++)
    {
      alpha_HE[j] = fabs(c[j])*Kconv*(j+1.0);
    }
#elif defined(HE_STOP_BINI)
  for (j=0; j <= 6; j++)
    {
      alpha_HE[j] = fabs(c[j]);
    }
#endif
#if 1
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) 
    {
      printf("[WARNING] Root must be bracketed in rtsafe\n");
      printf("fl=%.15G fh=%.15G\n", fl, fh);
      printf("xg=%.15G x1=%.15G x2=%.15G\n", xg, x1, x2);
      return xg;
    }
#endif
  if (fl == 0.0) 
    return x1;
  if (fh == 0.0) 
    return x2;
  if (fl < 0.0) 
    {
      xl=x1;
      xh=x2;
    } 
  else 
    {
      xh=x1;
      xl=x2; 
    }
  if (guess)
    rts = xg;
  else
    rts = 0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
#ifdef HE_USE_HALLEY
  err=polyalphadd(rts,c,&f,&df,&ddf);
#else
  err=polyalphad(rts,c,&f,&df);
#endif
  for (j=0;j<MAXIT;j++) 
    {
#if defined(HE_STOP_CAMERON)
      absp=fabs(f);
      if (absp <= DBL_EPSILON*err)// || (DBL_EPSILON*err <= minf && absp < minf))
        {
          return rts;
        }
#elif defined(HE_STOP_BINI)
      abx = fabs(rts);
      s=alpha_HE[6];
      for (j2=5; j2 >=0; j2--) 
        {
          s=abx*s+alpha_HE[j2];
        }
      if (fabs(f) <= 2.0*DBL_EPSILON*(4.0*6.0+1.0)*s) // stopping criterion of bini 
        {
          return rts;
        }
#endif 
      //Orient the search so that f .xl/ < 0.
      //Initialize the guess for root, the “stepsize before last,” and the last step.
      //Loop over allowed iterations.
      if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) 
	  || (fabs(2.0*f) > fabs(dxold*df)) || (df==0)) 
	{ 
	  //Bisect if Newton out of range, 
	  //or not decreasing fast enough.
	  dxold=dx;
	  dx=0.5*(xh-xl);
	  rts=xl+dx;
	  if (xl == rts) 
            {
              //printf("qui\n");
              return rts;
            }
	} 
      else 
	{
	  dxold=dx;
#ifdef HE_USE_HALLEY
          /* cubic convergence but we need to evaluate second derivative too */       
          corr = 1.0-f*ddf/(2.0*df*df);
          if (corr > 1.2)
            corr=1.2;
          if (corr < 0.8)
            corr=0.8;
          dx = f/df/corr;
          //dx = f/(df*(1.0-f*ddf/(2.0*df*df)));
#else
	  dx= f/df;
#endif
          temp=rts;
	  rts -= dx;
	  if (temp == rts) 
            {
              return rts;
            }
	}
#if !defined(HE_STOP_CAMERON) && !defined(HE_STOP_BINI)
      if (fabs(dx) < xacc) 
        {
          return rts;
        }
#endif
#ifdef HE_USE_HALLEY
      err=polyalphadd(rts,c,&f,&df,&ddf);
#else
      err=polyalphad(rts, c, &f, &df);
#endif

      //The one new function evaluation per iteration.
      if (f < 0.0) //Maintain the bracket on the root.
	xl=rts; 
      else
	xh=rts;
    }

  printf("[WARNING] Maximum number of iterations exceeded in rtsafe\n");
  *ok=0;
  return rts;
}
//#define MC_ALPHAGUESS
double pi_HE[3], pj_HE[3];
double check_overlap_polyell(int i, int j, double shift[3])
{
  const double GOLD=1.618034;
  double Rjp[3][3], r0jp[3], ri[3], rj[3], Di[3], Dj[3], Mjp[3][3], gradj[3], xppg[3], M2I[3][3], invM2I[3][3], oax[3], theta, norm;
  double RM[3][3], Mtmp[3][3], Mjpp[3][3], r0jpp[3], Mjp3[3][3], r0jp3[3], xp3g[3];
  double Mip[3][3], Mipp[3][3];
  double evec[3][3], eval[3], coeff[4], coeffpa[7], sx, sy, sz, x0, y0, z0;
  double sx2, sy2, sz2, sx4, sy4, sz4, x02, y02, z02, sai[3];
  double vt, at[3], Mi[3][3], Mj[3][3], Ri[3][3], Rj[3][3], dist, m00, m01, m02, m11, m22, m12, b, docc[3], doccn;
  double aR, aL, dR, dL, alpha, alg;
  int ig, it;
#ifdef MC_ALPHAGUESS
  double ng, gg[3], alq, solq[2], coeffq[3], xx[3];
  int numsol;
#endif
  int np, typei, typej;
  int kk1, kk2, kk3, k1, k2, k3, ok;
  complex double roots[6];
#define HE_OPT
#ifdef HE_OPT
  double m022, m012, m122, m002, m222, m112, m014, m024, m124, m123, m003, m015, m013, m023;
  double m113, m223;
#endif
#ifdef HE_FAST_DIST
  double gi[3], Mg[3], pmx0[3], gMg, Mp[3], pMp, gMp, cq_fd[3], solq_fd[2], lam_fd, dpij[3];
  int numsol_fd;
#endif 
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  
  /* apply affinity to reduce first ellipsoid to a sphere */
  for (k1=0; k1 < 3; k1++)
    {
      Di[k1]= 1.0/Sqr(typesArr[typei].sax[k1]);
      Dj[k1]= 1.0/Sqr(typesArr[typej].sax[k1]);
      for (k2=0; k2 < 3; k2++)
        {
          Ri[k1][k2] = R[i][k1][k2];
          Rj[k1][k2] = R[j][k1][k2];
        }
    }
  for (k1=0; k1 < 3; k1++)
    {
      sai[k1] = typesArr[typei].sax[k1];
    } 
  tRDiagRqe(i, Mi, Di, Ri);
  tRDiagRqe(j, Mj, Dj, Rj);
  ri[0] = rx[i];
  ri[1] = ry[i];
  ri[2] = rz[i];
  rj[0] = rx[j]+shift[0];
  rj[1] = ry[j]+shift[1];
  rj[2] = rz[j]+shift[2];

  /* verifico che il centro di i non appartenga a j e viceversa come check preliminare */
#if !defined(HE_FAST_DIST)
  if (calcfel(Mi,ri,rj) < 0.0)
    return -1.0;
  if (calcfel(Mj,rj,ri) < 0.0)
    return -1.0;
#endif
  //printf("coords i: %f %f %f j: %f %f %f\n", ri[0], ri[1], ri[2], rj[0], rj[1], rj[2]);
  /* switch to ellipsoid i reference system */
  for (kk1=0; kk1 < 3; kk1++)
    {
      r0jp[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
        {
          r0jp[kk1] += R[i][kk1][kk2]*(rj[kk2]-ri[kk2]);
          Rjp[kk1][kk2] = 0;
          //Aip[kk1] = 0;
          for (kk3=0; kk3 < 3; kk3++)
            {
              Rjp[kk1][kk2] += R[j][kk1][kk3]*R[i][kk2][kk3];
              //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
            } 
        }
    }
  //tRDiagRpw(i, Mi, DA, R[i]);
  tRDiagRqe(j, Mjp, Dj, Rjp);

  /* calculate matrix and position of ellipsoid j after application of affinity
   * which reduces ellipsoid i to a sphere */   
  Mjpp[0][0] = Mjp[0][0]*Sqr(sai[0]);
  Mjpp[0][1] = Mjp[0][1]*sai[0]*sai[1];
  Mjpp[0][2] = Mjp[0][2]*sai[0]*sai[2];
  Mjpp[1][0] = Mjp[1][0]*sai[0]*sai[1];
  Mjpp[1][1] = Mjp[1][1]*Sqr(sai[1]);
  Mjpp[1][2] = Mjp[1][2]*sai[1]*sai[2];
  Mjpp[2][0] = Mjp[2][0]*sai[0]*sai[2];
  Mjpp[2][1] = Mjp[2][1]*sai[1]*sai[2];
  Mjpp[2][2] = Mjp[2][2]*Sqr(sai[2]);
  r0jpp[0] = r0jp[0]/sai[0];
  r0jpp[1] = r0jp[1]/sai[1];
  r0jpp[2] = r0jp[2]/sai[2];
  x0 = r0jpp[0];
  y0 = r0jpp[1];
  z0 = r0jpp[2];
  m00 = Mjpp[0][0];
  m01 = Mjpp[0][1];
  m02 = Mjpp[0][2];
  m11 = Mjpp[1][1];
  m12 = Mjpp[1][2];
  m22 = Mjpp[2][2]; 
#ifdef HE_OPT
  x02  = Sqr(x0);
  y02  = Sqr(y0);
  z02  = Sqr(z0);
  m022 = Sqr(m02);
  m012 = Sqr(m01);
  m122 = Sqr(m12);
  m002 = Sqr(m00);
  m222 = Sqr(m22);
  m112 = Sqr(m11); 
  m014 = Sqr(m012);
  m024 = Sqr(m022);
  m124 = Sqr(m122); 
  m123 = m122*m12;
  m003 = m002*m00;
  m015 = m014*m01;
  m013 = m012*m01;
  m023 = m022*m02;
  m113 = m112*m11;
  m223 = m222*m22;
  
  b = -1 + m00*x02 + 2*m01*x0*y0 + m11*y02 + 2*m02*x0*z0 + 2*m12*y0*z0 + m22*z02;
  coeffpa[0] = -(PowerI(m022*m11 - 2*m01*m02*m12 + m012*m22 + 
       m00*(m122 - m11*m22),2)* (-b + m00*x02 + 2*m01*x0*y0 + m11*y02 + 2*m02*x0*z0 + 
       2*m12*y0*z0 + m22*z02));
  coeffpa[1] = -2*(m012 + m022 - m00*m11 + m122 - 
     (m00 + m11)*m22)*(m022*m11 - 2*m01*m02*m12 + 
     m012*m22 + m00*(m122 - m11*m22))*
   (-b + m00*x02 + 2*m01*x0*y0 + m11*y02 + 2*m02*x0*z0 + 
     2*m12*y0*z0 + m22*z02);
  coeffpa[2] = b*(m014 + m024 + m002*m112 - 2*m002*m122 - 4*m00*m11*m122 + 
      m124 + 4*m00*m11*(m00 + m11)*m22 - 2*(2*m00 + m11)*m122*m22 + 
      (m002 + 4*m00*m11 + m112)*m222 + 4*m01*m02*m12*(m00 + m11 + m22) - 
      2*m022*(2*m00*m11 + m112 - m122 + m00*m22 + 2*m11*m22) + 2*m012*
       (m022 - m00*m11 + m122 - 2*(m00 + m11)*m22 - m222)) + m022*m11*m122*x02 - 
   2*m01*m02*m123*x02 - m022*m112*m22*x02 + 2*m01*m02*m11*m12*m22*x02 + m012*m122*m22*x02 - 
   m012*m11*m222*x02 - m003*(m112 - 2*m122 + 4*m11*m22 + m222)* x02 - 2*m015*x0*y0 - 4*m013*m022*x0*y0 - 2*m01*m024*x0*y0 + 
   4*m01*m022*m112*x0*y0 - 8*m012*m02*m11*m12*x0*y0 - 2*m023*m11*m12*x0*y0 - 
   4*m013*m122*x0*y0 - 2*m01*m124*x0*y0 + 8*m013*m11*m22*x0*y0 + 10*m01*m022*m11*m22*x0*y0 - 
   14*m012*m02*m12*m22*x0*y0 + 4*m01*m11*m122*m22*x0*y0 + 6*m013*m222*x0*y0 - 2*m01*m112*m222*x0*y0 - m014*m11*y02 - 
   2*m012*m022*m11*y02 + 2*m022*m113*y02 - 2*m01*m023*m12*y02 - 4*m01*m02*m112*m12*y02 - 
   2*m012*m11*m122*y02 - 2*m022*m11*m122*y02 - m11*m124*y02 + m012*m022*m22*y02 + 4*m012*m112*m22*y02 + 
   4*m022*m112*m22*y02 - 4*m01*m02*m11*m12*m22*y02 + 2*m112*m122*m22*y02 + 2*m012*m11*m222*y02 - 
   m113*m222*y02 - 2*(m02*(m024 + m022* (-3*m112 + 2*m122 - 4*m11*m22) + Sqr(m122 - m11*m22))*x0 + 
      m12*(m024 + Sqr(m122 - m11*m22) - 2*m022*(m112 - m122 + 2*m11*m22))*y0 + m013*m22*(m12*x0 + m02*y0) + 
      m014*(m02*x0 + m12*y0) + m01*m02*(m02*m12*(7*m11 + 4*m22)*x0 + m022*m11*y0 + 4*m122*(m11 + m22)*y0) + 
      m012*(2*m023*x0 - m02*m22*(5*m11 + 2*m22)*x0 + 2*m12*(m122 - m22*(2*m11 + m22))*y0))*z0 - 
   (2*m013*m02*m12 + 4*m01*m02*m12*m22*(m11 + m22) + m22*(m024 + Sqr(m122 - m11*m22) - 
         2*m022*(m112 - m122 + 2*m11*m22)) - m012*(m022*(m11 - 2*m22) + 2*m22*(-m122 + 2*m11*m22 + m222)))*z02 + 
   m002*(4*m11*m122*x02 - 4*m112*m22*x02 + 4*m122*m22*x02 - 4*m11*m222*x02 + 2*m022*(2*m11 + m22)*x02 + 
      2*m012*(m11 + 2*m22)*x02 - m113*y02 + 2*m11*m122*y02 - 4*m112*m22*y02 - m122*m22*y02 - 2*m01*x0*(2*m02*m12*x0 + 
         (m112 - 2*m122 + 4*m11*m22 + m222)*y0) - 2*m02*(m112 - 2*m122 + 4*m11*m22 + m222)*x0*
       z0 - 2*m12*(m112 - 3*m122 + 5*m11*m22 + m222)*y0*z0 - (m11*m122 - 2*m122*m22 + 4*m11*m222 + 
         m223)*z02) + m00*(-(m014*x02) - m024*x02 + 4*m013*(m11 + 2*m22)*x0*y0 + 4*m023*(2*m11 + m22)*x0*z0 - 
      2*m02*(m122 - m11*m22)*x0*(m12*y0 - 5*m11*z0 - 4*m22*z0) - 
      4*(m11 + m22)*(-m122 + m11*m22)* (m11*y02 + z0*(2*m12*y0 + m22*z0)) + m012*(-2*m022*x02 + 4*m11*m22*x02 + 
         2*m222*x02 + 2*m112*y02 + 4*m11*m22*y02 - m222*y02 + 2*m12*(2*m11 + 5*m22)*y0*z0 + 4*m222*z02 + 
         4*m02*x0*(-2*m12*y0 + (m11 + 2*m22)*z0) + m122*(-2*x02 + z02)) + m022*(m122*(-2*x02 + y02) + 
         10*m11*m12*y0*z0 + 4*m12*m22*y0*z0 + 2*m222*z02 + m112*(2*x02 + 4*y02 - z02) + 4*m11*m22*(x02 + z02)) + 
      2*m01*(2*m022*x0*(2*m11*y0 + m22*y0 - 2*m12*z0) - (m122 - m11*m22)*x0*(-4*m11*y0 - 5*m22*y0 + m12*z0) + 
         m02*(-7*m122*y0*z0 + m11*m22*y0*z0 + m12*m22*(-2*x02 + y02 - 2*z02) + m11*m12*(-2*(x02 + y02) + z02))));
  coeffpa[3] = 2*(b*(2*m01*m02*m12 - m11*(2*m022 + m122) - (m022 - m112 + m122)*m22 + 
        m11*m222 + m002*(m11 + m22) - m012*(m11 + 2*m22) + m00*(-m012 - m022 + m112 - 
           2*m122 + 4*m11*m22 + m222)) - m022*m112*x02 + 2*m01*m02*m11*m12*x02 - m012*m11*m22*x02 - 
     m022*m11*m22*x02 + 2*m01*m02*m12*m22*x02 - m012*m222*x02 - m003*(m11 + m22)*x02 + 2*m013*m11*x0*y0 + 
     6*m01*m022*m11*x0*y0 - 8*m012*m02*m12*x0*y0 + 2*m01*m11*m122*x0*y0 + 6*m013*m22*x0*y0 + 
     2*m01*m022*m22*x0*y0 - 2*m01*m112*m22*x0*y0 + 2*m01*m122*m22*x0*y0 - 2*m01*m11*m222*x0*y0 + 
     m012*m112*y02 + 2*m022*m112*y02 - 2*m01*m02*m11*m12*y02 + m112*m122*y02 + 2*m012*m11*m22*y02 - m113*m22*y02 + 
     2*m01*m02*m12*m22*y02 + m11*m122*m22*y02 - m012*m222*y02 - m112*m222*y02 + 2*(-4*m01*m02*m12 + m022*(3*m11 + m22) + 
        m012*(m11 + 3*m22) - (m11 + m22)*(-m122 + m11*m22)) *(m02*x0 + m12*y0)*z0 + (2*m01*m02*m12*(m11 - m22) + 
        m022*(-m112 + 2*m11*m22 + m222) + m22*(2*m012*m22 - (m11 + m22)*(-m122 + m11*m22)))* z02 + m002*
      (m012*x02 + m022*x02 - m112*x02 + 2*m122*x02 - 4*m11*m22*x02 - m222*x02 - 2*m01*(m11 + m22)*x0*y0 - m112*y02 - 
        m122*y02 - 2*(m11 + m22)*(m02*x0 + m12*y0)*z0 - (m122 + m222)*z02) + m00*(2*m022*m11*x02 + m022*m22*x02 + 
        2*m013*x0*y0 - m113*y02 + 2*m11*m122*y02 - 4*m112*m22*y02 - m122*m22*y02 + 
        2*(m022 - m112 + 3*m122 - 5*m11*m22 - m222)*(m02*x0 + m12*y0)*z0 - (m11*m122 + m022*(m11 - m22) - 
           2*m122*m22 + 4*m11*m222 + m223)* z02 + m012* (2*m22*x02 - m22*y02 + 
           m11*(x02 + y02) + 2*m02*x0*z0 + 2*m12*y0*z0) + 2*m01*(m022*x0*y0 - 
           (m112 - 3*m122 + 5*m11*m22 + m222)*x0* y0 + m02*m12*(-x02 + y02 + z02))));
  coeffpa[4] = b*(m002 + m112 - 2*(m012 + m022 + m122) + 4*m11*m22 + 
      m222 + 4*m00*(m11 + m22)) - m003*x02 + 2*m00*m022*x02 - 4*m002*m11*x02 - 
   m012*m11*x02 - 4*m022*m11*x02 + 6*m01*m02*m12*x02 - 4*m002*m22*x02 - 4*m012*m22*x02 - m022*m22*x02 - 
   2*m002*m01*x0*y0 + 6*m013*x0*y0 + 6*m01*m022*x0*y0 - 10*m00*m01*m11*x0*y0 - 
   2*m01*m112*x0*y0 - 2*m00*m02*m12*x0*y0 - 2*m02*m11*m12*x0*y0 + 
   6*m01*m122*x0*y0 - 8*m00*m01*m22*x0*y0 - 8*m01*m11*m22*x0*y0 - 2*m02*m12*m22*x0*y0 + 2*m012*m11*y02 - 
   4*m00*m112*y02 - m113*y02 + 6*m01*m02*m12*y02 - 4*m00*m122*y02 + 2*m11*m122*y02 - 4*m012*m22*y02 - 
   4*m112*m22*y02 - m122*m22*y02 + m00*m012*(2*x02 - y02) - 2*((m002*m02 - 3*m012*m02 + m01*m12*(m11 + m22) + 
         m00*(4*m02*m11 + m01*m12 + 5*m02*m22) + m02*(-3*(m022 + m122) + 4*m11*m22 + m222))
        *x0 + (m00*m01*m02 - 3*m012*m12 + m01*m02*(m11 + m22) + 4*m00*m12*(m11 + m22) + 
         m12*(-3*m022 + m112 - 3*m122 + 5*m11*m22 + m222))*y0)*z0 - (m00*m022 + 4*m022*m11 - 6*m01*m02*m12 + 
      4*m00*m122 + m11*m122 - 2*(m022 + m122)*m22 + 4*(m00 + m11)*m222 + m223)*z02;
  coeffpa[5] = 2*b*(m00 + m11 + m22) - 2*(m002 + m012 + m022)* x02 - 4*(m01*(m00 + m11) + m02*m12)*x0*y0 - 
    2*(m012 + m112 + m122)*y02 - 4*(m00*m02*x0 + m02*m22*x0 + m12*(m11 + m22)*y0 + m01*(m12*x0 + m02*y0))*z0 - 
   2*(m022 + m122 + m222)*z02;
  coeffpa[6] = b;
#else
   b = -1 + m00*x0*x0 + 2*m01*x0*y0 + m11*y0*y0 + 
    2*m02*x0*z0 + 2*m12*y0*z0 + m22*z0*z0;

#if 0
  printf("M=%f %f %f %f %f %f\n", m00, m01, m02, m11, m12, m22);
  printf("x0=%f %f %f\n", x0, y0, z0);
  printf("boh=%f\n", x0*(m00*x0 + m01*y0 + m02*z0) + y0*(m01*x0 + m11*y0 + m12*z0) + 
    z0*(m02*x0 + m12*y0 + m22*z0));

  printf("m00*x0=%f\n", Mjpp[0][0]*x0);
  printf("m01*y0=%f\n", m01*y0);
  printf("m02*z0=%f\n", m02*z0);
#endif
 coeffpa[0] = -(PowerI(PowerI(m02,2)*m11 - 2*m01*m02*m12 + PowerI(m01,2)*m22 + 
       m00*(PowerI(m12,2) - m11*m22),2)*
     (-b + m00*PowerI(x0,2) + 2*m01*x0*y0 + m11*PowerI(y0,2) + 2*m02*x0*z0 + 
       2*m12*y0*z0 + m22*PowerI(z0,2)));
  coeffpa[1] = -2*(PowerI(m01,2) + PowerI(m02,2) - m00*m11 + PowerI(m12,2) - 
     (m00 + m11)*m22)*(PowerI(m02,2)*m11 - 2*m01*m02*m12 + 
     PowerI(m01,2)*m22 + m00*(PowerI(m12,2) - m11*m22))*
   (-b + m00*PowerI(x0,2) + 2*m01*x0*y0 + m11*PowerI(y0,2) + 2*m02*x0*z0 + 
     2*m12*y0*z0 + m22*PowerI(z0,2));
  coeffpa[2] = b*(PowerI(m01,4) + PowerI(m02,4) + PowerI(m00,2)*PowerI(m11,2) - 
      2*PowerI(m00,2)*PowerI(m12,2) - 4*m00*m11*PowerI(m12,2) + 
      PowerI(m12,4) + 4*m00*m11*(m00 + m11)*m22 - 
      2*(2*m00 + m11)*PowerI(m12,2)*m22 + 
      (PowerI(m00,2) + 4*m00*m11 + PowerI(m11,2))*PowerI(m22,2) + 
      4*m01*m02*m12*(m00 + m11 + m22) - 
      2*PowerI(m02,2)*(2*m00*m11 + PowerI(m11,2) - PowerI(m12,2) + m00*m22 + 
         2*m11*m22) + 2*PowerI(m01,2)*
       (PowerI(m02,2) - m00*m11 + PowerI(m12,2) - 2*(m00 + m11)*m22 - 
         PowerI(m22,2))) + PowerI(m02,2)*m11*PowerI(m12,2)*PowerI(x0,2) - 
   2*m01*m02*PowerI(m12,3)*PowerI(x0,2) - 
   PowerI(m02,2)*PowerI(m11,2)*m22*PowerI(x0,2) + 
   2*m01*m02*m11*m12*m22*PowerI(x0,2) + 
   PowerI(m01,2)*PowerI(m12,2)*m22*PowerI(x0,2) - 
   PowerI(m01,2)*m11*PowerI(m22,2)*PowerI(x0,2) - 
   PowerI(m00,3)*(PowerI(m11,2) - 2*PowerI(m12,2) + 4*m11*m22 + PowerI(m22,2))*
    PowerI(x0,2) - 2*PowerI(m01,5)*x0*y0 - 
   4*PowerI(m01,3)*PowerI(m02,2)*x0*y0 - 2*m01*PowerI(m02,4)*x0*y0 + 
   4*m01*PowerI(m02,2)*PowerI(m11,2)*x0*y0 - 
   8*PowerI(m01,2)*m02*m11*m12*x0*y0 - 2*PowerI(m02,3)*m11*m12*x0*y0 - 
   4*PowerI(m01,3)*PowerI(m12,2)*x0*y0 - 2*m01*PowerI(m12,4)*x0*y0 + 
   8*PowerI(m01,3)*m11*m22*x0*y0 + 10*m01*PowerI(m02,2)*m11*m22*x0*y0 - 
   14*PowerI(m01,2)*m02*m12*m22*x0*y0 + 4*m01*m11*PowerI(m12,2)*m22*x0*y0 + 
   6*PowerI(m01,3)*PowerI(m22,2)*x0*y0 - 
   2*m01*PowerI(m11,2)*PowerI(m22,2)*x0*y0 - PowerI(m01,4)*m11*PowerI(y0,2) - 
   2*PowerI(m01,2)*PowerI(m02,2)*m11*PowerI(y0,2) + 
   2*PowerI(m02,2)*PowerI(m11,3)*PowerI(y0,2) - 
   2*m01*PowerI(m02,3)*m12*PowerI(y0,2) - 
   4*m01*m02*PowerI(m11,2)*m12*PowerI(y0,2) - 
   2*PowerI(m01,2)*m11*PowerI(m12,2)*PowerI(y0,2) - 
   2*PowerI(m02,2)*m11*PowerI(m12,2)*PowerI(y0,2) - 
   m11*PowerI(m12,4)*PowerI(y0,2) + 
   PowerI(m01,2)*PowerI(m02,2)*m22*PowerI(y0,2) + 
   4*PowerI(m01,2)*PowerI(m11,2)*m22*PowerI(y0,2) + 
   4*PowerI(m02,2)*PowerI(m11,2)*m22*PowerI(y0,2) - 
   4*m01*m02*m11*m12*m22*PowerI(y0,2) + 
   2*PowerI(m11,2)*PowerI(m12,2)*m22*PowerI(y0,2) + 
   2*PowerI(m01,2)*m11*PowerI(m22,2)*PowerI(y0,2) - 
   PowerI(m11,3)*PowerI(m22,2)*PowerI(y0,2) - 
   2*(m02*(PowerI(m02,4) + PowerI(m02,2)*
          (-3*PowerI(m11,2) + 2*PowerI(m12,2) - 4*m11*m22) + 
         PowerI(PowerI(m12,2) - m11*m22,2))*x0 + 
      m12*(PowerI(m02,4) + PowerI(PowerI(m12,2) - m11*m22,2) - 
         2*PowerI(m02,2)*(PowerI(m11,2) - PowerI(m12,2) + 2*m11*m22))*y0 + 
      PowerI(m01,3)*m22*(m12*x0 + m02*y0) + 
      PowerI(m01,4)*(m02*x0 + m12*y0) + 
      m01*m02*(m02*m12*(7*m11 + 4*m22)*x0 + PowerI(m02,2)*m11*y0 + 
         4*PowerI(m12,2)*(m11 + m22)*y0) + 
      PowerI(m01,2)*(2*PowerI(m02,3)*x0 - m02*m22*(5*m11 + 2*m22)*x0 + 
         2*m12*(PowerI(m12,2) - m22*(2*m11 + m22))*y0))*z0 - 
   (2*PowerI(m01,3)*m02*m12 + 4*m01*m02*m12*m22*(m11 + m22) + 
      m22*(PowerI(m02,4) + PowerI(PowerI(m12,2) - m11*m22,2) - 
         2*PowerI(m02,2)*(PowerI(m11,2) - PowerI(m12,2) + 2*m11*m22)) - 
      PowerI(m01,2)*(PowerI(m02,2)*(m11 - 2*m22) + 
         2*m22*(-PowerI(m12,2) + 2*m11*m22 + PowerI(m22,2))))*PowerI(z0,2) + 
   PowerI(m00,2)*(4*m11*PowerI(m12,2)*PowerI(x0,2) - 
      4*PowerI(m11,2)*m22*PowerI(x0,2) + 4*PowerI(m12,2)*m22*PowerI(x0,2) - 
      4*m11*PowerI(m22,2)*PowerI(x0,2) + 
      2*PowerI(m02,2)*(2*m11 + m22)*PowerI(x0,2) + 
      2*PowerI(m01,2)*(m11 + 2*m22)*PowerI(x0,2) - 
      PowerI(m11,3)*PowerI(y0,2) + 2*m11*PowerI(m12,2)*PowerI(y0,2) - 
      4*PowerI(m11,2)*m22*PowerI(y0,2) - PowerI(m12,2)*m22*PowerI(y0,2) - 
      2*m01*x0*(2*m02*m12*x0 + 
         (PowerI(m11,2) - 2*PowerI(m12,2) + 4*m11*m22 + PowerI(m22,2))*y0) - 
      2*m02*(PowerI(m11,2) - 2*PowerI(m12,2) + 4*m11*m22 + PowerI(m22,2))*x0*
       z0 - 2*m12*(PowerI(m11,2) - 3*PowerI(m12,2) + 5*m11*m22 + 
         PowerI(m22,2))*y0*z0 - 
      (m11*PowerI(m12,2) - 2*PowerI(m12,2)*m22 + 4*m11*PowerI(m22,2) + 
         PowerI(m22,3))*PowerI(z0,2)) + 
   m00*(-(PowerI(m01,4)*PowerI(x0,2)) - PowerI(m02,4)*PowerI(x0,2) + 
      4*PowerI(m01,3)*(m11 + 2*m22)*x0*y0 + 
      4*PowerI(m02,3)*(2*m11 + m22)*x0*z0 - 
      2*m02*(PowerI(m12,2) - m11*m22)*x0*(m12*y0 - 5*m11*z0 - 4*m22*z0) - 
      4*(m11 + m22)*(-PowerI(m12,2) + m11*m22)*
       (m11*PowerI(y0,2) + z0*(2*m12*y0 + m22*z0)) + 
      PowerI(m01,2)*(-2*PowerI(m02,2)*PowerI(x0,2) + 4*m11*m22*PowerI(x0,2) + 
         2*PowerI(m22,2)*PowerI(x0,2) + 2*PowerI(m11,2)*PowerI(y0,2) + 
         4*m11*m22*PowerI(y0,2) - PowerI(m22,2)*PowerI(y0,2) + 
         2*m12*(2*m11 + 5*m22)*y0*z0 + 4*PowerI(m22,2)*PowerI(z0,2) + 
         4*m02*x0*(-2*m12*y0 + (m11 + 2*m22)*z0) + 
         PowerI(m12,2)*(-2*PowerI(x0,2) + PowerI(z0,2))) + 
      PowerI(m02,2)*(PowerI(m12,2)*(-2*PowerI(x0,2) + PowerI(y0,2)) + 
         10*m11*m12*y0*z0 + 4*m12*m22*y0*z0 + 2*PowerI(m22,2)*PowerI(z0,2) + 
         PowerI(m11,2)*(2*PowerI(x0,2) + 4*PowerI(y0,2) - PowerI(z0,2)) + 
         4*m11*m22*(PowerI(x0,2) + PowerI(z0,2))) + 
      2*m01*(2*PowerI(m02,2)*x0*(2*m11*y0 + m22*y0 - 2*m12*z0) - 
         (PowerI(m12,2) - m11*m22)*x0*(-4*m11*y0 - 5*m22*y0 + m12*z0) + 
         m02*(-7*PowerI(m12,2)*y0*z0 + m11*m22*y0*z0 + 
            m12*m22*(-2*PowerI(x0,2) + PowerI(y0,2) - 2*PowerI(z0,2)) + 
            m11*m12*(-2*(PowerI(x0,2) + PowerI(y0,2)) + PowerI(z0,2)))));
  coeffpa[3] = 2*(b*(2*m01*m02*m12 - m11*(2*PowerI(m02,2) + PowerI(m12,2)) - 
        (PowerI(m02,2) - PowerI(m11,2) + PowerI(m12,2))*m22 + 
        m11*PowerI(m22,2) + PowerI(m00,2)*(m11 + m22) - 
        PowerI(m01,2)*(m11 + 2*m22) + 
        m00*(-PowerI(m01,2) - PowerI(m02,2) + PowerI(m11,2) - 
           2*PowerI(m12,2) + 4*m11*m22 + PowerI(m22,2))) - 
     PowerI(m02,2)*PowerI(m11,2)*PowerI(x0,2) + 
     2*m01*m02*m11*m12*PowerI(x0,2) - PowerI(m01,2)*m11*m22*PowerI(x0,2) - 
     PowerI(m02,2)*m11*m22*PowerI(x0,2) + 2*m01*m02*m12*m22*PowerI(x0,2) - 
     PowerI(m01,2)*PowerI(m22,2)*PowerI(x0,2) - 
     PowerI(m00,3)*(m11 + m22)*PowerI(x0,2) + 2*PowerI(m01,3)*m11*x0*y0 + 
     6*m01*PowerI(m02,2)*m11*x0*y0 - 8*PowerI(m01,2)*m02*m12*x0*y0 + 
     2*m01*m11*PowerI(m12,2)*x0*y0 + 6*PowerI(m01,3)*m22*x0*y0 + 
     2*m01*PowerI(m02,2)*m22*x0*y0 - 2*m01*PowerI(m11,2)*m22*x0*y0 + 
     2*m01*PowerI(m12,2)*m22*x0*y0 - 2*m01*m11*PowerI(m22,2)*x0*y0 + 
     PowerI(m01,2)*PowerI(m11,2)*PowerI(y0,2) + 
     2*PowerI(m02,2)*PowerI(m11,2)*PowerI(y0,2) - 
     2*m01*m02*m11*m12*PowerI(y0,2) + 
     PowerI(m11,2)*PowerI(m12,2)*PowerI(y0,2) + 
     2*PowerI(m01,2)*m11*m22*PowerI(y0,2) - PowerI(m11,3)*m22*PowerI(y0,2) + 
     2*m01*m02*m12*m22*PowerI(y0,2) + m11*PowerI(m12,2)*m22*PowerI(y0,2) - 
     PowerI(m01,2)*PowerI(m22,2)*PowerI(y0,2) - 
     PowerI(m11,2)*PowerI(m22,2)*PowerI(y0,2) + 
     2*(-4*m01*m02*m12 + PowerI(m02,2)*(3*m11 + m22) + 
        PowerI(m01,2)*(m11 + 3*m22) - (m11 + m22)*(-PowerI(m12,2) + m11*m22))
       *(m02*x0 + m12*y0)*z0 + 
     (2*m01*m02*m12*(m11 - m22) + 
        PowerI(m02,2)*(-PowerI(m11,2) + 2*m11*m22 + PowerI(m22,2)) + 
        m22*(2*PowerI(m01,2)*m22 - (m11 + m22)*(-PowerI(m12,2) + m11*m22)))*
      PowerI(z0,2) + PowerI(m00,2)*
      (PowerI(m01,2)*PowerI(x0,2) + PowerI(m02,2)*PowerI(x0,2) - 
        PowerI(m11,2)*PowerI(x0,2) + 2*PowerI(m12,2)*PowerI(x0,2) - 
        4*m11*m22*PowerI(x0,2) - PowerI(m22,2)*PowerI(x0,2) - 
        2*m01*(m11 + m22)*x0*y0 - PowerI(m11,2)*PowerI(y0,2) - 
        PowerI(m12,2)*PowerI(y0,2) - 2*(m11 + m22)*(m02*x0 + m12*y0)*z0 - 
        (PowerI(m12,2) + PowerI(m22,2))*PowerI(z0,2)) + 
     m00*(2*PowerI(m02,2)*m11*PowerI(x0,2) + PowerI(m02,2)*m22*PowerI(x0,2) + 
        2*PowerI(m01,3)*x0*y0 - PowerI(m11,3)*PowerI(y0,2) + 
        2*m11*PowerI(m12,2)*PowerI(y0,2) - 4*PowerI(m11,2)*m22*PowerI(y0,2) - 
        PowerI(m12,2)*m22*PowerI(y0,2) + 
        2*(PowerI(m02,2) - PowerI(m11,2) + 3*PowerI(m12,2) - 5*m11*m22 - 
           PowerI(m22,2))*(m02*x0 + m12*y0)*z0 - 
        (m11*PowerI(m12,2) + PowerI(m02,2)*(m11 - m22) - 
           2*PowerI(m12,2)*m22 + 4*m11*PowerI(m22,2) + PowerI(m22,3))*
         PowerI(z0,2) + PowerI(m01,2)*
         (2*m22*PowerI(x0,2) - m22*PowerI(y0,2) + 
           m11*(PowerI(x0,2) + PowerI(y0,2)) + 2*m02*x0*z0 + 2*m12*y0*z0) + 
        2*m01*(PowerI(m02,2)*x0*y0 - 
           (PowerI(m11,2) - 3*PowerI(m12,2) + 5*m11*m22 + PowerI(m22,2))*x0*
            y0 + m02*m12*(-PowerI(x0,2) + PowerI(y0,2) + PowerI(z0,2)))));
  coeffpa[4] = b*(PowerI(m00,2) + PowerI(m11,2) - 
      2*(PowerI(m01,2) + PowerI(m02,2) + PowerI(m12,2)) + 4*m11*m22 + 
      PowerI(m22,2) + 4*m00*(m11 + m22)) - PowerI(m00,3)*PowerI(x0,2) + 
   2*m00*PowerI(m02,2)*PowerI(x0,2) - 4*PowerI(m00,2)*m11*PowerI(x0,2) - 
   PowerI(m01,2)*m11*PowerI(x0,2) - 4*PowerI(m02,2)*m11*PowerI(x0,2) + 
   6*m01*m02*m12*PowerI(x0,2) - 4*PowerI(m00,2)*m22*PowerI(x0,2) - 
   4*PowerI(m01,2)*m22*PowerI(x0,2) - PowerI(m02,2)*m22*PowerI(x0,2) - 
   2*PowerI(m00,2)*m01*x0*y0 + 6*PowerI(m01,3)*x0*y0 + 
   6*m01*PowerI(m02,2)*x0*y0 - 10*m00*m01*m11*x0*y0 - 
   2*m01*PowerI(m11,2)*x0*y0 - 2*m00*m02*m12*x0*y0 - 2*m02*m11*m12*x0*y0 + 
   6*m01*PowerI(m12,2)*x0*y0 - 8*m00*m01*m22*x0*y0 - 8*m01*m11*m22*x0*y0 - 
   2*m02*m12*m22*x0*y0 + 2*PowerI(m01,2)*m11*PowerI(y0,2) - 
   4*m00*PowerI(m11,2)*PowerI(y0,2) - PowerI(m11,3)*PowerI(y0,2) + 
   6*m01*m02*m12*PowerI(y0,2) - 4*m00*PowerI(m12,2)*PowerI(y0,2) + 
   2*m11*PowerI(m12,2)*PowerI(y0,2) - 4*PowerI(m01,2)*m22*PowerI(y0,2) - 
   4*PowerI(m11,2)*m22*PowerI(y0,2) - PowerI(m12,2)*m22*PowerI(y0,2) + 
   m00*PowerI(m01,2)*(2*PowerI(x0,2) - PowerI(y0,2)) - 
   2*((PowerI(m00,2)*m02 - 3*PowerI(m01,2)*m02 + m01*m12*(m11 + m22) + 
         m00*(4*m02*m11 + m01*m12 + 5*m02*m22) + 
         m02*(-3*(PowerI(m02,2) + PowerI(m12,2)) + 4*m11*m22 + PowerI(m22,2)))
        *x0 + (m00*m01*m02 - 3*PowerI(m01,2)*m12 + m01*m02*(m11 + m22) + 
         4*m00*m12*(m11 + m22) + 
         m12*(-3*PowerI(m02,2) + PowerI(m11,2) - 3*PowerI(m12,2) + 
            5*m11*m22 + PowerI(m22,2)))*y0)*z0 - 
   (m00*PowerI(m02,2) + 4*PowerI(m02,2)*m11 - 6*m01*m02*m12 + 
      4*m00*PowerI(m12,2) + m11*PowerI(m12,2) - 
      2*(PowerI(m02,2) + PowerI(m12,2))*m22 + 4*(m00 + m11)*PowerI(m22,2) + 
      PowerI(m22,3))*PowerI(z0,2);
  coeffpa[5] = 2*b*(m00 + m11 + m22) - 2*(PowerI(m00,2) + PowerI(m01,2) + PowerI(m02,2))*
    PowerI(x0,2) - 4*(m01*(m00 + m11) + m02*m12)*x0*y0 - 
   2*(PowerI(m01,2) + PowerI(m11,2) + PowerI(m12,2))*PowerI(y0,2) - 
   4*(m00*m02*x0 + m02*m22*x0 + m12*(m11 + m22)*y0 + 
      m01*(m12*x0 + m02*y0))*z0 - 
   2*(PowerI(m02,2) + PowerI(m12,2) + PowerI(m22,2))*PowerI(z0,2);
  coeffpa[6] = b;
#endif
#if 1
  /* ci deve essere sempre un solo zero positivo e quindi calcolo solo quello
   * usando zbrent con un opportuno bracketing (è la soluzione più veloce) */
#ifdef MC_ALPHAGUESS
  coeffq[0] = -1 + m00*x02+ 2*m01*x0*y0 + m11*y02 + 2*m02*x0*z0 + 2*m12*y0*z0 + m22*z02;
  coeffq[1] = -2*(m00*x02 + 2*m01*x0*y0 + m11*y02 + 2*m02*x0*z0 + 2*m12*y0*z0 + m22*z02);
  coeffq[2] = m00*x02 + 2*m01*x0*y0 + m11*y02 + 2*m02*x0*z0 + 2*m12*y0*z0 + m22*z02; 
  solve_quadratic(coeffq, &numsol, solq);
#if 0
  if (solq[0] < 0.0 || solq[1] < 0.0)
    printf("solq=%f %f\n", solq[0], solq[1]);
#endif
  if (fabs(solq[0]) < fabs(solq[1]))
    alq = solq[0];
  else
    alq = solq[1]; 
  xx[0] = alq*x0;
  xx[1] = alq*y0;
  xx[2] = alq*z0;
  for (kk1=0; kk1 < 3; kk1++)
    {
      gg[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
        gg[kk1] += Mjpp[kk1][kk2]*xx[kk2];
    } 
  //printf("x0=%f %f %f alg=%.15G norm rj0=%f %f\n",x0, y0, z0, alg, calc_norm(r0jpp), calc_norm(gg));
  ng = calc_norm(gg);
  if (ng > 0.0)
    alg=ng;
  else
    alg=0.5;
  
  //printf("alg=%f\n", alg);
  ig=1;
#else
  alg=0;
  ig=0;
#endif
  aR = (ig==1)?alg:1.0;
  aL= 0.0;
  dL = coeffpa[0];
  it=0;
  while (dL*polyalpha(coeffpa,aR) > 0.0)
    {
      aL = aR; 
      aR *= GOLD;
      it++;
    }
  //if (it > 0)
    //printf("it=%d aL=%f aR=%f\n", it, aL, aR);
  for (kk1=0; kk1 < 7; kk1++)
    coeffpa_HE[kk1] = coeffpa[kk1];
  /* N.B. con i criteri di terminazione di Bini o Cameron la tolleranza (5.0*DBL_EPSILON) 
   * passata come argomento non viene usata */
  alpha = rtsafe(coeffpa, alg, aL, aR, 5.0*DBL_EPSILON, ig, &ok); 

  if (!ok)
    {
      dist=distSq2origM(alpha, m00, m01, m02, m11, m12, m22, x0, y0, z0) - 1.0;
      printf("alpha=%.20G dist=%.16G dbleps=%.20G\n", alpha, dist, DBL_EPSILON);
      solve_numrec(coeffpa, roots, &ok);

      printf("\nri=%f %f %f rj=%f %f %f\n", rx[i], ry[i], rz[i], rx[j], ry[j], rz[j]);
      
      for (kk1=0; kk1 < 3; kk1++)
        for (kk2=0; kk2 < 3; kk2++)
          {
            printf("R[%d][%d]=%.15G\n", kk1, kk2, R[i][kk1][kk2]);
          }
      for (kk1=0; kk1 < 7; kk1++)
        {
          printf("coeff[%d]=%.15G\n", kk1, coeffpa[kk1]);
        }
      np=0;
      for (kk1=0; kk1 < 6; kk1++)
        {
          //printf("root[%d]=%.15G %.15G\n", kk1, creal(roots[kk1]), cimag(roots[kk1]));
          if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
            np++;
          if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
            {
              //dist=distSq2origM(creal(roots[kk1]), m00, m01, m02, m11, m12, m22, x0, y0, z0)-1.0;
              printf("dist=%.16G alpha=%.20G x=%f %f %f\n", dist, creal(roots[kk1]), x0, y0, z0);
              //exit(-1);
            }
        }
      printf("np=%d\n", np);
      store_bump(i, j);
      exit(-1);
    }
#if defined(HE_FAST_DIST)
  // il codice seguente non funziona...poiché l'idea è sbagliata
  intersectPoint(alpha, m00, m01, m02, m11, m12, m22, m012, m022, m122, m222, x0, y0, z0, sai, pi_HE, pj_HE);
  // ora calcolo il punto sull'altro ellissoide usando il gradiente 
#if 0
  gi[0] = 2.0*sai[0]*pi_HE[0];
  gi[1] = 2.0*sai[1]*pi_HE[1];
  gi[2] = 2.0*sai[2]*pi_HE[2];
  for (kk1=0; kk1 < 3; kk1++)
    {
      Mg[kk1]=0;
      for (kk2=0; kk2 < 3; kk2++)
        {
          Mg[kk1] += Mjp[kk1][kk2]*gi[kk2]; 
        }
    }
   gMg=scalProd(gi, Mg);// coefficiente del termine quadratico
   
   for (kk1=0; kk1 < 3; kk1++)
     pmx0[kk1] = pi_HE[kk1] - r0jp[kk1]; 
   for (kk1=0; kk1 < 3; kk1++)
    {
      Mp[kk1]=0;
      for (kk2=0; kk2 < 3; kk2++)
        {
          Mp[kk1] += Mjp[kk1][kk2]*pmx0[kk2]; 
        }
    }
   pMp=scalProd(pmx0, Mp)-1.0;// termine di ordine 0
   gMp=2.0*scalProd(gi, Mp);// termine lineare 

   cq_fd[0] = pMp;
   cq_fd[1] = gMp;
   cq_fd[2] = gMg;
   printf("coeff=%f %f %f\n", cq_fd[0], cq_fd[1], cq_fd[2]);
   solve_quadratic(cq_fd, &numsol_fd, solq_fd);
   if (fabs(solq_fd[0]) < fabs(solq_fd[1]))
     lam_fd = solq_fd[0];
   else
     lam_fd = solq_fd[1]; 
   printf("numsol=%d lam=%f %f\n", numsol_fd, solq_fd[0], solq_fd[1]);
   for (kk1=0; kk1 < 3; kk1++)
     {
       pj_HE[kk1] = pi_HE[kk1] + lam_fd*gi[kk1];
       dpij[kk1] = pj_HE[kk1] - pi_HE[kk1];
     }
   dist=calc_norm(dpij);
   
   //printf("x0=%f y0=%f z0=%f dist=%.15G sai=%f %f %f\n", x0, y0, z0, dist, sai[0], sai[1], sai[2]);
   return (lam_fd < 0.0)?-dist:dist;
#endif
#if 0
    { 
      double pi[3], pj[3], gi[3], gj[3], pij[3];
  printf("fi(pi)=%.15G\n", Sqr(pi_HE[0]/sai[0])+ Sqr(pi_HE[1]/sai[1])+ Sqr(pi_HE[2]/sai[2])-1.0);
  printf("fj(pj)=%.15G\n", calcfel(Mjp,r0jp, pj_HE));
  gi[0] = pi_HE[0]/Sqr(sai[0]);
  gi[1] = pi_HE[1]/Sqr(sai[1]);
  gi[2] = pi_HE[2]/Sqr(sai[2]);
#if 1
  for (kk1=0; kk1 < 3; kk1++)
    {
      gj[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
        {
          gj[kk1] += Mjp[kk1][kk2]*(pj_HE[kk2]-r0jp[kk2]);
        }
    }
#endif
  for (kk1=0; kk1 < 3; kk1++)
    pij[kk1] = pi_HE[kk1]-pj_HE[kk1];

  printf("NEW gi=%.15G %.15G %.15G gj=%.15G %.15G %.15G\n", gi[0], gi[1], gi[2], gj[0], gj[1], gj[2]);
  printf("NEW gi/gj=%.15G %.15G %.15G\n", gi[0]/gj[0], gi[1]/gj[1], gi[2]/gj[2]);
  printf("pij/gj =   %.15G %.15G %.15G\n", pij[0]/gi[0], pij[1]/gi[1],pij[2]/gi[2]);
    }
#endif
  //printf("x0=%f y0=%f z0=%f dist=%.15G sai=%f %f %f\n", x0, y0, z0, dist, sai[0], sai[1], sai[2]);
  return dist;
#else
#ifdef HE_OPT
  dist=distSq2origMopt(alpha, m00, m01, m02, m11, m12, m22, m012, m022, m122, m222, x0, y0, z0) - 1.0;
#ifdef DIST_OF_CLOSE_APPR
  intersectPoint(alpha, m00, m01, m02, m11, m12, m22, m012, m022, m122, m222, x0, y0, z0, sai, pi_HE, pj_HE);
  for (kk1=0; kk1 < 3; kk1++)
    {
      docc[kk1]= r0jp[kk1] + (pi_HE[kk1]-pj_HE[kk1]);  
    }
  doccn=calc_norm(docc); 
  //printf("dist of closest cont=%.15G actual dist=%.15G\n", doccn, calc_norm(r0jp));
#endif
#else
  dist=distSq2origM(alpha, m00, m01, m02, m11, m12, m22, x0, y0, z0) - 1.0;
#endif
#endif
  /* la distanza di closest-approach si puo' calcolare trasformando i due punti che determinano la distanza tra sfera ed ellissoide */
  if (dist < 0.0)
    {
      if (calc_norm(r0jp) > doccn)
        printf("[neg] dist of closest cont=%.15G actual dist=%.15G\n", doccn, calc_norm(r0jp));
      return -1.0;
    }
  else
    {

      if (calc_norm(r0jp) < doccn)
        printf("[pos] dist of closest cont=%.15G actual dist=%.15G\n", doccn, calc_norm(r0jp));
      return 1.0;
    }
#else
  solve_numrec(coeffpa, roots, &ok);
  if (!ok)
    {
      for (kk1=0; kk1 < 7; kk1++)
        {
          printf("c[%d]=%.15G\n", kk1, coeffpa[kk1]);
        }
      printf("M=%f %f %f %f %f %f\n", m00, m01, m02, m11, m12, m22);
      printf("x0= %f %f %f\n", x0, y0, z0);
      printf("b=%.15G\n", b);
      exit(-1);
    }
  for (kk1=0; kk1 < 6; kk1++)
    {
      //printf("root[%d]=%.15G %.15G\n", kk1, creal(roots[kk1]), cimag(roots[kk1]));
      if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
      np++;
      if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
        {
          dist=distSq2origM(creal(roots[kk1]), m00, m01, m02, m11, m12, m22, x0, y0, z0)-1.0;
          //printf("dist=%f alpha=%.15G x=%f %f %f\n", dist, creal(roots[kk1]), x0, y0, z0);
#if 0
          double distpw;
          distpw=check_overlap_pw(i, j, shift);
          if (dist*distpw < 0)
            {
              printf("evalpoly(alpha)=%.16G\n",polyalpha(coeffpa,roots[kk1]));
              printf("dist=%f distpw=%f\n", dist, distpw);
              printf("alpha=%.15G\n", creal(roots[kk1]));
              printf("sa=%f %f %f x=%f %f %f dist=%f\n", sx, sy, sz, x0, y0, z0, distSq2orig(creal(roots[kk1]), sx, sy, sz, x0, y0, z0));
#if 0
              rx[i] = ry[i] = rz[i] = 0.0;
              rx[j] = r0jpp[0];
              ry[j] = r0jpp[1];
              rz[j] = r0jpp[2];
              typesArr[typej].sax[0] = sx;
              typesArr[typej].sax[1] = sy;
              typesArr[typej].sax[2] = sz;
              for (kk1=0; kk1 < 3; kk1++)
                {
                  for (kk2=0; kk2 < 3; kk2++)
                    {
                      R[i][kk1][kk2] = (kk1==kk2)?1.0:0.0;
                      R[j][kk1][kk2] = Rjpp[kk1][kk2];
                    }
                }

#endif
#if 1
              rx[i] = ry[i] = rz[i] = 0.0;
              rx[j] = r0jp[0];
              ry[j] = r0jp[1];
              rz[j] = r0jp[2];
              for (kk1=0; kk1 < 3; kk1++)
                {
                  for (kk2=0; kk2 < 3; kk2++)
                    {
                      R[i][kk1][kk2] = (kk1==kk2)?1.0:0.0;
                      R[j][kk1][kk2] = Rjp[kk1][kk2];
                    }
                }

#endif
#if 0
              typesArr[typej].sax[0] = sx;
              typesArr[typej].sax[1] = sy;
              typesArr[typej].sax[2] = sz;
              for (kk1=0; kk1 < 3; kk1++)
                {
                  for (kk2=0; kk2 < 3; kk2++)
                    {
                      R[i][kk1][kk2] = R[j][kk1][kk2] = (kk1==kk2)?1.0:0.0;
                    }
                }
              rx[i] = 0.0;
              ry[i] = 0.0;
              rz[i] = 0.0;
              rx[j] = x0;
              ry[j] = y0;
              rz[j] = z0;
#endif
              store_bump(i, j);
              exit(-1);
            }
#endif
          if (dist < 0.0)
            return -1.0;
          else
            return 1.0;

        }
    }
#endif
#if 0
      if (np > 1)
        {
          printf("boh\n");
          exit(-1);
        }
#endif 
#if 0
  wrap_dsyev_he(Mjpp, evec, eval, &ok);
  /* eigenvalues are provided by zeros of the 3th order characteristic polynomial */
#if 0
  coeff[0] = -Sqr(Mjpp[0][2])*Mjpp[1][1] + 2.0*Mjpp[0][1]*Mjpp[0][2]*Mjpp[1][2] - Mjpp[0][0]*Sqr(Mjpp[1][2]) - Sqr(Mjpp[0][1])*Mjpp[2][2] + Mjpp[0][0]*Mjpp[1][1]*Mjpp[2][2];
  coeff[1] = Sqr(Mjpp[0][1]) + Sqr(Mjpp[0][2]) - Mjpp[0][0]*Mjpp[1][1] + Sqr(Mjpp[1][2]) - Mjpp[0][0]*Mjpp[2][2] 
    - Mjpp[1][1]*Mjpp[2][2];
  coeff[2] = Mjpp[0][0] + Mjpp[1][1] + Mjpp[2][2];
  coeff[3] = -1.0;
  solve_cubic_analytic(coeff, eval);
#endif
  sx = 1.0/sqrt(eval[0]);
  sy = 1.0/sqrt(eval[1]);
  sz = 1.0/sqrt(eval[2]);
  //printf("semi-axes=%.15G %.15G %.15G\n x0=%f %f %f", sx, sy, sz, x0, y0, z0);

  //since we move to reference frame of ellipsoid j we need to rotate the position of ellipsoid i (which is now a sphere)
  for (kk1=0; kk1 < 3; kk1++)
    {
      r0jp3[kk1]=0.0;
      for (kk2=0; kk2 < 3; kk2++)
        {
          r0jp3[kk1] += evec[kk1][kk2]*r0jpp[kk2]; 
        }
    }
  x0 = r0jp3[0];
  y0 = r0jp3[1];
  z0 = r0jp3[2];
  if (sx==sz)
    {
      /* swap y and z axes */
      vt = sy;
      sy = sz;
      sz = vt;
      vt = y0;
      y0 = z0;
      z0 = vt;
    }
  else if (sy==sz)
    {
      at[0] = sx;
      at[1] = sy;
      at[2] = sz;
      sx = at[1];
      sy = at[2];
      sz = at[0];
      at[0] = x0;
      at[1] = y0;
      at[2] = z0;
      x0 = at[1];
      y0 = at[2];
      z0 = at[0];
    }
  sx2 = sx*sx;
  sy2 = sy*sy;
  sz2 = sz*sz;
  x02 = x0*x0;
  y02 = y0*y0;
  z02 = z0*z0; 
  sx4 = sx2*sx2;
  sy4 = sy2*sy2;
  sz4 = sz2*sz2;

  /* maybe check for equal semi-axes an reorder semi-axes accordingly */
  if (sx==sy)
    {
      /* uniaxial ellipsoids */
      coeffpa[0] = -1;
      coeffpa[1] = -2.0*sx2-2.0*sz2;
      coeffpa[2] = -sx4 - 4.0*sx2*sz2 - sz4 + sx2*x02 + sx2*y02 + sz2*z02;
      coeffpa[3]=-2.0*sx2*sz4 - sx4*sz*(sz - z0) - sx4*sz*(sz + z0) + 
        2.0*sx2*sz2*(x02 + y02 + z02);
      coeffpa[4] = sx2*sz4*(x02 + y02) - sx4*sz2*(sz - z0)*(sz + z0);
      oqs_quartic_solver(coeffpa, roots);
      /* un'unica soluzione deve essere positiva */
      for (kk1=0; kk1 < 6; kk1++)
        {
          if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
            {
              dist=distSq2orig(creal(roots[kk1]), sx, sy, sz, x0, y0, z0);
              if (dist < 1.0)
                return -1.0;
              else
                return 1.0;
            }
        }
    }
  else
    {
      coeffpa[0] = -1; 
      coeffpa[1] = -2.0*sx2 - 2.0*sy2 - 2.0*sz2;
      coeffpa[2] = -sx4 - 4.0*sx2*sy2 - sy4 - 4.0*sx2*sz2 - 4.0*sy2*sz2 - sz4 + sx2*x02 + sy2*y02 + sz2*z02;
      coeffpa[3] = -2.0*sx4*sy2 - 2.0*sx2*sy4 - 2.0*sx4*sz2 - 8.0*sx2*sy2*sz2 - 2.0*sy4*sz2 - 2.0*sx2*sz4 - 
        2.0*sy2*sz4 + 2.0*sx2*sy2*x02 + 2.0*sx2*sz2*x02 + 2.0*sx2*sy2*y02 + 2*sy2*sz2*y02 + 2.0*sx2*sz2*z02 + 2.0*sy2*sz2*z02;
      coeffpa[4] = -sx4*sy4-4.0*sx4*sy2*sz2 - 4.0*sx2*sy4*sz2 - sx4*sz4 - 4.0*sx2*sy2*sz4 - sy4*sz4 + 
        sx2*sy4*x02 + 4.0*sx2*sy2*sz2*x02 + sx2*sz4*x02 + sx4*sy2*y02 + 4.0*sx2*sy2*sz2*y02 + sy2*sz4*y02 + 
        sx4*sz2*z02 + 4.0*sx2*sy2*sz2*z02 + sy4*sz2*z02;
      coeffpa[5] =-2.0*sx4*sy4*sz2 - 2.0*sx4*sy2*sz4 - 2.0*sx2*sy4*sz4 + 2.0*sx2*sy4*sz2*x02 + 2.0*sx2*sy2*sz4*x02 
        + 2.0*sx4*sy2*sz2*y02 + 2.0*sx2*sy2*sz4*y02 + 2.0*sx4*sy2*sz2*z02 + 2.0*sx2*sy4*sz2*z02;
      coeffpa[6] = -sx4*sy4*sz4 + sx2*sy4*sz4*x02 + sx4*sy2*sz4*y02 + sx4*sy4*sz2*z02;
      /* NOTA: se l'origine degli assi è dentro all'ellissoide o il centro dell'ellissoide (x0,y0,z0)
       * è dentro la circonferenza unitaria i due ellissoidi si sovrappongono e non va calcolato nulla */
      solve_numrec(coeffpa, roots, &ok);
      /* un'unica soluzione deve essere positiva */
      //np=0;
      for (kk1=0; kk1 < 6; kk1++)
        {
          //printf("root[%d]=%.15G %.15G\n", kk1, creal(roots[kk1]), cimag(roots[kk1]));
          //if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
            //np++;
          if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
            {
              dist=distSq2orig(creal(roots[kk1]), sx, sy, sz, x0, y0, z0)-1.0;
              //printf("dist=%f alpha=%.15G sa=%f %f %f x=%f %f %f\n", dist, creal(roots[kk1]), sx, sy, sz, x0, y0, z0);
#if 0
              double distpw;
              distpw=check_overlap_pw(i, j, shift);
              if (dist*distpw < 0)
                {
                  printf("evalpoly(alpha)=%.16G\n",polyalpha(coeffpa,roots[kk1]));
                  printf("dist=%f distpw=%f\n", dist, distpw);
                  printf("alpha=%.15G\n", creal(roots[kk1]));
                  printf("sa=%f %f %f x=%f %f %f dist=%f\n", sx, sy, sz, x0, y0, z0, distSq2orig(creal(roots[kk1]), sx, sy, sz, x0, y0, z0));
#if 0
                  rx[i] = ry[i] = rz[i] = 0.0;
                  rx[j] = r0jpp[0];
                  ry[j] = r0jpp[1];
                  rz[j] = r0jpp[2];
                  typesArr[typej].sax[0] = sx;
                  typesArr[typej].sax[1] = sy;
                  typesArr[typej].sax[2] = sz;
                  for (kk1=0; kk1 < 3; kk1++)
                    {
                      for (kk2=0; kk2 < 3; kk2++)
                        {
                          R[i][kk1][kk2] = (kk1==kk2)?1.0:0.0;
                          R[j][kk1][kk2] = Rjpp[kk1][kk2];
                        }
                    }
                  
#endif
#if 1
                  rx[i] = ry[i] = rz[i] = 0.0;
                  rx[j] = r0jp[0];
                  ry[j] = r0jp[1];
                  rz[j] = r0jp[2];
                  for (kk1=0; kk1 < 3; kk1++)
                    {
                      for (kk2=0; kk2 < 3; kk2++)
                        {
                          R[i][kk1][kk2] = (kk1==kk2)?1.0:0.0;
                          R[j][kk1][kk2] = Rjp[kk1][kk2];
                        }
                    }
                  
#endif
#if 0
                  typesArr[typej].sax[0] = sx;
                  typesArr[typej].sax[1] = sy;
                  typesArr[typej].sax[2] = sz;
                  for (kk1=0; kk1 < 3; kk1++)
                    {
                      for (kk2=0; kk2 < 3; kk2++)
                        {
                          R[i][kk1][kk2] = R[j][kk1][kk2] = (kk1==kk2)?1.0:0.0;
                        }
                    }
                  rx[i] = 0.0;
                  ry[i] = 0.0;
                  rz[i] = 0.0;
                  rx[j] = x0;
                  ry[j] = y0;
                  rz[j] = z0;
#endif
                  store_bump(i, j);
                  exit(-1);
                }
#endif
              if (dist < 0.0)
                return -1.0;
              else
                return 1.0;
              
            }
        }
#if 0
      if (np > 1)
        {
          printf("boh\n");
          exit(-1);
        }
#endif 
    }
#endif
}
void tRDiagRqe2d(int i, double M[2][2], double D[2], double Ri[2][2])
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
#if 1
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

#else
double check_overlap_polyell_2D(int i, int j, double shift[3])
{
  const double GOLD=1.618034;
  double Rjp[2][2], r0jp[2], ri[2], rj[2], Di[2], Dj[2], Mjp[2][2], gradj[2], xppg[2], M2I[2][2], invM2I[2][2], oax[2], theta, norm;
  double RM[2][2], Mtmp[2][2], Mjpp[2][2], r0jpp[2], Mjp3[2][2], r0jp3[2], xp3g[2];
  double Mip[2][2], Mipp[2][2];
  double evec[2][2], eval[2], coeffpa[5], sx, sy, x0, y0, docc[2], doccn;
  double sx2, sy2, sx4, sy4, sz4, x02, y02, sai[2];
  double vt, at[2], Mi[2][2], Mj[2][2], Ri[2][2], Rj[2][2], dist, m00, m01, m11, b2d;
  double alpha;
  int np, typei, typej;
  int kk1, kk2, kk3, k1, k2, k3, ok;
  complex double roots[4];
 
  typei = typeOfPart[i];
  typej = typeOfPart[j];
  
  /* apply affinity to reduce first ellipsoid to a sphere */
  for (k1=0; k1 < 2; k1++)
    {
      Di[k1]= 1.0/Sqr(typesArr[typei].sax[k1]);
      Dj[k1]= 1.0/Sqr(typesArr[typej].sax[k1]);
      for (k2=0; k2 < 2; k2++)
        {
          Ri[k1][k2] = R[i][k1][k2];
          Rj[k1][k2] = R[j][k1][k2];
        }
    }
  /* sai[0] e sai[1] sono i semiassi della particella i */
  for (k1=0; k1 < 2; k1++)
    {
      sai[k1] = typesArr[typei].sax[k1];
    } 
  tRDiagRqe2d(i, Mi, Di, Ri);
  tRDiagRqe2d(j, Mj, Dj, Rj);
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
  tRDiagRqe2d(j, Mjp, Dj, Rjp);

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
  m00 = Mjpp[0][0];
  m01 = Mjpp[0][1];
  m11 = Mjpp[1][1];
  b2d = -1.0 + m00*x0*x0 + 2.0*m01*x0*y0 + m11*y0*y0;

  coeffpa[0] = -Sqr(m01*m01 - m00*m11)*(-b2d + m00*x0*x0 + y0*(2*m01*x0 + m11*y0)); 
  coeffpa[1] = -2*(m00 + m11)*(-m01*m01 + m00*m11)*(-b2d + m00*x0*x0 + 
   y0*(2*m01*x0 + m11*y0));
  coeffpa[2] =b2d*(m00*m00 - 2*m01*m01 + 4*m00*m11 + m11*m11) - (m00*m00*m00 - 2*m00*m01*m01 + 
    4*m00*m00*m11 + m01*m01*m11)*x0*x0 - 
    2*m01*(m00*m00 - 3*m01*m01 + 5*m00*m11 + m11*m11)*x0*y0 - (-2*m01*m01*m11 + 
    m11*m11*m11 + m00*(m01*m01 + 4*m11*m11))*y0*y0;
  coeffpa[3] = 2*b2d*(m00 + m11) - 2*(m00*m00 + m01*m01)*x0*x0 - 
    4*m01*(m00 + m11)*x0*y0 - 2*(m01*m01 + m11*m11)*y0*y0;
  coeffpa[4] = b2d;
  oqs_quartic_solver(coeffpa, roots);
  for (kk1=0; kk1 < 4; kk1++)
    {
      if (cimag(roots[kk1])==0 && creal(roots[kk1]) > 0.0)
        {
          alpha=creal(roots[kk1]);
          break;
        }
    }
  dist=distSq2origM2d(alpha, m00, m01, m11, x0, y0) - 1.0;
  /* trasformando tramite l'affinità inversa i punti che individuano la distanza tra sfera ed ellissoide
   * si avrà la distanza tra i due ellissoidi che si può usare nella dinamica event-driven */
  if (dist < 0.0)
    return -1.0;
  else
    return 1.0;
}
#endif
#endif
