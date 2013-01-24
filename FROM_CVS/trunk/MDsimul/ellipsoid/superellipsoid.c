#if defined(MD_SUPERELLIPSOID) 
#include<mdsimul.h>
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG10(x)  
#define MD_DEBUG11(x) 
#define MD_DEBUG15(x) 
#define MD_DEBUG20(x) 
#define MD_DEBUG29(x) 
#define MD_DEBUG30(x)  //qui 
#define MD_DEBUG31(x) //qui 
#define MD_DEBUG32(x)    
#define MD_DEBUG33(x) 
#define MD_DEBUG34(x) 
#define MD_DEBUG36(x) 
#define MD_DEBUG37(x) 
#define MD_DEBUG38(x) 
#define MD_NEGPAIRS
#define MD_NO_STRICT_CHECK
#define MD_OPTDDIST
#undef MD_FDJAC_SYM
#ifdef MD_ASYM_ITENS
extern double *phi0, *psi0, *costheta0, *sintheta0, **REt, **RE0, *angM, ***RM, **REtA, **REtB, **Rdot;
extern double cosEulAng[2][3], sinEulAng[2][3];
#endif
extern double rxC, ryC, rzC, trefG;
extern double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
extern void symtop_evolve_orient(int i, double ti, double **Ro, double **REt, double cosea[3], double sinea[3], double *phir, double *psir);
#ifdef EDHE_FLEX
extern void set_angmom_to_zero(int i);
extern int *is_a_sphere_NNL;
#endif
extern int icg, jcg, doneryck;
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **REtA, **REtB;
extern double cosEulAng[2][3], sinEulAng[2][3];
extern long long int itsFNL, timesFNL, timesSNL, itsSNL;
extern int do_check_negpairs;

extern int isSymItens(int i);
#ifdef EDHE_FLEX
extern int *mapbondsaFlex, *mapbondsbFlex, nbondsFlex;
extern double *mapBheightFlex, *mapBhinFlex, *mapBhoutFlex, *mapSigmaFlex; 
extern double *t2arr, *distsOld, *dists, *distsOld2, *maxddoti;
extern int *crossed, *tocheck, *dorefine, *crossed, *negpairs;
#endif
#ifdef MD_ASYM_ITENS
extern double **Ia, **Ib, **invIa, **invIb, **Iatmp, **Ibtmp, *angM;
#else
extern double Ia, Ib, invIa, invIb;
#endif
extern double gradplane[3];
extern struct LastBumpS *lastbump;
extern double *axa, *axb, *axc;
extern int *scdone;
extern double *maxax;
/* Routines for LU decomposition from Numerical Recipe online */
void ludcmpR(double **a, int* indx, double* d, int n);
void lubksbR(double **a, int* indx, double *b, int n);
extern void update_MSDrot(int i);
extern void InvMatrix(double **a, double **b, int NB);
extern void calc_angmom(int i, double **I);
double min(double a, double b);
extern void upd_refsysM(int i);
extern double calc_maxddot(int i, int j);
double zbrent(double (*func)(double), double x1, double x2, double tol);
extern double invaSq[2], invbSq[2], invcSq[2];
extern double rxC, ryC, rzC;
extern int SolveLineq (double **a, double *x, int n); 
extern double calc_norm(double *vec);
extern int calcdist_retcheck;
extern double rA[3], rB[3];
extern double gradplane_all[6][3], rBall[6][3];
extern int polinterr, polinterrRyck;
extern double max3(double a, double b, double c);
extern double scalProd(double *A, double *B);
extern int fdjac_disterr;
extern double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **powdirs;
extern void tRDiagR(int i, double **M, double a, double b, double c, double **Ri);
extern void lab2body(int i, double x[], double xp[], double *rO, double **R);
int is_superellipse(int i)
{
  //return 1;
  MD_DEBUG37(return 1);
  if (typesArr[typeOfPart[i]].n[0]==2.0 && typesArr[typeOfPart[i]].n[1]==2.0 &&
      typesArr[typeOfPart[i]].n[2]==2.0)
    return 0;
  else 
    return 1;
}
#undef MD_TWOPARAMSE
#ifdef MD_TWOPARAMSE
double calcf(double *x,int i)
{
  double a,b,c,e,n;
  a = typesArr[typeOfPart[i]].sax[0];
  b = typesArr[typeOfPart[i]].sax[1];
  c = typesArr[typeOfPart[i]].sax[2];
  e = typesArr[typeOfPart[i]].n[0];
  n = typesArr[typeOfPart[i]].n[1];
  return pow(pow(fabs(x[0])/a,2.0/e)+pow(fabs(x[1])/b,2.0/e),e/n)+pow(fabs(x[2])/c,2.0/n)-1.0; 
}
#define Sign(x) ((x>=0)?1:-1)
void calcfx(double *fx, double x, double y, double z, int i)
{
/* calcola il gradiente della superficie nel riferimento del corpo rigido*/
  double a, b, c, n, e;
  double xa2e, yb2e, A, inve2;
  a = typesArr[typeOfPart[i]].sax[0];
  b = typesArr[typeOfPart[i]].sax[1];
  c = typesArr[typeOfPart[i]].sax[2];
  e = typesArr[typeOfPart[i]].n[0];
  n = typesArr[typeOfPart[i]].n[1];
  //printf("a=%f b=%f c=%f e=%f n=%f\n", a,b,c,e,n);
#if 0
  inve2 = 2.0/e;	
  xa2e = pow(fabs(x)/a,inve2);
  yb2e = pow(fabs(y)/b,inve2);
  A = pow(xa2e+yb2e,-1.0+e/n);
  fx[0] = (2.0*xa2e*A)/(n*x);
  fx[1] = (2.0*yb2e*A)/(n*y);
  fx[2] = (2.0*pow(fabs(z)/c,2.0/n))/(n*z);
#else
 A =  pow(pow(fabs(x)/a,2.0/e)+pow(fabs(y)/b,2.0/e),e/n-1.0)/n;
 fx[0] = 2.0*Sign(x)*pow(fabs(x)/a,2.0/e-1.0)*A/a;
 fx[1] = 2.0*Sign(y)*pow(fabs(y)/b,2.0/e-1.0)*A/b;
 fx[2] = 2.0*Sign(z)*pow(fabs(z)/c,2.0/n-1.0)/(c*n);
#endif

}
#define Power(x,y) pow(x,y)
void calcfxx(double df[3][3], double x, double y, double z, int i)
{
  /* calcola fxx nel riferimento del corpo rigido */
  double a, b, c, n, e;
  double xa2e, yb2e, A, xa2eyb2e, emn, nm2, em2, inve2, eSqrnxy, eSqrnxSq;
  a = typesArr[typeOfPart[i]].sax[0];
  b = typesArr[typeOfPart[i]].sax[1];
  c = typesArr[typeOfPart[i]].sax[2];
  e = typesArr[typeOfPart[i]].n[0];
  n = typesArr[typeOfPart[i]].n[1];
#if 0
  inve2 = 2.0/e;
  xa2e = pow(fabs(x)/a,inve2);
  yb2e = pow(fabs(y)/b,inve2);
  eSqrnxy = e*Sqr(n)*x*y;
  eSqrnxSq = e*Sqr(n)*x*x;
  A = pow(xa2e+yb2e, -2.0+e/n);
  xa2eyb2e = xa2e*yb2e;
  emn = e-n;
  nm2 = n-2.0;
  em2 = e-2.0;
  df[0][0] = (A*(-2.0*e*nm2*Sqr(xa2e) - 
       2.0*em2*n*xa2eyb2e))/ (e*Sqr(n)*Sqr(x));
  df[0][1] = (4*emn*xa2eyb2e*A)/ (eSqrnxy); 
  df[0][2] = 0.0;
  df[1][0] = 4*emn*xa2eyb2e*A/(eSqrnxy);
  df[1][1] = (A*(-2.0*em2*n*xa2e*yb2e- 2.0*e*nm2*Sqr(yb2e)))/(e*Sqr(n)*Sqr(y));
  df[1][2] = df[2][0] = df[2][1] = 0.0;
  df[2][2] = -(2.0*nm2*pow(fabs(z)/c,2.0/n))/(Sqr(n)*Sqr(z)); 
#else
  df[0][0] = (4.0*(-1.0 + e/n)*Power(fabs(x)/a,-2.0 + 4.0/e)*
      Power(Power(fabs(x)/a,2.0/e) + Power(fabs(y)/b,2/e),-2.0 + e/n))/
    (a*a*e*n) + (2.0*(-1.0 + 2.0/e)*Power(fabs(x)/a,-2 + 2.0/e)*
      Power(Power(fabs(x)/a,2.0/e) + Power(fabs(y)/b,2.0/e),-1.0 + e/n))/(a*a*n);
  df[0][1] = Sign(x)*Sign(y)*(4*(-1 + e/n)*Power(fabs(x)/a,-1 + 2/e)*Power(fabs(y)/b,-1 + 2/e)*
     Power(Power(fabs(x)/a,2/e) + Power(fabs(y)/b,2/e),-2 + e/n))/(a*b*e*n);
  df[0][2] = 0.0;
  df[1][0] = Sign(y)*Sign(x)*(4*(-1 + e/n)*Power(fabs(x)/a,-1 + 2/e)*Power(fabs(y)/b,-1 + 2/e)*
     Power(Power(fabs(x)/a,2/e) + Power(fabs(y)/b,2/e),-2 + e/n))/(a*b*e*n);
  df[1][1] = (4*(-1 + e/n)*Power(fabs(y)/b,-2 + 4/e)*
      Power(Power(fabs(x)/a,2/e) + Power(fabs(y)/b,2/e),-2 + e/n))/
    (b*b*e*n) + (2*(-1 + 2/e)*Power(fabs(y)/b,-2 + 2/e)*
      Power(Power(fabs(x)/a,2/e) + Power(fabs(y)/b,2/e),-1 + e/n))/(b*b*n);
  df[1][2] = 0.0;
  df[2][0] = 0.0;
  df[2][1] = 0.0;
  df[2][2] =(2*Power(1/c,2/n)*(-1 + 2/n)*Power(fabs(z),-2 + 2/n))/n;
#endif
}

#else
double calcf(double *x, int i)
{
  double a,b,c,n1,n2,n3;
  if (OprogStatus.targetPhi > 0.0)
    {
      a = axa[i];
      b = axb[i];
      c = axc[i];
    }
  else
    {
      a = typesArr[typeOfPart[i]].sax[0];
      b = typesArr[typeOfPart[i]].sax[1];
      c = typesArr[typeOfPart[i]].sax[2];
    }
  n1 = typesArr[typeOfPart[i]].n[0];
  n2 = typesArr[typeOfPart[i]].n[1];
  n3 = typesArr[typeOfPart[i]].n[2];
  return pow(fabs(x[0])/a, n1)+pow(fabs(x[1])/b, n2)+pow(fabs(x[2])/c, n3)-1.0;
}
#define Sign(x) ((x>=0)?1:-1)
void calcfx(double *fx, double x, double y, double z, int i)
{
/* calcola il gradiente della superficie nel riferimento del corpo rigido*/
  double a, b, c, n1, n2, n3;
  double xa, yb, A, inve2;
  double an1, bn2, cn3;
  if (OprogStatus.targetPhi > 0.0)
    {
      a = axa[i];
      b = axb[i];
      c = axc[i];
    }
  else
    {
      a = typesArr[typeOfPart[i]].sax[0];
      b = typesArr[typeOfPart[i]].sax[1];
      c = typesArr[typeOfPart[i]].sax[2];
    }
  n1 = typesArr[typeOfPart[i]].n[0];
  n2 = typesArr[typeOfPart[i]].n[1];
  n3 = typesArr[typeOfPart[i]].n[2];
  an1=n1/pow(a,n1);
  bn2=n2/pow(b,n2);
  cn3=n3/pow(c,n3);
  //printf("a=%f b=%f c=%f e=%f n=%f\n", a,b,c,e,n);
  fx[0] = Sign(x)*an1*pow(fabs(x),n1-1.0);
  fx[1] = Sign(y)*bn2*pow(fabs(y),n2-1.0);
  fx[2] = Sign(z)*cn3*pow(fabs(z),n3-1.0);
}

void calcfxx(double df[3][3], double x, double y, double z, int i)
{
  /* calcola fxx nel riferimento del corpo rigido */
  double a, b, c, n1, n2, n3;
  double an1, bn2, cn3;
  if (OprogStatus.targetPhi > 0.0)
    {
      a = axa[i];
      b = axb[i];
      c = axc[i];
    }
  else
    {
      a = typesArr[typeOfPart[i]].sax[0];
      b = typesArr[typeOfPart[i]].sax[1];
      c = typesArr[typeOfPart[i]].sax[2];
    }
  n1 = typesArr[typeOfPart[i]].n[0];
  n2 = typesArr[typeOfPart[i]].n[1];
  n3 = typesArr[typeOfPart[i]].n[2];
  an1=(n1)*(n1-1)/pow(a,n1);
  bn2=(n2)*(n2-1)/pow(b,n2);
  cn3=(n3)*(n3-1)/pow(c,n3);
  df[0][0]=an1*pow(fabs(x),n1-2.0);
  df[0][1]=0;
  df[0][2]=0;
  df[1][0]=0;
  df[1][1]=bn2*pow(fabs(y),n2-2.0);
  df[1][2]=0;
  df[2][0]=0;
  df[2][1]=0;
  df[2][2]=cn3*pow(fabs(z),n3-2.0);
#if 0
  printf("n2=%.15G y=%.15G bn2=%f\n", n2, y, bn2);
  if (i==121)
    {
      printf("pow():%.15G n2-2=%.15G bn2=%.15G\n", pow(fabs(y),n2-2.0), n2-2, bn2);
      printf("calcfxx: df= %f %f %f\n", df[0][0], df[1][1], df[2][2]);
      printf("xyz=%.15G %.15G %.15G\n", x, y, z);
    }
#endif
}
#endif
void body2lab_fx(int i, double xp[], double x[], double **R)
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
void body2lab_fxx(int i, double fxxp[3][3], double fxx[3][3], double **Ri)
{
  int k1, k2, k3;
  double Rtmp[3][3];
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    Rtmp[k1][k2] += fxxp[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	fxx[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    fxx[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
extern void calc_Rdot(int i, double cosea[3], double sinea[3], double **Ro);

void calcFxtFtSE(int i, double x[3], double **RM, double cosea[3], double sinea[3], double Omega[3][3], 
	       double **R, double pos[3], double vel[3], double fxp[3], double fxxp[3][3],
	       double Fxt[3], double *Ft)
{
  double tRfxxRdot[3][3], tRfxxR[3][3], tRdotfx[3], tRfxx[3][3]; 
  double dx[3];
  int k1, k2, k3;
  if (isSymItens(i))
    {
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      Rdot[k1][k2] = 0.0;
	      for (k3 = 0; k3 < 3; k3++)
		{
		  Rdot[k1][k2] += -R[k1][k3]*Omega[k3][k2];
		}
	    }
	}	
    }
  else
    calc_Rdot(i, cosea, sinea, Rdot);
#if 1
  for (k1 = 0; k1 < 3; k1++)
    {
      tRdotfx[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRdotfx[k1] += Rdot[k2][k1]*fxp[k2];
	}
    }

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRfxx[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (R[k3][k1] == 0.0 || fxxp[k3][k2] == 0.0)
		continue;
	      tRfxx[k1][k2] += R[k3][k1]*fxxp[k3][k2];
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRfxxRdot[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (tRfxx[k1][k3] == 0.0 || Rdot[k3][k2] == 0.0)
		continue;
	      tRfxxRdot[k1][k2] += tRfxx[k1][k3]*Rdot[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRfxxR[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (tRfxx[k1][k3] == 0.0 || R[k3][k2] == 0.0)
		continue;
	      tRfxxR[k1][k2] += tRfxx[k1][k3]*R[k3][k2]; 
	    }
	}
    }
#if 0
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      DtX[k1][k2] = tRfxxRdot[k1][k2] + tRfxxR[k1][k2];
#endif  
#if 1
   for (k1 = 0; k1 < 3; k1++)
     dx[k1] = x[k1]-pos[k1];
#endif
   for (k1 = 0; k1 < 3; k1++)
    {
      Fxt[k1] = tRdotfx[k1];
      for (k2 = 0; k2 < 3; k2++)
	Fxt[k1] += tRfxxRdot[k1][k2]*dx[k2] - tRfxxR[k1][k2]*vel[k2]; 
     } 
   *Ft = 0;
   for (k1 = 0; k1 < 3; k1++)
     {
       for (k2 = 0; k2 < 3; k2++)
	 {
	   /* 24/07/07: il primo addendo è corretto il secondo no! */
	   /* 04/03/09: CHECK THIS FORMULA!!! <======================================= !!!!!!!!!!! 
	      10/03/09: la formula sembra corretta, forse il commento del 24/07/07 non era stato
	                aggiornato.*/

	   *Ft += -fxp[k1]*R[k1][k2]*vel[k2] + fxp[k1]*Rdot[k1][k2]*dx[k2];
	 }
     }
#endif
}
extern void UpdateOrient(int i, double ti, double **Ro, double Omega[3][3]);
void fdjacSE(int n, double x[], double fvec[], double **df, 
	   void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  int na; 
#ifdef EDHE_FLEX
  int typei, typej;
#endif
  double  rA[3], rB[3], ti, vA[3], vB[3];
  double OmegaA[3][3], OmegaB[3][3];
  double fx[3], gx[3];
  double Fxt[3], Gxt[3], Ft=0.0, Gt=0.0;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  int k1, k2;
  double xpA[3], xpB[3], fxxp[3][3], gxxp[3][3], fxx[3][3], gxx[3][3], fxp[3], gxp[3];

  ti = x[4] + (trefG - atomTime[iA]);
  rA[0] = rx[iA] + vx[iA]*ti;
  rA[1] = ry[iA] + vy[iA]*ti;
  rA[2] = rz[iA] + vz[iA]*ti;
  vA[0] = vx[iA];
  vA[1] = vy[iA];
  vA[2] = vz[iA];
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  if (isSymItens(iA))
    UpdateOrient(iA, ti, RA, OmegaA);
  else
    symtop_evolve_orient(iA, ti, RA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(iA, ti, RA, OmegaA);
#endif
  MD_DEBUG2(printf("i=%d ti=%f", iA, ti));
  MD_DEBUG2(print_matrix(RA, 3));
  na = (iA < Oparams.parnumA)?0:1;
  typei = typeOfPart[iA];
  ti = x[4] + (trefG - atomTime[iB]);
  rB[0] = rx[iB] + vx[iB]*ti + shift[0];
  rB[1] = ry[iB] + vy[iB]*ti + shift[1];
  rB[2] = rz[iB] + vz[iB]*ti + shift[2];
  vB[0] = vx[iB];
  vB[1] = vy[iB];
  vB[2] = vz[iB];
#ifdef MD_ASYM_ITENS
  if (isSymItens(iB))
    UpdateOrient(iB, ti, RB, OmegaB);
  else
    symtop_evolve_orient(iB, ti, RB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  UpdateOrient(iB, ti, RB, OmegaB);
#endif
  na = (iB < Oparams.parnumA)?0:1;
  typej = typeOfPart[iB];

  lab2body(iA, &x[0], xpA, rA, RA);
  lab2body(iB, &x[0], xpB, rB, RB);  
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfx(gxp, xpB[0], xpB[1], xpB[2], iB);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(gxxp, xpB[0], xpB[1], xpB[2], iB);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RA);
  body2lab_fx(iB, gxp, gx, RB);  
  body2lab_fxx(iA, fxxp, fxx, RA);
  body2lab_fxx(iB, gxxp, gxx, RB);
  MD_DEBUG37(printf("1)SE fx=%.15G %.15G %.15G\n", fx[0], fx[1], fx[2]));
  MD_DEBUG37(printf("1)SE gx=%.15G %.15G %.15G\n", gx[0], gx[1], gx[2]));
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
  	{
	  df[k1][k2] = fxx[k1][k2] + Sqr(x[3])*gxx[k1][k2];
	}
    }

  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = gx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = 2.0*x[3]*gx[k1];
    } 
  df[3][3] = 0.0;
  df[4][3] = 0.0;
  /* questi vanno ricalcolati per i superellissoidi!*/
  calcFxtFtSE(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], OmegaA, RA, rA, vA, fxp, fxxp, Fxt, &Ft);
  calcFxtFtSE(iB, x, RM[iB], cosEulAng[1], sinEulAng[1], OmegaB, RB, rB, vB, gxp, gxxp, Gxt, &Gt);
  
  //printf("1)[SE] Ft=%.15G Fxt=%.15G %.15G %.15G\n", Ft, Fxt[0], Fxt[1], Fxt[2]);
  //printf("1)[SE] Gt=%.15G Gxt=%.15G %.15G %.15G\n", Gt, Gxt[0], Gxt[1], Gxt[2]);
  /* -------------------------------------------- */
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][4] = Fxt[k1]+Sqr(x[3])*Gxt[k1]; 
    } 
  df[3][4] = Ft;
  df[4][4] = Gt;
#ifndef MD_GLOBALNR
  /* and now evaluate fvec */
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  /* WARNING: fix this for SE!!!! */
#if 0
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;
#endif
  fvec[3] = calcf(xpA, iA);
  fvec[4] = calcf(xpB, iB);
#endif
}
void fdjacDistNegSE(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int, int, double []), 
	       int iA, int iB, double shift[3], double *fx, double *gx)
{
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
  double rDC[3];
  double xpA[3], xpB[3], fxxp[3][3], gxxp[3][3], fxp[3], gxp[3], fxx[3][3], gxx[3][3];
  int k1, k2;
  /* ci mettiamo nel riferimento del corpo rigido dove lo Jacobiano
     assume la forma più semplice */
  lab2body(iA, &x[0], xpA, rA, RtA);
  lab2body(iB, &x[3], xpB, rB, RtB);  
#if 0
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
#endif
  MD_DEBUG37(printf("fx=%.15G %.15G %.15G gx=%.15G %.15G %.15G\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfx(gxp, xpB[0], xpB[1], xpB[2], iB);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(gxxp, xpB[0], xpB[1], xpB[2], iB);
#if 0
  printf("fxxp=%f %f %f %f %f %f %f %f %f\n", fxxp[0][0], fxxp[0][1], fxxp[0][2],
	 fxxp[1][0], fxxp[1][1], fxxp[1][2], fxxp[2][0], fxxp[2][1], fxxp[2][2]);
#endif
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fx(iB, gxp, gx, RtB);  
  body2lab_fxx(iA, fxxp, fxx, RtA);
  body2lab_fxx(iB, gxxp, gxx, RtB);
 
  MD_DEBUG37(printf("NEW fx=%.15G %.15G %.15G gx=%.15G %.15G %.15G\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  //printf("fxx[%d][%d]=%.15G gxx=%.15G\n", k1, k2, fxx[k1][k2], gxx[k1][k2]);
	  df[k1][k2] = fxx[k1][k2];
	  df[k1][k2+3] = Sqr(x[6])*gxx[k1][k2];
	  MD_DEBUG37(printf("fxx[%d][%d]:%.15g Xa:%.15G\n", k1, k2, fxx[k1][k2], 2.0*Xa[k1][k2]));
	  MD_DEBUG37(printf("gxx[%d][%d]:%.15g Xb:%.15G\n", k1, k2, gxx[k1][k2], 2.0*Xb[k1][k2]));
	}
    }
#if 1
  if (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = x[k1+3] - x[k1];
#ifdef EDHE_FLEX
	if (OprogStatus.targetPhi > 0.0)
	  {
	    axi[0] = axa[iA];
	    axi[1] = axb[iA];
	    axi[2] = axc[iA];
	    axj[0] = axa[iB];
	    axj[1] = axb[iB];
	    axj[2] = axc[iB];
	  }
	else
	  {
	    for (kk=0; kk < 3; kk++)
	      {
		axi[kk] = typesArr[typeOfPart[iA]].sax[kk];
		axj[kk] = typesArr[typeOfPart[iB]].sax[kk];
	      }
	  }
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axi[0],axi[1],axi[2])+max3(axj[0],axj[1],axj[2])))
	{
	  fdjac_disterr = 1;	
	}
#else
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axa[iA],axb[iA],axc[iA])+max3(axa[iB],axb[iB],axc[iB])))
	{
	  fdjac_disterr = 1;	
	}
#endif 
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
	    df[k1+5][k2] = 1.0 + x[7]*fxx[k1][k2];
	  else 
	    df[k1+5][k2] = x[7]*fxx[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1.0;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = fx[k1];
  MD_DEBUG37(printf("SUPERELLIPSE:\n"));
  MD_DEBUG37(print_matrix(df, 8));
#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
#if 0
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;
#endif
 MD_DEBUG37(printf("OLD)f=%.15G g=%.15G\n", fvec[3], fvec[4]));
 for (k1 = 0; k1 < 3; k1++)
   {
     fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
   }
 fvec[3] = calcf(xpA,iA);
 fvec[4] = calcf(xpB,iB);
 MD_DEBUG37(printf("NEW)f=%.15G g=%.15G\n", fvec[3], fvec[4]));
  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
  MD_DEBUG37(printf("SUPERELLIPSE fvec=%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4], fvec[5],
	 fvec[6], fvec[7]));
}
void fdjacDistNeg5SE(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int, int, double []), 
		   int iA, int iB, double shift[3], double *fx, double *gx)
{
  double rDC[3], rD[3], A[3][3], b[3], c[3];
  int k1, k2, k3;
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
  double xpA[3], xpB[3], fxxp[3][3], gxxp[3][3], fxp[3], gxp[3], fxx[3][3], gxx[3][3];
  lab2body(iA, &x[0], xpA, rA, RtA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fxx(iA, fxxp, fxx, RtA);
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      rD[k1] = x[k1] + fx[k1]*x[4];
    } 
  lab2body(iB, rD, xpB, rB, RtB);  
  calcfx(gxp, xpB[0], xpB[1], xpB[2], iB);
  calcfxx(gxxp, xpB[0], xpB[1], xpB[2], iB);
  body2lab_fx(iB, gxp, gx, RtB);  
  body2lab_fxx(iB, gxxp, gxx, RtB);
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{/* 02/03/09: rivedere bene qui!!! */
	  A[k1][k2] = gxx[k1][k2];
	  for (k3 = 0; k3 < 3; k3++) 
	    A[k1][k2] += x[4]*gxx[k1][k3]*fxx[k3][k2];
	  A[k1][k2] *= Sqr(x[3]);
	}
    }	
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = fxx[k1][k2] + A[k1][k2];
	}
    }
    //printf("rC: %f %f %f rD: %f %f %f\n", x[0], x[1], x[2], rD[0], rD[1], rD[2]);
  //printf("fx: %f %f %f x[4]: %f\n", fx[0], fx[1], fx[2], x[4]);
#if 1
  if (OprogStatus.SDmethod==2 || OprogStatus.SDmethod==3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = rD[k1] - x[k1];
#ifdef EDHE_FLEX
      if (OprogStatus.targetPhi > 0.0)
	{
	  axi[0] = axa[iA];
	  axi[1] = axb[iA];
	  axi[2] = axc[iA];
	  axj[0] = axa[iB];
	  axj[1] = axb[iB];
	  axj[2] = axc[iB];
	}
      else
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      axi[kk] = typesArr[typeOfPart[iA]].sax[kk];
	      axj[kk] = typesArr[typeOfPart[iB]].sax[kk];
	    }
	}
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axi[0],axi[1],axi[2])+max3(axj[0],axj[1],axj[2])))
	{
	  fdjac_disterr = 1;	
	}

#else
      if (scalProd(rDC, fx) < 0.0 && calc_norm(rDC) > (max3(axa[iA],axb[iA],axc[iA])+max3(axa[iB],axb[iB],axc[iB])))
	{
	  fdjac_disterr = 1;	
	}
#endif
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      b[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  b[k1] += gxx[k1][k2]*fx[k2];
	}
      b[k1] *= Sqr(x[3]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      c[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  c[k1] += x[4]*gx[k2]*fxx[k2][k1];
	}
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

#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = calcf(xpA,iA);
  fvec[4] = calcf(xpB,iB);
#endif
}

void fdjacNeighPlaneSE(int n, double x[], double fvec[], double **df, 
		     void (*vecfunc)(int, double [], double [], int), int iA)
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  double  rA[3], ti, vA[3], vB[3], OmegaB[3][3];
  double DA[3][3], fx[3], invaSqN, invbSqN, invcSqN;
  double Fxt[3], Ft;
  double xpA[3], fxxp[3][3], fxp[3], fxx[3][3];
#ifdef EDHE_FLEX
  int typei;
#endif
  double psi, phi;
  double OmegaA[3][3];
  int k1, k2;
  ti = x[4] + (trefG - atomTime[iA]);
  rA[0] = rx[iA] + vx[iA]*ti;
  rA[1] = ry[iA] + vy[iA]*ti;
  rA[2] = rz[iA] + vz[iA]*ti;
  vA[0] = vx[iA];
  vA[1] = vy[iA];
  vA[2] = vz[iA];
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  if (isSymItens(iA))
    UpdateOrient(iA, ti, RA, OmegaA);
  else
    symtop_evolve_orient(iA, ti, RA, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(iA, ti, RA, OmegaA);
#endif
  MD_DEBUG2(printf("i=%d ti=%f", iA, ti));
  MD_DEBUG2(print_matrix(RA, 3));
#ifdef EDHE_FLEX
  typei = typeOfPart[iA];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1/Sqr(axa[iA]);
      invbSqN = 1/Sqr(axb[iA]);
      invcSqN = 1/Sqr(axc[iA]);
    }
  else
    {
      invaSqN = 1/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1/Sqr(axa[iA]);
  invbSqN = 1/Sqr(axb[iA]);
  invcSqN = 1/Sqr(axc[iA]);
#endif
  //tRDiagR(iA, Xa, invaSqN, invbSqN, invcSqN, RA);
  MD_DEBUG2(printf("invabc: (%f,%f,%f)\n", invaSq, invbSq, invcSq));
  MD_DEBUG2(print_matrix(Xa, 3));
  DA[0][1] = DA[0][2] = DA[1][0] = DA[1][2] = DA[2][0] = DA[2][1] = 0.0;
  DA[0][0] = invaSqN;
  DA[1][1] = invbSqN;
  DA[2][2] = invcSqN;

  /*N.B. il corpo rigido B (ossia il piano) in tale caso non evolve! */
  vB[0] = 0.0;
  vB[1] = 0.0;
  vB[2] = 0.0;
  OmegaB[0][0] = 0;
  OmegaB[0][1] = 0;
  OmegaB[0][2] = 0;
  OmegaB[1][0] = 0;
  OmegaB[1][1] = 0;
  OmegaB[1][2] = 0;
  OmegaB[2][0] = 0;
  OmegaB[2][1] = 0;
  OmegaB[2][2] = 0;
   
  lab2body(iA, &x[0], xpA, rA, RA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RA);
  body2lab_fxx(iA, fxxp, fxx, RA);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = fxx[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = gradplane[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = -2.0*x[3]*gradplane[k1];
    } 
  df[3][3] = 0.0;
  df[4][3] = 0.0;
  calcFxtFtSE(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], OmegaA, RA, rA, vA, fxp, fxxp, Fxt, &Ft);
#if 0
#ifdef MD_ASYM_ITENS
  if (isSymItens(iA))
    calcFxtFtSym(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
  else
    calcFxtFt(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fx, Fxt, &Ft);
#else
  calcFxtFtSym(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
#endif
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][4] = Fxt[k1];
    } 
 df[3][4] = Ft;
 df[4][4] = 0.0;
#ifndef MD_GLOBALNRNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
 fvec[3] = calcf(xpA, iA);
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[4] += (x[k1]-rB[k1])*gradplane[k1];
   }
 MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
#endif
}

void funcs2beZeroedNeighPlaneSE(int n, double x[], double fvec[], int i)
{
  int na, k1, k2; 
  double  rA[3], ti;
  double fx[3];
  double invaSqN, invbSqN, invcSqN;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#else
  double Omega[3][3];
#endif
#ifdef EDHE_FLEX
  int typei;
#endif
  double xpA[3], fxp[3];

  /* x = (r, alpha, t) */ 
  ti = x[4] + (trefG - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  symtop_evolve_orient(i, ti, Rt, REtA, cosEulAng[0], sinEulAng[0], &phi, &psi);
#else
  UpdateOrient(i, ti, Rt, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  typei = typeOfPart[i];  
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSqN = 1.0/Sqr(axa[i]);
      invbSqN = 1.0/Sqr(axb[i]);
      invcSqN = 1.0/Sqr(axc[i]);
    }
  else
    {
      invaSqN = 1/Sqr(typesArr[typei].sax[0]);
      invbSqN = 1/Sqr(typesArr[typei].sax[1]);
      invcSqN = 1/Sqr(typesArr[typei].sax[2]);
    }
#else
  invaSqN = 1.0/Sqr(axa[i]);
  invbSqN = 1.0/Sqr(axb[i]);
  invcSqN = 1.0/Sqr(axc[i]);
#endif
  //tRDiagR(i, Xa, invaSqN, invbSqN, invcSqN, Rt);

  lab2body(i, &x[0], xpA, rA, Rt);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, Rt);

  /* il secondo ellissoide resta fermo al tempo iniziale */
  /* il secondo corpo rigido (ossia il piano) non deve evolvere nel tempo */
  MD_DEBUG(print_matrix(Xb,3));
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
  fvec[3] = calcf(xpA, i);
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[4] += (x[k1]-rB[k1])*gradplane[k1];
    }
#if 0
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//[k1]*x[7]; 
#endif
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.
  MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
}
#ifdef MD_FDJAC_SYM
void fdjacDistNegNeighPlaneSE(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3];
  int k1, k2;
  double xpA[3], fxxp[3][3], fxp[3], fxx[3][3];

  /* ci mettiamo nel riferimento del corpo rigido dove lo Jacobiano
     assume la forma più semplice */
  lab2body(iA, &(x[0]), xpA, rA, RtA);
#if 0
  printf("in fdjac\n");
  print_matrix(RtA,3);
  printf("rA=%f %f %f xpA=%f %f %f %f %f %f\n", rA[0], rA[1], rA[2], xpA[0], xpA[1], xpA[2], x[0], x[1], x[2]);
#endif
#if 0
  printf("\n x-rA=%f %f %f\n", x[0]-rA[0], x[1]-rA[1], x[2]-rA[2]);
  print_matrix(RtA,3);
  printf("============\n");
#endif
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fxx(iA, fxxp, fxx, RtA);
#if 0
  printf("fxx0: %f %f %f\n", fxx[0][0], fxx[0][1], fxx[0][2]);
  printf("fxx1: %f %f %f\n", fxx[1][0], fxx[1][1], fxx[1][2]);
  printf("fxx2: %f %f %f\n", fxx[2][0], fxx[2][1], fxx[2][2]); 
#endif
#if 0
  if (fabs(fxx[0][0])> 0)
    {
      printf("x=%f %f %f\n", x[0], x[1], x[2]);
      printf("fxx: %f %f %f\n", fxx[0][0], fxx[1][1], fxx[2][2]);
      printf("fxx: %f %f %f\n", fxxp[0][0], fxxp[1][1], fxxp[2][2]);
      print_matrix(RtA, 3);
      printf("xpA=%f %f %f\n", xpA[0], xpA[1], xpA[2]);
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  if (k1==k2)
	    {
	      df[k1][k2]=-1.0;
	      df[k1][k2+3]=1.0;
	    }
	  else
	    {
	      df[k1][k2]=0.0;
	      df[k1][k2+3]=0.0;
	    }
	    df[k1][k2] += -x[6]*fxx[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  for (k1 = 0; k1 < 5; k1++)
    {
      df[3][k1+3] = 0.0;
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = 0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1+3] = gradplane[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = -fx[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1.0; //+ 2.0*x[7]*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 0.0;//2.0*x[7]*Xa[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1.0;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = gradplane[k1];//fx[k1];
  //print_matrix(df,8);
#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
#if 0 
/* ma che cavolo ho scritto qua?!? */
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = x[k1+3] - x[k1] - x[6]*fx[k1];
    }
#else
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gradplane[k1];
    }
#endif
 fvec[3] = calcf(xpA, iA);
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
     fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
   }
 //fvec[4] = 0.5*fvec[4]-1.0;
 /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif

}
#else
void fdjacDistNegNeighPlaneSE(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3];
  int k1, k2;
  double xpA[3], fxxp[3][3], fxp[3], fxx[3][3];

  /* ci mettiamo nel riferimento del corpo rigido dove lo Jacobiano
     assume la forma più semplice */
  lab2body(iA, &(x[0]), xpA, rA, RtA);
#if 0
  printf("in fdjac\n");
  print_matrix(RtA,3);
  printf("rA=%f %f %f xpA=%f %f %f %f %f %f\n", rA[0], rA[1], rA[2], xpA[0], xpA[1], xpA[2], x[0], x[1], x[2]);
#endif
#if 0
  printf("\n x-rA=%f %f %f\n", x[0]-rA[0], x[1]-rA[1], x[2]-rA[2]);
  print_matrix(RtA,3);
  printf("============\n");
#endif
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fxx(iA, fxxp, fxx, RtA);
#if 0
  printf("fxx0: %f %f %f\n", fxx[0][0], fxx[0][1], fxx[0][2]);
  printf("fxx1: %f %f %f\n", fxx[1][0], fxx[1][1], fxx[1][2]);
  printf("fxx2: %f %f %f\n", fxx[2][0], fxx[2][1], fxx[2][2]); 
#endif
#if 0
  if (fabs(fxx[0][0])> 0)
    {
      printf("x=%f %f %f\n", x[0], x[1], x[2]);
      printf("fxx: %f %f %f\n", fxx[0][0], fxx[1][1], fxx[2][2]);
      printf("fxx: %f %f %f\n", fxxp[0][0], fxxp[1][1], fxxp[2][2]);
      print_matrix(RtA, 3);
      printf("xpA=%f %f %f\n", xpA[0], xpA[1], xpA[2]);
    }
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = fxx[k1][k2];
	  df[k1][k2+3] = 0;
	}
    }
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
      df[4][k1+3] = gradplane[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = -2.0*x[6]*gradplane[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1.0; //+ 2.0*x[7]*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 0.0;//2.0*x[7]*Xa[k1][k2];
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
    df[k1+5][7] = gradplane[k1];//fx[k1];
  //print_matrix(df,8);
#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gradplane[k1];
    }
 fvec[3] = calcf(xpA, iA);
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
     fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
   }
 //fvec[4] = 0.5*fvec[4]-1.0;
 /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif

}
#endif
void fdjacDistNegNeighPlane5SE(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3], rD[3];
  int k1, k2;
  double xpA[3], fxxp[3][3], fxp[3], fxx[3][3];

   /* ci mettiamo nel riferimento del corpo rigido dove lo Jacobiano
     assume la forma più semplice */
  lab2body(iA, &x[0], xpA, rA, RtA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2] ,iA);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fxx(iA, fxxp, fxx, RtA);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = fxx[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      rD[k1] = x[k1] + gradplane[k1]*x[4];
    } 
  //printf("rC: %f %f %f rD: %f %f %f\n", x[0], x[1], x[2], rD[0], rD[1], rD[2]);
  //printf("fx: %f %f %f x[4]: %f\n", fx[0], fx[1], fx[2], x[4]);
  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  df[3][3] = 0.0;
  df[3][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = gradplane[k1];
    } 
  df[4][3] = 0.0;
  df[4][4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    df[4][4] += gradplane[k1]*gradplane[k1];
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = -2.0*x[3]*gradplane[k1];
      df[k1][4] = 0.0;
    } 

#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
 //fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      //fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gradplane[k1];
   }
 fvec[3] = calcf(xpA, iA);
#endif

}

void funcs2beZeroedDistNegSE(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1; 
  double fx[3], gx[3], xpA[3], xpB[3], fxp[3], gxp[3];
  /* x = (r, alpha, t) */ 

  lab2body(i, &x[0], xpA, rA, RtA);
  lab2body(j, &x[3], xpB, rB, RtB);  
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  calcfx(gxp, xpB[0], xpB[1], xpB[2], j);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, RtA);
  body2lab_fx(j, gxp, gx, RtB);  
 
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }

  fvec[3] = calcf(xpA,i);
  fvec[4] = calcf(xpB,j);

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}

void funcs2beZeroedDistNeg5SE(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1; 
  double fx[3], gx[3], rD[3];
  /* x = (r, alpha, t) */ 
  double xpA[3], xpB[3], fxp[3], gxp[3];
  /* ci mettiamo nel riferimento del corpo rigido dove lo Jacobiano
     assume la forma più semplice */
  lab2body(i, &x[0], xpA, rA, RtA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, RtA);

  for (k1 = 0; k1 < 3; k1++)
    {
      rD[k1] = x[k1] + fx[k1]*x[4];
    }
  lab2body(j, rD, xpB, rB, RtB);  
  calcfx(gxp, xpB[0], xpB[1], xpB[2], j);
  body2lab_fx(j, gxp, gx, RtB);  
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = calcf(xpA,i);
  fvec[4] = calcf(xpB,j);
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;
#endif
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
#ifdef MD_FDJAC_SYM
void funcs2beZeroedDistNegNeighPlaneSE(int n, double x[], double fvec[], int i)
{
  int k1; 
  double fx[3];
  double xpA[3], fxp[3];

  lab2body(i, &x[0], xpA, rA, RtA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, RtA);

  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = x[k1+3] - x[k1] - x[6]*fx[k1];
    }

  fvec[3] = calcf(xpA, i);
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
    }
  //fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//fx[k1]*x[7]; 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
#else
void funcs2beZeroedDistNegNeighPlaneSE(int n, double x[], double fvec[], int i)
{
  int k1; 
  double fx[3];
  double xpA[3], fxp[3];

  lab2body(i, &x[0], xpA, rA, RtA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, RtA);

  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gradplane[k1];
    }
  fvec[3] = calcf(xpA, i);
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
    }
  //fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//fx[k1]*x[7]; 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
#endif

void funcs2beZeroedDistNegNeighPlane5SE(int n, double x[], double fvec[], int i)
{
  int k1, k2; 
  double fx[3], rD[3];
  double xpA[3], fxp[3];

  lab2body(i, &x[0], xpA, rA, RtA);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, RtA);
  for (k1 = 0; k1 < 3; k1++)
    {
      rD[k1] = x[k1] + gradplane[k1]*x[4];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[4] += (rD[k1]-rB[k1])*gradplane[k1];
    }
  fvec[3] = calcf(xpA, i);
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif

}
/* 03/03/2009 TODO: 1) make NNL work 
                    2) steepest descent for super-ellipsoids! */
double calc_sign_SE(int i, double *r, double **R, double *x, double **X)
{
  double a,b,c,e,n;
  int k1, k2;
  double xp[3], segno;
#ifdef MD_TWOPARAMSE
  a = typesArr[typeOfPart[i]].sax[0];
  b = typesArr[typeOfPart[i]].sax[1];
  c = typesArr[typeOfPart[i]].sax[2];
  e = typesArr[typeOfPart[i]].n[0];
  n = typesArr[typeOfPart[i]].n[1];
#endif
  if (!is_superellipse(i))
    {
      segno = -1;
      /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++) 
	  segno += (x[k1]-r[k1])*X[k1][k2]*(x[k2]-r[k2]); 
      return segno;
      MD_DEBUG37(printf("HE segno=%.15G\n", segno));
    }
#if 0
  for (k1=0; k1 < 3; k1++)
    x[k1] += L*rint((x[k1] - r[k1])/L);
#endif
  lab2body(i, x, xp, r, R);
#ifdef MD_TWOPARAMSE
  segno = pow(pow(fabs(xp[0])/a,2.0/e)+pow(fabs(xp[1])/b,2.0/e),e/n)+pow(fabs(xp[2])/c,2.0/n)-1.0; 
#else
  segno = calcf(xp, i);
#endif
  MD_DEBUG37(printf("SE segno = %.15G\n" , segno));

  return segno;
}
extern double costolSDgrad;
void funcs2beZeroedSE(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int na, k1; 
  double  rA[3], rB[3], ti;
  double fx[3], gx[3], xpA[3], xpB[3], fxp[3], gxp[3];
#ifdef EDHE_FLEX
  int typei, typej;
#endif
  double OmegaA[3][3], OmegaB[3][3];
#ifdef MD_ASYM_ITENS
  double cosea[3], sinea[3], phi, psi;
#endif
  /* x = (r, alpha, t) */ 
  ti = x[4] + (trefG - atomTime[i]);
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
#ifdef MD_ASYM_ITENS
  if (isSymItens(i))
    UpdateOrient(i, ti, RA, OmegaA);
  else
    symtop_evolve_orient(i, ti, RA, REt, cosea, sinea, &phi, &psi);
#else
  UpdateOrient(i, ti, Rt, Omega);
#endif
  na = (i < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typei = typeOfPart[i];
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
  else
    {
      invaSq[na] = 1/Sqr(typesArr[typei].sax[0]);
      invbSq[na] = 1/Sqr(typesArr[typei].sax[1]);
      invcSq[na] = 1/Sqr(typesArr[typei].sax[2]);
    }
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[i]);
      invbSq[na] = 1/Sqr(axbP[i]);
      invcSq[na] = 1/Sqr(axcP[i]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[i]);
      invbSq[na] = 1/Sqr(axb[i]);
      invcSq[na] = 1/Sqr(axc[i]);
    }
#endif
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = x[4] + (trefG - atomTime[j]);
  MD_DEBUG(printf("x[4]:%.15f atomTime[%d]:%.15f\n",x[4], j, atomTime[j]));
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
#ifdef MD_ASYM_ITENS
  if (isSymItens(j))
    UpdateOrient(j, ti, RB, OmegaB);
  else
    symtop_evolve_orient(j, ti, RB, REt, cosea, sinea, &phi, &psi);
#else
  UpdateOrient(j, ti, Rt, Omega);
#endif
  na = (j < Oparams.parnumA)?0:1;
#ifdef EDHE_FLEX
  na = 0;
  typej = typeOfPart[j];
  if (OprogStatus.targetPhi > 0.0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
  else
    {
      invaSq[na] = 1/Sqr(typesArr[typej].sax[0]);
      invbSq[na] = 1/Sqr(typesArr[typej].sax[1]);
      invcSq[na] = 1/Sqr(typesArr[typej].sax[2]);
    }
#elif defined(MD_POLYDISP)
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
  else
    {
      invaSq[na] = 1/Sqr(axaP[j]);
      invbSq[na] = 1/Sqr(axbP[j]);
      invcSq[na] = 1/Sqr(axcP[j]);
    }
#else
  if (OprogStatus.targetPhi > 0)
    {
      invaSq[na] = 1/Sqr(axa[j]);
      invbSq[na] = 1/Sqr(axb[j]);
      invcSq[na] = 1/Sqr(axc[j]);
    }
#endif
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  
  lab2body(i, &x[0], xpA, rA, RA);
  lab2body(j, &x[0], xpB, rB, RB);  
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  calcfx(gxp, xpB[0], xpB[1], xpB[2], j);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, fx, RA);
  body2lab_fx(j, gxp, gx, RB);  

  MD_DEBUG(print_matrix(Xb,3));
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = calcf(xpA, i);
  fvec[4] = calcf(xpB, j);
}
void calc_norm_SE(int i, double *x, double *n, double *r, double **R, double **X)
{
  double xp[3], fxp[3];
  int a, b;
  if (!is_superellipse(i))
    {
      for (a=0; a < 3; a++)
	{
	  n[a] = 0;
	  for (b = 0; b < 3; b++)
	    {
	      n[a] += -X[a][b]*(r[b]-x[b]);
	    }
	}
      return;
    }
  lab2body(i, x, xp, r, R);
  calcfx(fxp, xp[0], xp[1], xp[2], i);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(i, fxp, n, R);
}
/* steepest descent per super-ellissoidi e piani*/
double gradcgfuncRyckNNLSE(double *vec, double *grad, double *fx, double *signA)
{
  int kk, k1; 
  double K1, K2, F, nf, ng, dd[3], normdd, ngA, ngB;
  double S=1.0, A=1.0, B, gradfx, gradgx, gx[3];
  doneryck = 0;
  calc_norm_SE(icg, vec, fx, rA, RtA, Xa);

  for (k1 = 0; k1 < 3; k1++)
    {
      dd[k1] = vec[k1+3]-vec[k1];
    }
  /* NOTA 25/03/09: calcdist ha come parametro forceguess che in alcuni casi puo' essere 1 ed in altri 0 
     usando quindi OprogStatus.forceguess non si tiene conto di ciò!! CONTROLLARE!!! */
  for (k1 = 0; k1 < 3; k1++)
    gx[k1] = gradplane[k1];
#if 0
  if (OprogStatus.forceguess)
    {
      A = calc_sign_SE(icg, rA, R[icg],&(vec[3]), Xa); 
      //B = calc_sign_SE(jcg, rB, R[jcg],&(vec[0]), Xb);
      if (A<0 && B<0)
	S = -1.0;
      else
	S = 1.0;
    }
#endif
  normdd = calc_norm(dd);
  if (normdd==0)
    {
      doneryck = 2;
      return 0;
    }
  /* la norma dei gradienti e' sempre stepSDA e stepSDB*/ 
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
#if 1
  gradfx = 0;
  gradgx = 0;
  for (k1=0; k1 < 3; k1++)
    {
      gradfx += grad[k1]*fx[k1]; 
      gradgx += grad[k1+3]*gx[k1];
    }
  //printf("gradfx=%.15G gradgx=%.15G nf=%.15G ng=%.15G\n", gradfx, gradgx, nf, ng);
  for (kk=0; kk < 3; kk++)
    {
      grad[kk] -= gradfx*fx[kk];
      grad[kk+3] -= gradgx*gx[kk];
    }
#endif
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
	  if (sqrt(ngA) < OprogStatus.tolSDgrad*(icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB))
	      //&& sqrt(ngB) < OprogStatus.tolSDgrad*(jcg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB))
	    {
	      //accngA++;
	      //accngB++;
	      doneryck = 1;
	    }
	}
      else
	{
	  if (fabs(scalProd(fx,dd)) > normdd*costolSDgrad && fabs(scalProd(gx,dd)) > normdd*costolSDgrad)
	    {
	      doneryck = 1;
	    }
	}
    }
  S *= OprogStatus.springkSD;
  F = S*Sqr(normdd);
  return F; 
}
/* =========================== >>> forces <<< ======================= */
double cgfuncRyckNNLSE(double *vec)
{
  int kk;
  double A, B, F;
  /* NOTA 25/03/09: calcdist ha come parametro forceguess che in alcuni casi puo' essere 1 ed in altri 0 
     usando quindi OprogStatus.forceguess non si tiene conto di ciò!! CONTROLLARE!!! */
#if 0
  if (OprogStatus.forceguess)
    {
      A = calc_sign_SE(icg, rA, R[icg],&(vec[3]), Xa); 
      B = calc_sign_SE(jcg, rB, R[jcg],&(vec[0]), Xb);
      if (A<0 && B<0)
	A = -OprogStatus.springkSD;
      else
	A = OprogStatus.springkSD;
    }
  else
#endif
    A = OprogStatus.springkSD;
 
  F = 0.0;
  for (kk=0; kk < 3; kk++)
    F += A*Sqr(vec[kk]-vec[kk+3]);
  return F;
}
/* steepest descent per super-ellissoidi */
double gradcgfuncRyckSE(double *vec, double *grad, double *fx, double *gx, double *signA, double *signB)
{
  int kk, k1; 
  double K1, K2, F, nf, ng, dd[3], normdd, ngA, ngB;
  double S=1.0, A=1.0, B, gradfx, gradgx;
  doneryck = 0;
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += Xa[k1][k2]*(vec[k2] - rA[k2]);
      fx[k1] *= 2.0;
    }
#endif
  calc_norm_SE(icg, vec, fx, rA, RtA, Xa);

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += Xb[k1][k2]*(vec[k2+3] - rB[k2]);
      gx[k1] *= 2.0;
#endif
      dd[k1] = vec[k1+3]-vec[k1];
    }
  calc_norm_SE(jcg, &(vec[3]), gx, rB, RtB, Xb);
  /* NOTA 25/03/09: calcdist ha come parametro forceguess che in alcuni casi puo' essere 1 ed in altri 0 
     usando quindi OprogStatus.forceguess non si tiene conto di ciò!! CONTROLLARE!!! */
  if (OprogStatus.forceguess)
    {
#if 0
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
#endif
      A = calc_sign_SE(icg, rA, R[icg],&(vec[3]), Xa); 
      B = calc_sign_SE(jcg, rB, R[jcg],&(vec[0]), Xb);
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
  /* la norma dei gradienti e' sempre stepSDA e stepSDB*/ 
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
#if 0
  gradfx = 0;
  gradgx = 0;
  for (k1=0; k1 < 3; k1++)
    {
      gradfx += grad[k1]*fx[k1]; 
      gradgx += grad[k1+3]*gx[k1];
    }
  //printf("gradfx=%.15G gradgx=%.15G nf=%.15G ng=%.15G\n", gradfx, gradgx, nf, ng);
  for (kk=0; kk < 3; kk++)
    {
      grad[kk] -= gradfx*fx[kk];
      grad[kk+3] -= gradgx*gx[kk];
    }
#endif
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
#if 0
	  printf("normdd=%.15G\n", normdd);
	  printf("rAB=%.15G\n", sqrt(Sqr(rA[0]-rB[0]) + Sqr(rA[1]-rB[1]) + Sqr(rA[2]-rB[2])));
#endif
	  if (fabs(scalProd(fx,dd)) > normdd*costolSDgrad && fabs(scalProd(gx,dd)) > normdd*costolSDgrad)
	    {
#if 0
	      printf("fx+qlpha^2*gx=%.15G\n", fx[0]*nf+Sqr(nf/ng)*gx[0]*ng);

	      printf("fx=%.15G %.15G %.15G gx=%.15G %.15G %.15G\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]);
	      printf("qui scalprod(fx,dd)=%.15G scalprod(gx,dd)=%.15G normdd=%.15G costolSDgrad=%.15G\n", fabs(scalProd(fx,dd))/normdd, fabs(scalProd(gx,dd))/normdd,normdd,costolSDgrad);
#endif
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
double cgfuncRyckSE(double *vec)
{
  int kk;
  double A, B, F;
  if (OprogStatus.forceguess)
    {
      A = calc_sign_SE(icg, rA, R[icg],&(vec[3]), Xa); 
      B = calc_sign_SE(jcg, rB, R[jcg],&(vec[0]), Xb);
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
extern long long int itsfrprmn, callsfrprmn,callsok, callsprojonto, itsprojonto;

extern double sfA, sfB;
extern double shiftcg[3], lambdacg, minaxicg, minaxjcg;
extern double gradfG[3], gradgG[3], dxG[6];
extern double max3(double a, double b, double c);
extern double min3(double a, double b, double c);
extern double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB;
extern int cghalfspring;
double calcfLab(int i, double *x, double *rA, double **Ri)
{
  double xp[3];
  lab2body(i, &x[0], xp, rA, Ri);
  //printf("rA = %.15G %.15G %.15G\n", rA[0], rA[1], rA[2]);
  //printf("xp=%f %f %f\n", xp[0], xp[1], xp[2]);
  //printf("Ri=%f %f %f %f %f %f %f %f %f\n", Ri[0][0], Ri[0][1], Ri[0][1], Ri[1][0], Ri[1][1], Ri[1][2],
  //  Ri[2][0], Ri[2][1], Ri[2][2]);
  return calcf(xp, i);
}
double func_to_zero(double chsi, int i, double *x, double *n, double *r, double **Ri)
{
  double x1[3], p[3], A;
  int kk;
  /* n è la normale al punto x relativamente al super-ellissoide
     i, che ha il centro di massa posizionato in r a ha orientazione Ri.
     Se si tratta di un ellissoide X non è altro che la matrice Xa "ruotata" di Ri */
  ///printf("f(x)=%.15G\n", calcfLab(i, ri, r, Ri));
  /* N.B. per ora ho scelto di cercare di portare il glider sulla superficie
     lungo la direzione individuata dal centro del super-ellissoide e dal punto x+dr.
     Per superfici convesse esiste sempre una soluzione di tale problema.
     In alternativa, come ho già fatto per gli ellissoidi, si puo' usare la direzione individuata
     dalla normale nel punto x tuttavia in tale caso la soluzione non esiste sempre, anche se
     con un passo sufficientemente piccolo esiste anche per superfici concave.  */

  for (kk=0; kk < 3; kk++)
    x1[kk] = x[kk] - chsi*n[kk];
    //x1[kk] = p[kk] - chsi*(p[kk]-r[kk]);
  
    ///printf("sf=%f x1=%f %f %f chsi=%f p=%f %f %f r=%f %f %f dr=%f %f %f\n", sf, x1[0], x1[1], x1[2], chsi, p[0], p[1], p[2], r[0], r[1], r[2], dr[0], dr[1], dr[2]);
  A = calcfLab(i, x1, r, Ri);	
  ////printf("chsi=%.15G f(x1)=%.15G f(x)=%.15G\n", chsi, A, calcfLab(i, x, r, Ri));
  //printf("f(r1)=%.15G \n", calcfLab(i, ri, r, Ri));
  return A;
}
int iSE;
double *xSE, *nSE, *rSE, **RiSE;
double func_to_zero_zb(double chsi)
{
  return func_to_zero(chsi, iSE, xSE, nSE, rSE, RiSE);
}
extern int zbrac(double (*func)(double), double *x1, double *x2);

int find_surf_sol(int i, double *x, double *n, double *r, double *sol, double **Ri, double *dr, double dx)
{
  int kk;
  double chsi1, chsi2;
  /* unidimensional root finding (NR or Brent? see Numerical Recipe to make a decision) */
  /* Evaluate the initial value of chsi according to sign of point x with respect to the surface of the SE*/
#if 0
  chsi = calc_sign_SE(i, r, Ri, x, X);
  if (chsi < 0)
    {
      chsi1 = chsi;
      chsi2 = -chsi;
    }
  else
    {
      chsi1 = chsi;
      chsi2 = -chsi;
    }
#else
  chsi1 = -dx;
  chsi2 = dx;
#endif
  iSE = i;
  rSE = r;
  RiSE = Ri;
  nSE = n;
  xSE = x;
#if 1
  if (!zbrac(func_to_zero_zb, &chsi1, &chsi2))
    {
      return 0;
    }
#endif
   //riSE= ri;
  //printf("func_to_zero_zb(0)=%.15G func_to_zero_zb(1):%.15G\n", func_to_zero_zb(0), func_to_zero_zb(1));
  *sol = zbrent(func_to_zero_zb, chsi1, chsi2, 1E-8);
  if (polinterr)
    return 0;
  else 
    return 1;
}
/* questa e' la routine piu' complicata da riscrivere per i super-ellissoidi */ 
extern int accngA, accngB;
void calc_constr_force(int i, double *x, double *n, double dt2, double *g, double *r, double **Ri, double dx)
{
  int kk, done=0, ret;
  const double GOLD=1.618034;
  double r1A[3];
  double sol=0.0;
 
  callsprojonto++;
  //printf("dr=%f %f %f\n", dr[0], dr[1], dr[2]);
  ret = find_surf_sol(i, x, n, r, &sol, Ri, g, dx);
  ////printf("i=%d SE sol=%.15G ret=%d\n", i, sol, ret);
  for (kk = 0; kk < 3; kk++)
    {
      g[kk] = sol*n[kk]/dt2; 
      //dr[kk] = sf*dr[kk]-sol*(r1[kk]-rA[kk]); 
    }
  /* commentando questa riga il valore di sf usato per rimanere "aderenti" alla superficie
   * non viene mantenuto.
   * In tal modo il passo non puo' decrescere in maniera irreversibile se non intorno al minimo. */
}
int check_doneSE(double fp, double fpold, double minax)
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
void check_are_on_surf(char* str, int i, int j, double *ri, double *rj)
{
  double A, B;
  A=  calc_sign_SE(i, rA, RtA, ri, Xa);
  B=  calc_sign_SE(j, rB, RtB, rj, Xb);
  printf("[%s] (%d,%d) calcf(ri)=%.15G calcf(rj)=%.15G\n",str, i, j, A, B);

}

#if 0
void projectgradSE(double *p, double *xi, double *gradf, double *gradg)
{
  int kk;
  double dist;
  dist = 0;
  for (kk=0; kk < 3; kk++)
    {
      dist+=Sqr(p[kk+3]-p[kk]);
    }
  dist = sqrt(dist);
  //check_are_on_surf("SD LOOP BEF P", icg, jcg, p, &(p[3]));
  projontoSE(icg, p, xi, rA, Xa, RtA, gradf, &sfA, dist);
  projontoSE(jcg, &p[3], &xi[3], rB, Xb, RtB, gradg, &sfB, dist);

#if 0
    {
      double pp[6];
      for (kk=0; kk < 3; kk++)
	{
	  pp[kk] = p[kk]+xi[kk];
	  pp[kk+3] = p[kk+3] + xi[kk+3];
	}
      check_are_on_surf("SD LOOP AFT P", icg, jcg, pp, &(pp[3]));
    }
#endif
}
#endif
void frprmnRyckNNLSE(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), double (*dfunc)(double [], double [], double [], double []))
  /*Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func,
   * using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is
   * input as ftol. Returned quantities are p (the location of the minimum), iter
   * (the number of iterations that were performed), and fret (the minimum value of the function).
   * The routine linmin is called to perform line minimizations. */
{ 
  int j,its;
  const int ITMAXFR = OprogStatus.maxitsSD;
  //const double GOLD=1.618034;
  double dx, dt2, fp, fpold=0.0, signA, signB;
  double minax, xi[6], xiold[6], g[6];
  double signAold, signBold, pold[6], poldold[6];
  //printf("primaprima p= %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2], p[3], p[4], p[5]);
 
  minax = min(minaxicg,minaxjcg);
  sfA = icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
  callsfrprmn++;
  /*Initializations.*/
  ////check_are_on_surf("SD BEGIN", icg, jcg, p, &(p[3]));
  fp = (*dfunc)(p,xi,gradfG, &signA); 
  if (doneryck==2)
    {
      callsok++;
      return;
    }
#if 1
  if ((OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 4) &&
      check_doneSE(fp, fpold, minax))
    {
      //printf("BOH\n");
      callsok++;
      return;
    }
#endif
  //printf("p=%f %f %f xi=%f %f %f gradfG=%f %f %f gradgG=%f %f %f\n", p[0], p[1], p[2],
  //	 xi[0], xi[1], xi[2], gradfG[0], gradfG[1], gradfG[2], gradgG[0], gradgG[1], gradgG[2]);

  //projectgradSE(p,xi,gradfG,gradgG);  
  //check_are_on_surf("SD BEGIN", icg, jcg, p, &(p[3]));
  dt2 = Sqr(sfA);
  for (j=0; j < n; j++)
    {
      pold[j] = p[j];
    }	  
  for (its=1;its<=ITMAXFR;its++)
    { 
      itsfrprmn++;      
      *iter=its;
      for (j=0; j < n; j++)
	{
	  poldold[j] = pold[j];
	  pold[j] = p[j];
	  xiold[j] = xi[j];
	  /* uncostrained move (verlet algorithm) */
	  //p[j] = 2.0*pold[j] - poldold[j] + dt2*xi[j];
	  p[j] = pold[j] + dt2*xi[j];
	}
      if (0)
	{
	  double ddd[3];
	  int kk;
	  for (kk=0; kk < 3; kk++)
	    ddd[kk] = p[kk]-rA[kk];
	  printf("n.ddd=%.15G\n", scalProd(gradfG, ddd));
	  printf("its=%d p %f %f %f pold %f %f %f\n", its, p[0], p[1], p[2], pold[0], pold[1], pold[2]);
	  printf("rA %f %f %f\n", rA[0], rA[1], rA[2]);
	  printf("n %f %f %f\n", gradfG[0], gradfG[1], gradfG[2]); 
	}
      dx = calc_norm(&(xi[0]))*dt2;
      calc_constr_force(icg, &(p[0]), gradfG, dt2, &(g[0]), rA, RtA, dx); 
      //dx = calc_norm(&(xi[3]))*dt2;
      //calc_constr_force(jcg, &(p[3]), gradgG, dt2, &(g[3]), rB, RtB, dx);
      /* NOTA 25/03/09: la funzione seguente non serve in quanto è sufficente proiettare lo spostamento
	 sul piano per far si che il punto appartenga ad esso */
      //calc_constr_force_onplane(&(p[3]), gradplane, dt2, &(g[3]), rB);
      for (j=0; j < 3; j++)
	p[j] += -g[j]*dt2;

      ////check_are_on_surf("SD LOOP", icg, jcg, p, &(p[3]));
      //printf("A) distance(its=%d)=%.15G\n", its, sqrt(Sqr(p[0]-p[3])+Sqr(p[1]-p[4])+Sqr(p[2]-p[5])));
      signAold = signA;
      fpold = fp; 
      fp = (*dfunc)(p,xi,gradfG, &signA);
      if (doneryck==2)
	{
	  callsok++;
	  return;
	 }
      if (check_doneSE(fp, fpold, minax))
	{
	  callsok++;
	  return;
	}
    } 
#ifndef MD_CALC_VBONDING
  printf("[Steepest Descent NNL LOOP] TOO MANY ITERATIONS\n");
#endif
  return; 
  //nrerror("Too many iterations in frprmn");
  
}
void frprmnRyckSE(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), double (*dfunc)(double [], double [], double [], double [], double*, double*))
  /*Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func,
   * using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is
   * input as ftol. Returned quantities are p (the location of the minimum), iter
   * (the number of iterations that were performed), and fret (the minimum value of the function).
   * The routine linmin is called to perform line minimizations. */
{ 
  int j,its;
  const int ITMAXFR = OprogStatus.maxitsSD;
  //const double GOLD=1.618034;
  double dx, dt2, fp, fpold=0.0, signA, signB;
  double minax, xi[6], xiold[6], g[6];
  double signAold, signBold, pold[6], poldold[6];
  //printf("primaprima p= %.15G %.15G %.15G %.15G %.15G %.15G\n", p[0], p[1], p[2], p[3], p[4], p[5]);
 
  minax = min(minaxicg,minaxjcg);
  sfA = icg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
  sfB = jcg<Oparams.parnumA?OprogStatus.stepSDA:OprogStatus.stepSDB;
  callsfrprmn++;
  /*Initializations.*/
  ////check_are_on_surf("SD BEGIN", icg, jcg, p, &(p[3]));
  fp = (*dfunc)(p,xi,gradfG,gradgG, &signA, &signB); 
  if (doneryck==2)
    {
      callsok++;
      return;
    }
#if 1
  if ((OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 4) &&
      check_doneSE(fp, fpold, minax))
    {
      callsok++;
      return;
    }
#endif
  //printf("p=%f %f %f xi=%f %f %f gradfG=%f %f %f gradgG=%f %f %f\n", p[0], p[1], p[2],
  //	 xi[0], xi[1], xi[2], gradfG[0], gradfG[1], gradfG[2], gradgG[0], gradgG[1], gradgG[2]);

  //projectgradSE(p,xi,gradfG,gradgG);  
  //check_are_on_surf("SD BEGIN", icg, jcg, p, &(p[3]));
  dt2 = Sqr(sfA);
  for (j=0; j < n; j++)
    {
      pold[j] = p[j];
    }	  
  for (its=1;its<=ITMAXFR;its++)
    { 
      itsfrprmn++;      
      *iter=its;
      for (j=0; j < n; j++)
	{
	  poldold[j] = pold[j];
	  pold[j] = p[j];
	  xiold[j] = xi[j];
	  /* uncostrained move (verlet algorithm) */
	  //p[j] = 2.0*pold[j] - poldold[j] + dt2*xi[j];
	  p[j] = pold[j] + dt2*xi[j];
	}
      if (0)
	{
	  double ddd[3];
	  int kk;
	  for (kk=0; kk < 3; kk++)
	    ddd[kk] = p[kk]-rA[kk];
	  printf("n.ddd=%.15G\n", scalProd(gradfG, ddd));
	  printf("its=%d p %f %f %f pold %f %f %f\n", its, p[0], p[1], p[2], pold[0], pold[1], pold[2]);
	  printf("rA %f %f %f\n", rA[0], rA[1], rA[2]);
	  printf("n %f %f %f\n", gradfG[0], gradfG[1], gradfG[2]); 
	}
      dx = calc_norm(&(xi[0]))*dt2;
      calc_constr_force(icg, &(p[0]), gradfG, dt2, &(g[0]), rA, RtA, dx); 
      dx = calc_norm(&(xi[3]))*dt2;
      calc_constr_force(jcg, &(p[3]), gradgG, dt2, &(g[3]), rB, RtB, dx);
      for (j=0; j < n; j++)
	p[j] += -g[j]*dt2;

      ////check_are_on_surf("SD LOOP", icg, jcg, p, &(p[3]));
      //printf("A) distance(its=%d)=%.15G\n", its, sqrt(Sqr(p[0]-p[3])+Sqr(p[1]-p[4])+Sqr(p[2]-p[5])));
      signAold = signA;
      signBold = signB;
      fpold = fp; 
      fp = (*dfunc)(p,xi,gradfG, gradgG, &signA, &signB);
      if (doneryck==2)
	{
	  //printf("QUIIII AAA\n");
	  callsok++;
	  return;
	 }
      if (check_doneSE(fp, fpold, minax))
	{
	  //printf("QUIIII BBB\n");
	  callsok++;
	  return;
	}
    } 
#ifndef MD_CALC_VBONDING
  printf("[Steepest Descent LOOP] TOO MANY ITERATIONS\n");
#endif
  return; 
  //nrerror("Too many iterations in frprmn");
  
}
void distSDSupEll(int i, int j, double shift[3], double *vecg, double lambda, int halfspring)
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
  if (OprogStatus.targetPhi > 0.0)
    {
      axaiF = axa[i];
      axbiF = axb[i];
      axciF = axc[i];
    }
  else
    {
      axaiF = typesArr[typei].sax[0];
      axbiF = typesArr[typei].sax[1];
      axciF = typesArr[typei].sax[2];
    }
  minaxicg = min3(axaiF, axbiF, axciF);
  if (OprogStatus.targetPhi > 0.0)
    {
      axajF = axa[j];
      axbjF = axb[j];
      axcjF = axc[j];
    }
  else
    {
      axajF = typesArr[typej].sax[0];
      axbjF = typesArr[typej].sax[1];
      axcjF = typesArr[typej].sax[2];
    }  
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
  cghalfspring = halfspring;
  for (kk=0; kk < 3; kk++)
    {
      shiftcg[kk] = shift[kk];
    }
  for (kk=0; kk < 6; kk++)
    {
      vec[kk] = vecg[kk];
    }
  
  frprmnRyckSE(vec, 6, OprogStatus.tolSD, &iter, &Fret, cgfuncRyckSE, gradcgfuncRyckSE);
  for (kk=0; kk < 6; kk++)
    {
      vecg[kk] = vec[kk];
    }
}
void distSD_NNL(int i, double *vecg, double lambda, int halfspring)
{
  int kk;
  double Fret;
  int iter;
  double vec[8];
#ifdef EDHE_FLEX
  int typei;
  double axaiF, axbiF, axciF, axajF, axbjF, axcjF;
#endif
  icg = i;
#ifdef EDHE_FLEX
  typei = typeOfPart[i];
  if (OprogStatus.targetPhi > 0.0)
    {
      axaiF = axa[i];
      axbiF = axb[i];
      axciF = axc[i];
    }
  else
    {
      axaiF = typesArr[typei].sax[0];
      axbiF = typesArr[typei].sax[1];
      axciF = typesArr[typei].sax[2];
    }
  minaxicg = min3(axaiF, axbiF, axciF);
#endif
  cghalfspring = halfspring;
  for (kk=0; kk < 6; kk++)
    {
      vec[kk] = vecg[kk];
    }
  frprmnRyckNNLSE(vec, 6, OprogStatus.tolSD, &iter, &Fret, cgfuncRyckNNLSE, gradcgfuncRyckNNLSE);
  for (kk=0; kk < 6; kk++)
    {
      vecg[kk] = vec[kk];
    }
}

double func_to_zero_intersec(double chsi, int i, double *x, double *r, double **Ri)
{
  double x1[3], p[3], A;
  int kk;
  for (kk=0; kk < 3; kk++)
    x1[kk] = r[kk] + chsi*(x[kk]-r[kk]);
  A = calcfLab(i, x1, r, Ri);	
  return A;
}
double func_to_zero_intersec_zb(double chsi)
{
  return func_to_zero_intersec(chsi, iSE, xSE, rSE, RiSE);
}

int calc_intersecSE(int i, double *rB, double *rA, double **Ri, double* rI)
{
  double chsi1, chsi2, chsi;
  int kk;
#if 0
  chsi = calc_sign_SE(i, r, Ri, x, X);
  if (chsi < 0)
    {
      chsi1 = chsi;
      chsi2 = -chsi;
    }
  else
    {
      chsi1 = chsi;
      chsi2 = -chsi;
    }
#else
  chsi1 = 0.0;
  chsi2 = 1.0;
#endif
  iSE = i;
  rSE = rA;  /* centro del SE */
  RiSE = Ri; /* orientazione del SE */
  xSE = rB;  /* a partire da rA, rArB è la direzione lungo la quale calcolare il punto d'intersezione 
	        con il SE*/
  //printf("func_to_zero_zb(0)=%.15G func_to_zero_zb(1):%.15G\n", func_to_zero_zb(0), func_to_zero_zb(1));
  chsi = zbrent(func_to_zero_intersec_zb, chsi1, chsi2, 1.0E-14);
  //printf("i=%d chsi=%.15G\n", i, chsi);
  for (kk=0; kk < 3; kk++)
    rI[kk] = rA[kk] + chsi*(rB[kk]-rA[kk]);
  //printf("polinterr=%d intersec: %.15G\n", polinterr, calcfLab(i, rI, rA, Ri));
  
  if (polinterr)
    return 0;
  else 
    return 1;
}
void calcfxLabSE(int i, double *x, double *r, double **Ri, double fx[3])
{
  /* r = centro SE, Ri = orientazione SE, x=punto sulla superficie di cui si deve calcolare il gradiente */ 
  double fxp[3], xpA[3];
  lab2body(i, &x[0], xpA, r, Ri);
  calcfx(fxp, xpA[0], xpA[1], xpA[2], i);
  body2lab_fx(i, fxp, fx, Ri);
}
/* ======================================== */
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define KRI 5
#define FUNC(x) ((*func)(x))
int polinterrSQ, growthSQ;
void polintSQ(double xain[KRI], double yain[KRI], int n, double x, double *y, double *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y,
 * and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xai) = yai, 
 * i = 1, . . . , n, then the returned value y = P(x).*/
{ 
  int i,m,ns=1; 
  double den,dif,dift,ho,hp,w, xa[KRI+1], ya[KRI+1];
  double c[KRI+1], d[KRI+1];
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
	      polinterrSQ=1;
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
}
double trapzdInt(double (*func)(double), double a, double b)
/* This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
   as a pointer to the function to be integrated between limits a and b, also input. When called with
   n=1, the routine returns the crudest estimate of b
   a f(x)dx. Subsequent calls with n=2,3,...
   (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points. */
{
  double sum,dx,x1,x2;
  //static double s;
  int it,j, n;
  n=100;
  dx = (a-b)/n;
  sum = 0;
  for (j = 0; j < n-1; j++)
    { 
      x1 = a+dx*j;
      x2 = x1 + dx; 
      sum+=dx*0.5*((*func)(x1)+(*func)(x2));
    }
  return sum;
}


double trapzd(double (*func)(double), double a, double b, int n)
/* This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
   as a pointer to the function to be integrated between limits a and b, also input. When called with
   n=1, the routine returns the crudest estimate of b
   a f(x)dx. Subsequent calls with n=2,3,...
   (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points. */
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  if (n == 1) 
    {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } 
  else 
    {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm; /* This is the spacing of the points to be added.*/
	x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
      s=0.5*(s+(b-a)*sum/tnm); /* This replaces s by its refined value. */
	return s;
    }
}
/*Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate;
  JMAX limits the total number of steps; K is the number of points used in the extrapolation.*/
double qromb(double (*func)(double), double a, double b)
/* Returns the integral of the function func from a to b. Integration is performed by Romberg’s
   method of order 2K, where, e.g., K=2 is Simpson’s rule. */
{
  double ss,dss=0.0;
  double s[JMAXP],h[JMAXP+1]; /* These store the successive trapezoidal approximations */
  int j; /*and their relative stepsizes.*/
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) 
    {
      s[j]=trapzd(func,a,b,j);
      if (j >= KRI) 
	{
	  polinterrSQ=0;
	  polintSQ(&h[j-KRI],&s[j-KRI],KRI,0.0,&ss,&dss);
	  if (fabs(dss) <= EPS*fabs(ss)) 
	    return ss;
	  if (polinterrSQ)
	    {
	      printf("interpolation error\n");
	      exit(-1);
	    }
	}
      h[j+1]=0.25*h[j];
      /* This is a key step: The factor is 0.25 even though the stepsize is decreased by only
	 0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1),
	 not just a polynomial in h.*/
    }
  printf("Too many steps in routine qromb");
  exit(-1);
  return 0.0; 
}
double *volSQ;
#define GAULEG_EPS 3.0e-15 /* EPS is the relative precision.*/
double *qgaus_w, *qgaus_x, *gauleg_w, *gauleg_x;
void gauleg(double x1, double x2, double *x, double* w, int n)
/* Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
   arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
   Legendre n-point quadrature formula.*/
{
  int m,j,i;
  double pi, z1,z,xm,xl,pp,p3,p2,p1; /*High precision is a good idea for this routine.*/
  pi = 2.0*acos(0.0);
  m=(n+1)/2; /* The roots are symmetric in the interval, so*/
  xm=0.5*(x2+x1); /* we only have to find half of them. */
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) { /* Loop over the desired roots.*/
    //z=cos(3.141592654*(i-0.25)/(n+0.5));
      z=cos(pi*(i-0.25)/(n+0.5)); 
    /* Starting with the above approximation to the ith root, we enter the main loop of
       refinement by Newton’s method. */
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) { /* Loop up the recurrence relation to get the*/
	p3=p2; /* Legendre polynomial evaluated at z.*/
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      /*p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
	by a standard relation involving also p2, the polynomial of one lower order.
       */
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp; /* Newton’s method.*/
    } while (fabs(z-z1) > GAULEG_EPS);
    x[i]=xm-xl*z; /* Scale the root to the desired interval,*/
    x[n+1-i]=xm+xl*z; /* and put in its symmetric counterpart.*/
    w[i]=2.0*xl/((1.0-z*z)*pp*pp); /* Compute the weight */
    w[n+1-i]=w[i]; /* and its symmetric counterpart.*/
  }
}
void init_gauleg_weights(void)
{
  int n, i;
  n = OprogStatus.n_gauleg;
  gauleg_x = malloc(sizeof(double)*(n+1));
  gauleg_w = malloc(sizeof(double)*(n+1));
  qgaus_x = malloc(sizeof(double)*(n/2+1));
  qgaus_w = malloc(sizeof(double)*(n/2+1));
  if (n%2 != 0)
    {
      printf("ERROR: n_gauleg must be a multiple of 2\n");
      exit(-1);
    }
  gauleg(-1.0,1.0,gauleg_x,gauleg_w, n);
  qgaus_x[0] = 0.0;
  qgaus_w[0] = 0.0;
  //printf("x gaus = \n");
  for (i = n/2+1; i <= n; i++)
    {
      qgaus_x[i-n/2] = gauleg_x[i];
      qgaus_w[i-n/2] = gauleg_w[i]; 
      //printf("qgaus_x[%d]=%.15G", i-n/2, qgaus_x[i-n/2]);
    }
  //printf("\n");
  free(gauleg_x);
  free(gauleg_w);
}
double qgaus(double (*func)(double), double a, double b)
/* Returns the integral of the function func between a and b, by ten-point Gauss-Legendre integration:
   the function is evaluated exactly ten times at interior points in the range of integration.*/
{
  int j;
  double  xr,xm,dx,s;
  /* The abscissas and weights.  First value of each array not used.*/
#if 0 
  static double x[]={0.0,0.1488743389,0.4333953941, 0.6794095682,0.8650633666,0.9739065285};
  static double w[]={0.0,0.2955242247,0.2692667193, 0.2190863625,0.1494513491,0.0666713443};
  int n=10;
  static double w[] = {0.0,0.0775059479784209,0.0770398181639281,0.0761103619006262,0.0747231690579683,0.0728865823958041,0.0706116473912868,0.0679120458152339,0.0648040134566011,0.061306242492929,0.0574397690993916,0.0532278469839368,0.0486958076350722,0.0438709081856731,0.0387821679744716,0.0334601952825464,0.0279370069800163,0.0222458491941251,0.0164210583815049,0.0104982845214056,0.0045212770985331};
  static double x[] = {0.0,0.0387724175060508,0.116084070675255,0.192697580701371,0.268152185007254,0.341994090825758,0.413779204371605,0.483075801686179,0.549467125095128,0.61255388966798,0.67195668461418,0.727318255189927,0.778305651426519,0.824612230833312,0.865959503212259,0.902098806968874,0.932812808278677,0.957916819213792,0.977259949983774,0.990726238699457,0.998237709710559};
/*
  static double w[] = {0.0,0.152753387130726,0.149172986472604,0.142096109318382,0.131688638449177,0.118194531961518,0.101930119817233,0.0832767415766291,0.0626720483330506,0.040601429767854,0.0176140071391506};
  static double x[] = {0.0,0.0765265211334973,0.227785851141645,0.37370608871542,0.510867001950827,0.636053680726515,0.746331906460151,0.839116971822219,0.912234428251326,0.963971927277914,0.993128599185095};*/
  int n = 40;
#endif
  int n;
  n = OprogStatus.n_gauleg;
  xm=0.5*(b+a);
  xr=0.5*(b-a);
  s=0; 
  /* Will be twice the average value of the function, since the
     ten weights (five numbers above each used twice) sum to 2.*/
  for (j=1;j<=n/2;j++) 
    {
      dx=xr*qgaus_x[j];
      s += qgaus_w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
  return s *= xr; /* Scale the answer to the range of integration.*/
}
static double xsav;
static double (*nrfunc)(double, double);
double yy2(double x);
double yy2_pf(double x);
double KK=0.999999;
double quad3d_pf(double (*func)(double, double), double x1, double x2)
/*Returns the integral of a user-supplied function func over a three-dimensional region specified
  by the limits x1, x2, and by the user-supplied functions yy1, yy2, z1, and z2, as defined in
  (4.6.2). (The functions y1 and y2 are here called yy1 and yy2 to avoid conflict with the names
  of Bessel functions in some C libraries). Integration is performed by calling qgaus recursively.*/
{
  double f1_pf(double x);
  nrfunc=func;
  return qgaus(f1_pf,x1,x2);
}
double quad3d(double (*func)(double, double), double x1, double x2)
/*Returns the integral of a user-supplied function func over a three-dimensional region specified
  by the limits x1, x2, and by the user-supplied functions yy1, yy2, z1, and z2, as defined in
  (4.6.2). (The functions y1 and y2 are here called yy1 and yy2 to avoid conflict with the names
  of Bessel functions in some C libraries). Integration is performed by calling qgaus recursively.*/
{
  double f1(double x);
  nrfunc=func;
  return qgaus(f1,x1,x2);
}
double f1_pf(double x) /*This is H of eq. (4.6.5).*/
{
  double f2(double y);
  xsav=x;
  return qgaus(f2,0.0,yy2_pf(x));
}
double f1(double x) /*This is H of eq. (4.6.5).*/
{
  double f2(double y);
  xsav=x;
  return qgaus(f2,0.0,yy2(x));
}
double f2(double y)
{
  return (*nrfunc)(xsav,y);
}
inline double funcSQ_pf(int t, double x, double y)
{
  double n1, n2, n3;
  n1 = typesArr[t].n[0];
  n2 = typesArr[t].n[1];
  n3 = typesArr[t].n[2];
  return 8.0*pow(1.0-pow(fabs(x), n1)-pow(fabs(y), n2),1.0/n3);
}

inline double funcSQ(int t, double x, double y)
{
  double a, b, c, n1, n2, n3;
  if (growthSQ >= 0.0)
    {
      t = typeOfPart[growthSQ];
      a = axa[growthSQ];
      b = axb[growthSQ];
      c = axc[growthSQ];
    }
  else
    {
      a = typesArr[t].sax[0];
      b = typesArr[t].sax[1];
      c = typesArr[t].sax[2];
    }
  n1 = typesArr[t].n[0];
  n2 = typesArr[t].n[1];
  n3 = typesArr[t].n[2];
  return 8.0*c*pow(1.0-pow(fabs(x)/a, n1)-pow(fabs(y)/b, n2),1.0/n3);
}
int typeSQ;
double yy2_pf(double x)
{
  double n1, n2;
  int t;
  t = typeSQ;
  n1 = typesArr[t].n[0];
  n2 = typesArr[t].n[1];
  return pow(1.0-pow(x,n1),1.0/n2);
}
double yy2(double x)
{
  double a, b, n1, n2;
  int t;
  if (growthSQ >= 0)
    {
      t = typeOfPart[growthSQ];
      a = axa[growthSQ];
      b = axb[growthSQ];
    }
  else
    {
      t = typeSQ;
      a = typesArr[t].sax[0];
      b = typesArr[t].sax[1];
    }
  n1 = typesArr[t].n[0];
  n2 = typesArr[t].n[1];
  return b*pow(1.0-pow(x/a,n1),1.0/n2);
}
double funcSQwrap(double x, double y)
{
  return funcSQ(typeSQ, x, y);
}
double funcSQwrap_pf(double x, double y)
{
  return funcSQ_pf(typeSQ, x, y);
}
double calcVolSQ(int type)
{
  double vol;
  typeSQ = type;
  growthSQ = -1;
  vol = quad3d(funcSQwrap, 0.0, typesArr[type].sax[0]); /* a semiaxis (along x) */
  return vol;
}
double calcVolSQ_growth(int i)
{
  double vol;
  growthSQ = i;
  typeSQ = typeOfPart[i];
  vol = quad3d(funcSQwrap, 0.0, axa[i]); /* a semiaxis (along x) */
  return vol;
}
double calcVolSQ_pf(int type)
{
  double vol;
  typeSQ = type;
  vol = quad3d_pf(funcSQwrap_pf, 0.0, 1.0); /* a semiaxis (along x) */
  return vol;
}
double calcPhiSQ(void);
double calc_SQ_volprefact(int type)
{
  double tpo, pf;
  tpo = OprogStatus.targetPhi;
  OprogStatus.targetPhi = 1.0;
  pf = calcVolSQ_pf(type);
  OprogStatus.targetPhi = tpo;
  //printf("targetPhi=%f\n", OprogStatus.targetPhi);
  return pf;
}
double calcPhiSQ(void)
{
  int typei, i, nt;
  double N=0;
  if (OprogStatus.targetPhi > 0.0)
    {
      for (i=0; i < Oparams.parnum; i++)
	{
	  N += calcVolSQ_growth(i);
	  //printf("particle=%d volume=%.15G\n", i, calcVolSQ_growth(typei));
	}	  
    }
  else
    {
      for (nt = 0; nt < Oparams.ntypes; nt++)
	{
	  volSQ[nt] = calcVolSQ(nt);
	  printf("type=%d volume=%.15G\n", nt, volSQ[nt]);
	}
      for (i=0; i < Oparams.parnum; i++)
	{
	  typei = typeOfPart[i];
	  N += volSQ[typei];
	}
    }
#ifdef MD_LXYZ
  return N / (L[0]*L[1]*L[2]);
#else
  return N / (L*L*L);
#endif
}
#endif
