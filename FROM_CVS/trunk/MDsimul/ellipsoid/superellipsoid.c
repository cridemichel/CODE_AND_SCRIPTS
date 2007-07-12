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
#define MD_NEGPAIRS
#define MD_NO_STRICT_CHECK
#define MD_OPTDDIST
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
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **REtA, **REtB;
extern double cosEulAng[2][3], sinEulAng[2][3];
extern long long int itsFNL, timesFNL, timesSNL, itsSNL;
extern int do_check_negpairs;

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
struct LastBumpS *lastbump;
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
int calcdist_retcheck;
extern double rA[3], rB[3];
extern double gradplane_all[6][3], rBall[6][3];
extern int polinterr, polinterrRyck;
extern double max3(double a, double b, double c);
extern double scalProd(double *A, double *B);
extern int fdjac_disterr;
extern double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **powdirs;

int is_superellipse(int i)
{
  if (typesArr[typeOfPart [i]].n[0]==1.0 && typesArr[typeOfPart [i]].n[1]==1.0)
    return 0;
  else 
    return 1;
}
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
  inve2 = 2.0/e;	 
  xa2e = pow(x/a,inve2);
  yb2e = pow(y/a,inve2);
  A = pow(xa2e+yb2e,-1.0+e/n);
  inve2 = 2.0/e;
  fx[0] = (2.0*xa2e*A)/(n*x);
  fx[1] = (2.0*yb2e*A)/(n*y);
  fx[2] = (2.0*pow(z/c,2.0/n))/(n*z);
}

void calcfxx(double df[3][3], double x, double y, double z, int i)
{
  /* calcola fxx nel riferimento del corpo rigido */
  double a, b, c, n, e;
  double xa2e, yb2e, A, xa2eyb2e, emn, nm2, em2, inve2, eSqrnxy;
  a = typesArr[typeOfPart[i]].sax[0];
  b = typesArr[typeOfPart[i]].sax[1];
  c = typesArr[typeOfPart[i]].sax[2];
  e = typesArr[typeOfPart[i]].n[0];
  n = typesArr[typeOfPart[i]].n[1];
  inve2 = 2.0/e;
  xa2e = pow(x/a,inve2);
  yb2e = pow(y/a,inve2);
  eSqrnxy = e*Sqr(n)*x*y;
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
  df[1][1] = (A*(-2*em2*n*xa2e*yb2e- 2*e*nm2*Sqr(yb2e)))/(e*Sqr(n)*Sqr(y));
  df[1][2] = df[2][0] = df[2][1] = 0.0;
  df[2][2] = (-2*nm2*pow(z/c,2/n))/(Sqr(n)*Sqr(z)); 
}
extern void lab2body(int i, double x[], double xp[], double *rO, double **R);
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

void calcFxtFtSE(int i, double x[3], double **RM, double cosea[3], double sinea[3], double **X,
	       double D[3][3], double **R, 
	       double pos[3], double vel[3], double fxp[3], double fxxp[3][3],
	       double Fxt[3], double *Ft)
{
  double tRfxRdot[3][3], tRdotfxR[3][3], fxRdot[3]; 
  double DtX[3][3], dx[3];
  int k1, k2, k3;
  calc_Rdot(i, cosea, sinea, Rdot);
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {
      fxRdot[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fxRdot[k1] += Rdot[k1][k2]*fxp[k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRfxxRdot[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (fxR[k3][k1] == 0.0 || Rdot[k3][k2] == 0.0)
		continue;
	      tRfxxRdot[k1][k2] += fxx[k3][k1]*Rdot[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tRdotDfx[k1][k2] = 0.0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (Rdot[k3][k1] == 0.0 || fxR[k3][k2] == 0.0)
		continue;
	      tRdotfxR[k1][k2] += Rdot[k3][k1]*fxR[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      DtX[k1][k2] = tRdotfxR[k1][k2] + tRfxRdot[k1][k2];
  for (k1 = 0; k1 < 3; k1++)
    {
      Fxt[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	Fxt[k1] += DtX[k1][k2]*(x[k2]-pos[k2]) - X[k1][k2]*vel[k2]; 
      Fxt[k1] *= 2.0;
     } 
   *Ft = 0;
   for (k1 = 0; k1 < 3; k1++)
     dx[k1] = x[k1]-pos[k1];
   for (k1 = 0; k1 < 3; k1++)
     {
       for (k2 = 0; k2 < 3; k2++)
	 {
	   *Ft += -vel[k1]*X[k1][k2]*dx[k2]+dx[k1]*DtX[k1][k2]*dx[k2]-dx[k1]*X[k1][k2]*vel[k2];
	 }
     }
#endif

}
void fdjacSE(int n, double x[], double fvec[], double **df, 
	   void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  int na; 
#ifdef EDHE_FLEX
  int typei, typej;
#endif
  double  rA[3], rB[3], ti, vA[3], vB[3];
#ifndef MD_ASYM_ITENS
  double OmegaA[3][3], OmegaB[3][3];
#endif
  double DA[3][3], DB[3][3], fx[3], gx[3];
  double Fxt[3], Gxt[3], Ft=0.0, Gt=0.0;
#ifdef MD_ASYM_ITENS
  double phi, psi;
#endif
  int k1;
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
  symtop_evolve_orient(iB, ti, RB, REtB, cosEulAng[1], sinEulAng[1], &phi, &psi);
#else
  UpdateOrient(iB, ti, RB, OmegaB);
#endif
  na = (iB < Oparams.parnumA)?0:1;
  typej = typeOfPart[iB];

  lab2body(iA, &x[0], xpA, rA, RtA);
  lab2body(iB, &x[3], xpB, rB, RtB);  
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfx(gxp, xpB[0], xpB[1], xpB[2], iB);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2], iA);
  calcfxx(gxxp, xpB[0], xpB[1], xpB[2], iB);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fx(iB, gxp, gx, RtA);  
  body2lab_fxx(iA, fxxp, fxx, RtA);
  body2lab_fxx(iB, gxxp, gxx, RtB);

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
  calcFxtFtSE(iA, x, RM[iA], cosEulAng[0], sinEulAng[0], Xa, DA, RA, rA, vA, fxp, fxxp, Fxt, &Ft);
  calcFxtFtSE(iB, x, RM[iB], cosEulAng[1], sinEulAng[1], Xb, DB, RB, rB, vB, gxp, gxxp, Gxt, &Gt);
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
  calcfx(fxp, xpA[0], xpA[1], xpA[2], iA);
  calcfx(gxp, xpB[0], xpB[1], xpB[2], iB);
  calcfxx(fxxp, xpA[0], xpA[1], xpA[2] ,iA);
  calcfxx(gxxp, xpB[0], xpB[1], xpB[2], iB);
  /* ...and now we have to go back to laboratory reference system */
  body2lab_fx(iA, fxp, fx, RtA);
  body2lab_fx(iB, gxp, gx, RtA);  
  body2lab_fxx(iA, fxxp, fxx, RtA);
  body2lab_fxx(iB, gxxp, gxx, RtB);
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*fxx[k1][k2];
	  df[k1][k2+3] = 2.0*Sqr(x[6])*gxx[k1][k2];
	}
    }
#if 1
  if (OprogStatus.SDmethod == 2 || OprogStatus.SDmethod == 3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = x[k1+3] - x[k1];
#ifdef EDHE_FLEX
      for (kk=0; kk < 3; kk++)
	{
	  axi[kk] = typesArr[typeOfPart[iA]].sax[kk];
	  axj[kk] = typesArr[typeOfPart[iB]].sax[kk];
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
	    df[k1+5][k2] = 1 + 2.0*x[7]*fxx[k1][k2];
	  else 
	    df[k1+5][k2] = 2.0*x[7]*fxx[k1][k2];
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
void fdjacDistNeg5SE(int n, double x[], double fvec[], double **df, 
		   void (*vecfunc)(int, double [], double [], int, int, double []), 
		   int iA, int iB, double shift[3], double *fx, double *gx)
{
  double rDC[3], rD[3], A[3][3], b[3], c[3];
  int k1, k2;
#ifdef EDHE_FLEX
  int kk;
  double axi[3], axj[3];
#endif
  double xpA[3], xpB[3], fxxp[3][3], gxxp[3][3], fxp[3], gxp[3], fxx[3][3], gxx[3][3];
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
	  A[k1][k2] = Sqr(x[3])*(gxx[k1][k2] + x[4]*fxx[k1][k2]*gxx[k1][k2]);
	}
    }	
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = fxx[k1][k2] + A[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      rD[k1] = x[k1] + fx[k1]*x[4];
    } 
  lab2body(iB, rD, xpB, rB, RtB);  
  calcfx(gxp, xpB[0], xpB[1], xpB[2], iB);
  calcfxx(gxxp, xpB[0], xpB[1], xpB[2], iB);
  body2lab_fx(iB, gxp, gx, RtA);  
  body2lab_fxx(iB, gxxp, gxx, RtB);
  //printf("rC: %f %f %f rD: %f %f %f\n", x[0], x[1], x[2], rD[0], rD[1], rD[2]);
  //printf("fx: %f %f %f x[4]: %f\n", fx[0], fx[1], fx[2], x[4]);
#if 1
  if (OprogStatus.SDmethod==2 || OprogStatus.SDmethod==3)
    {
      for (k1 = 0; k1 < 3; k1++)
	rDC[k1] = rD[k1] - x[k1];
#ifdef EDHE_FLEX
      for (kk=0; kk < 3; kk++)
	{
	  axi[kk] = typesArr[typeOfPart[iA]].sax[kk];
	  axj[kk] = typesArr[typeOfPart[iB]].sax[kk];
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
	  b[k1] += Xb[k1][k2]*fx[k2];
	}
      b[k1] *= 2.0*Sqr(x[3]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      c[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  c[k1] += gx[k2]*fxx[k2][k1];
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
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
 //print_matrix(df, 5);
 //printf("fx: %f %f %f gx: %f %f %f\n", fx[0]/calc_norm(fx), fx[1]/calc_norm(fx), 
 //	fx[2]/calc_norm(fx), gx[0]/calc_norm(gx), gx[1]/calc_norm(gx), gx[2]/calc_norm(gx));
 //printf("F2BZdistNeg5 fvec (%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]);
#endif
}

void fdjacDistNegNeighPlaneSE(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double [], int), int iA)
{
  double fx[3];
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
#ifndef MD_GLOBALNRDNL
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[6])*gradplane[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 //fvec[4] = 0.5*fvec[4]-1.0;
 /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + gradplane[k1]*x[7];//[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif

}
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
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gradplane[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
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
  body2lab_fx(j, gxp, gx, RtA);  
 
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
  body2lab_fx(j, gxp, gx, RtA);  
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}

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
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gradplane[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
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

void funcs2beZeroedDistNegNeighPlane5SE(int n, double x[], double fvec[], int i)
{
  int k1, k2; 
  double fx[3], rD[3];

  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
	}
      rD[k1] = x[k1] + gradplane[k1]*x[4];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] - Sqr(x[3])*gradplane[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (rD[k1]-rB[k1])*gradplane[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif

}
/* 06/07/2007 TODO: steepest descent for super-ellipsoids! */
#endif
