#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
/* posizioni e velocit� di tutti gli atomi in un dischetto
 * (cio� anche quello massless!) */ 
#if 0
extern double **rallx, **rally, **rallz, **Fallx, **Fally, **Fallz,
  **rallx_old, **rally_old, **rallz_old, **vxold, **vyold, **vzold;
#endif
extern double **vxold, **vyold, **vzold, **rx_old, **ry_old, **rz_old;
#ifdef MD_RESPA
extern double  **rx_oldLong, **ry_oldLong, **rz_oldLong;
extern double *rxi[NA], *ryi[NA], *rzi[NA];
#endif
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern void  links(int Nm, COORD_TYPE rcut);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut);
extern void BuildNebrList(int Nm, COORD_TYPE rCut);
extern void CoMB(int a, int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
extern void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
extern void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm);
extern void checkNebrRebuild(void);
extern void LJForce(int Nm, COORD_TYPE rcut);
extern void resetCM(int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void kinet(int Nm, COORD_TYPE** velx, COORD_TYPE** vely, 
		  COORD_TYPE** velz, COORD_TYPE VOL1);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Wm, Wmxx, Wmyy, Wmzz, 
  Wmxy, Wmyz, Wmzx, Pmxx, Pmyy, Pmzz, Pmxy, Pmyz, Pmzx, T1mxy, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
  T1myz, T1mzx, T1mxx, T1myy, T1mzz;  
#ifdef MD_RESPA
extern COORD_TYPE WLong, WxxLong, WyyLong, WzzLong,
  WxyLong, WyzLong, WzxLong, WmLong, WmxxLong, WmyyLong, WmzzLong, 
  WmxyLong, WmyzLong, WmzxLong, WmyxLong, WmzyLong, WmxzLong,
  WCxxLong, WCyyLong, WCzzLong, WCxyLomg, WCyzLong, WCzxLong, WShort, VShort, VcShort;
double WmShort, WmxxShort, WmyyShort, WmzzShort, WmxyShort, WmyzShort, WmzxShort,
       WShort, VcShort, VShort, WxxShort, WyyShort, WzzShort, WxyShort, WyzShort, WzxShort;
#endif
extern double DrSq,  Mtot;
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;
extern double Volold, Volo1, Volo2, Volot;
#ifdef MD_RESPA
extern double VololdLong;
#endif
/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
#ifdef MD_RESPA
extern int **nebrTabLong, nebrNowLong, nebrTabLenLong, nebrTabMaxLong;
#endif
/* ================================= */

extern COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

extern COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
extern COORD_TYPE  *Rmx, *Rmy, *Rmz;
/* ========================================================================= */
extern void check_distances(char* str);
extern double *atcharge;
#ifdef MD_RAPACONSTR
#define N_MAX 1000
/* Routines for LU decomposition from Numerical Recipe online */
#define TINY 1E-20
#define MD_NCMAX 200
void ludcmp(double a[NA][NA], int* indx, double* d)
{
  /* A[i][j] = Aij 
   * A x = b  
   * per semplicit� nel seguito si assume che l'ordine della matrice � 3 */
  int i,imax=0,j,k;
  int n = NA;
  double big,dum,sum,temp; 
  double vv[NA]; /* vv stores the implicit scaling of each row.*/
  
  /*vv = vector(1,n);*/
  *d=1.0; /* No row interchanges yet. */
  for (i=0;i<n;i++) 
    { 
      /* Loop over rows to get the implicit scaling information.*/ 
      big=0.0; 
      for (j=0;j<n;j++) 
	if ((temp=fabs(a[i][j])) > big) big=temp; 
      if (big == 0.0)
	{
	  printf("ERROR: Singular matrix in routine ludcmp\n"); 
	  exit(-1);
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

void lubksb(double a[NA][NA], int* indx, double b[NA])
{ 
  int i,ii=0,ip,j; 
  double sum; 
  const int n=NA;
  for (i=0;i<n;i++) 
    { 
      /* When ii is set to a positive value, it will become the index of the  
       * rst nonvanishing element of b. Wenow do the forward substitution,
       * equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go. */
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i]; 
      if (ii) 
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
void SolveLineq (double a[NA][NA], double x[NA], int n) 
{
  int indx[NA];
  double dd;
  ludcmp(a, indx, &dd);
  lubksb(a, indx, x);
}
void BuildConstraintMatrix (void) 
{
  int i, m, j;
  for (i = 0; i < NA; i ++) 
    for (j = 0; j < NA; j++)
      cMat[i][j] = 0.;
  for (i = 0; i < NA; i ++) 
    {
      m = i - 1;
      if (m >= 0) 
	cMat[m][i] = 2;
      m = m + 1;
      if (m < NA-1) 
	cMat[m][i] = -2;
    }
  for (m = 0; m < NA-1; m ++) 
    {
      cDistSq[m] = Sqr(Oparams.d);
      cAtom1[m] = m; 
      cAtom2[m] = m + 1;
    }
}
double cvMat[NA][NA], cDistSq[NA], vVec[NA], curBondLenSq[NA]; 
int cMat[NA][NA]; 
int cAtom1[NA], cAtom2[NA];
void ComputeConstraints(int RefSys)
{
  /* RefSys=1 allora calcola le forze del ref system */
  double **FxI,**FyI, **FzI, **FxO, **FyO, **FzO; 
  double dv, w;
  int NB = NA-1;
  double rp1[3], rp2[3];
  int i, k, m, mDif, m1, m2, n, nn;
  if (RefSys)
    {
       FxI = Fx;
       FyI = Fy;
       FzI = Fz;
       FxO = FxC;
       FyO = FyC;
       FzO = FzC;
    }
  else
    {
      FxI = FxLong;
      FyI = FyLong;
      FzI = FzLong;
      FxO = FxCL;
      FyO = FyCL;
      FzO = FzCL;
    }
  for (n = 0; n < Oparams.parnum; n ++)
    {
      //nn = n  * NA;
      for (m = 0; m < NB; m ++) 
	{
	rp1[0] = rx[cAtom1[m]][n];
	rp1[1] = ry[cAtom1[m]][n];
	rp1[2] = rz[cAtom1[m]][n];
	rp2[0] = rx[cAtom2[m]][n];
	rp2[1] = ry[cAtom2[m]][n];
	rp2[2] = rz[cAtom2[m]][n];
	
	for (k = 0; k < 3; k ++) 
	  {
	    cVec[k][m] = rp1[k] - rp2[k];
	    cVec[k][m] = cVec[k][m] - L*rint (cVec[k][m]/L);
	  } 
      }
    for (m1 = 0; m1 < NB; m1 ++) 
      {
  	for (m2 = 0; m2 < NB; m2 ++) 
	  {
	    mDif = cMat[m1][cAtom1[m2]]/Oparams.m[m1] - cMat[m1][cAtom2[m2]]/Oparams.m[m2];
	    cvMat[m1][m2] = 0.;
	    if (mDif != 0) cvMat[m1][m2] = mDif * (cVec[0][m1] * cVec[0][m2] +
			       		      cVec[1][m1] * cVec[1][m2] + cVec[2][m1] * cVec[2][m2]);
	  } 
      }
    for (m = 0; m < NB; m ++) 
      {
	vVec[m] = 0.;
	if (RefSys)
	  dv = vx[cAtom1[m]][n] - vx[cAtom2[m]][n];
	else
	  dv = 0;
	vVec[m] = vVec[m] - (FxI[cAtom1[m]][n]/Oparams.m[cAtom1[m]] -
			     FxI[cAtom2[m]][n]/Oparams.m[cAtom2[m]]) * cVec[k][m] - Sqr (dv);
	if (RefSys)
	  dv = vy[cAtom1[m]][n] - vy[cAtom2[m]][n];
	else
	  dv = 0;
	vVec[m] = vVec[m] - (FyI[cAtom1[m]][n]/Oparams.m[cAtom1[m]] -
			     FyI[cAtom2[m]][n]/Oparams.m[cAtom2[m]]) * cVec[k][m] - Sqr (dv);
	if (RefSys)
	  dv = vz[cAtom1[m]][n] - vz[cAtom2[m]][n];
	else
	  dv = 0;
	vVec[m] = vVec[m] - (FzI[cAtom1[m]][n]/Oparams.m[cAtom1[m]] -
			     FzI[cAtom2[m]][n]/Oparams.m[cAtom2[m]]) * cVec[k][m] - Sqr (dv);
      }
    SolveLineq (cvMat, vVec, NB);
    for (m = 0; m < NB; m ++)
      {
	for (i = 0; i < NA; i ++) 
	  {
	    w = cMat[m * NA + i];
	    if (w != 0.) 
	      {
		/* qui calcolo le forze vincolari dovute alle forze nel reference system */
		FxO[i][n] = w * vVec[m] * cVec[0][m];
		FyO[i][n] = w * vVec[m] * cVec[1][m];
		FzO[i][n] = w * vVec[m] * cVec[2][m];
	      } 
	  }
      }
  }
}

void  AnlzConstraintDevs (void) 
{
  double dr1[3], sumL;
  int i, k, n, ni;
  sumL = 0.;
  for (n = 0; n < Oparams.parnum; n ++) 
    {
      for (i = 0; i < NA - 1; i ++) 
	{
	  dr1[0] = rx[i + 1][n] - rx[i][n];
  	  dr1[0] = dr1[0] - L*rint(dr1[0]/L);
	  dr1[1] = ry[i + 1][n] - ry[i][n];
  	  dr1[1] = dr1[1] - L*rint(dr1[1]/L);
	  dr1[2] = rz[i + 1][n] - rz[i][n];
  	  dr1[2] = dr1[2] - L*rint(dr1[2]/L);
	  curBondLenSq[i] = Sqr (dr1[0]) + Sqr (dr1[1]) + Sqr (dr1[2]);
	  sumL = sumL + curBondLenSq[i];
	}
    }
  constraintDevL = sqrt (sumL / (Oparams.parnum * (NA-1))) - Oparams.d;
  if (constraintDevL > 1E-6)
    {
      /* do Shake!! */
    }
}
#endif
#ifdef MD_RESPA_NPT
void shakePosRespa(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB, COORD_TYPE d, 
	   COORD_TYPE m[NA], int Nm)
{
  COORD_TYPE dSq;
  int      done;
  int      moving[NA], moved[NA];
  COORD_TYPE  tol2, pxab, pyab, pzab, pabSq, dt2, dtSq2;
  COORD_TYPE  rabSq, diffSq, rxab, ryab, rzab, rpab, gab;
  COORD_TYPE  dx, dy, dz, rma, rmb;
  COORD_TYPE  axia, ayia, azia;
  COORD_TYPE  vxi[NA], vyi[NA], vzi[NA], pxi[NA], pyi[NA], pzi[NA];
  int i, a, b, it;
  const COORD_TYPE rptol = 1.0E-6;
  double RCMx, RCMy, RCMz;

  tol2   = 2.0 * tol;
  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;
  dSq = Sqr(d); /* In general you must supply a vector of bond lengths */
 /* ===== LOOP OVER MOLECULES ===== */
  L = cbrt(Vol);
  for (i=0; i < Nm; i++)
    {
      /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
      for(a=0; a < NA; a++)
	{
	  //rxi[a][i] = rx[a][i];
  	  //ryi[a][i] = ry[a][i];
  	  //rzi[a][i] = rz[a][i];

	  pxi[a] = rx[a][i];
	  pyi[a] = ry[a][i];
	  pzi[a] = rz[a][i];
	  vxi[a] = px[a][i]/m[a];
	  vyi[a] = py[a][i]/m[a];
	  vzi[a] = pz[a][i]/m[a];
	  moving[a] = 0;
	  moved[a]  = 1;
	}
      it = 0;
      done = 0;
      /* START OF ITERATIVE LOOP */
      while ( (done == 0) && (it <= maxIt) )
	{
	  done = 1;

	  for(a=0; a < NB; a++)
	    {
	      b = a + 1;
	      /* la catena deve essere aperta */
	      if ( (moved[a]==1) || (moved[b] == 1)) 
		{
		  pxab = pxi[a] - pxi[b];
		  pyab = pyi[a] - pyi[b];
		  pzab = pzi[a] - pzi[b];
		  pxab = pxab - L * rint(invL * pxab);
		  pyab = pyab - L * rint(invL * pyab);
		  pzab = pzab - L * rint(invL * pzab);
		  pabSq = Sqr(pxab) + Sqr(pyab) + Sqr(pzab);
		  rabSq = dSq;
		  diffSq = rabSq - pabSq;
#if 0
		  printf("pabSq: %.15G diffSq: %.15G i=%d (%d-%d) rabSq:%.15G d:%.15G (%.15G, %.15G, %.15G)\n", pabSq, diffSq,i, a, b, sqrt(rabSq), d, vx[a][i], vy[a][i], vz[a][i]);
		  printf("a=%d b=%d rabSq: %f pabSq: %f\n", a, b, rabSq, pabSq);
		  printf("a=%d (%f,%f,%f) b=%d (%f,%f,%f)\n", a, pxi[a], pyi[a], pzi[a],
			 b, pxi[b], pyi[b], pzi[b]);
#endif
		  if ( fabs(diffSq) > ( rabSq * tol2 ) ) 
		    {
		      rxab = rxi[a][i] - rxi[b][i];
		      ryab = ryi[a][i] - ryi[b][i];
		      rzab = rzi[a][i] - rzi[b][i];
		      rxab = rxab - L * rint(invL * rxab);
		      ryab = ryab - L * rint(invL * ryab);
		      rzab = rzab - L * rint(invL * rzab);
		      rpab = rxab * pxab + ryab * pyab + rzab * pzab;
		      
		      if ( rpab < (rabSq * rptol) )
			{
			
			  mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
				"CONSTRAINT FAILURE IN MOVEA",
				NULL);
			  printf("i=%d between %d and %d pab=%f rab=%f\n", i, a, b, pabSq,
				 sqrt(Sqr(rxab)+Sqr(ryab)+Sqr(rzab)));
			  printf("rpab=%.15f\n", rpab);
			  printf("rab(%f,%f,%f) pab(%f,%f,%f)\n", rxab, ryab, rzab,
				 pxab, pyab, pzab);
			   printf("[%d-%d] v=(%f,%f,%f) a=(%f,%f,%f) pr=(%f,%f,%f) r=(%f,%f,%f)\n", a, i, 
				  vx[a][i], vy[a][i], vz[a][i],
				  axia, ayia, azia, pxi[a], pyi[a], pzi[a],
				  rxi[a][i], ryi[a][i], rzi[a][i]);

			  exit(-1);
			}
			   
		      rma = 1.0 / m[a];
		      rmb = 1.0 / m[b];
		      gab = diffSq / ( 2.0 * ( rma + rmb ) * rpab );
		      dx  = rxab * gab;
		      dy  = ryab * gab;
		      dz  = rzab * gab;
		      pxi[a] = pxi[a] + rma * dx;
		      pyi[a] = pyi[a] + rma * dy;
		      pzi[a] = pzi[a] + rma * dz;
		      pxi[b] = pxi[b] - rmb * dx;
		      pyi[b] = pyi[b] - rmb * dy;
		      pzi[b] = pzi[b] - rmb * dz;
			   
#if 1
		      dx = dx / dt;
		      dy = dy / dt;
		      dz = dz / dt;
		      vxi[a] = vxi[a] + rma * dx;
		      vyi[a] = vyi[a] + rma * dy;
		      vzi[a] = vzi[a] + rma * dz;
		      vxi[b] = vxi[b] - rmb * dx;
		      vyi[b] = vyi[b] - rmb * dy;
		      vzi[b] = vzi[b] - rmb * dz;
#endif			   
		      moving[a] = 1;
		      moving[b] = 1;
		      done = 0;
			   
		    }
		       
		}
		   
	    }
	       
	  for(a=0; a < NA; a++)
	    {
	      moved[a]  = moving[a];
	      moving[a] = 0;
	    }
	  it = it + 1;
	}
      /* END OF ITERATIVE LOOP */
      if  (done == 0) 
	{
	  sprintf(msgStrA, "MOLECULE N. %d", i);
	  mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
		"TOO MANY CONSTRAINT ITERATIONS IN MOVEA",
		msgStrA,
		NULL);
	  exit(-1);
	}
      /* STORE AWAY NEW VALUES */
      for(a=0; a < NA; a++)
	{
	  rx[a][i] = pxi[a];
	  ry[a][i] = pyi[a];
	  rz[a][i] = pzi[a];
#if 1
	  px[a][i] = vxi[a]*m[a];
	  py[a][i] = vyi[a]*m[a];
	  pz[a][i] = vzi[a]*m[a];
#endif
	  //printf("v(%f,%f,%f) F(%f,%f,%f)\n", vx[a][i], vy[a][i], vz[a][i],
	//	 Fx[a][i], Fy[a][i], Fz[a][i]);
	}
    }

}
void shakeVelRespaNPT(int Nm, COORD_TYPE dt, COORD_TYPE m[NA], int maxIt, int NB, 
		      COORD_TYPE d, COORD_TYPE tol, COORD_TYPE **p2sx, 
		      COORD_TYPE** p2sy, COORD_TYPE** p2sz )
{
  COORD_TYPE DRx, DRy, DRz;
  COORD_TYPE rxi[NA], ryi[NA], rzi[NA], vxi[NA], vyi[NA], vzi[NA];
  COORD_TYPE rxab, ryab, rzab, rvab, gab, dSq;
  COORD_TYPE vxab, vyab, vzab;
  COORD_TYPE dx, dy, dz, dt2, rma, rmb;
  int i, a, b, it;
  int done;
  int moving[NA], moved[NA];
  
  dSq = Sqr(d);

  WC = 0.0;
#if !defined(MOLPTENS) || defined(ATPTENS)  
  WCxy = 0.0; /* Constraints-virial off-diagonal terms of pressure tensor */
  WCyz = 0.0;
  WCzx = 0.0;
  WCxx = 0.0; /* Constraints-virial off-diagonal terms of pressure tensor */
  WCyy = 0.0;
  WCzz = 0.0;
#endif
  dt2 = dt * 0.5;
  for (i=0; i < Nm; i++)
    {
      /* VELOCITY VERLET ALGORITHM PART B */
      for(a=0; a < NA; a++)
	{
	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];

	  vxi[a] = p2sx[a][i]/Oparams.m[a];
	  vyi[a] = p2sy[a][i]/Oparams.m[a];
	  vzi[a] = p2sz[a][i]/Oparams.m[a];
	  
	  moving[a] = 0;
	  moved[a] = 1;
	}
      
      /* START OF ITERATIVE LOOP */
          it = 0;
	  done = 0;
	  while ((done == 0 ) && ( it <= maxIt ))
	    {
	      done = 1;
	      for(a=0; a < NB; a++)
		{
		  b = a + 1;
		  /* la catena � aperta */
		  /*if (b >= NA-1) b = 0;
		  */
		  if ( (moved[a] == 1) || (moved[b] == 1) ) 
		    {
		      rxab = rxi[a] - rxi[b];
		      ryab = ryi[a] - ryi[b];
		      rzab = rzi[a] - rzi[b];
		      
		      /* (DRx, DRy, DRz) is the change in the distance to
			 employ minimum image conventions, note that also 
			 relative velocity change, because the velocity is 
			 dependent upon position */
		      DRx = - L * rint(invL*rxab);
		      DRy = - L * rint(invL*ryab);
		      DRz = - L * rint(invL*rzab);
		      
		      rxab = rxab + DRx;
		      ryab = ryab + DRy;
		      rzab = rzab + DRz;
		      
		      vxab = vxi[a] - vxi[b];// + DRx * dlnV;
		      vyab = vyi[a] - vyi[b];// + DRy * dlnV;
		      vzab = vzi[a] - vzi[b];// + DRz * dlnV;
		      
		      rvab = rxab * vxab + ryab * vyab + rzab * vzab;
		      //printf("rvab:%f\n", rvab);
		      rma  = 1.0 / m[a];
		      rmb  = 1.0 / m[b];
		      gab  = (-rvab) / ( ( rma + rmb ) * dSq );
		      if ( fabs(gab) > tol )
			{
			  WC = WC + gab * dSq;

			  dx = rxab * gab;
			  dy = ryab * gab;
			  dz = rzab * gab;
			  
			  /* Non-diagonal terms of stress tensor 
			     NOTE: gab * (rxab, ryab, rzab) = F * (dt/2) */
#if !defined(MOLPTENS) || defined(ATPTENS)
			  WCxx += rxab * dx;
			  WCyy += ryab * dy;
			  WCzz += rzab * dz;
			  WCxy = WCxy + rxab * dy;
			  WCyz = WCyz + ryab * dz;
			  WCzx = WCzx + rzab * dx;
#endif 
			  vxi[a] = vxi[a] + rma * dx;
			  vyi[a] = vyi[a] + rma * dy;
			  vzi[a] = vzi[a] + rma * dz;
			  vxi[b] = vxi[b] - rmb * dx;
			  vyi[b] = vyi[b] - rmb * dy;
			  vzi[b] = vzi[b] - rmb * dz;
			  moving[a] = 1;
			  moving[b] = 1;
			  done = 0;
			}
		    }
		}
	      for(a=0; a < NA; a++)
		{
		  moved[a]  = moving[a];
		  moving[a] = 0;
		}
	      it = it + 1;
	    } 
	  
	  /* END OF ITERATIVE LOOP */
	  
	  if (done == 0)
	    {
	      sprintf(msgStrA, "MOLECULE N. %d", i);
	      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
		    "TOO MANY CONSTRAINT ITERATIONS IN MOVEB",
		    msgStrA,
		    NULL);
	      exit(-1);
	      
	    }
	  
	  for(a=0; a < NA; a++)
	    {
	      p2sx[a][i] = vxi[a]*Oparams.m[a];
	      p2sy[a][i] = vyi[a]*Oparams.m[a];
	      p2sz[a][i] = vzi[a]*Oparams.m[a];
	    }
    }
  /* END OF LOOP OVER MOLECULES */

  WC = WC / dt2 / 3.0;
#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = WCxy / dt2; /* WCxy, ... are not exactly virial terms and the 3.0
			is not present */
  WCyz = WCyz / dt2;
  WCzx = WCzx / dt2;

  WCxx = WCxx / dt2;
  WCyy = WCyy / dt2;
  WCzz = WCzz / dt2;
#endif
}
#endif
#if defined(MD_RESPA)
void calcPressTens(void)
{
  Wm = Wm + WmLong;
#ifdef MOLPTENS
  Wmxy = Wmxy + WmyzLong;
  Wmyz = Wmyz + WmyzLong;
  Wmzx = Wmzx + WmzxLong;
#endif
#if defined(MOLPTENS)  || (!defined(ATPTENS) && !defined(ATPRESS))
  Wmxx = Wmxx + WmxxLong;
  Wmyy = Wmyy + WmyyLong;
  Wmzz = Wmzz + WmzzLong;
#endif
#if defined(ATPTENS)
  Wxx = Wxx + WCxx + WxxLong + WCxxLong;
  Wyy = Wyy + WCyy + WyyLong + WCyyLong;
  Wzz = Wzz + WCzz + WzzLong + WCzzLong;
  Wxy = Wxy + WCxy + WCxyLong + WxyLong;
  Wzx = Wzx + WCzx + WCzxLong + WzxLong;
  Wyz = Wyz + WCyz + WCyzLong + WyzLong;
#endif
  //VcR = Vc;
  Vc = Vc + VcLong;
  V = V + VLong;
  W = W + WLong;
#if defined(ATPRESS)
  W = W + WCLong;
#endif
#ifdef MOLPTENS
  /* Calculate molecular pressure */
#ifdef MD_RESPA_NPT
  calcPtensMolRespa(Oparams.parnum);
#else
  calcPtensMol(Oparams.parnum, Vol1);
#endif
  press_m =  T1mxx + T1myy + T1mzz + Wm;
  press_m /= Vol * 3.0;
#endif  
  
#if !defined(MOLPTENS) || defined(ATPTENS)
  /* Calculate the other element (off-diagonal) of the pressure tensor */
#ifdef MD_RESPA_NPT
  calcPtensAtRespa(Oparams.parnum);
#else
  calcPressTens(Oparams.parnum, Vol1);
#endif
#if 0
  Patxy = T1xy + Wxy + WCxy;
  Patyz = T1yz + Wyz + WCyz;
  Patzx = T1zx + Wzx + WCzx;  
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
#endif
  press_at = (T1xx + T1yy + T1zz)/3.0 + WC + W;
  press_at /= Vol;
#elif defined(ATPRESS)
  /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
  press_at = calcT1diagAt(Nm, 0.0) + (W + WC) / Vol;
  /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
#endif

#ifdef MOLPTENS
   /* Calculate the other element (off-diagonal) of the molecular 
      pressure tensor */
#if 0
  Pmxy = T1mxy + Wmxy;
  Pmyz = T1myz + Wmyz;
  Pmzx = T1mzx + Wmzx;  
  Pmxy /= Vol;
  Pmyz /= Vol;
  Pmzx /= Vol;
#endif
  /* Calculate molecular pressure */
  press_m =  T1mxx + T1myy + T1mzz + Wm;
  press_m /= Vol * 3.0;
#endif

#ifdef MOLPTENS
  Pxy = Pmxy;
  Pyz = Pmyz;
  Pzx = Pmzx;
  press = press_m;
#else
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;
  press = press_at;
#endif
  /* ===================================================*/
}

#ifdef MD_RESPA_NPT
void p2v(void)
{ 
  int i, a;
  double Rxl, Ryl, Rzl;
  Vol1 = Pv * Sqr(s) / OprogStatus.W;
  s1   = Ps * Sqr(s) / OprogStatus.Q; 
  for (i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &Rxl, &Ryl, &Rzl);
      for (a=0; a < NA; a++)
	{
	  vx[a][i] = px[a][i]/Oparams.m[a] + (Vol1 / Vol / 3.0)*Rxl;
	  vy[a][i] = py[a][i]/Oparams.m[a] + (Vol1 / Vol / 3.0)*Ryl;
	  vz[a][i] = pz[a][i]/Oparams.m[a] + (Vol1 / Vol / 3.0)*Rzl;
	}
    }
}
void v2p(void)
{
  int i, a;
  double Rxl, Ryl, Rzl;
  for (i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &Rxl, &Ryl, &Rzl);
      for (a=0; a < NA; a++)
	{
	  px[a][i] = Oparams.m[a]*(vx[a][i] - (Vol1 / Vol / 3.0)*Rxl);
	  py[a][i] = Oparams.m[a]*(vy[a][i] - (Vol1 / Vol / 3.0)*Ryl);
	  pz[a][i] = Oparams.m[a]*(vz[a][i] - (Vol1 / Vol / 3.0)*Rzl);
	}
    }
  Ps = OprogStatus.Q * s1 / Sqr(s);
  Pv = OprogStatus.W * Vol1 / Sqr(s);
}

void updImpLong(double dt, double c)
{
  int i, a;
  double cdt;
  cdt = c*dt;
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	px[a][i] += cdt * FxLong[a][i];
      	py[a][i] += cdt * FyLong[a][i];
	pz[a][i] += cdt * FzLong[a][i];
      }
}
void updImpNose(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], cost[NA], dlns, cdt2;
  cdt = c*dt;
  cdt2 = cdt / 2.0;
  dlns = s*Ps / OprogStatus.Q ;
  for (a = 0; a < NA; a++)
    {
      expdt[a] = exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlns;*/
    }
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
  //      px[a][i] += Fx[a][i] * cdt2;
//	py[a][i] += Fy[a][i] * cdt2;
//	pz[a][i] += Fz[a][i] * cdt2;
	px[a][i] = px[a][i]*expdt[a];
      	py[a][i] = py[a][i]*expdt[a];
	pz[a][i] = pz[a][i]*expdt[a];
	px[a][i] += Fx[a][i] * cdt;
	py[a][i] += Fy[a][i] * cdt;
	pz[a][i] += Fz[a][i] * cdt;
      }
#ifndef MD_FENE
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000001, px, py, pz);
#endif
}
void updImpNoseAft(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], cost[NA], dlns, cdt2;
  cdt = c*dt;
  cdt2 = cdt / 2.0;
  dlns = s*Ps / OprogStatus.Q ;
  for (a = 0; a < NA; a++)
    {
      expdt[a] = exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlns;*/
    }
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
        px[a][i] += Fx[a][i] * cdt;
	py[a][i] += Fy[a][i] * cdt;
	pz[a][i] += Fz[a][i] * cdt;
	px[a][i] = px[a][i]*expdt[a];
      	py[a][i] = py[a][i]*expdt[a];
	pz[a][i] = pz[a][i]*expdt[a];
      }
#if 0
#ifndef MD_FENE
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000000001, px, py, pz);
#endif
#endif
}

void updImpLongNoseSym(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], cost[NA], dlns, cdt2;
  cdt = c*dt;
  cdt2 = cdt / 2.0;
  dlns = s*Ps / OprogStatus.Q ;
  for (a = 0; a < NA; a++)
    {
      expdt[a] = exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlns;*/
    }
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	px[a][i] += FxLong[a][i] * cdt2;
	py[a][i] += FyLong[a][i] * cdt2;
	pz[a][i] += FzLong[a][i] * cdt2;
	px[a][i] = px[a][i]*expdt[a];
      	py[a][i] = py[a][i]*expdt[a];
	pz[a][i] = pz[a][i]*expdt[a];
	px[a][i] += FxLong[a][i] * cdt2;
	py[a][i] += FyLong[a][i] * cdt2;
	pz[a][i] += FzLong[a][i] * cdt2;

      }
#ifndef MD_FENE
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000000001, px, py, pz);
#ifdef ATPRESS 
  WCLong = WC;
#endif
#ifdef ATPTENS  
  WCxxLong = WCxx;
  WCyyLong = WCyy;
  WCzzLong = WCzz;
  WCxyLong = WCxy;
  WCzxLong = WCzx;
  WCyzLong = WCyz;
#endif 
#endif
}


void updImpLongNoseAft(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], cost[NA], dlns, cdt2;
  cdt = c*dt;
  cdt2 = cdt / 2.0;
  dlns = s*Ps / OprogStatus.Q ;
  for (a = 0; a < NA; a++)
    {
      expdt[a] = exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlns;*/
    }
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	px[a][i] += FxLong[a][i] * cdt;
	py[a][i] += FyLong[a][i] * cdt;
	pz[a][i] += FzLong[a][i] * cdt;
	px[a][i] = px[a][i]*expdt[a];
      	py[a][i] = py[a][i]*expdt[a];
	pz[a][i] = pz[a][i]*expdt[a];
      }
#if 0
#ifndef MD_FENE
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000000001, px, py, pz);
#ifdef ATPRESS 
  WCLong = WC;
#endif
#ifdef ATPTENS  
  WCxxLong = WCxx;
  WCyyLong = WCyy;
  WCzzLong = WCzz;
  WCxyLong = WCxy;
  WCzxLong = WCzx;
  WCyzLong = WCyz;
#endif 
#endif
#endif
}

void updImpLongNose(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], cost[NA], dlns, cdt2;
  cdt = c*dt;
  cdt2 = cdt / 2.0;
  dlns = s*Ps / OprogStatus.Q ;
  for (a = 0; a < NA; a++)
    {
      expdt[a] = exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlns;*/
    }
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	px[a][i] = px[a][i]*expdt[a];
      	py[a][i] = py[a][i]*expdt[a];
	pz[a][i] = pz[a][i]*expdt[a];
	px[a][i] += FxLong[a][i] * cdt;
	py[a][i] += FyLong[a][i] * cdt;
	pz[a][i] += FzLong[a][i] * cdt;
      }
#ifndef MD_FENE
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.0000000000001, px, py, pz);
#ifdef ATPRESS 
  WCLong = WC;
#endif
#ifdef ATPTENS  
  WCxxLong = WCxx;
  WCyyLong = WCyy;
  WCzzLong = WCzz;
  WCxyLong = WCxy;
  WCzxLong = WCzx;
  WCyzLong = WCyz;
#endif 
#endif
}
void updImp(double dt, double c)
{
  int i, a;
  double cdt;
  cdt = c*dt;
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a=0; a < NA; a++)
	{
	  px[a][i] = cdt*Fx[a][i] + px[a][i];
	  py[a][i] = cdt*Fy[a][i] + py[a][i];
	  pz[a][i] = cdt*Fz[a][i] + pz[a][i];
	}
    }
}
void calcImp(int i, double *PCMx, double *PCMy, double *PCMz)
{
  int a;
  *PCMx = 0.0;
  *PCMy = 0.0;
  *PCMz = 0.0;
  for (a = 0; a < NA; a++)
    {
      *PCMx += px[a][i];
      *PCMy += py[a][i];
      *PCMz += pz[a][i];
    }
}

void updImpNoseAndAft(double dt, double c)
{
  int i, a;
  double cdt, expdtN[NA], expdt[NA], mM[NA], cost[NA], dlnVmM[NA];
  double dlnV, dlns;
  double DPx, DPy, DPz, PCMx, PCMy, PCMz, cdt2, pxo, pyo, pzo;
  cdt = c*dt;
  cdt2 = cdt/2.0;
  dlnV = Pv*Sqr(s)/OprogStatus.W/3.0/Vol;
  dlns = s*Ps / OprogStatus.Q;
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      dlnVmM[a] = mM[a] * dlnV;
      expdt[a] = exp(-dlnVmM[a] * cdt);
      expdtN[a]= exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlnVmM[a];*/
    }

   for (i=0; i < Oparams.parnum; i++)
    {
      PCMx = 0.0;
      PCMy = 0.0;
      PCMz = 0.0;
      for (a = 0; a < NA; a++)
	{
	  PCMx += px[a][i];
	  PCMy += py[a][i];
	  PCMz += pz[a][i];
	}
      for (a=0; a < NA; a++)
	{
	  pxo = px[a][i];
	  pyo = py[a][i];
	  pzo = pz[a][i];
	  DPx = PCMx-px[a][i];
	  DPy = PCMy-py[a][i];
	  DPz = PCMz-pz[a][i];
	  //px[a][i] = px[a][i]*expdt2[a];
	  //py[a][i] = py[a][i]*expdt2[a];
	  //pz[a][i] = pz[a][i]*expdt2[a];
	  px[a][i] += (-dlnVmM[a]*DPx)*cdt;
	  py[a][i] += (-dlnVmM[a]*DPy)*cdt;
	  pz[a][i] += (-dlnVmM[a]*DPz)*cdt;
	  px[a][i] = px[a][i]*expdt[a];
	  py[a][i] = py[a][i]*expdt[a];
	  pz[a][i] = pz[a][i]*expdt[a];
	  PCMx += px[a][i] - pxo;
	  PCMy += py[a][i] - pyo;
	  PCMz += pz[a][i] - pzo;
	} 
    }
#ifndef MD_RESPA_NOSELONG 
   for (i=0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	 //px[a][i] += Fx[a][i]*cdt2;
	 //py[a][i] += Fy[a][i]*cdt2;
	 //pz[a][i] += Fz[a][i]*cdt2;
	  px[a][i] = px[a][i]*expdtN[a];
	  py[a][i] = py[a][i]*expdtN[a];
	  pz[a][i] = pz[a][i]*expdtN[a];
	  px[a][i] += Fx[a][i]*cdt;
	  py[a][i] += Fy[a][i]*cdt;
	  pz[a][i] += Fz[a][i]*cdt;
	 }
    }
#endif
#if 0
#ifndef MD_FENE  
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 
		   0.0000000000001, px, py, pz);
#endif
#endif
}
#if 0
void updImpNoseAndSym(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], expdt2[NA], mM[NA], cost[NA], dlnVmM[NA];
  double dlnV, dlns;
  double DPx, DPy, DPz, PCMx, PCMy, PCMz, cdt2, pxo, pyo, pzo, cdt4, expdt4[NA];
  cdt = c*dt;
  cdt2 = cdt/2.0;
  cdt4 = cdt2 /2.0;
  dlnV = Pv*Sqr(s)/OprogStatus.W/3.0/Vol;
  dlns = s*Ps / OprogStatus.Q;
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      dlnVmM[a] = mM[a] * dlnV;
      expdt[a] = exp(-(dlnVmM[a] + dlns) * cdt);
      expdt2[a] = exp(-(dlnVmM[a] + dlns) * cdt2);
      expdt4[a] = exp(-(dlnVmM[a] + dlns) * cdt4);
      /*cost[a] = (expdt[a] - 1.0) / dlnVmM[a];*/
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      PCMx = 0.0;
      PCMy = 0.0;
      PCMz = 0.0;
      for (a = 0; a < NA; a++)
	{
	  PCMx += px[a][i];
	  PCMy += py[a][i];
	  PCMz += pz[a][i];
	}
      for (a=NA-1; a > 0; a--)
	{
	  pxo = px[a][i];
	  pyo = py[a][i];
	  pzo = pz[a][i];
	  DPx = PCMx-px[a][i];
	  DPy = PCMy-py[a][i];
	  DPz = PCMz-pz[a][i];
	  px[a][i] += (Fx[a][i]-dlnVmM[a]*DPx)*cdt4;
	  py[a][i] += (Fy[a][i]-dlnVmM[a]*DPy)*cdt4;
	  pz[a][i] += (Fz[a][i]-dlnVmM[a]*DPz)*cdt4;
	  px[a][i] = px[a][i]*expdt2[a];
	  py[a][i] = py[a][i]*expdt2[a];
	  pz[a][i] = pz[a][i]*expdt2[a];
	  px[a][i] += (Fx[a][i]-dlnVmM[a]*DPx)*cdt4;
	  py[a][i] += (Fy[a][i]-dlnVmM[a]*DPy)*cdt4;
	  pz[a][i] += (Fz[a][i]-dlnVmM[a]*DPz)*cdt4;
	  //PCMx += px[a][i] - pxo;
	  //PCMy += py[a][i] - pyo;
	  //PCMz += pz[a][i] - pzo;
	}
      a = 0;
      pxo = px[a][i];
      pyo = py[a][i];
      pzo = pz[a][i];
      DPx = PCMx-px[a][i];
      DPy = PCMy-py[a][i];
      DPz = PCMz-pz[a][i];
      px[a][i] += (Fx[a][i]-dlnVmM[a]*DPx)*cdt2;
      py[a][i] += (Fy[a][i]-dlnVmM[a]*DPy)*cdt2;
      pz[a][i] += (Fz[a][i]-dlnVmM[a]*DPz)*cdt2;
      px[a][i] = px[a][i]*expdt[a];
      py[a][i] = py[a][i]*expdt[a];
      pz[a][i] = pz[a][i]*expdt[a];
      px[a][i] += (Fx[a][i]-dlnVmM[a]*DPx)*cdt2;
      py[a][i] += (Fy[a][i]-dlnVmM[a]*DPy)*cdt2;
      pz[a][i] += (Fz[a][i]-dlnVmM[a]*DPz)*cdt2;
      //PCMx += px[a][i] - pxo;
      //PCMy += py[a][i] - pyo;
      //PCMz += pz[a][i] - pzo;
      for (a=1; a < NA; a++)
	{
	  pxo = px[a][i];
	  pyo = py[a][i];
	  pzo = pz[a][i];
	  DPx = PCMx-px[a][i];
	  DPy = PCMy-py[a][i];
	  DPz = PCMz-pz[a][i];
	  px[a][i] += (Fx[a][i]-dlnVmM[a]*DPx)*cdt4;
	  py[a][i] += (Fy[a][i]-dlnVmM[a]*DPy)*cdt4;
	  pz[a][i] += (Fz[a][i]-dlnVmM[a]*DPz)*cdt4;
	  px[a][i] = px[a][i]*expdt2[a];
	  py[a][i] = py[a][i]*expdt2[a];
	  pz[a][i] = pz[a][i]*expdt2[a];
	  px[a][i] += (Fx[a][i]-dlnVmM[a]*DPx)*cdt4;
	  py[a][i] += (Fy[a][i]-dlnVmM[a]*DPy)*cdt4;
	  pz[a][i] += (Fz[a][i]-dlnVmM[a]*DPz)*cdt4;
	  //PCMx += px[a][i] - pxo;
	  //PCMy += py[a][i] - pyo;
	  //PCMz += pz[a][i] - pzo;
	}
    }
}
#endif
void updImpNoseAnd(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], expdtN[NA], mM[NA], cost[NA], dlnVmM[NA];
  double dlnV, dlns;
  double DPx, DPy, DPz, PCMx, PCMy, PCMz, cdt2, pxo, pyo, pzo;
  cdt = c*dt;
  cdt2 = cdt/2.0;
  dlnV = Pv*Sqr(s)/OprogStatus.W/3.0/Vol;
  dlns = s*Ps / OprogStatus.Q;
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      dlnVmM[a] = mM[a] * dlnV;
      expdt[a] = exp(-dlnVmM[a] * cdt);
      expdtN[a]= exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlnVmM[a];*/
    }
#ifndef MD_RESPA_NOSELONG
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	  px[a][i] += Fx[a][i]*cdt;
	  py[a][i] += Fy[a][i]*cdt;
	  pz[a][i] += Fz[a][i]*cdt;
	  px[a][i] = px[a][i]*expdtN[a];
	  py[a][i] = py[a][i]*expdtN[a];
	  pz[a][i] = pz[a][i]*expdtN[a];
	  //px[a][i] += Fx[a][i]*cdt2;
	  //py[a][i] += Fy[a][i]*cdt2;
	  //pz[a][i] += Fz[a][i]*cdt2;
	 }
    }
#endif
  for (i=0; i < Oparams.parnum; i++)
    {
      PCMx = 0.0;
      PCMy = 0.0;
      PCMz = 0.0;
      for (a = 0; a < NA; a++)
	{
	  PCMx += px[a][i];
	  PCMy += py[a][i];
	  PCMz += pz[a][i];
	}
      for (a=NA-1; a >= 0; a--)
	{
	  pxo = px[a][i];
	  pyo = py[a][i];
	  pzo = pz[a][i];
	  DPx = PCMx-px[a][i];
	  DPy = PCMy-py[a][i];
	  DPz = PCMz-pz[a][i];
	  px[a][i] = px[a][i]*expdt[a];
	  py[a][i] = py[a][i]*expdt[a];
	  pz[a][i] = pz[a][i]*expdt[a];
	  px[a][i] += (-dlnVmM[a]*DPx)*cdt;
	  py[a][i] += (-dlnVmM[a]*DPy)*cdt;
	  pz[a][i] += (-dlnVmM[a]*DPz)*cdt;
	  //px[a][i] = px[a][i]*expdt2[a];
	  //py[a][i] = py[a][i]*expdt2[a];
	  //pz[a][i] = pz[a][i]*expdt2[a];
	  PCMx += px[a][i] - pxo;
	  PCMy += py[a][i] - pyo;
	  PCMz += pz[a][i] - pzo;
	} 
    }
#if 1
#ifndef MD_FENE  
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 
		   0.000000001, px, py, pz);
#endif
#endif
}

void updImpAnd(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], expdtN[NA], mM[NA], cost[NA], dlnVmM[NA];
  double dlnV, dlns;
  double DPx, DPy, DPz, PCMx, PCMy, PCMz, cdt2, pxo, pyo, pzo;
  cdt = c*dt;
  cdt2 = cdt/2.0;
  dlnV = Pv*Sqr(s)/OprogStatus.W/3.0/Vol;
  dlns = s*Ps / OprogStatus.Q;
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      dlnVmM[a] = mM[a] * dlnV;
      expdt[a] = exp(-dlnVmM[a] * cdt);
      expdtN[a]= exp(-dlns * cdt);
      /*cost[a] = (expdt[a] - 1.0) / dlnVmM[a];*/
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	  px[a][i] += Fx[a][i]*cdt;
	  py[a][i] += Fy[a][i]*cdt;
	  pz[a][i] += Fz[a][i]*cdt;
	  px[a][i] = px[a][i]*expdtN[a];
	  py[a][i] = py[a][i]*expdtN[a];
	  pz[a][i] = pz[a][i]*expdtN[a];
	  //px[a][i] += Fx[a][i]*cdt2;
	  //py[a][i] += Fy[a][i]*cdt2;
	  //pz[a][i] += Fz[a][i]*cdt2;
	 }
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      PCMx = 0.0;
      PCMy = 0.0;
      PCMz = 0.0;
      for (a = 0; a < NA; a++)
	{
	  PCMx += px[a][i];
	  PCMy += py[a][i];
	  PCMz += pz[a][i];
	}
      for (a=NA-1; a >= 0; a--)
	{
	  pxo = px[a][i];
	  pyo = py[a][i];
	  pzo = pz[a][i];
	  DPx = PCMx-px[a][i];
	  DPy = PCMy-py[a][i];
	  DPz = PCMz-pz[a][i];
	  px[a][i] = px[a][i]*expdt[a];
	  py[a][i] = py[a][i]*expdt[a];
	  pz[a][i] = pz[a][i]*expdt[a];
	  px[a][i] += (-dlnVmM[a]*DPx)*cdt;
	  py[a][i] += (-dlnVmM[a]*DPy)*cdt;
	  pz[a][i] += (-dlnVmM[a]*DPz)*cdt;
	  //px[a][i] = px[a][i]*expdt2[a];
	  //py[a][i] = py[a][i]*expdt2[a];
	  //pz[a][i] = pz[a][i]*expdt2[a];
	  PCMx += px[a][i] - pxo;
	  PCMy += py[a][i] - pyo;
	  PCMz += pz[a][i] - pzo;
	} 
    }
}

void updPositions(double dt, double c)
{
  int i, a;
  double cdt, *m;
  cdt = c*dt;
  m = Oparams.m;
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a=0; a < NA; a++)
	{
	  rx[a][i] = cdt*px[a][i]/m[a] + rx[a][i];
	  ry[a][i] = cdt*py[a][i]/m[a] + ry[a][i];
	  rz[a][i] = cdt*pz[a][i]/m[a] + rz[a][i];
	}
    }
} 

void updPositionsNPT(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], mM[NA], *m, expdt2[NA];
  double cost2, cdt2, cost2mM[NA], cdt4;
  double RCMx, RCMy, RCMz, rxo, ryo, rzo, DRCMx, DRCMy, DRCMz;
  cdt = c*dt;
  cdt2 = cdt / 2.0;
  cdt4 = cdt2 / 2.0;
  m = Oparams.m;
  cost2 = Pv*Sqr(s)/(3.0*Vol*OprogStatus.W);
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      cost2mM[a] = cost2 * mM[a];
      expdt[a] = exp(cost2mM[a]*cdt);
      expdt2[a] = exp(cost2mM[a]*cdt2);
    }
#if 0
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a=0; a < NA; a++)
	{
	  rx[a][i] += cdt2*px[a][i]/m[a];
	  ry[a][i] += cdt2*py[a][i]/m[a];
	  rz[a][i] += cdt2*pz[a][i]/m[a];
	}
    }
#endif
  for (i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &RCMx, &RCMy, &RCMz);
      for (a=NA-1; a >= 0; a--)
	{
	  rxo = rx[a][i];
	  ryo = ry[a][i];
	  rzo = rz[a][i];
	  DRCMx = RCMx - rx[a][i]*mM[a];
          DRCMy = RCMy - ry[a][i]*mM[a];
	  DRCMz = RCMz - rz[a][i]*mM[a];
	  rx[a][i] += cdt2*cost2*DRCMx;
	  ry[a][i] += cdt2*cost2*DRCMy;
	  rz[a][i] += cdt2*cost2*DRCMz;
       	  rx[a][i] = rx[a][i]*expdt2[a];
	  ry[a][i] = ry[a][i]*expdt2[a];
	  rz[a][i] = rz[a][i]*expdt2[a];
	  //rx[a][i] += cdt4*(px[a][i]/m[a] + cost2*DRCMx);
	  //ry[a][i] += cdt4*(py[a][i]/m[a] + cost2*DRCMy);
	  //rz[a][i] += cdt4*(pz[a][i]/m[a] + cost2*DRCMz);
	  RCMx = RCMx + (-rxo + rx[a][i])*mM[a]; 
	  RCMy = RCMy + (-ryo + ry[a][i])*mM[a];
	  RCMz = RCMz + (-rzo + rz[a][i])*mM[a];
	}
    }
#if 0	
  a = 0;
  rxo = rx[a][i];
  ryo = ry[a][i];
  rzo = rz[a][i];
  DRCMx = RCMx - rx[a][i]*mM[a];
  DRCMy = RCMy - ry[a][i]*mM[a];
  DRCMz = RCMz - rz[a][i]*mM[a];
  rx[a][i] += cdt2*cost2*DRCMx;
  ry[a][i] += cdt2*cost2*DRCMy;
  rz[a][i] += cdt2*cost2*DRCMz;
  rx[a][i] = rx[a][i]*expdt[a];
  ry[a][i] = ry[a][i]*expdt[a];
  rz[a][i] = rz[a][i]*expdt[a];
  rx[a][i] += cdt2*cost2*DRCMx;
  ry[a][i] += cdt2*cost2*DRCMy;
  rz[a][i] += cdt2*cost2*DRCMz;
  RCMx = RCMx + (-rxo + rx[a][i])*mM[a]; 
  RCMy = RCMy + (-ryo + ry[a][i])*mM[a];
  RCMz = RCMz + (-rzo + rz[a][i])*mM[a];
#endif
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a=0; a < NA; a++)
	{
	  rx[a][i] += cdt*px[a][i]/m[a];
	  ry[a][i] += cdt*py[a][i]/m[a];
	  rz[a][i] += cdt*pz[a][i]/m[a];
	}
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &RCMx, &RCMy, &RCMz);
      for (a=0; a < NA; a++)
	{
	  rxo = rx[a][i];
	  ryo = ry[a][i];
	  rzo = rz[a][i];
	  DRCMx = RCMx - rx[a][i]*mM[a];
	  DRCMy = RCMy - ry[a][i]*mM[a];
	  DRCMz = RCMz - rz[a][i]*mM[a];
	  //rx[a][i] += cdt4*(px[a][i]/m[a] + cost2*DRCMx);
	  //ry[a][i] += cdt4*(py[a][i]/m[a] + cost2*DRCMy);
	  //rz[a][i] += cdt4*(pz[a][i]/m[a] + cost2*DRCMz);
	  rx[a][i] = rx[a][i]*expdt2[a];
	  ry[a][i] = ry[a][i]*expdt2[a];
	  rz[a][i] = rz[a][i]*expdt2[a];
	  rx[a][i] += cdt2*cost2*DRCMx;
	  ry[a][i] += cdt2*cost2*DRCMy;
	  rz[a][i] += cdt2*cost2*DRCMz;
	  RCMx = RCMx + (-rxo + rx[a][i])*mM[a]; 
	  RCMy = RCMy + (-ryo + ry[a][i])*mM[a];
	  RCMz = RCMz + (-rzo + rz[a][i])*mM[a];
     	}
    }
#if 0
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a=0; a < NA; a++)
	{
	  rx[a][i] += cdt2*px[a][i]/m[a];
	  ry[a][i] += cdt2*py[a][i]/m[a];
	  rz[a][i] += cdt2*pz[a][i]/m[a];
	}
    }
#endif

} 
void updPs(double dt, double c)
{
  double cdt, cdt2;
  double dof, Nm;
  double DT, Kin;
  int i, a;
  cdt = c*dt;
  Nm = Oparams.parnum;
  cdt2 = cdt / 2.0;
  Kin = 0;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	Kin += (Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i]))/Oparams.m[a];  
      }
  Kin *= 0.5;
#ifdef MD_FENE
  dof = 3*NA*Oparams.parnum;
#else
  dof = (2*NA + 1)*Oparams.parnum;
#endif
  /*printf("temp= %f\n", 2*Kin/dof);*/
  DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
  /*printf("DT: %f T: %f Oparams.T: %f s:%f\n", DT, 2.0*Kin/dof, Oparams.T, s);*/
  Ps += DT * cdt;
  Ps = Ps / (1 + Ps*cdt*s/OprogStatus.Q);
  //Ps += DT * cdt2;
}
void updPsAft(double dt, double c)
{
  double cdt, cdt2;
  double dof, Nm;
  double DT, Kin;
  int i, a;
  cdt = c*dt;
  Nm = Oparams.parnum;
  cdt2 = cdt / 2.0;
  Kin = 0;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	Kin += (Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i]))/Oparams.m[a];  
      }
  Kin *= 0.5;
#ifdef MD_FENE
  dof = 3*NA*Oparams.parnum;
#else
  dof = (2*NA + 1)*Oparams.parnum;
#endif
  /*printf("temp= %f\n", 2*Kin/dof);*/
  DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
  /*printf("DT: %f T: %f Oparams.T: %f s:%f\n", DT, 2.0*Kin/dof, Oparams.T, s);*/
  Ps = Ps / (1 + Ps*cdt*s/OprogStatus.Q);
  Ps += DT * cdt;
  //Ps += DT * cdt2;
}
void upds(double dt, double c)
{
  double cdt;
  cdt = c*dt;
  s = s / (1 - s*Ps*cdt/OprogStatus.Q);
}
void updLs(double dt, double c)
{
  double cdt, cdt2;
  double dof, Nm;
  double DT, Kin;
  int i, a;
  cdt = c*dt;
  Nm = Oparams.parnum;
  cdt2 = cdt / 2.0;
  Kin = 0;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	Kin += (Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i]))/Oparams.m[a];  
      }
  Kin *= 0.5;
#ifdef MD_FENE
  dof = 3*NA*Oparams.parnum;
#else
  dof = (2*NA + 1)*Oparams.parnum;
#endif
  /*printf("temp= %f\n", 2*Kin/dof);*/
  DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
  /*printf("DT: %f T: %f Oparams.T: %f s:%f\n", DT, 2.0*Kin/dof, Oparams.T, s);*/
  Ps = Ps / (1 + Ps*cdt*s/OprogStatus.Q);
  Ps += DT * cdt;
  s = s / (1 - s*Ps*cdt/OprogStatus.Q);
}
void updLsAft(double dt, double c)
{
  double cdt, cdt2;
  double dof, Nm;
  double DT, Kin;
  int i, a;
  cdt = c*dt;
  Nm = Oparams.parnum;
  cdt2 = cdt / 2.0;
  s = s / (1 - s*Ps*cdt/OprogStatus.Q);
  Kin = 0;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	Kin += (Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i]))/Oparams.m[a];  
      }
  Kin *= 0.5;
#ifdef MD_FENE
  dof = 3*NA*Oparams.parnum;
#else
  dof = (2*NA + 1)*Oparams.parnum;
#endif
  /*printf("temp= %f\n", 2*Kin/dof);*/
  DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
  /*printf("DT: %f T: %f Oparams.T: %f s:%f\n", DT, 2.0*Kin/dof, Oparams.T, s);*/
  Ps += DT * cdt;
  Ps = Ps / (1 + Ps*cdt*s/OprogStatus.Q);
}
/*========================== >>> calcPtensMol <<< ===========================*/
void calcPtensMolRespa(int Nm)
{
  int i, a;
  COORD_TYPE Px, Py, Pz, *m;
  
  T1mxy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1myz = 0.0;
  T1mzx = 0.0;
  T1mxx = 0.0; /* Kinetic diagonal terms of pressure tensor */
  T1myy = 0.0;
  T1mzz = 0.0;
  m = Oparams.m;
  for(i=0; i < Nm; i++)
    {
      Px = Py = Pz = 0.0;  /* Velocity of center of mass of molecule i */
      for(a=0; a < NA; a++)
	{
	  Px += px[a][i];
	  Py += py[a][i];
	  Pz += pz[a][i];
	}
      /* Kinetic component of pressure tensor (all terms) */
      T1mxy += Px * Py / Mtot; 
      T1myz += Py * Pz / Mtot;
      T1mzx += Pz * Px / Mtot;
      T1mxx += Px * Px / Mtot;
      T1myy += Py * Py / Mtot;
      T1mzz += Pz * Pz / Mtot;
    }
  /* Calculate the other element (off-diagonal) of the molecular 
     pressure tensor */
  Pmxy = T1mxy + Wmxy;
  Pmyz = T1myz + Wmyz;
  Pmzx = T1mzx + Wmzx;  
  Pmxx = T1mxx + Wmxx;
  Pmyy = T1myy + Wmyy;
  Pmzz = T1mzz + Wmzz;

  Pmxx /= Vol;
  Pmyy /= Vol;
  Pmzz /= Vol;
  Pmxy /= Vol;
  Pmyz /= Vol;
  Pmzx /= Vol;
}
/*============================ >>> calcPtensAt <<< =========================*/
void calcPtensAtRespa(int Nm, COORD_TYPE VOL1)
{
  /* Calculate all components of atomic pressure tensor */
  int i, a;
  COORD_TYPE *m;
  
  m = Oparams.m;
  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;

  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  /* Kinetic component of pressure tensor (all terms) */
	  T1xy += px[a][i] * py[a][i] / m[a]; 
	  T1yz += py[a][i] * pz[a][i] / m[a];
	  T1zx += pz[a][i] * px[a][i] / m[a] ;
	  T1xx += px[a][i] * px[a][i] / m[a];
	  T1yy += py[a][i] * py[a][i] / m[a];
	  T1zz += pz[a][i] * pz[a][i] / m[a];
	}
    }
  /* Calculate the other element (off-diagonal) of the molecular 
     pressure tensor */
  Patxy = T1xy + Wxy + WCxy;
  Patyz = T1yz + Wyz + WCyz;
  Patzx = T1zx + Wzx + WCzx;  
  Patxx = T1xx + Wxx + WCxx;
  Patyy = T1yy + Wyy + WCyy;
  Patzz = T1zz + Wzz + WCzz;
  Patxx /= Vol;
  Patyy /= Vol;
  Patzz /= Vol;
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
}

/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagMolRespa(int Nm)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of molecular pressure tensor */
  int i, a;
  COORD_TYPE kin;
  COORD_TYPE Px, Py, Pz;

  kin = 0.0;
  for(i=0; i < Nm; i++)
    {
      Px = 0.0;
      Py = 0.0;
      Pz = 0.0;
      for (a=0; a < NA; a++)
	{
	  Px += px[a][i];
	  Py += py[a][i];
	  Pz += pz[a][i];
	}
      kin +=  Sqr(Px) + Sqr(Py) + Sqr(Pz);
    }
  kin /= 3.0 * Vol * Mtot;
  return kin;

}
/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagAtRespa(int Nm)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of atomic pressure tensor */
  int i, a;
  COORD_TYPE kin, kina, *m;
  m = Oparams.m;
  kin = 0.0;
  for(i=0; i < Nm; i++)
    {
      kina = 0.0;
      for(a=0; a < NA; a++)
	{
	  kina +=  Sqr(px[a][i]) + Sqr(py[a][i]) + Sqr(pz[a][i]);
	}
      kina /= m[a];
      kin += kina;
    }
  kin /= 3.0 * Vol;
  return kin;
}
void updVol(double dt, double c)
{
  double cdt = c*dt;
  Vol += cdt*Sqr(s)*Pv/OprogStatus.W;
}
void updPvLong(double dt, double c)
{
  double press, cdt, cdt2, DP, Nm;
  Nm = Oparams.parnum;
  cdt = c * dt;
  cdt2 = c * dt / 2.0;
#ifdef MOLPTENS
  DP = WmLong  / 3.0 / Vol; /* press(t+dt) */
#else
  DP = ( WLong + WCLong) / Vol; /* press(t+dt) */
#endif
#if 0
    {
    FILE *f;
    f = fopen("press.dat", "a");
    fprintf(f, "%d %f\n", Oparams.curStep, press);
    fclose(f);
    }
#endif
  /*printf(">>>>>>>>> press: %f DP: %f WmShort: %f WmLong: %f\n", press, DP, WmShort, WmLong);*/
  Pv += DP  * cdt;
}
void updPv(double dt, double c)
{
  double press, cdt, cdt2, DP, Nm;
  Nm = Oparams.parnum;
  cdt = c * dt;
  cdt2 = c * dt / 2.0;
#ifdef MOLPTENS
  press = calcT1diagMolRespa(Nm) + WmShort / 3.0 / Vol; /* press(t+dt) */
#else
  press = calcT1diagAtRespa(Nm) + (WShort + WCShort ) / Vol; /* press(t+dt) */
#endif
  DP = press - Oparams.P;
#if 0
    {
    FILE *f;
    f = fopen("press.dat", "a");
    fprintf(f, "%d %f\n", Oparams.curStep, press);
    fclose(f);
    }
#endif
  /*printf(">>>>>>>>> press: %f DP: %f WmShort: %f WmLong: %f\n", press, DP, WmShort, WmLong);*/
  Pv += DP  * cdt;
  Pv *= exp(-cdt*Ps*s/OprogStatus.Q);
  ///Pv += DP  * cdt2;
}
void updPvAft(double dt, double c)
{
  double press, cdt, cdt2, DP, Nm;
  Nm = Oparams.parnum;
  cdt = c * dt;
  cdt2 = c * dt / 2.0;
#ifdef MOLPTENS
  press = calcT1diagMolRespa(Nm) + (WmShort) / 3.0 / Vol; /* press(t+dt) */
#else
  press = calcT1diagAtRespa(Nm) + (WShort + WCShort) / Vol; /* press(t+dt) */
#endif
  DP = press - Oparams.P;
#if 0
    {
    FILE *f;
    f = fopen("press.dat", "a");
    fprintf(f, "%d %f\n", Oparams.curStep, press);
    fclose(f);
    }
#endif
  /*printf(">>>>>>>>> press: %f DP: %f WmShort: %f WmLong: %f\n", press, DP, WmShort, WmLong);*/
  Pv *= exp(-cdt*Ps*s/OprogStatus.Q);
  Pv += DP  * cdt;
  ///Pv += DP  * cdt2;
}
/* =========================== >>> kinet <<< ============================== */
void kinetRespaNPT(int Nm, COORD_TYPE** px, COORD_TYPE** py, COORD_TYPE** pz)
{
  int i, a;
  K = 0.0;
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  K = K + (Sqr(px[a][i]) + Sqr(py[a][i]) + Sqr(pz[a][i]))/Oparams.m[a];
	}
    }
  K *= 0.5;
}
void movelongRespaNPTBef(double dt)
{
#ifdef MD_RESPA_NOSELONG
  if (OprogStatus.Nose == 0)
    updImpLong(dt, 0.5);
  else
    {
      /*printf("1) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);*/
      updLs(dt, 0.5);
      updImpLongNose(dt, 0.5);
      if (OprogStatus.Nose==1)
	updPvLong(dt, 0.5);
      /*printf("7) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);*/
    }
#else
   updImpLong(dt, 0.5);
   if (OprogStatus.Nose==1)
     updPvLong(dt, 0.5);
#if 1
#ifndef MD_FENE
  if (OprogStatus.rcutInner != Oparams.rcut)
    shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.0000000000001, px, py, pz);
#ifdef ATPRESS 
  WCLong = WC;
#endif
#ifdef ATPTENS  
  WCxxLong = WCxx;
  WCyyLong = WCyy;
  WCzzLong = WCzz;
  WCxyLong = WCxy;
  WCzxLong = WCzx;
  WCyzLong = WCyz;
#endif 
#endif
#endif
#endif
}

void movelongRespaNPTBefAlt(double dt)
{
#ifdef MD_RESPA_NOSELONG
  double cdt, cdt2;
  double dof, Nm;
  double DT, Kin;
  int i, a;
  double c;
  if (OprogStatus.Nose == 0)
    updImpLong(dt, 0.5);
  else
    {
      cdt = 0.5*dt;
      LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
      Nm = Oparams.parnum;
      cdt2 = cdt / 2.0;
      ///s = s / (1 - s*Ps*cdt/OprogStatus.Q);
      Kin = 0;
      for (i = 0; i < Oparams.parnum; i++)
    	for (a = 0; a < NA; a++)
	  {
	    Kin += (Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i]))/Oparams.m[a];  
	  }
      Kin *= 0.5;
#ifdef MD_FENE
      dof = 3*NA*Oparams.parnum;
#else
      dof = (2*NA + 1)*Oparams.parnum;
#endif
      /*printf("temp= %f\n", 2*Kin/dof);*/
      DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
      /*printf("DT: %f T: %f Oparams.T: %f s:%f\n", DT, 2.0*Kin/dof, Oparams.T, s);*/
      Ps += DT * cdt;
      Ps = Ps / (1 + Ps*cdt*s/OprogStatus.Q);
      //Ps += DT * cdt2;
      updImpLongNoseSym(dt, 0.5);
      if (OprogStatus.Nose==1)
	updPvLong(dt, 0.5);
      s = s / (1 - s*Ps*cdt/OprogStatus.Q);
    }
#else
  updImpLong(dt, 0.5);
  if (OprogStatus.Nose == 1)
    updPvLong(dt, 0.5);
#endif
}

void movelongRespaNPTAftAlt(double dt)
{
#ifdef MD_RESPA_NOSELONG
  double cdt, cdt2;
  double dof, Nm;
  double DT, Kin;
  int i, a;
  double c = 0.5;
  if (OprogStatus.Nose == 0)
    {
      LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
      updImpLong(dt, 0.5);
    }
  else
    {
      cdt = c*dt;
      LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
      Nm = Oparams.parnum;
      cdt2 = cdt / 2.0;
      s = s / (1 - s*Ps*cdt/OprogStatus.Q);
      if (OprogStatus.Nose==1)
	updPvLong(dt, 0.5);
      updImpLongNoseSym(dt, 0.5);
      Kin = 0;
      for (i = 0; i < Oparams.parnum; i++)
    	for (a = 0; a < NA; a++)
	  {
	    Kin += (Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i]))/Oparams.m[a];  
	  }
      Kin *= 0.5;
#ifdef MD_FENE
      dof = 3*NA*Oparams.parnum;
#else
      dof = (2*NA + 1)*Oparams.parnum;
#endif
      /*printf("temp= %f\n", 2*Kin/dof);*/
      DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
      /*printf("DT: %f T: %f Oparams.T: %f s:%f\n", DT, 2.0*Kin/dof, Oparams.T, s);*/
      Ps = Ps / (1 + Ps*cdt*s/OprogStatus.Q);
      Ps += DT * cdt;
      //Ps += DT * cdt2;
      //s = s / (1 - s*Ps*cdt2/OprogStatus.Q);
    }
#else
  LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
  if (OprogStatus.Nose==1)
    updPvLong(dt, 0.5);
  updImpLong(dt, 0.5);
#endif   
}

void movelongRespaNPTAft(double dt)
{
#ifdef MD_RESPA_NOSELONG
  if (OprogStatus.Nose == 0)
    {
      LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
      updImpLong(dt, 0.5);
    }
  else
    {
      LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
      if (OprogStatus.Nose==1)
	updPvLong(dt, 0.5);
      updImpLongNoseAft(dt, 0.5);
      updLsAft(dt, 0.5);
      /*printf("A1) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);*/
    }
#else
  LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
  if (OprogStatus.Nose==1)
    updPvLong(dt, 0.5);
  updImpLong(dt, 0.5);
#if 1
#ifndef MD_FENE
  if (OprogStatus.rcutInner != Oparams.rcut)
    shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000001, px, py, pz);
#ifdef ATPRESS 
  WCLong = WC;
#endif
#ifdef ATPTENS  
  WCxxLong = WCxx;
  WCyyLong = WCyy;
  WCzzLong = WCzz;
  WCxyLong = WCxy;
  WCzxLong = WCzx;
  WCyzLong = WCyz;
#endif 
#endif
#endif
#endif
}
/* ========================== >>> movea <<< =============================== */
void moveaRespa(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB, COORD_TYPE d, 
	   COORD_TYPE m[NA], int Nm)
{   
  int i, a;
  
  if ( ( NB != NA ) && ( NB != NA-1 ) ) 
    {
      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,  
	    "NB IN ERROR", 
	    NULL);
      exit(-1);
    }
  L = cbrt(Vol);
  invL = 1.0 / L;
#if 0
  printf("d=%.15f NB=%d NA=%d m[0]=%f m[1]=%f m[2]=%f\n",d,NB,NA,m[0], m[1], m[2]);
  printf("tol=%.20f\n", tol);
  printf("L=%f Vol=%.f dt=%f\n", L, Vol, dt);
#endif
#ifndef MD_FENE
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	rxi[a][i] = rx[a][i];
	ryi[a][i] = ry[a][i];
	rzi[a][i] = rz[a][i];
      }	 
#endif
#ifdef MD_RESPA_NOSELONG
  if (OprogStatus.Nose == 1)
    {
      updImpNoseAnd(dt, 0.5);
      updPv(dt, 0.5);
      updVol(dt, 0.5);
      updPositionsNPT(dt, 1.0);
      updVol(dt, 0.5);
    }
#else
  if (OprogStatus.Nose == 1)
    {
#ifdef MD_FENE
      updPv(dt, 0.5);
      updImpNoseAnd(dt, 0.5);
      updPs(dt, 0.5);
      updVol(dt, 0.5);
      upds(dt, 0.5);
      updPositionsNPT(dt, 1.0);
      upds(dt, 0.5);
      updVol(dt, 0.5);
#else
      updPv(dt, 0.5);
      updPs(dt, 0.5);
      updImpNoseAnd(dt, 0.5);
      updVol(dt, 0.5);
      upds(dt, 0.5);
      updPositionsNPT(dt, 1.0);
      upds(dt, 0.5);
      updVol(dt, 0.5);
#endif
    }
  else if (OprogStatus.Nose == 2)
    {
#ifdef MD_FENE
      updImpNose(dt, 0.5);
      updPs(dt, 0.5);
      upds(dt, 0.5);
      updPositions(dt, 1.0);
      upds(dt, 0.5);
#else
      updPs(dt, 0.5);
      updImpNose(dt, 0.5);
      upds(dt, 0.5);
      updPositions(dt, 1.0);
      upds(dt, 0.5);
#endif
    }
#endif
  else
    {
      updImp(dt, 0.5);
#ifndef MD_FENE  
      shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 
    		       0.000000001, px, py, pz);
#endif
      updPositions(dt, 1.0);
    }

#if !defined(MD_FENE)
  shakePosRespa(Oparams.steplength/OprogStatus.nrespa, 0.000000001, 150, NA-1, Oparams.d,
		Oparams.m, Oparams.parnum);
#endif

#if 0
  if (OprogStatus.Nose == 1)
    updVol(dt, 0.5);
#endif
}
#endif
#ifdef MD_RESPA_NPT
/* ========================= >>> moveb <<< =========================== */
void movebRespa(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	   COORD_TYPE m[NA], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE dx, dy, dz, dt2;
  /* ******************************************************************* */
  dt2 = dt*0.5;
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;
#ifdef MD_RESPA_NOSELONG
  if (OprogStatus.Nose == 1)
    {
      updPvAft(dt, 0.5);
      updImpNoseAndAft(dt, 0.5);
    }
#else
  if (OprogStatus.Nose == 1)
    {
#ifdef MD_FENE
      updPsAft(dt, 0.5);
      updImpNoseAndAft(dt, 0.5);
      updPvAft(dt, 0.5);
#else
      updImpNoseAndAft(dt, 0.5);
      updPsAft(dt, 0.5);
      updPvAft(dt, 0.5);
#endif
    }
  else if (OprogStatus.Nose == 2)
    {
#ifdef MD_FENE
      updPsAft(dt, 0.5);
      updImpNoseAft(dt, 0.5);
#else
      updImpNoseAft(dt, 0.5);
      updPsAft(dt, 0.5);
#endif
    }
#endif
  else
    {
      updImp(dt, 0.5);
#if 0
#ifndef MD_FENE  
      shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 
	  	       0.000000000001, px, py, pz);
#endif
#endif
    } 
#if 0
#ifndef MD_FENE  
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 
		   0.000000000001, px, py, pz);
#endif
#endif
}
#endif
#endif
#ifdef MD_RESPA
/* ============================ >>> move<<< =================================*/
void move(void)
{
  /* DESCRIPTION:
     Move the particles by one step */
  const int n = OprogStatus.nrespa;
  double distance;
  /* distanza fra i 3 atomi le cui coordinate evolvono nel tempo
   * Notare che i 3 atomi formano un triangolo equiliatero quindi
   * le 3 distanze sono uguali */
  double dt = Oparams.steplength;
  double VcR;
  int a, i, kk;
  distance = Oparams.d;
  /* calc predicted coords*/
  Mtot = 0;
  for (a = 0; a < NA; a++)
    Mtot += Oparams.m[a];
#ifdef MD_RESPA_NPT
  movelongRespaNPTBef(Oparams.steplength);
#else
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	vx[a][i] += 0.5 * dt * FxLong[a][i]/Oparams.m[a];
      	vy[a][i] += 0.5 * dt * FyLong[a][i]/Oparams.m[a];
	vz[a][i] += 0.5 * dt * FzLong[a][i]/Oparams.m[a];
      }
#ifndef MD_FENE
  shakeVel(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000000001, vx, vy, vz);
#endif
#endif
  for (kk=0; kk < n; kk++)
    {

#ifdef MD_RESPA_NPT
      moveaRespa(Oparams.steplength/n, 0.000000000001, 150, NA-1, distance, Oparams.m, 
	    Oparams.parnum);        
#else
      movea(Oparams.steplength/n, 0.000000000001, 150, NA-1, distance, Oparams.m, 
	    Oparams.parnum);        
#endif
      /* buildAtomsPositions();*/
      if (nebrNow)
	{
	  nebrNow = 0;
	  /* build up linked list on predicted 
	     coordinates (only father do it)*/
	  if (OprogStatus.noLinkedList)
	    {
	      BuildNebrListNoLinked(Oparams.parnum, OprogStatus.rcutInner);
	    }
	  else
	    {
	      links(Oparams.parnum, OprogStatus.rcutInner);
	      /* Build up neighbour list */  
	      BuildNebrList(Oparams.parnum, OprogStatus.rcutInner);
	    }
	}
#if 0 
      if (OprogStatus.Nose==1 && nebrNowLong)
	{
	  nebrNowLong = 0;
	  /* build up linked list on predicted 
	     coordinates (only father do it)*/
	  if (OprogStatus.noLinkedList)
	    {
	      BuildNebrListNoLinkedLong(Oparams.parnum, Oparams.rcut);
	      if (OprogStatus.Nose == 1)
		LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
	    }
	  else
	    {
	      links(Oparams.parnum, Oparams.rcut);
	      /* Build up neighbour list */  
	      BuildNebrListLong(Oparams.parnum, Oparams.rcut);
	      if (OprogStatus.Nose == 1)
		LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
	    }
	}
#endif 
      
      LJForce(Oparams.parnum, OprogStatus.rcutInner);
#ifdef MD_FENE
      FENEForce();
#endif
            /* kinet(Oparams.parnum, vx, vy, vz, Vol1); */
#ifdef MD_RESPA_NPT
      /* NVE ensemble o Dinamica Browniana */
      movebRespa(Oparams.steplength/n, 0.00000000001, 150, NA-1, Oparams.m, distance, 
		 Oparams.parnum); 
#else
      /* correct the coords */
      if (OprogStatus.Nose == 1)
	{  
	  /* NPT ensemble */
	  movebNPT(Oparams.steplength/n, 0.00000000001, 150, NA-1, Oparams.m, distance, 
		   Oparams.parnum);             
	}
      else if (OprogStatus.Nose == 2)
	{
	  /* NPT ensemble */
	  /*movebNTV(Oparams.steplength, 0.0000000001, 30, 3, Oparams.m, distance, 
	    Oparams.parnum); */            
	  movebNTV(Oparams.steplength/n, 0.00000000001, 150, NA-1, Oparams.m, distance, 
	       Oparams.parnum);  
	}
      else
	/* 0 = NVE, -1 = dinamica browniana */ 
	{
	  /* NVE ensemble o Dinamica Browniana */
	  moveb(Oparams.steplength/n, 0.00000000001, 150, NA-1, Oparams.m, distance, 
	    Oparams.parnum);             
	}
#endif

      if (OprogStatus.Nose==1)
	{
	  checkNebrRebuildNPT();
	  checkNebrRebuildNPTLong();
	}
      else
	{
	  checkNebrRebuild();
	  checkNebrRebuildLong();
	}
    }
  if (nebrNowLong)
    {
      nebrNowLong = 0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      if (OprogStatus.noLinkedList)
	{
	  BuildNebrListNoLinkedLong(Oparams.parnum, Oparams.rcut);
	}
      else
	{
	  links(Oparams.parnum, Oparams.rcut);
	  /* Build up neighbour list */  
	  BuildNebrListLong(Oparams.parnum, Oparams.rcut);
	}
    }
  
   //printf("Steps: %d VcR: %f VcL: %f\n",  Oparams.curStep, VcR, VcLong);
  
#ifdef MD_RESPA_NPT
  movelongRespaNPTAft(Oparams.steplength);
  p2v();
#else
  LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
  for (i=0; i < Oparams.parnum; i++)
    for (a=0; a < NA; a++)
      {
	vx[a][i] += 0.5 * dt * FxLong[a][i]/Oparams.m[a];
      	vy[a][i] += 0.5 * dt * FyLong[a][i]/Oparams.m[a];
	vz[a][i] += 0.5 * dt * FzLong[a][i]/Oparams.m[a];
      }
#ifndef MD_FENE  
  shakeVel(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000000001, vx, vy, vz);
#endif
#endif
  calcPressTens();
  kinetRespaNPT(Oparams.parnum, px, py, pz);
  /* Calculate the kinetic energy */
  //kinet(Oparams.parnum, vx, vy, vz, Vol1);

  if ( (OprogStatus.Nose == 1) || (OprogStatus.Nose == 2))
    {
      /*scalCor(Oparams.parnum);*/
      if ( ( (OprogStatus.sResetSteps > 0) &&
	     (Oparams.curStep == OprogStatus.sResetSteps) ) 
	   || ( OprogStatus.sResetSteps < 0 &&
		Oparams.curStep % abs(OprogStatus.sResetSteps) == 0) )

	{
	  /* NOTE:
	     Reset s to 1 and set to zero its time derivatives
	     This is a trick to reduce the energy drift during equilibration,
	     infact if s about unity the total energy drift less */
	  s = 1.0;
	  s1 = 0.0;
	  s2 = 0.0;
	}
      if (OprogStatus.Nose == 1)
	{
	  chksVol();
	}
      else 
	{
	  chks();
	}
    }
  if (OprogStatus.grow)
    {
      update_sigmas();

      if (check_sigmas())
	{
	  printf("Tutti i monomeri hanno raggiunto il raggio desiderato\n");
	  printf("Simulazione terminata...\n");
	  ENDSIM=1;
	}
    }
  if (  ( OprogStatus.snapSteps < 0 && (abs(OprogStatus.snapSteps) == Oparams.curStep) ) || 
	( OprogStatus.snapSteps > 0 && (Oparams.curStep % OprogStatus.snapSteps == 0) )  )
    savesnap();
  /* Update accumulators for calculating the angular diffusion coefficent */
  updateAng(Oparams.parnum);  

  /* Update the integral of the pressure tensor */
  updateDQ(Oparams.steplength);
  if ( ((OprogStatus.CMreset != 0) &&
	/* 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) */
       (Oparams.curStep == OprogStatus.CMreset)) 
    || ( OprogStatus.CMreset < 0 && (abs(Oparams.curStep % OprogStatus.CMreset) == 0) ) ) 
      resetCM(Oparams.parnum);
}
#endif
