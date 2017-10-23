#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
/* posizioni e velocità di tutti gli atomi in un dischetto
 * (cioè anche quello massless!) */ 
extern double **rallx, **rally, **rallz, **Fallx, **Fally, **Fallz,
  **rallx_old, **rally_old, **rallz_old, **vxold, **vyold, **vzold;

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
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Wm, Wmxx, Wmyy, Wmzz, 
  Wmxy, Wmyz, Wmzx, Pmxx, Pmyy, Pmzz, Pmxy, Pmyz, Pmzx, T1mxy, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
  T1myz, T1mzx, T1mxx, T1myy, T1mzz;  
COORD_TYPE DrSq = 0.0, Mtot;
/* used by linked list routines */
int *head, *list, *map;  /* arrays of integer */
int NCell, mapSize, M;

/* neighbour list method variables */
COORD_TYPE dispHi;
int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
COORD_TYPE  *Rmx, *Rmy, *Rmz;
/* ========================================================================= */
extern void check_distances(void);
extern double *atcharge;
/* ========================== >>> movea <<< =============================== */
void movea(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB, COORD_TYPE d, 
	   COORD_TYPE m[NA], int Nm)
{  
  COORD_TYPE dSq;
  int      done;
  int      moving[NA], moved[NA];
  COORD_TYPE  tol2, pxab, pyab, pzab, pabSq, dt2, dtSq2;
  COORD_TYPE  rabSq, diffSq, rxab, ryab, rzab, rpab, gab;
  COORD_TYPE  dx, dy, dz, rma, rmb;
  COORD_TYPE  axia, ayia, azia;
  COORD_TYPE  rxi[NA], ryi[NA], rzi[NA], pxi[NA], pyi[NA], pzi[NA],
    vxi[NA], vyi[NA], vzi[NA];
  int i, a, b, it;
  const COORD_TYPE rptol = 1.0E-6;

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
  tol2   = 2.0 * tol;
  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;
  dSq = Sqr(d); /* In general you must supply a vector of bond lengths */

  /* ===== LOOP OVER MOLECULES ===== */
  for (i=0; i < Nm; i++)
    {
      /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
      for(a=0; a < NA; a++)
	{
	  axia = Fx[a][i] / m[a];
	  ayia = Fy[a][i] / m[a];
	  azia = Fz[a][i] / m[a];
	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];
	  pxi[a] = rx[a][i] + dt * vx[a][i] + dtSq2 * axia;
	  pyi[a] = ry[a][i] + dt * vy[a][i] + dtSq2 * ayia;
	  pzi[a] = rz[a][i] + dt * vz[a][i] + dtSq2 * azia;
#if 0 
	  printf("[%d-%d] v=(%f,%f,%f) a=(%f,%f,%f) pr=(%f,%f,%f) r=(%f,%f,%f)\n", a, i, 
		 vx[a][i], vy[a][i], vz[a][i],
		 axia, ayia, azia, pxi[a], pyi[a], pzi[a],
		 rxi[a], ryi[a], rzi[a]);
#endif
	  vxt[a][i] = vx[a][i]; /* v(t) */
	  vyt[a][i] = vy[a][i];
	  vzt[a][i] = vz[a][i];
	  
	  vxi[a] = vx[a][i] + dt2 * axia;
	  vyi[a] = vy[a][i] + dt2 * ayia;
	  vzi[a] = vz[a][i] + dt2 * azia;
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
	      if (b >= NA) b = 0;
		   
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
		  printf("a=%d b=%d rabSq: %f pabSq: %f\n", a, b, rabSq, pabSq);
		  printf("a=%d (%f,%f,%f) b=%d (%f,%f,%f)\n", a, pxi[a], pyi[a], pzi[a],
			 b, pxi[b], pyi[b], pzi[b]);
#endif
		  if ( fabs(diffSq) > ( rabSq * tol2 ) ) 
		    {
		      rxab = rxi[a] - rxi[b];
		      ryab = ryi[a] - ryi[b];
		      rzab = rzi[a] - rzi[b];
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
				 Sqr(rxab)+Sqr(ryab)+Sqr(rzab));
			  printf("rab(%f,%f,%f) pab(%f,%f,%f)\n", rxab, ryab, rzab,
				 pxab, pyab, pzab);
			   printf("[%d-%d] v=(%f,%f,%f) a=(%f,%f,%f) pr=(%f,%f,%f) r=(%f,%f,%f)\n", a, i, 
				  vx[a][i], vy[a][i], vz[a][i],
				  axia, ayia, azia, pxi[a], pyi[a], pzi[a],
				  rxi[a], ryi[a], rzi[a]);

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
			   
		      dx = dx / dt;
		      dy = dy / dt;
		      dz = dz / dt;
		      
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

	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
	}
      
    }
  /* END OF LOOP OVER MOLECULES */

  if (OprogStatus.Nose == 0) return;
  /* Calculate the friction coefficent at time t and its derivative at time
     t + dt/2 */
  s = s + dt * s1 + dtSq2 * s2;
  s1t = s1;       /* s1(t) */
  s1 = s1 + dt2 * s2;
  
  if (OprogStatus.Nose == 2) return;
  Vol = Vol + dt * Vol1 + dtSq2 * Vol2;
  Vol1t = Vol1;   /* Vol1(t) */
  Vol1 = Vol1 + dt2 * Vol2;
 
}
double calc_norm(double *v)
{
  int a;
  double norm=0;
  for (a = 0; a < 3; a++)
   norm += Sqr(v[a]);
  return sqrt(norm); 
}
#ifdef MD_BROWN_BETTER
void build_ref_axes(int i, double u1[3], double u2[3], double u3[3])
{
  int a;
  double n;
  u1[0] = rx[1][i] - rx[0][i];
  u1[1] = ry[1][i] - ry[0][i];
  u1[2] = rz[1][i] - rz[0][i];
  n = calc_norm(u1);
  for (a = 0; a < 3; a++)
    u1[a] /= n;
  for (a=0; a < 3; a++)
    
  u2[0] = 2.0*rx[2][i]-(rx[0][i] + rx[1][i]);
  u2[1] = 2.0*ry[2][i]-(ry[0][i] + ry[1][i]);
  u2[2] = 2.0*rz[2][i]-(rz[0][i] + rz[1][i]);
  
  n = calc_norm(u2);
  for (a = 0; a < 3; a++)
    u2[a] /= n;
    
  vectProd(u1[0], u1[1], u1[2], u2[0], u2[1], u2[2], &u3[0], &u3[1], &u3[2]);  

}
void lab2body(double u1[3], double u2[3], double u3[3], double vx, double vy, double vz,
	      double *vxp, double *vyp, double *vzp)
{
  *vxp = vx*u1[0]+vy*u1[1]+vz*u1[2];
  *vyp = vx*u2[0]+vy*u2[1]+vz*u2[2];
  *vzp = vx*u3[0]+vy*u3[1]+vz*u3[2];
}
void body2lab(double u1[3], double u2[3], double u3[3], double vxp, double vyp, double vzp,
	      double *vx, double *vy, double *vz)
{
  *vx = vxp*u1[0]+vyp*u2[0]+vzp*u3[0];
  *vy = vxp*u1[1]+vyp*u2[1]+vzp*u3[1];
  *vz = vxp*u1[2]+vyp*u2[2]+vzp*u3[2];
}
extern double **u1, **u2, **u3;
/* ========================== >>> movea <<< =============================== */
void movea_Brownian(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB, COORD_TYPE d, 
	   COORD_TYPE m[NA], int Nm)
{  
  COORD_TYPE dSq;
  int      done;
  int      moving[NA], moved[NA];
  COORD_TYPE  tol2, pxab, pyab, pzab, pabSq, dtSq, dtSq2, dt2;
  COORD_TYPE  rabSq, diffSq, rxab, ryab, rzab, rpab, gab;
  COORD_TYPE  dx, dy, dz, rma, rmb;
  COORD_TYPE  axia, ayia, azia, sigmarN, sigmavN, sigmarP, sigmavP;
  COORD_TYPE  rxi[NA], ryi[NA], rzi[NA], pxi[NA], pyi[NA], pzi[NA],
    vxi[NA], vyi[NA], vzi[NA];
  double drx, dry, drz, dvx, dvy, dvz, kTm; 
  double echsidtP, e2chsidtP, chsiP, chsidtP, crvP, echsidtN, e2chsidtN, chsiN, chsidtN, crvN;
  double c0N, c1N, c2N, c1mc2N, c0P, c1P, c2P, c1mc2P;
  double vxp, vyp, vzp, axp, ayp, azp, rxp, ryp, rzp;
  int i, a, b, it;
  const COORD_TYPE rptol = 1.0E-6;

  if ( ( NB != NA ) && ( NB != NA-1 ) ) 
    {
      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,  
	    "NB IN ERROR", 
	    NULL);
      exit(-1);
    }
        
  L = cbrt(Vol);
  invL = 1.0 / L;

  tol2   = 2.0 * tol;
  dtSq  = dt * dt;
  dt2 = dt / 2.0;
  dtSq2 = dt * dt2;
 
  dSq = Sqr(d); /* In general you must supply a vector of bond lengths */

  kTm = Oparams.T / Oparams.m[0];
  chsidtN = Oparams.chsiN * Oparams.steplength;
  echsidtN  = exp(-Oparams.chsiN*Oparams.steplength);
  e2chsidtN = exp(-2.0*Oparams.chsiN*Oparams.steplength);
  dt = Oparams.steplength;
  chsiN = Oparams.chsiN;
  c0N = exp(-chsiN*dt);
  c1N = (1-c0N)/(chsiN*dt);
  c2N = (1-c1N)/(chsiN*dt);
  c1mc2N = c1N - c2N;
  /* qui si assume che le masse di tutti gli atomi del disco sono uguali! */
  sigmarN = sqrt( dt * (kTm / chsiN) *
		(2.0 -  ( 3.0 - 4.0 * echsidtN + e2chsidtN) / chsidtN) );
  sigmavN = sqrt(kTm * ( 1.0 - e2chsidtN ) );
  crvN = dt * kTm * Sqr( 1.0 - echsidtN) / chsidtN / sigmarN / sigmavN;
  
  chsidtP = Oparams.chsiP * Oparams.steplength;
  echsidtP  = exp(-Oparams.chsiP*Oparams.steplength);
  e2chsidtP = exp(-2.0*Oparams.chsiP*Oparams.steplength);
  chsiP = Oparams.chsiP;
  c0P = exp(-chsiP*dt);
  c1P = (1-c0P)/(chsiP*dt);
  c2P = (1-c1P)/(chsiP*dt);
  c1mc2P = c1P - c2P;
  /* qui si assume che le masse di tutti gli atomi del disco sono uguali! */
  sigmarP = sqrt( dt * (kTm / chsiP) *
		(2.0 -  ( 3.0 - 4.0 * echsidtP + e2chsidtP) / chsidtP) );
  sigmavP = sqrt(kTm * ( 1.0 - e2chsidtP ) );
  crvP = dt * kTm * Sqr( 1.0 - echsidtP) / chsidtP / sigmarP / sigmavP;

  /* ===== LOOP OVER MOLECULES ===== */
  for (i=0; i < Nm; i++)
    {
      build_ref_axes(i, u1[i], u2[i], u3[i]);
      /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
      for(a=0; a < NA; a++)
	{
	  axia = Fx[a][i] / m[a];
	  ayia = Fy[a][i] / m[a];
	  azia = Fz[a][i] / m[a];

	  bigauss(sigmarP, sigmavP, crvP, &drx, &dvx);
          bigauss(sigmarP, sigmavP, crvP, &dry, &dvy);
          bigauss(sigmarN, sigmavN, crvN, &drz, &dvz);

	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];

	  lab2body(u1[i], u2[i], u3[i], rxi[a], ryi[a], rzi[a], &rxp, &ryp, &rzp);
	  lab2body(u1[i], u2[i], u3[i], axia, ayia, azia, &axp, &ayp, &azp);
	  lab2body(u1[i], u2[i], u3[i], vx[a][i], vy[a][i], vz[a][i], &vxp, &vyp, &vzp);

	  rxp = rxp + c1P * dt * vxp + c2P * dtSq * axp + drx;
	  ryp = ryp + c1P * dt * vyp + c2P * dtSq * ayp + dry;
	  rzp = rzp + c1N * dt * vzp + c2N * dtSq * azp + drz;
	  
	  vxt[a][i] = vx[a][i]; /* v(t) */
	  vyt[a][i] = vy[a][i];
	  vzt[a][i] = vz[a][i];
	  
	  vxp = c0P*vxp + c1mc2P * dt * axp + dvx;
	  vyp = c0P*vyp + c1mc2P * dt * ayp + dvy;
	  vzp = c0N*vzp + c1mc2N * dt * azp + dvz;

	  body2lab(u1[i], u2[i], u3[i], rxp, ryp, rzp, &pxi[a], &pyi[a], &pzi[a]);
	  body2lab(u1[i], u2[i], u3[i], vxp, vyp, vzp, &vxi[a], &vyi[a], &vzi[a]);

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
	      if (b >= NA) b = 0;
		   
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
		  
		  if ( fabs(diffSq) > ( rabSq * tol2 ) ) 
		    {
		      rxab = rxi[a] - rxi[b];
		      ryab = ryi[a] - ryi[b];
		      rzab = rzi[a] - rzi[b];
		      rxab = rxab - L * rint(invL * rxab);
		      ryab = ryab - L * rint(invL * ryab);
		      rzab = rzab - L * rint(invL * rzab);
		      rpab = rxab * pxab + ryab * pyab + rzab * pzab;
		      
		      if ( rpab < (rabSq * rptol) )
			{
			
			  mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
				"CONSTRAINT FAILURE IN MOVEA",
				NULL);
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
			   
		      dx = dx / dt;
		      dy = dy / dt;
		      dz = dz / dt;
		      
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

	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
	}
      
    }
  /* END OF LOOP OVER MOLECULES */

  /* -1 = Brownian dyanmics  -2 = Brownian dynamics + Andersen */
  if (OprogStatus.Nose == -1) return;
  Vol = Vol + dt * Vol1 + dtSq2 * Vol2;
  Vol1t = Vol1;   /* Vol1(t) */
  Vol1 = Vol1 + dt2 * Vol2;
 
}
#else
/* ========================== >>> movea <<< =============================== */
void movea_Brownian(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB, COORD_TYPE d, 
	   COORD_TYPE m[NA], int Nm)
{  
  COORD_TYPE dSq;
  int      done;
  int      moving[NA], moved[NA];
  COORD_TYPE  tol2, pxab, pyab, pzab, pabSq, dtSq, dtSq2, dt2;
  COORD_TYPE  rabSq, diffSq, rxab, ryab, rzab, rpab, gab;
  COORD_TYPE  dx, dy, dz, rma, rmb;
  COORD_TYPE  axia, ayia, azia, sigmar, sigmav;
  COORD_TYPE  rxi[NA], ryi[NA], rzi[NA], pxi[NA], pyi[NA], pzi[NA],
    vxi[NA], vyi[NA], vzi[NA];
  double drx, dry, drz, dvx, dvy, dvz, kTm, echsidt, e2chsidt, chsi, chsidt, crv;
  double c0, c1, c2, c1mc2;
  int i, a, b, it;
  const COORD_TYPE rptol = 1.0E-6;

  if ( ( NB != NA ) && ( NB != NA-1 ) ) 
    {
      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,  
	    "NB IN ERROR", 
	    NULL);
      exit(-1);
    }
        
  L = cbrt(Vol);
  invL = 1.0 / L;

  tol2   = 2.0 * tol;
  dtSq  = dt * dt;
  dt2 = dt / 2.0;
  dtSq2 = dt * dt2;
  
  dSq = Sqr(d); /* In general you must supply a vector of bond lengths */

  kTm = Oparams.T / Oparams.m[0];
  chsidt = Oparams.chsi * Oparams.steplength;
  echsidt  = exp(-Oparams.chsi*Oparams.steplength);
  e2chsidt = exp(-2.0*Oparams.chsi*Oparams.steplength);
  dt = Oparams.steplength;
  chsi = Oparams.chsi;
  c0 = exp(-chsi*dt);
  c1 = (1-c0)/(chsi*dt);
  c2 = (1-c1)/(chsi*dt);
  c1mc2 = c1 - c2;
  /* qui si assume che le masse di tutti gli atomi del disco sono uguali! */
  sigmar = sqrt( dt * (kTm / chsi) *
		(2.0 -  ( 3.0 - 4.0 * echsidt + e2chsidt) / chsidt) );
  sigmav = sqrt(kTm * ( 1.0 - e2chsidt ) );
  crv = dt * kTm * Sqr( 1.0 - echsidt) / chsidt / sigmar / sigmav;

  /* ===== LOOP OVER MOLECULES ===== */
  for (i=0; i < Nm; i++)
    {
      /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
      for(a=0; a < NA; a++)
	{
	  axia = Fx[a][i] / m[a];
	  ayia = Fy[a][i] / m[a];
	  azia = Fz[a][i] / m[a];

	  bigauss(sigmar, sigmav, crv, &drx, &dvx);
          bigauss(sigmar, sigmav, crv, &dry, &dvy);
          bigauss(sigmar, sigmav, crv, &drz, &dvz);

	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];
	  pxi[a] = rx[a][i] + c1 * dt * vx[a][i] + c2 * dtSq * axia + drx;
	  pyi[a] = ry[a][i] + c1 * dt * vy[a][i] + c2 * dtSq * ayia + dry;
	  pzi[a] = rz[a][i] + c1 * dt * vz[a][i] + c2 * dtSq * azia + drz;
	  
	  vxt[a][i] = vx[a][i]; /* v(t) */
	  vyt[a][i] = vy[a][i];
	  vzt[a][i] = vz[a][i];
	  
	  vxi[a] = c0*vx[a][i] + c1mc2 * dt * axia + dvx;
	  vyi[a] = c0*vy[a][i] + c1mc2 * dt * ayia + dvy;
	  vzi[a] = c0*vz[a][i] + c1mc2 * dt * azia + dvz;

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
	      if (b >= NA) b = 0;
		   
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
		  
		  if ( fabs(diffSq) > ( rabSq * tol2 ) ) 
		    {
		      rxab = rxi[a] - rxi[b];
		      ryab = ryi[a] - ryi[b];
		      rzab = rzi[a] - rzi[b];
		      rxab = rxab - L * rint(invL * rxab);
		      ryab = ryab - L * rint(invL * ryab);
		      rzab = rzab - L * rint(invL * rzab);
		      rpab = rxab * pxab + ryab * pyab + rzab * pzab;
		      
		      if ( rpab < (rabSq * rptol) )
			{
			
			  mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
				"CONSTRAINT FAILURE IN MOVEA",
				NULL);
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
			   
		      dx = dx / dt;
		      dy = dy / dt;
		      dz = dz / dt;
		      
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

	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
	}
      
    }
  /* END OF LOOP OVER MOLECULES */

  /* -1 = Brownian dyanmics  -2 = Brownian dynamics + Andersen */
  if (OprogStatus.Nose == -1) return;
  Vol = Vol + dt * Vol1 + dtSq2 * Vol2;
  Vol1t = Vol1;   /* Vol1(t) */
  Vol1 = Vol1 + dt2 * Vol2;
 
}
#endif

/* ============================ >>> shakeVel <<< ========================== */
void shakeVel(int Nm, COORD_TYPE dt, COORD_TYPE m[NA], int maxIt, int NB, 
	      COORD_TYPE d, COORD_TYPE tol, COORD_TYPE **v2sx, 
	      COORD_TYPE** v2sy, COORD_TYPE** v2sz )
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

	  vxi[a] = v2sx[a][i];
	  vyi[a] = v2sy[a][i];
	  vzi[a] = v2sz[a][i];
	  
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
		  if (b >= NA) b = 0;
		  
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
	      v2sx[a][i] = vxi[a];
	      v2sy[a][i] = vyi[a];
	      v2sz[a][i] = vzi[a];
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

/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagMol(int Nm, COORD_TYPE VOL1)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of molecular pressure tensor */
  int i;
  COORD_TYPE kin, dlnV;
  COORD_TYPE Vx, Vy, Vz;

  kin = 0.0;
  dlnV = VOL1 / Vol / 3.0;
  for(i=0; i < Nm; i++)
    {
      CoMV(i, &Vx, &Vy, &Vz);
      kin +=  Sqr(Vx - dlnV * Rx[i]) + Sqr(Vy - dlnV * Ry[i]) + 
	Sqr(Vz - dlnV * Rz[i]);
    }
  
  kin *= Mtot;
  kin /= 3.0 * Vol;
  return kin;

}
/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagAt(int Nm, COORD_TYPE VOL1)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of atomic pressure tensor */
  int i, a;
  COORD_TYPE kin, dlnV, *m;
  
  m = Oparams.m;
  kin = 0.0;
  dlnV = VOL1 / Vol / 3.0;
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  kin +=  m[a]*(Sqr(vx[a][i] - dlnV * Rx[i]) + 
	    Sqr(vy[a][i] - dlnV * Ry[i]) + 
	    Sqr(vz[a][i] - dlnV * Rz[i]));
	}
    }
  kin /= 3.0 * Vol;
  return kin;
}
/*========================== >>> calcPtensMol <<< ===========================*/
void calcPtensMol(int Nm, COORD_TYPE VOL1)
{
  int i, a;
  COORD_TYPE vmx, vmy, vmz, dlnV, Px, Py, Pz, *m;
  
  T1mxy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1myz = 0.0;
  T1mzx = 0.0;
  T1mxx = 0.0; /* Kinetic diagonal terms of pressure tensor */
  T1myy = 0.0;
  T1mzz = 0.0;
  m = Oparams.m;
  dlnV = VOL1 / Vol / 3.0;
  for(i=0; i < Nm; i++)
    {
      vmx = vmy = vmz = 0.0;  /* Velocity of center of mass of molecule i */
      for(a=0; a < NA; a++)
	{
	  vmx += vx[a][i] * m[a];
	  vmy += vy[a][i] * m[a];
	  vmz += vz[a][i] * m[a];
	}
      /* Kinetic component of pressure tensor (all terms) */
      Px = vmx - dlnV * Rx[i];
      Py = vmy - dlnV * Ry[i];
      Pz = vmz - dlnV * Rz[i];

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
void calcPtensAt(int Nm, COORD_TYPE VOL1)
{
  /* Calculate all components of atomic pressure tensor */
  int i, a;
  COORD_TYPE dlnV, px, py, pz, *m;
  
  m = Oparams.m;
  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;

  dlnV = VOL1 / Vol / 3.0;
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  px = vx[a][i] - dlnV * Rx[i];
	  py = vy[a][i] - dlnV * Ry[i];
	  pz = vz[a][i] - dlnV * Rz[i];
	  /* Kinetic component of pressure tensor (all terms) */
	  T1xy += px * py * m[a]; 
	  T1yz += py * pz * m[a];
	  T1zx += pz * px * m[a] ;
	  T1xx += px * px * m[a];
	  T1yy += py * py * m[a];
	  T1zz += pz * pz * m[a];
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
  //printf("Wxx:%f WCxx:%f Wyy: %f WCyy: %f Wzz: %f WCzz: %f\n",
	//Wxx, WCxx, Wyy, WCyy, Wzz, WCzz);
  Patxx /= Vol;
  Patyy /= Vol;
  Patzz /= Vol;
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
}

const double ittol = 1E-20;
#ifdef MD_BROWN_BETTER
void movebBrownAnd(double dt, double tol, int maxIt, int NB, double m[3], double d, 
		    int Nm) 
{
  /* NTP fatto con dinamica browniana e pistone meccanico alla Andersen */ 
  COORD_TYPE FxNose, FyNose, FzNose, FxAnd, FyAnd, FzAnd, cfAnd, dlns, dlnV,
    dlnVSq; 
  COORD_TYPE Vol1g, s1i, Vol1i;
  COORD_TYPE DT, A, B, DP, c0N, c1N, c2N, chsiN, c0P, c1P, c2P, chsiP; 
  double vxp, vyp, vzp, Fxp, Fyp, Fzp;
  double u1[3], u2[3], u3[3];
  int i, a, k, numok;
  const int MAXNUMIT = 20;

  /* ******************************************************************* */
  chsiP = Oparams.chsiP;
  chsiN = Oparams.chsiN;
  c0N = exp(-chsiN*dt);
  c1N = (1-c0N)/(chsiN*dt);
  c2N = (1-c1N)/(chsiN*dt);
  c0P = exp(-chsiP*dt);
  c1P = (1-c0P)/(chsiP*dt);
  c2P = (1-c1P)/(chsiP*dt);

  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx[i], &Ry[i], &Rz[i]);
      for(a=0; a < NA; a++)
	{
	  vxt2[a][i] = vx[a][i]; /* v*t2 = v*(t+dt/2) */
	  vyt2[a][i] = vy[a][i];
	  vzt2[a][i] = vz[a][i];

	  FxLJ[a][i] = Fx[a][i];
	  FyLJ[a][i] = Fy[a][i];
	  FzLJ[a][i] = Fz[a][i];
	
	  /* predicted values of velocities at time t+dt */
	  /*vx[a][i] = 13.0 * vxt[a][i] / 6.0 - 4.0 * vxo1[a][i] / 3.0 + 
	    vxo2[a][i] / 6.0;
	    vy[a][i] = 13.0 * vyt[a][i] / 6.0 - 4.0 * vyo1[a][i] / 3.0 + 
	    vyo2[a][i] / 6.0;
	    vz[a][i] = 13.0 * vzt[a][i] / 6.0 - 4.0 * vzo1[a][i] / 3.0 + 
	    vzo2[a][i] / 6.0;
	  */
	  /*vx[a][i] = 2.0 * vxt[a][i] - vxo1[a][i]; */ /* Verlet */
	  /*vy[a][i] = 2.0 * vyt[a][i] - vyo1[a][i];
	    vz[a][i] = 2.0 * vzt[a][i] - vzo1[a][i];*/

	  vx[a][i] = 5.0 * vxt[a][i] / 2.0 - 2.0 * vxo1[a][i] + 
	    vxo2[a][i] / 2.0;
	  vy[a][i] = 5.0 * vyt[a][i] / 2.0 - 2.0 * vyo1[a][i] + 
	    vyo2[a][i] / 2.0;
	  vz[a][i] = 5.0 * vzt[a][i] / 2.0 - 2.0 * vzo1[a][i] + 
	    vzo2[a][i] / 2.0;
	}
    }
#if 0
  s1i = s1;    /* s1i = s1(t+dt/2) */
#endif
  Vol1i = Vol1;/* Vol1i = Vol1(t+dt/2)*/

  /* Initial guess for Vol1 */
  /* Vol1 = 13.0*Vol1t/6.0 - 4.0*Vol1o1/3.0 + Vol1o2 / 6.0; */ /* Vol1(t+dt) */
  
  /* Vol1 = 2.0 * Vol1t - Vol1o1; *//* Verlet */  
  
  Vol1 = 5.0 * Vol1t / 2.0 - 2.0 * Vol1o1 + Vol1o2 / 2.0;  

  shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
  for(k=0; k < MAXNUMIT; k++) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
#if 0
      kinet(Nm, vx, vy, vz, Vol1);  /* K(t+dt) */

      DT = s * (2.0 * K - (5.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
#endif 
      /* Calculate pressure, calcT1diagMol is a term proportional to the 
	 translational kinetic energy, see Ferrario and Ryckaert */
#ifdef MOLPTENS
      press_m = calcT1diagMol(Nm, Vol1) + Wm / 3.0 / Vol; /* press(t+dt) */
      
      /* Volume acceleration */
      DP = Sqr(s) * (press_m - Oparams.P) / OprogStatus.W;
#else
      press_at = calcT1diagAt(Nm, Vol1) + (W + WC) / Vol; /* press(t+dt) */
      /* Volume acceleration */
      DP = Sqr(s) * (press_at - Oparams.P) / OprogStatus.W;
#endif
      B = Vol1i + 0.5 * dt * DP;
      Vol1 = B / (1.0 - 0.5 * dt * s1 / s); /* Vol1(t+dt) */
      Vol2 = DP + s1 * Vol1 / s;            /* Vol2(t+dt) */
      
      /* Calculate the Andersen forces */
      dlnV   = Vol1 / Vol; 
      dlnVSq = Sqr(dlnV);
      cfAnd  = - 2.0 * dlnVSq / 9.0;
      cfAnd += Vol2 / Vol / 3.0;
#if 0
      cfAnd += dlns * dlnV / 3.0;
#endif	 
      for(i=0; i < Nm; i++)
	{
	  build_ref_axes(i, u1, u2, u3);
	  for(a=0; a < NA; a++)
	    {
	      /* Calculate center of mass position for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	    
	      FxAnd  = cfAnd * Rx[i] * m[a];
	      FyAnd  = cfAnd * Ry[i] * m[a];
	      FzAnd  = cfAnd * Rz[i] * m[a];
	      Fx[a][i] = FxLJ[a][i] + FxAnd;
	      Fy[a][i] = FyLJ[a][i] + FyAnd;
	      Fz[a][i] = FzLJ[a][i] + FzAnd;
	      /* vt2 = v(t+dt/2) */
	      vxold[a][i] = vx[a][i];
	      vyold[a][i] = vy[a][i];
	      vzold[a][i] = vz[a][i];
	      lab2body(u1, u2, u3, vxt2[a][i], vyt2[a][i], vzt2[a][i],
		       &vxp, &vyp, &vzp);
	      lab2body(u1, u2, u3, Fx[a][i], Fy[a][i], Fx[a][i],
		       &Fxp, &Fyp, &Fzp);
	      vxp = vxp + c2P * dt * Fxp / m[a];
	      vyp = vyp + c2P * dt * Fyp / m[a];
	      vzp = vzp + c2N * dt * Fzp / m[a];
	      body2lab(u1, u2, u3, vxp, vyp, vzp, &vx[a][i], &vy[a][i], &vz[a][i]);
#if 0
	      vx[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vy[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vz[a][i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[a][i] * m[a];
	      FyNose = - dlns * vy[a][i] * m[a];
	      FzNose = - dlns * vz[a][i] * m[a];
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[a][i] = Fx[a][i] + FxNose;
	      Fy[a][i] = Fy[a][i] + FyNose;
	      Fz[a][i] = Fz[a][i] + FzNose;
#endif
	    } 
	}
      shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
      /* Zero the velocity along bond and calculate the constraint virial, this
	 should be as before because the Andersen Forces don't act along 
	 bond.*/
      numok = 0;
      for (i=0; i < Oparams.parnum; i++)
	{
	  for (a = 0; a < NA; a++)
	    {
	      /* F = Flennard-jones + Fnose + Fandersen */
	      if (fabs(vx[a][i]- vxold[a][i]) < ittol &&
		  fabs(vy[a][i]- vyold[a][i]) < ittol &&
		  fabs(vz[a][i]- vzold[a][i]) < ittol)
		numok++;
	    }
	}
      if (numok == NA*Oparams.parnum)
	{
	  /* for all atoms diff between current velocities and prev vels is below tol! */ 
	  /*printf("Done %d iterations instead of %d\n", k, NUMIT);
	  */
	  break;
	}


    }

  /* Calculate kinetic energy */
  kinet(Nm, vx, vy, vz, Vol1);/* K(t+dt) */
  //printf("K exact: %.20f\n", K);

#ifdef MOLPTENS
  /* Calculate all components of molecular pressure tensor */
  calcPtensMol(Nm, Vol1);
  press_m = (Pmxx + Pmyy + Pmzz) / 3.0; /* press(t+dt) */
#endif

#if !defined(MOLPTENS) || defined(ATPTENS)
  /* Calculate all components of atomic pressure tensor */
  calcPtensAt(Nm, Vol1);
  press_at = (Patxx + Patyy + Patzz) / 3.0;
#endif 

#if !defined(MOLPTENS) && !defined(ATPTENS) && defined(ATPRESS)
  /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
  press_at = calcT1diagAt(Nm, Vol1) + (W + WC) / Vol;
  /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
#endif

#ifdef MOLPTENS
  /* Atomic pressure */
  press = press_m;
  Pxy = Pmxy;
  Pyz = Pmyz;
  Pzx = Pmzx;
#else
  press = press_at;
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;
#endif

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  Vol1o2 = Vol1o1;
#if 0
  s1o2  = s1o1;
#endif
  Vol1o1 = Vol1t;
#if 0
  s1o1 = s1t;
#endif
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }
  /* END OF ITERATION LOOP TO IMPROVE ANDERSEN-NOSE FORCE ESTIMATION */

}
#else
void movebBrownAnd(double dt, double tol, int maxIt, int NB, double m[3], double d, 
		    int Nm) 
{
  /* NTP fatto con dinamica browniana e pistone meccanico alla Andersen */ 
  COORD_TYPE FxNose, FyNose, FzNose, FxAnd, FyAnd, FzAnd, cfAnd, dlns, dlnV,
    dlnVSq; 
  COORD_TYPE Vol1g, s1i, Vol1i;
  COORD_TYPE DT, A, B, DP, c0, c1, c2, chsi; 
  int i, a, k, numok;
  const int MAXNUMIT = 20;

  /* ******************************************************************* */
  chsi = Oparams.chsi;
  c0 = exp(-chsi*dt);
  c1 = (1-c0)/(chsi*dt);
  c2 = (1-c1)/(chsi*dt);

  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx[i], &Ry[i], &Rz[i]);
      for(a=0; a < NA; a++)
	{
	  vxt2[a][i] = vx[a][i]; /* v*t2 = v*(t+dt/2) */
	  vyt2[a][i] = vy[a][i];
	  vzt2[a][i] = vz[a][i];

	  FxLJ[a][i] = Fx[a][i];
	  FyLJ[a][i] = Fy[a][i];
	  FzLJ[a][i] = Fz[a][i];
	
	  /* predicted values of velocities at time t+dt */
	  /*vx[a][i] = 13.0 * vxt[a][i] / 6.0 - 4.0 * vxo1[a][i] / 3.0 + 
	    vxo2[a][i] / 6.0;
	    vy[a][i] = 13.0 * vyt[a][i] / 6.0 - 4.0 * vyo1[a][i] / 3.0 + 
	    vyo2[a][i] / 6.0;
	    vz[a][i] = 13.0 * vzt[a][i] / 6.0 - 4.0 * vzo1[a][i] / 3.0 + 
	    vzo2[a][i] / 6.0;
	  */
	  /*vx[a][i] = 2.0 * vxt[a][i] - vxo1[a][i]; */ /* Verlet */
	  /*vy[a][i] = 2.0 * vyt[a][i] - vyo1[a][i];
	    vz[a][i] = 2.0 * vzt[a][i] - vzo1[a][i];*/

	  vx[a][i] = 5.0 * vxt[a][i] / 2.0 - 2.0 * vxo1[a][i] + 
	    vxo2[a][i] / 2.0;
	  vy[a][i] = 5.0 * vyt[a][i] / 2.0 - 2.0 * vyo1[a][i] + 
	    vyo2[a][i] / 2.0;
	  vz[a][i] = 5.0 * vzt[a][i] / 2.0 - 2.0 * vzo1[a][i] + 
	    vzo2[a][i] / 2.0;
	}
    }
#if 0
  s1i = s1;    /* s1i = s1(t+dt/2) */
#endif
  Vol1i = Vol1;/* Vol1i = Vol1(t+dt/2)*/

  /* Initial guess for Vol1 */
  /* Vol1 = 13.0*Vol1t/6.0 - 4.0*Vol1o1/3.0 + Vol1o2 / 6.0; */ /* Vol1(t+dt) */
  
  /* Vol1 = 2.0 * Vol1t - Vol1o1; *//* Verlet */  
  
  Vol1 = 5.0 * Vol1t / 2.0 - 2.0 * Vol1o1 + Vol1o2 / 2.0;  

  shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
  for(k=0; k < MAXNUMIT; k++) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
#if 0
      kinet(Nm, vx, vy, vz, Vol1);  /* K(t+dt) */

      DT = s * (2.0 * K - (5.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
#endif 
      /* Calculate pressure, calcT1diagMol is a term proportional to the 
	 translational kinetic energy, see Ferrario and Ryckaert */
#ifdef MOLPTENS
      press_m = calcT1diagMol(Nm, Vol1) + Wm / 3.0 / Vol; /* press(t+dt) */
      
      /* Volume acceleration */
      DP = Sqr(s) * (press_m - Oparams.P) / OprogStatus.W;
#else
      press_at = calcT1diagAt(Nm, Vol1) + (W + WC) / Vol; /* press(t+dt) */
      /* Volume acceleration */
      DP = Sqr(s) * (press_at - Oparams.P) / OprogStatus.W;
#endif
      B = Vol1i + 0.5 * dt * DP;
      Vol1 = B / (1.0 - 0.5 * dt * s1 / s); /* Vol1(t+dt) */
      Vol2 = DP + s1 * Vol1 / s;            /* Vol2(t+dt) */
      
      /* Calculate the Andersen forces */
      dlnV   = Vol1 / Vol; 
      dlnVSq = Sqr(dlnV);
      cfAnd  = - 2.0 * dlnVSq / 9.0;
      cfAnd += Vol2 / Vol / 3.0;
#if 0
      cfAnd += dlns * dlnV / 3.0;
#endif	 
      for(i=0; i < Nm; i++)
	{
	  for(a=0; a < NA; a++)
	    {
	      /* Calculate center of mass position for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	    
	      FxAnd  = cfAnd * Rx[i] * m[a];
	      FyAnd  = cfAnd * Ry[i] * m[a];
	      FzAnd  = cfAnd * Rz[i] * m[a];
	      Fx[a][i] = FxLJ[a][i] + FxAnd;
	      Fy[a][i] = FyLJ[a][i] + FyAnd;
	      Fz[a][i] = FzLJ[a][i] + FzAnd;
	      /* vt2 = v(t+dt/2) */
	      vxold[a][i] = vx[a][i];
	      vyold[a][i] = vy[a][i];
	      vzold[a][i] = vz[a][i];

	      vx[a][i] = vxt2[a][i] + c2 * dt * Fx[a][i] / m[a];
	      vy[a][i] = vyt2[a][i] + c2 * dt * Fy[a][i] / m[a];
	      vz[a][i] = vzt2[a][i] + c2 * dt * Fz[a][i] / m[a];
#if 0
	      vx[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vy[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vz[a][i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[a][i] * m[a];
	      FyNose = - dlns * vy[a][i] * m[a];
	      FzNose = - dlns * vz[a][i] * m[a];
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[a][i] = Fx[a][i] + FxNose;
	      Fy[a][i] = Fy[a][i] + FyNose;
	      Fz[a][i] = Fz[a][i] + FzNose;
#endif
	    } 
	}
      shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
      /* Zero the velocity along bond and calculate the constraint virial, this
	 should be as before because the Andersen Forces don't act along 
	 bond.*/
      numok = 0;
      for (i=0; i < Oparams.parnum; i++)
	{
	  for (a = 0; a < NA; a++)
	    {
	      /* F = Flennard-jones + Fnose + Fandersen */
	      if (fabs(vx[a][i]- vxold[a][i]) < ittol &&
		  fabs(vy[a][i]- vyold[a][i]) < ittol &&
		  fabs(vz[a][i]- vzold[a][i]) < ittol)
		numok++;
	    }
	}
      if (numok == NA*Oparams.parnum)
	{
	  /* for all atoms diff between current velocities and prev vels is below tol! */ 
	  /*printf("Done %d iterations instead of %d\n", k, NUMIT);
	  */
	  break;
	}


    }

  /* Calculate kinetic energy */
  kinet(Nm, vx, vy, vz, Vol1);/* K(t+dt) */
  //printf("K exact: %.20f\n", K);

#ifdef MOLPTENS
  /* Calculate all components of molecular pressure tensor */
  calcPtensMol(Nm, Vol1);
  press_m = (Pmxx + Pmyy + Pmzz) / 3.0; /* press(t+dt) */
#endif

#if !defined(MOLPTENS) || defined(ATPTENS)
  /* Calculate all components of atomic pressure tensor */
  calcPtensAt(Nm, Vol1);
  press_at = (Patxx + Patyy + Patzz) / 3.0;
#endif 

#if !defined(MOLPTENS) && !defined(ATPTENS) && defined(ATPRESS)
  /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
  press_at = calcT1diagAt(Nm, Vol1) + (W + WC) / Vol;
  /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
#endif

#ifdef MOLPTENS
  /* Atomic pressure */
  press = press_m;
  Pxy = Pmxy;
  Pyz = Pmyz;
  Pzx = Pmzx;
#else
  press = press_at;
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;
#endif

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  Vol1o2 = Vol1o1;
#if 0
  s1o2  = s1o1;
#endif
  Vol1o1 = Vol1t;
#if 0
  s1o1 = s1t;
#endif
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }
  /* END OF ITERATION LOOP TO IMPROVE ANDERSEN-NOSE FORCE ESTIMATION */

}
#endif
/* ========================= >>> moveb <<< =========================== */
void movebNTV(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	      COORD_TYPE m[3], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, dlns; 
  COORD_TYPE s1i;
  COORD_TYPE DT, A, dt2; 
  int i, a, k, numok;
  const int MAXNUMIT = 40;
  /* ******************************************************************* */
  dt2 = dt * 0.5;
  
  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx[i], &Ry[i], &Rz[i]);
      for(a=0; a < NA; a++)
	{
	  vxt2[a][i] = vx[a][i]; /* v*t2 = v*(t+dt/2) */
	  vyt2[a][i] = vy[a][i];
	  vzt2[a][i] = vz[a][i];

	  FxLJ[a][i] = Fx[a][i];
	  FyLJ[a][i] = Fy[a][i];
	  FzLJ[a][i] = Fz[a][i];
	
	  /* predicted values of velocities at time t+dt */
	  //vx[a][i] = 13.0 * vxt[a][i] / 6.0 - 4.0 * vxo1[a][i] / 3.0 + 
	    //vxo2[a][i] / 6.0;
	  //vy[a][i] = 13.0 * vyt[a][i] / 6.0 - 4.0 * vyo1[a][i] / 3.0 + 
	    //vyo2[a][i] / 6.0;
	  //vz[a][i] = 13.0 * vzt[a][i] / 6.0 - 4.0 * vzo1[a][i] / 3.0 + 
	    //vzo2[a][i] / 6.0;
	  
	  //vx[a][i] = 2.0 * vxt[a][i] - vxo1[a][i]; /* Verlet */
	  //vy[a][i] = 2.0 * vyt[a][i] - vyo1[a][i];
	  //vz[a][i] = 2.0 * vzt[a][i] - vzo1[a][i];
	
	  vx[a][i] = 5.0 * vxt[a][i] / 2.0 - 2.0 * vxo1[a][i] + 
	    vxo2[a][i] / 2.0;
	  vy[a][i] = 5.0 * vyt[a][i] / 2.0 - 2.0 * vyo1[a][i] + 
	    vyo2[a][i] / 2.0;
	  vz[a][i] = 5.0 * vzt[a][i] / 2.0 - 2.0 * vzo1[a][i] + 
	    vzo2[a][i] / 2.0;
	}
    }
  s1i = s1;    /* s1i = s1(t+dt/2) */

  shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
  //printf("Vol1: %f Vol1o1: %f Vol1o2: %f Vol1t:%f\n", 
	// Vol1, Vol1o1, Vol1o2, Vol1t);
  for(k=0; k < MAXNUMIT; k++) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm, vx, vy, vz, 0.0);  /* K(t+dt) */

      DT = s * (2.0 * K - (6.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;
#if 0
      printf("K guess: %.20f s1=%.15f s2=%.15f Oparams.T=%f DT=%f\n", K, s1, s2, Oparams.T, DT);
#endif
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      numok=0; 
      for(i=0; i < Nm; i++)
	{
	  for(a=0; a < NA; a++)
	    {
	      /* Calculate center of mass posiotin for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	      Fx[a][i] = FxLJ[a][i];
	      Fy[a][i] = FyLJ[a][i];
	      Fz[a][i] = FzLJ[a][i];
	     
	      /* vt2 = v(t+dt/2) */
	      vxold[a][i] = vx[a][i];
	      vyold[a][i] = vy[a][i];
	      vzold[a][i] = vz[a][i];

	      vx[a][i] = vxt2[a][i] + dt2 * Fx[a][i] / m[a];
	      vy[a][i] = vyt2[a][i] + dt2 * Fy[a][i] / m[a];
	      vz[a][i] = vzt2[a][i] + dt2 * Fz[a][i] / m[a];
	      vx[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vy[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vz[a][i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[a][i] * m[a];
	      FyNose = - dlns * vy[a][i] * m[a];
	      FzNose = - dlns * vz[a][i] * m[a];
      	      Fx[a][i] = FxLJ[a][i] + FxNose;
	      Fy[a][i] = FyLJ[a][i] + FyNose;
	      Fz[a][i] = FzLJ[a][i] + FzNose;
	    } 
	}
      shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
      for (i=0; i < Oparams.parnum; i++)
	{
	  for (a = 0; a < NA; a++)
	    {
	      /* F = Flennard-jones + Fnose + Fandersen */
	      if (fabs(vx[a][i]- vxold[a][i]) < ittol &&
		  fabs(vy[a][i]- vyold[a][i]) < ittol &&
		  fabs(vz[a][i]- vzold[a][i]) < ittol)
		numok++;
	      /*printf("[%d][%d] %f,%f,%f\n", a, i, vx[a][i], vy[a][i], vz[a][i]);*/
	    }
	}
      if (numok == NA*Oparams.parnum)
	{
	  /*printf("Done %d iterations instead of %d\n", k, NUMIT);
	  */
	  break;
	}
      /* Zero the velocity along bond and calculate the constraint virial, this
	 should be as before because the Andersen Forces don't act along 
	 bond.*/
    }
  kinet(Nm, vx, vy, vz, 0.0);/* K(t+dt) */

#ifdef MOLPTENS
  /* Calculate all components of molecular pressure tensor */
  calcPtensMol(Nm, 0.0);
  /* Molecular pressure */
  press_m = (Pmxx + Pmyy + Pmzz) / 3.0; /* press(t+dt) */
#endif

  /* Calculate all components of atomic pressure tensor */
#if !defined(MOLPTENS) || defined(ATPTENS)
  calcPtensAt(Nm, 0.0);
  press_at = (Patxx + Patyy + Patzz) / 3.0;
#endif
#if !defined(MOLPTENS) && !defined(ATPTENS) && defined(ATPRESS)
  /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
  press_at = calcT1diagAt(Nm, 0.0) + (W + WC) / Vol;
  /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
#endif

//printf("K exact: %.20f\n", K);

/* Atomic pressure */
#ifdef MOLPTENS
  press = press_m;
  Pxy = Pmxy;
  Pyz = Pmyz;
  Pzx = Pmzx;
#else
  press = press_at;
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;
#endif

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  s1o2  = s1o1;
  s1o1 = s1t;
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }
  /* END OF ITERATION LOOP TO IMPROVE ANDERSEN-NOSE FORCE ESTIMATION */
}

/* ========================= >>> moveb <<< =========================== */
void movebNPT(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	      COORD_TYPE m[NA], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, FxAnd, FyAnd, FzAnd, cfAnd, dlns, dlnV,
    dlnVSq; 
  COORD_TYPE Vol1g, s1i, Vol1i;
  COORD_TYPE DT, A, B, DP, dt2; 
  int i, a, k, numok;
  const int MAXNUMIT = 20;

  /* ******************************************************************* */
  dt2 = dt / 2.0;
  
  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx[i], &Ry[i], &Rz[i]);
      for(a=0; a < NA; a++)
	{
	  vxt2[a][i] = vx[a][i]; /* v*t2 = v*(t+dt/2) */
	  vyt2[a][i] = vy[a][i];
	  vzt2[a][i] = vz[a][i];

	  FxLJ[a][i] = Fx[a][i];
	  FyLJ[a][i] = Fy[a][i];
	  FzLJ[a][i] = Fz[a][i];
	
	  /* predicted values of velocities at time t+dt */
	  /*vx[a][i] = 13.0 * vxt[a][i] / 6.0 - 4.0 * vxo1[a][i] / 3.0 + 
	    vxo2[a][i] / 6.0;
	    vy[a][i] = 13.0 * vyt[a][i] / 6.0 - 4.0 * vyo1[a][i] / 3.0 + 
	    vyo2[a][i] / 6.0;
	    vz[a][i] = 13.0 * vzt[a][i] / 6.0 - 4.0 * vzo1[a][i] / 3.0 + 
	    vzo2[a][i] / 6.0;
	  */
	  /*vx[a][i] = 2.0 * vxt[a][i] - vxo1[a][i]; */ /* Verlet */
	  /*vy[a][i] = 2.0 * vyt[a][i] - vyo1[a][i];
	    vz[a][i] = 2.0 * vzt[a][i] - vzo1[a][i];*/

	  vx[a][i] = 5.0 * vxt[a][i] / 2.0 - 2.0 * vxo1[a][i] + 
	    vxo2[a][i] / 2.0;
	  vy[a][i] = 5.0 * vyt[a][i] / 2.0 - 2.0 * vyo1[a][i] + 
	    vyo2[a][i] / 2.0;
	  vz[a][i] = 5.0 * vzt[a][i] / 2.0 - 2.0 * vzo1[a][i] + 
	    vzo2[a][i] / 2.0;
	}
    }
  s1i = s1;    /* s1i = s1(t+dt/2) */
  Vol1i = Vol1;/* Vol1i = Vol1(t+dt/2)*/

  /* Initial guess for Vol1 */
  /* Vol1 = 13.0*Vol1t/6.0 - 4.0*Vol1o1/3.0 + Vol1o2 / 6.0; */ /* Vol1(t+dt) */
  
  /* Vol1 = 2.0 * Vol1t - Vol1o1; *//* Verlet */  
  
  Vol1 = 5.0 * Vol1t / 2.0 - 2.0 * Vol1o1 + Vol1o2 / 2.0;  

  shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
  for(k=0; k < MAXNUMIT; k++) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm, vx, vy, vz, Vol1);  /* K(t+dt) */

      DT = s * (2.0 * K - (6.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      
      /* Calculate pressure, calcT1diagMol is a term proportional to the 
	 translational kinetic energy, see Ferrario and Ryckaert */
#ifdef MOLPTENS
      press_m = calcT1diagMol(Nm, Vol1) + Wm / 3.0 / Vol; /* press(t+dt) */
      
      /* Volume acceleration */
      DP = Sqr(s) * (press_m - Oparams.P) / OprogStatus.W;
#else
      press_at = calcT1diagAt(Nm, Vol1) + (W + WC) / Vol; /* press(t+dt) */
      /* Volume acceleration */
      DP = Sqr(s) * (press_at - Oparams.P) / OprogStatus.W;
#endif
      B = Vol1i + 0.5 * dt * DP;
      Vol1 = B / (1.0 - 0.5 * dt * s1 / s); /* Vol1(t+dt) */
      Vol2 = DP + s1 * Vol1 / s;            /* Vol2(t+dt) */
      
      /* Calculate the Andersen forces */
      dlnV   = Vol1 / Vol; 
      dlnVSq = Sqr(dlnV);
      cfAnd  = - 2.0 * dlnVSq / 9.0;
      cfAnd += Vol2 / Vol / 3.0;
      cfAnd += dlns * dlnV / 3.0;
	 
      for(i=0; i < Nm; i++)
	{
	  for(a=0; a < NA; a++)
	    {
	      /* Calculate center of mass position for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	    
	      FxAnd  = cfAnd * Rx[i] * m[a];
	      FyAnd  = cfAnd * Ry[i] * m[a];
	      FzAnd  = cfAnd * Rz[i] * m[a];
	      Fx[a][i] = FxLJ[a][i] + FxAnd;
	      Fy[a][i] = FyLJ[a][i] + FyAnd;
	      Fz[a][i] = FzLJ[a][i] + FzAnd;
	      /* vt2 = v(t+dt/2) */
	      vxold[a][i] = vx[a][i];
	      vyold[a][i] = vy[a][i];
	      vzold[a][i] = vz[a][i];

	      vx[a][i] = vxt2[a][i] + dt2 * Fx[a][i] / m[a];
	      vy[a][i] = vyt2[a][i] + dt2 * Fy[a][i] / m[a];
	      vz[a][i] = vzt2[a][i] + dt2 * Fz[a][i] / m[a];
	      vx[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vy[a][i] /= 1.0 + 0.5 * dt * dlns;
	      vz[a][i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[a][i] * m[a];
	      FyNose = - dlns * vy[a][i] * m[a];
	      FzNose = - dlns * vz[a][i] * m[a];
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[a][i] = Fx[a][i] + FxNose;
	      Fy[a][i] = Fy[a][i] + FyNose;
	      Fz[a][i] = Fz[a][i] + FzNose;
	    } 
	}
      shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
      /* Zero the velocity along bond and calculate the constraint virial, this
	 should be as before because the Andersen Forces don't act along 
	 bond.*/
      numok=0;
      for (i=0; i < Oparams.parnum; i++)
	{
	  for (a = 0; a < NA; a++)
	    {
	      /* F = Flennard-jones + Fnose + Fandersen */
	      if (fabs(vx[a][i]- vxold[a][i]) < ittol &&
		  fabs(vy[a][i]- vyold[a][i]) < ittol &&
		  fabs(vz[a][i]- vzold[a][i]) < ittol)
		numok++;
	    }
	}
      if (numok == NA*Oparams.parnum)
	{
	  /* for all atoms diff between current velocities and prev vels is below tol! */ 
	  /*printf("Done %d iterations instead of %d\n", k, NUMIT);
	  */
	  break;
	}
    }

  /* Calculate kinetic energy */
  kinet(Nm, vx, vy, vz, Vol1);/* K(t+dt) */
  //printf("K exact: %.20f\n", K);

#ifdef MOLPTENS
  /* Calculate all components of molecular pressure tensor */
  calcPtensMol(Nm, Vol1);
  press_m = (Pmxx + Pmyy + Pmzz) / 3.0; /* press(t+dt) */
#endif

#if !defined(MOLPTENS) || defined(ATPTENS)
  /* Calculate all components of atomic pressure tensor */
  calcPtensAt(Nm, Vol1);
  press_at = (Patxx + Patyy + Patzz) / 3.0;
#endif 

#if !defined(MOLPTENS) && !defined(ATPTENS) && defined(ATPRESS)
  /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
  press_at = calcT1diagAt(Nm, Vol1) + (W + WC) / Vol;
  /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
#endif

#ifdef MOLPTENS
  /* Atomic pressure */
  press = press_m;
  Pxy = Pmxy;
  Pyz = Pmyz;
  Pzx = Pmzx;
#else
  press = press_at;
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;
#endif

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  Vol1o2 = Vol1o1;
  s1o2  = s1o1;
  Vol1o1 = Vol1t;
  s1o1 = s1t;
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }
  /* END OF ITERATION LOOP TO IMPROVE ANDERSEN-NOSE FORCE ESTIMATION */
}
#ifdef MD_BROWN_BETTER
/* ========================= >>> moveb <<< =========================== */
void moveb(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	   COORD_TYPE m[NA], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  int done;
  int moving[NA], moved[NA];
  COORD_TYPE Ka[NA], Mtot, vmx, vmy, vmz;
  COORD_TYPE rxab, ryab, rzab, rvab, gab, dSq;
  double vxp, vyp, vzp, Fxp, Fyp, Fzp;
  COORD_TYPE vxab, vyab, vzab;
  COORD_TYPE dx, dy, dz, dt2, rma, rmb, c2dtP, c2dtN, c2dt;
  COORD_TYPE rxi[NA], ryi[NA], rzi[NA], vxi[NA], vyi[NA], vzi[NA];
  double c0P, c1P, c2P, chsiP, c0, c1, chsi; 
  double c0N, c1N, c2N, chsiN;
  int i, a, b, it;

  
  /* ******************************************************************* */
  dt2 = dt*0.5;
  if (OprogStatus.Nose == -1)
    {
      chsiP = Oparams.chsiP;
      c0P = exp(-chsiP*dt);
      c1P = (1-c0P)/(chsiP*dt);
      c2P = (1-c1P)/(chsiP*dt);
      c2dtP = dt * c2P;
      chsiN = Oparams.chsiN;
      c0N = exp(-chsiP*dt);
      c1N = (1-c0N)/(chsiN*dt);
      c2N = (1-c1N)/(chsiN*dt);
      c2dtN = dt * c2N;
    }
  else
    {
      c2dt= dt / 2.0;
    }
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;
  Mtot = 0.0;
  for(a=0; a < NA; a++)
    {
      Mtot += m[a];
      Ka[a] = 0.0;
    }

  WC = 0.0;
#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = 0.0; /* Constraints-virial off-diagonal terms of pressure tensor */
  WCyz = 0.0;
  WCzx = 0.0;

  WCxx = 0.0; /* Constraints-virial off-diagonal terms of pressure tensor */
  WCyy = 0.0;
  WCzz = 0.0;
  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;
#endif

#ifdef MOLPTENS
  T1mxy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1myz = 0.0;
  T1mzx = 0.0;
  T1mxx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1myy = 0.0;
  T1mzz = 0.0;
#endif

  dSq = Sqr(d);
  /* LOOP OVER ALL MOLECULES */
  for(i=0; i < Nm; i++)
    {
      /* VELOCITY VERLET ALGORITHM PART B */
      for(a=0; a < NA; a++)
	{
	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];
	  if (OprogStatus.Nose == -1)
	    {
	      lab2body(u1[i], u2[i], u3[i], vx[a][i], vy[a][i], vz[a][i], &vxp, &vyp, &vzp);
	      lab2body(u1[i], u2[i], u3[i], Fx[a][i], Fy[a][i], Fz[a][i], &Fxp, &Fyp, &Fzp);
	      vxp = vxp + c2dtP * Fxp / m[a];
	      vyp = vyp + c2dtP * Fyp / m[a];
	      vzp = vzp + c2dtN * Fzp / m[a];
	      body2lab(u1[i], u2[i], u3[i], vxp, vyp, vzp, &vxi[a], &vyi[a], &vzi[a]);
	    }
	  else
	    {
	      vxi[a] = vx[a][i] + c2dt * Fx[a][i] / m[a];
	      vyi[a] = vy[a][i] + c2dt * Fy[a][i] / m[a];
	      vzi[a] = vz[a][i] + c2dt * Fz[a][i] / m[a];
	    }
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
	      if (b >= NA) b = 0;

	      if ( (moved[a] == 1) || (moved[b] == 1) ) 
		{
		  vxab = vxi[a] - vxi[b];
		  vyab = vyi[a] - vyi[b];
		  vzab = vzi[a] - vzi[b];
		  rxab = rxi[a] - rxi[b];
		  ryab = ryi[a] - ryi[b];
		  rzab = rzi[a] - rzi[b];
		  rxab = rxab - L*rint(invL*rxab);
		  ryab = ryab - L*rint(invL*ryab);
		  rzab = rzab - L*rint(invL*rzab);
		  rvab = rxab * vxab + ryab * vyab + rzab * vzab;
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

		      /* =================================== */
		      
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
#ifdef MOLPTENS
    vmx = vmy = vmz = 0.0;  /* Velocity of center of mass of molecule i */
#endif
    for(a=0; a < NA;  a++)
	{
	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
#ifdef MOLPTENS
	  vmx += vxi[a] * m[a];
	  vmy += vyi[a] * m[a];
	  vmz += vzi[a] * m[a];
#endif
#if !defined(MOLPTENS) || defined(ATPTENS)
	  /* Kinetic terms of the pressure-tensor */
	  T1xy += vx[a][i] * vy[a][i] * m[a]; 
	  T1yz += vy[a][i] * vz[a][i] * m[a];
	  T1zx += vz[a][i] * vx[a][i] * m[a];
	  T1xx += vx[a][i] * vx[a][i] * m[a]; 
	  T1yy += vy[a][i] * vy[a][i] * m[a];
	  T1zz += vz[a][i] * vz[a][i] * m[a];
#endif
	  Ka[a] = Ka[a] + Sqr(vxi[a]) + Sqr(vyi[a]) + Sqr(vzi[a]);
	}
#ifdef MOLPTENS
    /* Kinetic component of pressure tensor (all terms) */
    T1mxy += vmx * vmy / Mtot; 
    T1myz += vmy * vmz / Mtot;
    T1mzx += vmz * vmx / Mtot;
    T1mxx += vmx * vmx / Mtot;
    T1myy += vmy * vmy / Mtot;
    T1mzz += vmz * vmz / Mtot;
#endif    
    }
  /* END OF LOOP OVER MOLECULES */
  
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }

  for(a=0; a < NA; a++)
    {
      K  = K + Ka[a] * m[a] * 0.5;
    }
     
  WC = WC / dt2 / 3.0;

#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = WCxy / dt2; /* WCxy, ... are not exactly virial terms and the 3.0
			is not present */
  WCyz = WCyz / dt2;
  WCzx = WCzx / dt2;
#endif

  /* !!!!!##$$%%^^*/
#ifdef MOLPTENS
  /* Calculate molecular pressure */
  press_m =  T1mxx + T1myy + T1mzz + Wm;
  press_m /= Vol * 3.0;

#endif  
  
#if !defined(MOLPTENS) || defined(ATPTENS)
  /* Calculate the other element (off-diagonal) of the pressure tensor */
  Patxy = T1xy + Wxy + WCxy;
  Patyz = T1yz + Wyz + WCyz;
  Patzx = T1zx + Wzx + WCzx;  
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
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
  Pmxy = T1mxy + Wmxy;
  Pmyz = T1myz + Wmyz;
  Pmzx = T1mzx + Wmzx;  
  Pmxy /= Vol;
  Pmyz /= Vol;
  Pmzx /= Vol;
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
#else
/* ========================= >>> moveb <<< =========================== */
void moveb(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	   COORD_TYPE m[NA], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  int done;
  int moving[NA], moved[NA];
  COORD_TYPE Ka[NA], Mtot, vmx, vmy, vmz;
  COORD_TYPE rxab, ryab, rzab, rvab, gab, dSq;
  COORD_TYPE vxab, vyab, vzab;
  COORD_TYPE dx, dy, dz, dt2, rma, rmb, c2dt;
  COORD_TYPE rxi[NA], ryi[NA], rzi[NA], vxi[NA], vyi[NA], vzi[NA];
  double c0, c1, c2, chsi; 
  int i, a, b, it;

  
  /* ******************************************************************* */
  dt2 = dt*0.5;
  if (OprogStatus.Nose == -1)
    {
      chsi = Oparams.chsi;
      c0 = exp(-chsi*dt);
      c1 = (1-c0)/(chsi*dt);
      c2 = (1-c1)/(chsi*dt);
      c2dt = dt * c2;
    }
  else
    {
      c2dt= dt / 2.0;
    }
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;
  Mtot = 0.0;
  for(a=0; a < NA; a++)
    {
      Mtot += m[a];
      Ka[a] = 0.0;
    }

  WC = 0.0;
#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = 0.0; /* Constraints-virial off-diagonal terms of pressure tensor */
  WCyz = 0.0;
  WCzx = 0.0;

  WCxx = 0.0; /* Constraints-virial off-diagonal terms of pressure tensor */
  WCyy = 0.0;
  WCzz = 0.0;
  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;
#endif

#ifdef MOLPTENS
  T1mxy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1myz = 0.0;
  T1mzx = 0.0;
  T1mxx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1myy = 0.0;
  T1mzz = 0.0;
#endif

  dSq = Sqr(d);
  /* LOOP OVER ALL MOLECULES */
  for(i=0; i < Nm; i++)
    {
      /* VELOCITY VERLET ALGORITHM PART B */
      for(a=0; a < NA; a++)
	{
	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];
	  vxi[a] = vx[a][i] + c2dt * Fx[a][i] / m[a];
	  vyi[a] = vy[a][i] + c2dt * Fy[a][i] / m[a];
	  vzi[a] = vz[a][i] + c2dt * Fz[a][i] / m[a];
	
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
	      if (b >= NA) b = 0;

	      if ( (moved[a] == 1) || (moved[b] == 1) ) 
		{
		  vxab = vxi[a] - vxi[b];
		  vyab = vyi[a] - vyi[b];
		  vzab = vzi[a] - vzi[b];
		  rxab = rxi[a] - rxi[b];
		  ryab = ryi[a] - ryi[b];
		  rzab = rzi[a] - rzi[b];
		  rxab = rxab - L*rint(invL*rxab);
		  ryab = ryab - L*rint(invL*ryab);
		  rzab = rzab - L*rint(invL*rzab);
		  rvab = rxab * vxab + ryab * vyab + rzab * vzab;
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

		      /* =================================== */
		      
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
#ifdef MOLPTENS
    vmx = vmy = vmz = 0.0;  /* Velocity of center of mass of molecule i */
#endif
    for(a=0; a < NA;  a++)
	{
	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
#ifdef MOLPTENS
	  vmx += vxi[a] * m[a];
	  vmy += vyi[a] * m[a];
	  vmz += vzi[a] * m[a];
#endif
#if !defined(MOLPTENS) || defined(ATPTENS)
	  /* Kinetic terms of the pressure-tensor */
	  T1xy += vx[a][i] * vy[a][i] * m[a]; 
	  T1yz += vy[a][i] * vz[a][i] * m[a];
	  T1zx += vz[a][i] * vx[a][i] * m[a];
	  T1xx += vx[a][i] * vx[a][i] * m[a]; 
	  T1yy += vy[a][i] * vy[a][i] * m[a];
	  T1zz += vz[a][i] * vz[a][i] * m[a];
#endif
	  Ka[a] = Ka[a] + Sqr(vxi[a]) + Sqr(vyi[a]) + Sqr(vzi[a]);
	}
#ifdef MOLPTENS
    /* Kinetic component of pressure tensor (all terms) */
    T1mxy += vmx * vmy / Mtot; 
    T1myz += vmy * vmz / Mtot;
    T1mzx += vmz * vmx / Mtot;
    T1mxx += vmx * vmx / Mtot;
    T1myy += vmy * vmy / Mtot;
    T1mzz += vmz * vmz / Mtot;
#endif    
    }
  /* END OF LOOP OVER MOLECULES */
  
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }

  for(a=0; a < NA; a++)
    {
      K  = K + Ka[a] * m[a] * 0.5;
    }
     
  WC = WC / dt2 / 3.0;

#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = WCxy / dt2; /* WCxy, ... are not exactly virial terms and the 3.0
			is not present */
  WCyz = WCyz / dt2;
  WCzx = WCzx / dt2;
#endif

  /* !!!!!##$$%%^^*/
#ifdef MOLPTENS
  /* Calculate molecular pressure */
  press_m =  T1mxx + T1myy + T1mzz + Wm;
  press_m /= Vol * 3.0;

#endif  
  
#if !defined(MOLPTENS) || defined(ATPTENS)
  /* Calculate the other element (off-diagonal) of the pressure tensor */
  Patxy = T1xy + Wxy + WCxy;
  Patyz = T1yz + Wyz + WCyz;
  Patzx = T1zx + Wzx + WCzx;  
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
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
  Pmxy = T1mxy + Wmxy;
  Pmyz = T1myz + Wmyz;
  Pmzx = T1mzx + Wmzx;  
  Pmxy /= Vol;
  Pmyz /= Vol;
  Pmzx /= Vol;
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
#endif
/* DESCRIPTION:
   From Allen-Tildesley: 
   "When a molecule leaves the box by crossing one of the boundaries, 
   it is usual, to switch attention to the image molecule entering the box,
   by simply adding L, to or subtracting L from, the appropriate 
   coordinate." 
   
   NOTE: In the present case the box is cubic, but it is possible to 
   implement other geometries (see F.01)*/ 

/* ========================== >>> scalCor <<< ============================= */
void scalCor(int Nm)
{ 
  int i, a;
  COORD_TYPE cost, costo1, costo2, DRx, DRy, DRz, Rx, Ry, Rz;
  
  L = cbrt(Vol);
  invL = 1.0 / L;
  cost = Vol1 / Vol / 3.0;
  costo1 = Vol1o1 / Vol / 3.0;
  costo2 = Vol1o2 / Vol / 3.0;
  
  /* Reduced particles to first box */
  for(i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &Rx, &Ry, &Rz);
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * Rx);
      DRy = - L * rint(invL * Ry);
      DRz = - L * rint(invL * Rz);
      
      for(a=0; a < NA; a++)
	{
	  rx[a][i] += DRx;
	  ry[a][i] += DRy;
	  rz[a][i] += DRz;
	  
	  /* The velocity in the Andersen method depend upon the position,
	     so we must add a term also to velocity when switching attention
	     to the image particle. 
	     See Brown and Clarke Mol. Phys., Vol.51, No.5, 1243-1252, 1984 */ 
	  vx[a][i] += cost*DRx;
	  vy[a][i] += cost*DRy;
	  vz[a][i] += cost*DRz;
	  vxo1[a][i] += costo1*DRx;
	  vyo1[a][i] += costo1*DRy;
	  vzo1[a][i] += costo1*DRz;
	  vxo2[a][i] += costo2*DRx;
	  vyo2[a][i] += costo2*DRy;
	  vzo2[a][i] += costo2*DRz;
	}
    }
}

const COORD_TYPE bc1 = 14.0/45.0, bc2 = 64.0/45.0, bc3 = 24.0/45.0;
  
/* =========================== >>> BodeTerm <<< ============================*/
COORD_TYPE BodeTerm(COORD_TYPE dt, COORD_TYPE* fi)
{
  return dt * (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	       bc1 * fi[4]);
}
/* Routines for LU decomposition from Numerical Recipe online */
#define TINY 1E-20
void ludcmp(double a[3][3], int* indx, double* d, int *ok)
{
  /* A[i][j] = Aij 
   * A x = b  
   * per semplicità nel seguito si assume che l'ordine della matrice è 3 */
  int i,imax=-1,j,k;
  const int n = 3; 
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

void lubksb(double a[3][3], int* indx, double *b)
{ 
  int i,ii=0,ip,j; 
  double sum; 
  const int n = 3;
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
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);

/* ========================== >>> updateAng <<< ============================ */
void updateAng(int Nm)
{
  int ok, i, a, b, indx[3];
  double rcmx[3], rcmy[3], rcmz[3], norm, 
	 vcmx[3], vcmy[3], vcmz[3], RCMx, RCMy, RCMz, VCMx, VCMy, VCMz, dt;
  double L[3], rvx, rvy, rvz, Itens[3][3], d, u1[3], u2[3], u3[3], omCM[3],
	 omu3, norm1, norm2, u1u2, drCM[3];
#if 0 
  Itens = malloc(sizeof(double*)*3);
  for (a=0; a < 3; a++)
    Itens[a] = mall
#endif
  dt = Oparams.steplength;
  for(i=0; i < Oparams.parnum; i++)
    {
      L[0] = L[1] = L[2] = 0;
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  Itens[a][b] = 0.0;
      CoM(i, &RCMx, &RCMy, &RCMz);
      CoMV(i, &VCMx, &VCMy, &VCMz); 
      for (a=0; a < 3; a++)
	{
	  rcmx[a] = rx[a][i] - RCMx;
	  rcmy[a] = ry[a][i] - RCMy;
	  rcmz[a] = rz[a][i] - RCMz;
	  vcmx[a] = vx[a][i] - VCMx;
	  vcmy[a] = vy[a][i] - VCMy;
	  vcmz[a] = vz[a][i] - VCMz;
	  vectProd(rcmx[a], rcmy[a], rcmy[a], vcmx[a], vcmy[a], vcmz[a], &rvx, &rvy, &rvz);
	  L[0] += rvx * Oparams.m[a];
	  L[1] += rvy * Oparams.m[a];
	  L[2] += rvz * Oparams.m[a];

	  /* and now we calculate inertia tensor */
	  Itens[0][0] += Oparams.m[a]*(Sqr(rcmy[a]) + Sqr(rcmz[a]));   
	  Itens[1][1] += Oparams.m[a]*(Sqr(rcmz[a]) + Sqr(rcmx[a]));
	  Itens[2][2] += Oparams.m[a]*(Sqr(rcmx[a]) + Sqr(rcmy[a]));
	  Itens[0][1] += -Oparams.m[a]*rcmx[a]*rcmy[a];
	  Itens[0][2] += -Oparams.m[a]*rcmx[a]*rcmz[a];
	  Itens[1][2] += -Oparams.m[a]*rcmy[a]*rcmz[a];
	}
      //printf("L=%.15G %.15G %.15G\n", L[0], L[1], L[2]);
      Itens[1][0] = Itens[0][1];
      Itens[2][0] = Itens[0][2];
      Itens[2][1] = Itens[1][2];
      /* now we solve the linear system L = I omega using LU decomposition */
      /*gaussElimination(Itens, Lx, Ly, Lz);*/
      ludcmp(Itens,indx, &d, &ok);
      /* forward and backward substitution */
      lubksb(Itens, indx, L);
      ox[i] = L[0];
      oy[i] = L[1];
      oz[i] = L[2];
#if 0
      oz[a] = Lz / Itens[2][2];
      oy[a] = (Ly - Itens[1][2] * oz[a]) / Itens[1][1];
      ox[a] = (Lx - Itens[0][1] * oy[a] - Itens[0][2] * oz[a] ) / Itens[0][0]; 
#endif
      /*
      detI = Itens[0][0] * (Itens[1][1]*Itens[2][2]-Itens[1][2]*Itens[2][1]) -
       Itens[0][1] * (Itens[1][0]*Itens[2][2]-Itens[1][2]*Itens[2][1]) +
       Itens[0][2] * (Itens[1][0]*Itens[2][1] - Itens[1][1]*Itens[2][0]);
       */
      u1[0] = rx[1][i] - rx[0][i];
      u1[1] = ry[1][i] - ry[0][i];
      u1[2] = rz[1][i] - rz[0][i];
      u2[0] = rx[2][i] - rx[0][i];
      u2[1] = ry[2][i] - ry[0][i];
      u2[2] = rz[2][i] - rz[0][i];
      norm2 = calc_norm(u2);
      for (a=0;  a < 3; a++)
	u2[a] /= norm2;
      u1u2 = u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
      for (a=0;  a < 3; a++)
	{
	  u1[a] = u1[a] - u1u2*u2[a];
	}
      norm1 = calc_norm(u1);
      for (a=0;  a < 3; a++)
	{
	  u1[a] /= norm1;
	}
      vectProd(u1[0],u1[1],u1[2],u2[0],u2[1], u2[2],&u3[0], &u3[1], &u3[2]); 
#if 0
      printf("u1=%.15G %.15G %.15G\n", u1[0], u1[1], u1[2]);
      printf("u2=%.15G %.15G %.15G\n", u2[0], u2[1], u2[2]);
      printf("u3=%.15G %.15G %.15G\n", u3[0], u3[1], u3[2]);
      printf("normu1=%.15G normu2=%.15G normu3=%.15G\n", calc_norm(u1), 
	     calc_norm(u2), calc_norm(u3));
      printf("u1*u2=%.15G\n",u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2]);
#endif
 
      omu3 = u3[0]*ox[i] + u3[1]*oy[i] + u3[2]*oz[i];
      omCM[0] = u1[0]*ox[i] + u1[1]*oy[i] + u1[2]*oz[i];
      omCM[1] = u2[0]*ox[i] + u2[1]*oy[i] + u2[2]*oz[i];
      omCM[2] = u3[0]*ox[i] + u3[1]*oy[i] + u3[2]*oz[i];
      drCM[0] = u1[0]*VCMx + u1[1]*VCMy + u1[2]*VCMz;
      drCM[1] = u2[0]*VCMx + u2[1]*VCMy + u2[2]*VCMz;
      drCM[2] = u3[0]*VCMx + u3[1]*VCMy + u3[2]*VCMz;
      OprogStatus.sumox[i] += omCM[0] * dt; 
      OprogStatus.sumoy[i] += omCM[1] * dt;
      OprogStatus.sumoz[i] += omCM[2] * dt;
      
      OprogStatus.sumdrx[i] += drCM[0] * dt; 
      OprogStatus.sumdry[i] += drCM[1] * dt;
      OprogStatus.sumdrz[i] += drCM[2] * dt;

      /* Time integral of ang. vel. stored in xva variables */
      Dphix[i] = OprogStatus.sumox[i];
      Dphiy[i] = OprogStatus.sumoy[i];
      Dphiz[i] = OprogStatus.sumoz[i];
    }
}
/* ============================ >>> updateQ <<< =========================== */
void updateDQ(COORD_TYPE dt)
{
  /* Intgerate the pressure tensor using Bode's rule
     WARNING: the minimum effective update  for DQ(t) is 4 steps, because
     Bode rule needs 5 points at least */ 
  int curStep = Oparams.curStep, i1;
  COORD_TYPE *PxyArr = OprogStatus.PxyArr,
    *PyzArr = OprogStatus.PyzArr,
    *PzxArr = OprogStatus.PzxArr;
  
  //printf(" curStep: %d\n", curStep);
  /* Store points to use later for integration */
  
  if (curStep == 1)
    {
      /* WARNING: we assume here that the simulation starts from 1, 
	 that is the first step must be 1 <<<=================== !!!!!!!!!! */
      /* First time only HERE */
      PxyArr[4] = Pxy;
      PyzArr[4] = Pyz;
      PzxArr[4] = Pzx;
      return;
    }
 
  i1 = (curStep - 1) % 4;
  
  if (i1 != 0)
    {
      PxyArr[i1] = Pxy;
      PyzArr[i1] = Pyz;
      PzxArr[i1] = Pzx;
    }

  if (i1 == 0)
    {
      /* Last point of previuos set is the first point of the actual set */
      PxyArr[0] = PxyArr[4];
      PyzArr[0] = PyzArr[4];
      PzxArr[0] = PzxArr[4];
      /* 5th point is the actual value of pressure tensor */
      PxyArr[4] = Pxy;
      PyzArr[4] = Pyz;
      PzxArr[4] = Pzx;
      /* Here we use molecular pressure tensor */
      OprogStatus.DQxy += BodeTerm(dt, PxyArr);
      OprogStatus.DQyz += BodeTerm(dt, PyzArr);
      OprogStatus.DQzx += BodeTerm(dt, PzxArr);
    }

  /*         _ t
	    |
    DQab =  |  Pab dt
            |
	   - 0
	   
    con a, b = x,y,z 
  */
}


/* ============================ >>> chkVol <<< ============================= */
void chksVol(void)
{
  /* DESCRIPTION:
     This subroutine check if the actual volume is near the volume
     OprogStatus.avVol, if it is the case the program exit 
     This is useful if we want to perform a micronical production run 
     after the equilibration one (NPT) */
  COORD_TYPE DVol, rels, relVol, Ds;
  COORD_TYPE tolVol;
  COORD_TYPE tols; 
  
  tolVol  = OprogStatus.tolVol;
  tols    = OprogStatus.tols;
  if ( (OprogStatus.avVol > 0.0) && (OprogStatus.avs > 0.0))
    {
      DVol = fabs(Vol - OprogStatus.avVol);
      relVol = DVol / Vol;
      Ds   = fabs(s - OprogStatus.avs);
      rels   = Ds / s;
      if ( (relVol < tolVol) && (rels < tols) )
	{
	  /* This cancel the status file so the program doesn't try
	     to restart automatically */ 
	  printf("Volume and s within tolerance, exit...\n");
	  // printf("avVol1: %.10f relVol1: %.10f tolVol1: %.10f\n", OprogStatus.avVol1, relVol1, tolVol1);
	  // printf("Vol1: %f DVol1: %f\n", Vol1, DVol1);
	  ENDSIM = 1; /* End simulation */
	}
      
    }
}

void chks(void)
{
  /* DESCRIPTION:
     This subroutine check if the actual volume is near the volume
     OprogStatus.avVol, if it is the case the program exit 
     This is useful if we want to perform a micronical production run 
     after the equilibration one (NPT) */
  COORD_TYPE rels, Ds;
  COORD_TYPE tols, relT, DT, temp; 
  
  tols    = OprogStatus.tols;
  if (OprogStatus.avs > 0.0)
    {
      temp = 2.0 * K / (5.0 * Oparams.parnum - 3.0);
      DT = fabs(temp - Oparams.T);
      Ds   = fabs(s - OprogStatus.avs);
      relT = DT / Oparams.T;
      rels   = Ds / s;
      /* QUESTION: We must check also T or s1 or it is enough to check
         s? (actually it check only s) */
      if (rels < tols)// && (relT < tols) )
	{
	  /* This cancel the status file so the program doesn't try
	     to restart automatically */ 
	  printf("s within tolerance, exit...\n");
	  // printf("avVol1: %.10f relVol1: %.10f tolVol1: %.10f\n", OprogStatus.avVol1, relVol1, tolVol1);
	  // printf("Vol1: %f DVol1: %f\n", Vol1, DVol1);
	  ENDSIM = 1; /* End simulation */
	}
      
    }
}
#define HALFSQRT3 0.86602540378443864676 
static double HP[6][2] = {{-0.5,HALFSQRT3},{0.5,HALFSQRT3},{1.0,0.0},
    {0.5,-HALFSQRT3},{-0.5,-HALFSQRT3},{-1.0,0.0}};
void buildhexagon(double diam, int n, double *rx, double *ry, int *np)
{
  /* costruisce un esagono di diametro diam mettendo n punti per lato.
   * Le coordinate dei punti vengono salvate nei vettori rx e ry che 
   * quindi devono "ospitare" 6*(n+1) punti */
  int nvert, nint, nextvert;
  double displx, disply; 
  double halfdiam;
  /*FILE *f;*/
  
  halfdiam = 0.5*diam;
  /*f=fopen("coordhex.a","a+");*/
  for (nvert = 0; nvert < 6; nvert++)
    {
      nextvert = nvert+1;
      if (nextvert == 6) 
	nextvert = 0;
      displx = HP[nextvert][0] - HP[nvert][0];
      disply = HP[nextvert][1] - HP[nvert][1];
      for (nint = 0; nint <= n; nint++)
	{
	  rx[*np] = halfdiam*(HP[nvert][0] + ((double) nint)*displx / ((double)n+1));
	  ry[*np] = halfdiam*(HP[nvert][1] + ((double) nint)*disply / ((double)n+1));
	  /*fprintf(f,"%f %f\n", rx[*np], ry[*np]);*/
	  (*np)++;
	}
      
    }
  /*fclose(f);*/
}
void savesnap(void)
{
  int i, a;
  FILE *f;
  char fileop2[1024], fileop[1024];
  double L, rrx, rry, rrz, RCMx, RCMy, RCMz, nx, ny, nz;
  double d1x, d1y, d1z, d2x, d2y, d2z;
#ifdef MDLLINT
  sprintf(fileop2 ,"SnaT%.6G_%s_%lld", 
	  Oparams.T, 
	  OprogStatus.nRun, Oparams.curStep);
#else
  sprintf(fileop2 ,"SnaT%.6G_%s_%d", 
	  Oparams.T, 
	  OprogStatus.nRun, Oparams.curStep);
#endif
  strcpy(fileop, absTmpAsciiHD(fileop2));
  
  if ( (f = fopen(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  L = cbrt(Vol);
  fprintf(f, ".Vol: %f\n", Vol);
  for (i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &RCMx, &RCMy, &RCMz);
      if (OprogStatus.snapmode==1)
	{
	  RCMx = RCMx -  L * rint( RCMx / L);
	  RCMy = RCMy -  L * rint( RCMy / L);
	  RCMz = RCMz -  L * rint( RCMz / L);
	  d1x = rx[0][i] - rx[1][i];
	  d1y = ry[0][i] - ry[1][i];
	  d1z = rz[0][i] - rz[1][i];
	  d2x = rx[1][i] - rx[2][i];
	  d2y = ry[1][i] - ry[2][i];
	  d2z = rz[1][i] - rz[2][i];
	  vectProd(d1x, d1y, d1z, d2x, d2y, d2z, &nx, &ny, &nz);
	  fprintf(f, "%.15G %.15G %.15G %.15G %.15G %.15G @ %f %f\n", 
		  RCMx, RCMy, RCMz, nx, ny, nz, Oparams.Diam/2.0, 1.0); 
	}
      else
	{
	  for (a=0; a < Oparams.nsites; a++)
	    {
	      rrx = rallx[a][i] - L * rint( RCMx / L); 
	      rry = rally[a][i] - L * rint( RCMy / L);
	      rrz = rallz[a][i] - L * rint( RCMz / L);
	      fprintf(f, "%.15G %.15G %.15G\n", rrx, rry, rrz); 
	    }
	}
    }
  fclose(f);
}
void buildAtomsPositions(void)
{
  int i, a, n, spe, np=0, ip;
  double RCMx, RCMy, RCMz;
  double halfD, invhalfD, quartD, invquartD, invsRoot3, invsr3D, sr3D;
  double sRoot3, tmpx, tmpy, e1x, e1y, e1z, e2x, e2y, e2z, racmx, racmy;
  double *rhexagx, *rhexagy;
  FILE* f;
  sRoot3 = sqrt(3);
  halfD = Oparams.Diam * 0.5;
  quartD = Oparams.Diam *0.25;
  invhalfD = 1.0 / halfD;
  invquartD = 1.0 / quartD;
  invsRoot3 = 1.0 / sRoot3;  
  invsr3D = invsRoot3 * invquartD;
  sr3D = quartD * sRoot3;
#if 0
  if (Oparams.curStep==1)
    f=fopen("lapo.pos","w");
#endif
  rhexagx = malloc(sizeof(double)*Oparams.nsites);
  rhexagy = malloc(sizeof(double)*Oparams.nsites);
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      /* centro di massa */
      CoM(i, &RCMx, &RCMy, &RCMz);
 
      /* NOTA: e1 ed e2 sono scazzati! */
      /* costruisco i vettore di base unitari*/
      e1x = -(rx[2][i] - RCMx) * invhalfD;
      e1y = -(ry[2][i] - RCMy) * invhalfD;
      e1z = -(rz[2][i] - RCMz) * invhalfD;
      /* notare che si sottrae la proiezione di r[i] - rCM lungo e1 */
      e2x = (rx[0][i] - RCMx - quartD * e1x) * invsr3D;
      e2y = (ry[0][i] - RCMy - quartD * e1y) * invsr3D;
      e2z = (rz[0][i] - RCMz - quartD * e1z) * invsr3D;

      /* I primi tre atomi sono quelli dotati di massa */
            /* quindi costruisco tutti gli altri esagoni */
      np=0;
      if (Oparams.invsp >= 2)
	for (n=(Oparams.invsp/2)-1; n >= 0; n--)
	  buildhexagon(2*(n+1)*(Oparams.Diam / Oparams.invsp), n, rhexagx, rhexagy, &np);
      /* l'esagono e' stato costruito nel sistema di riferimento individuato
       * da e1 e e2 ora bisogna passare alle coordinate della scatola */
      spe = Oparams.invsp / 2;
      
      for (a = 0; a < 3; a++)
	{
	  /* i primi 3 atomi sono quelli "di base" */
#if 0
	  rallx[i][a] = rallx[i][1 + a*spe];
	  rally[i][a] = rally[i][1 + a*spe];
#endif
	  ip = (a*2+1)*spe;
	  tmpx = rhexagx[ip];
	  tmpy = rhexagy[ip];
#if 0
	  if (Oparams.curStep==1)
  	    fprintf(f, "%.15f %.15f %.15f\n", rhexagx[1+a*spe], rhexagy[1+a*spe], 0.0);
#endif
	  rhexagx[ip] = rhexagx[a];
	  rhexagy[ip] = rhexagy[a];
	  rhexagx[a] = tmpx;
	  rhexagy[a] = tmpy;
	  /*fprintf(f, "%.15f %.15f %.15f\n", rallx[a][i], rally[a][i], rallz[a][i]);*/
	}
      rallx[Oparams.nsites-1][i] = RCMx;
      rally[Oparams.nsites-1][i] = RCMy;
      rallz[Oparams.nsites-1][i] = RCMz;
#if 0
      if (Oparams.curStep==1)
	fprintf(f,"%.15f %.15f %.15f\n", 0.0, 0.0, 0.0);
#endif
#if 0
      for (a = 0; a < 3; a++)
	{
	  fprintf(f,"%.15f %.15f %.15f @ 2.0\n", rx[a][i], ry[a][i], rz[a][i]);
	}
#endif
      for (a = 0; a < Oparams.nsites-1; a++)
	{
#if 0
	  racmx = rallx[a][i];
	  racmy = rally[a][i]; 
#endif
	  rallx[a][i] = RCMx + e1x*rhexagx[a] + e2x*rhexagy[a];
	  rally[a][i] = RCMy + e1y*rhexagx[a] + e2y*rhexagy[a];
	  rallz[a][i] = RCMz + e1z*rhexagx[a] + e2z*rhexagy[a];
#if 0
	  if (Oparams.curStep==1)
	    {
	      if (a < 3)
		fprintf(f, "%.15f %.15f %.15f C[blue]\n", rhexagx[a], rhexagy[a], 0.0);
	      else if (a < 6*(Oparams.invsp / 2)) 
		fprintf(f, "%.15f %.15f %.15f C[green]\n", rhexagx[a], rhexagy[a], 0.0);
	      else
		fprintf(f, "%.15f %.15f %.15f\n", rhexagx[a], rhexagy[a], 0.0);
	    }
#endif
	  /*fprintf(f, "%.15f %.15f %.15f\n", rhexagx[a], rhexagy[a], 0.0);*/
	}
    }
#if 0
  if (Oparams.curStep==1)
    fclose(f);
#endif
  free(rhexagx);
  free(rhexagy);
}
extern int **ncut;
#ifdef MD_TSHAPED
#undef MD_HOC
void move(void)
{
  /* DESCRIPTION:
     Move the particles by one step */
  double distance;
  int i; 
  FILE* f;
  double dist;
  f = fopen("V.dat", "w");
  for (i=1; i < 3000; i++)
    {
      nebrNow=1;
      dist =0.05 + 0.005 * i;
#ifdef MD_HOC
      /* T-SHAPED */
      FCC2(Oparams.parnum, Oparams.Diam, Oparams.m,dist-1.67468);
#else
      FCC2(Oparams.parnum, Oparams.Diam*0.5, Oparams.m,dist-2*1.67468);
#endif
      /* distanza fra i 3 atomi le cui coordinate evolvono nel tempo
       * Notare che i 3 atomi formano un triangolo equiliatero quindi
       * le 3 distanze sono uguali */
      distance = sqrt(3)*0.5*Oparams.Diam;
      /* calc predicted coords*/
      /* -1 = brownian dynamics NTV 
       * -2 = brownian dynamics at fixed pressure NTP */
      buildAtomsPositions();
      if (nebrNow)
	{
	  nebrNow = 0;
	  dispHi = 0.0;
	  /* build up linked list on predicted 
	     coordinates (only father do it)*/
	  if (OprogStatus.noLinkedList)
	    {
	      BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut);
	    }
	  else
	    {
	      links(Oparams.parnum, Oparams.rcut);
	      /* Build up neighbour list */  
	      BuildNebrList(Oparams.parnum, Oparams.rcut);
	    }
	}

      LJForce(Oparams.parnum, Oparams.rcut);

      /* considera tutti i contributi alle forza agente sugli atomi "di base"
       * ossia somma anche le forze dovute agli atomi senza massa moltiplicate
       * per gli opportuni coefficienti "vincolari" 
       * */
      ForceOn123();
      checkNebrRebuild();
      if (  ( OprogStatus.snapSteps < 0 && (abs(OprogStatus.snapSteps) == Oparams.curStep) ) || 
	    ( OprogStatus.snapSteps > 0 && (Oparams.curStep % OprogStatus.snapSteps == 0) )  )
	savesnap();
      /* Update accumulators for calculating the angular diffusion coefficent */
      fprintf(f, "%f %f\n", dist, V);
    }
  printf("boh dist 0-1: %f\n", sqrt(Sqr(rallx[4][0]-rallx[5][0])+
			Sqr(rally[4][0]-rally[5][0])+Sqr(rallz[4][0]-rallz[5][0])) );
  fclose(f);
}
#elif defined(MD_EFFPOT)
void move(void)
{
  /* DESCRIPTION:
     Move the particles by one step */
  double distance, sumpot;
  int i, phiI, thetaI, phi2I, theta2I, chsi2I; 
  FILE* f;
  double PI;
  const int nr = 20;
  const double V0 = 0.0;
  double rend, rini, dr, dist, phi, theta, phi2, theta2, chsi2;
  double dtheta, dphi, dtheta2, dphi2, dchsi2;
  int nthetaI = 5, nphiI = 5, nphi2I = 10, ntheta2I = 10, nchsi2I = 10;
  f = fopen("Veff.dat", "w");
  PI = 2.0*acos(0);
  dphi = 2.0*PI/((double)nphiI);
  dtheta = PI/((double)nthetaI);
  dphi2 = 2.0*PI/((double)nphi2I);
  dtheta2 = PI/((double)ntheta2I);
  dchsi2 = 2.0*PI/((double)nchsi2I);
  /* build up linked list on predicted 
     coordinates (only father do it)*/
  if (OprogStatus.noLinkedList)
    {
      BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut);
    }
  else
    {
      links(Oparams.parnum, Oparams.rcut);
      /* Build up neighbour list */  
      BuildNebrList(Oparams.parnum, Oparams.rcut);
    }
  rini = Oparams.Diam - 2.0;
  rend = (Oparams.rcut + 5.0) * Oparams.lambdaD;
  dr = (rend - rini) / nr; 
  for (i=1; i < nr; i++)
    {
      dist = rini + i * dr;
      /* loop over all angles we're dealing with */
      sumpot = 0.0;
      printf("dist=%.15G\n", dist);
      for (phiI=0; phiI < nphiI; phiI++)
	{
	  phi = phiI * dphi;  
	  for (thetaI=0; thetaI < nthetaI; thetaI++)
	    {
	      printf("phiI=%d thetaI=%d\n", phiI, thetaI);
	      theta = thetaI * dtheta;
	      for (phi2I=0; phi2I < nphi2I; phi2I++)
		{
		  phi2 = phi2I * dphi2;
		  for (theta2I=0; theta2I < ntheta2I; theta2I++)
		    {
		      theta2 = theta2I * dtheta2;
		      for (chsi2I=0; chsi2I < nchsi2I; chsi2I++)
			{
			  chsi2 =  chsi2I * dchsi2;
			  setpos(Oparams.parnum, Oparams.Diam, Oparams.m, dist, phi, theta, 
				 phi2, theta2, chsi2);
			  /* distanza fra i 3 atomi le cui coordinate evolvono nel tempo
			   * Notare che i 3 atomi formano un triangolo equiliatero quindi
			   * le 3 distanze sono uguali */
			  distance = sqrt(3)*0.5*Oparams.Diam;
			  /* calc predicted coords*/
			  /* -1 = brownian dynamics NTV 
			   * -2 = brownian dynamics at fixed pressure NTP */
			  buildAtomsPositions();
 			  LJForce(Oparams.parnum, Oparams.rcut);
			  /* considera tutti i contributi alle forza agente sugli atomi "di base"
			   * ossia somma anche le forze dovute agli atomi senza massa moltiplicate
			   * per gli opportuni coefficienti "vincolari" 
			   * */

			  //printf("Vc:%.15G T:%.15G\n", Vc, Oparams.T);
			  sumpot += exp(-(Vc-V0)/Oparams.T)*sin(theta2)*sin(theta)*
			    dtheta*dphi*dtheta2*dphi2*dchsi2; 
			  //ForceOn123();
			  //checkNebrRebuild();
			
			}
		    }
		}
	    }
	}
      
      /* Update accumulators for calculating the angular diffusion coefficent */
      fprintf(f, "%.15G %.15G\n", dist, V0-log(sumpot));
    }
  /* printf("boh dist 0-1: %f\n", sqrt(Sqr(rallx[4][0]-rallx[5][0])+
     Sqr(rally[4][0]-rally[5][0])+Sqr(rallz[4][0]-rallz[5][0])) );
     */
  fclose(f);
  ENDSIM=1;
}
#else
/* ============================ >>> move<<< =================================*/
void move(void)
{
  /* DESCRIPTION:
     Move the particles by one step */
  double distance;
  /* distanza fra i 3 atomi le cui coordinate evolvono nel tempo
   * Notare che i 3 atomi formano un triangolo equiliatero quindi
   * le 3 distanze sono uguali */
  distance = sqrt(3)*0.5*Oparams.Diam;
  /* calc predicted coords*/
  /* -1 = brownian dynamics NTV 
   * -2 = brownian dynamics at fixed pressure NTP */
  if (OprogStatus.Nose < 0)
    {
      movea_Brownian(Oparams.steplength, 0.0000000001, 30, 3, distance, Oparams.m, 
	    Oparams.parnum);  
    }
  else
    {
      movea(Oparams.steplength, 0.000000000001, 30, 3, distance, Oparams.m, 
	    Oparams.parnum);        
    }
  buildAtomsPositions();

  if (nebrNow)
    {
      nebrNow = 0;
      dispHi = 0.0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      if (OprogStatus.noLinkedList)
	{
	  BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut);
	}
      else
	{
	  links(Oparams.parnum, Oparams.rcut);
	  /* Build up neighbour list */  
	  BuildNebrList(Oparams.parnum, Oparams.rcut);
	}
    }
    
  if (OprogStatus.switchForcesOff==0)
    {
      LJForce(Oparams.parnum, Oparams.rcut);
      
      /* considera tutti i contributi alle forza agente sugli atomi "di base"
       * ossia somma anche le forze dovute agli atomi senza massa moltiplicate
       * per gli opportuni coefficienti "vincolari" 
       * */
      ForceOn123();
    }
  kinet(Oparams.parnum, vx, vy, vz, Vol1);

  /* correct the coords */
  if (OprogStatus.Nose == 1)
    {  
      /* NPT ensemble */
      movebNPT(Oparams.steplength, 0.0000000001, 30, 3, Oparams.m, distance, 
	       Oparams.parnum);             
    }
  else if (OprogStatus.Nose == 2)
    {
      /* NPT ensemble */
      /*movebNTV(Oparams.steplength, 0.0000000001, 30, 3, Oparams.m, distance, 
	       Oparams.parnum); */            
      movebNTV(Oparams.steplength, 0.000000000001, 30, 3, Oparams.m, distance, 
	       Oparams.parnum);  }
  else if (OprogStatus.Nose == -2)
    {
      movebBrownAnd(Oparams.steplength, 0.0000000001, 30, 3, Oparams.m, distance, 
		    Oparams.parnum); 
    }
  else
    /* 0 = NVE, -1 = dinamica browniana */ 
    {
      /* NVE ensemble o Dinamica Browniana */
      moveb(Oparams.steplength, 0.0000000001, 30, 3, Oparams.m, distance, 
	    Oparams.parnum);             
    }
    
  /* Calculate the kinetic energy */
  kinet(Oparams.parnum, vx, vy, vz, Vol1);

  checkNebrRebuild();
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
  
    if (  ( OprogStatus.snapSteps < 0 && (abs(OprogStatus.snapSteps) == Oparams.curStep) ) || 
	  ( OprogStatus.snapSteps > 0 && (Oparams.curStep % OprogStatus.snapSteps == 0) )  )
      savesnap();
  /* Update accumulators for calculating the angular diffusion coefficent */
  updateAng(Oparams.parnum);  
  /* Update the integral of the pressure tensor */
  updateDQ(Oparams.steplength);
  if ( (OprogStatus.CMreset != 0) &&
      /* 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) */
      (Oparams.curStep == OprogStatus.CMreset) ) 
      resetCM(Oparams.parnum);
#if 0
  check_distances();
#endif
}
#endif
