#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

extern volatile int *sc_;           /* shared counter */
extern volatile int *semA_, *semB_; /* semaphores */

extern int process;          /* 0 = child 1 = father ( see mdsimul.c )*/ 

extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern void  links(int Nm, COORD_TYPE rcut, COORD_TYPE sigab[NA][NA]);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]);
extern void BuildNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]);
extern void CoMB(int a, int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
extern void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
extern void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm);
extern void checkNebrRebuild(void);
extern void LJForce(int Nm, COORD_TYPE epsab[NA][NA], 
		    COORD_TYPE sigab[NA][NA], COORD_TYPE rcut);
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
COORD_TYPE Vc, V, W, K, WC, T1xx, T1yy, T1zz,
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

COORD_TYPE *FxL[NA], *FyL[NA], *FzL[NA]; /* Local arrays of forces
					    (used by the LJForce routine 
					    to avoid interferences) */
COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
COORD_TYPE  *Rmx, *Rmy, *Rmz;
/* ========================================================================= */

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

  tol2   = 2.0 * tol;
  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;
  dSq = Sqr(d); /* In general you must supply a vector of bond lengths */

  /* ===== LOOP OVER MOLECULES ===== */
  loopShrC(i, 1, Nm)
    {
      /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
      loop(a, 1, NA)
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

	  loop(a, 1, NB)
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
	       
	  loop(a, 1, NA)
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

      loop(a, 1, NA)
	{
	  rx[a][i] = pxi[a];
	  ry[a][i] = pyi[a];
	  rz[a][i] = pzi[a];

	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
	}
      
    }
  poolShrC;
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
  loopShrC(i, 1, Nm)
    {
      /* VELOCITY VERLET ALGORITHM PART B */
      loop(a, 1, NA)
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
      // *  
	  it = 0;
	  done = 0;
	  while ((done == 0 ) && ( it <= maxIt ))
	    {
	      done = 1;
	      loop(a, 1, NB)
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
	      loop(a, 1, NA)
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
	  
	  loop(a, 1, NA)
	    {
	      v2sx[a][i] = vxi[a];
	      v2sy[a][i] = vyi[a];
	      v2sz[a][i] = vzi[a];
	    }
    }
  poolShrC;
  /* END OF LOOP OVER MOLECULES */

  ProcSync0();
  WC = WC / dt2 / 3.0;
#ifndef NO_PARALLEL_CODE  
  sumAllProc(FATHER, &WC);
#endif

#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = WCxy / dt2; /* WCxy, ... are not exactly virial terms and the 3.0
			is not present */
  WCyz = WCyz / dt2;
  WCzx = WCzx / dt2;

  WCxx = WCxx / dt2;
  WCyy = WCyy / dt2;
  WCzz = WCzz / dt2;
#ifndef NO_PARALLEL_CODE
  sumAllProcMA(&WCxy, &WCyz, &WCzx,
	       &WCxx, &WCyy, &WCzz, NULL);
#endif
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
  loopShrC(i, 1, Nm)
    {
      CoMV(i, &Vx, &Vy, &Vz);
      kin +=  Sqr(Vx - dlnV * Rx[i]) + Sqr(Vy - dlnV * Ry[i]) + 
	Sqr(Vz - dlnV * Rz[i]);
    }
  poolShrC;
#ifndef NO_PARALLEL_CODE
  sumAllProcMA(&kin, NULL);
#endif
  
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
  loopShrC(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  kin +=  m[a]*(Sqr(vx[a][i] - dlnV * Rx[i]) + 
	    Sqr(vy[a][i] - dlnV * Ry[i]) + 
	    Sqr(vz[a][i] - dlnV * Rz[i]));
	}
    }
  poolShrC;
#ifndef NO_PARALLEL_CODE
  sumAllProcMA(&kin, NULL);
#endif
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
  loopShrC(i, 1, Nm)
    {
      vmx = vmy = vmz = 0.0;  /* Velocity of center of mass of molecule i */
      loop(a, 1, NA)
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
  poolShrC;
#ifndef NO_PARALLEL_CODE  
  sumAllProcMA(&T1mxy, &T1myz, &T1mzx, 
	       &T1mxx, &T1myy, &T1mzz, NULL);
#endif
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
  loopShrC(i, 1, Nm)
    {
      loop(a, 1, NA)
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
  poolShrC;
#ifndef NO_PARALLEL_CODE  
  sumAllProcMA(&T1xx, &T1yy, &T1zz, 
	       &T1xy, &T1yz, &T1zx, NULL);
#endif
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

/* ========================= >>> moveb <<< =========================== */
void movebNTV(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	      COORD_TYPE m[NA], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, dlns; 
  COORD_TYPE s1i;
  COORD_TYPE DT, A, dt2; 
  int i, a, k;
  const int NUMIT = 8;

  /* ******************************************************************* */
  dt2 = dt * 0.5;
  
  loopShrC(i, 1, Nm)
    {
      CoM(i, &Rx[i], &Ry[i], &Rz[i]);
      loop(a, 1, NA)
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
  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm, vx, vy, vz, 0.0);  /* K(t+dt) */
      //printf("K guess: %.20f\n", K);

      DT = s * (2.0 * K - (5.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      
      loopShrC(i, 1, Nm)
	{
	  loop(a, 1, NA)
	    {
	      /* Calculate center of mass posiotin for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	      Fx[a][i] = FxLJ[a][i];
	      Fy[a][i] = FyLJ[a][i];
	      Fz[a][i] = FzLJ[a][i];
	     
	      /* vt2 = v(t+dt/2) */
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
	      Fx[a][i] = FxLJ[a][i] + FxNose;
	      Fy[a][i] = FyLJ[a][i] + FyNose;
	      Fz[a][i] = FzLJ[a][i] + FzNose;
	    } 
	}
      shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
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
  loopShrC(i, 1, Nm)
    {
      loop(a, 1, NA)
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
  int i, a, k;
  const int NUMIT = 8;

  /* ******************************************************************* */
  dt2 = dt / 2.0;
  
  loopShrC(i, 1, Nm)
    {
      CoM(i, &Rx[i], &Ry[i], &Rz[i]);
      loop(a, 1, NA)
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
  Vol1i = Vol1;/* Vol1i = Vol1(t+dt/2)*/

  /* Initial guess for Vol1 */
  //Vol1 = 13.0*Vol1t/6.0 - 4.0*Vol1o1/3.0 + Vol1o2 / 6.0;/* Vol1(t+dt) */
  
  //Vol1 = 2.0 * Vol1t - Vol1o1; /* Verlet */  
  
  Vol1 = 5.0 * Vol1t / 2.0 - 2.0 * Vol1o1 + Vol1o2 / 2.0;  

  shakeVel(Nm, dt, m, maxIt, NB, d, tol, vx, vy, vz);
  //printf("Vol1: %f Vol1o1: %f Vol1o2: %f Vol1t:%f\n", 
	// Vol1, Vol1o1, Vol1o2, Vol1t);
  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm, vx, vy, vz, Vol1);  /* K(t+dt) */
      //printf("K guess: %.20f\n", K);

      DT = s * (2.0 * K - (5.0 * Nm - 3.0) * Oparams.T) / 
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
      //      printf("k: %d press_at: %f press_m: %f\n", k, press_at, press_m);
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
	 
      loopShrC(i, 1, Nm)
	{
	  loop(a, 1, NA)
	    {
	      /* Calculate center of mass posiotin for molecule to which
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
  loopShrC(i, 1, Nm)
    {
      loop(a, 1, NA)
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
  COORD_TYPE dx, dy, dz, dt2, rma, rmb;
  COORD_TYPE rxi[NA], ryi[NA], rzi[NA], vxi[NA], vyi[NA], vzi[NA];
 
  int i, a, b, it;
  /* ******************************************************************* */
  dt2 = dt / 2.0;
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;
  Mtot = 0.0;
  loop(a, 1, NA)
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
  loopShrC(i, 1, Nm)
    {
      /* VELOCITY VERLET ALGORITHM PART B */
      loop(a, 1, NA)
	{
	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];
	  vxi[a] = vx[a][i] + dt2 * Fx[a][i] / m[a];
	  vyi[a] = vy[a][i] + dt2 * Fy[a][i] / m[a];
	  vzi[a] = vz[a][i] + dt2 * Fz[a][i] / m[a];
	
	  moving[a] = 0;
	  moved[a] = 1;
	}

      /* START OF ITERATIVE LOOP */
      it = 0;
      done = 0;
      while ((done == 0 ) && ( it <= maxIt ))
	{
	  done = 1;
	  loop(a, 1, NB)
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
	  loop(a, 1, NA)
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
    loop(a, 1, NA)
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
  poolShrC;
  /* END OF LOOP OVER MOLECULES */
  
  loopShrC(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }

  loop(a, 1, NA)
    {
      K  = K + Ka[a] * m[a] * 0.5;
    }
     
  ProcSync0();
  WC = WC / dt2 / 3.0;

#if !defined(MOLPTENS) || defined(ATPTENS)
  WCxy = WCxy / dt2; /* WCxy, ... are not exactly virial terms and the 3.0
			is not present */
  WCyz = WCyz / dt2;
  WCzx = WCzx / dt2;
#ifndef NO_PARALLEL_CODE
  sumAllProcMA(&WCxy, &WCyz, &WCzx, &WCxx, &WCyy, &WCzz, 
	       &T1xy, &T1yz, &T1zx, &T1xx, &T1yy, &T1zz, NULL);
#endif

#endif

#ifndef NO_PARALLEL_CODE  
  sumAllProcMA(&WC, &K, NULL);
#endif

  /* !!!!!##$$%%^^*/
#ifdef MOLPTENS
#ifndef NO_PARALLEL_CODE
  sumAllProcMA(&T1mxy, &T1myz, &T1mzx,
	       &T1mxx, &T1myy, &T1mzz, NULL);
#endif
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
  loopShrC(i, 1, Oparams.parnum)
    {
      CoM(i, &Rx, &Ry, &Rz);
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * Rx);
      DRy = - L * rint(invL * Ry);
      DRz = - L * rint(invL * Rz);
      
      loop(a, 1, NA)
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

/* ========================== >>> updateAng <<< ============================ */
void updateAng(int Nm)
{
  int i;
  COORD_TYPE r01x, r01y, r01z, dt;
  COORD_TYPE v01x, v01y, v01z, r01Sq, dSq, d;

  dSq = Sqr(Oparams.d);
  d = Oparams.d;
  
  /* Preliminar syncronization */
  ProcSync0();
 
  /* Only Father calculate mesures so it only updates angular accumulators */
  doFather{ 
    dt = Oparams.steplength;
    loop(i, 1, Oparams.parnum)
      {
	r01x = rx[0][i] - rx[1][i];
	r01y = ry[0][i] - ry[1][i];
	r01z = rz[0][i] - rz[1][i];

	/* Store molecular orientations */
	ux[i] = r01x / d;
	uy[i] = r01y / d;
	uz[i] = r01z / d;
	
	/*r01x = r01x - L * rint(invL * r01x);
	r01y = r01y - L * rint(invL * r01y);
	r01z = r01z - L * rint(invL * r01z);*/
	
	v01x = vx[0][i] - vx[1][i];
	v01y = vy[0][i] - vy[1][i];
	v01z = vz[0][i] - vz[1][i];
	
	/*r01Sq = Sqr(r01x) + Sqr(r01y) + Sqr(r01z);*/
	vectProd(r01x, r01y, r01z, v01x, v01y, v01z, &ox[i], &oy[i], &oz[i]);
	/* The angular velocities is a global array so now we have the actual
	   angular velocity */
	ox[i] = ox[i] / dSq;
	oy[i] = oy[i] / dSq;
	oz[i] = oz[i] / dSq;
	OprogStatus.sumox[i] += ox[i] * dt; 
	OprogStatus.sumoy[i] += oy[i] * dt;
	OprogStatus.sumoz[i] += oz[i] * dt;

	/* Time integral of ang. vel. stored in xva variables */
	Dphix[i] = OprogStatus.sumox[i];
	Dphiy[i] = OprogStatus.sumoy[i];
	Dphiz[i] = OprogStatus.sumoz[i];
      }
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



/* ============================ >>> move<<< =================================*/
void move(void)
{
  /* DESCRIPTION:
     Move the particles by one step; this procedure is PARALLELIZED !!!
     Valid shared counters are: scN where N = 0..9 
     WARNING: "outer loops" should be ever parallelized, that is you must use 
     loopShr() macro for them inside move() and its subroutine 
     NOTE:
     linked list building could not be parallelized because father and child 
     could disturb each other */
 
  /* calc predicted coords*/
  movea(Oparams.steplength, 0.0000000001, 30, 1, Oparams.d, Oparams.m, 
	Oparams.parnum);        

  if (nebrNow)
    {
      //printf("building neighbour listr step: %d\n", Oparams.curStep);
      nebrNow = 0;
      dispHi = 0.0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      if (OprogStatus.noLinkedList)
	{
	  BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut, Oparams.sigab);
	}
      else
	{
	  links(Oparams.parnum, Oparams.rcut, Oparams.sigab);
	  
	  /* Build up neighbour list */  
	  BuildNebrList(Oparams.parnum, Oparams.rcut, Oparams.sigab);
	}
    }
  
  LJForce(Oparams.parnum, Oparams.epsab, Oparams.sigab,/* calculate forces */
	  Oparams.rcut);

  /* correct the coords */
  if (OprogStatus.Nose == 1)
    {  
      /* NPT ensemble */
      movebNPT(Oparams.steplength, 0.0000000001, 30, 1, Oparams.m, Oparams.d, 
	       Oparams.parnum);             
    }
  else if (OprogStatus.Nose == 2)
    {
      /* NPT ensemble */
      movebNTV(Oparams.steplength, 0.0000000001, 30, 1, Oparams.m, Oparams.d, 
	       Oparams.parnum);             
    }
  else    
    {
      /* NVE ensemble */
      moveb(Oparams.steplength, 0.0000000001, 30, 1, Oparams.m, Oparams.d, 
	    Oparams.parnum);             
    }
  /* Calculate the kinetic energy */
  kinet(Oparams.parnum, vx, vy, vz, Vol1);
  
  checkNebrRebuild();
  if ( (OprogStatus.Nose == 1) || (OprogStatus.Nose == 2))
    {
      //scalCor(Oparams.parnum);
      if ( (OprogStatus.sResetSteps != 0) &&
	  (Oparams.curStep == OprogStatus.sResetSteps) )
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
  
  if ( (OprogStatus.CMreset != 0) &&
      // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
      (Oparams.curStep == OprogStatus.CMreset) ) 
      resetCM(Oparams.parnum);

  /* Update accumulators for calculating the angular diffusion coefficent */
  updateAng(Oparams.parnum);  
  
  /* Update the integral of the pressure tensor */
  updateDQ(Oparams.steplength);
}

