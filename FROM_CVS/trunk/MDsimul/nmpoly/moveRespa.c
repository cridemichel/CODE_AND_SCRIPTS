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
#if 0
extern double **rallx, **rally, **rallz, **Fallx, **Fally, **Fallz,
  **rallx_old, **rally_old, **rallz_old, **vxold, **vyold, **vzold;
#endif
extern double **vxold, **vyold, **vzold, **rx_old, **ry_old, **rz_old;
#ifdef MD_RESPA
extern double  **rx_oldLong, **ry_oldLong, **rz_oldLong;
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
  WmxyLong, WmyzLong, WmzxLong, WmyxLong, WmzyLong, WmxzLong;
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
#ifdef MD_RESPA_NPT
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
		  /* la catena è aperta */
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
#ifdef MD_RESPA_NPT
void p2v(void)
{ 
  int i, a;
  double Rxl, Ryl, Rzl;
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
  Vol1 = Pv * Sqr(s) / OprogStatus.W;
  s1   = Ps * Sqr(s) / OprogStatus.Q;  
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
#ifndef MD_FENE
  shakeVelRespaNPT(Oparams.parnum, Oparams.steplength, Oparams.m, 150, NA-1, Oparams.d, 0.000000000001, vx, vy, vz);
#endif
}
void updNoseAnd(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], mM[NA], cost[NA];
  double dlns, dlnV;
  double PCMx, PCMy, PCMz;
  cdt = c*dt;
  dlns = s*Ps / OprogStatus.Q ;
  dlnV = Pv*Sqr(s)/OprogStatus.W/3.0/Vol; 
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      expdt[a] = exp(-(dlns + dlnV * mM[a]) * cdt);
      cost[a] = (expdt[a] - 1.0) * dlnV * mM[a] / (dlns + dlnV*mM[a]);
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
	  px[a][i] = cost[a]*(PCMx - px[a][i]) + px[a][i]*expdt[a];
	  py[a][i] = cost[a]*(PCMy - py[a][i]) + py[a][i]*expdt[a];
	  pz[a][i] = cost[a]*(PCMz - pz[a][i]) + pz[a][i]*expdt[a];
	}
    }
}
void updNoseAndRef(double dt, double c)
{
  int i, a;
  double cdt, expdt[NA], mM[NA];
  double cost[NA];
  double RCMx, RCMy, RCMz;
  cdt = c*dt;
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      expdt[a] = exp(Pv*Sqr(s)*mM[a]*cdt/(3.0*Vol*OprogStatus.W));
      cost[a] = (expdt[a] - 1.0);
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      CoM(i, &RCMx, &RCMy, &RCMz);
      for (a=0; a < NA; a++)
	{
	  rx[a][i] = cost[a]*(RCMx - rx[a][i]*mM[a]) + rx[a][i]*expdt[a];
	  ry[a][i] = cost[a]*(RCMy - ry[a][i]*mM[a]) + ry[a][i]*expdt[a];
	  rz[a][i] = cost[a]*(RCMz - rz[a][i]*mM[a]) + rz[a][i]*expdt[a];
	}
    }

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
  s = s / (1 - s*Ps*cdt2/OprogStatus.Q);
  Kin = 0;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	Kin += Sqr(px[a][i])+Sqr(py[a][i])+Sqr(pz[a][i])/Oparams.m[a];  
      }
  Kin *= 0.5;
#ifdef MD_FENE
  dof = 3*NA*Oparams.parnum;
#else
  dof = (2*NA - 1)*Oparams.parnum;
#endif
  //printf("temp= %f\n", 2*Kin/dof);
  DT =  (2.0 * Kin - (dof - 3.0) * Oparams.T)/s;
  //printf("DT: %f Kin: %f\n", DT, Kin);
  Ps += DT * cdt2;
  Ps = Ps / (1 - Ps*cdt*s/OprogStatus.Q);
  Ps += DT * cdt2;
  s = s / (1 - s*Ps*cdt2/OprogStatus.Q);
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

void updLv(double dt, double c)
{
  double press, cdt, cdt2, DP, Nm;
  Nm = Oparams.parnum;
  cdt = c * dt;
  cdt2 = c * dt / 2.0;
  Vol += cdt2*Sqr(s)*Pv/OprogStatus.W;
#ifdef MOLPTENS
  press = calcT1diagMolRespa(Nm) + Wm / 3.0 / Vol; /* press(t+dt) */
#else
  press = calcT1diagAtRespa(Nm) + (W + WC) / Vol; /* press(t+dt) */
#endif
  printf("press: %f\n", press);
  DP = press - Oparams.P;
  Pv += DP  * cdt2;
  Pv *= exp(-Ps*s/OprogStatus.Q);
  Pv += DP  * cdt2;
  Vol += cdt2*Sqr(s)*Pv/OprogStatus.W;
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
  updImpLong(dt, 0.25);
  printf("1) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updNoseAnd(dt, 0.25);
  printf("2) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updLv(dt, 0.25);
  printf("3) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updLs(dt, 0.5);
  printf("4) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updLv(dt, 0.25);
  printf("5) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updNoseAnd(dt, 0.25);
  printf("6) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  LJForceLong(Oparams.parnum, Oparams.rcut, Oparams.rcut);
  updImpLong(dt, 0.25);
  updNoseAndRef(dt, 0.5); 
  printf("7) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
}

void movelongRespaNPTAft(double dt)
{
  updNoseAndRef(dt, 0.5); 
  printf("A1) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  LJForceLong(Oparams.parnum, Oparams.rcut, Oparams.rcut);
  updImpLong(dt, 0.25);
  printf("A2) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updNoseAnd(dt, 0.25);
  printf("A3) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updLv(dt, 0.25);
  printf("A4) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updLs(dt, 0.5);
  printf("A5) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updLv(dt, 0.25);
  printf("A6) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updNoseAnd(dt, 0.25);
  printf("A7) Pv: %f Ps: %f s: %f Vol: %f\n", Pv, Ps, s, Vol);
  updImpLong(dt, 0.25);
}
/* ========================== >>> movea <<< =============================== */
void moveaRespa(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB, COORD_TYPE d, 
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
  double expdt[NA], mM[NA];
  double cost[NA];
  double RCMx, RCMy, RCMz;

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
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      expdt[a] = exp(Pv*Sqr(s)*mM[a]*dt/(3.0*Vol*OprogStatus.W));
      cost[a] = (expdt[a] - 1.0);
    }
#ifdef MOLPTENS
  press = calcT1diagMolRespa(Nm) + (WmLong + WmShort) / 3.0 / Vol; /* press(t+dt) */
#else
  press = calcT1diagAtRespa(Nm) + (WLong + WShort + WC) / Vol; /* press(t+dt) */
#endif
  DP = press - Oparams.P;
  Pv += DP  * dt2;
  
  /* ===== LOOP OVER MOLECULES ===== */
  for (i=0; i < Nm; i++)
    {
      CoM(i, &RCMx, &RCMy, &RCMz);
      /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
      for(a=0; a < NA; a++)
	{
#if 0
	  axia = Fx[a][i] / m[a];
	  ayia = Fy[a][i] / m[a];
	  azia = Fz[a][i] / m[a];
#endif
	  rxi[a] = rx[a][i];
	  ryi[a] = ry[a][i];
	  rzi[a] = rz[a][i];
	  px[a][i] = px[a][i] + dt2 * Fx[a][i];
	  py[a][i] = py[a][i] + dt2 * Fy[a][i];
	  pz[a][i] = pz[a][i] + dt2 * Fz[a][i];
	  
	  pxi[a][i] = dt * px[a][i]/m[a] + cost[a]*(RCMx - rx[a][i]*mM[a]) + rx[a][i]*expdt[a];
	  pyi[a][i] = dt * py[a][i]/m[a] + cost[a]*(RCMy - ry[a][i]*mM[a]) + ry[a][i]*expdt[a];
	  pzi[a][i] = dt * pz[a][i]/m[a] + cost[a]*(RCMz - rz[a][i]*mM[a]) + rz[a][i]*expdt[a];

#if 0 
	  printf("[%d-%d] v=(%f,%f,%f) a=(%f,%f,%f) pr=(%f,%f,%f) r=(%f,%f,%f) (%f,%f,%f)\n", a, i, 
		 vx[a][i], vy[a][i], vz[a][i],
		 axia, ayia, azia, pxi[a], pyi[a], pzi[a],
		 rxi[a], ryi[a], rzi[a], Fx[a][i], Fy[a][i], Fz[a][i]);
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
#if !defined(MD_FENE)
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
				 sqrt(Sqr(rxab)+Sqr(ryab)+Sqr(rzab)));
			  printf("rpab=%.15f\n", rpab);
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
#endif
      for(a=0; a < NA; a++)
	{
	  rx[a][i] = pxi[a];
	  ry[a][i] = pyi[a];
	  rz[a][i] = pzi[a];

	  vx[a][i] = vxi[a];
	  vy[a][i] = vyi[a];
	  vz[a][i] = vzi[a];
	  //printf("v(%f,%f,%f) F(%f,%f,%f)\n", vx[a][i], vy[a][i], vz[a][i],
	//	 Fx[a][i], Fy[a][i], Fz[a][i]);
	}
      
    }
  
  Vol += dt * Sqr(s)*Pv/OprogStatus.W;

}

/* ========================= >>> moveb <<< =========================== */
void movebRespa(COORD_TYPE dt, COORD_TYPE tol, int maxIt, int NB,
	   COORD_TYPE m[NA], COORD_TYPE d, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  int done;
  int moving[NA], moved[NA];
  COORD_TYPE Ka[NA], Mtot, pmx, pmy, pmz;
  COORD_TYPE rxab, ryab, rzab, rvab, gab, dSq;
  COORD_TYPE vxab, vyab, vzab;
  COORD_TYPE dx, dy, dz, dt2, rma, rmb, c2dt;
  COORD_TYPE rxi[NA], ryi[NA], rzi[NA], pxi[NA], pyi[NA], pzi[NA];
  double c0, c1, c2,chsi; 
  int i, a, b, it;

  
  /* ******************************************************************* */
  dt2 = dt*0.5;
#if 0
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
#endif
  c2dt= dt / 2.0;
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
  for (a = 0; a < NA; a++)
    {
      mM[a] = Oparams.m[a] / Mtot;
      expdt[a] = exp(Pv*Sqr(s)*mM[a]*dt/(3.0*Vol*OprogStatus.W));
      cost[a] = (expdt[a] - 1.0);
    }
#ifdef MOLPTENS
  press = calcT1diagMolRespa(Nm) + (WmLong + WmShort) / 3.0 / Vol; /* press(t+dt) */
#else
  press = calcT1diagAtRespa(Nm) + (WLong + WShort + WC) / Vol; /* press(t+dt) */
#endif
  DP = press - Oparams.P;
  Pv += DP  * dt2;
 
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
	  vxi[a] = (px[a][i] + c2dt * Fx[a][i]) / m[a];
	  vyi[a] = (py[a][i] + c2dt * Fy[a][i]) / m[a];
	  vzi[a] = (pz[a][i] + c2dt * Fz[a][i]) / m[a];
	  moving[a] = 0;
	  moved[a] = 1;
	}
#if !defined(MD_FENE)
      /* START OF ITERATIVE LOOP */
      it = 0;
      done = 0;
      while ((done == 0 ) && ( it <= maxIt ))
	{
	  done = 1;
	  for(a=0; a < NB; a++)
	    {
	      b = a + 1;
	      //if (b >= NA) b = 0;

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
#endif
#ifdef MOLPTENS
    vmx = vmy = vmz = 0.0;  /* Velocity of center of mass of molecule i */
#endif
    for(a=0; a < NA;  a++)
	{
	  px[a][i] = vxi[a]*m[a];
	  py[a][i] = vyi[a]*m[a];
	  pz[a][i] = vzi[a]*m[a];
#ifdef MOLPTENS
	  pmx += pxi[a] * m[a];
	  pmy += pyi[a] * m[a];
	  pmz += pzi[a] * m[a];
#endif
#if !defined(MOLPTENS) || defined(ATPTENS)
	  /* Kinetic terms of the pressure-tensor */
	  T1xy += px[a][i] * py[a][i] / m[a]; 
	  T1yz += py[a][i] * pz[a][i] / m[a];
	  T1zx += pz[a][i] * px[a][i] / m[a];
	  T1xx += px[a][i] * px[a][i] / m[a]; 
	  T1yy += py[a][i] * py[a][i] / m[a];
	  T1zz += pz[a][i] * pz[a][i] / m[a];
#endif
	  Ka[a] = Ka[a] + Sqr(pxi[a]) + Sqr(pyi[a]) + Sqr(pzi[a]);
	}
#ifdef MOLPTENS
    /* Kinetic component of pressure tensor (all terms) */
    T1mxy += pmx * pmy / Mtot; 
    T1myz += pmy * pmz / Mtot;
    T1mzx += pmz * pmx / Mtot;
    T1mxx += pmx * pmx / Mtot;
    T1myy += pmy * pmy / Mtot;
    T1mzz += pmz * pmz / Mtot;
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

  if (Oparams.curStep == 1)  
    {
      BuildNebrListNoLinkedLong(Oparams.parnum, Oparams.rcut);
#if 0
      BuildNebrListNoLinked(Oparams.parnum, OprogStatus.rcutInner);
#endif
      LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
#if 0
      LJForce(Oparams.parnum, OprogStatus.rcutInner);
#ifdef MD_FENE
      FENEForce();
#endif
#endif
    } 
#ifdef MD_RESPA_NPT
  v2p();
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
      LJForce(Oparams.parnum, OprogStatus.rcutInner);
#ifdef MD_FENE
      FENEForce();
#endif
      /* considera tutti i contributi alle forza agente sugli atomi "di base"
       * ossia somma anche le forze dovute agli atomi senza massa moltiplicate
       * per gli opportuni coefficienti "vincolari" 
       * */
      Wm += WmLong;
#ifdef MOLPTENS
      Wmxy += WmyzLong;
      Wmyz += WmyzLong;
      Wmzx += WmzxLong;
#endif
#if defined(MOLPTENS)  || (!defined(ATPTENS) && !defined(ATPRESS))
      Wmxx += WmxxLong;
      Wmyy += WmyyLong;
      Wmzz += WmzzLong;
#endif
#if defined(ATPTENS)
      Wxx += WxxLong;
      Wyy += WyyLong;
      Wzz += WzzLong;
      Wxy += WxyLong;
      Wzx += WzxLong;
      Wyz += WyzLong;
#endif
      VcR = Vc;
      Vc += VcLong;
      V += VLong;
      W += WLong;
      kinet(Oparams.parnum, vx, vy, vz, Vol1);
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
  
#ifdef MD_RESPA_NPT
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	px[a][i] = vx[a][i] * Oparams.m[a];
	py[a][i] = vy[a][i] * Oparams.m[a];
	pz[a][i] = vz[a][i] * Oparams.m[a];
      }
#endif
 
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
  //LJForceLong(Oparams.parnum, Oparams.rcut, Oparams.rcut);
#ifdef MD_RESPA_NPT
  movelongRespaNPTAft(Oparams.steplength);
  p2v();
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
  /* Calculate the kinetic energy */
  kinet(Oparams.parnum, vx, vy, vz, Vol1);

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
