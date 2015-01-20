#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern void links(COORD_TYPE rcut, COORD_TYPE sigab[NA][NA]);
extern void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]);
extern void BuildNebrList(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]);

extern void checkNebrRebuild(void);
extern void LJForce(COORD_TYPE epsab[NA][NA], 
		    COORD_TYPE sigab[NA][NA], COORD_TYPE rcut);
extern void resetCM();

extern void kinet(COORD_TYPE** velx, COORD_TYPE** vely, 
		  COORD_TYPE** velz, COORD_TYPE VOL1);
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE V, W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx;  

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

/* ========================================================================= */

/* ========================== >>> movea <<< =============================== */
void movea(COORD_TYPE dt, COORD_TYPE m[NA])
{  

  COORD_TYPE  dt2, dtSq2;
  COORD_TYPE  axia, ayia, azia;
  int i, a;
        
  L = cbrt(Vol);
  invL = 1.0 / L;

  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;

  /* ===== LOOP OVER MOLECULES ===== */
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* ====== >>>> VELOCITY VERLET ALGORITHM PART A <<< ======= */
	  axia = Fx[a][i] / m[a];
	  ayia = Fy[a][i] / m[a];
	  azia = Fz[a][i] / m[a];
		  
	  rx[a][i] = rx[a][i] + dt * vx[a][i] + dtSq2 * axia;
	  ry[a][i] = ry[a][i] + dt * vy[a][i] + dtSq2 * ayia;
	  rz[a][i] = rz[a][i] + dt * vz[a][i] + dtSq2 * azia;
	  
	  vxt[a][i] = vx[a][i]; /* v(t) */
	  vyt[a][i] = vy[a][i];
	  vzt[a][i] = vz[a][i];
	  
	  vx[a][i] = vx[a][i] + dt2 * axia;
	  vy[a][i] = vy[a][i] + dt2 * ayia;
	  vz[a][i] = vz[a][i] + dt2 * azia;
	  
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

/* ======================= >>> calcT1diagAt <<< =========================== */
COORD_TYPE  calcT1diagAt(COORD_TYPE VOL1)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of atomic pressure tensor */
  int i, a;
  COORD_TYPE kin, dlnV, *m;
  
  m = Oparams.m;
  kin = 0.0;
  dlnV = VOL1 / Vol / 3.0;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  kin +=  m[a]*(Sqr(vx[a][i] - dlnV * rx[a][i]) + 
			Sqr(vy[a][i] - dlnV * ry[a][i]) + 
			Sqr(vz[a][i] - dlnV * rz[a][i]));
	}
    }
  
  kin /= 3.0 * Vol;
  return kin;
}


/*============================ >>> calcPtensAt <<< =========================*/
void calcPtensAt(COORD_TYPE VOL1)
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
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  px = vx[a][i] - dlnV * rx[a][i];
	  py = vy[a][i] - dlnV * ry[a][i];
	  pz = vz[a][i] - dlnV * rz[a][i];
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
  Pxy = T1xy + Wxy;
  Pyz = T1yz + Wyz;
  Pzx = T1zx + Wzx;  
  Pxx = T1xx + Wxx;
  Pyy = T1yy + Wyy;
  Pzz = T1zz + Wzz;
  //printf("Wxx:%f WCxx:%f Wyy: %f WCyy: %f Wzz: %f WCzz: %f\n",
	//Wxx, WCxx, Wyy, WCyy, Wzz, WCzz);
  Pxx /= Vol;
  Pyy /= Vol;
  Pzz /= Vol;
  Pxy /= Vol;
  Pyz /= Vol;
  Pzx /= Vol;
}
#ifdef MD_GRAVITY
/* ========================= >>> moveb <<< =========================== */
void movebdamped(COORD_TYPE dt, COORD_TYPE m[NA])
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE Fxdamp, Fydamp, Fzdamp, damp; 
  COORD_TYPE dt2; 
  int i, a, k;
  const int NUMIT = 4;
  double Nm;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];

  /* ******************************************************************* */
  dt2 = dt * 0.5;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vxt2[a][i] = vx[a][i]; /* v*t2 = v*(t+dt/2) */
	  vyt2[a][i] = vy[a][i];
	  vzt2[a][i] = vz[a][i];

	  FxLJ[a][i] = Fx[a][i];
	  FyLJ[a][i] = Fy[a][i];
	  FzLJ[a][i] = Fz[a][i];
	
	  vx[a][i] = 5.0 * vxt[a][i] / 2.0 - 2.0 * vxo1[a][i] + 
	    vxo2[a][i] / 2.0;
	  vy[a][i] = 5.0 * vyt[a][i] / 2.0 - 2.0 * vyo1[a][i] + 
	    vyo2[a][i] / 2.0;
	  vz[a][i] = 5.0 * vzt[a][i] / 2.0 - 2.0 * vzo1[a][i] + 
	    vzo2[a][i] / 2.0;
	}
    }

  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      for(a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
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
	      vx[a][i] /= 1.0 + 0.5 * dt * Oparams.partdamp;
	      vy[a][i] /= 1.0 + 0.5 * dt * Oparams.partdamp;
	      vz[a][i] /= 1.0 + 0.5 * dt * Oparams.partdamp;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      Fxdamp = - Oparams.partdamp * vx[a][i] * m[a];
	      Fydamp = - Oparams.partdamp * vy[a][i] * m[a];
	      Fzdamp = - Oparams.partdamp * vz[a][i] * m[a];
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[a][i] = FxLJ[a][i] + Fxdamp;
	      Fy[a][i] = FyLJ[a][i] + Fydamp;
	      Fz[a][i] = FzLJ[a][i] + Fzdamp;
	    } 
	}
      /* Zero the velocity along bond and calculate the constraint virial, this
	 should be as before because the Andersen Forces don't act along 
	 bond.*/
      

    }
  kinet(vx, vy, vz, 0.0);/* K(t+dt) */

#if 0
  calcPtensAt(0.0);
  press = (Pxx + Pyy + Pzz) / 3.0;
#endif
  press = calcT1diagAt(0.0) + W / Vol;

  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
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
void movebNTV(COORD_TYPE dt, COORD_TYPE m[NA])
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, dlns; 
  COORD_TYPE s1i;
  COORD_TYPE DT, A, dt2; 
  int i, a, k;
  const int NUMIT = 4;
  double Nm;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];

  /* ******************************************************************* */
  dt2 = dt * 0.5;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
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
  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(vx, vy, vz, 0.0);  /* K(t+dt) */
      /* CALCOLARE QUI CORRETTAMENTE LA TEMPERATURA DELLA MISTURA!!!!!!*/
      DT = s * (2.0 * K - (3.0 * Nm - 3.0) * Oparams.T) / OprogStatus.Q;
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
#ifdef MD_GRAVITY
      dlns = s1 /s + Oparams.partdamp ;
#else
      dlns = s1 / s ; 
#endif     
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      for(a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
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
      /* Zero the velocity along bond and calculate the constraint virial, this
	 should be as before because the Andersen Forces don't act along 
	 bond.*/
      

    }
  kinet(vx, vy, vz, 0.0);/* K(t+dt) */

#if 0
  calcPtensAt(0.0);
  press = (Pxx + Pyy + Pzz) / 3.0;
#endif
  press = calcT1diagAt(0.0) + W / Vol;

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  s1o2  = s1o1;
  s1o1 = s1t;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
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
void movebNPT(COORD_TYPE dt, COORD_TYPE m[NA])
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, FxAnd, FyAnd, FzAnd, cfAnd, dlns, dlnV,
    dlnVSq; 
  COORD_TYPE s1i, Vol1i;
  COORD_TYPE DT, A, B, DP, dt2; 
  int i, a, k;
  const int NUMIT = 4;
  double Nm;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  /* ******************************************************************* */
  dt2 = dt / 2.0;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
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

  //printf("Vol1: %f Vol1o1: %f Vol1o2: %f Vol1t:%f\n", 
	// Vol1, Vol1o1, Vol1o2, Vol1t);
  for(k = 0; k < NUMIT; k++) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(vx, vy, vz, Vol1);  /* K(t+dt) */
      //printf("K guess: %.20f\n", K);
      // SISTEMARE QUI BISOGNA CALCOLARE CORRETTAMENTE LA TEMPERATURA
      // DELLA MISTURA !!!!!!!!!!!!
      DT = s * (2.0 * K - (3.0 * Nm - 3.0) * Oparams.T) / OprogStatus.Q;
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      
      /* Calculate pressure, calcT1diagMol is a term proportional to the 
	 translational kinetic energy, see Ferrario and Ryckaert */
      press = calcT1diagAt(Vol1) + W  / Vol; /* press(t+dt) */
      /* Volume acceleration */
      DP = Sqr(s) * (press - Oparams.P) / OprogStatus.W;

      B = Vol1i + 0.5 * dt * DP;
      Vol1 = B / (1.0 - 0.5 * dt * s1 / s); /* Vol1(t+dt) */
      Vol2 = DP + s1 * Vol1 / s;            /* Vol2(t+dt) */
      
      /* Calculate the Andersen forces */
      dlnV   = Vol1 / Vol; 
      dlnVSq = Sqr(dlnV);
      cfAnd  = - 2.0 * dlnVSq / 9.0;
      cfAnd += Vol2 / Vol / 3.0;
      cfAnd += dlns * dlnV / 3.0;
      for(a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      /* Calculate center of mass posiotin for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	    
	      FxAnd  = cfAnd * rx[a][i] * m[a];
	      FyAnd  = cfAnd * ry[a][i] * m[a];
	      FzAnd  = cfAnd * rz[a][i] * m[a];
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
    }

  /* Calculate kinetic energy */
  kinet(vx, vy, vz, Vol1);/* K(t+dt) */
  //printf("K exact: %.20f\n", K);
  
  press = calcT1diagAt(Vol1) + W / Vol;

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  Vol1o2 = Vol1o1;
  s1o2  = s1o1;
  Vol1o1 = Vol1t;
  s1o1 = s1t;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
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


#if !defined(MD_GRAVITY)
/* ========================= >>> moveb <<< =========================== */
void moveb(COORD_TYPE dt, COORD_TYPE m[NA])
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE Ka[NA];
  COORD_TYPE dt2;
  int i, a;
  /* ******************************************************************* */
  dt2 = dt / 2.0;
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  for(a = 0; a < NA; a++)
    {
      Ka[a] = 0.0;
    }

  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;

  /* LOOP OVER ALL MOLECULES */
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vx[a][i] = vx[a][i] + dt2 * Fx[a][i] / m[a];
	  vy[a][i] = vy[a][i] + dt2 * Fy[a][i] / m[a];
	  vz[a][i] = vz[a][i] + dt2 * Fz[a][i] / m[a];

	  T1xy += vx[a][i] * vy[a][i] * m[a]; 
	  T1yz += vy[a][i] * vz[a][i] * m[a];
	  T1zx += vz[a][i] * vx[a][i] * m[a];
	  T1xx += vx[a][i] * vx[a][i] * m[a]; 
	  T1yy += vy[a][i] * vy[a][i] * m[a];
	  T1zz += vz[a][i] * vz[a][i] * m[a];

	  Ka[a] = Ka[a] + Sqr(vx[a][i]) + Sqr(vy[a][i]) + Sqr(vz[a][i]);
	}
    }
  /* END OF LOOP OVER MOLECULES */
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vxo2[a][i] = vxo1[a][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[a][i] = vyo1[a][i];
	  vzo2[a][i] = vzo1[a][i];
	  vxo1[a][i] = vxt[a][i];  /* vo1(t+dt) = v(t) */
	  vyo1[a][i] = vyt[a][i];
	  vzo1[a][i] = vzt[a][i];
	}
    }
  
  for(a = 0; a < NA; a++)
    {
      K  = K + Ka[a] * m[a] * 0.5;
    }
#if 0
  Pxy = T1xy + Wxy;
  Pyz = T1yz + Wyz;
  Pzx = T1zx + Wzx;  
  Pxy /= Vol;
  Pyz /= Vol;
  Pzx /= Vol;
  press = (T1xx + T1yy + T1zz)/3.0 + W;
  press /= Vol;
#endif
  
  press = (T1xx + T1yy + T1zz)/3.0 + W;
  press /= Vol;

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
void scalCor()
{ 
  int i, a;
  COORD_TYPE cost, costo1, costo2, DRx, DRy, DRz;
  
  L = cbrt(Vol);
  invL = 1.0 / L;
  cost = Vol1 / Vol / 3.0;
  costo1 = Vol1o1 / Vol / 3.0;
  costo2 = Vol1o2 / Vol / 3.0;
  
  /* Reduced particles to first box */
  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	     scale them */
	  DRx = - L * rint(invL * rx[a][i]);
	  DRy = - L * rint(invL * ry[a][i]);
	  DRz = - L * rint(invL * rz[a][i]);

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
      temp = 2.0 * K / (3.0 * (Oparams.parnum[0]+Oparams.parnum[1]) - 3.0);
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

void reassignvel(void)
{
  int a,i;
  
  s = 1.0; 
  s1 = 0.0;
  s2 = 0.0;

  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vxo1[a][i] = 0.0;
	  vyo1[a][i] = 0.0;
	  vzo1[a][i] = 0.0;
	  vxo2[a][i] = 0.0;
	  vyo2[a][i] = 0.0;
	  vzo2[a][i] = 0.0;
	}
    }

  comvel(Oparams.T, Oparams.m);
  nebrNow = 0;
  dispHi = 0.0;
  /* build up linked list on predicted 
     coordinates (only father do it)*/
  if (OprogStatus.noLinkedList)
    {
      BuildNebrListNoLinked(Oparams.rcut, Oparams.sigab);
    }
  else
    {
      links(Oparams.rcut, Oparams.sigab);
      
      /* Build up neighbour list */  
      BuildNebrList(Oparams.rcut, Oparams.sigab);
    }

  /* ricalcola le forze poiché ora il sistema è nel minimo locale */ 
  LJForce(Oparams.epsab, Oparams.sigab,/* calculate forces */
	  Oparams.rcut);
  printf("velocities reassinged\n");
}
extern void conjgrad(void);
#ifdef MD_GRAVITY
void awall(void)
{
  double L, L2;
  int a, i;
  const double wallFriction = 1.0;
  /* qui bisogna tener conto della fondo del recipiente */
  L = cbrt(Vol);
  L2 = L/2.0;
  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* urto con la parete */
	  if (rz[a][i] < -L2)
	    vz[a][i] = -wallFriction*vz[a][i];
	}
    }	
}
#endif
int check_equilibrat(void)
{
  int a, i;
  double Drx, Dry, Drz, DrSq[NA];
  
  for (a = 0; a < 2; a++)
    {
      DrSq[a] = 0.0;
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  Drx = rx[a][i] - OprogStatus.rxi[a][i];
	  Dry = ry[a][i] - OprogStatus.ryi[a][i];
	  Drz = rz[a][i] - OprogStatus.rzi[a][i];
	  DrSq[a] = DrSq[a] + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
	}
      DrSq[a] /= ((double) Oparams.parnum[a]);
    }
  /* Qui si è scelto come criterio di equilibratura che le particelle
   * si siano spostate in media più di 1 sigma */
  if (DrSq[0] > Oparams.sigab[0][0])
    return 1;
  else 
    return 0;
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
  int Nm; 
#ifdef MPI
  int equilib, sumeq, i, eqarray[128];
#endif
  /* calc predicted coords*/
  movea(Oparams.steplength, Oparams.m);
  if (nebrNow)
    {
      /*printf("building neighbour listr step: %d\n", Oparams.curStep);*/
      nebrNow = 0;
      dispHi = 0.0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      if (OprogStatus.noLinkedList)
	{
	  BuildNebrListNoLinked(Oparams.rcut, Oparams.sigab);
	}
      else
	{
	  links(Oparams.rcut, Oparams.sigab);
	  
	  /* Build up neighbour list */  
	  BuildNebrList(Oparams.rcut, Oparams.sigab);
	}
    }
  
  LJForce(Oparams.epsab, Oparams.sigab,/* calculate forces */
	  Oparams.rcut);
#if !defined MD_GRAVITY
  /* correct the coords */
  if (OprogStatus.Nose == 1)
    {  
      /* NPT ensemble */
      movebNPT(Oparams.steplength, Oparams.m);
    }
  else if (OprogStatus.Nose == 2)
    {
      /* NPT ensemble */
      movebNTV(Oparams.steplength, Oparams.m);
    }

  else    
    {
      /* NVE ensemble */
      moveb(Oparams.steplength, Oparams.m);
    }
#else
  if (OprogStatus.Nose == 1)
    {  
      /* NPT ensemble */
      movebNPT(Oparams.steplength, Oparams.m);
    }
  else if (OprogStatus.Nose == 2) 
    {
      /* NPT ensemble */
      movebNTV(Oparams.steplength, Oparams.m);
    }
  else
    {
      movebdamped(Oparams.steplength, Oparams.m);
    }
#endif
#ifdef MD_GRAVITY
  awall();
  
  if (OprogStatus.taptausteps > 0)
    {
      if (OprogStatus.quenchend == -1)
	{
	  Nm = Oparams.parnum[0] + Oparams.parnum[1];
	  if ( (2.0*K/(3.0*((double)Nm)-3.0)/Oparams.T) < 
	       OprogStatus.quenchtol)
	    {
	      printf("STEP: %d QUENCH DONE!!!\n", Oparams.curStep);
	      comvel(Oparams.parnum, Oparams.T, Oparams.m);
	    }
	}
      else
	{
	   if ((Oparams.curStep - OprogStatus.quenchend)  >= OprogStatus.taptausteps)
	     {
	       OprogStatus.quenchend = -1;
	     }
	}
    }
#endif
  /* Calculate the kinetic energy */
  kinet(vx, vy, vz, Vol1);
  
  checkNebrRebuild();
  if ( (OprogStatus.Nose == 1) || (OprogStatus.Nose == 2))
    {
      //scalCor();
      if ( ( (OprogStatus.sResetSteps > 0) &&
	     (Oparams.curStep == OprogStatus.sResetSteps) )
	   || ( (OprogStatus.sResetSteps < 0) &&
		(Oparams.curStep % (-OprogStatus.sResetSteps) == 0) )  )
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
  
  if (  ( (OprogStatus.CMreset > 0) &&
	  // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	  (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	     // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	     (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
    resetCM(Oparams.parnum);
#ifdef TAPPING 
  if ( OprogStatus.tapping && (Oparams.curStep % OprogStatus.taptau == 0))
    {
      conjgrad();
      reassignvel();
    } 
#endif
  
  if (OprogStatus.chkeqstps && 
      Oparams.curStep % OprogStatus.chkeqstps == 0)
    {
      if (!OprogStatus.equilibrated)
	OprogStatus.equilibrated = check_equilibrat();
#ifdef MPI
      equilib = OprogStatus.equilibrated;
      MPI_Allgather(&equilib, 1, MPI_INT, 
		    equilibrated, 1, MPI_INT, MPI_COMM_WORLD);
      sumeq = 0;
      for (i=0; i < numOfProcs; i++)
	sumeq += equilibrated[i];
      /* se sumeq = numOfProcs vuol dire che tutti i processi sono
       * equilibrati quindi la simulazione può terminare */
      if (sumeq == numOfProcs)
	{
	  mdPrintf(ALL,"All systems reached equilibrium, simulation completed\n", 
		   NULL);
	  printf("[MSDcheck] time= %.15G\n", Oparams.curStep*Oparams.steplength);
	  ENDSIM = 1;
	}
#if 1
      else
	printf("rank#%d not all systems reached equilibrium sumeq: %d\n", my_rank, sumeq);
#endif
#else
      mdPrintf(ALL, "All systems reached equilibrium, simulation completed",
	       NULL);
      ENDSIM = OprogStatus.equilibrated;
#endif
    }
  
  /* Update the integral of the pressure tensor */
#if 0
  updateDQ(Oparams.steplength);
#endif
}

