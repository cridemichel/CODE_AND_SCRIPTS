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
extern char TXT[MSG_LEN];
extern void  links(int Nm, COORD_TYPE rcut, COORD_TYPE sigma);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void BuildNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void checkNebrRebuild(void);
extern void LJForce(int Nm, COORD_TYPE epsilon, 
		    COORD_TYPE sigma, COORD_TYPE rcut);
extern void resetCM(int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void kinet(int Nm, COORD_TYPE* velx, COORD_TYPE* vely, 
		  COORD_TYPE* velz, COORD_TYPE VOL1);

extern zeroArrays(C_T *arrx, C_T* arry, C_T *arrz, int N);
extern sumFRContribs(void);
extern void buildMesh(C_T kmax, C_T Dk);
extern void ACForce(int Nm, COORD_TYPE kmax, COORD_TYPE Dk, C_T alpha, C_T S0);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE V, W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx,
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

/* ========================================================================= */

/* ========================== >>> movea <<< =============================== */
void movea(COORD_TYPE dt, COORD_TYPE m, int Nm)
{  
  COORD_TYPE  dt2, dtSq2;
  COORD_TYPE  axi, ayi, azi;
  
  int i;
 
  L = cbrt(Vol);
  invL = 1.0 / L;

  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;
 
  /* ===== LOOP OVER MOLECULES ===== */
  loop(i, 1, Nm)
    {
      axi = Fx[i] / m;
      ayi = Fy[i] / m;
      azi = Fz[i] / m;
      
      rx[i] = rx[i] + dt * vx[i] + dtSq2 * axi;
      ry[i] = ry[i] + dt * vy[i] + dtSq2 * ayi;
      rz[i] = rz[i] + dt * vz[i] + dtSq2 * azi;
	  
      vxt[i] = vx[i]; /* v(t) */
      vyt[i] = vy[i];
      vzt[i] = vz[i];
      
      vx[i] = vx[i] + dt2 * axi;
      vy[i] = vy[i] + dt2 * ayi;
      vz[i] = vz[i] + dt2 * azi;
    }
  pool;
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

/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagAt(int Nm, COORD_TYPE VOL1)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of atomic pressure tensor */
  int i;
  COORD_TYPE kin, dlnV, m;
  
  m = Oparams.m;
  kin = 0.0;
  dlnV = VOL1 / Vol / 3.0;
  loop(i, 1, Nm)
    {
	  kin +=  m*(Sqr(vx[i] - dlnV * rx[i]) + 
	    Sqr(vy[i] - dlnV * ry[i]) + 
	    Sqr(vz[i] - dlnV * rz[i]));
    }
  pool;
  /*
    #ifndef NO_PARALLEL_CODE
    sumAllProcMA(&kin, NULL);
    #endif
  */
  kin /= 3.0 * Vol;
  return kin;
}

/*============================ >>> calcPtensAt <<< =========================*/
void calcPtensAt(int Nm, COORD_TYPE VOL1)
{
  /* Calculate all components of atomic pressure tensor */
  int i;
  COORD_TYPE dlnV, px, py, pz, m;
  
  m = Oparams.m;
  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;

  dlnV = VOL1 / Vol / 3.0;
  loop(i, 1, Nm)
    {
      px = vx[i] - dlnV * rx[i];
      py = vy[i] - dlnV * ry[i];
      pz = vz[i] - dlnV * rz[i];
      /* Kinetic component of pressure tensor (all terms) */
      T1xy += px * py * m; 
      T1yz += py * pz * m;
      T1zx += pz * px * m ;
      T1xx += px * px * m;
      T1yy += py * py * m;
      T1zz += pz * pz * m;
    }
  pool;
  /*
    #ifndef NO_PARALLEL_CODE  
    sumAllProcMA(&T1xx, &T1yy, &T1zz, 
    &T1xy, &T1yz, &T1zx, NULL);
    #endif
  */
  
  /* Calculate the other element (off-diagonal) of the molecular 
     pressure tensor */
  Patxy = T1xy + Wxy;
  Patyz = T1yz + Wyz;
  Patzx = T1zx + Wzx;  
  Patxx = T1xx + Wxx;
  Patyy = T1yy + Wyy;
  Patzz = T1zz + Wzz;
  //printf("Wxx:%f WCxx:%f Wyy: %f WCyy: %f Wzz: %f WCzz: %f\n",
	//Wxx, WCxx, Wyy, WCyy, Wzz, WCzz);
  Patxx /= Vol;
  Patyy /= Vol;
  Patzz /= Vol;
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
}

void calcKtot();

/* ========================= >>> moveb <<< =========================== */
void movebNTV(COORD_TYPE dt, COORD_TYPE m, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, dlns; 
  COORD_TYPE s1i;
  COORD_TYPE DT, A, dt2; 
  int i, k;
  const int NUMIT = 4;

  /* ******************************************************************* */
  dt2 = dt * 0.5;
  
  loop(i, 1, Nm)
    {
      vxt2[i] = vx[i]; /* v*t2 = v*(t+dt/2) */
      vyt2[i] = vy[i];
      vzt2[i] = vz[i];
      
      FxLJ[i] = Fx[i];
      FyLJ[i] = Fy[i];
      FzLJ[i] = Fz[i];
      
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
      
      vx[i] = 5.0 * vxt[i] / 2.0 - 2.0 * vxo1[i] + 
	vxo2[i] / 2.0;
      vy[i] = 5.0 * vyt[i] / 2.0 - 2.0 * vyo1[i] + 
	vyo2[i] / 2.0;
      vz[i] = 5.0 * vzt[i] / 2.0 - 2.0 * vzo1[i] + 
	vzo2[i] / 2.0;
    }
  s1i = s1;    /* s1i = s1(t+dt/2) */

  //printf("Vol1: %f Vol1o1: %f Vol1o2: %f Vol1t:%f\n", 
	// Vol1, Vol1o1, Vol1o2, Vol1t);
  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm, vx, vy, vz, 0.0);  /* K(t+dt) */

      //printf("K guess: %.20f\n", K);

      DT = s * (2.0 * K - (3.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;

      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      
      loop(i, 1, Nm)
	{
	  /* Calculate center of mass posiotin for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	      Fx[i] = FxLJ[i];
	      Fy[i] = FyLJ[i];
	      Fz[i] = FzLJ[i];
	     
	      /* vt2 = v(t+dt/2) */
	      vx[i] = vxt2[i] + dt2 * Fx[i] / m;
	      vy[i] = vyt2[i] + dt2 * Fy[i] / m;
	      vz[i] = vzt2[i] + dt2 * Fz[i] / m;
	      vx[i] /= 1.0 + 0.5 * dt * dlns;
	      vy[i] /= 1.0 + 0.5 * dt * dlns;
	      vz[i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[i] * m;
	      FyNose = - dlns * vy[i] * m;
	      FzNose = - dlns * vz[i] * m;
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[i] = FxLJ[i] + FxNose;
	      Fy[i] = FyLJ[i] + FyNose;
	      Fz[i] = FzLJ[i] + FzNose;
	}
    }
  kinet(Nm, vx, vy, vz, 0.0);/* K(t+dt) */

  /* Calculate all components of atomic pressure tensor */
  calcPtensAt(Nm, 0.0);
  press_at = (Patxx + Patyy + Patzz)/3.0;

  /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
  //press_at = calcT1diagAt(Nm, 0.0) + W / Vol;
  /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
  //printf("K exact: %.20f\n", K);

  press = press_at;
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;

  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  s1o2  = s1o1;
  s1o1 = s1t;
  loop(i, 1, Nm)
    {
      vxo2[i] = vxo1[i]; /* vo2(t+dt) = v(t-dt) */
      vyo2[i] = vyo1[i];
      vzo2[i] = vzo1[i];
      vxo1[i] = vxt[i];  /* vo1(t+dt) = v(t) */
      vyo1[i] = vyt[i];
      vzo1[i] = vzt[i];
    }
  /* END OF ITERATION LOOP TO IMPROVE ANDERSEN-NOSE FORCE ESTIMATION */
}


/* ========================= >>> moveb <<< =========================== */
void movebNPT(COORD_TYPE dt, COORD_TYPE m, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, FxAnd, FyAnd, FzAnd, cfAnd, dlns, dlnV,
    dlnVSq; 
  COORD_TYPE s1i, Vol1i;
  COORD_TYPE DT, A, B, DP, dt2; 
  int i, k;
  const int NUMIT = 4;

  /* ******************************************************************* */
  dt2 = dt / 2.0;
  
  loop(i, 1, Nm)
    {
      vxt2[i] = vx[i]; /* v*t2 = v*(t+dt/2) */
      vyt2[i] = vy[i];
      vzt2[i] = vz[i];
      
      FxLJ[i] = Fx[i];
      FyLJ[i] = Fy[i];
      FzLJ[i] = Fz[i];
      
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
	
      vx[i] = 5.0 * vxt[i] / 2.0 - 2.0 * vxo1[i] + 
	    vxo2[i] / 2.0;
      vy[i] = 5.0 * vyt[i] / 2.0 - 2.0 * vyo1[i] + 
	vyo2[i] / 2.0;
      vz[i] = 5.0 * vzt[i] / 2.0 - 2.0 * vzo1[i] + 
	vzo2[i] / 2.0;
    }
  s1i = s1;    /* s1i = s1(t+dt/2) */
  Vol1i = Vol1;/* Vol1i = Vol1(t+dt/2)*/

  /* Initial guess for Vol1 */
  //Vol1 = 13.0*Vol1t/6.0 - 4.0*Vol1o1/3.0 + Vol1o2 / 6.0;/* Vol1(t+dt) */
  
  //Vol1 = 2.0 * Vol1t - Vol1o1; /* Verlet */  
  
  Vol1 = 5.0 * Vol1t / 2.0 - 2.0 * Vol1o1 + Vol1o2 / 2.0;  

  //printf("Vol1: %f Vol1o1: %f Vol1o2: %f Vol1t:%f\n", 
	// Vol1, Vol1o1, Vol1o2, Vol1t);
  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm, vx, vy, vz, Vol1);  /* K(t+dt) */
      //printf("K guess: %.20f\n", K);

      DT = s * (2.0 * K - (3.0 * Nm - 3.0) * Oparams.T) / 
	OprogStatus.Q;
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * Sqr(A) * A / Sqr(s);
      dlns = s1 / s ; 
      s2 = Sqr(s1) / s + DT; /* s2(t+dt) */
      
      /* Calculate pressure, calcT1diagMol is a term proportional to the 
	 translational kinetic energy, see Ferrario and Ryckaert */
      press_at = calcT1diagAt(Nm, Vol1) +  W / Vol; /* press(t+dt) */
      //      printf("k: %d press_at: %f press_m: %f\n", k, press_at, press_m);
      /* Volume acceleration */
      DP = Sqr(s) * (press_at - Oparams.P) / OprogStatus.W;
      B = Vol1i + 0.5 * dt * DP;
      //printf("press_at: %f\n", press_at);
      Vol1 = B / (1.0 - 0.5 * dt * s1 / s); /* Vol1(t+dt) */
      Vol2 = DP + s1 * Vol1 / s;            /* Vol2(t+dt) */
      
      /* Calculate the Andersen forces */
      dlnV   = Vol1 / Vol; 
      dlnVSq = Sqr(dlnV);
      cfAnd  = - 2.0 * dlnVSq / 9.0;
      cfAnd += Vol2 / Vol / 3.0;
      cfAnd += dlns * dlnV / 3.0;
	 
      loop(i, 1, Nm)
	{
	  /* Calculate center of mass posiotin for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	    
	      FxAnd  = cfAnd * rx[i] * m;
	      FyAnd  = cfAnd * ry[i] * m;
	      FzAnd  = cfAnd * rz[i] * m;
	      Fx[i] = FxLJ[i] + FxAnd;
	      Fy[i] = FyLJ[i] + FyAnd;
	      Fz[i] = FzLJ[i] + FzAnd;
	      /* vt2 = v(t+dt/2) */
	      vx[i] = vxt2[i] + dt2 * Fx[i] / m;
	      vy[i] = vyt2[i] + dt2 * Fy[i] / m;
	      vz[i] = vzt2[i] + dt2 * Fz[i] / m;
	      vx[i] /= 1.0 + 0.5 * dt * dlns;
	      vy[i] /= 1.0 + 0.5 * dt * dlns;
	      vz[i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[i] * m;
	      FyNose = - dlns * vy[i] * m;
	      FzNose = - dlns * vz[i] * m;
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[i] = Fx[i] + FxNose;
	      Fy[i] = Fy[i] + FyNose;
	      Fz[i] = Fz[i] + FzNose;
	}
    }

  /* Calculate kinetic energy */
  kinet(Nm, vx, vy, vz, Vol1);/* K(t+dt) */
  //printf("K exact: %.20f\n", K);

  /* Calculate all components of atomic pressure tensor */
  calcPtensAt(Nm, Vol1);
  press_at = (Patxx + Patyy + Patzz) / 3.0;
  //press_at = calcT1diagAt(Nm, Vol1) + W  / Vol;
  
  press = press_at;
  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;


  /* These are the values of velocities at t-dt and at t-2*dt, used to 
     estimate the temperature at time t+dt at begin of this procedure */
  Vol1o2 = Vol1o1;
  s1o2  = s1o1;
  Vol1o1 = Vol1t;
  s1o1 = s1t;
  loop(i, 1, Nm)
    {
      vxo2[i] = vxo1[i]; /* vo2(t+dt) = v(t-dt) */
      vyo2[i] = vyo1[i];
      vzo2[i] = vzo1[i];
      vxo1[i] = vxt[i];  /* vo1(t+dt) = v(t) */
      vyo1[i] = vyt[i];
      vzo1[i] = vzt[i];
    }
  /* END OF ITERATION LOOP TO IMPROVE ANDERSEN-NOSE FORCE ESTIMATION */
}



/* ========================= >>> moveb <<< =========================== */
void moveb(COORD_TYPE dt, COORD_TYPE m, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE dt2;
  int i;
  /* ******************************************************************* */

  dt2 = dt / 2.0;
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;

  /* LOOP OVER ALL MOLECULES */
  loop(i, 1, Nm)
    {
      FxLJ[i] = Fx[i];
      FyLJ[i] = Fy[i];
      FzLJ[i] = Fz[i];
      
      vx[i] = vx[i] + dt2 * Fx[i] / m;
      vy[i] = vy[i] + dt2 * Fy[i] / m;
      vz[i] = vz[i] + dt2 * Fz[i] / m;
      
      /* Kinetic terms of the pressure-tensor */
      T1xy += vx[i] * vy[i] * m; 
      T1yz += vy[i] * vz[i] * m;
      T1zx += vz[i] * vx[i] * m;
      T1xx += vx[i] * vx[i] * m; 
      T1yy += vy[i] * vy[i] * m;
      T1zz += vz[i] * vz[i] * m;
      K = K + Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]);
    }
  pool;
  /* END OF LOOP OVER MOLECULES */
  
  loop(i, 1, Nm)
    {
      vxo2[i] = vxo1[i]; /* vo2(t+dt) = v(t-dt) */
      vyo2[i] = vyo1[i];
      vzo2[i] = vzo1[i];
      vxo1[i] = vxt[i];  /* vo1(t+dt) = v(t) */
      vyo1[i] = vyt[i];
      vzo1[i] = vzt[i];
    }

  
  K  = K * m * 0.5;
         
  ProcSync0();

  /* 
     #ifndef NO_PARALLEL_CODE
     sumAllProcMA(&T1xy, &T1yz, &T1zx, &T1xx, &T1yy, &T1zz, NULL);
     #endif
     
     #ifndef NO_PARALLEL_CODE  
     sumAllProcMA(&K, NULL);
     #endif
  */

  Patxy = T1xy + Wxy;
  Patyz = T1yz + Wyz;
  Patzx = T1zx + Wzx;  
  Patxy /= Vol;
  Patyz /= Vol;
  Patzx /= Vol;
  press_at = (T1xx + T1yy + T1zz)/3.0 + W;
  press_at /= Vol;

  Pxy = Patxy;
  Pyz = Patyz;
  Pzx = Patzx;
  press = press_at;
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
  int i;
  COORD_TYPE cost, costo1, costo2, DRx, DRy, DRz;
  
  L = cbrt(Vol);
  invL = 1.0 / L;
  cost = Vol1 / Vol / 3.0;
  costo1 = Vol1o1 / Vol / 3.0;
  costo2 = Vol1o2 / Vol / 3.0;
  
  /* Reduced particles to first box */
  loop(i, 1, Oparams.parnum)
    {
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * rx[i]);
      DRy = - L * rint(invL * ry[i]);
      DRz = - L * rint(invL * rz[i]);
      
      rx[i] += DRx;
      ry[i] += DRy;
      rz[i] += DRz;
      
      /* The velocity in the Andersen method depend upon the position,
	 so we must add a term also to velocity when switching attention
	 to the image particle. 
	 See Brown and Clarke Mol. Phys., Vol.51, No.5, 1243-1252, 1984 */ 
      vx[i] += cost*DRx;
      vy[i] += cost*DRy;
      vz[i] += cost*DRz;
      vxo1[i] += costo1*DRx;
      vyo1[i] += costo1*DRy;
      vzo1[i] += costo1*DRz;
      vxo2[i] += costo2*DRx;
      vyo2[i] += costo2*DRy;
      vzo2[i] += costo2*DRz;
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
	  sprintf(TXT, "Volume and s within tolerance, exit...\n");
	  mdPrintf(STD, TXT, NULL);
	  // printf("avVol1: %.10f relVol1: %.10f tolVol1: %.10f\n", OprogStatus.avVol1, relVol1, tolVol1);
	  // printf("Vol1: %f DVol1: %f\n", Vol1, DVol1);
	  ENDSIM = 1; /* End simulation */
	}
      
    }
}

/* ============================= >>> chkT <<< ============================== */
void chkT(void)
{
  C_T DT, temp, relT;
  
  if (OprogStatus.tolT > 0.0)
    {
      temp = 2.0 * K / (3.0 * Oparams.parnum - 3.0);
      DT = fabs(temp - Oparams.T);
      relT = DT / Oparams.T;
      if (relT < OprogStatus.tolT)
	{
	  sprintf(TXT, "Temperature within tolerance, exit...\n");
	  mdPrintf(STD, TXT, NULL);
	  // printf("avVol1: %.10f relVol1: %.10f tolVol1: %.10f\n", OprogStatus.avVol1, relVol1, tolVol1);
	  // printf("Vol1: %f DVol1: %f\n", Vol1, DVol1);
	  ENDSIM = 1; /* End simulation */

	}

    }
}

/* ============================= >>> updatePE <<< ========================= */
void updatePE(int Nm)
{
  int iE;
  double ENmin, ENmax;
  
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  iE = (int) ((Vc - ((double) Nm)*ENmin) / 
	      ( ((double) Nm) * ((double) ENmax - ENmin)) * 
	      ((double) PE_POINTS));
  //  sprintf(TXT, "energia: %f iE:%d\n", EE[my_rank], iE);
  //mdPrintf(ALL, TXT, NULL);
  //printf("PEPOINTS: %d Vc:%f\n", PE_POINTS, Vc);
  if ( (iE >= 0) && (iE < PE_POINTS)) 
    {
      ++(OprogStatus.PE[iE]);
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
      temp = 2.0 * K / (3.0 * Oparams.parnum - 3.0);
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
	  sprintf(TXT, "s within tolerance, exit...\n");
	  mdPrintf(STD, TXT, NULL);
	  // printf("avVol1: %.10f relVol1: %.10f tolVol1: %.10f\n", OprogStatus.avVol1, relVol1, tolVol1);
	  // printf("Vol1: %f DVol1: %f\n", Vol1, DVol1);
	  ENDSIM = 1; /* End simulation */
	}
      
    }
}

extern C_T WLJ;

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
  movea(Oparams.steplength, Oparams.m, Oparams.parnum);        

  if (nebrNow)
    {
      //printf("building neighbour listr step: %d\n", Oparams.curStep);
      nebrNow = 0;
      dispHi = 0.0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      if (OprogStatus.noLinkedList)
	{
	  BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	}
      else
	{
	  links(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	  
	  /* Build up neighbour list */  
	  BuildNebrList(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	}
    }

  /* zero forces */
  zeroArrays(Fx, Fy, Fz, Oparams.parnum);
  /* Zero all components of pressure tensor */
  Wxy = 0.0;
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;

if ((OprogStatus.Nose == 1) && (OprogStatus.AntiCry))
  buildMesh(OprogStatus.kmax, OprogStatus.Dk);

if (OprogStatus.AntiCry)
  {
    ACForce(Oparams.parnum, OprogStatus.kmax, OprogStatus.Dk, OprogStatus.alpha,
	    OprogStatus.S0);
    /* calculate forces */
    LJForce(Oparams.parnum, Oparams.epsilon, Oparams.sigma,
	    Oparams.rcut);
    sumFRContribs();
  }
else
  {
    /* calculate forces */
    LJForce(Oparams.parnum, Oparams.epsilon, Oparams.sigma,
	    Oparams.rcut);
    W = WLJ;
  }

  /* correct the coords */
  if (OprogStatus.Nose == 1)
    {  
      /* NPT ensemble */
      movebNPT(Oparams.steplength, Oparams.m, Oparams.parnum);             
    }
  else if (OprogStatus.Nose == 2)
    {
      /* NTV ensemble */
      movebNTV(Oparams.steplength, Oparams.m, Oparams.parnum);             
    }
  else    
    {
      /* NVE ensemble */
      moveb(Oparams.steplength, Oparams.m, Oparams.parnum);             
    }
  /* Calculate the kinetic energy */
  kinet(Oparams.parnum, vx, vy, vz, Vol1);

  checkNebrRebuild();
  if ( (OprogStatus.Nose == 1) || (OprogStatus.Nose == 2))
    {
      //scalCor(Oparams.parnum);
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

  chkT();

  if (  ( (OprogStatus.CMreset > 0) &&
	  // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	  (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	     // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	     (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
      resetCM(Oparams.parnum);

  /* Update the integral of the pressure tensor */
  updateDQ(Oparams.steplength);
  updatePE(Oparams.parnum);

}

