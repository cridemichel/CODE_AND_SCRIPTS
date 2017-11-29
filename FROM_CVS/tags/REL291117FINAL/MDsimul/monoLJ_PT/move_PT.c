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
char  TXT1[MSG_LEN];
extern void links(int ss, int Nm, COORD_TYPE rcut, COORD_TYPE sigma);
extern void BuildNebrListNoLinked(int ss, int Nm, COORD_TYPE rCut, 
				  COORD_TYPE sigma);
extern void BuildNebrList(int ss, int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void checkNebrRebuild(int ss);
extern void LJForce(int ss, int Nm, COORD_TYPE epsilon, 
		    COORD_TYPE sigma, COORD_TYPE rcut);
extern void resetCM(int ss, int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void kinet(int ss, int Nm);

extern double kinetTot(int Nm);

extern void zeroArrays(C_T *arrx, C_T* arry, C_T *arrz, int N);
extern void sumFRContribs(int ss);
extern void buildMesh(C_T kmax, C_T Dk);
extern void ACForce(int ss, int Nm, COORD_TYPE kmax, COORD_TYPE Dk, C_T alpha, C_T S0);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, s1t[MAX_M], Vol1t, L, invL, Ktot;   
COORD_TYPE EE[MAX_M], Vc[MAX_M], V[MAX_M], W[MAX_M], K[MAX_M], T1xx[MAX_M],
  T1yy[MAX_M], T1zz[MAX_M],
  T1xx[MAX_M], T1yy[MAX_M], T1zz[MAX_M], T1xy[MAX_M], 
  T1yz[MAX_M], T1zx[MAX_M], Wxx[MAX_M], Wyy[MAX_M], Wzz[MAX_M],
  Wxy[MAX_M], Wyz[MAX_M], Wzx[MAX_M], Pxx[MAX_M], Pyy[MAX_M],
  Pzz[MAX_M], Pxy[MAX_M], Pyz[MAX_M], Pzx[MAX_M],
  Patxy[MAX_M], Patyz[MAX_M], Patzx[MAX_M], Patxx[MAX_M], Patyy[MAX_M],
  Patzz[MAX_M], T1myz[MAX_M], T1mzx[MAX_M], T1mxx[MAX_M], T1myy[MAX_M],
  T1mzz[MAX_M];  
COORD_TYPE DrSq = 0.0, Mtot;
/* used by linked list routines */
int **head, **list, *map;  /* arrays of integer */
int NCell, mapSize, M;

/* neighbour list method variables */
COORD_TYPE dispHi[MAX_M];
int ***nebrTab, *nebrNow, *nebrTabLen, nebrTabMax;
/* ================================= */


/* ========================== >>> movea <<< =============================== */
void movea(COORD_TYPE dt, COORD_TYPE m, int Nm)
{  
  COORD_TYPE  dt2, dtSq2;
  COORD_TYPE  axi, ayi, azi;
  int i, ss;
 
  L = cbrt(Vol);
  invL = 1.0 / L;

  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;
 
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      //printf("movea: %d\n", Oparams.PTM);
      /* ===== LOOP OVER MOLECULES ===== */
      loop(i, 1, Nm)
	{
	  axi = Fx[ss][i] / m;
	  ayi = Fy[ss][i] / m;
	  azi = Fz[ss][i] / m;
	 
	  rx[ss][i] = rx[ss][i] + dt * vx[ss][i] + dtSq2 * axi;
	  ry[ss][i] = ry[ss][i] + dt * vy[ss][i] + dtSq2 * ayi;
	  rz[ss][i] = rz[ss][i] + dt * vz[ss][i] + dtSq2 * azi;
	  
	  vxt[ss][i] = vx[ss][i]; /* v(t) */
	  vyt[ss][i] = vy[ss][i];
	  vzt[ss][i] = vz[ss][i];
	  
	  vx[ss][i] = vx[ss][i] + dt2 * axi;
	  vy[ss][i] = vy[ss][i] + dt2 * ayi;
	  vz[ss][i] = vz[ss][i] + dt2 * azi;
	}
      pool;
    }
}

/* ========================== >>> movea <<< =============================== */
void moveaNTV(COORD_TYPE dt, COORD_TYPE m, int Nm)
{  
  COORD_TYPE  dt2, dtSq2;
  COORD_TYPE  axi, ayi, azi;
  int i, ss;
 
  L = cbrt(Vol);
  invL = 1.0 / L;

  dt2    = dt / 2.0;
  dtSq2  = dt * dt2;
 
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      //printf("movea: %d\n", Oparams.PTM);
      /* ===== LOOP OVER MOLECULES ===== */
      loop(i, 1, Nm)
	{
	  axi = Fx[ss][i] / m;
	  ayi = Fy[ss][i] / m;
	  azi = Fz[ss][i] / m;
	 
	  rx[ss][i] = rx[ss][i] + dt * vx[ss][i] + dtSq2 * axi;
	  ry[ss][i] = ry[ss][i] + dt * vy[ss][i] + dtSq2 * ayi;
	  rz[ss][i] = rz[ss][i] + dt * vz[ss][i] + dtSq2 * azi;
	  
	  vxt[ss][i] = vx[ss][i]; /* v(t) */
	  vyt[ss][i] = vy[ss][i];
	  vzt[ss][i] = vz[ss][i];
	  
	  vx[ss][i] = vx[ss][i] + dt2 * axi;
	  vy[ss][i] = vy[ss][i] + dt2 * ayi;
	  vz[ss][i] = vz[ss][i] + dt2 * azi;
	  	 

	}
      pool;
      /* END OF LOOP OVER MOLECULES */
      s[ss] = s[ss] + dt * s1[ss] + dtSq2 * s2[ss];
      s1t[ss] = s1[ss];       /* s1(t) */
      s1[ss] = s1[ss] + dt2 * s2[ss];
    }
}

/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagAt(int ss, int Nm, COORD_TYPE VOL1)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of atomic pressure tensor */
  int i;
  COORD_TYPE kin, m;
  
  m = Oparams.m;
  kin = 0.0;
  loop(i, 1, Nm)
    {
	  kin +=  m*(Sqr(vx[ss][i]) + 
	    Sqr(vy[ss][i]) + 
	    Sqr(vz[ss][i]));
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
void calcPtensAt(int ss, int Nm)
{
  /* Calculate all components of atomic pressure tensor */
  int i;
  COORD_TYPE px, py, pz, m;
  
  m = Oparams.m;
  T1xy[ss] = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz[ss] = 0.0;
  T1zx[ss] = 0.0;
  T1xx[ss] = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy[ss] = 0.0;
  T1zz[ss] = 0.0;

  loop(i, 1, Nm)
    {
      px = vx[ss][i];
      py = vy[ss][i];
      pz = vz[ss][i];
      /* Kinetic component of pressure tensor (all terms) */
      T1xy[ss] += px * py * m; 
      T1yz[ss] += py * pz * m;
      T1zx[ss] += pz * px * m ;
      T1xx[ss] += px * px * m;
      T1yy[ss] += py * py * m;
      T1zz[ss] += pz * pz * m;
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
  Patxy[ss] = T1xy[ss] + Wxy[ss];
  Patyz[ss] = T1yz[ss] + Wyz[ss];
  Patzx[ss] = T1zx[ss] + Wzx[ss];  
  Patxx[ss] = T1xx[ss] + Wxx[ss];
  Patyy[ss] = T1yy[ss] + Wyy[ss];
  Patzz[ss] = T1zz[ss] + Wzz[ss];
  //printf("Wxx:%f WCxx:%f Wyy: %f WCyy: %f Wzz: %f WCzz: %f\n",
	//Wxx, WCxx, Wyy, WCyy, Wzz, WCzz);
  Patxx[ss] /= Vol;
  Patyy[ss] /= Vol;
  Patzz[ss] /= Vol;
  Patxy[ss] /= Vol;
  Patyz[ss] /= Vol;
  Patzx[ss] /= Vol;
}

void calcKtot();

/* ========================= >>> moveb <<< =========================== */
void movebNTV(COORD_TYPE dt, COORD_TYPE m, int Nm)
{
  /* *******************************************************************
     ** SECOND PART OF VELOCITY VERLET WITH CONSTRAINTS               **
     ******************************************************************* */
  COORD_TYPE FxNose, FyNose, FzNose, dlns; 
  COORD_TYPE s1i[MAX_M];
  COORD_TYPE DT, A, dt2; 
  int i, k, ss;
  const int NUMIT = 4;

  /* ******************************************************************* */
  dt2 = dt * 0.5;

  for (ss = 0; ss < Oparams.PTM; ss++)
    {    
      loop(i, 1, Nm)
	{
      
	  vxt2[ss][i] = vx[ss][i]; /* v*t2 = v*(t+dt/2) */
	  vyt2[ss][i] = vy[ss][i];
	  vzt2[ss][i] = vz[ss][i];
	  
	  FxLJ[ss][i] = Fx[ss][i];
	  FyLJ[ss][i] = Fy[ss][i];
	  FzLJ[ss][i] = Fz[ss][i];
	  
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
	  
	  vx[ss][i] = 5.0 * vxt[ss][i] / 2.0 - 2.0 * vxo1[ss][i] + 
	    vxo2[ss][i] / 2.0;
	  vy[ss][i] = 5.0 * vyt[ss][i] / 2.0 - 2.0 * vyo1[ss][i] + 
	    vyo2[ss][i] / 2.0;
	  vz[ss][i] = 5.0 * vzt[ss][i] / 2.0 - 2.0 * vzo1[ss][i] + 
	    vzo2[ss][i] / 2.0;
	}

      s1i[ss] = s1[ss];    /* s1i = s1(t+dt/2) */
    }
 
  //printf("Vol1: %f Vol1o1: %f Vol1o2: %f Vol1t:%f\n", 
	// Vol1, Vol1o1, Vol1o2, Vol1t);
  //Ktot = kinetTot(Nm);
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
	{
	  kinet(ss, Nm);
	  /* Ora il numero di gradi di liberta' e' numOfProcs(=M nella'art.
	     di Kob) superiore */
	  DT = s[ss] * (2.0 * K[ss] - (3.0 * Nm - 3.0) 
	  		* Oparams.T) / OprogStatus.Q;
	  //DT = s[ss] * (2.0 * Ktot - Oparams.PTM * (3.0 * Nm - 3.0) 
	  //	* Oparams.T) / OprogStatus.Q;
	  

	  A = s1i[ss] + 0.5 * dt * DT;
	  /* s1(t+dt) */
	  s1[ss] = A + 0.5 * Sqr(A) * dt / s[ss] + Sqr(dt) * 0.5 * 
	    Sqr(A) * A / Sqr(s[ss]);
	  dlns = s1[ss] / s[ss] ; 
	  s2[ss] = Sqr(s1[ss]) / s[ss] + DT; /* s2(t+dt) */
	  
	  loop(i, 1, Nm)
	    {
	      /* Calculate center of mass posiotin for molecule to which
		 belong atom (a, i) */
	      /* The forces due to interactions don't change inside this 
		 loop */
	      Fx[ss][i] = FxLJ[ss][i];
	      Fy[ss][i] = FyLJ[ss][i];
	      Fz[ss][i] = FzLJ[ss][i];
	      
	      /* vt2 = v(t+dt/2) */
	      vx[ss][i] = vxt2[ss][i] + dt2 * Fx[ss][i] / m;
	      vy[ss][i] = vyt2[ss][i] + dt2 * Fy[ss][i] / m;
	      vz[ss][i] = vzt2[ss][i] + dt2 * Fz[ss][i] / m;
	      vx[ss][i] /= 1.0 + 0.5 * dt * dlns;
	      vy[ss][i] /= 1.0 + 0.5 * dt * dlns;
	      vz[ss][i] /= 1.0 + 0.5 * dt * dlns;
	      /* Velocity calculated with contribution of Andersen and Nose 
		 Forces */
	      FxNose = - dlns * vx[ss][i] * m;
	      FyNose = - dlns * vy[ss][i] * m;
	      FzNose = - dlns * vz[ss][i] * m;
	      /* F = Flennard-jones + Fnose + Fandersen */
	      Fx[ss][i] = FxLJ[ss][i] + FxNose;
	      Fy[ss][i] = FyLJ[ss][i] + FyNose;
	      Fz[ss][i] = FzLJ[ss][i] + FzNose;
	    }
	}
    }
  

  Ktot = kinetTot(Nm);/* K(t+dt) */

  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      /* Calculate all components of atomic pressure tensor */
      calcPtensAt(ss, Nm);
      press_at[ss] = (Patxx[ss] + Patyy[ss] + Patzz[ss])/3.0;
      
      /* NOTE: Calculate Pressure (Eliminate in the future NOT USED only output)*/
      //press_at = calcT1diagAt(Nm, 0.0) + W / Vol;
      /* NOTE: Nm*5/3 = (6Nm - Nm) / 3.0 */ 
      //printf("K exact: %.20f\n", K);
      
      press[ss] = press_at[ss];
      Pxy[ss] = Patxy[ss];
      Pyz[ss] = Patyz[ss];
      Pzx[ss] = Patzx[ss];
      
      /* These are the values of velocities at t-dt and at t-2*dt, used to 
	 estimate the temperature at time t+dt at begin of this procedure */
    
    
      loop(i, 1, Nm)
	{
	  vxo2[ss][i] = vxo1[ss][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[ss][i] = vyo1[ss][i];
	  vzo2[ss][i] = vzo1[ss][i];
	  vxo1[ss][i] = vxt[ss][i];  /* vo1(t+dt) = v(t) */
	  vyo1[ss][i] = vyt[ss][i];
	  vzo1[ss][i] = vzt[ss][i];
	}
      
      s1o2[ss]  = s1o1[ss];
      s1o1[ss] = s1t[ss];

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
  int i, ss;
  /* ******************************************************************* */

  dt2 = dt / 2.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      K[ss] = 0.0;
      T1xy[ss] = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
      T1yz[ss] = 0.0;
      T1zx[ss] = 0.0;
      T1xx[ss] = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
      T1yy[ss] = 0.0;
      T1zz[ss] = 0.0;
      
      /* LOOP OVER ALL MOLECULES */
      loop(i, 1, Nm)
	{
	  FxLJ[ss][i] = Fx[ss][i];
	  FyLJ[ss][i] = Fy[ss][i];
	  FzLJ[ss][i] = Fz[ss][i];
	    
	  vx[ss][i] = vx[ss][i] + dt2 * Fx[ss][i] / m;
	  vy[ss][i] = vy[ss][i] + dt2 * Fy[ss][i] / m;
	  vz[ss][i] = vz[ss][i] + dt2 * Fz[ss][i] / m;
	  
	  /* Kinetic terms of the pressure-tensor */
	  T1xy[ss] += vx[ss][i] * vy[ss][i] * m; 
	  T1yz[ss] += vy[ss][i] * vz[ss][i] * m;
	  T1zx[ss] += vz[ss][i] * vx[ss][i] * m;
	  T1xx[ss] += vx[ss][i] * vx[ss][i] * m; 
	  T1yy[ss] += vy[ss][i] * vy[ss][i] * m;
	  T1zz[ss] += vz[ss][i] * vz[ss][i] * m;
	  K[ss] = K[ss] + Sqr(vx[ss][i]) + Sqr(vy[ss][i]) + Sqr(vz[ss][i]);
	}
      pool;
      /* END OF LOOP OVER MOLECULES */
      
      loop(i, 1, Nm)
	{
	  vxo2[ss][i] = vxo1[ss][i]; /* vo2(t+dt) = v(t-dt) */
	  vyo2[ss][i] = vyo1[ss][i];
	  vzo2[ss][i] = vzo1[ss][i];
	  vxo1[ss][i] = vxt[ss][i];  /* vo1(t+dt) = v(t) */
	  vyo1[ss][i] = vyt[ss][i];
	  vzo1[ss][i] = vzt[ss][i];
	}

      
      K[ss]  = K[ss] * m * 0.5;
      //printf("k[%d]:%f\n", ss, K[ss]);
      //  ProcSync0();
      
      /* 
	 #ifndef NO_PARALLEL_CODE
	 sumAllProcMA(&T1xy, &T1yz, &T1zx, &T1xx, &T1yy, &T1zz, NULL);
	 #endif
	 
	 #ifndef NO_PARALLEL_CODE  
	 sumAllProcMA(&K, NULL);
	 #endif
      */

      Patxy[ss] = T1xy[ss] + Wxy[ss];
      Patyz[ss] = T1yz[ss] + Wyz[ss];
      Patzx[ss] = T1zx[ss] + Wzx[ss];  
      Patxy[ss] /= Vol;
      Patyz[ss] /= Vol;
      Patzx[ss] /= Vol;
      press_at[ss] = (T1xx[ss] + T1yy[ss] + T1zz[ss])/3.0 + W[ss];
      press_at[ss] /= Vol;
      
      Pxy[ss] = Patxy[ss];
      Pyz[ss] = Patyz[ss];
      Pzx[ss] = Patzx[ss];
      press[ss] = press_at[ss];
    }
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

const COORD_TYPE bc1 = 14.0/45.0, bc2 = 64.0/45.0, bc3 = 24.0/45.0;
  
/* =========================== >>> BodeTerm <<< ============================*/
COORD_TYPE BodeTerm(COORD_TYPE dt, COORD_TYPE* fi)
{
  return dt * (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	       bc1 * fi[4]);
}

/* =========================== >>> BodeTerm <<< ============================*/
COORD_TYPE BodeTermB(COORD_TYPE* fi)
{
  return (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	  bc1 * fi[4]);
}


/* ============================ >>> updateQ <<< =========================== */
void updateDQ(int ss, COORD_TYPE dt)
{
  /* Intgerate the pressure tensor using Bode's rule
     WARNING: the minimum effective update  for DQ(t) is 4 steps, because
     Bode rule needs 5 points at least */ 
  int curStep = Oparams.curStep, i1;
  COORD_TYPE *PxyArr, *PyzArr, *PzxArr;  
  int* lt;
  PxyArr = OprogStatus.PxyArr[ss],
  PyzArr = OprogStatus.PyzArr[ss],
  PzxArr = OprogStatus.PzxArr[ss];

  lt = Oparams.lambdat;

  //printf(" curStep: %d\n", curStep);
  /* Store points to use later for integration */
  if (curStep == 1)
    {
      /* WARNING: we assume here that the simulation starts from 1, 
	 that is the first step must be 1 <<<========= !!!!!!!!!! */
      /* First time only HERE */
      PxyArr[4] = Pxy[ss];
      PyzArr[4] = Pyz[ss];
      PzxArr[4] = Pzx[ss];
      return;
    }
  
  i1 = (curStep - 1) % 4;
  
  if (i1 != 0)
    {
      PxyArr[i1] = Pxy[ss];
      PyzArr[i1] = Pyz[ss];
      PzxArr[i1] = Pzx[ss];
    }
  
  if (i1 == 0)
    {
      /* Last point of previuos set is the first point of the actual set */
      PxyArr[0] = PxyArr[4];
      PyzArr[0] = PyzArr[4];
      PzxArr[0] = PzxArr[4];
      /* 5th point is the actual value of pressure tensor */
      PxyArr[4] = Pxy[ss];
      PyzArr[4] = Pyz[ss];
      PzxArr[4] = Pzx[ss];
      /* Here we use molecular pressure tensor */
      OprogStatus.DQxy[lt[ss]] += BodeTerm(dt, PxyArr);
      OprogStatus.DQyz[lt[ss]] += BodeTerm(dt, PyzArr);
      OprogStatus.DQzx[lt[ss]] += BodeTerm(dt, PzxArr);
    }
  /*         _ t
	    |
    DQab =  |  Pab dt
            |
	   - 0
	   
    con a, b = x,y,z 
  */
}

/* =========================== >>> updateV <<< ========================== */
void updateV(int Nm, C_T dt)
{
  int i, *lt;
  int ss;
  /*
             _ t
	    |
    V_i  =  |  v_i dt
            |
       	   - 0
	
  */
  lt = Oparams.lambdat;
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      for(i = 0; i < Nm; i++)
	{
	  /* lt[ss] e' la temperatura a cui si trova la replica
	     ss */
	  OprogStatus.sumVx[lt[ss]][i] += vx[ss][i]*dt;
	  OprogStatus.sumVy[lt[ss]][i] += vy[ss][i]*dt;
	  OprogStatus.sumVz[lt[ss]][i] += vz[ss][i]*dt;
	}
    }
}


/* ============================ >>> chkPTEq <<< =============================*/
void chkPTEq(void)
{
  double PEij[MAX_M][PE_POINTS];
  double lambdai, lambdaj, norm, maxE;
  double ENmin, ENmax, beta0, dblE, refPE, actPE, diffTot, maxPE;
  double bot[5];
  const int refj = 0;
  int maxiE, iE, ss, Nm, nDiff, i;
    
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  Nm = Oparams.parnum;
  beta0 = 1.0/Oparams.T;

  /* Questa e' il livello al di sotto del quale i valori delle
     distribuzioni di energia vengono tagliati */
  /* Se non c'e' abbastanza statistica esci */
  maxPE = 0.0;
  maxiE = 0;
  /* Determina il massimo della curva di riferimento e usa il valore 
     dell'energia dell massimo come riferimento per l'energia */
  for (iE = 0; iE < PE_POINTS; iE++)
    {
      if (OprogStatus.PE[refj][iE] > maxPE)
	{
	  maxiE = iE;
	  maxPE = OprogStatus.PE[refj][iE];
	}
    }
  maxE = ((double)Nm) * (ENmin + ((double)maxiE)
			 * (ENmax - ENmin) / ((double)PE_POINTS));
	
  lambdaj = Oparams.lambda0[refj];
  for (ss = refj; ss < refj + OprogStatus.NUM_PEij; ss++)
   {
     /* Usa per il check dell'equilibratura i quattro 
	sistemi a piu' bassa temperatura */
     lambdai = Oparams.lambda0[ss];
     /* Fattore di normalizzazione */
     for (iE = 0; iE < PE_POINTS; iE++)
       {
	 /* Calcola valori non senza normalizzazione */
	 dblE = ((double)Nm) * (ENmin + ((double)iE) 
				* (ENmax - ENmin) / ((double)PE_POINTS));
	 dblE -= maxE;
	 PEij[ss][iE] = ((double) OprogStatus.PE[ss][iE])*
	   exp(beta0*(lambdai-lambdaj)*dblE);
       }
     /* Calcola il fattore di normalizzazione integrando con il metodo
	di Bode (ved. Numerical Recipe) */
     norm = 0.0;
     for (iE = 0; iE < PE_POINTS - 4; iE = iE + 4)
       {
	 for (i = 0; i < 5; i++)
	   {
	     /* Notare che l'energia e' opprtunamente shiftata per 
		avere nnumeri piu' piccoli a causa dell'esponeziale 
		(si rischia l'overflow!!!) */
	     bot[i] = PEij[ss][iE+i];
	   }
	 /*
	   dblE = ((double)Nm) * (EN_MIN + ((double)iE+1 - maxiE) 
	   * (EN_MAX - EN_MIN) / 
	   ((double)PE_POINTS));
	   PEij[ss][iE+1] = ((double) OprogStatus.PE[ss][iE+1])*
	   exp(beta0*(lambdai-lambdaj)*dblE);
	   norm += PEij[ss][iE] + PEij[ss][iE+1]; */
	 norm += BodeTermB(bot);
       }
     /* moltiplica per il deltaE */
     norm *= Nm * (ENmax - ENmin) / PE_POINTS;;

     for (iE = 0; iE < PE_POINTS - 1; iE++)
       {
	 if (norm != 0)
	 {
	   /* Normalizza! */
	   PEijM[ss*PE_POINTS + iE] = PEij[ss][iE] / norm;
	   PEij[ss][iE] /= norm;
	 }
       }
   }

  
  /* E ora confrontiamo le varie distribuzioni riscalate ottenute
     (ved. cond-mat/00001042 Kob e Yamamoto) */
  nDiff = 0;
  diffTot = 0.0;
  for (iE = 0; iE < PE_POINTS; iE++)
    {  
      refPE = PEij[refj][iE];
      for (ss = refj ; ss < refj + OprogStatus.NUM_PEij ; ss++)
	{
	  if (ss == refj)
	    continue;
	  actPE = PEij[ss][iE];
	  //if ( (refPE / maxPE[refj]) > minPerc && 
	  //   (actPE / maxPE[ss]) > minPerc )
	  if ((refPE > OprogStatus.sogliaPEij) && 
	      (actPE > OprogStatus.sogliaPEij))
	    {
	      /* Calcola la somma di tutte le differenze in valore assoluto 
		 delle curve riscalate in scala semilogaritmica (
		 ved. cond-mat/0001042) e se la media di tali differenze e' minore
		 di una soglia prestabilita allora il sistema e' equilibrato */
	      //printf("actPE: %f refPE: %f logaPE: %f logrefPE: %f\n"
	      //     , actPE, refPE, log10(actPE),log10(refPE));
	      if (OprogStatus.linearDiff == 0)
		diffTot += fabs(log10(refPE)-log10(actPE));
	      else
		diffTot += fabs(refPE-actPE);
	      nDiff++;
	    }
	}
    }

  diffMedia = diffTot / ((double)nDiff);
  if (isnan(diffMedia))
    diffMedia = -1.0;
  
  //printf("diffTot: %f nDiff: %d  media = %f\n", 
  //diffTot, nDiff, diffMedia);
  if (nDiff > 0 && diffMedia < SOGLIA_PTEQ)
    {
      OprogStatus.PTequilibrated = 1;
      //printf("EQUILIBRATO!!!\n");
    }
  else
    {
      OprogStatus.PTequilibrated = 0;
    }
}


/* ============================ >>> chkPTEq <<< =============================*/
void chkPTEqB(void)
{
  /* Questa routine differisce dalla precedente poiche' ogni curva e'
     confrontata con la curva riscalata successiva ossia vengono calcolate
     PE_i(E, beta0*lambda_i) e PE_i(beta0*lambda_i+1) dove i = 1...M 
     (ved. cond-mat/0001042) */
  double PEij[MAX_M][PE_POINTS];
  double lambdai, lambdaj, norm, maxE;
  double ENmin, ENmax, beta0, dblE, refPE, actPE, diffTot, maxPE;
  double bot[5];
  int jj, maxiE, iE, ss, Nm, nDiff, nstp, i;

  //  printf("eccomi qua!!!\n");
  ENmax = OprogStatus.ENmax;
  ENmin = OprogStatus.ENmin;
  Nm = Oparams.parnum;
  //  maxE = -6.0 * ((double)Nm);
  beta0 = 1.0/Oparams.T;
  //printf("chPTEqB!!!!!!!!!!!1\n");
  /* Questa e' il livello al di sotto del quale i valori delle
     distribuzioni di energia vengono tagliati */
  /* Se non c'e' abbastanza statistica esci */
  nstp = (PE_POINTS / Oparams.PTM);
  nstp *= rint(sqrt(nstp));
  //printf("nstp: %d\n", nstp);
  //if (Oparams.curStep < nstp ) return;
  /* Determina il massimo della curva di riferimento e usa il valore 
     dell'energia dell massimo come riferimento per l'energia */

  for (ss = 0; ss < Oparams.PTM - 1; ss++)
   {
     /* Riscala due curve attigue */
     lambdaj = Oparams.lambda0[ss+1];
     /* Usa per il check dell'equilibratura i quattro 
	sistemi a piu' bassa temperatura */
     lambdai = Oparams.lambda0[ss];
     /* Fattore di normalizzazione */
     
     maxPE = 0.0;
     maxiE = 0;
     for (iE = 0; iE < PE_POINTS; iE++)
       {
	 if (OprogStatus.PE[ss][iE] > maxPE)
	   {
	     maxiE = iE;
	     maxPE = OprogStatus.PE[ss][iE];
	   }
       }
     maxE = ((double)Nm) * (ENmin + ((double)maxiE)
			    * (ENmax - ENmin) / ((double)PE_POINTS));
       
     for (iE = 0; iE < PE_POINTS; iE++)
       {
	 /* Calcola valori non senza normalizzazione */
	 dblE = ((double)Nm) * (ENmin + ((double)iE) 
				* (ENmax - ENmin) / ((double)PE_POINTS));
	 dblE -= maxE;
	 PEij[2*ss][iE] = ((double) OprogStatus.PE[ss][iE])*
	   exp(beta0*(lambdai-lambdaj)*dblE);
	 PEij[2*ss+1][iE] = ((double) OprogStatus.PE[ss+1][iE]);
	 //printf("PEij: %f\n", PEij[2*ss][iE]);
       }
     /* Calcola il fattore di normalizzazione integrando con il metodo
	di Bode (ved. Numerical Recipe) */
     for (jj = 0; jj < 2; jj++) 
     {
       norm = 0.0;
       for (iE = 0; iE < PE_POINTS - 4; iE = iE + 4)
	 {
	   for (i = 0; i < 5; i++)
	     {
	       /* Notare che l'energia e' opprtunamente shiftata per 
		  avere nnumeri piu' piccoli a causa dell'esponeziale 
		  (si rischia l'overflow!!!) */
	       bot[i] = PEij[2*ss+jj][iE+i];
	     }
	   /*
	     dblE = ((double)Nm) * (EN_MIN + ((double)iE+1 - maxiE) 
	     * (EN_MAX - EN_MIN) / 
	     ((double)PE_POINTS));
	     PEij[ss][iE+1] = ((double) OprogStatus.PE[ss][iE+1])*
	     exp(beta0*(lambdai-lambdaj)*dblE);
	     norm += PEij[ss][iE] + PEij[ss][iE+1]; */
	   norm += BodeTermB(bot);
	 }
       /* moltiplica per il deltaE */
       norm *= ((double)Nm) * (ENmax - ENmin) / ((double)PE_POINTS);
       
       for (iE = 0; iE < PE_POINTS - 1; iE++)
	 {
	   if (norm != 0)
	     {
	       /* Normalizza! */
	       PEijM[(2*ss+jj)*PE_POINTS + iE] = PEij[2*ss+jj][iE] / norm;
	       //printf("PeijM: %f\n", PEijM[(2*ss+jj)*PE_POINTS + iE]);
	       PEij[2*ss+jj][iE] /= norm;
	     }
	 }
     }
   }
  
  /* E ora confrontiamo le varie distribuzioni riscalate ottenute
     (ved. cond-mat/00001042 Kob e Yamamoto) */
  nDiff = 0;
  diffTot = 0.0;
  for (iE = 0; iE < PE_POINTS; iE++)
    {  
      for (ss = 0 ; ss < Oparams.PTM - 1 ; ss++)
	{
	  /* confronta i e j=i+1 */
	  refPE = PEij[2*ss+1][iE];
	  actPE = PEij[2*ss][iE];
	  //if ( (refPE / maxPE[refj]) > minPerc && 
	  //   (actPE / maxPE[ss]) > minPerc )
	  if ((refPE > OprogStatus.sogliaPEij) && 
	      (actPE > OprogStatus.sogliaPEij))
	    {
	      /* Calcola la somma di tutte le differenze in valore assoluto 
		 delle curve riscalate in scala semilogaritmica (
		 ved. cond-mat/0001042) e se la media di tali differenze e' minore
		 di una soglia prestabilita allora il sistema e' equilibrato */
	      //printf("actPE: %f refPE: %f logaPE: %f logrefPE: %f\n"
	      //     , actPE, refPE, log10(actPE),log10(refPE));
 	      if (OprogStatus.linearDiff == 0)
		diffTot += fabs(log10(refPE)-log10(actPE));
	      else
		diffTot += fabs(refPE-actPE);
	      nDiff++;
	    }
	}
    }

  diffMedia = diffTot / ((double)nDiff);
  if (isnan(diffMedia))
    diffMedia = -1.0;
  //printf("diffTot: %f nDiff: %d  media = %f\n", 
  //diffTot, nDiff, diffMedia);
  if (nDiff > 0 && diffMedia < SOGLIA_PTEQ)
    {
      OprogStatus.PTequilibrated = 1;
      //printf("EQUILIBRATO!!!\n");
    }
  else
    {
      OprogStatus.PTequilibrated = 0;
    }
}

/* ============================ >>> saveSnaps <<< ========================= */
void saveSnapsPT()
{
  int ss;
  if (OprogStatus.PTequilibrated)
    {
      for (ss = 0; ss < Oparams.PTM; ss++)
	{
	  
	  //savedSteps[ss]++
	}

    }
}

extern double ranf(void);

/* ============================= >>> chkPT <<< ============================= */
/* ADDED 15/09/2000: Parallel Tempering */
void PTexchange()
{
  double soglia, w;//, Ffact0, Ffact1;
  double neighE, myE, Delta;
  int  myLambda, neighLambda;
  int primoSS, ss;
  int *lt;
  double *l0;
  /* da il sistema che corrisponde ad un certo Lambda 
   e' la funzione inversa di lambda[t] ossia:
   lambda[Mt[ss]] = ss 
  */
  int Mt[MAX_M];
  /* Oparams.lambdat e' un array di puntatori ai lambda iniziali, i
     n sostanza di ogni sottosistema dice qual'e' il lambda che gli 
     appartiene ad un dato istante */   
  lt = Oparams.lambdat;
  l0 = Oparams.lambda0;
 
  if (OprogStatus.PT == 0) 
    return;

  if ( (OprogStatus.RW == 0) && 
       (Oparams.lambdat[OprogStatus.refReplica] == Oparams.PTM) )
    {
      OprogStatus.RW = 1;
    }
  
  if ( (OprogStatus.RW == 1) && 
       (Oparams.lambdat[OprogStatus.refReplica] == 0) )
    {
      OprogStatus.RW = 2;
      /* Quest vuol dire che il sistema che all'inizio era a temperatura
	 piu' bassa ha completato il suo random walk nello spazio delle
	 temperature, ossia e' arrivato alla temperatura piu' alta
	 ed e' tornato indietro alla piu' bassa.
	 Se all'inizio si pone RW = 2 il programma comincia da subito
	 a calcolare le funzioni di distribuzione P(E,beta)
      */
    }
  
  /* Controlla se e' ora di fare gli scambi ed eventualmente li fa */
  if ( Oparams.curStep % OprogStatus.DtPT == 0) 
    {
      for (ss = 0; ss < Oparams.PTM; ss++)
	{
	  Mt[lt[ss]] = ss;
	}
      /*
	for (ss = 0; ss < Oparams.PTM; ss++)
	{
	printf("Mt[%d]: %d lt[%d]:%d\n", ss, Mt[ss], ss, lt[ss]);
	}
      */
      
      /* Raccoglie tutti i lambda e tutte le energie potenziali */
      if ( (Oparams.curStep % (2*OprogStatus.DtPT)) == 0) 
	{
	  /* scambia gli i = 0, 2, 4 con gli i+1 = 1, 3, 5 */
	  primoSS = 0;
	}
      else
	{
	  /* scambia gli i = 1, 3, 5 con gli i+1 = 2, 4, 6 */
	  primoSS = 1;
	}
      //printf("primoSS: %d\n", primoSS);
      if (OprogStatus.srescale != 0)
	{

	  for (ss = 0; ss < Oparams.PTM; ss++)
	    {
	      /* 03/10/2000 ADDED s reset in chkPT */
	      /* Resetta tutte le variabili s */
	      //s1[ss] /= s[ss];
	      /* forse cio' non serve !!!!!!!!!!!!!*/
	      //	  break;
	      
	      s[ss] = 1.0;
	      s1[ss] = 0.0;
	      s2[ss] = 0.0;
	      s1o1[ss] = 0.0;
	      s1o2[ss] = 0.0;
	    }
	}
      strcpy(TXT1, "Prob. Sc.:");
      for (ss = primoSS; ss + 1 < Oparams.PTM; ss+=2)
	{
	  soglia = ranf();
	  OprogStatus.attempted[ss]++;	
	  //printf("soglia: %f\n", soglia);
	  //printf("soglia[%d]: %f\n", ss, soglia);
	  myLambda = ss;
	  myE = EE[Mt[ss]];
	  
	  neighLambda = ss+1;
	  neighE = EE[Mt[ss+1]];
	  //printf("EEd[%d]: %f EE[%d]: %f diff: %f\n",
	  // ss, EE[ss], ss+1, EE[ss+1], EE[ss+1] - EE[ss]);
	  Delta = 1.0/Oparams.T * 
	    (l0[neighLambda] - l0[myLambda]) * (myE - neighE);
	  if (Delta <= 0) 
	    {
	      w = 1;
	    }
	  else 
	    {
	      w = exp(-Delta);
	    }
	  //printf("w: %f\n", w);
	  if (w >= soglia)
	    {
	      /* Questa e' la probabilita' di scambio
		 tra ss e ss+1 */
	      OprogStatus.scambi[ss]++;
	      /* Scambia anzitutto i lambda */
	      /*
	      Ffact0 = l0[ss+1] / l0[ss];
	      Ffact1 = l0[ss] / l0[ss+1];
	      for (i = 0; i < Oparams.parnum; i++)
		{
		  Fx[Mt[ss]][i] *= Ffact0;
		  Fy[Mt[ss]][i] *= Ffact0;
		  Fz[Mt[ss]][i] *= Ffact0;
		  Fx[Mt[ss+1]][i] *= Ffact1;
		  Fy[Mt[ss+1]][i] *= Ffact1;
		  Fz[Mt[ss+1]][i] *= Ffact1;
		  }*/
	      lt[Mt[ss]] = ss+1;
	      lt[Mt[ss + 1]] = myLambda;
	    }
	  
	  sprintf(TXT, "[%d<->%d:%.4f]", ss, ss+1,
		  ((double) OprogStatus.scambi[ss]) / 
		  ((double)OprogStatus.attempted[ss]) );
	  strcat(TXT1, TXT);
	}
      strcat(TXT1,"\n");
      /* Stampa tutte le probabilita' di scambio */
      mdPrintf(STD, TXT1, NULL);
    }
}
/* ============================= >>> updatePE <<< ========================= */
void updatePE(int Nm)
{
  int iE, *lt;
  double ENmin, ENmax;
  int ss;
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  lt = Oparams.lambdat;
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      //printf("lambda[%d]:%f", ss, Oparams.lambda0[Oparams.lambdat[ss]]);
      iE = (int) ((EE[ss] - ((double) Nm)*ENmin) / 
		  ( ((double) Nm) * ((double) ENmax - ENmin)) * 
		  ((double) PE_POINTS));
      //printf("EE[%d]: %f iE:%d\n", ss, EE[ss], iE);
      if ( (iE >= 0) && (iE < PE_POINTS)) 
	{
	  ++(OprogStatus.PE[lt[ss]][iE]);
	  //printf("step: %d opsPE[%d][%d]: %d\n", 
	  // Oparams.curStep, lt[ss], iE, OprogStatus.PE[lt[ss]][iE]);
	}
    }
}

/* ============================= >>> chkT <<< ============================== */
void chkT(void)
{
  C_T DT, temp, relT;
  double Ktot;

  if (OprogStatus.tolT > 0.0)
    {
      Ktot = kinetTot(Oparams.parnum);
      /* Calcola la temperatura "totale" */
      temp = 2.0 * Ktot / Oparams.PTM * (5.0 * Oparams.parnum - 3.0);
      DT = fabs(temp - Oparams.T);
      relT = DT / Oparams.T;
      if (relT < OprogStatus.tolT)
	{
	  sprintf(TXT, "Global temperature within tolerance, exit...\n");
	  mdPrintf(STD, TXT, NULL);
	  // printf("avVol1: %.10f relVol1: %.10f tolVol1: %.10f\n", OprogStatus.avVol1, relVol1, tolVol1);
	  // printf("Vol1: %f DVol1: %f\n", Vol1, DVol1);
	  ENDSIM = 1; /* End simulation */

	}

    }
}

/* =============================== >>> chs <<< =========================== */
/* RIVEDERE !!!!!!!!!!!!!!!!!!!!!!!! */
void chks(void)
{
  /* DESCRIPTION:
     This subroutine check if the actual volume is near the volume
     OprogStatus.avVol, if it is the case the program exit 
     This is useful if we want to perform a micronical production run 
     after the equilibration one (NPT) */
  COORD_TYPE rels, Ds;
  COORD_TYPE tols;// relT, DT, temp; 
  int ss;

  tols = OprogStatus.tols;
  if (OprogStatus.avs > 0.0)
    {
      //temp = 2.0 * K / (5.0 * Oparams.parnum - 3.0);
      //DT = fabs(temp - Oparams.T);
      rels = 0.0;
      for (ss = 0; ss < Oparams.PTM; ss++)
	{
	  Ds   = fabs(s[ss] - OprogStatus.avs);
	  //     relT = DT / Oparams.T;
	  rels   += Ds / s[ss];
	}
      /* QUESTION: We must check also T or s1 or it is enough to check
         s? (actually it check only s) */
      /* Controlla se il valore si rels medio e' inferiore ad una data
	 soglia */
      if (rels / ((double) Oparams.PTM) < tols)// && (relT < tols) )
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

extern C_T WLJ[MAX_M];

/* =========================== >>> scaleVelocities <<< ==================== */
void scaleVelocities(int ss, double fact)
{
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

}

/* ============================ >>> move<<< =================================*/
void move(void)
{
  int ss;
  /* DESCRIPTION:
     Move the particles by one step; this procedure is PARALLELIZED !!!
     Valid shared counters are: scN where N = 0..9 
     WARNING: "outer loops" should be ever parallelized, that is you must use 
     loopShr() macro for them inside move() and its subroutine 
     NOTE:
     linked list building could not be parallelized because father and child 
     could disturb each other */


  /* calc predicted coords*/
  if (OprogStatus.Nose != 0)
    moveaNTV(Oparams.steplength, Oparams.m, Oparams.parnum); 
  else
    movea(Oparams.steplength, Oparams.m, Oparams.parnum); 

  for(ss = 0; ss < Oparams.PTM; ss++)
    {
      if (nebrNow[ss])
	{
	  //printf("building neighbour listr step: %d\n", Oparams.curStep);
	  nebrNow[ss] = 0;
	  dispHi[ss] = 0.0;
	
	  /* build up linked list on predicted 
	     coordinates (only father do it)*/
	  if (OprogStatus.noLinkedList)
	    {
	      BuildNebrListNoLinked(ss, Oparams.parnum, Oparams.rcut, 
				    Oparams.sigma);
	    }
	  else
	    {
	      links(ss, Oparams.parnum, Oparams.rcut, Oparams.sigma);
	  
	      /* Build up neighbour list */  
	      BuildNebrList(ss, Oparams.parnum, Oparams.rcut, Oparams.sigma);
	    }
	}
      //printf("lambda0[%d]: %f lambdat: %d\n", ss, 
      //     Oparams.lambda0[ss], Oparams.lambdat[ss]);
      zeroArrays(Fx[ss], Fy[ss], Fz[ss], Oparams.parnum);
      /* Zero all components of pressure tensor */
      Wxy[ss] = 0.0;
      Wyz[ss] = 0.0;
      Wzx[ss] = 0.0;
      Wxx[ss] = Wyy[ss] = Wzz[ss] = 0.0;
    
      if ((OprogStatus.Nose == 1) && (OprogStatus.AntiCry))
	buildMesh(OprogStatus.kmax, OprogStatus.Dk);

      if (OprogStatus.AntiCry)
	{
	  ACForce(ss, Oparams.parnum, OprogStatus.kmax, OprogStatus.Dk, 
		  OprogStatus.alpha, OprogStatus.S0);
	  /* calculate forces */
	  LJForce(ss, Oparams.parnum, Oparams.epsilon, Oparams.sigma,
		  Oparams.rcut);
	  sumFRContribs(ss);
	  EE[ss] = Vc[ss] / Oparams.lambda0[Oparams.lambdat[ss]];
	}
      else
	{
	  /* calculate forces */
	  //printf("STEP: %d\n", Oparams.curStep);
	  LJForce(ss, Oparams.parnum, Oparams.epsilon, Oparams.sigma,
		  Oparams.rcut);
	  W[ss] = WLJ[ss];
	  EE[ss] = Vc[ss] / Oparams.lambda0[Oparams.lambdat[ss]];
	}
    }


  /* correct the coords */
  if (OprogStatus.Nose == 2)
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
  Ktot = kinetTot(Oparams.parnum);

  for (ss = 0; ss < Oparams.PTM; ss++)    
    {
      checkNebrRebuild(ss);
      
      if ( (OprogStatus.CMreset != 0) &&
	   // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	   (Oparams.curStep == OprogStatus.CMreset) ) 
	resetCM(ss, Oparams.parnum);
      
      /* Update the integral of the pressure tensor */
      updateDQ(ss, Oparams.steplength);
      
    }
  
  if (OprogStatus.Nose == 2)
    {
      //scalCor(Oparams.parnum);
      if ( (OprogStatus.sResetSteps != 0) &&
	   (Oparams.curStep == OprogStatus.sResetSteps) )
	{
	  /* NOTE:
	     Reset s to 1 and set to zero its time derivatives
	     This is a trick to reduce the energy drift during equilibration,
	     infact if s about unity the total energy drift less */
	  for (ss = 0; ss < Oparams.PTM; ss++)
	    {
	      s[ss] = 1.0;
	      s1[ss] = 0.0;
	      s2[ss] = 0.0;
	    }
	}
      chks();
    }
  
  chkT();
  
  /* ADDED 15/09/2000: Parallel Tempering
     controlla se e' il caso di fare lo scambio */
  updateV(Oparams.parnum, Oparams.steplength);
  if (OprogStatus.RW == 2)
    updatePE(Oparams.parnum);
  if ((Oparams.curStep % OprogStatus.DtPTEq == 0) && 
      (OprogStatus.RW == 2))
    {
      /* Controlla se il sistema ha equilibrato e genera le curve
	 PEij */ 
      chkPTEq();
    }
  PTexchange();
  saveSnapsPT();
}

