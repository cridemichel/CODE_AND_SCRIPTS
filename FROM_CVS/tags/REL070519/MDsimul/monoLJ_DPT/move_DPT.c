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
extern void links(int Nm, COORD_TYPE rcut, COORD_TYPE sigma);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, 
				  COORD_TYPE sigma);
extern void BuildNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void checkNebrRebuild();
extern void LJForce(int Nm, COORD_TYPE epsilon, 
		    COORD_TYPE sigma, COORD_TYPE rcut, double Lambda);
extern void resetCM(int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
		     COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
		     COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void kinet(int Nm);


extern void zeroArrays(C_T *arrx, C_T* arry, C_T *arrz, int N);
extern void sumFRContribs();
extern void buildMesh(C_T kmax, C_T Dk);
extern void ACForce(int Nm, COORD_TYPE kmax, COORD_TYPE Dk, C_T alpha, C_T S0, 
		    double Lambda);
extern int my_rank;
extern int numOfProcs;

/* ============ >>> move PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, s1t, Vol1t, L, invL, Ktot;   
COORD_TYPE EE[MAX_M],Vc, V, W, K, T1xx, T1yy, T1zz, T1xx, T1yy, T1zz, T1xy, 
  T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy,
  Pzz, Pxy, Pyz, Pzx,
  Patxy, Patyz, Patzx, Patxx, Patyy,
  Patzz, T1myz, T1mzx, T1mxx, T1myy,
  T1mzz;  
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
 
  
  //printf("movea: %d\n", Oparams.PTM);
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
}

/* ========================== >>> movea <<< =============================== */
void moveaNTV(COORD_TYPE dt, COORD_TYPE m, int Nm)
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
  s = s + dt * s1 + dtSq2 * s2;
  s1t = s1;       /* s1(t) */
  s1 = s1 + dt2 * s2;
}

/* ======================= >>> calcT1diagMol <<< =========================== */
COORD_TYPE  calcT1diagAt(int Nm, COORD_TYPE VOL1)
{
  /* calculate (T1mxx+T1myy+T1mzz) / 3.0, that is the average of diagonal
     terms of the kinetic part of atomic pressure tensor */
  int i;
  COORD_TYPE kin, m;
  
  m = Oparams.m;
  kin = 0.0;
  loop(i, 1, Nm)
    {
	  kin +=  m*(Sqr(vx[i]) + 
	    Sqr(vy[i]) + 
	    Sqr(vz[i]));
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
void calcPtensAt(int Nm)
{
  /* Calculate all components of atomic pressure tensor */
  int i;
  COORD_TYPE px, py, pz, m;
  
  m = Oparams.m;
  T1xy = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yz = 0.0;
  T1zx = 0.0;
  T1xx = 0.0; /* Kinetic off-diagonal terms of pressure tensor */
  T1yy = 0.0;
  T1zz = 0.0;

  loop(i, 1, Nm)
    {
      px = vx[i];
      py = vy[i];
      pz = vz[i];
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
  const int NUMIT = 6;

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

 
  loop(k, 1, NUMIT) /* Loop to convergence (NUMIT ~ 5 is enough)*/
    {
      kinet(Nm);
      /* Ora il numero di gradi di liberta' e' numOfProcs(=M nella'art.
	 di Kob) superiore */
      DT = s * (2.0 * K - (3.0 * Nm - 3.0) 
		* Oparams.T) / OprogStatus.Q;
      //DT = s * (2.0 * Ktot - Oparams.PTM * (3.0 * Nm - 3.0) 
      //	* Oparams.T) / OprogStatus.Q;
      
      //printf("[RANK:%d|iter:%d]DT:%.20f\n",my_rank, k, DT);
      A = s1i + 0.5 * dt * DT;
      /* s1(t+dt) */
      s1 = A + 0.5 * Sqr(A) * dt / s + Sqr(dt) * 0.5 * 
	    Sqr(A) * A / Sqr(s);
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

  kinet(Nm);/* K(t+dt) */

  /* Calculate all components of atomic pressure tensor */
  calcPtensAt(Nm);
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
  
  
  loop(i, 1, Nm)
    {
      vxo2[i] = vxo1[i]; /* vo2(t+dt) = v(t-dt) */
      vyo2[i] = vyo1[i];
      vzo2[i] = vzo1[i];
      vxo1[i] = vxt[i];  /* vo1(t+dt) = v(t) */
      vyo1[i] = vyt[i];
      vzo1[i] = vzt[i];
	}
  
  s1o2  = s1o1;
  s1o1 = s1t;
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
  L = cbrt(Vol);
  invL = 1.0 / L;

  
  K = 0.0;
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
  PxyArr = OprogStatus.PxyArr;
  PyzArr = OprogStatus.PyzArr;
  PzxArr = OprogStatus.PzxArr;

  lt = Oparams.lambdat;

  //printf(" curStep: %d\n", curStep);
  /* Store points to use later for integration */
  if (curStep == 1)
    {
      /* WARNING: we assume here that the simulation starts from 1, 
	 that is the first step must be 1 <<<========= !!!!!!!!!! */
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

/* =========================== >>> updateV <<< ========================== */
void updateV(int Nm, C_T dt)
{
  int i, *lt;
  /*
             _ t
	    |
    V_i  =  |  v_i dt
            |
       	   - 0
	
  */
  lt = Oparams.lambdat;
  for(i = 0; i < Nm; i++)
    {
      /* lt[my_rank] e' la temperatura a cui si trova la replica
	 my_rank */
      OprogStatus.sumVx[i] += vx[i]*dt;
      OprogStatus.sumVy[i] += vy[i]*dt;
      OprogStatus.sumVz[i] += vz[i]*dt;
    }
}


/* ============================ >>> chkPTEq <<< =============================*/
void chkPTEqB(void)
{
  double PEij[MAX_M][PE_POINTS];
  int dims[MAX_M], displs[MAX_M];
  
  double lambdai, lambdaj, norm, maxE;
  double ENmin, ENmax, beta0, dblE, refPE, actPE, diffTot, maxPE;
  double bot[5];
  const int refj = 0;
  int maxiE, iE, ss, Nm, nDiff, i;
    
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  Nm = Oparams.parnum;
  beta0 = 1.0/Oparams.T;

  /* Raccoglie tutte le distribuzioni di energia da tutti i processi */
  for (i = 0; i <  Oparams.PTM; i++)
    {
      dims[i] = PE_POINTS;
      /* Così l'array proveniente dal processo di rank 'my_rank' verra' messo
	 nella posizione lambdat[my_rank] e in tal modo le varie distribuzioni
	 di energia saranno ordinate in base alla temperatura */ 
      displs[i] = PE_POINTS * Oparams.lambdat[i];
   }
 
  /* Raccoglie tutte le distribuzioni ordinate in base alla temperatura
     cioe' l'i-esima distribuzione e' quella che relativa a lambda(t) = i */
  MPI_Allgatherv(OprogStatus.PE, PE_POINTS, MPI_INTEGER, PEint, dims, displs,
	      MPI_INTEGER, MPI_COMM_WORLD);

  
  /* Tutti i processi hanno tutte le distribuzioni e quindi fanno gli stessi
     conti giungendo alle stesso conclusioni */

  /* Questa e' il livello al di sotto del quale i valori delle
     distribuzioni di energia vengono tagliati */
  /* Se non c'e' abbastanza statistica esci */
  maxPE = 0.0;
  maxiE = 0;
  /* Determina il massimo della curva di riferimento e usa il valore 
     dell'energia dell massimo come riferimento per l'energia */
  for (iE = 0; iE < PE_POINTS; iE++)
    {
      if (PEint[refj][iE] > maxPE)
	{
	  maxiE = iE;
	  maxPE = PEint[refj][iE];
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
	  PEij[ss][iE] = ((double) PEint[ss][iE])*
	    exp(beta0*(lambdai-lambdaj)*dblE);
	}
      /* Calcola il fattore di normalizzazione integrando con il metodo
	 di Bode (ved. Numerical Recipe) */
      norm = 0.0;
      for (iE = 0; iE < PE_POINTS - 4; iE = iE + 4)
	{
	  for (i = 0; i < 5; i++)
	    {
	      /* Notare che l'energia e' opportunamente shiftata per 
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
void chkPTEq(void)
{
  /* Questa routine differisce dalla precedente poiche' ogni curva e'
     confrontata con la curva riscalata successiva ossia vengono calcolate
     PE_i(E, beta0*lambda_i) e PE_i(beta0*lambda_i+1) dove i = 1...M 
     (ved. cond-mat/0001042) */
  double PEij[2*MAX_M][PE_POINTS];
  int dims[MAX_M], displs[MAX_M];

  double lambdai, lambdaj, norm, maxE;
  double ENmin, ENmax, beta0, dblE, refPE, actPE, diffTot, maxPE;
  double bot[5];
  int jj, maxiE, iE, ss, Nm, nDiff, nstp, i;
  
  /* Raccoglie tutte le distribuzioni di energia da tutti i processi */
  for (i = 0; i < Oparams.PTM; i++)
    {
      dims[i] = PE_POINTS;
      /* Così l'array proveniente dal processo di rank 'my_rank' verra' messo
	 nella posizione lambdat[my_rank] e in tal modo le varie distribuzioni
	 di energia saranno ordinate in base alla temperatura */ 
      displs[i] = PE_POINTS * Oparams.lambdat[i];
   }
  
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
	 if (PEint[ss][iE] > maxPE)
	   {
	     maxiE = iE;
	     maxPE = PEint[ss][iE];
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
	 PEij[2*ss][iE] = ((double) PEint[ss][iE])*
	   exp(beta0*(lambdai-lambdaj)*dblE);
	 PEij[2*ss+1][iE] = ((double) PEint[ss+1][iE]);
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

/* ============================= >>> collectExchg <<< ====================== */
void collectExchg(int *ptrScam, int *ptrAtt)
{
  /* ptrAtt e ptrScam sono puntatori all'area di memoria in cui
     memorizzare la somma degli scambi e dei tentativi di scambio di tutti
     i processi, tali puntatori possono coincidere con OprogStatus.attmpted
     e OprogStatus.scambi solo se subito dopo tali arrays vengono azzerati
     infatti diversamente accumuleremmo valori sbagliati
  */
  int i;
  int att[MAX_M], scam[MAX_M];


  for (i = 0; i < Oparams.PTM; i++)
    {
      att[i] = OprogStatus.attempted[i];
      scam[i] = OprogStatus.scambi[i];
      //sprintf(TXT, "[%d] att: %d sca:%d\n", i, att[i], scam[i]);
      //mdPrintf(ALL, TXT, NULL);
    }
  // printf("porca troia\n");

  MPI_Allreduce(att, ptrAtt, Oparams.PTM, MPI_INTEGER, 
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(scam, ptrScam, Oparams.PTM, MPI_INTEGER, 
		MPI_SUM, MPI_COMM_WORLD);
  /*
    for (i = 0; i < Oparams.PTM; i++)
    {
    sprintf(TXT, "[%d] OS.att: %d OS.sca:%d\n", i, 
    OprogStatus.attempted[i], OprogStatus.scambi[i]);
    mdPrintf(ALL, TXT, NULL);
    }
  */

}

/* ============================ >>> collectlt <<< ============================== */
void collectlt(void)
{
  int lttmp;
  //double EEtmp;
  /* first of all get all lambda(t) */
  /* questo e' inefficiente ma tanto i dati sono una manciata */
  lttmp = Oparams.lambdat[my_rank];
  //EEtmp = EE[my_rank];
  MPI_Allgather(&lttmp, 1, MPI_INTEGER, 
		Oparams.lambdat, 1, MPI_INTEGER, MPI_COMM_WORLD);
  //MPI_Allgather(&EEtmp, 1, MPI_DOUBLE, 
  //		EE, 1, MPI_DOUBLE, MPI_COMM_WORLD);


}

void collectEE(void)
{
  double EEtmp;
  EEtmp = EE[my_rank];
  MPI_Allgather(&EEtmp, 1, MPI_DOUBLE, 
		EE, 1, MPI_DOUBLE, MPI_COMM_WORLD);

}
/* ============================= >>> zeroPE <<< ============================ */
void zeroPE(void)
{
  /* Azzera tutte le PE poiche' ormai il sistema dovrebbe aver equilibrato,
     avendo completato mezzo random walk */
  int i;
  for(i = 0; i < PE_POINTS; i++)
    {
      OprogStatus.PE[i] = 0;
    }
}
/* ============================= >>> chkPT <<< ============================= */
/* ADDED 15/09/2000: Parallel Tempering */
void PTexchange()
{
  double soglia, w;//, Ffact0, Ffact1;
  double neighE, myE, Delta, rate;
  int  mylt, neighlt, i;
  int primoSS;
  int *lt, ss;
  double *l0;
  MPI_Status stato;
  int exchg, tmpAtt[MAX_M], tmpScam[MAX_M];
  /* da il sistema che corrisponde ad un certo Lambda 
   e' la funzione inversa di lambda[t] ossia:
   lambda[rank[ss]] = ss 
  */
  int rank[MAX_M], itmp;

  /* Oparams.lambdat e' un array di puntatori ai lambda iniziali, i
     n sostanza di ogni sottosistema dice qual'e' il lambda che gli 
     appartiene ad un dato istante */   
  lt = Oparams.lambdat;
  l0 = Oparams.lambda0;
	
  if (OprogStatus.PT == 0) 
    return;
  
  if ( (OprogStatus.RW == 0) && 
       /* PATCH: le repliche sono Oparams.PTM ma lambdat varia tra 0 e 
	  Oparams.PTM-1 */ 
       (Oparams.lambdat[OprogStatus.refReplica] == (Oparams.PTM-1) ) )
    {
      OprogStatus.RW = 1;
    }
  
  if ( (OprogStatus.RW == 1) && 
       (Oparams.lambdat[OprogStatus.refReplica] == 0) )
    {
      OprogStatus.RW = 2;
      zeroPE();
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
      collectExchg(tmpScam, tmpAtt);
      /* BUG COLOSSALE: senza questo utilizzava l'energia di 1000 passi prima 
	 ma deve usare per valutare la probabilita' di scambio l'energia attuale!
       */
      collectEE();
      if (tmpAtt[0] == 0)
	{
	  sprintf(TXT, "[PTexchange(%d)] T = [ %.5f(att=0) ", 
		  Oparams.curStep,
		  Oparams.T / Oparams.lambda0[0]);
	}
      else
	{
	  sprintf(TXT, "[PTexchange(%d)] T = [ %.5f(%d/%d=%.5f) ", 
		  Oparams.curStep,
		  Oparams.T / Oparams.lambda0[0], 
		  tmpScam[0],
		  tmpAtt[0],
		  ((double)tmpScam[0])/
		  ((double)tmpAtt[0]));
	}
      for (i = 1; i < Oparams.PTM - 1; i++)
	{
	  if (tmpAtt[i] == 0)
	    {
	      sprintf(TXT1, ", %.5f(att=0)", 
		      Oparams.T / Oparams.lambda0[i]);
	      strcat(TXT, TXT1);
	    }
	  else
	    {
	      rate = 
		((double ) tmpScam[i]) / ((double) tmpAtt[i]);
	      sprintf(TXT1, ", %.5f(%d/%d=%.5f)", 
		      Oparams.T / Oparams.lambda0[i],
		      tmpScam[i], tmpAtt[i],
		      rate);
	      strcat(TXT, TXT1);
	    }
	}
      sprintf(TXT1, ", %.5f(0)", Oparams.T / Oparams.lambda0[Oparams.PTM-1]);
      strcat(TXT, TXT1);
      /* log out all new temperatures */
      strcat(TXT, " ]\n");
      mdPrintf(ALL, TXT, NULL);

      /* calcola la mappa inversa rispetto a lambda(t) */
      for (ss = 0; ss < Oparams.PTM; ss++)
	{
	  rank[lt[ss]] = ss;
	  //sprintf(TXT, "[%d]=%d\n", ss,rank[lt[ss]]);
	  //mdPrintf(ALL, TXT, NULL);
	}
      
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
	  /* 03/10/2000 ADDED s reset in chkPT */
	  /* Resetta le variabili s */
	  /* forse cio' non serve !!!!!!!!!!!!!*/
	  s = 1.0;
	  s1 = 0.0;
	  s2 = 0.0;
	  s1o1 = 0.0;
	  s1o2 = 0.0;
	}
      strcpy(TXT1, "Prob. Sc.:");
      /* Notare che devono inviare decidere di scambiare tutti i processi 
	 che hanno una temperatura pari o dispari non my_rank, infatti
	 se dato my_rank la temperatura che ha quel processo potrebbe essere
	 qualsiasi e cio' che conta e' che solo processi con temperature attigue
	 possono scambiarsi la temperatura.
	 Notare anche che se my_rank e' l'ultimo processo allora non deve scambiare 
	 perche' non c'e' un processo successivo. */
      if ( (lt[my_rank] - primoSS) % 2 == 0 ) 
	/* se primoSS = 0 => lt[my_rank] = 0 o 2 o 4 ...
	   se primoSS = 1 => lt[my_rank] = 1 o 3 o 5 ... */
	
	{
	  if (lt[my_rank] + 1 == Oparams.PTM)
	    goto fine;

	  //mdPrintf(ALL, "dentro", NULL);
	  soglia = ranf();

	  //printf("soglia: %f\n", soglia);
	  //printf("soglia[%d]: %f\n", ss, soglia);
	  /* NOTA BENE: notare che EE[i] e' l'energia del processo di rank i*/
	  mylt = lt[my_rank];
	  myE = EE[my_rank];
	  OprogStatus.attempted[mylt]++;	
	  neighlt = mylt+1;
	  neighE = EE[rank[mylt+1]];
	  printf("neighlt:%d mylt:%d rank[neighlt]: %d\n", 
		 neighlt, mylt, rank[neighlt]);
	  printf("my_rank: %d lt[my_rank]:%d rank[mylt]:%d\n", my_rank, lt[my_rank], 
		 rank[mylt]);
	  /* Mt[my_rank+1] da l'energia del sistema che ha la temperatura
             attigua a my_rank */
	
	  Delta = 1.0/Oparams.T * 
	    (l0[neighlt] - l0[mylt]) * (myE - neighE);
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
	      OprogStatus.scambi[mylt]++;
	      /* Scambia anzitutto i lambda (lo scambio effettivo lo fara'
		 poi con una MPI_Allgather) */
	      lt[my_rank] = mylt + 1;                                   
	      
	      //lt[rank[neighlt]]= mylt;
	      	      
	      /* e poi tutte le misure... */
	      exchg = my_rank;
	 
	      MPI_Ssend(&exchg, 1, MPI_INTEGER, rank[neighlt], 
			MD_MPI_EXCHG, MPI_COMM_WORLD);
	      
	      
	      MPI_Sendrecv_replace(OprogStatus.sumVx, Oparams.parnum, 
				   MPI_DOUBLE, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);
	      MPI_Sendrecv_replace(OprogStatus.sumVy, Oparams.parnum,
				   MPI_DOUBLE, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);
	      MPI_Sendrecv_replace(OprogStatus.sumVz, Oparams.parnum, 
				   MPI_DOUBLE, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);

	      MPI_Sendrecv_replace(OprogStatus.PE, 
				   PE_POINTS, MPI_INTEGER, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);

	    }
	  else
	    {
	      exchg = -1;
	      MPI_Ssend(&exchg, 1, MPI_INTEGER, rank[neighlt], 
			MD_MPI_EXCHG, MPI_COMM_WORLD);
	    }
	  
	  
	  sprintf(TXT, "[%.4f]",
		  ((double) OprogStatus.scambi[mylt]) / 
		  ((double)OprogStatus.attempted[mylt]) );
	  strcat(TXT1, TXT);
	  mdPrintf(ALL, TXT, NULL);
	}
      else
	{
	  /* se my_rank = 0 il processo non deve effettuare nessuno scambio, non 
	     esistendo un processo di rank minore */
	  
	  mylt = lt[my_rank];
	  
	  neighlt = mylt - 1; 
	  //sprintf(TXT,"my_rank: %d neighlt:%d mylt:%d rank[neighlt]: %d\n", 
	  //	  my_rank, neighlt, mylt, rank[neighlt]);
	  //mdPrintf(ALL, TXT, NULL);
	  
	  if (neighlt >= 0)
	    MPI_Recv(&exchg, 1, MPI_INTEGER, rank[neighlt], MD_MPI_EXCHG, 
		     MPI_COMM_WORLD, &stato);
	  else
	    exchg = -1;
	  
	  //mdPrintf(ALL, "ricevo", NULL);
	  /* exchg == -1 => lo scambio non va fatto 
	     exchg != -1 => exchg = {rank del processo con cui scambiare} */
	  if (exchg != -1)
	    {
	     
	      if (exchg != rank[neighlt])
		{
		  mdPrintf(STD, "['exchg != rank[neighlt]']: Che cazzo sta succedendo?!?!!!???\n", NULL);
		  MPI_Abort(MPI_COMM_WORLD, -1);
		}
	      //lt[exchg]= mylt;
	      
	      /* scambia i lambda e l'energia potenziale */
	      lt[my_rank] = neighlt;    

	      MPI_Sendrecv_replace(OprogStatus.sumVx, Oparams.parnum, 
				   MPI_DOUBLE, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);
	      MPI_Sendrecv_replace(OprogStatus.sumVy, Oparams.parnum, 
				   MPI_DOUBLE, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);
	      MPI_Sendrecv_replace(OprogStatus.sumVz, Oparams.parnum, 
				   MPI_DOUBLE, rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);
	      MPI_Sendrecv_replace(OprogStatus.PE, PE_POINTS, MPI_INTEGER, 
				   rank[neighlt], 
				   MD_MPI_EXCHG, rank[neighlt], MD_MPI_EXCHG, 
				   MPI_COMM_WORLD, &stato);
	      
	    }
	  
	}
    fine:
      collectlt();

      /*      for (i = 0; i < Oparams.PTM; i++)
	{
	sprintf(TXT, "DOPO lambdat[%d]:%d\n",i, Oparams.lambdat[i]);
	mdPrintf(ALL, TXT, NULL);
	}
    
	strcat(TXT1,"finito\n");*/
    }
}
/* ============================= >>> updatePE <<< ========================= */
void updatePE(int Nm)
{
  int iE, *lt;
  double ENmin, ENmax;
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  lt = Oparams.lambdat;
  iE = (int) ((EE[my_rank] - ((double) Nm)*ENmin) / 
	      ( ((double) Nm) * ((double) ENmax - ENmin)) * 
	      ((double) PE_POINTS));
  //  sprintf(TXT, "energia: %f iE:%d\n", EE[my_rank], iE);
  //mdPrintf(ALL, TXT, NULL);
  if ( (iE >= 0) && (iE < PE_POINTS)) 
    {
      ++(OprogStatus.PE[iE]);
    }

}

/* ============================= >>> chkT <<< ============================== */
void chkT(void)
{
  C_T DT, temp, relT;
 
  if (OprogStatus.tolT > 0.0)
    {
      kinet(Oparams.parnum);
      /* Calcola la temperatura "totale" */
      temp = 2.0 * K / Oparams.PTM * (5.0 * Oparams.parnum - 3.0);
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

  tols = OprogStatus.tols;
  if (OprogStatus.avs > 0.0)
    {
      //temp = 2.0 * K / (5.0 * Oparams.parnum - 3.0);
      //DT = fabs(temp - Oparams.T);
      rels = 0.0;
      Ds   = fabs(s - OprogStatus.avs);
	  //     relT = DT / Oparams.T;
      rels   = Ds / s;
    
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

extern C_T WLJ;

/* ========================= >>> funz <<< ================================ */
double funz(double c, double r)
{
  /* MODIFICARE !!!!!!!!!!!!!!!!!! */
  //return c*x*r;

  /* FUNZIONE ALTERNATIVA 
     -0.5 < r < 0.5 , x invece è la distanza fra due lambda consecutivi */
  double ret;

  if (r >= 0)
    {
      ret = 1.0 + OprogStatus.optFact * 2.0 * r; 
      /* aggiunge x cioe' raddoppia la distanza */
    }
  else
    {
      ret = 1.0 + OprogStatus.optFact * r;  /* toglie x/2 cioe' dimezza la distanza */ 
    }
return ret;
  
}

/* ========================== >>> optLamb <<< ============================ */
void optLamb()
{
  double DLtot, DL[MAX_M], Drate, rate, dtmp;
  int neighlt, mylt, stopOpt = 0, i, esciOpt = 0;
  int dims[MAX_M], displs[MAX_M];
  char Tstr[1024], tmpStr[255];
  int rank[MAX_M], *lt;
  double TARGET_RATE = 0.5;
  double RATE_TOL = 0.1;

  /* RATE_TOL è la tolleranza in percentuale entro la quale si puo' ritenere 
     di aver raggiunto un valore accettabile di rate di scambio */

  RATE_TOL = OprogStatus.RATE_TOL;
  TARGET_RATE = OprogStatus.TARGET_RATE;

  if (OprogStatus.optL == 0)
    {
      mdPrintf(ALL, "L'ottimizzazione è finita!\n", NULL);
      return;
    }
  /* Raccoglie tutte i lambda0 di energia da tutti i processi */
  
  collectExchg(OprogStatus.scambi, OprogStatus.attempted);
  mylt = Oparams.lambdat[my_rank];
  lt = Oparams.lambdat;

  //  if (mylt + 1 < Oparams.PTM )
  neighlt = mylt + 1;
  //  else
  //  neighlt = 0;

  /* Costruisce i Delta Lambda */
  if (lt[my_rank]+1 < Oparams.PTM) 
    {
      DL[lt[my_rank]] = Oparams.lambda0[lt[my_rank]+1] - 
	Oparams.lambda0[lt[my_rank]];
    }
  else
    {
      /* L'ultima temperatura non ha un temperatura vicina più grande! */
      DL[Oparams.PTM - 1] = 0;
    }
  
  for (i = 0; i < Oparams.PTM; i++)
    {
      rank[lt[i]] = i;
      //sprintf(TXT, "[%d]=%d\n", ss,rank[lt[ss]]);
      //mdPrintf(ALL, TXT, NULL);
    }

  for (i = 0; i < Oparams.PTM - 1; i++)
    {
      esciOpt |= (OprogStatus.attempted[i] == 0); 
    }

  if (esciOpt)
    {
      mdPrintf(ALL, "Non cambio i lambda, esco...\n", NULL);
      return;
    }
  /* controlla se tutti i rate sono entro la tolleranza, se è così => fine 
     ottimizzazione */
  if (neighlt < Oparams.PTM)
    {
      rate = 
	((double ) OprogStatus.scambi[0]) / ((double) OprogStatus.attempted[0]);
      stopOpt = (  ( fabs(rate - TARGET_RATE) / TARGET_RATE ) <= RATE_TOL  ); 
      sprintf(Tstr, "[optLamb(%d)] T = [ %.5f(%d/%d=%.4f) ", Oparams.curStep,
	      Oparams.T / Oparams.lambda0[0], OprogStatus.scambi[0],
	      OprogStatus.attempted[0], rate);
    }

  for (i = 1; i < Oparams.PTM - 1; i++)
    {
      rate = 
	((double ) OprogStatus.scambi[i]) / ((double) OprogStatus.attempted[i]);
      sprintf(tmpStr, ", %.5f(%d/%d=%.4f)", 
	      Oparams.T / Oparams.lambda0[i], 
	      OprogStatus.scambi[i], OprogStatus.attempted[i], rate);
      strcat(Tstr, tmpStr);
      
      stopOpt &= ( (fabs(rate - TARGET_RATE) / TARGET_RATE) <= RATE_TOL ); 
    }
  sprintf(tmpStr, ", %.5f(0)", Oparams.T / Oparams.lambda0[Oparams.PTM-1]);
  strcat(Tstr, tmpStr);
  /* log out all new temperatures */
  strcat(Tstr, " ]\n");
  mdPrintf(ALL, Tstr, NULL);
 
  if (stopOpt)
    {
      mdPrintf(ALL, "Optimization completed successfully!\n", NULL);
      OprogStatus.optL = 0; 
      /* disattiva l'ottimizzazione delle distanze in temperatura poiché
	 il target è stato raggiunto */
    }
  else
    {
      /* calcolo del rate del sistema attuale (my_rank) */
      rate = ((double ) OprogStatus.scambi[mylt]) / ((double) OprogStatus.attempted[mylt]);
      sprintf(TXT, "mylt:%d scambi:%d attemp:%d rate: %.5f\n", mylt,
	      OprogStatus.scambi[mylt],
	      OprogStatus.attempted[mylt], rate);
      mdPrintf(ALL, TXT, NULL);
      if (neighlt < Oparams.PTM)
	{
	  Drate = rate - TARGET_RATE;
	  /* Al massimo se Drate = 0.5 => Dlamb = |lambda(mylt) - lambda(neighlt)| */  
	  sprintf(TXT, "Drate:%.5f DL[%d]:%.5f\n", Drate, lt[my_rank],
	      funz(OprogStatus.optFact, Drate));
	  mdPrintf(ALL, TXT, NULL);
	  DL[lt[my_rank]] *= funz(OprogStatus.optFact, Drate);
	}
      
      for (i = 0; i < Oparams.PTM; i++)
	{
	  dims[i] = 1;
	  /* Così l'array proveniente dal processo di rank 'my_rank' verra' messo
	     nella posizione lambdat[my_rank] e in tal modo le varie distribuzioni
	     di energia saranno ordinate in base alla temperatura */ 
	  displs[i] = lt[i];
	  sprintf(TXT, "lambdat[%d]:%d\n", i, lt[i]);
	  mdPrintf(ALL, TXT, NULL);
	}
      
      dtmp = DL[mylt];
      MPI_Allgatherv(&dtmp, 1, MPI_DOUBLE, 
		     DL, dims, displs, MPI_DOUBLE, MPI_COMM_WORLD);
      
      
      /* Ricostruisce i lambda (PER ORA ASSUMO che lambda0[0] = 1!!!!!!!!!!!!) */
      Oparams.lambda0[0] = 1.0;
      DLtot = 0.0;
      for (i = 1; i < Oparams.PTM; i++)
	{
	  DLtot += DL[i-1];  
      sprintf(TXT, "lt[my_rank]:%d DL[%d]:%.5f\n", lt[my_rank], i-1, DL[i-1]);
      mdPrintf(ALL, TXT, NULL);
      Oparams.lambda0[i] = Oparams.lambda0[0] + DLtot; 
	}
      
      /* reset scambi and attempted to 0 */
    }
  for (i = 0; i < Oparams.PTM; i++) 
    {
      OprogStatus.scambi[i] = 0;
      OprogStatus.attempted[i] = 0;
    }
  zeroPE();
  sprintf(TXT, "=====>>>>>>[%d] scambi:%d attemp:%d\n",mylt,
	  OprogStatus.scambi[mylt], OprogStatus.attempted[mylt]);
  mdPrintf(ALL, TXT, NULL);
}

/* =========================== >>> scaleVelocities <<< ==================== */
void scaleVelocities(double fact)
{
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

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
  double Lambda;
  /*
  int i;
  for (i = 0; i < Oparams.PTM; i++)
    {
      sprintf(TXT, "lambda[%d]: %.5f\n", Oparams.lambdat[i], Oparams.lambda0[Oparams.lambdat[i]]);
      mdPrintf(ALL,TXT, NULL);
      }*/
  /* calc predicted coords*/
  if (OprogStatus.Nose != 0)
    moveaNTV(Oparams.steplength, Oparams.m, Oparams.parnum); 
  else
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
	  BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut, 
				Oparams.sigma);
	}
      else
	{
	  links(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	  
	  /* Build up neighbour list */  
	  BuildNebrList(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	}
    }

  zeroArrays(Fx, Fy, Fz, Oparams.parnum);
  /* Zero all components of pressure tensor */
  Wxy = 0.0;
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz  = 0.0;
  
  if ((OprogStatus.Nose == 1) && (OprogStatus.AntiCry))
    buildMesh(OprogStatus.kmax, OprogStatus.Dk);
  
  if (OprogStatus.AntiCry)
    {
      Lambda = Oparams.lambda0[Oparams.lambdat[my_rank]];
      ACForce(Oparams.parnum, OprogStatus.kmax, OprogStatus.Dk, 
	      OprogStatus.alpha, OprogStatus.S0, Lambda);
      /* calculate forces */
      LJForce(Oparams.parnum, Oparams.epsilon, Oparams.sigma,
	      Oparams.rcut, Lambda);
      sumFRContribs();
      EE[my_rank] = Vc / Oparams.lambda0[Oparams.lambdat[my_rank]];
    }
  else
    {
      /* calculate forces */
      Lambda = Oparams.lambda0[Oparams.lambdat[my_rank]];
      LJForce(Oparams.parnum, Oparams.epsilon, Oparams.sigma,
	      Oparams.rcut, Lambda );
      W = WLJ;
      EE[my_rank] = Vc / Oparams.lambda0[Oparams.lambdat[my_rank]];
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
  kinet(Oparams.parnum);

  checkNebrRebuild();
      
  if (  ( (OprogStatus.CMreset > 0) &&
	  // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	  (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	    // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	    (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
    resetCM(Oparams.parnum);
      
  /* Update the integral of the pressure tensor */
  updateDQ(my_rank, Oparams.steplength);
  
  if (OprogStatus.Nose == 2)
    {
      //scalCor(Oparams.parnum);
      /* se sResetSteps > 0  allora resetta una volta sola ma se
	 sResetSteps < 0 allora resetta s ogni sResetSteps */
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
      chks();
    }

  chkT();
  
  /* ADDED 15/09/2000: Parallel Tempering
     controlla se e' il caso di fare lo scambio */
  updateV(Oparams.parnum, Oparams.steplength);
  
  //if (OprogStatus.RW == 2)
  updatePE(Oparams.parnum);
  
  if ((Oparams.curStep % OprogStatus.DtPTEq == 0) && 
      (OprogStatus.RW == 2))
    {
      /* Controlla se il sistema ha equilibrato e genera le curve
	 PEij */ 
      chkPTEq();/* questa routine da problemi su Cray e Origin ovveo genera un 
                   floating point execption */
    }
  
  PTexchange();
  
  if (Oparams.curStep % (2 * OprogStatus.DtPT * OprogStatus.DtOptL) == 0)
    {
      optLamb();
      /* Aumenta progressivamente l'intervallo per avere risulatati
	 sempre più precisi */
      //OprogStatus.DtOptL = rint(OprogStatus.DtOptL * 1.2);
    }

  saveSnapsPT();
}

