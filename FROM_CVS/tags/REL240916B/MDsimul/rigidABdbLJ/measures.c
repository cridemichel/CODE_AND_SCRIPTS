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

/* See after for declaration */
void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Wm, Wmxx, Wmyy, Wmzz, 
  Wmxy, Wmyz, Wmzx, Pmxx, Pmyy, Pmzz, Pmxy, Pmyz, Pmzx, T1mxy, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
  T1myz, T1mzx, T1mxx, T1myy, T1mzz;  

extern COORD_TYPE Mtot;
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

extern COORD_TYPE *FxL[NA], *FyL[NA], *FzL[NA]; /* Local arrays of forces
					    (used by the LJForce routine 
					    to avoid interferences) */
extern COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

extern COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
extern COORD_TYPE  *Rmx, *Rmy, *Rmz;
/* ========================================================================= */
/*=========================== >>> vectProd <<< =========================== */
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);

/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
  COORD_TYPE dist, Px, Py, Pz, cost, Rx, Ry, Rz, RCMx, RCMy, RCMz, Mtot;
  int mol, Nm, i, a;
  COORD_TYPE px, py, pz, lx, ly, lz, Lx, Ly, Lz;
  COORD_TYPE Rx_scal, Ry_scal, Rz_scal;
  COORD_TYPE L, invL;

  Nm = Oparams.parnum;
  mdShow("SHIFTED POTENTIAL ENERGY: %.10f", Vc);
  printf("KINETIC ENERGY: %.10f\n", K);
  printf("s: %.10f s1: %.10f Vol: %.10f Vol1: %.10f\n", s, s1, Vol, Vol1);
  E = K + Vc;
  /* And now add the contribute due to the thermal bath */
  E = E + 0.5 * Sqr(s1) * OprogStatus.Q / Sqr(s) + 
    (5.0 * ((COORD_TYPE) Nm) - 3.0) * Oparams.T * log(s) + 
    0.5 * Sqr(Vol1) * OprogStatus.W / Sqr(s) + Oparams.P * Vol; 
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */

  printf("TOTAL ENERGY: %.10f\n", E);
  printf("Elrc:: %.6f Plrc: %.6f\n", Elrc, Plrc);
  mol = 10;
  dist = sqrt( Sqr(rx[0][mol] - rx[1][mol]) + Sqr(ry[0][mol] - ry[1][mol]) + 
	       Sqr(rz[0][mol] - rz[1][mol]) );
  printf("Atom position: %.15f\n", rx[0][mol]);
  printf("Molecule atoms distance: %.15f\n", dist);
  dist= (rx[0][mol] - rx[1][mol]) * (vx[0][mol] - vx[1][mol]) + 
    (ry[0][mol] - ry[1][mol]) * (vy[0][mol] - vy[1][mol]) + 
    (rz[0][mol] - rz[1][mol]) * (vz[0][mol] - vz[1][mol]);
  printf("Molecule relative velocity along bond: %.15f\n", dist); 

  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  Lx = 0.0;
  Ly = 0.0;
  Lz = 0.0;
  Mtot = Oparams.m[0] + Oparams.m[1];
  cost = Vol1 / Vol / 3.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  //L = cbrt(Vol);
  //invL = 1.0 / L;

  loop(i, 1, Nm)
    {
      CoM(i, &Rx, &Ry, &Rz);
      
      px = 0.0;
      py = 0.0;
      pz = 0.0;
      
      loop(a, 1, NA)
	{
	  px += Oparams.m[a] * (vx[a][i] - 
	    cost * Rx);
	  py += Oparams.m[a] * (vy[a][i] - 
	    cost * Ry);
	  pz += Oparams.m[a] * (vz[a][i] - 
	    cost * Rz);
	  //vectProd(rx[a][i], ry[a][i], rz[a][i], 
	       //     vx[a][i], vy[a][i], vz[a][i],
          //&lx, &ly, &lz);
	  //Lx += lx;
	  //Ly += ly;
	  //Lz += lz;
	}
      Rx_scal = Rx - L*rint(invL * Rx);
      Ry_scal = Ry - L*rint(invL * Ry);
      Rz_scal = Rz - L*rint(invL * Rz);
      vectProd(Rx_scal, Ry_scal, Rz_scal, px, py, pz, &lx, &ly, &lz);

      
      Px += px;
      Py += py;
      Pz += pz;
      Lx += lx;
      Ly += ly;
      Lz += lz;
    }
  printf("Px: %.15f Py: %.15f Pz: %.15f\n", Px, Py, Pz);
  printf("Lx: %.10f Ly: %.10f Lz: %.10f\n", Lx, Ly, Lz);
  
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  
  loop(i, 1, Nm)
    {
      CoM(i, &Rx, &Ry, &Rz);
      RCMx += Rx;
      RCMy += Ry;
      RCMz += Rz;
    }
  printf("Box CoM RCMx: %.15f RCMy: %.15f RCMz: %.15f\n", RCMx, RCMy, RCMz);
}

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  COORD_TYPE Drx, Dry, Drz, Mtot, invMtot, Dr4;
  COORD_TYPE m0, m1;
  COORD_TYPE rxCM, ryCM, rzCM; /* Center of mass position */
  int i;
  
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  Mtot = m0 + m1;
  invMtot = 1.0 / Mtot;
 
  DrSqTot = 0.0;
  Dr4 = 0.0;

  loop(i, 1, Oparams.parnum)
    {
      rxCM = (m0 * rx[0][i] +  m1 * rx[1][i]) * invMtot;
      ryCM = (m0 * ry[0][i] +  m1 * ry[1][i]) * invMtot;
      rzCM = (m0 * rz[0][i] +  m1 * rz[1][i]) * invMtot;
      Drx = rxCM - OprogStatus.rxCMi[i]; 
      Dry = ryCM - OprogStatus.ryCMi[i];
      Drz = rzCM - OprogStatus.rzCMi[i];
      if (OprogStatus.ipart == i)
	{
	  printf("i = %d\n", i);
	  /* Motion of the OprogStatus.ipart particle */
	  sqrtdr2 = sqrt(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
	}
      DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
      Dr4 += Sqr(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
   }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
 
  Dtrans = DrSqTot / ( 6.0 * ((COORD_TYPE) Oparams.steplength) *
		       ((COORD_TYPE) Oparams.curStep) * 
		       ((COORD_TYPE) Oparams.parnum ) );   
  Aa = ((COORD_TYPE) Oparams.parnum ) * 3.0 * 
    Dr4 / Sqr(DrSqTot) / 5.0 - 1.0; /* Non-Gaussian parameter */  
  
  DrSqTot /= ((COORD_TYPE) Oparams.parnum);
  
  printf("Diffusion coefficent: %.15f\n", Dtrans); 
}

/* ============================ >>> temperat <<< =========================== */
void temperat(void)
{
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
  int i;
  COORD_TYPE invMtot,Mtot, Ktransl, vxCM, vyCM, vzCM, m0, m1;
  temp = 2.0 * K / (5.0 * Oparams.parnum - 3.0);

  Ktransl = 0.0;
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  /*
  Mtot = m0 + m1;
  invMtot = 1.0 / Mtot;
  loop(i, 1, Oparams.parnum)
    {
      vxCM = (m0 * vx[0][i] +  m1 * vx[1][i]) * invMtot;
      vyCM = (m0 * vy[0][i] +  m1 * vy[1][i]) * invMtot;
      vzCM = (m0 * vz[0][i] +  m1 * vz[1][i]) * invMtot;
      Ktransl += 0.5 * Mtot * (Sqr(vxCM) + Sqr(vyCM) + Sqr(vzCM));
    }
  temp_transl = 2.0 * Ktransl  / (3.0 * Oparams.parnum - 3.0);
  */
  if (OprogStatus.avngTemp == 1)
    {
      OprogStatus.sumTemp += temp;
      temp = OprogStatus.sumTemp / NUMCALCS;
    }

  /* pressure */
  if (OprogStatus.avngPress == 1)
    {
      OprogStatus.sumPress += press;
      press = OprogStatus.sumPress / NUMCALCS;
    }
  printf("Instantaneous temperature:%.15f\n", temp);
  printf("Wm: %f WC: %.15f W: %.15f\n", Wm, W, WC);
  //printf("Molecular pressure: %f\n", press_m);


#ifdef MOLPTENS
  printf("MOLECULAR Pressure tensor\n");
  printf("Molecular pressure: %.15f\n", press);
#endif
#ifdef ATPTENS
  printf("ATOMIC Pressure tensor\n");
#endif
#if !defined(MOLPTENS) || defined(ATPRESS) || defined(ATPTENS)
  printf("Atomic pressure: %.15f\n", press_at);
#endif

  printf("Vol: %.15f\n", Vol);
}

/* ======================== >>> structFacts <<< =============================*/
void structFacts(void)
{
  int n, bin;
  COORD_TYPE sum, pi4, rhoAv, r;
  COORD_TYPE scalFact, kr, k, DELR;
  COORD_TYPE* sumS;
  //printf("calcolo fattore di struttura!!!\n");
  sumS = OprogStatus.sumS;
  DELR = (REND - RBEG) / MAXBIN;
  pi4 = pi * 4.0;
  rhoAv = (COORD_TYPE) Oparams.parnum / Vol; /* V = 1 here !!!!!!!!!!!!*/ 
  scalFact = (KEND - KBEG) / NUMK;
  
  loop(n, 1, NUMK)
    {
      k = KBEG + n * scalFact;
      sum = 0.0;
      /* Bode's rule */
      bin = 0;
     
      while(bin + 4 < MAXBIN)
	{
	  r = RBEG + bin * DELR;

	  kr = k*r;
	  sum += 14.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 64.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 24.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 64.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 14.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  //bin++;
	}

      S[n] = 1.0 + pi4 * sum * rhoAv * DELR / 45.0;

      if (OprogStatus.avngS == 1)
	{
	  sumS[n] += S[n];
	  S[n] = sumS[n] / NUMCALCS;
	}
    }
}

/* ============================ >>> structFacts <<< ======================== */
void structFactsB(void)
{
  /* DESCRIPTION:
     This mesuring function calculates the static structure factor */
  COORD_TYPE reRho, imRho, *cosL, *sinL, pi2, teta;
  COORD_TYPE m0, m1, invMtot, invNm;
  COORD_TYPE kSq, rxCM, ryCM, rzCM, rCMk, sumRho;
  COORD_TYPE rx0, ry0, rz0, rx1, ry1, rz1, scalFact, kx, ky, kz, siniTeta;
  int iTeta, iPhi, i, maxI, Nm, n;
  
  //return;
  /* useful quantities */

  L = cbrt(Vol);
  invL = 1.0 / L;
  m0 = Oparams.m[0];
  m1 = Oparams.m[1]; 
  scalFact = (KEND - KBEG) / NUMK;
  pi2 = 2.0 * pi;
  
  maxI = rint(sqrt(NUMK2AV)); 
  //printf("maxI : %d\n", maxI);
  /* We take 'maxI' values of Phi and 'maxI' values of Teta, so we have in this
     way NUMK2AV = maxI * maxI (about) values for 'k' over the 
     semi-sphere with same modulus */
  
  /* Allocate arrays for storin the values of Rho(k) for variuos directions */
  cosL  = malloc(sizeof(COORD_TYPE) * maxI);
  sinL  = malloc(sizeof(COORD_TYPE) * maxI);
  
  /* Build up look up table for sine and cosine */
  loop(iTeta, 1, maxI)
    {
      /* NOTICE:
	 Note that the values of teta for which we need the cosines and sines
	 are identical to the value of phi for which we need them, so we 
         calculate the cosines and sine for different teta values only*/
     
      teta = pi * ( iTeta / ((COORD_TYPE)maxI) ); 
      /* teat goes from 0 to a value next to 2*pi */ 
     
      cosL[iTeta] = cos(teta);
      sinL[iTeta] = sin(teta);
     
    }
  
  invMtot = 1.0 / (m0 + m1);
  Nm = Oparams.parnum;
  invNm = 1.0 / Nm;

  loop(n, 1, NUMK)
    {
      sumRho = 0.0;
      loop(iTeta, 1, maxI)
	{
	  loop(iPhi,  1, maxI) /* from 0 to NUMK-1 */
	    {
	      reRho = 0.0;
	      imRho = 0.0;
	      loop(i, 1, Oparams.parnum)
		{
		  /* Center of mass position */
		  rx0 = rx[0][i] - L * rint(invL * rx[0][i]); 
		  /*Reduced coordinates !!!!*/
		  ry0 = ry[0][i] - L * rint(invL * ry[0][i]);
		  rz0 = rz[0][i] - L * rint(invL * rz[0][i]);
		  rx1 = rx[1][i] - L * rint(invL * rx[1][i]);
		  ry1 = ry[1][i] - L * rint(invL * ry[1][i]);
		  rz1 = rz[1][i] - L * rint(invL * rz[1][i]);
		  rxCM = (m0 * rx0 + m1 * rx1) * invMtot;
		  ryCM = (m0 * ry0 + m1 * ry1) * invMtot;;
		  rzCM = (m0 * rz0 + m1 * rz1) * invMtot;;
		  /* Sx is calculated using k parallel to x-axis, Sy with k 
		     parallel to y-axis and Sz with k parallel to the z-axis, 
		     S, the structure factor is then computed averaging over 
		     Sx, Sy and Sz */
		  siniTeta = sinL[iTeta];
		  kx = siniTeta * cosL[iPhi];
		  ky = siniTeta * sinL[iPhi];
                  kz = cosL[iTeta];
		  //printf("kx: %.15f ky: %.15f kz:%.15f\n", kx, ky, kz); 
		  rCMk = (n * scalFact + KBEG) * 
		    (rxCM * kx + ryCM * ky + rzCM * kz);
		  
		  reRho = reRho + cos(rCMk) ; 
		  /* Real part of exp(i*k*r) for the actual molecule */
		  
		  imRho = imRho + sin(rCMk) ; 
		  /* Imaginary part of exp(i*k*r) for the actual molecule*/
		  
		}
	      sumRho = sumRho + Sqr(reRho) + Sqr(imRho)
		- Nm * reRho;//  + Sqr(Nm);
	      //printf("iTeta: %d iPhi:%d sumRho: %.15f\n", iTeta, iPhi, sumRho);
	    }
	}
      //sumS[n] = sumS[n] + invNm * 1.0 / 3.0; 
      //S[n] = sumS[n] * OprogStatus.measCalc[3] / Oparams.curStep;
      //kSq = Sqr(n * scalFact + KBEG);
      S[n] = 4.0 * pi * sumRho  * invNm / Sqr(maxI);  
      /*Fattore di struttura istantaneo */
      
      /*                 2 * pi
	  k = sqrt(n) * --------  because sqrt(n) = |(n, 0, 0)| for example
	                   L  
	  where L is the box length, in our particular case L = 1. */
	 
    }
  free(cosL);
  free(sinL);
}

/* ============================= >>> maxwell <<< =========================== */
void maxwell(void)
{
  int n, i, Nm, averages;
  int* histMB;
  
  COORD_TYPE vMod, vSq, vxCM, vyCM, vzCM, m0, m1,invMtot;

  Nm = Oparams.parnum;
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  invMtot = 1.0 / (m0 + m1);
  //averages = OprogStatus.measCalc[4];
  histMB = OprogStatus.histMB;
  
  if (OprogStatus.avngMB == 0)
    {
      loop(i, 1, NUMV)
	{
	  histMB[i] = 0;
	}
    } 
      
  loop(i, 1, Nm)
    {
      vxCM = (m0 * vx[0][i] + m1 * vx[1][i]) * invMtot;
      vyCM = (m0 * vy[0][i] + m1 * vy[1][i]) * invMtot;
      vzCM = (m0 * vz[0][i] + m1 * vz[1][i]) * invMtot;
      vSq = Sqr(vxCM) + Sqr(vyCM) + Sqr(vzCM);
      vMod = sqrt(vSq);
      n = ( NUMV - 1) * (vMod - VBEG)/ (VEND - VBEG); 
      /* VBEG must be always set to 0 */
      
      if( n < NUMV ) ++histMB[n];
    }
    
  loop(i, 1, NUMV)
	{
	  MB[i] = histMB[i];
	  if (OprogStatus.avngMB == 1)
	    MB[i] /= NUMCALCS; /* Mean if needed */
	}
}

/* ========================= >>> rotDiff <<< ===============================*/
void rotDiff(void)
{
  int i, Nm;
  COORD_TYPE DoSqTot, dt;
  Nm = Oparams.parnum;
  dt = Oparams.steplength;
  DoSqTot = 0.0;
  loop(i, 1, Nm)
    {
      DoSqTot = DoSqTot + Sqr(OprogStatus.sumox[i]) + 
	Sqr(OprogStatus.sumoy[i]) + Sqr(OprogStatus.sumoz[i]);
    }
  
  /* Drot should be multiplied by 1/<1/s> because in the denominator here
     we have the virtual time, for more details see S. Nose' Mol. Phys.
     1984, vol. 52, No.2, 255-268 */
  Drot = DoSqTot / ( 4.0 * dt * ((COORD_TYPE) Oparams.curStep) * 
			((COORD_TYPE) Nm) ); 
  DoSqTot /= ((COORD_TYPE) Nm);
  DphiSq = DoSqTot;
  mdShow("Drot: %.15f", Drot);
}

/* ======================= >>> viscosity <<< ============================== */
void viscosity(void)
{
  /* DESCRIPTION:
     Calculate the shear viscosity of the molecular liquid */
  COORD_TYPE dt;
  dt = Oparams.steplength;

  eta = (Sqr(OprogStatus.DQxy) + Sqr(OprogStatus.DQyz) + Sqr(OprogStatus.DQzx))
    / ( 3.0 * dt * 2.0 * ((COORD_TYPE) Oparams.curStep) );
  
  eta = Vol * eta / Oparams.T;
  
  if (OprogStatus.avngEta == 1)
    {
      OprogStatus.sumEta += eta;
      /* NOTE:
	 Now eta is the actual mean value, this operation smooth eta, that 
	 otherwise is very noisly */
      eta = OprogStatus.sumEta / NUMCALCS; 
      /* NUMCALCS is the number of calculations perfomed for 
	 the actual measure */
    }
  
  mdShow("eta: %.15f", eta);
}

/* =========================== >>> radDens <<< ============================= */
void radDens(void)
{
  int bin, Nm, i, j;
  int* hist;

  COORD_TYPE rij, rxij, ryij, rzij, rijSq, rxCMi, ryCMi, rzCMi, rxCMj, ryCMj,
    rzCMj;
  COORD_TYPE rlower, rupper, cost, nIdeal; 
  COORD_TYPE rhoAv, m0, m1, invMtot, DELR;
  //printf("calcolo gr!!!\n");
  L = cbrt(Vol);
  invL = 1.0 / L;
  Nm = Oparams.parnum;
  rhoAv = (COORD_TYPE) Nm;
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  invMtot = 1.0 / (m0 + m1);

  DELR = (REND - RBEG) / MAXBIN;

  hist = OprogStatus.hist;
  
  if (OprogStatus.avnggr == 0)
    {
      loop(i, 1, MAXBIN)
	{
	  hist[i] = 0;
	}
    }
  
  loop(i, 1, Nm - 1)
    {
      loop(j, i+2, Nm)
	{
	  rxCMi = (m0 * rx[0][i] +  m1 * rx[1][i]) * invMtot;
	  ryCMi = (m0 * ry[0][i] +  m1 * ry[1][i]) * invMtot;
	  rzCMi = (m0 * rz[0][i] +  m1 * rz[1][i]) * invMtot;
	  rxCMj = (m0 * rx[0][j] +  m1 * rx[1][j]) * invMtot;
	  ryCMj = (m0 * ry[0][j] +  m1 * ry[1][j]) * invMtot;
	  rzCMj = (m0 * rz[0][j] +  m1 * rz[1][j]) * invMtot;
	  rxij = rxCMi - rxCMj;
	  ryij = ryCMi - ryCMj;
	  rzij = rzCMi - rzCMj;
	  rxij = rxij - L * rint(invL * rxij);
	  ryij = ryij - L * rint(invL * ryij);
	  rzij = rzij - L * rint(invL * rzij);
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  rij =  sqrt(rijSq);
	  bin = ( (int) ( (rij - RBEG) / DELR) ); 
	                                                                                                if ( (bin < MAXBIN) && (bin >= 0) )
	    {
	      hist[bin] = hist[bin] + 2;
	    }
	}
    }
  /* Normalization */
  cost = 4.0 * pi * Nm / 3.0 / Vol;
  loop(bin, 1, MAXBIN)
    {
      rlower = RBEG + ( (COORD_TYPE) bin) * DELR;
      rupper = rlower + DELR;
      nIdeal = cost * (Sqr(rupper)*rupper - Sqr(rlower)*rlower);
      gr[bin] = ((COORD_TYPE) hist[bin]) / ((COORD_TYPE) Nm) / nIdeal;
      if (OprogStatus.avnggr == 1) 
	{
	  /* This a mean value of gr[bin] */
	  gr[bin] /= NUMCALCS;
	}
    }
}

/* ========================= >>> Ptensor <<< =============================== */
void Ptensor(void)
{
  /* DESCRIPTION:
     Store the three off-diagonal terms of the pressure tensor in an array 
     to save on disk like a single measure */
  Ptens[0] = Pxy;
  Ptens[1] = Pyz;
  Ptens[2] = Pzx;

}

/* ========================== >>> DQtensor <<< ============================= */
void DQtensor(void)
{
  DQtens[0] = OprogStatus.DQxy;
  DQtens[1] = OprogStatus.DQyz;
  DQtens[2] = OprogStatus.DQzx;
}

/* ======================= >>> corrFuncs <<< =============================== */
void corrFuncs(void)
{
  /* DESCRIPTION:
     Calculate auto-correlation functions of single particle quantity,
     averaged on all particles */
  COORD_TYPE o0, ott0, *ox0, *oy0, *oz0, dotu, Mtot, m0, m1;
  COORD_TYPE vcmx, vcmy,vcmz, *ux0, *uy0, *uz0, *vcmx0, *vcmy0, *vcmz0, doto;
  int i, Nm;
  
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  Mtot = m0 + m1;
  
  ox0 = OprogStatus.ox0;
  oy0 = OprogStatus.oy0;
  oz0 = OprogStatus.oy0;
  ux0 = OprogStatus.ux0;
  uy0 = OprogStatus.uy0;
  uz0 = OprogStatus.uz0;
  vcmx0 = OprogStatus.vcmx0;
  vcmy0 = OprogStatus.vcmy0;
  vcmz0 = OprogStatus.vcmz0;

  Nm = Oparams.parnum;

  /* Various correlators */
  psi1c = 0.0;
  psi2c = 0.0;
  C1c = 0.0;
  C2c = 0.0;
  C3c = 0.0;
  C4c = 0.0;
  velc = 0.0;
  
  loop(i, 1, Nm)
    {
      /* Correlator for angular velocities */
      o0 = sqrt(Sqr(ox0[i]) + Sqr(oy0[i]) + Sqr(oz0[i]));
      ott0 = sqrt(Sqr(ox[i]) + Sqr(oy[i]) + Sqr(oz[i]));
      doto = ox0[i] * ox[i] + oy0[i] * oy[i] +
	oz0[i] * oz[i];
      doto /= o0 * ott0; 
      psi1c += doto;
      psi2c += 0.5 * ( 3.0 * Sqr(doto) - 1.0);
      
      /* Correlator for orientatios */
      dotu = ux0[i] * ux[i] + uy0[i] * uy[i] +
	uz0[i] * uz[i];
      C1c += dotu;
      C2c += 0.5 * (3.0 * Sqr(dotu) - 1.0);
      C3c += 0.5 * (5.0 * Sqr(dotu) * dotu - 
		      3.0 * dotu);
      C4c += 0.125 * (35.0 * Sqr(dotu) * Sqr(dotu) - 
			30.0 * Sqr(dotu) + 3.0);
     
      /* Velocity correlation function */
      vcmx = (m0 * vx[0][i] + m1 * vx[0][i]) / Mtot;
      vcmy = (m0 * vy[0][i] + m1 * vy[0][i]) / Mtot;
      vcmz = (m0 * vz[0][i] + m1 * vz[0][i]) / Mtot;
      
      velc += vcmx0[i] * vcmx + vcmy0[i] * vcmy +
	vcmz0[i] * vcmz;
    }
  
  psi1c /= (COORD_TYPE) Nm;
  psi2c /= (COORD_TYPE) Nm;
  C1c /= (COORD_TYPE) Nm;
  C2c /= (COORD_TYPE) Nm;
  C3c /= (COORD_TYPE) Nm;
  C4c /= (COORD_TYPE) Nm;
  velc /= (COORD_TYPE) Nm;
}
/* ============================ >>> vanHove <<< ============================ */
void vanHoveSelf(void)
{
  int i, j, Nm;
  COORD_TYPE rxcm, rycm, rzcm, deltar, r;
  COORD_TYPE *rxcm0, *rycm0, *rzcm0; /* Initial center of mass positions */
  int Gsint[GSPOINT];

  rxcm0 = OprogStatus.rxCMi;
  rycm0 = OprogStatus.ryCMi;
  rzcm0 = OprogStatus.rzCMi;

  Nm = Oparams.parnum;
 
  /* Initilization */
  loop(j, 1, GSPOINT)
    {
      Gsint[j] = 0;
    }
  /* Gsint counts the number of particles, that have travelled a certain
     distance at the actual time */
  loop(i, 1, Nm)
    {
      CoM(i, &rxcm, &rycm, &rzcm);
      /* deltar is the distance travelled by the actual particles center of 
	 mass */
      deltar = sqrt(Sqr(rxcm - rxcm0[i]) + Sqr(rycm - rycm0[i]) +
	Sqr(rzcm - rzcm0[i]));
      j = (int) (GSPOINT * deltar / GSRMAX); 
      if (j < GSPOINT) ++Gsint[j];
    }

  /* Normalization of the van Hove function */
  loop(j, 1, GSPOINT)
    {
      Gs[j] = ((COORD_TYPE)Gsint[j]) / Nm; /* Self-Part of the Van Hove 
					      function */
      r = ((COORD_TYPE) j + 0.5) * (GSRMAX / GSPOINT);
      Gs[j] *= GSPOINT / GSRMAX;/* Gs is a density of probability */
      Gs[j] /= Sqr(r) * 4.0 * pi;
    }
}

/* ========================== >>> vanHoveAng << =========================== */
void vanHoveAng(void)
{
  COORD_TYPE dotu, teta;
  COORD_TYPE *ux0, *uy0, *uz0;
  int i, tetai, Nm;
  const COORD_TYPE dteta = 0.0000001; 

  ux0 = OprogStatus.ux0;
  uy0 = OprogStatus.uy0;
  uz0 = OprogStatus.uz0;
  Nm = Oparams.parnum;

  /* Initializa GsAng vector */
  loop(i, 1, GSANGPOINT)
    {
      GsAng[i] = 0.0;
    }

  loop(i, 1, Nm)
    {
      dotu = ux0[i] * ux[i] + uy0[i] * uy[i] +
	uz0[i] * uz[i];
      
      teta = acos(dotu);
      
      if ( (teta > dteta) && (teta < (pi - dteta)) )
	{
	  /* nTeta is the number of angles for which we 
	     calculate the van Hove function G = G(t, teta) */
	  tetai = (int) (GSANGPOINT * acos(dotu) / pi);
	  GsAng[tetai] += 1.0; 
	}
    }

  /* Normalization of the Van Hove function */
  for (i = 0; i < GSANGPOINT; ++i) 
    {
      /* The interval between 0 and pi is divided in nTeta subintervals,
	 so the for each interval the corrisponding teta is the point
	 in the middle */
      teta = ((COORD_TYPE)i + 0.5) * (pi / GSANGPOINT); 
      /* 0 < teta < pi */  
      GsAng[i] /= ((COORD_TYPE) Nm) * sin(teta);
      
      /* In this way when t -> +oo G(teta, t) = 1 for all teta */ 
      GsAng[i] *= 2 * GSANGPOINT / pi;
    }
}
/* ========================== >>> intermScat <<< ===========================*/
void intermScatt(void)
{
  int kInt, rInt;
  COORD_TYPE k, r;

  /* Fs is the fourier transform of the van Hove function */
  loop(kInt, 1, FSPOINT)
    {
      Fs[kInt] = 0.0;
      k = ( ((COORD_TYPE)kInt) + 0.5) * ( FSKMAX / ((COORD_TYPE)FSPOINT) );
      loop(rInt, 1, GSPOINT)
	{
	  r = ( ((COORD_TYPE) rInt) + 0.5) * (GSRMAX / ((COORD_TYPE)GSPOINT));
	  Fs[kInt] += 4.0 * pi * sin(k*r) * Gs[rInt] * r / k;
	}
    }
}
/* ============================ >>> GsvsGaus <<< ========================== */
void gaussapprox(void)
{
  /* Calculate
     (Gs - Gsguass) / GsGauss where GsGauss is 
     GsGauss
     DEPENDENCIES:
     - transdiff()
     - vanHoveSelf()
  */
  int i;
  COORD_TYPE Gsg, r, prefact;
  prefact = pow(3.0 / DrSqTot / pi / 2.0, 1.5);
  loop(i, 1, GSPOINT)
    {
      r = ((COORD_TYPE) i + 0.5) * (GSRMAX / GSPOINT); 
      Gsg = prefact * exp(-3.0*Sqr(r)/(2.0*DrSqTot));
      /* Difference from a gaussian */
      GsGsg[i] = (Gs[i] - Gsg) / Gsg;
    }
}
