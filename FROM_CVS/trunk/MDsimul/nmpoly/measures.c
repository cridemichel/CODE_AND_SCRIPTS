#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
char TXT[MSG_LEN], TXTA[10][MSG_LEN];
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
extern COORD_TYPE W, K, WC, T1xx, T1yy, T1zz,
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
extern void calcEnePot(int Nm, double rcut);

/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
  COORD_TYPE dist, Px, Py, Pz, cost, Rx, Ry, Rz, RCMx, RCMy, RCMz, Mtot;
  int n, mol, Nm, i, a, dof;
  COORD_TYPE px, py, pz, lx, ly, lz, Lx, Ly, Lz;
  COORD_TYPE Rx_scal, Ry_scal, Rz_scal;
  COORD_TYPE L, invL;

  Nm = Oparams.parnum;
  //calcEnePot(Oparams.parnum, Oparams.rcut);
  mdShow("SHIFTED POTENTIAL ENERGY: %.10f", Vc);
  mdShow("POTENTIAL ENERGY: %.10f", V);
  printf("KINETIC ENERGY: %.10f\n", K);
  printf("s: %.10f s1: %.10f Vol: %.10f Vol1: %.10f\n", s, s1, Vol, Vol1);
  E = K + Vc;
  /* And now add the contribute due to the thermal bath */
#ifdef MD_FENE
  printf("FENE energy Vfe:%f\n", Vfe);
  E += Vfe;
  dof = 3.0*NA;
#else
  dof = 2.0*NA+1;
#endif
  
  if (OprogStatus.Nose > 0)
    {
      if (OprogStatus.Nose == 2 || OprogStatus.Nose == 1)
	{
	  E = E + 0.5 * Sqr(s1) * OprogStatus.Q / Sqr(s) + 
	    ( dof * ((COORD_TYPE) Nm) - 3.0) * Oparams.T * log(s);
	}
      if (OprogStatus.Nose == 1)
	E = E + 
	  0.5 * Sqr(Vol1) * OprogStatus.W / Sqr(s) + Oparams.P * Vol; 
      /* So now E is the extended hamiltonian that should be an integral of 
	 motion */
    }
#if 0
  else if (OprogStatus.Nose == -2)
    {
      /* -2 = dinamica browniana + Andersen (s=1) */
      E = E + Oparams.P * Vol + 0.5 * Sqr(Vol1) * OprogStatus.W;
    }
#endif
  printf("TOTAL ENERGY: %.10f\n", E);
  printf("Elrc:: %.6f Plrc: %.6f\n", Elrc, Plrc);
#if 0
  mol = 0;
  dist = sqrt( Sqr(rx[0][mol] - rx[1][mol]) + Sqr(ry[0][mol] - ry[1][mol]) + 
	       Sqr(rz[0][mol] - rz[1][mol]) );
  printf("Atom position: %.15f\n", rx[0][mol]);
  printf("Molecule atoms distance: %.15f\n", dist);
  dist = ((rx[0][mol] - rx[1][mol]) * (vx[0][mol] - vx[1][mol]) + 
    (ry[0][mol] - ry[1][mol]) * (vy[0][mol] - vy[1][mol]) + 
    (rz[0][mol] - rz[1][mol]) * (vz[0][mol] - vz[1][mol]))/dist;
  printf("Molecule relative velocity along bond: %.15f\n", dist); 
#endif
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  Lx = 0.0;
  Ly = 0.0;
  Lz = 0.0;
  Mtot=0;
  for (n=0; n < NA; n++)
    Mtot += Oparams.m[n];
  cost = Vol1 / Vol / 3.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx, &Ry, &Rz);
      
      px = 0.0;
      py = 0.0;
      pz = 0.0;
      
      for(a=0; a < NA; a++)
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
  
  for(i=0; i < Nm; i++)
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
  COORD_TYPE *m;
  COORD_TYPE rxCM, ryCM, rzCM; /* Center of mass position */
  int i, a;
  
  m = Oparams.m;
  Mtot = 0;
  for (a=0; a < NA; a++)
    Mtot += m[a];
  invMtot = 1.0 / Mtot;
 
  DrSqTot = 0.0;
  Dr4 = 0.0;

  for(i=0; i < Oparams.parnum; i++)
    {
      rxCM = ryCM = rzCM = 0.0;
      for (a=0; a < NA; a++)
	{
	  rxCM += m[a] * rx[a][i];
	  ryCM += m[a] * ry[a][i];
	  rzCM += m[a] * rz[a][i];
	}
      rxCM *= invMtot;
      ryCM *= invMtot;
      rzCM *= invMtot;
      Drx = rxCM - OprogStatus.rxCMi[i] + OprogStatus.DR[i][0]; 
      Dry = ryCM - OprogStatus.ryCMi[i] + OprogStatus.DR[i][0];
      Drz = rzCM - OprogStatus.rzCMi[i] + OprogStatus.DR[i][0];
      /*printf("step: %d DR: (%f,%f,%f)\n", Oparams.curStep, Drx, Dry, Drz);*/
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
  COORD_TYPE invMtot,Mtot, Ktransl, vxCM, vyCM, vzCM, m0, m1, m2;
  int dof;

#ifdef MD_FENE
  dof = 3*NA;
#else
  dof = (NA*2.0+1); 
#endif
  if (OprogStatus.Nose >= 0)
    temp = 2.0 * K / (dof * Oparams.parnum - 3.0);
  else
    temp = 2.0 * K / (dof * Oparams.parnum);
  
  Ktransl = 0.0;
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  m2 = Oparams.m[2];
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
#ifdef MD_LOADMESH
extern int mesh[100][150][3];
extern int ntripl[100];
#endif

/* ======================== >>> structFacts <<< =============================*/
void SkCM(void)
{
  /* Riscrivere usando il codice del monoatomico!!!!!!!!!!!!!!!!!!!! */
 /* DESCRIPTION:
     This mesuring function calculates the static structure factor */
  COORD_TYPE reRho, imRho, pi2;
  COORD_TYPE invNm, scalFact;
  COORD_TYPE rCMk, sumRho, RCMx, RCMy, RCMz;
  COORD_TYPE *sumS, kbeg;
#ifndef MD_LOADMESH
  int mesh[][150][3]=
#include "../bimix/kmesh.dat"
  int ntripl[]=
#include "../bimix/ntripl.dat"
#endif
  int mp, i, Nm, n;
  const int numk = 98;
  /* useful quantities */
  sumS = OprogStatus.sumS;
  L = cbrt(Vol);
  invL = 1.0 / L;
  /* scalFact = (KEND - KBEG) / NUMK; */
  
  pi2 = 2.0 * pi;
  kbeg = 0.0; /*pi2 * invL;*/
  scalFact = pi2 * invL;
  Nm = Oparams.parnum;
  invNm = 1.0 / Nm;

  for(n = 0; n < numk; n++)
    {
      sumRho = 0.0;
      for(mp = 0; mp < ntripl[n]; mp++)
	{
	  reRho = 0.0;
	  imRho = 0.0;
	  for(i=0; i < Nm; i++)
	    {
	      /* il passo della mesh e' 0.5*pi2/L */
	      if (mesh[n][mp][0]==0 && mesh[n][mp][1] == 0 && 
		  mesh[n][mp][2] == 0)
		{
		  printf("ERRORE nella MESH!!!!!!!! n=%d mp=%d ntripl[n]:%d\n", n,
			 mp, ntripl[n]);
		  exit(-1);
		}
	      /* calculate structure factors*/ 
	      CoM(i, &RCMx, &RCMy, &RCMz);	
	      rCMk = kbeg + scalFact * 
		(RCMx * mesh[n][mp][0] + RCMy * mesh[n][mp][1] + 
		 RCMz * mesh[n][mp][2]);
	      reRho = reRho + cos(rCMk) ; 
	      imRho = imRho + sin(rCMk);
	      /* Imaginary part of exp(i*k*r) for the actual molecule*/
	    }
	  sumRho = sumRho + Sqr(reRho) + Sqr(imRho);
	}

      S[n] = sumRho  * invNm / ((COORD_TYPE) ntripl[n]);  
    }
  if (OprogStatus.avngS == 1)
    {
      for(n = 0; n < numk; n++)
	{
	  sumS[n] += S[n];
	  /*if (n == 10) printf("somma di S:%f\n", sumS[n]);*/
	  S[n] = sumS[n] / NUMCALCS;
	}
    }

  
}

/* ============================= >>> maxwell <<< =========================== */
void maxwell(void)
{
  int n, i, Nm, averages;
  int* histMB;
  
  COORD_TYPE vMod, vSq, vxCM, vyCM, vzCM, m0, m1, m2, invMtot;

  Nm = Oparams.parnum;
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  m2 = Oparams.m[2];
  invMtot = 1.0 / (m0 + m1 + m2);
  //averages = OprogStatus.measCalc[4];
  histMB = OprogStatus.histMB;
  
  if (OprogStatus.avngMB == 0)
    {
      for(i=0; i <  NUMV; i++)
	{
	  histMB[i] = 0;
	}
    } 
      
  for(i=0; i < Nm; i++)
    {
      vxCM = (vx[0][i]*m0 + vx[1][i]*m1 + vy[2][i]*m2) * invMtot;
      vyCM = (vy[0][i]*m0 + vy[1][i]*m1 + vy[2][i]*m2) * invMtot;
      vzCM = (vz[0][i]*m0 + vz[1][i]*m1 + vz[2][i]*m2) * invMtot;
      vSq = Sqr(vxCM) + Sqr(vyCM) + Sqr(vzCM);
      vMod = sqrt(vSq);
      n = ( NUMV - 1) * (vMod - VBEG)/ (VEND - VBEG); 
      /* VBEG must be always set to 0 */
      
      if( n < NUMV ) ++histMB[n];
    }
    
  for(i=0; i < NUMV; i++)
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
  for(i=0; i < Nm; i++)
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
  int bin, Nm, a, i, j;
  int* hist;
  COORD_TYPE rij, rxij, ryij, rzij, rijSq, rxCMi, ryCMi, rzCMi, rxCMj, ryCMj,
    rzCMj, *m;
  COORD_TYPE rlower, rupper, cost, nIdeal; 
  COORD_TYPE rhoAv, invMtot, DELR, Mtot;
  double rend, rbeg;
  L = cbrt(Vol);
  invL = 1.0 / L;
  Nm = Oparams.parnum;
  rhoAv = (COORD_TYPE) Nm;
  m = Oparams.m;
 

  Mtot = 0.0;
  for (a = 0; a < NA; a++)
    Mtot += m[a];
  
  invMtot = 1.0 / Mtot;
  
  rbeg = 0.2;
  rend = L * 0.5;
  DELR = (rend - rbeg) / MAXBIN;

  hist = OprogStatus.hist;
  
  if (OprogStatus.avnggr == 0)
    {
      for(i=0; i < MAXBIN; i++)
	{
	  hist[i] = 0;
	}
    }
  
  for(i=0; i < Nm - 1; i++)
    {
      for(j=i+1; j < Nm; j++)
	{
	  rxCMi = ryCMi = rzCMi = 0;
	  rxCMj = ryCMj = rzCMj = 0;
	  for (a=0; a < NA; a++)
	    {
	      rxCMi += m[a]*rx[a][i];
	      ryCMi += m[a]*ry[a][i];
	      rzCMi += m[a]*rz[a][i];
	      rxCMj += m[a]*rx[a][j];
	      ryCMj += m[a]*ry[a][j];
	      rzCMj += m[a]*rz[a][j];
	    }
	  rxCMi *= invMtot;
	  ryCMi *= invMtot;
	  rzCMi *= invMtot;
	  rxCMj *= invMtot;
	  ryCMj *= invMtot;
	  rzCMj *= invMtot;
	  rxij = rxCMi - rxCMj;
	  ryij = ryCMi - ryCMj;
	  rzij = rzCMi - rzCMj;
	  rxij = rxij - L * rint(invL * rxij);
	  ryij = ryij - L * rint(invL * ryij);
	  rzij = rzij - L * rint(invL * rzij);
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  rij =  sqrt(rijSq);
	  bin = ( (int) ( (rij - rbeg) / DELR) );
	  if ( (bin < MAXBIN) && (bin >= 0) )
	    {
	      hist[bin] = hist[bin] + 2;
	    }
	}
    }
  /* Normalization */
  cost = 4.0 * pi * Nm / 3.0 / Vol;
  for(bin=0; bin < MAXBIN; bin++)
    {
      rlower = rbeg + ( (COORD_TYPE) bin) * DELR;
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
void radDens23(void)
{
  /* fill in!! */ 
}
  
void radDens33(void)
{ 
  /* fill in!!*/
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

#if 0
/* ======================= >>> corrFuncs <<< =============================== */
void corrFuncs(void)
{
  /* DESCRIPTION:
     Calculate auto-correlation functions of single particle quantity,
     averaged on all particles */
  COORD_TYPE o0, ott0, *ox0, *oy0, *oz0, dotu, Mtot, m0, m1, m2;
  COORD_TYPE vcmx, vcmy,vcmz, *ux0, *uy0, *uz0, *vcmx0, *vcmy0, *vcmz0, doto;
  int i, Nm;
  
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  m2 = Oparams.m[2];
  Mtot = m0 + m1 + m2;

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
  
  for(i = 0; i < Nm; i++)
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
      vcmx = (m0*vx[0][i] + m1*vx[1][i] + m2*vx[2][i]) / Mtot;
      vcmy = (m0*vy[0][i] + m1*vy[1][i] + m2*vy[2][i]) / 3.0;
      vcmz = (m0*vz[0][i] + m1*vz[1][i] + m2*vz[2][i]) / 3.0;
      
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
  for(j=0; j < GSPOINT; j++)
    {
      Gsint[j] = 0;
    }
  /* Gsint counts the number of particles, that have travelled a certain
     distance at the actual time */
  for(i=0; i < Nm; i++)
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
  for(j=0; j < GSPOINT; j++)
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
  for(i=0; i < GSANGPOINT; i++)
    {
      GsAng[i] = 0.0;
    }

  for(i=0; i < Nm; i++)
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
  for(kInt=0; kInt < FSPOINT; kInt++)
    {
      Fs[kInt] = 0.0;
      k = ( ((COORD_TYPE)kInt) + 0.5) * ( FSKMAX / ((COORD_TYPE)FSPOINT) );
      for(rInt=0; rInt < GSPOINT; rInt++)
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
  for(i=0; i < GSPOINT; i++)
    {
      r = ((COORD_TYPE) i + 0.5) * (GSRMAX / GSPOINT); 
      Gsg = prefact * exp(-3.0*Sqr(r)/(2.0*DrSqTot));
      /* Difference from a gaussian */
      GsGsg[i] = (Gs[i] - Gsg) / Gsg;
    }
}
#endif
