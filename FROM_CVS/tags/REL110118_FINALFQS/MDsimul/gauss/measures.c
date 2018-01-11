#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

#ifdef MPI
extern int my_rank;
#endif
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
char TXTA[10][MSG_LEN];
char TXT[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz;  

/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

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
  COORD_TYPE Px, Py, Pz, cost, RCMx, RCMy, RCMz;
  int mol, Nm, i;
  COORD_TYPE px, py, pz;
  COORD_TYPE L, invL;
  FILE* mf;

  Nm = Oparams.parnum;
  sprintf(TXTA[0],"STEP %d [MEASURE]\n  s: %.10f s1: %.10f Vol: %.10f Vol1: %.10f\n", 
	  Oparams.curStep, s, s1, Vol, Vol1);
  

  E = K + Vc;
  /* And now add the contribute due to the thermal bath */
  E = E + 0.5 * Sqr(s1) * OprogStatus.Q / Sqr(s) + 
    (3.0 * ((COORD_TYPE) Nm) - 3.0) * Oparams.T * log(s) + 
    0.5 * Sqr(Vol1) * OprogStatus.W / Sqr(s) + Oparams.P * Vol; 
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */
  mol = 10;
  //printf("Elrc:: %.6f Plrc: %.6f\n", Elrc, Plrc);
  //printf("Atom position: %.15f\n", rx[mol]);
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  cost = Vol1 / Vol / 3.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  //L = cbrt(Vol);
  //invL = 1.0 / L;

  loop(i, 1, Nm)
    {
      
      // printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
     
      px = Oparams.m * (vx[i] - 
			    cost * rx[i]);
      py = Oparams.m * (vy[i] - 
			    cost * ry[i]);
      pz = Oparams.m * (vz[i] - 
			    cost * rz[i]);
      //vectProd(rx[a][i], ry[a][i], rz[a][i], 
      //     vx[a][i], vy[a][i], vz[a][i],
          //&lx, &ly, &lz);
      //Lx += lx;
      //Ly += ly;
      //Lz += lz;
      Px += px;
      Py += py;
      Pz += pz;
    }
  sprintf(TXTA[1], "  E=%.15f P=(%.15f,%.15f,%.15f)\n", E, Px, Py, Pz);
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  
  loop(i, 1, Nm)
    {
      RCMx += rx[i];
      RCMy += ry[i];
      RCMz += rz[i];
    }
 
  sprintf(TXTA[2],"  BOX CM=(%.15f,%.15f,%.15f)\n", RCMx, RCMy, RCMz);

  //printf("RANK: %d STEP: %d\n", my_rank, Oparams.curStep);
  //fflush(stdout);
  mdPrintf(ALL, TXTA[0], TXTA[1], TXTA[2], NULL);

  mf = fopen(absMisHD("Vc.a"),"a");
  fprintf(mf, "%d %.15G\n", Oparams.curStep, Vc);
  fclose(mf);
}

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  COORD_TYPE Drx, Dry, Drz, Dr4;
  int i;
  
  DrSqTot = 0.0;
  Dr4 = 0.0;

  loop(i, 1, Oparams.parnum)
    {
      Drx = rx[i] - OprogStatus.rxCMi[i]; 
      Dry = ry[i] - OprogStatus.ryCMi[i];
      Drz = rz[i] - OprogStatus.rzCMi[i];
      if (OprogStatus.ipart == i)
	{
	  //sprintf(TXT,"i = %d\n", i);
	  //mdPrintf(STD, TXT, NULL);
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
  
}

/* ============================ >>> temperat <<< =========================== */
void temperat(void)
{
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
  COORD_TYPE Ktransl, m;
  temp = 2.0 * K / (3.0 * Oparams.parnum - 3.0);

  Ktransl = 0.0;
  m = Oparams.m;
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
  sprintf(TXT, "P:%.10f T:%.10f W: %10f\n", press, temp, W);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);
}

/* ======================== >>> structFacts <<< =============================*/
void structFacts_OLD(void)
{
  int n, bin;
  COORD_TYPE sum, pi4, rhoAv, r;
  COORD_TYPE scalFact, kr, k, DELR;
  COORD_TYPE* sumS;
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

/* =========================== >>> structFacta <<< =========================*/
void structFacts(void)
{
  /* DESCRIPTION:
     This mesuring function calculates the static structure factor */
  COORD_TYPE reRho, imRho, pi2;
  COORD_TYPE invNm, scalFact;
  COORD_TYPE rCMk, sumRho;
  COORD_TYPE *sumS, kbeg;
  int mp, i, Nm, n;
  int mesh[][150][3]= 
#include "./kmesh.dat"
  int ntripl[]=
#include "./ntripl.dat"

  /* useful quantities */
  sumS = OprogStatus.sumS;
  L = cbrt(Vol);
  invL = 1.0 / L;
  //scalFact = (KEND - KBEG) / NUMK;
  
  pi2 = 2.0 * pi;
  kbeg = 0.0; //pi2 * invL;
  scalFact = pi2 * invL;
  //printf("maxI : %d\n", maxI);
  /* We take 'maxI' values of Phi and 'maxI' values of Teta, so we have in this
     way NUMK2AV = maxI * maxI (about) values for 'k' over the 
     semi-sphere with same modulus */
  
  Nm = Oparams.parnum;
  invNm = 1.0 / Nm;

  loop(n, 1, NUMK)
    {
      sumRho = 0.0;
      //printf("nummp:%d\n", ntripl[n]);      
      loop(mp, 1, ntripl[n])
	{
	  reRho = 0.0;
	  imRho = 0.0;
	  loop(i, 1, Oparams.parnum)
	    {
	      // il passo della mesh e' 0.5*pi2/L
	      if (mesh[n][mp][0]==0 && mesh[n][mp][1] == 0 && 
		  mesh[n][mp][2] == 0)
		{
		  printf("ERRORE nella MESH!!!!!!!!\n");
		  exit(-1);
		}
	      rCMk = kbeg + scalFact * 
		(rx[i] * mesh[n][mp][0] + ry[i] * mesh[n][mp][1] + 
		 rz[i] * mesh[n][mp][2]);
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
      loop(n, 1, NUMK)
	{
	  sumS[n] += S[n];
	  //if (n == 10) printf("somma di S:%f\n", sumS[n]);
	  S[n] = sumS[n] / NUMCALCS;
	}
    }
}


/* ============================= >>> maxwell <<< =========================== */
void maxwell(void)
{
  int n, i, Nm;
  int* histMB;
  
  COORD_TYPE vMod, vSq, m;

  Nm = Oparams.parnum;
  m = Oparams.m;
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
      vSq = Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]);
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

/* =========================== >>> radDens <<< ============================= */
void radDens(void)
{
  int bin, Nm, i, j;
  int* hist;

  COORD_TYPE rij, rxij, ryij, rzij, rijSq;
  COORD_TYPE m, rlower, rupper, cost, nIdeal; 
  COORD_TYPE rhoAv, DELR;
  L = cbrt(Vol);
  invL = 1.0 / L;
  Nm = Oparams.parnum;
  rhoAv = (COORD_TYPE) Nm;
  m= Oparams.m;
  
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
	  rxij = rx[i] - rx[j];
	  ryij = ry[i] - ry[j];
	  rzij = rz[i] - rz[j];
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


