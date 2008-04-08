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

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE V, W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx;  

extern COORD_TYPE Mtot;
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */


/* ========================================================================= */

/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
  COORD_TYPE Px, Py, Pz, cost, RCMx, RCMy, RCMz;
  int mol, Nm, i, a;
  COORD_TYPE L, invL, MTOT;
  double *m;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  mdShow("SHIFTED POTENTIAL ENERGY: %.10f", Vc);
  printf("KINETIC ENERGY: %.10f\n", K);
  printf("s: %.10f s1: %.10f Vol: %.10f Vol1: %.10f\n", s, s1, Vol, Vol1);
  E = K + Vc;
  /* And now add the contribute due to the thermal bath */
  E = E + 0.5 * Sqr(s1) * OprogStatus.Q / Sqr(s) + 
    (3.0 * ((COORD_TYPE) Nm) - 3.0) * Oparams.T * log(s) + 
    0.5 * Sqr(Vol1) * OprogStatus.W / Sqr(s) + Oparams.P * Vol; 
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */

  printf("TOTAL ENERGY: %.10f\n", E);
  printf("Elrc:: %.6f Plrc: %.6f\n", Elrc, Plrc);
  mol = 10;
  printf("Atom position: %.15f\n", rx[0][mol]);

  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  cost = Vol1 / Vol / 3.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  //L = cbrt(Vol);
  //invL = 1.0 / L;
  
  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  Px += Oparams.m[a] * (vx[a][i] - 
	    cost * rx[a][i]);
	  Py += Oparams.m[a] * (vy[a][i] - 
	    cost * ry[a][i]);
	  Pz += Oparams.m[a] * (vz[a][i] - 
	    cost * rz[a][i]);
	  //vectProd(rx[a][i], ry[a][i], rz[a][i], 
	       //     vx[a][i], vy[a][i], vz[a][i],
          //&lx, &ly, &lz);
	  //Lx += lx;
	  //Ly += ly;
	  //Lz += lz;
	}
    }
  printf("Px: %.15f Py: %.15f Pz: %.15f\n", Px, Py, Pz);
  
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  MTOT = 0.0;

  m = Oparams.m;

  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{  
	  /* ATTENZIONE: Se le masse delle due specie 
	     sono diverse cambiare il codice QUI!!!!!!!!!! */
	  RCMx += m[a]*rx[a][i];
	  RCMy += m[a]*ry[a][i];
	  RCMz += m[a]*rz[a][i];
	  MTOT += m[a];
	}
    }
  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;
  printf("Box CoM RCMx: %.15f RCMy: %.15f RCMz: %.15f\n", RCMx, RCMy, RCMz);
}

/* =========================== >>> press_ex <<< ============================ */
void press_ex(void)
{
  /* Questa routine calcola la pressione in eccesso rispetto al gas ideale */
  /*Pex = press - (Oparams.parnum[0] + Oparams.parnum[1]) * Oparams.T / Vol;*/
  Pex = W/Vol;
}

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  int a, i;
  double DrSq=0.0, Drx, Dry, Drz;
  FILE *f;  
  for (a = 0; a < 2; a++)
    {
      DrSqTot[a] = 0.0;
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  Drx = rx[a][i] - OprogStatus.rxi[a][i];
	  Dry = ry[a][i] - OprogStatus.ryi[a][i];
	  Drz = rz[a][i] - OprogStatus.rzi[a][i];
	  DrSqTot[a] = DrSqTot[a] + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
	}
      Dtrans[a] = DrSqTot[a] / ( 6.0 * ((double) Oparams.steplength) *
	     		    ((double) Oparams.curStep) * 
			    ((double) Oparams.parnum[a] ) );   
      DrSq += DrSqTot[a];
      DrSqTot[a] /= ((double) Oparams.parnum[a]);
    }
  f=fopen("msd.dat", "a");
  fprintf(f, "%.15G %.15G\n", Oparams.steplength*Oparams.curStep, 
	  DrSq/((double)Oparams.parnum[0]+Oparams.parnum[1]));
  fclose(f);
  if (Oparams.parnum[1]!=0)
    {
      f=fopen("msdA.dat", "a");
      fprintf(f, "%.15G %.15G\n", Oparams.steplength*Oparams.curStep, 
	      DrSqTot[0]);
      fclose(f);
      f=fopen("msdB.dat", "a");
      fprintf(f, "%.15G %.15G\n", Oparams.steplength*Oparams.curStep, 
	      DrSqTot[1]);
      fclose(f);
    }
  //printf("Dtr[0]: %.6f Dtr[1]: %.6f\n", Dtrans[0], Dtrans[1]);
}

/* ============================ >>> temperat <<< =========================== */
void temperat(void)
{
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
  temp = 2.0 * K / (3.0 * 
		    (((double)Oparams.parnum[0]+Oparams.parnum[1])) - 3.0);

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
#if 0
  printf("Instantaneous temperature:%.15f\n", temp);
  printf("W: %.15f\n", W);
  //printf("Molecular pressure: %f\n", press_m);


  printf("Pressure: %.15f\n", press);
  printf("Vol: %.15f\n", Vol);
#endif
}
#ifdef MD_LOADMESH
extern int mesh[100][150][3];
extern int ntripl[100];
#endif
/* ======================== >>> structFacts <<< =============================*/
void structFacts(void)
{
  /* Riscrivere usando il codice del monoatomico!!!!!!!!!!!!!!!!!!!! */
 /* DESCRIPTION:
     This mesuring function calculates the static structure factor */
  COORD_TYPE reRho, imRho, pi2;
  COORD_TYPE invNm, scalFact;
  COORD_TYPE rCMk, sumRho;
  COORD_TYPE *sumS, kbeg;
#ifndef MD_LOADMESH
  int mesh[][150][3]=
#include "./kmesh.dat"
  int ntripl[]=
#include "./ntripl.dat"
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
  Nm = Oparams.parnum[0];
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
	      rCMk = kbeg + scalFact * 
		(rx[0][i] * mesh[n][mp][0] + ry[0][i] * mesh[n][mp][1] + 
		 rz[0][i] * mesh[n][mp][2]);
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
  /* Distribuzione MB per gli atomi di tipo A (0) */
  int n, i, Nm, averages;
  int* histMB;
  
  COORD_TYPE vMod, vSq, m0;

  Nm = Oparams.parnum[0];
  m0 = Oparams.m[0];
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
      vSq = Sqr(vx[0][i]) + Sqr(vy[0][1]) + Sqr(vz[0][1]);
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
  /* Distribuzione radial per gli atomi di tipo A */
  int bin, Nm, i, j;
  int* hist;

  COORD_TYPE rij, rxij, ryij, rzij, rijSq;
  COORD_TYPE rlower, rupper, cost, nIdeal; 
  COORD_TYPE rhoAv, m0, DELR;
  //printf("calcolo gr!!!\n");
  L = cbrt(Vol);
  invL = 1.0 / L;
  Nm = Oparams.parnum[0];
  rhoAv = (COORD_TYPE) Nm;
  m0 = Oparams.m[0];

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
	  rxij = rx[0][i] - rx[0][j];
	  ryij = ry[0][i] - ry[0][j];
	  rzij = rz[0][i] - rz[0][j];
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

