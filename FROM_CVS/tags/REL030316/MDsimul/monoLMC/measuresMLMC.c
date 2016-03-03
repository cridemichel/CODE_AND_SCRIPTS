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
extern int *rxS, *ryS, *rzS;
extern int* pmap;
/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE K, Wxx, Wyy, Wzz, Wxy, Wyz, Wzx;  

extern COORD_TYPE Mtot;
/* used by linked list routines */

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax;
/* ================================= */


/* ========================================================================= */

/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
  COORD_TYPE RCMx, RCMy, RCMz, MTOT;
  int mol, Nm, i, a;
  COORD_TYPE m, la, L, invL;
  
  Nm = Oparams.parnum;
  mdShow("POTENTIAL ENERGY: %.10f", V);
  printf("VLJ: %f\n", VLJ);
  mol = 10;
  printf("temp: %.6f Atom position: %d\n", Oparams.T, rx[mol]);


  L = cbrt(Vol);
  invL = 1.0 / L;

  //L = cbrt(Vol);
  //invL = 1.0 / L;
  
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  la = Oparams.lattice_a;
  m = Oparams.m;

  MTOT = 0.0;

  for (i = 0; i < Oparams.parnum; i++)
    {  
      /* ATTENZIONE: Se le masse delle due specie 
	 sono diverse cambiare il codice QUI!!!!!!!!!! */
      RCMx += m*la * ((double)rx[i]);
      RCMy += m*la * ((double)ry[i]);
      RCMz += m*la * ((double)rz[i]);
      MTOT += m;
    }
    
  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;

  printf("Box CoM RCMx: %.15f RCMy: %.15f RCMz: %.15f\n", RCMx, RCMy, RCMz);
}

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  
}


/* ======================== >>> structFacts <<< =============================*/
void structFacts(void)
{
  /* Riscrivere usando il codice del monoatomico!!!!!!!!!!!!!!!!!!!! */
}

/* =========================== >>> radDens <<< ============================= */
void radDens(void)
{
  /* Distribuzione radial per gli atomi di tipo A */
  int bin, Nm, i, j;
  int* hist;
  double rij;
  int rijSqI, rxijI, ryijI, rzijI;
  COORD_TYPE rlower, rupper, cost, nIdeal; 
  COORD_TYPE rhoAv, m, DELR;
  double la;
  //printf("calcolo gr!!!\n");
  L = cbrt(Vol);
  invL = 1.0 / L;
  Nm = Oparams.parnum;
  rhoAv = (COORD_TYPE) Nm;
  m = Oparams.m;
  la = Oparams.lattice_a;

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
	  rxijI = pmap[rxS[i] - rxS[j]];
	  ryijI = pmap[ryS[i] - ryS[j]];
	  rzijI = pmap[rzS[i] - rzS[j]];
	  rijSqI = Sqr(rxijI) + Sqr(ryijI) + Sqr(rzijI);
	  rij =  Oparams.lattice_a*sqrt(((double)rijSqI));
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
      rlower = RBEG + ((COORD_TYPE)bin) * DELR;
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



