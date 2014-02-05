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

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, Elrc, Plrc;   

/* used by linked list routines */

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax, **nebrTaba, **nebrTabi;
extern int nebrSkTabLen, *nebrSkTab;

extern int cubeSize;
extern int *pmap, *pabs, *pSqr, *pgetj, *pgetb;
extern double *wRadLT, *vRadLT;

extern double *sumCos, *sumSin, *DsumSin, *DsumCos;
int *rxS, *ryS, *rzS;
int maxdelta;

/* ================================= */


double rcut, rcutSq, dvdr; 
double sigmaSq, epsilon4, epsilon24;
/* ================== >>> NO LINKED LIST IN MC <<< ===================== */

inline int minimumImageI(int pl);

/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigma) 
{
  int tot=0, i, j;
  COORD_TYPE rcut, rcutSq; 
  COORD_TYPE rrNebr;
  int rrNebrI;
  double la, rrNebrD;
  int Nm, rxi, ryi, rzi, rxijI, ryijI, rzijI, rijSqI;

  L = cbrt(Vol);
  invL = 1.0  / L;
  Nm = Oparams.parnum;
  la = Oparams.lattice_a;

  
  /* useful ab-constants inside OUTER LOOP below */
  rcut = rCut * sigma;
  rcutSq = Sqr(rcut);
  rrNebr = Sqr(rcut + OprogStatus.rNebrShell);
  rrNebrD = Sqr((rcut + OprogStatus.rNebrShell)/
		Oparams.lattice_a);
  
  rrNebrI = (int) rrNebrD;
  if (((double)rrNebrI) < rrNebrD)
    rrNebrI += 1;
  
  if (rrNebr > Sqr(L / 2.0))
    {
      printf("(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
	     sqrt(rrNebr), L/2.0);
      exit(-1);
    }
  //printf("sqrt(rrNebr[%d][%d]):%f\n", a, b, sqrt(rrNebr[a][b]));

  for (i = 0; i < Oparams.parnum; i++ )
    {
      rxi = rxS[i];
      ryi = ryS[i];
      rzi = rzS[i];
      
      nebrTabLen[i] = 0;
      for (j = 0; j < Oparams.parnum; j++)
	{
	  /* Bisogna assicurare solo il fatto che ogni particella
	     non puo' interagire con se stessa in un Monte Carlo */
	  if ( j == i )
	    continue;
	  /* INNER LOOP BEGINS */
	  rxijI = pmap[rxi - rxS[j]];
	  ryijI = pmap[ryi - ryS[j]];
	  rzijI = pmap[rzi - rzS[j]];
		  
	  /*
	    rxab = rxab - L * rint(rxab * invL);  
	    ryab = ryab - L * rint(ryab * invL);
	    rzab = rzab - L * rint(rzab * invL);*/
	  rijSqI = pSqr[rxijI] + pSqr[ryijI] + pSqr[rzijI];

	  if ( rijSqI < rrNebrI )/* 'rcut' is the cutoff for V */
	    {
	      nebrTab[i][nebrTabLen[i]] = j;
	      ++nebrTabLen[i];
	      //printf("nebrTabLen[%d]:%d\n", i, nebrTabLen[i]);
	    }
	}
      tot += nebrTabLen[i];

      /* INNER LOOP ENDS */
    }
  //printf("STEP:%d tot: %d\n", Oparams.curStep, tot);
}

/* ===================== >>> checkNebrRebuild <<< ======================== */
void checkNebrRebuild(void)
{
  /* Lo spostamento massimo lo conosco a priori poiché 
     so che ogni particella puo' andare in un sito primo 
     vicino e il passo reticolare è noto */
   
  dispHi = dispHi + Oparams.lattice_a;
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */

  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;

}

/* ========================== >>> minimumImageFB <<< =======================*/
inline int minimumImageFB(int pl)
{
  /* Calcola l'immagine minima assumendo che l'atomo si trovi
     nella prima scatola (quella che include l'origine */
  //  printf("PRIMA: %d\n", *pv);
  if (abs(pl) > (Oparams.lattice_M / 2))
    {
      if (pl > 0)
	pl -= Oparams.lattice_M;
      else
	pl += Oparams.lattice_M;
    }
  return pl;

}
/* ======================== >>> minimumImage <<< =========================== */
inline int minimumImageI(int pl)
{
  /* QUESTA E' LENTA DA FAR SCHIFO!!!!!!!!!!!!!!!! */
  /* calcola l'immagine minima di A cioè: 
     A - lM * rint(lM*A) non usando però rint ma solo
     funzioni intere per avere una maggiore efficienza */
  int lM, halfLM;

  lM = Oparams.lattice_M;
  halfLM = lM / 2;
  pl -= lM * (pl / lM);
  if (abs(pl % lM) > halfLM)
    {
      if (pl > 0)
	pl -= lM;
      else
	pl += lM;
    }
  return pl;
}

extern int csSq, NNC, KK;
extern int *kcx, *kcy, *kcz;
extern double *cosLU, *sinLU;
extern double **dcosLU, **dsinLU;
/* ============================ >>> energy <<< ==============================*/
void calcEnergy(int i, int rxx, int ryy, int rzz, 
		int dx, int dy, int dz,
		double *dWLJ, double *dVLJ, double *dWAC, double *dVAC)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  register int j;
  double ddSq, VACN;
  double sDs, sSs,dd, la;
  int drki, rki, rkio, rxi, ryi, rzi, rxijI, ryijI, rzijI;
  /* Local variables to implement linked list */
  register int  Nm, n, nk;
  int ind, DrijSqI, rijSqI, lM, numk;// csSq[NA][NA];
  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */
  
  Nm = Oparams.parnum;
  /* Loop over cells */
  
  rxi = rxx;
  ryi = ryy;
  rzi = rzz;

  la = Oparams.lattice_a;
  lM = Oparams.lattice_M;
  
  *dVLJ = 0.0;
  *dWLJ = 0.0;
  numk=0;
  //printf("nebrTabLen[%d]:%d\n", i, nebrTabLen[i]);
  for(n = 0; n < nebrTabLen[i]; n++)
    { 
     
      j = nebrTab[i][n]; 
      
      rxijI = pmap[rxi - rxS[j]];
      ryijI = pmap[ryi - ryS[j]];
      rzijI = pmap[rzi - rzS[j]];

      rijSqI = pSqr[rxijI] + pSqr[ryijI] + pSqr[rzijI];
      //DrabSqI =0;
      if (dx)
	{
	  DrijSqI = 1 + 2*dx*rxijI;
	}
      else if (dy)
	{
	  DrijSqI = 1 + 2*dy*ryijI;
	}
      else //if (dz) 
	{
	  DrijSqI = 1 + 2*dz*rzijI;
	}
      
      if (rijSqI < csSq)
	{
	  /* store old values */
	  *dVLJ -= vRadLT[rijSqI];
#ifdef PRESSURE
	  *dWLJ -= wRadLT[rijSqI]; 
#endif
	}

      rijSqI += DrijSqI;
      if  (rijSqI < csSq)
	{
	  *dVLJ += vRadLT[rijSqI];
#ifdef PRESSURE
	  *dWLJ += wRadLT[rijSqI]; 
#endif
	}
    }
  
  //  return;
  /* ======== E ORA IL POTENZIALE ANTI-CRISTALLINO ======= */
  *dVAC = 0.0;
  *dWAC = 0.0;
  VACN = 0.0;
  /* calcola le variazioni degli S(k) per tutti i k...*/
  for(nk = 0; nk < nebrSkTabLen; nk++)
    {
      ind = nebrSkTab[nk];
      
      //printf("NNC:%d\n", NNC);
      //printf("ind:%d nebrSkTabLen:%d\n", ind, nebrSkTabLen);
      /* il termine cosLU[rki] che corrisponde alla poiz. prec.
	 potrebbe essere memorizzato in un array (ma ci vuole tanta memoria) 
	 NOTA: Ora non calcola il viriale AC */
      rkio = rxx * kcx[ind] +  ryy * kcy[ind] + rzz * kcz[ind];
      //DsumCos[ind] = -cosLU[rki];
      //DsumSin[ind] = -sinLU[rki];
      if (dx)
	{
	  drki = dx*kcx[ind] + KK;
	}
      else if (dy)
	{
	  drki = dy*kcy[ind] + KK;
	}
      else //if(dz) 
	{
	  drki = dz*kcz[ind] + KK;
	}
      /*
	DsumCos[ind] = cosLU[rki] - cosLU[rkio];
	DsumSin[ind] = sinLU[rki] - sinLU[rkio];
      */

      /* sumCos[ind] è la somma sulle particelle di cos(kr) osssia
	 sumCos(k) = Sum_i cos(k r_i) che dipende ovviamente da k */ 
      DsumCos[ind] = dcosLU[drki][rkio];
      DsumSin[ind] = dsinLU[drki][rkio];
      
      /* ...ma controlla se è il caso di aggiustare il potenziale anticristallino
	 usando la nighbour list per Sk */
      sDs = sumCos[ind] + DsumCos[ind];
      sSs = sumSin[ind] + DsumSin[ind];
      //printf(": %f sSs:%f\n", sDs, sSs);

      dd = Sqr(sDs) + Sqr(sSs); //<------------ il peso è tutto qui!!!!!!!!!!
      if (dd > OprogStatus.S0) 
	{
	  numk++;
	  dd -= OprogStatus.S0;
	  //printf("step: %d VAC[%d]: %f\n", Oparams.curStep, ind, 
	  // 0.5*OprogStatus.alpha*Sqr(dd));
	  ddSq = Sqr(dd);
	  VACN += ddSq;
	  //VACk[ind] = ddSq;
	}
    }
  
  *dVAC = OprogStatus.alpha * VACN - VAC;
  //printf("i:%d numk:%d\n", i, numk);
#ifdef PRESSURE
#else
  *dWAC = 0.0;
#endif
}
  

