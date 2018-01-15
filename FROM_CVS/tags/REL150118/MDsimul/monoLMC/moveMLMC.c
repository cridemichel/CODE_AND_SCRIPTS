#include <mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

extern double rcut, rcutSq, dvdr; 
extern double sigmaSq, epsilon4, epsilon24;

extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern void links(COORD_TYPE rcut, COORD_TYPE sigma);
extern void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigma);
extern void BuildNebrList(COORD_TYPE rCut, COORD_TYPE sigma);

extern void checkNebrRebuild(void);
extern void resetCM();
extern char TXT[MSG_LEN], TXTA[10][MSG_LEN];
extern inline int minimumImageI(int);


/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE K, Wxx, Wyy, Wzz, Wxy, Wyz, Wzx;  

COORD_TYPE DrSq = 0.0, Mtot;
/* used by linked list routines */

/* neighbour list method variables */
COORD_TYPE dispHi;
int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax, **nebrTaba, **nebrTabi;
/* ================================= */
double *cosLU, *sinLU, *sumCosIni, *sumSinIni, *sumCos, *sumSin, *DsumSin, *DsumCos;

extern int *rxS, *ryS, *rzS; 
extern int *kcx, *kcy, *kcz;
int mosseAccettate = 0;
const int spost[6][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1}};
/*const int spost[26][3] = 
{{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1},
 {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},
 {0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
 {1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{-1,-1,1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}
 };*/

/* ========================= >>> randDisp <<< ========================== */
inline void randDisp(int *dx, int *dy, int *dz)
{
  int rn;
  /* sceglie uno dei primi vicini per fare lo spostamento */
  /* Sceglie un numero casuale intero tra 0...5 */
  rn =(int) (6.0*rand()/(RAND_MAX+1.0));
  //rn =(int) (26.0*rand()/(RAND_MAX+1.0));
  *dx = spost[rn][0]; 
  *dy = spost[rn][1];
  *dz = spost[rn][2];
}


extern void calcEnergy(int i, int rxx, int ryy, int rzz,
		       int dx, int dy, int dz, 
		       double *dWLJ, double *dVLJ, double *dWAC, double *dVAC);

extern int NNC;
extern double ranf(void);
int nSTLen, NnSTLen, NOTnebrSkTabLen, nebrSkTabLen, nebrSkNow=1;
int *nebrSkTab, *NOTnebrSkTab, *NnST, *nST;


/* ========================== >>> updateSk << ============================= */
inline double calcSk(int ind)
{
  int rkio, i;
 
  sumCosIni[ind] = sumCos[ind];
  sumSinIni[ind] = sumSin[ind];
  
  sumCos[ind] = 0.0;
  sumSin[ind] = 0.0;
  /* ricalcola sumCos e sumSin per tutti i k */ 
  for (i= 0; i < Oparams.parnum; i++)
    {
      rkio = rxS[i] * kcx[ind] +  ryS[i] * kcy[ind] + rzS[i] * kcz[ind];
      sumCos[ind] += cosLU[rkio];
      sumSin[ind] += sinLU[rkio];
    }

  return Sqr(sumCos[ind]) + Sqr(sumSin[ind]);
}

 
/* ========================== >>> buildSkNebr <<< ======================== */
void buildSkNebr(void)
{
  int ind, nk, i, rkio;
  double DSmax, Skold, Sk, S0, DS0, sqrt2;
  S0 = OprogStatus.S0;
  DS0 = OprogStatus.DS0;//OprogStatus.DS0;


  DSmax = 0.0;

  nSTLen = 0;
  NnSTLen = 0;

  //  printf("S: %d %d+ %d= %d\n", Oparams.curStep, NOTnebrSkTabLen,nebrSkTabLen, NOTnebrSkTabLen+nebrSkTabLen);
  if (NOTnebrSkTabLen+nebrSkTabLen == 0)
    {
      nebrSkTabLen = NNC;
      /* se è la prima volta allora considera tutte le componenti ddi k */
      for (nk = 0; nk < nebrSkTabLen; nk++)
	{
	  nebrSkTab[nk] = nk;
	}
    }

  /* scorre la vecchia lista dei primi non-vicini per S(k) e verifica
     che nessuna componente non-vicina abbia "percorso" piu' di DS0 
     poiché in tale caso la scelta dei parametre DS0 e chiamate potrebbe essere errata
  */
  for (nk = 0; nk < NOTnebrSkTabLen; nk++)
    {
      ind = NOTnebrSkTab[nk];
      Sk = calcSk(ind);
      //printf("NOTnebrSkTabLen: %d Sk:%f\n", NOTnebrSkTabLen, Sk);
      if (Sk > (S0 - DS0))
	{
	  nST[nSTLen] = ind;
	  nSTLen++;
	}
      else
	{
	  Skold = Sqr(sumCosIni[ind]) + Sqr(sumSinIni[ind]);
	  if (Sk-Skold > DSmax)
	    DSmax = Sk-Skold; 
	  /* questa è la tabella con le componenti scartate */
	  NnST[NnSTLen] = ind;
	  NnSTLen++;
	  /* uso la tabella temporanea poiché altrimenti aggiornerei la tabella
	     che sto usando in questo ciclo */
	}
    }

  /* scorre la lista dei primi vicini per S(k) */
  for (nk = 0; nk < nebrSkTabLen; nk++)
    {
      ind = nebrSkTab[nk];
      Sk = calcSk(ind);
      //printf("nk: %d Sk:%f\n", nk,Sk);
      if (Sk > (S0 - DS0))
	{
	  nST[nSTLen] = ind;
	  nSTLen++;
	  /* come nel caso precedente uso la tabella temporanea poichè
	     la tabella nebrSkTab la sto scorrendo in questo ciclo */
	}
      else
	{
	  /* questa è la tabella con le componenti scartate(non-vicini) */
	  NnST[NnSTLen] = ind;
	  NnSTLen++;
	}
    }
  

  nebrSkTabLen = nSTLen;
  NOTnebrSkTabLen = NnSTLen;
 
  /* assegna le componenti di non-vicini, memorizzate temporaneamente in 
     NnST, all'array di non-vicini NOTnebrSkTab */
  for (nk = 0; nk < NOTnebrSkTabLen; nk++)
    {
      NOTnebrSkTab[nk] = NnST[nk];
    }

  /* assegna le componenti di vicini, memorizzate temporaneamente in nST,
     all'array di non-vicini nebrSkTab */
  for (nk = 0; nk < nebrSkTabLen; nk++)
    {
      nebrSkTab[nk] = nST[nk];
    }

  
  //printf("nebrLen:%d\n", nebrSkTabLen);
  /* se un non-vicino si è spostato per piu' di DS0 allora sono cazzi */
  if (DSmax > DS0 && Oparams.curStep > 1) 
    {
      /* Qui controlla che le componenti non presenti nella neighbour ist non siano 
	 variate per piu' di S0-DS0 infatti se cosi' fosse il calcolo del potenziale
	 anticristallino potrebbe essere errato */
      sprintf(TXT, "ATTENZIONE: DSmax:%f maggiore di DS0 = %f!!!!!!!\n", DSmax, DS0);
      mdPrintf(ALL, TXT, NULL);
    }

}

/* ========================== >>> buildSkNebr <<< ======================== */
void buildSkNebrFast(void)
{
  int ind;
  double Sk, S0, DS0;
  S0 = OprogStatus.S0;
  DS0 = OprogStatus.DS0;//OprogStatus.DS0;

  nebrSkTabLen = 0;

  for (ind = 0; ind < NNC; ind++)
    {
      Sk = calcSk(ind);
      //printf("NOTnebrSkTabLen: %d Sk:%f\n", NOTnebrSkTabLen, Sk);
      if (Sk > (S0 - DS0))
	{
	  nebrSkTab[nebrSkTabLen] = ind;
	  nebrSkTabLen++;
	}
    }
}


/* ==================== >>> chkSkNebr <<< ============================= */
void chkSkNebr(void)
{
  static int chiamate=0;

  if (chiamate > 0)
    {
      chiamate--;
      nebrSkNow = 0;
    }
  else
    {
      chiamate = OprogStatus.chiamate - 1;
      /* ogni chiamate volte rigenera la lista per Sk */
      nebrSkNow = 1;
      //printf("STEP: %d Rigenera la neighbour list per Sk\n", Oparams.curStep);
    }

}

/* ========================== >>> movea <<< =============================== */
inline void moveMC(COORD_TYPE m)
{  
  int i, nk;
  int dx, dy, dz, ind;
  double Beta, deltV, deltW, deltVB, deltVLJ, deltWLJ, deltVAC, deltWAC;
  int rxold, ryold, rzold;
 
  /* ===== LOOP OVER MOLECULES ===== */
  //L = cbrt(Vol);
  //invL = 1.0 / L;

  Beta = 1.0 / Oparams.T;
  
  chkSkNebr();

  if (nebrSkNow)
    buildSkNebr();

  /*
    for (ind = 0; ind < NNC; ind++)
    {
    sumCosIni[ind] = sumCos[ind];
    sumSinIni[ind] = sumCos[ind];
    }
  */
  //if (Oparams.curStep%100 == 0)printf("VAC:%f\n", VAC);
  for(i = 0; i < Oparams.parnum; i++)
    {
      rxold = rxS[i];
      ryold = ryS[i];
      rzold = rzS[i];
      randDisp(&dx, &dy, &dz);
      
      calcEnergy(i, rxold, ryold, rzold, dx, dy, dz,
		 &deltWLJ, &deltVLJ, &deltWAC, &deltVAC);
      //deltVB = Beta * deltV;
      deltV = deltVLJ + deltVAC;

      //printf("deltVLJ:%f deltVAC: %f deltV: %f\n", deltVLJ, deltVAC, deltV);
      //printf("nebrSkTabLen: %d\n", nebrSkTabLen);
#ifdef PRESSURE
      deltW = deltWLJ + deltWAC;
#endif      
      deltVB = Beta * deltV;
      OprogStatus.mosseTentate += 1;
      /* NOTA: Vc (mentre V è quello non shiftato)
	 è il potenziale schiftato e qui uso quello
	 per stabilire se la mossa è eccettabile o meno */
      
      if ( deltVB < 75.0 ) 
	    {
	      if ( deltV <= 0.0 )  
                {
		  V += deltV;
		  VAC += deltVAC;
		  VLJ += deltVLJ;
#ifdef PRESSURE
		  W += deltW;
		  WAC += deltWAC;
		  WLJ += deltVLJ;
#endif
		  for (nk = 0; nk < nebrSkTabLen; nk++)
		    {
		      ind = nebrSkTab[nk];
		      sumCos[ind] += DsumCos[ind];
		      sumSin[ind] += DsumSin[ind];
		    }

		  rx[i] += dx;
		  ry[i] += dy;
		  rz[i] += dz;
		  rxS[i] = minimumImageFB(rxold + dx);
		  ryS[i] = minimumImageFB(ryold + dy);
		  rzS[i] = minimumImageFB(rzold + dz);

		
		  OprogStatus.mosseAccettate += 1;
		}
	      else if ( exp(-deltVB) > ranf() ) 
		{
		  V  += deltV;
		  VAC += deltVAC;
		  VLJ += deltVLJ;
		  for (nk = 0; nk < nebrSkTabLen; nk++)
		    {
		      ind = nebrSkTab[nk];
		      sumCos[ind] += DsumCos[ind];
		      sumSin[ind] += DsumSin[ind];
		    }

#ifdef PRESSURE
		  W  += deltW;
		  WAC += deltWAC;
		  WLJ += deltVLJ;

#endif
		  rx[i] += dx;
		  ry[i] += dy;
		  rz[i] += dz;
		  rxS[i] = minimumImageFB(rxold + dx);
		  ryS[i] = minimumImageFB(ryold + dy);
		  rzS[i] = minimumImageFB(rzold + dz);

		  OprogStatus.mosseAccettate += 1;
		}
	    }
    }
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
void scalCor()
{ 
  int i;
  COORD_TYPE cost, costo1, costo2, DRx, DRy, DRz;
  
  L = cbrt(Vol);
  invL = 1.0 / L;
  cost = Vol1 / Vol / 3.0;
  costo1 = Vol1o1 / Vol / 3.0;
  costo2 = Vol1o2 / Vol / 3.0;
  
  for(i = 0; i < Oparams.parnum; i++)
    {
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * rx[i]);
      DRy = - L * rint(invL * ry[i]);
      DRz = - L * rint(invL * rz[i]);
      
      rx[i] += DRx;
      ry[i] += DRy;
      rz[i] += DRz;
      rxS[i] = minimumImageI(rx[i]);
      ryS[i] = minimumImageI(ry[i]);
      rzS[i] = minimumImageI(rz[i]);
    }
}


/* ============================ >>> chkVol <<< ============================= */
void chksVol(void)
{
  /* DESCRIPTION:
     This subroutine check if the actual volume is near the volume
     OprogStatus.avVol, if it is the case the program exit 
     This is useful if we want to perform a micronical production run 
     after the equilibration one (NPT) */
  COORD_TYPE DVol, relVol;
  COORD_TYPE tolVol;
  
  tolVol  = OprogStatus.tolVol;

  if (OprogStatus.avVol > 0.0)
    {
      DVol = fabs(Vol - OprogStatus.avVol);
      relVol = DVol / Vol;
      
      if (relVol < tolVol)
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

/* =========================== >>> calcConst <<< ========================= */
void calcConst(void)
{
  
  double sr2, sr6, sr12;

  /* useful ab-constants inside OUTER LOOP below */
  rcut = Oparams.rcut * Oparams.sigma;
  
  if ((rcut + Oparams.lattice_a) > (cbrt(Vol)/2))  
    {
      printf("VALORI ERRATI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
      printf("il raggio di cutoff deve essere inferiore\n");
      printf("a mezza scatola di un passo reticolare ossia:%.10f",
	     Oparams.lattice_a);
    }

  rcutSq = Sqr(rcut);
  sigmaSq = Sqr(Oparams.sigma);
  epsilon4 = 4.0 * Oparams.epsilon;
  epsilon24 = 24.0 * Oparams.epsilon;
  sr2 = sigmaSq / rcutSq;
  sr6 = sr2 * sr2 * sr2;
  sr12 = Sqr(sr6);
  dvdr = epsilon24 * (sr6 - 2.0 * sr12) / rcut;
  /* initialize ab-variables */

  //V = 0.0; /* potential energy */
  //W = 0.0; /* virial function */
  //Vc = 0.0;
  //L = cbrt(Vol);
  //invL = 1.0  / L;

}

/* =========================== >>> adjQty <<< ========================= */


/* ============================ >>> move<<< =================================*/
inline void move(void)
{
  double moveRate;

  //printf("Vol: %f\n", Vol);
  if (nebrNow)
    {
      //printf("building neighbour listr la:%.5f step: %d\n", 
      //     Oparams.lattice_a, Oparams.curStep);
      nebrNow = 0;
      dispHi = 0.0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      BuildNebrListNoLinked(Oparams.rcut, Oparams.sigma);
    }

  calcConst();
 
  moveMC(Oparams.m);
  
  checkNebrRebuild();
  
  if ((OprogStatus.rateCheck != 0) && 
      (Oparams.curStep % OprogStatus.rateCheck == 0))
    {
      moveRate = ((double) OprogStatus.mosseAccettate) / 
	((double)OprogStatus.mosseTentate);
      if (moveRate > 0.5)
	{
	  sprintf(TXTA[0], "Il Rate è alto\n");
	}
      else
	{
	  sprintf(TXTA[0], "Il Rate è basso\n");
	}

      sprintf(TXTA[1], "acc: %d tent: %d moveRate: %.5f\n", 
	      OprogStatus.mosseAccettate,OprogStatus.mosseTentate, moveRate);
      mdPrintf(ALL, TXTA[0], TXTA[1], NULL);
      OprogStatus.mosseTentate = OprogStatus.mosseAccettate = 0;
      //srand((int)time(NULL));
    }
  
  if (  ( (OprogStatus.CMreset > 0) &&
	  // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	  (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	     // 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) 
	     (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
    resetCM(Oparams.parnum);
  
}

