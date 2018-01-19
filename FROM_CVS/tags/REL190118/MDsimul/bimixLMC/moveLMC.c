#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

extern double rcutab[NA][NA], rcutabSq[NA][NA], dvdr[NA][NA]; 
extern double sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA];

extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern void links(COORD_TYPE rcut, COORD_TYPE sigab[NA][NA]);
extern void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]);
extern void BuildNebrList(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]);

extern void checkNebrRebuild(void);
extern void resetCM();
extern char TXT[MSG_LEN], TXTA[10][MSG_LEN];
extern inline int minimumImageI(int);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE W, K, Wxx, Wyy, Wzz, Wxy, Wyz, Wzx;  

COORD_TYPE DrSq = 0.0, Mtot;
/* used by linked list routines */

/* neighbour list method variables */
COORD_TYPE dispHi;
int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax, **nebrTaba, **nebrTabi;
/* ================================= */

extern int *rxS[NA], *ryS[NA], *rzS[NA]; 
int mosseAccettate = 0;
// const int spost[6][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1}};
const int spost[26][3] = 
{{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1},
 {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},
 {0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
 {1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{-1,-1,1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}
 };
 
/* ========================= >>> randDisp <<< ========================== */
inline void randDisp(int *dx, int *dy, int *dz)
{
  int rn;
  /* sceglie uno dei primi vicini per fare lo spostamento */
  /* Sceglie un numero casuale intero tra 0...5 */
  //rn =(int) (6.0*rand()/(RAND_MAX+1.0));
  rn =(int) (26.0*rand()/(RAND_MAX+1.0));
  *dx = spost[rn][0]; 
  *dy = spost[rn][1];
  *dz = spost[rn][2];
}


extern void calcEnergy(int a, int i, int rxx, int ryy, int rzz,
		       int dx, int dy, int dz, 
		       double *dW, double *dV);

extern double ranf(void);

/* ========================== >>> movea <<< =============================== */
inline void moveMC(COORD_TYPE m[NA])
{  
  int i, a;
  int dx, dy, dz;
  double Beta, deltV, deltW, deltVB;
  int rxnew, rynew, rznew, rxold, ryold, rzold;
 
  /* ===== LOOP OVER MOLECULES ===== */
  //L = cbrt(Vol);
  //invL = 1.0 / L;

  Beta = 1.0 / Oparams.T;

  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
#ifdef LMC_FAST
	  rxold = rxS[a][i];
	  ryold = ryS[a][i];
	  rzold = rzS[a][i];
	  randDisp(&dx, &dy, &dz);
	  calcEnergy(a, i, rxold, ryold, rzold, dx, dy, dz, 
		     &deltW, &deltV);
	  //printf("deltVc_ %.10f\n", deltVc);
	  deltVB = Beta * deltV;
	  OprogStatus.mosseTentate += 1;
	  if ( deltVB < 75.0 ) 
	    {
	      if ( deltV <= 0.0 )  
                {
		  //Vc += deltVc;    
		
		  V += deltV;
#ifdef PRESSURE
		  W += deltW;
#endif
		  rx[a][i] += dx;
		  ry[a][i] += dy;
		  rz[a][i] += dz;
		  
		  rxS[a][i] = minimumImageFB(rxold + dx);
		  ryS[a][i] = minimumImageFB(ryold + dy);
		  rzS[a][i] = minimumImageFB(rzold + dz);
		
		  OprogStatus.mosseAccettate += 1;
		}
	      else if ( exp(-deltVcB) > ranf() ) 
		{
		  //Vc += deltVc;
		  V  += deltV;
#ifdef PRESSURE
		  W  += deltW;
#endif
		  rx[a][i] += dx;
		  ry[a][i] += dy;
		  rz[a][i] += dz;
		  rxS[a][i] = minimumImageFB(rxold + dx);
		  ryS[a][i] = minimumImageFB(ryold + dy);
		  rzS[a][i] = minimumImageFB(rzold + dz);

		  OprogStatus.mosseAccettate += 1;
		}
	    }
#else
	  rxold = rxS[a][i];
	  ryold = ryS[a][i];
	  rzold = rzS[a][i];
	  randDisp(&dx, &dy, &dz);
	  calcEnergy(a, i, rxold, ryold, rzold, dx, dy, dz,
		     &deltW, &deltV);
	  //deltVB = Beta * deltV;
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
#ifdef PRESSURE
		  W += deltW;
#endif
		  rx[a][i] += dx;
		  ry[a][i] += dy;
		  rz[a][i] += dz;
		  rxS[a][i] = minimumImageFB(rxold + dx);
		  ryS[a][i] = minimumImageFB(ryold + dy);
		  rzS[a][i] = minimumImageFB(rzold + dz);

		
		  OprogStatus.mosseAccettate += 1;
		}
	      else if ( exp(-deltVB) > ranf() ) 
		{
		  V  += deltV;
#ifdef PRESSURE
		  W  += deltW;
#endif
		  rx[a][i] += dx;
		  ry[a][i] += dy;
		  rz[a][i] += dz;
		  rxS[a][i] = minimumImageFB(rxold + dx);
		  ryS[a][i] = minimumImageFB(ryold + dy);
		  rzS[a][i] = minimumImageFB(rzold + dz);

		  OprogStatus.mosseAccettate += 1;
		}
	    }
#endif	      

	}
    }
  
  return;
for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (rxS[a][i] != rx[a][i] || ryS[a][i] != ry[a][i] || 
	      rzS[a][i] != rz[a][i])
	    {
	      printf("step: %d PORCA TROIA!!!!!!!!!!!!!!!!!\n", Oparams.curStep);
	      printf("[%d][%d] Scaled(%d,%d,%d) (%d,%d,%d)\n", a, i,
		     rxS[a][i], ryS[a][i], rzS[a][i], rx[a][i], 
		     ry[a][i], rz[a][i]);
	      //exit(-1);
	      
	      printf("MI[(%d,%d,%d)]=(%d,%d,%d)\n", 
		     rx[a][i], ry[a][i], rz[a][i], minimumImageI(rx[a][i]),
		     minimumImageI(ry[a][i]), minimumImageI(rz[a][i]));
	      exit(-1);
	    }
	}
    }

    exit(-1);
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
  int i, a;
  COORD_TYPE cost, costo1, costo2, DRx, DRy, DRz;
  
  L = cbrt(Vol);
  invL = 1.0 / L;
  cost = Vol1 / Vol / 3.0;
  costo1 = Vol1o1 / Vol / 3.0;
  costo2 = Vol1o2 / Vol / 3.0;
  
  /* Reduced particles to first box */
  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	     scale them */
	  DRx = - L * rint(invL * rx[a][i]);
	  DRy = - L * rint(invL * ry[a][i]);
	  DRz = - L * rint(invL * rz[a][i]);

    	  rx[a][i] += DRx;
	  ry[a][i] += DRy;
	  rz[a][i] += DRz;
	  rxS[a][i] = minimumImageI(rx[a][i]);
	  ryS[a][i] = minimumImageI(ry[a][i]);
	  rzS[a][i] = minimumImageI(rz[a][i]);
	}
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
  int a, b;
  double srab2, srab6, srab12;

  for(a = 0; a < NA; a++)
    {
      for(b = 0; b < NA; b++) /* b >= a because of symmetry */
	{
	  /* useful ab-constants inside OUTER LOOP below */
	  rcutab[a][b] = Oparams.rcut * Oparams.sigab[a][b];
	  if ((rcutab[a][b]+Oparams.lattice_a) > (cbrt(Vol)/2))  
	    {
	      printf("VALORI ERRATI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
	      printf("a=%d b=%d il raggio di cutoff deve essere inferiore\n",a,b);
	      printf("a mezza scatola di un passo reticolare ossia:%.10f",
		     Oparams.lattice_a);
	    }
	  rcutabSq[a][b] = Sqr(rcutab[a][b]);
	  sigabSq[a][b] = Sqr(Oparams.sigab[a][b]);
	  epsab4[a][b] = 4.0 * Oparams.epsab[a][b];
	  epsab24[a][b] = 24.0 * Oparams.epsab[a][b];
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = Sqr(srab6);
	  dvdr[a][b] = epsab24[a][b] * (srab6 - 2.0 * srab12) / rcutab[a][b];
	  /* initialize ab-variables */
	}
    }

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
 
  if (nebrNow)
    {
      //printf("building neighbour listr la:%.5f step: %d\n", 
      //     Oparams.lattice_a, Oparams.curStep);
      nebrNow = 0;
      dispHi = 0.0;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      BuildNebrListNoLinked(Oparams.rcut, Oparams.sigab);
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

