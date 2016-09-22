#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

#define LOOKUP
#define RADE
extern int ENDSIM;
extern char msgStrA[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W;  

/* used by linked list routines */

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax, **nebrTaba, **nebrTabi;

extern double ***vabLT[NA][NA];
extern double ***wabLT[NA][NA];
extern int cubeSize[NA][NA];
extern int *pmap, *pabs, *pSqr, *pgetj, *pgetb;
extern double *wabRadLT[NA][NA], *vabRadLT[NA][NA];
extern double **DwabRadLT[NA][NA], **DvabRadLT[NA][NA];

int *rxS[NA], *ryS[NA], *rzS[NA];
int maxdelta[NA][NA];

/* ================================= */


double rcutab[NA][NA], rcutabSq[NA][NA], dvdr[NA][NA]; 
double sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA];
/* ================== >>> NO LINKED LIST IN MC <<< ===================== */

inline int minimumImageI(int pl);

/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]) 
{
  int aiat, bjat, a, b, i, j;
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE rrNebr[NA][NA];
  int rrNebrI[NA][NA];
  COORD_TYPE rabSq, rxab, ryab, rzab;
  double la, rrNebrD;
  int Nm, rxa, rya, rza, rxabI, ryabI, rzabI, rabSqI;

  L = cbrt(Vol);
  invL = 1.0  / L;
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  la = Oparams.lattice_a;

  for(a = 0; a < NA; a++)
    {
      for(b = 0; b < NA; b++)
	{
	  /* useful ab-constants inside OUTER LOOP below */
	  rcutab[a][b] = rCut * sigab[a][b];
	  rcutabSq[a][b] = Sqr(rcutab[a][b]);
	  rrNebr[a][b] = Sqr(rcutab[a][b] + OprogStatus.rNebrShell);
	  rrNebrD = Sqr((rcutab[a][b] + OprogStatus.rNebrShell)/
	    Oparams.lattice_a);
	  
	  rrNebrI[a][b] = (int) rrNebrD;
	  if (((double)rrNebrI[a][b]) < rrNebrD)
	    rrNebrI[a][b] += 1;
	  
	  if (rrNebr[a][b] > Sqr(L / 2.0))
	    {
	      printf("(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
		     sqrt(rrNebr[a][b]), L/2.0);
	      exit(-1);
	    }
	  //printf("sqrt(rrNebr[%d][%d]):%f\n", a, b, sqrt(rrNebr[a][b]));
	}
    }

  
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++ )
      //loop(i, 1, Nm - 1)
	{
	  rxa = rxS[a][i];
	  rya = ryS[a][i];
	  rza = rzS[a][i];
	  aiat = i+Nm*a;
	  nebrTabLen[aiat] = 0;
	  for (b = 0; b < NA; b++)
	    {
	      for (j = 0; j < Oparams.parnum[b]; j++)
		//loop(b, 1, NA)  /* b >= a because of symmetry */
		{
		  /* Bisogna assicurare solo il fatto che ogni particella
		     non puo' interagire con se stessa in un Monte Carlo */
		  bjat = Nm*b+j;
		  if ( (bjat) == (aiat) )
		    continue;
		  /* INNER LOOP BEGINS */
		  //loop(j, i + 2, Nm) 
		  /* + 2 because you must remeber that really the all indices 
		     (a, b, i, j) start from 0 ... */
		  /*   i > j because of 3rd law */
		  
		  rxabI = pmap[rxa - rxS[b][j]];
		  ryabI = pmap[rya - ryS[b][j]];
		  rzabI = pmap[rza - rzS[b][j]];
		  
		  /*
		    rxab = rxab - L * rint(rxab * invL);  
		    ryab = ryab - L * rint(ryab * invL);
		    rzab = rzab - L * rint(rzab * invL);*/
		  rabSqI = pSqr[rxabI] + pSqr[ryabI] + pSqr[rzabI];
		  /*
		    printf("la:%.6f rabSq:%.6f rrNebr[%d][%d]:%.5f\n", la, rabSq, a, b,
		    rrNebr[a][b]);
		  */
		  if ( rabSqI < rrNebrI[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      /*
			if (nebrTabLen[aiat] >= nebrTabMax)
			{
			printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			nebrTabLen[i+a*Nm]);
			   printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
			   printf("ERROR: Neighbourlist overflow!\n");
			   exit(-1);
			   }
		      */
		      nebrTab[aiat][nebrTabLen[aiat]] = bjat;
		      //nebrTaba[aiat][nebrTabLen[aiat]] = a;
		      //nebrTabi[aiat][nebrTabLen[aiat]] = i;
		      //printf("(%d-%d)[(%d,%d)-(%d,%d)]\n ",
		      //     i + a*Nm, j + b*Nm, a, i, b, j);
		      
		      /* nebrTabLen[i+a*Nm] contiene il numero di primi
			 vicini dell'atomo (a,i) */
		      ++nebrTabLen[aiat];
		      /* Increment table element counter */
		      
		    }
		}
	      /* INNER LOOP ENDS */
	    }
	}
    }
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

extern int csSq[NA][NA];

#ifdef LMC_FAST
/* ============================ >>> energy <<< ==============================*/
void calcEnergy(int a, int i, int rxx, int ryy, int rzz, 
		int dx, int dy, int dz,
		double *dV, double *dW, double *dVc)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  register int b, j;
  double rxab, ryab, rzab, rabSq;
  double la, vab, vabCut, wab, srab2, srab6, srab12;
  int arxabI, aryabI, arzabI, rxa, rya, rza, rxabI, ryabI, rzabI;
  /* Local variables to implement linked list */
  register int  Nm, n, atia, nebrat;
  int rabSqI, cs, lM, aa, bb;//, csSq[NA][NA];
  int DrabSqI;
  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */
  
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  /* Loop over cells */
  
  atia = i+a*Nm;

  rxa = rxx;
  rya = ryy;
  rza = rzz;

  la = Oparams.lattice_a;
  lM = Oparams.lattice_M;


  *dV = 0.0;
  *dW = 0.0;
  *dVc = 0.0;
  
  //  return;
  for (aa = 0; aa < NA; aa++)
    {
      for (bb = 0; bb < NA; bb++)
	{
	  csSq[aa][bb] = Sqr(cubeSize[aa][bb]);
	}
    }
  
  for(n = 0; n < nebrTabLen[atia]; n++)
    { 
      nebrat = nebrTab[atia][n];

      if (nebrat >= Nm)
	{
	  b = 1;
	  j = nebrat - Nm;
	}
      else
	{
	  b = 0;
	  j = nebrat;
	}

      rxabI = pmap[rxa - rxS[b][j]];
      ryabI = pmap[rya - ryS[b][j]];
      rzabI = pmap[rza - rzS[b][j]];
      
      rabSqI = pSqr[rxabI] + pSqr[ryabI] + pSqr[rzabI];

      if (dx)
	{
	  DrabSqI = 1 + 2*dx*rxabI;
	}
      //else if (dy)
      if (dy)
	{
	  DrabSqI = 1 + 2*dy*ryabI;
	}
      //else
      if(dz) 
	{
	  DrabSqI = 1 + 2*dz*rzabI;
	}
      // printf("dx:%d dy:%d dz:%d DrabSqI: %d (%d,%d,%d)[%d]\n", dx, dy, dz, 
      //    DrabSqI,
      //    rxabI, ryabI, rzabI, maxdelta[a][b]);
      if (rabSqI < csSq[a][b])
	{
	  *dVc += DvabRadLT[a][b][DrabSqI+maxdelta[a][b]][rabSqI];
	  //if (DrabSqI == 0)
	  //deltWab[a][b] += wabRadLT[a][b][rabSqI]; 
	}
    }
  /* OUTER LOOP ENDS */
}  
#else
/* ============================ >>> energy <<< ==============================*/
void calcEnergy(int a, int i, int rxx, int ryy, int rzz, 
		int dx, int dy, int dz,
		double *dW, double *dV)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  register int b, j;
  double rxab, ryab, rzab, rabSq;
  double la, vab, vabCut, wab, srab2, srab6, srab12;
  int arxabI, aryabI, arzabI, rxa, rya, rza, rxabI, ryabI, rzabI;
  /* Local variables to implement linked list */
  register int  Nm, n, atia, nebrat;
  int DrabSqI, rabSqI, cs, lM, aa, bb;// csSq[NA][NA];

  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */
  
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  /* Loop over cells */
  
  atia = i+a*Nm;

  rxa = rxx;
  rya = ryy;
  rza = rzz;

  la = Oparams.lattice_a;
  lM = Oparams.lattice_M;
  
  *dV = 0.0;
  *dW = 0.0;

  /*
    for (aa = 0; aa < NA; aa++)
    {
    for (bb = 0; bb < NA; bb++)
    {
    csSq[aa][bb] = Sqr(cubeSize[aa][bb]);
    }
    
    }
  */
  for(n = 0; n < nebrTabLen[atia]; n++)
    { 
      nebrat = nebrTab[atia][n];

      if (nebrat >= Nm)
	{
	  b = 1;
	  j = nebrat - Nm;
	}
      else
	{
	  b = 0;
	  j = nebrat;
	}


      rxabI = pmap[rxa - rxS[b][j]];
      ryabI = pmap[rya - ryS[b][j]];
      rzabI = pmap[rza - rzS[b][j]];

#ifdef LOOKUP
      rabSqI = pSqr[rxabI] + pSqr[ryabI] + pSqr[rzabI];
      DrabSqI =0;
      if (dx)
	{
	  DrabSqI += 1 + 2*dx*rxabI;
	}
      //else if (dy)
      if (dy)
	{
	  DrabSqI += 1 + 2*dy*ryabI;
	}
      //else 
      if (dz) 
	{
	  DrabSqI += 1 + 2*dz*rzabI;
	}

      
      if (rabSqI < csSq[a][b])
	{
	  /* store old values */
	  *dV -= vabRadLT[a][b][rabSqI];
#ifdef PRESSURE
	  *dW -= wabRadLT[a][b][rabSqI]; 
#endif
	}

      rabSqI += DrabSqI;
      if  (rabSqI < csSq[a][b])
	{
	  *dV += vabRadLT[a][b][rabSqI];
#ifdef PRESSURE
	  *dW += wabRadLT[a][b][rabSqI]; 
#endif
	}
#else
      /* Moltiplica le coordinate intere per il passo del reticolo
	 in modo da ottenere le coordinate in unità ridotte */

      rxab = la * ((double) rxabI);
      ryab = la * ((double) ryabI);
      rzab = la * ((double) rzabI);
      
      /*
	rxab = rxab - L * rint(invL * rxab);     
	ryab = ryab - L * rint(invL * ryab);
	rzab = rzab - L * rint(invL * rzab);
      */
      
      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
	     
      if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
	{
	  srab2   = sigabSq[a][b] / rabSq;
	  srab6   = srab2 * srab2 * srab2;
          srab12  = Sqr(srab6);
	  
	  vab     = srab12 - srab6;
	  //vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);
	  wab     = vab + srab12;
	
	  Vab[a][b] = Vab[a][b] + vab;
	  //printf("(%d,%d)-(%d,%d):%.5f | ", a, i, b, j, vab);
	  /* total potential between all a-b atoms pairs */
	  Wab[a][b]   = Wab[a][b] + wab; 
	  /* Virial off-diagonal terms of atomic pressure tensor */
          ++ncut[a][b];
	}
#endif
    }

    /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
 /* OUTER LOOP ENDS */
}  
#endif
