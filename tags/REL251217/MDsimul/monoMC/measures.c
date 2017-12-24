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
extern COORD_TYPE pi, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, WC;

/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi, ratio;
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
  COORD_TYPE RCMx, RCMy;
  int mol, Nm, i;
  COORD_TYPE L, invL;

  Nm = Oparams.parnum;
  
  printf("STEP %d current acceptance ratio: %.15g DrMax:%.8g Taux: %.15G\n", 
	 Oparams.curStep,
	 ratio, OprogStatus.drMax,Oparams.Taux);
  Vtot =  Vlj+Vg+Vf;
  sprintf(TXTA[0], "  V=%.15f Vlj=%.15e Vg=%.15e Vf=%.15e Vw=%15e\n", Vtot, Vlj,Vg,Vf,Vw);
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  
  for(i = 0; i < Oparams.parnum; i++)
    {
      RCMx += rx[i];
      RCMy += ry[i];
      RCMz += rz[i];
    }
 
  sprintf(TXTA[1],"  BOX CM=(%.15f,%.15f,%.15f)\n", RCMx/Oparams.parnum, 
	  RCMy/Oparams.parnum, RCMz/Oparams.parnum);
  mdPrintf(ALL, TXTA[0], TXTA[1], NULL);
}

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  COORD_TYPE Drx, Dry, Drz;
  int i;
  
  DrSqTot = 0.0;

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
   }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
 
  Dtrans = DrSqTot / ( 6.0 * ((COORD_TYPE) Oparams.curStep) * 
		       ((COORD_TYPE) Oparams.parnum ) );   
  
  DrSqTot /= ((COORD_TYPE) Oparams.parnum);
  
}
#ifdef MD_LOADMESH
extern int mesh[100][150][3];
extern int ntripl[100];
#endif

/* =========================== >>> structFacta <<< =========================*/
void structFacts(void)
{
  /* DESCRIPTION:
     This mesuring function calculates the static structure factor */
  COORD_TYPE reRho, imRho, pi2;
  COORD_TYPE L, invNm, scalFact;
  COORD_TYPE rCMk, sumRho;
  COORD_TYPE *sumS, kbeg;
  int mp, i, Nm, n;
#ifndef MD_LOADMESH
  int mesh[][150][3]= 
#include "./kmesh.dat"
  int ntripl[]=
#include "./ntripl.dat"
#endif
  /* useful quantities */
  sumS = OprogStatus.sumS;
  L = Oparams.L;
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

/* =========================== >>> radDens <<< ============================= */
void radDens(void)
{
  int bin, Nm, i, j;
  int* hist;

  COORD_TYPE Vol, rij, rxij, ryij, rzij, rijSq;
  COORD_TYPE m, rlower, rupper, cost, nIdeal; 
  COORD_TYPE rhoAv, DELR;
  L = Oparams.L;
  invL = 1.0 / L;
  Vol = pow(L, 3.0);
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
