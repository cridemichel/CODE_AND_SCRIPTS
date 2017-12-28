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
extern char TXT[MSG_LEN];
extern void links(int Nm, COORD_TYPE rcut, COORD_TYPE sigma);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void BuildNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigma);
extern void checkNebrRebuild(void);
extern void resetCM(int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);

extern void zeroArrays(C_T *arrx, C_T* arry, C_T *arrz, int N);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
COORD_TYPE pi, L, invL, s1p, Elrc, Plrc;   
COORD_TYPE Vc, V, W, K;
#if 0
  T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx,
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
  T1myz, T1mzx, T1mxx, T1myy, T1mzz;  
#endif
COORD_TYPE DrSq = 0.0, Mtot;
/* used by linked list routines */
int *head, *list, *map;  /* arrays of integer */
int NCell, mapSize, M;

/* neighbour list method variables */
COORD_TYPE dispHi;
int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
int **nebrTabInteracting;
/* ================================= */

/* DESCRIPTION:
   From Allen-Tildesley: 
   "When a molecule leaves the box by crossing one of the boundaries, 
   it is usual, to switch attention to the image molecule entering the box,
   by simply adding L, to or subtracting L from, the appropriate 
   coordinate." 
   
   NOTE: In the present case the box is cubic, but it is possible to 
   implement other geometries (see F.01)*/ 

/* ========================== >>> scalCor <<< ============================= */
void scalCor(int Nm)
{ 
  int i;
  COORD_TYPE DRx, DRy, DRz;
  
  L = Oparams.L;
  invL = 1.0 / L;
  
  /* Reduced particles to first box */
  for(i=0; i < Oparams.parnum; i++)
    {
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * rx[i]);
      DRy = - L * rint(invL * ry[i]);
      DRz = - L * rint(invL * rz[i]);
      
      rx[i] += DRx;
      ry[i] += DRy;
      rz[i] += DRz;
    }
}

const COORD_TYPE bc1 = 14.0/45.0, bc2 = 64.0/45.0, bc3 = 24.0/45.0;
  
/* =========================== >>> BodeTerm <<< ============================*/
COORD_TYPE BodeTerm(COORD_TYPE dt, COORD_TYPE* fi)
{
  return dt * (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	       bc1 * fi[4]);
}

/* ============================= >>> updatePE <<< ========================= */
void updatePE(int Nm)
{
  int iE;
  double ENmin, ENmax;
  
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  iE = (int) ((Vc - ((double) Nm)*ENmin) / 
	      ( ((double) Nm) * ((double) ENmax - ENmin)) * 
	      ((double) PE_POINTS));
  if ( (iE >= 0) && (iE < PE_POINTS)) 
    {
      ++(OprogStatus.PE[iE]);
    }

}
double ratio;
extern double WLJ;
extern void energyi(int i, double rxi, double ryi, double rzi, COORD_TYPE epsilon, 
	    COORD_TYPE sigma, COORD_TYPE rcut, double *Vlji, double *Vgi, double* Vfi,
	    double *Vw, int recursive);
double ranf(void);
/* ============================ >>> move<<< =================================*/
void move(void)
{
  double VljOLD, VgOLD, VfOLD, deltVw, deltVlj, deltVg, deltVf, BetaEdw, BetaAux;
  double VwNEW, VwOLD, VljNEW, VgNEW, VfNEW, deltVB, L2, L;
  double rxinew, ryinew, rzinew, rxiold, ryiold, rziold;
#if 0
  double VfOOLD;
#endif
  int i;
  /* DESCRIPTION:
     Move the particles by one step; this procedure is PARALLELIZED !!!
     Valid shared counters are: scN where N = 0..9 
     WARNING: "outer loops" should be ever parallelized, that is you must use 
     loopShr() macro for them inside move() and its subroutine 
     NOTE:
     linked list building could not be parallelized because father and child 
     could disturb each other */
  BetaEdw = 1.0/Oparams.T;
  BetaAux = 1.0/Oparams.Taux;
  Vf = OprogStatus.Vf;
  Vg = OprogStatus.Vg;
  Vlj= OprogStatus.Vlj;
  Vw = OprogStatus.Vw;
  L = Oparams.L;
  L2 = L/2.0;
  if (nebrNow)
    {
      /*printf("step: %d rebuilding neighbour list\n", Oparams.curStep);*/
      nebrNow = 0;
      dispHi = sqrt(3.0)*OprogStatus.drMax;
      /* build up linked list on predicted 
	 coordinates (only father do it)*/
      if (OprogStatus.noLinkedList)
	{
	  BuildNebrListNoLinked(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	}
      else
	{
	  links(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	  
	  /* Build up neighbour list */  
	  BuildNebrList(Oparams.parnum, Oparams.rcut, Oparams.sigma);
	}
    }


  for (i=0; i < Oparams.parnum; i++)
    {
      rxiold = rx[i];
      ryiold = ry[i];
      rziold = rz[i];
      
      energyi(i, rxiold, ryiold, rziold, Oparams.epsilon, Oparams.sigma, 
	      Oparams.rcut, &VljOLD, &VgOLD, &VfOLD, &VwOLD, 1);
   
      /*sumup(Oparams.parnum,  Oparams.epsilon, Oparams.sigma, Oparams.rcut);
	VfOOLD = Vf;*/
#if 0
      if (0 && Vlj == 0.0  && VljOLD > 0.0)
	{
	  printf("VljOLD: %e\n", VljOLD);
	  printf("Vlj: %e\n", Vlj);
	  /*exit(-1);*/
	}
#endif
      rxinew = rxiold + ( 2.0 * ranf() - 1.0 ) * OprogStatus.drMax;
      ryinew = ryiold + ( 2.0 * ranf() - 1.0 ) * OprogStatus.drMax;
      rzinew = rziold + ( 2.0 * ranf() - 1.0 ) * OprogStatus.drMax;
      
      /*printf("drMax[%d]: %.10g opdrmax: %.10g\n", i, drMax[i], OprogStatus.drMax);
      */
      /*if (sqrt(Sqr(rxinew-rxiold)+Sqr(ryinew-ryiold)+Sqr(rzinew-rziold))
	  > OprogStatus.rNebrShell*0.5)
	{
	  printf("errore spostamente troppo grande: %f\n", sqrt(Sqr(rxinew-rxiold)+Sqr(ryinew-ryiold)+Sqr(rzinew-rziold)));
	  printf("i=%d (%f,%f,%f)\n", i, rxiold, ryiold, rziold);
	  printf("DrMax:%f\n", OprogStatus.drMax);
	  exit(-1);
	}*/
#if 0
      /* scale to first box */
      rxinew = rxinew + L * rint(rxinew * invL);
      ryinew = ryinew + L * rint(ryinew * invL);
#endif
      energyi(i, rxinew, ryinew, rzinew, Oparams.epsilon, Oparams.sigma, 
	      Oparams.rcut, &VljNEW, &VgNEW, &VfNEW, &VwNEW, 2);
      deltVlj = VljNEW - VljOLD;
      deltVg  = VgNEW  - VgOLD;
      deltVf  = VfNEW  - VfOLD; 
      deltVw  = VwNEW  - VwOLD;
      deltVB  = BetaEdw * (deltVg+deltVlj+deltVw) + BetaAux * deltVf;
      if (deltVf > 0.0)
	printf("Bah deltVf: %f\n", deltVf);
     /* printf("step: %d deltVB: %f deltVg: %f deltVlj: %f, deltVw:%f deltVf: %f\n", 
	     Oparams.curStep, deltVB, deltVg, deltVlj, deltVw, deltVf);*/
      if (deltVB < 75.0)
	{
	  if (deltVB <= 0.0)
	    {
	      Vg += deltVg;
	      Vf += deltVf;
	      Vlj += deltVlj;
	      Vw += deltVw;
	      rx[i] = rxinew;
	      ry[i] = ryinew;
	      rz[i] = rzinew;
	      OprogStatus.accepted++;
	    }
	  else if ( exp(-deltVB) > ranf())
	    {
	      Vg += deltVg;
	      Vf += deltVf;
	      Vlj += deltVlj;
	      Vw += deltVw;
	      rx[i] = rxinew;
	      ry[i] = ryinew;
	      rz[i] = rzinew;
	      OprogStatus.accepted++;
	    }
	 /* sumup(Oparams.parnum,  Oparams.epsilon, Oparams.sigma, Oparams.rcut);
	  printf("deltaVf1: %.15e deltaVf: %.15e\n", deltVf, Vf - VfOOLD);*/
#if 0
      printf("PRIMA sumup i=%d F=(%.15e,%.15e,%.15e)\n",i, Fx[i],Fy[i], Fz[i]);
      sumup(Oparams.parnum,  Oparams.epsilon, Oparams.sigma, Oparams.rcut);
      printf("DOPO sumup i=%d F=(%.15e,%.15e,%.15e)\n", i, Fx[i],Fy[i], Fz[i]);
#endif

#if 0
	  if (VwNEW > 0.0 || VwOLD > 0.0)
	    {
	      printf("VwNEW: %e VwOLD: %e\n", VwNEW, VwOLD);
	    }
	   /*if (fabs(deltVf) > 1E2|| deltVf <0)
	printf("passo: %d AZZZZZZZOO: %.10e\n", Oparams.curStep, deltVf);
	*/
	  if (0 && (Vf < 0.0||Vlj <0))
	    {
	      printf("==> Vf: %e Vlf:%e\n", Vf, Vlj);
	      printf("Vf: %e deltVf:%e, VfNEW: %e VfOLD:%e\n", Vf-deltVf, deltVf, VfNEW, VfOLD);
	      printf("Vlj: %e deltVlj: %e\n", Vlj-deltVlj, deltVlj);
	      /*printf("i=%d rz[i]:%e rzwall:%e\n", i, rz[i], 1.0-Oparams.Lz*0.5-Oparams.sigma*0.5);
	      */
	      sumup(Oparams.parnum,  Oparams.epsilon, Oparams.sigma, Oparams.rcut);
	      printf("after sumup Vlj:%e Vf: %e\n", Vlj, Vf);
	      printf("VljOLD:%e VljNEW: %e\n", VljOLD,  VljNEW);
	      energyi(i, rxiold, ryiold, rziold, Oparams.epsilon, Oparams.sigma, 
		      Oparams.rcut, &VljOLD, &VgOLD, &VfOLD, &VwOLD, 1);
	      printf("dopo energy i: VljOLD:%e\n", VljOLD);
	      printf("DrMax: %f Vlj: %e Vw: %e Vg:%e\n", OprogStatus.drMax, Vlj, Vw, Vg);
       	      printf("rzold: %e rznew: %e\n", rziold,  rzinew);
	      printf("F = (%e,%e,%e)\n", Fx[i], Fy[i], Fz[i]);
	      
	      exit(1);
	    }
#endif
	}
      OprogStatus.attempted++;
    }
  if (OprogStatus.scaleTauxSteps && (Oparams.curStep % OprogStatus.scaleTauxSteps == 0))
    {
      if (OprogStatus.scaleMode==0)
	Oparams.Taux -= OprogStatus.deltaTaux;
      else
	Oparams.Taux /= OprogStatus.deltaTaux;
    }
 
  
  if (OprogStatus.iratio!=0 && Oparams.curStep % abs(OprogStatus.iratio) == 0)
    {
      ratio = ((double)OprogStatus.accepted) / 
	((double)abs(OprogStatus.iratio)*Oparams.parnum);
      /*printf("current acceptance ratio: %f DrMax:%.8g\n", ratio, 
	OprogStatus.drMax);*/
      if (OprogStatus.iratio > 0)
	{
	  if (ratio > 0.5)
	    {
	      OprogStatus.drMax *= 1.05;
	      /*
		 printf("drMax: %f drmaxsqrt3: %f rnebrsh*0.5:%f\n", OprogStatus.drMax,
		 sqrt(3.0)*OprogStatus.drMax, OprogStatus.rNebrShell*0.5 );
		 */
	      if (sqrt(3.0)*OprogStatus.drMax > OprogStatus.rNebrShell*0.5)
		{
		  /*printf("rejecting new drMax...\n");*/
		  OprogStatus.drMax /= 1.05;
		}
	    }
	  else
	    {
	      OprogStatus.drMax *= 0.95;
	    }
	}
      OprogStatus.accepted = 0;
    }
    
    checkNebrRebuild();
    Vtot = Vg + Vf + Vlj + Vw;
    OprogStatus.Vf = Vf;
    OprogStatus.Vg = Vg;
    OprogStatus.Vlj = Vlj;
    OprogStatus.Vw = Vw;
#if 0
    if (Oparams.curStep%500==0)
      {
	sumup(Oparams.parnum,  Oparams.epsilon, Oparams.sigma, Oparams.rcut);
	printf("DVf: %.15e DVg: %.15e DVlj: %15e DVw: %15e\n", (OprogStatus.Vf - Vf)/Vf,
	     (OprogStatus.Vg-Vg)/Vg, (OprogStatus.Vlj-Vlj)/Vlj, 
	     (OprogStatus.Vw - Vw)/Vf);
      }
#endif 
    /*checkNebrList(Oparams.parnum, Oparams.rcut, Oparams.sigma);*/
    /*if (fabs(OprogStatus.Vf - Vf)> 1E-6)
      {
      printf("step: %d\n", Oparams.curStep);
      printf("Oprog.Vf:%.10e Vf:%.10e\n", OprogStatus.Vf, Vf);
      exit(-1);
      }*/
    if (  ( (OprogStatus.CMreset > 0) &&
	  (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	     (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
      resetCM(Oparams.parnum);

  /* Update the integral of the pressure tensor */
  /*updateDQ(Oparams.steplength);*/
  
  /*updatePE(Oparams.parnum);*/
}
