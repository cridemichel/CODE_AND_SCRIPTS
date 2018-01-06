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
extern C_T Ktot;
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
char TXTA[10][MSG_LEN];
char TXT[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL;   
extern COORD_TYPE Vc, V, W, K, 
  WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, 
  T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, 
  Wzz, Wxy, Wyz, Wzx, 
  Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, 
  Patxy, Patyz, Patzx, Patxx, Patyy, 
  Patzz;  

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
  COORD_TYPE Px, Py, Pz, RCMx, RCMy, RCMz;
  int mol, Nm, i;
  COORD_TYPE px, py, pz;
  COORD_TYPE L, invL;

  Nm = Oparams.parnum;

  E = K + Vc;
  //Etot += E;
  /*  And now add the contribute due to the thermal bath */
  E += 0.5 * Sqr(s1) * OprogStatus.Q / Sqr(s) + 
    (3.0 * ((COORD_TYPE) Nm) - 3.0) * Oparams.T * log(s); 
  //printf("(%.2f,%.4f) ", s, s1);
  
  //printf("\n");
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */
  mol = 10;
  //printf("Elrc:: %.6f Plrc: %.6f\n", Elrc, Plrc);
  //printf("Atom position: %.15f\n", rx[mol]);
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  //L = cbrt(Vol);
  //invL = 1.0 / L;

  loop(i, 1, Nm)
    {
      
      // printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
     
      px = Oparams.m * vx[i];
      py = Oparams.m * vy[i];
      pz = Oparams.m * vz[i];
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
  sprintf(TXTA[1], "STEP %d | Vc:%f E=%.10f P=(%.15f,%.15f,%.15f)\n", 
	  Oparams.curStep, Vc, E, Px, Py, Pz);
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
  mdPrintf(STD, TXTA[1], TXTA[2], NULL);
}

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  COORD_TYPE Drx, Dry, Drz;
  int i;
  int ss;
  
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      DrSqTot = 0.0;
      for(i = 0; i < Oparams.parnum; i++)
	{
	  Drx = OprogStatus.sumVx[i]; 
  	  Dry = OprogStatus.sumVy[i];
	  Drz = OprogStatus.sumVz[i];
	  //sprintf(TXT,"i = %d\n", i);
	  //mdPrintf(STD, TXT, NULL);
	  DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
	}
      /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
	 because DrSq is zero! */
      
      Dtrans = DrSqTot / ( 6.0 * ((COORD_TYPE) Oparams.steplength) *
			   ((COORD_TYPE) Oparams.curStep) * 
			   ((COORD_TYPE) Oparams.parnum ) );   
      //DrSqTot /= ((COORD_TYPE) Oparams.parnum);
    }
}

/* ============================== >>> PE <<< =============================== */
void probE(void)
{
  int ss, iE;
  double norm;

  /* Questa misura va fatta dopo un controllo dell'equilibratura poiche'
     diversamente non PE contiene valori aggiornati */
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      /* normalizzazione */
      norm = 0.0;
      for (iE = 0; iE < PE_POINTS; iE++)
	{
	  //printf("iE: %d\n", iE);
	  //printf("PE: %f\n",(double) OprogStatus.PE[ss][iE]);
	  norm = norm + ((double) PEint[ss][iE]);
	  //if (norm != 0) printf("opsPE: %d\n", norm);
	}
      //printf("norm: %f\n", norm);
      for (iE = 0; iE < PE_POINTS; iE++)
	{
	  if (norm != 0)
	    {
	      PE[ss*PE_POINTS + iE] = ((double) PEint[ss][iE]) / norm;
	      //printf("opsPE: %d norm: %d\n", OprogStatus.PE[ss][iE], norm);
	      //printf("PE[%d][%d]:%f\n", ss, iE, PE[ss*PE_POINTS + iE]);
	    }
	}
    }
}

/* ============================ >>> temperat <<< =========================== */
void temperat(void)
{
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
  COORD_TYPE m;
  
  //temp = 2.0 * Ktot / Oparams.PTM / (5.0 * Oparams.parnum - 3.0);
  
  temp = 2.0 * K / (3.0 * Oparams.parnum - 3.0);

  m = Oparams.m;
 
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
  sprintf(TXT, "s:%.6f s1: %.6f T:%.10f\n", s, s1, temp);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);
}

/* ======================= >>> media << ============================ */
void media(void)
{
  /* pressure */
  if (OprogStatus.avgMedia == 1)
    {
      OprogStatus.sumMedia += diffMedia;
      diffMedia = OprogStatus.sumMedia / NUMCALCS;
    }

  sprintf(TXT, "media:%.10f\n", diffMedia);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);

}

/* ============================= >>> probScMis <<< ==================== */
void probScMis(void)
{
  if (OprogStatus.attempted[Oparams.lambdat[my_rank]])
    probSc = ((double) OprogStatus.scambi[Oparams.lambdat[my_rank]])/
      ((double) OprogStatus.attempted[Oparams.lambdat[my_rank]]);
}
