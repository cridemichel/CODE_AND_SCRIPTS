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
extern COORD_TYPE Vc[MAX_M], V[MAX_M], W[MAX_M], K[MAX_M], 
  WC[MAX_M], T1xx[MAX_M], T1yy[MAX_M], T1zz[MAX_M],
  T1xx[MAX_M], T1yy[MAX_M], T1zz[MAX_M], T1xy[MAX_M], 
  T1yz[MAX_M], T1zx[MAX_M], WCxy[MAX_M], WCyz[MAX_M], WCzx[MAX_M], 
  WCxx[MAX_M], WCyy[MAX_M], WCzz[MAX_M], Wxx[MAX_M], Wyy[MAX_M], 
  Wzz[MAX_M], Wxy[MAX_M], Wyz[MAX_M], Wzx[MAX_M], 
  Pxx[MAX_M], Pyy[MAX_M], Pzz[MAX_M], Pxy[MAX_M], Pyz[MAX_M], Pzx[MAX_M], 
  Patxy[MAX_M], Patyz[MAX_M], Patzx[MAX_M], Patxx[MAX_M], Patyy[MAX_M], 
  Patzz[MAX_M];  

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
  int mol, Nm, i, ss;
  COORD_TYPE px, py, pz;
  COORD_TYPE L, invL;

  Nm = Oparams.parnum;
  Etot = 0;

  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      E[ss] = K[ss] + Vc[ss];
      Etot += E[ss];
      /*  And now add the contribute due to the thermal bath */
      Etot += 0.5 * Sqr(s1[ss]) * OprogStatus.Q / Sqr(s[ss]) + 
	(3.0 * ((COORD_TYPE) Nm) - 3.0) * Oparams.T * log(s[ss]); 
      printf("[%d]:(%.2f,%.4f) ", ss, s[ss], s1[ss]);
    }
  printf("\n");
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */
  mol = 10;
  ss = 0;
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
     
      px = Oparams.m * vx[ss][i];
      py = Oparams.m * vy[ss][i];
      pz = Oparams.m * vz[ss][i];
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
  sprintf(TXTA[1], "  Etot=%.10f P[%d]=(%.15f,%.15f,%.15f)\n", 
	  Etot, ss, Px, Py, Pz);
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  
  loop(i, 1, Nm)
    {
      RCMx += rx[ss][i];
      RCMy += ry[ss][i];
      RCMz += rz[ss][i];
    }
 
  sprintf(TXTA[2],"  BOX CM[%d]=(%.15f,%.15f,%.15f)\n", ss, RCMx, RCMy, RCMz);

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
	  Drx = OprogStatus.sumVx[ss][i]; 
	  Dry = OprogStatus.sumVy[ss][i];
	  Drz = OprogStatus.sumVz[ss][i];
	  //sprintf(TXT,"i = %d\n", i);
	  //mdPrintf(STD, TXT, NULL);
	  DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
	}
      /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
	 because DrSq is zero! */
      
      Dtrans[ss] = DrSqTot / ( 6.0 * ((COORD_TYPE) Oparams.steplength) *
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
  
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      /* normalizzazione */
      norm = 0.0;
      for (iE = 0; iE < PE_POINTS; iE++)
	{
	  //printf("iE: %d\n", iE);
	  //printf("PE: %f\n",(double) OprogStatus.PE[ss][iE]);
	  norm = norm + ((double) OprogStatus.PE[ss][iE]);
	  //if (norm != 0) printf("opsPE: %d\n", norm);
	}
      //printf("norm: %f\n", norm);
      for (iE = 0; iE < PE_POINTS; iE++)
	{
	  if (norm != 0)
	    {
	      PE[ss*PE_POINTS + iE] = ((double) OprogStatus.PE[ss][iE]) / norm;
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
  COORD_TYPE Ktransl, m;
  int ss;

  //temp = 2.0 * Ktot / Oparams.PTM / (5.0 * Oparams.parnum - 3.0);
  
  temp = 2.0 * K[Oparams.PTM-1] / (3.0 * Oparams.parnum - 3.0);

  Ktransl = 0.0;
  m = Oparams.m;
  ss = 0;
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
      OprogStatus.sumPress[ss] += press[ss];
      press[ss] = OprogStatus.sumPress[ss] / NUMCALCS;
    }
  sprintf(TXT, "T:%.10f\n", temp);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);
}

/* ========================= >>> misuras <<< ============================== */
void misuras(void)
{
  int ss;
  C_T smed;
  smed = 0.0;
  for( ss = 0; ss < Oparams.PTM; ss++)
    {
      smed += s[ss];
    }
  
  smed /= (double) Oparams.PTM;
  printf ("smed: %f\n", smed);
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
  int ss;
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      if (OprogStatus.attempted[ss])
	probSc[ss] = ((double) OprogStatus.scambi[ss])/
	  ((double) OprogStatus.attempted[ss]);
    }
}
