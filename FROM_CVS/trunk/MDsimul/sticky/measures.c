#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
#ifdef MD_AB41
extern int *bondscache, *numbonds, **bonds, *numbonds0, **bonds0;
#endif
#ifdef MPI
extern int my_rank;
#endif
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
char TXTA[10][MSG_LEN];
char TXT[MSG_LEN];
extern double Vz;
/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern double pi, s1t, Vol1t, invL, s1p, Elrc, Plrc;   
extern double W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz;  

/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern double dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
extern void UpdateSystem(void);
/* ================================= */

/* ========================================================================= */
/*=========================== >>> vectProd <<< =========================== */
extern void vectProd(double r1x, double r1y, COORD_TYPE r1z, 
	 double r2x, double r2y, COORD_TYPE r2z, 
	 double* r3x, double* r3y, COORD_TYPE* r3z);
extern int *inCell[3], cellsx, cellsy, cellsz;
extern int *cellList;
extern int bound(int na, int n);
extern int *numbonds;
FILE* mf;
double calcpotene(void)
{
  double Epot; 
  int na;
#ifdef MD_AB41
  int kk, j, numbondsAA, numbondsAB;
#endif
#if 0
  double shift[NDIM];
  int evCode, cellRangeEne[2*NDIM], signDir[NDIM];
  int iX, iY, iZ, jX, jY, jZ, k, n;
#endif
  Epot = 0;
#if 0
  for (k = 0; k < NDIM; k++)
    { 
      cellRangeEne[2*k]   = - 1;
      cellRangeEne[2*k+1] =   1;
    }
  for (na = 0; na < Oparams.parnum; na++)
    {
      for (iZ = cellRangeEne[4]; iZ <= cellRangeEne[5]; iZ++) 
	{
	  jZ = inCell[2][na] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellsz - 1;    
	      shift[2] = - L;
	    } 
	  else if (jZ == cellsz) 
	    {
	      jZ = 0;    
	      shift[2] = L;
	    }
	  for (iY = cellRangeEne[2]; iY <= cellRangeEne[3]; iY ++) 
	    {
	      jY = inCell[1][na] + iY;    
	      shift[1] = 0.0;
	      if (jY == -1) 
		{
		  jY = cellsy - 1;    
		  shift[1] = -L;
		} 
	      else if (jY == cellsy) 
		{
		  jY = 0;    
		  shift[1] = L;
		}
	      for (iX = cellRangeEne[0]; iX <= cellRangeEne[1]; iX ++) 
		{
		  jX = inCell[0][na] + iX;    
		  shift[0] = 0.0;
		  if (jX == -1) 
		    {
		      jX = cellsx - 1;    
		      shift[0] = - L;
		    } 
		  else if (jX == cellsx) 
		    {
		      jX = 0;   
		      shift[0] = L;
		    }
		  n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;

		  for (n = cellList[n]; n > -1; n = cellList[n]) 
		    {
		      if (n < na) 
			{
			  if (bound(n,na))
			    Epot -= Oparams.bheight;
			}
		    }
		}
	    }
	}
    }
#else
 for (na = 0; na < Oparams.parnum; na++)
   {
#ifdef MD_AB41
     numbondsAA = numbondsAB = 0;
     for (kk=0; kk < numbonds[na]; kk++)
       {
	 j=bonds[na][kk]/(NA*NA);
	 if (na < Oparams.parnumA && j < Oparams.parnumA)
	   numbondsAA++;
	 else
	   numbondsAB++; 
       }
     Epot -= Oparams.bheightAA*numbondsAA;
     Epot -= Oparams.bheightAB*numbondsAB;
#else
     Epot -= Oparams.bheight*numbonds[na];
#endif
   }
#endif
#ifdef MD_GRAVITY
 for (na = 0; na < Oparams.parnum; na++)
   {
     /* add gravitational potential energy */
     Epot += 2.0*Oparams.m[(na < Oparams.parnumA)?0:1]*Oparams.ggrav*rz[na];
   }
#endif
 return 0.5*Epot;
}
void calcV(void)
{
  double tref;
 V = calcpotene();
 mf = fopenMPI(absMisHD("energy.dat"),"a");
#ifdef MD_BIG_DT
 tref = OprogStatus.refTime;
#else
 tref = 0.0;
#endif
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
   fprintf(mf, "%15G %.15G\n", Oparams.time + tref, V);
#else
#ifdef MD_AB41
   fprintf(mf, "%15G %.15G\n", Oparams.time + tref, V/((double)Oparams.parnum));
#else
 if (Oparams.parnumA < Oparams.parnum)
   fprintf(mf, "%15G %.15G\n", Oparams.time + tref, V/((double)Oparams.parnum-Oparams.parnumA));
 else
   fprintf(mf, "%15G %.15G\n", Oparams.time + tref, V/((double)Oparams.parnum));
#endif
#endif
#else
 fprintf(mf, "%15G %.15G\n", Oparams.time + tref, V/((double)Oparams.parnum));
#endif
 fclose(mf);
}
void calc_energy(char *msg);

/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
  double Px, Py, Pz, RCMx, RCMy, RCMz;
  int mol, Nm, i;
  double px, py, pz;
  double invL;

  Nm = Oparams.parnum;
  sprintf(TXTA[0],"STEP %lld [MEASURE]\n  L: %.10f\n", 
	  (long long int)Oparams.curStep, L);
  calc_energy("[MEASURES]"); 
  V = calcpotene();
  E = K + V;
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */
  mol = 10;
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  invL = 1.0 / L;
  /* 05/06/08: questo sembra essere pericoloso e non serve! */
  //UpdateSystem();
  for(i = 0; i < Oparams.parnumA; i++)
    {
      px = Oparams.m[0] * vx[i];
      py = Oparams.m[0] * vy[i];
      pz = Oparams.m[0] * vz[i];
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
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      px = Oparams.m[1] * vx[i];
      py = Oparams.m[1] * vy[i];
      pz = Oparams.m[1] * vz[i];
      Px += px;
      Py += py;
      Pz += pz;
    }
#ifdef MD_BIG_DT
  sprintf(TXTA[1], "t=%f E=%.15f V=%.15f P=(%.14G,%.14G,%.14G)\n", Oparams.time + OprogStatus.refTime,
	  E, V, Px, Py, Pz);
#else
  sprintf(TXTA[1], "t=%f E=%.15f V=%.15f P=(%.14G,%.14G,%.14G)\n", Oparams.time,
	  E, V, Px, Py, Pz);
#endif
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  
  for(i = 0; i < Oparams.parnumA; i++)
    {
      RCMx += rx[i]*Oparams.m[0];
      RCMy += ry[i]*Oparams.m[0];
      RCMz += rz[i]*Oparams.m[0];
    }
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      RCMx += rx[i]*Oparams.m[1];
      RCMy += ry[i]*Oparams.m[1];
      RCMz += rz[i]*Oparams.m[1];
    }
  sprintf(TXTA[2],"  BOX CM=(%.15f,%.15f,%.15f)\n", RCMx, RCMy, RCMz);
  //printf("RANK: %d STEP: %d\n", my_rank, Oparams.curStep);
  //fflush(stdout);
  mdPrintf(ALL, TXTA[0], TXTA[1], TXTA[2], NULL);
}
/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  double Drx, Dry, Drz, Dr4;
  int i;
  DrSqTot = 0.0;
#if 0
  Dr4 = 0.0;
#endif
  for(i=0; i < Oparams.parnumA; i++)
    {
      Drx = rx[i] - OprogStatus.rxCMi[i] + L*OprogStatus.DR[i][0]; 
      Dry = ry[i] - OprogStatus.ryCMi[i] + L*OprogStatus.DR[i][1];
      Drz = rz[i] - OprogStatus.rzCMi[i] + L*OprogStatus.DR[i][2];
#if 0
      if (OprogStatus.ipart == i)
	{
	  //sprintf(TXT,"i = %d\n", i);
	  //mdPrintf(STD, TXT, NULL);
	  /* Motion of the OprogStatus.ipart particle */
	  sqrtdr2 = sqrt(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
	}
#endif
      DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
#if 0
      Dr4 += Sqr(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
#endif
   }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
#ifdef MD_BIG_DT
  Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time + OprogStatus.refTime) *
		       ((double) Oparams.parnumA ) );   
#else
  Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time) *
		       ((double) Oparams.parnumA ) );   
#endif
  //printf("Dtr: %f\n", Dtrans);
#if 0
  Aa = ((double) Oparams.parnumA ) * 3.0 * 
    Dr4 / Sqr(DrSqTot) / 5.0 - 1.0; /* Non-Gaussian parameter */  
#endif
  DrSqTot /= ((double) Oparams.parnumA);
  DrSqTotA = DrSqTot;
#if 1
  mf = fopenMPI(absMisHD("msdA.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, DrSqTot);
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, DrSqTot);
#endif
  fclose(mf);
  DrSqTot = 0.0;
  if (Oparams.parnum > Oparams.parnumA)
    {
      for(i=Oparams.parnumA; i < Oparams.parnum; i++)
	{
	  Drx = rx[i] - OprogStatus.rxCMi[i] + L*OprogStatus.DR[i][0]; 
	  Dry = ry[i] - OprogStatus.ryCMi[i] + L*OprogStatus.DR[i][1];
	  Drz = rz[i] - OprogStatus.rzCMi[i] + L*OprogStatus.DR[i][2];
	  if (OprogStatus.ipart == i)
	    {
	      //sprintf(TXT,"i = %d\n", i);
	      //mdPrintf(STD, TXT, NULL);
	      /* Motion of the OprogStatus.ipart particle */
	      sqrtdr2 = sqrt(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
	    }
	  DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
	}
      mf = fopenMPI(absMisHD("msdB.dat"),"a");
#ifdef MD_BIG_DT
      fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime,
	      DrSqTot / ((double)(Oparams.parnum-Oparams.parnumA)));
#else
      fprintf(mf, "%15G %.15G\n", Oparams.time, DrSqTot / ((double)(Oparams.parnum-Oparams.parnumA)));
#endif
      DrSqTotB = DrSqTot / ((double)(Oparams.parnum-Oparams.parnumA));
      fclose(mf);
    }
#endif
}

/* ============================ >>> temperat <<< =========================== */
void temperat(void)
{
  double tempRot, tempTra;
#ifdef MD_THREESPOTS
  double dogTot, dogTra, dogRot;
#endif
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
#if 0
  int i;
  double m;
  K = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      K += Oparams.m[0]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
    }
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      K += Oparams.m[1]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
    }
  K *= 0.5;
  if (OprogStatus.brownian==1)
    temp = 2.0 * K / (3.0 * Oparams.parnum);
  else
    temp = 2.0 * K / (3.0 * Oparams.parnum - 3.0);
#endif
  calc_energy(NULL);
#ifdef MD_THREESPOTS
  if (OprogStatus.brownian)
    {
      dogTra = 3.0*Oparams.parnum;
      dogRot = 2.0*Oparams.parnumA + 3.0*(Oparams.parnum-Oparams.parnumA);
      dogTot = dogTra + dogRot; 
    }
  else
    {
      dogTra = 3.0*Oparams.parnum-3;
      dogRot = 2.0*Oparams.parnumA+3.0*(Oparams.parnum-Oparams.parnumA);
      dogTot = dogTra + dogRot; 
    }
  //printf("dogTot: %.15G dogTra: %.15G dogRot:%.15G\n", dogTot, dogTra, dogRot);
  temp = 2.0 * K / dogTot;
  tempTra =  2.0 * Ktra / dogTra;
  tempRot =  2.0 * Krot / dogRot;
#else
#ifdef MD_SPOT_OFF
  temp = 2.0 * K / (3.0 * Oparams.parnum - 3.0);
  tempTra =  2.0 * Ktra / (3.0 * Oparams.parnum);
  tempRot =  0.0;
#else
  temp = 2.0 * K / (6.0 * Oparams.parnum - 3.0);
  tempTra =  2.0 * Ktra / (3.0 * Oparams.parnum);
  tempRot =  2.0 * Krot / (3.0 * Oparams.parnum);
#endif
#endif
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
#if 1
  mf = fopenMPI(absMisHD("temp.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, temp, tempTra, tempRot);
#else
  fprintf(mf, "%15G %.15G %.15G %.15G\n", Oparams.time, temp, tempTra, tempRot);
#endif
  fclose(mf);
#endif
#if 0
  sprintf(TXT, "P:%.10f T:%.10f W: %10f\n", press, temp, W);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);
#endif
}
#ifdef MD_ROTDIFF_MIS
void calcrotMSD(void)
{
  FILE *fA, *fB;
  int i, a;
  double DphiA[3], DphiB[3];
  DphiSqA = DphiSqB = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      DphiSqA += Sqr(OprogStatus.sumox[i])+Sqr(OprogStatus.sumoy[i])+
	Sqr(OprogStatus.sumoz[i]);
    }
  DphiSqA /= ((double)Oparams.parnumA);
  if (Oparams.parnumA < Oparams.parnum)
    {
      for (i = Oparams.parnumA; i < Oparams.parnum; i++)
	{
	  DphiSqB += Sqr(OprogStatus.sumox[i])+Sqr(OprogStatus.sumoy[i])+
	    Sqr(OprogStatus.sumoz[i]);
	}
      DphiSqB /= ((double)Oparams.parnum - Oparams.parnumA);
    }
  for (a = 0; a < 3; a++)
    DphiA[a] = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      DphiA[0] += Sqr(OprogStatus.sumox[i]);
      DphiA[1] += Sqr(OprogStatus.sumoy[i]);
      DphiA[2] += Sqr(OprogStatus.sumoz[i]);
    } 
  for (a = 0; a < 3; a++)
    DphiA[a] /= ((double) Oparams.parnumA);
  if (Oparams.parnumA < Oparams.parnum)
    {
      for (a = 0; a < 3; a++)
	DphiB[a] = 0.0;
      for (i = Oparams.parnumA; i < Oparams.parnum; i++)
	{
       	  DphiB[0] += Sqr(OprogStatus.sumox[i]);
	  DphiB[1] += Sqr(OprogStatus.sumoy[i]);
	  DphiB[2] += Sqr(OprogStatus.sumoz[i]);
	} 
      for (a = 0; a < 3; a++)
	DphiB[a] /= ((double) Oparams.parnum-Oparams.parnumA);
    }
  DphiSq = DphiSqA;
  fA = fopenMPI(absMisHD("rotMSDA.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(fA,"%.15G %.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, DphiSqA, 
	  DphiA[0], DphiA[1], DphiA[2]); 
#else
  fprintf(fA,"%.15G %.15G %.15G %.15G %.15G\n", Oparams.time, DphiSqA, 
	  DphiA[0], DphiA[1], DphiA[2]); 
#endif
  fclose(fA);
  if (Oparams.parnum > Oparams.parnumA)
    {
      fB = fopenMPI(absMisHD("rotMSDB.dat"),"a");
      
#ifdef MD_BIG_DT
      fprintf(fB,"%.15G %.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, DphiSqB,
	      DphiB[0], DphiB[1], DphiB[2]); 
#else
      fprintf(fB,"%.15G %.15G %.15G %.15G %.15G\n", Oparams.time, DphiSqB,
	      DphiB[0], DphiB[1], DphiB[2]); 
#endif
      fclose(fB);
    }
  
}
#endif
/* ======================== >>> structFacts <<< =============================*/
void structFacts_OLD(void)
{
  int n, bin;
  double sum, pi4, rhoAv, r, Vol;
  double scalFact, kr, k, DELR;
  double* sumS;
  sumS = OprogStatus.sumS;
  DELR = (REND - RBEG) / MAXBIN;
  pi4 = pi * 4.0;
  Vol = L*L*L;  
  rhoAv = (double) Oparams.parnum / Vol; /* V = 1 here !!!!!!!!!!!!*/ 
  scalFact = (KEND - KBEG) / NUMK;
  
  loop(n, 1, NUMK)
    {
      k = KBEG + n * scalFact;
      sum = 0.0;
      /* Bode's rule */
      bin = 0;
     
      while(bin + 4 < MAXBIN)
	{
	  r = RBEG + bin * DELR;

	  kr = k*r;
	  sum += 14.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 64.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 24.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 64.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  bin++;
	  r += DELR;
	  kr = k*r;
	  sum += 14.0*(gr[bin] - 1.0) * Sqr(r) * sin(kr) / kr;
	  //bin++;
	}

      S[n] = 1.0 + pi4 * sum * rhoAv * DELR / 45.0;

      if (OprogStatus.avngS == 1)
	{
	  sumS[n] += S[n];
	  S[n] = sumS[n] / NUMCALCS;
	}
    }
}

/* =========================== >>> structFacta <<< =========================*/
void structFacts(void)
{
  /* DESCRIPTION:
     This mesuring function calculates the static structure factor */
  double reRho, imRho, pi2;
  double invNmA, scalFact;
  double rCMk, sumRho;
  double *sumS, kbeg;
  int mp, i, NmA, n;
  int mesh[][150][3]= 
#include "../bimix/kmesh.dat"
  int ntripl[]=
#include "../bimix/ntripl.dat"
  const int numk = 98;
  /* useful quantities */
  sumS = OprogStatus.sumS;
  invL = 1.0 / L;
  //scalFact = (KEND - KBEG) / NUMK;
  
  pi2 = 2.0 * pi;
  kbeg = 0.0; //pi2 * invL;
  scalFact = pi2 * invL;
  //printf("maxI : %d\n", maxI);
  /* We take 'maxI' values of Phi and 'maxI' values of Teta, so we have in this
     way NUMK2AV = maxI * maxI (about) values for 'k' over the 
     semi-sphere with same modulus */
  
  NmA = Oparams.parnumA;
  invNmA = 1.0 / NmA;

  for(n = 0; n < numk; n++)
    {
      sumRho = 0.0;
      //printf("nummp:%d\n", ntripl[n]);      
      for(mp = 0; mp < ntripl[n]; mp++)
	{
	  reRho = 0.0;
	  imRho = 0.0;
	  for(i=0; i < Oparams.parnumA; i++)
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

      S[n] = sumRho  * invNmA / ((double) ntripl[n]);  
    }
  if (OprogStatus.avngS == 1)
    {
      for(n=0; n < numk; n++)
	{
	  sumS[n] += S[n];
	  //if (n == 10) printf("somma di S:%f\n", sumS[n]);
	  S[n] = sumS[n] / NUMCALCS;
	}
    }
}


/* ============================= >>> maxwell <<< =========================== */
void maxwell(void)
{
  int n, i;
  int* histMB;
  double vMod, vSq;

  histMB = OprogStatus.histMB;
  
  if (OprogStatus.avngMB == 0)
    {
      for(i = 0; i < NUMV; i++)
	{
	  histMB[i] = 0;
	}
    } 
      
  for(i = 0; i < Oparams.parnumA; i++)
    {
      vSq = Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]);
      vMod = sqrt(vSq);
      n = ( NUMV - 1) * (vMod - VBEG)/ (VEND - VBEG); 
      /* VBEG must be always set to 0 */
      
      if( n < NUMV ) ++histMB[n];
    }
    
  for(i = 0; i < NUMV; i++)
	{
	  MB[i] = histMB[i];
	  if (OprogStatus.avngMB == 1)
	    MB[i] /= NUMCALCS; /* Mean if needed */
	}
}

/* =========================== >>> radDens <<< ============================= */
void radDens(void)
{
  int bin, Nm, i, j;
  int* hist;
  double Vol, rij, rxij, ryij, rzij, rijSq;
  double m, rlower, rupper, cost, nIdeal; 
  double rhoAv, DELR;
  invL = 1.0 / L;
  Nm = Oparams.parnumA;
  rhoAv = (double) Nm;
  m = Oparams.m[0];
  Vol = L*L*L;
  DELR = (L/2.0 - RBEG) / MAXBIN;

  hist = OprogStatus.hist;
  
  if (OprogStatus.avnggr == 0)
    {
      loop(i, 1, MAXBIN)
	{
	  hist[i] = 0;
	}
    }
  
  for(i = 0; i < Nm - 1; i++)
    {
      for(j = i+1; j < Nm; j++)
	{
	  rxij = rx[i] - rx[j];
	  ryij = ry[i] - ry[j];
	  rzij = rz[i] - rz[j];
	  rxij = rxij - L * rint(invL * rxij);
	  ryij = ryij - L * rint(invL * ryij);
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
      rlower = RBEG + ( (double) bin) * DELR;
      rupper = rlower + DELR;
      nIdeal = cost * (Sqr(rupper)*rupper - Sqr(rlower)*rlower);
      gr[bin] = ((double) hist[bin]) / ((double) Nm) / nIdeal;
      if (OprogStatus.avnggr == 1) 
	{
	  /* This a mean value of gr[bin] */
	  gr[bin] /= NUMCALCS;
	}
    }
}
/* ========================= >>> Ptensor <<< =============================== */
void calcpress(void)
{
  /* DESCRIPTION:
     Store the three off-diagonal terms of the pressure tensor in an array 
     to save on disk like a single measure */
#ifdef MD_HSVISCO

  mf = fopenMPI(absMisHD("press.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%.15G %.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime,
	  press, pressKin, pressHS, pressST);
#else
  fprintf(mf, "%.15G %.15G %.15G %.15G %.15G\n", Oparams.time, press, pressKin, pressHS, pressST);
#endif
  fclose(mf);
#endif

}


/* ========================= >>> Ptensor <<< =============================== */
void Ptensor(void)
{
  /* DESCRIPTION:
     Store the three off-diagonal terms of the pressure tensor in an array 
     to save on disk like a single measure */
#ifdef MD_HSVISCO
  mf = fopenMPI(absMisHD("Ptens.dat"),"a");
#endif
  Ptens[0] = Pxy;
  Ptens[1] = Pyz;
  Ptens[2] = Pzx;
#ifdef MD_HSVISCO
#ifdef MD_BIG_DT
  fprintf(mf, "%.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, Pxy, Pyz, Pzx);
#else
  fprintf(mf, "%.15G %.15G %.15G %.15G\n", Oparams.time, Pxy, Pyz, Pzx);
#endif
  fclose(mf);
#endif

}

/* ========================== >>> DQtensor <<< ============================= */
void DQtensor(void)
{
#ifdef MD_HSVISCO
  mf = fopenMPI(absMisHD("DQtens.dat"),"a");
#endif
  DQtens[0] = OprogStatus.DQxy;
  DQtens[1] = OprogStatus.DQyz;
  DQtens[2] = OprogStatus.DQzx;
#ifdef MD_HSVISCO
#ifdef MD_BIG_DT
  fprintf(mf, "%.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, OprogStatus.DQxy, OprogStatus.DQyz,  OprogStatus.DQzx);
#else
  fprintf(mf, "%.15G %.15G %.15G %.15G\n", Oparams.time, OprogStatus.DQxy, OprogStatus.DQyz,  OprogStatus.DQzx);
#endif
  fclose(mf);
#endif


}


