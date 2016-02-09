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
extern double Vz;
FILE* mf;
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
extern int *bondscache, **bonds;

/* neighbour list method variables */
extern double dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

/* ========================================================================= */
/*=========================== >>> vectProd <<< =========================== */
extern void vectProd(double r1x, double r1y, COORD_TYPE r1z, 
	 double r2x, double r2y, COORD_TYPE r2z, 
	 double* r3x, double* r3y, COORD_TYPE* r3z);
#if defined(MD_SQWELL) || defined(MD_INFBARRIER) 
extern int *inCell[3], cellsx, cellsy, cellsz;
extern int *cellList;
extern int bound(int na, int n);
extern int *numbonds;
double calcpotene(void)
{
#ifdef MD_MIXWDEPTH
  int i, j, kk;  
#endif  
  double shift[NDIM], Epot; 
  int cellRangeEne[2 * NDIM], signDir[NDIM], evCode,
      iX, iY, iZ, jX, jY, jZ, k, n, na;
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
#ifdef MD_MIXWDEPTH
  for (i = 0; i < Oparams.parnum; i++)
    {
      for (kk=0; kk < numbonds[i]; kk++)
	{
	  j = bonds[i][kk];
	  if (i < Oparams.parnumA && j < Oparams.parnumA)
	    Epot -= Oparams.bheight[0][0];
	  else if (i >= Oparams.parnumA && j >= Oparams.parnumA)
	    Epot -= Oparams.bheight[1][1];
	  else
	    Epot -= Oparams.bheight[0][1];
	}
    }
#else	
 for (na = 0; na < Oparams.parnum; na++)
    Epot -= numbonds[na];
#endif
#endif
 return 0.5*Epot;
}
#endif
void calcV(void)
{
#ifdef MD_SQWELL
  V = calcpotene();
  mf = fopenMPI(absMisHD("energy.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, V);
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, V);
#endif
  fclose(mf);
#else
  V = 0;
#endif
}
/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
  double Px, Py, Pz, RCMx, RCMy, RCMz;
  int mol, Nm, i;
  double px, py, pz;
  double invL;

#if 0  
  FILE* mf;
#endif
  Nm = Oparams.parnum;
#ifdef MD_GRAVITY
  sprintf(TXTA[0],"STEP %lld [MEASURE]\n  L: %.10f Lz: %.10f\n", 
	  (long long int)Oparams.curStep, L, Lz);
#else
  sprintf(TXTA[0],"STEP %lld [MEASURE]\n  L: %.10f\n", 
	  (long long int)Oparams.curStep, L);
#endif
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
  
#ifdef MD_GRAVITY
  E = K + V;
#else
  E = K;
#endif
#ifdef MD_SQWELL
  V = calcpotene();
  E = K + V;
#endif
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */
  mol = 10;
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  invL = 1.0 / L;
  UpdateSystem();
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
#ifdef MD_GRAVITY
  calcKVz();
  sprintf(TXTA[1], "t=%f E=%.15f P=(%.14G,%.14G,%.14G) Vz=%f\n", Oparams.time,
	  E, Px, Py, Pz, Vz);
#else
#ifdef MD_SQWELL
  sprintf(TXTA[1], "t=%f E=%.15f V=%.15f P=(%.14G,%.14G,%.14G)\n", Oparams.time,
	  E, V, Px, Py, Pz);
#else
   sprintf(TXTA[1], "t=%f E=%.15f P=(%.14G,%.14G,%.14G)\n", Oparams.time,
	  E, Px, Py, Pz);
#endif
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
#if 0
  mf = fopen(absMisHD("V.a"),"a");
  fprintf(mf, "%d %.15G\n", Oparams.curStep, V);
  fclose(mf);
#endif
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
  Dr4 = 0.0;

  for(i=0; i < Oparams.parnumA; i++)
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
      Dr4 += Sqr(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
   }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
 
  Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time) *
		       ((double) Oparams.parnumA ) );   
  Aa = ((double) Oparams.parnumA ) * 3.0 * 
    Dr4 / Sqr(DrSqTot) / 5.0 - 1.0; /* Non-Gaussian parameter */  
  
  DrSqTot /= ((double) Oparams.parnumA);
  mf = fopenMPI(absMisHD("msdA.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, DrSqTot);
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, DrSqTot);
#endif
  fclose(mf);
  if (Oparams.parnum == Oparams.parnumA)
    return;
 
  DrSqTot = 0.0;
  Dr4 = 0.0;

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
      Dr4 += Sqr(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
   }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
 
  Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time) *
		       ((double) Oparams.parnumA ) );   
  Aa = ((double) Oparams.parnumA ) * 3.0 * 
    Dr4 / Sqr(DrSqTot) / 5.0 - 1.0; /* Non-Gaussian parameter */  
  
  DrSqTot /= ((double) Oparams.parnumA);
  mf = fopenMPI(absMisHD("msdB.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, DrSqTot);
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, DrSqTot);
#endif
  fclose(mf);
}

/* ============================ >>> temperat <<< =========================== */
void temperat(void)
{
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
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
  mf = fopenMPI(absMisHD("T.dat"),"a");
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, temp);
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, temp);
#endif
  fclose(mf);

#if 0
  sprintf(TXT, "P:%.10f T:%.10f W: %10f\n", press, temp, W);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);
#endif
}

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
#if defined(MD_GRAVITY)
  Vol = Lz*Sqr(L);
#else
  Vol = L*L*L;  
#endif
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
  int n, i, Nm;
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
#ifdef MD_GRAVITY
  Vol = Sqr(L)*Lz;
#else
  Vol = L*L*L;
#endif 
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
#if !defined(MD_GRAVITY)
	  rzij = rzij - L * rint(invL * rzij);
#endif
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
void Ptensor(void)
{
  /* DESCRIPTION:
     Store the three off-diagonal terms of the pressure tensor in an array 
     to save on disk like a single measure */
  mf = fopenMPI(absMisHD("Ptens.dat"),"a");

  Ptens[0] = Pxy;
  Ptens[1] = Pyz;
  Ptens[2] = Pzx;
  fprintf(mf, "%.15G %.15G %.15G %.15G\n", Oparams.time, Pxy, Pyz, Pzx);
  fclose(mf);
}

/* ========================== >>> DQtensor <<< ============================= */
void DQtensor(void)
{
  DQtens[0] = OprogStatus.DQxy;
  DQtens[1] = OprogStatus.DQyz;
  DQtens[2] = OprogStatus.DQzx;
}


