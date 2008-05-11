#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/
#define MD_DEBUG21(x)
/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
#ifdef MPI
extern int my_rank;
#endif
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
char TXTA[10][MSG_LEN];
char TXT[MSG_LEN];
extern double Vz;
FILE *mf;
#ifdef MD_INELASTIC
FILE *mf2;
#endif
extern double ***R;
void radius_of_gyration(void);
extern void UpdateSystem(void);
double DphiSqA=0.0, DphiSqB=0.0, DrSqTotA=0.0, DrSqTotB=0.0;
/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern double pi, s1t, Vol1t, s1p, Elrc, Plrc;   
#ifdef MD_LXYZ
extern double invL[3];
#else
extern double invL;
#endif
extern double W, K, WC, T1xx, T1yy, T1zz, Ktra, Krot,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz;  
#ifdef EDHE_FLEX
#ifdef MD_LL_BONDS
extern long long int *bondscache, **bonds;
#else
extern int *bondscache, **bonds;
#endif
#endif
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern double dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

/* ========================================================================= */
/*=========================== >>> vectProd <<< =========================== */
extern void vectProd(double r1x, double r1y, COORD_TYPE r1z, 
	 double r2x, double r2y, COORD_TYPE r2z, 
	 double* r3x, double* r3y, COORD_TYPE* r3z);
#ifdef MD_PATCHY_HE 
extern int *inCell[3], cellsx, cellsy, cellsz;
extern int *cellList;
extern int bound(int na, int n);
extern int *numbonds;
#ifdef EDHE_FLEX
extern int is_in_ranges(int A, int B, int nr, rangeStruct* r);
#endif
double calcpotene(void)
{
  double Epot; 
  int na;
#ifdef EDHE_FLEX
#ifdef MD_LL_BONDS
  long long int jj, jj2, aa, bb;
  int kk, kk2;
#else
  int kk2, jj, kk, jj2, aa, bb;
#endif
#endif
#if 0
  double shift[NDIM];
  int cellRangeEne[2 * NDIM];
  int iX, iY, iZ, jX, jY, jZ, k, n, signDir[NDIM], evCode;
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
  MD_DEBUG21(printf("BEGIN\n"));
  for (na = 0; na < Oparams.parnum; na++)
    {
#ifdef EDHE_FLEX
      for (kk=0; kk < numbonds[na]; kk++)
	{
	  jj = bonds[na][kk]/(NANA);
	  jj2 = bonds[na][kk]%(NANA);
	  aa = jj2 / NA;
	  bb = jj2 % NA;
	  //printf("numbonds[%d]=%d aa=%d bb=%d\n", na, numbonds[na], aa, bb);
	  for (kk2 = 0; kk2 < Oparams.ninters; kk2++)
	    {
	      if ( (is_in_ranges(typeOfPart[na], intersArr[kk2].type1, intersArr[kk2].nr1, intersArr[kk2].r1) && 
		    is_in_ranges(typeOfPart[jj], intersArr[kk2].type2, intersArr[kk2].nr2, intersArr[kk2].r2) &&
		    intersArr[kk2].spot1 == aa-1 && intersArr[kk2].spot2 == bb-1) || 
		   (is_in_ranges(typeOfPart[jj], intersArr[kk2].type1, intersArr[kk2].nr1, intersArr[kk2].r1) && 
		    is_in_ranges(typeOfPart[na], intersArr[kk2].type2, intersArr[kk2].nr2, intersArr[kk2].r2) &&
		    intersArr[kk2].spot1 == bb-1 && intersArr[kk2].spot2 == aa-1) )  
		{
		  MD_DEBUG21(printf("(%d,%d)-(%d,%d) height=%.15G\n", na, aa-1, jj, bb-1, intersArr[kk2].bheight));
		  Epot -= intersArr[kk2].bheight;
		}		 
	    }
	  if (Oparams.nintersIJ > 0)
	    {
	      for (kk2 = 0; kk2 < Oparams.nintersIJ; kk2++)
		{
		  if ( (is_in_ranges(na, intersArrIJ[kk2].i, intersArrIJ[kk2].nr1, intersArrIJ[kk2].r1) && 
			is_in_ranges(jj, intersArrIJ[kk2].j, intersArrIJ[kk2].nr2, intersArrIJ[kk2].r2) &&
		      	intersArrIJ[kk2].spot1 == aa-1 && intersArrIJ[kk2].spot2 == bb-1) || 
		       (is_in_ranges(jj, intersArrIJ[kk2].i, intersArrIJ[kk2].nr1, intersArrIJ[kk2].r1) && 
			is_in_ranges(na, intersArrIJ[kk2].j, intersArrIJ[kk2].nr2, intersArrIJ[kk2].r2) &&
			intersArrIJ[kk2].spot1 == bb-1 && intersArrIJ[kk2].spot2 == aa-1) )  
		    {
		      MD_DEBUG21(printf("(%d,%d)-(%d,%d) height=%.15G\n", na, aa-1, jj, bb-1, intersArrIJ[kk2].bheight));
		      Epot -= intersArrIJ[kk2].bheight;
		    }
		}	
	    }
	}
      //Epot -= numbonds[na];
#else
      Epot -= numbonds[na];
#endif
    }
  MD_DEBUG21(printf("END\n"));
#endif
 return 0.5*Epot;
}
#endif
#ifdef MD_RABBIT
extern int getnumbonds(int np, interStruct *ts, int inverted);
void get_bimono_bonds(int *bulk, int *mono, int *bi)
{
  int nb, i, ti;
  interStruct ts;
  *bulk = 0;
  *mono = 0;
  *bi = 0;
  nb=0;
  for (i=0; i < Oparams.parnum; i++)
    {
      ti = typeOfPart[i];
      if (ti==0 || ti==1)
	{
	  ts.type1 = ti;
	  ts.type2 = 5;
	  ts.spot1 = 1;
	  ts.spot2 = 0; 
	  if (ti==0)
	    nb = getnumbonds(i, &ts, 0);
	  else
	    nb += getnumbonds(i, &ts, 0);
	  if (ti==1)
	    {
	      if (nb==1) 
		(*mono)++;
	      else if (nb==2)
		(*mono)++;
	      else
		(*bulk)++;
	    }
	}
    }
}
#endif
void calcV(void)
{
#ifdef MD_RABBIT
  int bulk, mono, bi;
#endif
#ifdef MD_PATCHY_HE
#ifdef MD_FOUR_BEADS
  radius_of_gyration();
#endif
  V = calcpotene();
  mf = fopenMPI(absMisHD("energy.dat"),"a");
#if 0
  if (Oparams.parnumA < Oparams.parnum)
    fprintf(mf, "%15G %.15G\n", Oparams.time, V/((double)Oparams.parnum-Oparams.parnumA));
  else
    fprintf(mf, "%15G %.15G\n", Oparams.time, V/((double)Oparams.parnum));
#else
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, V/((double)Oparams.parnum));
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, V/((double)Oparams.parnum));
#endif
#endif
 fclose(mf);
#else
  V = 0;
#endif
#ifdef MD_RABBIT
  /* Salva il numero di legami monovalenti e bivalenti */
  mf = fopenMPI(absMisHD("bi-mono-bonds.dat"),"a");
  get_bimono_bonds(&bulk, &mono, &bi);
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %d %d %d\n", Oparams.time + OprogStatus.refTime, bulk, mono, bi);
#else
  fprintf(mf, "%15G %d %d %d\n", Oparams.time, bulk, mono, bi);
#endif
  fclose(mf);
#endif
}
void calc_energy(char *msg);
extern double *angM;
#ifdef EDHE_FLEX
void calc_momentum_filtered(double P[3], int filter)
{
  double mass;
  int i;
  double px, py, pz;
  UpdateSystem();
  P[0] = 0.0;
  P[1] = 0.0;
  P[2] = 0.0;
  //invL = 1.0 / L;
  for(i = 0; i < Oparams.parnum; i++)
    {
      if (filter != 0 && filter != typesArr[typeOfPart[i]].brownian)
	continue;
      mass = typesArr[typeOfPart[i]].m;
      px = mass * vx[i];
      py = mass * vy[i];
      pz = mass * vz[i];
      P[0] += px;
      P[1] += py;
      P[2] += pz;
    }
}
#endif
/* ============================== >>> Energy <<< ============================*/
void energy(void)
{
  /* DESCRIPTION:
     This measuring function calculate the total energy of the system */
#ifdef EDHE_FLEX
  double mass, Mtot;
  int cc=0;
#endif
  double Px, Py, Pz, RCMx, RCMy, RCMz;
  int mol, Nm, i;
  double px, py, pz;
//  double invL;
  double tref;
#ifdef MD_BIG_DT
  tref = OprogStatus.refTime;
#else
  tref = 0.0;
#endif
#if 0  
  FILE* mf;
#endif
  Nm = Oparams.parnum;
#ifdef MD_LXYZ
  sprintf(TXTA[0],"STEP %lld [MEASURE]\n  Lx: %.10f Ly: %.10f Lz:%.15f\n", 
	  (long long int)Oparams.curStep, L[0], L[1], L[2]);
#else
#ifdef MD_GRAVITY
  sprintf(TXTA[0],"STEP %lld [MEASURE]\n  L: %.10f Lz: %.10f\n", 
	  (long long int)Oparams.curStep, L, Lz);
#else
  sprintf(TXTA[0],"STEP %lld [MEASURE]\n  L: %.10f\n", 
	  (long long int)Oparams.curStep, L);
#endif
#endif
  calc_energy("[MEASURES]"); 
  //printf("angM[0]:%.15G angM[10]:%.15G\n", angM[0], angM[10]);
#ifdef MD_GRAVITY
  E = K + V;
#else
  E = K;
#endif
#ifdef MD_PATCHY_HE
  V = calcpotene();
  E = K + V;
#endif
  /* So now E is the extended hamiltonian that should be an integral of 
     motion */
  mol = 10;
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  //invL = 1.0 / L;
  UpdateSystem();
  for(i = 0; i < Oparams.parnumA; i++)
    {
#ifdef EDHE_FLEX
      mass = typesArr[typeOfPart[i]].m;
      px = mass * vx[i];
      py = mass * vy[i];
      pz = mass * vz[i];
#else
      px = Oparams.m[0] * vx[i];
      py = Oparams.m[0] * vy[i];
      pz = Oparams.m[0] * vz[i];
#endif
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
#ifdef EDHE_FLEX
      mass = typesArr[typeOfPart[i]].m;
      px = mass * vx[i];
      py = mass * vy[i];
      pz = mass * vz[i];
#else
      px = Oparams.m[1] * vx[i];
      py = Oparams.m[1] * vy[i];
      pz = Oparams.m[1] * vz[i];
#endif
      Px += px;
      Py += py;
      Pz += pz;
    }
#ifdef MD_GRAVITY
  calcKVz();
  sprintf(TXTA[1], "t=%f E=%.15f P=(%.14G,%.14G,%.14G) Vz=%f\n", Oparams.time + tref,
	  E, Px, Py, Pz, Vz);
#else
#ifdef MD_PATCHY_HE
  sprintf(TXTA[1], "t=%f E=%.15f V=%.15f P=(%.14G,%.14G,%.14G)\n", Oparams.time + tref,
	  E, V, Px, Py, Pz);
#else
   sprintf(TXTA[1], "t=%f E=%.15f P=(%.14G,%.14G,%.14G)\n", Oparams.time + tref,
	  E, Px, Py, Pz);
#endif
#endif
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
#ifdef EDHE_FLEX
  cc = 0; 
  Mtot = 0.0;	
#endif 
  for(i = 0; i < Oparams.parnumA; i++)
    {
#ifdef EDHE_FLEX
      mass = typesArr[typeOfPart[i]].m;
      if (mass > MD_INF_MASS)
	{
	  cc++;
	  continue;
	}
      Mtot += mass;
      RCMx += mass * rx[i];
      RCMy += mass * ry[i];
      RCMz += mass * rz[i];
#else 
      RCMx += rx[i]*Oparams.m[0];
      RCMy += ry[i]*Oparams.m[0];
      RCMz += rz[i]*Oparams.m[0];
#endif
    }
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
#ifdef EDHE_FLEX
      mass = typesArr[typeOfPart[i]].m;
      RCMx += mass * rx[i];
      RCMy += mass * ry[i];
      RCMz += mass * rz[i];
#else 
      RCMx += rx[i]*Oparams.m[1];
      RCMy += ry[i]*Oparams.m[1];
      RCMz += rz[i]*Oparams.m[1];
#endif
    }
#ifdef EDHE_FLEX
  RCMx /= Mtot;
  RCMy /= Mtot;
  RCMz /= Mtot;
#endif
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
#ifdef EDHE_FLEX
void radius_of_gyration(void)
{
  int i, kk, Namino;
  double rcx, rcy, rcz, dx, dy, dz, rmean[3], rgyr;
  FILE *f;
  f = fopenMPI(absMisHD("radius_of_gyration.dat"),"a");

  Namino = typeNP[0];
  /* evaluate center of mass of the protein */
  rmean[0] = rcx = rx[0];
  rmean[1] = rcy = ry[0];
  rmean[2] = rcz = rz[0];
  /* evaluate rad of gyr just for aminoacids not for plates */
  for (i = 0; i < Namino; i++)
    {
      dx = rx[i+1]-rx[i];
      dy = ry[i+1]-ry[i];
      dz = rz[i+1]-rz[i];
      /* minimum image */
#ifdef MD_LXYZ
      dx = dx - L[0]*rint(dx/L[0]);
      dy = dy - L[1]*rint(dy/L[1]);
      dz = dz - L[2]*rint(dz/L[2]);
#else
      dx = dx - L*rint(dx/L);
      dy = dy - L*rint(dy/L);
      dz = dz - L*rint(dz/L);
#endif
      rcx = rcx + dx;
      rcy = rcy + dy;
      rcz = rcz + dz;
      rmean[0] += rcx;
      rmean[1] += rcy;
      rmean[2] += rcz;
    }
  for (kk=0; kk < 3; kk++)
    rmean[kk] /= ((double)Namino);
  rgyr = 0.0;
  for (i = 0; i < Namino; i++)
    {
      dx = rx[i]-rmean[0];
      dy = ry[i]-rmean[1];
      dz = rz[i]-rmean[2];
      /* minimum image */
#ifdef MD_LXYZ
      dx = dx - L[0]*rint(dx/L[0]);
      dy = dy - L[1]*rint(dy/L[1]);
      dz = dz - L[2]*rint(dz/L[2]);
#else
      dx = dx - L*rint(dx/L);
      dy = dy - L*rint(dy/L);
      dz = dz - L*rint(dz/L);
#endif
      rgyr += Sqr(dx);
      rgyr += Sqr(dy);
      rgyr += Sqr(dz);
    }
  rgyr /= ((double)Namino);
  rgyr = sqrt(rgyr);
#ifdef MD_BIG_DT
  fprintf(f, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, rgyr);
#else
  fprintf(f, "%15G %.15G\n", Oparams.time, rgyr);
#endif
  fclose(f);
}
#endif

/* ========================== >>> transDiff <<< =============================*/
void transDiff(void)
{
  FILE *f;
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
#ifdef MD_CALC_DPP
  double MSDx, MSDy, MSDz;
#endif
  double Drx, Dry, Drz, Dr4;
  int i;
  DrSqTot = 0.0;
  Dr4 = 0.0;
  f = fopen(absMisHD("MSDA.dat"), "a");
  for(i=0; i < Oparams.parnumA; i++)
    {
#ifdef MD_LXYZ
      Drx = rx[i] - OprogStatus.rxCMi[i] + L[0]*OprogStatus.DR[i][0]; 
      Dry = ry[i] - OprogStatus.ryCMi[i] + L[1]*OprogStatus.DR[i][1];
      Drz = rz[i] - OprogStatus.rzCMi[i] + L[2]*OprogStatus.DR[i][2];
#else
      Drx = rx[i] - OprogStatus.rxCMi[i] + L*OprogStatus.DR[i][0]; 
      Dry = ry[i] - OprogStatus.ryCMi[i] + L*OprogStatus.DR[i][1];
      Drz = rz[i] - OprogStatus.rzCMi[i] + L*OprogStatus.DR[i][2];
#endif
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
      //Dr4 += Sqr(Sqr(Drx) + Sqr(Dry) + Sqr(Drz));
   }
  DrSqTotA =  DrSqTot / ((double) Oparams.parnumA);
#ifdef MD_BIG_DT
  fprintf(f, "%.15G %.15G\n", Oparams.time + OprogStatus.refTime, DrSqTotA);
#else
  fprintf(f, "%.15G %.15G\n", Oparams.time, DrSqTotA);
#endif
  fclose(f);
  if (Oparams.parnumA < Oparams.parnum)
    {
      f = fopen(absMisHD("MSDB.dat"), "a");
      DrSqTot = 0.0;
      for(i=Oparams.parnumA; i < Oparams.parnum; i++)
	{
#ifdef MD_LXYZ
	  Drx = rx[i] - OprogStatus.rxCMi[i] + L[0]*OprogStatus.DR[i][0]; 
	  Dry = ry[i] - OprogStatus.ryCMi[i] + L[1]*OprogStatus.DR[i][1];
	  Drz = rz[i] - OprogStatus.rzCMi[i] + L[2]*OprogStatus.DR[i][2];
#else
	  Drx = rx[i] - OprogStatus.rxCMi[i] + L*OprogStatus.DR[i][0]; 
	  Dry = ry[i] - OprogStatus.ryCMi[i] + L*OprogStatus.DR[i][1];
	  Drz = rz[i] - OprogStatus.rzCMi[i] + L*OprogStatus.DR[i][2];
#endif
	  DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
	}
      
      DrSqTotB = DrSqTot / ((double)Oparams.parnum - Oparams.parnumA);
#ifdef MD_BIG_DT
      fprintf(f, "%.15G %.15G\n", Oparams.time + OprogStatus.refTime, DrSqTotB);
#else
      fprintf(f, "%.15G %.15G\n", Oparams.time, DrSqTotB);
#endif
      fclose(f);
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
  DrSqTot = DrSqTotA;
#ifdef MD_CALC_DPP
  f = fopen(absMisHD("MSDAxyz.dat"), "a");
  MSDx = MSDy = MSDz = 0.0;
  for(i=0; i < Oparams.parnumA; i++)
    {
      MSDx += Sqr(OprogStatus.sumdx[i]);
      MSDy += Sqr(OprogStatus.sumdy[i]);
      MSDz += Sqr(OprogStatus.sumdz[i]);
    }
  MSDx /= Oparams.parnumA;
  MSDy /= Oparams.parnumA;
  MSDz /= Oparams.parnumA;
#ifdef MD_BIG_DT
  fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, MSDx+MSDy+MSDz, MSDx, MSDy, MSDz);
#else
  fprintf(f, "%.15G %.15G %.15G %.15G %.15G\n", Oparams.time, MSDx+MSDy+MSDz, MSDx, MSDy, MSDz);
#endif
  fclose(f);
  if (Oparams.parnumA < Oparams.parnum)
    {
      f = fopen(absMisHD("MSDBxyz.dat"), "a");
      for(i=Oparams.parnumA; i < Oparams.parnum; i++)
	{
	  MSDx += Sqr(OprogStatus.sumdx[i]);
	  MSDy += Sqr(OprogStatus.sumdy[i]);
	  MSDz += Sqr(OprogStatus.sumdz[i]);
	}
      MSDx /= (Oparams.parnum-Oparams.parnumA);
      MSDy /= (Oparams.parnum-Oparams.parnumA);
      MSDz /= (Oparams.parnum-Oparams.parnumA);
#ifdef MD_BIG_DT
      fprintf(f, "%.15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, MSDx, MSDy, MSDz);
#else
      fprintf(f, "%.15G %.15G %.15G %.15G\n", Oparams.time, MSDx, MSDy, MSDz);
#endif
      fclose(f);
    }
#endif
}
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
/* ============================ >>> temperat <<< =========================== */
#ifdef EDHE_FLEX
extern int get_dof_flex(int filter);
extern void calc_energy_filtered(int filter);
#endif
void temperat(void)
{
  double dof;
#ifdef EDHE_FLEX
  int kk;
  double P[3];
#endif
#if defined(MD_INELASTIC) || defined(MD_FOUR_BEADS)
  double dofTra, dofRot, tempRot, tempTra;
#endif
  /* DESCRIPTION:
     This the calculation of the instantaneous temperature */
#if 0
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
#endif
#ifdef EDHE_FLEX
  calc_energy_filtered(0);
#ifdef MD_FOUR_BEADS
  /* note that also angular momentum is conserved, hence we have 6 degrees of freedom less! */
  dof = ((double)Oparams.parnum)*6.0;
#else
  dof = get_dof_flex(0) - OprogStatus.frozenDOF;
#endif
#else
  calc_energy(NULL);
  dof = OprogStatus.dofA*((double)Oparams.parnumA) + 
    OprogStatus.dofB*((double) (Oparams.parnum-Oparams.parnumA));
#endif
  if (OprogStatus.brownian==1)
    temp = 2.0 * K / dof;
  else
#ifdef EDHE_FLEX
#ifdef MD_FOUR_BEADS
    temp = 2.0 * K / (dof-6.0);
#else
    temp = 2.0 * K / dof;
#endif
#else
    temp = 2.0 * K / (dof - 3.0);
#endif
  
  if (OprogStatus.avngTemp == 1)
    {
      OprogStatus.sumTemp += temp;
      temp = OprogStatus.sumTemp / NUMCALCS;
    }
#ifdef MD_FOUR_BEADS
  dofTra = 3.0*((double)Oparams.parnum);
  dofRot = dof - dofTra;
   if (OprogStatus.brownian==1)
    {
      tempRot = 2.0 * Krot / dofRot;
      tempTra = 2.0 * Ktra / dofTra;
    }
  else
    {
      tempRot = 2.0 * Ktra / dofRot;
      tempTra = 2.0 * Krot / (dofTra-6.0);
    }
#elif defined(MD_INELASTIC)
  dofTra = 3*((double)Oparams.parnum);
  dofRot = 2*((double)Oparams.parnum);
  if (OprogStatus.brownian==1)
    {
      tempRot = 2.0 * Krot / dofRot;
      tempTra = 2.0 * Ktra / dofTra;
    }
  else
    {
      tempRot = 2.0 * Ktra / dofRot;
      tempTra = 2.0 * Krot / (dofTra-3.0);
    }
#endif
  mf = fopenMPI(absMisHD("temp.dat"),"a");
#ifdef MD_INELASTIC
  mf2 =fopenMPI(absMisHD("temp_granular.dat"), "a"); 
#endif
#ifndef MD_FOUR_BEADS
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G\n", Oparams.time + OprogStatus.refTime, temp);
#else
  fprintf(mf, "%15G %.15G\n", Oparams.time, temp);
#endif
#endif
#ifdef MD_FOUR_BEADS
#ifdef MD_BIG_DT
  fprintf(mf, "%15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, temp, tempTra, tempRot);
#else
  fprintf(mf, "%15G %.15G %.15G %.15G\n", Oparams.time, temp, tempTra, tempRot);
#endif
#else
#ifdef MD_INELASTIC
#ifdef MD_BIG_DT
  fprintf(mf2, "%15G %.15G %.15G %.15G\n", Oparams.time + OprogStatus.refTime, ((double)OprogStatus.collCount), tempTra, tempRot);
#else
  fprintf(mf2, "%15G %.15G %.15G %.15G\n", Oparams.time, ((double)OprogStatus.collCount), tempTra, tempRot);
#endif
#endif
#endif
  fclose(mf);
#ifdef MD_INELASTIC
  fclose(mf2);
#endif
  /* pressure */
  if (OprogStatus.avngPress == 1)
    {
      OprogStatus.sumPress += press;
      press = OprogStatus.sumPress / NUMCALCS;
    }
#ifdef EDHE_FLEX
  sprintf(TXT, "DOF:%.15G T:%.15G", dof, temp);
  mdMsg(STD,NOSYS, NULL, "NOTICE", NULL,  TXT, NULL);
#endif
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
#ifdef MD_LXYZ
  Vol = L[0]*L[1]*L[2];  
#else
#if defined(MD_GRAVITY)
  Vol = Lz*Sqr(L);
#else
  Vol = L*L*L;  
#endif
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
#ifdef MD_LXYZ
  double scalFact[3];
  int kk;
#else
  double scalFact;
#endif
  double invNmA;
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
#ifdef MD_LXYZ
  for (kk=0; kk < 3; kk++)
    invL[kk] = 1.0 / L[kk];
#else
  invL = 1.0 / L;
#endif
  //scalFact = (KEND - KBEG) / NUMK;
  
  pi2 = 2.0 * pi;
  kbeg = 0.0; //pi2 * invL;
#ifdef MD_LXYZ
  for (kk=0; kk < 3; kk++)
    scalFact[kk] = pi2 * invL[kk];
#else
  scalFact = pi2 * invL;
#endif
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
#ifdef MD_LXYZ
	      rCMk = kbeg + scalFact[0] * rx[i] * mesh[n][mp][0] + scalFact[1] * ry[i] * mesh[n][mp][1] + 
		 scalFact[2] * rz[i] * mesh[n][mp][2];
#else
	      rCMk = kbeg + scalFact * 
		(rx[i] * mesh[n][mp][0] + ry[i] * mesh[n][mp][1] + 
		 rz[i] * mesh[n][mp][2]);
#endif
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
#ifdef MD_LXYZ
extern double min3(double a, double b, double c);
#endif
void radDens(void)
{
  int bin, Nm, i, j;
  int* hist;
  double Vol, rij, rxij, ryij, rzij, rijSq;
  double m, rlower, rupper, cost, nIdeal; 
  double rhoAv, DELR;

#ifdef MD_LXYZ
  for (i=0; i < 3; i++)
    invL[i] = 1.0 / L[i];
#else
  invL = 1.0 / L;
#endif
  Nm = Oparams.parnumA;
  rhoAv = (double) Nm;
  m = Oparams.m[0];
#ifdef MD_LXYZ
  Vol = L[0]*L[1]*L[2];
#else
#ifdef MD_GRAVITY
  Vol = Sqr(L)*Lz;
#else
  Vol = L*L*L;
#endif
#endif 

#ifdef MD_LXYZ
  DELR = (min3(L[0],L[1],L[2])/2.0 - RBEG) / MAXBIN;
#else
  DELR = (L/2.0 - RBEG) / MAXBIN;
#endif
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
#ifdef MD_LXYZ
	  rxij = rxij - L[0] * rint(invL[0] * rxij);
	  ryij = ryij - L[1] * rint(invL[1] * ryij);
	  rzij = rzij - L[2] * rint(invL[2] * rzij);
#else
	  rxij = rxij - L * rint(invL * rxij);
	  ryij = ryij - L * rint(invL * ryij);
#if !defined(MD_GRAVITY)
	  rzij = rzij - L * rint(invL * rzij);
#endif
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
  Ptens[0] = Pxy;
  Ptens[1] = Pyz;
  Ptens[2] = Pzx;

}

/* ========================== >>> DQtensor <<< ============================= */
void DQtensor(void)
{
  DQtens[0] = OprogStatus.DQxy;
  DQtens[1] = OprogStatus.DQyz;
  DQtens[2] = OprogStatus.DQzx;
}


