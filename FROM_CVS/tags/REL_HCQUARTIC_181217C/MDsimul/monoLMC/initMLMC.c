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
#define NUM_SPOST 6

extern int ENDSIM;
extern char msgStrA[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, Wxx, Wyy, Wzz,  Wxy, Wyz, Wzx;  

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax, **nebrTabi, **nebrTaba;

int cubeSize;

int *pmap, *pabs, *pSqr, *pgetb, *pgetj ;
double *wRadLT, *vRadLT, *sumCos, *sumSin, *DsumCos, *DsumSin;

extern int *rxS, *ryS, *rzS;
extern int maxdelta;

double **dcosLU, **dsinLU, *sinLU, *cosLU, *sumCosLU, *sumSinLU, *sumCosIni, *sumSinIni;
int *kcx, *kcy, *kcz;
int NNC, KK;
double *vx, *vy, *vz;

extern int minimumImageI(int);
extern int nebrSkNow;

/* ================================= */

extern inline int minimumImageI(int pl);

/* ============================= >>> FCC <<< ================================*/
void FCC(COORD_TYPE m)
{
  /*   DESCRIPTION:
       Sets up the alpha fcc lattice for n linear molecules.   
       The simulation box is a unit cube centred at the origin.
       N should be an integer of the form ( 4 * ( Nc ** 3 ) ),
       Where Nc is the number of FCC unit cells in each direction.  
       See figure 5.10 for a diagram of the lattice and a           
       definition of the four orientational sublattices.            
       PRINCIPAL VARIABLES:                                         
       COORD_TYPE    rxCm, ryCm, rzCm     Molecular Center of mass 
                                          positions             
       COORD_TYPE    rRoot3               1.0 / sqrt ( 3.0 ) */
  int Nc, Nm;
  COORD_TYPE rRoot3; // = 0.5773503;
  COORD_TYPE  Cell, Cell2, rxCm, ryCm, rzCm;
  int np,i, ix, iy, iz, iref, ii;
  COORD_TYPE bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  double la;

  //printf("FCC Vol: %f\n", Vol);
  L = cbrt(Vol);
  Nm = Oparams.parnum;
  Nc = ceil(  pow( ((COORD_TYPE)Nm)/4.0, 1.0/3.0 )  );
  //printf("Nc: %d\n", Nc);
  /* Calculate the side of the unit cell */
  Cell  = L / ((COORD_TYPE) Nc); /* unit cell length */
  Cell2 = 0.5 * Cell;              /* half unit cell length */

  /* Sublattice A */
  rRoot3 = 1.0 / sqrt(3.0);
  bx[0] =  0.0;
  by[0] =  0.0;
  bz[0] =  0.0;
  /*  Sublattice B */
  bx[1] =  Cell2;
  by[1] =  Cell2;
  bz[1] =  0.0;
  /* Sublattice C */
  bx[2] =  0.0;
  by[2] =  Cell2;
  bz[2] =  Cell2;
  /* Sublattice D */
  bx[3] =  Cell2;
  by[3] =  0.0;
  bz[3] =  Cell2;
  /* Construct the lattice from the unit cell */
  
  ii = 0;
  np = 0;
  la = Oparams.lattice_a;
  
  //printf("rx[1][50]: %f\n", rx[1][39]);
  //printf("NA: %d NB:%d\n", Oparams.parnum[0], Oparams.parnum[1]);
  //printf("modA: %d modB: %d mod: %d\n", modA, modB, mod);
  loop(iz, 1, Nc) /* loops over unit cells (that are simply cubes) */ 
    {
      loop(iy, 1, Nc)
	{
	  loop(ix, 1, Nc)
	    {
	      loop(iref, 1, 4) /* In each primitive cell there are four 
				  molecules */
		{

		  /* Questo vuol dire che ho piazzato tutte le particelle */
		  if (np >= Oparams.parnum) 
		      
		    break;
		  
		  /* nuova possibile posizione */
		  rxCm = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ryCm = by[iref] + Cell * ((COORD_TYPE) iy);
		  rzCm = bz[iref] + Cell * ((COORD_TYPE) iz);
		  //printf("CM:(%f,%f,%f)\n", rxCm, ryCm, rzCm);
		  
		  //printf("np[%d]=%d\n ", a, np[a]);
		  rx[np] = ((int)rint(rxCm/la));
		  ry[np] = ((int)rint(ryCm/la));
		  rz[np] = ((int)rint(rzCm/la));
		  //printf("(%d,%d,%d|%.5f)\n", rx[a][np[a]], 
		  // ry[a][np[a]], rz[a][np[a]], la);
		    
		  np++;
		}
	      
	      ii = ii + 4;

	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */


  for ( i = 0; i < Oparams.parnum; i++)
    {
      /* Initial position values are between -0.5 and 0.5 */
      rx[i] = minimumImageI(rx[i]); 
      ry[i] = minimumImageI(ry[i]);
      rz[i] = minimumImageI(rz[i]);
      //printf("(%d,%d,%d)\n", rx[a][i], 
      //  	 ry[a][i], rz[a][i]);
      
    }
    
  return;
}

/* ============================ >>> ranf <<< =============================== */
COORD_TYPE ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return rand() / ( (COORD_TYPE) RAND_MAX );
}

/* ============================= >>> gauss <<< ============================= */
COORD_TYPE gauss(void)
{
  
  /* 
     Random variate from the standard normal distribution.
     
     The distribution is gaussian with zero mean and unit variance.
     REFERENCE:                                                    
                                                                
     Knuth D, The art of computer programming, (2nd edition        
     Addison-Wesley), 1978                                      
                                                                
     ROUTINE REFERENCED:                                           
                                                                
     COORD_TYPE ranf()                                  
     Returns a uniform random variate on the range zero to one  
  */

  COORD_TYPE  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  COORD_TYPE sum, r, r2;
  int i;

  sum = 0.0;

  loop(i, 1, 12)
    {
      sum = sum + ranf();
    }

  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}

/* ========================= >>> resetCM <<< ==============================*/
void resetCM()
{
  COORD_TYPE RCMx, RCMy, RCMz;
  COORD_TYPE m, MTOT;
  double la;
  int i;

  m = Oparams.m;
  /* Remove net momentum, to have a total momentum equals to zero */
  la = Oparams.lattice_a;

  MTOT = 0.0;

  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  for(i = 0; i < Oparams.parnum; i++)
    {
      RCMx += m*la*((double)rx[i]); 
      /*Here RCM is the center of mass of the box */
      RCMy += m*la*((double)ry[i]);
      RCMz += m*la*((double)rz[i]);
      MTOT += m;
    }

  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;

  for(i = 0; i < Oparams.parnum; i++)
    {
      rx[i] -= rint(RCMx/la);
      ry[i] -= rint(RCMy/la);
      rz[i] -= rint(RCMz/la);
      rxS[i] = minimumImageI(rx[i]);
      ryS[i] = minimumImageI(ry[i]);
      rzS[i] = minimumImageI(rz[i]);
    }
    
}

/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  setToZero(Oparams.parnum, SAVE_LIST, 
	    NULL);  /* Set to zero all the coordinates */
  
  FCC(Oparams.m); 
  
  /* Put the baricenter of each molecule on a FCC lattice, and set 
     their orientations */  
  
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
  //angvel(Oparams.parnum, Oparams.T, Oparams.m, Oparams.d); 
}

/* =========================== >>> usrInitBef <<< ========================== */
void usrInitBef(void)
{
  /* DESCRIPTION:
     This function is called before any other initialization, put here 
     yours, for example initialize accumulators ! 
     NOTE: You should supply parameters value in parameters file, so this 
           initilization are quite fictitiuos for parameters, anyway 
	   accumulators initialization is crucial */
  
  /* ===================== >>> INIT PARAMETERS <<< ======================== 
   All the values set here for Oparams structure are taken as defaults if you
   don't specify corresponding parameters in the parameters file  */

  Dtrans = 0.0; /* DtransOld should become a field of OprogStatus */

  Vol = 1400.0;
  Vol1 = 0.0;
  Vol2 = 0.0;
  Vol1o1 = 0.0;
  Vol1o2 = 0.0;

  Oparams.T = 2.0;
  Oparams.P = 1.0;

  OprogStatus.avVol = 0.0;
  OprogStatus.tolVol = 0.0001;
  OprogStatus.tolVol1 = 0.1;
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  OprogStatus.avngS     = 0;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
  OprogStatus.bakStepsAscii = 0;
  OprogStatus.rateCheck = 10000;
  OprogStatus.mosseAccettate = 0;
  OprogStatus.mosseTentate = 0;
  Oparams.m = 1.0;

  OprogStatus.alpha = 0.83;
  OprogStatus.kmax = 7.12;
  OprogStatus.Dk = 0.34;
  OprogStatus.S0 = 10.0;

  Oparams.sigma = 1.0;
  Oparams.epsilon = 1.0;
  
  V = 0.0; /* potential energy */
  W = 0.0; /* virial function */
  VAC = 0.0;
  VLJ = 0.0;
  WAC = 0.0;
  WLJ = 0.0;

  srand((int)time(NULL));
  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  /* ======================================================================= */
}


extern void calcConst(void);
extern double rcut, rcutSq, dvdr; 
extern double sigmaSq, epsilon4, epsilon24;

/* ========================= >>> sumup <<< ================================ */
void sumupLJ(void)
{
  /* calcola i valori delle energie potenziali all'inizio della
     simulazione */
  int i, j, Nm, rxi, ryi, rzi, rxijI, ryijI, rzijI;
  double vij, wij, rSq, sr2, sr6, sr12, rxij, ryij, rzij, la;
  double vCut;
  int ncut;

  Nm = Oparams.parnum;

  calcConst();
  VLJ = 0.0;
  WLJ = 0.0;
  la = Oparams.lattice_a;
  printf("rcutSq: %f\n", rcutSq);
  ncut = 0;
    
  for (i = 0; i < Oparams.parnum; i++)
    {
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i];
      for (j = 0; j < Oparams.parnum; j++ )
	{
		  
	  if (j > i) 
	    {
	      rxijI = minimumImageI(rxi - rx[j]);
	      ryijI = minimumImageI(ryi - ry[j]);
	      rzijI = minimumImageI(rzi - rz[j]);
	      /* Moltiplica le coordinate intere per il passo 
		 del reticolo
		 in modo da ottenere le coordinate 
		 in unità ridotte */
	      rxij = la * ((double) rxijI);
	      ryij = la * ((double) ryijI);
	      rzij = la * ((double) rzijI);
	      
	      rSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	      
	      if ( rSq < rcutSq )
		/* 'rcut' is the cutoff for V */
		{
		  sr2   = sigmaSq / rSq;
		  sr6   = sr2 * sr2 * sr2;
		  sr12  = Sqr(sr6);
		  
		  vij     = sr12 - sr6;
		  //vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);
		  wij     = vij + sr12;
		  
		  VLJ = VLJ + vij;
		  WLJ = WLJ + wij; 
		  ++ncut;
		}
	    }
	}
    }

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */

  sr2 = sigmaSq / rcutSq;
  sr6 = sr2 * sr2 * sr2;
  sr12 = sr6 * sr6;
  vCut = sr12 - sr6;
  printf("ncut: %d vCut: %f ncut*vCut: %f\n", ncut, vCut, ((double)ncut)*vCut);
  VLJ = VLJ - ((double)ncut) * vCut;
  VLJ *= epsilon4;
  WLJ *= epsilon24 / 3.0;
  printf("epsilon4: %f INIZIO VLJ: %f rcutSq: %f\n", epsilon4, VLJ, rcutSq);
}

/* ========================== >>> rdiq <<< ============================= */
void rdiq(int *rx, int *ry, int *rz,  
	  double kx, double ky, double kz, double* Re_rh, double* Im_rh)
{
  double la, arg;//, sqrtNm;
  int Nm, i;
  
  *Re_rh = 0.0;
  *Im_rh = 0.0;
  Nm = Oparams.parnum;
  //sqrtNm = sqrt((C_T)Nm);

  la = Oparams.lattice_a;

  for(i = 0; i< Nm; i++)
    {
      arg = ((double)rx[i])*kx + ((double)ry[i])*ky + ((double)rz[i])*kz;
      arg *= la;
      *Re_rh = *Re_rh + cos(arg);
      *Im_rh = *Im_rh + sin(arg);
    }
  
  *Re_rh /= (sqrt((double)Nm));
  *Im_rh /= (sqrt((double)Nm));

}

/* ========================= >>> sumup <<< ================================ */
void sumupAC(void)
{
  int Nm, ind, indMax, i;
  double kx, ky, kz, pi2, sks, cost1, invNm, cost2;
  double Re_rh, Im_rh, dd, invL;
  double alpha, Dk, kmax, S0;

  kmax = OprogStatus.kmax;
  alpha = OprogStatus.alpha;
  Dk = OprogStatus.Dk;
  S0 = OprogStatus.S0;

  pi2 = 2.0*pi;

  invL = 1.0 / cbrt(Vol);
  cost2 = pi2*invL;

  Nm = Oparams.parnum;
  invNm = 1.0 / ((double)Nm);
  //cost1 = -8.0 * alpha / sqrt((C_T)Nm); 
  // N.B. 
  // The AC potential is 1/2*Sum[(S(k)-S0)^2]
  cost1 = -2.0 * alpha / sqrt((double)Nm); 
  VAC = 0.0;
  WAC = 0.0;
  printf("QUI NNC: %d\n", NNC);
  for(ind = 0; ind < NNC; ind++)
    {
      kx = cost2 * kcx[ind];
      ky = cost2 * kcy[ind];
      kz = cost2 * kcz[ind];
      //printf("cost2:%f kcx[%d]:%d\n", cost2, ind, kcx[ind]);
      rdiq(rx, ry, rz, kx, ky, kz, &Re_rh, &Im_rh);
      /* questo 2.0 viene dalla simmetria rispetto al piano yz */
      sumCos[ind] = Re_rh;
      sumSin[ind] = Im_rh; 
      sks = Sqr(Re_rh) + Sqr(Im_rh);
     
      dd = sks - S0;
      if (dd > 0.0)
	{
	  VAC = VAC + Sqr(dd);
	}
    }
  VAC *= alpha;/*lo 0.5 non c'è più per via della simmetria nello spazio k*/
  WAC = 0.0; /* per il momento non calcola il viriale AC!!!!!!!!!!!!!!! */
  printf("VAC:%f\n", VAC);
  //EE = VAC / Lambda;
  //VAC *= alpha;
}

  
/* ========================== >>> AllocCubeR <<< ===========================*/
double*** AllocCubeR(int size1, int size2, int size3)
{
  double*** v;
  void* buffer;
  int k1, k2;
  
  /* alloca la memoria per il cubo */
  buffer = malloc(size1 * size2 * size3 * sizeof(double));

  v = (double***) malloc(size1 * sizeof(double**));
  for (k1 = 0; k1 < size2; k1++)
    v[k1] = (double**) malloc(size2 * sizeof(double*));

  for (k1 = 0; k1 < size1; k1++)
    {
      v[k1][0] = ((double*) buffer) + k1 * size2 * size3;
      for (k2 = 1; k2 < size2; k2++)
	{
	  v[k1][k2] = v[k1][k2-1] + size3;
	}
    }

  return v;
}

/* ========================== >>> AllocCubeR <<< ===========================*/
int*** AllocCubeI(int size1, int size2, int size3)
{
  int*** v;
  void* buffer;
  int k1, k2;
  
  /* alloca la memoria per il cubo */
  buffer = malloc(size1 * size2 * size3 * sizeof(int));

  v = (int***) malloc(size1 * sizeof(int**));
  for (k1 = 0; k1 < size2; k1++)
    v[k1] = (int**) malloc(size2 * sizeof(int*));

  for (k1 = 0; k1 < size1; k1++)
    {
      v[k1][0] = ((int*) buffer) + k1 * size2 * size3;
      for (k2 = 1; k2 < size2; k2++)
	{
	  v[k1][k2] = v[k1][k2-1] + size3;
	}
    }

  return v;
}

int csSq;

/* ========================= >>> buildEnergiesLJ <<< ======================= */
void buildEnergiesLJ()
{
  int i; //csSq;
  double rSq, sr2, sr6, sr12;
  double vCut;
  double cutoff, cutoffSq;

  calcConst();

  /* MODIFICA 25/04/01:
     controllare qui perché prima c'era un baco enorme sul calcolo
     dei cutoff interi che dava delle energie potenziali differenti */ 
  cutoff = (Oparams.sigma * Oparams.rcut)/Oparams.lattice_a;
  cubeSize = (int)cutoff;
  if (((double)cubeSize) < cutoff )
    cubeSize += 1;
  
  cutoffSq = Sqr(cutoff);
  //csSq = Sqr(cubeSize);
  csSq = (int) cutoffSq;
  
  if (((double)csSq) < cutoffSq )
    csSq += 1;
	  
  printf("sigmaSq: %f espilon4: %f rcut*sigma: %.8f sigma: %.8f cutoff: %.6f\n",
	 sigmaSq, epsilon4, Oparams.sigma * Oparams.rcut, Oparams.sigma, 
	 cutoff);
  
  printf("la: %.10f csSq:%d\n",Oparams.lattice_a, 
	 csSq);
  
  vRadLT = malloc(sizeof(double)*csSq);
  wRadLT = malloc(sizeof(double)*csSq);
	  
  sr2 = sigmaSq / rcutSq;
  sr6 = sr2 * sr2 * sr2;
  sr12 = sr6 * sr6;
  vCut = sr12 - sr6;
  
  for (i = 0; i < csSq; i++)
    {
      
      rSq = Sqr((double)Oparams.lattice_a) * ((double)i);
      
      if ( rSq < rcutSq )/* 'rcut' is the cutoff for V */
	{
	  
	  sr2   = sigmaSq / rSq;
	  sr6   = sr2 * sr2 * sr2;
	  sr12  = Sqr(sr6);
	  vRadLT[i] = epsilon4*(sr12 - sr6 - vCut);
	  //wRadLT[i] = epsilon24*
	  // ((sr12 - sr6) + sr12)/3.0;
	}
    }
}


/* ========================= >>> buildEnergiesAC <<< ======================= */
void buildEnergiesAC()
{
  int Nm, ind, rix, riy, riz;
  double pi2, cost1, invNm, cost2;
  double invL, maxdcos, maxdsin;
  double fact, kpDkSq, kmDkSq, kSq;
  double radi2, *buf, S0, kmax, Dk;
  int numirk, kri, ikx, iky, ikz;
  int i, kk, lM;
  
  pi2 = 2.0*pi;
  invL = 1.0 / cbrt(Vol);
  cost2 = pi2*invL;

  Nm = Oparams.parnum;
  invNm = 1.0 / ((double)Nm);
  cost1 = -2.0 * OprogStatus.alpha / sqrt((double)Nm); 
  
  kmax = OprogStatus.kmax;
  S0 = OprogStatus.S0;
  Dk = OprogStatus.Dk;

  VAC = 0.0;
  WAC = 0.0;
  
  /* ======================== >>> buildMesh <<< ==================== */
  L = cbrt(Vol);
  fact = L/(2.0*pi);
  kpDkSq = Sqr((kmax+Dk)*fact);
  kmDkSq = Sqr((kmax-Dk)*fact);
  
  KK = ((int)sqrt(kpDkSq)) + 1;

  ind = 0;
  
  /* scopre quanti k ci sono */
  for (ikx = 0; ikx <= KK; ikx++)
    for(iky = -KK; iky <= KK; iky++)
      for(ikz = -KK; ikz <= KK; ikz++)
	{   
	  /* prende la metà dei k sfruttando la simmetria rispetto al piano 
	     yz */
	  if ((ikx == 0) && (iky < 0))
	    continue;
	  if ((ikx == 0) && (iky == 0) && (ikz < 0))
	    continue;

	  kSq = Sqr(ikx)+Sqr(iky)+Sqr(ikz);
	  if((kSq < kpDkSq) &&  (kSq > kmDkSq))  
	    {
	      ind++;
	    }
	}
  
  NNC = ind;

  sumCos = malloc(NNC*sizeof(double));
  sumSin = malloc(NNC*sizeof(double));
  sumCosIni = malloc(NNC*sizeof(double));
  sumSinIni = malloc(NNC*sizeof(double));

  DsumCos = malloc(NNC*sizeof(double));
  DsumSin = malloc(NNC*sizeof(double));
  
  kcx = malloc(NNC*sizeof(int));
  kcy = malloc(NNC*sizeof(int));
  kcz = malloc(NNC*sizeof(int));

  ind = 0;
  //printf("antiCry: %d NNC: %d\n", OprogStatus.AntiCry, NNC);
  
  for (ikx = 0; ikx <= KK; ikx++)
    for(iky = -KK; iky <= KK; iky++)
      for(ikz = -KK; ikz <= KK; ikz++)
	{ 
	  /* prende la metà dei k sfruttando la simmetria rispetto al piano 
	     yz */
	  if ((ikx == 0) && (iky < 0))
	    continue;
	  if ((ikx == 0) && (iky == 0) && (ikz < 0))
	    continue;

	  kSq = Sqr(ikx)+Sqr(iky)+Sqr(ikz);
	  if((kSq < kpDkSq) &&  (kSq > kmDkSq))  
	    {
	      kcx[ind]=ikx;
	      kcy[ind]=iky;
	      kcz[ind]=ikz;
	      //printf("(%d,%d,%d) |", kcx[ind], kcy[ind], kcz[ind]);
	      ind++;
	    }
	}

  lM = Oparams.lattice_M;
  numirk = 3*KK*lM;
  printf("NNC: %d lattice_M: %d numrik: %d\n", NNC, Oparams.lattice_M, numirk);
  buf = malloc(sizeof(double)*(2*numirk+1));
  cosLU = buf + numirk;
  buf = malloc(sizeof(double)*(2*numirk+1));
  sinLU = buf + numirk;
  /* ricordarsi si shiftare il primo indice di KK */
  dcosLU = (double**) malloc(sizeof(double*)*(2*KK+1));
  dsinLU = (double**) malloc(sizeof(double*)*(2*KK+1));

  dcosLU[0] = (double*) malloc(sizeof(double)*(2*KK+1)*(2*numirk+1));
  dsinLU[0] = (double*) malloc(sizeof(double)*(2*KK+1)*(2*numirk+1));
 
  for (i = -KK + 1; i <= KK; i++)
    {
      dcosLU[i+KK] = dcosLU[i+KK-1] + (2*numirk+1);
      dsinLU[i+KK] = dsinLU[i+KK-1] + (2*numirk+1);
    }

  printf("ind up to %d\n", NNC);
  
  //radi2 = sqrt(2.0);
  for (kri = -numirk; kri <= numirk; kri++)
    {
      cosLU[kri] = cos(kri*Oparams.lattice_a*cost2)/(sqrt((double)Nm));
      sinLU[kri] = sin(kri*Oparams.lattice_a*cost2)/(sqrt((double)Nm));
    }

  maxdcos = 0.0;
  maxdsin = 0.0;
  for (kri = -numirk; kri <= numirk; kri++)
    {
      for (kk = -KK; kk <= KK; kk++)
	{
	  if (abs(kri + kk) < numirk)
	    {
	      dcosLU[kk+KK][kri] = cosLU[kri+kk] - cosLU[kri];
	      dsinLU[kk+KK][kri] = sinLU[kri+kk] - sinLU[kri];
	      if (dcosLU[kk+KK][kri] > maxdcos)
		maxdcos = dcosLU[kk+KK][kri];
	      if (dsinLU[kk+KK][kri] > maxdsin)
		maxdsin = dsinLU[kk+KK][kri];
	    }
	}
    }
  printf("maxdcos: %f maxdsin: %f\n", maxdcos, maxdsin);
  /* In questo modo abbiamo tutti i seni e coseni che ci sevono senza doverli calcolare
   */

}

/* ========================== >>> builgPMap <<< ========================= */
void buildMaps(void)
{
  int i, lM, Nm;
  int *buf;
  lM = Oparams.lattice_M;
  Nm = Oparams.parnum;
  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pmap = buf + lM;

  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pabs = buf + lM;

  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pSqr = buf + lM;

  for (i = -lM; i <= lM; i++)
    {
      pmap[i] = minimumImageI(i);
      pabs[i] = abs(i); 
      pSqr[i] = Sqr(i);
     
    }
}

extern int *nebrSkTab, *NnST, *NOTnebrSkTab, *nST;

/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */

  int Nm, i, sct;
  COORD_TYPE m;

  //COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;

  /* initialize global varibales */
  pi = 2.0 * acos(0);
  
  mdMsg(ALL, NOSYS, "usrInitAft", "NOTICE", NULL,
	"Using Neighbour List method",
	NULL);
    
  /* Allocate arrays of integers, you must supply the array length and the
     address of the pointer to the array */
  Nm = Oparams.parnum;

  sct = sizeof(COORD_TYPE);
  
  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac;
  nebrNow = 1;
  nebrSkNow = 1;
  nebrTab = AllocMatI(Nm, nebrTabMax); 
  nebrTabi = AllocMatI(Nm, nebrTabMax); 
  nebrTaba = AllocMatI(Nm, nebrTabMax); 
 
  nebrTabLen = malloc(sizeof(int)*Nm);

  vx = malloc(sizeof(double)*Nm);
  vy = malloc(sizeof(double)*Nm);
  vz = malloc(sizeof(double)*Nm);
 
  /* =============================================== */
  
  for (i=0; i < Nm; i++)
    {
      nebrTabLen[i] = 0;
    }
  
  Oparams.lattice_a = cbrt(Vol) / Oparams.lattice_M;
  //  printf("lattice_a: %.5f Vol: %.5f\n", Oparams.lattice_a, Vol);
  /* Store the Center of Mass initial position for all particles */
  m = Oparams.m;

  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      for(i = 0; i < Oparams.parnum; i++)
	{
	  
	  /* store the initial positions of particles */
	  OprogStatus.rxi[i] = rx[i];
	  OprogStatus.ryi[i] = ry[i];
	  OprogStatus.rzi[i] = rz[i];
	}
	
      
      loop(i, 1, NUMK) 
	{
	  OprogStatus.sumS[i] = 0.0;
	}
      
    }
  printf("Vol: %.15f\n", Vol);
  
  
  /* coordinate scalate alla prima scatola */
  rxS = (int *)malloc(sizeof(int)*Oparams.parnum);
  ryS = (int *)malloc(sizeof(int)*Oparams.parnum);
  rzS = (int *)malloc(sizeof(int)*Oparams.parnum);
    

  for (i = 0; i < Oparams.parnum; i++)
    {
      /* coordinate scalate alla prima scatola */
      rxS[i] = minimumImageI(rx[i]);
      ryS[i] = minimumImageI(ry[i]);
      rzS[i] = minimumImageI(rz[i]);
    }

  /* calcola tutte le energie di interazione all'inizio */

  buildEnergiesLJ();
  
  buildEnergiesAC();
  
  buildMaps();
  
  nebrSkTab = malloc(sizeof(int)*NNC);
  NOTnebrSkTab = malloc(sizeof(int)*NNC);
  NnST = malloc(sizeof(int)*NNC);
  nST = malloc(sizeof(int)*NNC);
  
  /* calcola le energie potenziali e il viriale all'inizio della simulazione */
  sumupLJ();
  sumupAC();
  V = VAC + VLJ;
  printf("VAC:%f VLJ:%f V:%f\n", VAC, VLJ, V);
  W = WAC + WLJ;


}



/* ========================== >>> comvel <<< =============================== */
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE m)
{
  COORD_TYPE rTemp, sumx, sumy, sumz, RCMx, RCMy, RCMz;
  //COORD_TYPE Px, Py, Pz;
  int i;
  
  rTemp = sqrt(temp / Oparams.m);  
  /* variance of the velocities distribution function, we assume k = 1 */ 

  loop(i, 1, Nm)
    {
      /* Set the velocities of both atoms to the center of mass velocities,
         the exact velocities will be set in the angvel() routine, where we 
         will set:
	 Atom 1: v1  = Vcm + W^(d21 * m2/(m2+m1))
	 Atom 2: v2  = Vcm - W^(d21 * m1/(m1+m2))
	 where Vcm is the center of mass velocity (that is the actual 
	 velocity of both atoms), W is the angular velocity of the molecule,
	 d21 is the vector joining the two atoms (from 2 to 1) and 
	 m1 and m2 are the masses of two atoms 
      */
      
      vx[i] = rTemp * gauss(); 
      vy[i] = rTemp * gauss();
      vz[i] = rTemp * gauss();
      //printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
    }
  
    /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  loop(i, 1, Nm)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 sumx = sumx + vx[i];
	 sumy = sumy + vy[i];
	 sumz = sumz + vz[i];
	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  loop(i, 1, Nm)
    {
      vx[i] = vx[i] - sumx;
      vy[i] = vy[i] - sumy;
      vz[i] = vz[i] - sumz;
      //printf("ss=%d i=%d VVV(%.10f,%.10f,%.10f)\n", ss, i, 
      //     vx[ss][i], vy[ss][i], vz[ss][i]);
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }
}

/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i;
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", ((double)rx[i])*Oparams.lattice_a, 
	      ((double)ry[i])*Oparams.lattice_a, ((double)rz[i])*Oparams.lattice_a);
    }
  

  /* assegna le velocità in manierra gaussiana (Maxwell-Boltzman) in modo 
     che il file generato possa essere usato da una MD (monoLJDPT) */
  comvel(Oparams.parnum, Oparams.T, Oparams.m);

  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G\n", vx[i], vy[i], vz[i]);
    }

  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G %.15G",  Vol, 1.0, 0.0, 0.0, 0.0, 0.0);


}

/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i;
  double xx, yy, zz;

  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf\n", &xx, &yy, &zz) < 3)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
      rx[i] = rint(xx/Oparams.lattice_a);
      ry[i] = rint(yy/Oparams.lattice_a);
      rz[i] = rint(zz/Oparams.lattice_a);
    }
  
}
