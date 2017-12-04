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

double ***vabLT[NA][NA];
double ***wabLT[NA][NA];

int cubeSize[NA][NA];

int *pmap, *pabs, *pSqr, *pgetb, *pgetj ;
double *wabRadLT[NA][NA], *vabRadLT[NA][NA];
double **DwabRadLT[NA][NA], **DvabRadLT[NA][NA];

extern int *rxS[NA], *ryS[NA], *rzS[NA];
extern int maxdelta[NA][NA];

extern int minimumImageI(int);

/* ================================= */

extern inline int minimumImageI(int pl);

/* ============================= >>> FCC <<< ================================*/
void FCC(COORD_TYPE* m)
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
  int modA, modB, mod;
  COORD_TYPE rRoot3; // = 0.5773503;
  COORD_TYPE  Cell, Cell2, rxCm, ryCm, rzCm;
  int np[NA], a, i, ix, iy, iz, iref, ii;
  COORD_TYPE bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  double la;

  //printf("FCC Vol: %f\n", Vol);
  L = cbrt(Vol);
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
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
  np[0] = 0;
  np[1] = 0;

  la = Oparams.lattice_a;

  if (Oparams.parnum[0] > Oparams.parnum[1])
    {
      modA = rint(Oparams.parnum[0] / Oparams.parnum[1]);
      modB = 1;
      mod  = modA + modB; 
      /* Piazzera' modA atomi A(0) e poi 1 atomo B(1) */
    }
  else
    {
      modB = rint(Oparams.parnum[1] / Oparams.parnum[0]);
      modA = 1;
      mod  = modA + modB;
      /* Piazzera' modB atomi B(1) e poi 1 atomo A(0) */
    }

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
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
		  if ( (np[0] >= Oparams.parnum[0]) &&
		       (np[1] >= Oparams.parnum[1]) )
		    break;
		  
		  /* nuova possibile posizione */
		  rxCm = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ryCm = by[iref] + Cell * ((COORD_TYPE) iy);
		  rzCm = bz[iref] + Cell * ((COORD_TYPE) iz);
		  //printf("CM:(%f,%f,%f)\n", rxCm, ryCm, rzCm);
		  
		  if (modA == 1) 
		    {
		      /* questa condizione e' vera una volta ogni mod volte
			 in questo modo mette un atomo A e poi 
			 modB atomi B */
		      if ( (ii + iref) % mod == 0 )
			{
			  a = 0;
			}
		      else
			{
			  a = 1;
			}
		    }
		  else//qui se modB == 1
		    {
		      /* questa condizione e' vera una volta ogni mod volte
			 in questo modo mette un atomo B e poi 
			 modA atomi A */
		      if ( (ii + iref) % mod == 0 )
			{
			  a = 1;
			}
		      else
			{
			  
			  a = 0;
			}
		    }
		  /* se gli atomi di tipo a (0 o 1) sono stati tutti 
		     piazzati allora mette quelli dell'altro tipo */
		  if (np[a] == Oparams.parnum[a])
		    {
		      a = (~a) & 1;
		    }
		  //printf("np[%d]=%d\n ", a, np[a]);
		  rx[a][np[a]] = ((int)rint(rxCm/la));
		  ry[a][np[a]] = ((int)rint(ryCm/la));
		  rz[a][np[a]] = ((int)rint(rzCm/la));
		  //printf("(%d,%d,%d|%.5f)\n", rx[a][np[a]], 
		  // ry[a][np[a]], rz[a][np[a]], la);
		    
		  np[a]++;
		}
	      ii = ii + 4;
	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */

  for ( a = 0; a < NA; a++)
    {
      for ( i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* Initial position values are between -0.5 and 0.5 */
	  rx[a][i] = minimumImageI(rx[a][i]); 
	  ry[a][i] = minimumImageI(ry[a][i]);
	  rz[a][i] = minimumImageI(rz[a][i]);
	  //printf("(%d,%d,%d)\n", rx[a][i], 
	  //  	 ry[a][i], rz[a][i]);
		  
	}
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
  COORD_TYPE *m, MTOT;
  double la;
  int i, a;

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

  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  RCMx += m[a]*la*((double)rx[a][i]); 
	  /*Here RCM is the center of mass of the box */
	  RCMy += m[a]*la*((double)ry[a][i]);
	  RCMz += m[a]*la*((double)rz[a][i]);
	  MTOT += m[a];
	}
    }
  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;

  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  rx[a][i] -= rint(RCMx/la);
	  ry[a][i] -= rint(RCMy/la);
	  rz[a][i] -= rint(RCMz/la);
	  rxS[a][i] = minimumImageI(rx[a][i]);
	  ryS[a][i] = minimumImageI(ry[a][i]);
	  rzS[a][i] = minimumImageI(rz[a][i]);
	}
    }
}

/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  setToZero(Oparams.parnum[0], SAVE_LISTA, 
	    NULL);  /* Set to zero all the coordinates */
  
  setToZero(Oparams.parnum[1], SAVE_LISTB, 
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
  int i, a, b;
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
  OprogStatus.fstps = 1;	
  for (i=0; i<NA; ++i)
    {
      Oparams.m[i] = 1.0;
    }

  for(a=0; a<NA; ++a)
    {
      for(b=0; b<NA; ++b)
	{
	  Oparams.sigab[a][a] = 1.0;
	  Oparams.epsab[a][b] = 1.0;
	}
    }
  V = 0.0; /* potential energy */
  W = 0.0; /* virial function */

  srand((int)time(NULL));
  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  /* ======================================================================= */
}


extern void calcConst(void);
extern double rcutab[NA][NA], rcutabSq[NA][NA], dvdr[NA][NA]; 
extern double sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA];

/* ========================= >>> sumup <<< ================================*/
void sumup(void)
{
  /* calcola i valori delle energie potenziali all'inizio della
     simulazione */
  int aa, bb, a, i, ai, b, j, bj, Nm, rxa, rya, rza, rxabI, ryabI, rzabI;
  double vab, wab, rabSq, srab2, srab6, srab12, rxab, ryab, rzab, la;
  double Wab[NA][NA], Vab[NA][NA], ncut[NA][NA], Vcab[NA][NA], vabCut;
  Nm = Oparams.parnum[0]+Oparams.parnum[1];

  calcConst();
  V = 0.0;
  W = 0.0;
  la = Oparams.lattice_a;

  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  ncut[a][b] = 0;
	  Vab[a][b] = 0;
	  Wab[a][b] = 0;
	}
    }
  for(a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];

	  ai = i+Nm*a;
	  for(b = 0; b < NA; b++)
	    {
	      for (j = 0; j < Oparams.parnum[b]; j++ )
		{
		  bj = j+b*Nm;
		  if (bj > ai) 
		    {
		      rxabI = minimumImageI(rxa - rx[b][j]);
		      ryabI = minimumImageI(rya - ry[b][j]);
		      rzabI = minimumImageI(rza - rz[b][j]);
		      /* Moltiplica le coordinate intere per il passo 
			 del reticolo
			 in modo da ottenere le coordinate 
			 in unità ridotte */
		      rxab = la * ((double) rxabI);
		      ryab = la * ((double) ryabI);
		      rzab = la * ((double) rzabI);
      
		      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		      
		      if ( rabSq < rcutabSq[a][b] )
			/* 'rcut' is the cutoff for V */
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
			  /* Virial off-diagonal terms of atomic pressure 
			     tensor */
			  ++ncut[a][b];
			}
		    }
		}
	    }
	  
	}
    }

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
  
  for(aa = 0; aa < NA; aa++)
    {
      for(bb = 0; bb < NA; bb++) /* b >= a */
	{
	  srab2 = sigabSq[aa][bb] / rcutabSq[aa][bb];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;
	  Vcab[aa][bb] = Vab[aa][bb] - ncut[aa][bb] * vabCut;
	  /* ncut[a][b] is the number of atoms pairs a-b within 
	     rcutab[a][b] */ 
	}
    }

  /* Moltiplica per le costanti  */
  for (aa = 0; aa < NA; aa++)
    {
      for (bb = 0; bb < NA; bb++)
	{
	  V += Vcab[aa][bb] * epsab4[aa][bb];
	  W  += Wab[aa][bb]  * epsab24[aa][bb] / 3.0;
	}
    }

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

/* ======================= >>> calcolaEnergie <<< ======================== */
void calcolaCubi(int a, int b)
{
  int ix, iy, iz;
  double la, rxab, ryab, rzab;
  double rabSq, srab2, srab6, srab12;

  la = Oparams.lattice_a;
  printf("cubeSize[%d][%d]:%d\n", a, b, cubeSize[a][b]);
  printf("rcutabSq:%.10f\n", rcutabSq[a][b]);
  for (ix = 0; ix < cubeSize[a][b]; ix++)
    {
      for (iy = 0; iy < cubeSize[a][b]; iy++)
	{
	  for (iz = 0; iz < cubeSize[a][b]; iz++)
	    {


	      /* Moltiplica le coordinate intere per il passo del reticolo
		 in modooo da ottenere le coordinate in unità ridotte */
	      rxab = la * ((double) ix);
	      ryab = la * ((double) iy);
	      rzab = la * ((double) iz);
	      
	      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
	     
	      if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		{
		  srab2   = sigabSq[a][b] / rabSq;
		  srab6   = srab2 * srab2 * srab2;
		  srab12  = Sqr(srab6);
		  //printf("(%d,%d,%d) 12-6:%.20f\n", ix, iy, iz,srab12 - srab6);
		  vabLT[a][b][ix][iy][iz] = srab12 - srab6;
		  //vab[a][b][ix][iy][iz] -=     
		  // dvdr[a][b] * (rab - rcutab[a][b]);
		  wabLT[a][b][ix][iy][iz] = (srab12 - srab6) + srab12;
		}
	      else
		{
		  vabLT[a][b][ix][iy][iz] = -1E100;
		  wabLT[a][b][ix][iy][iz] = 0.0;
		}
	    }
	}
    }
  //printf(">>>>>>>>>>>>>>>>> 1-0 12 12 0: %.30f %.30f\n ", 
  // vabLT[1][0][12][12][0], wabLT[1][0][12][12][0]);
}

/* ======================== >>> buildEnergies <<< ==========================*/
void buildEnergies(void)
{
  int a, b;
   //for()
  
  calcConst();
  /* calcola tutte le energie per le varie possibili interazioni
     0-1 0-0 1-1 */
  /* ALLOCAZIONE */
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  cubeSize[a][b] = ((int)
	    ( (Oparams.sigab[a][b]*Oparams.rcut)/Oparams.lattice_a) )+1;
	  //printf("sigab: %.6f rcut: %.6f la: %.6f cubeSize[%d][%d]: %d\n", 
	  // Oparams.sigab[a][b], Oparams.rcut, Oparams.lattice_a, a, b, cubeSize[a][b]);
	  if ((a == 1) && (b == 0)) 
	    continue;
	    

	    vabLT[a][b] = 
	      AllocCubeR(cubeSize[a][b], cubeSize[a][b], cubeSize[a][b]);
	   
	    wabLT[a][b] = 
	      AllocCubeR(cubeSize[a][b], cubeSize[a][b], cubeSize[a][b]);
	   
	}
    }
  
  vabLT[1][0] = vabLT[0][1];
  wabLT[1][0] = wabLT[0][1];
	  
  /* CALCOLO */
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  if (a == 1 && b == 0)
	    continue;
	  calcolaCubi(a, b);
	}
    }
}

int csSq[NA][NA];

/* ========================= >>> buildEnergiesR <<< ======================= */
void buildEnergiesRad()
{
  int i, a, b; //csSq;
  double rabSq, srab2, srab6, srab12;
  double vabCut;
  double cutoff, cutoffSq;

  calcConst();


  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  /* MODIFICA 25/04/01:
	     controllare qui perché prima c'era un baco enorme sul calcolo
	     dei cutoff interi che dava delle energie potenziali differenti */ 
	  cutoff = (Oparams.sigab[a][b]*Oparams.rcut)/Oparams.lattice_a;
	  cubeSize[a][b] = (int)cutoff;
	  if (((double)cubeSize[a][b]) < cutoff )
	    cubeSize[a][b] += 1;
	  
	  cutoffSq = Sqr(cutoff);
	  //csSq[a][b] = Sqr(cubeSize[a][b]);
	  csSq[a][b] = (int) cutoffSq;

	  if (((double)csSq[a][b]) < cutoffSq )
	    csSq[a][b] += 1;
	  
	  printf("rcut*sigabSq: %.8f sigab: %.8f cutoff: %.6f\n",
		 Oparams.sigab[a][b]*Oparams.rcut, Oparams.sigab[a][b], 
		 cutoff);
	  
	  printf("la: %.10f csSq[%d][%d]:%d\n",Oparams.lattice_a, a, b, 
		 csSq[a][b]);

	  vabRadLT[a][b] = malloc(sizeof(double)*csSq[a][b]);
	  wabRadLT[a][b] = malloc(sizeof(double)*csSq[a][b]);
	}
    }


  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  //csSq = Sqr(cubeSize[a][b]);
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;

	  for (i = 0; i < csSq[a][b]; i++)
	    {
	      
	      rabSq = Sqr((double)Oparams.lattice_a) * ((double)i);
	     
	      if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		{

		  srab2   = sigabSq[a][b] / rabSq;
		  srab6   = srab2 * srab2 * srab2;
		  srab12  = Sqr(srab6);

		  //printf("(%d) 12-6:%.20f\n", i,srab12 - srab6);
		  vabRadLT[a][b][i] = epsab4[a][b]*(srab12 - srab6 - vabCut);
		  //printf("8--() DENTRO!!!!%.10f\n", vabRadLT[a][b][i]);
		  //vab[a][b][ix][iy][iz] -=     
		  // dvdr[a][b] * (rab - rcutab[a][b]);
		  //wabRadLT[a][b][i] = epsab24[a][b]*
		  // ((srab12 - srab6) + srab12)/3.0;
		}
	    }
	}
    }
  
}

/* ====================== >>> buildDeltEnergiesR <<< ======================= */
void buildDeltEnergiesRadORIG()
{
  int i, a, b, csSq, d;
  double DrSq, vabCut, rabSq, srab2, srab6, srab12;

  calcConst();
  
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  cubeSize[a][b] = 
	    ((int) ( (Oparams.sigab[a][b]*Oparams.rcut)/Oparams.lattice_a) )+1
	    + 1;
	  /* 
	     questo +1 serve poiché vengono costruite tabelle do lookup
	     per differenze di energie rlativi a siti che differiscono 
	     al piu' per un passo reticolare
	   */
	  maxdelta[a][b] = 1+2*cubeSize[a][b];
	  csSq = Sqr(cubeSize[a][b]);
	  printf("devo allocare: %d bytes!!!!\n", (2*maxdelta[a][b]+1)*sizeof(double)*csSq);
	  //	  exit(-1);

	  DvabRadLT[a][b] = AllocMatR(2*maxdelta[a][b]+1, sizeof(double)*csSq);
	  DwabRadLT[a][b] = AllocMatR(2*maxdelta[a][b]+1, sizeof(double)*csSq);
	}
    }
    
  
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;

	  for (d = 0; d < 2*maxdelta[a][b]+1; d++)
	    {
	      
	      csSq = Sqr(cubeSize[a][b]);
	      
	      for (i = 0; i < csSq; i++)
		{
		  
		  rabSq = Sqr((double)Oparams.lattice_a) * ((double)i);
		  DrSq = Sqr((double)Oparams.lattice_a) 
		    * ((double)(d - maxdelta[a][b])); 
		  if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      
		      srab2   = sigabSq[a][b] / rabSq;
		      srab6   = srab2 * srab2 * srab2;
		      srab12  = Sqr(srab6);
		      /* Calcola il potenziale e il viriale includendo
			 i fattori moltiplicativi */
		      //printf("(%d) 12-6:%.20f\n", i,srab12 - srab6);
		      DvabRadLT[a][b][d][i] = epsab4[a][b]*
			((srab12 - srab6) - vabCut);
		      DwabRadLT[a][b][d][i] = 
			epsab24[a][b]*((srab12 - srab6) + srab12)/3.0;
		    }
		  else
		    {
		      DvabRadLT[a][b][d][i] = 0.0;
		      DwabRadLT[a][b][d][i] = 0.0;
		    }

		  rabSq += DrSq; 
		  if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      srab2   = sigabSq[a][b] / rabSq;
		      srab6   = srab2 * srab2 * srab2;
		      srab12  = Sqr(srab6);
		      
		      DvabRadLT[a][b][d][i] = epsab4[a][b]*
			((srab12 - srab6) - vabCut) - DvabRadLT[a][b][d][i];
		      DwabRadLT[a][b][d][i] = 
			epsab24[a][b]*((srab12 - srab6) + srab12)/3.0 -
			DwabRadLT[a][b][d][i];
		    }
		  if (i == 0)printf("[%d,%d]\n", i, d);
		  if (i == 324 && d == 140)
		  printf("la: %.6f epsab4: %.10f DrSq: %.10f [%d,%d,%d,%d]: %.10f \n ", 
			 Oparams.lattice_a, epsab4[a][b], DrSq,a, b, d, i, DvabRadLT[a][b][d][i]);
		}
	    }
	}
    }
  exit(-1);
}



/* ====================== >>> buildDeltEnergiesR <<< ======================= */
void buildDeltEnergiesRad()
{
  int i, a, b, csSq, d;
  double DrSq, vabCut, rabSq, srab2, srab6, srab12;

  calcConst();
  
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  cubeSize[a][b] = 
	    ((int) ( (Oparams.sigab[a][b]*Oparams.rcut)/Oparams.lattice_a) )+1
	    + 1;
	  /* 
	     questo +1 serve poiché vengono costruite tabelle do lookup
	     per differenze di energie rlativi a siti che differiscono 
	     al piu' per un passo reticolare
	   */
	  maxdelta[a][b] = 1+2*cubeSize[a][b];
	  csSq = Sqr(cubeSize[a][b]);
	  printf("devo allocare: %d bytes!!!!\n", (2*maxdelta[a][b]+1)*sizeof(double)*csSq);
	  //	  exit(-1);

	  DvabRadLT[a][b] = AllocMatR(6*maxdelta[a][b]+1, sizeof(double)*csSq);
	  DwabRadLT[a][b] = AllocMatR(6*maxdelta[a][b]+1, sizeof(double)*csSq);
	}
    }
    
  
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;

	  for (d = 0; d < 6*maxdelta[a][b]+1; d++)
	    {
	      
	      csSq = Sqr(cubeSize[a][b]);
	      
	      for (i = 0; i < csSq; i++)
		{
		  
		  rabSq = Sqr((double)Oparams.lattice_a) * ((double)i);
		  DrSq = Sqr((double)Oparams.lattice_a) 
		    * ((double)(d - 3*maxdelta[a][b])); 
		  if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      
		      srab2   = sigabSq[a][b] / rabSq;
		      srab6   = srab2 * srab2 * srab2;
		      srab12  = Sqr(srab6);
		      /* Calcola il potenziale e il viriale includendo
			 i fattori moltiplicativi */
		      //printf("(%d) 12-6:%.20f\n", i,srab12 - srab6);
		      DvabRadLT[a][b][d][i] = epsab4[a][b]*
			((srab12 - srab6) - vabCut);
		      DwabRadLT[a][b][d][i] = 
			epsab24[a][b]*((srab12 - srab6) + srab12)/3.0;
		    }
		  else
		    {
		      DvabRadLT[a][b][d][i] = 0.0;
		      DwabRadLT[a][b][d][i] = 0.0;
		    }

		  rabSq += DrSq; 
		  if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      srab2   = sigabSq[a][b] / rabSq;
		      srab6   = srab2 * srab2 * srab2;
		      srab12  = Sqr(srab6);
		      
		      DvabRadLT[a][b][d][i] = epsab4[a][b]*
			((srab12 - srab6) - vabCut) - DvabRadLT[a][b][d][i];
		      DwabRadLT[a][b][d][i] = 
			epsab24[a][b]*((srab12 - srab6) + srab12)/3.0 -
			DwabRadLT[a][b][d][i];
		    }
		  if (i == 0)printf("[%d,%d]\n", i, d);
		  if (i == 324 && d == 140)
		  printf("la: %.6f epsab4: %.10f DrSq: %.10f [%d,%d,%d,%d]: %.10f \n ", 
			 Oparams.lattice_a, epsab4[a][b], DrSq,a, b, d, i, DvabRadLT[a][b][d][i]);
		}
	    }
	}
    }
  exit(-1);
}

/* ========================== >>> builgPMap <<< ========================= */
void buildMaps(void)
{
  int i, lM, Nm;
  int *buf;
  lM = Oparams.lattice_M;
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pmap = buf + lM;

  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pabs = buf + lM;

  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pSqr = buf + lM;

  pgetj = (int *) malloc((2*Nm)*sizeof(int));
  pgetb = (int *) malloc((2*Nm)*sizeof(int));

  for (i = -lM; i <= lM; i++)
    {
      pmap[i] = minimumImageI(i);
      pabs[i] = abs(i); 
      pSqr[i] = Sqr(i);
     
    }
  
  for(i = 0; i < NA*Nm; i++)
    {  
      if (i >= Nm)
	{
	  pgetb[i] = 1;
	  pgetj[i] = i - Nm;
	}
      else
	{
	  pgetb[i] = 0;
	  pgetj[i] = i;
	}
    }
}


/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */

  int Nm, i, sct;
  COORD_TYPE* m;
  int a;

  //COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;

  /* initialize global varibales */
  pi = 2.0 * acos(0);
  
  mdMsg(ALL, NOSYS, "usrInitAft", "NOTICE", NULL,
	"Using Neighbour List method",
	NULL);
    
  /* [0][1] elements are read from parameter file */
  Oparams.sigab[1][0] = Oparams.sigab[0][1];
  Oparams.epsab[1][0] = Oparams.epsab[0][1];
  
  /* Allocate arrays of integers, you must supply the array length and the
     address of the pointer to the array */
  Nm = Oparams.parnum[0] + Oparams.parnum[1];

  sct = sizeof(COORD_TYPE);
  
  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac;
  nebrNow = 1;
  nebrTab = AllocMatI(NA*Nm, nebrTabMax); 
  nebrTabi = AllocMatI(NA*Nm, nebrTabMax); 
  nebrTaba = AllocMatI(NA*Nm, nebrTabMax); 
 
  nebrTabLen = malloc(sizeof(int)*Nm*NA);
  /* =============================================== */
  
  for (i=0; i < NA*Nm; i++)
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
      for ( a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      
	      /* store the initial positions of particles */
	      OprogStatus.rxi[a][i] = rx[a][i];
	      OprogStatus.ryi[a][i] = ry[a][i];
	      OprogStatus.rzi[a][i] = rz[a][i];
	    }
	}
      
      loop(i, 1, NUMK) 
	{
	  OprogStatus.sumS[i] = 0.0;
	}
      
    }
  printf("Vol: %.15f\n", Vol);

  for (a = 0; a < NA; a++)
    {
      /* coordinate scalate alla prima scatola */
      rxS[a] = (int *)malloc(sizeof(int)*Oparams.parnum[a]);
      ryS[a] = (int *)malloc(sizeof(int)*Oparams.parnum[a]);
      rzS[a] = (int *)malloc(sizeof(int)*Oparams.parnum[a]);
    }

  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* coordinate scalate alla prima scatola */
	  //printf("(%d,%d,%d)\n", rx[a][i], ry[a][i], rz[a][i]);
	  rxS[a][i] = minimumImageI(rx[a][i]);
	  ryS[a][i] = minimumImageI(ry[a][i]);
	  rzS[a][i] = minimumImageI(rz[a][i]);
	}
    }
  /* calcola le energie potenziali e il viriale all'inizio della simulazione */
  sumup();
  
  /* calcola tutte le energie di interazione all'inizio */

#ifdef LMC_FAST
  buildDeltEnergiesRad();
#else
  buildEnergiesRad();
#endif

  //buildEnergies();
  buildMaps();

}


/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i, a;
  double la;
  la = Oparams.lattice_a;
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* salva le coordinate come numeri in virgola mobile */
	  fprintf(fs, "%.15G %.15G %.15G\n", la*((double)rx[a][i]), 
		  la*((double)ry[a][i]), la*((double)rz[a][i]));
	}
    }

}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i, a;
  double xx,yy,zz, la;
  
  la = Oparams.lattice_a;

  for (a = 0;  a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &xx, &yy, &zz) < 3)
	    {
	      mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	      exit(-1);
	    }
	  /* converte le coordinate negli interi corrisponednti */
	  rx[a][i] = rint(xx/la);
	  ry[a][i] = rint(yy/la);
	  rz[a][i] = rint(zz/la);
	}
    }
}
