#include <mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

extern char TXT[MSG_LEN];
extern int ENDSIM;
extern char msgStrA[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   

extern COORD_TYPE Vc[MAX_M], V[MAX_M], W[MAX_M], K[MAX_M], 
  WC[MAX_M], T1xx[MAX_M], T1yy[MAX_M], T1zz[MAX_M],
  T1xx[MAX_M], T1yy[MAX_M], T1zz[MAX_M], T1xy[MAX_M], 
  T1yz[MAX_M], T1zx[MAX_M], Wxx[MAX_M], Wyy[MAX_M], Wzz[MAX_M],
  Wxy[MAX_M], Wyz[MAX_M], Wzx[MAX_M], Pxx[MAX_M], Pyy[MAX_M], 
  Pzz[MAX_M], Pxy[MAX_M], Pyz[MAX_M], Pzx[MAX_M], 
  Patxy[MAX_M], Patyz[MAX_M], Patzx[MAX_M], Patxx[MAX_M], Patyy[MAX_M], 
  Patzz[MAX_M];  

/* used by linked list routines */
extern int **head, **list, *map;  /* arrays of integer */

extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int ***nebrTab, *nebrNow, *nebrTabLen, nebrTabMax;
int *kcx, *kcy, *kcz;
int NNC;
/* ================================= */

/* ========================================================================= */
/*=========================== >>> vectProd <<< =========================== */
void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z)
{
  /* DESCRIPTIOM:
     r3 = [ r1, r2 ] where [ , ] the vectorial product */
  *r3x = r1y * r2z - r1z * r2y; 
  *r3y = r1z * r2x - r1x * r2z;
  *r3z = r1x * r2y - r1y * r2x;
}


/* ============================= >>> FCC <<< ================================*/
void FCC(int ss, int Nm, COORD_TYPE m)
{
  /*   DESCRIPTION:
       Sets up the alpha fcc lattice for n linear molecules.   
       The simulation box is a unit cube centred at the origin.
       N should be an integer of the form ( 4 * ( Nc ** 3 ) ),
       Where Nc is the number of FCC unit cells in each direction.  
       See figure 5.10 for a diagram of the lattice and a           
       definition of the four orientational sublattices.            
       PRINCIPAL VARIABLES:                                         
       int           Nm                   Number of molecules             
       COORD_TYPE    d                    diatomic molecule length 
       COORD_TYPE    rxCm, ryCm, rzCm     Molecular Center of mass 
                                          positions             
       COORD_TYPE    ex, ey, ez           half of vector joining atom a and b 
                                          in a molecule 
       COORD_TYPE    rRoot3               1.0 / sqrt ( 3.0 ) */
  int Nc;
  COORD_TYPE rRoot3; // = 0.5773503;
  COORD_TYPE  Cell, Cell2;
  int i, ix, iy, iz, iref, ii;
  COORD_TYPE bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  
  //printf("FCC Vol: %f\n", Vol);
  L = cbrt(Vol);
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
  
  loop(iz, 1, Nc) /* loops over unit cells (that are simply cubes) */ 
    {
      loop(iy, 1, Nc)
	{
	  loop(ix, 1, Nc)
	    {
	      loop(iref, 1, 4) /* In each primitive cell there are four 
				  molecules */
		{
		  if ((ii + iref) >= Nm) break;
		  /* If Nm < 4 * Nc^3 the we have more lattice sites than 
		     particles so we stop if all the particles was positioned.
		     This condition in fact means: 'If all the particles was
		     positioned => end'
		  */

		  /* Center of Mass of the actual molecule (m + iref) */
		  rx[ss][ii+iref] = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ry[ss][ii+iref] = by[iref] + Cell * ((COORD_TYPE) iy);
		  rz[ss][ii+iref] = bz[iref] + Cell * ((COORD_TYPE) iz);
		}
	      
	      ii = ii + 4;
	      
	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */
  
  loop(i, 1, Nm)
    {
      /* Initial position values are between -0.5 and 0.5 */
      rx[ss][i] = rx[ss][i] - 0.5 * L; 
      ry[ss][i] = ry[ss][i] - 0.5 * L;
      rz[ss][i] = rz[ss][i] - 0.5 * L;
      //printf("ss=%d i=%d (%.10f,%.10f,%.10f)\n", ss, i, 
      //     rx[ss][i], ry[ss][i], rz[ss][i]);
    }
  return;
}

/* Centre of mass and angular velocities for linear molecules    
   PRINCIPAL VARIABLES:                                       
   int           Nm                   The number of molecules        
   COORD_TYPE    rx[Nm],ry[Nm],rz[Nm] Positions                      
   COORD_TYPE    vx[Nm],vy[Nm],vz[Nm] Velocities                     
   COORD_TYPE    ex[Nm],ey[Nm],ez[Nm] Orientations 
   COORD_TYPE    ox[Nm],oy[Nm],oz[Nm] Space-fixed angular velocities 
   COORD_TYPE    temp                 Reduced temperature            
   COORD_TYPE    inert                moment of inertia      
                                                            
 SUPPLIED ROUTINES:                                           
                                                              
  void comvel(Nm, temp, m)                                   
    Sets the centre of mass velocities for a configuration of 
    linear molecules at a given temperature.                  
  void angvel(Nm, temp, m, d)                            
    Sets the angular velocities for a configuration of linear 
    molecules at a given temperature.                         
  COORD_TYPE ranf(void)                                 
    Returns a uniform random variate on the range zero to one 
  COORD_TYPE gauss(void)                                
    Returns a uniform random normal variate from a            
    distribution with zero mean and unit variance.            
                                                              
 UNITS:                                                       
                                                              
 We employ Lennard-Jones units 
 PROPERTY                    UNITS                     
       rx, ry, rz           (epsilon/Mtot)**(1.0/2.0)             
       ox, oy, oz           (epsilon/Mtot*sigma**2)**(1.0/2.0)    
       inert                 Mtot*sigma**2                        
*/
#ifdef MPI
extern int my_rank;
#endif
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

void resetCM(int ss, int Nm)
{
  COORD_TYPE sumx, sumy, sumz, RCMx, RCMy, RCMz;
  COORD_TYPE m;
  int i;

  m = Oparams.m;
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  loop(i, 1, Nm)
       {
	 sumx = sumx + vx[ss][i];
	 sumy = sumy + vy[ss][i];
	 sumz = sumz + vz[ss][i];
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  loop(i, 1, Nm)
    {
      vx[ss][i] = vx[ss][i] - sumx;
      vy[ss][i] = vy[ss][i] - sumy;
      vz[ss][i] = vz[ss][i] - sumz;
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }
  
  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  loop(i, 1, Nm)
    {
      RCMx += rx[ss][i]; /* Here RCM is the center of mass of the box */
      RCMy += ry[ss][i];
      RCMz += rz[ss][i];
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  loop(i, 1, Nm)
    {
      rx[ss][i] -= RCMx;
      ry[ss][i] -= RCMy;
      rz[ss][i] -= RCMz;
    }
}

/* ========================== >>> comvel <<< =============================== */
void comvel (int ss, int Nm, COORD_TYPE temp, COORD_TYPE m)
{
  /*
    Translational velocities from maxwell-boltzmann distribution  
    The routine put in vx, vy, vz a velocity choosen from a M.-B. 
    distribution.
    
    The distribution is determined by temperature and (unit) mass.
    This routine is general, and can be used for atoms, linear    
    molecules, and non-linear molecules.                          
    
    ROUTINE REFERENCED:                                          
    
    COORD_TYPE gauss(void)
    Returns a uniform random normal variate from a           
    distribution with zero mean and unit variance.           
    
    VARIABLES 
    COORD_TYPE temp       Temperature 
    m[NA]                 Masses of atoms (NA is the number of atoms)
    int  Nm               Number of molecules  
  */
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
      
      vx[ss][i] = rTemp * gauss(); 
      vy[ss][i] = rTemp * gauss();
      vz[ss][i] = rTemp * gauss();
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
	 sumx = sumx + vx[ss][i];
	 sumy = sumy + vy[ss][i];
	 sumz = sumz + vz[ss][i];

	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  loop(i, 1, Nm)
    {
      vx[ss][i] = vx[ss][i] - sumx;
      vy[ss][i] = vy[ss][i] - sumy;
      vz[ss][i] = vz[ss][i] - sumz;
      //printf("ss=%d i=%d VVV(%.10f,%.10f,%.10f)\n", ss, i, 
      //     vx[ss][i], vy[ss][i], vz[ss][i]);
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }

  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  loop(i, 1, Nm)
    {
      RCMx += rx[ss][i]; /* Here RCM is the center of mass of the box */
      RCMy += ry[ss][i];
      RCMz += rz[ss][i];
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  loop(i, 1, Nm)
    {
      //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
      rx[ss][i] -= RCMx;
      ry[ss][i] -= RCMy;
      rz[ss][i] -= RCMz;
    }
  /* Now the center of mass of the box is in the origin */
}

extern void setToZero(COORD_TYPE** ptr, ...);

/* =========================== >>> initCoord <<< =========================*/
void initCoord(void)
{
  int ss;
  setToZero(SAVE_LIST, 
	    NULL);  /* Set to zero all the coordinates */

  for (ss = 0; ss < Oparams.PTM; ss++ )
    {
      FCC(ss, Oparams.parnum, Oparams.m); 
      /* Put the baricenter of each molecule on a FCC lattice, and set 
	 their orientations */  
      
      /* set both atoms velocity to the center of mass velocity */
      comvel(ss, Oparams.parnum, Oparams.T, Oparams.m); 
    }
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
  
  int ss, i;
  /* ===================== >>> INIT PARAMETERS <<< ======================== 
   All the values set here for Oparams structure are taken as defaults if you
   don't specify corresponding parameters in the parameters file  */
  

  Vol = 1400.0;
  Vol1 = 0.0;
  Vol2 = 0.0;
  Vol1o1 = 0.0;
  Vol1o2 = 0.0;

  /* initialize the Nose parameter */
  for (ss = 0; ss < MAX_M; ss++)
    {
      s[ss] = 1.0;
      s1[ss] = 0.0;
      s2[ss] = 0.0;
    }
  Oparams.T = 2.0;
  Oparams.P = 1.0;
  Oparams.M = 5; /* cells in each direction for linked lists */
  
  OprogStatus.Q = 1.0;  /* Default value of the parameter of the Nose method */
  OprogStatus.Nose = 1; /* Use nose method by default */ 
  OprogStatus.W = 0.0027;
  OprogStatus.avs   = 0.0;
  OprogStatus.tols  = 0.001;
  OprogStatus.tolT = 0.0;
  OprogStatus.sResetSteps = 0; /* Don't reset s by default */
  OprogStatus.zeroall_s = 0;
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  OprogStatus.noLinkedList = 0; /* Use Linked List */
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  OprogStatus.avngS     = 0;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
 
  Oparams.m = 1.0;
  
  Oparams.sigma = 1.0;
  Oparams.epsilon = 1.0;

  /* Parameters relative to Ruocco AC force
     See: cond-mat/00001311, Ruocco et al. */
  OprogStatus.alpha = 0.83;
  OprogStatus.kmax = 7.12;
  OprogStatus.Dk = 0.34;
  OprogStatus.S0 = 10.0;
  OprogStatus.ENmin = -4.0;
  OprogStatus.ENmax = -6.25;
  for (ss = 0; ss < MAX_M; ss++)
    {
      OprogStatus.scambi[ss] = 0;
      OprogStatus.attempted[ss] = 0;
      for (i = 0; i < MAXPAR; i++)
	{
	  OprogStatus.sumVx[ss][i] = 0.0;
	  OprogStatus.sumVy[ss][i] = 0.0;
	  OprogStatus.sumVz[ss][i] = 0.0;
	}
    }
  OprogStatus.RW = 0;
  OprogStatus.srescale = 0;
  OprogStatus.linearDiff = 0; 
  /* fa la differenza media fra i logaritmi delle distribuzioni riscalate 
     se 0 altrimenti fa la differenza fra le distr.*/
  
  OprogStatus.sogliaPEij = 1E-4; 
  /* per il calcolo della media confronta solo valori inferiori a tale 
     soglia evitando cosi' di fare differenze fra valori eccessivamente 
     improbabili e quindi rumorosi delle distribuzioni riscalate */
				    
  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  /* ======================================================================= */

  srand(time(NULL));
}

extern void buildMesh(C_T, C_T);

/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */

  int Nm, i, sct, ss;
  COORD_TYPE m;
  
  //COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;
 
  /* initialize global varibales */
  pi = 2.0 * acos(0);
  
  M = Oparams.M;

  if (M<=2) 
    {
      sprintf(TXT, "ERROR: cellNum must be > 2 !!!\n");
      mdPrintf(STD, TXT, NULL);
      exit(-1);
    }
  NCell = M * M * M;
  mapSize = 13 * NCell;
    
  if (OprogStatus.noLinkedList)
    {
      mdMsg(ALL, NOSYS, "usrInitAft", "NOTICE", NULL,
	    "Not using linked list method",
	    NULL);
    }
  else
    {
      mdMsg(ALL, NOSYS, "usrInitAft", "NOTICE", NULL,
	    "Using linked list method",
	    NULL);
    }
  
  /* Allocate arrays of integers, you must supply the array length and the
     address of the pointer to the array */
  head = AllocMatI(Oparams.PTM, NCell);
  list = AllocMatI(Oparams.PTM, Oparams.parnum);
  map  = malloc(sizeof(int)*mapSize);
  Nm = Oparams.parnum;
  sct = sizeof(COORD_TYPE);

  kcx = malloc(NK*sizeof(int));
  kcy = malloc(NK*sizeof(int));
  kcz = malloc(NK*sizeof(int));
  NNC = 0;

   /* NOTE:
     Inside shared loops some of the work is done by the Father, and some
     by the Child, so if you want something that depends upon the works of both
     processes you must put their calculations in the shared memory, and 
     not in the process local memory (Not shared) */ 
  maps();

  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac * Nm;
  nebrTab = (int ***) malloc(sizeof(int **)*Oparams.PTM);
  nebrNow = (int *) malloc(sizeof(int)*Oparams.PTM);
  nebrTabLen = (int *) malloc(sizeof(int)*Oparams.PTM);

  OprogStatus.NUM_PEij = ((int) (((double)Oparams.PTM) / 2.0));

  if (OprogStatus.PT == 0)
    {
      /* se non si sta usando il Parallel Tempering non 
	 si puo' usare il criteri ode lRandom Walk 
	 per stabilire se il sistema ha equilibrato */
      OprogStatus.RW = 2;
    }

  if ((OprogStatus.RW != 2) && (OprogStatus.RW != 0))
    {
      mdPrintf(STD, "Errore nel valore di RW: puo' essere solo 0 o 2!",
	       NULL);
      exit(-1);
    }

  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      if (Oparams.lambdat[ss] == 0)
	{
	  OprogStatus.refReplica = ss;
	}
      
      nebrNow[ss] = 1;
      nebrTab[ss] = AllocMatI(2, nebrTabMax); 
    }

  /* =============================================== */
  m = Oparams.m;

  buildMesh(OprogStatus.kmax, OprogStatus.Dk);

  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      for (ss = 0; ss < Oparams.PTM; ss++)
	{
	  OprogStatus.sumPress[ss] = 0.0;
	  loop(i, 1, NUMK) 
	    {
	      OprogStatus.sumS[ss][i] = 0.0;
	    }
	  
	  loop(i, 1, MAXBIN)
	    {
	      OprogStatus.hist[ss][i] = 0;
	    }

	  OprogStatus.sumEta[ss]   = 0.0;
	  OprogStatus.DQxy[ss] = 0.0;
	  OprogStatus.DQyz[ss] = 0.0;
	  OprogStatus.DQzx[ss] = 0.0;
        }

      if (OprogStatus.Nose == 0)
	{
	  /* If we begin a new simulation and we don't use Nose-Andersen method
	     the velocity and acceleration of Volume and s must be zero */
	  Vol1  = 0.0;
	  Vol2  = 0.0;
	  Vol1o1 = 0.0;
	  Vol1o2 = 0.0;

	  for (ss = 0; ss < Oparams.PTM; ss++)
	    {
	      s[ss] = 1.0;
	      s1[ss]    = 0.0;
	      s2[ss]    = 0.0;
	    }
	  /* 25/10/2000 ADDED 
	     Usa le forze LJ, cosi' se una precedente simulazione
	     e' stata fatta con il Nose' dovrebbero ridurre i casini */
	  for (i = 0; i < Oparams.parnum; i++)
	    {
	      Fx[ss][i] = FxLJ[ss][i];
	      Fy[ss][i] = FyLJ[ss][i];
	      Fz[ss][i] = FzLJ[ss][i];
	    }
	}
      
      if (OprogStatus.Nose == 2)/* NTV=ensemble */
	{
	  Vol1  = 0.0;
	  Vol2  = 0.0;
	  Vol1o1 = 0.0;
	  Vol1o2 = 0.0;
	  /*
	    if (!OprogStatus.sResetSteps)
	    {
	    s = 1.0;
	    s1 = 0.0;
	    }
	  */
	  /* setting zeroall_s it is possible to reset derivatives to 0
	     and s to 1, this is useful if you want ta make the temperature
	     piston mass large to perform quasi-NVE */
	  if (OprogStatus.zeroall_s)
	    {
	      for (ss = 0; ss < Oparams.PTM; ss++)
		{
		  s[ss] = 1.0;
		  s1[ss] = s1o1[ss] = s1o2[ss] = 0.0;
		}
	    }

	}
      
      OprogStatus.sumMedia = 0.0;
      OprogStatus.sumTemp  = 0.0;
      
    }
  // printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);
}

/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i, ss;
  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", rx[ss][i], ry[ss][i], rz[ss][i]);
	}
    }
  for (ss = 0; ss < Oparams.parnum; ss++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", vx[ss][i], vy[ss][i], vz[ss][i]);
	}
    }
  
  for (ss = 0; ss < Oparams.parnum; ss++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", Fx[ss][i], Fy[ss][i], Fz[ss][i]);
	}
    }

  for (ss = 0; ss < Oparams.PTM; ss++)
    {
      fprintf(fs, "%.15G %.15G %.15G %.15G %.15G\n",  
	      s[ss], s1[ss], s2[ss], s1o1[ss], s1o2[ss]);
    }
}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i,ss;

  for (ss = 0; ss < Oparams.parnum; ss++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &rx[ss][i], &ry[ss][i], &rz[ss][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
  
  for (ss = 0; ss < Oparams.parnum; ss++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &vx[ss][i], &vy[ss][i], &vz[ss][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
  
  for (ss = 0; ss < Oparams.parnum; ss++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &Fx[ss][i], &Fy[ss][i], &Fz[ss][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
  
  for (ss = 0; ss < Oparams.parnum; ss++)
    {
      if (fscanf(fs, "%lf %lf %lf %lf %lf", 
		 &s[ss], &s1[ss], &s2[ss], &s1o1[ss], &s1o2[ss]) < 3)
	{
	  mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
}
