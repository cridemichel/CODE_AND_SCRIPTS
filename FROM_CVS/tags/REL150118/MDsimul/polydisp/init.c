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

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE V, W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx;  

/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */


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
		  rx[a][np[a]] = rxCm;
		  ry[a][np[a]] = ryCm;
		  rz[a][np[a]] = rzCm;
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
	  rx[a][i] = rx[a][i] - 0.5 * L; 
	  ry[a][i] = ry[a][i] - 0.5 * L;
	  rz[a][i] = rz[a][i] - 0.5 * L;
	}
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

/* See after for declaration */
void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm);

void resetCM()
{
  COORD_TYPE sumx, sumy, sumz, RCMx, RCMy, RCMz;
  COORD_TYPE *m, MTOT;
  int i, a;

  m = Oparams.m;
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  MTOT = 0.0;

  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  sumx = sumx + m[a]*vx[a][i];
	  sumy = sumy + m[a]*vy[a][i];
	  sumz = sumz + m[a]*vz[a][i];
	  MTOT += m[a];
	}
      
    }
  sumx = sumx / MTOT; 
  sumy = sumy / MTOT;
  sumz = sumz / MTOT;

  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  vx[a][i] = vx[a][i] - sumx;
	  vy[a][i] = vy[a][i] - sumy;
	  vz[a][i] = vz[a][i] - sumz;
	  /* In this way the total (net) momentum of the system of 
	     molecules is zero */
	}
    }
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
	  RCMx += m[a]*rx[a][i]; /*Here RCM is the center of mass of the box */
	  RCMy += m[a]*ry[a][i];
	  RCMz += m[a]*rz[a][i];
	}
    }
  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;

  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  rx[a][i] -= RCMx;
	  ry[a][i] -= RCMy;
	  rz[a][i] -= RCMz;
	}
    }
}

/* ========================== >>> comvel <<< =============================== */
void comvel (COORD_TYPE temp, COORD_TYPE m[NA])
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
  COORD_TYPE rTemp[NA], sumx, sumy, sumz, RCMx, RCMy, RCMz;
  double MTOT;
  //COORD_TYPE Px, Py, Pz;
  int i, a;
  
  rTemp[0] = sqrt(temp / Oparams.m[0]);
  rTemp[1] = sqrt(temp / Oparams.m[1]);
  /* variance of the velocities distribution function, we assume k = 1 */ 

  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vx[a][i] = rTemp[a] * gauss(); 
	  vy[a][i] = rTemp[a] * gauss();
	  vz[a][i] = rTemp[a] * gauss();
	  /* gauss() is a gaussian variable with mean = 0 and variance = 1, 
	     that is
	     2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
	}
    }
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  MTOT = 0.0;
  for (a = 0;  a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* (sumx, sumy, sumz) is the total momentum */ 
	  sumx = sumx + m[a] * vx[a][i];
	  sumy = sumy + m[a] * vy[a][i];
	  sumz = sumz + m[a] * vz[a][i];
	  MTOT += m[a];
	}
    }
  sumx = sumx / MTOT; 
  sumy = sumy / MTOT;
  sumz = sumz / MTOT;

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  
  for (a = 0;  a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vx[a][i] = vx[a][i] - sumx;
	  vy[a][i] = vy[a][i] - sumy;
	  vz[a][i] = vz[a][i] - sumz;
	  //printf("([%d,%d]: %.5f %.5f %.5f)\n",a , i, vx[a][i], vy[a][i],
	  // vz[a][i]);
	  /* In this way the total (net) momentum of the system of 
	     molecules is zero */
	}
    }
  
  
  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  for (a = 0;  a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  
	  RCMx += m[a] * rx[a][i];/*Here RCM is the center of mass of the box*/
	  RCMy += m[a] * ry[a][i];
	  RCMz += m[a] * rz[a][i];
	}
    }
  
  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;

  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  rx[a][i] -= RCMx;
	  ry[a][i] -= RCMy;
	  rz[a][i] -= RCMz;
	}
    }
  /* Now the center of mass of the box is in the origin */
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
  
  /* set both atoms velocity to the center of mass velocity */
  comvel(Oparams.T, Oparams.m); 

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
  Dtrans[0] = Dtrans[1] = 0.0; /* DtransOld should become a field of OprogStatus */

  Vol = 1400.0;
  Vol1 = 0.0;
  Vol2 = 0.0;
  Vol1o1 = 0.0;
  Vol1o2 = 0.0;

  /* initialize the Nose parameter */
  s = 1.0;
  s1 = 0.0;
  s2 = 0.0;

  Oparams.T = 2.0;
  Oparams.P = 1.0;
  Oparams.M = 5; /* cells in each direction for linked lists */
  Oparams.d = 0.025;
  Oparams.tol = 0.0000001;
  
  OprogStatus.Q = 1.0;  /* Default value of the parameter of the Nose method */
  OprogStatus.Nose = 1; /* Use nose method by default */ 
  OprogStatus.W = 0.0027;
  OprogStatus.avVol = 0.0;
  OprogStatus.avs   = 0.0;
  OprogStatus.tolVol = 0.0001;
  OprogStatus.tols  = 0.001;
  OprogStatus.tolVol1 = 0.1;
  OprogStatus.sResetSteps = 0; /* Don't reset s by default */
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  OprogStatus.noLinkedList = 0; /* Use Linked List */
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  OprogStatus.avngS     = 0;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
  OprogStatus.avngEta   = 0;

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
  
  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  /* ======================================================================= */
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
  int a, b;

  //COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;

  /* initialize global varibales */
  pi = 2.0 * acos(0);
  
  M = Oparams.M;

  if (M<=2) 
    {
      printf("ERROR: cellNum must be > 2 !!!\n");
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
  /* [0][1] elements are read from parameter file */
  Oparams.sigab[1][0] = Oparams.sigab[0][1];
  Oparams.epsab[1][0] = Oparams.epsab[0][1];
  
  /* Allocate arrays of integers, you must supply the array length and the
     address of the pointer to the array */
  Nm = Oparams.parnum[0]+Oparams.parnum[1];
  head = malloc(sizeof(int)*NCell);
  list = malloc(sizeof(int)*NA*Nm);
  map  = malloc(sizeof(int)*mapSize);

  sct = sizeof(COORD_TYPE);
  
  /* NOTE:
     Inside shared loops some of the work is done by the Father, and some
     by the Child, so if you want something that depends upon the works of both
     processes you must put their calculations in the shared memory, and 
     not in the process local memory (Not shared) */ 
  maps();

  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac * Nm * NA;
  printf("INIT nebrTabMax: %d\n", nebrTabMax);
  nebrNow = 1;
  nebrTab = AllocMatI(2, nebrTabMax); 
  /* =============================================== */
  
  /* Store the Center of Mass initial position for all particles */
  m = Oparams.m;

  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      if (OprogStatus.Nose == 0)
	{
	  /* If we begin a new simulation and we don't use Nose-Andersen method
	     the velocity and acceleration of Volume and s must be zero */
	  Vol1  = 0.0;
	  Vol2  = 0.0;
	  Vol1o1 = 0.0;
	  Vol1o2 = 0.0;
	  s1    = 0.0;
	  s2    = 0.0;
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
	}
     
      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;
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
      OprogStatus.sumEta   = 0.0;
      OprogStatus.sumTemp  = 0.0;
      OprogStatus.sumPress = 0.0;
      
      loop(i, 1, NUMK) 
	{
	  OprogStatus.sumS[i] = 0.0;
	}
      
      loop(i, 1, MAXBIN)
	{
	  OprogStatus.hist[i] = 0;
	}
    }
  printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);
}


/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i, a;
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", rx[a][i], ry[a][i], rz[a][i]);
	}
    }
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", vx[a][i], vy[a][i], vz[a][i]);
	}
    }
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", Fx[a][i], Fy[a][i], Fz[a][i]);
	}
    }

  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G\n", Vol, Vol1, Vol2, Vol1o1, Vol1o2);
  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G",  s, s1, s2, s1o1, s1o2);

}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i, a;

  for (a = 0;  a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &rx[a][i], &ry[a][i], &rz[a][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &vx[a][i], &vy[a][i], &vz[a][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
  
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &Fx[a][i], &Fy[a][i], &Fz[a][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
  

  if (fscanf(fs, "%lf %lf %lf %lf %lf",  &Vol, &Vol1, &Vol2, &Vol1o1, &Vol1o2) < 5)
    {
      mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
      exit(-1);
    }
  
  if (fscanf(fs, "%lf %lf %lf %lf %lf",  &s, &s1, &s2, &s1o1, &s1o2) < 5)
    {
      mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
      exit(-1);
    }
      
}
