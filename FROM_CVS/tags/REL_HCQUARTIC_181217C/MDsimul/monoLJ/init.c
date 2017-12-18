#include<mdsimul.h>
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
extern COORD_TYPE Vc, V, W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz;  

/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
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
void FCC(int Nm, COORD_TYPE m)
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
		  rx[ii+iref] = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ry[ii+iref] = by[iref] + Cell * ((COORD_TYPE) iy);
		  rz[ii+iref] = bz[iref] + Cell * ((COORD_TYPE) iz);
		}
	      
	      ii = ii + 4;
	      
	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */
  
  loop(i, 1, Nm)
    {
      /* Initial position values are between -0.5 and 0.5 */
      rx[i] = rx[i] - 0.5 * L; 
      ry[i] = ry[i] - 0.5 * L;
      rz[i] = rz[i] - 0.5 * L;
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

void resetCM(int Nm)
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
	 sumx = sumx + vx[i];
	 sumy = sumy + vy[i];
	 sumz = sumz + vz[i];
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
      RCMx += rx[i]; /* Here RCM is the center of mass of the box */
      RCMy += ry[i];
      RCMz += rz[i];
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  loop(i, 1, Nm)
    {
      rx[i] -= RCMx;
      ry[i] -= RCMy;
      rz[i] -= RCMz;
    }
}

/* ========================== >>> comvel <<< =============================== */
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE m)
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
      RCMx += rx[i]; /* Here RCM is the center of mass of the box */
      RCMy += ry[i];
      RCMz += rz[i];
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  loop(i, 1, Nm)
    {
      //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
      rx[i] -= RCMx;
      ry[i] -= RCMy;
      rz[i] -= RCMz;
    }
  /* Now the center of mass of the box is in the origin */
}

/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  setToZero(SAVE_LIST, 
	    NULL);  /* Set to zero all the coordinates */

  FCC(Oparams.parnum, Oparams.m); 
  /* Put the baricenter of each molecule on a FCC lattice, and set 
     their orientations */  
  
  /* set both atoms velocity to the center of mass velocity */
  comvel(Oparams.parnum, Oparams.T, Oparams.m); 
  
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
  //angvel(Oparams.parnum, Oparams.T, Oparams.m, Oparams.d); 
}

/* =========================== >>> usrInitBef <<< ========================== */
void usrInitBef(void)
{
  int i;
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

  /* initialize the Nose parameter */
  s = 1.0;
  s1 = 0.0;
  s2 = 0.0;

  Oparams.T = 2.0;
  Oparams.P = 1.0;
  Oparams.M = 5; /* cells in each direction for linked lists */
  
  OprogStatus.Q = 1.0;  /* Default value of the parameter of the Nose method */
  OprogStatus.Nose = 1; /* Use nose method by default */ 
  OprogStatus.W = 0.0027;
  OprogStatus.avVol = 0.0;
  OprogStatus.avs   = 0.0;
  OprogStatus.tolVol = 0.0001;
  OprogStatus.tols  = 0.001;
  OprogStatus.tolVol1 = 0.1;
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

  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  for (i = 0; i < PE_POINTS; i++)
    OprogStatus.PE[i] = 0;
  
  /* ======================================================================= */
}

extern buildMesh(C_T, C_T);

/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */

  int Nm, i, sct;
  COORD_TYPE vcmx, vcmy, vcmz;
  COORD_TYPE m;
  char msgs[1024];
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
  head = malloc(sizeof(int)*NCell);
  list = malloc(sizeof(int)*Oparams.parnum);
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
  nebrNow = 1;
  nebrTab = AllocMatI(2, nebrTabMax); 
  /* =============================================== */
  
  m = Oparams.m;

  buildMesh(OprogStatus.kmax, OprogStatus.Dk);

  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
#if defined(SOFT_SPHERE)
#if defined(MD_STATIC_PP36)
      Oparams.PP = 36;
      mdPrintf(ALL, "WARNING: MD_STATIC_PP36 defined!\n", NULL);
#else
      sprintf(msgs, "WARNING: DYNAMIC PP=%d\n", Oparams.PP);
      mdPrintf(ALL, msgs, NULL);
#endif
#endif
      if (OprogStatus.Nose == 0)
	{
	  /* If we begin a new simulation and we don't use Nose-Andersen method
	     the velocity and acceleration of Volume and s must be zero */
	  Vol1  = 0.0;
	  Vol2  = 0.0;
	  Vol1o1 = 0.0;
	  Vol1o2 = 0.0;
	  s = 1.0;
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
	  /* setting zeroall_s it is possible to reset derivatives to 0
	     and s to 1, this is useful if you want ta make the temperature
	     piston mass large to perform quasi-NVE */
	  if (OprogStatus.zeroall_s)
	    {
	      s = 1.0;
	      s1 = s1o1 = s1o2 = 0.0;
	    }

	}
     
      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;
        
      loop(i, 1, Oparams.parnum)
	{
	  /* store the initial positions of particles */
	  OprogStatus.rxCMi[i] = rx[i];
	  OprogStatus.ryCMi[i] = ry[i];
	  OprogStatus.rzCMi[i] = rz[i];	 
	  
	  /* Center of mass velocities */
	  vcmx = vx[i];
	  vcmy = vy[i];
	  vcmz = vz[i];
	
	  OprogStatus.vcmx0[i] = vcmx;
	  OprogStatus.vcmy0[i] = vcmy;
	  OprogStatus.vcmz0[i] = vcmz;
	  
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
  // printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);
}


/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i;
  const char tipodat[] = "%.15G %.15G %.15G\n";
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat, rx[i], ry[i], rz[i]);
    }
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat, vx[i], vy[i], vz[i]);
    }

  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G %.15G",  Vol, s, s1, s2, s1o1, s1o2);

}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i;

  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf\n", &rx[i], &ry[i], &rz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf\n", &vx[i], &vy[i], &vz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
      
    }
  

  if (fscanf(fs, "%lf %lf %lf %lf %lf %lf",  &Vol, &s, &s1, &s2, &s1o1, &s1o2) < 3)
    {
      mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
      exit(-1);
    }
      
}
