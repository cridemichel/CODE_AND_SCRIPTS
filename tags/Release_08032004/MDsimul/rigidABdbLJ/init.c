#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

extern volatile int *sc_;           /* shared counter */
extern volatile int *semA_, *semB_; /* semaphores */

extern int process;          /* 0 = child 1 = father ( see mdsimul.c )*/ 

extern int ENDSIM;
extern char msgStrA[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, WC, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Wm, Wmxx, Wmyy, Wmzz, 
  Wmxy, Wmyz, Wmzx, Pmxx, Pmyy, Pmzz, Pmxy, Pmyz, Pmzx, T1mxy, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
  T1myz, T1mzx, T1mxx, T1myy, T1mzz;  
extern COORD_TYPE Mtot;
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

extern COORD_TYPE *FxL[NA], *FyL[NA], *FzL[NA]; /* Local arrays of forces
					    (used by the LJForce routine 
					    to avoid interferences) */
extern COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

extern COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
extern COORD_TYPE  *Rmx, *Rmy, *Rmz;
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
void FCC(int Nm, COORD_TYPE d, COORD_TYPE* m)
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
  COORD_TYPE dist;
  COORD_TYPE rRoot3; // = 0.5773503;
  COORD_TYPE  Cell, Cell2, rxCm, ryCm, rzCm, Mtot, fact;
  int a, i, ix, iy, iz, iref, ii;
  COORD_TYPE bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  COORD_TYPE ex[4], ey[4], ez[4]; /* orientations of each molecule in the 
				     primitive cell */
  COORD_TYPE d0x, d0y, d0z, d1x, d1y, d1z;

  //printf("FCC Vol: %f\n", Vol);
  L = cbrt(Vol);
  Nc = ceil(  pow( ((COORD_TYPE)Nm)/4.0, 1.0/3.0 )  );
  Mtot = m[0] + m[1];
  //printf("Nc: %d\n", Nc);
  /* Calculate the side of the unit cell */
  Cell  = L / ((COORD_TYPE) Nc); /* unit cell length */
  Cell2 = 0.5 * Cell;              /* half unit cell length */

  /* Sublattice A */
  rRoot3 = 1.0 / sqrt(3.0);
  bx[0] =  0.0;
  by[0] =  0.0;
  bz[0] =  0.0;
  ex[0] =  rRoot3 ;
  ey[0] =  rRoot3;
  ez[0] =  rRoot3;
  
  /*  Sublattice B */
  
  bx[1] =  Cell2;
  by[1] =  Cell2;
  bz[1] =  0.0;
  ex[1] =  rRoot3;
  ey[1] = -rRoot3;
  ez[1] = -rRoot3;
  
  /* Sublattice C */
  
  bx[2] =  0.0;
  by[2] =  Cell2;
  bz[2] =  Cell2;
  ex[2] = -rRoot3;
  ey[2] =  rRoot3;
  ez[2] = -rRoot3;
  
  /* Sublattice D */
  
  bx[3] =  Cell2;
  by[3] =  0.0;
  bz[3] =  Cell2;
  ex[3] = -rRoot3;
  ey[3] = -rRoot3;
  ez[3] =  rRoot3;	       
  
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
		  rxCm = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ryCm = by[iref] + Cell * ((COORD_TYPE) iy);
		  rzCm = bz[iref] + Cell * ((COORD_TYPE) iz);
		  //printf("CM:(%f,%f,%f)\n", rxCm, ryCm, rzCm);
		  /* (d0x, d0y, d0z) is the vector from the Ceneter of Mass
                     to the atom 0 */
		  fact = - d * m[1] / Mtot;
		  
		  d0x  = ex[iref] * fact;
		  d0y  = ey[iref] * fact;
		  d0z  = ez[iref] * fact;
 
		  /* (d1x, d1y, d1z) is the vector from the Ceneter of Mass
                     to the atom 1 */
		  fact = d * m[0] / Mtot;
		  
		  d1x  = ex[iref] * fact;
		  d1y  = ey[iref] * fact;
		  d1z  = ez[iref] * fact;
		  

		  /* The positions of two atoms are obtained by a displacement 
		     from the Center of Mass given by adding or subtracting the
		     d1 and d0 vectors (see above) */
		  rx[0][ii + iref] = rxCm + d0x;
		  rx[1][ii + iref] = rxCm + d1x;
		  ry[0][ii + iref] = ryCm + d0y;
		  ry[1][ii + iref] = ryCm + d1y;
		  rz[0][ii + iref] = rzCm + d0z;
		  rz[1][ii + iref] = rzCm + d1z;
		}
	      
	      ii = ii + 4;
	      
	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */
  
  loop(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  /* Initial position values are between -0.5 and 0.5 */
	  rx[a][i] = rx[a][i] - 0.5 * L; 
	  ry[a][i] = ry[a][i] - 0.5 * L;
	  rz[a][i] = rz[a][i] - 0.5 * L;
	}
    }
  return;
  /* DEBUG !!!!!!!!!!!! */
  loop(i, 1, Nm)
    {
      dist = sqrt( Sqr(rx[0][i] - rx[1][i]) + Sqr(ry[0][i] - ry[1][i]) + 
		 Sqr(rz[0][i] - rz[1][i]) );
      printf("dist: %f\n", dist);
    }
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

void resetCM(int Nm)
{
  COORD_TYPE sumx, sumy, sumz, RCMx, RCMy, RCMz, Rx, Ry, Rz;
  COORD_TYPE Vx, Vy, Vz, *m;
  int i, a;

  m = Oparams.m;
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  loop(i, 1, Nm)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 CoMV(i, &Vx, &Vy, &Vz);
	 sumx = sumx + Vx;
	 sumy = sumy + Vy;
	 sumz = sumz + Vz;
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  loopShrC(i, 1, Nm)
    {
      vx[0][i] = vx[0][i] - m[0]*sumx/Mtot;
      vx[1][i] = vx[1][i] - m[1]*sumx/Mtot;
      vy[0][i] = vy[0][i] - m[0]*sumy/Mtot;
      vy[1][i] = vy[1][i] - m[1]*sumy/Mtot;
      vz[0][i] = vz[0][i] - m[0]*sumz/Mtot;
      vz[1][i] = vz[1][i] - m[1]*sumz/Mtot;
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
      CoM(i, &Rx, &Ry, &Rz);
      RCMx += Rx; /* Here RCM is the center of mass of the box */
      RCMy += Ry;
      RCMz += Rz;
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  loopShrC(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  rx[a][i] -= RCMx;
	  ry[a][i] -= RCMy;
	  rz[a][i] -= RCMz;
	}
    }
}

/* ========================== >>> comvel <<< =============================== */
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE m[NA])
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
  COORD_TYPE rTemp, sumx, sumy, sumz, Mtot, RCMx, RCMy, RCMz, Rx, Ry, Rz;
  //COORD_TYPE Px, Py, Pz;
  int i, a;
  
  Mtot = m[0] + m[1];                       /* total mass of molecule */
  rTemp = sqrt(temp / Mtot);  
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
      
      vx[0][i] = rTemp * gauss(); 
      vx[1][i] = vx[0][i];
      vy[0][i] = rTemp * gauss();
      vy[1][i] = vy[0][i];
      vz[0][i] = rTemp * gauss();
      vz[1][i] = vz[0][i];
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
	 sumx = sumx + vx[0][i];
	 sumy = sumy + vy[0][i];
	 sumz = sumz + vz[0][i];
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  loop(i, 1, Nm)
    {
      vx[0][i] = vx[0][i] - sumx;
      vx[1][i] = vx[1][i] - sumx;
      vy[0][i] = vy[0][i] - sumy;
      vy[1][i] = vy[1][i] - sumy;
      vz[0][i] = vz[0][i] - sumz;
      vz[1][i] = vz[1][i] - sumz;
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
      CoM(i, &Rx, &Ry, &Rz);
      RCMx += Rx; /* Here RCM is the center of mass of the box */
      RCMy += Ry;
      RCMz += Rz;
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  loop(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  rx[a][i] -= RCMx;
	  ry[a][i] -= RCMy;
	  rz[a][i] -= RCMz;
	}
    }
  /* Now the center of mass of the box is in the origin */
}

/* ========================== >>> angvel <<< =============================== */
void angvel (int Nm, COORD_TYPE temp, COORD_TYPE m[NA], COORD_TYPE d)
{
  /* Angular velocities from the maxwell-boltzmann distribution.   
    
     The distribution is determined by temperature and inertia. 
     This routine is specific to linear molecules.              
     It chooses the direction of the angular velocity randomly but
     perpendicular to the molecular axis. The square of the       
     magnitude of the angular velocity is chosen from an          
     exponential distribution. There is no attempt to set the     
     total angular momentum to zero.                              
     
     ROUTINE REFERENCED:                                          
     
     COORD_TYPE ranf();                                 
     Returns a uniform random variate on the range zero to one 
  
     VARIABLES:
     
     d       Molecule length 
     temp    Desired temperature
     Nm      Number of molecules

  */ 
  
  COORD_TYPE inert;                 /* momentum of inertia of the molecule */
  
  COORD_TYPE  norm, dot, osq, o, mean;
  COORD_TYPE  xisq, xi1, xi2, xi;
  COORD_TYPE  d01x, d01y, d01z;    /* unit vector parallel to the molecule */
  COORD_TYPE d0x, d0y, d0z, d1x, d1y, d1z;                    /* see later */
  COORD_TYPE d0Sq, d1Sq, Mtot;
  COORD_TYPE oVd0_x, oVd0_y, oVd0_z, oVd1_x, oVd1_y, oVd1_z;
  COORD_TYPE ox, oy, oz;
  COORD_TYPE rab, rabSq, rx0, ry0, rz0;
  COORD_TYPE rx1, ry1, rz1;
  COORD_TYPE rx01, ry01, rz01;
  int i;

//  return;

  L = cbrt(Vol);
  invL = 1.0 / L;
  
  Mtot = m[0] + m[1]; /* total mass of molecule */

  /* square of distance between center of mass and atom 0 */
  d0Sq = Sqr(d * m[1] / Mtot); 
  
  /* square distance between center of mass and atom 1 */
  d1Sq = Sqr(d * m[0] / Mtot);

  inert = m[0] * d0Sq  +  m[1] * d1Sq; /* momentum of inertia */
  
  mean = 2.0 * temp / inert;

  loop(i, 1, Nm)
    {
      rx0 = rx[0][i];
      ry0 = ry[0][i];
      rz0 = rz[0][i];
      rx1 = rx[1][i];
      ry1 = ry[1][i];
      rz1 = rz[1][i];
     
      rx01 = rx0 - rx1;
      rx01 = rx01 - L * rint(invL * rx01); /* minimum image */
      ry01 = ry0 - ry1;
      ry01 = ry01 - L * rint(invL * ry01);
      rz01 = rz0 - rz1;
      rz01 = rz01 - L * rint(invL * rz01);

      rabSq = Sqr(rx01) + Sqr(ry01) + Sqr(rz01);
      rab = sqrt(rabSq);
      //printf("rab: %.15f\n", rab);    
      /* calculate the unit vector parallel to the molecule 
         (from 1 to 0)*/
      d01x = rx01 / rab;
      d01y = ry01 / rab;
      d01z = rz01 / rab;
      
      /* (d0x, d0y, d0z) is the vector joining the center of masse and the atom
	 0, (d1x, d1y, d1z) is the vector joining the center of mass and the
	 atom 1 ( both vectors originate from the CM ) */
      //printf("initCord m0: %.15f m1: %.15f\n", m[0], m[1]);
      d0x =   rx01 * m[1] / Mtot;
      d1x = - rx01 * m[0] / Mtot;
      d0y =   ry01 * m[1] / Mtot;
      d1y = - ry01 * m[0] / Mtot;  ;  
      d0z =   rz01 * m[1] / Mtot;
      d1z = - rz01 * m[0] / Mtot;
      //printf("d0: %.15f d1: %.15f\n", sqrt(Sqr(d0x)+Sqr(d0y)+Sqr(d0z)),
	//     sqrt( Sqr(d1x) + Sqr(d1y) + Sqr(d1z))); 

      /* Choose a random vector on the surface of a unit sphere
	 (see Allen-Tildesley pag.349-350)
         The method is suggested  by Marsaglia and is an improvment of 
         von Neumann method */
      
      xisq = 1.0;
      
      while (xisq >= 1.0)
	{
	  xi1  = ranf() * 2.0 - 1.0;
	  xi2  = ranf() * 2.0 - 1.0;
	  xisq = xi1 * xi1 + xi2 * xi2;
	}
      
      xi = sqrt (fabs(1.0 - xisq));
      ox = 2.0 * xi1 * xi;
      oy = 2.0 * xi2 * xi;
      oz = 1.0 - 2.0 * xisq;
      
      /* Constrain the vector to be perpendicular to the molecule */
      
      dot = ox * d01x + oy * d01y + oz * d01z;
      ox = ox - dot * d01x;
      oy = oy - dot * d01y;
      oz = oz - dot * d01z;
      
      /* Renormalize */
      
      osq   = ox * ox + oy * oy + oz * oz;
      norm  = sqrt(fabs(osq));
      ox    = ox / norm;
      oy    = oy / norm;
      oz    = oz / norm;
      
      /* Choose the magnitude of the angular velocity
         NOTE: consider that it is an exponential distribution 
	 (i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
      
      osq   = - mean * log(ranf());
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
     
      /* Components of vectorial product of angular velocity (ox, oy, oz) and
         (d0x, d0y, d0z) */
      oVd0_x = oy * d0z - oz * d0y; 
      oVd0_y = oz * d0x - ox * d0z;
      oVd0_z = ox * d0y - oy * d0x;

      /* Components of vectorial product of angular velocity (ox, oy, oz) and
         (d1x, d1y, d1z) */
      oVd1_x = oy * d1z - oz * d1y; 
      oVd1_y = oz * d1x - ox * d1z;
      oVd1_z = ox * d1y - oy * d1x;
      
      vx[0][i] = vx[0][i] + oVd0_x;
      vx[1][i] = vx[1][i] + oVd1_x;
      vy[0][i] = vy[0][i] + oVd0_y;
      vy[1][i] = vy[1][i] + oVd1_y;
      vz[0][i] = vz[0][i] + oVd0_z;
      vz[1][i] = vz[1][i] + oVd1_z;
      
    }
} 

/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  setToZero(SAVE_LIST, 
	    NULL);  /* Set to zero all the coordinates */

  FCC(Oparams.parnum, Oparams.d, Oparams.m); 
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
  COORD_TYPE invMtot, m0, m1, r01x, r01y, r01z, v01x, v01y, v01z, d, dSq;
  COORD_TYPE vcmx, vcmy, vcmz;
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
  d = Oparams.d;
  dSq = Sqr(d);
  
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
  
  /* =========== 
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;
  
  loop(i, 1, Oparams.parnum)
    {
      CoM(i, &Rx, &Ry, &Rz);
      RCMx += Rx;
      RCMy += Ry;
      RCMz += Rz;
    }
  printf("Box CoM RCMx: %f RCMy: %f RCMz: %f\n", RCMx, RCMy, RCMz);
  ======== */
  /*
  loop(a, 1, NA)
    {
      loop(b, 1, NA)
	{
	  printf("sigab[%d][%d]: %.10f\n", a, b, Oparams.sigab[a][b]);
	  printf("epsab[%d][%d]: %.10f\n", a, b, Oparams.epsab[a][b]);
	}
    }
  */
  
  /* Allocate arrays of integers, you must supply the array length and the
     address of the pointer to the array */
  head = malloc(sizeof(int)*NCell);
  list = malloc(sizeof(int)*NA*Oparams.parnum);
  map  = malloc(sizeof(int)*mapSize);

  Nm = Oparams.parnum;
  sct = sizeof(COORD_TYPE);
  
  FxL[0] = malloc(sct*Nm);
  FyL[0] = malloc(sct*Nm);
  FzL[0] = malloc(sct*Nm);
  FxL[1] = malloc(sct*Nm);
  FyL[1] = malloc(sct*Nm);
  FzL[1] = malloc(sct*Nm);
  Rmx    = malloc(sct*Nm);
  Rmy    = malloc(sct*Nm);
  Rmz    = malloc(sct*Nm);
  Dphix = malloc(sct*Nm);
  Dphiy = malloc(sct*Nm);
  Dphiz = malloc(sct*Nm);
  
  ox = malloc(sct*Nm);
  oy = malloc(sct*Nm);
  oz = malloc(sct*Nm);

  ux = malloc(sct*Nm);
  uy = malloc(sct*Nm);
  uz = malloc(sct*Nm);
  
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
  m0 = Oparams.m[0];
  m1 = Oparams.m[1];
  m = Oparams.m;

  Mtot = 0.0;
  loop(a, 1, NA)
    {
      Mtot += m[a];
    }
  invMtot = 1.0 / Mtot;

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
      loop(i, 1, Oparams.parnum)
	{
	  /* Calculate angular velocities and orientations 
	     of each particle */
	  r01x = rx[0][i] - rx[1][i];
	  r01y = ry[0][i] - ry[1][i];
	  r01z = rz[0][i] - rz[1][i];
	  
	  /* Store molecular orientations */
	  ux[i] = r01x / d;
	  uy[i] = r01y / d;
	  uz[i] = r01z / d;
	  
	  v01x = vx[0][i] - vx[1][i];
	  v01y = vy[0][i] - vy[1][i];
	  v01z = vz[0][i] - vz[1][i];
	  
	  vectProd(r01x, r01y, r01z, v01x, v01y, v01z, &ox[i], &oy[i], &oz[i]);
	  /* The angular velocities is a global array so now we have the actual
	     angular velocity */
	  ox[i] = ox[i] / dSq;
	  oy[i] = oy[i] / dSq;
	  oz[i] = oz[i] / dSq;
	}
   
      loop(i, 1, Oparams.parnum)
	{
	  
	  /* store the initial positions of particles */
	  OprogStatus.rxCMi[i] = (m0 * rx[0][i] + m1 * rx[1][i]) * invMtot;
	  OprogStatus.ryCMi[i] = (m0 * ry[0][i] + m1 * ry[1][i]) * invMtot;
	  OprogStatus.rzCMi[i] = (m0 * rz[0][i] + m1 * rz[1][i]) * invMtot;
	 
	  /* zero the accumulators that contains the temporal integral of 
	     the angular velocity for each particle */
	  OprogStatus.sumox[i] = 0.0;
	  OprogStatus.sumoy[i] = 0.0;
	  OprogStatus.sumoz[i] = 0.0;
	  
	  OprogStatus.ox0[i] = ox[i];
	  OprogStatus.oy0[i] = oy[i];
	  OprogStatus.oz0[i] = oz[i];
	  
	  OprogStatus.ux0[i] = ux[i];
	  OprogStatus.uy0[i] = uy[i];
	  OprogStatus.uz0[i] = uz[i];
	  
	  /* Center of mass velocities */
	  vcmx = (m0 * vx[0][i] + m1 * vx[0][i]) * invMtot;
	  vcmy = (m0 * vy[0][i] + m1 * vy[0][i]) * invMtot;
	  vcmz = (m0 * vz[0][i] + m1 * vz[0][i]) * invMtot;
	
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
  printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);
}
