#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
#define MD_DEBUG(x) x
extern char TXT[MSG_LEN];
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
void setToZero(COORD_TYPE* ptr, ...);
double maxax[2];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
double mgA, mgB, g2;
extern COORD_TYPE pi, L, invL, L2;   
#ifdef MD_GRAVITY
extern double Lz2;
#endif
extern COORD_TYPE W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB; 

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern const double timbig;
double **Xa, **Xb, **RA, **RB, ***R, **Rt;
#ifdef MD_ASYM_ITENS
double **Ia, **Ib, **invIa, **invIb;
#else
double Ia, Ib, invIa, invIb;
#endif
extern double **matrix(int n, int m);
int poolSize;
int parnumA, parnumB;
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
int *bondscache, *numbonds, **bonds, *numbonds0, **bonds0;
#endif
double invaSq[2], invbSq[2], invcSq[2];
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

void ScheduleEvent (int idA, int idB, double tEvent); 

/* ============================= >>> FCC <<< ================================*/
void FCC(int Nm, COORD_TYPE *m)
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
  int Nc, Ncz;
  double rRoot3; /* = 0.5773503; */
  double Cell, Cell2;
#ifdef MD_GRAVITY
  double Cellz, Cell2z;
#endif
  int i, ix, iy, iz, iref, ii;
  double bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  
#ifdef MD_GRAVITY
  /* NOTA:
     Ncz è tanto più grande quanto maggiore è Lz  rispetto a L,
     cioè così si tiene conto che la scatola può avere un'altezza 
     differente dagli altri due lati */
  Nc = ceil( pow ( (L/Lz)*((double)Nm) / 4.0, 1.0/3.0) ); 
  Ncz = ceil((Lz/L)*Nc);	      
  printf("Nc: %d Ncz:%d\n", Nc, Ncz);
#else
  Ncz = Nc = ceil(  pow( ((double)Nm)/4.0, 1.0/3.0 )  );
  printf("Nc: %d\n", Nc);
#endif
  /* Calculate the side of the unit cell */
  Cell  = L / ((double) Nc); /* unit cell length */
  Cell2 = 0.5 * Cell;              /* half unit cell length */
#ifdef MD_GRAVITY
  Cellz = Lz / ((double) Nc);
  Cell2z = 0.5 * Cellz;
#endif
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
#ifdef MD_GRAVITY
  bz[2] = Cell2z;
#else
  bz[2] = Cell2;
#endif  
  /* Sublattice D */
  bx[3] =  Cell2;
  by[3] =  0.0;
#ifdef MD_GRAVITY
  bz[3] = Cell2z;
#else
  bz[3] =  Cell2;
#endif  
  
  /* Construct the lattice from the unit cell */
  
  ii = 0;
  for(iz = 0; iz < Ncz; iz++) /* loops over unit cells (that are simply cubes) */ 
    {
      for(iy = 0; iy < Nc; iy++)
	{
	  for(ix = 0; ix < Nc; ix++)
	    {
	      for(iref = 0; iref < 4; iref++) /* In each primitive cell there are four 
						 molecules */
		{
		  if ((ii + iref) >= Nm) break;
		  /* If Nm < 4 * Nc^3 the we have more lattice sites than 
		     particles so we stop if all the particles was positioned.
		     This condition in fact means: 'If all the particles was
		     positioned => end'
		  */
		  
		  /* Center of Mass of the actual molecule (m + iref) */
		  rx[ii+iref] = bx[iref] + Cell * ((double) ix);
		  ry[ii+iref] = by[iref] + Cell * ((double) iy);
#ifdef MD_GRAVITY
		  rz[ii+iref] = bz[iref] + Cellz * ((double) iz);
#else
		  rz[ii+iref] = bz[iref] + Cell * ((double) iz);
#endif
#if 0
		  printf("#%d (%f,%f,%f) ix:%d iy:%d iz:%d\n", ii+iref, rx[ii+iref],
			 ry[ii+iref], rz[ii+iref], ix, iy, iz);
#endif
		}
	      ii = ii + 4;
	    }
	}
    }
  
  /* Shift centre of box to the origin */
  
  for(i = 0;i < Nm; i++)
    {
      /* Initial position values are between -0.5 and 0.5 */
      rx[i] = rx[i] - 0.5 * L; 
      ry[i] = ry[i] - 0.5 * L;
#ifdef MD_GRAVITY
      if (i < Oparams.parnumA)
	rz[i] = rz[i] - 0.5 * Lz + Oparams.sigma[0][0]*0.5 + 0.1;
      else
	rz[i] = rz[i] - 0.5 * Lz + Oparams.sigma[1][1]*0.5 + 0.1;
#else
      rz[i] = rz[i] - 0.5 * L;
#endif
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

void resetCM(int onlyz)
{
  COORD_TYPE sumx, sumy, sumz, RCMx, RCMy, RCMz;
  COORD_TYPE Mtot;
  int i;
  /* Remove net momentum, to have a total momentum equals to zero */
  Mtot = Oparams.parnumA*Oparams.m[0]+(Oparams.parnum-Oparams.parnumA)*Oparams.m[1];
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  for(i=0; i < Oparams.parnumA; i++)
    {
      sumx = sumx + vx[i]*Oparams.m[0];
      sumy = sumy + vy[i]*Oparams.m[0];
      sumz = sumz + vz[i]*Oparams.m[0];
    }
  for(i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      sumx = sumx + vx[i]*Oparams.m[1];
      sumy = sumy + vy[i]*Oparams.m[1];
      sumz = sumz + vz[i]*Oparams.m[1];
    }
     
  sumx = sumx / Mtot; 
  sumy = sumy / Mtot;
  sumz = sumz / Mtot;

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i = 0; i < Oparams.parnum; i++)
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

  for(i = 0; i < Oparams.parnumA; i++)
    {
      RCMx += rx[i]*Oparams.m[0]; /* Here RCM is the center of mass of the box */
      RCMy += ry[i]*Oparams.m[0];
      RCMz += rz[i]*Oparams.m[0];
    }
  
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      RCMx += rx[i]*Oparams.m[1]; /* Here RCM is the center of mass of the box */
      RCMy += ry[i]*Oparams.m[1];
      RCMz += rz[i]*Oparams.m[1];
    }
  
  RCMx /= Mtot;
  RCMy /= Mtot;
  RCMz /= Mtot;

  for(i = 0; i < Oparams.parnum; i++)
    {
      if (onlyz)
	 rz[i] -= RCMz;
      else
	{
       	  rx[i] -= RCMx;
	  ry[i] -= RCMy;
	  rz[i] -= RCMz;
	}
    }
}
void comvel_brown (COORD_TYPE temp, COORD_TYPE *m)
{
  COORD_TYPE rTemp[2] ;
  /*COORD_TYPE Px, Py, Pz;*/
  int i;
  rTemp[0] = sqrt(temp / m[0]);  
  rTemp[1] = sqrt(temp / m[1]);
  /* variance of the velocities distribution function, we assume k = 1 */ 
  for (i = 0; i < Oparams.parnumA; i++)
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
      
      vx[i] = rTemp[0] * gauss(); 
      vy[i] = rTemp[0] * gauss();
      vz[i] = rTemp[0] * gauss();
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
    }
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      vx[i] = rTemp[1] * gauss(); 
      vy[i] = rTemp[1] * gauss();
      vz[i] = rTemp[1] * gauss();
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
    }
 
}

/* ========================== >>> comvel <<< =============================== */
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE *m, int resetCM)
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
  COORD_TYPE rTemp[2], sumx, sumy, sumz, RCMx, RCMy, RCMz, Mtot;
  /*COORD_TYPE Px, Py, Pz;*/
  int i;
  Mtot = Oparams.parnumA*Oparams.m[0]+(Oparams.parnum-Oparams.parnumA)*Oparams.m[1];
  rTemp[0] = sqrt(temp / m[0]);  
  rTemp[1] = sqrt(temp / m[1]);
  /* variance of the velocities distribution function, we assume k = 1 */ 
  K = 0;
  for (i = 0; i < Oparams.parnumA; i++)
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
      
      vx[i] = rTemp[0] * gauss(); 
      vy[i] = rTemp[0] * gauss();
      vz[i] = rTemp[0] * gauss();
      //printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
      K = K + 0.5 * Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));
    }
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
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
      
      vx[i] = rTemp[1] * gauss(); 
      vy[i] = rTemp[1] * gauss();
      vz[i] = rTemp[1] * gauss();
      //printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
      K = K + 0.5 * Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));
    }
 
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  for(i = 0; i < Oparams.parnumA; i++)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 sumx = sumx + vx[i]*Oparams.m[0];
	 sumy = sumy + vy[i]*Oparams.m[0];
	 sumz = sumz + vz[i]*Oparams.m[0];
	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
   for(i = Oparams.parnumA; i < Oparams.parnum; i++)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 sumx = sumx + vx[i]*Oparams.m[1];
	 sumy = sumy + vy[i]*Oparams.m[1];
	 sumz = sumz + vz[i]*Oparams.m[1];
	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
 
  sumx = sumx / Mtot; 
  sumy = sumy / Mtot;
  sumz = sumz / Mtot;

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i = 0; i < Oparams.parnum; i++)
    {
      vx[i] = vx[i] - sumx;
      vy[i] = vy[i] - sumy;
      vz[i] = vz[i] - sumz;
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }


  if (!resetCM)
    return;
#ifndef MD_GRAVITY
  printf("temp: %f T: %f\n", temp, 2.0*K/(3.0*Oparams.parnum - 3.0));
  scalevels(temp, K);
#endif
  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  for(i = 0; i < Oparams.parnumA; i++)
    {
      RCMx += rx[i]*Oparams.m[0]; /* Here RCM is the center of mass of the box */
      RCMy += ry[i]*Oparams.m[0];
      RCMz += rz[i]*Oparams.m[0];
    }
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      RCMx += rx[i]*Oparams.m[1]; /* Here RCM is the center of mass of the box */
      RCMy += ry[i]*Oparams.m[1];
      RCMz += rz[i]*Oparams.m[1];
    }
  

  RCMx /= Mtot;
  RCMy /= Mtot;
  RCMz /= Mtot;
  for(i=0; i < Oparams.parnum; i++)
    {
      //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
      rx[i] -= RCMx;
      ry[i] -= RCMy;
#ifndef MD_GRAVITY
      rz[i] -= RCMz;
#endif
    }
  /* Now the center of mass of the box is in the origin */
}
void angvel(void)
{
  int i;

  for (i=0; i < Oparams.parnum; i++)
    {
      uxx[i] = 1.0;
      uyx[i] = 0.0;
      uzx[i] = 0.0;
      uxy[i] = 0.0;
      uyy[i] = 1.0;
      uzy[i] = 0.0;
      uxz[i] = 0.0;
      uyz[i] = 0.0;
      uzz[i] = 1.0;
      wx[i] = 0;
      wy[i] = 0;
      wz[i] = 0;
    }
  
}
void wrap_initCoord(void)
{
  /* A x(-0.603750000000000,4.226250000000000,-0.805000000000000) v(-0.099616130522196,-1.839280599669232,0.357754947051492f)-B x(-2.616250000000000,2.213750000000000,-0.805000000000000) v(1.011838511395152,0.876050550528104,-0.426995365917961)
   * */
  rx[0] = -2.0;
  ry[0] = 0.0;
  rz[0] =  0;
  
  vx[0] = 0.1;
  vy[0] = 0;
  vz[0] = 0;
  wx[0] = 1.0;
  wy[0] = -1.0;
  wz[0] = 0.5;
#if 0
  uxx[0] = 0.707;
  uyx[0] = -0.707;
  uzx[0] = 0.0;
  uxy[0] = 0.707;
  uyy[0] = 0.707;
  uzy[0] = 0.0;
#endif
  uxz[0] = 0.0;
  uyz[0] = 0.0;
  uzz[0] = 1.0;
  rx[1] = 2.0;
  ry[1] = 0.0;
  rz[1] = 0;
  vx[1] = -0.1;
  vy[1] = 0;
  vz[1] = 0;
  wx[1] =-2;
  wy[1] = -1;
  wz[1] = 0.1;
}

void adjust_norm(double **R);
/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  int i;
  setToZero(SAVE_LIST, 
	    NULL);  /* Set to zero all the coordinates */

  FCC(Oparams.parnum, Oparams.m); 
  /* Put the baricenter of each molecule on a FCC lattice, and set 
     their orientations */  
  
  /* set both atoms velocity to the center of mass velocity */
  comvel(Oparams.parnum, Oparams.T, Oparams.m, 1); 
#if 1
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
  printf("All'inizio T=%f\n", 2.0 * K / (3.0 * Oparams.parnum - 3.0));

#endif
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
  angvel(); 
  //wrap_initCoord();
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

    V = 0.0;
    L = 9.4;
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
    Oparams.delta[0][0] = Oparams.delta[1][1] = Oparams.delta[0][1] = Oparams.delta[1][0] = 0.0;
    Oparams.bheight = 0.0;
    OprogStatus.maxbonds = 20;
#endif
#ifdef MD_GRAVITY
    Lz = 9.4;
#endif
    Oparams.T = 2.0;
    Oparams.P = 1.0;
    Oparams.M = 5; /* cells in each direction for linked lists */
    
    Oparams.time = 0.0;
    OprogStatus.tolT = 0.0;
    OprogStatus.nebrTabFac = 150;
    OprogStatus.rNebrShell = 0.4;
    /* If 1 the program calculate of the corrisponding variable a mean from
       the begin of the run and not the instanteaneous value */
    OprogStatus.avnggr    = 0;
    Oparams.Dt = 0.01;
    OprogStatus.avngS     = 0;
    OprogStatus.avngPress = 0;
    OprogStatus.avngTemp  = 0;
    OprogStatus.scalevel = 0;
    OprogStatus.endtime = 0;
    OprogStatus.rescaleTime = 1.0;
    OprogStatus.brownian = 0;
#ifdef MD_GRAVITY
    OprogStatus.taptau = 0.0;
    OprogStatus.rzup = 0.0;
    OprogStatus.expandFact= 1.0;
    OprogStatus.quenchend = 0.0;
    OprogStatus.tc = 1E-5;
    OprogStatus.accrcmz = 0.0;
    OprogStatus.wallcollCount = 0;
    OprogStatus.checkquenchTime = 1.0;
    OprogStatus.numquench = 0;
    OprogStatus.maxquench = 0;
    OprogStatus.rescaleTime = 0.0;
    OprogStatus.extraLz = 10.0;
    OprogStatus.rhobh = 0.0;
    OprogStatus.vztap = 10.0;
#endif
#ifndef MD_ASYM_ITENS
    for (i = 0; i < 2; i++)
      Oparams.I[i] = 1.0;
#endif
    OprogStatus.eventMult = 100;
    OprogStatus.overlaptol = 0.0001;
    /* Il promo step inizia con un tapping a temperatura T */
    Oparams.m[0] = Oparams.m[1] = 1.0;
    //Oparams.sigma[0][0] = Oparams.sigma[1][1] = Oparams.sigma[1][0]= Oparams.sigma[0][1]=1.0;
    OprogStatus.collCount = 0;
    OprogStatus.crossCount = 0;
    OprogStatus.nextSumTime = 0.0;
    OprogStatus.nextcheckTime = 0.0;
    OprogStatus.intervalSum = 1.0;
    OprogStatus.storerate = 0.01;
    OprogStatus.KK = 0;
    OprogStatus.JJ = 0;
    /* Parameters relative to Ruocco AC force
       See: cond-mat/00001311, Ruocco et al. */
    OprogStatus.rNebrShell = 2.7; /* the radius of neighbour list shell */
    for (i = 0; i < PE_POINTS; i++)
      OprogStatus.PE[i] = 0;
    /* ======================================================================= */
  }
extern void check (int *overlap, double *K, double *V);
double *atomTime, *treeTime, *treeRxC, *treeRyC, *treeRzC;
int **tree, *inCell[3], *cellList, cellsx, cellsy, cellsz, cellRange[2*NDIM];
extern void PredictEvent(int, int);
extern void InitEventList(void);
void StartRun(void)
  {
    int j, k, n;

    for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
      cellList[j] = -1;
    /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
    for (n = 0; n < Oparams.parnum; n++)
      {
	atomTime[n] = Oparams.time;
	inCell[0][n] =  (rx[n] + L2) * cellsx / L;
	inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
	inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
	inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
#endif
#if 0
	if (inCell[0][n]>=cellsx ||inCell[1][n]>= cellsy||inCell[2][n]>= cellsz) 
	  {
	    printf("BOH?!?L:%f L2:%f n:%d rx[n]:%f\n", L, L2, n, rx[n]);
	    printf("(%d,%d,%d) (%d,%d,%d)\n",cellsx , cellsy,cellsz,
		   inCell[0][n],inCell[1][n], inCell[2][n]);
	  }
#endif	  
	j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	  inCell[0][n] + Oparams.parnum;
	cellList[n] = cellList[j];
	cellList[j] = n;
      }
    InitEventList();
    for (k = 0;  k < NDIM; k++)
      {
	cellRange[2*k]   = - 1;
	cellRange[2*k+1] =   1;
      }
    for (n = 0; n < Oparams.parnum; n++)
      PredictEvent(n, -2); 
  }

#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
extern void add_bond(int na, int n);
extern void remove_bond(int na, int n);
extern double calcpotene(void);
#endif
#if defined(MD_SQWELL) && defined(MD_BONDCORR) 
double corrini3, corrini0, corrini1, corrini2, corrnorm;
double *lastbreak1, *lastbreak2;
#endif
extern void print_matrix(double **M, int n);

void u2R(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      R[i][0][0] = uxx[i];
      R[i][0][1] = uxy[i];
      R[i][0][2] = uxz[i];
      R[i][1][0] = uyx[i];
      R[i][1][1] = uyy[i];
      R[i][1][2] = uyz[i];
      R[i][2][0] = uzx[i];
      R[i][2][1] = uzy[i];
      R[i][2][2] = uzz[i];
      MD_DEBUG2(print_matrix(R[i], 3));
    }
  
}
void R2u(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      uxx[i] = R[i][0][0];
      uxy[i] = R[i][0][1];
      uxz[i] = R[i][0][2];
      uyx[i] = R[i][1][0];
      uyy[i] = R[i][1][1];
      uyz[i] = R[i][1][2];
      uzx[i] = R[i][2][0];
      uzy[i] = R[i][2][1];
      uzz[i] = R[i][2][2];
    }
}
/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
  {
    /* DESCRIPTION:
       This function is called after the parameters were read from disk, put
       here all initialization that depends upon such parameters, and call 
       all your function for initialization, like maps() in this case */
#if defined(MD_SQWELL) && defined(MD_BONDCORR) 
    char fileop[1024], fileop2[1024], fileop3[1024];
    FILE *bof;
#endif    
    int Nm, i, sct, overlap;
    COORD_TYPE vcmx, vcmy, vcmz;
    COORD_TYPE *m;
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
    double sigDeltaSq, drx, dry, drz;
    int j;
#endif
    int a;
    /*COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;*/

    /* initialize global varibales */
    pi = 2.0 * acos(0);
    
    Nm = Oparams.parnumA;
    parnumA = Oparams.parnumA;
    parnumB = Oparams.parnum - Oparams.parnumA;
    sct = sizeof(COORD_TYPE);

    invL = 1.0/L;
    L2 = 0.5*L;
#ifdef MD_GRAVITY
    rcmz = -Lz*0.5;
    Lz2 = Lz*0.5;
#endif
    poolSize = OprogStatus.eventMult*Oparams.parnum;
    m = Oparams.m;
    Mtot = Oparams.m[0]*parnumA+Oparams.m[1]*parnumB;
    invmA = 1.0/Oparams.m[0];
    invmB = 1.0/Oparams.m[1];
#if 0
    Oparams.sigma[1][0] = Oparams.sigma[0][1];
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
    Oparams.delta[1][0] = Oparams.delta[0][1];
#endif
#endif
    Mred[0][0] = Mred[1][1] = 0.5;
    Mred[0][1] = Mred[1][0] = (Oparams.m[0]*Oparams.m[1])/(Oparams.m[0]+Oparams.m[1]);
    /* Calcoliam rcut assumendo che si abbian tante celle quante sono 
     * le particelle */
    if (Oparams.rcut <= 0.0)
      Oparams.rcut = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
    cellsx = L / Oparams.rcut;
    cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
    cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
    cellsz = L / Oparams.rcut;
#endif 
    printf("Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", Oparams.rcut,
	   cellsx, cellsy, cellsz);
#if 0
    printf("massA: %f massB: %f sigmaA:%f sigmaB:%f sigmaAB:%f\n", Oparams.m[0], Oparams.m[1],
	 Oparams.sigma[0][0], Oparams.sigma[1][1], Oparams.sigma[0][1]);
#endif
#if defined(MD_GRAVITY)
    g2 = 0.5*Oparams.ggrav;
    mgA = Oparams.m[0]*Oparams.ggrav; 
    mgB = Oparams.m[1]*Oparams.ggrav;
#endif
   
  /*    
     ** CHECK FOR PARTICLE OVERLAPS **
     ** CALCULATE ENERGY            ** */
    /*lastcol= malloc(sizeof(double)*Oparams.parnum);*/
    atomTime = malloc(sizeof(double)*Oparams.parnum);
    cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
    inCell[0] = malloc(sizeof(int)*Oparams.parnum);
    inCell[1]= malloc(sizeof(int)*Oparams.parnum);
    inCell[2] = malloc(sizeof(int)*Oparams.parnum);
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
    tree = AllocMatI(10, poolSize);
    bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
    bonds0 = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
    numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
    numbonds0 = (int *) malloc(Oparams.parnum*sizeof(int));
#ifdef MD_BONDCORR
    lastbreak1=(double*)malloc(Oparams.parnum*sizeof(double));
    lastbreak2=(double*)malloc(Oparams.parnum*sizeof(double));
#endif
    bondscache = (int *) malloc(sizeof(int)*OprogStatus.maxbonds);
#else
    tree = AllocMatI(9, poolSize);
#endif
    treeTime = malloc(sizeof(double)*poolSize);
    treeRxC  = malloc(sizeof(double)*poolSize);
    treeRyC  = malloc(sizeof(double)*poolSize);
    treeRzC  = malloc(sizeof(double)*poolSize);
    Xa = matrix(3, 3);
    Xb = matrix(3, 3);
#ifdef MD_ASYM_ITENS
    Ia = matrix(3, 3);
    Ib = matrix(3, 3);
    invIa = matrix(3, 3);
    invIb = matrix(3, 3);
#endif
    RA = matrix(3, 3);
    RB = matrix(3, 3);
    Rt = matrix(3, 3);
    R = malloc(sizeof(double**)*Oparams.parnum);
    for (i=0; i < Oparams.parnum; i++)
      {
	R[i] = matrix(3, 3);
      }
    u2R();
    if (OprogStatus.CMreset==-1)
      {
	comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
	resetCM(1);
      }
    else if (OprogStatus.CMreset==-2)
      {
	comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
      }

    if (Oparams.curStep == 1)
    {
      check (&overlap, &K, &V);
     
      if ( overlap ) 
     	{
 	  printf("ERROR: Particle overlap in initial configuration\n");
 	  exit(1);      
  	}
    }
#if 0
  for (i= 0; i < Oparams.parnum; i++)
    {
      /*printf("i=%d v(%.15G,%.15G,%.15G)\n", i, vx[i], vy[i], vz[i]);*/
      lastcol[i] = 0.0;
    }
#endif
  
  if (newSim)
    {
      Oparams.time=0.0;
      for (i=0; i < Oparams.parnum; i++)
	atomTime[i] = 0.0;
      OprogStatus.nextcheckTime += fabs(OprogStatus.rescaleTime);
      OprogStatus.nextSumTime += OprogStatus.intervalSum;
      if (OprogStatus.storerate > 0.0)
      	OprogStatus.nextStoreTime = OprogStatus.storerate;
      OprogStatus.nextDt += Oparams.Dt;
    }
  else
    {
      for (i=0; i < Oparams.parnum; i++)
	atomTime[i] = Oparams.time;
    }
  /* evaluation of principal inertia moments*/ 
  for (a = 0; a < 2; a++)
    {
      invaSq[a] = Sqr(1/Oparams.a[a]);
      invbSq[a] = Sqr(1/Oparams.b[a]);
      invcSq[a] = Sqr(1/Oparams.c[a]);
#ifdef MD_ASYM_ITENS
      ItensD[a][0] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.b[a])+Sqr(Oparams.c[a]));
      ItensD[a][1] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.a[a])+Sqr(Oparams.c[a]));
      ItensD[a][2] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.a[a])+Sqr(Oparams.b[a]));
#endif
    };
 
  /* maxax è il diametro del centroide */
  for (a = 0; a < 2; a++)
    {
      maxax[a] = 0.0;
      if (Oparams.a[a] > maxax[a])
	maxax[a] = Oparams.a[a];
      if (Oparams.b[a] > maxax[a])
	maxax[a] = Oparams.b[a];
      if (Oparams.c[a] > maxax[a])
	maxax[a] = Oparams.c[a];
      maxax[a] *= 2.0;
    }
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
#ifdef MD_BONDCORR
      lastbreak1[i] = 0.0;
      lastbreak2[i] = 0.0;
#endif
    }
  for ( i = 0; i < Oparams.parnum-1; i++)
    for ( j = i + 1; j < Oparams.parnum; j++)
      {
	if (i < parnumA && j < parnumA)
	  {
	    sigDeltaSq = Sqr(Oparams.sigma[0][0]+Oparams.delta[0][0]);
	  }
	else if (i >= parnumA && j >= parnumA)
	  {
	    sigDeltaSq = Sqr(Oparams.sigma[1][1]+Oparams.delta[1][1]);
	  }
       	else
	  {
	    sigDeltaSq = Sqr(Oparams.sigma[0][1]+Oparams.delta[0][1]);
	  }
	drx = rx[i] - rx[j];
	dry = ry[i] - ry[j];
	drz = rz[i] - rz[j];
	drx = drx - L * rint(drx / L);
	dry = dry - L * rint(dry / L);
	drz = drz - L * rint(drz / L);
	if (Sqr(drx)+Sqr(dry)+Sqr(drz) < sigDeltaSq) 
	  {
	    add_bond(i, j);
	    add_bond(j, i);
	  }
      }
#ifdef MD_BONDCORR
  corrnorm=0;
  for (i=0; i < Oparams.parnum; i++)
    { 
      numbonds0[i]=numbonds[i];
      corrnorm+=numbonds[i];
      for (j=0; j < numbonds0[i]; j++)
	bonds0[i][j]=bonds[i][j];
    }
  corrini3 = 0;
  corrini2 = 0;
  corrini1 = 0;
  corrini0 = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (numbonds0[i]==2)
	{
	  corrini2++;
	}
      if (numbonds0[i]==1 && numbonds0[bonds0[i][0]]==1)
	{
	  corrini1++;
	}
      if (numbonds0[i]==0)
	{
	  corrini0++;
	}
      if (numbonds0[i]==3)
	{
	  corrini3++;
	}
    }
  printf("------------> corrini0: %f\n", corrini0);
  printf("------------> corrini1: %f\n", corrini1);
  printf("------------> corrini2: %f\n", corrini2);
  printf("------------> corrini3: %f\n", corrini3);
  printf("------------> corrnorm: %f\n", corrnorm);
  sprintf(fileop2 ,"bondcorr.dat");
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bof = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  fclose(bof);
#if 1
 sprintf(fileop2 ,"BondCorrFuncB2.dat");
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bof = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  //fprintf(bof,"%.15f %.15f\n", Oparams.time+1E-5, corrini1);
  fclose(bof);
  sprintf(fileop2 ,"BondCorrFuncB3.dat");
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bof = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  //fprintf(bof,"%.15f %.15f\n", Oparams.time+1E-5, corrini2);
  fclose(bof);
#endif
#endif
  printf("Energia potenziale all'inizio: %.15f\n", calcpotene());
#endif
  StartRun(); 
  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      FILE *f;
      /* truncate file to zero lenght */
#ifdef MD_GRAVITY
      f = fopenMPI(MD_HD_MIS "T.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "Vz2.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rcmz.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rho.dat", "w");
      fclose(f);
#else
      f = fopenMPI(MD_HD_MIS "T.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "D.dat", "w");
      fclose(f);
#endif
      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;
#ifdef MD_GRAVITY
      if (!strcmp(OprogStatus.inifile,"*"))
	{
	  for (i= 0; i < Oparams.parnum; i++)
	    {
	      /*printf("i=%d v(%.15G,%.15G,%.15G)\n", i, vx[i], vy[i], vz[i]);*/
	      lastcol[i] = 0.0;
	    }
	}
#endif
  
      for(i = 0; i < Oparams.parnum; i++)
	{
	  /* store the initial positions of particles */
	  OprogStatus.rxCMi[i] = rx[i];
	  OprogStatus.ryCMi[i] = ry[i];
	  OprogStatus.rzCMi[i] = rz[i];	 
	  /* Center of mass velocities */
	  /*printf("(%f,%f,%f)\n", rx[i], ry[i], rx[i]);*/
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
      
      for(i = 0; i < NUMK; i++) 
	{
	  OprogStatus.sumS[i] = 0.0;
	}
      
      for(i = 0; i < MAXBIN; i++)
	{
	  OprogStatus.hist[i] = 0;
	}
    }
  /* printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);*/
}


/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i;
  const char tipodat[] = "%.15G %.15G %.15G %.15G %.15G %.15G\n";
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n";
   
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat2,rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	      uyz[i], uzx[i], uzy[i], uzz[i]);
    }
#ifndef MD_STOREMGL 
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat, vx[i], vy[i], vz[i], wx[i], wy[i], wz[i]);
    }
#ifdef MD_GRAVITY
  fprintf(fs, "%.15G %.15G\n", L, Lz);
#else
  fprintf(fs, "%.15G\n", L);
#endif
#endif
}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i;

  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf ", &rx[i], &ry[i], &rz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
      if (fscanf(fs, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		 &uxx[i], &uxy[i], &uxz[i], &uyx[i], &uyy[i], &uyz[i], &uzx[i], &uzy[i], &uzz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
  
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf ", &vx[i], &vy[i], &vz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
      if (fscanf(fs, "%lf %lf %lf\n", &wx[i], &wy[i], &wz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
  

#ifdef MD_GRAVITY
  if (fscanf(fs, "%lf %lf\n",  &L, &Lz) < 2)
    {
      mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
      exit(-1);
    }
#else
  if (fscanf(fs, "%lf\n",  &L) < 1)
    {
      mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
      exit(-1);
    }
#endif
      
}
