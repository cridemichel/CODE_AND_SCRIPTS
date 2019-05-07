#include<mdsimul.h>

#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
NOTE: The box edge length is unity, so every length must be referred to 
the this quantity.
 */
/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
extern double max(double a, double b);
#define MD_DEBUG10(x) 
extern char TXT[MSG_LEN];
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
void setToZero(COORD_TYPE* ptr, ...);
double **radat, **deltat, *maxax;
extern struct LastBumpS *lastbump;
extern double *lastcol;
double *axa,*axb,*axc;
double *a0I;
double **Aip;
#ifdef MD_CALENDAR_HYBRID
extern int *linearLists;
extern int numevPQ, totevHQ, overevHQ;
#endif
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
extern int mapbondsaAB[MD_PBONDS_AB];
extern int mapbondsbAB[MD_PBONDS_AB];
extern int mapbondsaAA[MD_PBONDS_AA];
extern int mapbondsbAA[MD_PBONDS_AA];
extern int mapbondsaBB[MD_PBONDS_BB];
extern int mapbondsbBB[MD_PBONDS_BB];
#elif MD_AB41
extern int mapbondsaAB[MD_PBONDS_AB];
extern int mapbondsbAB[MD_PBONDS_AB];
extern int mapbondsaAA[MD_PBONDS_AA];
extern int mapbondsbAA[MD_PBONDS_AA];
#else
extern int mapbondsaSiO[MD_PBONDS];
extern int mapbondsbSiO[MD_PBONDS];
#endif
extern int *mapbondsa;
extern int *mapbondsb;
int *crossevtodel;
#else
extern int mapbondsa[MD_PBONDS];
extern int mapbondsb[MD_PBONDS];
#endif
int *scdone;
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
double **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **powdirs;
#ifdef MD_ASYM_ITENS
extern double **Ia, **Ib, **invIa, **invIb;
#else
extern double Ia, Ib, invIa, invIb;
#endif
extern double **matrix(int n, int m);
int poolSize;
extern int parnumA, parnumB;
int *bondscache, *numbonds, **bonds, *numbonds0, **bonds0;
double invaSq[2], invbSq[2], invcSq[2];
double calcDistNeg(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, 
		   double dists[MD_PBONDS], int bondpair);

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
void check_coord(void)
{
  int i;
  for (i = 0; i < Oparams.parnum; i++)
    if (fabs(rx[i]) > L*0.5 || fabs(ry[i])>L*0.5 || fabs(rz[i]) > L*0.5)
      {
	printf("%d is out of box!\n", i);
	exit(-1);
      }
}
COORD_TYPE ranf(void);

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
      rx[i] = rx[i] - 0.5 * L + ranf()*1E-7; 
      ry[i] = ry[i] - 0.5 * L + ranf()*1E-7;
#ifdef MD_GRAVITY
      if (i < Oparams.parnumA)
	rz[i] = rz[i] - 0.5 * Lz + Oparams.sigma[0][0]*0.5 + 0.1;
      else
	rz[i] = rz[i] - 0.5 * Lz + Oparams.sigma[1][1]*0.5 + 0.1;
#else
      rz[i] = rz[i] - 0.5 * L + ranf()*1E-7;
#endif
      printf("%d = (%f,%f,%f)\n", i, rx[i], ry[i], rz[i]);
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

void resetCM(void)
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
      rx[i] -= RCMx;
      ry[i] -= RCMy;
      rz[i] -= RCMz;
    }
}
void angvel(void);
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
#ifndef MD_DOUBLE_DT
  angvel();
#endif
}
void scalevels(double temp, double K);

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
#ifdef MD_SILICA
/* nella silica l'ossigeno ha solo due sticky points quindi l'inizializzazione
 * è diversa */
#ifdef MD_THREESPOTS
void buildTetrahedras(void)
{
  int i;
  double Kl, Ktr, pi;
  //double Oangle;
  double radius; 
  /* NOTA: i < Oparams.parnumA => O
   *       i >= Oparams.parnumA => Si */
  /* Oxygen */
  //Oangle = acos(0) * 2.0 * 145.8 / 180.0;
  pi = acos(0.0)*2.0;
  Kl = cos(pi/6.0), Ktr = sin(pi/6.0);

  printf("pi=%.15G radius = %f\n", pi, Oparams.sigma[0][1]/2.0);
  for (i=0; i < Oparams.parnumA; i++)
    {
      /* il raggio è quello dell'interazione Si-O */
      radius = Oparams.sigma[0][1]/2.0;
      uxx[i] = -radius;
      uyx[i] = 0.0;
      uzx[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = radius;
      uyy[i] = 0.0;
      uzy[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      /* NOTA: l'ossigeno ha solo due sticky spots */
      uxz[i] = 0.0;
      uyz[i] = 0.0;
      uzz[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
    }

  /* Silicon */
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      /* il raggio è quello dell'interazione Si-O */
      radius = Oparams.sigma[0][1]/2.0;
      uxx[i] = Kl * radius;
      uyx[i] = -Ktr * radius;
      uzx[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = -Kl* radius;
      uyy[i] = -Ktr * radius;
      uzy[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      uxz[i] = 0.0;
      uyz[i] = radius;
      uzz[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
    }
  //printf("dist=%.15G\n", sqrt( Sqr(uxx[i-1]-uxy[i-1]) + Sqr(uyx[i]-uyy[i]) + Sqr(uzx[i]-uzz[i])));
}
#elif MD_AB41
void buildTetrahedras(void)
{
  int i;
  //double Kl, Ktr, pi;
  double Oangle;
  const double Kl = sqrt(8.0/3.0), Kdh = 1.0/3.0, Ktr = sqrt(8.0)/6.0;
  double radius; 
  /* NOTA: i < Oparams.parnumA => O
   *       i >= Oparams.parnumA => Si */
  /* Oxygen */
  //Oangle = acos(0) * 2.0 * 145.8 / 180.0;
  pi = acos(0.0)*2.0;
  //Kl = cos(pi/6.0), Ktr = sin(pi/6.0);
  Oangle = acos(0) * 2.0 * 145.8 / 180.0;

  printf("pi=%.15G radius = %f\n", pi, Oparams.sigma[0][1]/2.0);
  for (i=0; i < Oparams.parnumA; i++)
    {
      /* il raggio è quello dell'interazione Si-O */
      radius = Oparams.sigma[0][0]/2.0;
      uxx[i] = Kl * radius / 2.0;
      uyx[i] = -Ktr * radius;
      uzx[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = -Kl * radius / 2.0;
      uyy[i] = -Ktr * radius;
      uzy[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      uxz[i] = 0.0;
      uyz[i] = Ktr * 2.0 * radius;
      uzz[i] = -Kdh * radius;
    }

  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      /* il raggio è quello dell'interazione Si-O */
      radius = Oparams.sigma[1][1]/2.0;
      uxx[i] = radius * cos(Oangle);
      uyx[i] = radius * sin(Oangle);
      uzx[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = 0.0;
      uyy[i] = 0.0;
      uzy[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      /* NOTA: l'ossigeno ha solo due sticky spots */
      uxz[i] = 0.0;
      uyz[i] = 0.0;
      uzz[i] = 0.0;
    }
  //printf("dist=%.15G\n", sqrt( Sqr(uxx[i-1]-uxy[i-1]) + Sqr(uyx[i]-uyy[i]) + Sqr(uzx[i]-uzz[i])));
}
#else
void buildTetrahedras(void)
{
  int i;
  const double Kl = sqrt(8.0/3.0), Kdh = 1.0/3.0, Ktr = sqrt(8.0)/6.0;
  double Oangle;
  double radius; 
  /* NOTA: i < Oparams.parnumA => O
   *       i >= Oparams.parnumA => Si */
  /* Oxygen */
  Oangle = acos(0) * 2.0 * 145.8 / 180.0;
  printf("radius = %f\n", Oparams.sigma[0][1]);
  for (i=0; i < Oparams.parnumA; i++)
    {
      /* il raggio è quello dell'interazione Si-O */
      radius = Oparams.sigma[0][1]/2.0;
      uxx[i] = radius * cos(Oangle);
      uyx[i] = radius * sin(Oangle);
      uzx[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = radius;
      uyy[i] = 0.0;
      uzy[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      /* NOTA: l'ossigeno ha solo due sticky spots */
      uxz[i] = 0.0;
      uyz[i] = 0.0;
      uzz[i] = 0.0;
      //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
    }

  /* Silicon */
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      /* il raggio è quello dell'interazione Si-O */
      radius = Oparams.sigma[0][1]/2.0;
      uxx[i] = Kl * radius / 2.0;
      uyx[i] = -Ktr * radius;
      uzx[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = -Kl * radius / 2.0;
      uyy[i] = -Ktr * radius;
      uzy[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      uxz[i] = 0.0;
      uyz[i] = Ktr * 2.0 * radius;
      uzz[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
    }
  //printf("dist=%.15G\n", sqrt( Sqr(uxx[i-1]-uxy[i-1]) + Sqr(uyx[i]-uyy[i]) + Sqr(uzx[i]-uzz[i])));
}
#endif
#if 0
void angvel(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      /* N.B. ora i tre vettori ux uy e uz non sono altro che le coordinate
       * di due hydrogen sites e 1 electron sites. Quindi vanno 
       * scelti in modo da stare su tre spigoli di un tetraedro regolare.
       * In particolare formeranno un triangolo equilatero.
       * Qui il terzo vettore uXz è un electron site. */ 
      wx[i] = 0;
      wy[i] = 0;
      wz[i] = 0;
    }
}
#else
extern double scalProd(double *A, double *B);
extern double calc_norm(double *vec);
extern void vectProdVec(double *A, double *B, double *C);
/* ============================ >>> ranf <<< =============================== */
double ranfRandom(void)
{
  /*  Returns a uniform random variate in the range 0 to 1 (excluding 0).         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return (1.0-(((double)rand()) / (((double) RAND_MAX) + 1)));
}
#if 0
void angvel(void)
{
  int i, a;
  double pi, inert;                 /* momentum of inertia of the molecule */
  double norm, osq, o, mean, symax[3];
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, ww[3], wsz;
#ifdef MD_THREESPOTS
  double u1[3], u2[3], u3[3];
#endif
  //L = cbrt(Vol);
  invL = 1.0 / L;

  Mtot = Oparams.m[0]; /* total mass of molecule */

  inert = Oparams.I[0]; /* momentum of inertia */
  pi = acos(0)*2; 

  mean = sqrt(Oparams.T / inert);
  for (i = 0; i < Oparams.parnumA; i++)
    {
      wx[i] = mean * gauss(); 
      wy[i] = mean * gauss();
      wz[i] = mean * gauss();
    }

  Mtot = Oparams.m[1]; /* total mass of molecule */

  inert = Oparams.I[1]; /* momentum of inertia */

  mean = sqrt(Oparams.T / inert);

  for (i = Oparams.parnumA; i < Oparams.parnumA; i++)
    {
      wx[i] = mean * gauss(); 
      wy[i] = mean * gauss();
      wz[i] = mean * gauss();
    }
}
#else
void angvel(void)
{
  int i, a;
  double pi, inert;                 /* momentum of inertia of the molecule */
  double norm, osq, o, mean, symax[3];
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, ww[3], wsz;
#ifdef MD_THREESPOTS
  double u1[3], u2[3], u3[3];
#endif
  //L = cbrt(Vol);
  invL = 1.0 / L;

  Mtot = Oparams.m[0]; /* total mass of molecule */

  inert = Oparams.I[0]; /* momentum of inertia */
  pi = acos(0)*2; 

#ifdef MD_THREESPOTS
  mean = 2.0*Oparams.T / inert;
#else
  /* N.B. QUESTO E' SBAGLIATO PERCHE' LA DISTRIBUZIONE
     E' COME QUELLA TRASLAZIONALE (VEDI LANDAU) 
     CORREGGERE!!!! */
#if 1
  mean = sqrt(Oparams.T / inert);
#else
  mean = 3.0*Oparams.T / inert;
#endif
#endif
  for (i = 0; i < Oparams.parnumA; i++)
    {
#ifndef MD_THREESPOTS
#if 1
      wx[i] = mean * gauss(); 
      wy[i] = mean * gauss();
      wz[i] = mean * gauss();
#else
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
#if 0
      ww[0] = ox;
      ww[1] = oy;
      ww[2] = oz;
      for (a=0; a < 3; a++)
	symax[a] = R[i][a][0];
      norm = calc_norm(symax);
      for (a=0; a < 3; a++)
	symax[a] /= norm;
      wsz = scalProd(ww, symax);
      ox = ox-symax[0]*wsz;
      oy = oy-symax[1]*wsz;
      oz = oz-symax[2]*wsz;
#endif
      /* Renormalize */
      osq   = ox * ox + oy * oy + oz * oz;
      norm  = sqrt(fabs(osq));
      ox    = ox / norm;
      oy    = oy / norm;
      oz    = oz / norm;

      /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/

#if 1
      osq   = - mean * log(ranfRandom());
#else
      osq   = - mean * log(ranf());
#endif
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      wx[i] = ox;
      wy[i] = oy;
      wz[i] = oz;
#endif
#else
      xi1  = ranf()*2.0*pi;
      ox = cos(xi1);
      oy = sin(xi1);
      /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
#if 1
      osq   = - mean * log(ranfRandom());
#else
      do
	{
	  osq   = - mean * log(ranf());
	}
      while (isnan(osq)||isinf(osq));
#endif
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      for (a=0; a < 3; a++)
	u3[a] = R[i][a][0];
      norm = calc_norm(u3);
      for (a=0; a < 3; a++)
	u3[a] /= norm;
      u2[0] = 1;
      u2[1] = 1;
      u2[2] = 1;
      wsz = scalProd(u2, u3);
      for (a=0; a < 3; a++)
	u2[a] = u2[a]-u3[a]*wsz;
      norm=calc_norm(u2);
      for (a=0; a < 3; a++)
	u2[a] /= norm;
      vectProdVec(u2, u3, u1);
      wx[i] = u1[0]*ox+u2[0]*oy;
      wy[i] = u1[1]*ox+u2[1]*oy;
      wz[i] = u1[2]*ox+u2[2]*oy;
#endif
    }

  Mtot = Oparams.m[1]; /* total mass of molecule */

  inert = Oparams.I[1]; /* momentum of inertia */

#if 1
  mean = sqrt(Oparams.T / inert);
  for (i = Oparams.parnumA; i < Oparams.parnumA; i++)
    {

      wx[i] = mean * gauss(); 
      wy[i] = mean * gauss();
      wz[i] = mean * gauss();
    }
#else
  mean = 3.0*Oparams.T / inert;
  for (i = Oparams.parnumA; i < Oparams.parnumA; i++)
    {
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

      /* Renormalize */
      osq   = ox * ox + oy * oy + oz * oz;
      norm  = sqrt(fabs(osq));
      ox    = ox / norm;
      oy    = oy / norm;
      oz    = oz / norm;

      /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
#if 1 
      osq   = - mean * log(ranfRandom());
#else
      osq   = - mean * log(ranf());
#endif
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      wx[i] = ox;
      wy[i] = oy;
      wz[i] = oz;
    }
#endif
}
#endif
#endif
#else
void buildTetrahedras(void)
{
  int i;
  const double Kl = sqrt(8.0/3.0), Kdh = 1.0/3.0, Ktr = sqrt(8.0)/6.0;
  double radius; 
  for (i=0; i < Oparams.parnum; i++)
    {
      if (i < Oparams.parnumA)
	radius = Oparams.sigma[0][0]/2.0;
      else
	radius = Oparams.sigma[1][1]/2.0;
      uxx[i] = Kl * radius / 2.0;
      uyx[i] = -Ktr * radius;
      uzx[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[red]\n", uxx[i], uyx[i], uzx[i]);
      uxy[i] = -Kl * radius / 2.0;
      uyy[i] = -Ktr * radius;
      uzy[i] = -Kdh * radius;
      //printf("%f %f %f @ 0.075 C[red]\n", uxy[i], uyy[i], uzy[i]);
      uxz[i] = 0.0;
      uyz[i] = Ktr * 2.0 * radius * (MD_DIST_ELECTSITES * 2.0);
      uzz[i] = -Kdh * radius * (MD_DIST_ELECTSITES * 2.0);
      //printf("%f %f %f @ 0.075 C[green]\n", uxz[i], uyz[i], uzz[i]);
    }
}
#if 0
void angvel(void)
{
  int i;
  double inert;                 /* momentum of inertia of the molecule */
  double norm, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz;
  //L = cbrt(Vol);
  invL = 1.0 / L;

  Mtot = Oparams.m[0]; /* total mass of molecule */

  inert = Oparams.I[0]; /* momentum of inertia */

  mean = 3.0*Oparams.T / inert;

  for (i = 0; i < Oparams.parnum; i++)
    {
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
      wx[i] = ox;
      wy[i] = oy;
      wz[i] = oz;
    }
}
#else
void angvel(void)
{
  int i, a;
  double pi, inert;                 /* momentum of inertia of the molecule */
  double norm, osq, o, mean, symax[3];
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, ww[3], wsz;
#ifdef MD_THREESPOTS
  double u1[3], u2[3], u3[3];
#endif
  //L = cbrt(Vol);
  invL = 1.0 / L;

  Mtot = Oparams.m[0]; /* total mass of molecule */

  inert = Oparams.I[0]; /* momentum of inertia */
  pi = acos(0)*2; 

#ifdef MD_THREESPOTS
  mean = 2.0*Oparams.T / inert;
#else
  /* N.B. QUESTO E' SBAGLIATO PERCHE' LA DISTRIBUZIONE
     E' COME QUELLA TRASLAZIONALE (VEDI LANDAU) 
     CORREGGERE!!!! */
#if 1
  mean = sqrt(Oparams.T / inert);
#else
  mean = 3.0*Oparams.T / inert;
#endif
#endif
  for (i = 0; i < Oparams.parnumA; i++)
    {
#ifndef MD_THREESPOTS
#if 1
      wx[i] = mean * gauss(); 
      wy[i] = mean * gauss();
      wz[i] = mean * gauss();
#else
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
#if 0
      ww[0] = ox;
      ww[1] = oy;
      ww[2] = oz;
      for (a=0; a < 3; a++)
	symax[a] = R[i][a][0];
      norm = calc_norm(symax);
      for (a=0; a < 3; a++)
	symax[a] /= norm;
      wsz = scalProd(ww, symax);
      ox = ox-symax[0]*wsz;
      oy = oy-symax[1]*wsz;
      oz = oz-symax[2]*wsz;
#endif
      /* Renormalize */
      osq   = ox * ox + oy * oy + oz * oz;
      norm  = sqrt(fabs(osq));
      ox    = ox / norm;
      oy    = oy / norm;
      oz    = oz / norm;

      /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/

#if 1
      osq   = - mean * log(ranfRandom());
#else
      osq   = - mean * log(ranf());
#endif
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      wx[i] = ox;
      wy[i] = oy;
      wz[i] = oz;
#endif
#else
      xi1  = ranf()*2.0*pi;
      ox = cos(xi1);
      oy = sin(xi1);
      /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
#if 1
      osq   = - mean * log(ranfRandom());
#else
      do
	{
	  osq   = - mean * log(ranf());
	}
      while (isnan(osq)||isinf(osq));
#endif
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      for (a=0; a < 3; a++)
	u3[a] = R[i][a][0];
      norm = calc_norm(u3);
      for (a=0; a < 3; a++)
	u3[a] /= norm;
      u2[0] = 1;
      u2[1] = 1;
      u2[2] = 1;
      wsz = scalProd(u2, u3);
      for (a=0; a < 3; a++)
	u2[a] = u2[a]-u3[a]*wsz;
      norm=calc_norm(u2);
      for (a=0; a < 3; a++)
	u2[a] /= norm;
      vectProdVec(u2, u3, u1);
      wx[i] = u1[0]*ox+u2[0]*oy;
      wy[i] = u1[1]*ox+u2[1]*oy;
      wz[i] = u1[2]*ox+u2[2]*oy;
#endif
    }

  Mtot = Oparams.m[1]; /* total mass of molecule */

  inert = Oparams.I[1]; /* momentum of inertia */

#if 1
  mean = sqrt(Oparams.T / inert);
  for (i = Oparams.parnumA; i < Oparams.parnumA; i++)
    {

      wx[i] = mean * gauss(); 
      wy[i] = mean * gauss();
      wz[i] = mean * gauss();
    }
#else
  mean = 3.0*Oparams.T / inert;
  for (i = Oparams.parnumA; i < Oparams.parnumA; i++)
    {
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

      /* Renormalize */
      osq   = ox * ox + oy * oy + oz * oz;
      norm  = sqrt(fabs(osq));
      ox    = ox / norm;
      oy    = oy / norm;
      oz    = oz / norm;

      /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
#if 1 
      osq   = - mean * log(ranfRandom());
#else
      osq   = - mean * log(ranf());
#endif
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      wx[i] = ox;
      wy[i] = oy;
      wz[i] = oz;
    }
#endif
}
#endif
#endif
void wrap_initCoord(void)
{
  /* A x(-0.603750000000000,4.226250000000000,-0.805000000000000) v(-0.099616130522196,-1.839280599669232,0.357754947051492f)-B x(-2.616250000000000,2.213750000000000,-0.805000000000000) v(1.011838511395152,0.876050550528104,-0.426995365917961)
   * */
  rx[0] = -0.501;
  ry[0] = 0.0;
  rz[0] =  0;

  vx[0] = 0.0;
  vy[0] = 0;
  vz[0] = 0;
  /* -0.285316712824933 -0.182347469854598 -0.530547025349427*/

  wx[0] = 0.0;//-0.285312;// .003;
  wy[0] = 0.0;//-0.1823475;// -1.5;
  wz[0] = 0.0;//-0.530547;// -0.5;

  rx[1] = 0.501;
  ry[1] = 0.0;
  rz[1] = 0.0;
  vx[1] = 0.0;
  vy[1] = 0;
  vz[1] = 0;
  /* -0.102514772783053 -0.439677384690882 0.330913950385712*/
  wx[1] =0.0;//-0.102415;//-1;
  wy[1] =0.0;//-0.43968;//-0.3;
  wz[1] =1.0;//0.330914;// 0.1;
    { 
      double Rt[3][3],Omega[3][3],wSq, w, Ro[3][3], OmegaSq[3][3], M[3][3];
      double ti=4.5, sinw, cosw;
      int k1, k2, k3,i=1;
      wSq = 1.0;
      w = sqrt(wSq);
      if (w != 0.0) 
	{
#if 0
	  if (fabs(w*ti) < 1E-8)
	    {
	      sinw = ti*(1-Sqr(w*ti)/6.0);	  
	      cosw = Sqr(ti)*(0.5 - Sqr(w*ti)/24.0);
	    }
	  else 
	    {
	      sinw = sin(w*ti)/w;
	      cosw = (1.0 - cos(w*ti))/wSq;
	    }
#endif
	  sinw = sin(w*ti)/w;
	  cosw = (1.0 - cos(w*ti))/wSq;
	  Omega[0][0] = 0;
	  Omega[0][1] = -wz[i];
	  Omega[0][2] = wy[i];
	  Omega[1][0] = wz[i];
	  Omega[1][1] = 0;
	  Omega[1][2] = -wx[i];
	  Omega[2][0] = -wy[i];
	  Omega[2][1] = wx[i];
	  Omega[2][2] = 0;
	  OmegaSq[0][0] = -Sqr(wy[i]) - Sqr(wz[i]);
	  OmegaSq[0][1] = wx[i]*wy[i];
	  OmegaSq[0][2] = wx[i]*wz[i];
	  OmegaSq[1][0] = wx[i]*wy[i];
	  OmegaSq[1][1] = -Sqr(wx[i]) - Sqr(wz[i]);
	  OmegaSq[1][2] = wy[i]*wz[i];
	  OmegaSq[2][0] = wx[i]*wz[i];
	  OmegaSq[2][1] = wy[i]*wz[i];
	  OmegaSq[2][2] = -Sqr(wx[i]) - Sqr(wy[i]);

	  for (k1 = 0; k1 < 3; k1++)
	    {
	      for (k2 = 0; k2 < 3; k2++)
		{
		  //Omega[k1][k2] = -Omega[k1][k2];
		  M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
#ifdef MD_USE_CBLAS
		  Rtmp2[k1][k2] = Rtmp[k1][k2] = R[i][k1][k2];
#endif
#if 0
		  if (k1==k2)
		    M[k1][k1] += 1.0;
#endif
		}
	    }
#ifdef MD_USE_CBLAS
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		      3, 3, 3, 1.0, &Rtmp[0][0],
		      3, &M[0][0], 3,
		      1.0, &Rtmp2[0][0], 3);
#if 1
	  for (k1 = 0; k1 < 3; k1++)
	    for (k2 = 0; k2 < 3; k2++)
	      Ro[k1][k2] = Rtmp2[k1][k2];
#endif
#else 
	  Ro[0][0] = uxx[i];
	  Ro[0][1] = uxy[i];
	  Ro[0][2] = uxz[i];
	  Ro[1][0] = uyx[i];
	  Ro[1][1] = uyy[i];
	  Ro[1][2] = uyz[i];
	  Ro[2][0] = uzx[i];
	  Ro[2][1] = uzy[i];
	  Ro[2][2] = uzz[i];
	  Rt[0][0] = uxx[i];
	  Rt[0][1] = uxy[i];
	  Rt[0][2] = uxz[i];
	  Rt[1][0] = uyx[i];
	  Rt[1][1] = uyy[i];
	  Rt[1][2] = uyz[i];
	  Rt[2][0] = uzx[i];
	  Rt[2][1] = uzy[i];
	  Rt[2][2] = uzz[i];
	  for (k1 = 0; k1 < 3; k1++)
	    for (k2 = 0; k2 < 3; k2++)
	      {
		for (k3 = 0; k3 < 3; k3++)
		  Ro[k1][k2] += M[k1][k3]*Rt[k3][k2];
		//Ro[k1][k2] += R[i][k1][k3]*M[k3][k2];
	      }
#endif
	  uxx[i] = Ro[0][0];
	  uxy[i] = Ro[0][1];
	  uxz[i] = Ro[0][2];
	  uyx[i] = Ro[1][0];
	  uyy[i] = Ro[1][1];
	  uyz[i] = Ro[1][2];
	  uzx[i] = Ro[2][0];
	  uzy[i] = Ro[2][1];
	  uzz[i] = Ro[2][2];
	}
      else
	{
	  Omega[0][0] = 0;
	  Omega[0][1] = 0;
	  Omega[0][2] = 0;
	  Omega[1][0] = 0;
	  Omega[1][1] = 0;
	  Omega[1][2] = 0;
	  Omega[2][0] = 0;
	  Omega[2][1] = 0;
	  Omega[2][2] = 0;
	  for (k1 = 0; k1 < 3; k1++)
	    for (k2 = 0; k2 < 3; k2++)
	      {
		Ro[k1][k2] = R[i][k1][k2];
	      }
	}
    }
}

void adjust_norm(double **R);
#ifdef MD_HSVISCO
void calcT(void)
{
  double mass;
  int i;
  OprogStatus.Txy = 0.0;
  OprogStatus.Tyz = 0.0;
  OprogStatus.Tzx = 0.0;
  OprogStatus.Txx = 0.0;
  OprogStatus.Tyy = 0.0;
  OprogStatus.Tzz = 0.0;

  for (i=0; i < Oparams.parnum; i++)
    {
      if (i < Oparams.parnumA)
	mass = Oparams.m[0];
      else 
	mass = Oparams.m[1];
      OprogStatus.Txy += mass*vx[i]*vy[i];
      OprogStatus.Tyz += mass*vy[i]*vz[i];
      OprogStatus.Tzx += mass*vz[i]*vx[i];
      OprogStatus.Txx += mass*vx[i]*vx[i];
      OprogStatus.Tyy += mass*vy[i]*vy[i];
      OprogStatus.Tzz += mass*vz[i]*vz[i];
    }
} 
#endif
/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  int i;
  setToZero(SAVE_LIST, 
	    NULL);  /* Set to zero all the coordinates */

#if 0 
  rx[0] = 0.0;
  ry[0] = +0.6;
  rz[0] = 0.0;

  rx[1] = 0.0;
  ry[1] = -0.6;
  rz[1] = 0.0;
#endif

  FCC(Oparams.parnum, Oparams.m); 
  /* Put the baricenter of each molecule on a FCC lattice, and set 
     their orientations */  

  /* set both atoms velocity to the center of mass velocity */
  comvel(Oparams.parnum, Oparams.T, Oparams.m, 0); 
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
  printf("All'inizio T=%f\n", 2.0 * K / (6.0 * Oparams.parnum - 3.0));

#endif
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
#if 1
  for (i=0; i < Oparams.parnum; i++)
    {
      wx[i]=wy[i]=wz[i]=0.0;
    }
#else
  angvel(); 
#endif
  buildTetrahedras();
  //wrap_initCoord();
}
#ifdef MD_DYNAMIC_OPROG
int dyn_alloc_oprog(void);
void set_dyn_ascii(void);
#endif
/* =========================== >>> usrInitBef <<< ========================== */
void usrInitBef(void)
{
  int i, k;
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
#ifdef MD_DYNAMIC_OPROG
  OprogStatus.ptr = NULL;
  OprogStatus.len = 0;
#endif
#if MD_AB41
  Oparams.bheightAA = Oparams.bheightAB = 0.0;
#else
  Oparams.bheight = 0.0;
#endif
  OprogStatus.tmsd2endA = -1.0;
  OprogStatus.tmsd2endB = -1.0;
     
#ifdef MD_DYNAMIC_OPROG
  OprogStatus.dyn_alloc_oprog = dyn_alloc_oprog;
  OprogStatus.set_dyn_ascii = set_dyn_ascii;
#endif
  OprogStatus.maxbonds = 20;
  Oparams.T = 2.0;
  Oparams.P = 1.0;
  Oparams.M = 5; /* cells in each direction for linked lists */

  Oparams.time = 0.0;
  OprogStatus.tolT = 0.0;
  OprogStatus.targetPhi = 0.0;
  OprogStatus.scalfact = 0.8;
  OprogStatus.reducefact = 0.9;
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  Oparams.Dt = 0.01;
#ifdef MD_DOUBLE_DT
  Oparams.DtR = Oparams.Dt*(3.0/10.0);
#endif
  OprogStatus.avngS     = 0;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
  OprogStatus.scalevel = 0;
  OprogStatus.endtime = 0;
  //OprogStatus.tryharder = 0;
  OprogStatus.rescaleTime = 1.0;
  OprogStatus.brownian = 0;
  OprogStatus.checkGrazing = 0;
#ifdef MD_BIG_DT
  OprogStatus.refTime = 0.0;
  OprogStatus.bigDt = -1.0;
#endif
#ifndef MD_ASYM_ITENS
  for (i = 0; i < 2; i++)
    Oparams.I[i] = 1.0;
#endif

#if defined(MD_ROTDIFF_MIS) && !defined(MD_DYNAMIC_OPROG)
  for (i = 0; i < MAXPAR; i++)
    {
      OprogStatus.lastcolltime[i] = 0.0;
      OprogStatus.sumox[i] = 0.0;
      OprogStatus.sumoy[i] = 0.0;
      OprogStatus.sumoz[i] = 0.0;
    }
#endif
#ifndef MD_DYNAMIC_OPROG
  for (i = 0; i < MAXPAR; i++)
    {
#if 0
      OprogStatus.vcmx0[i] = 0.0;
      OprogStatus.vcmy0[i] = 0.0;
      OprogStatus.vcmz0[i] = 0.0;
#endif
      for (k=0; k < 3; k++)
	OprogStatus.DR[i][k] = 0.0;
    }
#endif
  OprogStatus.eventMult = 100;
  OprogStatus.overlaptol = 0.0001;
  /* Il promo step inizia con un tapping a temperatura T */
  Oparams.m[0] = Oparams.m[1] = 1.0;
  //Oparams.sigma[0][0] = Oparams.sigma[1][1] = Oparams.sigma[1][0]= Oparams.sigma[0][1]=1.0;
  OprogStatus.collCount = 0;
  OprogStatus.crossCount = 0;
  OprogStatus.h=1E-7;
  OprogStatus.assumeOneBond=0;
  OprogStatus.epsd = 0.0005;
  OprogStatus.epsdFast = 0.002;
  OprogStatus.epsdFastR = 0.0025;
  OprogStatus.epsdMax = 0.001;
  OprogStatus.zbrakn = 100;
  OprogStatus.zbrentTol = 0.00001;
  OprogStatus.forceguess = 1;
  OprogStatus.phitol = 1E-12;
  OprogStatus.axestol = 1E-8;
#ifdef MD_CALENDAR_HYBRID
  OprogStatus.scaleHQ = 50;
  OprogStatus.nlistsHQ = 50000;  
  OprogStatus.adjustHQ = 0;
  OprogStatus.baseIndex = 0;
  OprogStatus.curIndex = 0;
  OprogStatus.overthrHQ = 10;
#endif
  OprogStatus.nextSumTime = 0.0;
  OprogStatus.nextcheckTime = 0.0;
  OprogStatus.intervalSum = 1.0;
  OprogStatus.n1 = 160;
  OprogStatus.n2 = 60;
  OprogStatus.storerate = 0.01;
  OprogStatus.KK = 0;
  OprogStatus.JJ = 0;
  OprogStatus.rNebrShell = 2.7; /* the radius of neighbour list shell */
  for (i = 0; i < PE_POINTS; i++)
    OprogStatus.PE[i] = 0;
  /* ======================================================================= */
#ifdef MD_SAVE_REALLY_ALL
  OprogStatus.saveReallyAll=0;
  OprogStatus.readBinTree=0;
  strcpy(OprogStatus.iniTree,"SaveTree");
  strcpy(OprogStatus.iniBak,"BinaryBak");
#endif
  maxcoll=-1;
#ifdef MD_SURV_PROB
  OprogStatus.spdeltat = 10.0;
#endif
}
extern void check (int *overlap, double *K, double *V);
extern double *atomTime, *treeTime, *treeRxC, *treeRyC, *treeRzC;
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
extern int *inCell[2][3], *cellList[4], cellsx[4], cellsy[4], cellsz[4];
#else
extern int *inCell[3], *cellList, cellsx, cellsy, cellsz;
#endif
extern int **tree, cellRange[2*NDIM];
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
extern void PredictColl(int na, int nb, int nl);
extern void PredictCellCross(int na, int nc);
#else
extern void PredictEvent(int, int);
#endif
extern void InitEventList(void);
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
void StartRun(void)
{
  int j, k, n, nl, nc, iA, nl_ignore;

  for (nl=0; nl < 4; nl++)
    {
      for (j = 0; j < cellsx[nl]*cellsy[nl]*cellsz[nl] + Oparams.parnum; j++)
	cellList[nl][j] = -1;
    }
  for (nc=0; nc < 2; nc++)
    {
      /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
      for (n = 0; n < Oparams.parnum; n++)
	{
	  iA = (n < Oparams.parnumA)?0:1;
	  /*
	   * 0 è la lista della specie A per l'interazione A-A
	   * 1 è la lista della specie B per l'interazione B-B
	   * 2 è la lista costituita da molecole B per l'interazione A-B
	   * 3 è la lista costituita da molecole A per l'interazione B-A
	   * */
	  if (iA == 0 && nc == 0)
	    nl = 0;
	  else if (iA == 1 && nc == 0)
	    nl = 1;
	  else if (iA == 0 && nc == 1)
	    nl = 2;
	  else
	    nl = 3;
	  atomTime[n] = Oparams.time;
	  inCell[nc][0][n] =  (rx[n] + L2) * cellsx[nl] / L;
	  inCell[nc][1][n] =  (ry[n] + L2) * cellsy[nl] / L;
#ifdef MD_GRAVITY
	  inCell[nc][2][n] =  (rz[n] + Lz2) * cellsz[nl] / (Lz+OprogStatus.extraLz);
#else
	  inCell[nc][2][n] =  (rz[n] + L2)  * cellsz[nl] / L;
#endif
	  /*printf("nl=%d nc=%d n=%d inCell: %d %d %d cells: %d %d %d\n",
	    nl, nc, n, inCell[nc][0][n], inCell[nc][1][n], inCell[nc][2][n],
	    cellsx[nl], cellsy[nl], cellsz[nl]);
	   */
#if 0
	  if (inCell[0][n]>=cellsx ||inCell[1][n]>= cellsy||inCell[2][n]>= cellsz) 
	    {
	      printf("BOH?!?L:%f L2:%f n:%d rx[n]:%f\n", L, L2, n, rx[n]);
	      printf("(%d,%d,%d) (%d,%d,%d)\n",cellsx , cellsy,cellsz,
		     inCell[0][n],inCell[1][n], inCell[2][n]);
	    }
#endif	  
	  j = (inCell[nc][2][n]*cellsy[nl] + inCell[nc][1][n])*cellsx[nl] + 
	    inCell[nc][0][n] + Oparams.parnum;
	  cellList[nl][n] = cellList[nl][j];
	  cellList[nl][j] = n;
	}
    }
  InitEventList();
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (n = 0; n < Oparams.parnum; n++)
    {
      iA = (n<Oparams.parnumA)?0:1;
      nl_ignore = (n<Oparams.parnumA)?1:0;
      for (nc = 0; nc < 2; nc++)
	{
	  PredictCellCross(n, nc);
	}
      for (nl = 0; nl < 4; nl++)
	{
	  /* iA+2 è la lista Si-O (interazione fra specie diverse)
	   * ma con tutte e sole le molecole iA */
	  if (nl==nl_ignore || nl==iA+2)
	    continue;
	  /*TO BE REMOVED*/
	  //if (nl==3)
	  //continue;
	  //printf("======>qui nl=%d\n", nl);
	  PredictColl(n, -2, nl); 
	}
    }
  //exit(-1);
#if 0
    {
      double dist, rC[3], rD[3], shift[3];
      int i, j;
      for (i = 0; i < Oparams.parnum; i++)
	{
	  j=-1;
	  dist = get_min_dist(i, &j, rC, rD, shift);
	  printf("dist %d:%.8G\n", i, dist);
	}
    }
#endif
  //printf("L'albero e le liste ora dovrebbero essere popolate\n");
  //exit(-1);
}
#else
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
  //exit(-1);
#if 0
    {
      double dist, rC[3], rD[3], shift[3];
      int i, j;
      for (i = 0; i < Oparams.parnum; i++)
	{
	  j=-1;
	  dist = get_min_dist(i, &j, rC, rD, shift);
	  printf("dist %d:%.8G\n", i, dist);
	}
    }
#endif
}
#endif
extern void add_bond(int na, int n, int a, int b);
extern void remove_bond(int na, int n, int a, int b);
extern double calcpotene(void);
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
#define SIGN(X) ((X>0)?1.0:(-1.0)) 
typedef struct {
  double x,y,z;
} XYZ;
typedef struct {
  double point[3];
  //	double grad[3];
  struct {
    int i;	  
    int j;
  } neigh[4];
} MESHXYZ;

extern MESHXYZ **ellips_mesh[2];
void EvalSuperEllipse(double theta,double phi, double a, double b, double c, MESHXYZ *pm)
{
  double cth,cphi,sth,sphi;

  cth = cos(theta);
  cphi = cos(phi);
  sth= sin(theta);
  sphi = sin(phi);
  pm->point[0] = a * cphi * sth;
  pm->point[1] = b * sphi * sth;
  pm->point[2] = c * cth;
#if 1
    {
      FILE* f;
      f = fopen("mesh.dat", "a");
      fprintf(f,"%f %f %f\n", pm->point[0], pm->point[1], pm->point[2]);
      fclose(f);
    }
#endif
#if 0
  pm->grad[0] = 2.0*pm->point[0]/Sqr(a);
  pm->grad[1] = 2.0*pm->point[1]/Sqr(b);
  pm->grad[2] = 2.0*pm->point[2]/Sqr(c);
#endif
}
void add_neighbours(MESHXYZ** mesh, int i, int j)
{
  int n1, n2;
  n1 = OprogStatus.n1;
  n2 = OprogStatus.n2;

  if (i==1)
    {
      mesh[i][j].neigh[0].i = i;
      if (j==0)
	mesh[i][j].neigh[0].j = OprogStatus.n2-1;
      else
	mesh[i][j].neigh[0].j = j-1;
      mesh[i][j].neigh[1].i = i;
      if (j==OprogStatus.n2-1)
	mesh[i][j].neigh[1].j = 0;
      else
	mesh[i][j].neigh[1].j = j+1;
      mesh[i][j].neigh[2].i = i+1;
      mesh[i][j].neigh[2].j = j;
      mesh[i][j].neigh[3].i = -1;
      mesh[i][j].neigh[3].j = -1;
    }
  else if (i==n1/2-1)
    {
      mesh[i][j].neigh[0].i = i;
      if (j==0)
	mesh[i][j].neigh[0].j = OprogStatus.n2-1;
      else
	mesh[i][j].neigh[0].j = j-1;
      mesh[i][j].neigh[1].i = i;
      if (j == OprogStatus.n2-1)
	mesh[i][j].neigh[1].j = 0;
      else	
	mesh[i][j].neigh[1].j = j+1;
      mesh[i][j].neigh[2].i = i-1;
      mesh[i][j].neigh[2].j = j;
      mesh[i][j].neigh[3].i = -1;
      mesh[i][j].neigh[3].j = -1;
    }
  else
    {
      mesh[i][j].neigh[0].i = i;
      if (j==0)
	mesh[i][j].neigh[0].j = OprogStatus.n2-1;
      else
	mesh[i][j].neigh[0].j = j-1;
      mesh[i][j].neigh[1].i = i;
      if (j == OprogStatus.n2-1)
	mesh[i][j].neigh[1].j = 0;
      else	
	mesh[i][j].neigh[1].j = j+1;
      mesh[i][j].neigh[2].i = i-1;
      mesh[i][j].neigh[2].j = j;
      mesh[i][j].neigh[3].i = i+1;
      mesh[i][j].neigh[3].j = j;
    }
}
void build_mesh(MESHXYZ** mesh, double a, double b, double c)
{
  int i,j, n1, n2;
  double theta, phi;
  const double TWOPI=2.0*pi;
  /* n1 = stacks
   * n2 = slides */
  n1 = OprogStatus.n1;
  n2 = OprogStatus.n2;
  for (j=0;j<n2;j++)
    {
      phi = j * TWOPI / (double)n2;
      for (i=1;i<n1/2;i++) 
	{
	  theta = i * TWOPI / (double)n1;
	  EvalSuperEllipse(theta,phi,a,b,c,&mesh[i][j]);
	  add_neighbours(mesh, i, j); 
	}
    }
}
#ifdef MD_GROWTH_CODE
extern double calc_phi(void);
#endif
#ifdef MD_SILICA
extern void assign_bond_mapping(int i, int j);
#endif
extern void writeAllCor(FILE* fs);
extern void writeAsciiPars(FILE* fs, struct pascii strutt[]);

void save_init_conf(void)
{
#ifndef MD_STOREMGL
  const char sepStr[] = "@@@\n";
#endif
  FILE *bf;
  char fileop[512],fileop2[512];
#ifndef MD_STOREMGL
  char fileop3[512];
#endif
  sprintf(fileop2 ,"Store-Init");
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
#ifndef MD_STOREMGL
  writeAsciiPars(bf, opro_ascii);
  fprintf(bf, sepStr);
  writeAsciiPars(bf, opar_ascii);
  fprintf(bf, sepStr);
  printf("qui\n");
#endif
  fprintf(bf, ".Vol: %f\n", L*L*L);
  writeAllCor(bf);
  fclose(bf);
#ifndef MD_STOREMGL
#ifdef MPI
#ifdef MD_MAC
  sprintf(fileop3, "/usr/bin/gzip -f %s_R%d", fileop, my_rank);
#else
  sprintf(fileop3, "/bin/gzip -f %s_R%d", fileop, my_rank);
#endif
#else 
#ifdef MD_MAC
  sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
#else
   sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
#endif
#ifndef MD_NO_SYSTEM
  system(fileop3);
#endif
#endif
}
extern int bound(int na, int n, int a, int b);
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
extern int set_pbonds(int i, int j);
void check_all_bonds(void)
{
  int nl, nn, warn, amin, bmin, i, j, nb, iA, nc, nl_ignore, npbonds;
  double shift[3], dist, dists[MD_PBONDS];
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k;
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }

  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
  warn = 0;
  for ( i = 0; i < Oparams.parnum; i++)
    {
      if (warn)
	break;
      nb = 0;
      iA = (i<Oparams.parnumA)?0:1;
      nl_ignore = (i<Oparams.parnumA)?1:0;

      for (nl = 0; nl < 4; nl++)
	{
	  /* i legami possono essere solo tra Si e O!! */
	  if (nl < 2 )
	    continue;

	  // iA = (i < Oparams.parnumA)?0:1;
	  if (nl < 2)
	    nc = 0;
	  else
	    nc = 1; 
	  if (nl==nl_ignore || nl==iA+2)
	    continue;

	  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	    {
	      jZ = inCell[nc][2][i] + iZ;    
	      shift[2] = 0.;
	      /* apply periodico boundary condition along z if gravitational
	       * fiels is not present */
	      if (jZ == -1) 
		{
		  jZ = cellsz[nl] - 1;    
		  shift[2] = - L;
		} 
	      else if (jZ == cellsz[nl]) 
		{
		  jZ = 0;    
		  shift[2] = L;
		}
	      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
		{
		  jY = inCell[nc][1][i] + iY;    
		  shift[1] = 0.0;
		  if (jY == -1) 
		    {
		      jY = cellsy[nl] - 1;    
		      shift[1] = -L;
		    } 
		  else if (jY == cellsy[nl]) 
		    {
		      jY = 0;    
		      shift[1] = L;
		    }
		  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
		    {
		      jX = inCell[nc][0][i] + iX;    
		      shift[0] = 0.0;
		      if (jX == -1) 
			{
			  jX = cellsx[nl] - 1;    
			  shift[0] = - L;
			} 
		      else if (jX == cellsx[nl]) 
			{
			  jX = 0;   
			  shift[0] = L;
			}
		      j = (jZ *cellsy[nl] + jY) * cellsx[nl] + jX + Oparams.parnum;
		      for (j = cellList[nl][j]; j > -1; j = cellList[nl][j]) 
			{
			  if (i == j)
			    continue;
#if 0 
			  drx = rx[i] - rx[j];
			  shift2[0] = L*rint(drx/L);
			  dry = ry[i] - ry[j];
			  shift2[1] = L*rint(dry/L);
			  drz = rz[i] - rz[j]; 
			  shift2[2] = L*rint(drz/L);
#endif

			  assign_bond_mapping(i, j);
			  dist = calcDistNeg(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
			  npbonds = set_pbonds(i, j);
			  for (nn=0; nn < npbonds; nn++)
			    {
			      if (dists[nn]<0.0 && fabs(dists[nn])>OprogStatus.epsd 
				  && !bound(i,j,mapbondsa[nn], mapbondsb[nn]))
				// && fabs(dists[nn]-Oparams.sigmaSticky)>1E-4)
				{
				  warn=1;
				  printf("dists[%d]:%.15G i=%d j=%d\n", nn, dists[nn], i, j);
#if 0
				  aa = mapbondsa[nn];
				  bb = mapbondsb[nn];
				  wdist=dists[nn];
				  wnn = nn;
				  wj = j;
#endif
				  //nb++;
				}
			      else if (dists[nn]>0.0 && 
				       fabs(dists[nn])> OprogStatus.epsd && 
				       bound(i,j,mapbondsa[nn], mapbondsb[nn]))
				{
				  warn = 2;
				  printf("wrong number of bonds between %d and %d\n", i, j);
				  if (OprogStatus.checkGrazing==1)
				    {
				      remove_bond(i, j, mapbondsa[nn], mapbondsb[nn]);
				    }
				}
			    }
			}
		    }
		}
	    }
	}
      if (warn)
	{
	  mdPrintf(ALL, "[WARNING] wrong number of bonds\n", NULL);
	  sprintf(TXT,"[WARNING] Number of bonds for molecules %d incorrect\n", i);
	  mdPrintf(ALL, TXT, NULL);
	  if (warn==1)
	    mdPrintf(ALL,"Distance < 0 but not bonded, probably a grazing collision occurred\n",NULL);
	  else
	    mdPrintf(ALL,"Distance > 0 but bonded, probably a collision has been missed\n", NULL);
	  //printf("time=%.15G current value: %d real value: %d\n", Oparams.time,
	  //	 numbonds[i], nb);
	  //printf("I've adjusted the number of bonds\n");
	  //printf("Probably a grazing collisions occurred, try to reduce epsd...\n");
	  //store_bump(i,j);
	  if (warn==2)
	    {
	      if (OprogStatus.checkGrazing==2)
		exit(-1);
	      else
		mdPrintf(ALL,"I adjusted the number of bonds...energy won't conserve!", NULL);
	    }
	}
    }
}
#else
void check_all_bonds(void)
{
  int nn, warn, amin, bmin, i, j, nb;
  double shift[3], dist, dists[MD_PBONDS];
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k;
  int md_pbonds;
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }

  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  warn = 0;
  for ( i = 0; i < Oparams.parnum; i++)
    {
      if (warn)
	break;
      nb = 0;
      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCell[2][i] + iZ;    
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
	  for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	    {
	      jY = inCell[1][i] + iY;    
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
	      for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
		{
		  jX = inCell[0][i] + iX;    
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
		  j = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
		  for (j = cellList[j]; j > -1; j = cellList[j]) 
		    {
		      if (i == j)
			continue;
#if 0
		      drx = rx[i] - rx[j];
		      shift2[0] = L*rint(drx/L);
		      dry = ry[i] - ry[j];
		      shift2[1] = L*rint(dry/L);
		      drz = rz[i] - rz[j]; 
		      shift2[2] = L*rint(drz/L);
#endif
#ifdef MD_SILICA
		      assign_bond_mapping(i, j);
#endif
		      dist = calcDistNeg(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
#ifdef MD_THREESPOTS
		      if (i < Oparams.parnumA && j < Oparams.parnumA)
			md_pbonds = MD_PBONDS_AA;
		      else if (i >= Oparams.parnumA && j >= Oparams.parnumA)
			md_pbonds = MD_PBONDS_BB;
		      else
			md_pbonds = MD_PBONDS_AB;
#else
		      md_pbonds = MD_PBONDS;
#endif
		      for (nn=0; nn < md_pbonds; nn++)
			{
			  if (dists[nn]<0.0 && fabs(dists[nn])>OprogStatus.epsd 
			      && !bound(i,j,mapbondsa[nn], mapbondsb[nn]))
			    // && fabs(dists[nn]-Oparams.sigmaSticky)>1E-4)
			    {
			      warn=1;
			      printf("i=%d j=%d pos=%f %f %f vel=%f %f %f\n", i,j, rx[i], ry[i], rz[i], vx[i], vy[i],
		    		     vz[i]);
#if 0
			      aa = mapbondsa[nn];
			      bb = mapbondsb[nn];
			      wdist=dists[nn];
			      wnn = nn;
			      wj = j;
#endif
			      //nb++;
			    }
			  else if (dists[nn]>0.0 && 
				   fabs(dists[nn])> OprogStatus.epsd && 
				   bound(i,j,mapbondsa[nn], mapbondsb[nn]))
			    {
			      warn = 2;
			      printf("wrong number of bonds between %d and %d\n", i, j);
			      if (OprogStatus.checkGrazing==1)
				{
				  remove_bond(i, j, mapbondsa[nn], mapbondsb[nn]);
				}
			    }
			}
		    }
		}
	    }
	}
      if (warn)
	{
	  mdPrintf(ALL, "[WARNING] wrong number of bonds\n", NULL);
	  sprintf(TXT,"[WARNING] Number of bonds for molecules %d incorrect\n", i);
	  mdPrintf(ALL, TXT, NULL);
	  sprintf(TXT,"Step N. %d time=%.15G\n", Oparams.curStep, Oparams.time);
	  mdPrintf(ALL, TXT, NULL);
	  if (warn==1)
	    mdPrintf(ALL,"Distance < 0 but not bonded, probably a grazing collision occurred\n",NULL);
	  else
	    mdPrintf(ALL,"Distance > 0 but bonded, probably a collision has been missed\n", NULL);
	  //printf("time=%.15G current value: %d real value: %d\n", Oparams.time,
	  //	 numbonds[i], nb);
	  //printf("I've adjusted the number of bonds\n");
	  //printf("Probably a grazing collisions occurred, try to reduce epsd...\n");
	  //store_bump(i,j);
	  if (warn==2)
	    {
	      if (OprogStatus.checkGrazing==2)
		exit(-1);
	      else
		mdPrintf(ALL,"I adjusted the number of bonds...energy won't conserve!", NULL);
	    }
	}
    }
}
#endif
/* ======================== >>> usrInitAft <<< ==============================*/
extern void calc_energy(char *msg);
extern int set_pbonds(int i, int j);
#ifdef MD_SAVE_REALLY_ALL
extern void readTreeBondsLL(char *fn);
extern void readBinBak(char *fn);
#endif
#ifdef MD_CALENDAR_HYBRID
extern int *linearLists;
extern void rebuild_linked_list(void);
extern void rebuildCalendar(void);
void estimate_HQ_params(double phi)
{
  /* linear approximation for linear dependence on N */
  double scalevsNfact[4]={0.0877,0.04862,0.9408,7.532}; 
  double nlistsvsNfact[4]={48.67,97.34,244.7,564.182};
  double volfact[4] = {0.01,0.12,0.4,0.7}, msc, mnl, qsc, qnl, scf, nlf;
  int MAXINT = 500000000;
  long long int nlsize;
  int k, k1=-1, k2=-1;
  /* From Gerald Paul J. Comp. Phys. 221, 615 (2006) */

  if (phi < volfact[0])
    {
      /* do not choose values below those for phi=0.01 */
      OprogStatus.scaleHQ = Oparams.parnum*scalevsNfact[0];
      OprogStatus.nlistsHQ = Oparams.parnum*nlistsvsNfact[0];
      return;
    }

  if (phi > volfact[3])
    {
      k1 = 2;
      k2 = 3;
    }
  else
    {
      for (k = 0; k < 3; k++)
	if (phi > volfact[k] && phi < volfact[k+1])
	  {
	    k1 = k;
	    k2 = k+1;
	    break;
	  }
    }
  /* pendenze */
  msc = (scalevsNfact[k2] - scalevsNfact[k1])/(volfact[k2]-volfact[k1]);
  mnl  = (nlistsvsNfact[k2] - nlistsvsNfact[k1])/(volfact[k2]-volfact[k1]);
  /* ordinata all'origine */ 
  qsc = scalevsNfact[k1]-msc*volfact[k1];
  qnl = nlistsvsNfact[k1]-mnl*volfact[k1];
  scf = msc*phi+qsc; 
  nlf = mnl*phi+qnl;
  if (nlf <= 0.0 || scf <= 0.0)
    {
#if 0
      printf("phi=%.15G k1=%d k2=%d\n", phi, k1, k2);
      printf("msc=%.15G qsc=%.15G\n", msc, qsc);
#endif
      printf("[WARNING] estimate of HQ params using default values\n");
      printf("perfomance may be far from optimal, please check\n");
      OprogStatus.scaleHQ = 50;
      OprogStatus.nlistsHQ = 50000;
      printf("scaleHQ=%G nlistsHQ=%d\n", OprogStatus.scaleHQ, OprogStatus.nlistsHQ);
      //exit(-1);
    }
  else
    {
      OprogStatus.scaleHQ = scf*Oparams.parnum;
  
      nlsize = (long long int) (nlf*Oparams.parnum);
      /* evita allocazioni eccessive (oltre i 2Gb) */
      if (nlsize >= (long long int) MAXINT)
	{
	  printf("[WARNING] nlistsHQ will be limited to %d\n", MAXINT);
	  printf("check performance monitoring number of overflows\n");
	  OprogStatus.nlistsHQ = MAXINT;
	}
      else
	OprogStatus.nlistsHQ = (int) nlsize;
    }
#if 0
  printf("nlf=%G BOH nlsize=%lld nlistsHQ=%d\n", nlf, nlsize, OprogStatus.nlistsHQ);
  exit(-1);
#endif	
}
#if 1
void rebuild_linked_list();
void adjust_HQ_params(void)
{
  int targetNE = 15, del=5;
  int k, i, NAVG=10;
  double GOLD = 1.3;
  static int calls=0, sumNumevPQ=0, sumOverevHQ=0;
  double overevHQavg, numevPQavg;

  calls++;
  sumNumevPQ += numevPQ;
  sumOverevHQ += overevHQ;

  if (calls < NAVG)
    return;
  numevPQavg = ((double)sumNumevPQ) / calls;
  overevHQavg= ((double)sumOverevHQ) / calls; 

  /* reset accumulators */
  calls = 0;
  sumNumevPQ = 0;
  sumOverevHQ = 0;

  printf("average values over %d calls: numevPQ=%G overevHQ=%G\n", NAVG, numevPQavg, overevHQavg);
  if (targetNE - del <= numevPQavg && targetNE + del >= numevPQavg && 
      (overevHQavg <= OprogStatus.overthrHQ) )// && numovHQ < totevHQ/OprogStatus.nlistsHQ)
    {
      printf("Hybrid Calendar parameters adjusted!\n");
      OprogStatus.adjustHQ = 0;
    } 
  else
    {
      if (numevPQavg > targetNE)
	{
	  OprogStatus.scaleHQ *= GOLD;
	  //OprogStatus.nlistsHQ *=GOLD;
	 }
      else
	{
	  OprogStatus.scaleHQ /= GOLD;
	  //OprogStatus.nlistsHQ /=GOLD;
	}
      if (overevHQavg > OprogStatus.overthrHQ) 
	OprogStatus.nlistsHQ *= GOLD;
    }
  printf("Adjusting HQ params: scaleHQ=%G nlistsHQ=%d\n", OprogStatus.scaleHQ, OprogStatus.nlistsHQ);
  free(linearLists);
  linearLists = malloc(sizeof(int)*(OprogStatus.nlistsHQ+1));
  UpdateSystem();
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (i=0; i < Oparams.parnum; i++)
    crossevtodel[i] = -1;

  rebuild_linked_list();

  rebuildCalendar();
  if (OprogStatus.intervalSum > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_DOUBLE_DT
  if (OprogStatus.brownian)
    ScheduleEvent(-1, ATOM_LIMIT+12,OprogStatus.nextDtR);
#endif

}
#endif
#endif
#ifdef MD_DYNAMIC_OPROG
int dyn_alloc_oprog(void)
{
  int np, i;  
  void *last_ptr;
  if (OprogStatus.ptr)
    return OprogStatus.len;
  np = Oparams.parnum;
#ifdef MD_ROTDIFF_MIS
  OprogStatus.len = sizeof(double)*10*Oparams.parnum;
#else
  OprogStatus.len = sizeof(double)*6*Oparams.parnum;
#endif
  OprogStatus.ptr = malloc(OprogStatus.len);
  last_ptr = OprogStatus.ptr;
#ifdef MD_ROTDIFF_MIS
  OprogStatus.sumox = (double*)last_ptr;
  OprogStatus.sumoy = OprogStatus.sumox + np;
  OprogStatus.sumoz = OprogStatus.sumoy + np;
  OprogStatus.lastcolltime = OprogStatus.sumoz + np;
  last_ptr = (void*) (OprogStatus.lastcolltime + np);
#endif
  OprogStatus.rxCMi = ((double*)last_ptr);
  OprogStatus.ryCMi = OprogStatus.rxCMi + np;
  OprogStatus.rzCMi = OprogStatus.ryCMi + np;
  OprogStatus.DR = malloc(sizeof(double*)*np);
  for (i=0; i < np; i++)
    {
      OprogStatus.DR[i] = OprogStatus.rzCMi + np + i*3;
    }
#if 0
  OprogStatus.vcmx0 = OprogStatus.DR[np-1] + 3;
  OprogStatus.vcmy0 = OprogStatus.vcmx0 + np;
  OprogStatus.vcmz0 = OprogStatus.vcmy0 + np;
#endif
  OprogStatus.set_dyn_ascii();
  return OprogStatus.len;
}
void set_dyn_ascii(void)
{
  int k;
  k=0;
  do
    {
      if (!strcmp(opro_ascii[k].parName,"rxCMi"))
	opro_ascii[k].ptr = OprogStatus.rxCMi;
      if (!strcmp(opro_ascii[k].parName,"ryCMi"))
	opro_ascii[k].ptr = OprogStatus.ryCMi;
      if (!strcmp(opro_ascii[k].parName,"rzCMi"))
	opro_ascii[k].ptr = OprogStatus.rzCMi;
      if (!strcmp(opro_ascii[k].parName,"DR"))
	{
	  opro_ascii[k].ptr = OprogStatus.DR[0];
	}
#ifdef MD_ROTDIFF_MIS
      if (!strcmp(opro_ascii[k].parName,"sumox"))
	opro_ascii[k].ptr = OprogStatus.sumox;
      if (!strcmp(opro_ascii[k].parName,"sumoy"))
	opro_ascii[k].ptr = OprogStatus.sumoy;
      if (!strcmp(opro_ascii[k].parName,"sumoz"))
	opro_ascii[k].ptr = OprogStatus.sumoz;
      if (!strcmp(opro_ascii[k].parName,"lastcolltime"))
	opro_ascii[k].ptr = OprogStatus.lastcolltime;
#endif
      k++;
    }
  while (strcmp(opro_ascii[k].parName,""));

}
#endif
#ifdef MD_SURV_PROB
int *sp_has_collided, sp_equilib, *sp_coll_type;
double *sp_firstcolltime, sp_start_time;
int sp_tot_collisions;
void sp_reset_fct(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      sp_has_collided[i] = 0;
      sp_coll_type[i] = -1;
    }
  sp_tot_collisions = 0;
}
void save_sp(void)
{
  FILE *fAA, *fBB, *fAB;
  int i;
  fAA = fopen("surv_prob_AA.dat","a");
  fAB = fopen("surv_prob_AB.dat","a");
  fBB = fopen("surv_prob_BB.dat","a");

  for (i=0; i < Oparams.parnum; i++)
    {
      if (sp_coll_type[i] == 0)
 	fprintf(fAA, "%.15G\n", sp_firstcolltime[i]);
      else if (sp_coll_type[i] == 1)
	fprintf(fBB, "%.15G\n", sp_firstcolltime[i]);
      else if (sp_coll_type[i] == 2)
	fprintf(fAB, "%.15G\n", sp_firstcolltime[i]);
    }

  fclose(fAA);
  fclose(fBB);
  fclose(fAB);
}
#endif

void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */
  double dist, phiIni;
  int nn, aa, bb, Nm, k, i, sct, overlap, amin, bmin;
  COORD_TYPE vcmx, vcmy, vcmz;
  COORD_TYPE *m;
  double drx, dry, drz, shift[3], dists[MD_PBONDS];
  int j;
#ifdef MD_SAVE_REALLY_ALL
  int readBinTree=0;
#endif
#ifdef MD_STORE_BONDS
  FILE *fnb;
  char fname[128];
#endif
#ifdef MD_SILICA
  int nl, nc, npbonds;
#endif
  int a;
  /*COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;*/

 /* initialize global varibales */
  pi = 2.0 * acos(0);
#ifdef MD_SAVE_REALLY_ALL
  if (OprogStatus.readBinTree)
    {
      readBinTree = OprogStatus.readBinTree;
      readBinBak(OprogStatus.iniBak);
      newSim=0;
      //printf("OprogStatus.nextSumTime:%.15G", OprogStatus.nextSumTime);
    }
#endif
#if defined(MAXPAR) && !defined(MD_DYNAMIC_OPROG)
  if (Oparams.parnum >= MAXPAR)
    {
      printf("ERROR: Too many particles, increase MAXPAR in sticky.h and recompile\n");
      exit(-1);
    } 
#endif
#ifdef MD_DYNAMIC_OPROG
  OprogStatus.dyn_alloc_oprog();
#endif
#if defined(MD_ROTDIFF_MIS) && defined(MD_DYNAMIC_OPROG)
  for (i = 0; i < Oparams.parnum; i++)
    {
      OprogStatus.lastcolltime[i] = 0.0;
      OprogStatus.sumox[i] = 0.0;
      OprogStatus.sumoy[i] = 0.0;
      OprogStatus.sumoz[i] = 0.0;
    }
#endif
  Nm = Oparams.parnumA;
  parnumA = Oparams.parnumA;
  parnumB = Oparams.parnum - Oparams.parnumA;
  sct = sizeof(COORD_TYPE);
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt <= 0.0)
    OprogStatus.bigDt = 0.0;
#endif
#ifdef MD_GROWTH_CODE
  if (OprogStatus.targetPhi > 0.0)
    {
      printf("[GROWTH SIMULATION] WARNING: during growth spots will be disabled and sphere are additive\n");
      printf("sigma_ij = (sigma_i + sigma_j)*0.5 for every possible pair i,j\n");
    } 
#endif
  invL = 1.0/L;
  L2 = 0.5*L;
#ifdef MD_GRAVITY
  Lz2 = Lz*0.5;
#endif
  poolSize = OprogStatus.eventMult*Oparams.parnum;
  m = Oparams.m;
  Mtot = Oparams.m[0]*parnumA+Oparams.m[1]*parnumB;
#ifndef MD_SILICA
  Oparams.sigma[1][1] = Oparams.sigma[0][0];
#endif
  Oparams.sigma[1][0] = Oparams.sigma[0][1];
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  /* NOTA: nella lista nl=2 con rcut=rcutSiO ci sono tutte e sole le molecole O
   * mentre nella lista nl=3 con rcut=rcutSiO ci sono tutte e sole le molecole Si */
  Oparams.rcut[3] = Oparams.rcut[2];
#endif
  invmA = 1.0/Oparams.m[0];
  invmB = 1.0/Oparams.m[1];
  /* Calcoliamo rcut assumendo che si abbian tante celle quante sono 
   * le particelle */
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  for (nl = 0; nl < 4; nl++)
    {
      if (Oparams.rcut[nl] <= 0.0)
	Oparams.rcut[nl] = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
      cellsx[nl] = L / Oparams.rcut[nl];
      cellsy[nl] = L / Oparams.rcut[nl];
#ifdef MD_GRAVITY
      cellsz[nl] = (Lz+OprogStatus.extraLz) / Oparams.rcut[nl];
#else
      cellsz[nl] = L / Oparams.rcut[nl];
#endif
      printf("[%d] L=%.15G Oparams.rcut: %f %f %f cellsx:%d cellsy: %d cellsz:%d\n", nl, L,
	     Oparams.rcut[0], Oparams.rcut[1], Oparams.rcut[2],
	     cellsx[nl], cellsy[nl], cellsz[nl]);
    }
#else
  if (Oparams.rcut <= 0.0)
    Oparams.rcut = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
  cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
  cellsz = L / Oparams.rcut;
  printf("Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", Oparams.rcut,
	 cellsx, cellsy, cellsz);
#endif
#if defined(MD_GRAVITY)
  g2 = 0.5*Oparams.ggrav;
  mgA = Oparams.m[0]*Oparams.ggrav; 
  mgB = Oparams.m[1]*Oparams.ggrav;
#endif
#ifdef MD_SURV_PROB
  sp_firstcolltime = malloc(sizeof(double)*Oparams.parnum);
  sp_has_collided = malloc(sizeof(int)*Oparams.parnum);
  sp_coll_type = malloc(sizeof(int)*Oparams.parnum);
  sp_reset_fct();
#if 0
  if (Oparams.ntypes==2 && typeNP[0]==1)
    {
      /* if TRAPPING problem */
      for (i=0; i < Oparams.parnum; i++)
	if (typeOfPart[i]==1)
	  {
	    /* particles of type 1 are immobile */
	    vx[i]=vy[i]=vz[i]=0;
	  }
    }
#endif
  sp_equilib=1;
  sp_start_time = Oparams.time;
#endif
  lastcol= malloc(sizeof(double)*Oparams.parnum);
  atomTime = malloc(sizeof(double)*Oparams.parnum);
  lastbump = malloc(sizeof(struct LastBumpS)*Oparams.parnum);
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
  for (nl = 0; nl < 4; nl++)
    cellList[nl] = malloc(sizeof(int)*
			  (cellsx[nl]*cellsy[nl]*cellsz[nl]+Oparams.parnum));
  crossevtodel = malloc(sizeof(int)*Oparams.parnum);
  for (nc = 0; nc < 2; nc++)
    {
      inCell[nc][0] = malloc(sizeof(int)*Oparams.parnum);
      inCell[nc][1]= malloc(sizeof(int)*Oparams.parnum);
      inCell[nc][2] = malloc(sizeof(int)*Oparams.parnum);
    }
#else
  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  inCell[0] = malloc(sizeof(int)*Oparams.parnum);
  inCell[1]= malloc(sizeof(int)*Oparams.parnum);
  inCell[2] = malloc(sizeof(int)*Oparams.parnum);
#endif
#ifdef MD_CALENDAR_HYBRID
  tree = AllocMatI(16, poolSize);
#else
  tree = AllocMatI(12, poolSize);
#endif
  bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
  bonds0 = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
  numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
  numbonds0 = (int *) malloc(Oparams.parnum*sizeof(int));
  bondscache = (int *) malloc(sizeof(int)*OprogStatus.maxbonds);
  treeTime = malloc(sizeof(double)*poolSize);
#if 0
  treeRxC  = malloc(sizeof(double)*poolSize);
  treeRyC  = malloc(sizeof(double)*poolSize);
  treeRzC  = malloc(sizeof(double)*poolSize);
#endif
#ifdef MD_ASYM_ITENS
  Ia = matrix(3, 3);
  Ib = matrix(3, 3);
  invIa = matrix(3, 3);
  invIb = matrix(3, 3);
#endif
  RA = matrix(3, 3);
  RB = matrix(3, 3);
  Rt = matrix(3, 3);
  RtA = matrix(3, 3);
  RtB = matrix(3, 3);
  Aip = matrix(3,3);
  R = malloc(sizeof(double**)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    {
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
      crossevtodel[i] = -1;
#endif
      R[i] = matrix(3, 3);
      lastbump[i].mol = -1;
      lastbump[i].at = -1;
      lastcol[i] = 0.0;
    }
  u2R();
  if (Oparams.curStep == 1)
    {
      check (&overlap, &K, &V);

      if ( overlap ) 
	{
	  printf("ERROR: Particle overlap in initial configuration\n");
	  exit(1);      
	}
    }
  if (newSim)
    {
#ifdef MD_BIG_DT
      OprogStatus.refTime = 0.0;
#endif
      if (OprogStatus.CMreset==-1)
	{
	  comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
	  angvel(); 
	  calc_energy(NULL);
	  scalevels(Oparams.T, K);
	  resetCM();
	}
      else if (OprogStatus.CMreset==-2)
	{
	  comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
	  angvel(); 
	  calc_energy(NULL);
	  scalevels(Oparams.T, K);
	}
      else if (OprogStatus.CMreset==-3)
	{
	  /* 10/05/2010: assegna le velocità da una gaussiana a temperatura Oparams.T 
	     azzerando la velocità del centro di massa e riscalando poi le velocità
	     per avere una temperatura esattamente uguale a Oparams.T */
	  comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
#ifndef MD_SPOT_OFF
	  angvel();
#endif 
	  calc_energy(NULL);
	  scalevels(Oparams.T, K);
	}

      Oparams.time=0.0;
      for (i=0; i < Oparams.parnum; i++)
	atomTime[i] = 0.0;
      OprogStatus.nextcheckTime += fabs(OprogStatus.rescaleTime);
      OprogStatus.nextSumTime += OprogStatus.intervalSum;
      if (OprogStatus.storerate > 0.0)
	OprogStatus.nextStoreTime = OprogStatus.storerate;
      OprogStatus.nextDt += Oparams.Dt;
#ifdef MD_DOUBLE_DT
      OprogStatus.nextDtR += Oparams.DtR;	
#endif
    }
  else
    {
      for (i=0; i < Oparams.parnum; i++)
	atomTime[i] = Oparams.time;
    }
  radat = matrix(Oparams.parnum,NA);
  deltat = matrix(Oparams.parnum,NA);
  a0I = malloc(sizeof(double)*Oparams.parnum);

  maxax = malloc(sizeof(double)*Oparams.parnum);
#ifdef MD_GROWTH_CODE
  axa = malloc(sizeof(double)*Oparams.parnum);
  /* axb is used for interactions between A and B */
#if 0
  axb = malloc(sizeof(double)*Oparams.parnum);
#endif
  for (i=0; i < Oparams.parnumA; i++)
    {
      axa[i] = Oparams.sigma[0][0]*0.5;
#if 0
      axb[i] = Oparams.sigma[0][1]*0.5;
#endif
    } 
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      axa[i] = Oparams.sigma[1][1]*0.5;
#if 0
      axb[i] = Oparams.sigma[0][1]*0.5;
#endif
    } 

#endif
  scdone = malloc(sizeof(int)*Oparams.parnum);
  for (i=0; i < Oparams.parnumA; i++)
    {
      scdone[i] = 0;
      for (a = 0; a < NA; a++)
	{
	  radat[i][a] = Oparams.sigma[0][0];
	}
    }
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      scdone[i] = 0;
      for (a = 0; a < NA; a++)
	{
	  radat[i][a] = Oparams.sigma[1][1];
	}
    }
  /* evaluation of principal inertia moments*/ 
  for (a = 0; a < 2; a++)
    {
#ifdef MD_ASYM_ITENS
      ItensD[a][0] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.b[a])+Sqr(Oparams.c[a]));
      ItensD[a][1] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.a[a])+Sqr(Oparams.c[a]));
      ItensD[a][2] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.a[a])+Sqr(Oparams.b[a]));
#endif
    };
  if (OprogStatus.scalevel)
    {
      printf("CONSTANT TEMPERATURE T=%.15G MOLS=%d MOLA=%d\n", Oparams.T, Oparams.parnum, Oparams.parnum-Oparams.parnumA);
    }
  else
    {
      printf("CONSTANT ENERGY SIMULATION MOLS=%d MOLA=%d\n", Oparams.parnum, Oparams.parnum-Oparams.parnumA);
    }
#ifdef MD_GROWTH_CODE
  printf("INITIAL PHI=%.15G\n", phiIni=calc_phi());
#endif
#ifdef MD_CALENDAR_HYBRID
  /* se scaleHQ == 0 stima automaticamente i parametri 
     se scaleHQ < 0 disabilita il calendario O(1) */

  if (OprogStatus.scaleHQ == 0)
    {
      /* automagically estimate scaleHQ and nlistsHQ parameter */
      estimate_HQ_params(phiIni);
    }
  printf("Using Bounded Increasing Priority Queue, scale=%G nlists=%d\n",OprogStatus.scaleHQ, OprogStatus.nlistsHQ);
  linearLists = malloc(sizeof(int)*(OprogStatus.nlistsHQ+1));
#endif

#ifdef MD_SILICA
  /* write code for silica here!! */
  /* maxax è il diametro del centroide, notare che nel caso della
   * Silica l'interazione bonded è solo tra Si e O, per cui basta avere un solo maxax per 
   * particella!! */
  for (i = 0; i < Oparams.parnum; i++)
    {
#ifdef MD_AB41
      if (i < Oparams.parnumA)
	maxax[i] =Oparams.sigma[0][0] + max(Oparams.sigmaStickyAA,Oparams.sigmaStickyAB) + OprogStatus.epsd;
      else
	maxax[i] =Oparams.sigma[1][1] + Oparams.sigmaStickyAB + OprogStatus.epsd;
#else
      if (i < Oparams.parnumA)
	maxax[i] =(Oparams.sigma[0][0] + Oparams.sigmaSticky + OprogStatus.epsd);
      else
	maxax[i] =(Oparams.sigma[1][1] + Oparams.sigmaSticky + OprogStatus.epsd);
#endif
    }
#else
  /* maxax è il diametro del centroide */
  for (i = 0; i < Oparams.parnum; i++)
    {
      /* scegliere rigorosamente a seconda del modello!!*/
      if (i < Oparams.parnumA)
	{
	  maxax[i] =(Oparams.sigma[0][0] + Oparams.sigmaSticky + OprogStatus.epsd);
	}
      else
	{
	  maxax[i] = Oparams.sigma[1][1] + Oparams.sigmaSticky + OprogStatus.epsd;
	}
    }
#endif
#ifdef MD_SAVE_REALLY_ALL
  if (readBinTree)
    {
      printf("OprogStatus.nextDt:%.15G\n", OprogStatus.nextDt);
      printf("Oparams.time=%.15G\n", Oparams.time);
      readTreeBondsLL(OprogStatus.iniTree);
      return;
    }
#endif

  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
    }
#ifndef MD_SPOT_OFF
#ifdef MD_SILICA
  for ( i = 0; i < Oparams.parnum-1; i++)
    for ( j = i + 1; j < Oparams.parnum; j++)
      {
	/* l'interazione bonded è solo tra Si e O!! */
#if !defined(MD_THREESPOTS) && !defined(MD_AB41)
	if ( !((i < Oparams.parnumA && j >= Oparams.parnumA)||
	       (i >= Oparams.parnumA && j < Oparams.parnumA)) )
	  continue; 
#endif
#ifdef MD_AB41
	if (i >= Oparams.parnumA && j >= Oparams.parnumA)
	  continue;		  
#endif
	drx = rx[i] - rx[j];
	shift[0] = L*rint(drx/L);
	dry = ry[i] - ry[j];
	shift[1] = L*rint(dry/L);
	drz = rz[i] - rz[j]; 
	shift[2] = L*rint(drz/L);
#ifdef MD_SILICA
	assign_bond_mapping(i, j);
#endif
	dist = calcDistNeg(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
	npbonds = set_pbonds(i, j);
	for (nn=0; nn < npbonds; nn++)
	  {
#if 0
	    if (i==0) 
	      printf("i=0 j=%d dists[%d-%d]:%.15G\n", j, mapbondsa[nn], mapbondsb[nn], dists[nn]+Oparams.sigmaSticky);
#endif
	    if (dists[nn]<0.0)
	      {
		//printf("(%d,%d)-(%d,%d)\n", i, mapbondsa[nn], j, mapbondsb[nn]);
		aa = mapbondsa[nn];
		bb = mapbondsb[nn];
		add_bond(i, j, aa, bb);
		add_bond(j, i, bb, aa);
	      }
	  }
      }
#else
  for ( i = 0; i < Oparams.parnum-1; i++)
    for ( j = i + 1; j < Oparams.parnum; j++)
      {
	drx = rx[i] - rx[j];
	shift[0] = L*rint(drx/L);
	dry = ry[i] - ry[j];
	shift[1] = L*rint(dry/L);
	drz = rz[i] - rz[j]; 
	shift[2] = L*rint(drz/L);
	dist = calcDistNeg(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
	for (nn=0; nn < MD_PBONDS; nn++)
	  {
#if 0
	    if (i==0) 
	      printf("i=0 j=%d dists[%d-%d]:%.15G\n", j, mapbondsa[nn], mapbondsb[nn], dists[nn]+Oparams.sigmaSticky);
#endif
	    if (dists[nn]<0.0)
	      {
		//printf("(%d,%d)-(%d,%d)\n", i, mapbondsa[nn], j, mapbondsb[nn]);
		aa = mapbondsa[nn];
		bb = mapbondsb[nn];
		add_bond(i, j, aa, bb);
		add_bond(j, i, bb, aa);
	      }
	  }
      }
#endif
#endif
#ifndef MD_SPOT_OFF
  printf("Energia potenziale all'inizio: %.15f\n", calcpotene()/((double)Oparams.parnum));
#endif
  //exit(-1);
  StartRun(); 
  if (OprogStatus.intervalSum > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_DOUBLE_DT
  if (OprogStatus.brownian)
    ScheduleEvent(-1, ATOM_LIMIT+12,OprogStatus.nextDtR);
#endif
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT + 11, OprogStatus.bigDt);
#endif
  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  

  if (newSim == 1)
    {
      FILE *f;
      /* truncate file to zero lenght */
      f = fopenMPI(absMisHD("energy.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("msdA.dat"), "w+");
      fclose(f);
#ifdef MD_SURV_PROB
      f = fopenMPI(absMisHD("surv_prob_AA.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("surv_prob_BB.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("surv_prob_AB.dat"), "w+");
      fclose(f);
#endif

      if (Oparams.parnum > Oparams.parnumA)
	{
	  f = fopenMPI(absMisHD("msdB.dat"), "w+");
	  fclose(f);
	}
#ifdef MD_ROTDIFF_MIS
      f = fopenMPI(absMisHD("rotMSDA.dat"), "w+");
      fclose(f);
      if (Oparams.parnum > Oparams.parnumA)
	{
	  f = fopenMPI(absMisHD("rotMSDB.dat"), "w+");
	  fclose(f);
	}
#endif
      f = fopenMPI(absMisHD("temp.dat"), "w+");
      fclose(f);
#ifdef MD_HSVISCO
      f = fopenMPI(absMisHD("Ptens.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("press.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("DQtens.dat"), "w+");
      fclose(f);
      OprogStatus.DQTxy = 0.0;
      OprogStatus.DQTyz = 0.0;
      OprogStatus.DQTzx = 0.0;
      OprogStatus.DQTxx = 0.0;
      OprogStatus.DQTyy = 0.0;
      OprogStatus.DQTzz = 0.0;
      OprogStatus.DQWxy = 0.0;
      OprogStatus.DQWyz = 0.0;
      OprogStatus.DQWzx = 0.0;
      OprogStatus.DQWxx = 0.0;
      OprogStatus.DQWyy = 0.0;
      OprogStatus.DQWzz = 0.0;
      OprogStatus.DQWxxHS = 0.0;
      OprogStatus.DQWyyHS = 0.0;
      OprogStatus.DQWzzHS = 0.0;
      OprogStatus.DQWxxST = 0.0;
      OprogStatus.DQWyyST = 0.0;
      OprogStatus.DQWzzST = 0.0;
      OprogStatus.lastcoll = -1;
      OprogStatus.DQxx = 0.0;
      OprogStatus.DQyy = 0.0;
      OprogStatus.DQzz = 0.0;
      calcT();
#endif
      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;

      for(i = 0; i < Oparams.parnum; i++)
	{
#ifdef MD_STORE_BONDS
	  sprintf(fname, "numbonds-%d", i);
	  fnb = fopen(fname, "w+");
	  fclose(fnb);	
#endif
	  /* store the initial positions of particles */
	  OprogStatus.rxCMi[i] = rx[i];
	  OprogStatus.ryCMi[i] = ry[i];
	  OprogStatus.rzCMi[i] = rz[i];	 
	  /* Center of mass velocities */
	  /*printf("(%f,%f,%f)\n", rx[i], ry[i], rx[i]);*/
	  vcmx = vx[i];
	  vcmy = vy[i];
	  vcmz = vz[i];
	  for (k=0; k < 3; k++)
	    OprogStatus.DR[i][k] = 0.0;
#if 0
	  OprogStatus.vcmx0[i] = vcmx;
	  OprogStatus.vcmy0[i] = vcmy;
	  OprogStatus.vcmz0[i] = vcmz;
#endif
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

  /* save_init_conf()  puo' stare anche qui senza alcun problema (così 
     mi evito un po' d'unsulti da parte di valgrind :-) */
  save_init_conf();
  /* printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);*/
#ifdef MD_GRAZING_TRYHARDER
  printf("Grazing try harder code ENABLED!\n");
#endif

}

extern void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3]);
extern void free_matrix(double **M, int n);
/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i;
#ifndef MD_STOREMGL
  const char tipodat[] = "%.15G %.15G %.15G %.15G %.15G %.15G\n";
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n";
#endif
#ifdef MD_STOREMGL
#ifdef MD_SILICA
  int a;
  double rat[5][3], rO[3], **Rl;
  Rl = matrix(3,3);
  /* Oxygen */
  for (i = 0; i < Oparams.parnumA; i++)
    {
      rO[0] = rx[i];
      rO[1] = ry[i];
      rO[2] = rz[i];

      Rl[0][0] = uxx[i];
      Rl[0][1] = uxy[i];
      Rl[0][2] = uxz[i];
      Rl[1][0] = uyx[i];
      Rl[1][1] = uyy[i];
      Rl[1][2] = uyz[i];
      Rl[2][0] = uzx[i];
      Rl[2][1] = uzy[i];
      Rl[2][2] = uzz[i];

      BuildAtomPos(i, rO, Rl, rat);
      /* write coords */
      for (a = 0; a < 5; a++)
	{
	  if (a == 0)
	    {
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[red]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigma[0][1]/2.0);
	    }
	  else if (a < 3)
	    {
#ifdef MD_AB41
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[grey]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigmaStickyAA/2.0);
#else
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[grey]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigmaSticky/2.0);
#endif
	    }
	}
    }
  /* Silicon */
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      rO[0] = rx[i];
      rO[1] = ry[i];
      rO[2] = rz[i];

      Rl[0][0] = uxx[i];
      Rl[0][1] = uxy[i];
      Rl[0][2] = uxz[i];
      Rl[1][0] = uyx[i];
      Rl[1][1] = uyy[i];
      Rl[1][2] = uyz[i];
      Rl[2][0] = uzx[i];
      Rl[2][1] = uzy[i];
      Rl[2][2] = uzz[i];

      BuildAtomPos(i, rO, Rl, rat);
      /* write coords */
      for (a = 0; a < 5; a++)
	{
	  if (a == 0)
	    {
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[YellowGreen]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigma[0][1]/2.0);
	    }
	  else if (a < 5)
	    {
#ifdef MD_AB41
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[grey]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigmaStickyAB/2.0);
#else
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[grey]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigmaSticky/2.0);
#endif
	    }
	}
    }


#else
  int a;
  double rat[5][3], rO[3], **Rl;
  Rl = matrix(3,3);
  for (i = 0; i < Oparams.parnum; i++)
    {
      rO[0] = rx[i];
      rO[1] = ry[i];
      rO[2] = rz[i];

      Rl[0][0] = uxx[i];
      Rl[0][1] = uxy[i];
      Rl[0][2] = uxz[i];
      Rl[1][0] = uyx[i];
      Rl[1][1] = uyy[i];
      Rl[1][2] = uyz[i];
      Rl[2][0] = uzx[i];
      Rl[2][1] = uzy[i];
      Rl[2][2] = uzz[i];

      BuildAtomPos(i, rO, Rl, rat);
      /* write coords */
      for (a = 0; a < 5; a++)
	{
	  if (a == 0)
	    {
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[red]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigma[0][0]/2.0);
	    }
	  else if (a < 3)
	    {
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[blue]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigmaSticky/2.0);
	    }
	  else if (a < 5)
	    {
	      fprintf(fs, "%.15G %.15G %.15G @ %f C[green]\n", rat[a][0], rat[a][1], rat[a][2],
		      Oparams.sigmaSticky/2.0);
	    }
	}
    }
#endif
#else
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat2,rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	      uyz[i], uzx[i], uzy[i], uzz[i]);
    }
#endif
#ifndef MD_STOREMGL 
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat, vx[i], vy[i], vz[i], wx[i], wy[i], wz[i]);
    }
  fprintf(fs, "%.15G\n", L);
#else
  free_matrix(Rl, 3);
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
