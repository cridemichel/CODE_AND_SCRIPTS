#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/
#ifdef MD_LOADMESH
int mesh[100][150][3];
int ntripl[100];
#endif
void kinet(int Nm, COORD_TYPE** velx, COORD_TYPE** vely, COORD_TYPE** velz,
	   COORD_TYPE VOL1);

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
void forceCoeff(void);

extern int ENDSIM;
extern char msgStrA[MSG_LEN];
/* posizioni e velocità di tutti gli atomi in un dischetto
 * (cioè anche quello massless!) */ 
#if 0
double **rallx, **rally, **rallz, **Fallx, **Fally, **Fallz,
  **rallx_old, **rally_old, **rallz_old, *atcharge, 
#endif
double **rx_old, **ry_old, **rz_old, **vxold, **vyold, **vzold, **sigmag;
#ifdef MD_RESPA
double *rxi[NA], *ryi[NA], *rzi[NA];
double **rx_oldLong, **ry_oldLong, **rz_oldLong;
#endif
#ifdef MD_RAPACONSTR
extern double **cvMat, **cvMatInv, ***cvMatInvS, 
       *cDistSq, *vVec, *vVecLong, *curBondLenSq, *cVec[3]; 
extern double **lambvel;
extern int **cMat; 
extern int *cAtom1, *cAtom2;
extern int *laststep;
#endif
double *Fcoeff[3]; 
/* coefficienti delle forze dovute agli atomi senza massa
 * 0 = r1
 * 1 = r2
 * 2 = r3 */
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
extern double Volo1, Volo2, Volot;
/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
#ifdef MD_RESPA
extern int **nebrTabLong, nebrNowLong, nebrTabLenLong, nebrTabMaxLong;
#endif
/* ================================= */
#ifdef MD_RAPACONSTR
extern int *doshake;
#endif

extern COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

extern COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
extern COORD_TYPE  *Rmx, *Rmy, *Rmz;
/* ========================================================================= */
/*=========================== >>> vectProd <<< =========================== */
void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z)
{
  /* DESCRIPTION:
     r3 = [ r1, r2 ] where [ , ] the vectorial product */
  *r3x = r1y * r2z - r1z * r2y; 
  *r3y = r1z * r2x - r1x * r2z;
  *r3z = r1x * r2y - r1y * r2x;
}
void check_distances(char* str)
{
  double dist;
  int i, a, na, baddist=0;
  for(i=0; i < Oparams.parnum; i++)
    {
      for (a=0; a < NA-1;  a++)
	{
	  na = a+1; 
	  dist = sqrt( Sqr(rx[a][i] - rx[na][i]) + Sqr(ry[a][i] - ry[na][i]) + 
		       Sqr(rz[a][i] - rz[na][i]) );
	  if (fabs(dist - Oparams.d) > 1E-8)
	    {
	      printf("STEP: %d (%s) WRONG DISTANCE dist: %f i=%d between atoms %d and %d\n", 
	      Oparams.curStep, str,
	      dist, i, a, na);
	      baddist = 1;
	    }
	}
    }

  if (!baddist)
    printf("[Step N.%d] Tutte le distanze sono corrette\n", Oparams.curStep);
}
void Rand3(double *ex, double *ey, double *ez);

int checkdist(int ir, int ar)
{
  int i, a;
  double dist;
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	  if (abs(a-ar)<=1 && i==ir)
	    continue;
	  dist = sqrt(Sqr(rx[ar][ir]-rx[a][i])+Sqr(ry[ar][ir]-ry[a][i])+
		      Sqr(rz[ar][ir]-rz[a][i]));   
	  if (dist < Oparams.sigma)
	    return 0;
	}
    }
  return 1;
}
/* ============================= >>> FCC <<< ================================*/
void FCC(int Nm, COORD_TYPE D, COORD_TYPE* m)
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
       COORD_TYPE    D                    disk diameter of laponite molecule 
       COORD_TYPE    rxCm, ryCm, rzCm     Molecular Center of mass 
                                          positions             
       COORD_TYPE    ex, ey, ez           half of vector joining atom a and b 
                                          in a molecule 
       COORD_TYPE    rRoot3               1.0 / sqrt ( 3.0 ) */
  int Nc, n;
  const int maxit=100;
  COORD_TYPE dist, norm, dx, dy, dz;
  COORD_TYPE iRoot3, sRoot3; // = 0.5773503;
  COORD_TYPE  Cell, Cell2, rxCm, ryCm, rzCm, Mtot, fact, fact2;
  int a, i, ix, iy, iz, iref, ii, mid, notdone;
  COORD_TYPE bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  COORD_TYPE ex[4], ey[4], ez[4]; /* orientations of each molecule in the 
				     primitive cell */
  COORD_TYPE e2x[4], e2y[4], e2z[4]; /* orientations of each molecule in the 
				     primitive cell */
  double d;
  /*COORD_TYPE d0x, d0y, d0z, d1x, d1y, d1z, d2z, d2y, d2z;*/

  /*printf("FCC Vol: %f\n", Vol);*/
  L = cbrt(Vol);
  Nc = ceil(  pow( ((COORD_TYPE)Nm)/4.0, 1.0/3.0 )  );
  Mtot = m[0] + m[1];
  /*printf("Nc: %d\n", Nc);*/
  /* Calculate the side of the unit cell */
  Cell  = L / ((COORD_TYPE) Nc); /* unit cell length */
  Cell2 = 0.5 * Cell;              /* half unit cell length */
  d = Oparams.d;
  /* Sublattice A */
  iRoot3 = 1.0 / sqrt(3.0);
  sRoot3 = sqrt(3.0);
  /* NOTA: solo in un caso e and e2 sono ortogonali dunque sono stati 
   * costruiti male, risolvere tale problema!!!
   * 0 dovrebbe essere ok mentre 1 2 3 sono scazzati */
  bx[0] =  0.0;
  by[0] =  0.0;
  bz[0] =  0.0;
#if 0
  ex[0] =  iRoot3;
  ey[0] =  iRoot3;
  ez[0] =  iRoot3;
  /* questa è una direzione ortogonale per fissare l'orientazione
   * dei dischi di laponite */
  e2x[0] = iRoot3 * (Sqr(iRoot3) + 1);
  e2y[0] = iRoot3 * (Sqr(iRoot3) - 1);
  e2z[0] = iRoot3 * (Sqr(iRoot3) - 1);
  norm = sqrt(Sqr(e2x[0]) + Sqr(e2y[0]) + Sqr(e2z[0]));
  e2x[0] /= norm;
  e2y[0] /= norm;
  e2z[0] /= norm;
#endif
  /*  Sublattice B */
  bx[1] =  Cell2;
  by[1] =  Cell2;
  bz[1] =  0.0;
#if 0
  ex[1] =  iRoot3;
  ey[1] = -iRoot3;
  ez[1] = -iRoot3;
  /* e2 è ortogonale a e */
  e2x[1] = iRoot3*( Sqr(iRoot3) - 1); /* <== CHECK!!! */
  e2y[1] = iRoot3*(-Sqr(iRoot3) + 1);
  e2z[1] = iRoot3*(-Sqr(iRoot3) - 1);
  norm = sqrt(Sqr(e2x[1]) + Sqr(e2y[1]) + Sqr(e2z[1]));
  e2x[1] /= norm;
  e2y[1] /= norm;
  e2z[1] /= norm;
#endif
  /* Sublattice C */
  bx[2] =  0.0;
  by[2] =  Cell2;
  bz[2] =  Cell2;
#if 0
  ex[2] = -iRoot3;
  ey[2] =  iRoot3;
  ez[2] = -iRoot3;
  e2x[2] = -iRoot3*(Sqr(iRoot3)+1);
  e2y[2] =  iRoot3*(Sqr(iRoot3)-1);
  e2z[2] = iRoot3*(-Sqr(iRoot3)+1);
  norm = sqrt(Sqr(e2x[2]) + Sqr(e2y[2]) + Sqr(e2z[2]));
  e2x[2] /= norm;
  e2y[2] /= norm;
  e2z[2] /= norm;
#endif

  /* Sublattice D */
  bx[3] =  Cell2;
  by[3] =  0.0;
  bz[3] =  Cell2;
#if 0
  ex[3] = -iRoot3;
  ey[3] = -iRoot3;
  ez[3] =  iRoot3;	       
  e2x[3] = (1 - Sqr(iRoot3)) * iRoot3;
  e2y[3] = (1 - Sqr(iRoot3)) * iRoot3;
  e2z[3] = (1 + Sqr(iRoot3)) * iRoot3;
  norm = sqrt(Sqr(e2x[3]) + Sqr(e2y[3]) + Sqr(e2z[3]));
  e2x[3] /= norm;
  e2y[3] /= norm;
  e2z[3] /= norm;
#endif
  ii = 0;
  
  for(iz = 0; iz < Nc; iz++) /* loops over unit cells (that are simply cubes) */ 
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
		  rxCm = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ryCm = by[iref] + Cell * ((COORD_TYPE) iy);
		  rzCm = bz[iref] + Cell * ((COORD_TYPE) iz);
		  mid = NA/2;
		  /*printf("mid:%d NA:%d\n", mid, NA);*/
		  rx[mid][ii+iref] = rxCm;
		  ry[mid][ii+iref] = ryCm;
		  rz[mid][ii+iref] = rzCm;
		  for (n = mid+1; n < NA; n++)
		    {
		      notdone = maxit;
		      while (notdone)
			{
			  Rand3(&dx, &dy, &dz);
			  /* posiziona i vari monomeri */
			  rx[n][ii + iref] = rx[n-1][ii+iref] + d*dx;
			  ry[n][ii + iref] = ry[n-1][ii+iref] + d*dy;
			  rz[n][ii + iref] = rz[n-1][ii+iref] + d*dz;
			  if (checkdist(ii+iref,n))
			    notdone=0;
			  else
			    {
			      notdone--;
			      if (notdone==0)
				printf("FAILED! (%d,%d)\n", n, ii+iref);
			    }
			}
		    }
		   for (n = mid-1; n>=0; n--)
		    {
		      notdone = maxit;
		      while (notdone)
			{
			  Rand3(&dx, &dy, &dz);
			  /* posiziona i vari monomeri */
			  rx[n][ii + iref] = rx[n+1][ii+iref] + d*dx;
			  ry[n][ii + iref] = ry[n+1][ii+iref] + d*dy;
			  rz[n][ii + iref] = rz[n+1][ii+iref] + d*dz;
			  if (checkdist(ii+iref, n))
			    notdone=0;
			  else
			    notdone--;
			}
		    }
		    
		}
	      
	      ii = ii + 4;
	      
	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */
  
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  /* Initial position values are between -0.5 and 0.5 */
	  rx[a][i] = rx[a][i] - 0.5 * L; 
	  ry[a][i] = ry[a][i] - 0.5 * L;
	  rz[a][i] = rz[a][i] - 0.5 * L;
	}
    }
#if 1
  return;
#endif
  /* DEBUG !!!!!!!!!!!! */
  for(i=0; i < Nm; i++)
    {
      dist = sqrt( Sqr(rx[0][i] - rx[1][i]) + Sqr(ry[0][i] - ry[1][i]) + 
		 Sqr(rz[0][i] - rz[1][i]) );
      printf("dist1: %f\n", dist);
      dist = sqrt( Sqr(rx[1][i] - rx[2][i]) + Sqr(ry[1][i] - ry[2][i]) + 
		 Sqr(rz[1][i] - rz[2][i]) );
      printf("dist2: %f\n", dist);
      //dist = sqrt( Sqr(rx[2][i] - rx[0][i]) + Sqr(ry[2][i] - ry[0][i]) + 
	//	 Sqr(rz[2][i] - rz[0][i]) );
      //printf("dist3: %f\n", dist);
    }
  exit(1);
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

  for(i=0; i < 12; i++)
    {
      sum = sum + ranf();
    }

  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}

void bigauss(double sigma1, double sigma2, double c12, double* chsi1p, double* chsi2p)
{
  /* 
     Random bivariate from the standard normal distribution.
     To be used for Brownian Dynamics.
     */
  double chsi1, chsi2;

  chsi1 = gauss();
  chsi2 = gauss();
  *chsi1p = sigma1 * chsi1;
  *chsi2p = sigma2 * ( c12 * chsi1 + sqrt(1 - Sqr(c12))*chsi2);
}
/* See after for declaration */
void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm);
void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm);

void resetCM(int Nm)
{
  COORD_TYPE sumx, sumy, sumz, RCMx, RCMy, RCMz, Rx, Ry, Rz;
  COORD_TYPE Vx, Vy, Vz, *m;
  int i, a, n;

  m = Oparams.m;
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  for(i=0; i < Nm; i++)
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
#if 1
  //printf("SUM (Nm=%d %.15f,%.15f,%.15f\n)\n", Nm, sumx, sumy, sumz);
  //Px=0.0; Py=0.0; Pz=0.0;
#endif
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i=0; i <  Nm; i++)
    {
      for (n=0 ; n < NA; n++)
	{
	  vx[n][i] = vx[n][i] - sumx;
	  vy[n][i] = vy[n][i] - sumy;
	  vz[n][i] = vz[n][i] - sumz;
	}
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }
  
  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx, &Ry, &Rz);
      RCMx += Rx; /* Here RCM is the center of mass of the box */
      RCMy += Ry;
      RCMz += Rz;
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  rx[a][i] -= RCMx;
	  ry[a][i] -= RCMy;
	  rz[a][i] -= RCMz;
	}
    }
  nebrNow=1;
}

/* ========================== >>> comvel <<< =============================== */
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE m[NA])
{
  	int n;
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
  
  Mtot = m[0] + m[1] + m[2];                       /* total mass of molecule */
  rTemp = sqrt(temp / Mtot);  
  /* variance of the velocities distribution function, we assume k = 1 */ 

  for(i=0; i < Nm; i++)
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
      for (n = 0; n < NA; n++)
	{
	  vx[n][i] = rTemp * gauss(); 
	  vy[n][i] = vx[n][i];
	  vz[n][i] = vx[n][i];
	}
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
  
  for(i=0; i < Nm; i++)
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
  for(i=0; i < Nm; i++)
    {
      for (n = 0; n < NA; n++)
	{
	  vx[n][i] = vx[n][i] - sumx;
	  vy[n][i] = vy[n][i] - sumy;
	  vz[n][i] = vz[n][i] - sumz;
	}
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }

  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  for(i=0; i < Nm; i++)
    {
      CoM(i, &Rx, &Ry, &Rz);
      RCMx += Rx; /* Here RCM is the center of mass of the box */
      RCMy += Ry;
      RCMz += Rz;
    }
  
  RCMx /= (COORD_TYPE) Nm;
  RCMy /= (COORD_TYPE) Nm;
  RCMz /= (COORD_TYPE) Nm;

  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  rx[a][i] -= RCMx;
	  ry[a][i] -= RCMy;
	  rz[a][i] -= RCMz;
	}
    }
  kinet(Oparams.parnum, vx, vy, vz, 0);
  printf("temp:%f energia cinetica:%f temp_meas:%f\n",temp, K, 2.0 * K / (NA*3.0 * Oparams.parnum - 3.0));
;
  /* Now the center of mass of the box is in the origin */
}

void Rand3(double *ex, double *ey, double *ez)
{
  double  xisq, xi1, xi2, xi;

  xisq = 1.0;

  while (xisq >= 1.0)
    {
      xi1  = ranf() * 2.0 - 1.0;
      xi2  = ranf() * 2.0 - 1.0;
      xisq = xi1 * xi1 + xi2 * xi2;
    }

  xi = sqrt (fabs(1.0 - xisq));
  *ex = 2.0 * xi1 * xi;
  *ey = 2.0 * xi2 * xi;
  *ez = 1.0 - 2.0 * xisq;
#if 0
  osq   = ox * ox + oy * oy + oz * oz;
  norm  = sqrt(fabs(osq));
  ox    = ox / norm;
  oy    = oy / norm;
  oz    = oz / norm;
#endif
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
#if 0
  return;
#endif
  L = cbrt(Vol);
  invL = 1.0 / L;
  
  Mtot = m[0] + m[1] + m[2]; /* total mass of molecule */

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
void setToZeroPoly(COORD_TYPE** ptr, ...);

/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  setToZeroPoly(SAVE_LIST, 
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
  Oparams.NN = 12;
  Oparams.MM = 6;
  Oparams.d = 0.97;
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
#ifdef MD_RESPA
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  OprogStatus.nrespa = 3;
  OprogStatus.lambda = 0.1;
#endif
  OprogStatus.noLinkedList = 0; /* Use Linked List */
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  OprogStatus.avngS     = 0;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
  OprogStatus.avngEta   = 0;
  OprogStatus.grow = 0;
  OprogStatus.growth = 0.7; /* growth * sigma = minima distanza accettabile */
#ifdef MD_RAPACONSTR
  //OprogStatus.keepInvMat = 0;
  OprogStatus.rapatol = 1E-6;
#endif
#ifdef MD_FENE
  Oparams.kfe = 30.0;
  Oparams.R0 = 1.1;
#endif
  OprogStatus.snapSteps = 0;
  OprogStatus.snapmode = 0;
  for (i=0; i<NA; ++i)
    {
      Oparams.m[i] = 1.0;
    }

  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  /* ======================================================================= */
}
double mindist(int ar, int ir)
{
  int i, a;
  double L,dist, mindist, drx, dry, drz;
  mindist = 1E99;
  L=cbrt(Vol);
  for (i=0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	  if (abs(a-ar)<=1 && i==ir)
	    continue;
	  drx = rx[ar][ir]-rx[a][i];
	  dry = ry[ar][ir]-ry[a][i];
	  drz = rz[ar][ir]-rz[a][i];
	  drx = drx - L*rint(drx/L);
	  dry = dry - L*rint(dry/L);
	  drz = drz - L*rint(drz/L);
	  dist = sqrt(Sqr(drx)+Sqr(dry)+Sqr(drz));   
	  if (dist < mindist)
	    mindist = dist;
	}
    }
  return mindist;
}
void update_sigmas(void)
{
  int i, a, Nm, c;
  double sig;
  Nm = Oparams.parnum;
  c=0;
  for (i=0; i < Nm; i++)
    {
      for (a=0; a < NA; a++)
	{
	  if (sigmag[a][i] == Oparams.sigma)
	    { 
	      c++;
	      continue;
	    }
	  sig = mindist(a,i);
	  if (sig > Oparams.sigma*OprogStatus.growth)
	    {
	      sigmag[a][i] = Oparams.sigma;
	    }
	  else if (sig > sigmag[a][i])
	    sigmag[a][i] = sig;
	}
    }
    
  printf("[update_sigmas] Monomoeri restanti: %d\n", NA*Nm - c);
}
int check_sigmas(void)
{
  int i, a, Nm;
  double sig;
  Nm = Oparams.parnum;
  for (i=0; i < Nm; i++)
    {
      for (a=0; a < NA; a++)
	{
	  if (sigmag[a][i] != Oparams.sigma)
	    return 0;
	}
    }
  return 1;
}
extern void v2p(void);
extern void BuildNebrListNoLinkedLong(int Nm, double rCut);
extern void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut);
extern void LJForceLong(int Nm, double rcutI, double rcutO);
extern void LJForce(int Nm, double rcut);

/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */
#ifdef MD_LOADMESH
  int i1, i2;
  FILE* f;
  const int numk = 98;
  const int numk2av = 150;
#endif
  double sig;
  int Nm, i, sct, cc;
  COORD_TYPE invMtot, d, dSq;
  COORD_TYPE vcmx, vcmy, vcmz;
  COORD_TYPE* m, invrcut3, rcut;
#ifdef MD_RAPACONSTR
  int NB = NA-1;
#endif
  int a, b, n;
  /* le masse dei tre atomi "di base" sono uguali */
  
  /* initialize global varibales */
  pi = 2.0 * acos(0);
  M = Oparams.M;
  printf("NN: %d\n", Oparams.NN);
  /* calcola il numero di siti per disco usando invsp */
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
#ifdef MD_LOADMESH 
  printf("reading mesh files:%s\n",MD_MESHDIR "/ntripl.dat");
  f = fopen(MD_MESHDIR "/ntripl.dat", "r");
  if (!f)
    exit(-1);
  Volo1  = Volo2 = Vol;
  for (i1 = 0; i1 < numk; i1++)
    {
      ntripl[i1] = 0;
      for (i2 = 0; i2 < numk2av; i2++)
	{
	  mesh[i1][i2][0] = mesh[i1][i2][1] = mesh[i1][i2][2] = 0;
	}
    }

  for (i1 = 0; i1 < numk; i1++)
    {
      if (fscanf(f, "%d ", &ntripl[i1]) < 1)
	{
	  break; 
	}
      printf("[%d]=%d ", i1, ntripl[i1]);
    }
  fclose(f);
  f = fopen(MD_MESHDIR "/kmesh.dat", "r");
  for (i1 = 0; i1 < numk; i1++)
    {
      for (i2 = 0; i2 < numk2av; i2++)
	{
	  if (fscanf(f, "%d %d %d ", 
	    	     &mesh[i1][i2][0], &mesh[i1][i2][1], &mesh[i1][i2][2])<3)
	    break;
	  /*
	     printf("[%d][%d]: %d,%d,%d ", i1, i2, mesh[i1][i2][0], 
	     mesh[i1][i2][1], mesh[i1][i2][2]);
	     */
	}
    }
  fclose(f);
#endif 


#if 0
  Oparams.lambdaD = sqrt(  ( Oparams.T * Oparams.epsilon ) / 
    ( 4.0 * pi * (Oparams.nplus * Sqr(Oparams.zplus) + Oparams.nminus * Sqr(Oparams.zminus)) )  ); 
#endif
  vxold = malloc(sizeof(double*)*NA);
  vyold = malloc(sizeof(double*)*NA);
  vzold= malloc(sizeof(double*)*NA);
  rx_old = malloc(sizeof(double*)*NA);
  ry_old = malloc(sizeof(double*)*NA);
  rz_old= malloc(sizeof(double*)*NA);
#ifdef MD_RESPA
  rx_oldLong = malloc(sizeof(double*)*NA);
  ry_oldLong = malloc(sizeof(double*)*NA);
  rz_oldLong = malloc(sizeof(double*)*NA);
#endif
#ifdef MD_RAPACONSTR
  lambvel = malloc(sizeof(double*)*Oparams.parnum);
  cvMatInvS = malloc(sizeof(double**)*NB);
  cvMatInv  = malloc(sizeof(double*)*NB);
  cvMat     = malloc(sizeof(double*)*NB);
  cDistSq = malloc(sizeof(double)*NB);
  curBondLenSq = malloc(sizeof(double)*NB);
  vVec = malloc(sizeof(double)*NB);
  vVecLong = malloc(sizeof(double)*NB);
  cAtom1 = malloc(sizeof(int)*NB);
  cAtom2 = malloc(sizeof(int)*NB);
#endif
  sigmag = malloc(sizeof(double*)*NA);
#if 0
  Rx = malloc(sizeof(double)*Oparams.parnum);
  Ry = malloc(sizeof(double)*Oparams.parnum);
  Rz = malloc(sizeof(double)*Oparams.parnum);
#endif
#if 0
  Fcoeff[0] = malloc(sizeof(double)*NA);
  Fcoeff[1] = malloc(sizeof(double)*NA);
  Fcoeff[2] = malloc(sizeof(double)*NA);
#endif
#ifdef MD_RAPACONSTR
  cVec[0] = malloc(sizeof(double)*NB);
  cVec[1] = malloc(sizeof(double)*NB);
  cVec[2] = malloc(sizeof(double)*NB);
  cMat    = malloc(sizeof(int*)*NB);
  laststep = malloc(sizeof(int)*Oparams.parnum);
  doshake = malloc(sizeof(int)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    {
      lambvel[i] = malloc(sizeof(double)*NB);
      for (a = 0; a < NB; a++)
	lambvel[i][a] = 0;
      doshake[i] = 0;
      laststep[i] = 0;
    }
  for (a = 0; a < NB; a++)
    {
      cvMatInvS[a] = malloc(sizeof(double*)*NB);
      cvMatInv[a] = malloc(sizeof(double)*NB);
      cvMat[a] = malloc(sizeof(double)*NB);
#ifdef MD_RAPACONSTR
      cMat[a] = malloc(sizeof(int)*NA);
#endif
#if 0 
      if (OprogStatus.keepInvMat)
	{
	  for (b = 0; b < NB; b++)
	    cvMatInvS[a][b] = malloc(sizeof(double)*Oparams.parnum);
	}
#endif
    }
#endif
  for (a = 0; a < NA; a++)
    {
      rx_old[a] = malloc(sizeof(double)*Oparams.parnum);
      ry_old[a] = malloc(sizeof(double)*Oparams.parnum);
      rz_old[a] = malloc(sizeof(double)*Oparams.parnum);
#ifdef MD_RESPA
      rx_oldLong[a] = malloc(sizeof(double)*Oparams.parnum);
      ry_oldLong[a] = malloc(sizeof(double)*Oparams.parnum);
      rz_oldLong[a] = malloc(sizeof(double)*Oparams.parnum);
      rxi[a] = malloc(sizeof(double)*Oparams.parnum);
      ryi[a] = malloc(sizeof(double)*Oparams.parnum);
      rzi[a] = malloc(sizeof(double)*Oparams.parnum);
#endif
      vxold[a] = malloc(sizeof(double)*Oparams.parnum);
      vyold[a] = malloc(sizeof(double)*Oparams.parnum);
      vzold[a] = malloc(sizeof(double)*Oparams.parnum);
      sigmag[a] =  malloc(sizeof(double)*Oparams.parnum);
    }
  /* calculate force coefficients */
  /*forceCoeff();*/
 
  for (i=0; i < Nm; i++)
    {
      for (a=0; a < NA; a++)
	{
	  sig = mindist(a,i);
	  if (sig < Oparams.sigma*OprogStatus.growth)
	    sigmag[a][i] = sig;
	  else
	    sigmag[a][i] = Oparams.sigma;
	}
    }
  Rmx    = malloc(sct*Nm);
  Rmy    = malloc(sct*Nm);
  Rmz    = malloc(sct*Nm);
  Dphix = malloc(sct*Nm);
  Dphiy = malloc(sct*Nm);
  Dphiz = malloc(sct*Nm);
  rcut = Oparams.sigma*Oparams.rcut;
  invrcut3 = 1.0/(rcut*rcut*rcut);
  /* correzione alla pressione */
#if defined(SOFT_SPHERE)
  Plrc = 0;
#elif defined(NM_SPHERE)
  Plrc = 0;
#else
  Plrc = 16. * pi * invrcut3 * Sqr(Oparams.parnum * NA) * 
    (2. * Sqr(invrcut3) / 3. - 1.) / 3. / Sqr(Vol);
  printf("Plrc: %f\n", Plrc);
#endif
  ox = malloc(sct*Nm);
  oy = malloc(sct*Nm);
  oz = malloc(sct*Nm);
#if 0
  u1x = malloc(sct*Nm);
  u1y = malloc(sct*Nm);
  u1z = malloc(sct*Nm);
  u2x = malloc(sct*Nm);
  u2y = malloc(sct*Nm);
  u2z = malloc(sct*Nm);
#endif    
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
#ifdef MD_RESPA
  nebrTabMaxLong = OprogStatus.nebrTabFacLong * Nm * NA;
  printf("INIT nebrTabMax: %d\n", nebrTabMaxLong);
  nebrNowLong = 1;
  nebrTabLong = AllocMatI(2, nebrTabMaxLong); 
#endif
  /* =============================================== */
  
  /* Store the Center of Mass initial position for all particles */
  m = Oparams.m;
  
  Mtot = 0.0;
  for(a=0; a < NA; a++)
    {
      m[a] = m[0];
      printf("m:%f\n", Oparams.m[a]);
      Mtot += m[a];
    }
  invMtot = 1.0 / Mtot;

  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      if (OprogStatus.Nose == 0 || OprogStatus.Nose == -1)
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
      if (OprogStatus.Nose == -2)
	{
	  s     = 1.0;
	  s1    = 0.0;
	  s2    = 0.0;
	}

      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;
      for(i=0; i < Oparams.parnum; i++)
	{
	  
	  /* store the initial positions of particles */
	  OprogStatus.rxCMi[i] = OprogStatus.ryCMi[i] = OprogStatus.rzCMi[i] = 0.0;
	  for (a=0; a < NA; a++)
	    {
	      OprogStatus.rxCMi[i] += m[a] * rx[a][i];
	      OprogStatus.ryCMi[i] += m[a] * ry[a][i];
	      OprogStatus.rzCMi[i] += m[a] * rz[a][i]; 
	      for ( cc = 0; cc < 3; cc++)
		OprogStatus.DR[i][cc] = 0.0;
	    }
	  OprogStatus.rxCMi[i] /= Mtot;
	  OprogStatus.ryCMi[i] /= Mtot;
	  OprogStatus.rzCMi[i] /= Mtot;
	  /* zero the accumulators that contains the temporal integral of 
	     the angular velocity for each particle */
	  OprogStatus.sumox[i] = 0.0;
	  OprogStatus.sumoy[i] = 0.0;
	  OprogStatus.sumoz[i] = 0.0;
#if 0 
	  OprogStatus.ox0[i] = ox[i];
	  OprogStatus.oy0[i] = oy[i];
	  OprogStatus.oz0[i] = oz[i];
	  
	  OprogStatus.ux0[i] = ux[i];
	  OprogStatus.uy0[i] = uy[i];
	  OprogStatus.uz0[i] = uz[i];
#endif	  
	  /* Center of mass velocities */
#if 0	
	  vcmx = (m0 * vx[0][i] + m1 * vx[1][i] + m2 * vx[2][i]) * invMtot;
	  vcmy = (m0 * vy[0][i] + m1 * vy[1][i] + m2 * vy[2][i]) * invMtot;
	  vcmz = (m0 * vz[0][i] + m1 * vz[1][i] + m2 * vz[2][i]) * invMtot;
	  OprogStatus.vcmx0[i] = vcmx;
	  OprogStatus.vcmy0[i] = vcmy;
	  OprogStatus.vcmz0[i] = vcmz;
#endif	  
	}
      
      OprogStatus.sumEta   = 0.0;
      OprogStatus.sumTemp  = 0.0;
      OprogStatus.sumPress = 0.0;
      
      for(i=0; i < NUMK; i++) 
	{
	  OprogStatus.sumS[i] = 0.0;
	}
      
      for(i=0; i < MAXBIN; i++)
	{
	  OprogStatus.hist[i] = 0;
	}
    }
#ifdef MD_RESPA
#ifdef MD_RAPACONSTR
  BuildConstraintMatrix();
#endif
#ifdef MD_RESPA_NPT  
  v2p();
#endif
  BuildNebrListNoLinkedLong(Oparams.parnum, Oparams.rcut);
  BuildNebrListNoLinked(Oparams.parnum, OprogStatus.rcutInner);
  LJForceLong(Oparams.parnum, OprogStatus.rcutInner, Oparams.rcut);
  LJForce(Oparams.parnum, OprogStatus.rcutInner);
#endif
  printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);
}
/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i, a;
  for (i = 0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
 	  fprintf(fs, "%.15G %.15G %.15G\n", rx[a][i], ry[a][i], rz[a][i]);
	}
    }
  for (i = 0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", vx[a][i], vy[a][i], vz[a][i]);
	}
    }
#if 0
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", Fx[a][i], Fy[a][i], Fz[a][i]);
	}
    }
#endif
  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G\n", Vol, Vol1, Vol2, Vol1o1, Vol1o2);
  fprintf(fs, "%.15G %.15G %.15G %.15G %.15G",  s, s1, s2, s1o1, s1o2);

}

/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i, a;
  for (i = 0; i < Oparams.parnum; i++)
    {
      for (a = 0;  a < NA; a++)
    	{
	  if (fscanf(fs, "%lf %lf %lf\n", &rx[a][i], &ry[a][i], &rz[a][i]) < 3)
	    {
	      mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	      exit(-1);
	    }
	  //printf("pos i=%d a=%d rx=%.15G ry=%.15G rz=%.15G\n", i, a, rx[a][i], ry[a][i], rz[a][i]);
	}
    }
  for (i = 0; i < Oparams.parnum; i++)
    {
      for (a = 0; a < NA; a++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &vx[a][i], &vy[a][i], &vz[a][i]) < 3)
	    {
	      //printf("i=%d a=%d\n", i, a);;
	      mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	      exit(-1);
	    }
	}
    }
#if 0 
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
#endif

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
