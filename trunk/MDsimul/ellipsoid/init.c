#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/
#define MD_DEBUG31(x) 
/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
void writeAsciiPars(FILE* fs, struct pascii strutt[]);
void writeAllCor(FILE* fs);
void scalevels(double temp, double K);
extern void rebuildNNL(void);
extern char TXT[MSG_LEN];
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
void setToZero(COORD_TYPE* ptr, ...);
double *maxax;
double calc_norm(double *vec);
#ifdef MD_PATCHY_HE
extern struct LastBumpS *lastbump;
#else
extern int *lastbump;
#endif
extern double *lastcol;
double *axa, *axb, *axc;
double **Aip;
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
int **tree, *inCell[3], *cellList, cellsx, cellsy, cellsz, cellRange[2*NDIM];
/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern const double timbig;
double **XbXa, **Xa, **Xb, **RA, **RB, ***R, **Rt, **RtA, **RtB, **powdirs;
#ifdef MD_ASYM_ITENS
double **Ia, **Ib, **invIa, **invIb, **Iatmp, **Ibtmp;
#else
double Ia, Ib, invIa, invIb;
#endif
#ifdef MD_ASYM_ITENS
double *theta0, *phi0, *psi0, *costheta0, *sintheta0, **REt, **REtA, **REtB, *angM, ***RM, **RE0, **Rdot;
double cosEulAng[2][3], sinEulAng[2][3];
#endif
extern double **matrix(int n, int m);
extern int *ivector(int n);
extern double *vector(int n);
int poolSize;
int parnumA, parnumB;
#ifdef MD_PATCHY_HE
int *bondscache, *numbonds, **bonds, *numbonds0, **bonds0;
double *treeRxC, *treeRyC, *treeRzC;
extern int *mapbondsa;
extern int *mapbondsb;
extern int bound(int na, int n, int a, int b);
extern void remove_bond(int na, int n, int a, int b);
extern void assign_bond_mapping(int i, int j);
#endif
double invaSq[2], invbSq[2], invcSq[2];
extern double *fvec, *fvecG, *fvecD;
extern double **fjac,*g,*p,*xold;
extern int *indx;
#ifdef MD_PATCHY_HE
double calcDistNegSP(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, double dists[MD_PBONDS], int bondpair);
extern double spXYZ_A[MD_STSPOTS_A][3];
extern double spXYZ_B[MD_STSPOTS_B][3];
#endif
struct nebrTabStruct *nebrTab;
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

#ifdef MD_PATCHY_HE
extern void check_shift(int i, int j, double *shift);
void check_all_bonds(void)
{
  int nn, warn, amin, bmin, i, j, nb;
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
		      if (! ((i < Oparams.parnumA && j >= Oparams.parnumA)||
			     (i >= Oparams.parnumA && j < Oparams.parnumA)))
			  continue;
		      check_shift(i, j, shift);
		      assign_bond_mapping(i,j);
		      dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
		      for (nn=0; nn < MD_PBONDS; nn++)
			{
			  if (dists[nn]<0.0 && fabs(dists[nn])>OprogStatus.epsd 
			      && !bound(i,j,mapbondsa[nn], mapbondsb[nn]))
			  // && fabs(dists[nn]-Oparams.sigmaSticky)>1E-4)
			    {
			      warn=1;
			      MD_DEBUG31(
			      //printf("dists[1]:%.15G\n", dists[1]);
			      printf("[dist<0]dists[%d]:%.15G\n", nn, dists[nn]);
			      printf("i=%d j=%d %d %d\n", i, j, mapbondsa[nn], mapbondsb[nn]);
			      printf("NA*NA*i+a*NA+b=%d\n", NA*NA*i+mapbondsa[nn]*NA+mapbondsb[nn]);
			      )
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
			      //printf("[dist>0]dists[%d]:%.15G\n", nn, dists[nn]);
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
	  printf("numbonds[%d]:%d bonds[][]:%d\n", i, numbonds[i], bonds[i][0]);
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
      //printf("%d = (%f,%f,%f)\n", i, rx[i], ry[i], rz[i]);
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

void wrap_initCoord(void)
{
  /* A x(-0.603750000000000,4.226250000000000,-0.805000000000000) v(-0.099616130522196,-1.839280599669232,0.357754947051492f)-B x(-2.616250000000000,2.213750000000000,-0.805000000000000) v(1.011838511395152,0.876050550528104,-0.426995365917961)
   * */
  rx[0] = -3.2;
  ry[0] = 0.0;
  rz[0] =  0;
  
  vx[0] = 0.5;
  vy[0] = 0;
  vz[0] = 0;
  /* -0.285316712824933 -0.182347469854598 -0.530547025349427*/

#if 0
  wx[0] = -0.285312;// .003;
  wy[0] = -0.1823475;// -1.5;
  wz[0] = -0.530547;// -0.5;
#else
  wx[0] = -0.3;
  wy[0] = -0.8;
  wz[0] = 0.0;
#endif
#if 0
  uxx[0] = 0.707;
  uyx[0] = -0.707;
  uzx[0] = 0.0;
  uxy[0] = 0.707;
  uyy[0] = 0.707;
  uzy[0] = 0.0;
  uxz[0] = 0.0;
  uyz[0] = 0.0;
  uzz[0] = 1.0;
#endif
  rx[1] = 3.2;
  ry[1] = 0.0;
  
  rz[1] = 0.2;
  vx[1] = -0.5;
  vy[1] = 0.0;
  vz[1] = 0;
  /* -0.102514772783053 -0.439677384690882 0.330913950385712*/
#if 0
  wx[1] -0.102415;//-1;
  wy[1] =-0.43968;//-0.3;
  wz[1] =0.330914;// 0.1;
#else
  wx[1] = 0.3;
  wy[1] = 0.43;
  wz[1] = 0.0;
#endif
  }

void adjust_norm(double **R);
void PredictEventNNL(int na, int nb);
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
  printf("All'inizio T=%f\n", 2.0 * K / (3.0 * Oparams.parnum - 3.0));

#endif
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
  angvel(); 
#if 0
    {
     const char sepStr[] = "@@@\n";
     FILE *bf;
     char fileop[512],fileop2[512], fileop3[512];
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
          sprintf(fileop3, "/bin/gzip -f %s_R%d", fileop, my_rank);
#else 
          sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
	  system(fileop3);
#endif
    }
#endif
  //wrap_initCoord();
}
#ifdef MD_ASYM_ITENS
void calc_omega(int i);
#endif
/* =========================== >>> usrInitBef <<< ========================== */
void usrInitBef(void)
{
  int i;
#ifdef MD_ASYM_ITENS
  int n;
#endif
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
#ifdef MD_PATCHY_HE
    Oparams.sigmaSticky = 1.0;
    Oparams.bheight = 0.0;
    Oparams.bhin = 0.0;
    Oparams.bhout = 0.0;
    OprogStatus.assumeOneBond = 0;
    OprogStatus.checkGrazing = 0;
    OprogStatus.maxbonds = 100;
    OprogStatus.dofA = 6.0;
    OprogStatus.dofB = 5.0;
#else
    OprogStatus.dofA = 5.0;
    OprogStatus.dofB = 5.0;
#endif
#ifdef MD_GRAVITY
    Lz = 9.4;
#endif
    Oparams.T = 2.0;
    Oparams.P = 1.0;
    Oparams.M = 5; /* cells in each direction for linked lists */
    
    Oparams.time = 0.0;
    OprogStatus.tolT = 0.0;
    OprogStatus.targetPhi = 0.0;
    OprogStatus.scalfact = 0.8;
    OprogStatus.reducefact = 0.9;
    OprogStatus.nebrTabFac = 200;
    OprogStatus.rNebrShell = 1.0;
    OprogStatus.useNNL = 0;
    OprogStatus.dist5 = 0;
    OprogStatus.dist8stps = 0;
    OprogStatus.dist5NL = 0;
    OprogStatus.paralNNL = 1;
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
#else
    for (i = 0; i < 2; i++)
      for (n = 0; n < 3; n++)
	Oparams.I[i][n] = 1.0;
#endif
    for (i = 0; i < MAXPAR; i++)
      {
	OprogStatus.lastcolltime[i] = 0.0;
	OprogStatus.sumox[i] = 0.0;
	OprogStatus.sumoy[i] = 0.0;
	OprogStatus.sumoz[i] = 0.0;
      }
    OprogStatus.eventMult = 100;
    OprogStatus.overlaptol = 0.0001;
    /* Il promo step inizia con un tapping a temperatura T */
    Oparams.m[0] = Oparams.m[1] = 1.0;
    //Oparams.sigma[0][0] = Oparams.sigma[1][1] = Oparams.sigma[1][0]= Oparams.sigma[0][1]=1.0;
    OprogStatus.collCount = 0;
    OprogStatus.crossCount = 0;
    OprogStatus.epsd = 0.0005;
    OprogStatus.epsdNL = -1.0;
    OprogStatus.epsdSD = -1.0;
    OprogStatus.epsdGDO = -1.0;
    OprogStatus.h = 1E-10;
#ifdef MD_PATCHY_HE
    OprogStatus.epsdSP = -1.0;
    OprogStatus.epsdFastSP = -1.0;
    OprogStatus.epsdSPNL = -1.0;
    OprogStatus.epsdFastSPNL = -1.0;
#endif
    OprogStatus.epsdFast = 0.002;
    OprogStatus.epsdFastR = 0.0025;
    OprogStatus.epsdMax = 0.001;
    OprogStatus.epsdFastNL = -1.0;
    OprogStatus.epsdFastRNL = -1.0;
    OprogStatus.epsdMaxNL = -1.0;
    /* NOTA: gli epsd NL sono settati a -1.0 poiché
     * se tale valore resta vuol dire che non vengono settati nel file di parametri
     * e dunque assumeranno i valori degli epsd degli ellissoidi altrimenti
     * vengono usati i valori forniti dall'utente (ved. anche usrInitAft() */

    OprogStatus.guessDistOpt = 0;
    OprogStatus.tolSD = 0.01;
    OprogStatus.tolSDlong = -1.0;
    OprogStatus.tolSDconstr= 0.1;
    OprogStatus.tolSDgrad = 0.0;
    OprogStatus.tolAngSD = 0.0;
    OprogStatus.springkSD = 1.0;
    OprogStatus.toldxNR = 0.0;
    OprogStatus.tolAngNR = 0.0;
    OprogStatus.SDmethod = 0;
    OprogStatus.stepSDA = 1.0;
    OprogStatus.stepSDB = 1.0;
    OprogStatus.maxitsSD=200;
    OprogStatus.zbrakn = 100;
    OprogStatus.zbrentTol = 0.00001;
    OprogStatus.forceguess = 1;
    OprogStatus.phitol = 1E-12;
    OprogStatus.axestol = 1E-8;
    OprogStatus.minDist = 4E-8;
    OprogStatus.tmsd2end = -1.0;
    OprogStatus.rmsd2end = -1.0;
    OprogStatus.nextSumTime = 0.0;
    OprogStatus.nextcheckTime = 0.0;
    OprogStatus.intervalSum = 1.0;
    OprogStatus.n1 = 160;
    OprogStatus.n2 = 60;
    OprogStatus.storerate = 0.01;
    OprogStatus.KK = 0;
    OprogStatus.JJ = 0;
    /* Parameters relative to Ruocco AC force
       See: cond-mat/00001311, Ruocco et al. */
    for (i = 0; i < PE_POINTS; i++)
      OprogStatus.PE[i] = 0;
    /* ======================================================================= */
  }
extern void check (int *overlap, double *K, double *V);
double *atomTime, *treeTime, *treeRxC, *treeRyC, *treeRzC;
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
	//printf("inCell: %d, %d, %d\n", inCell[0][n], inCell[1][n], inCell[2][n]);
	//printf("n=%d(%f,%f,%f)\n",n,rx[n], ry[n], rz[n]);
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
    if (OprogStatus.useNNL)
      rebuildNNL();
    for (n = 0; n < Oparams.parnum; n++)
      {
	if (OprogStatus.useNNL)
	  PredictEventNNL(n, -2);
	else
	  PredictEvent(n, -2);
      }
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

#ifdef MD_PATCHY_HE
extern void build_atom_positions(void);
extern void add_bond(int na, int n, int a, int b);
extern double calcpotene(void);
int get_num_pbonds(int i, int j)
{
  return MD_PBONDS;
}
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

MESHXYZ **ellips_mesh[2];
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
#if 0
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
double calc_phi(void);
double costolSDgrad, costolAngSD;
extern double costhrNR;
#ifdef MD_ASYM_ITENS
extern void calc_euler_angles(int i, double **M, double *phi, double *theta, double *psi);
extern double scalProd(double *A, double *B);
double calc_norm(double *vec);
extern void vectProdVec(double *A, double *B, double *C);
void calc_angmom(int i, double **I)
{
  double wv[3], Mvec[3];
  int k1, k2;
  wv[0] = wx[i];
  wv[1] = wy[i];
  wv[2] = wz[i];
  for (k1 = 0; k1 < 3; k1++)
    {
      Mvec[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  Mvec[k1] += I[k1][k2]*wv[k2];
	}
    }
  angM[i] = calc_norm(Mvec);
  Mx[i] = Mvec[0];
  My[i] = Mvec[1];
  Mz[i] = Mvec[2];
}
void calc_RM(int i)
{
  double norm, Mvec[3];
  int k1;
  Mvec[0] = Mx[i];
  Mvec[1] = My[i];
  Mvec[2] = Mz[i];
  /* calcolo il prodotto vettore tra M e l'asse z */
  if (angM[i]==0.0)
    {
      RM[i][0][0] = 1.0;
      RM[i][0][1] = 0.0;
      RM[i][0][2] = 0.0;
      RM[i][1][0] = 0.0;
      RM[i][1][1] = 1.0;
      RM[i][1][2] = 0.0;
      RM[i][2][0] = 0.0;
      RM[i][2][1] = 0.0;
      RM[i][2][2] = 1.0;
      return;
    }
  for (k1 = 0; k1 < 3; k1++)
    RM[i][2][k1] = Mvec[k1]/angM[i];
  RM[i][0][0] = RM[i][2][1];
  RM[i][0][1] = -RM[i][2][0];
  RM[i][0][2] = 0.0;
  norm = calc_norm(RM[i][0]);
  if (norm == 0.0)
    {
      if (RM[i][2][2] >= 0.0)
	{
	  RM[i][0][0] = 1.0;
	  RM[i][0][1] = 0.0;
	  RM[i][0][2] = 0.0;
	  RM[i][1][0] = 0.0;
	  RM[i][1][1] = 1.0;
	  RM[i][1][2] = 0.0;
	  RM[i][2][0] = 0.0;
	  RM[i][2][1] = 0.0;
	  RM[i][2][2] = 1.0;
	}
      else
	{
  	  RM[i][0][0] = 1.0;
	  RM[i][0][1] = 0.0;
	  RM[i][0][2] = 0.0;
	  RM[i][1][0] = 0.0;
	  RM[i][1][1] = -1.0;
	  RM[i][1][2] = 0.0;
	  RM[i][2][0] = 0.0;
	  RM[i][2][1] = 0.0;
	  RM[i][2][2] = -1.0;
	}
      return;
    }
  for (k1 = 0; k1 < 3; k1++)
    RM[i][0][k1] /= norm;
  vectProdVec(RM[i][2],RM[i][0],RM[i][1]); 
#if 0
  //printf("1 scal 2: %.15G\n",  scalProd(RM[i][0], RM[i][2]));
  for (k1=0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      if (isnan(RM[i][k1][k2]))
	{
	  printf("angM routine RM[%d][%d][%d]:%.15G\n", i, k1,k2, RM[i][k1][k2]);
	  printf("angM[%d]:%.15G\n", i, angM[i]);
	  printf("Mvec=%f,%f,%f norm=%.15G\n", Mvec[0], Mvec[1], Mvec[2], norm);
	  exit(-1);
	}
#endif
}
#if 0
void calc_angmom(int i, double **I)
{
  double wv[3], Mvec[3], th, costh, sinth, phi, cosphi, sinphi, VP[3], Mu[3], VPN;
  int k1, k2;

  wv[0] = wx[i];
  wv[1] = wy[i];
  wv[2] = wz[i];
  for (k1 = 0; k1 < 3; k1++)
    {
      Mvec[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  Mvec[k1] += I[k1][k2]*wv[k2];
	}
    }
  angM[i] = calc_norm(Mvec);
  /* calcolo il prodotto vettore tra M e l'asse z */
  for (k1 = 0; k1 < 3; k1++)
    Mu[k1] = Mvec[k1]/angM[i];
  VP[0] = Mu[1];
  VP[1] = -Mu[0];
  VP[2] = 0.0;
  /* e ora calcolo RM */
  VPN = calc_norm(VP);
  for (k1 = 0; k1 < 3; k1++)
    VP[k1] = VP[k1]/VPN;
  th = acos(Mu[2]);//acos(Mvec[2]/angM[i]);
  costh = cos(th);
  sinth = sin(th);
  if (sinth==0.0)
    {
      RM[i][0][0] = 1.0;
      RM[i][0][1] = 0.0;
      RM[i][0][2] = 0.0;
      RM[i][1][0] = 0.0;
      RM[i][1][1] = 1.0;
      RM[i][1][2] = 0.0;
      RM[i][2][0] = 0.0;
      RM[i][2][1] = 0.0;
      RM[i][2][2] = 1.0;
      return;
    }
  phi = acos(VP[0]); 
  if (VP[1] < 0.0)
    phi = 2*pi - phi;
  cosphi = cos(phi);
  sinphi = sin(phi);
  RM[i][0][0] = cosphi;
  RM[i][0][1] = sinphi;
  RM[i][0][2] = 0.0;
  RM[i][1][0] = -costh*sinphi;
  RM[i][1][1] =  costh*cosphi;
  RM[i][1][2] = sinth;
  RM[i][2][0] = sinth*sinphi;
  RM[i][2][1] = -sinth*cosphi;
  RM[i][2][2] = costh;
  for (k1=0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      if (isnan(RM[i][k1][k2]))
	{
	  printf("angM routine RM[%d][%d][%d]:%.15G\n", i, k1,k2, RM[i][k1][k2]);
	  printf("sinth:%.15G cosphi: %.15G sinphi: %.15G costh: %.15G\n", sinth, cosphi,
		 sinphi, costh);
	  printf("VP=(%f,%f,%f) \n", VP[0], VP[1], VP[2]);
	  exit(-1);
	}
}
#endif
void tRDiagR(int i, double **M, double a, double b, double c, double **Ri);
void upd_refsysM(int i)
{
  int k1, k2, k3;
  calc_RM(i); 
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	RE0[k1][k2] = 0.0;
	for (k3 = 0; k3 < 3; k3++)
	  RE0[k1][k2] += R[i][k1][k3]*RM[i][k2][k3];
#if 0
	if (isnan(RM[i][k1][k2]))
	  {printf("RM[%d][%d][%d]:%.15G\n", i, k1, k2, RM[i][k1][k2]);
	  exit(-1);}
#endif
     }
#if 0
  printf("RE0[%d]=\n",i);
  print_matrix(RE0, 3);
#endif
  calc_euler_angles(i, RE0, &phi0[i], &theta0[i], &psi0[i]);
  //printf("RE0[2][2]: %.15G costheta=%.15G\n", RE0[2][2], cos(theta0[i]));
  costheta0[i] = cos(theta0[i]);
  sintheta0[i] = sin(theta0[i]);
}
#endif 
/* ======================== >>> usrInitAft <<< ==============================*/
void RDiagtR(int i, double **M, double a, double b, double c, double **Ri);
extern double max(double a, double b);
extern double max3(double a, double b, double c);
double calc_shell(void)
{
#ifdef MD_PATCHY_HE
  int aa, kk;
  double v[3], norm, delta;
#endif
  double deltamax = 0.0;
#ifdef MD_PATCHY_HE
  for (aa = 0; aa < MD_STSPOTS_A; aa++)
    {
       norm = calc_norm(spXYZ_A[aa]);
       for (kk=0; kk < 3; kk++)
	 v[kk] = (norm + Oparams.sigmaSticky*0.5) * spXYZ_A[aa][kk] /  norm;
       delta = max3(v[0]-Oparams.a[0],v[1]-Oparams.b[0],v[2]-Oparams.c[0]);
       if (aa == 0 || delta > deltamax)
	 deltamax  = delta;
    }
  for (aa = 0; aa < MD_STSPOTS_B; aa++)
    {
       norm = calc_norm(spXYZ_B[aa]);
       for (kk=0; kk < 3; kk++)
	 v[kk] = (norm + Oparams.sigmaSticky*0.5) * spXYZ_B[aa][kk] /  norm;
       delta = max3(v[0]-Oparams.a[1],v[1]-Oparams.b[1],v[2]-Oparams.c[1]);
       if (aa == 0 || delta > deltamax)
	 deltamax  = delta;
    }
#else
  deltamax = 0.0;
#endif
  //printf("deltamax=%.15G\n", deltamax);
  return deltamax;
}
double calc_nnl_rcut(void)
{
  double rcutA, rcutB;
  rcutA = 2.0*sqrt(Sqr(Oparams.a[0]+OprogStatus.rNebrShell)+
		   Sqr(Oparams.b[0]+OprogStatus.rNebrShell)+
		   Sqr(Oparams.c[0]+OprogStatus.rNebrShell));
  rcutB = 2.0*sqrt(Sqr(Oparams.a[1]+OprogStatus.rNebrShell)+
		   Sqr(Oparams.b[1]+OprogStatus.rNebrShell)+
		   Sqr(Oparams.c[1]+OprogStatus.rNebrShell));
  return 1.01*max(rcutA, rcutB);
}

void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */
  int Nm, i, sct, overlap;
  COORD_TYPE vcmx, vcmy, vcmz, MAXAX;
  COORD_TYPE *m;
#ifdef MD_PATCHY_HE
  double shift[3], drx, dry, drz, dist, dists[MD_PBONDS];
  int j, amin, bmin, nn, aa, bb, NPB;
  double distSPA, distSPB;
#endif
  int a;
  /*COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;*/

  /* initialize global varibales */
  pi = 2.0 * acos(0);

  Nm = Oparams.parnumA;
  parnumA = Oparams.parnumA;
  parnumB = Oparams.parnum - Oparams.parnumA;
  if (OprogStatus.epsdNL == -1.0)
    OprogStatus.epsdNL = OprogStatus.epsd;
  if (OprogStatus.epsdFastNL == -1.0)
    OprogStatus.epsdFastNL = OprogStatus.epsdFast;
  if (OprogStatus.epsdMaxNL == -1.0)
    OprogStatus.epsdMaxNL = OprogStatus.epsdMaxNL;
  if (OprogStatus.epsdFastRNL == -1.0)
    OprogStatus.epsdFastRNL = OprogStatus.epsdFastR;
#ifdef MD_PATCHY_HE
  if (OprogStatus.epsdSP == -1.0)
    OprogStatus.epsdSP = OprogStatus.epsd;
  if (OprogStatus.epsdFastSP == -1.0)
    OprogStatus.epsdFastSP = OprogStatus.epsdFast;
  if (OprogStatus.epsdSPNL == -1.0)
    OprogStatus.epsdSPNL = OprogStatus.epsdNL;
  if (OprogStatus.epsdFastSPNL == -1.0)
    OprogStatus.epsdFastSPNL = OprogStatus.epsdFastNL;
#endif
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt <= 0.0)
    OprogStatus.bigDt = 0.0;
#endif
  sct = sizeof(COORD_TYPE);
  costolSDgrad = cos(OprogStatus.tolSDgrad);
  costolAngSD =  fabs(cos(OprogStatus.tolAngSD) - 1.0);
  costhrNR = cos(OprogStatus.tolAngNR);
  if (OprogStatus.dist5==0)
    OprogStatus.dist8stps = 0;
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
  Mred[0][0] = Mred[1][1] = 0.5;
  Mred[0][1] = Mred[1][0] = (Oparams.m[0]*Oparams.m[1])/(Oparams.m[0]+Oparams.m[1]);
#if 0
  printf("massA: %f massB: %f sigmaA:%f sigmaB:%f sigmaAB:%f\n", Oparams.m[0], Oparams.m[1],
	 Oparams.sigma[0][0], Oparams.sigma[1][1], Oparams.sigma[0][1]);
#endif
#if defined(MD_GRAVITY)
  g2 = 0.5*Oparams.ggrav;
  mgA = Oparams.m[0]*Oparams.ggrav; 
  mgB = Oparams.m[1]*Oparams.ggrav;
#endif

  if (OprogStatus.epsdSD < 0.0)
    OprogStatus.epsdSD = Sqr(OprogStatus.epsd);
  if (OprogStatus.tolSDlong < 0.0)
    OprogStatus.tolSDlong = OprogStatus.tolSD;
  /*    
   ** CHECK FOR PARTICLE OVERLAPS **
   ** CALCULATE ENERGY            ** */
  lastcol= malloc(sizeof(double)*Oparams.parnum);
  atomTime = malloc(sizeof(double)*Oparams.parnum);
#ifdef MD_PATCHY_HE
  lastbump =  malloc(sizeof(struct LastBumpS)*Oparams.parnum);
#else
  lastbump = malloc(sizeof(int)*Oparams.parnum);
#endif
#ifdef MD_PATCHY_HE
  tree = AllocMatI(12, poolSize);
  bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
  bonds0 = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
  numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
  numbonds0 = (int *) malloc(Oparams.parnum*sizeof(int));
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
  XbXa = matrix(3, 3);
  indx=ivector(8); 
  fjac=matrix(8, 8);
  g=vector(8);
  p=vector(8); 
  xold=vector(8); 
  fvec=vector(8); 
  fvecD=vector(8);
  fvecG=vector(8);
#ifdef MD_ASYM_ITENS
  Oparams.I[0][1] = Oparams.I[0][0];
  Oparams.I[1][1] = Oparams.I[1][0];
  Ia = matrix(3, 3);
  Ib = matrix(3, 3);
  Iatmp = matrix(3,3);
  Ibtmp = matrix(3,3);
  invIa = matrix(3, 3);
  invIb = matrix(3, 3);
  angM = malloc(sizeof(double)*Oparams.parnum);
  phi0 = malloc(sizeof(double)*Oparams.parnum);
  psi0 = malloc(sizeof(double)*Oparams.parnum);
  costheta0 = malloc(sizeof(double)*Oparams.parnum);
  sintheta0 = malloc(sizeof(double)*Oparams.parnum);
  theta0 =    malloc(sizeof(double)*Oparams.parnum);
  psi0   =    malloc(sizeof(double)*Oparams.parnum);
  phi0   =    malloc(sizeof(double)*Oparams.parnum);
  angM   =    malloc(sizeof(double)*Oparams.parnum);
  REt = matrix(3,3);
  REtA = matrix(3,3);
  REtB = matrix(3,3);
  RE0 = matrix(3,3);
  Rdot = matrix(3,3);
  RM = malloc(sizeof(double**)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++) 
    RM[i] = matrix(3, 3);
#endif
  powdirs = matrix(6,6);
#if 0
  ellips_mesh[0]=malloc(sizeof(MESHXYZ*)*OprogStatus.n1*3);
  ellips_mesh[1]=malloc(sizeof(MESHXYZ*)*OprogStatus.n1*3);
  for (i = 0; i < OprogStatus.n1; i++)
    {
      ellips_mesh[0][i] = malloc(sizeof(MESHXYZ)*OprogStatus.n2*3);
      ellips_mesh[1][i] = malloc(sizeof(MESHXYZ)*OprogStatus.n2*3);
    }
  build_mesh(ellips_mesh[0], Oparams.a[0], Oparams.b[0], Oparams.c[0]);
  build_mesh(ellips_mesh[1], Oparams.a[1], Oparams.b[1], Oparams.c[1]);
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
      R[i] = matrix(3, 3);
#ifdef MD_PATCHY_HE
      lastbump[i].mol = -1;
      lastbump[i].at = -1;
#else
      lastbump[i] = -1;
#endif
      lastcol[i] = 0.0;
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
      FILE *f;
#ifdef MD_BIG_DT
      OprogStatus.refTime = 0.0;
#endif
      Oparams.time=0.0;
      /* truncate file to zero lenght */
#ifdef MD_PATCHY_HE
      f = fopenMPI(absMisHD("energy.dat"), "w+");
      fclose(f);
#endif
      f = fopenMPI(absMisHD("MSDA.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("rotMSDA.dat"), "w+");
      fclose(f);
      if (Oparams.parnum > Oparams.parnumA)
	{
	  f = fopenMPI(absMisHD("MSDB.dat"), "w+");
	  fclose(f);
	  f = fopenMPI(absMisHD("rotMSDB.dat"), "w+");
	  fclose(f);
	}
      f = fopenMPI(absMisHD("temp.dat"), "w+");
      fclose(f);
#ifdef MD_HSVISCO
      f = fopenMPI(absMisHD("Ptens.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("press.dat"), "w+");
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
      for (i=0; i < Oparams.parnum; i++)
	{
	  atomTime[i] = 0.0;
	  OprogStatus.sumox[i] = 0.0;
	  OprogStatus.sumoy[i] = 0.0;
	  OprogStatus.sumoz[i] = 0.0;
	}
      OprogStatus.nextcheckTime += fabs(OprogStatus.rescaleTime);
      OprogStatus.nextSumTime += OprogStatus.intervalSum;
      if (OprogStatus.storerate > 0.0)
	OprogStatus.nextStoreTime = OprogStatus.storerate;
#ifdef MD_STOREMGL
      if (OprogStatus.endtime < 0.0)
	OprogStatus.nextStoreTime = 0.0;
#endif
      OprogStatus.nextDt += Oparams.Dt;
    }
  else
    {
      for (i=0; i < Oparams.parnum; i++)
	atomTime[i] = Oparams.time;
    }

  axa = malloc(sizeof(double)*Oparams.parnum);
  axb = malloc(sizeof(double)*Oparams.parnum);
  axc = malloc(sizeof(double)*Oparams.parnum);
  maxax = malloc(sizeof(double)*Oparams.parnum);
  scdone = malloc(sizeof(int)*Oparams.parnum);
  if (OprogStatus.useNNL)
    {
      printf("I'm going to use NNL, good choice to go fast :)\n");
      nebrTab = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
    }
  for (i=0; i < Oparams.parnumA; i++)
    {
      if (OprogStatus.useNNL)
	{
	  nebrTab[i].len = 0;
	  nebrTab[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	  nebrTab[i].shift = matrix(OprogStatus.nebrTabFac, 3);
	  nebrTab[i].R = matrix(3, 3);
	}
      scdone[i] = 0;
      axa[i] = Oparams.a[0];
      axb[i] = Oparams.b[0];
      axc[i] = Oparams.c[0];
    }
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      if (OprogStatus.useNNL)
	{
	  nebrTab[i].len = 0;
	  nebrTab[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	  nebrTab[i].R = matrix(3, 3);
	}
      scdone[i] = 0;
      axa[i] = Oparams.a[1];
      axb[i] = Oparams.b[1];
      axc[i] = Oparams.c[1];
    }
  printf(">>>> phi=%.12G L=%f (%f,%f,%f)\n", calc_phi(), L, Oparams.a[0], Oparams.b[0], Oparams.c[0]); 
  if (Oparams.parnumA < Oparams.parnum)
    printf("semi-axes of B (%f, %f ,%f)\n",Oparams.a[1], Oparams.b[1], Oparams.c[1]);
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
#ifdef MD_PATCHY_HE
  build_atom_positions();
  distSPA = 0.0;
  for (aa = 0; aa < MD_STSPOTS_A; aa++)
    {
      dist = calc_norm(spXYZ_A[aa])+Oparams.sigmaSticky*0.5;
      //printf("calc_norm[%d]:%.15G\n", aa, calc_norm(spXYZ_A[aa]));
      if (dist > distSPA)
	distSPA = dist;
    }
  distSPB = 0.0;
  for (aa = 0; aa < MD_STSPOTS_B; aa++)
    {
      dist = calc_norm(spXYZ_B[aa])+Oparams.sigmaSticky*0.5;
      if (dist > distSPB)
	distSPB = dist;
    }
  //printf("distSPA: %.15G distSPB: %.15G\n", distSPA, distSPB);
#endif
  MAXAX = 0.0;
  for (i = 0; i < Oparams.parnum; i++)
    {
      maxax[i] = 0.0;
      a=(i<Oparams.parnumA)?0:1;
      if (Oparams.a[a] > maxax[i])
	maxax[i] = Oparams.a[a];
      if (Oparams.b[a] > maxax[i])
	maxax[i] = Oparams.b[a];
      if (Oparams.c[a] > maxax[i])
	maxax[i] = Oparams.c[a];
      //printf("distSPA=%.15G distSPB=%.15G\n", distSPA, distSPB);
#ifdef MD_PATCHY_HE
      //printf("maxax bef[%d]: %.15G\n", i, maxax[i]*2.0);
      if (i < Oparams.parnumA)
	{
	  if (distSPA > maxax[i])
	    maxax[i] = distSPA;
	}
      else
	{
	  if (distSPB > maxax[i])
	    maxax[i] = distSPB;
	}
#endif
      maxax[i] *= 2.0;
      if (maxax[i] > MAXAX)
	MAXAX = maxax[i];
      //printf("maxax aft[%d]: %.15G\n", i, maxax[i]);
    }
  /* Calcoliamo rcut assumendo che si abbian tante celle quante sono 
   * le particelle */
  if (OprogStatus.useNNL)
    {
      /* in questo modo rNebrShell è la "buccia" rispetto al minimo 
       * parallelepipedo che include l'ellissoide più gli sticky spots */
      if (OprogStatus.targetPhi <= 0.0)
	OprogStatus.rNebrShell += calc_shell();
      printf("[INFO] I've adjusted rNebrShell to %.15G\n", OprogStatus.rNebrShell);	  
      if (Oparams.rcut <= 0.0)
	{
	  Oparams.rcut = calc_nnl_rcut();
	  printf("[INFO] I've chosen rcut= %.15G\n", Oparams.rcut);
	}
    }
  else
    {
      if (Oparams.rcut <= 0.0)
	Oparams.rcut = MAXAX*1.01;
    }
  //Oparams.rcut = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
  cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
  cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
  cellsz = L / Oparams.rcut;
#endif 
  printf("Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", Oparams.rcut,
	 cellsx, cellsy, cellsz);
  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  inCell[0] = malloc(sizeof(int)*Oparams.parnum);
  inCell[1]= malloc(sizeof(int)*Oparams.parnum);
  inCell[2] = malloc(sizeof(int)*Oparams.parnum);
#ifdef MD_PATCHY_HE
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
    }
  printf("L=%f parnum: %d parnumA: %d\n", L, Oparams.parnum, Oparams.parnumA);
  printf("sigmaSticky=%.15G\n", Oparams.sigmaSticky);
  for ( i = 0; i < Oparams.parnum-1; i++)
    for ( j = i + 1; j < Oparams.parnum; j++)
      {
	/* l'interazione sticky è solo fra fra A e B! */
	if (!((i < Oparams.parnumA && j >= Oparams.parnumA)|| 
	      (i >= Oparams.parnumA && j < Oparams.parnumA)))
	  continue;
	drx = rx[i] - rx[j];
	shift[0] = L*rint(drx/L);
	dry = ry[i] - ry[j];
	shift[1] = L*rint(dry/L);
	drz = rz[i] - rz[j]; 
	shift[2] = L*rint(drz/L);
	assign_bond_mapping(i, j);
	dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
	NPB = get_num_pbonds(i, j);
	for (nn=0; nn < NPB; nn++)
	  {
   	    if (dists[nn] < 0.0)
	      {
		//printf("(%d,%d)-(%d,%d)\n", i, mapbondsa[nn], j, mapbondsb[nn]);
		aa = mapbondsa[nn];
		bb = mapbondsb[nn];
		add_bond(i, j, aa, bb);
		add_bond(j, i, bb, aa);
	      }
	  }
      }
  printf("Energia potenziale all'inizio: %.15f\n", calcpotene());
#endif
#ifdef MD_ASYM_ITENS
  for (i=0; i < Oparams.parnum; i++)
    {
      calc_omega(i);
      angM[i] = sqrt(Sqr(Mx[i])+Sqr(My[i])+Sqr(Mz[i]));
      upd_refsysM(i);
    }
#endif
  //exit(-1);
  StartRun(); 
  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
#ifdef MD_GRAVITY
  if (OprogStatus.scalevel || OprogStatus.taptau > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
#else
  if (OprogStatus.scalevel)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
#endif
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+11,OprogStatus.bigDt);
#endif
  MD_DEBUG(printf("scheduled rebuild at %.15G\n", nltime));
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
#if MD_BIG_DT
      OprogStatus.refTime = 0.0;
      OprogStatus.bigDt = -1.0;
#endif
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
extern double rA[3], rB[3];
#ifdef MD_PATCHY_HE
void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3]);
#endif
/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i;
#ifdef MD_STOREMGL
  int nn;
  double ratA[NA][3];
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G @ %.15G %.15G %.15G C[%s]\n";
  
#else
  const char tipodat[] = "%.15G %.15G %.15G %.15G %.15G %.15G\n";
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n";
#endif
#ifdef MD_STOREMGL
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (i < Oparams.parnumA)
	{
	  fprintf(fs, tipodat2,rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	  	  uyz[i], uzx[i], uzy[i], uzz[i], Oparams.a[0], Oparams.b[0], Oparams.c[0],
	  	  "red");
	}
      else
	{
	  fprintf(fs, tipodat2, rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	  	  uyz[i], uzx[i], uzy[i], uzz[i],
	  	  Oparams.a[1], Oparams.b[1], Oparams.c[1],
	  	  "green");
	}
#ifdef MD_PATCHY_HE
	  rA[0] = rx[i];
	  rA[1] = ry[i];
	  rA[2] = rz[i];
	  BuildAtomPos(i, rA, R[i], ratA);
	  for (nn = 1; nn < ((i < Oparams.parnumA)?MD_STSPOTS_A+1:MD_STSPOTS_B+1); nn++)
	    fprintf(fs,"%.15f %.15f %.15f @ %.15G C[orange]\n", 
		    ratA[nn][0], ratA[nn][1], ratA[nn][2], Oparams.sigmaSticky*0.5);
#endif
    }
#else
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat2, rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
	      uyz[i], uzx[i], uzy[i], uzz[i]);
    }
#endif
#ifndef MD_STOREMGL 
  for (i = 0; i < Oparams.parnum; i++)
    {
#ifdef MD_ASYM_ITENS
      fprintf(fs, tipodat, vx[i], vy[i], vz[i], Mx[i], My[i], Mz[i]);
#else
      fprintf(fs, tipodat, vx[i], vy[i], vz[i], wx[i], wy[i], wz[i]);
#endif
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
#ifdef MD_ASYM_ITENS
      if (fscanf(fs, "%lf %lf %lf\n", &Mx[i], &My[i], &Mz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
#else
      if (fscanf(fs, "%lf %lf %lf\n", &wx[i], &wy[i], &wz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
#endif
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
