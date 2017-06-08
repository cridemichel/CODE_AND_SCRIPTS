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
void setToZero(COORD_TYPE* ptr, ...);

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
double mg, g2;
extern COORD_TYPE pi, L, invL, L2;   
#ifdef MD_GRAVITY
extern double Lz2;
#endif
extern COORD_TYPE W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx; 

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern const double timbig;
int poolSize;
#ifndef MD_POLYDISP
double *radii;
#else
double *radiiINI;
#endif
int *scdone;
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
#ifdef MD_HSVISCO
void calcT(void)
{
  double mass;
  int i;
  OprogStatus.Txy = 0.0;
  OprogStatus.Tyz = 0.0;
  OprogStatus.Tzx = 0.0;
  mass = Oparams.m;
  for (i=0; i < Oparams.parnum; i++)
    {
      OprogStatus.Txy += mass*vx[i]*vy[i];
      OprogStatus.Tyz += mass*vy[i]*vz[i];
      OprogStatus.Tzx += mass*vz[i]*vx[i];
    }
} 
#endif

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
      rz[i] = rz[i] - 0.5 * Lz + Oparams.sigma*0.5 + 0.1;
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
#ifdef MD_FULL_LANG
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
#endif
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
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE m, int resetCM)
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
  /*COORD_TYPE Px, Py, Pz;*/
  int i;
  
  rTemp = sqrt(temp / Oparams.m);  
  /* variance of the velocities distribution function, we assume k = 1 */ 
  for (i = 0; i < Nm; i++)
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
#if defined(MD_FULL_LANG)
      v2x[i] = rTemp * gauss();
      v2y[i] = rTemp * gauss();
      v2z[i] = rTemp * gauss();
#elif defined(MD_MICRO_LANG)
      //v2x[i] = vx[i];
      //v2y[i] = vy[i];
      //v2z[i] = vz[i];
#endif
      //printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
    }
  if (OprogStatus.brownian)
    return;
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

  if (!resetCM)
    return;
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
#ifndef MD_GRAVITY
      rz[i] -= RCMz;
#endif
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
  comvel(Oparams.parnum, Oparams.T, Oparams.m, 1); 
  
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
  //angvel(Oparams.parnum, Oparams.T, Oparams.m, Oparams.d); 
}
#ifdef MD_LOADMESH
int mesh[100][150][3];
int ntripl[100];
#endif

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
#ifdef MD_GRAVITY
  Lz = 9.4;
#endif
  Oparams.T = 2.0;
  Oparams.P = 1.0;
  Oparams.M = 5; /* cells in each direction for linked lists */
  Oparams.wallDiss = 1.0;
  Oparams.partDiss = 1.0;
  OprogStatus.time = 0.0;
  OprogStatus.tolT = 0.0;
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  OprogStatus.targetPhi = 0.0;
#ifdef MD_POLYDISP
  OprogStatus.polydisp = 0.0;
  OprogStatus.polycutoff = 5.0;
#endif
  OprogStatus.phitol = 1E-12;
  OprogStatus.scalfact = 0.99;
  OprogStatus.axestol = 1E-8;
  OprogStatus.accrcmz = 0.0;
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  OprogStatus.avngS     = 0;
  OprogStatus.tc = 1E-5;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
  OprogStatus.taptau = 0.0;
  OprogStatus.rzup = 0.0;
  OprogStatus.expandFact= 1.0;
  OprogStatus.scalevel = 0;
  OprogStatus.brownian = 0;
  OprogStatus.eventMult = 100;
  OprogStatus.overlaptol = 0.0001;
  /* Il promo step inizia con un tapping a temperatura T */
  OprogStatus.quenchend = 0.0;
  Oparams.m = 1.0;
  Oparams.sigma = 1.0;
  OprogStatus.collCount = 0;
  OprogStatus.crossCount = 0;
  OprogStatus.wallcollCount = 0;
  OprogStatus.nextSumTime = 0.0;
  OprogStatus.nextcheckTime = 0.0;
  OprogStatus.intervalSum = 1.0;
  OprogStatus.storerate = 0.01;
  OprogStatus.KK = 0;
  OprogStatus.JJ = 0;
  OprogStatus.checkquenchTime = 1.0;
  OprogStatus.rescaleTime = 1.0;
  OprogStatus.brownian = 0;
  OprogStatus.numquench = 0;
  OprogStatus.maxquench = 0;
  OprogStatus.extraLz = 10.0;
  OprogStatus.rhobh = 0.0;
  OprogStatus.brownian = 0;
  Oparams.Dt = 0.01;
  #ifdef MD_FPBROWNIAN
  Oparams.xi = 100;
  #endif
  OprogStatus.vztap = 10.0;
  /* Parameters relative to Ruocco AC force
     See: cond-mat/00001311, Ruocco et al. */
  OprogStatus.rNebrShell = 2.7; /* the radius of neighbour list shell */
  for (i = 0; i < PE_POINTS; i++)
    OprogStatus.PE[i] = 0;
#ifdef MD_POLYDISP
  OprogStatus.ext_radii=0;
  strcpy(OprogStatus.radiiFile, "*");
#endif
  /* ======================================================================= */
}
extern void check ( double sigma, int *overlap, double *K, double *V);
extern double *treeTime, *atomTime;
extern int **tree, *inCell[3], *cellList, cellsx, cellsy, cellsz, cellRange[2*NDIM];
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
      inCell[0][n] =  (rx[n] + L2) * cellsx / L;
      inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
      inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
      inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
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
double *rhoz;
extern double calc_phi(void);
/* ======================== >>> usrInitAft <<< ==============================*/
#ifdef MD_POLYDISP
double PHI0;
#endif
#ifdef MD_POLYDISP
void read_radii(char *fn)
{
  FILE *fp;
  int i;
  fp=fopen(fn,"r");
  for (i=0; i < Oparams.parnum; i++)
    {
      fscanf(fp, "%lf ", &radii[i]);
      radii[i]/=2.0;
    }
  fclose(fp);
}
#endif
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
#ifdef MD_POLYDISP
  FILE *fp;
#endif
  int jZmax; 
  int Nm, i, sct, overlap;
  COORD_TYPE vcmx, vcmy, vcmz;
  COORD_TYPE m;
#ifndef MD_RAND_ALL  
  /*COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;*/
  if (OprogStatus.quenchend == -1 && OprogStatus.taptau > 0 && newSim &&
      !strcmp(OprogStatus.inifile,"*"))
    {
      srand(0);
      comvel(Oparams.parnum, 0.1, Oparams.m, 1); 
#if !defined(MPI)
      if (mdseed>=0) 
	srand(mdseed);
      else
	srand(((int)time(NULL)));
#else
      srand(my_rank*((int)time(NULL)));
#endif
    } 
#endif
  /* initialize global varibales */
  pi = 2.0 * acos(0);
  
  Nm = Oparams.parnum;
  sct = sizeof(COORD_TYPE);

  invL = 1.0/L;
  L2 = 0.5*L;
  rcmz = -Lz*0.5;
#ifdef MD_GRAVITY
  Lz2 = Lz*0.5;
#endif
  poolSize = OprogStatus.eventMult*Oparams.parnum;
  m = Oparams.m;
  /* Calcoliam rcut assumendo che si abbian tante celle quante sono 
   * le particelle */
#ifdef MD_POLYDISP
  if (OprogStatus.ext_radii)
    {
      if (strcmp(OprogStatus.radiiFile,"*"))
	read_radii(OprogStatus.radiiFile);
      else
	read_radii("radii.dat");
    }      
#endif
  if (Oparams.rcut <= 0.0)
    {
#ifdef MD_POLYDISP
      double maxrad=-1.0;
      if (!strcmp(OprogStatus.inifile,"*") && !OprogStatus.ext_radii)
	{
	  if (OprogStatus.polydisp > 0.0)
	    {
	      Oparams.rcut = 1.01*Oparams.sigma*(1.0 + OprogStatus.polydisp*OprogStatus.polycutoff);
	    }
	  else
	    {
	      Oparams.rcut = 1.01*Oparams.sigma;
	    }
	}
      else
	/*      else
		Oparams.rcut = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
		*/

    	{
	  for (i = 0; i < Oparams.parnum; i++)
	    {
	      if (radii[i] > maxrad)
		maxrad = radii[i];
	    }
	  Oparams.rcut = 1.01*maxrad*2.0;
	}
#else
      Oparams.rcut = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
#endif
    }
  cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
  cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
  cellsz = L / Oparams.rcut;
#endif 
  printf("Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", Oparams.rcut,
	 cellsx, cellsy, cellsz);
#if defined(MD_GRAVITY)
  g2 = 0.5*Oparams.ggrav;
  mg = Oparams.m*Oparams.ggrav; 
#endif
#ifdef MD_LOADMESH 
  printf("reading mesh files:%s\n",MD_MESHDIR "/ntripl.dat");
  f = fopen(MD_MESHDIR "/ntripl.dat", "r");
  if (!f)
    exit(-1);

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


/*    
   ** CHECK FOR PARTICLE OVERLAPS **
   ** CALCULATE ENERGY            ** */
#if 0
  lastcol= malloc(sizeof(double)*Oparams.parnum);
#endif
  atomTime = malloc(sizeof(double)*Oparams.parnum);
  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  inCell[0] = malloc(sizeof(int)*Oparams.parnum);
  inCell[1] = malloc(sizeof(int)*Oparams.parnum);
  inCell[2] = malloc(sizeof(int)*Oparams.parnum);
#ifdef MD_GRAVITY
  jZmax = (int) ((Lz / (OprogStatus.extraLz + Lz)) * cellsz)+1;  
  jZmax = (jZmax < cellsz)?jZmax:cellsz;
#else
  jZmax = cellsz;  
#endif
  rhoz = malloc(sizeof(double)*(jZmax+1));
  tree = AllocMatI(9, poolSize);
  treeTime = malloc(sizeof(double)*poolSize);
  if (OprogStatus.CMreset==-1)
    {
      comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
      resetCM(1);
    }
  else if (OprogStatus.CMreset==-2)
    {
      comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
    }

#ifndef MD_POLYDISP
  radii = malloc(sizeof(double)*Oparams.parnum);
#else
  if (OprogStatus.targetPhi > 0.0)
    radiiINI = malloc(sizeof(double)*Oparams.parnum);
#endif
  scdone = malloc(sizeof(int)*Oparams.parnum);
 for (i=0; i < Oparams.parnum; i++)
    {
      scdone[i] = 0;
#ifdef MD_POLYDISP
      if (OprogStatus.targetPhi > 0.0)
	radiiINI[i] = radii[i];

      if (newSim)
	{
	 /* if (OprogStatus.polydisp <= 0.0)
	    radii[i] = Oparams.sigma*0.5;
	  else*/
	  if (!strcmp(OprogStatus.inifile,"*") && !OprogStatus.ext_radii)
	    {
	      if (OprogStatus.polydisp > 0.0)
		{
		  do
		    {
		      radii[i] = (OprogStatus.polydisp*gauss() + 1.0)*Oparams.sigma*0.5; 
		    }
		  while (radii[i] < Oparams.sigma*0.5*(1.0 - OprogStatus.polycutoff*OprogStatus.polydisp) || radii[i] > Oparams.sigma*0.5*(1.0 + OprogStatus.polycutoff*OprogStatus.polydisp));
		  //printf("%.15G\n", radii[i]);
		}
	      else
		radii[i] = Oparams.sigma*0.5;
	    }
	}
#else
      radii[i] = Oparams.sigma*0.5;
#endif
    }
  if (Oparams.curStep == 1)
    {
      check ( Oparams.sigma, &overlap, &K, &V);
     
      if ( overlap ) 
     	{
 	  printf("ERROR: Particle overlap in initial configuration\n");
 	  exit(1);      
  	}
    }
  if (newSim)
    {
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
	atomTime[i] = OprogStatus.time;
    }
  StartRun(); 
  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel) 
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_POLYDISP
  printf(">>>> Volume fraction: %.15G\n", PHI0=calc_phi());
#else
  printf(">>>> Volume fraction: %.15G\n", calc_phi());
#endif
  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      FILE *f;
      /* truncate file to zero lenght */
      f = fopenMPI(MD_HD_MIS "T.dat", "w");
      fclose(f);
      f = fopenMPI(absMisHD("msd.dat"), "w+");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "Vz2.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rcmz.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rho.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rhoz.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rcmz_ist.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rho_ist.dat", "w");
      fclose(f);
#ifdef MD_HSVISCO
      OprogStatus.DQTxy = 0.0;
      OprogStatus.DQTyz = 0.0;
      OprogStatus.DQTzx = 0.0;
      OprogStatus.DQWxy = 0.0;
      OprogStatus.DQWyz = 0.0;
      OprogStatus.DQWzx = 0.0;
      OprogStatus.lastcoll = -1;
      calcT();
#endif
      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;
     
      if (!strcmp(OprogStatus.inifile,"*"))
	{
	  for (i= 0; i < Oparams.parnum; i++)
	    {
	      /*printf("i=%d v(%.15G,%.15G,%.15G)\n", i, vx[i], vy[i], vz[i]);*/
	      lastcol[i] = 0.0;
	    }
	}
      for(i = 0; i < Oparams.parnum; i++)
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
  const char tipodat[] = "%.15G %.15G %.15G\n";
#ifdef MD_POLYDISP
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, "%.15G %.15G %.15G %.15G\n", rx[i], ry[i], rz[i], radii[i]);
    }
#else
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat, rx[i], ry[i], rz[i]);
    }
#endif
  for (i = 0; i < Oparams.parnum; i++)
    {
      fprintf(fs, tipodat, vx[i], vy[i], vz[i]);
    }
#if defined(MD_MICRO_LANG) || defined(MD_FULL_LANG)
  for (i = 0; i < Oparams.parnum; i++)
    {
#ifdef MD_MICRO_LANG
      //v2x[i] = vx[i];
      //v2y[i] = vy[i];
      //v2z[i] = vz[i];
#endif
      //fprintf(fs, tipodat, v2x[i], v2y[i], v2z[i]);
    }
#endif
  fprintf(fs, "%.15G %.15G\n", L, Lz);
}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i;

#ifdef MD_POLYDISP
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf %lf\n", &rx[i], &ry[i], &rz[i], &radii[i]) < 4)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
#else
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf\n", &rx[i], &ry[i], &rz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
#endif 
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf\n", &vx[i], &vy[i], &vz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
      
    }
#if 0
#if defined(MD_MICRO_LANG) || defined(MD_FULL_LANG)
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf\n", &v2x[i], &v2y[i], &v2z[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[vel] reading ascii file\n", NULL);
	  exit(-1);
	}
      
    }

#endif
#endif
  if (fscanf(fs, "%lf %lf\n",  &L, &Lz) < 2)
    {
      mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
      exit(-1);
    }
      
}
