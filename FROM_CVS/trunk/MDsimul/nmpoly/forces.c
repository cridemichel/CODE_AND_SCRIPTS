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
#if 0
extern double **rallx, **rally, **rallz, **Fallx, **Fally, **Fallz,
  **rallx_old, **rally_old, **rallz_old, *atcharge; 
#endif
extern double **rx_old, **ry_old, **rz_old, **sigmag;
extern double *Fcoeff[3];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE W, K, WC, T1xx, T1yy, T1zz,
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

extern COORD_TYPE *ox, *oy, *oz; /* Angular velocities of each particle */

extern COORD_TYPE *ux, *uy, *uz; /* Molecular orientations */
extern COORD_TYPE  *Rmx, *Rmy, *Rmz;

/* =========================== >>> LINKED LIST <<< ========================= */
/* STATEMENT FUNCTION TO GIVE CELL INDEX */

int iCell(int ix, int iy, int iz) 
{
  return ((ix + M) % M) + ((iy + M) % M) * M + ( (iz + M) % M )* M * M;
}

/* =============================== >>> maps <<< ============================ */
void maps(void)
{      
  /*  DESCRIPTION:
      routine to set up a list of neighbouring cells               
                                                               
      PRINCIPAL VARIABLES:                                         
                                                               
      int M                 number of cells in each direction 
      int mapSize           size of cell-cell map             
      int map[mapSize]      list of neighbouring cells        
      
      USAGE:                                                       
      
      this subroutine sets up a list of the thirteen neighbouring  
      cells of each of the small cells in the central box. The     
      effects of the periodic boundary conditions are included.    
      The subroutine is called once at the beginning of the        
      simulation and the map is used in the force subroutine */
  
  int ix, iy, iz, iMap; 
  /* FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL  */
  loop(iz, 1, M)
    {
      loop(iy, 1, M)
	{
	  loop(ix, 1, M)
	    {  
	      iMap = iCell(ix, iy, iz) * 13;
	      /* iCell associate to the cell a unique number betwee 1 and 
		 M * M * M */ 
	      
	      /* Give a certain cell, the cells that have at least 
		 an edge common with that one are 26 */ 
	      map[iMap + 0]  = iCell(ix + 1, iy, iz); /*first neighbour cell*/ 
	      map[iMap + 1]  = iCell(ix + 1, iy + 1, iz);
	      map[iMap + 2]  = iCell(ix, iy + 1, iz);
	      map[iMap + 3]  = iCell(ix - 1, iy + 1, iz);
	      map[iMap + 4]  = iCell(ix + 1, iy, iz - 1);
	      map[iMap + 5]  = iCell(ix + 1, iy + 1, iz - 1);
	      map[iMap + 6]  = iCell(ix, iy + 1, iz - 1);
	      map[iMap + 7]  = iCell(ix - 1, iy + 1, iz - 1);
	      map[iMap + 8]  = iCell(ix + 1, iy, iz + 1);
	      map[iMap + 9]  = iCell(ix + 1, iy + 1, iz + 1);
	      map[iMap + 10] = iCell(ix, iy + 1, iz + 1);
	      map[iMap + 11] = iCell(ix - 1, iy + 1, iz + 1);
	      map[iMap + 12] = iCell(ix, iy, iz + 1);
	    }
	}
    }
}
/* ========================= >>> links <<< =================================*/
void  links(int Nm, COORD_TYPE rcut)
{
  /* DESCRIPTION:
     Routine to set up linked list, the head of chain arrays
     PRINCIPAL VARIABLES:                                       
     int     Nm                 number of molecules 
     int     M                  number of cells in each direction  
     int     NCell              total number of cells (M**3)       
     int     list[N]            linked list of atoms               
     int     head(NCell)        head of chain for each cell        
     COORD_TYPE  rx[N],ry[N],
              rz[N]             positions                          
     COORD_TYPE  rcut           the cutoff distance for the force  
     COORD_TYPE* sigma          potential length parameters (it is an array!)
     USAGE:                                                       
     
     Each atom is sorted into one of the M**3 small cells.        
     The first atom in each cell is placed in the head array.     
     Subsequent atoms are placed in the linked list array.        
     Atom coordinates are assumed to be between -0.5 and +0.5.    
     The routine is called every timestep before the force routine and
     in a Gear Predictor-Corrector after the predictor */
  COORD_TYPE  celli, cell, rxc, ryc, rzc;
  int iCell, i, a, b;
  
  L = cbrt(Vol);
  invL = 1.0 / L;
 
  for(iCell=0; iCell < NCell; iCell++)
    {
      head[iCell] = -1; /* -1 means end of list, that is no more particles */
    }
  pool;

  celli = (COORD_TYPE) M;
  cell  = 1.0 / celli;

  for(a=0; a < NA; a++)
    {
      for(b = 0; b < NA; b++)
	{
	  if (L * cell < (OprogStatus.rNebrShell + rcut * Oparams.sigma) ) 
	    /* avoid this check before, put outside !!!!!!!!!*/
	    {
	      printf("critical cutoff: %f\n", OprogStatus.rNebrShell +
		     rcut * Oparams.sigma);
	      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
		    "Cell size too small for cutoff",
		    NULL);
	      exit(-1);		   
	    }
	}
    }

  /* Sort all atoms */
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  /* rxc[][], ryc[][], rzc[] are the corrected coordinates, 
	     that is coordinates restricted to the central box */
	  
	  /* Corrected coordinates, that is restricted to the central box.
	     These coordinates are relative to all the particles actually 
	     in the central box.
	     Finally note that the rx, ry, rz coordinates are uncorrected 
	     because we are also interested in calculating transport 
	     coefficents */
	  rxc = rx[a][i] - L * rint(invL * rx[a][i]);
	  ryc = ry[a][i] - L * rint(invL * ry[a][i]);
	  rzc = rz[a][i] - L * rint(invL * rz[a][i]);
	  rxc /= L;
	  ryc /= L;
	  rzc /= L;
	  iCell = ( (int) (( rxc + 0.5 ) * celli) )
	    + ( (int) ((ryc + 0.5) * celli) ) * M
	    + ( (int) ((rzc + 0.5) * celli) ) * M * M;
	  /* At the atom (a,i) is associated the list element 'NA*i + a + 1',
	     and vice-versa if list[i] = j then j is the atoms:
	     (j % NA, j / NA) */ 

	  /* insert the atom (a,i) at the head of the linked list */
	  list[NA*i + a] = head[iCell]; /* Head[iCell] becomes the next
					   of 'NA*i + a' */ 
	  head[iCell] = NA*i + a;        
	  /* Now Head[iCell] is 'NA*i + a', in this way we have 
	     add 'NA*i + a' at the Head of the linked list.*/
	}
    }
}

/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut) 
{
  int a, b, i, j;
  COORD_TYPE rcutab, rcutabSq; 
  COORD_TYPE rrNebr;
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;
  L = cbrt(Vol);
  invL = 1.0  / L;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	rx_old[a][i] = rx[a][i];
	ry_old[a][i] = ry[a][i];
	rz_old[a][i] = rz[a][i];
      }

  /* useful ab-constants inside OUTER LOOP below */
  rcutab = rCut * Oparams.sigma;
  rcutabSq = Sqr(rcutab);
  rrNebr = Sqr(rcutab + OprogStatus.rNebrShell);
  if (rrNebr > Sqr(L / 2.0))
    {
      printf("(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
	     sqrt(rrNebr), L/2.0);
      exit(-1);
    }

  nebrTabLen = 0;
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];
	  /* INNER LOOP BEGINS */
	  for(j = i; j < Nm; j++) 
	    {
      	      for(b=0; b < NA; b++)  /* b >= a because of symmetry */
		/* + 2 because you must remeber that really the all indices 
		   (a, b, i, j) start from 0 ... */
		/*   i > j because of 3rd law */
		{
		  if (i==j && b < a+2)
		    continue;
		  rxab = rxa - rx[b][j]; /* distance between two atomes */
		  ryab = rya - ry[b][j];
		  rzab = rza - rz[b][j];
		  rxab = rxab - L * rint(rxab * invL);      /* minimum image */
		  ryab = ryab - L * rint(ryab * invL);
		  rzab = rzab - L * rint(rzab * invL);
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		  /*printf("rabSq: %f rrNebr: %f\n", rabSq, rrNebr);*/		  
		  if ( rabSq < rrNebr)/* 'rcut' is the cutoff for V */
		    {
		      if (nebrTabLen >= nebrTabMax)
			{
			  printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
				 nebrTabLen);
			  printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
			  printf("ERROR: Neighbourlist overflow!\n");
			  exit(-1);
			}
		      nebrTab[0][nebrTabLen] = NA*i + a;
		      nebrTab[1][nebrTabLen] = NA*j + b;
		      ++nebrTabLen; /* Increment table element counter */
		    }
		}
	      /* INNER LOOP ENDS */
	    }
	}
    }
  printf("step N. %d: nebrTabLen: %d\n", Oparams.curStep, nebrTabLen);
}

/* ======================== >>> BuildNebrList <<< ========================== */
void BuildNebrList(int Nm, COORD_TYPE rCut) 
{
  int a, b, i, j;
  COORD_TYPE rcutab, rcutabSq; 
  COORD_TYPE rrNebr;
  int  iCell, jCell0, jCell, nabor, hd, lst;
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;
  
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < NA; a++)
      {
	rx_old[a][i] = rx[a][i];
	ry_old[a][i] = ry[a][i];
	rz_old[a][i] = rz[a][i];
      }
  /* useful ab-constants inside OUTER LOOP below */
  rcutab = rCut * Oparams.sigma;
  rcutabSq = Sqr(rcutab);
  rrNebr = Sqr(rcutab + OprogStatus.rNebrShell);

  L = cbrt(Vol);
  invL = 1.0  / L;
  nebrTabLen = 0;
  /* Every process build its own neighbour list */
  for(iCell=0; iCell < NCell; iCell++)
    {
      /* LOOP OVER ALL MOLECULES IN THE CELL */
      hd = head[iCell]; 

      while(hd > -1)
	{
	  i = hd / NA;
	  a = hd % NA;
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];
	  
	  /* LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL */
	  lst = list[hd]; 
				
	  while(lst > -1)
	    /* Atoms in the same molecule don't interact */
	    { 
	      j = lst / NA;
	      b = lst % NA;
	 
	      if (j == i)
		/* atoms in the same molecule don't interact */
		{
		  lst = list[lst];
		  continue;
		}
	    
	      rxab = rxa - rx[b][j]; /* distance between two atomes */
	      ryab = rya - ry[b][j];
	      rzab = rza - rz[b][j];
	      
	      rxab = rxab - L * rint(invL * rxab);      /* minimum image */
	      ryab = ryab - L * rint(invL * ryab);
	      rzab = rzab - L * rint(invL * rzab);
	      
	      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
	     
	      if ( rabSq < rrNebr && abs(b - a) > 1) /* 'rcut' is the cutoff for V */
		{
		  if (nebrTabLen >= nebrTabMax)
		    {
		      printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			     nebrTabLen);
		      printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
		      printf("ERROR: Neighbourlist overflow!\n");
		      exit(-1);
		    }
		  nebrTab[0][nebrTabLen] = NA*i + a;
		  nebrTab[1][nebrTabLen] = NA*j + b;
		  ++nebrTabLen; /* Increment table element counter */
		}
	      
	      lst = list[lst]; /* next atom in the list */
	    }
	  /* LOOP OVER NEIGHBOURING CELLS */
	  jCell0 = 13 * iCell;
	  for(nabor=0; nabor < 13; nabor++)
	    {
	      jCell = map[jCell0 + nabor];
	      lst = head[jCell];
	      while(lst != -1) 
		{
		  j = lst / NA;
		  b = lst % NA;
		  if ( j == i) /* atoms in the same molecule don't interact */
		    {
		      lst = list[lst];
		      continue;
		    }
		 
		  rxab = rxa - rx[b][j]; /* distance between two atomes */
		  ryab = rya - ry[b][j];
		  rzab = rza - rz[b][j];
		  
		  rxab = rxab - L * rint(invL * rxab);    /* minimum image */
		  ryab = ryab - L * rint(invL * ryab);
		  rzab = rzab - L * rint(invL * rzab);
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		  
		  if ( rabSq < rrNebr && abs(b - a) > 1)/* 'rcut' is the cutoff for V */
		    {
		      if (nebrTabLen >= nebrTabMax)
			{
			  printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax,
				 nebrTabLen);
			  printf("ERROR: Neighbourlist overflow!\n");
			  printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
			  exit(-1);
			}
		      nebrTab[0][nebrTabLen] = NA*i + a;
		      nebrTab[1][nebrTabLen] = NA*j + b;
		      ++nebrTabLen;
		    }
		  
		  lst = list[lst];  /* next atom*/
		
		}
	    }
	  hd = list[hd]; /* next atom */
	}
    }
  /* OUTER LOOP ENDS */
}

/* ========================== >>> CoM <<< ================================== */
void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm)
{
  double *m, Mtot;
  int j, n;

  L = cbrt(Vol);
  invL = 1.0 / L;
  Mtot = 0.0;
  m = Oparams.m;
  for (j=0; j < NA; j++)
    Mtot += m[j];
  *rxcm = *rycm = *rzcm = 0;
  for (n = 0; n < NA; n++)
    { 
      *rxcm += rx[n][i]*m[n];
      *rycm += ry[n][i]*m[n];
      *rzcm += ry[n][i]*m[n];
    }
  *rxcm /= Mtot;
  *rycm /= Mtot;
  *rzcm /= Mtot;

}

/* ========================== >>> ComV <<< ================================ */
void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm)
{
  COORD_TYPE* m, Mtot;
  int j, n;

  /*L = cbrt(Vol);
  invL = 1.0 / L;*/

  m = Oparams.m;
  Mtot=0.0;
  for (j=0; j < NA; j++)
    Mtot += m[j];
  
  *vxcm = *vycm = *vzcm = 0; 
  for (n=0; n < NA; n++)
    {
      *vxcm += vx[n][i]*m[n];
      *vycm += vy[n][i]*m[n];
      *vzcm += vz[n][i]*m[n];
    }
  *vxcm /= Mtot;
  *vycm /= Mtot;
  *vzcm /= Mtot;
}
/* =========================== >>> kinet <<< ============================== */
void kinet(int Nm, COORD_TYPE** velx, COORD_TYPE** vely, COORD_TYPE** velz,
	   COORD_TYPE VOL1)
{
  /* DESCRIPTION:
     Calculate the kinetic energy of the sample, this is done for 
     correcting the Hoover friction coefficent with the appropriate value
     of the temperature (i.e. the predicted one) */
  int i, a;
  COORD_TYPE dlnV, px, py, pz;
#if 0
  double vmed;
#endif
  /* NOTE: K is not shared so this loop is done by father and child */
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  dlnV = VOL1 / Vol;
  dlnV /= 3.0;
#if 0
  vmed = 0;
#endif
  for(i=0; i < Nm; i++)
    {
      for(a=0; a < NA; a++)
	{
	  px = velx[a][i] - dlnV * Rx[i];
	  py = vely[a][i] - dlnV * Ry[i];
	  pz = velz[a][i] - dlnV * Rz[i];
	  K = K + Oparams.m[a] * (Sqr(px) + Sqr(py) + Sqr(pz));
	}
#if 0
      vmed = vmed + sqrt(Sqr((velx[0][i]+velx[1][i]+velx[2][i])/3.0) + 
			 Sqr(Sqr((vely[0][i]+vely[1][i]+vely[2][i])/3.0)) + 
			 Sqr(Sqr((velz[0][i]+velz[1][i]+velz[2][i])/3.0)) );
#endif
    }
#if 0
  printf("steplength: %f vmed=%f\n", Oparams.steplength, vmed / Oparams.parnum);
#endif
  /* Kp is the kinetic energy calculated using 'predicted' velocities at time
     t (v(t)) */
  K *= 0.5;
}

void checkNebrRebuild(void)
{
  int i, a, Nm = Oparams.parnum;
  double rNebrShellSq;
  /*double norm, vv, vvMax = 0.0;
  double RCMx, RCMy, RCMz, VCMx, VCMy, VCMz;*/

  rNebrShellSq = Sqr(0.5*OprogStatus.rNebrShell);
  for(i=0; i < Nm && !nebrNow ; i++)
    {
      /*CoM(i, &RCMx, &RCMy, &RCMz);
	CoMV(i, &VCMx, &VCMy, &VCMz);*/
      for(a=0; a < NA; a++)
	{
	  /* usando la velocità angolare calcolata sopra, 
	   * vengono calcolate le velocità di tutti gli atomi */
	  /* vectProd(omegax, omegay, omegaz, rx[a][i], ry[a][i], rz[a][i],
	     &vax, &vay, &vaz); */
	  if (Sqr(rx[a][i]-rx_old[a][i])+Sqr(ry[a][i]-ry_old[a][i])+
	      Sqr(rz[a][i]-rz_old[a][i]) > rNebrShellSq)
	    {
	      printf("STEP N. %d rebuilding neghbourlists!\n", Oparams.curStep);
	      printf("(%f,%f,%f)-(%f,%f,%f)\n", rx[a][i], ry[a][i],rz[a][i],
		     rx_old[a][i], ry_old[a][i],rz_old[a][i]);
	      nebrNow=1;
	      break;
	    } 
	    /*
	  vv = Sqr(VCMx) + Sqr(VCMy) + Sqr(VCMz);
	  if (vv > vvMax) 
	    vvMax = vv;*/
	}
    }
  /* dispHi = dispHi + sqrt(vvMax) * Oparams.steplength;*/
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  /*if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;*/

}
const double inv4pieps0=1.43999;

/* ============================ >>> force <<< ==============================*/
void LJForce(int Nm, double rcut)
{
  /* ======================== >>>LOCAL VARIABLES <<< ====================== */
  int a, b, i, j, ncut, grow;
  /*int ncut[NA][NA];*/
  COORD_TYPE rcutab, rcutabSq, epsab4, sigmaSq; 
  COORD_TYPE rxab, ryab, rzab, rabSq, fxab, fyab, fzab, rab;
  COORD_TYPE srab2, srab6, srab12, fab, vabhc, wab, vabyu, vab;
  COORD_TYPE Fxa, Fya, Fza, rxa, rya, rza;
  COORD_TYPE vabCut, Vcab;
  COORD_TYPE Vab, Wab, dvdr;
  COORD_TYPE rho, DRmx, DRmy, DRmz;
  /* Local variables to implement linked list */
  int  n, nebrTab0, nebrTab1;
  COORD_TYPE Wmyx, Wmzy, Wmxz, kD;
#ifdef NM_SPHERE
  double vabNN, vabMM;
#endif
  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to "ab-quantities" for all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */

  rcutab = rcut * Oparams.sigma;
  rcutabSq = Sqr(rcutab);
  sigmaSq = Sqr(Oparams.sigma);
  epsab4 = 4.0 * Oparams.epsilon;
  Vab = 0.0;
  Wab = 0.0;
  L = cbrt(Vol);
  invL = 1.0  / L;
  ncut = 0;
  /* initialize forces vector */
  for (i=0;i < Nm; i++) 
    {
#ifdef MOLPTENS      
      CoM(i, &Rmx[i], &Rmy[i], &Rmz[i]);
#endif
#if 0
      /* Calculate center of mass of all molecules */
      for(a=0; a < NA; a++)
	{
	  Fallx[a][i] = 0.0;
	  Fally[a][i] = 0.0;
	  Fallz[a][i] = 0.0;
	}
#endif
    }
  
  V = 0.0; /* potential energy */
  W = 0.0; /* virial function */
#ifdef ATPTENS
  Wxy = 0.0; /* virial off-diagonal terms of pressure tensor */
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;
#endif
#ifdef MOLPTENS
  Wm = 0.0;
  Wmxx = Wmyy = Wmzz = 0.0;
  Wmxy = Wmyz = Wmzx = 0.0;
  Wmyx = Wmzy = Wmxz = 0.0;
#endif
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
	  Fx[a][i] = 0.0;
	  Fy[a][i] = 0.0;
	  Fz[a][i] = 0.0;
	}
    }

  /* Loop over cells */
  Vc = 0.0;
  /*printf("nebrTabLen:%d\n", nebrTabLen);*/
  for (n=0; n < nebrTabLen; n++)
    {
      nebrTab0 = nebrTab[0][n]; 
      nebrTab1 = nebrTab[1][n];

      i = nebrTab0 / NA;
      a = nebrTab0 % NA;
      j = nebrTab1 / NA;
      b = nebrTab1 % NA;
      rxa = rx[a][i];
      rya = ry[a][i];
      rza = rz[a][i];
	  
      Fxa = Fx[a][i];
      Fya = Fy[a][i];
      Fza = Fz[a][i];
      
      rxab = rxa - rx[b][j]; /* distance between two atomes */
      ryab = rya - ry[b][j];
      rzab = rza - rz[b][j];
      
      rxab = rxab - L * rint(invL * rxab);      /* minimum image */
      ryab = ryab - L * rint(invL * ryab);
      rzab = rzab - L * rint(invL * rzab);
      
      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
      if (rabSq < rcutabSq )/* 'rcut' is the cutoff for V */
	{
	  /*rab   = sqrt(rabSq);*/
	  if (OprogStatus.grow)
	    srab2 = Sqr(sigmag[a][i])/rabSq; 
	  else
	    srab2 = sigmaSq / rabSq;
#if defined(SOFT_SPHERE)
	  vab = pow(srab2, ((double)Oparams.NN)/2.0);
  	  wab = ((double)Oparams.NN)*vab;
#elif defined(NM_SPHERE)
	  vabNN = pow(srab2, ((double)Oparams.NN)/2.0);
	  vabMM = -pow(srab2, ((double)Oparams.MM)/2.0);
	  vab = vabNN + vabMM;
  	  wab = ((double)Oparams.NN)*vabNN - ((double)Oparams.MM)*vabMM ;
#else
	  /* Lennard-Jones */
	  srab6   = srab2 * srab2 * srab2;
	  srab12  = Sqr(srab6);
	  vab     = srab12 - srab6;
	  /*vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);*/
	  wab     = 6.0*(vab + srab12);
#endif
#if 0
	  vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);
#endif
	  V = V + vab;
	  /* total potential between all a-b atoms pairs */
	  W = W + wab; 
	  /* NOTE: If you will use a shifted-force potential then 
	     calculate the force using that potential */
	  fab   = epsab4 * wab / rabSq;
#if 0
	  if (OprogStatus.grow && fabs(fab) > 1000)
	    fab = 1000;
#endif
	  /* force between two atoms */
	  fxab  = fab * rxab;         
	  fyab  = fab * ryab;
	  fzab  = fab * rzab;
	  /*printf("(%f,%f,%f)\n",fxab,fyab,fzab);*/
#ifdef ATPTENS
	  /* Virial off-diagonal terms of atomic pressure tensor */
	  Wxy += rxab * fyab;
	  Wyz += ryab * fzab;
	  Wzx += rzab * fxab;
	  Wxx += rxab * fxab;
	  Wyy += ryab * fyab;
	  Wzz += rzab * fzab;
#endif
	  /* Calculate all terms of molecular
	     pressure tensor */
#ifdef MOLPTENS	  
	  DRmx = (Rmx[i] - Rmx[j]);
	  DRmx = DRmx - L * rint(invL * DRmx);
	  DRmy = (Rmy[i] - Rmy[j]);
	  DRmy = DRmy - L * rint(invL * DRmy);
	  DRmz = (Rmz[i] - Rmz[j]);
	  DRmz = DRmz - L * rint(invL * DRmz);
	  
	  Wmxx += DRmx * fxab;
	  Wmyy += DRmy * fyab;
	  Wmzz += DRmz * fzab;
	  
	  Wmyx += DRmy * fxab;
	  Wmzy += DRmz * fyab;
	  Wmxz += DRmx * fzab;
	
	  Wmxy += DRmx * fyab;
	  Wmyz += DRmy * fzab;
	  Wmzx += DRmz * fxab;
#endif
	  Fxa   = Fxa + fxab;     /* total force acting on atom (a,i)*/
	  Fya   = Fya + fyab;
	  Fza   = Fza + fzab;
	  
	  Fx[b][j] = Fx[b][j] - fxab;  /* -fxab = fxba (3rd law) */
	  Fy[b][j] = Fy[b][j] - fyab;
	  Fz[b][j] = Fz[b][j] - fzab;
	  ++ncut;
	}
    
      Fx[a][i] = Fxa;
      Fy[a][i] = Fya;
      Fz[a][i] = Fza;
    }
  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
  rab   = rcut*Oparams.sigma;
  rabSq = Sqr(rab);
  srab2   = sigmaSq / rabSq;
#if defined(SOFT_SPHERE)
  vabCut = pow(srab2, Oparams.NN/2.0);
  //printf("NN: %d srab2:%f rab:%f vabCut: %f\n", Oparams.NN, srab2, rab, vabCut);
#elif defined(NM_SPHERE)
  vabCut = pow(srab2, Oparams.NN/2.0)-pow(srab2,Oparams.MM/2.0);
#else
  srab6 = srab2 * srab2 * srab2;
  srab12 = srab6 * srab6;
  vabCut = srab12 - srab6;
#endif
  /*vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);*/
  Vcab = ((double)ncut) * vabCut;
  /*printf("ncut[%d][%d]:%d\n Vcab:%f\n",a,b, ncut[a][b], vabCut);*/
  /* ncut[a][b] is the number of atoms pairs a-b within 
     rcutab[a][b] */ 
  W = epsab4 * W / 3.0;
  Vc = epsab4 * (V - Vcab); 
  V = epsab4 * V;
  /* MULTIPLY FOR ENERGY FACTORS */
#ifdef MOLPTENS
  Wm = Wmxx + Wmyy + Wmzz;
#endif
  /* NOTA: controllare se questo va effettivamente commentato!!!
     Wmxy = (Wmxy + Wmyx)/2.0;
     Wmyz = (Wmyz + Wmzy)/2.0;
     Wmzx = (Wmzx + Wmxz)/2.0;
     */
} 
