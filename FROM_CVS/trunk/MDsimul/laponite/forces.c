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
extern double **rallx, **rally, **rallz, **Fallx, **Fally, **Fallz,
  **rallx_old, **rally_old, **rallz_old, *atcharge; 
extern int **ncut;
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

  for(a=0; a < Oparams.nsites; a++)
    {
      for(b = 0; b < Oparams.nsites; b++)
	{
	  if (L * cell < (OprogStatus.rNebrShell + rcut * Oparams.lambdaD) ) 
	    /* avoid this check before, put outside !!!!!!!!!*/
	    {
	      printf("critical cutoff: %f\n", OprogStatus.rNebrShell +
		     rcut * Oparams.lambdaD);
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
      for(a=0; a < Oparams.nsites; a++)
	{
	  /* rxc[][], ryc[][], rzc[] are the corrected coordinates, 
	     that is coordinates restricted to the central box */
	  
	  /* Corrected coordinates, that is restricted to the central box.
	     These coordinates are relative to all the particles actually 
	     in the central box.
	     Finally note that the rx, ry, rz coordinates are uncorrected 
	     because we are also interested in calculating transport 
	     coefficents */
	  rxc = rallx[a][i] - L * rint(invL * rallx[a][i]);
	  ryc = rally[a][i] - L * rint(invL * rally[a][i]);
	  rzc = rallz[a][i] - L * rint(invL * rallz[a][i]);
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
	  list[Oparams.nsites*i + a] = head[iCell]; /* Head[iCell] becomes the next
					   of 'NA*i + a' */ 
	  head[iCell] = Oparams.nsites*i + a;        
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
    for (a = 0; a < Oparams.nsites; a++)
      {
	rallx_old[a][i] = rallx[a][i];
	rally_old[a][i] = rally[a][i];
	rallz_old[a][i] = rallz[a][i];
      }
 
  /* useful ab-constants inside OUTER LOOP below */
  rcutab = rCut * Oparams.lambdaD;
  rcutabSq = Sqr(rcutab);
  rrNebr = Sqr(rcutab + OprogStatus.rNebrShell);
#ifndef MD_EFFPOT
  if (rrNebr > Sqr(L / 2.0))
    {
      printf("(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
	     sqrt(rrNebr), L/2.0);
      exit(-1);
    }
#endif
#ifdef MD_EFFPOT
  rcutab = rCut * Oparams.lambdaD *100.0;
  rcutabSq = Sqr(rcutab);
  rrNebr = Sqr(rcutab + OprogStatus.rNebrShell);
#endif
  nebrTabLen = 0;
  for(i=0; i < Nm - 1; i++)
    {
      for(a=0; a < Oparams.nsites; a++)
	{
	  rxa = rallx[a][i];
	  rya = rally[a][i];
	  rza = rallz[a][i];
	  for(b=0; b < Oparams.nsites; b++)  /* b >= a because of symmetry */
	    {
	      /* INNER LOOP BEGINS */
	      for(j= i + 1; j < Nm; j++) 
		/* + 2 because you must remeber that really the all indices 
		   (a, b, i, j) start from 0 ... */
		/*   i > j because of 3rd law */
		{
		  rxab = rxa - rallx[b][j]; /* distance between two atomes */
		  ryab = rya - rally[b][j];
		  rzab = rza - rallz[b][j];
		  rxab = rxab - L * rint(rxab * invL);      /* minimum image */
		  ryab = ryab - L * rint(ryab * invL);
		  rzab = rzab - L * rint(rzab * invL);
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		  /*printf("rabSq: %f rrNebr: %f\n", rabSq, rrNebr);*/		  
		  if ( rabSq < rrNebr )/* 'rcut' is the cutoff for V */
		    {
		       if (nebrTabLen >= nebrTabMax)
			 {
			   printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
				  nebrTabLen);
			   printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
			   printf("ERROR: Neighbourlist overflow!\n");
			   exit(-1);
			 }
		       nebrTab[0][nebrTabLen] = Oparams.nsites*i + a;
		       nebrTab[1][nebrTabLen] = Oparams.nsites*j + b;
		       ++nebrTabLen; /* Increment table element counter */
		     
		    }
		}
	      /* INNER LOOP ENDS */
	    }
	}
    }
  //printf("step N. %d: nebrTabLen: %d\n", Oparams.curStep, nebrTabLen);
}

/* ======================== >>> BuildNebrList <<< ========================== */
void BuildNebrList(int Nm, COORD_TYPE rCut) 
{
  int a, b, i, j;
  COORD_TYPE rcutab, rcutabSq; 
  COORD_TYPE rrNebr;
  int  iCell, jCell0, jCell, nabor, hd, lst, nsites;
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;
  
  nsites= Oparams.nsites;
  for (i = 0; i < Oparams.parnum; i++)
    for (a = 0; a < Oparams.nsites; a++)
      {
	rallx_old[a][i] = rallx[a][i];
	rally_old[a][i] = rally[a][i];
	rallz_old[a][i] = rallz[a][i];
      }
  /* useful ab-constants inside OUTER LOOP below */
  rcutab = rCut * Oparams.lambdaD;
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
	  i = hd / nsites;
	  a = hd % nsites;
	  rxa = rallx[a][i];
	  rya = rally[a][i];
	  rza = rallz[a][i];
	  
	  /* LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL */
	  lst = list[hd]; 
				
	  while(lst > -1)
	    /* Atoms in the same molecule don't interact */
	    { 
	      j = lst / nsites;
	      b = lst % nsites;
	 
	      if (j == i)
		/* atoms in the same molecule don't interact */
		{
		  lst = list[lst];
		  continue;
		}
	    
	      rxab = rxa - rallx[b][j]; /* distance between two atomes */
	      ryab = rya - rally[b][j];
	      rzab = rza - rallz[b][j];
	      
	      rxab = rxab - L * rint(invL * rxab);      /* minimum image */
	      ryab = ryab - L * rint(invL * ryab);
	      rzab = rzab - L * rint(invL * rzab);
	      
	      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
	     
	      if ( rabSq < rrNebr )/* 'rcut' is the cutoff for V */
		{
		  if (nebrTabLen >= nebrTabMax)
		    {
		      printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			     nebrTabLen);
		      printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
		      printf("ERROR: Neighbourlist overflow!\n");
		      exit(-1);
		    }
		  nebrTab[0][nebrTabLen] = nsites*i + a;
		  nebrTab[1][nebrTabLen] = nsites*j + b;
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
		  j = lst / nsites;
		  b = lst % nsites;
		  if ( j == i) /* atoms in the same molecule don't interact */
		    {
		      lst = list[lst];
		      continue;
		    }
		 
		  rxab = rxa - rallx[b][j]; /* distance between two atomes */
		  ryab = rya - rally[b][j];
		  rzab = rza - rallz[b][j];
		  
		  rxab = rxab - L * rint(invL * rxab);    /* minimum image */
		  ryab = ryab - L * rint(invL * ryab);
		  rzab = rzab - L * rint(invL * rzab);
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		  
		  if ( rabSq < rrNebr )/* 'rcut' is the cutoff for V */
		    {
		      if (nebrTabLen >= nebrTabMax)
			{
			  printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax,
				 nebrTabLen);
			  printf("ERROR: Neighbourlist overflow!\n");
			  printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
			  exit(-1);
			}
		      nebrTab[0][nebrTabLen] = nsites*i + a;
		      nebrTab[1][nebrTabLen] = nsites*j + b;
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
  int j;

  L = cbrt(Vol);
  invL = 1.0 / L;
  Mtot = 0.0;
  m = Oparams.m;
  for (j=0; j < NA; j++)
    Mtot += m[j];
  *rxcm = (rx[0][i]*m[0] + rx[1][i]*m[1] + rx[2][i]*m[2]) / Mtot;
  *rycm = (ry[0][i]*m[0] + ry[1][i]*m[1] + ry[2][i]*m[2]) / Mtot;
  *rzcm = (rz[0][i]*m[0] + rz[1][i]*m[1] + rz[2][i]*m[2]) / Mtot;

}

/* ========================== >>> ComV <<< ================================ */
void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm)
{
  COORD_TYPE* m, Mtot;
  int j;

  /*L = cbrt(Vol);
  invL = 1.0 / L;*/

  m = Oparams.m;
  Mtot=0.0;
  for (j=0; j < NA; j++)
    Mtot += m[j];
    
  *vxcm = (vx[0][i]*m[0] + vx[1][i]*m[1] + vx[2][i]*m[2]) / Mtot;
  *vycm = (vy[0][i]*m[0] + vy[1][i]*m[1] + vy[2][i]*m[2]) / Mtot;
  *vzcm = (vz[0][i]*m[0] + vz[1][i]*m[1] + vz[2][i]*m[2]) / Mtot;
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
      for(a=0; a < Oparams.nsites; a++)
	{
	  /* usando la velocità angolare calcolata sopra, 
	   * vengono calcolate le velocità di tutti gli atomi */
	  /* vectProd(omegax, omegay, omegaz, rallx[a][i], rally[a][i], rallz[a][i],
	     &vax, &vay, &vaz); */
	  if (Sqr(rallx[a][i]-rallx_old[a][i])+Sqr(rally[a][i]-rally_old[a][i])+
	      Sqr(rallz[a][i]-rallz_old[a][i]) > rNebrShellSq)
	    {
	      printf("STEP N. %d rebuilding neghbourlists!\n", Oparams.curStep);
	      printf("(%f,%f,%f)-(%f,%f,%f)\n", rallx[a][i], rally[a][i],rallz[a][i],
		     rallx_old[a][i], rally_old[a][i],rallz_old[a][i]);
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
  int a, b, i, j;
  /*int ncut[NA][NA];*/
  COORD_TYPE rcutab, rcutabSq; 
  COORD_TYPE rxab, ryab, rzab, rabSq, fxab, fyab, fzab, rab;
  COORD_TYPE srab2, srab6, srab12, fab, vabhc, wab, vabyu, vab;
  COORD_TYPE Fxa, Fya, Fza, rxa, rya, rza;
  COORD_TYPE vabCut, Vcab;
  COORD_TYPE Vab, Wab, dvdr;
  COORD_TYPE rho, DRmx, DRmy, DRmz;
  /* Local variables to implement linked list */
  int  n, nebrTab0, nebrTab1;
  COORD_TYPE Wmyx, Wmzy, Wmxz, kD;

  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to "ab-quantities" for all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */

  rcutab = rcut * Oparams.lambdaD;
  rcutabSq = Sqr(rcutab);
  Vab = 0.0;
  Wab = 0.0;
  kD = 1.0/Oparams.lambdaD;
  L = cbrt(Vol);
  invL = 1.0  / L;
  for (a=0; a < Oparams.nsites; a++)
    for (b=0; b < Oparams.nsites; b++)
      ncut[a][b] = 0;
  /* initialize forces vector */
  for (i=0;i < Nm; i++) 
    {
#ifdef MOLPTENS      
      CoM(i, &Rmx[i], &Rmy[i], &Rmz[i]);
#endif
      /* Calculate center of mass of all molecules */
      for(a=0; a < Oparams.nsites; a++)
	{
	  Fallx[a][i] = 0.0;
	  Fally[a][i] = 0.0;
	  Fallz[a][i] = 0.0;
	}
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
  /* Loop over cells */
  Vc = 0.0;
  /*printf("nebrTabLen:%d\n", nebrTabLen);*/
  for (n=0; n < nebrTabLen; n++)
    {
      nebrTab0 = nebrTab[0][n]; 
      nebrTab1 = nebrTab[1][n];

      i = nebrTab0 / Oparams.nsites;
      a = nebrTab0 % Oparams.nsites;
      j = nebrTab1 / Oparams.nsites;
      b = nebrTab1 % Oparams.nsites;
      
      rxa = rallx[a][i];
      rya = rally[a][i];
      rza = rallz[a][i];
	  
      Fxa = Fallx[a][i];
      Fya = Fally[a][i];
      Fza = Fallz[a][i];
      
      rxab = rxa - rallx[b][j]; /* distance between two atomes */
      ryab = rya - rally[b][j];
      rzab = rza - rallz[b][j];
      
      rxab = rxab - L * rint(invL * rxab);      /* minimum image */
      ryab = ryab - L * rint(invL * ryab);
      rzab = rzab - L * rint(invL * rzab);
      
      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
      if ( rabSq < rcutabSq )/* 'rcut' is the cutoff for V */
	{
	  rab   = sqrt(rabSq);
  	  srab2   = 1.0 / rabSq;
	  srab6   = srab2 * srab2 * srab2;
	  /*srab12  = Sqr(srab6);*/
	  vabhc = Oparams.crep * srab6;
	  //printf("atcharge[%d]=%.15G atcharge[%d]=%.15G\n", a, atcharge[a], b, atcharge[b]);
	  vabyu = inv4pieps0*atcharge[a]*atcharge[b]*exp(-rab*kD)/(rab*Oparams.epsilon);
#if 0
	  if (rab< 0.96)
          printf("inv4pieps0:%f atcharge[a]:%f atcharge[b]:%f epsilon:%f 1/kD:%f crep:%f rab:%f (%d:%d)-(%d:%d)\n", 
		 inv4pieps0, atcharge[a],
		 atcharge[b], Oparams.epsilon, 1/kD, Oparams.crep, rab,
		 i,a,j,b);
//	  printf("(%d,%d)-(%d,%d) vabhc:%.15f vabyu:%.15f\n", i, a, j, b, vabhc, vabyu);
#endif
	  /* 0.1902 per avere eV */
	  /*vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);*/
	  vab = vabhc + vabyu;
	  wab     = 6.0*vabhc + vabyu * (1.0 + rab*kD);
#if 0
	  if (rab < 1.0)
	    {
	      int aa;
	      double raex, raey, raez;
	      printf("vabhc:%f vabyu:%f\n", vabhc, vabyu);
	      printf("(%d,%d)-(%d,%d) rab=%f\n", i, a, j, b, rab);
	      raex = rallx[a][i]- rallx[b][j]; 
	      raey = rally[a][i] - rally[b][j];
	      raez = rallz[a][i] - rallz[b][j];
	      printf("Distanza senaa MI: %f (%f,%f,%f)\n", sqrt(Sqr(raex)+Sqr(raey)+Sqr(raez)),
		     raex, raey, raez);
	      for (aa=0; aa < 3; aa++)
		{
		  printf("[vel] i=%d: (%f,%f,%f), j=%d (%f,%f,%f)\n", i,
			 vx[aa][i], vy[aa][i], vz[aa][i], 
			 j, vx[aa][j], vy[aa][j], vz[aa][j]);
		 printf("[pos] i=%d: (%f,%f,%f), j=%d (%f,%f,%f)\n", i,
			 rx[aa][i], ry[aa][i], rz[aa][i], 
			 j, rx[aa][j], ry[aa][j], rz[aa][j]);
		

		}
	    }
#endif
	  V = V + vab;
	  /* total potential between all a-b atoms pairs */
	  W = W + wab; 
	  /* NOTE: If you will use a shifted-force potential then 
	     calculate the force using that potential */
	  fab   = wab / rabSq;
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
	  
	  Fallx[b][j] = Fallx[b][j] - fxab;  /* -fxab = fxba (3rd law) */
	  Fally[b][j] = Fally[b][j] - fyab;
	  Fallz[b][j] = Fallz[b][j] - fzab;
	  ++ncut[a][b];
	}
    
      Fallx[a][i] = Fxa;
      Fally[a][i] = Fya;
      Fallz[a][i] = Fza;
    }
    if (Oparams.curStep==1)
  /* OUTER LOOP ENDS */
    {
	      FILE* f;
	      int ii;
	      double r;
	      f = fopen ("potentialRimFace.dat","w");
	      for (ii=0; ii < 300; ii++)
		{
		  r = 0.2 + ((2.5 - 0.2)/300)*ii;
		  fprintf(f, "%.15f %.15f\n", r, Oparams.crep / pow(r,6) +
			  inv4pieps0*atcharge[0]*atcharge[60]*exp(-r*kD)/(r*Oparams.epsilon));
		}
	      fclose(f);
	      f = fopen ("potentialRimRim.dat","w");
	      for (ii=0; ii < 300; ii++)
		{
		  r = 0.2 + ((5 - 0.2)/300)*ii;
		  fprintf(f, "%.15f %.15f\n", r, Oparams.crep / pow(r,6) +
			  inv4pieps0*atcharge[0]*atcharge[0]*exp(-r*kD)/(r*Oparams.epsilon));
		}
	      fclose(f);
	      f = fopen ("potentialFaceFace.dat","w");
	      for (ii=0; ii < 300; ii++)
		{
		  r = 0.2 + ((5 - 0.2)/300)*ii;
		  fprintf(f, "%.15f %.15f\n", r, Oparams.crep / pow(r,6) +
			  inv4pieps0*atcharge[60]*atcharge[60]*exp(-r*kD)/(r*Oparams.epsilon));
		}
	      fclose(f);

}

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
  Vcab = 0.0;
  for(a=0; a < Oparams.nsites; a++)
    {
      for(b=0; b < Oparams.nsites; b++) /* b >= a */
	{
	  rab   = rcut*Oparams.lambdaD;
	  rabSq = Sqr(rab);
	  srab2   = 1.0 / rabSq;
	  srab6   = srab2 * srab2 * srab2;
	  vabhc = Oparams.crep * srab6;
	  vabyu = inv4pieps0*atcharge[a]*atcharge[b]*exp(-rab*kD)/(rab*Oparams.epsilon);
	  /*vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);*/
	  vabCut = vabhc + vabyu;
	  Vcab += ncut[a][b] * vabCut;
	  /*printf("ncut[%d][%d]:%d\n Vcab:%f\n",a,b, ncut[a][b], vabCut);*/
	  /* ncut[a][b] is the number of atoms pairs a-b within 
	     rcutab[a][b] */ 
	}
    }
  Vc = V - Vcab; 
  /* MULTIPLY FOR ENERGY FACTORS */
  /*
     for(a = 0; a < NA; a++)
     {
     for (b = 0; b < NA; b++) 
     {
     V  = V + Vab[a][b]  * epsab4[a][b];
     Vc = Vc + Vcab[a][b] * epsab4[a][b];
     W  = W + Wab[a][b]  * epsab24[a][b] / 3.0;
     }
     }
    */
#ifdef MOLPTENS
  Wm = Wmxx + Wmyy + Wmzz;
#endif
  /* V, Vc and W  are local variables and Father and Child have their own 
     values, so we must sum these ones, putting the result in Father local
     variables *
  */ 
/* NOTA: controllare se questo va effettivamente commentato!!!
  Wmxy = (Wmxy + Wmyx)/2.0;
  Wmyz = (Wmyz + Wmzy)/2.0;
  Wmzx = (Wmzx + Wmxz)/2.0;
  */
} 


void ForceOn123(void)
{
  int i, a;
  for (i = 0; i < Oparams.parnum; i++)
    {
      Fx[0][i] = Fallx[0][i];
      Fy[0][i] = Fally[0][i];
      Fz[0][i] = Fallz[0][i];
      Fx[1][i] = Fallx[1][i];
      Fy[1][i] = Fally[1][i];
      Fz[1][i] = Fallz[1][i];
      Fx[2][i] = Fallx[2][i];
      Fy[2][i] = Fally[2][i];
      Fz[2][i] = Fallz[2][i];

      for (a = 3; a < Oparams.nsites; a++)
	{
	  Fx[0][i] += Fcoeff[0][a]*Fallx[a][i];
	  Fy[0][i] += Fcoeff[0][a]*Fally[a][i];
	  Fz[0][i] += Fcoeff[0][a]*Fallz[a][i];
	  Fx[1][i] += Fcoeff[1][a]*Fallx[a][i];
	  Fy[1][i] += Fcoeff[1][a]*Fally[a][i];
	  Fz[1][i] += Fcoeff[1][a]*Fallz[a][i];
	  Fx[2][i] += Fcoeff[2][a]*Fallx[a][i];
	  Fy[2][i] += Fcoeff[2][a]*Fally[a][i];
	  Fz[2][i] += Fcoeff[2][a]*Fallz[a][i];
	}
    }

}
