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
void  links(int Nm, COORD_TYPE rcut, COORD_TYPE sigab[NA][NA])
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
  
  /*  ZERO HEAD OF CHAIN ARRAY */
  ProcSync0();
  
  L = cbrt(Vol);
  invL = 1.0 / L;
 
  loop(iCell, 1, NCell)
    {
      head[iCell] = -1; /* -1 means end of list, that is no more particles */
    }
  pool;

  celli = (COORD_TYPE) M;
  cell  = 1.0 / celli;
  

  loop(a, 1, NA)
    {
      loop(b, 1, NA)
	{
	  //printf("L*cell: %f\n", L*cell);
	  if (L * cell < (OprogStatus.rNebrShell + rcut * sigab[a][b]) ) 
	    /* avoid this check before, put outside !!!!!!!!!*/
	    {
	      printf("critical cutoff: %f\n", OprogStatus.rNebrShell +
		     rcut * sigab[a][b]);
	      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
		    "Cell size too small for cutoff",
		    NULL);
	      exit(-1);		   
	    }
	}
    }

  /* Sort all atoms */
  loop(i, 1, Nm)
    {
      loop(a, 1, NA)
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
  pool;
}

/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]) 
{
  int a, b, i, j;
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE rrNebr[NA][NA];
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;

  ProcSync0();

  L = cbrt(Vol);
  invL = 1.0  / L;

  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a because of symmetry */
	{
	  /* useful ab-constants inside OUTER LOOP below */
	  rcutab[a][b] = rCut * sigab[a][b];
	  rcutabSq[a][b] = Sqr(rcutab[a][b]);
	  rrNebr[a][b] = Sqr(rcutab[a][b] + OprogStatus.rNebrShell);
	  if (rrNebr[a][b] > Sqr(L / 2.0))
	    {
	      printf("(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
		     sqrt(rrNebr[a][b]), L/2.0);
	      exit(-1);
	    }
	  //printf("sqrt(rrNebr[%d][%d]):%f\n", a, b, sqrt(rrNebr[a][b]));
	}
    }
  //doFather{printf("STEP: %d rebuilding nbr Lst\n",Oparams.curStep);};
  nebrTabLen = 0;
  loop(i, 1, Nm - 1)
    {
      loop(a, 1, NA)
	{
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];
	  loop(b, 1, NA)  /* b >= a because of symmetry */
	    {
	      /* INNER LOOP BEGINS */
	      loop(j, i + 2, Nm) 
		/* + 2 because you must remeber that really the all indices 
		   (a, b, i, j) start from 0 ... */
		/*   i > j because of 3rd law */
		{
		  rxab = rxa - rx[b][j]; /* distance between two atomes */
		  ryab = rya - ry[b][j];
		  rzab = rza - rz[b][j];
		  rxab = rxab - L * rint(rxab * invL);      /* minimum image */
		  ryab = ryab - L * rint(ryab * invL);
		  rzab = rzab - L * rint(rzab * invL);
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		  
		  if ( rabSq < rrNebr[a][b] )/* 'rcut' is the cutoff for V */
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
		       //doFather{printf("First Inter iCell: %d: (%d,%d)-(%d,%d)\n", 
				//       iCell,i,a,j,b);};
		       ++nebrTabLen; /* Increment table element counter */
		     
		    }
		}
	      /* INNER LOOP ENDS */
	    }
	}
    }
}

/* ======================== >>> BuildNebrList <<< ========================== */
void BuildNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]) 
{
  int a, b, i, j;
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE rrNebr[NA][NA];
  int  iCell, jCell0, jCell, nabor, hd, lst;
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;

  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a because of symmetry */
	{
	  /* useful ab-constants inside OUTER LOOP below */
	  rcutab[a][b] = rCut * sigab[a][b];
	  rcutabSq[a][b] = Sqr(rcutab[a][b]);
	  rrNebr[a][b] = Sqr(rcutab[a][b] + OprogStatus.rNebrShell);
	  //printf("sqrt(rrNebr[%d][%d]):%f\n", a, b, sqrt(rrNebr[a][b]));
	}
    }

  //doFather{printf("STEP: %d rebuilding nbr Lst\n",Oparams.curStep);};
  L = cbrt(Vol);
  invL = 1.0  / L;
  nebrTabLen = 0;
  /* Every process build its own neighbour list */
  loop(iCell, 1, NCell)
    {
      /* LOOP OVER ALL MOLECULES IN THE CELL */
      hd = head[iCell]; 

      while(hd > -1)
	{
	  i = hd / NA;
	  a = hd % NA;
	  //doFather{printf("cell: %d particle: (%d,%d)\n", iCell, i, a);};
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
	     
	      if ( rabSq < rrNebr[a][b] )/* 'rcut' is the cutoff for V */
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
		  //doFather{fprintf(stderr, "First Inter iCell: %d: (%d,%d)-(%d,%d)\n", 
		  //	  iCell,i,a,j,b);};
		  ++nebrTabLen; /* Increment table element counter */
		}
	      
	      lst = list[lst]; /* next atom in the list */
	    }
	  /* LOOP OVER NEIGHBOURING CELLS */
	  jCell0 = 13 * iCell;
	  loop(nabor, 1, 13)
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
		  
		  if ( rabSq < rrNebr[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      if (nebrTabLen >= nebrTabMax)
			{
			  printf("nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax,
				 nebrTabLen);
			  printf("ERROR: Neighbourlist overflow!\n");
			  printf("particles: (%d,%d)-(%d,%d)\n",i,a,j,b);
			  exit(-1);
			}
		      //doFather{fprintf(stderr, "jCell0: %d nabor:%d Second Inter cell:%d (%d,%d)-(%d,%d)\n",
		      //      jCell0, nabor, jCell, i,a,j,b);};
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
  pool;
  /* OUTER LOOP ENDS */
}




/* ========================== >>> Com <<< ================================== */
void CoMB(int a, int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm)
{
  /* DESCRIPTION:
     Calculate the center of mass position for the molecule to which belong
     the atom (a, i), this is used when we rescale not the molecule coordinates
     but the atoms ones */
  COORD_TYPE rx2, ry2, rz2; /* coords of the other particle */
  COORD_TYPE Dx, Dy, Dz;
  COORD_TYPE rxba, ryba, rzba;
  COORD_TYPE rxai, ryai, rzai, rxbi, rybi, rzbi;
  COORD_TYPE* m, Mtot, invMtot;
  int b, j;

  b = a + 1;
  if (b > 1) b = 0;
  L = cbrt(Vol);
  invL = 1.0 / L;
 
  m = Oparams.m;
  Mtot = 0.0;
  loop(j, 1, NA)
    {
      Mtot += m[j];
    }

  invMtot = 1.0 / Mtot;
  rxai = rx[a][i];
  ryai = ry[a][i];
  rzai = rz[a][i];
  rxbi = rx[b][i];
  rybi = ry[b][i];
  rzbi = rz[b][i];

  rxba = rxbi - rxai;
  ryba = rybi - ryai;
  rzba = rzbi - rzai;
  Dx = rxba - L * rint(invL * rxba);
  Dy = ryba - L * rint(invL * ryba); 
  Dz = rzba - L * rint(invL * rzba);
  rx2 = rxai + Dx;
  ry2 = ryai + Dy;
  rz2 = rzai + Dz;
  *rxcm = (m[a] * rxai + m[b] * rx2) * invMtot;
  *rycm = (m[a] * ryai + m[b] * ry2) * invMtot;
  *rzcm = (m[a] * rzai + m[b] * rz2) * invMtot;
}

/* ========================== >>> CoM <<< ================================== */
void CoM(int i, COORD_TYPE* rxcm, COORD_TYPE* rycm, COORD_TYPE* rzcm)
{
  COORD_TYPE* m, Mtot;

  L = cbrt(Vol);
  invL = 1.0 / L;
 
  m = Oparams.m;
  Mtot = 0.0;
  Mtot = m[0] + m[1];

  *rxcm = (m[0] * rx[0][i] + m[1] * rx[1][i]) / Mtot;
  *rycm = (m[0] * ry[0][i] + m[1] * ry[1][i]) / Mtot;
  *rzcm = (m[0] * rz[0][i] + m[1] * rz[1][i]) / Mtot;

}

/* ========================== >>> ComV <<< ================================ */
void CoMV(int i, COORD_TYPE* vxcm, COORD_TYPE* vycm, COORD_TYPE* vzcm)
{
  COORD_TYPE* m, Mtot;

  //L = cbrt(Vol);
  //invL = 1.0 / L;

  m = Oparams.m;
  Mtot = m[0] + m[1];
    

  *vxcm = (m[0] * vx[0][i] + m[1] * vx[1][i]) / Mtot;
  *vycm = (m[0] * vy[0][i] + m[1] * vy[1][i]) / Mtot;
  *vzcm = (m[0] * vz[0][i] + m[1] * vz[1][i]) / Mtot;

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

  ProcSync0();
  /* NOTE: K is not shared so this loop is done by father and child */
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  dlnV = VOL1 / Vol;
  dlnV /= 3.0;
 
  loop(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  px = velx[a][i] - dlnV * Rx[i];
	  py = vely[a][i] - dlnV * Ry[i];
	  pz = velz[a][i] - dlnV * Rz[i];
	  K = K + Oparams.m[a] * (Sqr(px) + Sqr(py) + Sqr(pz));
	}
    }

  /* Kp is the kinetic energy calculated using 'predicted' velocities at time
     t (v(t)) */
  K *= 0.5;
}

void checkNebrRebuild(void)
{
  int i, a, Nm = Oparams.parnum;
  COORD_TYPE vv, vvMax = 0.0;

  ProcSync0();
  loop(i, 1, Nm)
    {
      loop(a, 1, NA)
	{
	  vv = Sqr(vx[a][i]) + Sqr(vy[a][i]) + Sqr(vz[a][i]);
	}
      if (vv > vvMax) 
	vvMax = vv;
    }
  dispHi = dispHi + sqrt(vvMax) * Oparams.steplength;
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;

}
/* ======================= >>> sumForces <<< =============================== */
void sumForces(int Nm)
{
  /* DESCRIPTION:
     To avoid interferences the each process has its local Force array, in 
     which it stores its contribution to the total force acting on each atom.
     When the main loop in LJForce ends then we must add the two contributions,
     and this procedure do this */
  int i, a;

  ProcSync0(); /* Preliminar syncronization */

  doFather /* Put in the shared array of forces the father contribution */ 
    {
      loop(i, 1, Nm)
	{
	  loop(a, 1, NA)
	    {
	      Fx[a][i] = FxL[a][i]; /* F?S are shared variables */
	      Fy[a][i] = FyL[a][i];
	      Fz[a][i] = FzL[a][i];
	    }
	}
    }
  
  ProcSync0(); /* Syncronization, child waits father */
   
  doChild /* And now sum the Child contribution */
    { 
      loop(i, 1, Nm)
	{
	  loop(a, 1, NA)
	    {
	      Fx[a][i] += FxL[a][i]; /* F?S are shared variables */
	      Fy[a][i] += FyL[a][i];
	      Fz[a][i] += FzL[a][i];
	    }
	}
    }
}

/* ============================ >>> force <<< ==============================*/
void LJForce(int Nm, COORD_TYPE epsab[NA][NA], 
	     COORD_TYPE sigab[NA][NA], COORD_TYPE rcut)
{
  /* ======================== >>>LOCAL VARIABLES <<< ====================== */
  int a, b, i, j;
  int ncut[NA][NA];
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA], Epsab24;
  COORD_TYPE rxab, ryab, rzab, rabSq, fxab, fyab, fzab, rab;
  COORD_TYPE srab2, srab6, srab12, fab, vab, wab;
  COORD_TYPE Fxa, Fya, Fza, rxa, rya, rza;
  COORD_TYPE vabCut, Vcab[NA][NA];
  COORD_TYPE Vab[NA][NA], Wab[NA][NA], dvdr[NA][NA];
  COORD_TYPE rho, DRmx, DRmy, DRmz;
  /* Local variables to implement linked list */
  int  n, nebrTab0, nebrTab1;
  COORD_TYPE Wmyx, Wmzy, Wmxz;

  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */

  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a because of symmetry */
	{
	  /* useful ab-constants inside OUTER LOOP below */
	  rcutab[a][b] = rcut * sigab[a][b];
	  rcutabSq[a][b] = Sqr(rcutab[a][b]);
	  sigabSq[a][b] = Sqr(sigab[a][b]);
	  epsab4[a][b] = 4.0 * epsab[a][b];
	  epsab24[a][b] = 24.0 * epsab[a][b];
	  srab2   = sigabSq[a][b] / rcutabSq[a][b];
	  srab6   = srab2 * srab2 * srab2;
	  srab12  = Sqr(srab6);
	  dvdr[a][b] = epsab24[a][b] * (srab6 - 2.0 * srab12) / rcutab[a][b];
	  
	  /* initialize ab-variables */
	  ncut[a][b] = 0;
	  Vab[a][b] = 0.0;
	  Wab[a][b] = 0.0;
	  
	}
    }

  L = cbrt(Vol);
  invL = 1.0  / L;
  
  /* initialize forces vector */
  loop(i, 1, Nm) 
    {
#ifdef MOLPTENS      
      CoM(i, &Rmx[i], &Rmy[i], &Rmz[i]);
#endif
      /* Calculate center of mass of all molecules */
      loop(a, 1, NA)
	{
	  FxL[a][i] = 0.0;
	  FyL[a][i] = 0.0;
	  FzL[a][i] = 0.0;
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

  loopShrC(n, 1, nebrTabLen)
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
	  
      Fxa = FxL[a][i];
      Fya = FyL[a][i];
      Fza = FzL[a][i];
      
      Epsab24 = epsab24[a][b];
      rxab = rxa - rx[b][j]; /* distance between two atomes */
      ryab = rya - ry[b][j];
      rzab = rza - rz[b][j];
      
      rxab = rxab - L * rint(invL * rxab);      /* minimum image */
      ryab = ryab - L * rint(invL * ryab);
      rzab = rzab - L * rint(invL * rzab);
      
      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
	     
      if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
	{
	  //rab   = sqrt(rabSq);
	  
	  //printf("proc:%d step: %d (%d,%d)-(%d,%d) rab: %f\n", process,
	  // Oparams.curStep,
	  // i, a, j, b, sqrt(rabSq));
	 
	  srab2   = sigabSq[a][b] / rabSq;
	  srab6   = srab2 * srab2 * srab2;
	  srab12  = Sqr(srab6);
	  
	  vab     = srab12 - srab6;
	  //vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);
	  wab     = vab + srab12;
	  
	  Vab[a][b] = Vab[a][b] + vab;
	  
	  /* total potential between all a-b atoms pairs */
	  Wab[a][b]   = Wab[a][b] + wab; 
	  /* NOTE: If you will use a shifted-force potential then 
	     calculate the force using that potential */
	  fab   = wab / rabSq;
	  /* force between two atoms */
	  fxab  = fab * rxab * Epsab24;         
	  fyab  = fab * ryab * Epsab24;
	  fzab  = fab * rzab * Epsab24;
	  
#ifdef ATPTENS
	  /* Virial off-diagonal terms of atomic pressure tensor */
	  Wxy += rxab * fyab;
	  Wyz += ryab * fzab;
	  Wzx += rzab * fxab;
	  Wxx += rxab * fxab;
	  Wyy += ryab * fyab;
	  Wzz += rzab * fzab;
	  
#endif
	  /* Calculate all temrs of molecular
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
	  Fxa   = Fxa + fxab;     /* total force on an atom (a,i)*/
	  Fya   = Fya + fyab;
	  Fza   = Fza + fzab;
	  
	  FxL[b][j] = FxL[b][j] - fxab;  /* -fxab = fxba (3rd law) */
	  FyL[b][j] = FyL[b][j] - fyab;
	  FzL[b][j] = FzL[b][j] - fzab;
	  ++ncut[a][b];
	}
    
      FxL[a][i] = Fxa;
      FyL[a][i] = Fya;
      FzL[a][i] = Fza;
    }
  poolShrC;
  /*
  if ( (sc_[FATHER] >=  sc_[CHILD]) || (abs(sc_[FATHER] - sc_[CHILD]) > 1) ) 
    {
      sprintf(msgStrA, "sc_[FATHER]: %d  sc_[CHILD]:%d", sc_[FATHER], 
	      sc_[CHILD]); 
      mdMsg(ALL, NOSYS, NULL, "WARNING", NULL, 
	    "Processes overlap in the force routine.", 
	    msgStrA,
	    NULL);
    }
  */
 /* OUTER LOOP ENDS */

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a */
	{
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;
	  Vcab[a][b] = Vab[a][b] - ncut[a][b] * vabCut;
	  /* ncut[a][b] is the number of atoms pairs a-b within 
	     rcutab[a][b] */ 
	}
    }
  
  /* MULTIPLY FOR ENERGY FACTORS */
  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a and note that a varies from 0 to NA - 1 */
	{
	  V  = V + Vab[a][b]  * epsab4[a][b];
	  Vc = Vc + Vcab[a][b] * epsab4[a][b];
	  W  = W + Wab[a][b]  * epsab24[a][b] / 3.0;
	}
    }
#ifdef MOLPTENS
  Wm = Wmxx + Wmyy + Wmzz;
#endif
  /* V, Vc and W  are local variables and Father and Child have their own 
     values, so we must sum these ones, putting the result in Father local
     variables *
  */ 

#ifndef NO_PARALLEL_CODE
  sumAllProcMA(&V, &Vc, &W, NULL);
#ifdef ATPTENS
  /* Atomic virial contribute of pressure tensor */
  sumAllProcMA(&Wxy, &Wyz, &Wzx, &Wxx, &Wyy, &Wzz, NULL);

#endif
#ifdef MOLPTENS
  /* terms of molecular pressure tensor */ 
  sumAllProcMA(&Wm, &Wmxy, &Wmyz, &Wmzx, &Wmzx, &Wmyx,
	       &Wmzy, &Wmxz, &Wmxx, &Wmyy, &Wmzz, NULL);
#endif
  /* sumAllProc(FATHER, &Wm);
     sumAllProc(FATHER, &Wmxy);
     sumAllProc(FATHER, &Wmyz);
     sumAllProc(FATHER, &Wmzx);
     sumAllProc(FATHER, &Wmyx);
     sumAllProc(FATHER, &Wmzy);
     sumAllProc(FATHER, &Wmxz);*/
#endif
  Wmxy = (Wmxy + Wmyx)/2.0;
  Wmyz = (Wmyz + Wmzy)/2.0;
  Wmzx = (Wmzx + Wmxz)/2.0;

  /*sumAllProc(FATHER, &Wmxx);
    sumAllProc(FATHER, &Wmyy);
    sumAllProc(FATHER, &Wmzz);*/

  sumForces(Nm); /* sum the contributions of the two processes to forces acting
		    on atoms, and put the result in the shared forces arrays
		    Fx[][], Fy[][], Fz[][] */
}  
