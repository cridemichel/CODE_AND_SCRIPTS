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

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, WC; 
#if 0
  T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, WCxy, WCyz, WCzx, 
  WCxx, WCyy, WCzz, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Wm, Wmxx, Wmyy, Wmzz, 
  Wmxy, Wmyz, Wmzx, Pmxx, Pmyy, Pmzz, Pmxy, Pmyz, Pmzx, T1mxy, 
  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
  T1myz, T1mzx, T1mxx, T1myy, T1mzz; 
#endif
extern COORD_TYPE Mtot;
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax;
extern int **nebrTabInteracting;
/* ================================= */

extern char TXTA[10][MSG_LEN];
extern char TXT[MSG_LEN];

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
void  links(int Nm, COORD_TYPE rcut, COORD_TYPE sigma)
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
  int iCell, i;
  
  /*  ZERO HEAD OF CHAIN ARRAY */
  ProcSync0();
  
  L = Oparams.L;
  invL = 1.0 / L;
 
  loop(iCell, 1, NCell)
    {
      head[iCell] = -1; /* -1 means end of list, that is no more particles */
    }
  pool;

  celli = (COORD_TYPE) M;
  cell  = 1.0 / celli;
  

  if (L * cell < (OprogStatus.rNebrShell + rcut * sigma) ) 
    /* avoid this check before, put outside !!!!!!!!!*/
    {
     
      sprintf(TXT, "critical cutoff: %f\n", OprogStatus.rNebrShell +
	     rcut * sigma);
      mdPrintf(STD, TXT, NULL);
      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
	    "Cell size too small for cutoff",
	    NULL);
      exit(-1);		   
    }



  /* Sort all atoms */
  loop(i, 1, Nm)
    {
      /* rxc[][], ryc[][], rzc[] are the corrected coordinates, 
	 that is coordinates restricted to the central box */
      
      /* Corrected coordinates, that is restricted to the central box.
	 These coordinates are relative to all the particles actually 
	 in the central box.
	 Finally note that the rx, ry, rz coordinates are uncorrected 
	 because we are also interested in calculating transport 
	 coefficents */
      rxc = rx[i] - L * rint(invL * rx[i]);
      ryc = ry[i] - L * rint(invL * ry[i]);
      rzc = rz[i] - L * rint(invL * rz[i]);
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
      list[i] = head[iCell]; /* Head[iCell] becomes the next
				       of 'NA*i + a' */ 
      head[iCell] = i;        
      /* Now Head[iCell] is 'NA*i + a', in this way we have 
	 add 'NA*i + a' at the Head of the linked list.*/
    }
  pool;
}
/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void checkNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigma) 
{
  int i, j, n;
  COORD_TYPE rcut, rcutSq; 
  COORD_TYPE rrNebr;
  COORD_TYPE rxi, ryi, rzi, rijSq, rxij, ryij, rzij;

  L = Oparams.L;
  invL = 1.0  / L;

  /* useful ab-constants inside OUTER LOOP below */
  rcut = rCut * sigma;
  rcutSq = Sqr(rcut);
  rrNebr = Sqr(rcut + OprogStatus.rNebrShell);
  if (rrNebr > Sqr(L / 2.0))
    {
      sprintf(TXT, "(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
	     sqrt(rrNebr), L/2.0);
      mdPrintf(STD, TXT, NULL);
      exit(-1);
    }
  for(i=0; i < Nm; i++)
    {
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i];
      /* INNER LOOP BEGINS */
      for(j=0; j < Nm; j++) 
	{
	  if (j == i)
	    continue;
	  rxij = rxi - rx[j]; /* distance between two atomes */
	  ryij = ryi - ry[j];
	  rzij = rzi - rz[j];
	  rxij = rxij - L * rint(rxij * invL);      /* minimum image */
	  ryij = ryij - L * rint(ryij * invL);
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  
	  if ( rijSq < rrNebr )/* 'rcut' is the cutoff for V */
	    {
	      for (n=0; n < nebrTabLen[i]; n++)
		{
		  if (j == nebrTab[i][n])
		    {
		      return;
		    }
		}
	      printf ("particella che puo' urtare non nella lista!\n");
	      exit(-1);
	    }
	}
	      /* INNER LOOP ENDS */
    }
}

/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigma) 
{
  int i, j;
  COORD_TYPE rcut, rcutSq; 
  COORD_TYPE rrNebr;
  COORD_TYPE rxi, ryi, rzi, rijSq, rxij, ryij, rzij;

  L = Oparams.L;
  invL = 1.0  / L;

  /* useful ab-constants inside OUTER LOOP below */
  rcut = rCut * sigma;
  rcutSq = Sqr(rcut);
  rrNebr = Sqr(rcut + OprogStatus.rNebrShell);
  if (rrNebr > Sqr(L / 2.0))
    {
      sprintf(TXT, "(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
	     sqrt(rrNebr), L/2.0);
      mdPrintf(STD, TXT, NULL);
      exit(-1);
    }
  for(i=0; i < Nm; i++)
    {
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i];
      nebrTabLen[i] = 0;
      /* INNER LOOP BEGINS */
      for(j=0; j < Nm; j++) 
	{
	  if (j == i)
	    continue;
	  rxij = rxi - rx[j]; /* distance between two atomes */
	  ryij = ryi - ry[j];
	  rzij = rzi - rz[j];
	  rxij = rxij - L * rint(rxij * invL);      /* minimum image */
	  ryij = ryij - L * rint(ryij * invL);
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  if ( rijSq < rrNebr )/* 'rcut' is the cutoff for V */
	    {
	      if (nebrTabLen[i] >= OprogStatus.nebrTabFac)
		{
		  sprintf(TXTA[0], "nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			 nebrTabLen[i]);
		  sprintf(TXTA[1], "particles: (%d)-(%d)\n",i, j);
		  sprintf(TXTA[2], "ERROR: Neighbourlist overflow!\n");
		  mdPrintf(STD, TXTA[0], TXTA[1], TXTA[2], NULL);
		  exit(-1);
		}
	      nebrTab[i][nebrTabLen[i]] = j;
	      ++nebrTabLen[i];
	    }
	}
      /* INNER LOOP ENDS */
    }
}


/* ======================== >>> BuildNebrList <<< ========================== */
void BuildNebrList(int Nm, COORD_TYPE rCut, COORD_TYPE sigma) 
{
  int i, j;
  COORD_TYPE rcut, rcutSq; 
  COORD_TYPE rrNebr;
  int  iCell, jCell0, jCell, nabor, hd, lst;
  COORD_TYPE rxi, ryi, rzi, rijSq, rxij, ryij, rzij;

  /* useful ab-constants inside OUTER LOOP below */
  rcut = rCut * sigma;
  rcutSq = Sqr(rcut);
  rrNebr = Sqr(rcut + OprogStatus.rNebrShell);
  L = Oparams.L;
  invL = 1.0  / L;
  /* Every process build its own neighbour list */
  for(iCell=0; iCell < NCell; iCell++)
    {
      /* LOOP OVER ALL MOLECULES IN THE CELL */
      hd = head[iCell]; 

      while(hd > -1)
	{
	  i = hd;
	  nebrTabLen[i] = 0;
	  rxi = rx[i];
	  ryi = ry[i];
	  rzi = rz[i];
	  
	  /* LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL */
	  lst = list[hd]; 
				
	  while(lst > -1)
	    /* Atoms in the same molecule don't interact */
	    { 
	      j = lst;
	      	 
	      if (j == i)
		/* atoms in the same molecule don't interact */
		{
		  lst = list[lst];
		  continue;
		}
	    
	      rxij = rxi - rx[j]; /* distance between two atomes */
	      ryij = ryi - ry[j];
	      rzij = rzi - rz[j];
	      
	      rxij = rxij - L * rint(invL * rxij);      /* minimum image */
	      ryij = ryij - L * rint(invL * ryij);
	      rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	     
	      if ( rijSq < rrNebr)/* 'rcut' is the cutoff for V */
		{
		  if (nebrTabLen[i] >= OprogStatus.nebrTabFac)
		    {
		      sprintf(TXTA[0], "nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			     nebrTabLen[i]);
		      sprintf(TXTA[1], "particles: (%d)-(%d)\n", i, j);
		      sprintf(TXTA[2], "ERROR: Neighbourlist overflow!\n");
		      mdPrintf(STD, TXTA[0], TXTA[1], TXTA[2], NULL);
		      exit(-1);
		    }
		  nebrTab[i][nebrTabLen[i]] = j;
		  ++nebrTabLen[i]; /* Increment table element counter */
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
		  j = lst;
		  if ( j == i) /* atoms in the same molecule don't interact */
		    {
		      lst = list[lst];
		      continue;
		    }
		 
		  rxij = rxi - rx[j]; /* distance between two atomes */
		  ryij = ryi - ry[j];
		  rzij = rzi - rz[j];
		  
		  rxij = rxij - L * rint(invL * rxij);    /* minimum image */
		  ryij = ryij - L * rint(invL * ryij);
		  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
		  
		  if ( rijSq < rrNebr)/* 'rcut' is the cutoff for V */
		    {
		      if (nebrTabLen[i] >= OprogStatus.nebrTabFac)
			{
			  sprintf(TXTA[0], "nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax,
				 nebrTabLen[i]);
			  sprintf(TXTA[1], "ERROR: Neighbourlist overflow!\n");
			  sprintf(TXTA[2], "particles: (%d)-(%d)\n",i,j);
			  mdPrintf(STD, TXTA[0], TXTA[1], TXTA[2], NULL);
			  exit(-1);
			}
		      nebrTab[i][nebrTabLen[i]] = j;
		      ++nebrTabLen[i];
		    }
		  
		  lst = list[lst];  /* next atom*/
		
		}
	    }
	  hd = list[hd]; /* next atom */
	}
    }
}

void checkNebrRebuild(void)
{
  /*int i, Nm = Oparams.parnum;*/
  dispHi = dispHi + sqrt(3.0)*OprogStatus.drMax;
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  /*
     printf("0.5*fNebrS:%f dispHi:%f curStep:%d\n", 0.5*OprogStatus.rNebrShell, dispHi,
     Oparams.curStep);*/
  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;
}

/* =============================== >>> zeroForces <<< ====================== */
void zeroArrays(COORD_TYPE *arrx, COORD_TYPE *arry, COORD_TYPE *arrz, int N)
{
  int i;
  /* initialize forces vector */
  for(i = 0; i < N; i++) 
    {
      arrx[i] = 0.0;
      arry[i] = 0.0;
      arrz[i] = 0.0;
    }
}
#define MD_FCONT
#ifdef MPI
extern int my_rank;
#endif
/* ============================ >>> force <<< ==============================*/
void energyiF(int i, double rxi, double ryi, double rzi, COORD_TYPE epsilon, 
	      COORD_TYPE sigma, COORD_TYPE rcut, double* Vfi)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  int j, ncut;
  double sigmaSq, epsilon4, epsilon4_2N;
  double rxij, ryij, rzij, rijSq;
  double srij2, vij, L, invL;
  double rcutSig, rcutSq, A, B, C;
  /* Local variables to implement linked list */
  int  n;
  double rxic, ryic, rzic;
  double VfjOLD, VfjNEW, rzwall, mg, rziwSq, vw;
  double Fxi, Fyi, Fzi, vrij, Lz;
#ifdef MD_FCONT
  double dvdr, rij, rziw;
#endif
/* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */

  /* useful ab-constants inside OUTER LOOP below */
  rcutSig = rcut * sigma;
  rcutSq = Sqr(rcutSig);
  sigmaSq= Sqr(sigma);
  epsilon4 = 4.0 * epsilon;
  epsilon4_2N = epsilon4*2.0*((double)Oparams.N);
  L = Oparams.L;
  Lz =Oparams.Lz;
  invL = 1.0 / L;
  *Vfi = 0.0;
  Fxi = Fyi = Fzi = 0.0;
  mg = Oparams.m * Oparams.g;
#ifdef MD_FCONT
  dvdr = -pow(Sqr(sigmaSq/rcutSq), Oparams.N)/rcut;
#endif
#if 1
  for (n = 0; n < nebrTabLen[i]; n++)
    {
      j = nebrTab[i][n];
#else
  for (j=0; j < Oparams.parnum; j++)
    {
      if (j==i) 
	continue;
#endif
      rxij = rxi - rx[j]; /* distance between two atomes */
      ryij = ryi - ry[j];
      rzij = rzi - rz[j];
      
      rxij = rxij - L * rint(invL * rxij);      /* minimum image */
      ryij = ryij - L * rint(invL * ryij);
      rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
      
      if ( rijSq < rcutSq )/* 'rcut' is the cutoff for V */
	{
	  rij = sqrt(rijSq);
  	  srij2   = sigmaSq / rijSq;
	  vij = pow(srij2, ((double)Oparams.N));
	  /*printf("(%d,%d) rijSq: %f vij: %f\n", i, j, rijSq, vij);*/
#ifdef MD_FCONT
	  vrij = vij / rijSq + dvdr / rij;
#else
	  vrij = vij / rijSq;
#endif
	  Fxi += vrij * rxij;
	  Fyi += vrij * ryij;
	  Fzi += vrij * rzij;
	}
    }
  Fxi *= epsilon4_2N;
  Fyi *= epsilon4_2N;
  Fzi *= epsilon4_2N;
  Fzi -= mg;
  rzwall = -Lz*0.5 - Oparams.sigma*0.5;
  rziwSq = Sqr(rzi-rzwall);
#ifdef MD_FCONT
  dvdr = - pow(Sqr(sigmaSq/rcutSq), Oparams.Nw)/rcut;
#endif
  if (rziwSq < rcutSq)
    {
      /* -1.0 is a shift to have a zero potential when 
       * an atom touches the wall */
#ifdef MD_FCONT
      rziw = rzi-rzwall;
      vw = pow(sigmaSq/rziwSq, ((double)Oparams.Nw)); 
      /*printf("rziwSq: %f vw: %f Nw:%d\n", rziwSq, vw, Oparams.Nw);*/
      Fzi += 2.0*((double)Oparams.Nw)*(vw/rziwSq+dvdr/rziw)*(rzi-rzwall);
#else
      vw = pow(sigmaSq/rziwSq, ((double)Oparams.Nw)); 
      /*printf("rziwSq: %f vw: %f Nw:%d\n", rziwSq, vw, Oparams.Nw);*/
      Fzi += 2.0*((double)Oparams.Nw)*(vw/rziwSq)*(rzi-rzwall);
      *Vw = vw - pow(sigmaSq/rcutSq, Oparams.Nw);
#endif
    }	
  *Vfi += Oparams.alpha*(Sqr(Fxi)+Sqr(Fyi)+Sqr(Fzi));
  Fx[i] = Fxi;
  Fy[i] = Fyi;
  Fz[i] = Fzi;
  /*printf("DVlij: %f DVg: %f DVf: %f DVw:%f\n", *Vlji, *Vgi, *Vfi, *Vw );*/
} 

/* ============================ >>> force <<< ==============================*/
void energyi(int i, double rxi, double ryi, double rzi, COORD_TYPE epsilon, 
	    COORD_TYPE sigma, COORD_TYPE rcut, double *Vlji, double *Vgi, 
	    double* Vfi, double* Vw, int addnebr)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  int j, ncut;
  double sigmaSq, epsilon4, epsilon4_2N;
  double rxij, ryij, rzij, rijSq;
  double srij2, vij, L, invL;
  double rcutSig, rcutSq, A, B, C;
  /* Local variables to implement linked list */
  int  n;
  double rxic, ryic, rzic;
  double VfjOLD, VfjNEW, rzwall, mg, rziwSq, vw;
  double Fxi, Fyi, Fzi, vrij, Lz;
#ifdef MD_FCONT
  double dvdr, rij, rziw;
#endif
/* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */

  /* useful ab-constants inside OUTER LOOP below */
  rcutSig = rcut * sigma;
  rcutSq = Sqr(rcutSig);
  sigmaSq= Sqr(sigma);
  epsilon4 = 4.0 * epsilon;
  epsilon4_2N = epsilon4*2.0*((double)Oparams.N);
  L = Oparams.L;
  Lz =Oparams.Lz;
  invL = 1.0 / L;
  *Vlji = 0.0; /* potential energy */
  *Vgi = 0.0;
  *Vfi = 0.0;
  *Vw = 0.0;
  ncut = 0;
  Fxi = Fyi = Fzi = 0.0;
  mg = Oparams.m * Oparams.g;
#ifdef MD_FCONT
  dvdr = -pow(Sqr(sigmaSq/rcutSq), Oparams.N)/rcut;
#endif
#if 1
  for (n = 0; n < nebrTabLen[i]; n++)
    {
      j = nebrTab[i][n];
#else
  for (j=0; j < Oparams.parnum; j++)
    {
      if (j==i) 
	continue;
#endif
      rxij = rxi - rx[j]; /* distance between two atomes */
      ryij = ryi - ry[j];
      rzij = rzi - rz[j];
      
      rxij = rxij - L * rint(invL * rxij);      /* minimum image */
      ryij = ryij - L * rint(invL * ryij);
      rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
      
      if ( rijSq < rcutSq )/* 'rcut' is the cutoff for V */
	{
  	  srij2   = sigmaSq / rijSq;
	  vij = pow(srij2, ((double)Oparams.N));
	  /*printf("(%d,%d) rijSq: %f vij: %f\n", i, j, rijSq, vij);*/
#ifdef MD_FCONT
	  rij = sqrt(rijSq);
	  *Vlji = *Vlji + vij + dvdr * (rij -rcut);
#else
	  *Vlji = *Vlji + vij;
#endif
#ifdef MD_FCONT
	  vrij = vij / rijSq + dvdr / rij;
#else
	  vrij = vij / rijSq;
#endif
	  Fxi += vrij * rxij;
	  Fyi += vrij * ryij;
	  Fzi += vrij * rzij;
	  if (addnebr==1)
	    nebrTabInteracting[i][n] = 1;
	  ncut++;
	}
      else if (addnebr==1)
	nebrTabInteracting[i][n] = 0;
#if 1
      if (OprogStatus.block && addnebr==2 && (rijSq < rcutSq || nebrTabInteracting[i][n]))
	{
	  rxic = rx[i];
	  ryic = ry[i];
	  rzic = rz[i];
	  energyiF(j, rx[j], ry[j], rz[j], epsilon, sigma, rcut,
		  &VfjOLD);
	  rx[i] = rxi;
	  ry[i] = ryi;
	  rz[i] = rzi;
	  energyiF(j, rx[j], ry[j], rz[j], epsilon, sigma, rcut,
		  &VfjNEW);
	  *Vfi += VfjNEW - VfjOLD; 
	  rx[i] = rxic;
	  ry[i] = ryic;
	  rz[i] = rzic;
	}
#endif
    }
  Fxi *= epsilon4_2N;
  Fyi *= epsilon4_2N;
  Fzi *= epsilon4_2N;
  Fzi -= mg;
#if 0
  /* handle wall reaction though I think it's too rough
   * this solution...check */
  if (rz[i] + L*0.5 < 1E-15 && Fzi < 0.0)
    Fzi = 0.0;
#endif
  rzwall = -Lz*0.5 - Oparams.sigma*0.5;
  *Vlji  = epsilon4*(*Vlji); 
  *Vlji  -= epsilon4*((double)ncut)*pow(sigmaSq/rcutSq, Oparams.N);
  *Vgi   = mg*rzi;
  /*if (fabs(Fzi + mg)> 1E-15)
    printf("Fz[%d]=%.10g Vlji=%.10g\n", i, Fzi, *Vlji);*/
  rziwSq = Sqr(rzi-rzwall);
#ifdef MD_FCONT
  dvdr = - pow(Sqr(sigmaSq/rcutSq), Oparams.Nw)/rcut;
  rziw = rzi-rzwall;
#endif
  if (rziwSq < rcutSq)
    {
      /* -1.0 is a shift to have a zero potential when 
       * an atom touches the wall */
#ifdef MD_FCONT
      vw = pow(sigmaSq/rziwSq, ((double)Oparams.Nw)); 
      /*printf("rziwSq: %f vw: %f Nw:%d\n", rziwSq, vw, Oparams.Nw);*/
      Fzi += 2.0*((double)Oparams.Nw)*(vw/rziwSq+dvdr/rziw)*(rzi-rzwall);
      *Vw = vw - pow(sigmaSq/rcutSq, Oparams.Nw) + dvdr*(rziw-rcut);
#else
      vw = pow(sigmaSq/rziwSq, ((double)Oparams.Nw)); 
      /*printf("rziwSq: %f vw: %f Nw:%d\n", rziwSq, vw, Oparams.Nw);*/
      Fzi += 2.0*((double)Oparams.Nw)*(vw/rziwSq)*(rzi-rzwall);
      *Vw = vw - pow(sigmaSq/rcutSq, Oparams.Nw);
#endif
    }	
  *Vfi += Oparams.alpha*(Sqr(Fxi)+Sqr(Fyi)+Sqr(Fzi));
  if (!OprogStatus.block)
    *Vfi=0.0;
  Fx[i] = Fxi;
  Fy[i] = Fyi;
  Fz[i] = Fzi;
  /*printf("DVlij: %f DVg: %f DVf: %f DVw:%f\n", *Vlji, *Vgi, *Vfi, *Vw );*/
} 

void sumup(int Nm, COORD_TYPE epsilon, 
	   COORD_TYPE sigma, COORD_TYPE rcut)

{
  int i, j, ncut;
  double rxi, ryi, rzi, srij2, rijSq, rcutSq, sigmaSq;
  double vij, rxij, ryij, rzij, epsilon4, L, invL;
  double vw, rziwSq, rzwall, Lz2, epsilon4_2N, mg, Fxi, Fyi, Fzi, fxij, fyij, fzij, vrij, VwCut;
#ifdef MD_FCONT
  double dvdr, rij, rziw; 
#endif
  L = Oparams.L;
  invL = 1.0  / L;
  Vlj = 0.0;
  Vg = 0.0;
  Vf = 0.0;
  Vw = 0.0;
  Lz2 = Oparams.Lz * 0.5;
  epsilon4 = 4.0*epsilon;
  epsilon4_2N = 2.0*epsilon4*((double)Oparams.N);
  mg = Oparams.m * Oparams.g;
  rcutSq = Sqr(Oparams.rcut*Oparams.sigma);
  sigmaSq = Sqr(Oparams.sigma);
  ncut = 0;
#ifdef MD_FCONT
  dvdr = - pow(Sqr(sigmaSq/rcutSq), Oparams.N)/rcut;
#endif

  for (i=0; i < Oparams.parnum; i++)
    {
      Fx[i]=Fy[i]=Fz[i]=0.0;
      rx[i] -= L*rint(rx[i]/L);
      ry[i] -= L*rint(ry[i]/L);
    }
  for (i=0; i < Oparams.parnum-1; i++)
    {
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i]; 
      for (j=i+1; j < Oparams.parnum; j++)
	{
	  Fxi = Fx[i];
	  Fyi = Fy[i];
	  Fzi = Fz[i];
	  rxij = rxi - rx[j]; /* distance between two atomes */
	  ryij = ryi - ry[j];
	  rzij = rzi - rz[j];

	  rxij = rxij - L * rint(invL * rxij);      /* minimum image */
	  ryij = ryij - L * rint(invL * ryij);
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  if (rijSq < Sqr(0.7*Oparams.sigma))
	    {
	      printf("Overlap (%d,%d): %.15e!\n", i, j, rijSq);
	      exit(-1);
	    }
	  if ( rijSq < rcutSq )/* 'rcut' is the cutoff for V */
	    {
	      /*rab   = sqrt(rabSq);*/
	      /*printf("rijSq: %e\n", rijSq);*/
	      srij2   = sigmaSq / rijSq;
	      vij = pow(srij2, Oparams.N);
	      /*printf("vij= %e srij2=%e rijSq: %e rcutSq: %e\n", vij, srij2, rijSq, rcutSq);*/
#ifdef MD_FCONT
	      rij = sqrt(rijSq);
    	      vrij = vij / rijSq + dvdr / rij;
#else
	      vrij = vij / rijSq;
#endif
	      fxij = vrij * rxij;
	      fyij = vrij * ryij;
	      fzij = vrij * rzij;
#ifdef MD_FCONT
	      Vlj = Vlj + vij + dvdr * (rij - rcut);
#else
    	      Vlj = Vlj + vij;
#endif
	  /*printf("Vlj=%e\n", Vlj);*/
	      Fxi   = Fxi + fxij;     /* total force on an atom (a,i)*/
	      Fyi   = Fyi + fyij;
	      Fzi   = Fzi + fzij;

	      Fx[j] = Fx[j] - fxij;  /* -fxab = fxba (3rd law) */
	      Fy[j] = Fy[j] - fyij;
	      Fz[j] = Fz[j] - fzij;
	      ncut++;
	    }
	  Fx[i] = Fxi;
	  Fy[i] = Fyi;
	  Fz[i] = Fzi;
	}	
    }
  Vlj *= epsilon4;
  Vlj -= epsilon4*((double)ncut)*pow(sigmaSq/rcutSq, Oparams.N);
  VwCut = pow(sigmaSq/rcutSq, ((double)Oparams.Nw));
  rzwall = -Lz2-Oparams.sigma*0.5;
  /*printf("epsilon4_2N:%e Nw:%d rzwall:%f\n",epsilon4_2N, Oparams.Nw, rzwall); */
#ifdef MD_FCONT
  dvdr = - pow(Sqr(sigmaSq/rcutSq), Oparams.Nw)/rcut;
#endif
  for (i=0; i < Oparams.parnum; i++)
    {
      Fx[i] *= epsilon4_2N;
      Fy[i] *= epsilon4_2N;
      Fz[i] *= epsilon4_2N;
      Fz[i] -= mg;
      Vg += mg*rz[i];
      /* now we handle the wall through a steep "hardcore-like"
       * repulsive potential */
      rziwSq = Sqr(rz[i]-rzwall);
#ifdef MD_FCONT
      rziw = rz[i]-rzwall;
#endif
      if (rziwSq < rcutSq)
	{
#ifdef MD_FCONT
    	  vw = pow(sigmaSq/rziwSq, ((double)Oparams.Nw)); 
	  /*printf("rziwSq: %f vw: %f Nw:%d\n", rziwSq, vw, Oparams.Nw);*/
	  Fz[i] += 2.0*((double)Oparams.Nw)*(vw/rziwSq+dvdr/rziw)*(rz[i]-rzwall);
	  Vw += vw - VwCut + dvdr * (rziw - rcut);
#else
	  vw= pow(sigmaSq/rziwSq, Oparams.Nw);
	  Fz[i] += 2.0*((double)Oparams.Nw)*(vw/rziwSq)*(rz[i]-rzwall);
	  Vw += vw - VwCut;
#endif
	}	
      Vf += Oparams.alpha*(Sqr(Fx[i])+Sqr(Fy[i])+Sqr(Fz[i]));

    }
    if (!OprogStatus.block)
      Vf = 0.0;
    /*printf("Vlj: %f\n", Vlj); */
}
