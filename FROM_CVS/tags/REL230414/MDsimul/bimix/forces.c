#include<mdsimul.h>
#define SIMUL
/*#define MD_USE_RINT*/
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
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE V, W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx;  
/* used by linked list routines */
extern int *head, *list, *map;  /* arrays of integer */
extern int NCell, mapSize, M;

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, nebrTabLen, nebrTabMax;
/* ================================= */

/* ======================= >>> DSIGN <<< ================================ */
double DSIGN_(double a, double b)
{
  if (b > 0)
    return fabs(a);
  else
    return -fabs(a);
  //return (fabs(b)/b) * fabs(a);
}



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
void  links(COORD_TYPE rcut, COORD_TYPE sigab[NA][NA])
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
  int iCell, i, a, b, Nm;
  
  /*  ZERO HEAD OF CHAIN ARRAY */

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  
  L = cbrt(Vol);
  invL = 1.0 / L;

  loop(iCell, 1, NCell)
    {
      head[iCell] = -1; /* -1 means end of list, that is no more particles */
    }
  pool;

  celli = (COORD_TYPE) M;
  cell  = 1.0 / celli;
  
#ifndef MD_POLYDISP
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
	      printf("M=%d L*cell=%.15G\n", M, L*cell);
	      mdMsg(ALL, NOSYS, NULL, "ERROR", NULL,
		    "Cell size too small for cutoff",
		    NULL);
	      exit(-1);		   
	    }
	}
    }
#endif
  /* Sort all atoms */
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* rxc[][], ryc[][], rzc[] are the corrected coordinates, 
	     that is coordinates restricted to the central box */
	  
	  /* Corrected coordinates, that is restricted to the central box.
	     These coordinates are relative to all the particles actually 
	     in the central box.
	     Finally note that the rx, ry, rz coordinates are uncorrected 
	     because we are also interested in calculating transport 
	     coefficents */
#ifdef    MD_USE_RINT
	  rxc = rx[a][i] - L * rint(invL * rx[a][i]);
	  ryc = ry[a][i] - L * rint(invL * ry[a][i]);
	  rzc = rz[a][i] - L * rint(invL * rz[a][i]);
#else
	  rxc = rx[a][i] - L * RINT_(invL * rx[a][i]);
	  ryc = ry[a][i] - L * RINT_(invL * ry[a][i]);
	  rzc = rz[a][i] - L * RINT_(invL * rz[a][i]);
#endif	  
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
	  //list[NA*i + a] = head[iCell]; 
	  /* Head[iCell] becomes the next
	     of 'NA*i + a' */ 
	  //head[iCell] = NA*i + a;        
	  list[i+Nm*a] = head[iCell]; /* Head[iCell] becomes the next
					   of 'NA*i + a' */ 
	  head[iCell] = i+Nm*a;        
	  /* Now Head[iCell] is 'NA*i + a', in this way we have 
	     add 'NA*i + a' at the Head of the linked list.*/
	}
    }
  pool;
}

/* ====================== >>> BuildNebrListNoLinked <<< ==================== */
void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]) 
{
  int a, b, i, j;
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE rrNebr[NA][NA];
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;
#ifdef MD_POLYDISP
  double L2SQ;
#endif
  int Nm;

  L = cbrt(Vol);
  invL = 1.0  / L;
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
#ifdef MD_POLYDISP
  L2SQ = Sqr(L/2.0);
#endif
  loop(a, 1, NA)
    {
      loop(b, 1, NA)
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
  nebrTabLen = 0;
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++ )
      //loop(i, 1, Nm - 1)
	{
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];
	  for (b = 0; b < NA; b++)
	    {
	      for (j = 0; j < Oparams.parnum[b]; j++)
		//loop(b, 1, NA)  /* b >= a because of symmetry */
		{
		  /* Bisogna assicurare il fatto che ogni coppia
		     appaia una sola volta nella  neighbour list cosi' 
		     essendo a*Nm+i un numero univoco attribuito ad ogni atomo
		     le coppie vengono considerate una sola volta sfruttando
		     l'ordinamento
		  */
		  if ( (b*Nm+j) <= (a*Nm+i) )
		    continue;
		  /* INNER LOOP BEGINS */
		  //loop(j, i + 2, Nm) 
		  /* + 2 because you must remeber that really the all indices 
		     (a, b, i, j) start from 0 ... */
		  /*   i > j because of 3rd law */
		  
		  rxab = rxa - rx[b][j]; /* distance between two atomes */
		  ryab = rya - ry[b][j];
		  rzab = rza - rz[b][j];
#ifdef MD_USE_RINT
		  rxab = rxab - L * rint(rxab * invL);      /* minimum image */
		  ryab = ryab - L * rint(ryab * invL);
		  rzab = rzab - L * rint(rzab * invL);
#else
		  rxab = rxab - L * RINT_(rxab * invL);      /* minimum image */
		  ryab = ryab - L * RINT_(ryab * invL);
		  rzab = rzab - L * RINT_(rzab * invL);
#endif 
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
#ifdef MD_POLYDISP
		  rcutab[a][b] = rCut*(radii[0][i]+radii[0][j]);
		  rrNebr[a][b] = Sqr(rcutab[a][b] + OprogStatus.rNebrShell);
		  if (rrNebr[a][b] > L2SQ)
		    {
		      printf("(rcutoff + rNebrShell)=%f is  too large, it exceeds L/2 = %f\n",
			     sqrt(rrNebr[a][b]), L/2.0);
		      exit(-1);
		    }
		  
#endif
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
		       /*
			 nebrTab[0][nebrTabLen] = NA*i + a;
			 nebrTab[1][nebrTabLen] = NA*j + b;*/
		       nebrTab[0][nebrTabLen] = i + a*Nm;
		       nebrTab[1][nebrTabLen] = j + b*Nm;
		       //printf("(%d-%d)[(%d,%d)-(%d,%d)]\n ",
		       //     i + a*Nm, j + b*Nm, a, i, b, j);
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
void BuildNebrList(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]) 
{
  int Nm, a, b, i, j;
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE rrNebr[NA][NA];
  int  iCell, jCell0, jCell, nabor, hd, lst;
  COORD_TYPE rxa, rya, rza, rabSq, rxab, ryab, rzab;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];

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
	  /*
	    i = hd / NA;
	    a = hd % NA; */
	  a = hd / Nm;
	  i = hd % Nm; 
	  //doFather{printf("cell: %d particle: (%d,%d)\n", iCell, i, a);};
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];
	  
	  /* LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL */
	  lst = list[hd]; 
				
	  while(lst > -1)
	    /* Atoms in the same molecule don't interact */
	    { 
	      /*
		j = lst / NA;
		b = lst % NA;
	      */
	      b = lst / Nm;
	      j = lst % Nm;

	      if ((j == i) && (a == b))
		/* un atomo non puo' interagire con se stesso */
		{
		  lst = list[lst];
		  continue;
		}
	    
	      rxab = rxa - rx[b][j]; /* distance between two atomes */
	      ryab = rya - ry[b][j];
	      rzab = rza - rz[b][j];
	      
	      /* minimum image */
#ifdef        MD_USE_RINT
	      rxab = rxab - L * rint(invL * rxab);      
	      ryab = ryab - L * rint(invL * ryab);
	      rzab = rzab - L * rint(invL * rzab);
#else
	      rxab = rxab - L * RINT_(invL * rxab);      
	      ryab = ryab - L * RINT_(invL * ryab);
	      rzab = rzab - L * RINT_(invL * rzab);
#endif
	      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
#ifdef MD_POLYDISP
	      rcutab[a][b] = rCut*(radii[0][i]+radii[0][j]);
	      rrNebr[a][b] = Sqr(rcutab[a][b] + OprogStatus.rNebrShell);
#endif
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
		  /*
		    nebrTab[0][nebrTabLen] = NA*i + a;
		    nebrTab[1][nebrTabLen] = NA*j + b;*/
		  nebrTab[0][nebrTabLen] = i + a*Nm;
		  nebrTab[1][nebrTabLen] = j + b*Nm;
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
		  /*
		    j = lst / NA;
		    b = lst % NA;
		  */
		  b = lst / Nm;
		  j = lst % Nm;
		  if (  ( j == i ) && ( a == b )  ) 
		    /* Un atomo non puo' interagire con se stesso
		       !!!!!!!!!!!!!!!!!!!!!!!! */
		    {
		      lst = list[lst];
		      continue;
		    }
		 
		  rxab = rxa - rx[b][j]; /* distance between two atomes */
		  ryab = rya - ry[b][j];
		  rzab = rza - rz[b][j];
#ifdef  MD_USE_RINT 
		  rxab = rxab - L * rint(invL * rxab);    /* minimum image */
		  ryab = ryab - L * rint(invL * ryab);
		  rzab = rzab - L * rint(invL * rzab);
#else
		  rxab = rxab - L * RINT_(invL * rxab);    /* minimum image */
		  ryab = ryab - L * RINT_(invL * ryab);
		  rzab = rzab - L * RINT_(invL * rzab);

#endif
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
#ifdef MD_POLYDISP
    		  rcutab[a][b] = rCut*(radii[0][i]+radii[0][j]);
    		  rrNebr[a][b] = Sqr(rcutab[a][b] + OprogStatus.rNebrShell);
#endif
	
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
		      /*nebrTab[0][nebrTabLen] = NA*i + a;
			nebrTab[1][nebrTabLen] = NA*j + b;
		      */
		      nebrTab[0][nebrTabLen] = i + a*Nm;
		      nebrTab[1][nebrTabLen] = j + b*Nm;
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


/* =========================== >>> kinet <<< ============================== */
void kinet(COORD_TYPE** velx, COORD_TYPE** vely, COORD_TYPE** velz,
	   COORD_TYPE VOL1)
{
  /* DESCRIPTION:
     Calculate the kinetic energy of the sample, this is done for 
     correcting the Hoover friction coefficent with the appropriate value
     of the temperature (i.e. the predicted one) */
  int i, a;
  COORD_TYPE dlnV, px, py, pz;

  /* NOTE: K is not shared so this loop is done by father and child */
  K = 0.0;
  L = cbrt(Vol);
  invL = 1.0 / L;

  dlnV = VOL1 / Vol;
  dlnV /= 3.0;
 
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++) 
    	{
	  px = velx[a][i] - dlnV * rx[a][i];
	  py = vely[a][i] - dlnV * ry[a][i];
	  pz = velz[a][i] - dlnV * rz[a][i];
	  K = K + Oparams.m[a] * (Sqr(px) + Sqr(py) + Sqr(pz));
	}
    }

  /* Kp is the kinetic energy calculated using 'predicted' velocities at time
     t (v(t)) */
  K *= 0.5;
}

void checkNebrRebuild(void)
{
  int i, a;
  COORD_TYPE vv, vvMax = 0.0;

  for (a = 0; a < NA; a++)  
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  vv = Sqr(vx[a][i]) + Sqr(vy[a][i]) + Sqr(vz[a][i]);
	  if (vv > vvMax) 
	    vvMax = vv;
	}
    }
  dispHi = dispHi + sqrt(vvMax) * Oparams.steplength;
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;

}

#undef BASINS
#ifdef BASINS
/*extern double minimumImageFB(double c);*/
#endif

/* ============================ >>> force <<< ==============================*/
void LJForce(COORD_TYPE epsab[NA][NA], 
	     COORD_TYPE sigab[NA][NA], COORD_TYPE rcut)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  int a, b, i, j;
  int ncut[NA][NA];
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
#ifdef SOFT_SPHERE
#ifdef MD_POLYDISP
  double PP, rab, fcut;
#else
  int PP, rab=1.0, fcut=0.0;
#endif
  double Epsab4;
#elif defined(NM_SPHERE)
  double vabNN, vabMM, factor, sigmaFactorSq[NA][NA];
#endif
#ifdef MD_GRAVITY
  const double g = 9.81;/**(1E10/(0.5E12*0.5E12));*//* g in unità ridotte!*/
  double mg;
#endif
  COORD_TYPE sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA], Epsab24;
  COORD_TYPE rxab, ryab, rzab, rabSq, fxab, fyab, fzab;
  COORD_TYPE srab2, srab4, srab8, srab6, srab12, fab, vab, wab;
  COORD_TYPE Fxa, Fya, Fza, rxa, rya, rza;
  COORD_TYPE vabCut, Vcab[NA][NA];
  COORD_TYPE Vab[NA][NA], Wab[NA][NA], dvdr[NA][NA];
  /* Local variables to implement linked list */
  int  Nm, n, nebrTab0, nebrTab1, invLH;


  /* ======================================================================= */
  /*calculate useful quantities
   NOTE: We refer to ab-quantities as all bidimensional arrays 
         that depend upon atoms pair(e.g. 'Vab[a][b]' is a 'ab'-variable, 
	 instead 'sigab[a][b]' is an ab-constant */
#ifdef NM_SPHERE
  factor = pow(2.0,1.0/6.0);
#endif

  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a because of symmetry */
	{
#ifdef MD_POLYDISP
	  Vcab[a][b] = 0.0;
#endif
	  /* useful ab-constants inside OUTER LOOP below */
	  rcutab[a][b] = rcut * sigab[a][b];
	  rcutabSq[a][b] = Sqr(rcutab[a][b]);
	  sigabSq[a][b] = Sqr(sigab[a][b]);
#ifdef NM_SPHERE
	  epsab4[a][b] = epsab[a][b] / (((double) Oparams.NN) - ((double)Oparams.MM));
	  epsab24[a][b] = 6.0 * epsab[a][b] / (((double) Oparams.NN) - ((double)Oparams.MM));
	  sigmaFactorSq[a][b] = Sqr(sigab[a][b]*factor);
#else
	  epsab4[a][b] = 4.0 * epsab[a][b];
	  epsab24[a][b] = 24.0 * epsab[a][b];
#endif
	  srab2   = sigabSq[a][b] / rcutabSq[a][b];
	  srab6   = srab2 * srab2 * srab2;
	  srab12  = Sqr(srab6);
#ifdef SOFT_SPHERE
	  /* metterci il valore esatto eventualmente se si vuole la continuita'
	     nella forza al cutoff */
	  dvdr[a][b] = 0.0;
#elif defined(NM_SPHERE)
	  dvdr[a][b] = 0.0;
#else	  
	  dvdr[a][b] = epsab24[a][b] * (srab6 - 2.0 * srab12) / rcutab[a][b];
#endif	  
	  /* initialize ab-variables */
	  ncut[a][b] = 0;
	  Vab[a][b] = 0.0;
	  Wab[a][b] = 0.0;
	}
    }
#ifdef SOFT_SPHERE
  PP = Oparams.PP;
#endif
  L = cbrt(Vol);
  invL = 1.0  / L;
  invLH = 2.0 * invL;

  V = 0.0; /* potential energy */
  W = 0.0; /* virial function */
  Wxy = 0.0; /* virial off-diagonal terms of pressure tensor */
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  /* Loop over cells */
  Vc = 0.0;
  
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  Fx[a][i] = 0.0;
	  Fy[a][i] = 0.0;
	  Fz[a][i] = 0.0;
	}
    }

  for(n = 0; n < nebrTabLen; n++)
    { 
      nebrTab0 = nebrTab[0][n]; 
      nebrTab1 = nebrTab[1][n];

      /*
	i = nebrTab0 / NA;
	a = nebrTab0 % NA;
	j = nebrTab1 / NA;
	b = nebrTab1 % NA;
      */
#if 0
      a = nebrTab0 / Nm;
      i = nebrTab0 % Nm;
      b = nebrTab1 / Nm;
      j = nebrTab1 % Nm;
#else
      a = (nebrTab0 < Nm)?0:1;
      b = (nebrTab1 < Nm)?0:1;
      i = (a==1)?(nebrTab0 - Nm):nebrTab0;
      j = (b==1)?(nebrTab1 - Nm):nebrTab1; 
#endif 
#if 0
      rxa = rx[a][i];
      rya = ry[a][i];
      rza = rz[a][i];
	  
      Fxa = Fx[a][i];
      Fya = Fy[a][i];
      Fza = Fz[a][i];
#ifdef SOFT_SPHERE 
      Epsab4 = epsab4[a][b];
#else
      Epsab24 = epsab24[a][b];
#endif
#endif
      rxab = rx[a][i] - rx[b][j]; /* distance between two atomes */
      ryab = ry[a][i] - ry[b][j];
      rzab = rz[a][i] - rz[b][j];

#ifdef BASINS      
      /*rxab = minimumImageFB(rxab);
	ryab = minimumImageFB(ryab);
	rzab = minimumImageFB(rzab);*/
      rxab -= L*( (double) ((int) (rxab*invLH)) );
      ryab -= L*( (double) ((int) (ryab*invLH)) );
      rzab -= L*( (double) ((int) (rzab*invLH)) );
      /* quando si minimizza tutte le particelle sono
	 nella prima scatola quindi si può usare questa */
#else
#ifdef MD_USE_RINT 
      rxab = rxab - L * rint(invL * rxab); 
      ryab = ryab - L * rint(invL * ryab);
#if !defined(MD_GRAVITY)
      rzab = rzab - L * rint(invL * rzab);
#endif
#else
      rxab = rxab - L * RINT_(invL * rxab); 
      ryab = ryab - L * RINT_(invL * ryab);
#if !defined(MD_GRAVITY)
      rzab = rzab - L * RINT_(invL * rzab);
#endif
#endif      
#endif      
      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
#ifdef MD_POLYDISP
      sigabSq[a][b] = Sqr(radii[0][i]+radii[0][j]);
      rcutab[a][b] = rcut*(radii[0][i]+radii[0][j]);
      rcutabSq[a][b] = Sqr(rcutab[a][b]);
#endif
      if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
	{
#ifdef SOFT_SPHERE
	  srab2 = sigabSq[a][b] /rabSq;
#if   defined(MD_STATIC_PP36)
          srab4=srab2*srab2;
          srab8=srab4*srab4;
          vab = srab8*srab8*srab2;
	  vab = vab*vab;
          wab = 36.0*vab;
#elif   defined(MD_STATIC_PP18)
          srab4=srab2*srab2;
          srab8=srab4*srab4;
          vab = srab8*srab8*srab2;
          wab = 18.0*vab;
#elif   defined(MD_STATIC_PP12)
	  vab = srab2*srab2*srab2;
	  vab = Sqr(vab);
	  wab = 12.0*vab;
#elif defined(MD_STATIC_PP8)
	  vab = srab2*srab2;
	  vab = Sqr(vab);
	  wab = 8.0*vab;
#elif defined(MD_STATIC_PP6)
	  vab = srab2*srab2*srab2;
	  wab = 6.0*vab;   
#else
  	  vab = pow(srab2, ((double)PP)/2.0);
	  wab = ((double)PP)*vab;
#endif
#ifdef MD_POLYDISP
	  /* N.B. small PP means long range and with small cutoff
	     used here implies integrator instabilities, that is why
	     we emply here potential and force shift to fix this
	     (see Allen-Tildesley pag. 146) */
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  vabCut = pow(srab2,PP/2.0);
	  fcut = ((double)PP)*vabCut/rcutab[a][b];
	  rab = sqrt(rabSq);
#endif
	  Vab[a][b] += vab;
	  Wab[a][b] += wab;
	  fab   = epsab4[a][b] * (wab / rabSq - fcut / rab);
	  /* force between two atoms */
	  fxab  = fab * rxab;         
	  fyab  = fab * ryab;
	  fzab  = fab * rzab;
#ifdef MD_POLYDISP
	  Vcab[a][b] += vab - vabCut + (rab-rcutab[a][b])*fcut;
#endif
#elif defined(NM_SPHERE)
#ifdef MD_POLYDISP
	  sigmaFactorSq[a][b] = sigabSq[a][b]*Sqr(factor);
#endif
	  srab2 = sigmaFactorSq[a][b] / rabSq;
	  vabNN = ((double)Oparams.MM)*pow(srab2, ((double)Oparams.NN)/2.0);
	  vabMM = -((double)Oparams.NN)*pow(srab2, ((double)Oparams.MM)/2.0);
	  vab = vabNN + vabMM;
  	  wab = ((double)Oparams.NN)*vabNN + ((double)Oparams.MM)*vabMM ;
	  Vab[a][b] += vab;
	  Wab[a][b] += wab;
	  fab   = epsab4[a][b] * wab / rabSq;
	  /* force between two atoms */
	  fxab  = fab * rxab;         
	  fyab  = fab * ryab;
	  fzab  = fab * rzab;
#ifdef MD_POLYDISP
    	  srab2 = sigmaFactorSq[a][b] / rcutabSq[a][b];
       	  vabCut = ((double)Oparams.MM)*pow(srab2, Oparams.NN/2.0)-((double)Oparams.NN)*pow(srab2,Oparams.MM/2.0);
	  Vcab[a][b] += vab - vabCut;
#endif
/*
        # elif defined(SOFT_SPHERE)
	  srab2 = sigabSq[a][b] /rabSq;
	  vabNN = pow(srab2, ((double)Oparams.NN)/2.0);
	  vabMM = -pow(srab2, ((double)Oparams.MM)/2.0);
	  vab = vabNN + vabMM;
	  wab = ((double)Oparams.NN)*vabNN - ((double)Oparams.MM)*vabMM ;
	  Vab[a][b] += vab;
	  Wab[a][b] += wab;
	  fab   = epsab4[a][b] * wab / rabSq;
	  // force between two atoms 
	  fxab  = fab * rxab;         
	  fyab  = fab * ryab;
	  fzab  = fab * rzab;
	  */
#else
	  srab2   = sigabSq[a][b] / rabSq;
	  srab6   = srab2 * srab2 * srab2;
	  srab12  = Sqr(srab6);
	  
	  vab     = srab12 - srab6;
	  /*vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);*/
	  wab     = vab + srab12;
	  
	  Vab[a][b] = Vab[a][b] + vab;
	  
	  /* total potential between all a-b atoms pairs */
	  Wab[a][b]   = Wab[a][b] + wab; 
	  /* NOTE: If you will use a shifted-force potential then 
	     calculate the force using that potential */

	  fab   = epsab24[a][b] * wab / rabSq;
	  /* force between two atoms */
	  fxab  = fab * rxab;         
	  fyab  = fab * ryab;
	  fzab  = fab * rzab;

#ifdef MD_POLYDISP
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;
	  Vcab[a][b] += vab - vabCut;
#endif
#endif

	  /* Virial off-diagonal terms of atomic pressure tensor */
#if 0
	  Wxy += rxab * fyab;
	  Wyz += ryab * fzab;
	  Wzx += rzab * fxab;
	  Wxx += rxab * fxab;
	  Wyy += ryab * fyab;
	  Wzz += rzab * fzab;
	  
#endif
	  Fx[a][i]   = Fx[a][i] + fxab;     /* total force on an atom (a,i)*/
	  Fy[a][i]   = Fy[a][i] + fyab;
	  Fz[a][i]   = Fz[a][i] + fzab;
	  
	  Fx[b][j] = Fx[b][j] - fxab;  /* -fxab = fxba (3rd law) */
	  Fy[b][j] = Fy[b][j] - fyab;
	  Fz[b][j] = Fz[b][j] - fzab;
	  ++ncut[a][b];
	}
#if 0 
      Fx[a][i] = Fxa;
      Fy[a][i] = Fya;
      Fz[a][i] = Fza;
#endif
    }

  /* OUTER LOOP ENDS */

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
#ifdef SOFT_SPHERE
#ifndef MD_POLYDISP
  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a */
	{
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  vabCut = pow(srab2,PP/2.0);
	  Vcab[a][b] = Vab[a][b] - ((double)ncut[a][b]) * vabCut;
	  /* ncut[a][b] is the number of atoms pairs a-b within 
	     rcutab[a][b] */ 
	}
    }
#endif
  /* MULTIPLY FOR ENERGY FACTORS */
  loop(a, 1, NA)
    {
      loop(b, 1, NA) 
	{
	  V  = V + Vab[a][b]  * epsab4[a][b];
	  Vc = Vc + Vcab[a][b] * epsab4[a][b];
	  W  = W + Wab[a][b] * epsab4[a][b] / 3.0;
	}
    }
#elif defined(NM_SPHERE)
#ifndef MD_POLYDISP
 for (a = 0; a < NA; a++)
   for (b = 0; b < NA; b++)
     {
       srab2 = sigmaFactorSq[a][b] / rcutabSq[a][b];
       vabCut = ((double)Oparams.MM)*pow(srab2, Oparams.NN/2.0)-((double)Oparams.NN)*pow(srab2,Oparams.MM/2.0);
       Vcab[a][b] = Vab[a][b] - ((double)ncut[a][b]) * vabCut;
     }
#endif
 /* MULTIPLY FOR ENERGY FACTORS */
  loop(a, 1, NA)
    {
      loop(b, 1, NA) 
	{
	  V  = V + Vab[a][b]  * epsab4[a][b];
	  Vc = Vc + Vcab[a][b] * epsab4[a][b];
	  W  = W + Wab[a][b]  * epsab24[a][b] / 3.0;
	}
    }

#else
#ifndef MD_POLYDISP
  loop(a, 1, NA)
    {
      loop(b, 1, NA) /* b >= a */
	{
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;
	  Vcab[a][b] = Vab[a][b] - ((double)ncut[a][b]) * vabCut;
	  /* ncut[a][b] is the number of atoms pairs a-b within 
	     rcutab[a][b] */ 
	}
    }
#endif
  /* MULTIPLY FOR ENERGY FACTORS */
  loop(a, 1, NA)
    {
      loop(b, 1, NA) 
	{
	  V  = V + Vab[a][b]  * epsab4[a][b];
	  Vc = Vc + Vcab[a][b] * epsab4[a][b];
	  W  = W + Wab[a][b]  * epsab24[a][b] / 3.0;
	}
    }
#endif
#ifdef MD_GRAVITY
  for(a = 0; a < NA; a++)
    {
      mg = Oparams.m[a]*g;
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  Fz[a][i] -= mg; 
	  V  += mg*rz[a][i];
	  Vc += mg*rz[a][i];
	  /* e il viriale? */
	}
    }
#endif
}  
/* ============================ >>> force <<< ==============================*/
void LJForceNONL(COORD_TYPE epsab[NA][NA], 
	     COORD_TYPE sigab[NA][NA], COORD_TYPE rcut)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  int a, b, i, j;
  int ncut[NA][NA];
  COORD_TYPE rcutab[NA][NA], rcutabSq[NA][NA]; 
  COORD_TYPE sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA], Epsab24;
  COORD_TYPE rxab, ryab, rzab, rabSq, fxab, fyab, fzab;
  COORD_TYPE srab2, srab6, srab12, fab, vab, wab;
  COORD_TYPE Fxa, Fya, Fza, rxa, rya, rza;
  COORD_TYPE vabCut[NA][NA], Vcab[NA][NA];
  COORD_TYPE dvdr[NA][NA];
  /* Local variables to implement linked list */
  int  Nm, n, nebrTab0, nebrTab1;

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
	  vabCut[a][b] = epsab4[a][b]*(srab12 - srab6);
	  
	  /* initialize ab-variables */
	  ncut[a][b] = 0;
	  
	  
	}
    }

  L = cbrt(Vol);
  invL = 1.0  / L;
  V = 0.0; /* potential energy */
  W = 0.0; /* virial function */
  Wxy = 0.0; /* virial off-diagonal terms of pressure tensor */
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  /* Loop over cells */
  Vc = 0.0;

  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  Fx[a][i] = 0.0;
	  Fy[a][i] = 0.0;
	  Fz[a][i] = 0.0;
	}
    }

  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++ )
	//loop(i, 1, Nm - 1)
	{
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];
	  for (b = 0; b < NA; b++)
	    {
	      for (j = 0; j < Oparams.parnum[b]; j++)
		//loop(b, 1, NA)  /* b >= a because of symmetry */
		{
		  /* Bisogna assicurare il fatto che ogni coppia
		     appaia una sola volta nella  neighbour list cosi' 
		     essendo a*Nm+i un numero univoco attribuito ad ogni atomo
		     le coppie vengono considerate una sola volta sfruttando
		     l'ordinamento
		  */
		  if ( (b*Nm+j) <= (a*Nm+i) )
		    continue;

		  Fxa = Fx[a][i];
		  Fya = Fy[a][i];
		  Fza = Fz[a][i];
		  
		  Epsab24 = epsab24[a][b];
		  rxab = rxa - rx[b][j]; /* distance between two atomes */
		  ryab = rya - ry[b][j];
		  rzab = rza - rz[b][j];
#ifdef MD_USE_RINT 
		  rxab = rxab - L * rint(invL * rxab);   
		  ryab = ryab - L * rint(invL * ryab);
	      	  rzab = rzab - L * rint(invL * rzab);
#else
		  rxab = rxab - L * RINT_(invL * rxab);   
		  ryab = ryab - L * RINT_(invL * ryab);
	      	  rzab = rzab - L * RINT_(invL * rzab);

#endif
		  rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		  
		  if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		    {
		      srab2   = sigabSq[a][b] / rabSq;
		      srab6   = srab2 * srab2 * srab2;
		      srab12  = Sqr(srab6);
		      
		      vab     = srab12 - srab6;
		      //vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);
		      wab     = vab + srab12;
		      
		      V = V + epsab4[a][b]*vab;
		      Vc = Vc + epsab4[a][b]*vab - vabCut[a][b];
		      		      
		      /* NOTE: If you will use a shifted-force potential then 
			 calculate the force using that potential */
		      fab   = wab / rabSq;
		      /* force between two atoms */
		      fxab  = fab * rxab * Epsab24;         
		      fyab  = fab * ryab * Epsab24;
		      fzab  = fab * rzab * Epsab24;
		      
		      /* Virial off-diagonal terms of atomic pressure tensor */
#if 0
		      Wxy += rxab * fyab;
		      Wyz += ryab * fzab;
		      Wzx += rzab * fxab;
		      Wxx += rxab * fxab;
		      Wyy += ryab * fyab;
		      Wzz += rzab * fzab;
		      
#endif
		      Fxa   = Fxa + fxab;     /* total force on an atom (a,i)*/
		      Fya   = Fya + fyab;
		      Fza   = Fza + fzab;
		      
		      Fx[b][j] = Fx[b][j] - fxab;  /* -fxab = fxba (3rd law) */
		      Fy[b][j] = Fy[b][j] - fyab;
		      Fz[b][j] = Fz[b][j] - fzab;
		      ++ncut[a][b];
		    }
		  
		  Fx[a][i] = Fxa;
		  Fy[a][i] = Fya;
		  Fz[a][i] = Fza;
		}
	      
	      /* OUTER LOOP ENDS */
	    }
	}
    }
}  
