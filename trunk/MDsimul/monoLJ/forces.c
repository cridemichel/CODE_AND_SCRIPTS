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

extern char TXTA[10][MSG_LEN];
extern char TXT[MSG_LEN];

// Potential anti-crystalline energy
C_T VAC, WAC, VLJ, WLJ;
extern C_T *kcx, *kcy, *kcz;
extern int NNC;

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
  
  L = cbrt(Vol);
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
void BuildNebrListNoLinked(int Nm, COORD_TYPE rCut, COORD_TYPE sigma) 
{
  int i, j;
  COORD_TYPE rcut, rcutSq; 
  COORD_TYPE rrNebr;
  COORD_TYPE rxi, ryi, rzi, rijSq, rxij, ryij, rzij;

  ProcSync0();

  L = cbrt(Vol);
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
  //printf("sqrt(rrNebr[%d][%d]):%f\n", a, b, sqrt(rrNebr[a][b]));

  //doFather{printf("STEP: %d rebuilding nbr Lst\n",Oparams.curStep);};
  nebrTabLen = 0;
  loop(i, 1, Nm - 1)
    {
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i];
      /* INNER LOOP BEGINS */
      loop(j, i + 2, Nm) 
	/* + 2 because you must remeber that really the all indices 
	   (a, b, i, j) start from 0 ... */
	/*   i > j because of 3rd law */
	{
	  rxij = rxi - rx[j]; /* distance between two atomes */
	  ryij = ryi - ry[j];
	  rzij = rzi - rz[j];
	  rxij = rxij - L * rint(rxij * invL);      /* minimum image */
	  ryij = ryij - L * rint(ryij * invL);
	  rzij = rzij - L * rint(rzij * invL);
	  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	  
	  if ( rijSq < rrNebr )/* 'rcut' is the cutoff for V */
	    {
	      if (nebrTabLen >= nebrTabMax)
		{
		  sprintf(TXTA[0], "nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			 nebrTabLen);
		  sprintf(TXTA[1], "particles: (%d)-(%d)\n",i, j);
		  sprintf(TXTA[2], "ERROR: Neighbourlist overflow!\n");
		  mdPrintf(STD, TXTA[0], TXTA[1], TXTA[2], NULL);
		  exit(-1);
		}
	      nebrTab[0][nebrTabLen] = i;
	      nebrTab[1][nebrTabLen] = j;
	      //doFather{printf("First Inter iCell: %d: (%d,%d)-(%d,%d)\n", 
				//       iCell,i,a,j,b);};
	      ++nebrTabLen; /* Increment table element counter */
	      
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
  //printf("sqrt(rrNebr[%d][%d]):%f\n", a, b, sqrt(rrNebr[a][b]));

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
	  i = hd;
	  //doFather{printf("cell: %d particle: (%d,%d)\n", iCell, i, a);};
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
	      rzij = rzij - L * rint(invL * rzij);
	      
	      rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	     
	      if ( rijSq < rrNebr)/* 'rcut' is the cutoff for V */
		{
		  if (nebrTabLen >= nebrTabMax)
		    {
		      sprintf(TXTA[0], "nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax, 
			     nebrTabLen);
		      sprintf(TXTA[1], "particles: (%d)-(%d)\n", i, j);
		      sprintf(TXTA[2], "ERROR: Neighbourlist overflow!\n");
		      mdPrintf(STD, TXTA[0], TXTA[1], TXTA[2], NULL);
		      exit(-1);
		    }
		  nebrTab[0][nebrTabLen] = i;
		  nebrTab[1][nebrTabLen] = j;
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
		  rzij = rzij - L * rint(invL * rzij);
		  rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
		  
		  if ( rijSq < rrNebr)/* 'rcut' is the cutoff for V */
		    {
		      if (nebrTabLen >= nebrTabMax)
			{
			  sprintf(TXTA[0], "nebrTabMax: %d nebrTabLen: %d\n", nebrTabMax,
				 nebrTabLen);
			  sprintf(TXTA[1], "ERROR: Neighbourlist overflow!\n");
			  sprintf(TXTA[2], "particles: (%d)-(%d)\n",i,j);
			  mdPrintf(STD, TXTA[0], TXTA[1], TXTA[2], NULL);
			  exit(-1);
			}
		      //doFather{fprintf(stderr, "jCell0: %d nabor:%d Second Inter cell:%d (%d,%d)-(%d,%d)\n",
		      //      jCell0, nabor, jCell, i,a,j,b);};
		      nebrTab[0][nebrTabLen] = i;
		      nebrTab[1][nebrTabLen] = j;
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
void kinet(int Nm, COORD_TYPE* velx, COORD_TYPE* vely, COORD_TYPE* velz,
	   COORD_TYPE VOL1)
{
  /* DESCRIPTION:
     Calculate the kinetic energy of the sample, this is done for 
     correcting the Hoover friction coefficent with the appropriate value
     of the temperature (i.e. the predicted one) */
  int i;
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
      px = velx[i] - dlnV * rx[i];
      py = vely[i] - dlnV * ry[i];
      pz = velz[i] - dlnV * rz[i];
      K = K + Oparams.m * (Sqr(px) + Sqr(py) + Sqr(pz));
    }

  /* Kp is the kinetic energy calculated using 'predicted' velocities at time
     t (v(t)) */
  K *= 0.5;
}

void checkNebrRebuild(void)
{
  int i, Nm = Oparams.parnum;
  COORD_TYPE vv, vvMax = 0.0;

  ProcSync0();
  loop(i, 1, Nm)
    {
      vv = Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]);
      if (vv > vvMax) 
	vvMax = vv;
    }
 
 dispHi = dispHi + sqrt(vvMax) * Oparams.steplength;
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
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

/* =========================== >>> buildMesh <<< ========================= */
void buildMesh(COORD_TYPE kmax, COORD_TYPE Dk)
{
  /* Forza anti cristallina */
  int ikx, iky, ikz;
  int KK, ind;
  COORD_TYPE kpDkSq, kmDkSq, kSq;
  COORD_TYPE fact;
  /* ==== Ricerca terne k da controllare per potenziale anti-cristallino==== */

  L = cbrt(Vol);
  
  fact = L/(2.0*pi);
  kpDkSq = Sqr((kmax+Dk)*fact);
  kmDkSq = Sqr((kmax-Dk)*fact);
  
  KK = ((int)sqrt(kpDkSq)) + 1;
 
  ind=0;
  //for (ikx = -KK; ikx <= KK; ikx++)
  for (ikx = -KK; ikx <= KK; ikx++)
    for(iky = -KK; iky <= KK; iky++)
      for(ikz = -KK; ikz <= KK; ikz++)
	{   
	  kSq = Sqr(ikx)+Sqr(iky)+Sqr(ikz);
	  if((kSq < kpDkSq) &&  (kSq > kmDkSq))  
	    {
	      if(ind >= NK)
		{
		  mdPrintf(ALL, "Troppi Q da controllare !!!\n");
		  exit(-1);
		}
	      //kkx[ind]=ikx;
	      //kky[ind]=iky;
	      //kkz[ind]=ikz;
	      kcx[ind]=ikx;
	      kcy[ind]=iky;
	      kcz[ind]=ikz;
	      ind++;
	    }
	}
  //NNQ=ind;
  NNC=ind;
  //printf("NK: %d NNC: %d\n",NK,  NNC);
}

/* ========================== >>> rdiq <<< ============================= */
void rdiq(C_T *rx, C_T *ry, C_T *rz, C_T kx, C_T ky, C_T kz, C_T* Re_rh, C_T* Im_rh)
{
     
  C_T arg;//, sqrtNm;
  int Nm, i;
  
  *Re_rh = 0.0;
  *Im_rh = 0.0;
  Nm = Oparams.parnum;
  //sqrtNm = sqrt((C_T)Nm);

  for(i = 0; i< Nm; i++)
    {
      arg = rx[i]*kx + ry[i]*ky + rz[i]*kz;
      *Re_rh = *Re_rh + cos(arg);
      *Im_rh = *Im_rh + sin(arg);
    }
  
  *Re_rh /= (sqrt((C_T)Nm));
  *Im_rh /= (sqrt((C_T)Nm));

}

/* ======================== >>> sumContribs <<< ========================*/
void sumFRContribs(void)
{
  V = VAC + VLJ;
  W = WAC + WLJ;
  //printf("WAC: %f WLJ:%f\n", WAC, WLJ);
  Vc += VAC;
}

/* ========================== >>> ACForce <<< =========================== */
void ACForce(int Nm, COORD_TYPE kmax, COORD_TYPE Dk, C_T alpha, C_T S0)
{
  int ind, indMax, iq, i;
  COORD_TYPE wxx, wyy, wzz, kx, ky, kz, pi2, ss, Smax, cost1, invNm, cost2;
  C_T Re_rh, Im_rh, Im_ex, dd, gg, F, Fxx, Fyy, Fzz, arg, sm, invL;
  
  pi2 = 2.0*pi;
  sm = 0.0;
  invL = 1.0 / cbrt(Vol);
  cost2 = pi2*invL;

  invNm = 1.0 / ((C_T)Nm);
  //cost1 = -8.0 * alpha / sqrt((C_T)Nm); 
  // N.B. 
  // The AC potential is 1/2*Sum[(S(k)-S0)^2]
  cost1 = -2.0 * alpha / sqrt((C_T)Nm); 
  VAC = 0.0;
  WAC = 0.0;

  //printf("antiCry: %d NNC: %d\n", OprogStatus.AntiCry, NNC);
  for(ind = 0; ind < NNC; ind++)
    {
      kx = cost2 * kcx[ind];
      ky = cost2 * kcy[ind];
      kz = cost2 * kcz[ind];
      rdiq(rx, ry, rz, kx, ky, kz, &Re_rh, &Im_rh);
      ss = Sqr(Re_rh) + Sqr(Im_rh);
      if (ss > sm) 
	{
	  Smax = ss; /* max value reached by S(k) */
	  indMax = ind;/* k for which S(k) reach its maximum */
	}

      dd = ss - S0;
      if (dd > 0.0)
	{
	  //printf("STEP: %d ind: %d dd = %f ss:%f S0: %f|", Oparams.curStep, ind, dd, ss, S0);

	  gg =  dd * cost1;
	  for(i = 0; i < Nm; i++)
	   {
	     arg = kx*rx[i] + ky*ry[i] + kz*rz[i];
	     //ex = rh * cexp(-Im*arg);
	     //Re_ex = Re_rh * cos(-arg) - Im_rh * sin(-arg);
	     Im_ex = Re_rh * sin(-arg) + Im_rh * cos(-arg);
	     F = gg * Im_ex;
	     Fxx = F * kx;
	     Fyy = F * ky;
	     Fzz = F * kz;
	     //printf("Fxx: %f, Fyy: %f Fzz:%f\n ", Fxx, Fyy, Fzz);
	     Fx[i] += Fxx;
	     Fy[i] += Fyy;
	     Fz[i] += Fzz;
	     // Calcolca il contributo al viriale
	     // ma questa espressine non dipende dalla cella!!!!!!!
	     // boh !!!!!!!!!!!!!!!!
	     Wxy = rx[i] * Fyy;
	     Wyz = ry[i] * Fzz;
	     Wzx = rz[i] * Fxx;
	     wxx = rx[i] * Fxx;
	     wyy = ry[i] * Fyy;
	     wzz = rz[i] * Fzz;
	     Wxx += wxx;
	     Wyy += wyy;
	     Wzz += wzz;
	     WAC += wxx + wyy + wzz;//rx[i]*Fxx + ry[i]*Fyy + rz[i]*Fzz;
	   }
	  // Contribute to potential energy
	  VAC = VAC + Sqr(dd);
	}
    }
  VAC *= 0.5*alpha;
  //VAC *= alpha;
}

#ifdef MPI
extern int my_rank;
#endif

/* ============================ >>> force <<< ==============================*/
void LJForce(int Nm, COORD_TYPE epsilon, 
	     COORD_TYPE sigma, COORD_TYPE rcut)
{
  /* ======================== >>> LOCAL VARIABLES <<< ====================== */
  int i, j, ii;
  int ncut;
  COORD_TYPE sigmaSq, epsilon4, epsilon24;
  COORD_TYPE rxij, ryij, rzij, rijSq, fxij, fyij, fzij;
  COORD_TYPE srij2, srij4, srij8, srij6, srij12, fij, vij, wij;
  COORD_TYPE Fxi, Fyi, Fzi, rxi, ryi, rzi;
  COORD_TYPE vCut, rcutSig, rcutSq;
  COORD_TYPE dvdr, fict;
#ifdef SOFT_SPHERE
  int PP;
#endif
  /* Local variables to implement linked list */
  int  n, nebrTab0, nebrTab1;

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
  epsilon24 = 24.0 * epsilon;
  srij2   = sigma / rcutSq;
  srij6   = srij2 * srij2 * srij2;
  srij12  = Sqr(srij6);
#ifdef SOFT_SPHERE
  dvdr = 0;
#else
  dvdr = epsilon24 * (srij6 - 2.0 * srij12) / rcutSig;
#endif
  /* initialize ab-variables */
  ncut = 0;
#ifdef SOFT_SPHERE
  PP = Oparams.PP;
#endif
  L = cbrt(Vol);
  invL = 1.0  / L;
  VLJ = 0.0; /* potential energy */
  WLJ = 0.0; /* virial function */
  /* Loop over cells */
  Vc = 0.0;
  
  for (n = 0; n < nebrTabLen; n++)
    {
      i = nebrTab[0][n]; 
      j = nebrTab[1][n];
      
      rxi = rx[i];
      ryi = ry[i];
      rzi = rz[i];
	  
      Fxi = Fx[i];
      Fyi = Fy[i];
      Fzi = Fz[i];
     
      rxij = rxi - rx[j]; /* distance between two atomes */
      ryij = ryi - ry[j];
      rzij = rzi - rz[j];
      
      rxij = rxij - L * rint(invL * rxij);      /* minimum image */
      ryij = ryij - L * rint(invL * ryij);
      rzij = rzij - L * rint(invL * rzij);
      
      rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	     
      if ( rijSq < rcutSq )/* 'rcut' is the cutoff for V */
	{
	  //rab   = sqrt(rabSq);
#ifdef SOFT_SPHERE
	  srij2   = sigmaSq / rijSq;
#if defined(MD_STATIC_PP36)
	  srij4=srij2*srij2;
          srij8=srij4*srij4;
          vij = srij8*srij8*srij2;
	  vij= vij*vij;
          wij = 36.0*vij;
#else
	  vij = pow(srij2, ((double)PP) / 2.0);
	  wij = ((double) PP) * vij;
#endif
	  VLJ = VLJ + vij;
	  /* total potential between all a-b atoms pairs */
	  WLJ = WLJ + wij; 

	  /* NOTE: If you will use a shifted-force potential then 
	     calculate the force using that potential */
	  fij   = epsilon4 * wij / rijSq;
	  /* force between two atoms */
	  fxij  = fij * rxij;          
	  fyij  = fij * ryij;
	  fzij  = fij * rzij;
#else
	  srij2   = sigmaSq / rijSq;
	  srij6   = srij2 * srij2 * srij2;
	  srij12  = Sqr(srij6);

	  vij     = srij12 - srij6;
	  //vij     = vij -  dvdr * (rij - rcut);
	  wij     = vij + srij12;
	  
	  VLJ = VLJ + vij;
	  
	  /* total potential between all a-b atoms pairs */
	  WLJ = WLJ + wij; 

	  /* NOTE: If you will use a shifted-force potential then 
	     calculate the force using that potential */
	  fij   = wij / rijSq;
	  /* force between two atoms */
	  fxij  = fij * rxij * epsilon24;         
	  fyij  = fij * ryij * epsilon24;
	  fzij  = fij * rzij * epsilon24;
#endif
	  /* Virial off-diagonal terms of atomic pressure tensor */
	  Wxy += rxij * fyij;
	  Wyz += ryij * fzij;
	  Wzx += rzij * fxij;
	  Wxx += rxij * fxij;
	  Wyy += ryij * fyij;
	  Wzz += rzij * fzij;
	  
	  /* Calculate all temrs of molecular
	     pressure tensor */
	  Fxi   = Fxi + fxij;     /* total force on an atom (a,i)*/
	  Fyi   = Fyi + fyij;
	  Fzi   = Fzi + fzij;
	  
	  Fx[j] = Fx[j] - fxij;  /* -fxab = fxba (3rd law) */
	  Fy[j] = Fy[j] - fyij;
	  Fz[j] = Fz[j] - fzij;
	  ++ncut;
	}
    
      Fx[i] = Fxi;
      Fy[i] = Fyi;
      Fz[i] = Fzi;
    }
  /* OUTER LOOP ENDS */

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
#ifdef SOFT_SPHERE
  srij2 = sigmaSq / rcutSq;
  vCut = pow(srij2, ((double)PP) / 2.0);
  Vc = VLJ - ncut * vCut;
  VLJ  = VLJ * epsilon4;
  Vc = Vc * epsilon4;
  WLJ  = WLJ * epsilon4 / 3.0;
#else
  srij2 = sigmaSq / rcutSq;
  srij6 = srij2 * srij2 * srij2;
  srij12 = srij6 * srij6;
  vCut = srij12 - srij6;
  Vc = VLJ - ncut * vCut;
  VLJ  = VLJ * epsilon4;
  Vc = Vc * epsilon4;
  WLJ  = WLJ * epsilon24 / 3.0;
#endif
}  
