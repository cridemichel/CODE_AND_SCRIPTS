#include<mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/
#define MD_DEBUG31(x) 
#define MD_USE_CALLOC
/* Allocate memory for a matrix of integers */
#ifdef EDHE_FLEX
int MD_MAX_BOND_PER_PART = 3;
#endif
long long int** AllocMatLLI(int size1, int size2)
{
  long long int** v;
  int k;
  v = (long long int**) malloc(size1 * sizeof(long long int*));
  v[0] = (long long int*) malloc(size1 * size2 * sizeof(long long int));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */
void writeAsciiPars(FILE* fs, struct pascii strutt[]);
void writeAllCor(FILE* fs, int saveAll);
void scalevels(double temp, double K);
extern void rebuildNNL(void);
extern char TXT[MSG_LEN];
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
void setToZero(COORD_TYPE* ptr, ...);
double *maxax;
#ifdef MD_MULTIPLE_LL
extern void ProcessCellCrossingMLL(void);
extern void PredictEventMLL(void);
extern void PredictEventMLL_NLL(void);
extern void rebuildMultipleLL(void);
extern int ***inCellMLL;
extern int **cellListMLL;
extern double *rcutMLL;
extern void set_cells_size(void);
extern int *cellsxMLL, *cellsyMLL, *cellszMLL, *ignoreMLL;
#endif
double calc_norm(double *vec);
#ifdef MD_PATCHY_HE
extern struct LastBumpS *lastbump;
#else
extern int *lastbump;
#endif
#ifdef MD_PROTEIN_DESIGN
extern int nativeConfYN;
#endif
#ifdef MD_ABSORPTION
int *listtmp;
#endif
extern double *lastcol;
double *axa, *axb, *axc;
double **Aip;
#ifdef EDHE_FLEX
double *a0I;
extern int is_infinite_Itens(int i);
extern int is_infinite_mass(int i);
#endif
#ifdef MD_SUPERELLIPSOID
extern double *volSQ;
#endif
#ifdef MD_SPHERICAL_WALL
int sphWall, sphWallOuter;
#endif
int *scdone;
#ifdef MD_SPOT_GLOBAL_ALLOC
#if 0
extern double ratA[NA][3], ratB[NA][3];
extern double t2arrP[6][NA], distsP[6][NA], maxddotiP[6][NA], distsOldP[6][NA];
extern int crossedP[6][NA];
#ifndef MD_BASIC_DT
extern double distsOld2P[6][NA];
#endif
extern int tocheckP[6][NA], dorefineP[6][NA], crossedP[6][NA];
extern double distsSq[NA];
#else
extern double **ratA, **ratB;
extern double *t2arrP[6], *distsP[NA], *maxddotiP[NA], *distsOldP[NA];
extern int *crossedP[6];
#ifndef MD_BASIC_DT
extern double *distsOld2P[6];
#endif
extern int *tocheckP[6], *dorefineP[6], *crossedP[6];
extern double *distsSq;
#endif
#endif
/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
double mgA, mgB, g2;
#ifdef MD_LXYZ
extern COORD_TYPE pi, invL[3], L2[3];   
#else
extern COORD_TYPE pi, invL, L2;   
#endif
#ifdef MD_GRAVITY
extern double Lz2;
#endif
extern COORD_TYPE W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB; 
extern int **tree, *inCell[3], *cellList, cellsx, cellsy, cellsz, cellRange[2*NDIM];
#ifdef MD_EDHEFLEX_OPTNNL
int *inCell_NNL[3], *cellList_NNL;
double *rxNNL, *ryNNL, *rzNNL;
#endif
#ifdef EDHE_FLEX
int *is_a_sphere_NNL;
int *is2saveArr;
#ifdef MD_ABSORP_POLY
int *oldTypeOfPart;
#endif
#endif

#ifdef MD_CALENDAR_HYBRID
extern int numevPQ, totevHQ, overevHQ;
#endif
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
#if defined(MD_HE_PARALL)
int MPIpid;
int my_rank;
int numOfProcs; /* number of processeses in a communicator */
#endif 
#ifdef MD_MULTIPLE_LL
extern int **crossevtodel;
#endif
#ifdef EDHE_FLEX
char colsFlex[][256] = {"red","green","blue", "snow","gainsboro","OldLace","linen","PapayaWhip","blanched almond",
"BlanchedAlmond","bisque","peach puff","PeachPuff","navajo white","moccasin","cornsilk","ivory","lemon chiffon",
"LemonChiffon","seashell","honeydew","mint","MintCream","azure","alice","AliceBlue","lavender","lavender",
"LavenderBlush","misty","MistyRose","dark","DarkSlateGray","dark","DarkSlateGrey","dim",
"DimGray","dim","DimGrey","slate","SlateGray","slate","SlateGrey","light","LightSlateGray","light",
"LightSlateGrey","gray","grey","light","LightGrey","light","LightGray","midnight","MidnightBlue","navy",
"navy","NavyBlue","cornflower","CornflowerBlue","dark","DarkSlateBlue","slate","SlateBlue","medium","MediumSlateBlue",
"light","LightSlateBlue","medium","MediumBlue","royal","RoyalBlue","dodger","DodgerBlue","deep",
"DeepSkyBlue","sky","SkyBlue","light","LightSkyBlue","steel","SteelBlue","light","LightSteelBlue","light",
"LightBlue","powder","PowderBlue","pale","PaleTurquoise","dark","DarkTurquoise","medium","MediumTurquoise","turquoise",
"cyan","light","LightCyan","cadet","CadetBlue","medium","MediumAquamarine","aquamarine","dark","DarkGreen",
"dark","DarkOliveGreen","dark","DarkSeaGreen","sea","SeaGreen","medium","MediumSeaGreen","light","LightSeaGreen",
"pale","PaleGreen","spring","SpringGreen","lawn","LawnGreen","chartreuse","medium","MediumSpringGreen",
"GreenYellow","lime","LimeGreen","yellow","YellowGreen","forest","ForestGreen","olive","OliveDrab",
"dark","DarkKhaki","khaki","pale","PaleGoldenrod","light","LightGoldenrodYellow","light","LightYellow","yellow",
"gold","light","LightGoldenrod","goldenrod","dark","DarkGoldenrod","rosy","RosyBrown","indian","IndianRed",
"saddle","SaddleBrown","sienna","peru","burlywood","beige","wheat","sandy","SandyBrown","tan",
"chocolate","firebrick","brown","dark","DarkSalmon","salmon","light","LightSalmon","orange","dark",
"DarkOrange","coral","light","LightCoral","tomato","orange","OrangeRed","hot","HotPink",
"deep","DeepPink","pink","light","LightPink","pale","PaleVioletRed","maroon","medium","MediumVioletRed",
"violet","VioletRed","magenta","violet","plum","orchid","medium","MediumOrchid","dark","DarkOrchid",
"dark","DarkViolet","blue","BlueViolet","purple","medium","MediumPurple","thistle","snow1","snow2",
"snow3","snow4","seashell1","seashell2","seashell3","seashell4",
"bisque1","bisque2","bisque3","bisque4","PeachPuff1","PeachPuff2","PeachPuff3","PeachPuff4",
"LemonChiffon1","LemonChiffon2","LemonChiffon3","LemonChiffon4","cornsilk1","cornsilk2","cornsilk3","cornsilk4",
"ivory1","ivory2","ivory3","ivory4","honeydew1","honeydew2","honeydew3","honeydew4","LavenderBlush1","LavenderBlush2",
"LavenderBlush3","LavenderBlush4","MistyRose1","MistyRose2","MistyRose3","MistyRose4","azure1","azure2","azure3","azure4",
"SlateBlue1","SlateBlue2","SlateBlue3","SlateBlue4","RoyalBlue1","RoyalBlue2","RoyalBlue3","RoyalBlue4","blue1","blue2",
"blue3","blue4","DodgerBlue1","DodgerBlue2","DodgerBlue3","DodgerBlue4","SteelBlue1","SteelBlue2","SteelBlue3","SteelBlue4",
"DeepSkyBlue1","DeepSkyBlue2","DeepSkyBlue3","DeepSkyBlue4","SkyBlue1","SkyBlue2","SkyBlue3","SkyBlue4","LightSkyBlue1","LightSkyBlue2",
"LightSkyBlue3","LightSkyBlue4","SlateGray1","SlateGray2","SlateGray3","SlateGray4","LightSteelBlue1","LightSteelBlue2","LightSteelBlue3","LightSteelBlue4",
"LightBlue1","LightBlue2","LightBlue3","LightBlue4","LightCyan1","LightCyan2","LightCyan3","LightCyan4","PaleTurquoise1","PaleTurquoise2",
"PaleTurquoise3","PaleTurquoise4","CadetBlue1","CadetBlue2","CadetBlue3","CadetBlue4","turquoise1","turquoise2","turquoise3","turquoise4",
"cyan1","cyan2","cyan3","cyan4","DarkSlateGray1","DarkSlateGray2","DarkSlateGray3","DarkSlateGray4","aquamarine1","aquamarine2",
"aquamarine3","aquamarine4","DarkSeaGreen1","DarkSeaGreen2","DarkSeaGreen3","DarkSeaGreen4","SeaGreen1","SeaGreen2","SeaGreen3","SeaGreen4",
"PaleGreen1","PaleGreen2","PaleGreen3","PaleGreen4","SpringGreen1","SpringGreen2","SpringGreen3","SpringGreen4","green1","green2",
"green3","green4","chartreuse1","chartreuse2","chartreuse3","chartreuse4","OliveDrab1","OliveDrab2","OliveDrab3","OliveDrab4",
"DarkOliveGreen1","DarkOliveGreen2","DarkOliveGreen3","DarkOliveGreen4","khaki1","khaki2","khaki3","khaki4","LightGoldenrod1","LightGoldenrod2",
"LightGoldenrod3","LightGoldenrod4","LightYellow1","LightYellow2","LightYellow3","LightYellow4","yellow1","yellow2","yellow3","yellow4",
"gold1","gold2","gold3","gold4","goldenrod1","goldenrod2","goldenrod3","goldenrod4","DarkGoldenrod1","DarkGoldenrod2",
"DarkGoldenrod3","DarkGoldenrod4","RosyBrown1","RosyBrown2","RosyBrown3","RosyBrown4","IndianRed1","IndianRed2","IndianRed3","IndianRed4",
"sienna1","sienna2","sienna3","sienna4","burlywood1","burlywood2","burlywood3","burlywood4","wheat1","wheat2",
"wheat3","wheat4","tan1","tan2","tan3","tan4","chocolate1","chocolate2","chocolate3","chocolate4",
"firebrick1","firebrick2","firebrick3","firebrick4","brown1","brown2","brown3","brown4","salmon1","salmon2",
"salmon3","salmon4","LightSalmon1","LightSalmon2","LightSalmon3","LightSalmon4","orange1","orange2","orange3","orange4",
"DarkOrange1","DarkOrange2","DarkOrange3","DarkOrange4","coral1","coral2","coral3","coral4","tomato1","tomato2",
"tomato3","tomato4","OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","red1","red2","red3","red4",
"DeepPink1","DeepPink2","DeepPink3","DeepPink4","HotPink1","HotPink2","HotPink3","HotPink4","pink1","pink2",
"pink3","pink4","LightPink1","LightPink2","LightPink3","LightPink4","PaleVioletRed1","PaleVioletRed2","PaleVioletRed3","PaleVioletRed4",
"maroon1","maroon2","maroon3","maroon4","VioletRed1","VioletRed2","VioletRed3","VioletRed4","magenta1","magenta2",
"magenta3","magenta4","orchid1","orchid2","orchid3","orchid4","plum1","plum2","plum3","plum4",
"MediumOrchid1","MediumOrchid2","MediumOrchid3","MediumOrchid4","DarkOrchid1","DarkOrchid2","DarkOrchid3","DarkOrchid4","purple1","purple2",
"purple3","purple4","MediumPurple1","MediumPurple2","MediumPurple3","MediumPurple4","thistle1","thistle2","thistle3","thistle4","DarkBlue","dark","DarkCyan",
"dark","DarkMagenta","dark","DarkRed","light","LightGreen", ""};
int numcols;
int *mapbondsaFlex, *mapbondsbFlex, nbondsFlex;
double *mapBheightFlex, *mapBhinFlex, *mapBhoutFlex, *mapSigmaFlex; 
double *t2arr, *distsOld, *dists, *distsOld2, *maxddoti;
int *crossed, *tocheck, *dorefine, *negpairs, dofTot;
int **mapbondsaFlexS, **mapbondsbFlexS, *nbondsFlexS;
double **mapBheightFlexS, **mapBhinFlexS, **mapBhoutFlexS, **mapSigmaFlexS; 
#endif
extern double **matrix(int n, int m);
extern int *ivector(int n);
extern double *vector(int n);
int poolSize;
int parnumA, parnumB;
#ifdef MD_PATCHY_HE
#ifdef MD_LL_BONDS
long long int *bondscache, **bonds;
int *numbonds;
#else
int *bondscache, *numbonds, **bonds;
#endif
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
#ifndef EDHE_FLEX
extern double spXYZ_A[MD_STSPOTS_A][3];
extern double spXYZ_B[MD_STSPOTS_B][3];
#endif
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
#ifdef EDHE_FLEX
void check_these_bonds(int i, int j, double *shift, double t)
{
  int amin, warn, bmin, nbonds, nn;
  double dist;
  assign_bond_mapping(i,j);

  printf(">>>>>>>>>>QUI\n");
  dist = calcDistNegSP(t, 0.0, i, j, shift, &amin, &bmin, dists, -1);
#ifdef EDHE_FLEX
  nbonds = nbondsFlex;
#else
  nbonds = MD_PBONDS;
#endif
  for (nn=0; nn < nbonds; nn++)
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
		     printf("NA*NA*i+a*NA+b=%d\n", NANA*i+mapbondsa[nn]*NA+mapbondsb[nn]);
		    )
	}
      else if (dists[nn]>0.0 && 
	       fabs(dists[nn])> OprogStatus.epsd && 
	       bound(i,j,mapbondsa[nn], mapbondsb[nn]))
	{
	  warn = 2;
	  printf("[PredictEvent]wrong number of bonds between %d(%d) and %d(%d) nbonds=%d nn=%d\n",
		 i, mapbondsa[nn], j, mapbondsb[nn], nbonds, nn);
	  printf("sigmaSticky=%f dist=%.15G\n", 0.5*(typesArr[typeOfPart[i]].spots[mapbondsa[nn]-1].sigma+typesArr[typeOfPart[j]].spots[mapbondsb[nn]-1].sigma),
		 dists[nn]);
	  exit(-1);	
	}
    }
}
#endif
#ifdef EDHE_FLEX
#if 0
int all_bonds_of_same_type(int i, int ni)
{
  int kk, jj, jj2, aa, bb, cc;
  cc = 0;
  for (kk=0; kk < numbonds[i]; kk++)
    {
      jj = bonds[i][kk] / (NANA);
      jj2 = bonds[i][kk] % (NANA);
      aa = jj2 / NA;
      bb = jj2 % NA;
      if ((intersArr[ni].type1 == jj && intersArr[ni].spot1+1 == aa) ||
	  (intersArr[ni].type2 == jj && intersArr[ni].spot2+1 == bb))
	    cc++;
    }
  return cc;
}
int nmax_reached(int i, int j, int a, int b)
{
  int kk, type1, type2, ni;
  type1=typeOfPart[i];
  type2=typeOfPart[j];
  for (ni=0; ni < Oparams.ninters; ni++)
    {
      if (intersArr[ni].type1 == type1 && intersArr[ni].spot1+1 == a && intersArr[ni].type2 == type2 &&
	  intersArr[ni].spot2+1 == b)
	{
	  if (intersArr[ni].nmax == all_bonds_of_same_type(i,ni))
	    return 1;	    
	}
      else if (intersArr[ni].type1 == type2 && intersArr[ni].spot1+1 == b && intersArr[ni].type2 == type1 &&
	       intersArr[ni].spot2+1 == a )
	{
	  if (intersArr[ni].nmax == all_bonds_of_same_type(i,ni))
	    return 1;
	}
    }
  return 0;
}
#endif
#endif
#ifdef MD_CHECK_POINT
void read_check_point(void)
{
  if ((fs = fopenMPI(fn, "r")) == NULL)
    {
      sprintf(msgStrA, "Problem opening checkpoint file %s ", fn);
      mdMsg(ALL, NOSYS, "read_check_point", "ERROR", NULL,
	    msgStrA,
	    NULL);
      exit(-1);
    }

  readAsciiPars(fs, opro_ascii);
  readAsciiPars(fs, opar_ascii);
  /* Entrambe queste macro sono definite nel file mono_DPT.h */
  /* read up to coordinates begin */
  readAllCorCP(fs);

  fclose(fs);

}
void restart_from_check_point(void)
{
  read_check_point();
  StartRun();
  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
  if (OprogStatus.rescaleTime > 0)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  else
    OprogStatus.scalevel = 0;
  ScheduleEvent(-1, ATOM_LIMIT + 11, OprogStatus.bigDt);
}
#endif
#ifdef MD_GHOST_IGG
extern int areGhost(int i, int j);
#endif
#ifdef MD_MULTIPLE_LL
extern void check_all_bonds_MLL(void);
#endif
extern void check_all_bonds_NLL(void);
void check_all_bonds(void)
{
  int nn, warn, amin, bmin, i, j, nb, nbonds;
  double shift[3]={0.0,0.0,0.0}, dist;
#ifndef EDHE_FLEX
  double dists[MD_PBONDS];
#endif
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k;
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravit� � diretta lungo z negativo */ 

#if 1
  if (OprogStatus.useNNL)
    {
      check_all_bonds_NLL();
      return;
    }
#endif
#ifdef MD_MULTIPLE_LL
  /* N.B. 24/05/10: questa funzione va riscritta in multiple_ll.c */
  if (OprogStatus.multipleLL)
    {
      check_all_bonds_MLL();
      return;
    }
#endif
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  
  warn = 0;
  for ( i = 0; i < Oparams.parnum; i++)
    {
      nb = 0;
#if defined(MD_ABSORPTION) && defined(MD_SPHERICAL_WALL)
      if (i==sphWall || i==sphWallOuter)
	{
  	  for (j = 0; j < Oparams.parnum; j++)
	    {
	      if (j==sphWall || j==sphWallOuter)
		continue;
	      assign_bond_mapping(i,j);
	      if (nbondsFlex == 0)
		continue;
	     warn = check_bonds_ij(i, j, shift); 
	    }
	  continue;
	}
#endif
      for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_EDHEFLEX_WALL
      if (inCell[2][i] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (inCell[2][i] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
#endif 
      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCell[2][i] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      printf("BOHHHH\n");
	      jZ = cellsz - 1;    
#ifdef MD_LXYZ
	      shift[2] = - L[2];
#else
	      shift[2] = - L;
#endif
	    } 
	  else if (jZ == cellsz) 
	    {
	      jZ = 0;    
#ifdef MD_LXYZ
	      shift[2] = L[2];
#else
	      shift[2] = L;
#endif
	    }
	  for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	    {
	      jY = inCell[1][i] + iY;    
	      shift[1] = 0.0;
	      if (jY == -1) 
		{
		  jY = cellsy - 1;    
#ifdef MD_LXYZ
		  shift[1] = -L[1];
#else
		  shift[1] = -L;
#endif
		} 
	      else if (jY == cellsy) 
		{
		  jY = 0;    
#ifdef MD_LXYZ
		  shift[1] = L[1];
#else
		  shift[1] = L;
#endif
		}
	      for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
		{
		  jX = inCell[0][i] + iX;    
		  shift[0] = 0.0;
		  if (jX == -1) 
		    {
		      jX = cellsx - 1;    
#ifdef MD_LXYZ
		      shift[0] = - L[0];
#else
		      shift[0] = - L;
#endif
		    } 
		  else if (jX == cellsx) 
		    {
		      jX = 0;   
#ifdef MD_LXYZ
		      shift[0] = L[0];
#else
		      shift[0] = L;
#endif
		    }
		  j = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
		  for (j = cellList[j]; j > -1; j = cellList[j]) 
		    {
		      if (i == j)
			continue;
#ifdef MD_GHOST_IGG
		      /* do not check ghost particles, it is meaningless! */
		      if (Oparams.ghostsim && areGhost(i, j))
			continue;
#endif
#ifndef EDHE_FLEX
		      if (! ((i < Oparams.parnumA && j >= Oparams.parnumA)||
			     (i >= Oparams.parnumA && j < Oparams.parnumA)))
			  continue;
#endif
		      check_shift(i, j, shift);
		      assign_bond_mapping(i,j);
		      dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
#ifdef EDHE_FLEX
		      nbonds = nbondsFlex;
#else
		      nbonds = MD_PBONDS;
#endif
		      for (nn=0; nn < nbonds; nn++)
			{
			  if (dists[nn]<0.0 && fabs(dists[nn])>OprogStatus.epsd 
			      && !bound(i,j,mapbondsa[nn], mapbondsb[nn]) )
			  // && fabs(dists[nn]-Oparams.sigmaSticky)>1E-4)
			    {
			      warn=1;
		  	      MD_DEBUG31(
			      //printf("dists[1]:%.15G\n", dists[1]);
			      printf("[dist<0]dists[%d]:%.15G\n", nn, dists[nn]);
			      printf("i=%d j=%d %d %d\n", i, j, mapbondsa[nn], mapbondsb[nn]);
			      printf("NA*NA*i+a*NA+b=%lld\n", NANA*i+mapbondsa[nn]*NA+mapbondsb[nn]);
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
				   fabs(dists[nn]) > OprogStatus.epsd && 
				   bound(i,j,mapbondsa[nn], mapbondsb[nn]))
			    {
			      warn = 2;
			      printf("wrong number of bonds between %d(%d) and %d(%d) nbonds=%d nn=%d\n",
				     i, mapbondsa[nn], j, mapbondsb[nn], nbonds, nn);

			      printf("r=%f %f %f - %f %f %f\n", rx[i], ry[i], rz[i], rx[j], ry[j], rz[j]);
			      printf("[dist>0]dists[%d]:%.15G\n", nn, dists[nn]);
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
#ifdef MD_ALLOW_ONE_IGG_BOND
      if (warn==1 && get_igg_bonds(i, j)==1)
	continue;
#endif

      if (warn)
	{
	  mdPrintf(ALL, "[WARNING] wrong number of bonds\n", NULL);
	  sprintf(TXT,"[WARNING] Number of bonds for molecules %d incorrect\n", i);
	  mdPrintf(ALL, TXT, NULL);
	  sprintf(TXT,"Step N. %d time=%.15G\n", Oparams.curStep, Oparams.time);
	  mdPrintf(ALL, TXT, NULL);
#ifdef MD_LL_BONDS
	  printf("numbonds[%d]:%d bonds[][]:%lld\n", i, numbonds[i], bonds[i][0]);
#else
	  printf("numbonds[%d]:%d bonds[][]:%d\n", i, numbonds[i], bonds[i][0]);
#endif
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
    //  if (warn)
//	break;

    }
#ifdef MD_CHECK_POINT
  saveFullStore("ED_CHECK_POINT\n");
#endif
 
}
#endif
void ScheduleEvent (int idA, int idB, double tEvent); 
void check_coord(void)
{
  int i;
  for (i = 0; i < Oparams.parnum; i++)
#ifdef MD_LXYZ
    if (fabs(rx[i]) > L[0]*0.5 || fabs(ry[i])>L[1]*0.5 || fabs(rz[i]) > L[2]*0.5)
      {
	printf("%d is out of box!\n", i);
	exit(-1);
      }
#else
    if (fabs(rx[i]) > L*0.5 || fabs(ry[i])>L*0.5 || fabs(rz[i]) > L*0.5)
      {
	printf("%d is out of box!\n", i);
	exit(-1);
      }
#endif
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
#ifdef MD_LXYZ
  double Cell[3], Cell2[3];
#else
  double Cell, Cell2;
#endif
#ifdef MD_GRAVITY
  double Cellz, Cell2z;
#endif
  int i, ix, iy, iz, iref, ii;
  double bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
#ifdef MD_GRAVITY
  /* NOTA:
     Ncz � tanto pi� grande quanto maggiore � Lz  rispetto a L,
     cio� cos� si tiene conto che la scatola pu� avere un'altezza 
     differente dagli altri due lati */
  Nc = ceil( pow ( (L/Lz)*((double)Nm) / 4.0, 1.0/3.0) ); 
  Ncz = ceil((Lz/L)*Nc);	      
  printf("Nc: %d Ncz:%d\n", Nc, Ncz);
#else
  Ncz = Nc = ceil(  pow( ((double)Nm)/4.0, 1.0/3.0 )  );
  printf("Nc: %d\n", Nc);
#endif
  /* Calculate the side of the unit cell */
#ifdef MD_LXYZ
  for (ii=0; ii < 3; ii++)
    {
      Cell[ii] = L[ii] / ((double) Nc); /* unit cell length */
      Cell2[ii] = 0.5 * Cell[ii];              /* half unit cell length */
    }
#else
  Cell  = L / ((double) Nc); /* unit cell length */
  Cell2 = 0.5 * Cell;              /* half unit cell length */
#ifdef MD_GRAVITY
  Cellz = Lz / ((double) Nc);
  Cell2z = 0.5 * Cellz;
#endif
#endif
  /* Sublattice A */
  rRoot3 = 1.0 / sqrt(3.0);
  bx[0] =  0.0;
  by[0] =  0.0;
  bz[0] =  0.0;
  
  /*  Sublattice B */
#ifdef MD_LXYZ
  bx[1] =  Cell2[0];
  by[1] =  Cell2[1];
#else
  bx[1] =  Cell2;
  by[1] =  Cell2;
#endif
  bz[1] =  0.0;
  /* Sublattice C */
  
  bx[2] =  0.0;
#ifdef MD_LXYZ
  by[2] =  Cell2[1];
#else
  by[2] =  Cell2;
#endif
#ifdef MD_LXYZ
  bz[2] = Cell2[2];
#else
#ifdef MD_GRAVITY
  bz[2] = Cell2z;
#else
  bz[2] = Cell2;
#endif
#endif  
  /* Sublattice D */
#ifdef MD_LXYZ
  bx[3] =  Cell2[0];
#else
  bx[3] =  Cell2;
#endif
  by[3] =  0.0;
#ifdef MD_LXYZ
  bz[3] =  Cell2[2];
#else
#ifdef MD_GRAVITY
  bz[3] = Cell2z;
#else
  bz[3] =  Cell2;
#endif  
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
#ifdef MD_LXYZ
		  rx[ii+iref] = bx[iref] + Cell[0] * ((double) ix);
		  ry[ii+iref] = by[iref] + Cell[1] * ((double) iy);
		  rz[ii+iref] = bz[iref] + Cell[2] * ((double) iz);
#else
		  rx[ii+iref] = bx[iref] + Cell * ((double) ix);
		  ry[ii+iref] = by[iref] + Cell * ((double) iy);
#ifdef MD_GRAVITY
		  rz[ii+iref] = bz[iref] + Cellz * ((double) iz);
#else
		  rz[ii+iref] = bz[iref] + Cell * ((double) iz);
#endif
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
#ifdef MD_LXYZ
      rx[i] = rx[i] - 0.5 * L[0] + ranf()*1E-7; 
      ry[i] = ry[i] - 0.5 * L[1] + ranf()*1E-7;
      rz[i] = rz[i] - 0.5 * L[2] + ranf()*1E-7;
#else
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
#ifdef MD_RAND48
  return drand48();
#else
  return rand() / ( (COORD_TYPE) RAND_MAX );
#endif
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
#ifdef EDHE_FLEX
void resetCM(int onlyz)
{
  COORD_TYPE sumx, sumy, sumz, RCMx, RCMy, RCMz;
  COORD_TYPE Mtot, mass;
  int i;
  /* Remove net momentum, to have a total momentum equals to zero */
  Mtot = 0.0;
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  for(i=0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	{
	  continue;
	}
      mass = typesArr[typeOfPart[i]].m;
      Mtot += mass;
      sumx = sumx + vx[i]*mass;
      sumy = sumy + vy[i]*mass;
      sumz = sumz + vz[i]*mass;
    }
     
  sumx = sumx / Mtot; 
  sumy = sumy / Mtot;
  sumz = sumz / Mtot;

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i = 0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	{
	  continue;
	}
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

  for(i = 0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	{
	  continue;
	}
      mass = typesArr[typeOfPart[i]].m;
      RCMx += rx[i]*mass; /* Here RCM is the center of mass of the box */
      RCMy += ry[i]*mass;
      RCMz += rz[i]*mass;
    }
  
  RCMx /= Mtot;
  RCMy /= Mtot;
  RCMz /= Mtot;

  for(i = 0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	continue;
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
#else
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
#endif
void comvel_brown (COORD_TYPE temp, COORD_TYPE *m)
{
#ifdef EDHE_FLEX
  double rTemp, mass;
#else
  COORD_TYPE rTemp[2] ;
#endif
  /*COORD_TYPE Px, Py, Pz;*/
  int i;
#ifdef EDHE_FLEX
  for (i = 0; i < Oparams.parnum; i++)
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
      /* brownian = 0 non fa niente, 1 = brownian motion, 2 = rescale vels (see scalevels()) */
      if (typesArr[typeOfPart[i]].brownian==1)
	{
	  mass = typesArr[typeOfPart[i]].m;
	  rTemp = sqrt(temp / mass);  
	  vx[i] = rTemp * gauss(); 
	  vy[i] = rTemp * gauss();
#ifndef MD_EDHEFLEX_2D
	  vz[i] = rTemp * gauss();
#endif
	}	  
    }
#else
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
#endif 
}
#ifdef EDHE_FLEX
double ranfRandom(void)
{
  /*  Returns a uniform random variate in the range 0 to 1 (excluding 0).         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return (1.0-(((double)rand()) / (((double) RAND_MAX) + 1)));
}

void angvelMB(void)
{
  /* NOTA 16/02/11: questa routine ora presuppone che il corpo
     rigido sia simmetrico rispetto 
     all'asse x nel riferimento del corpo rigido
     cosicch� in tale rif. la vel. ang. intorno a tale asse sia zero 
     (l'ho scritta per le cokecan) */
  int i, a;
  double sp, pi, inert;                 /* momentum of inertia of the molecule */
  double norm, osq, o, mean, symax[3];
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, ww[3], wsz;

  pi = acos(0)*2; 

  /* N.B. QUESTO E' SBAGLIATO PERCHE' LA DISTRIBUZIONE
     E' COME QUELLA TRASLAZIONALE (VEDI LANDAU) 
     CORREGGERE!!!! */
  for (i = 0; i < Oparams.parnum; i++)
    {
      inert = typesArr[typeOfPart[i]].I[0];
      mean = sqrt(Oparams.T / inert);
      wx[i] = mean*gauss();
      wy[i] = mean*gauss();
      wz[i] = mean*gauss();
      sp = wx[i]*R[i][0][0]+wy[i]*R[i][0][1]+wz[i]*R[i][0][2];
      wx[i] -= sp*R[i][0][0];
      wy[i] -= sp*R[i][0][1];
      wz[i] -= sp*R[i][0][2];
      Mx[i] = wx[i]*inert;
      My[i] = wy[i]*inert;
      Mz[i] = wz[i]*inert;
      //printf("w=%f %f %f\n", wx[i], wy[i], wz[i]);
    }
}

/* ========================== >>> comvel <<< =============================== */
void comvel (int Nm, COORD_TYPE temp, COORD_TYPE *m, int resetCM)
{
  double rTemp, sumx, sumy, sumz, RCMx, RCMy, RCMz, Mtot, mass;
  /*COORD_TYPE Px, Py, Pz;*/
  int i;
  Mtot = 0.0;

  /* variance of the velocities distribution function, we assume k = 1 */ 
  K = 0;
  for (i = 0; i < Oparams.parnum; i++)
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

      if (is_infinite_mass(i))
	{
	  continue;
	}
      mass = typesArr[typeOfPart[i]].m;
      Mtot += mass;
      rTemp = sqrt(temp / mass);  
      vx[i] = rTemp * gauss(); 
      vy[i] = rTemp * gauss();
#ifndef MD_EDHEFLEX_2D
      vz[i] = rTemp * gauss();
#endif
      //printf("rank[%d] vx[%d]: %f\n", my_rank, i, vx[i]);
      /* gauss() is a gaussian variable with mean = 0 and variance = 1, that is
                               2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
      K = K + 0.5 * Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));
    }
 
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  for(i = 0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	continue;
      /* (sumx, sumy, sumz) is the total momentum */ 
      mass = typesArr[typeOfPart[i]].m;
      sumx = sumx + vx[i]*mass;
      sumy = sumy + vy[i]*mass;
      sumz = sumz + vz[i]*mass;
      //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
    }
 
  sumx = sumx / Mtot; 
  sumy = sumy / Mtot;
  sumz = sumz / Mtot;

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i = 0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	continue;
      vx[i] = vx[i] - sumx;
      vy[i] = vy[i] - sumy;
#ifndef MD_EDHEFLEX_2D
      vz[i] = vz[i] - sumz;
#endif
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
      if (is_infinite_mass(i))
	continue;
      mass = typesArr[typeOfPart[i]].m;
      RCMx += rx[i]*mass; /* Here RCM is the center of mass of the box */
      RCMy += ry[i]*mass;
      RCMz += rz[i]*mass;
    }
  

  RCMx /= Mtot;
  RCMy /= Mtot;
  RCMz /= Mtot;
  for(i=0; i < Oparams.parnum; i++)
    {
      if (is_infinite_mass(i))
	continue;
      //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
      rx[i] -= RCMx;
      ry[i] -= RCMy;
#ifndef MD_GRAVITY
#ifndef MD_EDHEFLEX_2D
      rz[i] -= RCMz;
#endif
#endif
    }
  /* Now the center of mass of the box is in the origin */
}
#else
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
#endif
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
#ifndef EDHE_FLEX
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
void calc_omega(int i, double *wwx, double *wwy, double *wwz);
#endif
/* =========================== >>> usrInitBef <<< ========================== */
#ifdef MD_PROTEIN_DESIGN
void read_native_conf(void);
#endif
#ifdef MD_DYNAMIC_OPROG
int dyn_alloc_oprog(void);
void set_dyn_ascii(void);
#endif

void usrInitBef(void)
{
  int i, k;
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
#ifdef MD_DYNAMIC_OPROG
  OprogStatus.ptr = NULL;
  OprogStatus.len = 0;
#endif
#ifdef MD_DYNAMIC_OPROG
  OprogStatus.dyn_alloc_oprog = dyn_alloc_oprog;
  OprogStatus.set_dyn_ascii = set_dyn_ascii;
#endif

#ifdef MD_LXYZ
    L[0] = L[1] = L[2] = 9.4;
#else
    L = 9.4;
#endif
#ifdef MD_PATCHY_HE
    OprogStatus.autocat = 0;
    OprogStatus.mac = 1.0;
    OprogStatus.k0 = 0.0;
    OprogStatus.k1 = 0.0;
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
#ifdef MD_POLYDISP
#ifdef MD_POLYDISP_XYZ
    OprogStatus.polydispX = OprogStatus.polydispY = OprogStatus.polydispZ = 0.0;
#else
    OprogStatus.polydisp = 0.0;
#endif
    OprogStatus.polycutoff = 5.0;
#endif
#ifdef EDHE_FLEX
#ifdef MD_SCALEPHI_STAGES
    /* 0 = grow all particles without any order, 1 = grow particles from type=0 to last type
       In future we could implement a method 2 that grow particles from bigger to smaller (according to type params) */ 
    OprogStatus.growthType = 0;
#endif
    OprogStatus.frozenDOF = 3;
    Oparams.saveBonds = 1;
    Oparams.maxbondsSaved = -1;
#endif
#ifdef MD_EDHEFLEX_WALL
    OprogStatus.hardwall = 0;
#endif
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
#ifdef EDHE_FLEX
    OprogStatus.optbm = 1;
#endif
#ifdef MD_EDHEFLEX_OPTNNL
    OprogStatus.optnnl = 0;
#endif
#ifdef EDHE_FLEX
    OprogStatus.stripStore = 0;
    strcpy(OprogStatus.par2save, "ALL"); 
    OprogStatus.xi = 0.00001;
    OprogStatus.Tf = 0.01;
#endif
#ifdef MD_GHOST_IGG
    Oparams.ghostsim = 0; /* 0 = no ghost IgG when they are bound, 1 = ghost */
#endif
#ifdef MD_RABBIT
    OprogStatus.time_limit = 100000000000.0;
    OprogStatus.first_time = 0.0;
    OprogStatus.rhozBinSize = 0.79/2.0+4.0+0.9+0.4+0.112+0.25;
#endif
#ifdef MD_BIG_DT
    OprogStatus.bigDt = -1.0;
#endif
    OprogStatus.avngS     = 0;
    OprogStatus.avngPress = 0;
    OprogStatus.avngTemp  = 0;
    OprogStatus.scalevel = 0;
    OprogStatus.endtime = 0;
    OprogStatus.rescaleTime = 1.0;
    OprogStatus.brownian = 0;
#if defined(MD_INELASTIC) || defined(MD_GRAVITY)
    Oparams.partDiss = 1.0;
    OprogStatus.tc = 1E-5;
#endif
#ifdef MD_GRAVITY
    OprogStatus.taptau = 0.0;
    OprogStatus.rzup = 0.0;
    OprogStatus.expandFact= 1.0;
    OprogStatus.quenchend = 0.0;
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
#ifndef MD_DYNAMIC_OPROG
    for (i = 0; i < MAXPAR; i++)
      {
	OprogStatus.lastcolltime[i] = 0.0;
	OprogStatus.sumox[i] = 0.0;
	OprogStatus.sumoy[i] = 0.0;
	OprogStatus.sumoz[i] = 0.0;
#ifdef MD_CALC_DPP
	OprogStatus.sumdx[i] = 0.0;
	OprogStatus.sumdy[i] = 0.0;
	OprogStatus.sumdz[i] = 0.0;
	OprogStatus.lastu1x[i] = 0.0;
	OprogStatus.lastu1y[i] = 0.0;
	OprogStatus.lastu1z[i] = 0.0;
	OprogStatus.lastu2x[i] = 0.0;
	OprogStatus.lastu2y[i] = 0.0;
	OprogStatus.lastu2z[i] = 0.0;
	OprogStatus.lastu3x[i] = 0.0;
	OprogStatus.lastu3y[i] = 0.0;
	OprogStatus.lastu3z[i] = 0.0;
#endif
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
    OprogStatus.epsd = 0.0005;
    OprogStatus.epsdNL = -1.0;
    OprogStatus.epsdSD = -1.0;
    OprogStatus.epsdGDO = -1.0;
    OprogStatus.h = 1E-10;
#ifdef MD_PATCHY_HE
    Oparams.nmax = 1;
    Oparams.Dr = 0.0;
    Oparams.theta = 0.54;
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
    /* NOTA: gli epsd NL sono settati a -1.0 poich�
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
#ifdef MD_CALENDAR_HYBRID
    OprogStatus.scaleHQ = 50;
    OprogStatus.nlistsHQ = 50000;  
    OprogStatus.adjustHQ = 0;
    OprogStatus.baseIndex = 0;
    OprogStatus.curIndex = 0;
    OprogStatus.overthrHQ = 20;
#endif
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
#ifdef MD_EDHEFLEX_WALL
    OprogStatus.epsdPlane=-1.0;
    OprogStatus.epsdFastPlane=-1.0;
    OprogStatus.n_gauleg = 40;
#endif
#ifdef EDHE_FLEX
    Oparams.ninters = 0;
    Oparams.nintersIJ = 0;
#endif   
#ifdef MD_PROTEIN_DESIGN
    strcpy(OprogStatus.nativeConf,"");
#endif
#ifdef MD_ABSORPTION
#ifdef MD_SPHERICAL_WALL
    OprogStatus.halfsolidangle = 0;
#endif
    OprogStatus.bufHeight = 1.2;
#endif
#ifdef MD_MULTIPLE_LL
    OprogStatus.multipleLL = 0;
    OprogStatus.rcutfactMLL = 1.01;
#endif
    maxcoll=-1;
#ifdef MC_SIMUL
    OprogStatus.adjuststepMC=1;
    OprogStatus.dthetaMC=0.1;
    OprogStatus.deltaMC=0.1;
    OprogStatus.ensembleMC=0; /* 0 = NVT 1=NPT */
#endif
}
extern void check (int *overlap, double *K, double *V);
double *atomTime, *treeTime, *treeRxC, *treeRyC, *treeRzC;
extern void PredictEvent(int, int);
extern void InitEventList(void);
void find_conciding_spots(void);
void find_bonds(void);

void StartRun(void)
{
  int j, k, n;
  
  find_conciding_spots();
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
     rebuildMultipleLL(); 
    }
  else
    {
      for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
	cellList[j] = -1;
      /* -1 vuol dire che non c'� nessuna particella nella cella j-esima */
      for (n = 0; n < Oparams.parnum; n++)
	{
#ifdef MD_SPHERICAL_WALL
	  if (n==sphWall)
	    {
	      cellList[sphWall]=-1;
	      continue;
	    }
	  if (n==sphWallOuter)
	    {
	      cellList[sphWallOuter]=-1;
	      continue;
	    }
#endif
	  atomTime[n] = Oparams.time;
	  //printf("qui n=%d %.15G %.15G %.15G\n", n, rx[n], ry[n], rz[n]);
#ifdef MD_LXYZ
	  inCell[0][n] =  (rx[n] + L2[0]) * cellsx / L[0];
	  inCell[1][n] =  (ry[n] + L2[1]) * cellsy / L[1];
	  inCell[2][n] =  (rz[n] + L2[2]) * cellsz / L[2];
#else
	  inCell[0][n] =  (rx[n] + L2) * cellsx / L;
	  inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
	  inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
	  inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
#endif
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
    }
#else
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList[j] = -1;
  /* -1 vuol dire che non c'� nessuna particella nella cella j-esima */
  for (n = 0; n < Oparams.parnum; n++)
    {
#ifdef MD_SPHERICAL_WALL
      if (n==sphWall)
	{
	  cellList[sphWall]=-1;
	  continue;
	}
      if (n==sphWallOuter)
	{
	  cellList[sphWallOuter]=-1;
	  continue;
	}
#endif
      atomTime[n] = Oparams.time;
      //printf("qui n=%d %.15G %.15G %.15G\n", n, rx[n], ry[n], rz[n]);
#ifdef MD_LXYZ
      inCell[0][n] =  (rx[n] + L2[0]) * cellsx / L[0];
      inCell[1][n] =  (ry[n] + L2[1]) * cellsy / L[1];
      inCell[2][n] =  (rz[n] + L2[2]) * cellsz / L[2];
#else
      inCell[0][n] =  (rx[n] + L2) * cellsx / L;
      inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
      inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
      inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
#endif
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
#endif
  InitEventList();
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  if (OprogStatus.useNNL)
    rebuildNNL();
#ifdef EDHE_FLEX
  /* N.B. 28/05/2010: notare che gli spot vanno cercati 
     dopo l'inizializzazione delle NNL e delle LL 
     ma prima di costruire il calendario (altrimenti
     mancherebbero tutti gli eventi relativi alle collisioni
     tra gli spot! */
  if (Oparams.maxbondsSaved==-1)
    {
      //printf("nbonds[1841]=%d\n",numbonds[1841]);
      printf("[INFO] finding spots\n");
      if (Oparams.ninters != 0)
	find_bonds();
   }

  if (Oparams.saveBonds && Oparams.maxbondsSaved==-1)
    Oparams.maxbondsSaved = OprogStatus.maxbonds;
#elif defined(MD_PATCHY_HE)
  find_bonds();
#endif

  for (n = 0; n < Oparams.parnum; n++)
    {
      if (OprogStatus.useNNL)
	PredictEventNNL(n, -2);
      else
	PredictEvent(n, -2);
    }
#if 0
    {
      
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
#ifdef EDHE_FLEX
int get_num_pbonds(int i, int j)
{
#if 1
  return nbondsFlex;
#else
  int ni, type1, type2, a;
  type1 = typeOfPart[i];
  type2 = typeOfPart[j];
  for (ni=0; ni < Oparams.ninters; ni++)
    {
      if ((intersArr[ni].type1 == type1 && intersArr[ni].type2 == type2) ||
	  (intersArr[ni].type1 == type2 && intersArr[ni].type1 == type1))
	{
	  a+=1;
	}	
    }
  for (ni=0; ni < Oparams.nintersIJ; ni++)
    {
      if ((intersArrIJ[ni].i == i && intersArrIJ[ni].j == j) ||
	  (intersArrIJ[ni].i == j && intersArrIJ[ni].j == i))
	{
	  a+=1;
	}	
    }

  return a;
#endif
}
#endif
#ifdef MD_PATCHY_HE
extern void build_atom_positions(void);
extern void add_bond(int na, int n, int a, int b);
extern double calcpotene(void);
#ifndef EDHE_FLEX
int get_num_pbonds(int i, int j)
{
  return MD_PBONDS;
}
#endif
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
  //printf("to ang mom ref sys\n");
  //print_matrix(RM[i],3);
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
  print_matrix(R[i],3);
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
#ifdef EDHE_FLEX
void calc_encpp(void)
{
  int kk;
  double v[3], norm;
  int pt, sp;
  double com[3]={0.0,0.0,0.0};

  for (pt = 0; pt < Oparams.ntypes; pt++)
    {
#ifdef MD_SPHERICAL_WALL
      if (pt == Oparams.ntypes-1)
	continue;
      if (pt == Oparams.ntypes-2)
	continue;
#endif
#ifdef MD_EDHEFLEX_OPTNNL
      /* N.B. 24/05/10: per usare l'rcut ottimizzato senza l'NNL bisognerebbe
	 predire il crossing del centro della NNL e non del centro di massa
	 dell'ellissoid/SQ.
	 Nel caso infatti in cui si usano le NNL le linked list ottimizzate
	 vengono ricostruite prima di ogni rebuild delle NNL tenendo conto dello 
	 shift del centro geometrico. Nel caso per� in cui non si usino le NNL
	 per tener conto di tale shift si dovrebbe predire l'attraversamento 
	 delle celle per tale centro e non per il centro degli ellissoidi/SQ, ma ancora 
	 questo non viene fatto.
       */
      if (OprogStatus.useNNL && OprogStatus.optnnl)
	{
	  com[0] = 0.0;
	  com[1] = 0.0;
	  com[2] = 0.0;
	  for (sp = 0; sp < typesArr[pt].nspots; sp++) 
	    {	
	      for (kk=0; kk < 3; kk++)
		com[kk] += typesArr[pt].spots[sp].x[kk];
	    }
	  for (kk=0; kk < 3; kk++)
	    com[kk] /= ((double)(typesArr[pt].nspots+1.0));
	  for (kk=0; kk < 3; kk++)
	    typesArr[pt].ppr[kk] = com[kk];
	  for (kk = 0; kk < 3; kk++)
	    typesArr[pt].ppsax[kk] = typesArr[pt].sax[kk]+fabs(com[kk]);
	  MD_DEBUG31(printf("typesArr.sax= %f %f %f .ppsax=%f %f %f\n", typesArr[pt].sax[0],typesArr[pt].sax[1],typesArr[pt].sax[2],
			    typesArr[pt].ppsax[0],typesArr[pt].ppsax[1],typesArr[pt].ppsax[2]));
	  //printf("pt=%d com=%f %f %f\n", pt, com[0], com[1], com[2]);
	}
      else
	{
	  for (kk = 0; kk < 3; kk++)
	    {
	      typesArr[pt].ppsax[kk] = typesArr[pt].sax[kk];
	      /* N.B. 15/06/10: inizializza anche ppr anche se in teoria non serve 
		 poich� non lo usa se optnnl=0 */
	      typesArr[pt].ppr[kk] = 0.0;
	    }
	}
#else
      for (kk = 0; kk < 3; kk++)
	typesArr[pt].ppsax[kk] = typesArr[pt].sax[kk];
#endif

      for (sp = 0; sp < typesArr[pt].nspots; sp++) 
	{
	  //norm = calc_norm(typesArr[pt].spots[sp].x);
#ifdef MD_EDHEFLEX_OPTNNL
	  if (OprogStatus.useNNL && OprogStatus.optnnl)
	    {
	      for (kk=0; kk < 3; kk++)
		v[kk] = typesArr[pt].spots[sp].sigma*0.5 + fabs(typesArr[pt].spots[sp].x[kk]-com[kk]);
	    }
	  else
	    {
	      for (kk=0; kk < 3; kk++)
		v[kk] = typesArr[pt].spots[sp].sigma*0.5 + fabs(typesArr[pt].spots[sp].x[kk]);
	    }
#else
	  for (kk=0; kk < 3; kk++)
    	    v[kk] = typesArr[pt].spots[sp].sigma*0.5 + fabs(typesArr[pt].spots[sp].x[kk]);
#endif
	  //printf("pt=%d sp=%d v=%.15G %.15G %.15G\n", pt, sp, v[0], v[1], v[2]);
	  for (kk = 0; kk < 3; kk++)
	    {
	      if (v[kk] > typesArr[pt].ppsax[kk]) 
		typesArr[pt].ppsax[kk] = v[kk];
	    }
	  //printf("pt=%d ppsax=%f %f %f\n", pt, typesArr[pt].ppsax[0],typesArr[pt].ppsax[1],typesArr[pt].ppsax[2]);
	}   
    }
}
#endif
double calc_shell(void)
{
#ifdef MD_PATCHY_HE
#ifndef EDHE_FLEX
  int aa;	
#endif
  int kk;
  double v[3], norm, delta;
#endif
  double deltamax = 0.0;
#ifdef EDHE_FLEX
  int pt, sp;
#endif
#ifdef MD_PATCHY_HE
#ifdef EDHE_FLEX
  for (pt = 0; pt < Oparams.ntypes; pt++)
    {
#ifdef MD_SPHERICAL_WALL
      if (pt == Oparams.ntypes-1)
	continue;
      if (pt == Oparams.ntypes-2)
	continue;
#endif
      for (sp = 0; sp < typesArr[pt].nspots; sp++) 
	{
	  norm = calc_norm(typesArr[pt].spots[sp].x);
	  if (norm!=0.0)
	    {
	      for (kk=0; kk < 3; kk++)
		v[kk] = (norm + typesArr[pt].spots[sp].sigma*0.5) * typesArr[pt].spots[sp].x[kk] / norm;
	    }
	  else
	    {
	      for (kk=0; kk < 3; kk++)
		v[kk] = typesArr[pt].spots[sp].sigma*0.5;
	    }	      
	  //printf("pt=%d sp=%d v=%.15G %.15G %.15G\n", pt, sp, v[0], v[1], v[2]);
	  delta = max3(v[0]-typesArr[pt].sax[0],v[1]-typesArr[pt].sax[1],v[2]-typesArr[pt].sax[2]);
	  if (( pt==0 && sp == 0) || delta > deltamax)
	    deltamax  = delta;
	  //printf("deltamax=%.15G\n", deltamax);
	}   
    }
#else
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
       if (delta > deltamax)
	 deltamax  = delta;
    }
#endif
#else
  deltamax = 0.0;
#endif
  //printf("deltamax=%.15G\n", deltamax);
  return deltamax;
}
double calc_nnl_rcut(void)
{
  double rcutA, rcutB, del;
#ifdef EDHE_FLEX
  int kk;
  double ax[3];
#endif
#if defined(MD_POLYDISP) || defined(EDHE_FLEX) 
  int i;
  double rcutMax=0.0;
#endif 
  /* nel caso in cui sono si usino le NNL � inutile aggiungere rNebrShell 
     per il calcolo di rcut. */
  if (OprogStatus.useNNL)
    {
      del = OprogStatus.rNebrShell;
    }
  else
    {
      del = 0.0;
    }
#ifdef EDHE_FLEX
  for (i = 0; i < Oparams.parnum; i++)
    {
#ifdef MD_SPHERICAL_WALL
      if (typeOfPart[i] == sphWall)
	continue;
      if (typeOfPart[i] == sphWallOuter)
	continue;
#endif
      for (kk=0; kk < 3; kk++)
	ax[kk] = typesArr[typeOfPart[i]].ppsax[kk];
      /* nel caso si tratti di un oggetto a simmetria sferica l'orientazione rimane della NNL (che � un cubo) rimane invariata
	 nel tempo per cui si pu� prendere rcut appena pi� grande del lato della NNL cubica*/
      if (is_a_sphere_NNL[i])
	rcutA = 2.0*max3(ax[0]+del,ax[1]+del,ax[2]+del);
      else
	rcutA = 2.0*sqrt(Sqr(ax[0]+del)+Sqr(ax[1]+del)+Sqr(ax[2]+del));
      if  (rcutA  > rcutMax)
	rcutMax = rcutA;
    }
  return 1.01*rcutMax;

#elif defined(MD_POLYDISP) 
  for (i = 0; i < Oparams.parnum; i++)
    {
      rcutA = 2.0*sqrt(Sqr(axaP[i]+del)+Sqr(axbP[i]+del)+Sqr(axcP[i]+del));
      if  (rcutA  > rcutMax)
	rcutMax = rcutA;
    }
  return 1.01*rcutMax;
#else
  rcutA = 2.0*sqrt(Sqr(Oparams.a[0]+del)+Sqr(Oparams.b[0]+del)+Sqr(Oparams.c[0]+del));
  rcutB = 2.0*sqrt(Sqr(Oparams.a[1]+del)+Sqr(Oparams.b[1]+del)+Sqr(Oparams.c[1]+del));
  return 1.01*max(rcutA, rcutB);
#endif
}
#ifdef MD_HE_PARALL
MPI_Datatype Particletype;
MPI_Datatype Eventtype;

void mpi_define_structs(void)
{
  MPI_Datatype type_pair[15]={MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
    MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE,
    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  int blocklen_pair[15]={2,6,18,12,8,6,6,4,2,2,2,2,2,2,18}, a;
  MPI_Aint displ_pair[15];
  MPI_Datatype type_ev[5]={MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT, MPI_INT};
  int blocklen_ev[5]={1,3,1,1,3};
  MPI_Aint displ_ev[5];
  MPI_Address(&parall_pair, &displ_pair[0]);
  MPI_Address(parall_pair.pos, &displ_pair[1]);
  MPI_Address(parall_pair.R, &displ_pair[2]);
  MPI_Address(parall_pair.vels, &displ_pair[3]);
  MPI_Address(parall_pair.axes, &displ_pair[4]);
  MPI_Address(parall_pair.cells, &displ_pair[5]);
  MPI_Address(parall_pair.lastbump, &displ_pair[6]);
  MPI_Address(parall_pair.time,  &displ_pair[7]);
  MPI_Address(parall_pair.atomTime, &displ_pair[8]);
#ifdef MD_ASYM_ITENS
  MPI_Address(parall_pair.angM,  &displ_pair[9]);
  MPI_Address(parall_pair.sintheta0, &displ_pair[10]);
  MPI_Address(parall_pair.costheta0, &displ_pair[11]);
  MPI_Address(parall_pair.phi0,      &displ_pair[12]);
  MPI_Address(parall_pair.psi0,      &displ_pair[13]);
  MPI_Address(parall_pair.RM,        &displ_pair[14]);
  for (a = 14; a >= 0; a--)
    displ_pair[a] -= displ_pair[0];
  MPI_Type_struct(15, blocklen_pair, displ_pair, type_pair, &Particletype);
#else
  for (a = 8; a >= 0; a--)
    displ_pair[a] -= displ_pair[0];
  MPI_Type_struct(9, blocklen_pair, displ_pair, type_pair, &Particletype);
#endif
  MPI_Type_commit(&Particletype);
  MPI_Address(&parall_event, &displ_ev[0]);
  MPI_Address(parall_event.rC, &displ_ev[1]);
  MPI_Address(&parall_event.a, &displ_ev[2]);
  MPI_Address(&parall_event.b, &displ_ev[3]);
#ifdef MD_PATCHY_HE
  MPI_Address(parall_event.sp, &displ_ev[4]);
  for (a = 4; a >= 0; a--)
    displ_ev[a] -= displ_ev[0];
  MPI_Type_struct(5, blocklen_ev, displ_ev, type_ev, &Eventtype);
#else
  for (a = 3; a >= 0; a--)
    displ_ev[a] -= displ_ev[0];
  MPI_Type_struct(4, blocklen_ev, displ_ev, type_ev, &Eventtype);
#endif
  MPI_Type_commit(&Eventtype);
}
void md_mpi_init(int *pargc, char***pargv)
{
  MPI_Init(pargc, pargv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs); 
  mpi_define_structs();
}
extern const int iwtagEvent, iwtagPair;
void md_mpi_finalize()
{
  int npr;
  char msgtype;
  for(npr = 1; npr < numOfProcs; npr++)
    {
      //parall_pair.p[0] = -2;
      //parall_pair.p[1] = -2;
      msgtype = 'T';
      MPI_Send(&msgtype, 1, MPI_CHAR, npr, iwtagPair, MPI_COMM_WORLD);
      //MPI_Send(&parall_pair, 1, Particletype, npr, iwtagPair, MPI_COMM_WORLD);
    }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
#endif
#ifdef MD_HE_PARALL
extern int parall_slave_get_data(parall_pair_struct *parall_pair);
extern void find_contact_parall(int na, int n, parall_event_struct *parall_event);
extern MPI_Status parall_status;
void slave_task(void)
{
  int njob, iriceve;
  char msgtype;

  if (my_rank == 0)
    return;
  while(1)
    {
      /* slaves processes here */
      for(njob = 0; njob >= 0; njob++)
	{
	  //printf("receiving new job rank=%d\n", my_rank);
	  MPI_Recv(&msgtype, 1, MPI_CHAR, 0, iwtagPair, MPI_COMM_WORLD, &parall_status);
	  if (msgtype == 'D')
	    {
	      break;
	    }
	  else if (msgtype == 'T')
	    {
	      //printf("FINE rank=%d\n", my_rank);
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Finalize();
	      exit(0);
	    }
	  MPI_Recv(&parall_pair, 1, Particletype, 0, iwtagPair, MPI_COMM_WORLD, &parall_status);
	  iriceve = parall_status.MPI_SOURCE;
	  //printf("received from %d rank=%d\n", iriceve, my_rank);
	  parall_slave_get_data(&parall_pair);
	  //printf("FINDING CONTACT TIME....rank=%d\n", my_rank);
	  find_contact_parall(parall_pair.p[0], parall_pair.p[1], &parall_event);
	  //printf("<<<<DONE my_rank=%d  i=%d j=%d\n", my_rank, parall_pair.p[0],
	  //     parall_pair.p[1]);
	  MPI_Send(&parall_event, 1, Eventtype, 0, iwtagEvent, MPI_COMM_WORLD);
	  //printf("mmmmmmaaa rank=%d\n", my_rank);
	  /* predict collision here */
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
}
#endif
#ifdef MD_CALC_DPP
extern void store_last_u(int i);
#endif
#ifdef EDHE_FLEX
int get_max_nbonds(void)
{
  int pt1, pt2, i, maxijbonds=0;
  int ni, type1, type2, a, maxpbonds=0;
  int *intersI;
  for (pt1=0; pt1 < Oparams.ntypes; pt1++)
    {
      for (pt2=pt1; pt2 < Oparams.ntypes; pt2++)
	{
	  type1 = pt1;
	  type2 = pt2;
	  a=0;
	  for (ni=0; ni < Oparams.ninters; ni++)
	    {
	      if ((intersArr[ni].type1 == type1 && intersArr[ni].type2 == type2) ||
		  (intersArr[ni].type1 == type2 && intersArr[ni].type2 == type1))
		{
		  a+=2;
		}	
	    }
	  if (a > maxpbonds)
	    maxpbonds = a;
	}
    }
  intersI = malloc(sizeof(int)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    intersI[i] = 0;
  for (ni=0; ni < Oparams.nintersIJ; ni++)
    {
      (intersI[intersArrIJ[ni].i])+=2; 
      (intersI[intersArrIJ[ni].j])+=2;
    }
  maxijbonds=0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (intersI[i] > maxijbonds)
	maxijbonds = intersI[i];
    }
  free(intersI);
  return maxpbonds+maxijbonds;
}
double eval_max_dist_for_spots(int pt)
{
  int ns;
  double dist, distMax=0.0;
  for (ns=0; ns < typesArr[pt].nspots; ns++)
    {
      dist = calc_norm(typesArr[pt].spots[ns].x) + typesArr[pt].spots[ns].sigma*0.5;
      if (dist > distMax)
	distMax = dist;
    }
  return distMax;
}
void assignPartTypes(void)
{
  int pt, i, ini, end, nt=0;
  typeOfPart = malloc(sizeof(int)*Oparams.parnum);
  ini = 0;
  end = typeNP[0];
  for (pt=0; pt < Oparams.ntypes; pt++)
    {
      for (i=ini; i < end; i++)
	{
	  typeOfPart[i] = nt; 
	}
      ini += typeNP[nt];
      end += typeNP[nt+1];
      nt++;
    }
}
void check_conf(void)
{
  int i, pt;
  int *typeNPL;
  typeNPL = malloc(sizeof(int)*Oparams.ntypes);
  for (i=0; i < Oparams.ntypes; i++)
    {
      typeNPL[i] = 0;
    }
  for (i=0; i < Oparams.parnum; i++)
    {
      typeNPL[typeOfPart[i]]++;
    }
  for (pt = 0; pt < Oparams.ntypes; pt++)
    {
      if (typeNPL[pt] != typeNP[pt])
	{
	  printf("WARNING: the number of particles of type %d (%d) is different from what expected (%d)!\n",
		 pt, typeNPL[pt], typeNP[pt]);
	}
    }
  free(typeNPL);
}
#endif
#ifdef EDHE_FLEX
extern int get_dof_flex(int filter);
extern int all_spots_in_CoM(int pt);
#ifdef MD_SUPERELLIPSOID
extern int is_superellipse(int i);
#endif
void find_spheres_NNL(void)
{
  int i, k1, k2, pt;
  for (i = 0; i < Oparams.parnum; i++)
    {
      is_a_sphere_NNL[i] = 1;
      pt = typeOfPart[i];
      if (!(typesArr[pt].sax[0] == typesArr[pt].sax[1] && 
	    typesArr[pt].sax[1] == typesArr[pt].sax[2] )) 
	{
	  is_a_sphere_NNL[i] = 0;
	  continue;
	}
#ifdef MD_SUPERELLIPSOID
      if (is_superellipse(i))
	{
	  is_a_sphere_NNL[i] = 0;
	  continue;
	}
#endif
      if (typesArr[pt].nspots > 0 && !all_spots_in_CoM(pt))
	{
	  is_a_sphere_NNL[i] = 0;
	  continue;
	}
      /* l'orientazione non cambia durante la simulazione
	 se si tratta di oggetti a simmetria sferica */
      for (k1 = 0; k1 < 3 ; k1++)
	for (k2 = 0; k2 < 3 ; k2++)
	  R[i][k1][k2] = (k1==k2)?1.0:0.0;
    }
}
#endif
#ifdef EDHE_FLEX
void par2saveArr(void);
extern void saveFullStore(char* fname);
void set_angmom_to_zero(int i)
{
  Mx[i]=My[i]=Mz[i]=0.0;
  wx[i]=wy[i]=wz[i]=0.0;
}
#endif
#if defined(EDHE_FLEX) || defined(MD_PATCHY_HE)
extern void find_bonds_one(int i);
extern void find_bonds_one_NLL(int i);

void find_bonds_flex_all(void)
{
  int i;
  printf("===========================>QUI\n");
  for (i=0; i < Oparams.parnum; i++)
    find_bonds_one(i); 
}
void find_bonds_flex_NNL(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    find_bonds_one_NLL(i); 
}
void find_bonds_flex(void)
{
#ifdef MD_SPHERICAL_WALL
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (typeOfPart[i]==1 || typeOfPart[i]==2)
	  {
	    if (typeOfPart[i]==1)
	      add_bond(i, sphWall, 1, 1);
	    add_bond(i, sphWallOuter, 1, 1);
	  }
    } 
#endif
  if (OprogStatus.useNNL)
    {
      find_bonds_flex_NNL();
    }
  else
    {
      find_bonds_flex_all();
    }
}
void find_bonds(void)
{
#ifndef EDHE_FLEX
  double dists[MD_PBONDS];
#endif
  int i, j, NPB, nn, aa, bb;
  double drx, dry, drz, dist, shift[3];
  int amin, bmin;
  //printf("INIZIO CERCO BOND\n");
#ifdef EDHE_FLEX
  find_bonds_flex();
  return;
#endif
  for ( i = 0; i < Oparams.parnum-1; i++)
    for ( j = i+1; j < Oparams.parnum; j++)
      {
#ifdef MD_SPHERICAL_WALL
	/* treat bonds with spherical walls separately */
	if (((j==sphWall||j==sphWallOuter) && typeOfPart[i]==1) || (typeOfPart[i]==2 && j==sphWallOuter))
	  {
	    //printf("qui\n");
	    //printf("SPHWALL i=%d j=%d numbonds[1841]=%d\n", i, j, numbonds[1841]);
	    
	    add_bond(i, j, 1, 1);
	    //printf("boh\n");
	    //add_bond(j, i, 1, 1);
	    continue;
	  }
#endif
	/* l'interazione sticky � solo fra fra A e B! */
#ifndef EDHE_FLEX
	if (!((i < Oparams.parnumA && j >= Oparams.parnumA)|| 
	      (i >= Oparams.parnumA && j < Oparams.parnumA)))
	  continue;
#endif
	drx = rx[i] - rx[j];
#ifdef MD_LXYZ
	shift[0] = L[0]*rint(drx/L[0]);
#else
	shift[0] = L*rint(drx/L);
#endif
	dry = ry[i] - ry[j];
#ifdef MD_LXYZ
	shift[1] = L[1]*rint(dry/L[1]);
#else
	shift[1] = L*rint(dry/L);
#endif
	drz = rz[i] - rz[j]; 
#ifdef MD_EDHEFLEX_WALL
	if (!OprogStatus.hardwall)
	  {
#ifdef MD_LXYZ
	    shift[2] = L[2]*rint(drz/L[2]);
#else
    	    shift[2] = L*rint(drz/L);
#endif
	  }
	else
	  shift[2] = 0.0;
#else
#ifdef MD_LXYZ
	shift[2] = L[2]*rint(drz/L[2]);
#else
	shift[2] = L*rint(drz/L);
#endif
#endif
	assign_bond_mapping(i, j);

	dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
	NPB = get_num_pbonds(i, j);
#if 0
	//if (NPB==0)
	  printf("i=%d j=%d NPB=%d\n",i, j, NPB);
#endif
	for (nn=0; nn < NPB; nn++)
	  {
	    //printf("Oparams.parnum=%d i=%d typei=%d j=%d typej=%d nn=%d dist=%.15G\n", Oparams.parnum, i, typeOfPart[i], j, typeOfPart[j], nn, dists[nn]);
	    //if (mapbondsa[nn] >= 15 && mapbondsb[nn] >= 15)
		//  printf("(%d,%d)-(%d,%d): %.15G\n", i, mapbondsa[nn]-1, j, mapbondsb[nn]-1, dists[nn]);
    	    if (dists[nn] < 0.0)
	      {
#if 0
		if (j==1001)
		  printf("Oparams.parnum=%d i=%d typei=%d j=%d typej=%d nn=%d dist=%.15G\n", Oparams.parnum, i, typeOfPart[i], j, typeOfPart[j], nn, dists[nn]);
#endif
		//printf("adding bond i=%d j=%d aa=%d bb=%d\n", i, j, aa, bb);
		aa = mapbondsa[nn];
		bb = mapbondsb[nn];
		//printf("adding bond (%d,%d)-(%d,%d)\n", i, aa, j, bb);
		add_bond(i, j, aa, bb);
		add_bond(j, i, bb, aa);
	      }
	  }
      }

  //printf("FINE CERCO BOND\n");
}
#endif
#ifdef EDHE_FLEX
int same_position(double x1[3], double x2[3])
{
  int kk;
  for (kk=0; kk < 3; kk++)
    if (x1[kk]!=x2[kk])
      return 0;
  return 1;
}
void find_conciding_spots(void)
{
  int nt, ns1, ns2;
  for (nt=0; nt < Oparams.ntypes; nt++)
    {
      /* 10/06/2010: qui comunque inizializzao anche gli spot utilizzati per le NNL
      ossia quelli che vanno da typesArr[nt].nspots a typesArr[nt].nspots + MD_SPNNL_NUMSP */
      for (ns1=0; ns1 < typesArr[nt].nspots + MD_SPNNL_NUMSP; ns1++)
	typesArr[nt].spots[ns1].same=ns1;
      for (ns1=0; ns1 < typesArr[nt].nspots; ns1++)
	{
	  /* se spots[ns1].same=ns1 vuol dire che tale spots non coincide con nessun altro */
	  if (typesArr[nt].spots[ns1].same != ns1)
	    continue;
	  for (ns2=ns1+1; ns2 < typesArr[nt].nspots; ns2++)
	    {
	      if (same_position(typesArr[nt].spots[ns1].x,typesArr[nt].spots[ns2].x) &&
		  typesArr[nt].spots[ns2].same == ns2)
		{
		  typesArr[nt].spots[ns2].same = ns1;
		}
	    }
	}
    }

}
#endif
#ifdef EDHE_FLEX
void boh(void)
{
  double dl, tini=0.0, t1=0.0, tl; 
  double shift[3]={0.0,0.0,0.0};
  double dist[1000];
  int amin, bmin; 
  tl = tini;
  dl = calcDistNegSP(tl, t1, 118, 58, shift, &amin, &bmin, dists, -1);
  printf("1)t1=%.15G t=%.15G d=%.15G\n", t1, tl, dl);
  tl = tini+1E-20;
  dl = calcDistNegSP(tl, t1, 118, 58, shift, &amin, &bmin, dists, -1);
  printf("2)t1=%.15G t=%.15G d=%.15G\n", t1, tl, dl);

  exit(-1);
}
#endif
#ifdef MD_SPHERICAL_WALL
void allocBondsSphWall(void)
{
  int i;
  for (i=0; i < Oparams.parnum; i++)
    {
      /* gli ultimi due tipi devono essere i "muri" sferici */
      if (typeOfPart[i]==Oparams.ntypes-1 || typeOfPart[i]==Oparams.ntypes-2)
	{
#ifdef MD_LL_BONDS
	  /* NOTA 21/04/2010: il 3 l'ho messo per tener conto del fatto che nel caso ad esempio 
	  con l'interazione SW i legami possono essere anche due per particella. */
	  bonds[i] = malloc(sizeof(long long int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#else
	  bonds[i] = malloc(sizeof(int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#endif
	  //break;
	}
    }
}
#endif
#ifdef MD_GHOST_IGG
extern int get_rabbit_bonds(int ifebA, int tA, int ifebB, int tB);
#endif
#ifdef MD_GHOST_IGG
void init_ghostArr(void)
{
  int i, iggStatus, nigg;
  nigg=0;
  iggStatus = 0;
  
  if (ghostInfoArr)
    return;/* if here it means that ghostInfoArr has been read from file */

  ghostInfoArr = malloc(sizeof(ghostInfo)*Oparams.parnum);
  /* here I assume that Igg HE are first Nigg */
  for (i=0; i < Oparams.parnum; i++)
    {
      if (typeOfPart[i] > 3)
	{
	  ghostInfoArr[i].iggnum = -1;
	  ghostInfoArr[i].ghost_status = -1; 
	  continue;
	}

      /* if ghostsim == 2 accept only transition 3->1, i.e. all antibody
	 will be in status 3 and only transition from 3 to 1 will be accepted.-
	 if ghostsim==2 all igg will be ghost and the will remain forever in state 3 (=ghost) */
      if (Oparams.ghostsim == 2 || Oparams.ghostsim == 3)
	{
	  iggStatus = 3;
	}
      else if (typeOfPart[i] == 0)
	{
	  iggStatus = (get_rabbit_bonds(i, 0, i+1, 1) > 0)?2:3;
	}
      /* set all IgG in state 3, so that even if IgG overlap 
     	 it is ok */
      ghostInfoArr[i].iggnum = nigg;
      ghostInfoArr[i].ghost_status = iggStatus; 
      //printf("i=%d status=%d nigg=%d\n", i, iggStatus, nigg);
      if (typeOfPart[i]==3)
	nigg++;
    }
}
#endif
#ifdef MD_PROTEIN_DESIGN
extern char nativeConf[512];
double *rxNat, *ryNat, *rzNat, *RNat[3][3];
#ifdef MD_LL_BONDS
long long int **bondsNat;
int *numbondsNat;
#else
int *numbondsNat, **bondsNat;
#endif
char dummy2[1024];
void read_native_conf(void)
{
  int i, j, parnum;
  FILE *fs;
  char sep[256];
  // increase OprogStatus.maxbonds in usrInitBef if needed!!!
  //printf("nativeConf:%s\n", nativeConf);
  
  fs = fopen(OprogStatus.nativeConf, "r");
  fscanf(fs, "%d ", &parnum);
  //printf("parnum=%d OprogStatus.maxbonds=%d\n", parnum, OprogStatus.maxbonds);
  rxNat = malloc(sizeof(double)*parnum);
  ryNat = malloc(sizeof(double)*parnum);
  rzNat = malloc(sizeof(double)*parnum);
  numbondsNat = malloc(sizeof(int)*parnum);
#ifdef MD_LL_BONDS
  bondsNat = AllocMatLLI(parnum, OprogStatus.maxbonds);
#else
  bondsNat = AllocMatI(parnum, OprogStatus.maxbonds);
#endif
  for (i = 0; i < parnum; i++)
    {
      fscanf(fs, "%d ", &numbondsNat[i]);
      for (j = 0; j < numbondsNat[i]; j++)
	{
#ifdef MD_LL_BONDS
	  fscanf(fs, "%lld ", &bondsNat[i][j]);
#else	      
	  fscanf(fs, "%d ", &bondsNat[i][j]);
#endif	    
	}
    }
  for (i=0; i < parnum; i++)
    fscanf(fs, "%lf %lf %lf %[^\n]", &(rxNat[i]), &(rxNat[i]), &(ryNat[i]), dummy2);
  fclose(fs);
}
int boundNat(int na, int n)
{
  int a;
  if (abs(na-n)==1)
    return 0;
  for (a=0; a < numbondsNat[na]; a++)
    {
#ifdef MD_LL_BONDS
    if (bondsNat[na][a] / (((long long int)NA)*NA) == n)
      return 1;
#else
    if (bondsNat[na][a] / (NANA) == n)
      return 1;
#endif
    }
  return 0;
}
double calc_order_param_native(void)
{
  int i, j;
  double RMSD, cc, DSQ, DSQN;
  RMSD=0.0;
  cc=0;
  for (i=0; i < Oparams.parnum; i++)
    for (j=i+1; j < Oparams.parnum; j++)
      {
	if (boundNat(i, j))
	  {
	    DSQ=0.0;
	    DSQ += Sqr((rx[i]-rx[j]));
	    DSQ += Sqr((ry[i]-ry[j]));
	    DSQ += Sqr((rz[i]-rz[j]));
	    DSQ = sqrt(DSQ);

	    DSQN=0.0;
	    DSQN += Sqr((rxNat[i]-rxNat[j]));
	    DSQN += Sqr((ryNat[i]-ryNat[j]));
	    DSQN += Sqr((rzNat[i]-rzNat[j]));
	    DSQN = sqrt(DSQN);

	    RMSD += Sqr(DSQ-DSQN);
	    //printf("qui!i=%d j=%d RMSDP=%.15G\n", i, j, Sqr((rxNat[i]-rxNat[j])-(rx[i]-rx[j])));
	    cc=cc+1.0;
	  }
      }
  RMSD = sqrt(RMSD/cc);
  //printf("RMSD=%.15G\n", RMSD);
  return RMSD;
}
#endif
#ifdef MD_ASYM_ITENS
extern int isSymItens(int i);
#endif
#ifdef MD_RABBIT
int nbins_rhoz;
double *rhoz;
extern double max3(double a, double b, double c);
#endif
#ifdef MD_SUPERELLIPSOID
void init_gauleg_weights(void);
double *SQvolPrefact;
double calc_SQ_volprefact(int);
#endif
#ifdef EDHE_FLEX
void buildSPNNL_spots_growth(int i)
{
  int pt, nsp, kk, k, vert[3];
  double signArr[MD_SPNNL_NUMSP][3] = {{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{-1,-1,1},{1,-1,-1},{-1,1,-1},{-1,-1,-1}};
  pt = typeOfPart[i];
  vert[0] = axa[i];
  vert[1] = axb[i];
  vert[2] = axc[i];
  nsp = typesArr[pt].nspots;
  for (k = 0; k < 8; k++)
    {
      /* NOTA 16/04/2010: in tal modo gli 8 spot coincidono con i vertici di un parallelepipedo
	 che hai i lati esattamente uguali al doppio dei semiassi della SQ */
      for (kk = 0; kk < 3; kk++)
	typesArr[pt].spots[nsp+k].x[kk] = vert[kk]*signArr[k][kk]; 
      typesArr[pt].spots[nsp+k].same=nsp+k;
      typesArr[pt].spots[nsp+k].sigma = MD_SPNNL_SIGMA; 
    }
}
void buildSPNNL_spots(void)
{
  int pt, nsp, kk, k, vert[3];
  double signArr[MD_SPNNL_NUMSP][3] = {{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{-1,-1,1},{1,-1,-1},{-1,1,-1},{-1,-1,-1}};
  for (pt=0; pt < Oparams.ntypes; pt++)
    {
      for (kk=0; kk < 3; kk++) 
	vert[kk] = typesArr[pt].sax[kk];
      nsp = typesArr[pt].nspots;
      for (k = 0; k < 8; k++)
	{
	/* NOTA 16/04/2010: in tal modo gli 8 spot coincidono con i vertici di un parallelepipedo
	   che hai i lati esattamente uguali al doppio dei semiassi della SQ */
	  for (kk = 0; kk < 3; kk++)
	    typesArr[pt].spots[nsp+k].x[kk] = vert[kk]*signArr[k][kk]; 
	  typesArr[pt].spots[nsp+k].same=nsp+k;
	  typesArr[pt].spots[nsp+k].sigma = MD_SPNNL_SIGMA; 
	}
    }
}
#endif
#ifdef MD_CALENDAR_HYBRID
extern int *linearLists;
void estimate_HQ_params(double phi)
{
  /* linear approximation for linear dependence on N */
  double scalevsNfact[4]={0.0877,0.04862,0.9408,7.532}; 
  double nlistsvsNfact[4]={48.67,97.34,244.7,564.182};
  double volfact[4] = {0.01,0.12,0.4,0.7}, msc, mnl, qsc, qnl, scf, nlf;
  int MAXINT=500000000;
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
      OprogStatus.scaleHQ = (int) (scf*Oparams.parnum);

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
}
#if 1
void rebuild_linked_list(void);
void rebuildCalendar(void);
void adjust_HQ_params(void)
{
  int targetNE = 15, del=5;
  int k, i, NAVG=10, k1, k2;
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
      ((double)overevHQavg) <= OprogStatus.overthrHQ)
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
      /* OprogStatus.overthrHQ � il numero massimo di eventi nella overflow list che viene considerato
	 accettabile */ 
      if (((double)overevHQavg) > OprogStatus.overthrHQ) 
	OprogStatus.nlistsHQ *= GOLD;
    }
  printf("ADJUSTING HQ PARAMS: scaleHQ=%G nlistsHQ=%d\n", OprogStatus.scaleHQ, OprogStatus.nlistsHQ);
  free(linearLists);
  linearLists = malloc(sizeof(int)*(OprogStatus.nlistsHQ+1));
  UpdateSystem();
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  rebuild_linked_list();

  if (OprogStatus.useNNL)
    rebuildNNL();
  rebuildCalendar();
  if (OprogStatus.intervalSum > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
#ifdef MD_BIG_DT
  if (OprogStatus.bigDt > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT + 11,OprogStatus.bigDt);
#endif
  if (OprogStatus.adjustHQ == 2)
    ENDSIM=1;
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
#ifdef MD_CALC_DPP
  OprogStatus.len = sizeof(double)*22*Oparams.parnum;
#else
  OprogStatus.len = sizeof(double)*10*Oparams.parnum;
#endif
  //printf("DYNAMIC ALLOCATION of %d bytes\n", OprogStatus.len);
  OprogStatus.ptr = malloc(OprogStatus.len);
  last_ptr = OprogStatus.ptr;
#ifdef MD_CALC_DPP
  OprogStatus.sumdx = (double*)last_ptr;
  OprogStatus.sumdy = OprogStatus.sumdx + np;
  OprogStatus.sumdz = OprogStatus.sumdy + np;
  OprogStatus.lastu1x = OprogStatus.sumdz + np;
  OprogStatus.lastu1y = OprogStatus.lastu1x + np;
  OprogStatus.lastu1z = OprogStatus.lastu1y + np;
  OprogStatus.lastu2x = OprogStatus.lastu1z + np;
  OprogStatus.lastu2y = OprogStatus.lastu2x + np;
  OprogStatus.lastu2z = OprogStatus.lastu2y + np;
  OprogStatus.lastu3x = OprogStatus.lastu2z + np;
  OprogStatus.lastu3y = OprogStatus.lastu3x + np;
  OprogStatus.lastu3z = OprogStatus.lastu3y + np;
  last_ptr = (void*) (OprogStatus.lastu3z + np);
#endif
  OprogStatus.sumox = (double*)last_ptr;
  OprogStatus.sumoy = OprogStatus.sumox + np;
  OprogStatus.sumoz = OprogStatus.sumoy + np;
  OprogStatus.lastcolltime = OprogStatus.sumoz + np;
  OprogStatus.rxCMi = OprogStatus.lastcolltime + np;
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
	opro_ascii[k].ptr = OprogStatus.DR[0];
      if (!strcmp(opro_ascii[k].parName,"sumox"))
	opro_ascii[k].ptr = OprogStatus.sumox;
      if (!strcmp(opro_ascii[k].parName,"sumoy"))
	opro_ascii[k].ptr = OprogStatus.sumoy;
      if (!strcmp(opro_ascii[k].parName,"sumoz"))
	opro_ascii[k].ptr = OprogStatus.sumoz;
      if (!strcmp(opro_ascii[k].parName,"lastcolltime"))
	opro_ascii[k].ptr = OprogStatus.lastcolltime;
#ifdef MD_CALC_DPP
      if (!strcmp(opro_ascii[k].parName,"sumdx"))
	opro_ascii[k].ptr = OprogStatus.sumdx;
      if (!strcmp(opro_ascii[k].parName,"sumdy"))
	opro_ascii[k].ptr = OprogStatus.sumdy;
      if (!strcmp(opro_ascii[k].parName,"sumdz"))
	opro_ascii[k].ptr = OprogStatus.sumdz;
      if (!strcmp(opro_ascii[k].parName,"lastu1x"))
	opro_ascii[k].ptr = OprogStatus.lastu1x;
      if (!strcmp(opro_ascii[k].parName,"lastu1y"))
	opro_ascii[k].ptr = OprogStatus.lastu1y;
      if (!strcmp(opro_ascii[k].parName,"lastu1z"))
	opro_ascii[k].ptr = OprogStatus.lastu1z;
      if (!strcmp(opro_ascii[k].parName,"lastu2x"))
	opro_ascii[k].ptr = OprogStatus.lastu2x;
      if (!strcmp(opro_ascii[k].parName,"lastu2y"))
	opro_ascii[k].ptr = OprogStatus.lastu2y;
      if (!strcmp(opro_ascii[k].parName,"lastu2z"))
	opro_ascii[k].ptr = OprogStatus.lastu2z;
      if (!strcmp(opro_ascii[k].parName,"lastu3x"))
	opro_ascii[k].ptr = OprogStatus.lastu3x;
      if (!strcmp(opro_ascii[k].parName,"lastu3y"))
	opro_ascii[k].ptr = OprogStatus.lastu3y;
      if (!strcmp(opro_ascii[k].parName,"lastu3z"))
	opro_ascii[k].ptr = OprogStatus.lastu3z;
#endif

      k++;
    }
  while (strcmp(opro_ascii[k].parName,""));
}
#endif
#ifdef MD_MULTIPLE_LL
extern void get_types_from_nl(int nl, int *t1, int *t2);
#endif
#define MD_MATRIX_CONTIGOUS
#if defined(MD_MULTIPLE_LL) && defined(MD_OPT_MULTLL)
extern int may_interact_all_type(int t1, int t2);
#endif
#ifdef MD_CALC_VBONDING
/* ============================ >>> ranf <<< =============================== */
double ranf_vb(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return random() / ( (double) RAND_MAX );
  //return random() / (2**31-1); 
}
double fons(double theta, double alpha)
{
  double pi;
  pi = acos(0.0)*2.0;
  return cosh(alpha*cos(theta))*alpha/(4.0*pi*sinh(alpha));
}
/* return an angle theta sampled from an Onsager angular distribution */
double theta_onsager(double alpha)
{
  /* sample orientation from an Onsager trial function (see Odijk macromol. (1986) )
     using rejection method */
  /* the comparison function g(theta) is just g(theta)=1 */ 
  double pi, y, f, theta, dtheta;
  pi = acos(0.0)*2.0;
  
  do 
    {
      /* uniform theta between 0 and pi */
      theta = pi*ranf_vb();
      /* uniform y between 0 and 1 */
      y = 1.01*fons(0.0,alpha)*ranf_vb();
      f = fons(theta,alpha);
      //printf("theta=%f y=%f\n", theta, y);
    }
  while (y >= f);
  return theta;
}
double *distro;
void orient_onsager(double *omx, double *omy, double* omz, double alpha)
{
  double thons;
  double pi, phi, verso;

  pi = acos(0.0)*2.0;
  /* random angle from onsager distribution */
  thons = theta_onsager(alpha);
  //printf("thos=%f\n", thons);
  distro[(int) (thons/(pi/1000.0))] += 1.0;
  phi = 2.0*pi*ranf_vb();
  //verso = (ranf_vb()<0.5)?1:-1;
  verso=1;
  *omx = verso*sin(thons)*cos(phi);
  *omy = verso*sin(thons)*sin(phi);
  *omz = verso*cos(thons); 
  //printf("norma=%f\n", sqrt(Sqr(*omx)+Sqr(*omy)+Sqr(*omz)));
}
void orient(double *omx, double *omy, double* omz)
{
  int i;
  //double inert;                 /* momentum of inertia of the molecule */
  //double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, osq, norm;
  
  //Mtot = m; /* total mass of molecule */

  //inert = I; /* momentum of inertia */
 
  //mean = 3.0*temp / inert;

  xisq = 1.0;

  while (xisq >= 1.0)
    {
      xi1  = ranf_vb() * 2.0 - 1.0;
      xi2  = ranf_vb() * 2.0 - 1.0;
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

  *omx = ox;
  *omy = oy;
  *omz = oz; 
#if 0
  /* Choose the magnitude of the angular velocity
NOTE: consider that it is an exponential distribution 
(i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/

  osq   = - mean * log(ranf());
  o     = sqrt(fabs(osq));
  ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      *wx = ox;
      *wy = oy;
      *wz = oz;
#endif 
}

int type;
void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;

  /* first row vector */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 1; u[1] = 1; u[2] = 1;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = -1; u[1] = -1; u[2] = 1;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[0][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[0][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
  if (typesArr[0].nspots==3 && type==0)
    {
      for (k=0; k < 3 ; k++)
	u[k] = R[1][k];
      vectProdVec(R[0], u, up);
      /* rotate randomly second axis */
      angle=4.0*acos(0.0)*ranf_vb();
      xx = cos(angle);
      yy = sin(angle);
      for (k=0; k < 3 ; k++)
	R[1][k] = u[k]*xx + up[k]*yy;
      //printf("calc_norm(R[1])=%.15G\n", calc_norm(R[1]));
    }
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];

  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
}
void set_semiaxes_vb(double fx, double fy, double fz)
{
  int ii;
  for (ii=0; ii < Oparams.parnum; ii++)
    {
      nebrTab[ii].axa = fx;//typesArr[typeOfPart[ii]].sax[0]+typesArr[typeOfPart[ii]].spots[0].sigma*0.5;
      nebrTab[ii].axb = fy;//typesArr[typeOfPart[ii]].sax[1];
      nebrTab[ii].axc = fz;//typesArr[typeOfPart[ii]].sax[2];
    }
}
extern double calcDistNegNNLoverlapPlane(double t, double t1, int i, int j, double shift[3]);
extern double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
extern int calcdist_retcheck;
extern int are_spheres(int i, int j);
extern void readAsciiPars(FILE* pfs, struct pascii strutt[]);
extern void AllocCoord(int size, COORD_TYPE** pointer, ...);
extern void Newsimul(char *argom);
extern void usrInitAft(void);
extern int iniCorFormat;
void initsa(double sax[3], double p[3])
{
  FILE* fs;
  int k1, ii;
  char fn[256];
  iniCorFormat = 1; /* 1 = ascii */
	 
  //usrInitBef();
  strcpy(fn,"./sq.par");
  //printf("qui\n");
  Newsimul(fn);
  usrInitAft();
  OprogStatus.optnnl = 0;
  assign_bond_mapping(0,1);
  nebrTab = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
  for (ii= 0; ii < Oparams.parnum; ii++)
    nebrTab[ii].R = matrix(3,3);
#if 0
  Oparams.parnum=2;

  typeOfPart = malloc(sizeof(int)*Oparams.parnum);
  typeOfPart[0]=typeOfPart[1]=0;
  typesArr = malloc(sizeof(partType)*Oparams.ntypes);
  typeNP= malloc(sizeof(int));
  typeNP[0]=2;
#endif
  for (k1=0; k1 < 3; k1++)
    {
      typesArr[0].sax[k1]=sax[k1];
      typesArr[1].sax[k1]=sax[k1];
      typesArr[0].n[k1]=p[k1];
      typesArr[1].n[k1]=p[k1];
      typesArr[0].I[k1]=1.0;
      typesArr[1].I[k1]=1.0;
    }
  typesArr[0].m=1.0;
  typesArr[1].m=1.0;
  typesArr[0].ignoreCore=0;
  typesArr[1].ignoreCore=0;
  typesArr[0].nspots=0;
  typesArr[1].nspots=0;
  typesArr[0].nhardobjs=0;
  typesArr[1].nhardobjs=0;

  Oparams.ninters=0;
  Oparams.nintersIJ=0;
  OprogStatus.maxbonds=0;
  numbonds= (int *) malloc(Oparams.parnum*sizeof(int));
  numbonds[0] = numbonds[1] = 0;

  AllocCoord(sizeof(double)*Oparams.parnum, ALLOC_LIST, NULL);
  //build_parallelepipeds();
}
void initsa_(double sax[3], double p[3])
{
  initsa(sax,p);
}
extern const double saxfactMC[3];
double calcdistsa(double ra[3], double rb[3], double uxa[3], double uxb[3], double *Lx, double *Ly, double *Lz,int *errchk)
{
  /* N.B. u1, u2 ed u3 sono i vettori del sistema di riferimento solidale con il corpo rigido 
     espressi nel riferminto del laboratorio */
  int k1, k2;
  double nn, vecg[8], vecgNeg[8], shift[3], Rla[3][3], Rlb[3][3];
  double d, r1[3], r2[3], alpha, d0;
  L[0]=*Lx;
  L[1]=*Ly;
  L[2]=*Lz;

  //printf("vec1=%f %f %f vec2=%f %f %f\n", ra[0], ra[1], ra[2], rb[0], rb[1], rb[2]);
   /* set positions and orientations */
  rx[0] = nebrTab[0].r[0] = ra[0];
  ry[0] = nebrTab[0].r[1] = ra[1];
  rz[0] = nebrTab[0].r[2] = ra[2];
  rx[1] = nebrTab[1].r[0] = rb[0];
  ry[1] = nebrTab[1].r[0] = rb[1];
  rz[1] = nebrTab[1].r[0] = rb[2];
  shift[0] = L[0]*rint((rx[0]-rx[1])/L[0]);
  shift[1] = L[1]*rint((ry[0]-ry[1])/L[1]);
  shift[2] = L[2]*rint((rz[0]-rz[1])/L[2]);

  nn = calc_norm(uxa);
  for (k1=0; k1 < 3; k1++)
    uxa[k1] /= nn;

  nn = calc_norm(uxb);
  for (k1=0; k1 < 3; k1++)
    uxb[k1] /= nn;
  type=-1;
  versor_to_R(uxa[0], uxa[1], uxa[2], Rla);
  versor_to_R(uxb[0], uxb[1], uxb[2], Rlb);
  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  R[0][k1][k2] = Rla[k1][k2];
	  R[1][k1][k2] = Rlb[k1][k2];
	}
    }

  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  nebrTab[0].R[k1][k2] = R[0][k1][k2];
	  nebrTab[1].R[k1][k2] = R[1][k1][k2];
	}
    }

  set_semiaxes_vb(1.01*(typesArr[typeOfPart[0]].sax[0]),
		  1.01*(typesArr[typeOfPart[0]].sax[1]), 
		  1.01*(typesArr[typeOfPart[0]].sax[2]));
  d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, 0, 1, shift);
  /* se d0 � positiva vuol dire che i due parallelepipedi non s'intersecano */
  if (d0 > 0.0)
    {
      return 1.0;
    }
  set_semiaxes_vb(saxfactMC[0]*typesArr[typeOfPart[0]].sax[0],
		  saxfactMC[1]*typesArr[typeOfPart[0]].sax[1], 
		  saxfactMC[2]*typesArr[typeOfPart[0]].sax[2]);

  d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, 0, 1, shift);
  /* se d0 � positiva vuol dire che i due parallelepipedi non s'intersecano */
  if (d0 < 0.0)
    {
      return -1.0;
    }
  OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
  calcdist_retcheck = 0;
  d=calcDistNeg(0.0, 0.0, 0, 1, shift, r1, r2, &alpha, vecg, 1);
  *errchk = calcdist_retcheck;
  //printf("QUI d=%f\n", d);
  return d;
}

double calcdistsa_(double ra[3], double rb[3], double uxa[3], double uxb[3], double *Lx, double *Ly, double *Lz, int *errchk)
{
  return calcdistsa(ra, rb, uxa, uxb, Lx, Ly, Lz, errchk);
}
double calcDistNeg_vb(int i, int j, double shift[3])
{
  double vecg[8], vecgNeg[8];
  double d, r1[3], r2[3], alpha, d0;
#ifdef MD_LXYZ
  shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
  shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
  shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#else
  shift[0] = L*rint((rx[i]-rx[j])/L);
  shift[1] = L*rint((ry[i]-ry[j])/L);
  shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
  //printf("semiax=%f %f %f\n", typesArr[typeOfPart[0]].sax[0],typesArr[typeOfPart[0]].sax[1], typesArr[typeOfPart[0]].sax[2]);
  if (!are_spheres(0,1))
    {
#if 1
      if (type==0||type==2||type==5)
	set_semiaxes_vb(1.01*(typesArr[typeOfPart[0]].sax[0]+typesArr[typeOfPart[0]].spots[0].sigma),
	  		1.01*(typesArr[typeOfPart[0]].sax[1]), 
	    		1.01*(typesArr[typeOfPart[0]].sax[2]));
      else
	set_semiaxes_vb(1.01*(typesArr[typeOfPart[0]].sax[0]),
			1.01*(typesArr[typeOfPart[0]].sax[1]), 
			1.01*(typesArr[typeOfPart[0]].sax[2]));
      d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
      /* se d0 � positiva vuol dire che i due parallelepipedi non s'intersecano */
      if (d0 > 0.0)
	{
	  if (type==0 || type==2 || type==5)
	    return -1.0;
	  else
	    return 1.0;
	}
#endif
#if 1
      set_semiaxes_vb(saxfactMC[0]*typesArr[typeOfPart[0]].sax[0],
		      saxfactMC[1]*typesArr[typeOfPart[0]].sax[1], 
	    	      saxfactMC[2]*typesArr[typeOfPart[0]].sax[2]);

      d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
      /* se d0 � positiva vuol dire che i due parallelepipedi non s'intersecano */
      if (d0 < 0.0)
	{
	  return -1.0;
	}
#if 0
      set_semiaxes_vb(0.9*typesArr[typeOfPart[0]].sax[0],
		      0.3*typesArr[typeOfPart[0]].sax[1], 
		      0.95*typesArr[typeOfPart[0]].sax[2]);
      d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
      /* se d0 � positiva vuol dire che i due parallelepipedi non s'intersecano */
      if (d0 < 0.0)
	{
	  return -1.0;
	}
      set_semiaxes_vb(0.9*typesArr[typeOfPart[0]].sax[0],
		      0.95*typesArr[typeOfPart[0]].sax[1], 
		      0.3*typesArr[typeOfPart[0]].sax[2]);
      d0 = calcDistNegNNLoverlapPlane(0.0, 0.0, i, j, shift);
      /* se d0 � positiva vuol dire che i due parallelepipedi non s'intersecano */
      if (d0 < 0.0)
	{
	  return -1.0;
	}
#endif
#endif
    }
  OprogStatus.targetPhi=1.0; /* valore fittizio dato solo per far si che non esca se calcDist fallisce */
  calcdist_retcheck = 0;
  d=calcDistNeg(0.0, 0.0, i, j, shift, r1, r2, &alpha, vecg, 1);
  if (calcdist_retcheck)
    printf("NR failure\n");
  //printf("QUI d=%f\n", d);
  return d;
}

#ifdef  MD_CHAIN_SIM
double chainlen;
void save_conf(int i)
{
  char fileop2[512], fileop[512];
  sprintf(fileop2 ,"Store-%d-%d", 
	  i, 0);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
  for (i=0; i < Oparams.parnum; i++)
    {
      update_MSDrot(i);
#ifdef MD_CALC_DPP
      update_MSD(i);
      store_last_u(i);
#endif
      OprogStatus.lastcolltime[i] = Oparams.time;
    }
  R2u();
#ifdef MD_SAVE_DISTANCE
  if (mgl_mode==1)
    {
      FILE *f;
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
      double dists[MD_PBONDS], d;
      int i;
#endif
      f = fopen("distance.dat", "a");
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
      d=calcJustDistNegSP(Oparams.time, 0, 1, dists);
      fprintf(f, "%.15G %.15G ", Oparams.time + OprogStatus.refTime, calcJustDistNeg(Oparams.time, 0, 1));
#if 1
      for (i=0; i < MD_PBONDS; i++)
	{
	  fprintf(f, "%.15G ", dists[i]);
	}
#else
      fprintf(f, "%.15G ", dists[1]);
      printf("mapbonds[0]:%d mapbondsb[1]:%d\n", mapbondsa[0], mapbondsb[1]);
#endif
      fprintf(f, "\n");
#else
      fprintf(f, "%.15G %.15G\n", Oparams.time + OprogStatus.refTime, calcJustDistNeg(Oparams.time, 0, 1));
#endif
      fclose(f);
    }
#endif

  if (mgl_mode==0)
    {
      writeAsciiPars(bf, opro_ascii);
      fprintf(bf, sepStr);
      writeAsciiPars(bf, opar_ascii);
      fprintf(bf, sepStr);
    }	      
  MD_DEBUG(printf("[Store event]: %.15G JJ=%d KK=%d\n", Oparams.time, OprogStatus.JJ, OprogStatus.KK));
  //fprintf(bf, ".semiAxes: %f %f %f, %f %f %f\n",
  //	  Oparams.a[0], Oparams.b[0], Oparams.c[0],
  //  Oparams.a[1], Oparams.b[1], Oparams.c[1]);
  writeAllCor(bf, 0);
  fclose(bf);
  if (mgl_mode==0)
    {
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
    }
}
void insert_sq(int j)
{
  double rp[3], rl[3], Rl[3][3];
  double ox, oy, oz, ene;
  int k1, k2, bonded=0;
  double rangexp[2][3];
  double rangexl[2][3];
#ifndef MD_SPOT_GLOBAL_ALLOC
  double ratA[NA][3], ratB[NA][3]; 
#endif
  if (j==0)
    {
      for (k1=0; k1 < 3; k1++)
	nebrTab[0].r[k1] = 0.0;
      rx[0] = ry[0] = rz[0] = 0.0;
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    R[0][k1][k2] = (k1==k2)?1:0;
	    nebrTab[0].R[k1][k2] = R[0][k1][k2];
	  }
    }
  else
    {
      rA[0] = rx[j-1];
      rA[1] = ry[j-1];
      rA[2] = rz[j-1];
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  RtA[k1][k2] = R[j-1][k1][k2];
      /* il semi-asse x � quello lungo */
      rangexp[0][0] = 2.0*typesArr[0].sax[0]-0.15;
      rangexp[1][0] = 2.0*typesArr[0].sax[0]+0.15;
      rangexp[0][1] = -1.0;
      rangexp[1][1] = 1.0;
      rangexp[0][2] = -1.0;
      rangexp[1][2] = 1.0;
      //body2lab(rangexp[0], rangexl[0], rA, RtA);
      //body2lab(rangexp[1], rangexl[1], rA, RtA); 
      bonded=0;
      do
	{
	  /* random position insider a box */
	  rp[0] = rangexp[0][0] + ranf_vb()*(rangexp[1][0]-rangexp[0][0]);
	  rp[1] = rangexp[0][1] + ranf_vb()*(rangexp[1][1]-rangexp[0][1]);
	  rp[2] = rangexp[0][2] + ranf_vb()*(rangexp[1][2]-rangexp[0][2]);
	  body2lab(rp, rl, rA, RtA);
	  rx[j] = rl[0];
	  ry[j] = rl[1];
	  rz[j] = rl[2];
	  orient(&ox, &oy, &oz);
	  versor_to_R(ox, oy, oz, Rl);
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		nebrTab[j].R[k1][k2] = Rl[k1][k2];
		R[j][k1][k2] = Rl[k1][k2];
	      }
	  assign_bond_mapping(j-1, j);
	  d = calcDistNeg_vb(j-1, j, shift);
	  if (calcdist_retcheck==0 && d >= 0.0)
	    {
	      dist = calcDistNegSP(0.0, 0.0, j-1, j, shift, &amin, &bmin, dists, -1);
	      //printf("dist=%f d=%f nbondsFlex=%d\n", dist, d, nbondsFlex);
	      for (nn=0; nn < nbondsFlex; nn++)
		{
		  //printf("dists[%d]=%f\n", nn, dists[nn]);
		  if (dists[nn]<0.0)
		    {	
		      return;
		    }
		}
	    }
	}	
      while (!bonded);
    }
}
void calc_persist_len(int maxtrials)
{
  int i=0, j;
  OprogStatus.optnnl = 0;
  while (i < maxtrials)
    {
      for (j=0; j < Oparams.parnum; j++)
	{
	  insert_sq(j);
	}
      i++;
      save_conf();
    }
}
#endif
void calc_vbonding(void)
{
  FILE *fi;
  int k1, k2, count, ii;
  long long int maxtrials, i;
  int amin, bmin, nn;
  double Lb, alpha, Rl[3][3], ox, oy, oz, d, ene, dist;
  double totene=0.0, shift[3];
  int n=1000;
  fi = fopen("vbonding.conf", "r");
  fscanf(fi, "%lld %d ", &maxtrials, &type);
  if (type==4 || type==5)
    fscanf(fi, " %lf", &alpha);
  /* type = 0 -> calculate bonding volume 
     type = 1 -> calculate co-volume 
     type = 2 -> calc. bonding volume for fixed orientation (aligned)
     type = 3 -> calc. co-volume for fixed orientation (aligned)
     type = 4 -> covolume using Onsager trial function
     type = 5 -> bonding volume using Onsager trial function
     type = 6 -> persistence length
   */
  
  fclose(fi);
#ifdef MD_MAC
  srandomdev();
#else
  srandom((int)(time(NULL)));
#endif
#ifdef MD_CHECK_POINT
  if (type==6)
    {
      calc_persist_len(maxtrials);
      exit(-1);
    }
#endif
  OprogStatus.optnnl = 0;
  assign_bond_mapping(0,1);
  i=0;
  nebrTab = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
  for (ii= 0; ii < Oparams.parnum; ii++)
    nebrTab[ii].R = matrix(3,3);

  if (type==4 || type==5)
    {
      distro=malloc(sizeof(double)*n);
      for (i=0; i < n; i++)
	distro[i] = 0.0;
    }
  for (k1=0; k1 < 3; k1++)
    nebrTab[0].r[k1] = 0.0;
  rx[0] = ry[0] = rz[0] = 0.0;
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
     	R[0][k1][k2] = (k1==k2)?1:0;
	nebrTab[0].R[k1][k2] = R[0][k1][k2];
	if (type==2 || type==3 || type == 4)
	  {
	    nebrTab[1].R[k1][k2] = R[0][k1][k2];
	    R[1][k1][k2] = R[0][k1][k2];
	  }
      }
  while (i < maxtrials)
    {
#ifdef MD_LXYZ
      rx[1] = L[0]*(ranf_vb()-0.5);
      ry[1] = L[1]*(ranf_vb()-0.5);
      rz[1] = L[2]*(ranf_vb()-0.5);
#else
      rx[1] = L*(ranf_vb()-0.5);
      ry[1] = L*(ranf_vb()-0.5);
      rz[1] = L*(ranf_vb()-0.5);
#endif      
      if (Sqr(rx[1])+Sqr(ry[1])+Sqr(rz[1]) >= Sqr(2.0*typesArr[typeOfPart[1]].sax[0]+typesArr[typeOfPart[1]].spots[0].sigma)) 
	{
	  //printf("1)i=%d\n", i);
	  i++;
	  continue;
	}
#if 0
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  R[k1][k2] = (k1==k2)?1:0;
#endif
//printf ("alpha=%f\n", alpha);
      if (type==4 || type==5)
	{
	  orient_onsager(&ox, &oy, &oz, alpha);
	}
      else
	orient(&ox, &oy, &oz);
#if 0
      rx=10; ry=0; rz=0;
      ox=-1.0/sqrt(2); oy=1.0/sqrt(2); oz=0;
#endif
      if (type==0 || type==1 || type==4 || type==5)
	{
	  versor_to_R(ox, oy, oz, Rl);
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		nebrTab[1].R[k1][k2] = Rl[k1][k2];
		R[1][k1][k2] = Rl[k1][k2];
	      }
	}

      if (type==4 || type==5)
	{
	  orient_onsager(&ox, &oy, &oz, alpha);
	  versor_to_R(ox, oy, oz, Rl);
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
	      	R[0][k1][k2] = Rl[k1][k2];
		nebrTab[0].R[k1][k2] = Rl[k1][k2];
	      }
	}
#if 0
      if (type==2)
	{
	  double s;
	  s=ranf_vb();
	  if (s < 0.5)
	    {ox=1;oy=0;oz=0;}
	  else
	    {ox=-1;oy=0;oz=0;}
	  versor_to_R(ox, oy, oz, Rl);

	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
	      	R[1][k1][k2] = Rl[k1][k2];
		nebrTab[1].R[k1][k2] = Rl[k1][k2];
	      }
	  //print_matrix(R[1],3);
	}
#endif
      nebrTab[1].r[0] = rx[1];
      nebrTab[1].r[1] = ry[1];
      nebrTab[1].r[2] = rz[1];
      //if (i%1000==0)
	//printf("prima i=%d\n", i);
      d = calcDistNeg_vb(0, 1, shift);
      if (i%100000==0)
	printf("i=%lld d=%f calcdist_retcheck=%d\n", i, d, calcdist_retcheck);
      if (calcdist_retcheck == 0 && d < 0 && ( type==1 || type == 3 || type==4))
	totene += 1.0;
      if (calcdist_retcheck==0 && d >= 0.0 && ( type==0 || type == 2 || type == 5))
	{
	  dist = calcDistNegSP(0.0, 0.0, 0, 1, shift, &amin, &bmin, dists, -1);
	  //printf("dist=%f d=%f nbondsFlex=%d\n", dist, d, nbondsFlex);
	  ene = 0;
	  for (nn=0; nn < nbondsFlex; nn++)
	    {
	      //printf("dists[%d]=%f\n", nn, dists[nn]);
	      if (dists[nn]<0.0)
		{	
		  ene=1.0;
		  break;
		}
	    }
	  totene += ene;
	}
      i++;
	  //printf("2)i=%d\n", i);
    }
#ifdef MD_LXYZ
  Lb = L[0];
#else
  Lb = L;
#endif
  if (type==0)
    printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)i))*(Lb*Lb*Lb)/Sqr(typesArr[0].nspots), totene);
  else if (type==2)
    printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)i))*(Lb*Lb*Lb)/typesArr[0].nspots, totene);
  else if (type==5)
    printf("Vbonding=%.10f (totene=%f)\n", (totene/((double)i))*(Lb*Lb*Lb)/typesArr[0].nspots, totene);
  else if (type==1 || type == 3 || type == 4)
    {
      printf("co-volume=%.10f (totene=%f)\n", (totene/((double)i))*(Lb*Lb*Lb), totene);
      printf("%.15G\n",(totene/((double)i))*(Lb*Lb*Lb));
    }
  if (type==4||type==5)
    {
      double norm, dtheta, pi;
      FILE *F;
      pi=2.0*acos(0.0);
      norm=0.0;
      dtheta = pi/n;
      for (i=0; i < n; i++)
	norm += sin(i*dtheta)*distro[i]*2*pi*dtheta;
      for (i=0; i < n; i++)
	distro[i]/= norm;

      F=fopen("fons.dat", "w");  
      for (i=0; i < n; i++)
	fprintf(F, "%.15G %.15G\n",i*dtheta,distro[i]); 
      fclose(F);

      F=fopen("fonsExact.dat", "w");  
      for (i=0; i < n; i++)
	fprintf(F, "%.15G %.15G\n",i*dtheta,fons(i*dtheta,alpha)); 
      fclose(F);
    }
  exit(-1);
}

#endif
#if (defined(MC_SIMUL) || defined(MD_STANDALONE)) && 0
double MC_funcSQ(double x, double z, double sa[3], double ee[3])
{
  return sa[1]*pow(1.0-(pow(x/a,ee[0])+pow(z/a,ee[2])),1/ee[1]);
}
struct mboxstr 
{
  double sax[3];
  double dr[3];
  int nbox;
} **mbox;
void build_parallelepipeds(void)
{
  double sa[3], dx;
  int tt, kk, k1, k2;

  mbox = malloc(sizeof(struct mboxstr*)*Oparams.ntypes);
  nmbox = 5;
  /* 2 multibox per tipo */
  for (tt=0; tt < Oparams.ntypes; tt++)
   mbox[tt] = malloc(sizeof(struct mboxstr)*2); 

  for (tt=0; tt < Oparams.ntypes; tt++)
    {
      /* il primo set di parallelepipedi � costituito 
	 da un solo parallelepipedo */
      mbox[tt][0].nbox=1;
      for (kk = 0; kk < 3; kk++)
	{
	  mbox[tt][0].dr[kk] = 0.0;
	}

      for (kk = 0; kk < 3; kk++)
	sa[kk] = typesArr[tt].sax[kk]; 
      mbox[tt][0].sa[0] = saxfactMC[0]*sa[0];
      mbox[tt][0].sa[1] = saxfactMC[1]*sa[1];
      mbox[tt][0].sa[2] = saxfactMC[2]*sa[2]; 
      mbox[tt][1]=nmbox;
      /* secondo set di parallelepipedi: approssimazione stepwise
	 della forma della superquadrica (molto pi� accurata) */
      lastx=x=0.0;
      while (!fine)
	{
	  dx=sa[0]/nmboxmax;
	  for (ix=0; ix < 1000; ix++)
	    {
	      if (fabs(MC_funcSQ(x)-MC_funcSQ(lastx)) > OprogStatus.MC_deltambox)
		{

		}
	      x+=dx;
	    }
	}
    }
}
#endif
void usrInitAft(void)
{
  long long int maxp;
  int numll, k1old, k2old, nl, numbm;
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */

#ifdef MD_LXYZ
  int kk;
#endif
  int k1, k2;
  double phiIni;
  int Nm, i, sct, overlap, k;
  COORD_TYPE vcmx, vcmy, vcmz, MAXAX;
  int mls=0, ls;
#ifndef EDHE_FLEX
  COORD_TYPE *m;
#else
  double rcut2;
#endif
  int maxnbonds;
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
  double shift[3], drx, dry, drz, dist;
  int pt;
#ifndef EDHE_FLEX
  double dists[MD_PBONDS];
#endif
  int j, amin, bmin, nn, aa, bb, NPB;
#ifndef EDHE_FLEX
  double distSPA, distSPB;
#endif
#endif
#ifdef EDHE_FLEX
  double maxSpots;
  int maxsp=-1, first=1;
#endif
#if defined(MD_POLYDISP) 
  double stocvar;
#endif
  int a;
  /*COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;*/

#if defined(MAXPAR) && !defined(MD_DYNAMIC_OPROG)
  if (Oparams.parnum >= MAXPAR)
    {
      printf("ERROR: Too many particles, increase MAXPAR in ellipsoid.h and recompile\n");
      exit(-1);
    } 
#endif

#ifdef MD_DYNAMIC_OPROG
  OprogStatus.dyn_alloc_oprog();
#endif
#ifdef EDHE_FLEX
  Oparams.parnumA = Oparams.parnum;
#endif
  /* initialize global varibales */
  pi = 2.0 * acos(0);
#ifdef EDHE_FLEX 
#if 0
 if (OprogStatus.targetPhi > 0.0)
   {
     printf("WARNING: You can not use growth with -DEDHE_FLEX!");
     printf("Exiting...");
     exit(-1);
   } 
#else
 /* NOTA 120410: la crescita anche nel caso del codice ED_HEFLEX ora � stata 
    implementata */
 if (OprogStatus.targetPhi > 0.0)
   {
     printf("WARNING: Growth simulation, spots interactions will be disabled\n");
   }
#ifdef MD_SUPERELLIPSOID
 if (OprogStatus.useNNL > 4  || OprogStatus.useNNL < 0)
   {
     printf("[-D MD_SUPERELLIPSOID] useNNL must be between 0 and 4 (3=SPNNL and 1=usual NNL)\n");
     exit(-1);
   }
#else
 if (OprogStatus.useNNL > 4 || OprogStatus.useNNL < 0)
   {
     printf("[HEFLEX] useNNL must be between 0 and 2 (1=usual NNL, 0=disabled)\n");
     exit(-1);
   }
#endif
#endif
#if defined(MD_ABSORPTION) && defined(MD_EDHEFLEX_OPTNNL)
 if (OprogStatus.optnnl)
   {
     printf("[WARNING] you can not use optimal NNL with -DMD_ABSORPTION, I disabled optimal NNL\n");
     OprogStatus.optnnl = 0;
   }
 listtmp = malloc(sizeof(int)*OprogStatus.nebrTabFac);
#endif
 numcols=0;
 for (i=0;;i++)
   {
     if (!strcmp(colsFlex[i],""))
       break;
     numcols++;
   }
 OprogStatus.stepSDB = OprogStatus.stepSDA;
#endif
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
#ifdef MD_EDHEFLEX_WALL
  if (OprogStatus.epsdPlane==-1.0)
    OprogStatus.epsdPlane = OprogStatus.epsd;
  if (OprogStatus.epsdFastPlane==-1.0)
    OprogStatus.epsdFastPlane = OprogStatus.epsdFast;
#endif
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
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
#ifdef MD_LXYZ
  for (kk=0; kk < 3; kk++)
    {
      invL[kk] = 1.0/L[kk];
      L2[kk] = 0.5*L[kk];
    }
#else
  invL = 1.0/L;
  L2 = 0.5*L;
#endif
#ifdef MD_GRAVITY
  rcmz = -Lz*0.5;
  Lz2 = Lz*0.5;
#endif
  if (OprogStatus.eventMult==1)
    {
      printf("[ERROR] eventMult parameter must be greater than 1, typical values are 20-30\n");
      exit(-1);
    }
#ifdef MD_SPHERICAL_WALL
  poolSize = OprogStatus.eventMult*Oparams.parnum+2*Oparams.parnum;
#else
  poolSize = OprogStatus.eventMult*Oparams.parnum;
#endif
#ifndef EDHE_FLEX
  m = Oparams.m;
  Mtot = Oparams.m[0]*parnumA+Oparams.m[1]*parnumB;
  invmA = 1.0/Oparams.m[0];
  invmB = 1.0/Oparams.m[1];
  Mred[0][0] = Mred[1][1] = 0.5;
  Mred[0][1] = Mred[1][0] = (Oparams.m[0]*Oparams.m[1])/(Oparams.m[0]+Oparams.m[1]);
#endif
#if 0
  printf("massA: %f massB: %f sigmaA:%f sigmaB:%f sigmaAB:%f\n", Oparams.m[0], Oparams.m[1],
	 Oparams.sigma[0][0], Oparams.sigma[1][1], Oparams.sigma[0][1]);
#endif
#if defined(MD_GRAVITY)
  g2 = 0.5*Oparams.ggrav;
  mgA = Oparams.m[0]*Oparams.ggrav; 
  mgB = Oparams.m[1]*Oparams.ggrav;
#endif
#if defined(MD_POLYDISP) 
  if (Oparams.parnumA < Oparams.parnum)
    {
      printf("ERROR: Oparams.parnum has to be equal to Oparams.parnumA with polydispersity!\n");
      exit(-1);
    }
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
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
#ifdef MD_HE_PARALL
  if (my_rank == 0)
    tree = AllocMatI(10, poolSize);
#else
#ifdef MD_CALENDAR_HYBRID
  tree = AllocMatI(16, poolSize);
#else
  tree = AllocMatI(13, poolSize);
#endif
#endif
#ifdef EDHE_FLEX
  //printf("maxbonds=%d\n", OprogStatus.maxbonds);
  if (Oparams.maxbondsSaved==-1)
    {
#ifdef MD_LL_BONDS
      bonds = AllocMatLLI(Oparams.parnum, OprogStatus.maxbonds);
#else
      bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
#endif
      numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
#ifdef MD_SPHERICAL_WALL
      allocBondsSphWall();
#endif
    }
#else
#ifdef MD_LL_BONDS
  bonds = AllocMatLLI(Oparams.parnum, OprogStatus.maxbonds);
#else
  bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
#endif  
  numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
#endif

#ifdef MD_SPHERICAL_WALL
#ifdef MD_LL_BONDS
  bondscache = malloc(sizeof(long long int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#else
  bondscache = malloc(sizeof(int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#endif
#else
#ifdef MD_LL_BONDS
  bondscache = (long long int *) malloc(sizeof(long long int)*OprogStatus.maxbonds);
#else
  bondscache = (int *) malloc(sizeof(int)*OprogStatus.maxbonds);
#endif
#endif  
#else
#ifdef MD_HE_PARALL
  if (my_rank == 0)
    tree = AllocMatI(10, poolSize);
#else
#ifdef MD_CALENDAR_HYBRID
  tree = AllocMatI(13, poolSize);
#else
  tree = AllocMatI(10, poolSize);
#endif
#endif
#endif
#ifdef MD_HE_PARALL
  if (my_rank == 0)
    {
      treeTime = malloc(sizeof(double)*poolSize);
      treeRxC  = malloc(sizeof(double)*poolSize);
      treeRyC  = malloc(sizeof(double)*poolSize);
      treeRzC  = malloc(sizeof(double)*poolSize);
    }
#else
  treeTime = malloc(sizeof(double)*poolSize);
  treeRxC  = malloc(sizeof(double)*poolSize);
  treeRyC  = malloc(sizeof(double)*poolSize);
  treeRzC  = malloc(sizeof(double)*poolSize);
#endif
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
  //angM = malloc(sizeof(double)*Oparams.parnum);
  //phi0 = malloc(sizeof(double)*Oparams.parnum);
  //psi0 = malloc(sizeof(double)*Oparams.parnum);
  costheta0 = malloc(sizeof(double)*Oparams.parnum);
  sintheta0 = malloc(sizeof(double)*Oparams.parnum);
  theta0 =    malloc(sizeof(double)*Oparams.parnum);
  psi0   =    malloc(sizeof(double)*Oparams.parnum);
  phi0   =    malloc(sizeof(double)*Oparams.parnum);
  angM   =    malloc(sizeof(double)*Oparams.parnum);
#ifdef MD_ABSORP_POLY
  oldTypeOfPart = malloc(sizeof(int)*Oparams.parnum);
#endif
  REt = matrix(3,3);
  REtA = matrix(3,3);
  REtB = matrix(3,3);
  RE0 = matrix(3,3);
  Rdot = matrix(3,3);
  RM = malloc(sizeof(double**)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++) 
    {
#ifdef MD_MATRIX_CONTIGOUS
      /* alloca R in maniera contigua */
      if (i==0)
	{
  	  RM[i] = malloc(sizeof(double*)*3);
	  RM[i][0] = malloc(sizeof(double)*Oparams.parnum*9);
	  RM[i][1] = RM[i][0] + 3;
	  RM[i][2] = RM[i][1] + 3;
	}
      else
	{
	  RM[i] = malloc(sizeof(double*)*3);
	  RM[i][0] = RM[i-1][2] + 3;
	  RM[i][1] = RM[i][0] + 3;
	  RM[i][2] = RM[i][1] + 3;
	}
#else
      RM[i] = matrix(3, 3);
#endif
    }
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
#if 0
      double dist;
      dist = sqrt(Sqr(rx[i]-rx[1])+Sqr(ry[i]-ry[1])+Sqr(rz[i]-rz[1]));
      if (i > 1 && sqrt(Sqr(rx[i]-rx[1])+Sqr(ry[i]-ry[1])+Sqr(rz[i]-rz[1])) > 19.5)
	  {	  
	  printf("i=%d BOH dist=%.15G\n",i,sqrt(Sqr(rx[i]-rx[1])+Sqr(ry[i]-ry[1])+Sqr(rz[i]-rz[1])));
			 
	  }
      //else if (i > 1 && dist > 19.0)
	//fprintf(stderr,"%.15G %.15G %.15G @ 0.5\n", rx[i], ry[i], rz[i]);
#endif
#ifdef MD_MATRIX_CONTIGOUS
      /* alloca R in maniera contigua */
      if (i==0)
	{
  	  R[i] = malloc(sizeof(double*)*3);
	  R[i][0] = malloc(sizeof(double)*Oparams.parnum*9);
	  R[i][1] = R[i][0] + 3;
	  R[i][2] = R[i][1] + 3;
	}
      else
	{
	  R[i] = malloc(sizeof(double*)*3);
	  R[i][0] = R[i-1][2] + 3;
	  R[i][1] = R[i][0] + 3;
	  R[i][2] = R[i][1] + 3;
	}
#else
      R[i] = matrix(3, 3);
#endif
#if defined(MD_PATCHY_HE)||defined(EDHE_FLEX)
      lastbump[i].mol = -1;
      lastbump[i].at = -1;
#else
      lastbump[i] = -1;
#endif
      lastcol[i] = 0.0;
    }
#ifdef EDHE_FLEX
#ifdef MD_SPOT_GLOBAL_ALLOC
  for (pt = 0; pt < Oparams.ntypes; pt++)
    {
      if (first || typesArr[pt].nspots > maxsp)
	{
	  maxsp = typesArr[pt].nspots;
	  first = 0;
	}
    }
  maxsp += MD_SPNNL_NUMSP+1;
  printf("============>maxsp=%d\n", maxsp);
  distsSq = malloc(sizeof(double)*maxsp);
  ratA = (double**)malloc(sizeof(double*)*maxsp);
  ratA[0] = (double*)malloc(sizeof(double)*3*maxsp);
  ratB = (double**)malloc(sizeof(double*)*maxsp);
  ratB[0] = (double*)malloc(sizeof(double)*3*maxsp);
  for (k = 1; k < maxsp; k++)
    {
      ratA[k] = ratA[k-1] + 3;
      ratB[k] = ratB[k-1] + 3;
    }
  for (k = 0; k < 6; k++)
    {
      tocheckP[k] = malloc(sizeof(int)*maxsp);
      dorefineP[k] = malloc(sizeof(int)*maxsp);
      crossedP[k]  = malloc(sizeof(int)*maxsp);
      distsOldP[k] = malloc(sizeof(double)*maxsp);
      distsP[k]    = malloc(sizeof(double)*maxsp);
      t2arrP[k]    = malloc(sizeof(double)*maxsp);
      maxddotiP[k]  = malloc(sizeof(double)*maxsp);
#ifndef MD_BASIC_DT
      distsOld2P[k] = malloc(sizeof(double)*maxsp);
#endif
    }
  /* inizializzo questi array globali per evitare (comunque inutili) lamentele di valgrind */
  printf("[INFO] GLOBAL ALLOCATION OF SPOTS VARIABLES\n");
  for (i=0; i < maxsp; i++)
    {
      distsSq[i] = 0.0;
      for (k=0; k < 3; k++)
	{
	  ratA[i][k] = ratB[i][k] = 0;

	}
      for (k = 0; k < 6; k++)
	{

	  t2arrP[k][i] = distsP[k][i] = maxddotiP[k][i] = distsOldP[k][i] = 0.0;
#ifndef MD_BASIC_DT
	  distsOld2P[k][i] = 0.0;
#endif
	  crossedP[k][i] = tocheckP[k][i] = dorefineP[k][i] = crossedP[k][i] = 0;
	}
    }
#endif
  /* check moments of inertia */
  for (pt=0; pt < Oparams.ntypes; pt++)
    {
      if (typesArr[pt].I[0]!=typesArr[pt].I[1])
	{
	  printf("ERROR: Moments of inertia along x-axis and y-axis must be equal!\n");
	  printf("Aborting...\n");
	  exit(-1);
	}
      MD_DEBUG50(printf("semi-assi del tipo %d=%f %f %f\n", pt, typesArr[pt].sax[0], typesArr[pt].sax[1],
			typesArr[pt].sax[2]);)
    } 
  maxnbonds = get_max_nbonds();
#ifdef MD_EDHEFLEX_OPTNNL
  if (OprogStatus.optnnl)
    {
      printf("[INFO] Enabling optimal choice of linked lists cutoff for NNL\n");
    }
#endif
#ifdef MD_SPOT_GLOBAL_ALLOC
  /* l'array bonds � al massimo long long int e quindi bisogna evitare overflow */
  maxp = ((long long int)MAX_ALLOWED_INT - (((long long int)NA)*maxsp + maxsp)) / ((long long int) NANA);
  //printf("MAX_ALLOWED_INT=%lld maxsp=%ld NA=%ld", MAX_ALLOWED_INT, maxsp, NA);
  //printf("NA*maxsp+maxsp=%d\n", NA*maxsp+maxsp);
  if (maxnbonds > 0 && Oparams.parnum > maxp)
    {
      printf("I am sorry but actually I can not simulate more than %lld particles\n", maxp);
      printf("If you want you can decrease NA=%lld in ellipsoid.h to increase maximum number of possible particles\n",
	     (long long int) NA);
      exit(-1);

    }
  else if (maxnbonds > 0)
    printf("[INFO] maximum number of allowed particles is %lld\n", maxp);
#endif
  //printf("maxnbonds:%d\n OprogStatus.maxbonds: %d\n", maxnbonds, OprogStatus.maxbonds);
  //if (OprogStatus.maxbonds > maxnbonds)
    //maxnbonds = OprogStatus.maxbonds;
  dofTot = get_dof_flex(0);
  mapbondsaFlex = (int*)malloc(sizeof(int)*maxnbonds);
  mapbondsbFlex = (int*)malloc(sizeof(int)*maxnbonds);
  mapBheightFlex = (double*) malloc(sizeof(double)*maxnbonds);
  mapBhinFlex    = (double*)malloc(sizeof(double)*maxnbonds);
  mapBhoutFlex   = (double*)malloc(sizeof(double)*maxnbonds);
  mapSigmaFlex   = (double*)malloc(sizeof(double)*maxnbonds);
  if (OprogStatus.optbm)
    {
      printf("[INFO] Optimizing assign_bond_mapping(i,j) function\n");
      numbm = Oparams.ntypes*Oparams.ntypes;
      nbondsFlexS = malloc(sizeof(int)*numbm);
      mapbondsaFlexS = malloc(sizeof(int*)*numbm);
      mapbondsbFlexS = malloc(sizeof(int*)*numbm);
      mapBheightFlexS = malloc(sizeof(double*)*numbm);
      mapBhinFlexS = malloc(sizeof(double*)*numbm);
      mapBhoutFlexS = malloc(sizeof(double*)*numbm);
      mapSigmaFlexS = malloc(sizeof(double*)*numbm);
      for (nl = 0; nl < numbm; nl++)
	{
	  nbondsFlexS[nl] = -1;
	  mapbondsaFlexS[nl] = malloc(sizeof(int)*maxnbonds);
	  mapbondsbFlexS[nl] = malloc(sizeof(int)*maxnbonds); 
	  mapBheightFlexS[nl] = malloc(sizeof(double)*maxnbonds);
	  mapBhinFlexS[nl] = malloc(sizeof(double)*maxnbonds);
	  mapBhoutFlexS[nl] = malloc(sizeof(double)*maxnbonds);
	  mapSigmaFlexS[nl] = malloc(sizeof(double)*maxnbonds);
	}
    }

  dists    = (double*)malloc(maxnbonds*sizeof(double));
  distsOld = (double*)malloc(maxnbonds*sizeof(double));
  distsOld2= (double*)malloc(maxnbonds*sizeof(double));
  t2arr    = (double*)malloc(maxnbonds*sizeof(double));
  maxddoti = (double*)malloc(maxnbonds*sizeof(double));
  //crossed  = (int*)malloc(maxnbonds*sizeof(int));
  tocheck  = (int*)malloc(maxnbonds*sizeof(int));
  dorefine = (int*)malloc(maxnbonds*sizeof(int));
  crossed  = (int*)malloc(maxnbonds*sizeof(int));
  negpairs = (int*)malloc(maxnbonds*sizeof(int));
  check_conf();
#endif
#ifdef MD_SUPERELLIPSOID
  volSQ = malloc(Oparams.ntypes*sizeof(double));
#endif
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
  else if (OprogStatus.CMreset==-3)
   {
      resetCM(0);
   }
  else if (OprogStatus.CMreset==-4)
    {
      printf("[INFO] setting translational and angular velocities\n");
      comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
      angvelMB();
    }


#if defined(EDHE_FLEX) && defined(MD_HANDLE_INFMASS)
#ifdef MD_SPHERICAL_WALL
  for (i=0; i < Oparams.parnum; i++)
    {
      /* N.B. il tipo numero Oparams.ntypes-1 deve essere il muro sferico (interno)! */
      if (typeOfPart[i] == Oparams.ntypes-2)
	{
	  sphWall = i;
	}
      /* ntypes-1 � il muro esterno */
      if (typeOfPart[i] == Oparams.ntypes-1)
	{
	  sphWallOuter = i;	  
	}
    }
#endif
  if (newSim)
    {  
      for (i=0; i < Oparams.parnum; i++)
	{
	  if (is_infinite_mass(i))
	    {
	      vx[i] = 0.0;
	      vy[i] = 0.0;
	      vz[i] = 0.0;
	    }
	  if (is_infinite_Itens(i))
	    {
	      Mx[i] = 0.0;
	      My[i] = 0.0;
	      Mz[i] = 0.0;
	      wx[i] = 0.0;
	      wy[i] = 0.0;
	      wz[i] = 0.0;
	    }
	}

    }
#endif
#ifdef MD_EDHEFLEX_2D
  for (i=0; i < Oparams.parnum; i++)
    {
      vz[i] = 0.0;
    }
#endif
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
#ifdef MD_PROTEIN_DESIGN
  if (newSim && nativeConfYN==1)
    strcpy(OprogStatus.nativeConf, nativeConf);
  if (strcmp(OprogStatus.nativeConf, ""))
    read_native_conf();
#endif
#ifdef MD_RABBIT
#ifdef MD_LXYZ
  nbins_rhoz = L[2] / OprogStatus.rhozBinSize; 
#else
  nbins_rhoz = L / OprogStatus.rhozBinSize; 
#endif
  rhoz = malloc(sizeof(double)*nbins_rhoz);
#endif
  if (newSim)
    {
      FILE *f;
#ifdef MD_BIG_DT
      OprogStatus.refTime = 0.0;
#endif
      Oparams.time=0.0;
      /* truncate file to zero lenght */
#ifdef MD_SAVE_DISTANCE
      f = fopenMPI(absMisHD("distance.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("distance-pred.dat"), "w+");
      fclose(f);
#endif
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
      f = fopenMPI(absMisHD("radius_of_gyration.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("energy.dat"), "w+");
      fclose(f);
#ifdef MD_PROTEIN_DESIGN
      if (strcmp(OprogStatus.nativeConf,""))
	{
	  f = fopenMPI(absMisHD("ordpara.dat"), "w+");
	  fclose(f);
	}
#endif
#ifdef MD_RABBIT
      f = fopenMPI(absMisHD("bi-mono-bonds.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("rates.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("rhoz.dat"), "w+");
      fclose(f);
#endif
#ifdef MD_ABSORPTION
      f = fopenMPI(absMisHD("absorption.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("buffer.dat"), "w+");
      fclose(f);
#endif
#endif
      f = fopenMPI(absMisHD("MSDA.dat"), "w+");
      fclose(f);
      f = fopenMPI(absMisHD("rotMSDA.dat"), "w+");
      fclose(f);
#ifdef MD_CALC_DPP
      f = fopenMPI(absMisHD("MSDAxyz.dat"), "w+");
      fclose(f);
#endif
      if (Oparams.parnum > Oparams.parnumA)
	{
	  f = fopenMPI(absMisHD("MSDB.dat"), "w+");
	  fclose(f);
	  f = fopenMPI(absMisHD("rotMSDB.dat"), "w+");
	  fclose(f);
#ifdef MD_CALC_DPP
    	  f = fopenMPI(absMisHD("MSDBxyz.dat"), "w+");
	  fclose(f);
#endif
 
	}
      f = fopenMPI(absMisHD("temp.dat"), "w+");
      fclose(f);
#ifdef MD_INELASTIC
      f = fopenMPI(absMisHD("temp_granular.dat"), "w+");
      fclose(f);
#endif
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
#ifdef MD_CALC_DPP
  	  OprogStatus.sumdx[i] = 0.0;
  	  OprogStatus.sumdy[i] = 0.0;
  	  OprogStatus.sumdz[i] = 0.0;
	  store_last_u(i);
#if 0
  	  OprogStatus.lastu1x[i] = 0.0;
  	  OprogStatus.lastu1y[i] = 0.0;
  	  OprogStatus.lastu1z[i] = 0.0;
  	  OprogStatus.lastu2x[i] = 0.0;
  	  OprogStatus.lastu2y[i] = 0.0;
  	  OprogStatus.lastu2z[i] = 0.0;
  	  OprogStatus.lastu3x[i] = 0.0;
  	  OprogStatus.lastu3y[i] = 0.0;
  	  OprogStatus.lastu3z[i] = 0.0;
#endif
#endif
	}
      OprogStatus.nextcheckTime += fabs(OprogStatus.rescaleTime);
      OprogStatus.nextSumTime += OprogStatus.intervalSum;
      if (OprogStatus.storerate > 0.0)
	OprogStatus.nextStoreTime = OprogStatus.storerate;
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
#ifdef EDHE_FLEX
#ifdef MD_SUPERELLIPSOID
  init_gauleg_weights();
  SQvolPrefact = malloc(sizeof(double)*Oparams.ntypes);
  for (i=0; i < Oparams.parnum; i++)
    {
      /* to calculate geometrical factor of SQ volume
	 we need unitary semiaxis, later (below) we will set them
	 to proper values */
      axa[i] = 1.0;
      axb[i] = 1.0;
      axc[i] = 1.0;
    }
  for (pt = 0; pt < Oparams.ntypes; pt++)
    {
      /* se non si tratta di una superquadrica alloca il prefattore � semlicemente
	 4*pi/3 */
      if (typesArr[pt].n[0]==2.0 && typesArr[pt].n[1]==2.0 &&
	      typesArr[pt].n[2]==2.0)
	SQvolPrefact[pt] = 8.0*acos(0.0)/3.0;
      else
	SQvolPrefact[pt] = calc_SQ_volprefact(pt);
      printf("[type=%d] Volume SQ prefactor = %.15G\n", pt, SQvolPrefact[pt]);
    }
#endif
	
  if (OprogStatus.useNNL == 3 || OprogStatus.useNNL == 4)  
    { 
      buildSPNNL_spots();
    }
  a0I = malloc(sizeof(double)*Oparams.parnum);
#endif
  maxax = malloc(sizeof(double)*Oparams.parnum);
  scdone = malloc(sizeof(int)*Oparams.parnum);
  if (OprogStatus.useNNL)
    {
      printf("I'm going to use NNL, good choice to go fast :)\n");
#ifdef MD_HE_PARALL
      if (my_rank == 0)	
 	nebrTab = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
#else
      nebrTab = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
#endif
    }
  for (i=0; i < Oparams.parnumA; i++)
    {
#ifdef MD_HE_PARALL
      if (my_rank == 0 && OprogStatus.useNNL)	
	{
	  nebrTab[i].len = 0;
	  nebrTab[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	  //nebrTab[i].shift = matrix(OprogStatus.nebrTabFac, 3);
	  nebrTab[i].shift = NULL;
	  nebrTab[i].R = matrix(3, 3);
	}
#else
      if (OprogStatus.useNNL)
	{
	  nebrTab[i].len = 0;
	  nebrTab[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	  //nebrTab[i].shift = matrix(OprogStatus.nebrTabFac, 3);
	  /* tanto shift non viene usata per cui � inutile sprecare memoria*/ 
	  nebrTab[i].shift = NULL;
	  nebrTab[i].R = matrix(3, 3);
	}
#endif
      scdone[i] = 0;
#if defined(MD_POLYDISP)
      /* assegna i semiassi usando una gaussiana con deviazione standard fissata da OprogStatus.polydisp (%) */
      if (newSim)
	{
#ifdef MD_POLYDISP_XYZ
	  if (OprogStatus.polydispX > 0.0 || OprogStatus.polydispY > 0.0 || OprogStatus.polydispZ > 0.0)
	    {
	      do
		{
		  /* N.B. i semiassi vengono scalati di un fattore casuale in maniera non isotropa */
		  stocvar = gauss();
		  axaP[i] = (OprogStatus.polydispX*stocvar + 1.0)* Oparams.a[0]; 
		  stocvar = gauss();
		  axbP[i] = (OprogStatus.polydispY*stocvar + 1.0)* Oparams.b[0];
		  stocvar = gauss();
		  axcP[i] = (OprogStatus.polydispZ*stocvar + 1.0)* Oparams.c[0];
		}
	      while ( axaP[i] < Oparams.a[0]*(1.0 - OprogStatus.polycutoff*OprogStatus.polydispX) ||
		      axaP[i] > Oparams.a[0]*(1.0 + OprogStatus.polycutoff*OprogStatus.polydispX) ||
		      axbP[i] < Oparams.b[0]*(1.0 - OprogStatus.polycutoff*OprogStatus.polydispY) ||
		      axbP[i] > Oparams.b[0]*(1.0 + OprogStatus.polycutoff*OprogStatus.polydispY) ||
		      axcP[i] < Oparams.c[0]*(1.0 - OprogStatus.polycutoff*OprogStatus.polydispZ) ||
		      axcP[i] > Oparams.c[0]*(1.0 + OprogStatus.polycutoff*OprogStatus.polydispZ) );

	      //printf("%.15G\n", radii[i]);
	      axa[i] = axaP[i];
	      axb[i] = axbP[i];
	      axc[i] = axcP[i];
	    }
	  else
	    {
	      axa[i] = axaP[i];
	      axb[i] = axbP[i];
	      axc[i] = axcP[i];
	    }
#else
	  if (OprogStatus.polydisp > 0.0)
	    {
	      /* notare che le seguenti condizioni non dipendono dai semiassi ma solo dal valore restituito
	       * da gauss() quindi basta controllare solo uno dei tre semiassi. */
#if 1
	     do
	       {
		 /* N.B. i semiassi vengono scalati di un fattore casuale ma in maniera
		  * isotropa. */
		 stocvar = gauss();
       		 axaP[i] = (OprogStatus.polydisp*stocvar + 1.0)* Oparams.a[0]; 
		 axbP[i] = (OprogStatus.polydisp*stocvar + 1.0)* Oparams.b[0];
		 axcP[i] = (OprogStatus.polydisp*stocvar + 1.0)* Oparams.c[0];
	       }
	     while ( axaP[i] < Oparams.a[0]*(1.0 - OprogStatus.polycutoff*OprogStatus.polydisp) ||
		    axaP[i] > Oparams.a[0]*(1.0 + OprogStatus.polycutoff*OprogStatus.polydisp) );
#else
	     /* this is just for testing purpose */
	     if (i < 128)
	       {
		 axaP[i] = Oparams.a[0];
		 axbP[i] = Oparams.b[0];
		 axcP[i] = Oparams.c[0];
	       }
	     else
	       {
	       	 axaP[i] = Oparams.a[1];
		 axbP[i] = Oparams.b[1];
		 axcP[i] = Oparams.c[1];
	       }
#endif

	      //printf("%.15G\n", radii[i]);
	     axa[i] = axaP[i];
	     axb[i] = axbP[i];
	     axc[i] = axcP[i];
	    }
	  else
	    {

	      axa[i] = axaP[i];
	      axb[i] = axbP[i];
	      axc[i] = axcP[i];
	    }
#endif
	}
      else
	{
	  axa[i] = axaP[i];
	  axb[i] = axbP[i];
	  axc[i] = axcP[i];
	}
      //printf("$$$ axes[%d]=(%f,%f,%f)\n", i, axaP[i], axbP[i], axcP[i]);
#else
#ifndef EDHE_FLEX
      axa[i] = Oparams.a[0];
      axb[i] = Oparams.b[0];
      axc[i] = Oparams.c[0];
#endif
#endif
    }
#ifdef EDHE_FLEX
  for (i=0; i < Oparams.parnum; i++)
    {
      axa[i] = typesArr[typeOfPart[i]].sax[0];
      axb[i] = typesArr[typeOfPart[i]].sax[1];
      axc[i] = typesArr[typeOfPart[i]].sax[2];
    }
#endif
#ifndef EDHE_FLEX
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
#ifdef MD_HE_PARALL
      if (my_rank == 0 && OprogStatus.useNNL)
	{
	  nebrTab[i].len = 0;
	  nebrTab[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	  nebrTab[i].R = matrix(3, 3);
	}
#else
      if (OprogStatus.useNNL)
	{
	  nebrTab[i].len = 0;
	  nebrTab[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	  nebrTab[i].R = matrix(3, 3);
	}
#endif
      scdone[i] = 0;
#ifndef EDHE_FLEX
      axa[i] = Oparams.a[1];
      axb[i] = Oparams.b[1];
      axc[i] = Oparams.c[1];
#endif
    }
#endif
#ifndef EDHE_FLEX
#ifdef MD_LXYZ
  printf(">>>> phi=%.12G L=%f %f %f (%f,%f,%f)\n", phiIni=calc_phi(), L[0], L[1], L[2], Oparams.a[0], Oparams.b[0], Oparams.c[0]); 
#else
  printf(">>>> phi=%.12G L=%f (%f,%f,%f)\n", phiIni=calc_phi(), L, Oparams.a[0], Oparams.b[0], Oparams.c[0]); 
#endif
  if (Oparams.parnumA < Oparams.parnum)
    printf("semi-axes of B (%f, %f ,%f)\n",Oparams.a[1], Oparams.b[1], Oparams.c[1]);
#else
#ifdef MD_LXYZ
  printf("[FLEX] phi=%.12G L=%f %f %f\n", phiIni=calc_phi(), L[0], L[1], L[2]); 
#else
  printf("[FLEX] phi=%.12G L=%f\n", phiIni=calc_phi(), L); 
#endif
#endif
  /* evaluation of principal inertia moments*/ 
  for (a = 0; a < 2; a++)
    {
#ifndef EDHE_FLEX
      invaSq[a] = Sqr(1/Oparams.a[a]);
      invbSq[a] = Sqr(1/Oparams.b[a]);
      invcSq[a] = Sqr(1/Oparams.c[a]);
#endif
#ifdef MD_ASYM_ITENS
      ItensD[a][0] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.b[a])+Sqr(Oparams.c[a]));
      ItensD[a][1] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.a[a])+Sqr(Oparams.c[a]));
      ItensD[a][2] = 1.0;//(1.0/5.0)*Oparams.m[a]*(Sqr(Oparams.a[a])+Sqr(Oparams.b[a]));
#endif
    };
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

  /* maxax e' il diametro del centroide */
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
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
#if defined(MD_POLYDISP) 
      if (axaP[i] > maxax[i])
	maxax[i] = axaP[i];
      if (axbP[i] > maxax[i])
	maxax[i] = axbP[i];
      if (axcP[i] > maxax[i])
	maxax[i] = axcP[i];
#elif defined(EDHE_FLEX)
#ifdef MD_SUPERELLIPSOID
#if 0
      if (typesArr[typeOfPart[i]].sax[0] > maxax[i])
	maxax[i] = typesArr[typeOfPart[i]].sax[0];
      if (typesArr[typeOfPart[i]].sax[1] > maxax[i])
	maxax[i] = typesArr[typeOfPart[i]].sax[1];
      if (typesArr[typeOfPart[i]].sax[2] > maxax[i])
	maxax[i] = typesArr[typeOfPart[i]].sax[2];
#else 
      /* per i superellissoidi il centroide non puo' essere uguale al raggio maggiore,
	 cos� lo sceglo pari alla met� della diagonale maggiore del parallelepipedo
	 con i lati pari al doppio dei semi-assi */
      maxax[i] = sqrt(Sqr(typesArr[typeOfPart[i]].sax[0])+Sqr(typesArr[typeOfPart[i]].sax[1])+
	Sqr(typesArr[typeOfPart[i]].sax[2]));
#endif
#else
      maxax[i] = sqrt(Sqr(typesArr[typeOfPart[i]].sax[0])+Sqr(typesArr[typeOfPart[i]].sax[1])+
	Sqr(typesArr[typeOfPart[i]].sax[2]));
#endif
      maxSpots = eval_max_dist_for_spots(typeOfPart[i]);
      if (maxSpots > maxax[i])
	maxax[i] = maxSpots;
      //printf("maxax[%d]:%f maxSpots:%f\n", i, maxax[i], maxSpots);
#else
      a=(i<Oparams.parnumA)?0:1;
      if (Oparams.a[a] > maxax[i])
	maxax[i] = Oparams.a[a];
      if (Oparams.b[a] > maxax[i])
	maxax[i] = Oparams.b[a];
      if (Oparams.c[a] > maxax[i])
	maxax[i] = Oparams.c[a];
#endif
      //printf("distSPA=%.15G distSPB=%.15G\n", distSPA, distSPB);
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
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

#ifdef MD_SPHERICAL_WALL
      if (maxax[i] > MAXAX && typeOfPart[i] != Oparams.ntypes-1 && typeOfPart[i] != Oparams.ntypes-2)
	MAXAX = maxax[i];
#else
      if (maxax[i] > MAXAX)
	MAXAX = maxax[i];
#endif
      //printf("maxax aft[%d]: %.15G\n", i, maxax[i]);
    }
#ifdef EDHE_FLEX
  is_a_sphere_NNL = malloc(sizeof(int)*Oparams.parnum);
  find_spheres_NNL();
  for (i=0; i < Oparams.parnum; i++)
    {
      if (is_a_sphere_NNL[i])
	set_angmom_to_zero(i);
    }
  is2saveArr = malloc(sizeof(int)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    is2saveArr[i] = 0;
  par2saveArr();
#endif/* Calcoliamo rcut assumendo che si abbian tante celle quante sono 
   * le particelle */
#ifdef EDHE_FLEX
  if (!OprogStatus.useNNL)
    {
      /* giusto per inizializzare tali vettore ed evitare che valgrind si lamenti */
      for (i=0; i < Oparams.ntypes; i++)
	{
	  for (k=0; k < 3; k++)
	    {
	      typesArr[i].ppsax[k] = 0.0;
	      typesArr[i].ppr[k] = 0.0;
	    }
	}
    }
#endif
  if (newSim)
    {
      if (OprogStatus.useNNL)
	{
	  /* in questo modo rNebrShell � la "buccia" rispetto al minimo 
	   * parallelepipedo che include l'ellissoide pi� gli sticky spots */
#ifndef EDHE_FLEX
	  if (OprogStatus.targetPhi <= 0.0)
	    OprogStatus.rNebrShell += calc_shell();
	  
	  printf("[INFO] I've adjusted rNebrShell to %.15G\n", OprogStatus.rNebrShell);	  
#else
	  calc_encpp();
#endif	  
	  if (Oparams.rcut <= 0.0)
	    {
#ifdef EDHE_FLEX
	      Oparams.rcut = calc_nnl_rcut();
#elif defined(MD_POLYDISP) 
#ifdef MD_POLYDISP_XYZ
	      if (OprogStatus.polydispX > 0.0||OprogStatus.polydispY > 0.0||OprogStatus.polydispZ > 0.0)
		Oparams.rcut = calc_nnl_rcut();//*(1.0+OprogStatus.polydisp*OprogStatus.polycutoff);
	      else
		Oparams.rcut = calc_nnl_rcut();
#else
	      if (OprogStatus.polydisp > 0.0)
		Oparams.rcut = calc_nnl_rcut();//*(1.0+OprogStatus.polydisp*OprogStatus.polycutoff);
	      else
		Oparams.rcut = calc_nnl_rcut();
#endif
#else
	      Oparams.rcut = calc_nnl_rcut();
#endif
	      printf("[INFO] I've chosen rcut= %.15G\n", Oparams.rcut);
	    }
	}
      else
	{
#if defined(MD_POLYDISP) 
	  if (Oparams.rcut <= 0.0)
	    {
	      Oparams.rcut = MAXAX*1.01;
	    }
#else
	  if (Oparams.rcut <= 0.0 
#ifdef MD_MULTIPLE_LL
	      || OprogStatus.multipleLL
#endif
 	     )
	    {
#ifdef EDHE_FLEX
	      calc_encpp();
	      rcut2 = calc_nnl_rcut();
	      if (rcut2 > MAXAX*1.01)
		Oparams.rcut = MAXAX*1.01;
	      else
		Oparams.rcut = rcut2;
#else
	      Oparams.rcut = MAXAX*1.01;
#endif
	    }
#endif
	}
    }
  else 
    {
#ifdef EDHE_FLEX
      if (OprogStatus.useNNL
#ifdef MD_MULTIPLE_LL
	  ||OprogStatus.multipleLL
#endif
	 )
	calc_encpp();
#endif
    }
  printf("MAXAX: %.15G rcut: %.15G\n", MAXAX, Oparams.rcut);
  //Oparams.rcut = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
      int t1, t2;
      numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
      rcutMLL = malloc(sizeof(double)*numll);
      cellListMLL = malloc(sizeof(int*)*numll);
      cellsxMLL = malloc(sizeof(int)*numll);
      cellsyMLL = malloc(sizeof(int)*numll);
      cellszMLL = malloc(sizeof(int)*numll);
      ignoreMLL = malloc(sizeof(int)*numll);
      set_cells_size();

      for (k1=0; k1 < numll; k1++)
	{
	  get_types_from_nl(k1, &t1, &t2);
#ifdef MD_SPHERICAL_WALL
	  if (t1==typeOfPart[sphWall] || t1==typeOfPart[sphWallOuter] ||
	      t2==typeOfPart[sphWall] || t2==typeOfPart[sphWallOuter])
	    {
	      ignoreMLL[k1] = 1;
	      printf("I ignore list related to spherical walls (types=%d-%d)\n", t1, t2);
	      continue;
	    }
#endif
	  if (!may_interact_all_type(t1, t2))
	    {
	      ignoreMLL[k1] = 1;
	      printf("Types %d and %d are different and they do not interact hence I ignore this linked list\n", t1, t2);
	      continue;
	    }
	  else
	    {
	      ignoreMLL[k1] = 0;
	    }
#ifdef MD_LXYZ
    	  printf("[%d] L=%.15G %.15G %.15G Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d t1=%d t2=%d\n", k1, L[0], L[1], L[2],
		 rcutMLL[k1], cellsxMLL[k1], cellsyMLL[k1], cellszMLL[k1], t1, t2);
#else
	  printf("[%d] L=%.15G Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d t1=%d t2=%d\n", k1, L,
		 rcutMLL[k1], cellsxMLL[k1], cellsyMLL[k1], cellszMLL[k1], t1, t2);
#endif
	}
      for (k1 = 0; k1 < numll; k1++)
	{
	  if (ignoreMLL[k1])
	    continue;
	  ls = cellsxMLL[k1]*cellsyMLL[k1]*cellszMLL[k1];
	  cellListMLL[k1] = malloc(sizeof(int)*(ls+Oparams.parnum));
	  if (k1==0|| ls > mls)
	    mls = ls;
	}
#if 1
      inCellMLL = malloc(sizeof(int**)*Oparams.ntypes);
      for (k1 = 0; k1 < Oparams.ntypes; k1++)
	{
	  inCellMLL[k1] = malloc(sizeof(int*)*3);
	}
      for (k1=0; k1 < Oparams.ntypes; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    inCellMLL[k1][k2] = malloc(sizeof(int)*Oparams.parnum);
	  }
#else
      inCellMLL = malloc(sizeof(int**)*Oparams.ntypes);
      for (k1 = 0; k1 < Oparams.ntypes; k1++)
	{
	  inCellMLL[k1] = malloc(sizeof(int*)*3);
	}
     inCellMLL[0][0] = malloc(sizeof(int)*Oparams.parnum*3*Oparams.ntypes); 
     for (k1 = 0; k1 < Oparams.ntypes; k1++)
       for (k2 = 0; k2 < 3; k2++)
	 {
	   if (k1==0 && k2==0)
	     continue;
	   k2old = k2-1;
	   k1old = k1;
	   if (k2old < 0)
	     {	
	       k2old = 2;
	       k1old = k1-1;
	       if (k1old < 0)
		 k1old = Oparams.ntypes-1;
	     }
	   inCellMLL[k1][k2] = inCellMLL[k1old][k2old] + Oparams.parnum;
	 }	 
#endif
      crossevtodel = malloc(sizeof(int*)*Oparams.ntypes);
      for (k1=0; k1 < Oparams.ntypes; k1++)
	{
	  crossevtodel[k1] = malloc(sizeof(int)*Oparams.parnum);
	}
      for (k1 = 0; k1 < Oparams.ntypes; k1++)
	{
	  for (k2 = 0; k2 < Oparams.parnum; k2++)
	    crossevtodel[k1][k2] = -1;
	}
#ifdef MD_EDHEFLEX_OPTNNL
      if (OprogStatus.optnnl)
	{
	  //printf("====>mls=%d\n", mls);
	  cellList_NNL = malloc(sizeof(int)*(mls+Oparams.parnum));
	  inCell_NNL[0] = malloc(sizeof(int)*Oparams.parnum);
	  inCell_NNL[1]= malloc(sizeof(int)*Oparams.parnum);
	  inCell_NNL[2] = malloc(sizeof(int)*Oparams.parnum);
	  rxNNL = malloc(sizeof(double)*Oparams.parnum);
	  ryNNL = malloc(sizeof(double)*Oparams.parnum);
	  rzNNL = malloc(sizeof(double)*Oparams.parnum);
	}
#endif
    }
  else
    {
#ifdef MD_LXYZ
      cellsx = L[0] / Oparams.rcut;
      cellsy = L[1] / Oparams.rcut;
      cellsz = L[2] / Oparams.rcut;
#else
      cellsx = L / Oparams.rcut;
      cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
      cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
      cellsz = L / Oparams.rcut;
#endif 
#endif
      printf("Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", Oparams.rcut,
	     cellsx, cellsy, cellsz);
      cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
      inCell[0] = malloc(sizeof(int)*Oparams.parnum);
      inCell[1]= malloc(sizeof(int)*Oparams.parnum);
      inCell[2] = malloc(sizeof(int)*Oparams.parnum);
#ifdef MD_EDHEFLEX_OPTNNL
      if (OprogStatus.optnnl)
	{
	  cellList_NNL = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
	  inCell_NNL[0] = malloc(sizeof(int)*Oparams.parnum);
	  inCell_NNL[1]= malloc(sizeof(int)*Oparams.parnum);
	  inCell_NNL[2] = malloc(sizeof(int)*Oparams.parnum);
	  rxNNL = malloc(sizeof(double)*Oparams.parnum);
	  ryNNL = malloc(sizeof(double)*Oparams.parnum);
	  rzNNL = malloc(sizeof(double)*Oparams.parnum);
	}
#endif
    }
#else
#ifdef MD_LXYZ
      cellsx = L[0] / Oparams.rcut;
      cellsy = L[1] / Oparams.rcut;
      cellsz = L[2] / Oparams.rcut;
#else
      cellsx = L / Oparams.rcut;
      cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
      cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
      cellsz = L / Oparams.rcut;
#endif 
#endif
#ifdef ED_PARALL_DD
      dd_init(); 
#endif
      printf("Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", Oparams.rcut,
	     cellsx, cellsy, cellsz);
      cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
      inCell[0] = malloc(sizeof(int)*Oparams.parnum);
      inCell[1] = malloc(sizeof(int)*Oparams.parnum);
      inCell[2] = malloc(sizeof(int)*Oparams.parnum);
#ifdef MD_EDHEFLEX_OPTNNL
      if (OprogStatus.optnnl)
	{
	  cellList_NNL = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
	  inCell_NNL[0] = malloc(sizeof(int)*Oparams.parnum);
	  inCell_NNL[1]= malloc(sizeof(int)*Oparams.parnum);
	  inCell_NNL[2] = malloc(sizeof(int)*Oparams.parnum);
	  rxNNL = malloc(sizeof(double)*Oparams.parnum);
	  ryNNL = malloc(sizeof(double)*Oparams.parnum);
	  rzNNL = malloc(sizeof(double)*Oparams.parnum);
	}
#endif

#endif
 
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
#ifdef EDHE_FLEX
  if (Oparams.maxbondsSaved==-1)
    {
      for (i=0; i < Oparams.parnum; i++)
	{
	  numbonds[i] = 0;
	}
    }
#else
  for (i=0; i < Oparams.parnum; i++)
    {
      numbonds[i] = 0;
    }
#endif
#endif
#ifdef MD_LXYZ
  printf("L=%f %f %f parnum: %d parnumA: %d\n", L[0], L[1], L[2], Oparams.parnum, Oparams.parnumA);
#else
  printf("L=%f parnum: %d parnumA: %d\n", L, Oparams.parnum, Oparams.parnumA);
#endif
#if defined(MD_PATCHY_HE) && !defined(EDHE_FLEX)
  printf("sigmaSticky=%.15G\n", Oparams.sigmaSticky);
#endif
#ifdef MD_ASYM_ITENS
  for (i=0; i < Oparams.parnum; i++)
    {
      if (isSymItens(i))
	{
	  double II;
#ifdef EDHE_FLEX
	  II=typesArr[typeOfPart[i]].I[0];
#else
	  if (i < Oparams.parnumA)
	    II=Oparams.I[0][0];
	  else
	    II=Oparams.I[1][0];
#endif
	  wx[i]=Mx[i]/II;
	  wy[i]=My[i]/II;
	  wz[i]=Mz[i]/II;
	}
      else
	calc_omega(i, &(wx[i]), &(wy[i]), &(wz[i]));
      angM[i] = sqrt(Sqr(Mx[i])+Sqr(My[i])+Sqr(Mz[i]));
      upd_refsysM(i);
    }
#endif
#ifdef MD_PATCHY_HE
  printf("Energia potenziale all'inizio: %.15f\n", calcpotene());
#endif
  //exit(-1);
#ifdef MD_HE_PARALL
  slave_task();
#endif
#ifdef MD_GHOST_IGG
  if (Oparams.ghostsim)
    {
      printf("[GHOST SIMULATION]: Initializing structures...\n");    
      init_ghostArr();
    }
#if 0
  /* comunque alloca l'array poich� writeBinCoord_heflex() lo usa per scrivere!
     uso la calloc per inizializzare tutto a 0 ed evitare insulti da parte di valgrind */
  if (ghostInfoArr==NULL)
    ghostInfoArr = calloc(Oparams.parnum, sizeof(ghostInfo));
#endif
#endif
#ifdef MD_GRAZING_TRYHARDER
  printf("[INFO] Grazing Try-Harder code ENABLED!\n");
#else
  printf("[INFO] Grazing Try-Harder code DISABLED!\n");
#endif

  StartRun();

#ifdef MD_SPHERICAL_WALL
  printf("[SPHERICAL WALL] checking bonds\n");
  if (OprogStatus.checkGrazing)
    check_all_bonds();
  printf("check done all bonds OK\n");
#endif
 
  if (mgl_mode != 2)
    ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0 && mgl_mode!=2)
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
      /* truncate file to zero lenght */
#ifdef MD_GRAVITY
      FILE *f;
      f = fopenMPI(MD_HD_MIS "Vz2.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rcmz.dat", "w");
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rho.dat", "w");
      fclose(f);
#endif
      OprogStatus.DQxy = 0.0;
      OprogStatus.DQyz = 0.0;
      OprogStatus.DQzx = 0.0;
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

#if 0
	  OprogStatus.vcmx0[i] = vcmx;
	  OprogStatus.vcmy0[i] = vcmy;
	  OprogStatus.vcmz0[i] = vcmz;
#endif
#ifdef MD_DYNAMIC_OPROG
	  OprogStatus.lastcolltime[i] = 0.0;
	  OprogStatus.sumox[i] = 0.0;
	  OprogStatus.sumoy[i] = 0.0;
	  OprogStatus.sumoz[i] = 0.0;
#ifdef MD_CALC_DPP
	  OprogStatus.sumdx[i] = 0.0;
	  OprogStatus.sumdy[i] = 0.0;
	  OprogStatus.sumdz[i] = 0.0;
	  OprogStatus.lastu1x[i] = 0.0;
	  OprogStatus.lastu1y[i] = 0.0;
	  OprogStatus.lastu1z[i] = 0.0;
	  OprogStatus.lastu2x[i] = 0.0;
	  OprogStatus.lastu2y[i] = 0.0;
	  OprogStatus.lastu2z[i] = 0.0;
	  OprogStatus.lastu3x[i] = 0.0;
	  OprogStatus.lastu3y[i] = 0.0;
	  OprogStatus.lastu3z[i] = 0.0;
#endif
	  for (k=0; k < 3; k++)
	    OprogStatus.DR[i][k] = 0.0;
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
#ifdef EDHE_FLEX
 if (!strcmp(OprogStatus.par2save,"ALL") || !strcmp(OprogStatus.par2save,"all"))
    globSaveAll = 1;
  else
    globSaveAll = 0;
  if (!globSaveAll || OprogStatus.stripStore)
    {
      saveFullStore("StoreInit");
    }
#endif
  /* printf("Vol: %.15f Vol1: %.15f s: %.15f s1: %.15f\n", Vol, Vol1, s, s1);*/
#if defined(MD_CALC_VBONDING) && !defined(MD_STANDALONE) && !defined(MC_SIMUL)
  calc_vbonding();
#endif
}
extern double rA[3], rB[3];
#ifdef EDHE_FLEX
void save_coordtmp_ascii(unsigned char hdWhich)
{
  char fileop[1024], fileop2[1024];
  FILE *bf;
  char fileop3[1024];
  const char sepStr[] = "@@@\n";
  int i;	
  sprintf(fileop2 ,"COORD_TMP_ASCII%d", (int)hdWhich);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  saveFullStore(fileop);
}
#endif
#ifdef EDHE_FLEX
int is_valid_parname_progStatus(char *pn)
{
  /* salva solo quello che pu� servire per l'analisi */
  if (!strcmp(pn, "DR") ||
      !strcmp(pn, "sumox") ||
      !strcmp(pn, "sumoy") ||
      !strcmp(pn, "sumoz") ||
      !strcmp(pn, "refTime")||
      !strcmp(pn, "KK")||
      !strcmp(pn, "JJ")
     )
    return 1; 
  else 
    return 0;
}
int is_valid_parname_params(char *pn)
{
  /* salva solo quello che pu� servire per l'analisi */
   if (!strcmp(pn, "totSteps") ||
       !strcmp(pn, "time") ||
       !strcmp(pn, "curStep") 
      )
    return 1; 
  else 
    return 0;
}
#endif
#ifdef MD_PATCHY_HE
#ifdef MD_SPOT_GLOBAL_ALLOC
void BuildAtomPos(int i, double *rO, double **R, double **rat);
#else
void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3]);
#endif
#endif
#ifdef EDHE_FLEX
void par2saveArr(void)
{
  char *ns;
  int first=1, i;
  char s1[32], s2[32];
  if (!strcmp(OprogStatus.par2save,"all")||!strcmp(OprogStatus.par2save,"ALL"))
    return;

  ns = strtok(OprogStatus.par2save, ",");
  while(ns || first)
    {
      first = 0;
      if (!ns && first)
	ns = OprogStatus.par2save;
      if (sscanf(ns, "%[^-]-%s", s1, s2)==2)
	{
	  //printf("s1=%s s2=%s\n",s1,s2);
	  for (i=atoi(s1); i <= atoi(s2); i++)
	    {
	      is2saveArr[i] = 1;
	    }
	}
      else
	is2saveArr[atoi(ns)]=1;
      ns = strtok(NULL, ",");
    }
  //for (i=0; i < Oparams.parnum; i++)
    //{
      //printf("is_to_save[%d]=%d\n", i, is2saveArr[i]);
    //}
}
int is_to_save(int i)
{
  if (is2saveArr[i])
    return 1;
  else
    return 0;
}
#endif
#ifdef EDHE_FLEX
void fprintf_ranges(FILE *f, int A, int nr, rangeStruct* r)
{
  int kk;
  if (A==-2)
    {
      for (kk=0; kk < nr; kk++)
	{
	  if (r[kk].min==r[kk].max)
	    fprintf(f, "%d", r[kk].min);
	  else
	    fprintf(f, "%d-%d", r[kk].min, r[kk].max);	
	  if (kk==nr-1)
	    fprintf(f," ");
	  else
	    fprintf(f,","); 
	}
    }
  else
    fprintf(f, "%d ", A);
}
#endif
/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs, int saveAll)
{
  int i;
  int nn;
#ifndef MD_SPOT_GLOBAL_ALLOC
  double ratA[NA][3];
#endif
#ifdef MD_FOUR_BEADS
  int nn2;
  char bcolor[32];
  char beadcol[2][4][32] = {{"red","green","blue","orange"},{"red","green","blue","orange"}};
#endif
#ifdef MD_SUPERELLIPSOID
  const char tipodat2_mgl[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G @ %.15G %.15G %.15G C[%s] Q %.15G %.15G %.15G\n";
#else
  const char tipodat2_mgl[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G @ %.15G %.15G %.15G C[%s]\n";
#endif
  const char tipodat[] = "%.15G %.15G %.15G %.15G %.15G %.15G\n";
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n";
#ifdef EDHE_FLEX
  const char tipodat2_flex[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n";
  int j;
#ifdef MC_SIMUL
  R2u();
#endif
  if (!mgl_mode && !OprogStatus.stripStore)
    {
#ifdef MD_MULTIPLE_LL
      fprintf(fs, "RF ");
      for (i=0; i < Oparams.ntypes; i++)
	fprintf(fs, "%.8G ", typesArr[i].rcutFact);
      fprintf(fs, "\n");
#endif
      for (i=0; i < Oparams.ntypes; i++)
	fprintf(fs, "%d ", typeNP[i]);
      fprintf(fs, "\n");
      for (i=0; i < Oparams.ntypes; i++)
	{
	  /* write particles parameters */
	  fprintf(fs, "%.15G %.15G %.15G\n", typesArr[i].sax[0], typesArr[i].sax[1], typesArr[i].sax[2]); 
	  fprintf(fs, "%.15G %.15G %.15G\n", typesArr[i].n[0], typesArr[i].n[1], typesArr[i].n[2]);
	  fprintf(fs, "%.15G %.15G %.15G %.15G %d %d\n", typesArr[i].m, typesArr[i].I[0], typesArr[i].I[1],
		  typesArr[i].I[2], typesArr[i].brownian, typesArr[i].ignoreCore);
#if 0
	  fprintf(fs, "%.15G %.15G %.15G\n", typesArr[i].xoff[0], typesArr[i].xoff[1], typesArr[i].xoff[2]); 
#endif
	  /* write sticky spots parameters */
	  fprintf(fs, "%d %d\n", typesArr[i].nspots, typesArr[i].nhardobjs);
	  for (j = 0; j < typesArr[i].nspots; j++)
	    fprintf(fs, "%.15G %.15G %.15G %.15G\n", typesArr[i].spots[j].x[0],typesArr[i].spots[j].x[1],
		    typesArr[i].spots[j].x[2], typesArr[i].spots[j].sigma);
	  for (j = 0; j < typesArr[i].nhardobjs; j++)
	    fprintf(fs, "%.15G %.15G %.15G % .15G %.15G %.15G %.15G %.15G %.15G\n", 
		   typesArr[i].hardobjs[j].x[0],typesArr[i].hardobjs[j].x[1],typesArr[i].hardobjs[j].x[2], 
		   typesArr[i].hardobjs[j].sax[0],typesArr[i].hardobjs[j].sax[1],typesArr[i].hardobjs[j].sax[2],
		   typesArr[i].hardobjs[j].n[0],typesArr[i].hardobjs[j].n[1],typesArr[i].hardobjs[j].n[2]);
	} 
      /* write interactions */
      for (i=0; i < Oparams.ninters; i++)
	{
#if 0
	  fprintf (fs, "%d %d %d %d %.15G %.15G %.15G %d\n", intersArr[i].type1, intersArr[i].spot1,
  		   intersArr[i].type2, 
  		   intersArr[i].spot2, 
  		   intersArr[i].bheight, intersArr[i].bhin, intersArr[i].bhout, intersArr[i].nmax);
#else
	  fprintf_ranges(fs, intersArr[i].type1, intersArr[i].nr1, intersArr[i].r1);
	  fprintf(fs, "%d ", intersArr[i].spot1);
	  fprintf_ranges(fs, intersArr[i].type2, intersArr[i].nr2, intersArr[i].r2);
	  fprintf(fs, "%d ", intersArr[i].spot2);
	  fprintf (fs, "%.15G %.15G %.15G %d\n", 
		   intersArr[i].bheight, intersArr[i].bhin, intersArr[i].bhout, intersArr[i].nmax);
#endif
	}
      if (Oparams.nintersIJ > 0)
	{
    	  /* write interactions */
	  for (i=0; i < Oparams.nintersIJ; i++)
	    {
#if 0
	      fprintf (fs, "%d %d %d %d %.15G %.15G %.15G\n", intersArrIJ[i].i, intersArrIJ[i].spot1,
		       intersArrIJ[i].j, 
	      	       intersArrIJ[i].spot2, 
	      	       intersArrIJ[i].bheight, intersArrIJ[i].bhin, intersArrIJ[i].bhout);
#else
	      fprintf_ranges(fs, intersArrIJ[i].i, intersArrIJ[i].nr1, intersArrIJ[i].r1);
	      fprintf(fs, "%d ", intersArrIJ[i].spot1);
	      fprintf_ranges(fs, intersArrIJ[i].j, intersArrIJ[i].nr2, intersArrIJ[i].r2);
	      fprintf(fs, "%d ", intersArrIJ[i].spot2);
	      fprintf (fs, "%.15G %.15G %.15G\n",   
		       intersArrIJ[i].bheight, intersArrIJ[i].bhin, intersArrIJ[i].bhout);
#endif
	    } 
	}
      //fprintf(fs, "\n");
      if (Oparams.saveBonds)
	{
	  for (i = 0; i < Oparams.parnum; i++)
	    {
	      fprintf(fs, "%d ", numbonds[i]);
	      for (j = 0; j < numbonds[i]; j++)
		{
#ifdef MD_LL_BONDS
		  fprintf(fs, "%lld ", bonds[i][j]);
#else
		  fprintf(fs, "%d ", bonds[i][j]);
#endif
		}
	      fprintf(fs, "\n");
	    }
	}
    }
#if 0
#ifdef MD_GHOST_IGG
  if (Oparams.ghostsim)
    {
      for (i=0; i < Oparams.parnum; i++)
	{
	  fprintf(fs, "%d %d\n", ghostInfoArr[i].iggnum, ghostInfoArr[i].ghost_status);
	}
    }
#endif
#endif
#endif
  if (mgl_mode)
    {
#ifdef MD_LXYZ
      fprintf(fs, ".Vol: %f\n", L[0]*L[1]*L[2]);
#else
      fprintf(fs, ".Vol: %f\n", L*L*L);
#endif
      for (i = 0; i < Oparams.parnum; i++)
	{
#ifdef EDHE_FLEX
	  if (!saveAll && !globSaveAll)
	    {
	      if (!is_to_save(i))
		continue;
	    }
#ifdef MD_FOUR_BEADS
	  if (i >= typeNP[0])
	    continue;
#endif
#endif

#ifdef EDHE_FLEX
#ifndef MD_FOUR_BEADS
#ifdef MD_SUPERELLIPSOID
	  fprintf(fs, tipodat2_mgl,rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
		  uyz[i], uzx[i], uzy[i], uzz[i], typesArr[typeOfPart[i]].sax[0], 
		  typesArr[typeOfPart[i]].sax[1], typesArr[typeOfPart[i]].sax[2],
		  colsFlex[typeOfPart[i]%numcols], typesArr[typeOfPart[i]].n[0],
		  typesArr[typeOfPart[i]].n[1], typesArr[typeOfPart[i]].n[2]);
#else
      	  fprintf(fs, tipodat2_mgl,rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
		  uyz[i], uzx[i], uzy[i], uzz[i], typesArr[typeOfPart[i]].sax[0], 
		  typesArr[typeOfPart[i]].sax[1], typesArr[typeOfPart[i]].sax[2],
		  colsFlex[typeOfPart[i]%numcols]);
#endif
#endif
#else
	  if (i < Oparams.parnumA)
	    {
	      fprintf(fs, tipodat2_mgl,rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
		      uyz[i], uzx[i], uzy[i], uzz[i], Oparams.a[0], Oparams.b[0], Oparams.c[0],
		      "red");
	    }
	  else
	    {
	      fprintf(fs, tipodat2_mgl, rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
		      uyz[i], uzx[i], uzy[i], uzz[i],
		      Oparams.a[1], Oparams.b[1], Oparams.c[1],
		      "green");
	    }
#endif
#ifdef MD_PATCHY_HE
	  rA[0] = rx[i];
	  rA[1] = ry[i];
	  rA[2] = rz[i];
	  BuildAtomPos(i, rA, R[i], ratA);
#ifdef EDHE_FLEX
#ifdef MD_FOUR_BEADS
	 for (nn = 0; nn < 16; nn++)
	   {
	     nn2 = nn % 4;
	     if (nn2 != 0)
	       continue;
	     if (i%2==0)
	       strcpy(bcolor,beadcol[0][nn/4]);
	     else
	       strcpy(bcolor,beadcol[1][nn/4]);
	     fprintf(fs,"%.15f %.15f %.15f @ %.15G C[%s]\n", 
	 	     ratA[nn+1][0], ratA[nn+1][1], ratA[nn+1][2], typesArr[typeOfPart[i]].spots[nn].sigma*0.5,bcolor);
	   }
	 fprintf(fs, ".Bonds: 0-1[0.1:green],0-2[0.1:green],0-3[0.1:green],1-2[0.1:green],1-3[0.1:green],2-3[0.1:green]\n");
#else
	  for (nn = 1; nn < typesArr[typeOfPart[i]].nspots+1; nn++)
	    fprintf(fs,"%.15f %.15f %.15f @ %.15G C[orange]\n", 
		    ratA[nn][0], ratA[nn][1], ratA[nn][2], typesArr[typeOfPart[i]].spots[nn-1].sigma*0.5);
#endif
#else
	  for (nn = 1; nn < ((i < Oparams.parnumA)?MD_STSPOTS_A+1:MD_STSPOTS_B+1); nn++)
	    fprintf(fs,"%.15f %.15f %.15f @ %.15G C[orange]\n", 
		    ratA[nn][0], ratA[nn][1], ratA[nn][2], Oparams.sigmaSticky*0.5);
#endif
#endif
#ifdef MD_FOUR_BEADS
	  if (i < Oparams.parnum-1)
	    fprintf(fs, ".newmol\n");
#endif
	}
    }
  else
    {
#ifdef EDHE_FLEX
      if (!OprogStatus.stripStore)
	fprintf(fs, "@@@\n");
#endif
      for (i = 0; i < Oparams.parnum; i++)
	{
#ifdef EDHE_FLEX
	  if (!saveAll && !globSaveAll)
	    {
	      if (!is_to_save(i))
		continue;
	    }
#endif
#ifdef EDHE_FLEX
	  fprintf(fs, tipodat2_flex, rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
		  uyz[i], uzx[i], uzy[i], uzz[i], typeOfPart[i]);
#else
	  fprintf(fs, tipodat2, rx[i], ry[i], rz[i], uxx[i], uxy[i], uxz[i], uyx[i], uyy[i], 
		  uyz[i], uzx[i], uzy[i], uzz[i]);
#endif
	}
    }
  if (mgl_mode==0)
    {
      for (i = 0; i < Oparams.parnum; i++)
	{
#ifdef EDHE_FLEX
	  if (!saveAll && !globSaveAll)
	    {
	      if (!is_to_save(i))
		continue;
	    }
#endif
#ifdef MD_ASYM_ITENS
	  fprintf(fs, tipodat, vx[i], vy[i], vz[i], Mx[i], My[i], Mz[i]);
#else
	  fprintf(fs, tipodat, vx[i], vy[i], vz[i], wx[i], wy[i], wz[i]);
#endif
	}
#ifdef MD_POLYDISP
      for (i = 0; i < Oparams.parnum; i++)
	{
	  fprintf(fs, "%.15G %.15G %.15G\n", axaP[i], axbP[i], axcP[i]);
	}
#endif
#ifdef MD_LXYZ
      fprintf(fs, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
#else
#ifdef MD_GRAVITY
      fprintf(fs, "%.15G %.15G\n", L, Lz);
#else
      fprintf(fs, "%.15G\n", L);
#endif
#endif
    }
}
#ifdef EDHE_FLEX
int readBinCoord_heflex(int cfd)
{
  int i;
  int size;
  unsigned char rerr = 0;
#ifdef EDHE_FLEX 
  int sizeSPNNL;
#endif
#ifdef MD_GHOST_IGG
  if (Oparams.ghostsim)
    {
      size = sizeof(ghostInfo)*Oparams.parnum;
#ifdef MD_USE_CALLOC
      ghostInfoArr = calloc(Oparams.parnum, sizeof(ghostInfo));
#else
      ghostInfoArr = malloc(sizeof(ghostInfo)*Oparams.parnum);
#endif
      rerr |= readSegs(cfd, "Init", "Error reading ghostInfoArr", CONT, size, ghostInfoArr, NULL);
    }
#endif
  //printf("status:%d\n", ghostInfoArr[Oparams.parnum-1].ghost_status);
  size = sizeof(int)*Oparams.parnum;
  typeOfPart = malloc(size);
  rerr |= -readSegs(cfd, "Init", "Error reading typeNP", CONT, size, typeOfPart, NULL);

  size = sizeof(int)*Oparams.ntypes;
  typeNP = malloc(size);
  rerr |= -readSegs(cfd, "Init", "Error reading typeNP", CONT, size, typeNP, NULL);
 

  size = sizeof(partType)*Oparams.ntypes; 
#ifdef MD_USE_CALLOC
  typesArr = calloc(Oparams.ntypes, sizeof(partType));
#else
  typesArr = malloc(size);
#endif
  rerr |= -readSegs(cfd, "Init", "Error reading typesArr", CONT, size, typesArr, NULL);
  for (i=0; i < Oparams.ntypes; i++)
    {
      size = sizeof(spotStruct)*typesArr[i].nspots;
      if (typesArr[i].nspots >= NA)
	{
	  printf("[ERROR] too many spots (%d) for type %d increase NA (actual value is %d) in ellipsoid.h and recompile\n",
		 typesArr[i].nspots, i, NA);
	  exit(-1);
	}

#ifdef MD_SUPERELLIPSOID
     /*  16/05/2010: notare che nel caso si abbiano spot per le NNL bisogna allocare pi� spazio
	ma questi non vengono salvati, quindi il numero di spot letti non ne deve tener conto */ 
      sizeSPNNL = sizeof(spotStruct)*MD_SPNNL_NUMSP;
#ifdef MD_USE_CALLOC
      typesArr[i].spots = calloc(MD_SPNNL_NUMSP+typesArr[i].nspots,sizeof(spotStruct));
#else
      typesArr[i].spots = malloc(size+sizeSPNNL);
#endif
      if (typesArr[i].nspots == 0)
	continue;
#else
#if 1
      sizeSPNNL = sizeof(spotStruct)*MD_SPNNL_NUMSP;
#ifdef MD_USE_CALLOC
      typesArr[i].spots = calloc(MD_SPNNL_NUMSP+typesArr[i].nspots,sizeof(spotStruct));
#else
      typesArr[i].spots = malloc(size+sizeSPNNL);
#endif
      if (typesArr[i].nspots == 0)
	continue;
#else
      sizeSPNNL = 0;
      if (typesArr[i].nspots == 0)
	continue;
#ifdef MD_USE_CALLOC
      typesArr[i].spots = calloc(typesArr[i].nspots,sizeof(spotStruct));
#else
      typesArr[i].spots = malloc(size+sizeSPNNL);
#endif
#endif
#endif
      rerr |= -readSegs(cfd, "Init", "Error reading spots", CONT, size, typesArr[i].spots, NULL);
    } 
  /* read interactions */
  if (Oparams.ninters > 0)
    {
#ifdef MD_USE_CALLOC
      intersArr = calloc(Oparams.ninters, sizeof(interStruct));
#else
      intersArr = malloc(sizeof(interStruct)*Oparams.ninters);
#endif
      size = sizeof(interStruct)*Oparams.ninters;
      rerr |= -readSegs(cfd, "Init", "Error reading intersArr", CONT, size, intersArr, NULL);
      for (i=0; i < Oparams.ninters; i++)
	{
	  if (intersArr[i].type1==-2)
	    {
	      size = sizeof(rangeStruct)*intersArr[i].nr1;
	      intersArr[i].r1 = malloc(size);
	      rerr |= readSegs(cfd, "Init", "Error reading ranges", CONT, size, intersArr[i].r1, NULL);
	    }
	  if (intersArr[i].type2==-2)
	    {
	      size = sizeof(rangeStruct)*intersArr[i].nr2;
	      intersArr[i].r2 = malloc(size);
	      rerr |=readSegs(cfd, "Init", "Error reading ranges", CONT, size, intersArr[i].r2, NULL);
	    }
	}
    }
  if (Oparams.nintersIJ > 0)
    {
#ifdef MD_USE_CALLOC
      intersArrIJ = calloc(Oparams.nintersIJ, sizeof(interStructIJ));
#else
      intersArrIJ = malloc(sizeof(interStructIJ)*Oparams.nintersIJ);
#endif
      size = sizeof(interStructIJ)*Oparams.nintersIJ;
      rerr |= -readSegs(cfd, "Init", "Error reading intersArrIJ", CONT, size, intersArrIJ, NULL);
      for (i=0; i < Oparams.nintersIJ; i++)
	{
	  if (intersArrIJ[i].i==-2)
	    {
	      size = sizeof(rangeStruct)*intersArrIJ[i].nr1;
	      intersArrIJ[i].r1 = malloc(size);
	      rerr |= readSegs(cfd, "Init", "Error reading ranges", CONT, size, intersArrIJ[i].r1, NULL);
	    }
	  if (intersArrIJ[i].j==-2)
	    {
	      size = sizeof(rangeStruct)*intersArrIJ[i].nr2;
	      intersArrIJ[i].r2 = malloc(size);
	      rerr |=readSegs(cfd, "Init", "Error reading ranges", CONT, size, intersArrIJ[i].r2, NULL);
	    }
	}
    }
  if (Oparams.saveBonds)
    { 
      /* N.B. i file binari d'ora (07/06/07) in poi  conterranno sempre i bonds se Oparams.saveBonds=1! */
      size = Oparams.parnum*sizeof(int);
      numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
      rerr |= -readSegs(cfd, "Init", "Error reading numbonds", CONT, size, numbonds, NULL);
     
      if (OprogStatus.maxbonds < Oparams.maxbondsSaved)
	OprogStatus.maxbonds = Oparams.maxbondsSaved;
#ifdef MD_LL_BONDS
      bonds = AllocMatLLI(Oparams.parnum, OprogStatus.maxbonds);
#else
      bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
#endif
#ifdef MD_SPHERICAL_WALL
      allocBondsSphWall();
#endif
      for (i=0; i < Oparams.parnum; i++)
	{
#ifdef MD_LL_BONDS
	  size = sizeof(long long int)*numbonds[i];
#else
	  size = sizeof(int)*numbonds[i];
#endif
	  //printf("maxbondsSaved=%d numbonds[%d]:%d\n", Oparams.maxbondsSaved, i, numbonds[i]);
	  rerr |= -readSegs(cfd, "Init", "Error reading bonds", CONT, size, bonds[i], NULL);
	}
    }

  return rerr;
}
void writeBinCoord_heflex(int cfd)
{
  int i;
  int size;

#ifdef MD_GHOST_IGG
  if (Oparams.ghostsim)
    {
      size = sizeof(ghostInfo)*Oparams.parnum;
      writeSegs(cfd, "Init", "Error writing ghostInfoArr", CONT, size, ghostInfoArr, NULL);
    }
#endif

  size = sizeof(int)*Oparams.parnum;
  writeSegs(cfd, "Init", "Error writing typeOfPart", CONT, size, typeOfPart, NULL);

  size = sizeof(int)*Oparams.ntypes;
  writeSegs(cfd, "Init", "Error writing typeNP", CONT, size, typeNP, NULL);

  size = sizeof(partType)*Oparams.ntypes; 
  writeSegs(cfd, "Init", "Error writing typesArr", CONT, size, typesArr, NULL);

  for (i=0; i < Oparams.ntypes; i++)
    {
      if (typesArr[i].nspots == 0)
	continue;
      size = sizeof(spotStruct)*typesArr[i].nspots;
      writeSegs(cfd, "Init", "Error writing spots", CONT, size, typesArr[i].spots, NULL);
    } 
  if (Oparams.ninters > 0)
    {
      size = sizeof(interStruct)*Oparams.ninters;
      /* write interactions */
      writeSegs(cfd, "Init", "Error writing intersArr", CONT, size, intersArr, NULL);
      /* writing all ranges here */
      for (i=0; i < Oparams.ninters; i++)
	{
	  if (intersArr[i].type1==-2)
	    {
	      size = sizeof(rangeStruct)*intersArr[i].nr1;
	      writeSegs(cfd, "Init", "Error writing ranges", CONT, size, intersArr[i].r1, NULL);
	    }
	  if (intersArr[i].type2==-2)
	    {
	      size = sizeof(rangeStruct)*intersArr[i].nr2;
	      writeSegs(cfd, "Init", "Error writing ranges", CONT, size, intersArr[i].r2, NULL);
	    }
	}
    }
  if (Oparams.nintersIJ > 0)
    {
      size = sizeof(interStructIJ)*Oparams.nintersIJ;
      /* write interactions */
      writeSegs(cfd, "Init", "Error writing intersArrIJ", CONT, size, intersArrIJ, NULL);
      /* writing all ranges here */
      for (i=0; i < Oparams.nintersIJ; i++)
	{
	  if (intersArrIJ[i].i==-2)
	    {
	      size = sizeof(rangeStruct)*intersArrIJ[i].nr1;
	      writeSegs(cfd, "Init", "Error writing ranges", CONT, size, intersArrIJ[i].r1, NULL);
	    }
	  if (intersArrIJ[i].j==-2)
	    {
	      size = sizeof(rangeStruct)*intersArrIJ[i].nr2;
	      writeSegs(cfd, "Init", "Error writing ranges", CONT, size, intersArrIJ[i].r2, NULL);
	    }
	}
    }
  /* se Oparams.maxbondsSaved � > 0 vuol dire che sono stati salvati nella presente
     configurazione anche i bond */
  if (Oparams.saveBonds)
    { 
      /* N.B. i file binari d'ora (07/06/07) in poi  conterranno sempre i bonds se Oparams.saveBonds=1! */
      size = Oparams.parnum*sizeof(int);
      writeSegs(cfd, "Init", "Error writing numbonds", CONT, size, numbonds, NULL);
      for (i=0; i < Oparams.parnum; i++)
	{
#ifdef MD_LL_BONDS
	  size = sizeof(long long int)*numbonds[i];
#else
	  size = sizeof(int)*numbonds[i];
#endif
	  //printf("writing bonds: %d\n", numbonds[i]);
	  writeSegs(cfd, "Init", "Error writing bonds", CONT, size, bonds[i], NULL);
	}
    }
}

#endif
#ifdef EDHE_FLEX
void parse_one_range(char *s, int *A, int *nr, rangeStruct **r)
{
  int m, M;
  if (sscanf(s, "%d-%d", &m, &M)==2)
    {
      *A = -2;
      if (*nr == 0)
	*r = malloc(sizeof(rangeStruct)*(*nr+1));
      else
	*r = realloc(*r, sizeof(rangeStruct)*(*nr+1));
      (*r)[*nr].min = m;
      (*r)[*nr].max = M;
      (*nr)++;
    }
  else
    {
      *A = atoi(s);
      *r = NULL;
      *nr = 0;
    }
}
void parse_ranges(char *s, int *A, int *nr, rangeStruct **r)
{
  char *ns;
  ns = strtok(s, ",");

  *nr = 0;
  if (!ns)
    {
      parse_one_range(s, A, nr, r);
      return;
    }
  while (ns)
    {
      MD_DEBUG31(printf("ns=%s\n", ns));
      parse_one_range(ns, A, nr, r);
      ns = strtok(NULL, ",");
    } 
}
#endif
char line[4096];
/* ========================== >>> readAllCor <<< ========================== */

void readAllCor(FILE* fs)
{
  int i;
#ifdef EDHE_FLEX
#ifdef MD_GHOST_IGG
  int size;
#endif
#ifdef MD_MULTIPLE_LL
  int beg;
#endif
  char sep[256];
  int j;
  char *s1, *s2;
  int sizeSPNNL;

  s1 = malloc(sizeof(char)*65535);
  s2 = malloc(sizeof(char)*65535);
  typeOfPart = malloc(sizeof(int)*Oparams.parnum);
  typeNP = malloc(sizeof(int)*Oparams.ntypes);

#ifdef MD_USE_CALLOC
  typesArr = calloc(Oparams.ntypes, sizeof(partType));
#else
  typesArr = malloc(sizeof(partType)*Oparams.ntypes);
#endif
#ifdef MD_MULTIPLE_LL
  fscanf(fs, "%s ", line);
  if (!strcmp(line, "RF"))
    { 
      //printf("line=%s\n", line);
      beg=0;
      for (i=0; i < Oparams.ntypes; i++)
	{
	  fscanf(fs, "%lf ", &(typesArr[i].rcutFact));
	  printf("typeArr[%d].rcutFact=%f\n", i, typesArr[i].rcutFact);
	}
    }
  else
    {
      sscanf(line, "%d", &typeNP[0]);
      beg=1;
      for (i=0; i < Oparams.ntypes; i++)
	typesArr[i].rcutFact = -1.0;
 
    }
  for (i=beg; i < Oparams.ntypes; i++)
    {
      fscanf(fs, "%d ", &typeNP[i]);
      //printf("typeNP[%d]=%d\n", i, typeNP[i]);
    }

#else
  for (i=0; i < Oparams.ntypes; i++)
    {
      fscanf(fs, "%d ", &typeNP[i]);
      //printf("typeNP[%d]=%d\n", i, typeNP[i]);
    }  
#endif
    
  for (i=0; i < Oparams.ntypes; i++)
    {
        /* read particles parameters */
      fscanf(fs, "%lf %lf %lf ", &typesArr[i].sax[0], &typesArr[i].sax[1], &typesArr[i].sax[2]); 
      fscanf(fs, "%lf %lf %lf ", &typesArr[i].n[0], &typesArr[i].n[1], &typesArr[i].n[2]);
      fscanf(fs, "%lf %lf %lf %lf %d %d ", &typesArr[i].m, &typesArr[i].I[0], &typesArr[i].I[1],
	   &typesArr[i].I[2], &typesArr[i].brownian, &typesArr[i].ignoreCore);
#if 0
      fscanf(fs, "%lf %lf %lf ", &typesArr[i].xoff[0], &typesArr[i].xoff[1], &typesArr[i].xoff[2]); 
#endif
      /* read sticky spots parameters */
      fscanf(fs, "%d %d ", &typesArr[i].nspots, &typesArr[i].nhardobjs);
      if (typesArr[i].nspots >= NA)
	{
	  printf("[ERROR] too many spots (%d) for type %d increase NA (actual value is %d) in ellipsoid.h and recompile\n",
		 typesArr[i].nspots, i, NA);
	  exit(-1);
	}
#ifdef MD_SUPERELLIPSOID
      sizeSPNNL = sizeof(spotStruct)*MD_SPNNL_NUMSP;
#else
      sizeSPNNL = sizeof(spotStruct)*MD_SPNNL_NUMSP;
      //sizeSPNNL = 0;
#endif
      if (sizeSPNNL > 0 || typesArr[i].nspots > 0)
	{
#ifdef MD_USE_CALLOC
	  typesArr[i].spots = calloc(typesArr[i].nspots+MD_SPNNL_NUMSP,sizeof(spotStruct));
#else
	  typesArr[i].spots = malloc(sizeof(spotStruct)*typesArr[i].nspots+sizeSPNNL);
#endif
	}
      else
	typesArr[i].spots = NULL;
      //printf("QUI nhardobjs=%d ntypes=%d\n",typesArr[i].nhardobjs, Oparams.ntypes );
      for (j = 0; j < typesArr[i].nspots; j++)
	fscanf(fs, "%lf %lf %lf %lf ", &typesArr[i].spots[j].x[0],&typesArr[i].spots[j].x[1],
	       &typesArr[i].spots[j].x[2], &typesArr[i].spots[j].sigma);
      /* hard objects (for now super-ellipsoids) */
      if (typesArr[i].nhardobjs > 0)
	{
#ifdef MD_USE_CALLOC
	  typesArr[i].hardobjs = calloc(typesArr[i].nhardobjs,sizeof(hardobjsStruct));
#else
	  typesArr[i].hardobjs = malloc(sizeof(hardobjsStruct)*typesArr[i].nhardobjs);
#endif
	}
      else
	typesArr[i].hardobjs = NULL;
	  
      for (j = 0; j < typesArr[i].nhardobjs; j++)
	fscanf(fs, "%lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
	       &typesArr[i].hardobjs[j].x[0],&typesArr[i].hardobjs[j].x[1], &typesArr[i].hardobjs[j].x[2], 
	       &typesArr[i].hardobjs[j].sax[0],&typesArr[i].hardobjs[j].sax[1], &typesArr[i].hardobjs[j].sax[2],
	       &typesArr[i].hardobjs[j].n[0], &typesArr[i].hardobjs[j].n[1], &typesArr[i].hardobjs[j].n[2]);
    } 
  //printf("qui nintersIJ=%d\n", Oparams.nintersIJ);
  /* read interactions */
  if (Oparams.ninters > 0)
    {
#ifdef MD_USE_CALLOC
      intersArr = calloc(Oparams.ninters, sizeof(interStruct));
#else
      intersArr = malloc(sizeof(interStruct)*Oparams.ninters);
#endif
    }
  else
    intersArr = NULL;
#if 0
  for (i=0; i < Oparams.ninters; i++)
   {
     fscanf(fs, "%d %d %d %d %lf %lf %lf %d ", &intersArr[i].type1, &intersArr[i].spot1, &intersArr[i].type2, 
	    &intersArr[i].spot2, 
	    &intersArr[i].bheight, &intersArr[i].bhin, &intersArr[i].bhout, &intersArr[i].nmax);
   } 
#else
    
  for (i=0; i < Oparams.ninters; i++)
    {
      fscanf(fs, "%s %d %s %d %lf %lf %lf %d ", s1, &intersArr[i].spot1, s2, 
 	     &intersArr[i].spot2, 
 	     &intersArr[i].bheight, &intersArr[i].bhin, &intersArr[i].bhout, &intersArr[i].nmax);
      parse_ranges(s1, &intersArr[i].type1, &intersArr[i].nr1, &intersArr[i].r1);
      parse_ranges(s2, &intersArr[i].type2, &intersArr[i].nr2, &intersArr[i].r2);
      MD_DEBUG31(printf("type1=%d type2=%d\n", intersArr[i].type1, intersArr[i].type2));
    } 
#endif
  
  /* read interactions */
  if (Oparams.nintersIJ > 0)
    {
#ifdef MD_USE_CALLOC
      intersArrIJ = calloc(Oparams.nintersIJ, sizeof(interStructIJ));
#else
      intersArrIJ = malloc(sizeof(interStructIJ)*Oparams.nintersIJ);
#endif
      for (i=0; i < Oparams.nintersIJ; i++)
	{
#if 0
	  fscanf(fs, "%d %d %d %d %lf %lf %lf ", &intersArrIJ[i].i, &intersArrIJ[i].spot1, &intersArrIJ[i].j, 
		 &intersArrIJ[i].spot2, 
		 &intersArrIJ[i].bheight, &intersArrIJ[i].bhin, &intersArrIJ[i].bhout);
#else
	  fscanf(fs, "%s %d %s %d %lf %lf %lf ", s1, &intersArrIJ[i].spot1, s2, 
		 &intersArrIJ[i].spot2, 
		 &intersArrIJ[i].bheight, &intersArrIJ[i].bhin, &intersArrIJ[i].bhout);
	  parse_ranges(s1, &intersArrIJ[i].i, &intersArrIJ[i].nr1, &intersArrIJ[i].r1);
	  parse_ranges(s2, &intersArrIJ[i].j, &intersArrIJ[i].nr2, &intersArrIJ[i].r2);
	  MD_DEBUG31(printf("i=%d j=%d\n", intersArrIJ[i].i, intersArrIJ[i].j));
#if 0
	  if (intersArrIJ[i].j==-2)
	    {
	      int kk;
	      for (kk = 0; kk < intersArrIJ[i].nr2; kk++)
		printf(">>> %d-%d\n", intersArrIJ[i].r2[kk].min, intersArrIJ[i].r2[kk].max);
	    }
#endif
#endif
	} 
    }
  fscanf(fs, "%s ", sep);
  /* se ci sono i bond li legge, altrimenti ci deve essere il separatore "@@@" */
  if (strcmp(sep, "@@@"))
    {
      /* Oparams.maxbondsSaved deve contenere il numero massimo di bond salvati 
	 nella configurazione che si sta leggendo */ 
      if (Oparams.maxbondsSaved < 0)
	Oparams.maxbondsSaved = OprogStatus.maxbonds;
      if (OprogStatus.maxbonds < Oparams.maxbondsSaved)
	OprogStatus.maxbonds = Oparams.maxbondsSaved;
#ifdef MD_LL_BONDS
      bonds = AllocMatLLI(Oparams.parnum, OprogStatus.maxbonds);
#else
      bonds = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
#endif
      numbonds = (int *) malloc(Oparams.parnum*sizeof(int));
#ifdef MD_SPHERICAL_WALL
      allocBondsSphWall();
#endif
 
      for (i = 0; i < Oparams.parnum; i++)
	{
	  if (i==0)
	    sscanf(sep, "%d", &numbonds[0]);
	  else
	    fscanf(fs, "%d ", &numbonds[i]);
	  //printf("i=%d, numbonds[]=%d\n", i, numbonds[i]);
	  for (j = 0; j < numbonds[i]; j++)
	    {
#ifdef MD_LL_BONDS
	      fscanf(fs, "%lld ", &bonds[i][j]);
#else	      
	      fscanf(fs, "%d ", &bonds[i][j]);
#endif	    
	    }
	}
      fscanf(fs, "@@@ ");
    }
#if 0
#ifdef MD_GHOST_IGG
  /* for now ascii restore file does not contain ghostInfoArr data, hence
     to restart a simulation preserving ghost status a binary COORD_TMPX file has to be used (ellipsoid -c)*/
  fscanf(fs, "%s ", sep);
  if (!strcmp(sep, "@@@ GHOST @@@"))
    {
      size = sizeof(ghostInfo)*Oparams.parnum;
      ghostInfoArr = malloc(sizeof(ghostInfo)*Oparams.parnum);
      for (i=0; i < Oparams.parnum; i++)
	fscanf(fs, "%d %d ", &(ghostInfoArr[i].iggnum), &(ghostInfoArr[i].ghost_status));
    }
#endif
#endif
#endif
  for (i = 0; i < Oparams.parnum; i++)
    {
      if (fscanf(fs, "%lf %lf %lf ", &rx[i], &ry[i], &rz[i]) < 3)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
#ifdef EDHE_FLEX
      if (fscanf(fs, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
		 &uxx[i], &uxy[i], &uxz[i], &uyx[i], &uyy[i], &uyz[i], &uzx[i], &uzy[i], &uzz[i],
		 &typeOfPart[i]) < 10)
	{
	  //printf("typeOfPart[%d]:%d\n", i, typeOfPart[i]);
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
      //printf("pos=%f %f %f type=%d\n", rx[i], ry[i], rz[i], typeOfPart[i]);
#else
      if (fscanf(fs, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		 &uxx[i], &uxy[i], &uxz[i], &uyx[i], &uyy[i], &uyz[i], &uzx[i], &uzy[i], &uzz[i]) < 9)
	{
	  mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	  exit(-1);
	}
#endif
      //printf("%d r=(%f,%f,%f) type=%d\n", i, rx[i], ry[i], rz[i], typeOfPart[i]);
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
      //printf("%d v=(%f,%f,%f)\n", i, vx[i], vy[i], vz[i]);
    }
#ifdef MD_POLYDISP
  for (i = 0; i < Oparams.parnum; i++)
    {
      fscanf(fs, "%lf %lf %lf\n", &axaP[i], &axbP[i], &axcP[i]);
    }
#endif
#ifdef MD_LXYZ
  fscanf(fs, "%[^\n]", line);

  if (sscanf(line, "%lf %lf %lf\n", &L[0], &L[1], &L[2]) < 3)
    {
      if (sscanf(line, "%lf\n",  &L[0]) == 1)
	{
	  L[1]=L[2]=L[0];
	}
      else
	{
	  mdPrintf(STD, "ERROR[extra] reading ascii file\n", NULL);
	  exit(-1);
	}
    }
#else
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
#endif      
#ifdef EDHE_FLEX
  free(s1);
  free(s2);
#endif
}
