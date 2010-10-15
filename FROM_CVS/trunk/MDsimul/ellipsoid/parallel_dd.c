#include<mdsimul.h>
#ifdef ED_PARALL_DD
/* ==== >>> routines to initialize structures <<< ==== */
extern int cellsx, cellsy, cellsz;
void calc_xyz(unsigned long long int cellnum, unsigned int *ix, unsigned int *iy, unsigned int *iz) 
{
  /* è la funzione inversa di calc_cellnum(...) */
  int i, k, km;
  const int maxk=60;
  unsigned int ni[3] = {0x0,0x0,0x0};
  unsigned long long int cn;
  cn = cellnum;
  /* km sarà il numero di bit non nulli di cellnum */
  for (k=0; k < maxk; k++)
    {
      if (cn==0x0)
	{
	  km = k;
	  break;
	}
      cn = cn >> 1;
    }
  for (k=0; k < km; k++)
    {
      ni[k % 3] |= (cellnum & 0x1) << (k/3);
      *iz |= (cellnum & 0x1) << k;
      cellnum = cellnum >> 1; 
    } 
  *ix = ni[0];
  *iy = ni[1];
  *iz = ni[2];
}
unsigned long long int calc_cellnum(int ix, int iy, int iz)
{
  /* N.B. interleaving cells coordinates in binary represantion a unique cell  
     number can be calculated (see Miller and Luding J. Comp. Phys. 2003) */
  const unsigned int maxk=20;
  unsigned int ixl, iyl, izl, km;
  unsigned long long int im=0x1, bitx, bity, bitz, res;
  /* max number of cells is 2^(3*20)*/
  ixl = ix;
  iyl = iy;
  izl = iz;
  if (dd_inCell[ix][iy][iz] != -1)
    return dd_inCell[ix][iy][iz];
  /* se n_bits(i) è il numero massimo di bits dell'unsigned int i
     allora di seguito si valuta max{n_bits(ix),n_bits(iy), n_bits(iz)}*/
  for (k=0; k < maxk; k++)
    {
      if (ixl==0x0 && iyl==0x0 && izl==0x0)
	{
	  km = k;
	  break;
	}
      ixl = ixl >> 1;
      iyl = iyl >> 1;
      izl = izl >> 1;
    }
  res = 0x0;
  for (k=0; k < km; k++)
    {
      bitx = ix & im;
      bity = iy & im;
      bitz = iz & im;
      res |= (bitz << (2*k+2)) | (bity << (2*k+1)) | (bitx << (2*k)); 
      im = im << 1;
    }
  dd_inCell[ix][iy][iz] = res;
  return res;
}
int cellToRegion(unsigned long long int cell)
{
  /* implementare una ricerca dicotomica */
  int i;
  for (i=0; i < numOfProcs; i++)
    {
      if (cell < regionsArr[i])
	return i;
    }
}
void rollback_init(void);
int *border_cells_head, *border_cells_ll, *dd_inCell;
int *neigh_regions_of_cell_list_head, *neigh_regions_of_cell_list;
enum {DD_REAL=0, DD_BORDER_ZONE, DD_VIRTUAL};
void calc_dd_inCell(int i, int ix, int iy, int iz)
{
  /* questa funzione va chiamata ogni volta che si costruiscono
     le linked cell list ed ogni volta che viene processato un 
     evento di cell-crossing */
  dd_inCell[i] = calc_cellnum(ix, iy, iz);
}  
void dd_init(void)
{
  /* queste dichiarazioni vanno poi rese globali */
  int *inRegion, dd_numreg, cc, iX, iY, iZ, ncp;
  /* 08/10/10: nel caso di codice parallelo è meglio non usare le multiple linked list poiche'
     il tutto si complica significativamente */
  /* numero totale di regioni (1 regione = 1 processo) */
  dd_numreg = numOfProcs; /* numOfProcs viene definita se si usa MPI */
  /* number of regions per process */
  ncp = (int) pow(2,(int)log2(cellsx*cellsy*cellsz/dd_numreg));
  if (OprogStatus.dd_numcells_per_proc == 0)
    {
      printf("ERROR: less than a cell per region, exiting...!\n");
      exit(-1);
    }
  /* inRegion[c] gives region associated to cell c */
  inRegion = malloc(sizeof(int)*cellsx*cellsy*cellsx);
  /* array contenente la fine di ogni regione nella sequenza associata alle celle
     tramite interleaving dei bit (vedi Sec. 3.1 S. Miller and S. Luding J. Comput. Phys. 193, 
     306-316 (2003) */
  regionsArr = malloc(sizeof(int)*dd_numreg);
  /* dd_inCell contiene la cella associata ad ogni particella (calcolata tramite interleaving come 
     spiegato da Luding and Miller, J. Comp. Phys. 2003) */
  dd_inCell = malloc(sizeof(int)*Oparams.parnum);

  border_cells_head = malloc(sizeof(int)*dd_numreg);
  border_cells_ll = malloc(sizeof(int)*cellsx*cellsy*cellsz);

  /* le regioni sono in generale dei parallelepipedi contenenti un numero intero 
     di celle ed individuati da 6 numeri (x1,x2), (y1,y2), (z1,z2) ossia x1i e x2i, all'inizio
     i parallelepipedi vengono scelti tutti uguali (a parte gli ultimi se il numero di regioni
     lungo un asse non è sottomultiplo del numero di celle. */
  for (i = 0; i < dd_numreg-1; i++)
    {
      /* regionsArr[i] contiene la fine della regione i-esima (in numero di celle), ossia
	 la prima cella che non gli appartiene. */
      regionsArr[i] = (i+1)*ncp;
    }
  /* all cells left assigned to last region */
  regionsArr[dd_numreg-1] = cellsx*cellsy*cellsz;
  if (regionsArr[dd_numreg-1] % 2 != 0)
    {
      printf("ERROR: total number of cell is %d, but it must be even in a parallel simulation\n", 
	     cellsx*cellsy*cellsz);
      exit(-1);
    }
  for (iZ = 0; iZ < cellsz; iZ++)
    for (iY = 0; iY < cellsy; iY++)
      for (iX = 0; iX < cellsx; iX++)
	{
	  cellnum = calc_cellnum(iX,iY,iZ);
	  inRegion[cellnum] = cellToRegion(cellnum);
	  cell_type[cellnum] = DD_REAL;
	}
  rollback_init();
}

/* ==== >>> rollback <<< ==== */
struct params rb_Oparams;
struct progStatus rb_OprogStatus;
extern double *lastcol, *atomTime, *lastbump;
extern int *inCell[3], **tree, *cellList;
extern double *treeRxC, *treeRyC, *treeRzC, *treeTime;
extern int *numbonds, **bonds, *oldTypeOfPart;
extern double **R, **RM;
extern struct nebrTabStruct *nebrTab;
extern ghostInfo *ghostInfoArr;
extern int *typeOfPart, *typeNP;
extern int *linearLists;
extern long long int itsF, timesF, itsS, timesS, numcoll, itsFNL, timesFNL, 
     timesSNL, itsSNL, numcalldist, numdisttryagain;
extern long long int itsfrprmn, callsfrprmn, callsok, callsprojonto, itsprojonto;
extern long long int accngA, accngB;
#ifdef MD_ASYM_ITENS
extern double *theta0, *phi0, *psi0, *costheta0, *sintheta0, *angM, ***RM;
#endif
extern double *axa, *axb, *axc;
extern double *a0I, *maxax;
extern int *scdone;

int dd_totBytes; 
void *dd_coord_ptr, *dd_coord_ptr_rb;
double *lastcol_rb, *atomTime_rb, *lastbump_rb;
int *inCell_rb[3], **tree_rb, *cellList_rb;
double *treeRxC_rb, *treeRyC_rb, *treeRzC_rb, *treeTime_rb;
int *numbonds_rb, **bonds_rb, *oldTypeOfPart_rb;
double **R_rb, **RM_rb;
struct nebrTabStruct *nebrTab_rb;
ghostInfo *ghostInfoArr_rb;
int *typeOfPart_rb, *typeNP_rb;
int *linearLists_rb;
long long int itsF_rb, timesF_rb, itsS_rb, timesS_rb, numcoll_rb, itsFNL_rb, timesFNL_rb, 
     timesSNL_rb, itsSNL_rb, numcalldist_rb, numdisttryagain_rb;
long long int itsfrprmn_rb, callsfrprmn_rb, callsok_rb, callsprojonto_rb, itsprojonto_rb;
long long int accngA_rb, accngB_rb;
#ifdef MD_ASYM_ITENS
double *theta0_rb, *phi0_rb, *psi0_rb, *costheta0_rb, *sintheta0_rb, *angM_rb, ***RM_rb;
#endif
double *axa_rb, *axb_rb, *axc_rb;
double *a0I_rb, *maxax_rb;
int *scdone_rb;
/* dd_coord_ptr viene inizializzato a seguito dell'allocazione
   della memoria per le coordinate, tramite la funzione mdarray.c:AllocCoord().
   Invece dd_coord_ptr_rb è un puntatore alla memoria usata per il memorizzare 
   i dati puntati da dd_coord_ptr. */
#ifdef MD_SPHERICAL_WALL
void allocBondsSphWall_rb(void)
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
	  bonds_rb[i] = malloc(sizeof(long long int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#else
	  bonds_rb[i] = malloc(sizeof(int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#endif
	  //break;
	}
    }
}
#endif

void rollback_init()
{
  int SEGSIZE, poolSize;
  /* allocate memory for saving rollback data */
#ifdef MD_DYNAMIC_OPROG
  rb_OprogStatus.ptr = malloc(OprogStatus.len);
#endif 
  dd_coord_ptr_rb = malloc(dd_totBytes); 
  lastcol_rb= malloc(sizeof(double)*Oparams.parnum);
  atomTime_rb = malloc(sizeof(double)*Oparams.parnum);
#ifdef MD_PATCHY_HE
  lastbump_rb =  malloc(sizeof(struct LastBumpS)*Oparams.parnum);
#else
  lastbump_rb = malloc(sizeof(int)*Oparams.parnum);
#endif
  cellList_rb = malloc(sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  inCell_rb[0] = malloc(sizeof(int)*Oparams.parnum);
  inCell_rb[1] = malloc(sizeof(int)*Oparams.parnum);
  inCell_rb[2] = malloc(sizeof(int)*Oparams.parnum);
#ifdef MD_LL_BONDS
  bonds_rb = AllocMatLLI(Oparams.parnum, OprogStatus.maxbonds);
#else
  bonds_rb = AllocMatI(Oparams.parnum, OprogStatus.maxbonds);
#endif
  numbonds_rb = (int *) malloc(Oparams.parnum*sizeof(int));
#ifdef MD_SPHERICAL_WALL
  allocBondsSphWall_rb();
#endif
#ifdef MD_SPHERICAL_WALL
  poolSize = OprogStatus.eventMult*Oparams.parnum+2*Oparams.parnum;
#else
  poolSize = OprogStatus.eventMult*Oparams.parnum;
#endif
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
#ifdef MD_CALENDAR_HYBRID
  tree_rb = AllocMatI(16, poolSize);
#else
  tree_rb = AllocMatI(13, poolSize);
#endif
#else
#ifdef MD_CALENDAR_HYBRID
  tree_rb = AllocMatI(13, poolSize);
#else
  tree_rb = AllocMatI(10, poolSize);
#endif
#endif
  treeTime_rb = malloc(sizeof(double)*poolSize);
  treeRxC_rb  = malloc(sizeof(double)*poolSize);
  treeRyC_rb  = malloc(sizeof(double)*poolSize);
  treeRzC_rb  = malloc(sizeof(double)*poolSize);
#ifdef MD_ABSORP_POLY
  oldTypeOfPart_rb = malloc(sizeof(int)*Oparams.parnum);
#endif
#ifdef MD_GHOST_IGG
  ghostInfoArr_rb = malloc(sizeof(ghostInfo)*Oparams.parnum);
#endif
  typeOfPart_rb = malloc(sizeof(int)*Oparams.parnum);
  typeNP_rb=malloc(sizeof(int)*Oparams.ntypes);
  if (OprogStatus.useNNL)
    {  
      nebrTab_rb = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
      for (i=0; i < Oparams.parnum; i++)
	{
    	  nebrTab_rb[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	}
    }
  RM_rb = malloc(sizeof(double**)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++) 
    {
#ifdef MD_MATRIX_CONTIGOUS
      /* alloca R in maniera contigua */
      if (i==0)
	{
  	  RM_rb[i] = malloc(sizeof(double*)*3);
	  RM_rb[i][0] = malloc(sizeof(double)*Oparams.parnum*9);
	  RM_rb[i][1] = RM[i][0] + 3;
	  RM_rb[i][2] = RM[i][1] + 3;
	}
      else
	{
	  RM_rb[i] = malloc(sizeof(double*)*3);
	  RM_rb[i][0] = RM_rb[i-1][2] + 3;
	  RM_rb[i][1] = RM_rb[i][0] + 3;
	  RM_rb[i][2] = RM_rb[i][1] + 3;
	}
#else
      RM_rb[i] = matrix(3, 3);
#endif
    }
#endif
  R_rb = malloc(sizeof(double**)*Oparams.parnum);
  
  for (i=0; i < Oparams.parnum; i++)
    {
#ifdef MD_MATRIX_CONTIGOUS
      /* alloca R in maniera contigua */
      if (i==0)
	{
  	  R_rb[i] = malloc(sizeof(double*)*3);
	  R_rb[i][0] = malloc(sizeof(double)*Oparams.parnum*9);
	  R_rb[i][1] = R_rb[i][0] + 3;
	  R_rb[i][2] = R_rb[i][1] + 3;
	}
      else
	{
	  R_rb[i] = malloc(sizeof(double*)*3);
	  R_rb[i][0] = R_rb[i-1][2] + 3;
	  R_rb[i][1] = R_rb[i][0] + 3;
	  R_rb[i][2] = R_rb[i][1] + 3;
	}
#else
      R_rb[i] = matrix(3, 3);
#endif
    }
#ifdef MD_ASYM_ITENS
  costheta0_rb = malloc(sizeof(double)*Oparams.parnum);
  sintheta0_rb = malloc(sizeof(double)*Oparams.parnum);
  theta0_rb =    malloc(sizeof(double)*Oparams.parnum);
  psi0_rb   =    malloc(sizeof(double)*Oparams.parnum);
  phi0_rb   =    malloc(sizeof(double)*Oparams.parnum);
  angM_rb   =    malloc(sizeof(double)*Oparams.parnum);
#endif
  axa_rb = malloc(sizeof(double)*Oparams.parnum);
  axb_rb = malloc(sizeof(double)*Oparams.parnum);
  axc_rb = malloc(sizeof(double)*Oparams.parnum);
  /* these array are for growth simulations and in this implementation
     growth is not allowed in parallel... (?!?) */
  a0I_rb = malloc(sizeof(double)*Oparams.parnum);
  maxax_rb = malloc(sizeof(double)*Oparams.parnum);
  scdone_rb = malloc(sizeof(int)*Oparams.parnum);
#ifdef MD_CALENDAR_HYBRID
  linearLists_rb = malloc(sizeof(int)*(OprogStatus.nlistsHQ+1));
#endif
}

void rollback_save(void)
{
  int dblparnum, intparnum, i;
  dblparnum = sizeof(double)*Oparams.parnum;
  intparnum = sizeof(int)*Oparams.parnum;
  memcpy(rb_Oparams,Oparams,sizeof(struct params));
  memcpy(rb_OprogStatus,OprogStatus,sizeof(struct progStatus));
  memcpy(rb_OprogStatus.ptr,OprogStatus.ptr, OprogStatus.len);
  memcpy(dd_coord_ptr_rb,dd_coord_ptr,dd_totBytes);
  memcpy(lastcol_rb,lastcol,dblparnum);
  memcpy(atomTime_rb,atomTime,dblparnum);
#ifdef MD_PATCHY_HE
  memcpy(lastbump_rb,lastbump,sizeof(struct LastBumpS)*Oparams.parnum);
#else
  memcpy(lastbump_rb,lastbump,intparnum);
#endif
  memcpy(cellList_rb, cellList, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  memcpy(inCell_rb[0], inCell[0], intparnum);
  memcpy(inCell_rb[1], inCell[1], intparnum);
  memcpy(inCell_rb[2], inCell[2], intparnum);
#ifdef MD_LL_BONDS
  memcpy(bonds_rb[0], bonds[0], sizeof(long long int)*Oparams.parnum*OprogStatus.maxbonds);
#else
  memcpy(bonds_rb[0], bonds[0], intparnum*OprogStatus.maxbonds);
#endif
  memcpy(numbonds_rb, numbonds, intparnum);
#ifdef MD_SPHERICAL_WALL
  for (i=0; i < Oparams.parnum; i++)
    {
      /* gli ultimi due tipi devono essere i "muri" sferici */
      if (typeOfPart[i]==Oparams.ntypes-1 || typeOfPart[i]==Oparams.ntypes-2)
	{
#ifdef MD_LL_BONDS
	  /* NOTA 21/04/2010: il 3 l'ho messo per tener conto del fatto che nel caso ad esempio 
	  con l'interazione SW i legami possono essere anche due per particella. */
	  memcpy(bonds_rb[i], bonds[i], sizeof(long long int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#else
	  memcpy(bonds_rb[i], bonds[i], intparnum*MD_MAX_BOND_PER_PART);
#endif
	}
    }
#endif
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
#ifdef MD_CALENDAR_HYBRID
  memcpy(tree_rb[0], tree[0], 16*poolSize*sizeof(int));
#else
  memcpy(tree_rb[0], tree[0], 13*poolSize*sizeof(int));
#endif
#else
#ifdef MD_CALENDAR_HYBRID
  memcpy(tree_rb[0], tree[0], 13*poolSize*sizeof(int));
#else
  memcpy(tree_rb[0], tree[0], 10*poolSize*sizeof(int));
#endif
#endif
  memcpy(treeTime_rb, treeTime, sizeof(double)*poolSize);
  memcpy(treeRxC_rb, treeRxC,   sizeof(double)*poolSize);
  memcpy(treeRyC_rb, treeRyC,   sizeof(double)*poolSize);
  memcpy(treeRzC_rb, treeRzC,   sizeof(double)*poolSize);
#ifdef MD_ABSORP_POLY
  memcpy(oldTypeOfPart_rb, oldTypeOfPart, intparnum);
#endif
#ifdef MD_GHOST_IGG
  memcpy(ghostInfoArr_rb, ghostInfoArr, sizeof(ghostInfo)*Oparams.parnum);
#endif
  memcpy(typeOfPart_rb, typeOfPart, intparnum);
  memcpy(typeNP_rb, typeNP, sizeof(int)*Oparams.ntypes);
  if (OprogStatus.useNNL)
    {  
      memcpy(nebrTab_rb, nebrTab, sizeof(struct nebrTabStruct)*Oparams.parnum);
      for (i=0; i < Oparams.parnum; i++)
	{
    	  memcpy(nebrTab_rb[i].list, nebrTab[i].list, sizeof(int)*nebrTab[i].len);
	}
    }
#ifdef MD_MATRIX_CONTIGOUS
  memcpy(RM_rb[0], RM[0], dblparnum*9);
#else
  for (i=0; i < Oparams.parnum; i++) 
    {
      memcpy(RM_rb[i], RM[i], 9*sizeof(double));
    }
#endif
#ifdef MD_MATRIX_CONTIGOUS
  memcpy(R_rb[0], R[0], dblparnum*9);
#else
  for (i=0; i < Oparams.parnum; i++) 
    {
      memcpy(R_rb[i], R[i], 9*sizeof(double));
    }
#endif
#ifdef MD_ASYM_ITENS
  memcpy(costheta0_rb, costhet0, dblparnum);
  memcpy(sintheta0_rb, sintheta0, dblparnum);
  memcpy(theta0_rb, theta0, dblparnum);
  memcpy(psi0_rb, psi0, dblparnum);
  memcpy(phi0_rb, phi0, dblparnum);
  memcpy(angM_rb, angM, dblparnum);
#endif
  memcpy(axa_rb, axa, dblparnum);
  memcpy(axb_rb, axb, dblparnum);
  memcpy(axc_rb, axc, dblparnum);
  if (OprogStatus.targetPhi > 0.0)
    {
      memcpy(a0I_rb, a0I, dblparnum);
      memcpy(maxax_rb, maxax, dblparnum);
      memcpy(scdone_rb, scdone, intparnum);
    }
#ifdef MD_CALENDAR_HYBRID
  memcpy(linearLists_rb, linearLists, sizeof(int)*(OprogStatus.nlistsHQ+1));
#endif
  itsF_rb = itsF;
  timeF_rb = timesF;
  itsS_rb = itsS;
  timesS_rb =  timesS;
  numcoll_rb = numcoll;
  itsFNL_rb =  itsFNL;
  timesFNL_rb = timesFNL; 
  timesSNL_rb = timesSNL;
  itsSNL_rb = itsSNL; 
  numcalldist_rb = numcalldist; 
  numdisttryagain_rb = numdisttryagain;
  itsfrprmn_rb = itsfrprmn;
  callsfrprmn_rb = callsfrprmn;
  callsok_rb = callsok;
  callsprojonto_rb = callsprojonto;
  itsprojonto_rb = itsprojonto;
  accngA_rb = accngA;
  accngB_rb = accngB;
}

void rollback_load(void)
{
  int dblparnum, intparnum, i;
  memcpy(Oparams,rb_Oparams,sizeof(struct params));
  memcpy(OprogStatus,rb_OprogStatus,sizeof(struct progStatus));
  memcpy(OprogStatus.ptr,rb_OprogStatus.ptr, OprogStatus_rb.len);
  memcpy(dd_coord_ptr,dd_coord_ptr_rb,dd_totBytes);
  dblparnum = sizeof(double)*Oparams.parnum;
  intparnum = sizeof(int)*Oparams.parnum; 
  memcpy(lastcol,lastcol_rb,dblparnum);
  memcpy(atomTime,atomTime_rb,dblparnum);
#ifdef MD_PATCHY_HE
  memcpy(lastbump,lastbump_rb,sizeof(struct LastBumpS)*Oparams.parnum);
#else
  memcpy(lastbump,lastbump_rb,intparnum);
#endif
  memcpy(cellList, cellList_rb, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
  memcpy(inCell[0], inCell_rb[0], intparnum);
  memcpy(inCell[1], inCell_rb[1], intparnum);
  memcpy(inCell[2], inCell_rb[2], intparnum);
#ifdef MD_LL_BONDS
  memcpy(bonds[0], bonds_rb[0], sizeof(long long int)*Oparams.parnum*OprogStatus.maxbonds);
#else
  memcpy(bonds[0], bonds_rb[0], intparnum*OprogStatus.maxbonds);
#endif
  memcpy(numbonds, numbonds_rb, intparnum);
#ifdef MD_SPHERICAL_WALL
  for (i=0; i < Oparams.parnum; i++)
    {
      /* gli ultimi due tipi devono essere i "muri" sferici */
      if (typeOfPart[i]==Oparams.ntypes-1 || typeOfPart[i]==Oparams.ntypes-2)
	{
#ifdef MD_LL_BONDS
	  /* NOTA 21/04/2010: il 3 l'ho messo per tener conto del fatto che nel caso ad esempio 
	  con l'interazione SW i legami possono essere anche due per particella. */
	  memcpy(bonds[i], bonds_rb[i], sizeof(long long int)*Oparams.parnum*MD_MAX_BOND_PER_PART);
#else
	  memcpy(bonds[i], bonds_rb[i], intparnum*MD_MAX_BOND_PER_PART);
#endif
	}
    }
#endif
#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
#ifdef MD_CALENDAR_HYBRID
  memcpy(tree[0], tree_rb[0], 16*poolSize*sizeof(int));
#else
  memcpy(tree[0], tree_rb[0], 13*poolSize*sizeof(int));
#endif
#else
#ifdef MD_CALENDAR_HYBRID
  memcpy(tree[0], tree_rb[0], 13*poolSize*sizeof(int));
#else
  memcpy(tree[0], tree_rb[0], 10*poolSize*sizeof(int));
#endif
#endif
  memcpy(treeTime, treeTime_rb, sizeof(double)*poolSize);
  memcpy(treeRxC, treeRxC_rb,   sizeof(double)*poolSize);
  memcpy(treeRyC, treeRyC_rb,   sizeof(double)*poolSize);
  memcpy(treeRzC, treeRzC_rb,   sizeof(double)*poolSize);
#ifdef MD_ABSORP_POLY
  memcpy(oldTypeOfPart, oldTypeOfPart_rb, intparnum);
#endif
#ifdef MD_GHOST_IGG
  memcpy(ghostInfoArr, ghostInfoArr_rb, sizeof(ghostInfo)*Oparams.parnum);
#endif
  memcpy(typeOfPart, typeOfPart_rb, intparnum);
  memcpy(typeNP, typeNP_rb, sizeof(int)*Oparams.ntypes);
  if (OprogStatus.useNNL)
    {  
      memcpy(nebrTab, nebrTab_rb, sizeof(struct nebrTabStruct)*Oparams.parnum);
      for (i=0; i < Oparams.parnum; i++)
	{
    	  memcpy(nebrTab[i].list, nebrTab_rb[i].list, sizeof(int)*nebrTab[i].len);
	}
    }
#ifdef MD_MATRIX_CONTIGOUS
  memcpy(RM[0], RM_rb[0], dblparnum*9);
#else
  for (i=0; i < Oparams.parnum; i++) 
    {
      memcpy(RM[i], RM_rb[i], 9*sizeof(double));
    }
#endif
#ifdef MD_MATRIX_CONTIGOUS
  memcpy(R[0], R_rb[0], dblparnum*9);
#else
  for (i=0; i < Oparams.parnum; i++) 
    {
      memcpy(R[i], R_rb[i], 9*sizeof(double));
    }
#endif
#ifdef MD_ASYM_ITENS
  memcpy(costheta0, costhet0_rb, dblparnum);
  memcpy(sintheta0, sintheta0_rb, dblparnum);
  memcpy(theta0, theta0_rb, dblparnum);
  memcpy(psi0, psi0_rb, dblparnum);
  memcpy(phi0, phi0_rb, dblparnum);
  memcpy(angM, angM_rb, dblparnum);
#endif
  memcpy(axa, axa_rb, dblparnum);
  memcpy(axb, axb_rb, dblparnum);
  memcpy(axc, axc_rb, dblparnum);
  if (OprogStatus.targetPhi > 0.0)
    {
      memcpy(a0I, a0I_rb, dblparnum);
      memcpy(maxax, maxax_rb, dblparnum);
      memcpy(scdone, scdone_rb, intparnum);
    }
#ifdef MD_CALENDAR_HYBRID
  memcpy(linearLists, linearLists_rb, sizeof(int)*(OprogStatus.nlistsHQ+1));
#endif
  itsF = itsF_rb;
  timeF = timesF_rb;
  itsS = itsS_rb;
  timesS =  timesS_rb;
  numcoll = numcoll_rb;
  itsFNL =  itsFNL_rb;
  timesFNL = timesFNL_rb; 
  timesSNL = timesSNL_rb;
  itsSNL = itsSNL_rb; 
  numcalldist = numcalldist_rb; 
  numdisttryagain = numdisttryagain_rb;
  itsfrprmn = itsfrprmn_rb;
  callsfrprmn = callsfrprmn_rb;
  callsok = callsok_rb;
  callsprojonto = callsprojonto_rb;
  itsprojonto = itsprojonto_rb;
  accngA = accngA_rb;
  accngB = accngB_rb;
}
/* ===== >>>> border zone <<< ===== */
/* max_neigh_regions contiene il numero di neighboring regions, mentre 
   l'array all_neighregions contiene tutte le "max_neigh_regions" neighboring regions */
int max_neigh_regions=0; 
int all_neighregions[26];
void add_to_neigh_regions_arr(int reg)
{
  int r, isnew=1;
 
  for (r = 0; r < max_neigh_regions; r++)
    {
      if (reg == all_neighregions[r])
      {
	isnew=0;
      }
    }
  if (isnew)
    {
      all_neighregions[max_neigh_regions++] = reg;
    }
}
void build_border_zone_list(int regnum)
{
  unsigned int ix, iy, iz, c, inicell, oldc;
  int ixp, iyp, izp, isbordercell;
  unsigned long long int nc;
  int relc[6][3] = {{0,0,1},{0,0,-1},{0,-1,0},{0,1,0},{1,0,0},{-1,0,0}}; 
  /* nc sono le sei neighbour cells lungo gli assi x, y, z */
  /* le border cells sono tutte quelle che hanno vicini in altre regioni */
  inicell = (regnum==0)?0:regionsArr[regnum-1];
  border_cells_head[regnum] = -1;
  for (c = inicell; c < regionsArr[regnum]; c++)
    {
      isbordercell = 0;
      calc_xyz(c, &ix, &iy, &iz);
      for (i = 0; i < 6; i++)
	{
	  ixp = ix + relc[i][0];
	  iyp = iy + relc[i][1];
	  izp = iz + relc[i][2];
	  if (ixp >= cellsz)
	    ixp = 0;
	  else if (ixp <= 0)
	    ixp = cellsx-1;
	  if (iyp >= cellsy)
	    iyp = 0;
	  else if (iyp <= 0)
	    iyp = cellsy-1;
#ifdef MD_EDHEFLEX_WALL
	  if (OprogStatus.hardwall)
	    {
	      /* se c'è il muro non ci sarà un processo adiacente che "scambia" particelle
		 con quello corrente quindi, possiamo evitare di considerare tale cella
		 come appartenente alla border zone */
	      if (izp >= cellsz || izp <= 0)
		izp = iz;
	    }  
	  else
	    {
	      if (izp >= cellsz)
		izp = 0;
	      else if (izp <= 0)
		izp = cellsz-1;
	    }
#else
	  if (izp >= cellsz)
	    izp = 0;
	  else if (izp <= 0)
	    izp = cellsz-1;
#endif
	  nc = calc_cellnum(ixp, iyp, izp); 
	  if (inRegion[nc] != numreg)
	    {
	      isbordercell = 1;
	      add_to_neigh_regions_arr(inRegion[nc]);
	      break;
	    }
	}	
      if (isbordercell == 1)
	{
	  /* add to border cells list */
	  oldc = border_cells_head[regnum];
	  border_cells_head[regnum] = c;
	  border_cells_ll[c] = oldc;
	  cell_type[c] = DD_BORDER_ZONE;
	}
    }
}
extern int cellRange[2 * NDIM];
int in_vborder_list(unsigned long long int cell, unsigned int regnum)
{
  int c;
  c = vborder_cells_head[regnum];
  while (c!=-1)
    {
      if (c == cell)
	return 1;
      c = vborder_cells_head[c];
    }
  return 0;
}
void build_vborder_zone_list(int regnum)
{
  unsigned int ix, iy, iz, c, oldc;
  int ixp, iyp, izp, isvbordercell, dx, dy, dz;
  unsigned long long int nc;
  int cellRangeT[2 * NDIM], k;
  /* nc sono le sei neighbour cells lungo gli assi x, y, z */
  /* le border cells sono tutte quelle che hanno vicini in altre regioni */
  vborder_cells_head[regnum] = -1;

  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  /* scorre tutte le celle della border zone e costruisce
     la lista delle celle della virtual border zone */
  c=vborder_cells_head[regnum];
  while (c!=-1)
    {
      isvbordercell = 0;
      calc_xyz(c, &ix, &iy, &iz);
#ifdef MD_EDHEFLEX_WALL
      if (OprogStatus.hardwall)
	{
	  if (iz + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
	  if (iz + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
	}
#endif
      /* check 26 neighbour cells */
      for (dz = cellRangeT[4]; dz <= cellRangeT[5]; dz++) 
	{
	  izp = iz + dz;  
	  if (izp == -1) 
	    {
	      izp = cellsz - 1;    
	    } 
	  else if (izp == cellsz) 
	    {
	      izp = 0;    
	    }
	  for (dy = cellRangeT[2]; dy <= cellRangeT[3]; dy++) 
	    {
	      iyp = iy + dy;    
	      if (iyp == -1) 
		{
		  iyp = cellsy - 1;    
		} 
	      else if (iyp == cellsy) 
		{
		  iyp = 0;    
		}
	      for (dx = cellRangeT[0]; dx <= cellRangeT[1]; dx++) 
		{
		  ixp = ix + dx;    

		  if (ixp == -1) 
		    {
		      ixp = cellsx - 1;    
		    } 
		  else if (ixp == cellsx) 
		    {
		      ixp = 0;   
		    }
		  nc = calc_cellnum(ixp, iyp, izp); 
		  if (inRegion[nc] != numreg)
		    {
		      if (!in_vborder_list(c, regnum))
			{
			  /* se non è stata già inserita inserisci la cella
			     nella virtual border zone */
			  oldc = vborder_cells_head[regnum];
			  vborder_cells_head[regnum] = c;
			  vborder_cells_ll[c] = oldc;
			  /* VIRTUAL_CELL */
			  cell_type[c] = DD_VIRTUAL;
			}
		    }
		}	
	    }
	}
      c = vborder_cells_head_ll[c];
    }
}
int is_border_zone_cell(unsigned int cn)
{
  if (cell_type[cn]==DD_BORDER_ZONE)
    return 1;
  else
    return 0;
}
void schedule_border_zone_event(int idA, int idB, double tEvent, unsigned int dest_cell)
{
  /* according to Luding and Miller we associate an event to each particle */
  unsigned int cn, ix, iy, iz;
  int idd;
 
  if (idB < ATOM_LIMIT) /* urto fra due particelle */
    {
      ix = inCell[idd][0];
      iy = inCell[idd][1];
      iz = inCell[idd][2];
      cn = calc_cellnum(ix, iy, iz);
      idd = cn+1;

      if (is_border_zone_cell(cn))
	{
	  if (tEvent < treeTimeBZ[idd])
	    ScheduleEventBZ(cn, tEvent); 
	}
      else
       	{
    	  ix = inCell[idd][0];
	  iy = inCell[idd][1];
	  iz = inCell[idd][2];
	  cn = calc_cellnum(ix, iy, iz);
	  idd = cn+1;
	  if (is_border_zone_cell(cn))
	    {
	      if (tEvent < treeTimeBZ[idd])
		ScheduleEventBZ(cn, tEvent);
	    }
	}
    }
  else if (idB < ATOM_LIMIT + 2*NDIM)
    {
      /* calc destination cell cn here*/
      idd = dest_cell+1;
      if (is_border_zone_cell(dest_cell))
	{
	  if (tEvent < treeTimeBZ[idd])
	    ScheduleEventBZ(dest_cell, tEvent);
	}
    }
}
void send_celltime_request_to_region(int regnum, int cellnum, double tEvent)
{
#ifdef MPI
  /* non-blocking send for request */
  MPI_Isend(MPI_COMM_WORLD);
#endif
}
unsigned int num_stored_replies;
struct struct_strep {
  double t_replied; 
  double t_received;
} stored_replies[26];
void store_adj_min_time(int p, double trep, double trec)
{
  stored_replies[p].t_replied=trep;
  stored_replies[p].t_received=trec;
}
double request_cell_time(unsigned int cellnum, double tEvent)
{
  unsigned int ix, iy, iz, c, oldc;
  int ixp, iyp, izp, isvbordercell, dx, dy, dz;
  unsigned long long int nc;
  int cellRangeT[2 * NDIM], k, newregion, completed=0;
  unsigned int neighregions[26], cr, r, numregions=0, cellnum_star;
  double min_time, cellnum_star_time, adj_min_time;

  min_time = tEvent;
  /* request lesser cell time for neighboring cells of cellnum in region regnum */
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
  calc_xyz(cellnum, &ix, &iy, &iz);
#ifdef MD_EDHEFLEX_WALL
  if (OprogStatus.hardwall)
    {
      if (iz + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (iz + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
    }
#endif
  /* check 26 neighbour cells */
  for (dz = cellRangeT[4]; dz <= cellRangeT[5]; dz++) 
    {
      izp = iz + dz;  
      if (izp == -1) 
	{
	  izp = cellsz - 1;    
	} 
      else if (izp == cellsz) 
	{
	  izp = 0;    
	}
      for (dy = cellRangeT[2]; dy <= cellRangeT[3]; dy++) 
	{
	  iyp = iy + dy;    
	  if (iyp == -1) 
	    {
	      iyp = cellsy - 1;    
	    } 
	  else if (iyp == cellsy) 
	    {
	      iyp = 0;    
	    }
	  for (dx = cellRangeT[0]; dx <= cellRangeT[1]; dx++) 
	    {
	      ixp = ix + dx;    

	      if (ixp == -1) 
		{
		  ixp = cellsx - 1;    
		} 
	      else if (ixp == cellsx) 
		{
		  ixp = 0;   
		}
	      nc = calc_cellnum(ixp, iyp, izp); 
	      /* if not current process */
	      if ((cr=inRegion[nc])!=my_rank)
		{
		  /* build a list of all neighboring processes */
		  newregion = 1;
		  for (r=0; r < numregions; r++)
		    {
		      if (cr == neighregions[r])
			{
			  newregion=0;
			  break;
			}
		    }
		  if (newregion)
		    {
		      neighregions[numregions++] = cr;
		    }
		}
	    }
	}
    }
  /* send requests to all neighboring processes */
  for (r = 0; r < numregions; r++)
    send_celltime_request_to_region(neighregions[r], cellnum, tEvent);
  /* each process send a non-blocking broadcast message to all neighboring processes to tell them 
     that there no more mesages to process */
  for (r=0; r < max_neigh_regions; r++)
    send_celltime_request_to_region(all_neighregions[r], -1, 0.0); /* -1 means: "no more messages" from me */
  completed = 0;
  num_stored_replies = 0;
  do
    {
#ifdef MPI
      /* receive all pending messages (requests of cell times)*/
      MPI_Receive(MPI_COMM_WORLD);
#endif
      if (cellnum_star != -1)
	{
	  /* memorizza tutti i reply alla query del tempo minimo nelle celle adiacenti
	     cellnum (notare che memorizza anche il tempo minimo inviato dal processo richiedente
	     cellnum_star_min_time poiche' questo tempo servirà poi per stabilire se e' necessario
	     un rollback (vedi articolo Comp. Phys. Comm.,I. Marin (1997)) */
	  adj_min_time = check_adjacent_cells(cellnum_star);
	  store_adj_min_time(num_stored_replies++, adj_min_time, cellnum_star_min_time);
#ifdef MPI
    	  MPI_Isend(MPI_COMM_WORLD); /* send adj_min_time */
#endif
	}
      else
	completed++;
    }
  while (completed==max_neigh_regions);

  for (r = 0; r < numregions; r++)
    {
#ifdef MPI
      MPI_Receive(MPI_COMM_WORLD);
#endif
      if (adj_min_time < min_time)
	min_time = adj_min_time;
    } 
  return min_time;
}
double dd_calc_tstep(void)
{
#if 0
#ifdef MPI
  /* syncronize all processes, now using MPI */
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  NextEventBZ();
  return request_cell_time(evCellBZ, evTimeBZ);
}
void send_vparticle_to_region(int i)
{
  unsigned int cn, r; 
  /* i is the index of particle to send to process regnum */
  
  if (is_border_zone_cell[cn=dd_inCell[i]])
    {
      /* send particle state to neighboring processes */
      r = neigh_regions_of_cell_list_head[cn];
      while (r!=-1)
	{
	  MPI_Isend(MPI_COMM_WORLD); /* send message to neighboring region r */
	  r = neigh_regions_of_cell_list[r];
	}    
    }
}
void receive_vparticle(void)
{



} 
void dd_syncronize(void)
{


}
#endif
