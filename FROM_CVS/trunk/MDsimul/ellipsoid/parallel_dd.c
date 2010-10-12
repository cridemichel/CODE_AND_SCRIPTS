#include<mdsimul.h>
#ifdef ED_PARALL_DD
/* ==== >>> routines to initialize structures <<< ==== */
extern int cellsx, cellsy, cellsz;
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
	}
  rollback_init();
}

/* ==== >>> rollback <<< ==== */
struct params rb_Oparams;
struct progStatus rb_OprogStatus;
int dd_totBytes; 
void *dd_coord_ptr, *dd_coord_ptr_rb;
double *lastcol_rb, *atomTime_rb, *lastbump_rb;
int *inCell_rb[3], **tree_rb, *cellList_rb;
double *treeRxC_rb, *treeRyC_rb, *treeRzC_rb, *treeTime_rb;
int *numbonds_rb, **bonds_rb, *oldTypeOfPart_rb;
double **R_rb, **RM_rb;
struct nebrTabStruct *nebrTab_rb;
ghostInfo *ghostInfoArr_rb;
int *typeOfPart_rb, *typeNP;
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
  memcpy(rb_Oparams,Oparams,sizeof(struct params));
  memcpy(rb_OprogStatus,OprogStatus,sizeof(struct progStatus));
  memcpy(rb_OprogStatus.ptr,OprogStatus.ptr, OprogStatus.len);
  memcpy(dd_coord_ptr_rb,dd_coord_ptr,dd_totBytes);
}

void rollback_load(void)
{



}

#endif
