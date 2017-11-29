#include<mdsimul.h>
#ifdef ED_PARALL_DD
/* ==== >>> routines to initialize structures <<< ==== */
extern int cellsx, cellsy, cellsz;
int *part_to_send_per_region[26], part_to_send_per_region_head, *part_region;
int **neigh_processes_per_cell, neigh_processes_per_cell_head[26];
double dd_tstep, rb_save_time, *treeTimeBZ;
int **treeBZ;
int num_particles; /* number of particles including virtual ones belonging to current region (process) */ 
int causality_error=0, global_causality_error=0;
int *part_global_idx;
#ifdef MPI
void mpi_check_status(MPI_Status *status)
{


}
#endif
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
inline unsigned long long int calc_cellnum(int ix, int iy, int iz)
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
  idx = ix*cellsy*cellsz + iy*cellsz + iz;
  if (dd_inCell[idx] != -1)
    return dd_inCell[idx];
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
  dd_inCell[idx] = res;
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
int *border_cells_head, *border_cells_ll;
int *dd_inCell, *global_part_idx, *local_part_idx, bz_old_cellnum;
int *neigh_regions_of_cell_list_head, *neigh_regions_of_cell_list;
enum {DD_REAL=0, DD_BORDER_ZONE, DD_VIRTUAL};
/* tags for mpi messages */
enum {MPI_END_OF_MSGS=0, DD_PART_STATE, DD_BZ_EVENT};
/* BZ event type to communicate to neighboring processes */
enum {BZ_PART_REMOVE=0, BZ_PART_ADD, BZ_REG_EV};
int num_part_to_send, head_part_to_send;
#ifdef MPI
MPI_Datatype MPI_PART_STATE;
struct pstate {int i; double x[3]; double v[3]; int inCell[3]; double time; int evtype} *part_state;
const int pstate_fields=6;
void create_mpi_pstate_datatype(struct pstate *ps)
{
  int block_lengths[pstate_fields];
  MPI_Aint displacements[pstate_fields];
  MPI_Aint addresses[pstate_fields+1];
  MPI_Datatype typelist[pstate_fields];
  typelist[0] = MPI_INT;
  typelist[1] = MPI_DOUBLE;
  typelist[2] = MPI_DOUBLE;
  typelist[3] = MPI_INT;
  typelist[4] = MPI_DOUBLE;
  typelist[5] = MPI_INT;
  block_lengths[0] = 1;
  block_lengths[1] = 3;
  block_lengths[2] = 3;
  block_lengths[3] = 3;
  block_lengths[4] = 1;
  block_lengths[5] = 1;
  MPI_Address(ps, &adresses[0]);  
  MPI_Address(&(ps->i),&addresses[1]);
  MPI_Address(&(ps->x),&addresses[2]);
  MPI_Address(&(ps->v),&addresses[3]);
  MPI_Address(&(ps->inCell),&addresses[4]);
  MPI_Address(&(ps->time),&addresses[5]);
  MPI_Address(&(ps->evtype),&addresses[6]);
  displacements[0] = addresses[1]-addresses[0];
  displacements[1] = addresses[2]-addresses[0];
  displacements[2] = addresses[3]-addresses[0];
  displacements[3] = addresses[4]-addresses[0]; 
  displacements[4] = addresses[5]-addresses[0];
  displacements[5] = addresses[6]-addresses[0];
  MPI_Create_type_struct(pstate_fields, block_lengths, displacements, typelist, &MPI_PART_STATE);
  MPI_Type_commit(MPI_PART_STATE);	
} 
MPI_Datatype MPI_CELL_TIME;
struct ct_struct {unsigned int cellnum; double mintime;} *cell_time_msgs;
const int ct_struct_fields=2;

void create_cell_time_event_datatype(struct ct_struct *cs)
{
  int block_lengths[ct_struct_fields];
  MPI_Aint displacements[ct_struct_fields];
  MPI_Aint addresses[ct_struct_fields+1];
  MPI_Datatype typelist[ct_struct_fields];
  typelist[0] = MPI_UINT;
  typelist[1] = MPI_DOUBLE;
  block_lengths[0] = 1;
  block_lengths[1] = 1;
  MPI_Address(ps, &adresses[0]);  
  MPI_Address(&(cs->cellnum),&addresses[1]);
  MPI_Address(&(cs->time),&addresses[2]);
  displacements[0] = addresses[1]-addresses[0];
  displacements[1] = addresses[2]-addresses[0];
  MPI_Create_type_struct(ct_struct_fields, block_lengths, displacements, typelist, &MPI_CELL_TIME);
  MPI_Type_commit(MPI_CELL_TIME);	
} 
void create_mpi_derived_datatypes(void)
{
  struct pstate ps;
  struct ct_struct cs;
  create_mpi_pstate_datatype(&ps);
  create_cell_time_event_datatype(&cs);
}
#endif
void dd_init(void)
{
  /* queste dichiarazioni vanno poi rese globali */
  int *inRegion, dd_numreg, cc, iX, iY, iZ, ncp, BZcalsize;
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
#ifdef MPI
  create_mpi_derived_datatypes();
#endif
  rb_since_save_changed = malloc(sizeof(int)*Oparams.parnum);
  rb_since_load_changed = malloc(sizeof(int)*Oparams.parnum);
  /* global_part_idx[i] da l'indice globale della particella il cui indice locale è i */
  global_part_idx = malloc(sizeof(int)*Oparams.parnum);
  /* local_part_idx[i]  da l'indice locale della particella il cui indice globale è i (quindi 
     è la funzione inversa di global_part_idx[...] */
  local_part_idx = malloc(sizeof(int)*Oparams.parnum);
  /* valore della cella all'inizio della parallel phase: tale 
     valore viene usato per capire quali sono le particelle rimaste nella border
     zone o quali quelle da aggiungere o togliere.
   */
  bz_old_cellnum =  malloc(sizeof(int)*Oparams.parnum);

  for (i=0; i < Oparams.parnum; i++)
   {
     rb_since_load_changed[i] = rb_since_save_changed[i] = 0;
     global_part_idx[i] = i;
     local_part_idx[i] = i;
   }
  /* inRegion[c] gives region associated to cell c */
  inRegion = malloc(sizeof(int)*cellsx*cellsy*cellsx);
  /* array contenente la fine di ogni regione nella sequenza associata alle celle
     tramite interleaving dei bit (vedi Sec. 3.1 S. Miller and S. Luding J. Comput. Phys. 193, 
     306-316 (2003) */
  regionsArr = malloc(sizeof(int)*dd_numreg);
  /* dd_inCell contiene l'indice di cella calcolato tramite interleaving (come 
     spiegato da Luding and Miller, J. Comp. Phys. 2003) */
  dd_inCell =   malloc(sizeof(int)*cellsx*cellsy*cellsz);
  for (i=0; i < cellsx*cellsy*cellsz; i++)
    dd_inCell[i] = -1;
  part_to_send = malloc(sizeof(int)*(Oparams.parnum+1));
  /* part_global_idx[i] gives the global index of particle i */
  part_global_idx = malloc(sizeof(int)*Oparams.parnum);
  num_part_to_send= 0;
  head_part_to_send = -1;
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
  cell_time_msgs = malloc(sizeof(ct_struct)*max_neigh_regions);
  /* part_region[i] = regione (processo) a cui appartiene la particella i-esima */
  part_region = malloc(sizeof(int)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    {
      cellnum = calc_cellnum(inCell[0][i], inCell[1][i], inCell[2][i]);
      part_region[i] = inRegion[cellnum];
    }
  for (k=0; k < 26; k++)
    {
      part_to_send_per_region[k] = malloc(sizeof(int)*Oparams.parnum);
      part_to_send_per_region_head[k] = -1;
    }
  part_state = malloc(sizeof(struct pstate)*Oparams.parnum);
  /* BZ event calendar allocations */
  BZcalsize = cellsx*cellsy*cellsz+1;
  treeTimeBZ = malloc(sizeof(double)*BZcalsize);
  treeBZ = AllocMatI(4,BZcalsize);	  
  /* ============================== */
  rollback_init();
  neigh_processes_per_cell_head = malloc(sizeof(int)*cellsx*cellsy*cellsz);
  neigh_processes_per_cell = malloc(sizeof(int*)*cellsx*cellsy*cellsz);   
  for (k=0; k < cellsx*cellsy*cellsz; k++)
    {
      neigh_processes_per_cell[k] = malloc(sizeof(int)*26);
      neigh_processes_per_cell_head[k] = -1;
    }
  inv_all_neighregions = malloc(sizeof(int)*numOfProcs);
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
#if 0
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
#else
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
#ifdef MD_GHOST_IGG
  ghostInfoArr_rb = malloc(sizeof(ghostInfo)*Oparams.parnum);
#endif
  typeOfPart_rb = malloc(sizeof(int)*Oparams.parnum);
  if (OprogStatus.useNNL)
    {  
      nebrTab_rb = malloc(sizeof(struct nebrTabStruct)*Oparams.parnum);
      for (i=0; i < Oparams.parnum; i++)
	{
    	  nebrTab_rb[i].list = malloc(sizeof(int)*OprogStatus.nebrTabFac);
	}
    }
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
  axa_rb = malloc(sizeof(double)*Oparams.parnum);
  axb_rb = malloc(sizeof(double)*Oparams.parnum);
  axc_rb = malloc(sizeof(double)*Oparams.parnum);
  /* these array are for growth simulations and in this implementation
     growth is not allowed in parallel... (?!?) */
}
#endif
void rb_particle_state_changed(int i)
{
  rb_since_save_changed[i] = rb_since_load_changed[i] = 1;
}
#if 0
void rollback_save(void)
{
  int dblparnum, intparnum, i;
#if 1
  parnum = num_particles;
#else
  parnum = Oparams.parnum;
#endif
  dblparnum = sizeof(double)*parnum;
  intparnum = sizeof(int)*parnum;
  memcpy(rb_Oparams,Oparams,sizeof(struct params));
  memcpy(rb_OprogStatus,OprogStatus,sizeof(struct progStatus));
  memcpy(rb_OprogStatus.ptr,OprogStatus.ptr, OprogStatus.len);
  memcpy(dd_coord_ptr_rb,dd_coord_ptr,dd_totBytes);
  memcpy(lastcol_rb,lastcol,dblparnum);
  memcpy(atomTime_rb,atomTime,dblparnum);
#ifdef MD_PATCHY_HE
  memcpy(lastbump_rb,lastbump,sizeof(struct LastBumpS)*parnum);
#else
  memcpy(lastbump_rb,lastbump,intparnum);
#endif
  memcpy(cellList_rb, cellList, sizeof(int)*(cellsx*cellsy*cellsz+parnum));
  memcpy(inCell_rb[0], inCell[0], intparnum);
  memcpy(inCell_rb[1], inCell[1], intparnum);
  memcpy(inCell_rb[2], inCell[2], intparnum);
#ifdef MD_LL_BONDS
  memcpy(bonds_rb[0], bonds[0], sizeof(long long int)*parnum*OprogStatus.maxbonds);
#else
  memcpy(bonds_rb[0], bonds[0], intparnum*OprogStatus.maxbonds);
#endif
  memcpy(numbonds_rb, numbonds, intparnum);
#ifdef MD_SPHERICAL_WALL
  for (i=0; i < parnum; i++)
    {
      /* gli ultimi due tipi devono essere i "muri" sferici */
      if (typeOfPart[i]==Oparams.ntypes-1 || typeOfPart[i]==Oparams.ntypes-2)
	{
#ifdef MD_LL_BONDS
	  /* NOTA 21/04/2010: il 3 l'ho messo per tener conto del fatto che nel caso ad esempio 
	  con l'interazione SW i legami possono essere anche due per particella. */
	  memcpy(bonds_rb[i], bonds[i], sizeof(long long int)*parnum*MD_MAX_BOND_PER_PART);
#else
	  memcpy(bonds_rb[i], bonds[i], intparnum*MD_MAX_BOND_PER_PART);
#endif
	}
    }
#endif
#ifdef MD_SPHERICAL_WALL
  poolSize = OprogStatus.eventMult*parnum+2*parnum;
#else
  poolSize = OprogStatus.eventMult*parnum;
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
  memcpy(ghostInfoArr_rb, ghostInfoArr, sizeof(ghostInfo)*parnum);
#endif
  memcpy(typeOfPart_rb, typeOfPart, intparnum);
  memcpy(typeNP_rb, typeNP, sizeof(int)*Oparams.ntypes);
  if (OprogStatus.useNNL)
    {  
      memcpy(nebrTab_rb, nebrTab, sizeof(struct nebrTabStruct)*parnum);
      for (i=0; i < parnum; i++)
	{
    	  memcpy(nebrTab_rb[i].list, nebrTab[i].list, sizeof(int)*nebrTab[i].len);
	}
    }
#ifdef MD_MATRIX_CONTIGOUS
  memcpy(RM_rb[0], RM[0], dblparnum*9);
#else
  for (i=0; i < parnum; i++) 
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
#else
void rollback_save(void)
{
  int i, k, k1, k2;
  rb_save_time = Oparams.time; 
  for (i=0; i < Oparams.parnum; i++)
    {
      /* salva lo stato solo delle particelle
	 aggiornate dall'ultimo salvataggio di stato */
      if (rb_since_save_changed[i])
	{
	  atomTime_rb[i] = atomTime[i];
	  rx_rb[i] = rx[i];
	  ry_rb[i] = ry[i];
	  rz_rb[i] = rz[i];
	  vx_rb[i] = vx[i];
	  vy_rb[i] = vy[i];
	  vz_rb[i] = vz[i];
	  wx_rb[i] = wx[i];
	  wy_rb[i] = wy[i];
	  wz_rb[i] = wz[i];
	  /* i semi assi durante una crescita possono cambiare */
	  axa_rb[i] = axa[i];
	  axb_rb[i] = axb[i];
	  axc_rb[i] = axc[i];
	  for (k=0; k < 3; k++)
	    inCell_rb[k][i] = inCell[k][i];
	  rb_since_save_changed[i] = 0;
	  DR_rb[i] = OprogStatus.DR[i]; 
	  Mx_rb[i] = Mx[i];
	  My_rb[i] = My[i];
	  Mz_rb[i] = Mz[i];
	  numbonds_rb[i] = numbonds[i];
	  for (k=0; k < numbonds[i]; k++)
	    bonds_rb[i][k] = bonds[i][k];
	  typeOfPart_rb[i] = typeOfPart[i];
	  oldTypeOfPart_rb[i] = oldTypeOfPart[i];
	  if (Oparams.ghostsim)
	    memcpy(&(ghostInfoArr_rb[i]), &(ghostInfoArr[i]), sizeof(ghostInfo));
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		R_rb[i][k1][k2] = R[i][k1][k2];  
	      }
	}
    }
  /* salva tutte le neighbor list */
  if (OprogStatus.useNNL)
    {  
      memcpy(nebrTab_rb, nebrTab, sizeof(struct nebrTabStruct)*parnum);
      for (i=0; i < parnum; i++)
	{
    	  memcpy(nebrTab_rb[i].list, nebrTab[i].list, sizeof(int)*nebrTab[i].len);
	}
    }

  /* accumulatori vari */
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
#endif
#if 0
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
#else
int doing_rollback_load=0;
void delete_all_events(int i)
{
  int id, idd;
  id=i+1;
#ifdef MD_ABSORPTION
  if (id==sphWall+1)
    continue; 
  if (id==sphWallOuter+1)
    continue; 
#endif
#if 0
  /* TODO: capire meglio come trattare questi eventi */
#ifdef MD_ABSORPTION
  /* NOTA: questo problema si pone solo se non si usa il buffer sferico */
  /* NOTA: gli id che vanno da 0...Oparams.parnum sono riservati ai cell crossing
     o agli urti con le pareti della scatola, nel caso di urto con la membrana 
     semi-permeabile tale evento va rimosso esplicitamente altrimentri rimane nel calendario 
	 degli eventi */
  if (evIdB==ATOM_LIMIT+50)
    DeleteEvent(idNow);
#endif
#endif
  /* delete cell-crossing event*/
  DeleteEvent(id);
  /* delete collisions from calendar and circular lists */
  for (idd = treeCircAL[id]; idd != id; idd = treeCircAL[idd]) 
    {
      /* il successivo (R) del precedente (L) diviene il successivo 
       * del nodo corrente poiché il nodo corrente è stato eliminato */
      treeCircBR[treeCircBL[idd]] = treeCircBR[idd];
      /* il precedente del successivo diviene il precedente del nodo
       * corrente */
      treeCircBL[treeCircBR[idd]] = treeCircBL[idd];
      DeleteEvent (idd);
    }
  /* treeIdA[0] punta al primo nodo della lista dei nodi liberi 
   * nel pool, quindi qui inserisce la lista di nodi appena liberati fra i nodi 
   * non utilizzati del pool */
  treeCircAR[treeCircAL[id]] = treeIdA[0];
  treeIdA[0] = treeCircAR[id];
  /* tutte le liste circolari vengono svuotate */
  treeCircAL[id] = treeCircAR[id] = id;
  for (idd = treeCircBL[id]; idd != id; idd = treeCircBL[idd]) 
    {
      /* vedere sopra infatti è lo stesso solo per la lista in cui
       * la particella è la prima della coppia (A) */
      treeCircAR[treeCircAL[idd]] = treeCircAR[idd];
      treeCircAL[treeCircAR[idd]] = treeCircAL[idd];
      DeleteEvent (idd);
      treeCircAR[idd] = treeIdA[0];    
      treeIdA[0] = idd;
    }
  treeCircBL[id] = treeCircBR[id] = id;
}
inline int dd_is_virtual(int i)
{
  return !(part_region[i] == my_rank);
}
void rollback_load(void)
{
  int i, k, k1, k2, na;
  Oparams.time =  rb_save_time; 
  for (i=0; i < Oparams.parnum; i++)
    {
      /* salva lo stato solo delle particelle
	 aggiornate dall'ultimo salvataggio di stato */
      if (rb_since_save_changed[i])
	{
	  atomTime[i] = atomTime_rb[i];
	  rx[i] = rx_rb[i];
	  ry[i] = ry_rb[i];
	  rz[i] = rz_rb[i];
	  vx[i] = vx_rb[i];
	  vy[i] = vy_rb[i];
	  vz[i] = vz_rb[i];
	  wx[i] = wx_rb[i];
	  wy[i] = wy_rb[i];
	  wz[i] = wz_rb[i];
	  /* i semi assi durante una crescita possono cambiare */
	  axa[i] = axa_rb[i];
	  axb[i] = axb_rb[i];
	  axc[i] = axc_rb[i];
	  for (k=0; k < 3; k++)
	    inCell[k][i] = inCell_rb[k][i];
	  OprogStatus.DR[i] = DR_rb[i]; 
	  Mx[i] = Mx_rb[i];
	  My[i] = My_rb[i];
	  Mz[i] = Mz_rb[i];
	  numbonds[i] = numbonds_rb[i];
	  for (k=0; k < numbonds[i]; k++)
	    bonds[i][k] = bonds_rb[i][k];
	  typeOfPart[i] = typeOfPart_rb[i];
	  oldTypeOfPart[i] = oldTypeOfPart_rb[i];
	  if (Oparams.ghostsim)
	    memcpy(&(ghostInfoArr[i]), &(ghostInfoArr_rb[i]), sizeof(ghostInfo));
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		R[i][k1][k2] = R_rb[i][k1][k2];  
	      }
	  /* TODO: setta scdone e maxax per la crescita */
	  upd_refsysM(i);
	  angM[i] = sqrt(Sqr(Mx[i])+Sqr(My[i])+Sqr(Mz[i]));
	}
    }

  doing_rollback_load=1;
  /* cancella tutti gli eventi in cui sono coinvolte le particelle aggiornate */
  for (i=0; i < Oparams.parnum; i++)
    {
      if (rb_since_load_changed[i] && !dd_is_virtual(i))
	delete_all_events(i);
    }

  /* predice i nuovi eventi per le particelle aggiornate */
  if (OprogStatus.useNNL)
    {  
      memcpy(nebrTab, nebrTab_rb, sizeof(struct nebrTabStruct)*Oparams.parnum);
      for (i=0; i < Oparams.parnum; i++)
	{
    	  memcpy(nebrTab[i].list, nebrTab_rb[i].list, sizeof(int)*nebrTab[i].len);
	}
    }

  for (i=0; i < Oparams.parnum; i++)
    {
      if (OprogStatus.useNNL)
	PredictEventNNL(n, -1);
      else
	PredictEvent(n, -1);
    }
  doing_rollback_load=0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (rb_since_load_changed[i])
	rb_since_load_changed[i] = 0;
    }
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
#endif
/* ===== >>>> border zone <<< ===== */
/* max_neigh_regions contiene il numero di neighboring regions, mentre 
   l'array all_neighregions contiene tutte le "max_neigh_regions" neighboring regions */
int max_neigh_regions=0; 
int all_neighregions[26];
int *inv_all_neighregions;
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
      /*questa di seguito è la funzione inversa */
      inv_all_neighregions[reg] = max_neigh_regions;
    }
}

void build_neigh_processes_list_of_cell(int ix, int iy, int iz, int numreg)
{
  unsigned int cn, k;
  int cellRangeT[2 * NDIM], doneReg[26];
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
  for (k = 0; k < 26; kèè)
   doneReg[k] = 0; 
#ifdef MD_EDHEFLEX_WALL
  if (OprogStatus.hardwall)
    {
      if (iz + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (iz + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
    }
#endif

  cn = calc_cellnum(ix,iy,iz);
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
	      nr = inRegion[nc]
	      if (!doneReg[inv_all_neighregions[nr]])
		{
		  doneReg[inv_all_neighregions[nr]] = 1;
		  /* all_inv_neighregions[] associa alle regioni vicine un numero tra 0 e max_neighregions */
		  neigh_processes_per_cell[cn][inv_all_neighregions[nr]] = neigh_processes_per_cell_head[cn];
		  neigh_processes_per_cell_head[cn] = inv_all_neighregions[nr];
		}	      
	    }
	}
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
	      build_neigh_processes_list_of_cell(ix, iy, iz, inRegion[nc]);
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
inline int ax_get_nc(int idA, int k, double velk, int cellsk)
{
  int nc;
  if (velk > 0.0)
    {
      nc = inCell[k][idA] + 1;
      if (nc == cellsk) 
	{
	  nc = 0;
	}
    }
  else
    { 
      nc = inCell[k][idA] - 1;
      if (nc == -1) 
	{
	  nc = cellsk - 1;
	}
    }
  return nc;
}
inline int calc_dest_cell(int idA, int idB)
{
  int k;
  k = idB - 100 - ATOM_LIMIT; 
  switch (k)
    {
    case 0: 
      return ax_get_nc(idA, 0, vx[idA], &(rx[idA]), cellsx);
      break;
    case 1: 
      return ax_get_nc(idA, 1, vy[idA], &(ry[idA]), cellsy);
      break;
    case 2:
      return ax_get_nc(idA, 2, vz[idA], &(rz[idA]), cellsz);
      break;
    }
}
void schedule_border_zone_event(int idA, int idB, double tEvent, unsigned int dest_cell)
{
  /* according to Luding and Miller we associate an event to each particle */
  unsigned int cn, ix, iy, iz;
  int idd;
  ix = inCell[0][idA];
  iy = inCell[1][idA];
  iz = inCell[2][idA];
  cn = calc_cellnum(ix, iy, iz);
  idd = cn+1;

  if (is_border_zone_cell(cn))
    {
      if (tEvent < treeTimeBZ[idd])
	ScheduleEventBZ(cn, tEvent); 
    }

  if (idB < ATOM_LIMIT) /* urto fra due particelle */
    {
      ix = inCell[0][idB];
      iy = inCell[1][idB];
      iz = inCell[2][idB];
      cn = calc_cellnum(ix, iy, iz);
      idd = cn+1;

      if (is_border_zone_cell(cn))
	{
	  if (tEvent < treeTimeBZ[idd])
	    ScheduleEventBZ(cn, tEvent); 
	}
    }
  else if (idB < ATOM_LIMIT + 2*NDIM)
    {
      /* TODO: calc destination cell cn here*/
      dest_cell = calc_dest_cell(idA, idB);
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
  MPI_Request request;
  struct ct_struct min_cell_time;
  /* non-blocking send for request */
  MPI_Isend(min_cell_time, 1, MPI_CELL_TIME, regnum, DD_BZ_EVENT, MPI_COMM_WORLD, &request);
#endif
}
unsigned int num_stored_replies;
struct struct_strep {
  unsigned int cellnum;
  double t_replied; 
  double t_received;
} stored_replies[26];

void store_adj_min_time(int p, double trep, double trec, unsigned int cn)
{
  stored_replies[p].cellnume = cn;
  stored_replies[p].t_replied=trep;
  stored_replies[p].t_received=trec;
}

double request_cell_time(unsigned int cellnum, double tEvent)
{
  unsigned int ix, iy, iz, c, oldc;
  int ixp, iyp, izp, isvbordercell, dx, dy, dz, source_process;
  unsigned long long int nc;
  int cellRangeT[2 * NDIM], k, newregion, completed=0;
  unsigned int neighregions[26], cr, r, numregions=0, cellnum_star;
  double min_time, cellnum_star_time, adj_min_time;
#ifdef MPI
  MPI_Status status;
  MPI_Request request;
#endif
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
  num_ct_received=0;
  do
    {
#ifdef MPI
      /* receive all pending messages (requests of cell times)*/
      MPI_Receive(&(cell_time_msgs[num_ct_received++]), MPI_BZ_EVENT, MPI_ANY_SOURCE, DD_BZ_EVENT, 
		  MPI_COMM_WORLD, &status);
      mpi_check_status(status);
      source_process=status.MPI_SOURCE;
#endif
      cellnum_star = cell_time_msgs[num_ct_received].cellnum;
      cellnum_star_min_time = cell_time_msgs[num_ct_received].mintime;
      if (cellnum_star != -1)
	{
	  /* memorizza tutti i reply alla query del tempo minimo nelle celle adiacenti
	     cellnum (notare che memorizza anche il tempo minimo inviato dal processo richiedente
	     cellnum_star_min_time poiche' questo tempo servirà poi per stabilire se e' necessario
	     un rollback (vedi articolo Comp. Phys. Comm.,I. Marin (1997)) */
	  adj_min_time = check_adjacent_cells(cellnum_star);
	  store_adj_min_time(num_stored_replies++, adj_min_time, cellnum_star_min_time, cellnum_star);
#ifdef MPI
	  MPI_Isend(min_cell_time, 1, DD_BZ_EVENT, source_process, 
		    MPI_COMM_WORLD, &request); /* send adj_min_time */
#endif
	}
      else
	completed++;
    }
  while (completed==max_neigh_regions);

  for (r = 0; r < numregions; r++)
    {
#ifdef MPI
      MPI_Receive(MPI_COMM_WORLD,status);
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
void bz_particle_to_send(int i)
{
  unsigned int cn, r; 
  /* i is the index of particle to send to process regnum,
     set all BZ particles to send at the end of parallel phase. */
  /* TODO: meglio sarebbe un flag anziché la lista altrimenti rischio di 
     avere più particelle uguali nella lista! */
  if (is_border_zone_cell[cn=calc_cellnum(inCell[0][i],inCell[1][i],inCell[2][i])])
    {
      /* update linked list of particles to send */
      part_to_send[i] = head_part_to_send; 
      head_part_to_send = i;
      num_part_to_send++;
    }
}

int npart_in_buf=0;
void add_part_to_mpi_buffer(int ipart)
{
  struct pstate *ps;
  unsigned int cn, oldcn;
  if (ipart==-1)
    {
      ps = &(part_state[npart_in_buf++]);
      ps->i = -1;
      return;
    }
  ps = &(part_state[npart_in_buf++]);
  /* send the local absolute identifier of particle ipart */
  ps->i = global_part_idx[ipart];
  ps->x[0] = rx[ipart];
  ps->x[1] = ry[ipart];
  ps->x[2] = rz[ipart];
  ps->v[0] = vx[ipart];
  ps->v[1] = vy[ipart];
  ps->v[2] = vz[ipart];
  ps->inCell[0] = inCell[0][ipart];
  ps->inCell[1] = inCell[1][ipart];
  ps->inCell[2] = inCell[2][ipart];
  ps->time = atomTime[ipart];
  cn = calc_cellnum(inCell[0][ipart], inCell[1][ipart], inCell[2][ipart]); 
  oldcn = bz_old_cellnum[ipart];
  if (is_border_zone_cell[cn] && is_border_zone_cell[oldcn])
    ps->evtype = BZ_REG_EV;
  else if (is_border_zone_cell[cn] && !is_border_zone_cell[oldcn])
    ps->evtype = BZ_PART_ADD;
  /* rimuove la particella se passa da una border zone cell del processo corrente
     ad una cella non border zone, poiché in tal caso non sarà più una particella
     virtuale per il neighboring process */
  else if (inRegion[cn] == my_rank && !is_border_zone_cell[cn] && is_border_zone_cell[oldcn])
    ps->evtype = BZ_PART_REMOVE; 
  else /* in quest'ultimo caso la particella è diventata reale in una neighboring region */
    ps->evtype = BZ_REG_EV;
}

void send_bz_particles(void)
{
  int i, k;
  unsigned int cn, p;
#ifdef MPI
  MPI_Request request;
#endif
  /* send all border zones particles updated during parallel phase 
     (use a unique mpi call for optimizing this stage) */
  i = part_to_send_head;
  while (i!=-1)
    {
      /* part_to_send_per_region_head[k] indica la prima particella della linked list da mandara
	 alla regione k */
      cn = calc_cellnum(inCell[0][i], inCell[1][i], inCell[2][i]);
      /* notare che p è un numero compreso tra 0 e max_neighregions */
      p = neigh_processes_per_cell_head[cn];
      while (p != -1)
	{
	  part_to_send_per_region[p][i] = part_to_send_per_region_head[p];
	  part_to_send_per_region_head[p] = i;
	  p = neigh_processes_per_cell[cn][p];
	}
      i = part_to_send[i];
    }
  for (p = 0; p < max_neighregions; p++)
    {
      npart_in_buf = 0;
      i = part_to_send_per_region_head[p];
      while (i!=-1)
	{
	  add_part_to_mpi_buffer(i);
	  i = part_to_send_per_region[p][i];
	}
      add_part_to_mpi_buffer(-1); /* to mark end of particles */
#ifdef MPI
      /* send particles to region (processor) all_neighregions[p] */
      MPI_Isend(part_state, npart_in_buf, MPI_PART_STATE, 
		all_neighregions[p], DD_PART_STATE, MPI_COMM_WORLD, &request);
#endif
    }
#if 0
  for (p = 0; p < 26; p++)
    {
      /* send a message to all neighbors processes with mpi tag=MPI_END_OF_MSGS
	 meaning "end of messages" */
      MPI_Isend(&eom, 1, MPI_CHAR, all_neighregions[p], MPI_END_OF_MSGS, MPI_COMM_WORLD, &request);
    }
#endif
}
double do_rollback(void)
{
  rollback_load();
}
void schedule_syncronization(double t)
{
  ScheduleEvent(-1, ATOM_LIMIT+20, t);
}

void process_causality_error(double t)
{
  do_rollback();
  schedule_syncronization(t);
}
void dd_updateCalendar(int i)
{
  int ipart;
  ipart = local_part_idx[part_state[i].i];
  if (OprogStatus.useNNL)
    {
      /* ricalcola i tempi di collisione con la NL */
      updrebuildNNL(ipart);
      PredictEventNNL(ipart, -1);
    }
  else
    {
      PredictEvent(ipart, -1);
    }
}

inline int check_causality_error_for_particle(int i)
{
  int ipart;
  /* verificare se basta questa condizione e se è giusta */
  ipart = local_part_idx[part_state[i].i];
  if (atomTime[ipart] < Oparams.time)
    return 1;
  else
    return 0;
}

void bz_remove_particle(int i)
{
  /* mark particle to remove during particles sorting (-1 means "remove") */
  local_part_idx[global_part_idx[i]] = -1;
  global_part_idx[i] = -1; 
}
inline int bz_add_particle(int glob_idx)
{
  Oparams.parnum++;
  local_part_idx[glob_idx] = Oparams.parnum-1;
  global_part_idx[Oparams.parnum-1] = glob_idx;
}
inline void dd_updateParticleState(int i)
{ 
  int ipart, k, evtype;
  unsigned int cellnum; 
  if (part_state[i].evtype == BZ_PART_REMOVE)
    return;
  ipart = local_part_idx[part_state[i].i];
  cellnum = calc_cellnum(part_state[i].inCell[0], part_state[i].inCell[1], part_state[i].inCell[2]);
  rx[ipart] = part_state[i].x[0];
  ry[ipart] = part_state[i].x[1];
  rz[ipart] = part_state[i].x[2];
  vx[ipart] = part_state[i].v[0];
  vy[ipart] = part_state[i].v[1];
  vz[ipart] = part_state[i].v[2];
  for (k = 0; k < 3; k++)
    inCell[k][ipart] = part_state[i].inCell[k];
  atomTime[ipart] = part_state[i].time;
}

void dd_update_particles_state(void)
{
  int i, j, ipart;

  for (i=0; part_state[i].i != -1; i++)
    {
      /* mark particles to remove */
      if (part_state[i].evtype == BZ_PART_REMOVE)
	{
	  ipart = local_part_idx[part_state[i].i];
	  /* elimina dal calendario tutti gli eventi in cui è coinvolta
	     la particella rimossa */
	  delete_all_events(ipart);
	  bz_remove_particle(ipart);
	  part_to_remove++;
	}
    }
  j=0;
  /* inserisce le nuove particelle negli "slot" lasciati liberi dalle particelle rimosse,
     se non sono abbastanza accoda le particelle nuove a quelle già presenti nella regione attuale */
  for (i=0; part_state[i].i != -1 && j < Oparams.parnum; i++)
    {
      if (part_state[i].evtype == BZ_PART_ADD)
	{
	  while (j < Oparams.parnum && global_part_idx[j] != -1)
	    j++; 
	  if (j < Oparams.parnum)
	    {
	      global_part_idx[j] = part_state[i].i;
	      local_part_idx[part_state[i].i] = j;
	    }
	}
    }
  /* add remaining particles received */
  for (; part_state[i].i != -1; i++)
    {
      if (part_state[i].evtype == BZ_PART_ADD)
	{
	  bz_add_particle(part_state[i].i);
	}
    }	
  /* aggiorna lo stato delle particelle nuove e di quelle da aggiornare soltanto */
  for (i=0; part_state[i].i != -1; i++)
    {
      dd_updateParticleState(i);
    }

  /* Predice i nuovi eventi relativi alle particelle ricevute (nuove ed esistenti) */
  for (i=0; part_state[i].i != -1; i++)
    {
      if (part_state[i].evtype == BZ_PART_REMOVE)
	continue;
      dd_updateCalendar(i);
      if (check_causality_error_for_particle(i))
	{
	  causality_error=1;
	}
    } 
}
void receive_vparticles(void)
{
  int completed = 0;
#ifdef MPI
  MPI_Status status;
#endif
  do
    {
#ifdef MPI
      MPI_Receive(part_state, Oparams.parnum, MPI_PART_STATE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
		  MPI_COMM_WORLD, &status);
      mpi_check_status(status);
#endif
      dd_update_particles_state();
      completed++;
    }
  while (completed < 26); /* se la particella è -1 vuol dire che i messaggi sono finiti */
  
} 
void check_causality_error_messages(void)
{ 
  int p;
#if 0
  for (p=0; p < numOfProcess; p++)
    {
      if (causality_error_arr[p])
	causality_error=1;
    }
#endif
  if (global_causality_error)
    {
      dd_syncronize();
      process_causality_error();
    }
}
void check_causality_error_bzevents(void)
{
  int kk, error_detected=0;
  double adj_min_time;
  /* TODO: qui deve controllare che i tempi di cella mandati in risposta
     ad altri processi siano corretti */
  for (k=0; k < num_stored_replies; k++)
    {
      adj_min_time = check_adjacent_cells(stored_replies[k].cn);
      if (stored_replied[k].trep >= stored_replies[k].trec && adj_min_time < stored_replies[k].trec) 
	{
	  error_detected = 1;
	  break;
	}
    }
  if (error_detected)
    causality_error=1;
  else
    causality_error=0;
}
int check_BZ_event(int i, int evIdB)
{
  int dest_cell;
  if (i < ATOM_LIMIT)
    {
      /* following conditions means that a BZ event has not been anticipated
	 hence stop parallel phase immediately in order to avoid a causality error */
      if (evIdB >= ATOM_LIMIT && evIdB < ATOM_LIMIT + 2*NDIM)
	{
	  /* TODO: calculate destination cell here */
	  dest_cell = calc_dest_cell(evIdA, evIdB);
	  if (is_border_zone_cell(dest_cell) && Oparams.time < dd_step)
	    return 1;
	}
      if (is_border_zone_cell[calc_cellnum(inCell[0][i],inCell[1][i],inCell[2][i])] && Oparams.time < dd_tstep)
	return 1;
    }
  return 0;
}
int check_tstep(double t)
{
  causality_error = 0;
  if (t >= dd_tstep)
    {
      /* ==== >>> END OF PARALLEL PHASE <<< ==== */
      send_bz_particles();
      receive_vparticles();
      /* la funzione check_causality_error_bzevents() setta la variabile globale causality_error */
      check_causality_error_bzevents();
#ifdef MPI
      //MPI_Allgather(MPI_COMM_WORLD, causality_error_arr);
      /* check if at least one regione raised a causality error */
      MPI_Allreduce(&causality_error, &global_causality_error, 1, MPI_INTEGER, MPI_LOR, MPI_COMM_WORLD);
#endif
      check_causality_error_messages();
      /* initialize list of BZ particles to send */
      head_part_to_send = -1;
      /* store current cell numbers associated to particles */
      for (i=0; i < Oparams.parnum; i++)
	{
	  bz_old_cellnum[i] = calc_cellnum(inCell[0][i], inCell[1][i], inCell[2][i]);
	}
      return 1;
    }
  return 0;
}
void dd_syncronize(void)
{
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  rollback_save();
}
/* ================ >>> BZ event calendar <<< ================= */
double timeBZ;
int evIdABZ;
void InitCalendarBZ(void)
{
  int id, i;
  treeLeftBZ[0] = treeRightBZ[0] = -1;
  //treeIdABZ[0] = 1;
  /* i nodi da Oparams.parnum + 1 (compreso) in poi sono il pool (cioè quelli dinamici) */
}
void DeleteEventBZ(int idd)
{
  int idp, idq, idr;
  idr = treeRightBZ[id];
  if (idr == -1)
    idq = treeLeftBZ[id];
  else 
    {
      if (treeLeftBZ[id] == -1) 
	idq = idr;
      else 
	{
	  if (treeLeftBZ[idr] == -1) 
	    idq = idr;
	  else 
	    {
	      idq = treeLeftBZ[idr];
	      while (treeLeftBZ[idq] > -1)
		{
		  idr = idq;    
		  idq = treeLeftBZ[idr];
		}
	      treeLeftBZ[idr] = treeRightBZ[idq];
	      if (treeRightBZ[idq]>-1)
	        treeUpBZ[treeRightBZ[idq]] = idr;
	      treeRightBZ[idq] = treeRightBZ[id];
	      if (treeRightBZ[id]>-1)
	        treeUpBZ[treeRightBZ[id]] = idq;
	    }
	  treeUpBZ[treeLeftBZ[id]] = idq;
	  treeLeftBZ[idq] = treeLeftBZ[id];
	} 
    }
  idp = treeUpBZ[id];    
  if (idq > -1)
    treeUpBZ[idq] = idp;
  if (treeRightBZ[idp] != id)
    treeLeftBZ[idp] = idq;
  else 
    treeRightBZ[idp] = idq;
}
void ScheduleEventBZ(int idA, double tEvent)
{
  int id, idNew, more;
  id = 0;
  /* N.B. il numero di eventi in tale calendario di BZ è pari 
     al numero di celle, cioè si ha un evento per cella per cui 
     non è necessaria alcuna gestione dinamica (ad es. liste circolari) */
  idNew = idA + 1;
  MD_DEBUG34(printf("idNew=%d\n", idNew));

  if (treeRightBZ[id] == -1) 
    treeRightBZ[id] = idNew;
  else 
    {
      /* Cerca la giusta collocazione nell'albero per l'evento da
       * schedulare */
      more = 1; 
      id = treeRightBZ[id];
      while (more) 
	{
	  if (tEvent <= treeTimeBZ[id]) 
	    {
	      if (treeLeftBZ[id] > -1) 
		id = treeLeftBZ[id];
	      else 
		{
		  more = 0;    
		  treeLeftBZ[id] = idNew;
		}
	    } 
	  else
	    {
	      if (treeRightBZ[id] > -1) 
		id = treeRightBZ[id];
	      else 
		{
		  more = 0;    
		  treeRightBZ[id] = idNew;
		} 
	    }
	} 
    }
    
  treeTimeBZ[idNew] = tEvent;
  treeIdABZ[idNew] = idA;    
  treeLeftBZ[idNew] = treeRightBZ[idNew] = -1;
  treeUpBZ[idNew] = id;
}
void NextEventBZ(void)
{
  int id, idAx, idBx, idd, idNow, idtx;
  /* Il nodo root (0), a cui non è associato alcun evento,
   * è linkato con il suo right pointer al primo nodo che contiene
   * un evento */
  idNow = treeRightBZ[0];  
  /* Cerca l'evento con tempo minore 
   * NOTA: l'albero è ordinato e ogni nodo sinistro ha un tempo inferiore */
  while (treeLeftBZ[idNow] > -1) 
    idNow = treeLeftBZ[idNow];
  timeBZ = treeTimeBZ[idNow];   
  evIdABZ = treeIdABZ[idNow];    
     
  /* L'evento è un cell-crossing o un evento generico (output o misura) */
  DeleteEvent (idNow);
  treeIdABZ[0] = idNow;
}
/* ============================================================ */
#endif
