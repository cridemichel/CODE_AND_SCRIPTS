#include<mdsimul.h>
#ifdef EDHE_FLEX
/* questa serve anche se non si usano le LL multiple
   per la funzione assign_bond_mapping ottimizzata */
extern char TXT[MSG_LEN];
#ifdef MD_GHOST_IGG
extern int areGhost(int i, int j);
#endif

#ifdef MD_LL_BONDS
extern long long int *bondscache, **bonds;
extern int *numbonds;
#else
extern int *bondscache, *numbonds, **bonds;
#endif
extern int *mapbondsa;
extern int *mapbondsb;

extern double *dists;
extern int **tree;
extern void assign_bond_mapping(int i, int j);
extern int *mapbondsaFlex, *mapbondsbFlex, nbondsFlex;
extern int bound(int na, int n, int a, int b);
extern double calcDistNegSP(double t, double t1, int i, int j, double shift[3], int *amin, int *bmin, double *dists, int bondpair);
extern void check_shift(int i, int j, double *shift);
int get_linked_list_type(int typena, int nc)
{
  int sum, t1, nc1;
  /* per N tipi le linked lists sono N*(N+1)/2. 
     typeOfPart[na]-nc -> nl
     consideriamo il caso ntypes=4
     0-0 -> 0
     1-1 -> 1
     2-2 -> 2
     3-3 -> 3
     0-0 -> X
     0-1 -> 4 = 4 + 1 - 1
     0-2 -> 5 = 4 + 2 - 1
     0-3 -> 6 = 4 + 3 - 1
     1-0 -> X
     1-1 -> X
     1-2 -> 7 = 4 + 4*1 + 2 - 3 
     1-3 -> 8 = 4 + 4*1 + 3 - 3  
     2-0 -> X
     2-1 -> X
     2-2 -> X
     2-3 -> 9 = 4 + 4*2 + 3 - 6
     in quest'ultimo caso il numero da sottrarre è pari alla somma delle X ossia
     6 = 1 + 2 + 3 = 3*(3+1)/2 = (typena+1)*(typena+2)/2
     inoltre 2=typne 4=ntypes e 3=nc e con questo si ottiene la formula riportata sotto
   */
  //typena = typeOfPart[na];
  if (nc==typena)
    return nc;
  else 
    {
      if (typena < nc)
	{
	  t1 = typena+1;
	  sum = t1*(t1+1)/2;
	  return Oparams.ntypes*t1 + (nc-sum);

	}
      else
	{
	  /* lo scambio serve poiché ad es. la lista 3-2 deve essere identica alla lista 2-3 */
	  nc1 = nc+1;
	  sum = nc1*(nc1+1)/2;
	  return Oparams.ntypes*nc1 + (typena-sum);
	}
    }
      /*
      typena=3 nc=2 => 3 + 3*(4 - (2+2)/2) / 2 = 9 OK  
      typena=1 nc=2 => 4*2 + 2 - 3 = 7 OK
   */
}
#endif
extern struct nebrTabStruct *nebrTab;
#ifdef MD_MULTIPLE_LL
#define MD_NNLPLANES
#ifdef MD_ABSORPTION
extern int *listtmp;
#ifdef MD_SPHERICAL_WALL
extern int sphWall, sphWallOuter;
extern void locate_spherical_wall(int na, int outer);
#endif
#endif
/* NOTA: nel caso della silica si avevano 4 linked lists in verità ossia:
   0 è la lista della specie A per l'interazione A-A
   1 è la lista della specie B per l'interazione B-B
   2 è la lista costituita da molecole B per l'interazione A-B
   3 è la lista costituita da molecole A per l'interazione B-A
   ossia per due tipi ho 2*2 linked lists.
   Qui invece vorrei implementare le linked in modo che 
   se ntypes è il numero di tipi il numero di linked lists  
   sia (ntypes+1)*ntypes/2 (nel caso della silica ad es. ne avremmo 3).
   Per ottenere questo le liste miste conterranno sia particelle di un tipo 
   che di un altro e quando si dovranno predire gli eventi si dovrà considerare
   che solo l'interazione mista va considerata scartando quella diretta.
*/
#ifdef MD_EDHEFLEX_OPTNNL
extern int *inCell_NNL[3], *cellList_NNL;
extern double *rxNNL, *ryNNL, *rzNNL;
#endif
extern double nextNNLrebuild;
extern int cellRange[2*NDIM];
int **crossevtodel;
extern const double timbig;
#ifdef MD_LXYZ
extern double L2[3];
#else
extern double L2;
#endif 
#ifdef MD_PATCHY_HE
extern int evIdC, evIdD, evIdE;
extern double *treeRxC, *treeRyC, *treeRzC;
extern double rxC, ryC, rzC;
#endif
double *rcutMLL;
int ***inCellMLL;
int **cellListMLL;
int *cellsxMLL, *cellsyMLL, *cellszMLL, *ignoreMLL;
/* con questo switch si puo' disabilitare l'ottimizzazione delle LL multiple
   con la quale le LL tra tipi non interagenti non vengono utilizzate. */

#if defined(EDHE_FLEX) && defined(MD_OPT_MULTLL)
extern int is_in_ranges(int A, int B, int nr, rangeStruct* r);
int may_interact_core_type(int typei, int typej)
{
  if (!typesArr[typei].ignoreCore && 
      !typesArr[typej].ignoreCore)
   return 1;
 else
   return 0; 
}
int may_interact_spots_type(int type1, int type2)
{
  int ni;
#if 0
  rangeStruct *r1, *r2;
  int nr1, nr2, inti, intj;
#endif
  if (typesArr[type1].nspots == 0 || typesArr[type2].nspots == 0)
    return 0;
  for (ni = 0; ni < Oparams.ninters; ni++)
    {
      if (is_in_ranges(type1, intersArr[ni].type1, intersArr[ni].nr1, intersArr[ni].r1) && 
	  is_in_ranges(type2, intersArr[ni].type2, intersArr[ni].nr2, intersArr[ni].r2))
	{
	  return 1;
	}	
      else if (is_in_ranges(type2, intersArr[ni].type1, intersArr[ni].nr1, intersArr[ni].r1) && 
	       is_in_ranges(type1, intersArr[ni].type2, intersArr[ni].nr2, intersArr[ni].r2))
	{
	  return 1;
	}
    }
  /* nintersIJ sono le interazioni specifiche tra particelle, se ci sono 
     non possiamo a priori scartare tutte le interazioni tra due tipi */
  if (Oparams.nintersIJ > 0)
    {
      return 1;
    }
  return 0;
}
int may_interact_all_type(int t1, int t2)
{
  if (t1==t2)
    return 1;
  if (may_interact_core_type(t1, t2))
    return 1;
  if (may_interact_spots_type(t1, t2))
    return 1;
	
  return 0;
}
#endif

int is_superellips_type(int pt)
{
  if (typesArr[pt].n[0]==2.0 && typesArr[pt].n[1]==2.0 &&
      typesArr[pt].n[2]==2.0)
    return 0;
  else 
    return 1;
}

int get_linked_list(int na, int nc)
{
  return get_linked_list_type(typeOfPart[na], nc);

}
extern int all_spots_in_CoM(int pt);
#ifdef MD_SUPERELLIPSOID
extern int is_superellipse(int i);
#endif
int is_a_sphere_NNL_type(int pt)
{
  if (!(typesArr[pt].sax[0] == typesArr[pt].sax[1] && 
	typesArr[pt].sax[1] == typesArr[pt].sax[2] )) 
    {
      return 0;
    }
#ifdef MD_SUPERELLIPSOID
  if (is_superellipse(pt))
    {
      return 0;
    }
#endif
  if (typesArr[pt].nspots > 0 && !all_spots_in_CoM(pt))
    {
      return 0;
    }
  return 1;
}
#if 0
void get_types_from_nl(int nl, int *t1, int *t2)
{
  int ta, tb, numll, l;
  /* funzione inversa di get_linked_list */
  static int *mat[2]={NULL,NULL};

  if (mat[0]==NULL)
    {
      numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
      mat[0] = malloc(sizeof(int)*numll);
      mat[1] = malloc(sizeof(int)*numll);
      for (l=0; l < numll; l++)
	{
	  for (ta = 0; ta < Oparams.ntypes; ta++)
	    for (tb = 0; tb < Oparams.ntypes; tb++)
	      {
		if (ta > tb)
		  continue;

		if (l==get_linked_list_type(ta, tb))
		  {
		    mat[0][l] = ta;
		    mat[1][l] = tb;
		    //*t1 = ta;
		    //*t2 = tb;
		    break;
		  }
	      }
	}
    }
  *t1 = mat[0][nl];
  *t2 = mat[1][nl];
}
#else
void get_types_from_nl(int nl, int *t1, int *t2)
{
  int ta, tb;
  /* funzione inversa di get_linked_list */
  for (ta = 0; ta < Oparams.ntypes; ta++)
    for (tb = 0; tb < Oparams.ntypes; tb++)
      {
	if (ta > tb)
	  continue;

	if (nl==get_linked_list_type(ta, tb))
	  {
	    *t1 = ta;
	    *t2 = tb;
	    return;
	  }
      }
}
#endif
extern double max3(double a, double b, double c);

double calc_rcut_type(int t)
{
  double rcutA;
  int kk;
  double ax[3], del;
  for (kk=0; kk < 3; kk++)
    ax[kk] = typesArr[t].ppsax[kk];

  if (OprogStatus.useNNL)  
    del = OprogStatus.rNebrShell;
  else
    del = 0.0;
  /* nel caso si tratti di un oggetto a simmetria sferica l'orientazione rimane della NNL (che è un cubo) rimane invariata
     nel tempo per cui si può prendere rcut appena più grande del lato della NNL cubica*/
  if (is_a_sphere_NNL_type(t))
    rcutA = 2.0*max3(ax[0]+del,ax[1]+del,ax[2]+del);
  else
    rcutA = 2.0*sqrt(Sqr(ax[0]+del)+Sqr(ax[1]+del)+Sqr(ax[2]+del));
  //printf("rcutFact[%d]=%f\n", t, typesArr[t].rcutFact);
  if (typesArr[t].rcutFact > 1.0)
    return typesArr[t].rcutFact*rcutA;
  else
    return OprogStatus.rcutfactMLL*rcutA;
}

double calc_rcut(int nl)
{
  int t1=-1, t2=-1;
  double rc, rc1, rc2;
  /* le celle liste per ora vengono sempre scelte automaticamente */
  get_types_from_nl(nl, &t1, &t2);
  rc1 = calc_rcut_type(t1);
  rc2 = calc_rcut_type(t2);
  rc = (rc1 + rc2)*0.5;
  return rc;
}
#if 0
void check_equal(int nl)
{
  int k;
  for (k=0; k < nl; k++)
    {
      if (cellsxMLL[nl]==cellsxMLL[k] &&
	  cellsyMLL[nl]==cellsyMLL[k] &&
	  cellszMLL[nl]==cellszMLL[k] && MLLremapto[k]==k)
	{
	  /* rimappa la linked list attual nl ad una esistente uguale
	     rimappando a linked list che non sono state già rimappate (condizione
	     MLLremapto[k]=k) */
	  MLLremapto[nl] = k;
	}
    } 
}
#endif
void set_cells_size(void)
{
  int nl, numll;
  double rcut;
  
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  for (nl = 0; nl < numll; nl++)
    {
      //if (Oparams.rcut[nl] <= 0.0)
	//Oparams.rcut[nl] = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
      rcut = rcutMLL[nl] = calc_rcut(nl);
#ifdef MD_LXYZ
      cellsxMLL[nl] = L[0] / rcut;
      cellsyMLL[nl] = L[1] / rcut;
      cellszMLL[nl] = L[2] / rcut;
#else
      cellsxMLL[nl] = L / rcut;
      cellsyMLL[nl] = L / rcut;
      cellszMLL[nl] = L / rcut;
#endif
   }

}
extern int evIdA, evIdB;
extern int cellRange[2*NDIM];
int check_boxwall(int k, int nc, int nl)
{
  int cellsk=0;
  double vel=0.0;
  switch (k)
    {
    case 0:
      cellsk = cellsxMLL[nl];
      vel = vx[evIdA];
      break;
    case 1:
      cellsk = cellsyMLL[nl];
      vel = vy[evIdA];
      break;
    case 2:
      cellsk = cellszMLL[nl];
      vel = vz[evIdA];
      break;
    }
  //printf("CHECK BOXWALL k=%d inCell[%d][%d][%d]:%d\n", k, nc, k, evIdA, inCell[nc][k][evIdA]);
  if ((vel < 0 && inCellMLL[nc][k][evIdA]==0) || (vel > 0 && inCellMLL[nc][k][evIdA]==cellsk-1))
    return 1;
  else 
    return 0;
}
extern void DeleteEvent(int id);
void docellcross2MLL(int k, double velk, int cellsk, int nc)
{
  if (velk > 0.0)
    {
      inCellMLL[nc][k][evIdA] = inCellMLL[nc][k][evIdA] + 1;
      cellRange[2 * k] = 1;
      if (inCellMLL[nc][k][evIdA] == cellsk) 
	inCellMLL[nc][k][evIdA] = 0;
    }
  else
    { 
      cellRange[2 * k + 1] = -1;
      inCellMLL[nc][k][evIdA] = inCellMLL[nc][k][evIdA] - 1;
      if (inCellMLL[nc][k][evIdA] == -1) 
	inCellMLL[nc][k][evIdA] = cellsk - 1;
    }
}
void docellcrossMLL(int k, double velk, double *rkptr, int cellsk, int nc)
{
  if (velk > 0.0)
    {
      inCellMLL[nc][k][evIdA] = inCellMLL[nc][k][evIdA] + 1;
#if 0
      if (evIdA==511)
	printf("QUI nc=%d k=%d inCellMLL=%d %d %d cellsk=%d\n", nc, k, inCellMLL[nc][0][evIdA], inCellMLL[nc][1][evIdA], inCellMLL[nc][2][evIdA],cellsk);
#endif
      cellRange[2 * k] = 1;
      if (inCellMLL[nc][k][evIdA] == cellsk) 
	{
	  inCellMLL[nc][k][evIdA] = 0;
#ifdef MD_LXYZ
	  *rkptr = -L2[k];
#else
	  *rkptr = -L2;
#endif
#ifdef MD_LXYZ
	  if (OprogStatus.useNNL)
	    nebrTab[evIdA].r[k] -= L[k];
#else
	  if (OprogStatus.useNNL)
	    nebrTab[evIdA].r[k] -= L;
#endif
	  OprogStatus.DR[evIdA][k]++;
	}

    }
  else
    { 
      cellRange[2 * k + 1] = -1;
      inCellMLL[nc][k][evIdA] = inCellMLL[nc][k][evIdA] - 1;
      if (inCellMLL[nc][k][evIdA] == -1) 
	{
	  inCellMLL[nc][k][evIdA] = cellsk - 1;
	  if (OprogStatus.useNNL)
#ifdef MD_LXYZ
	    nebrTab[evIdA].r[k] += L[k];
#else
	    nebrTab[evIdA].r[k] += L;
#endif

#ifdef MD_LXYZ
	  *rkptr = L2[k];
#else
	  *rkptr = L2;
#endif
	  OprogStatus.DR[evIdA][k]--;
	}
    }
}
extern void UpdateAtom(int i);
void PredictCellCross(int na, int nc);
void PredictCollMLL(int na, int nb, int nl);
void PredictCollMLL_NLL(int na, int nb);

void ProcessCellCrossingMLL(void)
{
  int kk, n, k, nl, numll;
  int nc, boxwall, nc_bw, nlcross_bw;//nlcoll, nlcoll_bw;
  int typei;

  UpdateAtom(evIdA);
  /* kk ci da la direzione lungo cui si sta realizzando il cell crossing */
  kk = evIdB - 100 - ATOM_LIMIT; 
  /* evIdC è semplicemente nc cioè ci dice il tipo di cella attraversata dalla particella evIdA, ossia
     nc = cella per interazine con tipo nc */
  nc = evIdC;

  typei = typeOfPart[evIdA];

  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;

  nl = get_linked_list_type(typei, nc);
 // printf("CELL CROSSING DI evIdA=%d nl=%d nc=%d typei=%d\n",evIdA, nl, nc, typei);
  //nc_bw = ;
#if 0
  if (iA == 0 && nc == 0)
    {
      nlcoll = 0;
      nlcross = 0;
      nlcoll_bw = 3;
      nlcross_bw = 2;
      nc_bw = 1;
    }
  else if (iA == 1 && nc == 0)
    { 
      nlcoll = 1;
      nlcross = 1;
      nlcoll_bw = 2;
      nlcross_bw = 3;
      nc_bw = 1; 
    }
  else if (iA == 0 && nc == 1)
    {
      nlcoll = 3;
      nlcross = 2;
      nlcoll_bw = 0;
      nlcross_bw = 0;
      nc_bw = 0;
    }
  else /* iA == 1 && nc == 1 */
    {
      nlcoll = 2;
      nlcross = 3;
      nlcoll_bw = 1;
      nlcross_bw = 1;
      nc_bw = 0;
    }
#endif 
  boxwall = check_boxwall(kk, nc, nl);
  /* questa condizione non si dovrebbe *mai* verificare, quindi 
   * in linea di principio le due linee che seguono potrebbero anche essere eliminate */
#if 0
  if (evIdA==870)
    {
      printf("[PROCCELLCROSS] type=%d inCell=%d %d %d\n", typeOfPart[evIdA], inCellMLL[nc][0][evIdA],inCellMLL[nc][1][evIdA],inCellMLL[nc][2][evIdA]);
      printf("[PROCESS CELL CROSS]kk=%d nc=%d nl=%d time=%.15G\n", kk, nc, nl, Oparams.time);
    }
#endif
  if (nc!=typei && boxwall)
    {
      
      printf("[PROCCELLCROSS] evIdA=%d\n",evIdA);
      printf("[PROCCELLCROSS] nc=%d and boxwall!!!! <===!!! typei=%d nl=%d kk=%d\n", nc, typei, nl, kk);
      exit(-1);
      return;
    }
  //printf("ProcellCellCross nl=%d nc=%d k=%d\n", nl, nc, k);
  /* NOTA: cellList[i] con 0 < i < Oparams.parnum è la cella in cui si trova la particella
   * i-esima mentre cellList[j] con 
   * Oparams.parnum <= j < cellsx*cellsy*cellsz+Oparams.parnum
   * è la prima particella che si trova nella cella j-esima
   */
  n = (inCellMLL[nc][2][evIdA] * cellsyMLL[nl] + inCellMLL[nc][1][evIdA])*cellsxMLL[nl] + 
    inCellMLL[nc][0][evIdA] + Oparams.parnum;
#if 0
  printf("nc=%d n=%d cellList[%d][%d]:%d\n",nc, n, nlcross, n, cellList[nlcross][n]);
  printf("vel=(%f,%f,%f) inCell= %d %d %d\n", vx[evIdA], vy[evIdA], vz[evIdA], inCell[nc][0][evIdA],inCell[nc][1][evIdA], inCell[nc][2][evIdA]);
#endif
  while (cellListMLL[nl][n] != evIdA) 
    n = cellListMLL[nl][n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellListMLL[nl][n] = cellListMLL[nl][evIdA];
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }

#if 0
  if (evIdA==511)
    printf("===> na=%d boxwall=%d time=%.15G\n",evIdA, boxwall, Oparams.time);
#endif
  if (boxwall)
    {
      /* se si è attraversata la scatola allora processa l'attraversamento per tutte le linked lists */ 
      for (nc_bw = 0; nc_bw < Oparams.ntypes; nc_bw++)
       	{
	  /* notare che se si tratta di un boxwall per forza di cose
	     deve essere nc = typeOfPart[na]. */
	  if (nc_bw == nc)
	    continue;
#ifdef MD_OPT_MULTLL
#if 0
	  if (!may_interact_all_type(nc, nc_bw))
	      continue;
#else
	  if (ignoreMLL[get_linked_list_type(nc, nc_bw)])
	    continue;
#endif
#endif 
	  nlcross_bw = get_linked_list_type(typei, nc_bw);
	  n = (inCellMLL[nc_bw][2][evIdA] * cellsyMLL[nlcross_bw] + inCellMLL[nc_bw][1][evIdA])*cellsxMLL[nlcross_bw] + 
	    inCellMLL[nc_bw][0][evIdA]
	    + Oparams.parnum;
	  while (cellListMLL[nlcross_bw][n] != evIdA) 
	    n = cellListMLL[nlcross_bw][n];
	  /* Eliminazione di evIdA dalla lista della cella n-esima della lista nl2 */
	  cellListMLL[nlcross_bw][n] = cellListMLL[nlcross_bw][evIdA];
	}	
    }
  switch (kk)
    {
    case 0: 
      docellcrossMLL(0, vx[evIdA], &(rx[evIdA]), cellsxMLL[nl], nc);
      break;
    case 1: 
      docellcrossMLL(1, vy[evIdA], &(ry[evIdA]), cellsyMLL[nl], nc);
      break;
    case 2:
      docellcrossMLL(2, vz[evIdA], &(rz[evIdA]), cellszMLL[nl], nc);
      break;
    }
  PredictCellCross(evIdA, nc);
  if (OprogStatus.useNNL==0)
    PredictCollMLL(evIdA, evIdB, nl);

  n = (inCellMLL[nc][2][evIdA] * cellsyMLL[nl] + inCellMLL[nc][1][evIdA])*cellsxMLL[nl] + 
    inCellMLL[nc][0][evIdA] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellListMLL[nl][evIdA] = cellListMLL[nl][n];
  cellListMLL[nl][n] = evIdA;
#if 0
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
#endif
#if 0
  printf("DOPO boxwall=%d nc=%d n=%d cellList[%d][%d]:%d\n",boxwall, nc, n, nlcross, n, cellList[nlcross][n]);
  printf("DOPO vel=(%f,%f,%f) inCell= %d %d %d\n", vx[evIdA], vy[evIdA], vz[evIdA], inCell[nc][0][evIdA],inCell[nc][1][evIdA], inCell[nc][2][evIdA]);
#endif
  if (boxwall)
    {
      for (nc_bw = 0; nc_bw < Oparams.ntypes; nc_bw++)
	{
	  if (nc_bw == nc)
	    continue;
#ifdef MD_OPT_MULTLL
#if 0
	  if (!may_interact_all_type(nc, nc_bw))
	      continue;
#else
	  if (ignoreMLL[get_linked_list_type(nc, nc_bw)])
	      continue;
#endif
#endif
	  nlcross_bw = get_linked_list(evIdA, nc_bw);
	  switch (kk)
	    {
	    case 0: 
	      docellcross2MLL(0, vx[evIdA], cellsxMLL[nlcross_bw], nc_bw);
	      break;
	    case 1: 
	      docellcross2MLL(1, vy[evIdA], cellsyMLL[nlcross_bw], nc_bw);
	      break;
	    case 2:
	      docellcross2MLL(2, vz[evIdA], cellszMLL[nlcross_bw], nc_bw);
	      break;
	    }
	  /* crossevtodel[0...Oparams.ntypes-1] deve contenere gli eventi di cell crossing per ogni tipo  
	     relativi alla particella evIdA, chiaramente crossevtode[typeOfPart[evIdA]][evIdA] sarà
	     sempre uguale al valore d'inizializzazione -1 */
	  if (crossevtodel[nc_bw][evIdA]!=-1)
	    {
	      //printf("DELETING CROSS EVENT evIdA=%d\n", evIdA);
	      DeleteEvent(crossevtodel[nc_bw][evIdA]);
	      crossevtodel[nc_bw][evIdA] = -1;
	    }
	  PredictCellCross(evIdA, nc_bw);
	  if (OprogStatus.useNNL==0)
	    PredictCollMLL(evIdA, evIdB, nlcross_bw);
	  n = (inCellMLL[nc_bw][2][evIdA] * cellsyMLL[nlcross_bw] + inCellMLL[nc_bw][1][evIdA])*cellsxMLL[nlcross_bw] + 
	    inCellMLL[nc_bw][0][evIdA] + Oparams.parnum;
	  /* Inserimento di evIdA nella nuova cella (head) */
	  cellListMLL[nlcross_bw][evIdA] = cellListMLL[nlcross_bw][n];
	  cellListMLL[nlcross_bw][n] = evIdA;
	}
    }
}
extern int locateHardWall(int na, int nplane, double tsup, double vecg[5], int ghw);
void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent); 

#if defined(MD_EDHEFLEX_WALL)
void PredictHardWall(int na, int evCode, double cctime)
{
  double vecg[5];
  int nl, typena, nplane=-1;
  /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  //for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  typena = typeOfPart[na];
  if (OprogStatus.hardwall)
    {
#if defined(MD_ABSORPTION) 
      if (OprogStatus.bufHeight > 0.0)
	{
#if !defined(MD_SPHERICAL_WALL)
	  if (vz[na] != 0.0)
	    {
#ifdef MD_LXYZ
	      hwcell = (L[2]-OprogStatus.bufHeight)*cellsz[typena]/L[2];
#else
	      hwcell = (L-OprogStatus.bufHeight)*cellsz[typena]/L;
#endif
#if 1
	      if (hwcell-inCell[typena][2][na] < 2)
		{
		  /* the semi-permeable plane is just one (nplane=0) */
		  if (locateHardWall(na, 0, cctime, vecg, 2))
		    {
		      rxC = vecg[0];
		      ryC = vecg[1];
		      rzC = vecg[2];
		      ScheduleEventBarr (na, ATOM_LIMIT+50, typena, 0, MD_WALL, vecg[4]);
		    }
		}
#endif
	    }
#endif
	}
#endif
      /* spostare in PredictColl con opportuno if (... )*/
      nl = get_linked_list(na, typena);
      //if (inCellMLL[typena][2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      //if (inCellMLL[typena][2][na] + cellRangeT[2 * 2 + 1] == cellszMLL[nl]) cellRangeT[2 * 2 + 1] = 0;
      /* ========= */
      if (inCellMLL[typena][2][na] == 0)
	nplane = 0;
      else if (inCellMLL[typena][2][na] == cellszMLL[nl]-1)
	nplane = 1;
      //if (na==955)
	//printf("na=%d cella:%d\n", na, inCell[2][na]);
      if (nplane!=-1 && locateHardWall(na, nplane, cctime, vecg, 1))
	{
	  rxC = vecg[0];
	  ryC = vecg[1];
	  rzC = vecg[2];
	  /* l'urto con il muro riguarda le LL relative all'interazione della particella con se stessa */
	  ScheduleEventBarr (na, ATOM_LIMIT + nplane, typena, 0, MD_WALL, vecg[4]);
	}
      else
	{
	  ScheduleEventBarr (na, ATOM_LIMIT + evCode, typena, 0, MD_EVENT_NONE, cctime);
	}
    }
  else
    {
      ScheduleEventBarr (na, ATOM_LIMIT + evCode, typena, 0, MD_EVENT_NONE, cctime);
    }
}
#endif
/* 21/05/2010: nc indica la linked lists da considerare ossia se nc=typeOfPart[na]=interazione con lo stesso tipo
   mentre valori diversi si riferiscono alle interazioni con gli altri tipi, 
   ad es. se ci sono 4 tipi la particella na è tipo 1 allora nc=1 sono le celle per l'interazione
   con lo stesso tipo mentre 0 = interazione con particelle di tipo 0, 2 interazione con particelle di tipo 2
   3 con quelle di tipo 3. */

void PredictCellCross(int na, int nc)
{
  int ignorecross[3], k, evCode, signDir[NDIM]={0,0,0}, nl;
  double tm[3], cctime=timbig;
  ignorecross[0] = ignorecross[1] = ignorecross[2] = 1;

  nl = get_linked_list(na, nc);
#ifdef MD_OPT_MULTLL
#if 0
  if (!may_interact_all_type(typeOfPart[na], nc) && typeOfPart[na] != nc)
    {
      //printf("ignoring typeOfPart[%d]=%d-nc=%d\n", na, typeOfPart[na], nc);    
      return;
    }
#else
  if (ignoreMLL[nl])
    return;
#endif
#endif
#if 0
  printf("[PredictCellCross ]time=%f nl=%d nc=%d n=%d inCell: %d %d %d cells: %d %d %d\n",
	 Oparams.time, nl, nc, na, inCell[nc][0][na], inCell[nc][1][na], inCell[nc][2][na],
	 cellsx[nl], cellsy[nl], cellsz[nl]);
#endif
  if (vz[na] != 0.0) 
    {
      if (vz[na] > 0.0) 
	{
	  signDir[2] = 0;/* direzione positiva */
	  if (nc != typeOfPart[na] && inCellMLL[nc][2][na]==cellszMLL[nl]-1)
	    ignorecross[2] = 1;
	  else
	    ignorecross[2] = 0;
	}
      else 
	{
	  signDir[2] = 1;/* direzione negativa */
	  if (nc != typeOfPart[na] && inCellMLL[nc][2][na]==0)
	    ignorecross[2] = 1;
	  else
	    ignorecross[2] = 0;
	}
#ifdef MD_LXYZ
#if defined(MD_EDHEFLEX_WALL) 
      if (ignorecross[2] || (OprogStatus.hardwall && ((signDir[2]==0 && inCellMLL[2][nc][na]==cellszMLL[nl]-1) || (signDir[2]==1 && inCellMLL[nc][2][na]==0))))
	tm[2] = timbig;
      else
	tm[2] = ((inCellMLL[nc][2][na] + 1 - signDir[2]) * L[2] /
		 cellszMLL[nl] - rz[na] - L2[2]) / vz[na];
#else
      if (ignorecross[2])
	tm[2] = timbig;
      else
	tm[2] = ((inCellMLL[nc][2][na] + 1 - signDir[2]) * L[2]/
	      	 cellszMLL[nc] - rz[na] - L2[2]) / vz[na];
#endif
#else
#if defined(MD_EDHEFLEX_WALL) 
      if (ignorecross[2] || (OprogStatus.hardwall && ((signDir[2]==0 && inCellMLL[nc][2][na]==cellszMLL[nl]-1) || (signDir[2]==1 && inCellMLL[nc][2][na]==0))))
	tm[2] = timbig;
      else
	tm[2] = ((inCellMLL[nc][2][na] + 1 - signDir[2]) * L /
		 cellszMLL[nl] - rz[na] - L2) / vz[na];
#else
      if (ignorecross[2])
	tm[2] = timbig;
      else
	tm[2] = ((inCellMLL[nc][2][na] + 1 - signDir[2]) * L /
	      	 cellszMLL[nl] - rz[na] - L2) / vz[na];
#endif
#endif

    } 
  else 
    tm[2] = timbig;
  /* end forcefield[k] != 0*/

  if (vx[na] != 0.0) 
    {
      if (vx[na] > 0.0) 
	{
	  if (nc != typeOfPart[na] && inCellMLL[nc][0][na]==cellsxMLL[nl]-1)
	    ignorecross[0] = 1;
	  else
	    ignorecross[0] = 0;
	  signDir[0] = 0;/* direzione positiva */
	}
      else
	{
	  if (nc != typeOfPart[na] && inCellMLL[nc][0][na]==0)
	    ignorecross[0] = 1;
	  else
	    ignorecross[0] = 0;
	  signDir[0] = 1;/* direzione negativa */
	}
      if (ignorecross[0])
	tm[0] = timbig;
      else
#ifdef MD_LXYZ
	tm[0] = ((inCellMLL[nc][0][na] + 1 - signDir[0]) * L[0] /
	      	 cellsxMLL[nl] - rx[na] - L2[0]) / vx[na];
#else
	tm[0] = ((inCellMLL[nc][0][na] + 1 - signDir[0]) * L /
	      	 cellsxMLL[nl] - rx[na] - L2) / vx[na];
#endif
    } 
  else 
    tm[0] = timbig;

  if (vy[na] != 0.0) 
    {
      if (vy[na] > 0.0) 
	{
	  if (nc != typeOfPart[na] && inCellMLL[nc][1][na]==cellsyMLL[nl]-1)
	    ignorecross[1] = 1;
	  else
	    ignorecross[1] = 0;
	  signDir[1] = 0;
	}
      else
	{
	  if (nc != typeOfPart[na] && inCellMLL[nc][1][na]==0)
	    ignorecross[1] = 1;
	  else
	    ignorecross[1] = 0;
	  signDir[1] = 1;
	}
      if (ignorecross[1])
	tm[1] = timbig;
      else
#ifdef MD_LXYZ
	tm[1] = ((inCellMLL[nc][1][na] + 1 - signDir[1]) * L[1] /
	      	 cellsyMLL[nl] - ry[na] - L2[1]) / vy[na];
#else
	tm[1] = ((inCellMLL[nc][1][na] + 1 - signDir[1]) * L /
	      	 cellsyMLL[nl] - ry[na] - L2) / vy[na];
#endif    
    } 
  else 
    tm[1] = timbig;
  /* ====== */
  /* Find minimum time */
  k = -1; /* giusto per dare un valore ed evitare una warning */
  if (tm[1] <= tm[2]) 
    {
      if (tm[0] <= tm[1]) 
	k = 0;
      else 
	k = 1;
    } 
  else
    {
      if (tm[0] <= tm[2]) 
	k = 0;
      else 
	k = 2;
    }
  /* Se un errore numerico fa si che tm[k] < 0 allora lo poniamo uguale a 0
   * (ved. articolo Lubachevsky) */
  if (tm[k]<0.0)
    {
      printf("tm[%d]: %.15G\n", k, tm[k]);
      tm[k] = 0.0;
      printf("k=%d nc=%d na=%d\n", k, nc, na);
#ifdef MD_LXYZ      
      printf("real cells: %d %d %d\n", (int)((rx[na] + L2[0]) * cellsxMLL[nl] / L[0]),
	     (int)((ry[na] + L2[1]) * cellsyMLL[nl] / L[1]), (int)((rz[na] + L2[2])  * cellszMLL[nl] / L[2]));
#else
      printf("real cells: %d %d %d\n", (int)((rx[na] + L2) * cellsxMLL[nl] / L),
	     (int)((ry[na] + L2) * cellsyMLL[nl] / L), (int)((rz[na] + L2)  * cellszMLL[nl] / L));
#endif
    }
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;// + 3*nc;
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 

  if (!ignorecross[k])
    {
#ifdef MD_EDHEFLEX_WALL 
      cctime = Oparams.time + tm[k];
      if (nc == typeOfPart[na] && OprogStatus.hardwall==1)
	{
	  PredictHardWall(na, evCode, cctime);
	}
      else
	{
#if 0
	  if (na==511)
	    {
	      printf("scheduling na=511 at t=%.15G k=%d curtime=%.15G\n", cctime, k, Oparams.time);
	      printf("r=%f %f %f v=%f %f %f\n", rx[na], ry[na], rz[na], vx[na], vy[na], vz[na]);
	      printf("inCell[511]==============%d %d %d\n", inCellMLL[0][0][511],inCellMLL[0][1][511],inCellMLL[0][2][511]);


	    }
#endif
#if 0
	  if (na==870 && nc!=typeOfPart[na])
	    {
	      printf("[PREDICT CELL CROSS]\n");
	      printf("nc=%d nl=%d typeOfPart[%d]=%d\n", nc, nl, na, typeOfPart[na]);
	      printf("cells = %d %d %d\n", cellsxMLL[nl], cellsyMLL[nl], cellszMLL[nl]);
	      printf("inCell=%d %d %d\n", inCellMLL[nc][0][na],inCellMLL[nc][1][na],inCellMLL[nc][2][na]);
	      printf("k=%d v=%f %f %f\n", k, vx[na], vy[na], vz[na]);
	      printf("ignorecross[k]=%d\n", ignorecross[k]);
	      printf("cctime=%.15G\n", cctime);
	      printf("[PREDICT CELL CROSS]-END\n");
	    }
#endif
	  ScheduleEventBarr (na, ATOM_LIMIT + evCode, nc, 0, MD_EVENT_NONE, cctime);
	}
#else
      ScheduleEventBarr (na, ATOM_LIMIT + evCode, nc, 0, MD_EVENT_NONE, Oparams.time + tm[k]);
#endif
    }
}
extern double *atomTime;
void rebuildMultipleLL(void)
{
  int nl, numll, maxnc, j, nc, n;
  
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  for (nl=0; nl < numll; nl++)
    {
#ifdef MD_OPT_MULTLL
      if (ignoreMLL[nl])
	continue;
#endif
      for (j = 0; j < cellsxMLL[nl]*cellsyMLL[nl]*cellszMLL[nl] + Oparams.parnum; j++)
	cellListMLL[nl][j] = -1;
    }
  //printf("rebuilding multiLL...");
  maxnc = Oparams.ntypes;
  for (nc=0; nc < maxnc; nc++)
    {
      /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
      for (n = 0; n < Oparams.parnum; n++)
	{
	  nl = get_linked_list(n, nc);
#ifdef MD_OPT_MULTLL
#if 0
	  if (!may_interact_all_type(typeOfPart[n], nc) && typeOfPart[n] != nc)
	    {
	      /* questa assegnazione non dovrebbe serivire...*/
	      atomTime[n] = Oparams.time;
	      continue;
	    }
#else
	  if (ignoreMLL[nl])
	    {
	      atomTime[n] = Oparams.time;
	      continue;
	    }
#endif
#endif

#ifdef MD_SPHERICAL_WALL
    	  if (n==sphWall)
    	    {
    	      cellListMLL[nl][sphWall]=-1;
    	      continue;
    	    }
	  if (n==sphWallOuter)
	    {
	      cellListMLL[nl][sphWallOuter]=-1;
	      continue;
	    }
#endif
	  atomTime[n] = Oparams.time;
#ifdef MD_LXYZ
	  inCellMLL[nc][0][n] =  (rx[n] + L2[0]) * cellsxMLL[nl] / L[0];
	  inCellMLL[nc][1][n] =  (ry[n] + L2[1]) * cellsyMLL[nl] / L[1];
	  inCellMLL[nc][2][n] =  (rz[n] + L2[2])  * cellszMLL[nl] / L[2];
#else
	  inCellMLL[nc][0][n] =  (rx[n] + L2) * cellsxMLL[nl] / L;
	  inCellMLL[nc][1][n] =  (ry[n] + L2) * cellsyMLL[nl] / L;
	  inCellMLL[nc][2][n] =  (rz[n] + L2) * cellszMLL[nl] / L;
#endif
#if 0
	  if (n==511)
	    {
	      printf("BOHBOH scheduling na=511 curtime=%.15G\n",Oparams.time);
	      printf("r=%f %f %f v=%f %f %f\n", rx[n], ry[n], rz[n], vx[n], vy[n], vz[n]);
	      printf("inCell[511]==============%d %d %d\n", inCellMLL[0][0][511],inCellMLL[0][1][511],inCellMLL[0][2][511]);


	    }
#endif
/*printf("nl=%d nc=%d n=%d inCell: %d %d %d cells: %d %d %d\n",
	    nl, nc, n, inCell[nc][0][n], inCell[nc][1][n], inCell[nc][2][n],
	    cellsx[nl], cellsy[nl], cellsz[nl]);
	   */
#if 0
	  if (inCell[0][n]>=cellsx ||inCell[1][n]>= cellsy||inCell[2][n]>= cellsz) 
	    {
	      printf("BOH?!?L:%f L2:%f n:%d rx[n]:%f\n", L, L2, n, rx[n]);
	      printf("(%d,%d,%d) (%d,%d,%d)\n",cellsx , cellsy,cellsz,
		     inCell[0][n],inCell[1][n], inCell[2][n]);
	    }
#endif	  
	  j = (inCellMLL[nc][2][n]*cellsyMLL[nl] + inCellMLL[nc][1][n])*cellsxMLL[nl] + 
	    inCellMLL[nc][0][n] + Oparams.parnum;
	  //if (n==511)
	    //printf("inCell==============%d %d %d\n", inCellMLL[nc][0][n],inCellMLL[nc][1][n],inCellMLL[nc][2][n]);
	  cellListMLL[nl][n] = cellListMLL[nl][j];
	  cellListMLL[nl][j] = n;
	}
    }
  //printf("...done\n")
}
extern int use_bounding_spheres(int na, int n);
extern int may_interact_all(int i, int j);
extern double *maxax;
extern double max_ax(int i);
extern int locate_contactSP(int i, int j, double shift[3], double t1, double t2, double *evtime, int *ata, int *atb, int *collCode);
extern int locate_contact(int i, int j, double shift[3], double t1, double t2, double vecg[5]);

void PredictCollMLL(int na, int nb, int nl) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq=0.0, dr[NDIM], dv[NDIM], shift[NDIM],  
	 b, d, t, tInt, vv, distSq, t1=0.0, t2=0.0;
  int ty1=-1, ty2=-1;
  int overlap, nc;
  int cellRangeT[2 * NDIM];
#ifdef MD_PATCHY_HE
  int ac, bc, collCode, collCodeOld, acHC, bcHC;
  double evtime, evtimeHC;
#endif
  double vecg[5];
  int iX, iY, iZ, jX, jY, jZ, k, n, typena;
#ifdef MD_SPHERICAL_WALL
  if (na==sphWall|| nb==sphWall)
    return;
  if (na==sphWallOuter|| nb==sphWallOuter)
    return;
#endif

  /* se il tipo della particella na è relativo alla linked lists nl 
     allora procedi nel predire le collisioni altrimenti non ha ovviamente
     senso ed esci da tale routine. */
  get_types_from_nl(nl, &ty1, &ty2);
  typena = typeOfPart[na]; 

  if (typena==ty1)
    nc = ty2;
  else if (typena==ty2)
    nc = ty1;
  else 
    return;
#ifdef MD_OPT_MULTLL
#if 0
  if (!may_interact_all_type(typena, nc) && typena != nc)
    {
      //printf("ignoring typeOfPart[%d]=%d-nc=%d\n", na, typeOfPart[na], nc);    
      return;
    }
#else
  if (ignoreMLL[nl])
    return;
#endif
#endif

#ifdef MD_EDHEFLEX_WALL
  if (OprogStatus.hardwall)
    {
      for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
      if (inCellMLL[nc][2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (inCellMLL[nc][2][na] + cellRangeT[2 * 2 + 1] == cellszMLL[nl]) cellRangeT[2 * 2 + 1] = 0;
    }
  else
    for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#else
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#endif
  //for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCellMLL[nc][2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
      if (jZ == -1) 
	{
	  jZ = cellszMLL[nl] - 1;    
#ifdef MD_LXYZ
	  shift[2] = - L[2];
#else
	  shift[2] = - L;
#endif
	} 
      else if (jZ == cellszMLL[nl]) 
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
	  jY = inCellMLL[nc][1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cellsyMLL[nl] - 1;    
#ifdef MD_LXYZ
	      shift[1] = -L[1];
#else
	      shift[1] = -L;
#endif
	    } 
	  else if (jY == cellsyMLL[nl]) 
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
	      jX = inCellMLL[nc][0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cellsxMLL[nl] - 1;    
#ifdef MD_LXYZ
		  shift[0] = - L[0];
#else
		  shift[0] = - L;
#endif
		} 
	      else if (jX == cellsxMLL[nl]) 
		{
		  jX = 0;   
#ifdef MD_LXYZ
		  shift[0] = L[0];
#else
		  shift[0] = L;
#endif
		}
	      n = (jZ *cellsyMLL[nl] + jY) * cellsxMLL[nl] + jX + Oparams.parnum;
	      for (n = cellListMLL[nl][n]; n > -1; n = cellListMLL[nl][n]) 
		{
		  /* se si tratta di una lista mista e il tipo della particella n è uguale 
		     a quello della particella na allora non considerare tale coppia
		     (questo è necessario poiché diversamente dalla silica qui si considera una solo linked list
		     mista e non due) */
		  if (nl >= Oparams.ntypes && typeOfPart[n] == typena)
		    {
		      continue;
		    }
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
#ifdef EDHE_FLEX
    		      if (!may_interact_all(na, n))
    			continue;
#endif
		      /* maxax[...] è il diametro dei centroidi dei due tipi
		       * di ellissoidi */

		      if (use_bounding_spheres(na, n))
			{
			  if (OprogStatus.targetPhi > 0)
			    {
			      sigSq = Sqr(max_ax(na)+max_ax(n)+OprogStatus.epsd);
			      //printf("maxax[n]=%.15G max_ax(n)=%.15G\n", maxax[n], max_ax(n));
			      //printf("sigSqHere=%.15G (%.15G)\n", sigSq, Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd));

			    }
			  else
			    {
#if defined(MD_POLYDISP) || defined(EDHE_FLEX)
			      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#else
			      if (na < parnumA && n < parnumA)
				sigSq = Sqr(maxax[na]+OprogStatus.epsd);
			      else if (na >= parnumA && n >= parnumA)
				sigSq = Sqr(maxax[na]+OprogStatus.epsd);
			      else
				sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#endif
			    }
			  tInt = Oparams.time - atomTime[n];
			  dr[0] = rx[na] - (rx[n] + vx[n] * tInt) - shift[0];	  
			  dv[0] = vx[na] - vx[n];
			  dr[1] = ry[na] - (ry[n] + vy[n] * tInt) - shift[1];
			  dv[1] = vy[na] - vy[n];
#ifdef MD_GRAVITY
			  dr[2] = rz[na] - 
			    (rz[n] + (vz[n] - 0.5 * Oparams.ggrav * tInt) * tInt) - shift[2];
			  dv[2] = vz[na] - (vz[n] - Oparams.ggrav * tInt);
#else
			  dr[2] = rz[na] - (rz[n] + vz[n] * tInt) - shift[2];
			  dv[2] = vz[na] - vz[n];

#endif
			  b = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
			  distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
			  vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
			  d = Sqr (b) - vv * (distSq - sigSq);

			  if (d < 0 || (b > 0.0 && distSq > sigSq)) 
			    {
			      /* i centroidi non collidono per cui non ci può essere
			       * nessun urto sotto tali condizioni */
			      continue;
			    }
			  if (vv==0.0)
			    {
			      if (distSq >= sigSq)
				continue;
			      /* la vel relativa è zero e i centroidi non si overlappano quindi
			       * non si possono urtare! */
			      t1 = t = 0;
			      t2 = 10.0;/* anche se sono fermi l'uno rispetto all'altro possono 
					   urtare ruotando */
			    }
			  else if (distSq >= sigSq)
			    {
			      t = t1 = - (sqrt (d) + b) / vv;
			      t2 = (sqrt (d) - b) / vv;
			      overlap = 0;
			    }
			  else 
			    {
			      t2 = t = (sqrt (d) - b) / vv;
			      t1 = 0.0; 
			      overlap = 1;
			      MD_DEBUG(printf("altro d=%f t=%.15f\n", d, (-sqrt (d) - b) / vv));
			      MD_DEBUG(printf("vv=%f dv[0]:%f\n", vv, dv[0]));
			    }
			  //t += Oparams.time; 
			  t2 += Oparams.time;
			  t1 += Oparams.time;
			}
		      else
			{
			  t1 = 0.0;
			  t2 = timbig;
			}
#if 0 
			  calcDist(Oparams.time, na, n, shift, r1, r2);
			  continue;
			  exit(-1);
#endif
#ifdef MD_PATCHY_HE
		      evtime = t2;
		      collCode = MD_EVENT_NONE;
		      rxC = ryC = rzC = 0.0;
		      collCodeOld = collCode;
		      evtimeHC = evtime;
		      acHC = ac = 0;
		      bcHC = bc = 0;
#ifdef EDHE_FLEX
		      if (OprogStatus.targetPhi <= 0.0)
			{
			  /* during growth disable interactions between spots */
			  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			    {
			      collCode = MD_EVENT_NONE;
			    }
			}
		      MD_DEBUG(if (collCode!=MD_EVENT_NONE) check_these_bonds(na, n, shift, evtime));
#else
		      if (OprogStatus.targetPhi <=0 && ((na < Oparams.parnumA && n >= Oparams.parnumA)|| 
							(na >= Oparams.parnumA && n < Oparams.parnumA)))
			{
			  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			    {
			      collCode = MD_EVENT_NONE;
			    }
			}
#endif
		      if (collCode!=MD_EVENT_NONE)
			t2 = evtime+1E-7;
#ifdef ED_HESW
		      /* per ora assumiamo un solo corpo rigido extra che serve per avere 
			 la buca quadrata, l'urto tra le buche ellissoidali viene considerato
			 come un urto tra gli spot 0 (per semplicità) */
		      if (typesArr[typeOfPart[na]].nhardobjs == 1 &&
			  typesArr[typeOfPart[n]].nhardobjs == 1)
			{
			  if (!locate_contact_hesw(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
			    { 
			      collCode = MD_EVENT_NONE;
			    }
			}
#endif

		      if (locate_contact(na, n, shift, t1, t2, vecg))
			{
			  if (collCode == MD_EVENT_NONE || (collCode!=MD_EVENT_NONE && vecg[4] <= evtime))
			    {
			      collCode = MD_CORE_BARRIER;
			      ac = bc = 0;
			      evtime = vecg[4];
			      rxC = vecg[0];
			      ryC = vecg[1];
			      rzC = vecg[2];
			    }
#ifdef MD_SAVE_DISTANCE
    			  printf("found collision exiting...\n");
    			  //exit(-1);
#endif
			}
		      else
			{
			  if (collCode == MD_EVENT_NONE)
			    continue;
			}

		      t = evtime;
#else
		      if (!locate_contact(na, n, shift, t1, t2, vecg))
		      	continue;

		      rxC = vecg[0];
		      ryC = vecg[1];
		      rzC = vecg[2];
		      MD_DEBUG(printf("A x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)-B x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)",
				      rx[na], ry[na], rz[na], vx[na], vy[na], vz[na],
				      rx[n], ry[n], rz[n], vx[n], vy[n], vz[n]));
		      t = vecg[4];

#endif
#ifdef MD_PATCHY_HE
		      if (t < Oparams.time)
			{
#if 1
			  printf("time:%.15f \n",Oparams.time);
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = Oparams.time;
			}

#else
		      if (t < 0)
			{
#if 1
			  printf("time:%.15f tInt:%.15f\n", Oparams.time,
				 tInt);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
	     					      Sqr(dr[2]))-1.0 );
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = 0;
			}
#endif
		      /* il tempo restituito da newt() è già un tempo assoluto */
#ifdef MD_PATCHY_HE
		      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
		      //printf("time: %f Adding collision (%d,%d)-(%d,%d)\n", t, na, ac, n, bc);

#else
		      ScheduleEvent (na, n, t);
#endif
		    }
		} 
	    }
	}
    }
}

void PredictEventMLL(int na, int nb)
{
  int nc, nl, numll, t1=-1, t2=-1;
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=0;
#endif
  for (nc = 0; nc < Oparams.ntypes; nc++)
    {
      PredictCellCross(na, nc);
    }
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;

#ifdef MD_SPHERICAL_WALL
  /* inner wall */
  locate_spherical_wall(na, 0);
  /* outer wall */
  locate_spherical_wall(na, 1);
#endif

  for (nl = 0; nl < numll; nl++)
    {
      get_types_from_nl(nl, &t1, &t2);
      if (typeOfPart[na]!=t1 && typeOfPart[na]!=t2)
	continue;
      PredictCollMLL(na, nb, nl);
    }
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=1;
#endif
}
#ifdef MD_EDHEFLEX_OPTNNL
void rebuildMultipleLL_NLL(int nl)
{
  /* N.B. Se si usano le NNL ottimizzate il centro di massa geometrico delle NNL non coincide con il 
     centro di massa degli ellissoidi, in tale caso quindi le linked lists degli ellissoidi vengono usate 
     per avere una sovrastima in caso di urti con la parete dura e per mantenere gli ellissoidi nella
     first box. */
  //double L2, invL;
#ifdef MD_LXYZ
  int kk;
#endif
  int t1=-1, t2=-1;
  int j, n, numll;
#ifdef MD_LXYZ
  double invL[NDIM], L2[NDIM];
#else
  double invL, L2;
#endif
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
#ifdef MD_LXYZ
  for (kk = 0; kk < 3; kk++)
    {
      L2[kk] = 0.5 * L[kk];
      invL[kk] = 1.0/L[kk];
    }
#else
  L2 = 0.5 * L;
  invL = 1.0/L;
#endif
  set_cells_size();
  for (j = 0; j < cellsxMLL[nl]*cellsyMLL[nl]*cellszMLL[nl] + Oparams.parnum; j++)
    cellList_NNL[j] = -1;
    
  /* NOTA: rcut delle LL per le NNL e' uguale a quello delle LL per gli ellissoidi
     ma quest'ultimo potrebbe essere scelto ad hoc. Inoltre le LL per gli ellissoidi 
     vengono solo usate per avere una upper limit per il tempo di collisione contro la pareti dure */  
  /* rebuild event calendar */

  get_types_from_nl(nl, &t1, &t2); 
  for (n = 0; n < Oparams.parnum; n++)
    {
      if (typeOfPart[n]!=t1 && typeOfPart[n]!=t2)
	continue;
#ifdef MD_SPHERICAL_WALL
      if (n==sphWall)
	{
	  cellList_NNL[sphWall]=-1;
	  continue;
	}
      if (n==sphWallOuter)
	{
	  cellList_NNL[sphWallOuter]=-1;
	  continue;
	}
#endif
      /* reduce to first box */
#ifdef MD_LXYZ
      rxNNL[n] = nebrTab[n].r[0] - L[0]*rint(nebrTab[n].r[0]*invL[0]);
      ryNNL[n] = nebrTab[n].r[1] - L[1]*rint(nebrTab[n].r[1]*invL[1]);
#ifdef MD_EDHEFLEX_WALL
      if (!OprogStatus.hardwall)
	rzNNL[n] = nebrTab[n].r[2] - L[2]*rint(nebrTab[n].r[2]*invL[2]);
      else
	rzNNL[n] = nebrTab[n].r[2];
#else
      rzNNL[n] = nebrTab[n].r[2] - L[2]*rint(nebrTab[n].r[2]*invL[2]);
#endif
      inCell_NNL[0][n] =  (rxNNL[n] + L2[0]) * cellsxMLL[nl] / L[0];
      inCell_NNL[1][n] =  (ryNNL[n] + L2[1]) * cellsyMLL[nl] / L[1];
      inCell_NNL[2][n] =  (rzNNL[n] + L2[2]) * cellszMLL[nl] / L[2];
#else
      rxNNL[n] = nebrTab[n].r[0] - L*rint(nebrTab[n].r[0]*invL);
      ryNNL[n] = nebrTab[n].r[1] - L*rint(nebrTab[n].r[1]*invL);
#ifdef MD_EDHEFLEX_WALL
      if (!OprogStatus.hardwall)
	rzNNL[n] = nebrTab[n].r[2] - L*rint(nebrTab[n].r[2]*invL);
      else
	rzNNL[n] = nebrTab[n].r[2];
#else
      rzNNL[n] = nebrTab[n].r[2] - L*rint(nebrTab[n].r[2]*invL);
#endif
      //printf("com %d = %f %f %f nebrtag.r=%f %f %f\n", n, rx[n], ry[n], rx[n], rxNNL[n], ryNNL[n],
	//    rzNNL[n]);
      inCell_NNL[0][n] =  (rxNNL[n] + L2) * cellsxMLL[nl] / L;
      inCell_NNL[1][n] =  (ryNNL[n] + L2) * cellsyMLL[nl] / L;
#ifdef MD_GRAVITY
      inCell_NNL[2][n] =  (rzNNL[n] + Lz2) * cellszMLL[nl] / (Lz+OprogStatus.extraLz);
#else
      inCell_NNL[2][n] =  (rzNNL[n] + L2)  * cellszMLL[nl] / L;
#endif
#endif
      j = (inCell_NNL[2][n]*cellsyMLL[nl] + inCell_NNL[1][n])*cellsxMLL[nl] + 
       inCell_NNL[0][n] + Oparams.parnum;
      cellList_NNL[n] = cellList_NNL[j];
      cellList_NNL[j] = n;
    }
}
#endif
extern void check_nnl_size(int na);
extern double calcDistNegNNLoverlap(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha, double *vecgsup, int calcguess);
extern double calcDistNegNNLoverlapPlane(double t, double t1, int i, int j, double shift[3]);

void BuildNNL_MLL(int na, int nl) 
{
  double shift[NDIM];
  int kk;
  double dist; 
  int cellsx, cellsy, cellsz;
  int *inCellL[3], *cellListL, t1=-1, t2=-1, nc;
#ifndef MD_NNLPLANES
  double vecgsup[8], alpha;
#endif
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cels[NDIM];*/
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
  /* reset list only on first call*/
  if (nl==0)
    nebrTab[na].len = 0;
#ifdef MD_SPHERICAL_WALL
  if (na==sphWall || na==sphWallOuter)
    return;
#endif
 for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (kk=0; kk < 3; kk++)
    shift[kk] = 0;
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

#ifdef MD_OPT_MULTLL
#if 0
  get_types_from_nl(nl, &t1, &t2);
  if (t1==typeOfPart[na])
    nc = t2;
  else if (t2==typeOfPart[na])
    nc = t1;
  else
    return;
  if (!may_interact_all_type(typeOfPart[na], nc) && nc!=typeOfPart[na])
    {
      return;
    }
#else
  if (ignoreMLL[nl])
    return;
#endif
#endif
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=0;
#endif
 
#ifdef MD_EDHEFLEX_OPTNNL
  if (OprogStatus.optnnl)
    {
      get_types_from_nl(nl, &t1, &t2);
      if (t1!=typeOfPart[na]&&t2!=typeOfPart[na])
        return;
      for (k=0; k < 3; k++)
	inCellL[k] = inCell_NNL[k];
      cellListL = cellList_NNL;
    }
  else
    {
      get_types_from_nl(nl, &t1, &t2);
      if (t1==typeOfPart[na])
	nc = t2;
      else if (t2==typeOfPart[na])
	nc = t1;
      else
	return;
      for (k=0; k < 3; k++)
	inCellL[k] = inCellMLL[nc][k];
      cellListL = cellListMLL[nl];
    }
#else
  get_types_from_nl(nl, &t1, &t2);
  if (t1==typeOfPart[na])
    nc = t2;
  else if (t2==typeOfPart[na])
    nc = t1;
  else
    return;
  for (k=0; k < 3; k++)
    inCellL[k] = inCellMLL[nc][k];
  cellListL = cellListMLL[nl];
#endif
  cellsx = cellsxMLL[nl];
  cellsy = cellsyMLL[nl];
  cellsz = cellszMLL[nl];
#if defined(MD_EDHEFLEX_WALL)
  /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (OprogStatus.hardwall)
    {
      if (inCellL[2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
      if (inCellL[2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
    }
#endif
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCellL[2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodic boundary condition along z if gravitational
       * fiels is not present */
      if (jZ == -1) 
	{
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
	  jY = inCellL[1][na] + iY;    
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
	      jX = inCellL[0][na] + iX;    
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
	      n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	      for (n = cellListL[n]; n > -1; n = cellListL[n]) 
		{
		  if (n != na)// && n != nb && (nb >= -1 || n < na)) 
		    {
#ifdef EDHE_FLEX
		      if (!may_interact_all(na, n))
			continue;
#endif
		      //assign_bond_mapping(na, n);
		      //dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
#ifdef MD_NNLPLANES
		      dist = calcDistNegNNLoverlapPlane(Oparams.time, 0.0, na, n, shift); 
#else
		      dist = calcDistNegNNLoverlap(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecgsup, 1); 
#endif
		      /* 0.1 è un buffer per evitare problemi, deve essere un parametro 
		       * in OprogStatus */
		      if (dist < 0)
			{
 			  nebrTab[na].list[nebrTab[na].len] = n;
			  //for (kk=0; kk < 3; kk++)
			  //nebrTab[na].shift[nebrTab[na].len][kk] = shift[kk];
			  nebrTab[na].len++;
			  check_nnl_size(na);
			}
		    }
		} 
	    }
	}
    }
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=1;
#endif
 
}
void BuildAllNNL_MLL_one(int i)
{
  int nl, numll;
  
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  for (nl = 0; nl < numll; nl++)
    {
#ifdef MD_EDHEFLEX_OPTNNL
      if (OprogStatus.optnnl)
	rebuildMultipleLL_NLL(nl);
#endif
      BuildNNL_MLL(i, nl);
    }
}
void BuildAllNNL_MLL(void)
{
  int nl, i, numll;
  
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;

  for (nl = 0; nl < numll; nl++)
    {
#ifdef MD_EDHEFLEX_OPTNNL
      if (OprogStatus.optnnl)
	rebuildMultipleLL_NLL(nl);
#endif
      for (i=0; i < Oparams.parnum; i++)
	{
	  BuildNNL_MLL(i, nl);
	}
    }
}

void PredictCollMLL_NLL(int na, int nb)
{
  int i, n;
  double vecg[5], shift[3], t1, t2, t;
  double sigSq, tInt, d, b, vv, dv[3], dr[3], distSq;
  int overlap;
#ifdef MD_PATCHY_HE
  int ac, bc, collCode, collCodeOld, acHC, bcHC;
  double evtime, evtimeHC;
#endif

#ifdef MD_SPHERICAL_WALL
  locate_spherical_wall(na, 0);
  locate_spherical_wall(na, 1);
#endif 
  /* NOTA: nel caso di attraversamento di una cella non deve predire le collisioni (visto che in tal caso stiamo 
     usando le NNL */
  if (nb >= ATOM_LIMIT+2*NDIM)
    return;

  for (i=0; i < nebrTab[na].len; i++)
    {
      n = nebrTab[na].list[i]; 
#if 0
      if (na==35 || na==16)
	printf("[PredictEventNNL] na=%d n=%d nb=%d\n",na,  n, nb);
#endif
      if (!(n != na && n!=nb && (nb >= -1 || n < na)))
	continue;
      	
      //
      // for (kk=0; kk < 3; kk++)
      //	shift[kk] = nebrTab[na].shift[i][kk];
#ifdef MD_LXYZ
      shift[0] = L[0]*rint((rx[na]-rx[n])/L[0]);
      shift[1] = L[1]*rint((ry[na]-ry[n])/L[1]);
#ifdef MD_EDHEFLEX_WALL
      if (!OprogStatus.hardwall)
	shift[2] = L[2]*rint((rz[na]-rz[n])/L[2]);
      else
	shift[2] = 0.0;
#else
      shift[2] = L[2]*rint((rz[na]-rz[n])/L[2]);
#endif
#else
      shift[0] = L*rint((rx[na]-rx[n])/L);
      shift[1] = L*rint((ry[na]-ry[n])/L);
#ifdef MD_EDHEFLEX_WALL
      if (!OprogStatus.hardwall)
	shift[2] = L*rint((rz[na]-rz[n])/L);
      else
	shift[2] = 0.0;
#else
      shift[2] = L*rint((rz[na]-rz[n])/L);
#endif
#endif
      /* maxax[...] è il diametro dei centroidi dei due tipi
       * di ellissoidi */
      if (use_bounding_spheres(na, n))
	{
	  if (OprogStatus.targetPhi > 0)
	    {
	      sigSq = Sqr(max_ax(na)+max_ax(n)+OprogStatus.epsd);
	    }
	  else
	    {
#ifdef MD_POLYDISP
	      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#else
#ifdef EDHE_FLEX
	      sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
	      //printf("max=(%d)%.15G (%d)%.15G\n", n, maxax[n], na, maxax[na]);
#else
	      if (na < parnumA && n < parnumA)
		sigSq = Sqr(maxax[na]+OprogStatus.epsd);
	      else if (na >= parnumA && n >= parnumA)
		sigSq = Sqr(maxax[na]+OprogStatus.epsd);
	      else
		sigSq = Sqr((maxax[n]+maxax[na])*0.5+OprogStatus.epsd);
#endif
#endif
	    }
	  tInt = Oparams.time - atomTime[n];
	  dr[0] = rx[na] - (rx[n] + vx[n] * tInt) - shift[0];	  
	  dv[0] = vx[na] - vx[n];
	  dr[1] = ry[na] - (ry[n] + vy[n] * tInt) - shift[1];
	  dv[1] = vy[na] - vy[n];
#ifdef MD_GRAVITY
	  dr[2] = rz[na] - 
	    (rz[n] + (vz[n] - 0.5 * Oparams.ggrav * tInt) * tInt) - shift[2];
	  dv[2] = vz[na] - (vz[n] - Oparams.ggrav * tInt);
#else
	  dr[2] = rz[na] - (rz[n] + vz[n] * tInt) - shift[2];
	  dv[2] = vz[na] - vz[n];
#endif
	  b = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
	  distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
	  vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
	  d = Sqr (b) - vv * (distSq - sigSq);
	  if (d < 0 || (b > 0.0 && distSq > sigSq)) 
	    {
	      /* i centroidi non collidono per cui non ci può essere
	       * nessun urto sotto tali condizioni */
	      continue;
	    }
	  if (vv==0.0)
	    {
	      if (distSq >= sigSq)
		{
		  continue;
		}
	      /* la vel relativa è zero e i centroidi non si overlappano quindi
	       * non si possono urtare! */
	      t1 = t = 0;
	      t2 = 10.0;/* anche se sono fermi l'uno rispetto all'altro possono 
			   urtare ruotando */
	    }
	  else if (distSq >= sigSq)
	    {
	      t = t1 = - (sqrt (d) + b) / vv;
	      t2 = (sqrt (d) - b) / vv;
	      overlap = 0;
	    }
	  else 
	    {
	      MD_DEBUG(printf("Centroids overlap!\n"));
	      t2 = t = (sqrt (d) - b) / vv;
	      t1 = 0.0; 
	      overlap = 1;
	      MD_DEBUG(printf("altro d=%f t=%.15f\n", d, (-sqrt (d) - b) / vv));
	      MD_DEBUG(printf("vv=%f dv[0]:%f\n", vv, dv[0]));
	    }
	  MD_DEBUG(printf("t=%f curtime: %f b=%f d=%f\n", t, Oparams.time, b ,d));
	  MD_DEBUG(printf("dr=(%f,%f,%f) sigSq: %f", dr[0], dr[1], dr[2], sigSq));
	  //t += Oparams.time; 
	  t2 += Oparams.time;
	  t1 += Oparams.time;
	}
     else
       {
	 t1 = 0.0;
	 t2 = timbig;
       } 
#if 0
      tnnl = min(nebrTab[na].nexttime,nebrTab[n].nexttime);
      if (tnnl < t2)
	t2 = tnnl;
#else      
      /* WARNING: OprogStatus.h è un buffer di sicurezza */
      if (nextNNLrebuild < t2)
	t2 = nextNNLrebuild + OprogStatus.h;
#endif
      // t1 = Oparams.time;
      // t2 = nebrTab[na].nexttime;//,nebrTab[n].nexttime);

#ifdef MD_PATCHY_HE
      evtime = t2;
      collCode = MD_EVENT_NONE;
      rxC = ryC = rzC = 0.0;
      collCodeOld = collCode;
      evtimeHC = evtime;
      acHC = ac = 0;
      bcHC = bc = 0;
#ifdef EDHE_FLEX
      if (OprogStatus.targetPhi <=0 && typesArr[typeOfPart[na]].nspots > 0 && typesArr[typeOfPart[n]].nspots > 0)
	{
	  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
	    {
	      collCode = MD_EVENT_NONE;
	    }
	}
#else
      if (OprogStatus.targetPhi <=0 && ((na < Oparams.parnumA && n >= Oparams.parnumA)|| 
					(na >= Oparams.parnumA && n < Oparams.parnumA)))
	{
	  if (!locate_contactSP(na, n, shift, t1, t2, &evtime, &ac, &bc, &collCode))
	    {
	      collCode = MD_EVENT_NONE;
	    }
	}
#endif
      if (collCode!=MD_EVENT_NONE)
	t2 = evtime+1E-7;
      if (locate_contact(na, n, shift, t1, t2, vecg))
	{
	  if (collCode == MD_EVENT_NONE || (collCode!=MD_EVENT_NONE && vecg[4] <= evtime))
	    {
	      collCode = MD_CORE_BARRIER;
	      ac = bc = 0;
	      evtime = vecg[4];
	      rxC = vecg[0];
	      ryC = vecg[1];
	      rzC = vecg[2];
	    }
#ifdef MD_SAVE_DISTANCE
	  printf("found collision exiting...\n");
	  //exit(-1);
#endif
	}
      else
	{
	  if (collCode == MD_EVENT_NONE)
	    continue;
	}

      t = evtime;
#else
      if (!locate_contact(na, n, shift, t1, t2, vecg))
	{
	  continue;
	}
      rxC = vecg[0];
      ryC = vecg[1];
      rzC = vecg[2];
      t = vecg[4];
#endif
#ifdef MD_PATCHY_HE
      //printf("Scheduling collision between %d and %d ac=%d bc=%d at t=%.15G\n", na, n, ac, bc, t);
      ScheduleEventBarr (na, n,  ac, bc, collCode, t);
#else
      ScheduleEvent (na, n, t);
#endif
    }
}
void PredictEventNNL_MLL(int na, int nb)
{
  int nc;
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=0;
#endif
  for (nc = 0; nc < Oparams.ntypes; nc++)
    {
      PredictCellCross(na, nc);
      //PredictCellCross(evIdB, nc);
    }
#if 0
  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  for (nl = 0; nl < numll; nl++)
    {
      get_types_from_nl(nl, &t1, &t2);
      if (typeOfPart[na]!=t1 && typeOfPart[na]!=t2)
    	continue;
      PredictCollMLL_NLL(na, nb);
    }
#endif     
  PredictCollMLL_NLL(na, nb);
#if 0
  for (nl = 0; nl < numll; nl++)
    {
      PredictCollMLL_NLL(evIdB, evIdA, nl);
    }
#endif
#ifdef MD_SOLVENT_NOHW
  if (typeOfPart[na]==4)
    OprogStatus.hardwall=1;
#endif
}
extern int calcdist_retcheck;
extern double calcDistNeg(double t, double t1, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
	      		  double *vecgsup, int calcguess);

double get_min_dist_MLL(int na, int *jmin, double *rCmin, double *rDmin, double *shiftmin) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  int cellsx, cellsy, cellsz, numll;
  int t1=-1, t2=-1, nc;
  double distMin=1E10,dist,vecg[8], alpha, shift[3];
  /*double cells[NDIM];*/
  int kk, nl;
  double r1[3], r2[3];
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;

  for (nl = 0; nl < numll; nl++)
    {
      if (ignoreMLL[nl])
	continue;

      calcdist_retcheck = 0;

      get_types_from_nl(nl, &t1, &t2);
      if (t1==typeOfPart[na])
	nc = t2;
      else if (t2==typeOfPart[na])
	nc = t1;
      else
	continue;

      cellsx = cellsxMLL[nl];
      cellsy = cellsyMLL[nl];
      cellsz = cellszMLL[nl];
#if defined(MD_EDHEFLEX_WALL)
      /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
      if (OprogStatus.hardwall)
	{
	  if (inCellMLL[nc][2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
	  if (inCellMLL[nc][2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
	}
#endif
      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCellMLL[nc][2][na] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
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
	      jY = inCellMLL[nc][1][na] + iY;    
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
		  jX = inCellMLL[nc][0][na] + iX;    
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
		  n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
		  for (n = cellListMLL[nl][n]; n > -1; n = cellListMLL[nl][n]) 
		    {
		      if (n!=na) 
			{
			  dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
			  if (calcdist_retcheck)
			    continue;
#if 0
			  if ((na==125||na==15) && (n==15||n==125))
			    printf("$$$$ dist: %.12G\n", dist);
#endif
			  if (*jmin == -1 || dist<distMin)
			    {
			      distMin = dist;
			      for (kk = 0; kk < 3; kk++)
				{
				  rCmin[kk] = r1[kk];
				  rDmin[kk] = r2[kk];
				  shiftmin[kk] = shift[kk];
				}
			      *jmin = n;
			    }
			}
		    } 
		}
	    }
	}
    }
  return distMin;
}


#ifdef MD_SPHERICAL_WALL
void reinsert_protein_MLL(int protein, int oldtype)
{
  int nl, n, ty1=-1, ty2=-1, typena, nc, numll;
  int j, kk, k, i, ii;

  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  typena = typeOfPart[protein];
  
  for (nl = 0; nl < numll; nl++)
    {
#ifdef MD_OPT_MULTLL
      if (ignoreMLL[nl])
	continue;
#endif
      get_types_from_nl(nl, &ty1, &ty2);
      //printf("[REINSERT PROTEIN MLL] resinserting protein %d typena=%d nl=%d type1=%d type2=%d\n", protein, typena, nl, ty1, ty2);
      /* N.B. 27/05/2010: al momento dell'assorbimento la particella è non-ghost ma appena dopo è ghost quindi 
	 va tolta con il vecchio tipo e reinserita con il nuovo, analogamente nel caso di urto
	 con la membrana semipermeabile la particella da ghost diventa non-ghost e questo richiede un aggiornamento 
	 delle linked lists. */
      if (ty1==oldtype)
	nc = ty2;
      else if (ty2==oldtype)
	nc = ty1;
      else 
	continue;
      n = (inCellMLL[nc][2][protein] * cellsyMLL[nl] + inCellMLL[nc][1][protein] )*cellsxMLL[nl] + inCellMLL[nc][0][protein]
	+ Oparams.parnum;

      while (cellListMLL[nl][n] != protein) 
	n = cellListMLL[nl][n];
      /* Eliminazione di protein dalla lista della cella n-esima */
      cellListMLL[nl][n] = cellListMLL[nl][protein];
    }
  /* inserimento nella nuova cella */
  for (nl = 0; nl < numll; nl++)
    {
      if (ignoreMLL[nl])
	continue;
      get_types_from_nl(nl, &ty1, &ty2);
      if (ty1==typena)
	nc = ty2;
      else if (ty2==typena)
	nc = ty1;
      else 
	continue;
#ifdef MD_LXYZ
      inCellMLL[nc][0][protein] =  (rx[protein] + L2[0]) * cellsxMLL[nl] / L[0];
      inCellMLL[nc][1][protein] =  (ry[protein] + L2[1]) * cellsyMLL[nl] / L[1];
      inCellMLL[nc][2][protein] =  (rz[protein] + L2[2]) * cellszMLL[nl] / L[2];
#else
      inCellMLL[nc][0][protein] =  (rx[protein] + L2) * cellsxMLL[nl] / L;
      inCellMLL[nc][1][protein] =  (ry[protein] + L2) * cellsyMLL[nl] / L;
#ifdef MD_GRAVITY
      inCellMLL[nc][2][protein] =  (rz[protein] + Lz2) * cellszMLL[nl] / (Lz+OprogStatus.extraLz);
#else
      inCellMLL[nc][2][protein] =  (rz[protein] + L2)  * cellszMLL[nl] / L;
#endif
#endif
      j = (inCellMLL[nc][2][protein]*cellsyMLL[nl] + inCellMLL[nc][1][protein])*cellsxMLL[nl] + 
	inCellMLL[nc][0][protein] + Oparams.parnum;
      cellListMLL[nl][protein] = cellListMLL[nl][j];
      cellListMLL[nl][j] = protein;
      //printf("QUI2\n");
    }
  if (OprogStatus.useNNL)
    {
      //listtmp = malloc(sizeof(int)*OprogStatus.nebrTabFac);
      /* rimuove protein da tutte le NNL delle particelle nella NNL 
	 di protein attuale */
      for (i = 0; i < Oparams.parnum; i++)
	{
	  kk=0;
	  for (k = 0; k < nebrTab[i].len; k++)
	    {
	      n = nebrTab[i].list[k];
	      if (n!=protein)
		{
		  listtmp[kk] = nebrTab[i].list[k]; 
		  kk++;
		}
	    }
	  nebrTab[i].len = kk;
	  for (k = 0; k < nebrTab[i].len; k++) 
	    {
	      nebrTab[i].list[k] = listtmp[k];
	    }
	}
      BuildAllNNL_MLL_one(protein);
      /* ricostruisce le NNL per tutte le particelle nella NNL di protein */
      for (ii=0; ii < nebrTab[protein].len; ii++)
	{
	  n = nebrTab[protein].list[ii]; 
	  BuildAllNNL_MLL_one(n);
	}
      //free(listtmp);
    }
}
#endif
extern void remove_bond(int na, int n, int a, int b);
void check_all_bonds_MLL(void)
{
  int warn, j;
  int i, cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k;
  int *cellListL, *inCellL[3], cellsx, cellsy, cellsz;
  int nc, t1=-1, t2=-1, nl, numll;
  double shift[3]={0.0,0.0,0.0};

  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  warn = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
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

      for (nl = 0; nl < numll; nl++)
	{
	  //printf("find bonds one i=%d nl=%d ignore=%d numll=%d\n",i,  nl, ignoreMLL[nl], numll);
	  if (ignoreMLL[nl])
	    continue;
	  get_types_from_nl(nl, &t1, &t2);
	  if (t1==typeOfPart[i])
	    nc = t2;
	  else if (t2==typeOfPart[i])
	    nc = t1;
	  else
	    continue;
	  for (k=0; k < 3; k++)
	    inCellL[k] = inCellMLL[nc][k];
	  cellListL = cellListMLL[nl];
	  cellsx = cellsxMLL[nl];
	  cellsy = cellsyMLL[nl];
	  cellsz = cellszMLL[nl];
	  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_EDHEFLEX_WALL
	  if (OprogStatus.hardwall==1)
	    {
	      if (inCellL[2][i] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
	      if (inCellL[2][i] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
	    }
#endif 
	  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	    {
	      jZ = inCellL[2][i] + iZ;    
	      shift[2] = 0.;
	      /* apply periodico boundary condition along z if gravitational
	       * fiels is not present */
	      if (jZ == -1) 
		{
		  //printf("BOHHHH\n");
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
		  jY = inCellL[1][i] + iY;    
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
		      jX = inCellL[0][i] + iX;    
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
		      for (j = cellListL[j]; j > -1; j = cellListL[j]) 
			{
			  if (i == j)
			    continue;
#ifdef MD_SPHERICAL_WALL
			  if (j==sphWall || j==sphWallOuter)
			    continue;
#endif
			  check_shift(i, j, shift);
			  assign_bond_mapping(i, j);
#ifdef EDHE_FLEX
			  if (nbondsFlex==0)
			    continue;
#endif 
			  warn = check_bonds_ij(i, j, shift);  
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
    }
}
#endif
#ifdef EDHE_FLEX
int check_bonds_ij(int i, int j, double shift[3])
{
  int nn, warn, amin, bmin, nbonds;
  double dist;
  dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
  nbonds = nbondsFlex;
  warn = 0;
#ifdef MD_GHOST_IGG
  if (Oparams.ghostsim)
    {
     if (areGhost(i,j))
       return warn;
    }
#endif

  for (nn=0; nn < nbonds; nn++)
    {
      if (dists[nn]<0.0 && fabs(dists[nn])>OprogStatus.epsd 
	  && !bound(i,j,mapbondsa[nn], mapbondsb[nn]) )
	{
	  warn=1;
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
  return warn;
}
#endif
#if defined(EDHE_FLEX) 
void check_all_bonds_NLL(void)
{
  int i, k, j, warn;
  double shift[3]={0.0,0.0,0.0};
  warn = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      //if (warn)
	//break;
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
      for (k=0; k < nebrTab[i].len; k++)
	{
	  j = nebrTab[i].list[k];
	  if (i == j)
	    continue;
#ifdef MD_SPHERICAL_WALL
	  if (j==sphWall || j==sphWallOuter)
	    continue;
#endif
	  /* calculate shift */
#ifdef MD_LXYZ
	  shift[0] = L[0]*rint((rx[i]-rx[j])/L[0]);
	  shift[1] = L[1]*rint((ry[i]-ry[j])/L[1]);
#ifdef MD_EDHEFLEX_WALL
	  if (!OprogStatus.hardwall)
	    shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
	  else
	    shift[2] = 0.0;
#else
	  shift[2] = L[2]*rint((rz[i]-rz[j])/L[2]);
#endif
#else
	  shift[0] = L*rint((rx[i]-rx[j])/L);
	  shift[1] = L*rint((ry[i]-ry[j])/L);
#ifdef MD_EDHEFLEX_WALL
	  if (!OprogStatus.hardwall)
	    shift[2] = L*rint((rz[i]-rz[j])/L);
	  else
	    shift[2] = 0.0;
#else
	  shift[2] = L*rint((rz[i]-rz[j])/L);
#endif
#endif
	  //check_shift(i, j, shift);
	  assign_bond_mapping(i,j);
	  if (nbondsFlex==0)
	    continue;
	  warn = check_bonds_ij(i,j,shift); 
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
    }
}

extern void add_bond(int na, int n, int a, int b);

#ifdef MD_MULTIPLE_LL
void find_bonds_one_MLL(int i)
{
  int nn,  amin, bmin, j, nbonds;
  double shift[3], dist;
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k;
  int *cellListL, *inCellL[3], cellsx, cellsy, cellsz;
  int nc, t1=-1, t2=-1, nl, numll;

  numll = Oparams.ntypes*(Oparams.ntypes+1)/2;
  for (nl = 0; nl < numll; nl++)
    {
      //printf("find bonds one i=%d nl=%d ignore=%d numll=%d\n",i,  nl, ignoreMLL[nl], numll);
      if (ignoreMLL[nl])
	continue;
      get_types_from_nl(nl, &t1, &t2);
      if (t1==typeOfPart[i])
	nc = t2;
      else if (t2==typeOfPart[i])
	nc = t1;
      else
	continue;
      for (k=0; k < 3; k++)
	inCellL[k] = inCellMLL[nc][k];
      cellListL = cellListMLL[nl];
      cellsx = cellsxMLL[nl];
      cellsy = cellsyMLL[nl];
      cellsz = cellszMLL[nl];
      for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_EDHEFLEX_WALL
      if (OprogStatus.hardwall==1)
	{
	  if (inCellL[2][i] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
	  if (inCellL[2][i] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
	}
#endif 
      for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
	{
	  jZ = inCellL[2][i] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      //printf("BOHHHH\n");
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
	      jY = inCellL[1][i] + iY;    
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
		  jX = inCellL[0][i] + iX;    
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
		  for (j = cellListL[j]; j > -1; j = cellListL[j]) 
		    {
		      if (i == j)
			continue;
#ifdef MD_SPHERICAL_WALL
		      if (j==sphWall || j==sphWallOuter)
			continue;
#endif
		      check_shift(i, j, shift);
		      assign_bond_mapping(i,j);
		      dist = calcDistNegSP(Oparams.time, 0.0, i, j, shift, &amin, &bmin, dists, -1);
		      nbonds = nbondsFlex;
		      //printf("nbondsFlex=%d checking i=%d j=%d\n", nbondsFlex, i, j);
		      for (nn=0; nn < nbonds; nn++)
			{
			  if (dists[nn]<0.0 && !bound(i, j, mapbondsaFlex[nn], mapbondsbFlex[nn]))
			    {
			      //printf("[find_bonds_one] found bond between ghost particles! i=%d j=%d typei=%d typej=%d\n",
			      //	 i, j, typeOfPart[i], typeOfPart[j]);
			      add_bond(i, j, mapbondsaFlex[nn], mapbondsbFlex[nn]);
			      add_bond(j, i, mapbondsbFlex[nn], mapbondsaFlex[nn]);
			    }
			}
		    }
		}
	    }
	}
    }
}
#endif
#endif
