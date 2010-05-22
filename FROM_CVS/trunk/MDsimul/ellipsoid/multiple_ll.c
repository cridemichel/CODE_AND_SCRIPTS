#include<mdsimul.h>
int is_superellips_type(int pt)
{
  if (typesArr[pt].n[0]==2.0 && typesArr[pt].n[1]==2.0 &&
      typesArr[pt].n[2]==2.0)
    return 0;
  else 
    return 1;
}

inline int get_linked_list(int na, int nc)
{
  int typena, sum, t1, nc1;
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
  typena = typeOfPart[na];
  if (nc==typeOfPart[na])
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
int is_a_sphere_NNL_type(int pt)
{
  int i, k1, k2, is_sph;
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
void get_types_from_nl(int nl, int *t1, int *t2)
{
  int ta, tb;
  /* funzione inversa di get_linked_list */
  for (ta = 0; ta < Oparams.ntypes; ta++)
    for (tb = 0; tb < Oparams.ntypes; tb++)
      {
	if (ta >= tb)
	  continue;
	if (nl==get_linked_list(ta, tb))
	  {
	    *t1 = ta;
	    *t2 = tb;
	    return;
	  }
      }
}
double calc_rcut_type(int t)
{
  double rcutA;
  int kk;
  double ax[3], del;
  int i;
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
  return OprogStatus.rcutfactMLL*rcutA;
}

double calc_rcut(int nl)
{
  int t1, t2;
  double rc, rc1, rc2;
  /* le celle liste per ora vengono sempre scelte automaticamente */
  get_types_from_nl(nl, &t1, &t2);
  rc1 = calc_rcut_type(t1);
  rc2 = calc_rcut_type(t2);
  rc = (rc1 + rc2)*0.5;
  return rc;
}
void set_cells_size(void)
{
  int nl, numll;
  double rcut;
  
  numll = Oparams.ntype*(Oparams.ntypes+1)/2;
  for (nl = 0; nl < numll; nl++)
    {
      //if (Oparams.rcut[nl] <= 0.0)
	//Oparams.rcut[nl] = pow(L*L*L / Oparams.parnum, 1.0/3.0); 
      rcut = calc_rcut(nl);
#ifdef MD_LXYZ
      cellsx[nl] = L[0] / rcut;
      cellsy[nl] = L[1] / rcut;
      cellsz[nl] = L[2] / rcut;
#else
      cellsx[nl] = L / rcut;
      cellsy[nl] = L / rcut;
      cellsz[nl] = L / rcut;
#endif
#ifdef MD_LXYZ
      printf("[%d] L=%.15G %.15G %.15G Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", nl, L[0], L[1], L[2],
	     rcut, cellsx[nl], cellsy[nl], cellsz[nl]);
#else
      printf("[%d] L=%.15G Oparams.rcut: %f cellsx:%d cellsy: %d cellsz:%d\n", nl, L,
	     rcut, cellsx[nl], cellsy[nl], cellsz[nl]);
#endif
    }

}
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
void ProcessCellCrossingMLL(void)
{
  int kk, n, k, nl;
  int nc, boxwall, nlcoll, nlcross, nc_bw, nlcross_bw, nlcoll_bw;
  int typei;

  UpdateAtom(evIdA);
  /* kk ci da la direzione lungo cui si sta realizzando il cell crossing */
  kk = evIdB - 100 - ATOM_LIMIT; 
  /* evIdC è semplicemente nc cioè ci dice il tipo di cella attraversata dalla particella evIdA, ossia
     nc = cella per interazine con tipo nc */
  nc = evIdC;

  typei = typeOfPart[evIdA];

  nl = get_linked_list(typei, nc);
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
  boxwall = check_boxwall(kk, nc, nlcross);
  /* questa condizione non si dovrebbe *mai* verificare, quindi 
   * in linea di principio le due linee che seguono potrebbero anche essere eliminate */
  if (nc==1 && boxwall)
    {
      printf("nc=1 and boxwall!!!! <===!!!\n");
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
  while (cellList[nlcross][n] != evIdA) 
    n = cellList[nlcross][n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellList[nlcross][n] = cellList[nlcross][evIdA];
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }

  if (boxwall)
    {
      //printf("BOXWALL nc=%d nc2=%d nl=%d nl2=%d evIdA=%d time=%.15G\n", nc, nc2, nl, nl2, evIdA, Oparams.time);
      n = (inCell[nc_bw][2][evIdA] * cellsy[nlcross_bw] + inCell[nc_bw][1][evIdA])*cellsx[nlcross_bw] + 
	inCell[nc_bw][0][evIdA]
	+ Oparams.parnum;
      while (cellList[nlcross_bw][n] != evIdA) 
	n = cellList[nlcross_bw][n];
      /* Eliminazione di evIdA dalla lista della cella n-esima della lista nl2 */
      cellList[nlcross_bw][n] = cellList[nlcross_bw][evIdA];
    }
  switch (kk)
    {
    case 0: 
      docellcross(0, vx[evIdA], &(rx[evIdA]), cellsx[nlcross], nc);
      break;
    case 1: 
      docellcross(1, vy[evIdA], &(ry[evIdA]), cellsy[nlcross], nc);
      break;
    case 2:
      docellcross(2, vz[evIdA], &(rz[evIdA]), cellsz[nlcross], nc);
      break;
    }
  PredictCellCross(evIdA, nc);
  PredictColl(evIdA, evIdB, nlcoll);
  n = (inCell[nc][2][evIdA] * cellsy[nlcross] + inCell[nc][1][evIdA])*cellsx[nlcross] + 
    inCell[nc][0][evIdA] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellList[nlcross][evIdA] = cellList[nlcross][n];
  cellList[nlcross][n] = evIdA;
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
#if 0
  printf("DOPO boxwall=%d nc=%d n=%d cellList[%d][%d]:%d\n",boxwall, nc, n, nlcross, n, cellList[nlcross][n]);
  printf("DOPO vel=(%f,%f,%f) inCell= %d %d %d\n", vx[evIdA], vy[evIdA], vz[evIdA], inCell[nc][0][evIdA],inCell[nc][1][evIdA], inCell[nc][2][evIdA]);
#endif
  if (boxwall)
    {
      switch (kk)
	{
	case 0: 
	  docellcross2(0, vx[evIdA], cellsx[nlcross_bw], nc_bw);
	  break;
	case 1: 
	  docellcross2(1, vy[evIdA], cellsy[nlcross_bw], nc_bw);
	  break;
      	case 2:
	  docellcross2(2, vz[evIdA], cellsz[nlcross_bw], nc_bw);
	  break;
	}
      if (crossevtodel[evIdA]!=-1)
	{
	  //printf("DELETING CROSS EVENT evIdA=%d\n", evIdA);
	  DeleteEvent(crossevtodel[evIdA]);
	  crossevtodel[evIdA] = -1;
	}
      PredictCellCross(evIdA, nc_bw);
      PredictColl(evIdA, evIdB, nlcoll_bw);
      n = (inCell[nc_bw][2][evIdA] * cellsy[nlcross_bw] + inCell[nc_bw][1][evIdA])*cellsx[nlcross_bw] + 
	inCell[nc_bw][0][evIdA] + Oparams.parnum;
      /* Inserimento di evIdA nella nuova cella (head) */
      cellList[nlcross_bw][evIdA] = cellList[nlcross_bw][n];
      cellList[nlcross_bw][n] = evIdA;
    }
}

/* 21/05/2010: nc indica la linked lists da considerare ossia se nc=typeOfPart[na]=interazione con lo stesso tipo
   mentre valori diversi si riferiscono alle interazioni con gli altri tipi, 
   ad es. se ci sono 4 tipi la particella na è tipo 1 allora nc=1 sono le celle per l'interazione
   con lo stesso tipo mentre 0 = interazione con particelle di tipo 0, 2 interazione con particelle di tipo 2
   3 con quelle di tipo 3. */
double PredictCellCross(int na, int nc)
{
  int ignorecross[3], k, evCode, signDir[NDIM]={0,0,0}, iA, nl;
  double tm[3], cctime=timbig;

  ignorecross[0] = ignorecross[1] = ignorecross[2] = 1;

  nl =  get_linked_list(na, nc);
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
#ifdef MD_EDHEFLEX_WALL
      /* il muro attualmente è solo lungo l'asse z */
      if (OprogStatus.hardwall && ((signDir[2]==0 && inCellMLL[2][na]==cellszMLL[nl]-1) || (signDir[2]==1 && inCellMLL[2][na]==0)))
	ignorecross[2] = 1;
#endif
      if (ignorecross[2])
	tm[2] = timbig;
      else
#ifdef MD_LXYZ
	tm[2] = ((inCellMLL[nc][2][na] + 1 - signDir[2]) * L[2] /
  		 cellszMLL[nl] - rz[na] - L2[2]) / vz[na];
#else
	tm[2] = ((inCellMLL[nc][2][na] + 1 - signDir[2]) * L /
  		 cellszMLL[nl] - rz[na] - L2) / vz[na];
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

  if (vy[na] != 0.) 
    {
      if (vy[na] > 0.) 
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
      if (tm[0] <= tm[1]) k = 0;
      else k = 1;
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
      printf("real cells: %d %d %d\n", (int)((rx[na] + L2) * cellsxMLL[nl] / L),
	     (int)((ry[na] + L2) * cellsyMLL[nl] / L), (int)((rz[na] + L2)  * cellszMLL[nl] / L));
    }
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;// + 3*nc;
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 
  MD_DEBUG15(printf("schedule event [WallCrossing](%d,%d) tm[%d]: %.8G\n", 
		    na, ATOM_LIMIT+evCode, k, tm[k]));

  if (!ignorecross[k])
    {
      ScheduleEventBarr (na, ATOM_LIMIT + evCode, nc, 0, MD_EVENT_NONE, Oparams.time + tm[k]);
      cctime = Oparams.time + tm[kk];
    }
  return cctime;
}

void rebuildMultipleLL(void)
{
  int nl, numll, maxnc;
  
  numll = Oparams.ntype*(Oparams.ntypes+1)/2;
  for (nl=0; nl < numll; nl++)
    {
      for (j = 0; j < cellsxMLL[nl]*cellsyMLL[nl]*cellszMLL[nl] + Oparams.parnum; j++)
	cellListMLL[nl][j] = -1;
    }

  maxnc = Oparams.ntypes;
  for (nc=0; nc < maxnc; nc++)
    {
      /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
      for (n = 0; n < Oparams.parnum; n++)
	{
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
 
	  nl = get_linked_list(n, nc);
	  atomTime[n] = Oparams.time;
#ifdef MD_LXYZ
	  inCellMLL[nc][0][n] =  (rx[n] + L2) * cellsxMLL[nl] / L[0];
	  inCellMLL[nc][1][n] =  (ry[n] + L2) * cellsyMLL[nl] / L[1];
	  inCellMLL[nc][2][n] =  (rz[n] + L2)  * cellszMLL[nl] / L[2];
#else
	  inCellMLL[nc][0][n] =  (rx[n] + L2) * cellsxMLL[nl] / L;
	  inCellMLL[nc][1][n] =  (ry[n] + L2) * cellsyMLL[nl] / L;
	  inCellMLL[nc][2][n] =  (rz[n] + L2)  * cellszMLL[nl] / L;
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
	  cellListMLL[nl][n] = cellListMLL[nl][j];
	  cellListMLL[nl][j] = n;
	}
    }
}


PredictEventMLL()
{

}

BuildNNLwithMLL()
{


}

