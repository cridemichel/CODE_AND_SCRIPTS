#include<mdsimul.h>

void set_cells_size(void)
{

}
ProcessCollision()
{}

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
	  if (nc > 0 && inCellMLL[nc][2][na]==cellszMLL[nl]-1)
	    ignorecross[2] = 1;
	  else
	    ignorecross[2] = 0;
	}
      else 
	{
	  signDir[2] = 1;/* direzione negativa */
	  if (nc > 0 && inCellMLL[nc][2][na]==0)
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
	  if (nc > 0 && inCellMLL[nc][0][na]==cellsxMLL[nl]-1)
	    ignorecross[0] = 1;
	  else
	    ignorecross[0] = 0;
	  signDir[0] = 0;/* direzione positiva */
	}
      else
	{
	  if (nc > 0 && inCellMLL[nc][0][na]==0)
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
	  if (nc > 0 && inCellMLL[nc][1][na]==cellsyMLL[nl]-1)
	    ignorecross[1] = 1;
	  else
	    ignorecross[1] = 0;
	  signDir[1] = 0;
	}
      else
	{
	  if (nc > 0 && inCellMLL[nc][1][na]==0)
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
#if 1
  if (tm[k]<0.0)
    {
      printf("tm[%d]: %.15G\n", k, tm[k]);
      tm[k] = 0.0;
      printf("real cells: %d %d %d\n", (int)((rx[na] + L2) * cellsxMLL[nl] / L),
	     (int)((ry[na] + L2) * cellsyMLL[nl] / L), (int)((rz[na] + L2)  * cellszMLL[nl] / L));
#if 0
      printf("idA=%d idB=%d treeQIndex[%d]=%d treeStatus[]=%d\n ", treeIdA[na+1], treeIdB[na+1], na+1, treeQIndex[na+1], treeStatus[na+1]);
      printf("currentIndex=%d\n", OprogStatus.curIndex);
      printf("nc=%d na=%d nl=%d\n",nc,na,nl);
      printf("tm[%d]<0 step %lld na=%d\n", k, (long long int)Oparams.curStep, na);
      printf("Cells(%d,%d,%d)\n", inCell[nc][0][na], inCell[nc][1][na], 
	     inCell[nc][2][na]);
      printf("cells= (%d,%d,%d)\n", cellsx[nl], cellsy[nl], cellsz[nl]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
      exit(-1);
      /*tm[k] = 0.0;*/
#endif
    }
#endif
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
#if 0
      printf("<<< NOT IGNORE >>> evIdA=%d nc=%d time=%.15G k=%d\n", na, nc, Oparams.time+tm[k],k);
      printf("inCell = %d %d %d <=>\n", inCell[nc][k][na], inCell[nc][k][na], inCell[nc][k][na]);
#endif
      ScheduleEventBarr (na, ATOM_LIMIT + evCode, nc, 0, MD_EVENT_NONE, Oparams.time + tm[k]);
      cctime = Oparams.time + tm[kk];
    }
  //printf("===>crossevtodel[%d]:%d\n", na, crossevtodel[na]);
  //printf("schedule event [WallCrossing](%d,%d) tm[%d]: %.16G time=%.15G evCode:%d\n", 
//	 na, ATOM_LIMIT+evCode, k, tm[k], tm[k]+Oparams.time, evCode);
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


PredictEventMulLL()
{}


PredictEventMulLL_NNL()
{}


