/* tree.c */
  
/********************************************************************
   This software is copyrighted material, reproduced from the book 
   "The Art of Molecular Dynamics Simulation" by D. C. Rapaport, 
   published by Cambridge University Press (1995).
*********************************************************************/
#include <mdsimul.h>
extern int **tree, evIdA, evIdB;
extern double *treeTime;
void DeleteEvent(int );
void NextEvent(void);
extern int poolSize;
#if 0
#define treeLeft   tree[0]
#define treeRight  tree[1]
#define treeUp     tree[2]
#define treeCircAL tree[3]
#define treeCircBL tree[4]
#define treeCircAR tree[5]
#define treeCircBR tree[6]
#define treeIdA    tree[7]
#define treeIdB    tree[8]
#endif
/* NOTA: treeIdA[0] punta al primo nodo della lista dei nodi liberi nel pool 
 */
#define check(X_,XS_) /*if (X_<0) {printf(XS_ "minore di 0!\n"); exit(-1);}*/
#define check2(X_,XS_) /*if (X_==0) {printf(XS_ "minore di 0!\n"); exit(-1);}*/
void ErrExit(char *str)
{
  printf(str);
  exit(-1);
}
int check_node(char* str, int id, int idNew, int idUp)
{
  /*int idd;*/
#if 0
  if (Oparams.curStep==1266)
    {
      printf("(id:%d,u:%dr:%d,l:%d) ",
	     id, treeUp[id], treeRight[id], treeLeft[id]);
    }
#endif
  if (treeUp[id] != idUp && idUp != 0)
    {
      printf("[%s] ERRORE\n", str);
      printf("treeUp: %d treeRight: %d treeLeft: %d\n", treeUp[id], treeRight[id],
	     treeLeft[id]);
	    printf("right[up]: %d left[up]: %d\n", treeRight[treeUp[id]],
		   treeLeft[treeUp[id]]);
      printf("STEP: %d idUp: %d id:%d\n", Oparams.curStep, idUp, id);
      return 1;
    }
   
  if (treeRight[id] != -1)
    {
      if (check_node(str, treeRight[id], idNew, id))
	return 1;
    }
  if (treeLeft[id] != -1 && id != 0)
    {
      if (check_node(str, treeLeft[id], idNew, id))
	return 1;
    }
  return 0;
}
void ScheduleEvent (int idA, int idB, double tEvent) 
{
  int id, idNew, more;
  id = 0;

  MD_DEBUG2(printf("#%d ScheduleEvent() idA:%d idB:%d evtime:%.15f\n", 
		   Oparams.curStep, idA, idB,
		  tEvent));
  if ((idB < ATOM_LIMIT ||
      idB >= ATOM_LIMIT + 2 * NDIM) && idB < ATOM_LIMIT + 100)
    {
      /* urto con altra particella o altro evento (misura o output)*/
      if (treeIdA[0] == -1)
	ErrExit ("empty event pool");
      /* treeIdA[0] è un puntatore ad un nodo utilizzabile nel pool */
      idNew = treeIdA[0];
      MD_DEBUG2(
      if (idB >= ATOM_LIMIT + 2 * NDIM)
	printf("idNew: %d tEvent: %.15f\n", idNew, tEvent));
      /* all'inizio treeCircAR[treeIdA[0]] = treeIdA[0]+1 quindi è un nodo 
       * non utilizzato nel pool */
      treeIdA[0] = treeCircAR[treeIdA[0]];
    }
  else 
    idNew = idA + 1;
  /* Se qui vuol dire che si tratta di un cell-crossing o 
     di un urto con parete
     NOTA: urto con parete e cell-crossing sono esclusivi, per cui basta un nodo 
     inoltre c'è sempre un evento di tale tipo associato con ogni particella 
     */
 check(idNew, "idNew");
  /* treeRight[id] == -1 => il calendario è vuoto */
  if (treeRight[id] == -1) 
    treeRight[id] = idNew;
  else 
    {
      /* Cerca la giusta collocazione nell'albero per l'evento da
       * schedulare */
      more = 1; 
      id = treeRight[id];
      while (more) 
	{
	  if (tEvent <= treeTime[id]) 
	    {
	      if (treeLeft[id] > -1) 
		id = treeLeft[id];
	      else 
		{
		  more = 0;    
		  treeLeft[id] = idNew;
		}
	    } 
	  else
	    {
	      if (treeRight[id] > -1) 
		id = treeRight[id];
	      else 
		{
		  more = 0;    
		  treeRight[id] = idNew;
		} 
	    }
	} 
    }
    
  if (idB < ATOM_LIMIT) 
    {
      /* Chiaramente ad idNew sono associate le particelle idA e idB
       * relative all'evento che si sta schedulando */
      /* inserisce idNew nella circular list della particelle idA */
      treeCircAR[idNew] = treeCircAR[idA + 1];
      treeCircAL[idNew] = idA + 1;
      treeCircAL[treeCircAR[idA + 1]] = idNew;
      treeCircAR[idA + 1] = idNew;
      /* inserisce idNew nella circular list di idB */
      treeCircBR[idNew] = treeCircBR[idB + 1];
      treeCircBL[idNew] = idB + 1;
      treeCircBL[treeCircBR[idB + 1]] = idNew;
      treeCircBR[idB + 1] = idNew;
    }
  treeTime[idNew] = tEvent;
  treeIdA[idNew] = idA;    
  treeIdB[idNew] = idB;
  treeLeft[idNew] = treeRight[idNew] = -1;
  treeUp[idNew] = id;
}

void NextEvent (void) 
{
  int id, idAx, idBx, idd, idNow, idtx;
  /* Il nodo root (0), a cui non è associato alcun evento,
   * è linkato con il suo right pointer al primo nodo che contiene
   * un evento */
  idNow = treeRight[0];  
  /* Cerca l'evento con tempo minore 
   * NOTA: l'albero è ordinato e ogni nodo sinistro ha un tempo inferiore */
  while (treeLeft[idNow] > -1) 
    idNow = treeLeft[idNow];
  OprogStatus.time = treeTime[idNow];   
  evIdA = treeIdA[idNow];    
  evIdB = treeIdB[idNow];
  MD_DEBUG2(printf("[ NextEvent ] #%d event(%d,%d) curtime:%f\n", 
		   Oparams.curStep, evIdA, evIdB, OprogStatus.time));
  if (evIdB < ATOM_LIMIT + 2 * NDIM) 
    /* Se si tratta di un urto fra particelle o con le pareti...
     * qui in sostanza si considerano solo gli eventi che cambiano lo stato 
     * della particella, cioè la sua velocità */
    {
      /* qui incrementa di 1 poiché il nodo root non contiene eventi */
      if (evIdA < evIdB) 
	{
	  idAx = evIdA + 1;    
	  idBx = evIdB + 1;
	} 
      else
	{
  	  idAx = evIdB + 1;    
  	  idBx = evIdA + 1;
    	}
      idtx = idBx - idAx;
      if (evIdB >= ATOM_LIMIT) 
	idBx = idAx;
      /* qui rimuove dal calendario tutti gli eventi in cui sono coinvolte A o B */
      for (id = idAx; id <= idBx; id += idtx) 
	{
	  DeleteEvent (id);
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
	  check(treeCircAL[id], "treeCircAL");
	  treeCircAR[treeCircAL[id]] = treeIdA[0];
	  treeIdA[0] = treeCircAR[id];
	  /* tutte le liste circolari vengono svuotate */
	  treeCircAL[id] = treeCircAR[id] = id;
	  for (idd = treeCircBL[id]; idd != id; idd = treeCircBL[idd]) 
	    {
	      /* vedere sopra infatti è lo stesso solo per la lista in cui
	       * la particella è la prima della coppia (A) */
	      check(treeCircAL[idd], "treeCircAL");
	      check(treeCircAR[idd], "treeCircAR");
	      treeCircAR[treeCircAL[idd]] = treeCircAR[idd];
	      treeCircAL[treeCircAR[idd]] = treeCircAL[idd];
	      DeleteEvent (idd);
	      treeCircAR[idd] = treeIdA[0];    
	      treeIdA[0] = idd;
	    }
	  treeCircBL[id] = treeCircBR[id] = id;
	}
    } 
  else 
    {
      /* L'evento è un cell-crossing o un evento generico (output o misura) */
      DeleteEvent (idNow);
      if (evIdB < ATOM_LIMIT + 100) 
	{
	  treeCircAR[idNow] = treeIdA[0];
	  treeIdA[0] = idNow;
	} 
    }
  /*check_node("next event", 0, -1, 0);*/
  if (OprogStatus.time < 0)  
    { 
       printf("#%d time = %.15f\n", Oparams.curStep, OprogStatus.time); 
       exit(-1);     
     } 
    
  /*printf("Next event evIdA: %d, evIdB:%d\n", evIdA, evIdB);*/
}

void DeleteEvent (int id)
{
  int idp, idq, idr;

  MD_DEBUG2(printf("[ DeleteEvent ] deleting node #%d\n", id));
  idr = treeRight[id];
  if (idr == -1)
    idq = treeLeft[id];
  else 
    {
      if (treeLeft[id] == -1) 
	idq = idr;
      else 
	{
	  if (treeLeft[idr] == -1) 
	    idq = idr;
	  else 
	    {
	      idq = treeLeft[idr];
	      while (treeLeft[idq] > -1)
		{
		  idr = idq;    
		  idq = treeLeft[idr];
		}
	      treeLeft[idr] = treeRight[idq];
	      check2(treeRight[idq], "treeRight[idq]");
	      check2(treeRight[id], "treeRight[id]");
	      check(idr, "idr");
	      check(idq, "idq");
	      check(id, "id");
	      if (treeRight[idq]>-1)
	        treeUp[treeRight[idq]] = idr;
	      treeRight[idq] = treeRight[id];
	      if (treeRight[id]>-1)
	        treeUp[treeRight[id]] = idq;
	    }
	  treeUp[treeLeft[id]] = idq;
	  treeLeft[idq] = treeLeft[id];
	} 
    }
  idp = treeUp[id];    
  if (idq > -1)
    treeUp[idq] = idp;
  if (treeRight[idp] != id)
    treeLeft[idp] = idq;
  else 
    treeRight[idp] = idq;
  MD_DEBUG2(printf("idq: %d idp: %d id:%d treeUp[idp]:%d Left[id]: %d Right[id]:%d LUp[idp]:%d RUp[idp]:%d\n", 
		   idq, idp, id,
		   treeUp[idp], treeLeft[id], treeRight[id],
		   treeLeft[treeUp[idp]], treeRight[treeUp[idp]]));
}

void InitEventList (void) 
{
  int id;
  treeLeft[0] = treeRight[0] = -1;
  treeIdA[0] = Oparams.parnum + 1;
  /* i nodi da Oparams.parnum + 1 (compreso) in poi sono il pool (cioè quelli dinamici) */
  for (id = treeIdA[0]; id <= poolSize - 2; id++) 
    treeCircAR[id] = id + 1;
  treeCircAR[poolSize-1] = -1;
  for (id = 1; id <= Oparams.parnum; id++) 
    {
      treeCircAL[id] = treeCircBL[id] = id;
      treeCircAR[id] = treeCircBR[id] = id;
    }
}
