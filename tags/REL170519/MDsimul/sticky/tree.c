/* tree.c */
  
/********************************************************************
   This software is copyrighted material, reproduced from the book 
   "The Art of Molecular Dynamics Simulation" by D. C. Rapaport, 
   published by Cambridge University Press (1995).
*********************************************************************/
#include <mdsimul.h>
extern int **tree, evIdA, evIdB, evIdC, evIdD, evIdE;
extern double *treeTime, *treeRxC, *treeRyC, *treeRzC;
void DeleteEvent(int );
void NextEvent(void);
extern int poolSize;
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
extern int *crossevtodel;
#endif
extern double rxC, ryC, rzC;
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
#ifdef MD_CALENDAR_HYBRID
int numevPQ=0; /* numero di eventi nella PQ (i.e. binary tree) */
int totevHQ=0;
int overevHQ=0;
//int linearLists[nlists+1];/*+1 for overflow*/ /* dynamically allocated */
int *linearLists; /* dynamically allocated */
//int currentIndex=0;
int populatePQ(void);
//double baseIndex=0;//in OprogStatus
void deleteFromEventQ(int e);
int insertInEventQ(int p);
#endif
void ErrExit(char *str)
{
  printf(str);
  exit(-1);
}
int check_node(char* str, int id, int idNew, int idUp)
{
  int idd;
  if (treeUp[id] != idUp && idUp != 0)
    {
      printf("[%s] ERRORE\n", str);
      printf("treeUp: %d treeRight: %d treeLeft: %d\n", treeUp[id], treeRight[id],
	     treeLeft[id]);
	    printf("right[up]: %d left[up]: %d\n", treeRight[treeUp[id]],
		   treeLeft[treeUp[id]]);
      printf("STEP: %lld idUp: %d id:%d\n",(long long int) Oparams.curStep, idUp, id);
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

#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
#ifdef MD_CALENDAR_HYBRID
int get_new_node(int idA, int idB, int idata)
{
  /* get_new_node(): questa funzione ottiene un nodo libero */
  int idNew;
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
	printf("idNew: %d tEvent: %.15f\n", idNew, treeTime[idNew]));
      /* all'inizio treeCircAR[treeIdA[0]] = treeIdA[0]+1 quindi è un nodo 
       * non utilizzato nel pool */
      treeIdA[0] = treeCircAR[treeIdA[0]];
      /* essendo qui si tratta di un urto quindi bisobna valutare se inserire
	 tale evento nella BSPQ (i.e. Hybrid Queue o nella PQ) */
    }
  else 
    { /* NOTA: siccome nel caso della silica 
	 ogni molecola appartiene a due linked list 
	 non basta un solo slot ma ne servono due qui! */
      /* Se qui vuol dire che si tratta di un cell-crossing o 
	 di un urto con parete
         NOTA: urto con parete e cell-crossing sono esclusivi, per cui basta un nodo 
         inoltre c'è sempre un evento di tale tipo associato con ogni particella 
	 */
       	idNew = 1 + (Oparams.parnum*idata + idA);
	if (idata)
	  {
	    crossevtodel[idA] = idNew;  
	  }
    }
  return idNew;

}
void setPQnode(int idNew, int idA, int idB, int idata, int idatb, int idcollcode, double tEvent)
{
  /* assegna i dati relativi all'evento al nodo idNew */
  treeTime[idNew] = tEvent;
  treeIdA[idNew] = idA;    
  treeIdB[idNew] = idB;
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
}
void InsertPQ(int idNew) 
{
  int more, id;
  double tEvent;
  numevPQ++;
  id = 0;
  tEvent = treeTime[idNew];
#ifdef MD_CALENDAR_HYBRID
  treeStatus[idNew] = 2; /* 2 = belonging to binary tree */
#endif
  MD_DEBUG2(printf("InsertPQ tEvent=%.15G idA=%d idB=%d\n", tEvent, treeIdA[idNew], treeIdB[idNew]);)
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
  treeUp[idNew] = id;
  treeLeft[idNew] = treeRight[idNew] = -1;
}
void insertInCircularLists(int idNew)
{
  int idA, idB;
  idA = treeIdA[idNew];
  idB = treeIdB[idNew];
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
}
void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent) 
{
  int idNew;
  idNew = get_new_node(idA, idB, idata);
  /* assegna i dati al nodo relativi all'evento */
  setPQnode(idNew, idA, idB, idata, idatb, idcollcode, tEvent);
  insertInEventQ(idNew);
  /* 07/05/2010: le liste circolari servono per eliminare tutti gli urti in cui è coinvolta
     una certa particella.
     Notare anche se l'evento è stato inserito nelle liste lineari va comunque inserito
     nelle liste circolari per poter essere rimosso dopo un urto.
   */
  insertInCircularLists(idNew);
}
#else
void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent) 
{
  int id, idNew, more;
  id = 0;

  MD_DEBUG2(printf("#%lld ScheduleEvent() idA:%d idB:%d evtime:%.15f\n", 
		   (long long int)Oparams.curStep, idA, idB,
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
      /* essendo qui si tratta di un urto quindi bisobna valutare se inserire
	 tale evento nella BSPQ (i.e. Hybrid Queue o nella PQ) */
    }
  else 
    { /* NOTA: siccome nel caso della silica 
	 ogni molecola appartiene a due linked list 
	 non basta un solo slot ma ne servono due qui! */
      /* Se qui vuol dire che si tratta di un cell-crossing o 
	 di un urto con parete
         NOTA: urto con parete e cell-crossing sono esclusivi, per cui basta un nodo 
         inoltre c'è sempre un evento di tale tipo associato con ogni particella 
	 */
       	idNew = 1 + (Oparams.parnum*idata + idA);
	if (idata)
	  {
	    crossevtodel[idA] = idNew;  
	  }
    }
  
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
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
  treeLeft[idNew] = treeRight[idNew] = -1;
  treeUp[idNew] = id;
}
#endif
#else
#ifdef MD_CALENDAR_HYBRID
int get_new_node(int idA, int idB, int idata)
{
  /* get_new_node(): questa funzione ottiene un nodo libero */
  int idNew;
  MD_DEBUG2(printf("#%lld ScheduleEvent() idA:%d idB:%d evtime:%.15f\n", 
		   (long long int)Oparams.curStep, idA, idB,
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
  return idNew;
}
void setPQnode(int idNew, int idA, int idB, int idata, int idatb, int idcollcode, double tEvent)
{
  /* assegna i dati relativi all'evento al nodo idNew */
  treeTime[idNew] = tEvent;
  treeIdA[idNew] = idA;    
  treeIdB[idNew] = idB;
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
}
void InsertPQ(int idNew) 
{
  int more, id;
  double tEvent;

  numevPQ++;
  id = 0;
  tEvent = treeTime[idNew];

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
  treeUp[idNew] = id;
  treeLeft[idNew] = treeRight[idNew] = -1;
} 
void insertInCircularLists(int idNew)
{
  int idA, idB;
  idA = treeIdA[idNew];
  idB = treeIdB[idNew];
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
} 
void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent) 
{
  int idNew;
  idNew = get_new_node(idA, idB, idata);
  /* assegna i dati al nodo relativi all'evento */
  setPQnode(idNew, idA, idB, idata, idatb, idcollcode, tEvent);
  insertInEventQ(idNew);
  /* 07/05/2010: le liste circolari servono per eliminare tutti gli urti in cui è coinvolta
     una certa particella.
     Notare anche se l'evento è stato inserito nelle liste lineari va comunque inserito
     nelle liste circolari per poter essere rimosso dopo un urto.
   */
  insertInCircularLists(idNew);
}
#else
void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent) 
{
  int id, idNew, more;

  id = 0;

  MD_DEBUG2(printf("#%lld ScheduleEvent() idA:%d idB:%d evtime:%.15f\n", 
		   (long long int)Oparams.curStep, idA, idB,
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
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
  treeLeft[idNew] = treeRight[idNew] = -1;
  treeUp[idNew] = id;
}
#endif
#endif
void ScheduleEvent(int IdA, int IdB, double tEvent)
{
  ScheduleEventBarr(IdA, IdB, 0, 0, MD_EVENT_NONE, tEvent);
}
#if 0
void ScheduleEvent (int idA, int idB, double tEvent) 
{
  int id, idNew, more;
  id = 0;

  MD_DEBUG2(printf("#%lld ScheduleEvent() idA:%d idB:%d evtime:%.15f\n", 
		   (long long int)Oparams.curStep, idA, idB,
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
  treeRxC[idNew] = rxC;
  treeRyC[idNew] = ryC;
  treeRzC[idNew] = rzC;
  treeIdA[idNew] = idA;    
  treeIdB[idNew] = idB;
  treeLeft[idNew] = treeRight[idNew] = -1;
  treeUp[idNew] = id;
}
#endif
void delete_events(int evIdA)
{
  int id, idd;
  
  id = evIdA + 1;
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
void NextEvent (void) 
{
  int id, idAx, idBx, idd, idNow, idtx;
  /* Il nodo root (0), a cui non è associato alcun evento,
   * è linkato con il suo right pointer al primo nodo che contiene
   * un evento */
#ifdef MD_CALENDAR_HYBRID
  populatePQ();
#endif
  idNow = treeRight[0];  

  /* Cerca l'evento con tempo minore 
   * NOTA: l'albero è ordinato e ogni nodo sinistro ha un tempo inferiore */
  while (treeLeft[idNow] > -1) 
    idNow = treeLeft[idNow];
  Oparams.time = treeTime[idNow];   
#if 0
  if (treeUp[idNow] == 0 && treeRight[idNow]==-1 && treeLeft[idNow]==-1)
    {
      printf("solo un evento!\n");
      exit(-1);
    }
#endif
#if 0
  rxC = treeRxC[idNow];
  ryC = treeRyC[idNow];
  rzC = treeRzC[idNow];
#endif
  //printf("INIZIO evtime = %.15G\n", Oparams.time);
  evIdA = treeIdA[idNow];    
  evIdB = treeIdB[idNow];
  evIdC = treeIdC[idNow];
  evIdD = treeIdD[idNow];
  evIdE = treeIdE[idNow];
  /*
     printf("[ NextEvent ] #%lld event(%d,%d) curtime:%f\n", 
     (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time);
  */
#if 0
  if (treeIdA[idNow] == 20446 && treeIdB[idNow] == ATOM_LIMIT+102)
    {
      printf("[ NextEvent ] #%lld event(%d,%d) curtime:%f treeRight[0]=%d(idA=%d idB=%d)\n", (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time, treeRight[0], treeIdA[treeRight[0]], treeIdB[treeRight[0]]);
      printf("idNow=%d idA=%d idB=%d time=%.15G\n",idNow, treeIdA[idNow], treeIdB[idNow], treeTime[idNow]);
	printf("treeQindex[%d]=%d\n", treeRight[0], treeQIndex[treeRight[0]]);
    }
#endif
#if 0
  if (evIdA==24||evIdB==24)
    {
  printf("[ NextEvent ] #%lld event(%d,%d) curtime:%.15G treeRight[0]=%d(idA=%d idB=%d)\n", (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time, treeRight[0], treeIdA[treeRight[0]], treeIdB[treeRight[0]]);
      printf("idNow=%d idA=%d idB=%d time=%.15G\n",idNow, treeIdA[idNow], treeIdB[idNow], treeTime[idNow]);
	printf("treeQindex[%d]=%d\n", treeRight[0], treeQIndex[treeRight[0]]);
    }

#endif
  MD_DEBUG2(printf("[ NextEvent ] #%lld event(%d,%d) curtime:%f\n", 
		   (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time));
  MD_DEBUG(printf("[ NextEvent ] #%lld event(%d,%d) curtime:%f\n", 
		   (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time));
  if (evIdB < ATOM_LIMIT + 2 * NDIM) 
    /* Se si tratta di un urto fra particelle o con le pareti...
     * qui in sostanza si considerano solo gli eventi che cambiano lo stato 
     * della particella, cioè la sua velocità */
    {
      //printf("QUI evtime = %.15G\n", Oparams.time);
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
	  /* qui elimina anche gli eventi relativi al cell crossing con le altre liste.
	   * Notare che se si hanno più di due specie si devono eliminare *tutti* gli eventi 
	   * di cell-crossing magari con un loop */
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
	  if (crossevtodel[id-1]!=-1)
	    {
#ifdef MD_GRAVITY
#if 0
	      if (evIdA==24)
		printf("part. N. %d time=%.15G idA=%d idB=%d deleting crossevtodel\n", id-1,Oparams.time, evIdA, evIdB);
#endif
#endif
	      DeleteEvent (id+Oparams.parnum);
	    }
	  /* notare che la condizione crossevtodel[id-1]!=-1 è necessaria in quanto 
	   * in quanto se una particella attraversa le pareti del box l'evento con nc==1
	   * non viene schedulato visto che ProcessCellCrossing si occupa di entrambi gli eventi
	   * (nc==0 e nc==1) */
	  /* Serve assolutamente porre crossevtodel a -1  per evitare che ProcessCellCross 
	   * tenti di rimuovere un evento già rimosso. */
	  crossevtodel[id-1] = -1;
#endif
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
      
    } 
  else 
    {
      /* L'evento è un cell-crossing o un evento generico (output o misura) */
      DeleteEvent (idNow);
      /* se evIdC == 1 vuol dire che si tratta di un cell crossing ma nc=1 */
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
      if (evIdB < ATOM_LIMIT + 100) 
	{
	  treeCircAR[idNow] = treeIdA[0];
	  treeIdA[0] = idNow;
	} 
      if (evIdA >=0 && evIdA < ATOM_LIMIT && evIdC==1)
	crossevtodel[evIdA] = -1;
#else
      if (evIdB < ATOM_LIMIT + 100) 
	{
	  treeCircAR[idNow] = treeIdA[0];
	  treeIdA[0] = idNow;
	}
#endif
    }

  /*check_node("next event", 0, -1, 0);*/
  if (Oparams.time < 0)  
    { 
       printf("#%lld time = %.15f\n",(long long int) Oparams.curStep, Oparams.time); 
       exit(-1);     
     } 
    
  /*printf("Next event evIdA: %d, evIdB:%d\n", evIdA, evIdB);*/
}
void DeletePQ (int id)
{
  int idp, idq, idr;

#if MD_CALENDAR_HYBRID
  numevPQ--;
#endif
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
#ifdef MD_BIG_DT
  treeUp[id] = -1;
#endif
  MD_DEBUG2(printf("idq: %d idp: %d id:%d treeUp[idp]:%d Left[id]: %d Right[id]:%d LUp[idp]:%d RUp[idp]:%d\n", 
		   idq, idp, id,
		   treeUp[idp], treeLeft[id], treeRight[id],
		   treeLeft[treeUp[idp]], treeRight[treeUp[idp]]));
}
#ifdef MD_CALENDAR_HYBRID
void DeleteEvent(int id)
{
  deleteFromEventQ(id);
}
#else
void DeleteEvent(int id)
{
  DeletePQ(id);
}
#endif
#ifdef MD_CALENDAR_HYBRID
void initHQlist(void)
{
  int i;
  /* base index non deve essere azzerato se non all'inizio della simulazione quando
     t = 0 o quando si fa un bigDt (08/05/10: Anche se di ciò ancora non sono sicurao al 100%) */
  //OprogStatus.baseIndex = 0;
  //currentIndex = 0;
  //OprogStatus.curIndex = 0;
  /* inizializzare anche le linked lists lineari? */
  for (i=0; i < OprogStatus.nlistsHQ+1; i++)
    linearLists[i] = -1;
  for (i=1; i < Oparams.parnum*OprogStatus.eventMult; i++) 
    treeStatus[i] = 0;/* 0 = free node */

  numevPQ = overevHQ = totevHQ = 0; 
  OprogStatus.curIndex=0;
  /* NOTA 11/10/2010: l'importante è che: baseIndex < OprogStatus.scaleHQ*Oparams.time */  
  OprogStatus.baseIndex = Oparams.time*OprogStatus.scaleHQ;

}
#endif
#if defined(MD_SILICA) && !defined(MD_USE_SINGLE_LL)
void InitEventList (void) 
{
  int id;
  treeLeft[0] = treeRight[0] = -1;
  treeIdA[0] = 2*Oparams.parnum + 1;
  /* i nodi da 2*Oparams.parnum + 1 (compreso) in poi sono il pool (cioè quelli dinamici) */
  for (id = treeIdA[0]; id <= poolSize - 2; id++) 
    treeCircAR[id] = id + 1;
#ifdef MD_BIG_DT
  for (id = 1; id < poolSize; id++) 
    treeUp[id] = -1;
#endif
  treeCircAR[poolSize-1] = -1;
  for (id = 1; id <= Oparams.parnum*2; id++) 
    {
      treeCircAL[id] = treeCircBL[id] = id;
      treeCircAR[id] = treeCircBR[id] = id;
#if 0
      treeCircAR[id + Oparams.parnum] = treeCircAR[id];
      treeCircAL[id+Oparams.parnum] = id;
      treeCircAL[treeCircAR[id]] = id + Oparams.parnum;
      treeCircAR[id] = id + Oparams.parnum ;
      treeCircBR[id+Oparams.parnum] = treeCircBR[id];
      treeCircBL[id+Oparams.parnum] = id;
      treeCircBL[treeCircBR[id]] = id + Oparams.parnum;
      treeCircBR[id] = id+Oparams.parnum;
#endif
    }
#ifdef MD_CALENDAR_HYBRID
  initHQlist();
#endif
}
#else
void InitEventList (void) 
{
  int id;
  treeLeft[0] = treeRight[0] = -1;
  treeIdA[0] = Oparams.parnum + 1;
  /* i nodi da Oparams.parnum + 1 (compreso) in poi sono il pool (cioè quelli dinamici) */
#ifdef MD_BIG_DT
  for (id = 1; id < poolSize; id++) 
    treeUp[id] = -1;
#endif
  for (id = treeIdA[0]; id <= poolSize - 2; id++) 
    treeCircAR[id] = id + 1;
  treeCircAR[poolSize-1] = -1;
  for (id = 1; id <= Oparams.parnum; id++) 
    {
      treeCircAL[id] = treeCircBL[id] = id;
      treeCircAR[id] = treeCircBR[id] = id;
    }
#ifdef MD_CALENDAR_HYBRID
  initHQlist();
#endif
}
#endif

#ifdef MD_CALENDAR_HYBRID
#if 0
extern eventQEntry *eventQEntries;

int *CBT; /* complete binary tree implemented in an array of
	     2*N integers */
int NP=0; /*current number of events */
int evIdHQ, noinsertHQ;
#endif

int insertInEventQ(int p)
{
  int oldFirst;
  double idbl;
  double IBIG = 1.0E18; /* il massimo long long int = +9,223,372,036,854,775,807 > IBIG */	
  long long int i;
  //eventQEntry * pt;
  //pt=eventQEntries+p; /* use pth entry */
  /* NOTA baseIndex va messo in OprogStatus! */
  idbl=(OprogStatus.scaleHQ*treeTime[p]-OprogStatus.baseIndex);

  /* very big numbers go to overflow list directly avoiding long long int overflows
     and segfaults */
  if (idbl > IBIG)
    {
      i=OprogStatus.nlistsHQ; /* store in overflow list */
      overevHQ++;
    }
  else
    {
      i = (long long int) idbl;
      /* N.B. se scaleHQ è grande un int qui puo' non bastare e ci si becca
	 un segfault! */
      /* O(1) is disabled forcing i to be currentIndex
	 i=currentIndex;*/
      if (OprogStatus.scaleHQ  < 0)
	i = OprogStatus.curIndex;
      //printf("baseIndex=%.15G p=%d i=%d\n", baseIndex, p, i);
      if(i>(OprogStatus.nlistsHQ-1)) /* account for wrap */
	{
	  i-=OprogStatus.nlistsHQ;
	  if(i>=OprogStatus.curIndex-1)
	    {
	      i=OprogStatus.nlistsHQ; /* store in overflow list */
	      overevHQ++;
	    }
	}
    }
  //pt->qIndex=i;
  totevHQ++;
  treeQIndex[p] = i;
  if(i==OprogStatus.curIndex)
    {
      InsertPQ(p); /* insert in PQ */
      return p;
    }
  else
    {
      //numevLL++; /* numero eventi nelle linked lists*/
      /* insert in linked list */
      treeStatus[p] = 1; /* 1 = belonging to linked lists of HQ */
    
      oldFirst=linearLists[i];
      MD_DEBUG2(printf("Inserting in linked lists oldFirst=%d p=%d idA=%d idB=%d\n", oldFirst,p,treeIdA[p],
		       treeIdB[p]));
      treePrev[p] = -1; /* treeLeft = previous */
      //pt->previous=-1;
      treeNext[p] = oldFirst; /* treeRight = next */
      //pt->next=oldFirst;
      linearLists[i]=p;
      if(oldFirst!=-1)
	treePrev[oldFirst] = p;
	//eventQEntries[oldFirst].previous=p;
    }
  return p;
}

void processOverflowList(void)
{
  int i,e,eNext;
  i=OprogStatus.nlistsHQ; /* overflow list */
  e=linearLists[i];
  linearLists[i]=-1; /* mark empty; we will treat all entries and may re-add some */

  overevHQ=0;
  while(e!=-1)
    {
      eNext = treeNext[e];
      //eNext=eventQEntries[e].next; /* save next */
      insertInEventQ(e); /* try add to regular list now */
      e=eNext;
    }
}
void deleteFromEventQ(int e)
{
  int prev,next,i;
  //eventQEntry *pt=eventQEntries+e;
  //i=pt->qIndex;
  totevHQ--;
  i = treeQIndex[e];
  if (i==OprogStatus.nlistsHQ)
    overevHQ--;
  treeStatus[e] = 0; /* free node */
  if(i==OprogStatus.curIndex)
    {
      MD_DEBUG2(printf("[delete] e=%d PQ node\n", e));
      DeletePQ(e); /* delete from pq */
    }
  else
    {
      /* remove from linked list */
      MD_DEBUG(printf("[delete] e=%d PQ node\n", e));
      prev = treePrev[e];
      next = treeNext[e];
      //prev=pt->previous;
      //next=pt->next;
      if(prev==-1)
	linearLists[i] = treeNext[e];
	//linearLists[i]=pt->next;
      else
	treeNext[prev]=next;
	//eventQEntries[prev].next=next;
      if(next!=-1)
	treePrev[next] = prev;
	//eventQEntries[next].previous=prev;
    }
}
//int deleteFirstFromEventQ()
int populatePQ(void)
{
  int e=-1;
  while(treeRight[0]==-1)
    /*if priority queue exhausted, i.e. if binary tree calendar is void */
    {
      /* change current index */
      OprogStatus.curIndex++;
#if 0
      if (OprogStatus.refTime > 0.0)
    	printf("currentIndex=%d feeding PQ linearLists[]=%d\n", currentIndex, linearLists[currentIndex]);
#endif
      MD_DEBUG2(printf("currentIndex=%d feeding PQ linearLists[]=%d\n", currentIndex, linearLists[currentIndex]));
      //printf("currentIndex=%d feeding PQ linearLists[]=%d baseIndex=%d\n", currentIndex, linearLists[currentIndex],
	//     OprogStatus.baseIndex);
      if(OprogStatus.curIndex==OprogStatus.nlistsHQ)
	{
	  OprogStatus.curIndex=0;
	  OprogStatus.baseIndex+=OprogStatus.nlistsHQ;
	  //printf("baseIndex=%f\n", baseIndex);
	  processOverflowList();
	}
      /* populate pq */
      e=linearLists[OprogStatus.curIndex];
      while(e!=-1)
	{
	  InsertPQ(e);
	  e = treeNext[e];
	  //e=eventQEntries[e].next;
	}
      /* 10/05/2010: notare che appena inserite gli eventi nell'albero
	 binario questi non sono più raggiungibili tramite le linked lists
	 anche se hanno treeQIndex[]=currentIndex */
      linearLists[OprogStatus.curIndex]=-1;
    }
  /* delete from binary tree here! */
#if 0
  e=CBT[1]; /* root contains shortest time entry */
  Delete(CBT[1]);
  /* a prendere l'evento successivo ci pensa NextEvent() 
     all'interno della quale all'inizio viene chiamata deleteFirstFromEventQ()
     per popolare il binary tree se vuoto */
#endif
  return e;
}
#endif
