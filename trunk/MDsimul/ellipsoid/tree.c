/* tree.c */
  
/********************************************************************
   This software is copyrighted material, reproduced from the book 
   "The Art of Molecular Dynamics Simulation" by D. C. Rapaport, 
   published by Cambridge University Press (1995).
*********************************************************************/
#define MD_DEBUG34(x) 
#include <mdsimul.h>
extern int **tree, evIdA, evIdB;
#ifdef MD_PATCHY_HE
extern int evIdC, evIdD, evIdE;
#endif
extern double *treeTime, *treeRxC, *treeRyC, *treeRzC;
void DeleteEvent(int );
void NextEvent(void);
extern int poolSize;
extern double rxC, ryC, rzC;
#ifdef MD_SPHERICAL_WALL
extern int sphWall, sphWallOuter;
#endif
#if defined(MD_MULTIPLE_LL) 
extern int **crossevtodel;
#endif

#define MD_DEBUG20(x) 
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
void enlarge_queue_size(void);

void ErrExit(char *str)
{
  printf(str);
  printf("\n");
  enlarge_queue_size();
  //exit(-1);
}
int check_node(char* str, int id, int idNew, int idUp)
{
  if (Oparams.curStep==1266)
    {
      printf("(id:%d,u:%dr:%d,l:%d) ",
	     id, treeUp[id], treeRight[id], treeLeft[id]);
    }
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
/* Allocate memory for a matrix of integers */
int** ReallocMatI(int **v, int size1, int size2, int size2old)
{
  int k, kk;
  int **vaux;
  vaux = (int **) malloc(sizeof(int*)*size1);
  vaux[0] = (int*) malloc(size1 * size2 * sizeof(int));
  if (!vaux || !vaux[0])
    {
      printf("[CRITICAL ERROR] I tried to enlarge event queue but malloc() failed\n");
      exit(-1);
    }

  for (k = 1; k < size1; k++)
    vaux[k] = vaux[k-1] + size2;
    /* copy old tree onto new one */
  for (k = 0; k < size1; k++)
    for (kk = 0; kk < size2old; kk++)
      vaux[k][kk] = v[k][kk];
  /* free old allocated memory */  
  free(v[0]);
  free(v);
  
  return vaux;
}

void enlarge_queue_size(void)
{
  double dblps;
  int oldps, id;
  oldps = poolSize;
  dblps = (double) poolSize;
  while (((int)dblps)==poolSize)
   dblps *= 1.10; /* 10% increment */ 
  poolSize = (int) dblps;
//  printf("treeCircAR[oldps-1]=%d treeCircAR[%d]=%d\n", treeCircAR[oldps-1],idd, treeCircAR[idd]);

#if defined(MD_PATCHY_HE) || defined(EDHE_FLEX)
#ifdef MD_HE_PARALL
  if (my_rank == 0)
    tree = ReallocMatI(9, poolSize, oldps);
#else
#ifdef MD_CALENDAR_HYBRID
  tree = ReallocMatI(tree, 16, poolSize, oldps);
#else
  tree = ReallocMatI(tree, 12, poolSize, oldps);
#endif
#endif
#else
#ifdef MD_HE_PARALL
  if (my_rank == 0)
    tree = RellocMatI(tree, 9, poolSize, oldps);
#else
#ifdef MD_CALENDAR_HYBRID
  tree = RellocMatI(tree, 13, poolSize, oldps);
#else
  tree = ReallocMatI(tree, 9, poolSize, oldps);
#endif
#endif
#endif  

#ifdef MD_HE_PARALL
  if (my_rank == 0)
    {
      treeTime = realloc(treeTime,sizeof(double)*poolSize);
      treeRxC  = realloc(treeRxC,sizeof(double)*poolSize);
      treeRyC  = realloc(treeRyC,sizeof(double)*poolSize);
      treeRzC  = realloc(treeRzC,sizeof(double)*poolSize);
    }
#else
  treeTime = realloc(treeTime,sizeof(double)*poolSize);
  treeRxC  = realloc(treeRxC,sizeof(double)*poolSize);
  treeRyC  = realloc(treeRyC,sizeof(double)*poolSize);
  treeRzC  = realloc(treeRzC,sizeof(double)*poolSize);
#endif
  if (treeTime==NULL || treeRxC==NULL || treeRyC == NULL || treeRzC == NULL)
    {
      printf("[CRITICAL ERROR] I tried to enlarge queue but realloc() failed\n");
      exit(-1);
    }
  else
    printf("[INFO] I have just enlarged event queue to %d nodes (old size was %d)\n", poolSize, oldps);
  /* inizializza opportunamente i nuovi nodi, ossia li inserisce tutti nel pool di quelli liberi */
  for (id=oldps; id <= poolSize-2; id++) 
    treeCircAR[id] = id + 1;
  treeCircAR[poolSize-1]=-1;
  /* oldps is a newly allocated free node, hence now treeIdA[0] is no more -1 but points 
     to a available node */
  treeIdA[0] = oldps;
  /* reset status of new nodes 0 = free (this is needed for timeshift_calendar()) */
  for (id=oldps; id < poolSize; id++)
    treeStatus[id] = 0;
}
#ifdef MD_CALENDAR_HYBRID
int get_new_node(int idA, int idB, int idata)
{
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
	printf("idNew: %d tEvent: %.15f\n", idNew, tEvent));
      /* all'inizio treeCircAR[treeIdA[0]] = treeIdA[0]+1 quindi è un nodo 
       * non utilizzato nel pool */
      treeIdA[0] = treeCircAR[treeIdA[0]];
    }
  else 
    {
#ifdef MD_MULTIPLE_LL
      /* notare comunque che se idata=typeOfPart[idA] allora
	 l'assegnazione è inutile */
      /* idata si riferisce al tipo d'interazione relativa al cell crossing attuale
	 ossia la cella attraversata è quella relativa all'interazione tra typeOfPart[idA] e idata. */
      idNew = idA + 1 + idata*Oparams.parnum;
      if (OprogStatus.multipleLL && idata != typeOfPart[idA])
	{
	  crossevtodel[idata][idA] = idNew;  
	}
#else
      idNew = idA + 1;
#endif
    }
  MD_DEBUG34(printf("idNew=%d\n", idNew));
  /* Se qui vuol dire che si tratta di un cell-crossing o 
     di un urto con parete
     NOTA: urto con parete e cell-crossing sono esclusivi, per cui basta un nodo 
     inoltre c'è sempre un evento di tale tipo associato con ogni particella 
     */
  return idNew;
}
void InsertPQ(int idNew)
{
  int more, id;
  double tEvent;
  numevPQ++;
  tEvent = treeTime[idNew];
  treeStatus[idNew] = 2; /* 2 = belonging to binary tree */
  id = 0;
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
    
  treeLeft[idNew] = treeRight[idNew] = -1;
  treeUp[idNew] = id;
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
#endif
#ifdef MD_PATCHY_HE
#ifdef MD_CALENDAR_HYBRID
void setPQnode(int idNew, int idA, int idB, int idata, int idatb, int idcollcode, double tEvent, double rxC, double ryC, double rzC)
{
  treeTime[idNew] = tEvent;
  treeRxC[idNew] = rxC;
  treeRyC[idNew] = ryC;
  treeRzC[idNew] = rzC;
  treeIdA[idNew] = idA;    
  treeIdB[idNew] = idB;
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
}
void ScheduleEventBarr (int idA, int idB, int idata, int idatb, int idcollcode, double tEvent)
{
  int idNew;
  idNew = get_new_node(idA, idB, idata);
  /* assegna i dati al nodo relativi all'evento */
  setPQnode(idNew, idA, idB, idata, idatb, idcollcode, tEvent, rxC, ryC, rzC);
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

  MD_DEBUG34(printf("#%lld ScheduleEvent() idA:%d idB:%d evtime:%.15f\n", 
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
    {
#ifdef MD_MULTIPLE_LL
      /* notare comunque che se idata=typeOfPart[idA] allora
	 l'assegnazione è inutile */
      /* idata si riferisce al tipo d'interazione relativa al cell crossing attuale
	 ossia la cella attraversata è quella relativa all'interazione tra typeOfPart[idA] e idata. */
      idNew = idA + 1 + idata*Oparams.parnum;
      if (OprogStatus.multipleLL && idata != typeOfPart[idA])
	{
	  crossevtodel[idata][idA] = idNew;  
	}
#else
      idNew = idA + 1;
#endif
    }
  MD_DEBUG34(printf("idNew=%d\n", idNew));
  /* Se qui vuol dire che si tratta di un cell-crossing o 
     di un urto con parete
     NOTA: urto con parete e cell-crossing sono esclusivi, per cui basta un nodo 
     inoltre c'è sempre un evento di tale tipo associato con ogni particella 
     */
  /* treeRight[id] == -1 => il calendario è vuoto */
  treeStatus[idNew] = 2; /* node belonging to binary tree */

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
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
  treeLeft[idNew] = treeRight[idNew] = -1;
  treeUp[idNew] = id;
}
#endif
void ScheduleEvent(int IdA, int IdB, double tEvent)
{
  ScheduleEventBarr(IdA, IdB, 0, 0, MD_EVENT_NONE, tEvent);
}
#else
#ifdef MD_CALENDAR_HYBRID
void setPQnode(int idNew, int idA, int idB, int idata, int idatb, int idcollcode, double tEvent, double rxC, double ryC, double rzC)
{
  treeTime[idNew] = tEvent;
  treeRxC[idNew] = rxC;
  treeRyC[idNew] = ryC;
  treeRzC[idNew] = rzC;
  treeIdA[idNew] = idA;    
  treeIdB[idNew] = idB;
  treeIdC[idNew] = idata;
  treeIdD[idNew] = idatb;
  treeIdE[idNew] = idcollcode;
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
#ifdef EDHE_FLEX
int check_next_event(void)
{
  int idNow, idLast, fine;
  double evtime;
  idNow = treeRight[0];  
  /* Cerca la prossima collisione con tempo minore 
   * NOTA: l'albero è ordinato e ogni nodo sinistro ha un tempo inferiore */

  while (!(treeIdB[idNow] >= ATOM_LIMIT + 100 || treeIdB[idNow] < ATOM_LIMIT + NDIM * 2))
    {
      fine = 0;
      idLast = idNow;
      idNow = treeRight[0];
      while (!fine)
	{
	  if (idNow == idLast)
	    idNow = treeRight[idNow];
	  else
	    idNow = treeLeft[idNow];

	  if (idNow == -1)
	    fine=1;
	}
    }
  /* se non è una collisione ignoralo (return 0)*/ 

  evtime = treeTime[idNow];
  /* se la collisione è molto ravvicinata attenzione! (return 1) */
  if (fabs(evtime-Oparams.time) < 1E-14)
    {
      return 1;
    }
  return 0;
}
#endif
#ifdef ED_PARALL_DD
extern int check_tstep(double time);
extern int check_BZ_event(int i);
#endif
void NextEvent (void) 
{
  int id, idAx, idBx, idd, idNow, idtx;
  /* Il nodo root (0), a cui non è associato alcun evento,
   * è linkato con il suo right pointer al primo nodo che contiene
   * un evento */
#ifdef MD_MULTIPLE_LL
  int nc;
#endif
#ifdef MD_CALENDAR_HYBRID
  populatePQ();
#endif
  idNow = treeRight[0];  
  /* Cerca l'evento con tempo minore 
   * NOTA: l'albero è ordinato e ogni nodo sinistro ha un tempo inferiore */
  while (treeLeft[idNow] > -1) 
    idNow = treeLeft[idNow];
  Oparams.time = treeTime[idNow];   
#ifdef ED_PARALL_DD
  if ( check_tstep(Oparams.time) || (check_BZ_event(evIdA)||check_BZ_event(evIdB)) )
    {
      /* after check_tstep calendar may be different, hence
	 reread next event from it */
#ifdef MD_CALENDAR_HYBRID
      populatePQ();
#endif
      idNow = treeRight[0];  
      /* Cerca l'evento con tempo minore 
       * NOTA: l'albero è ordinato e ogni nodo sinistro ha un tempo inferiore */
      while (treeLeft[idNow] > -1) 
	idNow = treeLeft[idNow];
      Oparams.time = treeTime[idNow]; 
    }
#endif
  rxC = treeRxC[idNow];
  ryC = treeRyC[idNow];
  rzC = treeRzC[idNow];
  evIdA = treeIdA[idNow];    
  evIdB = treeIdB[idNow];
#ifdef MD_PATCHY_HE
  evIdC = treeIdC[idNow];
  evIdD = treeIdD[idNow];
  evIdE = treeIdE[idNow];
#endif
#if 0
  if (evIdA==870)
    printf("evIdA=%d evIDB=%d evIdC=%d evIdD=%d time=%.15G\n", evIdA, evIdB, evIdC, evIdD, treeTime[idNow]);
#endif
  MD_DEBUG34(printf("[ NextEvent ] #%lld event(%d,%d) curtime:%.15G id=%d\n", (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time, idNow));
  MD_DEBUG20(printf("[ NextEvent ] #%lld event(%d,%d) curtime:%f\n", 
		   (long long int)Oparams.curStep, evIdA, evIdB, Oparams.time));
  if (evIdB < ATOM_LIMIT+2*NDIM || evIdB==ATOM_LIMIT+50)
    /* Se si tratta di un urto fra particelle o con le pareti...
     * qui in sostanza si considerano solo gli eventi che cambiano lo stato 
     * della particella, cioè la sua velocità */
    {
      /* qui incrementa di 1 poiché il nodo root non contiene eventi */
#ifdef MD_ABSORPTION
      /* NOTA: gli id che vanno da 0...Oparams.parnum sono riservati ai cell crossing
         o agli urti con le pareti della scatola, nel caso di urto con la membrana 
	 semi-permeabile tale evento va rimosso esplicitamente altrimentri rimane nel calendario 
	 degli eventi */
      if (evIdB==ATOM_LIMIT+50)
	DeleteEvent(idNow);
#endif
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
#ifdef MD_SPHERICAL_WALL
	  /* N.B. gli eventi in cui compare il muro sferico, a parte l'urto con evIdA non vanno cancellati  */
	  if (id==sphWall+1)
	    continue; 
	  if (id==sphWallOuter+1)
	    continue; 
#endif
	  MD_DEBUG34(printf("deleted event #%d\n", id));
#if defined(MD_MULTIPLE_LL) 
	  if (OprogStatus.multipleLL)
	    {
	      /* cancella tutti gli eventi di cell-crossing se si tratta di una particella */

	      /* prima cancella cell-crossing per le liste relative a particelle dello stesso tipo*/
#if 0
	      if (treeStatus[id + typeOfPart[id-1]*Oparams.parnum]==0)
		{
		  printf("evIdA=%d evIdB=%d id=%d I should delete event not in PQ?!?\n", evIdA, evIdB, id);
		  exit(-1);
		}
#endif
	      DeleteEvent(id + typeOfPart[id-1]*Oparams.parnum);
	      for (nc = 0; nc < Oparams.ntypes; nc++)
		{
		  if (nc == typeOfPart[id-1])
		    continue;
		  /* id-1 poiché idAx = evIdA (vedi sopra) */
		  /* avendo messo la condizione in delete event su treeStatus[id]
		     ossia se treeStatus[id]==0 non tenta di eliminare il nodo dalla PQ
		     visto che non è inserito, la condizione seguente è in realtà inutile */
		  if (crossevtodel[nc][id-1]!=-1)
		    {
		      DeleteEvent (id+nc*Oparams.parnum);
		    }
		  /* notare che la condizione crossevtodel[id-1]!=-1 è necessaria 
		   * in quanto se una particella attraversa le pareti del box l'evento con nc==1
		   * non viene schedulato visto che ProcessCellCrossing si occupa di entrambi gli eventi
		   * (nc==0 e nc==1) */
		  /* Serve assolutamente porre crossevtodel a -1  per evitare che ProcessCellCross 
		   * tenti di rimuovere un evento già rimosso. */
		      crossevtodel[nc][id-1] = -1;
		}
	    }
	  else
	    DeleteEvent(id);
#else
	  /* cancella l'evento di cell-crossing */
	  DeleteEvent (id);
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
#ifdef MD_MULTIPLE_LL
      if (evIdB < ATOM_LIMIT + 100) 
	{
	  treeCircAR[idNow] = treeIdA[0];
	  treeIdA[0] = idNow;
	} 
      //printf("evIdA=%d typeOfPart[]=%d\n", evIdA, typeOfPart[evIdA]);
      if (OprogStatus.multipleLL && evIdA >=0 && evIdA < ATOM_LIMIT && evIdC != typeOfPart[evIdA])
	{
	  crossevtodel[evIdC][evIdA] = -1;
	}
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
    
}

void DeletePQ (int id)
{
  int idp, idq, idr;
#ifdef MD_CALENDAR_HYBRID
  numevPQ--;
#endif
#ifdef MD_SPHERICAL_WALL
  /* N.B. sphWall+1 è l'evento di cell-crossing del 
     muro sferico ma tale evento non viene schedulato affatto
     e quindi non va neanche rimosso */
 if (id == sphWall+1)
  return; 
 if (id == sphWallOuter+1)
  return; 
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
  /* event not inserted in PQ hence it has not to be deleted this may happen
     for cell crossings of infinite mass particles!*/
  if (treeStatus[id]==0) /* 0 means not inserted in PQ */
    return;
  deleteFromEventQ(id);
}
#else
void DeleteEvent(int id)
{
  /* event not inserted in PQ hence it has not to be deleted, this may happen
     for cell crossings of infinite mass particles! */
  if (treeStatus[id] == 0) /* 0 means not inserted in PQ */
    return;
  treeStatus[id] = 0;
  DeletePQ(id);
}
#endif
#ifdef MD_CALENDAR_HYBRID
void initHQlist(void)
{
  int i;
  //OprogStatus.curIndex = 0;
  /* inizializzare anche le linked lists lineari? */
  for (i=0; i < OprogStatus.nlistsHQ+1; i++)
    linearLists[i] = -1;
#if 0
  for (i=1; i < poolSize; i++) 
    treeStatus[i] = 0;/* 0 = free node */
#endif
  numevPQ = overevHQ = totevHQ = 0; 
 
#ifdef MD_CALENDAR_HYBRID
  /* sembra che sia necessario per ricostruire il calendario, boh...*/
  OprogStatus.curIndex=0;
  /* NOTA 11/10/2010: l'importante è che: baseIndex < OprogStatus.scaleHQ*Oparams.time */  
  OprogStatus.baseIndex = Oparams.time*OprogStatus.scaleHQ;
  //printf("currentIndex=%d baseIndex=%d\n", OprogStatus.curIndex, OprogStatus.baseIndex);
#endif
	
}
#endif

void InitEventList (void) 
{
  int id, i;
  treeLeft[0] = treeRight[0] = -1;
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    treeIdA[0] = Oparams.ntypes*Oparams.parnum + 1;
  else
    treeIdA[0] = Oparams.parnum + 1;
#else
  treeIdA[0] = Oparams.parnum + 1;
#endif
  /* i nodi da Oparams.parnum + 1 (compreso) in poi sono il pool (cioè quelli dinamici) */
  for (id = treeIdA[0]; id <= poolSize - 2; id++) 
    treeCircAR[id] = id + 1;
#ifdef MD_BIG_DT
  for (id = 1; id < poolSize; id++) 
    treeUp[id] = -1;
#endif
  for (i=1; i < poolSize; i++) 
    treeStatus[i] = 0;/* 0 = free node */
  treeCircAR[poolSize-1] = -1;
#ifdef MD_MULTIPLE_LL
  if (OprogStatus.multipleLL)
    {
      for (id = 1; id <= Oparams.parnum*Oparams.ntypes; id++) 
	{
	  treeCircAL[id] = treeCircBL[id] = id;
	  treeCircAR[id] = treeCircBR[id] = id;
	}
    }
  else
    {
      for (id = 1; id <= Oparams.parnum; id++) 
	{
	  treeCircAL[id] = treeCircBL[id] = id;
	  treeCircAR[id] = treeCircBR[id] = id;
	}
    }

#else
  for (id = 1; id <= Oparams.parnum; id++) 
    {
      treeCircAL[id] = treeCircBL[id] = id;
      treeCircAR[id] = treeCircBR[id] = id;
    }
#endif
#ifdef MD_CALENDAR_HYBRID
  initHQlist();
#endif
}
#ifdef MD_CALENDAR_HYBRID
int insertInEventQ(int p)
{
  int oldFirst;
  double idbl;
  double IBIG = 1.0E18;
  long long i;
  //eventQEntry * pt;
  //pt=eventQEntries+p; /* use pth entry */
  //i=(int)(OprogStatus.scaleHQ*treeTime[p]-OprogStatus.baseIndex);
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
      /* se si scommenta questa riga di fatto si disattiva il calendario O(1) */
      //printf("curIndex=%d baseIndex=%d\n", OprogStatus.curIndex, OprogStatus.baseIndex);
      /* disable O(1) if scaleHQ < 0 */
      if (OprogStatus.scaleHQ < 0)
	i=OprogStatus.curIndex;

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
      treeStatus[p] = 1; /* 1 = belonging to linked lists of HQ */
      /* insert in linked list */
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
  overevHQ = 0;
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
  i = treeQIndex[e];
  if (i==OprogStatus.nlistsHQ)
    overevHQ--;
  treeStatus[e] = 0;
  totevHQ--;
#if 0
#ifdef MD_SPHERICAL_WALL
  /* N.B. sphWall+1 è l'evento di cell-crossing del 
     muro sferico ma tale evento non viene schedulato affatto
     e quindi non va neanche rimosso */
  if (id == sphWall+1)
    return; 
  if (id == sphWallOuter+1)
    return; 
#endif
#endif
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
      //printf("currentIndex=%d feeding PQ linearLists[]=%d nlists=%d\n", OprogStatus.curIndex, linearLists[OprogStatus.curIndex], OprogStatus.nlistsHQ);
      MD_DEBUG2(printf("currentIndex=%d feeding PQ linearLists[]=%d\n", currentIndex, linearLists[currentIndex]));
      if(OprogStatus.curIndex==OprogStatus.nlistsHQ)
	{
	  OprogStatus.curIndex=0;
	  OprogStatus.baseIndex+=OprogStatus.nlistsHQ;
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
