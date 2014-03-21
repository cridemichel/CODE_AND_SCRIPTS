/* Ogni elemento dell'array shmem_res_list contiene un puntatore ad un blocco
   di memoria condivisa allocato.
   Tale puntatore e' restituito in fase di allocazione dalla funzione shmget();
   se un elemento contiene NULL allora quel blocco non e' stato allocato. 
   Getshmem_r() e freeshmem_r() sono versioni rudimentali che necessitano 
   del passaggio di un valore intero (BlockA,BlockB,ecc.) compreso tra 
   0 e NUMSHAREBLOCK-1, l'intero passato in tale caso indetifica la risorsa 
   cosicche' getshmem_r(BlockA,1000) freeshmem_r(BlockA) allocano e deallocano
   1k di memoria condivisa.
   Per cio' che riguarda getshmem() e freeshmem() invece e' sufficiente
   notare che si comportano come malloc() e free() */
#include<ShrMem.h>

/* la seguente funzione quando invocata cancella tutte le risorse shared 
   memory  allocate fino al momento della sua invocazione */
void del_all(void)
{
  int i;
  for (i=0; i<NUMSHAREBLOCK; ++i)
    if (shmem_res_list[i]!=NULL) freeshmem(shmem_res_list[i]);
  printf("All shared memory resources eliminated\n");
}

void *getshmem_r(int nbloc,int size) 
{
key_t key;
void* segptr;
char keycar;
keycar = 64 + nbloc;
/* il processo che usa la memoria condivisa non dovrebbe cambiare la 
directory */ 

/* 19/3/99 CHANGED (OLD: key=ftok(KEY_PATH,keycar); doesn't work)*/
key = ftok(".", keycar);

/* crea la risorsa sharedmemory */
if ((shmid[nbloc] = shmget(key,size,IPC_CREAT| IPC_EXCL |0644)) == -1)
  {
    printf("Shared memory segment exists - trying as client\n");
    if ( (shmid[nbloc] = shmget(key,size,0)) == -1 )
      {
	perror("shmget");
	printf("Remove manually shared memory blocks with ipcrm.\n");
	printf("Use delshm utility\n");
	exit(1);
      }
  }
else 
  {
    printf("Creating new shared memory segment\n");
  }

/* crea un puntatore alla memoria condivisa */
if ((segptr = shmat(shmid[nbloc], 0, 0))== (void *)-1)
  {
    perror("shmat");
    exit(1);
  }
return segptr;
}
void *getshmem(int size)
{ 
  int i;
  /* se e' la prima volta nel programma che si chiama getshrmem inizializza 
     a false (0) tutti gli elementi dell'array shmem_res_list */ 
  if (FIRST_TIME==1) 
    {
      for (i=0; i<NUMSHAREBLOCK; ++i) shmem_res_list[i]=NULL;
      FIRST_TIME=0;
    }
  /* cerca una risorsa shmem libera nell'array */
  for(i=0; i<NUMSHAREBLOCK; ++i) 
      if (shmem_res_list[i]==NULL) /* NULL = blocco non allocato */
	{
	  free_shr=i; /* il blocco da allocare e' l'i-esimo */
	  break; /* esce dal for */
	}
  /* se si sono allocati gia' NUMSHAREBLOCK blocchi di memoria condivisa allora
     si incappa in un errore */
  if (i==NUMSHAREBLOCK)
    {
      printf("ERROR: You can't allocate more than %d shared memory blocks\n",NUMSHAREBLOCK);
      del_all(); /* rimuove tutte le risorse 'shared memory'
		    allocate fino ad ora */
      exit(-1);
    }
  else 
    {
      /* alloca la memoria condivisa richiesta */ 
      shmem_res_list[i] = getshmem_r(free_shr,size);
      /* restituisce il puntatore del blocco di memoria condivisa allocato */ 
      return shmem_res_list[i];
    }
}

void freeshmem_r(int nbloc)
{
/* rimuove la risorsa shared memory */
  shmctl(shmid[nbloc], IPC_RMID, 0);
  printf("Shared memory segment marked for deletion\n");
  fflush(stdout);
}
 
void delallshm(void)
{ 
  register int i; 
  for (i=0; i<NUMSHAREBLOCK; ++i)
    if (shmem_res_list[i]!=NULL) 
      freeshmem_r(i);
}

void freeshmem(void* punta)
{
  int i; 
  for (i=0; i<NUMSHAREBLOCK; i++)
    {
      /* se il puntatore passato come parametro e' nell'elemento i-esimo 
	 dell'array allora rimuovi il blocco i-esimo tramite 
	 freeshmem_private */
      if (shmem_res_list[i]==punta) 
	{
	  freeshmem_r(i);
	  break;
	}
    }

if (i==NUMSHAREBLOCK) 
  {
    /* se si arriva qui il puntatore indicato in freshmem() non appartiene
       all'array shmem_res_list, ovvero non si riferisce ad un blocco 
       di memoria condivisa precedentemente allocato */
    printf("ERROR: Not existing shared memory block\n");
    del_all(); /* rimuove tutte le risorse 'shared memory'
		  allocate fino ad ora */
    exit(-1);
  }
}

/* la seguente procedura crea un set di nesems semafori */
void createsem(int nsems)
{
  int ii;
  key_t key;
  key=ftok(".",'S');
  /* crea un set di due semafori */
  if ( (semid = semget(key, nsems, IPC_CREAT | 0644)) == -1 )
    {
      perror("semget");
      printf("Unable to get semaphore set id\n");
      exit(-1);
    }
  /* setta il valore iniziale dei semafori a zero */
  semopts.val = 0;
  for (ii=0; ii<nsems; ++ii) 
    { 
      semctl(semid, ii, SETVAL, semopts); /* semaforo del proc padre */
    }
}

/* Set the the value of a semaphore */
void setsem(int nsem, int val) 
{
  semopts.val = val;
  semctl(semid, nsem, SETVAL, semopts);
}

/* rimuove  il set di semafori creato con createsem() */
void removesem(void)
{
  /* rimuove la risorsa semaforo */
   semctl(semid, 0, IPC_RMID, 0); 
   printf("Semaphores set removed\n");
}

/* numersem = 0,1,2,..NUMERO_SEMAFORI 
dove NUMERO_SEMAFORI e' il numero di semafori creati con createsem(). 
Lock decrementa di 1 il semaforo specificato da numerosem, ma
se un certo semaforo vale zero e si chiama lock allora il processo dorme 
(sleep) finche' il suo valore non diviene positivo (unlock), cioe' finche' non 
viene fatta una unlock.
Questo perche' non si e' specificato il flag IPC_NOWAIT nella struttura 
sem_lock.
Tale ultima caratteristica permette di utilizzare i semafori per la 
sincronizzazione di due processi.
*/

void lock(ushort numerosem)
{
  sem_lock.sem_num = numerosem;
  /* lock del semaforo specificato da numerosem (decrementa di 1)*/
  if ( semop(semid, &sem_lock, 1) == -1 ) perror("semop"); 
} 

/* unlock() incrementa di 1 il semaforo specificato da numerosem */
void unlock(ushort numerosem)
{
  sem_unlock.sem_num = numerosem;
  /* unlock (incrementa di 1) del semaforo del padre */
  if ( semop(semid, &sem_unlock, 1) == -1) perror("semop");
}

int getsemval(int numerosem)
{
  return semctl(semid,numerosem,GETVAL,0);
}
