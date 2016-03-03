/* FATHER = semaforo del processo padre dopo una fork 
   CHILD = semaforo del processo figlio dopo una fork
   
   */
#define FATHER 1
#define CHILD  0

#define SPLIT pid = fork(); 
#define BEGPROC1  if (pid==0) {
#define ENDPROC1  }
#define BEGPROC2  else {
#define ENDPROC2  }
#define BlockA 0
#define BlockB 1
#define BlockC 2
#define BlockD 3
#define BlockE 4
#define BlockF 5
#define BlockG 6
#define BlockH 7
#define BlockI 8
#define BlockL 9
#define BlockM 10
#define BlockN 11
#define BlockO 12
#define BlockP 13
#define BlockQ 14
#define BlockR 15
#define BlockS 16
#define BlockT 17
#define BlockU 18
#define BlockV 19
#define BlockZ 20

/* active semaphores made by using shared memory, you must supply 
   a shared memory address of an array of two integers 
   SEM_ADDR = pointer to an array of two integers
   WHICH = 0 or 1, specify which semaphore to lock or unlock 
*/ 

#define lockSh(SEM_ADDR, WHICH) while(SEM_ADDR[WHICH] == 0);\
                                --(SEM_ADDR[WHICH])

#define unlockSh(SEM_ADDR, WHICH) ++SEM_ADDR[WHICH]
                                  

/* to avoid problems you must use these doubly, that is
   lock/unlock/lock/unlock, using two different semaphores for each 
   lock/unlock pair */ 
#define lockSh10(SEM_ADDR,WHICH) while(SEM_ADDR[WHICH] == 0);\
                                  SEM_ADDR[WHICH]=0
				   
#define unlockSh10(SEM_ADDR, WHICH) SEM_ADDR[WHICH]=1

#include<sys/types.h>
/* vedere il file ShrMem.c per i dettagli delle varie primitive */
void* getshmem(int size);
void* getshmem_r(int nbloc,int size);
void freeshmem(void *punta);
void freeshmem_r(int nbloc);
void createsem(int nsems);
void removesem(void);
void lock(ushort numerosem);
void unlock(ushort numerosem);
int getsemval(int numerosem);
void setsem(int nsem, int val);
void delallshm(void); 














