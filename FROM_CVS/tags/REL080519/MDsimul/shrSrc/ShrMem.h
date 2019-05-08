#include<sys/types.h>
#include<sys/ipc.h>
#include<sys/shm.h>
#include<sys/sem.h>
#include<stdio.h>
#include<unistd.h>
#define NUMSHAREBLOCK 90
#define KEY_PATH "/tmp" 
/* normalmente si usa . tuttavia cosi' la chiave e' assoluta. Cio' puo' 
   costituire un rischio se qualche altro processo usa /tmp come percorso 
   per la chiave */

/* dimensione della memoria condivisa 
   MAX 4Mb */
/* unione che serve per settare il valore iniziale del semaforo */ 
union semun semopts; 
/* crea delle strutture per il lock l'unlock del semaforo  */
struct sembuf sem_lock = {0, -1, 0};
struct sembuf sem_unlock = {0, 1, IPC_NOWAIT };
int shmid[NUMSHAREBLOCK];
int stato,semid; 
void *shmem_res_list[NUMSHAREBLOCK];
int FIRST_TIME = 2;
int free_shr;
/* prototipi di tutte le funzioni private e pubbliche della libreria share */
void del_all(void);
void *getshmem_r(int nbloc,int size);
void *getshmem(int size);
void createsem(int nsems);
void removesem(void);
void lock(ushort numerosem);
void unlock(ushort numerosem);
int getsemval(int numerosem);
void freeshmem_r(int nbloc);
void freeshmem(void* punta);
void delallshm(void);






