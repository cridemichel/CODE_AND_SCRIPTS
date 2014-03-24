/* ========================== >>> struct params <<< =========================*/

struct paramsBM
{
  /* DESCRIPTION:
     This is the header put at the top of a coordinate file. 
     It is strongly dependent upon simulation, besides to make 
     it more general as possible, it has some unused attribs */

  /* =================== >>> DON'T TOUCH THESE !!!! <<< =================== */
  int parnum[NA];        	/* particles number of
				 species 0 (A) and 1(B)*/
  COORD_TYPE steplength;		/* temporal step length */
  int totStep;	/* temporal step number that simulation 
				   must do */
  int curStep;	/* current step of simulation */
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  COORD_TYPE d;                 /* distance between atoms */
  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m[NA];             /* atoms masses */
  COORD_TYPE rcut;              /* cutoff for the pair potential */ 
  COORD_TYPE sigab[NA][NA];     /* pair potential length parameters */
  COORD_TYPE epsab[NA][NA];     /* pair potential energy parameters */
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
  
  COORD_TYPE tol;               /* Tolerance of the shake algoritm used 
				   by RATTLE */
  /* questi li ho messi per poter leggere i .cor file prodotti 
     dalla simulazione MC */
  int lattice_M;                /* Numero di celle su di un asse */
  double lattice_a;             /* passo del reticolo */
  int num_neigh;                
  /* Ogni volta che tenta una mossa nella MC deve scegliere dei vicini
     questo parametro indica quanti bisogna considerarne lungo un asse 
     (1,2,..) */
  
  /*======================================================================= */

} OparamsBM;

#define BM_SAVE_LISTA rxBM[0], ryBM[0], rzBM[0],\
                  vxBM[0], vyBM[0], vzBM[0],\
                  FxBM[0], FyBM[0], FzBM[0],\
                  vxo1BM[0], vyo1BM[0], vzo1BM[0],\
                  vxo2BM[0], vyo2BM[0], vzo2BM[0]
#define BM_SAVE_LISTB rxBM[1], ryBM[1], rzBM[1],\
                   vxBM[1], vyBM[1], vzBM[1],\
                   FxBM[1], FyBM[1], FzBM[1],\
                   vxo1BM[1], vyo1BM[1], vzo1BM[1],\
                   vxo2BM[1], vyo2BM[1], vzo2BM[1]


#define BM_EXT_SLST  &sBM, &s1BM, &s2BM, &VolBM, &Vol1BM, &Vol2BM, &Vol1o1BM, &s1o1BM, &Vol1o2BM,&s1o2BM

#define BM_ALLOC_LISTA  &rxBM[0], &ryBM[0], &rzBM[0],\
                    &vxBM[0], &vyBM[0], &vzBM[0],\
                    &FxBM[0], &FyBM[0], &FzBM[0],\
                    &vxo1BM[0], &vyo1BM[0], &vzo1BM[0],\
                    &vxo2BM[0], &vyo2BM[0], &vzo2BM[0]
                
#define BM_ALLOC_LISTB &rxBM[1], &ryBM[1], &rzBM[1],\
                    &vxBM[1], &vyBM[1], &vzBM[1],\
                    &FxBM[1], &FyBM[1], &FzBM[1],\
                    &vxo1BM[1], &vyo1BM[1], &vzo1BM[1],\
                    &vxo2BM[1], &vyo2BM[1], &vzo2BM[1]

#define BM_DECL_LIST *rxBM[NA], *ryBM[NA], *rzBM[NA],\
                    *vxBM[NA], *vyBM[NA], *vzBM[NA],\
                    *FxBM[NA], *FyBM[NA], *FzBM[NA],\
                    *vxo1BM[NA], *vyo1BM[NA], *vzo1BM[NA],\
                    *vxo2BM[NA], *vyo2BM[NA], *vzo2BM[NA]
            				   
#define BM_EXT_DLST    sBM, s1BM, s2BM, VolBM, Vol1BM, Vol2BM, Vol1o1BM, s1o1BM, Vol1o2BM, s1o2BM

