/* ========================== >>> struct params <<< =========================*/

struct paramsMO
{
  /* =================== >>> DON'T TOUCH THESE !!!! <<< =================== */
  int parnum;        	/* particles number */
  COORD_TYPE steplength;		/* temporal step length */
  int totStep;	/* temporal step number that simulation 
				   must do */
  int curStep;	/* current step of simulation */
  /* ======================================================================= */
 
  /* ==================== >>> PUT HERE YOUR PARAMS <<< ===================== */

  COORD_TYPE P;			/* pressure */
  COORD_TYPE T;			/* temperature */
  COORD_TYPE m;             /* atoms masses */
  COORD_TYPE rcut;              /* cutoff for the pair potential */ 
  COORD_TYPE sigma;     /* pair potential length parameters */
  COORD_TYPE epsilon;     /* pair potential energy parameters */
  int equilibrat;               /* != 0 if equilibrating */
  int M;                        /* number of cells in each direction 
				   (linked list) */   
  COORD_TYPE tol;               /* Tolerance of the shake algoritm used 
				   by RATTLE */

} OparamsMO;

#define MO_SAVE_LIST rxMO, ryMO, rzMO,\
                  vxMO, vyMO, vzMO,\
                  FxMO, FyMO, FzMO,\
                  vxo1MO, vyo1MO, vzo1MO,\
                  vxo2MO, vyo2MO, vzo2MO


#define MO_EXT_SLST  &sMO, &s1MO, &s2MO, &VolMO, &Vol1MO, &Vol2MO, &Vol1o1MO, &s1o1MO, &Vol1o2MO,&s1o2MO

#define MO_ALLOC_LIST  &rxMO, &ryMO, &rzMO,\
                    &vxMO, &vyMO, &vzMO,\
                    &FxMO, &FyMO, &FzMO,\
                    &vxo1MO, &vyo1MO, &vzo1MO,\
                    &vxo2MO, &vyo2MO, &vzo2MO

#define MO_DECL_LIST *rxMO, *ryMO, *rzMO,\
                    *vxMO, *vyMO, *vzMO,\
                    *FxMO, *FyMO, *FzMO,\
                    *vxo1MO, *vyo1MO, *vzo1MO,\
                    *vxo2MO, *vyo2MO, *vzo2MO
            				   
#define MO_EXT_DLST    sMO, s1MO, s2MO, VolMO, Vol1MO, Vol2MO, Vol1o1MO, s1o1MO, Vol1o2MO, s1o2MO

