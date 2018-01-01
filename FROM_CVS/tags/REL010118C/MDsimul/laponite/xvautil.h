/* GLOBAL VARIABLES */
char velFile[128], psi1File[128], psi2File[128], C1File[128], C2File[128], 
  C3File[128], C4File[128], dphiSqFile[128], drSqFile[128], DtFile[128], 
  DrFile[128], vhFile[128], ddtdrFile[128], ddtdphiFile[128],
  AalphaFile[128], GselfFile[128], GsGsgaussFile[128], FselfFile[128],
  wtdFile[128];
/* Corresponding flags */
int velFlag=1, psi1Flag=1, psi2Flag=1, vhFlag=1, C1Flag=1, C2Flag=1, C3Flag=1,
  phiFlag=1, C4Flag=1, dphiSqFlag=1, drSqFlag=1, DtFlag=1, DrFlag=1, 
  ddtdrFlag = 1, ddtdphiFlag = 1, AalphaFlag = 1, GselfFlag = 1,
  GsGsgaussFlag = 1, FselfFlag = 1, wtdFlag = 1;

char phiFile[128]; /* File containing the angular positions at all instants */
int tCor, printEvery, nTeta, tBeg, Gsnr; 
/* The van Hove function is calculated between dteta and pi - dteta, that is:
   dteta < teta < (pi - dteta) */
				       
COORD_TYPE m0, m1, d, Vol, GsrMax, wtdRmax, kMax;
COORD_TYPE Mtot, xvadt, dt, dteta, T; /* T = temperature */
COORD_TYPE *phix, *phiy, *phiz, *phi0x, *phi0y, *phi0z,
  *phitt0x, *phitt0y, *phitt0z;
/* GLOBAL VARIABLES */
char inputFile[NAME_LENGTH], xvaparsFile[NAME_LENGTH];
/* input file (measures file) */

char precision[64];
int Nm; /* Number of points for which calculate the time correlation
	   functions */
int tgap; /* Increment of t0 */
int vhgap; /* steps every which save the van Hove function */
COORD_TYPE *velacf, *psi1acf, *psi2acf, *C1acf, *C2acf, *C3acf, 
  *C4acf, *wtd;
COORD_TYPE **vanHove;
int *normv, *norm;
COORD_TYPE *drSq, *dphiSq, *Dr, *Dt, *ddtdrSq, *ddtdphiSq, *Aalpha, *dr4,
  **Gself, **GsGsgauss, *Fself;
COORD_TYPE pi;

/* Array containing the valid parameter in the parameter file the xvautil */
struct singlePar OxvaPar[] = { 
  /* =================== >>> DON'T TOUCH THESE !!!! <<< =================== 
     The first string of each array elemen ( singlePas.parName ), is the
     string you can use in the parameters file to set its value.
     For example considering the definition below of the array me could in the
     parameters file the line:
        parnum: 100 
     It then put at begin of the simulation the value 10 in Opramas.parnum/
     Every parameter could be of this types:
     - INT = integer ( Linux-> 4 bytes signed integer )
     - STR = string NAME_LENGTH bytes long (see mdsimul.h)
     - CT = COORD_TYPE 
     WARNING: the params types should be consistent with previous params 
       structure declaration.
  */
  {"inputFile", &inputFile,       STR},
  {"tCor",    &tCor,              INT},
  {"vel" ,    velFile,            STR},
  {"psi1",    psi1File,           STR},
  {"psi2",    psi2File,           STR},
  {"vanHove", vhFile,             STR},
  {"vhgap",   &vhgap,             INT},
  {"nTeta",   &nTeta,             INT},
  {"dteta",   &dteta,              CT},
  {"C1",      C1File,             STR},
  {"C2",      C2File,             STR},
  {"C3",      C3File,             STR},
  {"C4",      C4File,             STR},
  {"Dt",      DtFile,             STR},
  {"Dr",      DrFile,             STR},
  {"drSq",    drSqFile,           STR},
  {"Aalpha",  AalphaFile,         STR},
  {"dphiSq",  dphiSqFile,         STR},
  {"ddtdrSq", ddtdrFile,          STR},
  {"ddtdphiSq", ddtdphiFile,      STR},
  {"Gself",   GselfFile,          STR},
  {"wtd",     wtdFile,            STR},
  {"wtdRmax", &wtdRmax,            CT}, /* wtd gives the prob. that a part.
                                          moves not over wdtRmax */
  {"GsGsgauss",GsGsgaussFile,    STR},
  {"Gsnr",    &Gsnr,              INT},
  {"GsrMax",  &GsrMax,            CT},
  {"Fself",   FselfFile,          STR},
  {"kMax",    &kMax,              CT},
  {"prec",    precision,          STR},
  {"m0",      &m0,                CT},
  {"m1",      &m1,                CT},
  {"d",       &d,                 CT},
  {"dt",      &dt,                CT},
  {"Vol",     &Vol,               CT},
  {"printEvery", &printEvery,     INT},
  {"tgap",       &tgap,           INT},
  {"tBeg",       &tBeg,           INT},
  {"T",          &T,              CT},
  /* ======================================================================= */
  
  {"", NULL, 0} /* end of list, don't touch this !!! */
};
