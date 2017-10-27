/* GLOBAL VARIABLES */
char wtdAngFile[128], wtdTraFile[128];
/* Corresponding flags */
int wtdTraFlag = 1, wtdAngFlag = 1;

int printEvery, tRun, DtTra_star, DtAng_star, DtTra_mean = 0, DtAng_mean = 0; 

COORD_TYPE m0, m1, d, Vol, DR_star, Dphi_star;
COORD_TYPE Mtot, xvadt, dt, dteta, T; /* T = temperature */

/* GLOBAL VARIABLES */
char inputFile[NAME_LENGTH], xvaparsFile[NAME_LENGTH];
/* input file (measures file) */

int* lastJump;
char precision[64];
int Nm; /* Number of points for which calculate the time correlation
	   functions */
int *wtdTra, *wtdAng;
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
  {"wtdAng",    wtdAngFile,       STR},
  {"DtAng_star",&DtAng_star,      INT},
  {"DtTra_star",&DtTra_star,      INT},
  {"Dphi_star", &Dphi_star,       CT},
  {"wtdTra",    wtdTraFile,       STR},
  {"DR_star",   &DR_star,         CT}, /* wtd gives the prob. that a part.
                                          moves not over DR_star */
  {"DtTra_mean",   &DtTra_mean,   INT},
  {"DtAng_mean",   &DtAng_mean,   INT},
  {"prec",    precision,          STR},
  {"m0",      &m0,                CT},
  {"m1",      &m1,                CT},
  {"d",       &d,                 CT},
  {"dt",      &dt,                CT},
  {"Vol",     &Vol,               CT},
  {"printEvery", &printEvery,     INT},
  /* ======================================================================= */
  
  {"", NULL, 0} /* end of list, don't touch this !!! */
};
