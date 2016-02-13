/* GLOBAL VARIABLES */
char heteroAngFile[128], heteroTraFile[128];
/* Corresponding flags */
int heteroTraFlag = 1, heteroAngFlag = 1, mathFlag = 0;
/* if mathFlag is true then produce output for mathematica */

int baseGrey = 20, maxGrey = 225;
int greyLevel = 20, printEvery, tRun, DtTra, DtAng, DtTra_mean = 0, DtAng_mean = 0; 

int t0;
COORD_TYPE m0, m1, d, Vol, mobTraPerc, mobAngPerc;
COORD_TYPE Mtot, xvadt, dt, dteta, T; /* T = temperature */

COORD_TYPE RADFACT = 4.0, *rCMt0x, *rCMt0y,
  *rCMt0z, *rCMt1x, *rCMt1y, *rCMt1z;


/* GLOBAL VARIABLES */
char inputFile[NAME_LENGTH], xvaparsFile[NAME_LENGTH];
/* input file (measures file) */

int* lastJump;
char precision[64];
int Nm; /* Number of points for which calculate the time correlation
	   functions */
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
  {"baseGrey", &baseGrey, INT},
  {"maxGrey",  &maxGrey,  INT},
  {"radFact",   &RADFACT,         CT},
  {"greyLevel", &greyLevel,       INT},
  {"inputFile", &inputFile,       STR},
  {"heteroTra",  heteroTraFile,       STR},
  {"heteroAng",  heteroAngFile,       STR},
  {"t0",        &t0,          INT},
  {"DtAng",&DtAng,      INT},
  {"DtTra",&DtTra,      INT},
  {"mobTraPerc",      &mobTraPerc,      CT},
  {"mobAngPerc",      &mobAngPerc,      CT},
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
