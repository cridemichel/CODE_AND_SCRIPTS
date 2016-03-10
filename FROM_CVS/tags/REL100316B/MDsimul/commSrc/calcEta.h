/* GLOBAL VARIABLES */
char Pinput[128], DQinput[128], PFile[128], EtaFile[128], 
  DQFile[128], ddtDQFile[128];
 
/* flags */
int PFlag = 1, EtaFlag = 1, DQFlag = 1, ddtDQFlag = 1, ddtNormFlag = 1,
  PNormFlag = 1;
 
int tCor, printEvery, tBeg, tRun; 
/* The van Hove function is calculated between dteta and pi - dteta, that is:
   dteta < teta < (pi - dteta) */
				       
COORD_TYPE Vol;
COORD_TYPE xvadt, dt, T; /* T = temperature */
COORD_TYPE *Pacf, *Eta, DQxy0, DQyz0, DQzx0, DQxytt0, DQyztt0, DQzxtt0, 
  Pxy0, Pyz0, Pzx0, Pxytt0, Pyztt0, Pzxtt0, *DQ, *ddtDQ,
  *Pxy, *Pyz, *Pzx, *DQxy, *DQyz, *DQzx;
 
/* GLOBAL VARIABLES */
char xvaparsFile[NAME_LENGTH];
/* input file (measures file) */

char precision[64];
int Nm; /* Number of points for which calculate the time correlation
	   functions */
int tgap; /* Increment of t0 */
int *normv, *norm;
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
  {"Pinput",  Pinput,             STR},
  {"DQinput", DQinput,            STR},
  {"DQFile",  DQFile,             STR},
  {"ddtDQFile", ddtDQFile,        STR},
  {"ddtNorm",   &ddtNormFlag,     INT},
  {"PNorm",    &PNormFlag,      INT},
  {"Pacf",    PFile,              STR},
  {"Eta",     EtaFile,            STR},
  {"tCor",    &tCor,              INT},
  {"prec",    precision,          STR},
  {"dt",      &dt,                CT},
  {"Vol",     &Vol,               CT},
  {"printEvery", &printEvery,     INT},
  {"tgap",       &tgap,           INT},
  {"tBeg",       &tBeg,           INT},
  {"T",          &T,              CT},
  /* ======================================================================= */
  
  {"", NULL, 0} /* end of list, don't touch this !!! */
};
