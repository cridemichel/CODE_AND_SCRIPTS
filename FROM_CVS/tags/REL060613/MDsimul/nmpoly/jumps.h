#define MP 500
COORD_TYPE xvadt, dt;
COORD_TYPE m0, m1, d, Vol, Mtot;
char inputFile[NAME_LENGTH], xvaparsFile[NAME_LENGTH], sqrtdr2File[255], 
tetaFile[255];
COORD_TYPE pi;

/* input file (measures file) */

char wp[1024];
char precision[64];
int Nm; /* Number of points for which calculate the time correlation
	   functions */
int tTot, printEvery, tBeg = 0, tgap, numParts, np[MP],
  DtTra_mean = 0, DtAng_mean = 0; 

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
  {"inputFile",  &inputFile,      STR},
  {"sqrtdr2File",&sqrtdr2File,    STR},
  {"DtTra_mean",    &DtTra_mean,  INT},
  {"DtAng_mean",    &DtAng_mean,  INT},
  {"tetaFile",   &tetaFile,       STR},
  {"wp",        &wp,              STR},
  {"m0",        &m0,              CT},
  {"m1",        &m1,              CT},
  {"d",         &d,               CT},
  {"dt",        &dt,              CT},
  {"Vol",       &Vol,             CT},
  {"tgap",      &tgap,            INT},
  {"prec",    precision,          STR},
  {"tBeg",      &tBeg,            INT},
  {"tTot",      &tTot,            INT},
  /* ======================================================================= */
  
  {"", NULL, 0} /* end of list, don't touch this !!! */
};
