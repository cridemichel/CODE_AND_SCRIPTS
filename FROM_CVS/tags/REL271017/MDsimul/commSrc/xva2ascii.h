#define MP 500
COORD_TYPE xvadt;
COORD_TYPE m, d, Vol, Mtot;
int nRun, PTM;

char fileLista[NAME_LENGTH],
  inputFile[NAME_LENGTH], xvaparsFile[NAME_LENGTH];
COORD_TYPE pi;

/* input file (measures file) */

char wp[1024];
char precision[64];
int Nm; /* Number of points for which calculate the time correlation
	   functions */
int tTot, printEvery, tBeg = 0, tgap, numParts, np[MP]; 
double Temp[512];
char tempStr[1024]; /* Stringa che contiene una lista di temperature separate
		       da virgole, ad es. 0.475, 0.485, ... */

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
  {"inputFile", inputFile,        STR},
  {"fileLista", fileLista,        STR},
  {"m",         &m,               CT},
  {"Vol",       &Vol,             CT},
  {"tgap",      &tgap,            INT},
  {"prec",      precision,        STR},
  {"tBeg",      &tBeg,            INT},
  {"PTM",       &PTM,             INT},
  {"n",         &nRun,            INT},
  {"T",         tempStr,          STR},
  /* ======================================================================= */
  
  {"", NULL, 0} /* end of list, don't touch this !!! */
};
