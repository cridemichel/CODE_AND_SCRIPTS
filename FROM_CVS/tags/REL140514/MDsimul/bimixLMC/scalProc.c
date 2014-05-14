#include<mdsimul.h>
#include<scalProc.h>

/* GLOBAL VARIABLES */
extern char outFile[NAME_LENGTH], outFileBM[NAME_LENGTH];
/* output file (ascii file) name */
extern char inputFile[NAME_LENGTH];/* input file (measures file) */
extern char posFile[NAME_LENGTH];  /* name of the file with positions inside*/
extern int SEGSIZE, posBool, corBool, tBool, outBool, ihdr, ohdr, 
  infoBool, boxcmBool, bimixBool;

/* strings to store messages (use calling mdMsg()) */
extern COORD_TYPE T;
extern int DECL_LIST;
extern COORD_TYPE EXT_DLST;
double BM_EXT_DLST;
double BM_DECL_LIST; 

/* ============================ >>> help <<< ============================= */
void help(void)
{
  printf("Syntax: corutil [[-T: <temperature>] [-o <output file>]\n");
  printf(" [-Nm:<numb. of A atoms>,<numb. of B atoms>] [-nooh] [-noih] [-cor] [-infos] [-h]]\n");
  printf(" <coordinate file>\n");
}

/* =========================== >>> invalArg <<< ============================ */
void invalArg(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  help();
  printf("Invalid argument!\n");
  exit(-1);
}

/* =============================== >>> args <<< =============================*/
void args(int argc,char **argv)
{
  int i; /* fictitiuos counter */
  
  /* ==== set defaults =====*/
  
  T = 1.0; /* default temperature */

  strcpy(outFile, "scaled.cor"); 
  /* this is the defaukt output filename */

  /* =======================*/
  i = 0; 
  if (argc == 1) 
    invalArg();
  do /* scan all argument after the first one and before the last one */
    { 
      ++i;                            /* first i is 1, that is 
					 the second arg => OK */
      if (!strcmp(argv[i], "-h"))
	{
	  help();
	  exit(-1);
	}
      else if (!strcmp(argv[i], "-bimix:")) 
	{
	  bimixBool = 1;
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  
	  strcpy(outFileBM, argv[i]);

	}

      else if (!strcmp(argv[i], "-boxcm")) 
	{
	  boxcmBool = 1;
	}
      else if (!strcmp(argv[i], "-o")) 
	{
	  outBool = 1;
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  strcpy(outFile, argv[i]);
	} 
      else if (!strcmp(argv[i], "-Nm:"))
	{
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  sscanf(argv[i], "%d,%d", &Oparams.parnum[0], &Oparams.parnum[1]);
	  
	}

      else if (!strcmp(argv[i], "-pf:"))
	{
	  posBool = 1;
	  if ( ++i == (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg(); 
	    }
	  strcpy(posFile, argv[i]);
	}
      else if (!strcmp(argv[i], "-cor"))
	{
	  corBool = 1;
	}
      else if (!strcmp(argv[i], "-noih"))
	{
	  /* In this way the header is not read from the input file, that is 
	     the input file should be a raw coordinate file (raw = 
	     without header) */
	  ihdr = 0;
	}
      else if (!strcmp(argv[i], "-nooh"))
	{
	  /* omit the header in the output coordinate file */
	  ohdr = 0;
	}
      else if (!strcmp(argv[i], "-infos"))
	{
	  /* Print all the infos contained in the header of the file */
	  infoBool = 1;
	}
      else                            /* option not existant */ 
	{
	  invalArg();
	}
    } while ( i < (argc - 2) );
  /* and now take the input file */
  
  /* there is another arg ?  if N => ERROR */
  if (++i == argc)  
    {
      invalArg();
    }
 
  /* last arg that is input file name */
  strcpy(inputFile, argv[i]);
}

/* ========================= >>> AllocCoord <<< ========================= */
void allocCor(int size, int** pointer, ...)
{
  va_list ap;
  int** sptr;
  va_start(ap, pointer);
  *pointer = malloc(size);
  while((sptr = va_arg(ap, int**)) != NULL)
    {
      *sptr = malloc(size); 
    }
  va_end(ap);
}

/* ========================= >>> AllocCoord <<< ========================= */
void allocCorBM(int size, double** pointer, ...)
{
  va_list ap;
  double** sptr;
  va_start(ap, pointer);
  *pointer = malloc(size);
  while((sptr = va_arg(ap, double**)) != NULL)
    {
      *sptr = malloc(size); 
    }
  va_end(ap);
}

/* ======================== >>> loadSegs  <<< ===============================*/
void rSegs(int fdes, int size, void* pointer, ...)
{
  void* sptr;
  va_list ap;
  va_start(ap, pointer);
  
  /* Load first pointerof the list */
  if (read(fdes, pointer, size) < size) 
    {
      printf("ERROR: Too few bytes read!");
      exit(-1);
    }
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      if ( read(fdes, sptr, size) < size )
	{
	  printf("ERROR: Too few bytes read!");
	  exit(-1);
	}
    }
  va_end(ap);
}

/* =========================== >>> loadCor <<< ============================= */
void loadCor(int cfd)
{
  /* DESCRIPTION:
     Read structure params and coordinates from file whose file descriptor is 
     cfd.
     NOTE: nel salvataggio rispettare l'ordine di caricamento <----!!!!!
     Attualmente non viene effettuato nessun controllo sugli fread <----!!!!
   */
  if (ihdr == 1)
    {
      if (read(cfd, &Oparams,  sizeof(struct params)) < sizeof(struct params))
	{
	  printf("ERROR: Too few bytes read in the header!");
	  exit(-1);
	}
      //Nm = Oparams.parnum;
    }
  else if (Oparams.parnum[0] == 0 || Oparams.parnum[1] == 0)
    {
      printf("ERROR: You want to discard the header but you don't specify both particles numbers\n");
      scanf("Particles Number of type A? %d", &Oparams.parnum[0]);
      scanf("Particles Number of type B? %d", &Oparams.parnum[1]);
      printf("OK thanks\n");
    }
 
  /* allocate  coordinates needed by simulation (see. mdsimul_p) */
  SEGSIZE = sizeof(int) * Oparams.parnum[0];
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  allocCor(SEGSIZE, ALLOC_LISTA,
	   NULL);
  SEGSIZE = sizeof(int) * Oparams.parnum[1];
  allocCor(SEGSIZE, ALLOC_LISTB,
	   NULL);
  
  /* loads all arrays from the file associated with the fdes descriptor */
  SEGSIZE = sizeof(int) * Oparams.parnum[0];
  rSegs(cfd, SEGSIZE, SAVE_LISTA,
	NULL);    

  SEGSIZE = sizeof(int) * Oparams.parnum[1];
  rSegs(cfd, SEGSIZE, SAVE_LISTB,
	NULL);    

  /* NULL means: 'no more pointers to load' */
   
  /* loads extra variable such as the friction coefficents and its derivatives
     in the Hoover isotermal method */
  rSegs(cfd, sizeof(COORD_TYPE), EXT_SLST,
	NULL);    /* NULL means: 'no more pointers to load' */

}

/* ======================== >>> saveSegs  <<< ===============================*/
void wSegs(int fdes, int size, void* pointer, ...)
{
  void* sptr;
  va_list ap;
  va_start(ap, pointer);
  write(fdes, pointer, size); 
  /* at least one pointer should be present */
  
  /* if sptr = NULL => end of pointer list */ 
  while ( (sptr = va_arg(ap, void*)) != NULL )
    {
      /* if mode == EXIT if an error occurs then exit */
      write(fdes, sptr, size);
    }
  va_end(ap);
}


/* =========================== >>> saveCor <<< ============================*/
void saveCor(char* fileName)
{
  /* DESCRIPTION:
     save last coordinates on file fileName */
 
  int cfd; /* descriptor of coordinate file named fileName */ 
 
  cfd = creat(fileName, 0666);
  if (ohdr == 1)
    {
      write(cfd, &Oparams, sizeof(struct params));
    }
  /* writes all arrays to the disk physically making a sync() */
  //printf("SEGSIZE: %d\n", SEGSIZE);
  SEGSIZE = sizeof(int) * Oparams.parnum[0];
  wSegs(cfd, SEGSIZE, SAVE_LISTA,
	NULL);
  /* Write the extra variable (COORD_TYPE) */
  SEGSIZE = sizeof(int) * Oparams.parnum[1];
  wSegs(cfd, SEGSIZE, SAVE_LISTB,
	NULL);

  wSegs(cfd, sizeof(COORD_TYPE), EXT_SLST,
	NULL);

  close(cfd);

}

/* =========================== >>> setToZero <<< ========================== */
void setToZeroBM(int parnum, double* ptr, ...)
{
  /* DESCRIPTION:
     This procedure set to zero all the coordinates, passed as arguments */
  va_list ap;
  double* sptr;
  int i;
  
  va_start(ap, ptr);
  for(i = 0; i < parnum; ++i) ptr[i]=0.0;
  while ( (sptr = va_arg(ap, double*)) != NULL)
    {
       for(i = 0; i < parnum; ++i) sptr[i] = 0.0;
    }
  va_end(ap);
}

/* ============================ >>> ranf <<< =============================== */
COORD_TYPE ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return rand() / ( (COORD_TYPE) RAND_MAX );
}

/* ============================= >>> gauss <<< ============================= */
COORD_TYPE gauss(void)
{
  
  /* 
     Random variate from the standard normal distribution.
     
     The distribution is gaussian with zero mean and unit variance.
     REFERENCE:                                                    
                                                                
     Knuth D, The art of computer programming, (2nd edition        
     Addison-Wesley), 1978                                      
                                                                
     ROUTINE REFERENCED:                                           
                                                                
     COORD_TYPE ranf()                                  
     Returns a uniform random variate on the range zero to one  
  */

  COORD_TYPE  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  COORD_TYPE sum, r, r2;
  int i;

  sum = 0.0;

  loop(i, 1, 12)
    {
      sum = sum + ranf();
    }

  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}

/* ==================== >>> assignVelocities <<< =========================== */
void assignVelocities(double temp, double m[NA])
{
  COORD_TYPE rTemp[NA], sumx, sumy, sumz, RCMx, RCMy, RCMz;
  double MTOT;
  int i, a;
  
  rTemp[0] = sqrt(temp / Oparams.m[0]);
  rTemp[1] = sqrt(temp / Oparams.m[1]);
  /* variance of the velocities distribution function, we assume k = 1 */ 

  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vxBM[a][i] = rTemp[a] * gauss(); 
	  vyBM[a][i] = rTemp[a] * gauss();
	  vzBM[a][i] = rTemp[a] * gauss();
	  /* gauss() is a gaussian variable with mean = 0 and variance = 1, 
	     that is
	     2
	     1                X
        ----------- * exp( - --- )         
	 sqrt(2*PI)           2     */
	}
    }
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  MTOT = 0.0;
  for (a = 0;  a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* (sumx, sumy, sumz) is the total momentum */ 
	  sumx = sumx + m[a] * vxBM[a][i];
	  sumy = sumy + m[a] * vyBM[a][i];
	  sumz = sumz + m[a] * vzBM[a][i];
	  MTOT += m[a];
	}
    }
  sumx = sumx / MTOT; 
  sumy = sumy / MTOT;
  sumz = sumz / MTOT;

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  
  for (a = 0;  a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  vxBM[a][i] = vxBM[a][i] - sumx;
	  vyBM[a][i] = vyBM[a][i] - sumy;
	  vzBM[a][i] = vzBM[a][i] - sumz;
	  //printf("([%d,%d]: %.5f %.5f %.5f)\n",a , i, vxBM[a][i], vyBM[a][i],
	  // vzBM[a][i]);
	  /* In this way the total (net) momentum of the system of 
	     molecules is zero */
	}
    }
}

/* ======================= >>> saveCorBimix <<< ============================*/
void saveCorBimix(char* fileName)
{
  /* DESCRIPTION:
     save last coordinates on file fileName */
  
  /* per adesso il passo reticolare è costante ma poi deve diventare
     una parametro da riga di comando!*/
  double la = 0.0427272727272727;
  //const double la = 0.047;
  //const double la = 0.03133333333333;
  int a, b, i, SEGSIZE;
  
  int cfd; /* descriptor of coordinate file named fileName */ 
  
  /* Setta correttamente l'header file */
  OparamsBM.parnum[0] = Oparams.parnum[0]; 
  OparamsBM.parnum[1] = Oparams.parnum[1];
  OparamsBM.steplength = 0.01;
  OparamsBM.totStep = 100000;
  OparamsBM.curStep = 1;
  OparamsBM.d = 0;
  OparamsBM.P = Oparams.P;
  OparamsBM.T = Oparams.T;
  OparamsBM.m[0] = Oparams.m[0];
  OparamsBM.m[1] = Oparams.m[1];
  OparamsBM.rcut = Oparams.rcut;
  VolBM = Vol;
  sBM = 1.0;
  s2BM = Vol2BM = Vol1BM = s1BM = Vol1o1BM = Vol1o2BM = s1o1BM = s1o2BM = 0.0;
  la = cbrt(Vol) / Oparams.lattice_M; 
  printf("la: %.6f T: %.6f\n", la, Oparams.T);
  for (a = 0; a < NA; a++)
    {
      for (b = 0; b <NA; b++)
	{
	  OparamsBM.sigab[a][b] = Oparams.sigab[a][b];
	  OparamsBM.epsab[a][b] = Oparams.epsab[a][b];
	}
    }
  /* ===================== >>> FINE INIT <<<< ==================== */
  
  cfd = creat(fileName, 0666);
  if (ohdr == 1)
    {
      write(cfd, &OparamsBM, sizeof(struct paramsBM));
    }

  SEGSIZE = sizeof(double) * Oparams.parnum[0];
  /* SEGSIZE is the size in bytes of an array of coordinates */
  
  /* ALLOC_LIST is a macro defined in mdsimul.h and contains a list of 
     all addresses of the coordinates declared in the simulaiton
     (see that file) */
  allocCorBM(SEGSIZE, BM_ALLOC_LISTA,
	   NULL);
  SEGSIZE = sizeof(double) * Oparams.parnum[1];
  allocCorBM(SEGSIZE, BM_ALLOC_LISTB,
	   NULL);

  setToZeroBM(Oparams.parnum[0], BM_SAVE_LISTA, 
	    NULL);  /* Set to zero all the coordinates */
  
  setToZeroBM(Oparams.parnum[1], BM_SAVE_LISTB, 
	    NULL);  /* Set to zero all the coordinates */


  assignVelocities(Oparams.T, Oparams.m);

  /* Assegna le coordinate */
  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  rxBM[a][i] = la * ((double)rx[a][i]);
	  ryBM[a][i] = la * ((double)ry[a][i]);
	  rzBM[a][i] = la * ((double)rz[a][i]);
	  //printf("([%d,%d]: %.5f %.5f %.5f)\n",a , i, rxBM[a][i], ryBM[a][i],
	  // rzBM[a][i]);
		  
	}
    }

  /* writes all arrays to the disk physically making a sync() */
  //printf("SEGSIZE: %d\n", SEGSIZE);
  SEGSIZE = sizeof(double) * Oparams.parnum[0];
  wSegs(cfd, SEGSIZE, BM_SAVE_LISTA,
	NULL);

  /* Write the extra variable (COORD_TYPE) */
  SEGSIZE = sizeof(double) * Oparams.parnum[1];
  wSegs(cfd, SEGSIZE, BM_SAVE_LISTB,
	NULL);

  wSegs(cfd, sizeof(COORD_TYPE), BM_EXT_SLST,
	NULL);

  close(cfd);

}
