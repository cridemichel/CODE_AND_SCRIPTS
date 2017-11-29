#include<mdsimul.h>

/* GLOBAL VARIABLES */
char outFile[NAME_LENGTH];  /* output file (ascii file) name */
char outFileMO[NAME_LENGTH];
char inputFile[NAME_LENGTH];/* input file (measures file) */
char posFile[NAME_LENGTH];  /* name of the positions file */
int SEGSIZE, corBool = 0, tBool = 0, ihdr = 1, ohdr = 1, posBool = 0, Nm,
  boxcmBool = 0;
int outBool = 0; /* if true writes output coordinates file */
int infoBool = 0;
int monoBool = 0;

/* strings to store messages (use calling mdMsg()) */
COORD_TYPE T;
int DECL_LIST;
COORD_TYPE EXT_DLST;

/*  ======================= >>> Function prototypes <<< =====================*/
void args(int argc, char **argv);
void loadCor(int cfd);
void saveCor(char* fileName);
/* ========================================================================= */

/* ============================ >>> info <<< =============================== */
void info(void)
{
  printf("Informations about the coordinate file:\n");
  printf("Partcles number: %d\n",
	 Oparams.parnum);
  printf("External pressure P = %f\n", Oparams.P);
  printf("External temperature T = %f\n", Oparams.T);
  printf("mass of atom :%f\n", Oparams.m);
  printf("Total number of steps: %d steplength: %.10f\n", 
	 Oparams.totStep, Oparams.steplength);
  printf("Sigma:%f Epsilon;%f\n", Oparams.sigma, Oparams.epsilon);

  printf("curStep: %d Lattice_M:%d Lattice_a: %.6f\n", 
	 Oparams.curStep, Oparams.lattice_M, Oparams.lattice_a);
}

/* ============================== >>> main <<< ============================= */
void main(int argc, char* argv[])
{
  /* DESCRIPTION:
     This program scales the velocities, contained in the coordinate file
     specified as argument, with the option '-o' is possible to specify the
     name of the putput file.
     Beside it is possible to reduced all the coordinates to the first box
     and also to generate a raw coordinare file with the option '-nooh', or 
     to use a raw coordinate file (specifying the number of particles 
     with the option '-Nm: ')
     This programs should be recompile every time we use a different algorithm
     for the simulation because it works only on coordinates file of the 
     current simulation. */
  int cfd, i, Nm;
  COORD_TYPE K, Tcur, Rx, Ry, Rz, RCMx, RCMy, RCMz, Drx, Dry, Drz;
  COORD_TYPE m, scalFact;
  COORD_TYPE L, invL, Mtot; //r01;
  FILE* posf;
  double la;

  args(argc, argv);

  if ( (cfd = open(inputFile, O_RDONLY)) == -1)
    {
      printf("ERROR: Unable to open the coordinate file!\n");
      exit(-1);
    }
	       
  /* Read header and coordinates alocating the needed memory */
  printf("Loading coordinates...\n");
  loadCor(cfd);
  close(cfd);

  Nm = Oparams.parnum;

  L = cbrt(Vol);
  invL = 1.0 / L;
  /* Evaluate the instantenous temperature of the input coordiantes */
  K = 0.0;
  
  
  m = Oparams.m;
    
  

  /* coordinates rescaling */
  if ( corBool == 1)
    {
      printf("Reducing coordinates to first box...\n");
      for(i = 0; i < Oparams.parnum; i++)
	{
	  Rx =  rx[i];
	  Ry =  ry[i];
	  Rz =  rz[i];
	  Drx = - L * rint(Rx * invL);
	  Dry = - L * rint(Ry * invL);
	  Drz = - L * rint(Rz * invL);
	  
	  /* Scale coordinates to first cell */
	  rx[i] = rx[i] + Drx;
	  ry[i] = ry[i] + Dry;
	  rz[i] = rz[i] + Drz;
	}
	
    }

  /* ADD 27/1/1998:
     Put the center of mass of the box in the origin of axis */
  if (boxcmBool == 1)
    {
      Mtot = 0.0;
      RCMx = 0.0;
      RCMy = 0.0;
      RCMz = 0.0;
      
      for(i = 0;  i < Oparams.parnum; i++)
	{
	  RCMx += m * rx[i];
	  RCMy += m * ry[i];
	  RCMz += m * rz[i];
	  Mtot += m;
	}
	
      RCMx /= Mtot;
      RCMy /= Mtot;
      RCMz /= Mtot;
      printf("The actual box center of mass (%f, %f, %f) is now set to the origin\n", RCMx, RCMy, RCMz);
      
      
      for(i = 0; i < Oparams.parnum; i++)
	{
	  rx[i] -= RCMx;
	  ry[i] -= RCMy;
	  rz[i] -= RCMz;
	}
    
    }
  /* ============= >>> writing positions in an ascii form <<< ============== */
  if (posBool == 1)
    {
      la = cbrt(Vol) / Oparams.lattice_M; 
      printf("Writing positions...\n");
      posf = fopen(posFile, "w");
      for(i = 0; i < Oparams.parnum; i++)
	{
	  fprintf(posf, "%.8f %.8f %.8f\n", la*((double)rx[i]),
		  la*((double)ry[i]), 
		  la*((double)rz[i]));
	}
      
      fclose(posf);
    }

  /* ======================= >> saving coordinates <<< ====================*/
  if (outBool == 1)
    {
      /* And now save the new scaled coordinates in the output file */
      saveCor(outFile);
      printf("Coordinates file saved\n");
    }
  
  if (monoBool == 1)
    {
      /* And now save the new scaled coordinates in the output file */
      saveCorMono(outFileMO);
      printf("Coordinates file for mono simulation saved\n");
    }

  if (infoBool == 1)
    info();

}
