#include<mdsimul.h>

/* GLOBAL VARIABLES */
char outFile[NAME_LENGTH];  /* output file (ascii file) name */
char inputFile[NAME_LENGTH];/* input file (measures file) */
char posFile[NAME_LENGTH];  /* name of the positions file */
int SEGSIZE, corBool = 0, tBool = 0, ihdr = 1, ohdr = 1, posBool = 0, Nm,
  boxcmBool = 0;
int outBool = 0; /* if true writes output coordinates file */
int infoBool = 0;

/* strings to store messages (use calling mdMsg()) */
COORD_TYPE T;
COORD_TYPE DECL_LIST;
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
  printf("Nm = %d; P = %f; T= %f\n",
	 Oparams.parnum, Oparams.P, Oparams.T);
  printf("Total number of steps: %d\n", 
	 Oparams.totStep);
  printf("mass: %f sigma:%f epsilon:%f\n", Oparams.m, Oparams.sigma, Oparams.epsilon);
  printf("Linked Cell M = %d\n", Oparams.M);
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
  COORD_TYPE K, Tcur, RCMx, RCMy, RCMz, Drx, Dry, Drz;
  COORD_TYPE m, scalFact;
  COORD_TYPE L, invL; //r01;
  FILE* posf;
  
  args(argc, argv);
  Oparams.parnum = 0;

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
  m = Oparams.m;

  L = Oparams.L;
  invL = 1.0 / L;
  /* Evaluate the instantenous temperature of the input coordiantes */
  K = 0.0;
  
  if (tBool == 1)
    {
#if 0
      /* Scale velocity to obtain the desired temperature */
     for(i = 0; i < Nm; i++)
	{
	  K = K + 0.5 * m * 
	    ( Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]) );
	    
	}
     
     Tcur = 2.0 * K / (5.0 * Nm - 3.0); 
     /* Now scale the velocity */
     printf("Scaling velocities...\n");
     scalFact = sqrt(T / Tcur);
     for(i = 0; i< Nm; i++)
	{
	  vx[i] = vx[i] * scalFact;
	  vy[i] = vy[i] * scalFact;
	  vz[i] = vz[i] * scalFact;
	}
#endif
    }
  
  /* coordinates rescaling */
  if ( corBool == 1)
    {
      printf("Reducing coordinates to first box...\n");
      for(i = 0; i < Nm; i++)
	{
	  
	  Drx = - L * rint(rx[i] * invL);
	  Dry = - L * rint(ry[i] * invL);
	  Drz = - L * rint(rz[i] * invL);
	  
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
      
      RCMx = 0.0;
      RCMy = 0.0;
      RCMz = 0.0;
      
      for(i = 0; i < Nm; i++)
	{
	  RCMx += rx[i];
	  RCMy += ry[i];
	  RCMz += rz[i];
	}
      printf("The actual box center of mass (%f, %f, %f) is now set to the origin\n", RCMx, RCMy, RCMz);
      
      RCMx /= (COORD_TYPE) Nm;
      RCMy /= (COORD_TYPE) Nm;
      RCMz /= (COORD_TYPE) Nm;
      
      for(i = 0; i < Nm; i++)
	{
	      rx[i] -= RCMx;
	      ry[i] -= RCMy;
	      rz[i] -= RCMz;
	}
    }
  /* ============= >>> writing positions in an ascii form <<< ============== */
  if (posBool == 1)
    {
      printf("Writing positions...\n");
      posf = fopen(posFile, "w");
      for(i = 0; i < Nm; i++)
	{
	  /*if ((Sqr(Rxa)+Sqr(Rya)+Sqr(Rza))<10.0)*/
	  /*if (fabs(Rxa)<3.0 && fabs(Rya)<3.0 && fabs(Rza)<3.0)*/
	  //{
	  fprintf(posf, "%.8f %.8f %.8f\n", rx[i], ry[i], 
		  rz[i]);
	  //}
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
  
  if (infoBool == 1)
    info();

}
