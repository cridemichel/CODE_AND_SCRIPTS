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
  printf("Partcles number: (A:%d,B:%d\n",
	 Oparams.parnum[0], Oparams.parnum[1]);
  printf("External pressure P = %f\n", Oparams.P);
  printf("External temperature T = %f\n", Oparams.T);
  printf("mass of atom 0:%f of atom 1 %f\n", Oparams.m[0], Oparams.m[1]);
  printf("Total number of steps: %d steplength: %.10f\n", 
	 Oparams.totStep, Oparams.steplength);
  printf("Sigma00:%f Sigma11:%f Sigma10:%f\n", 
	 Oparams.sigab[0][0], Oparams.sigab[1][1], Oparams.sigab[1][0]);
  printf("Epsab00:%f Epsab11:%f Epsab10:%f\n", 
	 Oparams.epsab[0][0], Oparams.epsab[1][1], Oparams.epsab[1][0]);
  printf("Linked Cell M = %d\n", Oparams.M);
  printf("Curstep: %d\n", Oparams.curStep);
#ifdef SOFT_SPHERE
  printf("SOFT SPHERE PP=%d\n", Oparams.PP);
#endif
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
  int cfd, i, a, Nm;
  COORD_TYPE K, Tcur, Rx, Ry, Rz, RCMx, RCMy, RCMz, Drx, Dry, Drz;
  COORD_TYPE m[NA], scalFact;
  COORD_TYPE L, invL, Mtot; //r01;
  FILE* posf;
  
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

  Nm = Oparams.parnum[0] + Oparams.parnum[1];

  L = cbrt(Vol);
  invL = 1.0 / L;
  /* Evaluate the instantenous temperature of the input coordiantes */
  K = 0.0;
  
  loop(a, 1, NA)
    {
      m[a] = Oparams.m[a];
    }
  

  if (tBool == 1)
    {
      /* Scale velocity to obtain the desired temperature */
      for (a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a];  i++)
	    {
	      K = K + 0.5 * m[a] * 
		( Sqr(vx[a][i]) + Sqr(vy[a][i]) + Sqr(vz[a][i]) );
	    }
	}
      Tcur = 2.0 * K / (3.0 * Nm - 3.0); 
      /* Now scale the velocity */
      printf("Scaling velocities...\n");
      scalFact = sqrt(T / Tcur);
      for (a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      vx[a][i] = vx[a][i] * scalFact;
	      vy[a][i] = vy[a][i] * scalFact;
	      vz[a][i] = vz[a][i] * scalFact;
	    }
	}
    }
  
  /* coordinates rescaling */
  if ( corBool == 1)
    {
      printf("Reducing coordinates to first box...\n");
      for (a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      Rx =  rx[a][i];
	      Ry =  ry[a][i];
	      Rz =  rz[a][i];
	      Drx = - L * rint(Rx * invL);
	      Dry = - L * rint(Ry * invL);
	      Drz = - L * rint(Rz * invL);
	      
	      /* Scale coordinates to first cell */
	      rx[a][i] = rx[a][i] + Drx;
	      ry[a][i] = ry[a][i] + Dry;
	      rz[a][i] = rz[a][i] + Drz;
	    }
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
      
      for (a = 0; a < NA; a++)
	{
	  for(i = 0;  i < Oparams.parnum[a]; i++)
	    {
	      RCMx += m[a] * rx[a][i];
	      RCMy += m[a] * ry[a][i];
	      RCMz += m[a] * rz[a][i];
	      Mtot += m[a];
	    }
	}
      RCMx /= Mtot;
      RCMy /= Mtot;
      RCMz /= Mtot;
      printf("The actual box center of mass (%f, %f, %f) is now set to the origin\n", RCMx, RCMy, RCMz);
      

      for (a = 0;  a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      rx[a][i] -= RCMx;
	      ry[a][i] -= RCMy;
	      rz[a][i] -= RCMz;
	    }
	}
    }
  /* ============= >>> writing positions in an ascii form <<< ============== */
  if (posBool == 1)
    {
      printf("Writing positions...\n");
      posf = fopen(posFile, "w");
      for (a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      fprintf(posf, "%.8f %.8f %.8f\n", rx[a][i], ry[a][i], 
		      rz[a][i]);
	    }
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
