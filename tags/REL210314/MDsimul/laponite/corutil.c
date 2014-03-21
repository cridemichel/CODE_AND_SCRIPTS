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
  printf("Partcles number: %d Atoms distance d = %f\n",
	 Oparams.parnum, Oparams.d);
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
  COORD_TYPE m[NA], scalFact, Mtot, Rxc, Ryc, Rzc, Rxa, Rya, Rza;
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
  if(Oparams.parnum != 0)
    {
      Nm = Oparams.parnum;
    }
  else 
    {
      Nm =0.0;
    }

  L = cbrt(Vol);
  invL = 1.0 / L;
  /* Evaluate the instantenous temperature of the input coordiantes */
  K = 0.0;
  
  loop(a, 1, NA)
    {
      m[a] = Oparams.m[a];
    }
  
  Mtot = m[0] + m[1];

  if (tBool == 1)
    {
      /* Scale velocity to obtain the desired temperature */
      loop(i, 1, Nm)
	{
	  loop(a, 1, NA)
	    {
	      K = K + 0.5 * m[a] * 
		( Sqr(vx[a][i]) + Sqr(vy[a][i]) + Sqr(vz[a][i]) );
	    }
	}
      Tcur = 2.0 * K / (5.0 * Nm - 3.0); 
      /* Now scale the velocity */
      printf("Scaling velocities...\n");
      scalFact = sqrt(T / Tcur);
      loop(i, 1,Nm)
	{
	  loop(a, 1, NA)
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
      loop(i, 1, Nm)
	{
	  
	  Rx = (m[0] * rx[0][i] + m[1] * rx[1][i]) / Mtot;
	  Ry = (m[0] * ry[0][i] + m[1] * ry[1][i]) / Mtot;
	  Rz = (m[0] * rz[0][i] + m[1] * rz[1][i]) / Mtot;
	  Drx = - L * rint(Rx * invL);
	  Dry = - L * rint(Ry * invL);
	  Drz = - L * rint(Rz * invL);
	  
	  loop(a, 1, NA)
	    {
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
      
      RCMx = 0.0;
      RCMy = 0.0;
      RCMz = 0.0;
      
      loop(i, 1, Nm)
	{
	  Rx = (m[0] * rx[0][i] + m[1] * rx[1][i]) / Mtot;
	  Ry = (m[0] * ry[0][i] + m[1] * ry[1][i]) / Mtot;
	  Rz = (m[0] * rz[0][i] + m[1] * rz[1][i]) / Mtot;
	  RCMx += Rx;
	  RCMy += Ry;
	  RCMz += Rz;
	}
      printf("The actual box center of mass (%f, %f, %f) is now set to the origin\n", RCMx, RCMy, RCMz);
      
      RCMx /= (COORD_TYPE) Nm;
      RCMy /= (COORD_TYPE) Nm;
      RCMz /= (COORD_TYPE) Nm;
      
      loop(i, 1, Nm)
	{
	  loop(a, 1, NA)
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
	loop(i, 1, Nm)
	{
	  //i=i+2;
	  Rxa=(rx[0][i]+rx[1][i])/2;
	  Rya=(ry[0][i]+ry[1][i])/2;
	  Rza=(rz[0][i]+rz[1][i])/2;
	  /*if ((Sqr(Rxa)+Sqr(Rya)+Sqr(Rza))<10.0)*/
	  /*if (fabs(Rxa)<3.0 && fabs(Rya)<3.0 && fabs(Rza)<3.0)*/
	  //{
	  fprintf(posf, "%.8f %.8f %.8f\n", rx[0][i], ry[0][i], 
		  rz[0][i]);
	  fprintf(posf, "%.8f %.8f %.8f\n", rx[1][i], ry[1][i], 
		  rz[1][i]);
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
