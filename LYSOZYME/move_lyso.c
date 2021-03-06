#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
double REul[3][3];
#define MD_CALC_ITENS
#define NUM_ATOMS 1319
#define NUM_LYSO 2
#define Sqr(x) ((x)*(x))
#define INIFILEGRO "inifile.gro"

/* it is possible to fix some selected degrees of freedom */
const int fix_r=0, fix_phi_r=0, fix_theta_r=0, fix_psi=0, fix_phi=0, fix_theta=0;
const double boxfact=2.0, kb = 0.00831451, Tamb=293.16, tempFact=500.0;
const double boxlength=50.0;
const double selfenergy=-20230;
double orig[3], totq, totmass; 
double xx_old[3], r_old, phi_old, theta_old, psi_old, energy_old, maxdelene, mindelene;
double r, theta_r, phi_r, theta, psi, phi, deltheta, delpsi, deltheta_r, 
       delphi, delth_r, delphi_r, phi_r_old, theta_r_old;
double min_del_r, min_deltheta_r, min_delphi_r, min_deltheta, min_delphi, min_delpsi;
#ifdef MD_CALC_ITENS
double ItensTot[3][3], Itens[NUM_LYSO][3], RotMat[NUM_LYSO][3][3];
#endif
double pdb_charge, pdb_mass;
int pdb_nr, pdb_resnr, pdb_cgnr, checkene, last_changed;
char pdb_type[256], pdb_residue[256], pdb_atom[256];
double com[2][3];
const double RMIN=3.3, RMAX=7.5;
double del_r=0.1;
/* euler angles */
double phi, theta, psi;
double PI;

#define DEBUG(x) x
#define DEBUG2(x) x
#if 1
#define DEVNULL " > /dev/null 2>&1"
#else
#define DEVNULL ""
#endif
/* ==== >>> GROMACS INTERFACE ==== <<< */  
#define MYHOME "/Users/demichel/"
//#define GROPATH MYHOME "TEMPO_DET/GROMACS/gromacs/bin/"
#define GROPATH "/usr/local/bin/"
#define GROCONV GROPATH "grompp_d"
#define GROEXE GROPATH "mdrun_d"
#define GROEDITCONF GROPATH "editconf_d"
#define GROGENBOX GROPATH "genbox_d"
#define MPIPATH "/usr/bin/"
#define MPIEXE MPIPATH "mpirun -np 2 "
#define GROENEEXE GROPATH "g_energy_d"
#define GROFILE "minimized_box.gro"
#define GROFILE0 "minimized.gro"
#define GROMOVEDCONF "moved.gro"
#define GROTOPOLOGY "3A3Q.top"
#define GROBINFILE "gromd.tpr"
#define GROMDPFILE "fullmd.mdp"
#define GROOUTPUTGRO "md_final.gro"
#define GROPOTENERGY "md_energy.edr"
#define GROTRAJTRR "md_traj.trr"
#define GROTRAJXTC "md_traj.xtc"
#define GROENEOUT "md_ener.xvg"
#define TMPENEFILE "_lastene_"
/* ===================================== */

char fin_name[4096], fout_name[4096];
//char foutgro_name[4096];
char line[4096], blocktype[4096], line2[4096];
int resnum;
double xin[NUM_LYSO][NUM_ATOMS][3], xout[NUM_LYSO][NUM_ATOMS][3];
double mass[NUM_ATOMS], charge[NUM_ATOMS];
char dummy_str1[4096], dummy_str2[4096];
char resname[10], atomname[10];
double dummy_dbl1, dummy_dbl2, Lx, Ly, Lz;
double energy;
void move_to_origin(int nprot, double oXYZ[3])
{
  int n, kk;
  for (n=0; n < NUM_ATOMS; n++)
    {
      for (kk=0; kk < 3; kk++)
	xin[nprot][n][kk] -= oXYZ[kk];
    }
}
void read_masses_and_charge(void)
{
  FILE *f;
  int atomsblk=0, i, ni;
  f = fopen(GROTOPOLOGY,"r");
  while (!feof(f))
   {
     fscanf(f,"%[^\n]\n", line);
     //printf("line=%s\n", line);
     if (!atomsblk)
       {
	 sscanf(line, " [ %s ] ", blocktype);
	 if (!strcmp(blocktype, "atoms"))
	   {
	     atomsblk = 1;
	     continue;
	   }
       }
     if (atomsblk)
       {
	 sscanf(line, "%[^;]\n", line2);

	 //printf("line2=%s\n", line2);
	 pdb_mass = 0.0;
	 pdb_charge = 0.0;
	 ni = sscanf(line2, " %d %s %d %s %s %d %lf %lf ", 
		    &pdb_nr, pdb_type, &pdb_resnr, pdb_residue, pdb_atom, &pdb_cgnr, &pdb_charge, &pdb_mass); 
	 //printf("pdb_nr=%d ni=%d\n", pdb_nr, ni);
	 if (!strcmp(line2,"") || ni < 8 )
	   continue;
	 if (ni == 8)
	   {
	     mass[pdb_nr-1] = pdb_mass; 
	   }
	 else 
	   {
	     printf("Missing Mass in Topology File (atom=%d, residue=%d)\n",pdb_nr, pdb_resnr);
	     exit(-1);
	   }
     	 if (ni >= 7)
	   {
	     charge[pdb_nr-1] = pdb_charge;
	   }
	 else
	   { 
	     printf("Missing Charge in Topology File (atom=%d, residue=%d)\n",pdb_nr, pdb_resnr);
	     exit(-1);
	   }
       }
     if (pdb_nr == NUM_ATOMS)
       break;
   }
  fclose(f);
  totq = 0.0;
  totmass = 0.0;
  for (i=0; i < NUM_ATOMS; i++ )
    {
      totq += charge[i];
      totmass += mass[i];
    }
  printf("TOTAL CHARGE=%8.5f TOTAL MASS=%8.5f\n", totq, totmass);
}
void read_gro_coords(void)
{
  FILE *fin;
  int i, n, kk;
  fin = fopen(fin_name, "r"); 
  if (fin==NULL)
    {
      perror("exiting...");
      exit(-1);
    }
  /* read first two lines */
  fscanf(fin, "%[^\n] ", line);
  //printf("line1=%s\n", line);
  fscanf(fin, "%[^\n] ", line);
  //printf("line2=%s\n", line);
  /* read coords */
  for (i = 0; i < NUM_LYSO; i++)
    {
      for (n = 0; n < NUM_ATOMS; n++)
	{
	  fscanf(fin, "%s %s %lf %lf %lf %lf\n", dummy_str1, dummy_str2, &dummy_dbl1, 
		 &xin[i][n][0], &xin[i][n][1], &xin[i][n][2]);  
	  //if (i==1)
	   // printf("xin=%f %f %f\n", xin[i][n][0], xin[i][n][1], xin[i][n][2]);
	  //printf("dummstr1=%s dummy2=%s dummy dbl=%G\n", dummy_str1, dummy_str2, dummy_dbl1);
	}
    }
  fscanf(fin, "%lf %lf %lf\n", &Lx, &Ly, &Lz);
  //printf("Lx=%.15G Ly=%.15G Lz=%.15G\n", Lx, Ly, Lz);
  //exit(-1);
  fclose(fin);
}
void write_gro_coords(void)
{
  FILE *fout, *fin;
  int i, n, kk;  
  fout = fopen(GROMOVEDCONF, "w");
  fin = fopen(fin_name, "r");
  /* read and write first two lines to output file */
  fscanf(fin, "%[^\n] ", line);
  //printf("prima linea=%s\n", line);
  fprintf(fout, "%s\n", line);
  fscanf(fin, "%[^\n] ", line);
  //printf("seconda linea=%s\n", line);
  fprintf(fout, "%s\n", line);
  for (i = 0; i < NUM_LYSO; i++)
    {
      for (n = 0; n < NUM_ATOMS; n++)
	{
	  fscanf(fin, "%d %s %s %[^\n] ", &resnum, resname, atomname, dummy_str1);
	  //printf("i=%d n=%d resnum=%d resname=%s atomname=%s dummy=%s\n", i, n, resnum, resname, atomname, dummy_str1);
	  fprintf(fout, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resnum, resname, atomname, 
		  i*NUM_ATOMS+n+1, xout[i][n][0], xout[i][n][1], xout[i][n][2]);
	  //printf("xout=%f %f %f\n", xout[i][n][0], xout[i][n][1], xout[i][n][2]);
	}
    }	  
  fprintf(fout, "%.15G %.15G %.15G\n", Lx*boxfact, Ly*boxfact, Lz*boxfact);
  fclose(fin);
  fclose(fout);
}

void build_euler_matrix(double phi, double theta, double psi, double Reul[3][3])
{
  /* build euler matrix */
  double cosphi, sinphi, costheta, sintheta, cospsi, sinpsi;

  /* phi e psi variano tra -pi e pi sono definite modulo 2*pi mentre theta varia fra 0 e pi 
     (27/01/10: CONTROLLARE SU LANDAU!)*/
  cosphi = cos(phi);
  sinphi = sin(phi);
  costheta = cos(theta);
  sintheta = sin(theta);
  cospsi = cos(psi);
  sinpsi = sin(psi);
  Reul[0][0] = cospsi*cosphi-costheta*sinphi*sinpsi;
  Reul[0][1] = cospsi*sinphi+costheta*cosphi*sinpsi;
  Reul[0][2] = sinpsi*sintheta;
  Reul[1][0] = -sinpsi*cosphi-costheta*sinphi*cospsi;
  Reul[1][1] = -sinpsi*sinphi+costheta*cosphi*cospsi;
  Reul[1][2] = cospsi*sintheta;
  Reul[2][0] = sintheta*sinphi;
  Reul[2][1] = -sintheta*cosphi;
  Reul[2][2] = costheta;
}
void move_prot_copy(int nprot)
{
  int n, kk;
  for (n = 0; n < NUM_ATOMS; n++)
    {
      for (kk=0; kk < 3; kk++)
	xout[nprot][n][kk] = xin[nprot][n][kk];
    } 
}
void move_prot(int nprot, double xx[3], double psi, double phi, double theta)
{
  /* xx is a displacement for protein nprot while psi, phi and theta are euler
     angles to set its orientation */
  int n, kk, jj;
  double xp[3], xl[3];

  build_euler_matrix(psi, phi, theta, REul);
  for (n = 0; n < NUM_ATOMS; n++)
    {
#if 1
      /* rotation */
      for (jj=0; jj < 3; jj++)
	{
	  xout[nprot][n][jj] = 0.0;
	  for (kk=0; kk < 3; kk++)
	    xout[nprot][n][jj]+=REul[jj][kk]*xin[nprot][n][kk];
	}	 
#else
      for(jj=0; jj < 3; jj++)
	{
	  xout[nprot][n][jj] = xin[nprot][n][jj];
	}	
#endif
      /* translation */
      for (kk=0; kk < 3; kk++)
	xout[nprot][n][kk] += xx[kk];	  
      //printf("xx=%f %f %f xout[%d]=%f %f %f\n", xx[0], xx[1], xx[2], n, xout[nprot][n][0],
	//     xout[nprot][n][1], xout[nprot][n][2]);
    }
    
}

void calcCOM(int np, double com[3])
{
  int i, n, kk;

  for (kk=0; kk < 3; kk++)
    com[kk] = 0.0;

  for (n=0; n < NUM_ATOMS; n++)
    for (kk=0; kk < 3; kk++)
      com[kk] += mass[n]*xin[np][n][kk]; 
  
  for (kk=0; kk < 3; kk++)
    com[kk] /= totmass;
}


void calc_pos(double r, double theta_r, double phi_r, double xx[3])
{
  xx[0] = r*cos(theta_r)*cos(phi_r);
  xx[1] = r*cos(theta_r)*sin(phi_r);
  xx[2] = r*sin(theta_r);  
}
char groconvstr[4096];
char grorunstr[4096];
char groenestr[4096];
char groecstr[4096];
void grosimulate(void)
{
  /* delete old files */
  system("rm -f " GROOUTPUTGRO);
  system("rm -f " GROPOTENERGY);
  system("rm -f " GROTRAJTRR);
  system("rm -f " GROTRAJXTC);
  system("rm -f *mdout*");
  system("rm -f *step*");
  system("rm -f *md.log*");
  system("rm -f *md_ener.xvg*");
  system("rm -f state*cpt");
  sprintf(grorunstr, GROEXE " -nice 0 -v -s " GROBINFILE " -o " GROTRAJTRR " -x " GROTRAJXTC " -c " GROOUTPUTGRO " -e " GROPOTENERGY DEVNULL);
  DEBUG(fprintf(stderr, "grorunstr: %s\n", grorunstr));
  system(grorunstr);
}

void ascii2bin_grocoords()
{
  system("rm -f " GROBINFILE);
  sprintf(groconvstr, GROCONV " -f " GROMDPFILE " -c " GROMOVEDCONF " -p " GROTOPOLOGY " -o "  GROBINFILE DEVNULL);
  DEBUG(fprintf(stderr,"groconvstr: %s\n", groconvstr));
  system(groconvstr);
}

void run_gromacs(void)
{
  ascii2bin_grocoords(); 
  grosimulate();
}

double read_energy(void)
{
  FILE *f;
  double ene;
  sprintf(groenestr, "echo ""10"" | " GROENEEXE " -dp -f " GROPOTENERGY " -o " GROENEOUT DEVNULL);
  system(groenestr);
  DEBUG(fprintf(stderr,"%s\n", groenestr));
  system("tail -1 " GROENEOUT " | awk '{print $2}' > " TMPENEFILE);
  f = fopen (TMPENEFILE, "r");
  if (fscanf(f,"%lf\n", &ene) < 1)
    {
      perror("boh, problem reading energy calculated by gromacs\n");
      exit(-1);
    }

  fclose(f);
  system("rm -f " TMPENEFILE);	 
  return ene;
}

double calc_energy_gromacs(void)
{
  write_gro_coords();
  run_gromacs();
  return read_energy()-selfenergy;
}
/* =============================== */
void write_mesh_point_and_energy(double xx[3], double psi, double phi, double theta, double energy)
{
  FILE *f;
  f = fopen(fout_name, "a+");
  fprintf(f, "%.15G %.15G %.15G %.15G %.15G %.15G %.15G\n", xx[0], xx[1], xx[2], psi, phi, theta, energy);
  fclose(f);
}
char delstr[4096];
#if 0
void gro_genbox(void)
{
  sprintf(groecstr, "%s -f %s -o %s -box %f %f %f -bt cubic", GROEDITCONF, fin_name, GROFILE, boxlength,
	  boxlength, boxlength);
  system(groecstr);
  sprintf(groecstr, "%s -maxsol 0 -cp %s -o %s -p %s", GROGENBOX, GROFILE, "out.gro", GROTOPOLOGY);
  system(groecstr);
  system("cp -f out.gro "  GROFILE);
}
#endif
#ifdef MD_CALC_ITENS
void calcItensTot(int np, double I[3][3], double RCM[3])
{
  int i, j, k;
  double distSq, ri[3], rj[3];
  double Icom[3][3];
  for (j=0; j < 3; j++)
    for (k=0; k < 3; k++)
      I[j][k] = 0.0;
  /* moment of inertia with respect to center of mass */
  printf("CoM= %f %f %f\n", RCM[0], RCM[1], RCM[2]);
  for (j=0; j < 3; j++)
    for (k=0; k < 3; k++)
      {
	I[j][k] = 0.0;
	for (i=0; i < NUM_ATOMS; i++)
	  {
	    ri[0] = xin[np][i][0]-RCM[0];
	    ri[1] = xin[np][i][1]-RCM[1];
	    ri[2] = xin[np][i][2]-RCM[2];
	    distSq = Sqr(ri[0])+Sqr(ri[1])+Sqr(ri[2]);
	    I[j][k] += mass[i]*(((j==k)?distSq:0.0) - ri[j]*ri[k]);
	    //printf("np=%d mass=%.15G x=%f %f %f\n", np, mass[i], ri[0], ri[1], ri[2]);
	    //printf("xin= %f %f %f\n", xin[np][i][0], xin[np][i][1], xin[np][i][2]);
	  }
      }
  printf("I= {{%f, %f, %f},\n", I[0][0], I[0][1], I[0][2]);
  printf("    {%f, %f, %f},\n", I[1][0], I[1][1], I[1][2]);
  printf("    {%f, %f, %f}}\n", I[2][0], I[2][1], I[2][2]);
}

/* find eigenvectors and eigenvalues */
void wrap_dsyev(double a[3][3], double b[3][3], double x[3], int *ok)
{
  char JOBZ, UPLO;
  double AT[9], work[10];
  int i, j, c1, c2, c3;
  JOBZ='V';
  UPLO='U';
  extern void dsyev_(char* , char*, int*, double* , int*, double*, double*, int*, int*);
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) AT[j+3*i]=a[j][i];		
    }						
  c1 = 3;
  c2 = 1;
  c3 = 8;
  dsyev_(&JOBZ, &UPLO, &c1, AT, &c1, x, work, &c3, ok);      
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) b[i][j]=AT[j+3*i];		
    }	
  if (*ok != 0)
    printf("not ok (%d)\n", *ok);
}

void calcEigenVectVal(double I[3][3], double Itens[3], double R[3][3])
{
  int ok;
  /* è la matrice con gli autovalori (=matrice di orientazione) */
  wrap_dsyev(I, R, Itens, &ok);
  if (ok!=0)  
    {
      printf("Eigenvalues calculated but some problems arose (ok=%d)\n", ok);
      exit(-1);
    }
}

void calcIprincAxes(int np)
{
  int kk;
  calcItensTot(np, ItensTot, com[np]); 
  calcEigenVectVal(ItensTot, Itens[np], RotMat[np]);
  /* il vettore riga RotMat[m][0..2] di RotMat è l'autovettore di ItensTot 
     con autovalore Itens[m]  
     quindi se x' sono le coordinate nel sistema di riferimento del 
     corpo rigido e x quelle nel laboratorio abbiamo:
     x' = R x ovvero x = Trasposta(R) x' */
  printf("R= {{%f, %f, %f},\n", RotMat[np][0][0], RotMat[np][0][1], RotMat[np][0][2]);
  printf("    {%f, %f, %f},\n", RotMat[np][1][0], RotMat[np][1][1], RotMat[np][1][2]);
  printf("    {%f, %f, %f}}\n", RotMat[np][2][0], RotMat[np][2][1], RotMat[np][2][2]);

  printf("Moments of Intertia = %.15G %.15G %.15G\n", Itens[np][0], Itens[np][1], Itens[np][2]);
}
#endif
void rotate_to_princaxes(int np)
{
  int i, m, n;
  double xx[3];
  for (i=0; i < NUM_ATOMS; i++) 
    {
      for (m=0; m < 3; m++)
	{
	  xx[m] = 0.0;
	  for (n=0; n < 3; n++)
	    {
	      xx[m] += RotMat[np][m][n]*xin[np][i][n];
	    }
	}
      for (m=0; m < 3; m++)
	xin[np][i][m] = xx[m];
    }
}
int check_min_step(int lc)
{
  switch (lc)
    {
    case 0:
      if (del_r < min_del_r)
	{
	  del_r = min_del_r;
	  r = r_old + del_r;
	}
      break;
    case 1:
      if (delphi_r < min_delphi_r)
	{
	  delphi_r = min_delphi_r;
	  phi_r = phi_r_old + delphi_r;
	}	  
      break;
    case 2:
      if (deltheta_r < min_deltheta_r)
	{
	  deltheta_r = min_deltheta_r;
  	  theta_r = theta_r_old + deltheta_r;
	}
      break;
    case 3:
      if (delpsi < min_delpsi)
	{
	  delpsi = min_delpsi;
	  psi = psi_old + delpsi;
	}
      break;
    case 4:
      if (delphi < min_delphi)
	{
	  delphi = min_delphi;
	  phi = phi_old + delphi;
	}
      break;
    case 5:
	if (deltheta < min_deltheta)
	  {
	    deltheta = min_deltheta;
	    theta = theta_old + deltheta;
	  }  
      break;
    }

}
/* NOTA 26/04/2010: routine che cerca di adattare il passo della mesh alla 
   ripidezza del potenziale */
void adapt_step(int lc, int todo) /* todo=1 increase todo=0 decrease */ 
{
  const double CGOLD = 1.618034;
  static double sf = 1.618034; 
  double GOLD = 1.618034;
  static int lasttodo=-1, lastlc=-1;
  /* 0 = r, 1= phi_r 2=theta_r, 3=psi, 4=phi, 5=theta */
  if (lastlc!=-1 && lasttodo!=-1)
    {
      if (lastlc==lc && lasttodo!=todo)
	if (todo==1)
	  sf /= CGOLD;
	else
	  sf *= CGOLD;
    }
  else
    sf = CGOLD;
  lasttodo = todo;
  lastlc = lc;
  if (todo==1)
    GOLD = 1.0/sf;
  else
    GOLD = sf;
  printf("delphi=%.15G lastlc=%d lc=%d lasttodo=%d todo=%d GOLD=%.15G\n", delphi, lastlc, lc, lasttodo, todo, GOLD);
  switch (lc)
    {
    case 0:
      del_r /= GOLD;
      r = r_old + del_r;
      break;
    case 1:
      delphi_r /= GOLD;
      phi_r = phi_r_old + delphi_r;
      break;
    case 2:
      deltheta_r /= GOLD;
      theta_r = theta_r_old + deltheta_r;
      break;
    case 3:
      delpsi /= GOLD;
      psi = psi_old + delpsi;
      break;
    case 4:
      delphi /= GOLD;
      phi = phi_old + delphi;
      break;
    case 5:
      deltheta /=GOLD;
      theta = theta_old + deltheta;  
      break;
    }
}
void reset_steps(void)
{
  del_r = 0.1;
  delphi_r = 0.3;/* rad (0-2*PI) */
  deltheta_r = 0.3;
  delphi = 0.3;
  deltheta = 0.3;
  delpsi = 0.3;
}
double r_ini, theta_r_ini, phi_r_ini, psi_ini, theta_ini, phi_ini;
int check_changed_dof(void)
{
  /* check which steps has been updated last time */
  /* 0 = r, 1= phi_r 2=theta_r, 3=psi, 4=phi, 5=theta */
  
  if (r_old != r) 
    return 0;
  if (phi_r_old != phi_r)
    return 1;
  if (theta_r_old != theta_r)
    return 2;
  if(psi_old != psi)
    return 3;
  if (phi_old != phi)
    return 4;
  if (theta_old != theta)
    return 5;
}
int first_check_ene;

int check_ene_func(void)
{
  if (fabs(energy) > tempFact*kb*Tamb)
    return 2;
  if (first_check_ene)
    {
      return 0;
    }
  else if (fabs(energy-energy_old) > maxdelene)
    return 1;
  else if (fabs(energy-energy_old) < mindelene)
    /* if energy change is too small reset steps */
    return 3;
  else if (fabs(energy-energy_old) > mindelene && fabs(energy-energy_old) < maxdelene)
    return 0;
  else
    return 0;
}
int check_first_ene(int lc)
{
  switch (lc)
    {
    case 0:
      if (r == r_ini)
	return 1;
      else
	return 0;
      break;
    case 1:
      if (phi_r == phi_r_ini)
	return 1;
      else
	return 0;
      break;
    case 2:
      if (theta_r == theta_r_ini)
	return 1;
      else
	return 0;
      break;
    case 3:
      if (psi == psi_ini)
	return 1;
      else
	return 0; 
      break;
    case 4:
      if (phi == phi_ini)
	return 1;
      else
	return 0; 
      break;
    case 5:
      if (theta == theta_ini)
	return 1;
      else
	return 0; 
      break;
    }
}
void set_inivalues(void)
{
  r_ini=RMAX;
  theta_r_ini=-PI*0.5;
  phi_r_ini=0.0; 
  psi_ini=-PI; 
  theta_ini=0.0; 
  phi_ini=-PI;
}
void set_min_steps(void)
{
  double F = 100.0;
  min_del_r = del_r / F;
  min_delphi_r = delphi_r / F;
  min_deltheta_r = deltheta_r / F;
  min_deltheta = deltheta / F;
  min_delpsi = delpsi / F;
  min_delphi = delphi / F;  
}
int main(int argc, char **argv)
{
  double xx[3]; 
  int kk, first=1;
  PI = 2.0*acos(0.0);
  if (argc == 1)
    {
      fprintf(stderr, "ERROR: You must supply gromacs input and mesh output files!\n");
      exit(-1);    
    } 
  strcpy(fin_name, argv[1]);
  strcpy(fout_name, argv[2]);
  //strcpy(foutgro_name, GROMOVEDCONF);
  /* delete mesh file if it exists already */
  read_masses_and_charge();
  sprintf(delstr, "rm -f %s",  fout_name);
  system(delstr);
  printf("input: %s output: %s\n", fin_name, fout_name);
  //gro_genbox();
  read_gro_coords();  
  /* calculate geometrical center of mass of the two lysozyme proteins */
  calcCOM(0, com[0]);
  /* move to origin, i.e. subtract center of mass of protein #0 */
  printf("com[0]=%.15G %.15G %.15G\n", com[0][0], com[0][1], com[0][2]);
  calcCOM(1, com[1]);
  printf("com[1]=%.15G %.15G %.15G\n", com[1][0], com[1][1], com[1][2]);
#ifdef MD_CALC_ITENS
  /* calculate moments of inertia around principal axes */
  calcIprincAxes(0);  
  calcIprincAxes(1);
#endif
  move_to_origin(0, com[0]);
#if 0
  for (kk=0; kk < 3; kk++)
    orig[kk] = com[1]-com[0];
#else
  for (kk=0; kk < 3; kk++)
    orig[kk] = com[1][kk];
#endif
  move_to_origin(1, orig);
  calcCOM(0, com[0]);
  /* move to origin, i.e. subtract center of mass of protein #0 */
  printf("DOPO com[0]=%.15G %.15G %.15G\n", com[0][0], com[0][1], com[0][2]);
  calcCOM(1, com[1]);
  printf("DOPO com[1]=%.15G %.15G %.15G\n", com[1][0], com[1][1], com[1][2]);

 
  theta_r = 0.0; 
  phi_r = 0.0;
  theta=psi=0.0;
  phi=0.0;
  rotate_to_princaxes(0);
  rotate_to_princaxes(1);
  DEBUG(printf("Inertia tensors now should be diagonal:\n"));
  DEBUG(calcItensTot(0, ItensTot, com[0])); 
  DEBUG(calcItensTot(1, ItensTot, com[0])); 
  DEBUG(printf("...are they?\n"));
#if 0
  /* just in to out */
  move_prot_copy(0);
  move_prot_copy(1);  
  calcCOM(1, com[1]);
  printf("DOPO com[1]=%.15G %.15G %.15G\n", com[1][0], com[1][1], com[1][2]);
  write_gro_coords();
  exit(-1);
#endif 
  first=1;
#if 0
  delphi = 3.0;
  deltheta = 3.0;
  delpsi = 3.0;
  delphi_r = 3.0;
  delth_r = 3.0;
#endif
 /* assign default values to steps */
  reset_steps();
  set_inivalues();
  set_min_steps();
  maxdelene = 10.0;/* in kJ/mol (unità di gromacs) */
  mindelene = 0.1;
  for (r=RMAX; r >= RMIN && !fix_r; r -= del_r)
  //for (r=RMIN; r <= RMAX && !fix_r; r += del_r)
    {
      //calc_pos(r, theta_r, phi_r, xx);
      DEBUG2(printf("r=%.15G (MAX=%.15G MIN=%.15G)\n", r, RMAX, RMIN));
#if 1
      for (theta_r = theta_r_ini; theta_r < PI*0.5 && !fix_theta_r; theta_r += deltheta_r)
	{
	  for (phi_r = phi_r_ini; phi_r < 2.0*PI && !fix_phi_r; phi_r += delphi_r)
	    {
	      calc_pos(r, theta_r, phi_r, xx);
	      for (psi = psi_ini; psi < PI && !fix_psi; psi += delpsi)
		{
		  for (phi = phi_ini; phi < PI && !fix_phi; phi += delphi)
		    {
		      for (theta = theta_ini; theta < PI && !fix_theta; theta += deltheta)
			{
#endif
			  /* 0 = r, 1= phi_r 2=theta_r, 3=psi, 4=phi, 5=theta */
			  last_changed = check_changed_dof();
			  first_check_ene = check_first_ene(last_changed);
			  //printf("first_check_ene:%d last_changed: %d\n", first_check_ene, last_changed);
			  /* just in to out */
			  do 
			    {
			      move_prot_copy(0);
			      /* protein 0 is fixed while protein 1 moves */
			      move_prot(1, xx, psi, phi, theta);
			      /* calculate interaction energy between xout coordinates */
			      /* store old values */
			      energy = calc_energy_gromacs();
			      /* se l'energia d'interazione è troppo alta non cercare di adattare il paso */
			      checkene = check_ene_func();
			      printf("checkene: %d\n", checkene);
			      if (checkene==2)
				{
				  //reset_steps();
				  break;
				}
			      else if (checkene==1||checkene==3)
				adapt_step(last_changed, (checkene==3)?1:0);
			      if (check_min_step(last_changed))
				break;
			    }
			  while (checkene);
			  for (kk=0; kk < 3; kk++)
			    xx_old[kk] = xx[kk];
			  psi_old = psi;
			  phi_old = phi;
			  theta_old = theta;
			  phi_r_old = phi_r;
			  theta_r_old = theta_r;
			  r_old = r;
			  energy_old = energy;
			  first = 0;
			  /* ================ */
			  /* write new element of potential energy mesh */
			  write_mesh_point_and_energy(xx, psi, phi, theta, energy);
#if 1
			}
		    }
		}
	    }
	}
#endif
  }

  //write_gro_coords();
  //system("");
  return 0;
}
