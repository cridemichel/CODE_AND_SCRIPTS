#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
double REul[3][3];
#define NUM_ATOMS 1319
#define NUM_LYSO 2
#define INIFILEGRO "inifile.gro"
/* it is possible to fix some selected degrees of freedom */
const int fix_r=0, fix_phi_r=1, fix_theta_r=1, fix_psi=1, fix_phi=1, fix_theta=1;
const double boxfact=2.0;
const double boxlength=50.0;
const double selfenergy=-20714.2;
double orig[3];
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
#define GROENEEXE GROPATH "g_energy"
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
char line[4096];
int resnum;
double xin[NUM_LYSO][NUM_ATOMS][3], xout[NUM_LYSO][NUM_ATOMS][3];
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
      com[kk] += xin[np][n][kk]; 
  
  for (kk=0; kk < 3; kk++)
    com[kk] /= NUM_ATOMS;
}

double com[2][3];
const double RMIN=6.5, RMAX=6.5, delI=0.05;
/* euler angles */
double phi, theta, psi;
double PI;
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
int main(int argc, char **argv)
{
  double xx[3], r, theta_r, phi_r, theta, psi, phi, deltheta, delpsi, 
	 delphi, delth_r, delphi_r;
  int kk;
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
#if 0
  /* just in to out */
  move_prot_copy(0);
  move_prot_copy(1);  
  calcCOM(1, com[1]);
  printf("DOPO com[1]=%.15G %.15G %.15G\n", com[1][0], com[1][1], com[1][2]);
  write_gro_coords();
  exit(-1);
#endif 
  //for (r=RMAX; r >= RMIN && !fix_r; r -= delI)
  for (r=RMIN; r <= RMAX && !fix_r; r += delI)
    {
      calc_pos(r, theta_r, phi_r, xx);
      DEBUG2(printf("r=%.15G (MAX=%.15G MIN=%.15G)\n", r, RMAX, RMIN));
#if 0
      for (theta_r = -PI*0.5; theta_r < PI*0.5 && !fix_theta_r; theta_r += delth_r)
	{
	  for (phi_r = 0.0; phi_r < 2.0*PI && !fix_phi_r; theta_r += delphi_r)
	    {
	      calc_pos(r, theta_r, phi_r, xx);
	      for (psi = -PI; psi < PI && !fix_psi; psi = psi + delpsi)
		{
		  for (phi = -PI; phi < PI && !fix_phi; phi = phi + delphi)
		    {
		      for (theta = 0.0; theta < PI && !fix_theta; theta + deltheta)
			{
#endif
			  /* just in to out */
			  move_prot_copy(0);
			  /* protein 0 is fixed while protein 1 moves */
			  move_prot(1, xx, psi, phi, theta);
			  /* calculate interaction energy between xout coordinates */
			  energy = calc_energy_gromacs();
			  /* write new element of potential energy mesh */
			  write_mesh_point_and_energy(xx, psi, phi, theta, energy);
#if 0
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
