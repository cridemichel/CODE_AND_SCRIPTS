#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 1000
#define NA 6
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2
#define MD_SP_DELR 0.0
int MAXBONDS = 10;

#define Sqr(x) ((x)*(x))
char **fname; 
double *DR0[3], pi, time, *ti, *r0[3], *r1[3], L, refTime, *R[3][3], maxsax, maxax0, maxax1, maxsaxAA, 
       maxsaxAB, maxsaxBB, RCUT;
double sa[2]={-1.0,-1.0}, sb[2]={-1.0,-1.0}, sc[2]={-1.0,-1.0}, 
       Dr, theta, sigmaSticky, ratL[NA][3], *rat[NA][3], sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0;
int *inCell[3]={NULL,NULL,NULL}, *cellList=NULL, cellsx, cellsy, cellsz;
double *Fb[2], shift[3];
int *numbondst, **bondst, *numbonds0, **bonds0;
int points=-1, assez, NP, NPA, foundDRs=0;
int particles_type = 1;
char parname[128], parval[256000], line[256000];
char dummy[2048];
double A0, A1, B0, B1, C0, C1;
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
/* Allocate memory for a matrix of integers */
int** AllocMatI(int size1, int size2)
{
  int** v;
  int k;
  v = (int**) malloc(size1 * sizeof(int*));
  v[0] = (int*) malloc(size1 * size2 * sizeof(int));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}


void build_linked_list(void)
{
  double L2;
  int j, n;
  L2 = 0.5 * L;

  for (j = 0; j < cellsx*cellsy*cellsz + NP; j++)
    cellList[j] = -1;

  for (n = 0; n < NP; n++)
    {
      inCell[0][n] =  (rat[0][0][n] + L2) * cellsx / L;
      inCell[1][n] =  (rat[0][1][n] + L2) * cellsy / L;
      inCell[2][n] =  (rat[0][2][n] + L2) * cellsz / L;
      j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + NP;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}

double spApos[MD_STSPOTS_A][3] = {{MD_SP_DELR, 0.54, 0.0},{MD_SP_DELR, 0.54, 3.14159},{MD_SP_DELR, 2.60159,0.0},
    {MD_SP_DELR, 2.60159, 3.14159},{MD_SP_DELR, 1.5708, 0.0}};
double spBpos[MD_STSPOTS_B][3] = {{MD_SP_DELR, 0.0, 0.0},{MD_SP_DELR, 3.14159, 0.0}};

double spXYZ_A[MD_STSPOTS_A][3];
double spXYZ_B[MD_STSPOTS_B][3];

void build_atom_positions(void)
{
 /* N.B. le coordinate spXpos sono del tipo (Dr, theta, phi),
  * dove se Dr=0 la sfera sticky viene posizionata esattamente in 
  * maniera tangente e theta (0 <= theta <= Pi) e phi (0 <= phi < 2Pi)
  * sono gli angoli in coordinate sferiche che individuano il punto di contatto
  * tra sticky sphere ed ellissoide.
  * Tale routine converte le coordinate spXpos in coordinate cartesiane 
  * riferite al riferimento del corpo rigido. */  
  int kk, k1, aa;
  double x,y,z, grad[3], ng, dd[3];

  spApos[0][1] = theta;
  spApos[1][1] = theta;
  spApos[2][1] = pi - theta;
  spApos[3][1] = pi - theta;
  spApos[0][0] = Dr;
  spApos[1][0] = Dr;
  spApos[2][0] = Dr;
  spApos[3][0] = Dr;
  spApos[4][0] = Dr;
  spBpos[0][0] = Dr;
  spBpos[1][0] = Dr;
  //printf("Dr:%.15G pi=%.15G theta=%.15G\n", Dr, pi, theta);
  //printf("Dr: %.15G theta: %.15G pi=%.15G a=%.15G b=%.15G c=%.15G\n",
  //	 Dr, theta, pi, sa[1], sb[1], sc[1]);
  for (k1 = 0; k1 < MD_STSPOTS_A; k1++)
    {
      x = sa[0]*cos(spApos[k1][2])*sin(spApos[k1][1]);
      y = sb[0]*sin(spApos[k1][2])*sin(spApos[k1][1]);
      z = sc[0]*cos(spApos[k1][1]);
      //printf("xyz=%f %f %f\n", x, y, z);
      grad[0] = 2.0 * x / Sqr(sa[0]);
      grad[1] = 2.0 * y / Sqr(sb[0]);
      grad[2] = 2.0 * z / Sqr(sc[0]);
      ng = calc_norm(grad);
      for (aa = 0; aa < 3; aa++)
	grad[aa] /= ng;
      spXYZ_A[k1][0] = x + grad[0]*(sigmaSticky*0.5 + spApos[k1][0]);
      spXYZ_A[k1][1] = y + grad[1]*(sigmaSticky*0.5 + spApos[k1][0]);
      spXYZ_A[k1][2] = z + grad[2]*(sigmaSticky*0.5 + spApos[k1][0]);
	
	      //printf("k1=%d %f %f %f \n", k1,  spXYZ_A[k1][0] ,    spXYZ_A[k1][1] ,  spXYZ_A[k1][2]  );
    }
  for (kk=0; kk < 3; kk++)
    dd[kk] = spXYZ_A[0][kk] - spXYZ_A[1][kk];
  printf("Molecule A distance between Atoms 0 and 1: %.15G\n", calc_norm(dd));;
  for (kk=0; kk < 3; kk++)
    dd[kk] = spXYZ_A[2][kk] - spXYZ_A[3][kk];
  printf("Molecule A distance between Atoms 2 and 3: %.15G\n", calc_norm(dd));;

  for (k1 = 0; k1 < MD_STSPOTS_B; k1++)
    {
      x = sa[1]*cos(spBpos[k1][2])*sin(spBpos[k1][1]);
      y = sb[1]*sin(spBpos[k1][2])*sin(spBpos[k1][1]);
      z = sc[1]*cos(spBpos[k1][1]);
      grad[0] = 2.0 * x / Sqr(sa[1]);
      grad[1] = 2.0 * y / Sqr(sb[1]);
      grad[2] = 2.0 * z / Sqr(sc[1]);
      ng = calc_norm(grad);
      for (aa = 0; aa < 3; aa++)
	grad[aa] /= ng;
      spXYZ_B[k1][0] = x + grad[0]*(sigmaSticky*0.5 + spBpos[k1][0]);
      spXYZ_B[k1][1] = y + grad[1]*(sigmaSticky*0.5 + spBpos[k1][0]);
      spXYZ_B[k1][2] = z + grad[2]*(sigmaSticky*0.5 + spBpos[k1][0]) ;
    }
}
/* array con le posizioni degli atomi nel riferimento del corpo rigido 
 * nel caso dell'acqua i siti idrogeno ed elettroni sono disposti su 
 * di un tetraedro */
void BuildAtomPosAt(int i, int ata, double rO[3], double R[3][3], double rat[3])
{
  /* QUESTA VA RISCRITTA PER GLI ELLISSOIDI STICKY!!! */
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int k1, k2;
  double *spXYZ=NULL;
  //double radius; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */

  if (ata > 0)
    {
      if (i < NPA)
	spXYZ = spXYZ_A[ata-1];
      else  
	spXYZ = spXYZ_B[ata-1];
    }
  //radius = Oparams.sigma[0][1] / 2.0;
  if (ata == 0)
    {
      for (k1 = 0; k1 < 3; k1++)
	rat[k1] = rO[k1];
    }
  else 
    {
      for (k1 = 0; k1 < 3; k1++)
	{ 
	  rat[k1] = rO[k1];
	  for (k2 = 0; k2 < 3; k2++)
	    rat[k1] += R[k2][k1]*spXYZ[k2]; 
	}
    }
  
}

void BuildAtomPosAt32(int i, int ata, double rO[3], double R[3][3], double rat[3])
{
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  double r1[3], r2[3], r3[3], nr;
  double radius; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */

  /* NOTA: qui si assume che tutte e due le specie abbiamo lo stesso diametro!!! */
  radius = sigmaAA / 2.0;
  if (ata == 0)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk];
    }
  else if (ata <= 3)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] + R[kk][ata-1]; 
    }
  else
    {
      for (kk = 0; kk < 3; kk++)
	{
	  r1[kk] = R[kk][1]-R[kk][0];
	  r2[kk] = R[kk][2]-R[kk][0];
	}
      vectProdVec(r1, r2, r3);
      nr = calc_norm(r3);
      for (kk = 0; kk < 3; kk++)
	r3[kk] *= radius/nr;
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] - r3[kk]; 
    }
}

void BuildAtomPos32(int i, double rO[3], double R[3][3], double rat[5][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a, NUMAT;
  /* l'atomo zero si suppone nell'origine */
  if (i >= NPA)
    NUMAT = 4;
  else
    NUMAT = 3;
  for (a=0; a < NUMAT; a++)
    BuildAtomPosAt32(i, a, rO, R, rat[a]);
}

void BuildAtomPos(int i, double rO[3], double R[3][3], double rat[NA][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a;
  /* l'atomo zero si suppone nell'origine */
  if (i < NPA)
    {
      for (a=0; a < MD_STSPOTS_A+1; a++)
	BuildAtomPosAt(i, a, rO, R, rat[a]);
    }
  else
    {
      for (a=0; a < MD_STSPOTS_B+1; a++)
	BuildAtomPosAt(i, a, rO, R, rat[a]);
    }
}
#define Sqr(x) ((x)*(x))
int check_distance(int i, int j, double Dx, double Dy, double Dz)
{
  double DxL, DyL, DzL;
  double ma=0.0;
  DxL = fabs(Dx);
  DyL = fabs(Dy);
  DzL = fabs(Dz);

  if (particles_type == 1)
    {
      ma = maxsax;
    }
  else if (particles_type == 0)
    {
      if (i < NPA && j < NPA)
	ma = maxsaxAA;
      else if (i >= NPA && j >= NPA)
	ma = maxsaxBB;
      else
	ma = maxsaxAB;
    }
  if (DxL > ma || DyL > ma || DzL > ma)
    return 1;
  else 
    return 0;
}


double distance(int i, int j)
{
  int a, b;
  int maxa=0, maxb=0;
  double imgx, imgy, imgz;
  double Dx, Dy, Dz;

  Dx = rat[0][0][i] - rat[0][0][j];
  Dy = rat[0][1][i] - rat[0][1][j];
  Dz = rat[0][2][i] - rat[0][2][j];
  imgx = -L*rint(Dx/L);
  imgy = -L*rint(Dy/L);
  imgz = -L*rint(Dz/L);

  if (check_distance(i, j, Dx+imgx, Dy+imgy, Dz+imgz))
    return 1;

  if (particles_type == 1)
    {
      maxa = MD_STSPOTS_A;
      maxb = MD_STSPOTS_B;
    }
  else if (particles_type == 0)
    {
      if (i < NPA && j >= NPA)
	{
	  maxa = 2;
	  maxb = 3;
	}
      else if (i < NPA && j < NPA)
	{
	  maxa = 2;
	  maxb = 2;
	}
      else if (i >= NPA && j >= NPA)
	{
	  maxa = 3; 
	  maxb = 3;
	}
      else
	{
	  maxa = 3;
	  maxb = 2;
	}
    }
  for (a = 1; a < maxa+1; a++)
    {
      for (b = 1; b < maxb+1; b++)
	{

	  //printf("dist=%.14G\n", sqrt( Sqr(rat[a][0][i] + img*L -rat[a][0][j])+Sqr(rat[a][1][i] + img*L -rat[a][1][j])
	    //  +Sqr(rat[a][2][i] + img*L -rat[a][2][j])));
	  //printf("[DISTANCE] (%d,%d)-(%d,%d) dist=%.14G\n", i, a, j, b, sqrt( Sqr(rat[a][0][i] + imgx -rat[b][0][j])+Sqr(rat[a][1][i] + imgy -rat[b][1][j]) +Sqr(rat[a][2][i] + imgz -rat[b][2][j])));
	 
	  if (Sqr(rat[a][0][i] + imgx -rat[b][0][j])+Sqr(rat[a][1][i] + imgy -rat[b][1][j])
	      +Sqr(rat[a][2][i] + imgz -rat[b][2][j]) < Sqr(sigmaSticky))	  
		return -1;
	}
    }
  return 1;
}
#if 0
int check_distanceR(int i, int j, int imgix, int imgiy, int imgiz,
		  int imgjx, int imgjy, int imgjz, double imgx, double imgy, double imgz)
{
  double Dx, Dy, Dz;

  Dx = fabs(rat[0][0][i] + (imgix-imgjx)*L + imgx  - rat[0][0][j]);
  Dy = fabs(rat[0][1][i] + (imgiy-imgjy)*L + imgy - rat[0][1][j]);
  Dz = fabs(rat[0][2][i] + (imgiz-imgjz)*L + imgz - rat[0][2][j]);

  if (Dx > maxsax || Dy > maxsax || Dz > maxsax)
    {
      return 1;
    }
  else 
    {
      return 0;
    }
}
#endif
double distanceR(int i, int j, int imgix, int imgiy, int imgiz,
		  int imgjx, int imgjy, int imgjz, double Lbig)
{
  int a, b, maxa=0, maxb=0;
  double imgx, imgy, imgz;
  double Dx, Dy, Dz, dx, dy, dz;
  if (particles_type == 1)
    {
      if (i < NPA)
	{
	  maxa = MD_STSPOTS_A;
	  maxb = MD_STSPOTS_B;
	}
      else
	{
	  maxa = MD_STSPOTS_B;
	  maxb = MD_STSPOTS_A;
	}
    }
  else if (particles_type == 0)
    {
      if (i < NPA && j >= NPA)
	{
	  maxa = 2;
	  maxb = 3;
	}
      else if (i < NPA && j < NPA)
	{
	  maxa = 2;
	  maxb = 2;
	}
      else if (i >= NPA && j >= NPA)
	{
	  maxa = 3; 
	  maxb = 3;
	}
      else
	{
	  maxa = 3;
	  maxb = 2;
	}
    }

  dx = L*(imgix-imgjx);
  dy = L*(imgiy-imgjy);
  dz = L*(imgiz-imgjz);
  Dx = rat[0][0][i] - rat[0][0][j] + dx;
  Dy = rat[0][1][i] - rat[0][1][j] + dy;
  Dz = rat[0][2][i] - rat[0][2][j] + dz;
  imgx = -Lbig*rint(Dx/Lbig);
  imgy = -Lbig*rint(Dy/Lbig);
  imgz = -Lbig*rint(Dz/Lbig);

  if (check_distance(i, j, Dx + imgx, Dy + imgy, Dz + imgz))
    return 1;
  //printf("i=%d j=%d rat[][0][i]=%.15G,%.15G\n", i, j, rat[1][0][i], rat[2][0][i]);
  for (a = 1; a < maxa+1; a++)
    {
      for (b = 1; b < maxb+1; b++)
	{
	  //printf("[DISTANCER] (%d,%d)-(%d,%d) dist=%.14G\n", i, a, j, b, sqrt( Sqr(rat[a][0][i] + dx + imgx -rat[b][0][j])+Sqr(rat[a][1][i] + dy + imgy -rat[b][1][j]) +Sqr(rat[a][2][i] + dz + imgz -rat[b][2][j])));
	  if (Sqr(rat[a][0][i] + dx + imgx - rat[b][0][j])+Sqr(rat[a][1][i] + dy + imgy - rat[b][1][j])
	      + Sqr(rat[a][2][i] + dz + imgz - rat[b][2][j]) < Sqr(sigmaSticky))	  
		return -1;
	}
    }
  return 1;
}


int bond_found(int i, int j)
{
  if (distance(i, j) < 0.0)
    return 1;
  else
    return 0;
}
int bond_foundR(int i, int j, int imgix, int imgiy, int imgiz,
		int imgjx, int imgjy, int imgjz, double Lbig)
{
  if (distanceR(i, j, imgix, imgiy, imgiz, imgjx, imgjy, imgjz, Lbig) < 0.0)
    return 1;
  else
    return 0;
}

void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *DR[3], double *R[3][3])
{
  FILE *f;
  int nat=0, i, cpos;
  f = fopen(fname, "r");
  while (!feof(f) && nat < 2) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	}
      if (nat < 2)
	{
	  fseek(f, cpos, SEEK_SET);
	  fscanf(f, "%[^:]:", parname);
	  //printf("[%s] parname=%s\n", fname, parname);
	  if (!strcmp(parname,"DR"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf %lf %lf ", &DR[0][i], &DR[1][i], &DR[2][i]);
		}
	      foundDRs = 1;
	    }
#if 0
	  else if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	      foundrot = 1;
	    }
	  else if (!strcmp(parname,"sumoy"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[1][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoz"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[2][i]); 
		}
	    }
#endif
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else if (!strcmp(parname, "refTime"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *refTime = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else
	{
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%[^\n]\n", line); 
	      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %[^\n]\n", 
    		     &r[0][i], &r[1][i], &r[2][i], 
		     &R[0][0][i], &R[0][1][i], &R[0][2][i], &R[1][0][i], &R[1][1][i], &R[1][2][i],
		     &R[2][0][i], &R[2][1][i], &R[2][2][i], dummy); 
	      //printf("%.15G %.15G %.15G\n", R[2][0][i],R[2][1][i], R[2][2][i] );
	    
	    }
	  break; 
	}

    }
  fclose(f);
}

#define KMODMAX 600
#define NKSHELL 150
double qx[KMODMAX][NKSHELL], qy[KMODMAX][NKSHELL], qz[KMODMAX][NKSHELL];
double *sqReA[KMODMAX], *sqImA[KMODMAX], *sqReB[KMODMAX], *sqImB[KMODMAX];
double *cc;
char fname2[512];
char inputfile[1024];
double twopi;
void print_usage(void)
{
  printf("calcbondcorr [ --help/-h ] <lista_files> [points]\n");
  printf("where points is the number of points of the correlation function\n");
  exit(0);
}
double qavg[KMODMAX];

void parse_param(int argc, char** argv)
{
  int cc=1, extraparam=0;
  
  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      //printf("cc=%d extraparam=%d argc=%d\n", cc, extraparam, argc);
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (cc == argc || extraparam == 2)
	print_usage();
      else if (extraparam == 0)
	{ 
	  extraparam++;
	  //printf("qui1 extraparam:%d\n", extraparam);
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam == 1)
	{
	  extraparam++;
	  //printf("qui2 argv[%d]:%s\n",cc, argv[cc]);
	  points = atoi(argv[cc]);
	  //printf("points:%d\n", points);
	}
      else
	print_usage();
      cc++;
    }
}
void add_bond(int i, int j, int *numbonds, int **bonds)
{
  if (numbonds[i] > MAXBONDS)
    {
      printf("Too many bonds?!?\n");
      exit(-1);
    }
  bonds[i][numbonds[i]] = j;
  numbonds[i]++;
}

void find_bonds(int *numbonds, int **bonds, double *eneO)
{
  int i, iZ, jZ, iX, jX, iY, jY, j;
  double shift[3];
  double ene=0.0;
  for (i = 0; i < NP; i++)
    {
      for (iZ = -1; iZ <= +1; iZ++) 
	{
	  jZ = inCell[2][i] + iZ;    
	  shift[2] = 0.;
	  /* apply periodico boundary condition along z if gravitational
	   * fiels is not present */
	  if (jZ == -1) 
	    {
	      jZ = cellsz - 1;    
	      shift[2] = - L;
	    } 
	  else if (jZ == cellsz) 
	    {
	      jZ = 0;    
	      shift[2] = L;
	    }
	  for (iY = -1; iY <= +1; iY ++) 
	    {
	      jY = inCell[1][i] + iY;    
	      shift[1] = 0.0;
	      if (jY == -1) 
		{
		  jY = cellsy - 1;    
		  shift[1] = -L;
		} 
	      else if (jY == cellsy) 
		{
		  jY = 0;    
		  shift[1] = L;
		}
	      for (iX = -1; iX <= +1; iX ++) 
		{
		  jX = inCell[0][i] + iX;    
		  shift[0] = 0.0;
		  if (jX == -1) 
		    {
		      jX = cellsx - 1;    
		      shift[0] = - L;
		    } 
		  else if (jX == cellsx) 
		    {
		      jX = 0;   
		      shift[0] = L;
		    }
		  j = (jZ *cellsy + jY) * cellsx + jX + NP;
		  for (j = cellList[j]; j > -1; j = cellList[j]) 
		    {
		      switch (particles_type)
			{
			case 0:
			  if (j <= i) 
			    continue;
			  break;
			case 1:
			  if ((i < NPA && j < NPA) || ( i >= NPA && j >= NPA) ||
			      (i >= NPA && j < NPA))
			    continue;
			  break;
			}
		      if (bond_found(i, j))  
			{
			  //printf("qui1\n");
		      	  add_bond(i, j, numbonds, bonds);
	    		  add_bond(j, i, numbonds, bonds);
			  //printf("qui2\n");
			  ene=ene+1.0;
			}
		    }
		}
	    }
	}
    }
  *eneO = ene;
}
int exist_bond(int na, int n, int a, int b, int *numbonds, int **bonds)
{
  int i;
  for (i = 0; i < numbonds[na]; i++)
    if (bonds[na][i] == n*(NA*NA)+a*NA+b)
      return 1;
  return 0;
}
void build_spots_positions(double *r0[3], double *R[3][3])
{
  int i, a, b;
  double r0L[3], RL[3][3], ratL[NA][3];
  for (i = 0; i < NP; i++)
    {
      /* qui va il codice per individuare i cluster */
      for (a = 0; a < 3; a++)
	{
	  r0L[a] = r0[a][i];
	  for (b = 0; b < 3; b++)
	    RL[a][b] = R[a][b][i];
	}
      //printf("r0L[%d]=%.15G %.15G %.15G\n", i, r0L[0], r0L[1], r0L[2]);
      if (particles_type == 1)
	BuildAtomPos(i, r0L, RL, ratL);
      else if (particles_type == 0)
	BuildAtomPos32(i, r0L, RL, ratL);
      for (a = 0; a < NA; a++)
	for (b = 0; b < 3; b++)
	  rat[a][b][i] = ratL[a][b];
      //printf("rat[]=%.15G %.15G %.15G\n", rat[0][0][i], rat[0][1][i], rat[0][2][i]);
    }
}
int main(int argc, char **argv)
{
  FILE *f, *f2;
  int first=1, firstp=1, c1, c2, c3, i, ii, nr1, nr2, a, b;
  int iq, NN, fine, JJ, maxl, nfiles, nat, np, maxnp, jj2;
  int nb, j; 
  double invL, rxdummy, sumImA, sumReA, sumImB, sumReB, scalFact, ene0, enet;

  twopi = acos(0)*4.0;	  
  pi = twopi / 2.0;	
#if 0
  if (argc <= 1)
    {
      printf("Usage: calcfqself <lista_file> [points] [qmin] [qmax]\n");
      printf("where points is the number of points of the correlation function\n");
      exit(-1);
    }
#endif 
  parse_param(argc, argv);
  //printf(">>qmin=%d\n", qmin);
  c2 = 0;
  if (!(f2 = fopen(inputfile, "r")))
    {
      printf("ERROR: I can not open file %s\n", inputfile);
      exit(-1);
    }
  maxl = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", dummy); 
      if (strlen(dummy)+1 > maxl)
	maxl = strlen(dummy)+1;
      c2++;
    }	
  nfiles = c2;
  rewind(f2);
  fname = malloc(sizeof(char*)*nfiles);
  for (ii=0; ii < nfiles; ii++)
    {
      fname[ii] = malloc(sizeof(char)*maxl);
      fscanf(f2, "%[^\n]\n", fname[ii]); 
    }

  if (!(f = fopen(fname[0], "r")))
    {
      printf("ERROR: I can not open file %s\n", fname[0]);
      exit(-1);
    }
  nat = 0;
  while (!feof(f) && nat < 2) 
    {
      fscanf(f, "%[^\n]\n)", line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  if (nat==2)
	    {
	      for (i=0; i < 2*NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      fscanf(f, "%lf\n", &L);
	      break;
	    }
	  continue;
	}
      sscanf(line, "%[^:]:%[^\n]\n", parname, parval); 
      if (!strcmp(parname,"parnum"))
	NP = atoi(parval);
      else if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      else if (!strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (!strcmp(parname, "a"))
       	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &A0, &A1);
	}
      else if (!strcmp(parname, "b"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &B0, &B1);
	}
      else if (!strcmp(parname, "c"))
	{
	  fscanf(f, "%[^\n]\n", parval);
	  sscanf(parval, "%lf %lf ", &C0, &C1);
	}
      else if (nat==1 && !strcmp(parname,"sigma"))
	sscanf(parval, "%lf %lf %lf %lf\n", &sigmaAA, &sigmaAB, &sigmaAB, &sigmaBB);	
      else if (nat==1 && !strcmp(parname,"sigmaSticky"))
	sigmaSticky = atof(parval);
      else if (nat==1 && !strcmp(parname,"theta"))
	theta = atof(parval);
      else if (nat==1 && !strcmp(parname,"Dr"))
	Dr = atof(parval);

    }
  fclose(f);
  invL = 1.0/L;
#if 0
  if (argc >= 3)
    points = atoi(argv[2]);
  else
    points = NN;
  if (argc >= 4)
    qmin = atoi(argv[3]);
  if (argc == 5)
    qmax = atoi(argv[4]);
#endif
  if (points == -1)
    points = NN;
  scalFact = twopi * invL;
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;

  printf("invL=%.15G\n", invL);
  //printf("maxnp=%d points=%d\n",maxnp, points);
  if ((A0 > B0 && A0 > C0) || (A0 < B0 && A0 < C0))
    assez = 0;
  else if ((B0 > A0 && B0 > C0) || (B0 < A0 && B0 < C0))
    assez = 1;
  else if ((C0 > A0 && C0 > B0) || (C0 < A0 && C0 < B0))
    assez = 2;
  if (sigmaAA != -1.0)
    particles_type = 0;

  if (NPA == -1)
    NPA = NP;
  //fprintf(stderr, "allocating %d items NN=%d NP=%d num files=%d maxnp=%d\n", points, NN, NP, nfiles, maxnp);
  cc = malloc(sizeof(double)*points);
  ti = malloc(sizeof(double)*points);
  numbonds0  = malloc(sizeof(int)*NP); 
  bonds0 = AllocMatI(NP, MAXBONDS);
  numbondst  = malloc(sizeof(int)*NP); 
  bondst = AllocMatI(NP, MAXBONDS);
  Fb[0] = malloc(sizeof(double)*points);
  Fb[1] = malloc(sizeof(double)*points);
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      r1[a] = malloc(sizeof(double)*NP);
    }

  for (ii=0; ii < points; ii++)
    {
      ti[ii] = -1.0;
      Fb[0][ii] = 0.0;
      Fb[1][ii] = 0.0;
    }
  first = 0;
  fclose(f2);
  for (ii=0; ii < points; ii++)
    cc[ii] = 0.0;
  if (particles_type == 1)
    {
      maxax0 = sa[0];
      if (sb[0] > maxax0)
	maxax0 = sb[0];
      if (sc[0] > maxax0)
	maxax0 = sc[0];
      maxax1 = sa[1];
      if (sb[1] > maxax1)
	maxax1 = sb[1];
      if (sc[0] > maxax1)
	maxax1 = sc[1];
      maxsax = fabs(maxax1)+fabs(maxax0)+2.0*sigmaSticky;
      //printf("maxsax=%.15G\n", maxsax);
    }
  else
    {
      maxsaxAA = fabs(sigmaAA)+2.0*sigmaSticky;
      maxsaxAB = fabs(sigmaAB)+2.0*sigmaSticky;
      maxsaxBB = fabs(sigmaBB)+2.0*sigmaSticky;
    }
  /* WARNING: se i diametri sono diversi va cambiato qua!! */ 
  if (particles_type == 1)
    RCUT = maxsax;
  else if (particles_type == 0)
    RCUT = maxsaxAA*1.01;

  for (a = 0; a < 3; a++)
    {
      for (b = 0; b < NA; b++)
	rat[b][a] = malloc(sizeof(double)*NP);
      r0[a] = malloc(sizeof(double)*NP);
      DR0[a] = malloc(sizeof(double)*NP);
      for (b = 0; b < 3; b++)
	{
	  R[a][b] = malloc(sizeof(double)*NP);
	}
    }
  if (particles_type==1)
    {
      build_atom_positions();
      printf("SYSTEM: ELLISPOIDS - DGEBA\n");
    }
  else
    {
      printf("SYSTEM: SPHERES 3-2\n");
    }
  if (NPA != NP)
    printf("[MIXTURE] files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", nfiles, NP, NPA, L, NN, maxl);
  else
    printf("[MONODISPERE] files=%d NP = %d L=%.15G NN=%d maxl=%d\n", nfiles, NP, L, NN, maxl);

  c2 = 0;
  JJ = 0;
  cellsx = L / RCUT;
  cellsy = L / RCUT;
  cellsz = L / RCUT;
  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+NP));
  inCell[0] = malloc(sizeof(int)*NP);
  inCell[1] = malloc(sizeof(int)*NP);
  inCell[2] = malloc(sizeof(int)*NP);

  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN)
    {	
      readconf(fname[nr1], &time, &refTime, NP, r0, DR0, R);
      /* costruisce la posizione di tutti gli sticky spots */
      build_spots_positions(r0, R);
      build_linked_list();
      for ( i = 0; i < NP; i++)
	numbonds0[i] = 0;
      find_bonds(numbonds0, bonds0, &ene0);
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      /* N.B. considera NN punti in maniera logaritmica e poi calcola i punti in maniera lineare 
	       * distanziati di NN punti. */
              np = (JJ == 0)?nr2-nr1:NN-1+JJ;	      
	      if (nr2 >= nfiles || np >= points)
		{
		  fine = 1;
		  break;
		}
	      if (JJ > 0 && (nr2 - nr1) % NN != 0)
		continue;
	      readconf(fname[nr2], &time, &refTime, NP, r1, DR0, R);
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("np=%d time=%.15G\n", np, ti[np]);
		}
	      /* accumulate bond correlation function */ 
	      build_spots_positions(r0, R);
	      build_linked_list();
	      for ( i = 0; i < NP; i++)
		numbondst[i] = 0;
	      find_bonds(numbondst, bondst, &enet);
	      for (i = 0; i < NP; i++)
		{
		  for (nb = 0; nb < numbondst[i]; nb++)
		    {
		      j = bondst[nb][i] / NA;
		      jj2 = bondst[nb][i] % (NA*NA);
		      a = jj2 / NA; 
		      b = jj2 % NA;
		      if (i < NPA)
			{
			  if (exist_bond(i, a, j, b, numbonds0, bonds0))
			    Fb[0][np] += 1.0;  
			}
		      else
			{
			  if (exist_bond(i, a, j, b, numbonds0, bonds0))
			    Fb[1][np] += 1.0;
			}
		    }
		}
	      cc[np] += 1.0;
	    }
	}
    }
  /* print correlation function to a file */
  if (NPA < NP)
    {
      f = fopen("FbondsA.dat", "w");
    }
  else
    {
      f = fopen("Fbonds.dat", "w");
    }
  for (ii=0; ii < points; ii++)
    {
      Fb[0][ii] /= ((double)cc[ii]);
    }
  for (ii=0; ii < points; ii++)
    {
      Fb[0][ii] /=  Fb[0][0];
    }
  for (ii=0; ii < points; ii++)
    {
      printf("%.15G %.15G\n", ti[ii], Fb[0][ii]);
    }
  fclose(f);
  if (NPA < NP)
    {
      f = fopen("FbondsB.dat", "w");
      for (ii=0; ii < points; ii++)
	{
	  Fb[1][ii] /= ((double)cc[ii]);
	}
      for (ii=0; ii < points; ii++)
	{
	  Fb[1][ii] /=  Fb[1][0];
	}
      for (ii=0; ii < points; ii++)
	{
	  fprintf(f, "%.15G %.15G\n", ti[ii], Fb[1][ii]);
	}
      fclose(f);
    }
  return 0;
}
