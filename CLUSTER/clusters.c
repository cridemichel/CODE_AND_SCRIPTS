#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 5000
#define NA 6
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2
#define MD_PBONDS 10
#define Sqr(x) ((x)*(x))
char **fname; 
double L, time, *ti, *R[3][3], *r0[3], r0L[3], RL[3][3], *DR0[3], maxsax, maxax0, maxax1,
       maxsaxAA, maxsaxAB, maxsaxBB;
double pi, sa[2]={-1.0,-1.0}, sb[2]={-1.0,-1.0}, sc[2]={-1.0,-1.0}, 
       Dr, theta, sigmaSticky, ratL[NA][3], *rat[NA][3], sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0;
int *dupcluster, RCUT, shift[3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int NP, NPA=-1, ncNV, ncNV2;
int check_percolation = 1, *nspots, particles_type=1;
/* particles_type= 0 (sphere3-2), 1 (ellipsoidsDGEBA) */ 
char inputfile[1024];
int foundDRs=0, foundrot=0, *color, *color2, *clsdim, *clsdim2, *clsdimNV, *clscolNV, *clscol, 
    *clsdimsort, *clssizedst, *clssizedstAVG, *percola;
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
		     &R[2][0][i], &R[2][1][i], &R[2][2][i]); 
	      //printf("%.15G %.15G %.15G\n", R[2][0][i],R[2][1][i], R[2][2][i] );
	    
	    }
	  break; 
	}

    }
  fclose(f);
}
#define MD_SP_DELR 0.0


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
  double ma;
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
  int maxa, maxb;
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
  int a, b, maxa, maxb;
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

void change_all_colors(int NP, int* color, int colorsrc, int colordst)
{
  int ii;
  for (ii = 0; ii < NP; ii++)
    {
      if (color[ii] == colorsrc)
	color[ii] = colordst;
    }
}
char fncls[1024];
int findmaxColor(int NP, int *color)
{
  int i, maxc=-1;
  for (i = 0; i < NP; i++) 
    {
      if (color[i] > maxc)
	maxc = color[i];
    }
  return maxc;
}

struct cluster_sort_struct { 
  int dim;
  int color;
};
struct cluster_sort_struct *cluster_sort;
int compare_func (const void *aa, const void *bb)
{
  int ai, bi;
  int temp;
  struct cluster_sort_struct *a, *b;
  a = (struct cluster_sort_struct*) aa;
  b = (struct cluster_sort_struct*) bb;
  ai = a->dim;
  bi = b->dim;
  temp = ai - bi;
  if (temp < 0)
    return 1;
  else if (temp > 0)
    return -1;
  else
    return 0;
}
#if 0
const int images_array[27][3]={{0,0,0},
{1,0,0},{0,1,0},{0,0,1},
{-1,0,0},{0,-1,0},{0,0,-1},
{1,1,0}, {0,1,1}, {1,0,1},
{-1,-1,0},{0,-1,-1},{-1,0,-1},
{-1,+1,0},{0,-1,+1},{-1,0,+1},
{+1,-1,0},{0,+1,-1},{+1,0,-1},
{1,1,1},{-1,-1,-1},
{-1,1,1},{1,-1,1},{1,1,-1},
{-1,-1,1},{1,-1,-1},{-1,1,-1}};
#else
const int images_array[8][3]={{0,0,0},
{1,0,0},{0,1,0},{0,0,1},
{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
#endif
void choose_image(int img, int *dix, int *diy, int *diz)
{
  *dix = images_array[img][0];
  *diy = images_array[img][1];
  *diz = images_array[img][2];
}
void print_usage(void)
{
  printf("Usage: clusters <listafile>\n");
  exit(0);
}

void parse_params(int argc, char** argv)
{
  int cc=1;
  if (argc <= 1)
    {
      print_usage();
    }
  while (cc < argc)
    {
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"--noperc") || !strcmp(argv[cc],"-np" ))
	{
	 check_percolation = 0;
	} 
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
}
int *inCell[3], *cellList, cellsx, cellsy, cellsz;
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

int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int  NN, fine, JJ, nat, maxl, maxnp, np, nc2, nc, dix, diy, diz, djx,djy,djz,imgi2, imgj2, jbeg, ifin;
  int jX, jY, jZ, iX, iY, iZ;
  //int coppie;
  double refTime=0.0, ti, ene=0.0;
  const int NUMREP = 8;
  int curcolor, ncls, b, j, almenouno, na, c, i2, j2, ncls2;
  pi = acos(0.0)*2.0;
    /* parse arguments */
  parse_params(argc, argv);

  f2 = fopen(inputfile, "r");
  c2 = 0;
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
  fclose(f2);
  f = fopen(fname[0], "r");
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
      else if (nat == 0 && !strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (nat==1 && !strcmp(parname,"a"))
	sscanf(parval, "%lf %lf\n", &sa[0], &sa[1]);	
      else if (nat==1 && !strcmp(parname,"b"))
	sscanf(parval, "%lf %lf\n", &sb[0], &sb[1]);	
      else if (nat==1 && !strcmp(parname,"c"))
	sscanf(parval, "%lf %lf\n", &sc[0], &sc[1]);	
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
  /* default = ellipsoids */
  //printf("sigmaAA=%.15G\n", sigmaAA);
  if (sigmaAA != -1.0)
    particles_type = 0;
  
  if (NPA == -1)
    NPA = NP;
  color = malloc(sizeof(int)*NP);
  color2= malloc(sizeof(int)*NP*NUMREP);
  clsdim2=malloc(sizeof(int)*NP*NUMREP);
  nspots = malloc(sizeof(int)*NP);
  clsdim = malloc(sizeof(int)*NP);
  clsdimNV = malloc(sizeof(int)*NP);
  clscolNV = malloc(sizeof(int)*NP);
  clscol   = malloc(sizeof(int)*NP);
  cluster_sort = malloc(sizeof(struct cluster_sort_struct)*NP);
  clssizedst = malloc(sizeof(int)*NP);
  clssizedstAVG = malloc(sizeof(int)*NP);
  dupcluster = malloc(sizeof(int)*NP*NUMREP); 
  
  percola = malloc(sizeof(int)*NP);
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
  cellsx = L / RCUT;
  cellsy = L / RCUT;
  cellsz = L / RCUT;
  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+NP));
  inCell[0] = malloc(sizeof(int)*NP);
  inCell[1] = malloc(sizeof(int)*NP);
  inCell[2] = malloc(sizeof(int)*NP);

  for (i = 0; i < NP; i++)
    {
      clssizedstAVG[i] = 0;
      clssizedst[i] = 0; 
      percola[i] = 0; 
    }
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
  //printf("sigmaSticky=%.15G\n", sigmaSticky);
  for (nr1 = 0; nr1 < nfiles; nr1++)
    {	
      readconf(fname[nr1], &time, &refTime, NP, r0, DR0, R);
      ti = time + refTime;
      /* costruisce la posizione di tutti gli sticky spots */
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
      for (i = 0; i < NP; i++)
	{
	  color[i] = -1;	  
	  clssizedst[i] = 0;
	}
      curcolor = 0;
      ene=0;
      //coppie = 0;
      build_linked_list();

      if (particles_type == 1)
	{
	  jbeg = NPA;
	  ifin = NPA;
	}
      else 
	{
	  jbeg = 0; 
	  ifin = NP;
	}
      for (i = 0; i < ifin; i++)
	{
	  if (particles_type == 1)
	    color[i] = curcolor;
	  else if (particles_type == 0)
	    { 
	      if (color[i] == -1)
		color[i] = curcolor;
	    }
	  for (iZ = -1; iZ <= 1; iZ++) 
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
	      for (iY = -1; iY <= 1; iY ++) 
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
			      if ((i < NPA && j < NPA) || ( i >= NPA && j >= NPA))
				continue;
			      break;
			    }
			  if (bond_found(i, j))  
			    {
			      ene=ene+1.0;
			      if (color[j] == -1)
				color[j] = color[i];
			      else
				{
				  if (color[i] < color[j])
				    change_all_colors(NP, color, color[j], color[i]);
				  else if (color[i] > color[j])
				    change_all_colors(NP, color, color[i], color[j]);
				}
			    }
			}
		    }
		}
	    }
#if 0 
	  for (j = jbeg; j < NP; j++)
	    {
	      //coppie++;
	      if (particles_type == 0 && j <= i) 
		continue;
      	      if (bond_found(i, j))  
		{
		  ene=ene+1.0;
		  if (color[j] == -1)
		    color[j] = color[i];
		  else
		    {
		      if (color[i] < color[j])
			change_all_colors(NP, color, color[j], color[i]);
		      else if (color[i] > color[j])
			change_all_colors(NP, color, color[i], color[j]);
		    }
		  
		}
	    }
#endif
	  curcolor = findmaxColor(NP, color)+1;
	}
      /* considera la particelle singole come cluster da 1 */
      for (i = 0; i < NP; i++)
	{
	  if (color[i]==-1)
	    {	    
	      color[i] = curcolor;
	      curcolor++;
	    } 
	  //printf("color[%d]=%d\n", i, color[i]);
	}
      ncls = curcolor;
      //printf("curcolor:%d\n", curcolor);
      sprintf(fncls, "%s.clusters", fname[nr1]);
      f = fopen(fncls, "w+");
      for (nc = 0; nc < ncls; nc++)
	{
	  clsdim[nc] = 0; 
	}
      for (nc = 0; nc < ncls; nc++)
	{
	  for (a = 0; a < NP; a++)
	    if (color[a] == nc)
	      {
		clsdim[color[a]]++;
		clscol[nc] = color[a];
	      }
	}
      //printf("NP=%d ncls=%d\n", NP, ncls);
      /*  ==== >>> REMOVE VOIDS <<< ==== */
      ncNV=0;
      for (nc = 0; nc < ncls; nc++)
	{
	  if (clsdim[nc] != 0)
	    {
	      clsdimNV[ncNV] = clsdim[nc];
	      clscolNV[ncNV] = clscol[nc]; 
	      ncNV++;
	    }

	}
      ncls = ncNV;
      printf("E/N = %.15G\n", ene/((double)NP));
      //printf("coppie PERC=%d\n", coppie);
      for (nc = 0; nc < ncls; nc++)
	{
	  //printf("clsdimNV[%d]=%d\n",nc ,clsdimNV[nc]);
	  cluster_sort[nc].dim = clsdimNV[nc];
	  cluster_sort[nc].color = clscolNV[nc];
	}
      qsort(cluster_sort, ncls, sizeof(struct cluster_sort_struct), compare_func);
      /* ============== >>> PERCOLATION <<< ================== */
      if (check_percolation)
	{

	  for (i=0; i < NP; i++)
	    {
	      if (particles_type == 1)
		{
		  if (i < NPA)
		    nspots[i] = MD_STSPOTS_A;
		  else
		    nspots[i] = MD_STSPOTS_B;		
		}
	      else
		{
		  if (i < NPA)
		    nspots[i] = 2;
		  else
		    nspots[i] = 3;		
		}
	    }	
	  ene=0;
	  //coppie = 0;
	  for (nc = 0; nc < ncls; nc++)
	    {
	      if (cluster_sort[nc].dim==1)
		continue;
	      //printf("Analysing cluster #%d of #%d\n", nc+1, ncls);
	      /* N.B per verificare la percolazione ogni cluster va "duplicato"
	       * in tutte le direzioni e se alla fine risulta comunque un unico 
	       * cluster allora tale cluster è percolante.*/
	      na = 0;
	      //printf("i=1011 j=277 rat=%.15G %.15G\n", rat[0][0][1011], rat[0][0][377]);
	      for (i=0; i < NP; i++)
		{
		  if (color[i]==cluster_sort[nc].color)
		    {
		      for (c = 0; c < NUMREP; c++)
			{
			  dupcluster[c*cluster_sort[nc].dim+na] = i;
			}
		      na++;
		    }
		}
	      //printf("i=1011 j=277 rat=%.15G %.15G\n", rat[0][0][1011], rat[0][0][377]);
	      //printf("NP=%d NPA=%d na=%d,clsdim[%d]=%d\n", NP, NPA,na,nc,cluster_sort[nc].dim);
	      curcolor = 0;
	      for (i2 = 0; i2 < na*NUMREP; i2++)
		{
		  color2[i2] = -1;	  
		}
	      for (i2 = 0; i2 < na*NUMREP; i2++)
		{
		  if (color2[i2]==-1)
		    color2[i2] = curcolor;
		  //printf("nc=%d na*NUMREP=%d i2=%d\n",nc, na*NUMREP, i2);
		  //printf("curcolor:%d\n", curcolor);
		  for (j2 = 0; j2 < na*NUMREP; j2++)
		    {
		      i = dupcluster[i2];
		      j = dupcluster[j2];
		      if (i2 >= j2)
			continue;
		      if (particles_type == 1)
			{
			  if ((nspots[i]==MD_STSPOTS_A && nspots[j]==MD_STSPOTS_A) ||
			      (nspots[i]==MD_STSPOTS_B && nspots[j]==MD_STSPOTS_B))
			    continue;
			}
		      //coppie++;
		      dix = diy = diz = 0;
		      djx = djy = djz = 0;
		      imgi2 = i2 / na;
		      imgj2 = j2 / na;
		      if (i2==j2)
			continue;
		      //printf("i2=%d j2=%d imgi2=%d imgj2=%d i=%d j=%d\n", i2, j2, imgi2, imgj2, i, j);
		      choose_image(imgi2, &dix, &diy, &diz);
		      choose_image(imgj2, &djx, &djy, &djz);
		      //if (dix!=0||diy!=0||diz!=0||djx!=0||djy!=0||djz!=0)
			//printf("(%d,%d,%d)-(%d,%d,%d)\n", dix, diy, diz, djx, djy, djz);
		      if ( bond_foundR(i, j, dix, diy, diz, djx, djy, djz, 2.0*L) )
			{
			  ene=ene+1.0;
			  //printf("qui!!!\n");
			  if (color2[j2] == -1)
			    {
			      color2[j2] = color2[i2];
			      //printf("color2[j2]=color2[i2]=%d\n", color2[i2]);
			    }
			  else
			    {
			      if (color2[i2] < color2[j2])
				change_all_colors(na*NUMREP, color2, color2[j2], color2[i2]);
			      else if (color2[i2] > color2[j2])
				change_all_colors(na*NUMREP, color2, color2[i2], color2[j2]);
			    }
			}
		    }
		  curcolor = findmaxColor(na*NUMREP, color2)+1;
		  //printf("curcolor2=%d\n", curcolor);
		}
	      ncls2 = curcolor;
	      for (nc2 = 0; nc2 < ncls2; nc2++)
		{
		  clsdim2[nc2] = 0; 
		}
	      for (nc2 = 0; nc2 < ncls2; nc2++)
		{
		  for (a = 0; a < na*NUMREP; a++)
		    if (color2[a] == nc2)
		      {
			clsdim2[color2[a]]++;
		      }
		}

	      /* ==== >>> REMOVE_VOIDS <<< ==== */
	      ncNV2=0;
	      for (nc2 = 0; nc2 < ncls2; nc2++)
		{
		  if (clsdim2[nc2] != 0)
		    {
		      ncNV2++;
		    }
		}
	      //printf("ncls2=%d\n", ncNV2);
	      if (ncNV2 < NUMREP)
		percola[nc] = 1;
	    }
	  printf("E/N (PERCOLATION) = %.15G\n", ene/((double)(NUMREP))/((double)NP));
	}
      //printf("coppie PERC=%d\n", coppie);
      almenouno = 0;
      for (nc = 0; nc < ncls; nc++)
	{
	  //if (cluster_sort[nc].dim >= 2)
	    //almenouno = 1;
	  if (percola[cluster_sort[nc].color])
	    fprintf(f, "1 ");
	  else
	    fprintf(f, "0 ");
	      
	  for (i = 0; i < NP; i++)
	    {
	      if (color[i]==cluster_sort[nc].color)
		{
		  fprintf(f, "%d ", i);
		}
	    }
	  fprintf(f, "\n");
	}
      //if (almenouno==0)
	//fprintf(f, "WARNING: No clusters found!\n");
      fclose(f);
      for (nc = 0; nc < ncls; nc++)
	{
	  //printf("cluster_sort[%d].dim=%d color=%d\n", nc, cluster_sort[nc].dim, cluster_sort[nc].color);
	  clssizedst[cluster_sort[nc].dim]++;
	  clssizedstAVG[cluster_sort[nc].dim]++;
	}
      sprintf(fncls, "%s.clsdst", fname[nr1]);
      f = fopen(fncls, "w+");
      for (i = 1; i < NP; i++)
	{
	  if (clssizedst[i] != 0)
	    fprintf(f, "%d %d\n", i, clssizedst[i]);
	}
      fclose(f);
    }
  f = fopen("avg_cluster_size_distr.dat", "w+");
  for (i = 1; i < NP; i++)
    {
      if (clssizedstAVG[i] != 0)
	fprintf(f, "%d %d\n", i, clssizedstAVG[i]/nfiles);
    }
  fclose(f);
  return 0;
}
