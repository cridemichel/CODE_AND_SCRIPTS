/* NOTA 09/11/10: questo programma per ora non supporta le interazioni specifiche tra particelle intersIJ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 10000
#define NA 100
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2
#define MD_PBONDS 10
#define Sqr(x) ((x)*(x))
/* get this from ellipsoid.h */
#define EDHE_FLEX
#ifdef EDHE_FLEX
#define NAFLX 1000
#define LINES_TO_SKIP 10
#define NANA (((long long int)NAFLX)*((long long int)NAFLX))
#endif
/* NOTA: 
 * particles_type == 0 ( DGEBA - sticky ellipsoid), 1 (sticky 2-3), 2 (bimixhs) */
int *mapbondsa, *mapbondsb, only_average_clsdistro=0;
int *mapbondsaFlex, *mapbondsbFlex, nbondsFlex;
double *mapBheightFlex, *mapBhinFlex, *mapBhoutFlex, *mapSigmaFlex; 
int **mapbondsaFlexS, **mapbondsbFlexS, *nbondsFlexS;
double **mapBheightFlexS, **mapBhinFlexS, **mapBhoutFlexS, **mapSigmaFlexS; 

double totdist=0.0, distcc=0.0;
typedef struct 
{
  int min;
  int max;
} rangeStruct;

typedef struct {
  double x[3];
  double sax[3];
  double ppsax[3]; /* semi-lati di un parallelepipedo che circoscrive l'ellissoide con tutti i suoi sticky spots */
  double ppr[3];
  double n[3];	
#if 0
  double R[3][3]; /* orientazione relativa al sistema di riferimento del corpo rigido */
#endif
} hardobjsStruct;
typedef struct {
  double x[3];
  double sigma;
  int same;
} spotStruct;
typedef struct 
{
  double sax[3];
  double ppsax[3]; /* semi-lati di un parallelepipedo che circoscrive l'ellissoide con tutti i suoi sticky spots */
  double ppr[3];
  double n[3];/*super-ellipsoids double parameters (as in Povray, i.e. only n[0] and n[1] are used, n[2] is a place holder 
		for now) */
  double m;
  double I[3];
#if 0
  double xoff[3];
#endif
  int ignoreCore;/* if 1 ignore collisions of the core object */
  int brownian;
  int nspots;
  spotStruct* spots; 
  int nhardobjs; /* one can add other sub-objects to build a super-object made of 
		   several super-ellipsoids with their spots */
  double rcutFact;
  hardobjsStruct* hardobjs;
} partType;
typedef struct 
{
  int type1; /* se type1 == -2 usa i range r1 */
  int spot1; 
  int type2; /* se type2 == -2 usa i range r2 */
  int spot2;
  double bheight;
  double bhin;
  double bhout;
  int nmax;
  int nr1; /* numero di range */
  rangeStruct* r1;
  int nr2; /* numero di range */
  rangeStruct* r2;
} interStruct;


interStruct* intersArr;
double MAXAX, maxax, maxSpots;
char **fname; 
int *typeNP=NULL;
partType *typesArr=NULL;
int *typeOfPart;
const int NUMREP = 8;
int MAXBONDS = 100;
double L, time, *ti, *R[3][3], *r0[3], r0LA[3], r0LB[3], RLA[3][3], RLB[3][3], *DR0[3], maxsax, maxax0, maxax1,
       maxsaxAA, maxsaxAB, maxsaxBB, RCUT;
double pi, sa[2]={-1.0,-1.0}, sb[2]={-1.0,-1.0}, sc[2]={-1.0,-1.0}, 
       Dr, theta, sigmaSticky=-1.0, ratLA[NA][3], ratLB[NA][3], *rat[NA][3], sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0;
double deltaAA=-1.0, deltaAB=-1.0, deltaBB=-1.0;
int *dupcluster, shift[3];
#ifdef EDHE_FLEX
long long int **bonds;
int *numbonds;
#else
int *numbonds, **bonds;
#endif
char parname[1024], parval[10000000], line[10000000];
char dummy[20480000];
int NP, NPA=-1, ncNV, ncNV2, START, END, NT, NI;
int check_percolation = 1, *nspots, output_bonds=0, mix_type=-1, saveBonds=-1;
int avgbondlen=1, calcordparam=0;
/* particles_type= 0 (sphere3-2), 1 (ellipsoidsDGEBA) */ 
char inputfile[1024];
int foundDRs=0, foundrot=0, *color, *color2, *clsdim, *clsdim2, *clsdimNV, *clscolNV, *clscol, 
    *clsdimsort, *clssizedst, *percola;
double *clssizedstAVG;
int media_log=0, mc_sim=0;
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
  static first=1;
  int nat=0, i, cpos, j;
  f = fopen(fname, "r");
  while (!feof(f)) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  continue;
	}

      if (nat==2)
	{
	  for (i = 0; i < typeNP[0]*4; i++) 
	    {
	      if (i > 0)
		fscanf(f, "%[^\n]\n", line); 
	      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %[^\n]\n", 
      		     &r[0][i], &r[1][i], &r[2][i], 
		     &R[0][0][i], &R[0][1][i], &R[0][2][i], &R[1][0][i], &R[1][1][i], &R[1][2][i],
		     &R[2][0][i], &R[2][1][i], &R[2][2][i], &(typeOfPart[i]), dummy); 
	      //printf("%.15G %.15G %.15G\n", R[2][0][i],R[2][1][i], R[2][2][i] );
	    
	    }
	  break; 
	}

    }
  fclose(f);
}
#if 0
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
#endif
/* array con le posizioni degli atomi nel riferimento del corpo rigido 
 * nel caso dell'acqua i siti idrogeno ed elettroni sono disposti su 
 * di un tetraedro */
#if 0
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
#endif
#if 0
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
#endif
void BuildAtomPosAtDNAD(int i, int ata, double *rO, double R[3][3], double rat[3])
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
      spXYZ = typesArr[typeOfPart[i]].spots[ata-1].x;
    }
  //radius = Oparams.sigma[0][1] / 2.0;
  if (ata == 0)
    {
      for (k1 = 0; k1 < 3; k1++)
	rat[k1] = rO[k1];
      //printf("%f %f %f @ 0.5 C[red]\n", rat[0], rat[1], rat[2]);
    }
  else 
    {
      for (k1 = 0; k1 < 3; k1++)
	{ 
	  rat[k1] = rO[k1];
	  for (k2 = 0; k2 < 3; k2++)
	    rat[k1] += R[k2][k1]*spXYZ[k2]; 
	}
#if 0
      printf("ata= %d rat= %f %f %f\n", ata, rat[0], rat[1], rat[2]);
      printf("rO = %f %f %f \n", rO[0], rO[1], rO[2]);
      printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
      printf("ata=%d %f %f %f @ 0.075 C[blue]\n", ata, R[0][ata-1], R[1][ata-1], R[2][ata-1]);
#endif
    }
  
}
void BuildAtomPosSQ(int i, double rO[3], double R[3][3], double rat[NA][3])
{
  int a;
#if 0
  for (a = 0; a < 3; a++)
    rat[1][a] = rO[a];
#endif   
  for (a=0; a < typesArr[typeOfPart[i]].nspots+1; a++)
    BuildAtomPosAtDNAD(i, a, rO, R, rat[a]);
}
#if 0
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
#endif
#define Sqr(x) ((x)*(x))

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
int is_in_ranges(int A, int B, int nr, rangeStruct* r)
{
  int kk;
  if (A==-1 && B==-1)
    return 0;
  if (A==B)
    return 1;
  if (B==-2)
    {
      for (kk=0; kk < nr; kk++)  
	{
	  if (A >= r[kk].min && A <= r[kk].max)
	    {
	      return 1;
	    }
	}
    }
  return 0;
}

void assign_bond_mapping(int i, int j)
{
  int ni, type1, type2, a, k, nl;
  type1 = typeOfPart[i];
  type2 = typeOfPart[j];
  a=0;
  
  /* nl is a unique number assigned to each type1-type2 pair */
  for (ni=0; ni < NI; ni++)
    {
      if (is_in_ranges(type1, intersArr[ni].type1, intersArr[ni].nr1, intersArr[ni].r1) && 
	  is_in_ranges(type2, intersArr[ni].type2, intersArr[ni].nr2, intersArr[ni].r2))
	{
	  /* N.B. il +1 c'è poiché mapbondsa[]=0 si riferisce all'atomo centrato nell'origine (il core) */
	  mapbondsaFlex[a] = intersArr[ni].spot1+1;
	  mapbondsbFlex[a] = intersArr[ni].spot2+1;
	  mapBheightFlex[a] = intersArr[ni].bheight;
	  mapBhinFlex[a] = intersArr[ni].bhin;
	  mapBhoutFlex[a] = intersArr[ni].bhout;
	  mapSigmaFlex[a] = 0.5*(typesArr[type1].spots[intersArr[ni].spot1].sigma
				 + typesArr[type2].spots[intersArr[ni].spot2].sigma);
	  a++;
	  if (type1 == type2 && intersArr[ni].spot1 != intersArr[ni].spot2)
	    {
	      mapbondsaFlex[a] = intersArr[ni].spot2+1;
	      mapbondsbFlex[a] = intersArr[ni].spot1+1;
	      mapBheightFlex[a] = mapBheightFlex[a-1];
	      mapBhinFlex[a] = mapBhinFlex[a-1];
	      mapBhoutFlex[a] = mapBhoutFlex[a-1];
	      mapSigmaFlex[a] = mapSigmaFlex[a-1];
	      a++;
	    }
	}	
      else if (is_in_ranges(type2, intersArr[ni].type1, intersArr[ni].nr1, intersArr[ni].r1) && 
	       is_in_ranges(type1, intersArr[ni].type2, intersArr[ni].nr2, intersArr[ni].r2))
	{
	  /* N.B. il +1 c'è poiché mapbondsa[]=0 si riferisce all'atomo centrato nell'origine (il core) */
	  mapbondsaFlex[a] = intersArr[ni].spot2+1;
	  mapbondsbFlex[a] = intersArr[ni].spot1+1;
	  mapBheightFlex[a] = intersArr[ni].bheight;
	  mapBhinFlex[a] = intersArr[ni].bhin;
	  mapBhoutFlex[a] = intersArr[ni].bhout;
	  mapSigmaFlex[a] = 0.5*(typesArr[type2].spots[intersArr[ni].spot1].sigma
				 + typesArr[type1].spots[intersArr[ni].spot2].sigma);
	  a++;
	}
    }
#if 0
  /* qui assenga le interazioni specifiche tra particelle i-j 
     (queste non le ottimizziamo) */
  for (ni=0; ni < Oparams.nintersIJ; ni++)
    {
      if (is_in_ranges(i, intersArrIJ[ni].i, intersArrIJ[ni].nr1, intersArrIJ[ni].r1) && 
	  is_in_ranges(j, intersArrIJ[ni].j, intersArrIJ[ni].nr2, intersArrIJ[ni].r2))
	{
	  mapbondsaFlex[a] = intersArrIJ[ni].spot1+1;
	  mapbondsbFlex[a] = intersArrIJ[ni].spot2+1;
	  mapBheightFlex[a] = intersArrIJ[ni].bheight;
	  mapBhinFlex[a] = intersArrIJ[ni].bhin;
          mapBhoutFlex[a] = intersArrIJ[ni].bhout;
	  mapSigmaFlex[a] = 0.5*(typesArr[type1].spots[intersArrIJ[ni].spot1].sigma
				 + typesArr[type2].spots[intersArrIJ[ni].spot2].sigma);
	  a++;
       	}	 
     else if (is_in_ranges(j, intersArrIJ[ni].i, intersArrIJ[ni].nr1, intersArrIJ[ni].r1) && 
	      is_in_ranges(i, intersArrIJ[ni].j, intersArrIJ[ni].nr2, intersArrIJ[ni].r2)) 
       {
	 mapbondsaFlex[a] = intersArrIJ[ni].spot2+1;
	 mapbondsbFlex[a] = intersArrIJ[ni].spot1+1;
	 mapBheightFlex[a] = intersArrIJ[ni].bheight;
	 mapBhinFlex[a] = intersArrIJ[ni].bhin;
	 mapBhoutFlex[a] = intersArrIJ[ni].bhout;
	 mapSigmaFlex[a] = 0.5*(typesArr[type2].spots[intersArrIJ[ni].spot1].sigma
				+ typesArr[type1].spots[intersArrIJ[ni].spot2].sigma);
	 a++;	 
       } 
    }
#endif
  nbondsFlex = a;
  mapbondsa = mapbondsaFlex;
  mapbondsb = mapbondsbFlex;
} 

char fncls[1024];
char fn[1024];
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
  printf("Usage: clusters [-mc] <listafile>\n");
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
      else if (!strcmp(argv[cc],"--montecarlo") || !strcmp(argv[cc],"-mc" ))
	{
	   mc_sim = 1;
	}
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
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
/* Allocate memory for a matrix of integers */
long long int** AllocMatLLI(int size1, int size2)
{
  long long int** v;
  int k;
  v = (long long int**) malloc(size1 * sizeof(long long int*));
  v[0] = (long long int*) malloc(size1 * size2 * sizeof(long long int));
  for (k = 1; k < size1; k++)
    v[k] = v[k-1] + size2;
  return v;
}
double eval_max_dist_for_spots(int pt)
{
  int ns;
  double dist, distMax=0.0;
  for (ns=0; ns < typesArr[pt].nspots; ns++)
    {
      dist = calc_norm(typesArr[pt].spots[ns].x) + typesArr[pt].spots[ns].sigma*0.5;
      if (dist > distMax)
	distMax = dist;
    }
  return distMax;
}
void parse_one_range(char *s, int *A, int *nr, rangeStruct **r)
{
  int m, M;
  if (sscanf(s, "%d-%d", &m, &M)==2)
    {
      *A = -2;
      if (*nr == 0)
	*r = malloc(sizeof(rangeStruct)*(*nr+1));
      else
	*r = realloc(*r, sizeof(rangeStruct)*(*nr+1));
      (*r)[*nr].min = m;
      (*r)[*nr].max = M;
      (*nr)++;
    }
  else
    {
      *A = atoi(s);
      *r = NULL;
      *nr = 0;
    }
}
void parse_ranges(char *s, int *A, int *nr, rangeStruct **r)
{
  char *ns;
  ns = strtok(s, ",");

  *nr = 0;
  if (!ns)
    {
      parse_one_range(s, A, nr, r);
      return;
    }
  while (ns)
    {
      parse_one_range(ns, A, nr, r);
      ns = strtok(NULL, ",");
    } 
}
int get_max_nbonds(void)
{
  int pt1, pt2, i, maxijbonds=0;
  int ni, type1, type2, a, maxpbonds=0;
  int *intersI;
  for (pt1=0; pt1 < NT; pt1++)
    {
      for (pt2=pt1; pt2 < NT; pt2++)
	{
	  type1 = pt1;
	  type2 = pt2;
	  a=0;
	  for (ni=0; ni < NI; ni++)
	    {
	      if ((intersArr[ni].type1 == type1 && intersArr[ni].type2 == type2) ||
		  (intersArr[ni].type1 == type2 && intersArr[ni].type2 == type1))
		{
		  a+=2;
		}	
	    }
	  if (a > maxpbonds)
	    maxpbonds = a;
	}
    }
#if 0
  intersI = malloc(sizeof(int)*Oparams.parnum);
  for (i=0; i < Oparams.parnum; i++)
    intersI[i] = 0;
  for (ni=0; ni < Oparams.nintersIJ; ni++)
    {
      (intersI[intersArrIJ[ni].i])+=2; 
      (intersI[intersArrIJ[ni].j])+=2;
    }
  maxijbonds=0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (intersI[i] > maxijbonds)
	maxijbonds = intersI[i];
    }
  free(intersI);
#else
  maxijbonds = 0;
#endif
  return maxpbonds+maxijbonds;
}

#define npmax 10001
const int nlin=20;
int l1[npmax], l2[npmax];
double dlog[npmax], xlog[npmax];
int maxnbonds;
double antpos[3], ndb;
double pos4dist[2][3];
double *g0, delr;
void readCorIni(char* cif)
{
  FILE *f;
  char *s1, *s2;
  int nat, i, a, b, j, beg;
  f = fopen(cif, "r");
  nat = 0;
  printf("reading: %s\n", cif);
  while (!feof(f) && nat < 2) 
    {
      fscanf(f, "%[^\n]\n)", line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  if (nat==1)
	    {
	      double rcutFact;
	      /* read spots in HRB ref. system */
	      if (!typeNP)
		{
		  typeNP = malloc(sizeof(int)*NP);
		  typesArr = malloc(sizeof(partType)*NT);
		}
	      printf("NT=%d\n", NT);
	      fscanf(f, "%s ", line);
	      if (!strcmp(line, "RF"))
		{ 
		  //printf("line=%s\n", line);
		  beg=0;
		  for (i=0; i < NT; i++)
		    {
		      fscanf(f, "%lf ", &(typesArr[i].rcutFact));
		      printf("typeArr[%d].rcutFact=%f\n", i, typesArr[i].rcutFact);
		    }
		}
	      else
		{
		  sscanf(line, "%d", &typeNP[0]);
		  beg=1;
		  for (i=0; i < NT; i++)
		    typesArr[i].rcutFact = -1.0;

		}
	      for (i=beg; i < NT; i++)
		{
		  fscanf(f, "%d ", &typeNP[i]);
		  //printf("typeNP[%d]=%d\n", i, typeNP[i]);
		}

#if 0
	      for (i=0; i < NT; i++)
		{
		  sscanf(f, "%d ", &typeNP[i]);
		}
#endif
	      for (i=0; i < NT; i++)
		{
		  /* read particles parameters */
		  fscanf(f, "%lf %lf %lf ", &typesArr[i].sax[0], &typesArr[i].sax[1], &typesArr[i].sax[2]); 
		  printf("sax[%d]= %f %f %f\n", i, typesArr[i].sax[0], typesArr[i].sax[1], typesArr[i].sax[2]);
		  fscanf(f, "%lf %lf %lf ", &typesArr[i].n[0], &typesArr[i].n[1], &typesArr[i].n[2]);
		  fscanf(f, "%lf %lf %lf %lf %d %d ", &typesArr[i].m, &typesArr[i].I[0], &typesArr[i].I[1],
			 &typesArr[i].I[2], &typesArr[i].brownian, &typesArr[i].ignoreCore);
		  /* read sticky spots parameters */
		  fscanf(f, "%d %d ", &typesArr[i].nspots, &typesArr[i].nhardobjs);
		  if (typesArr[i].nspots >= NA)
		    {
		      printf("[ERROR] too many spots (%d) for type %d increase NA (actual value is %d) in ellipsoid.h and recompile\n",
			     typesArr[i].nspots, i, NA);
		      exit(-1);
		    }
		  if (typesArr[i].nspots > 0)
		    {
		      typesArr[i].spots = malloc(sizeof(spotStruct)*typesArr[i].nspots);
		    }
		  else
		    typesArr[i].spots = NULL;
		  for (j = 0; j < typesArr[i].nspots; j++)
		    fscanf(f, "%lf %lf %lf %lf ", &typesArr[i].spots[j].x[0],&typesArr[i].spots[j].x[1],
			   &typesArr[i].spots[j].x[2], &typesArr[i].spots[j].sigma);

		}
	      if (NI > 0)
		{
		  intersArr = malloc(sizeof(interStruct)*NI);
		  maxnbonds = get_max_nbonds();
		  mapbondsaFlex = (int*)malloc(sizeof(int)*maxnbonds);
		  mapbondsbFlex = (int*)malloc(sizeof(int)*maxnbonds);
		  mapBheightFlex = (double*) malloc(sizeof(double)*maxnbonds);
		  mapBhinFlex    = (double*)malloc(sizeof(double)*maxnbonds);
		  mapBhoutFlex   = (double*)malloc(sizeof(double)*maxnbonds);
		  mapSigmaFlex   = (double*)malloc(sizeof(double)*maxnbonds);
		}
	      else
		intersArr = NULL;
	      s1 = malloc(sizeof(char)*65535);
	      s2 = malloc(sizeof(char)*65535);

	      for (i=0; i < NI; i++)
		{
		  fscanf(f, "%s %d %s %d %lf %lf %lf %d ", s1, &intersArr[i].spot1, s2, 
			 &intersArr[i].spot2, 
			 &intersArr[i].bheight, &intersArr[i].bhin, &intersArr[i].bhout, &intersArr[i].nmax);
		  parse_ranges(s1, &intersArr[i].type1, &intersArr[i].nr1, &intersArr[i].r1);
		  parse_ranges(s2, &intersArr[i].type2, &intersArr[i].nr2, &intersArr[i].r2);
		} 
	      free(s1);
	      free(s2);
	    }
	  if (nat==2)
	    {
	      int l2s;
	      double dummy1, dummy2;
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
	      typeOfPart = malloc(sizeof(int)*NP);
	      for (i = 0; i < NP; i++) 
		{
		  fscanf(f, "%[^\n]\n", line); 
		  //printf("NP=%d i=%d line=%s\n", NP, i, line);
		  sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %[^\n]\n", 
			 &r0[0][i], &r0[1][i], &r0[2][i], 
			 &R[0][0][i], &R[0][1][i], &R[0][2][i], &R[1][0][i], &R[1][1][i], &R[1][2][i],
			 &R[2][0][i], &R[2][1][i], &R[2][2][i], &(typeOfPart[i]), dummy); 
		  //printf("%.15G %.15G %.15G\n", R[2][0][i],R[2][1][i], R[2][2][i] );

		}
	      for (i = 0; i < NP; i++) 
		{
		  fscanf(f, "%[^\n]\n", line); 
		}	  
	      fscanf(f, "%lf\n", &L);
	      printf("====> L=%f\n", L);
	      break;
	    }
	  continue;
	}
      sscanf(line, "%[^:]:%[^\n]\n", parname, parval); 
      if (!strcmp(parname,"parnum"))
	NP = atoi(parval);
      else if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      //else if (nat == 0 && !strcmp(parname,"NN"))
      //NN = atoi(parval);
      else if (nat==0 && !strcmp(parname,"saveBonds"))
	{
	  sscanf(parval, "%d\n", &saveBonds);	
	  //particles_type = 3;
	}
      else if (nat==0 && !strcmp(parname,"a"))
	sscanf(parval, "%lf %lf\n", &sa[0], &sa[1]);	
      else if (nat==0 && !strcmp(parname,"b"))
	sscanf(parval, "%lf %lf\n", &sb[0], &sb[1]);	
      else if (nat==0 && !strcmp(parname,"c"))
	sscanf(parval, "%lf %lf\n", &sc[0], &sc[1]);	
      else if (nat==0 && !strcmp(parname,"sigma"))
	sscanf(parval, "%lf %lf %lf %lf\n", &sigmaAA, &sigmaAB, &sigmaAB, &sigmaBB);	
      else if (nat==0 && !strcmp(parname,"delta"))
	sscanf(parval, "%lf %lf %lf %lf\n", &deltaAA, &deltaAB, &deltaAB, &deltaBB);	
      else if (nat==0 && !strcmp(parname,"sigmaSticky"))
	sigmaSticky = atof(parval);
      else if (nat==0 && !strcmp(parname,"theta"))
	theta = atof(parval);
      else if (nat==0 && !strcmp(parname,"Dr"))
	Dr = atof(parval);
      else if (nat==0 && !strcmp(parname,"ntypes"))
	NT = atoi(parval);
      else if (nat==0 && !strcmp(parname,"ninters"))
	NI = atoi(parval);
    }
  fclose(f);
}

int main(int argc, char **argv)
{
  int kk, kmax, kj, i3, j, nbonds, bin, points;
  long long int jj2, aa, bb;
  double am, xmed, dist, distSq, rlower, rupper, nIdeal, g0m, r, cost, dx, dy;
  FILE *f, *f2, *f3;
  char *s1, *s2;
  int beg, c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int  NN, fine, JJ, nat, maxl, maxnp, np, nc2, nc, dix, diy, diz, djx,djy,djz,imgi2, imgj2, jbeg, ifin;
  int jX, jY, jZ, iX, iY, iZ, jj;
  //int coppie;
  double refTime=0.0, ti, ene=0.0;
  int curcolor, ncls, b, almenouno, na, c, i2, j2, ncls2;
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
  readCorIni("CorIni");
  printf("NP=%d\n", NP);
  nspots = malloc(sizeof(int)*NP);
  points=100;
  delr= 15.0/points;

  /* CALCOLA I MAXAX QUI !!! servono per la costruzione delle LL */
  /* WARNING: se i diametri sono diversi va cambiato qua!! */ 
  g0 = malloc(sizeof(double)*points);
  for (ii=0; ii < points; ii++)
    g0[ii] = 0.0;
  ndb=0;
  //printf("r0[1003]= %f %f %f\n", r0[0][1003], r0[1][1003], r0[2][1003]);
  for (nr1=0; nr1 < nfiles; nr1++)
    {
      
      if (fname[nr1][strlen(fname[nr1])-1]=='i' &&
	  fname[nr1][strlen(fname[nr1])-3]=='I')
	{
	  readCorIni(fname[nr1]);
	  continue;
	}
      printf("reading conf: %s\n", fname[nr1]);
      readconf(fname[nr1], &time, &refTime, NP, r0, DR0, R);
      
      for (i=0; i < typeNP[0]; i++)
	{
      	  for (a = 0; a < 3; a++)
	    {
	      r0LA[a] = r0[a][i*4];
	      for (b = 0; b < 3; b++)
		RLA[a][b] = R[a][b][i*4];
	      //printf("r0=%f %f %f R=%f %f %f\n", r0L[0], r0L[1], r0L[2], RL[0][0], RL[0][1], RL[0][2]);
	    }
	  for (a = 0; a < 3; a++)
	    {
	      r0LB[a] = r0[a][i*4+1];
	      for (b = 0; b < 3; b++)
		RLB[a][b] = R[a][b][i*4+1];
	      //printf("r0=%f %f %f R=%f %f %f\n", r0L[0], r0L[1], r0L[2], RL[0][0], RL[0][1], RL[0][2]);
	    }

	  /* calc spot positions of two FABs */
       	  BuildAtomPosSQ(i*4, r0LA, RLA, ratLA);
	  BuildAtomPosSQ(i*4+1, r0LB, RLB, ratLB);
	  //printf("A typeOfPart[%d]=%d\n", i*4, typeOfPart[i*4]);
	  //printf("B typeOfPart[%d]=%d\n", i*4+1, typeOfPart[i*4+1]);
	  /* check double bonding here */
	  nbonds=0;
	  //printf("A i=%d %f %f %f\n", i*4, RLA[0][2],RLA[1][2],RLA[2][2]);
	  //printf("B i=%d %f %f %f\n", i*4+1, RLB[0][2],RLB[1][2],RLB[2][2]);

	  //printf("typeof A=%d i=%d %f %f %f\n", typeOfPart[i*4], i*4, ratLA[2][0],ratLA[2][1],ratLA[2][2]);
	  //printf("typeof B=%d i=%d %f %f %f\n", typeOfPart[i*4+1], i*4+1, ratLB[2][0],ratLB[2][1],ratLB[2][2]);
	  for (j=typeNP[0]*4; j < NP; j++)
	    {
	      //printf("r0[1003]= %f %f %f\n", r0[0][1003], r0[1][1003], r0[2][1003]);
	      for (kk=0; kk < 3; kk++)
		antpos[kk] = r0[kk][j];
	      //printf("j=%d %f %f %f\n", j, antpos[0], antpos[1], antpos[2]);
	      distSq = 0.0;
	      dx = ratLA[2][0]-antpos[0];
	      dy = ratLA[2][1]-antpos[1];
	      dx = dx - L*rint(dx/L);
	      dy = dy - L*rint(dy/L);
	      distSq = Sqr(dx)+Sqr(dy)+Sqr(ratLA[2][2]-antpos[2]);
	      if (distSq < Sqr((0.612+0.79)*0.5))
		{
		  for (kk=0; kk < 3; kk++)
		    pos4dist[nbonds][kk] = antpos[kk];
		  nbonds+=1;
		  if (nbonds < 2)
		    continue;	 
		 else 
		   break; 
		}
	      distSq = 0.0;
	      dx = ratLB[2][0]-antpos[0];
	      dy = ratLB[2][1]-antpos[1];
	      dx = dx - L*rint(dx/L);
	      dy = dy - L*rint(dy/L);
	      distSq = Sqr(dx)+Sqr(dy)+Sqr(ratLB[2][2]-antpos[2]);
	      //printf("dist=%f\n", sqrt(distSq)); 
	      if (distSq < Sqr((0.612+0.79)*0.5))
		{
		  for (kk=0; kk < 3; kk++)
		    pos4dist[nbonds][kk] = antpos[kk];
		  nbonds+=1;
		}
	      if (nbonds==2)
		{
		  break;
		}
	    }
	  if (nbonds==2)
	    {
	      distSq=0.0;
	      for (kk=0; kk < 3; kk++)
		distSq += Sqr(pos4dist[0][kk]-pos4dist[1][kk]);
	      dist = sqrt(distSq);
	      bin = ((int) (dist / delr)); 
              ndb++;
	      if (bin < points && bin >= 0)
  		{
  		  g0[bin] += 2.0;
  		  //printf("g0[%d]=%.15G\n", bin, g0[bin]);
  		}
	    }
	}
    }
    
  f2 = fopen("grBIV.dat","w+");
  cost = 3.14159265359;//*ndb/(L*L); 
  r=delr*0.5; 
  printf("ndb=%f\n", ndb);
  for (i=0; i < points; i++)
    {
      rlower = ( (double) i )*delr;
      rupper = rlower + delr;
      nIdeal= cost * (Sqr(rupper)-Sqr(rlower));
      g0m = g0[i] /((double)ndb)/nIdeal;
      fprintf(f2, "%.15G %.15G\n", r, g0m);
      r += delr;
    }
  fclose(f2);

  return 0;
}
