#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 10000
#define NA 6
#define MD_STSPOTS_A 5
#define MD_STSPOTS_B 2
#define MD_PBONDS 10
#define Sqr(x) ((x)*(x))

#define npmax 10000001
const int nlin=20;
int l1[npmax], l2[npmax];
double dlog[npmax], xlog[npmax];
int kk, kmax, kj, i3, block, dummyint;
double am, xmed;
int *ip;
/* NOTA: 
 * particles_type == 0 ( DGEBA - sticky ellipsoid), 1 (sticky 2-3), 2 (bimixhs) */
char **fname; 

const int NUMREP = 8;
int MAXBONDS = 10;
double wellWidth;
double Lx, Ly, Lz, L, time, *ti, *R[3][3], *r0[3], r0L[3], RL[3][3], *DR0[3], maxsax, maxax0, maxax1,
       maxsaxAA, maxsaxAB, maxsaxBB, RCUT;
double pi, sa[2]={-1.0,-1.0}, sb[2]={-1.0,-1.0}, sc[2]={-1.0,-1.0}, 
       Dr, theta, sigmaSticky=-1.0, sigmaAA=-1.0, sigmaAB=-1.0, sigmaBB=-1.0;
double deltaAA=-1.0, deltaAB=-1.0, deltaBB=-1.0;
int *dupcluster, shift[3], *numbonds, **bonds;
char parname[128], parval[256000], line[256000];
char dummy[2048];
int NP, NPA=-1, ncNV, ncNV2, START, END;
int check_percolation = 1, *nspots, output_bonds=0, mix_type=-1, media_log=0;
/* particles_type= 0 (sphere3-2), 1 (ellipsoidsDGEBA) */ 
char inputfile[1024];
int foundDRs=0, foundrot=0, *color, *color2, *clsdim, *clsdim2, *clsdimNV, *clscolNV, *clscol, 
    *clsdimsort, *clssizedst, *percola;
double *clssizedstAVG;
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
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3])
{
  FILE *f;
  int nat=0, i, cpos;
  f = fopen(fname, "r");
  fscanf(f,"%[^\n]\n", line);
  fscanf(f,"%[^\n]\n", line);
  fscanf(f,"%lf %lf %lf ", &Lx, &Ly, &Lz);
  fscanf(f,"%[^\n]\n", line);
  fscanf(f,"%[^\n]\n", line);
  //cpos = ftell(f);
  //printf("cpos=%d\n", cpos);
  //fscanf(f, "%[^\n]\n",line);
  for (i = 0; i < NP; i++) 
    {
      fscanf(f, "%lf %lf %lf\n", &(r[0][i]), &(r[1][i]), &(r[2][i]), &(ip[i])); 
      //printf("%.15G %.15G %.15G\n", R[2][0][i],R[2][1][i], R[2][2][i] );
      //printf("%f, %f, %f\n", r[0][i], r[1][i], r[2][i]);
      r[0][i] -= Lx*0.5;
      r[1][i] -= Ly*0.5;
      r[2][i] -= Lz*0.5;
    }
    
  fclose(f);
}
#define MD_SP_DELR 0.0


double spApos[MD_STSPOTS_A][3] = {{MD_SP_DELR, 0.54, 0.0},{MD_SP_DELR, 0.54, 3.14159},{MD_SP_DELR, 2.60159,0.0},
    {MD_SP_DELR, 2.60159, 3.14159},{MD_SP_DELR, 1.5708, 0.0}};
double spBpos[MD_STSPOTS_B][3] = {{MD_SP_DELR, 0.0, 0.0},{MD_SP_DELR, 3.14159, 0.0}};

double spXYZ_A[MD_STSPOTS_A][3];
double spXYZ_B[MD_STSPOTS_B][3];

#define Sqr(x) ((x)*(x))
#if 0
int check_distance(int i, int j, double Dx, double Dy, double Dz)
{
  double DxL, DyL, DzL;
  double ma=0.0;
  DxL = fabs(Dx);
  DyL = fabs(Dy);
  DzL = fabs(Dz);

  ma = maxsaxBB;

  if (DxL > ma || DyL > ma || DzL > ma)
    return 1;
  else 
    return 0;
}

#endif
double distance(int i, int j)
{
  int a, b;
  int maxa=0, maxb=0;
  double imgx, imgy, imgz;
  double Dx, Dy, Dz;
  //double wellWidth;

  Dx = r0[0][i] - r0[0][j];
  Dy = r0[1][i] - r0[1][j];
  Dz = r0[2][i] - r0[2][j];
  imgx = -L*rint(Dx/L);
  imgy = -L*rint(Dy/L);
  imgz = -L*rint(Dz/L);

#if 0
  if (check_distance(i, j, Dx+imgx, Dy+imgy, Dz+imgz))
    return 1;
#endif
  //wellWidth = sigmaBB;
  if (Sqr(r0[0][i] + imgx - r0[0][j])+Sqr(r0[1][i] + imgy - r0[1][j])
      +Sqr(r0[2][i] + imgz - r0[2][j]) < Sqr(wellWidth))	  
    return -1;
	
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
  
  //wellWidth = sigmaBB;
    
  dx = L*(imgix-imgjx);
  dy = L*(imgiy-imgjy);
  dz = L*(imgiz-imgjz);
  Dx = r0[0][i] - r0[0][j] + dx;
  Dy = r0[1][i] - r0[1][j] + dy;
  Dz = r0[2][i] - r0[2][j] + dz;
  imgx = -Lbig*rint(Dx/Lbig);
  imgy = -Lbig*rint(Dy/Lbig);
  imgz = -Lbig*rint(Dz/Lbig);

#if 0
  if (check_distance(i, j, Dx + imgx, Dy + imgy, Dz + imgz))
    return 1;
#endif
  if (Sqr(r0[0][i] + dx + imgx - r0[0][j])+Sqr(r0[1][i] + dy + imgy - r0[1][j])
      + Sqr(r0[2][i] + dz + imgz - r0[2][j]) < Sqr(wellWidth))	  
    return -1;

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
char fn[1024];
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
{1,0,0},{0,1,0},{0,0,1},{1,1,1},{1,1,0},{0,1,1},{1,0,1}};
#endif
void choose_image(int img, int *dix, int *diy, int *diz)
{
  *dix = images_array[img][0];
  *diy = images_array[img][1];
  *diz = images_array[img][2];
}
void print_usage(void)
{
  printf("Usage: clusters [--ptype/-pt] [--noperc/-np] [ --medialog/-ml ] [--bonds/-b] [--maxbonds] [--wellwidth/-ww <value>] <listafile>\n");
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
      else if (!strcmp(argv[cc],"--ptype") || !strcmp(argv[cc],"-pt" ))
	{
	  cc++;
          if (cc == argc)
	     print_usage();
	  mix_type = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--noperc") || !strcmp(argv[cc],"-np" ))
	{
	  check_percolation = 0;
	} 
      else if (!strcmp(argv[cc],"--medialog") || !strcmp(argv[cc],"-ml" ))
	{
	   media_log = 1;
	}
      else if (!strcmp(argv[cc],"--bonds") || !strcmp(argv[cc],"-b" ))
	{
	  output_bonds = 1;
	} 
      else if (!strcmp(argv[cc],"-ww") || !strcmp(argv[cc],"--wellwidth" ))
	{
	  cc++;
	  wellWidth = atof(argv[cc]);
	} 
      else if (!strcmp(argv[cc],"--maxbonds") || !strcmp(argv[cc],"-mb" ))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  MAXBONDS = atoi(argv[cc]);
 	  output_bonds = 1;
	} 
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
}
int *inCell[3]={NULL,NULL,NULL}, *cellList=NULL, cellsx, cellsy, cellsz;
void build_linked_list(void)
{
  double L2;
  int j, n;
  L2 = 0.5 * L;

  for (j = 0; j < cellsx*cellsy*cellsz + NP; j++)
    cellList[j] = -1;

  for (n = START; n < END; n++)
    {
      inCell[0][n] =  (cylinders[n].r[0] + L2) * cellsx / L;
      inCell[1][n] =  (cylinders[n].r[1] + L2) * cellsy / L;
      inCell[2][n] =  (cylinders[n].r[2] + L2) * cellsz / L;
      if (inCell[0][n] == cellsx)
	inCell[0][n]=cellsx-1;
      if (inCell[1][n] == cellsy)
	inCell[1][n]=cellsy-1;
      if (inCell[2][n] == cellsz)
	inCell[2][n]=cellsz-1;
      if (inCell[0][n] == -1)
	inCell[0][n]=0;
      if (inCell[1][n] == -1)
	inCell[1][n]=0;
      if (inCell[2][n] == -1)
	inCell[2][n]=0;
      j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + NP;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}
void build_linked_list_perc(int clsdim, double Lbig)
{
  double L2;
  int img, j, n, np, dix, diy, diz;
  L2 = 0.5 * L;

  for (j = 0; j < cellsx*cellsy*cellsz + NP*NUMREP; j++)
    cellList[j] = -1;
  //printf("cells=%d %d %d Lbig=%.15G L=%.15G\n", cellsx, cellsy, cellsz, Lbig, L);
  for (n = 0; n < clsdim*NUMREP; n++)
    {
      img = n / clsdim;
      choose_image(img, &dix, &diy, &diz);
      np = dupcluster[n];
      inCell[0][n] =  (r0[0][np] + dix*L + L2) * cellsx / Lbig;
      inCell[1][n] =  (r0[1][np] + diy*L + L2) * cellsy / Lbig;
      inCell[2][n] =  (r0[2][np] + diz*L + L2) * cellsz / Lbig;
      if (inCell[0][n] == cellsx)
	inCell[0][n]=cellsx-1;
      if (inCell[1][n] == cellsy)
	inCell[1][n]=cellsy-1;
      if (inCell[2][n] == cellsz)
	inCell[2][n]=cellsz-1;
      if (inCell[0][n] == -1)
	inCell[0][n]=0;
      if (inCell[1][n] == -1)
	inCell[1][n]=0;
      if (inCell[2][n] == -1)
	inCell[2][n]=0;
   j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + NP*NUMREP;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}
void add_bond(int i, int j)
{
  bonds[i][numbonds[i]] = j;
  numbonds[i]++;
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
double calcDistNegHC(int i, int j, double shift[3], int* retchk)
{
  const int MAX_ITERATIONS = 1000000;
#ifdef MC_HC_SPHERO_OPT
  int rim;
  double sphov;
#endif
  int it, k2;
  double normNSq, ViVj[3], lambdai, lambdaj;
  double sp, Q1, Q2, normPiDi, normPjDj, normN, L, D, DiN, DjN, niN[3], njN[3], Djni, Djnj;
  double PiPj[3], N[3], Pi[3], Pj[3], VV[3], Di[2][3], Dj[2][3], ni[3], nj[3], Ci[3], Cj[3];
  double normPiPj, Ui[3], DiCi[3], DiCini, normDiCi, DjCi[3], normDjCi;
  double PiDi[3], PjDj[3], Ai[3], Tj[3], Tjp[3], Tjm[3], TjpCi[3], TjmCi[3], TjpCini, TjmCini;
  double DjUini, DjUi[3], normDjUi, AiDj[3], AiDjnj, AiDjnjvec[3], TjNew[3], TjNewCi[3], TjNewCini;
  double TjOld[3], ninj, CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3], TipCjnj, TimCjnj;
  double Aj[3], AjDini, AjDinivec[3], AjDi[3], Tip[3], Tim[3], TipCj[3], TimCj[3], Dini;
  double DiCj[3], normDiCj, DiCjnj, Uj[3], DiUj[3], normDiUj, DiUjnj;
  double Tim_perp[3], Tip_perp[3], Tim_para[3], Tip_para[3], normTim_perp, DjCini;
  double Tjm_perp[3], Tjp_perp[3], Tjm_para[3], Tjp_para[3], normTjm_perp;
  double TiOld[3], TiNew[3], TiNewCj[3], TiNewCjnj;	
  double normCiCj;	
  double DjTmp[2][3], CiTmp[3], niTmp[3], njTmp[3];
  int kk, j1, j2;

  *retchk = 0; 

  for (kk=0; kk < 3; kk++)
    {
      ni[kk] = R[i][0][kk];
      nj[kk] = R[j][0][kk];
    }
  Ci[0] = rx[i];
  Ci[1] = ry[i];
  Ci[2] = rz[i]; 
  Cj[0] = rx[j] + shift[0];
  Cj[1] = ry[j] + shift[1];
  Cj[2] = rz[j] + shift[2]; 
  L = 2.0*typesArr[typeOfPart[i]].sax[0];
  D = 2.0*typesArr[typeOfPart[i]].sax[1];
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }

  for (kk=0; kk < 3; kk++)
    {
      /* centers of mass of disks */
      Di[0][kk]=Ci[kk]+0.5*L*ni[kk];
      Di[1][kk]=Ci[kk]-0.5*L*ni[kk];
      Dj[0][kk]=Cj[kk]+0.5*L*nj[kk];
      Dj[1][kk]=Cj[kk]-0.5*L*nj[kk];
    }
#ifdef MC_HC_SPHERO_OPT
  if ((sphov=check_spherocyl(CiCj, D, L, Di, Ci, ni, Dj, Cj, nj, &rim)) > 0.0)
    return 1;
#endif
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= L)
	    return -1;
	  else
	    return 1;
	}

      /* parallel disks */
      for (j1=0; j1 < 2; j1++)
	for (j2=j1; j2 < 2; j2++)
	  {
	    sp=0.0;
	    for (kk=0; kk < 3; kk++)
	      {
		VV[kk] = Di[j1][kk]-Dj[j2][kk];
		sp += ni[kk]*VV[kk];
	      }
	    if (sp == 0 && calc_norm(VV) < D)
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
	    if (normPiDi <= 0.5*D && normPjDj <= 0.5*D)
	      {
		Q1 = sqrt(Sqr(D)/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D)/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  /* case A.2 overlap of rim and disk */

  /* =================================== >>> Part A <<< ========================= */
  for (j1=0; j1 < 2; j1++)
    {
      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	    }
	}
      for (j2=0; j2 < 2; j2++)
	{
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);
#if 0
	  if (dostorebump)
	    {
	      printf("normDjUi=%.15G DjUini=%.15G\n", normDjUi, DjUini);
	      printf("Ci=%f %f %f Dj=%f %f %f\n", Ci[0], Ci[1], Ci[2], Dj[0], Dj[1], Dj[2]);
	      printf("DjUi=%.15G %.15G %.15G\n", DjUi[0], DjUi[1], DjUi[2]); 
	      printf("Uj=%.15G %.15G %.15G\n", Ui[0], Ui[1], Ui[2]); 
	      printf("nj=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
	      printf("DjCini= %.15G\n", DjCini);
	    }
#endif 
	  if (normDjUi > D)
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Ui| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

	   */
	  if (normDjUi < D*0.5 && fabs(DjCini) > L*0.5)
	    continue;

	  if (normDjUi < D*0.5 && fabs(DjCini) <= L*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }
#if 1
	  find_initial_guess(Ai, Ci, ni, Dj[j2], nj, D);

#else
	  for (kk=0; kk < 3; kk++)
	    {
	      //Ai[kk] = Ci[kk];
	      Ai[kk] = Ui[kk];  
	    }
#endif
	  for (it = 0; it < MAX_ITERATIONS; it++)
	    {
	      for (kk=0; kk < 3; kk++)
		{
		  AiDj[kk] = Ai[kk] - Dj[j2][kk];
		}
	      AiDjnj = scalProd(AiDj,nj);
	      vectProdVec(AiDj,nj,AiDjnjvec);
	      for (kk=0; kk < 3; kk++)
		VV[kk] =  0.5*D*(AiDj[kk]-AiDjnj*nj[kk])/calc_norm(AiDjnjvec);
	      for (kk=0; kk < 3; kk++)
		{
		  Tjp[kk] = Dj[j2][kk] + VV[kk];
		  Tjm[kk] = Dj[j2][kk] - VV[kk];
		  TjpCi[kk] = Tjp[kk] - Ci[kk];
		  TjmCi[kk] = Tjm[kk] - Ci[kk];
		}
	      TjpCini = scalProd(TjpCi,ni);  
	      TjmCini = scalProd(TjmCi,ni);
	      for (kk=0; kk < 3; kk++)
		{
		  Tjp_perp[kk] = TjpCi[kk]-TjpCini*ni[kk];
		  Tjp_para[kk] = TjpCini*ni[kk];
		  Tjm_perp[kk] = TjmCi[kk]-TjmCini*ni[kk];
		  Tjm_para[kk] = TjmCini*ni[kk];
		} 
	      normTjm_perp = calc_norm(Tjp_perp);
	      for (kk=0; kk < 3; kk++)
		TjOld[kk] = TjNew[kk];
	      if (calc_norm(Tjm_perp) < calc_norm(Tjp_perp))
		{
		  for (kk=0; kk < 3; kk++)
		    TjNew[kk] = Tjm[kk];
		}	  
	      else
		{
		  for (kk=0; kk < 3; kk++)
		    TjNew[kk] = Tjp[kk];
		}

	      for (kk=0; kk < 3; kk++)
		TjNewCi[kk] = TjNew[kk] - Ci[kk];
	      TjNewCini = scalProd(TjNewCi,ni);

#ifdef DEBUG_HCMC
	      printf("j1=%d A it=%d Aiold=%.15G %.15G %.15G\n", j1, it, Ai[0], Ai[1], Ai[2]);
#endif
	      for (kk=0; kk < 3; kk++)
		Ai[kk] = TjNewCini*ni[kk] + Ci[kk]; 
#ifdef DEBUG_HCMC
	      printf("A it=%d Ainew=%.15G %.15G %.15G TjNewCini=%.15G\n", it, Ai[0], Ai[1], Ai[2], TjNewCini);
	      printf("A Ci=%.15G %.15G %.15G\n", Ci[0], Ci[1], Ci[2]);
	      printf("A ni=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
#endif
	      if ( it > 0 && check_convergence(TjOld,TjNew) ) 
		break;
	    }
	  totitsHC += it;
#ifdef DEBUG_HCMC
	  printf("A #1 number of iterations=%d Tjold=%.15G %.15G %.15G Tjnew=%.15G %.15G %.15G\n",it, 
		 TjOld[0], TjOld[1], TjOld[2], TjNew[0], TjNew[1], TjNew[2]);
#endif
	  if (it >= MAX_ITERATIONS)
	    {
	      printf("MAX ITERATIONS REACHED in A!\n");
	      *retchk=1;
	      return -1;
	    }
	  if ( (calc_norm(Tjp_para) <= L*0.5 && calc_norm(Tjp_perp) <= D*0.5)||
	       (calc_norm(Tjm_para) <= L*0.5 && calc_norm(Tjm_perp) <= D*0.5) )
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #2 disk-rim\n");
#endif	   
	      return -1;
	    }
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	    }
	}

    }
  /* =================================== >>> Part B <<< ========================= */
#if 0
  for (j1=0; j1 < 2; j1++)
    {
      for (kk=0; kk < 3; kk++)
	DiCj[kk] = Di[j1][kk] - Cj[kk];
      normDiCj = calc_norm(DiCj);
      DiCjnj = scalProd(DiCj,nj);
      for (kk=0; kk < 3; kk++)
	{
	  Uj[kk] = Cj[kk] + DiCjnj*nj[kk];
	  DiUj[kk] = Di[j1][kk] - Uj[kk];
	}

      DiUjnj = scalProd(DiUj,nj);
      normDiUj = calc_norm(DiUj);
#ifdef DEBUG_HCMC
      if (dostorebump)
	{
	  printf("B normDiUj=%.15G DiUjnj=%.15G\n", normDiUj, DiUjnj);
	  printf("B Cj=%f %f %f Di=%f %f %f\n", Cj[0], Cj[1], Cj[2], Di[j1][0], Di[j1][1], Di[j1][2]);
	  printf("B DiUj=%.15G %.15G %.15G\n", DiUj[0], DiUj[1], DiUj[2]); 
	  printf("B Uj=%.15G %.15G %.15G\n", Uj[0], Uj[1], Uj[2]); 
	  printf("B nj=%.15G %.15G %.15G\n", nj[0], nj[1], nj[2]);
	  printf("DjCini= %.15G\n", DjCini);
	}
#endif 
 
      if (normDiUj > D)
	continue;

      if (normDiUj < D*0.5 && fabs(DiCjnj) > L*0.5)
	continue;

      if (normDiUj < D*0.5 && fabs(DiCjnj) <= L*0.5)
	{
#ifdef DEBUG_HCMC
	  if (dostorebump)
	    printf("B #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	  return -1;
	}      
      for (kk=0; kk < 3; kk++)
	{
	  Aj[kk] = Cj[kk];
	}
      for (it = 0; it < MAX_ITERATIONS; it++)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      AjDi[kk] = Aj[kk] - Di[j1][kk];
	    }
	  AjDini = scalProd(AjDi,ni);
	  vectProdVec(AjDi,ni,AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    VV[kk] =  0.5*D*(AjDi[kk]-AjDini*ni[kk])/calc_norm(AjDinivec);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tip[kk] = Di[j1][kk] + VV[kk];
	      Tim[kk] = Di[j1][kk] - VV[kk];
	      TipCj[kk] = Tip[kk] - Cj[kk];
	      TimCj[kk] = Tim[kk] - Cj[kk];
	    }
	  TipCjnj = scalProd(TipCj,nj);  
	  TimCjnj = scalProd(TimCj,nj);
	  for (kk=0; kk < 3; kk++)
	    {
	      Tip_perp[kk] = TipCj[kk]-TipCjnj*nj[kk];
	      Tip_para[kk] = TipCjnj*nj[kk];
	      Tim_perp[kk] = TimCj[kk]-TimCjnj*nj[kk];
	      Tim_para[kk] = TimCjnj*nj[kk];
	    } 
	  normTim_perp = calc_norm(Tip_perp);
	  for (kk=0; kk < 3; kk++)
	    TiOld[kk] = TiNew[kk];
	  if (calc_norm(Tim_perp) < calc_norm(Tip_perp))
	    {
	      for (kk=0; kk < 3; kk++)
		TiNew[kk] = Tim[kk];
	    }	  
	  else
	    {
	      for (kk=0; kk < 3; kk++)
		TiNew[kk] = Tip[kk];
	    }

	  for (kk=0; kk < 3; kk++)
	    TiNewCj[kk] = TiNew[kk] - Cj[kk];
	  TiNewCjnj = scalProd(TiNewCj,nj);
#ifdef DEBUG_HCMC
	  printf("B it=%d Ajold=%.15G %.15G %.15G\n", it, Aj[0], Aj[1], Aj[2]);
#endif
	  for (kk=0; kk < 3; kk++)
	    Aj[kk] = TiNewCjnj*nj[kk] + Cj[kk]; 
#ifdef DEBUG_HCMC
	  printf("B it=%d Ajnew=%.15G %.15G %.15G TiNewCjnj=%.15G\n", it, Aj[0], Aj[1], Aj[2], TiNewCjnj);
	  printf("B Ci=%.15G %.15G %.15G\n", Cj[0], Cj[1], Cj[2]);
	  printf("B ni=%.15G %.15G %.15G\n", nj[0], nj[1], nj[2]);
	  printf("B #1 number of iterations=%d Tiold=%.15G %.15G %.15G Tinew=%.15G %.15G %.15G\n",it, 
	     TiOld[0], TiOld[1], TiOld[2], TiNew[0], TiNew[1], TiNew[2]);

#endif
	
	  if ( it > 0 && check_convergence(TiOld,TiNew) ) 
	    {
	      break;
	    }
	} 
      totitsHC += it;
#ifdef DEBUG_HCMC
      printf("B #1 number of iterations=%d Tiold=%.15G %.15G %.15G Tinew=%.15G %.15G %.15G\n",it, 
	     TiOld[0], TiOld[1], TiOld[2], TiNew[0], TiNew[1], TiNew[2]);
#endif
 
      if (it >= MAX_ITERATIONS)
       	{
 	  printf("MAX ITERATIONS REACHED IN B\n");
	  *retchk=1;
#ifdef DEBUG_HCMC
	  //exit(-1);
#endif
 	  return -1;
  	}
      
     // printf("#2 number of iterations=%d\n",it);
      if ( (calc_norm(Tip_para) <= L*0.5 && calc_norm(Tip_perp) <= D*0.5)||
	   (calc_norm(Tim_para) <= L*0.5 && calc_norm(Tim_perp) <= D*0.5) )
	{
#ifdef DEBUG_HCMC
	  if (dostorebump)
	    printf("B #2 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	  return -1;
	}
    }
#endif
  numcallsHC += 4.0; 

  /* case A.3 rim-rim overlap */
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < D && fabs(lambdai) < 0.5*L && fabs(lambdaj) < 0.5*L)
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 1;
}
struct cylstr 
{
  double u[3];
  double r[3];
  double L;
  double D;
  int i1;
  int i2;
}
*cylinders;
int allocated_cyls=100;
int numcyls=0;
void add_cylinder(double r[], doulbe u[], double L, int i, int j)
{
  int k;
  numcyls++; 
  if (numcyls >= allocated_cyls)
    {
      allocated_cyls += 100;
      cylinders = realloc(cylinders, allocated_cyls);
    }
  for (k=0; k < 3; k++)
    {
      cylinders[numcyls-1].r[k] = r[k];
      cylinders[numcyls-1].u[k] = u[k];  
      cylinders[numcyls-1].i1 = i;
      cylinders[numcyls-1].i2 = j;
      cylinders[numcyls-1].D = wellWidth;
      cylinders[numcyls-1].L = L;
    }
}
int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int  NN=-1, fine, JJ, nat, maxl, maxnp, np, nc2, nc, dix, diy, diz, djx,djy,djz,imgi2, imgj2, jbeg, ifin;
  int jX, jY, jZ, iX, iY, iZ, iold;
  //int coppie;
  double refTime=0.0, ti, ene=0.0, ucyl[3], rcm[3], r0old[3];
  int curcolor, ncls, b, j, almenouno, na, c, i2, j2, ncls2;
  wellWidth=-1.0;
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
  fscanf(f,"%d %d\n", &block, &dummyint);
  //fscanf(f,"%[^\n]\n", line);
  fscanf(f,"%d\n", &NP);
  fscanf(f,"%lf %lf %lf ", &Lx, &Ly, &Lz);
  fscanf(f,"%[^\n]\n", line);
  fscanf(f,"%[^\n]\n", line);
  fclose(f);

  L=Lx;
#if 0 
  f = fopen("mols.dat", "r");
  fscanf(f, "%s %d ", line, &NP);
  fscanf(f, "%s %d ", line, &NPA);
  fclose(f);
#endif
  f = fopen("sigma.dat", "r");
  fscanf(f, "%lf %lf ", &sigmaAA, &sigmaBB);
  sigmaAB=0.5*(sigmaAA+sigmaBB);
  fclose(f);

  cylinders = malloc(sizeof(struct cylstr)*allocated_cyls);
 
  printf("NP=%d NPA=%d sigmaBB=%f\n", NP, NPA, sigmaBB);  
  NPA=0;
  mix_type=-1;
  if (mix_type==-1 || NPA==NP)
    {
    	START=0;
        END=NP;
    }
   else if (mix_type==0)
    {
      START=0;
      END=NPA;
    }	 
  else
    {
      START=NPA;
      END=NP;
    }
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
  clssizedstAVG = malloc(sizeof(double)*NP);
  dupcluster = malloc(sizeof(int)*NP*NUMREP); 
  percola = malloc(sizeof(int)*NP);
  ip = malloc(sizeof(int)*NP);
  if (output_bonds)
    {
      numbonds  = malloc(sizeof(int)*NP); 
      bonds = AllocMatI(NP, MAXBONDS);
    }
  if (wellWidth==-1)
    wellWidth=sigmaBB;
  maxsaxAA = fabs(sigmaAA);
  maxsaxAB = fabs(sigmaAB);
  maxsaxBB = fabs(wellWidth);
    /* le AA sono le grandi quindi usiamo quelle per RCUT */
  RCUT = maxsaxBB*1.01;
  printf("maxsaxBB=%f RCUT=%f\n", maxsaxBB, RCUT);	
  for (a = 0; a < 3; a++)
    {
#if 0
      for (b = 0; b < NA; b++)
	rat[b][a] = malloc(sizeof(double)*NP);
#endif
      r0[a] = malloc(sizeof(double)*NP);
      //DR0[a] = malloc(sizeof(double)*NP);
#if 0
      for (b = 0; b < 3; b++)
	{
	  R[a][b] = malloc(sizeof(double)*NP);
	}
#endif
    }
  printf("[MIXTURE] files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", nfiles, NP, NPA, L, NN, maxl);
  //printf("sigmaSticky=%.15G\n", sigmaSticky);
  for (i = 0; i < NP; i++)
    {
      clssizedstAVG[i] = 0.0;
    }      
  for (nr1 = 0; nr1 < nfiles; nr1++)
    {	
      if (cellList)
	{
	  free(cellList);
	  free(inCell[0]);
	  free(inCell[1]);
	  free(inCell[2]);
	}
     readconf(fname[nr1], &time, &refTime, NP, r0);
      /* =================== >>> CREATE CYLINDERS <<< ===============
	 se la stessa particella appartiene a n frame allora il cluster di n particelle
	 corrispondente conta 1 */
      for (i=0; i < NP/block; i++)
	{
	  /* considera tutti i frame */
	  for (j=i; j < NP; j=j+block)
	    {
	      if (ip[j] != i)
		{
		  printf("boh \n");
		  exit(-1);
		}
	      if (j > i)
		{
		  for (k=0; k < 3; k++)
		    ucyl[0] = r0[k][j] - r0old[k];
		  norm = calc_norm(ucyl);
		  for (k=0; k < 3; k++)
		    ucyl[k] /= norm;
		  for (k=0; k < 3; k++)
		    rcm[k] = (r0[k][j] + r0old[k])*0.5;
		  add_cylinder(rcm, ucyl, norm, iold, i);
		}
	      for (k=0; k < 3; k++)
		r0old[k] = r0[k][i];
	      iold = i;
	    }
	}

      printf("created %d cylinders\n", numcyls);
      NP = numcyls;
      ti = time + refTime;
      for (i = 0; i < NP; i++)
	{
	  percola[i] = 0; 
	}
      for (i=0; i < NP; i++)
	{

	}
      cellsx = L / RCUT;
      cellsy = L / RCUT;
      cellsz = L / RCUT;
      cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+NP));
      inCell[0] = malloc(sizeof(int)*NP);
      inCell[1] = malloc(sizeof(int)*NP);
      inCell[2] = malloc(sizeof(int)*NP);
       /* costruisce la posizione di tutti gli sticky spots */
      for (i = 0; i < NP; i++)
	{
	  color[i] = -1;	  
	  clssizedst[i] = 0;
	}
      curcolor = 0;
      ene=0;
      //coppie = 0;
      build_linked_list();

      jbeg = 0; 
      ifin = NP;
      for (i = START; i < END; i++)
	{
    	  if (color[i] == -1)
	    color[i] = curcolor;
	    
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
		      	  if (j <= i) 
	    		    continue;
			  if (bond_found(i, j))  
			    {
			      if (output_bonds)
				{
				  add_bond(i, j);
				  add_bond(j, i);
				}
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
      for (i = START; i < END; i++)
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

      /* stima volumi cluster */
      for (nc = 0; nc < ncls; nc++)
	{

	}	
#if 0
      /* ============== >>> PERCOLATION <<< ================== */
      if (check_percolation)
	{
	  free(cellList);
	  free(inCell[0]);
	  free(inCell[1]);
	  free(inCell[2]);
	  cellsx = 2.0*L / RCUT;
	  cellsy = 2.0*L / RCUT;
	  cellsz = 2.0*L / RCUT;
	  cellList = malloc(sizeof(int)*(cellsx*cellsy*cellsz+NP*NUMREP));
	  inCell[0] = malloc(sizeof(int)*NP*NUMREP);
	  inCell[1] = malloc(sizeof(int)*NP*NUMREP);
	  inCell[2] = malloc(sizeof(int)*NP*NUMREP);

	  ene=0;
	  //coppie = 0;
	  for (nc = ncls-1; nc >= 0; nc--)
	    {
	      if (cluster_sort[nc].dim==1)
		continue;

	      //printf("Analysing cluster #%d of #%d\n", nc+1, ncls);
	      /* N.B per verificare la percolazione ogni cluster va "duplicato"
	       * in tutte le direzioni e se alla fine risulta comunque un unico 
	       * cluster allora tale cluster è percolante.*/
	      na = 0;
	      //printf("i=1011 j=277 rat=%.15G %.15G\n", rat[0][0][1011], rat[0][0][377]);
	      for (i=START; i < END; i++)
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

	      build_linked_list_perc(cluster_sort[nc].dim, 2.0*L);
	      //printf("i=1011 j=277 rat=%.15G %.15G\n", rat[0][0][1011], rat[0][0][377]);
	      //printf("NP=%d NPA=%d na=%d,clsdim[%d]=%d\n", NP, NPA,na,nc,cluster_sort[nc].dim);
	      curcolor = 0;
	      for (i2 = 0; i2 < na*NUMREP; i2++)
		{
		  color2[i2] = -1;	  
		}
	      //printf("na=%d clsdim=%d\n", na, cluster_sort[nc].dim);
	      for (i2 = 0; i2 < na*NUMREP; i2++)
		{
		  if (color2[i2]==-1)
		    color2[i2] = curcolor;
		  //printf("nc=%d na*NUMREP=%d i2=%d\n",nc, na*NUMREP, i2);
		  //printf("curcolor:%d\n", curcolor);
		  i = dupcluster[i2];
		  for (iZ = -1; iZ <= 1; iZ++) 
		    {
		      jZ = inCell[2][i2] + iZ;    
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
			  jY = inCell[1][i2] + iY;    
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
			      jX = inCell[0][i2] + iX;    
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
			      j2 = (jZ *cellsy + jY) * cellsx + jX + NP*NUMREP;
			      for (j2 = cellList[j2]; j2 > -1; j2 = cellList[j2]) 
				{
				  //if (color[j] != cluster_sort[nc].color)
				    //continue;
				  //coppie++;
				  j = dupcluster[j2];
			 	  if (j2 <= i2) 
				    continue;
				  //dix = diy = diz = 0;
				  //djx = djy = djz = 0;
				  imgi2 = i2 / na;
				  imgj2 = j2 / na;
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
			    }
			}
		    }
#if 0
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
#endif
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
		{
		  percola[nc] = 1;
		  printf("#clusters in the replicated system: %d (of %d replicas)\n", ncNV2, NUMREP);
		  break;
		}
	    }
	  printf("E/N (PERCOLATION) = %.15G\n", ene/((double)(NUMREP))/((double)NP));
	}
#endif
      //printf("coppie PERC=%d\n", coppie);
      almenouno = 0;

      sprintf(fn, "perc%s.dat", fname[nr1]);
      f2 = fopen(fn, "w");
      fclose(f2);
      for (nc = 0; nc < ncls; nc++)
	{
	  //if (cluster_sort[nc].dim >= 2)
	    //almenouno = 1;
	  if (percola[nc])
	    {
	      sprintf(fn, "perc%s.dat", fname[nr1]);
	      f2 = fopen(fn, "a");
	      fprintf(f2, "%d %d\n", nc, cluster_sort[nc].dim);
	      fclose(f2);
	    }
	  if (percola[nc])
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
	  clssizedstAVG[cluster_sort[nc].dim] += 1.0;
	}
      sprintf(fncls, "%s.clsdst", fname[nr1]);
      f = fopen(fncls, "w+");
      for (i = 1; i < NP; i++)
	{
	  if (clssizedst[i] != 0)
	    fprintf(f, "%d %d\n", i, clssizedst[i]);
	}
      fclose(f);
      if (output_bonds)
	{
	  sprintf(fn, "%s.bonds", fname[nr1]);
	  f = fopen(fn, "w+");
	  fprintf(f, "%d %.15G\n", START-END, L);
	  for (i = START; i < END; i++)
	    {
	      fprintf(f,"%.15G %.15G %.15G\n", r0[0][i], r0[1][i], r0[2][i]);
	    }	  
	  for (i = START; i < END; i++)
	    {
	      fprintf(f,"%d %d\n", i+1, numbonds[i]);
	      for (c = 0; c < numbonds[i]-1; c++)
    		fprintf(f, "%d ", bonds[i][c]+1);
	      fprintf(f, "%d\n", bonds[i][numbonds[i]-1]+1);
	    }
	  fclose(f);
	}
    }

  f = fopen("avg_cluster_size_distr.dat", "w+");
  if (media_log)
    {
      for (kk=1; kk <= 51; kk++)
	l1[kk]=(int) nlin*pow(1.25,kk-1);

      for(kk=1; kk <= 50; kk++)
	{
	  l2[kk]=l1[kk+1]-1;
	  if (l2[kk] < npmax) 
	    kmax=kk;
	}
      for(kk=1; kk <= kmax; kk++)
	{
	  dlog[kk]=0.0;
	  xlog[kk]=0.0;
	  for (kj=l1[kk]; kj <= l2[kk]; kj++)
	    {
	      if (kj < NP && clssizedstAVG[kj] !=0)
		{   
		  dlog[kk]=dlog[kk]+clssizedstAVG[kj];
		}		  
	      xlog[kk]=xlog[kk]+kj;
	    }
	}
      for (i3=0; i3 < nlin; i3++)
	{
	  if (clssizedstAVG[i3] != 0) fprintf(f,"%d %.15G\n", i3,((double)clssizedstAVG[i3])/((double)(nfiles)));
	}
      for (kk=1; kk <= kmax; kk++)
	{
	  if (dlog[kk]!=0) 
	    {
	      am=l2[kk]-l1[kk]+1; 
	      xmed=xlog[kk]/am;
	      dlog[kk]=dlog[kk]/am; 
	      fprintf(f,"%.15G %.15G\n", xmed,((double) dlog[kk])/((double)nfiles));
	    }
	}
    }
  else
    {
      for (i = 1; i < NP; i++)
	{
	  if (clssizedstAVG[i] != 0.0)
	    fprintf(f, "%d %.15G\n", i, ((double)clssizedstAVG[i])/((double)nfiles));
	}
    }
  fclose(f);
  return 0;
}
