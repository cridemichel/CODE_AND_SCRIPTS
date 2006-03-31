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
double L, time, *ti, *R[3][3], *cc, *r0[3], r0L[3], RL[3][3], *DR0[3];
double pi, sa[2], sb[2], sc[2], Dr, theta, sigmaSticky, ratL[NA][3], *rat[NA][3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int NP, NPA=-1;
int points, foundDRs=0, foundrot=0, *color, *clsdim, *clsdimsort, *clssizedst, *clssizedstAVG, *percola;
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
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
	      if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
		{
		  sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], 
			 &R[0][0][i], &R[0][1][i], &R[0][2][i], &R[1][0][i], &R[1][1][i], &R[1][2][i],
			 &R[2][0][i], &R[2][1][i], &R[2][2][i]); 
		}
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
void BuildAtomPos(int i, double rO[3], double R[3][3], double rat[NA][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a;
  /* l'atomo zero si suppone nell'origine */
  if (i < NP)
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
double distance(int i, int j, int img)
{
  int a, b;

  for (a = 0; a < MD_STSPOTS_A; a++)
    {
      for (b = 0; b < MD_STSPOTS_B; b++)
	{
	  if (Sqr(rat[a][0][i] + img*L -rat[a][0][j])+Sqr(rat[a][1][i] + img*L -rat[a][1][j])
	      +Sqr(rat[a][2][i] + img*L -rat[a][2][j]) < Sqr(sigmaSticky))	  
		return -1;
	}
    }
  return 1;
}

int bond_found(int i, int j, int img)
{
  if (distance(i, j, img) < 0.0)
    return 1;
  else
    return 0;
}
void change_all_colors(int colorsrc, int colordst)
{
  int ii;
  for (ii = 0; ii < NP; ii++)
    {
      if (color[ii] == colorsrc)
	color[ii] = colordst;
    }
}
char fncls[1024];
int findmaxColor(int *color)
{
  int i, maxc=-1;
  for (i = 0; i < NP; i++) 
    {
      if (color[i] > maxc)
	maxc = color[i];
    }
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

int main(int argc, char **argv)
{
  FILE *f, *f2, *f3;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int  NN, fine, JJ, nat, maxl, maxnp, np, nc;
  double refTime=0.0, ti;
  int curcolor, ncls, b, j;
  pi = acos(0.0)*2.0;
  if (argc <= 1)
    {
      printf("Usage: clusters <listafile>\n");
      //printf("where NN il the lenght of the logarithmic block\n");
      exit(-1);
    }
  f2 = fopen(argv[1], "r");
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
      else if (nat == 1 && !strcmp(parname,"NN"))
	NN = atoi(parval);
      else if (nat==1 && !strcmp(parname,"a"))
	sscanf(parval, "%lf %lf\n", &sa[0], &sa[1]);	
      else if (nat==1 && !strcmp(parname,"b"))
	sscanf(parval, "%lf %lf\n", &sb[0], &sb[1]);	
      else if (nat==1 && !strcmp(parname,"c"))
	sscanf(parval, "%lf %lf\n", &sc[0], &sc[1]);	
      else if (nat==1 && !strcmp(parname,"sigmaSticky"))
	sigmaSticky = atof(parval);
      else if (nat==1 && !strcmp(parname,"theta"))
	theta = atof(parval);
      else if (nat==1 && !strcmp(parname,"Dr"))
	Dr = atof(parval);
    }
  fclose(f);
  if (NPA == -1)
    NPA = NP;
  if (argc == 3)
    points=atoi(argv[2]);
  else
    points=NN;
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  //ti = malloc(sizeof(double)*points);
  cc = malloc(sizeof(double)*points);
  color = malloc(sizeof(int)*NP);
  clsdim = malloc(sizeof(int)*NP);
  cluster_sort = malloc(sizeof(struct cluster_sort_struct)*NP);
  clssizedst = malloc(sizeof(int)*NP);
  clssizedstAVG = malloc(sizeof(int)*NP);
  percola = malloc(sizeof(int)*NP);
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
  for (ii=0; ii < points; ii++)
    {
      cc[ii]=0.0;
    }
  build_atom_positions();

  if (NPA != NP)
    printf("[MIXTURE] files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", nfiles, NP, NPA, L, NN, maxl);
  else
    printf("[MONODISPERE] files=%d NP = %d L=%.15G NN=%d maxl=%d\n", nfiles, NP, L, NN, maxl);
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
	  BuildAtomPos(i, r0L, RL, ratL);
	  for (a = 0; a < NA; a++)
	    for (b = 0; b < 3; b++)
	      rat[a][b][i] = ratL[a][b];
	}
      for (i = 0; i < NP; i++)
	{
	  color[i] = -1;	  
	  clssizedst[i] = 0;
	}
      curcolor = 0;
      for (i = 0; i < NPA; i++)
	{
	  color[i] = curcolor;
	  for (j = NPA; j < NP; j++)
	    {
      	      if (bond_found(i, j, 0))
		{
		  if (color[j] == -1)
		    color[j] = color[i];
		  else
		    {
		      if (color[i] < color[j])
			change_all_colors(color[j], color[i]);
		      else if (color[i] > color[j])
			change_all_colors(color[i], color[j]);
		    }
		}
	      if (bond_found(i, j, -1) || bond_found(i, j, +1))
		percola[color[i]] = 1;

	    }
	  curcolor = findmaxColor(color)+1;
	}	  
      sprintf(fncls, "%s.clusters", fname[nr1]);
      f = fopen(fncls, "w+");
      ncls = curcolor;
      for (nc = 0; nc < ncls; nc++)
	{
	  for (a = 0; a < NP; a++)
	    if (color[a] == nc)
	      clsdim[color[a]]++;
	}
      for (nc = 0; nc < ncls; nc++)
	{
	  cluster_sort[nc].dim = clsdim[nc];
	  cluster_sort[nc].color = nc;
	}
      qsort(cluster_sort, ncls, sizeof(struct cluster_sort_struct), compare_func);
      for (nc = 0; nc < ncls; nc++)
	{
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
      fclose(f);
      for (nc = 0; nc < ncls; nc++)
	{
	  clssizedst[clsdim[nc]]++;
	  clssizedstAVG[clsdim[nc]]++;
	}
      sprintf(fncls, "%s.clsdst", fname[nr1]);
      f = fopen(fncls, "w+");
      for (i = 2; i < NP; i++)
	{
	  if (clssizedst[i] != 0)
	    fprintf(f, "%d %d\n", i, clssizedst[i]);
	}
      fclose(f);
    }
  f = fopen("avg_cluster_size_distr.dat", "w+");
  for (i = 2; i < NP; i++)
    {
      if (clssizedstAVG[i] != 0)
	fprintf(f, "%d %d\n", i, clssizedstAVG[i]/nfiles);
    }
  fclose(f);
  return 0;
}
