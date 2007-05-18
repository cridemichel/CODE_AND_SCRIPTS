#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 5000
#define Sqr(x) ((x)*(x))
int NP, NPA=-1;
char **fname; 
double L, time, *ti, *rotMSD, *MSD, *rotMSDA, *MSDA, *rotMSDB, *MSDB, *cc, *rotMSDcls[2], *MSDcls[2], *cc_cls[2];
double **DR, deltat, ***R;
double *r0[3], *w0[3], *rt[3], *wt[3], *rtold[3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int points=-1, foundDRs=0, foundrot=0, eventDriven=0, skip=1, clusters=0;
int *isPercPart, curStep;
char *pnum;
char inputfile[2048], cluststr[2048];
double storerate = -1;
int bakSaveMode = -1;
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
double **matrix(int n, int m)
{
  double **M;
  int i;
  M = malloc(sizeof(double*)*n);
  for (i=0; i < n; i++)
    M[i] = calloc(m, sizeof(double));
  return M;
}

#define MD_SILICA
#define MD_THREESPOTS
void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3])
{
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  double r1[3], r2[3], r3[3], nr;
  double radius; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */
  /* WARNING: se le particelle che interagiscono hanno diametro diverso da sigma[0][1] qui va cambiato!!!! */ 
#ifdef MD_AB41
  if (i < Oparams.parnumA)
    radius = Oparams.sigma[0][0] / 2.0;
  else
    radius = Oparams.sigma[1][1] / 2.0;
  //printf("i=%d sigma=%f\n", i, radius*2.0);
#else
  radius =  0.5;
#endif
  if (ata == 0)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk];
      //printf("%f %f %f @ 0.5 C[red]\n", rat[0], rat[1], rat[2]);
    }
  else if (ata <= 3)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk] + R[kk][ata-1]; 
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
      //printf("ata=%d %f %f %f @ 0.075 C[blue]\n", ata, R[0][ata-1], R[1][ata-1], R[2][ata-1]);
    }
  else
    {
      /* l'atomo restante è un electron site */
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
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
    }

}
void BuildAtomPos(int i, double *rO, double **R, double rat[5][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a, NUMAT;
  /* l'atomo zero si suppone nell'origine */
#ifdef MD_SILICA
#ifdef MD_THREESPOTS
  if (i >= NPA)
    NUMAT = 4;
  else
    NUMAT = 3;
#elif defined(MD_AB41)
  if (i >= Oparams.parnumA)
    NUMAT = 2;
  else
    NUMAT = 5;
#else
  if (i >= Oparams.parnumA)
    NUMAT = 5;
  else
    NUMAT = 3;
#endif
#else
  NUMAT = 5;
#endif
  for (a=0; a < NUMAT; a++)
    BuildAtomPosAt(i, a, rO, R, rat[a]);
}
double **Rt;
#define NA 3
void save_fra(char *fname)
{
  char fileop[1024], fileop2[1024], fileop3[1024];
  double rat[NA][3];
  double rcm[3];	
  int i, k1, k2;
  FILE* f;
  sprintf(fileop2 ,"%s.fra", fname);
  /* store conf */
  f = fopen(fileop2,"w");
  fprintf(f, "%d 0 %d %d 0\n", curStep, NP, NP-NPA);
  fprintf(f, "%.15G %.15G %.15G 0 0 %.15G\n", L, L, L, deltat);
  for (i = NPA; i < NP; i++)
    {
      //printf("i=%d\n",i);
      rcm[0] = r0[0][i];
      rcm[1] = r0[1][i];
      rcm[2] = r0[2][i];
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  Rt[k1][k2] = R[i][k1][k2];
      BuildAtomPos(i, rcm, Rt, rat);
      fprintf(f, "%.15G %.15G %.15G\n", rat[1][0], rat[1][1], rat[1][2]);
      fprintf(f, "%.15G %.15G %.15G\n", rat[2][0], rat[2][1], rat[2][2]);
      fprintf(f, "%.15G %.15G %.15G\n", rat[0][0], rat[0][1], rat[0][2]);
    }
  for (i = 0; i < NPA; i++)
    {
      rcm[0] = r0[0][i];
      rcm[1] = r0[1][i];
      rcm[2] = r0[2][i];
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  Rt[k1][k2] = R[i][k1][k2];
      BuildAtomPos(i, rcm, Rt, rat);
      fprintf(f, "%.15G %.15G %.15G\n", rat[1][0], rat[1][1], rat[1][2]);
      fprintf(f, "%.15G %.15G %.15G\n", rat[2][0], rat[2][1], rat[2][2]);
      fprintf(f, "%.15G %.15G %.15G\n", rat[0][0], rat[0][1], rat[0][2]);
    }
#if 0
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      //printf("i=%d\n",i);
      fprintf(f, "%.15G %.15G %.15G\n", v0[1][0], v0[1][1], v0[1][2]);
    }
  for (i = 0; i < Oparams.parnumA; i++)
    {
      fprintf(f, "%.15G %.15G %.15G\n", v0[1][0], v0[1][1], v0[1][2]);
    }
#endif
  fclose(f);
}

void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double **DR, 
	      double ***R)
{
  FILE *f;
  int nat=0, i, cpos;
  double dt=-1;
  int curstp=-1;

  *ti = -1.0;
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
		  fscanf(f, " %lf %lf %lf ", &DR[i][0], &DR[i][1], &DR[i][2]);
		}
	      foundDRs = 1;
	    }
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
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else if (!strcmp(parname, "curStep"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      curstp = atoi(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else if (!strcmp(parname, "steplength"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      dt = atof(parval);
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
	      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
    		     &r[0][i], &r[1][i], &r[2][i], &R[i][0][0], 
		     &R[i][0][1], &R[i][0][2], &R[i][1][0], 
		     &R[i][1][1], &R[i][1][2], &R[i][2][0], 
		     &R[i][2][1], &R[i][2][2]); 
	    }
	  break; 
	}
    }
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
  if (*ti == -1)
    *ti = ((double)curstp)*dt;
  fclose(f);
}
void print_usage(void)
{
  printf("Usage: cri2fra <lista_file_da_convertire>\n");
  exit(0);
}
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
	  points = atoi(argv[cc]);
	}
      else
	print_usage();
      cc++;
    }
}

int main(int argc, char **argv)
{
  FILE *f, *f2, *f3, *fA, *fB, *f2A, *f2B;
  double *adjDr[3], Dr, Dw, A1, A2, A3, dr;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int NN=-1, fine, JJ, nat, maxl, maxnp, np, NP1, NP2;
  double refTime=0.0, tmpdbl;
  int isperc, kk;
  parse_param(argc, argv);
  f2 = fopen(inputfile, "r");
  c2 = 0;
  maxl = 0;
  Rt = malloc(sizeof(double*)*3);
  for (kk=0; kk < 3; kk++)
    Rt[kk] = malloc(sizeof(double)*3);
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
      if (!strcmp(parname,"storerate"))
	{
	  eventDriven = 1;
	  storerate = atof(parval);
	}
      if (!strcmp(parname,"curStep"))
	{
	  curStep = atoi(parval);
	}
      if (!strcmp(parname,"Dt"))
	{
	  deltat = atof(parval);
	}
      if (!strcmp(parname,"bakSavedMode"))
	bakSaveMode = atoi(parval);
      if (!strcmp(parname,"parnum"))
	{
	  if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
	    {
	      NP = atoi(parval);
	    }
	  else
	    {
	      NP = NP1+NP2;
	      NPA = NP1;
	    }
	}
      if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      if (nat==0 && !strcmp(parname,"NN"))
	NN = atoi(parval);
    }
  fclose(f);
  if (NPA == -1)
    NPA = NP;
  if (points == -1)
    points=NN;
  if ((eventDriven==1 && storerate <= 0.0 && bakSaveMode <= 0)
      || (eventDriven==0 && bakSaveMode <= 0)) 
    NN = 1;
  if (NN!=1)
    skip = 0;
  
  if (eventDriven==0)
    L = cbrt(L);
  if (NPA != NP)
    printf("[MIXTURE] points=%d files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", points, nfiles, NP, NPA, L, NN, maxl);
  else
    printf("[MONODISPERSE] points=%d files=%d NP = %d L=%.15G NN=%d maxl=%d\n", points, nfiles, NP, L, NN, maxl);
  if (eventDriven)
    printf("[ED] Event-Driven simulation\n");
  else
    printf("[MD] Time-Driven simulation\n");
  DR = malloc(sizeof(double*)*NP);
  for (ii = 0; ii < NP; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*3);
    }

  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      w0[a] = malloc(sizeof(double)*NP);
    }
  R = malloc(sizeof(double**)*NP);
  for (i=0; i < NP; i++)
    {
      R[i] = matrix(3, 3);
    }

  for (nr1 = 0; nr1 < nfiles; nr1++)
    {	
      readconf(fname[nr1], &time, &refTime, NP, r0, w0, DR, R);
      save_fra(fname[nr1]);
    }
  return 0;
}
