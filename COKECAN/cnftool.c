#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define NA 10
#define MAXFILES 5000
#define Sqr(x) ((x)*(x))
int NP, NPA=-1, MCsim=0, dummyint, reorder=-1, rempar=0, partial=0;
char **fname; 
double sax, say, saz;
double X0, L, Lx, Ly, Lz, time, *ti, *rotMSD, *MSD, *rotMSDA, *MSDA, *rotMSDB, *MSDB, *cc, *rotMSDcls[2], *MSDcls[2], *cc_cls[2];
double **DR, deltat, ***R;
double *r0[3], *w0[3], *rt[3], *wt[3], *rtold[3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int partType, points=-1, foundDRs=0, foundrot=0, eventDriven=0, skip=1, clusters=0;
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

void BuildAtomPosAt(int i, int ata, double *rO, double **R, double rat[3])
{
  /* calcola le coordinate nel laboratorio di uno specifico atomo */
  int kk;
  double r1[3], r2[3], r3[3], nr;
  double radius, spdist; 
  /* l'atomo zero si suppone nell'origine 
   * la matrice di orientazione ha per vettori colonna le coordinate nel riferimento
   * del corpo rigido di tre sticky point. Il quarto sticky point viene ricostruito
   * a partire da questi. */
  /* WARNING: se le particelle che interagiscono hanno diametro diverso da sigma[0][1] qui va cambiato!!!! */ 
  spdist = 0.5;/*X0-(2.0-1.541665);*/
  radius=1.21667;
  if (ata == 0)
    {
      for (kk = 0; kk < 3; kk++)
	rat[kk] = rO[kk];
      //printf("%f %f %f @ 0.5 C[red]\n", rat[0], rat[1], rat[2]);
    }
  else if (ata <= 3)
    {
      for (kk=0; kk < 3; kk++) 
	rat[kk] = rO[kk] + ((ata==1)?+1.0:-1.0)*spdist*R[0][kk]; 
      //printf("%f %f %f @ 0.075 C[blue]\n", rat[0], rat[1], rat[2]);
      //printf("ata=%d %f %f %f @ 0.075 C[blue]\n", ata, R[0][ata-1], R[1][ata-1], R[2][ata-1]);
    }
}
void BuildAtomPos(int i, double *rO, double **R, double rat[NA][3])
{
  /* calcola le posizioni nel laboratorio di tutti gli atomi della molecola data */
  int a, NUMAT;
  /* l'atomo zero si suppone nell'origine */
  NUMAT=3;
  for (a=0; a < NUMAT; a++)
    BuildAtomPosAt(i, a, rO, R, rat[a]);
}
double **Rt;
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
  fprintf(f, "%d 0 %d %d 2\n", curStep, NP, NP);
  fprintf(f, "%.15G %.15G %.15G 0 0 0\n", Lx, Ly, Lz);
  if (MCsim)
    {
      for (i = 0; i < NP; i++)
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

    }
  else
    {
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
void outconf(char *fname)
{
  FILE *f;
  int nat=0, i, cpos;
  double dt=-1;
  int curstp=-1;

  f = fopen(fname, "r");
  while (!feof(f)) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      //fprintf(stderr, "line=%s nat=%d\n", line, nat);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  fprintf(stdout, "%s\n", line);
	  //cpos =ftell(f);
	  if (nat==1)
	    {
	      fscanf(f, "%[^\n]\n",line);
	      fprintf(stdout, "%d\n", NP-rempar);
	    }
	}
      else 
	{
	  if (nat == 0)
	    {
	      sscanf(line, "%[^:]:%s\n", parname, parval);
	      //printf("[%s] parname=%s\n", fname, parname);
	      if (!strcmp(parname,"parnum"))
		{
		  fprintf(stdout, "parnum: %d\n", NP-rempar);
		}
	      else
	    	fprintf(stdout, "%s\n",line);
	    }
	  // printf("nat=%d line=%s\n", nat, line);
	  else if (nat==2)
	    {
	      fseek(f, cpos, SEEK_SET);
	      for (i = 0; i < NP; i++) 
		{
		  fscanf(f, "%[^\n]\n", line); 
		  if (i < rempar)
		    continue;
		  else
		    fprintf(stdout, "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %d\n",
			    r0[0][i], r0[1][i], r0[2][i], R[i][0][0], 
			    R[i][0][1], R[i][0][2], R[i][1][0], 
			    R[i][1][1], R[i][1][2], R[i][2][0], 
			    R[i][2][1], R[i][2][2], 0);
		}		
	      fscanf(f, "%[^\n]\n", line); 
	      fprintf(stdout, "%s\n", line);
	    }
	  else
	    fprintf(stdout, "%s\n", line);
	}
    }
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
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
  while (!feof(f)) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  cpos =ftell(f);
	}
     // printf("nat=%d line=%s\n", nat, line);
      if (nat < 1)
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
      else if (nat==2)
	{
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%[^\n]\n", line); 
	      if (MCsim)
		{
		  sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
			 &r[0][i], &r[1][i], &r[2][i], &R[i][0][0], 
			 &R[i][0][1], &R[i][0][2], &R[i][1][0], 
			 &R[i][1][1], &R[i][1][2], &R[i][2][0], 
			 &R[i][2][1], &R[i][2][2], &partType); 

		}
	      else
		{
		  sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
			 &r[0][i], &r[1][i], &r[2][i], &R[i][0][0], 
			 &R[i][0][1], &R[i][0][2], &R[i][1][0], 
			 &R[i][1][1], &R[i][1][2], &R[i][2][0], 
			 &R[i][2][1], &R[i][2][2]); 
		}
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
  fprintf(stderr,"Usage: cnftool --partial/-p [--removepart||-rp] <val> [--reordconf||-rc <val> ] <file_da_convertire>\n");
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
      else if (!strcmp(argv[cc],"--partial")||!strcmp(argv[cc],"-p"))
	{
	  partial = 1;
	}
      else if (!strcmp(argv[cc],"--removepart")||!strcmp(argv[cc],"-rp"))
	{
	  cc++;
	  if (cc >= argc)
	    {
  	      print_usage();
	      exit(-1);
	    }
	  rempar=atoi(argv[cc]);
	}

      else if (!strcmp(argv[cc],"--reordconf")||!strcmp(argv[cc],"-rc"))
	{
	  cc++;
	  if (cc >= argc)
	    {
  	      print_usage();
	      exit(-1);
	    }
	  reorder=atoi(argv[cc]);
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
void full_reorder_conf(int reord)
{
  int ifree, cc, i, k1, k2, j;
  double L2, r0tmp[3], Rtmp[3][3];
  /* reord=0->x-axis 1->y-axis 2->z-axis */
#if 0
  switch (reord)
    {
    case 0:
      L2 = Lx*0.5;
      break;
    case 1:
      L2 = Ly*0.5;
      break;
    case 2:
      L2 = Lz*0.5;
      break;
    }
#endif
  cc=0;
  for (i=0; i < NP-1; i++)
    {
      for (j=i+1; j > 0; j--)
	{
	    
	  if (r0[reord][j] < r0[reord][j-1]) 
	    {
	      /* fai in modo che le particelle con i=0...N/2 siano nella metà (-L2,0] della scatola */ 
	      for (k1=0; k1 < 3; k1++)
		{
		  r0tmp[k1] = r0[k1][j];
		  for (k2=0; k2 < 3; k2++)
		    Rtmp[k1][k2] = R[j][k1][k2];
		}
	      for (k1=0; k1 < 3; k1++)
		{
		  r0[k1][j] = r0[k1][j-1];
		  for (k2=0; k2 < 3; k2++)
		    R[j][k1][k2] = R[j-1][k1][k2];
		}
	      for (k1=0; k1 < 3; k1++)
		{
		  for (k2=0; k2 < 3; k2++)
		    R[j-1][k1][k2] = Rtmp[k1][k2];
		  r0[k1][j-1] = r0tmp[k1];
		}

	    }
	  else 
	    break;
	}
    }
}
void reorder_conf(int reord)
{
  int ifree, cc, i, k1, k2;
  double L2, r0tmp[3], Rtmp[3][3];
  /* reord=0->x-axis 1->y-axis 2->z-axis */
#if 0
  switch (reord)
    {
    case 0:
      L2 = Lx*0.5;
      break;
    case 1:
      L2 = Ly*0.5;
      break;
    case 2:
      L2 = Lz*0.5;
      break;
    }
#endif
  cc=0;
  for (i=0; i < NP; i++)
    {
      if (r0[reord][i] < 0.0)
	{
	  ifree = 0;
	  while (r0[reord][ifree] < 0.0)
	    {
	      ifree++;
	      if (ifree >= i)
		break;
	    }
	  if (ifree >= i)
	    continue;
	  /* fai in modo che le particelle con i=0...N/2 siano nella metà (-L2,0] della scatola */ 
	  for (k1=0; k1 < 3; k1++)
	    {
	      r0tmp[k1] = r0[k1][ifree];
	      for (k2=0; k2 < 3; k2++)
		Rtmp[k1][k2] = R[ifree][k1][k2];
	    }
	  for (k1=0; k1 < 3; k1++)
	    {
	      r0[k1][ifree] = r0[k1][i];
	      for (k2=0; k2 < 3; k2++)
		R[ifree][k1][k2] = R[i][k1][k2];
		    }
	  for (k1=0; k1 < 3; k1++)
	    {
	      for (k2=0; k2 < 3; k2++)
		R[i][k1][k2] = Rtmp[k1][k2];
	      r0[k1][i] = r0tmp[k1];
	    }
	  cc++;
	}
	if (cc==NP/2)
	  break;
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
  f = fopen(inputfile, "r");
  c2 = 0;
  maxl = 0;
  Rt = malloc(sizeof(double*)*3);
  for (kk=0; kk < 3; kk++)
    Rt[kk] = malloc(sizeof(double)*3);
#if 0
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
#endif
  nat = 0;
  //while (!feof(f) && nat < 4) 
  //{
  if (MCsim && nat==2)
    {
      /* <parnum> 
	 2 1 1
	 16 2 2
	 1 1 1 1 2 0
	 2 0
	 1.541665 0 0 1.21667
	 -1.541665 0 0 1.21667
	 0 0 0 0 1 0 0 100000
	 0 0 0 1 1 0 0 100000
	 0 1 0 1 1 0 0 100000  */
      fscanf(f, "%d\n", &dummyint);
      fscanf(f, "%lf %lf %lf\n", &sax, &say, &saz);
      for (i=0; i < 8; i++)
	fscanf(f, "%[^\n]\n", line);
      X0 = sax/say;
      //printf("X0=%f\n", X0);
      fscanf(f, "%[^\n]\n)", line);
    }
  else
    fscanf(f, "%[^\n]\n)", line);

  if (!strcmp(line,"@@@"))
    {
      nat++;
      if (nat==3)
	{
	  //printf("NP=%d\n", NP);
	  if (MCsim)
	    {
	      for (i=0; i < NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      fscanf(f, "%lf %lf %lf\n", &Lx, &Ly, &Lz);
	    }
	  else
	    {
	      for (i=0; i < 2*NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      fscanf(f, "%lf\n", &L);
	    }
	}
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
  if (!strcmp(parname, "ensembleMC"))
    {
      MCsim=1;
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
  //}
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
  
  if (eventDriven==0 && MCsim==0)
    L = cbrt(L);
  if (NPA != NP)
    fprintf(stderr,"[MIXTURE] points=%d files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", points, nfiles, NP, NPA, L, NN, maxl);
  else
    fprintf(stderr,"[MONODISPERSE] points=%d files=%d NP = %d L=%.15G NN=%d maxl=%d\n", points, nfiles, NP, L, NN, maxl);
  if (MCsim==1)
    fprintf(stderr,"[MC] Monte Carlo Simulation X0=%f reorder=%d\n", X0, reorder);
  else if (eventDriven)
    fprintf(stderr,"[ED] Event-Driven simulation\n");
  else
    fprintf(stderr,"[MD] Time-Driven simulation\n");
  DR = malloc(sizeof(double*)*NP);
  for (ii = 0; ii < NP; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*3);
    }

  fprintf(stderr, "reord. axis=%d, num. of particle to remove: %d\n", reorder, rempar); 
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

  readconf(inputfile, &time, &refTime, NP, r0, w0, DR, R);
  if (reorder!=-1)
    {
      if (partial)
	reorder_conf(reorder);
      else
	full_reorder_conf(reorder);
    }
  //save_fra(fname[nr1]);
  outconf(inputfile);
  return 0;
}
