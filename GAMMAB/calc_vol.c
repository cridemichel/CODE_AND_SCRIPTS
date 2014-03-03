#include<stdio.h>
#include<stdlib.h>

char inputfile[4096];
int numat;
double *pos[3];
double *rad, proberad, cells[3];
double com[3], maxdist[3], L[3], maxrad;
double *cellList, *inCell[3];
long long int outits=100000;
void rebuildLinkedList(void)
{
  int j, n;

  for (j = 0; j < cells[0]*cells[1]*cells[2] + numat; j++)
    cellList[j] = -1;
  /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
 
  for (n = 0; n < numat; n++)
    {
      inCell[0][n] =  (pos[n] + L2[0]) * cells[0] / L[0];
      inCell[1][n] =  (pos[n] + L2[1]) * cells[1] / L[1];
      inCell[2][n] =  (pos[n] + L2[2]) * cells[2] / L[2];
      j = (inCell[2][n]*cells[1] + inCell[1][n])*cells[0] + 
	inCell[0][n] + numat;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}

void print_usage(void)
{
  printf("calc_vol [--outits|-oi <output iterations>] [--max-trials|-mt <number of trials>]  [--radius|-r <probe_radius (default=1.5)] <coord_file>\n");
  exit(0);
}

void parse_param(int argc, char** argv)
{
  int cc=1;
  int extraparam=0;
  if (argc==1)
    {
      print_usage();
      exit(1);
    }
  
  while (cc < argc)
    {
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
        {
          print_usage();
        }
      else if (!strcmp(argv[cc],"--radius") || !strcmp(argv[cc],"-r" ))
        {
          cc++;
	  if (cc == argc)
	    print_usage();
	  proberad = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--max-trials") || !strcmp(argv[cc],"-mt" ))
        {
          cc++;
	  if (cc == argc)
	    print_usage();
	  tt = atoi(argv[cc]);
        }
      else if (!strcmp(argv[cc],"--outits") || !strcmp(argv[cc],"-oi" ))
     	{
     	  cc++;
	  if (cc == argc)
	    print_usage();
	  outits = atoi(argv[cc]);
        }
      else if (extraparam == 0)
	{ 
	  extraparam = 1;
	  strcpy(inputfile,argv[cc]);
	}
      else
        print_usage();
   cc++;
  };
#if 0  
  strcpy(inputfile,argv[cc]);
  if (argc == 3)
    points=atoi(argv[2]);
  else
    points=100;
#endif
}
void read_file(void)
{ 
  FILE* f;
  int cc=0;
  f = fopen(inputfile, "r");
  /* count atoms */
  while (!feof(f))
    {
      cc++;
    }
  numat=cc+1;
  rewind(f);
  for (kk=0; kk < 3; kk++)
    {
      /* l'ultima particella è il probe */
      pos[kk] = malloc(sizeof(double)*numat);
      rad = malloc(sizeof(double)*numat);
    }
  for (cc=0; cc < numat; cc++)
    {
      fscanf(f, "%lf %lf %lf %lf\n", &(pos[0][cc]),&(pos[1][cc]),&(pos[2][cc]), &(rad[cc]));
    }
  fclose(f);
}

int main(int argc, char **argv)
{
  long long int tt;
  double totov, vol;
  int MAXTRIALS, kk, i, j;

  read_file();

  for (kk=0; kk < 3; kk++)
    com[kk] = 0.0;
  for (i = 0; i < numat; i++)
    {
      for (kk=0; kk < 3; kk++)
	com[kk] + = pos[kk];
    }
  for (kk=0; kk < 3; kk++)
    com[kk] /= ((double)numat);

  maxrad = 0.0;
  for (i = 0; i < numat; i++)
    {
      for (kk=0; kk < 3; kk++)
	pos[kk][i] -= com[kk];
      if (maxrad < rad[i])
       maxrad = rad[i];	
    }
  for (kk=0; kk < 3; kk++)
    maxdist[kk] = 0.0;
  for (i = 0; i < numat; i++)
    {
      for (j = i+1; j < numat; j++)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      if (maxdist[kk] < fabs(pos[kk][i]-pos[kk][j]))
		maxdist[kk] = fabs(pos[kk][i]-pos[kk][j]);
	    }
	}
    }
  for (kk=0; kk < 3; kk++)
    {
      L[kk] = maxdist[kk] + 2.1*proberad;
      cells[kk] = L[kk]/maxrad;
    }
  cellList = malloc(sizeof(int)*(cells[0]*cells[1]*cells[2]+numat));
  inCell[0] = malloc(sizeof(int)*numat);
  inCell[1] = malloc(sizeof(int)*numat);
  inCell[2] = malloc(sizeof(int)*numar);
  f = fopen("volume.dat","w+");
  fclose(f); 
  totov=0.0;
  for (tt=0; tt < MAXTRIALS; tt++)
    {
      for (kk = 0; kk < 3; kk++)
	pos[kk][numat-1] = L[kk]*(drand48()-0.5); 
      build_linked_lists();
      if (check_overlap())
	totov++;
      if (tt % outits==0)
	{
	   if (tt!=0)
	    {
	      vol = (totov/((double)tt))*(L[0]*L[1]*L[2]);
	      f=fopen("volume.dat", "a");
	      printf("volume=%.10f (totov=%f/%lld)\n", vol, totov, tt);
	      fprintf(f, "%lld %.15G %.15G\n", tt, vol, totene);
	      fclose(f);
	      sync();
	    }
	}
    }
  printf("DONE volume=%.12f\n",(totov/((double)tt))*(L[0]*L[1]*L[2]));
}
