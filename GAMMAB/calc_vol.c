#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define Sqr(x) ((x)*(x))

char inputfile[4096];
int numat;
double *pos[3];
double *rad, proberad, cells[3];
double com[3], maxdist[3], L[3], L2[3], maxrad;
double *cellList, *inCell[3];
long long int outits=100000, maxtrials=1000000;
void build_linked_lists(void)
{
  int j, n;

  for (j = 0; j < cells[0]*cells[1]*cells[2] + numat; j++)
    cellList[j] = -1;
  /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
 
  for (n = 0; n < numat; n++)
    {
      inCell[0][n] =  (pos[0][n] + L2[0]) * cells[0] / L[0];
      inCell[1][n] =  (pos[1][n] + L2[1]) * cells[1] / L[1];
      inCell[2][n] =  (pos[2][n] + L2[2]) * cells[2] / L[2];
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
	  maxtrials = atoi(argv[cc]);
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
  int cc=0, kk;
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

int check_overlap_ij(int i, int j)
{
  double distSq, radSq, avrad;
  int kk;

  avrad = 0.5*(rad[i]+rad[j]);

  radSq = Sqr(avrad);
  distSq = 0.0;
  for (kk=0; kk < 3; kk++)
    distSq += Sqr(pos[kk][i] - pos[kk][j]);
  if ( distSq < radSq)
    return 1;  
  else
    return 0;
}

int check_overlap(int ip)
{
  int kk, nb, k, iZ, jZ, iX, jX, iY, jY, n, na;
  int cellRangeT[6];
  int cellRange[6];
  double shift[3];
  na=ip;
  nb=-1;
  for (kk = 0;  kk < 3; kk++)
    {
      cellRange[2*kk]   = - 1;
      cellRange[2*kk+1] =   1;
    }

  for (k = 0; k < 6; k++) cellRangeT[k] = cellRange[k];

  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
      if (jZ == -1) 
	{
	  jZ = cells[2] - 1;    
	  shift[2] = - L[2];
	} 
      else if (jZ == cells[2]) 
	{
	  jZ = 0;    
	  shift[2] = L[2];
	}
      
      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	{
	  jY = inCell[1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cells[1] - 1;    
	      shift[1] = -L[1];
	    } 
	  else if (jY == cells[1]) 
	    {
	      jY = 0;    
	      shift[1] = L[1];
	    }
	  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
	    {
	      jX = inCell[0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cells[0] - 1;    
		  shift[0] = - L[0];
		} 
	      else if (jX == cells[0]) 
		{
		  jX = 0;   
		  shift[0] = L[0];
		}
	      n = (jZ *cells[1] + jY) * cells[0] + jX + numat;
	      for (n = cellList[n]; n > -1; n = cellList[n]) 
		{
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
		      if (check_overlap_ij(na, n)<0.0)
			{
			  return 1;
			}
		    }
		} 
	    }
	}
    }
  return 0;
}
int main(int argc, char **argv)
{
  FILE *f;
  long long int tt;
  double totov, vol;
  int kk, i, j;

  parse_param(argc, argv);

  read_file();

  for (kk=0; kk < 3; kk++)
    com[kk] = 0.0;
  for (i = 0; i < numat; i++)
    {
      for (kk=0; kk < 3; kk++)
	com[kk] += pos[kk][i];
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
      L2[kk] = 0.2*L[kk];
      cells[kk] = L[kk]/maxrad;
    }
  /* set probe radius */
  rad[numat-1] = proberad;
  cellList = malloc(sizeof(int)*(cells[0]*cells[1]*cells[2]+numat));
  inCell[0] = malloc(sizeof(int)*numat);
  inCell[1] = malloc(sizeof(int)*numat);
  inCell[2] = malloc(sizeof(int)*numat);
  f = fopen("volume.dat","w+");
  fclose(f); 
  totov=0.0;
  for (tt=0; tt < maxtrials; tt++)
    {
      for (kk = 0; kk < 3; kk++)
	pos[kk][numat-1] = L[kk]*(drand48()-0.5); 
      build_linked_lists();
      if (check_overlap(numat-1))
	totov++;
      if (tt % outits==0)
	{
	   if (tt!=0)
	    {
	      vol = (totov/((double)tt))*(L[0]*L[1]*L[2]);
	      f=fopen("volume.dat", "a");
	      printf("volume=%.10f (totov=%f/%lld)\n", vol, totov, tt);
	      fprintf(f, "%lld %.15G %.15G\n", tt, vol, totov);
	      fclose(f);
	      sync();
	    }
	}
    }
  printf("DONE volume=%.12f\n",(totov/((double)tt))*(L[0]*L[1]*L[2]));
}
