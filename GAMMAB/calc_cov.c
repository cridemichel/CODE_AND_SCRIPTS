#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define Sqr(x) ((x)*(x))

char inputfile[4096];
int numat, cells[3];
double *pos[3], *posBody[3];
double *rad, proberad;
double com[3], maxdist[3], L[3], L2[3], maxrad;
double R[3][3];
int *cellList;
int *inCell[3];
long long int outits=100000, maxtrials=1000000;

/* apply a random rotation around the supplied axis because 
   bent cylinders do not have azimuthal symmetry */
double thetaGlobalBondangle;
void orient(double *omx, double *omy, double* omz)
{
  int i;
  //double inert;                 /* momentum of inertia of the molecule */
  //double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, osq, norm;
  
  //Mtot = m; /* total mass of molecule */

  //inert = I; /* momentum of inertia */
 
  //mean = 3.0*temp / inert;

  xisq = 1.0;

  while (xisq >= 1.0)
    {
      xi1  = drand48() * 2.0 - 1.0;
      xi2  = drand48() * 2.0 - 1.0;
      xisq = xi1 * xi1 + xi2 * xi2;
    }

  xi = sqrt (fabs(1.0 - xisq));
  ox = 2.0 * xi1 * xi;
  oy = 2.0 * xi2 * xi;
  oz = 1.0 - 2.0 * xisq;

  /* Renormalize */
  osq   = ox * ox + oy * oy + oz * oz;
  norm  = sqrt(fabs(osq));
  ox    = ox / norm;
  oy    = oy / norm;
  oz    = oz / norm;
  *omx = ox;
  *omy = oy;
  *omz = oz; 
}

void add_rotation_around_axis(double ox, double oy, double oz, double Rin[3][3], double Rout[3][3])
{
  double theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];
  int k1, k2, k3;
  /* pick a random rotation angle between 0 and 2*pi*/
  theta = 4.0*acos(0.0)*drand48();
  /* set to be used in az. angle distro calculation */
  thetaGlobalBondangle = theta;

  thetaSq=Sqr(theta);
  sinw = sin(theta);
  cosw = (1.0 - cos(theta));
  Omega[0][0] = 0;
  Omega[0][1] = -oz;
  Omega[0][2] = oy;
  Omega[1][0] = oz;
  Omega[1][1] = 0;
  Omega[1][2] = -ox;
  Omega[2][0] = -oy;
  Omega[2][1] = ox;
  Omega[2][2] = 0;
  OmegaSq[0][0] = -Sqr(oy) - Sqr(oz);
  OmegaSq[0][1] = ox*oy;
  OmegaSq[0][2] = ox*oz;
  OmegaSq[1][0] = ox*oy;
  OmegaSq[1][1] = -Sqr(ox) - Sqr(oz);
  OmegaSq[1][2] = oy*oz;
  OmegaSq[2][0] = ox*oz;
  OmegaSq[2][1] = oy*oz;
  OmegaSq[2][2] = -Sqr(ox) - Sqr(oy);

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = Rin[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += Rin[k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     Rout[k1][k2] = Ro[k1][k2]; 
}

void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
  double Rout[3][3];
  int k1, k2;
  /* first row vector */
  R[0][0] = ox;
  R[0][1] = oy;
  R[0][2] = oz;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[0][0] && u[1]==R[0][1] && u[2]==R[0][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[0][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[0][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
  if (typesArr[0].nspots==3 && type==0)
    {
      for (k=0; k < 3 ; k++)
	u[k] = R[1][k];
      vectProdVec(R[0], u, up);
      /* rotate randomly second axis */
      angle=4.0*acos(0.0)*ranf_vb();
      xx = cos(angle);
      yy = sin(angle);
      for (k=0; k < 3 ; k++)
	R[1][k] = u[k]*xx + up[k]*yy;
      //printf("calc_norm(R[1])=%.15G\n", calc_norm(R[1]));
    }
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
  /* add a random rotation around the axis (ox, oy, oz) */
  add_rotation_around_axis(ox, oy, oz, R, Rout);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R[k1][k2] = Rout[k1][k2];
}

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
      // printf("n=%d inCell=%d %d %d\n", n, inCell[0][n], inCell[1][n],inCell[2][n]);
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
  double a, b, c, d;
  int cc=0, kk;
  f = fopen(inputfile, "r");
  /* count atoms */
  while (!feof(f))
    {
      fscanf(f, "%lf %lf %lf %lf ", &a, &b, &c, &d);
      cc++;
    }
  numat=2*cc;
  rewind(f);
  for (kk=0; kk < 3; kk++)
    {
      /* l'ultima particella è il probe */
      pos[kk] = malloc(sizeof(double)*numat);
      posBody[kk] = malloc(sizeof(double)*numat);
      rad = malloc(sizeof(double)*numat);
    }
  for (cc=0; cc < numat-1; cc++)
    {
      fscanf(f, "%lf %lf %lf %lf\n", &(pos[0][cc]),&(pos[1][cc]),&(pos[2][cc]), &(rad[cc]));
    }
  fclose(f);
}

int check_overlap_ij(int i, int j)
{
  double distSq, sigmaSq, avsigma;
  int kk;

  avsigma = rad[i]+rad[j];
  //printf("rad[%d]=%f rad[%d]=%f\n", i, rad[i], j, rad[j]);
  sigmaSq = Sqr(avsigma);
  distSq = 0.0;
  for (kk=0; kk < 3; kk++)
    distSq += Sqr(pos[kk][i] - pos[kk][j]);
  //printf("dist=%f avrad=%f\n", sqrt(distSq), avrad);
  if ( distSq < sigmaSq )
    return -1;  
  else
    return 1;
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
		      /* gli atomi devono appartenere a proteine diverse */
		      if ( !((n < numat && na >= numat) || (n >= numat && na < numat)) )
			continue;
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
int check_overlap_12(void)
{
  int i;
  for (i=0; i < numat/2; i++)
    {
      if(check_overlap(i))
	return 1;
    }
  return 0;
}
int main(int argc, char **argv)
{
  FILE *f;
  long long int tt;
  double totov, vol;
  int kk, i, j, k1, k2;

  parse_param(argc, argv);

  read_file();
  for (kk=0; kk < 3; kk++)
    com[kk] = 0.0;
  /*  la prima proteina sarà fissa in (0,0,0) con il suo centro di massa */
  for (i = 0; i < numat/2; i++)
    {
      for (kk=0; kk < 3; kk++)
	com[kk] += pos[kk][i];
    }
  for (kk=0; kk < 3; kk++)
    com[kk] /= ((double)numat)/2.0;

  maxrad = 0.0;
  for (i = 0; i < numat; i++)
    {
      for (kk=0; kk < 3; kk++)
	pos[kk][i] -= com[kk];
      if (maxrad < rad[i])
       maxrad = rad[i];	
    }

  for (i = numat/2; i < numat; i++)
    {
      for (kk=0; kk < 3; kk++)
	com[kk] += pos[kk][i];
    }
  for (kk=0; kk < 3; kk++)
    com[kk] /= ((double)numat)/2.0;

  for (i = numat/2; i < numat; i++)
    {
      for (kk=0; kk < 3; kk++)
	pos[kk][i] -= com[kk];
    }

  for (kk=0; kk < 3; kk++)
    maxdist[kk] = 0.0;

  for (i = 0; i < numat/2; i++)
    {
      for (j = i+1; j < numat/2; j++)
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
      L[kk] = 3.1*maxdist[kk];
      L2[kk] = 0.5*L[kk];
      cells[kk] = L[kk]/(2.0*maxrad);
    }
  /* set probe radius */
  rad[numat-1] = proberad;
  cellList = malloc(sizeof(int)*(cells[0]*cells[1]*cells[2]+numat));
  inCell[0] = malloc(sizeof(int)*numat);
  inCell[1] = malloc(sizeof(int)*numat);
  inCell[2] = malloc(sizeof(int)*numat);
  printf("cells = %d %d %d  atoms=%d \n", cells[0], cells[1], cells[2], numat);
  f = fopen("covolume.dat","w+");
  fclose(f); 
  totov=0.0;
  for (i=numat/2; i < numat; i++)
    {
      for (kk = 0; kk < 3; kk++)
	posBody[kk][i] = pos[kk][i];
    }      
  for (tt=0; tt < maxtrials; tt++)
    {
      orient(&ox, &oy, &oz); 
      versor_to_R(ox, oy, oz, Rl);
      for (kk = 0; kk < 3; kk++)
	com[kk] = L[kk]*(drand48()-0.5); 
      for (i=numat/2; i < numat; i++)
	{
	  for (k1 = 0; k1 < 3; k1++)
	    {
	      pos[k1][i] = com[k1];
	      for (k2 = 0; k2 < 3; k2++)
		pos[k1][i] += R[k2][k1]*posBody[k2][i];
	    }
	}	

      build_linked_lists();

	
      if (check_overlap_12())
	totov++;
      if (tt % outits==0)
	{
	   if (tt!=0)
	    {
	      vol = (totov/((double)tt))*(L[0]*L[1]*L[2]);
	      f=fopen("covolume.dat", "a");
	      printf("covolume=%.10f (totov=%f/%lld)\n", vol, totov, tt);
	      fprintf(f, "%lld %.15G %.15G\n", tt, vol, totov);
	      fclose(f);
	      sync();
	    }
	}
    }
  printf("DONE covolume=%.12f\n",(totov/((double)tt))*(L[0]*L[1]*L[2]));
}
