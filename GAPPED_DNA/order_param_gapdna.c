#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <lapack.h>

extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );


double max3(double a, double b, double c)
{
  double m;
  m = a;
  if (b > m)
    m = b;
  if (c > m)
    m = c;
  return m;
}
char line[100000], parname[124], parval[256000];
int N, mcsim=0;
double x[3], R[3][3], Q[3][3], *Qn[3][3], *Nn;
double r1[3], r2[3], r3[3], u1[3], u2[3], u3[3], dt, L[3];
double *ev_n[3];
double ppos;
int plane=-1; /* -1=none 0 = x 1=y 2=z */
char fname[1024], inputfile[1024];
int readCnf = 0, timeEvol = 0, ordmatrix=0, eigenvectors=0, curStep, readLapo=0, heflex=0;
int nslabs=2;
#define Sqr(x) ((x)*(x))
void vectProd(double r1x, double r1y, double r1z, 
	 double r2x, double r2y, double r2z, 
	 double* r3x, double* r3y, double* r3z)
{
  /* DESCRIPTIOM:
     r3 = [ r1, r2 ] where [ , ] the vectorial product */
  *r3x = r1y * r2z - r1z * r2y; 
  *r3y = r1z * r2x - r1x * r2z;
  *r3z = r1x * r2y - r1y * r2x;
}
void get_L(FILE* f, double L[3])
{
  int i;
  if (readCnf)
    {
      do 
	{
	  fscanf(f,"%[^\n]\n",line);
	}
      while (strcmp(line,"@@@"));
    }
  do 
    {
      fscanf(f,"%[^\n]\n",line);
    }
  while (strcmp(line,"@@@"));
  /* BUG FIX 21/01/13 */
  if (heflex)
    {
      do 
	{
	  fscanf(f,"%[^\n]\n",line);
	}
      while (strcmp(line,"@@@"));
    }

  for (i=0; i < N; i++)
    fscanf(f,"%[^\n]\n",line);
  if (!mcsim)
    {
      for (i=0; i < N; i++)
	fscanf(f,"%[^\n]\n",line);
    }
  if (fscanf(f, "%lf %lf %lf\n", &L[0], &L[1], &L[2]) == 1)
    {
      L[2] = L[0];
      L[1] = L[0];
    }
  printf("L=%f %f %f\n", L[0], L[1], L[2]); 
}
double eigvec[3][3], *eigvec_n[3][3], eigvec_t[3][3];
void diagonalize(double M[3][3], double ev[3])
{
  double a[9], work[45];
  char jobz, uplo;
  int info, i, j, lda, lwork;
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) a[j+3*i]=M[j][i];		
      //for(j=0; j<3; j++) a[j][i]=M[j][i];		
    }	
  lda = 3;
  if (eigenvectors)
    jobz='V';
  else
    jobz='N';
  uplo='U';
  lwork = 45;
  dsyev_(&jobz, &uplo, &lda, a, &lda, ev, work, &lwork,  &info);  
  if (!eigenvectors)
    return;
  for (i=0; i<3; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
      for(j=0; j<3; j++) eigvec[i][j]=a[j+3*i];		
    }	
}
void print_usage(void)
{
  printf("order_param [ --nslabs/-n | --plane/-p | --ppos/-x | -mc | --lapo/-l | --heflex/-fl |--cnf/-c | --time/-t | --ordmatrix/-Q ] [--eigenvec/-ev] <confs_file>\n");
  exit(0);
}
void parse_param(int argc, char** argv)
{
  int cc=1;
  
  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"--nslabs") || !strcmp(argv[cc],"-n"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  nslabs = atoi(argv[cc]);
	}

      else if (!strcmp(argv[cc],"--plane") || !strcmp(argv[cc],"-p"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  plane = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--ppos") || !strcmp(argv[cc],"-x"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  ppos = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"-mc"))
	{
	  mcsim = 1;
	} 
      else if (!strcmp(argv[cc],"--cnf") || !strcmp(argv[cc],"-c" ))
	{
	  readCnf = 1;
	} 
      else if (!strcmp(argv[cc],"--heflex") || !strcmp(argv[cc],"-fl" ))
	{
	  heflex = 1;
	} 
      else if (!strcmp(argv[cc],"--lapo") || !strcmp(argv[cc],"-l" ))
	{
	  readLapo = 1;
	} 
      else if (!strcmp(argv[cc],"--time") || !strcmp(argv[cc],"-t"))
	{
	  timeEvol = 1;
	}
      else if (!strcmp(argv[cc],"--ordmatrix") || !strcmp(argv[cc],"-Q"))
	{
	  ordmatrix = 1;
	}
      else if (!strcmp(argv[cc],"--eigenvec") || !strcmp(argv[cc],"-ev"))
	{
	  eigenvectors = 1;
	}
      else if (cc == argc)
	print_usage();
      else
	strcpy(inputfile,argv[cc]);
      cc++;
    }
}
double calc_norm(double *v)
{
  int a;
  double norm=0;
  for (a = 0; a < 3; a++)
   norm += Sqr(v[a]);
  return sqrt(norm); 
}
void build_ref_axes2(double u1[3], double u2[3], double u3[3])
{
  int a;
  double norm1, u1u2, norm2;
  u1[0] = r2[0] - r1[0];
  u1[1] = r2[1] - r1[1];
  u1[2] = r2[2] - r1[2];
  u2[0] = r3[0] - r1[0];
  u2[1] = r3[1] - r1[1];
  u2[2] = r3[2] - r1[2];
  norm2 = calc_norm(u2);
  for (a=0;  a < 3; a++)
    u2[a] /= norm2;
  u1u2 = u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
  for (a=0;  a < 3; a++)
    {
      u1[a] = u1[a] - u1u2*u2[a];
    }
  norm1 = calc_norm(u1);
  for (a=0;  a < 3; a++)
    {
      u1[a] /= norm1;
    }
  vectProd(u1[0],u1[1],u1[2],u2[0],u2[1], u2[2],&u3[0], &u3[1], &u3[2]); 
}

void build_ref_axes(double u1[3], double u2[3], double u3[3])
{
  int a;
  double n;
  u1[0] = r2[0] - r1[0];
  u1[1] = r2[1] - r1[1];
  u1[2] = r2[2] - r1[2];
  n = calc_norm(u1);
  for (a = 0; a < 3; a++)
    u1[a] /= n;
    
  u2[0] = 2.0*r3[0]-(r1[0] + r2[0]);
  u2[1] = 2.0*r3[1]-(r1[1] + r2[1]);
  u2[2] = 2.0*r3[2]-(r1[2] + r2[2]);
  
  n = calc_norm(u2);
  for (a = 0; a < 3; a++)
    u2[a] /= n;
   
  //fprintf(stderr,"u1.u2=%.15G\n", u1[0]*u2[0]+u1[1]*u2[1]+u1[2]*u2[2]);
  vectProd(u1[0], u1[1], u1[2], u2[0], u2[1], u2[2], &u3[0], &u3[1], &u3[2]);  

}

void calc_I2(double theta, double phi, double *reI2, double *imI2)
{
  double A;
  double sin2th, cos2th, cosphi, sinphi, cos2phi, sin2phi, costh, sinth;
  costh= cos(theta);
  sinth= sin(theta);
  cos2th = Sqr(cos(theta));
  //cos4th = Sqr(cos2th);
  sin2th = Sqr(sin(theta));
  //sin4th = Sqr(sin2th);
  cosphi = cos(phi);
  sinphi = sin(phi);
  cos2phi = cos(2.0*phi);
  sin2phi = sin(2.0*phi);
  A = (1.0/4.0)*sqrt(15.0/2.0/M_PI)*sin2th;
  /* 0, 1... are m=-2, -1...*/
  reI2[0] =  cos2phi*A;
  imI2[0] =  -sin2phi*A;
  reI2[4] =  cos2phi*A;
  imI2[4] =  sin2phi*A;
  A = (1.0/2.0)*sqrt(15.0/2.0/M_PI)*costh*sinth;
  reI2[1] = cosphi*A;
  imI2[1] = -sinphi*A;
  reI2[3] = cosphi*A;
  imI2[3] = sinphi*A;
  A = (1.0/4.0)*sqrt(5.0/M_PI)*(-1.0+3.0*cos2th);
  reI2[2] = A;
  imI2[2] = 0;
}

void calc_I4(double theta, double phi, double *reI2, double *imI2)
{
  double A;
  double sin2th, sin4th, cos2th, cos4th, cosphi, sinphi, cos2phi, sin2phi, cos3phi, sin3phi;
  double costh, sinth, cos4phi, sin4phi;
  costh = cos(theta);
  sinth = sin(theta);
  cos2th = Sqr(cos(theta));
  cos4th = Sqr(cos2th);
  sin2th = Sqr(sin(theta));
  sin4th = Sqr(sin2th);
  cosphi = cos(phi);
  sinphi = sin(phi);
  cos2phi = cos(2.0*phi);
  sin2phi = sin(2.0*phi);
  cos3phi = cos(3.0*phi);
  sin3phi = sin(3.0*phi);
  cos4phi = cos(4.0*phi);
  sin4phi = sin(4.0*phi); 
  /* here 0,1,... are m=-4,3,... */
  A = (3.0/16.0)*sqrt(35.0/2.0/M_PI)*sin4th; 
  reI2[0] = cos4phi*A;
  imI2[0] = -sin4phi*A;
  reI2[8] = cos4phi*A;
  imI2[8] = sin4phi*A;

  A = (3.0/8.0)*sqrt(35.0/M_PI)*costh*sin2th*sinth; 
  reI2[1] = cos3phi*A;
  imI2[1] = -sin3phi*A;
  reI2[7] = cos3phi*A;
  imI2[7] = sin3phi*A;

  A = (3.0/8.0)*sqrt(5.0/2.0/M_PI)*(-1.0+7.0*cos2th)*sin2th;
  reI2[2] = cos2phi*A;
  imI2[2] = -sin2phi*A;
  reI2[6] = cos2phi*A;
  imI2[6] = sin2phi*A;

  A = (3.0/8.0)*sqrt(5.0/2.0/M_PI)*(-3.0+7.0*cos2th)*sinth; 
  reI2[3] = cosphi*A;
  imI2[3] = -sinphi*A;
  reI2[5] = cosphi*A;
  imI2[5] = sinphi*A;

  reI2[4] = 3.0*(3.0-30.0*cos2th+35*cos4th)/16.0/sqrt(M_PI);
  imI2[4] = 0.0;
}

int main(int argc, char** argv)
{
  FILE *f, *f2, *f_n;
  int nf, i, a, b, dummyint, k, nbin, mm;
  double sum_reI2[5], sum_imI2[5], sum_reI4[9], sum_imI4[9], I2, I4, theta, phi, reK2[5], imK2[5], reK4[9], imK4[9];
  double ev[3], ti=-1, S, tref=0.0, Dx, Qt[3][3], ev_t[3];
#if 0
  if (argc == 1)
    {
      printf("file with all confs to read as input\n");
      exit(-1);
    }
#endif
  parse_param(argc, argv);
  f2 = fopen(inputfile,"r");
  nf = 0;
  if (plane >= 0)
    {
      f_n = fopen("Sslabs.dat","w+");
    }
  if (plane >= 0)
    {
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  {
  	    Qn[a][b] = malloc(sizeof(double)*nslabs);
	    eigvec_n[a][b] = malloc(sizeof(double)*nslabs);
	  }
      Nn = malloc(sizeof(double)*nslabs);
      for (k=0; k < 3; k++)
	{
	  ev_n[k] = malloc(sizeof(double)*nslabs);
	}
    }
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Q[a][b] = 0.0;      
	if (plane >= 0)
	  {
	    for (k=0; k < nslabs; k++)
	      {
  		Qn[k][a][b] = 0.0;
		Nn[k] = 0.0;
	      }
	  }
      }
  I2=I4=0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", fname);
      //printf("fname=%s argv[2]=%s\n",fname, argv[2]);
      f = fopen(fname,"r");
      /* BUG FIX 21/01/13 */
      if (N!=-1 && plane>=0)
	{
	  get_L(f, L); 
	  rewind(f);
	}
      nf++;
      tref=0.0;
      if (readCnf)
	{
	  do 
	    {
	      fscanf(f,"%[^\n]\n",line);
	      sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	      if (!strcmp(parname,"refTime"))
		tref = atof(parval);
	    }
	  while (strcmp(line,"@@@"));
	}
      do 
	{
	  fscanf(f,"%[^\n]\n",line);
	  sscanf(line, "%[^:\n ]:%[^\n]\n", parname, parval);
	  if (!strcmp(parname,"parnum"))
	    N = atoi(parval);
	  if (!strcmp(parname,"time"))
	    ti = atof(parval);
	  if ((readLapo||mcsim) && !strcmp(parname,"curStep"))
	    curStep = atoi(parval);
	  if (readLapo && !strcmp(parname,"steplength"))
	    dt = atof(parval);
  
  	}
      while (strcmp(line,"@@@"));
      if (heflex)
	{
	  do 
	    {
	      fscanf(f,"%[^\n]\n",line);
	    }
	  while (strcmp(line,"@@@"));
	}
      //printf("fname=%s %d ellipsoids...\n", fname, N);
      if (timeEvol)
	{
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
	    	Q[a][b] = 0.0;      
		if (plane >= 0)
		  {
		    for (k=0; k < nslabs; k++)		    
		      Qn[k][a][b] = 0.0;
		    Nn[k] = 0;
		  }
	      }
	}

      for (mm=0; mm < 5; mm++)
	sum_reI2[mm]=sum_imI2[mm]=0;
      for (mm=0; mm < 9; mm++)
	sum_reI4[mm]=sum_imI4[mm]=0;
      for (i=0; i < N; i++)
	{
	  if (readLapo)
	    {
	      fscanf (f, "%lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n",
	  	    &(r1[0]), &(r1[1]), &(r1[2]), 
	  	    &(r2[0]), &(r2[1]), &(r2[2]), 
		    &(r3[0]), &(r3[1]), &(r3[2]));
	      build_ref_axes(u1, u2, u3);
	      for (a=0; a < 3; a++)
		{
		  R[0][a] = u1[a];
		  R[1][a] = u2[a];
		  R[2][a] = u3[a];
		}
	    }
	  else if (heflex)
	    fscanf (f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
	  	    &(x[0]), &(x[1]), &(x[2]), 
	  	    &(R[0][0]), &(R[0][1]), &(R[0][2]), &(R[1][0]), &(R[1][1]), 
	  	    &(R[1][2]), &(R[2][0]), &(R[2][1]), &(R[2][2]), &dummyint);
	  else
	    fscanf (f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
	  	    &(x[0]), &(x[1]), &(x[2]), 
	  	    &(R[0][0]), &(R[0][1]), &(R[0][2]), &(R[1][0]), &(R[1][1]), 
	  	    &(R[1][2]), &(R[2][0]), &(R[2][1]), &(R[2][2]));

	  theta = acos(R[0][2]);
	  phi = asin(R[0][1]/sin(theta));
	  calc_I2(theta, phi, reK2, imK2); 
	  for (mm=0; mm < 5; mm++)
	    {
	      sum_reI2[mm] += reK2[mm];
	      sum_imI2[mm] += imK2[mm];
	    }
	    
	  calc_I4(theta, phi, reK4, imK4);
	   for (mm=0; mm < 9; mm++)
	    {
	      sum_reI4[mm] += reK4[mm];
	      sum_imI4[mm] += imK4[mm];
	    }
	   //printf("R=%.15G %.15G %.15G\n", R[0][1], R[0][1], R[0][2]);
	  /* BUG FIX 21/01/13 */
	  if (plane>=0)
	    {
	      if (nslabs==2)
		{
		  if (x[plane] >= ppos)
		    Nn[1] += 1.0;
		  else
		    Nn[0] += 1.0;
		}
	      else
		{
		  Dx = L[plane]/nslabs; 
		  nbin = (int) ((x[plane]+L[plane]*0.5)/Dx);
		  Nn[nbin] += 1.0;
		}
	    }
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
		if (readLapo)
		  Q[a][b] += 1.5 * R[2][a]*R[2][b];
		else
		  Q[a][b] += 1.5 * R[0][a]*R[0][b];
		if (a==b)
		  Q[a][a] -= 0.5;
		if (plane >= 0)
		  {
		    if (nslabs==2)
		      {
			if (x[plane] >= ppos)
			  {
			    if (readLapo)
			      Qn[a][b][1] += 1.5 * R[2][a]*R[2][b];
			    else
			      Qn[a][b][1] += 1.5 * R[0][a]*R[0][b];
			    if (a==b)
			      Qn[a][a][1] -= 0.5;
			  }
			else
			  { 
			    if (readLapo)
			      Qn[a][b][0] += 1.5 * R[2][a]*R[2][b];
			    else
			      Qn[a][b][0] += 1.5 * R[0][a]*R[0][b];
			    if (a==b)
			      Qn[a][a][0] -= 0.5;
			  }
		      }
	    	    else
		      {
		      	if (readLapo)
	    		  Qn[nbin][a][b] += 1.5 * R[2][a]*R[2][b];
    			else
			  Qn[nbin][a][b] += 1.5 * R[0][a]*R[0][b];
			Qn[nbin][a][a] -= 0.5;
		      }
		  }
	      }
	}
      for (mm=0; mm < 5; mm++)
	{
  	  sum_reI2[mm] /= ((double)N);
	  sum_imI2[mm] /= ((double)N);
	}
      for (mm=0; mm < 9; mm++)
	{
   	  sum_reI4[mm] /= ((double)N);
       	  sum_imI4[mm] /= ((double)N);
	}
      if (timeEvol)
	{
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
	    	Q[a][b] /= ((double)N);
		if (plane >= 0)
		  {
		    for (k=0; k < nslabs; k++)
	      	      {
      			if (Nn[k] > 0)
			  Qn[a][b][k] /= ((double)Nn[k]);
		      }
		  }
	      } 
	  diagonalize(Q, ev);
	  if (plane >=0)
	    {
	      for (k=0; k < nslabs; k++)
		{
		  for (a=0; a < 3; a++)
		    {
		      for (b=0; b < 3; b++)
			Qt[a][b] = Qn[a][b][k];
		    }
		  diagonalize(Qt, ev_t);
		  for (a=0; a < 3; a++)
		    {
		      ev_n[a][k] = ev_t[a];
		      for (b=0; b < 3; b++)
			Qn[a][b][k] = Qt[a][b];
		    }
		}
	    }

	  if (mcsim) 
	    {
	      printf("%d %.15G %.15G %.15G\n", curStep, ev[0], ev[1], ev[2]); 
	      if (plane >= 0)
		{
		  fprintf(f_n, "%d ", curStep); 
		  for (k=0; k < nslabs; k++)
		    fprintf(f_n, "%.15G ", max3(ev_n[0][k], ev_n[1][k], ev_n[2][k]));
		  fprintf(f_n, "\n");
		}
	    }
	  else if (readLapo)
	    printf("%.15G %.15G %.15G %.15G\n", curStep*dt, ev[0], ev[1], ev[2]); 
	  else
	    printf("%.15G %.15G %.15G %.15G\n", ti+tref, ev[0], ev[1], ev[2]); 
	}

      fclose(f);
      for (mm = 0; mm < 5; mm++)
	{
	  I2 += (4.0*M_PI/5.0)*(Sqr(sum_reI2[mm]) + Sqr(sum_imI2[mm])) / ((double)nf);
	}
      for (mm = 0; mm < 9; mm++)
	{
 	  I4 += (4.0*M_PI/9.0)*(Sqr(sum_reI4[mm]) + Sqr(sum_imI4[mm])) / ((double)nf);
	}    
    }
  fclose(f2); 
  
  if (timeEvol)
    {
      if (plane >= 0)
	{	
	  fclose(f_n);
	}
      return 0;
    }
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      Q[a][b] /= ((double)N)*((double)nf); 
  //printf("N*nf=%f Nr=%f Nl=%f\n", ((double)N)*nf, Nr, Nl);
  I2 /= ((double)nf);
  I4 /= ((double)nf);
  printf("I2 = %f I4 = %f\n", I2, I4);
  if (plane >= 0)
    {
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  {
	    for (k=0; k < nslabs; k++)
	      {
	       	if (Nn[k] > 0)
		  Qn[a][b][k] /= Nn[k];
	      }
	  }
      for (k=0; k < nslabs; k++)
	{
	  for (a=0; a < 3; a++)
	    {
	      for (b=0; b < 3; b++)
		Qt[a][b] = Qn[a][b][k];
	    }
	  diagonalize(Qt, ev_t);
	  for (a=0; a < 3; a++)
	    {
	      ev_n[a][k] = ev_t[a];
	      for (b=0; b < 3; b++)
		eigvec_n[a][b][k] = eigvec[a][b];	  
	    }
	}
    }
  diagonalize(Q, ev);
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      eigvec_t[a][b] = eigvec[a][b];	  

  if (fabs(ev[0]) > fabs(ev[1]))
    S = ev[0];
  else
    S = ev[1];  
  if (fabs(ev[2]) > S)
    S = ev[2];
  printf("Order Parameter: %.8G (of %.8G %.8G %.8G)\n", S, 
	 ev[0], ev[1], ev[2]);
  if (plane >= 0)
    {
      for (k=0; k < nslabs; k++)
	{
	  if (fabs(ev_n[0][k]) > fabs(ev_n[1][k]))
	    S = ev_n[0][k];
	  else
	    S = ev_n[1][k];  
	  if (fabs(ev_n[2][k]) > S)
	    S = ev_n[2][k];

	  printf("[slab %d/%d] Order Parameter (plane = %f): %.8G (of %.8G %.8G %.8G)\n",
		 k, nslabs, ppos, S, 
		 ev_n[0][k], ev_n[1][k], ev_n[2][k]);
	}
    }
#if 1
  if (!timeEvol && ordmatrix)
    {
      printf("Ordering matrix Q:\n");
      printf("{");
      for (a=0; a < 3; a++)
	{
	  printf("{");
	  for (b=0; b < 3; b++)
	    {
	      if (b < 2)
		printf("%.15f,", Q[a][b]);	  
	      else
		printf("%.15f",Q[a][b]);
	    }
	  if (a < 2)
	    printf("},\n");
	  else if (b < 2)
	    printf("}\n");
	  else 
	    printf("}");
	}
      printf("}\n");
      if (plane >= 0)
	{
	  for (k=0; k < nslabs; k++)
	    {
	      printf("[slab %d/%d] Ordering matrix Qn:\n", k, nslabs);
	      printf("{");

	      for (a=0; a < 3; a++)
		{
		  printf("{");
		  for (b=0; b < 3; b++)
		    {
		      if (b < 2)
			printf("%.15f,", Qn[a][b][k]);	  
		      else
		    printf("%.15f",Qn[a][b][k]);
		    }
		  if (a < 2)
		    printf("},\n");
		  else if (b < 2)
		    printf("}\n");
		  else 
		    printf("}");
		}
	      printf("}\n");
	    }
	}
    }
  if (!timeEvol && eigenvectors)
    {
      printf("Eigenvectors matrix:\n");
      for (a=0; a < 3; a++)
	{
	  printf("lambda=%.15G {%.15f,%.15G,%.15G}\n",ev[a], eigvec_t[a][0], eigvec_t[a][1], eigvec_t[a][2]);	  
	  if (plane >= 0)
	    {
	      for (k=0; k < nslabs; k++)
		printf("[left] lambda=%.15G {%.15f,%.15G,%.15G}\n",
		       ev_n[a][k], eigvec_n[a][0][k], eigvec_n[a][1][k], eigvec_n[a][2][k]);	  
	    }
	}
    }
#endif
}
