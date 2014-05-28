#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#define Sqr(x) ((x)*(x))


struct vector 
{
  double x;
  double y;
  double z;
};

char infile[1024], e2efile[1024]="e2e.dat", bcparfile[1024]="bcpar.dat";
char string[1024], dummystr[256], a[1024];
struct vector PA1_1, PA2_1, PB2_1, PB1_1;
struct vector bar1, bar2, bar1_1, bar2_1;
struct vector *P;
double Lx=-1., Ly=-1., Lz=-1., ltotmin=38.0, ltotmax=42.0, l1l2min=0.8, l1l2max=1.2, frame, dist_ist, dist_aver;
int fixbroken=0, outframes=1000, nl=20, nlt=20, nphi=20; 
double Lhc, Dhc=20.0, Lhc1, Lhc2;
double PI, Prad=1.8;
int mglout=0, firstframe=-1, lastframe=-1; /* lastframe/firstfram=-1 means do not check */
FILE *mglfile;
char mglfn[1024]="out.mgl"; 
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


double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}

double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
void body2lab(double xp[3], double x[3], double rO[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
      x[k1] += rO[k1];
    }
}
void body2labR(double xp[3], double x[3], double R[3][3])
{
  int k1, k2;
  for (k1=0; k1 < 3; k1++)
    {
      x[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  x[k1] += R[k2][k1]*xp[k2];
       	} 
    }
}
#define DISTMAX 1E100
double dist_func(int i0, double l1, double l2, double phi, double Ro[3][3], double b1[3], double b2[3], 
		 double pos1[3], double n1[3], double pos2[3], double n2[3])
{
  /* params = {L_1, L_2, theta} */
  int k, i, dontcheck, kk;
  double fact, disttot, delx, angle, pos1B[3], pos2B[3];
  double vp[3], rp[3], drp[3], b12[3], distbb;
  double dist, norm, n1B[3], n2B[3], sp, dist1Sq, dist2Sq, dist3Sq, dist4Sq;
  double ltot, th1, th2, db12[3];
  
  for (kk=0; kk < 3; kk++)
    b12[kk] = b2[kk]-b1[kk];

  distbb = calc_norm(b12);
  if (l1+l2 < distbb)
    return DISTMAX;
  //angle = calc_angle(params[0], params[1]);
  th1 = acos((-Sqr(l2)+Sqr(l1)+Sqr(distbb))/(2.0*l1*distbb));
  th2 = acos((-Sqr(l1)+Sqr(l2)+Sqr(distbb))/(2.0*l2*distbb));
  angle = 2.0*acos(0.0)-th1-th2;
#if 0
  printf("distbb=%f l1=%f l2=%f\n", distbb, l1, l2);
  printf("th1=%f th2=%f angle=%f\n", th1*180/3.14, th2*180/3.14, angle*180/3.14);
#endif
  ////delx = tan(angle/2.0)*Dhc/2.0; 
  delx = 0.0;
  th2=2.0*acos(0.0)-th2;
  Lhc = (Dhc + delx);
 
  Lhc1 = l1 + delx;
  Lhc2 = l2 + delx; 
   //printf("l1=%f l2=%f th1=%f th2=%f phi=%f\n", l1, l2, th1, th2, phi);
  n1B[0] = sin(th1)*cos(phi);
  n1B[1] = sin(th1)*sin(phi);
  n1B[2] = cos(th1);
  body2labR(n1B, n1, Ro);

  n2B[0] = sin(th2)*cos(phi);
  n2B[1] = sin(th2)*sin(phi);
  n2B[2] = cos(th2);
  body2labR(n2B, n2, Ro);
#if 0
  n2[0] = -n1[0];
  n2[1] = -n1[1];
  n2[2] = -n1[2];
#endif  
  fact = Lhc1*0.5 - delx*0.5;
#if 1
  pos1[0] = b1[0]+n1[0]*fact;
  pos1[1] = b1[1]+n1[1]*fact;
  pos1[2] = b1[2]+n1[2]*fact;
#else
  pos1B[0]=n1B[0]*fact;
  pos1B[1]=n1B[1]*fact;
  pos1B[2]=n1B[2]*fact;
  body2lab(pos1B, pos1, b1, Ro);
#endif
  //printf("n=%f %f %f pos1=%f %f %f (norm=%f)\n", n1[0], n1[1], n1[2], pos1[0], pos1[1], pos1[2], calc_norm(pos1));
  /* com2 = {-(X0 D /2 - delx*0.5), 0, 0};*/
  //fact = Dhc/2.0 - delx*0.5;
  fact = Lhc2*0.5 - delx*0.5;
  
#if 1
  pos2[0] = b2[0]+n2[0]*fact;
  pos2[1] = b2[1]+n2[1]*fact;
  pos2[2] = b2[2]+n2[2]*fact;
#else
  pos2B[0]=n2B[0]*fact;
  pos2B[1]=n2B[1]*fact;
  pos2B[2]=n2B[2]*fact;
  body2lab(pos2B, pos2, b2, Ro);
#endif
#if 0
  printf("P1=P2 P1:%f %f %f P2:%f %f %f\n",  b1[0]+n1[0]*Lhc1, b1[1]+n1[1]*Lhc1, b1[2]+n1[2]*Lhc1,
                                             b2[0]+n2[0]*Lhc2, b2[1]+n2[1]*Lhc2, b2[2]+n2[2]*Lhc2);
#endif
  dist=DISTMAX;
  disttot=0.0;

  for(k=0; k<22; k++)
    {
      /* check overlap here */
      drp[0] = P[i0+k].x - pos1[0];
      drp[1] = P[i0+k].y - pos1[1];
      drp[2] = P[i0+k].z - pos1[2];
      dontcheck=0;
      sp = scalProd(drp, n1);
      if (sp < 0.0)
	{
	  dist1Sq = Sqr(fabs(sp) - Lhc1 * 0.5);
	  if (dist1Sq < dist)
	    dist=dist1Sq;
	}
      vectProdVec(drp, n1, vp);
      norm=calc_norm(vp);
      dist2Sq = Sqr(norm - Dhc*0.5);
      if (dist2Sq < dist)
	dist = dist2Sq;
      //printf("dist2Sq=%f\n", dist2Sq);
      drp[0] = P[i0+k].x - pos2[0];
      drp[1] = P[i0+k].y - pos2[1];
      drp[2] = P[i0+k].z - pos2[2];
      sp = scalProd(drp, n2);
      if (sp < 0.0)
	{
  	  dist3Sq = Sqr(fabs(sp) - Lhc2 * 0.5); 
	  if (dist3Sq < dist)
	    dist = dist3Sq;
	}
      vectProdVec(drp, n2, vp);
      norm=calc_norm(vp);
      dist4Sq = Sqr(norm - Dhc*0.5);   
      if (dist < dist4Sq)
	dist = dist4Sq;
      disttot+=dist;
    }
  if (disttot == DISTMAX)
    {
      printf("We have a problem, dist=-1\n");
      exit(-1);
    }
  return disttot;
}
#if 0
void conjgrad_grad(double *angs, double *grad)
{
  double r1[3], r2[3], xp[6], dd[3], jac[6][4];
  double sin0, sin1, cos0, cos1, sin2, sin3, cos2, cos3;
  int kk, k1, k2;

  sin0 = sin(angs[0]);
  sin1 = sin(angs[1]);
  cos0 = cos(angs[0]);
  cos1 = cos(angs[1]);
  sin2 = sin(angs[2]);
  sin3 = sin(angs[3]);
  cos2 = cos(angs[2]);
  cos3 = cos(angs[3]);
  xp[0] = axa[icg]*cos0*sin1;
  xp[1] = axb[icg]*sin0*sin1;
  xp[2] = axc[icg]*cos1;
  xp[3] = axa[jcg]*cos2*sin3;
  xp[4] = axb[jcg]*sin2*sin3;
  xp[5] = axc[jcg]*cos3;

  body2lab(icg, xp, r1, rA, RtA);
  body2lab(jcg, &xp[3], r2, rB, RtB);

  for (kk = 0; kk < 3; kk++)
    dd[kk] = -2.0*(r2[kk] - r1[kk]);   
  printf("[CG grad] i=%d j=%d dist=%.15G\n", icg, jcg, calc_norm(dd)/2.0);
  /* NOTA: jac[][] è lo jacobiano del cambio di coordinate 
   * (theta1,phi1,theta2,phi2)->(x1,y1,z1,x2,y2,z2) cioè 
   * \frac{\delta(x1,y1,z1,x2,y2,z2)}{\delta(theta1,phi1,theta2,phi2)}*/
  jac[0][0] = -axa[icg]*sin0*sin1;
  jac[0][1] = axa[icg]*cos0*cos1;
  jac[0][2] = 0.0;
  jac[0][3] = 0.0;
  jac[1][0] = axb[icg]*cos0*sin1;
  jac[1][1] = axb[icg]*sin0*cos1; 
  jac[1][2] = 0.0;
  jac[1][3] = 0.0;
  jac[2][0] = 0.0;
  jac[2][1] = -axc[icg]*sin1; 
  jac[2][2] = 0.0;
  jac[2][3] = 0.0;

  jac[3][0] = 0.0;
  jac[3][1] = 0.0;
  jac[3][2] = -axa[jcg]*sin2*sin3;
  jac[3][3] = axa[jcg]*cos2*cos3;
  jac[4][0] = 0.0;
  jac[4][1] = 0.0;
  jac[4][2] = axb[jcg]*cos2*sin3;
  jac[4][3] = axb[jcg]*sin2*cos3;
  jac[5][0] = 0.0;
  jac[5][1] = 0.0;
  jac[5][2] = 0.0;
  jac[5][3] = -axc[jcg]*sin3;

  for (kk = 0; kk < 2; kk++)
    {
      grad[kk] = 0.0;
      grad[kk+2] = 0.0;
      for (k1 = 0; k1 < 3; k1++)
	{
	  grad[kk] += -dd[k1]*jac[k1][kk];
	  grad[kk+2] += dd[k1]*jac[k1+3][kk+2];
	}
    }
  if (angs[1]+grad[1] > pi)
    grad[1] = 2.0*pi - grad[1];  
  if (angs[1]+grad[1] < 0.0)
    grad[1] = -grad[1];
  if (angs[3]+grad[3] > pi)
    grad[3] = 2.0*pi - grad[3];  
  if (angs[3]+grad[3] < 0.0)
    grad[3] = -grad[3];
   //printf("grad=(%.15G,%.15G,%.15G,%.15G\n", grad[0], grad[1], grad[2], grad[3]);
}
double rDcg[3];
int conjgrad(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double []))
{
  const int ITMAX=10;
  const double EPS = 1.0E-12;
  int j,its;
  double gg,gam,fp,dgg, g[8],h[8],xi[8], dist, distold;
  /* Initializations.*/
  fp=(*func)(p); 
  (*dfunc)(p,xi);
  //distold = Sqr(rDcg[0]-p[0])+Sqr(rDcg[1]-p[1])+Sqr(rDcg[2]-p[2]);

  for (j=0; j < n; j++) 
    {
      g[j] = -xi[j];
      xi[j] = h[j] = g[j];
    }
  for (its = 0; its < ITMAX; its++) 
    { 
      *iter=its;
      linmin(p, xi, n, fret, func); /* Next statement is the normal return: */
      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) 
	{
	  return 0;
	}
      fp= *fret;
      (*dfunc)(p,xi);
      dgg = gg = 0.0;
      for (j = 0; j < n; j++)
	{
	  gg += g[j]*g[j];
	  /* dgg += xi[j]*xi[j]; */ /* This statement for Fletcher-Reeves.*/
	  dgg += (xi[j]+g[j])*xi[j]; /* This statement for Polak-Ribiere.*/
	}
      if (gg == 0.0) 
	{ 
	  /* Unlikely. If gradient is exactly zero then
	     we are already done. */
	  return 0;
	}
      gam=dgg/gg;
      for (j = 0; j < n; j++) 
	{
	  g[j] = -xi[j];
	  xi[j]=h[j]=g[j]+gam*h[j];
	}
    }
  //printf("ERROR: Too many iterations in frprmn *fret=%.15G\n", *fret);
  //exit(-1);
  return 0;
}
#endif
void print_matrix(double M[3][3], int n)
{
  int k1, k2;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}
void print_usage(void)
{
  printf("BCAparam [-o <params_file> | -e <end2end_file> | --mglmode/-m | --firstframe/-f <first_frame> \n");
  printf("   | --lastframe/-l <last_frame> | --mglfn|-mf <mgl_file_name> | --Prad/-Pr <phospate radius>\n"); 
  printf("   | --l12min/-l1m <l1_over_l2_min> | --l12max|-l12M <l1_over_l2_max> | --ltotmin/-ltm <ltotmin>\n");
  printf("   | --ltotmax|-ltM <ltotmax> | -nl <mesh_points_for_l> | -nlt <mesh_points_for_ltot>\n");
  printf("   | -nphi <mesh_points_for_phi> | --outframes/-of <output frames> | -Lx/y/z <box_size_along_x/y/z>\n");
  printf(" | --fixbroken|-fb ] <pdb_file>\n");
  exit(0);
}
void parse_param(int argc, char** argv)
{
  int cc=1;
  int extraparam = 0;  

  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"-o" ))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  strcpy(bcparfile,argv[cc]);
	} 
      else if (!strcmp(argv[cc],"-e"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  strcpy(e2efile,argv[cc]);
	}
      else if (!strcmp(argv[cc],"--mglfn")|!strcmp(argv[cc],"-mf"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  strcpy(mglfn,argv[cc]);
	}
      else if (!strcmp(argv[cc],"--mglmode") || !strcmp(argv[cc],"-m"))
	{
	  mglout=1;
	}
      else if (!strcmp(argv[cc],"--firstframe") || !strcmp(argv[cc],"-f"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  firstframe = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--lastframe") || !strcmp(argv[cc],"-l"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  lastframe = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--l1l2min") || !strcmp(argv[cc],"-l12m"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  l1l2min = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--l1l2max") || !strcmp(argv[cc],"-l12M"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  l1l2max = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--diameter") || !strcmp(argv[cc],"-D"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  Dhc = atof(argv[cc]);
	}

      else if (!strcmp(argv[cc],"--ltotmin") || !strcmp(argv[cc],"-ltm"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  ltotmin = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--ltotmax") || !strcmp(argv[cc],"-ltM"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  ltotmax = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"-Lx"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  Lx = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"-Ly"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  Ly = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"-Lz"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  Lz = atof(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--fixbroken")||!strcmp(argv[cc],"-fb"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  fixbroken = 1;
	}
      else if (!strcmp(argv[cc],"-nl"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  nl = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"-nlt"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  nlt = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"-nphi"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  nphi = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--outframes") || !strcmp(argv[cc],"-of") )
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  outframes = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--Prad") || !strcmp(argv[cc],"-Pr"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  Prad = atof(argv[cc]);
	}
      else if (cc == argc || extraparam == 1)
	print_usage();
      else if (extraparam == 0)
	{ 
	  extraparam = 1;
	  strcpy(infile,argv[cc]);
	}
      else
	print_usage();
      cc++;
    }
}

int main(int argc, char *argv[])
{
  int ibeg, iend, opt, res, numP, P_count, i, k, found_one=0, first, kk;
  double pi, x, y, z, l, m, norm, comx, comy, comz, distbest, l1best, l2best, phi, dphi, l1, l2;
  double dl, ltot, l1min, l1max, del_l1, del_ltot, sp;
  double e2eav, xv[3], yv[3], zv[3], Ro[3][3], b1[3], b2[3], n1[3], n2[3], pos1[3], pos2[3];
  double Lgx, Lgy, Lgz, n1best[3], n2best[3], pos1best[3], pos2best[3], Lmgl=0.0;
  double cc=0, angle, angleav, distav, l1av, l2av, dist, anglebest, delb[3], e2ebest, Lhc1best, Lhc2best,
	 mglcomx, mglcomy, mglcomz;
  pi = 2.0*acos(0.0);
#if 0
  angle=PI*(180.0 - atof(argv[1]))/180.0;
  n1[0] = cos(angle);
  n1[1] = sin(angle);
  n1[2] = 0.0;
  n2[0] = -1.0;
  n2[1] = 0.0;
  n2[2] = 0.0;
#endif
  FILE *buffer, *in, *e2e;
#if 0
  while ((opt = getopt (argc, argv, "i:e:o:")) != -1)
    {
    switch (opt)
      {
      case 'i':
	printf ("Input file: \"%s\"\n", optarg);
	strcpy(infile, optarg);
	break;
      case 'e':
	printf ("Output file: \"%s\"\n", optarg);
	strcpy(e2efile,optarg);
	break;
      case 'o':
	printf ("bcparfile: \"%s\"\n", optarg);
	strcpy(bcparfile,optarg);
	break;
	//	case 'h':
	//	printf("./nameprogram -i input_trajectori.pdb -e output_end2end -o output_angle");
      }
    }
#endif
  dl=1./6.;
  ltotmin = 36.0-dl;
  ltotmax = 44.0-dl;
     
  parse_param(argc, argv);
  //buffering P positions into a file
  in= fopen(infile, "r");
  buffer = fopen("buffer.pdb", "w");
  if (mglout)
    mglfile = fopen(mglfn, "w+");
  while ( fgets(string, 100, in) != NULL )
    {
      if(string[13] == 'P'  )
	{
	  fprintf(buffer, "%s", string);
	}
    }
  fclose(buffer);
  fclose(in);
  frame = 0;
  buffer = fopen("buffer.pdb", "r");
  e2e = fopen(e2efile, "w+");  
  numP=0;
  dist_aver=0.0;
  while ( !feof(buffer) )
    {
    fscanf(buffer, "%22c %d %lf %lf %lf %lf %lf\n", dummystr, &res, &x, &y, &z, &l, &m);
    //printf("x=%f %f %f\n", x, y, z);
    if ((res-1)%24==2-1)       
      { 
	PA1_1.x=x; PA1_1.y=y; PA1_1.z=z;
      }
    else if ((res-1)%24==12-1) 
      { 
	PA2_1.x=x; PA2_1.y=y; PA2_1.z=z; 
      }
    else if ((res-1)%24==14-1) 
      {
       	PB2_1.x=x; PB2_1.y=y; PB2_1.z=z; 
      }
    else if ((res-1)%24==24-1) 
      { 
	PB1_1.x=x; PB1_1.y=y; PB1_1.z=z;
      }

    if((res-1)%24==24-1)
      {
	bar1_1.x=(PA1_1.x+PB1_1.x)/2.; bar1_1.y=(PA1_1.y+PB1_1.y)/2.; bar1_1.z=(PA1_1.z+PB1_1.z)/2.;
	bar2_1.x=(PA2_1.x+PB2_1.x)/2.; bar2_1.y=(PA2_1.y+PB2_1.y)/2.; bar2_1.z=(PA2_1.z+PB2_1.z)/2.;

	dist_ist=pow ( (bar1_1.x-bar2_1.x)*(bar1_1.x-bar2_1.x) + (bar1_1.y-bar2_1.y)*(bar1_1.y-bar2_1.y) + (bar1_1.z-bar2_1.z)*(bar1_1.z-bar2_1.z) , 0.5);

	frame+=1.0;
	//printf("res=%d mod=%d\n", res,(res-1)%24);
	dist_aver = dist_aver+dist_ist;
	fprintf(e2e, "%f %f %f\n", frame, dist_ist, dist_aver/frame);
      }
    numP++;
    }
  e2eav=dist_aver/frame;
  printf("\nafter %f steps the average e2e is:  %f\n", frame, dist_aver/frame);
  fclose(e2e);
  rewind(buffer);

  P = malloc(sizeof(struct vector)*numP); 
  frame = 0;
  while ( !feof(buffer))
    {
      fscanf(buffer, "%22c %d %lf %lf %lf %lf %lf\n", a, &res, &x, &y, &z, &l, &m);
      P[P_count].x = x;
      P[P_count].y = y;
      P[P_count].z = z;
      P_count++;
      //if (P_count%10000==0) printf("P_count=%d\n", P_count);
      //printf("%f %f %f \n", P_x[P_count], P_y[P_count], P_z[P_count]);
    }
  fclose(buffer);
  srand48(145); /* seeding */
  distav=angleav=l1av=l2av=0.0;
  if (firstframe > 0)
    ibeg = 22*firstframe;
  else
    ibeg = 0;

  /* un frame è un duplex */
  if (lastframe > 0 && lastframe > firstframe)
    {
      iend = 22*lastframe;
      if (iend > numP)
	iend = numP;
    }
  else
    iend=numP;
  if (mglout)
    printf("[mglmode] ibeg=%d iend=%d\n", ibeg, iend);

  printf("Finding best parameters for BCA model\n");
  printf("# duplexes is %d\n", (iend-ibeg)/22);
  printf("l1l2 range=(%f-%f) ltot range=(%f,%f) dl=%d dlt=%d dphi=%d\n",
	 l1l2min, l1l2max, ltotmin, ltotmax, nl, nlt, nphi);

  if (mglout)
    {
      cc = 0;
      mglcomx = mglcomy = mglcomz = 0.0;
      for (i=ibeg; i < iend; i=i+1)
	{
	  mglcomx += P[i].x;
	  mglcomy += P[i].y;
	  mglcomz += P[i].z;
	  cc++;
	}
      /* quest è il centro di massa di tutti i fosfori che si considerano */
      mglcomx /= cc;
      mglcomy /= cc;
      mglcomz /= cc;
      Lmgl = 0.0;
      Lgx=Lgy=Lgz=0.0;
      for (i=ibeg; i < iend; i=i+1)
	{
	  if (fabs(P[i].x - mglcomx) > Lgx)
	    Lgx = fabs(P[i].x - mglcomx);
	  if (fabs(P[i].y - mglcomy) > Lgy)
	    Lgy = fabs(P[i].y - mglcomy);
	  if (fabs(P[i].z - mglcomz) > Lgz)
	    Lgz = fabs(P[i].z - mglcomz);
	}
      
      Lmgl = 2.0*max3(Lgx,Lgy,Lgz);
      if (fixbroken && (Lx<=0 || Ly<=0 || Lz<=0))
	{
	  Lx = Lgx*2.0;
	  Ly = Lgy*2.0;
	  Lz = Lgz*2.0;
	  printf("[WARNING] You want me to fix broken duplexes but you did not supply box dimensions, hence\n");
	  printf("I tried to guess box dimensions and they turn out to be: %f %f %f\n", Lx, Ly, Lz);
	}
      fprintf(mglfile, ".Vol: %f\n", pow(Lmgl,3.0));
    }

  cc=0;
  for (i=ibeg; i < iend; i=i+22)
    {
      if ((i/22)%outframes==0 && i > ibeg) 
	printf("# duplex=%d/%d\n", i/22, numP/22);
      //center of mass of terminal Phosphate pairs 
      bar1.x = (P[i].x+P[i+21].x)*0.5;
      bar1.y = (P[i].y+P[i+21].y)*0.5;
      bar1.z = (P[i].z+P[i+21].z)*0.5;

      bar2.x = (P[i+10].x+P[i+11].x)*0.5;
      bar2.y = (P[i+10].y+P[i+11].y)*0.5;
      bar2.z = (P[i+10].z+P[i+11].z)*0.5;

      /* center of mass */
      found_one=0;
      comx=comy=comz=0.0;	
      for(k=0; k<22; k++)
	{
	  comx += P[i+k].x;
	  comy += P[i+k].y;
	  comz += P[i+k].z;
	  if (mglout==1)
	    {
	      fprintf(mglfile, "%f %f %f @ %f C[green]\n", P[i+k].x-mglcomx, P[i+k].y-mglcomy, P[i+k].z-mglcomz, Prad);
	    }
	}
      comx /= 22.;
      comy /= 22.;
      comz /= 22.;
      
      pi = 2.0*acos(0.0);
      /* swearch among all possible values of l_1, l_2 and phi! */
      b1[0] = bar1.x;
      b1[1] = bar1.y;
      b1[2] = bar1.z;
      b2[0] = bar2.x;
      b2[1] = bar2.y;
      b2[2] = bar2.z;

      zv[0] = bar2.x-bar1.x;
      zv[1] = bar2.y-bar1.y;
      zv[2] = bar2.z-bar1.z;
      norm = calc_norm(zv);
      for (k=0; k < 3; k++)
	zv[k] /= norm;

      xv[0] = bar2.x - comx;
      xv[1] = bar2.y - comy;
      xv[2] = bar2.z - comz;

      sp = scalProd(xv, zv);
      for (k=0; k < 3; k++)
	xv[k] = xv[k] - zv[k]*sp;

      norm = calc_norm(xv);
      for (k=0; k < 3; k++)
	xv[k] /= norm;

      vectProdVec(zv, xv, yv);
      Ro[0][0] = xv[0];
      Ro[0][1] = xv[1];
      Ro[0][2] = xv[2];
      Ro[1][0] = yv[0];
      Ro[1][1] = yv[1];
      Ro[1][2] = yv[2];
      Ro[2][0] = zv[0];
      Ro[2][1] = zv[1];
      Ro[2][2] = zv[2];
      //print_matrix(Ro,3);
      del_ltot = (ltotmax-ltotmin)/((double)nlt);
      //del_l1 = (l1max-l1min)/20.;
      dphi = 2.0*pi/((double)nphi);
      /* le lunghezze sono in angstrom */
      //l1min = 12;
      //l1max = 20;

      first=1;
      for (phi=0; phi < 2.0*pi; phi += dphi)
	{
	  for (ltot=ltotmin; ltot < ltotmax; ltot += del_ltot)
	    {
	      l1min = l1l2min*ltot/2.0;
	      l1max = l1l2max*ltot/2.0;
	      del_l1 =  (l1max-l1min)/((double)nl);
	      for (l1 = l1min; l1 < l1max; l1 += del_l1)
		{
		  dist=dist_func(i, l1, ltot-l1, phi, Ro, b1, b2, pos1, n1, pos2, n2);
		  //printf("l1=%f del_l1=%f\n", l1, del_l1);
		  if (first || dist < distbest)
		    {
		      first=0;
		      l1best = l1;
		      l2best = ltot-l1best;
		      for (kk=0; kk < 3; kk++)
			delb[kk] = b2[kk]-b1[kk];
		      e2ebest = calc_norm(delb);
		      distbest = dist;
		      for (kk=0; kk < 3; kk++)
			{
			  n1best[kk] = n1[kk];
			  n2best[kk] = n2[kk];
			  pos1best[kk] = pos1[kk];
			  pos2best[kk] = pos2[kk];
			}
		      Lhc1best = Lhc1;
		      Lhc2best = Lhc2;
		    }
		}
	    }
	}
      anglebest = 180.*acos((-Sqr(e2ebest)+Sqr(l1best)+Sqr(l2best))/(2.0*l1best*l2best))/acos(0.0)/2.0;
      distav += distbest;
      if (l1best < l2best)
	{
	  l1av += l1best;
	  l2av += l2best;
	}
      else
	{
	  l1av += l2best;
	  l2av += l1best;
	}
      if (mglout==1)
	{
	  fprintf(mglfile, "%f %f %f %f %f %f @ %f %f C[red]\n", pos1best[0]-mglcomx, pos1best[1]-mglcomy, pos1best[2]-mglcomz, n1best[0], n1best[1], n1best[2], Dhc/2.0, Lhc1best);
	  fprintf(mglfile, "%f %f %f %f %f %f @ %f %f C[red]\n", pos2best[0]-mglcomx, pos2best[1]-mglcomy, pos2best[2]-mglcomz, n2best[0], n2best[1], n2best[2], Dhc/2.0, Lhc2best);
	}

      //printf("anglebest=%f\n", anglebest);
      angleav += anglebest;
      cc++;
    }
  if (mglout)
    fclose(mglfile);
  printf("l1=%.15G l2=%.15G theta_b=%.15G\n", l1av/cc, l2av/cc, angleav/cc);

}

