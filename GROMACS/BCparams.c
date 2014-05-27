#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#define Sqr(x) ((x)*(x))
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
double frame, dist_ist, dist_aver;
double Lhc, Dhc=2.0;
double PI;

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

double dist_func(int i0, double l1, double l2, double phi, double Ro[3][3], double b1[3], double b2[3])
{
  /* params = {L_1, L_2, theta} */
  int k, i, dontcheck;
  double fact, pos1[3], pos2[3], delx, angle, pos1B[3], pos2B[3];
  double vp[3], rp[3], drp[3];
  double dist, norm, n1B[3], n2B[3], n1[3], n2[3], sp, dist1Sq, dist2Sq, dist3Sq, dist4Sq;
  double ltot, th1, th2, db12[3];
  
  ltot = l1+l2;
  /* use cosine theorem to calculate angle */
  //angle = calc_angle(params[0], params[1]);
  delx = tan(angle/2.0)*Dhc/2.0; 
  Lhc = (Dhc + delx);

  th1 = acos((-Sqr(l2)+Sqr(l1)+Sqr(ltot))/(2.0*l1*ltot));
  th2 = 2.0*acos(0.0)-acos((-Sqr(l1)+Sqr(l2)+Sqr(ltot))/(2.0*l2*ltot));

  //printf("l1=%f l2=%f th1=%f th2=%f phi=%f\n", l1, l2, th1, th2, phi);
  n1B[0] = sin(th1)*cos(phi);
  n1B[1] = sin(th1)*sin(phi);
  n1B[2] = cos(th1);
  body2labR(n1B, n1, Ro);

  n2B[0] = sin(th2)*cos(phi);
  n2B[1] = sin(th2)*sin(phi);
  n2B[2] = cos(th2);
  body2labR(n2B, n2, Ro);

  fact= Dhc*0.5 - delx*0.5;

  pos1B[0]=n1B[0]*fact;
  pos1B[1]=n1B[1]*fact;
  pos1B[2]=n1B[2]*fact;
  body2lab(pos1B, pos1, b1, Ro);

  //printf("n=%f %f %f pos1=%f %f %f (norm=%f)\n", n1[0], n1[1], n1[2], pos1[0], pos1[1], pos1[2], calc_norm(pos1));
  /* com2 = {-(X0 D /2 - delx*0.5), 0, 0};*/
  fact = -(Dhc/2.0 - delx*0.5);
  pos2B[0]=n2B[0]*fact;
  pos2B[1]=n2B[1]*fact;
  pos2B[2]=n2B[2]*fact;
  body2lab(pos2B, pos2, b2, Ro);

  dist=1000000;

  for(k=0; k<22; k++)
    {
      /* check overlap here */
      drp[0] = P[i0+k].x - pos1[0];
      drp[1] = P[i0+k].y - pos1[1];
      drp[2] = P[i0+k].z - pos1[2];
      dontcheck=0;
      sp = scalProd(drp, n1);
      if (sp > 0.0)
	{
	  dist1Sq = Sqr(fabs(sp) - Lhc * 0.5);
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
      if (sp > 0.0)
	{
  	  dist3Sq = Sqr(fabs(sp) - Lhc * 0.5); 
	  if (dist3Sq < dist)
	    dist = dist3Sq;
	}
      vectProdVec(drp, n2, vp);
      norm=calc_norm(vp);
      dist4Sq = Sqr(norm - Dhc*0.5);   
      if (dist < dist4Sq)
	dist = dist4Sq;
    }
  if (dist == 1000000)
    {
      printf("We have a problem, dist=-1\n");
      exit(-1);
    }
  return dist;
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

int main(int argc, char *argv[])
{
  int opt, res, numP, P_count, i, k, found_one=0, first, kk;
  double pi, x, y, z, l, m, norm, comx, comy, comz, distbest, l1best, l2best, phi, dphi, l1, l2;
  double dl, ltot, l1min, l1max, ltotmin, ltotmax, del_l1, del_ltot, sp;
  double e2eav, xv[3], yv[3], zv[3], Ro[3][3], b1[3], b2[3];
  double cc, angle, angleav, distav, l1av, l2av, dist, anglebest, delb[3], e2ebest;
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
  //buffering P positions into a file
  in= fopen(infile, "r");
  buffer = fopen("buffer.pdb", "w");
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
      if (P_count%10000==0) printf("P_count=%d\n", P_count);
      //printf("%f %f %f \n", P_x[P_count], P_y[P_count], P_z[P_count]);
    }
  fclose(buffer);
  srand48(145); /* seeding */
  distav=angleav=l1av=l2av=0.0;
  for (i=0; i < numP; i=i+22)
    {
      if (i%100==0) 
	printf("i=%d/%d\n", i, numP);
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
	}
      comx /= 22.;
      comy /= 22.;
      comz /= 22.;

      pi = 2.0*acos(0.0);
      /* swearch among all possible values of l_1, l_2 and phi! */
      dl = 1./6.;
      /* le lunghezze sono in angstrom */
      ltotmin = 38.0-dl;
      ltotmax = 42.0-dl;
      l1min = 18.0;
      l1max = 22.0;

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
      Ro[1][0] = xv[1];
      Ro[2][0] = xv[2];
      Ro[0][1] = yv[0];
      Ro[1][1] = yv[1];
      Ro[2][1] = yv[2];
      Ro[0][2] = zv[0];
      Ro[1][2] = zv[1];
      Ro[2][2] = zv[2];
      del_ltot = (ltotmax-ltotmin)/10.;
      del_l1 = (l1max-l1min)/20.;
      dphi = 2.0*pi/20;
      first=1;
      for (phi=0; phi < 2.0*pi; phi += dphi)
	{
	  for (ltot=ltotmin; ltot < ltotmax; ltot +=del_ltot)
	    {
	      for (l1 = l1min; l1 < l1max; l1 += del_l1)
		{
		  dist=dist_func(i, l1, ltot-l1, phi, Ro, b1, b2);
		  if (first || dist < distbest)
		    {
		      first=0;
		      l1best = l1;
		      l2best = ltot-l1best;
		      for (kk=0; kk < 3; kk++)
			delb[kk] = b2[kk]-b1[kk];
		      e2ebest = calc_norm(delb);
		      distbest = dist;
		    }
		}
	    }
	}
      anglebest = 180.*acos((-Sqr(e2ebest)+Sqr(l1best)+Sqr(l2best))/(2.0*l1best*l2best))/acos(0.0)/2.0;
      distav += distbest;
      l1av += l1best;
      l2av += l2best;
      angleav += anglebest;
      cc++;
    }
  printf("l1=%.15G l2=%.15G theta_b=%.15G\n", l1av/cc, l2av/cc, angleav/cc);

}

