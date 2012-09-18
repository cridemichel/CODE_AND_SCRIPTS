#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define Sqr(X) ((X)*(X))

/* <<< COSTANTI DA CAMBIARE >>> */
const double m=1.0; // massa 
double temp = 0.1;// temperatura
const double I=0.1;//momentum of inertia

double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
/* ------------------------------ */
void vectProdVec(double *A, double *B, double *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1]; 
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
}
void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;

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
#if 0
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
#endif
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
#if 0
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    Rt[k1][k2]=R[k2][k1];
  for (k1=0; k1 < 3 ; k1++)
    for (k2=0; k2 < 3 ; k2++)
    R[k1][k2]=Rt[k1][k2];
#endif

  //printf("calc_norm R[2]=%f vp=%f\n", calc_norm(R[2]), scalProd(R[1],R[2]));
}

/* ============================ >>> ranf <<< =============================== */
double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return rand() / ( (double) RAND_MAX );
}

double gauss(void)
{
  
  /* 
     Random variate from the standard normal distribution.
     
     The distribution is gaussian with zero mean and unit variance.
     REFERENCE:                                                    
                                                                
     Knuth D, The art of computer programming, (2nd edition        
     Addison-Wesley), 1978                                      
                                                                
     ROUTINE REFERENCED:                                           
                                                                
     COORD_TYPE ranf()                                  
     Returns a uniform random variate on the range zero to one  
  */

  double  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0;i < 12; i++)
    {
      sum = sum + ranf();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}

void angvel(double *wx, double *wy, double* wz, int Nm)
{
  int i;
  double inert;                 /* momentum of inertia of the molecule */
  double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, Mtot;
  
  Mtot = m; /* total mass of molecule */

  inert = I; /* momentum of inertia */
 
  mean = 3.0*temp / inert;

  for (i = 0; i < Nm; i++)
    {
      xisq = 1.0;
      
      while (xisq >= 1.0)
	{
	  xi1  = ranf() * 2.0 - 1.0;
	  xi2  = ranf() * 2.0 - 1.0;
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
      
      /* Choose the magnitude of the angular velocity
         NOTE: consider that it is an exponential distribution 
	 (i.e. Maxwell-Boltzmann, see Allen-Tildesley pag. 348-349)*/
      
      osq   = - mean * log(ranf());
      o     = sqrt(fabs(osq));
      ox    = o * ox;
      oy    = o * oy;
      oz    = o * oz;
      *wx = ox;
      *wy = oy;
      *wz = oz;
    }
}
void evalfourth(double *r1, double *r2, double *r3, double sigma)
{
  double pi;
  double etaaa;
  double ankat, gekat, phi;
  int iv[3]={1,2,0}, jv[3]={2,0,1};
  double d, el, em, dn, dip[3], vn[3], vecn, vec;
  int kk, lc, mc;
  etaaa = 109.471/2.0;
  pi = acos(0)*2.0;
  phi = etaaa*pi/180.0;
  //fprintf(stderr,"radius:%f phi=%f\n", radius, phi);
  ankat=(sigma/2.0)*cos(phi);
  gekat=(sigma/2.0)*sin(phi);
  el = ankat;
  em = gekat;
  dn = 0.0;
  for (kk=0; kk < 3; kk++)
    {
      dip[kk] = r1[kk]+r2[kk];
      dn += dip[kk]*dip[kk];
    }
  dn = sqrt(dn);
  //fprintf(stderr,"dn=%f\n", dn);
  for (kk=0; kk < 3; kk++)
    dip[kk]=dip[kk]/dn;
  d = 0.0;
  for (kk=0; kk < 3; kk++)
    {
      lc = iv[kk];
      mc = jv[kk];
      vn[kk] = r1[lc]*r2[mc] - r1[mc]*r2[lc];
      d = d + vn[kk]*vn[kk];
    }
  d = sqrt(d);
  for (kk=0; kk < 3; kk++)
    {
      vn[kk] /= d;
     
      vecn = vn[kk]*em;
      vec = -dip[kk]*el;
      //fprintf(stderr,"em: %f el:%f vec: %f vecn: %f\n", em,el,vec, vecn);
      //r3[kk] = vec - vecn; 
      r3[kk] = vec+vecn;
    }
}
double *vx, *vy, *vz, *wx, *wy, *wz;
double calc_energy(int nummol)
{
  int i, k1;
  double K,wt[3];
  K = 0;
  for (i=0; i < nummol; i++)
    {
	 /* calcola tensore d'inerzia e le matrici delle due quadriche */
	  K += m*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    K += Sqr(wt[k1])*I;
    }
  K *= 0.5;
  return K;
}
int nummolmax;
double R[3][3];
char line[4096];
int main(int argc, char **argv)
{
  FILE* f, *f2;
  int i, kk, nummol, junki1, junki2, junki3, passi;
  double sf, radius, junkd1, junkd2, junkd3, trad, patchRad, patchPos, tlen;
  double r1[3], r2[3], r3[3], r21[3], r31[3], nn[3], norm, rO[3];
  double K, r4[3], bx, by, bz, r41[3], sumx, sumy, sumz;
  //double  rTemp;
  //printf("Immetti la  temperatura: ");
  //scanf("%lf",&temp);
  f2 = fopen(argv[2],"w");
  f = fopen(argv[1], "r");
  fscanf(f, "%[^\n]\n", line);
  fscanf(f, "%d %lf %lf %lf %lf %lf\n", &nummol, &bx, &trad, &tlen, &patchRad, &patchPos);
  by = bz = bx;
  printf("NUMMOL: %d\n", nummol);
#if 0
  vx = malloc(sizeof(double)*nummol);
  vy = malloc(sizeof(double)*nummol);
  vz = malloc(sizeof(double)*nummol);
  wx = malloc(sizeof(double)*nummol);
  wy = malloc(sizeof(double)*nummol);
  wz = malloc(sizeof(double)*nummol);
#endif		
  if (argc==4)	
   nummol=atoi(argv[3]);
#ifdef MKCNF
  fprintf(f2,"ensembleMC:0\n@@@\n");
#endif
  fprintf(f2,"parnum: %d\n", nummol);
  //fprintf(f2,"parnumA: %d\n", nummol);
  //fprintf(f2,"T:%.15G\n",temp);
/* HARD CYLINDERS WITH PARTCHES:
tol: 0
ninters: 3
nintersIJ: 0
ntypes: 1
saveBonds: 0
maxbondsSaved: -1
@@@
850
1 1 1
2 2 2
1 1 1 1 2 0
2 0
1.15 0 0 0.5
-1.15 0 0 0.5
0 0 0 0 1 0 0 100000
0 0 0 1 1 0 0 100000
0 1 0 1 1 0 0 100000
   
*/
  fprintf(f2,"ninters: 3\nnintersIJ: 0\nntypes: 1\nsaveBonds: 0\nmaxbondsSaved: -1\n@@@\n");
  fprintf(f2,"%d\n", nummol);	
  fprintf(f2,"%G %G %G\n", tlen*0.5, trad, trad);	
  fprintf(f2,"2 2 2\n");
  fprintf(f2,"1 1 1 1 2 0\n");
  fprintf(f2,"2 0\n");
  fprintf(f2,"%.15G 0 0 %.15G\n", tlen*0.5+patchPos, patchRad*2.0);
  fprintf(f2,"%.15G 0 0 %.15G\n", -tlen*0.5-patchPos, patchRad*2.0);
  fprintf(f2,"0 0 0 0 1 0 0 100000\n");
  fprintf(f2,"0 0 0 1 1 0 0 100000\n");
  fprintf(f2,"0 1 0 1 1 0 0 100000\n");
  fprintf(f2,"@@@\n");
  for (i=0; i < nummol; i++)
    {
      fscanf(f, "%lf %lf %lf\n", &r1[0], &r1[1], &r1[2]);
      fscanf(f, "%lf %lf %lf\n", &r2[0], &r2[1], &r2[2]);
      for (kk=0; kk < 3; kk++)
	{
	  rO[kk] = r1[kk];
	  r21[kk] = r2[kk];
	}
      norm = calc_norm(r21);
      for (kk=0; kk < 3; kk++)
         r21[kk]/=norm;	
      versor_to_R(r21[0], r21[1], r21[2], R);
      //evalfourth(r21,r31,r4,0.9*radius*2.0);
#if 0
      if (i==0)
	fprintf(stderr, "radius21: %f radius31: %f\n", calc_norm(r21), calc_norm(r31));
      for (kk=0; kk < 3; kk++)
	{
	  r4[kk] += rO[kk];
	  //fprintf(stderr,"r1: %f r4: %f\n", r1[kk], r4[kk]);
	}
      for (kk=0; kk < 3; kk++)
	r41[kk] = r4[kk] - r1[kk];

      if (i==0)
	{
	  fprintf(stderr,"radius41=%.15G \n", calc_norm(r41));
	}
#endif
       fprintf(f2,"%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G 0\n",
	     rO[0], rO[1], rO[2], R[0][0], R[0][1], R[0][2], R[1][0], R[1][1], R[1][2],
             R[2][0], R[2][1], R[2][2]);
    }
  fclose(f);
#if 0
  rTemp = sqrt(temp / m); 
   /* assegna le velocita'! */
  sumx = sumy = sumz = 0;
  for (i = 0; i < nummol; i++)
    {
      vx[i] = rTemp * gauss(); 
      vy[i] = rTemp * gauss();
      vz[i] = rTemp * gauss();
      sumx += vx[i]*m;
      sumy += vy[i]*m;
      sumz += vz[i]*m;
      angvel(&wx[i], &wy[i], &wz[i], nummol);
    }
   sumx /= nummol*m;
   sumy /= nummol*m;
   sumz /= nummol*m;
   for (i = 0; i < nummol; i++)
    {
      vx[i] -= sumx;
      vy[i] -= sumy;
      vz[i] -= sumz;
    }	
  K = calc_energy(nummol);
  sf = sqrt( ( (6.0*((double)nummol)-3.0) * temp ) / (2.0*K) );
  for (i = 0; i < nummol; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      wx[i] *= sf;
      wy[i] *= sf;
      wz[i] *= sf;
      /* scala anche i tempi di collisione! */
    } 
 
   for (i = 0; i < nummol; i++)
    {
      fprintf(f2,"%.15G %.15G %.15G %.15G %.15G %.15G\n", vx[i], vy[i], vz[i], 
      	wx[i], wy[i], wz[i]);
    }	
#endif
  fprintf(f2,"%f %f %f\n", bx, by, bz); // lato scatola
  fclose(f2); 
#if 0
  free(vx);
  free(vy);
  free(vz);
  free(wx);
  free(wy);
  free(wz);
#endif
  return 0;
}


