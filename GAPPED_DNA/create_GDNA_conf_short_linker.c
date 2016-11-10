#include<stdlib.h>
#include<stdio.h>
#include<math.h>
double nx, ny, nz, L[3], *rx, *ry, *rz, extradel;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
double *rxCM, *ryCM, *rzCM;
double *Rc[3][3], *Ri[3][3];
int full, ibeg, numpoly;
#define maxpolylen 10000;
double R0[3][3];
#define Sqr(x) ((x)*(x))
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
void versor_to_R(double ox, double oy, double oz, double R[3][3])
{
  int k;
  double angle, u[3], sp, norm, up[3], xx, yy;
#ifdef MC_BENT_DBLCYL
  double Rout[3][3];
  int k1, k2;
#endif
  /* first row vector */
  norm = sqrt(Sqr(ox)+Sqr(oy)+Sqr(oz));
  R[0][0] = ox/norm;
  R[0][1] = oy/norm;
  R[0][2] = oz/norm;
  printf("orient=%f %f %f\n", ox, oy, oz);
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
  /* third row vector */
  vectProdVec(R[0], R[1], u);
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}


double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
void genorient(double *omx, double *omy, double* omz)
{
  int i;
  //double inert;                 /* momentum of inertia of the molecule */
  //double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, osq, norm;
  
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

  *omx = ox;
  *omy = oy;
  *omz = oz; 
}
double min3(double a, double b, double c)
{
  double m;
  m = a;
  if (b < m)
    m = b;
  if (c < m)
    m = c;
  return m;
}


int main(int argc, char **argv)
{
  FILE *f;
  double orient, theta0, theta0rad, Diam, del0, del0x, del0y, del0z, maxL, pi, norm;
  double vol, permdiam, thmax, del, sigb, delfb1, delfb2, delfb3, delfb4, Len, inclth, dd, uu[4][3];
  double del00x, del00y, del00z, *rxCM, *ryCM, *rzCM, bs[3], factor[3], delta;
  double phi, targetphi=0.25, xtrafact, rp1[3], rp2[3], rp3[3], rp4[3];
  int k1, k2, numpoly, parnum=1000, i, j, polylen=4, a, b;
  int nx, ny, nz, nxmax, nymax, nzmax, idx, kk;
  del=0.5;
  /* permanent spots diameter */
  permdiam=2.205; /* 20 T * 0.63 nm dove 0.63 nm è la lunghezza per base stimata per ssDNA in BiophysJ 86, 2630 (2004) */  
  //permdiam=0.5;
  Len=16.0; /* 48 bp equal roughly to 16 nm */
  if (argc == 1)
   {
     printf("create_GDNA_conf <conf_file_name> <nxmax> <nymax> <nzmax> <phi> <diam> <permdiam> \n"); 
     exit(-1);
   }
  f = fopen(argv[1], "w+");
  if (f==NULL)
    {
      printf("boh...f is null!\n");
      exit(-1);
    }
  Diam=2.0;
  pi = acos(0.0)*2.0;

  nxmax = 5;
  if (argc > 2) 
    nxmax = atoi(argv[2]);
  
  nymax = 5;
  if (argc > 3)
    nymax = atoi(argv[3]);

  nzmax = 5;
  if (argc > 4)
    nzmax = atoi(argv[4]);

  if (argc > 5)
    targetphi = atof(argv[5]);

  if (argc > 6)
    Diam = atof(argv[6]);

  if (argc > 7)
    permdiam = atof(argv[7]);

  parnum = 4*nxmax*nymax*nzmax;
  numpoly = nxmax*nymax*nzmax;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxCM = malloc(sizeof(double)*numpoly);
  ryCM = malloc(sizeof(double)*numpoly);
  rzCM = malloc(sizeof(double)*numpoly);
  rxc =malloc(sizeof(double)*polylen);
  ryc =malloc(sizeof(double)*polylen);
  rzc =malloc(sizeof(double)*polylen);
  /* NOTA: calcola l'angolo d'inclinazione dei due cilindri per far si che le patch
   * permamenti siano sovrapposte */
  //dd = sqrt(Sqr((Diam - permdiam)*0.5) + Sqr(permdiam*0.5));  
  inclth = M_PI*0.5-2.0*atan(permdiam/Diam);
  inclth *= 1.1;
  printf("inclination angle=%f deg\n", 180*inclth/M_PI);
  /* ----------------------------- */
  for (a=0; a < 3; a++)
    {
      for (b=0; b < 3; b++)
	{
	  Rc[a][b] = malloc(sizeof(double)*polylen);
	  Ri[a][b] = malloc(sizeof(double)*parnum);
	}
    }
  
#if 0
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Ri[k1][k2] = malloc(sizeof(double)*parnum);     
	R[k1][k2] =  malloc(sizeof(double)*polylen);
      }
#endif
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R0[k1][k2]=0.0;
  R0[0][0]=R0[1][1]=R0[2][2]=1.0;
  delta = 0.1;
  /* building the dimer... */
#if 0
  rxc[0] = Len/2.0+delta;
  ryc[0] = -Diam/2.0-delta;
  rzc[0] = 0.0;
#endif
  //dd = permdiam
  rp1[0] = 0.0;
  rp1[1] = permdiam*0.5-0.0001;
  rp1[2] = 0.0;
  rp2[0] = 0.0;
  rp2[1] = -permdiam*0.5+0.0001;
  rp2[2] = 0.0;

  rxc[0] = rp1[0] - (permdiam*0.5+Len*0.5)*cos(inclth);
  ryc[0] = rp1[1] + (permdiam*0.5+Len*0.5)*sin(inclth);
  rzc[0] = 0.0;
  printf("cos %f sin %f\n", cos(inclth), sin(inclth));
  uu[0][0] = rp1[0] - rxc[0];
  uu[0][1] = rp1[1] - ryc[0];
  uu[0][2] = rp1[2] - rzc[0];

   printf("rcos=%f rsin=%f #0 uu=%f %f %f\n", (permdiam*0.5+Len*0.5)*sin(inclth), (permdiam*0.5+Len*0.5)*cos(inclth), uu[0][0], uu[0][1], uu[0][2]);
   norm = calc_norm(uu[0]);
   for (kk=0; kk < 3; kk++)
     uu[0][kk] /= norm;
   versor_to_R(uu[0][0], uu[0][1], uu[0][2], R0);
#if 1
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][0] = R0[a][b];
      }
#endif

#if 0
  rxc[1] = Len/2.0+delta;
  ryc[1] = Diam/2.0+delta;
  rzc[1] = 0.0;
#endif
  rxc[1] = rp2[0] - (permdiam*0.5+Len*0.5)*cos(inclth);
  ryc[1] = rp2[1] - (permdiam*0.5+Len*0.5)*sin(inclth);
  rzc[1] = 0.0;

  uu[1][0] = rp2[0] - rxc[1];
  uu[1][1] = rp2[1] - ryc[1];
  uu[1][2] = rp2[2] - rzc[1];

  norm = calc_norm(uu[1]);
  for (kk=0; kk < 3; kk++)
    uu[1][kk] /= norm;

  versor_to_R(uu[1][0], uu[1][1], uu[1][2], R0);

  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][1] = R0[a][b];
      }

  dd = Diam*0.5*sin(inclth) + (Len+permdiam*0.5)*cos(inclth);
  /* set origin in the center of mass of the dimer */
  rxc[0] += dd;
  rxc[1] += dd;
  /* ------------------------------------------- */

  rxc[2] = -rxc[0];
  ryc[2] =  ryc[0];
  rzc[2] =  rzc[0];
  for (kk=1; kk < 3; kk++)
    uu[2][kk] = uu[0][kk];
  uu[2][0] = -uu[0][0]; 
  versor_to_R(uu[2][0], uu[2][1], uu[2][2], R0);
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][2] = R0[a][b];
      }

  rxc[3] = -rxc[1];
  ryc[3] =  ryc[1];
  rzc[3] =  rzc[1];

  for (kk=1; kk < 3; kk++)
    uu[3][kk] = uu[1][kk];
  uu[3][0] = -uu[1][0];
  versor_to_R(uu[3][0], uu[3][1], uu[3][2], R0);
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][3] = R0[a][b];
      }
  bs[0] = 2.0*dd + permdiam + 2.0*delta;
  bs[1] = 2.0*Len*sin(inclth)+2.0*Diam*cos(inclth)+2.0*delta;
  bs[2] = Diam;
  printf("bs=%f %f %f\n", bs[0], bs[1], bs[2]);
	printf("A=%f B=%f\n",2.0*Len*sin(inclth),2.0*Diam*cos(inclth));
#if 0
  for (i=0; i < polylen; i++)
    {
      rxc[i] = 0.0;
      ryc[i] = i*1.0001*Diam;
      rzc[i] = 0.0;
#if 0
      /*apply theta0 rotation around x-axis */
      theta0rad=pi*(i*theta0)/180.0;
      R0[1][1]=cos(theta0rad);
      R0[1][2]=-sin(theta0rad);
      R0[2][1]=sin(theta0rad);
      R0[2][2]=cos(theta0rad);
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  R[k1][k2][i] = R0[k1][k2];
#endif
    }
#endif
  fprintf(f, "parnum: %d\n", parnum);
  fprintf(f,"ninters: 2\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"%f %f %f\n", Len/2.0, Diam/2.0, Diam/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"%f 0 0 %f\n", permdiam*0.5+Len/2.0, permdiam);/* 0: along x axis (permanent) 0.05 means lp=20 */
  fprintf(f,"%f 0 0 %f\n", -Len/2.0-0.15, 0.5);
  fprintf(f,"0 0 0 0 1 0 0 1\n");
  fprintf(f,"0 1 0 1 1 0 0 100000\n");
  fprintf(f, "@@@\n");
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  //maxL = L;
  //numpoly=parnum/polylen;
  //printf("numpoly=%d numpoly*polylen=%d\n", numpoly, numpoly*polylen);
#if 0
  del00x = Len*0.5;
  del00y = Diam;
  del00z = Diam;
  del0x = 0.01;
  del0y=del0z=Diam*0.1;
#endif
  drx = bs[2]; //0.500000000001;
  dry = bs[2];
  drz = bs[2];
#if 0
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
#endif
  /* place on a cubic lattice */ 
  for (nx=0; nx < nxmax; nx++)
    for (ny=0; ny < nymax; ny++)
      for (nz=0; nz < nzmax; nz++)
	{
	  idx = nz+nzmax*ny+nzmax*nymax*nx;
      	  rxCM[idx] = (nx+0.5001)*drx;
	  ryCM[idx] = (ny+0.5001)*dry;
	  rzCM[idx] = (nz+0.5001)*drz;
       	}
  L[0] = (nxmax)*drx;
  L[1] = (nymax)*dry;
  L[2] = (nzmax)*drz;

  printf("step #1 L=%f %f %f\n", L[0], L[1], L[2]);
  /* expand system according to tetramer size bs[3] */
  factor[0] = 1.0001*(bs[0]/bs[2]);
  factor[1] = 1.0001*(bs[1]/bs[2]);
  factor[2] = 1.0001;
  for (i=0; i < numpoly; i++)
    {
      rxCM[i] *= factor[0];
      ryCM[i] *= factor[1];
      rzCM[i] *= factor[2];
    } 
  for (a=0; a < 3; a++)
    L[a] *= factor[a];
  printf("step #2 [factor] L=%f %f %f\n", L[0], L[1], L[2]);
  vol = acos(0.0)*2.0*(Diam/2)*(Diam/2)*Len;
  phi=parnum*vol/(L[0]*L[1]*L[2]);
  xtrafact = pow(phi/targetphi, 1.0/3.0);
  printf("vol=%f targetphi=%f phi=%f xtrafact=%f\n", vol, targetphi, phi, xtrafact);

  if (xtrafact < 1.0)
    xtrafact = 1.0001;
  for (i=0; i < numpoly; i++)
    {
      rxCM[i] *= xtrafact;
      ryCM[i] *= xtrafact;
      rzCM[i] *= xtrafact;
    } 
  for (a=0; a < 3; a++)
    L[a] *= xtrafact;
 
  printf("parnum=%d nx=%d ny=%d nz=%d argc=%d L=%f %f %f\n", parnum, nxmax, nymax, nzmax, argc, L[0], L[1], L[2]);
#if 0
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx;
      if (rxl+Len*0.5 > L)
	{
	  ny++;
	  nx=0;
	  rxl = rxc[polylen-1]+nx*drx;
	}
      ryl = ryc[polylen-1]+ny*dry;
      if (ryl+Diam*0.5 > L)
	{
	  nz++;
	  ny=0;
	  ryl = ryc[polylen-1]+ny*dry;
      	}
      rzl = rzc[polylen-1]+nz*drz;
      if (rzl+Diam*0.5 > L)
	{
	  full=1;
	}

      for (j = 0; j < polylen; j++)
	{
	  rx[i*polylen+j] = rxc[j]+nx*drx;
	  //printf("nx=%f i=%d j=%d x=%f\n", nx, i, j,  rx[i*polylen+j]);
	  ry[i*polylen+j] = ryc[j]+ny*dry;
	  rz[i*polylen+j] = rzc[j]+nz*drz+del00z;
#if 0
      	  for (k1=0; k1 < 3; k1++)
    	    for (k2=0; k2 < 3; k2++)
    	      Ri[k1][k2][i*polylen+j]= R[k1][k2][j];
#endif
	  //printf("np=%d\n", i*polylen+j);
	}
      if (full==1)
	{
	  printf("I could only place %d particles!\n", i);
	  break;
	}
      nx++;
    }
#endif
  //  if (full)
  //   exit(-1);
#if 1
  for (i=0; i < numpoly; i++)
    {
      for (j=0; j < polylen; j++)
	{
	  idx = i*polylen + j;
	  rx[idx] = rxCM[i] + rxc[j];
	  ry[idx] = ryCM[i] + ryc[j];
	  rz[idx] = rzCM[i] + rzc[j];
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      Ri[a][b][idx] = Rc[a][b][j];
	  //printf("qui idx=%d\n", idx);
	}
    }
  for (i=0; i < parnum; i++)
    {
      rx[i] -= L[0]*0.5;
      ry[i] -= L[1]*0.5;
      rz[i] -= L[2]*0.5;
      //printf("qui i=%d\n", i);
      //orient=(i%2==0)?1.0:-1.0;
      fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f 0\n", rx[i], ry[i], rz[i], Ri[0][0][i], Ri[0][1][i], Ri[0][2][i], Ri[1][0][i], Ri[1][1][i], Ri[1][2][i], Ri[2][0][i], Ri[2][1][i], Ri[2][2][i]);
      //fprintf(f, "%f %f %f  0\n", rx[i], ry[i], rz[i]);
      //printf("qui2\n");
    }
#endif	
  fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  printf("phi=%f\n", parnum*vol/(L[0]*L[1]*L[2]));
  fclose(f);
  return 1;
} 
