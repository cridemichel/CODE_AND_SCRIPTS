#include<stdlib.h>
#include<stdio.h>
#include<math.h>
double nx, ny, nz, L[3], *rx, *ry, *rz, extradel;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
double *rxCM, *ryCM, *rzCM;
double *Rc[3][3], *Ri[3][3];
int full, ibeg, numpoly;
#define maxpolylen 10000;
double R0[3][3], R1[3][3];
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
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

double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
#ifdef MULTIARM
/* a seconda del numero di braccia le patch vengono disposte in una maniera regolare opportuna,
 * ossia triangolo (n=3), tetraedro (n=4), triangolo + 2 patch sull'asse ad esso perpendicolare (n=5).
 * Notare che sono versori unitari. */
double patchGeomN3[3][3] = {{0,0,-1},{0.866025,0,0.5},{-0.866025,0,0.5}};
double patchGeomN4[4][3] = {{0.816497, -1.1547, -0.333333},{-0.816497, -1.1547, -0.333333},{0., 2.3094, -0.333333},{1.,0,0}};
double patchGeomN5[5][3] = {{0,0,-1},{0.866025,0.0,0.5},{-0.866025,0.0,0.5},{0,-1,0},{0,1,0}};
double patchGeom[5][3];
#endif
#ifdef MULTIARM
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
  R[2][0] = ox/norm;
  R[2][1] = oy/norm;
  R[2][2] = oz/norm;
  //printf("orient=%f %f %f\n", ox, oy, oz);
  u[0] = 0.0; u[1] = 1.0; u[2] = 0.0;
  if (u[0]==R[2][0] && u[1]==R[2][1] && u[2]==R[2][2])
    {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;
    }
  /* second row vector */
  sp = 0;
  for (k=0; k < 3 ; k++)
    sp+=u[k]*R[2][k];
  for (k=0; k < 3 ; k++)
    u[k] -= sp*R[2][k];
  norm = calc_norm(u);
  //printf("norm=%f u=%f %f %f\n", norm, u[0], u[1], u[2]);
  for (k=0; k < 3 ; k++)
    R[1][k] = u[k]/norm;
  /* third row vector */
  vectProdVec(R[1], R[2], u);
 
  for (k=0; k < 3 ; k++)
    R[0][k] = u[k];
}
#endif
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}

double calcDistBox(double rbi[3], double rbj[3], double saxi[3], double saxj[3], double Ri[3][3], double Rj[3][3])
{
  double RR, R0, R1, cij[3][3], fabscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k, k1, k2, existsParallelPair = 0;
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  for (k=0; k < 3; k++)
    {
      rA[k] = rbi[k];
      rB[k] = rbj[k];
      EA[k] = saxi[k];
      EB[k] = saxj[k];
    }
#if 0
  /* riportare qua anche l'analogo routin per sfere se servirà */
  if (is_a_sphere_NNL[i] && is_a_sphere_NNL[j])
    return calcDistNegNNLoverlapPlaneHS(i, j, rA, rB);
#endif
  
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  AA[k1][k2] = Ri[k1][k2];
	  BB[k1][k2] = Rj[k1][k2];
	}
    	DD[k1] = rA[k1] - rB[k1];
    }
  /* axis C0+s*A0 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[0][k1] =  scalProd(AA[0], BB[k1]);
      fabscij[0][k1] = fabs(cij[0][k1]);
      if ( fabscij[0][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[0] = scalProd(AA[0],DD);
  RR = fabs(AD[0]);
  R1 = EB[0]*fabscij[0][0]+EB[1]*fabscij[0][1]+EB[2]*fabscij[0][2];
  R01 = EA[0] + R1;
  if ( RR > R01 )
    return 1.0; /* non si intersecano */
  /* axis C0+s*A1 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[1][k1] = scalProd(AA[1],BB[k1]);
      fabscij[1][k1] = fabs(cij[1][k1]);
      if ( fabscij[1][k1] == 1.0  )
	existsParallelPair = 1;
    }
  AD[1] = scalProd(AA[1],DD);
  RR = fabs(AD[1]);
  R1 = EB[0]*fabscij[1][0]+EB[1]*fabscij[1][1]+EB[2]*fabscij[1][2];
  R01 = EA[1] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*A2 */
  for (k1= 0; k1 < 3; k1++)
    {
      cij[2][k1] = scalProd(AA[2], BB[k1]);
      fabscij[2][k1] = fabs(cij[2][k1]);
      if ( fabscij[2][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[2] = scalProd(AA[2],DD);
  RR = fabs(AD[2]);
  R1 = EB[0]*fabscij[2][0]+EB[1]*fabscij[2][1]+EB[2]*fabscij[2][2];
  R01 = EA[2] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*B0 */
  RR = fabs(scalProd(BB[0],DD));
  R0 = EA[0]*fabscij[0][0]+EA[1]*fabscij[1][0]+EA[2]*fabscij[2][0];
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B1 */
  RR = fabs(scalProd(BB[1],DD));
  R0 = EA[0]*fabscij[0][1]+EA[1]*fabscij[1][1]+EA[2]*fabscij[2][1];
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;
  
  /* axis C0+s*B2 */
  RR = fabs(scalProd(BB[2],DD));
  R0 = EA[0]*fabscij[0][2]+EA[1]*fabscij[1][2]+EA[2]*fabscij[2][2];
  R01 = R0 + EB[2];
  if ( RR > R01 )
    return 1.0;

  /* At least one pair of box axes was parallel, therefore the separation is
   * effectively in 2D, i.e. checking the "edge" normals is sufficient for
   * the separation of the boxes. 
   */
  if ( existsParallelPair )
    return -1.0;

  /* axis C0+s*A0xB0 */
  RR = fabs(AD[2]*cij[1][0]-AD[1]*cij[2][0]);
  R0 = EA[1]*fabscij[2][0] + EA[2]*fabscij[1][0];
  R1 = EB[1]*fabscij[0][2] + EB[2]*fabscij[0][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB1 */
  RR = fabs(AD[2]*cij[1][1]-AD[1]*cij[2][1]);
  R0 = EA[1]*fabscij[2][1] + EA[2]*fabscij[1][1];
  R1 = EB[0]*fabscij[0][2] + EB[2]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB2 */
  RR = fabs(AD[2]*cij[1][2]-AD[1]*cij[2][2]);
  R0 = EA[1]*fabscij[2][2] + EA[2]*fabscij[1][2];
  R1 = EB[0]*fabscij[0][1] + EB[1]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB0 */
  RR = fabs(AD[0]*cij[2][0]-AD[2]*cij[0][0]);
  R0 = EA[0]*fabscij[2][0] + EA[2]*fabscij[0][0];
  R1 = EB[1]*fabscij[1][2] + EB[2]*fabscij[1][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB1 */
  RR = fabs(AD[0]*cij[2][1]-AD[2]*cij[0][1]);
  R0 = EA[0]*fabscij[2][1] + EA[2]*fabscij[0][1];
  R1 = EB[0]*fabscij[1][2] + EB[2]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB2 */
  RR = fabs(AD[0]*cij[2][2]-AD[2]*cij[0][2]);
  R0 = EA[0]*fabscij[2][2] + EA[2]*fabscij[0][2];
  R1 = EB[0]*fabscij[1][1] + EB[1]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB0 */
  RR = fabs(AD[1]*cij[0][0]-AD[0]*cij[1][0]);
  R0 = EA[0]*fabscij[1][0] + EA[1]*fabscij[0][0];
  R1 = EB[1]*fabscij[2][2] + EB[2]*fabscij[2][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB1 */
  RR = fabs(AD[1]*cij[0][1]-AD[0]*cij[1][1]);
  R0 = EA[0]*fabscij[1][1] + EA[1]*fabscij[0][1];
  R1 = EB[0]*fabscij[2][2] + EB[2]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB2 */
  RR = fabs(AD[1]*cij[0][2]-AD[0]*cij[1][2]);
  R0 = EA[0]*fabscij[1][2] + EA[1]*fabscij[0][2];
  R1 = EB[0]*fabscij[2][1] + EB[1]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  return -1.0;
}

void orient(double *omx, double *omy, double* omz)
{
  int i;
  //double inert;                 /* momentum of inertia of the molecule */
  //double norm, dot, osq, o, mean;
  double  xisq, xi1, xi2, xi;
  double ox, oy, oz, osq, norm;
  
  xisq = 1.0;

  while (xisq >= 1.0)
    {
      xi1  = ranf_vb() * 2.0 - 1.0;
      xi2  = ranf_vb() * 2.0 - 1.0;
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

int main(int argc, char **argv)
{
  FILE *f;
  double ox, oy, oz;
  double sigChain, DiamSph, orient, theta0, theta0rad, Diam, del0, del0x, del0y, del0z, maxL, pi;
  double vol, permdiam, thmax, del, sigb, delfb1, delfb2, delfb3, delfb4, Len, sigmaAntigens, sigSph;
  double del00x, del00y, del00z, *rxCM, *ryCM, *rzCM, bs[3], factor[3], deltax, deltay, deltaz, DiamAntigen;
  double Rcmx, Rcmy, Rcmz, rxa, rya, rza, dist, dx, dy, dz, dxMax, dyMax, dzMax, distRevPatch, bigAntigenSurfDiam=0.0;
  double phi, targetphi=0.25, xtrafact, Lx=10.0, Ly=10.0, Lz=10.0, nanorevpatchDiam, rsp1[3], rsp2[3];
  double rA[3], rB[3], saA[3], saB[3];
  int numpolyignored=0;
  int k1, k2, numpoly, parnum=1000, i, j, polylen, a, b, numSpheres=3, idx1, idx2;
  int type, kk, k, overlap, nx, ny, nz, nxmax, nymax, nzmax, idx, numantigens, *ignorepoly, numpolyeff;
#ifdef MULTIARM
  int numarms=2, *typesnb;
#endif
  del=0.5;
  /* nanobody (Fab) permanent spots diameter */
  permdiam=0.8; 

  /* diametro antigeni */
  DiamAntigen = 0.79;
  
  //permdiam=0.5;
  if (argc == 1)
    {
#ifdef MULTIARM
      printf("create_NANOBODY_conf_nematic_unfolded <conf_file_name> <nxmax> <nymax> <nzmax> <Lx> <Ly> <Lz> <DensSuperfAntigens> <numspheres-per-arm> <QFab-diam> <QFab-len> <QFab-diam-permpatch> <QFab-diam-revpatch> <QFab-dist-revpatch> <sphere-diam> <sphere-revpatch-diam> <DiametroAntigene> <bigAntigenSurfDiam> <numarms>\n"); 
#else
     printf("create_NANOBODY_conf_nematic_unfolded <conf_file_name> <nxmax> <nymax> <nzmax> <Lx> <Ly> <Lz> <DensSuperfAntigens> <numspheres> <QFab-diam> <QFab-len> <QFab-diam-permpatch> <QFab-diam-revpatch> <QFab-dist-revpatch> <sphere-diam> <sphere-revpatch-diam> <DiametroAntigene> <bigAntigenSurfDiam>\n"); 
#endif
     exit(-1);
   }
  f = fopen(argv[1], "w+");
  if (f==NULL)
    {
      printf("boh...f is null!\n");
      exit(-1);
    }
  /* diametro e lunghezza nanobody Fab (ellissoide prolato) */
  Diam=2.0;
  Len = 4.0;
  nanorevpatchDiam = 0.45; 
  /* diametro della sfera della catena */
  DiamSph=1.0; 

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
    Lx = atof(argv[5]);

  if (argc > 6)
    Ly = atof(argv[6]);

  if (argc > 7)
    Lz = atof(argv[7]);

  if (argc > 8)
    sigmaAntigens = atof(argv[8]);

  if (argc > 9)
    numSpheres = atoi(argv[9]);

  if (argc > 10 && atof(argv[10])!=-1)
    Diam = atof(argv[10]);

  if (argc > 11 && atof(argv[11])!=-1)
    Len = atof(argv[11]);

  if (argc > 12 && atof(argv[12])!=-1)
    permdiam = atof(argv[12]);
  
  if (argc > 13 && atof(argv[13])!=-1)
    nanorevpatchDiam = atof(argv[13]);

  if (argc > 14 && atof(argv[14])!=-1)
    distRevPatch = atof(argv[14]);
  else
    distRevPatch = Len/2.0; 

  if (argc > 15 && atof(argv[15])!=-1)
    DiamSph = atof(argv[15]);

  if (argc > 16 && atof(argv[16])!=-1)
    sigSph = atof(argv[16]);
  else
    sigSph = 0.119*DiamSph;

  if (argc > 17 && atof(argv[17])!=-1)
    DiamAntigen = atof(argv[17]);

  /* se bigAntigenSurfDiam = 0 allora mette gli antigeni sulla faccia in basso */
  if (argc > 18 && atof(argv[18])!=-1)
    {
      bigAntigenSurfDiam = atof(argv[19]);
    }
  else
    bigAntigenSurfDiam = 120.0;
   
#ifdef MULTIARM
  if (argc > 19)
    numarms = atoi(argv[19]);

  if (numarms > 5)
    {
      printf("ERROR: maximum supported branches is 5!\n");
      exit(-1);	
    }
  if (numarms < 3)
    {
      printf("Compile without -DMULTI_ARM flag to generate bi-bodies\n");
      exit(-1);
    }
#endif

#ifdef MULTIARM
  /* ogni braccio ha numSphere sfere ed in più c'è la branching sphere */
  polylen = numarms*numSpheres + 1 + numarms;
  parnum = polylen*nxmax*nymax*nzmax;
  typesnb = malloc(sizeof(int)*polylen);
#else
  polylen = numSpheres + 2;
  parnum = polylen*nxmax*nymax*nzmax;
#endif
  numpoly = nxmax*nymax*nzmax;
  ignorepoly = malloc(sizeof(int)*numpoly);
  for (i=0; i < numpoly; i++)
    ignorepoly[i] = 0;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxCM = malloc(sizeof(double)*numpoly);
  ryCM = malloc(sizeof(double)*numpoly);
  rzCM = malloc(sizeof(double)*numpoly);
  rxc =malloc(sizeof(double)*polylen);
  ryc =malloc(sizeof(double)*polylen);
  rzc =malloc(sizeof(double)*polylen);
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
  deltax = 0.1;
  deltay = 0.1;
  deltaz = 0.1;
  /* building the nanobody... */
#ifdef MULTIARM
  if (numarms==3)
    {
      for (k1=0; k1 < 3; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      patchGeom[k1][k2] = patchGeomN3[k1][k2];
	    }
	}
    }
  else if (numarms == 4)
    {
      for (k1=0; k1 < 4; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      patchGeom[k1][k2] = patchGeomN4[k1][k2];
	    }
	}
    }
  else if (numarms == 5)
    {
      for (k1=0; k1 < 5; k1++)
	{
	  for (k2=0; k2 < 3; k2++)
	    {
	      patchGeom[k1][k2] = patchGeomN5[k1][k2];
	    }
	}
    }

  /* first place the branching sphere... */
  rxc[0] = 0.0;
  ryc[0] = 0.0;
  rzc[0] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][0] = R0[a][b];
      }
  typesnb[0] = 2;
  idx = 0;
  for (k1=0; k1 < numarms; k1++)
    {
      for (k2 = 0; k2 < numSpheres; k2++)
	{
	  dist = (DiamSph+deltaz)*(k2+1);
       	  idx++;
	  rxc[idx] = dist*patchGeom[k1][0];
	  ryc[idx] = dist*patchGeom[k1][1];
	  rzc[idx] = dist*patchGeom[k1][2];
	  typesnb[idx] = 1;
	  versor_to_R(-patchGeom[k1][0], -patchGeom[k1][1], -patchGeom[k1][2], R1);
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
		Rc[a][b][idx] = R1[a][b];
	      }
 	}
      /* place nanobody at the end of the arm */
      idx++;
      dist = (DiamSph+deltaz)*numSpheres + (DiamSph+deltaz)*0.5 + (Len+deltaz)*0.5;
      typesnb[idx] = 0;
      rxc[idx] = dist*patchGeom[k1][0];
      ryc[idx] = dist*patchGeom[k1][1];
      rzc[idx] = dist*patchGeom[k1][2];
      /* il meno ci vuole poiché solo una della due patch del nanobody è irreversibile */
      versor_to_R(-patchGeom[k1][0], -patchGeom[k1][1], -patchGeom[k1][2], R1);
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  {
	    Rc[a][b][idx] = R1[a][b];
	  }
    }
#else
   /* ellipsoid */	
  rxc[0] = 0.0;
  ryc[0] = 0.0;
  rzc[0] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][0] = R0[a][b];
      }
 /* spheres */
  for (k=0; k < numSpheres; k++)
    {
      rxc[k+1] = 0.0;
      ryc[k+1] = 0.0;
      rzc[k+1] = (Len+deltaz)*0.5+(DiamSph+deltaz)*0.5 + (DiamSph+deltaz)*k;
      for (a=0; a < 3; a++)
	for (b=0; b < 3; b++)
	  {
	    Rc[a][b][k+1] = R0[a][b];
	  }
    }
  /* ellipsoid */
  rxc[polylen-1] = 0.0;
  ryc[polylen-1] = 0.0;
  rzc[polylen-1] = (Len+deltaz)+(DiamSph+deltaz)*numSpheres;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][polylen-1] = R0[a][b];
      }
#endif
  /* center the nanobody */
#ifdef MULTIARM
  Rcmx=Rcmy=Rcmz=0.0;
  for (k=0; k < polylen; k++)
    {
      Rcmx += rxc[k];
      Rcmy += ryc[k];
      Rcmz += rzc[k];
    }
  Rcmx /= ((double) polylen);
  Rcmy /= ((double) polylen);
  Rcmz /= ((double) polylen);
  for (k=0; k < polylen; k++)
    {
      rxc[k]-=Rcmx;
      ryc[k]-=Rcmy;
      rzc[k]-=Rcmz;
    }
#else
  Rcmz=0.0;
  for (k=0; k < polylen; k++)
    {
      Rcmz += rzc[k];
    }
  Rcmz /= (2.0+numSpheres);
  for (k=0; k < polylen; k++)
    {
      rzc[k]-=Rcmz;
    }
#endif
#ifdef MULTIARM
  /* calcola il box che racchiude il nanobody usando solo i Fab,
   * nel caso numarms=3 il box lungo y deve avere larghezza pari al diamtro 
   * delle o se più grande a quello del semiasse minore degli ellissoidi */
  for (k1 = 0; k1 < numarms; k1++)
    {
      for (k2 = k1+1; k2 < numarms; k2++)
	{
	  idx1 = (numSpheres+1)*(k1+1); 
	  idx2 = (numSpheres+1)*(k2+1);
	  printf("type[%d]=%d type[%d]=%d\n", idx1, typesnb[idx1], idx2, typesnb[idx2]);
      	  rsp1[0] = rxc[idx1] - distRevPatch*Rc[2][0][idx1];
	  rsp2[0] = rxc[idx2] - distRevPatch*Rc[2][0][idx2];
	  rsp1[1] = ryc[idx1] - distRevPatch*Rc[2][1][idx1];
	  rsp2[1] = ryc[idx2] - distRevPatch*Rc[2][1][idx2];
	  rsp1[2] = rzc[idx1] - distRevPatch*Rc[2][2][idx1];
	  rsp2[2] = rzc[idx2] - distRevPatch*Rc[2][2][idx2];
	 
	  //printf("orient=%f %f %f  - %f %f %f\n", Rc[2][0][idx1], Rc[2][1][idx1], Rc[2][2][idx1],
	    //  Rc[2][0][idx2], Rc[2][1][idx2], Rc[2][2][idx2]);
	  dx = rsp1[0]-rsp2[0];
	  dy = rsp1[1]-rsp2[1];
	  dz = rsp1[2]-rsp2[2];
	  printf("dy=%f\n", dy);
	  if (k1==0 && k2 ==1)
	    {
	      dxMax = fabs(dx);
	      dyMax = fabs(dy);
	      dzMax = fabs(dz);
	    }
	  if (fabs(dx) > dxMax) 
	    dxMax = fabs(dx);
	  if (fabs(dy) > dyMax) 
	    dyMax = fabs(dy);
	  if (fabs(dz) > dzMax) 
	    dzMax = fabs(dz);
	}
    } 
  //dxMax += Len;
  //dyMax += Len;
  //dzMax += Len;

  bs[0] = dxMax;
  if (numarms==3)
    bs[1] = (DiamSph > Diam)?DiamSph:Diam;
  else
    bs[1] = dyMax;
 
  bs[2] = dzMax;
  printf("calculated bs=%f %f %f\n", bs[0], bs[1], bs[2]);
#else
  bs[0] = Diam>DiamSph?Diam:DiamSph;
  bs[1] = Diam>DiamSph?Diam:DiamSph; 
  bs[2] = 2.0*(Len+deltaz)+numSpheres*(DiamSph+deltaz);
#endif
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
  drx = bs[0]; //0.500000000001;
  dry = bs[0];
  drz = bs[0];
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
  factor[0] = 1.0001;
  factor[1] = 1.0001;
  factor[2] = 1.0001*(bs[2]/bs[0]);
  for (i=0; i < numpoly; i++)
    {
      rxCM[i] *= factor[0];
      ryCM[i] *= factor[1];
      rzCM[i] *= factor[2];
    } 
  for (a=0; a < 3; a++)
    L[a] *= factor[a];
  printf("step #2 [factor] L=%f %f %f\n", L[0], L[1], L[2]);
  //vol = acos(0.0)*2.0*(Diam/2)*(Diam/2)*Len;
  //phi=parnum*vol/(L[0]*L[1]*L[2]);
  //xtrafact = pow(phi/targetphi, 1.0/3.0);

  //printf("vol=%f targetphi=%f xtrafact=%f\n", vol, targetphi);

  for (i=0; i < numpoly; i++)
    {
      rxCM[i] *= Lx/L[0];
      ryCM[i] *= Ly/L[1];
      rzCM[i] *= Lz/L[2];
    } 
  L[0] = Lx;
  L[1] = Ly;
  L[2] = Lz;
  /* gli antigeni vengono messi nel piano xy (cioè z = -L[2]*0.5) */	
  numantigens = ((int)(sigmaAntigens*(L[0]*L[1]))); 
  printf("parnum=%d nx=%d ny=%d nz=%d argc=%d L=%f %f %f\n", parnum, nxmax, nymax, nzmax, argc, L[0], L[1], L[2]);
#if 1
  for (i=0; i < numpoly; i++)
    {
      if (bigAntigenSurfDiam > 0.0)
	{
	  /* roughly check whether the nanobody is inside
	   * the big sphere */
	  rA[0] = 0;
	  rA[1] = 0;
	  rA[2] = 0;
	  saA[0] = bigAntigenSurfDiam*0.5;
	  saA[1] = bigAntigenSurfDiam*0.5;
	  saA[2] = bigAntigenSurfDiam*0.5;

	  rB[0] = rxCM[i];
	  rB[1] = ryCM[i];
	  rB[2] = rzCM[i];
	  saB[0] = bs[0]*0.5;
	  saB[1] = bs[1]*0.5;
	  saB[2] = bs[2]*0.5;

      	  if (calcDistBox(rA, rB, saA, saB, R0, R0) < 1.0)
	    {
	      ignorepoly[i] = 1;
	      numpolyignored++;
	    }
	  if (ignorepoly[i])
	    continue;	
	}

      for (j=0; j < polylen; j++)
	{
	  idx = i*polylen + j;
	  rx[idx] = rxCM[i] + rxc[j];
	  ry[idx] = ryCM[i] + ryc[j];
	  rz[idx] = rzCM[i] + rzc[j];
  	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      Ri[a][b][idx] = Rc[a][b][j];
	  //printf("qui idx=%d R[0][0]=%f\n", idx, Ri[0][0][idx]);
	}
    }

  /* HEADER */
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "parnum: %d\n", parnum - numpolyignored*polylen + numantigens + 1);
  else
    fprintf(f, "parnum: %d\n", parnum - numpolyignored*polylen + numantigens);
#ifdef MULTIARM
  fprintf(f,"ninters: %d\n", 6+numarms*2);
#else
  fprintf(f,"ninters: 9\n");
#endif
  fprintf(f,"nintersIJ: 0\n");
#ifdef MULTIARM
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "ntypes: 5\n");
  else
    fprintf(f, "ntypes: 4\n");
#else
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "ntypes: 5\n");
  else
    fprintf(f, "ntypes: 4\n");
#endif
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
#ifdef MULTIARM
  numpolyeff = numpoly - numpolyignored;
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "%d %d %d %d 1\n", numpolyeff*numarms, numarms*numSpheres*numpolyeff, numpolyeff, numantigens);
  else
    fprintf(f, "%d %d %d %d\n", numpolyeff*numarms, numarms*numSpheres*numpolyeff, numpolyeff, numantigens);
#else
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "%d %d %d %d 1\n", numpolyeff, numpolyeff, numSpheres*numpolyeff, numantigens);
  else
    fprintf(f, "%d %d %d %d\n", numpolyeff, numpolyeff, numSpheres*numpolyeff, numantigens);
#endif
  /* first Fab */
  fprintf(f,"%f %f %f\n", Diam/2.0, Diam/2.0, Len/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n", Len/2.0, permdiam);/* 0: along z axis (permanent) 0 */
  fprintf(f,"0 0 %f %f\n", -distRevPatch, nanorevpatchDiam); /* 1: along z axis patch which will form bonds with antigens */

#ifndef MULTIARM
  /* second Fab */
  fprintf(f,"%f %f %f\n", Diam/2.0, Diam/2.0, Len/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n", -Len/2.0, permdiam);/* 0: along z axis (permanent) 0 */
  fprintf(f,"0 0 %f %f\n", (Len/2.0-0.5), 0.612); /* 1: along z axis patch which will form bonds with antigens */
#endif

  /* bi-sphere */
  fprintf(f,"%f %f %f\n", DiamSph/2.0, DiamSph/2.0, DiamSph/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n", DiamSph/2.0, sigSph);/* 0: along z axis (permanent) 0 */
  fprintf(f,"0 0 %f %f\n", -DiamSph/2.0, sigSph); /* 1: along z axis patch which will form bonds with antigens */

#ifdef MULTIARM
  /* branching sphere */
  fprintf(f,"%f %f %f\n", DiamSph/2.0, DiamSph/2.0, DiamSph/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"%d 0\n", numarms);
  for (k1=0; k1 < numarms; k1++)
    {
      fprintf(f,"%f %f %f %f\n", patchGeom[k1][0]*(DiamSph/2.0), patchGeom[k1][1]*(DiamSph/2.0),
	      patchGeom[k1][2]*(DiamSph/2.0), sigSph);/* 0: along z axis (permanent) 0 */
    }
#endif

  /* antigen */
  fprintf(f,"%f %f %f\n", DiamAntigen/2.0, DiamAntigen/2.0, DiamAntigen/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1e+200 1 1 1 0 1\n");
  fprintf(f,"1 0\n");
  fprintf(f,"0 0 0 %f\n", DiamAntigen);

  /* big sphere */
  if (bigAntigenSurfDiam > 0.0)
    {
      fprintf(f,"%f %f %f\n", bigAntigenSurfDiam/2.0, bigAntigenSurfDiam/2.0, bigAntigenSurfDiam/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
      fprintf(f,"1 1 1\n");
      fprintf(f, "1e+200 1 1 1 0 0\n");
      fprintf(f,"0 0\n");
    }
  /* interactions */
#ifdef MULTIARM
  fprintf(f,"0 0 1 0 1 1000000 1000000 1\n");
  fprintf(f,"0 0 1 1 1 1000000 1000000 1\n");
  fprintf(f,"1 0 1 0 1 1000000 1000000 1\n");
  fprintf(f,"1 0 1 1 1 1000000 1000000 1\n");
  fprintf(f,"1 1 1 1 1 1000000 1000000 1\n");
  for (k1 = 0; k1 < numarms; k1++)
    {
      fprintf(f, "1 0 2 %d 1 1000000 1000000 1\n", k1);
      fprintf(f, "1 1 2 %d 1 1000000 1000000 1\n", k1);
    }
  fprintf(f,"0 1 3 0 1 0 0 1\n");
#else
  fprintf(f,"0 0 2 0 1 1000000 1000000 1\n");
  fprintf(f,"1 0 2 0 1 1000000 1000000 1\n");
  fprintf(f,"0 0 2 1 1 1000000 1000000 1\n");
  fprintf(f,"1 0 2 1 1 1000000 1000000 1\n");
  fprintf(f,"2 0 2 0 1 1000000 1000000 1\n");
  fprintf(f,"2 0 2 1 1 1000000 1000000 1\n");
  fprintf(f,"2 1 2 1 1 1000000 1000000 1\n");
  fprintf(f,"0 1 3 0 1 0 0 1\n");
  fprintf(f,"1 1 3 0 1 0 0 1\n");
#endif
  fprintf(f, "@@@\n");

  /* place nanobodies */
  for (i=0; i < parnum; i++)
    {
      rx[i] -= L[0]*0.5;
      ry[i] -= L[1]*0.5;
      rz[i] -= L[2]*0.5;
      //printf("qui i=%d\n", i);
      //orient=(i%2==0)?1.0:-1.0;
      k = parnum / polylen;
      if (ignorepoly[k])
	continue;
#ifdef MULTIARM
      kk = i % polylen;
      type = typesnb[kk];
#else
      kk = i % polylen;
      if (kk==0) 
	type = 0;
      else if (kk==polylen-1) 
	type = 1;
      else
	type = 2;
#endif
      fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %d\n", rx[i], ry[i], rz[i], Ri[0][0][i], Ri[0][1][i], Ri[0][2][i], Ri[1][0][i], Ri[1][1][i], Ri[1][2][i], Ri[2][0][i], Ri[2][1][i], Ri[2][2][i], type);
      //fprintf(f, "%f %f %f  0\n", rx[i], ry[i], rz[i]);
      //printf("qui2\n");
    }
#endif	
  /* place antigens */
  rx = realloc(rx,sizeof(double)*(parnum+numantigens));
  ry = realloc(ry,sizeof(double)*(parnum+numantigens));
  rz = realloc(rz,sizeof(double)*(parnum+numantigens));
  if (bigAntigenSurfDiam > 0.0)
    {
      /* setup a non flat geometry here...to start a sphere? */
      for (i=0; i < numantigens; i++)
	{
	  overlap=1;
	  while (overlap)
	    {
	      orient(&ox, &oy, &oz);
	      ox *= bigAntigenSurfDiam*0.5;
	      oy *= bigAntigenSurfDiam*0.5;
	      oz *= bigAntigenSurfDiam*0.5;

	      /* check overlaps between antigens */
	      overlap = 0;
	      for (k = 0; k < i-1; k++)
		{
		  if (Sqr(rx[parnum+i]-rxa) + Sqr(ry[parnum+i]-rya) < Sqr(DiamAntigen))
		    {
		      overlap=1;
		      break;
		    }
		}	  
	    }
	  fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 3\n", ox, oy, oz); 
	}
    }
  else
    {
      /* setup a non flat geometry here...to start a sphere? */
      for (i=0; i < numantigens; i++)
	{
	  overlap = 1;
	  rza = -L[2]*0.5;      
	  while (overlap)
	    {
	      rxa = (ranf()-0.5)*L[0];
	      rya = (ranf()-0.5)*L[1];
	      /* check overlaps between antigens */
	      overlap = 0;
	      for (k = 0; k < i-1; k++)
		{
		  if (Sqr(rx[parnum+i]-rxa) + Sqr(ry[parnum+i]-rya) < Sqr(DiamAntigen))
		    {
		      overlap=1;
		      break;
		    }
		}	  
	    }
	  //rx[i+parnum] = rxa;
	  ry[i+parnum] = rya;
	  rz[i+parnum] = rza;
	  fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 3\n", rxa, rya, rza); 
	} 
    }
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "0 0 0 1 0 0 0 1 0 0 0 1 4\n"); 

  /* velocities */
  for (i=0; i < parnum; i++)
    {
      fprintf(f, "%f %f %f 0 0 0\n", ranf(), ranf(), ranf());
    }
  for (i=0; i < numantigens; i++)
    {
      fprintf(f, "0 0 0 0 0 0\n");
    }
  if (bigAntigenSurfDiam > 0.0)
    fprintf(f, "0 0 0 0 0 0\n");

  fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  printf("Number of Nanobodies=%d - Number of Antigens=%d\n", numpoly, numantigens);
  printf("polylen=%d\n", polylen);
  //printf("phi=%f\n", parnum*vol/(L[0]*L[1]*L[2]));
  fclose(f);
  return 1;
} 
