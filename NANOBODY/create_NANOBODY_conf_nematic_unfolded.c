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
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */

double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}


int main(int argc, char **argv)
{
  FILE *f;
  double sigChain, DiamSph, orient, theta0, theta0rad, Diam, del0, del0x, del0y, del0z, maxL, pi;
  double vol, permdiam, thmax, del, sigb, delfb1, delfb2, delfb3, delfb4, Len, sigmaAntigens, sigSph;
  double del00x, del00y, del00z, *rxCM, *ryCM, *rzCM, bs[3], factor[3], deltax, deltay, deltaz, DiamAntigen;
  double Rcmx, Rcmy, Rcmz, rxa, rya, rza;
  double phi, targetphi=0.25, xtrafact, Lx=10.0, Ly=10.0, Lz=10.0;
  int k1, k2, numpoly, parnum=1000, i, j, polylen, a, b, numSpheres=3;
  int type, kk, k, overlap, nx, ny, nz, nxmax, nymax, nzmax, idx, numantigens;
  del=0.5;
  /* nanobody (Fab) permanent spots diameter */
  permdiam=0.8; 
  sigSph = 0.119;
  DiamAntigen = 0.395;
  //permdiam=0.5;
  if (argc == 1)
   {
     printf("create_NANOBODY_conf_nematic_unfolded <conf_file_name> <nxmax> <nymax> <nzmax> <Lx> <Ly> <Lz> <DensSuperfAntigens> <numspheres> <Fab-patch-diam> <sphere-patch-diam> <DiametroAntigene>\n"); 
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

  if (argc > 10)
    permdiam = atof(argv[10]);

  if (argc > 11)
    sigSph = atof(argv[11]);
  
  if (argc > 12)
    DiamAntigen = atof(argv[12]);

  polylen = numSpheres + 2;
  parnum = polylen*nxmax*nymax*nzmax;
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
  /* center the nanobody */
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

  bs[0] = Diam>DiamSph?Diam:DiamSph;
  bs[1] = Diam>DiamSph?Diam:DiamSph; 
  bs[2] = 2.0*(Len+deltaz)+numSpheres*(DiamSph+deltaz);
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
  fprintf(f, "parnum: %d\n", parnum + numantigens);
  fprintf(f,"ninters: 9\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 4\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d %d %d %d\n", numpoly, numpoly, numSpheres*numpoly, numantigens);
  /* first Fab */
  fprintf(f,"%f %f %f\n", Diam/2.0, Diam/2.0, Len/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n", Len/2.0, permdiam);/* 0: along z axis (permanent) 0 */
  fprintf(f,"0 0 %f %f\n", -(Len/2.0-0.5), 0.612); /* 1: along z axis patch which will form bonds with antigens */

  /* second Fab */
  fprintf(f,"%f %f %f\n", Diam/2.0, Diam/2.0, Len/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n", -Len/2.0, permdiam);/* 0: along z axis (permanent) 0 */
  fprintf(f,"0 0 %f %f\n", (Len/2.0-0.5), 0.612); /* 1: along z axis patch which will form bonds with antigens */

  /* sphere */
  fprintf(f,"%f %f %f\n", DiamSph/2.0, DiamSph/2.0, DiamSph/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1 1 1 1 1 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"0 0 %f %f\n", DiamSph/2.0, sigSph);/* 0: along z axis (permanent) 0 */
  fprintf(f,"0 0 %f %f\n", -DiamSph/2.0, sigSph); /* 1: along z axis patch which will form bonds with antigens */

  /* antigen */
  fprintf(f,"%f %f %f\n", DiamAntigen/2.0, DiamAntigen/2.0, DiamAntigen/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"1 1 1\n");
  fprintf(f, "1e+200 1 1 1 0 1\n");
  fprintf(f,"1 0\n");
  fprintf(f,"0 0 0 %f\n", DiamAntigen);

  /* interactions */
  fprintf(f,"0 0 2 0 1 1000000 100000 1\n");
  fprintf(f,"1 0 2 0 1 1000000 100000 1\n");
  fprintf(f,"0 0 2 1 1 1000000 100000 1\n");
  fprintf(f,"1 0 2 1 1 1000000 100000 1\n");
  fprintf(f,"2 0 2 0 1 1000000 100000 1\n");
  fprintf(f,"2 0 2 1 1 1000000 100000 1\n");
  fprintf(f,"2 1 2 1 1 1000000 100000 1\n");
  fprintf(f,"0 1 3 0 1 0 0 1\n");
  fprintf(f,"1 1 3 0 1 0 0 1\n");
  fprintf(f, "@@@\n");

  /* place nanobodies */
  for (i=0; i < parnum; i++)
    {
      rx[i] -= L[0]*0.5;
      ry[i] -= L[1]*0.5;
      rz[i] -= L[2]*0.5;
      //printf("qui i=%d\n", i);
      //orient=(i%2==0)?1.0:-1.0;
      kk = i % polylen;
      if (kk==0) 
	type = 0;
      else if (kk==polylen-1) 
	type = 1;
      else
	type = 2;
      fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %d\n", rx[i], ry[i], rz[i], Ri[0][0][i], Ri[0][1][i], Ri[0][2][i], Ri[1][0][i], Ri[1][1][i], Ri[1][2][i], Ri[2][0][i], Ri[2][1][i], Ri[2][2][i], type);
      //fprintf(f, "%f %f %f  0\n", rx[i], ry[i], rz[i]);
      //printf("qui2\n");
    }
#endif	
  /* place antigens */
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
      rx[i+parnum] = rxa;
      ry[i+parnum] = rya;
      rz[i+parnum] = rza;
      fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 3\n", rxa, rya, rza); 
    } 
  /* velocities */
  for (i=0; i < parnum; i++)
    {
      fprintf(f, "%f %f %f 0 0 0\n", ranf(), ranf(), ranf());
    }
  for (i=0; i < numantigens; i++)
    {
      fprintf(f, "0 0 0 0 0 0\n");
    }
  fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  printf("Number of Nanobodies=%d - Number of Antigens=%d\n", numpoly, numantigens);
  //printf("phi=%f\n", parnum*vol/(L[0]*L[1]*L[2]));
  fclose(f);
  return 1;
} 
