#include<stdlib.h>
#include<stdio.h>
#include<math.h>
double nx, ny, nz, L[3], *rx, *ry, *rz, extradel;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
double *rxCM, *ryCM, *rzCM;
double *Rc[3][3], *Ri[3][3];
int full, ibeg, numpoly;
#define maxpolylen 10000;
double R0[3][3], sigma;

int main(int argc, char **argv)
{
  FILE *f;
  double Diam, del0x, del0y, del0z, maxL, pi;
  double vol, permdiam, del, Len;
  double *rxCM, *ryCM, *rzCM, bs[3], factor[3], deltax, deltay;
  double phi, targetphi=0.25, xtrafact;
  int k1, k2, numpoly, parnum=1000, i, j, polylen=2, a, b;
  int nx, ny, nz, nxmax, nymax, nzmax, idx;
  del=0.5;
  /* permanent spots diameter */
  // we want to simulate 1T and 4 T 
  permdiam=0.63/2.0; /* this parameter has to be half the length of the spacer (since the center of patch is at perdiam/2 
                    from the surface of the base, for example
                    for 20T spacesr is  20 T * 0.63 nm / 2 = 6.3 nm 
                    dove 0.63 nm è la lunghezza per base stimata per ssDNA in BiophysJ 86, 2630 (2004) */  
  //permdiam=0.5;
  Len=8.0; /* 24 bp equal roughly to 8 nm */
  Diam=2.4;
  if (argc == 1)
   {

     printf("create_gapdblhairpin_conf <conf_file_name> <nxmax> <nymax> <nzmax> <phi> <diam> <permdiam>\n"); 
     exit(-1);
   }
  f = fopen(argv[1], "w+");
  if (f==NULL)
    {
      printf("boh...f is null!\n");
      exit(-1);
    }
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

  sigma = sqrt(2.0*Diam*Diam+Len*Len);

  parnum = 2*nxmax*nymax*nzmax;
  numpoly = nxmax*nymax*nzmax;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxCM = malloc(sizeof(double)*numpoly);
  ryCM = malloc(sizeof(double)*numpoly);
  rzCM = malloc(sizeof(double)*numpoly);
  rxc = malloc(sizeof(double)*polylen);
  ryc = malloc(sizeof(double)*polylen);
  rzc = malloc(sizeof(double)*polylen);
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
  deltax = 0.05;
  /* building the dimer... */
  rxc[0] = Len/2.0+deltax;
  ryc[0] = 0.0;
  rzc[0] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][0] = R0[a][b];
      }

  //Rc[0][0][0] = -Rc[0][0][0];
  rxc[1] = -Len/2.0-deltax;
  ryc[1] = 0.0;
  rzc[1] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][1] = R0[a][b];
      }

  Rc[0][0][1] = -Rc[0][0][1];
#if 0
  rxc[2] = -Len/2.0-deltax;
  ryc[2] = Diam/2.0+deltay;
  rzc[2] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][2] = R0[a][b];
      }

  Rc[0][0][2] = -Rc[0][0][2];
  rxc[3] = Len/2.0+deltax;
  ryc[3] = Diam/2.0+deltay;
  rzc[3] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][3] = R0[a][b];
      }
  //Rc[0][0][3] = -Rc[0][0][3];
  //Rc[0][0][3] = -Rc[0][0][3];
#endif
  bs[0] = 2.0*Len+2.0*deltax;
  bs[1] = Diam;
  bs[2] = Diam;
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
  fprintf(f, "%d %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G 0\n", parnum, L[0], L[1], L[2], sigma, sigma, Len/2.0, Diam/2.0, Diam/2.0);
  fprintf(f, "%.15G 1\n", permdiam);
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
      fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f\n", rx[i], ry[i], rz[i], Ri[0][0][i], Ri[0][1][i], Ri[0][2][i], Ri[1][0][i], Ri[1][1][i], Ri[1][2][i], Ri[2][0][i], Ri[2][1][i], Ri[2][2][i]);
      //fprintf(f, "%f %f %f  0\n", rx[i], ry[i], rz[i]);
      //printf("qui2\n");
    }
#endif	
  //fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  printf("phi=%f\n", parnum*vol/(L[0]*L[1]*L[2]));
  fclose(f);
  return 1;
} 
