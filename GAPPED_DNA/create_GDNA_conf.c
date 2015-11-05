#include<stdlib.h>
#include<stdio.h>
#include<math.h>
double nx, ny, nz, L, *rx, *ry, *rz, extradel;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz;
int full, ibeg, numpoly;
#define maxpolylen 10000;
double *R[3][3], R0[3][3], *Ri[3][3];
int main(int argc, char **argv)
{
  FILE *f;
  double orient, theta0, theta0rad, phi, Diam, del0, del0x, del0y, del0z, maxL, pi;
  double vol, permdiam, thmax, del, sigb, delfb1, delfb2, delfb3, delfb4, Len;
  double del00x, del00y, del00z;
  int k1, k2, numpoly, parnum=1000, i, j, polylen=2;

  del=0.5;
  /* permanent spots diameter */
  permdiam=6.3; /* 20 T * 0.63 nm dove 0.63 nm è la lunghezza per base stimata per ssDNA in BiophysJ 86, 2630 (2004) */  
  Len=32.0;
  if (argc == 1)
   {

     printf("create_GDNA_conf <conf_file_name> <num. polymers> <volume fraction> <polymer length>\n"); 
     exit(-1);
   }
  f = fopen(argv[1], "w+");
  Diam=2.0;
  pi = acos(0.0)*2.0;

  if (argc > 4) 
    polylen = atoi(argv[4]);

  numpoly=500;
  if (argc > 2)
    {
      numpoly = atoi(argv[2]);
      parnum=numpoly*polylen;
    }
  else
    parnum =500*polylen; 
  if (argc > 3)
    {
      phi = atof(argv[3]);
      vol=3.14159*(Diam*0.5)*(Diam*0.5)*Len;
      printf("parnum=%d vol=%f phi=%f\n", parnum, vol, phi);
      L=pow(((double)parnum)*vol/phi,1.0/3.0);
      //L= pow(((double)2.0),2.0/3.0)*pow(((double)numpoly),1.0/3.0)*pow(pi/3.0,1.0/3.0)*pow(((double)polylen),1.0/3.0)/pow(phi,1.0/3.0);
    }
  else
    L=50.0;
  extradel=0.5;
  if (L < Diam*polylen*1.01)
    {
      L = Diam*polylen*1.01;
    }
  printf("L=%f phi=%f argc=%d polylen=%d\n", L, phi, argc, polylen);
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  rxc =malloc(sizeof(double)*polylen);
  ryc =malloc(sizeof(double)*polylen);
  rzc =malloc(sizeof(double)*polylen);
#if 0
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Ri[k1][k2] = malloc(sizeof(double)*parnum);     
	R[k1][k2] =  malloc(sizeof(double)*polylen);
      }
#endif
  del00x = Len*0.5+0.0001;
  del00y=del00z=Diam*0.5+0.0001;
  del0x = 0.01;
  del0y=del0z=Diam*0.1+0.0001;
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R0[k1][k2]=0.0;
  R0[0][0]=R0[1][1]=R0[2][2]=1.0;
  for (i=0; i < polylen; i++)
    {
      rxc[i] = i*1.0001*Len;
      ryc[i] = 0.0;
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
  fprintf(f, "parnum: %d\n", parnum);
  fprintf(f,"ninters: 2\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"8 1 1\n"); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"2 0\n");
  fprintf(f,"%f 0 0 %f\n", permdiam*0.5+8.0, permdiam);/* 0: along x axis (permanent) 0.05 means lp=20 */
  fprintf(f,"%f 0 0 %f\n", -8.15, 0.5);
  fprintf(f,"0 0 0 0 1 0 0 1\n");
  fprintf(f,"0 1 0 1 1 0 0 100000\n");
  fprintf(f, "@@@\n");
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  maxL = L;
  numpoly=parnum/polylen;
  printf("numpoly=%d numpoly*polylen=%d\n", numpoly, numpoly*polylen);
  drx = Len*((double)polylen)+del0x; //0.500000000001;
  dry = Diam+del0y;
  drz = Diam+del0z;
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx+del00x;
      if (rxl+Len*0.5 > L)
	{
	  ny++;
	  nx=0;
	  rxl = rxc[polylen-1]+nx*drx+del00x;
	}
      ryl = ryc[polylen-1]+ny*dry+del00y;
      if (ryl+Diam*0.5 > L)
	{
	  nz++;
	  ny=0;
	  ryl = ryc[polylen-1]+ny*dry+del00y;
      	}
      rzl = rzc[polylen-1]+nz*drz+del00z;
      if (rzl+Diam*0.5 > L)
	{
	  full=1;
	}

      for (j = 0; j < polylen; j++)
	{
	  rx[i*polylen+j] = rxc[j]+nx*drx+del00x;
	  //printf("nx=%f i=%d j=%d x=%f\n", nx, i, j,  rx[i*polylen+j]);
	  ry[i*polylen+j] = ryc[j]+ny*dry+del00y;
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
  //  if (full)
  //   exit(-1);
  for (i=0; i < parnum; i++)
    {
      rx[i] -= L*0.5;
      ry[i] -= L*0.5;
      rz[i] -= L*0.5;
      orient=(i%2==0)?1.0:-1.0;
      fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f 0\n", rx[i], ry[i], rz[i],
	      orient*R0[0][0], R0[0][1], R0[0][2], R0[1][0], R0[1][1], R0[1][2],
	      R0[2][0], R0[2][1], R0[2][2]);
    }
  fprintf(f, "%.15G %.15G %.15G\n", L, L, L);
  return 1;
} 
