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
  double theta0, theta0rad, phi, Diam, del0, del0x, del0y, del0z, maxL, pi;
  double thmax, del, sigb, delfb1, delfb2, delfb3, delfb4;
  int k1, k2, numpoly, parnum=2800, i, j, polylen=20;

  del=0.5;
  thmax=(30.0/180.0)*3.14159265358979;
  if (argc == 1)
   {

     printf("create_KFpoly_conf <conf_file_name> <num. polymers> <volume fraction> <polymer length>\n"); 
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
      L= pow(((double)2.0),2.0/3.0)*pow(((double)numpoly),1.0/3.0)*pow(pi/3.0,1.0/3.0)*pow(((double)polylen),1.0/3.0)/pow(phi,1.0/3.0);
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
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Ri[k1][k2] = malloc(sizeof(double)*parnum);     
	R[k1][k2] =  malloc(sizeof(double)*polylen);
      }
  del0 = Diam*0.5+0.0001;
  del0x = 0.0001;
  del0y=del0z=Diam*2;
  theta0=10.0; /* in gradi */
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      R0[k1][k2]=0.0;
  R0[0][0]=1.0;
  for (i=0; i < polylen; i++)
    {
      rxc[i] = i*Diam;
      ryc[i] = 0.0;
      rzc[i] = 0.0;
      /*apply theta0 rotation around x-axis */
      theta0rad=pi*(i*theta0)/180.0;
      R0[1][1]=cos(theta0rad);
      R0[1][2]=-sin(theta0rad);
      R0[2][1]=sin(theta0rad);
      R0[2][2]=cos(theta0rad);
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  R[k1][k2][i] = R0[k1][k2];
    }
  fprintf(f, "parnum: %d\n", parnum);
  fprintf(f,"ninters: 11\n");
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 1\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d\n", parnum);
  fprintf(f,"1 1 1\n");
  fprintf(f,"2 2 2\n");
  fprintf(f, "1 1 1 1 2 0\n");
  fprintf(f,"6 0\n");
  fprintf(f,"1 0 0 0.119\n");/* 0: along x axis (permanent) 0.05 means lp=20 */
  fprintf(f,"-1 0 0 0.119\n");/* 1: along x axis (permanent) */
  fprintf(f,"0  1  0 0.5\n");/* 2: along y axis */ 
  fprintf(f,"0 -1  0 0.5\n");/* 3: along y axis */
  fprintf(f,"0  0  1 0.5\n");/* 4: along z axis */
  fprintf(f,"0  0 -1 0.5\n");/* 5: along z axis */
  fprintf(f,"0 0 0 0 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 0 0 1 0.0001 1000000000 10000000000 100000\n");
  fprintf(f,"0 1 0 1 0.0001 1000000000 10000000000 100000\n");
  sigb=log(2.0*(pow(1.0+del/2.0,3.0)-1.0)*(1.0-cos(thmax)));
  delfb1=3.1-sigb;	
  delfb2=0.6-sigb;
  delfb3=0.215-sigb;
  delfb4=1.25-sigb;
  fprintf(f,"0 2 0 2 %.15G 0 0 100000\n",delfb1); /* y-axis kern-frenkel patches */
  fprintf(f,"0 2 0 3 %.15G 0 0 100000\n",delfb1); /* y-axis kern-frenkel patches */
  fprintf(f,"0 3 0 3 %.15G 0 0 100000\n",delfb1); /* y-axis kern-frenkel patches */
  fprintf(f,"0 4 0 4 %.15G 0 0 100000\n",delfb2); /* z-axis KF patch */ 
  fprintf(f,"0 5 0 5 %.15G 0 0 100000\n",delfb2);  /* z-axis KF patch */
  fprintf(f,"0 4 0 5 %.15G 0 0 100000\n",delfb3);  /* z-axis KF patch */
  fprintf(f,"0 3 0 4 %.15G 0 0 100000\n",delfb4);  /* perp conf */
  fprintf(f,"0 2 0 5 %.15G 0 0 100000\n",delfb4);  /* perp conf */
  fprintf(f,"@@@\n");  
  nx=ny=nz=0;
  full=0;
  ibeg = 0;
  maxL = L;
  numpoly=parnum/polylen;
  printf("numpoly=%d numpoly*polylen=%d\n", numpoly, numpoly*polylen);
  drx = Diam*((double)polylen)+del0x; //0.500000000001;
  dry = Diam+del0y;
  drz = Diam+del0z;
  if (parnum%polylen != 0)
    {
      printf("number of particles must a multiple of %d\n", polylen);
      exit(-1);	
    }
  for (i=0; i < numpoly; i++)
    {
      rxl = rxc[polylen-1]+nx*drx+del0;
      if (rxl+Diam*0.5 > L)
	{
	  ny++;
	  nx=0;
	  rxl = rxc[polylen-1]+nx*drx+del0;
	}
      ryl = ryc[polylen-1]+ny*dry+del0;
      if (ryl+Diam*0.5 > L)
	{
	  nz++;
	  ny=0;
	  ryl = ryc[polylen-1]+ny*dry+del0;
      	}
      rzl = rzc[polylen-1]+nz*drz+del0;
      if (rzl+Diam*0.5 > L)
	{
	  full=1;
	}

      for (j = 0; j < polylen; j++)
	{
	  rx[i*polylen+j] = rxc[j]+nx*drx+del0;
	  //printf("nx=%f i=%d j=%d x=%f\n", nx, i, j,  rx[i*polylen+j]);
	  ry[i*polylen+j] = ryc[j]+ny*dry+del0;
	  rz[i*polylen+j] = rzc[j]+nz*drz+del0;
      	  for (k1=0; k1 < 3; k1++)
    	    for (k2=0; k2 < 3; k2++)
    	      Ri[k1][k2][i*polylen+j]= R[k1][k2][j];
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
      fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f 0\n", rx[i], ry[i], rz[i],
	      Ri[0][0][i], Ri[0][1][i], Ri[0][2][i], Ri[1][0][i], Ri[1][1][i], Ri[1][2][i],
	      Ri[2][0][i], Ri[2][1][i], Ri[2][2][i]);
    }
  fprintf(f, "%.15G %.15G %.15G\n", L, L, L);
  return 1;
} 
