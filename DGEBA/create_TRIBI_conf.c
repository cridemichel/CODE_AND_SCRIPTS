#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//#define HARD_SPHERES
double nx, ny, nz, L[3], *rx, *ry, *rz, extradel;
double *rxc, *ryc, *rzc, rxl, ryl, rzl, drx, dry, drz, *vx, *vy, *vz, *typeofpar;
double *rxCM, *ryCM, *rzCM;
double *Rc[3][3], *Ri[3][3];
int full, ibeg, numpoly;
#define maxpolylen 10000;
double R0[3][3];
double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
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

  double a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0; i < 12; i++)
    {
      sum = sum + ranf();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}

int main(int argc, char **argv)
{
  FILE *f;
  double orient, theta0, theta0rad, Diam, DiamA, DiamB, volA, volB, del0, del0x, del0y, del0z, maxL, pi;
  double vol, permdiam, thmax, del, sigb, delfb1, delfb2, delfb3, delfb4, Len;
  double del00x, del00y, del00z, *rxCM, *ryCM, *rzCM, bs[3], factor[3], delta, MoI, MoIA, MoIB;
  double phi, targetphi=0.25, xtrafact, pD, temp, rTemp, rTempA, rTempB, wx, wy, wz, mean, mass, totvx, totvy, totvz;
  int k1, k2, numpoly, parnum=1000, i, j, polylen=1, a, b, parnumA, parnumB;
  int nx, ny, nz, nxmax, nymax, nzmax, idx, type, numtypeA, numtypeB;
  del=0.5;
  /* permanent spots diameter */
  permdiam=6.3; /* 20 T * 0.63 nm dove 0.63 nm è la lunghezza per base stimata per ssDNA in BiophysJ 86, 2630 (2004) */  
  Len=16.0; /* 48 bp equal roughly to 16 nm */
  if (argc == 1)
   {

     printf("create_GDNA_conf <conf_file_name> <nxmax> <nymax> <nzmax> <phi> <diam>\n"); 
     exit(-1);
   }
  srand48(-1);
  f = fopen(argv[1], "w+");
  if (f==NULL)
    {
      printf("boh...f is null!\n");
      exit(-1);
    }
  DiamA=1.0;
  DiamB=2.0;
  Diam=DiamB;
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

  parnum = nxmax*nymax*nzmax;
  parnumA = 2*parnum/5;
  parnumB = 3*parnum/5;
  numpoly = nxmax*nymax*nzmax;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  vx = malloc(sizeof(double)*parnum);
  vy = malloc(sizeof(double)*parnum);
  vz = malloc(sizeof(double)*parnum);
  typeofpar = malloc(sizeof(int)*parnum);
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
  delta = 0.1;
  /* building the dimer... */
  rxc[0] = 0.0;
  ryc[0] = 0.0;
  rzc[0] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][0] = R0[a][b];
      }

  bs[0] = DiamB;
  bs[1] = DiamB;
  bs[2] = DiamB;
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
#ifdef HARD_SPHERES
  fprintf(f,"ninters: 0\n");
#else
  fprintf(f,"ninters: 6\n");
#endif
  fprintf(f,"nintersIJ: 0\n");
  fprintf(f,"ntypes: 2\n");
  fprintf(f,"saveBonds: 0\n");
  fprintf(f, "@@@\n");
  fprintf(f, "%d %d\n", parnumA, parnumB);
  fprintf(f,"%f %f %f\n", DiamA/2.0, DiamA/2.0, DiamA/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"2 2 2\n");
  MoIA=MoI=2.0/3.0*(DiamA/2.0)*(DiamA/2.0);
  pD = 0.11965683746373801*(DiamA+DiamB)/2.;
  fprintf(f, "1 %f %f %f 2 0\n", MoI, MoI, MoI);
#ifdef HARD_SPHERES
  fprintf(f, "0 0\n");
#else
  fprintf(f,"3 0\n");
  fprintf(f,"%f 0 0 %f\n", DiamA/2.0, pD);/* 0: along x axis (permanent) 0.05 means lp=20 */
  fprintf(f,"%f %f 0 %f\n", -DiamA*cos(pi/6.0), -DiamA*sin(pi/6.0), pD);
  fprintf(f,"%f %f 0 %f\n", DiamA*cos(pi/6.0), -DiamA*sin(pi/6.0), pD);
  fprintf(f,"%f %f %f\n", DiamB/2.0, DiamB/2.0, DiamB/2.0); /* each dsDNA of 48 bp which is roughly equal to 48 / 3 nm = 16 nm (D=2 nm in our case) */ 
  fprintf(f,"2 2 2\n");
  MoIB=MoI=2.0/3.0*(DiamB/2.0)*(DiamB/2.0);
  fprintf(f, "8 %f %f %f 2 0\n", MoI, MoI, MoI);
  fprintf(f,"2 0\n");
  fprintf(f,"%f 0 0 %f\n", DiamB/2.0, pD);/* 0: along x axis (permanent) 0.05 means lp=20 */
  fprintf(f,"%f 0 0 %f\n", -DiamB/2.0, pD);
  fprintf(f,"0 0 1 0 1 25 1000000 1\n");
  fprintf(f,"0 0 1 1 1 25 1000000 1\n");
  fprintf(f,"0 1 1 0 1 25 1000000 1\n");
  fprintf(f,"0 1 1 1 1 25 1000000 1\n");
  fprintf(f,"0 2 1 0 1 25 1000000 1\n");
  fprintf(f,"0 2 1 1 1 25 1000000 1\n");
#endif
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
  L[0] = (nxmax)*dry;
  L[1] = (nymax)*drx;
  L[2] = (nzmax)*drz;
#if 0
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
#endif
  //vol = (4.0*pi/3.0)*(Diam/2)*(Diam/2)*(Diam/2.0);
  volA=(4.0*pi/3.0)*(DiamA/2)*(DiamA/2)*(DiamA/2.0);
  volB=(4.0*pi/3.0)*(DiamB/2)*(DiamB/2)*(DiamB/2.0); 
  phi=(parnumA*volA + parnumB*volB)/(L[0]*L[1]*L[2]);
  xtrafact = pow(phi/targetphi, 1.0/3.0);
  printf("volA=%f targetphi=%f phi=%f xtrafact=%f\n", volA, targetphi, phi, xtrafact);

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
  numtypeA=numtypeB=0;
  for (i=0; i < numpoly; i++)
    typeofpar[i] = -1;
   
  for (i=0; i < parnum; i++)
    {
      rx[i] -= L[0]*0.5;
      ry[i] -= L[1]*0.5;
      rz[i] -= L[2]*0.5;
      //printf("qui i=%d\n", i);
      //orient=(i%2==0)?1.0:-1.0;
      type = (drand48() > 0.4)?1:0;
      if (type==0)
	{
	  numtypeA++;
	  if (numtypeA > parnumA)
	    type=1;
	}
      else
	{	
	  numtypeB++;
	  if (numtypeB > parnumB)
	    type=0;
	}
      typeofpar[i] = type;
      fprintf(f, "%.15G %.15G %.15G %f %f %f %f %f %f %f %f %f %d\n", rx[i], ry[i], rz[i], Ri[0][0][i], Ri[0][1][i], Ri[0][2][i], Ri[1][0][i], Ri[1][1][i], Ri[1][2][i], Ri[2][0][i], Ri[2][1][i], Ri[2][2][i], type);
      //fprintf(f, "%f %f %f  0\n", rx[i], ry[i], rz[i]);
      //printf("qui2\n");
    }
#endif	
  temp=1.0;
  rTempA=sqrt(temp);
  rTempB=sqrt(temp/8.);
  
  totvx = totvy = totvz = 0;
  for (i=0; i < parnum; i++)
    {
      if (typeofpar[i] == 0)
	{
	  rTemp = rTempA;
	  mass = 1.0;
	}
      else
	{
	  mass = 8.0;
	  rTemp = rTempB;
	}
      vx[i] = rTemp*gauss();
      vy[i] = rTemp*gauss();
      vz[i] = rTemp*gauss();
      totvx += mass*vx[i];
      totvy += mass*vy[i];
      totvz += mass*vz[i];
    }
  printf("Mtot=%f totvx=%f totvy=%f totvz=%f\n", parnumA + parnumB*8., totvx, totvy, totvz);
  totvx /= (parnumA + parnumB*8.0);
  totvy /= (parnumA + parnumB*8.0);
  totvz /= (parnumA + parnumB*8.0);
  for (i=0; i < parnum; i++)
    {
      vx[i] -= totvx;
      vy[i] -= totvy;
      vz[i] -= totvz;
    }
  totvx=totvy=totvz=0;
  for (i=0; i < parnum; i++)
    {
      if (typeofpar[i] == 0)
	{
	  rTemp = rTempA;
	  mass = 1.0;
	}
      else
	{
	  mass = 8.0;
	  rTemp = rTempB;
	} 
      totvx += mass*vx[i];
      totvy += mass*vy[i];
      totvz += mass*vz[i];
    } 
  printf("parnum=%d parnumA=%d parnumB=%d\n", parnum, parnumA, parnumB);
  totvx /= (parnumA + parnumB*8.0);
  totvy /= (parnumA + parnumB*8.0);
  totvz /= (parnumA + parnumB*8.0);
  printf("totv= %.15G %.15G %.15G\n", totvx, totvy, totvz); 
  for (i=0; i < parnum; i++)
    {
      if (typeofpar[i]==0)
	mean = sqrt(temp / MoIA);
      else
	mean = sqrt(temp / MoIB);
      wx = mean*gauss();
      wy = mean*gauss();
      wz = mean*gauss();
      fprintf(f, "%.15G %.15G %.15G %.15G %.15G %.15G\n", vx[i], vy[i], vz[i],
	      wx, wy, wz);
    }
  fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);
  printf("phi=%f\n", (parnumA*volA+parnumB*volB)/(L[0]*L[1]*L[2]));
  fclose(f);
  return 1;
} 
