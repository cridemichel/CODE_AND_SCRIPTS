#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#define MD_INF_ITENS 1E199
#define MD_INF_MASS  1E199
double nx, ny, nz, L[3], *rx, *ry, *rz, extradel;
//double *rxc, *ryc, *rzc; 
double rxl, ryl, rzl, drx, dry, drz;
double *rxCM, *ryCM, *rzCM, *vx, *vy, *vz, *wx, *wy, *wz, uinner, uouter;
double *Rc[3][3], *Ri[3][3];
int full, ibeg, numpoly, *top, nspots, nlines;
double DiamSpot=10.0;
char spotfile[1024], line[4096];
#define maxpolylen 10000;
#define Sqr(VAL_) ( (VAL_) * (VAL_) ) /* Sqr(x) = x^2 */
double R0[3][3];
struct onepart {
   double rx; 
   double ry;
   double rz;
   double R[3][3];
   double vx;
   double vy;
   double vz;
   double wx;
   double wy;
   double wz;
   int t;	
};
struct onespot
{
  double rx; 
  double ry;
  double rz;
  double sig;
};
struct onespot *allspots;
struct onepart *allpart;
double ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return drand48();
}
int compare_func (const void *aa, const void *bb)
{
  struct onepart *a, *b;
  int temp;
  a = (struct onepart*) aa;
  b = (struct onepart*) bb;
  temp = b->t - a->t;
  if (temp < 0)
    return 1;
  else if (temp > 0)
    return -1;
  else
    return 0;
}


/* ============================= >>> gauss <<< ============================= */
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

  for (i=0; i < 12; i++)
    {
      sum = sum + ranf();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}
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

double max2(double a, double b)
{
  double m;
  m = a;
  if (b > m)
    m = b;
  return m;
}


int main(int argc, char **argv)
{
  FILE *f, *sf;
  double MoI1, MoI2, mass, massE, massS, massP, massC, rTemp, temp=1.0;
  double vxCM, vyCM, vzCM, mCM=0, mean;
  double orient, theta0, theta0rad, Diam, DiamE, DiamS, DiamP, DiamC, LenE, LenS, LenP, LenC,
	 del0, del0x, del0y, del0z, maxL, pi;
  double vol, permdiam, thmax, del, sigb, delfb1, delfb2, delfb3, delfb4, Len, volE, volS, volC;
  double del00x, del00y, del00z, *rxCM, *ryCM, *rzCM, bs[3], factor[3], delta;
  double phi, targetphi=0.25, xtrafact, LBOX=-1.0;
  int k1, k2, numpoly, parnum=1000, i, j, polylen=1, a, b, parnumE=0, parnumS=0, typeP, parnumC=0;
  int nx, ny, nz, nxmax, nymax, nzmax, idx, nE, nS, nC, done, xi;
  del=0.5;
  /* permanent spots diameter */
  permdiam=6.3; /* 20 T * 0.63 nm dove 0.63 nm è la lunghezza per base stimata per ssDNA in BiophysJ 86, 2630 (2004) */  
  if (argc == 1)
   {
     printf("create_SURFANTIG_conf <conf_file_name> <nxmax> <nymax> <nzmax> <phi> <diamE> <diamS> <diamC> <LenE> <LenS> <LenC> <# spots> <# substrates> <diamspots> <spot file>\n"); 
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
  
  if (targetphi > 1.0)
    LBOX=targetphi;

  DiamE=DiamS=DiamP=DiamC=1.0;
  if (argc > 6)
    DiamE = atof(argv[6]);

  if (argc > 7)
    DiamS = atof(argv[7]);

#if 0
  if (argc > 8)
    DiamP = atof(argv[8]);
#endif
  
  if (argc > 8)
    DiamC = atof(argv[8]);

  LenE=DiamE;
  LenS=DiamS;
  //LenP=DiamP;
  LenC=DiamC;
  
  if (argc > 9)
    LenE = atof(argv[9]);

  if (argc > 10)
    LenS = atof(argv[10]);

#if 0
  if (argc > 12)
    LenP = atof(argv[12]);
#endif

  if (argc > 11)
    LenC = atof(argv[11]);

  //LenP = LenS;
  //DiamP = DiamS;
  parnumE=1;
  if (argc > 12)
    nspots = atoi(argv[12]);
  if (argc > 13)
    parnumS = atoi(argv[13]);

  uinner=1.0;
  uouter=1.0;
  if (argc > 14)
    DiamSpot = atof(argv[14]);
  if (argc > 15)
    strcpy(spotfile, argv[15]);
  else 
    strcpy(spotfile, "spots.xyz");

  sf = fopen(spotfile, "r");
  nlines=0;
  while (!feof(sf))
    {
      fscanf(sf, "%[!\n]\n", line);
      nlines++;
    }
  fclose(sf);
  nspots = nlines;
  allspots = malloc(sizeof(struct onespot)*nspots);
  sf = fopen(spotfile, "r");
  for (i=0; i < nspots; i++)
    {
      fscanf(sf, "%lf %lf %lf\n", &(allspots[i].rx), &(allspots[i].ry), &(allspots[i].rz));
      allspots[i].sig = DiamSpot;
    }
  fclose(sf);
  printf("uinner=%f uouter=%f\n", uinner, uouter);	
  Len=max2(LenS,LenC);
  Diam=max2(DiamS,DiamC);

  printf("Len=%f Diam=%f LenE=%f DiamE=%f LenC=%f DiamC=%f LenS=%f DiamS=%f\n", Len, Diam, LenE, DiamE, LenC, DiamC, LenS, DiamS);
  parnum = polylen*nxmax*nymax*nzmax;
  numpoly = nxmax*nymax*nzmax;
  rx = malloc(sizeof(double)*parnum);
  ry = malloc(sizeof(double)*parnum);
  rz = malloc(sizeof(double)*parnum);
  top = malloc(sizeof(int)*parnum);
  rxCM = malloc(sizeof(double)*numpoly);
  ryCM = malloc(sizeof(double)*numpoly);
  rzCM = malloc(sizeof(double)*numpoly);
  //rxc =malloc(sizeof(double)*polylen);
  //ryc =malloc(sizeof(double)*polylen);
  //rzc =malloc(sizeof(double)*polylen);
  vx = malloc(sizeof(double)*parnum);
  vy = malloc(sizeof(double)*parnum);
  vz = malloc(sizeof(double)*parnum);
  wx = malloc(sizeof(double)*parnum);
  wy = malloc(sizeof(double)*parnum);
  wz = malloc(sizeof(double)*parnum);
 	
  allpart = malloc(sizeof(struct onepart)*parnum);

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
  delta = 0.01;
  /* building the dimer... */
#if 0
  rxc[0] = 0.0;
  ryc[0] = 0.0;
  rzc[0] = 0.0;
#endif
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][0] = R0[a][b];
      }

#if 0
  rxc[1] = Len/2.0+delta;
  ryc[1] = Diam/2.0+delta;
  rzc[1] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][1] = R0[a][b];
      }

  rxc[2] = -Len/2.0-delta;
  ryc[2] = -Diam/2.0-delta;
  rzc[2] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][2] = R0[a][b];
      }
  Rc[0][0][2] = -Rc[0][0][2];
  rxc[3] = -Len/2.0-delta;
  ryc[3] = Diam/2.0+delta;
  rzc[3] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Rc[a][b][3] = R0[a][b];
      }
  Rc[0][0][3] = -Rc[0][0][3];
#endif
  bs[0] = Len+2.0*delta;
  bs[1] = Diam+2.0*delta;
  bs[2] = Diam+2.0*delta;
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
  fprintf(f, "ninters: 2\n");
  fprintf(f, "nintersIJ: 0\n");
  /* #Enzymes #Substrates #Products at t=0 */
  if (parnumS > 0 && parnumS < parnum-parnumE)
    {
#if 0
      parnumS=parnum-parnumE;
      parnumS = (int)(((double)parnumS)/(1.+fractC));
      parnumC = parnum-parnumE-parnumS;
#else
      parnumC=parnum-parnumE-parnumS;
#endif
    }
  else
    {
      parnumS=parnum-parnumE;
    }
  if (parnumC > 0)
    fprintf(f, "ntypes: 3\n");
  else
    fprintf(f, "ntypes: 2\n");
  fprintf(f, "saveBonds: 0\n");
  fprintf(f, "@@@\n");
    printf("parnumE=%d parnumC=%d \n", parnumE, parnumC);
  if (parnumC > 0)
    fprintf(f, "%d %d %d\n", parnumE, parnumS, parnumC);
  else
    fprintf(f, "%d %d\n", parnumE, parnumS);
  fprintf(f,"%f %f %f\n", LenE/2.0, DiamE/2.0, DiamE/2.0); 
  fprintf(f,"2 2 2\n");
  /* set here moment of inertia of a uniaxial ellipsoid */
  mass = massE = MD_INF_MASS;// 1.0*pow(DiamE/DiamS,3.0);
  MoI1 = MD_INF_ITENS;
  MoI2 = MD_INF_ITENS;
  fprintf(f, "%f %f %f %f 1 0\n", mass, MoI1, MoI2, MoI2);
// spots of giant sphere
  fprintf(f, "%d 0\n", nspots);
  for (i=0; i < nspots; i++)
    {
      fprintf(f, "0 0 0 %f\n", DiamSpot);
    }
#if 0
  fprintf(f,"2 0\n");
  fprintf(f, "0 0 0 %f\n", DiamE+0.001);
  fprintf(f, "0 0 0 %f\n", DiamE*1.1);
#endif

  fprintf(f,"%f %f %f\n", LenS/2.0, DiamS/2.0, DiamS/2.0);
  fprintf(f,"2 2 2\n");
  /* set here moment of inertia of spheres */
  mass = massS = 1.0;
  MoI1 = 1.0;
  MoI2 = 1.0;
  fprintf(f, "%f %f %f %f 1 0\n", mass, MoI1, MoI2, MoI2);
  fprintf(f,"1 0\n");
  //fprintf(f, "0 0 0 %f\n", DiamS+0.001);
  fprintf(f, "0 0 0 %f\n", DiamS*1.1);

  if (parnumC > 0)
    {
      fprintf(f,"%f %f %f\n", LenC/2.0, DiamC/2.0, DiamC/2.0); 
      fprintf(f,"2 2 2\n");
      /* set here moment of inertia of spheres */
      mass = massC = 1.0*pow(DiamC/DiamS,3.0);
      MoI1 = 1.0;
      MoI2 = 1.0;
      fprintf(f, "%f %f %f %f 1 0\n", mass, MoI1, MoI2, MoI2);
      fprintf(f,"0 0\n");
    }

#if 0
  fprintf(f,"%f 0 0 %f\n", permdiam*0.5+Len/2.0, permdiam);/* 0: along x axis (permanent) 0.05 means lp=20 */
  fprintf(f,"%f 0 0 %f\n", -Len/2.0-0.15, 0.5);
  fprintf(f,"0 0 0 0 1 0 0 1\n");
  fprintf(f,"0 1 0 1 1 0 0 100000\n");
#endif 
// interactions between spots and ligands here
  for (i=0; i < nspots; i++)
    {
      fprintf(f, "0 %d 1 0 0 %f 0 1\n", i, uouter);
    }
//  fprintf(f, "0 0 1 0 0 %f 0 1\n", uinner);
// fprintf(f, "0 1 1 1 %f 0 0 1\n", uouter);
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
  drx = 1.0;//bs[0]; //0.500000000001;
  dry = 1.0;//bs[1];
  drz = 1.0;//bs[2];
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
  printf("nx=%d %d %d dr=%f %f %f\n", nxmax, nymax, nzmax, drx, dry, drz);
  printf("step #1 L=%f %f %f\n", L[0], L[1], L[2]);
  /* expand system according to tetramer size bs[3] */
#if 1
  factor[0] = 1.0001*(bs[0]);
  factor[1] = 1.0001*(bs[1]);
  factor[2] = 1.0001*(bs[2]);
  for (i=0; i < numpoly; i++)
    {
      rxCM[i] *= factor[0];
      ryCM[i] *= factor[1];
      rzCM[i] *= factor[2];
    } 
  for (a=0; a < 3; a++)
    L[a] *= factor[a];
  printf("step #2 [factor] L=%f %f %f\n", L[0], L[1], L[2]);
#endif
  volE = M_PI*4.0*(DiamE/2.)*(DiamE/2.)*(LenE/2.)/3.0;
  volS = M_PI*4.0*(DiamS/2.)*(DiamS/2.)*(LenS/2.)/3.0;
  volC = M_PI*4.0*(DiamC/2.)*(DiamC/2.)*(LenC/2.)/3.0;
  printf("Diam=%f Len=%f\n", Diam, Len);
  if (parnumC == 0)
    {
      if (LBOX > 0.0)
	targetphi=(parnumE*volE+parnumS*volS)/(LBOX*LBOX*LBOX);
      phi=(parnumE*volE+parnumS*volS)/(L[0]*L[1]*L[2]);
    }
  else
    {
      if (LBOX > 0.0)
	targetphi=(parnumE*volE+parnumS*volS+parnumC*volC)/(LBOX*LBOX*LBOX);
      phi=(parnumE*volE+parnumS*volS+parnumC*volC)/(L[0]*L[1]*L[2]);
    }
  printf("LBOX=%f targetphi=%f\n", LBOX, targetphi);
  xtrafact = pow(phi/targetphi, 1.0/3.0);
  printf("volE=%f volS=%f volC=%f targetphi=%f phi=%f xtrafact=%f\n", volE, volS, volC, targetphi, phi, xtrafact);

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
  /* immobile particle */
  rx[0] = ry[0] = rz[0] = 0.0;
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	Ri[a][b][0] = (a==b)?1.0:0.0;
	allpart[i].R[a][b] = Ri[a][b][0];
      }
	 
  idx = 1;
  while (idx < parnum)
    {
      i = (int) (numpoly*ranf());
      /* l'enzima gigante è in (0, 0, 0) */
      if (Sqr(rxCM[i])+Sqr(ryCM[i])+Sqr(rzCM[i]) > Sqr((DiamE + DiamS)*0.5) )
	{
	  //printf("idx=%d i=%d\n", idx, i);
	  rx[idx] = rxCM[i];
	  ry[idx] = ryCM[i];
	  rz[idx] = rzCM[i];
	  for (a=0; a < 3; a++)
	    for (b=0; b < 3; b++)
	      {
		Ri[a][b][idx] = Rc[a][b][0];
		allpart[idx].R[a][b] = Rc[a][b][0];
	      }
	  idx++;
	} 
      //printf("qui idx=%d\n", idx);
    }
  nE = nS = nC = 0;
  top[0] = 0;
  allpart[0].t = 0;
       
  for (i=1; i < parnum; i++)
    {
      rx[i] -= L[0]*0.5;
      ry[i] -= L[1]*0.5;
      rz[i] -= L[2]*0.5;
      //printf("qui i=%d\n", i);
      //orient=(i%2==0)?1.0:-1.0;
      if (parnumC == 0 && parnumE==1)
	{
	  if (i < parnumE)
	    typeP = 0;
	  else
	    typeP = 1;
	}
      else
	{
	  done = 0;
	  while (!done)
	    {
	      xi = drand48()*parnum;
	      if (nE < parnumE && xi < parnumE)
		{
		  typeP = 0;
		  nE++;
		  done = 1;
		}
	      else if (nS < parnumS && xi < parnumS+parnumE)
		{
		  typeP = 1;
		  nS++;
		  done = 1;
		}
	      else if (nC < parnumC)
		{
		  nC++;
		  typeP = 3;
		  done = 1;
		}
	      else if (nE+nS+nC==parnum)
		done = 1;
	    }
	}
      top[i] = typeP;
      allpart[i].t = top[i];
        //fprintf(f, "%f %f %f  0\n", rx[i], ry[i], rz[i]);
      //printf("qui2\n");
    }
#endif	
  /* velocities */
  for (i=0; i < parnum; i++)
    {
      //rx[i] -= L[0]*0.5;
      //ry[i] -= L[1]*0.5;
      //rz[i] -= L[2]*0.5;
      //printf("qui i=%d\n", i);
      //orient=(i%2==0)?1.0:-1.0;
      allpart[i].rx = rx[i];
      allpart[i].ry = ry[i];
      allpart[i].rz = rz[i];
      if (parnumE==1 && parnumC == 0.0)
	{
	  if (i < parnumE)
	    {
	      mass = massE;
	    }
	  else
	    {
	      mass = massS;
	    }
	}
      else
	{
	  if (top[i]==0)
	    mass = massE;
	  else if (top[i]==1)
	    mass = massS;
	  else if (top[i]==3)
	    mass = massC;
	}
      if (i==0)
	{
	  vx[i] = vy[i] = vz[i] = 0.0;
	}
      else
	{
	  rTemp = sqrt(temp / mass);  
	  vx[i] = rTemp * gauss(); 
	  vy[i] = rTemp * gauss();
	  vz[i] = rTemp * gauss();
	}
    }
  vxCM=vyCM=vzCM=0.0;
  for (i=0; i < parnum; i++)
    {
      if (parnumE==1 && parnumC ==0.0)
	{
	  if (i < parnumE)
    	    mass = massE;
	  else 
	    mass = massS;
	}
      else
	{
	  if (top[i]==0)
	    mass = massE;
	  else if (top[i]==1)
	    mass = massS;
	  else if (top[i]==3)
	    mass = massC;
	}
      mCM += mass;
      if (i > 0)
	{
	  vxCM+=mass*vx[i];
	  vyCM+=mass*vy[i];
	  vzCM+=mass*vz[i];
	}
    }
  for (i=0; i < parnum; i++)
    {
      if (i > 0)
	{
	  vx[i] -= vxCM/mCM;
	  vy[i] -= vyCM/mCM;
	  vz[i] -= vzCM/mCM;
	}
      allpart[i].vx = vx[i];
      allpart[i].vy = vy[i];
      allpart[i].vz = vz[i];
#if 1
      if (LenE==DiamE && LenS==DiamS && LenC==DiamC)
	wx[i]=wy[i]=wz[i]=0;
      else
	{
	  if (i==0)
	    {
	      wx[i] = wy[i] = wz[i] = 0.0;
	    }
	  else
	    {
	      if (i < parnumE)
		mean = sqrt(temp / MoI1);
	      else if (i < parnumS+parnumE)
		mean = sqrt(temp / MoI1);
	      else 
		mean = sqrt(temp / MoI1);
	      wx[i] = mean*gauss();
	      wy[i] = mean*gauss();
	      wz[i] = mean*gauss();
	    }
	}
#else
      wx[i]=wy[i]=wz[i]=0;
#endif
      allpart[i].wx = wx[i];
      allpart[i].wy = wy[i];
      allpart[i].wz = wz[i];
     }  
   /* sort particles here (enzyme first!) */
   qsort(allpart, parnum, sizeof(struct onepart), compare_func);
   for (i=0; i < parnum; i++)
     {
       fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %d\n", allpart[i].rx, allpart[i].ry, allpart[i].rz, 
	allpart[i].R[0][0], allpart[i].R[0][1], allpart[i].R[0][2], allpart[i].R[1][0], 
	allpart[i].R[1][1], allpart[i].R[1][2], 
	allpart[i].R[2][0], allpart[i].R[2][1], allpart[i].R[2][2], allpart[i].t);
     }
   for (i=0; i < parnum; i++)
     {
 
      fprintf(f, "%f %f %f %f %f %f\n", allpart[i].vx, allpart[i].vy, allpart[i].vz, 
	allpart[i].wx, allpart[i].wy, allpart[i].wz);
     }
  fprintf(f, "%.15G %.15G %.15G\n", L[0], L[1], L[2]);

  printf("phi=%f\n", (parnumE*volE+parnumS*volS+parnumC*volC)/(L[0]*L[1]*L[2]));
  fclose(f);
  return 1;
} 
