#include<mdsimul.h>
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#define MD_DEBUG(x) x 
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
extern double **Xa, **Xb, **RA, **RB, ***R, **Rt;
#ifdef MD_ASYM_ITENS
double **Ia, **Ib, **invIa, **invIb;
#else
double Ia, Ib, invIa, invIb;
#endif
int *lastbump;
extern double maxax[2];
/* Routines for LU decomposition from Numerical Recipe online */
void ludcmpR(double **a, int* indx, double* d, int n);
void lubksbR(double **a, int* indx, double *b, int n);
void SolveLineq (double **a, double *x, int n);
void InvMatrix(double **a, double **b, int NB);
extern double invaSq[2], invbSq[2], invcSq[2];
double rxC, ryC, rzC;
void print_matrix(double **M, int n)
{
  int k1, k2;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}
void print_matrixArr(double M[3][3])
{
  int k1, k2;
  const int n = 3;
  printf("{");
  for (k1 = 0; k1 < n; k1++)
    {
      printf("{");
      for (k2 = 0; k2 < n; k2++)
	{
	  printf("%.15G", M[k1][k2]);
	  if (k2 < n - 1)
	    printf(", ");
	}
      printf("}");
      if (k1 < n-1)
	printf(",\n");
    }
  printf("}\n");
}
#ifdef MD_GRAVITY
int checkz(char *msg)
{
  int ii;
  double sig;
  for (ii=0; ii < Oparams.parnum; ii++)
    {
      if (ii < Oparams.parnumA)
	sig = Oparams.sigma[0][0];
      else
        sig = Oparams.sigma[1][1];
      if ((rz[ii]+Lz*0.5-sig/2.0)<0. &&  
	fabs(rz[ii]+Lz*0.5-sig/2.0) > 0.1)
	{
	  printf("[%s]ii=%d diff:%.15f\n", msg,
		 ii, rz[ii]+Lz*0.5-sig/2.0);
	  printf("sotto more in checkz!!!\n");
	  return 1;
	}
    }
  return 0;
}
#endif
extern int ENDSIM;
extern char msgStrA[MSG_LEN];
extern char TXT[MSG_LEN];
extern void resetCM(int Nm);
extern void vectProd(COORD_TYPE r1x, COORD_TYPE r1y, COORD_TYPE r1z, 
	 COORD_TYPE r2x, COORD_TYPE r2y, COORD_TYPE r2z, 
	 COORD_TYPE* r3x, COORD_TYPE* r3y, COORD_TYPE* r3z);
extern void kinet(int Nm, COORD_TYPE* velx, COORD_TYPE* vely, 
		  COORD_TYPE* velz, COORD_TYPE VOL1);
extern void ScheduleEvent (int idA, int idB, double tEvent);
extern void NextEvent (void);
void distanza(int ia, int ib);
double pi, invL, L2, Vz;   
#ifdef MD_GRAVITY
double Lz2;
#endif
double W, K, T1xx, T1yy, T1zz,
  T1xx, T1yy, T1zz, T1xy, T1yz, T1zx, Wxx, Wyy, Wzz,
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB;
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
double DrSq = 0.0; 
const double timbig = 1E12;
/* used by linked list routines */
#ifdef MD_GRAVITY
extern double g2, mgA, mgB;
#endif
double *lastcol;
double *treetime, *atomTime, *rCx, *rCy, *rCz; /* rC è la coordinata del punto di contatto */
int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
  cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
int evIdA, evIdB, parnumB, parnumA;
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
int evIdC;
extern int *bondscache, *numbonds, **bonds;
#endif
/* ========================== >>> scalCor <<< ============================= */
void scalCor(int Nm)
{ 
  int i;
  double DRx, DRy, DRz, invL;
  
  invL = 1.0 / L;
  /* Reduced particles to first box */
  for(i=0; i < Oparams.parnum; i++)
    {
      /* (DRx, DRy, DRz) is the quantity to add to the positions to 
	 scale them */
      DRx = - L * rint(invL * rx[i]);
      DRy = - L * rint(invL * ry[i]);
      DRz = - L * rint(invL * rz[i]);
      rx[i] += DRx;
      ry[i] += DRy;
      rz[i] += DRz;
    }
}
#ifdef MD_GRAVITY
void calcKVz(void)
{
  int i;
  double dd;
  Vz = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      Vz += vz[i];
    }
  Vz *= Oparams.m[0];
  dd = 0.0;
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += vz[i];
    }
  dd *= Oparams.m[1];
  Vz += dd;
  Vz /= Mtot;
  K = 0.0;
  for (i=0; i < Oparams.parnumA; i++)
    {
      K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
    }
  K *= Oparams.m[0]*0.5;
  dd = 0.0;
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
    }
  dd *= Oparams.m[1]*0.5;
  K += dd;
}
#endif

void outputSummary(void)
{
  FILE *f;
  int i;
  /* mettere qualcosa qui */
#ifdef MD_GRAVITY
  printf("K= %.15f V=%.15f T=%.15f Vz: %f\n", K, V, 
	 (2.0*K/(3.0*Oparams.parnum-3.0)), Vz);
#else
#if 0
  printf("time=%.15f\n", OprogStatus.time);
  printf("K= %.15f T=%.15f\n", K, 
	 (2.0*K/(3.0*Oparams.parnum-3.0)));
#endif
#endif
#if 0
  f = fopenMPI(MD_HD_MIS "T.dat", "a");
  if (OprogStatus.brownian==1)
    fprintf(f, "%.15f %.15f\n", Oparams.time, (2.0*K/(3.0*Oparams.parnum)));
  else
    fprintf(f, "%.15f %.15f\n", Oparams.time, (2.0*K/(3.0*Oparams.parnum-3.0)));
  fclose(f);
#endif
#ifdef MD_GRAVITY
  f = fopenMPI(MD_HD_MIS "Vz2.dat", "a");
  fprintf(f, "%.15f %.15f\n", Oparams.time, Sqr(Vz));
  fclose(f);
#endif
}
#ifdef MD_GRAVITY
void scalevels(double temp, double K, double Vz)
{
  int i; 
  double sf, VVx, VVy, VVz, ddx, ddy, ddz;
  sf = sqrt( ( (3.0*((double)Oparams.parnum)-3.0) * temp ) / (2.0*K) );

  VVx = VVy = VVz = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] = (vz[i] - Vz)*sf + Vz;
      VVx += vx[i];
      VVy += vy[i];
      VVz += vz[i];
      /* scala anche i tempi di collisione! */
    } 
  VVx *= Oparams.m[0];
  VVy *= Oparams.m[0];
  VVz *= Oparams.m[0];
  ddx = ddy = ddz = 0.0;
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] = (vz[i] - Vz)*sf + Vz;
      ddx += vx[i];
      ddy += vy[i];
      ddz += vz[i];
      /* scala anche i tempi di collisione! */
    } 
  ddx *= Oparams.m[1];
  ddy *= Oparams.m[1];
  ddz *= Oparams.m[1];
  VVx += ddx;
  VVy += ddy;
  VVz += ddz;
  VVx = VVx / Mtot;
  VVy = VVy / Mtot;
  VVz = VVz / Mtot;
  for (i = 0; i < Oparams.parnum; i++)
    {
      vx[i] -= VVx;
      vy[i] -= VVy;
      vz[i] = vz[i] - VVz + Vz;
    }
  MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}
#else
void scalevels(double temp, double K)
{
  int i; 
  double sf;
    
  sf = sqrt( ( (3.0*((double)Oparams.parnum)-3.0) * temp ) / (2.0*K) );
  for (i = 0; i < Oparams.parnumA; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      /* scala anche i tempi di collisione! */
    } 
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] *= sf;
      /* scala anche i tempi di collisione! */
    } 
  MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}
#endif
/* ============================ >>> updateQ <<< =========================== */
void updateDQ(COORD_TYPE dt)
{
  /* Here we use molecular pressure tensor */
  OprogStatus.DQxy += dt*Pxy;
  OprogStatus.DQyz += dt*Pyz;
  OprogStatus.DQzx += dt*Pzx;

  /*         _ t
	    |
    DQab =  |  Pab dt
            |
	   - 0
	   
    con a, b = x,y,z 
  */
}

/* ============================= >>> updatePE <<< ========================= */
void updatePE(int Nm)
{
  int iE;
  double ENmin, ENmax;
  
  ENmin = OprogStatus.ENmin;
  ENmax = OprogStatus.ENmax;
  iE = (int) ((V - ((double) Nm)*ENmin) / 
	      ( ((double) Nm) * ((double) ENmax - ENmin)) * 
	      ((double) PE_POINTS));
  /*  sprintf(TXT, "energia: %f iE:%d\n", EE[my_rank], iE);
      mdPrintf(ALL, TXT, NULL);
      printf("PEPOINTS: %d Vc:%f\n", PE_POINTS, Vc);*/
  if ( (iE >= 0) && (iE < PE_POINTS)) 
    {
      ++(OprogStatus.PE[iE]);
    }

}

extern double WLJ;
#if 1
void check (int *overlap, double *K, double *V)
{
  *overlap = 0;
}
#else
void check (int *overlap, double *K, double *V)
{
  /* *******************************************************************
   ** TESTS FOR PAIR OVERLAPS AND CALCULATES KINETIC ENERGY.        **
   *******************************************************************
   */
   int  i, j;
   double rxi, ryi, rzi, rxij, ryij, rzij, rijSq, sigSq, rij;
   double tol=OprogStatus.overlaptol;
   
   printf("overlaptol: %f\n", OprogStatus.overlaptol);
   *overlap = 0;
   *V = 0.0;
   *K = 0.0;   
   for (i = 0; i < Oparams.parnum-1; i++)
     {     
       rxi = rx[i];
       ryi = ry[i];
       rzi = rz[i];
       for (j = i+1; j < Oparams.parnum; j++)
	 {
	   rxij = rxi - rx[j];
	   ryij = ryi - ry[j];
	   rzij = rzi - rz[j];
	   rxij = rxij - L*rint(invL*rxij);
	   ryij = ryij - L*rint(invL*ryij);
#if !defined(MD_GRAVITY)
	   rzij = rzij - L*rint(invL*rzij);
#endif
	   rijSq = Sqr(rxij) + Sqr(ryij) + Sqr(rzij);
	   if (i < parnumA && j < parnumA)
	     sigSq = Sqr(Oparams.sigma[0][0]);
	   else if (i >= parnumA && j >= parnumA)
	     sigSq = Sqr(Oparams.sigma[1][1]);
	   else
	     sigSq = Sqr(Oparams.sigma[0][1]);
	   if ( rijSq < sigSq ) 
	     {
	       rij = sqrt(rijSq / sigSq);
	       if ( ( 1.0 - rij ) > tol ) 
		 {
		   *overlap = 1;
#if 1
    		   printf("#%d,#%d,rij/sigma =%f15.8 | ", i, j, rij);
    		   printf("(%f,%f,%f)-(%f,%f,%f)\n", rxi, ryi, rzi, 
			  rx[j], ry[j], rz[j]);
#endif

		 }
	     }
         }
     }
}
#endif
/* effettua la seguente moltiplicazione tra matrici:
 * 
 *                 | a 0 0 |               
 *  Trasposta(R) * | 0 b 0 | * R  dove R = {matrix unitaria relativa all'orientazione dell'ell.}
 *                 | 0 0 c |             
 * */
void tRDiagR(int i, double **M, double a, double b, double c, double **Ri)
{
  int na;
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  na = (i < Oparams.parnumA)?0:1;
  Di[0][0] = a;
  Di[1][1] = b;
  Di[2][2] = c;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  MD_DEBUG2(printf("a=%f b=%f c=%f Di=\n", a, b, c));
  MD_DEBUG2(print_matrixArr(Di));
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  MD_DEBUG2(print_matrixArr(Rtmp));
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
void RDiagtR(int i, double **M, double a, double b, double c, double **Ri)
{
  int na;
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  na = (i < Oparams.parnumA)?0:1;
  Di[0][0] = a;
  Di[1][1] = b;
  Di[2][2] = c;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  MD_DEBUG2(printf("a=%f b=%f c=%f Di=\n", a, b, c));
  MD_DEBUG2(print_matrixArr(Di));
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k2][k3];
	  }
      }
  MD_DEBUG2(print_matrixArr(Rtmp));
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k1][k3]*Rtmp[k3][k2];
	  }
      }
}

#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
void add_bond(int na, int n);
void remove_bond(int na, int n);
void bump (int i, int j, double* W, int bt)
{
  /*
   *******************************************************************
   ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
   **                                                               **
   ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
   ** THE ROUTINE ALSO COMPUTES COLLISIONAL VIRIAL W.               **
   *******************************************************************
   */
  double rxij, ryij, rzij, factor, invmi, invmj;
  double delpx, delpy, delpz;
  double mredl, ene;
  double sigSq, sigDeltaSq, intdistSq;
  double distSq;
  double vxij, vyij, vzij, b;

  if (i < parnumA && j < parnumA)
    {
      sigSq = Sqr(Oparams.sigma[0][0]);
      sigDeltaSq = Sqr(Oparams.sigma[0][0]+Oparams.delta[0][0]);
      mredl = Mred[0][0];
    }
  else if (i >= parnumA && j >= parnumA)
    {
      sigSq = Sqr(Oparams.sigma[1][1]);
      sigDeltaSq = Sqr(Oparams.sigma[1][1]+Oparams.delta[1][1]);
      mredl = Mred[1][1];
    }
  else
    {
      sigSq = Sqr(Oparams.sigma[0][1]);
      sigDeltaSq = Sqr(Oparams.sigma[0][1]+Oparams.delta[0][1]);
      mredl = Mred[0][1]; 
    }
  vxij = vx[i] - vx[j];
  vyij = vy[i] - vy[j];
  vzij = vz[i] - vz[j];
  /*printf("(i:%d,j:%d sigSq:%f\n", i, j, sigSq);*/
  /*printf("mredl: %f\n", mredl);*/
  rxij = rx[i] - rx[j];
  if (fabs (rxij) > L2)
    rxij = rxij - SignR(L, rxij);
  ryij = ry[i] - ry[j];
  if (fabs (ryij) > L2)
    ryij = ryij - SignR(L, ryij);
  rzij = rz[i] - rz[j];
  if (fabs (rzij) > L2)
    rzij = rzij - SignR(L, rzij);
  /* Nel caso di gravita' e' intuile implementare il TC-model di Luding
   * per evitare il collasso inelastico.
   * Gli urti in tale caso sono tutti elastici. */ 
  /* SQUARE WELL: modify here */
  distSq = Sqr(rxij)+Sqr(ryij)+Sqr(rzij);
  /*printf("distSq:%.20f\n",distSq);*/
  b = rxij * vxij + ryij * vyij + rzij * vzij;
  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;
  factor = 0.0;
  //printf("bump(%d,%d):%d distSq: %.20f b:%f\n", i, j, bt, distSq, b);
  switch (bt)
    {
    /* N.B.
     * Notare che Oparams.bheight è la profondità della buca ed 
     * è una quantità positiva!!*/
    case MD_CORE_BARRIER:
      factor = -2.0*b;
      factor *= mredl / sigSq;
      break;
    case MD_INOUT_BARRIER:
#ifdef MD_INFBARRIER
      factor = -2.0*b;
#elif defined(MD_SQWELL)
      if (Sqr(b) < 2.0*sigDeltaSq*Oparams.bheight/mredl)
	{
	  factor = -2.0*b;
	}
      else
	{
	  factor = -b + sqrt(Sqr(b) - 2.0*sigDeltaSq*Oparams.bheight/mredl);
	  remove_bond(i, j);
	  remove_bond(j, i);
	}
#endif
      factor *= mredl / sigDeltaSq;
#if 0
      if (fabs(distSq - sigDeltaSq)>1E-12)    
	printf("[bump]dist:%.20f\n",sqrt(distSq));
#endif
      break;
    case MD_OUTIN_BARRIER:
#ifdef MD_INFBARRIER
      factor = -2.0*b;
#elif defined(MD_SQWELL)
      add_bond(i, j);
      add_bond(j, i);
      factor = -b - sqrt(Sqr(b) + 2.0*sigDeltaSq*Oparams.bheight/mredl);
#endif
      factor *= mredl / sigDeltaSq;
      break;
    }
  
  delpx = factor * rxij;
  delpy = factor * ryij;
  delpz = factor * rzij;
#if 0
  ene= (Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])+
	Sqr(vx[j])+Sqr(vy[j])+Sqr(vz[j])); 
#endif
  vx[i] = vx[i] + delpx*invmi;
  vx[j] = vx[j] - delpx*invmj;
  vy[i] = vy[i] + delpy*invmi;
  vy[j] = vy[j] - delpy*invmj;
  vz[i] = vz[i] + delpz*invmi;
  vz[j] = vz[j] - delpz*invmj;
}
#else
void check_contact(int i, int j, double** Xa, double **Xb, double *rAC, double *rBC)
{
  int k1, k2;
  double f, g, del[3];
  f = g = -1; 
  for (k1=0; k1 < 3; k1++)
    del[0] = -L*rint((rAC[k1] - rBC[k1])/L);
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	f += rAC[k1]*Xa[k1][k2]*rAC[k2];
	g += rBC[k1]*Xb[k1][k2]*rBC[k2];
      } 
  printf("f(rC)=%.15G g(rC)=%.15G\n", f, g);
  //if (fabs(f) > 1E-5||fabs(g) > 1E-5)
   // exit(-1);
}
extern double **matrix(int n, int m);
extern void free_matrix(double **M, int n);
double calcDist(double t, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
     		double *vecgsup, int calcguess);
double calcDistNeg(double t, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
     		double *vecgsup, int calcguess);

void bump (int i, int j, double rCx, double rCy, double rCz, double* W)
{
  /*
   *******************************************************************
   ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
   **                                                               **
   ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
   ** THE ROUTINE ALSO COMPUTES COLLISIONAL VIRIAL W.               **
   *******************************************************************
   */
  double rxij, ryij, rzij, factor, invmi, invmj;
  double delpx, delpy, delpz, sigSq, wrx, wry, wrz, rACn[3], rBCn[3], rnI[3];
  double rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double norm[3], norm2[3];
  double modn, denom;
#ifndef MD_ASYM_ITENS
  double factorinvIa, factorinvIb;
#endif
  double shift[3];
  int na, a, b;
#if 0
  if (i < parnumA && j < parnumA)
    {
      sigSq = Sqr(Oparams.sigma[0][0]);
      //mredl = Mred[0][0];
    }
  else if (i >= parnumA && j >= parnumA)
    {
      sigSq = Sqr(Oparams.sigma[1][1]);
      //mredl = Mred[1][1];
    }
  else
    {
      sigSq = Sqr(Oparams.sigma[0][1]);
      mredl = Mred[0][1]; 
    }
#endif
  /*printf("(i:%d,j:%d sigSq:%f\n", i, j, sigSq);*/
  /*printf("mredl: %f\n", mredl);*/
  //MD_DEBUG(calc_energy("dentro bump1"));
  MD_DEBUG(printf("[bump] t=%f contact point: %f,%f,%f \n", Oparams.time, rxC, ryC, rzC));
  rAC[0] = rx[i] - rCx;
  rAC[1] = ry[i] - rCy;
  rAC[2] = rz[i] - rCz;
#if 1
  for (a=0; a < 3; a++)
    if (fabs (rAC[a]) > L2)
      rAC[a] -= SignR(L, rAC[a]);
#endif
  rBC[0] = rx[j] - rCx;
  rBC[1] = ry[j] - rCy;
  rBC[2] = rz[j] - rCz;
#if 1
#if MD_DEBUG(x) == x
    {
      double shift[3], r1[3], r2[3], alpha, vecgd[8], r12[3];
      shift[0] = L*rint((rx[i]-rx[j])/L);
      shift[1] = L*rint((ry[i]-ry[j])/L);
      shift[2] = L*rint((rz[i]-rz[j])/L);
      printf("shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]);
      printf("[bump] distance between %d-%d: %.15f\n", i, j, calcDistNeg(Oparams.time, i, j, shift, r1, r2, &alpha, vecgd, 1));
    }
#endif
  for (a=0; a < 3; a++)
    {
      MD_DEBUG(printf("P rBC[%d]:%.15f ", a, rBC[a]));
      if (fabs (rBC[a]) > L2)
	rBC[a] -= SignR(L, rBC[a]);
      MD_DEBUG(printf("D rBC[%d]:%.15f ", a, rBC[a]));
    }
  MD_DEBUG(printf("\n"));
#endif 
  /* calcola tensore d'inerzia e le matrici delle due quadriche */
  na = (i < Oparams.parnumA)?0:1;
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], R[i]);
#ifdef MD_ASYM_ITENS
  RDiagtR(i, Ia, ItensD[na][0], ItensD[na][1], ItensD[na][2], R[i]);
#else
  Ia = Oparams.I[na];
#endif
  na = (j < Oparams.parnumA)?0:1;
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], R[j]);
#ifdef MD_ASYM_ITENS
  RDiagtR(j, Ib, ItensD[na][0], ItensD[na][1], ItensD[na][2], R[j]);
#else
  Ib = Oparams.I[na];
#endif
  //MD_DEBUG(calc_energy("dentro bump2"));
  MD_DEBUG(check_contact(evIdA, evIdB, Xa, Xb, rAC, rBC));
  
  //MD_DEBUG(calc_energy("dentro bump3"));
  /* calcola le matrici inverse del tensore d'inerzia */
#ifdef MD_ASYM_ITENS
  InvMatrix(Ia, invIa, 3);
  InvMatrix(Ib, invIb, 3);
#else
  invIa = 1/Ia;
  invIb = 1/Ib;
  MD_DEBUG(printf("Ia=%f Ib=%f\n", Ia, Ib));
#endif
  for (a=0; a < 3; a++)
    {
      norm[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  norm[a] += -Xa[a][b]*rAC[b];
	}
    }
  modn = 0.0;
  for (a = 0; a < 3; a++)
    modn += Sqr(norm[a]);
  modn = sqrt(modn);
  for (a=0; a < 3; a++)
    norm[a] /= modn;
  MD_DEBUG(printf("CYL %f %f %f %f %f %f\n", rCx, rCy, rCz, norm[0], norm[1], norm[2]));
  /* calcola le velocità nel punto di contatto */
  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[0] = vx[i] + wrx;
  vCA[1] = vy[i] + wry;
  vCA[2] = vz[i] + wrz;
  vectProd(wx[j], wy[j], wz[j], -rBC[0], -rBC[1], -rBC[2], &wrx, &wry, &wrz);
  vCB[0] = vx[j] + wrx;
  vCB[1] = vy[j] + wry;
  vCB[2] = vz[j] + wrz;
 
  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;
  
  denom = invmi + invmj; 
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  MD_DEBUG(printf("[bump] before bump vc=%.15G\n", vc));
  if (vc < 0)// && fabs(vc) > 1E-10)
    {
      MD_DEBUG(printf("norm = (%f,%f,%f)\n", norm[0], norm[1],norm[2]));
      MD_DEBUG(printf("vel  = (%f,%f,%f)\n", vx[i], vy[i], vz[i]));
      MD_DEBUG(printf("i=%d r = (%f,%f,%f)\n", i, rx[i], ry[i], rz[i]));
      printf("[ERROR] maybe second collision has been wrongly predicted\n");
      printf("relative velocity (vc=%.15G) at contact point is negative! I ignore this event...\n", vc);
      return;
    }
  vectProd(rAC[0], rAC[1], rAC[2], norm[0], norm[1], norm[2], &rACn[0], &rACn[1], &rACn[2]);
#ifdef MD_ASYM_ITENS 
  for (a=0; a < 3; a++)
    {
      rnI[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  rnI[a] += invIa[a][b]*rACn[b]; 
	}
    }
  for (a = 0; a < 3; a++)
    denom += rnI[a]*rACn[a];
#else
  for (a = 0; a < 3; a++)
    denom += invIa*Sqr(rACn[a]);
#endif
  vectProd(rBC[0], rBC[1], rBC[2], norm[0], norm[1], norm[2], &rBCn[0], &rBCn[1], &rBCn[2]);
#ifdef MD_ASYM_ITENS  
  for (a=0; a < 3; a++)
    {
      rnI[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  rnI[a] += invIb[a][b]*rBCn[b]; 
	}
    }
  for (a = 0; a < 3; a++)
    denom += rnI[a]*rBCn[a];
#else
  for (a = 0; a < 3; a++)
    denom += invIb*Sqr(rBCn[a]);
#endif
#ifdef MD_GRAVITY
  factor =2.0*vc/denom;
  /* Dissipation */
#if 0
  /* se si vuole avere una dissipazione nell'urto questo va sistemato!!! */
  if (!((Oparams.time - lastcol[i] < OprogStatus.tc)||
  	(Oparams.time - lastcol[j] < OprogStatus.tc)))
    factor *= mredl*(1+Oparams.partDiss);
#endif
#else
  /* Nel caso di gravita' e' intuile implementare il TC-model di Luding
   * per evitare il collasso inelastico.
   * Gli urti in tale caso sono tutti elastici. */ 
  /* SQUARE WELL: modify here */
  factor =2.0*vc/denom;
  MD_DEBUG(printf("factor=%f denom=%f\n", factor, denom));
#endif
  delpx = - factor * norm[0];
  delpy = - factor * norm[1];
  delpz = - factor * norm[2];
#if 0
  ene= (Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])+
	    Sqr(vx[j])+Sqr(vy[j])+Sqr(vz[j])); 
#endif
  vx[i] = vx[i] + delpx*invmi;
  vx[j] = vx[j] - delpx*invmj;
  vy[i] = vy[i] + delpy*invmi;
  vy[j] = vy[j] - delpy*invmj;
  vz[i] = vz[i] + delpz*invmi;
  vz[j] = vz[j] - delpz*invmj;
  MD_DEBUG(printf("delp=(%f,%f,%f)\n", delpx, delpy, delpz));
#ifdef MD_ASYM_ITENS
  for (a=0; a < 3; a++)
    {
      wx[i] += factor*invIa[0][a]*rACn[a];
      wx[j] -= factor*invIb[0][a]*rBCn[a];
      wy[i] += factor*invIa[1][a]*rACn[a];
      wy[j] -= factor*invIb[1][a]*rBCn[a];
      wz[i] += factor*invIa[2][a]*rACn[a];
      wz[j] -= factor*invIb[2][a]*rBCn[a];
    }
#else
  factorinvIa = factor*invIa;
  factorinvIb = -factor*invIb;
  wx[i] += factorinvIa*rACn[0];
  wx[j] -= factorinvIb*rBCn[0];
  wy[i] += factorinvIa*rACn[1];
  wy[j] -= factorinvIb*rBCn[1];
  wz[i] += factorinvIa*rACn[2];
  wz[j] -= factorinvIb*rBCn[2];
#endif
#if MD_DEBUG(x) == x
  for (a=0; a < 3; a++)
    {
      norm2[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  norm2[a] += -Xb[a][b]*rBC[b];
	}
    }
  modn = 0.0;
  for (a = 0; a < 3; a++)
    modn += Sqr(norm2[a]);
  modn = sqrt(modn);
  for (a=0; a < 3; a++)
    norm2[a] /= modn;
  modn = 0;
  for (a=0; a < 3; a++)
    modn += norm2[a]*norm[a];
  printf("norm.norm2=%.15G\n", modn);
  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[0] = vx[i] + wrx;
  vCA[1] = vy[i] + wry;
  vCA[2] = vz[i] + wrz;
  vectProd(wx[j], wy[j], wz[j], -rBC[0], -rBC[1], -rBC[2], &wrx, &wry, &wrz);
  vCB[0] = vx[j] + wrx;
  vCB[1] = vy[j] + wry;
  vCB[2] = vz[j] + wrz;
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  printf("[bump] after bump vc=%.15G\n", vc); 
  //Oparams.time += 1.0; 
  //UpdateAtom(0);
  //UpdateAtom(1);
    {
      double shift[3], r1[3], r2[3], alpha, vecgd[8], r12[3];
      shift[0] = L*rint((rx[i]-rx[j])/L);
      shift[1] = L*rint((ry[i]-ry[j])/L);
      shift[2] = L*rint((rz[i]-rz[j])/L);
      printf("shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]);
      printf("[bump] distance between %d-%d: %.15f\n", i, j, calcDistNeg(Oparams.time, i, j, shift, r1, r2, &alpha, vecgd, 1));
    }
//exit(-1); 
#endif
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  i, vx[i],vy[i],vz[i], j, vx[j],vy[j],vz[j]));
  MD_DEBUG(printf("after bump %d-(%.10f,%.10f,%.10f) %d-(%.10f,%.10f,%.10f)\n", 
		  
		  i, wx[i],wy[i],wz[i], j, wx[j],wy[j],wz[j]));
/* TO CHECK: il viriale ha senso solo se non c'è la gravità */
#if 0
  *W = delpx * rxij + delpy * ryij + delpz * rzij;
#endif
}
#endif
#ifdef MD_GRAVITY
void calccmz(void)
{
  int i;
  double dd;
  rcmz = 0.0;
  
  for(i = 0; i < Oparams.parnumA; i++)
    {
      rcmz += rz[i];
    }
  rcmz *= Oparams.m[0];
  dd = 0.0;
  for(i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += rz[i];
    }
  dd *= Oparams.m[1];
  rcmz += dd;
  rcmz /= Mtot;

}
void calcRho(void)
{
  int jZ, jY, jX, n, npart, jZmax;
  double rhohcp, hhcp, Lz2, sig;
  double dia;
  Lz2 = Lz*0.5;
  /* hhcp è l'altezza delle particelle con diametro piu' piccolo 
   * se fossero close-packed */
  dia = (Oparams.sigma[0][0] < Oparams.sigma[1][1])?Oparams.sigma[0][0]:Oparams.sigma[1][1];
  rhohcp =0.7405*24/(4*pi*dia*dia*dia) ; 
  hhcp = Oparams.parnum/(rhohcp*Sqr(L));

  if (OprogStatus.rhobh <= 0)
    hhcp = Oparams.parnum/(rhohcp*Sqr(L))+OprogStatus.rhobh;
  else
    hhcp = (Oparams.parnum/(rhohcp*Sqr(L)))*OprogStatus.rhobh;
  /* Se rhobh > 0 allora l'altezza per il calcolo della densità è:
   * h_closepacking * rhobh */
  jZmax = (int) ((Lz / (OprogStatus.extraLz + Lz)) * cellsz)+1;  
  jZmax = (jZmax < cellsz)?jZmax:cellsz;
  MD_DEBUG2(printf("cellsz: %d iZmax: %d hhcp: %.15f\n", cellsz, jZmax, hhcp));
  npart = 0;
  for (jZ = 0; jZ < jZmax; jZ++)
    for (jX = 0; jX < cellsx; jX++)
      for (jY = 0; jY < cellsy; jY++)
	{
	  n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	  for (n = cellList[n]; n > -1; n = cellList[n]) 
	    {
	      if (n < parnumA)
		sig = Oparams.sigma[0][0];
	      else
		sig = Oparams.sigma[1][1];
	      if (rz[n] + Lz2 + sig*0.5 < hhcp)
		npart++;
	    }
	}
  rho = ((double)npart)/(L*L*hhcp);
}

void save_rho(void)
{
  FILE *f;
  f = fopenMPI(MD_HD_MIS "rho.dat", "a");
  fprintf(f, "%d %.15f\n", OprogStatus.numquench, rho);
  fclose(f);
}
void save_rzcm(void)
{
  FILE* f;
  f = fopenMPI(MD_HD_MIS "rcmz.dat", "a");
  fprintf(f, "%d %.15f\n", OprogStatus.numquench, rcmz);
  fclose(f);
}
#endif
#ifdef MD_GRAVITY
void calcObserv(void)
{
  int i;
  double dd1, dd;
  double dd2;
  K = 0.0;
  V = 0.0;
  Vz = 0.0;
  /* Bisogna considerare le velocità rispetto al centro di massa! */
  for (i=0; i < Oparams.parnumA; i++)
    {
      Vz += vz[i]; 
    }
  Vz *= Oparams.m[0];
  dd = 0.0;
  for (i=Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd += vz[i];
    }
  dd *= Oparams.m[1];
  Vz += dd;
  Vz /= (Oparams.parnumA*Oparams.m[0] + (Oparams.parnum-Oparams.parnumA)*Oparams.m[1]);
  calccmz();

  for (i = 0; i < Oparams.parnumA; i++)
    {
      K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
      V += rz[i];
    }
  K *= Oparams.m[0] * 0.5;
  V *= mgA;
  dd1 = 0.0;
  dd2 = 0.0;
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      dd1 += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
      dd2 += rz[i];
    }
  dd1 *= Oparams.m[1] * 0.5;
  dd2 *= mgB;
  K += dd1;
  V += dd2;
}
#else
void calcObserv(void)
{
  /* DESCRIPTION:
     This mesuring functions calculates the Translational Diffusion 
     coefficent */
  FILE *f;
  double Drx, Dry, Drz;
  int i;
#ifdef MPI
  int equilib, sumeq;
#endif
  
  DrSqTot = 0.0;
  K = 0.0;
  for (i = 0; i < Oparams.parnumA; i++)
    {
      K += Oparams.m[0]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
    }
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      K += Oparams.m[1]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
    }
  K *= 0.5;
  for(i=0; i < Oparams.parnumA; i++)
    {
      Drx = rx[i] - OprogStatus.rxCMi[i] + L*OprogStatus.DR[i][0]; 
      Dry = ry[i] - OprogStatus.ryCMi[i] + L*OprogStatus.DR[i][1];
      Drz = rz[i] - OprogStatus.rzCMi[i] + L*OprogStatus.DR[i][2];
      DrSqTot = DrSqTot + Sqr(Drx) + Sqr(Dry) + Sqr(Drz);
   }
  /* NOTE: The first Dtrans(first simulation step) is not meaningful, 
     because DrSq is zero! */
  if (Oparams.time>0)
    Dtrans = DrSqTot / ( 6.0 * ((double) Oparams.time) *
		   	 ((double) Oparams.parnumA ) );   
  else 
    Dtrans = 0;

  DrSqTot /= ((double) Oparams.parnumA);
  if (OprogStatus.eqlevel > 0.0)
    {
      if (!OprogStatus.equilibrated)
	OprogStatus.equilibrated = (DrSqTot>OprogStatus.eqlevel*(2.0*maxax[0])?1:0);
#ifdef MPI
      equilib = OprogStatus.equilibrated;
      MPI_Allgather(&equilib, 1, MPI_INT, 
		    equilibrated, 1, MPI_INT, MPI_COMM_WORLD);
      sumeq = 0;
      for (i=0; i < numOfProcs; i++)
	sumeq += equilibrated[i];
      /* se sumeq = numOfProcs vuol dire che tutti i processi sono
       * equilibrati quindi la simulazione può terminare */
      if (sumeq == numOfProcs)
	{
	  mdPrintf(ALL,"All systems reached equilibrium, simulation completed", 
		   NULL);
	  ENDSIM = 1;
	}
#else
      mdPrintf(ALL, "All systems reached equilibrium, simulation completed",
	   NULL);
      ENDSIM = OprogStatus.equilibrated;
#endif
    }
  if (Oparams.time>0)
    {
      f = fopenMPI(MD_HD_MIS "D.dat", "a");
      fprintf(f, "%.15f %.15f\n", Oparams.time,  Dtrans);
      fclose(f);
    }
}
#endif
extern double *treeTime;
void adjust_norm(double **R)
{
  int k1, k2; 
  double n[3];
  for (k1 = 0; k1 < 3; k1++)
    {
      n[k1]=0;
      for(k2 = 0; k2 < 3; k2++)
	n[k1] += Sqr(R[k2][k1]);
      n[k1] = sqrt(n[k1]);
      if (fabs((n[k1])-1.0)>1E-10)
	{
	  MD_DEBUG(printf("Adjusting norm of orientations time=%.15f\n", Oparams.time));
	  MD_DEBUG(printf("delta = %.15f\n", fabs(n[k1]-1.0)));
	  for(k2 = 0; k2 < 3; k2++)
	    R[k2][k1] /= n[k1];
	}
    }
  
}

void UpdateAtom(int i)
{
  double ti;
  double wSq, w, sinw, cosw;
  double Omega[3][3], OmegaSq[3][3], Rtmp[3][3], M[3][3];
  double w1[3],w2[3];
  int k1, k2, k3;
  ti = Oparams.time - atomTime[i];
  
  rx[i] += vx[i]*ti;
  ry[i] += vy[i]*ti;
#if defined(MD_GRAVITY)
  rz[i] += vz[i]*ti - g2*Sqr(ti);
  vz[i] += -Oparams.ggrav*ti;
#else
  rz[i] += vz[i]*ti;
#endif
  /* ...and now orientations */
  wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
  w = sqrt(wSq);
  if (w != 0.0) 
    {
#if 0
      if (fabs(w*ti) < 1E-8)
	{
	  sinw = ti*(1-Sqr(w*ti)/6.0);	  
	  cosw = Sqr(ti)*(0.5 - Sqr(w*ti)/24.0);
	}
#endif
      sinw = sin(w*ti)/w;
      cosw = (1.0 - cos(w*ti))/wSq;
      Omega[0][0] = 0;
      Omega[0][1] = -wz[i];
      Omega[0][2] = wy[i];
      Omega[1][0] = wz[i];
      Omega[1][1] = 0;
      Omega[1][2] = -wx[i];
      Omega[2][0] = -wy[i];
      Omega[2][1] = wx[i];
      Omega[2][2] = 0;
      OmegaSq[0][0] = -Sqr(wy[i]) - Sqr(wz[i]);
      OmegaSq[0][1] = wx[i]*wy[i];
      OmegaSq[0][2] = wx[i]*wz[i];
      OmegaSq[1][0] = wx[i]*wy[i];
      OmegaSq[1][1] = -Sqr(wx[i]) - Sqr(wz[i]);
      OmegaSq[1][2] = wy[i]*wz[i];
      OmegaSq[2][0] = wx[i]*wz[i];
      OmegaSq[2][1] = wy[i]*wz[i];
      OmegaSq[2][2] = -Sqr(wx[i]) - Sqr(wy[i]);
      
      for (k1 = 0; k1 < 3; k1++)
	{
	  
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      Omega[k1][k2] = -Omega[k1][k2];
	      Rtmp[k1][k2] = R[i][k1][k2];
	      M[k1][k2] = sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
#if 0
	      if (k1==k2)
		M[k1][k1] += 1.0;
#endif
	    }
	}
#if MD_DEBUG(x)==x
      w1[0] = wx[i];
      w1[1] = wy[i];
      w1[2] = wz[i];
      for (k1 = 0; k1 < 3; k1++)
	{
	  w2[k1] = 0;
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      w2[k1] += (sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2])*w1[k2];   
	    }	 
	   
	}
      if (fabs(w2[0])>1E-12 || fabs(w2[1])>1E-12
	  || fabs(w2[2])>1E-12)
	printf("ti=%.15G w=%.15G cosw=%.15G sinw=%.15G w = (%.15f,%.15f,%.15f) Exp(Omega t) w = (%.15f,%.15f,%.15f)\n",
	       ti, w, cosw, sinw, w1[0],w1[1],w1[2],w2[0],w2[1],w2[2]);	
#endif
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
#if 0
	    R[i][k1][k2] = 0.0;
#endif
	    for (k3 = 0; k3 < 3; k3++)
	      R[i][k1][k2] += M[k1][k3]*Rtmp[k3][k2];
	  }
      adjust_norm(R[i]);
#if 0
      if (rz[i]+Lz*0.5-Oparams.sigma/2.0 < 0. && OprogStatus.quenchend > 0.0)
	{
	  int no;
	  printf("rz[i](t-ti):%.15f rz[i]:%.15f ti:%.15f", 
		 rz[i] - vz[i]*ti + g2*Sqr(ti), rz[i],
		 ti);
	  printf("vz[i]:%.15f\n", vz[i]);
	  printf("i=%d SOTTO MURO PT rz+Lz*0.5=%.30f\n", i, rz[i]+Lz*0.5-Oparams.sigma/2.0);
	  printf("*******************************************************************\n");
	  if (vz[i] < 0 && evIdA != i && evIdB != ATOM_LIMIT +4)
	    {
 
	      exit(-1);
	    }
	}
#endif
    }
  atomTime[i] = Oparams.time;
}
void UpdateSystem(void)
{
  int i;
  /* porta tutte le particelle allo stesso tempo */
  for (i=0; i < Oparams.parnum; i++)
    {
      UpdateAtom(i);
#if 0
      if (ry[i] > L2 || ry[i] < -L2)
	{
	  printf("Porca Troia!!!\n");
	  exit(-1);
	}
#endif
    }
}

#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
void remove_bond(int na, int n)
{
  int i, nb, ii;
  nb = numbonds[na];
  if (!nb)
    return;
  ii = 0;
  memcpy(bondscache, bonds[na], sizeof(int)*numbonds[na]);
  for (i = 0; i < nb; i++)
    if (bondscache[i] != n)
      {
	bonds[na][ii++] = bondscache[i];
      } 
    else
      numbonds[na]--;
  if (nb==numbonds[na])
    printf("nessun bond rimosso fra %d,%d\n", n, na);
}
int bound(int na, int n);

void add_bond(int na, int n)
{
  if (bound(na, n))
    {
      printf("il bond %d,%d eiste già!\n", na, n);
      return;
    }
  bonds[na][numbonds[na]] = n;
  numbonds[na]++;
}

int bound(int na, int n)
{
  int i;
  for (i = 0; i < numbonds[na]; i++)
    if (bonds[na][i] == n)
      return 1;
  return 0;
}
#endif
UpdateOrient(int i, double ti, double **Ro, double Omega[3][3])
{ 
  double wSq, w, OmegaSq[3][3], M[3][3];
  double sinw, cosw;
  int k1, k2, k3;
  wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
  w = sqrt(wSq);
  if (w != 0.0) 
    {
#if 0
      if (fabs(w*ti) < 1E-8)
	{
	  sinw = ti*(1-Sqr(w*ti)/6.0);	  
	  cosw = Sqr(ti)*(0.5 - Sqr(w*ti)/24.0);
	}
      else 
	{
	  sinw = sin(w*ti)/w;
    	  cosw = (1.0 - cos(w*ti))/wSq;
	}
#endif
      sinw = sin(w*ti)/w;
      cosw = (1.0 - cos(w*ti))/wSq;
      Omega[0][0] = 0;
      Omega[0][1] = -wz[i];
      Omega[0][2] = wy[i];
      Omega[1][0] = wz[i];
      Omega[1][1] = 0;
      Omega[1][2] = -wx[i];
      Omega[2][0] = -wy[i];
      Omega[2][1] = wx[i];
      Omega[2][2] = 0;
      OmegaSq[0][0] = -Sqr(wy[i]) - Sqr(wz[i]);
      OmegaSq[0][1] = wx[i]*wy[i];
      OmegaSq[0][2] = wx[i]*wz[i];
      OmegaSq[1][0] = wx[i]*wy[i];
      OmegaSq[1][1] = -Sqr(wx[i]) - Sqr(wz[i]);
      OmegaSq[1][2] = wy[i]*wz[i];
      OmegaSq[2][0] = wx[i]*wz[i];
      OmegaSq[2][1] = wy[i]*wz[i];
      OmegaSq[2][2] = -Sqr(wx[i]) - Sqr(wy[i]);
 
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      Omega[k1][k2] = -Omega[k1][k2];
	      M[k1][k2] = sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
#if 0
	      if (k1==k2)
	      	M[k1][k1] += 1.0;
#endif
	    }
	}
      
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    Ro[k1][k2] = R[i][k1][k2];
	    for (k3 = 0; k3 < 3; k3++)
	      Ro[k1][k2] += M[k1][k3]*R[i][k3][k2];
	  }
      adjust_norm(Ro);
    }
  else
    {
      Omega[0][0] = 0;
      Omega[0][1] = 0;
      Omega[0][2] = 0;
      Omega[1][0] = 0;
      Omega[1][1] = 0;
      Omega[1][2] = 0;
      Omega[2][0] = 0;
      Omega[2][1] = 0;
      Omega[2][2] = 0;
      for (k1 = 0; k1 < 3; k1++)
	for (k2 = 0; k2 < 3; k2++)
	  {
	    Ro[k1][k2] = R[i][k1][k2];
	  }
    }
}

#ifndef MD_APPROX_JACOB
extern double **matrix(int n, int m);
extern void free_matrix(double **M, int n);
void calcFxtFt(double x[3], double **X,
	       double D[3][3], double Omega[3][3], double **R, 
	       double pos[3], double vel[3], double gradf[3],
	       double Fxt[3], double *Ft)
{
  double tOmegaD[3][3], DOmega[3][3];
  double sumDOmega[3][3], Mtmp[3][3], DtX[3][3], dx[3];
  int k1, k2, k3;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  DOmega[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (D[k1][k3] == 0.0 || Omega[k3][k2] == 0.0)
		continue;
	      DOmega[k1][k2] += D[k1][k3]*Omega[k3][k2];
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  tOmegaD[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    {
	      if (D[k3][k2] == 0.0 || Omega[k3][k1] == 0.0)
		continue;
	      tOmegaD[k1][k2] += Omega[k3][k1]*D[k3][k2]; 
	    }
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      sumDOmega[k1][k2] = tOmegaD[k1][k2] + DOmega[k1][k2];

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  Mtmp[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    Mtmp[k1][k2] += sumDOmega[k1][k3]*R[k3][k2]; 
	}
    }
   for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  DtX[k1][k2] = 0;
	  for (k3 = 0; k3 < 3; k3++)
	    DtX[k1][k2] += R[k3][k1]*Mtmp[k3][k2]; 
	}
    }
   for (k1 = 0; k1 < 3; k1++)
     {
       Fxt[k1] = 0;
       for (k2 = 0; k2 < 3; k2++)
	 Fxt[k1] += DtX[k1][k2]*(x[k2]-pos[k2]) - Sqr(x[3])*X[k1][k2]*vel[k2]; 
       Fxt[k1] *= 2.0;
     } 
   *Ft = 0;
   for (k1 = 0; k1 < 3; k1++)
     dx[k1] = x[k1]-pos[k1];
   for (k1 = 0; k1 < 3; k1++)
     {
       for (k2 = 0; k2 < 3; k2++)
	 {
	   *Ft += -vel[k1]*X[k1][k2]*dx[k2]+dx[k1]*DtX[k1][k2]*dx[k2]-dx[k1]*X[k1][k2]*vel[k2];
	 }
     }
}
#undef MD_GLOBALNR
#undef MD_GLOBALNR2
double rA[3], rB[3];
void fdjacGuess(int n, double x[], double fvec[], double **df, 
	   void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  double fx[3], gx[3], tmp;
  double Fxt[3], Gxt[3], Ft, Gt;
  int k1, k2, k3;

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*(Xa[k1][k2] + Sqr(x[3])*Xb[k1][k2]);
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]);
	}
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][3] = 2.0*x[3]*gx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1] - gx[k1];
      printf("[fdjacGuess] fx[%d]:%.15G gx=%.15G\n",k1 , fx[k1],gx[k1]);
    } 

  df[3][3] = 0.0;
#ifndef MD_GLOBALNR2
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
  fvec[3] = 0.0;
  tmp = 0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      tmp += (x[k1]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0 - (0.5*tmp-1.0);
#endif
}

/* funzione che calcola lo Jacobiano */
void fdjac(int n, double x[], double fvec[], double **df, 
	   void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  /* N.B. QUESTA ROUTINE VA OTTIMIZZATA! ad es. calcolando una sola volta i gradienti di A e B...*/
  int na; 
  double  rA[3], rB[3], ti, vA[3], vB[3], OmegaA[3][3], OmegaB[3][3];
  double DA[3][3], DB[3][3], fx[3], gx[3];
  double Fxt[3], Gxt[3], Ft, Gt;
  int k1, k2, k3;
  ti = x[4] - atomTime[iA];
  rA[0] = rx[iA] + vx[iA]*ti;
  rA[1] = ry[iA] + vy[iA]*ti;
  rA[2] = rz[iA] + vz[iA]*ti;
  vA[0] = vx[iA];
  vA[1] = vy[iA];
  vA[2] = vz[iA];
  /* ...and now orientations */
  UpdateOrient(iA, ti, RA, OmegaA);
  MD_DEBUG2(printf("i=%d ti=%f", iA, ti));
  MD_DEBUG2(print_matrix(RA, 3));
  na = (iA < Oparams.parnumA)?0:1;
  tRDiagR(iA, Xa, invaSq[na], invbSq[na], invcSq[na], RA);
  MD_DEBUG2(printf("invabc: (%f,%f,%f)\n", invaSq[na], invbSq[na], invcSq[na]));
  MD_DEBUG2(print_matrix(Xa, 3));
  DA[0][1] = DA[0][2] = DA[1][0] = DA[1][2] = DA[2][0] = DA[2][1] = 0.0;
  DA[0][0] = invaSq[na];
  DA[1][1] = invbSq[na];
  DA[2][2] = invcSq[na];
  ti = x[4] - atomTime[iB];
  rB[0] = rx[iB] + vx[iB]*ti + shift[0];
  rB[1] = ry[iB] + vy[iB]*ti + shift[1];
  rB[2] = rz[iB] + vz[iB]*ti + shift[2];
  vB[0] = vx[iB];
  vB[1] = vy[iB];
  vB[2] = vz[iB];
  UpdateOrient(iB, ti, RB, OmegaB);
  na = (iB < Oparams.parnumA)?0:1;
  tRDiagR(iB, Xb, invaSq[na], invbSq[na], invcSq[na], RB);
  DB[0][1] = DB[0][2] = DB[1][0] = DB[1][2] = DB[2][0] = DB[2][1] = 0.0;
  DB[0][0] = invaSq[na];
  DB[1][1] = invbSq[na];
  DB[2][2] = invcSq[na];
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*(Xa[k1][k2] + Sqr(x[3])*Xb[k1][k2]);
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]);
	}
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[3][k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[3][k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]); 
#endif
      df[3][k1] = fx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[4][k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[4][k1] += 2.0*Xb[k1][k2]*(x[k2]-rB[k2]); 
#endif
      df[4][k1] = gx[k1];
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      df[k1][3] = 0;
      for (k2 = 0; k2 < 3; k2++)
	df[k1][3] += 4.0*x[3]*Xb[k1][k2]*(x[k2]-rB[k2]); 
#endif
      df[k1][3] = 2.0*x[3]*gx[k1];
    } 
  df[3][3] = 0.0;
  df[4][3] = 0.0;
  calcFxtFt(x, Xa, DA, OmegaA, RA, rA, vA, fx, Fxt, &Ft);
  calcFxtFt(x, Xb, DB, OmegaB, RB, rB, vB, gx, Gxt, &Gt);
  for (k1 = 0; k1 < 3; k1++)
    {
      //df[k1][4] = 0;
      //for (k2 = 0; k2 < 3; k2++)
      df[k1][4] = Fxt[k1]+Sqr(x[3])*Gxt[k1]; 
    } 
 df[3][4] = Ft;
 df[4][4] = Gt;
#ifndef MD_GLOBALNR
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
 MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
#endif
}
#endif
void upd2tGuess(int i, int j, double shift[3], double tGuess)
{
  double ti;
  int na;
  double Omega[3][3];
  ti = tGuess - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
  UpdateOrient(i, ti, Rt, Omega);
  na = (i < Oparams.parnumA)?0:1;
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = tGuess - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, Rt, Omega);
  na = (j < Oparams.parnumA)?0:1;
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);

}
void funcs2beZeroedGuess(int n, double x[], double fvecG[], int i, int j, double shift[3])
{
  int na, k1, k2; 
  double fx[3], gx[3], tmp;
  /* x = (r, alpha, t) */ 
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvecG[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
#if 0
  fvec[3] = -1.0;
  fvec[4] = -1.0;
#endif
#if 0
  MD_DEBUG(printf("fx+Sqr(alpha)*gx=(%f,%f,%f) fx=(%f,%f,%f) gx=(%f,%f,%f)\n", fvec[0], fvec[1], fvec[2],
	 fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
#endif
  fvecG[3] = 0.0;
  tmp = 0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvecG[3] += (x[k1]-rA[k1])*fx[k1];
      tmp += (x[k1]-rB[k1])*gx[k1];
    }
  fvecG[3] = 0.5*fvecG[3]-1.0 - (0.5*tmp-1.0);
}
void funcs2beZeroed(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int na, k1, k2; 
  double  rA[3], rB[3], ti;
  double fx[3], gx[3];
  double Omega[3][3];
  /* x = (r, alpha, t) */ 
  
  ti = x[4] - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  /* ...and now orientations */
  UpdateOrient(i, ti, Rt, Omega);
  na = (i < Oparams.parnumA)?0:1;
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = x[4] - atomTime[j];
  MD_DEBUG(printf("x[4]:%.15f atomTime[%d]:%.15f\n",x[4], j, atomTime[j]));
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, Rt, Omega);
  na = (j < Oparams.parnumA)?0:1;
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
      MD_DEBUG(printf("[FUNC2BEZ]x[%d]:%.15f rA[%d]:%f fx:%.15f\n", k1, x[k1], k1, rA[k1],fx[k1]));

    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2] - rB[k2]);
      MD_DEBUG(printf("[FUNC2BEZ]x[%d]:%.15f rB[%d]:%f gx:%.15f\n", k1, x[k1], k1, rB[k1],gx[k1]));
    }

  MD_DEBUG(print_matrix(Xb,3));
  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      fvec[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fvec[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]) + 2.0*Sqr(x[3])*Xb[k1][k2]*(x[k2] - rB[k2]);
#endif
      fvec[k1] = fx[k1] + Sqr(x[3])*gx[k1];
    }
#if 0
  fvec[3] = -1.0;
  fvec[4] = -1.0;
#endif
#if 0
  MD_DEBUG(printf("fx+Sqr(alpha)*gx=(%f,%f,%f) fx=(%f,%f,%f) gx=(%f,%f,%f)\n", fvec[0], fvec[1], fvec[2],
	 fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
#endif
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
#if 0
      for (k2 = 0; k2 < 3; k2++)
	{
	  fvec[3] += (x[k1]-rA[k1])*Xa[k1][k2]*(x[k2]-rA[k2]);
	  fvec[4] += (x[k1]-rB[k1])*Xb[k1][k2]*(x[k2]-rB[k2]);
	}
#endif
#if 1
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1]-rB[k1])*gx[k1];
#endif
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;
  MD_DEBUG(printf("F2BZ fvec (%.12f,%.12f,%.12f,%.12f,%.13f)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]));
}

double tdist;
double rA[3], rB[3];
#undef MD_GLOBALNRD
void fdjacDistNeg(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  int na; 
  double fx[3], gx[3];
  int k1, k2, k3;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = 2.0*Sqr(x[6])*Xb[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3]-rB[k2]);
	}
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  for (k1 = 0; k1 < 5; k1++)
    {
      df[3][k1+3] = 0;
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = 0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1+3] = gx[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = 2.0*x[6]*gx[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1 + 2.0*x[7]*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 2.0*x[7]*Xa[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = fx[k1];
#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
  //MD_DEBUG(printf("F2BZdistNeg fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
}

/* funzione che calcola lo Jacobiano */
void fdjacDist(int n, double x[], double fvec[], double **df, 
    	       void (*vecfunc)(int, double [], double []), int iA, int iB, double shift[3])
{
  int na; 
  double fx[3], gx[3];
  int k1, k2, k3;
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
       	{
	  df[k1][k2] = 2.0*Xa[k1][k2];
	  df[k1][k2+3] = 2.0*Sqr(x[6])*Xb[k1][k2];
	}
    }
  /* calc fx e gx */
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  fx[k1] += 2.0*Xa[k1][k2]*(x[k2]-rA[k2]);
	  gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3]-rB[k2]);
	}
    } 

  for (k1 = 0; k1 < 3; k1++)
    {
      df[3][k1] = fx[k1];
    } 
  for (k1 = 0; k1 < 5; k1++)
    {
      df[3][k1+3] = 0;
    } 
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1] = 0;
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      df[4][k1+3] = gx[k1];
    } 
  df[4][6] = df[4][7] = 0;

  for (k1 = 0; k1 < 3; k1++)
    {
      df[k1][6] = 2.0*x[6]*gx[k1];
      df[k1][7] = 0.0;
    } 

  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2] = 1 + 2.0*Sqr(x[7])*Xa[k1][k2];
	  else 
	    df[k1+5][k2] = 2.0*Sqr(x[7])*Xa[k1][k2];
	}
    }
  for (k1=0; k1<3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  if (k1==k2)
	    df[k1+5][k2+3] = -1;
	  else 
	    df[k1+5][k2+3] = 0;
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][6] = 0;
  for (k1 = 0; k1 < 3; k1++)
    df[k1+5][7] = 2.0*x[7]*fx[k1];
#ifndef MD_GLOBALNRD
 /* and now evaluate fvec */
 for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
 fvec[3] = 0.0;
 fvec[4] = 0.0;
 for (k1 = 0; k1 < 3; k1++)
   {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
   }
 fvec[3] = 0.5*fvec[3]-1.0;
 fvec[4] = 0.5*fvec[4]-1.0;
  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*Sqr(x[7]); 
  //MD_DEBUG(printf("F2BZdist fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G,%.12G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
#endif
}
void funcs2beZeroedDistNeg(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3];
  /* x = (r, alpha, t) */ 
  
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*x[7]; 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}

void funcs2beZeroedDist(int n, double x[], double fvec[], int i, int j, double shift[3])
{
  int k1, k2; 
  double fx[3], gx[3];
  /* x = (r, alpha, t) */ 
  
#if 0
  printf("Xa=\n");
  print_matrix(Xa, 3);
  printf("Xb=\n");
  print_matrix(Xb, 3);
#endif
  
  for (k1 = 0; k1 < 3; k1++)
    {
      fx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	fx[k1] += 2.0*Xa[k1][k2]*(x[k2] - rA[k2]);
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      gx[k1] = 0;
      for (k2 = 0; k2 < 3; k2++)
	gx[k1] += 2.0*Xb[k1][k2]*(x[k2+3] - rB[k2]);
    }

   for (k1 = 0; k1 < 3; k1++)
    {
      fvec[k1] = fx[k1] + Sqr(x[6])*gx[k1];
    }
  fvec[3] = 0.0;
  fvec[4] = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    {
      fvec[3] += (x[k1]-rA[k1])*fx[k1];
      fvec[4] += (x[k1+3]-rB[k1])*gx[k1];
    }
  fvec[3] = 0.5*fvec[3]-1.0;
  fvec[4] = 0.5*fvec[4]-1.0;

  /* N.B. beta=x[7] non è al quadrato poichè in questo modo la distanza puo' 
   * essere anche negativa! */
  for (k1=0; k1 < 3; k1++)
    fvec[k1+5] = x[k1] - x[k1+3] + fx[k1]*Sqr(x[7]); 
#if 0
  MD_DEBUG(printf("fx: (%f,%f,%f) gx (%f,%f,%f)\n", fx[0], fx[1], fx[2], gx[0], gx[1], gx[2]));
  MD_DEBUG(printf("fvec (%.12G,%.12G,%.12G,%.12G,%.12G,%.15G,%.15G,%.15G)\n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4],fvec[5],fvec[6],fvec[7]));
  MD_DEBUG(printf("x (%f,%f,%f,%f,%f,%f,%f)\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]));
#endif
}
void calc_intersec(double *rB, double *rA, double **Xa, double* rI)
{
  double A, B, C, D, tt;
  double rBA[3], XarA[3];
  int k1, k2, k3;
  for (k1=0; k1 < 3; k1++)
    rBA[k1] = rB[k1] - rA[k1];
  MD_DEBUG(printf("rBA=(%f,%f,%f)\n", rBA[0], rBA[1], rBA[2])); 
  MD_DEBUG(printf("rB= (%f,%f,%f)\n", rB[0], rB[1], rB[2]));
  A = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      A += rBA[k1]*Xa[k1][k2]*rBA[k2];
 
#if 0
  for (k1 = 0; k1 < 3; k1++)
    {   
      XarA[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	XarA[k1] += Xa[k1][k2]*rA[k2];
    }
  B = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    B += rBA[k1]*XarA[k1];
  B *= 2.0;
  C = 0.0;
  for (k1 = 0; k1 < 3; k1++)
    C +=  rA[k1]*XarA[k1];
  C = C - 1.0;
  D = Sqr(B) - 4.0*A*C;
  if (D < 0 || A == 0)
    {
      printf("[calc_intersec] Serious problem guessing distance, aborting...\n");
      printf("tt = %f D=%f A=%f B=%f C=%f\n", tt, D, A, B, C);
      printf("distance: %f\n", sqrt(Sqr(rBA[0])+Sqr(rBA[1])+Sqr(rBA[2])));
      print_matrix(Xa,3);
      print_matrix(Xb,3);
      exit(-1);
    }
  tt = (-B + sqrt(D))/(2.0*A); 
#endif
  if (A <= 0)
    {
      printf("[calc_intersec] Serious problem guessing distance, aborting...\n");
      printf("tt = %f D=%f A=%f B=%f C=%f\n", tt, D, A, B, C);
      printf("distance: %f\n", sqrt(Sqr(rBA[0])+Sqr(rBA[1])+Sqr(rBA[2])));
      print_matrix(Xa,3);
      print_matrix(Xb,3);
      exit(-1);
    }
  tt = sqrt(1 / A); 
  for (k1 = 0; k1 < 3; k1++)
    {
      rI[k1] = rA[k1] + tt*rBA[k1];  
    }
}
void calc_grad(double *rC, double *rA, double **Xa, double *grad)
{
  int k1, k2;
  for (k1 = 0; k1 < 3; k1++)
    {
      grad[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	grad[k1] += 2.0*Xa[k1][k2]*(rC[k2]-rA[k2]); 
    }
}
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
double calcDistNeg(double t, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
     		double *vecgsup, int calcguess)
{
  double vecg[8], rC[3], rD[3], rDC[3], r12[3];
  double ti, segno;
  int retcheck;
  double Omega[3][3], nf, ng, gradf[3], gradg[3];
  int k1, k2, na;
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  ti = t - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
  UpdateOrient(i, ti, Rt, Omega);
  na = (i < Oparams.parnumA)?0:1;
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = t - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, Rt, Omega);
  na = (j < Oparams.parnumA)?0:1;
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
  if (calcguess)
    {
      calc_intersec(rB, rA, Xa, rC);
      calc_intersec(rA, rB, Xb, rD);
      MD_DEBUG(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		      rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
      calc_grad(rC, rA, Xa, gradf);
      calc_grad(rD, rB, Xb, gradg);
      MD_DEBUG(printf("gradf=(%f,%f,%f) gradg=(%f,%f,%f)\n",
		      gradf[0], gradf[1], gradf[2], gradg[0], gradg[1], gradg[2]));
      nf = calc_norm(gradf);
      ng = calc_norm(gradg);
      
      vecg[6] = sqrt(nf/ng);
      for (k1=0; k1 < 3; k1++)
	{
	  vecg[k1] = rC[k1];
	  vecg[k1+3] = rD[k1];
	  rDC[k1] = rD[k1] - rC[k1];
	}
      vecg[7] = sqrt(calc_norm(rDC)/nf);  
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  MD_DEBUG(printf("alpha: %f beta: %f\n", vecg[6], vecg[7]));
  newtDistNeg(vecg, 8, &retcheck, funcs2beZeroedDistNeg, i, j, shift); 
  if (retcheck != 0)
    {
      printf("I couldn't calculate distance between %d and %d\n, exiting....\n", i, j);
      Oparams.time = t;
      store_bump(i, j);
      exit(-1);
    }
  for (k1 = 0; k1 < 8; k1++)
    {
      vecgsup[k1] = vecg[k1]; 
    }  
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      r2[k1] = vecg[k1+3];
      r12[k1] = r1[k1] - r2[k1];
    }
  *alpha = vecg[6];
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (r2[k1]-rA[k1])*Xa[k1][k2]*(r2[k2]-rA[k2]); 
#if 1
  if (segno > 0)
    return calc_norm(r12);
  else
    return -calc_norm(r12);
#else
#if MD_DEBUG(x) == x
  for (k1 = 0; k1 < 3; k1++)
    rDC[k1] = r1[k1] - r2[k1];
    {
      FILE *f;
      f =fopen("distNeg.dat","a");
      fprintf(f,"%.15f %.15f\n", Oparams.time, (segno>0)?calc_norm(rDC):-calc_norm(rDC));  
      fclose(f);
    }
  return calc_norm(rDC);
#endif
#endif
}
double calcDist(double t, int i, int j, double shift[3], double *r1, double *r2, double *alpha,
     		double *vecgsup, int calcguess)
{
  double vecg[8], rC[3], rD[3], rDC[3], r12[3];
  double ti, segno;
  int retcheck;
  double Omega[3][3], nf, ng, gradf[3], gradg[3];
  int k1, k2, na;
  MD_DEBUG(printf("t=%f tai=%f taj=%f i=%d j=%d\n", t, t-atomTime[i],t-atomTime[j],i,j));
  ti = t - atomTime[i];
  rA[0] = rx[i] + vx[i]*ti;
  rA[1] = ry[i] + vy[i]*ti;
  rA[2] = rz[i] + vz[i]*ti;
  MD_DEBUG(printf("rA (%f,%f,%f)\n", rA[0], rA[1], rA[2]));
  /* ...and now orientations */
  UpdateOrient(i, ti, Rt, Omega);
  na = (i < Oparams.parnumA)?0:1;
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], Rt);

  ti = t - atomTime[j];
  rB[0] = rx[j] + vx[j]*ti + shift[0];
  rB[1] = ry[j] + vy[j]*ti + shift[1];
  rB[2] = rz[j] + vz[j]*ti + shift[2];
  UpdateOrient(j, ti, Rt, Omega);
  na = (j < Oparams.parnumA)?0:1;
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], Rt);
  if (calcguess)
    {
      calc_intersec(rB, rA, Xa, rC);
      calc_intersec(rA, rB, Xb, rD);
      MD_DEBUG(printf("rC=(%f,%f,%f) rD=(%f,%f,%f)\n",
		      rC[0], rC[1], rC[2], rD[0], rD[1], rD[2]));
      calc_grad(rC, rA, Xa, gradf);
      calc_grad(rD, rB, Xb, gradg);
      MD_DEBUG(printf("gradf=(%f,%f,%f) gradg=(%f,%f,%f)\n",
		      gradf[0], gradf[1], gradf[2], gradg[0], gradg[1], gradg[2]));
      nf = calc_norm(gradf);
      ng = calc_norm(gradg);
      
      vecg[6] = sqrt(nf/ng);
      for (k1=0; k1 < 3; k1++)
	{
	  vecg[k1] = rC[k1];
	  vecg[k1+3] = rD[k1];
	  rDC[k1] = rD[k1] - rC[k1];
	}
      vecg[7] = calc_norm(rDC)/nf;  
    }
  else
    {
      for (k1 = 0; k1 < 8; k1++)
	vecg[k1] = vecgsup[k1];
    }
  MD_DEBUG(printf("alpha: %f beta: %f\n", vecg[6], vecg[7]));
  newtDist(vecg, 8, &retcheck, funcs2beZeroedDist, i, j, shift); 
  if (retcheck != 0)
    {
      printf("[calcDist] I couldn't calculate distance between %d and %d\n, exiting....\n", i, j);
      //Oparams.time = t;
      //store_bump(i, j);
      //exit(-1);
    }
  for (k1 = 0; k1 < 8; k1++)
    {
      vecgsup[k1] = vecg[k1]; 
    }  
  for (k1 = 0; k1 < 3; k1++)
    {
      r1[k1] = vecg[k1];
      r2[k1] = vecg[k1+3];
      r12[k1] = r1[k1] - r2[k1];
    }
  *alpha = vecg[6];
#if 0
  segno = -1;
  /* se rC è all'interno dell'ellissoide A allora restituisce una distanza negativa*/
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++) 
      segno += (r2[k1]-rA[k1])*Xa[k1][k2]*(r2[k2]-rA[k2]); 
#endif
#if 1
  return calc_norm(r12);
#if 0
  if (segno > 0)
    return calc_norm(r12);
  else
    return -calc_norm(r12);
#endif
#else
#if MD_DEBUG(x) == x
  for (k1 = 0; k1 < 3; k1++)
    rDC[k1] = r1[k1] - r2[k1];
    {
      FILE *f;
      f =fopen("dist.dat","a");
      fprintf(f,"%.15f %.15f\n", Oparams.time, (segno>0)?calc_norm(rDC):-calc_norm(rDC));  
      fclose(f);
    }
  return calc_norm(rDC);
#endif
#endif
}
void rebuildCalendar(void);
int vc_is_pos(int i, int j, double rCx, double rCy, double rCz,
	       double t)
{ 
  double rAC[3], rBC[3], vCA[3], vCB[3], vc;
  double norm[3], wrx, wry, wrz, OmegaA[3][3], OmegaB[3][3];
  double modn;
  double shift[3];
  int na, a, b;
  MD_DEBUG(printf("[bump] t=%f contact point: %f,%f,%f \n", Oparams.time, rxC, ryC, rzC));
  rAC[0] = rx[i] + vx[i]*(t - atomTime[i]) - rCx;
  rAC[1] = ry[i] + vy[i]*(t - atomTime[i]) - rCy;
  rAC[2] = rz[i] + vz[i]*(t - atomTime[i]) - rCz;
#if 1
  for (a=0; a < 3; a++)
    if (fabs (rAC[a]) > L2)
      rAC[a] -= SignR(L, rAC[a]);
#endif
  rBC[0] = rx[j] + vx[j]*(t-atomTime[j]) - rCx;
  rBC[1] = ry[j] + vy[j]*(t-atomTime[j]) - rCy;
  rBC[2] = rz[j] + vz[j]*(t-atomTime[j]) - rCz;
#if 1
  for (a=0; a < 3; a++)
    {
      MD_DEBUG(printf("P rBC[%d]:%.15f ", a, rBC[a]));
      if (fabs (rBC[a]) > L2)
	rBC[a] -= SignR(L, rBC[a]);
      MD_DEBUG(printf("D rBC[%d]:%.15f ", a, rBC[a]));
    }
  MD_DEBUG(printf("\n"));
#endif 
  /* calcola tensore d'inerzia e le matrici delle due quadriche */
  na = (i < Oparams.parnumA)?0:1;
  
  UpdateOrient(i, t-atomTime[i], RA, OmegaA);
  tRDiagR(i, Xa, invaSq[na], invbSq[na], invcSq[na], RA);

  na = (j < Oparams.parnumA)?0:1;
  UpdateOrient(j, t-atomTime[j], RB, OmegaB);
  tRDiagR(j, Xb, invaSq[na], invbSq[na], invcSq[na], RB);
  for (a=0; a < 3; a++)
    {
      norm[a] = 0;
      for (b = 0; b < 3; b++)
	{
	  norm[a] += -Xa[a][b]*rAC[b];
	}
    }
  modn = 0.0;
  for (a = 0; a < 3; a++)
    modn += Sqr(norm[a]);
  modn = sqrt(modn);
  for (a=0; a < 3; a++)
    norm[a] /= modn;
  /* calcola le velocità nel punto di contatto */
  vectProd(wx[i], wy[i], wz[i], -rAC[0], -rAC[1], -rAC[2], &wrx, &wry, &wrz);
  vCA[0] = vx[i] + wrx;
  vCA[1] = vy[i] + wry;
  vCA[2] = vz[i] + wrz;
  vectProd(wx[j], wy[j], wz[j], -rBC[0], -rBC[1], -rBC[2], &wrx, &wry, &wrz);
  vCB[0] = vx[j] + wrx;
  vCB[1] = vy[j] + wry;
  vCB[2] = vz[j] + wrz;
  vc = 0;
  for (a=0; a < 3; a++)
    vc += (vCA[a]-vCB[a])*norm[a];
  MD_DEBUG(printf("VCPOS vc=%.15f\n", vc));
  return (vc > 0);
}
#define EPS 1e-5
void evolveVec(int i, double ti, double *vecout, double *vecin)
{
  double wSq, w, OmegaSq[3][3], M[3][3], Omega[3][3];
  double sinw, cosw;
  int k1, k2;
  wSq = Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]);
  w = sqrt(wSq);
  if (w != 0.0) 
    {
      sinw = sin(w*ti)/w;
      cosw = (1.0 - cos(w*ti))/wSq;
      Omega[0][0] = 0;
      Omega[0][1] = -wz[i];
      Omega[0][2] = wy[i];
      Omega[1][0] = wz[i];
      Omega[1][1] = 0;
      Omega[1][2] = -wx[i];
      Omega[2][0] = -wy[i];
      Omega[2][1] = wx[i];
      Omega[2][2] = 0;
      OmegaSq[0][0] = -Sqr(wy[i]) - Sqr(wz[i]);
      OmegaSq[0][1] = wx[i]*wy[i];
      OmegaSq[0][2] = wx[i]*wz[i];
      OmegaSq[1][0] = wx[i]*wy[i];
      OmegaSq[1][1] = -Sqr(wx[i]) - Sqr(wz[i]);
      OmegaSq[1][2] = wy[i]*wz[i];
      OmegaSq[2][0] = wx[i]*wz[i];
      OmegaSq[2][1] = wy[i]*wz[i];
      OmegaSq[2][2] = -Sqr(wx[i]) - Sqr(wy[i]);
 
      for (k1 = 0; k1 < 3; k1++)
	{
	  for (k2 = 0; k2 < 3; k2++)
	    {
	      Omega[k1][k2] = -Omega[k1][k2];
	      M[k1][k2] = sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	    }
	}
      
      for (k1 = 0; k1 < 3; k1++)
	{
      	  vecout[k1] = vecin[k1];
	  for (k2 = 0; k2 < 3; k2++)
	    vecout[k1] += M[k1][k2]*vecin[k1];
	}
    }
  else
    {
      for (k1 = 0; k1 < 3; k1++)
	{
      	  vecout[k1] = vecin[k1];
	}
    }

}
double calcvecF(int i, int j, double t1, double t, double *r1, double *r2, double* ddot)
{
  int kk;
  double rcat[3], rdbt[3], wra[3], wrb[3];
  evolveVec(i, t-t1, rcat, r1);
  evolveVec(j, t-t1, rdbt, r2);
  ddot[0] = vx[i] - vx[j];
  ddot[1] = vy[i] - vy[j];
  ddot[2] = vx[i] - vz[j];
  vectProd(wx[i], wy[i], wz[i], rcat[0], rcat[1], rcat[2], &wra[0], &wra[1], &wra[2]);
  vectProd(wx[j], wy[j], wz[j], rdbt[0], rdbt[1], rdbt[2], &wrb[0], &wrb[1], &wrb[2]);
  for (kk=0; kk < 3; kk++)
    ddot[kk] += wra[kk] - wrb[kk];
  return calc_norm(ddot);
}


const COORD_TYPE bc1 = 14.0/45.0, bc2 = 64.0/45.0, bc3 = 24.0/45.0;
/* =========================== >>> BodeTerm <<< ============================*/
double BodeTerm(double dt, double* fi)
{
  return dt * (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	       bc1 * fi[4]);
}
int refine_contact(int i, int j, double t, double vecgd[8], double shift[3],double  vecg[5])
{
  int kk, retcheck;

  for (kk = 0; kk < 3; kk++)
    {
      vecg[kk] = (vecgd[kk]+vecgd[kk+3])*0.5; 
    }
  vecg[3] = vecgd[6];
  vecg[4] = t;
  newt(vecg, 5, &retcheck, funcs2beZeroed, i, j, shift); 
  if (retcheck==2)
    {
      MD_DEBUG(printf("newt did not find any contact point!\n"));
      return 0;
    }
  else if (vecg[4] < Oparams.time ||
	   fabs(vecg[4] - Oparams.time)<1E-12 )
    {
      /* se l'urto è nel passato chiaramente va scartato
       * tuttavia se t è minore di zero per errori di roundoff? */
      /* Notare che i centroidi si possono overlappare e quindi t può
       * essere tranquillamente negativo */
      MD_DEBUG(printf("i=%d j=%d <<< vecg[4]=%.15f time:%.15f\n",
		      i, j, vecg[4], Oparams.time));
      return 0;
    }
  else
    {
#if 0
      if (!vc_is_pos(i, j, vecg[0], vecg[1], vecg[2], vecg[4]))
	{
	  MD_DEBUG(printf("vc is positive!\n"));
	  MD_DEBUG(printf("t=%.15f collision predicted %d-%d\n",
			  Oparams.time, i, j));
	  return 0;
	}
#endif
      return 1; 
    }
}
int locate_contact_trivial(int i, int j, double shift[3], double t1, double t2, double vecg[5])
{
  double h, d1, d2, d1Neg, d1Pos, alpha, vecgd1[8], vecgd2[8], vecgd3[8], t, r1[3], r2[3]; 
  double vd, normddot, ddot[3], maxddot, delt;
  const double epsd = 0.001; 
  int foundrc, retcheck, kk;
  t = t1;
  MD_DEBUG(printf("[locate_contact] t1=%f t2=%f shift=(%f,%f,%f)\n", t1, t2, shift[0], shift[1], shift[2]));
  h = 0.01;//EPS*(t2-t1);
  if (lastbump[i]==j && lastbump[j]==i)
    {
      t += h;
      MD_DEBUG(printf("last collision was between %d-%d\n",i,j));
      MD_DEBUG(printf("atomTime[%d]:%.15G atomTime[%d]:%.15G\n", i, atomTime[i], j, atomTime[j])); 
    }
  d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 1);
  MD_DEBUG(printf("distances d1=%.12G", d1));
 #if 1
  d1Pos = calcDist(t, i, j, shift, r1, r2, &alpha, vecgd3, 1);
  MD_DEBUG(printf("distances d1=%.12G d1Neg: %.12G t1=%.15G curtime=%.15G\n", d1Pos, d1,
		  t1, Oparams.time));
  
  if (fabs(d1-d1Pos)>1E-3)
    {
      printf("distances differ d1=%f d1Neg: %f\n", d1Pos, d1);
      //exit(-1);
    }
#endif
#if 1
  if ((lastbump[0]==i || lastbump[0]==j) &&
      (lastbump[1]==i || lastbump[1]==j))
    {
      MD_DEBUG(printf("last collision was between (%d-%d)\n", i, j));
      while (d1 < 0)
	{
	  t += h;
	  if (t > t2)
	    return 0;
	  d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 0);
	}
    }
  else if (d1<0&&fabs(d1)>1E-2)
    {
      printf("d1 < 0 i=%d j=%d\n",i, j);
      return 0;
    }
#endif
  MD_DEBUG(printf(">>>>d1:%f\n", d1));
  foundrc = 0;
#if 0
  if (d1 < 0)
    {
      return 0;
    }
#endif
#if 0
  if (d1 < 0)
    {
      if (refine_contact(i, j, t, vecgd, shift, vecg))
	{
	  MD_DEBUG(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
	  return 1;
	}
      else 
	{
	  MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
	  if (d2 < 0)
	    {
	      MD_DEBUG(printf("d2 < 0 and I did not find contact point, boh...\n"));
	      return 0;
	    }
	}
    }
#endif
  while (t < t2)
    {
#if 0
      normddot = calcvecF(i, j, t1, t, r1, r2, ddot);
      if (normddot!=0)
	t += epsd/normddot;
      
#endif
      t += h;
      d2 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd1, 0);
      MD_DEBUG(printf(">>>> t = %f d1:%f d2:%f\n", t, d1, d2));
      if (d1 > 0 && d2 < 0)
	{
#if 0
	  if (d2 <0)
	    {
	      t -= h;
	      d2 = calcDist(t, i, j, shift, r1, r2, &alpha, vecgd, 0);
	    }
#endif
	  if (refine_contact(i, j, t-h, vecgd2, shift, vecg))
	    {
	      MD_DEBUG(printf("[locate_contact] Adding collision between %d-%d\n", i, j));
	      MD_DEBUG(printf("collision will occur at time %.15G\n", vecg[4])); 
	      if (vecg[4]>t2 || vecg[4]< t1 || 
		  (lastbump[i] == j && lastbump[j]==i && fabs(t - lastcol[i])<1E-12))
		return 0;
	      else
		return 1;
	    }
	  else 
	    {
	      MD_DEBUG(printf("[locate_contact] can't find contact point!\n"));
	      if (d2 < 0)
		{
		  MD_DEBUG(printf("d2 < 0 and I did not find contact point, boh...\n"));
		  return 0;
		  //exit(-1);
		}
	      else
		continue;
	      
	    }
	}
      d1 = d2;
      for (kk = 0; kk < 8; kk++)
	vecgd2[kk] = vecgd1[kk];
    }
  MD_DEBUG(  
  if (foundrc==0)
    printf("%d-%d t=%.12G > t2=%.12G I did not find any contact point!\n", i, j, t, t2);
  )
  return foundrc;
}
int locate_contact(int i, int j, double shift[3], double t1, double t2, double vecg[8])
{
  double h, d1, d2, alpha, vecgd[8], t, r1[3], r2[3]; 
  double normddot, ddot[3], maxddot, delt;
  const double epsd = 0.0001; 
  int foundrc, retcheck, kk;
  t = t1;
  MD_DEBUG(printf("[locate_contact] t1=%f t2=%f shift=(%f,%f,%f)\n", t1, t2, shift[0], shift[1], shift[2]));
  d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd, 1);
  maxddot = sqrt(Sqr(vx[i]-vx[j])+Sqr(vy[i]-vy[j])+Sqr(vz[i]-vz[j])) +
    sqrt(Sqr(wx[i])+Sqr(wy[i])+Sqr(wz[i]))*maxax[i<Oparams.parnumA?0:1]
	 + sqrt(Sqr(wx[j])+Sqr(wy[j])+Sqr(wz[j]))*maxax[j<Oparams.parnumA?0:1];
  if (maxddot==0)
    return 0;
  MD_DEBUG(printf(">>>>d1:%f\n", d1));
  h = EPS*(t2-t1);
  foundrc = 0;
  while (d1 < 0)
    {
      t+=h;
      if (t > t2)
	return 0;
      d1 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd, 0);
    }

  
  while (t < t2)
    {
      delt = d1 / maxddot;
      if ((delt < h)||(d1<10*epsd))
	{
	  while (d1 < 10*epsd)
	    {
	      normddot = calcvecF(i, j, t1, t, r1, r2, ddot);
	      if (normddot!=0)
		t += epsd/normddot;
	      else
		t += h;
	      d2 = calcDistNeg(t, i, j, shift, r1, r2, &alpha, vecgd, 0);
	      MD_DEBUG(printf(">>>> t = %f d1:%f d2:%f\n", t, d1, d2));
	      if (d1 > 0 && d2 < 0)
		{
		  for (kk = 0; kk < 3; kk++)
		    {
		      vecg[kk] = (vecgd[kk]+vecgd[kk+3])*0.5; 
		    }
		  vecg[3] = vecgd[6];
		  vecg[4] = t;
		  newt(vecg, 5, &retcheck, funcs2beZeroed, i, j, shift); 
		  if (retcheck==2 || vecg[4] < Oparams.time ||
		      fabs(vecg[4] - Oparams.time)<1E-12 )
		    {
		      /* se l'urto è nel passato chiaramente va scartato
		       * tuttavia se t è minore di zero per errori di roundoff? */
		      /* Notare che i centroidi si possono overlappare e quindi t può
		       * essere tranquillamente negativo */
		      MD_DEBUG(printf("<<< vecg[4]=%.15f time:%.15f\n",
				      vecg[4], Oparams.time));
		      continue;
		    }
		  else
		    return 1; 
		}
	      d1 = d2;
	    }
	}
      else
	{
	  t += delt;
	  d1 = calcDist(t, i, j, shift, r1, r2, &alpha, vecgd, 1);
	}
    }
  return foundrc;
}
void PredictEvent (int na, int nb) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq, dr[NDIM], dv[NDIM], shift[NDIM], tm[NDIM], vecgd[8],
	 b, d, t, tInt, vv, distSq, t1, t2;
  int et, kk, ii, overlap, mm;
  double ncong, cong[3], pos[3], vecg[5], pos2[3], r1[3], r2[3];
  /*N.B. questo deve diventare un paramtetro in OprogStatus da settare nel file .par!*/
  /*double cells[NDIM];*/
#ifdef MD_GRAVITY
  double Lzx, h1, h2, sig, hh1;
#endif
#if defined(MD_SQWELL) || defined(MD_INFBARRIER) 
  int collCode;
  double sigDeltaSq, intdistSq, distSq, s;
  const double EPSILON = 1E-10;
  double mredl;
#endif
  int cellRangeT[2 * NDIM], signDir[NDIM], evCode,
  iX, iY, iZ, jX, jY, jZ, k, n;

  MD_DEBUG(printf("PredictEvent: %d,%d\n", na, nb));
  MD_DEBUG(calc_energy("PredEv"));
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
#ifdef MD_GRAVITY
  Lzx = OprogStatus.extraLz + Lz;
  if (na < parnumA)
    sig = Oparams.sigma[0][0];
  else
    sig = Oparams.sigma[1][1];
  /* NOTA: Il muro inferiore è posto a Lz / 2 cioè sul fondo del box,
   * in questo modo la base di ogni sfera poggia esattamente 
   * sul fondo della scatola */
  if (inCell[2][na] == 0)
    {
      hh1 =  vz[na] * vz[na] + 2.0 * Oparams.ggrav *
	(rz[na] + Lz2);
      h1 = hh1 -  Oparams.ggrav * sig;
    }
  else
    h1 = hh1 = vz[na] * vz[na] + 2.0 * Oparams.ggrav *
      (rz[na] + Lz2 - inCell[2][na] * (Lzx) / cellsz);
#if 0
  if (evIdB < ATOM_LIMIT)
    {
  if (rx[na] < inCell[0][na]*L/cellsx - L2 || rx[na] > (inCell[0][na]+1)*L/cellsx-L2)
    {
      printf("ERRORE nelle CELLE LUNGO X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("na=%d (%f,%f,%f) cells(%d,%d,%d)\n",
	     na, rx[na], ry[na], rz[na], inCell[0][na], inCell[1][na], inCell[2][na]);
      exit(-1);
    }
   else if (ry[na] < inCell[1][na]*L/cellsy -L2 || ry[na] > (inCell[1][na]+1)*L/cellsy-L2)
     {
       printf("ERRORE nelle CELLE LUNGO Y!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
       printf("na=%d (%f,%f,%f) cells(%d,%d,%d)\n",
 	      na, rx[na], ry[na], rz[na], inCell[0][na], inCell[1][na], inCell[2][na]);
       exit(-1);
     }
   else if (rz[na] < inCell[2][na]*(Lz+OprogStatus.extraLz)/cellsz - Lz2 || rz[na] > (inCell[2][na]+1)*(Lz+OprogStatus.extraLz)/cellsz -Lz2) 
     {
       printf("ERRORE nelle CELLE LUNGO Z!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
       printf("na=%d (%f,%f,%f) cells(%d,%d,%d)\n",
	      na, rx[na], ry[na], rz[na], inCell[0][na], inCell[1][na], inCell[2][na]);
       printf("cellz0: %f cellz1: %f\n", 
	      inCell[2][na]*(Lz+OprogStatus.extraLz)/cellsz -Lz2,
	     (inCell[2][na]+1)*(Lz+OprogStatus.extraLz)/cellsz -Lz2 );
       printf("evIdA: %d evIdB: %d\n", evIdA, evIdB);
       printf("vz[na]:%f\n", vz[na]);
       exit(-1);
     } } 
#endif
  if (vz[na] > 0.0) 
    {
      /* h1 è il Discriminante dell'intersezione con la faccia 
       * inferiore della cella lungo z, h2 di quella superiore,
       * per cui se vz > 0 e h2 > 0 allora si la particella attraverserà
       * la faccia superiore poiché h2 < h1 => t_h2 < t_h1
       * (si noti che la soluzione con tempo positivo se vz > 0 è quella
       * con il + nella formula per la risoluz. di un'eq. di secondo grado */
      h2 = hh1 - 2.0 * Oparams.ggrav * Lzx  / cellsz;
      if (h2 > 0.0) 
	{
	  tm[2] =  (vz[na] - sqrt (h2)) / Oparams.ggrav;
	  signDir[2] = 0;/* signDir = 0 vuol dire che la direzione è
	  		    positiva (faccia superiore in questo caso) */
	} 
      else 
	{
	  tm[2] =  (vz[na] + sqrt (h1)) / Oparams.ggrav;
	  signDir[2] = 1;/* direzione negativa (faccia inferiore in questo caso) */
	}
    } 
  else 
    {
      tm[2] =  (vz[na] + sqrt (h1)) / Oparams.ggrav;
      signDir[2] = 1;
    }
#else
   if (vz[na] != 0.0) 
    {
      if (vz[na] > 0.0) 
	signDir[2] = 0;/* direzione positiva */
      else 
	signDir[2] = 1;/* direzione negativa */
      tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L /
	       cellsz - rz[na] - L2) / vz[na];
    } 
  else 
    tm[2] = timbig;
#endif
  /* end forcefield[k] != 0*/
  
  if (vx[na] != 0.0) 
    {
      if (vx[na] > 0.0) 
	signDir[0] = 0;/* direzione positiva */
      else 
	signDir[0] = 1;/* direzione negativa */
      tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L /
	       cellsx - rx[na] - L2) / vx[na];
    } 
  else 
    tm[0] = timbig;
  
  if (vy[na] != 0.) 
    {
      if (vy[na] > 0.) 
	signDir[1] = 0;
      else 
	signDir[1] = 1;
      tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L /
	       cellsy - ry[na] - L2) / vy[na];
    } 
  else 
    tm[1] = timbig;
  /* ====== */
  /* Find minimum time */
  k = -1; /* giusto per dare un valore ed evitare una warning */
  if (tm[1] <= tm[2]) {
    if (tm[0] <= tm[1]) k = 0;
    else k = 1;
  } else {
    if (tm[0] <= tm[2]) k = 0;
    else k = 2;
  }
  /* Se un errore numerico fa si che tm[k] < 0 allora lo poniamo uguale a 0
   * (ved. articolo Lubachevsky) */
#if 1
  if (tm[k]<0)
    {
      tm[k] = 0.0;
#if 1
      printf("tm[%d]<0 step %lld na=%d\n", k, (long long int)Oparams.curStep, na);
#ifdef MD_GRAVITY
      printf("rz:%f diff:%f\n", rz[na], rz[na]+Lz2);
      printf("h1:%f hh1:%f vz:%f cellz:%d\n", h1, hh1, vz[na], inCell[2][na]);
#endif
      printf("Cells(%d,%d,%d)\n", inCell[0][na], inCell[1][na], inCell[2][na]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
      /*exit(-1);*/
      /*tm[k] = 0.0;*/
#endif
    }
#endif
  /* 100+0 = attraversamento cella lungo x
   * 100+1 =       "           "     "   y
   * 100+2 =       "           "     "   z */
  evCode = 100 + k;
  /* urto con le pareti, il che vuol dire:
   * se lungo z e rz = -L/2 => urto con parete */ 
#ifdef MD_GRAVITY
  if (k == 2 && inCell[2][na] == 0 && signDir[2] == 1) 
    {
      evCode = 4;/* sarebbe 2*k con k=2 (z) e per me vuol dire urto con parete in basso
  		    che è anche l'unica nel presente caso */
      MD_DEBUG2(printf("wall!!! (evIdA: %d)\n", na));
    }
#endif
  MD_DEBUG(printf("schedule event [WallCrossing](%d,%d)\n", na, ATOM_LIMIT+evCode));
  ScheduleEvent (na, ATOM_LIMIT + evCode, Oparams.time + tm[k]);
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_GRAVITY
  /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (inCell[2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
  if (inCell[2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
#endif

  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
#ifndef MD_GRAVITY
      /* apply periodico boundary condition along z if gravitational
       * fiels is not present */
      if (jZ == -1) 
	{
	  jZ = cellsz - 1;    
	  shift[2] = - L;
	} 
      else if (jZ == cellsz) 
	{
	  jZ = 0;    
	  shift[2] = L;
	}
#endif
      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
	{
	  jY = inCell[1][na] + iY;    
	  shift[1] = 0.0;
	  if (jY == -1) 
	    {
	      jY = cellsy - 1;    
	      shift[1] = -L;
	    } 
	  else if (jY == cellsy) 
	    {
	      jY = 0;    
	      shift[1] = L;
	    }
	  for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
	    {
	      jX = inCell[0][na] + iX;    
	      shift[0] = 0.0;
	      if (jX == -1) 
		{
		  jX = cellsx - 1;    
		  shift[0] = - L;
		} 
	      else if (jX == cellsx) 
		{
		  jX = 0;   
		  shift[0] = L;
		}
	      n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	      for (n = cellList[n]; n > -1; n = cellList[n]) 
		{
		  if (n != na && n != nb && (nb >= -1 || n < na)) 
		    {
#if defined(MD_SQWELL) || defined(MD_INFBARRIER)
		      if (na < parnumA && n < parnumA)
			{
			  sigSq = Sqr(Oarams.sigma[0][0]);
			  sigDeltaSq = Sqr(Oparams.sigma[0][0]+Oparams.delta[0][0]);
			  mredl = Mred[0][0];
#if 0
			  inthreshold =  Sqr(Oparams.sigma[0][0]-Oparams.delta[0][0]/2.0);
		          outthreshold = Sqr(Oparams.sigma[0][0]+Oparams.delta[0][0]/2.0);	
#endif
			}
		      else if (na >= parnumA && n >= parnumA)
			{
			  sigSq = Sqr(Oparams.sigma[1][1]);
			  sigDeltaSq = Sqr(Oparams.sigma[1][1]+Oparams.delta[1][1]);
			  mredl = Mred[1][1]; 
#if 0
			  inthreshold =  Sqr(Oparams.sigma[1][1]-Oparams.delta[1][1]/2.0);
		          outthreshold = Sqr(Oparams.sigma[1][1]+Oparams.delta[1][1]/2.0);
#endif
			}
		      else
			{
			  sigSq = Sqr(Oparams.sigma[0][1]);
			  sigDeltaSq = Sqr(Oparams.sigma[0][1]+Oparams.delta[0][1]);
			  mredl = Mred[0][1]; 
#if 0
			  inthreshold =  Sqr(Oparams.sigma[0][1]-Oparams.delta[0][1]/2.0);
		          outthreshold = Sqr(Oparams.sigma[0][1]+Oparams.delta[0][1]/2.0);
#endif
			}
#else
		      /* maxax[...] è il diametro dei centroidi dei due tipi
		       * di ellissoidi */
		      if (na < parnumA && n < parnumA)
			sigSq = Sqr(maxax[0]);
		      else if (na >= parnumA && n >= parnumA)
			sigSq = Sqr(maxax[1]);
		      else
			sigSq = Sqr((maxax[0]+maxax[1])*0.5);
		      MD_DEBUG2(printf("sigSq: %f\n", sigSq));
#endif
		      tInt = Oparams.time - atomTime[n];
		      dr[0] = rx[na] - (rx[n] + vx[n] * tInt) - shift[0];	  
		      dv[0] = vx[na] - vx[n];
		      dr[1] = ry[na] - (ry[n] + vy[n] * tInt) - shift[1];
		      dv[1] = vy[na] - vy[n];
#ifdef MD_GRAVITY
		      dr[2] = rz[na] - 
			(rz[n] + (vz[n] - 0.5 * Oparams.ggrav * tInt) * tInt) - shift[2];
		      dv[2] = vz[na] - (vz[n] - Oparams.ggrav * tInt);
#else
		      dr[2] = rz[na] - (rz[n] + vz[n] * tInt) - shift[2];
		      dv[2] = vz[na] - vz[n];

#endif
     		      b = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
#if defined(MD_SQWELL)|| defined(MD_INFBARRIER)
		      distSq = Sqr(dr[0]) + Sqr(dr[1]) + Sqr(dr[2]);
		      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
		      collCode = MD_EVENT_NONE;
		      if (!bound(n, na))
			{
			  if ( b < 0.0 ) 
			    {
			      /* la piccola correzione serve poichè a causa
			       * di errori numerici dopo l'evento la particella
			       * potrebbe essere ancora fuori dalla buca (distSq > sigDeltaSq)*/
			      d = Sqr (b) - vv * (distSq - sigDeltaSq);
			      if (d > 0.0)
				{
				  t = (-sqrt (d) - b) / vv;
				  if (t > 0 || (t < 0 && distSq < sigDeltaSq))
				    collCode = MD_OUTIN_BARRIER;
				}
			    }
			}
		      else
			{
#ifdef MD_INFBARRIER
			  if (sigDeltaSq == 0)
			    goto no_core_bump;
			  
#endif
		  	    
			  if (b < 0.0)
			    { 
		      	      d= Sqr (b) - vv * (distSq - sigSq);
	      		      if (d > 0.0)
      				{
				  t = (-sqrt (d) - b) / vv;
				  if (t > 0 || (t < 0 && distSq < sigSq))
				    collCode = MD_CORE_BARRIER;
				}
			    }
#ifdef MD_INFBARRIER
no_core_bump:
#endif
			  if (collCode == MD_EVENT_NONE)
			    {
			      d = Sqr (b) - vv * (distSq - sigDeltaSq);
			      if (d > 0.0)
				{
				  t = ( sqrt (d) - b) / vv;
				  if (t > 0 || (t < 0 && distSq > sigDeltaSq))
				    collCode = MD_INOUT_BARRIER;
				}
			    }
			}
		      if (t < 0 && collCode!= MD_EVENT_NONE)
			{
#if 1
			  printf("time:%.15f tInt:%.15f t:%.20f\n", Oparams.time,
				 tInt, t);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
	     					      Sqr(dr[2])));
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
			  printf("collCode: %d\n", collCode);
			  //exit(-1);
#endif
			  t = 0;
			}

		      if (collCode != MD_EVENT_NONE)
			{
			  ScheduleEventBarr (na, n, collCode, Oparams.time + t);
			  MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
		      	} 

#else
		      distSq = Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]);
		      vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
	    	      d = Sqr (b) - vv * (distSq - sigSq);
		
		      if (d < 0 || (b > 0.0 && distSq > sigSq)) 
			{
			  /* i centroidi non collidono per cui non ci può essere
			   * nessun urto sotto tali condizioni */
			  continue;
			}
#if 1
		      MD_DEBUG(printf("PREDICTING na=%d n=%d\n", na , n));
		      if (vv==0.0)
			{
			  if (distSq >= sigSq)
			    continue;
			    /* la vel relativa è zero e i centroidi non si overlappano quindi
			     * non si possono urtare! */
			  t1 = t = 0;
			  t2 = 10.0;/* anche se sono fermi l'uno rispetto all'altro possono 
				       urtare ruotando */
			}
		      else if (distSq >= sigSq)
			{
			  t = t1 = - (sqrt (d) + b) / vv;
			  t2 = (sqrt (d) - b) / vv;
			  overlap = 0;
			}
		      else 
			{
			  MD_DEBUG(printf("Centroids overlap!\n"));
			  t2 = t = (sqrt (d) - b) / vv;
			  t1 = 0.0; 
			  overlap = 1;
			  MD_DEBUG(printf("altro d=%f t=%.15f\n", d, (-sqrt (d) - b) / vv));
			  MD_DEBUG(printf("vv=%f dv[0]:%f\n", vv, dv[0]));
			}
		      MD_DEBUG(printf("t=%f curtime: %f b=%f d=%f\n", t, Oparams.time, b ,d));
		      MD_DEBUG(printf("dr=(%f,%f,%f) sigSq: %f", dr[0], dr[1], dr[2], sigSq));
		      //t += Oparams.time; 
		      t2 += Oparams.time;
		      t1 += Oparams.time;
#if 0
		      /* t è il guess per il newton-raphson */
		      /* come guess per x possiamo usare il punto di contatto 
		       * fra i centroidi */
		      /* vecg è un guess per il vettore a 5 dimensioni (x, alpha ,t) */
		      /* NOTA: qui va calcolato il vettore guess vecg 
		       * che è il punto di contatto dei due centroidi */
		      cong[0] = rx[na] + vx[na] * (t-atomTime[na]) 
			- (rx[n] + vx[n] * (t-atomTime[n])) - shift[0];
		      cong[1] = ry[na] + vy[na] * (t-atomTime[na]) 
			- (ry[n] + vy[n] * (t-atomTime[n])) - shift[1];
		      cong[2] = rz[na] + vz[na] * (t-atomTime[na]) 
			- (rz[n] + vz[n] * (t-atomTime[n])) - shift[2];
		      ncong=0.0;
		      for (kk=0; kk < 3; kk++)
			ncong +=  Sqr(cong[kk]);
		      ncong = sqrt(ncong);
		      for (kk=0; kk < 3; kk++)
			{
			  cong[kk] /= ncong;
			}
		      pos[0] =  rx[na] + vx[na] * (t-atomTime[na]);
		      pos[1] =  ry[na] + vy[na] * (t-atomTime[na]);
		      pos[2] =  rz[na] + vz[na] * (t-atomTime[na]);
		      for (kk=0; kk < 3; kk++)
			{
			  vecg[kk] = pos[kk] - 0.5*cong[kk]*maxax[na<Oparams.parnumA?0:1];
			  pos[kk] = vecg[kk];
			}
		      MD_DEBUG(printf("shift (%f, %f, %f) vecg (%f, %f, %f)\n", shift[0], shift[1], shift[2], vecg[0], vecg[1], vecg[2]));
		      MD_DEBUG2(printf("r[%d](%f,%f,%f)-r[%d](%f,%f,%f)\n",
				       na, rx[na], ry[na], rz[na], n, rx[n], ry[n], rz[n]));
		      
#if 0
		      t = Oparams.time;
		      vecg[0] = (rx[na] + rx[n] + shift[0])*0.5;
		      vecg[1] = (ry[na] + ry[n] + shift[1])*0.5;
		      vecg[2] = (rz[na] + rz[n] + shift[2])*0.5;
#endif
		      vecg[3] = 1.0; /* questa stima di alpha andrebbe fatta meglio!*/
		      vecg[4] = t;
		      
		      MD_DEBUG(printf("time=%.15f vecguess: %f,%f,%f alpha=%f t=%f\n",Oparams.time, vecg[0], vecg[1], vecg[2], vecg[3],vecg[4]));
#endif
#endif		      
		      //calcDist(Oparams.time, na, n, shift, r1, r2);
		      //continue;
		      //exit(-1);
		      if (!locate_contact(na, n, shift, t1, t2, vecg))
			continue;
		     			  
#if 0
	    	      gd = (d2 - d1)/h; 
    		      if (fabs(gd) > 1E-12)
			delt = (epsd*(maxax[0]+maxax[1])*0.5)/gd;
		      d1o = d1;
		      d2o = d2;
#endif
#if 0
		      for (kk=0; kk < 5; kk++)
		      	vecgold[kk] = vecg[kk];
	     	      for(mm=0; mm < 10000; mm++)
			{
			  MD_DEBUG(printf(">>> mm=%d\n", mm));
			  newt(vecg, 5, &retcheck, funcs2beZeroed, na, n, shift); 
			  if (retcheck==2 || vecg[4] < Oparams.time ||
			      fabs(vecg[4] - Oparams.time)<1E-12 )
			    {
			      /* se l'urto è nel passato chiaramente va scartato
			       * tuttavia se t è minore di zero per errori di roundoff? */
			      /* Notare che i centroidi si possono overlappare e quindi t può
			       * essere tranquillamente negativo */
			      MD_DEBUG(printf("<<< vecg[4]=%.15f time:%.15f\n",
					      vecg[4], Oparams.time));
			      break;
			    }
	  		  if (vc_is_pos(na, n, vecg[0], vecg[1], vecg[2], vecg[4]))
			    {
			      MD_DEBUG(printf("vc is positive!\n"));
			      MD_DEBUG(printf("t=%.15f collision predicted %d-%d\n",
					      Oparams.time, na, n));
			      break;
			    }
			  else if (mm == 0)
			    {
			      MD_DEBUG(printf("vc is negative!\n"));
			      for (kk=0; kk < 5; kk++)
				vecgf[kk] = vecg[kk];
			    }
			  for (kk=0; kk < 5; kk++)
			    vecg[kk] = vecgold[kk] + (vecgf[kk] - vecgold[kk])/(10002-mm);
			
			}
		      if (mm==10000)
			{
			  printf("MAX VCPOS!!!!\n");
			  printf("%d-%d t=%.15G\n", na, n, Oparams.time);
			  continue;
			}
		      else
			{
			  printf("OK maybe got contact point mm=%d\n", mm);
			}
#endif
#if 0
		      if (retcheck==2 || vecg[4] < Oparams.time ||
	    		  fabs(vecg[4] - Oparams.time)<1E-12 )
    			{
			  /* se l'urto è nel passato chiaramente va scartato
			   * tuttavia se t è minore di zero per errori di roundoff? */
			  /* Notare che i centroidi si possono overlappare e quindi t può
			   * essere tranquillamente negativo */
			  MD_DEBUG(printf("<<< vecg[4]=%.15f time:%.15f\n",
					  vecg[4], Oparams.time));
			  continue;
			}
#endif		      
		      rxC = vecg[0];
		      ryC = vecg[1];
		      rzC = vecg[2];
		      MD_DEBUG(printf("A x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)-B x(%.15f,%.15f,%.15f) v(%.15f,%.15f,%.15f)",
				      rx[na], ry[na], rz[na], vx[na], vy[na], vz[na],
				      rx[n], ry[n], rz[n], vx[n], vy[n], vz[n]));
		      t = vecg[4];
#if 1
		      if (t < 0)
			{
#if 1
			  printf("time:%.15f tInt:%.15f\n", Oparams.time,
				 tInt);
			  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
	     					      Sqr(dr[2]))-1.0 );
			  printf("STEP: %lld\n", (long long int)Oparams.curStep);
			  printf("atomTime: %.10f \n", atomTime[n]);
			  printf("n:%d na:%d\n", n, na);
			  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
#endif
			  t = 0;
			}
#endif
		      /* il tempo restituito da newt() è già un tempo assoluto */
		      MD_DEBUG(printf("time: %f Adding collision %d-%d\n", Oparams.time, na, n));
		      ScheduleEvent (na, n, t);
		      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
#endif
		    }
		} 
	    }
	}
    }
}

#ifdef MD_GRAVITY
void ProcessCollWall(void)
{
  /* Dissipation */
  MD_DEBUG2(printf("Collision with wall evIdA: %d vz: %.15f\n", evIdA, vz[evIdA]));
#if 0
  printf("timeNow-lastcol[%d]:%.15f\n", evIdA, Oparams.time - lastcol[evIdA]);
  printf("walldiss:%.15f vz:%.10f\n", Oparams.wallDiss, vz[evIdA]);
  printf("timeNow:%.10f lastcol:%.10f\n", Oparams.time, lastcol[evIdA]);
#endif
  if (Oparams.time - lastcol[evIdA] < OprogStatus.tc)
    {
      vz[evIdA] = -vz[evIdA];
    }
  else
    {
      vz[evIdA] = -Oparams.wallDiss*vz[evIdA];
    }	 
  lastcol[evIdA] = Oparams.time;
}
#endif
void calc_energynew(char *msg)
{
  int i, k1, k2;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#else
  double Ia, Ib;
#endif
  K = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (i<Oparams.parnumA)
	{
	  /* calcola tensore d'inerzia e le matrici delle due quadriche */
	  K += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    {
#ifdef MD_ASYM_ITENS
      	      K += wt[k1]*ItensD[0][k1]*wt[k1];
#else
	      K += Oparams.I[0]*Sqr(wt[k1]);
#endif
	    }
	}
      else
	{
	  K += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
	  for (k1=0; k1 < 3; k1++)
	    {
#ifdef MD_ASYM_ITENS
      	      K += wt[k1]*ItensD[1][k1]*wt[k1];
#else
      	      K += Oparams.I[1]*Sqr(wt[k1]);
#endif
	    }
	}
    }
  K *= 0.5;
#ifdef MD_ASYM_ITENS
  free_matrix(Ia,3);
  free_matrix(Ib,3);
#endif
  printf("[%s] Kinetic Energy: %f\n", msg, K);
}
void calc_energy_i(char *msg, int i)
{
  int k1, k2;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#else
  double Ia, Ib;
#endif
  K = 0;
  if (i<Oparams.parnumA)
    {
      /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
      RDiagtR(i, Ia, ItensD[0][0], ItensD[0][1], ItensD[0][2], R[i]);
#endif
      K += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
      wt[0] = wx[i];
      wt[1] = wy[i];
      wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    K += wt[k1]*Ia[k1][k2]*wt[k2];
	  }
#else
      for (k1=0; k1 < 3; k1++)
	K += Sqr(wt[k1])*Oparams.I[0];
#endif
    }
  else
    {
#ifdef MD_ASYM_ITENS
      RDiagtR(i, Ib, ItensD[1][0], ItensD[1][1], ItensD[1][2], R[i]);
#endif
      K += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
      wt[0] = wx[i];
      wt[1] = wy[i];
      wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
      for (k1=0; k1 < 3; k1++)
	for (k2=0; k2 < 3; k2++)
	  {
	    K += wt[k1]*Ib[k1][k2]*wt[k2];
	  }
#else
      for (k1=0; k1 < 3; k1++)
	K += Sqr(wt[k1])*Oparams.I[1];
#endif
    }
  K *= 0.5;
#ifdef MD_ASYM_ITENS
  free_matrix(Ia,3);
  free_matrix(Ib,3);
#endif
  printf("[%s] Kinetic Energy of %d: %f\n", msg, i, K);
  
}
void calc_energy(char *msg)
{
  int i, k1, k2;
  double wt[3];
#ifdef MD_ASYM_ITENS
  double **Ia, **Ib;
  Ia = matrix(3,3); 
  Ib = matrix(3,3);
#else
  double Ia, Ib;
#endif

  K = 0;
  for (i=0; i < Oparams.parnum; i++)
    {
      if (i<Oparams.parnumA)
	{
	  /* calcola tensore d'inerzia e le matrici delle due quadriche */
#ifdef MD_ASYM_ITENS
	  RDiagtR(i, Ia, ItensD[0][0], ItensD[0][1], ItensD[0][2], R[i]);
#endif
	  K += Oparams.m[0]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		K += wt[k1]*Ia[k1][k2]*wt[k2];
	      }
#else
	  for (k1=0; k1 < 3; k1++)
	    K += Sqr(wt[k1])*Oparams.I[0];
#endif
	}
      else
	{
#ifdef MD_ASYM_ITENS
	  RDiagtR(i, Ib, ItensD[1][0], ItensD[1][1], ItensD[1][2], R[i]);
#endif
	  K += Oparams.m[1]*(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i]));  
	  wt[0] = wx[i];
	  wt[1] = wy[i];
	  wt[2] = wz[i];
#ifdef MD_ASYM_ITENS
	  for (k1=0; k1 < 3; k1++)
	    for (k2=0; k2 < 3; k2++)
	      {
		K += wt[k1]*Ib[k1][k2]*wt[k2];
	      }
#else
	  for (k1=0; k1 < 3; k1++)
	    K += Sqr(wt[k1])*Oparams.I[1];
#endif
	}
    }
  K *= 0.5;
#ifdef MD_ASYM_ITENS
  free_matrix(Ia,3);
  free_matrix(Ib,3);
#endif
  printf("[%s] Kinetic Energy: %f\n", msg, K);
}
void store_bump(int i, int j)
{
  char fileop2[512], fileop[512];
  int ii;
  FILE *bf;
  const char tipodat2[]= "%.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G %.15G\n";
  sprintf(fileop2 ,"StoreBump-%d-%d-t%.8f", i, j, Oparams.time);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (bf = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
      exit(-1);
    }
  UpdateSystem();
  R2u();
  fprintf(bf, ".Vol: %f\n", L*L*L);
  MD_DEBUG(printf("[Store bump]: %.15G\n", Oparams.time));
  for (ii = 0; ii < Oparams.parnum; ii++)
    {
      if (ii==i || ii==j)
	fprintf(bf, tipodat2,rx[ii], ry[ii], rz[ii], uxx[ii], uxy[ii], uxz[ii], uyx[ii], uyy[ii], 
		uyz[ii], uzx[ii], uzy[ii], uzz[ii]);
    }
  //writeAllCor(bf);
  fprintf(bf,"%.15f %.15f %.15f @ 0.1 C[green]\n", rxC, ryC, rzC);
  fclose(bf);

}
void ProcessCollision(void)
{
  int k;
  UpdateAtom(evIdA);
  UpdateAtom(evIdB);
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  MD_DEBUG(calc_energy("prima"));
#if defined(MD_SQWELL)||defined(MD_INFBARRIER)
  /* i primi due bit sono il tipo di event (uscit buca, entrata buca, collisione con core 
   * mentre nei bit restanti c'e' la particella con cui tale evento e' avvenuto */
  bump(evIdA, evIdB, &W, evIdC);
#else
  bump(evIdA, evIdB, rxC, ryC, rzC, &W);
#endif
  MD_DEBUG(calc_energy("dopo"));
  MD_DEBUG(store_bump(evIdA, evIdB));
    //ENDSIM=1;
  /*printf("qui time: %.15f\n", Oparams.time);*/
#ifdef MD_GRAVITY
  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
#else
  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
  lastbump[evIdA]=evIdB;
  lastbump[evIdB]=evIdA;
#endif
  PredictEvent(evIdA, -1);
  PredictEvent(evIdB, evIdA);
}
void docellcross(int k, double velk, double *rkptr, int cellsk)
{
#if 0
  if (inCell[0][evIdA]+1> cellsx ||inCell[1][evIdA]+1> cellsy||inCell[2][evIdA]+1> cellsz) 
    {printf("PRIMAin cell cross ?!?\n");
    printf("velk: %f (%d,%d,%d) (%d,%d,%d) k=%d cellsk:%d\n",velk,  cellsx , cellsy,cellsz,
    inCell[0][evIdA],inCell[1][evIdA], inCell[2][evIdA], k, cellsk );}
#endif
    if (velk > 0.0)
    {
      inCell[k][evIdA] = inCell[k][evIdA] + 1;
      cellRange[2 * k] = 1;
      if (inCell[k][evIdA] == cellsk) 
	{
	  inCell[k][evIdA] = 0;
#ifdef MD_GRAVITY
	  if (k==2)
	    {
	      printf("Un particella ha superato la massima altezza consentita (%.15f)\n",
		     Lz2 + OprogStatus.extraLz);
	      printf("Aumentare il parametro extraLz e rilanciare la simulazione\n");
	      exit(-1);
	    }
#endif
	  *rkptr = -L2;
	  OprogStatus.DR[evIdA][k]++;
	}

   }
  else
    { 
      cellRange[2 * k + 1] = -1;
      inCell[k][evIdA] = inCell[k][evIdA] - 1;
      if (inCell[k][evIdA] == -1) 
	{
	  inCell[k][evIdA] = cellsk - 1;
	  *rkptr = L2;
	  OprogStatus.DR[evIdA][k]--;
	}
    }
#if 0
  if (inCell[0][evIdA]> cellsx ||inCell[1][evIdA]> cellsy||inCell[2][evIdA]> cellsz) 
    {printf("in cell cross ?!?\n");
    printf("velk: %f(%d,%d,%d) (%d,%d,%d) k=%d cellsk:%d\n",  velk,cellsx , cellsy,cellsz,
    inCell[0][evIdA],inCell[1][evIdA], inCell[2][evIdA], k, cellsk );}
#endif
}
void ProcessCellCrossing(void)
{
#ifdef MD_GRAVITY
  int j; 
#endif
  int k, n;

  UpdateAtom(evIdA);
  /* NOTA: cellList[i] con 0 < i < Oparams.parnum è la cella in cui si trova la particella
   * i-esima mentre cellList[j] con 
   * Oparams.parnum <= j < cellsx*cellsy*cellsz+Oparams.parnum
   * è la prima particella che si trova nella cella j-esima
    */
  n = (inCell[2][evIdA] * cellsy + inCell[1][evIdA] )*cellsx + inCell[0][evIdA]
    + Oparams.parnum;
  
  while (cellList[n] != evIdA) 
    n = cellList[n];
  /* Eliminazione di evIdA dalla lista della cella n-esima */
  cellList[n] = cellList[evIdA];
  for (k = 0; k < NDIM; k++)
    { 
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
#ifdef MD_GRAVITY
  j = evIdB - ATOM_LIMIT;
  if (j >= 100)
    {
      k = j - 100; 
      switch (k)
	{
	case 0: 
	  docellcross(0, vx[evIdA], &(rx[evIdA]), cellsx);
	  break;
	case 1: 
	  docellcross(1, vy[evIdA], &(ry[evIdA]), cellsy);
	  break;
	case 2:
	  docellcross(2, vz[evIdA], &(rz[evIdA]), cellsz);
	  break;
	}
    }
  else
    {
      k = j / 2;
      cellRange[j] = 0;
      ProcessCollWall();
    }
#else
  k = evIdB - 100 - ATOM_LIMIT; 
#if 0
  if (inCell[0][evIdA]> cellsx ||inCell[1][evIdA]> cellsy||inCell[2][evIdA]> cellsz) 
  printf("Cells(%d,%d,%d)\n", inCell[0][evIdA],inCell[1][evIdA],inCell[2][evIdA]);
#endif
  switch (k)
    {
    case 0: 
      docellcross(0, vx[evIdA], &(rx[evIdA]), cellsx);
      break;
    case 1: 
      docellcross(1, vy[evIdA], &(ry[evIdA]), cellsy);
      break;
    case 2:
      docellcross(2, vz[evIdA], &(rz[evIdA]), cellsz);
      break;
    }
#endif
  PredictEvent(evIdA, evIdB);
  n = (inCell[2][evIdA] * cellsy + inCell[1][evIdA])*cellsx + 
    inCell[0][evIdA] + Oparams.parnum;
  /* Inserimento di evIdA nella nuova cella (head) */
  cellList[evIdA] = cellList[n];
  cellList[n] = evIdA;
}
void velsBrown(double T)
{
  comvel_brown(T, Oparams.m); 
}

void rebuildLinkedList(void)
{
  int j, n;
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList[j] = -1;
  /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
  for (n = 0; n < Oparams.parnum; n++)
    {
      atomTime[n] = Oparams.time;
      inCell[0][n] =  (rx[n] + L2) * cellsx / L;
      inCell[1][n] =  (ry[n] + L2) * cellsy / L;
#ifdef MD_GRAVITY
      inCell[2][n] =  (rz[n] + Lz2) * cellsz / (Lz+OprogStatus.extraLz);
#else
      inCell[2][n] =  (rz[n] + L2)  * cellsz / L;
#endif
      j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
	inCell[0][n] + Oparams.parnum;
      cellList[n] = cellList[j];
      cellList[j] = n;
    }
}
void rebuildCalendar(void)
{
  int k, n;
  
  InitEventList();
  for (k = 0;  k < NDIM; k++)
    {
      cellRange[2*k]   = - 1;
      cellRange[2*k+1] =   1;
    }
  for (n = 0; n < Oparams.parnum; n++)
    PredictEvent(n, -2); 
}
void distanza(int ia, int ib)
{
  double dx, dy, dz;
  dx = rx[ia]-rx[ib];
  dy = ry[ia]-ry[ib];
  dz = rz[ia]-rz[ib];
  dx = dx - L*rint(dx/L);
  dy = dx - L*rint(dy/L);
  printf("dist(%d,%d): %f\n", ia, ib, sqrt(Sqr(dx)+Sqr(dy)+Sqr(dz)));
}
void rebuildLinkedList(void);
#if defined(MD_SQWELL) && defined(MD_BONDCORR)
extern int **bonds0, *numbonds0; 
int bondhist[3]={0,0,0};
extern double corrnorm, corrini3, corrini1, corrini2, *lastbreak1, *lastbreak2;
#endif
/* ============================ >>> move<<< =================================*/
void move(void)
{
  char fileop[1024], fileop2[1024], fileop3[1024];
  FILE *bf;
#if defined(MD_SQWELL) && defined(MD_BONDCORR)
  FILE *bof;
  int j, cc;
  double cb, corr1, corr2, corr3;
#endif
  const char sepStr[] = "@@@\n";
  int db, i, innerstep=0;
#ifdef MD_GRAVITY
  int ii;
  double rzmax, zfact;
#endif

  /* Zero all components of pressure tensor */
#if 0
  Wxy = 0.0;
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;
#endif
  /* get next event */
  while (!ENDSIM)
    {
      innerstep++;
      NextEvent();
      /* Descrizione Eventi:
       * 0 <= evIdB < ATOM_LIMIT: 
       *        urto fra evIdA e evIdB 
       *
       * ATOM_LIMIT <= evIdB <= ATOM_LIMIT + 5:
       *        ATOM_LIMIT   -> urto con parete lungo x nella direzione negativa
       *        ATOM_LIMIT+1 -> urto con parete lungo x nella direzione positiva
       *
       * ATOM_LIMIT + 5 < evIdB < ATOM_LIMT + 100:
       *        eventi per usi vari (misure, output e altro)
       *
       * evIdB >= ATOM_LIMIT+100: 
       *        attraversamento della cella (cell-crossing) */        
      /* PROCESS EVENTS */ 
      if (evIdB < ATOM_LIMIT)
	{
	  MD_DEBUG(printf("collision (evIdA: %d evIdB:%d)\n", evIdA, evIdB));
	  ProcessCollision();
	  OprogStatus.collCount++;
	}

#ifdef MD_GRAVITY
      else if (evIdB >= ATOM_LIMIT + 100 || evIdB < ATOM_LIMIT + NDIM * 2)
	{
	  ProcessCellCrossing();
	  OprogStatus.crossCount++;
	}
#else
      else if ( evIdB >= ATOM_LIMIT + 100 )
	{
	  ProcessCellCrossing();
	  OprogStatus.crossCount++;
	}
#endif
      /* ATOM_LIMIT +6 <= evIdB < ATOM_LIMIT+100 eventi che si possono usare 
       * liberamente */
      else if (evIdB == ATOM_LIMIT + 7)
	{
	  UpdateSystem();
	  OprogStatus.nextSumTime += OprogStatus.intervalSum;
	  ScheduleEvent(-1, ATOM_LIMIT + 7, OprogStatus.nextSumTime);
#if defined(MD_SQWELL) && defined(MD_BONDCORR)
	    {
	      FILE *ppf;
	      sprintf(fileop2 ,"population.dat");
	      /* store conf */
	      strcpy(fileop, absTmpAsciiHD(fileop2));
	      if ( (ppf = fopenMPI(fileop, "w")) == NULL)
		{
		  mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
		  exit(-1);
		}
	      fprintf(ppf,"%d %d %d\n", bondhist[0], bondhist[1], bondhist[2]);
	      fclose(ppf); 
	    }
#endif
	  /*calcObserv();*/
	  outputSummary(); 
	}
      else if (evIdB == ATOM_LIMIT + 8)
	{
#if defined(MD_BONDCORR) && defined(MD_SQWELL)  
	  cb = 0;
	  corr1 = corr2 = corr3 = 0;
	  for (i=0; i < Oparams.parnum; i++)
	    {
	      for (j=0; j < numbonds0[i]; j++)
		{
		  if (bound(i,bonds0[i][j]))
		    cb++;
		}
	      cc = 0;
    	      for (j=0; j < numbonds0[i]; j++)
		{
		  if (bound(i,bonds0[i][j]))
		    cc++;
		}
	      if (cc > 0)
		{
		  if (numbonds0[i]==2)
	      	    corr2++;
		  if (numbonds0[i]==3)
		    corr3++;
		}
	      if (numbonds0[i]==1 && numbonds0[bonds0[i][0]]==1)
		{
		  if (numbonds[i]==1)	
		    corr1++;
		}	
	    }
	  sprintf(fileop2 ,"BondCorrFuncB2.dat");
	  /* store conf */
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	  if ( (bof = fopenMPI(fileop, "a")) == NULL)
	    {
	      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	      exit(-1);
	    }
	  fprintf(bof,"%.15f %.15f\n", Oparams.time, corr2/corrini2);
	  fclose(bof);
	  sprintf(fileop2 ,"BondCorrFuncB3.dat");
	  /* store conf */
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	  if ( (bof = fopenMPI(fileop, "a")) == NULL)
	    {
	      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	      exit(-1);
	    }
	  fprintf(bof,"%.15f %.15f\n", Oparams.time, corr3/corrini3);
	  fclose(bof);
  
	  sprintf(fileop2 ,"bondcorr.dat");
	  /* store conf */
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	  if ( (bf = fopenMPI(fileop, "a")) == NULL)
	    {
	      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	      exit(-1);
	    }
	  if (corrnorm)
	    fprintf(bf, "%.15f %.8f\n", Oparams.time, cb/corrnorm);
	  fclose(bf);
#else
	  sprintf(fileop2 ,"Store-%d-%d", 
		  OprogStatus.KK, OprogStatus.JJ);
	  /* store conf */
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	  if ( (bf = fopenMPI(fileop, "w")) == NULL)
	    {
	      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	      exit(-1);
	    }
	  UpdateSystem();
	  R2u();
#ifndef MD_STOREMGL
	  writeAsciiPars(bf, opro_ascii);
	  fprintf(bf, sepStr);
	  writeAsciiPars(bf, opar_ascii);
	  fprintf(bf, sepStr);
	  printf("qui\n");
#endif
	  MD_DEBUG(printf("[Store event]: %.15G JJ=%d KK=%d\n", Oparams.time, OprogStatus.JJ, OprogStatus.KK));
	  fprintf(bf, ".Vol: %f\n", L*L*L);
	  writeAllCor(bf);
	  fclose(bf);
#ifndef MD_STOREMGL
#ifdef MPI
          sprintf(fileop3, "/bin/gzip -f %s_R%d", fileop, my_rank);
#else 
          sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
	  system(fileop3);
#endif
#endif
	  OprogStatus.JJ++;
	  if (OprogStatus.JJ == OprogStatus.NN)
	    {
	      OprogStatus.JJ = 0;
	      OprogStatus.KK++;
	    }
          OprogStatus.nextStoreTime = OprogStatus.storerate *
	    (pow(OprogStatus.base,OprogStatus.NN)*OprogStatus.KK+pow(OprogStatus.base,OprogStatus.JJ));
	  ScheduleEvent(-1, ATOM_LIMIT + 8, OprogStatus.nextStoreTime);
	}
      else if (evIdB == ATOM_LIMIT + 10)
	{
	  UpdateSystem();
	  R2u();
#if 1
	    {
	      static double shift[3] = {0,0,0}, vecg[8], vecgNeg[8];
	      double d,r1[3], r2[3], alpha;
	      FILE* f;
	      static int first = 1;
	      shift[0] = L*rint((rx[0]-rx[1])/L);
	      shift[1] = L*rint((ry[0]-ry[1])/L);
	      shift[2] = L*rint((rz[0]-rz[1])/L);
	      MD_DEBUG(printf("[EVENT10] shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]));
	      d=calcDist(Oparams.time, 0, 1, shift, r1, r2, &alpha, vecg, 1);
	      if (first)
		f = fopen("distPos.dat","w");
	      else
		f = fopen("distPos.dat","a");
	      fprintf(f,"%.15G %.15G\n", Oparams.time, d);
	      fclose(f);
	      d=calcDistNeg(Oparams.time, 0, 1, shift, r1, r2, &alpha, vecgNeg, 1);
	      if (first)
		f = fopen("distNeg.dat","w");
	      else
		f = fopen("distNeg.dat","a");
	      fprintf(f,"%.15G %.15G\n", Oparams.time, d);
	      fclose(f);
	      if (first == 1)
		first = 0;
	    }
#endif
	  if (OprogStatus.brownian)
	    {
	      velsBrown(Oparams.T);
	      rebuildCalendar();
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	    }
	  OprogStatus.nextDt += Oparams.Dt;
	  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	  break;
	}
#ifdef MD_GRAVITY
      else if (evIdB == ATOM_LIMIT + 9)
	{
	  UpdateSystem();
	  if (OprogStatus.taptau > 0.0)
	    {
	      if (OprogStatus.quenchend < 0.0)
		{
#if 0
		  if ((V - Vold)/V < OprogStatus.quenchtol)
		    {
		      printf("QUENCH DONE! %d\n", Oparams.curStep);
		      /* se l'energia potenziale è ormai stabile considera il quench finito */
		      OprogStatus.quenchend = Oparams.curStep;
		    }
#else
		  calcKVz();
		  OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
		  if ( (2.0*K/(3.0*((double)Oparams.parnum)-3.0)) < 
		       OprogStatus.quenchtol)
		    {
#endif
		      printf("QUENCH DONE! %lld\n", (long long int) Oparams.curStep);
		      OprogStatus.numquench++;
		      /* calcola e salva le misure quando finisce il quench */
		      calcRho();
		      save_rho();
		      calccmz();
		      save_rzcm();
		      OprogStatus.quenchend = Oparams.time;
		      comvel(Oparams.parnum, Oparams.T, Oparams.m, 0);
#if 1
		      calcKVz();
		      scalevels(Oparams.T, K, Vz);
#endif
		      MD_DEBUG3(printf("rzmax:%f\n", rzmax));
		      rzmax = -Lz2;
		      for (ii=0; ii < Oparams.parnum; ii++)
			{
			  if (rz[ii] > rzmax)
			    rzmax = rz[ii];
			}
		      if (Lz / (rzmax+Lz2) < OprogStatus.expandFact)
			zfact = Lz/(rzmax+Lz2);
		      else
			zfact = OprogStatus.expandFact;
		      for (ii=0; ii < Oparams.parnum; ii++)
			{
			  rz[ii] = zfact*(rz[ii]+Lz2)-Lz2;
			  rz[ii] += OprogStatus.rzup;
			  vz[ii] += OprogStatus.vztap; 
			}
		      rebuildLinkedList();
		      MD_DEBUG3(distanza(996, 798));
		      rebuildCalendar();
		      if (OprogStatus.storerate > 0.0)
			ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
		      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
		      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
		    }
		}
	      else if ((Oparams.time - OprogStatus.quenchend)  < OprogStatus.taptau)
		{
		  /* se scalevelsteps = 0 allora scala ogni passo se si sta facendo il 
		     tapping */
		  OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
		  calcKVz();
		  MD_DEBUG4(printf("SCALVEL #%lld Vz: %.15f\n", (long long int) Oparams.curStep,Vz));
		  scalevels(Oparams.T, K, Vz);
		  rebuildLinkedList();
		  rebuildCalendar();
		  if (OprogStatus.storerate > 0.0)
		    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
		  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
		  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
		}
	      else
		{
		  /* start quench (-1  significa che il quench è iniziato) */
		  OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
		  OprogStatus.quenchend = -1;
		}
	      ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	    }
	  else if (OprogStatus.scalevel)
	    {
	      OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	      calcKVz();
	      MD_DEBUG2(printf("[TAPTAU < 0] SCALVEL #%lld Vz: %.15f\n", (long long int)Oparams.curStep,Vz));
	      scalevels(Oparams.T, K, Vz);
	      rebuildCalendar();
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	    }
#if 0
	  else if (2.0*K/(3.0*Oparams.parnum-3.0)>Oparams.T)
	    {
	      UpdateSystem();
	      calcKVz();
	      scalevels(Oparams.T, K, Vz);
	    }
#endif
	}
      if (OprogStatus.maxquench && OprogStatus.numquench == OprogStatus.maxquench)
	ENDSIM = 1;
#else
      else if (evIdB == ATOM_LIMIT + 9)
	{
	  if (OprogStatus.scalevel)
	    {
	      UpdateSystem();
	      OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	      MD_DEBUG2(printf("[TAPTAU < 0] SCALVEL #%lld Vz: %.15f\n", 
			       (long long int)Oparams.curStep,Vz));
	      K = 0.0;
	      for (i = 0; i < Oparams.parnumA; i++)
		{
		  K += Oparams.m[0]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
		}
	      for (i = Oparams.parnumA; i < Oparams.parnum; i++)
		{
		  K += Oparams.m[1]*(Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]));
		}
	      K *= 0.5;
	      scalevels(Oparams.T, K);
	      rebuildCalendar();
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	      if (OprogStatus.rescaleTime > 0)
		ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	      else
		OprogStatus.scalevel = 0;
	    }
	}
#endif
      if (OprogStatus.endtime > 0 && Oparams.time > OprogStatus.endtime)
	ENDSIM = 1;
#if 1 && defined(MD_SQWELL) && defined(MD_BONDCORR)
      corr3 = corr1 = corr2 = 0;
      for (i=0; i < Oparams.parnum; i++)
	{
	  if (numbonds0[i]==2)
	    {
	      cc = 0;
	      for (j=0; j < numbonds0[i]; j++)
		    {
		      if (bound(i,bonds0[i][j]))
			cc++;
		    }
	      if (cc>0)
		{
		  corr2++;
    		}
	      if (cc==0 && lastbreak1[i]>0.0)
		{
		  if (lastbreak1[i]>0.0 
		      && Oparams.time - lastbreak1[i] <
		      10.0*Sqr(Oparams.delta[0][0])/(Oparams.T*Oparams.Dt)) 
		    {
		      bondhist[1]++;
		    }
		  else
		    bondhist[0]+=2;
		  lastbreak1[i]=-1.0;
		}
	      if (cc==1 && lastbreak1[i] != -1.0)
		{
		  lastbreak1[i]=Oparams.time;
		}
	    }
	  if (numbonds[i]==3)
	    {
	      cc = 0;
	      for (j=0; j < numbonds0[i]; j++)
		    {
		      if (bound(i,bonds0[i][j]))
			cc++;
		    }
	      if (cc==0 && lastbreak2[i]>0.0)
		{
		  if (lastbreak2[i]>0.0 
		      && Oparams.time - lastbreak2[i] >
		      2.0*exp(2/Oparams.T)*Sqr(Oparams.delta[0][0])/(Oparams.T*Oparams.Dt)) 
		    {
		      bondhist[2]++;
		    }
		  else
		    bondhist[0]+=3;
		  lastbreak2[i]=-1.0;
		}
	      if (cc>0)
		{
		  corr3++;
		}
	      
	      if (cc==2 && lastbreak2[i] != -1.0)
		{
		  lastbreak2[i]=Oparams.time;
		}
	    }
	  if (numbonds0[i]==1 && numbonds0[bonds0[i][0]]==1)
	    {
	      if (numbonds[i]==1)	
		corr1++;
	    }	
	}
#if 0
      sprintf(fileop2 ,"BondCorrFuncB1.dat");
      /* store conf */
      strcpy(fileop, absTmpAsciiHD(fileop2));
      if ( (bof = fopenMPI(fileop, "a")) == NULL)
	{
	  mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	  exit(-1);
	}
      fprintf(bof,"%.15f %.15f\n", Oparams.time, corr1);
      fclose(bof);
      sprintf(fileop2 ,"BondCorrFuncB2.dat");
      /* store conf */
      strcpy(fileop, absTmpAsciiHD(fileop2));
      if ( (bof = fopenMPI(fileop, "a")) == NULL)
	{
	  mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	  exit(-1);
	}
      fprintf(bof,"%.15f %.15f\n", Oparams.time, corr2);
      fclose(bof);
#endif
#endif
#if 0
      if (Oparams.curStep == Oparams.totStep)
	{
	  printf(" **** end of dynamics **** \n");
	  printf(" final colliding pair (%d,%d)\n", i, j);
	  /* ** CHECK FOR PARTICLE OVERLAPS ** */
	  UpdateSystem();
	  /*check (Oparams.sigma, &overlap, &K, &V);*/
	  if ( overlap ) 
	    {
	      printf(" particle overlap in final configuration\n");
	    }
	}
#endif
#if 0
      if ( ( (OprogStatus.CMreset > 0) &&
	     /* 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) */
	   (Oparams.curStep == OprogStatus.CMreset) )
	|| ( (OprogStatus.CMreset < 0) &&
	     /* 24/3/99 CHG:((Oparams.curStep % OprogStatus.CMreset) == 0)) */
	  (Oparams.curStep % (-OprogStatus.CMreset) == 0) )  ) 
	  resetCM(Oparams.parnum);

      /* Update the integral of the pressure tensor */
      updateDQ(tij);
      updatePE(Oparams.parnum);
#endif
    }
}
