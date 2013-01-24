#include<mdsimul.h>
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
#if defined(MPI)
extern int my_rank;
extern int numOfProcs; /* number of processeses in a communicator */
extern int *equilibrated;
#endif 
#ifdef MD_HSVISCO
void calcT(void);
#endif
long long int numcoll=0;
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
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx, Mtot, Mred[2][2], invmA, invmB, DQxyOld, DQyzOld, DQzxOld;
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
double DrSq = 0.0; 
const double timbig = 1E12;
/* used by linked list routines */
#ifdef MD_GRAVITY
extern double g2, mgA, mgB;
#endif
/* double *lastcol;*/
double *treetime, *atomTime;
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
  printf(">>> numcoll:%lld\n", numcoll);
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
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, Txyold, Tyzold, Tzxold;
#endif
  double rxij, ryij, rzij, factor, invmi, invmj;
  double delpx, delpy, delpz;
  double mredl, ene;
  double sigSq, sigDeltaSq, intdistSq;
  double distSq;
  double vxij, vyij, vzij, b;
#ifdef MD_HSVISCO
  double taus;
#endif
   double bheight;
#ifdef MD_FINBARRIER
   double bhout, bhin;
#endif
#ifndef MD_MIXWDEPTH
   bheight = Oparams.bheight;
#endif
   if (i < parnumA && j < parnumA)
    {
      sigSq = Sqr(Oparams.sigma[0][0]);
#ifdef MD_MIXWDEPTH
      bheight = Oparams.bheight[0][0];
#endif
      sigDeltaSq = Sqr(Oparams.sigma[0][0]+Oparams.delta[0][0]);
      mredl = Mred[0][0];
    }
  else if (i >= parnumA && j >= parnumA)
    {
      sigSq = Sqr(Oparams.sigma[1][1]);
      sigDeltaSq = Sqr(Oparams.sigma[1][1]+Oparams.delta[1][1]);
#ifdef MD_MIXWDEPTH
      bheight = Oparams.bheight[1][1];
#endif
       mredl = Mred[1][1];
    }
  else
    {
      sigSq = Sqr(Oparams.sigma[0][1]);
      sigDeltaSq = Sqr(Oparams.sigma[0][1]+Oparams.delta[0][1]);
#ifdef MD_MIXWDEPTH
      bheight = Oparams.bheight[0][1];
#endif
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
#if defined(MD_FINBARRIER)
#ifdef MD_MIXWDEPTH
  if (i < Oparams.parnumA && j < Oparams.parnumA)  
    {
      bhin = Oparams.bhin[0][0];
      bhout= Oparams.bhout[0][0];
    }
  else if (i >= Oparams.parnumA && j >= Oparams.parnumA)
    {
      bhin = Oparams.bhin[1][1];
      bhout= Oparams.bhout[1][1];
    }
  else 
    {
      bhin = Oparams.bhin[0][1];
      bhout= Oparams.bhout[0][1];
    } 
#else
  bhout = Oparams.bhout;
  bhin  = Oparams.bhin;
#endif
#endif

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
#if defined(MD_FINBARRIER)
      if (bheight < 0)
	{
	  if (bhout >= 0.0 && Sqr(b) < 2.0*sigDeltaSq*bhout/mredl)
	    {
	      factor = -2.0*b;
	      //printf("BOND NOT BROKEN bump i=%d j=%d ata=%d atb=%d\n", i, j, ata, atb);
	    }
	  else
	    {
	      factor = -b + sqrt(Sqr(b) - 2.0*bheight/mredl);
	      //printf("BOND BROKEN bump i=%d j=%d ata=%d atb=%d\n", i, j, ata, atb);
	      remove_bond(i, j);
	      remove_bond(j, i);
	    }
	}
      else
	{
	  if (Sqr(b) < 2.0*sigDeltaSq*(bheight+bhout)/mredl)
	    {
	      factor = -2.0*b;
	    }
	  else
	    {
	      factor = -b + sqrt(Sqr(b) - 2.0*sigDeltaSq*bheight/mredl);
	      remove_bond(i, j);
	      remove_bond(j, i);
	    }
	}
#else
      if (Sqr(b) < 2.0*sigDeltaSq*bheight/mredl)
	{
	  factor = -2.0*b;
	}
      else
	{
	  factor = -b + sqrt(Sqr(b) - 2.0*sigDeltaSq*bheight/mredl);
	  remove_bond(i, j);
	  remove_bond(j, i);
	}
#endif
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
#ifdef MD_FINBARRIER
      if (bheight < 0)
	{
	  if (Sqr(b) < 2.0*sigDeltaSq*(-bheight+bhin)/mredl)
	    {
	      factor = -2.0*b;
	    }
	  else 
	    {
	      factor = -b - sqrt(Sqr(b) + 2.0*sigDeltaSq*bheight/mredl);
	      add_bond(i, j);
	      add_bond(j, i);
	    }
	}
      else
	{
	  if (bhin >= 0.0 && Sqr(b) < 2.0*sigDeltaSq*bhin/mredl)
	    {
	      factor = -2.0*b;
    	    }
	  else
	    {
	      add_bond(i, j);
	      add_bond(j, i);
	      factor = -b - sqrt(Sqr(b) + 2.0*sigDeltaSq*bheight/mredl);
	    }
	}
#else
      add_bond(i, j);
      add_bond(j, i);
      factor = -b - sqrt(Sqr(b) + 2.0*sigDeltaSq*bheight/mredl);
#endif
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
#ifdef MD_HSVISCO
  DTxy = delpx*delpy*invmi + vx[i]*delpy + delpx*vy[i];
  DTxy += delpx*delpy*invmj - vx[j]*delpy - delpx*vy[j]; 
  DTyz = delpy*delpz*invmi + vy[i]*delpz + delpy*vz[i];
  DTyz += delpy*delpz*invmj - vy[j]*delpz - delpy*vz[j];
  DTzx = delpz*delpx*invmi + vz[i]*delpx + delpz*vx[i];
  DTzx += delpz*delpx*invmj - vz[j]*delpx - delpz*vx[j];
#endif
  vx[i] = vx[i] + delpx*invmi;
  vx[j] = vx[j] - delpx*invmj;
  vy[i] = vy[i] + delpy*invmi;
  vy[j] = vy[j] - delpy*invmj;
  vz[i] = vz[i] + delpz*invmi;
  vz[j] = vz[j] - delpz*invmj;
  //calcT();
#ifdef MD_HSVISCO 
  if (OprogStatus.lastcoll!=-1)
    {
      taus = Oparams.time - OprogStatus.lastcoll;
      OprogStatus.DQTxy += OprogStatus.Txy*taus; 
      OprogStatus.DQTyz += OprogStatus.Tyz*taus;
      OprogStatus.DQTzx += OprogStatus.Tzx*taus;
      //taus = Oparams.time - OprogStatus.lastcoll;
      //printf("DQT= %f %f %f\n", OprogStatus.DQTxy, OprogStatus.DQTyz, OprogStatus.DQTzx);
      OprogStatus.DQWxy += rxij*delpy;
      OprogStatus.DQWyz += ryij*delpz;
      OprogStatus.DQWzx += rzij*delpx;
      //printf("DQW= %f %f %f\n", OprogStatus.DQWxy, OprogStatus.DQWyz, OprogStatus.DQWzx);
    }
#if 0
  Txyold = OprogStatus.Txy;
  Tyzold = OprogStatus.Tyz;
  Tzxold = OprogStatus.Tzx;
#endif
  OprogStatus.Txy += DTxy; 
  OprogStatus.Tyz += DTyz;
  OprogStatus.Tzx += DTzx;
#if 0
  printf("P STEP #%d T= %f %f %f\n", Oparams.curStep, OprogStatus.Txy, OprogStatus.Tyz, OprogStatus.Tzx);
  printf("P STEP #%d DT=%f %f %f\n", Oparams.curStep, DTxy, DTyz, DTzx);
  calcT();
  printf("D STEP #%d T= %f %f %f\n", Oparams.curStep, OprogStatus.Txy, OprogStatus.Tyz, OprogStatus.Tzx);
  printf("DT STEP #%d T= %f %f %f\n", Oparams.curStep, OprogStatus.Txy-Txyold, OprogStatus.Tyz-Tyzold, OprogStatus.Tzx-Tzxold);
#endif
#endif
#if 0
  *W = delpx * rxij + delpy * ryij + delpz * rzij;
#endif

}
#else
void bump (int i, int j, double* W)
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
  double delpx, delpy, delpz, sigSq;
  double mredl, ene;
  if (i < parnumA && j < parnumA)
    {
      sigSq = Sqr(Oparams.sigma[0][0]);
      mredl = Mred[0][0];
    }
  else if (i >= parnumA && j >= parnumA)
    {
      sigSq = Sqr(Oparams.sigma[1][1]);
      mredl = Mred[1][1];
    }
  else
    {
      sigSq = Sqr(Oparams.sigma[0][1]);
      mredl = Mred[0][1]; 
    }
  /*printf("(i:%d,j:%d sigSq:%f\n", i, j, sigSq);*/
  /*printf("mredl: %f\n", mredl);*/
  rxij = rx[i] - rx[j];
  if (fabs (rxij) > L2)
    rxij = rxij - SignR(L, rxij);
  ryij = ry[i] - ry[j];
  if (fabs (ryij) > L2)
    ryij = ryij - SignR(L, ryij);
  rzij = rz[i] - rz[j];
#if !defined(MD_GRAVITY)
  if (fabs (rzij) > L2)
    rzij = rzij - SignR(L, rzij);
#endif
#ifdef MD_GRAVITY
  factor = ( rxij * ( vx[i] - vx[j] ) +
	     ryij * ( vy[i] - vy[j] ) +
	     rzij * ( vz[i] - vz[j] ) ) / sigSq;
  /* Dissipation */
  if (!((Oparams.time - lastcol[i] < OprogStatus.tc)||
  	(Oparams.time - lastcol[j] < OprogStatus.tc)))
    factor *= mredl*(1+Oparams.partDiss);
#else
  /* Nel caso di gravita' e' intuile implementare il TC-model di Luding
   * per evitare il collasso inelastico.
   * Gli urti in tale caso sono tutti elastici. */ 
  /* SQUARE WELL: modify here */
  factor = ( rxij * ( vx[i] - vx[j] ) +
	     ryij * ( vy[i] - vy[j] ) +
	     rzij * ( vz[i] - vz[j] ) ) / sigSq;
  factor *= mredl*2; /*(1+Oparams.partDiss);*/
#endif
  delpx = - factor * rxij;
  delpy = - factor * ryij;
  delpz = - factor * rzij;
  invmi = (i<Oparams.parnumA)?invmA:invmB;
  invmj = (j<Oparams.parnumA)?invmA:invmB;
  ene= (Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])+
	    Sqr(vx[j])+Sqr(vy[j])+Sqr(vz[j])); 
  vx[i] = vx[i] + delpx*invmi;
  vx[j] = vx[j] - delpx*invmj;
  vy[i] = vy[i] + delpy*invmi;
  vy[j] = vy[j] - delpy*invmj;
  vz[i] = vz[i] + delpz*invmi;
  vz[j] = vz[j] - delpz*invmj;
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
	OprogStatus.equilibrated = (DrSqTot>OprogStatus.eqlevel*Oparams.sigma[0][0]?1:0);
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
void UpdateAtom(int i)
{
  double ti;
  ti = Oparams.time - atomTime[i];
  rx[i] += vx[i]*ti;
  ry[i] += vy[i]*ti;
#if defined(MD_GRAVITY)
  rz[i] += vz[i]*ti - g2*Sqr(ti);
  vz[i] += -Oparams.ggrav*ti;
#else
  rz[i] += vz[i]*ti;
#endif
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
void rebuildCalendar(void);
extern void ScheduleEventBarr (int idA, int idB, int idC, double tEvent); 

void PredictEvent (int na, int nb) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq, dr[NDIM], dv[NDIM], shift[NDIM], tm[NDIM],
  b, d, t, tInt, vv;
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
			  sigSq = Sqr(Oparams.sigma[0][0]);
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
		      if (na < parnumA && n < parnumA)
			sigSq = Sqr(Oparams.sigma[0][0]);
		      else if (na >= parnumA && n >= parnumA)
			sigSq = Sqr(Oparams.sigma[1][1]);
		      else
			sigSq = Sqr(Oparams.sigma[0][1]);
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
		      if (!bound(n, na) && sigSq != sigDeltaSq)
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
			  //printf("sigDeltaSq: %f sigSq: %f\n", sigDeltaSq, sigSq);
			  if (sigDeltaSq == sigSq)
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
			  if (collCode == MD_EVENT_NONE && sigSq != sigDeltaSq)
			    {
			      d = Sqr (b) - vv * (distSq - sigDeltaSq);
			      if (d > 0.0)
				{
				  t = ( sqrt (d) - b) / vv;
				  if (t > 0) //|| (t < 0 && distSq > sigDeltaSq))
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
		      if (b < 0.0) 
			{
			  vv = Sqr(dv[0]) + Sqr (dv[1]) + Sqr (dv[2]);
			  d = Sqr (b) - vv * 
			    (Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2]) - sigSq);
#if 0
			  if (OprogStatus.quenchend > 0.0 && Oparams.curStep > 700000)
			    printf("dist:%.15G\n", 
				   (Sqr (dr[0]) + Sqr (dr[1]) + Sqr(dr[2])) - sigSq );
#endif
    			  if (d >= 0.) 
			    {
			      t = - (sqrt (d) + b) / vv;
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
			      ScheduleEvent (na, n, Oparams.time + t);
			      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
			    } 
			}
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
#if defined(MD_SQWELL)||defined(MD_INFBARRIER)
  /* i primi due bit sono il tipo di event (uscit buca, entrata buca, collisione con core 
   * mentre nei bit restanti c'e' la particella con cui tale evento e' avvenuto */
  bump(evIdA, evIdB, &W, evIdC);
#else
  bump(evIdA, evIdB, &W);
#endif
  numcoll++;
#ifdef MD_HSVISCO
  OprogStatus.lastcoll = Oparams.time;
#endif
  /*printf("qui time: %.15f\n", Oparams.time);*/
#ifdef MD_GRAVITY
  lastcol[evIdA] = lastcol[evIdB] = Oparams.time;
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
extern void comvel_brown (COORD_TYPE temp, COORD_TYPE *m);

void velsBrown(double T)
{
  comvel_brown(T, Oparams.m); 
}
extern double gauss(void);
extern double ranf(void);
#ifdef MD_MICRO_LANG
void random_direction(double *nx, double *ny, double *nz)
{
  double xisq, xi1, xi2, xi, nsq, norm;
  xisq = 1.0;
  while (xisq >= 1.0)
    {
      xi1  = ranf() * 2.0 - 1.0;
      xi2  = ranf() * 2.0 - 1.0;
      xisq = xi1 * xi1 + xi2 * xi2;
    }
  xi = sqrt (fabs(1.0 - xisq));
  *nx = 2.0 * xi1 * xi;
  *ny = 2.0 * xi2 * xi;
  *nz = 1.0 - 2.0 * xisq;

  /* Renormalize */
  nsq   = (*nx)*(*nx) + (*ny)*(*ny) + (*nz)*(*nz);
  norm  = sqrt(fabs(nsq));
  *nx    = *nx / norm;
  *ny    = *ny / norm;
  *nz    = *nz / norm;
}


void velsMicroLang(double T, double xi)
{
   double c1, c2, M, n, gam, vpx, vpy, vpz, m, kTm, nx, ny, nz;
   double mredl, b, vxij, vyij, vzij, delpx, delpy, delpz, factor;	
   int i;
   n = 1.0/Oparams.Dt;
   //printf("c1=%f c2=%f T=%.15G m=%.15G\n", c1, c2, T, m);
   //c2 = kTm*2*m / (M+m);
   for (i = 0; i < Oparams.parnum; i++)
    {
       if (i == 0)
	 {
	   M = Oparams.m[0];
	   gam = Oparams.xi*M;
	   /* il 3 deriva dal fatto che bisogna mediare su metà angolo solido!*/
	   m = Oparams.Dt*gam/2.0;// *(3.0 / 2.0); 
	   mredl = m*M/(m+M);
	   kTm = sqrt(Oparams.T / m);
	   c1 = (M - m)/(M+m);
	   c2 = kTm*2*m / (M+m);
	 }
       else if (i == Oparams.parnumA)
	 {
	   M = Oparams.m[1];
	   /* il 3 deriva dal fatto che bisogna mediare su metà angolo solido!*/
	   //mredl = m*M/(m+M);
	   gam = Oparams.xi*M;
	   m = Oparams.Dt*gam/2.0;// *(3.0 / 2.0); 
	   kTm = sqrt(Oparams.T / m);
	   c1 = (M - m)/(M+m);
	   c2 = kTm*2*m / (M+m);
	 }
#if 1
       vpx = gauss();
       vpy = gauss();
       vpz = gauss();
#else
       vpx = kTm*gauss();
       vpy = kTm*gauss();
       vpz = kTm*gauss();
#endif
       //random_direction(&nx, &ny, &nz);	
#if 1
       vx[i] = c1*vx[i] + c2*vpx;
       vy[i] = c1*vy[i] + c2*vpy;
       vz[i] = c1*vz[i] + c2*vpz;
#else
       vxij = vx[i] - vpx;
       vyij = vy[i] - vpy;
       vzij = vz[i] - vpz;
#if 0
       /* se n è sempre parallelo alla velocità
       * relativa allora torniamo all'implementazione
       * originaria errata */
       {
         double norm;
	 norm = sqrt(Sqr(vxij)+Sqr(vyij)+Sqr(vzij));
	 nx = vxij/norm;
	 ny = vyij/norm;
	 nz = vzij/norm;
       }
#endif
#if 0
       b = nx*vxij + ny*vyij + nz*vzij; 
#if 0
       if (b > 0.0)
	 {
	    nx = -nx;
	    ny = -ny;
	    nz = -nz;
	    b = -b;
	 }
#endif
       //b = sqrt(Sqr(vxij)+Sqr(vyij)+Sqr(vzij));
       factor = -2.0*b;
       factor *= mredl;
       delpx = factor*nx;
       delpy = factor*ny;
       delpz = factor*nz;
#if 0
       vx[i] += delpx/M;
       vy[i] += delpy/M;
       vz[i] += delpz/M;
        
#if 0
       if (fabs(delpx) > 1.0 || fabs(delpy) > 1.0 || fabs(delpz) > 1.0)
	 {
	   printf("delp=(%.15G,%.15G,%.15G)\n", delpx, delpy, delpz);
	   printf("M=%.15G mredl=%.15G m=%.15G\n", M, mredl, m);
	   printf("v=(%.15G,%.15G,%.15G)\n", vx[i], vy[i], vz[i]);
	   printf("vp=(%.15G,%.15G,%.15G) n=%.15G,%.15G,%.15G\n", vpx, vpy, vpz,nx,ny,nz);
	   //printf("kTm=%.15G b=%.15G\n", kTm, b);
	   //printf("b=%.15G factor=%.15G\n", b,factor);
	 }
#endif
#else
	/* cosi' è come nel caso 1D, ma essendo qui in 3D cio' equivale ad assumere che il versore
	* n sia sempre parallelo alla velocità relativa vij! */
	vx[i] -= vxij*2.0*mredl/M;
	vy[i] -= vyij*2.0*mredl/M;
	vz[i] -= vzij*2.0*mredl/M;
#endif
       //printf("delp=(%.15G, %.15G, %.15G)\n", delpx, delpy, delpz);
#endif
       //printf("vpx=%.15G,%.15G,%.15G kTm=%.15G\n", vpx, vpy, vpz, kTm);
       //printf("(%.15G, %.15G, %.15G)\n", vx[i], vy[i], vz[i]);
#endif
    }
}
#elif defined(MD_FULL_LANG)
void velsFullLang(double T, double xi)
{
   double c1, c2, M, n, gam, vxp, vyp, vzp, m, kTm, nx, ny, nz;
   double mredl, b, vxij, vyij, vzij, delpx, delpy, delpz, factor;	
   double exp2gamt, expgamt;
   int i;
   //printf("c1=%f c2=%f T=%.15G m=%.15G\n", c1, c2, T, m);
   //c2 = kTm*2*m / (M+m);
   for (i = 0; i < Oparams.parnum; i++)
    {
       if (i == 0)
	 {
	   c1=exp(-Oparams.xi*Oparams.Dt);
	   kTm = Oparams.T / Oparams.m[0];
	   c2 = sqrt(kTm*(1 - Sqr(c1)));
           //printf("kTm=%.15G c1=%.15Gi c2:%.15G\n", kTm, c1, c2);
	 }
       else if (i == Oparams.parnumA)
	 {
	   c1=exp(-Oparams.xi*Oparams.Dt);
	   kTm = Oparams.T / Oparams.m[1];
	   c2 = sqrt(kTm*(1 - Sqr(c1)));
	 }
       //random_direction(&nx, &ny, &nz);	
       vx[i] = c1*vx[i] + c2*gauss();
       vy[i] = c1*vy[i] + c2*gauss();
       vz[i] = c1*vz[i] + c2*gauss();
}
} 
#endif
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
extern void InitEventList (void);

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
#ifdef MD_SAVEFRA
extern void error_on_writing(FILE* f, char* filename, char *where, char *command);
void save_fra(void)
{
  char fileop[1024], fileop2[1024], fileop3[1024];
  double rat[5][3];
  double rcm[3];	
  int i, ret, rerr;
  FILE* f;
  sprintf(fileop2 ,"Cnf-%d-%d", 
	  OprogStatus.KK, OprogStatus.JJ);
  /* store conf */
  strcpy(fileop, absTmpAsciiHD(fileop2));
  if ( (f = fopenMPI(fileop, "w")) == NULL)
    {
      mdPrintf(STD, "Error saving store file!\n", NULL);
      exit(-1);
    }
  rerr = 0;
#ifdef MD_BIG_DT
  ret = fprintf(f, "%lld 0 %d %d 0\n", 
		(long long int)((Oparams.time+OprogStatus.refTime)*1000.0/Oparams.Dt), 
		Oparams.parnum, Oparams.parnum-Oparams.parnumA);
#else
  ret = fprintf(f, "%lld 0 %d %d 0\n", (long long int)(Oparams.time*1000.0/Oparams.Dt), 
		Oparams.parnum, Oparams.parnum-Oparams.parnumA); 
#endif
  rerr |= (ret < 0)?1:0;  
  ret = fprintf(f, "%.15G %.15G %.15G 0 0 %.15G\n", L, L, L, Oparams.Dt/1000.0);
  rerr |= (ret < 0)?1:0;  
  for (i = Oparams.parnumA; i < Oparams.parnum; i++)
    {
      //printf("i=%d\n",i);
      rcm[0] = rx[i];//+L*OprogStatus.DR[i][0];
      rcm[1] = ry[i];//+L*OprogStatus.DR[i][1];
      rcm[2] = rz[i];//+L*OprogStatus.DR[i][2];
#if 0
      BuildAtomPos(i, rcm, R[i], rat);
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[1][0], rat[1][1], rat[1][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[2][0], rat[2][1], rat[2][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[0][0], rat[0][1], rat[0][2]);
#endif
      ret = fprintf(f, "%.15G %.15G %.15G\n", rcm[0], rcm[1], rcm[2]);
      rerr |= (ret < 0)?1:0;  
    }
  for (i = 0; i < Oparams.parnumA; i++)
    {
      rcm[0] = rx[i];//+L*OprogStatus.DR[i][0];
      rcm[1] = ry[i];//+L*OprogStatus.DR[i][1];
      rcm[2] = rz[i];//+L*OprogStatus.DR[i][2];
#if 0
      BuildAtomPos(i, rcm, R[i], rat);
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[1][0], rat[1][1], rat[1][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[2][0], rat[2][1], rat[2][2]);
      rerr |= (ret < 0)?1:0;  
      ret = fprintf(f, "%.15G %.15G %.15G\n", rat[0][0], rat[0][1], rat[0][2]);
#endif
      ret = fprintf(f, "%.15G %.15G %.15G\n", rcm[0], rcm[1], rcm[2]);
      rerr |= (ret < 0)?1:0;  
    }
  for (i = 0; i < Oparams.parnum; i++)
    {
      ret = fprintf(f, "%.15G %.15G %.15G\n", vx[i], vy[i], vz[i]);
      rerr |= (ret < 0)?1:0;  
    }

#ifdef MD_SAVEFRA_COMPRESSED
#ifdef MD_MAC
  sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
#else
  sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
#ifndef MD_NO_SYSTEM
  system(fileop3);
#endif
#endif
  if (rerr)
    error_on_writing(f, fileop, "save_fra", "fprintf");
  fclose(f);
}
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
  while (1)
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
#ifndef MD_SAVEFRA
	  sprintf(fileop2 ,"Store-%d-%d", 
		  OprogStatus.KK, OprogStatus.JJ);
	  /* store conf */
	  strcpy(fileop, absTmpAsciiHD(fileop2));
	  if ( (bf = fopenMPI(fileop, "w")) == NULL)
	    {
	      mdPrintf(STD, "Errore nella fopen in saveBakAscii!\n", NULL);
	      exit(-1);
	    }
#endif
	  UpdateSystem();
#ifdef MD_SAVEFRA
	  save_fra();
#else
	  writeAsciiPars(bf, opro_ascii);
	  fprintf(bf, sepStr);
	  writeAsciiPars(bf, opar_ascii);
	  fprintf(bf, sepStr);
	  writeAllCor(bf);
	  fclose(bf);
#ifdef MPI
#ifdef MD_MAC
          sprintf(fileop3, "/usr/bin/gzip -f %s_R%d", fileop, my_rank);
#else
          sprintf(fileop3, "/bin/gzip -f %s_R%d", fileop, my_rank);
#endif
#else 
#ifdef MD_MAC
          sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
#else
          sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
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
	  if (OprogStatus.brownian)
	    {
#ifdef MD_HSVISCO
	      double taus, Vol;
	      if (OprogStatus.lastcoll!=-1)
		{
		  /* notare che nel caso di dinamica browniana
		   * lastcoll è in generale l'ultima collisione o tra due particelle
		   * o tra le particelle e il fluido (reset delle velocità)*/
		  taus = Oparams.time - OprogStatus.lastcoll; 
		  OprogStatus.DQTxy += taus * OprogStatus.Txy;
		  OprogStatus.DQTyz += taus * OprogStatus.Tyz;
		  OprogStatus.DQTzx += taus * OprogStatus.Tzx;
#if 1
		  DQxyOld = OprogStatus.DQxy;
		  DQyzOld = OprogStatus.DQyz;
		  DQzxOld = OprogStatus.DQzx;
#endif
		  OprogStatus.DQxy = OprogStatus.DQTxy + OprogStatus.DQWxy;
		  OprogStatus.DQyz = OprogStatus.DQTyz + OprogStatus.DQWyz;
		  OprogStatus.DQzx = OprogStatus.DQTzx + OprogStatus.DQWzx;
		  Vol = L*L*L;
		  OprogStatus.DQxy /= Vol;
		  OprogStatus.DQyz /= Vol;
		  OprogStatus.DQzx /= Vol;
#if 1
		  Pxy = (OprogStatus.DQxy - DQxyOld)/Oparams.Dt;
		  Pyz = (OprogStatus.DQyz - DQyzOld)/Oparams.Dt;
		  Pzx = (OprogStatus.DQzx - DQzxOld)/Oparams.Dt;
#endif
	      	  //printf("DQ= %f %f %f\n", OprogStatus.DQxy, OprogStatus.DQyz, OprogStatus.DQzx);
		  OprogStatus.lastcoll = Oparams.time;
		}
#endif
#if defined(MD_MICRO_LANG)
	      velsMicroLang(Oparams.T, Oparams.xi);
#elif defined(MD_FULL_LANG)
	      velsFullLang(Oparams.T, Oparams.xi);

#else
	      velsBrown(Oparams.T);
#endif
#ifdef MD_HSVISCO
	      calcT();
#endif
	      rebuildCalendar();
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	    }
#ifdef MD_HSVISCO
	  else
	    {
	      double Vol;
#if 1
	      DQxyOld = OprogStatus.DQxy;
	      DQyzOld = OprogStatus.DQyz;
	      DQzxOld = OprogStatus.DQzx;
#endif
    	      OprogStatus.DQxy = OprogStatus.DQTxy + OprogStatus.DQWxy;
	      OprogStatus.DQyz = OprogStatus.DQTyz + OprogStatus.DQWyz;
	      OprogStatus.DQzx = OprogStatus.DQTzx + OprogStatus.DQWzx;
	      Vol = L*L*L;
	      OprogStatus.DQxy /= Vol;
	      OprogStatus.DQyz /= Vol;
	      OprogStatus.DQzx /= Vol;
#if 1
    	      Pxy = (OprogStatus.DQxy - DQxyOld)/Oparams.Dt;
	      Pyz = (OprogStatus.DQyz - DQyzOld)/Oparams.Dt;
	      Pzx = (OprogStatus.DQzx - DQzxOld)/Oparams.Dt;
#endif
	      //printf("STEP #%d DQ= %f %f %f\n", Oparams.curStep, OprogStatus.DQxy, OprogStatus.DQyz, OprogStatus.DQzx);
	    }
#endif
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
