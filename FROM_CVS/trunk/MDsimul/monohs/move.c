#include<mdsimul.h>
#define SIMUL
#define SignR(x,y) (((y) >= 0) ? (x) : (- (x)))
int checkz(char *msg)
{
  int ii;
  for (ii=0; ii < Oparams.parnum; ii++)
    {
      if ((rz[ii]+Lz*0.5-Oparams.sigma/2.0)<0. &&  
	fabs(rz[ii]+Lz*0.5-Oparams.sigma/2.0) > 0.1)
	{
	  printf("[%s]ii=%d diff:%.15f\n", msg,
		 ii, rz[ii]+Lz*0.5-Oparams.sigma/2.0);
	  printf("sotto more in checkz!!!\n");
	  return 1;
	}
    }
  return 0;
}

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
  Wxy, Wyz, Wzx, Pxx, Pyy, Pzz, Pxy, Pyz, Pzx;
/*  Patxy, Patyz, Patzx, Patxx, Patyy, Patzz,
    T1myz, T1mzx, T1mxx, T1myy, T1mzz;  */
double DrSq = 0.0; 
const double timbig = 1E12;
/* used by linked list routines */
#ifdef MD_GRAVITY
extern double g2, mg;
#endif
double *treetime, *atomTime;
int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
  cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
int evIdA, evIdB;
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
void calcKVz(void)
{
  int i;
  Vz = 0.0;
  for (i = 0; i < Oparams.parnum; i++)
    {
      Vz += vz[i];
    }
  Vz /= ((double)Oparams.parnum);
  K = 0.0;
  for (i=0; i < Oparams.parnum; i++)
    {
      K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
    }
  K *= Oparams.m * 0.5;
}
void calccmz(void);
void calcRho(void);

void outputSummary(void)
{
  FILE *f;
  /* mettere qualcosa qui */
#if 1
  printf("time=%.15f\n", OprogStatus.time);
  printf("K= %.15f V=%.15f T=%.15f Vz: %f\n", K, V, 
	 (2.0*K/(3.0*Oparams.parnum-3.0)), Vz);
#endif
  f = fopenMPI(MD_HD_MIS "T.dat", "a");
  fprintf(f, "%.15f %.15f\n", OprogStatus.time, (2.0*K/(3.0*Oparams.parnum-3.0)));
  fclose(f);
  f = fopenMPI(MD_HD_MIS "Vz2.dat", "a");
  fprintf(f, "%.15f %.15f\n", OprogStatus.time, Sqr(Vz));
  fclose(f);
  if (OprogStatus.numquench==0)
    {
      calcRho();
      calccmz();
      f = fopenMPI(MD_HD_MIS "rho_ist.dat", "a");
      fprintf(f, "%.15f %.15f\n", OprogStatus.time, rho);
      fclose(f);
      f = fopenMPI(MD_HD_MIS "rcmz_ist.dat", "a");
      fprintf(f, "%.15f %.15f\n", OprogStatus.time, rcmz);
      fclose(f);
    }
}
void scalevels(double temp, double K, double Vz)
{
  int i; 
  double sf, VVx, VVy, VVz;

  sf = sqrt( ( (3.0*((double)Oparams.parnum)-3.0) * temp ) / (2.0*K) );

  VVx = VVy = VVz = 0.0;
  for (i = 0; i < Oparams.parnum; i++)
    {
      vx[i] *= sf;
      vy[i] *= sf;
      vz[i] = (vz[i] - Vz)*sf + Vz;
      VVx += vx[i];
      VVy += vy[i];
      VVz += vz[i];
      /* scala anche i tempi di collisione! */
    } 
  VVx = VVx / (double) Oparams.parnum;
  VVy = VVy / (double) Oparams.parnum;
  VVz = VVz / (double) Oparams.parnum;
  for (i = 0; i < Oparams.parnum; i++)
    {
      vx[i] -= VVx;
      vy[i] -= VVy;
      vz[i] = vz[i] - VVz + Vz;
    }
  MD_DEBUG2(printf("sf: %.15f temp: %f K: %f Vz: %.15f minvz:%.15G\n", sf, temp, K, Vz));
}

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
void check ( double sigma, int *overlap, double *K, double *V)
{
  /* *******************************************************************
   ** TESTS FOR PAIR OVERLAPS AND CALCULATES KINETIC ENERGY.        **
   *******************************************************************
   */
   int  i, j;
   double mg, rxi, ryi, rzi, rxij, ryij, rzij, rijSq, sigSq, rij;
   double tol=OprogStatus.overlaptol;
   
   printf("overlaptol: %f\n", OprogStatus.overlaptol);
   sigSq  = Sqr(Oparams.sigma);
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
#ifdef MD_GRAVITY
   mg = Oparams.m * Oparams.ggrav;
#endif
   for (i = 0; i < Oparams.parnum; i++)
     {
       *K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]);
#ifdef MD_GRAVITY
       *V += rz[i];
#endif
     }
   
   *K = 0.5 * (*K) * Oparams.m;
#ifdef MD_GRAVITY
   *V *= mg; 
#endif
}

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
  double rxij, ryij, rzij, factor;
  double delvx, delvy, delvz, sigSq;

  sigSq = Sqr(Oparams.sigma);
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
  factor = ( rxij * ( vx[i] - vx[j] ) +
	     ryij * ( vy[i] - vy[j] ) +
	     rzij * ( vz[i] - vz[j] ) ) / sigSq;
  /* Dissipation */
  if (!((OprogStatus.time - lastcol[i] < OprogStatus.tc)||
      (OprogStatus.time - lastcol[j] < OprogStatus.tc)))
    factor *= (1+Oparams.partDiss)*0.5;
  delvx = - factor * rxij;
  delvy = - factor * ryij;
  delvz = - factor * rzij;
  vx[i] = vx[i] + delvx;
  vx[j] = vx[j] - delvx;
  vy[i] = vy[i] + delvy;
  vy[j] = vy[j] - delvy;
  vz[i] = vz[i] + delvz;
  vz[j] = vz[j] - delvz;
  /* TO CHECK: il viriale ha senso solo se non c'� la gravit� */
  *W = delvx * rxij + delvy * ryij + delvz * rzij;
}
void calccmz(void)
{
  int i;
  rcmz = 0.0;
  
  for(i = 0; i < Oparams.parnum; i++)
    {
      rcmz += rz[i];
    }
  
  rcmz /= Oparams.parnum;

}
extern double *rhoz;
void calcRho(void)
{
  int jZ, jY, jX, n, npart, jZmax, i;
  double rhohcp, hhcp, Lz2;
  double dia = Oparams.sigma;
  Lz2 = Lz*0.5;
  /* hhcp � l'altezza delle particelle se fossero close-packed */
  rhohcp =0.7405*24/(4*pi*dia*dia*dia) ; 
  hhcp = Oparams.parnum/(rhohcp*Sqr(L));

  if (OprogStatus.rhobh <= 0)
    hhcp = Oparams.parnum/(rhohcp*Sqr(L))+OprogStatus.rhobh;
  else
    hhcp = (Oparams.parnum/(rhohcp*Sqr(L)))*OprogStatus.rhobh;
  /* Se rhobh > 0 allora l'altezza per il calcolo della densit� �:
   * h_closepacking * rhobh */
  jZmax = (int) ((Lz / (OprogStatus.extraLz + Lz)) * cellsz)+1;  
  jZmax = (jZmax < cellsz)?jZmax:cellsz;
  MD_DEBUG2(printf("cellsz: %d iZmax: %d hhcp: %.15f\n", cellsz, jZmax, hhcp));
  npart = 0;
  for (i = 0; i < jZmax; i++)
    {
      rhoz[i] = 0.0;
    }
  for (jZ = 0; jZ < jZmax; jZ++)
    for (jX = 0; jX < cellsx; jX++)
      for (jY = 0; jY < cellsy; jY++)
	{
	  n = (jZ *cellsy + jY) * cellsx + jX + Oparams.parnum;
	  for (n = cellList[n]; n > -1; n = cellList[n]) 
	    {
	      if (rz[n] + Lz2 + Oparams.sigma*0.5 < hhcp)
		npart++;
	      if (jZ < jZmax)
		rhoz[jZ/2]+= 1.0;
	    }
	}
  rho = 0.523598776*((double)npart)/(L*L*hhcp);
  for (i = 0; i < jZmax/2; i++)
    {
      rhoz[i] *= 0.523598776/(2.0*L*L*Oparams.rcut);
    }
}

void save_rho(void)
{
  int jZmax, i;
  FILE *f;
  f = fopenMPI(MD_HD_MIS "rho.dat", "a");
  fprintf(f, "%d %.15f\n", OprogStatus.numquench, rho);
  jZmax = (int) ((Lz / (OprogStatus.extraLz + Lz)) * cellsz)+1;  
  jZmax = (jZmax < cellsz)?jZmax:cellsz;
  fclose(f);
  f = fopenMPI(MD_HD_MIS "rhoz.dat", "a");
  fprintf(f, "%d", OprogStatus.numquench);
  for (i = 0; i < jZmax; i++)
    fprintf(f, " %.15f", rhoz[i]);
  fprintf(f, "\n");
  fclose(f);
}
void save_rzcm(void)
{
  FILE* f;
  f = fopenMPI(MD_HD_MIS "rcmz.dat", "a");
  fprintf(f, "%d %.15f\n", OprogStatus.numquench, rcmz);
  fclose(f);
}
void calcObserv(void)
{
  int i;
  K = 0.0;
  V = 0.0;
  Vz = 0.0;
 
   /* Bisogna considerare le velocit� rispetto al centro di massa! */
  for (i=0; i < Oparams.parnum; i++)
    {
      Vz += vz[i]; 
    }
  Vz /= Oparams.parnum;
  calccmz();

  for (i = 0; i < Oparams.parnum; i++)
    {
      K += Sqr(vx[i]) + Sqr(vy[i]) + Sqr(vz[i]-Vz);
#ifdef MD_GRAVITY
      V += rz[i];
#endif
    }
  K *= Oparams.m * 0.5;
#ifdef MD_GRAVITY
  V *= mg;
#endif
}
extern double *treeTime;
void UpdateAtom(int i)
{
  double ti;
  ti = OprogStatus.time - atomTime[i];
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
  atomTime[i] = OprogStatus.time;
}
void UpdateSystem(void)
{
  int i;
  /* porta tutte le particelle allo stesso tempo */
  for (i=0; i < Oparams.parnum; i++)
    {
      UpdateAtom(i);
    }
}
void AdjustLastcol(void)
{
  int i;
  for (i=0; i <  Oparams.parnum; i++)
    lastcol[i] -= OprogStatus.time;
}
void rebuildCalendar(void);
void PredictEvent (int na, int nb) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double Lzx, sigSq, dr[NDIM], dv[NDIM], shift[NDIM], tm[NDIM], cells[NDIM],
  b, d, t, tInt, vv, h1, h2, hh1;
  int cellRangeT[2 * NDIM], signDir[NDIM], evCode,
  iX, iY, iZ, jX, jY, jZ, k, n;

  MD_DEBUG(printf("PredictEvent: %d,%d\n", na, nb));
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravit� � diretta lungo z negativo */ 
#ifdef MD_GRAVITY
  Lzx = OprogStatus.extraLz + Lz;

  /* NOTA: Il muro inferiore � posto a Lz / 2 cio� sul fondo del box,
   * in questo modo la base di ogni sfera poggia esattamente 
   * sul fondo della scatola */
  if (inCell[2][na] == 0)
    {
      hh1 =  vz[na] * vz[na] + 2.0 * Oparams.ggrav *
	(rz[na] + Lz2);
      h1 = hh1 -  Oparams.ggrav * Oparams.sigma;
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
      /* h1 � il Discriminante dell'intersezione con la faccia 
       * inferiore della cella lungo z, h2 di quella superiore,
       * per cui se vz > 0 e h2 > 0 allora si la particella attraverser�
       * la faccia superiore poich� h2 < h1 => t_h2 < t_h1
       * (si noti che la soluzione con tempo positivo se vz > 0 � quella
       * con il + nella formula per la risoluz. di un'eq. di secondo grado */
      h2 = hh1 - 2.0 * Oparams.ggrav * Lzx  / cellsz;
      if (h2 > 0.0) 
	{
	  tm[2] =  (vz[na] - sqrt (h2)) / Oparams.ggrav;
	  signDir[2] = 0;/* signDir = 0 vuol dire che la direzione �
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
      printf("tm[%d]<0 step %d na=%d\n", k, Oparams.curStep, na);
      printf("rz:%f diff:%f\n", rz[na], rz[na]+Lz2);
      printf("h1:%f hh1:%f vz:%f cellz:%d\n", h1, hh1, vz[na], inCell[2][na]);
      printf("Cells(%d,%d,%d)\n", inCell[0][na], inCell[1][na], inCell[2][na]);
      printf("signDir[0]:%d signDir[1]: %d signDir[2]: %d\n", signDir[0], signDir[1],
	     signDir[2]);
      /*exit(-1);*/
      tm[k] = 0.0;
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
  		    che � anche l'unica nel presente caso */
      MD_DEBUG2(printf("wall!!! (evIdA: %d)\n", na));
    }
#endif
  MD_DEBUG(printf("schedule event [WallCrossing](%d,%d)\n", na, ATOM_LIMIT+evCode));
  ScheduleEvent (na, ATOM_LIMIT + evCode, OprogStatus.time + tm[k]);
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];
#ifdef MD_GRAVITY
  /* k = 2 : lungo z con la gravita' non ci sono condizioni periodiche */
  if (inCell[2][na] + cellRangeT[2 * 2] < 0) cellRangeT[2 * 2] = 0;
  if (inCell[2][na] + cellRangeT[2 * 2 + 1] == cellsz) cellRangeT[2 * 2 + 1] = 0;
#endif
  sigSq = Sqr(Oparams.sigma);

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
		      tInt = OprogStatus.time - atomTime[n];
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
				  printf("time:%.15f tInt:%.15f\n", OprogStatus.time,
					 tInt);
				  printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
					 Sqr(dr[2]))-1.0 );
				  printf("STEP: %d\n", Oparams.curStep);
				  printf("atomTime: %.10f \n", atomTime[n]);
				  printf("n:%d na:%d\n", n, na);
				  printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
				  /*exit(-1);*/
				  t = 0;
				}
			      ScheduleEvent (na, n, OprogStatus.time + t);
			      MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
			    } 
			}
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
  printf("timeNow-lastcol[%d]:%.15f\n", evIdA, OprogStatus.time - lastcol[evIdA]);
  printf("walldiss:%.15f vz:%.10f\n", Oparams.wallDiss, vz[evIdA]);
  printf("timeNow:%.10f lastcol:%.10f\n", OprogStatus.time, lastcol[evIdA]);
#endif
  if (OprogStatus.time - lastcol[evIdA] < OprogStatus.tc)
    {
      vz[evIdA] = -vz[evIdA];
    }
  else
    {
      vz[evIdA] = -Oparams.wallDiss*vz[evIdA];
    }	 
  lastcol[evIdA] = OprogStatus.time;
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
  bump(evIdA, evIdB, &W);
  /*printf("qui time: %.15f\n", OprogStatus.time);*/
  lastcol[evIdA] = lastcol[evIdB] = OprogStatus.time;

  PredictEvent(evIdA, -1);
  PredictEvent(evIdB, evIdA);
}
void docellcross(int k, double velk, double *rkptr, int cellsk)
{
  if (velk > 0.0)
    {
      cellRange[2 * k] = 1;
      inCell[k][evIdA] = inCell[k][evIdA] + 1;
      if (inCell[k][evIdA] == cellsk) 
	{
#ifdef MD_GRAVITY
	  if (k==2)
	    {
	      printf("Un particella ha superato la massima altezza consentita (%.15f)\n",
		     Lz2 + OprogStatus.extraLz);
	      printf("Aumentare il parametro extraLz e rilanciare la simulazione\n");
	      exit(-1);
	    }
#endif
	  inCell[k][evIdA] = 0;
	  *rkptr = -L2;
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
	}
    }
}
void ProcessCellCrossing(void)
{
#ifdef MD_GRAVITY
  int j; 
#endif
  int k, n;

  UpdateAtom(evIdA);
  /* NOTA: cellList[i] con 0 < i < Oparams.parnum � la cella in cui si trova la particella
   * i-esima mentre cellList[j] con 
   * Oparams.parnum <= j < cellsx*cellsy*cellsz+Oparams.parnum
   * � la prima particella che si trova nella cella j-esima
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
void rebuildLinkedList(void)
{
  int j, n;
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList[j] = -1;
  /* -1 vuol dire che non c'� nessuna particella nella cella j-esima */
  for (n = 0; n < Oparams.parnum; n++)
    {
      atomTime[n] = OprogStatus.time;
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
/* ============================ >>> move<<< =================================*/
void move(void)
{
  int ii;
  double rzmax, zfact;
  /* Zero all components of pressure tensor */
#if 0
  Wxy = 0.0;
  Wyz = 0.0;
  Wzx = 0.0;
  Wxx = Wyy = Wzz = 0.0;
#endif
  /* get next event */
 
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
      calcObserv();
      outputSummary(); 
    }
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
    		  /* se l'energia potenziale � ormai stabile considera il quench finito */
    		  OprogStatus.quenchend = Oparams.curStep;
		}
#else
	      calcKVz();
	      OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
    	      if ( (2.0*K/(3.0*((double)Oparams.parnum)-3.0)) < 
		   OprogStatus.quenchtol)//*Oparams.T)
    		{
#endif
    		  printf("QUENCH DONE! %d\n", Oparams.curStep);
		  OprogStatus.numquench++;
		  /* calcola e salva le misure quando finisce il quench */
		  calcRho();
		  save_rho();
		  calccmz();
		  save_rzcm();
    		  OprogStatus.quenchend = OprogStatus.time;
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
   		  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
		  
		}
 	    }
	  else if ((OprogStatus.time - OprogStatus.quenchend)  < OprogStatus.taptau)
    	    {
	      /* se scalevelsteps = 0 allora scala ogni passo se si sta facendo il 
    		 tapping */
	      OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	      calcKVz();
	      MD_DEBUG4(printf("SCALVEL #%d Vz: %.15f\n", Oparams.curStep,Vz));
	      scalevels(Oparams.T, K, Vz);
	      rebuildLinkedList();
	      rebuildCalendar();
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	    }
	  else
	    {
	      /* start quench (-1  significa che il quench � iniziato) */
	      OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
	      OprogStatus.quenchend = -1;
	    }
	  ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	}
      else if (OprogStatus.scalevel)
	{
	  OprogStatus.nextcheckTime += OprogStatus.rescaleTime;
	  calcKVz();
	  MD_DEBUG2(printf("[TAPTAU < 0] SCALVEL #%d Vz: %.15f\n", Oparams.curStep,Vz));
	  scalevels(Oparams.T, K, Vz);
	  rebuildCalendar();
	  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	  ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	}
#if 0
      else if (2.0*K/(3.0*Oparams.parnum-3.0)>Oparams.T)
	{
	  UpdateSystem();
	  calcKVz();
	  scalevels(Oparams.T, K, Vz);
	}
#endif
#if 0 
	  for (ii= 0; ii < Oparams.parnum; ii++)
	    {
	      if (rz[ii]+Lz*0.5 < 0.0)
		{
		  printf("*******************************************************************\n");
		  printf("STEP: %d\n", Oparams.curStep);
		  printf("evIdA=%d SOTTO MURO PT rz+Lz*0.5=%.30f\n", ii, rz[ii]+Lz*0.5);
		  printf("*******************************************************************\n");
		}
	    }
#endif
    }
  if (OprogStatus.maxquench && OprogStatus.numquench == OprogStatus.maxquench)
    {
      ENDSIM = 1;
    }
 
#if 0 
  if (evIdA > -1 && rz[evIdA] + Lz*0.5 < - 0.5)
    { 
      printf("*******************************************************************\n");
      printf("STEP: %d\n", Oparams.curStep);
      printf("evIdA=%d SOTTO MURO PT rz+Lz*0.5=%.30f\n", evIdA, rz[evIdA]+Lz*0.5);
      printf("*******************************************************************\n");
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