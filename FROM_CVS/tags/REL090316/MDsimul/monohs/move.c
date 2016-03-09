#include<mdsimul.h>
#define SIMUL
#ifdef MD_HSVISCO
void calcT(void);
#endif
void rebuildCalendar(void);
extern void comvel (int, COORD_TYPE, COORD_TYPE, int);
extern void InitEventList (void);
extern void writeAllCor (FILE*);
extern void writeAsciiPars (FILE*,struct pascii[]);
#ifdef MD_FPBROWNIAN
extern COORD_TYPE gauss (void);
#endif

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
#ifdef MD_FULL_LANG
extern void bigauss(double sigma1, double sigma2, double c12, double* chsi1p, double* chsi2p);
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
double W=0.0, K, T1xx, T1yy, T1zz,
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
double *atomTime;
int *inCell[3], **tree, *cellList, cellRange[2*NDIM], 
    cellsx, cellsy, cellsz, initUcellx, initUcelly, initUcellz;
int evIdA, evIdB;
extern int poolSize;
#ifndef MD_POLYDISP
extern double *radii;
#else
extern double *radiiINI;
#endif
extern int *scdone;
double presst1=-1.0, presst2=-1.0;
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
double rcutL, radiusL;
double calc_norm(double *vec)
{
  int k1;
  double norm=0.0;
  for (k1 = 0; k1 < 3; k1++)
    norm += Sqr(vec[k1]);
  return sqrt(norm);
}
double get_min_dist (int na, int *jmin, double *shiftmin) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double distMin=1E10,dist,shift[3], ti;
  /*double cells[NDIM];*/
  int kk;
  double r1[3], r2[3];
  int cellRangeT[2 * NDIM], iX, iY, iZ, jX, jY, jZ, k, n;
 /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
  for (k = 0; k < 2 * NDIM; k++) cellRangeT[k] = cellRange[k];

  //calcdist_retcheck = 0;
  for (iZ = cellRangeT[4]; iZ <= cellRangeT[5]; iZ++) 
    {
      jZ = inCell[2][na] + iZ;    
      shift[2] = 0.;
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
		  if (n!=na) 
		    {
		      //dist = calcDistNeg(Oparams.time, 0.0, na, n, shift, r1, r2, &alpha, vecg, 1);
		      ti = OprogStatus.time-atomTime[na];
		      r1[0] = rx[na] + vx[na]*ti;
		      r1[1] = ry[na] + vy[na]*ti;
		      r1[2] = rz[na] + vz[na]*ti;
		      ti = OprogStatus.time-atomTime[n];
		      r2[0] = rx[n] + vx[n]*ti + shift[0];
		      r2[1] = ry[n] + vy[n]*ti + shift[1];
		      r2[2] = rz[n] + vz[n]*ti + shift[2];
		      dist = 0.0;
		      for (kk = 0; kk < 3; kk++)
			dist += Sqr(r1[kk]-r2[kk]);		      
		      dist = sqrt(dist);
		      //printf(">>>> dist (%d-%d) = %.15G\n", na, n, dist);
		      dist -= radii[na] + radii[n];
		      //printf("radii[%d]:%f radii[%d]: %f\n", na, radii[na],n, radii[n]);
		      //if (calcdist_retcheck)
			//continue;
		      if (*jmin == -1 || dist<distMin)
			{
			  distMin = dist;
			  //printf("distmin = %.15G\n", distMin);
			  for (kk = 0; kk < 3; kk++)
			    {
			      //rCmin[kk] = r1[kk];
			      //rDmin[kk] = r2[kk];
			      shiftmin[kk] = shift[kk];
			    }
			  *jmin = n;
			}
		    }
		} 
	    }
	}
    }
  return distMin;
}

double calc_phi(void)
{
  double N = 0;
  //const double pi = acos(0)*2;
  int i ;
  for (i=0; i < Oparams.parnum; i++)
    {
      N += radii[i]*radii[i]*radii[i];
    }
  N *= 4.0*pi/3.0;
  return N / (L*L*Lz);
}

void store_values(int i)
{
  rcutL = Oparams.rcut;
  radiusL = radii[i];
}
void restore_values(int i)
{
  Oparams.rcut = rcutL;
  radii[i] = radiusL;
}
double max_ax(int i)
{
  double ma;
  ma = 0;
  if (radii[i]>ma)
    ma = radii[i];
  return ma;
}
void scale_coords(double sf)
{
  int i; 
  L *= sf;
  Lz*= sf;
  for (i = 0; i < Oparams.parnum; i++)
    {
      rx[i] *= sf;
      ry[i] *= sf;
      rz[i] *= sf;
      radii[i] *= sf;
    }
}
#ifdef MD_POLYDISP
extern double PHI0;
#endif
double scale_radius(int i, int j, double rA[3], double rB[3], double shift[3], double *factor)
{
  int kk;
  double C, Ccur, F, phi0, phi, fact, L2, rAB[3], fact1, fact2, nrAB;

  L2 = 0.5 * L;
  phi = calc_phi();
#ifdef MD_POLYDISP
  phi0 = PHI0;
#else
  phi0 =((double)Oparams.parnum)*Sqr(Oparams.sigma*0.5)*(Oparams.sigma*0.5);
  phi0 *= 4.0*pi/3.0;
  phi0 /= L*L*Lz;
#endif
  C = cbrt(OprogStatus.targetPhi/phi0);
#ifdef MD_POLYDISP
  Ccur = radii[i]/radiiINI[i]; 
#else
  Ccur = radii[i]/(Oparams.sigma*0.5); 
#endif
  F = C / Ccur;
  for (kk=0; kk < 3; kk++)
    {
      rAB[kk] = rA[kk] - rB[kk] - shift[kk];
    }
  nrAB = calc_norm(rAB);
  /* 0.99 serve per evitare che si tocchino */
  if (F < 1)
    {
      //printf("qui...\n");
      fact = F;
    }
  else
    {
      fact1 = (nrAB-radii[j])/radii[i];
      fact2 = F;
      //printf("fact1: %15G fact2: %.15G\n", fact1, fact2);
      if (fact2 < fact1)
	fact = fact2;
      else
	{
	  fact1 *= OprogStatus.scalfact;
	  if (fact1 < 1.0)
	    fact1 = 1.0;
	  fact = fact1;
	}
    }
  //printf("phi=%f fact1=%.8G fact2=%.8G scaling factor: %.8G\n", phi, fact1, fact2, fact);
  radii[i] *= fact;
  *factor = fact;
  if (2.0*radii[i] > Oparams.rcut)
    Oparams.rcut = 2.0*radii[i]*1.01;
  return calc_phi();
}

void rebuild_linked_list()
{
  double L2;
  int j, n;
  L2 = 0.5 * L;
   cellsx = L / Oparams.rcut;
  cellsy = L / Oparams.rcut;
#ifdef MD_GRAVITY
  cellsz = (Lz+OprogStatus.extraLz) / Oparams.rcut;
#else
  cellsz = L / Oparams.rcut;
#endif 
  for (j = 0; j < cellsx*cellsy*cellsz + Oparams.parnum; j++)
    cellList[j] = -1;
  
  /* rebuild event calendar */
  for (n = 0; n < Oparams.parnum; n++)
    {
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
double check_dist_min(int i, char *msg)
{
  int j;
  double distMin=1E60, shift[3];
  
  j = -1;
  distMin = get_min_dist(i, &j, shift);
    
  if (msg)
    printf("[check_dist_min] %s distMin: %.12G\n", msg, distMin);

  return distMin;
} 
double check_alldist_min(char *msg)
{
  int j, i;
  double distMin=1E60, dist;
  double /*rC[3], rD[3],*/ shift[3];
  for (i=0; i < Oparams.parnum; i++)
    {
      j = -1;
      dist = get_min_dist(i, &j, shift);
      if (j != -1 && dist < distMin)
	distMin = dist;
    }
  if (msg)
    printf("[dist all] %s: %.10G\n", msg, distMin);
  return distMin;
  
}
void scale_Phi(void)
{
  int i, j, imin, kk, /*its,*/ done=0;
  static int first = 1;
  static double rad0I, target;
  double distMinT, distMin=1E60, /*rCmin[3], rDmin[3],*/ rAmin[3], rBmin[3]/*, rC[3], rD[3]*/;
  double L2, shift[3], /*shiftmin[3],*/ phi, factor, radai;
  if (OprogStatus.targetPhi <= 0.0)
    return;

  phi=calc_phi();
  if (first)
    {
      first = 0;
#ifdef MD_POLYDISP
      rad0I = radiiINI[0];
#else
      rad0I = Oparams.sigma*0.5;
#endif
      target = cbrt(OprogStatus.targetPhi/calc_phi());
    }
  //printf("[scale_Phi] >>>> Volume fraction: %.15G target: %.15G\n", phi, target);
  //UpdateSystem();   
  L2 = 0.5 * L;
  /* get the minimum distance in the system */
  phi = calc_phi();
  for (kk = 0;  kk < 3; kk++)
    {
      cellRange[2*kk]   = - 1;
      cellRange[2*kk+1] =   1;
    }
  imin = -1;
  for (i = 0; i < Oparams.parnum; i++)
    {
      j = -1;
      if (scdone[i]==1)
	{
	  done++;
	  continue;
	}
      distMin = get_min_dist(i, &j, shift);
      //if (calcdist_retcheck)
	//continue;
      if (j == -1)
	continue;
      //printf("i=%d j=%d distmin=%.10G\n", i, j, distMin);
      rAmin[0] = rx[i];
      rAmin[1] = ry[i];
      rAmin[2] = rz[i];
      rBmin[0] = rx[j];
      rBmin[1] = ry[j];
      rBmin[2] = rz[j];
      //scalfact = OprogStatus.scalfact;
      store_values(i);
      if (distMin < 1E-13)
	continue;
      phi = scale_radius(i, j, rAmin, rBmin, shift, &factor);
      rebuild_linked_list();
      distMinT = check_dist_min(i, NULL);
      if (distMinT  < 0 && fabs(distMinT) > 1E-10)
	{
	  printf("[scale_Phi] distanza minima < 0!\n");
	  exit(-1);
	}
#if 0
      if (calcdist_retcheck)
	{
	  restore_values(i);
	  rebuild_linked_list();
	  continue;
	}
#endif
#if 0
      its = 0;
      while (distMinT < 0)
	{
	  restore_values(i);
	  phi = scale_radius(i, j, rAmin, rBmin, shiftmin, &factor);
	  rebuild_linked_list();
	  distMinT = check_dist_min(i, "Alla fine di calc_Phi()");
	  its++;
	}
#endif
#ifdef MD_POLYDISP
      radai = radiiINI[i];//Oparams.sigma*0.5;
#else
      radai = Oparams.sigma*0.5;
#endif
      /*
      printf("radai: %.15G radii[i]: %.15G radii/radai=%.15G target: %.15G\n",
	     radai, radii[i], radii[i]/radai, target);
      */
      if (fabs(radii[i] / radai - target) < OprogStatus.axestol)
	{
	  done++;
	  scdone[i] = 1;
	  if (done == Oparams.parnum)
	    break;
	}
     
    }

  printf("Scaled axes succesfully phi=%.8G\n", phi);
  rebuild_linked_list();
  rebuildCalendar();
  ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
  if (OprogStatus.storerate > 0.0)
    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
  if (OprogStatus.scalevel)
    ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
  printf("Scaled successfully %d/%d spheres \n", done, Oparams.parnum);
  if (done == Oparams.parnum || fabs(phi - OprogStatus.targetPhi)<OprogStatus.phitol)
    {
      ENDSIM = 1;
      /* riduce gli ellissoidi alle dimensioni iniziali e comprime il volume */
      factor = rad0I/radii[0];
      Oparams.rcut *= factor;
      scale_coords(factor);
    }
}

void outputSummary(void)
{
  FILE *f;
  
  scale_Phi();
  
  /* mettere qualcosa qui */
#if 0
  printf("time=%.15f\n", OprogStatus.time);
  printf("K= %.15f V=%.15f T=%.15f Vz: %f\n", K, V, 
	 (2.0*K/(3.0*Oparams.parnum-3.0)), Vz);
#endif
#ifndef MD_FULL_LANG
  f = fopenMPI(MD_HD_MIS "T.dat", "a");
  fprintf(f, "%.15f %.15f\n", OprogStatus.time, (2.0*K/(3.0*Oparams.parnum-3.0)));
  fclose(f);
  f = fopenMPI(MD_HD_MIS "Vz2.dat", "a");
  fprintf(f, "%.15f %.15f\n", OprogStatus.time, Sqr(Vz));
  fclose(f);
#endif
  printf("Number of collisions: %lld\n", OprogStatus.collCount);
#if 0
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
#endif
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
#ifdef MD_GRAVITY
  double mg; 
#endif
  double rxi, ryi, rzi, rxij, ryij, rzij, rijSq, sigSq, rij;
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
#ifdef MD_POLYDISP
	  sigSq = Sqr(radii[i]+radii[j]);
#else
	  if (OprogStatus.targetPhi > 0.0)
	    sigSq = Sqr(radii[i]+radii[j]);
#endif
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
#ifdef MD_HSVISCO
  double  DTxy, DTyz, DTzx, Txyold, Tyzold, Tzxold;
#endif
#ifdef MD_HSVISCO
  double taus, invm, delpx, delpy, delpz;
#endif
#ifdef MD_FPBROWNIAN_DEBUG
  double check;
  check = sqrt(Sqr(vx[i]-vx[j]) + Sqr(vy[i]-vy[j]) + Sqr(vz[i]-vz[j]));
  if (check<.5) printf ("bump: velocity difference %d,%d %.15G\n", i,j,check);
#endif
#ifdef MD_POLYDISP
  sigSq = Sqr(radii[i]+radii[j]);
#else
  if (OprogStatus.targetPhi > 0.0)
    {
      sigSq = Sqr(radii[i]+radii[j]);
    }
  else
    sigSq = Sqr(Oparams.sigma);
#endif
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
#ifdef MD_HSVISCO
  invm = 1.0/Oparams.m;
  delpx = delvx*Oparams.m;
  delpy = delvy*Oparams.m;
  delpz = delvz*Oparams.m;
  DTxy = delpx*delpy*invm + vx[i]*delpy + delpx*vy[i];
  DTxy += delpx*delpy*invm - vx[j]*delpy - delpx*vy[j]; 
  DTyz = delpy*delpz*invm + vy[i]*delpz + delpy*vz[i];
  DTyz += delpy*delpz*invm - vy[j]*delpz - delpy*vz[j];
  DTzx = delpz*delpx*invm + vz[i]*delpx + delpz*vx[i];
  DTzx += delpz*delpx*invm - vz[j]*delpx - delpz*vx[j];
#endif
  vx[i] = vx[i] + delvx;
  vx[j] = vx[j] - delvx;
  vy[i] = vy[i] + delvy;
  vy[j] = vy[j] - delvy;
  vz[i] = vz[i] + delvz;
  vz[j] = vz[j] - delvz;
#ifdef MD_FULL_LANG
  factor = ( rxij * ( v2x[i] - v2x[j] ) +
  	ryij * ( v2y[i] - v2y[j] ) +
        rzij * ( v2z[i] - v2z[j] ) ) / sigSq;
  if (factor < 0.0)
    {
      delvx = - factor * rxij;
      delvy = - factor * ryij;
      delvz = - factor * rzij;
      v2x[i] = v2x[i] + delvx;
      v2x[j] = v2x[j] - delvx;
      v2y[i] = v2y[i] + delvy;
      v2y[j] = v2y[j] - delvy;
      v2z[i] = v2z[i] + delvz;
      v2z[j] = v2z[j] - delvz;
     }
#endif
#ifdef MD_HSVISCO 
  if (OprogStatus.lastcoll!=-1)
    {
      taus = OprogStatus.time - OprogStatus.lastcoll;
      OprogStatus.DQTxy += OprogStatus.Txy*taus; 
      OprogStatus.DQTyz += OprogStatus.Tyz*taus;
      OprogStatus.DQTzx += OprogStatus.Tzx*taus;
      //taus = Oparams.time - OprogStatus.lastcoll;
      //printf("DQT= %f %f %f\n", OprogStatus.DQTxy, OprogStatus.DQTyz, OprogStatus.DQTzx);
      OprogStatus.DQWxy += rxij*delpy;
      OprogStatus.DQWyz += ryij*delpz;
      OprogStatus.DQWzx += rzij*delpx;
      OprogStatus.DQW += rxij*delpx + ryij*delpy + rzij*delpz;
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
  /* TO CHECK: il viriale ha senso solo se non c'è la gravità */
  if (*W == 0.0)
    presst1 = presst2;
  *W += Oparams.m*(delvx * rxij + delvy * ryij + delvz * rzij);
  presst2 = OprogStatus.time;
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
  /* hhcp è l'altezza delle particelle se fossero close-packed */
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
#ifdef MD_POLYDISP
	      if (rz[n] + Lz2 + radii[n] < hhcp)
    		npart++;
#else
	      if (OprogStatus.targetPhi > 0.0)
		{
		  if (rz[n] + Lz2 + radii[n] < hhcp)
		    npart++;
		}
	      else
		{
		  if (rz[n] + Lz2 + Oparams.sigma*0.5 < hhcp)
		    npart++;
		}
#endif
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

  /* Bisogna considerare le velocità rispetto al centro di massa! */
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
double *treeTime;
void UpdateAtom(int i)
{
#ifdef MD_FPBROWNIAN
  double expt;
#endif
  double ti;
  ti = OprogStatus.time - atomTime[i];
  if (ti<0)
    MD_DEBUG(printf("WARNING: UpdateAtom[%d] with negative ti\n",i));
#ifdef MD_FPBROWNIAN
  expt = exp(-Oparams.xi*ti);
  if (OprogStatus.brownian) {
    #ifdef MD_FPBROWNIAN_DEBUG
    if (ti>0) {
    #endif
    rx[i] += vx[i]*(1-expt)/Oparams.xi;
    ry[i] += vy[i]*(1-expt)/Oparams.xi;
    rz[i] += vz[i]*(1-expt)/Oparams.xi;
    vx[i] *= expt;
    vy[i] *= expt;
    vz[i] *= expt;
    #ifdef MD_FPBROWNIAN_DEBUG
      printf("velval %d %.15G %.15G\n", i, OprogStatus.time, sqrt(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])));
    }
    #endif
  } else {
#endif
  rx[i] += vx[i]*ti;
  ry[i] += vy[i]*ti;
#if defined(MD_GRAVITY)
  rz[i] += vz[i]*ti - g2*Sqr(ti);
  vz[i] += -Oparams.ggrav*ti;
#else
  rz[i] += vz[i]*ti;
#endif
#ifdef MD_FPBROWNIAN
  }
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
if (ti>0)
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
void PredictEvent (int na, int nb) 
{
  /* na = atomo da esaminare 0 < na < Oparams.parnum 
   * nb = -2,-1, 0 ... (Oparams.parnum - 1)
   *      -2 = controlla solo cell crossing e urti con pareti 
   *      -1 = controlla urti con tutti gli atomi nelle celle vicine e in quella attuale 
   *      0 < nb < Oparams.parnum = controlla urto tra na e n < na 
   *      */
  double sigSq, dr[NDIM], dv[NDIM], shift[NDIM], tm[NDIM], b, d, t, tInt, vv;
#ifdef MD_GRAVITY
  double Lzx, cells[NDIM]; 
  double h1, h2, hh1;
#endif
  int cellRangeT[2 * NDIM], signDir[NDIM], evCode,
      iX, iY, iZ, jX, jY, jZ, k, n;
#ifdef MD_FPBROWNIAN
  double expt, Ltmp;
  expt=0;
#endif

  MD_DEBUG(printf("PredictEvent: %d,%d\n", na, nb));
  /* Attraversamento cella inferiore, notare che h1 > 0 nel nostro caso
   * in cui la forza di gravità è diretta lungo z negativo */ 
#ifdef MD_GRAVITY
  Lzx = OprogStatus.extraLz + Lz;

  /* NOTA: Il muro inferiore è posto a Lz / 2 cioè sul fondo del box,
   * in questo modo la base di ogni sfera poggia esattamente 
   * sul fondo della scatola */
  if (inCell[2][na] == 0)
    {
      hh1 =  vz[na] * vz[na] + 2.0 * Oparams.ggrav *
	(rz[na] + Lz2);
#ifdef MD_POLYDISP
      h1 = hh1 -  Oparams.ggrav * 2.0 * radii[na];
#else
      if (OprogStatus.targetPhi > 0.0)
	h1 = hh1 -  Oparams.ggrav * 2.0 * radii[na];
      else
	h1 = hh1 -  Oparams.ggrav * Oparams.sigma;
#endif
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
#ifdef MD_FPBROWNIAN
          Ltmp = ((inCell[2][na] + 1 - signDir[2]) * L / cellsz - rz[na] - L2);
          Ltmp = 1 - Oparams.xi * Ltmp / vz[na];
          tm[2] = ((Ltmp > 0 && Ltmp < 1) ? -log(Ltmp)/Oparams.xi : timbig);
#else
	  tm[2] = ((inCell[2][na] + 1 - signDir[2]) * L / cellsz - rz[na] - L2) / vz[na];
#endif
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
#ifdef MD_FPBROWNIAN
          Ltmp = ((inCell[0][na] + 1 - signDir[0]) * L / cellsx - rx[na] - L2);
          Ltmp = 1 - Oparams.xi * Ltmp / vx[na];
          tm[0] = ((Ltmp > 0 && Ltmp < 1) ? -log(Ltmp)/Oparams.xi : timbig);
#else
	  tm[0] = ((inCell[0][na] + 1 - signDir[0]) * L / cellsx - rx[na] - L2) / vx[na];
#endif
	} 
      else 
	tm[0] = timbig;

      if (vy[na] != 0.) 
	{
	  if (vy[na] > 0.) 
	    signDir[1] = 0;
	  else 
	    signDir[1] = 1;
#ifdef MD_FPBROWNIAN
          Ltmp = ((inCell[1][na] + 1 - signDir[1]) * L / cellsy - ry[na] - L2);
          Ltmp = 1 - Oparams.xi * Ltmp / vy[na];
          tm[1] = ((Ltmp > 0 && Ltmp < 1) ? -log(Ltmp)/Oparams.xi : timbig);
#else
	  tm[1] = ((inCell[1][na] + 1 - signDir[1]) * L / cellsy - ry[na] - L2) / vy[na];
#endif
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
#ifdef MD_GRAVITY
	  printf("rz:%f diff:%f\n", rz[na], rz[na]+Lz2);
	  printf("h1:%f hh1:%f vz:%f cellz:%d\n", h1, hh1, vz[na], inCell[2][na]);
#endif
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
			che è anche l'unica nel presente caso */
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
#ifdef MD_POLYDISP
			  sigSq = Sqr(radii[na]+radii[n]);
#else
		      	  if (OprogStatus.targetPhi > 0.0)
			    {
			      sigSq = Sqr(radii[na]+radii[n]);
			    }
#endif
			  tInt = OprogStatus.time - atomTime[n];
#ifdef MD_FPBROWNIAN
                          if (OprogStatus.brownian) {
                            expt = exp(-Oparams.xi*tInt);
                            dr[0] = rx[na] - (rx[n] + vx[n] * (1-expt)/Oparams.xi) - shift[0];
                            dv[0] = vx[na] - vx[n] * expt;
                            dr[1] = ry[na] - (ry[n] + vy[n] * (1-expt)/Oparams.xi) - shift[1];
                            dv[1] = vy[na] - vy[n] * expt;
                            dr[2] = rz[na] - (rz[n] + vz[n] * (1-expt)/Oparams.xi) - shift[2];
                            dv[2] = vz[na] - vz[n] * expt;
                          } else {
#endif
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
#ifdef MD_FPBROWNIAN
                          }
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
#ifdef MD_FPBROWNIAN
                                  if (OprogStatus.brownian) {
				    expt = 1.0 + Oparams.xi * (sqrt (d) + b) / vv;
                                    t = (expt > 0. ? -log(expt)/Oparams.xi : 0);
                                    #ifdef MD_FPBROWNIAN_DEBUG
                                      printf("%.15G %d %d expt=%.15G t=%.15G tcoll=%.15G\n",
                                      OprogStatus.time, n, na, expt, t, OprogStatus.time+t);
                                      printf("velvol %d %.15G %.15G %.15G\n", na, OprogStatus.time,
                                             sqrt(Sqr(vx[na])+Sqr(vy[na])+Sqr(vz[na])), atomTime[na]);
                                    #endif
                                  } else {
#endif
				  t = - (sqrt (d) + b) / vv;
#ifdef MD_FPBROWNIAN
                                  }
                                  if (expt > 0 && t < 0)
#else
				  if (t < 0)
#endif
				    {
#if 1
				      printf("time:%.15f tInt:%.15f\n", OprogStatus.time,
					     tInt);
				      printf("dist:%.15f\n", sqrt(Sqr(dr[0])+Sqr(dr[1])+
								  Sqr(dr[2]))-1.0 );
				      printf("STEP: %d\n", Oparams.curStep);
				      printf("atomTime: %.10f \n", atomTime[n]);
				      printf("n:%d na:%d\n", n, na);
				      printf("jZ: %d jY:%d jX: %d n:%d\n", jZ, jY, jX, n);
				      exit(-1);
#endif
				      t = 0;
				    }
#ifdef MD_FPBROWNIAN
                                  if (expt > 0) {
#endif
				  ScheduleEvent (na, n, OprogStatus.time + t);
				  MD_DEBUG(printf("schedule event [collision](%d,%d)\n", na, ATOM_LIMIT+evCode));
#ifdef MD_FPBROWNIAN
                                  }
#endif
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
#ifdef MD_HSVISCO
  OprogStatus.lastcoll = OprogStatus.time;
#endif
  /* questa invero è la pressione in eccesso */
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
  comvel(Oparams.parnum, T, Oparams.m, 0); 
}

#ifdef MD_FPBROWNIAN
/* give all velocities random "kicks", with variance sqrt(2kT.xi), where kT.xi is given */
void velsKick (double Txi)
{
  COORD_TYPE rTemp, sumx, sumy, sumz;
  int i;
  int Nm;
  Nm = Oparams.parnum; /* since the code copied from comvel() uses Nm */

  rTemp = sqrt (2*Txi);
  for (i = 0; i < Oparams.parnum; i++)
    {
      vx[i] += rTemp * gauss();
      vy[i] += rTemp * gauss();
      vz[i] += rTemp * gauss();
#ifdef MD_FPBROWNIAN_DEBUG
       printf("velvil %d %.15G %.15G\n", i, OprogStatus.time, sqrt(Sqr(vx[i])+Sqr(vy[i])+Sqr(vz[i])));
#endif
    }

  /* remove net momentum */
  /* TODO: this is copied from init.c:comvel(), and should be a subroutine of its own */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  loop(i, 1, Nm)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 sumx = sumx + vx[i];
	 sumy = sumy + vy[i];
	 sumz = sumz + vz[i];
	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  loop(i, 1, Nm)
    {
      vx[i] = vx[i] - sumx;
      vy[i] = vy[i] - sumy;
      vz[i] = vz[i] - sumz;
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }

  /* TODO: do we need to move the center of mass to the origin? see comvel() routine */
}
#endif
#ifdef MD_MICRO_LANG
extern double gauss(void);
extern double ranf(void);
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
   M = Oparams.m;
   gam = Oparams.xi*M;
   /* il 3 deriva dal fatto che bisogna mediare su metà angolo solido!*/
   m = Oparams.Dt*gam/2.0;//;(3.0 / 2.0); 
   c1 = (M - m)/(M+m);
   kTm = sqrt(Oparams.T / m);
   //printf("c1=%f c2=%f T=%.15G m=%.15G\n", c1, c2, T, m);
   c2 = kTm*2*m / (M+m);
   mredl = m*M/(m+M);
   for (i = 0; i < Oparams.parnum; i++)
    {
#if 1
       vpx = gauss();
       vpy = gauss();
       vpz = gauss();
#else
       vpx = kTm*gauss();
       vpy = kTm*gauss();
       vpz = kTm*gauss();
#endif
       random_direction(&nx, &ny, &nz);	
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
#if 1
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
#if 1
       vx[i] += delpx/M;
       vy[i] += delpy/M;
       vz[i] += delpz/M;
#else
	/* cosi' è come nel caso 1D, ma essendo qui in 3D cio' equivale ad assumere che il versore
	* n sia sempre parallelo alla velocità relativa vij! */
	vx[i] -= vxij*2.0*mredl/M;
	vy[i] -= vyij*2.0*mredl/M;
	vz[i] -= vzij*2.0*mredl/M;
#endif
#endif
       //printf("delp=(%.15G, %.15G, %.15G)\n", delpx, delpy, delpz);
#endif
       //printf("vpx=%.15G,%.15G,%.15G kTm=%.15G\n", vpx, vpy, vpz, kTm);
       //printf("(%.15G, %.15G, %.15G)\n", vx[i], vy[i], vz[i]);
     }  
    
} 
#endif
#ifdef MD_FULL_LANG
void velsFullLang(double T, double xi)
{
 /*
    Translational velocities from maxwell-boltzmann distribution  
    The routine put in vx, vy, vz a velocity choosen from a M.-B. 
    distribution.
    
    The distribution is determined by temperature and (unit) mass.
    This routine is general, and can be used for atoms, linear    
    molecules, and non-linear molecules.                          
    
    ROUTINE REFERENCED:                                          
    
    COORD_TYPE gauss(void)
    Returns a uniform random normal variate from a           
    distribution with zero mean and unit variance.           
    
    VARIABLES 
    COORD_TYPE temp       Temperature 
    m[NA]                 Masses of atoms (NA is the number of atoms)
    int  Nm               Number of molecules  
  */
  double rTemp, sumx, sumy, sumz, RCMx, RCMy, RCMz;
  /*COORD_TYPE Px, Py, Pz;*/
  int i, Nm;
  double kTm, chsidt, echsidt, e2chsidt, chsi, c0, c1, c2, c1mc2, sigmar, sigmav, crv;	 
  double drx, dry, drz, dt, dvx, dvy, dvz;

  Nm = Oparams.parnum;
  //rTemp = sqrt(temp / Oparams.m);  
  kTm = Oparams.T / Oparams.m;
  chsidt = Oparams.xi * Oparams.Dt;
  echsidt  = exp(-Oparams.xi*Oparams.Dt);
  e2chsidt = exp(-2.0*Oparams.xi*Oparams.Dt);
  dt = Oparams.Dt;
  chsi = Oparams.xi;
  c0 = exp(-chsi*dt);
  c1 = (1-c0)/(chsi*dt);
  c2 = (1-c1)/(chsi*dt);
  c1mc2 = c1 - c2;
  /* qui si assume che le masse di tutti gli atomi del disco sono uguali! */
  sigmar = sqrt( dt * (kTm / chsi) *
		(2.0 -  ( 3.0 - 4.0 * echsidt + e2chsidt) / chsidt) );
  sigmav = sqrt(kTm * ( 1.0 - e2chsidt ) );
  crv = dt * kTm * Sqr( 1.0 - echsidt) / chsidt / sigmar / sigmav;
  
  /* variance of the velocities distribution function, we assume k = 1 */ 
  for (i = 0; i < Nm; i++)
    {
      /* Set the velocities of both atoms to the center of mass velocities,
         the exact velocities will be set in the angvel() routine, where we 
         will set:
	 Atom 1: v1  = Vcm + W^(d21 * m2/(m2+m1))
	 Atom 2: v2  = Vcm - W^(d21 * m1/(m1+m2))
	 where Vcm is the center of mass velocity (that is the actual 
	 velocity of both atoms), W is the angular velocity of the molecule,
	 d21 is the vector joining the two atoms (from 2 to 1) and 
	 m1 and m2 are the masses of two atoms 
      */
      bigauss(sigmar, sigmav, crv, &drx, &dvx);
      bigauss(sigmar, sigmav, crv, &dry, &dvy);
      bigauss(sigmar, sigmav, crv, &drz, &dvz);

      /* 14/07/06: CHECK WHAT FOLLOWS!!! */
      vx[i] = v2x[i]*(1.0 - c0)/chsi + drx;
      vy[i] = v2y[i]*(1.0 - c0)/chsi + dry;
      vz[i] = v2z[i]*(1.0 - c0)/chsi + drz;
      vx[i] /= dt;
      vy[i] /= dt;
      vz[i] /= dt;
      v2x[i] = c0*v2x[i] + dvx;
      v2y[i] = c0*v2y[i] + dvy;
      v2z[i] = c0*v2z[i] + dvz;
      //vx[i] = rTemp * gauss();
      //vy[i] = rTemp * gauss();
      //vz[i] = rTemp * gauss();
    }
  /* Remove net momentum, to have a total momentum equals to zero */
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  return; 
  for (i=0; i < Nm; i++)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 sumx = sumx + vx[i];
	 sumy = sumy + vy[i];
	 sumz = sumz + vz[i];
	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i = 0;i <  Nm; i++)
    {
      vx[i] = vx[i] - sumx;
      vy[i] = vy[i] - sumy;
      vz[i] = vz[i] - sumz;
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
    }
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;
  
  for (i=0; i < Nm; i++)
       {
	 /* (sumx, sumy, sumz) is the total momentum */ 
	 sumx = sumx + v2x[i];
	 sumy = sumy + v2y[i];
	 sumz = sumz + v2z[i];
	 //printf("rank[%d] vx[%d]: %.20f\n", my_rank, i, vx[i]);
       }
     
  sumx = sumx / ((COORD_TYPE) Nm ); 
  sumy = sumy / ((COORD_TYPE) Nm );
  sumz = sumz / ((COORD_TYPE) Nm );

  //Px=0.0; Py=0.0; Pz=0.0;
  /* Now (sumx, sumy, sumz) is the total momentum per atom (Ptot/(2*Nm)) */
  for(i = 0;i <  Nm; i++)
    {
      v2x[i] = v2x[i] - sumx;
      v2y[i] = v2y[i] - sumy;
      v2z[i] = v2z[i] - sumz;
      /* In this way the total (net) momentum of the system of 
	 molecules is zero */
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
  int /*ii,*/ innerstep;
  char fileop[1024], fileop2[1024], fileop3[1024];
  FILE *bf;
  const char sepStr[] = "@@@\n";
  /*double rzmax, zfact;*/
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
	  calcObserv();
	  outputSummary(); 
	}
      else if (evIdB == ATOM_LIMIT + 8)
	{
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
	  writeAsciiPars(bf, opro_ascii);
	  fprintf(bf, sepStr);
	  writeAsciiPars(bf, opar_ascii);
	  fprintf(bf, sepStr);
	  writeAllCor(bf);
	  fclose(bf);
#ifdef MD_MAC
#ifdef MPI
          sprintf(fileop3, "/usr/bin/gzip -f %s_R%d", fileop, my_rank);
#else 
          sprintf(fileop3, "/usr/bin/gzip -f %s", fileop);
#endif
#else
#ifdef MPI
          sprintf(fileop3, "/bin/gzip -f %s_R%d", fileop, my_rank);
#else 
          sprintf(fileop3, "/bin/gzip -f %s", fileop);
#endif
#endif
	  system(fileop3);
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
		  taus = OprogStatus.time - OprogStatus.lastcoll; 
		  OprogStatus.DQTxy += taus * OprogStatus.Txy;
		  OprogStatus.DQTyz += taus * OprogStatus.Tyz;
		  OprogStatus.DQTzx += taus * OprogStatus.Tzx;
		  OprogStatus.DQxy = OprogStatus.DQTxy + OprogStatus.DQWxy;
		  OprogStatus.DQyz = OprogStatus.DQTyz + OprogStatus.DQWyz;
		  OprogStatus.DQzx = OprogStatus.DQTzx + OprogStatus.DQWzx;
		  Vol = L*L*Lz;
		  OprogStatus.DQxy /= Vol;
		  OprogStatus.DQyz /= Vol;
		  OprogStatus.DQzx /= Vol;
	      	  //printf("DQ= %f %f %f\n", OprogStatus.DQxy, OprogStatus.DQyz, OprogStatus.DQzx);
		  OprogStatus.lastcoll = OprogStatus.time;
		}
#endif
              #if defined(MD_FPBROWNIAN)
                /* in the FP-Brownian code, instead of re-randomizing velocities,
                   we need to evolve the velocities and then give them a "kick" */
                velsKick(Oparams.T*Oparams.xi);
	      #elif defined(MD_FULL_LANG)
	        velsFullLang(Oparams.T,Oparams.xi);
	      #elif defined(MD_MICRO_LANG)
		velsMicroLang(Oparams.T,Oparams.xi);
	      #else
	        velsBrown(Oparams.T);
              #endif
#ifdef MD_HSVISCO
	      calcT();
#endif
	      rebuildCalendar();
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
	      if (OprogStatus.scalevel)
	        ScheduleEvent(-1, ATOM_LIMIT+9, OprogStatus.nextcheckTime);
	    }
#ifdef MD_HSVISCO
	  else
	    {
	      double Vol;
    	      OprogStatus.DQxy = OprogStatus.DQTxy + OprogStatus.DQWxy;
	      OprogStatus.DQyz = OprogStatus.DQTyz + OprogStatus.DQWyz;
	      OprogStatus.DQzx = OprogStatus.DQTzx + OprogStatus.DQWzx;
	      Vol = L*L*Lz;
	      OprogStatus.DQxy /= Vol;
	      OprogStatus.DQyz /= Vol;
	      OprogStatus.DQzx /= Vol;
	      //printf("STEP #%d DQ= %f %f %f\n", Oparams.curStep, OprogStatus.DQxy, OprogStatus.DQyz, OprogStatus.DQzx);
	    }
#endif
	  if (presst2 > 0.0 && presst1 > 0.0)
	    press = W / (presst2 - presst1) / (L*L*Lz) / 3.0;
	  else 
	    press = 0.0;
	  W = 0;
	  OprogStatus.nextDt += Oparams.Dt;
	  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
	  break;
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
		      /* se l'energia potenziale è ormai stabile considera il quench finito */
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
#ifdef MD_GRAVITY
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
#endif
			   rebuildLinkedList();
			   MD_DEBUG3(distanza(996, 798));
			   rebuildCalendar();
			   ScheduleEvent(-1, ATOM_LIMIT+7, OprogStatus.nextSumTime);
			   if (OprogStatus.storerate > 0.0)
			     ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
			   ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
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
	      	  if (OprogStatus.storerate > 0.0)
	    	    ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
		  ScheduleEvent(-1, ATOM_LIMIT+10,OprogStatus.nextDt);
		}
	      else
		{
		  /* start quench (-1  significa che il quench è iniziato) */
		  OprogStatus.nextcheckTime += OprogStatus.checkquenchTime;
		  OprogStatus.quenchend = -1;
		}
	      if (OprogStatus.checkquenchTime > 0.0)
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
	      if (OprogStatus.storerate > 0.0)
		ScheduleEvent(-1, ATOM_LIMIT+8, OprogStatus.nextStoreTime);
	      if (OprogStatus.scalevel)
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
	{
	  ENDSIM = 1;
	}
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
