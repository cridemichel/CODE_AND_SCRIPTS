#include<mdsimul.h>

#define ITER_ELAPSED 10
/*=== Alcune di queste sono variabili fittizie definite per poter linkare il 
  file mdarray.c === */
extern double L;
static double W, dispHi;
static int nebrNow;
double Fret;
extern char TXT[MSG_LEN], TXTA[10][MSG_LEN];

void mnbrak(double *AX, double *BX, double *CX, double *FA, double *FB, double *FC,
	     double (*func)(double));

extern void BuildNebrListNoLinked(COORD_TYPE rCut, COORD_TYPE sigab[NA][NA]); 

/* ======================== >>> minimumImageFB <<< ======================== */
double minimumImageFB(double c)
{
  if (abs(c) > (L / 2.0))
    {
      if (c > 0)
	c -= L;
      else
	c += L;
    }
  return c;
}

/* ======================= >>> DSIGN <<< ================================ */
double DSIGN(double a, double b)
{
  if (b > 0)
    return fabs(a);
  else
    return -fabs(a);
  //return (fabs(b)/b) * fabs(a);
}

#ifdef TAPPING

#define DABS fabs
#ifdef MPI
MPI_Status status;
#endif
int dietag;

double Ftol, Epoten, EpotIni, Emin, L, invL, fnorm;

/*double Vc, V, W, VAC, WAC, VLJ;*/

#if 0
extern double *xicomx[NA], *xicomy[NA], *xicomz[NA], *pcomx[NA], 
  *pcomy[NA], *pcomz[NA], *xix[NA], *xiy[NA], 
  *xiz[NA], *Gx[NA], *Gy[NA], *Gz[NA], *Hx[NA], *Hy[NA], *Hz[NA];
#endif

void writeEnergy(char filewalter[132], double T, int steps,
		 double VpotIni, double VpotEnd);

void frprmn(int Nm[NA], double Ftol, double *Fret);

void  Forces(int Nm[NA]);

extern void saveBakAscii(char* fn);
extern void LJForce(COORD_TYPE epsab[NA][NA], 
       		    COORD_TYPE sigab[NA][NA], COORD_TYPE rcut);

void conjgrad(void)
{
  Forces(Oparams.parnum);
  
  EpotIni = Epoten;
	
  /*
     minimization configuration
   */
  /*sprintf(TXT, "Using tol= 1E-15\n");
    mdPrintf(ALL, TXT, NULL);*/

  /* Forse e' bene far diventare la tolleranza un parametro in OprogStatus */
  Ftol = OprogStatus.taptol;
  printf("Ftol:%f\n", Ftol);
  frprmn(Oparams.parnum, Ftol, &Fret);
  Forces(Oparams.parnum);
  Emin = Epoten;
  /*sprintf(TXT, "DONE  WORK N.%d\n", numwrk);
    mdPrintf(ALL, TXT, NULL);*/
 
  /* anche il nome del file su cui salvare puo' essere un parametro!!! */ 
  writeEnergy("tapping.dat", Oparams.T, Oparams.curStep, EpotIni, Emin);

  /* E qui dobbiamo rimettere il sistema alla temperatura Oparams.T, anche se ci pensa
     il pistone termico volendo... */
}

/*
    *****************************************************************
    Substitute the first three characters in the file with Qnc
    ******************************************************************  
    */
void writeconf(char filewalter[132])
{
  char fina[132];
  char fina2[512];
  char percorso[512];

  strcpy(fina, filewalter);
  strcpy(percorso, MD_HD_MIS);  
  strcpy(fina2, "Qnc_");
  strcat(fina2, percorso);
  strcat(fina2, fina);
  /*  
      sprintf(TXT, "Scritto il file quenchato:%s\n", fina2);
      mdPrintf(ALL, TXT, NULL);
   */
  saveBakAscii(fina2);
}
#ifdef MPI
extern int my_rank;
#endif
/* ========================= >>> writeEnergy <<< ========================= */
void writeEnergy(char filewalter[132], double T, int steps,
		 double VpotIni, double VpotEnd)
{
  char fina[132], fina2[512], percorso[512];
  FILE* ndat;
#ifdef MPI
  char fina3[32];
#endif  
  strcpy(fina, filewalter);
  strcpy(percorso, OprogStatus.tmpPath); 
  strcpy(fina2, percorso);
  strcat(fina2, "Ene_");
  strcat(fina2, fina);
#ifdef MPI
  sprintf(fina3, "_R%d", my_rank);
  strcat(fina2, fina3);
#endif
 
 printf("fina2:%s\n", fina2);
  ndat = fopen(fina2, "a");
  fprintf(ndat,"%.5G %d %.15G %.15G\n", T, steps, VpotIni, VpotEnd);
  /*
   sprintf(TXT, "Scritto il file Ene:%s\n", fina2);
   mdPrintf(ALL, TXT, NULL);
 */
  fclose(ndat);

}

/* ========================= >>> cleanstring <<< ========================== */
void cleanstring(char *filei)
{
  char fileo[132];
  char recname[132];
  int ex[132], kk, ns;

  //printf("qui\n");
  strcpy(fileo,filei);
  /*
    clean from empty spaces
  */
  
  for(kk = 0; kk < 132; kk++)
    {
      ex[kk] = 1; 
      if (recname[kk] == ' ') 
	ex[kk] = 0; 
    }  
  
  ns=0;
  for(kk = 0; kk < 132; kk++)
    {
      if (ex[kk]) 
	{
	  ns = ns+1;
	  recname[ns] = recname[kk];
	}
    }
  
  for(kk = ns+1; kk < 132; kk++)
    { 
      recname[kk]=' ';
    }

  strcpy(filei,fileo);
  printf("fine\n");
}

/* =========================== >>> min <<< =============================== */
double min(double a, double b)
{
  if (a >= b)
    {
      return b;
    }
  else
    {
      return a;
    }
}

/* =========================== >>> max <<< ================================= */
double max(double a, double b)
{
  if (a >= b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


/* ============================ >>> brent <<< ============================ */
double brent(double AX, double BX, double CX, double (*FFF)(double),
	     double tol, double *xmin)
{
  const int ITMAX = 100;
  const double CGOLD = 0.3819660, ZEPS = 1.0E-15;
  double FU, FV, FW, XM, X, U, D, ETEMP, FX, E, A, B, V, tol1, tol2, R, Q, P;
  int iter;

  A = min(AX, CX);
  B = max(AX, CX);

  D = 0.0; /*
	     verificare che questa inizializzazione sia giusta!!!!!!!!!!!!!!!!!!!!!
	     perché nel codice fortran non anche se viene fatta automaticamente 
	   */
  
  V = BX;
  W = V;
  X = V;
  E = 0.0;
  FX = FFF(X);
  FV = FX;
  FW = FX;

  for(iter = 0; iter < ITMAX; iter++)
    {
      XM = 0.5 * (A+B);
      tol1 = tol*fabs(X) + ZEPS;
      tol2 = 2.0*tol1;
        
      if ( fabs(X-XM) <= (tol2 - 0.5 * (B-A)) ) 
	goto tre;
      
      if(fabs(E) > tol1) 
	{
          R = (X-W) * (FX-FV);
          Q = (X-V) * (FX-FW);
          P = (X-V) * Q - (X-W) * R;
          Q = 2.0 * (Q-R);
          if(Q > 0.0) 
	    {
	      P = -P;
	    }
          Q = fabs(Q);
          ETEMP = E;
          E = D;
          if(fabs(P) >= fabs(0.5 * Q * ETEMP)|| 
	     P <= Q * (A-X) || P >= Q * (B-X)) 
	    goto uno;
	  D = P/Q;
	  U = X+D;
          if(U-A < tol2 || B-U < tol2) 
	    D = DSIGN(tol1,XM-X);
	 
	  goto due;
	}
    uno: 
      if(X >= XM)
	E = A-X;
      else
	E = B-X;
      D = CGOLD * E;
    due:
      if(fabs(D) >= tol1)
	U = X + D;
      else
	U = X + DSIGN(tol1,D);
      
      FU = FFF(U);
      if(FU <= FX)
	{
          if(U >= X)
            A = X;
	  else
            B = X;
	  
          V = W;
          FV = FW;
          W = X;
          FW = FX;
	  X = U;
          FX = FU;
	}
      else
	{
          if (U < X)
            A = U;
          else
            B = U;
          
          if (FU <= FW || W == X)
            {
	      V = W;
	      FV = FW;
	      W = U;
	      FW = FU;
	    }
          else if (FU <= FV || V == X || V == W)
            {
	      V = U;
	      FV = FU;
	    }
	  
	}
    }
  sprintf(TXT, "Brent exceed maximum iterations.\n");
  mdPrintf(ALL, TXT, NULL);
tre:     
  *xmin = X;
  return FX;
}

/* ============================ >>> F1dim <<< ========================== */
double F1dim(double X)
{
  //const int Nmax=1000;
  int j, i, a;
  for (a = 0; a < NA; a++)
    {
      for( j = 0; j < Oparams.parnum[a]; j++)
	{
	  rx[a][j] = pcomx[a][j] + X*xicomx[a][j];
	  ry[a][j] = pcomy[a][j] + X*xicomy[a][j];
	  rz[a][j] = pcomz[a][j] + X*xicomz[a][j];
	  //printf("X:%f xicomx:%f pcomx:%f\n", X, xicomx[j], pcomx[j]);
	}
    }
  /*
    check that all molecules are in the original box
    check of periodic boundary conditions on the read data
  */
  
  /* porta tutte le particelle nel primo box */
  for (a = 0;  a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  rx[a][i] = rx[a][i] - L * rint(invL * rx[a][i]);
	  ry[a][i] = ry[a][i] - L * rint(invL * ry[a][i]);
	  rz[a][i] = rz[a][i] - L * rint(invL * rz[a][i]);
	}   
    }
  //ENERGY
  Forces(Oparams.parnum);
  
  return Epoten;
  //  printf("x: %.6f f1dim: %.6f ncom:%d\n",x,f1dim,ncom)
}

/* ======================= >>> SDchkRebuild <<< ========================= */
void SDchkRebuild(void)
{
  int i, a;
  COORD_TYPE vv, vvMax = 0.0;

  for (a = 0; a < NA; a++)  
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  vv = Sqr(xix[a][i]) + Sqr(xiy[a][i]) + Sqr(xiz[a][i]);
	  if (vv > vvMax) 
	    vvMax = vv;
	}
    }
  dispHi = dispHi + sqrt(vvMax);
  /* If the maximum displacement is too high rebuild Neighbour List
     see Rapaport pag .54 */
  if (dispHi > 0.5 * OprogStatus.rNebrShell)
    nebrNow = 1;

}

/* ========================= >>> linmin <<< ======== ======================= */
void LinMin(int Nm[NA], double *Fret)
{
  int j, i, a;
  //const nmax = 1000;
  const double tol = 1.0E-10;
  double AX, XX, BX, xmin, FX, FB, FA;
  //int Ncom;
 
  //Ncom = N;

  for (a = 0; a < NA; a++)
    {
      for (j = 0; j < Nm[a]; j++) 
	{
	  pcomx[a][j]  = rx[a][j];
	  pcomy[a][j]  = ry[a][j];
	  pcomz[a][j]  = rz[a][j];
	  xicomx[a][j] = xix[a][j];
	  xicomy[a][j] = xiy[a][j];
	  xicomz[a][j] = xiz[a][j];
	}
    }


  Forces(Nm);
  //sprintf(TXT, "Prima di MNBRAK=%.6f\n", Epoten); 
  //mdPrintf(ALL, TXT, NULL);
  //printf("bx et al',box,boy,boz,boxhx,boxhy,boxhz
  AX = 0.0;
  XX = 0.001 / fnorm;
  BX = 0.002 / fnorm;
  
  mnbrak(&AX, &XX, &BX, &FA, &FX, &FB, F1dim);
  *Fret = brent(AX, XX, BX, F1dim, tol, &xmin);
  //sprintf(TXT, "IN LINMIN: MIN: %.6f\n", xmin);
  //mdPrintf(ALL, TXT, NULL);
  /*
    move the system
  */
  for (a = 0; a < NA; a++)
    {
      for(j = 0; j < Nm[a]; j++)
	{
	  xix[a][j] = xmin * xix[a][j];
	  xiy[a][j] = xmin * xiy[a][j];
	  xiz[a][j] = xmin * xiz[a][j];
	  //if ((xmin == 0.0) && (xi[j] != 0.0)
	  //  ) printf *,j,k,xi(j,k)
	  rx[a][j] = pcomx[a][j] + xix[a][j];
	  ry[a][j] = pcomy[a][j] + xiy[a][j];
	  rz[a][j] = pcomz[a][j] + xiz[a][j];
	}
    }
  /* porta tutte le particelle nel primo box */
  for (a = 0; a < NA; a++)
    {
      for (i = 1; i < Nm[a]; i++)
	{
	  rx[a][i]= rx[a][i] - L * rint(invL * rx[a][i]);
	  ry[a][i]= ry[a][i] - L * rint(invL * ry[a][i]);
	  rz[a][i]= rz[a][i] - L * rint(invL * rz[a][i]);
	}   
    }
  //sprintf(TXT, " XMIN= %.6f EN=%.6f\n", xmin, *Fret);
  //mdPrintf(ALL, TXT, NULL);

  SDchkRebuild();

  Forces(Nm);
  //sprintf(TXT, "Dopo di MNBRAK=%.6f\n", Epoten); 
  //mdPrintf(ALL, TXT, NULL);
  if (xmin == 0.0) 
    {
      sprintf(TXT, "ATTENZIONE X MIN= 0\n");
      mdPrintf(ALL, TXT, NULL);
    }

}

/* ======================= >>> frprmn <<< =========================== */
void frprmn(int Nm[NA], double Ftol, double *Fret)
{
  const int ITMAX = 20000;
  const double eps = 1E-14;
  double FP, GG, DGG, GAM;
  int a, j, iter;
  
  //BuildNebrListNoLinked(Nm, Oparams.rcut, Oparams.sigma);
  Forces(Nm);
  FP = Epoten;
  /*
    calculate force modulus
  */
  fnorm = 0.0;

  for (a = 0; a < NA; a++)
    {
      for (j = 0; j < Nm[a]; j++)
	{
	  fnorm = fnorm + Fx[a][j]*Fx[a][j]+Fy[a][j]*Fy[a][j]+Fz[a][j]*Fz[a][j];
	}
    }
  fnorm = sqrt(fnorm);
  
  //write(10,*) 0,fp,fnorm
  
  
  for (a = 0; a < NA; a++)
    {
      for(j = 0; j < Nm[a]; j++)
	{
	  Gx[a][j]  = Fx[a][j];
	  Hx[a][j]  = Gx[a][j];
	  xix[a][j] = Hx[a][j];
	  Gy[a][j]  = Fy[a][j];
	  Hy[a][j]  = Gy[a][j];
	  xiy[a][j] = Hy[a][j];
	  Gz[a][j]  = Fz[a][j];
	  Hz[a][j]  = Gz[a][j];
	  xiz[a][j] = Hz[a][j];
	}
    }
  
  for(iter = 0; iter < ITMAX; iter++)
    {
      LinMin(Nm, Fret);
      if (iter % ITER_ELAPSED == 0)
	{
	  sprintf(TXTA[0],"[#%d] ",iter);
	}
      Forces(Nm);
      
      if (iter % ITER_ELAPSED == 0)
	{
	  sprintf(TXTA[1], "E=%.15G\n", Epoten);
	  mdPrintf(ALL, TXTA[0], TXTA[1], NULL);
	}
      /*
	calculate force modulus
      */
      fnorm=0.0;
      for (a = 0; a < NA; a++)
	{
	  for(j = 0; j < Nm[a]; j++)
	    {
	      fnorm = fnorm + Fx[a][j]*Fx[a][j]+Fy[a][j]*Fy[a][j]+Fz[a][j]*Fz[a][j];
	    }
	}
      fnorm = sqrt(fnorm);
      // SISTEMARE !!!!!!!!!!!!!!!!!!!!!!!!!!
      //write(10,*) iter,*Fret,fnorm;
      
      //print *,iter,*Fret,fnorm
      if (2.0 * fabs(*Fret-FP) <= Ftol*(fabs(*Fret)+fabs(FP) + eps))
	return;
      
      FP = Epoten;
      
      GG = 0.0;
      DGG = 0.0;
      for (a = 0; a < NA; a++) 
	{
	  for(j = 0; j < Oparams.parnum[a]; j++)
	    {
	      GG = GG + Sqr(Gx[a][j]) + Sqr(Gy[a][j]) + Sqr(Gz[a][j]);
	      //DGG=DGG+F(J,K)**2
	      DGG = DGG + (-Fx[a][j] + Gx[a][j])*(-Fx[a][j]) + 
		(-Fy[a][j] + Gy[a][j])*(-Fy[a][j]) +
		(-Fz[a][j] + Gz[a][j])*(-Fz[a][j]);
	    }
	}
      if (GG == 0)
	return;
      
      GAM = DGG / GG;
      for (a = 0;  a < NA; a++)
	{
	  for(j = 0; j < Nm[a]; j++)
	    {
	      Gx[a][j]  = Fx[a][j];
	      Gy[a][j]  = Fy[a][j];
	      Gz[a][j]  = Fz[a][j];
	      Hx[a][j]  = Gx[a][j] + GAM * Hx[a][j];
	      Hy[a][j]  = Gy[a][j]+ GAM * Hy[a][j];
	      Hz[a][j]  = Gz[a][j]+ GAM * Hz[a][j];
	      xix[a][j] = Hx[a][j];
	      xiy[a][j] = Hy[a][j];
	      xiz[a][j] = Hz[a][j];    
	      //         per steepest discent
	    }
	}
    }
  sprintf(TXT, "FRPR maximum iterations exceeded\n");
  mdPrintf(ALL, TXT, NULL);
}


/* ==================== >>> mnbrak <<< =================================== */
void mnbrak(double *AX, double *BX, double *CX, double *FA, double *FB, double *FC,
	     double (*func)(double))
{
  const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0E-15;
  double Q, FU, U, ULIM, DUM, R;

  *FA = func(*AX);
  *FB = func(*BX);
  
  if( *FB > *FA)
    {
      DUM = *AX;
      *AX = *BX;
      *BX = DUM;
      DUM = *FB;
      *FB = *FA;
      *FA = DUM;
    }
  *CX = *BX + GOLD*(*BX-*AX);
  
  *FC = func(*CX);
 inizio:
  if (*FB >= *FC)
    {
      R = (*BX-*AX) * (*FB-*FC);
      Q = (*BX-*CX) * (*FB-*FA);
      U = *BX-((*BX-*CX)*Q - (*BX-*AX)*R)/(2.0*DSIGN(max(DABS(Q-R),TINY),Q-R));
      ULIM = *BX + GLIMIT*(*CX-*BX);
      if((*BX-U)*(U-*CX) > 0)
	{
	  FU = func(U);
          if (FU < *FC)
	    {
	      *AX = *BX;
	      *FA = *FB;
	      *BX = U;
	      *FB = FU;
	      goto inizio;
	    }
          else if(FU > *FB)
	    {
	      *CX = U;
	      *FC = FU;
	      goto inizio;
	    }
	  
	  U = *CX + GOLD*(*CX-*BX);
	  FU = func(U);
	}
      else if ((*CX-U)*(U-ULIM) > 0) 
	{
	  FU = func(U);
      
	  if (FU < *FC)
	    {
	      *BX = *CX;
	      *CX = U;
	      U = *CX + GOLD*(*CX-*BX);
	      *FB = *FC;
	      *FC = FU;
	      FU = func(U);
	    }
	}
      else if((U-ULIM)*(ULIM-*CX) >= 0)
	{
	  U = ULIM;
	  FU = func(U);
	}
      else
	{
	  U = *CX + GOLD*(*CX-*BX);
	  FU = func(U);
	}
      
      *AX = *BX;
      *BX = *CX;
      *CX = U;
      *FA = *FB;
      *FB = *FC;
      *FC = FU;
      goto inizio;
    }
  
  //sprintf(TXT, "fine MNBRAK\n");
  //mdPrintf(ALL, TXT, NULL);
}


/* =========================== >>> forces <<< ======================= */
void  Forces(int Nm[NA])
{
  int i, a;
  double L, invL;
  
  L = cbrt(Vol);
  invL = 1.0/L;

  //printf("Vol:%f\n", Vol);

  /* porta tutte le particelle nel primo box */
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Nm[a]; i++)
	{
	  rx[a][i]=rx[a][i] - L*rint(invL*rx[a][i]);
	  ry[a][i]=ry[a][i] - L*rint(invL*ry[a][i]);
	  rz[a][i]=rz[a][i] - L*rint(invL*rz[a][i]);
	}   
    }

  
  /* azzera tutte le forze */
  //setToZero(Oparams.parnum[0], Fx[0], Fy[0], Fz[0], 
  //	    NULL);  /* Set to zero all the coordinates */
  
  // setToZero(Oparams.parnum[1], Fx[0], Fy[0], Fz[0], 
  //    NULL);  /* Set to zero all the coordinates */

  /* Costruisce le neighbour list */
  
  if (nebrNow)
    {
      nebrNow = 0;
      dispHi = 0.0;
      BuildNebrListNoLinked(Oparams.rcut, Oparams.sigab);
    }

  LJForce(Oparams.epsab, Oparams.sigab, Oparams.rcut);
 
  /* Energia potenziale totale */
  Epoten = Vc;

  //sprintf(TXT,"Epoten IN Forces: %.6f\n",Epoten);
  //mdPrintf(ALL, TXT, NULL);
}
#endif
