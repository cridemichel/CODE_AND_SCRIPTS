#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#define Sqr(x) ((x)*(x))
double *r[3], *R[3][3], L[3], *sax[3], *msax[3], Lx, Ly, Lz;
int *frozen;
const int numprot = 64;
int mcsim=0;
void tRDiagRpw(int i, double M[3][3], double D[3], double Ri[3][3])
{
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = D[0];
  Di[1][1] = D[1];
  Di[2][2] = D[2];
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
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
void xlambda(double lambda, double rA[3], double A[3][3], double rB[3], double B[3][3], double x[3])
{
  double lamA[3][3], onemlamB[3][3], ABL[3][3], invABL[3][3];
  double x1[3], x2[3], x3[3], detinvABL;
  int k1, k2;
  /* calcola xlambda, vedi L. Paramonov and S. N. Yaliraki J. Chem. Phys. 123, 194111 (2005) */
  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  lamA[k1][k2] = lambda*A[k1][k2];
	  onemlamB[k1][k2] = (1.0-lambda)*B[k1][k2];
 	  ABL[k1][k2] = lamA[k1][k2] + onemlamB[k1][k2];
	}
    }
  for (k1=0; k1 < 3; k1++)
    {
      x1[k1]=0;
      x2[k1]=0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  x1[k1] += lamA[k1][k2]*rA[k2];
	  x2[k1] += onemlamB[k1][k2]*rB[k2];
	}
      x3[k1] = x1[k1] + x2[k1];
    }
  detinvABL=-ABL[0][2]*ABL[1][1]*ABL[2][0] + ABL[0][1]*ABL[1][2]*ABL[2][0] + 
    ABL[0][2]*ABL[1][0]*ABL[2][1] - ABL[0][0]*ABL[1][2]*ABL[2][1] - 
    ABL[0][1]*ABL[1][0]*ABL[2][2] + ABL[0][0]*ABL[1][1]*ABL[2][2]; 

  invABL[0][0] = -ABL[1][2]*ABL[2][1] + ABL[1][1]*ABL[2][2];
  invABL[0][1] =  ABL[0][2]*ABL[2][1] - ABL[0][1]*ABL[2][2];
  invABL[0][2] = -ABL[0][2]*ABL[1][1] + ABL[0][1]*ABL[1][2];
  invABL[1][0] =  ABL[1][2]*ABL[2][0] - ABL[1][0]*ABL[2][2]; /* a12 a20 - a10 a22 */
  invABL[1][1] = -ABL[0][2]*ABL[2][0] + ABL[0][0]*ABL[2][2]; /* -a02 a20 + a00 a22 */
  invABL[1][2] =  ABL[0][2]*ABL[1][0] - ABL[0][0]*ABL[1][2]; /* a02 a10 - a00 a12 */
  invABL[2][0] = -ABL[1][1]*ABL[2][0] + ABL[1][0]*ABL[2][1]; /* -a11 a20 + a10 a21 */
  invABL[2][1] =  ABL[0][1]*ABL[2][0] - ABL[0][0]*ABL[2][1]; /* a01 a20 - a00 a21 */
  invABL[2][2] = -ABL[0][1]*ABL[1][0] + ABL[0][0]*ABL[1][1]; /* -a01 a10 + a00 a11 */

  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      invABL[k1][k2] /= detinvABL;

  for (k1 = 0; k1 < 3; k1++)
    {
      x[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  x[k1] += invABL[k1][k2]*x3[k2];
	}
    }
}

double Slam(double lambda, double rA[3], double A[3][3], double rB[3], double B[3][3])
{
  int k1, k2;
  double xlam[3], fA[3], fB[3], SA, SB;

  xlambda(lambda, rA, A, rB, B, xlam);

  for (k1=0; k1 < 3; k1++)
    {
      fA[k1] = 0;
      fB[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  fA[k1] += A[k1][k2]*(xlam[k2]-rA[k2]);
	  fB[k1] += B[k1][k2]*(xlam[k2]-rB[k2]);
	}
    }

  SA = SB = 0.0;
  for (k1=0; k1 < 3; k1++)
    {
      SA += lambda*(xlam[k1]-rA[k1])*fA[k1];
      SB += (1.0-lambda)*(xlam[k1]-rB[k1])*fB[k1];
    }
  /* ho messo un - così la funzione ha un minimo invece
     che un massimo e questo minimo viene trovato dalla funzione brentPW */
  return -(SA+SB);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
int brentPWTooManyIter=0;
double brentPW(double ax, double bx, double cx, double tol, double *xmin, double rA[3], double A[3][3], double rB[3], double B[3][3])
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
 * (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)),
 * this routine isolates the minimum to a fractional precision of about tol using Brent's
 * method. The abscissa of the minimum is returned as xmin, and the minimum function value 
 * is returned as brent, the returned function value. */
{ 
  int iter, ITMAXBR=100;
  const double CGOLD=0.3819660;
  const double ZEPSBR=1E-20;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0, fuold;
  brentPWTooManyIter=0;
  /* This will be the distance moved on the step before last.*/
  a=(ax < cx ? ax : cx); /*a and b must be in ascending order, 
			   but input abscissas need not be.*/
  b=(ax > cx ? ax : cx);
  x=w=v=bx; /*Initializations...*/
  fw=fv=fx=Slam(x, rA, A, rB, B); 
  if (fw < -1.0)
    {
      /* non-overlap! */
      *xmin=x;
      return -100.0;
    }
  fuold = fv;
  for (iter=1;iter<=ITMAXBR;iter++)
    { 
      /*Main program loop.*/
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPSBR); 
      if (fabs(x-xm) <= (tol2-0.5*(b-a)))
	{ /*Test for done here.*/
	  *xmin=x;
	  return fx;
	} 
      if (fabs(e) > tol1) 
	{ /*Construct a trial parabolic fit.*/
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if (q > 0.0)
	    p = -p; 
	  q=fabs(q);
	  etemp=e; 
	  e=d; 
	  if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=CGOLD*(e=(x >= xm ? a-x : b-x)); 
	    /*The above conditions determine the acceptability of the parabolic fit.
	     * Here we take the golden section step into the larger of the two segments.*/
	  else
	    {
	      d=p/q; /* Take the parabolic step.*/
	      u=x+d; 
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x); 
	    }
	}
      else
	{
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	} 
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=Slam(u, rA, A, rB, B); /*This is the one function evaluation per iteration.*/
      if (fu < -1.0)
	{
	  /* non overlap! */
	  *xmin=x;
	  return -100.0;
	}
#if 0
      if (2.0*fabs(fuold-fu) <= tol*(fabs(fuold)+fabs(fu)+ZEPSBR)) 
	{ 
	  *xmin=u;
	  return fu;
	}
#endif
      fuold = fu;//
      if (fu <= fx)
	{ /*Now decide what to do with our function evaluation.*/
	  if (u >= x) 
	    a=x;
	  else
	    b=x;
	  SHFT(v,w,x,u); /* Housekeeping follows:*/
	  SHFT(fv,fw,fx,fu); 
	} 
      else
	{ 
	  if (u < x) 
	    a=u; 
	  else 
	    b=u; 
	  if (fu <= fw || w == x)
	    {
	      v=w; w=u; fv=fw; fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w)
	    { 
	      v=u; fv=fu;
	    }
	} /* Done with housekeeping. Back for another iteration.*/
    }
  printf("Too many iterations in brent!\n");
  brentPWTooManyIter=1;
  //nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}

double check_overlap_pw(int i, int j, double shift[3])
{
  const double tolPW=1.0E-12;
  double res, A[3][3], B[3][3], xmin, RMi[3][3], RMj[3][3]; 
  int k1, k2, a, b;
  double  DA[3], DB[3], rA[3], rB[3];


  rA[0] = r[0][i];
  rA[1] = r[1][i];
  rA[2] = r[2][i];

  rB[0] = r[0][j]+shift[0];
  rB[1] = r[1][j]+shift[1];
  rB[2] = r[2][j]+shift[2];

  for (k1=0; k1 < 3; k1++)
    {
      DA[k1]= 1.0/Sqr(sax[k1][i]);
      DB[k1]= 1.0/Sqr(sax[k1][j]);
    }
  for (a=0; a < 3; a++)
    for (b=0; b < 3; b++)
      {
	RMi[a][b] = R[a][b][i];
	RMj[a][b] = R[a][b][j];
      }
  tRDiagRpw(i, A, DA, RMi);
  tRDiagRpw(j, B, DB, RMj);

  res =  - brentPW(0, 0.5, 1.0, tolPW, &xmin, rA, A, rB, B);
  if (brentPWTooManyIter)
    {
      printf("res=%f xmin=%f\n", res, xmin);
      exit(-1);
    }
  //printf("res=%f\n", res);
  return res - 1.0;
}
char line[4096], parname[1024], parval[4096];
int NP, NPA;
void readconf(char *fname, double *r[3], double *R[3][3], int *NP, int *NPA)
{
  FILE *f;
  int nat=0, i, cpos, dummyint;
  double dt=-1;
  int curstp=-1, NP1, NP2;
  //  *ti = -1.0;
  //printf("fn=%s\n", fname);
  f = fopen(fname, "r");
  while (!feof(f)) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n] ",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  //printf("qui nat=%d\n", nat);
	  continue;
	}
	//printf("line=%s\n", line);
      if (nat < 2)
	{
	  fseek(f, cpos, SEEK_SET);
	  fscanf(f, "%[^:]:", parname);
	  //printf("[%s] parname=%s\n", fname, parname);
	  if (!strcmp(parname,"parnum"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
		{
		  *NP = atoi(parval);
		}
	      else
		{
		  *NP = NP1+NP2;
		  *NPA = NP1;
		}
	    }
	  else if (!strcmp(parname,"parnumA"))
	    *NPA = atoi(parval);
#if 0
	  else if (!strcmp(parname,"DR"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf %lf %lf ", &DR[0][i], &DR[1][i], &DR[2][i]);
		}
	      foundDRs = 1;
	    }
	  else if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	      foundrot = 1;
	    }
	  else if (!strcmp(parname,"sumoy"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf ", &w[1][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoz"))
	    {
	      for (i=0; i < *NP; i++)
		{
		  fscanf(f, " %lf ", &w[2][i]); 
		}
	    }
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      //*ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
#endif
	  else if (!strcmp(parname, "curStep"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      curstp = atoi(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else if (!strcmp(parname, "steplength"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      dt = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }
#if 0
	  else if (!strcmp(parname, "refTime"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *refTime = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
#endif
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else if (nat==2)
	{
	  if (*NPA==-1)
	    *NPA = *NP;
	  fseek(f, cpos, SEEK_SET);
	  for (i = 0; i < *NP; i++) 
	    {
	      fscanf(f, "%[^\n]\n", line); 
#if 0	      
	      if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
		{
		  sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
		}
#endif
	      sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    	      &(r[0][i]), &(r[1][i]), &(r[2][i]), 
	    	      &(R[0][0][i]), &(R[0][1][i]), &(R[0][2][i]), &(R[1][0][i]), &(R[1][1][i]), 
	    	      &(R[1][2][i]), &(R[2][0][i]), &(R[2][1][i]), &(R[2][2][i]), &sax[0][i], &sax[1][i], &sax[2][i],
		      &msax[0][i], &msax[1][i], &msax[2][i]);

	      //printf("r[%d]=%f %f %f sax=%f %f %f\n", i, r[0][i], r[1][i], r[2][i], sax[0][i], sax[1][i], sax[2][i]);
	    }
	  if (mcsim==1)
	    {
	      fscanf(f, "%lf %lf %lf\n", &L[0], &L[1], &L[2]);
	      break;
	    }
	  else
	    {
	      /* if MD read vels */	
	      for (i=0; i < *NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      // fscanf(f, "%lf\n", &L);
	      fscanf(f, "%lf %lf %lf\n", &L[0], &L[1], &L[2]);
	      break;
	    }

	  break; 
	}
    }
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
#if 0
  if (*ti == -1)
    *ti = ((double)curstp)*dt;
#endif 
  fclose(f);
}
char fname[1024], inputfile[1024];
void main(int argc, char **argv)
{
  FILE *f;
  double delsf, nf=0.0, shift[3], sf, avsax[3], saxl[3];
  int i, j, a, b, done=0, overlap=0;
  strcpy(inputfile, argv[1]);

  frozen = malloc(sizeof(int)*numprot);
  for (a=0; a < 3; a++)
    {
      sax[a] = malloc(sizeof(double)*numprot);
      msax[a] = malloc(sizeof(double)*numprot);
      r[a] = malloc(sizeof(double)*numprot);
      for (b=0; b < 3; b++)
	{
	  R[a][b] = malloc(sizeof(double)*numprot);
	}
    }
  f = fopen(inputfile,"r");
  for  (a=0; a <3; a++)
    {
      avsax[a] = 0.0;
    }
  while (!feof(f)) 
    {
      fscanf(f, "%[^\n]\n", fname);
      nf++; 
      printf("fname=%s\n", fname);
      readconf(fname, r, R, &NP, &NPA);
      if (nf==1)
	printf("BOX: %f %f %f\n", L[0], L[1], L[2]);
      delsf=0.05;
      for (i=0; i < numprot; i++)
	frozen[i] = 0;
      for (sf=delsf; sf < 100. && done==0; sf+=delsf)
	{
	  done = 1;
	  for (i=0; i < numprot; i++)
	    {
	      if (!frozen[i])
		done=0;
	    }

	  for (i=0; i < numprot; i++)
	    {
	      if (!frozen[i])
		{
		  sax[0][i] = msax[0][i]*sf;
		  sax[1][i] = msax[1][i]*sf;
		  sax[2][i] = msax[2][i]*sf;
		}
	      for (i=0; i < numprot; i++)
		{
#if 0
		  done=0;
		  saxl[0] = msax[0][i]*sf;
		  saxl[1] = msax[1][i]*sf;
		  saxl[2] = msax[2][i]*sf;
		  printf("i=%d saxl= %f %f %f\n", i, saxl[0], saxl[1], saxl[2]);
		  if (saxl[0] > msax[0][i] && saxl[1] > msax[1][i] && saxl[2] > msax[2][i])
		    {
		      done=1;
		      for (a=0; a < 3; a++)
			avsax[a] += msax[a][i];
		      break;
		    }
#endif
		  overlap = 0;
		  for (j=0; j < numprot; j++)
		    {
		      shift[0] = L[0]*rint((r[0][i]-r[0][j])/L[0]);
		      shift[1] = L[1]*rint((r[1][i]-r[1][j])/L[1]);
		      shift[2] = L[2]*rint((r[2][i]-r[2][j])/L[2]);
		      if (check_overlap_pw(i, j, shift) < 0.0)
			{
			  overlap = 1;
			  for (a=0; a < 3; a++)
			    { 
			      if (frozen[i]==0)
				{
				  sax[a][i]=(sf-delsf)*msax[a][i];
				  frozen[i] = 1;
				}
			      if (frozen[j]==0)
				{
				  sax[a][i]=(sf-delsf)*msax[a][i];
				  frozen[i] = 1;
				}
			    }
			  break;
			}
		    }
#if 0
		  if (overlap==1)
		    {
		      /* overlap not found */
		      saxl[0] = msax[0][i]*(sf-delsf);
		      saxl[1] = msax[1][i]*(sf-delsf);
		      saxl[2] = msax[2][i]*(sf-delsf);
		      for (a=0; a < 3; a++)
			avsax[a] += saxl[a];
		      done=1;
		    }
#endif
		}
	    }
	}
      /* calc average semi-axes */
      for (i=0; i < numprot; i++)
	{
	  for (a=0; a < 3; a++)
	    avsax[a] += sax[a][i];
	}
    }
  printf("Average Semi-axes: %f %f %f\n", avsax[0]/((double)nf)/((double)numprot),
	 avsax[1]/((double)nf)/((double)numprot), avsax[2]/((double)nf)/((double)numprot));
  fclose(f);
}
