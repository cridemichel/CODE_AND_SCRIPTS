#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 5000
char **fname; 
double L, time, *ti, *rotMSD, *MSD, *rotMSDA, *MSDA, *rotMSDB, *MSDB, *cc, *rotMSDcls[2], *MSDcls[2], *cc_cls[2];
double **DR, **DR0;
double *r0[3], *w0[3], *rt[3], *wt[3], *rtold[3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int points=-1, foundDRs=0, foundrot=0, eventDriven=0, skip=1, clusters=0;
int *isPercPart;
char *pnum;
char inputfile[2048], cluststr[2048];
double storerate = -1;
int bakSaveMode = -1;
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double **DR)
{
  FILE *f;
  int nat=0, i, cpos;
  double dt=-1;
  int curstp=-1;

  *ti = -1.0;
  f = fopen(fname, "r");
  while (!feof(f) && nat < 2) 
    {
      cpos = ftell(f);
      //printf("cpos=%d\n", cpos);
      fscanf(f, "%[^\n]\n",line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	}
      if (nat < 2)
	{
	  fseek(f, cpos, SEEK_SET);
	  fscanf(f, "%[^:]:", parname);
	  //printf("[%s] parname=%s\n", fname, parname);
	  if (!strcmp(parname,"DR"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf %lf %lf ", &DR[i][0], &DR[i][1], &DR[i][2]);
		}
	      foundDRs = 1;
	    }
	  else if (!strcmp(parname,"sumox"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[0][i]); 
		}
	      foundrot = 1;
	    }
	  else if (!strcmp(parname,"sumoy"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[1][i]); 
		}
	    }
	  else if (!strcmp(parname,"sumoz"))
	    {
	      for (i=0; i < NP; i++)
		{
		  fscanf(f, " %lf ", &w[2][i]); 
		}
	    }
	  else if (!strcmp(parname, "time"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *ti = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
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
	  else if (!strcmp(parname, "refTime"))
	    {
	      fscanf(f, "%[^\n]\n", parval);
	      *refTime = atof(parval);
	      //printf("[%s] TIME=%.15G %s\n",fname,*ti, parval);
	    }	
	  else
	    fscanf(f, " %[^\n]\n", parval);
	}
      else
	{
	  for (i = 0; i < NP; i++) 
	    {
	      fscanf(f, "%[^\n]\n", line); 
	      if (!sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3)
		{
		  sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
		}
	    }
	  break; 
	}
    }
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
  if (*ti == -1)
    *ti = ((double)curstp)*dt;
  fclose(f);
}
void print_usage(void)
{
  printf("Usage: calcmsd [ --skip/-s | --clusters/-c ] <listafile> [number of points]\n");
  exit(0);
}
void parse_param(int argc, char** argv)
{
  int cc=1, extraparam=0;
  
  if (argc==1)
    print_usage();
  while (cc < argc)
    {
      //printf("cc=%d extraparam=%d argc=%d\n", cc, extraparam, argc);
      if (!strcmp(argv[cc],"--help")||!strcmp(argv[cc],"-h"))
	{
	  print_usage();
	}
      else if (!strcmp(argv[cc],"--skip") || !strcmp(argv[cc],"-s"))
	{
	  cc++;
	  if (cc == argc)
	    print_usage();
	  skip = atoi(argv[cc]);
	}
      else if (!strcmp(argv[cc],"--clusters") || !strcmp(argv[cc],"-c"))
	{
	  clusters=1;
	}
      else if (cc == argc || extraparam == 2)
	print_usage();
      else if (extraparam == 0)
	{ 
	  extraparam++;
	  //printf("qui1 extraparam:%d\n", extraparam);
	  strcpy(inputfile,argv[cc]);
	}
      else if (extraparam == 1)
	{
	  extraparam++;
	  points = atoi(argv[cc]);
	}
      else
	print_usage();
      cc++;
    }
}

int main(int argc, char **argv)
{
  FILE *f, *f2, *f3, *fA, *fB, *f2A, *f2B;
  double *adjDr[3], Dr, Dw, A1, A2, A3, dr;
  int c1, c2, c3, i, nfiles, nf, ii, nlines, nr1, nr2, a;
  int NP, NPA=-1, NN=-1, fine, JJ, nat, maxl, maxnp, np, NP1, NP2;
  double refTime=0.0, tmpdbl;
  int isperc, kk;
#if 0
  if (argc <= 1)
    {
      printf("Usage: calcmsd <listafile> [number of points]\n");
      //printf("where NN il the lenght of the logarithmic block\n");
      exit(-1);
    }
#endif
  parse_param(argc, argv);
  f2 = fopen(inputfile, "r");
  c2 = 0;
  maxl = 0;
  while (!feof(f2))
    {
      fscanf(f2, "%[^\n]\n", dummy); 
      if (strlen(dummy)+1 > maxl)
	maxl = strlen(dummy)+1;
      c2++;
    }	
  nfiles = c2;
  rewind(f2);
  fname = malloc(sizeof(char*)*nfiles);
  for (ii=0; ii < nfiles; ii++)
    {
      fname[ii] = malloc(sizeof(char)*maxl);
      fscanf(f2, "%[^\n]\n", fname[ii]); 
    }
  fclose(f2);
  f = fopen(fname[0], "r");
  nat = 0;
  while (!feof(f) && nat < 2) 
    {
      fscanf(f, "%[^\n]\n)", line);
      if (!strcmp(line,"@@@"))
	{
	  nat++;
	  if (nat==2)
	    {
	      for (i=0; i < 2*NP; i++)
		fscanf(f, "%[^\n]\n", line);
	      fscanf(f, "%lf\n", &L);
	      break;
	    }
	  continue;
	}
      sscanf(line, "%[^:]:%[^\n]\n", parname, parval); 
      if (!strcmp(parname,"storerate"))
	{
	  eventDriven = 1;
	  storerate = atof(parval);
	}
      if (!strcmp(parname,"bakSavedMode"))
	bakSaveMode = atoi(parval);
      if (!strcmp(parname,"parnum"))
	{
	  if (sscanf(parval, "%d %d ", &NP1, &NP2) < 2)
	    {
	      NP = atoi(parval);
	    }
	  else
	    {
	      NP = NP1+NP2;
	      NPA = NP1;
	    }
	}
      if (!strcmp(parname,"parnumA"))
	NPA = atoi(parval);
      if (nat==0 && !strcmp(parname,"NN"))
	NN = atoi(parval);
    }
  fclose(f);
  if (NPA == -1)
    NPA = NP;
  if (points == -1)
    points=NN;
  if ((eventDriven==1 && storerate <= 0.0 && bakSaveMode <= 0)
      || (eventDriven==0 && bakSaveMode <= 0)) 
    NN = 1;
  if (NN!=1)
    skip = 0;
  
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  if (eventDriven==0)
    L = cbrt(L);
  ti = malloc(sizeof(double)*points);
  rotMSD = malloc(sizeof(double)*points);
  MSD = malloc(sizeof(double)*points);
  MSDA = malloc(sizeof(double)*points);
  MSDB = malloc(sizeof(double)*points);
  rotMSDA = malloc(sizeof(double)*points);
  rotMSDB = malloc(sizeof(double)*points);
  cc = malloc(sizeof(double)*points);
  if (clusters)
    {
      for (kk=0; kk < 2; kk++)
	{
	  cc_cls[kk] = malloc(sizeof(double)*points);
	  MSDcls[kk] = malloc(sizeof(double)*points);
	  rotMSDcls[kk] = malloc(sizeof(double)*points);
	}
    }
  DR = malloc(sizeof(double*)*NP);
  DR0= malloc(sizeof(double*)*NP);
  for (ii = 0; ii < NP; ii++)
    {
      DR[ii]  = malloc(sizeof(double)*3);
      DR0[ii] = malloc(sizeof(double)*3);
    }
  for (ii=0; ii < points; ii++)
    {
      ti[ii] = -1.0;
      rotMSD[ii] = 0.0;
      MSD[ii] = 0.0;
      rotMSDA[ii] = 0.0;
      MSDA[ii] = 0.0;
      rotMSDB[ii] = 0.0;
      MSDB[ii] = 0.0;
      cc[ii]=0.0;
      if (clusters)
	{
	  MSDcls[0][ii]=MSDcls[1][ii]=0.0;
	  cc_cls[0][ii]=cc_cls[1][ii]= 0;
	}
    }
    
  if (NPA != NP)
    printf("[MIXTURE] points=%d files=%d NP = %d NPA=%d L=%.15G NN=%d maxl=%d\n", points, nfiles, NP, NPA, L, NN, maxl);
  else
    printf("[MONODISPERSE] points=%d files=%d NP = %d L=%.15G NN=%d maxl=%d\n", points, nfiles, NP, L, NN, maxl);
  if (eventDriven)
    printf("[ED] Event-Driven simulation\n");
  else
    printf("[MD] Time-Driven simulation\n");
  for (a=0; a < 3; a++)
    {
      r0[a] = malloc(sizeof(double)*NP);
      w0[a] = malloc(sizeof(double)*NP);
      rt[a] = malloc(sizeof(double)*NP);
      rtold[a] = malloc(sizeof(double)*NP);
      wt[a] = malloc(sizeof(double)*NP);
      adjDr[a] = malloc(sizeof(double)*NP); 
    }
  if (clusters)
    {
      sprintf(cluststr,"%s.clusters",fname[nr1]);
      isPercPart = malloc(sizeof(int)*NP);
      for (i=0; i < NP; i++)
	isPercPart[i] = -1;
      f = fopen(cluststr, "r");
      while (!feof(f))
	{
	  fscanf(f, "%[^\n]\n", line);
	  isperc = atoi(strtok(line," "));
	  while (!(pnum=strtok(NULL," ")))
	    {
	      if (isperc)
		isPercPart[atoi(pnum)] = 1; 
	      else
		isPercPart[atoi(pnum)] = 0; 
	    } 
	}
      fclose(f);
    } 
  for (nr1 = 0; nr1 < nfiles; nr1=nr1+NN+skip)
    {	
      for (i=0; i < NP; i++)
	for (a=0; a < 3; a++)
	  adjDr[a][i] = 0.0;

      readconf(fname[nr1], &time, &refTime, NP, r0, w0, DR0);
      fine = 0;
      for (JJ = 0; fine == 0; JJ++)
	{
	  for (nr2 = nr1 + JJ*NN; nr2-nr1-JJ*NN < NN; nr2++)
	    {
	      /* N.B. considera NN punti in maniera logaritmica e poi calcola i punti in maniera lineare 
	       * distanziati di NN punti. */
	      np = (JJ == 0)?nr2-nr1:NN-1+JJ;	      
	      if (nr2 >= nfiles || np >= points)
		{
		  fine = 1;
		  break;
		}
	      if (JJ > 0 && (nr2 - nr1) % NN != 0)
		continue;
		
	      if (nr2==nr1)
		{
		  for (i=0; i < NP; i++)
		    for (a=0; a < 3; a++)
		      rtold[a][i] = r0[a][i];
		}
	      else
		{
		  for (i=0; i < NP; i++)
		    for (a=0; a < 3; a++)
		      rtold[a][i] = rt[a][i];
		}
	      //printf("fname[%d]:%s\n", nr2, fname[nr2]);
	      readconf(fname[nr2], &time, &refTime, NP, rt, wt, DR);
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("nr1=%d time=%.15G\n", np, ti[np]);
		}
  
	      if (nr2 == nr1)
		continue;
	      for (i = 0; i < NP; i++)
		for (a = 0; a < 3; a++)
		  {
		    Dw = wt[a][i] - w0[a][i];
		    Dr = rt[a][i] - r0[a][i];
		    dr = rt[a][i] - rtold[a][i];
		    if (foundDRs)
		      {
			adjDr[a][i] = L*(DR[i][a]-DR0[i][a]); 
		      }
		    else
		      {
			if (eventDriven)
			  {
			    if (nr2 > nr1 && fabs(dr) > L*0.5)
			      if (dr > 0.0)
				adjDr[a][i] -= L;
			      else
				adjDr[a][i] += L;
			  }
		      }
		    //printf("adjDr[%d][%d]:%f\n", a, i, adjDr[a][i]);

		    tmpdbl = (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
		    MSD[np] += tmpdbl;
		    if (clusters)
		      {
			if (isPercPart[i]==1)
			  MSDcls[0][np] += tmpdbl;
			else if (isPercPart[i] == 0)
			  MSDcls[1][np] += tmpdbl;
		      }	  
		    if (foundrot)
		      {
			tmpdbl = Dw*Dw;
			rotMSD[np] += tmpdbl;
			if (clusters)
			  {
			    if (isPercPart[i]==1)
			      rotMSDcls[0][np] += tmpdbl;
			    else if (isPercPart[i] == 0)
			      rotMSDcls[1][np] += tmpdbl;
			  }	  
		      }
	    	    if (NP != NPA)
		      {
			if (i < NPA)
			  {
			    MSDA[np] += (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
			    if (foundrot)
			      rotMSDA[np] += Dw*Dw;
			  }
			else
			  {
		    	    MSDB[np] += (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
	    		    if (foundrot)
			      rotMSDB[np] += Dw*Dw;
			  }
		      }
		    if (clusters)
		      {
			if (isPercPart[i]==1)
			  cc_cls[0][np] += 1.0;  
			else if (isPercPart[i] == 0)
			  cc_cls[1][np] += 1.0;
		      }		      
		  }
	      cc[np] += 1.0;
	      //printf("cc[%d]:%f\n", nr2-nr1, cc[nr2-nr1]);
	    }
	}
    }
  f = fopen("MSDcnf.dat", "w+");
  if (foundrot)
    f2 = fopen("rotMSDcnf.dat", "w+");
  if (NP != NPA)
    {
      fA = fopen("MSDAcnf.dat", "w+");
      if (foundrot)
	f2A = fopen("rotMSDAcnf.dat", "w+");
      fB = fopen("MSDBcnf.dat", "w+");
      if (foundrot)
	f2B = fopen("rotMSDBcnf.dat", "w+");
    }
  for (ii=1; ii < points; ii++)
    {
      //printf("cc[%d]=%f ti=%f\n", ii, cc[ii], ti[ii]);
      if (cc[ii] > 0 && ti[ii] > -1.0)
	{
	  fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], MSD[ii]/cc[ii]/((double)NP), cc[ii]);
	  if (foundrot)
	    fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSD[ii]/cc[ii]/((double)NP), cc[ii]);
	}
      if (NP != NPA)
	{
	  fprintf(fA, "%.15G %.15G %f\n", ti[ii]-ti[0], MSDA[ii]/cc[ii]/((double)NPA), cc[ii]);
	  if (foundrot)
	    fprintf(f2A, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSD[ii]/cc[ii]/((double)NPA), cc[ii]);
	  fprintf(fB, "%.15G %.15G %f\n", ti[ii]-ti[0], MSDB[ii]/cc[ii]/((double)NP-NPA), cc[ii]);
	  if (foundrot)
	    fprintf(f2B, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSDB[ii]/cc[ii]/((double)NP-NPA), cc[ii]);
	}
    }
  fclose(f);
  if (foundrot)
    fclose(f2);
  if (NP != NPA)
    {
      fclose(fA);
      if (foundrot)
	fclose(f2A);
      fclose(fB);
      if (foundrot)
	fclose(f2B);
    }
  if (clusters)
    {
      f = fopen("MSDcnfClsPerc.dat", "w+");
      f2 = fopen("rotMSDcnfClsPerc.dat","w+");
      for (ii=1; ii < points; ii++)
	{
	  //printf("cc[%d]=%f ti=%f\n", ii, cc[ii], ti[ii]);
	  if (cc_cls[0][ii] > 0 && ti[ii] > -1.0)
	    {
	      fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], MSDcls[0][ii]/cc_cls[0][ii], cc_cls[0][ii]);
	      if (foundrot)
		fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSDcls[0][ii]/cc_cls[0][ii], cc_cls[0][ii]);
	    }
	}
      fclose(f);
      fclose(f2);
      f = fopen("MSDcnfClsNotPerc.dat", "w+");
      f2 = fopen("rotMSDcnfClsNotPerc.dat","w+");
      for (ii=1; ii < points; ii++)
	{
	  //printf("cc[%d]=%f ti=%f\n", ii, cc[ii], ti[ii]);
	  if (cc_cls[1][ii] > 0 && ti[ii] > -1.0)
	    {
	      fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], MSDcls[1][ii]/cc_cls[1][ii], cc_cls[1][ii]);
	      if (foundrot)
		fprintf(f2, "%.15G %.15G %f\n", ti[ii]-ti[0], rotMSDcls[1][ii]/cc_cls[1][ii], cc_cls[1][ii]);
	    }
	}
      fclose(f);
      fclose(f2);
    }
  return 0;
}
