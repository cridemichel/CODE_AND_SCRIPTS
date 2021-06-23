#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPTS 10000
#define MAXFILES 5000
char **fname; 
double L[3], time, *ti, *rotMSD, *alpha2, *alpha2A, *alpha2B, *MSD, *rotMSDA, *MSDA, *rotMSDB, *MSDB, *cc, *rotMSDcls[2], *MSDcls[2], *cc_cls[2];
double **DR, **DR0;
double *r0[3], *w0[3], *rt[3], *wt[3], *rtold[3];
char parname[128], parval[256000], line[256000];
char dummy[2048];
int points=-1, foundDRs=0, foundrot=0, eventDriven=0, skip=1, clusters=0, nongauss=0;
int *isPercPart;
char *pnum;
char inputfile[2048], cluststr[2048];
double storerate = -1;
int bakSaveMode = -1;
void readconf(char *fname, double *ti, double *refTime, int NP, double *r[3], double *w[3], double **DR)
{
  FILE *f;
  int nat=0, i;
  double dt=-1;
  int curstp=-1;
  static int first = 1;

  *ti = -1.0;
  f = fopen(fname, "r");
  while (!feof(f) && nat < 2) 
    {
      fscanf(f, "%d %lf %lf %lf\n", &NP, &L[0], &L[1], &L[2]); 
      for (i = 0; i < NP; i++) 
        {
          fscanf(f, "%[^\n]\n", line); 
          if (!(sscanf(line, "%lf %lf %lf\n", &r[0][i], &r[1][i], &r[2][i])==3))
            {
              sscanf(line, "%lf %lf %lf %[^\n]\n", &r[0][i], &r[1][i], &r[2][i], dummy); 
            }
        }
    }
  /* N.B.nei codici non-event-driven non esiste il parametro time negli store 
   * files ascii, quindi il tempo lo calcolo usando i passi correnti e il passo
   * d'integrazione. */ 
  if (*ti == -1)
    *ti = ((double)curstp)*dt;
  fclose(f);
  first = 0;
}
void print_usage(void)
{
  printf("Usage: calcmsd [ --nongauss/-ng | --skip/-s | --clusters/-c ] <listafile> [number of points]\n");
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
      else if (!strcmp(argv[cc],"--nongauss") || !strcmp(argv[cc],"-ng"))
	{
	  nongauss=1;
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
  FILE *fngA, *fngB, *fng, *f, *f2, *fA, *fB, *f2A, *f2B;
  double *adjDr[3], Dr, Dw, dr, DR2, DR2DR2;
  int c2, i, nfiles, ii, nr1, nr2, a;
  int NP, NPA=-1, NN=-1, fine, JJ, nat, maxl, maxnp, np;
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
  fscanf(f, "%d %lf %lf %lf\n", &NP, &L[0], &L[1], &L[2]); 
  fclose(f);
  if (NPA == -1)
    NPA = NP;
  if (points == -1)
    points=NN;
#if 1
  if ((eventDriven==1 && storerate <= 0.0 && bakSaveMode <= 0)
      || (eventDriven==0 && bakSaveMode <= 0)) 
    NN = 1;
#endif
  if (NN!=1)
    skip = 0;
  
  maxnp = NN + (nfiles-NN)/NN;
  if (points > maxnp)
    points = maxnp;
  ti = malloc(sizeof(double)*points);
  rotMSD = malloc(sizeof(double)*points);
  MSD = malloc(sizeof(double)*points);
  if (nongauss)
    {
      alpha2 = malloc(sizeof(double)*points);
      alpha2A = malloc(sizeof(double)*points);
      alpha2B = malloc(sizeof(double)*points);
    }  
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
      if (nongauss)
	{
	  alpha2[ii] = 0.0;
	  alpha2A[ii] = 0.0;
	  alpha2B[ii] = 0.0;
	}
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
  NPA=NP; 
  if (NPA != NP)
    printf("[MIXTURE] points=%d files=%d NP = %d NPA=%d L=%.15G %.15G %.15G NN=%d maxl=%d\n", points, nfiles, NP, NPA, L[0], L[1], L[2], NN, maxl);
  else
    printf("[MONODISPERSE] points=%d files=%d NP = %d L=%.15G %.15G %.15G NN=%d maxl=%d\n", points, nfiles, NP, L[0], L[1], L[2], NN, maxl);
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
      sprintf(cluststr,"%s.clusters",fname[0]);
      isPercPart = malloc(sizeof(int)*NP);
      for (i=0; i < NP; i++)
	isPercPart[i] = -1;
      if (!(f = fopen(cluststr, "r")))
	{
	  printf("I can not open %s file...\n", cluststr);
	  exit(-1);
	}
      while (!feof(f))
	{
	  fscanf(f, "%[^\n]\n", line);
	  isperc = atoi(strtok(line," "));
	  while ((pnum=strtok(NULL," "))!=0)
	    {
	      if (isperc)
		isPercPart[atoi(pnum)] = 1; 
	      else
		isPercPart[atoi(pnum)] = 0; 
	      //printf("atopi(pnum):%d\n", atoi(pnum));
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
      time = nr1; 
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
	      //if (JJ > 0)
		//printf("nr2=%d nr1=%d nr2-nr1=%d NN=%d %d\n", nr2, nr2, nr2-nr1, (nr2-nr1)%NN);
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
              time=nr2;
	      if (np < points && ti[np] == -1.0)
		{
		  ti[np] = time + refTime;
		  //printf("nr1=%d np=%d time=%.15G\n", nr1, np, ti[np]);
		}
  
	      if (nr2 == nr1)
		continue;
	      for (i = 0; i < NP; i++)
		{
		  if (nongauss)
		    DR2=0;
		  for (a = 0; a < 3; a++)
		    {
		      Dw = wt[a][i] - w0[a][i];
		      Dr = rt[a][i] - r0[a][i];
		      dr = rt[a][i] - rtold[a][i];
		      if (foundDRs)
			{
			  adjDr[a][i] = L[a]*(DR[i][a]-DR0[i][a]); 
			}
		      else
			{
			  if (eventDriven)
			    {
			      if (nr2 > nr1 && fabs(dr) > L[a]*0.5)
                                {
                                  if (dr > 0.0)
                                    adjDr[a][i] -= L[a];
                                  else
                                    adjDr[a][i] += L[a];
                                }
			    }
			  else
			    {
			      if (nr2 > nr1 && fabs(dr) > L[a]*0.5)
                                {
                                  if (dr > 0.0)
                                    adjDr[a][i] -= L[a];
                                  else
                                    adjDr[a][i] += L[a];
                                }
			    }
			}
		      //printf("adjDr[%d][%d]:%f\n", a, i, adjDr[a][i]);

		      tmpdbl = (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
		      MSD[np] += tmpdbl;
		      if (nongauss)
			DR2 += tmpdbl;
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
			      tmpdbl = (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
			      MSDA[np] += tmpdbl;
			      if (foundrot)
				rotMSDA[np] += Dw*Dw;
			    }
			  else
			    {
			      tmpdbl = (Dr+adjDr[a][i])*(Dr+adjDr[a][i]);
			      MSDB[np] += tmpdbl;
			      if (foundrot)
				rotMSDB[np] += Dw*Dw;
			    }
			}
		    }
		  if (nongauss)
	    	    {
		      DR2DR2 = DR2*DR2;
	    	      alpha2[np] += DR2DR2;
		      if (i < NPA)
			alpha2A[np] += DR2DR2;
		      else
			alpha2B[np] += DR2DR2;
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
  if (nongauss)
    fng = fopen("alpha2.dat", "w+");

  if (foundrot)
    f2 = fopen("rotMSDcnf.dat", "w+");
  if (NP != NPA)
    {
      fA = fopen("MSDAcnf.dat", "w+");
      if (foundrot)
	f2A = fopen("rotMSDAcnf.dat", "w+");
      if (nongauss)
	{
  	  fngA = fopen("alpha2A.dat", "w+");
  	  fngB = fopen("alpha2B.dat", "w+");
	}
      fB = fopen("MSDBcnf.dat", "w+");
      if (foundrot)
	f2B = fopen("rotMSDBcnf.dat", "w+");
    }
  for (ii=1; ii < points; ii++)
    {
      //printf("cc[%d]=%f ti=%f\n", ii, cc[ii], ti[ii]);
      if (cc[ii] > 0 && ti[ii] > -1.0)
	{
	  tmpdbl = MSD[ii]/cc[ii]/((double)NP);
	  tmpdbl = tmpdbl*tmpdbl;
	  fprintf(f, "%.15G %.15G %f\n", ti[ii]-ti[0], MSD[ii]/cc[ii]/((double)NP), cc[ii]);
	  if (nongauss)
	    fprintf(fng, "%.15G %.15G %f\n", ti[ii]-ti[0], (3.0/5.0)*alpha2[ii]/cc[ii]/((double)NP)/tmpdbl-1.0, cc[ii]);
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
	  if (nongauss)
	    {
	      tmpdbl = MSDA[ii]/cc[ii]/((double)NPA);
	      tmpdbl = tmpdbl*tmpdbl;
	      fprintf(fngA, "%.15G %.15G %f\n", ti[ii]-ti[0], (3.0/5.0)*alpha2A[ii]/cc[ii]/((double)NPA)/tmpdbl-1.0, cc[ii]);
	      tmpdbl = MSDB[ii]/cc[ii]/((double)NP-NPA);
	      tmpdbl = tmpdbl*tmpdbl;
	      fprintf(fngB, "%.15G %.15G %f\n", ti[ii]-ti[0], (3.0/5.0)*alpha2B[ii]/cc[ii]/((double)NP-NPA)/tmpdbl-1.0, cc[ii]);
	    }
	}
    }
  fclose(f);
  if (nongauss)
    fclose(fng);
  if (foundrot)
    fclose(f2);
  if (NP != NPA)
    {
      if (nongauss)
	{
	  fclose(fngA);
	  fclose(fngB);
	}
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
