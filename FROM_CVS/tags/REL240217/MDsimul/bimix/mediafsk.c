#include <stdio.h>
#include <stdlib.h>
#define maxtime 200
const int npt = maxtime*maxtime;

double sqIm[maxtime*maxtime],sqRe[maxtime*maxtime], nacc[maxtime*maxtime];

int main(void)
{
  int time1[maxtime*maxtime],time2[maxtime*maxtime];
  double sssRe, sssIm, dt;
  char filename[256];
  int np, it1, it2, esiste, kk, ik, i, tt;
  FILE *fsl, *fsd;
  np=0;
    
  fsl = fopen("filesselfCnf22.list", "r");
  
  while (!feof(fsl))
    {
      fscanf(fsl, "%[^\n]\n", filename);
      fsd = fopen(filename, "r");
      fprintf(stderr,"processo file: %s\n", filename);
      /* la prima riga contiene il passo d'integrazione */
      fscanf(fsd, "%lf", &dt);
      while (!feof(fsd))
	{
	  fscanf(fsd, "%d %d %lf %lf\n", &it1, &it2, &sssRe, &sssIm);
	  fprintf(stderr, "nuovi dati: it2=%d\n", it2);
	  /* identificazione di it1,it2 */
	  esiste = 0;
	  for (kk=0; kk < np; kk++)
	    {
	      if (np&&(it1==time1[kk])&&(it2==time2[kk]))
		{
		  /* gia' esiste */
		  ik=kk;
		  sqRe[ik]=sqRe[ik]+sssRe;
		  sqIm[ik]=sqIm[ik]+sssIm;
		  nacc[ik]=nacc[ik]+1;
		  esiste = 1;
		}
	    }
	  if (!esiste)
	    {
	      /* se arriva qui non esiste! */
	      fprintf(stderr," New time: %d %d %d\n",np,it1,it2);
	      if (np > npt) 
		{
		  fprintf(stderr,"Errore dimensionamento");
		  exit(-1);
		}
	      time1[np]=it1;
	      time2[np]=it2;
	      ik=np;
	      sqRe[ik]=sqRe[ik]+sssRe;
      	      sqIm[ik]=sqIm[ik]+sssIm;
	      nacc[ik]=nacc[ik]+1;
	      np=np+1;
	    }
	}
      fclose(fsd);
      /* media dei dati */
    }
  fclose(fsl);
  /* Write results */

  //fsl = fopen("SelfCnfq22.dat", "w");
  
  for (i=0; i < np; i++)
    {
      tt=nacc[i];
      fprintf(stdout,"%.12f %f %d\n", dt*(time2[i]-time1[i]), sqRe[i]/tt, tt);
    }

  //fclose(fsl);
}
