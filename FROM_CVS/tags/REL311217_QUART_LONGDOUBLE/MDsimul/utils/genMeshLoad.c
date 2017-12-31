#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define Sqr(x) (x)*(x)

int esiste_wt(int nct, int wt, int* ct)
{
  int i;
  for (i=0; (nct > 0) && (i < nct); i++)
    if (wt == ct[i])
      return 1;
  return 0;
}
int triplets[1000000][3];
int main (void)
{
  FILE *fs, *fs2;
  int n1,n2,n3, nt;
  int wt, maxn = 150, ntripl, ntripl2, maxik;
  const int maxtripl = 150;
  /* qui k e' in unita' di 2pi/L */
  double Dk = 0.5, maxk = 50.0, k = 0.0, ktry = 0.0;
  int ik;
  int i, nct, ct[1000];
  
  fs = fopen("kmesh.dat", "w");
  fs2 = fopen("ntripl.dat", "w");

  maxik = ceil(maxk / Dk); 
  printf("k=");
  for (ik=2; ik < maxik; ++ik) 
    { 
      k = ik * Dk;
      printf("%f,", k);
      
      /*printf("%d,", ik);*/
      ntripl = 0;

      for (n1=-maxn; n1<=maxn;n1++)
	{
	  for (n2=-maxn;n2<=maxn;n2++)
	    for (n3=-maxn; n3<=maxn; n3++)
	      {
		ktry = sqrt(Sqr(n1) + Sqr(n2) + Sqr(n3));
		if ( (ktry >= k) && (ktry < k+Dk) )
		  {
		    triplets[ntripl][0] = n1;
		    triplets[ntripl][1] = n2;
		    triplets[ntripl][2] = n3;
		    ntripl++;
		  }
	      }
	  }
      /* zero the array with chosen triplets */
      for (i=0; i < maxtripl; i++) 
	{
	  ct[i] = 0;
	}
      printf("triplette: %d\n", ntripl);
      for (nct = 0; (nct < ntripl) && (nct < maxtripl);)
	{
	  /* printf("nct: %d\n", nct);*/
	  wt = (int)((rand()/((double)RAND_MAX+1.0)) * ntripl);
	  if (esiste_wt(nct, wt, ct))  continue;
	  ct[nct] = wt;
	  nct++;
	
	  fprintf(fs, "%d %d %d ", triplets[wt][0], 
		  triplets[wt][1], triplets[wt][2]);
	}
      
      if (ntripl < maxtripl)
	{
	  fprintf(fs2, "%d\n", ntripl);
	  for (i=ntripl; i<maxtripl; i++)
	    fprintf(fs, "0 0 0 ");
	}
      else 
	{
	  fprintf(fs2, "%d\n", maxtripl);
	}
      /* printf("nct: %d, ntripl = %d, maxtripl: %d\n", nct, ntripl,
 	 maxtripl); */
      /*fseek(fs, -1, SEEK_CUR);*/
      fprintf(fs, "\n");
   }
  printf("\n");

  /*fseek(fs2, -1, SEEK_CUR);
  fseek(fs, -1, SEEK_CUR);*/
  fclose(fs);
  fclose(fs2);
  return 0;
}
