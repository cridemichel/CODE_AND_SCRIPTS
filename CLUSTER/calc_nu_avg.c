#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define maxlen 1000
double avgcls[maxlen], count[maxlen];
void main(int argc, char** argv)
{
  FILE *f, *f2;
  char fn[256];
  int i, len, eml=0;
  double val, val1;
	
  for (i=0; i < maxlen; i++)
    {	
      avgcls[i]=count[i]=0.0;
    }
  f = fopen(argv[1],"r");
  while (!feof(f)) 
   {
      fscanf(f,"%s\n", fn);
      f2=fopen (fn, "r");
      while (!feof(f2))
	{
          fscanf(f2,"%lf %lf\n", &val1, &val);
	  len = (int) val1;
	  if (len < maxlen-1)
            {
              if (len > eml)
		eml=len;
              avgcls[len]+=val;
              count[len]++; 
            }
         }
      fclose(f2);

   }


  fclose(f);
  for (i=0; i <= eml; i++)
    if (count[i]!=0)
       printf("%i %.15G\n", i, avgcls[i]/count[i]);
}
