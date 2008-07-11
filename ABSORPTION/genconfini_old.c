#include<stdlib.h>
#include<stdio.h>
#include<math.h>
int Nrecept;
int Nprot, brownian;
double L, sigRecept, massRecept, sigProt, massProt, T;
double ranf(void)
{
 return rand() / ( (double) RAND_MAX );
}

double gauss(void)
{
  double  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0; i < 12; i++)
    {
      sum = sum + ranf();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;

}
int main(int argc, char** argv)
{
 FILE *f;
 double K;
 int i;
 T = 1.0;
 Nrecept = 1;
 Nprot=1000;
 sigRecept = 1.0;
 sigProt = 1.0;
 massProt = 1.0;
 massRecept = 1.0E200;
 brownian = 1;
 L=15.0;
 f=fopen("absorb.cnf","w");
 fprintf(f,"parnum:%d\n", Nrecept+Nprot); 
 fprintf(f,"totStep:%d\n", 20000);
 fprintf(f,"time:%.15G\n",0.0); 
 fprintf(f,"curStep:%d\n", 0); 
 fprintf(f,"P:%.15G\n", 1.0);
 fprintf(f,"T:%.15G\n", 1.0);
 fprintf(f,"rcut:%.15G\n", 1.05);
 fprintf(f,"equilibrat: %d\n", 0); 
 fprintf(f,"Dt:%.15G\n", 0.05); 
 fprintf(f,"ninters:%d\n", 1);
 fprintf(f,"ntypes:%d\n", 2); 
 fprintf(f,"@@@\n");
 fprintf(f,"%d %d\n", Nrecept, Nprot);
 /* Receptor type */
 fprintf(f, "%.15G %.15G %.15G\n", sigRecept/3.0, sigRecept/3.0, sigRecept/3.0);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massRecept, 1.0, 1.0, 1.0, 0, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigRecept);
 /* protein type */ 
 fprintf(f, "%.15G %.15G %.15G\n", sigProt/3.0, sigProt/3.0, sigProt/3.0);
 fprintf(f, "%.15G %.15G %.15G\n", 1.0, 1.0, 1.0);
 fprintf(f, "%.15G %.15G %.15G %.15G %d %d\n", massProt, 1.0, 1.0, 1.0, brownian, 1);
 fprintf(f, "%d %d\n", 1, 0);
 fprintf(f, "%.15G %.15G %.15G %.15G\n", 0.0, 0.0, 0.0, sigProt);
 /* interactions */ 
 fprintf(f, "%d %d %d %d %.15G %.15G %.15G %d\n", 0, 0, 1, 0, 1.0, 0.0, 0.0, 10);
 fprintf(f, "@@@\n");
 /* positions */
 fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 0\n", 0.0, 0.0, -L*0.5);
 for (i=0; i < Nprot; i++)
   {
     fprintf(f, "%.15G %.15G %.15G 1 0 0 0 1 0 0 0 1 1\n", (L-1.0)*(ranf()-0.5)+0.5, (L-1.0)*(ranf()-0.5)+0.5,
	     (L-1.0)*(ranf()-0.5)+0.5); 
   }
 /* velocities */
 fprintf(f, "0.0 0.0 0.0 0.0 0.0 0.0\n");
 K = sqrt(T/massProt);
 for (i=0; i < Nprot; i++)
   {
     fprintf(f, "%.15G %.15G %.15G 0.0 0.0 0.0\n", K*gauss(), K*gauss(), K*gauss());
   } 
 fprintf(f, "%.15G\n", L);
 fclose(f);
}
