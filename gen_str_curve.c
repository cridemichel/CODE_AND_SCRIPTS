#include<stdlib.h>
#include<stdio.h>
#include<math.h>
const int NN=30;
const double base=1.3, storerate=0.01;
double beta, t=0.0, t_beg=0.0, t_end=0.0, A=1.0, tau=1.0;
int main(int argc, char** argv)
{
  int i,j;
  if (argc<5)
    {
      printf("gen_str_curve <tau> <A> <beta> <t_beg> <t_end>\n");
      exit(-1);
    }
  t_beg=atof(argv[4]);
  t_end=atof(argv[5]);
  tau=atof(argv[1]);
  A=atof(argv[2]);
  beta=atof(argv[3]);
  t=t_beg;
  j=0;
  while (t < t_end)
    {
      for (i=0; i < NN; i++)
	{
	  printf("%.15G %.15G\n", t, A*exp(-pow(t/tau,beta)));
	  t=t_beg+storerate*pow(base,i)+j*storerate*pow(base,NN);
	}
      j++;
    }
  return 0;
}
