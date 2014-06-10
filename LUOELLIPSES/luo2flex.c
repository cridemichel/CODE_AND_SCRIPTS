#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void main(int argc, char **argv)
{
  FILE *outf, *inpf;
  int parnum;
  double k, k2, Lx, Ly, Lz;

  inpf = fopen(argv[1], "r");
  outf = fopen(argv[2], "w+");

  k2=k/2.0;
  fprintf(outf,"totStep: 50000\n");
  fprintf(outf,"equilibrat: 0\n");
  fprintf(outf,"curStep: 4400\n");
  fprintf(outf,"M: 5\n");
  fprintf(outf,"rcut: 16.766\n");
  fprintf(outf,"P: 1\n");
  fprintf(outf,"ninters: 0\n");
  fprintf(outf,"T: 1\n");
  fprintf(outf,"tol: 0\n");
  fprintf(outf,"time: 0\n");
  fprintf(outf,"Dt: 0.05\n");
  fprintf(outf,"parnum: %d\n", parnum);
  fprintf(outf,"ntypes: 1\n");
  fprintf(outf,"@@@\n");
  fprintf(outf,"%d\n",parnum); 
  fprintf(outf,"%f 0.5 0.5\n", k2); 
  fprintf(outf,"1 1 1\n"); 
  fprintf(outf,"%f %f %f %f 0 0\n", I, mass); 
  fprintf(outf,"0 0\n");
  fprintf("@@@\n");
  for (i=0; i < parnum; i++)
    {
      fprintf(outf, "%f %f %f %f %f %f %f %f %f %f %f %f 0\n", rx[i], ry[i], rz[i], R[0][0][i],R[0][0][i],R[0][0][i] 

	     );
    }
  for (i=0; i < parnum; i++)
    {
      fprintf(outf, "");
    }
  fprintf(outf, "%f %f %f\n", Lx, Ly, Lz);

  fclose(outf);
  fclose(inpf);
}
