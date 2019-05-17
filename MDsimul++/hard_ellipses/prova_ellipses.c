#include "./hardell.H"
extern double R[2][2][2], rx[2], ry[2], rz[2], sax[2][2];
#define DEBUG_HE
#ifdef DEBUG_HE
extern "C" {
  double check_overlap_pw2d(int i, int j, double *shift);
};
#endif
int main (int argc, char **argv)
{
  hardell<double> A, B;
  long long int tt;
  double ov=0.0, ovc, ovcpp, shift[2]={0.0,0.0};
  if (argc==1)
    {
      printf("ok argc=0\n");
    }
  for (tt=0; tt < atof(argv[1]); tt++)
    {
      if ((ovcpp=overlap(A,B)) < 0.0)
        ov+=1.0;
#ifdef DEBUG_HE
      ovc=check_overlap_pw2d(0,1,shift);
      if (ovc*ovcpp < 0)
        {
          B.r.show("rB=");
          B.na.show("naB=");
          B.nb.show("nbB=");
          exit(0);
        }
#endif
    }
}
