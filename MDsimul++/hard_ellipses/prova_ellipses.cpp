#include "./hardell.H"
extern double R[2][2][2], rx[2], ry[2], rz[2], sax[2][2];
//#define DEBUG_HE
#ifdef DEBUG_HE
extern "C" {
  double check_overlap_pw2d(int i, int j, double *shift);
  double check_overlap_polyell_2D(int i, int j, double *shift);
};
#endif
int main (int argc, char **argv)
{
  hardell<double> A, B;
  long long int tt;
  double Lbox, ov=0.0, ovc, ovcpp, shift[2]={0.0,0.0};
  if (argc==1)
    {
      printf("ok argc=0\n");
    }
  srand48(2);
  A.a = 1.0;
  A.b = 2.0;
  B.a = 1.0;
  B.b = 2.0;
  Lbox = 10.0;
  for (tt=0; tt < atof(argv[1]); tt++)
    {
      B.random_orient();
      B.random_box(Lbox);
      if ((ovcpp=overlap(A,B)) < 0.0)
        ov+=1.0;
#ifdef DEBUG_HE
      rx[0] = A.r[0];
      ry[0] = A.r[1];
      R[0][0][0] = A.na[0];
      R[0][0][1] = A.na[1];
      R[0][1][0] = A.nb[0];
      R[0][1][1] = A.nb[1];
      sax[0][0] = A.a;
      sax[0][1] = A.b;

      rx[1] = B.r[0];
      ry[1] = B.r[1];
      R[1][0][0] = B.na[0];
      R[1][0][1] = B.na[1];
      R[1][1][0] = B.nb[0];
      R[1][1][1] = B.nb[1];
      sax[1][0] = B.a;
      sax[1][1] = B.b;

      //ovc=check_overlap_pw2d(0,1,shift);
      ovc=check_overlap_polyell_2D(0,1,shift);
#if 1
      if (ovc*ovcpp < 0)
        {
          printf("ovc=%f ovcpp=%f\n", ovc, ovcpp);
          B.r.show("rB=");
          B.na.show("naB=");
          B.nb.show("nbB=");
          exit(0);
        }
#endif
#endif
    }
}
