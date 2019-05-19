#include "./hardell.H"
#include "../mdlib/boxes.H"
extern double R[2][2][2], rx[2], ry[2], rz[2], sax[2][2];
//#define DEBUG_HE
#define USE_BBOX
#ifdef DEBUG_HE
extern "C" {
  double check_overlap_pw2d(int i, int j, double *shift);
  double check_overlap_polyell_2D(int i, int j, double *shift);
};
#endif
int main (int argc, char **argv)
{
#ifdef USE_BBOX
  rectangle<double> rectA, rectB; 
#endif
  hardell<double> A, B;
  long long int tt;
#if DEBUG_HE
  double ovc, shift[2]={0.0,0.0};
#endif
  double Lbox, ov=0.0, ovcpp;
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

#ifdef USE_BBOX
  rectA.r = A.r;
  rectA.sax[0] = A.a;
  rectA.sax[1] = A.b;
  rectA.R.set_row(0,A.na);
  rectA.R.set_row(1,A.nb);
  rectB.sax[0] = B.a;
  rectB.sax[1] = B.b;
#endif
  for (tt=0; tt < atof(argv[1]); tt++)
    {
      B.random_orient();
      B.random_box(Lbox);
#ifdef USE_BBOX
      rectB.r = B.r;
      rectB.R.set_row(0,B.na);
      rectB.R.set_row(1,B.nb);
      if (overlap(rectA,rectB) < 0.0)
        {
          if ((ovcpp=overlap(A,B)) < 0.0)
            ov+=1.0;
        }
#else
      if ((ovcpp=overlap(A,B)) < 0.0)
        ov+=1.0;
#endif
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
  cout << setprecision(16) << "excluded surface=" << Lbox*Lbox*ov/((double)tt) << "\n";
}
