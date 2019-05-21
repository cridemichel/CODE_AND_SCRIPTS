
#include "./mcsim_hardell.H"
#include <time.h>
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
  mcsim<double> heMC;
  if (argc==1)
    {
      printf("ok argc=0\n");
    }
  heMC.read_pars("hepars.asc");
  heMC.init();
  heMC.run();

#if 0
  for (tt=0; tt < he; tt++)
    {
      B.random_orient();
      //B.na << 1,0;
      //B.nb << 0,1;
#if 1
      theta=2.0*M_PI*drand48();
      //dr = 2.0*(A.b-A.a); 

      B.r << r*cos(theta), r*sin(theta);
#else
      B.random_box(Lbox);
#endif
      rectB.r = B.r;
      rectB.R.set_row(0,B.na);
      rectB.R.set_row(1,B.nb);

#ifdef USE_BBOX
      if (overlap(rectA,rectB) < 0.0)
        {
          if (overlap(A,B) > 0.0)
            ov+=1.0;
        }
      else 
        ov+=1.0;
#else
      if ((ovcpp=overlap(A,B)) < 0.0)
        ov+=1.0;
#endif
#endif
}
