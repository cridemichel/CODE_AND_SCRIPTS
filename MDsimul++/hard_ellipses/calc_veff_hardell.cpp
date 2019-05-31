#include<iostream>
#include<fstream>
#include "./hardell.H"
#include "../mdlib/boxes.H"
extern double R[2][2][2], rx[2], ry[2], rz[2], sax[2][2];
//#define DEBUG_HE
#define USE_BBOX
//#define FLUSH_BUF
#ifdef DEBUG_HE
extern "C" {
  double check_overlap_pw2d(int i, int j, double *shift);
  double check_overlap_polyell_2D(int i, int j, double *shift);
};
#endif
int main (int argc, char **argv)
{
  ofstream of;
  const char *ofn="veff_vs_tt.dat";
#ifdef USE_BBOX
  rectangle<double> rectA, rectB; 
#endif
  hardell<double> A, B;
  long long int tt, ttmax, savett=10000, ttini;
  int tr, NUMR=30;
#ifdef DEBUG_HE
  double shift[2]={0.0,0.0};
#endif
  double X0, ov=0.0, theta, dr, r, r0, ovini;
  //double Lbox;
  bool justoner=false;
  int ir=0; 
  vector<double> veff;
  double veff1;
  if (argc==1)
    {
      printf("calc_veff_hardell <trials> <X0> <NUMR> <ir> <savett>\n");
      exit(1);
    }
  if (argc >= 2)
    ttmax=atoll(argv[1]);
  else
    ttmax = 1000000;

   if (argc>=3)
    X0 = atof(argv[2]);
  else 
    X0 = 2.0;
  
  if (argc>=4)
    NUMR = atoi(argv[3]);
 
  if (argc >=5)
    {
      justoner=true;
      ir = atoi(argv[4]);
    }
  else
    {
      ir=0;
      justoner =  false;
    }
  if (ir < 0 || ir >= NUMR)
    {
      cout << "ir[=" << ir << "] must be in the range [0," << NUMR << ")\n";
      exit(1);
    }
  if (argc >= 6)
    {
      savett=atoll(argv[5]);
    }
  fstream checkpoint;
  checkpoint.open("calcveff.chk");
  if (!checkpoint.fail())
    {
      cout << "Checkpoint file 'calcveff.chk' found, I am gonna use it\n";
      checkpoint >> ir;
      checkpoint >> ttini;
      checkpoint >> ovini;
      checkpoint >> X0;
      checkpoint >> NUMR;
      checkpoint >> savett;
      checkpoint >> ttmax;
      checkpoint >> justoner;
      checkpoint.close(); 
    }
  else
    {
      ttini=0;
      ovini=0;
    }

  veff.resize(NUMR); 
  srand48(time(0));
  A.a = 0.5;
  A.b = X0*A.a;
  A.na << 1,0;
  A.nb << 0,1;
  B.a = 0.5;
  B.b = X0*A.a;
  //Lbox = 1.0001*(max(B.a,B.b)+max(A.a,A.b))*2.0;

#ifdef USE_BBOX
  rectA.r = A.r;
  rectA.sax[0] = A.a;
  rectA.sax[1] = A.b;
  rectA.R.set_row(0,A.na);
  rectA.R.set_row(1,A.nb);
  rectB.sax[0] = B.a;
  rectB.sax[1] = B.b;
#endif
  dr = 2.0*(A.b-A.a)/((double) NUMR);
  cout << "dr = " << dr << "\n";
  cout << "NUMR=" << NUMR << " doing just ir=" << ir << "\n"; 
  cout << "ttmax=" << ttmax << " savett=" << savett << "\n";
  r = r0 = 2.0*A.a+dr*0.5+dr*ir;
  cout << setprecision(16) << "r=" << r << "\n"; 
  if (justoner)
    {
      of.open(ofn, ios::trunc);
      of.close();  
    }
  
  for (tr=ir; (tr < NUMR)&&(justoner==true && tr < ir+1); tr++)
    {
      //cout << "r=" << r << "\n";
      cout << "doing just ir=" << ir << "\n";
      ov=ovini;
      for (tt=ttini; tt < ttmax; tt++)
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
#else
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
#ifdef DEBUG_HE
          //ovc=check_overlap_pw2d(0,1,shift);
          if (check_overlap_pw2d(0,1,shift) > 0.0)
            {
              ov+=1.0;
            }
          //else
#if 0
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
          if (justoner==true && (tt % savett == 0))
            {
              if (tt > 0)
                {
                  veff1 = ov/((double)tt);
                  of.open(ofn,ios::app);
                  of << tt << " " << setprecision(16) << -log(veff1) << "\n";
                  checkpoint.open("calcveff.chk",ios::trunc|ios::out);
                  if (justoner==true && !checkpoint.fail())
                    {
                      checkpoint << ir << " ";
                      checkpoint << tt << " ";
                      checkpoint << setprecision(20) << ov << " ";
                      checkpoint << X0 << " ";
                      checkpoint << NUMR << " ";
                      checkpoint << savett << " ";
                      checkpoint << ttmax << " ";
                      checkpoint << justoner << "\n";
                      checkpoint.close();
                    }
                  of.close();
#ifdef FLUSH_BUF
                  of.flush();
#endif
                 }
            }
        }
      if (justoner==false)
        {
          veff[tr] = ov/((double)tt);
          r += dr;
        }
    }
  //cout << "veff=\n";
  r = r0;
  if (justoner==false)
    {
      for (auto i=0; i < NUMR; i++)
        {
          veff[i] = -log(veff[i]);
          cout << setprecision(16) << r << " " << veff[i] << "\n";
          r+=dr;
        } 
    }
  //cout << setprecision(16) << "excluded surface=" << Lbox*Lbox*ov/((double)tt) << "\n";
  remove("./calcveff.chk");
}
