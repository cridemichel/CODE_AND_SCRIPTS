#include "./mcsim_hardell.H"
int main (int argc, char **argv)
{
  mcsim<double> heMC;
  //heMC.read_pars("hepars.asc");
  auto simtype=atoi(argv[1]);
  if (argc > 1 || simtype!=0)
    {
      if (simtype==1)
        heMC.restart();
      else 
        heMC.readconf(heMC.pars.iniconf);
    }
  else
    heMC.createconf(0);

  heMC.run();
#if 0
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
