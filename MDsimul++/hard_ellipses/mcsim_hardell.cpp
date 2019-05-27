#include "./mcsim_hardell.H"
int main (int argc, char **argv)
{
  char *rf;
  mcsim<double> heMC;
  //heMC.read_pars("hepars.asc");
  int simtype;
  if (argc >= 1)
    simtype = atoi(argv[1]);
  else simtype=0;
  
  if (simtype==1 && argc==2)
    {
      cout << "You have to supply the filename of the restart file\n";
      exit(0);
    }
  if (argc >= 2)
    rf = argv[2];
  else
    rf = nullptr;

  switch (simtype)
    {
    case 1:
        {
          heMC.restart(rf);
          break;
        }
    case 2:
        {
          heMC.readconf(heMC.pars.iniconf);
          break;
        }
    default:
        {
          //cout << "N=" << heMC.pars.N << "\n";
          heMC.createconf(0);
          break;
        }
     }
  //cout << "qui #steps="<< heMC.pars.steps <<"\n";
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
