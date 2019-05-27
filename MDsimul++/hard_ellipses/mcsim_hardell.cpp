#include "./mcsim_hardell.H"
int main (int argc, char **argv)
{
  char *rf;
  mcsim<double> heMC;
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
          cout << "Restarting from file " << string(rf) << "\n";;
          heMC.restart(rf);
          cout << " step = " << heMC.iniStep << "\n";
          break;
        }
    case 2:
        {
          cout << "Using snapshot " << heMC.pars.iniconf << " as initial configuration\n";
          heMC.readconf(heMC.pars.iniconf);
          break;
        }
    default:
        {
          //cout << "N=" << heMC.pars.N << "\n";
          cout << "Creating initial conf\n";
          heMC.createconf(0);
          break;
        }
     }
  cout << "initial phi=" << heMC.calc_phi() << "\n";
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
