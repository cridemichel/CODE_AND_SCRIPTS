#include "hcsim.H"
sim hcsim;
int main(int argc, char **argv)
{
  hcsim.parseargs(argc, argv);
  hcsim.run();
  return 0;
}
