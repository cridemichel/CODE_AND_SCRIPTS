#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
template <int n>
class simulation: public params, simstate
{
  void randomgen(void);
  void onestep(void);
  void geninitialconf(void);
  void linkedlist(void);
  void neighbourlist(void);
public:
  void setup(void);
  void iterate(void);
  void measure(void);
  void savesnap(void);
  void savestate(void);
}
