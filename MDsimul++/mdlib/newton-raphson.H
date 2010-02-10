#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// Class implementing a generic Newton-Raphson method to solve systems of equations
template <int n>
class NR: public linear-sys
{
  inline void SWAP(double a, double b) 
    {
      double temp;
      temp=a;
      a=b;
      b=temp;
    }
  /* costruttore e distruttore */
  void NR(class NR-sys);
  void ~NR(void);
  Vector x;
  Vector fvec;
public:
  Vector<int n> FindZero(class NR-sys);
}
