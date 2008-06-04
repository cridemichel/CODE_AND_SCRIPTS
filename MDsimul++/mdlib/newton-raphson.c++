#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// Class implementing a generic Newton-Raphson method to solve systems of equations
template <int n, class F> 
class NR: public linear-sys
{
  inline void SWAP(double a, double b) 
    {
      double temp;
      temp=a;
      a=b;
      b=temp;
    }
 
public:
  Vector<int n> FindZero(SysEq F);

}
