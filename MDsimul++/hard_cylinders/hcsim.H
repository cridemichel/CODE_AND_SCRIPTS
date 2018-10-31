#include "../mdlib/simulation.H"
#include "hard_cylinders.H"
class hcsim: public sim<double>
{
public:
  virtual void onestep(void);
  virtual void saveobserv(void);
  virtual void termiante(void);
}
