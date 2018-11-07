#include "../mdlib/pmatrix.H"
#include "../mdlib/simulation.H"

template <class ntype>
class hard_cylinder: public particle
{

public:
  ntype overlap(vector<particle<ntype>>& parB, pvector<ntype,3>& shift);
  
};
