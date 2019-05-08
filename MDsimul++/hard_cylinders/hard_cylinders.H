#include "../mdlib/pmatrix.H"
#include "../mdlib/simulation.H"
#include "../mdlib/rpoly.H"
template <class ntype>
class hard_cylinder: public particle
{
  ntype rimdiskone_hybrid(ntype D, ntype L, pvector<ntype,3>& Ci, pvector<ntype,3> ni, pvector<ntype,3> Dj, pvector<ntype,3> nj);
public:
  ntype overlap(vector<particle<ntype>>& parB, pvector<ntype,3>& shift);
  
};
