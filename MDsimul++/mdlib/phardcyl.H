#include "./pmatrix.H"
#include "./rpoly.H"
#include "./pmatrix.H"
#include "./hardcyl.H"
template <class ntype>
using pvec3=pvector<ntype,3>;
template <class ntype>
using pmat3=pmatrixq<ntype,3>;
#define Sqr(x) ((x)*(x))
template <class ntype>
class phardcyl: public hardcyl<ntype> 
{
public:
  using hardcyl<ntype>::L;
  using hardcyl<ntype>::D;
  using hardcyl<ntype>::r;
  using hardcyl<ntype>::disk;
  using hardcyl<ntype>::n;
  using hardcyl<ntype>::shift;
  //using hardcyl<ntype>::overlap(class hardcyl<ntype> A, class hardcyl<ntype> B, pvector<ntype,3> shift);
  static const int NP;
  //pvector<ntype,3> rpaB[NP];// position of patches (body reference system)
  pvector<ntype,3> rpaL[NP];// position of patches (lab reference system)
  ntype sig[NP][NP]; // range of site-site interactions
  ntype u0[NP][NP];  // well depth of site-site interactions
  phardcyl(): hardcyl<ntype>()
    {
      update_patches();
    }
  phardcyl(ntype l, ntype d): hardcyl<ntype>(l,d)
    {
      update_patches();
    }
  void update_patches(void)
    {
      ntype L2=0.5*L;
      rpaL[0] = r + L2*n;
      rpaL[1] = r - L2*n;
    }
};
