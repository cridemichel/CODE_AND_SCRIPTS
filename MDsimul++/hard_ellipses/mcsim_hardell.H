#include "../mdlib/pvector.H"
template<class ntype>
class params
{

public:
  int N; // number of particles
  pvector<ntype,2> L; // box
  pvector<ntype,-1> Veff; //effective potential
  int Neff; // points of effective potential to use for interpolation
  pvector<ntype,2> M; // linked list cells 
  void setveff(int n, ntype* v)
    {
      Veff.allocate(n);
      Neff = n;
      for (auto i=0; i < n; i++)
        {
          Veff[i] = v[i];
        }
    }
  params()
    {
      N=1000;
      L[0]=L[1]=L[2]=10.0;
    }
};
