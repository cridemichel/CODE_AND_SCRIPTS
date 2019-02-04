// base container class for particles
#include "./pvector.H"
template <class ntype>
class particle 
{
public:
 pvector<ntype,3> r; // center of mass position
 virtual ntype overlap(particle<ntype> A, particle<ntype> B);
}
