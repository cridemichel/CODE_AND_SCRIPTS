typedef long long int llint;
template <class ntype>
class params
{
  llint curstep, maxstep;
public:
  ntype D, L, u0, sigR, sigP;
  ntype rcutpa, rcutho;
  ntype Lbox[3];
  ntype pbias;
  ntype T; // temperature
  int npav, nexc;
  int ensemble; 
  params()
    {
    }
  ~params()
    {

    }
};
