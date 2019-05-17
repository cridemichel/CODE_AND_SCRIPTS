#include "../mdlib/pmatrix.H"
#include "../mdlib/rpoly.H"
template <class ntype>
using pvec3=pvector<ntype,3>;
template <class ntype>
using pmat3=pmatrixq<ntype,3>;
#define Sqr(x) ((x)*(x))
template <class ntype>
class hardell 
{
public:
  ntype a, b; //semiaxes
  pvector<ntype,2> na, nb, r;//na nb versori degli assi e r centro di massa 
  pvector<ntype,2> shift;
  template <class U>
    friend U overlap (class hardell<U>& A, class hardell<U>& B, pvector<U,2>& shift);
  void random_orient(void)
    {
      na.random_orient();
      nb << 0,1;
      if (nb==na)
        {
          nb << 1,0;
        }
      nb = nb - (nb*na)*na;
      nb = (1.0/nb.norm())*nb;
    }
  void random_box(ntype Lbox)
    {
      r.random_box();
      r = Lbox*r;
    }
  hardell()
    {
      na << 1,0;
      nb << 0,1;
      shift << 0,0;
      a =1;
      b =2;
      r << 0,0;
    }
};

template <class ntype>
ntype overlap (class hardell<ntype>& A, class hardell<ntype>& B)
{
  ntype x0, y0, alpha;
  pvector <ntype,5> coeffpa;
  ntype x02, y02, sai[2];
  ntype vt, at[2], dist, m00, m01, m11, b2d, m002, m012, m112;
  pvector<ntype,2> r0jpp, r0jp, ri, rj, x;
  pmatrixq<ntype,2> Di, Dj, Mi, Mj, Ri, Rj, Rjp, Mjp, Mjpp, Mtmp, RM, invM2I;
  pvector<complex<ntype>,4> roots;
  rpoly<ntype,4> quartic;
  int kk1, k1;
 
  /* apply affinity to reduce first ellipsoid to a sphere */
  Di << 1.0/Sqr(A.a), 0, 0, 1.0/Sqr(A.b);
  Dj << 1.0/Sqr(B.a), 0, 0, 1.0/Sqr(B.b);
  for (k1=0; k1 < 2; k1++)
    {
      Ri[0][k1] = A.na[k1];
      Rj[1][k1] = A.nb[k1];
    }
  /* sai[0] e sai[1] sono i semiassi della particella i */
  sai[0] = A.a;
  sai[1] = A.b;
  Mi = Ri.transp()*Di*Ri;
  Mj = Rj.transp()*Dj*Rj;
  ri = A.r;
  rj = B.r+B.shift;

  /* verifico che il centro di i non appartenga a j e viceversa come check preliminare */
  auto rji = rj-ri;
  auto rij = (-1.0)*rji;
  if ((rji*Mi)*rji-1.0 < 0.0)
    return -1.0;
  if ((rij*Mj)*rij-1.0 < 0.0)
    return -1.0;

  /* switch to ellipsoid i reference system */
  r0jp = Ri*rji;
  Rjp = Rj*Ri.transp();
  Mjp = Rjp.transp()*Dj*Rjp;

  /* calculate matrix and position of ellipsoid j after application of affinity
   * which reduces ellipsoid i to a sphere */   
  Mjpp[0][0] = Mjp[0][0]*Sqr(sai[0]);
  Mjpp[0][1] = Mjp[0][1]*sai[0]*sai[1];
  //Mjpp[1][0] = Mjp[1][0]*sai[0]*sai[1];
  Mjpp[1][1] = Mjp[1][1]*Sqr(sai[1]);
  r0jpp[0] = r0jp[0]/sai[0];
  r0jpp[1] = r0jp[1]/sai[1];
  x0 = r0jpp[0];
  y0 = r0jpp[1];
  x02 = Sqr(x0);
  y02 = Sqr(y0);
  m00 = Mjpp[0][0];
  m01 = Mjpp[0][1];
  m11 = Mjpp[1][1];
  m002 = Sqr(m00);
  m012 = Sqr(m01);
  m112 = Sqr(m11);
  b2d = -1.0 + m00*x02 + 2.0*m01*x0*y0 + m11*y02;

  coeffpa[0] = -Sqr(m012 - m00*m11)*(-b2d + m00*x02 + y0*(2*m01*x0 + m11*y0)); 
  coeffpa[1] = -2*(m00 + m11)*(-m01*m01 + m00*m11)*(-b2d + m00*x02 + 
   y0*(2*m01*x0 + m11*y0));
  coeffpa[2] =b2d*(m002 - 2*m012 + 4*m00*m11 + m112) - (m002*m00 - 2*m00*m012 + 
    4*m002*m11 + m012*m11)*x02 - 
    2*m01*(m002 - 3*m012 + 5*m00*m11 + m112)*x0*y0 - (-2*m012*m11 + 
    m11*m112 + m00*(m012 + 4*m112))*y02;
  coeffpa[3] = 2*b2d*(m00 + m11) - 2*(m002 + m012)*x02 - 
    4*m01*(m00 + m11)*x0*y0 - 2*(m012 + m112)*y02;
  coeffpa[4] = b2d;

  quartic.set_coeff(coeffpa);
  quartic.find_roots(roots);
  for (kk1=0; kk1 < 4; kk1++)
    {
      if (imag(roots[kk1])==0 && real(roots[kk1]) > 0.0)
        {
          alpha=real(roots[kk1]);
          break;
        }
    }
  auto detMa =-m01*m01 + (m00 + alpha)*(m11 + alpha);
  x[0] = x0*(-m01*m01 + m00*(m11 + alpha)) + y0*(-m01*m11 + m01*(m11 + alpha));
  x[1] = x0*(-m00*m01 + m01*(m00 + alpha)) + y0*(-m01*m01 + m11*(m00 + alpha));
  x = (1.0/detMa)*x;
  dist = x*x-1.0; 
  /* trasformando tramite l'affinità inversa i punti che individuano la distanza tra sfera ed ellissoide
   * si avrà la distanza tra i due ellissoidi che si può usare nella dinamica event-driven */
  if (dist < 0.0)
    return -1.0;
  else
    return 1.0;
}
