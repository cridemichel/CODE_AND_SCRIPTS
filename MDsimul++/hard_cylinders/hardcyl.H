#include "../mdlib/pmatrix.H"
#include "../mdlib/rpoly.H"
template <class ntype>
using pvec3=pvector<ntype,3>;
template <class ntype>
using pmat3=pmatrixq<ntype,3>;
#define Sqr(x) ((x)*(x))

#if 0
template <class ntype>
class particle
{

};
#endif
#if 0
template <class ntype> class hardcyl;

template <class ntype>
ntype overlap(class hardcyl<ntype> A, class hardcyl<ntype> B, pvector<ntype,3> shift);
#endif
template <class U>
int test_for_fallback(pvec3<U> P, pvec3<U> Cip, pvec3<U> nip, U D2i, U D2j, U& diff);
template <class ntype>
void versor_to_R_alt_fb(pvector<ntype,3> Ci,  pvector<ntype,3> ni, 
                        pvector <ntype,3> Dj, pvector<ntype,3> nj, 
                        pmatrixq<ntype,3>& R, ntype D, pvector<ntype,3> Tj, 
                        int MAXITS)
{
  int k, its;
  ntype DjCini, AiDjnj, ragg;
  pvector<ntype,3> u, dsc, DjCi, Ai, AiDj, Tnew, VV, dscperp, dscpara, TnCi;
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

  /* N.B. Questa viene utilizzata in fallback e faccio un ulteriore step dell'algorito di Ibarra semplificato per
   * determinare l'asse y del riferimenti del disco (l'asse x è l'asse perpendicolare
   * al disco e l'asse z si ottiene con il prodotto vettore dell'asse x e y) */
  for (its = 0; its < MAXITS; its++)
    {
      DjCi = Dj - Ci;
      //normDjCi = DjCi.norm();
      DjCini = DjCi.dot(ni);

      Ai = Ci + DjCini*ni;
      AiDj = Ai - Dj; 
      AiDjnj = AiDj.dot(nj);
      VV = AiDj - AiDjnj*nj;
      //for (kk1=0; kk1 < 3; kk1++)
      //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
      //VV.show("VV=");
      ragg = 1.0/VV.norm();
      //for (kk1=0; kk1 < 3; kk1++)
      VV = ragg*VV;
      //cout << "boh=" << VV.norm() << "\n";
      Tj = (D*0.5)*VV;
    }

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k];
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
  u=nj^VV;
  //cout << "norm(VV)=" << VV.norm() << "\n";
  //cout << "norm(u)=" << u.norm() << "nj.u=" << nj.dot(u) << "\n";
 
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}

template <class U>
void rotate_axes_on_plane(pmat3<U>& RR);

template <class U> 
U perpcomp(pvec3<U> V, pvec3<U> C, pvec3<U> n);
 
template <class ntype>
class hardcyl
{
  int ii;
  template <class U>
  friend U rimrim(class hardcyl<U>& A, class hardcyl<U>& B);
  template <class U>
  friend U diskdisk(class hardcyl<U>& A, class hardcyl<U>& B);
  template <class U>
  friend U rimdisk(class hardcyl<U>& A, class hardcyl<U>& B);

  template <class U>
  friend U test_overlap_parall_cyl(class hardcyl<U>& A, class hardcyl<U>& B, int wd)
    {
      U DjCini, normDjUi;
      pvec3<U> DjCi, Ui, DjUi;
      DjCi = (B.disk[wd]+B.shift) - (A.r+A.shift);
      //normDjCi = calc_norm(DjCi);
      DjCini = DjCi*A.n;

      Ui = (A.r+A.shift) + DjCini*A.n;
      DjUi = B.disk[wd]+B.shift - Ui;
      normDjUi = abs(DjUi);

      if (normDjUi <= 0.5*(A.D+B.D) && abs(DjCini) <= A.L*0.5)
        return -1.0;
      else
        return 1.0;
    }
  template <class U>
  friend int test_for_fallback(pvec3<U> P, pvec3<U> Cip, pvec3<U> nip, U D2i, U D2j, 
                        U& diff)
    {
      const U DIST_THR=5E-13;// old value 5.0E-13 -- NOTA 22/02/18 con il nuovo quartic solver ci si può spingere a 5E-14 volendo!
      U diff1, diff2;
      diff1=abs(perpcomp(P, Cip, nip)-D2i); // qui D2 è il diametro del rim
      diff2=abs(sqrt(Sqr(P[1])+Sqr(P[2]))-D2j);// qui D2 è il diametro del disco

      diff=diff1+diff2;
      if (diff1 > DIST_THR*D2i || diff2 > DIST_THR*D2j)
        return 1;
      else
        return 0;
    }
  template <class U> 
  friend U perpcomp(pvec3<U> V, pvec3<U> C, pvec3<U> n)
    {
      pvec3<U> dsc, dscperp;
      U sp;
      dsc = V - C; 
      sp = dsc*n;
      dscperp = dsc-sp*n;
      return abs(dscperp);
    }
  template <class U>
  friend void rotate_axes_on_plane(pmat3<U>& RR)
    {
      pmat3<U> Rin, OmegaSq, Omega, M;
      U ox, oy, oz, theta, sinw, cosw;

      ox = RR[0][0];
      oy = RR[0][1];
      oz = RR[0][2];
      Rin=RR;
      theta = M_PI/4.0;
      //thetaSq=Sqr(theta);
      sinw = sin(theta);
      cosw = (1.0 - cos(theta));
      Omega[0][0] = 0;
      Omega[0][1] = -oz;
      Omega[0][2] = oy;
      Omega[1][0] = oz;
      Omega[1][1] = 0;
      Omega[1][2] = -ox;
      Omega[2][0] = -oy;
      Omega[2][1] = ox;
      Omega[2][2] = 0;
      OmegaSq[0][0] = -Sqr(oy) - Sqr(oz);
      OmegaSq[0][1] = ox*oy;
      OmegaSq[0][2] = ox*oz;
      OmegaSq[1][0] = ox*oy;
      OmegaSq[1][1] = -Sqr(ox) - Sqr(oz);
      OmegaSq[1][2] = oy*oz;
      OmegaSq[2][0] = ox*oz;
      OmegaSq[2][1] = oy*oz;
      OmegaSq[2][2] = -Sqr(ox) - Sqr(oy);
      M = -sinw*Omega+cosw*OmegaSq;
      RR = Rin*M;
    }

public:
  ntype L, D;
  pvector<ntype,3> disk[2];
  pvector<ntype,3> n, r;
  pvector<ntype,3> shift;
  //templated function overlap is a friend of class hardcyl
  //note that we if we use the following syntax
  //friend ntype overlap (class hardcyl<ntype> A, class hardcyl<ntype> B, pvector<ntype,3> shift);
  //the compiler will not instantiate overlap when instatiating the class and we will get a linker error
  //since it is a friend function and not a method!
  template <class U>
  friend U overlap (class hardcyl<U>& A, class hardcyl<U>& B, pvector<U,3>& shift);
  void random_orient(void)
    {
      n.random_orient();
      //update_disks();
    }
  void random_box(ntype Lbox)
    {
      r.random_box();
      r = Lbox*r;
    }
  void update_disks(void)
    {
      ntype Lh=0.5*L;
      disk[0] = r+Lh*n;
      disk[1] = r-Lh*n;
    }
  hardcyl()
    {
      L=1.0;
      D=1.0;
      n << 1,0,0;
      r << 0,0,0;
      update_disks();
    }
  hardcyl(ntype l, ntype d)
    {
      L=l;
      D=d;
      n << 1,0,0;
      r << 0,0,0;
      update_disks();
     }
};
template <class ntype>
ntype rimrim(hardcyl<ntype>& A, hardcyl<ntype>& B)
{
  ntype lambdai, lambdaj, ninj;
  ntype CiCjni, CiCjnj, detA; 
  pvector<ntype,3> CiCj, ViVj, Vi, Vj;
  /* case A.3 rim-rim overlap */
  CiCj = A.r - (B.r+B.shift);
  CiCjni = CiCj*A.n;
  CiCjnj = CiCj*B.n;
  ninj = A.n*B.n;
  detA = ninj*ninj-1.0;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  Vi = A.r + lambdai*A.n;
  Vj = (B.r+B.shift) + lambdaj*B.n;
  ViVj = Vi-Vj;
  if (ViVj.norm() < 0.5*(A.D+B.D) && abs(lambdai) < 0.5*A.L && abs(lambdaj) < 0.5*B.L)
    {
      return -1;
    }
  return 0;
}
template <class ntype>
ntype diskdisk(hardcyl<ntype>& A, hardcyl<ntype>& B)
{
  int j1, j2, kk; 
  ntype sp, normCiCj, Q1;
  pvector <ntype,3> CiCj, VV, niN, njN, PjDj, assex, Pi, N, Pj, PiPj, PiDi;
  ntype DiN, DjN, Djnj, Dini, normN;
  ntype normNSq, normPiPj, Q2, normPiDi, normPjDj;
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  CiCj = A.r - (B.r + B.shift);
  if (A.n==B.n)
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = CiCj.norm();
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (VV*A.n==1.0)
	{
	  if (normCiCj <= 0.5*(A.L+B.L))
	    return -1;
	  else
	    return 1;
	}

      /* parallel disks */
      for (j1=0; j1 < 2; j1++)
	for (j2=j1; j2 < 2; j2++)
	  {
            VV = A.disk[j1]-(B.disk[j2]+B.shift);
            sp=A.n*VV;
	    if (sp == 0 && VV.norm() < 0.5*(A.D+B.D))
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      N=A.n^B.n;// ^ is the inner product
      niN=A.n^N;
      njN=B.n^N;
      normN=N.norm();
      normNSq=normN*normN;
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = A.disk[j1]*N;
	    DjN = B.disk[j2]*N;
	    Dini = A.disk[j1]*A.n;
	    Djnj = B.disk[j2]*B.n;
            Pi = (DiN*N+Dini*njN-Djnj*niN)/normNSq;
            Pj = (DjN*N+Dini*njN-Djnj*niN)/normNSq;
            PiDi = Pi - A.disk[j1];
            PjDj = Pj - (B.disk[j2]+B.shift); 
	    normPiDi = PiDi.norm();
	    normPjDj = PjDj.norm();
	    if (normPiDi <= 0.5*A.D && normPjDj <= 0.5*B.D)
	      {
		Q1 = sqrt((A.D*A.D)/4.0-(normPiDi*normPiDi));
		Q2 = sqrt((B.D*B.D)/4.0-(normPjDj*normPjDj));
                PiPj = Pi-Pj;
		normPiPj = PiPj.norm();
		if (normPiPj <= Q1 + Q2)
		  {
		    return -1;
		  }
		//else 
		//return 1;
	      }
	    //else 
	    //return 1;
	  }
    }
  return 0;
}
template <class ntype>
void discard_spurious(ntype *solqua, int *numsol)
{
  /* each solution x must be such that |x| <= 1 */
  int k, nsol;
  const ntype EPS=1E-9;
  ntype solL[4];
  nsol=0;
  for (k=0; k < *numsol; k++)
    solL[k] = solqua[k];
  for (k=0; k < *numsol; k++)
    {
      if (fabsl(solL[k]) < 1.0+EPS)
	{
	  solqua[nsol] = solL[k];
	  nsol++;
	}
    }
  //printf("numsol=%d nsol=%d\n", *numsol, nsol);
  *numsol = nsol;
}
template <class ntype>
void solve_quartic(pvector<ntype,5>& coeff, int *numsol, ntype *solqua)
{
  int k;
  rpoly<ntype,4> rp;
  pvector<complex<ntype>,4> roots;
  rp.set_coeff(coeff);
  rp.find_roots(roots);

  *numsol=0;
  for (k=0; k < 4; k++)
    {
      //cout << "k=" << k << "\n";
      if (imag(roots[k])==0)
	{
	  solqua[*numsol] = real(roots[k]);
	  (*numsol)++;
	}
    }
}
template <class ntype>
ntype rimdiskone(hardcyl<ntype>& A, hardcyl<ntype>& B, int wd)
{
  int kk1, numsol[2]={0,0}, fallback, solset;
  ntype tmp, sp, solec[4][2], solqua[4];
  ntype c0, c1, c2, c3, c02, c12, c22, coeffEr[6], c32, c42, c52, c4, c5;  
  ntype diff[2][4], maxdiff[2], sumdiff[2];
  ntype nip02,nip12,nip22,Cip02,Cip12,Cip22, temp, normv;
  ntype omnip02, omnip12, omnip22;
  ntype D2isq, D2i, D2jsq, D2j, Cip0, Cip1, Cip2, nip0, nip1 , nip2, nip1nip2, nip0nip2, nip0nip1; 
  pvector<ntype,3> nip[2], Cip[2], Ci;
  pmatrixq<ntype,3> Rl;
  pvector<ntype,3> dsc, dscperp, uy;
  pvector<ntype,5> coeff;
  pvec3<ntype> solarr[2][4];
  /* rim of A and disk[wd] of B */
  D2i = A.D*0.5;
  D2j = B.D*0.5; 
  D2isq = D2i*D2i;
  D2jsq = D2j*D2j;
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  uy = B.disk[wd]+B.shift;
  //uy.show("uy");
  // verificare se D è di B o A!!! 
  versor_to_R_alt_fb(A.r+A.shift, A.n, B.disk[wd]+B.shift, B.n, Rl, B.D, uy, 1); 
  //versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  nip[0] = Rl*A.n;
  //(Rl*Rl.transp()).show("R");
  Cip[0] = Rl*((A.r+A.shift)-(B.disk[wd]+B.shift));
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  normv = nip[0].norm();
  nip0 = nip[0][0]/normv;
  nip1 = nip[0][1]/normv;
  nip2 = nip[0][2]/normv;
  Cip0 = Cip[0][0];
  Cip1 = Cip[0][1];
  Cip2 = Cip[0][2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);
  /* with some simplifications we save a bunch of FLOPS... */
  omnip02 = 1.0 - nip02;
  omnip12 = 1.0 - nip12;
  omnip22 = 1.0 - nip22;
  nip1nip2 = nip1*nip2;
  nip0nip2 = nip0*nip2;
  nip0nip1 = nip0*nip1;
  coeffEr[0] = omnip12;
  coeffEr[1] = omnip22;
  coeffEr[2] = -2.0*nip1nip2;  
  coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
								    Cip1*Cip2*nip1nip2) - D2jsq;
  coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
  coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  
  /* check ellipse */
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2jsq;
  coeffEr[1] *= D2jsq; 
  coeffEr[2] *= D2jsq;
  coeffEr[4] *= D2j;
  coeffEr[5] *= D2j;
  c0 = coeffEr[0];
  c1 = coeffEr[1];
  c2 = coeffEr[2];
  c3 = coeffEr[3];
  c4 = coeffEr[4];
  c5 = coeffEr[5];
  c02 = Sqr(c0);
  c12 = Sqr(c1);
  c22 = Sqr(c2);
  c32 = Sqr(c3);
  c42 = Sqr(c4);
  c52 = Sqr(c5);
#ifndef MC_EXCHG_QUART_SOL
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
#else
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
  coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
  coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
  coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
  //coeff.show("[C++] quartic coeff=");
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(A, B, wd);
    }
  else
    {
      solve_quartic(coeff, &(numsol[0]), solqua);
    }
  //discard_spurious(solqua, &(numsol[0]));
  fallback = 0;

  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      //  cout << "solqua[" << kk1 << "]=" << solqua[kk1] << "\n";
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
#if 0
      temp = c5 + c2*solqua[kk1];
      solec[kk1][0] = solqua[kk1];
      solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif     
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
      if (temp==0.0) 
	{
	  fallback=1;
	}
    }
  
  //cout << "numsol[0]=" << numsol[0] << "\n";
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2j*solec[kk1][0];
      solarr[0][kk1][2] = D2j*solec[kk1][1];
      if (test_for_fallback(solarr[0][kk1], Cip[0], nip[0], D2i, D2j, diff[0][kk1]))
	{
	  fallback=1;
#if 0
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[0][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[0][kk1], Cip[0], nip[0]));

	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
    }
#if 0
  if (tinyimagGBL)
    {
      //printf("qui\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
#endif
  solset=0;
  if (fallback)
    {
      rotate_axes_on_plane(Rl);
      nip[1] = Rl*A.n;
      Cip[1] = Rl*(A.r+A.shift-(B.disk[wd]+B.shift));
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
      ntype normv= nip[1].norm();
      nip0 = nip[1][0]/normv;
      nip1 = nip[1][1]/normv;
      nip2 = nip[1][2]/normv;
      Cip0 = Cip[1][0];
      Cip1 = Cip[1][1];
      Cip2 = Cip[1][2];
      nip02=Sqr(nip0);
      nip12=Sqr(nip1);
      nip22=Sqr(nip2);
      //nip04=Sqr(nip02);
      //nip14=Sqr(nip12);
      //nip24=Sqr(nip22);
      //nip03=nip02*nip0;
      //nip13=nip12*nip1;
      //nip23=nip22*nip2;
      Cip02=Sqr(Cip0);
      Cip12=Sqr(Cip1);
      Cip22=Sqr(Cip2);   
#if 1
      omnip02 = 1.0 - nip02;
      omnip12 = 1.0 - nip12;
      omnip22 = 1.0 - nip22;
      nip1nip2 = nip1*nip2;
      nip0nip2 = nip0*nip2;
      nip0nip1 = nip0*nip1;
      coeffEr[0] = omnip12;
      coeffEr[1] = omnip22;
      coeffEr[2] = -2.0*nip1nip2;  
      coeffEr[3] = Cip02*omnip02 + Cip12*omnip12 + Cip22*omnip22 - 2.0*(Cip0*Cip1*nip0nip1 + Cip0*Cip2*nip0nip2 +
									Cip1*Cip2*nip1nip2) - D2isq;
      coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
      coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  

#elif 0
      coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
	nip12*nip22;
      coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
	nip12*nip22 + nip24;
      coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
	2*nip1*nip23;
      coeffEr[3] = Cip02 + Cip12 + Cip22 - D2sq - 
	2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1 - 
	2*Cip12*nip12 + Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 + Cip12*nip14 - 
	4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2 - 4*Cip1*Cip2*nip1*nip2 + 2*Cip1*Cip2*nip02*nip1*nip2 + 
	2*Cip0*Cip2*nip0*nip12*nip2 + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22 + Cip02*nip02*nip22 + 
	Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + Cip12*nip12*nip22 + Cip22*nip12*nip22 + 
	2*Cip0*Cip2*nip0*nip23 + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;
      coeffEr[4] = -2*Cip1 + 4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 + 
	4*Cip1*nip12 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13 - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 
	2*Cip2*nip02*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip0*nip0*nip1*nip22 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;
      coeffEr[5] = -2*Cip2 + 4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 + 
	4*Cip1*nip1*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22 - 2*Cip2*nip02*nip22 - 2*Cip2*nip12*nip22 - 2*Cip0*nip0*nip23 - 2*Cip1*nip1*nip23 
	- 2*Cip2*nip24;
#else
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/

      coeffEr[0] = 1.0 + ( -2*nip12 + nip14 + nip12*nip22) + nip02*nip12;
      coeffEr[1] = 1.0 + ( -2*nip22 + nip12*nip22 + nip24) + nip02*nip22;
      coeffEr[2] = 2*nip02*nip1*nip2 + (- 4*nip1*nip2 + 2*nip13*nip2 + 
					2*nip1*nip23);

      coeffEr[3] = 
	(- 2*Cip02*nip02 + Cip02*nip04 - 4*Cip0*Cip1*nip0*nip1 + 2*Cip0*Cip1*nip03*nip1+ 
	 Cip02*nip02*nip12 + Cip12*nip02*nip12 + 2*Cip0*Cip1*nip0*nip13 - 4*Cip0*Cip2*nip0*nip2 + 2*Cip0*Cip2*nip03*nip2
	 + 2*Cip1*Cip2*nip02*nip1*nip2 + 2*Cip0*Cip2*nip0*nip12*nip2 + Cip02*nip02*nip22 + 
	 Cip22*nip02*nip22 + 2*Cip0*Cip1*nip0*nip1*nip22 + 2*Cip0*Cip2*nip0*nip23 ) 
	+ Cip02 + Cip12 + Cip22 - Sqr(D2)  - 
	2*Cip12*nip12  + Cip12*nip14 - 4*Cip1*Cip2*nip1*nip2  + 2*Cip1*Cip2*nip13*nip2 - 2*Cip22*nip22  + Cip12*nip12*nip22 
	+ Cip22*nip12*nip22  + 2*Cip1*Cip2*nip1*nip23 + Cip22*nip24;

      coeffEr[4] =
	(4*Cip0*nip0*nip1 - 2*Cip0*nip03*nip1 +  
	 - 2*Cip1*nip02*nip12 - 2*Cip0*nip0*nip13
	 - 2*Cip2*nip02*nip1*nip2 - 2*Cip0*nip0*nip1*nip22 ) 
	- 2*Cip1 + 4*Cip1*nip12  - 2*Cip1*nip14 + 4*Cip2*nip1*nip2 - 2*Cip2*nip13*nip2 - 2*Cip1*nip12*nip22 - 
	2*Cip2*nip1*nip23;

      coeffEr[5] = 
	(4*Cip0*nip0*nip2 - 2*Cip0*nip03*nip2 - 2*Cip1*nip02*nip1*nip2 - 2*Cip0*nip0*nip12*nip2 - 2*Cip2*nip02*nip22
	 - 2*Cip0*nip0*nip23 ) -2*Cip2 + 4*Cip1*nip1*nip2  - 2*Cip1*nip13*nip2 + 
	4*Cip2*nip22  - 2*Cip2*nip12*nip22 - 2*Cip1*nip1*nip23 - 2*Cip2*nip24;

#endif
      /* check ellipse */

      /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
      coeffEr[0] *= D2jsq;
      coeffEr[1] *= D2jsq; 
      coeffEr[2] *= D2jsq;
      coeffEr[4] *= D2j;
      coeffEr[5] *= D2j;
      c0 = coeffEr[0];
      c1 = coeffEr[1];
      c2 = coeffEr[2];
      c3 = coeffEr[3];
      c4 = coeffEr[4];
      c5 = coeffEr[5];
      c02 = Sqr(c0);
      c12 = Sqr(c1);
      c22 = Sqr(c2);
      c32 = Sqr(c3);
      c42 = Sqr(c4);
      c52 = Sqr(c5);

#if 0
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
#endif
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
      coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
      coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
      coeff[0] = c02 + 2*c0*c3 + c32 - c42;
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(A, B, wd);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      //discard_spurious(solqua, &(numsol[1]));
#ifdef MC_QUART_VERBOSE
      printf("falling back [#%ld] type=%d numsol=%d %d\n", numfb++,fallback, numsol[0], numsol[1]);
#endif
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0 
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
#endif
	  temp = c4 + c2*solqua[kk1];
	  if (temp==0)
	    {
	      printf("[WARNING] temp is 0 in fallback hybrid numsol=%d %d\n", numsol[0], numsol[1]);
	    }
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2j*solec[kk1][0];
	  solarr[1][kk1][2] = D2j*solec[kk1][1];
	  test_for_fallback(solarr[1][kk1], Cip[1], nip[1], D2i, D2j, diff[1][kk1]);
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  
	}
      if (fallback==2)
	solset=1;
      else if (numsol[1]==0 && numsol[0] > 0)
	solset=0;
      else
	{
	  if (maxdiff[1] < maxdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
	}
    }
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
      dsc = solarr[solset][kk1] - Cip[solset];
      if (test_for_fallback(solarr[solset][kk1], Cip[solset], nip[solset], D2i, D2j, tmp)) 
	{
	  printf("[C++] # %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("[C++] distanza punto-centro disk: %.15G\n", solarr[solset][kk1].norm());
	  printf("[C++] distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip[solset], nip[solset]));

	  if (kk1 < numsol[1-solset])
	    {
	      printf("[C++] DISCARDED SOLSET [%d]\n", 1-solset);
	      printf("[C++] distanza punto-centro disk: %.15G\n", abs(solarr[1-solset][kk1]));
	      printf("[C++] distanza punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip[1-solset], nip[1-solset]));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("[C++] distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("[C++] numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("[C++] ni=",ni);
	  print_vec("[C++] nj=",nj);
	  printf("[C++] c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("[C++] c4=%.15G c5=%.15G\n", c4, c5);
	  printf("[C++] solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("[C++] coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("[C++] solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("[C++] ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("[C++] (%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("[C++]{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("[C++] quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("[C++] temp=%.15G\n", temp);
#endif
	  printf("[C++] # %d >>>> =====================  \n", kk1);
	}
      sp = dsc*nip[solset];
      if (abs(sp) < A.L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
template <class ntype>
ntype rimdisk(hardcyl<ntype>& A, hardcyl<ntype>& B)
{
  ntype sp[2], perpdist[2];
  pvec3<ntype> DjCi;
  ntype DjCini;
  int whichdisk;
  /* We check overlap between the two disks of B and rim of A */
  /* choose the disk of B closest to axis of cylinder A */
  DjCi = B.disk[0] - A.r;
  sp[0] = DjCi*A.n;
  perpdist[0]=(DjCi - sp[0]*A.n).norm();
  DjCi = B.disk[1] - A.r;
  sp[1] = DjCi*A.n;
  perpdist[1]=(DjCi - sp[1]*A.n).norm();
  whichdisk=(perpdist[0]<perpdist[1])?0:1;
  DjCini=sp[whichdisk];
  //Ui = A.r + DjCini*A.n;
  //normDjUi = perpdist

  if (perpdist[whichdisk] < 0.5*A.D && abs(DjCini) <= 0.5*A.L)
    return -1;

  if ( !(perpdist[whichdisk] > 0.5*(A.D+B.D)) &&
       !(perpdist[whichdisk] < 0.5*A.D && abs(DjCini) > 0.5*A.L) )
    {
    
      if (rimdiskone(A, B, whichdisk) < 0.0)
        return -1;
    }
  /* Swapped A and B, i.e. we check overlap between the two disks of A and rim of B */
  /* choose the disk of A closest to axis of cylinder B */
  DjCi = A.disk[0] - B.r;
  sp[0] = DjCi*B.n;
  perpdist[0]=(DjCi - sp[0]*B.n).norm();
  DjCi = A.disk[1] - B.r;
  sp[1] = DjCi*B.n;
  perpdist[1]=(DjCi - sp[1]*B.n).norm();
  whichdisk=(perpdist[0]<perpdist[1])?0:1;
  DjCini = sp[whichdisk];

  if (perpdist[whichdisk] < 0.5*B.D && abs(DjCini) <= 0.5*B.L)
    return -1;

  if ( !(perpdist[whichdisk] > 0.5*(A.D+B.D)) &&
       !(perpdist[whichdisk] < 0.5*B.D && abs(DjCini) > 0.5*B.L) )
    {
      if (rimdiskone(B, A, whichdisk) < 0.0)
        return -1;
    }
  return 0;
}
template <class ntype>
ntype overlap(class hardcyl<ntype>& A, class hardcyl<ntype>& B, pvector<ntype,3>& shift)
{
  pvec3<ntype> CiCj, pv;
  ntype sp, perpc, parac;
  A.shift << 0,0,0;
  ntype ret;
  B.shift = shift;
  A.ii =0;
  //cout << "checking overlap\n";
  //numcallsRD++;
  /* if cyls are parallel is simpler... */
  if (1.0-abs(A.n*B.n) == 0)
    //if (1.0-fabs(scalProd(ni,nj)) < 1.0E-8)
    {
      CiCj = (A.r+A.shift) - (B.r+B.shift);
      sp = CiCj*A.n;
      pv = sp*A.n;
      parac = abs(pv);
      perpc = abs(CiCj - pv);
      if (perpc <= (A.D+B.D)*0.5 && parac <= (A.L+B.L)*0.5)
        return -1.0;
      else 
        return 1.0;
    }
  if ((ret=diskdisk(A,B))!=0)
    {
      return ret;
    } 
  if ((ret=rimrim(A,B))!=0)
    {
      return ret;
    }

  if ((ret=rimdisk(A,B))!=0)
    {
      return ret;
    } 
  return 1.0;
}
