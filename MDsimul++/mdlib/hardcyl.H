#include "./pmatrix.H"
#include "./rpoly.H"
template <class ntype>
using pvec3=pvector<ntype,3>;
#define Sqr(x) ((x)*(x))

template <class ntype>
class hardcyl
{
  friend double rimrim(class hardcyl<ntype> A, class hardcyl<ntype> B);
  friend double diskdisk(class hardcyl<ntype> A, class hardcyl<ntype> B);
  friend double rimdisk(class hardcyl<ntype> A, class hardcyl<ntype> B);
public:
  ntype L, D;
  pvector<ntype,3> disk[2];
  pvector<ntype,3> n;
  pvector<ntype,3> r;
  pvector<ntype,3> shift;
  friend double overlap(class hardcyl<ntype> A, class hardcyl<ntype> B, pvector<ntype,3> shift);
  pvector<ntype,3> update_disks(void)
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
double rimrim(hardcyl<ntype> A, hardcyl<ntype> B)
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
double diskdisk(hardcyl<ntype> A, hardcyl<ntype> B)
{
  int j1, j2, kk; 
  ntype sp, normCiCj, Q1;
  pvector <ntype,3> CiCj, VV, niN, njN, PjDj, assex, Pi, N, Pj, PiPj, PiDi;
  ntype DiN, DjN, Djni, Djnj, Dini, normN;
  ntype normNSq, normPiPj, Q2, normPiDi, normPjDj;
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  CiCj = A.r - (B.r + B.shift);
  if (A.n==B.n)
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
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
void solve_quartic(pvector<ntype,5> coeff, int *numsol, ntype *solqua)
{
  int k;
  rpoly<ntype,4> rp;
  pvector<complex<ntype>,4> roots;
  rp.set_coeff(coeff);
  rp.find_roots(roots);
  for (k=0; k < 4; k++)
    {
      if (imag(roots[k])==0)
	{
	  solqua[*numsol] = real(roots[k]);
	  (*numsol)++;
	}
    }

}
template <class ntype>
ntype rimdiskone(hardcyl<ntype> A, hardcyl<ntype> B, int wd)
{
  int kk1, kk2, numsol[2], nsc, fallback, solset;
  const ntype FALLBACK_THR = 1E-4;
  ntype tmp, sp, solec[4][2], solqua[4], solquasort[4], solquad[2];
  ntype dsc[3], dscperp[3], c0, c1, c2, c3, c02, c12, c22, nipp[3], Cipp[3], coeffEr[6], rErpp1sq, rErpp2sq, c32, c42, c52, c4, c5;  
  ntype diff[2][4], maxdiff[2], sumdiff[2], diffxy[2][4];
  ntype nip02,nip12,nip22,Cip02,Cip12,Cip22, temp;
  ntype omnip02, omnip12, omnip22;
  ntype D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2, nip1nip2, nip0nip2, nip0nip1; 
  pvector<ntype,3> nip[2], Cip[2], uy;
  pmatrixq<ntype,3> Rl;
  pvector<ntype,5> coeff;
  pvec3<ntype> solarr[2][4];
  /* rim of A and disk[wd] of B */
  D2 = B.D*0.5; 
  D2sq = D2*D2;
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  uy = B.disk[wd];
  //versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 1); 
  //versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  nip[0] = Rl*A.n;
  Cip[0] = Rl*((A.r+A.shift)-(B.disk[wd]+B.shift));

  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = abs(nip[0]);
  nip0 = nip[0][0]/norm;
  nip1 = nip[0][1]/norm;
  nip2 = nip[0][2]/norm;
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
								    Cip1*Cip2*nip1nip2) - D2sq;
  coeffEr[4] = 2.0*(Cip2*nip1nip2 + Cip0*nip0nip1 - Cip1*omnip12);
  coeffEr[5] = 2.0*(Cip0*nip0nip2 + Cip1*nip1nip2 - Cip2*omnip22);  
  /* check ellipse */
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= D2sq;
  coeffEr[1] *= D2sq; 
  coeffEr[2] *= D2sq;
  coeffEr[4] *= D2;
  coeffEr[5] *= D2;
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
  if (coeff[4]==0)
    {
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
    }
  else
    {
      solve_quartic(coeff, &(numsol[0]), solqua);
    }
  //discard_spurious(solqua, &(numsol[0]));
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
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
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = D2*solec[kk1][0];
      solarr[0][kk1][2] = D2*solec[kk1][1];
      if (test_for_fallback(solarr[0][kk1], Cip[0], nip[0], D2, &(diff[0][kk1])))
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
      for (kk1=0; kk1 < 3; kk1++)
	{
	  nip[1][kk1] = 0;
	  //Aip[kk1] = 0;
	  Cip[1][kk1] = 0;
	  for (kk2=0; kk2 < 3; kk2++)
	    {
	      nip[1][kk1] += Rl[kk1][kk2]*ni[kk2];
	      Cip[1][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	      //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	    } 
	}
      /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
      norm = calc_norm(nip[1]);
      nip0 = nip[1][0]/norm;
      nip1 = nip[1][1]/norm;
      nip2 = nip[1][2]/norm;
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
									Cip1*Cip2*nip1nip2) - D2sq;
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
      coeffEr[0] *= D2sq;
      coeffEr[1] *= D2sq; 
      coeffEr[2] *= D2sq;
      coeffEr[4] *= D2;
      coeffEr[5] *= D2;
      //printf("coeffEr=%.15G %.15G\n", coeffEr[0], coeffEr[1]);
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
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, L, D, D);
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
	  //printf("tempnew=%.16G\n", temp);
	  if (temp==0)
	    {
	      printf("[WARNING] temp is 0 in fallback hybrid numsol=%d %d\n", numsol[0], numsol[1]);
	    }
	  solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
	  solec[kk1][1] = solqua[kk1];
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = D2*solec[kk1][0];
	  solarr[1][kk1][2] = D2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
	  test_for_fallback(solarr[1][kk1], Cip[1], nip[1], D2, &(diff[1][kk1]));
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
	    //if (sumdiff[1] < sumdiff[0])
	    solset = 1;
	  else 
	    solset = 0;
	}
    }
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
      dsc = solarr[solset][kk1] - Cip[solset];
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallback(solarr[solset][kk1], Cip[solset], nip[solset], D2, &tmp)) 
	{
	  printf("# %d numsol=%d %d ===================== <<<< \n", kk1, numsol[0], numsol[1]);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
	  printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[solset][kk1], Cip[solset], nip[solset]));

	  if (kk1 < numsol[1-solset])
	    {
	      printf("DISCARDED SOLSET [%d]\n", 1-solset);
	      printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[1-solset][kk1]));
	      printf("distanz punto-asse rim=%.15G\n", perpcomp(solarr[1-solset][kk1], Cip[1-solset], nip[1-solset]));
	    }
#ifdef MC_QUART_VERBOSE
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(D2));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[solset][1], Cip[solset][2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("(%.18G)*x^4+(%.18G)*x^3+(%.18G)*x^2+(%.18G)*x+(%.18G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.18G,%.18G,%.18G,%.18G,%.18G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
#endif
	  printf("# %d >>>> =====================  \n", kk1);
	  //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
	  //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#if 0
	  if (coeff[4] < 1E-10) 
	    {
	      for (kk1=0; kk1 < numsol; kk1++)
		printf("sol=%.20G\n", solqua[kk1]);
	      exit(-1);
	    }
#endif
	}
#endif
      sp = dsc*nip[solset];
      if (fabs(sp) < B.L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
template <class ntype>
ntype rimdisk(hardcyl<ntype> A, hardcyl<ntype> B)
{
  ntype sp, perpdist[2];
  pvec3<ntype> DjCi;
  int whichdisk;
  /* choose the disk closest to cyl axis */
  DjCi = B.disk[0] - A.r;
  sp = DjCi*A.n;
  perpdist[0]=abs(DjCi - sp*A.n);
  DjCi = B.disk[1] - A.r;
  sp = DjCi*A.n;
  perpdist[1]=abs(DjCi - sp*A.n);
  whichdisk=(perpdist[0]<perpdist[1])?0:1;
  if (perpdist[whichdisk] <= 0.5*(A.D+B.D) && rimdiskone(A, B, whichdisk) < 0.0)
    return -1;
  /* choose the disk closest to cyl axis */
  DjCi = A.disk[0] - B.r;
  sp = DjCi*B.n;
  perpdist[0]=abs(DjCi - sp*B.n);
  DjCi = A.disk[1] - B.r;
  sp = DjCi*B.n;
  perpdist[1]=abs(DjCi - sp*B.n);
  whichdisk=(perpdist[0]<perpdist[1])?0:1;
  if (perpdist[whichdisk] <= 0.5*(A.D+B.D) && rimdiskone(B, A, whichdisk) < 0.0)
    return -1;
  return 0;
}
template <class ntype>
double overlap(class hardcyl<ntype> A, class hardcyl<ntype> B, pvector<ntype,3> shift)
{
  pvec3<ntype> CiCj, pv;
  ntype sp, perpc, parac;
  A.shift << 0,0,0;
  bool ret;
  B.shift = shift;
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
   return ret; 
  if ((ret=rimrim(A,B))!=0)
   return ret; 
  if ((ret=rimdisk(A,B))!=0)
   return ret; 
  return 1.0;
}
