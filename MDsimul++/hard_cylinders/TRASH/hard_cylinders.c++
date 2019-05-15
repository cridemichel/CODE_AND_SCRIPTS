#include<cmath>
#include"hard_cylinders.H"
#define Sqr(x) ((x)*(x))
int test_for_fallback(double *P, double *Cip, double *nip, double D2, double *diff)
{
#ifdef MC_QUART_USE_ANALYTIC
  const double DIST_THR=1E-4;
#else
#ifdef MC_QUART_HYBRID
#ifdef FAST_QUARTIC_SOLVER
  const double DIST_THR=5E-13;// old value 5.0E-13 -- NOTA 22/02/18 con il nuovo quartic solver ci si può spingere a 5E-14 volendo!
#else
  const double DIST_THR=5E-12;
#endif
#else
  const double DIST_THR=5E-12;
#endif
#endif
  double diff1, diff2;
  diff1=abs(perpcomp(P, Cip, nip)-D2); // qui D2 è il diametro del rim
  diff2=abs(sqrt(Sqr(P[1])+Sqr(P[2]))-D2);// qui D2 è il diametro del disco

  *diff=diff1+diff2;
  if (diff1 > DIST_THR*D2 || diff2 > DIST_THR*D2)
    return 1;
  else
    return 0;
#if 0
  const double DIST_THR=5E-8;
  if ((*diff=fabs(perpcomp(P, Cip, nip)-D2)) > DIST_THR*D2)
    return 1;
  else 
    return 0;
#endif
}
template <class ntype> ntype hard_cylinder<ntype>::rimdiskone_hybrid(ntype D, ntype L, pvector<ntype,3>& Ci, pvector<ntype,3> ni, 
                                                                     pvector<ntype,3> Dj, pvector<ntype,3> nj)
{
  int kk1, kk2, numsol[2], nsc, fallback, solset;
#ifdef MC_QUART_VERBOSE
  static long int numfb=0;
#endif
  const double FALLBACK_THR = 1E-4;
  pvector<ntype,3> nip[2], Cip[2];
  pvector<complex<ntype>,4> roots;
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4], solquasort[4], solquad[2], uy[3];
  double dsc[3], dscperp[3], c0, c1, c2, c3, c02, c12, c22, nipp[3], Cipp[3], coeffEr[6], rErpp1sq, rErpp2sq, c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2], diffxy[2][4];
  double norm, Rl[3][3];
  double nip02,nip12,nip22,Cip02,Cip12,Cip22, temp;
  double omnip02, omnip12, omnip22;
  double D2sq, D2, Cip0, Cip1, Cip2, nip0, nip1 , nip2, nip1nip2, nip0nip2, nip0nip1; 
  rpoly<ntype,4> oqs;
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  D2 = D*0.5; 
  D2sq = Sqr(D2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  for (kk1=0; kk1 < 3; kk1++)
    uy[kk1]=Dj[kk1];
  versor_to_R_alt_fb(Ci, ni, Dj, nj, Rl, D, uy, 1); 

  //versor_to_R_alt(Ci, ni, Dj, nj, Rl, D); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[0][kk1] = 0;
      //Aip[kk1] = 0;
      Cip[0][kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[0][kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[0][kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = nip.norm();
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
  //xC=yC=0;
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
      oqs.set_coeff(coeff);
      oqs.find_roots(roots);
      numsol[0]=0;
      for (kk1=0; kk1 < 4; kk1++)
        {
          if (imag(roots[kk1])==0)
            { 
              solqua[numsol[0]]=real(roots[kk1]);
              numsol[0]++;
            }
        }
      //solve_quartic(coeff, &(numsol[0]), solqua);
    }
  // 08/05/2019: SONO ARRIVATO QUI <===============================================================================================
  discard_spurious(solqua, &(numsol[0]));

  /* ora assegno a solec[][] e calcolo x */
  /* use bisection newton-raphson to refine solutions */
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
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
#if defined(MC_QUART_VERBOSE) && 0
      printf("numsol=%d,", numsol[0]);
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  printf("sol=%.16G ", solqua[kk1]);
	}
      printf("\n");
      for (kk1=0; kk1 < numsol[0]; kk1++)
	{
	  temp = c4 + c2*solqua[kk1];
	  printf("temp[%d]=%.16G\n", kk1, temp);
	}
      printf("c2=%.16G c4=%.16G c5=%.16G\n", c2, c4, c5);
      store_bump(iGbl,jGbl);
#endif

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
        {
          solve_quartic(coeff, &(numsol[1]), solqua);
        }
      discard_spurious(solqua, &(numsol[1]));
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
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[solset][kk2];
	}
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
      sp = scalProd(dsc, nip[solset]);
      if (fabs(sp) < L*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}
void versor_to_R_alt(double *Ci, double *ni, double *Dj, double *nj, double R[3][3], double D)
{
  int k, kk1, kk2, kk;
  double u[3], norm, sp, dsc[3]; 
  double normDjCi, DjCi[3], DjCini, Ai[3], AiDjnj, AiDjni, AiDj[3], Tnew[3], VV[3], dscperp[3], dscpara[3], ragg, TnCi[3];
  /* first row vector */
  for (k=0; k < 3; k++)
    R[0][k] = nj[k];

  /* N.B. qui faccio uno step dell'algorito di Ibarra semplificato per
   * determinare l'asse y del riferimenti del disco (l'asse x è l'asse perpendicolare
   * al disco e l'asse z si ottiene con il prodotto vettore dell'asse x e y) */
#if 1
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  normDjCi = calc_norm(DjCi);
  DjCini = scalProd(DjCi,ni);

  for (kk1 = 0; kk1 < 3; kk1++)
    Ai[kk1] = Ci[kk1] + DjCini*ni[kk1];
  for (kk1=0; kk1 < 3; kk1++)
    AiDj[kk1] = Ai[kk1] - Dj[kk1]; 
  AiDjnj = scalProd(AiDj, nj);
  for (kk1=0; kk1 < 3; kk1++)
    {
      VV[kk1] = AiDj[kk1] - AiDjnj*nj[kk1];
    }
  //for (kk1=0; kk1 < 3; kk1++)
    //dscpara[kk1] = dscperp[kk1] - Dj[kk1];
  ragg = calc_norm(VV);

  for(k=0;k<3;k++)
    {
      R[1][k] = VV[k]/ragg;
      //R[1][k] = VV[k];
      //TnCi[k] = Tnew[k]-Ci[k];
    }
#if 0
  ragg = scalProd(TnCi,ni);
  for (k=0;k<3;k++)
    Ai[k] = Ci[k] + ragg*ni[k];
#endif
#else

 for (k=0; k < 3; k++)
    dsc[k] = Ci[k] - Dj[k]; 
  sp = scalProd(dsc, nj);
  for (k=0; k < 3; k++)
    R[1][k] = dsc[k] - sp*nj[k];
#endif  
  //printf("scalProd=%.15G\n", scalProd(R[1],R[0]));
  vectProdVec(R[0], R[1], u);
  for (k=0; k < 3 ; k++)
    R[2][k] = u[k];
}
double test_overlap_parall_cyll(long double *Ci, long double *ni, long double *Dj, long double *nj, long double Li, 
				long double Diami, long double Diamj)
{
  int kk;
  long double DjCi[3], DjCini, Ui[3], DjUi[3], normDjUi;
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  //normDjCi = calc_norm(DjCi);
  DjCini = scalProdl(DjCi,ni);

  for (kk=0; kk < 3; kk++)
    {
      Ui[kk] = Ci[kk] + DjCini*ni[kk];
      DjUi[kk] = Dj[kk] - Ui[kk];
    }
  normDjUi = calc_norml(DjUi);

  if (normDjUi <= 0.5*(Diamj+Diami) && fabsl(DjCini) <= Li*0.5)
    return -1.0;
  else
    return 1.0;
}
double test_overlap_parall_cyl(double *Ci, double *ni, double *Dj, double *nj, double Li, double Diami, double Diamj)
{
  int kk;
  double DjCi[3], DjCini, Ui[3], DjUi[3], normDjUi;
  for (kk=0; kk < 3; kk++)
    DjCi[kk] = Dj[kk] - Ci[kk];
  //normDjCi = calc_norm(DjCi);
  DjCini = scalProd(DjCi,ni);

  for (kk=0; kk < 3; kk++)
    {
      Ui[kk] = Ci[kk] + DjCini*ni[kk];
      DjUi[kk] = Dj[kk] - Ui[kk];
    }
  normDjUi = calc_norm(DjUi);

  if (normDjUi <= 0.5*(Diamj+Diami) && fabs(DjCini) <= Li*0.5)
    return -1.0;
  else
    return 1.0;
}

void rotate_axes_on_plane(double RR[3][3])
{
  double Rin[3][3];
  int k1, k2, k3;
  double ox, oy, oz, theta, thetaSq, sinw, cosw;
  double OmegaSq[3][3],Omega[3][3], M[3][3], Ro[3][3];

  ox = RR[0][0];
  oy = RR[0][1];
  oz = RR[0][2];
  for (k1=0; k1 < 3; k1++)
    for (k2=0; k2 < 3; k2++)
      {
	Rin[k1][k2]=RR[k1][k2];
      }
  //NOTA 15/01/2018
  //se uso ranf altero la sequenza casuale e perdo il confronto la simulazione da 50x10^6 già fatta con il metodo di Alberto 
  //theta = (ranf()>0.5?1.:-1.)*M_PI/4.0;
  theta = M_PI/4.0;
  thetaSq=Sqr(theta);
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

  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  M[k1][k2] = -sinw*Omega[k1][k2]+cosw*OmegaSq[k1][k2];
	}
    }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Ro[k1][k2] = Rin[k1][k2];
	for (k3 = 0; k3 < 3; k3++)
	  Ro[k1][k2] += Rin[k1][k3]*M[k3][k2];
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
     RR[k1][k2] = Ro[k1][k2];
}

double rimdiskonediff(double Diami, double Diamj, double Li, double Lj, double Ci[3], double ni[3], double Dj[3], double nj[3], double DjCini)
{
  int kk1, kk2, numsol[2], nsc, fallback, solset;
  const double FALLBACK_THR = 1E-4;
  double tmp, sp, coeff[5], solarr[2][4][3], solec[4][2], solqua[4], solquasort[4], solquad[2];
  double dsc[3], dscperp[3], c0, c1, c2, c3, c02, c12, c22, nipp[3], Cipp[3], coeffEr[6], rErpp1sq, rErpp2sq, c32, c42, c52, c4, c5;  
  double diff[2][4], maxdiff[2], sumdiff[2], diffxy[2][4];
  double Cip[3], nip[3], norm, Rl[3][3];
  double nip02,nip12,nip22,nip03,nip13,nip23,nip04,nip14,nip24,Cip02,Cip12,Cip22, temp;
  //long double c0l, c1l, c2l, c3l, c4l, c5l, templ, solqual;
  //double aErcut, bErcut, nErcutx[3], nErcuty[3], nErcutz[3], rErcut[3], m00, m01, m10, m11, m002, m112, AA, BB, invm10, ev0, ev1, AA0, BB0;
  //double fact,nErcutxp[3], nErcutyp[3], nErcutzp[3], rErcutp[3], aErcut2, bErcut2, nErcutyp12, nErcutyp22, nErcutzp12, nErcutzp22;
  //double ia00, ia01, ia10, ia11, ia002, ia102, ia012, ia112, delta;
  double Di2sq, Dj2sq, Di2, Dj2, Cip0, Cip1, Cip2, nip0, nip1 , nip2; 
/* LAST ATTEMPT */
  /* se asse del rim e asse del disco sono paralleli si deve considerare un caso a parte */
  Dj2 = Diamj*0.5; 
  Di2 = Diami*0.5;
  Dj2sq = Sqr(Dj2);
  Di2sq = Sqr(Di2);
  /* mi metto nel riferimento del disco (p) */
#if 0
  versor_to_R(nj[0], nj[1], nj[2], Rl);
#else
  versor_to_R_alt(Ci, ni, Dj, nj, Rl, Diamj); 
#endif
  for (kk1=0; kk1 < 3; kk1++)
    {
      nip[kk1] = 0;
      //Aip[kk1] = 0;
      Cip[kk1] = 0;
      for (kk2=0; kk2 < 3; kk2++)
	{
	  nip[kk1] += Rl[kk1][kk2]*ni[kk2];
	  Cip[kk1] += Rl[kk1][kk2]*(Ci[kk2]-Dj[kk2]);
	  //Aip[kk1] += Rl[kk1][kk2]*(Ai[kk2]-Dj[j2][kk2]);
	} 
    }
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
  norm = calc_norm(nip);
  nip0 = nip[0]/norm;
  nip1 = nip[1]/norm;
  nip2 = nip[2]/norm;
  Cip0 = Cip[0];
  Cip1 = Cip[1];
  Cip2 = Cip[2];
  nip02=Sqr(nip0);
  nip12=Sqr(nip1);
  nip22=Sqr(nip2);
  nip04=Sqr(nip02);
  nip14=Sqr(nip12);
  nip24=Sqr(nip22);
  nip03=nip02*nip0;
  nip13=nip12*nip1;
  nip23=nip22*nip2;
  Cip02=Sqr(Cip0);
  Cip12=Sqr(Cip1);
  Cip22=Sqr(Cip2);   
#if 1
  coeffEr[0] = 1 - 2*nip12 + nip02*nip12 + nip14 + 
    nip12*nip22;
  coeffEr[1] = 1 - 2*nip22 + nip02*nip22 + 
    nip12*nip22 + nip24;
  coeffEr[2] = -4*nip1*nip2 + 2*nip02*nip1*nip2 + 2*nip13*nip2 + 
    2*nip1*nip23;
  /* sistemare quest con i giusti Dj2sq e Di2sq */
  coeffEr[3] = Cip02 + Cip12 + Cip22 - Di2sq - 
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
#if 0
  {
  /* ora trovo i 6 coefficienti dell'ellisse del rim (c0*x^2 + c1*y^2 + c2*xy + c3 + c4*x + c5*y=0)*/
    double cq[3], x, lam, solq[2], p[3];
    int numsol;
    lam = -Cip0/nip0;
    x = Cip1 + lam*nip1;
    cq[0] = coeffEr[3]+coeffEr[4]*x+coeffEr[0]*x*x;
    cq[1] = coeffEr[5]+coeffEr[2]*x;
    cq[2] = coeffEr[1];
    solve_quadratic(cq, &numsol, solq);
    p[0] = 0.0;
    p[1] = x;
    p[2] = solq[0];
    if (fabs(perpcomp(p, Cip, nip) - D2) > 3E-8)
      {
	printf("coeff quad=%.16G %.16G %.16G\n", cq[2], cq[1], cq[0]);
      	printf("distance punto ellipse axis=%.16G\n", perpcomp(p, Cip, nip));
	printf("nip.njp=%.15G lam=%.15G\n", nip0, lam);
      }
  }
#endif
  /* applico un'omotetia per ridurre la circonferenza del disco a quella unitaria */	
  coeffEr[0] *= Dj2sq;
  coeffEr[1] *= Dj2sq; 
  coeffEr[2] *= Dj2sq;
  coeffEr[4] *= Dj2;
  coeffEr[5] *= Dj2;
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
  //xC=yC=0;
  coeff[4] = c02 - 2*c0*c1 + c12 + c22;
  coeff[3] = 2*c2*c4 - 2*c0*c5 + 2*c1*c5;
  coeff[2] = -2*c02 + 2*c0*c1 - c22 - 2*c0*c3 + 2*c1*c3 + c42 + c52;
  coeff[1] = -2*c2*c4 + 2*c0*c5 + 2*c3*c5;
  coeff[0] = c02 + 2*c0*c3 + c32 - c42;
  if (coeff[4]==0)
    {
      /* N.B. 08/01/18 forse così è troppo restrittiva e dovrò cambiare in una condizione del tipo 
       * fabs(coeff[4]) < EPSILON, devo fare delle prove per stabilirlo... */
      /* cilindri paralleli */
      return test_overlap_parall_cyl(Ci, ni, Dj, nj, Li, Diami, Diamj);
    }
  else
    solve_quartic(coeff, &(numsol[0]), solqua);
#if 0
  if (numsol==1)
    {
      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
      printf("sol=%.15G\n", solqua[0]);
      printf("BOH\n");
    }
#endif
  discard_spurious(solqua, &(numsol[0]));

  //solve_fourth_deg(coeff, &numsol, solqua);
  /* ora assegno a solec[][] e calcolo x */
#if 0
  if (numsol > 1)
    {
      printf("PRIMA solqua=%.15G %.15G\n", solqua[0], solqua[1]);
      qsort(solqua, numsol, sizeof(double), compare_func);
      //printf("numsol=%d\n", numsol);
      //printf("DOPO solqua=%.15G %.15G\n", solqua[0], solqua[1]);
    }
#endif
  /* use bisection newton-raphson to refine solutions */
#if 0
  if (numsol > 2)
    {
      printf("PRIMA solqua(sorted)= ");
      for (kk1=0; kk1 < numsol; kk1++)
	printf(" %.15G ", solqua[kk1]);
      printf("\n");
    }
#endif
#if 0
  for (kk1=0; kk1 < numsol; kk1++)
    {
      double xg;

      if (kk1==0)
	x1b = -1.1; /* le soluzioni devono essere tra -1 e 1 */
      else
	x1b = (solqua[kk1-1]+solqua[kk1])*0.5;
      if (kk1==numsol-1)
	x2b = 1.1;
      else 
	x2b = (solqua[kk1+1]+solqua[kk1])*0.5;
      xg=solqua[kk1];
#if 0
      if ((kk1 == 0 && xg < -1)
	  ||(kk1==numsol-1 && xg > 1))
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 0);
      else
#endif
	solqua[kk1]=rtsafe(coeff, xg, x1b, x2b, 1E-12, 1);
    }
#endif
#if 0
  printf("DOPO solqua(sorted)= ");
  for (kk1=0; kk1 < numsol; kk1++)
    printf(" %.15G ", solqua[kk1]);
  printf("\n");
#endif
  //if (numsol > 0)
  //printf("numsol=%d\n", numsol);
  fallback = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
#if 0
      c0l = c0;
      c1l = c1;
      c2l = c2;
      c3l = c3;
      c4l = c4;
      c5l = c5;
      solqual=solqua[kk1];
      templ = c4l + c2l*solqual;
      solec[kk1][0]=((double)((-c0l - c3l - c5l*solqual + (c0l - c1l)*Sqr(solqual))/templ));
      solec[kk1][1] = solqua[kk1];
#else
      temp = c4 + c2*solqua[kk1];
      solec[kk1][0] = (-c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]))/temp;
      solec[kk1][1] = solqua[kk1];
      //printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
#if 0
      if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	{
  	  printf("solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	  printf("coeff=%.15G %.15G %.15G %.15G %.15G %.15G\n", c0, c1, c2, c3, c4, c5);
	  printf("numeratore=%.16G\n", -c0 - c3 - c5*solqua[kk1] + (c0 - c1)*Sqr(solqua[kk1]));
	}
#endif
#endif
      /* NOTA: siccome le solzuioni sono tali che |x| < 1 e |y| < 1 se temp è molto minore di 1 vuole dire 
       * anche il denominatore lo è quindi sto dividendo due numeri piccoli con conseguenti errori numerici 
       * per cui meglio se risolvo la quartica in x. */
#if 0
      if (test_solution_xy(solec[kk1], &(diffxy[0][kk1])))
	  fallback=1;
#endif
      if (temp==0.0) 
	{
	  fallback=1;
	}
#if 0
      if (fabs(temp) < FALLBACK_THR)
	{
	  fallback = 1;
	  break;
	}
#endif
    }
#if 0
      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
	     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
	     coeff[1]*solqua[kk1]+coeff[0]);
      //printf("semiaxes=%f %f %f %f\n", aEd, bEd, aEr, bEr);
      //printf("ellips(sol)=%.15G\n", Sqr(solec[kk1][0]/a)+Sqr(solec[kk1][1]/b)-1.0);
#endif
  /* ora trovo i 5 coefficienti della quartica c4*x^4+c3*x^3....*/
#if 0
  if (fallback)
    {
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      solve_quartic(coeff, &numsol, solqua);
      for (kk1=0; kk1 < numsol; kk1++)
	{
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	}
    }
#endif
  sumdiff[0] = maxdiff[0] = 0;
  for (kk1=0; kk1 < numsol[0]; kk1++)
    {
      /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
       * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
      solarr[0][kk1][0] = 0.0;
      solarr[0][kk1][1] = Dj2*solec[kk1][0];
      solarr[0][kk1][2] = Dj2*solec[kk1][1];
#if 1
      if (test_for_fallbackdiff(solarr[0][kk1], Cip, nip, Di2, Dj2, &(diff[0][kk1])))
	{
	  fallback=1;
#if 0
	  if (numsol==4)
	    {
	      printf("%d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		     coeff[1]*solqua[kk1]+coeff[0]);
	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[0][kk1]);
	      printf(">>>> =====================\n");
	    }
#endif
	}
      sumdiff[0] += diff[0][kk1];
      if (diff[0][kk1] > maxdiff[0] || kk1==0)
	maxdiff[0] = diff[0][kk1];  
#endif
    }
  if (tinyimagGBL)
    {
      //printf("BOHHHH\n");
      fallback=2;// 2 vuol dire che solset=0 non ha soluzioni reali quindi se ci sono soluzioni usa il fallback e basta
    }
  solset=0;
#if 1
  if (fallback)
    {
      //printf("falling back\n");
      coeff[4] = c02 - 2*c0*c1 + c12 + c22;
      coeff[3] = 2*c0*c4 - 2*c1*c4 + 2*c2*c5;
      coeff[2] = 2*c0*c1 - 2*c12 - c22 + 2*c0*c3 - 2*c1*c3 + c42 + c52;
      coeff[1] = 2*c1*c4 + 2*c3*c4 - 2*c2*c5;
      coeff[0] = c12 + 2*c1*c3 + c32 - c52;
      if (coeff[4]==0)
	{
	  /* cilindri paralleli */
	  return test_overlap_parall_cyl(Ci, ni, Dj, nj, Li, Diami, Diamj);
	}
      else
	solve_quartic(coeff, &(numsol[1]), solqua);
      discard_spurious(solqua, &(numsol[1]));

      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
#if 0
    	  c0l = c0;
	  c1l = c1;
	  c2l = c2;
	  c3l = c3;
	  c4l = c4;
	  c5l = c5;
	  solqual=solqua[kk1];
	  templ = c5l + c2l*solqual;
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = ((double)((-c1l - c3l - c4l*solqual + (c1l - c0l)*Sqr(solqual))/templ)); 
#else
	  temp = c5 + c2*solqua[kk1];
	  solec[kk1][0] = solqua[kk1];
      	  solec[kk1][1] = (-c1 - c3 - c4*solqua[kk1] + (c1 - c0)*Sqr(solqua[kk1]))/temp; 
	  //printf("fallback:");
	  //test_solution_xy(solec[kk1], &(diffxy[1][kk1]));
#if 0
	  if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
	    {
	      printf("[fallback] solec[%d]=%.16G %.16G temp=%.15G\n", kk1, solec[kk1][0], solec[kk1][1], temp);
	    }
#endif
#endif
	}
      sumdiff[1] = maxdiff[1]=0;
      for (kk1=0; kk1 < numsol[1]; kk1++)
	{
	  /* rimoltiplico le coordinate per D2 per riportarmi alla circonferenza di raggio D2 
	   * (ossia faccio l'omotetia inversa rispetto a quella precedente) */	
	  solarr[1][kk1][0] = 0.0;
	  solarr[1][kk1][1] = Dj2*solec[kk1][0];
	  solarr[1][kk1][2] = Dj2*solec[kk1][1];
#if 0
	  printf("[fallback] solarr[%d]=%.16G %.16G\n", kk1, solarr[0][kk1][1], solarr[0][kk1][2]);
	  printf("[fbprevsol]solarr[%d]=%.16G %.16G\n", kk1, solarr[1][kk1][1], solarr[1][kk1][2]);
#endif
#if 1
	  test_for_fallbackdiff(solarr[1][kk1], Cip, nip, Di2, Dj2, &(diff[1][kk1]));
	  sumdiff[1] += diff[1][kk1];
	  if (diff[1][kk1] > maxdiff[1] || kk1==0)
	    maxdiff[1] = diff[1][kk1];  

#endif
#if 0
  	  if (numsol==4)
  	    {
  	      printf("FALLBACK %d [solset=0] numsol=%d ===================== <<<< \n", kk1, numsol);
  	      printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
  	      printf("ni.nj=%.15G\n", scalProd(ni,nj));
  	      printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
  	      printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
  	      printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
  		     coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
  		     coeff[1]*solqua[kk1]+coeff[0]);
  	      printf("temp=%.15G\n", temp);
	      printf("diff=%.16G\n", diff[1][kk1]);
	      printf(">>>> =====================\n");
	  }
#endif
	
	}
#if 1
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
#if 0
	  if (fallback==3 && solset != 1)
	    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
	}
#endif
    }
#endif
#if 0
  if (fallback && numsol==4)
    printf("CHOSEN SOLSET IS N. %d\n", solset);
#endif
#if 0
  construct_inner_points(solarr, Ci, ni, Dj, nj, D);
#endif
  for (kk1=0; kk1 < numsol[solset]; kk1++)
    {
#if 0
      printf("solarr[%d]=(%f,%f,%f)\n", kk1, solarr[kk1][0],solarr[kk1][1],solarr[kk1][2]);
      printf("norm solarr=%.15G\n", calc_norm(solarr[kk1]));
#endif
      for (kk2=0; kk2 < 3; kk2++)
	{
	  dsc[kk2] = solarr[solset][kk1][kk2] - Cip[kk2];
	}
      //printf("dist centro-punto=%.15G\n", calc_distance(Cjpp,solarr[kk1]));

#if 0
      if (calc_normsq(solarr[kk1])-Sqr(D2) > NEWT_THR)
	{
	  newt2Dquartic(coeffEr, solarr[kk1], D2);
	}
#endif
#if 1
      //if (fabs(perpcomp(solarr[kk1], Cip, nip)-D2) > 1E-11)
      if (test_for_fallbackdiff(solarr[solset][kk1], Cip, nip, Di2, Dj2, &tmp)) 
	{
	  printf("# %d ===================== <<<< \n", kk1);
	  printf("distanza punto-centro disk: %.15G\n", calc_norm(solarr[solset][kk1]));
#if 1
	  printf("distanza punto-centro disksq: %.15G D2^2=%.15G\n", calc_norm(solarr[solset][kk1]), Sqr(Dj2));
	  printf("BOH2BOH2 perpcom=%.15G\n", perpcomp(solarr[solset][kk1], Cip, nip));
	  printf("Cip1=%15G Cip2=%.15G\n", Cip[1], Cip[2]);
	  printf("numsol=%d fallback=%d\n", numsol[solset], fallback);
	  print_vec("ni=",ni);
	  print_vec("nj=",nj);
	  printf("c02=%.15G c0=%.15G c1=%.15G c12=%.15G c22=%.15G\n", c02, c0, c1, c12, c22);
	  printf("c4=%.15G c5=%.15G\n", c4, c5);
	  printf("solec[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("coeffEr=%.16G %.16G %.16G %.16G %.16G %.16G\n", coeffEr[0], coeffEr[1], coeffEr[2], coeffEr[3], coeffEr[4],
		 coeffEr[5]);
#endif
	  //solve_quadratic(coeff, &numsol2, solquad);
	  //if (numsol2> 0)
	  //printf("solqua=%.15G %.15G\n", solquad[0], solquad[1]); 
	  printf("solqua[%d]=%.15G\n", kk1, solqua[kk1]);
	  printf("ni.nj=%.15G\n", scalProd(ni,nj));
	  printf("(%.15G)*x^4+(%.15G)*x^3+(%.15G)*x^2+(%.15G)*x+(%.15G)\n", coeff[4], coeff[3], coeff[2], coeff[1], coeff[0]);
	  printf("{%.15G,%.15G,%.15G,%.15G,%.15G}\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
	  printf("quart(sol)=%.15G\n", coeff[4]*Sqr(solqua[kk1])*Sqr(solqua[kk1])+
		 coeff[3]*Sqr(solqua[kk1])*solqua[kk1] + coeff[2]*Sqr(solqua[kk1])+
		 coeff[1]*solqua[kk1]+coeff[0]);
	  printf("temp=%.15G\n", temp);
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
      sp = scalProd(dsc, nip);
      if (fabs(sp) < Li*0.5)
	{
	  return -1;
	}
    }
  return 1;  
}

double rimdiskdiff(double *D, double *L, double Ci[3], double ni[3], double Di[2][3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, kk, j2, k2, ignore[2];
  char fileop2[512];
  double DjUini, DjTmp[2][3], DjCi[3], DjUi[3], niTmp[3], njTmp[3], perpdist[2];
  double normDjUi, normDjCi, DjCini, Ui[3], CiTmp[3], CjTmp[3];
  double LiTmp, LjTmp, DiamiTmp, DiamjTmp, Li, Lj, Diami, Diamj;
  Diami=D[0];
  Diamj=D[1];
  Li=L[0];
  Lj=L[1];
  for (j1=0; j1 < 2; j1++)
    {
      if (j1==1)
	{
	  //break;
	  for (kk=0; kk < 3; kk++)
	    {
	      for (k2=0; k2 < 2; k2++)
		DjTmp[k2][kk] = Dj[k2][kk];
	      CiTmp[kk] = Ci[kk];
	      niTmp[kk] = ni[kk];
	      njTmp[kk] = nj[kk];
	      DiamiTmp = Diami;
	      DiamjTmp = Diamj;
	      LiTmp = Li;
	      LjTmp = Lj;
	      /* exhange the two particles */	
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = Di[k2][kk];
	      Ci[kk] = Cj[kk];
	      ni[kk] = nj[kk];
	      nj[kk] = niTmp[kk];
	       Diami = Diamj;
	      Diamj = DiamiTmp;
	      Li = Lj;
	      Lj = LiTmp;

	    }
	}
      
      for (k2=0; k2 < 2; k2++)
	{
  	  perpdist[k2]=perpcomp(Dj[k2], Ci, ni);
	  ignore[k2] = 0;
	}

      if (perpdist[0] < perpdist[1])
	{ 
    	  ignore[1] = 1;
	}
      else
	{
	  ignore[0] = 1;
	}

      for (j2=0; j2 < 2; j2++)
	{
	  if (ignore[j2])
	    continue;
	  for (kk=0; kk < 3; kk++)
	    DjCi[kk] = Dj[j2][kk] - Ci[kk];
	  normDjCi = calc_norm(DjCi);
	  DjCini = scalProd(DjCi,ni);
	  for (kk=0; kk < 3; kk++)
	    {
	      Ui[kk] = Ci[kk] + DjCini*ni[kk];
	      DjUi[kk] = Dj[j2][kk] - Ui[kk];
	    }

	  DjUini = scalProd(DjUi,ni);
	  normDjUi = calc_norm(DjUi);
#if 0
	  if (dostorebump)
	    {
	      printf("normDjUi=%.15G DjUini=%.15G\n", normDjUi, DjUini);
	      printf("Ci=%f %f %f Dj=%f %f %f\n", Ci[0], Ci[1], Ci[2], Dj[0], Dj[1], Dj[2]);
	      printf("DjUi=%.15G %.15G %.15G\n", DjUi[0], DjUi[1], DjUi[2]); 
	      printf("Uj=%.15G %.15G %.15G\n", Ui[0], Ui[1], Ui[2]); 
	      printf("nj=%.15G %.15G %.15G\n", ni[0], ni[1], ni[2]);
	      printf("DjCini= %.15G\n", DjCini);
	    }
#endif 
	  if (normDjUi > 0.5*(Diami+Diamj))
	    continue;

	  /* NOTE: in Ibarra et al. Mol. Phys. 33, 505 (2007) 
	     there is some mess about following conditions:
	     The second and third condition on right column of page 514 
	     should read (D=sigma):
	     |Di-Ui| < D/2  && |(Dj-Ci).ni| > L/2

	     |Dj-Ui| < D/2  && |(Dj-Ci).ni| <= L/2

*/
#ifndef MC_IBARRA_SIMPLER
	  /* se sono quasi paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) < 1.0E-8)
	    {
	      if (normDjUi <= 0.5*(Diami+Diamj) && fabs(DjCini) <= Li*0.5)
		return -1;
	      else
		continue;
	    }
#else
	  /* se sono quasi paralleli... */
	  if (1.0-fabs(scalProd(ni,nj)) < 3E-16)
	    {
	      if (normDjUi <= 0.5*(Diami+Diamj) && fabs(DjCini) <= Li*0.5)
		return -1;
	      else
		continue;
	    }
#endif
	  if (normDjUi < Diami*0.5 && fabs(DjCini) > Li*0.5)
	    continue;

	  if (normDjUi < Diami*0.5 && fabs(DjCini) <= Li*0.5)
	    {
#ifdef DEBUG_HCMC
	      if (dostorebump)
		printf("A #1 disk-rim NP=%d\n", Oparams.parnum);
#endif	
	      return -1;
	    }

#ifdef MC_IBARRA_SIMPLER
  	  if (rimdiskone_ibarradiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;
#else
#ifdef MC_QUART_LONG_DOUBLE
#ifdef MC_DEBUG_HCALGO
	    {
	      double alg1, alg2; 
	      alg1 = rimdiskone_ibarradiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      alg2 = rimdiskoneldiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
		  //printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
	      if (alg1!=alg2)
		{
		  store_bump(iGbl, jGbl);
		  printf("Discrepancy between i=%d and j=%d!!\n", iGbl, jGbl);
		  printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
		  //saveCorAscii();
		  sprintf(fileop2 ,"coord-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveCoord(fileop2);
		  sprintf(fileop2 ,"coorbakascii-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveBakAscii(fileop2);
		  //exit(-1);
		}
	      if (alg2 < 0)
		return -1;
	    }
#else
	  /* N.B. NON ANCORA IMPLEMENTATA */
	  if (rimdiskonediffl(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;
#endif
#else
#ifdef MC_DEBUG_HCALGO
	    {
	      double alg1, alg2; 
	      alg1 = rimdiskone_ibarradiff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      alg2 = rimdiskonediff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini);
	      //if ((iGbl==469 || iGbl==38) && (jGbl==469 || jGbl==38))
		  //printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
	      if (alg1!=alg2)
		{
		  store_bump(iGbl, jGbl);
		  printf("Discrepancy between i=%d and j=%d!!\n", iGbl, jGbl);
		  printf("IBARRA=%f QUARTIC=%f\n", alg1, alg2);
		  //saveCorAscii();
		  sprintf(fileop2 ,"coord-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveCoord(fileop2);
		  sprintf(fileop2 ,"coorbakascii-%d-%d-s%d.cor", iGbl, jGbl, Oparams.curStep);
		  saveBakAscii(fileop2);
		  //exit(-1);
		}
	      if (alg2 < 0)
		return -1;
	    }
#else
	  if (rimdiskonediff(Diami, Diamj, Li, Lj, Ci, ni, Dj[j2], nj, DjCini) < 0.0)
	    return -1;

#endif
#endif
#endif
	}
      if (j1==1)
	{
	  for (kk=0; kk < 3; kk++)
	    {
	      /* restore particles*/
	      for (k2=0; k2 < 2; k2++)
		Dj[k2][kk] = DjTmp[k2][kk];
	      Ci[kk] = CiTmp[kk];
	      ni[kk] = niTmp[kk];
	      nj[kk] = njTmp[kk];
	    }
	}

    }
  return 0;
}
double diskdiskdiff(double *D, double *L, double Di[2][3], double Ci[3], double ni[3], double Dj[2][3], double Cj[3], double nj[3])
{
  int j1, j2, kk; 
  double sp, normCiCj, CiCj[3], VV[3], Q1;
  double DiN, DjN, niN[3], njN[3], Djni, Djnj, assex[3], Dini, Pi[3], N[3], Pj[3], normN;
  double normNSq, PiPj[3], normPiPj, Q2, PjDj[3], normPiDi, normPjDj, PiDi[3];
  /* case A.1 (see Appendix of Mol. Sim. 33 505-515 (2007) */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
  if (ni[0]==nj[0] && ni[1]==nj[1] && ni[2]==nj[2])
    {
      /* special case of collinear cylinders (parallel disks) */
      normCiCj = calc_norm(CiCj);
      for (kk=0; kk < 3; kk++)
	VV[kk] = CiCj[kk]/normCiCj;

      if (scalProd(VV,ni)==1.0)
	{
	  if (normCiCj <= 0.5*(L[0]+L[1]))
	    return -1;
	  else
	    return 1;
	}

      /* parallel disks */
      for (j1=0; j1 < 2; j1++)
	for (j2=j1; j2 < 2; j2++)
	  {
	    sp=0.0;
	    for (kk=0; kk < 3; kk++)
	      {
		VV[kk] = Di[j1][kk]-Dj[j2][kk];
		sp += ni[kk]*VV[kk];
	      }
	    if (sp == 0 && calc_norm(VV) < 0.5*(D[0]+D[1]))
	      {
		return -1;
	      }
	  }
    }
  else 
    {
      /* loop over all disk pairs (they are 4) */
      vectProdVec(ni, nj, N);
      vectProdVec(ni,N,niN);
      vectProdVec(nj,N,njN);
      normN=calc_norm(N);
      normNSq=Sqr(normN);
      for (j1=0; j1 < 2; j1++)
	for (j2=0; j2 < 2; j2++)
	  {
	    DiN = scalProd(Di[j1],N);
	    DjN = scalProd(Dj[j2],N);
	    Dini = scalProd(Di[j1],ni);
	    Djnj = scalProd(Dj[j2],nj);
	    for (kk=0; kk < 3; kk++)
	      { 
		Pi[kk] = (DiN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
		Pj[kk] = (DjN*N[kk] + Dini*njN[kk]-Djnj*niN[kk])/normNSq;
	      }
	    for (kk=0; kk < 3; kk++)
	      {
		PiDi[kk] = Pi[kk] - Di[j1][kk];
		PjDj[kk] = Pj[kk] - Dj[j2][kk];
	      }
	    normPiDi = calc_norm(PiDi);
	    normPjDj = calc_norm(PjDj);
	    if (normPiDi <= 0.5*D[0] && normPjDj <= 0.5*D[1])
	      {
		Q1 = sqrt(Sqr(D[0])/4.0-Sqr(normPiDi));
		Q2 = sqrt(Sqr(D[1])/4.0-Sqr(normPjDj));
		for (kk=0; kk < 3; kk++)
		  {
		    PiPj[kk] = Pi[kk] - Pj[kk];
		  }
		normPiPj = calc_norm(PiPj);
		if (normPiPj <= Q1 + Q2)
		  {
#ifdef DEBUG_HCMC
		    if (dostorebump)
		      printf("disk-disk\n");
#endif
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
double rimrimdiff(double *D, double *L, double Ci[3],double ni[3], double Cj[3], double nj[3])
{
  int kk;
  double ViVj[3], lambdai, lambdaj, ninj;
  double CiCj[3], CiCjni, CiCjnj, detA, Vi[3], Vj[3]; 
  /* case A.3 rim-rim overlap */
  for (kk=0; kk < 3; kk++)
    {
      CiCj[kk] = Ci[kk] - Cj[kk];
    }
  CiCjni = scalProd(CiCj,ni);
  CiCjnj = scalProd(CiCj,nj);
  ninj = scalProd(ni, nj);
  detA = Sqr(ninj)-1;

  /* WARNING: solution given in Ibarra et al. Mol. Sim. 33,505 (2007) is wrong */
  lambdai = ( CiCjni - CiCjnj*ninj)/detA;
  lambdaj = (-CiCjnj + CiCjni*ninj)/detA;

  for (kk=0; kk < 3; kk++)
    {
      Vi[kk] = Ci[kk] + lambdai*ni[kk];   
      Vj[kk] = Cj[kk] + lambdaj*nj[kk];
      ViVj[kk] = Vi[kk] - Vj[kk];
    }
  if (calc_norm(ViVj) < 0.5*(D[0]+D[1]) && fabs(lambdai) < 0.5*L[0] && fabs(lambdaj) < 0.5*L[1])
    {
#ifdef DEBUG_HCMC
      if (dostorebump)
	printf("rim-rim NP=%d\n", Oparams.parnum);
#endif	
//      if (sphov > 0.0)
//	printf("boh\n");
      return -1;
    }
  return 0;
}
template <class ntype> ntype HC::overlap(vector<particle<ntype>>& parA, vector<particle<ntype>>& parB, pvector<ntype,3>& shift)
{
  static int firstcall=1;
  const int MAX_ITERATIONS = 1000000;
  ntype L[2], D[2], ret;
  pvector<ntype,3> Ci, Di[2], Dj[2], ni, nj;

  ni = parA.R.extract_col_vec(0);
  nj = parB.R.extract_col_vec(0);
  Ci = parA.r + shift;
  Cj = parB.r + shift;
  L[0] = parA.L; 
  D[0] = parA.D; 
  L[1] = parB.L;
  D[1] = parB.D;
  
  /* centers of mass of disks */
  Di[0]=Ci+0.5*L[0]*ni;
  Di[1]=Ci-0.5*L[0]*ni;
  Dj[0]=Cj+0.5*L[1]*nj;
  Dj[1]=Cj-0.5*L[1]*nj;
  
  if ((ret=rimrimdiff(D, L, Ci, ni, Cj, nj)) != 0.0)
    return ret;

  if ((ret=diskdiskdiff(D, L, Di, Ci, ni, Dj, Cj, nj)) != 0.0)
    return ret;

  if ((ret=rimdiskdiff(D, L, Ci, ni, Di, Dj, Cj, nj)) != 0.0)
    return ret;
  /* case A.2 overlap of rim and disk */
  /* =================================== >>> Part A <<< ========================= */
 /* =================================== >>> Part B <<< ========================= */
  
  //numcallsHC += 4.0; 
  return 1;
}
