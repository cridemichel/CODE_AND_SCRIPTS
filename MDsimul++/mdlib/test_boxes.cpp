#include "./boxes.H"
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}
double calcDistBox(double rcmA[3], double saxA[3], double RA[3][3],double rcmB[3], double saxB[3], double RB[3][3])
{
  double RR, R0, R1, cij[3][3], abscij[3][3], AD[3], R01, DD[3];
  double AA[3][3], BB[3][3], EA[3], EB[3], rA[3], rB[3];
  int k, k1, k2, existsParallelPair = 0;
  /* N.B. Trattandosi di parallelepipedi la loro interesezione si puo' calcolare in 
   * maniera molto efficiente */ 
  //return -1;
  for (k=0; k < 3; k++)
    {
      rA[k] = rcmA[k];
      rB[k] = rcmB[k];
      EA[k] = saxA[k];
      EB[k] = saxB[k];
    }
  for (k1 = 0; k1 < 3; k1++)
    {
      for (k2 = 0; k2 < 3; k2++)
	{
	  AA[k1][k2] = RA[k1][k2];
	  BB[k1][k2] = RB[k1][k2];
	}
    	DD[k1] = rA[k1] - rB[k1];
    }
  /* axis C0+s*A0 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[0][k1] =  scalProd(AA[0], BB[k1]);
      abscij[0][k1] = abs(cij[0][k1]);
      if ( abscij[0][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[0] = scalProd(AA[0],DD);
  RR = abs(AD[0]);
  R1 = EB[0]*abscij[0][0]+EB[1]*abscij[0][1]+EB[2]*abscij[0][2];
  R01 = EA[0] + R1;
  if ( RR > R01 )
    return 1.0; /* non si intersecano */
  /* axis C0+s*A1 */
  for (k1 = 0; k1 < 3; k1++)
    {
      cij[1][k1] = scalProd(AA[1],BB[k1]);
      abscij[1][k1] = abs(cij[1][k1]);
      if ( abscij[1][k1] == 1.0  )
	existsParallelPair = 1;
    }
  AD[1] = scalProd(AA[1],DD);
  RR = abs(AD[1]);
  R1 = EB[0]*abscij[1][0]+EB[1]*abscij[1][1]+EB[2]*abscij[1][2];
  R01 = EA[1] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*A2 */
  for (k1= 0; k1 < 3; k1++)
    {
      cij[2][k1] = scalProd(AA[2], BB[k1]);
      abscij[2][k1] = abs(cij[2][k1]);
      if ( abscij[2][k1] == 1.0 )
	existsParallelPair = 1;
    }
  AD[2] = scalProd(AA[2],DD);
  RR = abs(AD[2]);
  R1 = EB[0]*abscij[2][0]+EB[1]*abscij[2][1]+EB[2]*abscij[2][2];
  R01 = EA[2] + R1;
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*B0 */
  RR = abs(scalProd(BB[0],DD));
  R0 = EA[0]*abscij[0][0]+EA[1]*abscij[1][0]+EA[2]*abscij[2][0];
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B1 */
  RR = abs(scalProd(BB[1],DD));
  R0 = EA[0]*abscij[0][1]+EA[1]*abscij[1][1]+EA[2]*abscij[2][1];
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;
  
  /* axis C0+s*B2 */
  RR = abs(scalProd(BB[2],DD));
  R0 = EA[0]*abscij[0][2]+EA[1]*abscij[1][2]+EA[2]*abscij[2][2];
  R01 = R0 + EB[2];
  if ( RR > R01 )
    return 1.0;

  /* At least one pair of box axes was parallel, therefore the separation is
   * effectively in 2D, i.e. checking the "edge" normals is sufficient for
   * the separation of the boxes. 
   */
  if ( existsParallelPair )
    return -1.0;

  /* axis C0+s*A0xB0 */
  RR = abs(AD[2]*cij[1][0]-AD[1]*cij[2][0]);
  R0 = EA[1]*abscij[2][0] + EA[2]*abscij[1][0];
  R1 = EB[1]*abscij[0][2] + EB[2]*abscij[0][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB1 */
  RR = abs(AD[2]*cij[1][1]-AD[1]*cij[2][1]);
  R0 = EA[1]*abscij[2][1] + EA[2]*abscij[1][1];
  R1 = EB[0]*abscij[0][2] + EB[2]*abscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB2 */
  RR = abs(AD[2]*cij[1][2]-AD[1]*cij[2][2]);
  R0 = EA[1]*abscij[2][2] + EA[2]*abscij[1][2];
  R1 = EB[0]*abscij[0][1] + EB[1]*abscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB0 */
  RR = abs(AD[0]*cij[2][0]-AD[2]*cij[0][0]);
  R0 = EA[0]*abscij[2][0] + EA[2]*abscij[0][0];
  R1 = EB[1]*abscij[1][2] + EB[2]*abscij[1][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB1 */
  RR = abs(AD[0]*cij[2][1]-AD[2]*cij[0][1]);
  R0 = EA[0]*abscij[2][1] + EA[2]*abscij[0][1];
  R1 = EB[0]*abscij[1][2] + EB[2]*abscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB2 */
  RR = abs(AD[0]*cij[2][2]-AD[2]*cij[0][2]);
  R0 = EA[0]*abscij[2][2] + EA[2]*abscij[0][2];
  R1 = EB[0]*abscij[1][1] + EB[1]*abscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB0 */
  RR = abs(AD[1]*cij[0][0]-AD[0]*cij[1][0]);
  R0 = EA[0]*abscij[1][0] + EA[1]*abscij[0][0];
  R1 = EB[1]*abscij[2][2] + EB[2]*abscij[2][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB1 */
  RR = abs(AD[1]*cij[0][1]-AD[0]*cij[1][1]);
  R0 = EA[0]*abscij[1][1] + EA[1]*abscij[0][1];
  R1 = EB[0]*abscij[2][2] + EB[2]*abscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB2 */
  RR = abs(AD[1]*cij[0][2]-AD[0]*cij[1][2]);
  R0 = EA[0]*abscij[1][2] + EA[1]*abscij[0][2];
  R1 = EB[0]*abscij[2][1] + EB[1]*abscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  return -1.0;
}

int main(int argc, char **argv)
{
  double cppov, cov;
  if (argc==0)
    printf("ok argv[0]=%s\n", argv[0]);
  pvector<long double,3> shift;
  box<double> B1, B2; 
  shift << 0,0,0;
  B1.sax[0] = 2.0;
  B1.sax[1] = 1.0;
  B1.sax[1] = 1.0;
  B1.R << 1,0,0,0,1,0,0,0,1;
  B1.r << 0,0,0;
  
  B2.sax[0] = 2.0;
  B2.sax[1] = 1.0;
  B2.sax[1] = 1.0;
  
  double Lbox=2.0*(B2.sax[0]+B1.sax[0]);
  double ov=0.0;
  long int tt;
  srand48(2);
#if 0
  C2.r << 2.4293795260408330705104162916541,-0.95339470844345441946643404662609,-0.37793181441483625349064823240042;
  C2.n << 0.69494497380222397531213118782034,0.11621765917114777744156839389689,-0.70960900437057961021025676018326;
  C2.update_disks();
  
  overlap(C1,C2,shift);
  //overlap_hc_c(C1.r.v, C1.n.v, C1.L, C1.D, C2.r.v, C2.n.v, C2.L, C2.D);
  exit(-1);
#endif
  for (tt=0; tt < atof(argv[1]); tt++)
    {
      B2.random_box(Lbox);
      B2.random_orient();
      //C2.r.show("r=");
      //C2.n.show("n=");
#if 1
      if ((cppov=overlap(B1,B2)) < 0.0)
        {
          ov += 1.0;
        }
#endif
#if defined(DEBUG_ALGO) 
      //cout << "C++ overlap tested\n";

      if (calcDistBox(B1.r.v, B1.sax.v, B1.R.m, B2.r.v, B2.sax.v, B2.R.m) < 0.0)
        ov+=1.0;
#if 0
      if (cov*cppov < 0)
        {
          cout << "we have a problem\n";
          exit(-1);
        }
#endif
      //cout << "C overlap tested\n";
#endif
    }
  cout << setprecision(20) << "excluded volume="  << Lbox*Lbox*Lbox*ov/((double)tt) << "\n";
  //printf("ov=%.15G\n", overlap(C1,C2,shift));
  //cout << "overlap=" << overlap(C1,C2,shift) << "\n"; 
  return 0;
}
