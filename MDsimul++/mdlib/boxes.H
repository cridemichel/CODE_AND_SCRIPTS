#include "./pmatrix.H"


template <class ntype> class box; 
template <class ntype> class rectangle; 

template <class ntype> ntype overlap(box<ntype>& A, box<ntype>& B);
template <class ntype> ntype overlap(rectangle<ntype>& A, rectangle<ntype>& B);

template <class ntype>
class box 
{
public:
  void random_orient(void)
    {
      R.random();
    }
  void random_box(ntype Lbox)
    {
      r.random_box();
      r = Lbox*r;
    }

  // < > means that we are a specialization of overlap is a friend function 
  // it is needed othwewise compiler consider following overlap not related to the overlap
  // function defined below
  friend ntype overlap<> (box<ntype>& A, box<ntype>& B);
  pvector<ntype,3> r, sax, shift;
  pmatrixq<ntype,3> R;
  box()
    {
      r << 0,0,0;
      R << 1,0,0,0,1,0,0,0,1;
      sax << 1,1,1;
      shift << 0,0,0;
    }

};

template <class ntype>
class rectangle 
{
public:
  void random_orient(void)
    {
      R.random();
    }
  void random_box(ntype Lbox)
    {
      r.random_box();
      r = Lbox*r;
    }

  // < > means that we are a specialization of overlap is a friend function 
  // it is needed othwewise compiler consider following overlap not related to the overlap
  // function defined below
  friend ntype overlap<> (rectangle<ntype>& A, rectangle<ntype>& B);
  pvector<ntype,2> r, sax, shift;
  pmatrixq<ntype,2> R;
  rectangle()
    {
      r << 0,0;
      R << 1,0,0,1;
      sax << 1,1;
      shift << 0,0;
    }

};

template <class ntype> ntype overlap(rectangle<ntype>& A, rectangle<ntype>& B)
{
  int k;
  pvector<ntype,2> AA[2], BB[2], DD, cij[2], fabscij[2], AD, EA, EB;
  ntype RR, R0, R1, R01;
  //bool existParallelPair=false;
  for (k=0; k < 2; k++)
    {
      AA[k]=A.R.row(k);  
      BB[k]=B.R.row(k);
    }
  EA = A.sax;
  EB = B.sax;
  DD = (A.r+A.shift) - (B.r+B.shift);
  /* axis C0+s*A0 */
  for (k=0; k < 2; k++) 
    {
      cij[0][k] = AA[0]*BB[k];
      fabscij[0][k] = abs(cij[0][k]);
      //if (fabscij[0][k]==1.0)
        //existParallelPair=true;
    }

  AD[0] = AA[0]*DD;
  RR= abs(AD[0]);
  R1 = EB*fabscij[0]; 
  R01 = EA[0] + R1;
  if (RR > R01)
    return 1.0;
  /* axis C0+s*A1 */
  for (k=0; k < 3; k++) 
    {
      cij[1][k] = AA[1]*BB[k];
      fabscij[1][k] = abs(cij[1][k]);
      //if (fabscij[1][k]==1.0)
        //existParallelPair=true;
    }
  AD[1] = AA[1]*DD;
  RR = abs(AD[1]);
  R1 = EB*fabscij[1];
  R01 = EA[1] + R1;

  if (RR > R01)
    return 1.0;
  /* axis C0+s*B0 */
  RR = abs(BB[0]*DD);
  R0 = EA[0]*fabscij[0][0]+EA[1]*fabscij[1][0];
  // R0 = EA*fabscij.col(0);
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B1 */
  RR = abs(BB[1]*DD);
  R0 = EA[0]*fabscij[0][1]+EA[1]*fabscij[1][1];
  // R0 = EA*fabscij.col(1);
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;

  return -1.0;
} 

template <class ntype>
ntype overlap(box<ntype>& A, box<ntype>& B)
{
  int k;
  pvector<ntype,3> AA[3], BB[3], DD, cij[3], fabscij[3], AD, EA, EB;
  ntype RR, R0, R1, R01;
  bool existParallelPair=false;
  for (k=0; k < 3; k++)
    {
      AA[k]=A.R.row(k);  
      BB[k]=B.R.row(k);
    }
  EA = A.sax;
  EB = B.sax;
  DD = (A.r+A.shift) - (B.r+B.shift);
  /* axis C0+s*A0 */
  for (k=0; k < 3; k++) 
    {
      cij[0][k] = AA[0]*BB[k];
      fabscij[0][k] = abs(cij[0][k]);
      if (fabscij[0][k]==1.0)
        existParallelPair=true;
    }

  AD[0] = AA[0]*DD;
  RR= abs(AD[0]);
  R1 = EB*fabscij[0]; 
  R01 = EA[0] + R1;
  if (RR > R01)
    return 1.0;
  /* axis C0+s*A1 */
  for (k=0; k < 3; k++) 
    {
      cij[1][k] = AA[1]*BB[k];
      fabscij[1][k] = abs(cij[1][k]);
      if (fabscij[1][k]==1.0)
        existParallelPair=true;
    }
  AD[1] = AA[1]*DD;
  RR = abs(AD[1]);
  R1 = EB*fabscij[1];
  R01 = EA[1] + R1;

  if (RR > R01)
    return 1.0;
  /* axis C0+s*A2 */
  for (k=0; k < 3; k++) 
    {
      cij[2][k] = AA[2]*BB[k];
      fabscij[2][k] = abs(cij[2][k]);
      if (fabscij[2][k]==1.0)
        existParallelPair=true;
    }
  AD[2] = AA[2]*DD;
  RR = abs(AD[2]);
  R1 = EB*fabscij[2];
  R01 = EA[2] + R1;
  if (RR > R01)
    return 1.0;
  /* axis C0+s*B0 */
  RR = abs(BB[0]*DD);
  R0 = EA[0]*fabscij[0][0]+EA[1]*fabscij[1][0]+EA[2]*fabscij[2][0];
  // R0 = EA*fabscij.col(0);
  R01 = R0 + EB[0];
  if ( RR > R01 )
    return 1.0;
  /* axis C0+s*B1 */
  RR = abs(BB[1]*DD);
  R0 = EA[0]*fabscij[0][1]+EA[1]*fabscij[1][1]+EA[2]*fabscij[2][1];
  // R0 = EA*fabscij.col(1);
  R01 = R0 + EB[1];
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*B2 */
  RR = abs(BB[2]*DD);
  R0 = EA[0]*fabscij[0][2]+EA[1]*fabscij[1][2]+EA[2]*fabscij[2][2];
  // R0 = EA*fabscij.col(2);
  R01 = R0 + EB[2];
  if ( RR > R01 )
    return 1.0;

  if (existParallelPair)
    return -1.0;

  /* axis C0+s*A0xB0 */
  RR = abs(AD[2]*cij[1][0]-AD[1]*cij[2][0]);
  R0 = EA[1]*fabscij[2][0] + EA[2]*fabscij[1][0];
  R1 = EB[1]*fabscij[0][2] + EB[2]*fabscij[0][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB1 */
  RR = abs(AD[2]*cij[1][1]-AD[1]*cij[2][1]);
  R0 = EA[1]*fabscij[2][1] + EA[2]*fabscij[1][1];
  R1 = EB[0]*fabscij[0][2] + EB[2]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A0xB2 */
  RR = abs(AD[2]*cij[1][2]-AD[1]*cij[2][2]);
  R0 = EA[1]*fabscij[2][2] + EA[2]*fabscij[1][2];
  R1 = EB[0]*fabscij[0][1] + EB[1]*fabscij[0][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB0 */
  RR = abs(AD[0]*cij[2][0]-AD[2]*cij[0][0]);
  R0 = EA[0]*fabscij[2][0] + EA[2]*fabscij[0][0];
  R1 = EB[1]*fabscij[1][2] + EB[2]*fabscij[1][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB1 */
  RR = abs(AD[0]*cij[2][1]-AD[2]*cij[0][1]);
  R0 = EA[0]*fabscij[2][1] + EA[2]*fabscij[0][1];
  R1 = EB[0]*fabscij[1][2] + EB[2]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A1xB2 */
  RR = abs(AD[0]*cij[2][2]-AD[2]*cij[0][2]);
  R0 = EA[0]*fabscij[2][2] + EA[2]*fabscij[0][2];
  R1 = EB[0]*fabscij[1][1] + EB[1]*fabscij[1][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB0 */
  RR = abs(AD[1]*cij[0][0]-AD[0]*cij[1][0]);
  R0 = EA[0]*fabscij[1][0] + EA[1]*fabscij[0][0];
  R1 = EB[1]*fabscij[2][2] + EB[2]*fabscij[2][1];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB1 */
  RR = abs(AD[1]*cij[0][1]-AD[0]*cij[1][1]);
  R0 = EA[0]*fabscij[1][1] + EA[1]*fabscij[0][1];
  R1 = EB[0]*fabscij[2][2] + EB[2]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  /* axis C0+s*A2xB2 */
  RR = abs(AD[1]*cij[0][2]-AD[0]*cij[1][2]);
  R0 = EA[0]*fabscij[1][2] + EA[1]*fabscij[0][2];
  R1 = EB[0]*fabscij[2][1] + EB[1]*fabscij[2][0];
  R01 = R0 + R1;
  if ( RR > R01 )
    return 1.0;

  return -1.0;

}  
