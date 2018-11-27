#include <stdlib.h>
#include "pmatrix.H"
#include "rpoly.H"
#define TINY 1E-20
#define MD_NBMAX 4
#define ARMA
#define GLM
#define EIGEN
#ifdef ARMA
#include<armadillo>
#ifdef GLM
#include <glm/ext/matrix_double4x4_precision.hpp> 
#endif
using namespace arma;
#endif
#ifdef EIGEN
#include "/usr/local/include/eigen3/Eigen/Dense"
#endif
double ranf()
{
  return drand48();
}
double ranf3()
{
return drand48();
}
#ifndef NMAT
#define NMAT 4
#endif
#ifdef GLM
#if NMAT==2
  glm::dmat2 mmg, mmg2, mmg3, mmg4;
#elif NMAT==3
  glm::dmat3 mmg, mmg2, mmg3, mmg4;
#else
  glm::dmat4 mmg, mmg2, mmg3, mmg4;
#endif
  glm::vec4 vg1, vg2, vg3, vg4;
#endif
#ifdef ARMA
  mat mma(NMAT,NMAT), mma2(NMAT,NMAT), mma3(NMAT,NMAT), mma4(NMAT,NMAT);
  vec va1(NMAT), va2(NMAT), va3(NMAT), va4(NMAT); 
#endif
  pmatrixq<double,NMAT> mm, mm2, mm3, mm4;
  //pmatrixq<double,NMAT> mm, mm2, mm3, mm4;
  pmatrixq<long double, NMAT> mml, mml2;
  pvector<double, NMAT> v1, v2, v3, v4;
  pvector<long double, NMAT> vl1, vl2, vl3, vl4;
#ifdef EIGEN
#if NMAT > 120
  Eigen::MatrixXd mme(NMAT,NMAT), mme2(NMAT,NMAT), mme3(NMAT,NMAT), mme4(NMAT,NMAT);
#else
  Eigen::Matrix<double, NMAT, NMAT> mme, mme2, mme3, mme4;
#endif
#endif
int main(int argc, char**argv)
{
  long double x1, y1;
  int i, j, k; 
  long long int cct, ccmax, ccmy, ccarm, ccglm, cceigen, numtent;
  complex<long double> x1c, x2c, x3c, x4c;
  pmatrixq<double,3> m, m2, m3;
  pvector<double,3> v;
  int t;
  double sig=1.0, fact;    
  rpoly<double,4> rp;
  pvector<complex<double>,4> r;
  pvector<double,5> c;
  pvector<complex<long double>,5> cc;
  double A[2][2], B[2][2];
  //pmatrixq<double,NMAT> m22= mm.I();
  //m22.show();
#if 0
  for (i=0; i < NMAT; i++)
    {
      v1[i] = drand48()*1.0-0.5;
      v2[i] = drand48()*1.0-0.5;
      v3[i] = drand48()*1.0-0.5;
      for (j=0; j < NMAT; j++)
	{
	  //mml[i][j]= drand48()*10.0-5;
	  mm[i][j] = drand48()*1.0-0.5;//mml[i][j];
	  mm2[i][j] =drand48()*1.0-0.5;
	  mm3[i][j] =drand48()*1.0-0.5;
	  mm4[i][j] =drand48()*1.0-0.5;

#ifdef GLM
	  mmg[i][j] = mm[i][j];
	  mmg2[i][j] = mm2[i][j];
	  mmg3[i][j] = mm3[i][j];
	  mmg4[i][j] = mm4[i][j];
#endif
	  mma(i,j)=mm[i][j];
	  mma2(i,j)=mm2[i][j];
	  mma3(i,j)=mm3[i][j];
	  mma4(i,j)=mm4[i][j];
	}	  
    }
  v1.show("v1=");
  v1.muladd(mm,v2);
  v2.show("v2=");
  mm.show("mm=");
  v1.show("mm*v2+v1=");

  exit(-1);
  printf("inizio\n");
  //mm.show("prima");
  mm /= (mm2+mm3);
  mm += 3.0*(mm2*3.0+mm3+1.3*mm4)/3.0+mm/3.0+ (mm4+3.0*mm2)/(mm2-mm*mm4+mm*2);
  mm = mm+mm2/2.0+2.0*mm+(mm+mm2)/3.0+3.0*(mm+mm2);
  v1 = v1+v2/3.0+2.0*v1+(v1+v2)/3.0+2.0*(v3+v1);
  mm.show("dopo");
  printf("fine\n");
  exit(-1);
  mm[0][0]=1;
  mm[0][1]=2;
  mm[0][2]=3;
  mm[1][0]=4;
  mm[1][1]=5;
  mm[1][2]=6;
  mm[2][0]=7;
  mm[2][1]=8;
  mm[2][2]=9;
  v1[0] = 1;
  v1[1] = 2;
  v1[2] = 3;
  v2=mm*v1;

  mm.show("mm=");
  v1.show("v1=");
  v2.show("mm2=");

  exit(-1);
#endif
#if 0
  //mm = (mm+mm2+mm3)*(mm2+mm3+mm)*(mm*mm);
  mm = mm+mm2+mm3; 
  // v1 = (mm+mm2+mm3)*v2;
  //mm.muladd(mm,mm2+mm3+mm4);
  exit(-1);
#endif
#if 0
  c[0] = 24.0;
  c[1] = -50.0;
  c[2] = 35.0;
  c[3] = -10.0;
  c[4] = 1.0;
  rp.set_coeff(c);
  //rp.forcehqr();
  rp.find_roots(r);
  rp.show("0==");
  for (i=0; i < 4; i++)
    {
      printf("[i=%d] %.15G+I*%.15G\n", i, r[i].real(), r[i].imag());
    }
  exit(-1);
#endif
  for (i=0; i < NMAT; i++)
    {
      for (j=0; j < NMAT; j++)
	{
	  mm[i][j]= drand48()*10.0-5;
	  mm2[i][j] =drand48()*10.0-5; 
#ifdef GLM
	  mmg[i][j] = mm[i][j];
	  mmg2[i][j] =mm2[i][j];
#endif
	}	  
    }
#if 0
  mm.show("mm=");
  mm2.show("mm2=");
  mm=mm2+mm;
  mm.show("mm+mm2");
  mmg = mmg + mmg2;
#ifdef GLM
  for (int i=0; i < 4; i++)
    {
      for (int j=0; j < 4; j++) 
	printf("%.15LG ", (long double)mmg[i][j]);
      printf("\n");
    }
#endif

  exit(-1);
#endif
#if 0 //////////////////////
  //printf("%Lf\n",mml[i][j]);
  srand48(time(0));
  ccmax=100;
  ccmy=ccarm=ccglm=cceigen=0;
  for (cct=0; cct < ccmax; cct++)
    {
      for (i=0; i < NMAT; i++)
	{
	  va1(i)= v1[i]= vl1[i] = drand48()*10-5.0;
	  va2(i)=v2[i]= vl2[i] = drand48()*10-5.0;
	  va3(i)=v3[i] = vl3[i] = drand48()*10-5.0;
	  va4(i)=v4[i] = vl4[i] = drand48()*10-5.0;
	  for (j=0; j < NMAT; j++)
	    {
	      mml[i][j]= drand48()*10.0-5;
	      mml2[i][j] =drand48()*10.0-5; 
	      mm[i][j] = mml[i][j];
	      mm2[i][j] = mml2[i][j];
	    }	  
	}
      //mm.show("mm=");
      //mm2.show("mm2"); 
#ifdef GLM
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      mmg[i][j] = mml[i][j];
	      mmg2[i][j] = mml2[i][j];
	    }
	}
#endif
#ifdef ARMA
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      mma(i,j) = mml[i][j];
	      mma2(i,j) = mml2[i][j];
	    }
	}
#endif
      //vl1 = mml*vl2 + mml.inv()*vl3;
      //v1 =  mm*v2 + mm.inv()*v3;
#ifdef EIGEN
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      mme(i,j) = mml[i][j];
	      mme2(i,j) = mml2[i][j];
	    }
	}
#endif
 
      mml = (mml-mml2)*mml2.inv()+(mml2+mml)*mml.inv(); 
      mm = (mm-mm2)*mm2.inv()+(mm2+mm)*mm.inv(); 
#ifdef ARMA
      mma = (mma-mma2)*mma2.i()+(mma2+mma)*mma.i(); 
#endif
#ifdef EIGEN
      mme = (mme-mme2)*mme2.inverse()+(mme2+mme)*mme.inverse(); 
#endif
#ifdef GLM
      mmg = (mmg-mmg2)*inverse(mmg2)+(mmg2+mmg)*inverse(mmg); 
#endif
      //mml.show("mml*mml*mml*mml");
     //printf("\nMYRES: ");
      long double res, resmax=1E-13;
#if 0
      //long double mmldet = mml.det();
      //printf("det of mml=%.15LG\n", mml.det());
      //printf("det of mm=%.15G\n", mm.det());
      //res = (mm.det()-mmldet)/mmldet;
      //if (fabs(res) > resmax)
	//printf("detres=%.15LG", res); 
      for (i=0; i < NMAT; i++)
	{
	  res = (((long double)v1[i])-vl1[i])/vl1[i]; 
	  if (fabs(res) > resmax)
	    {
	      ccmy++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
      	}
	
#else
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mm[i][j])-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  ccmy++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }
	  //printf("\n");
	}
#endif
#if 0
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    printf(" %.15G ",mmg[i][j]);
	  printf("\n");
	}
#endif
#ifdef GLM
#if 0
      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
    	  res = (((long double)vg1[i][j])-vl1[i][j])/mml[i][j]; 
	  if (fabs(res) > resmax)
	    {
	      ccglm++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
	}
#else

      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mmg[i][j])-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  ccglm++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }    
	}
#endif
#endif
#ifdef ARMA
      //printf("\nARMARES: ");
#if 1
      //va1 =  mma*va2 + mma.i()*va3;
      //mma = mma+mma2;
      //mma = (mma-mma2)*mma2.i()+(mma2+mma)*mma.i(); 
#else
      mma2 = mma;
      mma *= mma2;
      mma *= mma2;
      mma *= mma2.i();
#endif
#if 0
      for (i=0; i < NMAT; i++)
	{
      	  res = (((long double)va1(i))-vl1[i])/vl1[i]; 
	  if (fabs(res) > resmax)
	    {
    	      ccarm++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
	}
      //printf("det=%.15G\n", det(mma));
      //res = (det(mma)-mmldet)/mmldet;
      // if (fabs(res) > resmax)
	//printf("detres=%.15LG\n", res);
#else
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mma(i,j))-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  ccarm++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }
	  //printf("\n");
	}
      //printf("\n");
#endif
#ifdef EIGEN
#if 0
      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
    	  res = (((long double)vg1[i][j])-vl1[i][j])/mml[i][j]; 
	  if (fabs(res) > resmax)
	    {
	      ccglm++;
	      //printf(" res(%d,%d)=%.20LG ", i, j, res);
	    }
	}
#else

      //printf("\nGLMRES: ");
      for (i=0; i < NMAT; i++)
	{
	  for (j=0; j < NMAT; j++)
	    {
	      res = (((long double)mme(i,j))-mml[i][j])/mml[i][j]; 
	      if (fabs(res) > resmax)
		{
		  cceigen++;
		  //printf(" res(%d,%d)=%.20LG ", i, j, res);
		}
	    }    
	}
#endif
#endif
    }
  printf("ccmy=%lld ccglm=%lld ccarm=%lld cceigen=%lld\n", ccmy, ccglm , ccarm, cceigen);
  //mm.show("mm");
  cout.precision(16);
  cout.setf(ios::fixed);
  //mma.raw_print();
  exit(-1); 
#endif
#endif
  srand48(0);
  for (i=0; i < NMAT; i++)
    {
      v1[i] = drand48()*10.0-5; 
      v2[i] = drand48()*10.0-5;
      v3[i] = drand48()*10.0-5;
      v4[i] = drand48()*10.0-5;
#ifdef ARMA
      va1(i) = v1[i];
      va2(i) = v2[i];
      va3(i) = v3[i];
      va4(i) = v4[i];
#endif
#ifdef GLM
      vg1[i] = v1[i];
      vg2[i] =v2[i];
      vg3[i] =v3[i];
      vg4[i] = v4[i];
#endif
      for (j=0; j < NMAT; j++)
	{
	  //mml[i][j]= drand48()*10.0-5;
	  mm[i][j] = drand48()*10.0-5;//mml[i][j];
#ifdef GLM
	  mmg[i][j] = mm[i][j];
#endif
#ifdef ARMA
	  mma(i,j)=mm[i][j];
#endif
	}	  
    }
  //mm.show("mm=");
  srand48(0);
  if (argc==2)
    numtent=atoi(argv[1]);
  else
    numtent = 200;
  for (t=0; t < numtent; t++)
    {
#if 0
      x1 = sig*(ranf()-0.5);
      y1 = sig*(ranf()-0.5);
      x1c = complex<long double>(x1,0);
      x2c = complex<long double>(y1,0);
      x1 = sig*(ranf()-0.5);
      y1 = sig*(ranf()-0.5);
      x3c = complex<long double>(x1,y1);
      x4c = complex<long double>(x1,-y1);

      cc[4] = 1.0;
      cc[3] = -(x1c+x2c+x3c+x4c);
      cc[2] = (x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
      cc[1] = (-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
      cc[0] = (x1c*x2c*x3c*x4c);
      for (j=0; j < 5; j++)
	c[j] = cc[j].real();
      
      rp.set_coeff(c);
      //rp.forcehqr();
      rp.find_roots(r);
#endif 
  for (i=0; i < NMAT; i++)
    {
      for (j=0; j < NMAT; j++)
	{
	  //mml[i][j]= drand48()*10.0-5;
	  mm[i][j] = drand48()*1.0-0.5;//mml[i][j];
	  mm2[i][j] =drand48()*1.0-0.5;
	  mm3[i][j] =drand48()*1.0-0.5;
	  mm4[i][j] =drand48()*1.0-0.5;

#ifdef GLM
	  mmg[i][j] = mm[j][i];
	  mmg2[i][j] = mm2[j][i];
	  mmg3[i][j] = mm3[j][i];
	  mmg4[i][j] = mm4[j][i];
#endif
#ifdef ARMA
	  mma(i,j)=mm[i][j];
	  mma2(i,j)=mm2[i][j];
	  mma3(i,j)=mm3[i][j];
	  mma4(i,j)=mm4[i][j];
#endif
#ifdef EIGEN
	  mme(i,j)=mm[i][j];
	  mme2(i,j)=mm2[i][j];
	  mme3(i,j)=mm3[i][j];
	  mme4(i,j)=mm4[i][j];
#endif
	}	  
    }
#ifndef SEL
#define SEL 1
#endif
#if SEL==1
      //mm[0][0] += sin(t)/10000.0;
      // mm2 = mm;
      // mm2 *= mm;
      // mm2 *= mm;
      // mm2 *= mm.inv();
      //mm2 = (mm*mm2)*(mm*mm2);
      for (k=0; k < 100; k++)
	{
	  //v3=(v1+v2)+(v3+v4);
#ifdef TEST_SUM
	  mm += mm+mm2+mm3+mm4;
#elif defined(TEST_MUL)
	  mm += mm*mm2*mm3*mm4;
#elif defined(TEST_INV)
	  mm += mm.inv()+mm2.inv()+mm3.inv()+mm4.inv();
#else
	  mm += mm2*mm3*mm4.inv()+(mm4+mm3+mm2);
#endif
#if 1
	  //mm += mm+mm2+mm3;
	  //mm = (mm+mm2+mm3)*(mm2+mm3+mm)*(mm*mm);
	  //mm = (mm+mm2)+(mm3+mm4);
	  //mm.muladd(mm,mm2+mm3+mm4);
	  //mm = mm+mm4+mm2+mm3;
	  //mm.muladd(mm2, mm3);
	  //mm += mm2;
	  //mm += mm3;
	  //mm.fastadd(mm2,mm3);
#else
	  for (i=0; i < NMAT; i++)
	    for (j=0; j < NMAT; j++)
	      mm[i][j] += (mm[i][j] + mm3[i][j]+mm2[i][j]);//(mm[i][j]+mm4[i][j]);
#endif
	}
#elif SEL==2
      for (k=0; k < 100; k++)
	{
#ifdef TEST_SUM
	  mmg += mmg+mmg2+mmg3+mmg4;
#elif defined(TEST_MUL)
	  mmg += mmg*mmg2*mmg3*mmg4;
#elif defined(TEST_INV)
	  mmg += inverse(mmg)+inverse(mmg2)+inverse(mmg3)+inverse(mmg4);
#else
	  mmg += mmg2*mmg3*inverse(mmg4)+(mmg4+mmg3+mmg2);
#endif
	  //mmg += mmg2*mmg3+inverse(mmg4);
	 //mmg =mmg2*mmg3*mmg4;
	  //mmg += inverse(mmg3)+mmg*mmg3;
	  //mmg += inverse(mmg4)+mmg*mmg4;
	  //mmg = (mmg+mmg4)+(mmg2+mmg3);
	}
#if 1
      //mmg2 += ((mmg+mmg)+(mmg+mmg))/mmg;
      //mmg2 += ((mmg*inverse(mmg))*(mmg*inverse(mmg)));
#else
      mmg2 = mmg;
      mmg2 *= mmg;
      mmg2 *= mmg;
      mmg2 *= inverse(mmg);
#endif
#elif SEL==3
      //mma(0,0) += sin(t)/(10000.0);
      //mma2 = ((mma*mma2)*(mma*mma2));
      for (k=0; k < 100; k++)
	{
	  //va3 = cross(cross(va1,va2),cross(va3,va4));
	  //va1 += mma2*va2+mma3*va3+mma4*va4;
	  //mma += mma*mma2*mma3.i()+mma4.i();
	  //mma = (mma3+mma2)+(mma+mma4);
	  //va3=(va1+va2)*3.0+(va3+va4)*5.0;
#ifdef TEST_SUM
	  mma += mma+mma2+mma3+mma4;
#elif defined(TEST_MUL)
	  mma += mma*mma2*mma3*mma4;
#elif defined(TEST_INV)
	  mma += mma.i()+mma2.i()+mma3.i()+mma4.i();
#else
	  mma += mma2*mma3*mma4.i()+(mma4+mma3+mma2);
#endif
	  //mma += mma2*mma3+mma4.i();  
	}
      //mma2 = ((mma*mma.i())*(mma*mma.i()));
#elif SEL==4
      for (k=0; k < 100; k++)
	{
	  //va3 = cross(cross(va1,va2),cross(va3,va4));
	  //va1 += mma2*va2+mma3*va3+mma4*va4;
	  //mma += mma*mma2*mma3.i()+mma4.i();
	  //mma = (mma3+mma2)+(mma+mma4);
	  //va3=(va1+va2)*3.0+(va3+va4)*5.0;
#ifdef TEST_SUM
	  mme += mme*mme2*mme3*mme4;
#elif defined(TEST_MUL)
	  mme += mme*mme2*mme3*mme4;
#elif defined(TEST_INV)
	  mme += mme.inverse()+mme2.inverse()+mme3.inverse()+mme4.inverse();
#else
	  mme += mme2*mme3*mme4.inverse()+(mme4+mme3+mme2);
#endif
	  //mme += mme2*mme3+mme4.inverse();  
	}

#endif
    }
#ifdef GLM
#if 0
  printf("GLM last\n");
  for (int i=0; i < 3; i++)
    {
      for (int j=0; j < 3; j++) 
	printf("%.15G ", mmg2[i][j]);
      printf("\n");
    }
#endif
#endif
#if 0
#if SEL==1
  mm.show("mm=");
#elif SEL==3
#ifdef ARMA
  printf("lastarma=\n");
  mma.raw_print();
#endif
#endif
  mm=mm2+mm3;
  mm.show("mmfinal=");
  mma = mma2+mma3;
  printf("lastarmafinal=\n");
  mma.raw_print();
  //mm2 = ((mm*mm2)*(mm*mm2));
  //mm2.show("mm2last=");
#endif
  return 0;
}
