#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
//#define KAMENY
//#define MIGNOTTE
//#define WILK2
#ifdef KAMENY
#define NDEG 9
#elif defined(MIGNOTTE)
#define NDEG 6
#elif defined(WILKINSON)
#define NDEG 10
#elif defined(WILK2)
#define NDEG 10
#else
#ifndef NDEG
#define NDEG 6
#endif
#endif
// The main function that prints all combinations of  
// size r in arr[] of size n. This function mainly 
// uses combinationUtil() 
#if 0
void printCombination(int arr[], int n, int r) 
{ 
    // A temporary array to store all combination 
    // one by one 
    int data[r]; 
  
    // Print all combination using temprary array 'data[]' 
    combinationUtil(arr, n, r, 0, data, 0); 
} 
#endif
/* arr[]  ---> Input Array 
   n      ---> Size of input array 
   r      ---> Size of a combination to be printed 
   index  ---> Current index in data[] 
   data[] ---> Temporary array to store current combination 
   i      ---> index of current element in arr[]     */
void combinationUtil(int arr[], int n, int r, int index, 
                     int data[], int i) 
{ 
    // Current cobination is ready, print it 
    if (index == r) { 
        for (int j = 0; j < r; j++) 
            printf("%d ", data[j]); 
        printf("\n"); 
        return; 
    } 
  
    // When no more elements are there to put in data[] 
    if (i >= n) 
        return; 
  
    // current is included, put next at next location 
    data[index] = arr[i]; 
    combinationUtil(arr, n, r, index + 1, data, i + 1); 
  
    // current is excluded, replace it with next 
    // (Note that i+1 is passed, but index is not 
    // changed) 
    combinationUtil(arr, n, r, index, data, i + 1); 
} 
  
int main()
{
  using numty=double;
  int i, j, k, k2; 
  rpoly<numty,NDEG> rp;
  rpoly<numty,NDEG,true> rphqr;
  pvector<numty, NDEG+1> c;
  pvector<complex<numty>, NDEG> roots;
  pvector<complex<long double>, NDEG> er;
  pvector<complex<long double>, NDEG+1> cc;
  complex<long double> term, segno;
#if 1
  er[0]=1.0;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = 1.0;//er[i-1]/10.0L;
    }
  //cout << "ntype size=" << sizeof(numty) << "\n";
  cc[NDEG] = 1.0;
  static int subsets[NDEG][NDEG];
  int ns;
  // use Vieta's formulas to calculate polynomial coefficients
  // TO DO!! 
  //
  int inparr[NDEG];

  for (i=0; i < NDEG; i++)
    inparr[i] = i;
  for (i=0; i < NDEG; i++)
    {

      combinationUtil(inparr, NDEG, NDEG-i, data, i);
      cc[i]=0;

      for (j=0; j < ns; j++)
        {
          k2=0;
          term=1.0;
          while (subsets[ns][k2]!=-1)
            {
              term *= subsets[ns][k2];
              k2++;
            }
          if (i%2==0)
            segno=1;
          else
            segno =-1;
          cc[j] += segno*term;
        };
    }
#if 0
  c[5] = (-x1c - x2c - x3c - x4c - x5c - x6c).real();
  c[4] = (x1c*x2c + x1c*x3c + x2c*x3c + x1c*x4c + x2c*x4c + x3c*x4c + x1c*x5c + x2c*x5c + x3c*x5c + x4c*x5c + 
   x1c*x6c + x2c*x6c + x3c*x6c + x4c*x6c + x5c*x6c).real(); 
  c[3] = (-(x1c*x2c*x3c) - x1c*x2c*x4c - x1c*x3c*x4c - x2c*x3c*x4c - x1c*x2c*x5c - x1c*x3c*x5c - x2c*x3c*x5c - 
   x1c*x4c*x5c - x2c*x4c*x5c - x3c*x4c*x5c - x1c*x2c*x6c - x1c*x3c*x6c - x2c*x3c*x6c - x1c*x4c*x6c - 
   x2c*x4c*x6c - x3c*x4c*x6c - x1c*x5c*x6c - x2c*x5c*x6c - x3c*x5c*x6c - x4c*x5c*x6c).real();
  c[2] = (x1c*x2c*x3c*x4c + x1c*x2c*x3c*x5c + x1c*x2c*x4c*x5c + x1c*x3c*x4c*x5c + x2c*x3c*x4c*x5c + 
   x1c*x2c*x3c*x6c + x1c*x2c*x4c*x6c + x1c*x3c*x4c*x6c + x2c*x3c*x4c*x6c + x1c*x2c*x5c*x6c + 
   x1c*x3c*x5c*x6c + x2c*x3c*x5c*x6c + x1c*x4c*x5c*x6c + x2c*x4c*x5c*x6c + x3c*x4c*x5c*x6c).real();
  c[1] = (-(x1c*x2c*x3c*x4c*x5c) - x1c*x2c*x3c*x4c*x6c - x1c*x2c*x3c*x5c*x6c - x1c*x2c*x4c*x5c*x6c - 
   x1c*x3c*x4c*x5c*x6c - x2c*x3c*x4c*x5c*x6c).real(); 
  c[0] = (x1c*x2c*x3c*x4c*x5c*x6c).real();
#endif
#else
  //c << -1685011.48498632,162908947.86425, 464224691.56007,20202439788.8921,31022318528.7487, -352127.445371458,1;
   //c << -1685011.4849863185081630945205688,162908947.864249706268310546875,3464224691.560069561004638671875,20202439788.892055511474609375,31022318528.748729705810546875,-352127.44537145795766264200210571,1;
#ifdef KAMENY
//Polynomials with few very clustered roots.
//Kameny  
  numty K = 1E50;
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[0] = 9.0;
  c[2] = -6.0*K*K;
  c[4] = K*K*K*K;
  c[9] = K*K; 
#elif defined(MIGNOTTE)
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[NDEG] = 1.0;
  c[0]=1.0;
  c[1] = -300.0;
  c[2] = 30000.0;
  c[3] = -1E6;
#elif defined(WILKINSON)
c << 3628800, - 10628640, +12753576, - 8409500, 3416930, -902055, 157773, -18150,1320,-55, 1.0;
#elif defined(WILK2)
c << -1E22, 2E21, -1E20, 0, 0, 0, 0, 0, 0, 0, 1.0;
#else

#endif
#endif
  rp.set_coeff(c);
  rp.show();
  rp.find_roots(roots);
  roots.show();  
  exit(-1);
  return 0;
}
