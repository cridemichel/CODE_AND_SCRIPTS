#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
#include<list>
#include<string>
#define SET_ROOTS
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
#define NDEG 5
#endif
#endif
using numty = double;
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
void sort_sol_opt(int N, complex<numty> *sol, complex<numty> *exsol)
{
  int k1, k2, kk;
  double v, vmin;
  complex<numty> *solt = new complex<numty>[N];
  int *perm = new int[N];
  int *perm_min = new int[N];
  for (k1=0; k1 < N; k1++)
    {
      perm[k1] = k1;
    }
  do {
    //  std::cout << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';
    cout << perm << "\n"; 
    v = 0;
    for (k2=0; k2 < 5; k2++)
      {
        v += (exsol[k2]==complex<double>(0,0))?abs(sol[perm[k2]]-exsol[k2]):abs((sol[perm[k2]]-exsol[k2])/exsol[k2]);
      }
    if (k1==0 || v < vmin)
      {
        for (kk=0; kk < N; kk++)
          perm_min[kk]=perm[kk];
        vmin = v;
      }
  } 
  while (std::next_permutation(perm,perm+N));

  for (k2=0; k2 < N; k2++)
    solt[k2] = sol[k2];

  for (k2=0; k2 < N; k2++)
    sol[k2] = solt[perm_min[k2]];
  delete [] solt;
  delete [] perm; 
  delete [] perm_min;
}
numty print_accuracy_at(int N, char *str, complex<numty> *csol, complex<numty> *exsol)
{
  /* we follow FLocke here */
  int k1;
  numty relerr, relerrmax;
  for (k1=0; k1 < N; k1++)
    {
      relerr=abs((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16G\n", str, relerrmax);
  return relerrmax;
}

void print_roots(int N, char *str, complex<numty> *er)
{
  printf("CASE %s\n", str);
  for (auto i=0; i < N; i++)
    cout << "root #" << i << " "<< er[i] << setprecision(16) << "\n";
}

void print( list<int> l){
    for(list<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}

void subset(list<list<int>>& L, int arr[], int size, int left, int index, list<int> &l)
{
    if(left==0)
      {
        L.push_back(l);
        //print(L.back());
        return;
      }
    for(int i=index; i<size;i++)
      {
        l.push_back(arr[i]);
        subset(L,arr,size, left-1,i+1,l);
        l.pop_back();
      }
} 
int main()
{
  int i, j; 
  rpoly<numty,NDEG> rp;
  rpoly<numty,NDEG,true> rphqr;
  pvector<numty, NDEG+1> c;
  pvector<complex<numty>, NDEG> roots;
  pvector<complex<long double>, NDEG> er;
  pvector<complex<long double>, NDEG+1> cc;
  complex<long double> term, segno;
#ifdef SET_ROOTS 
  er[0]=1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/10.0L;
    }
  //cout << "ntype size=" << sizeof(numty) << "\n";
  // use Vieta's formulas to calculate polynomial coefficients
  //
  std::list<std::list<int>> subsets;// list containing all subsets of a given length
  // list with all list of subsets
  list<list<int>> lt;   
  // a subset of integers
  list<int> l;
  int array[NDEG];
  for (j=0; j < NDEG; j++)
    array[j] = j;
  cc[NDEG] = 1.0;
  for (i=0; i < NDEG; i++)
    {
      // all subsets of NDEG-i elements are stored in list lt
      subset(lt,array,NDEG,NDEG-i,0,l);
      cc[i] = 0;
      // loop over all subsets
      for(list<list<int>>::iterator it=lt.begin(); it!=lt.end(); ++it)
        {
          term=1.0;
          //print(*it);
          // loop over one subset
          for(list<int>::iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
            {
              term *= er[*it2]; 
              //cout << "er[" << *it2 << "]=" << er[*it2] << "\n";
            }
          if (i%2==0)
            segno=-1;
          else
            segno=1;
          cc[i] += term*segno;
        }
      lt.clear();
   }
  for (i=0; i < NDEG+1; i++)
    {
      c[i] = real(cc[i]);
    }
   //c << -1685011.4849863185081630945205688,162908947.864249706268310546875,3464224691.560069561004638671875,20202439788.892055511474609375,31022318528.748729705810546875,-352127.44537145795766264200210571,1;
#elif defined(KAMENY)
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
#endif
  rp.set_coeff(c);
  rp.show();
  rp.find_roots(roots);
  roots.show();  
  exit(-1);
  return 0;
}
