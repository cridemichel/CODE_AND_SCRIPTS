#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
#include<list>
#include<string>
#ifndef CASO
#define CASO 1
#endif
using numty = double;
void calc_coeff(numty* c, complex<long double> er[]);
void calc_coeff_dep_on_case(numty* c, complex<long double> **r)
{
  int i;
#if CASO==1 
#define NDEG 20
  static complex <long double> er[NDEG];
  er[0]=0.1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/(10.0L);
    }
  calc_coeff(c, er);
#if 0
  for (i=0; i < NDEG; i++)
    cout << "[CALC] c[" << i << "]=" << c[i] << setprecision(20) << "\n";
  c[0] = 1.00000000000000053791679056187E-55;
  c[1] = 1.11111111100000055584244679722E-45;
  c[2] = 1.12233445443322162787800313739E-36;
  c[3] = 1.12345790111098789834840761654E-28;
  c[4] = 1.12357014577977569473069186308E-21;
  c[5] = 1.12358025801221010756462539505E-15;
  c[6] = 1.12357014577977572218174173241E-10;
  c[7] = 1.12345790111098777138700591932E-6;
  c[8] = 0.0011223344544332213099796513589;
  c[9] = 0.111111111100000006790544659907;
  c[10] = 1.0;
#endif
  for (i=0; i < NDEG; i++)
    cout << " c[" << i << "]=" << c[i] << setprecision(20) << "\n";

  *r = er;
#elif CASO==2 
  //Polynomials with few very clustered roots.
  //Kameny  
#define NDEG 9
  numty K = 1E50;
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[0] = 9.0;
  c[2] = -6.0*K*K;
  c[4] = K*K*K*K;
  c[9] = K*K; 
  r = NULL;
#elif CASO==3
#define NDEG 10
  for (auto i=0; i <= NDEG; i++)
    c[i]=0.0;
  c[NDEG] = 1.0;
  c[0]=1.0;
  c[1] = -300.0;
  c[2] = 30000.0;
  c[3] = -1E6;
  r=NULL;
#elif CASO==5
#define NDEG 11
  static numty ct ={-1E22, 2E21, -1E20, 0, 0, 0, 0, 0, 0, 0, 1.0};
  for (i=0; i < 11; i++)
   c[i] = ct[i]; 
  r=NULL;
#elif CASO==6
#define NDEG 20
  static complex <long double> er[NDEG];
  er[0]=1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/10.0L;
    }
  calc_coeff(c, er);
  *r = er;
#endif
} 
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
void sort_sol_opt(complex<numty> sol[NDEG], complex<numty> exsol[NDEG])
{
  int k1, k2, kk;
  double v, vmin;
  complex<numty> solt[NDEG];
  int perm[NDEG];
  int perm_min[NDEG];
  for (k1=0; k1 < NDEG; k1++)
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
        for (kk=0; kk < NDEG; kk++)
          perm_min[kk]=perm[kk];
        vmin = v;
      }
  } 
  while (std::next_permutation(perm,perm+NDEG));

  for (k2=0; k2 < NDEG; k2++)
    solt[k2] = sol[k2];

  for (k2=0; k2 < NDEG; k2++)
    sol[k2] = solt[perm_min[k2]];
}
numty print_accuracy_at(char *str, complex<numty> csol[NDEG], complex<numty> exsol[NDEG])
{
  /* we follow FLocke here */
  int k1;
  numty relerr, relerrmax;
  for (k1=0; k1 < NDEG; k1++)
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

void print_roots(char *str, complex<long double> *er)
{
  printf("CASE %s\n", str);
  for (auto i=0; i < NDEG; i++)
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

void calc_coeff(numty* c, complex<long double> er[NDEG])
{
  int i, j;
  std::list<std::list<int>> subsets;// list containing all subsets of a given length
  complex<long double> cc[NDEG+1];
  complex<long double> term, segno;
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
          if (NDEG % 2==1)
            {
              if (i%2==0)
                segno=-1;
              else
                segno=1;
            }
          else
            {
              if (i%2==0)
                segno=1;
              else
                segno=-1;
            }
          cc[i] += term*segno;
        }
      lt.clear();
    }
  for (i=0; i < NDEG+1; i++)
    {
      c[i] = real(cc[i]);
    }
}
int main(int argc, char *argv[])
{
  rpoly<numty,NDEG> rp;
  rpoly<numty,NDEG,true> rphqr;
  pvector<complex<numty>,NDEG> roots;
  complex<long double> cr[NDEG], *er;
  pvector<numty,NDEG+1> c;
  int algo, i;
  numty ca[NDEG+1];
  if (argc == 2)
    {
      algo = atoi(argv[1]);
    }
  else
    {
      algo = 1;
    }
  if (algo < 1 || algo > 24)
    {
      printf("algorithm must be between 1 and XX algo=%d\n", algo);
      exit(-1);
    }
  calc_coeff_dep_on_case(ca, &er);
  char testo[] = "Case CASO";
  if (er!=nullptr)
    print_roots(testo, er); 
  for (i=0; i < NDEG+1; i++)
    c[i]=ca[i];
  rp.set_coeff(c);
  rp.show();
  rp.find_roots(roots);
  roots.show();  
  for (i=0; i < NDEG; i++)
    cr[i] = roots[i];
  // sort roots and calculate relative error
  return 0;
}
