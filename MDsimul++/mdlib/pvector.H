#ifndef _PVECTOR_
#define _PVECTOR_
#include <limits>
#include<cstdlib>
#include<iostream>
#include<iomanip>
using namespace std;
template <class ntype, int N> class pvector {
    int vsize;
    int maxdigits;
public:
    ntype v[N];
    inline pvector<ntype,N> operator +(const pvector<ntype,N>& param)
    {
	int i;
	static pvector<ntype,N> lv;
	for (i=0; i < N; i++)
	  lv.v[i] = v[i] + param.v[i];
	return lv;
    }
 
    inline pvector<ntype,N> operator += (const pvector<ntype,N>& param)
      {
       	int i;
	for (i=0; i < N; i++)
	  v[i] += param.v[i];
	return *this;
      }

    inline pvector<ntype,N> operator *(const double& param) 
    {
        int i;
	static pvector<ntype,N> lv;
	for (i=0; i < N; i++)
	  lv.v[i] = v[i]*param;
	return lv;
    }

    friend pvector<ntype,N> operator *(const double& param, const pvector<ntype,N>& v) 
      {
	int i;
	static pvector<ntype,N> lv;
	for (i=0; i < N; i++)
	  lv.v[i] = v.v[i]*param;
	return lv;
      }

    inline ntype& operator[](const int& i)
      {
	return v[i];
      }
    inline double dot(const pvector<ntype,N>& v1)
      {
	int i;
	double sum=0;
	for (i=0; i < N; i++)
	  {
	    sum+=v[i]*v1.v[i]; 
	  }
	return sum;
      }
    inline double norm()
      {
	return (*this).dot(*this);
      }
    // cross product
    inline double cross(const pvector<ntype,N>& v1)
      {
	static pvector<ntype,N> v2;
	v2.v[0] = v[1]*v1.v[2] - v[2]*v1.v[1];
	v2.v[1] = v[2]*v1.v[0] - v[0]*v1.v[2];
	v2.v[2] = v[0]*v1.v[1] - v[1]*v1.v[0];
	return v2;
      }
    //^ = cross product 
    inline pvector<ntype,N> operator ^(const pvector<ntype,N> &vv)
      {
	return (*this).cross(vv);
      }
    //* dot product
    inline double operator *(const pvector<ntype,N> &vv)
      {
	return (*this).dot(vv);
      }
    //friend double scalprod(const pvector<ntype,N>& v1, const pvector<ntype,N>& v2);
    void show(void)
      {
	show(NULL);
      }
    void show(const char* str)
      {
	int i;
	if (str!=NULL)
	  cout << str;
	cout << "{";
	for (i=0; i < N; i++)
	  {
	    cout << setprecision(maxdigits) << v[i];
	    if (i < N-1)
	      cout << ",";
	  }
       	cout << "}\n";
      }
    int size()
      {
	return N;
      } 
    pvector()
      {
	vsize=N;
	maxdigits=numeric_limits<ntype>::digits10;
      }
   
    ~pvector()
      {
      }
};
template <class ntype, int N>
double scalprod(const pvector<ntype,N>& v1, const pvector<ntype,N>& v2)
{
  int i;
  double sum=0;
  for (i=0; i < N; i++)
    {
      sum+=v1.v[i]*v2.v[i]; 
    }
  return sum;
}
/* some predefined vectors useful for simulations */
typedef pvector<double,2> vecd2;
typedef pvector<double,3> vecd3;
typedef pvector<double,4> vecd4;
#endif
