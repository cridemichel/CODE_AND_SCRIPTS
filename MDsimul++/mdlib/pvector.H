template <int N> class pvector {
    int vsize;
public:
    double v[N];
    inline pvector<N> operator +(const pvector<N>& param)
    {
	int i;
	pvector<N> lv;
	for (i=0; i < N; i++)
	  lv.v[i] = v[i] + param.v[i];
	return lv;
    }
 
    inline pvector<N> operator += (const pvector<N>& param)
      {
       	int i;
	for (i=0; i < N; i++)
	  v[i] += param.v[i];
	return *this;
      }

    inline pvector<N> operator *(const double& param) 
    {
        int i;
	pvector<N> lv;
	for (i=0; i < N; i++)
	  lv.v[i] = v[i]*param;
	return lv;
    }
    inline double& operator[](const int& i)
      {
	return v[i];
      }
    inline double scalprod(const pvector<N>& v1, const pvector<N> &v2)
      {
	int i;
	double sum=0;
	for (i=0; i < N; i++)
	  {
	    sum+=v1.v[i]*v2.v[i]; 
	  }
	return sum;
      }
    inline pvector<N> vectprod(const pvector<N>& v1, const pvector<N> &v2)
      {
	pvector<N> v3;
	v3.v[0] = v1.v[1]*v2.v[2] - v1.v[2]*v2.v[1];
	v3.v[1] = v1.v[2]*v2.v[0] - v1.v[0]*v2.v[2];
	v3.v[2] = v1.v[0]*v2.v[1] - v1.v[1]*v2.v[0];
	return v3;
      }
    void show(void)
      {
	int i;
	for (i=0; i < N; i++)
	  printf("%.15G ", v[i]);
	printf("\n");
      }
    void show(const char* str)
      {
	int i;
	printf("[%s] ", str);
	for (i=0; i < N; i++)
	  printf("%.15G ", v[i]);
	printf("\n");
      }
    int size()
      {
	return N;
      } 
    pvector()
      {
	vsize=N;
      }
   
    ~pvector()
      {
      }
};
typedef pvector<8> pvec;
// Lots of other methods here
pvec operator*(const double& s, const pvec& v)
{
  int i;
  pvec lv;
  for (i=0; i < lv.size(); i++)
    lv.v[i] = v.v[i]*s;
  return lv;
}
