template <int n> class linear-sys
{
  Matrix A[n][n];
  Matrix B[n][n];
  void ludcmp(double **a, int n,  int* indx, double* d, int *ok);
  void lubksb(double **a, int n, int* indx, double *b);
  void wrap_dgesv(double **a, double *x, int n, int *ok);
  void InvMatrix(double **a, double **b, int NB);
public:
  int SolveLineq (double **a, double *x, int n);
}
