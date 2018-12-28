#include <mps/mps.h>
mps_context *s_mps; 
mps_monomial_poly * p_mps;
int N_mps, digits_mps;
void mpsolve_find_roots(double *c)
{
  int i;
 for (i = 0; i <= n_mps; i++)
   mps_monomial_poly_set_coefficient_d (s_mps, p_mps, i, c[i], 0);
 mps_context_set_input_poly (s_mps, MPS_POLYNOMIAL (p_mps));

 /* Set the output precision to DBL_EPSILON and the default goal
  * to approximate. Try to find all the possible digits representable 
  * in floating point. */
 mps_context_set_output_prec (s_mps, digits_mps);
 mps_context_set_output_goal (s_mps, MPS_OUTPUT_GOAL_APPROXIMATE);
 mps_mpsolve (s_mps);
} 
void mpsolve_init(int N, int digits)
{
  s_mps = mps_context_new ();
  p_mps = mps_monomial_poly_new (s_mps, N);
  N_mps=N;
  digits_mps=digits;
}
void mps_free(void)
{
  mps_context_free (s_mps);
}
