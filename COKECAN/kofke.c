#include <kofke.h>
void krk4::evaluate_by_sim(int eqstps, int prstps, double yv, double kval)
{
  /* prepare dirs only here simulations are performed from shell */ 
 
}
int main (int argc, char **argv)
{
  krk4 *k1, *k2, *k3, *k4;
  int type;
  /* we use a fourth order runge-kutta */
  double y1[2], y2[2];
  type=atoi(argv[1]);

  k1 = new krk4();
  k2 = new krk4();
  k3 = new krk4();
  k4 = new krk4();
  /* get value of y(x) at current integration step (i.e. y<1,2>[0]) */
  if (type==1)
    {
      y1[0] = atof(argv[2]);
    }
  else if (type==2)
    {
      y1[0] = atof(argv[2]);
      y2[0] = atof(argv[3]);
    }
  if (type==1)
    {
      /* lambda=T: we consider here that temperature is varying 
       in this case we need a scalar y, i.e. we need only y1 */
      k1->evaluate_by_sim(simstr.eq_steps, simstr.pr_steps, yv[0], 0.0);
      k2->evaluate_by_sim(simstr.eq_steps, simstr.pr_steps, yv[0], 0.5*dT*k1.value);
      k3->evaluate_by_sim(simstr.eq_steps, simstr.pr_steps, yv[0], 0.5*dT*k2.value);
      k4->evaluate_by_sim(simstr.eq_steps, simstr.pr_steps, yv[0], dT*k3.value); 
      yv[1] = yv[0]+dT*k4.value;
      T += dT;
      /* store new coexistance T and pressure here */ 
    }
  else if (type==2)
    {
      /* we consider here the generalized clapeyron equation as discussed
	 in Vega et al. J. Phys.: Condens. Matter 20, 153101 (2008) */
    }
}
