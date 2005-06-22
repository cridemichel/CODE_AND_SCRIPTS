typedef struct {
  char nomefile[1000];   //
 
  int parnum;       // number of particles
  int parnumA;      // number of particles of type A
  int parnumB;      // number of particles of type B
  
  double time;      //
  
  double mA;     // mass for type A
  double mB;     // mass for type B
  
  double kbT;       // temperature

  double aA;     // x-axis for type A 
  double bA;     // y-axis for type A
  double cA;     // z-axis for type A
  double IA;     // momentum of inertia for type A

  double aB;     // x-axis for type B 
  double bB;     // y-axis for type B
  double cB;     // z-axis for type B
  double IB;     // momentum of inertia for type B

  double box;       // side of the box

  double *rx;
  double *ry;
  double *rz;

  double** R;        // R[i][0..8] = Rxx Rxy Rxz Ryx Ryy Ryz Rzx Rzy Rzz 

  double *vx;
  double *vy;
  double *vz;

  double *wx;
  double *wy;
  double *wz;

} conf_t;
