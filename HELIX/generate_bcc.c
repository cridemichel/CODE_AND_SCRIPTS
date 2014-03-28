#include <stdlib.h>
#include <stdio.h>
void main(int argc, char** argv)
{
   double Lx,Ly,Lz,len, D, rx0, ry0, rz0, rx, ry, rz;
   double Dx, Dy, Dz;
   int i;
   len=11.01;
   D=1.401;
   Lx=Ly=Lz=27.4;
   rz = rz0 = -Lx/2.0+len/2.0;
   ry = ry0 = -Ly/2.0+D/2.0;
   rx = rx0 = -Lz/2.0+D/2.0;
   Dz = 13.6;//Lx/((int) Lx / len);
   Dy = Ly/((int) Ly / D);
   Dx = Lz/((int) Lz / D);
   fprintf(stderr, "Dx=%f Dy=%f  Dz=%f\n", Dx, Dy, Dz);    
   for (i=0; i < 700; i++)
     {	
	printf("%f %f %f 1 0 0 0 1 0 0 0 1 0\n", rx, ry, rz);     
	ry += Dy;
        if (ry > Ly/2.0)
          {
            rx=rx+Dx;
            ry=ry0;
            if (rx > Lx/2.0)
              {
		rz += Dz;
                rx = rx0;
	        if (rz > Lz/2.0)
                  {
                    fprintf(stderr,"I cannot accomodate all particles! (i=%d)\n", i);
                    exit(-1);
		  }
               }
          }	
     }
  printf("%f %f %f\n", Lx, Ly, Lz);
}
