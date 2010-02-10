#include<stdlib.h>
#include<stdio.h>
#include<math.h>
double R, G, B;
double Hue_2_RGB( double v1, double v2, double vH )             //Function Hue_2_RGB
{
   if ( vH < 0.0 ) vH += 1.0;
   if ( vH > 1.0 ) vH -= 1.0;
   if ( ( 6.0 * vH ) < 1.0 ) return ( v1 + ( v2 - v1 ) * 6.0 * vH );
   if ( ( 2.0 * vH ) < 1.0 ) return ( v2 );
   if ( ( 3.0 * vH ) < 2.0 ) return ( v1 + ( v2 - v1 ) * ( ( 2.0 / 3.0 ) - vH ) * 6.0 );
   return v1;
}
void hsl2rgb_net(double H, double S, double L)
{
  double var_1, var_2;
  //fprintf(stderr,"DEFAULT ALGO\n");
  if ( S == 0 )                       //HSL from 0 to 1
    {
      R = L;                      //RGB results from 0 to 255
      G = L;
      B = L;
    }
  else
    {
      if ( L < 0.5 ) 
	var_2 = L * ( 1.0 + S );
      else           
	var_2 = ( L + S ) - ( S * L );

      var_1 = 2 * L - var_2;

      R = Hue_2_RGB( var_1, var_2, H + ( 1.0 / 3.0 ) );
      G = Hue_2_RGB( var_1, var_2, H );
      B = Hue_2_RGB( var_1, var_2, H - ( 1.0 / 3.0 ) );
    }
}
double ConvertHueToRGB(double m1, double m2,double hue)
{
  if (hue < 0.0)
    hue+=1.0;
  if (hue > 1.0)
    hue-=1.0;
  if ((6.0*hue) < 1.0)
    return(m1+6.0*(m2-m1)*hue);
  if ((2.0*hue) < 1.0)
    return(m2);
  if ((3.0*hue) < 2.0)
    return(m1+6.0*(m2-m1)*(2.0/3.0-hue));
  return(m1);
}
void ConvertHSLToRGB(const double hue,const double saturation,
  const double lightness)
{
  double b, g, r, m1, m2; 
  /*
    Convert HSL to RGB colorspace.
  */
  if (saturation == 0)
    {
      R=lightness;
      G=R;
      B=R;
      return;
    }
  if (lightness <= 0.5)
  m2=lightness*(saturation+1.0);
  else
    m2=(lightness+saturation)-(lightness*saturation);
  m1=2.0*lightness-m2;
  R=ConvertHueToRGB(m1,m2,hue+1.0/3.0);
  G=ConvertHueToRGB(m1,m2,hue);
  B=ConvertHueToRGB(m1,m2,hue-1.0/3.0);
}
void hsl2rgb_wikipedia(double H, double S, double L)
{
  /* WARNING: questo algoritmo da risultati non conformi a quello preso da ImageMagick
     quindi o c'è qualcosa che non va nell'algoritmo stesso o nella mia implementazione */
  double C, X, HP, R1, G1, B1, m;
  //fprintf(stderr,"WIKIPEDIA ALGO\n");
  if (L <= 0.5)
    C = 2.0*L*S;
  else
    C = (2.0-2.0*L)*S;
  HP = H*360/60.0; /* H is assumed to be between 0 and 1 here */
  X = C*(1 - abs(((int)HP) % 2 - 1));
  if (HP < 1)
    {
      R1 = C;
      G1 = X;
      B1 = 0.0; 
    }
  else if (HP < 2)
    {
      R1 = X;
      G1 = C;
      B1 = 0.0; 
    }
  else if (HP < 3)
    {
      R1 = 0.0;
      G1 = C;
      B1 = X; 
    }
  else if (HP < 4)
    {
      R1 = 0.0;
      G1 = X;
      B1 = C; 
    }
  else if (HP < 5)
    {
      R1 = X;
      G1 = 0.0;
      B1 = C; 
    }
  else if (HP < 6)
    {
      R1 = C;
      G1 = 0.0;
      B1 = X; 
    }
  else 
    R1=G1=B1=0.0; /* if H undefined */
  m = L - 0.5*C;
  R = R1 + m;
  G = G1 + m;
  B = B1 + m;
}
void ConvertHSBToRGB(const double hue,const double saturation, const double brightness)
{
  double f, h, p, q, t;
  /*
    Convert HSB to RGB colorspace.
    HSB=HSV
  */
  if (saturation == 0.0)
    {
      R=brightness;
      G=R;
      B=R;
      return;
    }
  h=6.0*(hue-floor(hue));
  f=h-floor((double) h);
  p=brightness*(1.0-saturation);
  q=brightness*(1.0-saturation*f);
  t=brightness*(1.0-(saturation*(1.0-f)));
  switch ((int) h)
  {
    case 0:
    default:
    {
      R=brightness;
      G=t;
      B=p;
      break;
    }
    case 1:
    {
      R=q;
      G=brightness;
      B=p;
      break;
    }
    case 2:
    {
      R=p;
      G=brightness;
      B=t;
      break;
    }
    case 3:
    {
      R=p;
      G=q;
      B=brightness;
      break;
    }
    case 4:
    {
      R=t;
      G=p;
      B=brightness;
      break;
    }
    case 5:
    {
      R=brightness;
      G=p;
      B=q;
      break;
    }
  }
}

int main(int argc, char** argv)
{
  double  H,S,L;
  if (argc < 4) 
    {
      printf("SYNTAX: hsl2rgb  <H> <S> <L> <R|G|B|A=all> <algorithm (0=default=ImageMagick|1=wikipedia|2=net|3=HSB to RGB from IM)>\n");
      exit(-1);
    }
  H = atof(argv[1]);
  S = atof(argv[2]);
  L = atof(argv[3]);
  
  if (argc==6)
    switch (atoi(argv[5]))
      { 
      case 0:
      default:
	ConvertHSLToRGB(H,S,L);
	break;
      case 1:
	hsl2rgb_wikipedia(H,S,L);
	break;
      case 2:
	hsl2rgb_net(H,S,L);
	break;
      case 3:
	ConvertHSBToRGB(H, S, L);
      } 
  else
    ConvertHSLToRGB(H,S,L);
  if (!strcmp(argv[4],"A")||!strcmp(argv[4],"All")||!strcmp(argv[4],"all"))
    printf("%f %f %f\n",R, G, B);
  else if (!strcmp(argv[4],"R"))
    printf("%f\n",R);
  else if (!strcmp(argv[4],"G"))
    printf("%f\n",G);
  else if (!strcmp(argv[4],"B"))
    printf("%f\n",B);
  else
    {   
      printf("Unrecognized output\n");
      exit(-1);
    }
  return 0;
}
