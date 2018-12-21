#include <stdlib.h>
#include <stdio.h>
#include "pmatrix.H"
#include "./rpoly.H"
#include<complex>
#include<list>
#include<string>

bool allreal=false, doswap=false;
using ldbl=long double;
#ifndef CASO
#define CASO 4
#endif
#undef M_PI
#define M_PI 3.1415926535897932384626433832795029L
#define Complex(x,y) (ldbl(x)+1il*ldbl(y))
using numty = double;//mpf_float_500;
void calc_coeff(long double* c, complex<long double> er[]);
void calc_coeff_dep_on_case(long double* c, complex<long double> **r)
{
  int i;
#if CASO==1
  // wilkinson
  //
#define NDEG 10
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==2
  // wilkinson
  //
#define NDEG 15
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==3
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = i+1;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==4
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  er[0] = -2.1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]+0.2L;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==5
#define NDEG 10
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = 1.0L/i;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==6
#define NDEG 15
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = 1.0L/i;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==7
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=1; i < NDEG+1; i++)
    { 
      er[i-1] = 1.0L/i;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==8
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = 1.0L/pow(2,NDEG/2-i);
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==9
#define NDEG 20
  static complex <long double> er[NDEG];
  allreal=true;
  for (i=0; i < NDEG; i++)
    { 
      er[i] = 1.0L/pow(2,NDEG/2-i)-3.0L;
    }
  calc_coeff(c, er);
  *r = er;
#elif CASO==10
#define NDEG 20
#if 0
  static complex <long double> er[NDEG];

  for (i=0; i < NDEG; i++)
    { 
      er[i] = cos(M_PI*(2.0L*(i+1L)-1L)/40.0L);
    }
#endif
  allreal=true;
#if 1
  doswap=true;
  static complex <long double> er[NDEG]=
    {-0.98883082622512854506974288293400861L - 
      0.14904226617617444692935471527721756L*1il, \
        -0.98883082622512854506974288293400861L + 
        0.14904226617617444692935471527721756L*1il, \
        -0.90096886790241912623610231950744505L - 
        0.43388373911755812047576833284835875L*1il, \
        -0.90096886790241912623610231950744505L + 
        0.43388373911755812047576833284835875L*1il, \
        -0.73305187182982632852243148927067191L - 
        0.68017273777091939018735870103374024L*1il, \
        -0.73305187182982632852243148927067191L + 
        0.68017273777091939018735870103374024L*1il, \
        -0.50000000000000000000000000000000000L - 
        0.86602540378443864676372317075293618L*1il, \
        -0.50000000000000000000000000000000000L + 
        0.86602540378443864676372317075293618L*1il, \
        -0.22252093395631440428890256449679476L - 
        0.97492791218182360701813168299393122L*1il, \
        -0.22252093395631440428890256449679476L + 
        0.97492791218182360701813168299393122L*1il, 
      0.07473009358642425429093974573476665L - 
        0.99720379718118014822502987087811927L*1il, 
      0.07473009358642425429093974573476665L + 
        0.99720379718118014822502987087811927L*1il, 
      0.36534102436639501454473799892976880L - 
        0.93087374864420425563779924195127531L*1il, 
      0.36534102436639501454473799892976880L + 
        0.93087374864420425563779924195127531L*1il, 
      0.62348980185873353052500488400423981L - 
        0.78183148246802980870844452667405775L*1il, 
      0.62348980185873353052500488400423981L + 
        0.78183148246802980870844452667405775L*1il, 
      0.82623877431599487194516257377267840L - 
        0.56332005806362202774926153802976051L*1il, 
      0.82623877431599487194516257377267840L + 
        0.56332005806362202774926153802976051L*1il, 
      0.95557280578614073281133405376746667L - 
        0.29475517441090421683077298196019097L*1il, 
      0.95557280578614073281133405376746667L + 
        0.29475517441090421683077298196019097L*1il};  
#endif
  calc_coeff(c, er);
  *r = er;
#elif CASO==11
#define NDEG 20
  static complex <long double> er[NDEG]=
    {-0.98883082622512854506974288293400861L- 
      0.14904226617617444692935471527721756L*1il, 
      -0.98883082622512854506974288293400861L+ 
        0.14904226617617444692935471527721756L*1il, 
      -0.90096886790241912623610231950744505L- 
        0.43388373911755812047576833284835875L*1il, 
      -0.90096886790241912623610231950744505L+ 
        0.43388373911755812047576833284835875L*1il, 
      -0.73305187182982632852243148927067191L- 
        0.68017273777091939018735870103374024L*1il, 
      -0.73305187182982632852243148927067191L+ 
        0.68017273777091939018735870103374024L*1il, 
      -0.50000000000000000000000000000000000L- 
        0.86602540378443864676372317075293618L*1il, 
      -0.50000000000000000000000000000000000L+ 
        0.86602540378443864676372317075293618L*1il, 
      -0.22252093395631440428890256449679476L- 
        0.97492791218182360701813168299393122L*1il, 
      -0.22252093395631440428890256449679476L+ 
        0.97492791218182360701813168299393122L*1il, 
      0.07473009358642425429093974573476665L- 
        0.99720379718118014822502987087811927L*1il, 
      0.07473009358642425429093974573476665L+ 
        0.99720379718118014822502987087811927L*1il, 
      0.36534102436639501454473799892976880L- 
        0.93087374864420425563779924195127531L*1il, 
      0.36534102436639501454473799892976880L+ 
        0.93087374864420425563779924195127531L*1il, 
      0.62348980185873353052500488400423981L- 
        0.78183148246802980870844452667405775L*1il, 
      0.62348980185873353052500488400423981L+ 
        0.78183148246802980870844452667405775L*1il, 
      0.82623877431599487194516257377267840L- 
        0.56332005806362202774926153802976051L*1il, 
      0.82623877431599487194516257377267840L+ 
        0.56332005806362202774926153802976051L*1il, 
      0.95557280578614073281133405376746667L- 
        0.29475517441090421683077298196019097L*1il, 
      0.95557280578614073281133405376746667L+ 
        0.29475517441090421683077298196019097L*1il};
  doswap=true;
  allreal=true;
  calc_coeff(c, er);
  *r = er;
#elif CASO==12
#define NDEG 24
  doswap=true;
  allreal=true;
 
  static complex <long double> er[NDEG];

  er[0]=Complex(-3.52E2L, 0);
  er[1]=Complex(-3.52E2L, 0);
  er[2]=Complex(-2.8371450777E2L, -2.9920517772E2L);
  er[3]=Complex(-2.8371450777E2L,  2.9920517772E2L);
  er[4]=Complex(-2.7867414048E2L,  6.1005469197E2L);
  er[5]=Complex(-2.7867414048E2L, -6.1005469197E2L);
  er[6]=Complex(-2.74892372E2L, 0L);
  er[7]=Complex(-2.014171531E2L, 0L);
  er[8]=Complex(-1.255366582E2L, 0L);
  er[9]=Complex(-9.599999999E1L, 0L);
  er[10]=Complex(-8.8692435121E1L,  5.5009607430E2L);
  er[11]=Complex(-8.869243512E1L, -5.5009607430E2L);
  er[12]=Complex(-1.6000000000E1L, 0L);
  er[13]=Complex( 8.23178509855E1L, 0L);
  er[14]=Complex( 8.8692435121E1L, -5.50096074303E2L);
  er[15]=Complex( 8.8692435121E1L,  5.5009607430E2L);
  er[16]=Complex( 1.9293739373E2L,  1.60865921259E3L);
  er[17]=Complex( 1.929373937E2L, -1.6086592125E3L);
  er[18]=Complex( 2.0141715312E2L, 0L);
  er[19]=Complex( 2.7489237213E2L, 0L);
  er[20]=Complex( 7.52E2L, 0L);
  er[21]=Complex( 7.52E2L, 0L);
  er[22]=Complex( 9.1106065E2L,  1.5722L);
  er[23]=Complex( 9.1106065E2L, -1.5722L);
  static long double cs[NDEG+1];
  cs[0]=-54765291428198020791747503747742749163073958404455022926495744.L;
  cs[1]=-4052135566767965847649766745769409681058667331648450681896960.L;
  cs[2]=-31969984081155943263834965670035075493639295858977076674560.L;
  cs[3]=575060225471570237690073740639182419333523437771848417280.L;
  cs[4]=7337981286595499156409929740830030318565357725459415040.L;
  cs[5]=6611223380089859336490797585290455483968982077145088.L;
  cs[6]=-195514288747757987122118583800597358656801082441728.L;
  cs[7]=-726907419403715013562762609680450059293446635520.L;
  cs[8]=197178719520196724204974332265013056299335680.L;
  cs[9]=5968852409133617129605588058090797893943296.L;
  cs[10]=16576506891508825500182005531742679597056.L;
  cs[11]=23375026506968330494765978581548924928.L;
  cs[12]=2206941937668751746514177591607296.L;
  cs[13]=-75617855277818001758431020580864.L;
  cs[14]=-204797687173976372829472423936.L;
  cs[15]= -143150263927579584306872320.L;
  cs[16]=  20214880144364480233472.L;
  cs[17]=  453786251090072698880.L;
  cs[18]=  1265052493274939392.L;
  cs[19]= -968887355572224.L;
  cs[20]=  1015406084096.L;
  cs[21]= -3949133824.L;
  cs[22]=  3284992.L;
  cs[23]= -1728.L;

  cs[24]=1.0L;

  for (i=0; i < NDEG+1; i++)
    {
      c[i] = cs[i];
    }
  *r = er;

#elif CASO==13


#elif CASO==14


#elif CASO==15
// Noferini
#define NDEG 12
  //roots and coefficients were calculated by Wolfram Mathematica with a precision of 1000 digits
  static complex <long double> er[NDEG]=
    {4.6783738689934637477470956503276e-14,3.593464181785079996420231643529e-10,
   0.71315932408901343233326922076029,-0.37460208174361406896066147763122,
   -0.72465168650045961190528692375604,0.73436066385675994997813804359153,
   0.75188917254981479758025308786281,1.3149964483387182814468080288333,
   0.62018121440825773751625731899171,-0.2348464075172120733811932571832,
   0.58416493058065636942158421397624,0.37071900029579493839310950837148};
  doswap=false;
  allreal=true;

  static long double cs[NDEG+1]={-7.4535809976555373129332791846744e-26,1.5934065119636233867683029604694e-12,
   -0.0044336022632989170089290997184383,0.0113879060225037412685170930364,
   0.087649800127665380958117305226893,-0.34172448886511450921576879764911,
   -0.1028434743793867531666720875538,2.0354824813590990435428325566779,
   -2.7765173990337953174050229883067,-0.74359721627240353964458460555424,
   4.5882338934369912760408592287635,-3.7553705787171229543394756980964,1.};

  for (i=0; i < NDEG+1; i++)
    {
      c[i] = cs[i];
    }
  *r = er;

#elif CASO==16
// Noferini
#define NDEG 35
  //roots and coefficients were calculated by Wolfram Mathematica with a precision of 1000 digits
  static complex <long double> er[NDEG]=
    {-1.2411014081169485749720046374712e-12,7.0414870784977958688104322165598e-11,
      -1.1009740691481262157709498052086,1.5099233137877368012385266380292,
      0.55325068571508822704738370366753,0.69822314886842319800567269068431,
      -0.39358122503696295358605458141417,0.92332581983459989341184681633702,
      0.0059157564026236025538034712028596,1.1993653017593120359390784380673,
      0.13629876113291881637131010537191,0.20488160759439544502968669948569,
      -0.11460746969200482916925947542847,-1.6927217660973512540579128019357,
      0.75019149770224904354365907207291,-0.70115292492435923728000908743125,
      0.33883130859046094724111648060249,0.74559022969330081544795377010652,
      -0.19925630407256960572549429611744,-0.19237419926621068066044345659583,
      0.96331179897280306263229146847686,0.50124728967473494688367329217709,
      -0.55557517369328925640364214928378,0.30102070930459680076348711791703,
      -3.01438299249866633422729745287,0.80925481894941304376993377443945,
      0.3413061066689615606686938162135,-0.27739316894258263576601307428711,
      1.4085416925483222067114209422981,-1.6864235615968405940826585651041,
      1.0577528780235342793204416601796,-0.87506945052129112363103003063229,
      1.5063591474908893183951161464415,-1.280008562142808787985258106032,
      -0.38663280512671670128482646047247};
  doswap=false;
  allreal=true;

  static long double cs[NDEG+1]={
    1.2206094960416340906834578680468e-31,9.6615438846790090923275525883592e-20,
    -1.3967063075503082065104296976722e-9,2.3785023237495565148223572868933e-7,
    -9.5300697944611868773994541074766e-8,-0.000034357952282030902005958581364579,
    0.000063628525666157172917759826411451,0.001766113547232277432011182234,
    -0.0049065161060369240035921703313614,-0.042474493857995891272956598697765,
    0.15787384728811660156539693601506,0.49266441526149102367125428105999,
    -2.5831914502149769268168924115978,-2.0890795971123032167677487371407,
    23.098676136323451754813980831842,-9.2318384904714597893797501179091,
    -114.50543995664166100992737517431,140.07544441105945977820146202565,
    296.19506352344378824080481925454,-636.11140719012898663241156002659,
    -279.05119261041554927295990161395,1458.6728965493940763310836619597,
    -386.12147270753405940471651139033,-1806.65305143170590063751606141,
    1308.5725639089255521500945479174,1134.3695853685003449774250740953,
    -1453.3654194955223302264452691014,-217.56788392523778495861322361917,
    815.62492183226671931294225426009,-128.04614992159652741651777995596,
    -238.89353750332911089274852049973,80.274876406730634107375813405455,
    32.99101215120731845199001445294,-15.775892157461408972046263156223,
    -1.4844382000237576047211077710708,1.};

  for (i=0; i < NDEG+1; i++)
    {
      c[i] = cs[i];
    }
  *r = er;
#elif CASO==30
#define NDEG 20
  allreal=true;
  er[0]=0.1;
  for (i=1; i < NDEG; i++)
    { 
      er[i] = er[i-1]/(10.0L);
    }
  calc_coeff(c, er);
  *r = er;
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
  //for (i=0; i < NDEG; i++)
  //cout << " c[" << i << "]=" << c[i] << setprecision(20) << "\n";
#elif CASO==31
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
#elif CASO==32
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
#elif CASO==33
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
void sort_sol_opt(complex<long double> *csol, complex<long double>* exsol, long double* allrelerr)
{
  int k1, k2, k2min;
  int perm[NDEG];
  long double relerr, relerrmin, relerrmax;
  complex <long double> diff, solt[NDEG];
  bool used_exsol[NDEG];
  for (k1=0; k1 < NDEG; k1++)
    used_exsol[k1]=false;
  for (k1=0; k1 < NDEG; k1++)
    {
      bool ini = true;
      for (k2=0; k2 < NDEG; k2++)
        {
          if (used_exsol[k2]==true)
            continue;
          diff = csol[k1] - exsol[k2];
          relerr = (exsol[k2]==complex<long double>(0.0+0.0*1i))?abs(diff):abs(diff/exsol[k2]);
          if (ini==true || relerr <= relerrmin)
           {
             ini=false;
             k2min=k2;
             relerrmin = relerr;
           } 
        }
      perm[k1] = k2min;
      //cout << "perm[" << k1 << "]=" << k2min << "\n";
      allrelerr[k2min] = relerrmin;
      used_exsol[k2min]=true;
    }

  for (k1=0; k1 < NDEG; k1++)
    solt[k1] = csol[k1];

  for (k1=0; k1 < NDEG; k1++)
    csol[perm[k1]] = solt[k1];
}
numty print_accuracy_at(char *str, complex<long double>* csol, complex<long double> *exsol, long double *allrelerr)
{
  /* we follow FLocke here */
  int k1;
  long double relerrmax;
  for (k1=0; k1 < NDEG; k1++)
    {
      if (k1==0 || allrelerr[k1] > relerrmax)
        {
          relerrmax=allrelerr[k1];
        }
    }
  printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
  return relerrmax;
}

void print_roots(char *str, complex<long double> *er, complex<long double> *cr, long double *allrelerr)
{
  printf("CASE %s\n", str);
  for (auto i=0; i < NDEG; i++)
    {
      cout << setprecision(16) << "root #" << i << " EX: "<< er[i] << " C:" << cr[i];
      cout << setprecision(3) << " [ eps: " << allrelerr[i] << " ]\n"; 
    }
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

void calc_coeff(long double* c, complex<long double> er[NDEG])
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
void print_backward_err(char *str, long double* c, complex<long double> *cr)
{
  /* we follow FLocke here */
  int k1;
  long double relerr, relerrmax, cc[NDEG+1];
  calc_coeff(cc,  cr);
  for (k1=0; k1 < NDEG+1; k1++)
    {
      relerr=abs((c[k1]==0)?(cc[k1] - c[k1]):(cc[k1] - c[k1])/c[k1]); 
      cout << "k1=" << k1 << " relerr" << relerr << " " << c[k1] <<" " << cc[k1] << "\n" ;
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((c[k1]==0)?(cc[k1] - c[k1]):(cc[k1] - c[k1])/c[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
}
int main(int argc, char *argv[])
{
  pvector<complex<numty>,NDEG> roots(NDEG);
  char testo2[256];
  complex<long double> cr[NDEG], *er;
  pvector<numty,NDEG+1> c(NDEG+1);
  int algo, i;
  long double ca[NDEG+1];
  long double allrelerr[NDEG];
  if (argc == 2)
    {
      algo = atoi(argv[1]);
    }
  else
    {
      algo = 1;
    }
  if (algo < 0 || algo > 1)
    {
      printf("algorithm must be between 1 and XX algo=%d\n", algo);
      exit(-1);
    }
  calc_coeff_dep_on_case(ca, &er);
  for (i=0; i < NDEG+1; i++)
    c[i]=ca[i];
  c.show("boh");
  cout << "coeff=" << c[NDEG] << "\n";
  if (algo==0)
    {
      rpoly<numty,NDEG> rp(NDEG);
      rp.set_coeff(c);
      rp.show();
      //rp.zroots(roots,false);
      rp.find_roots(roots);
      //roots.show();  
      sprintf(testo2, "OPS");
    }
  else if (algo==1)
    {
      rpoly<numty,NDEG,true> rphqr;
      rphqr.set_coeff(c);
      rphqr.show("p(x)=");
      rphqr.find_roots(roots);
      sprintf(testo2, "HQR");
      //roots.show();  
    }
  for (i=0; i < NDEG; i++)
    cr[i] = roots[i];
  // sort roots and calculate relative error
#if 0
  if (allreal==true)
    {
      struct rerot { long double re; long double im;};
      std::array<struct rerot,NDEG> rero;
      for (i=0; i < NDEG; i++)
        {
          rero[i].re = er[i].real();
          rero[i].im = er[i].imag();
        }
      std::sort(rero.begin(),rero.end(), [&] (struct rerot a, struct rerot b)-> bool {return a.re < b.re;});
      if (doswap==true)
        {
          for (i=0; i+1 < NDEG; i+=2)
            {
              if (rero[i].im > rero[i+1].im)
                swap(rero[i],rero[i+1]);
            }
        }
      for (i=0; i < NDEG; i++)
        {
          er[i] = rero[i].re + 1il*rero[i].im;
          cout << "EX root= " << er[i] << "\n";
        }
      for (i=0; i < NDEG; i++)
        {
          rero[i].re = cr[i].real();
          rero[i].im = cr[i].imag();
        }
      std::sort(rero.begin(),rero.end(), [&] (struct rerot a, struct rerot b)-> bool {return a.re < b.re;});
      if (doswap==true)
        {
          for (i=0; i+1 < NDEG; i+=2)
            {
              if (rero[i].im > rero[i+1].im)
                swap(rero[i],rero[i+1]);
            }
        }
      for (i=0; i < NDEG; i++)
        {
          cr[i] = rero[i].re + 1il*rero[i].im;
          cout << "CALC root= " << cr[i] << "\n";
        }
    }
  else
    {
      sort_sol_opt(cr, er);
    };
#endif
  sort_sol_opt(cr, er, allrelerr);
  print_roots(testo2, er, cr, allrelerr);
  cout << "Forward relarive error:\n";
  print_accuracy_at(testo2, cr, er, allrelerr);
  //cout << "Backward relarive error:\n";
  //print_backward_err(testo2, ca, cr); 
  return 0;
}
