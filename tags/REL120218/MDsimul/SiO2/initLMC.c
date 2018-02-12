#include <mdsimul.h>
#define SIMUL
/* CONVENTION: 
   indices a,b are used for atoms inside a molecule , while
   indices i,j ares used for molecules 
   NOTE: The box edge length is unity, so every length must be referred to 
         the this quantity.
*/

/* ==============>>> SHARED COUNTERS (DON'T TOUCH THESE)<<< ================ */

#define LOOKUP
#define RADE
#define NUM_SPOST 6

extern int ENDSIM;
extern char msgStrA[MSG_LEN];

/* ============ >>> MOVE PROCEDURE AND MEASURING FUNCTIONS VARS <<< =========
 Here you can put all the variable that you use only in this file, that is 
 in the move function and in the measuring functions, note that the variables 
 to measures have to be put in the 'mdsimdep.h' file (see that) */
extern COORD_TYPE pi, s1t, Vol1t, L, invL, s1p, Elrc, Plrc;   
extern COORD_TYPE Vc, V, W, K, Wxx, Wyy, Wzz,  Wxy, Wyz, Wzx;  

/* neighbour list method variables */
extern COORD_TYPE dispHi;
extern int **nebrTab, nebrNow, *nebrTabLen, nebrTabMax, **nebrTabi, **nebrTaba;

double ***vabLT[NA][NA];
double ***wabLT[NA][NA];

int cubeSize[NA][NA];

int *pmap, *pabs, *pSqr, *pgetb, *pgetj ;
double *wabRadLT[NA][NA], *vabRadLT[NA][NA];
double **DwabRadLT[NA][NA], **DvabRadLT[NA][NA];

extern int *rxS[NA], *ryS[NA], *rzS[NA];
extern int maxdelta[NA][NA];
 
IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                  
      INTEGER *4 NFI,NACC,KMAX,KSQMAX
      PARAMETER(NMAXM=1000,NLISTMAX=250000)
      PARAMETER ( KMAX = 6, KSQMAX = 36 )
C                                                                               
      LOGICAL  KNEBG,KDISK,KTAPE,KTEMP,KPRES,KGVRO,KINPO                        
C                                                                               
      COMMON /CMSTEK/ KNEBG,KDISK,KTAPE,KTEMP,KPRES,KGVRO,KINPO                
      COMMON /CMSTEL/ LNEBG,LDISK,LTAPE,LCRTL,LPROT,LRSTA,LIOST                
      COMMON /CMZAHL/ NFI  ,NACC ,NATU ,NZAHL,MOLS ,NGMAX,NLIST,KMOLS          
      COMMON /CMSCAL/ LDTEM,LDPRE,DTEM ,DPRE ,ATEMP,APRES,ARHO ,AEGES          
      COMMON /CMREAD/ ETEMP,ERHO ,EPRES,EEGES,DELTA,TAUT ,TAUP ,COMP           
      COMMON /CMPOPA/ AVOGA,EPSON,SIGMA,BOLTK,WAITS,DELTR,PI,WAIT    
      COMMON /SHIFT/  EPSSHIFT
      common /lewis/  aene,bene,epsonr,sigmar,rswitch2

C                                                                   
      COMMON /CMVALU/ VIRIA,VIRIRA,VIRIQA,EKINA,EKINRA,EKINQA,ELECA, 
     1                EGESA,EGESRA,EGESQA,EPOTA,EPOTRA,EPOTQA,EKCMA,           
     2                TEMPA,TEMPRA,TEMPQA,PRESA,PRESRA,PRESQA,FREEA,   
     3                RHOA ,RHORA ,RHOQA ,WORKA,WORKRA,WORKQA,TDPQA,      
     4                TROTA,TCMAA      
C                                                                               
      COMMON /CMOUTP/ EGES,EGESP,EGESRP,EGESF,TEMP,TEMPP,TEMPRP,TEMPF,          
     1                EPOT,EPOTP,EPOTRP,EPOTF,VIRI,VIRIP,VIRIRP,VIRIF,          
     2                EKIN,EKINP,EKINRP,EKINF,PRES,PRESP,PRESRP,PRESF,          
     3                RHO ,RHOP ,RHORP ,RHOF ,WORK,WORKP,WORKRP,WORKF,          
     4                TROT,TROTP,TCMA  ,TCMAP,EKCM,EKCMP,ELEC  ,ELECP,          
     5                FREE,FREEP,TDPQ  ,TDPQP,QDIP                              
C                                                                               
      COMMON /CMCUTO/ VOLUME,FACGVR,FACPOT,WWx,wwy,wwz                          
       COMMON /CMFCON/ CONREP,CONATP,CONREF,CONATF,FMOLS,RONEQ ,RTWOQ,          
     1                FACPRE,FACVIR,FACDEN,FACWOR,FACRF,RANG2Q,RANG2,           
     2                PERMIT,RANGE3,RANG3H,RANGEQ,RANGE,RANTQ ,RANT ,           
     3                R2MR  ,R3MR  ,TR2   ,TR3   ,EM   ,                        
     4                FACTMP(2),FACEKI
C                                                                               
      COMMON /CMLIST/ RNLAST(nmaxm,3)
c
      COMMON /CMKORM/ RM(nmaxm,3)
      COMMON /CMRSTA/ RSTART(nmaxm,3)
 
      COMMON /CMPARK/ PARKE,VE(nmaxm,3),RN(nmaxm,3)                        
      COMMON /CMFORC/ F(nmaxm,3),VIRIAL,EPOTEN,ELEPOT                          
      COMMON /CMALFI/ DIP(nmaxm,3),RCM(nmaxm,3),RNQ(nmaxm,3),TDIP(4)             
      common /boxs/ bx,by,bz,boxhx,boxhy,boxhz                                  
C                                                                              
      COMMON /BINFRAC/MOLA                            
      COMMON /BINMIX/ SIGAA,SIGBB,SIGAB,EPSAA,EPSAB,
     *                                EPSBB,AMASSA,AMASSB      
      COMMON /BINENE/CONREPAA,CONREFAA,CONREPBB,CONREFBB,
     *                                CONREPAB,CONREFAB,
     *               CONATPAA,CONATFAA,CONATPBB,CONATFBB,
     *                                CONATPAB,CONATFAB,
     *   AENEAA,AENEBB,AENEAB,BENEAA,BENEBB,BENEAB
      COMMON /BINRANGES/ RANGEAA,RANGEBB,RANGEAB,RANGEQAA,
     *      RANGEQBB,RANGEQAB
      
      COMMON /DIFFMASSES/ WAITA,WAITB,PARKEA,PARKEB,FACEKIA,FACEKIB

      COMMON /CMLISTAA/ LISTAA(NLISTMAX),LASTAA(nmaxm)
      COMMON /CMLISTBB/ LISTBB(NLISTMAX),LASTBB(nmaxm)
      COMMON /CMLISTAB/ LISTAB(NLISTMAX),LASTAB(nmaxm)
      COMMON /LOGTIME/ NSS(100),BASE,NPC

      COMMON /BKSPARAA/ QBKA,QBKAA,AMPAA,VVVAA,CCCAA,FOEXPAA,FORR6AA
      COMMON /BKSPARBB/ QBKB,QBKBB,AMPBB,VVVBB,CCCBB,FOEXPBB,FORR6BB
      COMMON /BKSPARAB/      QBKAB,AMPAB,VVVAB,CCCAB,FOEXPAB,FORR6AB
c
c
      COMMON /BKSPOTEN/ EPOTENCHARGES,VIRIALCHARGES,EPOTENEWK
      COMMON /BKSFORCES/  FCHARGES(nmaxm,3),FEWK(nmaxm,3)    
      COMMON /LJFORCES/ FLJ(NMAXM,3),EPOTENLJ,VIRIALLJ
      PARAMETER ( MAXK = 1000 )
      REAL*8      KVEC(MAXK), KAPPA
      COMMON /EWALD1/ KVEC,KAPPA  
      COMMON /BKSRANGES/ RANGECAA,RANGECBB,RANGECAB,RANGEQCAA,
     *      RANGEQCBB,RANGEQCAB
      COMMON /HARMRANGE/  RHARMBB,RHARMQBB, AMPHARMBB,EOFFSETBB,
     *    RHARMAB,RHARMQAB, AMPHARMAB,EOFFSETAB
      COMMON /LAMBDA/ ALAMBDA,DELTAEPOT
      COMMON /PATCHES/ SAA,SBB,SAB,EAA,EBB,EAB,EBB6TOT,EBB30TOT,
     *  FORR30BB,EAB6TOT,EAB30TOT,FORR30AB
      character*60 pwd
      character*20 pn
      COMMON /PERMPI/ PWD,PN
      character*50 filelabel
      character*132 fileinput
      COMMON /SWITCH / RS,RS2,RC,RC2,
     *                 VSWAAA,VSWBAA,VSWCAA,VSWAAB,VSWBAB,VSWCAB,
     *                 VSWABB,VSWBBB,VSWCBB,              
     *                 FSWAAA,FSWBAA,FSWCAA,FSWAAB,FSWBAB,FSWCAB,
     *                 FSWABB,FSWBBB,FSWCBB  


/* ============================= >>> fixswitch <<< ============================ */
void fixswitch(void)
{
  PI = 4.0*DATAN(1.0);
  TOTPI = 2.0/sqrt(PI);
  //  print *,' RS=',RS
  //  print *,' RC=',RC
  DR = RS - RC;
  /*
     this subroutine calculate the switchin parameters for VSW
     
     once RC AND RS ARE KNOW, FROM THE CONTINUITY OF V,V' and V''
     one GET THE THREE COEFFICIENT (ONE FOR EACH PAIR I,J)
     
     VSW= vswaij*(r-rc)**5 + vswbij*(r-rc)**4 + vswcij*(r-rc)**3
     
     AA interactions
   */
  POT = QBKAA* erfc(KAPPA*RS)/RS;
  DER= -QBKAA*
    (TOTPI*KAPPA*RS* exp(-KAPPA*RS*KAPPA*RS)+ erfc(KAPPA*RS))/RS/RS;
  SDE=QBKAA* (4.d0*KAPPA**3/sqrt(PI)* exp(-(KAPPA**2)
	*(RS**2)) +4.d0*KAPPA/sqrt(PI)*exp(-(KAPPA**2)
	*(RS**2))/RS**2 +2.d0*erfc(KAPPA*RS)/RS**3);

  //        print *,' POT,DER,SDE',POT,DER,SDE

  VSWAAA = (12.d0*POT-6.d0*DER*DR+SDE*DR**2)/2.d0/DR**5;
  VSWBAA = (-30.d0*POT+14.d0*DER*DR-2.d0*SDE*DR**2)/2.d0/DR**4;
  VSWCAA = (20.d0*POT-8.d0*DER*DR+SDE*DR**2)/2.d0/DR**3;
  /*
     AB interaction
   */
  RSINV= 1.0 / RS;
  POT= QBKAB* erfc(KAPPA*RS)/RS + AMPAB*exp(-VVVAB*RS)      
    + RSINV**6*(CCCAB + EAB6TOT) 
    + RSINV**30*EAB30TOT;
  DER = -QBKAB*
    (TOTPI*KAPPA*RS*exp(-KAPPA*RS*KAPPA*RS)+ erfc(KAPPA*RS))/RS/RS
     -   FOEXPAB*exp(-VVVAB*RS) - FORR6AB*RSINV**7
     - FORR30AB*RSINV**31; 
    
  SDE = QBKAB* (4.d0*KAPPA**3/sqrt(PI)*exp(-(KAPPA**2)
    *(RS**2)) +4.d0*KAPPA/sqrt(PI)*exp(-(KAPPA**2)
    *(RS**2))/RS**2 +2.d0*erfc(KAPPA*RS)/RS**3)+
    VVVAB**2*AMPAB*exp(-VVVAB*RS)
    + 42.0d0*RSINV**8*(CCCAB + EAB6TOT)
    + 930.0d0*RSINV**32*EAB30TOT;


  // print *,' POT,DER,SDE',POT,DER,SDE
  
  VSWAAB = (12.0*POT-6.0*DER*DR+SDE*DR**2)/2.0/DR**5;
  VSWBAB = (-30.0*POT+14.0*DER*DR
	    - 2.0*SDE*DR**2)/2.0/DR**4;
  VSWCAB = (20.0*POT-8.0*DER*DR+SDE*DR**2)/2.0/DR**3;

  /*
     BB INTERACTION
   */ 
  RSINV=1.0/RS;
  POT= QBKBB* erfc(KAPPA*RS)/RS + AMPBB*exp(-VVVBB*RS)      
    + RSINV**6*(CCCBB + EBB6TOT) 
    + RSINV**30*EBB30TOT;
  DER= -QBKBB*
    (TOTPI*KAPPA*RS*exp(-KAPPA*RS*KAPPA*RS)+erfc(KAPPA*RS))/RS/RS
    -   FOEXPBB*exp(-VVVBB*RS) - FORR6BB*RSINV**7
    - FORR30BB*RSINV**31;

  SDE = QBKBB* (4.0*KAPPA**3/sqrt(PI)*exp(-(KAPPA**2)
 	*(RS**2)) +4.0*KAPPA/sqrt(PI)*exp(-(KAPPA**2)
        *(RS**2))/RS**2 +2.d0*erfc(KAPPA*RS)/RS**3)+
        + VVVBB**2*AMPBB*exp(-VVVBB*RS)
        + 42.0*RSINV**8*(CCCBB + EBB6TOT)
        +930.0*RSINV**32*EBB30TOT;




  //        print *,' POT,DER,SDE',POT,DER,SDE

  VSWABB = (12.0*POT-6.0*DER*DR+SDE*DR**2)/2.0/DR**5;
  VSWBBB = (-30.0*POT+14.0*DER*DR
   	    -2.0*SDE*DR**2)/2.0/DR**4;
  VSWCBB = (20.0*POT-8.0*DER*DR+SDE*DR**2)/2.0/DR**3;

}

/* ========================== >>> setupewald <<< ========================== */
void setupewald(void)
{

  /*
   *******************************************************************
   ** ROUTINE TO SET UP THE WAVE-VECTORS FOR THE EWALD SUM.         **
   **                                                               **
   ** THE WAVEVECTORS MUST FIT INTO A BOX OF UNIT LENGTH.           **
   ** IN THIS EXAMPLE WE ALLOW A MAXIMUM OF 1000 WAVEVECTORS.       **
   *******************************************************************
   
   K IN UNITS OF TWOPI/BX
   */
  int TOTK;


  TWOPI=2.0*4.0*atan(1.0);
  B = 1.0 / 4.0 / KAPPA / KAPPA;
  TWOPIOVERBOX=TWOPI/BX;
  //    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **
  //    print *,' KAPPA=',KAPPA

  TOTK = 0;
  for (KX = 0; KX < KMAX; KX++)
    {
      RKX = TWOPIOVERBOX * KX;
      for (KY = -KMAX; KY < KMAX; KY++)
	{
     	  RKY = TWOPIOVERBOX * KY;
	  for (KZ = -KMAX; KZ < KMAX; -KZ++)
	    {
	      RKZ = TWOPIOVERBOX * KZ;
	      KSQ = KX * KX + KY * KY + KZ * KZ;
	      if ( ( KSQ < KSQMAX ) && ( KSQ != 0 ) ) 
		{
	       	  TOTK = TOTK + 1;
		  if ( TOTK > MAXK ) 
		    {
		      printf("KVEC is too small\n");
		      exit(-1);
		      //STOP 'KVEC IS TOO SMALL';
		    }
		  
		  RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ;
		  KVEC(TOTK) = TWOPI * exp(-B*RKSQ)/RKSQ /BX**3;
		}
	    }
	}
    }
  //WRITE( *, ' ( '' EWALD SUM SETUP COMPLETE ''     ) ' )
  //  WRITE( *, ' ( '' NUMBER OF WAVEVECTORS IS '', I5 ) ' ) TOTK
}
 
/* ============================= >>> initalcharges <<< ===========================*/
/* 
 ******************************************************************       

 Data for the potential  v(r)=q*q/r+amp*exp(-r*vvv)+ccc/r**8
 
 +4*eab*[ (Sab/r)^30 - (Sab/r)^6]

 (information on masses are in INITAL)
 
 Charges are in a.u.  
 vvv in nm
 amp in ?   --- 
 */

void  INITALCHARGES(void)
{                                  
  
  /*
     Il potenziale e' zero for r>rc
     e' switching function r>rs
     
   */
  PI = 4.0*atan(1.0);
  RS=0.77476;
  RC=1.0;        
  RS2=RS*RS;
  RC2=RC*RC;
  //                             RECIPROCAL SPACE K
  KAPPA=2.5;

  //        print *,'KAPPA',kappa

  //                      to go to  q*q/r in kJ/mol with r in nm

  QFACTOR = 11.7870740;
           
  AVN = 6.02205E23;
  //              1ev in coulomb
  EVINC = 1.60219E-19;
  //                      to go to  q*q/r in kJ/mol with r in nm
  qpieps0 = 4.0 * PI * 8.854188E-12;

//               da eV a kJ/mol        
  EFACTOR = EVINC*AVN/1000.0;
  Qfactor = (1.0/qpieps0*evinc*evinc/(1.0E-9)/ 1000.0 * avn)**0.5;
  //  print *,' Qfactor',qfactor       
        
  /*                                                              
	  SPECIE # 1  Si-Si
   */
  QBKA  = 2.4;
  QBKA  = QBKA*QFACTOR;
  QBKAA = QBKA*QBKA;
  AMPAA = 0.0;
  VVVAA = 0.0;
  CCCAA = 0.0;
  EAA   = 0.0;
  SAA   = 0.0;
  RANGECAA = RC;
  RANGEQCAA = RANGECAA*RANGECAA;

  /*
     SPECIE # 2  (BB) Oxygen-Oxygen
   */
  QBKB = -1.2;
  QBKB = QBKB*QFACTOR;
  QBKBB= QBKB*QBKB;
  //        print *,' ampbb da peter s code (J) ',2.224814e-16
  //            questa e' in eV
  AMPBB = 1388.7730;
  AMPBB = AMPBB*EFACTOR;
  /*        print *,'Amp OO kJ/mol',ampbb
	  print *,'Amp in J',ampbb/avn*1000.0d0 */
  VVVBB = 27.60;
  /*        print *,' C O-O da peters code', -2.80350e-23*/
  CCCBB = -175.0000E-6;
  CCCBB = CCCBB*EFACTOR;
  //                                  ebb in J (as given by Peter's student)
  EBB = 1.6839685E-22;
  SBB = 0.1779239;
  //                        E o-o in kJ/mol
  EBB=EBB*AVN/1000.0;
  /*
     print *,' CCC in kJ/mol nm6',  CCCBB
     print *,' CCC in J nm6',  CCCBB/avn*1000.0d0*/
  FOEXPBB = AMPBB*VVVBB;
  FORR6BB = 6.0*(CCCBB- 4.0*EBB*SBB**6);
  FORR30BB = 30.0*(4.0*EBB*SBB**30);
  EBB30TOT=4.0*EBB*SBB**30;
  EBB6TOT=-4.0*EBB*SBB**6;
  RANGECBB=RC;
  if (RANGECBB > BX*0.5) 
    {  
      printf("ERROR RANGECBB rinormalizzo ?\n");
      RANGECBB=BX*0.5;
    }

  RANGEQCBB=RANGECBB*RANGECBB;
  /* 
     AB interactions Si-Oxygen
   */
  QBKAB=QBKA*QBKB;
  /* OLD:        AMPAB=2.884201e-15
     print *,' ampbb da peter s code',ampab
     questa e' in eV */
  AMPAB = 18003.7572;
  AMPAB=AMPAB*EFACTOR;
  /*        print *,'Amp SiO kJ/mol',ampab
       	    print *,'Amp AB in J',ampab/avn*1000.0d0 */
  VVVAB = 48.7318;
  /* print *,' C Si-O da peters code', -2.13928e-23 */
  CCCAB = - 133.5381E-6;
  CCCAB = CCCAB*EFACTOR;	
  /*  E si-o in J */
  EAB = 4.9634598E-22;
  /*                         E si-o in kJ/mol */
  EAB = EAB*AVN / 1000.0;
  SAB= 0.1313635;
  /*  
      print *,' CCC in kJ/mol nm6',  CCCAB
      print *,' CCC in J nm6',  CCCAB/avn*1000.0d0 */
  FOEXPAB = AMPAB*VVVAB;
  FORR6AB = 6.0*(CCCAB-4.0*EAB*SAB**6);

  FORR30AB = 30.0 * (4.02*EAB*SAB**30);
  EAB30TOT = 4.0*EAB*SAB**30;
  EAB6TOT  = - 4.0*EAB*SAB**6;

  RANGECAB = RC;
  if (RANGECAB > BX*0.5) 
    {
      printf(" ERROR RANGECAB - RINORMALIZZO RANGE\n");
      RANGECAB=BX*0.5;
    }
  RANGEQCAB=RANGECAB*RANGECAB;
  /*
     Now fix the swithcing function
     VSW= vswaij*(r-rc)**5 + vswbij*(r-rc)**4 + vswcij*(r-rc)**3
   */
  FIXSWITCH(void);
  /*
     per le forze servono
   */
  FSWAAA= - 5.0*VSWAAA;
  FSWBAA= - 4.0*VSWBAA;	
  FSWCAA= - 3.0*VSWCAA;
  FSWAAB= - 5.0*VSWAAB;	
  FSWBAB= - 4.0*VSWBAB;
  FSWCAB= - 3.0*VSWCAB;
  FSWABB= - 5.0*VSWABB;
  FSWBBB= - 4.0*VSWBBB;
  FSWCBB= - 3.0*VSWCBB;

  SETUPEWALD(void);
     
}      


/* =========================== >>> inital <<< ============================= */
void inital(void)                                  
{
  const double AVOGA = 6.02205E+02, PI = 3.1415926535898, BOLTK= 8.314E-03;

  NLIST = NLISTMAX;
  /*
     This is for 33% 66% mixture
   */
  MOLA = MOLS/3;
  //      print *,' #A,#B,#(A+B)',MOLA,MOLS-MOLA,MOLS
       
  //     print *,'BOX SIZE=',bx,by,bz
  bmin = 95001.0;
  if (bx < bmin) 
    bmin = bx;
  if (by < bmin) 
    bmin = by;
  if (bz < bmin) 
    bmin = bz; 
  bmin=bmin/2.0;      
  boxhx = 0.5*bx;
  boxhy = 0.5*by; 
  boxhz = 0.5*bz;
/*                                                              
  ARGON EPS=0.99605728D0,SIGMA=0.341D0, MASSA=39.95

  SPECIE # 1

*/                 
  EPSON = 23.0;
  SIGMA = 0.33;
  WAITA = 28.085500; 
  SIGAA = SIGMA;
  EPSAA = EPSON;

  //  print ("Time Units=%e\n",dsqrt(wait*sigaa**2/48.0d0/epsaa);

  CL = 2.50;
  RSWITCHAA = CL*SIGMA;
  RSWITCH2AA=RSWITCHAA*RSWITCHAA;
  /*
     EPSONR=EPSON RINORMALIZZATO
     SIGMAR=SIGMA RINORMALIZZATA 
   */
  EPSONR=EPSON;
  SIGMAR=SIGMA;
  CONREPAA= 4.00*EPSONR*SIGMAR**12;
  CONATPAA= 4.00*EPSONR*SIGMAR**6;
  CONREFAA=12.00*CONREPAA;
  CONATFAA= 6.00*CONATPAA;
  RANGEAA=   RSWITCHAA;
  RANGEQAA=RANGEAA*RANGEAA;
  AENEAA= -(CONREPAA/RSWITCHAA**12-CONATPAA/RSWITCHAA**6);
  //  print *,' Energy Shift AA=',AENEAA
  /*
       SPECIE # 2  (BB)
   */
                 
  EPSON = 23.0;
  SIGMA = 0.28;
  WAITB = 15.999400;
  EPSBB = EPSON;
  SIGBB=SIGMA;

  CL=2.50;
  EPSONR=EPSON;
  SIGMAR=SIGMA;
  RSWITCHBB=CL*SIGMA;
  RSWITCH2BB=RSWITCHBB*RSWITCHBB;
  /*
     EPSONR=EPSON RINORMALIZZATO
     SIGMAR=SIGMA RINORMALIZZATA 
   */
  CONREPBB = 4.00*EPSONR*pow(SIGMAR,12);
  CONATPBB = 4.00*EPSONR*SIGMAR**6;
  CONREFBB = 12.00*CONREPBB;
  CONATFBB = 6.00*CONATPBB;
  RANGEBB  =  RSWITCHBB;
  RANGEQBB = RANGEBB*RANGEBB;
  AENEBB=-(CONREPBB/RSWITCHBB**12-CONATPBB/RSWITCHBB**6);
  //print *,' Energy Shift BB=',AENEBB;
  /*
     AB interactions 
   */
             
  EPSON = 32.0;
  SIGMA = 0.160;
  EPSONR = EPSON;
  SIGMAR = SIGMA;
  EPSAB = EPSON;
  SIGAB=SIGMA;
  CL = 2.50;
  RSWITCHAB = CL*SIGMA;
  RSWITCH2AB = RSWITCHAB*RSWITCHAB;
/*
   EPSONR=EPSON RINORMALIZZATO
   SIGMAR=SIGMA RINORMALIZZATA */

  CONREPAB = 4.00*EPSONR*pow(SIGMAR,12);
  CONATPAB = 4.00*EPSONR*pow(SIGMAR,6);
  CONREFAB = 12.00*CONREPAB;
  CONATFAB = 6.00*CONATPAB;
  RANGEAB  = RSWITCHAB;
  RANGEQAB = RANGEAB*RANGEAB;
  AENEAB=-(CONREPAB/pow(RSWITCHAB,12)-CONATPAB/pow(RSWITCHAB,6));
  //      print *,' Energy Shift AB=',AENEAB

  range=max(rangeaa,rangeab,rangebb);

  FMOLS = ( (double) MOLS);
  if (range > bmin) 
    range = bmin;
  printf("Using range =%.f\n",range);                                
  RANT  = 0.95*RANGE;
  // era 1.1
  RANG2  = 1.150*RANGE;
  RANGE3 = 1.00/pow(RANGE,3);
  RANG3H = 0.50*RANGE3; 
  PERMIT = 0.25*pow((RANG2-RANGE),2);
          
  RANTQ = RANT**2;
  RANGEQ= RANGE**2;
  RANG2Q= RANG2**2;
  VOLUME= BX*BY*BZ;
                                                        
  FACTMP(1)= 2.0/( 3.0*BOLTK );
  FACTMP(2)= 1.0/( 3.0*BOLTK );                                   
  FACEKIA= 0.5*WAITA/FMOLS;                                 
  FACEKIB= 0.5*WAITB/FMOLS;
/*                                                  
    TO be controlled                  
 */
  FACPRE= 2.0E+03*FMOLS/( 3.0*AVOGA );                     
  FACWOR= 1.0E-03*AVOGA/FMOLS;
  FACVIR= 0.5/FMOLS;
  FACDENA= WAITA*FMOLS/AVOGA;
  FACDENB= WAITB*FMOLS/AVOGA;                                        
  FACGVR= DELTR*VOLUME/( 2.0*PI*Sqr(FMOLS) );  
}
/*
  This calculate the real space part of the Ewald sum
  the exp and the LJ terms

 ****************************************************************** */
void ewrforcescharges(void)
{                                                              
      double RIJ[3],BBR[3],CCEL[3];//,FO(nmaxm,3)           
      
      PI = 4.0 * atan(1.0);
      TOTPI = 2.0 / sqrt(PI);

      dminab=1E10;
      dminbb=1E10;

      EPOTENCHARGES=0.0;
      VIRIALCHARGES=0.0;                      
                
        
      /*  
	  AA interactions - ALL OF THEM !
	  
       */
    for (I = 1; i < Oparams.parnum[0]; i++)   
      for (J = I+1; j < Oparams.parnum[1]; j++ )
	{
	  
	  rij[1] = rx[0][i] - rx[0][j];
	  rij[2] = ry[0][i] - ry[0][j];
	  rij[3] = rz[0][i] - rz[0][j];
	  
	  /* riduzione al first box */
	  BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx));
	  BBR(2)=By*DBLE(INT(RIJ(2)/BOXHy));  
	  BBR(3)=Bz*DBLE(INT(RIJ(3)/BOXHz));                                       
	  RIJ(1)=RIJ(1)-BBR(1);                
	  RIJ(2)=RIJ(2)-BBR(2);
	  RIJ(3)=RIJ(3)-BBR(3);
	  RO = Sqr(RIJ(1))+Sqr(RIJ(2))+Sqr(RIJ(3));   


	  if (RO <= RC2)  
	    {

	      RQINV =1.0D0/RO;                                           
	      RIJP=DSQRT(RO);
	      RIJINV= 1.0D0/RIJP;
		      
	      if (RO.LT.RS2) 
	       	{	  
		  ALPHAR= KAPPA*RIJP;             
		  ERC = erfc(ALPHAR); 
		  POTDM =  QBKAA*ERC/RIJP;
		  CCEL3 =  QBKAA*
		    * (TOTPI*ALPHAR*exp(-ALPHAR*ALPHAR)+ERC)*RQINV*RIJINV;
		  
		}
	      else
		{		
		  
		    DR1=RIJP-RC;
	       	    DR2=DR1*DR1;
      		    DR3=DR2*DR1;
		    DR4=DR3*DR1;
		    DR5=DR4*DR1;
		    
		    POTDM=VSWAAA*DR5+VSWBAA*DR4+VSWCAA*DR3;
		    CCEL3=(FSWAAA*DR4+FSWBAA*DR3+FSWCAA*DR2)*RIJINV;

		}

      EPOTENCHARGES=EPOTENCHARGES+POTDM;

      CCEL(1)=CCEL3*RIJ(1);
      CCEL(2)=CCEL3*RIJ(2);
      CCEL(3)=CCEL3*RIJ(3);

      //      FO(I,1)=FO(I,1)+CCEL(1)
      //      FO(I,2)=FO(I,2)+CCEL(2)
      //      FO(I,3)=FO(I,3)+CCEL(3)
                                                       
      //FO(J,1)=FO(J,1)-CCEL(1)
      //FO(J,2)=FO(J,2)-CCEL(2)
      //FO(J,3)=FO(J,3)-CCEL(3)
                                 
      VIRIALCHARGES=VIRIALCHARGES-CCEL(1)*RIJ(1) 
	-CCEL(2)*RIJ(2)
	-CCEL(3)*RIJ(3);
	    }
	}
      
    /* Ora BB */  
      for (i=MOLA+1; i < MOLS-1;i++)
	for (j = i+1, j < mols; j++)
	  {
	    RIJ(1)=RN(I,1)-RN(J,1);
	    RIJ(2)=RN(I,2)-RN(J,2);                                               
	    RIJ(3)=RN(I,3)-RN(J,3);                                                
                                                        
      	    BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHX));
	    BBR(2)=By*DBLE(INT(RIJ(2)/BOXHY));
	    BBR(3)=Bz*DBLE(INT(RIJ(3)/BOXHZ));
	    
      	    RIJ(1)=RIJ(1)-BBR(1);
	    RIJ(2)=RIJ(2)-BBR(2);
	    RIJ(3)=RIJ(3)-BBR(3);
	    
      	    RO=Sqr(RIJ(1))+Sqr(RIJ(2))+Sqr(RIJ(3));
	    if (RO > RC2) 
	      goto 1340

		RQINV =1.0/RO;                                           
		RIJP=sqrt(RO);
		RIJINV= 1.0/RIJP;
		
		if (RO < RS2) 
		  {
		    
		    ALPHAR= KAPPA*RIJP;
		      ERC = DERFC(ALPHAR);
		      POTDM =  QBKBB*ERC/RIJP;
		      CCEL3 =  QBKBB*
			* (TOTPI*ALPHAR*exp(-ALPHAR*ALPHAR)+ERC)*RQINV*RIJINV;
		      TSN   =  pow(RQINV,3);
		      BBEXP =  exp(-VVVBB*RIJP);
		      POTDM =  POTDM  + AMPBB*BBEXP   + TSN*(CCCBB + EBB6TOT) +
                                         pow(RQINV,15)*EBB30TOT;
		      CCEL3 =  CCEL3  + FOEXPBB*BBEXP*RIJINV 
      			+ FORR6BB*TSN*RQINV
       			+ FORR30BB*pow(RQINV,16);
	     	  }
		else
		  {
		    DR1=RIJP-RC;
		    DR2=DR1*DR1;	
		    DR3=DR2*DR1;
		    DR4=DR3*DR1;
		    DR5=DR4*DR1;

	    	    POTDM=VSWABB*DR5+VSWBBB*DR4+VSWCBB*DR3;
		    CCEL3=(FSWABB*DR4+FSWBBB*DR3+FSWCBB*DR2)*RIJINV;

		  }

	 	EPOTENCHARGES=EPOTENCHARGES+POTDM;
		
		CCEL(1)=CCEL3*RIJ(1);                                                    
		CCEL(2)=CCEL3*RIJ(2);
		CCEL(3)=CCEL3*RIJ(3);                                          
                                                                               
		//FO(I,1)=FO(I,1)+CCEL(1);
		//FO(I,2)=FO(I,2)+CCEL(2);
		//FO(I,3)=FO(I,3)+CCEL(3);
                                                                               
	  	//FO(J,1)=FO(J,1)-CCEL(1);
		//FO(J,2)=FO(J,2)-CCEL(2);
		  
		//FO(J,3)=FO(J,3)-CCEL(3);                                     
                                                                               
	  	VIRIALCHARGES=VIRIALCHARGES-CCEL(1)*RIJ(1) 
		  -CCEL(2)*RIJ(2)
		  -CCEL(3)*RIJ(3);
               
	  }

      /*
	 ORA AB
       */
      for (i=1; i < MOLA; i++)
	for (j = mola+1; j < mols; j++)
	  {
	    RIJ(1)=RN(I,1)-RN(J,1);
	    RIJ(2)=RN(I,2)-RN(J,2);                                                
	    RIJ(3)=RN(I,3)-RN(J,3);                                                
                                                                               
      	    BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx));
	    BBR(2)=BY*DBLE(INT(RIJ(2)/BOXHy));
	    BBR(3)=BZ*DBLE(INT(RIJ(3)/BOXHz));                                  

                                                                               
      	    RIJ(1)=RIJ(1)-BBR(1);
	    RIJ(2)=RIJ(2)-BBR(2);                                                      
	    RIJ(3)=RIJ(3)-BBR(3);                                                     
                                                                               
      	    RO=Sqr(RIJ(1))+Sqr(RIJ(2))+Sqr(RIJ(3));
	    if (RO > RC2) 
	      goto 2340;

	    RQINV =1.0/RO;                                      
	    RIJP=sqrt(RO);
	    RIJINV= 1.0/RIJP;
	    if (RO < RS2) 
	      {
		ALPHAR= KAPPA*RIJP             
		  ERC = DERFC(ALPHAR) 
		  POTDM =   QBKAB*ERC/RIJP
		  CCEL3 =  QBKAB*
		  * (TOTPI*ALPHAR*DEXP(-ALPHAR*ALPHAR)+ERC)*RQINV*RIJINV 
		  c
		  c  exp and r-n parts of the potential
		  c 
		  TSN   =  pow(RQINV,3) 
		  ABEXP =  DEXP(-VVVAB*RIJP)
		  POTDM =  POTDM + AMPAB*ABEXP + TSN*(CCCAB + EAB6TOT) +
		  *                                    pow(RQINV,15)*EAB30TOT
		  CCEL3 =  CCEL3 + FOEXPAB*ABEXP*RIJINV + FORR6AB*TSN*RQINV
		  *                  + FORR30AB*pow(RQINV,16)
	      }
	    else
	      {
		DR1=RIJP-RC;
		DR2=DR1*DR1;
		DR3=DR2*DR1;
		DR4=DR3*DR1;
		DR5=DR4*DR1;

	      	POTDM=VSWAAB*DR5+VSWBAB*DR4+VSWCAB*DR3;
		CCEL3=(FSWAAB*DR4+FSWBAB*DR3+FSWCAB*DR2)*RIJINV;
	      }
	    

      	    EPOTENCHARGES=EPOTENCHARGES+POTDM;

      	    CCEL(1)=CCEL3*RIJ(1);
	    CCEL(2)=CCEL3*RIJ(2);                                                      
	    CCEL(3)=CCEL3*RIJ(3);                                                      
                                                                               
	    // FO(I,1)=FO(I,1)+CCEL(1)                                                   
	    // FO(I,2)=FO(I,2)+CCEL(2) 
	    // FO(I,3)=FO(I,3)+CCEL(3)                                                   
                                                                               
      	    //FO(J,1)=FO(J,1)-CCEL(1)                                                   
	    //FO(J,2)=FO(J,2)-CCEL(2)                                                   
	    //FO(J,3)=FO(J,3)-CCEL(3)                                                   
                                                                               
      	    VIRIALCHARGES=VIRIALCHARGES-CCEL(1)*RIJ(1) 
	      -CCEL(2)*RIJ(2)
	      -CCEL(3)*RIJ(3);
	  }

      /*
	 DO 180 J=1,3                                                              
	 DO 180 I=1,MOLS                                                           
	 C                     
	 FCHARGES(I,J)= FO(I,J)                                               
	 C                                                                               
	 180 CONTINUE 
       */
}


/* 
     *******************************************************************
     ** CALCULATES K-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
     **                                                               **
     ** THE SELF TERM IS SUBTRACTED.                                  **
     ** IN ONE COORDINATE DIRECTION (X), SYMMETRY IS USED TO REDUCE   **
     ** THE SUM TO INCLUDE ONLY POSITIVE K-VECTORS.                   **
     ** THE NEGATIVE VECTORS IN THIS DIRECTION ARE INCLUDED BY USE    **
     ** OF THE MULTIPLICATIVE VARIABLE 'FACTOR'.                      **
     **                                                               **
     ** PRINCIPAL VARIABLES:                                          **
     **                                                               **
     ** INTEGER N                   NUMBER OF IONS                    **
     ** REAL    RX(N),RY(N),RZ(N)   POSITIONS OF IONS                 **
     ** REAL    Z(N)                IONIC CHARGES                     **
     ** REAL    VK                  K-SPACE POTENTIAL ENERGY          **
     ** REAL    VKS                 SELF PART OF K-SPACE SUM          **
     ******************************************************************* 
*/

void ewkforcescharges(void)    
{
      int     TOTK;
      int     KX, KY, KZ, I, KSQ;

      COMPLEX*16     EIKX(1:NMAXM, 0:KMAX);
      COMPLEX*16     EIKY(1:NMAXM, -KMAX:KMAX);
      COMPLEX*16     EIKZ(1:NMAXM, -KMAX:KMAX);
      COMPLEX*16     EIKR(NMAXM), SUM, SUMA, SUMB;
      // DIMENSION   FO(nmaxm,3)           
      
      PI=4.0*atan(1.0);
      TWOPI=2.0*PI;
      TWOPIOVERBOX=TWOPI/BX;
      RSQPI=1.0d0/sqrt(PI);
      /*
	 DO 50 J=1,3  
	 DO 50 I=1,MOLS                                
	 FO(I,J)  =0.0D0                                
	 50 CONTINUE   */
      /*

       *******************************************************************
       
       ** CONSTRUCT EXP(IK.R) FOR ALL IONS AND K-VECTORS **
       ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY **
       */
      for (i=1,i < mols; i++)
	{
     
  	  EIKX(I, 0) = (1.0d0, 0.0d0);
	  EIKY(I, 0) = (1.0d0, 0.0d0);
	  EIKZ(I, 0) = (1.0d0, 0.0d0);

	  EIKX(I, 1) = DCMPLX ( cos ( TWOPIOVERBOX * RN(I,1) ) ,
		sin ( TWOPIOVERBOX * RN(I,1) ) );
	  EIKY(I, 1) = DCMPLX ( DCOS ( TWOPIOVERBOX * RN(I,2) ) ,
				DSIN ( TWOPIOVERBOX * RN(I,2) ) );
	  EIKZ(I, 1) = DCMPLX ( DCOS ( TWOPIOVERBOX * RN(I,3) ) ,
				DSIN ( TWOPIOVERBOX * RN(I,3) ) );

	  EIKY(I, -1) = DCONJG ( EIKY(I, 1) );
       	  EIKZ(I, -1) = DCONJG ( EIKZ(I, 1) );
	}

      //    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **

      for (kx = 2; kx < kmax; kx++)
	for (i = 1; i < mols; i++)
	  {
	    EIKX(I, KX) = EIKX(I, KX-1) * EIKX(I, 1)
	  }

        for  (ky = 2; ky < KMAX; ky++)
	  for (i = 1; i < mols; i++)
	    { 
	      EIKY(I,  KY) = EIKY(I, KY-1) * EIKY(I, 1);
	      EIKY(I, -KY) = DCONJG ( EIKY(I, KY) );
	    }
        
	 for (kz = 2; kz < KMAX; kz++)
	   for (i = 1; i < mols; i++)
	     { 
	       EIKZ(I,  KZ) = EIKZ(I, KZ-1) * EIKZ(I, 1);
	       EIKZ(I, -KZ) = DCONJG ( EIKZ(I, KZ) );
	     }
	 /*    ** SUM OVER ALL VECTORS **
	       */

        VD   = 0.0;
        TOTK = 0;
	
        for(kx = 0; kx < KMAX; kx++)
	  {
	    // FACTOR tiene conto della riflession k-> -k 
	    if ( kx == 0 ) 
              FACTOR = 1.0;  
	    else
	      FACTOR = 2.0;
	    
           for(ky = -KMAX; ky <= KMAX; ky++)
	     for(kz = -KMAX; kz < KMAX; kz++)
	       {
		 KSQ = kx * kx + ky * ky + kz * kz;
		 if ( ( KSQ < KSQMAX ) && ( KSQ != 0 ) ) 
		   {
		     TOTK = TOTK + 1;
		     SUMA  = (0.0, 0.0);
		     for (i = 1; i < MOLA; i++)
		       {
			 EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ);
			 SUMA     = SUMA + EIKR(I);
		       }
		     SUMA=SUMA*QBKA;
		     SUMB = (0.0, 0.0);
		     for (i = MOLA + 1; i < MOLS; i++ )
		       {
    			 EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ);
			 SUMB     = SUMB +  EIKR(I);
		       }
		     SUMB=SUMB*QBKB;
		     SUM = SUMA+SUMB;
		     VD = VD + FACTOR * KVEC(TOTK) * DCONJG (SUM)*SUM;
		     /*
			Here forces  (see notes on NoteBook)
			
			DO  I = 1, MOLA 
			QFORCE= 2.0d0*
			c        FACTOR*KVEC(TOTK)*QBKA*DIMAG(EIKR(I)*DCONJG(SUM)) 
			FO(I,1)=FO(I,1)+QFORCE* TWOPIOVERBOX* KX 
			FO(I,2)=FO(I,2)+QFORCE* TWOPIOVERBOX* KY 
			FO(I,3)=FO(I,3)+QFORCE* TWOPIOVERBOX* KZ 
			END DO
			DO  I = MOLA+1,MOLS 
			QFORCE= 2.0d0*
			c        FACTOR*KVEC(TOTK)*QBKB*DIMAG( EIKR(I)*DCONJG(SUM)) 
			FO(I,1)=FO(I,1)+QFORCE* TWOPIOVERBOX* KX 
			FO(I,2)=FO(I,2)+QFORCE* TWOPIOVERBOX* KY 
			FO(I,3)=FO(I,3)+QFORCE* TWOPIOVERBOX* KZ 
			END DO    
			
			ENDIF
			
		      */
		   }
	       }
	  }

	/*    ** CALCULATES SELF PART OF K-SPACE SUM **
	             VS= SUM Q*Q
	      */
	VS= MOLA*QBKAA + (MOLS-MOLA)*QBKBB;
        VS = RSQPI * KAPPA * VS;

	//    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **

        VK = VD - VS;

	/*        write(88,*) rn(mols,1),vk,-fo(mols,1)
	 */
	/*       print *,' EPOTEN EWK ',VK
		 
		 Da decidere dove sia meglio sommare le F sulle cariche !
	 */
	/*
	   DO 180 J=1,3                                                              
 	   DO 180 I=1,MOLS                                                           
 	   
 	   FEWK(I,J)= FO(I,J)                                               
 	   
 	   180 CONTINUE                                                                  
	 */
        EPOTENEWK=VK;
}

extern int minimumImageI(int);

/* ================================= */

extern inline int minimumImageI(int pl);

/* ============================= >>> FCC <<< ================================*/
void FCC(COORD_TYPE* m)
{
  /*   DESCRIPTION:
       Sets up the alpha fcc lattice for n linear molecules.   
       The simulation box is a unit cube centred at the origin.
       N should be an integer of the form ( 4 * ( Nc ** 3 ) ),
       Where Nc is the number of FCC unit cells in each direction.  
       See figure 5.10 for a diagram of the lattice and a           
       definition of the four orientational sublattices.            
       PRINCIPAL VARIABLES:                                         
       COORD_TYPE    rxCm, ryCm, rzCm     Molecular Center of mass 
                                          positions             
       COORD_TYPE    rRoot3               1.0 / sqrt ( 3.0 ) */
  int Nc, Nm;
  int modA, modB, mod;
  COORD_TYPE rRoot3; // = 0.5773503;
  COORD_TYPE  Cell, Cell2, rxCm, ryCm, rzCm;
  int np[NA], a, i, ix, iy, iz, iref, ii;
  COORD_TYPE bx[4], by[4], bz[4]; /* base vectors for FCC lattice */
  double la;

  //printf("FCC Vol: %f\n", Vol);
  L = cbrt(Vol);
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  Nc = ceil(  pow( ((COORD_TYPE)Nm)/4.0, 1.0/3.0 )  );
  //printf("Nc: %d\n", Nc);
  /* Calculate the side of the unit cell */
  Cell  = L / ((COORD_TYPE) Nc); /* unit cell length */
  Cell2 = 0.5 * Cell;              /* half unit cell length */

  /* Sublattice A */
  rRoot3 = 1.0 / sqrt(3.0);
  bx[0] =  0.0;
  by[0] =  0.0;
  bz[0] =  0.0;
  /*  Sublattice B */
  bx[1] =  Cell2;
  by[1] =  Cell2;
  bz[1] =  0.0;
  /* Sublattice C */
  bx[2] =  0.0;
  by[2] =  Cell2;
  bz[2] =  Cell2;
  /* Sublattice D */
  bx[3] =  Cell2;
  by[3] =  0.0;
  bz[3] =  Cell2;
  /* Construct the lattice from the unit cell */
  
  ii = 0;
  np[0] = 0;
  np[1] = 0;

  la = Oparams.lattice_a;

  if (Oparams.parnum[0] > Oparams.parnum[1])
    {
      modA = rint(Oparams.parnum[0] / Oparams.parnum[1]);
      modB = 1;
      mod  = modA + modB; 
      /* Piazzera' modA atomi A(0) e poi 1 atomo B(1) */
    }
  else
    {
      modB = rint(Oparams.parnum[1] / Oparams.parnum[0]);
      modA = 1;
      mod  = modA + modB;
      /* Piazzera' modB atomi B(1) e poi 1 atomo A(0) */
    }

  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  //printf("rx[1][50]: %f\n", rx[1][39]);
  //printf("NA: %d NB:%d\n", Oparams.parnum[0], Oparams.parnum[1]);
  //printf("modA: %d modB: %d mod: %d\n", modA, modB, mod);
  loop(iz, 1, Nc) /* loops over unit cells (that are simply cubes) */ 
    {
      loop(iy, 1, Nc)
	{
	  loop(ix, 1, Nc)
	    {
	      loop(iref, 1, 4) /* In each primitive cell there are four 
				  molecules */
		{

		  /* Questo vuol dire che ho piazzato tutte le particelle */
		  if ( (np[0] >= Oparams.parnum[0]) &&
		       (np[1] >= Oparams.parnum[1]) )
		    break;
		  
		  /* nuova possibile posizione */
		  rxCm = bx[iref] + Cell * ((COORD_TYPE) ix);
		  ryCm = by[iref] + Cell * ((COORD_TYPE) iy);
		  rzCm = bz[iref] + Cell * ((COORD_TYPE) iz);
		  //printf("CM:(%f,%f,%f)\n", rxCm, ryCm, rzCm);
		  
		  if (modA == 1) 
		    {
		      /* questa condizione e' vera una volta ogni mod volte
			 in questo modo mette un atomo A e poi 
			 modB atomi B */
		      if ( (ii + iref) % mod == 0 )
			{
			  a = 0;
			}
		      else
			{
			  a = 1;
			}
		    }
		  else//qui se modB == 1
		    {
		      /* questa condizione e' vera una volta ogni mod volte
			 in questo modo mette un atomo B e poi 
			 modA atomi A */
		      if ( (ii + iref) % mod == 0 )
			{
			  a = 1;
			}
		      else
			{
			  
			  a = 0;
			}
		    }
		  /* se gli atomi di tipo a (0 o 1) sono stati tutti 
		     piazzati allora mette quelli dell'altro tipo */
		  if (np[a] == Oparams.parnum[a])
		    {
		      a = (~a) & 1;
		    }
		  //printf("np[%d]=%d\n ", a, np[a]);
		  rx[a][np[a]] = ((int)rint(rxCm/la));
		  ry[a][np[a]] = ((int)rint(ryCm/la));
		  rz[a][np[a]] = ((int)rint(rzCm/la));
		  //printf("(%d,%d,%d|%.5f)\n", rx[a][np[a]], 
		  // ry[a][np[a]], rz[a][np[a]], la);
		    
		  np[a]++;
		}
	      ii = ii + 4;
	    }
	  
	}
      
    }
  
  /* Shift centre of box to the origin */

  for ( a = 0; a < NA; a++)
    {
      for ( i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* Initial position values are between -0.5 and 0.5 */
	  rx[a][i] = minimumImageI(rx[a][i]); 
	  ry[a][i] = minimumImageI(ry[a][i]);
	  rz[a][i] = minimumImageI(rz[a][i]);
	  //printf("(%d,%d,%d)\n", rx[a][i], 
	  //  	 ry[a][i], rz[a][i]);
		  
	}
    }
  return;
}

/* ============================ >>> ranf <<< =============================== */
COORD_TYPE ranf(void)
{
  /*  Returns a uniform random variate in the range 0 to 1.         
      Good random number generators are machine specific.
      please use the one recommended for your machine. */
  return rand() / ( (COORD_TYPE) RAND_MAX );
}

/* ============================= >>> gauss <<< ============================= */
COORD_TYPE gauss(void)
{
  
  /* 
     Random variate from the standard normal distribution.
     
     The distribution is gaussian with zero mean and unit variance.
     REFERENCE:                                                    
                                                                
     Knuth D, The art of computer programming, (2nd edition        
     Addison-Wesley), 1978                                      
                                                                
     ROUTINE REFERENCED:                                           
                                                                
     COORD_TYPE ranf()                                  
     Returns a uniform random variate on the range zero to one  
  */

  COORD_TYPE  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  COORD_TYPE sum, r, r2;
  int i;

  sum = 0.0;

  loop(i, 1, 12)
    {
      sum = sum + ranf();
    }

  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;
}

/* ========================= >>> resetCM <<< ==============================*/
void resetCM()
{
  COORD_TYPE RCMx, RCMy, RCMz;
  COORD_TYPE *m, MTOT;
  double la;
  int i, a;

  m = Oparams.m;
  /* Remove net momentum, to have a total momentum equals to zero */
  la = Oparams.lattice_a;

  MTOT = 0.0;

  /* ADD 27/1/1998:
     And Now we put the center of mass of the box in the origin of axis
     because otherwise int NPT method the total momentum is not zero */
  RCMx = 0.0;
  RCMy = 0.0;
  RCMz = 0.0;

  for (a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  RCMx += m[a]*la*((double)rx[a][i]); 
	  /*Here RCM is the center of mass of the box */
	  RCMy += m[a]*la*((double)ry[a][i]);
	  RCMz += m[a]*la*((double)rz[a][i]);
	  MTOT += m[a];
	}
    }
  RCMx /= MTOT;
  RCMy /= MTOT;
  RCMz /= MTOT;

  for(a = 0; a < NA; a++)
    {
      for(i = 0; i < Oparams.parnum[a]; i++)
	{
	  rx[a][i] -= rint(RCMx/la);
	  ry[a][i] -= rint(RCMy/la);
	  rz[a][i] -= rint(RCMz/la);
	  rxS[a][i] = minimumImageI(rx[a][i]);
	  ryS[a][i] = minimumImageI(ry[a][i]);
	  rzS[a][i] = minimumImageI(rz[a][i]);
	}
    }
}

/* =========================== >>> initCoord <<< ============================*/
void initCoord(void)
{
  setToZero(Oparams.parnum[0], SAVE_LISTA, 
	    NULL);  /* Set to zero all the coordinates */
  
  setToZero(Oparams.parnum[1], SAVE_LISTB, 
	    NULL);  /* Set to zero all the coordinates */

  FCC(Oparams.m); 
  
  /* Put the baricenter of each molecule on a FCC lattice, and set 
     their orientations */  
  
  /* set the exact velocity of both atoms, considering the rotational motion 
     of the molecule, too. */
  //angvel(Oparams.parnum, Oparams.T, Oparams.m, Oparams.d); 
}

/* =========================== >>> usrInitBef <<< ========================== */
void usrInitBef(void)
{
  /* DESCRIPTION:
     This function is called before any other initialization, put here 
     yours, for example initialize accumulators ! 
     NOTE: You should supply parameters value in parameters file, so this 
           initilization are quite fictitiuos for parameters, anyway 
	   accumulators initialization is crucial */
  
  /* ===================== >>> INIT PARAMETERS <<< ======================== 
   All the values set here for Oparams structure are taken as defaults if you
   don't specify corresponding parameters in the parameters file  */
  int i, a, b;
  Dtrans = 0.0; /* DtransOld should become a field of OprogStatus */

  Vol = 1400.0;
  Vol1 = 0.0;
  Vol2 = 0.0;
  Vol1o1 = 0.0;
  Vol1o2 = 0.0;

  Oparams.T = 2.0;
  Oparams.P = 1.0;

  OprogStatus.avVol = 0.0;
  OprogStatus.tolVol = 0.0001;
  OprogStatus.tolVol1 = 0.1;
  OprogStatus.nebrTabFac = 150;
  OprogStatus.rNebrShell = 0.4;
  /* If 1 the program calculate of the corrisponding variable a mean from
     the begin of the run and not the instanteaneous value */
  OprogStatus.avnggr    = 0;
  OprogStatus.avngS     = 0;
  OprogStatus.avngPress = 0;
  OprogStatus.avngTemp  = 0;
  OprogStatus.bakStepsAscii = 0;
  OprogStatus.rateCheck = 10000;
  OprogStatus.mosseAccettate = 0;
  OprogStatus.mosseTentate = 0;
  OprogStatus.fstps = 1;	
  for (i=0; i<NA; ++i)
    {
      Oparams.m[i] = 1.0;
    }

  for(a=0; a<NA; ++a)
    {
      for(b=0; b<NA; ++b)
	{
	  Oparams.sigab[a][a] = 1.0;
	  Oparams.epsab[a][b] = 1.0;
	}
    }
  V = 0.0; /* potential energy */
  W = 0.0; /* virial function */

  srand((int)time(NULL));
  Oparams.rcut = 2.7; /* the cutoff is for each atoms 'a' in each molecule:
			 Oparams.rcut * Oparams.sigma[a] */
  /* ======================================================================= */
}


extern void calcConst(void);
extern double rcutab[NA][NA], rcutabSq[NA][NA], dvdr[NA][NA]; 
extern double sigabSq[NA][NA], epsab4[NA][NA], epsab24[NA][NA];

/* ========================= >>> sumup <<< ================================*/
void sumup(void)
{
  /* calcola i valori delle energie potenziali all'inizio della
     simulazione */
  int aa, bb, a, i, ai, b, j, bj, Nm, rxa, rya, rza, rxabI, ryabI, rzabI;
  double vab, wab, rabSq, srab2, srab6, srab12, rxab, ryab, rzab, la;
  double Wab[NA][NA], Vab[NA][NA], ncut[NA][NA], Vcab[NA][NA], vabCut;
  Nm = Oparams.parnum[0]+Oparams.parnum[1];

  calcConst();
  V = 0.0;
  W = 0.0;
  la = Oparams.lattice_a;

  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  ncut[a][b] = 0;
	  Vab[a][b] = 0;
	  Wab[a][b] = 0;
	}
    }
  for(a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  rxa = rx[a][i];
	  rya = ry[a][i];
	  rza = rz[a][i];

	  ai = i+Nm*a;
	  for(b = 0; b < NA; b++)
	    {
	      for (j = 0; j < Oparams.parnum[b]; j++ )
		{
		  bj = j+b*Nm;
		  if (bj > ai) 
		    {
		      rxabI = minimumImageI(rxa - rx[b][j]);
		      ryabI = minimumImageI(rya - ry[b][j]);
		      rzabI = minimumImageI(rza - rz[b][j]);
		      /* Moltiplica le coordinate intere per il passo 
			 del reticolo
			 in modo da ottenere le coordinate 
			 in unità ridotte */
		      rxab = la * ((double) rxabI);
		      ryab = la * ((double) ryabI);
		      rzab = la * ((double) rzabI);
      
		      rabSq = Sqr(rxab) + Sqr(ryab) + Sqr(rzab);
		      
		      if ( rabSq < rcutabSq[a][b] )
			/* 'rcut' is the cutoff for V */
			{
			  srab2   = sigabSq[a][b] / rabSq;
			  srab6   = srab2 * srab2 * srab2;
			  srab12  = Sqr(srab6);
			  
			  vab     = srab12 - srab6;
			  //vab     = vab -  dvdr[a][b] * (rab - rcutab[a][b]);
			  wab     = vab + srab12;
			  
			  Vab[a][b] = Vab[a][b] + vab;
			  //printf("(%d,%d)-(%d,%d):%.5f | ", a, i, b, j, vab);
			  /* total potential between all a-b atoms pairs */
			  Wab[a][b]   = Wab[a][b] + wab; 
			  /* Virial off-diagonal terms of atomic pressure 
			     tensor */
			  ++ncut[a][b];
			}
		    }
		}
	    }
	  
	}
    }

  /* CALCULATE SHIFTED POTENTIAL
     shifted potential, for each atoms within rcut 
     subtracts Vcut = V(rcut) 
     (see pag 145 A.T.) */
  
  for(aa = 0; aa < NA; aa++)
    {
      for(bb = 0; bb < NA; bb++) /* b >= a */
	{
	  srab2 = sigabSq[aa][bb] / rcutabSq[aa][bb];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;
	  Vcab[aa][bb] = Vab[aa][bb] - ncut[aa][bb] * vabCut;
	  /* ncut[a][b] is the number of atoms pairs a-b within 
	     rcutab[a][b] */ 
	}
    }

  /* Moltiplica per le costanti  */
  for (aa = 0; aa < NA; aa++)
    {
      for (bb = 0; bb < NA; bb++)
	{
	  V += Vcab[aa][bb] * epsab4[aa][bb];
	  W  += Wab[aa][bb]  * epsab24[aa][bb] / 3.0;
	}
    }

}

/* ========================== >>> AllocCubeR <<< ===========================*/
double*** AllocCubeR(int size1, int size2, int size3)
{
  double*** v;
  void* buffer;
  int k1, k2;
  
  /* alloca la memoria per il cubo */
  buffer = malloc(size1 * size2 * size3 * sizeof(double));

  v = (double***) malloc(size1 * sizeof(double**));
  for (k1 = 0; k1 < size2; k1++)
    v[k1] = (double**) malloc(size2 * sizeof(double*));

  for (k1 = 0; k1 < size1; k1++)
    {
      v[k1][0] = ((double*) buffer) + k1 * size2 * size3;
      for (k2 = 1; k2 < size2; k2++)
	{
	  v[k1][k2] = v[k1][k2-1] + size3;
	}
    }

  return v;
}

/* ========================== >>> AllocCubeR <<< ===========================*/
int*** AllocCubeI(int size1, int size2, int size3)
{
  int*** v;
  void* buffer;
  int k1, k2;
  
  /* alloca la memoria per il cubo */
  buffer = malloc(size1 * size2 * size3 * sizeof(int));

  v = (int***) malloc(size1 * sizeof(int**));
  for (k1 = 0; k1 < size2; k1++)
    v[k1] = (int**) malloc(size2 * sizeof(int*));

  for (k1 = 0; k1 < size1; k1++)
    {
      v[k1][0] = ((int*) buffer) + k1 * size2 * size3;
      for (k2 = 1; k2 < size2; k2++)
	{
	  v[k1][k2] = v[k1][k2-1] + size3;
	}
    }

  return v;
}

/* ========================= >>> buildEnergiesR <<< ======================= */
void buildEnergiesRad()
{
  int i, a, b; //csSq;
  double rabSq, srab2, srab6, srab12;
  double vabCut;
  double cutoff, cutoffSq;

  calcConst();


  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  /* MODIFICA 25/04/01:
	     controllare qui perché prima c'era un baco enorme sul calcolo
	     dei cutoff interi che dava delle energie potenziali differenti */ 
	  cutoff = (Oparams.sigab[a][b]*Oparams.rcut)/Oparams.lattice_a;
	  cubeSize[a][b] = (int)cutoff;
	  if (((double)cubeSize[a][b]) < cutoff )
	    cubeSize[a][b] += 1;
	  
	  cutoffSq = Sqr(cutoff);
	  //csSq[a][b] = Sqr(cubeSize[a][b]);
	  csSq[a][b] = (int) cutoffSq;

	  if (((double)csSq[a][b]) < cutoffSq )
	    csSq[a][b] += 1;
	  
	  printf("rcut*sigabSq: %.8f sigab: %.8f cutoff: %.6f\n",
		 Oparams.sigab[a][b]*Oparams.rcut, Oparams.sigab[a][b], 
		 cutoff);
	  
	  printf("la: %.10f csSq[%d][%d]:%d\n",Oparams.lattice_a, a, b, 
		 csSq[a][b]);

	  vabRadLT[a][b] = malloc(sizeof(double)*csSq[a][b]);
	  wabRadLT[a][b] = malloc(sizeof(double)*csSq[a][b]);
	}
    }


  for (a = 0; a < NA; a++)
    {
      for (b = 0; b < NA; b++)
	{
	  //csSq = Sqr(cubeSize[a][b]);
	  srab2 = sigabSq[a][b] / rcutabSq[a][b];
	  srab6 = srab2 * srab2 * srab2;
	  srab12 = srab6 * srab6;
	  vabCut = srab12 - srab6;

	  for (i = 0; i < csSq[a][b]; i++)
	    {
	      
	      rabSq = Sqr((double)Oparams.lattice_a) * ((double)i);
	     
	      if ( rabSq < rcutabSq[a][b] )/* 'rcut' is the cutoff for V */
		{

		  srab2   = sigabSq[a][b] / rabSq;
		  srab6   = srab2 * srab2 * srab2;
		  srab12  = Sqr(srab6);

		  //printf("(%d) 12-6:%.20f\n", i,srab12 - srab6);
		  vabRadLT[a][b][i] = epsab4[a][b]*(srab12 - srab6 - vabCut);
		  //printf("8--() DENTRO!!!!%.10f\n", vabRadLT[a][b][i]);
		  //vab[a][b][ix][iy][iz] -=     
		  // dvdr[a][b] * (rab - rcutab[a][b]);
		  //wabRadLT[a][b][i] = epsab24[a][b]*
		  // ((srab12 - srab6) + srab12)/3.0;
		}
	    }
	}
    }
  
}


/* ========================== >>> builgPMap <<< ========================= */
void buildMaps(void)
{
  int i, lM, Nm;
  int *buf;
  lM = Oparams.lattice_M;
  Nm = Oparams.parnum[0] + Oparams.parnum[1];
  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pmap = buf + lM;

  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pabs = buf + lM;

  buf = (int *) malloc((2*lM + 1)*sizeof(int));
  pSqr = buf + lM;

  pgetj = (int *) malloc((2*Nm)*sizeof(int));
  pgetb = (int *) malloc((2*Nm)*sizeof(int));

  for (i = -lM; i <= lM; i++)
    {
      pmap[i] = minimumImageI(i);
      pabs[i] = abs(i); 
      pSqr[i] = Sqr(i);
     
    }
  
  for(i = 0; i < NA*Nm; i++)
    {  
      if (i >= Nm)
	{
	  pgetb[i] = 1;
	  pgetj[i] = i - Nm;
	}
      else
	{
	  pgetb[i] = 0;
	  pgetj[i] = i;
	}
    }
}


/* ======================== >>> usrInitAft <<< ==============================*/
void usrInitAft(void)
{
  /* DESCRIPTION:
     This function is called after the parameters were read from disk, put
     here all initialization that depends upon such parameters, and call 
     all your function for initialization, like maps() in this case */

  int Nm, i, sct;
  COORD_TYPE* m;
  int a;

  //COORD_TYPE RCMx, RCMy, RCMz, Rx, Ry, Rz;

  /* initialize global varibales */
  pi = 2.0 * acos(0);
  
  mdMsg(ALL, NOSYS, "usrInitAft", "NOTICE", NULL,
	"Using Neighbour List method",
	NULL);
    
  /* [0][1] elements are read from parameter file */
  Oparams.sigab[1][0] = Oparams.sigab[0][1];
  Oparams.epsab[1][0] = Oparams.epsab[0][1];
  
  /* Allocate arrays of integers, you must supply the array length and the
     address of the pointer to the array */
  Nm = Oparams.parnum[0] + Oparams.parnum[1];

  sct = sizeof(COORD_TYPE);
  
  /* Initialize variables for neighbour list method */
  nebrTabMax = OprogStatus.nebrTabFac;
  nebrNow = 1;
  nebrTab = AllocMatI(NA*Nm, nebrTabMax); 
  nebrTabi = AllocMatI(NA*Nm, nebrTabMax); 
  nebrTaba = AllocMatI(NA*Nm, nebrTabMax); 
 
  nebrTabLen = malloc(sizeof(int)*Nm*NA);
  /* =============================================== */
  
  for (i=0; i < NA*Nm; i++)
    {
      nebrTabLen[i] = 0;
    }
  
  Oparams.lattice_a = cbrt(Vol) / Oparams.lattice_M;
  //  printf("lattice_a: %.5f Vol: %.5f\n", Oparams.lattice_a, Vol);
  /* Store the Center of Mass initial position for all particles */
  m = Oparams.m;

  /* The fields rxCMi, ... of OprogStatus must contain the centers of mass 
     positions, so wwe must initialize them! */  
  if (newSim == 1)
    {
      for ( a = 0; a < NA; a++)
	{
	  for(i = 0; i < Oparams.parnum[a]; i++)
	    {
	      
	      /* store the initial positions of particles */
	      OprogStatus.rxi[a][i] = rx[a][i];
	      OprogStatus.ryi[a][i] = ry[a][i];
	      OprogStatus.rzi[a][i] = rz[a][i];
	    }
	}
      
      loop(i, 1, NUMK) 
	{
	  OprogStatus.sumS[i] = 0.0;
	}
      
    }
  printf("Vol: %.15f\n", Vol);

  for (a = 0; a < NA; a++)
    {
      /* coordinate scalate alla prima scatola */
      rxS[a] = (int *)malloc(sizeof(int)*Oparams.parnum[a]);
      ryS[a] = (int *)malloc(sizeof(int)*Oparams.parnum[a]);
      rzS[a] = (int *)malloc(sizeof(int)*Oparams.parnum[a]);
    }

  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* coordinate scalate alla prima scatola */
	  //printf("(%d,%d,%d)\n", rx[a][i], ry[a][i], rz[a][i]);
	  rxS[a][i] = minimumImageI(rx[a][i]);
	  ryS[a][i] = minimumImageI(ry[a][i]);
	  rzS[a][i] = minimumImageI(rz[a][i]);
	}
    }
  /* calcola le energie potenziali e il viriale all'inizio della simulazione */
  sumup();
  
  /* calcola tutte le energie di interazione all'inizio */

  buildEnergiesRad();
  buildMaps();

}


/* ========================== >>> writeAllCor <<< ========================== */
void writeAllCor(FILE* fs)
{
  int i, a;
  double la;
  la = Oparams.lattice_a;
  for (a = 0; a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  /* salva le coordinate come numeri in virgola mobile */
	  fprintf(fs, "%.15G %.15G %.15G\n", la*((double)rx[a][i]), 
		  la*((double)ry[a][i]), la*((double)rz[a][i]));
	}
    }

}


/* ========================== >>> readAllCor <<< ========================== */
void readAllCor(FILE* fs)
{
  int i, a;
  double xx,yy,zz, la;
  
  la = Oparams.lattice_a;

  for (a = 0;  a < NA; a++)
    {
      for (i = 0; i < Oparams.parnum[a]; i++)
	{
	  if (fscanf(fs, "%lf %lf %lf\n", &xx, &yy, &zz) < 3)
	    {
	      mdPrintf(STD, "ERROR[pos] reading ascii file\n", NULL);
	      exit(-1);
	    }
	  /* converte le coordinate negli interi corrisponednti */
	  rx[a][i] = rint(xx/la);
	  ry[a][i] = rint(yy/la);
	  rz[a][i] = rint(zz/la);
	}
    }
}
            