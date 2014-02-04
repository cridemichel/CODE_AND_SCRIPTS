      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                  
      INTEGER *4 NFI,NACC                                                     
      PARAMETER(NMAXM=1000,NLISTMAX=100000)
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
      PARAMETER ( MAXK = 1000 )
      REAL*8      KVEC(MAXK), KAPPA
      COMMON /EWALD1/ KVEC,KAPPA  
      COMMON /BKSRANGES/ RANGECAA,RANGECBB,RANGECAB,RANGEQCAA,
     *      RANGEQCBB,RANGEQCAB
      COMMON /HARMRANGE/  RHARMBB,RHARMQBB, AMPHARMBB,EOFFSETBB,
     *    RHARMAB,RHARMQAB, AMPHARMAB,EOFFSETAB

