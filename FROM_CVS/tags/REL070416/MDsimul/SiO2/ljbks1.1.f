c
c  Test Function - Check Forces Energies
c
c  VERSIONE PER BINARY MIXTURES ---- DIFFERENT MASSES ALLOWED
c
c   Non V=Vlj+A+B*r/sigma bensi Vlj-V(2.5sigma)
c
c   Version for Free energy calculations
c
c   Switch between a BMLJ and a Charged BM  (with Ewald)
c
c
      INCLUDE 'ljbm.com'   
      open(unit=11,file='fort.101',status='unknown')
      open(unit=10,file='fort.100',status='unknown')
c      
c   read the potential parameters
c
      CALL READRESTART
      CALL INITAL     
      CALL INITALCHARGES 
c
      NACC=0
c
c
      DO NSTEP=1,NZAHL
       NACC=NACC+1
       NFI=NFI+1
c      CALL NELIST       
      
        CALL EWRFORCESCHARGES  
        CALL EWKFORCESCHARGES 
c
       CALL MOTION 
       IF (MOD(NSTEP,LPROT).EQ.0) CALL VALUES         
       IF (KDISK.AND.(MOD(NSTEP,LDISK)).EQ.0) CALL WRITERESTART
c
c  save also configuration if.... 
c
      IF (KTAPE.AND.(NFI.EQ.LTAPE)) THEN
         CALL WRITESTORY
         CALL WRITERESTART
      END IF 
      END DO                                 
c
C
      IF (KDISK) CALL WRITERESTART 
      STOP                                                                     
      END              
c                                                         
C     ******************************************************************       C
       SUBROUTINE INITAL                                  
       INCLUDE 'ljbks.com'  
       AVOGA=6.02205D+02
       PI=3.1415926535898d0
       BOLTK= 8.314D-03
       NLIST=NLISTMAX
c
c  This is for 33% 66% mixture
c
       MOLA=MOLS/3
       print *,' #A,#B,#(A+B)',MOLA,MOLS-MOLA,MOLS
c
       print *,'BOX SIZE=',bx,by,bz
       bmin=95001. 
       if (bx.Lt.bmin) bmin=bx     
       if (by.Lt.bmin) bmin=by 
       if (bz.Lt.bmin) bmin=bz 
       bmin=bmin/2.0d0      
       BOXHX = 0.5d0*BX 
       boxhy = 0.5d0*by 
       boxhz = 0.5d0*bz  
C                                                              
C ARGON EPS=0.99605728D0,SIGMA=0.341D0, MASSA=39.95
C
C
c
c  SPECIE # 1
c
c                 
       EPSON=1.0D0
       SIGMA=0.1D0
       WAITA= 28.085500D0 
       SIGAA=SIGMA
       EPSAA=EPSON
c
       print *,' Time Units=',dsqrt(wait*sigaa**2/48.0d0/epsaa)
c
       CL=2.50D0
       RSWITCHAA=CL*SIGMA
       RSWITCH2AA=RSWITCHAA*RSWITCHAA
c
c   EPSONR=EPSON RINORMALIZZATO
c   SIGMAR=SIGMA RINORMALIZZATA 
c
      EPSONR=EPSON
      SIGMAR=SIGMA
      CONREPAA= 4.00D0*EPSONR*SIGMAR**12
      CONATPAA= 4.00D0*EPSONR*SIGMAR**6
      CONREFAA=12.00D0*CONREPAA
      CONATFAA= 6.00D0*CONATPAA
      RANGEAA=   RSWITCHAA 
      RANGEQAA=RANGEAA*RANGEAA 
      AENEAA= -(CONREPAA/RSWITCHAA**12-CONATPAA/RSWITCHAA**6)
      print *,' Energy Shift AA=',AENEAA
c
c
c  SPECIE # 2  (BB)
c
c                 
       EPSON=0.5D0
       SIGMA=0.088D0
       WAITB= 15.999400D0
       EPSBB=EPSON
       SIGBB=SIGMA
c
       CL=2.50D0
       EPSONR=EPSON
       SIGMAR=SIGMA
       RSWITCHBB=CL*SIGMA
       RSWITCH2BB=RSWITCHBB*RSWITCHBB
c
c   EPSONR=EPSON RINORMALIZZATO
c   SIGMAR=SIGMA RINORMALIZZATA 
c
      CONREPBB= 4.00D0*EPSONR*SIGMAR**12
      CONATPBB= 4.00D0*EPSONR*SIGMAR**6
      CONREFBB=12.00D0*CONREPBB
      CONATFBB= 6.00D0*CONATPBB
      RANGEBB=   RSWITCHBB 
      RANGEQBB=RANGEBB*RANGEBB 
      AENEBB=-(CONREPBB/RSWITCHBB**12-CONATPBB/RSWITCHBB**6)
      print *,' Energy Shift BB=',AENEBB
c
c  AB interactions 
c
c                 
       EPSON=1.50D0
       SIGMA=0.080D0
       EPSONR=EPSON
       SIGMAR=SIGMA
       EPSAB=EPSON
       SIGAB=SIGMA
       CL=2.50D0
       RSWITCHAB=CL*SIGMA
       RSWITCH2AB=RSWITCHAB*RSWITCHAB
c
c   EPSONR=EPSON RINORMALIZZATO
c   SIGMAR=SIGMA RINORMALIZZATA 
c
      CONREPAB= 4.00D0*EPSONR*SIGMAR**12
      CONATPAB= 4.00D0*EPSONR*SIGMAR**6
      CONREFAB=12.00D0*CONREPAB
      CONATFAB= 6.00D0*CONATPAB
      RANGEAB=   RSWITCHAB 
      RANGEQAB=RANGEAB*RANGEAB 
      AENEAB=-(CONREPAB/RSWITCHAB**12-CONATPAB/RSWITCHAB**6)
      print *,' Energy Shift AB=',AENEAB
c
       range=max(rangeaa,rangeab,rangebb)   
c
       FMOLS=DFLOAT(MOLS)
       if (range.gt.bmin) range=bmin
       print *,' Using range =',range                                
       RANT  = 0.95D0*RANGE 
c era 1.1
       RANG2 = 1.150D0*RANGE 
       RANGE3= 1.00D0/RANGE**3
       RANG3H= 0.50D0*RANGE3 
       PERMIT= 0.25D0*(RANG2-RANGE)**2  
C                  
       RANTQ = RANT**2 
       RANGEQ= RANGE**2 
       RANG2Q= RANG2**2 
       VOLUME= BX*BY*BZ
C                                                        
      FACTMP(1)= 2.0d0/(3.0d0*BOLTK)     
      FACTMP(2)= 1.0d0/(3.0d0*BOLTK)                                           
      FACEKIA= 0.5d0*WAITA/FMOLS                                    
      FACEKIB= 0.5d0*WAITB/FMOLS
C                                                     
c
C  TO be controlled                  
c
      FACPRE= 2.0D+03*FMOLS/(3.0d0*AVOGA)                     
      FACWOR= 1.0D-03*AVOGA/FMOLS      
      FACVIR= 0.5d0/FMOLS            
c 
      FACDENA= WAITA*FMOLS/AVOGA
      FACDENB= WAITB*FMOLS/AVOGA                                                 
      FACGVR= DELTR*VOLUME/(2.0d0*PI*FMOLS**2)  
C                    
      RETURN
      END                                                       
C     ******************************************************************  
      SUBROUTINE READRESTART
      INCLUDE 'ljbks.com'  
c
c  read information on the type of run
c     
      NINP=1
      OPEN(UNIT=NINP,FILE='ljbks.inp',STATUS='OLD')  
      READ (NINP,*) KNEBG,KDISK,KTAPE,KTEMP,KPRES
      READ (NINP,*) LNEBG,LDISK,LTAPE,LDTEM,LDPRE
      READ (NINP,*) MOLS ,LCRTL,LPROT,NZAHL
      READ (NINP,*) ETEMP,ERHO ,EPRES
      READ (NINP,*) DELTA,TAUT ,TAUP ,COMP ,TIMA
      READ (NINP,*) BX,BY,BZ
      CLOSE(NINP)
c
c  read restart configuration
c
      IRR=0
      NDAT=1
      OPEN(UNIT=NDAT,FILE='ljbks.res',STATUS='OLD')
      if (IRR.NE.0) GOTO 777        
      READ (NDAT,*,IOSTAT=IRR)     NFI ,NACC , IMOLS, BX, BY, BZ
       if (IRR.NE.0) GOTO 777  
      READ (NDAT,*,IOSTAT=IRR)     AEGES,
     -                  ATEMP,APRES,ARHO,DELTAI,RHOP,EEGES,LIOST

      IF (IMOLS.NE.MOLS) PAUSE ' # MOLECOLE DIVERSO RESTART E INPUT'
      IF (DELTA.NE.DELTAI) PRINT *,' USO DELTA=',DELTA
      DO I=1,MOLS
       READ(NDAT,*,IOSTAT=IRR ) RN(I,1),RN(I,2),RN(I,3)
        if (IRR.NE.0) GOTO 777 
       DO K=1,3
         RNLAST(I,K)=RN(I,K)
       END DO
      END DO
      DO I=1,MOLS
       READ(NDAT,*,IOSTAT=IRR ) VE(I,1),VE(I,2),VE(I,3)
        if (IRR.NE.0) GOTO 777 
      END DO
      CLOSE(NDAT)
c
 777  IF (IRR.NE.0) THEN
         CALL CONFIG
          
      END IF  
c
c
c  check of periodic boundary conditions on the read data
c
      BOXHX=BX*0.5d0
      BOXHY=BY*0.5d0
      BOXHZ=BZ*0.5d0
c  
      DO 4800 I=1,MOLS                                                           
      NPERX =INT(ABS(RN(I,1))/BOXHX)                                          
      ANPERX=DSIGN((REAL(NPERX)*Bx),RN(I,1))                                  
      NPERy =INT(ABS(RN(I,2))/BOXHy)                                          
      ANPERy=DSIGN((REAL(NPERy)*By),RN(I,2))                                  
      NPERz =INT(ABS(RN(I,3))/BOXHZ)                                          
      ANPERz=DSIGN((REAL(NPERZ)*Bz),RN(I,3))                                  
c     if (abs(anperx)+abs(anpery)+abs(anperz).ne.0) then
c     print *, 'Fuori scatola in restart',i,nperx,npery,nperz
c     end if
C                                                                               
      RN(I,1)=RN(I,1)-ANPERX
      RN(I,2)=RN(I,2)-ANPERY                                                
      RN(I,3)=RN(I,3)-ANPERZ                                                
C                                                                               
 4800 CONTINUE  
C
C  it is a continuation of a already started simulation. Read other stuff
C     
      IF (KTAPE.AND.(.NOT.KNEBG)) THEN
      NDAT=1 
      OPEN(UNIT=NDAT,FILE='ljbks.cont',
     -            STATUS='OLD',form='unformatted')        
C
      READ (NDAT)    RSTART     
      READ (NDAT)
     -           VIRIA,VIRIRA,VIRIQA,EKINA,EKINRA,EKINQA,TROTA,
     -           EGESA,EGESRA,EGESQA,EPOTA,EPOTRA,EPOTQA,TCMAA,
     -           PRESA,PRESRA,PRESQA,TEMPA,TEMPRA,TEMPQA,EKCMA,DQM,
     -           RHOA ,RHORA ,RHOQA ,WORKA,WORKRA,WORKQA,FREEA,GKR,SQA,
     -           ELECA,TDPQA
      READ(ndat) base,npc,ltape 
      print *,' BASE NPC LTAPE',base,npc,ltape
      CLOSE(NDAT)
      END IF
C                       
C  Read the information on the number or records previosly written
C                                                        
c   knebg=t (if start conf)
c
      IF (.NOT.KTAPE) goto 333
      IF (KTAPE.AND.KNEBG)               THEN                    
       LIOST=0
       NFI=0
       PRINT *,' I WILL START TO SAVE FROM THIS RUN'
       print *,' INSERISCI BASE E  NPUNTI'
       read(5,*) base,npc
       call tempilog 

      ELSE
       PRINT *,' I CONTINUE TO SAVE FROM REC N.',LIOST,ltape
C
c calculate next time 
c
      call tempilog
      ncicli=liost/npc
c     print *,' N cicli completi=',ncicli
c      
      nexts=liost-ncicli*npc+1
c 
c     print *,' Prossimo stato (tra 1 e npc)',nexts
      ltape=ncicli*nss(npc)+nss(nexts)
c
c     print *,' Next time I will save when nfi=',ltape
      END IF 
333   CONTINUE
C                                                                               
      RETURN                                                                    
      END                                                                       
c            
C     ******************************************************************        
      SUBROUTINE WRITERESTART
      INCLUDE 'ljbks.com'  
c
c
c  write restart configuration
c
      NDAT=1
      OPEN(UNIT=NDAT,FILE='ljbks.res',STATUS='UNKNOWN')        
      WRITE (NDAT,*)     NFI ,NACC , MOLS, BX, BY, BZ
      WRITE (NDAT,*)     AEGES,
     -                  ATEMP,APRES,ARHO,DELTA,RHOP,EEGES,LIOST
c
      DO I=1,MOLS
       WRITE(NDAT,*) RN(I,1),RN(I,2),RN(I,3)
      END DO
      DO I=1,MOLS
       WRITE(NDAT,*) VE(I,1),VE(I,2),VE(I,3)
      END DO
      CLOSE(NDAT)
C
      OPEN(UNIT=NDAT,FILE='ljbks.cont',STATUS='UNKNOWN',
     -  FORM='UNFORMATTED')        
C
      WRITE (NDAT)    RSTART     
      WRITE (NDAT)
     -           VIRIA,VIRIRA,VIRIQA,EKINA,EKINRA,EKINQA,TROTA,
     -           EGESA,EGESRA,EGESQA,EPOTA,EPOTRA,EPOTQA,TCMAA,
     -           PRESA,PRESRA,PRESQA,TEMPA,TEMPRA,TEMPQA,EKCMA,DQM,
     -           RHOA ,RHORA ,RHOQA ,WORKA,WORKRA,WORKQA,FREEA,GKR,SQA,
     -           ELECA,TDPQA
      write(ndat) base,npc,ltape 
      CLOSE(NDAT)
C                
      RETURN                                                                    
      END                                                                       
c            
C     ******************************************************************        
      SUBROUTINE WRITESTORY
      INCLUDE 'ljbks.com'  
      character*10 nrecord
      character*1 recname(10)
      equivalence (nrecord, recname(1))
c
c  write configuration number liost
c
      liost=liost+1
      write(nrecord,'(i10)') liost
c
c  clean from empty spaces
c

 77    continue
       if (recname(1).eq.' ') then
          do k=2,10
            recname(k-1)=recname(k)
          end do
          recname(10)=' '
          goto 77
        end if
        print *,' Start writing File =','t5.'//nrecord

c
      NDAT=17
      OPEN(UNIT=NDAT,FILE='t043.'//nrecord,STATUS='UNKNOWN')        
      WRITE (NDAT,*)     NFI ,NACC , MOLS, BX, BY, BZ
      WRITE (NDAT,*)     AEGES,
     -                  ATEMP,APRES,ARHO,DELTA,RHOP,EPOTEN,LIOST
c
      DO I=1,MOLS
       WRITE(NDAT,*) RN(I,1),RN(I,2),RN(I,3)
      END DO
      DO I=1,MOLS
       WRITE(NDAT,*) VE(I,1),VE(I,2),VE(I,3)
      END DO

      CLOSE(NDAT)
      print *,' End writing' 
C
c calculate next time 
c
      ncicli=liost/npc
      nexts=liost-ncicli*npc+1
      ltape=ncicli*nss(npc)+nss(nexts)
c     print *,' Next time I will save when nfi=',ltape
      RETURN
      END                                                                       
c          
      SUBROUTINE NELIST
      INCLUDE 'ljbks.com'  
c
      DIMENSION RIJ(3),INDEX(nmaxm)                                              
      DATA INDEX /nmaxm*0/                                                       
      DO 10 I=1,MOLS                                                            
C                                                                               
      DIFFX=ABS(RN(I,1)-RNLAST(I,1))                                          
      DIFFY=ABS(RN(I,2)-RNLAST(I,2))                                          
      DIFFZ=ABS(RN(I,3)-RNLAST(I,3))                                          
C                                                                               
      IF (DIFFX.GT.BOXHx) DIFFX=DIFFX-BX                                        
      IF (DIFFY.GT.BOXHy) DIFFY=DIFFY-BY                                        
      IF (DIFFZ.GT.BOXHz) DIFFZ=DIFFZ-BZ                                        
C                                                                               
      DIFFSQ=DIFFX**2+DIFFY**2+DIFFZ**2                                         
C                                                                               
      IF((DIFFSQ.GT.PERMIT).OR.NACC.EQ.1) GO TO 15                                             
   10 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
C
 15   CONTINUE
C
c     print *,'Nelist MOL,DR2,NACC',i,diffsq,nacc           
c
      DO 20 I=1,MOLS                                                            
C                                                                               
      RNLAST(I,1)=RN(I,1)                                                     
      RNLAST(I,2)=RN(I,2)                                                     
      RNLAST(I,3)=RN(I,3)                                                     
C                                                                               
   20 CONTINUE    
c
c
      CALL NELISTAA
      CALL NELISTBB
      CALL NELISTAB
c
      RETURN
      END 
c
      SUBROUTINE NELISTAA
      INCLUDE 'ljbks.com'  
C
      DIMENSION RIJ(3),INDEX(nmaxm)                                              
      DATA INDEX /nmaxm*0/    
      INL=1                                                                     
C                                                                               
      DO 30 I=  1,MOLA-1                                                        
      DO 28 J=I+1,MOLA                                                          
C                                                                               
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      RIJ(1)=RIJ(1)-BX*REAL(INT(RIJ(1)/BOXHX))                                  
      RIJ(2)=RIJ(2)-By*REAL(INT(RIJ(2)/BOXHY))                                  
      RIJ(3)=RIJ(3)-Bz*REAL(INT(RIJ(3)/BOXHZ))                                  
C                                                                               
      ROO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2                                         
      IND=1                                                                     
      IF (ROO.GE.RANG2Q) IND=0                                                  
C                                                                               
      INDEX(J)=IND                                                              
   28 CONTINUE                                                                  
C                                                                               
      DO 29 J=I+1,MOLS                                                          
C                                                                               
      LISTAA(INL)=J                                                               
      INL=INL+INDEX(J)                                                          
C                                                                               
   29 CONTINUE                                                                  
C                                                                               
      LASTAA(I)=INL-1                                                             
C                                                                               
   30 CONTINUE                                                                  
C                                                                               
c     print *,'N interazioni AA ',INL
      IF (INL.LE.NLIST) RETURN                                                  
C                                                                               
      WRITE(6,*) ' ERROR'
C     
      PRINT *,' NLIST=',NLIST                                                                          
      STOP                                                                      
      END 
C     ******************************************************************        
      SUBROUTINE NELISTBB
      INCLUDE 'ljbks.com'  
C

      DIMENSION RIJ(3),INDEX(nmaxm)                                              
      DATA INDEX /nmaxm*0/                                                       
      INL=1                                                                     
C                                                                               
      DO 30 I=  MOLA+1,MOLS-1  
      DO 28 J=I+1,MOLS                                                          
C                                                                               
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      RIJ(1)=RIJ(1)-BX*REAL(INT(RIJ(1)/BOXHX))                                  
      RIJ(2)=RIJ(2)-By*REAL(INT(RIJ(2)/BOXHY))                                  
      RIJ(3)=RIJ(3)-Bz*REAL(INT(RIJ(3)/BOXHZ))                                  
C                                                                               
      ROO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2                                         
      IND=1                                                                     
      IF (ROO.GE.RANG2Q) IND=0                                                  
C                                                                               
      INDEX(J)=IND                                                              
   28 CONTINUE                                                                  
C                                                                               
      DO 29 J=I+1,MOLS                                                          
C                                                                               
      LISTBB(INL)=J                                                               
      INL=INL+INDEX(J)                                                          
C                                                                               
   29 CONTINUE                                                                  
C                                                                               
      LASTBB(I)=INL-1                                                             
C                                                                               
   30 CONTINUE                                                                  
C                                                                               
c     print *,'N interazioni BB ',INL
      IF (INL.LE.NLIST) RETURN                                                  
C                                                                               
      WRITE(6,*) 'ERROR'
C     
       PRINT *,' NLIST=',NLIST
       STOP
       END                                                                          

C     ******************************************************************        
      SUBROUTINE NELISTAB
      INCLUDE 'ljbks.com'  
C                                                                               
      DIMENSION RIJ(3),INDEX(nmaxm)                                              
      DATA INDEX /nmaxm*0/                                                       
C                                                                               
      INL=1                                                                     
C                                                                               
      DO 30 I=  1,MOLA                                                        
      DO 28 J=MOLA+1,MOLS                                                          
C                                                                               
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      RIJ(1)=RIJ(1)-BX*REAL(INT(RIJ(1)/BOXHX))                                  
      RIJ(2)=RIJ(2)-By*REAL(INT(RIJ(2)/BOXHY))                                  
      RIJ(3)=RIJ(3)-Bz*REAL(INT(RIJ(3)/BOXHZ))                                  
C                                                                               
      ROO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2                                         
      IND=1                                                                     
      IF (ROO.GE.RANG2Q) IND=0                                                  
C                                                                               
      INDEX(J)=IND                                                              
   28 CONTINUE                                                                  
C                                                                               
      DO 29 J=I+1,MOLS                                                          
C                                                                               
      LISTAB(INL)=J                                                               
      INL=INL+INDEX(J)                                                          
C                                                                               
   29 CONTINUE                                                                  
C                                                                               
      LASTAB(I)=INL-1                                                             
C                                                                               
   30 CONTINUE                                                                  
C                                                                               
c     print *,'N interazioni AB ',INL

      IF (INL.LE.NLIST) RETURN                                                  
C                                                                               
      WRITE(6,*) 'ERROR'
      PRINT *,' NLIST=',NLIST                   
      STOP    
      END 
C     ****************************************************************** 
      SUBROUTINE FORCES                                                 
      INCLUDE 'ljbks.com'  
C                                                              
      DIMENSION RIJ(3),BBR(3),VIRR(3),CCEL(3),FO(nmaxm,3)           
C                         
      EPOTEN=0.0D0                         
      VIRIAL=0.0D0                      
C                
      DO 50 J=1,3  
      DO 50 I=1,MOLS                                
      FO(I,J)  =0.0D0                                
   50 CONTINUE                                     
C  
c AA interactions
c                
c       print *,CONREPAA,CONATPAA,AENEAA,BENEAA,CONREFAA,CONATFAA,mola
c       pause
c 
      MFIRST=1            
      DO 300 I=1,MOLA-1                          
      IF(MFIRST.GT.LASTAA(I)) GO TO 300          
      DO 330 JM=MFIRST,LASTAA(I)
      J=LISTAA(JM)                                                        
C                                                                               
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx))                                         
      BBR(2)=By*DBLE(INT(RIJ(2)/BOXHy))                                         
      BBR(3)=Bz*DBLE(INT(RIJ(3)/BOXHz))                                         
c
C                                                                               
      RIJ(1)=RIJ(1)-BBR(1)                                                      
      RIJ(2)=RIJ(2)-BBR(2)                                                      
      RIJ(3)=RIJ(3)-BBR(3)                                                      
C                                                                               
      RO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2   
      IF (RO.GE.RANGEQAA) GO TO 330                                               

      RQ =1.0D0/RO                                           
      RIJP=DSQRT(RO)
C                                                                               
      VIRR(1)=0.0D0                                                             
      VIRR(2)=0.0D0                                                             
      VIRR(3)=0.0D0                                                             
C                    
c   lw LJ !
c                                                           
       TSN   =RQ**3 
c       POTDM =    TSN*(CONREPAA*TSN-CONATPAA)+AENEAA+BENEAA*RIJP
        POTDM =    TSN*(CONREPAA*TSN-CONATPAA)+AENEAA
c       CCEL3 = RQ*TSN*(CONREFAA*TSN-CONATFAA)-BENEAA/RIJP
        CCEL3 = RQ*TSN*(CONREFAA*TSN-CONATFAA) 
c
      EPOTEN=EPOTEN+POTDM                                                       

c       if (i.lt.10) write(98,98) i,j,potdm
c 98    format(i4,1x,i4,1x,e18.12)

      CCEL(1)=CCEL3*RIJ(1)                                                      
      CCEL(2)=CCEL3*RIJ(2)                                                      
      CCEL(3)=CCEL3*RIJ(3)                                                      
C                                                                               
      FO(I,1)=FO(I,1)+CCEL(1)                                                   
      FO(I,2)=FO(I,2)+CCEL(2)                                                   
      FO(I,3)=FO(I,3)+CCEL(3)                                                   
C                                                                               
      FO(J,1)=FO(J,1)-CCEL(1)                                                   
      FO(J,2)=FO(J,2)-CCEL(2)                                                   
      FO(J,3)=FO(J,3)-CCEL(3)                                                   
C                                                                               
      VIRIAL=VIRIAL-CCEL(1)*RIJ(1) 
     *             -CCEL(2)*RIJ(2)
     *             -CCEL(3)*RIJ(3)
C               
  340 CONTINUE                                                                  
  330 CONTINUE                                                                  
C                                                                               
       MFIRST=LASTAA(I)+1                                                          
  300 CONTINUE                                                                  
C                     
C             
c  ora BB
c   
C  
      MFIRST=1            
      DO 1300 I=MOLA+1,MOLS-1                          
      IF(MFIRST.GT.LASTBB(I)) GO TO 1300          
      DO 1330 JM=MFIRST,LASTBB(I)
      J=LISTBB(JM)                                                        
C                                                                               
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx))                                         
      BBR(2)=By*DBLE(INT(RIJ(2)/BOXHy))                                         
      BBR(3)=Bz*DBLE(INT(RIJ(3)/BOXHz))                                         
c
C                                                                               
      RIJ(1)=RIJ(1)-BBR(1)                                                      
      RIJ(2)=RIJ(2)-BBR(2)                                                      
      RIJ(3)=RIJ(3)-BBR(3)                                                      
C                                                                               
      RO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2   
      IF (RO.GE.RANGEQBB) GO TO 1330                                               

      RQ =1.0D0/RO                                           
      RIJP=DSQRT(RO)
C                                                                               
      VIRR(1)=0.0D0                                                             
      VIRR(2)=0.0D0                                                             
      VIRR(3)=0.0D0                                                             
C                    
c   lw LJ !
c                                                           
        TSN   =RQ**3 
c       POTDM =    TSN*(CONREPBB*TSN-CONATPBB)+AENEBB+BENEBB*RIJP
c       CCEL3 = RQ*TSN*(CONREFBB*TSN-CONATFBB)-BENEBB/RIJP
       POTDM =    TSN*(CONREPBB*TSN-CONATPBB)+AENEBB
       CCEL3 = RQ*TSN*(CONREFBB*TSN-CONATFBB)
c
       EPOTEN=EPOTEN+POTDM                                                      

c       if (i.lt.10) write(98,98) i,j,potdm
c 98    format(i4,1x,i4,1x,e18.12)

      CCEL(1)=CCEL3*RIJ(1)                                                      
      CCEL(2)=CCEL3*RIJ(2)                                                      
      CCEL(3)=CCEL3*RIJ(3)                                                      
C                                                                               
      FO(I,1)=FO(I,1)+CCEL(1)                                                   
      FO(I,2)=FO(I,2)+CCEL(2)                                                   
      FO(I,3)=FO(I,3)+CCEL(3)                                                   
C                                                                               
      FO(J,1)=FO(J,1)-CCEL(1)                                                   
      FO(J,2)=FO(J,2)-CCEL(2)                                                   
      FO(J,3)=FO(J,3)-CCEL(3)                                                   
C                                                                               
      VIRIAL=VIRIAL-CCEL(1)*RIJ(1) 
     *             -CCEL(2)*RIJ(2)
     *             -CCEL(3)*RIJ(3)
C               
1340  CONTINUE                                                                  
1330  CONTINUE                                                                  
C                                                                               
       MFIRST=LASTBB(I)+1
1300       CONTINUE       
c
c  ORA AB
c
      MFIRST=1                               
      DO 2300 I=1,MOLA                          
      IF(MFIRST.GT.LASTAB(I)) GO TO 2300          
      DO 2330 JM=MFIRST,LASTAB(I)                                                  
      J=LISTAB(JM)   
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx))                                         
      BBR(2)=BY*DBLE(INT(RIJ(2)/BOXHy))                                         
      BBR(3)=BZ*DBLE(INT(RIJ(3)/BOXHz))                                         
c
C                                                                               
      RIJ(1)=RIJ(1)-BBR(1)                                                      
      RIJ(2)=RIJ(2)-BBR(2)                                                      
      RIJ(3)=RIJ(3)-BBR(3)                                                      
C                                                                               
      RO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2   
      IF (RO.GE.RANGEQAB) GO TO 2330                                               

      RQ =1.0D0/RO                                           
      RIJP=DSQRT(RO)
C                                                                               
      VIRR(1)=0.0D0                                                             
      VIRR(2)=0.0D0                                                             
      VIRR(3)=0.0D0                                                             
C                    
c   lw LJ !
c                                                           
       TSN   =RQ**3 
c       POTDM =    TSN*(CONREPAB*TSN-CONATPAB)+AENEAB+BENEAB*RIJP
c       CCEL3 = RQ*TSN*(CONREFAB*TSN-CONATFAB)-BENEAB/RIJP
        POTDM =    TSN*(CONREPAB*TSN-CONATPAB)+AENEAB
        CCEL3 = RQ*TSN*(CONREFAB*TSN-CONATFAB)
c
      EPOTEN=EPOTEN+POTDM                                                       

c       if (i.lt.10) write(98,98) i,j,potdm
c 98    format(i4,1x,i4,1x,e18.12)

      CCEL(1)=CCEL3*RIJ(1)                                                      
      CCEL(2)=CCEL3*RIJ(2)                                                      
      CCEL(3)=CCEL3*RIJ(3)                                                      
C                                                                               
      FO(I,1)=FO(I,1)+CCEL(1)                                                   
      FO(I,2)=FO(I,2)+CCEL(2)                                                   
      FO(I,3)=FO(I,3)+CCEL(3)                                                   
C                                                                               
      FO(J,1)=FO(J,1)-CCEL(1)                                                   
      FO(J,2)=FO(J,2)-CCEL(2)                                                   
      FO(J,3)=FO(J,3)-CCEL(3)                                                   
C                                                                               
      VIRIAL=VIRIAL-CCEL(1)*RIJ(1) 
     *             -CCEL(2)*RIJ(2)
     *             -CCEL(3)*RIJ(3)
C               
 2340  CONTINUE                                                                  
 2330  CONTINUE                                                                  
C                                                                               
       MFIRST=LASTAB(I)+1                                                          
 2300  CONTINUE       
c
      DO 180 J=1,3                                                              
      DO 180 I=1,MOLS                                                           
C                                                                               
      F(I,J)= FO(I,J)                                               
C                                                                               
  180 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C     ******************************************************************        
      SUBROUTINE MOTION                                                         
      INCLUDE 'ljbks.com'  
c
      ALFA    =1.0D0                                                            
      PARKEA   =0.0D0                                                            
      PARKEB   =0.0D0
C                                                                               
C     *********************** Tempraturkontrolle ***********************        
C                                                                               
      IF (KTEMP) CALL VSCALE(ALFA)                                              
C      
      DO 450 J=1,3                                                              
      DO  I=1,MOLA                                                           
      PARKEA=PARKEA+VE(I,J)**2
      END DO
      DO  I=MOLA+1,MOLS                                                           
      PARKEB=PARKEB+VE(I,J)**2
      END DO
  450 CONTINUE              
C                    
c Sum of all forces, energies and virial
c
      DO J=1,3
        DO I=1,MOLS
          F(I,J)=FCHARGES(I,J)+FEWK(I,J)
        END DO
       END DO
c
       EPOTEN=EPOTENCHARGES+EPOTENEWK
       VIRIAL=VIRIALCHARGES
c
      DO 460 J=1,3                                                              
      DO I=1,MOLA                                                           
      VE(I,J)=(VE(I,J)+F(I,J)*DELTA/WAITA)*ALFA                        
      RM(I,J)= RN(I,J)                                                      
      RN(I,J)= RM(I,J)+VE(I,J)*DELTA                                      
      END DO
      DO I=MOLA+1,MOLS                                                           
      VE(I,J)=(VE(I,J)+F(I,J)*DELTA/WAITB)*ALFA                        
      RM(I,J)= RN(I,J)                                                      
      RN(I,J)= RM(I,J)+VE(I,J)*DELTA                                      
      END DO 
  460 CONTINUE                                                                  
C                                                                               
C            
      DO 480 I=1,MOLS                                                           
      NPERX =INT(ABS(RN(I,1))/BOXHX)                                          
      ANPERX=DSIGN((REAL(NPERX)*Bx),RN(I,1))                                  
      NPERy =INT(ABS(RN(I,2))/BOXHy)                                          
      ANPERy=DSIGN((REAL(NPERy)*By),RN(I,2))                                  
      NPERz =INT(ABS(RN(I,3))/BOXHZ)                                          
      ANPERz=DSIGN((REAL(NPERZ)*Bz),RN(I,3))                                  
C                                                                               
      RN(I,1)=RN(I,1)-ANPERX
      RN(I,2)=RN(I,2)-ANPERY                                                
      RN(I,3)=RN(I,3)-ANPERZ                                                
C                                                                               
  480 CONTINUE                                                                  
C                EVERY 10 STEPS CHECK DRIFT IN MOMENTUM                         
c      IF (MOD(NACC,10).eq.1) CALL DRIFT                                         
      RETURN                                                                    
      END                                                                       
C     ******************************************************************        
c
C     ******************************************************************        
      SUBROUTINE VSCALE(ALFA)                                                   
      INCLUDE 'ljbks.com'  

      TEMP=0.0D0                                                                
C                                                                               
      DO 1 J=1,3                                                                
      DO I=1,MOLA                                                             
      TEMP=TEMP+FACEKIA*VE(I,J)**2                                          
      END DO
      DO I=MOLA+1,MOLS                                                             
      TEMP=TEMP+FACEKIB*VE(I,J)**2                                          
      END DO
    1 CONTINUE                                                                  
C                                                                               
      TEMP=FACTMP(1)*TEMP
C                      
      xxx=1.0+DELTA*(ATEMP/TEMP-1.0)/TAUT
      if (xxx.le.0) then
        alfa=.1
        print *,' ALFA FIXED TO 0.1'
       else
      ALFA=DSQRT(1.0+DELTA*(ATEMP/TEMP-1.0)/TAUT)                               
c      print *,alfa,' questa T=',temp,' prevista=',atemp
      end if 
C                                                                               
      RETURN                                                                    
      END                                                                       
C     ******************************************************************        
      SUBROUTINE VALUES                                                         
      INCLUDE 'ljbks.com'  
C                                                                              
      FLPROT=REAL(LPROT)                                                        
      FNACC =REAL(NACC)                                                         
C                                                                               
C     ************ Berechnung der potentiellen Energie *****************        
C                                                                               
      EPOT=EPOTEN/FMOLS-FACPOT                                            
C                                                                               
      EPOTA =EPOTA +EPOT                                                        
      EPOTRA=EPOTRA+EPOT                                                        
      EPOTQA=EPOTQA+EPOT**2                                                     
C                                                                               
      EPOTP =EPOTA /FNACC                                                       
      EPOTRP=EPOTRA/FLPROT                                                      
      EPOTQP=EPOTQA/FNACC                                                       
      EPOTF =DSQRT(ABS(EPOTQP-EPOTP**2))                                        
C                                                                               
C     **** Berechnung der kinetischen Energie des Molekuels ************        
C                                                                               
      EKIN=FACEKIA*PARKEA+FACEKIB*PARKEB
C                                                                               
      EKINA =EKINA +EKIN                                                        
      EKINRA=EKINRA+EKIN                                                        
      EKINQA=EKINQA+EKIN**2                                                     
C                                                                               
      EKINP =EKINA /FNACC                                                       
      EKINRP=EKINRA/FLPROT                                                      
      EKINQP=EKINQA/FNACC                                                       
      EKINF =DSQRT(ABS(EKINQP-EKINP**2))                                        
C                                                                               
C     ***************** Berechnung der Tempratur ***********************        
C                                                                               
      TEMP=FACTMP(1)*EKIN                                                       
C                       
c    QUI TEMP IN KELVIN OK   PRINT *,'TEMP=',temp 
c
      TEMPA =TEMPA +TEMP                                                        
      TEMPRA=TEMPRA+TEMP                                                        
      TEMPQA=TEMPQA+TEMP**2                                                     
C                                                                               
      TEMPP =TEMPA /FNACC                                                       
      TEMPRP=TEMPRA/FLPROT                                                      
      TEMPQP=TEMPQA/FNACC                                                       
      TEMPF =DSQRT(ABS(TEMPQP-TEMPP**2))                                        
C                                                                               
C     ******************** VIRIALE PER BM LJ WK
C                          (NO CORRECTIONS)                       
c                                           
      VIRI=FACVIR*VIRIAL                                           
C                                                                               
      VIRIA =VIRIA +VIRI                                                        
      VIRIRA=VIRIRA+VIRI                                                        
      VIRIQA=VIRIQA+VIRI**2                                                     
C                                                                               
      VIRIP =VIRIA /FNACC                                                       
      VIRIRP=VIRIRA/FLPROT                                                      
      VIRIQP=VIRIQA/FNACC                                                       
      VIRIF =DSQRT(ABS(VIRIQP-VIRIP**2))                                        
C                                                                               
c  pressure according to Eq. 2.2.9 in Hansen and Mc Donalds
c
c                    PIG (units of Kj/mol divide nm^3)
c
       PIG=MOLS*BOLTK*TEMP/VOLUME
c
c  VIRIAL e' in Kj/mol ed e' ri-rj*Fij sum i sum j>i
c 
       PINT=-1.0d0/3.0d0*VIRIAL/VOLUME
c
      PRES=PIG+PINT
c

C                                                                               
      PRESA =PRESA +PRES                                                        
      PRESRA=PRESRA+PRES                                                        
      PRESQA=PRESQA+PRES**2                                                     
C                                                                               
      PRESP =PRESA /FNACC                                                       
      PRESRP=PRESRA/FLPROT                                                      
      PRESQP=PRESQA/FNACC                                                       
      PRESF =DSQRT(ABS(PRESQP-PRESP**2)) 
c

      write(10,*) nacc,real(pig+pint),real(pig),real(pint)        
      write(11,*) nacc,real(pres*volume)                                       
c      write(102,*)  nacc,real(PRESP),fnacc
C                  
C     ****************** Berechnung der Volumenarbeit ******************        
C                                                                               
      WORK=FACWOR*APRES*VOLUME                                                  
C                                                                               
      WORKA =WORKA +WORK                                                        
      WORKRA=WORKRA+WORK                                                        
      WORKQA=WORKQA+WORK**2                                                     
C                                                                               
      WORKP =WORKA /FNACC                                                       
      WORKRP=WORKRA/FLPROT                                                      
      WORKQP=WORKQA/FNACC                                                       
      WORKF =DSQRT(ABS(WORKQP-WORKP**2))                                        
C                                                                               
C     **********   Energia totale per molecola al passo NACC (non mediata dunque)
C                                                                               
      EGES=EKIN+EPOT   
      IF (KPRES) EGES=EKIN+EPOT+WORK                                            
C                                                                               
      EGESA =EGESA +EGES                                                        
      EGESRA=EGESRA+EGES                                                        
      EGESQA=EGESQA+EGES**2                                                     
C                                                                               
      EGESP =EGESA /FNACC                                                       
      EGESRP=EGESRA/FLPROT                                                      
      EGESQP=EGESQA/FNACC                                                       
      EGESF =DSQRT(ABS(EGESQP-EGESP**2))                                        
c
      write(98,99) nacc,eges,ekin,epot
 99   format(i15,3(2x,e20.12)) 
c      write(99,*) nacc, EGESP,FNACC             
c
C     ***********************   Densita'
C                                                                               
      RHO=FACDEN/VOLUME                                                         
C                                                                               
      RHOA =RHOA +RHO                                                           
      RHORA=RHORA+RHO                                                           
      RHOQA=RHOQA+RHO**2                                                        
C                                                                               
      RHOP =RHOA /FNACC                                                         
      RHORP=RHORA/FLPROT                                                        
      RHOQP=RHOQA/FNACC                                                         
      RHOF =DSQRT(ABS(RHOQP-RHOP**2))                                           
C                                                                               
C     *******************   Energia Libera
C                                                                               
      FREE=EGES                                                                 
C                                                                               
      IF (TEMP.EQ.0) THEN                                                       
      print *,' ATTENTION T=0 K'                                                
      T=1.                                                                      
      END IF                                                                    
      FREEH=FREE/(BOLTK*TEMP)                                                   
      IF (DABS(FREEH).GT.170.0D0) FREEH=DSIGN(170.0D0,FREEH)                    
C                                                                               
      FREEA= FREEA+EXP(-FREEH)                                                  
      FREEP=-BOLTK*TEMPP*LOG(FREEA/FNACC)                                       
C                                                                               
C     ********************** Controllo Pressione
C                                                                               
      IF (KPRES) CALL DSCALE                                              
C                                                                               
      IF (MOD(NACC,LCRTL).EQ.0)                        THEN                     
C                                                                               
      IF ((ABS(DTEM).GE.ABS(ETEMP-ATEMP)).AND.(LDTEM.NE.0)) LDTEM=0             
      IF ((ABS(DPRE).GE.ABS(EPRES-APRES)).AND.(LDPRE.NE.0)) LDPRE=0             
C                                                                               
      IF (KTEMP.AND.(LDTEM.NE.0)) ATEMP=ATEMP+DTEM                              
      IF (KPRES.AND.(LDPRE.NE.0)) APRES=APRES+DPRE                              
C                                                                               
      IF (KTEMP.AND.(LDTEM.EQ.0)) ATEMP=ETEMP                                   
      IF (KPRES.AND.(LDPRE.EQ.0)) APRES=EPRES                                   
C                                                                               
      END IF                                                                    
C                                                                               
      IF (MOD(NACC,LPROT).NE.0) RETURN                                          
C                                                                               
      PRESRA=0.0D0                                                              
      VIRIRA=0.0D0                                                              
      TEMPRA=0.0D0                                                              
      EPOTRA=0.0D0                                                              
      EKINRA=0.0D0                                                              
      EGESRA=0.0D0                                                              
      RHORA =0.0D0                                                              
      WORKRA=0.0D0                                                              
C                                                                               
      RETURN                                                                    
      END                                                                       
C     ******************************************************************        
C
      SUBROUTINE DSCALE                                                   
      INCLUDE 'ljbks.com'  
c
      BETA=(1.0-COMP*DELTA*(APRES-PRES)/TAUP)**(1.0/3.0)    
C                                                          
      BX =BETA*BX 
      by=beta*by 
      bz=beta*bz 
      BOXHx=0.50*Bx
      boxhy=.5*by  
      boxhz=.5*bz  
      WWx  =2.00*PI /BX
      WWy  =2.00*PI /By 
      WWz  =2.00*PI /Bz 
CC                       
      VOLUME=Bx*by*bz
      FACGVR=DELTR*VOLUME/(2.0*PI*FMOLS**2)
      FACPOT=8.0*PI*RANGE3*EPSON*FMOLS*SIGMA**6/(3.0*VOLUME)   
C
      DO 2 J=1,3 
      DO 2 I=1,MOLS
      RN(I,J)=RN(I,J)*BETA                                  
    2 CONTINUE 
C                                                                               
      RETURN                                                                    
      END                                                                       

      SUBROUTINE STRUCT
      INCLUDE 'ljbks.com'  

C                                                                               
C     ********************* QUADRATIC DISPLACEMENT *********************        
C                                                                               
      DQM=0.0d0    
      DO 6 J=1,MOLS                                                             
      DIF=ABS(RN(J,1)-RSTART(J,1))                           
      IF (DIF.GT.BOXHx) DIF=DIF-BX                                              
      DQM=DQM+DIF**2                                                            
      DIF=ABS(RN(J,2)-RSTART(J,2))                           
      IF (DIF.GT.BOXHy) DIF=DIF-BY                                              
      DQM=DQM+DIF**2                                                            
      DIF=ABS(RN(J,3)-RSTART(J,3))                           
      IF (DIF.GT.BOXHz) DIF=DIF-BZ                                              
      DQM=DQM+DIF**2                                                            
    6 CONTINUE                                                                  
C                                                                               
      DQM=DQM/FMOLS                                                             
C                                                                               
      RETURN                                                                    
      END              


                                                                 
C     ******************************************************************        
      SUBROUTINE DRIFT                                                          
      INCLUDE 'ljbks.com'  
c
      dimension p(3)                                                            
      p(1)=0.0d0                                                                
      p(2)=0.0d0                                                                
      p(3)=0.0d0                                                                
      do k=1,mola                                                            
       do l=1,3                                                              
        p(l)=p(l)+waita*ve(k,l)                                               
       end do
      end do
      do k=mola+1,mols                                                            
       do l=1,3                                                              
        p(l)=p(l)+waitb*ve(k,l)                                               
       end do
      end do

       totm= mola*waita+(mols-mola)*waitb

      do 34 l=1,3
34    p(l)=p(l)/totm                                               
      do 30 k=1,mols  
      do 30 l=1,3 
      ve(k,l)=ve(k,l)-p(l)                                                  
30    continue              
      return               
      end                                                                      
C     ****************************************************************** 
      SUBROUTINE CONFIG
      INCLUDE 'ljbks.com'  
      CALL INITAL
c
c  FIx the box size from the requested density
c
      totm= mola*waita+(mols-mola)*waitb
      arho=erho
      volume=bx*by*bz
      print *,' arho=',arho
      print *,'avoga',avoga
      VOLUME  =totm/(ARHO*AVOGA)
      vfat=(volume/bx/by/bz)**0.333333333333333333333333333
      bx=bx*vfat
      by=by*vfat
      bz=bz*vfat
      volume=bx*by*bz
c
c  calculate number density
c
      adn=float(mols)/volume
      print *,'Resulting Number density (units of sigmaAA)',adn
      print *, 'V per mol =',1.0d0/adn
      print *, 'Side per mol ',(1.0d0/adn)**0.333333333333333
      print *,' BX=',BX
      pause 
c
c  random insertion of particles (to avoid melting a crystall !!!)
c
      ih=113343454
      do i=1,mols
        do j=1,3
 666     xxx=rand(ih)
         if (xxx.eq.0) goto 666
         rn(i,j)=(xxx-0.5d0)*bx
         ve(i,j)=(rand(ih)-0.5d0)*1e-10
        end do
      end do
c
c   
c        

      DELTA=DELTA*1e-4
      print *,' DELTA=',DELTA
      KTAPE=.FALSE.
      KNEBG=.FALSE.
      KTEMP=.TRUE.
      TAUT=DELTA*3.5
      NSTEP=10000      
      LPROT=1
c
      RETURN
      END
       

      subroutine tempilog 
      INCLUDE 'ljbks.com'
      nss(1)=1
      LTAPE=1
      do kk=2,npc      
      nss(kk)=int(base**(kk-1))
      if (nss(kk).le.nss(kk-1)) nss(kk)=nss(kk-1)+1
      end do
c
      print *,' I will save conf at the following times '
      do kk=1,npc
       print *,kk,nss(kk),base**(kk-1)
      end do
      pause
      return
      end

C 
c    ******************************************************************       C
c
c  Data for the potential  v(r)=q*q/r+amp*exp(-r*vvv)+ccc/r**8
c
c    cut off 1/2 A (r-rcut)^2
c
c (information on masses are in INITAL)
c
c   Charges are in a.u.  
c                  vvv in nm
c                  amp in ?   --- 
c
       SUBROUTINE INITALCHARGES                                  
       INCLUDE 'ljbks.com'  
C
c                      to go to  q*q/r in kJ/mol with r in nm
        QFACTOR=11.7870740
c           
c               da eV a kJ/mol
        EFACTOR=96.484683
        AVN=6.02205d23
c              1ev in coulomb
        EVINC=1.60219d-19
c                      to go to  q*q/r in kJ/mol with r in nm
        qpieps0=4.0d0*3.14159*8.854188d-12
        Qfactor=(1.0d0/qpieps0*evinc*evinc/(1.0d-9)/1000.0d0*avn)**0.5
        print *,' Qfactor',qfactor        
        
C                                                              
c  SPECIE # 1  Si-Si
c 
        QBKA=2.4d0
        QBKA=QBKA*QFACTOR
        QBKAA=QBKA*QBKA
        AMPAA=0.0d0
        VVVAA=0.0d0
        CCCAA=0.0d0
        RANGECAA=0.0d0
        RANGEQCAA=RANGECAA*RANGECAA

c
c  SPECIE # 2  (BB) Oxygen-Oxygen
c
         QBKB=-1.2d0
         QBKB=QBKB*QFACTOR
         QBKBB=QBKB*QBKB
c        print *,' ampbb da peter s code (J) ',2.224814e-16
c             questa e' in eV
         AMPBB=1388.7730d0
         AMPBB=AMPBB*EFACTOR
c        print *,'Amp OO kJ/mol',ampbb
c        print *,'Amp in J',ampbb/avn*1000.0d0
c   
        VVVBB=27.60d0
c        print *,' C O-O da peters code', -2.80350e-23
        CCCBB= - 175.0000d-6
        CCCBB=CCCBB*EFACTOR
c        print *,' CCC in kJ/mol nm6',  CCCBB
c        print *,' CCC in J nm6',  CCCBB/avn*1000.0d0
        FOEXPBB= AMPBB*VVVBB  
        FORR6BB= 6.0d0*CCCBB
        RANGECBB=0.7d0
        RANGEQCBB=RANGECBB*RANGECBB
        IF (RANGECBB.GT.BX*0.5d0) pause' ERROR RANGECBB'
         RHARMBB=0.1439d0
         RHARMQBB=RHARMBB*RHARMBB
         AMPHARMBB=1.0d13


c
c  AB interactions Si-Oxygen
c
         QBKAB=QBKA*QBKB
c        AMPAB=2.884201e-15
c        print *,' ampbb da peter s code',ampab
c             questa e' in eV
        AMPAB=18003.7572d0
        AMPAB=AMPAB*EFACTOR
c        print *,'Amp SiO kJ/mol',ampab
c        print *,'Amp AB in J',ampab/avn*1000.0d0
         VVVAB= 48.7318d0
c        print *,' C Si-O da peters code', -2.13928e-23
        CCCAB= -133.5381d-6
        CCCAB=CCCAB*EFACTOR
c        print *,' CCC in kJ/mol nm6',  CCCAB
c        print *,' CCC in J nm6',  CCCAB/avn*1000.0d0
         FOEXPAB= AMPAB*VVVAB
         FORR6AB= 6.0d0*CCCAB
         RANGECAB=0.7d0
         RANGEQCAB=RANGECAB*RANGECAB
         IF (RANGECAB.GT.BX*0.5d0) pause' ERROR RANGECAB'
         RHARMAB=0.11936d0
         RHARMQAB=RHARMAB*RHARMAB
         AMPHARMAB=1.0d17
c
c  Offset per harmonic cat
c
        TSN   =  1.0d0/RHARMQAB**3 
        ABEXP =  DEXP(-VVVAB*RHARMAB)
        EOFFSETAB= qbkab/RHARMAB+ AMPAB*ABEXP   + TSN*CCCAB
C
        TSN   =  1.0d0/RHARMQBB**3 
        BBEXP =  DEXP(-VVVBB*RHARMBB)
        EOFFSETBB= qbkbb/RHARMBB + AMPBB*BBEXP   + TSN*CCCBB
        print *, 'OFFSET AB BB', EOFFSETAB,EOFFSETBB
c
c  Plot potential
c
       open(unit=33,file='potSiSi.dat',status='unknown')
       open(unit=34,file='potOO.dat',status='unknown')
       open(unit=35,file='potOSi.dat',status='unknown')
       do i=1,1000
         r=0.001*i
           write(33,*) r, qbkaa/r+ampaa*dexp(-vvvaa*r)+cccaa/r**6
          if (r.lt.RHARMBB) then
           write(34,*) r, AMPHARMBB*(r-RHARMBB)**2/2.0+eoffsetbb
           else
           write(34,*) r, qbkbb/r+ampbb*dexp(-vvvbb*r)+cccbb/r**6
          end if
          if (r.lt.RHARMAB) then
           write(35,*) r, AMPHARMAB*(r-RHARMAB)**2/2.0+eoffsetab
           else 
           write(35,*) r, qbkab/r+ampab*dexp(-vvvab*r)+cccab/r**6      
          end if
 
c           write(36,99) r, qbkaa/r,ampaa*dexp(-vvvaa*r),cccaa/r**6
c           write(37,99) r, qbkbb/r,ampbb*dexp(-vvvbb*r),cccbb/r**6
c           write(38,99) r, qbkab/r,ampab*dexp(-vvvab*r),cccab/r**6      
        end do
        close(33)
        close(34)
        close(35)
 99     format(4(e20.10,1x))
c  
        print *,' Calling SETUPEWALD'
        KAPPA = 5.0D0/BX
        print *,' Using KAPPA=',KAPPA
        CALL SETUPEWALD                   
C                    
      RETURN
      END                                                       
c
c
c  This calculate the real space part of the Ewald sum
c  the exp and the LJ terms
c
C     ****************************************************************** 
      SUBROUTINE EWRFORCESCHARGES    
      INCLUDE 'ljbks.com'  
C                                                              
      DIMENSION RIJ(3),BBR(3),CCEL(3),FO(nmaxm,3)           
C                         
      PI=4.0d0*DATAN(1.0d0)
      TOTPI=2.0d0/DSQRT(PI)
c      RHARMQBB=0.d0
      RHARMQAB=0.d0 
c
      EPOTENCHARGES=0.0D0                         
      VIRIALCHARGES=0.0D0                      
C                
      DO 50 J=1,3  
      DO 50 I=1,MOLS                                
      FO(I,J)  =0.0D0       
   50 CONTINUE                                     
C  
c AA interactions - ALL OF THEM !
c 
c
      DO 300 I=1,MOLA-1   
      DO 300 J=I+1,MOLA
C                                            
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C 
      BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx))                        
      BBR(2)=By*DBLE(INT(RIJ(2)/BOXHy))  
      BBR(3)=Bz*DBLE(INT(RIJ(3)/BOXHz))                                       
      RIJ(1)=RIJ(1)-BBR(1)                
      RIJ(2)=RIJ(2)-BBR(2)
      RIJ(3)=RIJ(3)-BBR(3)
      RO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2   

      RQINV =1.0D0/RO                                           
      RIJP=DSQRT(RO)
      RIJINV= 1.0D0/RIJP
C         
        ALPHAR= KAPPA*RIJP             
        ERC = DERFC(ALPHAR) 
        POTDM =  QBKAA*ERC/RIJP
        CCEL3 =  QBKAA*
     * (TOTPI*ALPHAR*DEXP(-ALPHAR*ALPHAR)+ERC)*RQINV*RIJINV
c
      EPOTENCHARGES=EPOTENCHARGES+POTDM

      CCEL(1)=CCEL3*RIJ(1)
      CCEL(2)=CCEL3*RIJ(2) 
      CCEL(3)=CCEL3*RIJ(3)
C 
      FO(I,1)=FO(I,1)+CCEL(1)
      FO(I,2)=FO(I,2)+CCEL(2)
      FO(I,3)=FO(I,3)+CCEL(3)
C                                                       
      FO(J,1)=FO(J,1)-CCEL(1)
      FO(J,2)=FO(J,2)-CCEL(2)
      FO(J,3)=FO(J,3)-CCEL(3)
C                                 
      VIRIALCHARGES=VIRIALCHARGES-CCEL(1)*RIJ(1) 
     *                           -CCEL(2)*RIJ(2)
     *                           -CCEL(3)*RIJ(3)
C               
  340 CONTINUE                                            
  330 CONTINUE                                                               
C                                                      
  300 CONTINUE                                                
 301  continue
C                     
c  ora BB
c   
C  
      DO 1300 I=MOLA+1,MOLS-1                          
      DO 1330 J=I+1,MOLS
C                                                                               
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHX))                                         
      BBR(2)=By*DBLE(INT(RIJ(2)/BOXHY))                                         
      BBR(3)=Bz*DBLE(INT(RIJ(3)/BOXHZ))                                         
c
C                                                                               
      RIJ(1)=RIJ(1)-BBR(1)                                                      
      RIJ(2)=RIJ(2)-BBR(2)                                                      
      RIJ(3)=RIJ(3)-BBR(3)                                                      
C                                                                               
      RO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2   


      RQINV =1.0D0/RO                                           
      RIJP=DSQRT(RO)
      RIJINV= 1.0D0/RIJP
C         
        ALPHAR= KAPPA*RIJP             
        ERC = DERFC(ALPHAR) 
        POTDM =  QBKBB*ERC/RIJP
        CCEL3 =  QBKBB*
     * (TOTPI*ALPHAR*DEXP(-ALPHAR*ALPHAR)+ERC)*RQINV*RIJINV 
c
c  exp and r-9 parts of the potential
c 
      IF (RO.LE.RANGEQCBB) THEN  

       IF (RO.GT.RHARMQBB) THEN
        TSN   =  RQINV**3 
        BBEXP =  DEXP(-VVVBB*RIJP)
        POTDM =  POTDM  + AMPBB*BBEXP   + TSN*CCCBB
        CCEL3 =  CCEL3  + FOEXPBB*BBEXP*RIJINV 
     *                  + FORR6BB*TSN*RQINV
        ELSE
         POTDM =  POTDM  + 
     *      0.5d0* AMPHARMBB*( RIJP - RHARMBB)**2 + EOFFSETBB
            
         CCEL3 =  CCEL3  - AMPHARMBB*(RIJP - RHARMBB)*RIJINV
       END IF
      END IF
c
       EPOTENCHARGES=EPOTENCHARGES+POTDM                   
C
      CCEL(1)=CCEL3*RIJ(1)                                                      
      CCEL(2)=CCEL3*RIJ(2)                                                      
      CCEL(3)=CCEL3*RIJ(3)                                                      
C                                                                               
      FO(I,1)=FO(I,1)+CCEL(1)                                                   
      FO(I,2)=FO(I,2)+CCEL(2)                                                   
      FO(I,3)=FO(I,3)+CCEL(3)                                                   
C                                                                               
      FO(J,1)=FO(J,1)-CCEL(1)                                                   
      FO(J,2)=FO(J,2)-CCEL(2)                                                   
      FO(J,3)=FO(J,3)-CCEL(3)                                                   
C                                                                               
      VIRIALCHARGES=VIRIALCHARGES-CCEL(1)*RIJ(1) 
     *                           -CCEL(2)*RIJ(2)
     *                           -CCEL(3)*RIJ(3)
C               
1340  CONTINUE                                                                  
1330  CONTINUE         
1300  CONTINUE       

c
c  ORA AB
c
      DO 2300 I=1,MOLA                          
      DO 2330 J=MOLA+1,MOLS 
      RIJ(1)=RN(I,1)-RN(J,1)                                                
      RIJ(2)=RN(I,2)-RN(J,2)                                                
      RIJ(3)=RN(I,3)-RN(J,3)                                                
C                                                                               
      BBR(1)=BX*DBLE(INT(RIJ(1)/BOXHx))                                         
      BBR(2)=BY*DBLE(INT(RIJ(2)/BOXHy))                                         
      BBR(3)=BZ*DBLE(INT(RIJ(3)/BOXHz))                                         
c
C                                                                               
      RIJ(1)=RIJ(1)-BBR(1)                                                      
      RIJ(2)=RIJ(2)-BBR(2)                                                      
      RIJ(3)=RIJ(3)-BBR(3)                                                      
C                                                                               
      RO=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2   
      RQINV =1.0D0/RO                                           
      RIJP=DSQRT(RO)
      RIJINV= 1.0D0/RIJP
C         
        ALPHAR= KAPPA*RIJP             
        ERC = DERFC(ALPHAR) 
        POTDM =   QBKAB*ERC/RIJP
        CCEL3 =  QBKAB*
     * (TOTPI*ALPHAR*DEXP(-ALPHAR*ALPHAR)+ERC)*RQINV*RIJINV 
c
c  exp and r-9 parts of the potential
c 
      IF (RO.LE.RANGEQCAB) THEN  
       IF (RO.GE.RHARMQAB) THEN
        TSN   =  RQINV**3 
        ABEXP =  DEXP(-VVVAB*RIJP)
        POTDM =  POTDM + AMPAB*ABEXP + TSN*CCCAB
        CCEL3 =  CCEL3 + FOEXPAB*ABEXP*RIJINV + FORR6AB*TSN*RQINV
        ELSE
         POTDM =  POTDM  + 
     *      0.5d0* AMPHARMAB*( RIJP - RHARMAB)**2 + EOFFSETAB
         CCEL3 =  CCEL3  - AMPHARMAB*( RIJP - RHARMAB)/RIJP
       END IF
      END IF
c
      EPOTENCHARGES=EPOTENCHARGES+POTDM 
c
      CCEL(1)=CCEL3*RIJ(1)                                                      
      CCEL(2)=CCEL3*RIJ(2)                                                      
      CCEL(3)=CCEL3*RIJ(3)                                                      
C                                                                               
      FO(I,1)=FO(I,1)+CCEL(1)                                                   
      FO(I,2)=FO(I,2)+CCEL(2) 
      FO(I,3)=FO(I,3)+CCEL(3)                                                   
C                                                                               
      FO(J,1)=FO(J,1)-CCEL(1)                                                   
      FO(J,2)=FO(J,2)-CCEL(2)                                                   
      FO(J,3)=FO(J,3)-CCEL(3)                                                   
C                                                                               
      VIRIALCHARGES=VIRIALCHARGES-CCEL(1)*RIJ(1) 
     *                           -CCEL(2)*RIJ(2)
     *                           -CCEL(3)*RIJ(3)
C               
 2340  CONTINUE                                                                  
 2330  CONTINUE                                                                  
C                                                                               
 2300  CONTINUE       
c
      DO 180 J=1,3                                                              
      DO 180 I=1,MOLS                                                           
C                     
      FCHARGES(I,J)= FO(I,J)                                               
C                                                                               
  180 CONTINUE                                                                  
c
      print *,' EPOTEN IN EWR',EPOTENCHARGES
c
      RETURN                                                                    
      END    

        SUBROUTINE SETUPEWALD
        INCLUDE 'ljbks.com'  

C    *******************************************************************
C    ** ROUTINE TO SET UP THE WAVE-VECTORS FOR THE EWALD SUM.         **
C    **                                                               **
C    ** THE WAVEVECTORS MUST FIT INTO A BOX OF UNIT LENGTH.           **
C    ** IN THIS EXAMPLE WE ALLOW A MAXIMUM OF 1000 WAVEVECTORS.       **
C    *******************************************************************
C  
C           K IN UNITS OF TWOPI/BX
C
        INTEGER     TOTK
        PARAMETER ( KMAX = 5, KSQMAX = 27  )

C    *******************************************************************

        TWOPI=2.0d0*4.0d0*DATAN(1.0d0) 
        B = 1.0D0 / 4.0D0 / KAPPA / KAPPA
        TWOPIOVERBOX=TWOPI/BX
C    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **


        TOTK = 0
        DO 100 KX = 0, KMAX
           RKX = TWOPIOVERBOX * KX 
           DO 99 KY = -KMAX, KMAX
              RKY = TWOPIOVERBOX * KY 
              DO 98 KZ = -KMAX, KMAX
                 RKZ = TWOPIOVERBOX * KZ 
                 KSQ = KX * KX + KY * KY + KZ * KZ
                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN
                    TOTK = TOTK + 1
                    IF ( TOTK .GT. MAXK ) STOP 'KVEC IS TOO SMALL'
                    RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ
                    KVEC(TOTK) = TWOPI * DEXP(-B*RKSQ)/RKSQ /BX**3
                 ENDIF
98            CONTINUE
99         CONTINUE

100     CONTINUE

        WRITE( *, ' ( '' EWALD SUM SETUP COMPLETE ''     ) ' )
        WRITE( *, ' ( '' NUMBER OF WAVEVECTORS IS '', I5 ) ' ) TOTK

        RETURN
        END

      SUBROUTINE EWKFORCESCHARGES    
      INCLUDE 'ljbks.com'  


C    *******************************************************************
C    ** CALCULATES K-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
C    **                                                               **
C    ** THE SELF TERM IS SUBTRACTED.                                  **
C    ** IN ONE COORDINATE DIRECTION (X), SYMMETRY IS USED TO REDUCE   **
C    ** THE SUM TO INCLUDE ONLY POSITIVE K-VECTORS.                   **
C    ** THE NEGATIVE VECTORS IN THIS DIRECTION ARE INCLUDED BY USE    **
C    ** OF THE MULTIPLICATIVE VARIABLE 'FACTOR'.                      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF IONS                    **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS OF IONS                 **
C    ** REAL    Z(N)                IONIC CHARGES                     **
C    ** REAL    VK                  K-SPACE POTENTIAL ENERGY          **
C    ** REAL    VKS                 SELF PART OF K-SPACE SUM          **
C    *******************************************************************

        INTEGER     TOTK

        INTEGER     KMAX, KX, KY, KZ, I, KSQMAX, KSQ
        PARAMETER ( KMAX = 5, KSQMAX = 27 )

        COMPLEX     EIKX(1:NMAXM, 0:KMAX)
        COMPLEX     EIKY(1:NMAXM, -KMAX:KMAX)
        COMPLEX     EIKZ(1:NMAXM, -KMAX:KMAX)
        COMPLEX     EIKR(NMAXM), SUM, SUMA, SUMB
        DIMENSION   FO(nmaxm,3)           
C
       PI=4.0d0*DATAN(1.0d0)
       TWOPI=2.0d0*PI
       TWOPIOVERBOX=TWOPI/BX
       RSQPI=1.0d0/DSQRT(PI)
C
      DO 50 J=1,3  
      DO 50 I=1,MOLS                                
      FO(I,J)  =0.0D0                                
   50 CONTINUE   
C
C
C    *******************************************************************
C
C    ** CONSTRUCT EXP(IK.R) FOR ALL IONS AND K-VECTORS **
C    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY **
C
        DO 10 I = 1, MOLS

           EIKX(I, 0) = (1.0, 0.0)
           EIKY(I, 0) = (1.0, 0.0)
           EIKZ(I, 0) = (1.0, 0.0)

           EIKX(I, 1) = CMPLX ( COS ( TWOPIOVERBOX * RN(I,1) ) ,
     :                          SIN ( TWOPIOVERBOX * RN(I,1) ) )
           EIKY(I, 1) = CMPLX ( COS ( TWOPIOVERBOX * RN(I,2) ) ,
     :                          SIN ( TWOPIOVERBOX * RN(I,2) ) )
           EIKZ(I, 1) = CMPLX ( COS ( TWOPIOVERBOX * RN(I,3) ) ,
     :                          SIN ( TWOPIOVERBOX * RN(I,3) ) )

           EIKY(I, -1) = CONJG ( EIKY(I, 1) )
           EIKZ(I, -1) = CONJG ( EIKZ(I, 1) )

10      CONTINUE

C    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **

        DO 12 KX = 2, KMAX
           DO 11 I = 1, MOLS
              EIKX(I, KX) = EIKX(I, KX-1) * EIKX(I, 1)
11         CONTINUE
12      CONTINUE

        DO 14 KY = 2, KMAX
           DO 13 I = 1, MOLS
              EIKY(I,  KY) = EIKY(I, KY-1) * EIKY(I, 1)
              EIKY(I, -KY) = CONJG ( EIKY(I, KY) )
13         CONTINUE
14      CONTINUE
c
        DO 16 KZ = 2, KMAX
           DO 15 I = 1, MOLS
              EIKZ(I,  KZ) = EIKZ(I, KZ-1) * EIKZ(I, 1)
              EIKZ(I, -KZ) = CONJG ( EIKZ(I, KZ) )
15         CONTINUE
16      CONTINUE

C    ** SUM OVER ALL VECTORS **
        print *,qbka,qbkb
C
        VD   = 0.0D0
        TOTK = 0
        DO 24 KX = 0, KMAX
c                             FACTOR tiene conto della riflession k-> -k 
           IF ( KX .EQ. 0 ) THEN
              FACTOR = 1.0D0  
           ELSE
              FACTOR = 2.0D0
           ENDIF
           DO 23 KY = -KMAX, KMAX
              DO 22 KZ = -KMAX, KMAX
                 KSQ = KX * KX + KY * KY + KZ * KZ
                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN
                    TOTK = TOTK + 1
                    SUMA  = (0.0, 0.0)
                    DO  I = 1, MOLA
                       EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ)
                       SUMA     = SUMA + EIKR(I)
                    END DO
                    SUMA=SUMA*QBKA
                    SUMB = (0.0, 0.0)
                    DO  I = MOLA+1, MOLS
                       EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ)
                       SUMB     = SUMB +  EIKR(I)
                    END DO
                    SUMB=SUMB*QBKB
                    SUM = SUMA+SUMB
                    VD = VD + FACTOR * KVEC(TOTK) * CONJG ( SUM ) * SUM

c
c Here forces  (see notes on NoteBook)
c
         DO  I = 1, MOLA 
           QFORCE= 
     c        FACTOR*KVEC(TOTK)*QBKA*AIMAG(EIKR(I)*CONJG(SUM)) 
               
              FO(I,1)=FO(I,1)+QFORCE* TWOPIOVERBOX* KX 
              FO(I,2)=FO(I,2)+QFORCE* TWOPIOVERBOX* KY 
              FO(I,3)=FO(I,3)+QFORCE* TWOPIOVERBOX* KZ 
          END DO
         DO  I = MOLA+1,MOLS 
           QFORCE= 2.0d0*
     c        FACTOR*KVEC(TOTK)*QBKB*AIMAG( EIKR(I)*CONJG(SUM)) 
              FO(I,1)=FO(I,1)+QFORCE* TWOPIOVERBOX* KX 
              FO(I,2)=FO(I,2)+QFORCE* TWOPIOVERBOX* KY 
              FO(I,3)=FO(I,3)+QFORCE* TWOPIOVERBOX* KZ 
          END DO    

                 ENDIF
                
22            CONTINUE
23         CONTINUE
24      CONTINUE

C    ** CALCULATES SELF PART OF K-SPACE SUM **
C       VS= SUM Q*Q
C
        VS= MOLA*QBKAA + (MOLS-MOLA)*QBKBB
        VS = RSQPI * KAPPA * VS

C    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **

        VK = VD - VS

c        write(88,*) rn(mols,1),vk,-fo(mols,1)
c
        print *,' EPOTEN EWK ',VK
c
c  Da decidere dove sia meglio sommare le F sulle cariche !
c
      DO 180 J=1,3                                                              
      DO 180 I=1,MOLS                                                           
                     
      FEWK(I,J)= FO(I,J)                                               
                                                                               
  180 CONTINUE                                                                  
c
       EPOTENEWK=VK
c
        RETURN
        END


 




