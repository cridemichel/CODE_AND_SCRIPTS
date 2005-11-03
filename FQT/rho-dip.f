c
c   calcola rho(l,m,k) per un fissato modulo - scrive i rho in files
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     *  (MOLS=512, LMAXX=2,KMAX=50,MAXK=KMAX)
      DIMENSION TH(MOLS),PHI(MOLS),CAI(MOLS)
      COMPLEX*16 RHO(0:LMAXX,-LMAXX:LMAXX,1:KMAX)
      COMPLEX*16 IMAG,RCE,RCEK,DL,DLCONJ
      DIMENSION KX(MAXK),KY(MAXK),KZ(MAXK)  
      COMPLEX*16 DXYZ(MOLS,0:LMAXX,-LMAXX:LMAXX),
     *      R(0:LMAXX,-LMAXX:LMAXX,-LMAXX:LMAXX,KMAX)    
      DIMENSION QROT(3,3) 
      character*2 label2
      character*3 label3,label
      character*80 nomefileconf
      external system
      common /sistema/  VROT(MOLS,3,3),CM(MOLS,3),BX
      common /AAA/ label3
      common /BBB/label2
c      
         open(unit=99,file='comando.sh',status='unknown')
         write (99,'(a)') 'mkdir RHOTMP'
         close(99)
        call system ('sh comando.sh')    
c
       pi=4.0d0*datan(1.0d0)
       twopi=2.0d0*pi  
       imag=(0.0d0,1.0d0) 

c
c read q-vectors with same modulus (file produced with program makeqrand)
c
         print *,' inserisci modulo '
         read(5,*) kmod
c        print *,' input lmin and lmax'
c        read(5,*) lmin,lmax
c
c         kmod=2
         lmin=0
         lmax=2
        write(label3,'(i3)') kmod
        call empty
        print *,' leggo dal file ','qvector.'//label3
        qmod=0.0d0
        open(unit=1,
     *    file='../../FQT/QVECTOR/qvector.'//label3,
     *             status='unknown')
        do i=1,maxk
         read(1,*,iostat=irr) kx(i),ky(i),kz(i)
         qmod=dsqrt(1.0d0*(kx(i)**2+ky(i)**2+kz(i)**2))+qmod
         if (irr.ne.0) goto 55
        end do
 55     nkv=i-1
        print *,' I am going to study # k-vectors=',nkv
        qmod=qmod/float(nkv)
        print *,' average modulus (units of twopi/box)',qmod
        close(1)
c
c
c
        do 20 k=1,nkv
c
c      calculate rotation that brings K along Z 
c
c        print *,k,KX(K), KY(K) , KZ(K)

         QNORM=DSQRT(1.0d0*(KX(K)*KX(K)+ KY(K)*KY(K)+ KZ(K)*KZ(K)))

c     Calcolo gli angoli di Eulero

         CAIRR=0.0d0

         THRR=dacos((1.0d0*kz(k))/qnorm)
         if (thrr.lt.0) thrr=-thrr
         if (thrr.gt.pi)thrr=thrr-pi 
c
c         if ((THRR.EQ.0).OR.(THRR.EQ.PI)) then
c
         if (dsin(thrr).eq.0) then
            PHIRR=0.0d0 
         else
            senphi=1.0d0*ky(k)/qnorm/dsin(THRR)
            if (senphi.ge.0.0d0) then
               xdum=(1.0d0*kx(k))/qnorm/dsin(THRR)
               if (xdum.ge.1.0d0) then
                  print *,' fisso phirr a zero poiche xdum=',xdum
                  phirr=0.0d0
               else
                  if (xdum.le.-1.0d0) then
                    print *,' fisso phirr a PI poiche xdum=',xdum
                    phirr=pi
                  else
                     PHIRR=dacos(xdum)
                  end if
               end if 
            else
               PHIRR=twopi-dacos((1.0d0*kx(k))/qnorm/dsin(THRR))
            endif   
         end if
c     
        do l=lmin,lmax
         do m=-l,l
          do n=-l,l
           call dlmn(l,m,n,thrr,res)
           dl=res*cdexp(-imag*m*phirr)  ! *cdexp(-imag*n*cairr)
           r(l,m,n,k)=dl
           end do
         end do
        end do   
 20     continue
c
c
c
          open(unit=8,file='listaconf',status='old')
          nconf=0          
c
 11       read(8,'(a)',end=22) nomefileconf
         
c
c  read 1 configuration  -- restituisce rm e vrot(i,3,3)
c
         nconf=nconf+1
c
         CALL READCONF(nomefileconf,ti)
         q0=twopi/bx
c        print *,' q0=',q0,mols
c
c  VROT(I,3,J)= 
c  VROT(I,2,J)= 
c  VROT(I,1,J)= 
c  
      DO 50 I=1,MOLS
      VMD1=0.0D0
      VMD2=0.0D0
      VMD3=0.0d0
      DO J=1,3
C
      VMD1=VMD1+ VROT(I,1,J)* VROT(I,1,J)
      VMD2=VMD2+ VROT(I,2,J)* VROT(I,2,J)      
      VMD3=VMD3+ VROT(I,3,J)* VROT(I,3,J)
C
      END DO
c
      VMD1=DSQRT( VMD1)
      VMD2=DSQRT( VMD2)
      VMD3=DSQRT( VMD3)
c
      if((dabs(vmd1-1.0d0).gt.1d-6).or.(dabs(vmd2-1.0d0).gt.1d-6).or.
     x      (dabs(vmd3-1.0d0).gt.1d-6) ) then
        print *,' MODULI VETTORI non unitari ', i,VMD1,vmd2,vmd3 
        stop
         end if
c
   50 CONTINUE
c
c
c    now....  euler angles for each molecule  phi,th,cai
c
       DO I=1,MOLS

        CTH   =VROT(I,3,3)
        IF (CTH.EQ.-1.0) THEN 
          TH(I) = 2.0*DACOS(0.0d0)
        ELSE IF (CTH.EQ.1.0) THEN
          TH(I) = 0.0
        ELSE 
          TH(I) = DACOS(CTH)
        END IF
c

        CCAI  =-VROT(I,1,3)/(1.0d0-CTH*CTH)**0.5   
c
        IF (CCAI.EQ.-1.0) THEN 
          CAI(I) = 2.0*DACOS(0.0d0)
        ELSE IF (CCAI.EQ.1.0) THEN
          CAI(I) = 0.0
        ELSE 
          CAI(I)= DACOS(CCAI)
        END IF
        IF (VROT(I,2,3).LT.0) CAI(I)=-CAI(I)
        IF (CAI(I).LT.0) CAI(I)=CAI(I)+TWOPI
c        
c  check  cx=sinth*sinphi
c         cy=-sinth*cosphi
c 
c         print *,VROT(I,2,3), sin(th(i))*sin(cai(i))
c         print *,VROT(I,1,3), -sin(th(i))*cos(cai(i))
c
        CPHI  =VROT(I,3,1)/(1.0d0-CTH*CTH)**0.5 
        IF (CPHI.EQ.-1.0) THEN 
          PHI(I) = 2.0*DACOS(0.0d0)
        ELSE IF (CPHI.EQ.1.0) THEN
          PHI(I) = 0.0
        ELSE 
          PHI(I)= DACOS(CPHI)
        END IF
        IF (VROT(I,3,2).LT.0) PHI(I)=-PHI(I)
        IF (PHI(I).LT.0) PHI(I)=PHI(I)+TWOPI 
c        
c  check  az=sinpsi*sinth
c         bz=cospsi*sinth
c 
c        print *,VROT(I,3,2),sin(th(i))*sin(phi(i))
c        print *,VROT(I,3,1),sin(th(i))*cos(phi(i))

c        print *, I, real(PHI(I)), real(th(i)), real(CAI(I))
       END DO
c
c
c     calculate all spherical harmonics related to the XYZ system  
c     DXYZ
c
       do i=1,mols
        do l=lmin,lmax
         do m=-l,l
           call dlmn(l,m,0,th(i),res)
           dl=res*cdexp(-imag*m*phi(i))
           DXYZ(i,l,m)=dl
          end do
        end do
       end do
c
        DO K=1,NKV
c
c   cleaning
c
        do l=lmin,lmax
         do m=-l,l
           rho(l,m,k)=(0.0d0,0.0d0)    
          end do
         end do
c  
c     print *,' starting',kx(k),ky(k),kz(k)

      do i=1,mols

       rcek=cdexp(imag*q0*(kx(k)*cm(i,1)+ky(k)*cm(i,2)+kz(k)*cm(i,3)))    

        do l=lmin,lmax
         do m=-l,l
c
c  calculate the rotated generalized harmonics
c
       dlconj=0.0d0   
        
       do mp=-l,l
         dlconj=dlconj+conjg(dxyz(i,l,mp))*r(l,mp,m,k)
         end do
         rho(l,m,k)= rho(l,m,k)+ rcek*dlconj
         end do
       end do
       end do
c
c   phase factor - according to Schilling's fax
c
        do l=lmin,lmax
         do m=-l,l
           rho(l,m,k)= rho(l,m,k)*dsqrt(1.0d0*
     *                      (2.0d0*l+1.0d0))*imag**l   
         end do
        end do
c
        end do   ! complete loop over k for the selected configuration
c        
c
c   save data in files named according to the harmonics 
c
c      print *,'nconf', nconf 
      if (nconf.eq.1) then
        write(label3,'(i3)') kmod
        call empty
        label=label3
        do l=lmin,lmax
         do m=-l,l
          ifn=l*10+m
          write(label2,'(i2)') ifn
          call empty2
          open(unit=10+ifn,
     *   file='RHOTMP/ro.'//label2//'.k='//label,
c                                         form='unformatted',
     *          status='unknown')
          print *,' Opened','RHOTMP/ro.'//label2//'.k='//label,l,m
          end do
        end do
       end if
c
        do l=lmin,lmax
         do m=-l,l
          ifn=l*10+m
           write(10+ifn,*) ti,nkv,(rho(l,m,ks),ks=1,nkv)
c           print *, nfi,nkt
c           print *,(rho(l,m,k),k=1,nkv)
          end do
        end do

       goto 11

c      here... go to next conf .....
c
c  close files
c
 22     continue
        do l=lmin,lmax
         do m=-l,l
          ifn=l*10+m
          close(unit=100+ifn)
         end do
        end do

       stop
       end


      subroutine dlmn(l,m,n,th,res)
      implicit double precision(a-h,o-z)
      dimension fact(0:10)
      logical first
      save first,fact
      data first /.TRUE./
      if (first) then
        fact(0)=1.0d0
         do k=1,10
          fact(k)=fact(k-1)*float(k)
         end do 
        first=.false.
      end if
c
      if (th.eq.0.0d0) then
         res=0.0d0
         if (m.eq.n) res=1.0d0
         return
      end if
c
       amp= dsqrt(fact(l+m)*fact(l-m)*fact(l+n)*fact(l-n))
c
       sth=dsin(th/2.0d0)
       cth=dcos(th/2.0d0)
       res=0.0d0
c        
       do 100 k=0,2*l
       i1=l+m-k
       if (i1.lt.0) return 
       i2=l-n-k
       if (i2.lt.0) return 
       i3=k-m+n
       if (i3.lt.0) goto 100
       res=res+amp*(-1.0d0)**k*cth**(2.0d0*l+m-n-2.0d0*k)*
     *           sth**(2.0d0*k-m+n)/(
     *  fact(l+m-k)*fact(l-n-k)*fact(k)*fact(k-m+n) )
 100   continue
       return
       end       
c
c
       subroutine empty 
       character*1 a(3)      
       common /AAA/ a
       do kk=1,3
         if (a(kk).eq.' ') a(kk)='0'
       end do  
       return            
       end
c
       subroutine empty2 
       character*1 a(2)      
       common /BBB/ a
       do kk=1,2
         if (a(kk).eq.' ') a(kk)='0'
       end do  
       return            
       end
c
       SUBROUTINE READCONF(nomefile,ti)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER(MOLS=512)
       character*80 nomefile
       double precision ti
       common /sistema/  VROT(MOLS,3,3),CM(MOLS,3),BX
       character*7 parnum
       character*5 time
       character*80 dummy
c
c        print *,' REading conf',nomefile
c
       open(unit=1,file=nomefile,status='old')
       read(1,10) parnum,nellips 
 10    format(a7,i5)
c       print *,' Number of ellipsoids=',nellips
        if (nellips.ne.mols) then
             print *,' Set mols to ',nellips
             stop
        end if
        do i=1,2
          read(1,*) dummy
        end do
        read(1,11) time,ti
 11     format(a5,D20.12)
        do i=1,14
          read(1,*) dummy
c        print *,dummy
        end do
        do i=1,nellips
            read(1,*) cm(i,1),cm(i,2),cm(i,3),
     x          vrot(i,1,1),vrot(i,1,2),vrot(i,1,3),
     x          vrot(i,2,1),vrot(i,2,2),vrot(i,2,3),
     x          vrot(i,3,1),vrot(i,3,2),vrot(i,3,3)
         end do
         do i=1,nellips
            read(1,*) dummy  ! velocities
         end do
         read(1,*) bx
c        print *,' BX=',bx
         close(1)
c
        return
        end
