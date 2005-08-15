c
c   calcola sq(k,l,m,lp,mp) per un fissato modulo - scrive i risultati 
c                 in files

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     *   (MOLS=512,LMAXX=1,KMAX=50,MAXK=KMAX,NPMAX=4000)
      DIMENSION TIME1(NPMAX)
c      DIMENSION SQLMN(0:LMAXX,0:LMAXX,0:LMAXX,1:KMAX)
      COMPLEX*16 RHO(0:LMAXX,-LMAXX:LMAXX,-LMAXX:LMAXX,1:KMAX)
      COMPLEX*16 IMAG,RCE,RCEK
      DIMENSION KX(MAXK),KY(MAXK),KZ(MAXK)  
      COMPLEX*16  R(0:LMAXX,-LMAXX:LMAXX,-LMAXX:LMAXX)    
      DIMENSION QROT(3,3)
      COMPLEX*16 RHO1(MAXK,NPMAX),RHO2(MAXK,NPMAX),
     *         DATA1(NPMAX),DATA2(NPMAX)
      COMPLEX*8 CORACCLIN(NPMAX),CROSSCORLIN(NPMAX)
      COMPLEX*8 CORACCLOG(NPMAX),CROSSCORLOG(NPMAX)
      character*3 label3,label
      character*6 label6
      character*2 label2

       complex*8 sqlmn
       common /AAA/ label3 
       common /BBB/ label6
       common /CCC/ label2
c   
       pi=4.0d0*datan(1.0d0)
       twopi=2.0d0*pi  
       imag=(0.0d0,1.0d0) 
c
         print *,' inserisci kmod min max'
         read(5,*) kmodmin,kmodmax
c        print *,' inserisci lmin and lmax'
c        read(5,*) lmin,lmax
c         kmodmin=2
         lmin=0
         lmax=2
c
c read q-vectors with same modulus (file produced with program makeq)
c
        do kmod= kmodmin, kmodmax
        write(label3,'(i3)') kmod
        call empty
        label=label3
        print *,' leggo dal file ','qvector.'//label3
        qmod=0.0d0

        open(unit=1,
     *    file='../../FQT/QVECTOR/qvector.'//label3,
     *             status='unknown')

        do i=1,maxk
         read(1,*,iostat=irr) kx(i),ky(i),kz(i)
         qmod=sqrt(real(kx(i)**2+ky(i)**2+kz(i)**2))+qmod
         if (irr.ne.0) goto 55
        end do
 55     nkv=i-1
        print *,' I am going to study # k-vectors=',nkv
        qmod=qmod/float(nkv)
        print *,' average modulus (units of twopi/box)',qmod
        close(1)
c
c             open files rho(l,m,n)
c
        do l=lmin,lmax
         do m=-l,l
          ifn=l*10+m
          write(label2,'(i2)') ifn
          call empty2
          open(unit=10+ifn,
     *   file='RHOTMP/ro.'//label2//'.k='//label,
c              form='unformatted',
     *          status='old')
         end do
        end do
c
c      open files sq-statico
c
          open(unit=90,
     *   file='sq0.all.dat'//label,
     *          status='unknown')
c
        do l=lmin,lmax
         do m=-l,l
          ifn=l*10+m
          print *,' generalized harmonic l,m = ',l,m
c
c   reads all data in a matrix  first l,m,n
c
           do inp=1,npmax
           read(10+ifn,*,iostat=irr) ti,nkv1,(rho1(k,inp),k=1,nkv)
           if (irr.ne.0) goto 555
           if (nkv1.ne.nkv) print *,' error in # q-vectors',inp
           time1(inp)=ti
c           print *,'time=',ti
           end do
 555       npt=inp-1
c           print *,' NPT, rho1=',npt,rho1(1,1)
           rewind(10+ifn)
c
c   now reads lp,mp
c
         do lp=lmin,lmax
          do mp=-lp,lp
           if (mp.ne.m) goto 1111
             ifn=lp*10+mp
c
c   reads all data in a matrix  now lp,mp
c
           do inp=1,npt
           read(10+ifn,*,iostat=irr) ti,nkv1,(rho2(k,inp),k=1,nkv)
           if (nkv1.ne.nkv) print *,' error in # q-vectors'
c           if (nfi.ne.time1(inp)) print *,' errore in time'
           if (irr.ne.0) print *,' some problems with file',ifn
           end do
c
           rewind(10+ifn)
c           print *,'rho2=',rho2(1,1)
c
c  the two matrices rho1(k,np) and rho2(k,np) contain all data
c  needed for the calculation of the autocorrelation functions
c 
c           call spacingtime(time1,npt,npc,ncic,base)
c
c
           print *,' lp,mp',lp,mp
         
          sqlmn=0.0d0
          nacc=0
          do k=1,nkv
          do nn=1,npt
            sqlmn=sqlmn+conjg(rho1(k,nn))*rho2(k,nn)
            nacc=nacc+1
c            print *,'sqklm ',k,sqlmn
          end do
         end do
c
c       at this point the correlation for l,m,lp,mp is ready
c
c        print *,'nacc,mols',nacc,mols
        adiv=float(nacc*mols)
c        print *,' adiv=',adiv
c
c        print *,'nacc=', nacc,' npt=',npt,' nkv=',nkv
             write(90,5555) l,m,lp,mp,
     *            kmod, real(sqlmn)/adiv, aimag(sqlmn)/adiv
5555             format(4i2,1x,i4,2e20.12)
c          print *,' SAVING',real(sqlmn)/adiv, aimag(sqlmn)/adiv
c
c       now.... next l,m,n,lp,mp,np
c
 1111     continue          
         end do
        end do
c
          end do
        end do

        end do
c
c  close files
c
        do l=lmin,lmax
         do m=-l,l
          ifn=l*10+m
          close(unit=10+ifn)
         end do
        end do
        close(90)
       stop
       end
c

c*************************************************************************
c
c   
c 
c*************************************************************************
c
       subroutine empty6
       character*1 a(6)
       common /BBB/ a
       do kk=1,6
         if (a(kk).eq.' ') a(kk)='0'
       end do
       return
       end
c
       subroutine empty2
       character*1 a(2)
       common /CCC/ a
       do kk=1,2
         if (a(kk).eq.' ') a(kk)='0'
       end do
       return
       end

       subroutine empty
       character*1 a(3)
       common /AAA/ a
       do kk=1,3
         if (a(kk).eq.' ') a(kk)='0'
       end do
       return
       end
