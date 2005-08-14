c
c   calcola sq(k,l,m,m,l,m,np) per un fissato modulo - scrive i risultati 
c                 in files  
c               NOTA l=lp,m=mp 

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     *  (MOLS=9261,
     *   LMAXX=1,KMAX=50,MAXK=KMAX,NPMAX=4000,MAXSLOG=60)
      DIMENSION TIME1(NPMAX)
c      DIMENSION SQLMN(0:LMAXX,0:LMAXX,0:LMAXX,1:KMAX)
      COMPLEX*16 RHO(0:LMAXX,-LMAXX:LMAXX,-LMAXX:LMAXX,1:KMAX)
      COMPLEX*16 IMAG,RCE,RCEK
      DIMENSION KX(MAXK),KY(MAXK),KZ(MAXK)  
      COMPLEX*16 R(0:LMAXX,-LMAXX:LMAXX,-LMAXX:LMAXX)    
      DIMENSION QROT(3,3)
      COMPLEX*16 RHO1(MAXK,NPMAX),RHO2(MAXK,NPMAX)
      COMPLEX*8  DATA1(NPMAX),DATA2(NPMAX)
      COMPLEX*8 CORACCLIN(NPMAX),CROSSCORLIN(NPMAX)
      COMPLEX*8 CORACCLOG(MAXSLOG,MAXSLOG),CROSSCORLOG(MAXSLOG,MAXSLOG)
      character*3 label3,label
      character*6 label6
      common /AAA/ label3
      common /BBB/ label6
      character*1 lab(-3:3)

       lab(0)='Z'
       lab(1)='U'
       lab(2)='D'
       lab(-1)='M'
       lab(-2)='N'
       lab(3)='T'
       lab(-3)='O'
c          
c
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi  
      imag=(0.0d0,1.0d0) 
      amh=1.008d0
      amo=16.00d0
      amwt=amo+2.0d0*amh
      lmin=0
      lmax=2
c
c read q-vectors with same modulus (file produced with program makeq)
c
        print *,' inserisci modulo '
        read(5,*) kmod
        write(label3,'(i3)') kmod
        call empty
        label=label3
        print *,' leggo dal file ','qvector.'//label3
        qmod=0.0d0
c
        open(unit=1,
     *    file='/usr1/users/fs/MD/QVECTOR/qvector.'//label3,
     *             status='old')
c

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
          do n=-l,l
          ifn=l*10+m
          write(label3,'(i3)') ifn
          call empty
          open(unit=1000+ifn,
     *   file='ro.'//label3//'.k='//label,form='unformatted',
     *          status='unknown')
          end do
         end do
        end do
c
        do l=lmin,lmax
         do m=-l,l
c          do n=-l,l
          ifn=l*100+m*10+n
c         print *,' generalized harmonic = ',ifn
c
c   reads all data in a matrix  first l,m,n
c
           do inp=1,npmax
           read(1000+ifn,iostat=irr) ti,nkv1,(rho1(k,inp),k=1,nkv)
           if (irr.ne.0) goto 555
           if (nkv1.ne.nkv) print *,' error in # q-vectors'
           time1(inp)=ti
           end do
 555       npt=inp-1
           rewind(1000+ifn)
c
c   now reads lp,mp,np
c
                     
         do lp=l,l
          do mp=m,m
c           do np=-lp,lp
             ifn=lp*10+mp
c            write(5,5555) ' doing',l,m,n,lp,mp,np
 5555         format(a6,6i4)
c
c   reads all data in a matrix  now lp,mp,np
c
           do inp=1,npt
           read(1000+ifn,iostat=irr) ti,nkv1,(rho2(k,inp),k=1,nkv)
           if (nkv1.ne.nkv) print *,' error in # q-vectors'
c           if (nfi.ne.time1(inp)) print *,' errore in time'
           if (irr.ne.0) print *,' some problems with file',ifn
           end do
           rewind(1000+ifn)
c
c  the two matrices rho1(k,np) and rho2(k,np) contain all data
c  needed for the calculation of the autocorrelation functions
c 
c           print *,' before time'
           call spacingtime(time1,npt,npc,ncic,base)
c            print *,' after time'
c
c  cleaning
c
            do j=1,ncic
              coracclin(j)=0.0d0
            end do
c     
            do j=1,npc
              do k=1,npc
              coracclog(j,k)=0.0d0
              end do
             end do
c
             nacc=0
c
          do k=1,nkv
c
c Note: croscorr subroutines do not perform any conmplex conjugation
c
c
          do nn=1,npt
            data1(nn)=conjg(rho1(k,nn))
            data2(nn)=rho2(k,nn)
          end do
c
c  calcola correlazione lineare e log -- not normalized
c
c       print *,' calling lin'
        CALL CROSSLIN (NPT,NCIC,NPC,DATA1,DATA2,CROSSCORLIN)
c       print *,' calling log'
        CALL CROSSLOG (NPT,NPC,NPC,DATA1,DATA2,CROSSCORLOG)
c       print *,' done'
c
c    media the autocorrelation function
c
         do j=1,NCIC
         CORACCLIN(J)=CORACCLIN(J)+CROSSCORLIN(J)
         end do
c


        do kk=1,NPC
          do ii=1,NPC-kk+1
          CORACCLOG(kk,ii)=CORACCLOG(kk,ii)+CROSSCORLOG(kk,ii)
          end do         
         end do


c
        NACC=NACC+1          
c       print *,' nacc=',nacc
c
        end do
c
c       at this point the correlation for l,m,n,lp,mp,np is ready
c
          open(unit=100,
     *   file='sqt.'//
     *   lab(l)//lab(m)//lab(n)//lab(lp)//lab(mp)//
     *   lab(np)//'.k='//label,
     *          status='unknown')
c
c     normalizza only on the # of molecules
c
        adiv=float(nacc*mols)
c       print *,' ADIV=',ADIV

c        do kk=1,NPC
         do kk=1,1 
        do i=2,NPC-kk+1
          CORACCLOG(kk,i)=CORACCLOG(kk,i)/adiv
          write(100,210)  (time1(i)-time1(1))*0.01,
     *        real(CORACCLOG(kk,i)),aimag(CORACCLOG(kk,i))
         end do
         end do
c
         do j=2,NCIC
          CORACCLIN(J)=CORACCLIN(J)/adiv
          write(100,210)    (time1(npc)-time1(1)+1)*0.01*(J-1), 
     *        real(CORACCLIN(J)),aimag(CORACCLIN(J))
         end do
        close(100)

210     format(f26.6,3x,f15.6,3x,f15.6)

c
c       now.... next l,m,n,lp,mp,np
c
c          end do
         end do
        end do
c
c          end do
         end do
        end do


c
c  close files
c
        do l=lmin,lmax
         do m=-l,l
          do n=-l,l
          ifn=l*100+m*10+n
          close(unit=1000+ifn)
          end do
         end do
        end do

       stop
       end
c
c**********************************************************************       
c
c calculate # points (npc*ncic),, # points per log cycle, # cycles
c
c**********************************************************************
c
        subroutine spacingtime(datac,ntime,npc,ncic,base)
        implicit double precision (a-h,o-z)
        dimension datac(ntime)
        logical first
        save first
        data first /.TRUE./
        save ncskip,isexp
        data ncskip,isexp /0,0/
c
c Questa diventa la base.
c
c
c
         if (first) then 
         print *,' PUNTI INVIATI ALLA SUB',NTIME
         do kkk=1,ntime
         print *,kkk,datac(kkk)
         end do
         open(unit=111,file='base-exp.dat',status='old')
         read(111,*) base,isexp
         close(111)
         print *,' LETTO BASE ED ISEXP',base,isexp
         if (ntime.lt.isexp) stop
         pause
         end if
c



66       npc=isexp+1
         NC=npc
         kk=npc
         if (first) print *,'Max Step=',datac(nc),' Conf #',kk
         NCIC=ntime/kk
         if (first) print *,' Numero cicli completi',ncic
         NTIME=KK*NCIC
         if (first) then
         print *,' Numero punti usati (multiplo di ciclo)',ntime
         print *,' Quanti cicli vuoi skippare per cattiva equil?'
         read(5,*) ncskip
         end if
         ncic=ncic-ncskip
         if (first) print *,' Uso allora solo ',ncic,' cicli'
         ntime=ntime-npc*ncskip
         if (first) print *,' Uso allora solo ',ntime,' punti '
c
c
         first=.false.
         return
         end
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
c     For complex DATA  - no normalization
c
      SUBROUTINE CROSSLIN(N,NC,NS,DATA1,DATA2,COR)
      parameter (NMAX=4000)
      real *4 tt, pcor(nmax)

      COMPLEX *8 DATA1(Nmax),COR(nmax),DATA2(nmax)
c                  
       do i=1,nc
         cor(i)=0.0
         pcor(i)=0.0
       end do
c
        do i=1,nc
         do j=1,n-(i-1)*ns,1
          cor(i)=cor(i)+data1(j)*data2(j+(i-1)*ns)
          pcor(i)=pcor(i)+1
        end do
       end do
c  Normalization  - - - sul numero dei punti
       do i=1,nc
         cor(i)=cor(i)/pcor(i)
       end do
      RETURN    
      END
c
c
c
      SUBROUTINE CROSSLOG (N,NC,NS,DATAX,DATAY,COR)
c
c
c
      parameter (NROW=4000, nmax=NROW, max=60)
      real *4 tt, pcor(max,max)
      COMPLEX *8 DATAX(Nmax),COR(max,max),DATAY(NMAX)
c
       do ik=1,nc
       do i=1,nc
         cor(ik,i)=0.0
         pcor(ik,i)=0.0
       end do
       end do
        do kk=1,nc
        do j=kk,n,ns
         do i= 1,nc-kk+1
         cor(kk,i)= cor(kk,i)+ (datax(j)*datay(j+i-1))
         pcor(kk,i)=pcor(kk,i)+1
        end do
       end do
       end do
c  Normalization       

      do kk=1,nc
       do i=1,nc-kk+1
         if (pcor(kk,i).eq.0) goto 6666
         cor(kk,i)=cor(kk,i)/pcor(kk,i)
6666     continue
       end do
       end do
      RETURN
      END
c

c*************************************************************************
c
       subroutine empty
       character*1 a(3)
       common /AAA/ a
       do kk=1,3
         if (a(kk).eq.' ') a(kk)='0'
       end do
       return
       end

