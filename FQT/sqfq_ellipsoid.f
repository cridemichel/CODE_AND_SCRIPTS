
c  Sq per  sistema  tridimensionale -- legge da NSCATT --- scrive su 
c                     non formattato
c
c     For cm 
c
c   NQM massimo vettore d'onda   -nqm < kx < nqm
c   NBIN numero di bin (circa nqm/qmesh)
c
c	Questa versione legge i parametri Qmin, Qmax, skip
c	dal file 'fqt.inp'
c
c         1/sqrt(n)* sum_i exp(iqr_i)
c
c      include 'mols'       

      parameter(maxconf=1000)
      PARAMETER(MOLS=30000)
      PARAMETER(NQM=100,NBIN=300,NPM=MOLS)
      PARAMETER(NROWMAX=1000,NTMAX=NROWMAX,NQMAX=300)
      dimension qex(npm,0:nqm),
     *     qey(npm,0:nqm),qez(npm,0:nqm)
      complex *16  rq1,rq2,imag,qex,qey,qez,cdumx,cdumy,cdumz,
     *     ex,ey,ez
      complex *8 DATARS1(NBIN),datars2(nbin)
      complex *8 rho1(nbin),rho2(nbin)
      character*3 label
      common /AAA/ label
      common /DDD/ label2
      integer calc
      REAL *8 sofq11(NQMAX),sofq12(NQMAX),sofq22(NQMAX)
      DIMENSION tconf(ntmax),nkv1(NQMAX)
      REAL*4 DATAR(MOLS,NTMAX,3)
      double precision anp,bx,by,bz,ener,time,press,temp
      double precision x(mols),y(mols),z(mols),vx,vy,vz
c    a,b,c sono i versori degli assi dell'ellissoide 	
      double precision uax(mols), uay(mols), uaz(mols)
      double precision ubx(mols), uby(mols), ubz(mols)
      double precision ucx(mols), ucy(mols), ucz(mols), FDF(mols)
      double precision pi,twopi,xxx,qa,qb,qc,qqx,qqy,qqz
      double precision a1,a2,b1,b2,c1,c2,dens1,dens2
      character*8 dummy   
      character*80 riga,dummil
      character*4  inizioriga 
      dimension itype(mols)
      character*2 nconfstring,label2,ncstring2
      character*3 ncstring
      
      character*132 nomefile,nomefile2,nomeqnc(maxconf)
      character*5 serie
      integer confskip

      ILONG=8
      IDOUBLE=8
C     parametri del modello
      a1 = 1.0d0
      b1 = 1.0d0
      c1 = 5.0d0
      a2 = 1.7d0
      b2 = 1.7d0
      c2 = 5.0d0
      pi=4.0d0*datan(1.0d0)
      dens1 = 58.0d0/(4*pi*a1*b1*c1/3)
      dens2 = 182.0d0/(4*pi*a2*b2*c2/3)

C ----   Legge i nomi di tutti i file di configurazione

      nf=0
      open(unit=71,file='filename.list',status='old')
      do ii=1,10000000
      read(71,'(a)',end=3742) nomefile
      nf=nf+1
      if (nf.gt.maxconf) pause 'Error in # of conf Qnc'
      nomeqnc(nf)=nomefile
      end do
 3742 continue
      close(71)

c 20   close(2)
c 5    read(71,'(a)',end=333) nomefile
    
      print *,'ci sono # files=',nf



      print *,' Insert modulus -> min & max & skip'
      open (unit=11, file='sqfq.inp', status='old')
      read(11,*) kmodmin,kmodmax,kskip
     

      print *,' PI=',pi
      twopi=2.0d0*pi
      imag=(0.0,1.0)


C     INIZIA Ciclo sulle Configurazioni

      nacc=0
      
      do kmod=kmodmin,kmodmax,kskip
      sofq11(kmod)=0.0d0
      sofq22(kmod)=0.0d0
      sofq12(kmod)=0.0d0
      end do

      do 100 kkkk=1,nf
            print *,nomeqnc(kkkk)

c     Legge coordinate

 
            open(1,file=nomeqnc(kkkk),status='old')
            do kkk=1,5000
               read(1,10) riga
 10            format(a80)
               inizioriga=riga
               if ('@@@ '.eq.inizioriga) then
                  read(1,*) dummy, nmols
                  read(1,*) dummy, nmola
            print *, nmols, nmola
               end if
               if ('nmax'.eq.inizioriga) then
                  read(1,*) dummy
c     print *, dummy
                  do i=1,nmols
                     read(1,*) x(i), y(i), z(i), 
     *       uax(i),uay(i),uaz(i),ubx(i),uby(i),ubz(i),
     *       ucx(i),ucy(i),ucz(i)
                  end do
                  do i=1,nmols
                     read(1,*)  aax, aay, aaz, ccx, ccy, ccz
                  end do
                  read(1,*) boxsize
                  goto 222
               end if
         
            end do

 222     continue
         bx=boxsize
         by=boxsize
         bz=boxsize
         
         boxhx=bx/2.0d0
         boxhy=by/2.0d0
         boxhz=bz/2.0d0

         print *, NMOLS, NMOLA, BX, BY, BZ

         np=nmols
         npt1=nmola
         npt2=nmols-nmola
c     Aggiorna Contatore sui file che legge
            nacc=nacc+1
            nqt=nqm
            q0=2.0*pi/bx
c            if (nacc.eq.1) then
c               print *,' q0=',q0,' nqt=',nqt
c            end if
            do i=1,np
               cdumx= -imag*q0*x(i)
               cdumy= -imag*q0*y(i)
               cdumz= -imag*q0*z(i)
               ex=cdexp(cdumx)
               ey=cdexp(cdumy)
               ez=cdexp(cdumz)
               do k=0,nqt
                  qex(i,k)=ex**k
                  qey(i,k)=ey**k
                  qez(i,k)=ez**k
               end do
            end do

             
C     Per ogni configurazione calcola le rho giuste
            do 1000 kmod=kmodmin,kmodmax,kskip
               write(label,'(i3)') kmod
               call empty
c     print *,' leggo dal file ','qvector.'//label       
               open(unit=33,
     *  file='/home/corezzi/QVECTOR/qvector.'//label,
     *              status='unknown')
c     ciclo che media sulle varie triplette di vettori d'onda
               do 444 inq=1,nbin 
                  read(33,*,iostat=irr) i,j,k
                  qmod=sqrt(real(i**2+j**2+k**2))+qmod
                  if (irr.ne.0) go to 55
c    calcola il fattore di forma e lo moltiplica 
                  qmodulo=sqrt(real(i**2+j**2+k**2))*q0
				qqx=real(i)*q0
				qqy=real(j)*q0
				qqz=real(k)*q0
				do l=1,np
				if (l.le.nmola) then
				 qa = a1*(qqx*uax(l)+qqy*uay(l)+qqz*uaz(l))
				 qb = b1*(qqx*ubx(l)+qqy*uby(l)+qqz*ubz(l))
				 qc = c1*(qqx*ucx(l)+qqy*ucy(l)+qqz*ucz(l))
				 xxx = sqrt(qa**2+qb**2+qc**2)  
	             FDF(l) = dsin(xxx)/xxx**3 - dcos(xxx)/xxx**2 
				 FDF(l) = FDF(l)*dens1*4*pi*a1*b1*c1
				else
				 qa = a2*(qqx*uax(l)+qqy*uay(l)+qqz*uaz(l))
				 qb = b2*(qqx*ubx(l)+qqy*uby(l)+qqz*ubz(l))
				 qc = c2*(qqx*ucx(l)+qqy*ucy(l)+qqz*ucz(l))
				 xxx = sqrt(qa**2+qb**2+qc**2)  
	             FDF(l) = dsin(xxx)/xxx**3 - dcos(xxx)/xxx**2  
				 FDF(l) = FDF(l)*dens2*4*pi*a2*b2*c2
				end if
				end do           
                  rq1=(0.0,0.0)
                  rq2=(0.0,0.0)
                  calc=5
                  if ((i.lt.0).and.(j.lt.0)) calc=1
                  if ((i.ge.0).and.(j.lt.0)) calc=2
                  if ((i.lt.0).and.(j.ge.0)) calc=3
                  if ((i.ge.0).and.(j.ge.0)) calc=4
               
			 goto (201,202,203,204) calc
               print *,' Some error',calc,i,j
               stop
 201           i=-i
               j=-j
               do l=1,np
c                  print *, itype(l)
                  if (l.le.nmola) then
                     rq1=rq1+FDF(l)*qez(l,k)/qex(l,i)/qey(l,j)
                  else
                     rq2=rq2+FDF(l)*qez(l,k)/qex(l,i)/qey(l,j)
                  end if   
               end do
               goto 300
 202           j=-j
               do l=1,np
                  if (l.le.nmola) then
                     rq1=rq1+FDF(l)*qex(l,i)/qey(l,j)*qez(l,k)
                  else
                     rq2=rq2+FDF(l)*qex(l,i)/qey(l,j)*qez(l,k) 
                  end if  
               end do
               goto 300
 203           i=-i
               do l=1,np
                  if (l.le.nmola) then
                     rq1=rq1+FDF(l)*qey(l,j)/qex(l,i)*qez(l,k)
                  else
                     rq2=rq2+FDF(l)*qey(l,j)/qex(l,i)*qez(l,k)
                  end if  
               end do
               goto 300
 204           continue
               do l=1,np
                  if (l.le.nmola) then
                     rq1=rq1+FDF(l)*qex(l,i)*qey(l,j)*qez(l,k)
                  else
                     rq2=rq2+FDF(l)*qex(l,i)*qey(l,j)*qez(l,k)
                  end if  
               end do
               goto 300
 300           continue
               anorm=sqrt(float(np))
               anorm1=sqrt(float(npt1))
               anorm2=sqrt(float(npt2))
			 rho1(inq)=rq1/anorm
               rho2(inq)=rq2/anorm
 444        end do
            
       close(33)
            
 55         continue
            nkv1(kmod)=inq-1
            print *,' I am going to study # k-vectors=',nkv1(kmod)
            qmod=qmod/float(nkv1(kmod))
            print *,' Average modulus (units of twopi/box)',qmod

            

         
         do k=1,nkv1(kmod)
               sofq11(kmod)=sofq11(kmod)+rho1(k)*conjg(rho1(k))
               sofq22(kmod)=sofq22(kmod)+rho2(k)*conjg(rho2(k))
               sofq12(kmod)=sofq12(kmod)+rho1(k)*conjg(rho2(k))
         end do
c     
c     save structure factor  - units of 2pi/box
c     


         

 1000 end do       
            
	    close(1)
c            call system('gzip '//filename)  
c            if (irr.ne.0) then
c               print *,' Tento di leggere senza successo '
c               stop
c            end if
c
c     ACCUMULA TEMPO E COORDINATE 
c
           
c            nt= nt+1
c            tconf(nt)=time
c            do kk=1,mols
c               datar(kk,nt,1)=x(kk)
c               datar(kk,nt,2)=y(kk)
c               datar(kk,nt,3)=z(kk)
c            end do
c            print *,nt
            
c            Continua ciclo sulle config.
 100     continue

         open(9,file='sqfq-partials.dat',status='unknown')
         do kmod=kmodmin,kmodmax,kskip
          
            dummy2=float(nacc)*nkv1(kmod)
            print *,' Normalization over ',dummy2
            write(9,*) kmod, real(sofq11(kmod)/dummy2),
     *       real(sofq22(kmod)/dummy2),real(sofq12(kmod)/dummy2),
     *       real((sofq11(kmod)+2.0d0*sofq12(kmod)+sofq22(kmod))/dummy2)
            
         end do
         close(9)

      stop
      end
      

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
      common /DDD/ a
      do kk=1,2
         if (a(kk).eq.' ') a(kk)='0'
      end do
      return
      end
      
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
c     
c
c
c     eliminate empty spaces in strings
c
      
      subroutine cleanstring(filei)
      character*132 filei,fileo
      character*1 recname(132)
      logical ex(132)
      equivalence (recname(1),fileo)
c     
      fileo=filei
c     
c     clean from empty spaces
c     
      do kk=1,132
         ex(kk)=.true.
         if (recname(kk).eq.' ') ex(kk)=.false.
      end do
c     
      ns=0
      do kk=1,132
         if (ex(kk)) then
            ns=ns+1
            recname(ns)=recname(kk)
         end if
      end do
c     
      do kk=ns+1,132
         recname(kk)=' '
      end do
      
c     
      filei=fileo
      return      
      end

