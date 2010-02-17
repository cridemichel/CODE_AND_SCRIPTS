       REAL*8 xmin, bret, TOL, LAM0, LAM1, LAM2
       PARAMETER(TOL=1.0D-15, LAM0=0.0, LAM1=0.5, LAM2=1.0)
       REAL*8 saA(3),saB(3),RA(3,3),RB(3,3)
       REAL*8 COMA(3),COMB(3)
       INTEGER ncall, i, j
!      NOTA sulle matrici di orientazione RA e RB:
!      RA(i,j) e RB(i,j) sono le matrici di orientazione dei 
!      dei due ellissoidi A e B ossia sono composte 
!      dai 3 vettori riga che rappresentano gli assi principali
!      dell'ellissoide. In particolare per passare dalle coordinate
!      del sistema di rif. del laboratorio (x) a quello del corpo rigido
!      (xp) si deve effettuare la trasformazione: xp = R*(x-xCM)
!      dove xCM è il centro di massa del corpo rigido.
       common / SIMPAR / saA,saB,RA,RB,COMA,COMB
!      USEDBRENT=.TRUE. to use brent with derivatives (5 times faster)
!      USEDIRINV=.TRUE. to use direct matrix inversion instead 
!      of LU decomposition (2 times faster!)
       common / LOGICPAR /  USEDIRINV, USEDBRENT
       LOGICAL USEDBRENT,USEDIRINV, READFROMFILE
!       PARAMETER(USEDBRENT=.TRUE.,USEDIRINV=.FALSE.)
       EXTERNAL Spw_wrap
       EXTERNAL SpwDer_wrap
!      set ellipsoids positions, orientations and parameters (i.e. semiaxes) 
       USEDBRENT=.TRUE.
       USEDIRINV=.TRUE.
       READFROMFILE=.TRUE.
       if (USEDBRENT) then
          print *, 'Using DBRENT'
       else
          print *, 'Using BRENT' 
       end if
       if (READFROMFILE) then
       OPEN(UNIT=11,FILE='ellips.pos', STATUS='OLD', iostat=irr) 
       if (irr.ne.0) then
       print *, 'BOH...'
       pause
       end if
       READ (11,*,IOSTAT=irr) COMA(1), COMA(2), COMA(3), 
     *  RA(1,1), RA(1,2), RA(1,3), RA(2,1), RA(2,2), RA(2,3),
     *  RA(3,1), RA(3,2), RA(3,3), saA(1), saA(2), saA(3)
       READ (11,*,IOSTAT=irr) COMB(1), COMB(2), COMB(3), 
     *  RB(1,1), RB(1,2), RB(1,3), RB(2,1), RB(2,2), RB(2,3),
     *  RB(3,1), RB(3,2), RB(3,3), saB(1), saB(2), saB(3)
       CLOSE(11)
!      it should be the 6th element read 
       print *, 'RA(1,...)=', RA(1,1), RA(1,2), RA(1,3)
       print *, 'RA(2,...)=', RA(2,1), RA(2,2), RA(2,3)
       print *, 'RA(3,...)=', RA(3,1), RA(3,2), RA(3,3)
       print *, 'semiaxes=',saA(1), saA(2), saA(3)
       else 
        saA(1)=2.
        saA(2)=2.
        saA(3)=1.
        saB(1)=2.
        saB(2)=2.
        saB(3)=1.
        COMA(1)=-1.99
        COMA(2)=0.
        COMA(3)=0.
        COMB(1)=1.99
        COMB(2)=0.
        COMB(3)=0.
        RA(1,1)=1.
        RA(1,2)=0.
        RA(1,3)=0.
        RA(2,1)=0.
        RA(2,2)=1.
        RA(2,3)=0.
        RA(3,1)=0.
        RA(3,2)=0.
        RA(3,3)=1.
        RB(1,1)=1.
        RB(1,2)=0.
        RB(1,3)=0.
        RB(2,1)=0.
        RB(2,2)=1.
        RB(2,3)=0.
        RB(3,1)=0.
        RB(3,2)=0.
        RB(3,3)=1.
       end if 
       do ncall=1,1000000
       if (USEDBRENT) then
        bret=dbrent(LAM0,LAM1,LAM2,SpwDer_wrap,TOL,xmin)
       else
        bret=brent(LAM0,LAM1,LAM2,Spw_wrap,TOL,xmin)
       end if
       end do
       PRINT *, 'F(A,B)=', -bret 
       end
       SUBROUTINE SpwDer_wrap(x,Sl,SlP,CALCSl,CALCSlP)
       EXTERNAL SpwDer
       REAL*8 x, Sl, SlP
       LOGICAL CALCSl, CALCSlP
       REAL*8 saA(3),saB(3),RA(3,3),RB(3,3)
       REAL*8 COMA(3),COMB(3)
       common / SIMPAR / saA,saB,RA,RB,COMA,COMB
       CALL SpwDer(x,Sl,SlP,CALCSL,CALCSlP) 
       RETURN
       END
       REAL*8 function Spw_wrap(x)
       REAL*8 saA(3),saB(3),RA(3,3),RB(3,3)
       REAL*8 COMA(3),COMB(3),boh, x
       common / SIMPAR / saA,saB,RA,RB,COMA,COMB
       Spw_wrap = Spw(x)
!      print *,"S1=", Spw_wrap
!       CALL SpwDer(x,Spw_wrap,boh,.TRUE.,.FALSE.) 
!       ,saA, COMA, RA, saB, COMB, RB) 
!      print *,'S2=', Spw_wrap
       RETURN
       end
       SUBROUTINE SpwDer(lambda,Sl,SlP,CALCSl,CALCSlP)
       EXTERNAL MATINV
       EXTERNAL MATINV33
       LOGICAL USEDBRENT,USEDIRINV
       common / LOGICPAR / USEDIRINV,USEDBRENT
       LOGICAL CALCSl,CALCSlP
       INTEGER a,i,j,n
       REAL*8 GinvR(3),saA(3), saB(3),B(3),H(3,3)
       REAL*8 Ginv(3,3),G(3,3),Ainv(3,3),Binv(3,3)
       REAL*8 R(3),RA(3,3),RB(3,3),COMA(3),COMB(3) 
       REAL*8 lambda,Sl,SlP
       common / SIMPAR / saA,saB,RA,RB,COMA,COMB
       do a=1,3
         R(a) = COMB(a)-COMA(a)
       end do  
       do i=1,3
        do j=1,3
         Ainv(i,j)=0.
         Binv(i,j)=0.
         do n=1,3
          Ainv(i,j)=Ainv(i,j)+saA(n)**2*RA(n,i)*RA(n,j) 
          Binv(i,j)=Binv(i,j)+saB(n)**2*RB(n,i)*RB(n,j)
         end do
         G(i,j)=(1.0D0 - lambda)*Ainv(i,j)+lambda*Binv(i,j) 
         if (CALCSlP) then
          H(i,j)=(1.0D0 - lambda)**2*Ainv(i,j)-lambda**2*Binv(i,j)
         end if
        end do
       end do
       if (USEDIRINV) then
         CALL MATINV33(G,Ginv)
       else
         CALL MATINV(G,Ginv)
       end if
       do i=1,3
        GinvR(i)=0.0d0
        do j=1,3 
          GinvR(i)=GinvR(i)+Ginv(i,j)*R(j)
        end do
       end do
!      se il valore di Slp passato è maggiore di 0 calcola la
!      derivata di S(lambda)
!      print *,'GinR=', GinvR(1), GinvR(2), GinvR(3)
       if (CALCSlP) then
        SlP = 0.0d0
        do i=1,3
         do j=1,3
          SlP = SlP + GinvR(i)*H(i,j)*GinvR(j) 
         end do
        end do
        SlP = -SlP
       end if
!      calculate also S(lambda) if Sl > 0.0
!      print *, 'R=', R(1), R(2), R(3)
       if (CALCSl) then
        Sl = 0.0d0
        do i=1,3
          Sl = Sl + R(i)*GinvR(i)  
        end do
        Sl = lambda*(1.0D0-lambda)*Sl
        Sl = -Sl
       end if
!       print *,'sl=', Sl, ' Slp=', SlP, ' l=',lambda 
       RETURN
       END
       REAL*8 function Spw(lambda)
!       ,saA,COMA,RA,saB,COMB,RB)
       LOGICAL USEDBRENT,USEDIRINV
       common / LOGICPAR / USEDIRINV, USEDBRENT
       INTEGER a,i,j,n
       EXTERNAL MATINV
       EXTERNAL MATINV33
       REAL*8 saA(3),saB(3),B(3)
       REAL*8 AB(3,3),ABinv(3,3),Ainv(3,3),Binv(3,3)
       REAL*8 R(3),RA(3,3),RB(3,3),COMA(3),COMB(3) 
       common / SIMPAR / saA,saB,RA,RB,COMA,COMB
       REAL*8 lambda
       do a=1,3
       R(a) = COMB(a)-COMA(a)
       end do  
       do i=1,3
       do j=1,3
       Ainv(i,j)=0.
       Binv(i,j)=0.
       do n=1,3
       Ainv(i,j)=Ainv(i,j)+saA(n)**2*RA(n,i)*RA(n,j) 
       Binv(i,j)=Binv(i,j)+saB(n)**2*RB(n,i)*RB(n,j)
       end do
       ABinv(i,j)=(1.0D0 - lambda)*Ainv(i,j)+lambda*Binv(i,j)
       end do
       end do
!      INVERT ABinv(3,3) matrix
       if (USEDIRINV) then
         CALL MATINV33(ABinv,AB)
       else
         CALL MATINV(ABinv,AB)
       end if
       Spw=0.0D0
       do i=1,3
       do j=1,3
       Spw=Spw+R(i)*AB(i,j)*R(j)
       end do
       end do
       Spw = Spw*lambda*(1.0D0-lambda)
       Spw = -Spw
       RETURN
       END 
       REAL*8 FUNCTION brent(ax,bx,cx,f,tol,xmin)
       INTEGER ITMAX
       REAL*8 ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
       EXTERNAL f
       PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
       INTEGER iter
       REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
       a=min(ax,cx) 
       b=max(ax,cx) 
       v=bx 
       w=v
       x=v
       e=0. 
       fx=f(x)
       fv=fx
       fw=fx
       do iter=1,ITMAX 
!      print *,'iter #', iter
       xm=0.5*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.*tol1
       if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3 
       if(abs(e).gt.tol1) then 
         r=(x-w)*(fx-fv)
         q=(x-v)*(fx-fw)
         p=(x-v)*q-(x-w)*r
         q=2.*(q-r)
       if(q.gt.0.) p=-p
         q=abs(q)
         etemp=e
         e=d
         if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.
     *    p.ge.q*(b-x)) goto 1
         d=p/q 
         u=x+d
         if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
         goto 2 
       endif
1      if(x.ge.xm) then
          e=a-x 
       else
         e=b-x
       endif
         d=CGOLD*e 
2      if(abs(d).ge.tol1) then 
         u=x+d 
       else
         u=x+sign(tol1,d)
       endif
       fu=f(u) 
       if(fu.le.fx) then 
         if(u.ge.x) then 
          a=x
         else
          b=x
         endif
       v=w
       fv=fw
       w=x
       fw=fx
       x=u
       fx=fu
       else
        if(u.lt.x) then
          a=u
        else
          b=u
        endif
       if(fu.le.fw .or. w.eq.x) then
         v=w
         fv=fw
         w=u
         fw=fu
       else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
         v=u
         fv=fu
       endif
       endif 
       end do 
       pause 'brent exceed maximum iterations'
3      xmin=x 
       brent=fx
       RETURN
       END 
       SUBROUTINE ludcmp(a,n,np,indx,d)
       INTEGER n,np,indx(n),NMAX
       REAL*8 d,a(np,np),TINY
       PARAMETER (NMAX=500,TINY=1.0e-20)       
       INTEGER i,imax,j,k
       REAL*8 aamax,dum,sum,vv(NMAX) 
       d=1.0D0 
       do i=1,n 
       aamax=0.0d0
       do j=1,n
       if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       end do 
       if (aamax.eq.0.) pause 'singular matrix in ludcmp' 
       vv(i)=1./aamax 
       end do
       do j=1,n 
       do i=1,j-1 
       sum=a(i,j)
       do k=1,i-1
       sum=sum-a(i,k)*a(k,j)
       end do 
       a(i,j)=sum
       end do 
       aamax=0.0 
       do i=j,n 
       sum=a(i,j) 
       do k=1,j-1
       sum=sum-a(i,k)*a(k,j)
       end do 
       a(i,j)=sum
       dum=vv(i)*abs(sum) 
       if (dum.ge.aamax) then 
        imax=i
        aamax=dum
       endif
       end do 
       if (j.ne.imax)then 
        do k=1,n 
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
       end do 
       d=-d 
       vv(imax)=vv(j) 
       endif
       indx(j)=imax
       if(a(j,j).eq.0.)a(j,j)=TINY
       if(j.ne.n)then 
         dum=1./a(j,j)
       do  i=j+1,n
       a(i,j)=a(i,j)*dum
       end do 
       endif
       end do  
       return
       END
       SUBROUTINE lubksb(a,n,np,indx,b)
       INTEGER n,np,indx(n)
       REAL*8 a(np,np),b(n) 
       INTEGER i,ii,j,ll
       REAL*8 sum
       ii=0 
       do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
         do j=ii,i-1
         sum=sum-a(i,j)*b(j)
         end do 
        else if (sum.ne.0.) then
         ii=i 
        endif 
        b(i)=sum
        end do 
        do i=n,1,-1 
         sum=b(i)
         do j=i+1,n
          sum=sum-a(i,j)*b(j)
         end do
         b(i)=sum/a(i,i) 
         end do
        return 
        END
       SUBROUTINE MATINV(A,Ainv)
       REAL*8 A(3,3),D,Ainv(3,3)
       INTEGER i,j,INDX(3)
       do i=1,3 
       do j=1,3
         Ainv(i,j)=0.0d0
       end do 
       Ainv(i,i)=1.0d0
       end do 
       CALL ludcmp(A,3,3,INDX,D) 
       do j=1,3
        CALL lubksb(A,3,3,INDX,Ainv(1,j))
!      Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
!      address of the jth column of y.
       end do
       RETURN
       END
       SUBROUTINE MATINV33(A,Ainv)
       REAL*8 A(3,3),Ainv(3,3)
       REAL*8 DET
       DET = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)
     * *A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)
     * *A(3,1)
       Ainv(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
       Ainv(1,2)=A(1,3)*A(3,2)-A(1,2)*A(3,3)
       Ainv(1,3)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
       Ainv(2,1)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
       Ainv(2,2)=A(1,1)*A(3,3)-A(3,1)*A(1,3)
       Ainv(2,3)=A(1,3)*A(2,1)-A(1,1)*A(2,3)
       Ainv(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
       Ainv(3,2)=A(1,2)*A(3,1)-A(1,1)*A(3,2)
       Ainv(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
       do i=1, 3
       do j=1, 3
          Ainv(i,j)=Ainv(i,j)/DET
       end do
       end do
       RETURN
       END
       FUNCTION dbrent(ax,bx,cx,f,tol,xmin)
!      Given a function f and its derivative function df, and given a bracketing triplet of abscissas
!      ax, bx, cx [such that bx is between ax and cx, and f(bx) is less than both f(ax) and
!      f(cx)], this routine isolates the minimum to a fractional precision of about tol using
!      a modification of Brent’s method that uses derivatives. The abscissa of the minimum is
!      returned as xmin, and the minimum function value is returned as dbrent, the returned
!      function value.
       INTEGER ITMAX
       REAL*8 dbrent,ax,bx,cx,tol,xmin,df,ZEPS
       EXTERNAL f
       PARAMETER (ITMAX=100,ZEPS=1.0e-10)
       INTEGER iter
       REAL*8 a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,
     * u,u1,u2,v,w,x,xm
!      Comments following will point out onlydi fferences from the routine brent. Read that
!      routine first.
       LOGICAL ok1,ok2 
       a=min(ax,cx) 
       b=max(ax,cx)
       v=bx
       w=v
       x=v
       e=0.
!      fx=f(x)
!      f now calculates both f and its derivative
       CALL f(x,fx,dx,.TRUE.,.TRUE.)
       fv=fx
       fw=fx
!      dx=df(x) 
       dv=dx
       dw=dx
       do iter=1,ITMAX
!      print *, 'iter=', iter
       xm=0.5*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.*tol1
       if(abs(x-xm).le.(tol2-.5*(b-a))) goto 13
       if(abs(e).gt.tol1) then
       d1=2.*(b-a) 
       d2=d1
       if(dw.ne.dx) d1=(w-x)*dx/(dx-dw) 
       if(dv.ne.dx) d2=(v-x)*dx/(dx-dv) 
!Which of these two estimates of d shall we take? We will insist that theyb e within
!the bracket, and on the side pointed to bythe derivative at x:
       u1=x+d1
       u2=x+d2
       ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
       ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
       olde=e 
!Movement on the step before last.
       e=d
       if(.not.(ok1.or.ok2))then
! Take onlya n acceptable d, and if both
!are acceptable, then take the smallest
!one.
       goto 11
       else if (ok1.and.ok2)then
       if(abs(d1).lt.abs(d2))then
        d=d1
       else
        d=d2
       endif
       else if (ok1)then
        d=d1
       else
        d=d2
       end if
       if(abs(d).gt.abs(0.5*olde))goto 11
        u=x+d
       if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
        goto 12
       endif
11      if(dx.ge.0.) then 
! Decide which segment bythe sign of the derivative.
        e=a-x
       else
        e=b-x
       end if 
       d=0.5*e 
!Bisect, not golden section.
12     if(abs(d).ge.tol1) then
        u=x+d
!       fu=f(u)
        CALL f(u,fu,du,.TRUE.,.TRUE.)
       else
        u=x+sign(tol1,d)
!       fu=f(u)
        CALL f(u,fu,du,.TRUE.,.TRUE.)
        if(fu.gt.fx) goto 13 
! If the minimum step in the downhill direction takes us uphill,
       endif 
!then we are done.
!       du=df(u) ! now we already have derivative calculated by f (see above) 
!Now all the housekeeping, sigh.
       if(fu.le.fx) then
        if(u.ge.x) then
         a=x
        else
         b=x
       end if
        v=w
        fv=fw
        dv=dw
        w=x
        fw=fx
        dw=dx
        x=u
        fx=fu
        dx=du
       else
        if(u.lt.x) then
         a=u
        else
         b=u
       end if
       if(fu.le.fw .or. w.eq.x) then
        v=w
        fv=fw
        dv=dw
        w=u
        fw=fu
        dw=du
       else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
        v=u
        fv=fu
        dv=du
       endif
       endif
       end do 
       pause 'dbrent exceeded maximum iterations'
13     xmin=x
       dbrent=fx
       return
       END
