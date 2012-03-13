      SUBROUTINE qromo(func,a,b,ss,choose)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.d-4, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      print*,'a,b,fa,fb',a,b,func(a),func(b)
      print*,'s(j)',s(1),s(2),s(3),s(4),s(5),s(6)
      pause 'too many steps in qromo'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE midpnt(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs


      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.d-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
!	  if(i.ge.7)then
!	     yscal(i)=max(0.01,y(i))
!	  end if
12      continue

        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
	write(*,*)'yscal ',(yscal(ijk),ijk=1,nvar)
	write(*,*)'y     ',(y(ijk),ijk=1,nvar)
	write(*,*)'dydx  ',(dydx(ijk),ijk=1,nvar)
	write(*,*)h,eps
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
	write(*,*)hdid,hnext
	
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,
     *              derivs,ierror)
      INTEGER n,NMAX,ierror
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=511)
CU    USES derivs,rkck
      INTEGER i
      REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
	if(abs(yerr(i)/yscal(i)).ge.errmax)then
	   igreatest=i
	end if
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
!	write(*,*)' reducing step because of i= ',igreatest,
!     *    y(igreatest),yerr(igreatest),yscal(igreatest),
!     *    ' h=',h,' t=',x,' z=',y(3)
!        write(*,*)'f=',(y(ij),ij=1,n)
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
*        print*,'h,yscale/er',h,yscal(1),yscal(2),yerr(1),yerr(2)
        xnew=x+h
        if(xnew.eq.x)then
           ierror=1
           return
        end if
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=511)
CU    USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL*8 zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.d-6)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
     *'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.

      SUBROUTINE chebft(a,b,c,n,func)
      INTEGER n,NMAX
      REAL*8 a,b,c(n),func,PI
      EXTERNAL func
      PARAMETER (NMAX=50, PI=3.141592653589793d0)
      INTEGER j,k
      REAL*8 bma,bpa,fac,y,f(NMAX)
      DOUBLE PRECISION sum
      bma=0.5*(b-a)
      bpa=0.5*(b+a)
      do 11 k=1,n
        y=cos(PI*(k-0.5)/n)
        f(k)=func(y*bma+bpa)
11    continue
      fac=2./n
      do 13 j=1,n
        sum=0.d0
        do 12 k=1,n
          sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
12      continue
        c(j)=fac*sum
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL*8 chebev,a,b,x,c(m)
      INTEGER j
      REAL*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
*
**********************************************************
*
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=5,NUSE2=5)
CU    USES chebev
      REAL*8 xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.,1.,c1,NUSE1,xx)
      gam2=chebev(-1.,1.,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      INTEGER MAXIT
      REAL*8 ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-10,FPMIN=1.d-30,MAXIT=10000,XMIN=2.,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,
     *gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,
     *ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=1.d0/(b+d)
        c=b+1.d0/c
        del=c*d
        h=del*h
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessik; try asymptotic expansion'
1     continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
        ritemp=fact*ril+ripl
        fact=fact-xi
        ripl=fact*ritemp+ril
        ril=ritemp
12    continue
      f=ripl/ril
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
13      continue
        pause 'bessk series failed to converge'
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
14      continue
        pause 'bessik: failure to converge in cf2'
3       continue
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do 15 i=1,nl
        rktemp=(xmu+i)*xi2*rk1+rkmu
        rkmu=rk1
        rk1=rktemp
15    continue
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      RECURSIVE SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-2, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      write(*,*)'a,b,ss,func =',a,b,ss,func(a),func(b)
      pause 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     *REDMAX=1.d-5,REDMIN=.7,TINY=1.d-30,SCALMX=.1)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL*8 h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE rzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k
      REAL*8 b,b1,c,ddy,v,yy,d(NMAX,IMAX),fx(IMAX),x(IMAX)
      SAVE d,x
      x(iest)=xest
      if(iest.eq.1) then
        do 11 j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
11      continue
      else
        do 12 k=1,iest-1
          fx(k+1)=x(iest-k)/xest
12      continue
        do 14 j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.) then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
13        continue
          dy(j)=ddy
          yz(j)=yy
14      continue
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software -0)+.
