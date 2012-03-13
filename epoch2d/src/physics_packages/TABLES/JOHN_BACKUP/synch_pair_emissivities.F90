      function fsynchqmexact(eta,chi)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)

!       The synchrotron emissivity according to Erber ("exact" formula)

      real (kind=real8) , dimension(500) :: ylogtable,j1logtable,j2logtable,j3logtable
      real (kind=real8)  :: eta,chi,x,y,ylog,j1log,j2log,j3log,j1,j2,j3,m1,m2,m3
      real (kind=real8)  :: fsynchqmexact,delta
      integer :: npoints,n,nlow
      save iread,npoints,nlow,ylogtable,j1logtable,j2logtable,j3logtable
      logical iread /.true./

!      x is the photon energy as a fraction of the electron energy

      x=2.*chi/eta

      if(x.ge.1.)then
         fsynchqmexact=0.
         return
      end if

!      y is the quantum version of yclassical = \nu/\nu_{crit}

      y=x/3./eta/(1.-x)
      ylog=log10(y)

      if(ylog.gt.1..or.ylog.le.-5.)then
         fsynchqmexact=0.
         return
      end if

      m1=1.+1./(1.-x)**2
      m2=2./(1.-x)
      m3=x**2/(1.-x)**2
      
!      Read in the tabulated values of the j integrals

      if(iread)then
         open(unit=11,file='j1j2j3.table',status='unknown')
         read(11,*)npoints
        do n=1,npoints
            read(11,*)ylogtable(n),j1logtable(n),j2logtable(n),j3logtable(n)
         end do
         nlow=1
         close(11)
      end if
      iread=.false.
      
      call hunt(ylogtable,npoints,ylog,nlow)
      delta=ylogtable(nlow+1)-ylogtable(nlow)

      j1log=(j1logtable(nlow+1)*(ylog-ylogtable(nlow))+ &
            j1logtable(nlow)*(ylogtable(nlow+1)-ylog))/delta

      j2log=(j2logtable(nlow+1)*(ylog-ylogtable(nlow))+ &
            j2logtable(nlow)*(ylogtable(nlow+1)-ylog))/delta

      j3log=(j3logtable(nlow+1)*(ylog-ylogtable(nlow))+ &
            j3logtable(nlow)*(ylogtable(nlow+1)-ylog))/delta

      j1=10.**j1log
      j2=10.**j2log
      j3=10.**j3log

      fsynchqmexact=2.*sqrt(3.)/9./atan(1.)*chi**2/eta**4* &
                   (m1*j1+m2*j2+m3*j3)
!      write(*,*)' eta,chi,fsynchqmexact',eta,chi,fsynchqmexact
end function fsynchqmexact

     function gsokolov(eta)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)

      real (kind=real8) , dimension(500) :: etalogtable,glogtable
      real (kind=real8)  :: eta,etalog,gsoklog,gsokolov,delta
      integer :: npoints,n,nlow
      save iread,npoints,nlow,etalogtable,glogtable
      logical iread /.true./

      etalog=log10(eta)
      if(eta.le.1.e-5)then
         gsokolov=1.
         return
      else if(eta.ge.10.)then
         gsokolov=0.
         return
      end if

!      Read in the tabulated values of log(gsokolov)

      if(iread)then
         open(unit=11,file='gsokolov.table',status='unknown')
         read(11,*)npoints
        do n=1,npoints
            read(11,*)etalogtable(n),glogtable(n)
         end do
         nlow=1
         close(11)
      end if
      iread=.false.
      
      call hunt(etalogtable,npoints,etalog,nlow)
      delta=etalogtable(nlow+1)-etalogtable(nlow)

      gsoklog=(glogtable(nlow+1)*(etalog-etalogtable(nlow))+ &
            glogtable(nlow)*(etalogtable(nlow+1)-etalog))/delta

      gsokolov=10.**gsoklog
      return
end function gsokolov

     function hsokolov(eta)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)

      real (kind=real8) , dimension(500) :: etalogtable,hlogtable
      real (kind=real8)  :: eta,etalog,hsoklog,hsokolov,delta
      integer :: npoints,n,nlow
      save iread,npoints,nlow,etalogtable,hlogtable
      logical iread /.true./

      etalog=log10(eta)
      if(eta.le.1.e-5)then
         hsokolov=5.23599
         return
      else if(eta.ge.10.)then
         hsokolov=0.
         return
      end if

!      Read in the tabulated values of log(hsokolov)

      if(iread)then
         open(unit=11,file='hsokolov.table',status='unknown')
         read(11,*)npoints
        do n=1,npoints
            read(11,*)etalogtable(n),hlogtable(n)
         end do
         nlow=1
         close(11)
      end if
      iread=.false.
      
      call hunt(etalogtable,npoints,etalog,nlow)
      delta=etalogtable(nlow+1)-etalogtable(nlow)

      hsoklog=(hlogtable(nlow+1)*(etalog-etalogtable(nlow))+ &
            hlogtable(nlow)*(etalogtable(nlow+1)-etalog))/delta

      hsokolov=10.**hsoklog
      return
end function hsokolov

     function fairy(x)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)
     real (kind=real8)  :: fairy,x,xlog,zlow,zhigh,answer
      common/rfc/xlog
      external midpnt,s1func,s2func
      xlog=log10(x)
      if(x.le.1.)then
         zlow=0.
         zhigh=1.
         call qromo(s1func,zlow,zhigh,answer,midpnt)
         fairy=3./2.*x**2*answer
         return
      else
         zlow=0.
         zhigh=1.
         call qromo(s2func,zlow,zhigh,answer,midpnt)
         fairy=x*answer
         return
      end if
end function fairy

      function s1func(z)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)
      real (kind=real8)  :: s1func,z,xlog,y,xarg,xnu,rk,ri,rip,rkp

      common/rfc/xlog
      y=z**(-3./2.)
      xarg=10.**xlog*y
      xnu=5./3.

      if(xarg.lt.0.001)then
         rk=0.902745/2.*(xarg/2.)**(-5./3.)
      else if(xarg.le.10.)then
         call bessik(xarg,xnu,ri,rk,rip,rkp)
      else
         rk=sqrt(2.*atan(1.)/xarg)*exp(-xarg)
      end if
      s1func=z**(-5./2.)*rk
      return
end function s1func

      function s2func(z)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)
      real (kind=real8)  :: s2func,z,xlog,y,x,xarg,xnu,rk,ri,rip,rkp
      common/rfc/xlog
      y=-log(z)
      x=10.**xlog
      xarg=x+y
      xnu=5./3.

      if(xarg.lt.0.001)then
         rk=0.902745/2.*(xarg/2.)**(-5./3.)
      else if(xarg.le.10.)then
         call bessik(xarg,xnu,ri,rk,rip,rkp)
      else
         rk=sqrt(2.*atan(1.)/xarg)*exp(-xarg)
      end if
      s2func=rk/z
      return
end function s2func

        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function omegahat(eta)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)

!       The function defined by Erber for trident 
!       pair production (via a virtual photon)

      real (kind=real8) , dimension(1500) :: etalogtable,ologtable,tlogtable
      real (kind=real8)  :: eta,olog,delta
      real (kind=real8)  :: omegahat,xlogmin,xlogmax,etalog
      integer :: npoints,n,nlow
      save iread,npoints,nlow,etalogtable,ologtable
      logical iread /.true./

!      Read in the tabulated values of the log(omegahat)

      if(eta.le.1.e-5)then
         omegahat=0.
         return
      end if
      etalog=log10(eta)

      if(iread)then
         open(unit=11,file='pairprod.table',status='unknown')
         read(11,*)npoints,xlogmin,xlogmax
        do n=1,npoints
            read(11,*)etalogtable(n),ologtable(n),tlogtable(n)
         end do
         nlow=1
         close(11)
      end if
      iread=.false.
      
      call hunt(etalogtable,npoints,etalog,nlow)
      delta=etalogtable(nlow+1)-etalogtable(nlow)

      olog=(ologtable(nlow+1)*(etalog-etalogtable(nlow))+ &
            ologtable(nlow)*(etalogtable(nlow+1)-etalog))/delta

      omegahat=10.**olog

end function omegahat
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function tpair(chi)
      implicit none
      integer, parameter :: real8=selected_real_kind(15,40)

!       The function defined by Erber for  
!       pair production by a (real) photon

      real (kind=real8) , dimension(1500) :: chilogtable,ologtable,tlogtable
      real (kind=real8)  :: chi,tlog,delta
      real (kind=real8)  :: tpair,xlogmin,xlogmax,chilog
      integer :: npoints,n,nlow
      save iread,npoints,nlow,chilogtable,tlogtable
      logical iread /.true./

!      Read in the tabulated values of the log(omegahat)

      chilog=log10(chi)

      if(iread)then
         open(unit=11,file='pairprod.table',status='unknown')
         read(11,*)npoints,xlogmin,xlogmax
        do n=1,npoints
            read(11,*)chilogtable(n),ologtable(n),tlogtable(n)
         end do
         nlow=1
         close(11)
      end if
      iread=.false.
      
      call hunt(chilogtable,npoints,chilog,nlow)
      if(nlow.eq.0)then
        tpair=0.
      else
      delta=chilogtable(nlow+1)-chilogtable(nlow)

      tlog=(tlogtable(nlow+1)*(chilog-chilogtable(nlow))+ &
            tlogtable(nlow)*(chilogtable(nlow+1)-chilog))/delta

      tpair=10.**tlog
      end if
      return
end function tpair
         
