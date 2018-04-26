c
c
      subroutine dskin(initial,energy,pitch,rho,fdist)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      data initialized/0/

c..................................................................
c
c     This f77 subroutine reads an ASCII file "diskf" produced by 
c     and described in the subroutine dskout in the 
c     CQL3D Fokker-Planck code.
c
c     Interpolate f for value at energy,pitch,rho ==> fdist
c
c     initial=1  ==> read file diskf,  else already open.
c
c     The parameters (ending in "a" below, and dimensions of the 
c     variables to be read into must be set in the file dskin.h,
c     in accord with the dimensions (iy,jx,lrz,lrzmax) in "diskf":
c
c       dskin.h:
c       iya=iy(1), jxa=jx, lrza=lrz, ngena=ngen
c
c     Pointered storage could be used instead, if this procedure
c       becomes too tedious.
c
c     Bob Harvey, June, 2000
c
c..................................................................

      include 'dskin.h'

      character*8  idskf
      character*80 line

      IF (initial .eq. 1) THEN

         initialized=1

c      write(*,*) 'dskin: initial = ',initial

      idskf="diskf"
      open(unit=4,file=idskf,status='old')

c..................................................................
c     The input namelist file has been  transcribed onto the beginning
c     of file idskf and we space over it to the FP data.
c..................................................................

 1    read(4,1003) line
      if (line(1:27).eq."***Begin computed output***") go to 2
      go to 1
 2    continue 

c..................................................................
c
c    CQL3D has a facility to computer FP solutions on a subset of
c      the full radial mesh on which plasma parameters are specified.
c      
c    In the following, we are only interested in the given data
c    which is at flux surfaces numbered ll=1:lrz
c    (thus we can consider lrindx(ll)=ll).
c
c
c
c     FROM CQL3D:dskout.f
c     In the following disk write to file named idskf:
c     This subroutine is called to write data for each FP'd
c          flux surface.
c     ll=  FP flux surface number 
c          (ll=1:lrz, lrz.le.lrzmax, see cqlinput_help))
c          (lrindx(ll) gives flux surface number on the full CQL3D
c                     radial mesh. (lrindx(ll)=ll if lrz=lrzmax,
c                     and using cql3d mode(cqlpmod="disabled"))).
c     iy(ll),jx= dimensions in theta and u(momentum/mass)
c           (In the case where iy varies with ll, iy(1) will be greatest.)
c     lrz= number of flux surfaces FP'd.
c     lrzmax= number of flux surfaces, including any not FP'd.
c     x = momentum-per-mass(nomalized to maximum 1.)
c         (same for each flux surface lrindx(ll))
c     y = theta(radians) mesh at  each flux surface lrindx(ll)
c     rovera= normalized radius (ll)
c             (rho, see Hinton and Haseltine for non-circ).
c             (generally ~sqrt(tor. flux), but other coords are
c              available for setup in the CQL3D input file.)
c     elecfld = toroidal electric field (volts/cm)
c     bthr - the poloidal magnetic field at theta-poloidal = pi/2.
c     btoru - the toroidal magnetic field at the same position. 
c     bthr0, btor0, are the poloidal and  toroidal
c         magnetic fields at the outer midplane of the flux surface. 
c     reden= electron density at minimum B point on flux surface.
c     temp= initial electron temperature (keV)
c     radmin= plasma minor radius (cms).
c     vnorm= normalization momentum-per-mass (maximum on grid) (cm/sec)
c     vmaxdvt= vnorm/(temp/mass)**0.5
c     eovedd= electric  field, normalized to Driecer field 
c             (calc'd in sub restvty).
c
c     distribution function normalised so that
c         integral( (dx)**3 f) = density at minumum B point.
c
c..................................................................
c
      l=1
c     Loop back to 10 reading data for each flux surfaces

 10   read(4,1004)  ll, iy(l),jx,lrz,lrzmax,ngen
      l=ll+1

c      write(*,*) ll, iy(ll),jx,lrz,lrzmax,ngen

c..................................................................
c       Check dimensions from dskin.h
c..................................................................

      if (iy(1).ne.iya .or. jx.ne.jxa .or. lrz.ne.lrza
     +     .or. ngen.ne.ngena) stop 'Check dskin.h parameters'
      
      
      read(4,1004)  itl(ll),itu(ll)
      read(4,1005)  (x(j),j=1,jx)
      read(4,1005)  (y(i,ll),i=1,iy(ll))
      do 20 k=1,ngen
         read(4,1005)  bnumb(k),fmass(k)
         read(4,1005)  rovera(ll),elecfld(ll),
     +        bthr(ll),btoru(ll)
         read(4,1005)  bthr0(ll),btor0(ll),
     +        reden(ll,k),temp(ll,k)
         read(4,1005)  radmin,vnorm,vmaxdvt,eovedd
         read(4,1005)  ((f(i,j,ll,k),i=1,iy(ll)),j=1,jx)
 20   continue
      if (ll.lt.lrz) go to 10
      
      
 1003 format(a80)
 1004 format(16i5)
 1005 format(5e16.8)
 1006 format(a27)
      close(unit=4)

      ELSE     !  END OF INITIALIZATION

c..................................................................
c     BEGIN: Interpolate f for value at energy,pitch,rho
c     Single species, although easy to generalize
c..................................................................


c      write(*,*) 'dskin: initial = ',initial
      if (initialized.ne.1) stop 'Stop in dskin: Data not initialized.'
         
      kk=1
      clight=2.99792458d10
      vel=sqrt(2.*energy*1.d3*1.6e-12/fmass(1))
      xvel=vel/vnorm
c     First, determine if the particle at (energy,pitch,radius)
c       is trapped.
c     Use lookup to find weights for distribution at neighboring
c       tabulated points.  Use closest pitch angle mesh to
c       determine whether trapped.

      if (rho .le. rovera(1)) then
         ll=1
      elseif (rho .ge. rovera(lrz)) then
         ll=lrz
      else
         call lookup(rho,rovera,lrzmax,weightu,weightl,ll)
         if (weightl.ge.weightu) ll=ll-1
      endif
               
      lpitch=luf(pitch,y(1,ll),iy(ll))
      if (lpitch.gt.(iy(ll)/2)) lpitch=lpitch-1 ! Gives symmetry.
      itrap=0                                   ! Trapping flag.
      if (lpitch.gt.itl(ll) .and. lpitch.lt.itu(ll)) itrap=1

c     Shift radial location of trapped particles by a half banana width
c     This will be an outward/inward shift for pos/neg bthr*bnumb
c     (bnumb is charge number, including sign).
c     Transiting particles are unshifted.

      rr=rho
      charge=4.8032d-10
      if (itrap.eq.1) then
         qb_mc=bnumb(kk)*charge*bthr(ll)/(fmass(kk)*clight)
         rban=xvel*cos(pitch)*vnorm/qb_mc
c     Shift radial location of trapped particles by a half banana width
c     This will be an outward/inward shift for pos/neg bthr*bnumb
c     (bnumb is charge number, including sign).
         rrr=rho+rban/radmin
         if (rr.lt. 0.d0) rr=0.d0
         if (rr.gt. 1.d0)  rr=1.d0
      endif

c..................................................................
c     Use lookup to find weights for distribution at neighboring
c     tabulated points.  Interpolate to velocity and pitch angle
c     to the specified (xvel,pitch) values.
c..................................................................

c     Stay on the grid:
      call lookup(rr,rovera,lrz,rweightu,rweightl,lrr)
      if (lrr.gt.lrz) then
         lrr=lrz
         rweightl=0.d0
         rweightu=1.d0
      elseif (lrr.le.1) then
         llr=2
          rweightl=1.d0
          rweightu=0.d0
      endif
      call lookup(pitch,y(1,ll),iy(ll),pweightu,pweightl,iyy)
      if (iyy.ge.iy(ll)) then
         iyy=iy(ll)
         pweightu=1.d0
         pweightl=0.d0
      elseif (iyy.le.1) then
         iyy=2
         pweightu=0.d0
         pweightl=1.d0
      endif
      call lookup(xvel,x,jx,xweightu,xweightl,jxx)
      if (jxx.ge.jx) then
         jxx=jx
         xweightu=1.d0
         xweightl=0.d0
      elseif (jxx.le.1) then
         jxx=2
         xweightu=0.d0
         xweightl=1.d0
      endif

c     Tri-linear interpolate (Eliminate f<0 values).

      fdist=rweightl*(pweightl*xweightl*max(0d0,f(iyy-1,jxx-1,lrr-1,1))
     +               +pweightu*xweightl*max(0d0,f(iyy,jxx-1,lrr-1,1))
     +               +pweightl*xweightu*max(0d0,f(iyy-1,jxx,lrr-1,1))
     +               +pweightu*xweightu*max(0d0,f(iyy,jxx,lrr-1,1)))
     +     +rweightu*(pweightl*xweightl*max(0d0,f(iyy-1,jxx-1,lrr,1))
     +               +pweightu*xweightl*max(0d0,f(iyy,jxx-1,lrr,1))
     +               +pweightl*xweightu*max(0d0,f(iyy-1,jxx,lrr,1))
     +               +pweightu*xweightu*max(0d0,f(iyy,jxx,lrr,1)))

c      Following write statments for code checkout.
c      write(*,*) 'dskin: fdist, f.....................:', fdist
c      write(*,*) 'dskin: rweightl,pweightl,xweightl'
c      write(*,*) '           ',rweightl,pweightl,xweightl
c      write(*,*) f(iyy-1,jxx-1,lrr-1,1),f(iyy,jxx-1,lrr-1,1)
c      write(*,*) f(iyy-1,jxx,lrr-1,1),f(iyy,jxx,lrr-1,1)
c      write(*,*) 'dskin: rweightu', rweightu
c      write(*,*) f(iyy-1,jxx-1,lrr,1),f(iyy,jxx-1,lrr,1)
c      write(*,*) f(iyy-1,jxx,lrr,1),f(iyy,jxx,lrr,1)
c      write(*,*) 'dskin..................................'


c     END of evaluation of f(energy,pitch,rho)

      ENDIF          
      
      return
      end

c
c
      subroutine lookup(x,xarray,length,weightu,weightl,lement)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine uses luf to do a table look up. Then it interpolates
c     to gain a bit of accuracy.
c     x is the argument; xarray is the monotonic array; length
c     is the length of the array. lement is the first index such
c     that xarray(lement).gt.x.
c..................................................................

      save
      dimension xarray(*)
      lement=luf(x,xarray,length)
      weightl=(xarray(lement)-x)/(xarray(lement)-xarray(lement-1))
      weightu=1.-weightl
      return
      end



      integer function luf(x,table,n)
      implicit integer (i-n), real*8 (a-h,o-z)
c
c     THIS ROUTINE SHOULD BE A BINARY SEARCH.  IT NEEDS
C        WORK!
c     luf(x,table,n) (MATHLIB) which is a function returning the index
c        of the first element in the table that is greater than x.
c        Elements must be strictly increasing. x.gt.table(n)==>n+1.
c
c     find first index such that table(luf).gt.x
c
c
      dimension table(n)
c
      do i=1,n
        if (table(i) .gt. x) go to 10
      end do
 10   continue
c     luf = 1 if x.lt.table(1) and luf=n+1 if x>ge.table(n)
      luf = i
c
      return
      end

