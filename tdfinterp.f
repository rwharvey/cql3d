c
c
      subroutine tdfinterp(u,pitch,rho,polang,fdist)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c
c     Interpolate f for value at u,pitch,rho,polang ==> fdist
c          u=momentum-per-rest-mass (cm/sec)
c          pitch=pitch angle at polang [0,pi]
c          rho=normalized radial coord
c          polang=poloidal angle, in [0,2*pi) 
c                                    [or [0,pi] if eqsym.ne.'none']
c
c     f will be in code units (i.e., f_cgs*vnorm**3). where
c     f_cgs is in cgs units 1/(volume*velocity**3).
c     Normalization is thus such that integral{f*d**3x}=ne(#/cm**3).
c
c     Zero-orbit-width (readily generalized).
c
c     Bob Harvey, June, 2000/April, 2010
c
c..................................................................

      include 'param.h'
      include 'comm.h'
      real*8, dimension(:), allocatable ::  rovera_ 
      data ifirst/1/, icount/0/


c.......................................................................
c
c    CQL3D has a facility to computer FP solutions on a subset of
c      the full radial mesh on which plasma parameters are specified.
c      
c    In the following, we are only interested in the given data
c    which is at flux surfaces numbered ll=1:lrz
c    (thus we can consider lrindx(ll)=ll).
c
c
c     distribution function normalized so that
c         integral( (dx)**3 f) = density at minumum B point.
c
c.......................................................................
c
c.......................................................................
c     Fill rovera_, containing FP'd flux surfaces (lrz.le.lrzmax).
c     f is given on this restricted set of surfaces f(,,,1:lrz).
c.......................................................................

      if (ifirst.eq.1 ) then
         allocate (rovera_(lrz), STAT=istat)
         do ll=1,lrz
            rovera_(ll)=rovera(lrindx(ll))
         enddo
         ifirst=0
      endif
c
c.......................................................................
c     For given poloidal angle and pitch, determine equatorial 
c     plane pitch angle.
c     ll= radial index in rovera_(1:lrz) grid
c     lr= radial index in rovera(1:lrzmax) grid
c     lookup() adjusts weights for endpoints.
c.......................................................................

      call lookup(rho,rovera_,lrz,rweightu_,rweightl_,ll)
      call lookup(rho,rovera,lrzmax,rweightu,rweightl,lr)

c     Lookup polang on surrounding pol[1:lz,1:lrzmax] meshes, 
c     interpolate for B/B0, and obtain equatorial plane pitch angle.

      if (eqsym.ne."none" .and. polang.gt.pi) 
     +     polang=2*pi-polang   !assumes up-down symm

      call lookup(polang,pol(1,lr-1),lz,pweightu1,pweightl1,lt1)
      call lookup(polang,pol(1,lr),lz,pweightu2,pweightl2,lt2)
      bdb01=pweightl1*bbpsi(lt1,lr-1)+pweightu1*bbpsi(lt1,lr-1)
      bdb02=pweightl2*bbpsi(lt2,lr)+pweightu2*bbpsi(lt2,lr)
      bdb0=rweightl*bdb01+rweightu*bdb02

cBH110817      sin_pitch02=sin(pitch)**2/bdb0
cBH110817      pitch0=asin(sqrt(sin_pitch02))  ! Corresponding midplane pitch
cBH110817: Fix for pitch.gt.pio2.  Since tdfinterp only used so far
cBH110817: for perp viewing NPA in C-Mod, this hasn't caused a problem.
      sin_pitch0=sin(pitch)/sqrt(bdb0)
      pitch0=asin(sin_pitch0)  ! Corresponding midplane pitch
      if (pitch.gt.pio2) pitch0=pi-pitch0


c.......................................................................
c     BEGIN: Interpolate f for value at u,pitch0,rho
c     Single species, although easy to generalize
c.......................................................................

      k=1
      xvel=u/vnorm

c     First, determine if the particle at (u,pitch0,radius)
c       is trapped.
c     Use lookup to find weights for distribution at neighboring
c       tabulated points.  Use closest pitch angle mesh to
c       determine whether trapped.
               
c$$$      if (rweightl.ge.rweightu) llp=lr-1
c$$$      lpitch=luf(pitch,y(1,llp),iy_(llp))
c$$$      if (lpitch.gt.(iy_(llp)/2)) lpitch=lpitch-1 ! Gives symmetry.
c$$$      itrap=0                                   ! Trapping flag.
c$$$      if (lpitch.gt.itl_(llp) .and. lpitch.lt.itu_(llp)) itrap=1
c$$$
c$$$c     Shift radial location of trapped particles by a half banana width
c$$$c     This will be an outward/inward shift for pos/neg bthr*bnumb
c$$$c     (bnumb is charge number, including sign).
c$$$c     Transiting particles are unshifted.
c$$$
c$$$      rr=rho
c$$$      charge=4.8032d-10
c$$$      if (itrap.eq.1) then
c$$$         qb_mc=bnumb(k)*charge*bthr(lr)/(fmass(k)*clight)
c$$$         rban=xvel*cos(pitch)*vnorm/qb_mc
c$$$c     Shift radial location of trapped particles by a half banana width
c$$$c     This will be an outward/inward shift for pos/neg bthr*bnumb
c$$$c     (bnumb is charge number, including sign).
c$$$         rrr=rho+rban/radmin
c$$$         if (rr.lt. 0.d0) rr=0.d0
c$$$         if (rr.gt. 1.d0)  rr=1.d0
c$$$      endif

c..................................................................
c     Use lookup to find weights for distribution at neighboring
c     tabulated points.  Interpolate to velocity and pitch angle
c     to the specified (xvel,pitch) values.
c..................................................................

c     lookup pitch0 on y-grid corresponding to neighboring f()
      lrm=lrindx(ll-1)
      call lookup(pitch0,y(1,lr),iy_(lr),pweightu,pweightl,iyy)
      call lookup(pitch0,y(1,lrm),iy_(lrm),pweightum,pweightlm,iyym)
      call lookup(xvel,x,jx,xweightu,xweightl,jxx)

c     Tri-linear interpolate (Eliminate f<0 values).

c$$$      fdist=
c$$$     +     rweightl_*(pweightlm*xweightl*max(0d0,f(iyym-1,jxx-1,k,ll-1))
c$$$     +               +pweightum*xweightl*max(0d0,f(iyym,jxx-1,k,ll-1))
c$$$     +               +pweightlm*xweightu*max(0d0,f(iyym-1,jxx,k,ll-1))
c$$$     +               +pweightum*xweightu*max(0d0,f(iyym,jxx,k,ll-1)))
c$$$     +     +rweightu_*(pweightl*xweightl*max(0d0,f(iyy-1,jxx-1,k,ll))
c$$$     +               +pweightu*xweightl*max(0d0,f(iyy,jxx-1,k,ll))
c$$$     +               +pweightl*xweightu*max(0d0,f(iyy-1,jxx,k,ll))
c$$$     +               +pweightu*xweightu*max(0d0,f(iyy,jxx,k,ll)))

c     Tri-linear interpolate log(f) (and eliminate log(f)<-200 values).
c     (The log of constant temperature distributions are linear vs vel,
c      making linear interpolation appropriate.)
      fmin= 1.383896526737250d-87  ! log(1.383896526737250d-87)= -200.
c      fmin= 3.720075976020916d-44  ! log(3.720075976020916d-44)= -100.


c     Adapt to calc of fdist from either favg or f:
      if (f4d_out.eq."tavg" .and. n.eq.nstop) then
         fmmm=favg(iyym-1,jxx-1,k,ll-1)
         f0mm=favg(iyym,  jxx-1,k,ll-1)
         fm0m=favg(iyym-1,jxx,  k,ll-1)
         f00m=favg(iyym,  jxx,  k,ll-1)
         fmm0=favg(iyym-1,jxx-1,k,ll)
         f0m0=favg(iyym,  jxx-1,k,ll)
         fm00=favg(iyym-1,jxx,  k,ll)
         f000=favg(iyym,  jxx,  k,ll)
      else
         fmmm=f(iyym-1,jxx-1,k,ll-1)
         f0mm=f(iyym,  jxx-1,k,ll-1)
         fm0m=f(iyym-1,jxx,  k,ll-1)
         f00m=f(iyym,  jxx,  k,ll-1)
         fmm0=f(iyym-1,jxx-1,k,ll)
         f0m0=f(iyym,  jxx-1,k,ll)
         fm00=f(iyym-1,jxx,  k,ll)
         f000=f(iyym,  jxx,  k,ll)
      endif

      fdist= rweightl_*
     +    (  pweightlm*xweightl*log(max(fmmm,fmin))
     +      +pweightum*xweightl*log(max(f0mm,fmin))
     +      +pweightlm*xweightu*log(max(fm0m,fmin))
     +      +pweightum*xweightu*log(max(f00m,fmin)) )
     +      +rweightu_*
     +    (  pweightl*xweightl*log(max(fmm0,fmin))
     +      +pweightu*xweightl*log(max(f0m0,fmin))
     +      +pweightl*xweightu*log(max(fm00,fmin))
     +      +pweightu*xweightu*log(max(f000,fmin)) )
      fdist=exp(fdist)
      
      return
      end

c
