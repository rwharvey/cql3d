
      subroutine eqcoord
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
      character*8 generate


c..................................................................
c     This routine controls the translation of the equilibrium psi
c     and f =(R*Btor) data into information specific to a flux surface
c     characterized by radial coordinate rovera(lr_). We have 
c     erhocon(lr_)=rovera(lr_)*rhomax where rhomax=eqrho(nconteqn).
c     For a circular flux surface rhomax would correspond to radmin.
C     !!!! for small aspect ratio only  !!!!
c..................................................................

c..................................................................
c     Generate an array eqrho(j), j=1,..., nconteqn. This array
c     contains the toroidal radial coordinate, rho, corresponding
c     to the psi array eqpsi(j). This will be used later as a
c     basis for spline interpolation.
c..................................................................
      generate="disabled"
      if (lr_.eq.lrzmax) then
        generate="enabled"
        call eqrhopsi(generate)
      endif

c..................................................................
c     Next determine the radial coordinate of interest.
c..................................................................

      if (rovera(lr_).ge.0) then
        erhocon(lr_)=rovera(lr_)*rhomax
      endif

c..................................................................
c     Determine the psi value (epsicon(lr_)) corresponding to 
c     erhocon(lr_) (erhocon(lr_) is passed in a common block).
c..................................................................

      call eqfndpsi(psides,areades,volum)
      epsicon(lr_)=psides
      areacon(lr_)=areades
      volcon(lr_)=volum
C-----------------------------------------------------------------------
c     Determine some geometrical factors needed to calculate the aspect
c     ratio and the elongation
C-----------------------------------------------------------------------
      rgeom(lr_) =0.5 * (rpcon(lr_) - rmcon(lr_))
      r0geom(lr_)=0.5 * (rpcon(lr_) + rmcon(lr_))
      call aminmx(solz_,1,lorbit(lr_),1,zmincon1,zmaxcon1,kmin,kmax)
      zgeom(lr_)=zmaxcon1

c..................................................................
c     If the density and temperature profiles are to have a strict
c     psi dependence, do that here. This will overwrite the default
c     rho profile determined in subroutine tdxinit
c..................................................................

      if (profpsi.eq."enabled".and. lrzmax.gt.1) then
        do 1000 k=1,ntotal
          reden(k,lr_)=reden(k,0)*(epsicon(lr_)/psimag)**.5
          temp(k,lr_)=(epsicon(lr_)/psimag)*temp(k,0)
 1000   continue
      endif
      return
      end
