c
c
      subroutine eqvolpsi(epsicon_,volum,areac)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'


c..................................................................
c     This routine returns the volume and area of a flux surface given
c     that the psi label of the flux surface is epsi. If we
c     can assume that eqorbit has been previously called for
c     this flux surface, set eqorb="disabled" before the call, otherwise
c     set eqorb="enabled". Note if eqorb is "disabled" then it is
c     assumed that the field line data exists in arrays such as
c     solr_, solz_, etc.
c..................................................................

      if (eqorb.eq."enabled") then
        rstart=0.d0 ! OUTPUT, when kopt=1
        zstart=zmag ! OUTPUT, when kopt=1
        call eqorbit(1,epsicon_,rstart,zstart) ! kopt=1
        !YuP[2020-06-30] Added kopt=2, which allows tracing surface
        ! directly from point (rstart,zstart) when it is given in INPUT
        ! (in this case value of epsicon_ is not needed).
        ! For the original design, use kopt=1,
        ! which means: find the starting point from knowledge of epsicon_ 
        ! and zmag coordinate (stored in comm.h).
      endif


c......................................................................
c     The volume (2*pi*Z*R*dR) is computed by a simple mid-point rule 
c     For the non-updown-symm case, measure volume and area to the
c     the contour from the diagonal from (rmcon_,zmcon) to (rpcon_,zpcon_).
c......................................................................

      volum=0.
      areac=0.
         do  l=2,lorbit_
            ravg=0.5*(solr_(l)+solr_(l-1))
            zdiag=zmcon_+(ravg-rmcon_)*(zpcon_-zmcon_)/(rpcon_-rmcon_)
            volum=volum+abs((solr_(l)-solr_(l-1))*
     1           (solz_(l)+solz_(l-1)-2.*zdiag)*ravg)
            areac=areac+abs((solr_(l)-solr_(l-1))*
     1           0.5*(solz_(l)+solz_(l-1)-2.*zdiag))
         enddo
c      else
c         do l=2,lorbit_
c            volum=volum+(solr_(l)-solr_(l-1))*(solz_(l)+solz_(l-1))*
c     1           (solr_(l)+solr_(l-1))
c            areac=areac+(solr_(l)-solr_(l-1))*(solz_(l)+solz_(l-1))
c         enddo
c      endif
      if (eqsym.ne."none") then ! up-down symmetry
          !-YuP: same as non-sym case, 
          ! but only half of flux surface was traced
          volum=volum*2.
          areac=areac*2.
      endif

      volum=abs(volum*pi)
      areac=abs(areac)
      
      !YuP[2015/05/02] Added a more accurate def. of volume and area
      ! near magnetic axis.
      !This definition is not dependent/sensitive to tracing the surface.
      ! Gives better convergence in eqfndpsi in cases for 
      ! radcoord.eq."sqtorflx", when 'rhonew' is proportional to 
      ! sqrt of 'volum' (near magn.axis).
      if(abs(epsicon_-psimag).le.psimag*1.d-3)then
         ! Very close to magn.axis: assume circular flux surfaces
         areac= pi*(solr_(1)-rmag)**2 ! solr_(1)-rmag == minor radius
         volum= 2*pi*rmag*areac
      endif

      return
      end
