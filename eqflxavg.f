!
!
      subroutine eqflxavg(epsicon_,a,flxavg_,flxavgd_)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

      dimension a(*)
!cdir$ nobounds
!..................................................................
!     This routine returns the flux surface average of a.
!     It works for both updown and non-up-down symmetric cases.
!..................................................................

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
      sum1=0.
      sum2=0.

!..................................................................
!     eqdell=dl; eqbpol=B-poloidal; both defined on constant phi
!     flux surface.
!..................................................................
      do l=2,lorbit_
        sum1=sum1+eqdell_(l)/eqbpol_(l)
        sum2=sum2+(a(l)+a(l-1))*.5*eqdell_(l)/eqbpol_(l)
      enddo
      flxavg_=sum2/sum1
      flxavgd_=sum1
      return
      end subroutine eqflxavg


!==================================================================
!==================================================================
! YuP[2019-12-12] migrated this subroutine from CQL3D-FOW version

      subroutine eqflxavg_lz(ir,a,flxavg_,sum_a_ds_bp)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
      include 'name.h'

      dimension a(*)
!..................................................................
!     This routine returns the flux surface average of a.
!     It works for both updown and non-up-down symmetric cases.
!..................................................................
      ! output: <a> = sum(a*ds/Bpol)/sum(ds/Bpol) 
      ! Note: If no 1/Bpol in averaging, it results in ~1% increase  
      ! of <a> at rho~0.5, 7% at rho~0.8, 16% at rho~1.0, 
      ! comparing to sum(a*ds/Bpol)/sum(ds/Bpol).
      ! Note also: the number of points for averaging is lz here,
      ! typically 30. This is sufficient: using 150 points yields 
      ! about same result (up to 4th digit).
      sum1=0.
      sum2=0.
      
      if(eqmod.eq."enabled")then
      
        do l=2,lz
        ! arc-length in poloidal direction:
        dsl=sqrt((solrz(l,ir)-solrz(l-1,ir))**2+
     &           (solzz(l,ir)-solzz(l-1,ir))**2 )  
        ! Include 1/Bpol factor:
        bpol= 0.5*(bpolz(l,ir)+bpolz(l-1,ir))
        dsl=dsl/bpol
        sum1=sum1+dsl       
        sum2=sum2+(a(l)+a(l-1))*.5*dsl
        enddo
 
      else ! eqmod.eq."disabled" (circular-shaped plasma eq. model)
        ! bpolz is not available in this case
        
        do l=2,lz
        ! arc-length in poloidal direction:
        dsl=sqrt((solrz(l,ir)-solrz(l-1,ir))**2+
     &           (solzz(l,ir)-solzz(l-1,ir))**2 )  
        ! Include 1/Bpol factor:
        rloc=0.5*(solrz(l,ir)+solrz(l-1,ir)) ! major radius at local point
        ! bthr(lr) is the pol.field without Rmag/R factor
        bpol= bthr(ir)*radmaj/rloc 
        dsl=dsl/bpol
        sum1=sum1+dsl       
        sum2=sum2+(a(l)+a(l-1))*.5*dsl
        enddo
        
      endif
 
      flxavg_=sum2/sum1
      sum_a_ds_bp=sum2
      return
      end subroutine eqflxavg_lz
