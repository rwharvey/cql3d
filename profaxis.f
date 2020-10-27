!
!
      subroutine profaxis(rn,expn1,expm1,dratio,rova)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!---------------------------------------------------------------------
!     Expands "parabolic" profiles by computing the ratio rn
!     of the local parameter value at normalized radius rova to the central 
!     value, given exponents expn1 and expm1 of the parabola, and
!     the ratio dratio of edge to central parameter value.
!---------------------------------------------------------------------
      include 'param.h'
      include 'comm.h'

      if(abs(dratio-1.).lt.1.e-12) then
        rn=1.0
      else
        rn=dratio+(1.0-dratio)*(1.-rova**expn1)**expm1
      endif
      return
      end subroutine profaxis


!======================================================================
!======================================================================
!YuP[2019-12-29] Another version, 
! where instead of dratio=e(1)/e(0) value, we use 
! both of these values (e(1)==e1 is for the plasma edge, 
! and e(0)==e0 is for the plasma center).
! This modification allows cases of e0=0.0.
! On output, use, e.g.,  a(ll)=e_out  
! [so, no need to have a(ll)=a(0)*rn ]
      subroutine profaxis1(e_out,expn1,expm1,e0,e1,rova)
      implicit none !integer (i-n), real*8 (a-h,o-z)
      save
!---------------------------------------------------------------------
!     Expands "parabolic" profiles by computing the value 
!     of the local parameter value at normalized radius rova, 
!     given exponents expn1 and expm1 of the parabola, and
!     the input values e1 at plasma edge and e0 at plasma center.
!---------------------------------------------------------------------
      !include 'param.h' !no need
      !include 'comm.h'  !no need
      real*8 expn1,expm1,e0,e1,rova ! input
      real*8 e_out ! output
      e_out= e1 + (e0-e1)*(1.d0-rova**expn1)**expm1
      !Note: in case of e0=e1, we get e_out=e1
      return
      end subroutine profaxis1
