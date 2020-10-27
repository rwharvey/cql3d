!
!
      subroutine tdxin33d(a,rya,klrz,expn1,expm1)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     Fill in input arrays between center and edge parabolically.
!..................................................................

      !include 'param.h'

      dimension rya(0:klrz), a(0:klrz)
      !YuP data em90 /1.d-90/  !YuP[2019-12-29] no need anymore
      !YuP if (abs(a(0)) .le. em90) a(0)=em90 !YuP[2019-12-29] no need anymore
      !YuP dratio=a(1)/a(0)
      e0=a(0) !YuP[2019-12-29]
      e1=a(1) !YuP[2019-12-29]
      do ll=1,klrz
        !YuP call profaxis(rn,expn1,expm1,dratio,rya(ll))
        !YuP a(ll)=a(0)*rn
        call profaxis1(e_out,expn1,expm1,e0,e1,rya(ll)) !YuP[2019-12-29]
        a(ll)=e_out !YuP[2019-12-29]      
      enddo
      return
      end subroutine tdxin33d
