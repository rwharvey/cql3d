!
!
      subroutine tdxin23d(a,rya,klrz,ngn,nso,k,kk,expn1,expm1)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     Fill in input arrays between center and edge parabolically.
!..................................................................

      !include 'param.h'

      dimension rya(0:klrz), a(ngn,nso,0:klrz)
      ! data em90 /1.d-90/  !YuP[2019-12-29] no need anymore
      ! if (abs(a(k,kk,0)) .le. em90) a(k,kk,0)=em90 !YuP[2019-12-29] no need anymore
      !YuP dratio=a(k,kk,1)/a(k,kk,0)
      e0=a(k,kk,0) !YuP[2019-12-29]
      e1=a(k,kk,1) !YuP[2019-12-29]
      do ll=1,klrz
        !YuP call profaxis(rn,expn1,expm1,dratio,rya(ll))
        !YuP a(k,kk,ll)=a(k,kk,0)*rn
        call profaxis1(e_out,expn1,expm1,e0,e1,rya(ll)) !YuP[2019-12-29]
        a(k,kk,ll)=e_out !YuP[2019-12-29]
      enddo
      return
      end subroutine tdxin23d
