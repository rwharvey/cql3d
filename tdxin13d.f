
!
!
      subroutine tdxin13d(a,rya,klrz,m,k,expn1,expm1)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!.............................................................
!     this is a utility parabolic "fill in" routine
!.............................................................

      !include 'param.h'
      dimension rya(0:klrz), a(m,0:klrz)
      ! data em90 /1.d-90/   !YuP[2019-12-29] no need anymore
      ! if (abs(a(k,0)).le. em90) a(k,0)=em90  !YuP[2019-12-29] no need anymore
      !YuP dratio=a(k,1)/a(k,0)
      e0=a(k,0) !YuP[2019-12-29]
      e1=a(k,1) !YuP[2019-12-29]
      do ll=1,klrz
        !YuP call profaxis(rn,expn1,expm1,dratio,rya(ll))
        !YuP a(k,ll)=a(k,0)*rn
        call profaxis1(e_out,expn1,expm1,e0,e1,rya(ll)) !YuP[2019-12-29]
        a(k,ll)=e_out !YuP[2019-12-29]
      enddo
      return
      end subroutine tdxin13d
