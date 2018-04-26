c
c
      subroutine tdxin33d(a,rya,klrz,expn1,expm1)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     Fill in input arrays between center and edge parabolically.
c..................................................................

      include 'param.h'


      dimension rya(0:klrz), a(0:klrz)
      data em90 /1.d-90/ 
      if (abs(a(0)) .le. em90) a(0)=em90
      dratio=a(1)/a(0)
      do 1 ll=1,klrz
        call profaxis(rn,expn1,expm1,dratio,rya(ll))
        a(ll)=a(0)*rn
 1    continue
      return
      end
