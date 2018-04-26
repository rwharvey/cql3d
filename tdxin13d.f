
c
c
      subroutine tdxin13d(a,rya,klrz,m,k,expn1,expm1)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c.............................................................
c     this is a utility parabolic "fill in" routine
c.............................................................

      include 'param.h'
      dimension rya(0:klrz), a(m,0:klrz)
      data em90 /1.d-90/ 
      if (abs(a(k,0)).le. em90) a(k,0)=em90
      dratio=a(k,1)/a(k,0)
      do 1 ll=1,klrz
        call profaxis(rn,expn1,expm1,dratio,rya(ll))
        a(k,ll)=a(k,0)*rn
 1    continue
      return
      end
