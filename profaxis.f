c
c
      subroutine profaxis(rn,expn1,expm1,dratio,rova)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c---------------------------------------------------------------------
c     Expands "parabolic" profiles by computing the ratio rn
c     of the local parameter value at normalized radius rova to the central 
c     value, given exponents expn1 and expm1 of the parabola, and
c     the ratio dratio of edge to central parameter value.
c---------------------------------------------------------------------
      include 'param.h'
      include 'comm.h'

      if(abs(dratio-1.).lt.1.e-12) then
        rn=1.0
      else
        rn=dratio+(1.0-dratio)*(1.-rova**expn1)**expm1
      endif
      return
      end
