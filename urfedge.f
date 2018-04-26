c
c
      subroutine urfedge(i1,i2,c1,c2,lp,lrad,jval)
      implicit integer (i-n), real*8 (a-h,o-z)

c..........................................................
c     Computes weights for diffusion coefficients if spectral
c     width of ray element is small.
c..........................................................

      include 'param.h'
      include 'comm.h'

      if (c2 .lt. -1.) c2=-1.
      if (c2 .gt. 1.) c2=1.
      if (c1 .lt. -1.) c1=-1.
      if (c1 .gt. 1.)  c1=1.  

      if (i1.ne.i2) then

        if (i1.ne.1) then
          fact1=abs(cosmz(i1-1,lp,lrad)-cosmz(i1,lp,lrad))
          del1=abs(c1-cosmz(i1-1,lp,lrad))
        else           !Not possible
          fact1=abs(-1.-cosmz(i1,lp,lrad))
          del1=abs(-1.-c1)
        endif
        if (i2.ne.1) then
          del2=abs(c2-cosmz(i2,lp,lrad))
          fact2=abs(cosmz(i2,lp,lrad)-cosmz(i2-1,lp,lrad))
        else
c*bh*940227del2=abs(c2-cosmz(1,lp,lrad))
c*bh*940227fact2=abs(-1.-cosmz(1,lp,lrad))!closetozero.
          del2=1.0
          fact2=1.0
        endif
        scal1=del1/fact1
        scal2=del2/fact2
      else
        del1=abs(c2-c1)
        if (i1.ne.1) then
          fact1=abs(cosmz(i1-1,lp,lrad)-cosmz(i1,lp,lrad))
        else             !Not possible
          fact1=abs(-1.-cosmz(i1,lp,lrad))
        endif
        scal1=del1/fact1
        scal2=scal1
      endif
      ifct1(jval)=65535*scal1+.5
      ifct2(jval)=65535*scal2+.5
      return
      end
