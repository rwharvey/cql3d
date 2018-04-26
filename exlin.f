c
c
      real*8 function exlin(y1,y2,x1,x2,xx)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
cmnt  linear interpolation.
c
      exlin=(y2-y1)/(x2-x1)*(xx-x1)+y1
      return
      end
