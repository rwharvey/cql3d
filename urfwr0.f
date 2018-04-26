c
c
      subroutine urfwr0(x,nrayelt,nray,nrayelts)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension x(nrayelts,*),nrayelt(*)
      save em300
      data em300 /1.d-300/
c....................................................................
c     For formatted o/p purposes, values of abs(x) .lt. 1.e-300
c     are set equal to 0. (Formats ok to 1.e-999, but real*8
c     restricted to .gt. 1.e-327 in some fortrans).
c....................................................................

      do 10  iray=1,nray
        do 11  is=1,nrayelt(iray)
          if (abs(x(is,iray)).lt.em300)  x(is,iray)=0.0
 11     continue
 10   continue

      return
      end
