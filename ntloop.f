c
c

c
      subroutine ntloop
      implicit integer (i-n), real*8 (a-h,o-z)
      save

      include 'param.h'
      include 'comm.h'
      include 'name.h'

      zdttot=1.0
      if (adimeth.eq."enabled" .and. transp.eq."enabled"
     +  .and. n.ge.nonadi) zdttot=2.0
      timet=timet+zdttot*dtreff
      time_(l_)=timet
      
      return
      end
