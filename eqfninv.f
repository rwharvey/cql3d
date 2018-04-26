
c
c
      real*8 function eqfninv(psival,scalfct)
      implicit integer (i-n), real*8 (a-h,o-z)

      include 'param.h'
      include 'comm.h'


c..................................................................
c     This routine returns E given psival..
c..................................................................

      if (eqmodel.eq."power") then
        eqfninv=(psival/scalfct)**(1./eqpower)
      else
        call eqwrng(6)
      endif
      return
      end
