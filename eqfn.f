c
c
      real*8 function eqfn(e,scalfct)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'


c..................................................................
c     This routine returns ad-hoc psi as a function of e.
c..................................................................

      if (eqmodel.eq."power") then
        if (e.eq.zero) then
          eqfn=one
        else
          eqfn=scalfct*e**eqpower
        endif
      else
        call eqwrng(5)
      endif
      return
      end
