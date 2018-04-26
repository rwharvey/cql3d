c
c
      subroutine eqflxavg(epsicon_,a,flxavg_,flxavgd_)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

      dimension a(*)
cdir$ nobounds
c..................................................................
c     This routine returns the flux surface average of a.
c     It works for both updown and non-up-down symmetric cases.
c..................................................................

      if (eqorb.eq."enabled") then
        call eqorbit(epsicon_)
      endif
      sum1=0.
      sum2=0.

c..................................................................
c     eqdell=dl; eqbpol=B-poloidal; both defined on constant phi
c     flux surface.
c..................................................................
      do 10 l=2,lorbit_
        sum1=sum1+eqdell_(l)/eqbpol_(l)
        sum2=sum2+(a(l)+a(l-1))*.5*eqdell_(l)/eqbpol_(l)
 10   continue
      flxavg_=sum2/sum1
      flxavgd_=sum1
      return
      end
