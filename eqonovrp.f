c
c
      subroutine eqonovrp(epsicon_,onovrp1,onovrp2)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'


c..................................................................
c     This routine computes <1./R**2>, <1./R> for flux surface psival
c..................................................................

      if (eqorb.eq."enabled") then
        call eqorbit(epsicon_)
      endif
      do 20 ipower=1,2
        do 10 l=1,lorbit_
          tlorb1(l)=1./(solr_(l)**ipower)
 10     continue
        call eqflxavg(epsicon_,tlorb1,onovrs,flxavgd_)
        if (ipower.eq.1) onovrp1=onovrs
        if (ipower.eq.2) onovrp2=onovrs
 20   continue
      return
      end
