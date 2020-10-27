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
        rstart=0.d0 ! OUTPUT, when kopt=1
        zstart=zmag ! OUTPUT, when kopt=1
        call eqorbit(1,epsicon_,rstart,zstart) ! kopt=1
        !YuP[2020-06-30] Added kopt=2, which allows tracing surface
        ! directly from point (rstart,zstart) when it is given in INPUT
        ! (in this case value of epsicon_ is not needed).
        ! For the original design, use kopt=1,
        ! which means: find the starting point from knowledge of epsicon_ 
        ! and zmag coordinate (stored in comm.h).
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
