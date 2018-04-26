c
c
      subroutine tdtrcon
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c..............................................................
c     Compute conservation constant.
c..............................................................
      include 'param.h'
      include 'comm.h'
      data sgainr /0.0/
c jakub urban 110708: commented out for g95 compiler
c the above blanket save statement should do the job
c      save sgainr

c..............................................................
c     Compute original number of particles in tokamak.
c..............................................................
      if (n.eq.1) then
        total0=0.
        do 10 l=1,lrz
          ilr=lrindx(l)
          total0=xlndn0(ilr)/zmaxpsi(ilr)*dvol(ilr)+total0
 10     continue
      endif

c..............................................................
c     Compute total number of particles in device now.
c...............................................................

      total=0.
      do 30 l=1,lrz
        ilr=lrindx(l)
        do 40 k=1,ngen
          total=xlndn(k,ilr)/zmaxpsi(ilr)*dvol(ilr)+total
 40     continue
 30   continue

c................................................................
c     Call routine to compute number of particles lost at limiter
c     this time step.
c................................................................

      call tdtrflx
      sgainr=sgainr+flxout
      conserv=(total-total0-sgainr)/(total*.5+total0*.5)
      return
      end
