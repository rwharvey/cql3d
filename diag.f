c
c
      subroutine diag
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This routine controls the calculation of various
c     density and power diagnostics.
c     In particular, it determines density gain and loss due to
c     various physical and numerical processes.
c..................................................................

      include 'param.h'
      include 'comm.h'

      if (n .gt. 0) go to 20

c..................................................................
c     compute the line averaged density at time=0 summed over all
c     general species
c..................................................................

      xlndn0(lr_)=0.
      do 10 k=1,ngen
        xlndn00(k,lr_)=xlndn(k,lr_)
        xlndn0(lr_)=xlndn0(lr_)+xlndn(k,lr_)
        do 31 iq=1,8
          sgaint(iq,k,l_)=0.
 31     continue
 10   continue
      sgaint1(l_)=0.
      consn0(l_)=xlndn0(lr_)
      return
      
      
 20   continue ! n>0
      dentot=0.

c..................................................................
c     compute cumulative in time particle density diagnostics.
c     see subs diagxswx and diagimpd.
c..................................................................

      do 40 k=1,ngen
        do 41 iq=1,8
          sgaint(iq,k,l_)=sgaint(iq,k,l_)+sgain(iq,k)
          sgaint1(l_)=sgaint1(l_)+sgain(iq,k)
          !YuP[2020] Normally, only these values are calculated:
          !sgain(3,k) from NBI or KO source
          !sgain(4,k) from flux of particles through j=jx boundary
          !sgain(5,k) from bad orbits (krook operator)
          !sgain(6,k) from toroidal losses
          !sgain(8,k) from setting negative parts of distr.func. to 0.
          ![other iq can be used in runs with implct.ne.enabled]
 41     continue
        dentot=dentot+xlndn(k,lr_)
 40   continue

c..................................................................
c     compute density conservation diagnostic - the closer to
c     1.e-14 the better
c..................................................................

      consn(l_)=(dentot-xlndn0(lr_)-sgaint1(l_)) /
     /  (.5*(xlndn0(lr_)+dentot))
     
ccc      if(l_.eq.1)then
ccc       write(*,*)'diag ',
ccc     + dentot-xlndn0(lr_),dentot-xlndn0(lr_)-sgaint1(l_),consn(l_)
ccc       write(*,'(i5,e13.5)') l_,sgain(4,1)
ccc      endif
c..................................................................
c     compute power transfer diagnostics (see sub diagentr)
c..................................................................

      do 60 k=1,ngen
        pelec=pelec+entr(k,2,l_)
        psou=psou+entr(k,5,l_)
        pwrf=pwrf+entr(k,3,l_)
 60   continue

c..................................................................
c     compute the total power being added to the system through
c     d.c. electric fields, RF power, and beam sources. (this does
c     not include power required to sustain backbround species
c     at specified fixed temperatures).
c..................................................................

      pinput=pelec+pwrf+psou
      return
      end
