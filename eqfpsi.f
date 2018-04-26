c
c
      subroutine eqfpsi(psval,fpsi__,fppsi__)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

c..................................................................
c     This routine provides f(psi) to model the toroidal
c     magnetic field. For cases that eqsource="ellipse"
c     the f is ad-hoc and is determined through the namelist
c     model, fpsimodel. In the case that eqsource="filename", then
c     a file exists on disk which provides f and the equilibrium
c     psi. As of 9/21/88 filename=eqdsk or topeol.
c     Also provided is the derivative df/dpsi, fppsi.
c..................................................................

      if (eqsource.eq."ellipse") then
        if (fpsimodl.eq."constant") then
          fpsi__=btor*radmaj
          fppsi__=0.
        else
          call eqwrng(7)
        endif
      else
        itab(1)=1
        itab(2)=1
        itab(3)=0
        call terp1(nfp,psiar,fpsiar,d2fpsiar,psval,1,tab,itab)
        fpsi__=tab(1)
        fppsi__=tab(2)
      endif
      return
      end


c
c
      subroutine eqppsi(psval,ppsi__,pppsi__)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

c..................................................................
c     This routine provides p(psi) to model the plasma
c     pressure. For cases that eqsource="ellipse"
c     the p i zero. In the case that eqsource="filename", then
c     a file exists on disk which provides f and the equilibrium
c     psi. As of 9/21/88 filename=eqdsk or topeol.
c     Also provided is the derivative dp/dpsi, pppsi.
c..................................................................

      if (eqsource.eq."ellipse") then
         ppsi__=0.
         pppsi__=0.
      else
        itab(1)=1
        itab(2)=1
        itab(3)=0
        call terp1(nfp,psiar,prar,d2prar,psval,1,tab,itab)
        ppsi__=tab(1)
        pppsi__=tab(2)
      endif
      return
      end
