c
c
      subroutine eqfpsi(psval,fpsi__,fppsi__)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
      character*8 eqmirror

c..................................................................
c     This routine provides f(psi) to model the toroidal
c     magnetic field. For cases that eqsource="ellipse"
c     the f is ad-hoc and is determined through the namelist
c     model, fpsimodel. In the case that eqsource="eqdsk", then a file
c     (eqdskin) exists on disk which provides f(psi) and the 
c     equilibrium psi. As of 9/21/88 eqsource=eqdsk or topeol.
c     As if March/2016, eqsource has five possible values.
!     (eqsource.eq."mirror1" is not accessible in this CQL3D version)
c     Also provided is the derivative df/dpsi, fppsi.
c..................................................................

      eqmirror="disabled"
      if( eqsource.eq."eqdsk".and.machine.eq."mirror") then
         eqmirror="enabled"
      else
         eqmirror="disabled"
      endif
      if ( (eqsource.eq."ellipse")  .or. 
     +     (eqsource.eq."miller")   .or.
     +     (eqsource.eq."mirror1" .or. eqmirror.eq."enabled") ) then
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
      character*8 eqmirror

c..................................................................
c     This routine provides p(psi) to model the plasma
c     pressure. For cases that eqsource="ellipse"
c     the p i zero. In the case that eqsource="filename", then
c     a file exists on disk which provides f and the equilibrium
c     psi. As of 9/21/88 filename=eqdsk or topeol.
c     Also provided is the derivative dp/dpsi, pppsi.
c..................................................................

      eqmirror="disabled"
      if( eqsource.eq."eqdsk".and.machine.eq."mirror") then
         eqmirror="enabled"
      else
         eqmirror="disabled"
      endif
      if ( (eqsource.eq."ellipse")  .or. 
     +     (eqsource.eq."miller")   .or.
     +     (eqsource.eq."mirror1" .or. eqmirror.eq."enabled") ) then
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
