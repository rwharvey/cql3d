

      subroutine ainitial
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c............................................................
c     initializes the driver.
c............................................................

      include 'param.h'
      include 'comm.h'
      include 'name.h'
CMPIINSERT_INCLUDE     
c.......................................................................

      if (n.gt.0) return

c...........................................................
c     call a routine which calculates mesh integration
c     coefficients
c YuP[Sept-2014] moved this part outside of ainitial: 
c needed for proper work of FOW-logic.
c See notes in tdinitl, after if(fow.ne.'disabled')...
cyup      call micxinil
c...........................................................

c.....................................................................
c     call orbit bounce time routine.  Also orbit shift calc.
c YuP[Sept-2014] moved this part outside of ainitial: 
c needed for proper work of FOW-logic.
c See notes in tdinitl, after if(fow.ne.'disabled')...
c      if (l_ .eq. lmdpln_) then
c         if(taunew.eq."enabled") then
cyup            call baviorbt
c         else
cyup            call baviorbto
c         endif
c      endif
c.....................................................................

c.....................................................................
c     bounce average various geometrical quantities..
c.....................................................................

      if (l_ .eq. lmdpln_) call bavgmax

c.......................................................................
c     Redefine vptb=1.0 if CQLP
c.......................................................................

      if (cqlpmod.eq."enabled" .and. l_.eq.lmdpln_) call wpvptb

c..............................................................
c     Evaluate an additional integration coefficient.
c..............................................................

      if (l_ .eq. lmdpln_) call micxinim

c...........................................................
c     call routine which initializes distribution functions.
c...........................................................

      call finit
c.....................................................................
c     compute the orbit loss krook operator and the toroidal
c     loss operator
c.....................................................................
cyup [07-14-2014] Moved losscone in front of diaggnde, so that 
c                 gone() array is available in diaggnde, at n=0.
      if (l_ .eq. lmdpln_) then
         call losscone
      endif

c...........................................................
c     call routine to calculate energies and densities.
c...........................................................

      call diaggnde

c.....................................................................
c     compute the orbit loss krook operator and the toroidal
c     loss operator
c.....................................................................
cyup      if (l_ .eq. lmdpln_) call losscone

c.....................................................................
c     Compute coefficients for synchrotron radiation loss
c.....................................................................

      if (l_ .eq. lmdpln_) call synchrad

c...................................................................
c     Compute coefficients for phenomenological energy loss
c...................................................................

      if (l_ .eq. lmdpln_) call lossegy

c.....................................................................
c     d.c. ohmic electric field coefficient set-up
c.....................................................................

      call coefefld

c.....................................................................
c     Call plasma resistivity diagnostic...
c.....................................................................

      call efield
      call restvty

c...........................................................
c     call driver routine which initializes analytic sources.
c...........................................................

      if (l_ .eq. lmdpln_) call sourcee

c.......................................................................
c     Call RF module if this is a 2-D calculation and if the RF
c     module is appended.
c.......................................................................

      if (lrzmax.eq.1 .or.(cqlpmod.eq."enabled".and.l_.eq.lmdpln_)) then
         if (nrf.ge.1) then
            if (vlfmod.eq."enabled") then
               call vlf("setup")
            elseif (vlhmod.eq."enabled") then
               call vlh("setup")
            else
c     call rf("setup")
               do 850 ku=1,ngen
                  call diagentr(3,ku)
 850           continue
            endif
         endif
      endif

      if (lrzmax.gt.1 .and. nrf.ge.1 .and. vlhmod.eq."enabled") 
     +     call vlh("setup")
      
      if (lrzmax.gt.1 .and. nrf.ge.1 .and. vlfmod.eq."enabled") 
     +     call vlf("setup") ! Added YuP[03-2016] Calling vlf for lrz>1
     
c.....................................................................
c     Call conservation diagnostic routine..
c.....................................................................

      call diag

c.....................................................................
c     Call time plot storage routine...
c.....................................................................

      call ntdstore

c.....................................................................
c     Call plotting routine...
c.....................................................................

      call pltmain

c.....................................................................
c     print out some parameters
c.......................................................................

      if (lrzmax .eq. 1) call tdoutput(1)
C

      return
      end
