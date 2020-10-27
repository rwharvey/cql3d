
c
c
      subroutine urfchief
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This routine controls the urf module
c     The following control variables are used in determining the
c     action of this subroutine:
c     nrfstep1, nrfstep2, nrfpwr, nrfitr1, nrfitr2, nrfitr3, urfncoef
c
c     nurf=a counter for the number of calls to urfchief which have 
c     resulted in calculation or recalculation of the diffusion 
c     coefficients. (Initialized to 0).
c     The nurf variable is incremented (at end of urfchief) each time
c     n/integer(urfncoef*ncoef)*integer(urfncoef*ncoef).eq.n.
c
c     The sequence of actions (as a function of
c     increasing nurf.ge.0 at each call) is given by the following steps
c     (after each diffusion coeff calc, control returns to the calling
c     subroutine):
c     1. For nurf=0, calc or read ray data for nrfstep1 steps along ray.
c     2. Calc. damping of ray data, and then resulting quasilinear
c     diffusion coeffs, using a fraction of the input power
c     = (1/2)**nrfpwr.
c     return.
c     3. Repeat step 2 for next nrfpwr calls, but with fractional input
c     power (1/2)**(nrfpwr-1), (1/2)**(nrfpwr-2),.... (1/2)**0.
c     (This step is a no-op if nrfpwr=0).
c     4. Iterate step 2 with full input power for next nrfitr1 calls.
c     (This step is a no-op if nrfitr1=0).
c     5. Extend extendable rays by nrfstep2 steps.
c     6. Re-calc damping from ray data and then quasilinear diffusion
c     coeffs.  Iterate this step nrfitr2 additional calls.
c     7. Steps 5 and 6 are carried out nrfitr3 times.
c
c     Thus choose
c     nstop=
c     (nrfpwr+1+nrfitr1+nrfitr3*(nrfitr2+1))*integer[urfncoef*ncoef]
c     if the above sequence is to be completed.
c
c     This urf module treats multiple wave types, multiple cyclotron 
c     harmonics for each wave type, and there may be multiple
c     rf diffused general (FP'd) species.   We refer to the sum
c     of the wave harmonics over the wave types, as the wave modes,
c     that is, one wave mode for each harmonic of each wave type.
c     The following variables are used in the urf*.f routines to
c     direct the logic:
c        mrf= number of rf "types"
c        mrfn= number of rf "modes" (sum over mrf of the nharms())
c        irfn(1:mrf)=  rf mode index (in 1:mrfn) of the lowest
c                      harmonic for each wave type
c        krfn(1:mrfn)=  wave type index (in 1:mrf) for each rf mode
c        nharm1(1:mrf)= lowest cyclotron harmonic, for each rf type
c        nharms(1:mrf)= number of cyclotron harmonics, for each rf type
c        nharm(1:mrfn)= harmonic number for each "mode".
c        nrfspecies(1:mrf) = general species index which each rf
c                            type is applied to (a nml variable)
c
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

c...................................................................
c     Control logic
c...................................................................
      
      imprf=0
      iurfncf=urfncoef*ncoef
      if(n/iurfncf*iurfncf.ne.n)  return

c.......................................................................
c     urfrstrt.eq."enabled"  => "do not update delpwr" option
c     For convergence studies with previously calculated ray data files...
c.......................................................................

      iopt=1
      if (urfrstrt.eq."enabled") iopt=2
      
cpu      call cpu_time(t_urf1) !-.-.-.-.-.-.-.-.-.-.-.-.-.

c...................................................................
c     Get a time average of the distribution function for purposes
c     of computing wave deposition and damping.
c...................................................................
      call urfavg
      
      
      if(nurf.eq.0)  then

c...................................................................
c     Initialize counters for control logic
c...................................................................

        irfpwr=nrfpwr
        irfitr1=nrfitr1
        irfitr2=0
        irfitr3=nrfitr3

c...................................................................
c     Calculate the ray data, depending on call_lh,call_ech,call_fw
c...................................................................

        initrf=0
cFollowing needs reconfiguring for arbitrary nmodsa [BH, 060314]
c$$$        nraypts1=nrfstep1(1)
c$$$        nraypts2=nrfstep1(2)
c$$$        nraypts3=nrfstep1(3)
c$$$        call urfrays(initrf,nraypts1,nraypts2,nraypts3)

c...................................................................
c     Read in ray data
c...................................................................
        if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
        write(*,*)'call urfread: n, nurf=',n,nurf
        endif
        call urfread

c...................................................................
c     Obtain segregation of ray data by flux surface
c...................................................................

        call urffflx

c...................................................................
c     Save initial power in the rays
c...................................................................

        do 10 krf=1,mrf
          kk=irfn(krf)       !irfn setup in urfread
          do 11  iray=1,nray(kk)
            delpwr0(iray,kk)=delpwr(1,iray,kk)
 11       continue
 10     continue

c...................................................................
c     Initialize Bessel function table
c...................................................................

        call urfbes

c...................................................................
c     Initialize packed data and do damping
c...................................................................

        call urfpack
             
        call urfdamp0(irfpwr,iopt)       
       
      elseif (iopt .eq. 2) then

c.......................................................................
c     compute damping and power but do not update delpwr according
c     to what has been absorbed before => compatible with lrz<lrzmax
c.......................................................................

        call urfdamp0(irfpwr,iopt)

      elseif (irfpwr.gt.0)  then
        irfpwr=irfpwr-1

c...................................................................
c     do damping
c...................................................................

        call urfdamp0(irfpwr,iopt)

      elseif (irfitr1.gt.0)  then
        irfitr1=irfitr1-1

c...................................................................
c     do damping
c...................................................................

cpu      call cpu_time(t_urf2) !-.-.-.-.-.-.-.-.-.-.-.-.-.
cpu      write(*,*)'urfchief t2=', t_urf2-t_urf1
      
        call urfdamp0(0,iopt)
        
cpu      call cpu_time(t_urf3) !-.-.-.-.-.-.-.-.-.-.-.-.-.
cpu      write(*,*)'urfchief t3=', t_urf3-t_urf1

      elseif (irfitr3.gt.0 .and. irfitr2.eq.0)  then
        irfitr3=irfitr3-1
        irfitr2=nrfitr2

c...................................................................
c     Extend rays
c...................................................................

        initrf=1
        nraypts1=nraypts1+nrfstep2
        nraypts2=nraypts2+nrfstep2
        nraypts3=nraypts3+nrfstep2
        call urfwrite
        call urfrays(initrf,nraypts1,nraypts2,nraypts3)
        if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
           write(*,*)'call urfread: n, nurf,irfitr2,irfitr3='
     +         ,n,nurf,irfitr2,irfitr3
        endif
        call urfread

c...................................................................
c     Obtain re-segregation by flux surface of ray data
c...................................................................

        call urffflx

c...................................................................
c     Extend packed data and do damping
c...................................................................

        call urfpack
        call urfdamp0(0,iopt)

      elseif(irfitr3.ge.0 .and. nrfitr3.gt.0 .and. irfitr2.gt.0)  then
        irfitr2=irfitr2-1

c...................................................................
c     do damping
c...................................................................

        call urfdamp0(0,iopt)
      endif


c..................................................................
c     Call the coefficient generator. 
c     Radial surfaces: lrz (and not lrors)
c     Generate B,C,E and F
c..................................................................
      call cpu_time(t_urf1) !-.-.-.-.-.-.-.-.-.-.-.-.-.
      call urfb0  ! YuP-110222: Now includes all lr_ internally 
      call cpu_time(t_urf2) !-.-.-.-.-.-.-.-.-.-.-.-.-.
      
CMPIINSERT_IF_RANK_EQ_0
      WRITE(*,'(a,f10.3)')'tcpu_URFB0=',t_urf2-t_urf1
CMPIINSERT_ENDIF_RANK

      imprf=1


      do 50 ll=1,lrz
c..................................................................
c     Read flux surface dependent data from disk.
c     for midplane parameter values
c..................................................................
        call tdnflxs(lmdpln(ll))
c%os  
        if (indxlr_ .ne. ll) then
CMPIINSERT_IF_RANK_EQ_0
          PRINT *,' error in defining l_,indxlr_,... in urfchief'
CMPIINSERT_ENDIF_RANK  
          stop 'urfchief: error in defining l_,indxlr_,...'
        endif
c%os  
c..................................................................
c     Plot out diffusion coefficients
c..................................................................
        if (noplots.ne."enabled1") call urfbplt !only when plturfb="enabled"
                                                !or  plturfb='color'
 50   continue

      nurf=nurf+1

      return
      end
