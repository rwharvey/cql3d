c
c
      subroutine achiefn(kopt)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     kopt=0: This routine advances the  equations in momemtum space 
c             over time-step n (or if cqlpmod.eq.'disabled',
c             transp='enabled',soln_method.eq.'it3drv', then only
c             calc coeffs until at last flux surface.
c     kopt=1: Compute plasma density and energy transfer
c             based on new f (after solution is obtained); 
c             when kopt=1, impavnc0 is not called.
c     kopt=2: Calculated h and g functions and associated fluxes,
c             to obtain electrostatic electric field for
c             cqlpmod="enabled",sbdry="periodic",esfld="enabled".
c             (BH:Evidently, have not fully implemented this yet. 
c             No kopt reference is in impanvc0.  It was
c             successfully implemented in the STELLA code (a derivative
c             of Kupfer's FPET code).  See achiefn in stella.f. 
c             Also check file Harvey_rf99.pdf, or related paper, from 
c             Radiofrequency Power in Plasmas meeting, Anapolis, 1999.)
c     kopt=3: Calculate h and g functions, for time-advancing of electric
c             field according to Ampere-Faraday equations.
c             Iteration on Kupfer's g and h functions is used in
c             obtaining the implicit toroidal electric field.
c             BH131104.
c             
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE     
 
      character cptline*80
      integer getpid
      real :: cputime
            

      if (kopt.eq.0  .or.  kopt.eq.2) then
            
c..................................................
c     bring background profiles up to time step n
c..................................................
ccc      if(nefiter.eq.1) call profiles !-YuP: moved to tdchief

c.....................................................................
c     compute electric field either according to direct specification, 
c     or to maintain a target toroidal current.
c.....................................................................
      if (eseswtch.ne."enabled") then
         call efield
      endif

        ! Start time advancement:
      if(nefiter.eq.1) then
      if (kopt.eq.0) then
           n=n+1
           n_(l_)=n ! new time-step for this flux surface
           ! for 2-d (v_par,v_perp) calculation ntloop controls
           ! end of run or restart.
           ! Also updates time.
           call ntloop
      endif
      endif
      
c..................................................
c     determine coulomb log...
c..................................................
      if(nefiter.eq.1) call cfpgamma

c....................................................
c     Apply additional preloading to f, if n=nfpld
c     (If only initial loading is used, then nfpld=0)
c....................................................
      if (n.eq.nfpld .and. nefiter.eq.1) then
          call finit ! Note: at n=0 finit is called from ainitial.
      endif          ! First call to achiefn is with n=0, then advanced.

c..................................................
c     determine particle sources...
c..................................................
      if (n.ne.1 .and. nefiter.eq.1) call sourcee
      ! Note: at n=0 sourcee is called from ainitial
      ! First call to achiefn is with n=0, then advanced in tdchief.
      

CMPIINSERT_IF_RANK_EQ_MPIWORKER
c..................................................
c     generate collisional F.P. coefficients...
c..................................................
       call cfpcoefc ! time-consuming calculations
c..................................................................
c     call implicit time advancement routine if implct.eq."enabled"
c..................................................................
      if (implct .eq. "enabled") call impavnc0(kopt)
CMPIINSERT_ENDIF_RANK

c...........................................................
c     call splitting time advancement if implct = "disabled"
c...........................................................
      if (implct .ne. "enabled") call exsweep
      
      go to 999 !-> return/end

      endif ! kopt=0 or 2
      
c..................................................................
c     Simultaneous FP and Ampere-Faraday Eqns solution, kopt=3 
c     Following above coding for kopt=0.   
c..................................................................
      if (kopt.eq.3) then
         n=n+1  ! n is used by logics in cfpcoefn, impavnc0, etc.
         zdttot=one
         call cfpgamma
         if (adimeth.eq."enabled" .and. transp.eq."enabled"
     +        .and. n.ge.nonadi) zdttot=2.0
         timet=timet+zdttot*dtreff
        ! n and time are not really updated yet - 
        ! n and time are updated after n_ and time_ are updated.
        call cfpgamma  !determine coulomb log...   But called 6 lines above?
cBH131107:  NEED to check this n=0 issue.
        ! Apply additional preloading to f, if n=nfpld
        ! (If only initial loading is used, then nfpld=0)
        if (n.eq.nfpld) then
               call finit
        endif
        ! Note: at n=0 finit is called from ainitial.
        ! First call to achiefn is with n=0, then advanced.
        if (n.ne.1) call sourcee !determine particle sources...
        ! Note: at n=0 sourcee is called from ainitial
        ! First call to achiefn is with n=0, then advanced in tdchief.

c For ZOW, parallelization (in lr) is done here.
CMPIINSERT_IF_RANK_EQ_MPIWORKER
          !collisional F.P. coefficients:
          call cfpcoefc ! time-consuming calculations
          !call implicit time advancement:
          if (implct .eq. "enabled") call impavnc0(kopt)
CMPIINSERT_ENDIF_RANK
cBH110309:  Restore n
cBH110309:        n=n-1
      endif  !On kopt.eq.3




      



      if (kopt.eq.1) then 
c..................................................
c     Compute plasma energy and density and energy transfer
c     rate terms...
c     Diagnostics ONLY COMPUTED for cases where transport
c     is not considered (transp="disabled")
c..................................................

      if (transp.eq."disabled" .or. 
     +  (cqlpmod.eq."enabled".and.n.ge.nontran)) then
        call diaggnde
c..................................................
c     Compute plasma resistivity for electron runs..
c..................................................
        call restvty
c..................................................
c     Compute various diagnostics...
c..................................................
        call diag
c..................................................
c     Obtain data for time dependent plots
c..................................................
        call ntdstore
      endif

c..................................................
c     plotting logic for 2-D (v-theta) code...
c..................................................
      iplot=0
      if (noplots.ne."enabled1") then
         do i=1,nplota
            if(n.eq.nplot(i)) then
               iplot=i
               tplot(i)=timet
            endif
         enddo
      endif

      if (lrzmax.eq.1) then
        if (n.eq.nstop .or. iplot.ne.0) call tdoutput(2)
        if (n.eq.nstop .or. iplot.ne.0) call pltmain
        if (n.eq.nstop.and.pltra.eq."enabled") then
           call pltrun
        endif
ccc        if (netcdfnm.ne."disabled" .and. nstop.ne.0) then
ccc           call netcdfrw2(1) !YuP: Why here? Also called from tdchief
ccc        endif
        if(n.eq.nstop) then
           call dskout(l_)
           call pgend
cBH080106           cputime=etime(tarray)
           call cpu_time(cputime)
           write (*,*) 'achiefn/lrzmax=1: CPU time (seconds)', cputime
           if (ichkpnt.ne."disabled") then
             write(cptline,100) i,ichkpnt
 100         format("chkpnt -p ",i5," -f ",a8)
             print *,cptline
           endif
        endif ! n.eq.nstop
      endif ! lrzmax=1

c..................................................
c     plotting (and netcdf store) logic for 3-D code ..
c     Used if transp="disabled".
c..................................................
      if (transp .eq. "disabled" .and. noplots.ne."enabled1") then
c$$$      if (transp .eq. "disabled") then
         if (lrzmax.gt.1) then
            if (n.ge.nstop) then
               call pltmain
               call netcdfmain
               go to 999 !-> return/end
            endif
            if (iplot.ne.0) call pltmain
         endif
      endif
      
      go to 999 !-> return/end
      
      endif ! if kopt=1
      
      
     

 999  return

      end
