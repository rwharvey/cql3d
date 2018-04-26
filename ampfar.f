
      double precision function ampfarl(distn,ll)
      implicit integer (i-n), real*8 (a-h,o-z)

c     Returns result of an integral of the distribution function
c     distn, arising in the Ampere-Faraday eqns.
c     distn = combinations of h and g,  where h and g are similar
c     to the iterative solution introduced in 
c     ``Kinetic Modeling of SOL Plasmas'', K. Kupfer, R.W. Harvey, 
c     O. Sauter, G. Staebler, and M.J. Schaffer,  Physics of Plasmas
c     3,  3644 (1996).

      include 'param.h'
      include 'comm.h'
      double precision distn(0:iy+1,0:jx+1)

c     following diaggnde integration
      call bcast(tam3,zero,jx)
      do j=1,jx
         do i=1,iy
            tam3(j)=tam3(j)+distn(i,j)*cynt2(i,ll)*coss(i,ll)
         enddo
      enddo
      
      cn=zero
      do j=1,jx
         tam3(j)=tam3(j)*x(j)*gammi(j)*cint2(j)*dxi(j)
         cn=cn+tam3(j)*dx(j) !we only need cn here; no need to keep tam3
      enddo

c     Adjust poloidal length dlpgpsii for eqsym.ne."none" cases
      if (eqsym.eq."none") then
         symm=one
      else
         symm=two
      endif

      ampfarl= (4.d0*pi*(-charge*vnorm))/
cBH131110     +        vnorm*drpmconz(ll)*dlpsii(ll)*cn
cBH131110, BUT think vnorm* is correct, cancels vnorm in impanvc0:dfdvz
     +        (clight*dtr) *cn    ! now [statA/cm^2 /cm]
        !YuP[02-2017] removed drpmconz(ll) factor.
        !YuP[03-2017] Removed symm factor next to dlpgpsii(ll) :
        !Both dlpgpsii(ll) and dlpsii(ll) are calc-ed along same 
        !poloidal arclength, so either both should have symm* in front
        !or both should not have.
        !
        !YuP[03-2017] removed dlpsii(ll)/dlpgpsii(ll) factor here
        ! and added it in ampsoln as dlpsii(llp)/dlpgpsii(ll) 
        ! [note that dlpsii is at llp==l' index now, not ll]

      return
      end function ampfarl



      subroutine ampfalloc
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     Allocates arrays used in Ampere-Faraday routines.
c     COULD move alloc of some other A-F specific arrays here.
c..............................................................

      include 'param.h'
      include 'comm.h'
c.......................................................................


c     NOTE: fh is alternatively allocated (differently) in wpalloc.f
c           when cqlpmod.eq.enabled.
c           fh and fg are allocated here for Kupfer functions.
c           Only needed for electrons.
      allocate(fh(0:iy+1,0:jx+1,1,1:lrz),STAT=istat)
      call bcast(fh,zero,SIZE(fh))

      allocate(fg(0:iy+1,0:jx+1,1,1:lrz),STAT=istat)
      call bcast(fg,zero,SIZE(fg))

      allocate(ampfln(1:lrz),STAT=istat)
      call bcast(ampfln,zero,SIZE(ampfln))

      allocate(ampflh(1:lrz),STAT=istat)
      call bcast(ampflh,zero,SIZE(ampflh))

      allocate(ampflg(1:lrz),STAT=istat)
      call bcast(ampflg,zero,SIZE(ampflg))

      allocate(ampfa(1:lrz,1:lrz),STAT=istat)
      call bcast(ampfa,zero,SIZE(ampfa))

      allocate(ampfb(1:lrz,1:lrz),STAT=istat)
      call bcast(ampfb,zero,SIZE(ampfb))

      allocate(ampfaa(1:lrz,1:lrz),STAT=istat)
      call bcast(ampfaa,zero,SIZE(ampfaa))

      allocate(ampfc(1:lrz),STAT=istat)
      call bcast(ampfc,zero,SIZE(ampfc))

      allocate(ampf2ebar(1:lrz),STAT=istat)
      call bcast(ampf2ebar,zero,SIZE(ampf2ebar))

      return
      end subroutine ampfalloc


      subroutine ampfinit 
      ! n is not updated yet, n=0,1,...,nstop-1
      ! Called from tdchief when n+1.eq.nonampf
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     Initialize arrays used in Ampere-Faraday routines.
c..............................................................

      include 'param.h'
      include 'comm.h'
c.......................................................................
c     Electric field elecfld(0:lrza) units are the mixed Volts/cm.  
c     (1 cgs statvolt=300 volts).  elecfld is used for FP eqn solutions.
c
c     efflag="toroidal" is ensured for ampfmod, giving elecfld is toroidal.
c     On a flux surface, toroidal electric field varies strictly as 1/R,
c     due to tor symmetry. Therefore, picking a point on a flux surface
c     to evaluate the electric field gives it everywhere on the flux surface.
c
c     For ampf: elecfldn,delecfld0,delecfld0 (cgs) are used, 
c     evaluated on each flux surface at major radius R=rmag posn on a
c     flux surface, as is the main toroidal electric field variable, elecfld.
c     elecfldn(,,) are the values of elecfld used at each radius,time step,
c     and iteration, but in statV/cm, evaluated at radial bin centers (as
c     is elecfld.
c     elecfldn(0:lrz+1,0:nstop,0:nampfmax)=elecfldc,elecfld(1:lrz),elecfldb
c                                 for each n,iteration (converted to cgs).
c     delecfld0(1:lrz,1:nstop)= elecfld(1:lrz,n+1)-elecfld(1:lrz,n)
c                              i.e., total change in elecfld over time step
c                                    evaluated at bin centers.
c     delecfld0n(1:lrz,1:nstop,1:nampfmax)=change in electric field a each
c        step n=1:nstop, and each iteration, evaluated at bin centers.
c.......................................................................

      do nn=nonampf,nstop
        ! Do not set elecfldn(.,nn,.) to 0.d0 for nn<nonampf
        ! They were saved as the background el.fld at earlier steps.
        ! See tdchief, L~639
         do ll=1,lrz
            psi0bar(ll,nn)=one  !Not presently recalc'd or used.
         enddo  !On ll
      enddo  !On nn
      do nn=1,nstop
         do ll=1,lrz
            delecfld0(ll,nn)=zero
            do it=1,nampfmax
               delecfld0n(ll,nn,it)=zero !actually, already done in ainalloc
            enddo
         enddo
      enddo

c     Set initial, n=nonampf profiles
      elecfld(0)=elecfldc !Check elecfldc is set properly

cBH170329
      do ll=0,lrz
         elecfldn(ll,nonampf,0)=elecfld(ll)/300.d0
      enddo

c     Fill in elecfldn for n.lt.nonampf
c     with the n=nonampf values, for plotting purposes.
c YuP[21-08-2017] This is already done in tdchief, line ~700,
c and it is done properly: elecfldn(ll,n+1,niter)=elecfld(ll)/300.d0
c for all niter=0,nampfmax.
c So commenting it here:
c      if (nonampf.gt.0) then
c         do nn=0,nonampf-1
c            do ll=0,lrz
c               do it=0,nampfmax
c                  elecfldn(ll,nn,it)=elecfldn(ll,nonampf,1) ! for all it
c               enddo
c            write(*,*)'nn,ll,elecfldn*300=',nn,ll,elecfldn(ll,nn,1)*300
c            enddo
c         enddo
c      endif

      write(*,*)'-------------- done with ampfinit ------------'

      return
      end subroutine ampfinit




      subroutine ampfefldb(nn,timett)
      ! n is not updated yet, n=0,1,...,nstop-1
      ! Called from tdchief when n+1.eq.nonampf, as ampfefldb(n+1,time+dtr)
      ! So, nn= nonampf,nonampf+1,...,nstop
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     Get specified time-advanced boundary toroidal electric field.
c     Follow methods in subroutines tdxinitl and profiles.
c     BH: Probably useful to add one-toroidal turn voltage as function 
c     of time (or constant) to the namelist input.
c     Another input method, more in keeping with experiments, would 
c     specify total plasma toroidal current as a function of time.
c     Subroutine ampfsol would iterate the boundary voltage to find
c     efld giving the specified  current.
c     nn=time-advanced step
c     timett=time-advanced time
c.......................................................................

      include 'param.h'
      include 'comm.h'

cBH131107: For starters, only use elec field specified by fixed parabola
cBH131107: input, or prbola-t, with efswtch=method1, tmdmeth=method1.

      if (nbctime.le.0) then !Time-indep edge bndry value elecfldb
      
c     Follow subroutine tdxinitl
        if (iproelec.eq."parabola") then
          elecfldn(lrz+1,nn,0)=elecfldb/300.d0 !here: ampfefldb(nn.ge.nonampf)
          ! nn= nonampf,nonampf+1,...,nstop
        elseif(iproelec.eq."spline") then ! YuP[03-2017] added option
          ! Was set in tdxinitl:
          !call tdinterp("zero","free",ryain,elecin,njene,rya(1),
     +    !       elecfld(1),lrzmax)
          !elecfldc=elecin(1)
          !elecfldb=elecin(njene)
          elecfldn(lrz+1,nn,0)=elecfldb/300.d0 !here: ampfefldb(nn.ge.nonampf)
        else
          STOP 
     +   'ampfefldb: only setup for iproelec.eq."parabola" or "spline" '
        endif

      else  !That is, time-dependent, nbctime.gt.0
c     Follow subroutine profiles
c     [But maybe just use eledfldb from sub profiles?]
      itab(1)=1
      itab(2)=0
      itab(3)=0

      if (tmdmeth.eq."method1") then
         itme=0
         do jtm=1,nbctime
            if (timett.ge.bctime(jtm)) itme=jtm
         enddo
         itme1=itme+1
      endif  !On tmdmeth
      
      if (efswtch.eq."method1") then
         if (iproelec.eq."prbola-t") then
            
            if (tmdmeth.eq."method1".and.elecc(1).ne.zero) then
               
               if (itme.eq.0) then
                  elecfldn(lrz+1,nn,0)=elecb(1)/300.d0
               elseif (itme.lt.nbctime) then
                  elecfldn(lrz+1,nn,0)=(elecb(itme)+(elecb(itme1)
     1                 -elecb(itme))/(bctime(itme1)-bctime(itme))
     1                 *(timett-bctime(itme)))/300.d0
               else
                  elecfldn(lrz+1,nn,0)=elecb(nbctime)/300.d0
               endif
            endif
         endif
      endif  !On efswtch
      
      endif  !On nbctime

c     Scale boundary electric field
      elecfldn(lrz+1,nn,0)=elecscal*elecfldn(lrz+1,nn,0) !here: ampfefldb(nn.ge.nonampf)
      
      !YuP[21-08-2017] Added: copy the boundary value for iteration=0
      ! to all other iterations:
      do niter=0,nampfmax ! or up to nefitera
         elecfldn(lrz+1,nn,niter)=elecfldn(lrz+1,nn,0) !here: ampfefldb(nn.ge.nonampf)
      enddo
      

c     Set zero iteration of elecfldn equal to previous time
c     step radial profile, or previous iteration.
      !it_prev=0 ! YuP[21-08-2017] was 0
      it_prev=nampfmax ! But logically, should be the last? YuP[21-08-2017]
      do ll=0,lrz
         elecfldn(ll,nn,0)=elecfldn(ll,nn-1,it_prev)
         !YuP: Now: previous time step, but the LAST iteration at that step
      enddo

c      do ll=0,lrz
c      write(*,*)'ampfefldb: n,ll,elecfldn(ll,n+1,0)*300=',
c     +                      n,ll,elecfldn(ll,nn,0)*300
c      enddo

c     Zero delecfld0
      do ll=1,lrz
         delecfld0(ll,nn)=zero
         delecfld0n(ll,nn,1)=zero
      enddo

     
      return
      end subroutine ampfefldb


      subroutine ampfdiff(iflag)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

c     No real test here, to start with.
c     Future: based on (elecfldn(ll,nn,it+1)-elecfldn(ll,nn,it+1))/
c             elecfldn(ll,nn,it+1).

      iflag=1  !returning test not satisfied, so all iterations carried out.

      return
      end subroutine ampfdiff



      subroutine ampfsoln(it,nn) ! Called from tdchief 
                                 ! here nn=nonampf,...,nstop
                                 ! it=1,nampfmax
      implicit integer (i-n), real*8 (a-h,o-z)
      integer, allocatable:: ipiv(:)

c.......................................................................
c     Updates the toroidal electric field one iteration, using 
c     Ampere-Faraday eqns.  The Kupfer h and g functions have
c     been stored in fh and fg, for this iteration.
c     it=iteration number, starting at 1
c     nn=(advanced) time step, starting at 1
c.......................................................................

      include 'param.h'
      include 'comm.h'


      do ll=1,lrz
         ampfln(ll)=ampfarl(f_(0,0,kelec,ll),ll)
         ampflh(ll)=ampfarl(fh(0,0,kelec,ll),ll)
         ampflg(ll)=ampfarl(fg(0,0,kelec,ll),ll)
      enddo

      do ll=1,lrz ! == l index in notes
c        BH140302: adjusting - to +, and back
         ampftmp=-(one/(clight*rmag))*drpmconz(ll)*rpconz(ll)
         do llp=1,lrz ! == l' index in notes
            !YuP[03-2017] removed dlpsii(ll)/dlpgpsii(ll) factor in ampfarl 
            !and added it here as dlpsii(llp)/dlpgpsii(ll) 
            ![note that dlpsii is at llp==l' index now, not ll]
            dldl= dlpsii(llp)/dlpgpsii(ll) ! YuP added here.
            ! For a cylinder, this is approximately r(l')/r(l),
            ! where r(l) is the minor radius of flux surf. #l.
            ampfa(ll,llp)=dldl*ampftmp*drpmconz(llp)*ampflg(llp)
            ampfb(ll,llp)=dldl*ampftmp*drpmconz(llp)*
     +                         (ampflh(llp)-ampfln(llp))
         enddo
      enddo

c     Set up matrix ampfaa for the delecfld0n(1:lrz,nn,it)
ctmp      Set up ampfaa (coefficient next to deltaE(llp,n+1)
      call bcast(ampfaa,zero,lrz*lrz)
      do ll=1,lrz     !== l  in notes 
         do llp=1,lrz !== l' in notes
            if (llp.eq.ll) then ! l'=l
               ampfaa(ll,llp)= two -ampfa(ll,llp)
            ! Summing a_{l,l'} over l'  -  Up to l or l-1 ?  
            ! (if - up to l-1, comment ampfa(ll,llp) in the above 
            ! Tests: a very little difference
            ! (summing up to l-1 gives a little bit faster decay of I current)
            elseif (llp.gt.ll) then ! l'>l
               ! This should be only present if(ll.le.lrz-1)
               ! But from condition ll < llp (<= lrz)
               ! it is clear that here we can only find ll.le.lrz-1
               ampfaa(ll,llp)= four*(-one)**(llp-ll)
            else ! llp < ll  (l' < l in notes)
               ampfaa(ll,llp)= -ampfa(ll,llp)
            endif
         enddo
      enddo

      write(*,'(a,2i4)')'ampfsoln: nn,it-1=', nn,it-1
      do ll=0,lrz+1
        write(*,'(a,i4,2e15.6)')
     +'ampfsoln before soln: ll,elecfldn(ll,nn,it-1)*300,elecfld(ll)=',
     +                       ll,elecfldn(ll,nn,it-1)*300,elecfld(ll)
      enddo

c     Set up the rhs vector
ctmp      Setup ampf2ebar,ampfc  (1:lrz), Check elecfldn has it=0 dimensioning
ctmp      and set up.  elecfldn(,,0) is set up from elecfld at previous
ctmp      time step.
      do ll=1,lrz
         ! This is 2E^bar (l,n) in notes
cBH170329 ampf2ebar(ll)=elecfld0n(ll,nn,it-1)+elecfld0n(ll-1,nn,it-1)![statV/cm]
         ampf2ebar(ll)=two*elecfldn(ll,nn,it-1) ![statV/cm]
         !here nn is from nonampf,...,nstop  range
      enddo

      call bcast(ampfc,zero,lrz)
      do ll=1,lrz
         ! This is -2E^bar(l) -2E(lrz)
         ampfc(ll)= -ampf2ebar(ll) 
     +     +two*(-one)**(lrz-ll)*elecfldn(lrz+1,nn,it-1) ![statV/cm]
         ! YuP: the edge term is corrected
         if(ll.le.lrz-1)then
         do llp=ll+1,lrz ! (l < l' < = lrz in notes)
            ampfc(ll)=ampfc(ll)-two*(-one)**(llp-ll)*ampf2ebar(llp) ![statV/cm]
            ! YuP[02-2017] corrected (-one)**(llp-l) to (-one)**(llp-ll)
         enddo
         endif
         do llp=1,ll !ll-1 !Add the SUM[b(l,l')], with l' from 1 to l
            ! Up to ll or ll-1 ?   Tests: a very little difference
            ! (summing up to ll-1 gives a little bit faster decay of I current)
            ampfc(ll)=ampfc(ll)+ampfb(ll,llp) ![statV/cm]
         enddo
      enddo

c     Solve linear equations for delecfld0n(ll,nn,it) using
c     Gaussian elimination with pivoting.

      if (.not. allocated(ipiv)) then
         allocate(ipiv(lrz))
      endif
      call dgesv(lrz,1,ampfaa,lrz,ipiv,ampfc,lrz,info)
      !    DGESV( N,NRHS, A,  LDA,IPIV,  B,  LDB,INFO ) 
      !computes the solution to a system of linear equations A * X = B
      !So, the matrix A(N,N) is ampfaa, and the vector B is ampfc.
      !On exit, B contains the solution vector X, which corresponds to
      !the values of delecfld0n(ll,*,*) for all ll=1:lrz,
      !which is the time-step (or iteration) change in electric field 
      !evaluated at bin centers ll=1:lrz.
      
      if (info .ne. 0) then
         print *,' warning after dgesv in ampfsoln: info= ',info
         stop 'ampfsoln'
      endif

c     Update elecfld (tor elec fld) to present iteration
c     The ampfc soln is based on statvolt/cm elecfld.
c     (Alternatively, consider code mod to work entirely in V/cm?)
      do ll=1,lrz
         delecfld0n(ll,nn,it)=ampfc(ll) ! statvolt/cm
         elecfldn(ll,nn,it)=elecfldn(ll,nn,it-1)+delecfld0n(ll,nn,it)
         elecfld(ll)=elecfldn(ll,nn,it)*300.d0 ! V/cm
      enddo
      

      do ll=0,lrz+1
        write(*,'(a,i4,2e15.6)')
     +'ampfsoln  after soln: ll,elecfldn(ll,nn,it)*300,elecfld(ll)=',
     +                       ll,elecfldn(ll,nn,it)*300,elecfld(ll)
      enddo
      !here nn is from nonampf,...,nstop  range
      
c     Update total delecfld0 for this time step nn
      do ll=1,lrz
         delecfld0(ll,nn)=delecfld0(ll,nn)+delecfld0n(ll,nn,it) ![statV/cm]
      enddo

      return
      end subroutine ampfsoln
