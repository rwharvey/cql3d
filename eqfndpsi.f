c
c
      subroutine eqfndpsi(psides,areades,volum)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE     


c     Set parameter for trapped particle frac calc (150 in ONETWO).
      parameter (nlam=150)
      dimension suml(nlam)


c..................................................................
c     This subroutine is called from subroutine eqcoord.
c
c     If 0 .le. rovera(lr_) .le. 1. then:
c     This routine does a Newton's iteration to determine the
c     value of the equilibrium psi, psides, associated with radial
c     coordinate (passed in common) erhocon(lr_)=rovera(lr_)*rhomax:
c     rhomax=max[rya,1.e-8], rya() is namelist input.
c     rhomax is determined in eqrhopsi for the radial coord
c     choice specified by namelist variable radcoord.
cBH090811:  Actually, from subroutine eqrhopsi, rhomax is obtained
cBH090811:  from a linear extrapolation of eqrho(eqpsi) to the
cBH090811:  eqdsk value psilim (eqsource='eqdsk').  This gives
cBH090811:  a more accurate value of rhomax, evidently for increased
cBH090811:  accuracy.
c     
c     The arrays eqrho(j), eqpsi(j) and eqfopsi(j), j=1:nconteqn, 
c     have been calculated in eqrhopsi.
c     eqpsi,eqrho are corresponding radial psi, and coord values 
c     (in accord with radcoord), and  eqfopsi=f=R*B_phi.
c
c     If (rovera(lr_).lt.0.) then:
c     Set psides=povdelp*delp and find the contour such that
c     psi=psides directly (no iterations required).
c..................................................................
      !YuP[2020-06-30] Added a new option for method of searching 
      !for the flux surface that corresponds 
      !to the target value of "rhodes" (rho value from rya() list).
      ! The old method is accessed with kopt=1.
      ! It is based on narrowing the window [psi1;psi2]
      ! (and concurrently,  narrowing of [rho1;rho2] window)
      ! that converges to psinew (while [rho1;rho2] converges to rhonew),
      ! maybe from one side only, or from both sides of the window.
      ! The iterations are done for  
      !    psinew= psi1 + (psi2-psi1)*
      !         ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
      ! where psi1 and rho1 (or psi2 and rho2) are adjusted 
      ! at each iteration.
      ! The problem with this method is that in the region 
      ! close to magnetic axis the profile of psi (poloidal flux)
      ! is nearly flat, and so the values of psi1 and psi2 could be
      ! very close to each other, while rho1 and rho2 are distinct 
      ! (especially when radcoord is related to enclosed volume, etc).
      ! The new option kopt=2 below uses a different method - 
      ! instead of "psi", the values of major radius R 
      ! for the start of a flux surface are iterated.
      ! The iterative formula is 
      !    rnew= rp1 + (rp2-rp1)*
      !         ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
      ! where [rp1;rp2] is the window in R that contains the starting
      ! point (R,Z)==(rnew,zmag) for the flux surface.
      ! At each iteration the surface is traced from (rnew,zmag),
      ! and the value of rhonew is calculated along such surface, 
      ! and then compared to the target value of rhodes.
      ! Then, the window [rho1,rho2], together with [rp1,rp2],
      ! is adjusted to make it narrower,
      ! and rnew is recalculated from the above formula.
      ! This new method shows a better convergence near magnetic axis,
      ! although not always it can resolve the problems.
      ! It may happen that the window [rp1;rp2] is narrow to almost a point,
      ! and yet, from one iteration to another the value of rhonew 
      ! changes in 1st digit (hence, "POOR CONVERG." messages).
      ! Such sensitivity of rho value to a tiny variation
      ! of starting point in R could be caused by a poor quality 
      ! of eqdsk data, particularly near the magnetic axis.
      ! It is recommended to recompute eqdsk file with a refined (R,Z) grid,
      ! or simply avoid the region close to magnetic axis in cqlinput setup.
c..................................................................

CMPIINSERT_IF_RANK_EQ_0      
c      write(*,*)'eqfndpsi: radcoord,rhomax,lr_,erhocon(lr_)= ',
c     +                     radcoord,rhomax,lr_,erhocon(lr_)
CMPIINSERT_ENDIF_RANK


      if (rovera(lr_).ge.0.) then
c..................................................................
c     Begin by finding the first index jval such that eqrho(jval) is
c     larger than rhodes.
c..................................................................

        rhodes=erhocon(lr_)
        if (rhodes.gt.rhomax) call eqwrng(8)
        do 10 j=2,nconteqn
          if (rhodes.le.eqrho(j)) go to 11
 10     continue
 11     continue
        jval=j
        
CMPIINSERT_IF_RANK_EQ_0      
        if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
        WRITE(*,*)
        WRITE(*,*)'eqfndpsi: lr_,rhodes,jval,eqrho(jval-1),eqrho(jval)',
     +                       lr_,rhodes,jval,eqrho(jval-1),eqrho(jval)
        WRITE(*,*)'eqfndpsi: eqrpcon(jval-1),eqrpcon(jval)',
     +                       eqrpcon(jval-1),eqrpcon(jval)
        WRITE(*,*)'eqfndpsi: eqrmcon(jval-1),eqrmcon(jval)',
     +                       eqrmcon(jval-1),eqrmcon(jval)
c        write(*,*)'eqfndpsi: eqrho(j),j=1,nconteqn',
c     +                     (eqrho(j),j=1,nconteqn)
        endif
CMPIINSERT_ENDIF_RANK

c..................................................................
c     Begin iteration loop
c..................................................................
        kopt=1 !2 ! Need to Add it into INPUT list !!!!!!!!!!!!!!!!

        !---> For kopt=1 (original design):
        psi2=eqpsi(jval) ! Initial guess
        rho2=eqrho(jval) ! Initial guess
        psi1=eqpsi(jval-1) ! Initial guess
	  rho1=eqrho(jval-1) ! Initial guess: 
	                     ! the target surface is between rho1 and rho2
        !---> For kopt=2 ( YuP[2020-06-30] design):
        rp2=eqrpcon(jval)   ! Initial guess (major radius, lower limit)
        rp1=eqrpcon(jval-1) ! Initial guess (major radius, upper limit)
	  ! The target surface (with rhodes) is between rp1 and rp2

	  !psimag=0 in a mirror machine.
        if(jval.le.2) then !-YuP: for better convergence near m.axis 
           psi1=psimag
           rho1=0.d0
        endif
        
        iter=0
 20     continue !-> iteration loop (through Line~175) ------------------
 
        if(kopt.eq.1)then
          !Bi-linear interpolation -> iterate psinew:
          psinew= psi1 + (psi2-psi1)*
     &          ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
          epsicon(lr_)=psinew !Saved; Also an INPUT for eqorbit, when kopt=1
          eqcall="disabled"
          rstart=0.d0 ! OUTPUT, when kopt=1
          zstart=zmag ! OUTPUT, when kopt=1
          !--------------------
          call eqorbit(1,psinew,rstart,zstart) ! Get (solr_(l),solz_(l)) tracing flux surface
          !--------------------
          rnew=rstart ! OUTPUT, when kopt=1
          !zstart is  OUTPUT, when kopt=1
        endif ! kopt=1
        
        if(kopt.eq.2)then
          !YuP[2020-06-30] Added kopt=2, which allows tracing surface
          ! directly from point (rstart,zstart) when it is given in INPUT
          ! (in this case value of epsicon_ is not needed).
          ! For the original design, use kopt=1,
          ! which means: find the starting point from knowledge of epsicon_ 
          ! and zmag coordinate (stored in comm.h).
          ! Iterate for rnew (rp1 and rp2 are adjusted after each iteration):
          rnew= rp1 + (rp2-rp1)*
     &          ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
          eqcall="disabled"
          rstart=rnew ! INPUT, when kopt=2
          zstart=zmag ! INPUT, when kopt=2
          !For kopt=2, psinew is not an input, but it is found in eqorbit().
          !--------------------
          call eqorbit(2,psinew,rstart,zstart) ! Get (solr_(l),solz_(l)) tracing flux surface
          !--------------------
          epsicon(lr_)=psinew ! Found in eqorbit(), when kopt=2.
        endif ! kopt=2
        
        do 70 l1=1,lorbit_
          ! The r.h.s. are values from eqorbit:
          solr(l1,lr_)=solr_(l1)
          solz(l1,lr_)=solz_(l1)
          es(l1,lr_)=es_(l1)
          eqbpol(l1,lr_)=eqbpol_(l1) !=Bpol along field line (pol.plane)
          bpsi(l1,lr_)=bpsi_(l1)
          thtpol(l1,lr_)=thtpol_(l1)
          eqdell(l1,lr_)=eqdell_(l1)
CMPIINSERT_IF_RANK_EQ_0      
          if(eqbpol_(l1).eq.0. )then
             WRITE(*,'(a,i5,2e17.10)')
     +       'eqfndpsi: l1,eqdell_(l1),eqbpol_(l1) ',
     +                  l1,eqdell_(l1),eqbpol_(l1)
          endif
CMPIINSERT_ENDIF_RANK
 70     continue
        ! The r.h.s. (*_) are values from eqorbit:
        eqdells(lr_)=eqdells_ ! =dl along field line (pol.plane)
        lorbit(lr_)=lorbit_
        rmcon(lr_)=rmcon_
        rpcon(lr_)=rpcon_
        zmcon(lr_)=zmcon_
        zpcon(lr_)=zpcon_
        es_bmax(lr_)=es_bmax_
        bpsi_max(lr_)=bpsi_max_
        bpsi_min(lr_)=bpsi_min_
        lbpsi_max(lr_)=lbpsi_max_
        lbpsi_min(lr_)=lbpsi_min_
        bthr(lr_)=bthr_
        btoru(lr_)=btoru_
        fpsi(lr_)=fpsi_
        zmax(lr_)=zmax_
        btor0(lr_)=btor0_
        bthr0(lr_)=bthr0_
        bmidplne(lr_)=bmidplne_  !At min bpsi_ point, not necessarily
                                 !the midplane, for eqsym.eq."none"
        !Get average values over surface epsicon(lr_) (corr.to rhodes)                      
        eqorb="disabled"
        call eqvolpsi(epsicon(lr_),volum,areac)
        call eqonovrp(epsicon(lr_),onovrp1,onovrp2) !get <1/R>, <1/R**2>
        
        ! Find flux-surf. averaged values of tlorb1(l),
        ! save into bpolsqa_,flxavgd_
        ! Note: epsicon_ is INPUT, but only needed when eqorb="enabled"
        ! which is not the case here.
        do 60 l=1,lorbit_
          tlorb1(l)=eqbpol_(l)**2 !=Bpol^2 along given surface (iterated)
 60     continue
        call eqflxavg(epsicon_,tlorb1,bpolsqa_,flxavgd_)
        !OUT= bpolsqa_= SUM(tlorb1*dl/Bpol)/SUM(dl/Bpol)
        
        ! Similarly, Find flux-surf. averaged values of tlorb1(l),
        ! save into psiovr_,flxavgd_
        do 40 l=1,lorbit_
          tlorb1(l)=bpsi_(l)/solr_(l) !=(B(l)/B0)/R
 40     continue
        call eqflxavg(epsicon_,tlorb1,psiovr_,flxavgd_)
        !OUT= psiovr_= SUM(tlorb1*dl/Bpol)/SUM(dl/Bpol)
        !(and again, epsicon_ is not used here)
        
        bpolsqa(lr_)=bpolsqa_
        psiovr(lr_)=psiovr_
        flxavgd(lr_)=flxavgd_ !=SUM(dl/Bpol)
        onovrp(1,lr_)=onovrp1 ! <1/R> over lr_ surface
        onovrp(2,lr_)=onovrp2 ! <1/R^2> over lr_ surface
        
        ! Get toroidal "f(psi)" (=R*B)
        ! for psi=epsicon(lr_) surface; based on interpolation:
        call eqfpsi(epsicon(lr_),fpsi_,fppsi_)

        ! Now the value of rhonew (iterated value) can be evaluated:
        if (radcoord.eq."sqtorflx") then
        !YuP[2019] Need to inspect this part in detail.
        !YuP: Try to use Bi-linear interpolation, 
        !similar to initial guess for psinew:
        ! psinew= psi1 + (psi2-psi1)*
        !         ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
        ! where psi1=eqpsi(jval-1), psi2=eqpsi(jval)
        ! and similarly for rho1,rho2
           dvolum=(volum-eqvol(jval-1))
           fpsih=(fpsi_+eqfopsi(jval-1))*.5 !similar to eqfh below
           onok=.5*(onovrp(1,lr_)+eqovrp(jval-1,1)) ! similar to eqrph
           tem=eqrho(jval-1)**2*pi*btor ! value at previous radial point
           onoh=(onovrp(2,lr_)+eqovrp(jval-1,2))*.5
           rhonew=tem+onoh*dvolum*fpsih/pi*0.5
           !Definition from eqrhopsi:
           ! dvolum=(eqvol(j)-eqvol(j-1))
           ! eqrph=(eqovrp(j,2)+eqovrp(j-1,2))*.5
           ! eqfh=(eqfopsi(j)+eqfopsi(j-1))*.5
           ! eqrho(j)=eqrho(j-1)+eqrph*dvolum*eqfh/pi*0.5
           !And then, take sqrt, as here below:
           areanew=areac
           rhonew=sqrt(rhonew/pi/btor)
        elseif (radcoord.eq."sqarea") then
           rhonew=sqrt(areac/pi)
           areanew=areac
        elseif (radcoord.eq."sqvol") then
           rhonew=sqrt(volum/(2.*pi**2*rmag))
           areanew=areac
        elseif (radcoord.eq."rminmax") then
           rhonew=0.5*(rpcon_-rmcon_)
           areanew=areac
        elseif (radcoord.eq."polflx") then
           rhonew=(psinew-psimag)/(psilim-psimag)
           areanew=areac
        elseif (radcoord.eq."sqpolflx") then
           rhonew=sqrt((psinew-psimag)/(psilim-psimag))
           areanew=areac
        endif
         
        !Compare the given rhonew (iterated) with target value of rhodes:    
        err=abs(rhonew-rhodes)/rhodes ! rhodes is the target.
        iter=iter+1 ! count iterations; usually 2-4 is sufficient

CMPIINSERT_IF_RANK_EQ_0      
      if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
!      WRITE(*,'(a,2i5,2f16.10,e14.6,e10.3)')
!     +'eqfndpsi: iter,lorbit_,rhonew,rhodes,psinew-psimag,err',
!     +        iter,lorbit_,rhonew,rhodes,psinew-psimag,err
      WRITE(*,'(a,i4,3f14.7)')'iter,rho1,rhonew,rho2=',
     +                         iter,rho1,rhonew,rho2
      WRITE(*,'(a,i4,3f14.7)')'iter,rp1, rnew,  rp2= ',
     +                         iter,rp1, rnew,  rp2
      endif
CMPIINSERT_ENDIF_RANK

        if (err.gt.1.e-5 .and. iter.lt.35) then !max number of iter: was 25
          if (rhonew.gt.rhodes) then
            rho2=rhonew 
            psi2=psinew ! For kopt=1 only
            rp2= rnew   ! For kopt=2 only
          else
            rho1=rhonew
            psi1=psinew ! For kopt=1 only
            rp1= rnew   ! For kopt=2 only
            !if(iter.gt.25)then
              !Poor convergence (usually near rho=0):
              ! try to reset rho2 and psi2 :
            !  rho2=0.d0
            !  psi2=psimag
            !endif
          endif
          go to 20  !  GO TO NEXT ITERATION
        endif
        
        if (err.gt.1.e-2) then

CMPIINSERT_IF_RANK_EQ_0      
           WRITE(*,'(a,i4,e12.3,i4,e12.3)')
     +      'eqfndpsi/WARNING: POOR CONVERG. lr,rhonew,iter,err=',
     +                                       lr_,rhonew,iter,err
CMPIINSERT_ENDIF_RANK

           !!!YuP call eqwrng(4) !-> will stop the job
           !Sometimes the convergence is poor near magn.axis,
           !presumably because PSI(rho) is nearly flat.
        endif
        psides=psinew
        areades=areanew

        fppsi(lr_)=fppsi_
        call eqppsi(epsicon(lr_),ppsi_,pppsi_)
        pppsi(lr_)=pppsi_
        

c.......................................................................
c     Now for the case that rovera(lr_)=.lt.0.
c.......................................................................

      else  !rovera(lr_).lt.0
        delp=(psimag-psilim)
        psides=psimag-povdelp*delp
        do 30 j=2,nconteqn
          if (eqpsi(j).lt.psides) go to 31
 30     continue
 31     continue
        jval=j
        epsicon(lr_)=psides
        rstart=0.d0 ! OUTPUT, when kopt=1
        zstart=zmag ! OUTPUT, when kopt=1
        call eqorbit(1,psides,rstart,zstart) ! Get (solr_(l),solz_(l)) tracing flux surface
        !YuP[2020-06-30] Added kopt=2, which allows tracing surface
        ! directly from point (rstart,zstart) when it is given in INPUT
        ! (in this case value of epsicon_ is not needed).
        ! For the original design, use kopt=1,
        ! which means: find the starting point from knowledge of epsicon_ 
        ! and zmag coordinate (stored in comm.h).
        do 90 l1=1,lorbit_
          solr(l1,lr_)=solr_(l1)
          solz(l1,lr_)=solz_(l1)
          es(l1,lr_)=es_(l1)
          eqbpol(l1,lr_)=eqbpol_(l1)
          bpsi(l1,lr_)=bpsi_(l1)
          thtpol(l1,lr_)=thtpol_(l1)
          eqdell(l1,lr_)=eqdell_(l1)
          bpsi(l1,lr_)=bpsi_(l1)
 90     continue
        eqdells(lr_)=eqdells_
        lorbit(lr_)=lorbit_
        rmcon(lr_)=rmcon_
        rpcon(lr_)=rpcon_
        bthr(lr_)=bthr_
        btoru(lr_)=btoru_
        fpsi(lr_)=fpsi_
        zmax(lr_)=zmax_
        btor0(lr_)=btor0_
        bthr0(lr_)=bthr0_
        eqorb="disabled"
        call eqvolpsi(epsicon(lr_),volum,areac)
        call eqonovrp(epsicon(lr_),onovrp1,onovrp2)
        do 50 l=1,lorbit_
          tlorb1(l)=bpsi_(l)/solr_(l)
 50     continue
        call eqflxavg(epsicon_,tlorb1,psiovr_,flxavgd_)
        do 80 l=1,lorbit_
          tlorb1(l)=eqbpol_(l)**2
 80     continue
        call eqflxavg(epsicon_,tlorb1,bpolsqa_,flxavgd_)
        psiovr(lr_)=psiovr_
        flxavgd(lr_)=flxavgd_
        onovrp(1,lr_)=onovrp1
        onovrp(2,lr_)=onovrp2
        call eqfpsi(epsicon(lr_),fpsi_,fppsi_)

        if (radcoord.eq."sqtorflx") then
           fpsih=(fpsi_+eqfopsi(jval-1))*.5
           onok=.5*(onovrp(1,lr_)+eqovrp(jval-1,1))
           tem=eqrho(jval-1)**2*pi*btor
           onoh=(onovrp(2,lr_)+eqovrp(jval-1,2))*.5
           rhonew=tem+onoh*(volum-eqvol(jval-1))*fpsih/pi*0.5
           areanew=areac
           rhonew=sqrt(rhonew/pi/btor)
        elseif (radcoord.eq."sqarea") then
           rhonew=sqrt(areac/pi)
           areanew=areac
        elseif (radcoord.eq."sqvol") then
           rhonew=sqrt(volum/(2.*pi**2*rmag))
           areanew=areac
        elseif (radcoord.eq."rminmax") then
           rhonew=0.5*(rpcon_-rmcon_)
           areanew=areac
        elseif (radcoord.eq."polflx") then
           rhonew=psides
           areanew=areac
        elseif (radcoord.eq."sqpolflx") then
           rhonew=sqrt(psides)
           areanew=areac
        endif
              
        areades=areanew
        psinew=psides
        erhocon(lr_)=rhonew
        rovera(lr_)=erhocon(lr_)/rhomax

        fppsi(lr_)=fppsi_
        call eqppsi(epsicon(lr_),ppsi_,pppsi_)
        pppsi(lr_)=pppsi_

      endif  !on rovera(lr_).ge./.le. 0.
c
c.......................................................................
c     compute <bpsi> and <bpsi**2>,<1/(bpsi*R**3>
c.......................................................................
c
      do 41 l=1,lorbit_
         tlorb1(l)= bpsi_(l)
         tlorb2(l)= bpsi_(l)**2
 41   continue
      call eqflxavg(epsicon_,tlorb1,zpsiavg,flxavgd_)
      call eqflxavg(epsicon_,tlorb2,zpsi2av,flxavgd_)
      psiavg(1,lr_)=zpsiavg
      psiavg(2,lr_)=zpsi2av
      do l=1,lorbit_
         tlorb1(l)=1./(bpsi_(l)*solr_(l)**3)
      enddo
      call eqflxavg(epsicon_,tlorb1,zpsiavg,flxavgd_)
      onovpsir3(lr_)=zpsiavg

c.......................................................................
c     Calculate effective trapped particle fraction 
c     (see e.g., Hirshman and Sigmar, Nucl. Fus. 21, 1079 (1981), 
c      Eq. 4.54)
c     trapfrac=1.-0.75*<B**2>*integral[0,Bmax**-1]{lambda*dlambda/
c                                        <(1-lambda*B)**0.5>}
c     Manipulated (BH) this to following (with zeta**0.5=Bmax*lambda)
c     to be close to a ONETWO expression, and integrated as below:
c     trapfrac=1-0.75*<h**2>*integral[0,1]{d_zeta/
c                                         <(2.*sqrt(1-zeta**.5*h)>
c      where <...> is flux surface avg, h=B/Bmax.
c.......................................................................

      !bmaxbmin=bpsi(lorbit_,lr_) ! YuP: maybe bpsi(lbpsi_max(lr_),lr_) ?
      ! For a general case of eqsym:
      bmaxbmin=bpsi_max(lr_) ! YuP [July 2014] ==Bmax/Bmin
      
      h2fsa=psiavg(2,lr_)/bmaxbmin**2
      
      dlam=1./(nlam-1.)
      do ilam=1,nlam
         rtlam=sqrt((ilam-1)*dlam)
         do l=1,lorbit_
            hlam=bpsi_(l)/bmaxbmin ! = B(l)/Bmax
            val=abs(1.0-rtlam*hlam)
            tlorb1(l)=sqrt(val)
         enddo
         call eqflxavg(epsicon_,tlorb1,suml(ilam),flxavgd_)
      enddo
         
      xi0=0.
      do ilam=1,nlam-1
         xi0=xi0+0.25*(1.0/suml(ilam) + 1.0/suml(ilam+1))*dlam
      enddo

      trapfrac(lr_)=1.-0.75*h2fsa*xi0
      
CMPIINSERT_IF_RANK_EQ_0      
      if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
      WRITE(*,*)'eqfndpsi/END: lr_,iter,rhonew,rhodes,solr_(1)=', 
     +  lr_,iter,rhonew,rhodes, solr_(1)
      endif
CMPIINSERT_ENDIF_RANK

      return
      end
