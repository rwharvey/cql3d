c
c
      subroutine losscone
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This routine computes the location of the loss orbits and puts
c     the information into gone(i,j,k,indxlr_).
c     [Subroutine called from ainitial. indxlr_ is set by tdnflxs(ll)
c      in a loop ll=lrors,1,-1 in tdinitl.] 
c     Several loss models are available and are invoked by setting
c     lossmode(k) to the chosen "string" in the input:
c     "snk0"= orbits whose energy is less than esink
c     are lost from the system with characteristic time equal to
c     the bounce time.
c     "snk0accl"=same as "snk0" except the loss mechanism is acceler-
c     ated by factor xsinkm (input).
c     "electron" forces a loss cone for general species electrons
c     for perpendicular energies greater that eperc and parallel
c     energies greater than eparc.
c     "mirrorcc"=potential square well model for tandem mirror central
c     cell calculations; orbits passing through a potential jump
c     at the throat of magnitude ephicc(kev) are lost.
c     "mirrsnk"=combination of mirrorcc and snk0.
c     "simplban"=
c     The array gone is utilized in subroutine coefload, which
c     determines suitable Krook operator loss times, depending on the
c     local value of gone.
c..................................................................

      include 'param.h'
      include 'comm.h'

c     Pointers for dynamic memory allocation, local usage:
      real*8, dimension(:), pointer :: upar,uprp,rho_a
      integer, dimension(:,:,:), pointer :: notlost
      real*8, dimension(:,:,:), pointer :: dnotlost

      do 100 k=1,ngen  !Loops down to end of the subroutine

         !call bcast(gone(0,0,k,indxlr_),zero,iyjx2) ! YuP commented-out
         !YuP[09-10-2014] do not initialize gone array here.
         ! it was already done in ainalloc.
         xmul=1.
         
         !===========================================================
         if (lossmode(k).eq."snk0" .or. lossmode(k).eq."snk0accl" .or.
     1        lossmode(k) .eq. "mirrsnk") then
            if (lossmode(k) .eq. "snk0accl") xmul=xsinkm
            if (esink.gt.em90) then
               xsink=sqrt(esink/fions(k))
c            else
c               xsink as specified in namelist
            endif
            
c..................................................................
c     Find the first mesh point x "outside" the sink hole.
c..................................................................
            
            do 10 j=2,jx
               if (x(j) .ge. xsink) go to 20
 10         continue
            j=jx
 20         continue
            do 30 jj=1,j
               do 31 i=1,iy
                  gone(i,jj,k,indxlr_)=-xmul
 31            continue
 30         continue
            
            
c..................................................................
c     Electron losses defined by eparc and eperc
c..................................................................
            
         elseif (lossmode(k).eq."electron" .and. k.eq.kelecg) then
         
            v2=vnorm*vnorm/(clight*clight)
            do 40  i=1,iy
               do 41  j=1,jx
                  gampar=sqrt(1.+(coss(i,l_)*x(j))**2*v2)
                  gamper=sqrt(1.+(sinn(i,l_)*x(j))**2*v2)
                  epar=(gampar-1.)*restmkev !YuP[01-2011] was 512.
                  eper=(gamper-1.)*restmkev !YuP[01-2011] was 512.
                  if(epar.gt.eparc(k,lr_).or.eper.gt.eperc(k,lr_)) then
                     gone(i,j,k,indxlr_)=-1.
                  else
                     gone(i,j,k,indxlr_)=0.
                  endif
 41            continue
 40         continue
            
c.......................................................................
c     Simple Banana + Gyro-Orbit loss model: 
c     -trapped particle bananas plus gyro-orbit width (circ plasma) 
c     -greater than dist to plasma edge are lost.
c     -Circulating particles with gyro-orbit width greater than
c      distance to the edge of the plasma are lost.
c     Banana width uses poloidal magnetic field which is an average
c     of Bpol at each radius and of the edge Bpol.
c     bthr0 is at outer midplane.  Might be that bth ---at thetapol
c       =pi/2--- would be more appropriate.  Might want to check this
c       by orbit integrations.  (BobH, 011125).
c.......................................................................
            
         elseif (     lossmode(k).eq."simplban" 
     +           .or. lossmode(k).eq."simplbn1"
     +	         .or. lossmode(k).eq."simplbn2") then
c$$$            if (ionce.eq.0) then
c$$$                 ionce=1
c$$$                 write(*,*)'losscone:bthr(ll),ll=1,lrzmax',
c$$$     +                               (bthr(ll),ll=1,lrzmax)
c$$$                 write(*,*)'losscone:bthr0(ll),ll=1,lrzmax',
c$$$     +                               (bthr0(ll),ll=1,lrzmax)
c$$$                 write(*,*)'Losscone:ll,eqbpol(1:lorbit,ll=1,lrzmax)'
c$$$                 do ll=1,lrzmax
c$$$                    write(*,*) ll,(eqbpol(lll,ll),lll=1,lorbit(lr_))
c$$$                 enddo
c$$$            endif

c         Project to (Last Closed Flux Surface)*simpbfac from flux 
c         surface radii tabulated on the rya() mesh  
            rpconz_max=rpconz(lrzmax)+(simpbfac-rya(lrzmax))*
     +                 (rpconz(lrzmax)-rpconz(lrzmax-1))/
     +                 (rya(lrzmax)-rya(lrzmax-1))
            bthr0edge=abs(bthr0(lrzmax)+(simpbfac-rya(lrzmax))*
     +                 (bthr0(lrzmax)-bthr0(lrzmax-1))/
     +                 (rya(lrzmax)-rya(lrzmax-1)))
            bmod0edge=abs(bmod0(lrzmax)+(simpbfac-rya(lrzmax))*
     +                 (bmod0(lrzmax)-bmod0(lrzmax-1))/
     +                 (rya(lrzmax)-rya(lrzmax-1)))
            bthr0avg=0.5*(abs(bthr0(lr_))+bthr0edge)
            x_loss=abs(bnumb(k)*charge)
     +           *(rpconz_max-rpconz(lr_))/(fmass(k)*vnorm*clight)
            r_loss=rpconz_max-rpconz(lr_)
            omcnst=abs(bnumb(k)*charge)/(fmass(k)*vnorm*clight)
c$$$            write(*,*)'Losscone: lr_,rpconz_max,x_loss',
c$$$     +                           lr_,rpconz_max,x_loss
c$$$            write(*,*)'losscone: bthr0avg,bmod0edge',
c$$$     +                           bthr0avg,bmod0edge

cBH160403: Rethinking gyro-orbit loss.  Should evaluate bmag0 at outside
cBH160403: of banana orbit.  Banana is centered on lr_ flux surface.
cBH160403: Linearly interpolate for bmod0 at outside of the the banana.
           
            drpconz0=(rpconz(2)-rpconz(1))/(rya(2)-rya(1))
            drpconz1=(rpconz(lrzmax)-rpconz(lrzmax-1))/
     +               (rya(lrzmax)-rya(lrzmax-1))
            dbmod00=(bmod0(2)-bmod0(1))/(rya(2)-rya(1))
            dbmod01=(bmod0(lrzmax)-bmod0(lrzmax-1))/
     +               (rya(lrzmax)-rya(lrzmax-1))
            do 50  i=1,iy
               do 51  j=3,jx    !was j=1,jx ! YuP
                  ! YuP[7-22-2014] skip j=1:2 (do not apply a loss cone)
                  !Sometimes a loss cone at v~0 gives unstable solution
                  !and issues for boundary condition at v=0.
                  !Other than that, the effect from skipping j=1:2
                  !is very small - difference in ~6th digit.
                  if (i.gt.itl .and. i.lt.itu) then  !ZOW trapped particles
cBH160529
cBH160529
cBH160529
cBH160529
cBH160529
cBH160529
                  ! BH: We take the abs() for delb1, since the orbit
                  !     shift is from BA radial flux surface, and outside
                  !     scrapeoff will occur on first or second part of
                  !     the bounce.
c YuP[01-2017] Commented with ! by YuP: does not work properly.
c YuP  One of problems: Failure in ngen=2 multiURF tests
c YuP  because subr. deltar is setup/called for one general species only,
c YuP  and so all deltar* arrays are saved for k=1 gen.species only.
!                     delb1=x(j)*abs(coss(i,indxlr_))/(omcnst*bthr0(lr_))
!                     call lookup(rpconz(lr_)+delb1,rpconz,lrzmax,
!     +                           wtu,wtl,linterp)
!                     rpconz1=wtl*rpconz(linterp-1)+wtu*rpconz(linterp)
!                     bmod01=wtl*bmod0(linterp-1)+wtu*bmod0(linterp)
!                     gyrorad1=x(j)*abs(sinn(i,indxlr_))/(omcnst*bmod01)
!                     gone1=0.
!                     if (delb1+gyrorad1 .gt. r_loss) gone1=-1.

!                     if (lossmode(k).eq.'simplbn1')then !YuP[2017-11-21]corr.lossmode(1)
!                        if (delb1+gyrorad1 .gt. r_loss) then
!                           gone(i,j,k,indxlr_)=-1.
!                        else
!                           gone(i,j,k,indxlr_)=0.
!                        endif   ! on delb1
!                     endif  ! on lossmode(k).eq.'simplbn1'

!                     gyrorad2=x(j)*abs(sinn(i,indxlr_))/bmod0edge/omcnst
!                     delb2=x(j)*abs(coss(i,indxlr_))/bthr0avg/omcnst
!                     gone2=0.
!                     if ((gyrorad2 + delb2).gt. x_loss/omcnst) gone2=-1.

!                     if (lossmode(k).eq.'simplbn2') then !YuP[2017-11-21]corr.lossmode(1)
!                        if ((gyrorad2 + delb2)
!     +                       .gt. x_loss/omcnst) then
!                           gone(i,j,k,indxlr_)=-1.
!                        else
!                           gone(i,j,k,indxlr_)=0.
!                        endif
!                     endif  ! on lossmode(k).eq.'simplbn2'

cBH160529 There are problems with application of NB for these models.
cBH160529 Need to first have NB birth points shifted inwards to their
cBH160529 BA radius, if they are going to be scraped off at a banana width.
cBH160529 First, compare delb with deltaRZ calc, which we will subsequently
cBH160529 implement for scrapeoff.
cBH160529 Calc orbit shift with deltarho, to be consistent with calc
cBH160529 of NB birth points (with ndeltarho='freya') 
cBH160529 Note: cql3d-FOW convention for deltarho shift.
cBH170708 Note: deltarho calc in subroutine deltar enabled for ngen=1 at this time.
cBH170708       See email YuP: 2017-01-04 12:27
cYuP: rho1 is not used, commenting out 
cYuP: (can give a problem in case of lossmode(1)=disabled; lossmode(2)=simplban)
cYuP                     rho1=rya(lr_)+x(j)*vnorm*abs(deltarho(i,1,lr_))
                  ! BH: We take the abs() for deltarho, since the orbit
                  !     shift is from BA radial flux surface, and outside
                  !     scrapeoff will occur on first or second part of
                  !     the bounce. deltarho is anti-symmetric about pi/2.
cYuP                     if (rho1.lt.zero) rho1=zero
cYuP                     if (rho1.gt.1.1d0) rho1=1.1d0
cYuP                     call lookup(rho1,rya(1),lrzmax,wtu,wtl,linterp)
cYuP                    !Interpolate, or extrapolate a bit.
cYuP                     if (linterp.eq. lrzmax .and. wtu.eq.1d0) then
cYuP                        R1=rpconz(linterp)+(rho1-rya(linterp))*drpconz1
cYuP                        b1=bmod0(linterp)+(rho1-rya(linterp))*dbmod01
cYuP                     elseif(linterp.eq.2 .and. wtl.eq.1d0) then
cYuP                        R1=rpconz(linterp)+(rho1-rya(linterp))*drpconz1
cYuP                        b1=bmod0(1)+(rho1-rya(1))*dbmod00
cYuP                        if (R1.lt.rmag) then
cYuP                           R1=rmag
cYuP                           b1=bmag
cYuP                        endif
cYuP                     else
cYuP                        R1=wtl*rpconz(linterp-1)+wtu*rpconz(linterp)
cYuP                        b1=wtl*bmod0(linterp-1)+wtu*bmod0(linterp)
cYuP                     endif
cYuP                     gyrorad=x(j)*abs(sinn(i,indxlr_))/(omcnst*bmod01)
cYuP                     delb_deltarho=R1-rpconz(lr_)
cYuP                     gone0=0.
cYuP                     if (delb_deltarho+gyrorad.gt.r_loss) gone0=-1.
cBH170708  For backward compatibility, reset  gyrorad and delb_deltarho
                    gyrorad=x(j)*abs(sinn(i,indxlr_))/(omcnst*bmod0edge)
                     delb_deltarho=x(j)*abs(coss(i,indxlr_))/(omcnst*
     + 		                   bthr0avg)
                     if (lossmode(k).eq.'simplban') then
                        if (delb_deltarho+gyrorad .gt. r_loss) then
                           gone(i,j,k,indxlr_)=-1.
                        else
                           gone(i,j,k,indxlr_)=0.
                        endif
                     endif
cBH160529
cBH160529
cBH160529      if (i.eq.(itl+1) .and. j.eq.3) then
cBH160529      write(*,*)'rpconz_max,r_loss',rpconz_max,r_loss
cBH160529      write(*,*)'i,j,lr_,rya,rho1,R,R1,delb_deltarho,'//
cBH160529     +  'delb1,delb2,gyrorad,gyrorad1,gyrorad2,gone0,gone1,gone2'
cBH160529      endif
cBH160529      write(*,*) i,j,lr_,rya(lr_),rho1,rpconz(lr_),R1,delb_deltarho,
cBH160529     +     delb1,delb2,gyrorad,gyrorad1,gyrorad2,gone0,gone1,gone2
cBH160529
cBH160529  Conclusion from comparison of above results, particularly gone0
cBH160529  (simplban, new), gone1(simplbn1, 1st adj), gone2(simplbn2, old):
cBH160529  All about the same, although simplbn2 least severe, simplbn1 is
cBH160529  most, and new method is intermediate. See log.1.simplban_extended.
cBH160529  Application of ndeltarho='freya' (a relevant correction) will 
cBH160529  reduce the co-inj losses.
cBH160529

                  else  ! On i
                  
cBH-YuP140807: Should have used bmod0() for transitting particles.
cBH-YuP140807              if (x(j)*abs(sinn(i,indxlr_))/bmod0edge
cBH160529 Could do more for passing particles with deltarho....
                     if (x(j)*abs(sinn(i,indxlr_))/bmod0(indxlr_)
     +                  .gt. x_loss) then
                        gone(i,j,k,indxlr_)=-1.
                     else
                        gone(i,j,k,indxlr_)=0.
                     endif
                  endif  ! on i
 51            continue
 50         continue
c..................................................................
c     Now for mirror cases.
c..................................................................
            
         else if (lossmode(k).eq."mirrorcc") then
         
            call lossorbm(ephicc,k)
            
c..................................................................
c     Now for file specification of prompt losses
c..................................................................
            
         else if (lossmode(k).eq."frm_file") then

c     Allocate local pointered memory for reading lossfile(k)
            allocate (uprp(n_uprp),STAT = istat)
            if(istat .ne. 0)
     .           call allocate_error("uprp, sub rdc_multi",0,istat)
            call bcast(uprp,zero,SIZE(uprp))
            
            allocate (upar(n_upar),STAT = istat)
            if(istat .ne. 0)
     .           call allocate_error("upar, sub rdc_multi",0,istat)
            call bcast(upar,zero,SIZE(upar))
            
            allocate (rho_a(n_psi),STAT = istat)
            if(istat .ne. 0)
     .           call allocate_error("rho_a, sub rdc_multi",0,istat)
            call bcast(rho_a,zero,SIZE(rho_a))
            
            allocate (notlost(n_uprp, n_upar, n_psi),STAT = istat)
            if(istat .ne. 0) 
     .           call allocate_error("notlost, sub rdc_multi",0,istat)
            call ibcast(notlost,1,SIZE(notlost))
            
            allocate (dnotlost(n_uprp, n_upar, n_psi),STAT = istat)
            if(istat .ne. 0) 
     .           call allocate_error("dnotlost, sub rdc_multi",0,istat)
            call bcast(dnotlost,one,SIZE(notlost))
            
            
c Open file to read lost regions, specified on an equispaced
c uprp,upar grid at each radial flux surface.
            
            iunit=14
            open(unit=iunit,file='lossfile(k)',status='unknown')
            read(iunit,309) n_uprp !number of uprp points
            read(iunit,309) n_upar !number of upar points
            read(iunit,309) n_R !number of radial (PSI) points
            read(iunit, 3310) vc_cgs !Max vel on the grid
            read(iunit, 3310) upar_min,upar_max !Generally, -1. +1.
                                                ! Max uprp is assumed
                                                ! to = max upar.
                                                ! min uprp=0.
            read(iunit,310) (((notlost(i_prp, i_par, i_psi),
     +            i_prp =0,n_uprp-1), i_par =0,n_upar-1), i_psi =1,n_R)
! notlost=0 for loss cone (if at least one i_phase or i_posn is lost),
! notlost=1 otherwise
            close(iunit)        ! Close file 
!-----------------------------------------------------------------
 3310       format(1p6e18.10)
 309        format(10i10)
c BH remove 1p descriptor for gfortran 4..5.1 (20101208) 
c[Nonsense for i decriptor].          
c 310        format(1p6i2)                                                              
 310        format(6i2)                                                                 


            uprp_min=0.d0
            uprp_max=upar_max

c     Check enorm for notlost data:
      
            vc_cgs2=vc_cgs*vc_cgs
            gammac=sqrt(1.d0+vc_cgs2/clite2)
            if ((gammac-1.d0).lt.1.e-6) then
               enormc=0.5d0*fmass(1)*(vc_cgs)**2/ergtkev
            else
               enormc=(gammac-1.d0)*fmass(1)*clite2/ergtkev
            endif
            if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
            write(*,*)'Enorm in losscone, gammac =',enormc,gammac
            endif

c.................................................................
c     The rho_a radial mesh and associated coeffc will be reduced
c     by a factor of 2, if lrz.eq.n_psi/2, enabling cql3d to run
c     on half the number of flux surfaces used in the full-wave code.
c     Check if n_psi.eq.lrz.   
c     If not, check n_psi/2.eq.lrz.  If so, omit every 2nd radial point
c       of the diffusion coeff grid and coeffs, and reset n_psi.
c       This enables factor of 2 reduction of the cql3d radial mesh.
c       The cql3d and du0u0_input radial meshes are assumed to be
c         the same (or close).
c       Future modification:  Interpolate the du0u0_input radial mesh
c                             to the cql3d radial mesh.
c     Check if n_psi.eq.lrz:  if not, STOP
     
      if (n_psi.ne.lrz .and.
     +    int((n_psi+1.1)/2) .eq. lrz) then

c.................................................................
c     Adjust the radial mesh to use every second one ==> 32 radial
c     points from 64, or 64 radial points from 128.  
c.................................................................

         n_psi=n_psi/2
         if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
         write(*,*)'losscone: n_psi,lrz=',n_psi,lrz
         endif
         do i_psi=1,n_psi
            rho_a(i_psi)=rho_a(2*i_psi)
            do i_upar=1,n_upar
               do i_uprp=1,n_uprp
                  notlost(i_uprp,i_upar,i_psi)=
     &                notlost(i_uprp,i_upar,2*i_psi)
               enddo
            enddo
         enddo
         
         if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
         write(*,*)'losscone: rho_a, no. elements/2:',
     +        (rho_a(i),i=1,n_psi)
         endif
                 
      endif

      if (n_psi.ne.lrz) then
         WRITE(*,*)'losscone:  n_psi.ne.lrz'
         STOP
      endif


c.................................................................
c  The velocity grids are assumed to be equi-spaced.  Renormalize
c  to account for possible vc_cgs.ne.vnorm

c  Setup arrays of initial upar0,uprp0, momentum-per-rest-mass at the
c  minimum B-field point.
c.................................................................

      dupar=(upar_max-upar_min)/(n_upar-1)*vc_cgs/vnorm
      do i=1,n_upar-1
         upar(i)=upar_min+(i-1)*dupar
      enddo
      
      duprp=(uprp_max-uprp_min)/(n_uprp-1)*vc_cgs/vnorm
      do i=1,n_uprp-1
         uprp(i)=uprp_min+(i-1)*duprp
      enddo
      dupdup=dupar*duprp

c.................................................................
c  Bilinear interpolate notlost onto the cql3d gone(i,j,k,ll) grid.
c  [Coding modified from rdc_multi.]
c  If any quantity remains on the grid which is less than 1.d0
c  then set it to zero.
c  If vc_cgs.lt.vnorm, and values of gone(i,j=jc,k,ll)=0.d0, then
c  set gone (i,j=jc+1:jx,k,ll)=0.d0.
c.................................................................

      do i_psi=1,n_psi
         do i_upar=1,n_upar
            do i_uprp=1,n_uprp
               if (notlost(i_uprp,i_upar,i_psi).eq.0) then
                  dnotlost(i_uprp,i_upar,i_psi)=-one !i.e., 
                                                     !gone-like lost 
               else
                  dnotlost(i_uprp,i_upar,i_psi)=zero !i.e., 
                                                     !gone-like not lost
               endif
            enddo
         enddo
      enddo

      do ll=1,lrz
      do j=1,jx
         do i=1,iy
            xupar=x(j)*coss(i,ll)
            xuprp=x(j)*sinn(i,ll)
            if (xupar.lt.upar_max .and. xupar.gt.upar_min) then
               if (xuprp.lt.uprp_max .and. xuprp.ge.uprp_min) then
                  ipar=(xupar-upar_min)/dupar+1
                  iprp=(xuprp-uprp_min)/duprp+1
                  w00=(upar(ipar+1)-xupar)*(uprp(iprp+1)-xuprp)/dupdup
                  w01=(upar(ipar+1)-xupar)*(xuprp-uprp(iprp))/dupdup
                  w10=(xupar-upar(ipar))*(uprp(iprp+1)-xuprp)/dupdup
                  w11=(xupar-upar(ipar))*(xuprp-uprp(iprp))/dupdup
                  wsum=w00+w01+w10+w11
                  gone(i,j,k,ll)=(w00*dnotlost(iprp,ipar,ll)
     +                     +w01*dnotlost(iprp+1,ipar,ll)
     +                     +w10*dnotlost(iprp,ipar+1,ll)
     +                     +w11*dnotlost(iprp+1,ipar+1,ll))
                  if (gone(i,j,k,ll).ne.zero) gone(i,j,k,ll)=-one
               endif
            endif
         enddo
      enddo
      enddo

c  If vc_cgs.lt.vnorm, haven't covered the whole iy,jx grid.
c  If last point ga vc_cgs indicates prompt loss, then extend
c    this at constant pitch angle to the jx edge of the mesh.
         jc=luf(vc_cgs/vnorm,x,jx)-1 !Will be less than jx if
                                  !vc_cgs/vnorm.lt.x(jx-1)

         if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
         write(*,*)'losscone: jx,jc=',jx,jc
         endif
         
         do i=1,iy
         if (jc.lt.jx .and. gone(i,jc,k,ll).eq.-one) then
            do j=jc,jx
               gone(i,j,k,ll)=-one
            enddo
         endif
         enddo
         



         endif      ! on lossmode(k) 
         !===========================================================

         if (lossmode(k) .eq. "mirrsnk") call lossorbm(ephicc,k)

 100  continue ! k  species
 
cBH160529      if (lr_.eq.1) STOP  !This is for stop after extended o/p
      return
      end
      
      
