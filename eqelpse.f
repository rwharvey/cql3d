
!
!
      subroutine eqelpse
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!
      include 'param.h'
      include 'comm.h'
      character*8 ifirst


!..................................................................
!     This routine generates an ad-hoc psi (an "solution" to
!     the equilibrium calculation)  with elliptical contours of psi
!     centered at rmag with ellipticity = ellptcty. The funtional form
!     of psi(E,lr_) where E**2=z**2+(r-rmag)**2/(1.-ellptcty**2) is
!     arbitrary and is handled by function eqfn. The psi function
!     is normalized so that the poloidal field at r=rmag, z=radmin
!     is bth.
!..................................................................

!..................................................................
!     Create the z,r meshes..
!..................................................................

      data ifirst /"first"/
      if (lr_.ne.lrzmax) return
      if (ifirst.eq."first") then
        zst=zbox/(nnz-1)
        rst=rbox/(nnr-1)
        er(1)=rboxdst
        ez(1)=-zbox*.5
        do 10 nn=2,nnr
          er(nn)=er(nn-1)+rst
 10     continue
        do 11 nn=2,nnz
          ez(nn)=ez(nn-1)+zst
 11     continue
        zmag=zero
        ifirst="notfirst"
      endif

!..................................................................
!     Set up the psi array. First determine the normalization
!     constant from the poloidal field constraint.
!..................................................................

      fhp=eqfn(radmin*1.0001,one)
      fhm=eqfn(radmin*.9999,one)
      deriv=(fhp-fhm)/(radmin*2.e-4)
      scalfct=-bth*rmag/deriv

!..................................................................
!     scalfct is the factor, now for the psi array..
!..................................................................

      do 20 ir=1,nnr
        do 21 iz=1,nnz
          epsi(ir,iz)=sqrt(ez(iz)**2+(er(ir)-rmag)**2/(1.-ellptcty**2))
 21     continue
 20   continue
      do 30 ir=1,nnr
        do 31 iz=1,nnz
          epsi(ir,iz)=eqfn(epsi(ir,iz),scalfct)
 31     continue
 30   continue
      return
      end subroutine eqelpse




!======================================================================
!======================================================================

      subroutine eq_miller_set
!     Setup (R,Z) grids and pol.flux array over (R,Z) for Miller equilibrium.
!     See equations (36-37) of
!     R.L. Miller et al., "Noncircular, finite aspect ratio, local
!     equilibrium model", Phys. Plasmas, Vol. 5, No. 4, April 1998.

!     INPUT: rbox,zbox,rboxdst   from namelist/setup/
!     (size of R,Z area for setting the grids, and innermost R)
!     (Note: zmag is assumed =0)

!     OUTPUT: er() and ez() arrays/grids, 
!     and epsi(ir,iz) array for pol.flux
!     Store them in comm.h
     
      implicit integer (i-n), real*8 (a-h,o-z)
!      save
!
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE
      
      parameter(nworka=3*nnra+1)
      character*8 ntitle,dat   ! Added for g95 compiler, Urban 110708
      dimension ntitle(5),workk(nworka)

      character*20 eq_miller_eqdskin ! generated data is saved to this file
      eq_miller_eqdskin='eqdskin_miller_saved' 

      !--- Define commonly used values through eq_miller_*** namelist:
      !......................... Stored in comm.h:
      rmag=   eq_miller_rmag ! Magnetic axis: major radius coord [cm]
      zmag=   eq_miller_zmag ! Magnetic axis: vertical coord [cm]
      cursign=eq_miller_cursign ! Sign of Plasma Current [+1. or -1.]
      psimag= eq_miller_psimag  ! Pol.flux at magn.axis [cgs] 
      psilim= eq_miller_psilim  ! Pol.flux at LCFS [cgs] 
      psi_n=  eq_miller_psi_n ! n and m powers for PSI(r) profile as in 
      psi_m=  eq_miller_psi_m ! for PSI(r)= psilim + (psimag-psilim)*(1-(r/a)^n)^m
      !......................... Stored in name_decl.h:
      btor=   eq_miller_btor    ! Tor field at Geom. center of LCFS [Gauss]
      bsign=  sign(1.d0,btor)   ! Sign of btor
      radmin= eq_miller_radmin  ! Plasma minor radius [cm]
      !......................... local/input:
      delta_e=eq_miller_deltaedge ! Triangularity of LCFS (at r=radmin)
      akappa= eq_miller_kappa  ! Vertical elongation (const for all surfaces)
      drr0=   eq_miller_drr0   ! dR0/dr  we assume Shafr.shift=-drr0*r
      
      !--- Derive other values, based on model (Miller+our_assumptions):
      !......................... Stored in name_decl.h:
      radmaj= rmag +drr0*radmin ! R at Geom. center of LCFS [cm]
                                ! Note: here, we assume Shafr.shift=-drr0*r
      !......................... local/input:
      fpsieq= btor*radmaj       ! B*R assumed constant (no depen. on r)
      
      
      rbox= 2*radmin*1.5  ! +-25% outside of LCFS
      zbox= 2*radmin*akappa*1.5  ! +-25% outside of LCFS
      ! Innermost/lowest R for the box:
      rboxdst=  max( radmaj-radmin*1.25, (radmaj-radmin)/2 ) 
      ymideqd=0.d0  !is the vertical shift of the rectangular box
    
CMPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)' '
      WRITE(*,*)'--------------------------'
      WRITE(*,*)'eq_miller_set: rmag, zmag, radmaj, radmin,rboxdst',
     &           rmag, zmag, radmaj, radmin,rboxdst
      WRITE(*,*)'eq_miller_set: cursign, bsign, psimag, psilim',
     &           cursign, bsign, psimag, psilim
CMPIINSERT_ENDIF_RANK

      ! See subr. eq_miller() for definition of surfaces and fields.
      !------------------------------------------------------------------

!..................................................................
!     Create the (R,Z) grids.
!..................................................................
      zst=zbox/(nnz-1)  ! nnr and nnz are in comm.h
      rst=rbox/(nnr-1)
      zmag=zero
      er(1)=rboxdst   ! Innermost/lowest R for the box.
      ez(1)=-zbox*.5  ! Lowest Z for the box; Assuming zmag=0.
      do 10 nn=2,nnr
          er(nn)=er(nn-1)+rst  ! R grid [cm]
 10   continue
      do 11 nn=2,nnz
          ez(nn)=ez(nn-1)+zst  ! Z grid [cm]
 11   continue
      ezmin=ez(1)   ! cgs
      ezmax=ez(nnz)
      ermin=er(1)
      ermax=er(nnr)

!..................................................................
!     Set up the pol.flux array. 
!..................................................................
      do 20 ir=1,nnr
          rr=er(ir) ! R
      do 21 iz=1,nnz
          zz=ez(iz) ! Z
          ! Convert (R,Z)->(rminor,polang)
         !Early version, only valid for ellipse (not D-shape):
         !call eq_miller_RZ_to_rpol(rmag,zmag, delta_e,akappa,drr0, 
!     &                              rr,zz, rminor,polang) 
          !New version, valid for D-shape:
          call eq_miller_RZ_to_rpol_fnd(rmag,zmag, radmin,
     &                               delta_e,akappa,drr0,rr,zz, 
     &                               rminor,polang) 
          ! Get PSI value:
          call eq_miller(rminor, polang, 
     &        rmag, zmag, cursign, bsign, radmin, psimag,psilim,
     &           psi_n,psi_m,
     &        delta_e, akappa, drr0, fpsieq,
     &        dpsidr, Rs,Zs, PSI, BR, BPHI, BZ)
          epsi(ir,iz)=PSI  ![cgs]
 21   continue
 20   continue

      
      ! Set up minor-r radial grid (nnv) and psiar(nnv), fpsiar(nnv), etc.
      nnv= int(max(nnr,nnz)/2) ! r-grid is about half-size of total R-grid
      do irho=1,nnv
          rminor= radmin*(irho-1)/(nnv-1) ! minor r= [0; radmin]  [cm]
          polang= 0.d0
          ! Get PSI value:
            call eq_miller(rminor, polang, 
     &           rmag, zmag, cursign, bsign, radmin, psimag,psilim,
     &           psi_n,psi_m,
     &           delta_e, akappa, drr0, fpsieq,
     &           dpsidr, Rs,Zs, PSI, BR, BPHI, BZ)
          psiar(irho)= PSI    ! Assumed linear dep. on r (for now)
          fpsiar(irho)=fpsieq ! Just a const (no dep. on r)
          prar(irho)=  0.d0   ! Undefined value, for now
          ffpar(irho)= 0.d0   ! Undefined value, for now
          ppar(irho)=  0.d0   ! Undefined value, for now
          qar(irho)=   0.d0   ! Undefined value, for now
      enddo

      ! Set LCFS : go along surface at r=radmin
      ncontr=100 ! 100 points along surface is good enough
      allocate(rcontr(ncontr),STAT=istat) 
      allocate(zcontr(ncontr),STAT=istat) 
      do ipol=1,ncontr
            rminor= radmin
            polang= twopi*(ipol-1)/(ncontr-1)
            ! Get Rs,Zs coords at this surface:
            call eq_miller(rminor, polang, 
     &           rmag, zmag, cursign, bsign, radmin, psimag,psilim,
     &           psi_n,psi_m,
     &           delta_e, akappa, drr0, fpsieq,
     &           dpsidr, Rs,Zs, PSI, BR, BPHI, BZ)
           rcontr(ipol)= Rs
           zcontr(ipol)= Zs
           !write(*,*)'polang,R,Z=',polang,Rs,Zs
      enddo

      ! Set limiter = LCFS+10%
      nlimiter=100
      allocate(rlimiter(nlimiter),STAT=istat) 
      allocate(zlimiter(nlimiter),STAT=istat) 
      do ipol=1,ncontr
            rminor= radmin*1.10  ! 10% outside of LCFS
            polang= twopi*(ipol-1)/(ncontr-1)
            ! Get Rs,Zs coords at this surface:
            call eq_miller(rminor, polang, 
     &           rmag, zmag, cursign, bsign, radmin, psimag,psilim,
     &           psi_n,psi_m,
     &           delta_e, akappa, drr0, fpsieq,
     &           dpsidr, Rs,Zs, PSI, BR, BPHI, BZ)
           Rlim= min(Rs,er(nnr)) ! not to exceed the grid
           Rlim= max(Rlim,er(1)) ! not smaller than left border of R-grid
           Zlim= min(Zs,ez(nnz)) ! not to exceed the grid
           Zlim= max(Zlim,ez(1)) ! not smaller than lower border of Z-grid
           rlimiter(ipol)= Rlim
           zlimiter(ipol)= Zlim
           !write(*,*)'polang,R,Z=',polang,Rs,Zs
      enddo
         

      ! Supposed to be plasma current, but not really used, only for sign of current:
      toteqd=cursign  
      
      ! Not used, but fill-in the lines in data file, to keep the structure:
      psimx1=0.d0
      psimx2=0.d0
      xax1=0.d0
      xax2=0.d0
      zax1=0.d0
      zax2=0.d0
      psisep=0.d0
      xsep=0.d0
      ysep=0.d0
      zsep=0.d0
      
      ntitle(1)='Miller  '
      ntitle(2)='Equilib '
      ntitle(3)='cursign,'
      ntitle(4)='nnr,nnz,'
      ntitle(5)='nnv '
      dat= '='
      ipestg= int(cursign) !?

CMPIINSERT_IF_RANK_EQ_0
      open(unit=10,file=eq_miller_eqdskin, STATUS='unknown')
        
      write(10,110)(ntitle(i),i=1,5),dat,ipestg,nnr,nnz,nnv ! (6a8,4i4)
      write(10,120) rbox/100,zbox/100,radmaj/100,rboxdst/100,
     &                ymideqd/100   ! All converted to [m]
      write(10,120) rmag/100,zmag/100,psimag/1.e8,psilim/1.e8,btor/1.e4
      write(10,120) toteqd,psimx1,psimx2,xax1,xax2
      write(10,120) zax1,zax2,psisep,xsep,zsep
      ! psimx1,psixm2,xax1,xax2,zax1,zax2,psisep,xsep,ysep - OBSOLETE.
      write(10,120) (fpsiar(i)/1.e6,i=1,nnv) ! converted to mks
      write(10,120) (prar(i),i=1,nnv)
      write(10,120) (ffpar(i),i=1,nnv)
      write(10,120) (ppar(i),i=1,nnv)
      write(10,120) ((epsi(i,j)/1.e8,i=1,nnr),j=1,nnz)
      write(10,120) (qar(i),i=1,nnv)        
      write(10,8210) ncontr, nlimiter
      write(10, 8200)  
     .     (rcontr(i)/100, zcontr(i)/100, i = 1,ncontr) ! [m]
      write(10, 8200)
     .     (rlimiter(i)/100, zlimiter(i)/100, i = 1,nlimiter) ![m]
 8210   format (2i5)
 8200   format(5e16.9)
 110    format(6a8,4i4)
 120    format(5e16.9)
 250    format( (5(e21.14)) )
 251    format(5i5)

      close(unit=10) ! data is saved
      
      WRITE(*,*)'eq_miller_set: Miller eqdskin is generated'
      WRITE(*,*)'--------------------------'
      WRITE(*,*)' '
CMPIINSERT_ENDIF_RANK
      
      ! This part follows subr.equilib
      ! Set up spline array for the f =(R*BTOR) subroutine,etc.
      i1p(1)=4
      i1p(2)=4
      call coeff1(nnv,psiar,fpsiar,d2fpsiar,i1p,1,workk)
      call coeff1(nnv,psiar,ffpar,d2ffpar,i1p,1,workk)
      call coeff1(nnv,psiar,prar,d2prar,i1p,1,workk)
      call coeff1(nnv,psiar,ppar,d2ppar,i1p,1,workk)
      call coeff1(nnv,psiar,qar,d2qar,i1p,1,workk)
      nfp=nnv ! Save into comm.h, for usage in terp1()
      ibd(1)=4
      ibd(2)=4
      ibd(3)=4
      ibd(4)=4
      call coeff2(nnr,er,nnz,ez,epsi,epsirr,epsizz,epsirz,nnra,ibd,
     &      wkepsi)

      return
      end subroutine eq_miller_set


!======================================================================
!======================================================================
      subroutine eq_miller_RZ_to_rpol(rmag,zmag, delta_e,akappa,drr0, 
     &                                R,Z, rminor,polang) 
      !NOT USED ANYMORE: Early version, only valid for ellipse (not D-shape)
      ! Convert (R,Z)->(r,polang)
      !!! Only valid for delta_e=0 (triangularity=0), 
      ! Instead, use subr. eq_miller_RZ_to_rpol_fnd()
      implicit integer (i-n), real*8 (a-h,o-z)
CMPIINSERT_INCLUDE
      
      ! Assumptions:
      ! Major radius [cm] along magnetic surface:
      ! Rs= rmag + drr0*r + r*cos(polang)  ! Assumes Shafr.shift ~ r
      ! Zs= zmag + akappa*r*sin(polang) ! Vertical coord [cm] along magnetic surface
      
      if(delta_e.ne.0.d0) then
         WRITE(*,*)'eq_miller_RZ_to_rpol: Not ready for delta_e > 0.'
         stop
      endif
      
      rmr= R-rmag
      zmz= (Z-zmag)/akappa
      drr2= drr0**2  ! Supposed to be <<1
      
      DDD= rmr**2*drr2 + (1.d0-drr2)*(rmr**2+zmz**2)
      ! If no Shafranov shift then drr0=0, then DDD=(rmr**2+zmz**2)
      if(DDD.lt.0.d0)then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'eq_miller_RZ_to_rpol: DDD<0'
CMPIINSERT_ENDIF_RANK
         stop
      endif
      
      rminor= ( -rmr*drr0 +sqrt(DDD) )/(1.d0-drr2) ! positive root
      
      if(rminor.lt.0.d0)then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'eq_miller_RZ_to_rpol: rminor<0',rminor
CMPIINSERT_ENDIF_RANK
         rminor=max(rminor,1.d-8)
         !stop
      endif
      
      ! Find polang from tg(polang) = (Z-zmag)/kappa/(R-rmag-drr0*r)
      polang= atan2(zmz,rmr-drr0*rminor)  ! [-pi,pi]
      if (polang.lt.0.d0) polang=polang+6.28318530717959 ! [0,2pi)

      return
      end subroutine eq_miller_RZ_to_rpol


!======================================================================
!======================================================================
      subroutine eq_miller_RZ_to_rpol_fnd(rmag,zmag, a,
     &                               delta_e,akappa,drr0,R,Z, 
     &                               rminor,polang) 
      !YuP[2020-07-13] For each given (R,Z) point,
      !find (rminor,polang) coordinate that can be used as input for 
      !subr.eq_miller(r, polang,...)
      !This subroutine combines the definition for Rs(r,polang)
      ! and for Zs(r,polang) to eliminate r variable and to 
      ! make a function [see function eq_miller_rr(pol)] 
      ! of only polang variable. The root of eq_miller_rr(pol)=0
      ! for each given (R,Z) point is saved as pol_fnd
      ! (see below). The value of r_fnd is readily evaluated
      ! from Zs(r,pol_fnd)=Z  equation 
      ![basically, r_fnd= (Z-zmag)/(kappa*sin(pol_fnd)) ]
      implicit none
CMPIINSERT_INCLUDE

      ! INPUT:
      real*8 rmag, zmag ! Major R of magn.axis, Vertical Z of magn.axis
      real*8 a ! minor radius at LCFS (rminor=a)
      real*8 delta_e ! Triangularity at plasma edge (r=a). 
                     ! At r, we assume delta(r)= delta_e*(r/a)
      real*8 akappa  ! Vertical elongation. Const. in COGENT? Then s_kappa==0
      real*8 drr0    ! dR0/dr !we assume Shafr.shift=-drr0*r
      real*8 R,Z  ! (R,Z) coordinates of a given point [cm]
      
      ! OUTPUT:
      real*8 rminor,polang ! Minor radius [cm] and generalized pol.angle [rad]

!     local variables:
      real*8 sint,cost, tpbs, ra
      real*8 pol_fnd,r_fnd,sint_fnd,delta_fnd,beta_fnd, Rs_fnd,Zs_fnd
      real*8 pi, costpbs, sintpbs, det,det1
      real*8 rr_target,zz_target
      real*8 pol_mn, pol_mx, pol_eps
      integer iter, itermx
      real*8 err0, err
      integer ier, nsig
      
      external eq_miller_rr ! external F(pol.angle),
                     ! for subr.zbrent(F,err,nsig,pol_mn,pol_mx,itermx,ier)
                     
      !To transfer some parameters 
      !between subr.eq_miller_RZ_to_rpol_fnd() and func.eq_miller_rr()
      common/eq_miller_f/rsin,a00,delta00,drr00,rmag00,r00
      real*8 rsin,a00,delta00,drr00,rmag00,r00
      ! For transfer through common :
      a00=a
      delta00=delta_e
      drr00=drr0
      rmag00=rmag
      r00=R
      
      pi=3.141592653589793d0
      
      itermx= 100 ! Max number of allowed iterations in zbrent;
                  ! Usually it takes 6-10 iterations.
      err=  1.d-9*a ! For zbrent(), conversion criteria (can set to 0.)
      nsig= 8  ! For zbrent(), conversion criteria
      ier=0 !in/out diagnostics
                        
      rr_target= R-rmag
      zz_target= Z-zmag
      
      !-1-> Set min/max values (boundaries) for pol.angle range.
      !     The external function eq_miller_rr(pol.angle) should
      !     have opposite values at these boundaries
      pol_eps= 1.d-8 ![rad]
      
      if(Z-zmag.gt. 1.d-8*a)then    ! Z>zmag
        if(R.gt.rmag)then
           pol_mn=   -pol_eps
           pol_mx= pi-pol_eps
        else ! R<rmag
           pol_mn=    pol_eps
           pol_mx= pi+pol_eps
        endif
      elseif(Z-zmag.lt.-1.d-8*a)then ! Z<zmag
        if(R.gt.rmag)then
           pol_mn=   pi+pol_eps
           pol_mx= 2*pi+pol_eps
        else ! R<rmag
           pol_mn=   pi-pol_eps
           pol_mx= 2*pi-pol_eps
        endif
      else ! Z=zmag (with machine accuracy)
        if(R.gt.rmag)then
          rminor= (R-rmag)/(drr0+1.d0)
          polang= 0.d0
          return ! Done
        elseif(R.lt.rmag)then
          rminor= (R-rmag)/(drr0-1.d0) !Note: typically drr0=-0.1...-0.3
          polang= pi
          return ! Done
        else ! R=rmag (and also here Z=zmag)
          rminor=1.d-8 ! Should be 0., but let's stay away from 0
          polang= 0.d0 ! just any value.
          return ! Done
        endif
      endif ! Z
      
      !-2-> Find zero of eq_miller_rr(pol.angle)
      iter=itermx
      rsin= (Z-zmag)/akappa  != r*sin(pol.angle) [r and pol.angle 
      !are unknown here, but this combination is a known value]
      call zbrent(eq_miller_rr,err,nsig,pol_mn,pol_mx,iter,ier)
      !On output, the root is put into pol_mx, 
      !Also, on output, iter contains the actual number of iterations
      pol_fnd= pol_mx
      sint_fnd= sin(pol_fnd) 
      if(abs(sint_fnd).gt.0.1*abs(rsin/a))then
        r_fnd= rsin/sint_fnd ! Can be r_fnd>>a, but it's ok.
      else ! sin(pol_fnd)~0
        r_fnd= 10.0*a
      endif
      rminor= r_fnd   !r is FOUND for given (R,Z) point
      polang= pol_fnd !polang FOUND
      
      !The rest is optional [could be commented]
      !-3-> Verify that (R_fnd,Z_fnd) is close to the target point (R,Z)
      delta_fnd= min(delta_e*(r_fnd/a),delta_e) !delta(r) Triangularity -
      ! limited by value at r=a (in case r gets larger than a)
      beta_fnd= asin(delta_fnd)   ! = 'x' in Miller paper 
      Rs_fnd= rmag + drr0*r_fnd + r_fnd*cos(pol_fnd + beta_fnd*sint_fnd)
      Zs_fnd= zmag + akappa*r_fnd*sint_fnd 
      err0= sqrt( (Rs_fnd-R)**2 + (Zs_fnd-Z)**2)/a
CMPIINSERT_IF_RANK_EQ_0
      if((ier.ne.0).or.(iter.gt.10).or.(err0.gt.1.d-6))then
        WRITE(*,*)' '
        WRITE(*,*)'eq_miller_RZ_to_rpol_fnd: TARGET (R,Z)=',
!BH     &  (R,Z)
     &  R,Z
        WRITE(*,*)'iter,ier=',iter,ier
        WRITE(*,*)'err0=',err0
        pause
      endif
CMPIINSERT_ENDIF_RANK
      
      return
      end subroutine eq_miller_RZ_to_rpol_fnd
      
!======================================================================
!======================================================================
      real*8 function eq_miller_rr(pol)
      !YuP[2020-07-13] See comments in the beginning 
      !of subr.eq_miller_RZ_to_rpol_fnd
      implicit none
      real*8 pol ! INPUT (to be iterated by subr.zbrent)
      !
      real*8 sint, r, delta, beta, costpbs, Rsint !local
      real*8 a, delta_e, drr0, rmag, Rtarget
      
      !To transfer some parameters 
      !between subr.eq_miller_RZ_to_rpol_fnd() and func.eq_miller_rr()
      common/eq_miller_f/rsin,a00,delta00,drr00,rmag00,r00 !INPUT here
      real*8 rsin,a00,delta00,drr00,rmag00,r00

      a=a00
      delta_e=delta00
      drr0=drr00
      rmag=rmag00
      Rtarget=r00
      
      sint= sin(pol) 
      if(abs(sint).gt.0.1*abs(rsin/a))then
        r= rsin/sint ! Can be r_fnd>>a, but it's ok.
      else ! sin(pol_fnd)~0
        r= 10.0*a
      endif
      delta= min(delta_e*(r/a),delta_e) !delta(r) Triangularity -
      ! limited by value at r=a (in case r gets larger than a)
      beta= asin(delta)   ! = 'x' in Miller paper 
      !Note: Rs= rmag + drr0*r + r*cos(pol + beta*sint)
      !Note: Zs= zmag + akappa*r*sint
      costpbs= cos(pol + beta*sint) 
      ! Rn*sin(pol)==
      Rsint= rmag*sint +rsin*(drr0+costpbs) !we assume Shafr.shift=-drr0*r
      eq_miller_rr= Rsint - Rtarget*sint 
      
      return
      end function eq_miller_rr

!======================================================================
!======================================================================
      subroutine eq_miller(r, polang, 
     &           rmag, zmag, cursign, bsign, a, psimag,psilim,
     &           psi_n,psi_m,
     &           delta_e, akappa, drr0, fpsieq,
     &           dpsidr, Rs,Zs, PSI, BR, BPHI, BZ)
!     Returns: Magnetic field components
!     given by equations (36-37) of
!     R.L. Miller et al., "Noncircular, finite aspect ratio, local
!     equilibrium model", Phys. Plasmas, Vol. 5, No. 4, April 1998.
!  Almost exact copy of COGENT version: MillerBlockCoordSysF.ChF  
!  The difference is in units: COGENT uses [Tesla, meters],
!  while CQL3D uses [Gauss, cm].
!  Also: Here, we are using factor cursign, 
!  which is the direction of plasma current in CQL3D.
!  Direction of Bpol is defined as
!  Bpol= -cursign*|grad(psi)|/R
!  If plasma current is in positive-phi direction 
!  (counter-clockwise if viewed from top at tokamak; cursign=+1.0)
!  then Bpol is directed down at the outermost point of flux surface
!  and upward at the innermost point.
!  Similarly, the sign of tor.field is defined by bsign.

      implicit none !integer (i-n), real*8 (a-h,o-z)
      
!     OUTPUT:
      real*8 Rs,Zs   ! Coordinates along pol.flux surface
      real*8 dpsidr  ! |d(psi)/dr| such that 
                     !            Bpol*R= -cursign*|d(psi)/dr|*|grad(r)|
      real*8 PSI ! poloidal flux [maxwell] (includes factor 2*pi)
      real*8 BR, BZ, BPHI   ! Components of magnetic field  [Gauss]
             ! along major radius R, vertical coord. Z and tor.angle phi
!     INPUT:
      real*8 r,polang  ! Minor radius [cm] and Poloidal angle [rad].
      real*8 delta_e ! Triangularity at plasma edge (r=a). 
                     ! At r, we assume delta(r)= delta_e*(r/a)
      real*8 akappa  ! Vertical elongation. Const. in COGENT? Then s_kappa==0
      real*8 drr0    ! dR0/dr !we assume Shafr.shift=-drr0*r
      real*8 fpsieq  ! Btor*R [Gauss*cm]
      real*8 rmag, zmag, cursign, bsign, a, psimag,psilim, psi_n, psi_m
      
!     Other INPUT:      
!     rmag, zmag   (coords of magnetic axis [cm])
!     cursign  (= +1. or -1.) Sign of plasma current in phi-direction
!     bsign    (= +1. or -1.) Sign of tor.field (+1 for positive in phi-direction)
!     a = plasma minor radius (edge/LCFS)  [cm]
!     psimag = pol.flux at magnetic axis [cgs]
!     psilim = Pol.flux at LCFS [cgs]
!     psi_n,psi_m =   n and m powers for PSI(r) profile as in 
!                   ! PSI(r)= psilim + (psimag-psilim)*(1-(r/a)^n)^m

!     Example of input data in COGENT for a circular equilibrium [SI units]:
!     gksystem.magnetic_geometry_mapping.miller.kappa   = 1.
!     gksystem.magnetic_geometry_mapping.miller.delta   = 0.
!     gksystem.magnetic_geometry_mapping.miller.dpsidr  = 1.71 (-> 1.71e6 cgs)
!     gksystem.magnetic_geometry_mapping.miller.drR0    = 0.
!     gksystem.magnetic_geometry_mapping.miller.s_kappa = 0.0
!     gksystem.magnetic_geometry_mapping.miller.s_delta = 0.0
!     # Coordinates (in units.length) of the magnetic axis (CQL3D: rmag,zmag):
!     gksystem.magnetic_geometry_mapping.miller.origin  = 1.7 0.
!     # "Btor_scale" sets B_tor (units.magnetic_field ) * R (units.length):
!     gksystem.magnetic_geometry_mapping.miller.Btor_scale= 25.65 (-> 25.65e6 cgs)

!     local variables:
      real*8 sint, cost, tpbs, ra,denom,one_ran,delta,beta,dr_dr,dr_dz
      real*8 psi_2
      
      real*8 sk  ! = s_kappa= (r/kappa)*d(kappa)/dr   
                 !Here, it's set to 0, since kappa=const
                 
      real*8 sd  ! = s_delta= (r/sqrt(1-delta^2))*d(delta)/dr  !const in COGENT?
                 !Here, it is set as delta/sqrt(1-delta^2)
                 !to be consistent with delta = delta_e*(r/a)

      sint = dsin(polang)
      cost = dcos(polang)
      
      ra=r/a
      delta= delta_e*(ra)  ! delta(r)    Triangularity
                            ! Note: In COGENT, delta is constant 
                            !(not consistent with def.of S_delta)
      beta= asin(delta)  ! = 'x' in Miller paper 
      tpbs= polang + beta*sint
      
      sd= delta/sqrt(1.d0-delta**2) 
      ! Follows from S_delta(r)= [d(delta)/dr]*r/sqrt(1-delta^2)
      ! and our definition of delta(r)= delta_e*(r/a) see above.
      
      sk=0.d0 ! S_kappa= (r/kappa)*d(kappa)/dr
      
      !dpsidr= -2.0*(psimag-psilim)*r/a**2  !assuming parabolic psi(r) !!!
      !Note: in the above line, using r/a**2 vs ra/a gives different accuracy.
!      if(ra/a-r/a**2.ne.0.d0)
!     &   write(*,*)'ra/a - r/a**2=',ra/a - r/a**2 !~ e-18 difference

      !Now generalized to PSI(r)= psilim +(psimag-psilim)*(1-(r/a)^n)^m
      one_ran= (1.d0-ra**psi_n) !== (1-(r/a)^n)
      dpsidr= -psi_n*psi_m*(psimag-psilim)*one_ran**(psi_m-1)
     &                    *(ra**(psi_n-1))/a
      
      ! Major radius [cm] along magnetic surface:
      Rs= rmag + drr0*r + r*dcos(tpbs)  !we assume Shafr.shift=-drr0*r
      Zs= zmag + akappa*r*sint ! Vertical coord [cm] along magnetic surface
            
      ! Denominator in Eq.(37), Miller:
      denom=  dcos(beta*sint) + drr0*cost 
     &      + (sk - sd*cost + (1.d0 + sk)*beta*cost)*sint*dsin(tpbs) 
     
      ! dr/dR   
      ![checked this for general case dR0/dr#0, d(delta)/dr#0, d(kappa)/dr#0  --YuP]
      dr_dR= cost/denom   !-> cos(polang) for circular 
      
      ! dr/dZ
      dr_dZ= dsin(tpbs)*(1.d0 + beta*cost)/(akappa*denom)  
      !-> sin(polang) for circular
      
      ! Vector in tangential to m.surf. direction, along pos. pol.angle:
      ! (-dr/dZ, +dr/dR)  -> (-sin(polang), +cos(polang)) for circular m.surf.

      !This is basically Eq.(37), Miller :
      !Bp= (dpsidr/Rs)*dsqrt( dr_dZ**2 + dr_dR**2 ) 
      
      !PSI_2=  psimag -(psimag-psilim)*(r/a)**2 !ASSUMING parabola !!!
      !Now generalized to PSI(r)= psilim +(psimag-psilim)*(1-(r/a)^n)^m
      PSI=  psilim +(psimag-psilim)*one_ran**psi_m
      ! Note: one_ran= (1.d0-ra**psi_n)
!      if(PSI_2-PSI.ne.0.d0)
!     &       write(*,*)'PSI_2-PSI=',PSI_2-PSI ! ~ 1e-10 difference
      ! In the above, we forced psi(r) to have max at magn.axis;
      ! but the proper direction of Bpol field is determined by cursign:
      BR=   +cursign*abs(dpsidr/Rs)*dr_dZ 
      BZ=   -cursign*abs(dpsidr/Rs)*dr_dR
      BPHI=  bsign*abs(fpsieq)/Rs

      return
      end subroutine eq_miller

!======================================================================
!======================================================================

