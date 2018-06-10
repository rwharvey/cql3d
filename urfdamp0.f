c
c
      subroutine urfdamp0(irfpwr,kopt)
      implicit integer (i-n), real*8 (a-h,o-z)
cyup      save

c..................................................................
c     This routine calculates the power flowing along each ray
c     accounting for the damping, starting from an edge value
c     equal to (1/2)**irfpwr*(initialized input value).
c
c     kopt  .ne. 2: Modify local power according to power already
c     damped along the ray
c     =   2: Keep delpwr as initially given. Used for diagnostic rerun
c
c     Reminder:
c        mrf= number of rf "types"
c        mrfn= number of rf "modes" (sum over mrf of the nharms())
c        irfn(1:mrf)=  rf mode index (in 1:mrfn) of the lowest
c                      harmonic for each wave type
c        irfm(1:mrf)=  number of modes for each distinct wave
c                      type (where there may be more than one
c                      identical wave file applied to multiple
c                      species.
c        irffile(1:mrf)= "separate", "combine1", or "combine2",
c                      indicating separate, first of series of
c                      identical files applied to different
c                      species, continuation of the identical files
c                      applied to different species.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

cyup      character*8 ifirst
cyup      data ifirst/"first"/

      call bcast(powrft(1),zero,lrzmax)
      call bcast(powurf(0),zero,nmodsa+1)
      call bcast(powurfc(0),zero,nmodsa+1)
      call bcast(powurfl(0),zero,nmodsa+1)
      call bcast(powrf(1,1),zero,lrza*nmodsa)
      call bcast(powrfc(1,1),zero,lrza*nmodsa)
      call bcast(powrfl(1,1),zero,lrza*nmodsa)
      call bcast(powurfi(0,0),zero,(lrza+1)*(nmodsa+1))

      call bcast(urfpwr(1,1,1),zero,nrayelts*nrayn*mrfn)
cyup      if (ifirst.eq."first") then
         call bcast(urfpwrc(1,1,1),zero,nrayelts*nrayn*mrfn)
         call bcast(urfpwrl(1,1,1),zero,nrayelts*nrayn*mrfn)
cyup      endif
      call bcast(scalurf,one, nrayelts*nrayn*mrfn)
cBH120815      call bcast(salphac,zero,nrayelts*nrayn*mrfn)

      do krfmode=1,mrfn ! loop over wave modes.
c..................................................................
c     Apply collisional or added linear damping to first harmonic
c     of separate wave types.
c     Indicate whether this wave mode, krfmode, corresponds to the 
c     lowest harmonic for a ray type (lk=1, if yes).
c..................................................................
        lk=0
        do ktype=1,mrf ! loop over wave types.    
           if (krfmode.eq.irfn(ktype)) lk=1
        enddo
c..................................................................
c     Get damping coefficients for all the ray elements
c..................................................................
        if (lk.eq.1) call urfdampa(krfmode) !collisional or added linear
             ! damping: urfpwrc(is,iray,krfmode) urfpwrl(is,iray,krfmode)
      enddo

ccc        call bcast(da(1,0),zero,iyjxp1) ! YuP: Why needed? No effect.
      if (n.ge.nondamp) then
ccc      do 6 ll=1,lrz ! Loop over flux surfaces
      do 5 krfmode=1,mrfn ! loop over wave modes.
ccc          call tdnflxs(lmdpln(ll)) !-> Determine l_, lr_, etc.
          if (n.gt.0 .and. urfdmp.eq."secondd") then
             call urfdamp2(krfmode)
          else
             call urfdamp1(krfmode)
          endif
 5    continue ! krfmode=1,mrfn ! loop over wave modes.
ccc 6    continue ! ll
      endif

CMPIINSERT_BARRIER 
CMPIINSERT_BCAST_URFPWR
 

cyup      ifirst="notfirst"

c..................................................................
c
c     If nharms(1:mrf).gt.0, then power flow along delpwr for
c     the given wave type is calculated only for the lowest harmonic.
c     The resulting delpwr is copied to the harmonics
c     at end of subroutine.
c
c     Explanation of some of following variables:
c
c     delpwr(is=1:nrayelt(iray,krf),iray=1:nray(krf),krf=1:mrfn), where
c     mrfn=sum_over_wave_types(nharms), is the power flowing
c     along each ray (ergs/sec), for each mode krf.
c
c     urfpwr(is,iray,krf) is the fractional power absorbed by each 
c     ray element, i.e., fraction of power flowing into the ray
c     element which is absorbed in the ray element. 
c     In nharms=0,1 case, this is for each mode krf.
c     In nharms.gt.1 case, this is for each harmonic, krf=1:nharms.
c
c     urfpwrc(is,iray,krf) fractional collisional power absorption,
c     passed from ray tracing code.
c
c     urfpwrl(is,iray,krf) additional linear power absorption passed
c     from ray tracing code.
c
c     powrf(lr,krf) is the sum of urf power deposited from each mode
c     (or harmonic for nharms.gt.1) in each radial bin, divided
c     by bin volume. (watts/cm**3).  
c
c     powrfc(lr,krf) collisional power/volume (watts/cm**3).      
c
c     powrfl(lr,krf) additional linear power/volume (watts/cm**3).
c
c     powrft(lr) power per volume (watts/cm**3) summed over modes
c     or harmonics, dur to urf, collisional and add. linear abs.
c
c     powurf(0) total rf power summed over modes and radius
c     powurf(1:mrfn) total power summed over radius, for each
c     mode or harmonic.  (watts)
c
c     powurfi(0:lr,0:mrfn) partial powers integrated up to rho(lr).
c     powurf(krf)=powurfi(lrzmax,krf).    (watts)
c
c     powurfc(0) total coll. urf power summed over modes and radius 
c     powurfc(1:mrfn) total coll. urf power summed over radius, for
c     each mode or harmonic. (watts)
c
c     powurfl(0:mrfn) ditto for urf linear power. 
c..................................................................


c..................................................................
c     Calculate time-dependent pwrscalet, based on 
c     nurftime,urftime(),pwrscale1()
c..................................................................

         if (nurftime.le.0) then
            pwrscalet=one
         else
            itme=1
            do jtm=1,nurftime
               if (timet.ge.urftime(jtm)) itme=jtm
            enddo
            itme1=itme+1
            if (itme.lt.nurftime) then
               pwrscalet=pwrscale1(itme)
     1              +(pwrscale1(itme1)-pwrscale1(itme))
     1              /(urftime(itme1)-urftime(itme))
     1              *(timet-urftime(itme))
            else
               pwrscalet=pwrscale1(nurftime)
            endif
         endif



      do 500 ktype=1,mrf ! loop over wave types.
         krf=irfn(ktype) !index of lowest harmonic in 1:mrfn (all modes)
                         !for each ktype
c         write(*,*)'urfdamp0: ktype, krf, pwrscale(ktype)',
c     +                        ktype, krf, pwrscale(ktype)
c.......................................................................
c     Update power flowing along the rays using run indicator irffile,
c     set in urfread.
c
c     If irffile()="separate" and nharms(krf).gt.0,
c     then sum power absorption over harmonics
c     for the given wave type.
c
c     If irffile()="combine1"
c     then sum power absorption over all modes,
c     and set pwr absorption data for irffile()="combine2" files
c     equal to the "combined1" data.
c.......................................................................

c        maxmodes=nharms(ktype)           
c        if (irffile.eq.'combined') then
c           maxmodes=mrfn
c           if (ktype.ge.2) goto 501 !I.E., finished with the loop
c        endif ! if (irffile.eq.'combined')

        if (irffile(ktype).eq."separate") then
           maxmodes=nharms(ktype)
        elseif (irffile(ktype).eq."combine1") then
              maxmodes=irfm(ktype)
        elseif (irffile(ktype).eq."combine2") then
           goto 500  !i.e., finished with this distinct wave type
        endif


        if(kopt.ne.2)then ! modify local power, and determine nrayelt0
        pwrfrac=0.5**irfpwr*pwrscale(ktype)*pwrscalet
        powray=0.d0
        do iray=1,nray(krf)
          nrayelt0(iray,krf)=nrayelt(iray,krf)
          delpwr(1,iray,krf)=pwrfrac*delpwr0(iray,krf) ! is=1
          powray=powray+delpwr(1,iray,krf) ! Not used???
          if(nrayelt(iray,krf).ge.2)  then ! A normal long ray
            do is=2,nrayelt(iray,krf)
               ism=is-1
               urfpwrt=0.d0
               do kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
                  kkk=krf+kk-1 !mode index
                  urfpwrt= urfpwrt + urfpwr(ism,iray,kkk) 
     +                    +urfpwrc(ism,iray,kkk)+urfpwrl(ism,iray,kkk)
               enddo ! kk=1,MAX(1,maxmodes)
               ! Redefine delpwr in this element through the previous:
               delpwr(is,iray,krf)=delpwr(ism,iray,krf)*(1.d0-urfpwrt)
               !!!!!!! here some of delpwr may become <0, because urfpwr>1.
               !(but normally, urfpwr<<1.)
               !write(*,*)'kopt=1: iray,is, urfpwrt=', iray,is, urfpwrt
               !if(urfpwrt.ge.1.d0) pause
               if(delpwr(is,iray,krf).le.1.d-6*delpwr(1,iray,krf))then
                 ! The left-over power in the ray element is negligible
                 delpwr(is:nrayelt(iray,krf),iray,krf)=em100
                 nrayelt0(iray,krf)=ism ! The last element where power is still not small
                 if(istarts(iray,krf).eq.1)  istarts(iray,krf)=-1
                 ! Finished with this ray, go to next ray:
                 go to 22  !!!!!!! this jump also happens when delpwr<0
               endif
            enddo ! is
 22       continue   ! for skipping all the rest ray-elements      
          else ! A short(one point) ray :
            delpwr(1,iray,krf)=em100
          endif
        enddo ! iray
        endif ! kopt=1
        
        
        if(kopt.eq.2)then !do NOT modify delpwr, only determine nrayelt0
        do iray=1,nray(krf)
          nrayelt0(iray,krf)=nrayelt(iray,krf)
          if(nrayelt(iray,krf).ge.2)  then ! A normal long ray
            do is=2,nrayelt(iray,krf)
               ism=is-1
               if(delpwr(is,iray,krf).le.1.d-6*delpwr(1,iray,krf))then
                 ! The left-over power in the ray element is negligible
                 nrayelt0(iray,krf)=ism ! The last element where power is still not small
                 if(istarts(iray,krf).eq.1)  istarts(iray,krf)=-1
                 go to 23  ! Finished with this ray, go to next ray:
               endif
            enddo ! is
 23         continue   ! for skipping all the rest ray-elements      
          endif ! if nrayelt(iray,krf)=1, do nothing
        enddo ! iray
        endif ! kopt=2
        
        
        ! Now sum-up deposited RF powers, sum over rays and ray-elements
        do 20  iray=1,nray(krf)
            do 15  is=1,nrayelt0(iray,krf) ! up to the element where delpwr is still not too small
              lrl=lloc(is,iray,krf) ! Only for ZOW: 
              !ZOW: one ray-element contributes to only one radial bin.
              do 17 kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
                kkk=krf+kk-1
                powrf(lrl,kkk)=powrf(lrl,kkk)
     1            +delpwr(is,iray,krf)*urfpwr(is,iray,kkk)*1.d-7
                !YuP[02/08/2015] Correction for powrfc and powrfl: is->is-1
                !Should there be a factor of scalurf present? 
                powrfc(lrl,kkk)=powrfc(lrl,kkk)
     1            +delpwr(is,iray,krf)*urfpwrc(is,iray,kkk)*1.d-7
                powrfl(lrl,kkk)=powrfl(lrl,kkk)
     1            +delpwr(is,iray,krf)*urfpwrl(is,iray,kkk)*1.d-7
 17           continue
 15         continue ! is=1,nrayelt0(iray,krf)
 20     continue ! do iray=1,nray(krf)

        ! Total powers (sum over radial index)
        do 30 l=1,lrzmax
            do 19 kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
              kkk=krf+kk-1
              powurfi(l,kkk)=powurfi(l-1,kkk)+powrf(l,kkk)
              powurfi(l,0)=powurfi(l,0)+powurfi(l,kkk)
              powrf(l,kkk)=powrf(l,kkk)/dvol(l)  !Converted to pwr dens
              powrft(l)=powrft(l)+powrf(l,kkk) !Pwr dens summed on modes
              !-> Coll.damping, if any:
              ! Total coll.power (integral over rho and sum over modes):
              powurfc(0)=powurfc(0)+powrfc(l,kkk) 
              ! Integral over rho, for each mode (Watts):
              powurfc(kkk)=powurfc(kkk)+powrfc(l,kkk)
              powrfc(l,kkk)=powrfc(l,kkk)/dvol(l) ! watt/cm^3
              !-> Linear damping, if any:
              powurfl(0)=powurfl(0)+powrfl(l,kkk)
              powurfl(kkk)=powurfl(kkk)+powrfl(l,kkk)
              powrfl(l,kkk)=powrfl(l,kkk)/dvol(l)
 19         continue ! kk=1,maxmodes
 30     continue ! l=1,lrzmax
        do 29 kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
          kkk=krf+kk-1
          powurf(kkk)=powurfi(lrzmax,kkk)
 29     continue ! kk=1,maxmodes
      powurf(0)=powurfi(lrzmax,0)

      if (maxmodes.gt.1) then
        do 600 kk=2,maxmodes
          kkk=krf+kk-1
          call dcopy(nrayn*nrayelts,delpwr(1,1,krf),1,delpwr(1,1,kkk),1)
          do iray=1,nrayn
             nrayelt0(iray,kkk)=nrayelt0(iray,krf)
          enddo
 600    continue ! kk=1,maxmodes
      endif

      
 500  continue ! do ktype=1,mrf
 
 
 
 
 501  continue ! After end of do 500 ! (Not used)

c..................................................................
c     Print out min/max , for diagnostic purposes 
c      do kkk=1,5
c      urfb_max=zero
c      urfb_min=ep100
c      do iray=1,nrayn
c         do is=1,nrayelts
c            yuri= urfpwr(is,iray,kkk) 
c            urfb_max1=max(urfb_max,yuri)
c            urfb_min1=min(urfb_min,yuri)
c            if (urfb_max1.gt.urfb_max) then
c               urfb_max=urfb_max1
c            endif
c            if (urfb_min1.lt.urfb_min) then
c               urfb_min=urfb_min1
c            endif
c         enddo
c      enddo
c      write(*,*)'powurf-e->:',kkk,powurf(kkk),powurfl(kkk)
c      enddo
c...................................................................

      return
      end
