c
c
      subroutine urffflx
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c..................................................................
c     Takes ray data from ray tracing codes and segregates
c     information by flux surface. 
c     Three arrays are calculated:
c     -ncontrib(l) will hold the number of ray elements that contribute
c      to the diffusion coefficient for flux surface l, i.e., in vol(l).
c      lloc(is,iray,krf) is the index of the flux surface to which ray
c      element (is,iray,krf) contributes.  The counting is restricted
c      to wave types, and does not include harmonics.
c     -psiloc(,,) is the value of pol. flux psi associated with the 
c      ray element;     krf is the index of the excitation mode.
c     -lloc(,,) is the radial mesh rho bin.
c     NOTE: psivalm(l) demarks the outer edge of the flux volume associated
c      with flux surface l (vol(l)).
c..................................................................

      include 'param.h'
      include 'comm.h'


      do 1 l=1,lrzmax
 1    tr2(l)=psimag-psivalm(l)

c..................................................................
c     According to the psi value determined for a particular ray element
c     assign it to a given flux surface.
c..................................................................

      call ibcast(ncontrib(1),0,lrzmax)
      icount_outside_lim=0 ! only for printout
      icount_outside_ez=0  ! only for printout
      icount_outside_er=0  ! only for printout
      do 100 krf=1,mrf
        do 20 iray=1,nray(irfn(krf))
          do 30 is=lrayelt(iray,irfn(krf))+1,nrayelt(iray,irfn(krf))
            !-------------------------
            zray=wz(is,iray,irfn(krf))
            if(zray.gt.ez(nnz))then 
              !For a mirror machine: ray can get outside of zbox
              !which defines the border of ez() equilibrium grid
              ![so that zbox=ez(nnz)-ez(1)]
              if (icount_outside_ez.eq.0) then
                if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
                write(*,*)'urffflx: Ray elements outside of ez grid'
                endif
              endif
              icount_outside_ez=icount_outside_ez+1 !for a printout
              if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
                write(*,'(a,i4,2i7)')
     +            'urffflx: zray>ez; iray,is,icount_outside_ez',
     +                               iray,is,icount_outside_ez
              endif
              ! Make an adjustment:
              zray=ez(nnz)
              !This correction is ok for a tokamak, too,
              !although not likely to happen.
            endif
            if(zray.lt.ez(1))then 
              if (icount_outside_ez.eq.0) then
                if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
                write(*,*)'urffflx: Ray elements outside of ez grid'
                endif
              endif
              icount_outside_ez=icount_outside_ez+1 !for a printout
              if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
                write(*,'(a,i4,2i7)')
     +            'urffflx: zray<ez; iray,is,icount_outside_ez',
     +                               iray,is,icount_outside_ez
              endif
              ! Similarly, Make an adjustment:
              zray=ez(1)
            endif
            !-------------------------
            rray=wr(is,iray,irfn(krf))
            if(rray.gt.er(nnr))then 
              !For a mirror machine: ray can get outside of 
              !er() equilibrium grid
              ![so that zbox=ez(nnz)-ez(1)]
              ! Make an adjustment:
              rray=er(nnr)
              !This correction is ok for a tokamak, too,
              !although not likely to happen.
            endif
            if(rray.lt.er(1))then ! this cannot happen, but ok to add.
              ! Similarly, Make an adjustment:
              rray=er(1)
            endif
            !-------------------------
            psiloc(is,iray,irfn(krf))= terp2(rray,zray,nnr,er,nnz,ez,
     1                               epsi,epsirr,epsizz,epsirz,nnra,0,0)
            apsi=psimag-psiloc(is,iray,irfn(krf))
            l=luf(apsi,tr2(1),lrzmax)
cBH090602   Ray elements outside LCFS (rho=1 surface) will be attributed to lrzmax
            if (l.gt.lrzmax) then
               if (icount_outside_lim.eq.0) then
                if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
                 write(*,*)'urffflx: Ray elements outside of rho=1'
                endif
               endif
               icount_outside_lim=icount_outside_lim+1 !for a printout
               if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
               write(*,'(a,i4,2i7)')
     +             'urffflx:l>lrzmax; iray,is,icount_outside_lim',
     +                                iray,is,icount_outside_lim
               endif
               l=lrzmax ! Adjusted
            endif
c$$$            if (l.le.0) then
c$$$               write(*,*)'urffflx:l,lrzmax,k,iray,is',l,lrzmax,k,iray,is
c$$$               go to 30
c$$$            endif
            lloc(is,iray,irfn(krf))=l
cBH090602            if (l.gt.lrzmax) go to 30
            ncontrib(l)=ncontrib(l)+1
 30       continue
 20     continue
 100  continue !  krf=1,mrf

c     Duplicate data for psiloc and lloc into multi-harmonic
c     related arrays.
      do 120 krf=1,mrf
      if (nharms(krf).gt.1) then
        do 110  i=1,nharms(krf)-1
          do 21  iray=1,nray(irfn(krf))
            do 31  is=lrayelt(iray,irfn(krf))+1,nrayelt(iray,irfn(krf))
              psiloc(is,iray,irfn(krf)+i)=psiloc(is,iray,irfn(krf))
              lloc(is,iray,irfn(krf)+i)=lloc(is,iray,irfn(krf))
 31         continue
 21       continue
 110    continue
      endif
 120  continue !  krf=1,mrf

c     Temporary print out checking number of elements at each rad bin:
c      write(*,*)'urffflx:ncontrib(1:lrzmax):',(ncontrib(l),l=1,lrzmax)
      
      return
      end
