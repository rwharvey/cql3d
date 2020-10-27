!
!
      real*8 function bsl(jj,kk,ll)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      include 'comm.h'
!      dimension f_(0:iy+1,0:jx+1,ngen,lrz)
!      f_ is passed in common
      allocatable :: bsl_s(:,:,:)

!.........................................................................
!     This routine is used to compute the skewing effect at the p/t bndry
!     due to the bootstrap effect. It is not used if bootcalc="disabled",
!     or if advnce="explicit" or if lrz=1
!     Subroutine bsl is for the lower theta tp-bndry at i=itl.  
!     There is also a separate subroutine bsu for the upper itu bndry.
!.........................................................................
! YuP[2019-12-12] migrated corrections from CQL3D-FOW version
!     YuP [2014-05-22] Corrected dfdr definition, and also rban.
!     See comment in bsu.f.

      bsl=0.
      if (bootcalc.eq."disabled") return
      if (implct.ne."enabled" .or. lrz.eq.1 .or. n.lt.nonboot) return
      if (jj.le.2) return ! YuP[07-24-2014] added: 
                          ! df/dtheta should be 0 for j=1 and 2,
                          ! so do not apply bootstrap modifications.
                          ! Tests: difference in 4th-5th digit

      if (.NOT. ALLOCATED(bsl_s)) then
        allocate( bsl_s(0:jx+1,lrz,ngen) )
        call bcast(bsl_s,zero,(jx+2)*lrz*ngen)
      endif
      
      qb_mc=bnumb(kk)*charge*bthr(ll)/(fmass(kk)*clight)

      jjj=max(jj,1)   ! to limit from below, for x(jjj)
      jjj=min(jjj,jx) ! to limit from above, for x(jjj)
      rban= 0.45*cursign*x(jjj)*coss(itl_(ll),ll)*vnorm/qb_mc
      ! YuP [140522] factor 0.45 added in rban

      if (bootcalc.eq."method1") then

        if (n.eq.nonboot.or.bootupdt.eq."enabled") then
          if (ll.eq.1) then
!            dfdr=(f_(itl_(ll+1),jj,kk,ll+1)-f_(itl_(ll),jj,kk,ll))/
!     &      (rz(ll+1)-rz(ll))
             p1=rpcon(ll+1)-rpcon(ll)
             p2=rpcon(ll+2)-rpcon(ll+1)
             p3=rpcon(ll+2)-rpcon(ll)
!             p1= (rpcon(ll+1)-rpcon(ll)   + rmcon(ll)-rmcon(ll+1)  )*0.5
!             p2= (rpcon(ll+2)-rpcon(ll+1) + rmcon(ll+1)-rmcon(ll+2))*0.5
!             p3= (rpcon(ll+2)-rpcon(ll)   + rmcon(ll)-rmcon(ll+2)  )*0.5
             dfdr=-(p1+p3)/(p1*p3)*f_(itl_(ll),jj,kk,ll)
     &            +p3/(p1*p2)*f_(itl_(ll+1),jj,kk,ll+1)
     &            -p1/(p2*p3)*f_(itl_(ll+2),jj,kk,ll+2)
          elseif (ll.eq.lrz) then
!            dfdr=(f_(itl_(ll),jj,kk,ll)-f_(itl_(ll-1),jj,kk,ll-1))/
!     &      (rz(ll)-rz(ll-1))
             p1=rpcon(ll-1)-rpcon(ll-2)
             p2=rpcon(ll)-rpcon(ll-1)
             p3=rpcon(ll)-rpcon(ll-2)
!             p1= (rpcon(ll-1)-rpcon(ll-2) + rmcon(ll-2)-rmcon(ll-1))*0.5
!             p2= (rpcon(ll)-rpcon(ll-1)   + rmcon(ll-1)-rmcon(ll)  )*0.5
!             p3= (rpcon(ll)-rpcon(ll-2)   + rmcon(ll-2)-rmcon(ll)  )*0.5
             dfdr=+p2/(p1*p3)*f_(itl_(ll-2),jj,kk,ll-2)
     &            -p3/(p1*p2)*f_(itl_(ll-1),jj,kk,ll-1)
     &            +(p2+p3)/(p2*p3)*f_(itl_(ll),jj,kk,ll)
          else ! 1<ll<lrz
             p1=rpcon(ll)-rpcon(ll-1)
             p2=rpcon(ll+1)-rpcon(ll)
             p3=rpcon(ll+1)-rpcon(ll-1)
!             p1= (rpcon(ll)-rpcon(ll-1)   + rmcon(ll-1)-rmcon(ll)  )*0.5
!             p2= (rpcon(ll+1)-rpcon(ll)   + rmcon(ll)-rmcon(ll+1)  )*0.5
!             p3= (rpcon(ll+1)-rpcon(ll-1) + rmcon(ll-1)-rmcon(ll+1))*0.5
!!             dfdr=-p2/(p1*p3)*f_(itl_(ll-1),jj,kk,ll-1)
!!     &            -(p1-p2)/(p1*p2)*f_(itl_(ll),jj,kk,ll)
!!     &            +p1/(p2*p3)*f_(itl_(ll+1),jj,kk,ll+1)
            dfdr=(f_(itl_(ll+1),jj,kk,ll+1)-f_(itl_(ll-1),jj,kk,ll-1))/
     &      p3 !-> same results as above dfdr, up to 4th digit
          endif
          bsl_=-bootsign*dfdr*rban
          bsl_s(jj,ll,kk)=bsl_  ! save into array
        else  ! n>nonboot
          bsl_=bsl_s(jj,ll,kk)  ! use saved array
        endif ! n=nonboot

!     Limit the jump at trapped-passing boundary to 0.2*(distn functn).
!     If getting larger values, should probably at least consider
!     method2.  (BH).
!     YuP: usually the jump is smaller: small change if skipped.
        bsum=abs(0.2*f_(itl_(ll),jj,kk,ll))
        bsl=sign(one,bsl_)*min(bsum,abs(bsl_))
          bsl=bsl_ ! uncomment to skip the above.
          ! For fast ions, the jump could be larger than 0.2*f,
          ! so, restricting the jump to 0.2f results in under-estimation
          ! of current.  On the other hand, at small rho and large energies,
          ! there should be no jump at all - pinch orbits disappear.

      elseif (bootcalc.eq."method2") then

!       Assuming positive current here (should be generalized).
        rrr=rpcon(ll)-rban ! rban carries the sign of vpar
        !write(*,*)'bsl: rmcon(1),rpcon(ll)=',ll,rmcon(1),rpcon(ll)
        if(rrr.gt.MAXVAL(rlimiter)) then 
           ! out of plasma -> lost orbit
           f_irrr= 1.383896526737250d-87
           !f_irrr= f_(itl_(lrz),jj,kk,lrz) ! YuP[2019-12-13] Not too diff.
           ! can be a large jump because of loss cone:
           bsl= f_irrr - f_(itl_(ll),jj,kk,ll)
        elseif(rpcon(ll)-abs(rban).lt.rmcon(1)) then 
           f_irrr= 1.383896526737250d-87
           ! No jump, for now. Orbits are too large 
           !(maybe no pinch orbits at all at small rho & large energy):
           bsl= 0.d0 ! f_irrr - f_(itl_(ll),jj,kk,ll)
        else ! rmag<rrr<rpcon(lrzmax)
           rrr=max(rpcon(1),rrr) ! in case rmag<rrr<rpcon(1)
           call lookup(rrr,rpcon(1),lrzmax,weightu,weightl,irrr)
           f_irrr= weightl*f_(itl_(irrr-1),jj,kk,irrr-1)
     &            +weightu*f_(itl_(irrr),jj,kk,irrr)
           bsl= f_irrr - f_(itl_(ll),jj,kk,ll)
        endif
        bsl_s(jj,ll,kk)=bsl ! save into array

      endif

      return     
      end function bsl

 