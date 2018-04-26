c
c
      real*8 function bsu(jj,kk,ll)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      include 'comm.h'
C      dimension f_(0:iy+1,0:jx+1,ngen,lrz)
C      f_ is passed in common
      allocatable :: bsu_s(:,:)
      
c.........................................................................
c     This routine is used to compute the skewing effect at the p/t bndry
c     due to the bootstrap effect. It is not used if bootcalc="disabled",
c     or if advnce="explicit" or if lrz=1
c     Subroutine bsu is for the upper theta tp-bndry at i=itu.  
c     There is also a separate subroutine bsl for the lower itl bndry.
c.........................................................................

      bsu=0.
      if (bootcalc.eq."disabled") return
      if (implct.ne."enabled" .or. lrz.eq.1 .or. n.lt.nonboot) return

      if (.NOT. ALLOCATED(bsu_s)) then
        allocate( bsu_s(0:jx+1,lrz) )
        call bcast(bsu_s,zero,(jx+2)*lrz)
      endif

      qb_mc=bnumb(kk)*charge*bthr(ll)/(fmass(kk)*clight)

      jjj=max(jj,1)   ! to limit from below, for x(jjj)
      jjj=min(jjj,jx) ! to limit from above, for x(jjj)
      rban=cursign*x(jjj)*coss(itu_(ll),ll)*vnorm/qb_mc

      if (bootcalc.eq."method1") then

        if (n.eq.nonboot.or.bootupdt.eq."enabled") then
          if (ll.eq.1) then
c            dfdr=(f_(itu_(ll+1),jj,kk,ll+1)-f_(itu_(ll),jj,kk,ll))/
c     1      (rz(ll+1)-rz(ll))
             p1=rz(ll+1)-rz(ll)
             p2=rz(ll+2)-rz(ll+1)
             p3=rz(ll+2)-rz(ll)
             dfdr=-(p1+p3)/(p1*p3)*f_(itu_(ll),jj,kk,ll)
     +            +p3/(p1*p2)*f_(itu_(ll+1),jj,kk,ll+1)
     +            -p1/(p2*p3)*f_(itu_(ll+2),jj,kk,ll+2)
          elseif (ll.eq.lrz) then
c            dfdr=(f_(itu_(ll),jj,kk,ll)-f_(itu_(ll-1),jj,kk,ll-1))/
c     1      (rz(ll)-rz(ll-1))
             p1=rz(ll-1)-rz(ll-2)
             p2=rz(ll)-rz(ll-1)
             p3=rz(ll)-rz(ll-2)
             dfdr=+p2/(p1*p3)*f_(itu_(ll-2),jj,kk,ll-2)
     +            -p3/(p1*p2)*f_(itu_(ll-1),jj,kk,ll-1)
     +            +(p2+p3)/(p2*p3)*f_(itu_(ll),jj,kk,ll)
          else
c            dfdr=(f_(itu_(ll+1),jj,kk,ll+1)-f_(itu_(ll-1),jj,kk,ll-1))/
c     1      (rz(ll+1)-rz(ll-1))
             p1=rz(ll)-rz(ll-1)
             p2=rz(ll+1)-rz(ll)
             p3=rz(ll+1)-rz(ll-1)
             dfdr=-p2/(p1*p3)*f_(itu_(ll-1),jj,kk,ll-1)
     +            -(p1-p2)/(p1*p2)*f_(itu_(ll),jj,kk,ll)
     +            +p1/(p2*p3)*f_(itu_(ll+1),jj,kk,ll+1)
          endif
          bsu_=-bootsign*dfdr*rban
          bsu_s(jj,ll)=bsu_
c          if(jj.eq.30 .and. (n.eq.nonboot .or. n.eq.4) )
c     +      write(*,'(a,2i3,e11.3,4e12.4)')
c     +     'kk,ll, dfdr, f_(ll=2,4,6,8)=',
c     +      kk,ll, dfdr, f_(itu_(2),jj,kk,2),f_(itu_(4),jj,kk,4),
c     +                   f_(itu_(6),jj,kk,6),f_(itu_(8),jj,kk,8)
        else
          bsu_=bsu_s(jj,ll)
        endif

c     Limit the jump at trapped-passing boundary to 0.2*(distn functn).
c     If getting larger values, should probably at least consider
c     method2.  (BH).
        bsum=abs(0.2*f_(itu_(ll),jj,kk,ll))
        bsu=sign(one,bsu_)*min(bsum,abs(bsu_))

      elseif (bootcalc.eq."method2") then

c       Assuming positive current here (should be generalized).
        rrr=rz(ll)-rban
        if (rrr.gt.rz(lrzmax)) rrr=rz(lrzmax)
cBH080714:  rz dimensioned 0:lrza
        call lookup(rrr,rz(1),lrzmax,weightu,weightl,irrr)
cBH051101        if (irrr.eq.lzrmax+1) then
c$$$        if (irrr.ge.lrzmax+1) then
c$$$          irrr=lrzmax+1
c$$$          bsu=f_(itu_(irrr-1),jj,kk,irrr-1)-f_(itu_(ll),jj,kk,ll)
c$$$        else
          bsu=weightl*f_(itu_(irrr-1),jj,kk,irrr-1)
     1        +weightu*f_(itu_(irrr),jj,kk,irrr)
     2        -f_(itu_(ll),jj,kk,ll)
c$$$        endif

        bsu_s(jj,ll)=bsu

      endif

      return     
      end

