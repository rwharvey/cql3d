


c
      subroutine vlh(action)
      implicit integer (i-n), real*8 (a-h,o-z)
      character*(*) action
      save

c.................................................................
cmnt  subroutine vlh is a fixed diffusion magnitude Landau
cmnt  damping model designed with electron heating in the
cmnt  lower hybrid range of frequencies in mind.
cmnt  subroutine vlh requires: nrf.eq.1 .and. vlhmod.eq."enabled"
c.................................................................
c
      include 'param.h'
      include 'comm.h'

      if (action.ne."setup") return
      
      mrfn=vlhmodes
      
      
      if(ASSOCIATED(cqlb)) then
        ! cqlb-cqlf are already allocated => do nothing 
      else ! Not allocated yet
        allocate(cqlb(iy,jx,lrz,mrfn),STAT=istat)
        allocate(cqlc(iy,jx,lrz,mrfn),STAT=istat)
        allocate(cqle(iy,jx,lrz,mrfn),STAT=istat)
        allocate(cqlf(iy,jx,lrz,mrfn),STAT=istat)
      endif
      cqlb_size=size(cqlb)
      call bcast(cqlb,zero,cqlb_size)
      call bcast(cqlc,zero,cqlb_size)
      call bcast(cqle,zero,cqlb_size)
      call bcast(cqlf,zero,cqlb_size)
      
      
      
      do 35  nmod=1,mrfn
        call bcast (bqlm,zero,iyjx)
        do 30 j=1,jx
          do 20 l=1,lz
            do 10 i=1,iy
              vll=x(j)*cosz(i,l,lr_)*vnorm/gamma(j)
              vprp=x(j)*sinz(i,l,lr_)*vnorm/gamma(j)
cBH091031                 if (l.eq.lz.or.lmax(i,lr_).ne.l) then
cBH091031                    ax=dtau(i,l,lr_)/tau(i,lr_)
cBH091031                 else
cBH091031                    ax=(dtau(i,l,lr_)+dtau(i,l+1,lr_))/tau(i,lr_)
cBH091031                 endif
              if (eqsym.ne."none") then !i.e. up-down symm
                 !if not bounce interval
                 if(l.eq.lz .or. l.ne.lmax(i,lr_)) then
                    ax=dtau(i,l,lr_)/tau(i,lr_)
                 else !bounce interval: additional contribution
                    ax=(dtau(i,l,lr_)+dtau(i,l+1,lr_))/tau(i,lr_)
                 endif
              else    !eqsym="none"
                 if (l.lt.lz_bmax(lr_) .and. l.eq.lmax(i,lr_))
     +                then
                      !trapped, with tips between l and l+1 (above midplane)
                      ax=(dtau(i,l,lr_)+dtau(i,l+1,lr_))/tau(i,lr_)
                      !-YuP  Note: dtau(i,l+1,lr_)=0
                 elseif (l.gt.lz_bmax(lr_) .and. l.eq.lmax(i+iyh,lr_))
     +                then
                      !trapped, with tips between l and l-1 (below midplane)     
                      ax=(dtau(i,l,lr_)+dtau(i,l-1,lr_))/tau(i,lr_) !NB:l-1
                      !-YuP  Note: dtau(i,l-1,lr_)=0
                 else
                      !passing (i<itl), or trapped but with tips at other l;
                      !also, at l=lz_bmax, includes last trapped particle i=itl
                      !(for such particle, lmax(itl)=lz_bmax; see micxinil)
                      ax=dtau(i,l,lr_)/tau(i,lr_)
                 endif
              endif
              bqlm(i,j)=ax*vptb(i,lr_)
     +          *vlhd(vll,vprp,pol(l,lr_),nmod)*cosz(i,l,lr_)**2

              if (bqlm(i,j).eq.zero) go to 10 

              if (vlhprprp(nmod).eq."parallel") then

                cqlb(i,j,indxlr_,1)=cqlb(i,j,indxlr_,1)+bqlm(i,j)*xsq(j)
                cqlc(i,j,indxlr_,1)=cqlc(i,j,indxlr_,1)-bqlm(i,j)*x(j)
     1            *sinn(i,l_)/coss(i,l_)
                cqle(i,j,indxlr_,1)=cqle(i,j,indxlr_,1)-bqlm(i,j)*x(j)
     1            *sinn(i,l_)**2/coss(i,l_)
                cqlf(i,j,indxlr_,1)=cqlf(i,j,indxlr_,1)+bqlm(i,j)
     1            *sinn(i,l_)**3/coss(i,l_)**2

              elseif (vlhprprp(nmod).eq."perp") then

                cqlb(i,j,indxlr_,1)=cqlb(i,j,indxlr_,1)+bqlm(i,j)*xsq(j)
     +            *sinz(i,l,lr_)**2/cosz(i,l,lr_)**2
                cqlc(i,j,indxlr_,1)=cqlb(i,j,indxlr_,1)+bqlm(i,j)*x(j)
     +            *sinn(i,l_)/coss(i,l_)
                cqle(i,j,indxlr_,1)=cqlb(i,j,indxlr_,1)+bqlm(i,j)*x(j)
     +            *sinn(i,l_)**2/coss(i,l_)
                cqlf(i,j,indxlr_,1)=cqlb(i,j,indxlr_,1)+bqlm(i,j)
     +            *cosz(i,l,lr_)**2*sinn(i,l_)
     +            /(bbpsi(l,lr_)*coss(i,l_)**2)

              endif

 10         continue
 20       continue
 30     continue

 35   continue

      do 70 j=1,jx
        do 60 i=itl,iyh
          ii=iy+1-i
          cqlb(i,j,indxlr_,1)=   cqlb(i,j,indxlr_,1)*.5
          cqlc(i,j,indxlr_,1)=   cqlc(i,j,indxlr_,1)*.5
          cqle(i,j,indxlr_,1)=   cqle(i,j,indxlr_,1)*.5
          cqlf(i,j,indxlr_,1)=   cqlf(i,j,indxlr_,1)*.5
          cqlb(ii,j,indxlr_,1)=  cqlb(i,j,indxlr_,1)
          cqlc(ii,j,indxlr_,1)= -cqlc(i,j,indxlr_,1)
          cqle(ii,j,indxlr_,1)= -cqle(i,j,indxlr_,1)
          cqlf(ii,j,indxlr_,1)=  cqlf(i,j,indxlr_,1)
 60     continue
 70   continue

c..................................................................
c     Taper diffusion over last 10 point of velocity, 
c     if ineg="trunc_d"
c..................................................................

      if (ineg.eq."trunc_d") then
        if (jx.le.11) stop 'vlh: Need to set jx>11'
        do 90 j=jx-11,jx
          do 91 i=1,iy
            cqlb(i,j,indxlr_,1)=truncd(j)*cqlb(i,j,indxlr_,1)
            cqlc(i,j,indxlr_,1)=truncd(j)*cqlc(i,j,indxlr_,1)
            cqle(i,j,indxlr_,1)=truncd(j)*cqle(i,j,indxlr_,1)
            cqlf(i,j,indxlr_,1)=truncd(j)*cqlf(i,j,indxlr_,1)
 91       continue
 90     continue
      endif

c.................................................................
c     Plotting diffusion coefficient
c.................................................................
      call vlhbplt

      return
      end
