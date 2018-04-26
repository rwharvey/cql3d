c
c
      subroutine lossorbm(ephi,ksp)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Mirror loss cone calculation.
c..................................................................

      save
      include 'param.h'
      include 'comm.h'

      vp2=2.*ephi*ergtkev/fmass(ksp)
      do 20 i=1,iy
        do 10 j=1,jx
          v02=(vnorm*x(j))**2
          vp02=v02*(coss(i,l_)**2+(1.-psimx(lr_))*sinn(i,l_)**2)
          vd2=vp2-vp02
          if(vd2.lt.0.d0) then
             temp1(i,j)=-one
          else
             temp1(i,j)=zero
          endif
 10     continue
 20   continue

c..................................................................
c     For case "mirrsnk" f_ in loss hole near x=0.
c..................................................................

      if (lossmode(ksp) .eq. "mirrorcc") then
        call dcopy(iyjx2,temp1(0,0),1,gone(0,0,ksp,indxlr_),1)
      elseif (lossmode(ksp) .eq. "mirrsnk") then
        do 30 i=1,iy
          do 40 j=1,jx
            if(gone(i,j,ksp,lr_).ge.-em90  .and. temp1(i,j).eq.0.) then
               gone(i,j,ksp,indxlr_)=zero
            else
               gone(i,j,ksp,indxlr_)=-one
            endif
 40       continue
 30     continue
      elseif (lossmode(ksp) .eq. "mirrsnk1") then
        do j=2,jx
           if (x(j).ge.xsink) go to 50
        enddo
 50     continue
        jsink=min(j,jx)
        do i=1,iy
          do j=1,jx
            if(gone(i,j,ksp,lr_).ge.-em90  .and. temp1(i,j).eq.0.) then
               gone(i,j,ksp,indxlr_)=zero
            elseif(gone(i,j,ksp,lr_).ge.-em90  .and. temp1(i,j).eq.-one 
     1             .and. j.lt.jsink ) then
               gone(i,j,ksp,indxlr_)=zero
            else
               gone(i,j,ksp,indxlr_)=-one
            endif
         enddo
      enddo

      endif
      return
      end
