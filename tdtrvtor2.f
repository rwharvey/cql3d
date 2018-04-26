c
c
      subroutine tdtrvtor2(f1,f2,vp,vp_,kopt)
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine takes the function f1 defined on the full
c     velocity mesh and interpolates it onto the velocity mesh
c     on which the transport is done. The new function is f2
c     and it has the same integral as f1 (particle conserving),
c     with vp and vp_ being the weights for f1 and f2 respectively.
c     That is: sum(f1(i)*vp(i)) = sum(f2(i)*vp_(i)), i=1,iy_(l)
c
c     This new version of tdtrvtor is such that two successive calls
c     to tdtrvtor2(f1,f2,.,.,iopt) and then tdtrrtov2(f2,f1b,.,.,iopt)
c     gives f1b = f1 at each point i=1,iy_(l)
c
c     kopt = 1: transformation of distribution function
c     2:       "        of quasilinear term (velsou) (ADI)
c     3:       "        of transport term (spasou)   (ADI)
c
c     Assumes that vp is vp(lrindx(l)) if kopt = 1 and vp(l) otherwise,
c     where l=1,lrors
c..............................................................

      include 'param.h'
      include 'comm.h'

      dimension f1(0:iyp1,0:jxp1,ngen,*)
      dimension f2(0:iyp1,0:jxp1,ngen,*)
      dimension vp(iy,lrza), vp_(iy,lrza)
      dimension vpeff(iy,lrza), vpeff_(iy,lrza)
c.......................................................................
c     copy weight according to option

      do 5 l=1,lrors
        if (kopt .eq. 1) then
          ilr=lrindx(l)
          do 6 i=1,iy_(l)
            vpeff(i,l)=vp(i,ilr)
            vpeff_(i,l)=vp_(i,ilr)
 6        continue
        else
          do 7 i=1,iy_(l)
            vpeff(i,l)=vp(i,l)
            vpeff_(i,l)=vp_(i,l)
 7        continue
        endif
 5    continue

c.......................................................................
c     kopt=3: loop over lrz instead of lrors
c.......................................................................
      leff=lrors
      if (kopt .eq. 3) leff=lrz

      call dcopy (iyjx2*ngen*leff,f1,1,f2,1)
      do 10 k=1,ngen
        do 20 l=1,leff
          itl=itl_(l)
          itu=itu_(l)
          iyy=iy_(l) !-YuP-101215: Don't use iy=; it's in common /params/
                     ! Don't let overwrite the cqlinput value!

          do 30 j=1,jx
            f2(itl-2,j,k,l)=(f1(itl-2,j,k,l)*vpeff(itl-2,l)+
     1        f1(itl-1,j,k,l)*vpeff(itl-1,l)+
     1        f1(itl,j,k,l)*vpeff(itl,l)*.5)/vpeff_(itl-2,l)
            f2(itl+2,j,k,l)=(f1(itl+2,j,k,l)*vpeff(itl+2,l)+
     1        f1(itl+1,j,k,l)*vpeff(itl+1,l)+
     1        f1(itl,j,k,l)*vpeff(itl,l)*.5)/vpeff_(itl+2,l)
            f2(itu-2,j,k,l)=f2(itl+2,j,k,l)
            f2(itu+2,j,k,l)=(f1(itu+2,j,k,l)*vpeff(itu+2,l)+
     1        f1(itu+1,j,k,l)*vpeff(itu+1,l)+
     1        f1(itu,j,k,l)*vpeff(itu,l)*.5)/vpeff_(itu+2,l)

            if (kopt .eq. 1) then
              if (f2(itl-2,j,k,l).ne.zero) then
                f_vtor(j,k,(kopt-1)*6+1,l)=
     +            f1(itl-2,j,k,l)/f2(itl-2,j,k,l)
              else
                f_vtor(j,k,(kopt-1)*6+1,l)=1.0
              endif
              if(f2(itl+2,j,k,l).ne.zero) then
                f_vtor(j,k,(kopt-1)*6+2,l)=
     +            f1(itl+2,j,k,l)/f2(itl+2,j,k,l)
              else
                f_vtor(j,k,(kopt-1)*6+2,l)=1.0
              endif
              if(f2(itu+2,j,k,l).ne.zero) then
                f_vtor(j,k,(kopt-1)*6+3,l)=
     +            f1(itu+2,j,k,l)/f2(itu+2,j,k,l)
              else
                f_vtor(j,k,(kopt-1)*6+3,l)=1.0
              endif
              if (f1(itl-2,j,k,l).ne.zero) then
                f_vtor(j,k,(kopt-1)*6+4,l)=f1(itl,j,k,l)/f1(itl-2,j,k,l)
              else
                f_vtor(j,k,(kopt-1)*6+4,l)=1.0
              endif
              if (f1(itl+2,j,k,l).ne.zero) then
                f_vtor(j,k,(kopt-1)*6+5,l)=f1(itl,j,k,l)/f1(itl+2,j,k,l)
              else
                f_vtor(j,k,(kopt-1)*6+5,l)=1.0
              endif
              if (f1(itu+2,j,k,l).ne.zero) then
                f_vtor(j,k,(kopt-1)*6+6,l)=f1(itl,j,k,l)/f1(itu+2,j,k,l)
              else
                f_vtor(j,k,(kopt-1)*6+6,l)=1.0
              endif
            endif

 30       continue

          do 40 i=itl-1,itl+1
            do 41 j=1,jx
              i_=iyy+1-i
              f2(i,j,k,l)=0.
              f2(i_,j,k,l)=0.
 41         continue
 40       continue

 20     continue
 10   continue

      return
      end
