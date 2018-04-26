c
c
      subroutine tdtrrtov2(f1,f2,vp,vp_,kopt)
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine interpolates f1, defined on the transport velocity
c     mesh, onto f2, defined on the full velocity mesh. See comments
c     of tdtrvtor2 for more detail
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
          do 30 j=1,jx
            f2(itl-2,j,k,l)=f1(itl-2,j,k,l) * f_vtor(j,k,(kopt-1)*6+1,l)
            f2(itl+2,j,k,l)=f1(itl+2,j,k,l) * f_vtor(j,k,(kopt-1)*6+2,l)
            f2(itu-2,j,k,l)=f2(itl+2,j,k,l)
            f2(itu+2,j,k,l)=f1(itu+2,j,k,l) * f_vtor(j,k,(kopt-1)*6+3,l)
            f2(itl,j,k,l)=(f2(itl-2,j,k,l)*f_vtor(j,k,(kopt-1)*6+4,l) +
     +        2.*f2(itl+2,j,k,l)*f_vtor(j,k,(kopt-1)*6+5,l) +
     +        f2(itu+2,j,k,l)*f_vtor(j,k,(kopt-1)*6+6,l))*0.25
            f2(itu,j,k,l)=f2(itl,j,k,l)
            f2(itl-1,j,k,l)=(f1(itl-2,j,k,l)*vpeff_(itl-2,l) -
     -        f2(itl-2,j,k,l)*vpeff(itl-2,l) -
     -        0.5*f2(itl,j,k,l)*vpeff(itl,l))/vpeff(itl-1,l)
            f2(itl+1,j,k,l)=(f1(itl+2,j,k,l)*vpeff_(itl+2,l) -
     -        f2(itl+2,j,k,l)*vpeff(itl+2,l) -
     -        0.5*f2(itl,j,k,l)*vpeff(itl,l))/vpeff(itl+1,l)
            f2(itu-1,j,k,l)=f2(itl+1,j,k,l)
            f2(itu+1,j,k,l)=(f1(itu+2,j,k,l)*vpeff_(itu+2,l) -
     -        f2(itu+2,j,k,l)*vpeff(itu+2,l) -
     -        0.5*f2(itl,j,k,l)*vpeff(itl,l))/vpeff(itu+1,l)
 30       continue
 20     continue
 10   continue

      return
      end
