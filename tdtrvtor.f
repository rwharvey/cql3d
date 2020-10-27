c
c
      subroutine tdtrvtor(f1,f2)
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine takes the distribution f1 defined on the full
c     velocity mesh and interpolates it onto the velocity mesh
c     on which the transport is done. The new distribution is f2
c     and it has the same density as f1 (particle conserving).
c..............................................................

      include 'param.h'
      include 'comm.h'

      dimension f1(0:iyp1,0:jxp1,ngen,0:*)
      dimension f2(0:iyp1,0:jxp1,ngen,0:*)
      dimension denrad(ngen,lrorsa),denvel(ngen,lrorsa)

      call dcopy (iyjx2*ngen*(lrors+1),f1,1,f2,1)

      call tdtrchkd(f1,vpint,denvel)

cBH070419:   removing special itl,itu treatment for ipacktp=0
c      write(*,*)'tdtrvtor: ipacktp', ipacktp
      if (ipacktp.eq.3) then

      do 10 k=1,ngen
        do 20 l=1,lrors
          ilr=lrindx(l)
          itl=itl_(l)
          itu=itu_(l)
          iyy=iy_(l) !-YuP-101215: Don't use iy=; it's in common /params/
                     ! Don't let overwrite the cqlinput value!
          do 30 j=1,jx
            do 40 i=itl-1,itl+1
              i_=iyy+1-i
              f2(i,j,k,l)=0.
              f2(i_,j,k,l)=0.
 40         continue

            if (f1(itl,j,k,l).ne.zero) then
              f_lm(j,k,l)=f1(itl-1,j,k,l)/f1(itl,j,k,l)
              f_lp(j,k,l)=f1(itl+1,j,k,l)/f1(itl,j,k,l)
              f_up(j,k,l)=f1(itu+1,j,k,l)/f1(itl,j,k,l)
            else
              f_lm(j,k,l)=1.0
              f_lp(j,k,l)=1.0
              f_up(j,k,l)=1.0
            endif

            f2(itl-2,j,k,l)=(f1(itl-2,j,k,l)*vpint(itl-2,ilr)+
     1        f1(itl-1,j,k,l)*vpint(itl-1,ilr)+
     1        f1(itl,j,k,l)*vpint(itl,ilr)*.5)/vpint_(itl-2,ilr)
            f2(itl+2,j,k,l)=(f1(itl+2,j,k,l)*vpint(itl+2,ilr)+
     1        f1(itl+1,j,k,l)*vpint(itl+1,ilr)+
     1        f1(itl,j,k,l)*vpint(itl,ilr)*.5)/vpint_(itl+2,ilr)
            f2(itu+2,j,k,l)=(f1(itu+2,j,k,l)*vpint(itu+2,ilr)+
     1        f1(itu+1,j,k,l)*vpint(itu+1,ilr)+
     1        f1(itu,j,k,l)*vpint(itu,ilr)*.5)/vpint_(itu+2,ilr)
            f2(itu-2,j,k,l)=f2(itl+2,j,k,l)
 30       continue
 20     continue
 10   continue

      elseif (ipacktp.ne.0) then
         stop 'STOP in tdtrvtor:  Check ipacktp'
      endif

      call tdtrchkd(f2,vpint_,denrad)
      return
      end
