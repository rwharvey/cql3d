c
c
      subroutine tdtrrtov(f1)
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine interpolates onto the full velocity mesh
c     from the transport velocity mesh.
c     f1 ======> f1,    THAT IS, just fix up itl+/-1,itu+/-1
c..............................................................
      include 'param.h'
      include 'comm.h'

      dimension f1(0:iyp1,0:jxp1,ngen,0:*),denrad(ngen,lrorsa)
      dimension denvel(ngen,lrorsa)

      call tdtrchkd(f1,vpint_,denrad)

cBH070419:   removing special itl,itu treatment for ipacktp=0
        if (ipacktp.eq.3) then

      do 10 k=1,ngen
        do 11 l=1,lrors
          ilr=lrindx(l)
          itl=itl_(l)
          itu=itu_(l)
          do 12 j=1,jx
            fact1=(vpint_(itl-2,ilr)-vpint(itl-2,ilr))*f1(itl-2,j,k,l)
     1        +2.*(vpint_(itl+2,ilr)-vpint(itl+2,ilr))*f1(itl+2,j,k,l)
     1        +(vpint_(itu+2,ilr)-vpint(itu+2,ilr))*f1(itu+2,j,k,l)
            fact2=vpint(itl-1,ilr)*f_lm(j,k,l)+
     1        2.*vpint(itl,ilr)+2.*vpint(itl+1,ilr)*
     1        f_lp(j,k,l)+vpint(itu+1,ilr)*f_up(j,k,l)
            f1(itl,j,k,l)=fact1/fact2
            f1(itu,j,k,l)=f1(itl,j,k,l)
            f1(itl-1,j,k,l)=f_lm(j,k,l)*f1(itl,j,k,l)
            f1(itu+1,j,k,l)=f_up(j,k,l)*f1(itl,j,k,l)
            f1(itl+1,j,k,l)=f_lp(j,k,l)*f1(itl,j,k,l)
            f1(itu-1,j,k,l)=f1(itl+1,j,k,l)
 12       continue
 11     continue
 10   continue

        elseif (ipacktp.ne.0) then
           stop 'STOP in tdtrrtov:  Check ipacktp'
        endif

      call tdtrchkd(f1,vpint,denvel)
      return
      end
