c
c
      subroutine soup0
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      include 'comm.h'

      if (soucoord.eq."cart") return
      do 10 kk=1,ngen
        do 11 m=1,nso
          do 12 j=1,jx
            tam1(j)=-(x(j)-xem1(kk,m,lr_))**2/xem2(kk,m,lr_)
            sovt(j,kk,m,lr_)=exp(tam1(j))
 12       continue
 11     continue
 10   continue
      return
      end
