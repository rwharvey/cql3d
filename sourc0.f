c
c
      subroutine sourc0
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     define source uniquely at v=0
c..................................................................

      include 'param.h'
      include 'comm.h'

      do 10 k=1,ngen
        s=0.
        u=0.
        do 11 i=1,iy
          u=u+source(i,1,k,indxlr_)*cynt2(i,l_)*vptb(i,lr_)
          s=s+cynt2(i,l_)*vptb(i,lr_)
 11     continue
        do 12 i=1,iy
          source(i,1,k,indxlr_)=u/s
 12     continue
 10   continue
      return
      end
