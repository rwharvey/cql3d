c
c
      subroutine tdtrwtl
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     Computes the Chang-Cooper weights for the transport
c     calculation.
c..............................................................

      include 'param.h'
      include 'comm.h'

      do 1 k=1,ngen
        do 15 j=1,jx
          do 10 l=0,lrors-1
            do 20 i=1,iymax
              dl(i,j,k,l)=.5
 20         continue
 10       continue
 15     continue
 1    continue

      return
      end
