c
c
      subroutine impnorm(xnorm,a,rhs,nn)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine normalizes the matrix a so that the maximum 
c     coefficient for each equation is of order 1.
c..................................................................
      dimension a(nn)
             
      xnorm=0.d0
      do 1 i=1,nn
        xnorm=xnorm+dabs(a(i))
 1    continue
 
      if (xnorm.gt.0.d0) then
      rhs=rhs/xnorm
      do 2 i=1,nn
        a(i)=a(i)/xnorm
 2    continue
      endif
      
      return
      end
