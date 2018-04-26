c
c
      subroutine urfavg
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

      alpha=1.
      if (n.gt.2) then
        alpha3=.35
        alpha=(1.-alpha3)/(nstop)*(n-3)+alpha3
      endif
      
      if (n.eq.3) then
ccc        call dcopy(iyjx2*ngen*lrors,f,1,g_,1)
      do l=1,lrors
         do k=1,ngen
            do j=0,jxp1
               do i=0,iyp1
                  g_(i,j,k,l)= f(i,j,k,l) 
               enddo
            enddo
         enddo
      enddo
      endif

      do k=1,ngen
         do 10 l=1,lrors
            do 20 j=1,jx
               do 30 i=1,iy
                  g_(i,j,k,l)=alpha*f(i,j,k,l)+(1.-alpha)*g_(i,j,k,l)
 30            continue
 20         continue
 10      continue
      enddo

      return
      end
