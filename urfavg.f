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
        !YuP[2019-07-15] Bad idea - to make alpha dependent on nstop.
        !Then, in two test runs - 
        ! one with nstop=10 and another with nstop=15 (example) -
        !the results will be different during first 10 steps already.
        ! This subr. is only needed for g_() below,
        ! which is only used for urfdmp='secondd'.
        ! Better not to use it, for now. So, better use urfdmp='firstd'
      endif
      
      !YuP[2019-07-15] In general, this subr.urfavg is not suitable
      ! for a restart run, because it requires averaging of f()
      ! over 3 time steps. So, if we want a restart run, 
      ! we need to save data from those 3 time steps,
      ! or at least save the g_() function (set below) into *.nc file.
      ! So, in case of restart it is better to use urfdmp='firstd' ;
      ! then this subr. is not called.
      
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
