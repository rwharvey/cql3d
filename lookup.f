c
c
      subroutine lookup(x,xarray,length,weightu,weightl,lement)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine uses luf_bin to do a table look up. 
c     Then it interpolates to gain a bit of accuracy.
c     x is the argument; xarray is the monotonic array; length
c     is the length of the array. lement is the first index such
c     that xarray(lement).gt.x.
c     weightu/weightl are the upper/lower weights for xarray
c     in  x(lement:lement+1).
c     If x falls outside bounds of the array, weightl/weightu set
c     to give constant value at relevant bounds of the array.
c..................................................................

      save
      dimension xarray(*)
      
      if(x.lt.xarray(1)) then
         ! YuP: this case happens quite often:
         ! orbit gets to R smaller than R of 1st flux surface.
         ! Suppressing print-out.
cyup         if ((abs(xarray(1)-x)/
cyup     1        max(abs(xarray(1)),abs(x))).gt.em12) then
c           write if outside roundoff limits:
cyup            write(*,*)'WARNING in lookup lement.lt.1'
cyup            write(*,*)'the code will set lement=2'
cyup            write(*,*)'x,xarray(1)=',x,xarray(1)
cyup         endif
         lement=2
         weightl=1.d0
         weightu=0.d0
         goto 10
      endif

      if(x.ge.xarray(length)) then 
         ! YuP: this case:
         ! orbit gets outside of equilpsi(lrz), the surface #lrz.
         if ((abs(xarray(length)-x)/
     1        max(abs(xarray(length)),abs(x))).gt.em12) then
c           write if outside roundoff limits:
            !write(*,*)'WARNING in lookup lement.gt.length'
            !write(*,*)'the code will set lement=length'
            !write(*,*)'x,xarray(length)=',x,xarray(length)
            ! This may happen if an orbit is out of LCFS (or R-grid)
            ! Suppress printout.
         endif
         lement=length
         weightl=0.d0
         weightu=1.d0
         goto 10
      endif
      
      lement=luf_bin(x,xarray,length)
      if(lement.le.1) stop 'lookup: lement.le.1' ! should never happen
      weightl=(xarray(lement)-x)/(xarray(lement)-xarray(lement-1))
      weightu=1.-weightl
      
 10   return
      end
c
c
      subroutine lookup_tdf(x,xarray,length,weightu,weightl,lement)
      implicit integer (i-n), real*8 (a-h,o-z)

c        lookup_tdf is special version of subroutine lookup, with
c        printout of warning about exceeding table limits turned
c        off, avoiding excess printout under normal conditions
c        for tdfinterp.

c..................................................................
c     This routine uses luf_bin to do a table look up. 
c     Then it interpolates to gain a bit of accuracy.
c     x is the argument; xarray is the monotonic array; length
c     is the length of the array. lement is the first index such
c     that xarray(lement).gt.x.
c     weightu/weightl are the upper/lower weights for xarray
c     in  x(lement:lement+1).
c     If x falls outside bounds of the array, weightl/weightu set
c     to give constant value at relevant bounds of the array.
c..................................................................

      save
      dimension xarray(*)
      
      if(x.lt.xarray(1)) then
         ! YuP: this case happens quite often:
         ! orbit gets to R smaller than R of 1st flux surface.
         ! Suppressing print-out.
cyup         if ((abs(xarray(1)-x)/
cyup     1        max(abs(xarray(1)),abs(x))).gt.em12) then
c           write if outside roundoff limits:
cyup            write(*,*)'WARNING in lookup lement.lt.1'
cyup            write(*,*)'the code will set lement=2'
cyup            write(*,*)'x,xarray(1)=',x,xarray(1)
cyup         endif
         lement=2
         weightl=1.d0
         weightu=0.d0
         goto 10
      endif

      if(x.ge.xarray(length)) then 
         ! YuP: this case:
         ! orbit gets outside of equilpsi(lrz), the surface #lrz.
c         if ((abs(xarray(length)-x)/
c     1        max(abs(xarray(length)),abs(x))).gt.em12) then
c           write if outside roundoff limits:
c            write(*,*)'WARNING in lookup lement.gt.length'
c            write(*,*)'the code will set lement=length'
c            write(*,*)'x,xarray(length)=',x,xarray(length)
c         endif
         lement=length
         weightl=0.d0
         weightu=1.d0
         goto 10
      endif
      
      lement=luf_bin(x,xarray,length)
      if(lement.le.1) stop 'lookup: lement.le.1' ! should never happen
      weightl=(xarray(lement)-x)/(xarray(lement)-xarray(lement-1))
      weightu=1.-weightl
      
 10   return
      end
