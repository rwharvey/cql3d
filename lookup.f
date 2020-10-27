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
      !implicit integer (i-n), real*8 (a-h,o-z)

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

      !save
      integer length,lement
      real*8 x,weightu,weightl
      real*8 xarray(length)
      
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
      
c      write(*,*)'lookup_tdf[befor.luf_bin] x,xarray=',x,xarray(1:length)
      lement=luf_bin(x,xarray,length)
c      write(*,*)
c     + 'lookup_tdf[aft.luf_bin] lement,xarray(lement-1),xarray(lement)',
c     +                          lement,xarray(lement-1),xarray(lement)
      if(lement.le.1)then !YuP[2018-06-29] modification for x<xarray(1) case
         !YuP/was:   stop 'lookup: lement.le.1'
         !YuP[2018-06-29] New version:
         !Note: luf_bin=1 happens only when px.lt.parray(1)  [see luf_bin]
         !Set lement to 2, and attribute all weight to lement-1 point
         lement=2
         weightl=1.d0 ! weight factor for (lement-1) point
         weightu=0.d0 ! weight factor for (lement) point
         return
      endif
      
      ! weight factor for (lement-1) point: 
      weightl=(xarray(lement)-x)/(xarray(lement)-xarray(lement-1))
      ! weight factor for (lement) point: 
      weightu=1.-weightl
      !Example: lement=2 means that x.ge.xarray(1)
      !   and x<xarray(2) (strictly less)
      
 10   return
      end
