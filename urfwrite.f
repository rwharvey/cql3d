c
c
      subroutine urfwrite
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine performs the converse operation on urfread:
c     i.e., a disk file "rayop" in written or updated with ray tracing
c     data.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

CMPIINSERT_IF_RANK_NE_0_RETURN

      krf=0
      if (lh.eq."enabled") then
        krf=krf+1
        open(unit=20,file='raylh',status='old')
        call urfwrite_(krf,20)
      endif
      if (ech.eq."enabled") then
        krf=krf+1
        open(unit=23,file='rayech',status='old')
        call urfwrite_(krf,23)
      endif
      if (fw.eq."enabled") then
        krf=krf+1
        open(unit=24,file='rayfw',status='old')
        call urfwrite_(krf,24)
      endif
      return
      end
