      subroutine rdc_bplt(krf)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     plots rf rdcb coefficient as a contour plot on cql3d grid.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*8 pltvlhb
      character*8 pltovlp

      REAL*4 :: R410=10.,R40=0.

      data pltvlhb /'enabled'/
      data pltovlp /'enabled'/


CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return
      if (pltvlhb.ne."enabled") return
c$$$      if (pltovlp.eq."enabled".and. mrfn.gt.1) then
ccc      call bcast(temp1(1,0),zero,iy*(jx+1)) ! YuP-101215: error?
      call bcast(temp1(0,0),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)
      
      write(*,*)'rdc_bplt(krf): mrfn =',mrfn,' krf=',krf
      
c$$   do 560 k=1,mrfn
      do 561 j=1,jx
         do 562 i=1,iy
            temp1(i,j)=rdcb(i,j,lr_,krf)
 562     continue
 561  continue
      CALL PGPAGE
      itype=7 ! means: plots are made for rdcb
      call pltcont(nrdcspecies(krf),1,
     +     'Contours of RdcB vs. v_parallel,v_perp',7)
      write(t_,552) lr_
 552  format(" Flux surface number",i3,";   all modes, krf=",i2)
      CALL PGMTXT('B',R410,R40,R40,t_)
      
c$$$  560    continue
      
c$$$  endif


      return
      end


