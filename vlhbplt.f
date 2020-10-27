      subroutine vlhbplt
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     plots rf cqlb coefficient as a contour plot.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*8 pltvlhb
      character*8 pltovlp
      REAL*4 :: R47=7.,R48=8.,R49=9.,R40=0.,R410=10.,R411=11.

      save pltvlhb,pltovlp
      data pltvlhb /'enabled'/
      data pltovlp /'enabled'/


CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return
      if (pltvlhb.ne."enabled") return
      if (pltovlp.eq."enabled".and. mrfn.gt.1) then
ccc        call bcast(temp1(1,0),zero,iy*(jx+1))  ! YuP-101215: error?
        call bcast(temp1(0,0),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)


        do 560 k=1,mrfn
          do 561 j=1,jx
            do 562 i=1,iy
              temp1(i,j)=cqlb(i,j,indxlr_,k)
 562        continue
 561      continue
          CALL PGPAGE
          itype=6 ! means: plots are made for vlhb
          call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp',
     +                 itype)
          write(t_,552) lr_
 552      format(" Flux surface number",i3,";   all modes")
          CALL PGMTXT('B',R410,R40,R40,t_)

 560    continue

      endif

      do 680 k=1,mrfn
c..................................................................
c     Compute vpar21/vte and vpar11/vte normalized velocities
c     Following vlh.f
c..................................................................
        if (vprprop .eq. "enabled") then
          xvpr=1.
        else
          xvpr=0.
        endif

c
c     determine the vparallel interval at R=R0
c
        vpmax=vparmax(k)*clight
        vpmin=vparmin(k)*clight
c
c     determine the vparallel range of nonzero D at outside
c     of flux surface.
c
        rovr0=(radmaj+xvpr*radmin*rovera(lr_))/radmaj
        vmin=vpmin*rovr0
        vmax=vpmax*rovr0

        vpar21dv=vmin/vth(1,lr_)
        vpar11dv=vmax/vth(1,lr_)


        do  j=1,jx
           do  i=1,iy
              temp1(i,j)=cqlb(i,j,indxlr_,k)
           enddo
        enddo

        CALL PGPAGE
        itype=6 ! means: plots are made for vlhb
        call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp',itype)
        write(t_,660) lr_,k
        CALL PGMTXT('B',R410,R40,R40,t_)
        write(t_,661) vpar21dv,vpar11dv
        CALL PGMTXT('B',R411,R40,R40,t_)

 680  continue
 660  format("Flux surface number",i3," mode=",i1)
 661  format("vpar21/vth=",1pe15.7,"   vpar11/vth=",1pe15.7)

      return
      end


