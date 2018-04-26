c
c
      subroutine dsk_gr
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     At the end of the run, this routine at the option of the user 
c     writes out to disk file 'idskf' the meshes on which the
c     distribution function is calculated, and 
c     (((f(i,j,k,l),i=1,iy),j=1,jx),l=1,lrors), for each species k.
c..................................................................

cBH000506:   De-activated by following goto.

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

CMPIINSERT_IF_RANK_NE_0_RETURN

      go to 999

      if (lrzmax.le.1) then
        if(idskf.eq. "disabled". or. n.ne.nstop+1)  return
      else
        if(idskf.eq. "disabled". or. n.ne.nstop)  return
      endif
      if(l_.ne.lrors)  return
      open(unit=41,file='graph',status='unknown')
ccc      close(unit=2) ! YuP: Why here?
      ilen=0
      write(41,1004)  iy,jx,lrors,lrzmax
      write(41,1005)  vnorm
      write(41,1005)  ((y(i,l),i=1,iy),l=1,lrors)
      write(41,1005)  (x(j),j=1,jx)
      write(41,1005)  (rz(l),l=1,lrzmax)
      do 1000 k=1,ngen
        write(41,1005)  (((f(i,j,k,l),i=1,iy),j=1,jx),l=1,lrors)
 1000 continue
 1003 format(a80)
 1004 format(5i16)
 1005 format(5e16.8)
      close(unit=41)
 999  return
      end
