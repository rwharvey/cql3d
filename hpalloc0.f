c
c
      subroutine hpalloc0(iptr,ln)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'

c     NO_OP:  Not needed after 081216, after changing from
c             CRAY pointers to f90 pointers.

c..................................................................
c     Allocates memory AND zeroes out the array.
c..................................................................
c990131      integer err
c$$$CPTR>>>DELETE PTR-HPALLOC
c$$$      pointer(iptr,a(1))
c$$$
c$$$c_cray      call hpalloc(iptr,ln,err,1)
c$$$c_pc Unix libU77.a library.  malloc allocates bytes.
c$$$c_pc      iptr=malloc(8*ln)
c$$$c_pc      if (iptr.eq.0) stop 'Problem with malloc in hpalloc0'
c$$$      iptr=malloc(8*ln)
c$$$      if (iptr.eq.0) then
c$$$         write(*,*)'hpalloc0: ln=',ln
c$$$         stop 'Problem with malloc in hpalloc0'
c$$$      endif
c$$$
c$$$      zero=0.d0
c$$$      call bcast(a,zero,ln)
c$$$CPTR<<<END PTR-HPALLOC
      return
      end
