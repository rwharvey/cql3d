c
c
      real*8 function siggy(ee)
      implicit integer (i-n), real*8 (a-h,o-z)
c---------------------------------------------------------------------
c    this routine is used in conjunction with array svtab defined
c    in subroutine sigsetup. it is the table look-up.
c---------------------------------------------------------------------
      include 'param.h'
      include 'comm.h'

      els=ee+em90
      inum=(els-elmin)/delegy+1.5
      inum=max0(inum,1) 
      inum=min0(inum,mtab)
      val=svtab(inum)
      siggy=val
      return
      end
