c
c
      subroutine frwrong(kerr)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE     

! print error messages - on mpirank.eq.0 only
CMPIINSERT_IF_RANK_EQ_0
      if (kerr.eq.1 .or. kerr.eq.2) then
        WRITE(*,10010)
      elseif (kerr.eq.3) then
        WRITE(*,10030)
      elseif (kerr.eq.4) then
        WRITE(*,10040)
      elseif (kerr.eq.5) then
        WRITE(*,10050)
      else
        WRITE(*,10990)
      endif
CMPIINSERT_ENDIF_RANK
      
10010 format("kspeci not properly input for use with fr-module")
10030 format("subroutine freyasou - loop 110")
10040 format("subroutine freyasou - rcon iteration -meth2")
10050 format("subroutine freyasou - luf failure i.gt.iy")
10990 format( "fr module-unspecified")
      stop 'frwrong:' ! stop at all MPI cores
      end
