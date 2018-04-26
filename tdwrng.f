

c-----------------------------------------------------------------
      subroutine tdwrng(kerr)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c...  
c
cmnt  diagnostic error messages
c
c...  

      include 'param.h'
CMPIINSERT_INCLUDE     
      
! print error messages - on mpirank.eq.0 only
CMPIINSERT_IF_RANK_EQ_0
      if (kerr.eq.1) then
        WRITE(*,10010)
      elseif (kerr.eq.2) then
        WRITE(*,10020)
      elseif (kerr.eq.3) then
        WRITE(*,10030)
      elseif (kerr.eq.4) then
        WRITE(*,10040)
      elseif (kerr.eq.5) then
        WRITE(*,10050)
      elseif (kerr.eq.6) then
        WRITE(*,10060)
      elseif (kerr.eq.7) then
        WRITE(*,10070)
      elseif (kerr.eq.8) then
        WRITE(*,10080)
      elseif (kerr.eq.9) then
        WRITE(*,10090)
      elseif (kerr.eq.10) then
        WRITE(*,10100)
      else
        WRITE(*,10990)
      endif
CMPIINSERT_ENDIF_RANK
      
      stop ! stop at all MPI cores
      
c      if (kerr.eq.0) return
10010 format("machine = toroidal only  for 3-d calc. in cqlinput")
10020 format("wdmodel =ech1 only for 3-d calc. in cqlinput")
10030 format("rya(ll).lt.rya(ll-1) for an ll; rzset error (cqlinput)")
10040 format ("subroutine tdsxr - thpol out of bounds")
10050 format("subroutine tdsxr - field line does not intersect viewing",
     1  " angle")
10060 format("subroutine tdtodskr - hpdealc returns error message")
10070 format("subroutine ???? - iostatus returns error message")
10080 format("transport model- meshy=fixed_mu.  Insufficient resolution"
     1  "in the trapped region for some l. Set tfac larger")
10090 format("transport model - transport mesh falls EXACTLTY on p/t",/,
     1  "boundary. Perturb tfac")
10100 format("STOP in tdtry: iactst=abort forces stop. error>1.e-8")
10990 format("unspecified")
      return
      end
