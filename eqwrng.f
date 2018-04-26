c
c
      subroutine eqwrng(kerr)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
CMPIINSERT_INCLUDE     

! print error messages - on mpirank.eq.0 only
CMPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'eqwrng(kerr) stopped with kerr=',kerr
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
      elseif (kerr.eq.11) then
        WRITE(*,10110)
      elseif (kerr.eq.12) then
        WRITE(*,10120)
      else
        WRITE(*,10990)
      endif
CMPIINSERT_ENDIF_RANK
      
10010 format("subroutine eqorbit - do loop 10 error")
10020 format("subroutine eqorbit - labels 15 to 20 - no convergence")
10030 format("subroutine eqorbit - zstep changed too often")
10040 format("subroutine eqfndpsi - no convergence")
10050 format("function eqfn - only psi model allowed is power model.")
10060 format("subroutine eqrfw")
10070 format("subroutine eqfpsi - model chosen is unavailable")
10080 format("subroutine eqfndpsi: erhocon(lr_) outside device")
10090 format("subroutine eqrhopsi - radmin too big")
10100 format("subroutine equilib-options are eqsource=eqdsk or topeol")
10110 format("subroutine equilib-nnv > nnr")
10120 format("sub equilib; nnr or nnz .gt. nnra or nnza")
10990 format("unspecified")

      stop 'eqwrng: error' ! stop at all MPI cores
      end
