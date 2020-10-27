      SUBROUTINE ADCERC(PTE,PNE)
C
C**********************************************
C     CALL ELECTRON RATE COEFFICIENT ROUTINES
C**********************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      CALL ADCECI(PTE,PNE)
      CALL ADCDREC(PTE,PNE)
      CALL ADCRREC(PTE,PNE)
      CALL ADCECEX(PTE,PNE)
      CALL ADCBREM(PTE,PNE)
C
      RETURN
      END
