
      SUBROUTINE ADCALC(INUCZ,TE,DENSE,TA,DENSA)
C
C
C**********************************************
C      PERFORM MAIN CALCULATIONS
C**********************************************
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
C----sequre the input
C
      ZTE=TE/1000.
      ZNE=DENSE
      ZTA=TA/1000.
      ZNA=DENSA
C
C----ADCINIT initialize the code
C
      CALL ADCINIT
C
C-----ADCSET sets up the energy levels
C
      CALL ADCSET(INUCZ)
C
C-----ADCERC computes the electron rates
C
      CALL ADCERC(ZTE,ZNE)
C
C-----ADCXR  computes the charge exchange rates
C
      CALL ADCXR(ZNA,ZTA)
C
C-----ADCE computes coronal equilibrium rates and charge state fractions
C
      CALL ADCE(ZTE,ZNE,ZTA,ZNA)
C
      RETURN
      END



      SUBROUTINE ADCDO(INUCZ,TE,DENSE,TA,DENSA,
     &                 ZS,Z2S, EZRAD, EZLOSS, ZDISTR)
C
C*******************************************************
C    CALCULATES AND PUT OURPUT DATA INTO USER ARRAYS
C*******************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      DIMENSION ZDISTR(*)
C
C-----DO THE CALCULATIONS
C
      CALL ADCALC(INUCZ,TE,DENSE,TA,DENSA)
C
C-----ADCOUT FILLS UP THE USER OUTPUT DATA
C
      CALL ADCOUT(ZS,Z2S,EZRAD,EZLOSS,ZDISTR)
C
      RETURN
      END




      SUBROUTINE ADCDU(INUCZ,TE,DENSE,TA,DENSA)
C
C*******************************************
C  CALCULATES AND PRINT I/O DATA INTO FILE
C*******************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
C-----DO THE CALCULATIONS
C
      CALL ADCALC(INUCZ,TE,DENSE,TA,DENSA)
C
C-----ADCAUT prints detailed output
C
      CALL ADCAUT(TE,DENSE,TA,DENSA)
C
      RETURN
      END




      SUBROUTINE ADCDU1(INUCZ,TE,DENSE,TA,DENSA)
C
C*******************************************
C  CALCULATES AND PRINT I/O DATA INTO FILE
C*******************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
C-----DO THE CALCULATIONS
C
      CALL ADCALC(INUCZ,TE,DENSE,TA,DENSA)
C
C-----ADCAUT prints detailed output
C
      CALL ADCAUT1(TE,DENSE,TA,DENSA)
C
      RETURN
      END

