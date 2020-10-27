      REAL*8 FUNCTION ADCFCHI(PCHI, PA, PB, PC, PD)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      Z1 = PA + PB * (1.0 + 1.0 / PCHI)
      Z2 = PC - (PA + 2.0*PB + PB/PCHI)/PCHI
      Z3 = PD / PCHI

      ZCHISQ = PCHI * PCHI
      ZCHICB = ZCHISQ * PCHI

      ZALPHA = (.001193 + .9764*PCHI +.6604*ZCHISQ +
     &          .02590*ZCHICB) / (1.0 + 1.488*PCHI +
     &          .2972*ZCHISQ + .004925*ZCHICB)
      ZBETA = (-.0005725 + .01345*PCHI +.8691*ZCHISQ +
     &         .03404*ZCHICB) / (1.0 + 2.197*PCHI +
     &         .2454*ZCHISQ + .002053*ZCHICB)

      ADCFCHI = (3.0E13 / PCHI) * (Z1 + Z2 * ZALPHA + Z3 * ZBETA)
      RETURN

      END



      REAL*8 FUNCTION ADCEXPE1(PX)
C
C****************************************************************
C     CALCULATE EXP(PX) * E1(PX), WHERE THE EXPONENTIAL INTEGRAL
C     E1(X) = INTEGRAL FROM X TO INFINITY OF (EXP(-T)/T)DT.
C
C     REFERENCE: 'HANDBOOK OF MATHEMATICAL FUNCTIONS', ABRAMOWITZ
C     AND STEGUN, NATIONAL BUREAU OF STANDARDS
C****************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      IF(PX.GE.1.) GO TO 2
      IF(PX.GT.1.d-4) GO TO 1
C
      ADCEXPE1 = DEXP(PX) * ( - DLOG(PX) - .57721566 + PX )
      RETURN
C
C
    1 ADCEXPE1 = DEXP(PX) * ( - DLOG(PX) - .57721566 + .99999193*PX
     &         - .24991055*PX**2 + .05519968*PX**3 - .00976004*PX**4
     &         + .00107857*PX**5 )
C
      RETURN
C
C
    2 ADCEXPE1 = (PX + 2.334733 + .250621/PX) /
     &         (PX**2 + 3.330657*PX + 1.681534)
C
      RETURN
      END





      REAL*8 FUNCTION ADCEUND(PX,PCHK)
C
C**********************************************************
C     PREVENT UNDERFLOWS IN RATE CALCULATIONS BY RETURNING
C     0. FOR EXP(PX) IF PX < PCHK.
C**********************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      IF( PX - PCHK ) 1 , 1 , 2
    1 ADCEUND = 0.
      RETURN
C
    2 ADCEUND = DEXP(PX)
      RETURN
C
      END