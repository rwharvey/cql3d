      SUBROUTINE ADCTIP(KTIPOK)
C
C*****************************************************************
C     SUBSTITUTE TABULATED IONIZATION POTENTIALS FOR THOSE
C     CALCULATED BY SHIELDING CONSTANTS IN SETENM.  NOTE THAT THE
C     'NUMBER OF EQUIVALENT ELECTRONS' IS NOT CORRECTLY RE-DEFINED
C     TO ACCOMPANY THIS PROCEDURE, BUT THE APPROXIMATION IS STILL
C     AN IMPROVEMENT FOR NEAR-NEUTRAL SPECIES.
C*****************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      COMMON / ADSDAT / NUCZ, NSPC, APN(100,10), AQN(100,10),
     &                  NVALNC(100), SIGMA(10,10), QN(100,10),
     &                  EN(100,10), EIN(100,5), ENM(100,5,10),
     &                  FNM(100,5,10), FNN(100), ENN(100)
C
      DIMENSION ZTIP(100)
C
      DO 1 JSPC = 1, 100
      ZTIP(JSPC) = 0.0
    1 CONTINUE
C
      IF(NUCZ .NE. 2) GO TO 3
C
C     HELIUM
C
      ZTIP(1) = 24.58E-3
      ZTIP(2) = 54.4E-3
C
      GO TO 100
    3 IF(NUCZ .NE. 3) GO TO 4
C
C     LITHIUM
C
      ZTIP(1) = 5.4E-3
      ZTIP(2) = 75.6E-3
      ZTIP(3) = 122.4E-3
C
      GO TO 100
C
    4 IF(NUCZ .NE. 4) GO TO 5
C
C     BERYLLIUM
C
      ZTIP(1) = 9.32E-3
      ZTIP(2) = 18.2E-3
      ZTIP(3) = 153.9E-3
      ZTIP(4) = 217.6E-3
C
      GO TO 100
C
    5 IF(NUCZ .NE. 5) GO TO 6
C
C     BORON
C
      ZTIP(1) = 8.296E-3
      ZTIP(2) = 25.15E-3
      ZTIP(3) = 37.9E-3
      ZTIP(4) = 259.3E-3
      ZTIP(5) = 340.0E-3
C
      GO TO 100
C
    6 IF(NUCZ .NE. 6) GO TO 8
C
C     CARBON
C
      ZTIP(1) = 11.3E-3
      ZTIP(2) = 24.4E-3
      ZTIP(3) = 47.9E-3
      ZTIP(4) = 64.5E-3
      ZTIP(5) = 0.392
      ZTIP(6) = 0.490
C
      GO TO 100
    8 IF(NUCZ .NE. 8) GO TO 10
C
C     OXYGEN
C
      ZTIP(1) = 13.6E-3
      ZTIP(2) = 35.1E-3
      ZTIP(3) = 54.9E-3
      ZTIP(4) = 77.4E-3
      ZTIP(5) = 0.114
      ZTIP(6) = 0.138
      ZTIP(7) = 0.739
      ZTIP(8) = 0.871
C
      GO TO 100
   10 IF(NUCZ .NE. 10) GO TO 13
C
C     NEON
C
      ZTIP(1) = 21.6E-3
      ZTIP(2) = 41.0E-3
      ZTIP(3) = 63.5E-3
      ZTIP(4) = 97.1E-3
      ZTIP(5) = 0.126
      ZTIP(6) = 0.158
      ZTIP(7) = 0.207
      ZTIP(8) = 0.239
      ZTIP(9) = 1.196
      ZTIP(10) = 1.362
C
      GO TO 100
   13 IF(NUCZ .NE. 13) GO TO 18
C
C     ALUMINUM
C
      ZTIP(1) = 5.99E-3
      ZTIP(2) = 18.8E-3
      ZTIP(3) = 28.4E-3
      ZTIP(4) = 0.120
      ZTIP(5) = 0.154
      ZTIP(6) = 0.190
      ZTIP(7) = 0.241
      ZTIP(8) = 0.285
      ZTIP(9) = 0.330
      ZTIP(10) = 0.399
      ZTIP(11) = 0.442
      ZTIP(12) = 2.086
      ZTIP(13) = 2.304
C
      GO TO 100
   18 IF(NUCZ .NE. 18) GO TO 22
C
C     ARGON
C
      ZTIP(1) = 15.8E-3
      ZTIP(2) = 27.6E-3
      ZTIP(3) = 40.7E-3
      ZTIP(4) = 59.8E-3
      ZTIP(5) = 75.0E-3
      ZTIP(6) = 91.0E-3
      ZTIP(7) = 0.124
      ZTIP(8) = 0.143
      ZTIP(9) = 0.422
      ZTIP(10) = 0.479
      ZTIP(11) = 0.539
      ZTIP(12) = 0.618
      ZTIP(13) = 0.686
      ZTIP(14) = 0.756
      ZTIP(15) = 0.855
      ZTIP(16) = 0.918
      ZTIP(17) = 4.121
      ZTIP(18) = 4.426
C
      GO TO 100
   22 IF(NUCZ .NE. 22) GO TO 26
C
C     TITANIUM
C
      ZTIP(1) = 6.82E-3
      ZTIP(2) = 13.6E-3
      ZTIP(3) = 27.5E-3
      ZTIP(4) = 43.3E-3
      ZTIP(5) = 99.2E-3
      ZTIP(6) = 0.119
      ZTIP(7) = 0.141
      ZTIP(8) = 0.169
      ZTIP(9) = 0.193
      ZTIP(10) = 0.216
      ZTIP(11) = 0.265
      ZTIP(12) = 0.292
      ZTIP(13) = 0.787
      ZTIP(14) = 0.861
      ZTIP(15) = 0.940
      ZTIP(16) = 1.042
      ZTIP(17) = 1.131
      ZTIP(18) = 1.220
      ZTIP(19) = 1.342
      ZTIP(20) = 1.425
      ZTIP(21) = 6.249
      ZTIP(22) = 6.626
C
      GO TO 100
   26 IF(NUCZ .NE. 26) GO TO 100
C
C     IRON
C
      ZTIP(1) = 7.87E-3
      ZTIP(2) = 16.2E-3
      ZTIP(3) = 30.7E-3
      ZTIP(4) = 54.8E-3
      ZTIP(5) = 75.0E-3
      ZTIP(6) = 99.0E-3
      ZTIP(7) = 0.125
      ZTIP(8) = 0.151
      ZTIP(9) = 0.235
      ZTIP(10) = 0.262
      ZTIP(11) = 0.290
      ZTIP(12) = 0.331
      ZTIP(13) = 0.361
      ZTIP(14) = 0.392
      ZTIP(15) = 0.457
      ZTIP(16) = 0.490
      ZTIP(17) = 1.265
      ZTIP(18) = 1.358
      ZTIP(19) = 1.456
      ZTIP(20) = 1.582
      ZTIP(21) = 1.689
      ZTIP(22) = 1.799
      ZTIP(23) = 1.950
      ZTIP(24) = 2.045
      ZTIP(25) = 8.828
      ZTIP(26) = 9.278
C
C     END OF CURRENT TABULATION; RETURN IF SELECTED
C     ELEMENT NOT PRESENT AND SET FLAG.
C
  100 KTIPOK = 0
      IF(ZTIP(1) .EQ. 0.0) RETURN
C
      KTIPOK = 1
C
      DO 150 JSPC = 1, NUCZ
      IVALNC = NVALNC(JSPC)
      JSPCP1 = JSPC + 1
C
      EIN(JSPC,IVALNC) = ZTIP(JSPC)
      EN(JSPCP1,IVALNC) = ZTIP(JSPC)
  150 CONTINUE
C
      RETURN
      END

