      SUBROUTINE ADCECEX(PTE, PNE)
C
C**************************************************************
C     COMPUTE RADIATION RATES RADCLX(JQ) (JOULE/SEC) DUE TO
C     COLLISIONAL EXCITATION FOR EACH IONIC SPECIES JQ = (IONIC
C     CHARGE + 1).  ALL EXCITATIONS ARE ASSUMED TO RESULT IN
C     RADIATIVE DECAYS.
C
C     REFERENCES:
C     1) POST ET AL. PPPL-1352 (CORRECTED)
C**************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      COMMON / ADSDAT / NUCZ, NSPC, APN(100,10), AQN(100,10),
     &                  NVALNC(100), SIGMA(10,10), QN(100,10),
     &                  EN(100,10), EIN(100,5), ENM(100,5,10),
     &                  FNM(100,5,10), FNN(100), ENN(100)
C
      COMMON / ADELEC / RCLION(100), RRAREC(100), RDIREC(100),
     &                  CIZLOS(100), RADRRC(100), RADDRC(100),
     &                  RADCLX(100), RADBRM(100), TOTREC(100),
     &                  TOTRAD(100)
C
      COMMON / CEXDAT / ZDNN(28), ZGAMNN(28)
C
      EXTERNAL ADCEUND, ADCEXPE1
C
      DATA (ZDNN(J),J=1,28) / 0., 0.,
     & 1.2, .91, .89, .92, .86, .87, .85, 0.,
     & .72, .78, .76, .78, .75, .70, .68, .65, .65,
     & 8 * .70, 0. /
C
      DATA (ZGAMNN(J),J=1,28) / 0., 0.,
     & .54, .77, .58, .57, .41, .63, .66, 0.,
     & .97, .96, .88, .86, .87, .85, .83, .81,
     & 9 * .80, 0. /
C
      ESML=-45.d0
      ZJLKEV = 1.6021E-16
C
C
C
      ZTE = PTE
      ZNUC = DBLE(NUCZ)
C
C     ******** LOOP FOR ALL SPECIES EXCEPT FULLY IONIZED ********     
C
      DO 100 JQ = 1 , NUCZ
C
      ZCHG = DBLE(JQ - 1)
      IVALNC = NVALNC(JQ)
      IBOUND = NUCZ - JQ + 1
C
C     ***** N - N RATE *****
C
      ZZNN = 0.
      IF(FNN(JQ) .EQ. 0.) GO TO 30
C
      ZXNN = ENN(JQ) / ZTE
C
      IF(JQ .NE. 1) GO TO 20
      ZTEMP = .06 * (DSQRT(ZXNN) - 2.) / (1. + ZXNN)
      GO TO 22
C
   20 IF(IBOUND .LE. 28) GO TO 21
      ZTEMP = .7 * (1. - .8 / ZCHG)
      GO TO 22
C
   21 ZTEMP = ZDNN(IBOUND) * (1. - ZGAMNN(IBOUND) / ZCHG)
C
   22 ZGNN = ZTEMP + .276 * ADCEXPE1(ZXNN)
C
      ZZNN = FNN(JQ) * ADCEUND(-ZXNN , ESML) * ZGNN
C
C     ***** SUM N - M RATES FROM ALL OCCUPIED SHELLS *****
C
   30 ZNMSUM = 0.
C
C
      DO 31 JN = 1 , IVALNC
C
      ZJN = DBLE(JN)
      IIJM = IVALNC
      IF(JN .EQ. IVALNC .OR.
     &        AQN(JQ,IVALNC) .EQ. 0.) IIJM = IVALNC + 1
C
      DO 31 JM = IIJM , 10
      ZJM = DBLE(JM)
C
      ZXNM = ENM(JQ,JN,JM) / ZTE
      ZTEMP = ZJM * (ZJM - ZJN) * (1. + (1. - 2./ZNUC) * ZXNM) / 20.
      ZGNM = .19 * (1. + .9 * (1. + ZTEMP) * ADCEXPE1(ZXNM))
C
C     NOTE THAT FNM HAS APN AND AQN FACTORS ALREADY INCLUDED
C
      ZZNM = FNM(JQ,JN,JM) * ADCEUND(-ZXNM , ESML) * ZGNM
      ZNMSUM = ZNMSUM + ZZNM
   31 CONTINUE
C
C
C
      ZTEMP = (ZZNN + ZNMSUM) / DSQRT(ZTE)
      IF(ZTEMP .LT. 1.E-10) ZTEMP = 0.
      RADCLX(JQ) = PNE * ZJLKEV * 4.99E-10 * ZTEMP
C
  100 CONTINUE
C
C     ***************END OF JQ SPECIES LOOP ***************
C
      RADCLX(NSPC) = 0.
C
      RETURN
      END


