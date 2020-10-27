      SUBROUTINE ADCBREM(PTE, PNE)
C
C***************************************************************
C     INPUT : PTE, NUCZ, NSPC
C     OUTPUT : RADBRM
C
C     COMPUTE BREMSSTRAHLUNG RADIATION RADBRM(JQ) (JOULE/SEC)
C     FOR EACH IONIC SPECIES JQ = (IONIC CHARGE + 1).
C
C     REFERENCES:
C     1) POST ET AL. PPPL-1352 (CORRECTED)
C     2) KARZAS AND LATTER
C     3) TUCKER
C***************************************************************
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
C
      ZTE = PTE
C
C     USE FIXED FREE-FREE GAUNT FACTOR ZGFF = 1.2
C
      ZGFF = 1.2
C
      DO 100 JQ = 1 , NSPC
      ZCHG = DBLE(JQ - 1)
      RADBRM(JQ) = PNE * 4.85E-31 * ZCHG * ZCHG * DSQRT(ZTE) * ZGFF
  100 CONTINUE
C
      RETURN
      END
