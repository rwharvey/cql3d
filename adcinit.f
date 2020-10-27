      SUBROUTINE ADCINIT
C
C**********************************************
C    Initializations for atomic data sets
C**********************************************
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
     &                  RADCLX(100), RADBRM(100)
C
      COMMON / ADNEUT / RCXREC(100), RADBCX(100)
C
      COMMON / ADCEQ / CEAVGZ, CEAVGZ2, CERAD, CELOS,
     &                 CEFRAC(100),TOTREC(100),TOTION(100), 
     &                 TOTRAD(100)
C
      COMMON / ADPARM / LADEN, LADTIP, LTIPOK, LECI, LECIOK
C
      COMMON / ADDIEL / CDNN(100), CDNM(100), 
     &                  YHNNA, YHNNB, YHNNC, YHNMA, YHNMB, YHNMC,
     &                  LDRMLT
C
      COMMON / ADCXXC / ANNEUT(5), VNEUT(5), 
     &                  NCXB, NCXOPT, NCXERR, IVUNIT
C
C---DIELECTRONIC MULTIPLIER (1=1.0, 2=HAHN 3=INPUT HAHN)
      LDRMLT = 2
C
C---HAHN CONSTANTS (NN A,B,C; NM A,B,C)
      IF(LDRMLT-2) 1,2,3
C
   1  YHNNA = 1.0
      YHNNB = 1.0
      YHNNC = 1.0
      YHNMA = 1.0
      YHNMB = 1.0
      YHNMC = 1.0
      GO TO 4
C
   2  YHNNA = 0.2
      YHNNB = 0.7
      YHNNC = 1.4
      YHNMA = 0.1
      YHNMB = 0.55
      YHNMC = 1.0
      GO TO 4
C
   3  YHNNA = 1.0
      YHNNB = 1.0
      YHNNC = 1.0
      YHNMA = 1.0
      YHNMB = 1.0
      YHNMC = 1.0
C
C----MAYER(0) OR MORE(1) ENERGIES; TAB IP(YES=1)
   4  LADEN = 1
      LADTIP = 1
C
C----IONIZATION RATES (XSNQ=1, BELFAST=2, YOUNGER=3)
      LECI = 2
C
C----X-S OPTION (0-NONE 1=OSAS,2=GJ,3=OSCT)
      NCXOPT=2
C
C----ENERGY UNITS(1=CM/S,2=KEV/AMU)
      IVUNIT=2
C
C----NUMBER OF NEUTRAL BEAM COMPONENTS
      NCXB=1
C
C----MISC
      NCXERR=0
C
      DO 100 I=1,100
      CDNN(I)=0.
      CDNM(I)=0.
100   CONTINUE
C
      DO 101 J=1,10
      NVALNC(I)=0
      FNN(I)=0.0 
      ENN(I)=0.0
      DO 102 I=1,100
        APN(I,J)=0.0
        AQN(I,J)=0.0
        SIGMA(I,J)=0.0 
        QN(I,J)=0.0
        EN(I,J)=0.0
 102  CONTINUE
 101  CONTINUE
C
      DO 103 I=1,100
      DO 104 K=1,5
      EIN(I,K)=0.0
      DO 105 J=1,10
      ENM(I,K,J)=0.0
      FNM(I,K,J)=0.0
 105  CONTINUE
 104  CONTINUE 
 103  CONTINUE
C
      DO 110 I=1,100
      RCLION(I)=0.0 
      RRAREC(I)=0.0
      RDIREC(I)=0.0
      CIZLOS(I)=0.0 
      RADRRC(I)=0.0
      RADDRC(I)=0.0
      RADCLX(I)=0.0 
      RADBRM(I)=0.0 
 110  CONTINUE
C
      DO 120 I=1,100
      RCXREC(I)=0.0
      RADBCX(I)=0.0
 120  CONTINUE
C
      CEAVGZ=0.0
      CEAVGZ2=0.0
      CERAD=0.0
      CELOS=0.0
      DO 130 I=1,100
      TOTION(I)=0.0
      CEFRAC(I)=0.0
      TOTREC(I)=0.0
      TOTRAD(I)=0.0
 130  CONTINUE
C
      DO 140 I=1,5
      ANNEUT(I)=0.0
      VNEUT(I)=0.0
 140  CONTINUE
C
      RETURN
      END
