      SUBROUTINE ADCOUT(ZS,Z2S,EZRAD,EZLOSS,ZDISTR)
C
C**********************************************
C       PASS THE OUPUT DATA TO URER ARRAYS
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
      COMMON / ADPARM / LADEN, LADTIP, LTIPOK, LECI, LECIOK
C
      COMMON / ADDIEL / CDNN(100), CDNM(100), 
     &                  YHNNA, YHNNB, YHNNC, YHNMA, YHNMB, YHNMC,
     &                  LDRMLT
C
      COMMON / ADCEQ / CEAVGZ, CEAVGZ2, CERAD, CELOS,
     &                 CEFRAC(100),TOTREC(100),TOTION(100), 
     &                 TOTRAD(100)
C
      COMMON / ADCXXC / ANNEUT(5), VNEUT(5), 
     &                  NCXB, NCXOPT, NCXERR, IVUNIT
C
      DIMENSION ZDISTR(*)
C
      ZS=CEAVGZ
      Z2S=CEAVGZ2
      EZRAD=CERAD 
      EZLOSS=CELOS

C      print 99, zs,z2s,ezrad,ezloss
C  99  format(1x,'adcout: ',4(1p,E9.2,1x))

      DO 100 I=1,NSPC
      ZDISTR(I)=CEFRAC(I)
 100  CONTINUE
C
      RETURN
      END




      SUBROUTINE ADCOUTP(ITYP,ZTE,ZNE,
     &           ZS,Z2S,EZRAD,EZLOSS,ZDISTR)
C
C**********************************************
C  PRINT THE OUPUT DATA INTO FILE 'adcauxE.txt'
C**********************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      DIMENSION ZDISTR(*)
      INTEGER IOUT
C
      OPEN(UNIT=IOUT,FILE='adcauxE.txt')
C
      CALL ADCAUTH(IOUT)
      IF(ITYP.GT.0) CALL ADCAUTI(IOUT)
      CALL ADCAUTP(IOUT,ZTE,ZNE)

      IF(ITYP.GT.0) CALL ADCAUTR(IOUT,ZTE,ZNE)

      CALL ADCAUTS(IOUT,ZS,Z2S,EZRAD,EZLOSS,ZDISTR)
C
      CLOSE(UNIT=IOUT)
C

      RETURN
      END
