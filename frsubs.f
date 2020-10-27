c***** SUBROUTINE EBALAF (A,N,IA,D,K,L)
C   IMSL ROUTINE NAME   - EBALAF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBALAF (A,N,IA,D,K,L)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,K,L
      REAL*8             A(IA,*),D(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L1,K1,K1P1,K11,JJ,J,I,LL,NOCONV
      REAL*8             R,C,F,G,B,S,B2,ONE,ZERO,P95
      DATA               B/16.0/,B2/256.0/
      DATA               ZERO/0.0/,ONE/1.0/,P95/.95/
C                                  REDUCE NORM A BY DIAGONAL SIMILARITY
C                                  TRANSFORMATION STORED IN D
C                                  FIRST EXECUTABLE STATEMENT
      L1 = 1
      K1 = N
C                                  SEARCH FOR ROWS ISOLATING AN EIGEN-
C                                    VALUE AND PUSH THEM DOWN
    5 K1P1 = K1+1
      IF (K1.LT.1) GO TO 35
      K11=K1
      DO 30 JJ=1,K11
         J = K1P1-JJ
         R = ZERO
         DO 10 I=1,K1
            IF (I.EQ.J) GO TO 10
            R=R+ABS(A(J,I))
   10    CONTINUE
         IF (R.NE.ZERO) GO TO 30
         D(K1) = J
         IF (J.EQ.K1) GO TO 25
         DO 15 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,K1)
            A(I,K1) = F
   15    CONTINUE
         DO 20 I=L1,N
            F = A(J,I)
            A(J,I) = A(K1,I)
            A(K1,I) = F
   20    CONTINUE
   25    K1 = K1-1
         GO TO 5
   30 CONTINUE
C                                  SEARCH FOR COLUMNS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM LEFT
   35 IF (K1.LT.L1) GO TO 65
      LL = L1
      DO 60 J=LL,K1
         C = ZERO
         DO 40 I=L1,K1
            IF (I.EQ.J) GO TO 40
            C = C+ABS(A(I,J))
   40    CONTINUE
         IF (C.NE.ZERO) GO TO 60
         D(L1) = J
         IF (J.EQ.L1) GO TO 55
         DO 45 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,L1)
            A(I,L1) = F
   45    CONTINUE
         DO 50  I=L1,N
            F = A(J,I)
            A(J,I) = A(L1,I)
            A(L1,I) = F
   50    CONTINUE
   55    L1 = L1+1
         GO TO 35
   60 CONTINUE
C                                  NOW BALANCE THE SUBMATRIX IN ROWS
C                                    L1 THROUGH K1
   65 K = L1
      L = K1
      IF (K1.LT.L1) GO TO 75
      DO 70  I=L1,K1
         D(I) = ONE
   70 CONTINUE
   75 NOCONV = 0
      IF (K1.LT.L1) GO TO 120
      DO 115 I=L1,K1
         C = ZERO
         R = ZERO
         DO 80 J=L1,K1
            IF (J.EQ.I) GO TO 80
            C = C+ABS(A(J,I))
            R = R+ABS(A(I,J))
   80    CONTINUE
         G = R/B
         F = ONE
         S = C+R
   85    IF (C.GE.G) GO TO 90
         F = F * B
         C = C*B2
         GO TO 85
   90    G = R*B
   95    IF (C.LT.G) GO TO 100
         F = F/B
         C = C/B2
         GO TO 95
C                                  NOW BALANCE
  100    IF ((C+R)/F.GE.P95*S) GO TO 115
         G = ONE/F
         D(I) = D(I)*F
         NOCONV = 1
         DO 105 J=L1,N
            A(I,J) = A(I,J)*G
  105    CONTINUE
         DO 110 J=1,K1
            A(J,I) = A(J,I)*F
  110    CONTINUE
  115 CONTINUE
  120 IF (NOCONV.EQ.1) GO TO 75
      RETURN
      END

c***** SUBROUTINE EBBCKF (D,Z,K,L,MM,N,IZ)
C   IMSL ROUTINE NAME   - EBBCKF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - DOUBLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBBCKF (D,Z,K,L,MM,N,IZ)
      implicit integer (i-n), real*8 (a-h,o-z)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,MM,N,IZ
      REAL*8             Z(IZ,*),D(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KM1,II,JJ,LP1
      REAL*8             S
C                                  COLUMN SCALE Z BY APPROPRIATE D VALUE
C                                  FIRST EXECUTABLE STATEMENT
      IF (L.EQ.0) GO TO 15
      DO 10 I=K,L
         S = D(I)
         DO 5 J=1,MM
            Z(I,J) = Z(I,J)*S
    5    CONTINUE
   10 CONTINUE
C                                  INTERCHANGE ROWS IF PERMUTATIONS
C                                    OCCURRED IN EBALAF
   15 IF (K.EQ.1) GO TO 30
      KM1 = K-1
      DO 25 I=1,KM1
         II = K-I
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 25
         DO 20 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   20    CONTINUE
   25 CONTINUE
   30 IF (L.EQ.N) GO TO 45
      LP1 = L+1
      DO 40 II=LP1,N
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 40
         DO 35 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   35    CONTINUE
   40 CONTINUE
   45 RETURN
      END

c***** SUBROUTINE EHBCKF (Z,H,D,N,MM,IZH,K,L)
C   IMSL ROUTINE NAME   - EHBCKF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - DOUBLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHBCKF (Z,H,D,N,MM,IZH,K,L)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MM,IZH,K,L
      REAL*8             Z(IZH,*),H(IZH,*),D(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LM2,KI,LTEMP,M,MA,MP2,I,J
      REAL*8             T,G,TINV,ZERO,ONE
      DATA               ZERO,ONE/0.0,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      LM2=L-2
      IF(LM2.LT.K) GO TO 9005
      LTEMP=LM2+K
      DO 30 KI=K,LM2
         M=LTEMP-KI
         MA=M+1
         T=H(MA,M)
         IF(T.EQ.ZERO) GO TO 30
         T=T*D(MA)
         MP2=M+2
         IF(MP2.GT.L) GO TO 10
         DO 5 I=MP2,L
            D(I)=H(I,M)
    5    CONTINUE
   10    IF(MA.GT.L) GO TO 30
         TINV = ONE / T
         DO 25 J=1,MM
            G=ZERO
            DO 15 I=MA,L
               G=G+D(I)*Z(I,J)
   15       CONTINUE
            G = G*TINV
            DO 20 I=MA,L
               Z(I,J)=Z(I,J)+G*D(I)
   20       CONTINUE
   25    CONTINUE
   30 CONTINUE
 9005 RETURN
      END

c***** SUBROUTINE EHESSF (A,K,L,N,IA,D)
C   IMSL ROUTINE NAME   - EHESSF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHESSF (A,K,L,N,IA,D)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IA
      REAL*8             A(IA,N),D(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LA,KP1,M,I,MP,II,J,JJ
      REAL*8             F,G,H,SCALE,ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      LA = L - 1
      KP1 = K + 1
      IF (LA .LT. KP1) GO TO 50
      DO 45 M = KP1, LA
         H = ZERO
         D(M) = ZERO
         SCALE = ZERO
C                                  SCALE COLUMN
         DO 5 I = M, L
            SCALE = SCALE + ABS(A(I,M-1))
    5    CONTINUE
         IF (SCALE .EQ. ZERO ) GO TO 45
         MP = M + L
C                                  DO 10 I=L,M,-1
         DO 10 II = M, L
            I = MP - II
            D(I) = A(I,M-1) / SCALE
            H = H + D(I) * D(I)
   10    CONTINUE
         G = -SIGN(SQRT(H),D(M))
         H = H - D(M) * G
         D(M) = D(M) - G
C                                  FORM (I-(U*UT)/H) * A
         DO 25 J = M,N
            F = ZERO
C                                  DO 15 I=L,M,-1
            DO 15 II = M, L
               I = MP - II
               F = F + D(I) * A(I,J)
   15       CONTINUE
            F = F / H
            DO 20 I = M, L
               A(I,J) = A(I,J) - F * D(I)
   20       CONTINUE
   25    CONTINUE
C                                  FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 40 I = 1,L
            F = ZERO
C                                  DO 30 J=L,M,-1
            DO 30 JJ = M, L
               J = MP - JJ
               F = F + D(J) * A(I,J)
   30       CONTINUE
            F = F / H
            DO 35 J = M, L
               A(I,J) = A(I,J) - F * D(J)
   35       CONTINUE
   40    CONTINUE
         D(M) = SCALE * D(M)
         A(M,M-1) = SCALE * G
   45 CONTINUE
   50 RETURN
      END
c
c
c***** SUBROUTINE EQRH3F (H,N,IH,K,L,WR,WI,Z,IZ,IER)
C   IMSL ROUTINE NAME   - EQRH3F
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - DOUBLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRH3F (H,N,IH,K,L,WR,WI,Z,IZ,IER)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IH,K,L,IZ,IER
      REAL*8             H(IH,N),WR(N),WI(N),Z(IZ,N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEN,ITS,IENM2,NPL,LL,LB,NAML,MM,M,MP2,KA,NA,
     *                   J,JJ
      REAL*8             T3(2),RDELP,P4,P5,P7,ZERO,ONE,T,X,Y,W,S,ZZ,R,P,
     *                   Q,RNORM,RA,SA,VR,VI
      COMPLEX*16         Z3
      LOGICAL            NOTLAS
      EQUIVALENCE        (Z3,T3(1))
      DATA               RDELP/0.710543E-14/
      DATA               P4 /0.4375/,P5 /0.5/,P7 /0.75/,ZERO /0.0/,ONE
     *                   /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  STORE ROOTS ISOLATED BY EBALAF
      RNORM = 0.0
      KA = 1
      DO 10 I=1,N
         DO 5 J=KA,N
    5    RNORM = RNORM+ABS(H(I,J))
         KA = I
         IF (I.GE.K .AND. I.LE.L) GO TO 10
         WR(I) = H(I,I)
         WI(I) = ZERO
   10 CONTINUE
      IEN = L
      T = ZERO
C                                  SEARCH FOR NEXT EIGENVALUES
   15 IF (IEN.LT.K) GO TO 145
      ITS = 0
      NA = IEN-1
      IENM2 = NA-1
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                  ELEMENT
   20 NPL = IEN+K
      DO 25 LL=K,IEN
         LB = NPL-LL
         IF (LB.EQ.K) GO TO 30
         S = ABS(H(LB-1,LB-1))+ABS(H(LB,LB))
         IF (S.EQ.0.0) S = RNORM
         IF (ABS(H(LB,LB-1)).LE.RDELP*S) GO TO 30
   25 CONTINUE
C
   30 X = H(IEN,IEN)
      IF (LB.EQ.IEN) GO TO 110
      Y = H(NA,NA)
      W = H(IEN,NA)*H(NA,IEN)
      IF (LB.EQ.NA) GO TO 115
      IF (ITS.EQ.30) GO TO 250
C                                  FORM SHIFT
      IF (ITS.NE.10 .AND. ITS.NE.20) GO TO 40
      T = T+X
      DO 35 I=K,IEN
         H(I,I) = H(I,I)-X
   35 CONTINUE
      S = ABS(H(IEN,NA))+ABS(H(NA,IENM2))
      X = P7*S
      Y = X
      W = -P4*S*S
   40 ITS = ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                  SUB-DIAGONAL ELEMENTS
      NAML = IENM2+LB
      DO 45 MM=LB,IENM2
         M = NAML-MM
         ZZ = H(M,M)
         R = X-ZZ
         S = Y-ZZ
         P = (R*S-W)/H(M+1,M)+H(M,M+1)
         Q = H(M+1,M+1)-ZZ-R-S
         R = H(M+2,M+1)
         S = ABS(P)+ABS(Q)+ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF (M.EQ.LB) GO TO 50
         IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.RDELP*ABS(P)*(ABS(H(M-1,
     *   M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GO TO 50
   45 CONTINUE
   50 MP2 = M+2
      DO 55 I=MP2,IEN
         H(I,I-2) = ZERO
         IF (I.EQ.MP2) GO TO 55
         H(I,I-3) = ZERO
   55 CONTINUE
C                                  DOUBLE QR STEP INVOLVING ROWS
C                                  L TO EN AND COLUMNS M TO EN
      DO 105 KA=M,NA
         NOTLAS = KA.NE.NA
         IF (KA.EQ.M) GO TO 60
         P = H(KA,KA-1)
         Q = H(KA+1,KA-1)
         R = ZERO
         IF (NOTLAS) R = H(KA+2,KA-1)
         X = ABS(P)+ABS(Q)+ABS(R)
         IF (X.EQ.ZERO) GO TO 105
         P = P/X
         Q = Q/X
         R = R/X
   60    CONTINUE
         S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (KA.EQ.M) GO TO 65
         H(KA,KA-1) = -S*X
         GO TO 70
   65    IF (LB.NE.M) H(KA,KA-1) = -H(KA,KA-1)
   70    P = P+S
         X = P/S
         Y = Q/S
         ZZ = R/S
         Q = Q/P
         R = R/P
C                                  ROW MODIFICATION
         DO 80 J=KA,N
            P = H(KA,J)+Q*H(KA+1,J)
            IF (.NOT.NOTLAS) GO TO 75
            P = P+R*H(KA+2,J)
            H(KA+2,J) = H(KA+2,J)-P*ZZ
   75       H(KA+1,J) = H(KA+1,J)-P*Y
            H(KA,J) = H(KA,J)-P*X
   80    CONTINUE
         J = MIN0(IEN,KA+3)
C                                  COLUMN MODIFICATION
         DO 90 I=1,J
            P = X*H(I,KA)+Y*H(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 85
            P = P+ZZ*H(I,KA+2)
            H(I,KA+2) = H(I,KA+2)-P*R
   85       H(I,KA+1) = H(I,KA+1)-P*Q
            H(I,KA) = H(I,KA)-P
   90    CONTINUE
         IF (IZ.LT.N) GO TO 105
C                                  ACCUMULATE TRANSFORMATIONS
         DO 100 I=K,L
            P = X*Z(I,KA)+Y*Z(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 95
            P = P+ZZ*Z(I,KA+2)
            Z(I,KA+2) = Z(I,KA+2)-P*R
   95       Z(I,KA+1) = Z(I,KA+1)-P*Q
            Z(I,KA) = Z(I,KA)-P
  100    CONTINUE
  105 CONTINUE
      GO TO 20
C                                  ONE ROOT FOUND
  110 H(IEN,IEN) = X+T
      WR(IEN) = H(IEN,IEN)
      WI(IEN) = ZERO
      IEN = NA
      GO TO 15
C                                  TWO ROOTS FOUND
  115 P = (Y-X)*P5
      Q = P*P+W
      ZZ = SQRT(ABS(Q))
      H(IEN,IEN) = X+T
      X = H(IEN,IEN)
      H(NA,NA) = Y+T
      IF (Q.LT.ZERO) GO TO 135
C                                  REAL PAIR
      ZZ = P+SIGN(ZZ,P)
      WR(NA) = X+ZZ
      WR(IEN) = WR(NA)
      IF (ZZ.NE.ZERO) WR(IEN) = X-W/ZZ
      WI(NA) = ZERO
      WI(IEN) = ZERO
      X = H(IEN,NA)
      R = SQRT(X*X+ZZ*ZZ)
      P = X/R
      Q = ZZ/R
C                                  ROW MODIFICATION
      DO 120 J=NA,N
         ZZ = H(NA,J)
         H(NA,J) = Q*ZZ+P*H(IEN,J)
         H(IEN,J) = Q*H(IEN,J)-P*ZZ
  120 CONTINUE
C                                  COLUMN MODIFICATION
      DO 125 I=1,IEN
         ZZ = H(I,NA)
         H(I,NA) = Q*ZZ+P*H(I,IEN)
         H(I,IEN) = Q*H(I,IEN)-P*ZZ
  125 CONTINUE
      IF (IZ.LT.N) GO TO 140
C                                  ACCUMULATE TRANSFORMATIONS
      DO 130 I=K,L
         ZZ = Z(I,NA)
         Z(I,NA) = Q*ZZ+P*Z(I,IEN)
         Z(I,IEN) = Q*Z(I,IEN)-P*ZZ
  130 CONTINUE
      GO TO 140
C                                  COMPLEX PAIR
  135 WR(NA) = X+P
      WR(IEN) = X+P
      WI(NA) = ZZ
      WI(IEN) = -ZZ
  140 IEN = IENM2
      GO TO 15
C                                  ALL ROOTS FOUND, NOW
C                                  BACKSUBSTITUTE
  145 IF (IZ.LT.N) GO TO 9005
      IF (RNORM.EQ.ZERO) GO TO 9005
      DO 220 NN=1,N
         IEN = N+1-NN
         P = WR(IEN)
         Q = WI(IEN)
         NA = IEN-1
         IF (Q.GT.ZERO) GO TO 220
         IF (Q.LT.ZERO) GO TO 180
C                                  REAL VECTOR
         M = IEN
         H(IEN,IEN) = ONE
         IF (NA.EQ.0) GO TO 220
         DO 175 II=1,NA
            I = IEN-II
            W = H(I,I)-P
            R = H(I,IEN)
            IF (M.GT.NA) GO TO 155
            DO 150 J=M,NA
               R = R+H(I,J)*H(J,IEN)
  150       CONTINUE
  155       IF (WI(I).GE.ZERO) GO TO 160
            ZZ = W
            S = R
            GO TO 175
  160       M = I
            IF (WI(I).NE.ZERO) GO TO 165
            T = W
            IF (W.EQ.ZERO) T = RDELP*RNORM
            H(I,IEN) = -R/T
            GO TO 175
C                                  SOLVE REAL EQUATIONS
  165       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
            T = (X*S-ZZ*R)/Q
            H(I,IEN) = T
            IF (ABS(X).LE.ABS(ZZ)) GO TO 170
            H(I+1,IEN) = (-R-W*T)/X
            GO TO 175
  170       H(I+1,IEN) = (-S-Y*T)/ZZ
  175    CONTINUE
C                                  END REAL VECTOR
         GO TO 220
C                                  LAST VECTOR COMPONENT CHOSEN
C                                    IMAGINARY SO THAT EIGENVECTOR
C                                    MATRIX IS TRIANGULAR
  180    M = NA
C                                  COMPLEX VECTOR
         IF (ABS(H(IEN,NA)).LE.ABS(H(NA,IEN))) GO TO 185
         H(NA,NA) = Q/H(IEN,NA)
         H(NA,IEN) = -(H(IEN,IEN)-P)/H(IEN,NA)
         GO TO 190
  185    CONTINUE
         Z3 = CMPLX(ZERO,-H(NA,IEN))/CMPLX(H(NA,NA)-P,Q)
         H(NA,NA) = T3(1)
         H(NA,IEN) = T3(2)
  190    H(IEN,NA) = ZERO
         H(IEN,IEN) = ONE
         IENM2 = NA-1
         IF (IENM2.EQ.0) GO TO 220
         DO 215 II=1,IENM2
            I = NA-II
            W = H(I,I)-P
            RA = ZERO
            SA = H(I,IEN)
            DO 195 J=M,NA
               RA = RA+H(I,J)*H(J,NA)
               SA = SA+H(I,J)*H(J,IEN)
  195       CONTINUE
            IF (WI(I).GE.ZERO) GO TO 200
            ZZ = W
            R = RA
            S = SA
            GO TO 215
  200       M = I
            IF (WI(I).NE.ZERO) GO TO 205
            Z3 = CMPLX(-RA,-SA)/CMPLX(W,Q)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            GO TO 215
C                                  SOLVE COMPLEX EQUATIONS
  205       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
            VI = (WR(I)-P)*Q
            VI = VI+VI
            IF (VR.EQ.ZERO .AND. VI.EQ.ZERO) VR = RDELP*RNORM*(ABS(W)
     *      +ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
            Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/CMPLX(VR,VI)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            IF (ABS(X).LE.ABS(ZZ)+ABS(Q)) GO TO 210
            H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,IEN))/X
            H(I+1,IEN) = (-SA-W*H(I,IEN)-Q*H(I,NA))/X
            GO TO 215
  210       CONTINUE
            Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,IEN))/CMPLX(ZZ,Q)
            H(I+1,NA) = T3(1)
            H(I+1,IEN) = T3(2)
  215    CONTINUE
C                                  END COMPLEX VECTOR
  220 CONTINUE
C                                  END BACKSUBSTITUTION
C                                  VECTORS OF ISOLATED ROOTS
      DO 230 I=1,N
         IF (I.GE.K .AND. I.LE.L) GO TO 230
         DO 225 J=I,N
            Z(I,J) = H(I,J)
  225    CONTINUE
  230 CONTINUE
      IF (L.EQ.0) GO TO 9005
C                                  MULTIPLY BY TRANSFORMATION MATRIX
      DO 245 JJ=K,N
         J = N+K-JJ
         M = MIN0(J,L)
         DO 240 I=K,L
            ZZ = ZERO
            DO 235 KA=K,M
               ZZ = ZZ+Z(I,KA)*H(KA,J)
  235       CONTINUE
            Z(I,J) = ZZ
  240    CONTINUE
  245 CONTINUE
      GO TO 9005
C                                  NO CONVERGENCE AFTER 30 ITERATIONS
C                                  SET ERROR INDICATOR  TO THE INDEX
C                                  OF THE CURRENT EIGENVALUE
  250 IER = 128+IEN
      DO 255 I=1,IEN
         WR(I) = ZERO
         WI(I) = ZERO
  255 CONTINUE
      IF (IZ.LT.N) GO TO 9000
      DO 265 I=1,N
         DO 260 J=1,N
            Z(I,J) = ZERO
  260    CONTINUE
  265 CONTINUE
 9000 CONTINUE
      CALL UERTST (IER,6HEQRH3F)
 9005 RETURN
      END

c
c
c***** SUBROUTINE EIGRF  (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C   IMSL ROUTINE NAME   - EIGRF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - AUGUST 1, 1981
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A REAL GENERAL MATRIX IN FULL STORAGE MODE
C
C   USAGE               - CALL EIGRF (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - THE INPUT REAL GENERAL MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                N      - THE INPUT ORDER OF THE MATRIX A.
C                IA     - THE INPUT ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IJOB   - THE INPUT OPTION PARAMETER. WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                W      - THE OUTPUT COMPLEX VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                         NOTE - THE ROUTINE TREATS W AS A REAL VECTOR
C                           OF LENGTH 2*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE W(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. IZ MUST BE GREATER
C                           THAN OR EQUAL TO N IF IJOB IS NOT EQUAL TO
C                           ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             (2+N)N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRH3F FAILED
C                           TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                           J+1,J+2,...,N HAVE BEEN COMPUTED CORRECTLY.
C                           EIGENVALUES 1,...,J ARE SET TO ZERO.
C                           IF IJOB = 1 OR 2 EIGENVECTORS ARE SET TO
C                           ZERO. THE PERFORMANCE INDEX IS SET TO 1000.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - EBALAF,EBBCKF,EHBCKF,EHESSF,EQRH3F,UERTST1,
C                           UGETIO1
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGRF  (A,N,IA,IJOB,W,Z,IZ,WK,IER)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IJOB,IZ,IER
      REAL*8             A(IA,*),WK(N,*),W(*),Z(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,IZ2,K,L,I,N1,N2,II,JJ,NP1,IIZ,NPI,JW,J,
     *                   IS,IG,IGZ,LW,LLZ,KKZ,LZ,KZ
      REAL*8             ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     *                   ZERO,ONE,THOUS,AN,Z11
      DATA               RDELP/0.710543E-14/
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IZ2 = IZ+IZ
      IF (IJOB .GE. 0 .AND. IJOB .LE. 3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 10
    5 IF (IJOB .EQ. 0) GO TO 16
   10 IF (IZ .GE. N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   15 IF (IJOB .EQ. 3) GO TO 95
C                                  PACK A INTO AN N BY N ARRAY
   16 K = 1
      L = 1
      DO 20 J=1,N
         DO 20 I=1,N
            A(K,L) = A(I,J)
C                                  SAVE INPUT A IF IJOB = 2
            IF (IJOB .EQ. 2) WK(I,J)=A(I,J)
            K = K+1
            IF (K .GT. IA) K = 1
            IF (K .EQ. 1) L = L+1
   20 CONTINUE
      N1 = 1
      IF (IJOB .EQ. 2) N1 = N+1
      N2 = N1+1
      IF (IJOB .EQ. 0) N2 = 1
C                                  BALANCE THE INPUT A
      CALL EBALAF (A,N,N,WK(1,N1),K,L)
      IF (IJOB .EQ. 0 .AND. L .EQ. 0) GO TO 35
C                                  IF L = 0, A IS ALREADY IN HESSENBERG
C                                    FORM
      CALL EHESSF (A,K,L,N,N,WK(1,N2))
      IF (IJOB .EQ. 0) GO TO 35
C                                  SET Z IDENTITY MATRIX
      II = 1
      JJ = 1
      NP1 = N+1
      DO 30 I=1,N
         DO 25 J=1,N
            Z(II) = ZERO
            II = II+1
   25    CONTINUE
         Z(JJ) = ONE
         JJ = JJ+NP1
   30 CONTINUE
      CALL EHBCKF (Z,A,WK(1,N2),N,N,N,K,L)
      IIZ = N
   35 IF (IJOB .EQ. 0) IIZ = 1
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z11 = Z(1)
      CALL EQRH3F (A,N,N,K,L,W(1),W(N+1),Z,IIZ,JER)
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z(1) = Z11
      IF (JER .GT. 128 .OR. IJOB .EQ. 0) GO TO 40
      CALL EBBCKF (WK(1,N1),Z,K,L,N,N,N)
C                                  CONVERT W (EIGENVALUES) TO COMPLEX
C                                    FORMAT
   40 DO 45 I=1,N
         NPI = N+I
         WK(I,N1) = W(NPI)
   45 CONTINUE
      JW = N+N
      J = N
      DO 50 I=1,N
         W(JW-1) = W(J)
         W(JW) = WK(J,N1)
         JW = JW-2
         J = J-1
   50 CONTINUE
      IF (IJOB .EQ. 0) GO TO 9000
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      J = N
   60 IF (J .LT. 1) GO TO 85
      IF (W(J+J) .EQ. ZERO) GO TO 75
C                                  MOVE PAIR OF COMPLEX CONJUGATE
C                                    EIGENVECTORS
      IS = IZ2*(J-1)+1
      IG = N*(J-2)+1
      IGZ = IG+N
C                                  MOVE COMPLEX CONJUGATE EIGENVECTOR
      DO 65 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IGZ)
         IS = IS+2
         IG = IG+1
         IGZ = IGZ+1
   65 CONTINUE
C                                  MOVE COMPLEX EIGENVECTOR
      IS = IZ2*(J-2)+1
      IG = IS+IZ2
      DO 70 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IG+1)
         IS = IS+2
         IG = IG+2
   70 CONTINUE
      J = J-2
      GO TO 60
C                                  MOVE REAL EIGENVECTOR
   75 IS = IZ2*(J-1)+N+N
      IG = N*J
      DO 80 I=1,N
         Z(IS-1) = Z(IG)
         Z(IS) = ZERO
         IS = IS-2
         IG = IG-1
   80 CONTINUE
      J = J-1
      GO TO 60
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
C                                    NEXT, MOVE ORIGINAL MATRIX BACK
C                                    TO A
   85 IF (IJOB .LE. 1) GO TO 9000
      DO 90 I=1,N
         DO 90 J=1,N
            A(I,J) = WK(I,J)
   90 CONTINUE
      WK(1,1) = THOUS
      IF (JER .NE. 0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
   95 ANORM = ZERO
      DO 105 J=1,N
         ASUM = ZERO
         DO 100 I=1,N
            ASUM = ASUM+ABS(A(I,J))
  100    CONTINUE
         ANORM = MAX (ANORM,ASUM)
  105 CONTINUE
      IF (ANORM .EQ. ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      LLZ = 0
      KKZ = 0
      DO 120 J=1,N
         S = ZERO
         SUMZ = ZERO
         LZ = LLZ+1
         KZ = KKZ+1
         LW = J+J-1
         DO 115 L=1,N
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            KZ = KKZ+1
            SUMR = -W(LW)*Z(LZ)+W(LW+1)*Z(LZ+1)
            SUMI = -W(LW)*Z(LZ+1)-W(LW+1)*Z(LZ)
            DO 110 K=1,N
               SUMR =SUMR+A(L,K)*Z(KZ)
               SUMI = SUMI+A(L,K)*Z(KZ+1)
               KZ = KZ+2
  110       CONTINUE
            S = S+CABS(CMPLX(SUMR,SUMI))
            LZ = LZ+2
  115    CONTINUE
         PI = MAX (PI,S/SUMZ)
         KKZ = KKZ+IZ2
         LLZ = LLZ+IZ2
  120 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1,1) = PI
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST1 (IER,6HEIGRF )
      IF (JER .EQ. 0) GO TO 9005
      IER = JER
      CALL UERTST1 (IER,6HEIGRF )
 9005 RETURN
      END

      subroutine leqt1fl(a, m, n, ia, b, idgt, wkarea, ier)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c-leqt1f--------s/d-----library 2--------------------------------------
c
c   function            - linear equation solution - full storage
c                           mode - space economizer solution.
c   usage               - call leqt1f (a,m,n,ia,b,idgt,wkarea,ier)
c   parameters   a      - input matrix of dimension n by n containing
c                           the coefficient matrix of the equation
c                           ax = b.
c                         on output, a is replaced by the lu
c                           decomposition of a rowwise permutation of
c                           a.
c                m      - number of right-hand sides.(input)
c                n      - order of a and number of rows in b.(input)
c                ia     - number of rows in the dimension statement
c                           for a and b in the calling program. (input)
c                b      - input matrix of dimension n by m containing
c                           right-hand sides of the equation ax = b.
c                         on output, the n by m solution x replaces b.
c                idgt   - input option.
c                         if idgt is greater than 0, the elements of
c                           a and b are assumed to be correct to idgt
c                           decimal digits and the routine performs
c                           an accuracy test.
c                         if idgt equals zero, the accuracy test is
c                           bypassed.
c                wkarea - work area of dimension greater than or equal
c                           to n.
c                ier    - error parameter
c                         terminal error = 128+n.
c                           n = 1 indicates that a is algorithmically
c                             singular. (see the chapter l prelude).
c                         warning error = 32+n.
c                           n = 2 indicates that the accuracy test
c                                 failed.
c                                 the computed solution may be in error
c                                 by more than can be accounted for by
c                                 the uncertainty of the data.
c                                 this warning can be produced only if
c                                 idgt is greater than 0 on input.
c                                 see chapter l prelude for further
c                                 discussion.
c   precision           - single/double
c   req'd IMSL routines - ludatf1,luelmf1,uertst1
c   language            - fortran
c ----------------------------------------------------------------------
c   latest revision     - march 22,1974
c                                  dec
c
      dimension          a(ia,*),b(ia,*),wkarea(*)
      dimension ipvt(n) ! YuP: added  ipvt - permutation indices.
****  real*8   a,b,wkarea,d1,d2,wa
c
      ier = 0
c                                  decompose a
c-YuP call ludatf1 (a,a,n,ia,idgt,d1,d2,wkarea,wkarea,wa,ier)
      call ludatf1 (a,a,n,ia,idgt,d1,d2,ipvt,wkarea,wa,ier) ! YuP: ipvt - out
      if (ier .gt. 128)  go to 9000
c                                  call routine luelmf1 (forward and
c                                  backward substitutions)
      do 10 j=1,m
c-YuP    call luelmf1 (a,b(1,j),wkarea,n,ia,b(1,j))
         call luelmf1 (a,b(1,j),ipvt,n,ia,b(1,j)) ! YuP: ipvt - input
   10 continue
      if (ier .eq. 0)  go to 9005
 9000 continue
      call uertst1 (ier,'leqt1fl')
 9005 return
c
      end


      subroutine ludatf1 (a, lu, n, ia, idgt, d1, d2, ipvt,
     .                    equil, wa, ier)
      implicit integer (i-n), real*8 (a-h,o-z)
c
c-ludatf1--------s/d-----library 2--------------------------------------
c
c   function            - l-u decomposition by the crout algorithm
c                           with optional accuracy test.
c   usage               - call ludatf1 (a,lu,n,ia,idgt,d1,d2,ipvt,
c                                       equil,wa,ier)
c   parameters   a      - input matrix of dimension n by n containing
c                           the matrix to be decomposed
c                lu     - real output matrix of dimension n by n
c                           containing the l-u decomposition of a
c                           rowwise permutation of the input matrix.
c                           for a description of the format of lu, see
c                           example.
c                n      - input scalar containing the order of the
c                           matrix a.
c                ia     - input scalar containing the row dimension of
c                           matrices a and lu in the calling program.
c                idgt   - input option.
c                           if idgt is greater than zero, the non-zero
c                           elements of a are assumed to be correct to
c                           idgt decimal places.  ludatf1 performs an
c                           accuracy test to determine if the computed
c                           decomposition is the exact decomposition
c                           of a matrix which differs from the given on
c                           by less than its uncertainty.
c                         if idgt is equal to zero, the accuracy test i
c                           bypassed.
c                d1     - output scalar containing one of the two
c                           components of the determinant. see
c                           description of parameter d2, below.
c                d2     - output scalar containing one of the
c                           two components of the determinant. the
c                           determinant may be evaluated as (d1)(2**d2)
c                ipvt   - output vector of length n containing the
c                           permutation indices. see document
c                           (algorithm).
c                equil  - output vector of length n containing
c                           reciprocals of the absolute values of
c                           the largest (in absolute value) element
c                           in each row.
c                wa     - accuracy test parameter, output only if
c                           idgt is greater than zero.
c                           see element documentation for details.
c                ier    - error parameter
c                         terminal error = 128+n
c                           n = 1 indicates that matrix a is
c                                 algorithmically singular. (see the
c                                 chapter l prelude).
c                         warning error = 32+n
c                           n = 2 indicates that the accuracy test
c                                 failed.
c                                 the computed solution may be in error
c                                 by more than can be accounted for by
c                                 the uncertainty of the data.
c                                 this warning can be produced only if
c                                 idgt is greater than 0 on input.
c                                 see chapter l prelude for further
c                                 discussion.
c   precision           - single/double
c   req'd IMSL routines - uertst1
c   language            - fortran
c ----------------------------------------------------------------------
c
c   latest revision     - march 22,1974 (dec)
c
      dimension a(ia,*), lu(ia,*), ipvt(*), equil(*)
      real*8    lu
      data      zero, one, four, sixtn, sixth
     .         /0.0 , 1.0, 4.0 , 16.0 , 0.0625/
c
      ier  = 0
      rn   = n
      wrel = zero
      d1   = one
      d2   = zero
      biga = zero
      do 10 i=1,n
         bigg = zero
         do 5 j=1,n
            p       = a(i,j)
            lu(i,j) = p
****        p       = DABS (p)
            p       =  ABS (p)
            if (p .gt. bigg)  bigg = p
    5    continue
         if (bigg .gt.  biga) biga = bigg
         if (bigg .eq.  zero)  go to 110
         equil(i) = one/bigg
   10 continue
      do 105 j=1,n
         jm1 = j-1
         if (jm1 .lt. 1)  go to 40
c
c compute u(i,j), i = 1,...,j-1
c
         do 35 i=1,jm1
            sum = lu(i,j)
            im1 = i-1
            if (idgt .eq. 0)  go to 25
c
c with accuracy test
c
            ai = ABS (sum)
            wi = zero
            if (im1 .lt. 1)  go to 20
            do 15 k=1,im1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi  = wi + ABS (t)
   15       continue
            lu(i,j) = sum
   20       wi = wi + ABS (sum)
            if (ai .eq. zero) ai = biga
            test = wi/ai
            if (test .gt. wrel) wrel = test
            go to 35
c
c without accuracy
c
   25       if (im1 .lt. 1)  go to 35
            do 30 k=1,im1
               sum = sum-lu(i,k)*lu(k,j)
   30       continue
            lu(i,j) = sum
   35    continue
   40    p = zero
c
c compute u(j,j) and l(i,j), i = j+1,...
c
         do 70 i=j,n
            sum = lu(i,j)
            if (idgt .eq. 0)  go to 55
c
c with accuracy test
c
            ai = ABS (sum)
            wi = zero
            if (jm1 .lt. 1)  go to 50
            do 45 k=1,jm1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi = wi + ABS (t)
   45       continue
            lu(i,j) = sum
   50       wi = wi + ABS (sum)
            if (ai .eq. zero) ai = biga
            test = wi/ai
            if (test .gt. wrel) wrel = test
            go to 65
c
c without accuracy test
c
   55       if (jm1 .lt. 1)  go to 65
            do 60 k=1,jm1
               sum = sum-lu(i,k)*lu(k,j)
   60       continue
            lu(i,j) = sum
   65       q = equil(i) * ABS (sum)
            if (p .ge. q)  go to 70
            p = q
            imax = i
   70    continue
c
c test for algorithmic singularity
c
         q = rn+p
         if (q .eq. rn)  go to 110
         if (j .eq. imax)  go to 80
c
c interchange rows j and imax
c
         d1 = -d1
         do 75 k=1,n
            p = lu(imax,k)
            lu(imax,k) = lu(j,k)
            lu(j,k) = p
   75    continue
         equil(imax) = equil(j)
   80    ipvt(j) = imax
         d1 = d1*lu(j,j)
   85    if (ABS (d1) .le. one)  go to 90
         d1 = d1*sixth
         d2 = d2+four
         go to 85
   90    if (ABS (d1) .ge. sixth)  go to 95
         d1 = d1*sixtn
         d2 = d2-four
         go to 90
   95    continue
         jp1 = j+1
         if (jp1 .gt. n)  go to 105
c
c divide by pivot element u(j,j)
c
         p = lu(j,j)
         do 100 i=jp1,n
            lu(i,j) = lu(i,j)/p
  100    continue
  105 continue
c
c perform accuracy test
c
      if (idgt .eq. 0)  go to 9005
      p = 3*n+3
      wa = p*wrel
****  q = wa+10.0d0**(-idgt)
      q = wa+10.0**(-idgt)
      if (q .ne. wa)  go to 9005
      ier = 34
      go to 9000
c
c algorithmic singularity
c
  110 ier = 129
      d1 = zero
      d2 = zero
c
c print error
c
 9000 call uertst1 (ier, 'ludatf1')
 9005 return
c
      end

      subroutine luelmf1 (a, b, ipvt, n, ia, x)
      implicit integer (i-n), real*8 (a-h,o-z)
c
c
c-luelmf1--------s/d-----library 2--------------------------------------
c
c   function            - elimination part of solution of ax = b -
c                           full storage mode
c   usage               - call luelmf1 (a,b,ipvt,n,ia,x)
c   parameters   a      - the result, lu, computed in the subroutine
c                           'ludatf1', where l is a lower triangular
c                           matrix with ones on the main diagonal. u is
c                           upper triangular. l and u are stored as a
c                           single matrix a, and the unit diagonal of
c                           l is not stored
c                b      - b is a vector of length n on the right hand
c                           side of the equation ax = b
c                ipvt   - the permutation matrix returned from the
c                           subroutine 'ludatf1', stored as an n length
c                           vector
c                n      - order of a and number of rows in b
c                ia     - number of rows in the dimension statement
c                           for a in the calling program.
c                x      - the result x
c   precision           - single/double
c   language            - fortran
c ----------------------------------------------------------------------
c   latest revision     - april 11,1975
c                                  dec
c
      dimension          a(ia,*),b(*),ipvt(*),x(*)
****  real*8   a,b,x,sum
c                                  solve ly = b for y
      do 5 i=1,n
    5 x(i) = b(i)
      iw = 0
      do 20 i=1,n
         ip = ipvt(i)
         sum = x(ip)
         x(ip) = x(i)
         if (iw .eq. 0)  go to 15
         im1 = i-1
         do 10 j=iw,im1
            sum = sum-a(i,j)*x(j)
   10    continue
         go to 20
   15    if (sum .ne. 0.0) iw = i
   20 x(i) = sum
c                                  solve ux = y for x
      do 30 ib=1,n
         i = n+1-ib
         ip1 = i+1
         sum = x(i)
         if (ip1 .gt. n)  go to 30
         do 25 j=ip1,n
            sum = sum-a(i,j)*x(j)
   25   continue
   30 x(i) = sum/a(i,i)
      return
c
      end


      subroutine uertst1 (ier, obsolete)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c --- changed name from UERTST to UERTST1 so that only IMSL routines
c --- with sources in ONETWO will call UERTST1. Other IMSL routines used
c --- in ONETWO which are linked from the IMSL library will thus not
c --- be involved with this subroutine.
c
c --- uertst1 ---------------- library 2 -------------------------------
c
c   function            - error message generation
c   usage               - call uertst1 (ier, obsolete)
c   parameters   ier    - error parameter. type + n  where
c                           type = 128 implies terminal error
c                                   64 implies warning with fix
c                                   32 implies warning
c                              n = error code relevant to calling routine
c              obsolete - input scalar (double precision on dec)
c                         containing the name of the calling routine
c                         as a 6-character literal string. --- OBSOLETE
c   language            - Fortran
c                         DEC
c
c ----------------------------------------------------------------------
c
c   latest revision     - october 1, 1975
c
c990131      include 'param.i'
c990131      include 'imsl.i'
c990131      include 'io.i'
c
      integer        warn, warf, term
      dimension      ibit(4)
      character*(*)  obsolete
      equivalence   (ibit(1), warn), (ibit(2), warf), (ibit(3), term)
      data           ibit / 32, 64, 128, 0 /
c
      ier2 = ier
      if (ier2 .ge. warn)  go to 5
c
c     non-defined
c
      ier1 = 4
      go to 20
    5 if (ier2 .lt. term)  go to 10
c
c     terminal
c
      ier1 = 3
      go to 20
   10 if (ier2 .lt. warf)  go to 15
c
c     warning (with fix)
c
      ier1 = 2
      go to 20
c
c     warning
c
   15 ier1 = 1
c
c     extract 'n'
c
   20 ier2 = ier2 - ibit(ier1)
c
c     print error message
c     disable undefined error output (it's not an error)
c
      if (ier1 .eq. 4)  return
c
c      write  (nout, 26)  imslmd, ier
c      write  (nqik, 26)  imslmd, ier
c      write  (ncrt, 26)  imslmd, ier
      write  (*, 26)  'imslmd', ier
      write  (*, 26)  'imslmd', ier
      write  (*, 26)  'imslmd', ier
   26 format (/ ' IMSL (uertst1) ERROR message from ', a,
     .          ':  ier (error code) =', i5)
      return
c
      end
c
c
      subroutine icsevu1 (x, y, nx, c, ic, u, s, m, ier)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c --- modified version of IMSL subroutine named ICSEVU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - evaluation of a cubic spline
c
c   usage               - call icsevu1 (x, y, nx, c, ic, u, s, m, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx (input). x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points
c                           (input).
c                nx     - number of elements in x and y (input).
c                           nx must be .ge. 2.
c                c      - spline coefficients (input). c is an nx-1 by
c                           3 matrix.
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement
c                           in the calling program (input).
c                           ic must be .ge. nx-1
c                u      - vector of length m containing the abscissae
c                           of the m points at which the cubic spline
c                           is to be evaluated (input).
c                s      - vector of length m (output).
c                           the value of the spline approximation at
c                           u(i) is
c                           s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c                           where x(j) .le. u(i) .lt. x(j+1) and
c                           d = u(i)-x(j).
c                m      - number of elements in u and s (input).
c                ier    - error parameter (output).
c                         warning error
c                           ier = 33, u(i) is less than x(1).
c                           ier = 34, u(i) is greater than x(nx).
c
c                           ********************************************
c                           output of warning errors has been suppressed
c                           ********************************************
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio1
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   remarks  1.  the routine assumes that the abscissae of the nx
c                data points are ordered such that x(i) is less than
c                x(i+1) for i = 1,...,nx-1. no check of this condition
c                is made in the routine. unordered abscissae will cause
c                the algorithm to produce incorrect results.
c            2.  the routine generates two warning errors. one error
c                occurs if u(i) is less than x(1), for some i in the
c                the interval (1,m) inclusively. the other error occurs
c                if u(i) is greater than x(nx), for some i in the
c                interval (1,m) inclusively.
c            3.  the ordinate y(nx) is not used by the routine. for
c                u(k) .gt. x(nx-1), the value of the spline, s(k), is
c                given by
c                 s(k) = ((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)
c                where d = u(k)-x(nx-1).
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,m,ier
      real*8             x(nx),y(nx),c(ic,3),u(m),s(m)
c
c     specifications for local variables
c
      integer            i,jer,ker,nxm1,k
      real*8             d,dd,zero
      data               i/1/, zero/0.0/
c
c     first executable statement
c
      jer = 0
      ker = 0
      if (m .le. 0)  go to 9005
      nxm1 = nx-1
      if (i .gt. nxm1)  i = 1
c
c     evaluate spline at m points
c
      do 40 k=1,m
c
c        find the proper interval
c
         d = u(k)-x(i)
         if (d) 5, 25, 15
    5    if (i .eq. 1)  go to 30
         i = i-1
         d = u(k)-x(i)
         if (d) 5, 25, 20
   10    i = i+1
         d = dd
   15    if (i .ge. nx)  go to 35
         dd = u(k)-x(i+1)
         if (dd .ge. zero)  go to 10
         if ( d .eq. zero)  go to 25
c
c        perform evaluation
c
   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
         go to 40
   25    s(k) = y(i)
         go to 40
c
c        u(k) < x(1)
c
   30    jer = 33
         go to 20
c
c        u(k) > x(nx)
c
   35    if (dd .gt. zero)  ker = 34
         d = u(k) - x(nxm1)
         i = nxm1
         go to 20
c
   40 continue
c
      ier = MAX0 (jer, ker)
c
****  if (jer .gt. 0)  call uertst1 (jer, 'icsevu1')
****  if (ker .gt. 0)  call uertst1 (ker, 'icsevu1')
c
 9005 return
c
      end
c
c
      subroutine icsicu1 (x, y, nx, bpar, c, ic, ier)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c --- modified version of IMSL subroutine named ICSICU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - interpolatory approximation by cubic splines
c                           with arbitrary second derivative end
c                           conditions.
c
c   usage               - call icsicu1 (x, y, nx, bpar, c, ic, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx. (input) x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points.
c                           (input)
c                nx     - number of elements in x and y. (input) nx
c                           must be .ge. 2.
c                bpar   - vector of length 4 containing the end
c                           condition parameters. (input)
c                           2.0*spp(1)+bpar(1)*spp(2) = bpar(2),
c                           bpar(3)*spp(nx-1)+2.0*spp(nx) = bpar(4),
c                           where spp(i) = second derivative of the
c                           cubic spline function s evaluated at x(i).
c                c      - spline coefficients. (output) c is an nx-1 by
c                           3 matrix. the value of the spline
c                           approximation at t is
c                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c                           where x(i) .le. t .lt. x(i+1) and
c                           d = t-x(i).
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement in
c                           the calling program. (input)
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129, ic is less than nx-1
c                           ier = 130, nx is less than 2.
c                           ier = 131, input abscissa are not ordered
c                             so that x(1) .lt. x(2) ... .lt. x(nx).
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio1
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,ier
      real*8             x(nx),y(nx),bpar(4),c(ic,3)
c
c     specifications for local variables
c
      integer            i,j,nxm1
      real*8             dx,dxj,dxjp1,dxp,dyj,dyjp1,half,one,pj,
     .                   six,sixi,two,yppa,yppb,zero
      equivalence        (dxj,yppb),(pj,sixi),(dxjp1,yppa)
      data               zero/0.0/,half/0.5/,one/1.0/,
     .                   two/2.0/,six/6.0/
c
      ier = 0
c
c     check error conditions
c
      nxm1 = nx-1
      if (ic .lt. nxm1)  go to 30
      if (nx .lt. 2   )  go to 35
      if (nx .eq. 2   )  go to 10
c
c     compute coefficients and right hand side of the tridiagonal
c     system defining the second derivatives of the spline interpolant for (x,y)
c
c     c(j,1) = lambda(j)
c     c(j,2) = mu(j)
c     c(j,3) = d(j)
c
      dxj = x(2)-x(1)
      if (dxj .le. zero)  go to 40
      dyj = y(2)-y(1)
      do 5 j=2,nxm1
         dxjp1 = x(j+1)-x(j)
         if (dxjp1 .le. zero)  go to 40
         dyjp1 = y(j+1)-y(j)
         dxp = dxj+dxjp1
         c(j,1) = dxjp1/dxp
         c(j,2) = one-c(j,1)
         c(j,3) = six*(dyjp1/dxjp1-dyj/dxj)/dxp
         dxj = dxjp1
         dyj = dyjp1
    5 continue
c
c     factor the tridiagonal matrix and solve for u
c
c     c(j,2)  = u(j)
c     c(j,1)  = q(j)
c     bpar(1) = lambda(1)
c     bpar(2) = d(1)
c     bpar(3) = mu(nx)
c     bpar(4) = d(nx)
c
   10 c(1,1) = -bpar(1)*half
      c(1,2) = bpar(2)*half
      if (nx .eq. 2)  go to 20
      do 15 j=2,nxm1
         pj = c(j,2)*c(j-1,1)+two
         c(j,1) = -c(j,1)/pj
         c(j,2) = (c(j,3)-c(j,2)*c(j-1,2))/pj
   15 continue
c
c     solve for cubic coefficients of spline interpolant
c     c(j,1), c(j,2), and c(j,3)
c
   20 yppb = (bpar(4)-bpar(3)*c(nxm1,2))/(bpar(3)*c(nxm1,1)+two)
      sixi = one/six
      do 25 i=1,nxm1
         j = nx-i
         yppa = c(j,1)*yppb+c(j,2)
         dx = x(j+1)-x(j)
         c(j,3) = sixi*(yppb-yppa)/dx
         c(j,2) = half*yppa
         c(j,1) = (y(j+1)-y(j))/dx-(c(j,2)+c(j,3)*dx)*dx
         yppb = yppa
   25 continue
      go to 9005
   30 ier = 129
      go to 9000
   35 ier = 130
      go to 9000
   40 ier = 131
c
 9000 call uertst1 (ier, 'icsicu1')
 9005 return
c
      end

c***** SUBROUTINE IBCIEU1 (F,IFD,X,NX,Y,NY,XL,NXL,YL,NYL,FL,IFLD,WK,IER)
C   IMSL ROUTINE NAME   - IBCIEU  - Modified for real*8
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BICUBIC SPLINE TWO-DIMENSIONAL INTERPOLATOR
C
C   USAGE               - CALL IBCIEU (F,IFD,X,NX,Y,NY,XL,NXL,YL,NYL,
C                           FL,IFLD,WK,IER)
C
C   ARGUMENTS    F      - NX BY NY MATRIX CONTAINING THE FUNCTION
C                           VALUES. (INPUT) F(I,J) IS THE FUNCTION VALUE
C                           AT THE POINT (X(I),Y(J)) FOR I=1,...,NX AND
C                           J=1,...,NY.
C                IFD    - ROW DIMENSION OF THE MATRIX F EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 2.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 2.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                XL     - VECTOR OF LENGTH NXL. (INPUT)
C                NXL    - NUMBER OF ELEMENTS IN XL. (INPUT)
C                YL     - VECTOR OF LENGTH NYL. (INPUT)
C                NYL    - NUMBER OF ELEMENTS IN YL. (INPUT)
C                         NOTE - THE COORDINATE PAIRS (XL(I),YL(J)),
C                           FOR I=1,...,NXL AND J=1,...,NYL, GIVE THE
C                           POINTS AT WHICH THE INTERPOLATORY BICUBIC
C                           SPLINE IS TO BE EVALUATED.
C                FL     - NXL BY NYL MATRIX CONTAINING THE INTERPOLATORY
C                           BICUBIC SPLINE VALUES. (OUTPUT) FL(I,J) IS
C                           SET TO THE VALUE OF THE INTERPOLATORY
C                           BICUBIC SPLINE AT (XL(I),YL(J)) FOR
C                           I=1,...,NXL AND J=1,...,NYL. NOTE THAT THE
C                           NUMBER OF COLUMNS IN FL MUST BE .GE.
C                           MAX(NYL,NY) SINCE FL IS ALSO USED AS
C                           WORKING STORAGE (OF SIZE NXL BY NY) DURING
C                           THE COMPUTATION.
C                IFLD   - ROW DIMENSION OF THE MATRIX FL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                WK     - WORK VECTOR OF LENGTH
C                           MAX((NX-1)*3,(NY-1)*3+NY).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IFD IS LESS THAN NX.
C                           IER = 130, IFLD IS LESS THAN NXL.
C                           IER = 131, NX IS LESS THAN 2.
C                           IER = 132, NY IS LESS THAN 2.
C                           IER = 133, VECTOR X IS NOT ORDERED SO THAT
C                             X(1) .LT. X(2) ... .LT. X(NX).
C                           IER = 134, VECTOR Y IS NOT ORDERED SO THAT
C                             Y(1) .LT. Y(2) ... .LT. Y(NY).
C                         WARNING ERROR
C                           IER = 37, XL(I) IS LESS THAN X(1) OR
C                             GREATER THAN X(NX).
C                           IER = 38, YL(I) IS LESS THAN Y(1) OR
C                             GREATER THAN Y(NY).
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - ICSEVU1,ICSCCU1,UERSET1,UERTST1,UGETIO1
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCIEU1 (F,IFD,X,NX,Y,NY,XL,NXL,YL,NYL,FL,IFLD,WK,IER)
      implicit integer (i-n), real*8 (a-h,o-z)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFD,NX,NY,NXL,NYL,IFLD,IER
      REAL*8               F(IFD,NY),X(NX),Y(NY),XL(NXL),YL(NYL),
     1                   FL(IFLD,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            KER,LER,MER,NXM1,NYM1,IY,JER,KYL,KYLP1,IXL,
     1                   IYL,LEVEL,LEVOLD
C                                  FIRST EXECUTABLE STATEMENT
      KER = 0
      LER = 0
C                                  CHECK IFD .GE. NX
      IF (IFD .LT. NX) GO TO 30
C                                  CHECK IFLD .GE. NX.
      IF (IFLD .LT. NXL) GO TO 35
C                                  CHECK NX .GE. 2
      IF (NX .LT. 2) GO TO 36
C                                  CHECK NY .GE. 2
      IF (NY .LT. 2) GO TO 37
C                                  CALL UERSET1 TO SILENCE WARNING
C                                    MESSAGES FROM ICSEVU
      LEVEL = 2
      CALL UERSET1 (LEVEL,LEVOLD)
      MER = 0
      NXM1 = NX-1
      NYM1 = NY-1
C                                  INTERPOLATE IN THE X-DIRECTION
      DO 10 IY=1,NY
C                                  CALCULATE THE COEFFICIENTS
         CALL ICSCCU1 (X,F(1,IY),NX,WK(1),NXM1,JER)
C                                  CHECK FOR ERROR IN ICSCCU
         IF (JER .NE. 0) GO TO 40
C                                  EVALUATE
         CALL ICSEVU1 (X,F(1,IY),NX,WK(1),NXM1,XL,FL(1,IY),NXL,JER)
C                                  CHECK FOR ERROR IN ICSEVU
         IF (JER .NE. 0) KER = 37
   10 CONTINUE
      KYL = NYM1*3
      KYLP1 = KYL+1
C                                  INTERPOLATE IN THE Y-DIRECTION
      DO 25 IXL=1,NXL
         DO 15 IY=1,NY
            WK(KYL+IY) = FL(IXL,IY)
   15    CONTINUE
C                                  CALCULATE THE COEFFICIENTS
         CALL ICSCCU1 (Y,WK(KYLP1),NY,WK(1),NYM1,JER)
C                                  CHECK FOR ERROR IN ICSCCU
         IF (JER .NE. 0) GO TO 45
C                                  EVALUATE
         DO 20 IYL=1,NYL
            CALL ICSEVU1 (Y,WK(KYLP1),NY,WK(1),NYM1,YL(IYL),FL(IXL,IYL),
     1                   1,JER)
C                                  CHECK FOR ERROR IN ICSEVU
            IF (JER .NE. 0) LER = 38
   20    CONTINUE
   25 CONTINUE
      GO TO 46
C                                  HANDLE ERRORS
   30 MER = 129
      GO TO 50
   35 MER = 130
      GO TO 50
   36 MER = 131
      GO TO 50
   37 MER = 132
      GO TO 50
   40 MER = 133
      GO TO 46
   45 MER = 134
   46 CALL UERSET1 (LEVOLD,LEVEL)
   50 IER = MAX0(MER,KER,LER)
 9000 CONTINUE
      IF (MER .NE. 0) CALL UERTST1(MER,6HIBCIEU)
      IF (KER .NE. 0) CALL UERTST1(KER,6HIBCIEU)
      IF (LER .NE. 0) CALL UERTST1(LER,6HIBCIEU)
 9005 RETURN
      END
c
c
c***** SUBROUTINE UGETIO1(IOPT,NIN,NOUT)
C   IMSL ROUTINE NAME   - UGETIO
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - AUGUST 1, 1981
C
C   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
C                           VALUES FOR INPUT AND OUTPUT UNIT
C                           IDENTIFIERS.
C
C   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
C
C   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
C                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
C                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
C                           AND NOUT, RESPECTIVELY.
C                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
C                           RESET FOR SUBSEQUENT USE.
C                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
C                           RESET FOR SUBSEQUENT USE.
C                NIN    - INPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
C                NOUT   - OUTPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
C                OPERATIONS CALLS UGETIO1 TO OBTAIN THE CURRENT UNIT
C                IDENTIFIER VALUES. IF UGETIO1 IS CALLED WITH IOPT=2 OR
C                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
C                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UGETIO1(IOPT,NIN,NOUT)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END
c
c
c***** SUBROUTINE UERSET1 (LEVEL,LEVOLD)
C   IMSL ROUTINE NAME   - UERSET1
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - CRAY/SINGLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SET MESSAGE LEVEL FOR IMSL ROUTINE UERTST1
C
C   USAGE               - CALL UERSET1 (LEVEL,LEVOLD)
C
C   ARGUMENTS    LEVEL  - NEW VALUE FOR MESSAGE LEVEL. (INPUT)
C                           OUTPUT FROM IMSL ROUTINE UERTST1 IS
C                           CONTROLLED SELECTIVELY AS FOLLOWS,
C                             LEVEL = 4 CAUSES ALL MESSAGES TO BE
C                                       PRINTED,
C                             LEVEL = 3 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 32,
C                             LEVEL = 2 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 64,
C                             LEVEL = 1 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 128,
C                             LEVEL = 0 ALL MESSAGE PRINTING IS
C                                       SUPPRESSED.
C                LEVOLD - PREVIOUS MESSAGE LEVEL. (OUTPUT)
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - UERTST1,UGETIO1
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UERSET1 (LEVEL,LEVOLD)
      implicit integer (i-n), real*8 (a-h,o-z)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LEVEL,LEVOLD
C                                  FIRST EXECUTABLE STATEMENT
      LEVOLD = LEVEL
      CALL UERTST1 (LEVOLD,6HUERSET)
      RETURN
      END
c
c
      real*8 function asimp (a1, b, ep, m, n, FUN)
      implicit integer (i-n), real*8 (a-h,o-z)
c
c --- author:  k. hillstrom (argonne national laboratory, chicago, illinois)
c
      external FUN
c
      real*8    a1, b, ep, FUN, a, eps, absar, est, fa, fm, fb, dx, sx,
     .          f1, f2, fbp, est2, nrtr, est1, sum, daft, esum, tsum, da
      real*8    aest2, ftst, fmax, aest1, delta, aest
      dimension f2(30), fbp(30), est2(30), nrtr(30), aest2(30), ftst(3)
c
c     the parameter setup for the initial call
c
      if (n .le. 0) then
        write (7,fmt='(a)')  ' **** ASIMP error return:  n .le. 0 ****'
        asimp = 0.0
        return
      end if
c
      if (n .gt. 3) then
        write (7, fmt='(a)')  ' **** ASIMP error return:  n .gt. 3 ****'
        asimp = 0.0
        return
      end if
c
      a       = a1
      eps     = ep*15.0
      esum    = 0.0
      tsum    = 0.0
      lvl     = 1
      da      = b-a
      fa      = FUN(a)
      fm      = FUN((a+b)*0.5)
      fb      = FUN(b)
      m       = 3
      fmax    = ABS (fa)
      ftst(1) = fmax
      ftst(2) = ABS (fm)
      ftst(3) = ABS (fb)
      do 10 i=2,3
        if (fmax .ge. ftst(i))  go to 10
        fmax = ftst(i)
   10 continue
      est   = (fa+4.0*fm+fb)*da/6.0
      absar = (ftst(1)+4.0*ftst(2)+ftst(3))*da/6.0
      aest  = absar
c
c 1 = recur
c
   20 dx         = da/(2.0**lvl)
      sx         = dx/6.0
      f1         = FUN(a+0.5*dx)
      f2(lvl)    = FUN(a+1.5*dx)
      est1       = sx*(fa+4.0*f1+fm)
      fbp(lvl)   = fb
      est2(lvl)  = sx*(fm+4.0*f2(lvl)+fb)
      sum        = est1+est2(lvl)
      ftst(1)    = ABS (f1)
      ftst(2)    = ABS (f2(lvl))
      ftst(3)    = ABS (fm)
      aest1      = sx*(ABS (fa)+4.0*ftst(1)+ftst(3))
      aest2(lvl) = sx*(ftst(3) +4.0*ftst(2) + ABS (fb))
      absar      = absar-aest+aest1+aest2(lvl)
      m          = m+2
      go to (60,30,70),n
   30 delta = absar
      go to 90
   60 delta = 1.0
      go to 90
   70 do 80 i=1,2
        if (fmax .ge. ftst(i))  go to 80
        fmax = ftst(i)
   80 continue
      delta = fmax
   90 diff  = ABS (est-sum)
      daft  = (est-sum)/15.0
      if (diff-eps*delta) 110,110,100
  100 if (lvl-30) 140,120,120
  110 if (lvl- 1) 120,140,120
c
c 2 = up
c
  120 a    = a+2.0*dx
  130 lvl  = lvl-1
      esum = esum+daft
      l    = nrtr(lvl)
      tsum = tsum+sum
      go to (160, 170), l
c
c 11 = r1,12=r2
c
  140 nrtr(lvl) = 1
      est = est1
      aest = aest1
      fb = fm
      fm = f1
      eps = eps/2.0
  150 lvl = lvl+1
      go to 20
  160 nrtr(lvl) = 2
      fa = fb
      fm = f2(lvl)
      fb = fbp(lvl)
      est = est2(lvl)
      aest = aest2(lvl)
      go to 150
  170 eps = 2.0 * eps
      sum = 0.0
      if (lvl-1) 180,180,130
  180 asimp = tsum
      a = ABS (esum)
      ep = diff/delta
      if (a .ge. ep)  go to 190
      asimp = asimp - esum
  190 return
c
      end


C---- IMSL ROUTINE -------------------------------------------------------------

c     subroutine icsicu1 (x, y, nx, bpar, c, ic, ier)
c --- modified version of IMSL subroutine named ICSICU (REAL->REAL*8)

c called from:
c frsubs.f(1922):  CALL ICSCCU1 (X,F(1,IY),NX,WK(1),NXM1,JER)
c frsubs.f(1938):  CALL ICSCCU1 (Y,WK(KYLP1),NY,WK(1),NYM1,JER)

C----------------------------------------------------------------------
C   IMSL ROUTINE NAME   - ICSCCU                                        
c   CUBIC SPLINE INTERPOLATION (EASY-TO-USE VERSION)
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/SINGLE                                    
C                                                                       
C   LATEST REVISION     - JUNE 1, 1980                                  
C                                                                       
C   PURPOSE             - CUBIC SPLINE INTERPOLATION                    
C                           (EASY-TO-USE VERSION)                       
C                                                                       
C   USAGE               - CALL ICSCCU (X,Y,NX,C,IC,IER)                 
C                                                                       
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE  
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,  
C                           NX. (INPUT) X MUST BE ORDERED SO THAT       
C                           X(I) .LT. X(I+1).                           
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES  
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS. 
C                           (INPUT)                                     
C                NX     - NUMBER OF ELEMENTS IN X AND Y. (INPUT) NX     !210
C                           MUST BE .GE. 2.                             !220
C                C      - SPLINE COEFFICIENTS. (OUTPUT) C IS AN NX-1 BY !230
C                           3 MATRIX. THE VALUE OF THE SPLINE           !240
C                           APPROXIMATION AT T IS                       !250
C                           S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)  !260
C                           WHERE X(I) .LE. T .LT. X(I+1) AND           !270
C                           D = T-X(I).                                 !280
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS          !290
C                           SPECIFIED IN THE DIMENSION STATEMENT IN     !300
C                           THE CALLING PROGRAM. (INPUT)                !310
C                IER    - ERROR PARAMETER. (OUTPUT)                     !320
C                         TERMINAL ERROR                                !330
C                           IER = 129, IC IS LESS THAN NX-1.            !340
C                           IER = 130, NX IS LESS THAN 2.               !350
C                           IER = 131, INPUT ABSCISSA ARE NOT ORDERED   !360
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX).    !370
C                                                                       !380
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         !390
C                       - SINGLE/H36,H48,H60                            !400
C                                                                       !410
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 !420
C                                                                       !430
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           !440
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      !450
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  !460
C                                                                       !470
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.       !480
C                                                                       !490
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN !500
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    !510
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        !520
C                                                                       !530
C-----------------------------------------------------------------------!540
C                                                                      
      SUBROUTINE ICSCCU1 (X,Y,NX,C,IC,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS         !570
      INTEGER            NX,IC,IER                                      !580
      REAL*8             X(NX),Y(NX),C(IC,3) !-YuP: changed REAL->REAL*8
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   !600
      INTEGER            IM1,I,JJ,J,MM1,MP1,M,NM1,NM2                   !610
      REAL*8             DIVDF1,DIVDF3,DTAU,G,CNX(3) !-YuP: changed REAL->REAL*8
C                                  FIRST EXECUTABLE STATEMENT           !630
      NM1 = NX-1                                                        !640
      IER = 129                                                         !650
      IF (IC .LT. NM1) GO TO 9000                                       !660
      IER = 130                                                         !670
      IF (NX .LT. 2) GO TO 9000                                         !680
      IER = 131                                                         !690
      IF (NX .EQ. 2) GO TO 45                                           !700
C                                  COMPUTE NOT-A-KNOT SPLINE            !710
      DO 5 M = 2,NM1                                                    !720
         MM1=M-1                                                        !730
         C(M,2) = X(M)-X(MM1)                                           !740
         IF (C(M,2).LE.0.0) GO TO 9000                                  !750
         C(M,3) = (Y(M)-Y(MM1))/C(M,2)                                  !760
    5 CONTINUE                                                          !770
      CNX(2) = X(NX)-X(NM1)                                             !780
      IF (CNX(2).LE.0.0) GO TO 9000                                     !790
      CNX(3) = (Y(NX)-Y(NM1))/CNX(2)                                    !800
      IER = 0                                                           !810
      NM2 = NX-2                                                        !820
      IF (NX .GT. 3) GO TO 10                                           !830
      C(1,3) = CNX(2)                                                   !840
      C(1,2) = C(2,2)+CNX(2)                                            !850
      C(1,1) = ((C(2,2)+2.*C(1,2))*C(2,3)*CNX(2)+C(2,2)**2*CNX(3))      !860
     1/C(1,2)                                                           !870
      GO TO 20                                                          !880
   10 C(1,3) = C(3,2)                                                   !890
      C(1,2) = C(2,2)+C(3,2)                                            !900
      C(1,1) = ((C(2,2)+2.*C(1,2))*C(2,3)*C(3,2)+C(2,2)**2*C(3,3))      !910
     1/C(1,2)                                                           !920
      DO 15 M=2,NM2                                                     !930
         MP1=M+1                                                        !940
         MM1=M-1                                                        !950
         G = -C(MP1,2)/C(MM1,3)                                         !960
         C(M,1) = G*C(MM1,1)+3.*C(M,2)*C(MP1,3)+3.*C(MP1,2)*C(M,3)      !970
         C(M,3) = G*C(MM1,2)+2.*C(M,2)+2.*C(MP1,2)                      !980
   15 CONTINUE                                                          !990
   20 G = -CNX(2)/C(NM2,3)                                              !1000
      C(NM1,1) = G*C(NM2,1)+3.*C(NM1,2)*CNX(3)+3.*CNX(2)*C(NM1,3)       !1010
      C(NM1,3) = G*C(NM2,2)+2.*C(NM1,2)+2.*CNX(2)                       !1020
      IF (NX.GT.3) GO TO 25                                             !1030
      CNX(1)=2.*CNX(3)                                                  !1040
      CNX(3)=1.                                                         !1050
      G=-1./C(NM1,3)                                                    !1060
      GO TO 30                                                          !1070
   25 G = C(NM1,2)+CNX(2)                                               !1080
      CNX(1) = ((CNX(2)+2.*G)*CNX(3)*C(NM1,2)+CNX(2)**2*                !1090
     1(Y(NM1)-Y(NX-2))/C(NM1,2))/G                                      !1100
      G = -G/C(NM1,3)                                                   !1110
      CNX(3) = C(NM1,2)                                                 !1120
   30 CNX(3) = G*C(NM1,2)+CNX(3)                                        !1130
      CNX(1) = (G*C(NM1,1)+CNX(1))/CNX(3)                               !1140
      C(NM1,1) = (C(NM1,1)-C(NM1,2)*CNX(1))/C(NM1,3)                    !1150
      DO 35 JJ=1,NM2                                                    !1160
         J = NM1-JJ                                                     !1170
         C(J,1) = (C(J,1)-C(J,2)*C(J+1,1))/C(J,3)                       !1180
   35 CONTINUE                                                          !1190
      DO 40 I=2,NM1                                                     !1200
         IM1 = I-1                                                      !1210
         DTAU = C(I,2)                                                  !1220
         DIVDF1 = (Y(I)-Y(IM1))/DTAU                                    !1230
         DIVDF3 = C(IM1,1)+C(I,1)-2.*DIVDF1                             !1240
         C(IM1,2) = (DIVDF1-C(IM1,1)-DIVDF3)/DTAU                       !1250
         C(IM1,3) = DIVDF3/DTAU**2                                      !1260
   40 CONTINUE                                                          !1270
      DTAU = CNX(2)                                                     !1280
      DIVDF1 = (Y(NX)-Y(NM1))/DTAU                                      !1290
      DIVDF3 = C(NM1,1)+CNX(1)-2.*DIVDF1                                !1300
      C(NM1,2) = (DIVDF1-C(NM1,1)-DIVDF3)/DTAU                          !1310
      C(NM1,3) = DIVDF3/DTAU**2                                         !1320
      GO TO 9005                                                        !1330
   45 IF (X(1) .GE. X(2)) GO TO 9000                                    !1340
      IER = 0                                                           !1350
      C(1,1) = (Y(2)-Y(1))/(X(2)-X(1))                                  !1360
      C(1,2) = 0.0                                                      !1370
      C(1,3) = 0.0                                                      !1380
      GO TO 9005                                                        !1390
 9000 CONTINUE                                                          !1400
      CALL UERTST(IER,8HICSCCU  )                                       !1410
 9005 RETURN                                                            !1420
      END                                                               !1430



c---IMSL---UERTST------------------------------------------------IMSL---UERTST--    
C   USAGE               - CALL UERTST (IER,NAME)
C   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
C                           IER = I+J WHERE
C                             I = 128 IMPLIES TERMINAL ERROR,
C                             I =  64 IMPLIES WARNING WITH FIX, AND
C                             I =  32 IMPLIES WARNING.
C                             J = ERROR CODE RELEVANT TO CALLING
C                                 ROUTINE.
C                NAME   - A SIX CHARACTER LITERAL STRING GIVING THE
C                           NAME OF THE CALLING ROUTINE. (INPUT)
C   PRECISION/HARDWARE  - SINGLE/ALL
C   REQD. IMSL ROUTINES - UGETIO
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
C                ONTO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
C                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
C                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
C                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
C                UGETIO AS FOLLOWS..
C                                NIN = 0
C                                NOUT = NEW OUTPUT UNIT NUMBER
C                                CALL UGETIO(3,NIN,NOUT)
C                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C-----------------------------------------------------------------------
      SUBROUTINE UERTST (IER,NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      INTEGER*2          NAME(3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER*2          NAMSET(3),NAMEQ(3)
	integer*4 ieq, i, ieqdf, levold, level
      DATA               NAMSET/2HUE,2HRS,2HET/
      DATA               NAMEQ/2H  ,2H  ,2H  /
C                                  FIRST EXECUTABLE STATEMENT
      DATA               LEVEL/4/,IEQDF/0/,IEQ/1H=/
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(*,35) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(*,35) IER,NAME
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C                                  PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(*,40) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(*,40) IER,NAME
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
C                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(*,45) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(*,45) IER,NAME
      GO TO 30
   15 CONTINUE
C                                  CHECK FOR UERSET CALL
      DO 20 I=1,3
         IF (NAME(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C                                  PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(*,50) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(*,50) IER,NAME
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
   40 FORMAT(36H *** WARNING WITH FIX ERROR  (IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     1       20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
C                                  SAVE P FOR P = R CASE
C                                    P IS THE PAGE NAME
C                                    R IS THE ROUTINE NAME
   55 IEQDF = 1
      DO 60 I=1,3
   60 NAMEQ(I) = NAME(I)
   65 RETURN
      END

    
!======================================================================
!   IMSL ROUTINE NAME   - ZBRENT
! YuP[2020-07-13] Imported by Yu.Petrov, changed to real*8
!-----------------------------------------------------------------------
!
!   COMPUTER            - DEC11/SINGLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
!                           GIVEN INTERVAL (BRENT ALGORITHM)
!
!   USAGE               - CALL ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)
!
!   ARGUMENTS    F      - AN EXTERNAL FUNCTION SUBPROGRAM F(X)
!                           PROVIDED BY THE USER WHICH COMPUTES F FOR
!                           ANY X IN THE INTERVAL (A,B). (INPUT)
!                           F MUST APPEAR IN AN EXTERNAL STATEMENT IN
!                           THE CALLING PROGRAM
!                EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT,
!                           B, IS ACCEPTED IF ABS(F(B)) IS LESS THAN OR
!                           EQUAL TO EPS.  EPS MAY BE SET TO ZERO.
!                NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
!                           B, IS ACCEPTED IF THE CURRENT APPROXIMATION
!                           AGREES WITH THE TRUE SOLUTION TO NSIG
!                           SIGNIFICANT DIGITS.
!                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
!                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE
!                           IN SIGN.
!                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B
!                           WILL CONTAIN THE BEST APPROXIMATION TO THE
!                           ROOT OF F. SEE REMARK 1.
!                MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND
!                           ON THE NUMBER OF FUNCTION EVALUATIONS
!                           REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN
!                           WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
!                           EVALUATIONS USED.
!                IER    - ERROR PARAMETER. (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129 INDICATES THE ALGORITHM FAILED TO
!                             CONVERGE IN MAXFN EVALUATIONS.
!                           IER = 130 INDICATES F(A) AND F(B) HAVE THE
!                             SAME SIGN.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS  1.  ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE
!                FOLLOWING,
!                F(A)*F(B) .LE.0,
!                ABS(F(B)) .LE. ABS(F(A)), AND
!                EITHER ABS(F(B)) .LE. EPS OR
!                ABS(A-B) .LE. MAX(ABS(B),0.1)*10.0**(-NSIG).
!                THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES
!                LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE
!                COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED
!                IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL
!                MAGNITUDE.
!            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN
!                K = (ALOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
!                  D=MIN(OVER X IN (A,B) OF
!                    MAX(ABS(X),0.1)*10.0**(-NSIG)).
!                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS.
!                RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY
!                ZBRENT EXCEED SQRT(K). D CAN BE COMPUTED AS FOLLOWS,
!                  P = AMIN1(ABS(A),ABS(B))
!                  P = AMAX1(0.1,P)
!                  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1
!                  D = P*10.0**(-NSIG)
!
!   COPYRIGHT           - 1977 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZBRENT(F,EPS,NSIG,A,B,MAXFN,IER)
      implicit none
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,MAXFN,IER
      REAL*8             F,EPS,A,B
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IC
      REAL*8             ZERO,HALF,ONE,THREE,TEN,
     &                   T,FA,FB,C,FC,D,E,TOL,RM,S,P,Q,R,RONE,TEMP

      DATA               ZERO/0.d0/,HALF/.5d0/,ONE/1.d0/,THREE/3.d0/,

     &                   TEN/10.d0/

!                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      T = TEN**(-NSIG)
      IC = 2
      S = A
      FA = F(S)
      S = B
      FB = F(S)
!                                  TEST FOR SAME SIGN
      IF (FA*FB.GT.ZERO) GO TO 50
    5 C = A
      FC = FA
      D = B-C
      E = D
   10 IF (ABS(FC).GE.ABS(FB)) GO TO 15
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
   15 CONTINUE
      TOL = T*MAX(ABS(B),0.1d0)
      RM = (C-B)*HALF
!                                  TEST FOR FIRST CONVERGENCE CRITERIA
      IF (ABS(FB).LE.EPS) GO TO 40
!                                  TEST FOR SECOND CONVERGENCE CRITERIA
      IF (ABS(C-B).LE.TOL) GO TO 40
!                                  CHECK EVALUATION COUNTER
      IF (IC.GE.MAXFN) GO TO 45
!                                  IS BISECTION FORCED
      IF (ABS(E).LT.TOL) GO TO 30
      IF (ABS(FA).LE.ABS(FB)) GO TO 30
      S = FB/FA
      IF (A.NE.C) GO TO 20
!                                  LINEAR INTERPOLATION
      P = (C-B)*S
      Q = ONE-S
      GO TO 25
!                                  INVERSE QUADRATIC INTERPOLATION
   20 Q = FA/FC
      R = FB/FC
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
   25 IF (P.GT.ZERO) Q = -Q
      IF (P.LT.ZERO) P = -P
      S = E
      E = D
!                                  IF ABS(P/Q).GE.75*ABS(C-B) THEN
!                                     FORCE BISECTION
      IF (P+P.GE.THREE*RM*Q) GO TO 30
!                                  IF ABS(P/Q).GE..5*ABS(S) THEN FORCE
!                                     BISECTION. S = THE VALUE OF P/Q
!                                     ON THE STEP BEFORE THE LAST ONE
      IF (P+P.GE.ABS(S*Q)) GO TO 30
      D = P/Q
      GO TO 35
!                                  BISECTION
   30 E = RM
      D = E
!                                  INCREMENT B
   35 A = B
      FA = FB
      TEMP = D
      IF (ABS(TEMP).LE.HALF*TOL) TEMP = SIGN(HALF*TOL,RM)
      B = B+TEMP
      S = B
      FB = F(S)
      IC = IC+1
      IF (FB*FC.LE.ZERO) GO TO 10
      GO TO 5
!                                  CONVERGENCE OF B
   40 A = C
      MAXFN = IC
      GO TO 9005
!                                  MAXFN EVALUATIONS
   45 IER = 129
      A = C
      MAXFN = IC
      GO TO 9000
!                                  TERMINAL ERROR - F(A) AND F(B) HAVE
!                                  THE SAME SIGN
   50 IER = 130
      MAXFN = IC
 9000 CONTINUE
      CALL UERTST1(IER,'ZBRENT')
 9005 RETURN
      END SUBROUTINE ZBRENT



     