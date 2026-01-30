      DOUBLE PRECISION FUNCTION ZABS1(ZR, ZI)
C***BEGIN PROLOGUE  ZABS1
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
C     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZABS1
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
      U = DABS(ZR)
      V = DABS(ZI)
      S = U + V
C-----------------------------------------------------------------------
C     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
C     TRUE FLOATING ZERO
C-----------------------------------------------------------------------
      S = S*1.0D+0
      IF (S.EQ.0.0D+0) GO TO 20
      IF (U.GT.V) GO TO 10
      Q = U/V
      ZABS1 = V*DSQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZABS1 = U*DSQRT(1.D+0+Q*Q)
      RETURN
   20 ZABS1 = 0.0D+0
      RETURN
      END
