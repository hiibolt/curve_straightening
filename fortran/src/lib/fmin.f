      DOUBLE PRECISION FUNCTION FMIN (AX,BX,TOL,F, IFLG )
      DOUBLE PRECISION AX, BX, TOL, F
      INTEGER IFLG
C
C***********************************************************
C
C                                                  From GLRH
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/23/95
C
C   This function returns an approximation to the point at
C which a real-valued function F attains a minimum on the
C interval (AX,BX).  It is based on the function by the same
C name by R. Brent (01/01/73), modified to use reverse
C communication.
C
C   The method used is a combination of golden section
C search and successive parabolic interpolation.  Conver-
C gence is never much slower than that for a Fibonacci
C search.  If F has a continuous second derivative which is
C positive at the minimum (which is not at AX or BX), then
C convergence is superlinear, and usually of the order of
C about 1.324....
C
C   The function F is never evaluated at two points closer
C together than EPS*Abs(FMIN) + (TOL/3), where EPS is
C approximately the square root of the relative machine
C precision.  If F is a unimodal function and the computed
C values of F are always unimodal when separated by at least
C EPS*Abs(X*) + (TOL/3), then FMIN approximates the abscissa
C of the global minimum of F on the interval (AX,BX) with
C an error less than 3*EPS*Abs(FMIN) + TOL.  If F is not
C unimodal, then FMIN may approximate a local, but perhaps
C non-global, minimum to the same accuracy.
C
C On input:
C
C       AX,BX = Endpoints of the initial interval -- the
C               interval over which F is to be minimized.
C
C       TOL = Desired length of the interval of uncertainty
C             of the final result.  TOL .GE. 0.
C
C       F = Function value F(FMIN) if IFLG .NE. 0, or unused
C           dummy parameter on the first call.
C
C The above parameters are not altered by this function.
C
C       IFLG = Reverse communication flag:
C              IFLG = 0 if this is the first call to FMIN
C                       for a given minimization problem.
C              IFLG > 0 if F contains a function value
C                       at the point FMIN returned by the
C                       previously call.  The value of
C                       IFLG must be the value returned
C                       by the previous call.
C
C On Output:
C
C       IFLG = Reverse communication flag:
C              IFLG = 0 if FMIN is the solution.
C              IFLG > 0 if FMIN contains a point at which
C                       F is to be evaluated.  FMIN must
C                       be called again.
C              IFLG < 0 if FMIN was called with an invalid
C                       value of IFLG (neither 0 nor the
C                       returned value).
C
C       FMIN = Approximation to the point at which F is
C              minimized (IFLG = 0), or point at which F
C              is to be evaluated (IFLG > 0).
C
C Module required by FMIN:  DSTORE
C
C Intrinsic functions called by FMIN:  ABS, SIGN, SQRT
C
C Reference:  Richard Brent, Algorithms for Minimization
C             Without Derivatives, Prentice-Hall, Inc.
C             (1973).
C
C***********************************************************
C
      DOUBLE PRECISION A, B, C, D, E, FU, FV, FW, FX, EPS,
     .                 P, Q, R, TOL1, TOL2, U, V, W, X, XM
      DOUBLE PRECISION DSTORE
      SAVE
C
C Test IFLG.
C
      IF (IFLG .EQ. 1) GO TO 2
      IF (IFLG .EQ. 2) GO TO 6
      IF (IFLG .NE. 0) THEN
        IFLG = -ABS(IFLG)
        RETURN
      ENDIF
C
C C is the squared inverse of the golden ratio.
C
      C = 0.5D0*(3.D0 - SQRT(5.D0))
C
C EPS is approximately the square root of the relative
C   machine precision.
C
      EPS = 1.D0
    1 EPS = EPS/2.D0
        TOL1 = DSTORE(1.D0 + EPS)
        IF (TOL1 .GT. 1.D0) GO TO 1
      EPS = SQRT(EPS)
C
C Initialization:
C
      A = AX
      B = BX
      V = A + C*(B - A)
      W = V
      X = V
      E = 0.D0
C
C Get F(X):  IFLG = 1.
C
      IFLG = 1
      FMIN = X
      RETURN
    2 FX = F
      FV = FX
      FW = FX
C
C Main loop:
C
    3 XM = 0.5D0*(A + B)
        TOL1 = EPS*ABS(X) + TOL/3.D0
        TOL2 = 2.D0*TOL1
C
C Test for termination.
C
        IF (ABS(X - XM) .LE. (TOL2 - 0.5D0*(B - A)))
     .    GO TO 7
C
C Test for golden-section necessary.
C
        IF (ABS(E) .LE. TOL1) GO TO 4
C
C Fit a parabola.
C
        R = (X - W)*(FX - FV)
        Q = (X - V)*(FX - FW)
        P = (X - V)*Q - (X - W)*R
        Q = 2.D0*(Q - R)
        IF (Q .GT. 0.D0) P = -P
        Q =  ABS(Q)
        R = E
        E = D
C
C Test for parabola acceptable.
C
        IF ( ABS(P) .GE. ABS(0.5D0*Q*R)  .OR.
     .       P .LE. Q*(A - X) .OR.  P .GE. Q*(B - X) )
     .    GO TO 4
C
C Take a parabolic interpolation step.
C
        D = P/Q
        U = X + D
C
C F must not be evaluated too close to AX or BX.
C
        IF ((U - A) .LT. TOL2) D = SIGN(TOL1, XM - X)
        IF ((B - U) .LT. TOL2) D = SIGN(TOL1, XM - X)
        GO TO 5
C
C Take a golden-section step.
C
    4   IF (X .GE. XM) E = A - X
        IF (X .LT. XM) E = B - X
        D = C*E
C
C F must not be evaluated too close to X.
C
    5   IF (ABS(D) .GE. TOL1) THEN
          U = X + D
        ELSE
          U = X + SIGN(TOL1, D)
        ENDIF
C
C Get F(U).  IFLG = 2.
C
        IFLG = 2
        FMIN = U
        RETURN
    6   FU = F
C
C Update A, B, V, W, and X.
C
        IF (FU .LE. FX) THEN
          IF (U .GE. X) A = X
          IF (U .LT. X) B = X
          V = W
          FV = FW
          W = X
          FW = FX
          X = U
          FX = FU
          GO TO 3
        ENDIF
C
        IF (U .LT. X) THEN
          A = U
        ELSE
          B = U
        ENDIF
        IF (FU .LE. FW  .OR.  W .EQ. X) THEN
          V = W
          FV = FW
          W = U
          FW = FU
          GO TO 3
        ENDIF
        IF (FU .LE. FV  .OR.  V .EQ. X  .OR.  V .EQ. W) THEN
          V = U
          FV = FU
        ENDIF
        GO TO 3
C
C Return the solution.
C
    7 IFLG = 0
      FMIN = X
      RETURN
      END