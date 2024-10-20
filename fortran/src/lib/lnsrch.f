      SUBROUTINE LNSRCH (M,C0,W,PHI0,DC,FMTOL,OPT, S, C,DS,
     .                   D1C,D2C,U,F,DF,HU,CL,E,FNS,PHIC,
     .                   NEVAL,CMAX,DELC,IER)
      INTEGER M, NEVAL, IER
      LOGICAL OPT
      DOUBLE PRECISION C0(3,0:M), W, PHI0, DC(3,M), FMTOL,
     .                 S, C(3,0:M), DS(M), D1C(3,M),
     .                 D2C(3,M), U(3,M), F(M), DF(3,M),
     .                 HU(3,M), CL, E, FNS, PHIC, CMAX, DELC
C
C***********************************************************
C
C                                                  From GEO2
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/19/03
C
C   This subroutine minimizes phi(C0+S*DC) over positive
C step-sizes S (unless OPT = FALSE), where phi is the
C functional defined by Function PHI.
C
C On input:
C
C       M = Number of vertices and curve segments.  M >= 3.
C
C       C0 = Array dimensioned 3 by 0:M containing the
C            current estimate of the minimizer of phi:
C            sequence of vertices with the Cartesian
C            coordinates of c(j) in column j.
C
C       W = Weight defining the penalty term.  W >= 0.
C
C       PHI0 = phi(C0).
C
C       DC = Array dimensioned 3 by M defining the search
C            direction for the descent method.
C
C       FMTOL = Nonnegative tolerance for Function FMIN:
C               the desired length of the interval of
C               uncertainty for the optimal step-size.
C
C       OPT = Flag with value TRUE iff the optimal step-
C             size is to be computed.
C
C The above parameters are not altered by this routine.
C
C       S = Step-size to be used if OPT = FALSE, or initial
C           estimate of the optimal step-size if OPT = TRUE.
C           The output from a previous call is a reasonable
C           choice for the initial estimate.  S > 0.
C
C       C = Array dimensioned 3 by 0:M.
C
C       DS = Array of length at least M.
C
C       D1C,D2C,U = Arrays dimensioned 3 by M.
C
C       F = Array of length at least M.
C
C       DF,HU = Array dimensioned 3 by M.
C
C On output:
C
C       S = Estimate of the optimal step-size (not altered
C           if OPT = FALSE).
C
C       C = Minimizer of phi:  C0 + S*DC
C
C       DS,D1C,D2C,U,F,DF,HU,CL,E,FNS = Terms computed by
C                                       Function PHI (and
C                                       associated with C).
C
C       PHIC = Objective function value phi(C).
C
C       NEVAL = Number of evaluations of phi.
C
C       CMAX = Max-norm of C.
C
C       DELC = Relative change in C:  Max-norm(C-C0)/(1+CMAX)
C            = S*Max-norm(DC)/(1+CMAX).
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M, W, FMTOL, or S is outside its
C                     valid range on input.
C             IER = 2 if duplicate vertices are found in
C                     Function PHI.
C             IER = 3 if an invalid gradient grad(f) = 0 is
C                     returned by FEVAL (called by PHI).
C             IER = 4 if DC = 0.
C             Output parameters are not altered if IER = 1.
C
C Modules required by LNSRCH:  FEVAL, FMIN, PHI
C
C Intrinsic functions called by LNSRCH:  ABS, MAX
C
C***********************************************************
C
      DOUBLE PRECISION PS2, PSS, S1, S2
      DOUBLE PRECISION FMIN, PHI
      INTEGER IFLG, J, NEV
C
C Test for invalid input.
C
      IF (M .LT. 3  .OR.  W .LT. 0.D0  .OR.  FMTOL .LT. 0.D0
     .    .OR.  S .LE. 0.D0) THEN
        NEVAL = 0
        IER = 1
        RETURN
      ENDIF
C
C Compute DELC = Max-norm(DC), and test for DC = 0.
C
      DELC = 0.
      DO 1 J = 1,M
        IF (ABS(DC(1,J)) .GT. DELC) DELC = ABS(DC(1,J))
        IF (ABS(DC(2,J)) .GT. DELC) DELC = ABS(DC(2,J))
        IF (ABS(DC(3,J)) .GT. DELC) DELC = ABS(DC(3,J))
    1   CONTINUE
      IF (DELC .EQ. 0.) THEN
        NEVAL = 0
        IER = 4
        RETURN
      ENDIF
      NEV = 0
C
      IF (.NOT. OPT) GO TO 6
C
C Find a bracketing interval [S1,S2] = [0,S2] such that
C   phi(C0+S2*DC) > PHI0 = phi(C0).  S2 is initialized to
C   2*S and doubled at each step as necessary.
C
      S1 = 0.D0
      S2 = S
    2 S2 = 2.D0*S2
        DO 3 J = 1,M
          C(1,J) = C0(1,J) + S2*DC(1,J)
          C(2,J) = C0(2,J) + S2*DC(2,J)
          C(3,J) = C0(3,J) + S2*DC(3,J)
    3     CONTINUE
        C(1,0) = C(1,M)
        C(2,0) = C(2,M)
        C(3,0) = C(3,M)
        PS2 = PHI (M,C,W, DS,D1C,D2C,U,F,DF,HU,CL,E,FNS,IER)
        NEV = NEV + 1
        IF (IER .NE. 0) RETURN
        IF (PS2 .LE. PHI0) GO TO 2
C
C Compute the optimal step-size S.
C
      IFLG = 0
      PSS = 0.D0
    4 S = FMIN (S1,S2,FMTOL,PSS, IFLG )
        IF (IFLG .EQ. 0) GO TO 6
        DO 5 J = 1,M
          C(1,J) = C0(1,J) + S*DC(1,J)
          C(2,J) = C0(2,J) + S*DC(2,J)
          C(3,J) = C0(3,J) + S*DC(3,J)
    5     CONTINUE
        C(1,0) = C(1,M)
        C(2,0) = C(2,M)
        C(3,0) = C(3,M)
        PSS = PHI (M,C,W, DS,D1C,D2C,U,F,DF,HU,CL,E,FNS,IER)
        NEV = NEV + 1
        IF (IER .NE. 0) RETURN
        GO TO 4
C
C Update C and compute the final value PHIC = phi(C), CMAX,
C  and DELC = Max-norm(C-C0)/(1+CMAX) for C = C0+S*DC.
C
    6 CMAX = MAX( ABS(C(1,0)), ABS(C(2,0)), ABS(C(3,0)) )
      DO 7 J = 1,M
        C(1,J) = C0(1,J) + S*DC(1,J)
        C(2,J) = C0(2,J) + S*DC(2,J)
        C(3,J) = C0(3,J) + S*DC(3,J)
        IF (ABS(C(1,J)) .GT. CMAX) CMAX = ABS(C(1,J))
        IF (ABS(C(2,J)) .GT. CMAX) CMAX = ABS(C(2,J))
        IF (ABS(C(3,J)) .GT. CMAX) CMAX = ABS(C(3,J))
    7   CONTINUE
      C(1,0) = C(1,M)
      C(2,0) = C(2,M)
      C(3,0) = C(3,M)
      PHIC = PHI (M,C,W, DS,D1C,D2C,U,F,DF,HU,CL,E,FNS,IER)
      DELC = S*DELC/(1.D0+CMAX)
      NEVAL = NEV + 1
      RETURN
      END
