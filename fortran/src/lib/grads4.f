      SUBROUTINE GRADS4 (M,DS,A0,RCMIN, AP,G,DA,WK, RCOND,
     .                   IER)
      INTEGER M, IER
      DOUBLE PRECISION DS(M), A0, RCMIN, AP(*), G(3,M),
     .                 DA(M), WK(M,3), RCOND
C
C***********************************************************
C
C                                                  From GEO2
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/30/03
C
C   This subroutine applies a 4-th order smoother DTD**(-1)
C to the L2 gradient grad(phi), returning a Sobolev gradient
C G of the functional phi associated with a periodic geo-
C desic curve (Function PHI).  The linear system is solved
C by Linpack Subroutines DPPCO and DPPSL.  Refer to the
C header comments for the definition of DTD = D**T*D.
C
C On input:
C
C       M = Number of vertices and curve segments.  M >= 5.
C
C       DS = Array of length M containing segment lengths
C            computed by Function PHI (and associated with
C            the current sequence of vertices C).
C
C       A0 = Positive weight defining the Sobolev inner
C            product <g,h>_c and the differential operator
C            D.  A small value places more weight on the
C            second derivatives resulting in a more effec-
C            tive smoother DTD**(-1), but a larger value may
C            be necessary to improve the condition number,
C            especially for large M.
C
C       RCMIN = Minimum value of RCOND (reciprocal of the
C               estimated condition number of DTD) for which
C               the linear system is to be solved.  A
C               reasonable value is 1.e3*EPS, where EPS is
C               the machine precision, guaranteeing at least
C               3 correct digits in the solution.
C               RCMIN > 0.
C
C The above parameters are not altered by this routine.
C
C       AP = Array of length at least M*(M+1)/2 used to
C            store the symmetric positive definite order-M
C            matrix D**T*D.
C
C       G = Array dimensioned 3 by M containing the compo-
C           nents of the L2-gradient grad(phi) (partial
C           derivatives of phi with respect to the vertices
C           c(j)).  G must be an element of S0 and may be
C           computed by Subroutine GRADL2.
C
C       DA = Work array of length at least M.
C
C       WK = Work array of length at least 3*M.
C
C On output:
C
C       AP = R**T*R factorization of the upper triangle of
C            D**T*D in Linpack packed storage format
C            (incomplete if IER = 2).
C
C       G = Sobolev gradient (an element of S0) unless
C           IER > 0, in which case G is not altered.
C
C       DA = Array of arc lengths:  DA(i) = 2*da(i) =
C            ds(i)+ds(i+1) for i = 1 to M.
C
C       WK = Garbage.
C
C       RCOND = Estimated reciprocal of the condition number
C               of DTD (0 if IER = 2) unless IER = 1, in
C               which case RCOND is not altered.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M, A0, or RCMIN is outside its
C                     valid range on input.  Output parame-
C                     ters are not altered in this case.
C             IER = 2 if DTD is not positive definite
C                     (singular to working precision).
C             IER = 3 if RCOND < RCMIN.
C
C Linpack subprograms required by GRADS4:  DASUM, DAXPY,
C                                        DDOT, DPPCO, DPPFA,
C                                        DPPSL, DSCAL
C
C***********************************************************
C
      DOUBLE PRECISION W0
      INTEGER I, J, K
C
C Test for invalid input parameters.
C
      IF (M .LT. 5  .OR.  A0 .LE. 0.D0  .OR.
     .    RCMIN .LE. 0.D0) GO TO 21
C
C  Store da(i)*2 = ds(i)+ds(i+1), and the coefficients of
C
C       D2c(i) = w1(i)*c(i-1) - w2(i)*c(i) + w3(i)*c(i+1),
C
C where
C
C       w1(i) = 2/[ds(i)*(ds(i)+ds(i+1))],
C
C       w2(i) = 2/[ds(i)*ds(i+1)],
C
C       w3(i) = 2/[ds(i+1)*(ds(i)+ds(i+1))]
C
C for i = 1 to m.
C
      DO 1 I = 1,M-1
        DA(I) = DS(I) + DS(I+1)
        WK(I,1) = 2.D0/(DS(I)*DA(I))
        WK(I,2) = 2.D0/(DS(I)*DS(I+1))
        WK(I,3) = 2.D0/(DS(I+1)*DA(I))
    1   CONTINUE
      DA(M) = DS(M) + DS(1)
      WK(M,1) = 2.D0/(DS(M)*DA(M))
      WK(M,2) = 2.D0/(DS(M)*DS(1))
      WK(M,3) = 2.D0/(DS(1)*DA(M))
C
C Store the upper triangle of D**T*D in AP in Linpack packed
C   storage format:  AP(k) = D**T*D(i,j) for j = 1 to m,
C   i = 1 to j.  Note that D2 is tridiagonal with nonzeros
C   in the upper right and lower left corners.  Hence,
C   D**T*D is pentadiagonal except for three nonzero
C   elements in the corners.
C
      W0 = A0*A0
      K = 0
C
C j = 1:
C
      K = K + 1
      AP(K) = WK(M,3)*WK(M,3)*DA(M) + WK(1,2)*WK(1,2)*
     .        DA(1) + WK(2,1)*WK(2,1)*DA(2) + W0*DA(1)
C
C j = 2:
C
      K = K + 1
      AP(K) = -WK(1,2)*WK(1,3)*DA(1) - WK(2,1)*WK(2,2)*DA(2)
      K = K + 1
      AP(K) = WK(1,3)*WK(1,3)*DA(1) + WK(2,2)*WK(2,2)*DA(2)
     .        + WK(3,1)*WK(3,1)*DA(3) + W0*DA(2)
C
C j = 3 to m-2:
C
      DO 3 J = 3,M-2
        DO 2 I = 1,J-3
          K = K + 1
          AP(K) = 0.D0
    2     CONTINUE
        K = K + 1
        AP(K) = WK(J-1,1)*WK(J-1,3)*DA(J-1)
        K = K + 1
        AP(K) = -WK(J-1,2)*WK(J-1,3)*DA(J-1) -
     .          WK(J,1)*WK(J,2)*DA(J)
        K = K + 1
        AP(K) = WK(J-1,3)*WK(J-1,3)*DA(J-1) +
     .          WK(J,2)*WK(J,2)*DA(J) +
     .          WK(J+1,1)*WK(J+1,1)*DA(J+1) + W0*DA(J)
    3   CONTINUE
C
C j = m-1:
C
      J = M - 1
      K = K + 1
      AP(K) = WK(M,1)*WK(M,3)*DA(M)
      DO 4 I = 2,J-3
        K = K + 1
        AP(K) = 0.D0
    4   CONTINUE
      K = K + 1
      AP(K) = WK(J-1,1)*WK(J-1,3)*DA(J-1)
      K = K + 1
      AP(K) = -WK(J-1,2)*WK(J-1,3)*DA(J-1) -
     .        WK(J,1)*WK(J,2)*DA(J)
      K = K + 1
      AP(K) = WK(J-1,3)*WK(J-1,3)*DA(J-1) +
     .        WK(J,2)*WK(J,2)*DA(J) +
     .        WK(J+1,1)*WK(J+1,1)*DA(J+1) + W0*DA(J)
C
C j = m:
C
      K = K + 1
      AP(K) = -WK(1,1)*WK(1,2)*DA(1) - WK(M,2)*WK(M,3)*DA(M)
      K = K + 1
      AP(K) = WK(1,1)*WK(1,3)*DA(1)
      DO 5 I = 3,M-3
        K = K + 1
        AP(K) = 0.D0
    5   CONTINUE
      K = K + 1
      AP(K) = WK(M-1,1)*WK(M-1,3)*DA(M-1)
      K = K + 1
      AP(K) = -WK(M-1,2)*WK(M-1,3)*DA(M-1) -
     .        WK(M,1)*WK(M,2)*DA(M)
      K = K + 1
      AP(K) = WK(M-1,3)*WK(M-1,3)*DA(M-1) +
     .        WK(M,2)*WK(M,2)*DA(M) +
     .        WK(1,1)*WK(1,1)*DA(1) + W0*DA(M)
C
C Overwrite AP with an R**T*R factorization (for upper
C   triangular matrix R) stored in packed format.
C
      CALL DPPCO (AP,M,RCOND,WK,IER)
      IF (IER .NE. 0) THEN
        RCOND = 0.D0
        IER = 2
        RETURN
      ENDIF
      IF (RCOND .LT. RCMIN) THEN
        IER = 3
        RETURN
      ENDIF
C
C Solve the three linear systems corresponding to the
C   x, y, and z components of G.
C
      DO 6 K = 1,M
        WK(K,1) = G(1,K)
    6   CONTINUE
      CALL DPPSL (AP,M,WK)
      DO 7 K = 1,M
        G(1,K) = WK(K,1)
    7   CONTINUE
      DO 8 K = 1,M
        WK(K,1) = G(2,K)
    8   CONTINUE
      CALL DPPSL (AP,M,WK)
      DO 9 K = 1,M
        G(2,K) = WK(K,1)
    9   CONTINUE
      DO 10 K = 1,M
        WK(K,1) = G(3,K)
   10   CONTINUE
      CALL DPPSL (AP,M,WK)
      DO 11 K = 1,M
        G(3,K) = WK(K,1)
   11   CONTINUE
      RETURN
C
C Invalid input parameter.
C
   21 IER = 1
      RETURN
      END