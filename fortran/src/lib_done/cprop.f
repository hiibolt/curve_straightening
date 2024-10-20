      SUBROUTINE CPROP (M,DS,U, S,CRV,CNRM2,CMAX,IER)
      INTEGER M, IER
      DOUBLE PRECISION DS(M), U(3,M), S(M), CRV(M), CNRM2,
     .                 CMAX
C
C***********************************************************
C
C                                                  From GEO1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/19/03
C
C   Given a polygonal approximation to a closed curve c in
C a level set {(x,y,z):  f(x,y,z) = 0}, along with segment
C lengths ds(i) and discrete approximations u(i) to vertex
C values of geodesic curvature Kg, this subroutine computes
C vertex values |u(i)| and norms of the total geodesic
C curvature.
C
C On input:
C
C       M = Number of vertices and polygonal curve segments
C           in the approximation of c.  M >= 1.
C
C       DS = Array of length M containing segment lengths:
C            ds(i) = |c(i)-c(i-1)| for i = 1 to m.
C
C       U = Array dimensioned 3 by M containing vertex
C           values of [I - fn(i)*fn(i)**T]*D2c(i), where
C           fn(i) is the unit normal to the surface and
C           D2c(i) is an approximation to c'' at c(i).
C
C The above parameters are not altered by this routine.
C
C       S,CRV = Arrays of length at least M.
C
C On output:
C
C       S = Cumulative segment lengths:  s(j) = Sum[ds(i)],
C           where the sum is over i = 1 to j for j = 1 to m.
C           For the natural parameterization by arc length,
C           s(j) is the parameter value associated with ver-
C           tex j, s(j)-ds(j)/2 is the parameter value
C           associated with midpoint j, and s(m) is the
C           total curve length.
C
C       CRV = Vertex values of the geodesic curvature
C             |Kg(s)| = |u(i)|.
C
C       CNRM2 = Approximation to the integral with respect
C               to arc length s of Kg(s)**2.
C
C       CMAX = Approximation to the max-norm of |Kg(s)|.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M < 1.
C
C Modules required by CPROP:  None
C
C Intrinsic functions called by CPROP:  ABS, MAX, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION CJ
      INTEGER J
C
C Test for error 1:
C
      IF (M .LT. 1) THEN
        IER = 1
        RETURN
      ENDIF
      IER = 0
C
C Compute cumulative arc length values s(j), j = 1 to m.
C
      S(1) = DS(1)
      DO 1 J = 2,M
        S(J) = S(J-1) + DS(J)
    1   CONTINUE
C
C Compute vertex geodesic curvature values Kg(s).
C
      DO 2 J = 1,M
        CRV(J) = SQRT(U(1,J)**2 + U(2,J)**2 + U(3,J)**2)
    2   CONTINUE
C
C Accumulate norms using the composite midpoint rule for
C   CNRM2.  CJ is the midpoint value of Kg(s).
C
      CNRM2 = 0.D0
      CMAX = 0.D0
      CJ = (CRV(M) + CRV(1))/2.D0
      DO 3 J = 1,M
        CNRM2 = CNRM2 + CJ*CJ*DS(J)
        CMAX = MAX(CMAX,ABS(CJ))
        IF (J .LT. M) CJ = (CRV(J) + CRV(J+1))/2.D0
    3   CONTINUE
      RETURN
      END
