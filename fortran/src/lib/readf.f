      SUBROUTINE READF (LIN,M,M1, C,C1,IER)
      INTEGER LIN, M, M1, IER
      DOUBLE PRECISION C(3,0:M), C1(3,0:M1)
C
C***********************************************************
C
C                                                  From GEO1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/08/03
C
C   This subroutine reads a data set representing a closed
C curve c1 stored in the format used by WRITF:  the number
C of polygon vertices m0+1 (format I4) followed by the se-
C quence of vertices (ordered triples of coordinates):
C (c1(i,j), i = 1,3), j = 0,m0), three values per line
C (format 3E23.15).  A polygonal curve c consisting of m+1
C vertices is then computed by taking m+1 points uniformly
C distributed on c1, with c(0) = c1(0) and c(m) = c1(m).
C Vertices c1(0) and c1(m) should be identical, and adjacent
C vertices should be distinct, but the data set is not
C tested for validity.
C
C On input:
C
C       LIN = Logical unit number for input.  0 .LE. LIN
C             .LE. 99.
C
C       M = Number of curve segments in c.  M >= 3.
C
C       M1 = Maximum number of curve segments M0 in c1.
C            M1 >= 3.
C
C The above parameters are not altered by this function.
C
C       C = Array dimensioned 3 by 0:M.
C
C       C1 = Array dimensioned 3 by 0:M1.
C
C On output:
C
C       C = Sequence of M+1 vertex coordinates unless
C           IER > 0.
C
C       C1 = Sequence of M0+1 input vertex coordinates
C            unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LIN, M, or M1 is outside its valid
C                     range on input.
C             IER = 2 if a read error occurred.
C             IER = 3 if M0 > M1, where M0 is the value read
C                     from the data set.
C
C Modules required by READF:  None
C
C Intrinsic function called by READF:  DBLE
C
C***********************************************************
C
      DOUBLE PRECISION B1, B2, R, R0
      INTEGER I, J, J0, M0
C
C Input formats:
C
  100 FORMAT (I4)
  110 FORMAT (3E23.15)
C
C Test for invalid parameters.
C
      IF (LIN .LT. 0  .OR.  LIN .GT. 99  .OR.  M .LT. 3
     .    .OR.  M1 .LT. 3) GO TO 11
C
C Read the data set:  M0 and C1.
C
      READ (LIN,100,ERR=12) M0
      IF (M0 .GT. M1) GO TO 13
      DO 1 J = 0,M0
        READ (LIN,110,ERR=12) (C1(I,J), I = 1,3)
    1   CONTINUE
C
C Compute M+1 uniformly distributed vertices c on c1.
C
      R = DBLE(M0)/DBLE(M)
      C(1,0) = C1(1,0)
      C(2,0) = C1(2,0)
      C(3,0) = C1(3,0)
      DO 2 J = 1,M
        R0 = DBLE(J)*R
        J0 = R0
        IF (J0 .EQ. M0) J0 = J0 - 1
        B2 = R0 - DBLE(J0)
        B1 = 1.D0 - B2
        C(1,J) = B1*C1(1,J0) + B2*C1(1,J0+1)
        C(2,J) = B1*C1(2,J0) + B2*C1(2,J0+1)
        C(3,J) = B1*C1(3,J0) + B2*C1(3,J0+1)
    2   CONTINUE
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
C
C Read error.
C
   12 IER = 2
      RETURN
C
C Invalid value of M0.
C
   13 IER = 3
      RETURN
      END
