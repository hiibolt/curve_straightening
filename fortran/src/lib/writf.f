      SUBROUTINE WRITF (LOUT,M,F)
      INTEGER LOUT, M
      DOUBLE PRECISION F(3,0:M)
C
C***********************************************************
C
C                                                 From DMVC3
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/18/03
C
C   This subroutine creates an output data set consisting of
C the number of polygon vertices M+1 (format I4) followed by
C the sequence of vertices (ordered triples of coordinates):
C (f(i,j), i = 1,3), j = 0,m), three values per line (format
C 3E23.15).
C
C On input:
C
C       LOUT = Logical unit number for writing the solution.
C              0 .LE. LOUT .LE. 99.
C
C       M = Number of curve segments.  M > 0.
C
C       F = Array dimensioned 3 by 0:M containing the
C           vertices with the Cartesian coordinates of
C           f(j) in column j.
C
C Input parameters are not altered by this routine.
C
C An error message is written to the standard output unit
C if an input parameter is invalid or a write error is
C encountered writing to unit LOUT.
C
C Modules required by WRITF:  None
C
C***********************************************************
C
      INTEGER I, J
C
C Output formats:
C
  100 FORMAT (I4)
  110 FORMAT (3E23.15)
C
C Test for invalid parameters.
C
      IF (LOUT .LT. 0  .OR.  LOUT .GT. 99  .OR.  M .LE. 0)
     .  GO TO 11
C
C Create the data set.
C
      WRITE (LOUT,100,ERR=12) M+1
      WRITE (LOUT,110,ERR=12) ((F(I,J), I = 1,3), J = 0,M)
      RETURN
C
C Invalid input parameter.
C
   11 WRITE (*,210)
  210 FORMAT (///10X,'*** Error in WRITF:  invalid input ',
     .               'parameter ***')
      RETURN
C
C Error writing to unit LOUT.
C
   12 WRITE (*,220)
  220 FORMAT (///10X,'*** Error in WRITF writing to unit ',
     .               'LOUT ***')
      RETURN
      END
