      DOUBLE PRECISION FUNCTION DSTORE (X)
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                                 From GLRH2
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/25/96
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C On input:
C
C       X = Double precision value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       DSTORE = Value of X after it has been stored and
C                possibly truncated or rounded to the double
C                precision word length.
C
C Modules required by DSTORE:  None
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
      Y = X
      DSTORE = Y
      RETURN
      END