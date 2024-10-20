      PROGRAM TEST_PHI
      IMPLICIT NONE
      
      INTEGER M, IER
      INTEGER I, J, LPRT
      DOUBLE PRECISION W,CL,E,FNS,TS,PHI
      DOUBLE PRECISION C(3,0:400),DS(400),D1C(3,400),D2C(3,400), 
     .                 U(3,400),F(400),G(3,400),HU(3,400)

      DOUBLE PRECISION FUNCTION PHI
      
      ! Initialize test data
      M = 400
      W = 0

      CALL INITC (M,C)

      WRITE (LPRT,122) ((C(I,J), I = 1,3), J = 0,M)
C Prints with each vertex on a separate line.
  122 FORMAT (5X,':33 Vertices of initial curve:  C(j) = ',
     .        3(1P,D35.15))
      
      ! Call PHI function
      TS = PHI(M,C,W, DS,D1C,D2C,U,F,G,HU,CL,E,FNS,IER)
      
      ! Check for errors
      IF (IER .NE. 0) THEN
         PRINT *, 'Error occurred. IER = ', IER
         STOP
      END IF
      
      ! Print results
      PRINT *, 'Curve Length (CL) = ', CL
      PRINT *, 'Total Geodesic Curvature (E) = ', E
      PRINT *, 'Half Sum of Squared Function Values (FNS) = ', FNS
      PRINT *, 'PHI = ', TS

      END
