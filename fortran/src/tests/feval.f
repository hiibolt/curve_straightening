      PROGRAM TEST_FEVAL
      DOUBLE PRECISION X,Y,Z,F,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,FZZ
      INTEGER ND

      EXTERNAL FEVAL

C Test Case 1: ND = 0 (Only the function value is calculated)
      X = 1.0D0
      Y = 1.0D0
      Z = 1.0D0
      ND = 0
      CALL FEVAL(X,Y,Z,ND, F,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,FZZ)
      PRINT *, 'Test Case 1 (ND = 0):'
      PRINT *, 'F = ', F

C Test Case 2: ND = 1 (First derivatives are also calculated)
      X = 1.0D0
      Y = 1.0D0
      Z = 1.0D0
      ND = 1
      CALL FEVAL(X,Y,Z,ND, F,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,FZZ)
      PRINT *, 'Test Case 2 (ND = 1):'
      PRINT *, 'F  = ', F
      PRINT *, 'FX = ', FX, ' FY = ', FY, ' FZ = ', FZ

C Test Case 3: ND = 2 (Second derivatives are also calculated)
      X = 1.0D0
      Y = 1.0D0
      Z = 1.0D0
      ND = 2
      CALL FEVAL(X,Y,Z,ND, F,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,FZZ)
      PRINT *, 'Test Case 3 (ND = 2):'
      PRINT *, 'F   = ', F
      PRINT *, 'FX  = ', FX, ' FY  = ', FY, ' FZ  = ', FZ
      PRINT *, 'FXX = ', FXX, ' FXY = ', FXY, ' FXZ = ', FXZ
      PRINT *, 'FYY = ', FYY, ' FYZ = ', FYZ, ' FZZ = ', FZZ

      STOP
      END
