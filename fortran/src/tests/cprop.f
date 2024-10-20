      PROGRAM TEST_CPROP
      INTEGER M, IER
      DOUBLE PRECISION CNRM2, CMAX
      PARAMETER (M=4) ! Set the number of vertices to 4 for this test case
      DOUBLE PRECISION DS(M), U(3, M), S(M), CRV(M)

C Initialize input values for DS (segment lengths)
      DATA DS / 1.0D0, 1.5D0, 2.0D0, 2.5D0 /

C Initialize input values for U (3 by M matrix)
      DATA U / 0.1D0, 0.2D0, 0.3D0,
     .         0.4D0, 0.5D0, 0.6D0,
     .         0.7D0, 0.8D0, 0.9D0,
     .         1.0D0, 1.1D0, 1.2D0 /

C Call the CPROP subroutine
      CALL CPROP(M, DS, U, S, CRV, CNRM2, CMAX, IER)

C Output the results
      PRINT *, 'IER = ', IER
      PRINT *, 'S = '
      DO 10 I = 1, M
            PRINT *, S(I)
   10 CONTINUE

      PRINT *, 'CRV = '
      DO 20 I = 1, M
            PRINT *, CRV(I)
   20 CONTINUE

      PRINT *, 'CNRM2 = ', CNRM2
      PRINT *, 'CMAX = ', CMAX

      END
