      SUBROUTINE INITC(M, C)
      INTEGER M
      DOUBLE PRECISION C(3,0:M)
      DOUBLE PRECISION THETA, PHI, PI
      PARAMETER (PI = 3.141592653589793D0)
      INTEGER J
      DOUBLE PRECISION AX, AY, AZ, EPS

C Initialize ellipsoid semi-axes values (same as in FEVAL)
      DATA EPS/0.09D0/
      AX = 1.D0
      AY = 1.D0 - EPS
      AZ = 0.25D0

C Distribute M points on the surface of the ellipsoid
C using spherical coordinates and map to Cartesian coordinates
      DO J = 1, M
         THETA = 2.0D0 * PI * (J-1) / M
         PHI = PI * (J-1) / M

C Cartesian coordinates of points on the ellipsoid
         C(1,J) = AX * COS(THETA) * SIN(PHI)
         C(2,J) = AY * SIN(THETA) * SIN(PHI)
         C(3,J) = AZ * COS(PHI)
      END DO

C Copy the first point to the last column (C(:,0) = C(:,M))
      C(1,0) = C(1,M)
      C(2,0) = C(2,M)
      C(3,0) = C(3,M)

      RETURN
      END
