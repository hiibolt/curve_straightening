      SUBROUTINE FEVAL (X,Y,Z,ND, F,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,FZZ)
      INTEGER ND
      DOUBLE PRECISION X, Y, Z, F, FX, FY, FZ, FXX, FXY,
     .                 FXZ, FYY, FYZ, FZZ
C
C***********************************************************
C
C                                              From GEODESIC
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/19/03
C
C   This subroutine returns the value and, optionally, first
C and second partial derivatives of a trivariate function f.
C
C On input:
C
C       X,Y,Z = Evaluation point.
C
C       ND = Order of derivatives to be returned:
C            ND = 0 if only F is to be computed.
C            ND = 1 if only F, FX, FY, and FZ are to be
C                   computed.
C            ND = 2 if F, first partials and second
C                   partials are to be computed.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Function value f(x,y,z).
C
C       FX,FY,FZ = Components of the gradient of f at
C                  (x,y,z) if ND >= 1.
C
C       FXX,...,FZZ = Components of the upper triangle of
C                     the Hessian of f at (x,y,z) if ND > 1.
C
C***********************************************************
C
      DOUBLE PRECISION AX, AY, AZ, EPS
C
C Ellipsoid:  x**2/Ax + y**2/Ay + z**2/Az - 1 = 0.
C
C The Gaussian curvature is K in [Kmin,Kmax], where
C Kmin = Az/(Ax*Ay) and Kmax = Ax/(Ay*Az) so that
C Kmin/Kmax = (Az/Ax)**2.
C
C If Az < Ay < Ax and Kmin/Kmax >= 1/16, then there are no
C simple closed (periodic) geodesics other than the principal
C ellipses (x = 0, y = 0, and z = 0).  For the
C following parameters, there are no non-standard simple
C closed geodesics if EPS > 0.
C
C
      F = 2*(X**2)+3*(Y**2)+5*(Z**2)+(X**3)*(Z)+(Y)*(Z**3)-1.D0
      IF (ND .LT. 1) RETURN
C
C Gradient of f.
C
      FX = 4.D0*X + 3.D0*(X**2)*Z
      FY = 6.D0*Y + (Z**3)
      FZ = 10.D0*Z + (X**3) + 3.D0*(Y)*(Z**2)

      IF (ND .LT. 2) RETURN
C
C Hessian of f.
C
      FXX = 4.D0 + 6.D0*X*Z
      FXY = 0.D0
      FXZ = 3.D0*(X**2)
      FYY = 6.D0
      FYZ = 3.D0*(Z**2)
      FZZ = 10.D0 + 6.D0*(Y)*(Z)
      RETURN
      END
