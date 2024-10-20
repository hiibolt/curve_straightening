      DOUBLE PRECISION FUNCTION PHI (M,C,W, DS,D1C,D2C,U,F,
     .                               G,HU,CL,E,FNS,IER)
      INTEGER M, IER, COUNT
      DOUBLE PRECISION C(3,0:M), W, DS(M), D1C(3,M),
     .                 D2C(3,M), U(3,M), F(M), G(3,M),
     .                 HU(3,M), CL, E, FNS
C
C***********************************************************
C
C                                                  From GEO2
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/19/03
C
C   Given a polygonal curve defined by a sequence of verti-
C ces c(i), i = 0 to m, along with a weight w associated
C with the constraint penalty, this function returns the
C following:
C
C       ds(i) = |c(i)-c(i-1)|,
C
C       D1c(i) = [c(i)-c(i-1)]/ds(i),
C
C       D2c(i) = [D1c(i+1)-D1c(i)]/([ds(i)+ds(i+1)]/2)
C
C for i = 1 to m, where ds(m+1) = ds(1) and D1c(m+1) =
C D1c(1), and the total geodesic curvature
C
C       E(c) = 0.5*Sum[ (|u(i)|**2)*da(i) ]
C
C where
C
C       da(i) = [ds(i)+ds(i+1)],
C
C       u(i) = [I - fn(i)*fn(i)**T]*D2c(i)
C            = D2c(i) - <fn(i),D2c(i)>*fn(i),
C
C       fn(i) = grad(f)(c(i))/|grad(f)(c(i))|,
C
C and Sum[_] denotes the sum over i = 1 to m.  This function
C also returns
C
C       FNS = 0.5*Sum[ F(i)**2 ]
C
C and
C
C       phi(c) = E(c) + w*FNS
C
C where F(i) = f(c(i)).
C
C Refer to the comments at the beginning of the main program
C for more details.
C
C   Subroutine FEVAL is called to compute F(i) = f(c(i)) and
C the gradients of f at c(i) for i = 1 to m.
C
C On input:
C
C       M = Number of vertices and curve segments.  M >= 3.
C
C       C = Array dimensioned 3 by 0:M containing the
C           vertices with the Cartesian coordinates of
C           c(j) in column j.  Adjacent vertices c(i-1)
C           and c(i) must be distinct for i = 1 to m.
C           It is assumed that c(0) = c(m).
C
C       W = Weight defining the penalty term.  W >= 0.
C
C The above parameters are not altered by this function.
C
C       DS = Array of length at least M.
C
C       D1C,D2C,U = Arrays dimensioned 3 by M.
C
C       F = Array of length at least M.
C
C       G,HU = Arrays dimensioned 3 by M.
C
C On output:
C
C       DS,D1C,D2C,U = Arrays containing the terms defined
C                      above:  arc lengths, unit tangent
C                      vectors, curvature vectors, and cur-
C                      vature vectors projected onto tangent
C                      planes.
C
C       F,G = Arrays containing function values f(c(i)) and
C             gradient vectors grad(f)(c(i)), respectively,
C             computed by Subroutine FEVAL.
C
C       HU = Array containing terms needed for the gradient
C            of phi:  s*H(i)*u(i), for i = 1 to m, where
C            s = -da(i)*<fn(i),D2c(i)>/|G( ,i)| and H(i) is
C            the Hessian of f at c(i).
C
C       CL = Total curve length:  Sum[ ds(i) ].
C
C       E = Value of the functional representing total
C           geodesic curvature.
C
C       FNS = Half the sum of squared components of F.
C
C       PHI = E(c) + w*FNS.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M or W is outside its valid range
C                     on input.  Output parameters are not
C                     altered in this case.
C             IER = 2 if c(i-1) = c(i) for some i (1 to m).
C             IER = 3 if grad(f)(c(i)) = 0 for some i
C                     (1 to m).
C
C Module required by PHI:  FEVAL
C
C Intrinsic function called by PHI:  SQRT
C
C***********************************************************
C
      DOUBLE PRECISION CSUM, DA, F11, F12, F13, F22, F23,
     .                 F33, FSUM, S
      INTEGER I, ND
      DATA ND/2/
      IF (M .LT. 3  .OR.  W .LT. 0.D0) GO TO 11
C
C Compute ds(i) and D1c(i) for i = 1 to m.  Sum[ds(i)] is
C   accumulated in CL.
C
      CL = 0.D0
      DO 1 I = 1,M
        DS(I) = (C(1,I)-C(1,I-1))**2 +
     .          (C(2,I)-C(2,I-1))**2 +
     .          (C(3,I)-C(3,I-1))**2
        IF (DS(I) .EQ. 0.D0) GO TO 12
        DS(I) = SQRT(DS(I))
        CL = CL + DS(I)
        D1C(1,I) = (C(1,I)-C(1,I-1))/DS(I)
        D1C(2,I) = (C(2,I)-C(2,I-1))/DS(I)
        D1C(3,I) = (C(3,I)-C(3,I-1))/DS(I)
    1   CONTINUE
C
C Compute D2c(i),
C         F(i) = f(c(i)),
C         G( ,i) = grad(f)(c(i)),
C         u(i) = D2c(i) - <fn(i),D2c(i)>*fn(i), and
C         HU( ,i) = s*H(i)*u(i),
C
C for i = 1 to m-1, where s = -da(i)*<fn(i),D2c(i)>/|G( ,i)|
C and H(i) is the Hessian of f at c(i).
C
C Sum[ (|u(i)|**2)*da(i) ] is accumulated in CSUM, and
C Sum[ F(i)**2 ] is accumulated in FSUM.
C
      CSUM = 0.D0
      FSUM = 0.D0
      COUNT = 0
      DO 2 I = 1,M-1
        DA = DS(I) + DS(I+1)
        S = 2.D0/DA
        D2C(1,I) = S*(D1C(1,I+1)-D1C(1,I))
        D2C(2,I) = S*(D1C(2,I+1)-D1C(2,I))
        D2C(3,I) = S*(D1C(3,I+1)-D1C(3,I))
        CALL FEVAL (C(1,I),C(2,I),C(3,I),ND, F(I),G(1,I),
     .              G(2,I),G(3,I),F11,F12,F13,F22,F23,F33)
        FSUM = FSUM + F(I)**2
        S = G(1,I)**2 + G(2,I)**2 + G(3,I)**2
        IF (S .EQ. 0.D0) GO TO 13
        S = (G(1,I)*D2C(1,I) + G(2,I)*D2C(2,I) +
     .       G(3,I)*D2C(3,I))/S
        U(1,I) = D2C(1,I) - S*G(1,I)
        U(2,I) = D2C(2,I) - S*G(2,I)
        U(3,I) = D2C(3,I) - S*G(3,I)
        CSUM = CSUM + DA*(U(1,I)**2 + U(2,I)**2 + U(3,I)**2)
        COUNT = COUNT + 1

        S = -DA*S
        HU(1,I) = S*(F11*U(1,I) + F12*U(2,I) + F13*U(3,I))
        HU(2,I) = S*(F12*U(1,I) + F22*U(2,I) + F23*U(3,I))
        HU(3,I) = S*(F13*U(1,I) + F23*U(2,I) + F33*U(3,I))
    2   CONTINUE
C
C Compute D2c(m), F(m), G( ,m), u(m), and HU( ,m), and add
C   the final contribution to FSUM and CSUM.
C
      DA = DS(M) + DS(1)
      S = 2.D0/DA
      D2C(1,M) = S*(D1C(1,1)-D1C(1,M))
      D2C(2,M) = S*(D1C(2,1)-D1C(2,M))
      D2C(3,M) = S*(D1C(3,1)-D1C(3,M))

      CALL FEVAL (C(1,M),C(2,M),C(3,M),ND, F(M),G(1,M),
     .            G(2,M),G(3,M),F11,F12,F13,F22,F23,F33)
      FSUM = FSUM + F(M)**2

      S = G(1,M)**2 + G(2,M)**2 + G(3,M)**2
      IF (S .EQ. 0.D0) GO TO 13
      S = (G(1,M)*D2C(1,M) + G(2,M)*D2C(2,M) +
     .     G(3,M)*D2C(3,M))/S
      U(1,M) = D2C(1,M) - S*G(1,M)
      U(2,M) = D2C(2,M) - S*G(2,M)
      U(3,M) = D2C(3,M) - S*G(3,M)
      CSUM = CSUM + DA*(U(1,M)**2 + U(2,M)**2 + U(3,M)**2)

      S = -DA*S
      HU(1,M) = S*(F11*U(1,M) + F12*U(2,M) + F13*U(3,M))
      HU(2,M) = S*(F12*U(1,M) + F22*U(2,M) + F23*U(3,M))
      HU(3,M) = S*(F13*U(1,M) + F23*U(2,M) + F33*U(3,M))
C
C Compute E, FNS, and PHI.
C
      E = 0.5D0*CSUM
      FNS = 0.5D0*FSUM
      PHI = E + W*FNS
C
C No error encountered.
C
      IER = 0
      RETURN
C
C M < 3 or W < 0.
C
   11 IER = 1
      RETURN
C
C DS(i) = 0 for 1 <= i <= m.
C
   12 IER = 2
      RETURN
C
C grad(f)(c(i)) = 0 for 1 <= i <= m.
C
   13 IER = 3
      RETURN
      END
