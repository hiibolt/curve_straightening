      SUBROUTINE GRADL2 (M,W,DS,D1C,D2C,U,F,HU, G, IER)
      INTEGER M, IER
      DOUBLE PRECISION W, DS(M), D1C(3,M), D2C(3,M), U(3,M),
     .                 F(M), HU(3,M), G(3,M)
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
C   This subroutine computes the L2 gradient grad(phi) of
C the functional phi associated with a periodic geodesic
C curve approximation (Function PHI).
C
C Refer to the header comments and Function PHI for further
C details.
C
C On input:
C
C       M = Number of vertices and curve segments.  M >= 3.
C
C       W = Weight defining the penalty term.  W >= 0.
C
C       DS,D1C,D2C,U,F,HU = Arrays containing curve segment
C                           lengths, unit tangent vectors,
C                           etc. computed by Function PHI
C                           for the current curve c.
C
C The above parameters are not altered by this routine.
C
C       G = Array dimensioned 3 by M containing gradient
C           vectors grad(f)(c(i)).
C
C On output:
C
C       G = Components of the discretized L2-gradient
C           grad(phi) -- partial derivatives of phi with
C           respect to the vertices c(j).
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M or W is outside its valid range.
C             G is not altered if IER = 1.
C
C Modules required by GRADL2:  None
C
C Intrinsic functions called by GRADL2:  None
C
C***********************************************************
C
      DOUBLE PRECISION D1I, D1IP1, D2I, D2IM1, D2IP1, S, S1,
     .                 S2, SI, SIP1, UI, UIM1, UIP1
      INTEGER I
C
C Test for invalid input.
C
      IF (M .LT. 3  .OR.  W .LT. 0.D0) GO TO 11
C
C Overwrite G with grad(phi) = grad(E) + w*F'**T*F.
C
C
C I = 1:
C
      SI = 2.D0/DS(1)
      SIP1 = 2.D0/DS(2)
      S = SI + SIP1
      D1I = D1C(1,1)*(U(1,1)-U(1,M)) +
     .      D1C(2,1)*(U(2,1)-U(2,M)) +
     .      D1C(3,1)*(U(3,1)-U(3,M))
      D1IP1 = D1C(1,2)*(U(1,2)-U(1,1)) +
     .        D1C(2,2)*(U(2,2)-U(2,1)) +
     .        D1C(3,2)*(U(3,2)-U(3,1))
      D2IM1 = D2C(1,M)*U(1,M) + D2C(2,M)*U(2,M) +
     .        D2C(3,M)*U(3,M)
      D2I = D2C(1,1)*U(1,1) + D2C(2,1)*U(2,1) +
     .      D2C(3,1)*U(3,1)
      D2IP1 = D2C(1,2)*U(1,2) + D2C(2,2)*U(2,2) +
     .        D2C(3,2)*U(3,2)
      UIM1 = (U(1,M)**2 + U(2,M)**2 + U(3,M)**2)/2.D0
      UI = (U(1,1)**2 + U(2,1)**2 + U(3,1)**2)/2.D0
      UIP1 = (U(1,2)**2 + U(2,2)**2 + U(3,2)**2)/2.D0
      S1 = UIM1 + UI + SI*D1I - D2IM1 - D2I
      S2 = UI + UIP1 + SIP1*D1IP1 - D2I - D2IP1
      G(1,1) = SI*U(1,M) - S*U(1,1) + SIP1*U(1,2) +
     .         S1*D1C(1,1) - S2*D1C(1,2) + HU(1,1) +
     .         W*F(1)*G(1,1)
      G(2,1) = SI*U(2,M) - S*U(2,1) + SIP1*U(2,2) +
     .         S1*D1C(2,1) - S2*D1C(2,2) + HU(2,1) +
     .         W*F(1)*G(2,1)
      G(3,1) = SI*U(3,M) - S*U(3,1) + SIP1*U(3,2) +
     .         S1*D1C(3,1) - S2*D1C(3,2) + HU(3,1) +
     .         W*F(1)*G(3,1)
C
C I = 2 to M-1:
C
      DO 1 I = 2,M-1
        SI = SIP1
        SIP1 = 2.D0/DS(I+1)
        S = SI + SIP1
        D1IP1 = D1C(1,I+1)*(U(1,I+1)-U(1,I)) +
     .          D1C(2,I+1)*(U(2,I+1)-U(2,I)) +
     .          D1C(3,I+1)*(U(3,I+1)-U(3,I))
        D2I = D2IP1
        D2IP1 = D2C(1,I+1)*U(1,I+1) + D2C(2,I+1)*U(2,I+1) +
     .          D2C(3,I+1)*U(3,I+1)
        UI = UIP1
        UIP1 = (U(1,I+1)**2 + U(2,I+1)**2 +
     .          U(3,I+1)**2)/2.D0
        S1 = S2
        S2 = UI + UIP1 + SIP1*D1IP1 - D2I - D2IP1
        G(1,I) = SI*U(1,I-1) - S*U(1,I) + SIP1*U(1,I+1) +
     .           S1*D1C(1,I) - S2*D1C(1,I+1) + HU(1,I) +
     .           W*F(I)*G(1,I)
        G(2,I) = SI*U(2,I-1) - S*U(2,I) + SIP1*U(2,I+1) +
     .           S1*D1C(2,I) - S2*D1C(2,I+1) + HU(2,I) +
     .           W*F(I)*G(2,I)
        G(3,I) = SI*U(3,I-1) - S*U(3,I) + SIP1*U(3,I+1) +
     .           S1*D1C(3,I) - S2*D1C(3,I+1) + HU(3,I) +
     .           W*F(I)*G(3,I)
    1   CONTINUE
C
C I = M:
C
      SI = 2.D0/DS(M)
      SIP1 = 2.D0/DS(1)
      S = SI + SIP1
      D1IP1 = D1C(1,1)*(U(1,1)-U(1,M)) +
     .        D1C(2,1)*(U(2,1)-U(2,M)) +
     .        D1C(3,1)*(U(3,1)-U(3,M))
      D2I = D2IP1
      D2IP1 = D2C(1,1)*U(1,1) + D2C(2,1)*U(2,1) +
     .        D2C(3,1)*U(3,1)
      UI = UIP1
      UIP1 = (U(1,1)**2 + U(2,1)**2 + U(3,1)**2)/2.D0
      S1 = S2
      S2 = UI + UIP1 + SIP1*D1IP1 - D2I - D2IP1
      G(1,M) = SI*U(1,M-1) - S*U(1,M) + SIP1*U(1,1) +
     .         S1*D1C(1,M) - S2*D1C(1,1) + HU(1,M) +
     .         W*F(M)*G(1,M)
      G(2,M) = SI*U(2,M-1) - S*U(2,M) + SIP1*U(2,1) +
     .         S1*D1C(2,M) - S2*D1C(2,1) + HU(2,M) +
     .         W*F(M)*G(2,M)
      G(3,M) = SI*U(3,M-1) - S*U(3,M) + SIP1*U(3,1) +
     .         S1*D1C(3,M) - S2*D1C(3,1) + HU(3,M) +
     .         W*F(M)*G(3,M)
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
      END
