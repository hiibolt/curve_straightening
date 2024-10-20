      SUBROUTINE PLTCRV (LUN,N,IND,M,F,I1,I2,PLTSIZ, IER)
      INTEGER LUN, N, IND(0:N), M, I1, I2, IER
      DOUBLE PRECISION F(3,0:M), PLTSIZ
C
C***********************************************************
C
C                                                 From DMVC3
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   12/04/03
C
C   Given a sequence of 3-D vertices f(j), j = 0 to m, this
C subroutine creates a file containing a level-2 Encapsula-
C ted PostScript (EPS) plot of the polygonal curve
C corresponding to the orthogonal projection of the sequence
C onto one of the Cartesian coordinate planes.  A subset of
C the projected vertices may be distinguished by marker
C symbols centered at the points.
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       N = k-1, where k is the number of vertices to be
C           marked.  -1 <= N <= M.
C
C       IND = Integer array of length N+1 containing a
C             strictly increasing sequence of indexes
C             (0 to M) of the vertices to be marked.
C
C       M = Number of polygonal curve segments.  M > 0.
C
C       F = Array dimensioned 3 by 0:M containing the
C           vertices with the Cartesian coordinates of
C           f(j) in column j.
C
C       I1,I2 = Row indexes for F (distinct values in the
C               range 1 to 3) specifying a pair of Cartesian
C               coordinates and defining the projection
C               plane as the I1-I2 plane.  I1 and I2 are
C               the horizontal and vertical components,
C               respectively.
C
C       PLTSIZ = Plot size in inches.  The view volume
C                (bounding box defined by the projected
C                curve) is mapped to a rectangular viewport
C                with maximum side-length equal to .88*
C                PLTSIZ.  The containing rectangle is cen-
C                tered on the 8.5 by 11 inch page, and its
C                boundary is drawn.  Labels below, to the
C                left of, and to the right of the rectangle
C                extend its dimensions by 1/4 inch, 2/3
C                inches, and 1/6 inch, respectively, for
C                FSIZ = 12 pts.  1.0 <= PLTSIZ <= 6.5.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, N, M, I1, I2, or PLTSIZ is
C                     outside its valid range.
C             IER = 2 if the elements of IND are invalid or
C                     not strictly increasing.
C             IER = 3 if all vertices are identical.
C             IER = 4 if an error was encountered in writing
C                     to unit LUN.
C
C   Various plotting options can be controlled by altering
C the data statement below.
C
C Modules required by PLTCRV:  None
C
C Intrinsic functions called by PLTCRV:  CHAR, MAX, MIN,
C                                        NINT, REAL
C
C***********************************************************
C
      INTEGER I, IH, IPX1, IPX2, IPY1, IPY2, IW, J
      LOGICAL ASPECT, AXES, DOTS, SQUARE
      REAL    BSIZ, DSIZ, DX, DY, FSIZ, R, SFX, SFY, T, TX,
     .        TY, WX1, WX2, WY1, WY2, XM, YM
C
      DATA    ASPECT/.TRUE./, AXES/.TRUE./,  BSIZ/6.0/,
     .        DOTS/.FALSE./,  DSIZ/3.0/,     FSIZ/12.0/,
     .        SQUARE/.FALSE./
C
C Local parameters:
C
C ASPECT =    Logical variable with value TRUE if the
C               aspect ratio of the view volume is to be
C               preserved in the viewport, FALSE if the
C               viewport is to be square.  Note the limit-
C               ation on the size of the aspect ratio R
C AXES =      Logical variable with value TRUE if the x and
C               y axes (intersected with the view volume)
C               are to be drawn
C BSIZ =      Box size in points for marker symbols
C DOTS =      Logical variable with value TRUE if the curve
C             is to be rendered as a set of small dots at
C             the vertices instead of a sequence of line
C             segments connecting the vertices
C DSIZ =      Dot size in points used if DOTS = TRUE
C DX =        View volume width WX2-WX1
C DY =        View volume height WY2-WY1
C FSIZ =      Font size in points for labels
C I,J =       Array indexes
C IH =        Height of the viewport in points
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the containing rectangle
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the containing rectangle
C IW =        Width of the viewport in points
C R =         Aspect ratio of the viewport:
C               Max(1/2,Min(2,DX/DY)) if ASPECT = TRUE,
C               1 otherwise
C SFX,SFY =   Scale factors for mapping world coordinates
C               (window coordinates in [WX1,WX2] X [WY1,WY2])
C               to viewport coordinates
C SQUARE =    Logical variable with value TRUE if markers
C               are to be square boxes, FALSE if markers are
C               to be X's
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C WX1,WY1 =   X and y (world) coordinates of the lower left
C               corner of the window (view volume)
C WX2,WY2 =   X and y (world) coordinates of the upper right
C               corner of the window (view volume)
C XM,YM =     Location in world coordinates of the center of
C               a marker symbol or beginning of a label
C
C
C Test for error 1.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.  N .LT. -1
     .    .OR.  N .GT. M  .OR.  M .LE. 0  .OR.  I1 .LT. 1
     .    .OR.  I1 .GT. 3  .OR.  I2 .LT. 1  .OR.  I2 .GT. 3
     .    .OR.  I1 .EQ. I2  .OR.  PLTSIZ .LT. 1.D0  .OR.
     .    PLTSIZ .GT. 6.5D0) GO TO 11
C
C Test for error 2.
C
      IF (N .GE. 0) THEN
        IF (IND(0) .LT. 0  .OR.  IND(0) .GT. M) GO TO 12
        DO 1 I = 1,N
          IF (IND(I) .LT. 0  .OR.  IND(I) .GT. M) GO TO 12
          IF (IND(I) .LE. IND(I-1)) GO TO 12
    1     CONTINUE
      ENDIF
C
C Compute the window (view volume) corner coordinates
C   (WX1,WY1) and WX2,WY2).
C
      WX1 = REAL(F(I1,0))
      WX2 = WX1
      WY1 = REAL(F(I2,0))
      WY2 = WY1
      DO 2 J = 1,M
        WX1 = MIN(WX1,REAL(F(I1,J)))
        WX2 = MAX(WX2,REAL(F(I1,J)))
        WY1 = MIN(WY1,REAL(F(I2,J)))
        WY2 = MAX(WY2,REAL(F(I2,J)))
    2   CONTINUE
C
C Compute the dimensions and aspect ratio of the view
C   volume.
C
      DX = WX2 - WX1
      DY = WY2 - WY1
      IF (DX .LT. 1.E-6  .AND.  DY .LT. 1.E-6) GO TO 13
      IF (DY .LT. 1.E-6) THEN
        WY1 = WY1 - DX/4.0
        WY2 = WY2 + DX/4.0
        DY = WY2 - WY1
      ELSEIF (DX .LT. 1.E-6) THEN
        WX1 = WX1 - DY/4.0
        WX2 = WX2 + DY/4.0
        DX = WX2 - WX1
      ENDIF
      R = 1.0
      IF (ASPECT) R = MAX(0.5,MIN(2.0,DX/DY))
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the containing
C   rectangle.  The coordinates, specified in default user
C   space units (points, at 72 points/inch with origin at
C   the lower left corner of the page), are chosen to have
C   aspect ratio R, and to center the plot on the 8.5 by 11
C   inch page.  The center of the page is (306,396), and
C   T = PLTSIZ/2 in points.
C
      T = 36.0*PLTSIZ
      IF (R .GE. 1.0) THEN
        IPX1 = 306 - NINT(T)
        IPX2 = 306 + NINT(T)
        IPY1 = 396 - NINT(T/R)
        IPY2 = 396 + NINT(T/R)
      ELSE
        IPX1 = 306 - NINT(T*R)
        IPX2 = 306 + NINT(T*R)
        IPY1 = 396 - NINT(T)
        IPY2 = 396 + NINT(T)
      ENDIF
C
C Output header comments.  The bounding box corner coordi-
C   nates are obtained by extending the containing rectangle
C   to include its boundary and the labels.
C
      WRITE (LUN,100,ERR=14) IPX1-NINT(4.0*FSIZ),
     .                       IPY1-NINT(1.5*FSIZ),
     .                       IPX2+NINT(FSIZ), IPY2+1
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Minimum Variation Curve'/
     .        '%%Creator:  DMVC3'/
     .        '%%EndComments')
C
C Set the line thickness to 2 points, and draw the boundary
C   of the containing rectangle.
C
      T = 2.0
      WRITE (LUN,110,ERR=14) T
      WRITE (LUN,120,ERR=14) IPX1, IPY1
      WRITE (LUN,130,ERR=14) IPX1, IPY2
      WRITE (LUN,130,ERR=14) IPX2, IPY2
      WRITE (LUN,130,ERR=14) IPX2, IPY1
      WRITE (LUN,140,ERR=14)
      WRITE (LUN,150,ERR=14)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT (2I4,' moveto')
  130 FORMAT (2I4,' lineto')
  140 FORMAT ('closepath')
  150 FORMAT ('stroke')
C
C Set IW and IH to the width and height of a viewport
C   obtained by shrinking the containing rectangle by 12%
C   in each direction.
C
      IW = NINT(0.88*REAL(IPX2-IPX1))
      IH = NINT(0.88*REAL(IPY2-IPY1))
C
C Set up a viewport mapping.
C
      SFX = REAL(IW)/DX
      SFY = REAL(IH)/DY
      TX = REAL(306-IW/2) - SFX*WX1
      TY = REAL(396-IH/2) - SFY*WY1
      IF (DOTS) THEN
C
C Set the line thickness to DSIZ, and use round line caps.
C
        WRITE (LUN,110,ERR=14) DSIZ
        WRITE (LUN,160,ERR=14)
  160   FORMAT ('1 setlinecap')
C
C Render the curve as a set of dots.
C
        DO 3 J = 0,M
          WRITE (LUN,170,ERR=14) SFX*REAL(F(I1,J))+TX,
     .                           SFY*REAL(F(I2,J))+TY
          WRITE (LUN,180,ERR=14) SFX*REAL(F(I1,J))+TX,
     .                           SFY*REAL(F(I2,J))+TY
    3     CONTINUE
C
C Paint the path and reset the line thickness to 1 point.
C
        WRITE (LUN,150,ERR=14)
        T = 1.0
        WRITE (LUN,110,ERR=14) T
      ELSE
C
C Set the line thickness to 1 point.
C
        T = 1.0
        WRITE (LUN,110,ERR=14) T
C
C Render the curve as a sequence of line segments.
C
        WRITE (LUN,170,ERR=14) SFX*REAL(F(I1,0))+TX,
     .                         SFY*REAL(F(I2,0))+TY
        DO 4 J = 1,M
          WRITE (LUN,180,ERR=14) SFX*REAL(F(I1,J))+TX,
     .                           SFY*REAL(F(I2,J))+TY
    4   CONTINUE
C
C Paint the path.
C
        WRITE (LUN,150,ERR=14)
      ENDIF
  170 FORMAT (2F12.6,' moveto')
  180 FORMAT (2F12.6,' lineto')
      IF (AXES) THEN
C
C Draw x and y axes (if contained in the view volume).
C
        IF (WX1 .LT. 0.0  .AND.  WX2 .GT. 0.0) THEN
          WRITE (LUN,170,ERR=14) TX, SFY*WY1+TY
          WRITE (LUN,180,ERR=14) TX, SFY*WY2+TY
        ENDIF
        IF (WY1 .LT. 0.0  .AND.  WY2 .GT. 0.0) THEN
          WRITE (LUN,170,ERR=14) SFX*WX1+TX, TY
          WRITE (LUN,180,ERR=14) SFX*WX2+TX, TY
        ENDIF
      ENDIF
C
C Draw marker symbols (squares or X's) at the vertices to
C   be distinguished.
C
C   T = BSIZ/2,
C   (XM,YM) = Location of symbol center.
C
      T = BSIZ/2.0
      DO 5 I = 0,N
        J = IND(I)
        XM = SFX*REAL(F(I1,J))+TX
        YM = SFY*REAL(F(I2,J))+TY
        WRITE (LUN,170,ERR=14) XM-T, YM-T
        IF (SQUARE) THEN
          WRITE (LUN,180,ERR=14) XM+T, YM-T
          WRITE (LUN,180,ERR=14) XM+T, YM+T
          WRITE (LUN,180,ERR=14) XM-T, YM+T
          WRITE (LUN,140,ERR=14)
        ELSE
          WRITE (LUN,180,ERR=14) XM+T, YM+T
          WRITE (LUN,170,ERR=14) XM+T, YM-T
          WRITE (LUN,180,ERR=14) XM-T, YM+T
        ENDIF
    5   CONTINUE
C
C Paint the path.
C
      WRITE (LUN,150,ERR=14)
C
C Select a font and scale it.
C
      T = FSIZ
      WRITE (LUN,200,ERR=14) T
  200 FORMAT ('/Helvetica findfont'/
     .        F12.6,' scalefont setfont')
C
C Draw tick marks and labels at the extreme points of the
C   domain and range.  The tick mark length is T/2.
C
      XM = SFX*WX1+TX
      YM = REAL(IPY1)
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,180,ERR=14) XM, YM-T/2.0
      XM = XM - 1.5*T
      YM = YM - 1.5*T
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,210,ERR=14) WX1
  210 FORMAT ('(',F6.2,') show')
C
      XM = SFX*WX2+TX
      YM = REAL(IPY1)
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,180,ERR=14) XM, YM-T/2.0
      XM = XM - 1.5*T
      YM = YM - 1.5*T
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,210,ERR=14) WX2
C
      XM = REAL(IPX1)
      YM = SFY*WY1+TY
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,180,ERR=14) XM-T/2.0, YM
      XM = XM - 4.0*T
      YM = YM - 0.25*T
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,210,ERR=14) WY1
C
      XM = REAL(IPX1)
      YM = SFY*WY2+TY
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,180,ERR=14) XM-T/2.0, YM
      XM = XM - 4.0*T
      YM = YM - 0.25*T
      WRITE (LUN,170,ERR=14) XM, YM
      WRITE (LUN,210,ERR=14) WY2
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,220,ERR=14)
  220 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,230,ERR=14) CHAR(4)
  230 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter LUN, N, M, I1, I2, or PLTSIZ.
C
   11 IER = 1
      RETURN
C
C IND is not strictly increasing or contains an invalid
C   index.
C
   12 IER = 2
      RETURN
C
C DX = DY = 0:  all vertices identical.
C
   13 IER = 3
      RETURN
C
C Error writing to unit LUN.
C
   14 IER = 4
      RETURN
      END
