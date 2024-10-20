C
C     GEO2:  Discrete Periodic Geodesics on a Surface
C
C                         12/19/03
C
C
C   Given a smooth function f defining a level surface
C S = {p = (x,y,z): f(p) = 0} and an initial curve c in S,
C this program constructs a discrete approximation to the
C periodic geodesic of S that is closest to c.  The surface
C is assumed to be regular so that the unit normal fn(p) =
C grad(f)(p)/|grad(f)(p)| exists at every point p in S,
C where |_| denotes Euclidean norm in 3-space.
C
C   Denote by c:[0,L] -> S a parametric representation of a
C unit-speed curve in S (one parameterized by arc length)
C with total length L:  f(c(s)) = 0 and |c'(s)| = 1 for all
C s in [0,L].  The total geodesic curvature of c is the
C L2-norm squared of the tangential curvature component:
C
C       J(c) = Int[ |Kg(s)|**2 ds ],
C
C where Int[ ] denotes the integral with respect to arc
C length s from 0 to L, and
C
C       Kg(s) = [I-fn(c(s))*fn(c(s))**T]*c''(s)
C
C is the projection of the curvature vector c''(s) onto the
C orthogonal complement of the normal fn at c(s).  A geo-
C desic is characterized by J(c) = 0.  The periodic end
C conditions are c(0) = c(L) and c'(0) = c'(L).
C
C
C DISCRETIZATION:
C
C   Let c now denote the polygonal curve defined by a
C cyclically-ordered sequence of m vertices c(i), i = 1 to
C m.  We denote segment lengths by ds(i), midpoint unit
C tangent vectors by D1c(i), and vertex curvature vectors
C by D2c(i):
C
C       ds(i) = |c(i)-c(i-1)|,
C
C       D1c(i) = [c(i)-c(i-1)]/ds(i), and
C
C       D2c(i) = [D1c(i+1)-D1c(i)]/([ds(i)+ds(i+1)]/2)
C
C for i = 1 to m, where c(0) = c(m), ds(m+1) = ds(1), and
C D1c(m+1) = D1c(1).
C
C   The integral defining geodesic curvature is approximated
C by the composite trapezoidal rule giving the discretized
C energy functional
C
C       E(c) = 0.5*Sum[ (|u(i)|**2)*da(i) ],
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
C and Sum[_] denotes the sum over i = 1 to m.
C
C
C MINIMIZATION PROBLEM:
C
C   We have the constrained optimization problem of minimi-
C zing E(c) subject to the end conditions and the m nonlin-
C ear constraints F(c) = 0, where the i-th component of F(c)
C is f(c(i)), i = 1 to m.  The periodic end conditions are
C implicit in the expressions for ds(1), D1c(1), and D2c(m).
C
C   For the nonlinear constraints, we add a penalty function
C to the energy functional, giving
C
C       phi(c) = E(c) + 0.5*w*F(c)**T*F(c),
C
C where w is a positive penalty weight.  As w increases,
C the solution to the unconstrained problem of minimizing
C phi approaches the solution to the constrained problem.
C Unfortunately, the Hessian of phi also becomes more nearly
C singular, and the optimization problem becomes increasing-
C ly ill-conditioned, as w increases.  We therefore solve a
C sequence of problems with increasing values of w, using
C each solution as the initial estimate for the subsequent
C problem.
C
C
C METHOD:
C
C   The functional phi is minimized by the Polak-Ribiere
C nonlinear conjugate gradient method.  In place of the
C standard gradient (vector of partial derivatives), how-
C ever, we use a discretized Sobolev gradient.
C
C   Denote by S0 the set of perturbations for c which pre-
C serve the end conditions.  This is the space of all m-
C vectors of vertices.
C
C   A Sobolev inner product associated with a curve (the
C current approximation to the solution) c is
C
C       <g,h>_c = Int[ <g''(s),h''(s)> +
C                      a0**2*<g(s),h(s)> ds ]
C
C where s is the arc length associated with c, and a0 is a
C positive weight.  Note that this is positive (defines
C an inner product) on functions with square integrable
C second derivatives, and is intrinsic to the curve
C (independent of the parametrization).  The purpose of
C the second term is both to force positivity and to allow
C some control over the conditioning of the smoothing
C operator defined below.  The discretized inner product on
C S0 is
C
C       <g,h>_c = 0.5*Sum[ <D2g(i),D2h(i)>*da(i) +
C                          a0**2*<g(i),h(i)>*da(i) ]
C
C               = <Dg,Dh>_(2*m),
C
C where the sum is over i = 1 to m, D = (D2',D0') is the
C discrete differential operator defined by
C
C       D2'g(i) = Sqrt(0.5*da(i))*D2g(i),
C
C       D0'g(i) = a0*Sqrt(0.5*da(i))*g(i),
C
C and, for 3-vectors r and s, the discretized L2 inner
C product is
C
C       <r,s>_m = Sum(i=1,m)[ <r(i),s(i)> ].
C
C Note that D2g(i) involves second divided difference of g
C values and values of ds(i) = |c(i)-c(i-1)| -- not
C |g(i)-g(i-1)|, which would involve g in a nonlinear
C fashion.
C
C   We have <g,h>_c = <Dg,Dh>_(2*m) = <D**T*Dg,h>_m =
C <D**T*Dg,h>_S0.  Now let g be the Sobolev gradient (c-
C gradient) of phi at c.  Then by the Rietz Representation
C Theorem,
C
C       phi'(c)h = <grad(phi),h>_S0 = <g,h>_c
C                = <D**T*Dg,h>_S0
C
C for all h in S0.  Thus the Sobolev gradient g is defined
C by
C
C       DTD*g = grad(phi),
C
C for DTD = D**T*D.  Note that the 4-th order smoothing
C operator DTD**(-1) is a preconditioner for the gradient
C descent method.
C
C The operator D is represented by a 2*m by m matrix which
C is applied to each of the three components of an element
C of S0.  Note also that D**T*D is a symmetric positive
C definite matrix of order m.  The system is solved
C by a direct method (using an R**T*R factorization for
C upper triangular matrix R).
C
C
C ADDITIONAL SOFTWARE REQUIRED:
C
C   This program must be linked with the following level-1
C BLAS subprograms:  DASUM, DAXPY, DDOT, DSCAL.  These are
C often available in a pre-compiled library fine-tuned to a
C particular computer architecture.
C
C   Two additional subroutines are required:
C
C FEVAL (x,y,z,ND, F,Fx,Fy,Fz,Fxx,Fxy,Fxz,Fyy,Fyz,Fzz)
C defines the function f, returning a value F = f(x,y,z),
C first partial derivatives (Fx,Fy,Fz) if ND > 0, and second
C partial derivatives Fxx, ..., Fzz if ND > 1.  All parame-
C ters are double precision scalars.
C
C INITC (M, C) must return the sequence of vertices defining
C the initial curve.  For j = 1 to M, the Cartesian coordi-
C nates of c(j) must be stored in column j of array C which
C must be declared as follows:  'DOUBLE PRECISION C(3,0:M)'.
C (The coordinates of c(0) will be copied into the last
C column and therefore need not be stored.)  The vertices
C should be zeros of f (but are not tested).
C
C
C SOFTWARE USAGE INSTRUCTIONS:
C
C   In addition to selecting the surface f and initial curve
C c through Subroutines FEVAL and INITC, the number of ver-
C tices M (and hence the discretization error) may be chosen
C by altering the PARAMETER statement below.  Further
C options are provided in the DATA statement below.  Note,
C in particular, that the choice of descent method is con-
C trolled by L2G, MAXCG, OPT, SSIZ0, and SSIZ1.  Also, the
C penalty weight is controlled by W0, WSF, and WMAX, and
C creation of postscript files containing plots is con-
C trolled by PLT0, NPLOT, and PLTC.  These options can be
C selected or altered dynamically in the case of interactive
C execution mode (INTER = TRUE).  A run may be restarted by
C reading in the solution from a previous run as the initial
C estimate instead of calling INITC.  Refer to READC0.
C
C
C ARRAY length parameters:
C
C   LA = Length of array AP used to store D**T*D in Linpack
C        packed storage format for the linear solver DPPCO.
C
C   M = Number of vertices and curve segments.  The array
C       storage requirement is approximately 4*M**2 + 244*M
C       + 24*M1 bytes for 8-byte double precision numbers.
C       M >= 5.
C
C   M1 = Maximum number of vertices in an initial estimate
C        curve input by Subroutine READF (in array C1) if
C        READC0 = TRUE.
C
      INTEGER LA, M, M1
      PARAMETER (M=400, M1=400, LA=(M*(M+1))/2)
C
C Arrays:
C
C Note that C and C0 provide for vertex c(0) = c(m).  This
C   is required by Subroutine PLTCRV in order to connect
C   the first and last vertices.
C
      DOUBLE PRECISION AP(LA), C(3,0:M), C0(3,0:M),
     .                 C1(3,0:M1), CRV(M), D1C(3,M),
     .                 D2C(3,M), DC(3,M), DS(M), F(M),
     .                 G(3,M), G0(3,M), HU(3,M), U(3,M)
C
C Scalars and functions:
C
      CHARACTER CHR
      CHARACTER*9      FNAMO
      CHARACTER*8      FNAMP
      CHARACTER*80     TEXT
      DOUBLE PRECISION PHI
      DOUBLE PRECISION A0, ALPHA, CL, CMAX, CNRM2, DELC,
     .                 DG, DGG, DPHI, DSMAX, DSMIN, E,
     .                 FMTOL, FNS, GG, PLTSIZ, PHI0, PHIC,
     .                 PHTOL, Q, RCMIN, RCOND, SSIZ0, SSIZ1,
     .                 T, W, W0, WMAX, WSF
      INTEGER I, IER, J, J1, J2, K, LIN, LOUT, LPLT, LPRT,
     .        MAXCG, NCG, NG, NGMAX, NGMAX0, NPLOT, NPRNT,
     .        NPHSM, NPHVAL, NS
      LOGICAL INTER, L2G, OPT, PLOT, PLT0, PLTC, PRNT,
     .        READC0
C
C Parameters:
C
C   A0 = Positive weight defining the Sobolev inner product
C        <g,h>_c and the differential operator D.  A small
C        value places more weight on the second derivatives
C        resulting in a more effective smoother DTD**(-1),
C        but a larger value may be necessary to improve the
C        condition number, especially for large M.
C
C   FMTOL = Nonnegative tolerance for the line search:
C           length of the interval of uncertainty.  Refer
C           to Function FMIN.
C
C   INTER = Flag with value TRUE if and only if the program
C           is to be run interactively, in which case output
C           is written to the screen as well as to a file.
C
C   L2G = Flag with value TRUE iff the L2-gradient is to be
C         used in place of the Sobolev gradient.
C
C   MAXCG = Number of Polak-Ribiere conjugate gradient steps
C           between restarts with a steepest descent step.
C
C   NGMAX0 = Maximum number of descent steps allowed for
C            each w value if INTER = FALSE.  NGMAX0 > 0.
C
C   NPLOT = Number of descent steps between plot steps if
C           INTER = FALSE.  For each plot step, a postscript
C           file containing a plot of the current approxima-
C           tion c is created with the file name geoxx.ps,
C           where xx is the plot number.  A text file with
C           file name geoxx.out is also created.  The first
C           and last descent steps are always plot steps.
C
C   NPRNT = Number of descent steps between print steps.
C           For each print step, several lines of statistics
C           describing the step are written to unit LPRT.
C           The first and last descent steps are necessarily
C           print steps, and output is written to the screen
C           on every step if INTER = TRUE.  NPRNT .GE. 1.
C
C   OPT = Flag with value TRUE iff the optimal step-size
C         (computed by a line search) is to be used at each
C         descent step.
C
C   PHTOL = Nonnegative tolerance which defines convergence
C           of the descent method if INTER = FALSE:  bound
C           on Max(DPHI,DC**2,DG**3), where
C             DPHI = abs(change in phi)/(1+phi),
C             DC = max-norm(change in c)/(1+max-norm(c)),
C             DG = max-norm(S-gradient)/(1+max-norm(c)).
C
C   PLT0 = Flag with value TRUE iff the initial curve is to
C          be plotted and written to disk (files geo00.ps
C          and geo00.out are to be created).
C
C   PLTC = Flag with value TRUE iff the geodesic curvature
C          |Kg(s)| for the final curve c is to be plotted
C          and written to disk (files geocr.ps and geocr.out
C          are to be created).
C
C   PLTSIZ = Plot size in inches for Postscript plots.
C            1.0 <= PLTSIZ <= 6.5.
C
C   RCMIN = Minimum value of RCOND (reciprocal of the esti-
C           mated condition number of DTD) for which the
C           linear system is to be solved for the Sobolev
C           gradient.  A reasonable value is 1.D3*EPS, where
C           EPS is the machine precision, guaranteeing at
C           least 3 correct digits in the solution.  If
C           GRADS4 returns RCOND < RCMIN, the L2 gradient
C           is used.  RCMIN > 0.
C
C   READC0 = Flag with value TRUE iff the initial curve is
C            to be read from file geo2.in rather than
C            obtained by a call to Subroutine INITC.
C
C   SSIZ0 = Descent step size used if OPT = FALSE.
C
C   SSIZ1 = Initial descent step size used if OPT = TRUE.
C           Note that if this value is too large, the method
C           can fail with an increase in phi.
C
C   W0 = Initial value of the penalty weight w if INTER =
C        FALSE.  W0 >= 0.
C
C   WMAX = Final value of the weight w if INTER = FALSE.
C          WMAX >= W0.
C
C   WSF = Scale factor for increasing w if INTER = FALSE.
C         WSF > 1.
C
      DATA A0/1.0D1/,      FMTOL/1.D-6/,   INTER/.FALSE./,
     .     L2G/.FALSE./,   MAXCG/3/,       NGMAX0/20000/,
     .     NPLOT/20000/,   NPRNT/200/,     OPT/.TRUE./,
     .     PHTOL/1.D-11/,  PLT0/.TRUE./,   PLTC/.TRUE./,
     .     PLTSIZ/3.5/,    RCMIN/1.D-13/,  READC0/.FALSE./,
     .     SSIZ0/1.D-3/,   SSIZ1/1.D-5/,   W0/1.D0/,
     .     WMAX/1.D3/,     WSF/2.0D0/
C
C Logical unit numbers:
C
C   LIN = Logical unit number for reading an initial curve
C          M, ((C(I,J), I = 1,3), J = 0,M) if READC0 = TRUE:
C          Subroutine READF.
C
C   LOUT = Logical unit number for writing curves
C          M, ((C(I,J), I = 1,3), J = 0,M):  Subroutine
C          WRITF.
C
C   LPLT = Logical unit number for files containing
C          Postscript plots:  initial estimate, geodesic
C          curve, curvature, etc.  Refer to Subroutine
C          PLTCRV.
C
C   LPRT = Logical unit number for output other than the
C          solution.  If INTER = TRUE, output is also
C          written to the screen.
C
      DATA LIN/1/,  LOUT/2/,  LPLT/3/,  LPRT/4/
C
C File names for postscript files containing plots of
C   curves and for output text files.  Each time a new
C   plot is interactively requested (or after each NPLOT
C   descent steps if INTER = FALSE), K is incremented and
C   inserted into the file names (in place of 00).
C
      DATA K/0/,  FNAMO/'geo00.out'/,  FNAMP/'geo00.ps'/
C
C Open print file, delete it, and reopen it.  This is
C   necessary on some systems (Microsoft) which ignore
C   end-of-file records.
C
      OPEN (LPRT,FILE='geo2.prt')
      CLOSE (LPRT,STATUS='DELETE')
      OPEN (LPRT,FILE='geo2.prt')
C
C Interactive input formats:
C
  500 FORMAT (I6)
  510 FORMAT (D22.15)
  520 FORMAT (A80)
  530 FORMAT (A1)
      IF (INTER) THEN
C
C Interactive mode:  initialize NGMAX to 1.
C
        NGMAX = 1
C
C Print a message and prompt for a line of text to be
C   written to the print file.
C
        WRITE (*,100)
  100   FORMAT (///10X,'GEO2:  Periodic Geodesic Curve',
     .                 ' Package'/
     .          /10X,'Respond to each of the following ',
     .              'prompts with one of:'//
     .          15X,'a)  the requested value followed by ',
     .              'Carriage Return,'/
     .          15X,'b)  Carriage Return only to specify ',
     .              'a default value, or'/
     .          15X,'c)  an invalid character (such as a ',
     .              'letter when numeric'/
     .          19X,'input is requested) to back up to the',
     .              ' previous option.'///
     .          10X,'Specify a line of text with at most ',
     .              '80 characters (the'/
     .          10X,'current date for example) to be ',
     .              'included in geo2.prt:'/)
        READ (*,520) TEXT
      ELSE
C
C Non-interactive mode:  initialize NGMAX = NGMAX0.
C
        NGMAX = NGMAX0
      ENDIF
C
C Print a heading with parameter values.
C
      WRITE (LPRT,110)
  110 FORMAT (///26X,'GEO2 Output'/)
      IF (INTER) WRITE (LPRT,112) TEXT
  112 FORMAT (5X,A80)
      WRITE (LPRT,114) M, W0, WSF, WMAX, L2G, A0, RCMIN,
     .                 MAXCG, OPT, SSIZ0, SSIZ1, FMTOL
  114 FORMAT (//5X,'Number of curve segments:  M = ',I4//1P,
     .        5X,'Penalty weight W parameters:'/
     .        10X,'Initial value:  W0 = ',D9.3/
     .        10X,'Scale factor:  WSF = ',D9.3/
     .        10X,'Final value:  WMAX = ',D9.3//
     .        5X,'Descent method parameters:'/
     .        10X,'L2 gradient flag:  L2G = ',L1/
     .        10X,'Weight on vertex values in D:  A0 = ',D9.3/
     .        10X,'Minimum condition # of DTD:  RCMIN = ',D9.3/
     .        10X,'# CG steps between restarts:  MAXCG = ',I3/
     .        10X,'Optimal step-size flag:  OPT = ',L1/
     .        10X,'Constant step-size if OPT = F:  SSIZ0 = ',D9.3/
     .        10X,'Initial step-size if OPT = T:  SSIZ1 = ',D9.3/
     .        10X,'Line search tolerance:  FMTOL = ',D9.3)
      IF (INTER) THEN
        WRITE (*,110)
        WRITE (*,112) TEXT
        WRITE (*,114) M, W0, WSF, WMAX, L2G, A0, RCMIN,
     .                MAXCG, OPT, SSIZ0, SSIZ1, FMTOL
      ELSE
        WRITE (LPRT,116) PHTOL, NGMAX0
  116   FORMAT (/5X,'Convergence parameters:  '/
     .          10X,'Tolerance:  PHTOL = ',1P,D9.3/
     .          10X,'Max # of descent steps per w value:  ',
     .              'NGMAX0 = ',I7)
      ENDIF
C
C Get initial vertices C.
C
      IF (READC0) THEN
        OPEN (LIN,FILE='geo2.in',ERR=40,STATUS='OLD')
        CALL READF (LIN,M,M1, C,C1,IER)
        IF (IER .NE. 0) GO TO 41
      ELSE
        CALL INITC (M, C)
      ENDIF
C
C Optional plot of initial estimate:  PLTCRV.
C
C J1 and J2 define orthogonal projections as follows:
C
C   J1 = 1, J2 = 2  ==>  x-y plane,
C   J1 = 2, J2 = 3  ==>  y-z plane,
C   J1 = 3, J2 = 1  ==>  z-x plane.
C
C Set the second parameter N to -1 to suppress marker
C   symbols.
C
      IF (PLT0) THEN
        J1 = 1
        J2 = 2
        OPEN (LPLT,FILE=FNAMP)
        CLOSE (LPLT,STATUS='DELETE')
        OPEN (LPLT,FILE=FNAMP)
        CALL PLTCRV (LPLT,0,0,M,C,J1,J2,PLTSIZ, IER)
        CLOSE (LPLT)
        IF (IER .NE. 0) GO TO 46
C
C Write the initial curve to disk.
C
        OPEN (LOUT,FILE=FNAMO)
        CLOSE (LOUT,STATUS='DELETE')
        OPEN (LOUT,FILE=FNAMO)
        CALL WRITF (LOUT,M,C)
        CLOSE (LOUT)
      ENDIF
C
C Initialize iteration counts and weight w.
C
C   NG = Number of iterations (gradient evaluations).
C
C   NGMAX = Maximum allowable value of NG with the current
C           w value.
C
C   NPHSM = Total number of evaluations of phi in the line
C           searches.
C
C   W = Weight defining penalty term.
C
      NG = 0
      NPHSM = 0
      W = W0
C
C DEBUG ONLY - Prints M, C, and W.
C
      WRITE (LPRT,118) M
  118 FORMAT (//5X,':3 Number of curve segments:  M = ',I4)
      WRITE (LPRT,120) W
  199 FORMAT (5X,':3 Penalty weight:  w = ',1P,D9.3)
      WRITE (LPRT,122) ((C(I,J), I = 1,3), J = 0,M)
C Prints with each vertex on a separate line.
  122 FORMAT (5X,':33 Vertices of initial curve:  C(j) = ',
     .        3(1P,D35.15))

C
C Top of outer loop on weights w.
C
C Compute ds(i), D1c(i), D2c(i), u(i), F(i) = f(c(i)),
C   G(i) = grad(f)(c(i)), curve length CL, E = E(c),
C   FNS = .5*||F||**2, and PHI0 = phi(c).
C
    2 PHI0 = PHI (M,C,W, DS,D1C,D2C,U,F,G,HU,CL,E,FNS,IER)
      IF (IER .NE. 0) GO TO 43
C
C Initialization for inner loop:
C
C   NCG = Number of conjugate gradient steps since the
C         previous restart with a steepest descent step.
C
      NCG = MAXCG
      WRITE (LPRT,120) W
      IF (INTER) WRITE (*,120) W
  120 FORMAT (//5X,'Penalty weight:  w = ',1P,D9.3)
C
C Top of inner loop:  overwrite G with the L2-gradient
C                     grad(phi).
C
    4 CALL GRADL2 (M,W,DS,D1C,D2C,U,F,HU, G, IER)
      IF (IER .NE. 0) GO TO 44
C
C Overwrite G with the Sobolev gradient unless L2G = TRUE.
C   Arrays CRV and U are used as work space.
C
      RCOND = 0.D0
      IF (.NOT. L2G) THEN
        CALL GRADS4 (M,DS,A0,RCMIN, AP,G,CRV,U, RCOND,IER)
        IF (IER .GT. 0) THEN
          WRITE (LPRT,125) RCOND
          IF (INTER) WRITE (*,125) RCOND
  125     FORMAT (//5X,'GRADS4 Failure:  RCOND = ',1P,D9.3)
        ENDIF
      ENDIF
C
C Update iteration count.
C
      NG = NG + 1
C
C Print statistics associated with the step.
C
C   PRNT = TRUE iff statistics describing the iteration are
C          to be printed.
C
      PRNT = (NPRNT*(NG/NPRNT) .EQ. NG)  .OR.  NG .EQ. 1
      IF (PRNT) THEN
        IF (NCG .GE. MAXCG) THEN
          WRITE (LPRT,130) NG
        ELSE
          WRITE (LPRT,135) NG
        ENDIF
        WRITE (LPRT,140) PHI0, CL, E, FNS, RCOND
      ENDIF
      IF (INTER) THEN
        IF (NCG .GE. MAXCG) THEN
          WRITE (*,130) NG
        ELSE
          WRITE (*,135) NG
        ENDIF
        WRITE (*,140) PHI0, CL, E, FNS, RCOND
      ENDIF
  130 FORMAT (//15X,'Iteration ',I7,' (Steepest descent)')
  135 FORMAT (//15X,'Iteration ',I7,' (Conjugate gradient)')
  140 FORMAT (5X,'Functional:                  phi(c) = ',
     .           1P,D21.15/
     .        5X,'              Total curve length CL = ',
     .           D21.15/
     .        5X,'                        Energy E(c) = ',
     .           D21.15/
     .        5X,'    Constraint residual .5*||F||**2 = ',
     .           D9.3/
     .        5X,'Reciprocal of Condition No.:  RCOND = ',
     .           D9.3)
C
      IF (NCG .GE. MAXCG) THEN
C
C Steepest descent step.
C
        Q = 0.0
        NCG = 0
      ELSE
C
C Conjugate gradient step.
C
        GG = 0.D0
        DGG = 0.D0
        DO 5 J = 1,M
          GG = GG + G0(1,J)*G0(1,J) + G0(2,J)*G0(2,J) +
     .              G0(3,J)*G0(3,J)
          DGG = DGG + G(1,J)*(G(1,J)-G0(1,J)) +
     .                G(2,J)*(G(2,J)-G0(2,J)) +
     .                G(3,J)*(G(3,J)-G0(3,J))
    5     CONTINUE
        Q = DGG/GG
        NCG = NCG + 1
      ENDIF
      IF (NCG .LT. MAXCG) THEN
C
C Save a copy of the gradient G in G0 for computing the
C   search direction in the next conjugate gradient step.
C
        DO 6 J = 1,M
          G0(1,J) = G(1,J)
          G0(2,J) = G(2,J)
          G0(3,J) = G(3,J)
    6     CONTINUE
      ENDIF
C
C Compute the search direction DC = -G + Q*DC, and
C   compute the max-norm DG of the gradient G.
C
      DG = 0.D0
      DO 7 J = 1,M
        DC(1,J) = -G(1,J) + Q*DC(1,J)
        DC(2,J) = -G(2,J) + Q*DC(2,J)
        DC(3,J) = -G(3,J) + Q*DC(3,J)
        DG = MAX(DG,ABS(G(1,J)),ABS(G(2,J)),ABS(G(3,J)))
    7   CONTINUE
C
C Copy C into C0, compute the optimal step-size ALPHA, and
C   update C to C0+ALPHA*DC.  The Max-norm of C is returned
C   in CMAX, and DELC is the relative change in C:
C   Max-norm(C-C0)/(1+CMAX) = ALPHA*Max-norm(DC)/(1+CMAX).
C   G is overwritten with grad(f)(c(i)).
C
C The initial estimate for ALPHA should be small.  Using
C   the step-size computed at the previously iteration can
C   result in an incorrect step-size caused by finding a
C   local minimum which is not a global minimum.
C
      C0(1,0) = C(1,0)
      C0(2,0) = C(2,0)
      C0(3,0) = C(3,0)
      DO 8 J = 1,M
        C0(1,J) = C(1,J)
        C0(2,J) = C(2,J)
        C0(3,J) = C(3,J)
    8   CONTINUE
      IF (OPT) THEN
        ALPHA = SSIZ1
      ELSE
C
C   Use small constant step-size.
C
        ALPHA = SSIZ0
      ENDIF
      CALL LNSRCH (M,C0,W,PHI0,DC,FMTOL,OPT, ALPHA, C,DS,
     .             D1C,D2C,U,F,G,HU,CL,E,FNS,PHIC,NPHVAL,
     .             CMAX,DELC,IER)
      IF (IER .NE. 0) GO TO 45
      NPHSM = NPHSM + NPHVAL
C
C Adjust DELC to the relative change in C squared.
C
      DELC = DELC**2
C
C Adjust DG to the cube of the max-norm of the gradient
C   relative to the solution c.
C
      DG = (DG/(1.D0+CMAX))**3
C
C Compute the relative change DPHI in phi and update PHI0.
C   Terminate on an increase in phi.
C
      DPHI = (PHI0-PHIC)/(1.D0+PHI0)
      PHI0 = PHIC
      IF (DPHI .LT. 0) THEN
        NGMAX = NG
        IF (.NOT. INTER) W = WMAX
      ENDIF
C
C Print statistics.
C
      IF (PRNT) WRITE (LPRT,150) NPHVAL, ALPHA, PHIC, DPHI,
     .                           DELC, DG
      IF (INTER) WRITE (*,150) NPHVAL, ALPHA, PHIC, DPHI,
     .                         DELC, DG
  150 FORMAT (/5X,'Line search:   No. phi evals NPHVAL = ',
     .            I3/
     .         5X,'                        Step-size S = ',
     .            1P,D9.3/
     .         5X,'                  Functional phi(c) = ',
     .            D21.15/
     .         5X,'     Relative change in phi(c) DPHI = ',
     .            D9.2/
     .         5X,'  Squared relative change in c:  DC = ',
     .            D9.2/
     .         5X,'    Cubed S-gradient rel. to c:  DG = ',
     .            D9.2/
     .         1X,'______________________________________',
     .            '___________________________')
C
C Test for a plot step.
C
      IF (.NOT. INTER) PLOT = NPLOT*(NG/NPLOT) .EQ. NG
      IF (PLOT  .AND.  K .EQ. 99) THEN
        WRITE (LPRT,155)
        IF (INTER) WRITE (*,155)
  155   FORMAT (/5X,'*** No more .ps or .out files can be',
     .              ' created. ***')
      ENDIF
      IF (PLOT  .AND.  K .LT. 99) THEN
C
C Create a Postscript file for the current approximation.
C
        J1 = 1
        J2 = 2
        K = K + 1
        I = K/10
        FNAMP(4:4) = CHAR(I+48)
        FNAMO(4:4) = CHAR(I+48)
        I = K - 10*I
        FNAMP(5:5) = CHAR(I+48)
        FNAMO(5:5) = CHAR(I+48)
        OPEN (LPLT,FILE=FNAMP)
        CLOSE (LPLT,STATUS='DELETE')
        OPEN (LPLT,FILE=FNAMP)
        CALL PLTCRV (LPLT,0,0,M,C,J1,J2,PLTSIZ, IER)
        CLOSE (LPLT)
        IF (IER .NE. 0) GO TO 46
C
C Create a text file:
C
        OPEN (LOUT,FILE=FNAMO)
        CLOSE (LOUT,STATUS='DELETE')
        OPEN (LOUT,FILE=FNAMO)
        CALL WRITF (LOUT,M,C)
        CLOSE (LOUT)
        PLOT = .FALSE.
      ENDIF
C
C Test for termination of inner loop.
C
      IF (.NOT. INTER) THEN
C
C   Noninteractive mode.
C
        IF ((DPHI .LT. PHTOL  .AND.  DELC .LT. PHTOL  .AND.
     .      DG .LT. PHTOL)  .OR.  NG .GE. NGMAX) GO TO 15
      ELSE
C
C   Interactive mode.
C
        IF (NG .GE. NGMAX) THEN
    9     WRITE (*,160)
  160     FORMAT (/10X,'Change options? (y/n)'/)
          READ (*,530,ERR=9) CHR
          IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
     .        CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 9
          IF (CHR .EQ. 'N'  .OR.  CHR .EQ. 'n') GO TO 14
C
C   Change options.
C
   10     WRITE (*,161)
  161     FORMAT (10X,'Use the L2-gradient (in place of ',
     .                'the Sobolev gradient)? (y/n)'/)
          READ (*,530,ERR=9) CHR
          IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
     .        CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 10
          L2G = CHR .EQ. 'Y'  .OR.  CHR .EQ. 'y'
C
   11     WRITE (*,162)
  162     FORMAT (10X,'Use a line search in each descent',
     .                ' step? (y/n)'/)
          READ (*,530,ERR=10) CHR
          IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
     .        CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 11
          OPT = CHR .EQ. 'Y'  .OR.  CHR .EQ. 'y'
C
          IF (.NOT. OPT) THEN
   12       WRITE (*,163) SSIZ0
  163       FORMAT (10X,'Descent step-size SSIZ0 = ',1P,D9.3
     .             /10X,'Specify a new (positive) value:'/)
            READ (*,510,ERR=11) SSIZ0
            IF (SSIZ0 .LE. 0.D0) GO TO 12
          ENDIF
C
   13     WRITE (*,164) A0
  164     FORMAT (10X,'D0 weight a0 = ',1P,D9.3/
     .            10X,'Specify a new (positive) value:'/)
          T = A0
          READ (*,510,ERR=11) A0
          IF (A0 .LE. 0.D0) GO TO 13
          IF (A0 .NE. T) WRITE (LPRT,165) A0
          WRITE (*,165) A0
  165     FORMAT (//5X,'D0 weight:  a0 = ',1P,D9.3/)
C
   14     WRITE (*,166)
  166     FORMAT (10X,'Create a .ps file containing a ',
     .                'plot of the next curve? (y/n)'/)
          READ (*,530,ERR=13) CHR
          IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
     .        CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 14
          PLOT = CHR .EQ. 'Y'  .OR.  CHR .EQ. 'y'
C
          WRITE (*,167)
  167     FORMAT (/10X,'Take another NS descent steps'/10X,
     .             '(terminate step if NS < 0) for NS = '/)
          READ (*,500,ERR=14) NS
          IF (NS .LT. 0) GO TO 15
          IF (NS .EQ. 0) NS = 1
          NGMAX = NGMAX + NS
        ENDIF
      ENDIF
      GO TO 4
C
C End of inner loop.
C
   15 IF (.NOT. PRNT) THEN
C
C Print parameter values.
C
        IF (NCG .EQ. 0) THEN
          WRITE (LPRT,130) NG
        ELSE
          WRITE (LPRT,135) NG
        ENDIF
        WRITE (LPRT,140) PHI0, CL, E, FNS, RCOND
        WRITE (LPRT,150) NPHVAL, ALPHA, PHIC, DPHI, DELC, DG
      ENDIF
      IF (INTER) THEN
C
C Interactive mode:  prompt for termination or new weight w.
C
        W0 = W
   16   WRITE (*,170) W
  170   FORMAT (//5X,'Current weight:  W = ',F8.3//
     .          10X,'Terminate program? (y/n)'/)
        READ (*,530,ERR=16) CHR
        IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
     .      CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 16
        IF (CHR .EQ. 'Y'  .OR.  CHR .EQ. 'y') GO TO 30
   17   WRITE (*,172)
  172   FORMAT (10X,'Specify a new weight:  W ='/)
        READ (*,510,ERR=17) W
   18   WRITE (*,174)
  174   FORMAT (10X,'Take NS descent steps for NS = '/)
        READ (*,500,ERR=18) NS
        IF (NS .LE. 0) NS = 1
        NGMAX = NG + NS
        GO TO 2
      ELSE
C
C Noninteractive mode:  test for termination of outer loop
C   on weights w.
C
        IF (W .LT. WMAX) THEN
          W = WSF*W
          IF (W .GT. WMAX) W = WMAX
          NGMAX = NG + NGMAX0
          GO TO 2
        ENDIF
      ENDIF
C
C End of outer loop on weights.  Print iteration counts and
C   statistics.
C
   30 DSMIN = DS(1)
      DSMAX = DS(1)
      DO 31 J = 2,M
        DSMIN = MIN(DSMIN,DS(J))
        DSMAX = MAX(DSMAX,DS(J))
   31   CONTINUE
      WRITE (LPRT,180) DSMIN, DSMAX
      IF (INTER) WRITE (*,180) DSMIN, DSMAX
  180 FORMAT (//5X,1P,'Minimum ds(j) = ',D9.3,
     .             ', Maximum ds(j) = ',D9.3)
      WRITE (LPRT,182) NG
      IF (INTER) WRITE (*,182) NG
  182 FORMAT (/5X,'Number of descent ',
     .        'iterations: ',I6)
      WRITE (LPRT,184) NPHSM
      IF (INTER) WRITE (*,184) NPHSM
  184 FORMAT (/5X,'Total number of phi evaluations ',
     .        'in the line searches: ',I7)
C
C Plot the solution curve.
C
      J1 = 1
      J2 = 2
      OPEN (LPLT,FILE='geo100.ps')
      CLOSE (LPLT,STATUS='DELETE')
      OPEN (LPLT,FILE='geo100.ps')
      CALL PLTCRV (LPLT,0,0,M,C,J1,J2,PLTSIZ, IER)
      CLOSE (LPLT)
      IF (IER .NE. 0) GO TO 46
C
C Write the solution to disk.
C
      OPEN (LOUT,FILE='geo100.out')
      CLOSE (LOUT,STATUS='DELETE')
      OPEN (LOUT,FILE='geo100.out')
      CALL WRITF (LOUT,M,C)
      CLOSE (LOUT)
C
C Compute intrinsic properties of the curve:
C
C   F = Vertex values of cumulative arc length s.
C   CRV = Vertex values of geodesic curvature |Kg(s)| =
C         |u(i)|, i = 1 to m.
C
C   CNRM2 = J(c) = Integral with respect to s of Kg(s)**2.
C   CMAX = Maximum of |Kg(s)|.
C
      CALL CPROP (M,DS,U, F,CRV,CNRM2,CMAX,IER)
      IF (PLTC) THEN
C
C Copy F and CRV into the rows of C0.
C
        DO 32 J = 1,M
          C0(1,J) = F(J)
          C0(2,J) = CRV(J)
   32     CONTINUE
        C0(1,0) = 0.D0
        C0(2,0) = C0(2,M)
C
C Create a postscript file containing a plot of geodesic
C   curvature |Kg(s)| as a function of s.
C
        J1 = 1
        J2 = 2
        OPEN (LPLT,FILE='geocr.ps')
        CLOSE (LPLT,STATUS='DELETE')
        OPEN (LPLT,FILE='geocr.ps')
        CALL PLTCRV (LPLT,0,0,M,C0,J1,J2,PLTSIZ, IER)
        CLOSE (LPLT)
        IF (IER .NE. 0) GO TO 46
      ENDIF
C
C Print norms and terminate the program.
C
      WRITE (LPRT,190) F(M), CNRM2, CMAX
      IF (INTER) WRITE (*,190) F(M), CNRM2, CMAX
  190 FORMAT (/5X,1P,'Total curve length = ',D10.4/
     .         5X,'Integral of geodesic curvature Kg(s)**2',
     .            ' = ',D10.4/
     .         5X,'Maximum of |Kg(s)| = ',D10.4)
      STOP
C
C Error opening file geo2.in.
C
   40 WRITE (LPRT,420)
      IF (INTER) WRITE (*,420)
  420 FORMAT (///10X,'*** Error opening file geo2.in ***')
      STOP
C
C Error in Subroutine READF.
C
   41 WRITE (LPRT,421) IER
      IF (INTER) WRITE (*,421) IER
  421 FORMAT (///10X,'*** Error in READF:  IER = ',I1,
     .        ' ***')
      STOP
C
C Error encountered in Function PHI.
C
   43 WRITE (LPRT,423) IER
      IF (INTER) WRITE (*,423) IER
  423 FORMAT (///10X,'*** Error in PHI:  IER = ',I1,' ***')
C
C Error encountered in Subroutine GRADL2.
C
   44 WRITE (LPRT,424) IER
      IF (INTER) WRITE (*,424) IER
  424 FORMAT (///10X,'*** Error in GRADL2:  ',
     .        'IER = ',I1,' ***')
      STOP
C
C Error encountered in Subroutine LNSRCH.
C
   45 WRITE (LPRT,425) IER
      IF (INTER) WRITE (*,425) IER
  425 FORMAT (///10X,'*** Error in LNSRCH:  IER = ',
     .        I1,' ***')
      STOP
C
C Error encountered in Subroutine PLTCRV.
C
   46 WRITE (LPRT,426) IER
      IF (INTER) WRITE (*,426) IER
  426 FORMAT (///10X,'*** Error in PLTCRV:  IER = ',
     .        I1,' ***')
      STOP
      END