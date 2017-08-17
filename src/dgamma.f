cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION DGAMMA(X)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c this routine calculates the gamma function for a real argument x.
c   computation is based on an algorithm outlined in reference 1.
c   the program uses rational functions that approximate the gamma
c   function to at least 20 significant decimal digits.  coefficients
c   for the approximation over the interval (1,2) are unpublished.
c   those for the approximation for x .ge. 12 are from reference 2.
c   the accuracy achieved depends on the arithmetic system, the
c   compiler, the intrinsic functions, and proper selection of the
c   machine-dependent constants.
c
c
c*******************************************************************
c*******************************************************************
c
c explanation of machine-dependent constants
c
c beta   - radix for the floating-point representation
c maxexp - the smallest positive power of beta that overflows
c xbig   - the largest argument for which gamma(x) is representable
c          in the machine, i.e., the solution to the equation
c                  gamma(xbig) = beta**maxexp
c xinf   - the largest machine representable floating-point number;
c          approximately beta**maxexp
c eps    - the smallest positive floating-point number such that
c          1.0+eps .gt. 1.0
c xminin - the smallest positive floating-point number such that
c          1/xminin is machine representable
c
c     approximate values for some important machines are:
c
c                            beta       maxexp        xbig
c
c cray-1         (s.p.)        2         8191        966.961
c cyber 180/855
c   under nos    (s.p.)        2         1070        177.803
c ieee (ibm/xt,
c   sun, etc.)   (s.p.)        2          128        35.040
c ieee (ibm/xt,
c   sun, etc.)   (d.p.)        2         1024        171.624
c ibm 3033       (d.p.)       16           63        57.574
c vax d-format   (d.p.)        2          127        34.844
c vax g-format   (d.p.)        2         1023        171.489
c
c                            xinf         eps        xminin
c
c cray-1         (s.p.)   5.45e+2465   7.11e-15    1.84e-2466
c cyber 180/855
c   under nos    (s.p.)   1.26e+322    3.55e-15    3.14e-294
c ieee (ibm/xt,
c   sun, etc.)   (s.p.)   3.40e+38     1.19e-7     1.18e-38
c ieee (ibm/xt,
c   sun, etc.)   (d.p.)   1.79d+308    2.22d-16    2.23d-308
c ibm 3033       (d.p.)   7.23d+75     2.22d-16    1.39d-76
c vax d-format   (d.p.)   1.70d+38     1.39d-17    5.88d-39
c vax g-format   (d.p.)   8.98d+307    1.11d-16    1.12d-308
c
c*******************************************************************
c*******************************************************************
c
c error returns
c
c  the program returns the value xinf for singularities or
c     when overflow would occur.  the computation is believed
c     to be free of underflow and overflow.
c
c
c  intrinsic functions required are:
c
c     int, dble, exp, log, real, sin
c
c
c references: "an overview of software development for special
c              functions", w. j. cody, lecture notes in mathematics,
c              506, numerical analysis dundee, 1975, g. a. watson
c              (ed.), springer verlag, berlin, 1976.
c
c              computer approximations, hart, et. al., wiley and
c              sons, new york, 1968.
c
c  latest modification: october 12, 1989
c
c  authors: w. j. cody and l. stoltz
c           applied mathematics division
c           argonne national laboratory
c           argonne, il 60439
c
c----------------------------------------------------------------------
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 I,N
      LOGICAL PARITY
      DOUBLE PRECISION
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
c----------------------------------------------------------------------
c  mathematical constants
c----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
c----------------------------------------------------------------------
c  machine dependent parameters
c----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
c----------------------------------------------------------------------
c  numerator and denominator coefficients for rational minimax
c     approximation over (1,2).
c----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
c----------------------------------------------------------------------
c  coefficients for minimax approximation over (12, inf).
c----------------------------------------------------------------------
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
c----------------------------------------------------------------------
c  statement functions for conversion between integer and float
c----------------------------------------------------------------------
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
c----------------------------------------------------------------------
c  argument is negative
c----------------------------------------------------------------------
            Y = -X
            Y1 = DINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. DINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
c----------------------------------------------------------------------
c  argument is positive
c----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
c----------------------------------------------------------------------
c  argument .lt. eps
c----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
c----------------------------------------------------------------------
c  0.0 .lt. argument .lt. 1.0
c----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
c----------------------------------------------------------------------
c  1.0 .lt. argument .lt. 12.0, reduce argument if necessary
c----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
c----------------------------------------------------------------------
c  evaluate approximation for 1.0 .lt. argument .lt. 2.0
c----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
c----------------------------------------------------------------------
c  adjust result for case  0.0 .lt. argument .lt. 1.0
c----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
c----------------------------------------------------------------------
c  adjust result for case  2.0 .lt. argument .lt. 12.0
c----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
c----------------------------------------------------------------------
c  evaluate for argument .ge. 12.0,
c----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*DLOG(Y)
                  RES = DEXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
c----------------------------------------------------------------------
c  final adjustments and return
c----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
  900 DGAMMA = RES
      RETURN
c ---------- last line of gamma ----------
      END
c
