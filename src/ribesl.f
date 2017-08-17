cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE RIBESL(X,ALPHA,NB,IZE,B,NCALC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    this routine calculates bessel functions i sub(n+alpha) (x)
c    for non-negative argument x, and non-negative order n+alpha,
c    with or without exponential scaling.
c
c  on input:
c    x: working precision non-negative real argument for which
c       i's or exponentially scaled i's (i*exp(-x))
c       are to be calculated.  if i's are to be calculated,
c       x must be less than exparg (see below).
c    alpha: working precision fractional part of order for which
c           i's or exponentially scaled i's (i*exp(-x)) are
c           to be calculated.  0 .le. alpha .lt. 1.0.
c    nb: integer number of functions to be calculated, nb .gt. 0.
c        the first function calculated is of order alpha, and the
c        last is of order (nb - 1 + alpha).
c    ize: integer type.  ize = 1 if unscaled i's are to calculated,
c         and 2 if exponentially scaled i's are to be calculated.
c    b: working precision output vector of length nb.  if the routine
c       terminates normally (ncalc=nb), the vector b contains the
c       functions i(alpha,x) through i(nb-1+alpha,x), or the
c       corresponding exponentially scaled functions.
c    ncalc: integer output variable indicating possible errors.
c           before using the vector b, the user should check that
c           ncalc=nb, i.e., all orders have been calculated to
c           the desired accuracy.  see error returns below.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c explanation of machine-dependent constants
c
c   beta   = radix for the floating-point system
c   minexp = smallest representable power of beta
c   maxexp = smallest power of beta that overflows
c   it     = number of bits in the mantissa of a working precision
c            variable
c   nsig   = decimal significance desired.  should be set to
c            int(log10(2)*it+1).  setting nsig lower will result
c            in decreased accuracy while setting nsig higher will
c            increase cpu time without increasing accuracy.  the
c            truncation error is limited to a relative error of
c            t=.5*10**(-nsig).
c   enten  = 10.0 ** k, where k is the largest integer such that
c            enten is machine-representable in working precision
c   ensig  = 10.0 ** nsig
c   rtnsig = 10.0 ** (-k) for the smallest integer k such that
c            k .ge. nsig/4
c   enmten = smallest abs(x) such that x/4 does not underflow
c   xlarge = upper limit on the magnitude of x when ize=2.  bear
c            in mind that if abs(x)=n, then at least n iterations
c            of the backward recursion will be executed.  the value
c            of 10.0 ** 4 is used on every machine.
c   exparg = largest working precision argument that the library
c            exp routine can handle and upper limit on the
c            magnitude of x when ize=1; approximately
c            log(beta**maxexp)
c
c
c     approximate values for some important machines are:
c
c                        beta       minexp      maxexp       it
c
c  cray-1        (s.p.)    2        -8193        8191        48
c  cyber 180/855
c    under nos   (s.p.)    2         -975        1070        48
c  ieee (ibm/xt,
c    sun, etc.)  (s.p.)    2         -126         128        24
c  ieee (ibm/xt,
c    sun, etc.)  (d.p.)    2        -1022        1024        53
c  ibm 3033      (d.p.)   16          -65          63        14
c  vax           (s.p.)    2         -128         127        24
c  vax d-format  (d.p.)    2         -128         127        56
c  vax g-format  (d.p.)    2        -1024        1023        53
c
c
c                        nsig       enten       ensig      rtnsig
c
c cray-1        (s.p.)    15       1.0e+2465   1.0e+15     1.0e-4
c cyber 180/855
c   under nos   (s.p.)    15       1.0e+322    1.0e+15     1.0e-4
c ieee (ibm/xt,
c   sun, etc.)  (s.p.)     8       1.0e+38     1.0e+8      1.0e-2
c ieee (ibm/xt,
c   sun, etc.)  (d.p.)    16       1.0d+308    1.0d+16     1.0d-4
c ibm 3033      (d.p.)     5       1.0d+75     1.0d+5      1.0d-2
c vax           (s.p.)     8       1.0e+38     1.0e+8      1.0e-2
c vax d-format  (d.p.)    17       1.0d+38     1.0d+17     1.0d-5
c vax g-format  (d.p.)    16       1.0d+307    1.0d+16     1.0d-4
c
c
c                         enmten      xlarge   exparg
c
c cray-1        (s.p.)   1.84e-2466   1.0e+4    5677
c cyber 180/855
c   under nos   (s.p.)   1.25e-293    1.0e+4     741
c ieee (ibm/xt,
c   sun, etc.)  (s.p.)   4.70e-38     1.0e+4      88
c ieee (ibm/xt,
c   sun, etc.)  (d.p.)   8.90d-308    1.0d+4     709
c ibm 3033      (d.p.)   2.16d-78     1.0d+4     174
c vax           (s.p.)   1.17e-38     1.0e+4      88
c vax d-format  (d.p.)   1.17d-38     1.0d+4      88
c vax g-format  (d.p.)   2.22d-308    1.0d+4     709
c
c*******************************************************************
c*******************************************************************
c
c error returns
c
c  in case of an error,  ncalc .ne. nb, and not all i's are
c  calculated to the desired accuracy.
c
c  ncalc .lt. 0:  an argument is out of range. for example,
c     nb .le. 0, ize is not 1 or 2, or ize=1 and abs(x) .ge. exparg.
c     in this case, the b-vector is not calculated, and ncalc is
c     set to min0(nb,0)-1 so that ncalc .ne. nb.
c
c  nb .gt. ncalc .gt. 0: not all requested function values could
c     be calculated accurately.  this usually occurs because nb is
c     much larger than abs(x).  in this case, b(n) is calculated
c     to the desired accuracy for n .le. ncalc, but precision
c     is lost for ncalc .lt. n .le. nb.  if b(n) does not vanish
c     for n .gt. ncalc (because it is too small to be represented),
c     and b(n)/b(ncalc) = 10**(-k), then only the first nsig-k
c     significant figures of b(n) can be trusted.
c
c
c intrinsic functions required are:
c
c     dble, exp, dgamma, gamma, int, max, min, real, sqrt
c
c
c acknowledgement
c
c  this program is based on a program written by david j.
c  sookne (2) that computes values of the bessel functions j or
c  i of real argument and integer order.  modifications include
c  the restriction of the computation to the i bessel function
c  of non-negative real argument, the extension of the computation
c  to arbitrary positive order, the inclusion of optional
c  exponential scaling, and the elimination of most underflow.
c  an earlier version was published in (3).
c
c references: "a note on backward recurrence algorithms," olver,
c              f. w. j., and sookne, d. j., math. comp. 26, 1972,
c              pp 941-947.
c
c             "bessel functions of real argument and integer order,"
c              sookne, d. j., nbs jour. of res. b. 77b, 1973, pp
c              125-132.
c
c             "algorithm 597, sequence of modified bessel functions
c              of the first kind," cody, w. j., trans. math. soft.,
c              1983, pp. 242-245.
c
c  latest modification: may 30, 1989
c
c  modified by: w. j. cody and l. stoltz
c               applied mathematics division
c               argonne national laboratory
c               argonne, il  60439
c
c-------------------------------------------------------------------
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 IZE,K,L,MAGX,N,NB,NBMX,NCALC,NEND,NSIG,NSTART
      DOUBLE PRECISION DGAMMA,
     1 ALPHA,B,CONST,CONV,EM,EMPAL,EMP2AL,EN,ENMTEN,ENSIG,
     2 ENTEN,EXPARG,FUNC,HALF,HALFX,ONE,P,PLAST,POLD,PSAVE,PSAVEL,
     3 RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,X,XLARGE,ZERO
      DIMENSION B(NB)
c-------------------------------------------------------------------
c  mathematical constants
c-------------------------------------------------------------------
      DATA ONE,TWO,ZERO,HALF,CONST/1.0D0,2.0D0,0.0D0,0.5D0,1.585D0/
c-------------------------------------------------------------------
c  machine-dependent parameters
c-------------------------------------------------------------------
      DATA NSIG,XLARGE,EXPARG /16,1.0D4,709.0D0/
      DATA ENTEN,ENSIG,RTNSIG/1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN/8.9D-308/
c-------------------------------------------------------------------
c  statement functions for conversion
c-------------------------------------------------------------------
      CONV(N) = DBLE(N)
      FUNC(X) = DGAMMA(X)
c-------------------------------------------------------------------
c check for x, nb, or ize out of range.
c-------------------------------------------------------------------
      IF ((NB.GT.0) .AND. (X .GE. ZERO) .AND.
     1    (ALPHA .GE. ZERO) .AND. (ALPHA .LT. ONE) .AND.
     2    (((IZE .EQ. 1) .AND. (X .LE. EXPARG)) .OR.
     3     ((IZE .EQ. 2) .AND. (X .LE. XLARGE)))) THEN
c-------------------------------------------------------------------
c use 2-term ascending series for small x
c-------------------------------------------------------------------
            NCALC = NB
            MAGX = INT(X)
            IF (X .GE. RTNSIG) THEN
c-------------------------------------------------------------------
c initialize the forward sweep, the p-sequence of olver
c-------------------------------------------------------------------
                  NBMX = NB-MAGX
                  N = MAGX+1
                  EN = CONV(N+N) + (ALPHA+ALPHA)
                  PLAST = ONE
                  P = EN / X
c-------------------------------------------------------------------
c calculate general significance test
c-------------------------------------------------------------------
                  TEST = ENSIG + ENSIG
                  IF (2*MAGX .GT. 5*NSIG) THEN
                        TEST = DSQRT(TEST*P)
                     ELSE
                        TEST = TEST / CONST**MAGX
                  END IF
                  IF (NBMX .GE. 3) THEN
c-------------------------------------------------------------------
c calculate p-sequence until n = nb-1.  check for possible overflow.
c-------------------------------------------------------------------
                     TOVER = ENTEN / ENSIG
                     NSTART = MAGX+2
                     NEND = NB - 1
                     DO 100 K = NSTART, NEND
                        N = K
                        EN = EN + TWO
                        POLD = PLAST
                        PLAST = P
                        P = EN * PLAST/X + POLD
                        IF (P .GT. TOVER) THEN
c-------------------------------------------------------------------
c to avoid overflow, divide p-sequence by tover.  calculate
c p-sequence until abs(p) .gt. 1.
c-------------------------------------------------------------------
                           TOVER = ENTEN
                           P = P / TOVER
                           PLAST = PLAST / TOVER
                           PSAVE = P
                           PSAVEL = PLAST
                           NSTART = N + 1
   60                      N = N + 1
                              EN = EN + TWO
                              POLD = PLAST
                              PLAST = P
                              P = EN * PLAST/X + POLD
                           IF (P .LE. ONE) GO TO 60
                           TEMPB = EN / X
c-------------------------------------------------------------------
c calculate backward test, and find ncalc, the highest n
c such that the test is passed.
c-------------------------------------------------------------------
                           TEST = POLD*PLAST / ENSIG
                           TEST = TEST*(HALF-HALF/(TEMPB*TEMPB))
                           P = PLAST * TOVER
                           N = N - 1
                           EN = EN - TWO
                           NEND = MIN0(NB,N)
                           DO 80 L = NSTART, NEND
                              NCALC = L
                              POLD = PSAVEL
                              PSAVEL = PSAVE
                              PSAVE = EN * PSAVEL/X + POLD
                              IF (PSAVE*PSAVEL .GT. TEST) GO TO 90
   80                      CONTINUE
                           NCALC = NEND + 1
   90                      NCALC = NCALC - 1
                           GO TO 120
                        END IF
  100                CONTINUE
                     N = NEND
                     EN = CONV(N+N) + (ALPHA+ALPHA)
c-------------------------------------------------------------------
c calculate special significance test for nbmx .gt. 2.
c-------------------------------------------------------------------
                     TEST = MAX(TEST,DSQRT(PLAST*ENSIG)*DSQRT(P+P))
                  END IF
c-------------------------------------------------------------------
c calculate p-sequence until significance test passed.
c-------------------------------------------------------------------
  110             N = N + 1
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN * PLAST/X + POLD
                  IF (P .LT. TEST) GO TO 110
c-------------------------------------------------------------------
c initialize the backward recursion and the normalization sum.
c-------------------------------------------------------------------
  120             N = N + 1
                  EN = EN + TWO
                  TEMPB = ZERO
                  TEMPA = ONE / P
                  EM = CONV(N) - ONE
                  EMPAL = EM + ALPHA
                  EMP2AL = (EM - ONE) + (ALPHA + ALPHA)
                  SUM = TEMPA * EMPAL * EMP2AL / EM
                  NEND = N - NB
                  IF (NEND .LT. 0) THEN
c-------------------------------------------------------------------
c n .lt. nb, so store b(n) and set higher orders to zero.
c-------------------------------------------------------------------
                        B(N) = TEMPA
                        NEND = -NEND
                        DO 130 L = 1, NEND
  130                      B(N+L) = ZERO
                     ELSE
                        IF (NEND .GT. 0) THEN
c-------------------------------------------------------------------
c recur backward via difference equation, calculating (but
c not storing) b(n), until n = nb.
c-------------------------------------------------------------------
                           DO 140 L = 1, NEND
                              N = N - 1
                              EN = EN - TWO
                              TEMPC = TEMPB
                              TEMPB = TEMPA
                              TEMPA = (EN*TEMPB) / X + TEMPC
                              EM = EM - ONE
                              EMP2AL = EMP2AL - ONE
                              IF (N .EQ. 1) GO TO 150
                              IF (N .EQ. 2) EMP2AL = ONE
                              EMPAL = EMPAL - ONE
                              SUM = (SUM + TEMPA*EMPAL) * EMP2AL / EM
  140                      CONTINUE
                        END IF
c-------------------------------------------------------------------
c store b(nb)
c-------------------------------------------------------------------
  150                   B(N) = TEMPA
                        IF (NB .LE. 1) THEN
                           SUM = (SUM + SUM) + TEMPA
                           GO TO 230
                        END IF
c-------------------------------------------------------------------
c calculate and store b(nb-1)
c-------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N)  = (EN*TEMPA) / X + TEMPB
                        IF (N .EQ. 1) GO TO 220
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
                  END IF
                  NEND = N - 2
                  IF (NEND .GT. 0) THEN
c-------------------------------------------------------------------
c calculate via difference equation and store b(n), until n = 2.
c-------------------------------------------------------------------
                     DO 200 L = 1, NEND
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*B(N+1)) / X +B(N+2)
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
  200                CONTINUE
                  END IF
c-------------------------------------------------------------------
c calculate b(1)
c-------------------------------------------------------------------
                  B(1) = TWO*EMPAL*B(2) / X + B(3)
  220             SUM = (SUM + SUM) + B(1)
c-------------------------------------------------------------------
c normalize.  divide all b(n) by sum.
c-------------------------------------------------------------------
  230             IF (ALPHA .NE. ZERO)
     1               SUM = SUM * FUNC(ONE+ALPHA) * (X*HALF)**(-ALPHA)
                  IF (IZE .EQ. 1) SUM = SUM * DEXP(-X)
                  TEMPA = ENMTEN
                  IF (SUM .GT. ONE) TEMPA = TEMPA * SUM
                  DO 260 N = 1, NB
                     IF (B(N) .LT. TEMPA) B(N) = ZERO
                     B(N) = B(N) / SUM
  260             CONTINUE
                  RETURN
c-------------------------------------------------------------------
c two-term ascending series for small x.
c-------------------------------------------------------------------
               ELSE
                  TEMPA = ONE
                  EMPAL = ONE + ALPHA
                  HALFX = ZERO
                  IF (X .GT. ENMTEN) HALFX = HALF * X
                  IF (ALPHA .NE. ZERO) TEMPA = HALFX**ALPHA /FUNC(EMPAL)
                  IF (IZE .EQ. 2) TEMPA = TEMPA * DEXP(-X)
                  TEMPB = ZERO
                  IF ((X+ONE) .GT. ONE) TEMPB = HALFX * HALFX
                  B(1) = TEMPA + TEMPA*TEMPB / EMPAL
                  IF ((X .NE. ZERO) .AND. (B(1) .EQ. ZERO)) NCALC = 0
                  IF (NB .GT. 1) THEN
                     IF (X .EQ. ZERO) THEN
                           DO 310 N = 2, NB
                              B(N) = ZERO
  310                      CONTINUE
                        ELSE
c-------------------------------------------------------------------
c calculate higher-order functions.
c-------------------------------------------------------------------
                           TEMPC = HALFX
                           TOVER = (ENMTEN + ENMTEN) / X
                           IF (TEMPB .NE. ZERO) TOVER = ENMTEN / TEMPB
                           DO 340 N = 2, NB
                              TEMPA = TEMPA / EMPAL
                              EMPAL = EMPAL + ONE
                              TEMPA = TEMPA * TEMPC
                              IF (TEMPA .LE. TOVER*EMPAL) TEMPA = ZERO
                              B(N) = TEMPA + TEMPA*TEMPB / EMPAL
                              IF ((B(N) .EQ. ZERO) .AND. (NCALC .GT. N))
     1                             NCALC = N-1
  340                      CONTINUE
                     END IF
                  END IF
            END IF
         ELSE
            NCALC = MIN0(NB,0)-1
      END IF
      RETURN
c---------- last line of ribesl ----------
      END
c
