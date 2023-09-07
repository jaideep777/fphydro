module md_sofunutils_plus
  use md_precision
  implicit none

  contains 

  function zero ( a, b, f, t ) result(val)

    !*****************************************************************************80
    !
    !! ZERO seeks the root of a function F(X) in an interval [A,B].
    !
    !  Discussion:
    !
    !    The interval [A,B] must be a change of sign interval for F.
    !    That is, F(A) and F(B) must be of opposite signs.  Then
    !    assuming that F is continuous implies the existence of at least
    !    one value C between A and B for which F(C) = 0.
    !
    !    The location of the zero is determined to within an accuracy
    !    of 6 * MACHEPS * abs ( C ) + 2 * T.
    !
    !    Thanks to Thomas Secretin for pointing out a transcription error in the
    !    setting of the value of P, 11 February 2013.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    11 February 2013
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Richard Brent.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Richard Brent,
    !    Algorithms for Minimization Without Derivatives,
    !    Dover, 2002,
    !    ISBN: 0-486-41998-3,
    !    LC: QA402.5.B74.
    !
    !  Parameters:
    !
    !    Input, real (kind = dbl8) A, B, the endpoints of the change of 
    !    sign interval.
    !
    !    Input, real (kind = dbl8) MACHEP, an estimate for the relative machine
    !    precision.
    !
    !    Input, real (kind = dbl8) T, a positive error tolerance.
    !
    !    Input, external real (kind = dbl8) F, the name of a user-supplied
    !    function, of the form "FUNCTION F ( X )", which evaluates the
    !    function whose zero is being sought.
    !
    !    Output, real (kind = dbl8) ZERO, the estimated value of a zero of
    !    the function F.
    !
    implicit none

    real (kind = dbl8)  ::  a, b, c, d, e
    real (kind = dbl8)  ::  f
    real (kind = dbl8)  ::  fa, fb, fc
    real (kind = dbl8)  ::  m
    real (kind = dbl8)  ::  machep
    real (kind = dbl8)  ::  p, q, r, s, sa, sb
    real (kind = dbl8)  ::  t
    real (kind = dbl8)  ::  tol
    real (kind = dbl8)  ::  val

    machep = epsilon ( 1D+00 )
    !
    !  Make local copies of A and B.
    !
    sa = a
    sb = b
    fa = f ( sa )
    fb = f ( sb )

    c = sa
    fc = fa
    e = sb - sa
    d = e

    do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

      else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
      q = - q
      else
      p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
      else
      e = m
      d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
        ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

    end do

    val = sb

  end function zero

  function alngam ( xvalue, ifault )

    !*****************************************************************************80
    !
    !! ALNGAM computes the logarithm of the gamma function.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Allan Macleod.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Allan Macleod,
    !    Algorithm AS 245,
    !    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
    !    Applied Statistics,
    !    Volume 38, Number 2, 1989, pages 397-402.
    !
    !  Parameters:
    !
    !    Input, real (kind = dbl8) XVALUE, the argument of the Gamma function.
    !
    !    Output, integer (kind = int4) IFAULT, error flag.
    !    0, no error occurred.
    !    1, XVALUE is less than or equal to 0.
    !    2, XVALUE is too big.
    !
    !    Output, real (kind = dbl8) ALNGAM, the logarithm of the gamma function of X.
    !
    implicit none
    
    real (kind = dbl8) alngam
    real (kind = dbl8), parameter :: alr2pi = 0.918938533204673D+00
    integer (kind = int4) ifault
    real (kind = dbl8), dimension ( 9 ) :: r1 = (/ &
      -2.66685511495D+00, &
      -24.4387534237D+00, &
      -21.9698958928D+00, &
      11.1667541262D+00, &
      3.13060547623D+00, &
      0.607771387771D+00, &
      11.9400905721D+00, &
      31.4690115749D+00, &
      15.2346874070D+00 /)
    real (kind = dbl8), dimension ( 9 ) :: r2 = (/ &
      -78.3359299449D+00, &
      -142.046296688D+00, &
      137.519416416D+00, &
      78.6994924154D+00, &
      4.16438922228D+00, &
      47.0668766060D+00, &
      313.399215894D+00, &
      263.505074721D+00, &
      43.3400022514D+00 /)
    real (kind = dbl8), dimension ( 9 ) :: r3 = (/ &
      -2.12159572323D+05, &
      2.30661510616D+05, &
      2.74647644705D+04, &
      -4.02621119975D+04, &
      -2.29660729780D+03, &
      -1.16328495004D+05, &
      -1.46025937511D+05, &
      -2.42357409629D+04, &
      -5.70691009324D+02 /)
    real (kind = dbl8), dimension ( 5 ) :: r4 = (/ &
      0.279195317918525D+00, &
      0.4917317610505968D+00, &
      0.0692910599291889D+00, &
      3.350343815022304D+00, &
      6.012459259764103D+00 /)
    real (kind = dbl8) ::  x
    real (kind = dbl8) ::  x1
    real (kind = dbl8) ::  x2
    real (kind = dbl8), parameter :: xlge = 5.10D+05
    real (kind = dbl8), parameter :: xlgst = 1.0D+30
    real (kind = dbl8) xvalue
    real (kind = dbl8) y
    
    x = xvalue
    alngam = 0.0D+00
    !
    !  Check the input.
    !
    if ( xlgst <= x ) then
      ifault = 2
      return
    end if
    
    if ( x <= 0.0D+00 ) then
      ifault = 1
      return
    end if
    
    ifault = 0
    !
    !  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
    !
    if ( x < 1.5D+00 ) then
    
      if ( x < 0.5D+00 ) then
    
      alngam = - log ( x )
      y = x + 1.0D+00
    !
    !  Test whether X < machine epsilon.
    !
      if ( y == 1.0D+00 ) then
        return
      end if
    
      else
    
      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00
    
      end if
    
      alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
            y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )
    
      return
    
    end if
    !
    !  Calculation for 1.5 <= X < 4.0.
    !
    if ( x < 4.0D+00 ) then
    
      y = ( x - 1.0D+00 ) - 1.0D+00
    
      alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
            x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
    !
    !  Calculation for 4.0 <= X < 12.0.
    !
    else if ( x < 12.0D+00 ) then
    
      alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
            x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
    !
    !  Calculation for 12.0 <= X.
    !
    else
    
      y = log ( x )
      alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi
    
      if ( x <= xlge ) then
    
      x1 = 1.0D+00 / x
      x2 = x1 * x1
    
      alngam = alngam + x1 * ( ( &
          r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )
    
      end if
    
    end if
    
    return
    end

    
    function alnorm ( x, upper )
    
    !*****************************************************************************80
    !
    !! ALNORM computes the cumulative density of the standard normal distribution.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by David Hill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    David Hill,
    !    Algorithm AS 66:
    !    The Normal Integral,
    !    Applied Statistics,
    !    Volume 22, Number 3, 1973, pages 424-427.
    !
    !  Parameters:
    !
    !    Input, real (kind = dbl8) X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity;
    !    .FALSE. => integrate from - Infinity to X.
    !
    !    Output, real (kind = dbl8) ALNORM, the integral of the standard normal
    !    distribution over the desired interval.
    !
    implicit none
    
    real (kind = dbl8), parameter :: a1 = 5.75885480458D+00
    real (kind = dbl8), parameter :: a2 = 2.62433121679D+00
    real (kind = dbl8), parameter :: a3 = 5.92885724438D+00
    real (kind = dbl8) alnorm
    real (kind = dbl8), parameter :: b1 = -29.8213557807D+00
    real (kind = dbl8), parameter :: b2 = 48.6959930692D+00
    real (kind = dbl8), parameter :: c1 = -0.000000038052D+00
    real (kind = dbl8), parameter :: c2 = 0.000398064794D+00
    real (kind = dbl8), parameter :: c3 = -0.151679116635D+00
    real (kind = dbl8), parameter :: c4 = 4.8385912808D+00
    real (kind = dbl8), parameter :: c5 = 0.742380924027D+00
    real (kind = dbl8), parameter :: c6 = 3.99019417011D+00
    real (kind = dbl8), parameter :: con = 1.28D+00
    real (kind = dbl8), parameter :: d1 = 1.00000615302D+00
    real (kind = dbl8), parameter :: d2 = 1.98615381364D+00
    real (kind = dbl8), parameter :: d3 = 5.29330324926D+00
    real (kind = dbl8), parameter :: d4 = -15.1508972451D+00
    real (kind = dbl8), parameter :: d5 = 30.789933034D+00
    real (kind = dbl8), parameter :: ltone = 7.0D+00
    real (kind = dbl8), parameter :: p = 0.398942280444D+00
    real (kind = dbl8), parameter :: q = 0.39990348504D+00
    real (kind = dbl8), parameter :: r = 0.398942280385D+00
    logical up
    logical upper
    real (kind = dbl8), parameter :: utzero = 18.66D+00
    real (kind = dbl8) x
    real (kind = dbl8) y
    real (kind = dbl8) z
    
    up = upper
    z = x
    
    if ( z < 0.0D+00 ) then
      up = .not. up
      z = - z
    end if
    
    if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then
    
      if ( up ) then
      alnorm = 0.0D+00
      else
      alnorm = 1.0D+00
      end if
    
      return
    
    end if
    
    y = 0.5D+00 * z * z
    
    if ( z <= con ) then
    
      alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 & 
      / ( y + a3 ))))
    
    else
    
      alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))
    
    end if
    
    if ( .not. up ) then
      alnorm = 1.0D+00 - alnorm
    end if
    
    return
    end


    function gammad ( x, p, ifault )
    
    !*****************************************************************************80
    !
    !! GAMMAD computes the Lower Incomplete Gamma Integral y(a,x)/G(a)
    !
    !  Auxiliary functions:
    !
    !    ALOGAM = logarithm of the gamma function, 
    !    ALNORM = algorithm AS66
    !
    !  Modified:
    !
    !    20 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by B Shea.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    B Shea,
    !    Algorithm AS 239:
    !    Chi-squared and Incomplete Gamma Integral,
    !    Applied Statistics,
    !    Volume 37, Number 3, 1988, pages 466-473.
    !
    !  Parameters:
    !
    !    Input, real (kind = dbl8) X, P, the parameters of the incomplete 
    !    gamma ratio.  0 <= X, and 0 < P.
    !
    !    Output, integer (kind = int4) IFAULT, error flag.
    !    0, no error.
    !    1, X < 0 or P <= 0.
    !
    !    Output, real (kind = dbl8) GAMMAD, the value of the incomplete 
    !    Gamma integral.
    !
    implicit none
    
    real (kind = dbl8) a
    ! real (kind = dbl8) alnorm
    ! real (kind = dbl8) alngam
    real (kind = dbl8) an
    real (kind = dbl8) arg
    real (kind = dbl8) b
    real (kind = dbl8) c
    real (kind = dbl8), parameter :: elimit = - 88.0D+00
    real (kind = dbl8) gammad
    integer (kind = int4) ifault
    real (kind = dbl8), parameter :: oflo = 1.0D+37
    real (kind = dbl8) p
    real (kind = dbl8), parameter :: plimit = 1000.0D+00
    real (kind = dbl8) pn1
    real (kind = dbl8) pn2
    real (kind = dbl8) pn3
    real (kind = dbl8) pn4
    real (kind = dbl8) pn5
    real (kind = dbl8) pn6
    real (kind = dbl8) rn
    real (kind = dbl8), parameter :: tol = 1.0D-14
    logical upper
    real (kind = dbl8) x
    real (kind = dbl8), parameter :: xbig = 1.0D+08
    
    gammad = 0.0D+00
    !
    !  Check the input.
    !
    if ( x < 0.0D+00 ) then
      ifault = 1
      return
    end if
    
    if ( p <= 0.0D+00 ) then
      ifault = 1
      return
    end if
    
    ifault = 0
    
    if ( x == 0.0D+00 ) then
      gammad = 0.0D+00
      return
    end if
    !
    !  If P is large, use a normal approximation.
    !
    if ( plimit < p ) then
    
      pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) &
      + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
    
      upper = .false.
      gammad = alnorm ( pn1, upper )
      return
    
    end if
    !
    !  If X is large set GAMMAD = 1.
    !
    if ( xbig < x ) then
      gammad = 1.0D+00
      return
    end if
    !
    !  Use Pearson's series expansion.
    !  (Note that P is not large enough to force overflow in ALOGAM).
    !  No need to test IFAULT on exit since P > 0.
    !
    if ( x <= 1.0D+00 .or. x < p ) then
    
      arg = p * log ( x ) - x - alngam ( p + 1.0D+00, ifault )
      c = 1.0D+00
      gammad = 1.0D+00
      a = p
    
      do
    
      a = a + 1.0D+00
      c = c * x / a
      gammad = gammad + c
    
      if ( c <= tol ) then
        exit
      end if
    
      end do
    
      arg = arg + log ( gammad )
    
      if ( elimit <= arg ) then
      gammad = exp ( arg )
      else
      gammad = 0.0D+00
      end if
    !
    !  Use a continued fraction expansion.
    !
    else 
    
      arg = p * log ( x ) - x - alngam ( p, ifault )
      a = 1.0D+00 - p
      b = a + x + 1.0D+00
      c = 0.0D+00
      pn1 = 1.0D+00
      pn2 = x
      pn3 = x + 1.0D+00
      pn4 = x * b
      gammad = pn3 / pn4
    
      do
    
      a = a + 1.0D+00
      b = b + 2.0D+00
      c = c + 1.0D+00
      an = a * c
      pn5 = b * pn3 - an * pn1
      pn6 = b * pn4 - an * pn2
    
      if ( pn6 /= 0.0D+00 ) then
    
        rn = pn5 / pn6
    
        if ( abs ( gammad - rn ) <= min ( tol, tol * rn ) ) then
        exit
        end if
    
        gammad = rn
    
      end if
    
      pn1 = pn3
      pn2 = pn4
      pn3 = pn5
      pn4 = pn6
    !
    !  Re-scale terms in continued fraction if terms are large.
    !
      if ( oflo <= abs ( pn5 ) ) then
        pn1 = pn1 / oflo
        pn2 = pn2 / oflo
        pn3 = pn3 / oflo
        pn4 = pn4 / oflo
      end if
    
      end do
    
      arg = arg + log ( gammad )
    
      if ( elimit <= arg ) then
      gammad = 1.0D+00 - exp ( arg )
      else
      gammad = 1.0D+00
      end if
    
    end if
    
    return
    end
    
  end module md_sofunutils_plus
  