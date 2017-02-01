module pdf_fcns

! Module cotains subroutines which calculate functions used in pdfs
!
! gamma_logfcn
! gammafcn
! betafcn


contains
subroutine gamma_logfcn ( x, gamma_log )
!
!*******************************************************************************
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.  
!    The program uses rational functions that theoretically approximate 
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The 
!    approximation for X > 12 is from reference 3, while approximations 
!    for X < 12.0E+00 are similar to those in reference 1, but are unpublished.  
!    The accuracy achieved depend subroutines on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-depend subroutineent 
!    constants.
!
!  Modified:
!
!    16 June 1999
!
!  Authors: 
!
!    W. J. Cody and L. Stoltz
!    Argonne National Laboratory
!
!  References:
!
!    # 1) 
!    W. J. Cody and K. E. Hillstrom, 
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp. 
!    Volume 21, 1967, pages 198-203.
!
!    # 2) 
!    K. E. Hillstrom, 
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, 
!    May 1969.
! 
!    # 3) 
!    Hart, Et. Al., 
!    Computer Approximations, 
!    Wiley and sons, New York, 1968.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.  X must be positive.
!
!    Output, real GAMMA_LOG, the logarithm of the Gamma function of X.
!
!*******************************************************************************
!
!  Explanation of machine-depend subroutineent constants
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62E+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72E+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08E+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                           FRTBIG
!
!  CRAY-1        (S.P.)   3.13E+615
!  Cyber 180/855
!    under NOS   (S.P.)   6.44E+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   1.42E+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   2.25D+76
!  IBM 3033      (D.P.)   2.56D+18
!  VAX D-Format  (D.P.)   1.20D+9
!  VAX G-Format  (D.P.)   1.89D+76
!
  implicit none
!
  real(8), intent(in) :: x
  real(8), intent(out) :: gamma_log
  
  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728E-03, &
     8.4171387781295E-04, &
    -5.952379913043012E-04, &
     7.93650793500350248E-04, &
    -2.777777777777681622553E-03, &
     8.333333333333333331554247E-02, &
     5.7083835261e-03 /)
  real(8) :: corr
  real(8), parameter :: d1 = - 5.772156649015328605195174E-01
  real(8), parameter :: d2 =   4.227843350984671393993777E-01
  real(8), parameter :: d4 =   1.791759469228055000094023E+00
  integer i
  real(8), parameter :: frtbig = 1.42E+09

  real(8), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888E+00, &
    2.018112620856775083915565E+02, &
    2.290838373831346393026739E+03, &
    1.131967205903380828685045E+04, &
    2.855724635671635335736389E+04, &
    3.848496228443793359990269E+04, &
    2.637748787624195437963534E+04, &
    7.225813979700288197698961E+03 /)
  real(8), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064E+00, &
    5.424138599891070494101986E+02, &
    1.550693864978364947665077E+04, &
    1.847932904445632425417223E+05, &
    1.088204769468828767498470E+06, &
    3.338152967987029735917223E+06, &
    5.106661678927352456275255E+06, &
    3.074109054850539556250927E+06 /)
  real(8), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062E+04, &
    2.426813369486704502836312E+06, &
    1.214755574045093227939592E+08, &
    2.663432449630976949898078E+09, &
    2.940378956634553899906876E+10, &
    1.702665737765398868392998E+11, &
    4.926125793377430887588120E+11, &
    5.606251856223951465078242E+11 /)
  real(8), parameter :: pnt68 = 0.6796875E+00
  real(8), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036E+01, &
    1.113332393857199323513008E+03, &
    7.738757056935398733233834E+03, &
    2.763987074403340708898585E+04, &
    5.499310206226157329794414E+04, &
    6.161122180066002127833352E+04, &
    3.635127591501940507276287E+04, &
    8.785536302431013170870835E+03 /)
  real(8), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942E+02, &
    7.765049321445005871323047E+03, &
    1.331903827966074194402448E+05, &
    1.136705821321969608938755E+06, &
    5.267964117437946917577538E+06, &
    1.346701454311101692290052E+07, &
    1.782736530353274213975932E+07, &
    9.533095591844353613395747E+06 /)
  real(8), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843E+03, &
    6.393885654300092398984238E+05, &
    4.135599930241388052042842E+07, &
    1.120872109616147941376570E+09, &
    1.488613728678813811542398E+10, &
    1.016803586272438228077304E+11, &
    3.417476345507377132798597E+11, &
    4.463158187419713286462081E+11 /)
  real(8)::  res
  real(8), parameter :: sqrtpi = 0.9189385332046727417803297E+00
  real(8), parameter :: xbig = 4.08E+36
  real(8) ::xden
  real(8) ::xm1
  real(8) ::xm2
  real(8) ::xm4
  real(8) ::xnum
  real(8) :: xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0E+00 .or. x > xbig ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = - log ( x )

  else if ( x <= 1.5E+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0E+00
      xm1 = ( x - 0.5E+00 ) - 0.5E+00
    end if

    if ( x <= 0.5E+00 .or. x >= pnt68 ) then

      xden = 1.0E+00
      xnum = 0.0E+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5E+00 ) - 0.5E+00
      xden = 1.0E+00
      xnum = 0.0E+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0E+00 ) then

    xm2 = x - 2.0E+00
    xden = 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0E+00 ) then

    xm4 = x - 4.0E+00
    xden = - 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0E+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5E+00 * corr
    res = res + x * ( corr - 1.0E+00 )

  end if

  gamma_log = res

  return
end subroutine gamma_logfcn


subroutine gammafcn ( x, gamma )
!
!*******************************************************************************
!
!! GAMMA calculates the Gamma function for a real argument X.
!
!
!  Definition:
!
!    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
!
!  Recursion:
!
!    GAMMA(X+1) = X * GAMMA(X)
!
!  Special values:
!
!    GAMMA(0.5) = SQRT(PI)
!    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for X .GE. 12 are from reference 2.
!    The accuracy achieved depend subroutines on the arithmetic system, the
!    compiler, the intrinsic functions, and proper selection of the
!    machine-depend subroutineent constants.
!
!  Machine-depend subroutineent constants:
!
!    BETA: radix for the floating-point representation.
!    MAXEXP: the smallest positive power of BETA that overflows.
!    XBIG: the largest argument for which GAMMA(X) is representable
!      in the machine, i.e., the solution to the equation
!      GAMMA(XBIG) = BETA**MAXEXP.
!    XMININ: the smallest positive floating-point number such that
!      1/XMININ is machine representable.
!
!    Approximate values for some important machines are:
!
!                               BETA       MAXEXP        XBIG
!
!    CRAY-1         (S.P.)        2         8191        966.961
!    Cyber 180/855
!      under NOS    (S.P.)        2         1070        177.803
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)        2          128        35.040
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)        2         1024        171.624
!    IBM 3033       (D.P.)       16           63        57.574
!    VAX D-Format   (D.P.)        2          127        34.844
!    VAX G-Format   (D.P.)        2         1023        171.489
!
!                              XMININ
!
!    CRAY-1         (S.P.)   1.84E-2466
!    Cyber 180/855
!      under NOS    (S.P.)   3.14E-294
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)   1.18E-38
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)   2.23D-308
!    IBM 3033       (D.P.)   1.39D-76
!    VAX D-Format   (D.P.)   5.88D-39
!    VAX G-Format   (D.P.)   1.12D-308
!
!  Reference: 
!
!    W J Cody,
!    "An Overview of Software Development for Special Functions", 
!    Lecture Notes in Mathematics, 506, 
!    Numerical Analysis Dundee, 1975, 
!    G. A. Watson (ed.),
!    Springer Verlag, Berlin, 1976.
!
!    Hart et al,
!    Computer Approximations, 
!    Wiley and sons, New York, 1968.
!
!  Author: 
!
!    W. J. Cody and L. Stoltz,
!    Applied Mathematics Division,
!    Argonne National Laboratory,
!    Argonne, Illinois, 60439.
!
!  Parameters:
!
!    Input, real X, the argument of the function.
!
!    Output, real GAMMA, the value of the function. 
!    The computation is believed to be free of underflow and overflow.
!
  implicit none
  
  real(8), intent(in) :: x
  real(8), intent(out) :: gamma
  
  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728E-03, &
     8.4171387781295E-04, &
    -5.952379913043012E-04, &
     7.93650793500350248E-04, &
    -2.777777777777681622553E-03, &
     8.333333333333333331554247E-02, &
     5.7083835261E-03 /)
  real(8) :: fact
  
  integer :: i
  integer :: n
  real, parameter, dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811E+00, &
     2.47656508055759199108314E+01, &
    -3.79804256470945635097577E+02, &
     6.29331155312818442661052E+02, &
     8.66966202790413211295064E+02, &
    -3.14512729688483675254357E+04, &
    -3.61444134186911729807069E+04, &
     6.64561438202405440627855E+04 /)
  logical parity
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real, parameter, dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353E+01, &
     3.15350626979604161529144E+02, &
    -1.01515636749021914166146E+03, &
    -3.10777167157231109440444E+03, &
     2.25381184209801510330112E+04, &
     4.75584627752788110767815E+03, &
    -1.34659959864969306392456E+05, &
    -1.15132259675553483497211E+05 /)
  real, parameter :: sqrtpi = 0.9189385332046727417803297E+00
  real(8) :: sum2
  real, parameter :: xbig = 35.040E+00
  real(8) :: xden
  real, parameter :: xminin = 1.18E-38
  real(8) :: xnum
  real(8) :: y
  real(8) :: y1
  real(8) :: ysq
  real(8) :: z
!
  parity = .false.
  fact = 1.0E+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0E+00 ) then

    y = - x
    y1 = aint ( y )
    gamma = y - y1

    if ( gamma /= 0.0E+00 ) then

      if ( y1 /= aint ( y1 * 0.5E+00 ) * 2.0E+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * gamma )
      y = y + 1.0E+00

    else

      gamma = huge ( gamma )
      return

    end if

  end if
!
!  Argument < EPS
!
  if ( y < epsilon ( y ) ) then

    if ( y >= xminin ) then
      gamma = 1.0E+00 / y
    else
      gamma = huge ( gamma )
      return
    end if

  else if ( y < 12.0E+00 ) then

    y1 = y
!
!  0.0E+00 < argument < 1.0E+00
!
    if ( y < 1.0E+00 ) then
      z = y
      y = y + 1.0E+00
!
!  1.0E+00 < argument < 12.0E+00, reduce argument if necessary.
!
    else
      n = int ( y ) - 1
      y = y - real ( n )
      z = y - 1.0E+00
    end if
!
!  Evaluate approximation for 1.0E+00 < argument < 2.0.
!
    xnum = 0.0E+00
    xden = 1.0E+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    gamma = xnum / xden + 1.0E+00
!
!  Adjust result for case  0.0E+00 < argument < 1.0.
!
    if ( y1 < y ) then
      gamma = gamma / y1
!
!  Adjust result for case  2.0E+00 < argument < 12.0.
!
    else if ( y1 > y ) then

      do i = 1, n
        gamma = gamma * y
        y = y + 1.0E+00
      end do

    end if
!
!  Evaluate for 12 <= argument.
!
  else

    if ( y <= xbig ) then

      ysq = y**2
      sum2 = c(7)
      do i = 1, 6
        sum2 = sum2 / ysq + c(i)
      end do
      sum2 = sum2 / y - y + sqrtpi
      sum2 = sum2 + ( y - 0.5E+00 ) * log ( y )
      gamma = exp ( sum2 )

    else

      gamma = huge ( gamma )
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    gamma = - gamma
  end if

  if ( fact /= 1.0E+00 ) then
    gamma = fact / gamma
  end if

  return
end subroutine gammafcn

subroutine betafcn( a, b, beta )
!
!*******************************************************************************
!
!! BETA returns the value of the Beta function.
!
!
!  Formula:
!
!    BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
!              = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
!
!  Modified:
!
!    10 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the function.
!    0.0E+00 < A,
!    0.0E+00 < B.
!
!    Output, real BETA, the value of the function.
!
  implicit none
  real(8), intent(in) :: a
  real(8), intent(in) :: b
  real(8), intent(out) :: beta
  real(8) :: gamma_loga
  real(8) :: gamma_logb
  real(8) :: gamma_logapb
  if ( a <= 0.0E+00 .or. b <= 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA - Fatal error!'
    write ( *, * ) '  Both A and B must be greater than 0.'
    stop
  end if
  call gamma_logfcn(a,gamma_loga)
  call gamma_logfcn(b,gamma_logb)
  call gamma_logfcn((a+b),gamma_logapb)
  beta = exp ( gamma_loga + gamma_logb - gamma_logapb )
  return
end subroutine betafcn
end module pdf_fcns