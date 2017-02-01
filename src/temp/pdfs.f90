module pdfs

! Module cotains the following pdfs
!
! beta_pdf ( x, a, b, pdf )
! gamma_pdf ( x, a, b, pdf )
! invgamma_pdf( x, a, b, pdf )
! gamma_pdf2 ( x, a, b, pdf )
! invgamma_pdf2( x, a, b, pdf )
! normal_pdf( x, a, b, pdf )
! uniform_pdf( x, a, b, pdf )

contains

subroutine gamma_converter(mu, sigma, a, b)
implicit none
real(8), intent(in) :: mu
real(8), intent(in) :: sigma
real(8), intent(out) :: a
real(8), intent(out) :: b
a = (mu/sigma)**2
b = mu/a
end subroutine

subroutine gamma_converter_mode(mode, sigma, a, b)
implicit none
real(8), intent(in) :: mode
real(8), intent(in) :: sigma
real(8), intent(out) :: a
real(8), intent(out) :: b
a = (1.0-mode/sigma)**(-1)
b = mode/(a-1.d0)
end subroutine

subroutine beta_converter(mu,sigma,a,b)
  implicit none

  real(8), intent(in) :: mu, sigma
  real(8), intent(out) :: a, b

  a = (1-mu)*mu**2/sigma**2 - mu
  b = a*(1/mu - 1)
end subroutine 


subroutine beta_pdf ( x, a, b, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: x
real(8), intent(in) :: a
real(8), intent(in) :: b
real(8), intent(out) :: pdf
real(8) :: beta
if ( x < 0.0E+00 .or. x > 1.0E+00 ) then
    pdf = 0.0E+00
else
    call betafcn(a,b,beta)
    pdf = x**( a - 1.0E+00 ) * ( 1.0E+00 - x )**( b - 1.0E+00 ) / beta
end if
end subroutine

subroutine gamma_pdf( x, a, b, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: x
real(8), intent(in) :: a
real(8), intent(in) :: b
real(8), intent(out) :: pdf
real(8) :: gamma_a
if ( x <= 0 ) then
    pdf = 0.0E+00
else
    call gammafcn(a,gamma_a)
    pdf = ( 1.0d0/ (gamma_a*(b**a) ) ) * x**( a - 1.0d0 ) * exp ( -x/b )
end if
end subroutine


subroutine gamma_pdf_2( x, v, s, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: x
real(8), intent(in) :: v
real(8), intent(in) :: s
real(8), intent(out) :: pdf
real(8) :: gamma_a
real(8) :: a
real(8) :: b
a = v/2.0d0
b = 2.0d0/s
if ( x <= 0 ) then
    pdf = 0.0E+00
else
    call gammafcn(a,gamma_a)
    pdf = ( 1.0d0/ (gamma_a*(b**a) ) ) * x**( a - 1.0d0 ) * exp ( -x/b )
end if
end subroutine

subroutine gamma_pdf_1( x, v, s, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: x
real(8), intent(in) :: v
real(8), intent(in) :: s
real(8), intent(out) :: pdf
real(8) :: gamma_a
real(8) :: a
real(8) :: b
a = v/2.0d0
b = 2.0d0/s

if ( x <= 0 ) then
    pdf = 0.0E+00
else
    call gammafcn(a,gamma_a)
    pdf = ( 2.0d0/ (gamma_a*(b**a) ) ) * x**( v - 1.0d0 ) * exp ( -0.5d0*s*x**2 )
end if
end subroutine


subroutine invgamma_pdf( x, a, b, pdf )
use pdf_fcns
implicit none
real(8),intent(in) :: a
real(8),intent(in) :: b
real(8),intent(in) :: x
real(8),intent(out) :: pdf
real(8) :: gamma_a
if ( x <= 0.0d0 ) then
    pdf = 0.0d0
  else
    call gammafcn(a,gamma_a)
    pdf = (b**a) / gamma_a * (1.0d0/x)**( a + 1.d0 )* exp ( -b/x )
end if
end subroutine


subroutine invgamma_pdf_dynare( x, a, b, pdf )
use pdf_fcns
implicit none
real(8),intent(in) :: a
real(8),intent(in) :: b
real(8),intent(in) :: x
real(8),intent(out) :: pdf
real(8) :: gamma_a
if ( x <= 0.0d0 ) then
    pdf = 0.0d0
  else
    call gammafcn(a,gamma_a)
    pdf = 1.0d0 / (gamma_a * b**a) * (x)**(-a -1.d0 ) * exp ( -1.d0/(x*b) )
end if
end subroutine

subroutine invgamma_pdf_2( x, v, s, pdf )
use pdf_fcns
implicit none
real(8),intent(in) :: v
real(8),intent(in) :: s
real(8),intent(in) :: x
real(8),intent(out) :: pdf
real(8) :: gamma_a
real(8) :: a
real(8) :: b
a = v/2.0d0
b = 2.0d0/s
if ( x <= 0.0d0 ) then
    pdf = 0.0d0
  else
    call gammafcn(a,gamma_a)
    pdf = 1 / (gamma_a * b**a) * (x)**(-a - 1.d0 )* exp ( 1.0d0/(x*b) )
end if
end subroutine

subroutine invgamma_pdf_1( x, v, s, pdf )
use pdf_fcns
implicit none
real(8),intent(in) :: v
real(8),intent(in) :: s
real(8),intent(in) :: x
real(8),intent(out) :: pdf
real(8) :: gamma_a
real(8) :: a
real(8) :: b
a = v/2.0d0
b = 2.0d0/s
pdf = 0.0d0
if ( x <= 0.0d0 ) then
    pdf = 0.0d0
  else
    call gammafcn(a,gamma_a)
    pdf = 2.0d0 / (gamma_a * b**a) * (x)**( -v - 1.d0 )* exp (-0.5d0*s/x**2 )
end if
end subroutine

subroutine normal_pdf ( x, a, b, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: a
real(8), intent(in) :: b
real(8), intent(in) :: x
real(8), intent(out) :: pdf
real(8), parameter :: pi = 3.14159265358979323846264338327950288419716939937510
real(8) :: y
y = ( x - a ) / b
pdf = exp ( - 0.5d0 * y**2 )  / ( b * sqrt ( 2.0d0 * pi ) )
end subroutine


subroutine uniform_pdf ( x, a, b, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: a
real(8), intent(in) :: b
real(8), intent(in) :: x
real(8), intent(out) :: pdf
pdf = 1.0d0/(b-a)
end subroutine


subroutine rayleigh_pdf ( x, a, pdf )
use pdf_fcns
implicit none
real(8), intent(in) :: a
real(8), intent(in) :: x
real(8), intent(out) :: pdf
real(8) :: y
y = x/a
pdf = exp ( - 0.5d0 * (y)**2 ) * y / a
end subroutine

function logigpdf(x,a,b)
  use pdf_fcns
  real(8), intent(in) :: x, a, b
  real(8) :: logigpdf

  real(8) :: lg
  
  call gamma_logfcn(b/2.0d0,lg)
  logigpdf = log(2.0d0) - lg + b/2.0d0*log(b*a**2/2.0d0) &
       -(b+1.0d0)/2.0d0*log(x**2) - b*a**2.0/(2.0d0*x**2)

end function logigpdf


end module pdfs
