module fortress_prior_t

  use RandomNumber_class, only : RandomNumber
  use logbeta, only : betaln, gamln

  implicit none

  integer, parameter :: PARA_BETA = 1
  integer, parameter :: PARA_GAMMA = 2
  integer, parameter :: PARA_NORMAL = 3
  integer, parameter :: PARA_INVGAMMA = 4
  integer, parameter :: PARA_UNIFORM = 5
  integer, parameter :: PARA_FIXED = 6

  integer, parameter :: PRIOR_FILE_UNIT = 26

  double precision, parameter :: M_PI = 3.14159265358979d0

  type, abstract :: fortress_abstract_prior
     integer :: npara
   contains
     !procedure(rvs_interface), deferred :: rvs
     procedure(logpdf_interface), deferred :: logpdf
  end type fortress_abstract_prior

  interface fortress_abstract_prior

     function rvs_interface(self,nsim) result(rvs)
       import fortress_abstract_prior
       class(fortress_abstract_prior), intent(inout) :: self
       integer, intent(in) :: nsim
       double precision :: rvs(self%npara,nsim)
     end function rvs_interface

     function logpdf_interface(self, para) result(lpdf)
       import fortress_abstract_prior
       class(fortress_abstract_prior), intent(inout) :: self
       double precision, intent(in) :: para(self%npara)
       double precision :: lpdf
     end function logpdf_interface

  end interface fortress_abstract_prior


  type, extends(fortress_abstract_prior) :: prior

     integer :: seed = 1848

     !integer :: npara

     integer, allocatable :: ptype(:), pfix(:)
     double precision, allocatable :: pmean(:), pstdd(:), pval(:)
     double precision, allocatable :: plower(:), pupper(:)

     type(RandomNumber) :: rn

   contains
     procedure :: logpdf
     !procedure :: rvs
  end type prior

  interface prior
     module procedure new_prior
  end interface prior

contains

  type(prior) function new_prior(frankfile) result(pr)

    integer :: nlines, unit, io, i

    character(len=*), intent(in) :: frankfile


    open(PRIOR_FILE_UNIT, file=frankfile)
    nlines = 0

    do
       read(PRIOR_FILE_UNIT,*,iostat=io)
       if (io/=0) exit
       nlines = nlines + 1
    end do

    rewind(PRIOR_FILE_UNIT)

    pr%npara = nlines
    allocate(pr%ptype(pr%npara), pr%pmean(pr%npara), pr%pstdd(pr%npara), &
         pr%pfix(pr%npara), pr%pval(pr%npara), &
         pr%plower(pr%npara), pr%pupper(pr%npara))

    do i = 1, pr%npara
       read(PRIOR_FILE_UNIT,*) pr%ptype(i), pr%pmean(i), pr%pstdd(i), pr%pfix(i), pr%pval(i)
       select case( pr%ptype(i) )
       case( PARA_BETA )
          pr%plower(i) = 0.000001d0
          pr%pupper(i) = 0.999999d0
       case( PARA_GAMMA )
          pr%plower(i) = 0.000001d0
          pr%pupper(i) = 50.00000d0
       case( PARA_NORMAL )
          pr%plower(i) = -9999999.0d0
          pr%pupper(i) =  9999999.0d0
       case( PARA_INVGAMMA )
          pr%plower(i) = 0.0000001d0
          pr%pupper(i) = 50.000000d0
       case( PARA_UNIFORM )
          pr%plower(i) = pr%pmean(i)
          pr%pupper(i) = pr%pstdd(i)
       case( PARA_FIXED )
          pr%plower(i) = pr%pval(i)
          pr%pupper(i) = pr%pval(i)
       case default
          print*,'ERROR: in prior(.), parameter', i, 'has a misspecified prior'
          stop
       end select
    end do

    close(PRIOR_FILE_UNIT)
  end function new_prior

  double precision function logpdf(self, para) result(logprior)

    class(prior), intent(inout) :: self
    double precision, intent(in) :: para(self%npara)

    double precision :: a, b
    integer :: i

    logprior = 0.0d0

    associate(pmean => self%pmean, pstdd => self%pstdd )
      do i = 1, self%npara

         select case ( self%ptype(i) )
         case( PARA_BETA )
            a = (1-pmean(i))*pmean(i)**2/pstdd(i)**2 - pmean(i)
            b = a*(1/pmean(i) - 1)
            logprior = logprior + logbetapdf(para(i),a,b)
         case ( PARA_GAMMA )
            b = pstdd(i)**2/pmean(i) !theta
            a = pmean(i)/b           !k
            logprior = logprior + loggampdf(para(i),a,b)
         case ( PARA_NORMAL )
            a = pmean(i)
            b = pstdd(i)
            logprior = logprior + lognorpdf(para(i),a,b)
         case ( PARA_INVGAMMA )
            a = pmean(i)
            b = pstdd(i)
            logprior = logprior + logigpdf(para(i),a,b)
          case ( PARA_UNIFORM )
             a = pmean(i)
             b = pstdd(i)
             logprior = logprior + log(1.0d0/(b - a))
         end select
      end do
    end associate

  end function logpdf

  ! function rvs(self, size) result(parasim)

  !   class(prior), intent(inout) :: self
  !   integer, intent(in) :: size

  !   double precision :: parasim(self%npara, nsim)



  ! end function rvs

  function logbetapdf(x, a, b)

    double precision, intent(in) :: x, a, b
    double precision :: logbetapdf

    logbetapdf = -betaln(a,b) + (a-1.0d0)*log(x) + (b-1.0d0)*log(1.0-x)

  end function logbetapdf

  function loggampdf(x, a, b)

    double precision, intent(in) :: x, a, b
    double precision :: loggampdf

    loggampdf = -gamln(a) -a*log(b) + (a-1.0d0)*log(x) - x/b

  end function loggampdf

  function lognorpdf(x, a, b)

    double precision, intent(in) :: x, a, b
    double precision :: lognorpdf

    lognorpdf = -0.5d0*log(2.0d0*M_PI) - log(b) - 0.5d0*(x-a)**2/b**2

  end function lognorpdf

  function logigpdf(x,a,b)

    double precision, intent(in) :: x, a, b
    double precision :: logigpdf

    logigpdf = log(2.0d0) - gamln(b/2.0d0) + b/2.0d0*log(b*a**2/2.0d0) &
         -(b+1.0d0)/2.0d0*log(x**2) - b*a**2.0/(2.0d0*x**2)

  end function logigpdf

end module fortress_prior_t
