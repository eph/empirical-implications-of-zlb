!------------------------------------------------------------------------------------------------------------------
! MODULE: class_model
!> @version 5.2
!> @brief Specifies solution details of NK model along the lines of CEE or SW.
!-----------------------------------------------------------------------------------------------------------------
module class_model

use model_details, only: initializesolution, decr, solutiondetails, soldeallocate, model_details_ss => steadystate
use get_decisionrule, only: nonlinearsolver
use get_decisionrule_parallel, only: nonlinearsolver_parallel
implicit none

double precision, parameter :: M_PI = 3.14159265358979323846d0

type model
   character(len=144) :: name = 'ghlss_v5_2'

   character(len=144) :: datafile = 'data/glss_data.txt'
   integer :: nobs = 7, T = 125, neps = 6
   double precision, allocatable :: yy(:,:) 
   double precision, allocatable :: HH(:,:)


   type(solutiondetails) :: solution
   double precision, allocatable :: params(:)

   integer :: npara, nvars, nexog

 contains
   procedure :: describe_params
   procedure :: describe
   generic :: solve => solve_serial, solve_parallel
   procedure :: solve_serial
   procedure :: solve_parallel
   procedure :: simulate_modeldata
   procedure :: simulate_modelirfs
   procedure :: cleanup


   procedure :: load_data
   procedure :: steadystate
   procedure :: g 
   procedure :: logpdfy_kernel
   procedure :: pdfy

   procedure :: pr


end type model

interface model
   module procedure new_model
end interface model

contains

  subroutine load_data(m)
    class(Model) :: m
    integer :: i

    open(1, file=m%datafile, action='read')

    do i = 1, m%T
       read(1,*) m%yy(:,i)
    end do

    close(1)


  end subroutine load_data

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: describe
!> @date 4-14-16
!> @brief Report some of the details about the solution method.  
!---------------------------------------------------------------------------------------------------------  
subroutine describe(m)
  implicit none
  class(model) :: m

  write(*,*) '- - - - - - - - - - - - '
  print *, 'Model Name: ', trim(m%name)
  print *, 'number of parameters: ', m%solution%poly%nparams
  write(*,*) 'number of potential shocks = ', m%solution%poly%nexog
  write(*,*) 'number of active shocks = ', m%solution%poly%nexogshock+m%solution%poly%nexogcont
  write(*,*) 'number of shocks in polynomial = ', m%solution%poly%nexogcont
  write(*,*) 'number of grid points used for finite-element part of decr = ', m%solution%poly%ns
  write(*,*) 'number of grid points used for poly part of decr = ', m%solution%poly%ngrid
  write(*,*) 'number of quadrature points = ', m%solution%poly%nquad
  write(*,*) '- - - - - - - - - - - - '
  
end subroutine describe

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: describe_params
!> @date 4-14-16
!> @brief Report parameter values.  
!---------------------------------------------------------------------------------------------------------  
subroutine describe_params(m)
  implicit none
  class(model) :: m

  integer :: i

  print *, 'model Name: ', trim(m%name)
  print *, 'number of parameters: ', m%solution%poly%nparams

  write(*,*) '- - - - - - - - - - - - '
  do i = 1,m%solution%poly%nparams 
     write(*,*) 'm%params(',i,') = ', m%params(i)
  end do
  write(*,*) '---------------------------------'

end subroutine describe_params

!--------------------------------------------------------------------------------------------------------
! FUNCTION: new_model
!> @date 11-16-16
!> @brief Set some parameter values for the model and choose details about the solution.  
!---------------------------------------------------------------------------------------------------------  
type(model) function new_model(zlbswitch) result(m)

  use utils, only: read_matrix
  implicit none
  
  logical, intent(in), optional :: zlbswitch

 !local variables
  character (LEN = 200) :: inputfile
  
  if (present(zlbswitch)) then
     m%solution%poly%zlbswitch = zlbswitch
  else
     m%solution%poly%zlbswitch = .true.
  end if

  if (zlbswitch .eqv. .false.) then
     m%name = 'ghlss_v5_2_unc'
  end if

  m%solution%poly%nparams = 43
  m%solution%poly%nexog = 6
  m%solution%poly%nexogcont = 0 
  !m%solution%poly%nexogcont = 1  !adds level technology shock into polynomial 
  m%solution%poly%nvars = 22
  m%solution%poly%nmsv = 7
  m%solution%poly%nfunc = 7
  m%solution%poly%nindplus = 1  !adds higher order terms for investment in polynomial 
  
  m%npara = m%solution%poly%nparams
  m%nvars = m%solution%poly%nvars+m%solution%poly%nexog
  m%nexog = m%solution%poly%nexog


  !allocate array for model's parameters
  allocate(m%params(m%solution%poly%nparams))

  allocate(m%solution%poly%indplus(m%solution%poly%nindplus))
  if (m%solution%poly%nindplus .eq. 1) m%solution%poly%indplus(1) = 3  !go higher order with respect to lagged investment in each polynomial

  !set number of grid points per shock
   allocate(m%solution%poly%nshockgrid(m%solution%poly%nexog))
   m%solution%poly%nshockgrid(1) = 7 !shock to safe asset
   m%solution%poly%nshockgrid(2) = 3 !MEI shock
   m%solution%poly%nshockgrid(3) = 3 !tech
   m%solution%poly%nshockgrid(4) = 3 !r shock
   m%solution%poly%nshockgrid(5) = 3 !g shock

   ! m%solution%poly%nshockgrid(2) = 1 !MEI shock
   ! m%solution%poly%nshockgrid(3) = 1 !tech
   ! m%solution%poly%nshockgrid(4) = 1 !r shock
   ! m%solution%poly%nshockgrid(5) = 1 !g shock
   if (m%solution%poly%nexogcont .eq. 0) m%solution%poly%nshockgrid(6) = 1  !level technology shock

   call  initializesolution(m%solution)

    allocate(m%yy(m%nobs,m%T), m%HH(m%nobs,m%nobs))
    
    call m%load_data()

end function new_model

!--------------------------------------------------------------------------------------------------------
! FUNCTION: solve
!> @date 4-14-16
!> @brief Solve the model using the serial code and check for convergence.
!---------------------------------------------------------------------------------------------------------  
function solve_serial(m, params) result(convergence)

  implicit none
  class(model) :: m

  double precision, intent(in) :: params(m%solution%poly%nparams)
  logical :: convergence

  call nonlinearsolver(params, m%solution, convergence)
  m%params = params

end function solve_serial

function steadystate(self, params) result(ss)
  class(Model) :: self

  double precision, intent(in) :: params(self%npara)
  double precision :: ss(self%nvars)

  integer :: i 
  ss = 0.0d0

  ss = model_details_ss(params,self%nvars,self%npara)

  do i = 1,self%nvars
     if (i <= self%solution%poly%nmsv) then 
        ss(i) = exp(ss(i))
     else 
        ss(i) = ss(i)
     end if


  end do

end function steadystate


!--------------------------------------------------------------------------------------------------------
! FUNCTION: solve_parallel
!> @date 4-18-16
!> @brief Solve the model using the parallelized code and check for convergence.
!---------------------------------------------------------------------------------------------------------  
function solve_parallel(m, params, nproc, rank) result(convergence)

  implicit none
  class(Model) :: m

  double precision, intent(in) :: params(m%solution%poly%nparams)
  integer, intent(in) :: nproc, rank
  logical :: convergence

  logical :: inbounds
  if (.not. inbounds(m%npara, params)) then
     convergence = .false.
     return
  end if
  call nonlinearsolver_parallel(nproc, rank, params, m%solution, convergence)

  m%params = params

end function solve_parallel

!--------------------------------------------------------------------------------------------------------
!SUBROUTINE: simulate_modeldata
!> @date 11-16-16
!> @brief Simulate data from the nonlinear or linear model.
!> param[in] nonlinearswitch if true, simulate from nonlinear model; otherwise linear model
!---------------------------------------------------------------------------------------------------------  
subroutine simulate_modeldata(m,capt,nonlinearswitch,modeldata,seed)

  use simulate_model, only: simulate_data
  implicit none
  class(Model) :: m

  integer, intent(in) :: capt
  logical, intent(in) :: nonlinearswitch
  integer, intent(in) :: seed
  double precision, intent(out) :: modeldata(m%solution%poly%nvars+2*m%solution%poly%nexog,capt)
  
  call simulate_data(capt,m%params,m%solution%poly,m%solution%linsol,m%solution%alphacoeff,modeldata,nonlinearswitch,seed)
  
end subroutine simulate_modeldata

!--------------------------------------------------------------------------------------------------------
!SUBROUTINE: simulate_modelirfs
!> @date 11-17-16
!> @brief Simulate irfs from the nonlinear and linear models.
!---------------------------------------------------------------------------------------------------------  
subroutine simulate_modelirfs(m,capt,nsim,shockindex,endogirf,linirf,euler_errors,neulererrors,scaleshock_opt)

  use simulate_model, only: simulate_irfs
  implicit none
  class(Model) :: m

  integer, intent(in) :: capt
  integer, intent(in) :: nsim
  integer, intent(in) :: shockindex
  integer, intent(in) :: neulererrors
  double precision, intent(in), optional :: scaleshock_opt
  double precision, intent(out) :: endogirf(m%solution%poly%nvars+m%solution%poly%nexog+2,capt)
  double precision, intent(out) :: linirf(m%solution%poly%nvars+m%solution%poly%nexog+2,capt)
  double precision, intent(out) :: euler_errors(2*neulererrors,capt)

  double precision :: scaleshock
  double precision :: endogvarbas0(m%solution%poly%nvars+m%solution%poly%nexog)
  double precision :: endogvarshk0(m%solution%poly%nvars+m%solution%poly%nexog)
  double precision :: innov0(m%solution%poly%nexog)
  double precision :: premiumirf(2,capt)
  
  if (present(scaleshock_opt)) then 
     scaleshock = scaleshock_opt
  else
     scaleshock = 1.d0
  end if


  endogvarbas0 = 0.0d0
  endogvarbas0(1:m%solution%poly%nvars) = exp(m%solution%poly%endogsteady(1:m%solution%poly%nvars))
  endogvarshk0 = endogvarbas0
  innov0 = 0.0d0
  
  !shock to liquidity    
  if (shockindex .eq. 1) then
     endogvarshk0(m%solution%poly%nvars+shockindex) = scaleshock*m%params(30)/(1.0d0-m%params(29)**2)    
  end if

  !shock to MEI -- note this is negative
  if (shockindex .eq. 2) then
     endogvarshk0(m%solution%poly%nvars+shockindex) = -scaleshock*m%params(28)/(1.0d0-m%params(27)**2)      
  end if
  endogvarbas0 = endogvarshk0  !change initial conditions of shock and baseline dset
  innov0(shockindex) = scaleshock

  call simulate_irfs(capt,nsim,m%params,m%solution%poly,m%solution%alphacoeff,m%solution%linsol,endogvarbas0,endogvarshk0,innov0,endogirf,linirf,&
       euler_errors,neulererrors,premiumirf)
    
end subroutine simulate_modelirfs



  function g(m, states_old, shocks_curr,regime_curr) result(states_curr)
    class(Model) :: m

    double precision, intent(in) :: states_old(m%nvars), shocks_curr(m%nexog)
    integer :: regime_curr(2)
    double precision :: states_curr(m%nvars)

    call decr(states_old, shocks_curr, m%params, &
         m%solution%poly, m%solution%alphacoeff, states_curr)

  end function g

  function pdfy(m, t, states_new, states_old, para) result(pdf)
    class(Model) :: m 

    integer, intent(in) :: t

    double precision, intent(in) :: states_new(m%nvars), states_old(m%nvars)
    double precision, intent(in) :: para(m%npara)
    double precision :: pdf

    double precision :: gz, zhat, observables(m%nobs)

    integer :: ngdp, ncons, ninv, nlab, nwage, ninfl, nnomr,ntechshk
    integer :: nobsgdp, nobscons, nobsinv, nobslab, nobswage, nobsinfl, nobsnomr
    integer :: nmegdp, nmecons, nmeinv, nmelab, nmewage, nmeinfl, nmenomr
    integer :: j
    real(8), parameter :: const_pi = 3.14159265358979323846d0
    real(8) :: temp
    integer :: nmestart

    ! order of observables: dy  dc  di  h dwp pi  ffr
    nobsgdp     = 1
    nobscons    = 2
    nobsinv     = 3
    nobsinfl    = 6
    nobsnomr    = 7

    nmestart = 36
    nmegdp  = nmestart + nobsgdp
    nmecons = nmestart + nobscons
    nmeinv  = nmestart + nobsinv
    !nmelab  = nmestart + nobslab
    !nmewage = nmestart + nobswage
    nmeinfl = nmestart + nobsinfl
    nmenomr = nmestart + nobsnomr

    observables = m%yy(:,t)

    !all variables are in levels except for shocks
    !shocks are in log-level deviations from SS
    !(1) cap, (2) cc, (3) inv, (4) rw, (5) notr, (6) dp, (7) gdp, (8) xhp, (9) nomr, (10) lam, (11) qq, (12) lab, (13) util, (14) mc, (15) rentalk, (16) muc
    !(17) vi, (18) vp, (19) vw, (20) dw, (21) bc, (22) bi, (23) liqshk, (24) invshk, (25) techshk, (26) intshk, (27) gshk, (28) ashk


    ngdp    = 7
    ncons   = 2
    ninv    = 3
    nlab    = 12
    nwage   = 20
    ninfl   = 6
    nnomr    = 9

    ntechshk = 25
    gz = log(para(3))  + states_new(ntechshk)



    !Weight update
    pdf = 1.0d0

    !Model specific (some of these for short run code clarity)

    ! log gdp growth
    !temp = states_new(ngdp) - states_old(ngdp) + gz 
    temp = log(states_new(ngdp)/states_old(ngdp))   + gz 
    temp = 1.0d0/(para(nmegdp))  * ( observables(nobsgdp) - temp )
    pdf = pdf*1.0d0/(para(nmegdp) )*exp(-0.5d0*temp*temp)
    
    ! log consumption growth
    temp = log(states_new(ncons)/states_old(ncons)) + gz 
    temp = 1.0d0/(para(nmecons))  * ( observables(nobscons) - temp )
    pdf = pdf*1.0d0/(para(nmecons) )*exp(-0.5d0*temp*temp)
    
    !! log investment growth
    temp = log(states_new(ninv)/states_old(ninv))  + gz 
    temp = 1.0d0/(para(nmeinv))  * ( observables(nobsinv) - temp )
    pdf = pdf*1.0d0/(para(nmeinv) )*exp(-0.5d0*temp*temp)
    
    ! log inflation
    temp = log(states_new(ninfl))
    temp = 1.0d0/para(nmeinfl)  * ( observables(nobsinfl) - temp )
    pdf = pdf*1.0d0/para(nmeinfl)*exp(-0.5d0*temp*temp)
    
    ! log nominal rate
    temp =  log(states_new(nnomr))
    temp = 1.0d0/para(nmenomr)  * ( observables(nobsnomr) - temp )
    pdf = pdf*1.0d0/para(nmenomr)*exp(-0.5d0*temp*temp)
    
    pdf = pdf*1.0d0/sqrt(2.0d0*const_pi)**(5)


  end function pdfy

  function logpdfy_kernel(m, t, states_new, states_old, para) result(pdf)
    class(Model) :: m 

    integer, intent(in) :: t

    double precision, intent(in) :: states_new(m%nvars), states_old(m%nvars)
    double precision, intent(in) :: para(m%npara)
    double precision :: pdf

    double precision :: gz, zhat, observables(m%nobs)

    integer :: ngdp, ncons, ninv, nlab, nwage, ninfl, nnomr,ntechshk
    integer :: nobsgdp, nobscons, nobsinv, nobslab, nobswage, nobsinfl, nobsnomr
    integer :: nmegdp, nmecons, nmeinv, nmelab, nmewage, nmeinfl, nmenomr
    integer :: j
    real(8), parameter :: const_pi = 3.14159265358979323846d0
    real(8) :: temp
    integer :: nmestart

    ! order of observables: dy  dc  di  h dwp pi  ffr
    nobsgdp     = 1
    nobscons    = 2
    nobsinv     = 3
    nobsinfl    = 6
    nobsnomr    = 7

    nmestart = 36
    nmegdp  = nmestart + nobsgdp
    nmecons = nmestart + nobscons
    nmeinv  = nmestart + nobsinv
    !nmelab  = nmestart + nobslab
    !nmewage = nmestart + nobswage
    nmeinfl = nmestart + nobsinfl
    nmenomr = nmestart + nobsnomr

    observables = m%yy(:,t)

    !all variables are in levels except for shocks
    !shocks are in log-level deviations from SS
    !(1) cap, (2) cc, (3) inv, (4) rw, (5) notr, (6) dp, (7) gdp, (8) xhp, (9) nomr, (10) lam, (11) qq, 
    !(12) lab, (13) util, (14) mc, (15) rentalk, (16) muc, (17) vi, (18) vp, (19) vw, (20) dw, (21) bc, (22) bi,
    !(23) liqshk, (24) invshk, (25) techshk, (26) rrshk, (27) gshk, (28) elastshk, (29) elastwshk
    !endogvar = (/cap, cc, inv, rw, notr, dp, gdp, xhp, nomr, lam, &
    !     qq, lab, util, mc, rentalk, muc, vi, vp, vw, dw, bc, bi, currentshockvalues/) 

    ngdp    = 7
    ncons   = 2
    ninv    = 3
    nlab    = 12
    nwage   = 20
    ninfl   = 6
    nnomr    = 9

    ntechshk = 25
    gz = log(para(3))  + states_new(ntechshk)



    !Weight update
    pdf = 0.0d0

    !Model specific (some of these for short run code clarity)

    ! log gdp growth
    !temp = states_new(ngdp) - states_old(ngdp) + gz 
    temp = log(states_new(ngdp)/states_old(ngdp))   + gz 
    temp = 1.0d0/(para(nmegdp))  * ( observables(nobsgdp) - temp )
    !pdf = pdf*1.0d0/(para(nmegdp) )*exp(-0.5d0*temp*temp)
    pdf = pdf - 0.5d0*temp**2
    ! log consumption growth
    temp = log(states_new(ncons)/states_old(ncons)) + gz 
    temp = 1.0d0/(para(nmecons))  * ( observables(nobscons) - temp )
    !pdf = pdf*1.0d0/(para(nmecons) )*exp(-0.5d0*temp*temp)
    pdf = pdf - 0.5d0*temp**2
    !! log investment growth
    temp = log(states_new(ninv)/states_old(ninv))  + gz 
    temp = 1.0d0/(para(nmeinv))  * ( observables(nobsinv) - temp )
    !pdf = pdf*1.0d0/(para(nmeinv) )*exp(-0.5d0*temp*temp)
    pdf = pdf - 0.5d0*temp**2
    ! log inflation
    temp = log(states_new(ninfl))
    temp = 1.0d0/para(nmeinfl)  * ( observables(nobsinfl) - temp )
    !pdf = pdf*1.0d0/para(nmeinfl)*exp(-0.5d0*temp*temp)
    pdf = pdf - 0.5d0*temp**2
    ! log nominal rate
    temp =  log(states_new(nnomr))
    temp = 1.0d0/para(nmenomr)  * ( observables(nobsnomr) - temp )
    !pdf = pdf*1.0d0/para(nmenomr)*exp(-0.5d0*temp*temp)
    pdf = pdf - 0.5d0*temp**2
    !pdf = pdf*1.0d0/sqrt(2.0d0*const_pi)**(5)
    !pdf = log(pdf)

    if ((isnan(pdf)) .or. (pdf < -10000000000.0d0)) then
       print*,states_new,states_old
       print*,pdf
       print*,observables
       print*,'end'
       stop
    end if

  end function logpdfy_kernel


  double precision function pr(m, para) result(pr0)

    class(Model) :: m
    double precision, intent(in) :: para(m%npara)

    call priorfcn(m%npara, para, pr0)

  end function pr


!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: cleanup
!> @date 4-14-16
!> @brief Deallocate arrays.
!---------------------------------------------------------------------------------------------------------  
subroutine cleanup(m)

  implicit none
  class(model) :: m

  deallocate(m%params)
  call soldeallocate(m%solution)
end subroutine cleanup

end module class_model
