!------------------------------------------------------------------------------------------------------------------
! MODULE: polydef
!> @version 5.0
!> @date 4-1-16
!> @author Christopher Gust
!
!DESCRIPTION: 
!> Defines a derived data types (polydetails and solutiondetails) that contains the polynomial and solution details.
!> Also, contains functions and subroutines involved in initializing polydetails. 
!-----------------------------------------------------------------------------------------------------------------

module polydef
use linear_solution
implicit none  

!------------------------------------------------------------------------------------------------------------------
! Derived Data Type: polydetails
!
!DESCRIPTION: 
!>Polydetails contains most of the Smolyak polynomial details needed to simulate decision rule.
!>  @param nfunc Number of functions approximated.
!>  @param nmsv Number of variables in the minimum state vector (MSV).  
!>  @param nvars Number of endogenous variables.
!>  @param ngrid  Number of grid points used to solve model.
!>  @param nparams  Number of structural model parameters.
!>  @param nindplus Number of variables in which polynomial goes to fourth order.
!>  @param nexog Total number of exogenous variables.
!>  @param nexogshock Number of active exogenous shocks modelled using spline functions.
!>  @param ns Number of exogenous grid points and regime-specific states.
!>  @param ninter Number of points used for interpolating the shocks.
!>  @param indplus Array that indicates the state variable which go to the fourth order.  
!>  @param nshockgrid Vector with number of grid points per shock.!>  @param nstate_reg Vector with number of states per regime-specific variable.
!>  @param interpolmat Matrix indexing grid values of shock used to interpolate between them.
!>  @param endogsteady Model steady state.
!>  @param slopeconmsv Slope coefficients and constants for function going from msv to xx.
!>  @param shockbounds Matrix with upper and lower bounds for each shock.
!>  @param shockdistance Distance between shock values along the grid.
!>  @param slopeconmsv Slope coefficients and constants for function going from msv to xx.
!>  @param exoggrid Matrix indexing grid values of the shocks.
!>  @param ghnodes Matrix of gauss-hermite nodes used to integrate over shocks.
!>  @param ghweights Vector of gauss-hermite weights used to integrate over shocks.
!-----------------------------------------------------------------------------------------------------------------
type :: polydetails
     integer :: nfunc  
     integer :: nmsv   
     integer :: nvars 
     integer :: ngrid   
     integer :: nparams
     integer :: nindplus 
     integer :: nexog
     integer :: nexogshock
     integer :: nexogcont
     integer :: ns
     integer :: nsexog
     integer :: ninter
     integer :: nquad
     logical :: zlbswitch
     integer, allocatable :: indplus(:)  
     integer, allocatable :: nshockgrid(:)
     integer, allocatable :: interpolatemat(:,:)   
     double precision, allocatable :: endogsteady(:)
     double precision, allocatable :: slopeconmsv(:)
     double precision, allocatable :: shockbounds(:,:)
     double precision, allocatable :: shockdistance(:)
     double precision, allocatable :: exoggrid(:,:)    
     double precision, allocatable :: ghnodes(:,:)
     double precision, allocatable :: ghweights(:)
end type polydetails


!------------------------------------------------------------------------------------------------------------------
! Derived Data Type: solutiondetails
!
!DESCRIPTION: 
!> Solutiondetails contains most of the details needed to solve the model.
!>  @param poly Polynomial details.
!>  @param linsol Linear solution details.
!>  @param number_shock_values Number of shock values.
!>  @param startingguess Logical variable indicating whether a starting guess is available for polynomial coefficients.
!>  @param exogvarinfo Matrix indexing the grid position of each shock.
!>  @param bbt Matrix with Smolyak basis functions along its columns evaluated at each grid point (rows).  
!>  @param bbtinv Matrix inverse of bbt. 
!>  @param xgrid Matrix of grid points along its columns. Each grid point is NmsvX1 
!>  @param alphacoeff Matrix of polynomial coefficients.   
!-----------------------------------------------------------------------------------------------------------------
type :: solutiondetails
     type (polydetails) :: poly
     type (linsoldetails)  :: linsol
     integer                       :: number_shock_values
     logical                       :: startingguess
     integer, allocatable          :: exogvarinfo(:,:)
     double precision, allocatable :: bbt(:,:) 
     double precision, allocatable :: bbtinv(:,:)
     double precision, allocatable :: xgrid(:,:) 
     double precision, allocatable :: alphacoeff(:,:)  
end type solutiondetails

public :: initializesolution,msv2xx,sparsegrid,smolyakpoly,exogposition,finite_grid,initialalphas,&
     simulate_linear,soldeallocate
private :: setgridsize,ghquadrature,exoggridindex

contains

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: initializesolution
!> @author Christopher Gust
!> @date 4-1-16
!> @brief Initializes the derived datatype, solution. Solution is mostly an output. The inputs are startingguess,nexog,nregvar,nshockgrid,nstate_reg,
!> nmsv,nvars,nfunc,nindplus,indplus, and zlbswitch.
!> @param[in] params Model parameters.
!> @param solution Solution details.  
!----------------------------------------------------------------------------------------------------------
subroutine initializesolution(solution)

implicit none 
type (solutiondetails), intent(inout) :: solution

integer, parameter :: nquadsingle = 3
integer :: i,j,k,blocksize,ncall
integer :: nexogadj,nquadadj,nmsvadj

!put shocks into polynomial approximation if necessary
nexogadj = solution%poly%nexog-solution%poly%nexogcont
nmsvadj = solution%poly%nmsv+solution%poly%nexogcont

solution%poly%ngrid = 2*(solution%poly%nmsv+solution%poly%nexogcont)+2*solution%poly%nindplus+1

!set nexogshock,ns, and number_shock_values
call setgridsize(nexogadj,solution%poly%nshockgrid,solution%poly%nexogshock,&
     solution%poly%ninter,solution%poly%ns,solution%number_shock_values)
solution%poly%nquad = nquadsingle**(solution%poly%nexogshock+solution%poly%nexogcont)

!allocate solution matrices
allocate(solution%poly%shockbounds(solution%poly%nexogshock,2))
allocate(solution%poly%shockdistance(solution%poly%nexogshock))  
allocate(solution%exogvarinfo(nexogadj,solution%poly%ns))
allocate(solution%poly%exoggrid(nexogadj,solution%poly%ns))
allocate(solution%poly%interpolatemat(solution%poly%nexogshock,2**solution%poly%nexogshock))
allocate(solution%poly%ghnodes(solution%poly%nexogshock+solution%poly%nexogcont,solution%poly%nquad))
allocate(solution%poly%ghweights(solution%poly%nquad))
allocate(solution%poly%endogsteady(solution%poly%nvars+solution%poly%nexog))
allocate(solution%xgrid(nmsvadj,solution%poly%ngrid))
allocate(solution%bbt(solution%poly%ngrid,solution%poly%ngrid))
allocate(solution%bbtinv(solution%poly%ngrid,solution%poly%ngrid))
allocate(solution%poly%slopeconmsv(2*nmsvadj))
allocate(solution%alphacoeff(solution%poly%nfunc*solution%poly%ngrid,2*solution%poly%ns))

solution%exogvarinfo(1:solution%poly%nexogshock,:) = exoggridindex(solution%poly%nshockgrid,solution%poly%nexogshock,&
     solution%poly%ns)

solution%exogvarinfo(solution%poly%nexogshock+1:nexogadj,:) = 1

!get matrix used for interpolating the shocks
solution%poly%interpolatemat = 0
do i = solution%poly%nexogshock,1,-1
   if (i == solution%poly%nexogshock) then
      blocksize = 1
   else
      blocksize = 2*blocksize
   end if
   ncall = 2**(solution%poly%nexogshock-1)/blocksize
   do j = 1,ncall
      do k = 1,2
         solution%poly%interpolatemat(i,2*blocksize*(j-1)+blocksize*(k-1)+1:2*blocksize*(j-1)+blocksize*k) = k-1 
      end do
   end do
end do

!get quadrature nodes and weights
nquadadj = solution%poly%nexogshock+solution%poly%nexogcont
call ghquadrature(nquadsingle,nquadadj,solution%poly%nquad,solution%poly%ghnodes,solution%poly%ghweights)

!construct sparse grid, bb matrix and its inverse
call sparsegrid(nmsvadj,solution%poly%nindplus,solution%poly%ngrid,solution%poly%indplus,&
     solution%xgrid,solution%bbt,solution%bbtinv)

solution%startingguess = .false.
solution%alphacoeff = 0.0d0

!initialize linear solution and kalman matrices
call initializelinearsolution(solution%poly%nparams,solution%poly%nvars,solution%poly%nexog,&
     solution%poly%nexogshock,solution%poly%nexogcont,solution%linsol)

end subroutine initializesolution

!--------------------------------------------------------------------------------------------------------
! FUNCTION: msv2xx
!> @author Christopher Gust
!> @date 10-30-14
!> @brief Returns endogenous state variables in [-1,1] domain or the inverse function depending on slopeconmsv.
!> @param[in] msv Endogenous state variables.
!> @param[in] nmsv Minimum number of endogenous state variables.
!> @param[in] slopeconmsv Slope and constants that depend on the upper and lower bounds of minimum state variables.
!> @return Value of endogenous state variables in [-1,1] domain or the inverse function.  
!----------------------------------------------------------------------------------------------------------
function msv2xx(msv,nmsv,slopeconmsv)

  implicit none 
  integer, intent(in)             :: nmsv
  double precision, intent(in)    :: msv(nmsv)
  double precision, intent(in)    :: slopeconmsv(2*nmsv)
  double precision                :: msv2xx(nmsv)

  msv2xx = slopeconmsv(1:nmsv)*msv + slopeconmsv(nmsv+1:2*nmsv)

end function msv2xx

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: sparsegrid
!> @author Christopher Gust
!> 9-13-14
!> @brief Returns Smolyak grid points and Smolyak matrix bb and its inverse. 
!>    The bb matrix (in JMM 2013) contains the basis functions evaluated at one point along its rows.
!>    We compute its transpose (bbt) since we will use it to grab points in column major form later. 
!>    We also need the inverse of bbt to update the polynomial coefficients.
!> @param[in] nmsv Minimum number of state variables.
!> @param[in] nindplus Number of variables to include up to fourth order in polynomial  
!> @param[in] ngrid Number of grid points.
!> @param[in] indplus Indicator array for variables that go to fourth order.
!> @param[out] xgrid Matrix of grid points
!> @param[out] bbt Matrix of Smolyak basis functions
!> @param[out] bbtinv Inverse matrix of bbt.
!----------------------------------------------------------------------------------------------------------
subroutine sparsegrid(nmsv,nindplus,ngrid,indplus,xgrid,bbt,bbtinv)

implicit none 
integer, intent(in)             :: nmsv,nindplus,ngrid
integer, intent(in)             :: indplus(nindplus)
double precision, intent(out)   :: xgrid(nmsv,ngrid)
double precision, intent(out)   :: bbt(ngrid,ngrid)
double precision, intent(out)   :: bbtinv(ngrid,ngrid)

integer                         :: i,info
integer                         :: ipiv(ngrid)
double precision                :: work(ngrid)

!form smolyak grid
xgrid = 0.0d0
do i = 1,nmsv
   xgrid(i,2*i) = -1.0d0
   xgrid(i,2*i+1) = 1.0d0
end do

do i = 1,nindplus
   xgrid(indplus(i),2*nmsv+2*(i-1)+2) = -1.0d0/sqrt(2.0d0)
   xgrid(indplus(i),2*nmsv+2*(i-1)+3) = 1.0d0/sqrt(2.0d0)
end do

!form bbt matrix 
bbt = 0.0d0
do i = 1,ngrid
   bbt(:,i) = smolyakpoly(nmsv,ngrid,nindplus,indplus,xgrid(:,i))
end do

bbtinv = bbt
call dgetrf(ngrid,ngrid,bbtinv,ngrid,ipiv,info)
if (info .eq. 0) then
   call dgetri(ngrid,bbtinv,ngrid,ipiv,work,ngrid,info)
else
   write(*,*) 'something went wrong with dgetrf (sparsegrid)'
   write(*,*) 'info = ', info
end if

end subroutine sparsegrid

!--------------------------------------------------------------------------------------------------------
! FUNCTION: smolyakpoly
!> @author Christopher Gust
!> @date 9-13-14
!> @brief Returns vector of Smolyak polynomials.
!> @param[in] nmsv Minimum number of state variables.  
!> @param[in] ngrid Number of grid points.
!> @param[in] nindplus Number of variables to include up to fourth order in polynomial
!> @param[in] indplus Indicator array for variables that go to fourth order.
!> @param[in] xx Endogenous state variables defined over [-1,1] domain.
!> @return Vector of Smolyak polynomials.  
!----------------------------------------------------------------------------------------------------------
function smolyakpoly(nmsv,ngrid,nindplus,indplus,xx)

  implicit none 
  integer, intent(in)             :: nmsv,ngrid,nindplus
  integer, intent(in)             :: indplus(nindplus)
  double precision, intent(in)    :: xx(nmsv)
  double precision                :: smolyakpoly(ngrid)

  integer                         :: i

  smolyakpoly(1) = 1.0d0
  do i = 1,nmsv
     smolyakpoly(2*i) = xx(i)
     smolyakpoly(2*i+1) = 2.0d0*xx(i)**2-1.0d0
  end do

  do i = 1,nindplus
     smolyakpoly(2*nmsv+2*(i-1)+2) = 4.0d0*xx(indplus(i))**3-3.0d0*xx(indplus(i))
     smolyakpoly(2*nmsv+2*(i-1)+3) = 8.0d0*xx(indplus(i))**4-8.0d0*xx(indplus(i))**2+1.0d0
  end do

end function smolyakpoly

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: setgridsize
!> @author Christopher Gust
!> @date 4-20-15
!> @brief Sets grid size for exogenous shocks.
!> @param[in] nexog Number of exogenous variables including those shocks held constant.  
!> @param[in] nshockgrid Vector containing grid size for each shock.
!> @param[out] nexogshock Number of active shocks (with states greater than one.)
!> @param[out] ns Total number of grid points.
!> @param[out] number_shock_values Values of shocks at the gridpoints.
!----------------------------------------------------------------------------------------------------------
subroutine setgridsize(nexog,nshockgrid,nexogshock,ninter,ns,number_shock_values)

implicit none 
integer, intent(in)             :: nexog
integer, intent(in)             :: nshockgrid(nexog)
integer, intent(out)            :: nexogshock
integer, intent(out)            :: ninter
integer, intent(out)            :: ns
integer, intent(out)            :: number_shock_values

integer                         :: i

nexogshock = 0
do i = 1,nexog
   if (nshockgrid(i) > 1) then
      nexogshock = nexogshock + 1
   else
      exit
   end if
end do
ninter = 2**nexogshock

!allocate matrices for shock processes
ns = 1
number_shock_values = 0
do i = 1,nexogshock
   ns = ns*nshockgrid(i)
   number_shock_values = number_shock_values + nshockgrid(i)   
end do

end subroutine setgridsize

!--------------------------------------------------------------------------------------------------------
! FUNCTION: exoggridindex
!> @author Christopher Gust
!> @date 2-19-15

!> @brief Returns the matrix of index values for each shock in the grid. The grid has a total of ns points of 
!> dimension nexogX1.
!> @param[in] ngrid Array containing number of shock values in the grid for each of the nexogshockX1 shocks.
!> @param[in] nexogshock Number of exogenous shocks.
!> @param[in] nexog Number of exogenous variables including shocks -- must be greater than nexogshock.  
!> @param[in] ns Total number of grid points.
!> @return Matrix of all possible states for the exogenous variables.
!----------------------------------------------------------------------------------------------------------
function exoggridindex(ngrid,nexog,ns)
             
implicit none 
integer, intent(in) :: nexog,ns,ngrid(nexog)
integer :: exoggridindex(nexog,ns)                 

integer ie,ic,ib,blocksize,ncall

exoggridindex = 0
do ie = nexog,1,-1
   if (ie == nexog) then
      blocksize = 1
   else
      blocksize = blocksize*ngrid(ie+1)
   end if
   ncall = ns/(blocksize*ngrid(ie))
   do ic = 1,ncall
      do ib = 1,ngrid(ie)
         exoggridindex(ie,ngrid(ie)*blocksize*(ic-1)+blocksize*(ib-1)+1:ngrid(ie)*blocksize*(ic-1)+blocksize*ib) = ib
      end do
   end do
end do

end function exoggridindex

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: finite_grid
!> @author Christopher Gust
!> @date 2-17-14

!> @brief Returns an array containing the values of shock at which model is solved.
!>  
!> @param[in] n Number of shock realizations.
!> @param[in] rho AR(1) coefficient.
!> @param[in] sigmaep Standard deviation of innovation.
!> @param[out] zstep Distance between shock values. 
!> @param[out] finite_grid Grid of shock values.
!---------------------------------------------------------------------------------------------------------------
subroutine finite_grid(n,rho,sigmaep,zstep,shockgrid)

implicit none
integer, intent(in) :: n     
double precision, intent(in) :: rho      
double precision, intent(in) :: sigmaep  
double precision, intent(out) :: shockgrid(n) 
double precision, intent(out) :: zstep

double precision, parameter :: maxgridstd = 3.0d0 
integer :: i             
double precision :: nu 
 

shockgrid = 0.0d0
nu = sqrt(1.0d0/(1.0d0-rho**2))*maxgridstd*sigmaep
shockgrid(n)  = nu
shockgrid(1)  = -shockgrid(n)
zstep = (shockgrid(n) - shockgrid(1)) / dble(n - 1)

do i = 2,n-1
    shockgrid(i) = shockgrid(1) + zstep * dble(i - 1)
end do

end subroutine finite_grid

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: qhquadrature
!> @author Christopher Gust
!> @date 8-27-15

!> @brief Produces gauss-hermite quadrature weights and nodes for mulivariate case from univariate case.
!>  
!> @param[in] nquadsingle Number of univariate quadrature points.  (must be 3,5 or 7)
!> @param[in] nexog Number of shocks.
!> @param[out] nquad Total number of quadrature nodes.
!> @param[out] qhnodes Multivariate quadrature nodes
!> @param[out] ghweights Multivariate quadrature weights. 
!---------------------------------------------------------------------------------------------------------------
subroutine ghquadrature(nquadsingle,nexog,nquad,ghnodes,ghweights) 

implicit none
integer, intent(in) :: nquadsingle
integer, intent(in) :: nexog
integer, intent(out) :: nquad  
double precision, intent(out) :: ghnodes(nexog,nquadsingle**nexog)
double precision, intent(out) :: ghweights(nquadsingle**nexog)

integer :: middle_integer
integer, parameter, dimension(3) :: nquadsingleset = (/3,5,7/) 
double precision, parameter :: const_pi = 3.14159265358979323846d0
integer :: ie,blocksize,ncall,ic,ib
double precision :: quadnodes_s(nquadsingle), quadweights_s(nquadsingle)
double precision :: ghweights_mat(nexog,nquadsingle**nexog)

if (any(nquadsingle .eq. nquadsingleset) .eq. .false.) then
   write(*,*) 'WARNING: nquadsingle must equal 3,5, or 7 (ghquadrature in polydef) and you set nquadsingle = ', nquadsingle
end if 

quadnodes_s = 0.0d0
quadweights_s = 0.0d0
if (nquadsingle .eq. 3) then
   quadnodes_s(1) = -sqrt(6.0d0)/2.0d0
   quadnodes_s(2) = 0.0d0
   quadweights_s(1) = sqrt(const_pi)/6.0d0
   quadweights_s(2) = 2.0d0*sqrt(const_pi)/3.0d0
elseif (nquadsingle .eq. 5) then
   quadnodes_s(1) =  -2.0201870456d0
   quadnodes_s(2) =  -0.9585724646d0 
   quadnodes_s(3) =   0.0d0
   quadweights_s(1) = 0.019953242059d0 
   quadweights_s(2) = 0.39361932315d0
   quadweights_s(3) = 0.94530872048d0
else
   quadnodes_s(1) =   -2.6519613568d0
   quadnodes_s(2) =   -1.6735516287d0
   quadnodes_s(3) =   -0.8162878828d0
   quadnodes_s(4) =   0.0d0
   quadweights_s(1) = 9.7117245d-04
   quadweights_s(2) = 0.054515582819d0
   quadweights_s(3) = 0.4256072526d0
   quadweights_s(4) = 0.8102646175568d0
end if

middle_integer = (nquadsingle+1)/2
do ie = 1,middle_integer-1
   quadnodes_s(nquadsingle-ie+1) = -quadnodes_s(ie)
   quadweights_s(nquadsingle-ie+1) = quadweights_s(ie)
end do

nquad = nquadsingle**nexog
ghnodes = 0.0d0
ghweights_mat = 0.0d0
do ie = nexog,1,-1
   if (ie == nexog) then
      blocksize = 1
   else
      blocksize = nquadsingle*blocksize
   end if
   ncall = nquad/(nquadsingle*blocksize)
   do ic = 1,ncall
      do ib = 1,nquadsingle
         ghnodes(ie,nquadsingle*blocksize*(ic-1)+blocksize*(ib-1)+1:nquadsingle*blocksize*(ic-1)+blocksize*ib) = sqrt(2.0d0)*quadnodes_s(ib)
         ghweights_mat(ie,nquadsingle*blocksize*(ic-1)+blocksize*(ib-1)+1:nquadsingle*blocksize*(ic-1)+blocksize*ib) = quadweights_s(ib)
      end do
   end do
end do

ghweights = (1.0d0/const_pi)**(nexog/2.0d0)*product(ghweights_mat,dim=1)

end subroutine ghquadrature

!--------------------------------------------------------------------------------------------------------
! FUNCTION: exogposition
!> @author Christopher Gust
!> @date 4-13-16

!> @brief  Computes the position of the regime and position on exogenous grid.  Hard-wired for 6 or 7 variable for speed.
!> @param[in] exogvec Array with index for each shock value.
!> @param[in] nrvec Number of realizations for each shock.
!> @param[in] nlength Length of exogenous processes.
!> @return Index for the exogenous grid.
!
!---------------------------------------------------------------------------------------------------------------
function exogposition(exogvec,nrvec,nlength)

implicit none
integer, intent(in)           :: nlength
integer, intent(in)           :: nrvec(nlength)
integer, intent(in)           :: exogvec(nlength)
integer                       :: exogposition

if (nlength .eq. 6) then
   exogposition = nrvec(6)*nrvec(5)*nrvec(4)*nrvec(3)*nrvec(2)*(exogvec(1)-1) + nrvec(6)*nrvec(5)*nrvec(4)*nrvec(3)*(exogvec(2)-1) + &
        nrvec(6)*nrvec(5)*nrvec(4)*(exogvec(3)-1) + nrvec(6)*nrvec(5)*(exogvec(4)-1) + nrvec(6)*(exogvec(5)-1) + exogvec(6)
elseif (nlength .eq. 5) then
   exogposition = nrvec(5)*nrvec(4)*nrvec(3)*nrvec(2)*(exogvec(1)-1) + nrvec(5)*nrvec(4)*nrvec(3)*(exogvec(2)-1) + &
        nrvec(5)*nrvec(4)*(exogvec(3)-1) + nrvec(5)*(exogvec(4)-1) + exogvec(5)
elseif (nlength .eq. 4) then
   exogposition = nrvec(4)*nrvec(3)*nrvec(2)*(exogvec(1)-1) + nrvec(4)*nrvec(3)*(exogvec(2)-1) + &
        nrvec(4)*(exogvec(3)-1) + exogvec(4)
elseif (nlength .eq. 3) then
   exogposition = nrvec(3)*nrvec(2)*(exogvec(1)-1) + nrvec(3)*(exogvec(2)-1) + exogvec(3)
elseif (nlength .eq. 2) then
   exogposition = nrvec(2)*(exogvec(1)-1) + exogvec(2)
elseif (nlength .eq. 1) then
   exogposition = exogvec(1)
else
   stop 'There can only be six shocks. You need to modify exogposition. (exogposition - module polydef).'
end if

end function exogposition

!--------------------------------------------------------------------------------------------
! SUBROUTINE: initialalphas
!> @author Christopher Gust
!> @date 4-16-16
!> @brief Uses linear decision rule to construct initial guess of polynomial coefficients.  
!>  @param[in] nfunc Number of functions approximated
!>  @param[in] ngrid  Number of grid points used to solve mode
!>  @param[in] ns Number of exogenous states.
!>  @param[in] nvars Number of endogenous variables used in nonlinear solution.
!>  @param[in] nexog Number of exogenous shocks.
!>  @param[in] nmsv Number of variables in the minimum state vector (MSV)  
!>  @param[in] shockinfo Matrix indexing each shock value for each state.
!>  @param[in] slopeconxx Slope coefficients and constants used to go from xx [-1,1] space to msv space.
!>  @param[in] xgrid Matrix of grid points along its columns. Each grid point is NmsvX1. 
!>  @param[in] bbtinv Matrix inverse of bbt. 
!>  @param[in] aalin Feedback part of linear decision rule that excludes shock processes.  
!>  @param[in] bblin Linear decision rule matrix that maps effect of shocks (not innovations) on endogenous variables.
!>  @param[in] endogsteady Steady state values of model's variables.
!>  @return Polynomial coefficients. 
!---------------------------------------------------------------------------------------------
function initialalphas(nfunc,ngrid,ns,nvars,nexog,nexogshock,nexogcont,nmsv,exoggrid,slopeconxx,&
     xgrid,bbtinv,aalin,bblin,endogsteady)

implicit none
integer, intent(in)            :: nfunc,ngrid,ns,nvars,nexog,nmsv,nexogshock,nexogcont
double precision, intent(in)   :: exoggrid(nexog-nexogcont,ns)
double precision, intent(in)   :: slopeconxx(2*(nmsv+nexogcont))
double precision, intent(in)   :: xgrid(nmsv+nexogcont,ngrid)
double precision, intent(in)   :: bbtinv(ngrid,ngrid)
double precision, intent(in)   :: aalin(nvars,nvars)
double precision, intent(in)   :: bblin(nvars,nexog)
double precision, intent(in)   :: endogsteady(nvars+nexog)
double precision               :: initialalphas(nfunc*ngrid,2*ns)

integer            :: ss,i,ifunc
integer            :: nmsvplus
double precision   :: yy(nfunc,ngrid)
double precision   :: endogvarm1(nvars,ngrid)
double precision   :: endogvar(nvars),exogpart(nvars)
double precision   :: exogval(nexog)
double precision   :: alphass(nfunc,ngrid)
double precision   :: slopeconxxmsv(2*nmsv)
double precision   :: slopeconcont(2*nexogcont)

initialalphas = 0.0d0
alphass = 0.0d0

nmsvplus = nmsv + nexogcont
slopeconxxmsv(1:nmsv) = slopeconxx(1:nmsv)
slopeconxxmsv(nmsv+1:2*nmsv) = slopeconxx(nmsvplus+1:nmsvplus+nmsv)
if (nexogcont .gt. 0) then
   slopeconcont(1:nexogcont) = slopeconxx(nmsv+1:nmsvplus)
   slopeconcont(nexogcont+1:2*nexogcont) = slopeconxx(nmsvplus+nmsv+1:2*nmsvplus)   
end if

endogvarm1 = 0.0d0
do i = 1,ngrid
   endogvarm1(1:nmsv,i) = msv2xx(xgrid(1:nmsv,i),nmsv,slopeconxxmsv)-endogsteady(1:nmsv)
end do

do ss = 1,ns
   yy = 0.0d0
   exogval = 0.0d0
   do i = 1,ngrid
      exogval(1:nexogshock) = exoggrid(1:nexogshock,ss)   !update shocks (in deviation from ss)
      if (nexogcont .gt. 0) then
         exogval(nexog-nexogcont+1:nexog) = msv2xx(xgrid(nmsv+1:nmsv+nexogcont,i),nexogcont,slopeconcont)-endogsteady(nvars+nexog-nexogcont+1:nvars+nexog)
      end if

      !get linear solution
      call dgemv('N', nvars, nvars, 1.0d0, aalin, nvars, endogvarm1(:,i), 1, &
          0.0d0, endogvar, 1)
      call dgemv('N', nvars, nexog, 1.0d0, bblin, nvars, exogval, 1, &
          0.0d0, exogpart, 1)
      endogvar = endogsteady + endogvar + exogpart
      yy(:,i) = endogvar( (/10,11,18,19,21,22,13/) )   
   end do
  
   call dgemm('N','N',nfunc,ngrid,ngrid,1.0d0,yy,nfunc,bbtinv,ngrid,0.0d0,alphass,nfunc)

   do i = 1,ngrid
      do ifunc = 1,nfunc
         !if (alphass(ifunc,i) < 1.0e-8) alphass(ifunc,i) = 0.0d0
         initialalphas((ifunc-1)*ngrid+i,ss) = alphass(ifunc,i)
         !initial guess for ZLB polynomials
         initialalphas((ifunc-1)*ngrid+i,ns+ss) = alphass(ifunc,i)
         !print*,initialalphas((ifunc-1)*ngrid+i,ns+ss)
      end do
   end do
end do

end function initialalphas

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: simulate_linear
!> @author Christopher Gust
!> @date 3-22-15

!> @brief Simulate data to compute ergodic means of economy's variables.
!> @param[in] pp Part of the linear decision rule (feedback part).  
!> @param[in] sigma Part of the linear decision rule (innovation part).
!> @param[in] endogsteady Steady state values of model variables
!> @param[in] nexog Number of exogenous variables (shocks and fixed values).
!> @param[in] nexogshock Number of exogenous shocks.
!> @param[in] nvars Number of endogenous variables.
!> @param[in] nmsv Number of minimum state endogenous variables.
!> @param[out] endog_emean Ergodic means of variables.
!> @param[out] zlbfrequency Frequency at ZLB.  
!> @param[out] msvbounds Bounds on minimum state variables (first nmsv are lower; next are upper). 
!> @param[out] statezlbinfo Indicator for whether 1 or 2 polynomials used to approximate solution at particular state.
!---------------------------------------------------------------------------------------------------------
subroutine simulate_linear(linsol,ns,nmsv,shockbounds,shockdistance,nshockgrid,&
     endog_emean,zlbfrequency,msvbounds,statezlbinfo,convergence)

use rng_serial
implicit none
type(linsoldetails), intent(in) :: linsol
integer, intent(in) :: nmsv,ns
double precision, intent(in)    :: shockbounds(linsol%nexogshock,2)
double precision, intent(in)    :: shockdistance(linsol%nexogshock)
integer,          intent(in)    :: nshockgrid(linsol%nexog-linsol%nexogcont)
integer,          intent(out)   :: statezlbinfo(ns)
double precision, intent(out)   :: endog_emean(linsol%nvars)
double precision, intent(out)   :: zlbfrequency
double precision, intent(out)   :: msvbounds(2*(nmsv+linsol%nexogcont))
logical, intent(out)   :: convergence

integer, parameter  :: capt = 100000
integer            :: counter,displacement,llim,ulim,ttsim,countzlb,i
integer            :: count_exploded
integer            :: nmsvplus
integer            :: shockindex(linsol%nexog-linsol%nexogcont),stateindex
integer            :: countzlbstates(ns)
double precision   :: endogvar(linsol%nvars,0:capt)
double precision   :: msvhigh(nmsv),msvlow(nmsv)
double precision   :: innovations(linsol%nexog)
logical            :: explosiveerror,non_explosive
double precision   :: xrandn(linsol%nexog,capt)
double precision   :: scalebd,msv_std(nmsv+linsol%nexogcont)
integer :: brng
integer :: method
integer :: errcode
integer :: iseed
type (vsl_stream_state) :: stream     

!get random normals
iseed = 101294
brng = vsl_brng_mt19937
method = vsl_rng_method_gaussian_boxmuller 
errcode = vslnewstream( stream,   brng,  iseed )
errcode = vdrnggaussian( method, stream, linsol%nexog*capt, xrandn, 0.0d0, 1.0d0)

nmsvplus = nmsv + linsol%nexogcont
convergence = .true.
endogvar = 0.0d0
counter = 0
countzlb = 0
countzlbstates = 0
statezlbinfo = 0
displacement = 400
endogvar(1:nmsv,0) = linsol%endogsteady(1:nmsv)
msvhigh = log(exp(linsol%endogsteady(1:nmsv))*2.0d0)
msvlow = log(exp(linsol%endogsteady(1:nmsv))*0.01d0)
innovations = 0.0d0
count_exploded = 0

counterloop: do 
   explosiveerror = .false.
   llim = counter+1
   ulim = min(capt,counter+displacement)
   do ttsim = llim,ulim
      innovations(1:linsol%nexogshock) = xrandn(1:linsol%nexogshock,ttsim)
      if (linsol%nexogcont > 0) innovations(linsol%nexog-linsol%nexogcont+1:linsol%nexog) = xrandn(linsol%nexog-linsol%nexogcont+1:linsol%nexog,ttsim)
      call decrlin(endogvar(:,ttsim-1),innovations,linsol,endogvar(:,ttsim))
      if (endogvar(5,ttsim) .lt. 0.0d0) then
         countzlb = countzlb + 1
         shockindex = 1
         do i = 1,linsol%nexogshock
            if (endogvar(linsol%nvars-linsol%nexog+i,ttsim) < shockbounds(i,1)) then 
               shockindex(i) = 1
            else  !interpolation case
               shockindex(i) = min( nshockgrid(i), int(1.0d0+(endogvar(linsol%nvars-linsol%nexog+i,ttsim)-&
                    shockbounds(i,1))/shockdistance(i)) )  
            end if
         end do
         stateindex = exogposition(shockindex,nshockgrid,linsol%nexog-linsol%nexogcont)
         countzlbstates(stateindex) = countzlbstates(stateindex) + 1
         if (countzlbstates(stateindex) > 5) then
            statezlbinfo(stateindex) = 1 
         end if
      end if
   end do

   checkloop: do ttsim = llim,ulim
      non_explosive = ( all(endogvar(1:nmsv,ttsim) .le. msvhigh) .and.  &
           all(endogvar(1:nmsv,ttsim) .ge. msvlow) ) 
      explosiveerror = ( (non_explosive .eq. .false.) .or. (isnan(endogvar(1,ttsim)) .eq. .true.) )
      if (explosiveerror .eq. .true.)  then
         counter = max(ttsim-200,0)
         errcode = vdrnggaussian( method, stream, linsol%nexog*capt, xrandn, 0.0d0, 1.0d0)
         if (explosiveerror == .true.) then 
            !write(*,*) 'solution exploded at ', ttsim, ' vs ulim ', ulim
	    count_exploded = count_exploded + 1
         end if
          
         if (counter .eq. 0) then
	    convergence = .false.
            write(*,*) 'degenerated back to the beginning'
         end if
         exit checkloop
      end if
   end do checkloop

   if (count_exploded .eq. 50) then
      write(*,*) 'linear solution exploded for the 50th time'
      convergence = .false.
   end if 

   if (convergence .eq. .false.) then
      return
   end if

   if (explosiveerror .eq. .false.) then
      counter = counter + displacement
   end if

   if (counter .ge. capt) then
      exit counterloop
   end if

end do counterloop

!calculate mean of all variables; std. dev of minimum state variables and exogenous variables
!that we include in polynomial part of approximated decision rule
zlbfrequency = 100*dble(countzlb)/dble(capt)
scalebd = 3.0d0
endog_emean = sum(endogvar(:,1:capt),dim=2)/dble(capt)

do counter = 1,nmsv
   msv_std(counter) = sqrt( sum( (endogvar(counter,1:capt)-endog_emean(counter))**2 )/dble(capt-1) )
end do
msvbounds(1:nmsv) = endog_emean(1:nmsv)-scalebd*msv_std
msvbounds(nmsvplus+1:nmsvplus+nmsv) = endog_emean(1:nmsv)+scalebd*msv_std
do counter = 1,linsol%nexogcont
   msv_std(nmsv+counter) = sqrt( sum( (endogvar(linsol%nvars-counter+1,1:capt)-endog_emean(linsol%nvars-counter+1))**2 )/dble(capt-1) )
   msvbounds(nmsv+counter) = endog_emean(linsol%nvars-counter+1)-scalebd*msv_std(nmsv+counter)
   msvbounds(nmsvplus+nmsv+counter) = endog_emean(linsol%nvars-counter+1)+scalebd*msv_std(nmsv+counter)
end do

end subroutine simulate_linear

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: soldeallocate
!> @author Christopher Gust
!> @date 4-13-16

!> @brief Deallocate the solution matrices.
!>  
!> @param solution Solution details (in\out).
!---------------------------------------------------------------------------------------------------------------
subroutine soldeallocate(solution) 

implicit none
type (solutiondetails), intent(inout) :: solution

deallocate(solution%poly%nshockgrid) 
deallocate(solution%poly%shockbounds)
deallocate(solution%poly%shockdistance)  
deallocate(solution%exogvarinfo)
deallocate(solution%poly%exoggrid)
deallocate(solution%poly%interpolatemat)
deallocate(solution%poly%ghnodes)
deallocate(solution%poly%ghweights)
deallocate(solution%poly%indplus)
deallocate(solution%xgrid)
deallocate(solution%bbt)
deallocate(solution%bbtinv)
deallocate(solution%poly%slopeconmsv)
deallocate(solution%poly%endogsteady)
deallocate(solution%alphacoeff)
call linsoldeallocate(solution%linsol)

end subroutine soldeallocate

end module polydef





