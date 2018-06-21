!------------------------------------------------------------------------------------------------------------------
! MODULE: get_decisionrule
!> @version 5.0
!> @author Christopher Gust
!
!DESCRIPTION
!> Functions and subroutines used to define and find model's decision rule.
!-----------------------------------------------------------------------------------------------------------------
module get_decisionrule

use model_details, only: get_shockdetails, decr_euler, steadystate, simulate_linear, initialalphas, solutiondetails, polydetails 
use linear_solution, only: get_aimsolution,lindecrule_markov
  
implicit none
public :: nonlinearsolver
private :: fixedpoint
  
contains

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: nonlinearsolver
!> @author Christopher Gust
!> @date 6-28-16

!> @brief Takes parameters and initial guess to solve nonlinear model using fixed point iterations.  
!> @param[in] params Model parameters.  
!> @param solution Solution details is inout.  Most data is input.  Exceptions: alphacoeff (inout), 
!>    slopeconmsv (out), exoggrid(out), shockbounds(out), shockdistance(out), linsolution (out), endogsteady (out).
!> @param[out] convergence Indicator variable reporting back whether polynomial coefficients converged.
!---------------------------------------------------------------------------------------------------------
subroutine nonlinearsolver(params,solution,convergence)

implicit none
type(solutiondetails), intent(inout) :: solution
double precision, intent(in)    :: params(solution%poly%nparams)
logical, intent(out)            :: convergence

integer          :: nmsvplus,i
integer          :: statezlbinfo(solution%poly%ns)
double precision :: avgerror
double precision :: aalin(solution%poly%nvars,solution%poly%nvars)
double precision :: bblin(solution%poly%nvars,solution%poly%nexog)
double precision :: alphacoeff0(solution%poly%nfunc*solution%poly%ngrid,2*solution%poly%ns)
double precision :: alphacoeffstar(solution%poly%nfunc*solution%poly%ngrid,2*solution%poly%ns)
double precision :: msvbounds(2*(solution%poly%nmsv+solution%poly%nexogcont))
double precision :: slopeconxx(2*(solution%poly%nmsv+solution%poly%nexogcont))
double precision :: zlbfrequency,endog_emean(solution%poly%nvars+solution%poly%nexog)

!get shock details
call get_shockdetails(solution%poly%nparams,solution%poly%nexog,solution%poly%nexogshock,solution%poly%nexogcont,&
     solution%poly%ns,solution%number_shock_values,solution%poly%nshockgrid,params,solution%exogvarinfo,&
     solution%poly%exoggrid,solution%poly%shockbounds,solution%poly%shockdistance)

!get bounds after computing variances from linear solution; 
!also decide which exogenous grid points need 2 polynomials instead of 1 based on whether ZLB is reached using linear solution
solution%poly%endogsteady = steadystate(params,solution%linsol%nvars,solution%poly%nparams)
solution%linsol%endogsteady = solution%poly%endogsteady
call get_aimsolution(params,solution%linsol)

call simulate_linear(solution%linsol,solution%poly%ns,solution%poly%nmsv,solution%poly%shockbounds,solution%poly%shockdistance,&
     solution%poly%nshockgrid,endog_emean,zlbfrequency,msvbounds,statezlbinfo,convergence)

!check to see how many times the linear model wants a  "nearly-explosive" path
if (convergence .eqv. .false.) then
   return
end if 

nmsvplus = solution%poly%nmsv+solution%poly%nexogcont
solution%poly%slopeconmsv(1:nmsvplus) = 2.0d0/(msvbounds(nmsvplus+1:nmsvplus)-msvbounds(1:nmsvplus))
solution%poly%slopeconmsv(nmsvplus+1:2*nmsvplus) = -2.0d0*msvbounds(1:nmsvplus)/&
     (msvbounds(nmsvplus+1:2*nmsvplus)-msvbounds(1:nmsvplus))-1.0d0
slopeconxx(1:nmsvplus) = 0.5d0*(msvbounds(nmsvplus+1:2*nmsvplus)-&
     msvbounds(1:nmsvplus))
slopeconxx(nmsvplus+1:2*nmsvplus) = msvbounds(1:nmsvplus) + 0.5d0*&
  (msvbounds(nmsvplus+1:2*nmsvplus)-msvbounds(1:nmsvplus))

!if no starting guess, construct one from linear solution
if (solution%startingguess .eqv. .false.) then
   call lindecrule_markov(solution%linsol,aalin,bblin)
   
   alphacoeff0 = initialalphas(solution%poly%nfunc,solution%poly%ngrid,&
        solution%poly%ns,solution%poly%nvars,solution%poly%nexog,solution%poly%nexogshock,&
        solution%poly%nexogcont,solution%poly%nmsv,solution%poly%exoggrid,slopeconxx,solution%xgrid,solution%bbtinv,&
        aalin,bblin,solution%poly%endogsteady)
else
   alphacoeff0 = solution%alphacoeff
end if

!find metaparameters that solve model 
call fixedpoint(params,solution%poly,alphacoeff0,slopeconxx,solution%bbt,&
     solution%bbtinv,solution%xgrid,statezlbinfo,alphacoeffstar,convergence,avgerror)

solution%alphacoeff = alphacoeffstar

! do i = 1,5
!    write(*,*) 'alphacoeff0(1:19,i) = ', alphacoeff0(1:19,i)
!    write(*,*) '-------------------------------'
!    write(*,*) 'alphacoeff(1:19,i) = ', solution%alphacoeff(1:19,i)
!    write(*,*) '-------------------------------'
!    write(*,*) '-------------------------------'
! end do

end subroutine nonlinearsolver 

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: fixedpoint
!> @author Christopher Gust
!> @date 11-3-14

!> @brief Iterates on a fixed point for polynomial coefficients.
!> @todo Fix errsum_grid and eliminate matmul.
!> @param[in] params Model parameters.  
!> @param[in]  poly Polynomial details.
!> @param[in] alphacoeff0 Initial guess at polynomial coefficients.
!> @param[in] slopeconxx Slope coefficients and constants for converting xx to msv space. 
!> @param[in] bbt Matrix with Smolyak basis functions along its columns evaluated at each grid point (rows).  
!> @param[in] bbtinv Matrix inverse of bbt. 
!> @param[in] xgrid Matrix of grid points along its columns. Each grid point is NmsvX1 
!> @param[in] statezlbinfo Indicator variable for states with two polynomials used to handle ZLB constraint.
!> @param[out] alphacoeffstar Values of polynomial coefficients at their solution.
!> @param[out] convergence Indicator variable reporting back whether polynomial coefficients converged.
!> @param[out] avgerror Average absolute error of polynomial values.
!> @param Niter Number of iterations (set to a 'large' number.) 
!> @param tolfun Tolerance criteria that determines whether algorithm converged 
!> @param step Step size for updating polynomial coefficients (set between 0 (small step) and 1 (large step))
!---------------------------------------------------------------------------------------------------------
subroutine fixedpoint(params,poly,alphacoeff0,slopeconxx,bbt,bbtinv,&
     xgrid,statezlbinfo,alphacoeffstar,convergence,avgerror)

implicit none
type(polydetails), intent(in)   :: poly
double precision, intent(in)    :: bbt(poly%ngrid,poly%ngrid)
double precision, intent(in)    :: bbtinv(poly%ngrid,poly%ngrid)
double precision, intent(in)    :: xgrid(poly%nmsv,poly%ngrid)
double precision, intent(in)    :: params(poly%nparams)
double precision, intent(in)    :: alphacoeff0(poly%nfunc*poly%ngrid,2*poly%ns)
double precision, intent(in)    :: slopeconxx(2*poly%nmsv)
integer,          intent(in)    :: statezlbinfo(poly%ns)
double precision, intent(out)   :: alphacoeffstar(poly%nfunc*poly%ngrid,2*poly%ns)
logical, intent(out)            :: convergence
double precision, intent(out)   :: avgerror

integer, parameter :: niter = 150
double precision, parameter :: tolfun = 1.0d-04
double precision, parameter :: step = 7.0d-01
integer :: ii,ss,igrid,ifunc,nrow
double precision :: alphanew(poly%nfunc*poly%ngrid,2*poly%ns)
double precision :: alphacur(poly%nfunc*poly%ngrid,2*poly%ns)
double precision :: polyappnew(2*poly%nfunc,poly%ngrid)
double precision :: errsum_grid,errsum,errmax
double precision :: alphass(2*poly%nfunc,poly%ngrid)

alphacur = alphacoeff0
convergence = .false.
mainloop: do ii = 1,niter

   avgerror = 0.0d0
   do ss = 1,poly%ns
      polyappnew = 0.0d0
      errsum_grid = 0.0d0
      do igrid = 1,poly%ngrid
         call decr_euler(igrid,ss,params,poly,alphacur,slopeconxx,&
              bbt,xgrid,statezlbinfo(ss),polyappnew(:,igrid),errsum,errmax) 
         errsum_grid = errsum_grid + errsum 
      end do

      nrow = 2*poly%nfunc
      call dgemm('N','N',nrow,poly%ngrid,poly%ngrid,1.0d0,polyappnew,nrow,bbtinv,poly%ngrid,0.0d0,alphass,nrow)
      do igrid = 1,poly%ngrid
         do ifunc = 1,poly%nfunc
            alphanew((ifunc-1)*poly%ngrid+igrid,ss) = alphass(ifunc,igrid)
            alphanew((ifunc-1)*poly%ngrid+igrid,poly%ns+ss) = alphass(poly%nfunc+ifunc,igrid)
         end do
      end do
      avgerror = avgerror + errsum_grid/poly%ngrid
   end do
   
   avgerror = avgerror/dble(2*poly%ns) 

   if (any(isnan(alphanew))) then
      convergence = .false.
      alphacoeffstar = alphacur
      exit mainloop
   end if
   
   ! if ( mod(ii,1) == 0) then
   !   write(*,*) 'ii = ', ii
   !   write(*,*) 'avgerror = ', avgerror  
   ! end if

   if (avgerror < tolfun) then 
      convergence = .true.
      alphacoeffstar = alphacur
      !write(*,*) 'avgerror = ', avgerror
      !write(*,*) 'iteration = ', ii
      exit mainloop
   end if
   alphacur = (1.0d0-step)*alphacur + step*alphanew
 
end do mainloop

end subroutine fixedpoint

end module get_decisionrule
