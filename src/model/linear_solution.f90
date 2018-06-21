!------------------------------------------------------------------------------------------------------------------
! MODULE: linear_solution
!
!> @author Christopher Gust
!
!DESCRIPTION: 
!> Procedures for computing the linear solution to the model. 
!-----------------------------------------------------------------------------------------------------------------
module linear_solution

use fortress_random_t, only: fortress_random
use fortress_util, only: write_array_to_file
use gensys, only: do_gensys

implicit none

!------------------------------------------------------------------------------------------------------------------
! Derived Data Type: linsoldetails
!
!DESCRIPTION: 
!>Polydetails contains most of the Smolyak polynomial details needed to simulate decision rule.
!>  @param nparams Number of functions approximated.
!>  @param nvars Number of variables in the minimum state vector (MSV).  
!>  @param nexog Total number of potential shocks.
!>  @param nexogshock Shocks on finite-element part of approximated nonlinear decision rule.
!>  @param nexogcont Shocks on smooth part of approximated nonlinear decision rule.
!>  @param pp  Part of linear decision rule that feeds back on to lagged variables. 
!>  @param sigma Part of linear decision rule that feed back on to innovations.
!>  @param endogsteady Model steady state.
!-----------------------------------------------------------------------------------------------------------------
type :: linsoldetails
     integer :: nparams  
     integer :: nvars
     integer :: nexog 
     integer :: nexogshock
     integer :: nexogcont
     double precision, allocatable :: pp(:,:) 
     double precision, allocatable :: sigma(:,:)
     double precision, allocatable :: endogsteady(:)
end type linsoldetails

public :: initializelinearsolution,get_aimsolution,get_kalmanmatrices,lindecrule_markov,&
     decrlin,decrlin2,simulate_kalman,linsoldeallocate

contains

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: initializelinearsolution
!> @author Christopher Gust
!> @date 10-8-15
!> @brief Initializes the derived datatype for the linear solution.
!> @param[in] nvars Number of endogenous variables in nonlinear decision rule.  linsol%nvars = poly%nvars+nexog. 
!> @param[out] linsol Linear solution details.  
!----------------------------------------------------------------------------------------------------------
subroutine initializelinearsolution(nparams,nvars,nexog,nexogshock,nexogcont,linsol)  

implicit none 
integer, intent(in)             :: nparams
integer, intent(in)             :: nvars
integer, intent(in)             :: nexog
integer, intent(in)             :: nexogshock
integer, intent(in)             :: nexogcont
type (linsoldetails), intent(out) :: linsol

linsol%nparams = nparams
linsol%nvars = nvars+nexog  !nvars is equal to poly%nvars which excludes exogenous shocks
linsol%nexog = nexog
linsol%nexogshock = nexogshock
linsol%nexogcont = nexogcont
allocate(linsol%pp(linsol%nvars,linsol%nvars)) 
allocate(linsol%sigma(linsol%nvars,linsol%nexog)) 
allocate(linsol%endogsteady(linsol%nvars))

end subroutine initializelinearsolution

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: get_aimsolution
!> @author William Gamber and Christopher Gust
!> @date 6-28-16
!> @brief Uses sparse AIM code to solve the linear model in file modelv81.
!> Generates the matrices, $P$ and $\Sigma$ defined by the system:
!> $x_t = P x_{t-1} + \Sigma \epsilon_{t}$.
!> The endogenous variables, $x_t$ are defined in linsol%endognames and 
!> the innovations, $\epsilon_t$, are defined in linsol%innovnames.
!> Part of this code is generated automatically and copied into this module.
!> Note that with -check all option aPointerToVoid generates a runtime error.
!> @param[in] params Model's parameters
!> @param[in] linsol%nparams Number of model parameters
!> @param[in] linsol%nvars Number of endogenous variables in linear model.
!> @param[in] linsol%nexog Number of exogenous (fixed) variables in linear model. 
!> @param[in] linsol%nexogshock Number of exogenous shocks in linear model. 
!> @param[out] linsol%pp Endogenous feedback matrix of linear-solution.
!> @param[out] linsol%sigma Exogenous linear-solution matrix (effect of innovations on endogenous variables).
!----------------------------------------------------------------------------------------------------------
subroutine get_aimsolution(params,linsol)

implicit none
type(linsoldetails), intent(inout) :: linsol
double precision, intent(in) :: params(linsol%nparams)

logical :: convergence

double precision :: beta_tr,pibar_tr,gz_tr,gamxhp,gamma,sigmal,phi,phiw,ap_tr,aw_tr,alpha,phii,sigmaa,gam_rs,gam_dp,gamdy,sdevtech_tr,rhog,sdevg_tr,rhoinv,sdevinv_tr,sdevliq_tr,sdevint_tr,sdeva_tr
double precision :: rbar,gz,rholiq,ep_tr,epw_tr,shrgy,delta,lamhp,psil,beta,pibar,ep,epw,ap,aw,bw,sdevtech,sdevg,sdevinv,sdevliq,sdevint,sdeva,gg,gamtil,mc,k2yrat,shriy,shrcy,labss,kappaw,kappap,kss,gdpss,invss,phii_jpt,css,rwss,mucss,lamss,rss,rkss,rhoint,rhoa

double precision :: rhoelast,sdevelast_tr,rhoelastw,sdevelastw_tr,sdevelastw,sdevelast

! gensys
double precision :: fmat, fwt, ywt, gev, loose, DIV
integer :: eu(2)

double precision :: TT(42, 42), RR(42, 8), QQ(8, 8), DD(5), ZZ(5, 42), HH(8,8)
double precision :: GAM0(42, 42), GAM1(42, 42), C(42), PSI(42, 8), PPI(42, 8), CC(42)     

! !get aim parameters
beta = params(1) 
pibar = params(2) 
gz = params(3)
psil = params(4)
gamma = params(5) 
sigmal = params(6) 
phi = params(7)
phiw = params(8)
ep = params(9) 
epw = params(10) 
ap = params(11)
aw = params(12)
bw = params(13) 
lamhp = params(14) 
alpha = params(15)
delta = params(16)
phii = params(17) 
sigmaa = params(18) 
gam_rs = params(19)
gam_dp = params(20)
gamxhp = params(21) 
gamdy = params(22) 
shrgy = params(23)
sdevtech = params(24)
rhog = params(25) 
sdevg = params(26) 
rhoinv = params(27)
sdevinv = params(28)
rholiq = params(29) 
sdevliq = params(30) 
rhoint = params(31)
sdevint = params(32)
rhoa = params(33) 
sdeva = params(34) 

! !fix bw = aw

gg=1.0/(1.0-shrgy)
gamtil=gamma/gz
mc=(ep-1.0)/ep
k2yrat=((mc*alpha)/(gz/beta-(1.0-delta)))*gz
shriy=(1.0-(1.0-delta)/gz)*k2yrat
shrcy=1.0-shrgy-shriy
labss=( ((epw-1.0)/epw)*(1.0-alpha)*(1.0-beta*gamtil)*((ep-1.0)/ep)*(1.0/(psil*(1.0-gamtil)))*(1.0/shrcy) )**(1.0/(sigmal+1.0))
kappaw=((1.0-gamtil)/(1.0-beta*gamtil))*epw*psil*labss**(1.0+sigmal)/phiw
kappap=(ep-1.0)/(phi*(1.0+beta*(1.0-ap)))
kss=labss*(gz**(alpha/(alpha-1.0)))*k2yrat**(1.0/(1.0-alpha))
gdpss=(kss/gz)**alpha*labss**(1.0-alpha)
invss=shriy*gdpss
phii_jpt=phii/invss
css=shrcy*gdpss
rwss=(1.0-alpha)*mc*gdpss/labss
mucss=(1.0/css)*(1.0/(1.0-gamtil))
lamss=mucss*(1.0-beta*gamtil)
rss=gz*pibar/beta
rkss=gz/beta-1.0+delta

GAM0 = 0.0d0
GAM1 = 0.0d0
PSI = 0.0d0
PPI = 0.0d0
C = 0.0d0

GAM0(3, 1) = 1.0d0
GAM0(2, 2) = -1.0d0*beta*gamtil**2.0d0 - 1.0d0
GAM0(7, 2) = gg*shrcy
GAM0(16, 2) = -1.0d0
GAM0(21, 2) = gamtil*1.0/(-1.0d0*gamtil + 1.0d0)
GAM0(37, 2) = -1.0d0
GAM0(3, 3) = 1.0/gz*(-1.0d0*delta + 1.0d0) - 1.0d0
GAM0(7, 3) = gg*shriy
GAM0(11, 3) = phii_jpt*(beta + 1.0d0)
GAM0(17, 3) = 1.0d0
GAM0(22, 3) = -1.0d0
GAM0(35, 3) = -1.0d0
GAM0(4, 4) = -1.0d0
GAM0(14, 4) = 1.0d0
GAM0(15, 4) = 1.0d0
GAM0(19, 4) = -1.0d0*kappaw
GAM0(5, 5) = -1.0d0
GAM0(9, 5) = 1.0d0
GAM0(10, 5) = 1.0d0
GAM0(4, 6) = -1.0d0
GAM0(5, 6) = gam_dp*(-1.0d0*gam_rs + 1.0d0)
GAM0(6, 6) = -1.0d0
GAM0(18, 6) = 1.0d0
GAM0(40, 6) = -1.0d0
GAM0(5, 7) = gamdy*(-1.0d0*gam_rs + 1.0d0)
GAM0(7, 7) = -1.0d0
GAM0(12, 7) = 1.0d0
GAM0(14, 7) = -1.0d0
GAM0(5, 8) = gamxhp*(-1.0d0*gam_rs + 1.0d0)
GAM0(8, 8) = -1.0d0
GAM0(9, 9) = -1.0d0
GAM0(1, 10) = -1.0d0
GAM0(2, 10) = -1.0d0*(-1.0d0*gamtil + 1.0d0)*(-1.0d0*beta*gamtil + 1.0d0 &
   )
GAM0(10, 10) = -1.0d0
GAM0(19, 10) = -1.0d0*kappaw
GAM0(42, 10) = -1.0d0
GAM0(1, 11) = -1.0d0
GAM0(11, 11) = -1.0d0
GAM0(41, 11) = -1.0d0
GAM0(8, 12) = -1.0d0*alpha + 1.0d0
GAM0(12, 12) = alpha - 1.0d0
GAM0(14, 12) = 1.0d0
GAM0(15, 12) = 1.0d0
GAM0(19, 12) = kappaw*sigmal
GAM0(7, 13) = alpha*1.0/ep*gg*(ep - 1.0d0)
GAM0(8, 13) = alpha
GAM0(12, 13) = -1.0d0*alpha
GAM0(13, 13) = -1.0d0
GAM0(15, 13) = -1.0d0
GAM0(6, 14) = kappap
GAM0(14, 14) = -1.0d0
GAM0(13, 15) = 1.0/sigmaa
GAM0(15, 15) = -1.0d0
GAM0(36, 15) = -1.0d0
GAM0(16, 16) = gamtil - 1.0d0
GAM0(17, 17) = -1.0d0
GAM0(18, 18) = -1.0d0
GAM0(19, 19) = -1.0d0
GAM0(20, 19) = 1.0d0
GAM0(38, 19) = -1.0d0
GAM0(4, 20) = 1.0d0
GAM0(20, 20) = -1.0d0
GAM0(21, 21) = -1.0d0
GAM0(22, 22) = -1.0d0
GAM0(10, 23) = 1.0d0
GAM0(27, 23) = -1.0d0
GAM0(3, 24) = 1.0/gz*(-1.0d0*delta + 1.0d0) - 1.0d0
GAM0(11, 24) = -1.0d0
GAM0(28, 24) = -1.0d0
GAM0(2, 25) = -1.0d0*gamtil
GAM0(3, 25) = 1.0/gz*(-1.0d0*delta + 1.0d0)
GAM0(4, 25) = -1.0d0
GAM0(5, 25) = gamdy*(-1.0d0*gam_rs + 1.0d0)
GAM0(11, 25) = phii_jpt
GAM0(12, 25) = alpha
GAM0(15, 25) = 1.0d0
GAM0(16, 25) = -1.0d0*gamtil
GAM0(17, 25) = 1.0d0
GAM0(20, 25) = -1.0d0*bw + 1.0d0
GAM0(29, 25) = -1.0d0
GAM0(39, 25) = -1.0d0
GAM0(5, 26) = 1.0d0
GAM0(30, 26) = -1.0d0
GAM0(7, 27) = 1.0d0
GAM0(31, 27) = -1.0d0
GAM0(32, 28) = -1.0d0
GAM0(33, 29) = -1.0d0
GAM0(12, 30) = alpha - 1.0d0
GAM0(34, 30) = -1.0d0
GAM0(23, 31) = -1.0d0
GAM0(24, 32) = -1.0d0
GAM0(25, 33) = -1.0d0
GAM0(26, 34) = -1.0d0
GAM0(11, 35) = -1.0d0*beta*phii_jpt
GAM0(22, 35) = 1.0d0
GAM0(1, 36) = -1.0d0*beta*1.0/gz*(-1.0d0*delta + 1.0d0) + 1.0d0
GAM0(2, 37) = beta*gamtil
GAM0(21, 37) = -1.0d0*1.0/(-1.0d0*gamtil + 1.0d0)
GAM0(19, 38) = beta
GAM0(1, 39) = -1.0d0
GAM0(2, 39) = beta*gamtil
GAM0(10, 39) = -1.0d0
GAM0(11, 39) = -1.0d0*beta*phii_jpt
GAM0(21, 39) = -1.0d0*gamtil*1.0/(-1.0d0*gamtil + 1.0d0)
GAM0(22, 39) = 1.0d0
GAM0(6, 40) = beta*1.0/(beta*(-1.0d0*ap + 1.0d0) + 1.0d0)
GAM0(10, 40) = -1.0d0
GAM0(1, 41) = beta*1.0/gz*(-1.0d0*delta + 1.0d0)
GAM0(1, 42) = 1.0d0
GAM0(10, 42) = 1.0d0

GAM1(3, 1) = 1.0/gz*(-1.0d0*delta + 1.0d0)
GAM1(12, 1) = alpha
GAM1(15, 1) = 1.0d0
GAM1(2, 2) = -1.0d0*gamtil
GAM1(16, 2) = -1.0d0*gamtil
GAM1(24, 2) = -1.0d0
GAM1(11, 3) = phii_jpt
GAM1(17, 3) = 1.0d0
GAM1(25, 3) = -1.0d0
GAM1(4, 4) = -1.0d0
GAM1(26, 4) = -1.0d0
GAM1(5, 5) = -1.0d0*gam_rs
GAM1(6, 6) = -1.0d0*(-1.0d0*ap + 1.0d0)*1.0/(beta*(-1.0d0*ap + 1.0d0) + &
   1.0d0)
GAM1(18, 6) = -1.0d0*ap + 1.0d0
GAM1(20, 6) = aw - 1.0d0
GAM1(5, 7) = gamdy*(-1.0d0*gam_rs + 1.0d0)
GAM1(23, 7) = -1.0d0
GAM1(27, 23) = -0.85d0
GAM1(28, 24) = -1.0d0*rhoinv
GAM1(30, 26) = -1.0d0*rhoint
GAM1(31, 27) = -1.0d0*rhog
GAM1(32, 28) = -1.0d0*rhoelast
GAM1(33, 29) = -1.0d0*rhoelastw
GAM1(35, 35) = -1.0d0
GAM1(36, 36) = -1.0d0
GAM1(37, 37) = -1.0d0
GAM1(38, 38) = -1.0d0
GAM1(39, 39) = -1.0d0
GAM1(40, 40) = -1.0d0
GAM1(41, 41) = -1.0d0
GAM1(42, 42) = -1.0d0

PSI(27, 1) = -1.0d0
PSI(28, 2) = -1.0d0
PSI(29, 3) = -1.0d0
PSI(30, 4) = -1.0d0
PSI(31, 5) = -1.0d0
PSI(32, 6) = -1.0d0
PSI(33, 7) = -1.0d0
PSI(34, 8) = -1.0d0

PPI(35, 1) = -1.0d0
PPI(36, 2) = -1.0d0
PPI(37, 3) = -1.0d0
PPI(38, 4) = -1.0d0
PPI(39, 5) = -1.0d0
PPI(40, 6) = -1.0d0
PPI(41, 7) = -1.0d0
PPI(42, 8) = -1.0d0


call do_gensys(TT, CC, RR, fmat, fwt, ywt, gev, eu, loose, GAM0, GAM1, C, PSI, PPI, DIV)

!call write_array_to_file('TT.txt', TT)

QQ = 0.0d0
QQ(1, 1) = sdevliq**2.0d0
QQ(2, 2) = sdevinv**2.0d0
QQ(3, 3) = sdevtech**2.0d0
QQ(4, 4) = sdevint**2.0d0
QQ(5, 5) = sdevg**2.0d0
QQ(6, 6) = sdevelast**2.0d0
QQ(7, 7) = sdevelastw**2.0d0
QQ(8, 8) = sdeva**2.0d0

RR = matmul(RR, sqrt(QQ))
!call write_array_to_file('RR.txt', RR)

linsol%pp = TT(1:28,1:28)
linsol%sigma = RR(1:28,1:6)

!for shocks not included in nonlinear model, set innovation part of solution matrix to zero
if (linsol%nexogshock+linsol%nexogcont < linsol%nexog) then
   linsol%sigma(:,linsol%nexogshock+1:linsol%nexog-linsol%nexogcont) = 0.0d0
end if

end subroutine get_aimsolution

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: get_kalmanmatrices
!> @author Christopher Gust
!> @date 10-8-15
!> @brief Transform linear solution into Kalman matrices defined as:
!> $X_t = A X_{t-1} + Q \epsilon_{t}$.
!> $Y_t = G X_{t},$
!> where $X_t =(x_t,gdp_{t-1},c_{t-1},inv_{t-1},lab_{t-1},realw_{t-1})'$ and $Y_t$ are the obervables. 
!> Y_t = (dgdp_t,dc_t,dinv_t,dlab_t,drealw_t,dp_t,nomr_t)' 
!> @param[in] linsol Linear solution details.
!> @param[in] nobs Number of observables (hard-wired to 7)
!> @param[out] aa The "A" matrix used by the Kalman filter.
!> @param[out] qq The "Q" matrix used by the Kalman filter.
!> @param[out] gg The "G" matrix used by the Kalman filter.
!----------------------------------------------------------------------------------------------------------
subroutine get_kalmanmatrices(linsol,nobs,aa,qq,gg)

implicit none
type(linsoldetails), intent(in) :: linsol
integer, intent(in)           :: nobs
double precision, intent(out) :: aa(linsol%nvars+5,linsol%nvars+5)
double precision, intent(out) :: qq(linsol%nvars+5,linsol%nexog)
double precision, intent(out) :: gg(nobs,linsol%nvars+5)

aa = 0.0d0
qq = 0.0d0
gg = 0.0d0

aa(1:linsol%nvars,1:linsol%nvars) = linsol%pp
aa(linsol%nvars+1,7) = 1.0d0  !7 = gdp
aa(linsol%nvars+2,2) = 1.0d0  !2 = cc
aa(linsol%nvars+3,3) = 1.0d0  !3 = inv
aa(linsol%nvars+4,12) = 1.0d0 !12 = lab
aa(linsol%nvars+5,4) = 1.0d0  !4 = realw

qq(1:linsol%nvars,:) = linsol%sigma

gg(1,7) = 1.0d0
gg(1,linsol%nvars+1) = -1.0d0
gg(1,25) = 1.0d0
gg(2,2) = 1.0d0
gg(2,linsol%nvars+2) = -1.0d0
gg(2,25) = 1.0d0
gg(3,3) = 1.0d0
gg(3,linsol%nvars+3) = -1.0d0
gg(3,25) = 1.0d0
gg(4,12) = 1.0d0
gg(4,linsol%nvars+4) = -1.0d0
gg(5,4) = 1.0d0
gg(5,linsol%nvars+5) = -1.0d0
gg(5,25) = 1.0d0
gg(6,6) = 1.0d0
gg(7,9) = 1.0d0

end subroutine get_kalmanmatrices

!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: lindecrule_markov
!> @author Christopher Gust
!> @date 10-8-15
!> @brief Transform linear solution into representation that can be used to simulate shocks rather than innovations directly. 
!>        Used to get initial guess for nonlinear solution.
!> The transformed system has the form:
!> $y_t = A y_{t-1} + B s_{t}$,
!> where y_t are the endogenous variables excluding the shocks and s_t are the shocks.
!> @param[in] linsol Linear solution details.
!> @param[out] aalin Feedback part of linear decision rule that excludes shock processes.  
!> @param[out] bblin Linear decision rule matrix that maps effect of shocks (not innovations) on endogenous variables.
!---------------------------------------------------------------------------------------------
subroutine lindecrule_markov(linsol,aalin,bblin)

implicit none
type(linsoldetails), intent(in) :: linsol
double precision, intent(out)  :: aalin(linsol%nvars-linsol%nexog,linsol%nvars-linsol%nexog)
double precision, intent(out)  :: bblin(linsol%nvars-linsol%nexog,linsol%nexog)

integer :: i

aalin = 0.0d0
bblin = 0.0d0

aalin = linsol%pp(1:linsol%nvars-linsol%nexog,1:linsol%nvars-linsol%nexog)
do i = 1,linsol%nexogshock
   bblin(:,i) = linsol%sigma(1:linsol%nvars-linsol%nexog,i)/linsol%sigma(linsol%nvars-linsol%nexog+i,i)
end do

do i = 1,linsol%nexogcont
   bblin(:,linsol%nexog-i+1) = linsol%sigma(1:linsol%nvars-linsol%nexog,linsol%nexog-i+1)/linsol%sigma(linsol%nvars-i+1,linsol%nexog-i+1)
end do

end subroutine lindecrule_markov

!--------------------------------------------------------------------------------------------
! SUBROUTINE: decrlin
!> @author Christopher Gust
!> @date 10-8-15
!> @brief Linear decision rule -- returns endogenous variables and shocks given lagged endogenous values and innovations.
!> @param[in] endogvarm1 Lagged endogenous variables and shock values.
!> @param[in] innovations Innovations to the shocks.  
!> @param[in] linsol Linear solution details.
!> @param[out] endogvar Endogenous variables and shock values.
!---------------------------------------------------------------------------------------------
subroutine decrlin(endogvarm1,innovations,linsol,endogvar)

implicit none
type(linsoldetails), intent(in) :: linsol
double precision, intent(in)    ::  innovations(linsol%nexog)
double precision, intent(in)    :: endogvarm1(linsol%nvars)
double precision, intent(out)   :: endogvar(linsol%nvars)  

integer :: nvars
double precision :: xxm1(linsol%nvars), xx(linsol%nvars)
double precision :: exogpart(linsol%nvars)

xxm1 = 0.0d0
xx = 0.0d0 

xxm1(1:linsol%nvars-linsol%nexog) = endogvarm1(1:linsol%nvars-linsol%nexog)-linsol%endogsteady(1:linsol%nvars-linsol%nexog)
xxm1(linsol%nvars-linsol%nexog+1:linsol%nvars) = endogvarm1(linsol%nvars-linsol%nexog+1:linsol%nvars)

call dgemv('N', linsol%nvars, linsol%nvars, 1.0d0, linsol%pp, linsol%nvars, xxm1, 1, &
          0.0d0, xx, 1)
call dgemv('N', linsol%nvars, linsol%nexog, 1.0d0, linsol%sigma, linsol%nvars, innovations, 1, &
          0.0d0, exogpart, 1)
xx = xx + exogpart

endogvar(1:linsol%nvars-linsol%nexog) = xx(1:linsol%nvars-linsol%nexog)+linsol%endogsteady(1:linsol%nvars-linsol%nexog)
endogvar(linsol%nvars-linsol%nexog+1:linsol%nvars) = xx(linsol%nvars-linsol%nexog+1:linsol%nvars)

end subroutine decrlin

!--------------------------------------------------------------------------------------------
! SUBROUTINE: decrlin2
!> @author Christopher Gust
!> @date 10-30-15
!> @brief Linear decision rule -- returns endogenous variables and shocks in exact some format as nonlinear decision rule (decr). 
!> @param[in] endogvarm1 Lagged endogenous variables and shock values.
!> @param[in] innovations Innovations to the shocks.  
!> @param[in] linsol Linear solution details.
!> @param[out] endogvar Endogenous variables and shock values.
!---------------------------------------------------------------------------------------------
subroutine decrlin2(endogvarm1,innovations,linsol,endogvar)

implicit none
type(linsoldetails), intent(in) :: linsol
double precision, intent(in)    ::  innovations(linsol%nexog)
double precision, intent(in)    :: endogvarm1(linsol%nvars)
double precision, intent(out)   :: endogvar(linsol%nvars)  

integer :: nvars
double precision :: xxm1(linsol%nvars), xx(linsol%nvars)
double precision :: exogpart(linsol%nvars)

xxm1 = 0.0d0
xx = 0.0d0 

xxm1(1:linsol%nvars-linsol%nexog) = log(endogvarm1(1:linsol%nvars-linsol%nexog))-linsol%endogsteady(1:linsol%nvars-linsol%nexog)
xxm1(linsol%nvars-linsol%nexog+1:linsol%nvars) = endogvarm1(linsol%nvars-linsol%nexog+1:linsol%nvars)

call dgemv('N', linsol%nvars, linsol%nvars, 1.0d0, linsol%pp, linsol%nvars, xxm1, 1, &
          0.0d0, xx, 1)
call dgemv('N', linsol%nvars, linsol%nexog, 1.0d0, linsol%sigma, linsol%nvars, innovations, 1, &
          0.0d0, exogpart, 1)
xx = xx + exogpart

endogvar(1:linsol%nvars-linsol%nexog) = exp(xx(1:linsol%nvars-linsol%nexog)+linsol%endogsteady(1:linsol%nvars-linsol%nexog))
endogvar(linsol%nvars-linsol%nexog+1:linsol%nvars) = xx(linsol%nvars-linsol%nexog+1:linsol%nvars)

end subroutine decrlin2

!--------------------------------------------------------------------------------------------
! SUBROUTINE: simulate_kalman
!> @author Christopher Gust
!> @date 10-8-15
!> @brief Uses AIM decision rule to simulate model observables using kalman matrices.
!> @param[in] linsol Linear solution details.
!> @param[in] aa The "A" matrix used by the Kalman filter.
!> @param[in] qq The "Q" matrix used by the Kalman filter.
!> @param[in] gg The "G" matrix used by the Kalman filter.
!> @param[out] obsmean Mean of the observable variables.
!> @param[out] obsvar Std. deviation of the observable variables. 
!> @param[out] obscorr Correlations of the obervable variables. 
!---------------------------------------------------------------------------------------------
subroutine simulate_kalman(linsol,nobs,aa,qq,gg,obsmean,obsstd,obscorr)

use rng_serial
implicit none

type(linsoldetails), intent(in) :: linsol
integer, intent(in)            :: nobs
double precision, intent(in)   :: aa(linsol%nvars+5,linsol%nvars+5)
double precision, intent(in)   :: qq(linsol%nvars+5,linsol%nexog)
double precision, intent(in)   :: gg(nobs,linsol%nvars+5)
double precision, intent(out) :: obsmean(nobs)
double precision, intent(out) :: obsstd(nobs)
double precision, intent(out) :: obscorr(nobs*(nobs-1)/2)

integer, parameter :: capt = 1000 !Length of simulated dataset.    
integer            :: i,counter,displacement,llim,ulim,ttsim
integer            :: ic,ivar,jvar
double precision   :: modelvar(linsol%nvars+5,0:capt)
double precision   :: observables(nobs,capt)
logical            :: explosiveerror
double precision   :: xrandn(linsol%nexog,capt)
double precision   :: exogpart(linsol%nvars)

double precision   :: varsum(nobs),covsum
integer :: brng
integer :: method
integer :: errcode
integer :: iseed
type(fortress_random) :: rv
!type (vsl_stream_state) :: stream     

!get random uniforms
!iseed = 101293
!brng = vsl_brng_mt19937
!method = vsl_rng_method_gaussian_boxmuller 
!errcode = vslnewstream( stream,   brng,  iseed )
!errcode = vdrnggaussian( method, stream, linsol%nexog*capt, xrandn, 0.0d0, 1.0d0)
rv = fortress_random(seed=101293)
xrandn = rv%norm_rvs(linsol%nexog,capt)


!initialize endogenous variables and others 
modelvar = 0.0d0
counter = 0
displacement = 200

counterloop: do 
   explosiveerror = .false.
   llim = counter+1
   ulim = min(capt,counter+displacement)
   do ttsim = llim,ulim
       call dgemv('N', linsol%nvars+5, linsol%nvars+5, 1.0d0, aa, linsol%nvars+5, modelvar(:,ttsim-1), 1, &
           0.0d0, modelvar(:,ttsim), 1)
       call dgemv('N', linsol%nvars+5, linsol%nexog, 1.0d0, qq, linsol%nvars+5, xrandn(:,ttsim), 1, &
           0.0d0, exogpart, 1)
       modelvar(:,ttsim) = modelvar(:,ttsim) + exogpart

       call dgemv('N', nobs, linsol%nvars+5, 1.0d0, gg, nobs, modelvar(:,ttsim), 1, &
           0.0d0, observables(:,ttsim), 1)

   end do

   checkloop: do ttsim = llim,ulim
      explosiveerror = (isnan(modelvar(1,ttsim)) .eqv. .true.) 
      if (explosiveerror .eqv. .true.)  then
         counter = max(ttsim-50,0)
         !errcode = vdrnguniform( method, stream, linsol%nexog*capt, xrandn, 0.0d0, 1.0d0)
         xrandn = rv%norm_rvs(linsol%nexog,capt)
         if (explosiveerror .eqv. .true.) then 
            write(*,*) 'solution exploded at ', ttsim, ' vs ulim ', ulim
         end if
        
         if (counter == 0) then
            write(*,*) 'degenerated back to the beginning'
         end if
         exit checkloop
      end if
   end do checkloop

   if (explosiveerror .eqv. .false.) then
      counter = counter + displacement
   end if

   if (counter >= capt) then
      exit counterloop
   end if

end do counterloop

obsmean = sum(observables(:,1:capt),dim=2)/dble(capt)
varsum = 0.0d0
do i = 1,capt
   varsum = varsum + (observables(:,i)-obsmean)**2
end do
varsum = varsum/dble(capt-1)
obsstd = sqrt(varsum)

!could rewrite so that capt loop is outside the ivar and jvar loops
ic = 0
do ivar = 1,nobs
   do jvar = 1,nobs      
      if (jvar > ivar) then
         ic = ic + 1
         covsum = 0.0d0
         do i = 1,capt
            covsum = covsum + (observables(ivar,i)-obsmean(ivar))*(observables(jvar,i)-obsmean(jvar))
         end do
         covsum = covsum/dble(capt-1)
         obscorr(ic) = covsum/(obsstd(ivar)*obsstd(jvar))
      end if
   end do
end do

end subroutine simulate_kalman

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: linsoldeallocate
!> @author Christopher Gust
!> @date 10-8-15

!> @brief Deallocate the linear solution matrices.  
!> @param linsol Linear solution details (in\out).
!---------------------------------------------------------------------------------------------------------------

subroutine linsoldeallocate(linsol)  

implicit none 
type (linsoldetails), intent(inout) :: linsol

deallocate(linsol%pp) 
deallocate(linsol%sigma) 
deallocate(linsol%endogsteady)

end subroutine linsoldeallocate

end module linear_solution
