!------------------------------------------------------------------------------------------------------------------
! MODULE: simulate_model
!> @version 5.0
!> @author Christopher Gust
!
!DESCRIPTION: 
!> Subroutines for simulating the nonlinear model's dynamics.  
!-----------------------------------------------------------------------------------------------------------------
module simulate_model

use rng_serial
use model_details, only: decr, polydetails
use linear_solution, only: decrlin, linsoldetails

use fortress_random_t, only: fortress_random

implicit none

public :: simulate_data, simulate_irfs, compute_zlbstats
private :: euler_errorsfcn

contains

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: simulate_data
!> @author Christopher Gust
!> @date 11-30-15

!> @brief Simulate innovations and endogenous variables from the nonlinear model.
!> @param[in] capt Length of dataset.
!> @param[in] params Model parameters.  
!> @param[in]  poly Polynomial details.
!> @param[in] alphacoeff Values of polynomial coefficients at their solution.
!> @param[out] modeldata Endogenous variables and innovations.  
!> @param[in] seed Optional argument specifying the seed for the random draws.
!---------------------------------------------------------------------------------------------------------
subroutine simulate_data(capt,params,poly,linsol,alphacoeff,modeldata,nonlinearswitch,seed)

implicit none
integer, intent(in)             :: capt
type (polydetails), intent(in)  :: poly
type(linsoldetails), intent(in) :: linsol
logical, intent(in)             :: nonlinearswitch
double precision, intent(in)    :: alphacoeff(poly%nfunc*poly%ngrid,2*poly%ns)
double precision, intent(in)    :: params(poly%nparams)
double precision, intent(out)   :: modeldata(poly%nvars+2*poly%nexog,capt)
integer, intent(in), optional :: seed

integer            :: counter,displacement,llim,ulim,ttsim
double precision   :: endogvar(poly%nvars+poly%nexog,0:capt)
double precision   :: msvhigh(poly%nmsv),msvlow(poly%nmsv)
double precision   :: innovations(poly%nexog)
logical            :: explosiveerror,non_explosive
double precision   :: xrandn(poly%nexog,capt)
integer :: brng
integer :: method
integer :: errcode
integer :: iseed
type(fortress_random) :: rng
!type (vsl_stream_state) :: stream   

!get random normals
iseed = 101293


if (present(seed)) iseed = seed
rng = fortress_random(seed=iseed)
xrandn = rng%norm_rvs(poly%nexog, capt)


! brng = vsl_brng_mt19937
! method = vsl_rng_method_gaussian_boxmuller 
! errcode = vslnewstream( stream,   brng,  iseed )
! errcode = vdrnggaussian( method, stream, poly%nexog*capt, xrandn, 0.0d0, 1.0d0)

modeldata = 0.0d0
endogvar = 0.0d0
counter = 0
displacement = 400
if (nonlinearswitch .eqv. .true.) then
   endogvar(1:poly%nmsv,0) = exp(poly%endogsteady(1:poly%nmsv))
   msvhigh = exp(poly%endogsteady(1:poly%nmsv))*2.0d0
   msvlow = exp(poly%endogsteady(1:poly%nmsv))*0.01d0
else
   endogvar(1:poly%nmsv,0) = poly%endogsteady(1:poly%nmsv)
   msvhigh = log( exp(poly%endogsteady(1:poly%nmsv))*2.0d0 )
   msvlow = log( exp(poly%endogsteady(1:poly%nmsv))*0.01d0 )
end if

innovations = 0.0d0

counterloop: do 
   explosiveerror = .false.
   llim = counter+1
   ulim = min(capt,counter+displacement)
   do ttsim = llim,ulim
      innovations(1:poly%nexogshock) = xrandn(1:poly%nexogshock,ttsim)
      if (poly%nexogcont > 0) innovations(poly%nexog-poly%nexogcont+1:poly%nexog) = xrandn(poly%nexog-poly%nexogcont+1:poly%nexog,ttsim)
      if (nonlinearswitch .eqv. .true.) then
         call decr(endogvar(:,ttsim-1),innovations,params,poly,alphacoeff,endogvar(:,ttsim))
         modeldata(1:poly%nvars+poly%nexog,ttsim) = endogvar(:,ttsim)
      else
         call decrlin(endogvar(:,ttsim-1),innovations,linsol,endogvar(:,ttsim))
         modeldata(1:16,ttsim) = exp(endogvar(1:16,ttsim)+poly%endogsteady(1:16))
         modeldata(20,ttsim) = exp(endogvar(20,ttsim)+poly%endogsteady(20))
         modeldata((/17,18,19,21,22/),ttsim) = endogvar((/17,18,19,21,22/),ttsim)
         modeldata(poly%nvars+1:poly%nvars+poly%nexog,ttsim) = endogvar(poly%nvars+1:poly%nvars+poly%nexog,ttsim)
      end if 
      
      modeldata(poly%nvars+poly%nexog+1:poly%nvars+2*poly%nexog,ttsim) = innovations
   end do

   checkloop: do ttsim = llim,ulim
      non_explosive = ( all(endogvar(1:poly%nmsv,ttsim) .le. msvhigh) .and.  &
           all(endogvar(1:poly%nmsv,ttsim) .ge. msvlow) ) 
      explosiveerror = ( (non_explosive .eqv. .false.) .or. (isnan(endogvar(1,ttsim)) .eqv. .true.) )
      if (explosiveerror .eqv. .true.)  then
         counter = max(ttsim-350,0)
         xrandn = rng%norm_rvs(poly%nexog, capt)
         !errcode = vdrnggaussian( method, stream, poly%nexog*capt, xrandn, 0.0d0, 1.0d0)
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

end subroutine simulate_data

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: compute_zlbstats
!> @author Christopher Gust
!> @date 11-15-15

!> @brief Computes zlb probability and zlb duration data. 
!> @param[in] nomrdata Annual nominal rate series in percent.
!> @param[in] bigcapt Length of nomrdata.  
!> @param[in] capt_sample Sample length used to compute zlb probability.
!> @param[in] number_samples Number of samples of length capt_sample in bigcapt.
!> @param[out] number_zlbspells Number of zlb spells.
!> @param[out] zlbduration Array containing duration of each zlb spell.
!> @param[out] zlbfrequency Array containing zlb probability in each sample.
!---------------------------------------------------------------------------------------------------------
subroutine compute_zlbstats(nomrdata,bigcapt,capt_sample,number_samples,nspellmax,number_zlbspells,zlbduration,zlbfrequency)

implicit none
integer, intent(in)             :: bigcapt,capt_sample,number_samples,nspellmax
double precision, intent(in)    :: nomrdata(bigcapt)
integer, intent(out)            :: number_zlbspells
integer, intent(out)            :: zlbduration(nspellmax)
double precision, intent(out)   :: zlbfrequency(number_samples)

double precision, parameter :: zlbub = 0.25d0
integer :: ii,TTzlb,tt,zlbentrydate,zlbexitdate
double precision :: nomr_sample(capt_sample)

!compute zlb frequency statistics
zlbfrequency = 0.0d0
do ii = 1,number_samples
   nomr_sample = nomrdata(capt_sample*(ii-1)+1:capt_sample*ii)
   TTzlb = count(nomr_sample < zlbub)
   zlbfrequency(ii) = 100.0d0*dble(TTzlb)/dble(capt_sample)
end do

!count zlb spells and compute duration of each one
number_zlbspells = 0
zlbduration = 0
do tt = 2,bigcapt
   if ((nomrdata(tt-1) .ge. zlbub) .and. (nomrdata(tt) .lt. zlbub)) then
      number_zlbspells = number_zlbspells + 1
      zlbentrydate = tt
   end if
   

   if (number_zlbspells > 0) then
      if ((nomrdata(tt-1) .lt. zlbub) .and. (nomrdata(tt) .ge. zlbub)) zlbexitdate = tt
      zlbduration(number_zlbspells) = zlbexitdate-zlbentrydate
   end if
   
end do

end subroutine compute_zlbstats

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: simulate_irfs
!> @author Christopher Gust
!> @date 11-17-16

!> @brief Simulate IRFs by simulating an unshocked and shocked path which differ 
!> in their initial conditions (e.g. lagged shock values). 
!> Computes average deviation from unshocked path except for nominal and notional rates.
!> @param[in] capt Length of IRFs.
!> @param[in] ndraw Number of simulations used to do integration for IRFs.
!> @param[in] params Model parameters.
!> @param[in] poly Polynomial details for solution to the model.  
!> @param[in]  endogvarbas0 Initial conditions for endogenous and exogenous variables for unshocked path.
!> @param[in]  endogvarshk0 Initial conditions for endogenous and exogenous variables for shocked path.
!> @param[out] endogirf Nonlinear model's IRFs.
!---------------------------------------------------------------------------------------------------------
subroutine simulate_irfs(capt,ndraw,params,poly,alphacoeff,linsol,endogvarbas0,endogvarshk0,innov0,endogirf,linirf,euler_errors,neulererrors,&
     premirf)

implicit none 
integer, intent(in) :: capt,ndraw,neulererrors
type (polydetails), intent(in) :: poly
type(linsoldetails), intent(in) :: linsol
double precision, intent(in)    :: alphacoeff(poly%nfunc*poly%ngrid,2*poly%ns)
double precision, intent(in)   :: params(poly%nparams)
double precision, intent(in) :: endogvarbas0(poly%nvars+poly%nexog)
double precision, intent(in) :: endogvarshk0(poly%nvars+poly%nexog)
double precision, intent(in) :: innov0(poly%nexog)
double precision, intent(out) :: endogirf(poly%nvars+poly%nexog+2,capt)  !+2 should match nwish_level
double precision, intent(out) :: linirf(poly%nvars+poly%nexog+2,capt)  !+2 shoud match nwish_level
double precision, intent(out) :: euler_errors(2*neulererrors,capt)  
double precision, intent(out) :: premirf(2,capt)

integer, parameter :: nwish_level = 2
logical, parameter :: eulererrorswitch = .true.
integer            :: i,j,ttsim,j_used
integer            :: wishlist_level(nwish_level),med_rank,info
logical            :: nlerrorswitch
double precision   :: endogvarbas(poly%nvars+poly%nexog,0:capt)
double precision   :: endogvarshk(poly%nvars+poly%nexog,0:capt)
double precision   :: endoglinbas(poly%nvars+poly%nexog,0:capt)
double precision   :: endoglinshk(poly%nvars+poly%nexog,0:capt)
double precision   :: endoglinshkm1_exp(poly%nvars+poly%nexog)
double precision   :: endoglinshk_exp(poly%nvars+poly%nexog)
double precision   :: varsort(ndraw)
double precision   :: levelresponse(nwish_level*capt,ndraw)
double precision   :: levellin(nwish_level*capt,ndraw)
double precision   :: innovations(poly%nexog),innovplus(poly%nexog)
double precision   :: xrandn(poly%nexog,capt*ndraw)
double precision   :: errors_nl(neulererrors),errors_lin(neulererrors)
double precision   :: premiumbas(capt),premiumshk(capt)
!double precision   :: nomrdata(capt*ndraw)
!character (LEN = 200) :: outputfile  
integer :: brng
integer :: method
integer :: errcode
integer :: iseed
type(fortress_random) :: rng
!type (vsl_stream_state) :: stream    

iseed = 101293
rng = fortress_random(seed=iseed)
xrandn = rng%norm_rvs(linsol%nexog, capt*ndraw)

! brng = vsl_brng_mt19937
! method = vsl_rng_method_gaussian_boxmuller 
! errcode = vslnewstream( stream,   brng,  iseed )
! errcode = vdrnggaussian( method, stream, poly%nexog*capt*ndraw, xrandn, 0.0d0, 1.0d0)
 
endogirf = 0.0d0
linirf = 0.0d0
euler_errors = 0.0d0
levelresponse = 0.0d0
levellin = 0.0d0
innovations = 0.0d0
premirf = 0.0d0
premiumshk = 1.0d0
premiumbas = 1.0d0
!nomrdata = 0.0d0
!9=nominal rate, 5=notional rate
wishlist_level = (/9,5/)
j_used = 0
do j = 1,ndraw
   endogvarbas = 0.0d0
   endogvarbas(:,0) = endogvarbas0
   endogvarshk = 0.0d0
   endogvarshk(:,0) = endogvarshk0
   endoglinbas = 0.0d0
   endoglinbas(1:poly%nvars,0) = log(endogvarbas0(1:poly%nvars))
   endoglinbas(poly%nvars+1:poly%nvars+poly%nexog,0) = endogvarbas0(poly%nvars+1:poly%nvars+poly%nexog)
   endoglinshk = 0.0d0
   endoglinshk(1:poly%nvars,0) = log(endogvarshk0(1:poly%nvars))
   endoglinshk(poly%nvars+1:poly%nvars+poly%nexog,0) = endogvarshk0(poly%nvars+1:poly%nvars+poly%nexog)
   simloop: do ttsim = 1,capt   
      innovations(1:poly%nexogshock) = xrandn(1:poly%nexogshock, capt*(j-1)+ttsim)
      if (poly%nexogcont > 0) innovations(poly%nexog-poly%nexogcont+1:poly%nexog) = xrandn(poly%nexog-poly%nexogcont+1:poly%nexog,ttsim)
      call decr(endogvarbas(:,ttsim-1),innovations,params,poly,alphacoeff,endogvarbas(:,ttsim))
      call decrlin(endoglinbas(:,ttsim-1),innovations,linsol,endoglinbas(:,ttsim))
      if (ttsim .eq. 1) then
          innovplus = innovations + innov0
      else
          innovplus = innovations
      end if
      call decr(endogvarshk(:,ttsim-1),innovplus,params,poly,alphacoeff,endogvarshk(:,ttsim)) 
      call decrlin(endoglinshk(:,ttsim-1),innovplus,linsol,endoglinshk(:,ttsim))

      !call calc_premium(endogvarshk(:,ttsim),params,poly,alphacoeff,premiumshk(ttsim))
      !call calc_premium(endogvarbas(:,ttsim),params,poly,alphacoeff,premiumbas(ttsim))
      
      !IRFs
      if ( any(isnan(endogvarshk(:,ttsim))) .or. any(isnan(endogvarbas(:,ttsim)))  ) then
         !write(*,*) 'failure at draw = ', j
         exit simloop 
      else
         
         !compute absolute value of euler errors
         if (eulererrorswitch .eqv. .true.) then
            nlerrorswitch = .true.
            errors_nl = euler_errorsfcn(neulererrors,endogvarshk(:,ttsim-1),endogvarshk(:,ttsim),params,poly,alphacoeff,nlerrorswitch,linsol)
            if (any(isnan(errors_nl))) then
               exit simloop
            else 
               euler_errors(1:neulererrors,ttsim) = euler_errors(1:neulererrors,ttsim) + errors_nl
            end if
            

            nlerrorswitch = .false.
            endoglinshkm1_exp(1:poly%nvars) = exp(endoglinshk(1:poly%nvars,ttsim-1))
            endoglinshkm1_exp(poly%nvars+1:poly%nvars+poly%nexog) = endoglinshk(poly%nvars+1:poly%nvars+poly%nexog,ttsim-1)
            endoglinshk_exp(1:poly%nvars) = exp(endoglinshk(1:poly%nvars,ttsim))
            endoglinshk_exp(poly%nvars+1:poly%nvars+poly%nexog) = endoglinshk(poly%nvars+1:poly%nvars+poly%nexog,ttsim)
            errors_lin = euler_errorsfcn(neulererrors,endoglinshkm1_exp,endoglinshk_exp,params,poly,alphacoeff,nlerrorswitch,linsol)
            if (any(isnan(errors_lin))) then
               exit simloop
            else 
               euler_errors(neulererrors+1:2*neulererrors,ttsim) = euler_errors(neulererrors+1:2*neulererrors,ttsim) + errors_lin   
            end if
         end if

         if (ttsim .eq. capt) then
            j_used = j_used + 1
         end if
         
         endogirf(1:poly%nvars+poly%nexog,ttsim) = endogirf(1:poly%nvars+poly%nexog,ttsim) + 100.0d0*&
              log(endogvarshk(:,ttsim)/endogvarbas(:,ttsim))  
         linirf(1:poly%nvars+poly%nexog,ttsim) = linirf(1:poly%nvars+poly%nexog,ttsim) + 100.0d0*&
              (endoglinshk(:,ttsim)-endoglinbas(:,ttsim)) 
         premirf(1,ttsim) = premirf(1,ttsim) + 100.0d0*log(premiumshk(ttsim)/premiumbas(ttsim))

         do i = 1,nwish_level   
            levelresponse((i-1)*capt+ttsim,j_used) = 400.0d0*log(endogvarshk(wishlist_level(i),ttsim))
            levellin((i-1)*capt+ttsim,j_used) = 400.0d0*endoglinshk(wishlist_level(i),ttsim)
         end do
      end if
   end do simloop
end do

endogirf = endogirf/dble(j_used)
linirf = linirf/dble(j_used)
premirf = premirf/dble(j_used)
euler_errors = euler_errors/dble(j_used)

!med_rank = int(0.5d0*dble(j_used)+1.0d0)
! do ttsim = 1,capt 
!    do i = 1,nwish_level
!       varsort(1:j_used) = levelresponse((i-1)*capt+ttsim,1:j_used)
!       call dlasrt('I',ndraw,varsort,info)
!       endogirf(poly%nvars+poly%nexog+i,ttsim) = varsort(med_rank)    
!       varsort = levellin((i-1)*capt+ttsim,:)
!       call dlasrt('I',ndraw,varsort,info)
!       linirf(poly%nvars+poly%nexog+i,ttsim) = varsort(med_rank)  
!    end do
! end do

end subroutine simulate_irfs

!--------------------------------------------------------------------------------------------------------
!FUNCTION: euler_errors
!> @author Christopher J. Gust
!> @date 12-4-15

!> @brief Compute euler errors given the lagged and current endogenous variables. 
!> @param[in] endogvar Current values of endogenous variables. 
!> @param[in] param Model parameters
!> @param[in] polysolution Polynomial details.  
!> @param[in] currentshocks Index of current value of shocks.  
!> @return Value of Euler errors. 
!------------------------------------------------------------------------------------------------------------------
function euler_errorsfcn(neulererrors,endogvarm1,endogvar,params,poly,alphacoeff,nlerrorswitch,linsol)

implicit none
type(polydetails), intent(in) :: poly
type(linsoldetails), intent(in) :: linsol
integer, intent(in) :: neulererrors
double precision, intent(in)    :: params(poly%nparams)
double precision, intent(in) :: endogvarm1(poly%nvars+poly%nexog),endogvar(poly%nvars+poly%nexog)
double precision, intent(in) :: alphacoeff(poly%nfunc*poly%ngrid,2*poly%ns)
logical, intent(in) :: nlerrorswitch
double precision :: euler_errorsfcn(neulererrors)

integer :: ss,i
double precision :: beta,pibar,gz,gamma,psil,sigmal,phi,phiw,ep,epw,ap,aw,bw,lamhp,alpha,delta,phii
double precision :: sigmaa,shrgy,gam_rs,gam_dp,gam_dy,gam_xhp,gss,rss,rkss
double precision :: utilcostp,dwtildem1,gzwage,dw,xhp,vp,dptildem1,gshk,aayy,utilcost
double precision :: techshk,techshkp,invshk,invshkp,liqshk
double precision :: endogvarp(poly%nvars+poly%nexog)
double precision :: lendogvar(poly%nvars+poly%nexog),lendogvarp(poly%nvars+poly%nexog)
double precision :: ev(6),exp_eul(6),exp_var(6)
double precision :: innovations(poly%nexog)

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
gam_xhp = params(21) 
gam_dy = params(22)
shrgy = params(23)

rkss = params(3)/params(1)-1.0d0+params(16)
exp_eul = 0.0d0
exp_var = 0.0d0
innovations = 0.0d0
do ss = 1,poly%nquad
   innovations(1:poly%nexogshock) = poly%ghnodes(:,ss) 
   if (nlerrorswitch .eqv. .true.) then
      call decr(endogvar,innovations,params,poly,alphacoeff,endogvarp)
   else
      lendogvar(1:poly%nvars) = log(endogvar(1:poly%nvars))
      lendogvar(poly%nvars+1:poly%nvars+poly%nexog) = endogvar(poly%nvars+1:poly%nvars+poly%nexog)
      call decrlin(lendogvar,innovations,linsol,lendogvarp)
      endogvarp(1:poly%nvars) = exp(lendogvarp(1:poly%nvars))
      endogvarp(poly%nvars+1:poly%nvars+poly%nexog) = lendogvarp(poly%nvars+1:poly%nvars+poly%nexog)
   end if
   techshkp = exp(endogvarp(25))
   invshkp = exp(endogvarp(24))
   ev(1) = endogvarp(10)/(endogvarp(6)*techshkp)
   utilcostp = (rkss/params(18))*(exp(params(18)*(endogvarp(13)-1.0d0))-1.0d0)
   ev(2) = (endogvarp(10)/techshkp)*( (endogvarp(15)*endogvarp(13))-utilcostp+(1.0d0-params(16))*endogvarp(11) )
   ev(3) = endogvarp(10)*(endogvarp(18)-1.0d0)*endogvarp(18)*endogvarp(7)
   ev(4) = (endogvarp(19)-1.0d0)*endogvarp(19) 
   ev(5) = endogvarp(16)/techshkp 
   ev(6) = endogvarp(10)*endogvarp(11)*invshkp*(endogvarp(17)-1.0d0)*endogvarp(17)*endogvarp(17)
   exp_var = exp_var + poly%ghweights(ss)*ev   
end do

!ep and epw shock not incorporated, need to update if you want markup shocks
techshk = exp(endogvar(25))
liqshk = exp(endogvar(23))
invshk = exp(endogvar(24))
exp_eul(1) =  (params(1)/params(3))*liqshk*endogvar(9)*exp_var(1)
exp_eul(2) = (params(1)/params(3))*exp_var(2)/endogvar(10)
exp_eul(3) =  params(1)*exp_var(3)/(endogvar(10)*endogvar(7))+&
     (params(9)/params(7))*( endogvar(14)-(params(9)-1.0d0)/params(9) ) 
exp_eul(4) = params(1)*exp_var(4)+(params(10)/params(8))*endogvar(10)*endogvar(12)*&
           ( (params(4)*endogvar(12)**params(6)/endogvar(10))-((params(10)-1.0d0)/&
           params(10))*endogvar(4) ) 
exp_eul(5) = exp_var(5)
exp_eul(6) = (params(1)/params(3))*exp_var(6)/endogvar(10) - 0.5d0*endogvar(11)*invshk*(endogvar(17)-1.0d0)*(endogvar(17)-1.0d0)

euler_errorsfcn = 0.0d0
euler_errorsfcn(1) = log(exp_eul(1)/endogvar(10)) 
euler_errorsfcn(2) = log(exp_eul(2)/endogvar(11))
euler_errorsfcn(3) = exp_eul(3)-(endogvar(18)-1.0d0)*endogvar(18)
euler_errorsfcn(4) = exp_eul(4)-(endogvar(19)-1.0d0)*endogvar(19)
euler_errorsfcn(5) = log(exp_eul(5)/endogvar(21))

euler_errorsfcn(6) = log( 1.0d0 + (1/params(18))*log(endogvar(15)/rkss) )-log(endogvar(13))

euler_errorsfcn(7) = endogvar(1)-(1.0d0-delta)*(endogvarm1(1)/(gz*techshk)) -&
     invshk*endogvar(3)*(1.0d0-(phii/2.0d0)*(endogvar(17)-1.0d0)**2) !capital acc
euler_errorsfcn(8) = endogvar(2)-(gamma*endogvarm1(2)/(gz*techshk)+1.0d0/endogvar(16))  !defines cc using muc
euler_errorsfcn(9) = endogvar(11)*invshk*(endogvar(17)-1.0d0)*endogvar(17) + (1.0d0/phii)*(1.0d0-endogvar(11)*invshk)-exp_eul(6) !investment demand

dwtildem1 = (pibar**aw)*(endogvarm1(6)**(1.0d0-aw))
gzwage = gz*techshk**(1.0d0-bw)
dw = endogvar(19)*dwtildem1*gzwage
euler_errorsfcn(10) = endogvar(4)-(endogvarm1(4)*dw)/(gz*techshk*endogvar(6))

xhp = alpha*log(endogvar(13)) + (1-alpha)*(log(endogvar(12))-poly%endogsteady(12))
rss = gz*pibar/beta
euler_errorsfcn(11) = log(rss)+gam_rs*(log(endogvarm1(5))-log(rss)) + (1.0d0-gam_rs)*( gam_dp*(log(endogvar(6)/pibar)) +&
     gam_dy*log(endogvar(7)*techshk/endogvarm1(7)) + gam_xhp*xhp ) + endogvar(26) - log(endogvar(5))

vp = endogvar(18)  
dptildem1 = (pibar**ap)*(endogvarm1(6)**(1.0d0-ap))
euler_errorsfcn(12) = vp*dptildem1 - endogvar(6)

gss = 1.0d0/(1.0d0-shrgy)
gshk = exp( log(gss) + endogvar(27) )
aayy = 1.0d0/gshk-(phi/2.0d0)*(vp-1.0d0)**2
utilcost = (rkss/sigmaa)*(exp(sigmaa*(endogvar(13)-1.0d0))-1.0d0)
euler_errorsfcn(13) = (1.0d0/aayy)*( endogvar(2)+endogvar(3) + utilcost*(endogvarm1(1)/(gz*techshk)) ) - endogvar(7)

euler_errorsfcn(14) = endogvar(7)**(1.0d0/(1.0d0-alpha))*(endogvar(13)*endogvarm1(1)/(gz*techshk))**(alpha/(alpha-1.0d0)) - endogvar(12) !labor
euler_errorsfcn(15) = endogvar(14)-endogvar(4)*endogvar(12)/((1.0d0-alpha)*endogvar(7))  !mc
euler_errorsfcn(16) = (alpha/(1.0d0-alpha))*(endogvar(4)*endogvar(12)*gz*techshk/(endogvar(13)*endogvarm1(1))) - endogvar(15) !rentalk
euler_errorsfcn(17) =  endogvar(10) + (gamma/gz)*beta*endogvar(21)  - endogvar(16) !muc

!average absolute value of euler errors
do i = 1,neulererrors-1
   euler_errorsfcn(i) = abs(euler_errorsfcn(i))
   euler_errorsfcn(neulererrors) = euler_errorsfcn(neulererrors) + euler_errorsfcn(i)/dble(17)
end do

end function euler_errorsfcn

end module simulate_model
