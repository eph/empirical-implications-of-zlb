!------------------------------------------------------------------------------------------------------------------
! MODULE: model_details
!> @version 5.2
!> @author Christopher Gust
!> @brief NK model similar to CEE or SW.
!-----------------------------------------------------------------------------------------------------------------

module model_details

use polydef

public :: decr,steadystate,decr_euler,get_shockdetails, calc_premium
private :: intermediatedec

contains

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: decr
!> @author Christopher Gust
!> @date 6-28-16

!> @brief The decision rule -- returns endogenous variables and shocks given lagged endogenous values and innovations.
!> @param[in] endogvarm1 Lagged endogenous variables and shock values.
!> @param[in] innovations Innovations to the shocks.  
!> @param[in] params Model parameters.  
!> @param[in] poly Polynomial details.
!> @param[in] alphacoeff Polynomial coefficients.
!> @param[out] endogvar Endogenous variables and shock values
!----------------------------------------------------------------------------------------------------------
subroutine decr(endogvarm1,innovations,params,poly,alphacoeff,endogvar)   

implicit none
type(polydetails), intent(in)   :: poly
double precision, intent(in)    ::  innovations(poly%nexog)
double precision, intent(in)    :: params(poly%nparams)
double precision, intent(in)    :: endogvarm1(poly%nvars+poly%nexog)
double precision, intent(in)    :: alphacoeff(poly%nfunc*poly%ngrid,2*poly%ns)
double precision, intent(out)   :: endogvar(poly%nvars+poly%nexog)  

integer          :: i,stateindex,stateindexplus,stateindex0,stateindex1,ifunc,nmsvplus
integer          :: shockindexall(poly%nexog-poly%nexogcont),shockindex(poly%nexogshock)
integer          :: shockindex_inter(poly%nexogshock)
double precision :: sdevtech,rhog,sdevg,rhoinv,sdevinv,rholiq,sdevliq,rhoint,sdevint,lmsv(poly%nmsv+poly%nexogcont)
double precision :: rhoa,sdeva
double precision :: currentshockvalues(poly%nexog)
double precision :: funcmat(poly%nfunc,poly%ninter)
double precision :: funcmatplus(poly%nfunc,poly%ninter)
double precision :: weightvec(poly%ninter),weighttemp(poly%nexogshock)
double precision :: funcapp(poly%nfunc),funcapp_plus(poly%nfunc)
double precision :: xx(poly%nmsv),polyvec(poly%ngrid)
double precision :: omegapoly,prod_sd
double precision :: llabss
double precision, parameter :: omegaweight = 100000.0d0
logical          :: zlbintermediate

nmsvplus = poly%nmsv+poly%nexogcont

!shock parameters
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

!update shocks (all shocks in deviation from SS)
currentshockvalues(1) = rholiq*endogvarm1(23) + sdevliq*innovations(1)
currentshockvalues(2) = rhoinv*endogvarm1(24) + sdevinv*innovations(2)
currentshockvalues(3) = sdevtech*innovations(3)
currentshockvalues(4) = rhoint*endogvarm1(26) + sdevint*innovations(4)
currentshockvalues(5) =  rhog*endogvarm1(27) + sdevg*innovations(5)
currentshockvalues(6) = rhoa*endogvarm1(28) + sdeva*innovations(6)

!find position of shocks for interpolation
do i = 1,poly%nexogshock
   if (currentshockvalues(i) .lt. poly%shockbounds(i,1)) then  !extrapolation below
      shockindex(i) = 1
   else if (currentshockvalues(i) .ge. poly%shockbounds(i,2)) then  !extrapolation above
      shockindex(i) = int( (poly%shockbounds(i,2)-poly%shockbounds(i,1))/poly%shockdistance(i) )
   else  
      shockindex(i) = int(1.0d0+(currentshockvalues(i)-poly%shockbounds(i,1))/poly%shockdistance(i))  !interpolation case
   end if
end do 

!loop for interpolating between shocks
shockindexall = 1
shockindexall(1:poly%nexogshock) = shockindex
stateindex0 = exogposition(shockindexall,poly%nshockgrid,poly%nexog-poly%nexogcont)
shockindexall(1:poly%nexogshock) = shockindex + poly%interpolatemat(:,poly%ninter)
stateindex1 = exogposition(shockindexall,poly%nshockgrid,poly%nexog-poly%nexogcont)
funcmat = 0.0d0
weightvec = 0.0d0
lmsv(1:poly%nmsv) = log(endogvarm1(1:poly%nmsv))
if (poly%nexogcont > 0) lmsv(poly%nmsv+1:poly%nmsv+poly%nexogcont) = currentshockvalues(poly%nexog-poly%nexogcont+1:poly%nexog)
xx = msv2xx(lmsv,nmsvplus,poly%slopeconmsv)

polyvec = smolyakpoly(nmsvplus,poly%ngrid,poly%nindplus,poly%indplus,xx)
prod_sd = product(poly%shockdistance)
do i = 1,poly%ninter
   shockindex_inter = shockindex + poly%interpolatemat(:,i)
   shockindexall(1:poly%nexogshock) = shockindex_inter
   stateindex = exogposition(shockindexall,poly%nshockgrid,poly%nexog-poly%nexogcont)
   stateindexplus = stateindex+poly%ns  
   do ifunc = 1,poly%nfunc
      funcmat(ifunc,i) = dot_product(alphacoeff((ifunc-1)*poly%ngrid+1:ifunc*poly%ngrid,stateindex),polyvec)
      funcmatplus(ifunc,i) = dot_product(alphacoeff((ifunc-1)*poly%ngrid+1:ifunc*poly%ngrid,stateindexplus),polyvec)
   end do
   weighttemp = dble(1-poly%interpolatemat(:,i))*(currentshockvalues(1:poly%nexogshock)-poly%exoggrid(1:poly%nexogshock,stateindex0)) +&
        dble(poly%interpolatemat(:,i))*(poly%exoggrid(1:poly%nexogshock,stateindex1)-currentshockvalues(1:poly%nexogshock))
   weightvec(poly%ninter-i+1) = product(weighttemp)/prod_sd
end do
call dgemv('N', poly%nfunc, poly%ninter, 1.0d0, funcmat, poly%nfunc, weightvec, 1, &
          0.0d0, funcapp, 1)

llabss = poly%endogsteady(12)
zlbintermediate = .false.  !start with evaluation of 1 poly case (omegapoly and 2nd funcapp irrelevant)
call intermediatedec(poly%nparams,poly%nvars,poly%nexog,poly%nfunc,endogvarm1,&
     currentshockvalues,params,funcapp,llabss,endogvar,omegapoly,funcapp,zlbintermediate)

if ( (endogvar(5) .lt. 1.0d0) .and. (poly%zlbswitch .eqv. .true.) ) then  !zlb case
   zlbintermediate = .true.
   omegapoly = exp(omegaweight*log(endogvar(5))) !now omegapoly and funcapp_plus relevant
   call dgemv('N', poly%nfunc, poly%ninter, 1.0d0, funcmatplus, poly%nfunc, weightvec, 1, &
          0.0d0, funcapp_plus, 1)

   call intermediatedec(poly%nparams,poly%nvars,poly%nexog,poly%nfunc,&
           endogvarm1,currentshockvalues,params,funcapp,llabss,endogvar,omegapoly,funcapp_plus,zlbintermediate)
   endogvar(9) = 1.0d0
end if

end subroutine decr

!--------------------------------------------------------------------------------------------------------
! FUNCTION: steadystate
!> @author Christopher Gust
!> @date 6-28-15
!> @brief Returns the parameters used by AIM and computes the steady state of the NK model.
!> @param[in] params Model parameters.
!> @param[in] nvars Number of steady state variables.
!> @param[in] nparams Number of model parameters.  
!> @return Steady state values of model's variables.
!----------------------------------------------------------------------------------------------------------
function steadystate(params,nvars,nparams)  

implicit none
integer, intent(in) :: nvars,nparams
double precision, intent(in) :: params(nparams)
double precision :: steadystate(nvars)

double precision :: beta,pibar,gz,psil,gamma,sigmal,phi,phiw,ep,epw,ap,aw,bw,lamhp,alpha,delta,phii,sigmaa,gam_rs,gam_dp
double precision :: gamxhp,gamdy,shrgy
double precision :: gg,gamtil,mc,k2yrat,shriy,shrcy,labss,kappaw,kappap,kss,gdpss,invss,phii_jpt,css,rwss,mucss,lamss
double precision :: rss,rkss

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

!fix bw = aw (doesn't matter for SS but done anyway for consistency)
bw = aw

!composite parameters and steady state calculations
gg = 1.0d0/(1.0d0-shrgy)
gamtil = gamma/gz
mc = (ep-1.0d0)/ep
k2yrat = ((mc*alpha)/(gz/beta-(1.0d0-delta)))*gz
shriy = (1.0d0-(1.0d0-delta)/gz)*k2yrat
shrcy = 1.0d0-shrgy-shriy
labss = ( ((epw-1.0d0)/epw)*(1.0d0-alpha)*(1.0d0-beta*gamtil)*((ep-1.0d0)/ep)*(1.0d0/(psil*(1.0d0-gamtil)))*&
     (1.0d0/shrcy) )**(1.0d0/(sigmal+1.0d0))
kappaw = ((1.0d0-gamtil)/(1.0d0-beta*gamtil))*epw*psil*labss**(1.0d0+sigmal)/phiw
kappap = (ep-1.0d0)/(phi*(1.0d0+beta*(1.0d0-ap)))
!print*,'kappap',kappap, kappap/(ep-1.0d0)
kss = labss*(gz**(alpha/(alpha-1.0d0)))*k2yrat**(1.0d0/(1.0d0-alpha))
gdpss = (kss/gz)**alpha*labss**(1.0d0-alpha)
invss = shriy*gdpss
phii_jpt = phii/invss
css = shrcy*gdpss
rwss = (1.0d0-alpha)*mc*gdpss/labss
mucss = (1.0d0/css)*(1.0d0/(1.0d0-gamtil))
lamss = mucss*(1.0d0-beta*gamtil)
rss = gz*pibar/beta
rkss = gz/beta-1.0d0+delta

!(1) cap, (2) cc, (3) inv, (4) rw, (5) notr, (6) dp, (7) gdp, (8) xhp, (9) nomr, (10) lam, (11) qq, (12) lab, (13) util, (14) mc, (15) rentalk, (16) muc
!(17) vi, (18) vp, (19) vw, (20) dw, (21) bc, (22) bi, (23) liqshk, (24) invshk, (25) techshk, (26) intshk, (27) gshk, (28) ashk
steadystate = (/log(kss),log(css),log(invss),log(rwss),log(rss),log(pibar),log(gdpss),0.0d0,log(rss),log(lamss),&
     0.0d0,log(labss),0.0d0,log(mc),log(rkss),log(mucss),0.0d0,0.0d0,0.0d0,log(pibar*gz),log(mucss),0.0d0,0.0d0,0.0d0,0.0d0,log(gg),0.0d0/)

end function steadystate


!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: intermediatedec
!> @author Christopher Gust
!> @date 6-28-15

!> @brief Returns endogenous variables and shock values given their lagged values and polynomial approximation.
!> @param[in] nparams  Number of structural model parameters. 
!> @param[in] nvars Number of endogenous variables.
!> @param[in] nexog Number of exogenous shocks.
!> @param[in] nfunc Number of functions approximated.
!> @param[in] endogvarm1 Lagged endogenous variables and shock values.
!> @param[in] currentshockvalues Current shock values.
!> @param[in] params Model parameters.  
!> @param[in] polyvar Values of polynomial functions in normal times.
!> @param[in] llabss Steady state value of labor.
!> @param[out] endogvar Endogenous variables and shock values
!> @param[in] omegapoly Weight on ZLB polynomials.
!> @param[in] polyvarplus Values of polynomial functions in zlb times.
!> @param[in] zlbintermediate Indicator for whether there are two regimes for the polynomial or not.
!----------------------------------------------------------------------------------------------------------
subroutine intermediatedec(nparams,nvars,nexog,nfunc,endogvarm1,currentshockvalues,params,polyvar,llabss,endogvar,omegapoly,&
     polyvarplus,zlbintermediate)

implicit none 
integer, intent(in)             :: nparams,nvars,nexog,nfunc
double precision, intent(in)    :: params(nparams)
double precision, intent(in)    :: endogvarm1(nvars+nexog)  
double precision, intent(in)    :: currentshockvalues(nexog)
double precision, intent(in)    :: polyvar(nfunc)
double precision, intent(in)    :: llabss
double precision, intent(out)   :: endogvar(nvars+nexog)  
double precision, intent(in)    :: omegapoly
double precision, intent(in)    :: polyvarplus(nfunc)
logical, intent(in)             :: zlbintermediate

!Declare local variables
double precision :: beta,pibar,gz,gamma,psil,sigmal,phi,phiw,ep,epw,ap,aw,bw,lamhp,alpha,delta,phii
double precision :: sigmaa,shrgy,gam_rs,gam_dp,gam_dy,gam_xhp,gss,lrss,rkss
double precision :: dpm1,rwm1,capm1,invm1,ccm1,gdpm1,notrm1
double precision :: dptildem1,dwtildem1,gzwage
double precision :: invshk,gshk,techshk,rrshk,util,utilcost,ashk
double precision :: lam,qq,bp,bww,bc,bi,vp,vw,dp,dw,muc,cc,vi,inv,aayy,gdp,notr,rw,lab,mc,rentalk,cap,nomr
double precision :: xhp,cquad_vi

!parameters
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

!fix bw = aw
bw = aw

invshk = exp( currentshockvalues(2) )
techshk = exp( currentshockvalues(3) ) 
rrshk = currentshockvalues(4)
gss = 1.0d0/(1.0d0-shrgy)
gshk = exp( log(gss) + currentshockvalues(5) )
ashk = exp( currentshockvalues(6) )

capm1 = endogvarm1(1)
ccm1 = endogvarm1(2)
invm1 = endogvarm1(3)
rwm1 = endogvarm1(4)
notrm1 = endogvarm1(5)
dpm1 = endogvarm1(6)
gdpm1 = endogvarm1(7)

if (zlbintermediate .eqv. .true.) then
   lam = omegapoly*exp(polyvar(1)) + (1.0d0-omegapoly)*exp(polyvarplus(1))
   qq = omegapoly*exp(polyvar(2)) + (1.0d0-omegapoly)*exp(polyvarplus(2))
   bp = omegapoly*polyvar(3) + (1.0d0-omegapoly)*polyvarplus(3)
   bww = omegapoly*polyvar(4) + (1.0d0-omegapoly)*polyvarplus(4)
   bc = omegapoly*exp(polyvar(5)) + (1.0d0-omegapoly)*exp(polyvarplus(5))
   bi = omegapoly*polyvar(6) + (1.0d0-omegapoly)*polyvarplus(6)
   util = omegapoly*exp(polyvar(7)) + (1.0d0-omegapoly)*exp(polyvarplus(7))
else
   lam = exp(polyvar(1)) 
   qq = exp(polyvar(2)) 
   bp = polyvar(3)
   bww = polyvar(4)
   bc = exp(polyvar(5)) 
   bi = polyvar(6)
   util = exp(polyvar(7))
end if
   
vp = (sqrt(1.0d0+4.0d0*bp)+1.0d0)/2.0d0
vw = (sqrt(1.0d0+4.0d0*bww)+1.0d0)/2.0d0
dptildem1 = (pibar**ap)*(dpm1**(1.0d0-ap))
dwtildem1 = (pibar**aw)*(dpm1**(1.0d0-aw))
gzwage = gz*techshk**(1.0d0-bw)
dp = vp*dptildem1
dw = vw*dwtildem1*gzwage
muc = lam + (gamma/gz)*beta*bc
cc = gamma*ccm1/(gz*techshk)+1.0d0/muc
cquad_vi = bi/(qq*invshk)-(1.0d0-qq*invshk)/(phii*qq*invshk)

vi = 0.5d0*(1.0d0+sqrt(1.0d0+4.0d0*cquad_vi))


inv = vi*invm1/techshk
aayy = 1.0d0/gshk-(phi/2.0d0)*(vp-1.0d0)*(vp-1.0d0)
rkss = gz/beta-1.0d0+delta
utilcost = (rkss/sigmaa)*(exp(sigmaa*(util-1.0d0))-1.0d0)
gdp = (1.0d0/aayy)*( cc+inv + utilcost*(capm1/(gz*techshk)) )
rw = rwm1*dw/(gz*techshk*dp)
lab = gdp**(1.0d0/(1.0d0-alpha))*(util*capm1/(gz*techshk))**(alpha/(alpha-1.0d0))/ashk
xhp = alpha*log(util) + (1-alpha)*(log(lab)-llabss)


lrss = log(gz*pibar/beta)
notr = exp( lrss+gam_rs*(log(notrm1)-lrss) + (1.0d0-gam_rs)*( gam_dp*log(dp/pibar)  +&
     gam_dy*log(gdp*techshk/gdpm1) + gam_xhp*xhp ) + rrshk )


mc = rw*lab/((1.0d0-alpha)*gdp)
rentalk = (alpha/(1.0d0-alpha))*(rw*lab*gz*techshk/(util*capm1))
cap = (1.0d0-delta)*(capm1/(gz*techshk)) + invshk*inv*(1.0d0- (phii/2.0d0)*(vi-1.0d0)*(vi-1.0d0) )
nomr = notr


!all variables are in levels except for shocks
!shocks are in log-level deviations from SS
!(1) cap, (2) cc, (3) inv, (4) rw, (5) notr, (6) dp, (7) gdp, (8) xhp, (9) nomr, (10) lam, (11) qq, 
!(12) lab, (13) util, (14) mc, (15) rentalk, (16) muc, (17) vi, (18) vp, (19) vw, (20) dw, (21) bc, (22) bi,
!(23) liqshk, (24) invshk, (25) techshk, (26) rrshk, (27) gshk, (28) ashk
endogvar = (/cap, cc, inv, rw, notr, dp, gdp, xhp, nomr, lam, &
     qq, lab, util, mc, rentalk, muc, vi, vp, vw, dw, bc, bi, currentshockvalues/) 
     
end subroutine intermediatedec


!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: decr_euler
!> @author Christopher Gust
!> @date 6-28-16

!> @brief For a given collocation point, return associated errors.
!> @param[in] gridindex Index of one of the collocation points.  
!> @param[in] shockpos Index of current exogenous state. 
!> @param[in] params Model parameters  
!> @param[in] poly Polynomial details.  
!> @param[in] alphacoeff Polynomial coefficients.
!> @param[in] slopeconxx Slope coefficients and constants for converting xx to msv space. 
!> @param[in] bbt Matrix with Smolyak basis functions along its columns evaluated at each grid point (rows).  
!> @param[in] xgrid Matrix of grid points along its columns. Each grid point is NmsvX1 
!> @param[in] zlbinfo Indicator equal to 0 if using zlb polynomials; 1 for normal polynomials
!> @param[out] polyappnew Updated value of polynomials at collocation point.
!> @param[out] errsum Sum of squared errors of difference between current and new polynomials.
!> @param[out] errmax Max absolute difference between current and new polynomials.
!----------------------------------------------------------------------------------------------------------
subroutine decr_euler(gridindex,shockpos,params,poly,alphacoeff,slopeconxx,&
     bbt,xgrid,zlbinfo,polyappnew,errsum,errmax)   

implicit none
type(polydetails), intent(in)   :: poly
integer,          intent(in)    :: shockpos
integer,          intent(in)    :: gridindex
double precision, intent(in)    :: params(poly%nparams)
double precision, intent(in)    :: bbt(poly%ngrid,poly%ngrid)
double precision, intent(in)    :: xgrid(poly%nmsv+poly%nexogcont,poly%ngrid)
double precision, intent(in)    :: alphacoeff(poly%nfunc*poly%ngrid,2*poly%ns)
double precision, intent(in)    :: slopeconxx(2*(poly%nmsv+poly%nexogcont))
integer,          intent(in)    :: zlbinfo
double precision, intent(out)   :: polyappnew(2*poly%nfunc)  
double precision, intent(out)   :: errsum
double precision, intent(out)   :: errmax  

integer          :: ifunc,ss,shockpospoly,imax,nmsvplus
double precision :: currentshockvalues(poly%nexog),polyapp(2*poly%nfunc)
double precision :: innovations(poly%nexog)
double precision :: endogvarm1(poly%nvars+poly%nexog)
double precision :: endogvar(poly%nvars+poly%nexog),endogvarzlb(poly%nvars+poly%nexog)
double precision :: endogvarp(poly%nvars+poly%nexog),endogvarzlbp(poly%nvars+poly%nexog)
double precision :: slopeconcont(2*poly%nexogcont),xgridshock(poly%nexogcont)
double precision :: slopeconxxmsv(2*poly%nmsv),xgridmsv(poly%nmsv)
double precision :: ev(12),exp_eul(12),exp_var(12)
double precision :: abserror(2*poly%nfunc)
double precision :: rkss,utilcostp,techshkp,invshkp,invshk,liqshk,ep
double precision :: llabss
double precision :: omegapoly
logical          :: zlbintermediate

nmsvplus = poly%nmsv+poly%nexogcont
currentshockvalues = 0.0d0
if (poly%nexogcont > 0) then
   slopeconcont(1:poly%nexogcont) = slopeconxx(poly%nmsv+1:nmsvplus)
   slopeconcont(poly%nexogcont+1:2*poly%nexogcont) = slopeconxx(nmsvplus+poly%nmsv+1:2*nmsvplus)
   xgridshock = xgrid(poly%nmsv+1:nmsvplus,gridindex)
   currentshockvalues(poly%nexog-poly%nexogcont+1:poly%nexog) = msv2xx(xgridshock,poly%nexogcont,slopeconcont)
end if

currentshockvalues(1:poly%nexogshock) = poly%exoggrid(1:poly%nexogshock,shockpos)
shockpospoly = shockpos + poly%ns
polyapp = 0.0d0
do ifunc = 1,poly%nfunc
  polyapp(ifunc) = dot_product(alphacoeff((ifunc-1)*poly%ngrid+1:ifunc*poly%ngrid,&
       shockpos),bbt(:,gridindex))
  polyapp(poly%nfunc+ifunc) = dot_product(alphacoeff((ifunc-1)*poly%ngrid+1:ifunc*poly%ngrid,&
       shockpospoly),bbt(:,gridindex))
end do
!print*,shockpos
!print*,'polyapp'
! do ifunc = 1,poly%ngrid
!    print*,alphacoeff(ifunc,1)
! end do
! do ifunc = 1,poly%ngrid
!    print*,bbt(ifunc,gridindex)
! end do
!print*,'polyapp',polyapp(1)

endogvarm1 = 0.0d0
xgridmsv = xgrid(1:poly%nmsv,gridindex)
slopeconxxmsv(1:poly%nmsv) = slopeconxx(1:poly%nmsv)
slopeconxxmsv(poly%nmsv+1:2*poly%nmsv) = slopeconxx(nmsvplus+1:nmsvplus+poly%nmsv)
endogvarm1(1:poly%nmsv) = exp( msv2xx(xgridmsv,poly%nmsv,slopeconxxmsv) )

endogvarm1(1:poly%nmsv) = exp( msv2xx(xgridmsv,poly%nmsv,slopeconxxmsv) )
zlbintermediate = .false.  !force evaluation of 1 poly case at date t
llabss = poly%endogsteady(12)

call intermediatedec(poly%nparams,poly%nvars,poly%nexog,poly%nfunc,endogvarm1,&
     currentshockvalues,params,polyapp(1:poly%nfunc),llabss,endogvar,omegapoly,polyapp(1:poly%nfunc),zlbintermediate)
if ((zlbinfo .ne. 0) .and. (poly%zlbswitch .eqv. .true.)) then 
   call intermediatedec(poly%nparams,poly%nvars,poly%nexog,poly%nfunc,endogvarm1,&
     currentshockvalues,params,polyapp(poly%nfunc+1:2*poly%nfunc),llabss,endogvarzlb,omegapoly,&
     polyapp(poly%nfunc+1:2*poly%nfunc),zlbintermediate)
   endogvarzlb(9) = 1.0d0
end if
rkss = params(3)/params(1)-1.0d0+params(16)


exp_var = 0.0d0
innovations = 0.0d0

do ss = 1,poly%nquad
   innovations(1:poly%nexogshock) = poly%ghnodes(:,ss) 
   if (poly%nexogcont > 0) innovations(poly%nexog-poly%nexogcont+1:poly%nexog) =  poly%ghnodes(poly%nexogshock+1:poly%nexogshock+poly%nexogcont,ss) 
   call decr(endogvar,innovations,params,poly,alphacoeff,endogvarp)
   !if((shockpos==1).and.(ss==1)) print*,'endogvarp',endogvarp(1)

   techshkp = exp(endogvarp(25))
   invshkp = exp(endogvarp(24))
   ev(1) = endogvarp(10)/(endogvarp(6)*techshkp)


   utilcostp = (rkss/params(18))*(exp(params(18)*(endogvarp(13)-1.0d0))-1.0d0)
   ev(2) = (endogvarp(10)/techshkp)*( (endogvarp(15)*endogvarp(13))-utilcostp+(1.0d0-params(16))*endogvarp(11) )
   ev(3) = endogvarp(10)*(endogvarp(18)-1.0d0)*endogvarp(18)*endogvarp(7)
   ev(4) = (endogvarp(19)-1.0d0)*endogvarp(19) 
   ev(5) = endogvarp(16)/techshkp 
   ev(6) = endogvarp(10)*endogvarp(11)*invshkp*(endogvarp(17)-1.0d0)*endogvarp(17)*endogvarp(17)

   if ((zlbinfo .ne. 0) .and. (poly%zlbswitch .eqv. .true.)) then
      call decr(endogvarzlb,innovations,params,poly,alphacoeff,endogvarzlbp)
      ev(7) = endogvarzlbp(10)/(endogvarzlbp(6)*techshkp) 
      utilcostp = (rkss/params(18))*(exp(params(18)*(endogvarzlbp(13)-1.0d0))-1.0d0)
      ev(8) = (endogvarzlbp(10)/techshkp)*( (endogvarzlbp(15)*endogvarzlbp(13))-utilcostp+(1-params(16))*&
           endogvarzlbp(11) )
      ev(9) = endogvarzlbp(10)*(endogvarzlbp(18)-1.0d0)*endogvarzlbp(18)*endogvarzlbp(7)
      ev(10) = (endogvarzlbp(19)-1.0d0)*endogvarzlbp(19) 
      ev(11) = endogvarzlbp(16)/techshkp 
      ev(12) = endogvarzlbp(10)*endogvarzlbp(11)*invshkp*(endogvarzlbp(17)-1.0d0)*endogvarzlbp(17)*endogvarzlbp(17)
   else
      ev(7:12) = ev(1:6)
   end if
   exp_var = exp_var + poly%ghweights(ss)*ev   
end do

liqshk = exp(endogvar(23))
invshk = exp(endogvar(24))
ep = params(9)
exp_eul(1) =  (params(1)/params(3))*liqshk*endogvar(9)*exp_var(1)
exp_eul(2) = (params(1)/params(3))*exp_var(2)/endogvar(10)
exp_eul(3) =  params(1)*exp_var(3)/(endogvar(10)*endogvar(7))+&
     (ep/params(7))*( endogvar(14)-(ep-1.0d0)/ep ) 
exp_eul(4) = params(1)*exp_var(4)+(params(10)/params(8))*endogvar(10)*endogvar(12)*&
           ( (params(4)*endogvar(12)**params(6)/endogvar(10))-((params(10)-1.0d0)/&
           params(10))*endogvar(4) ) 
exp_eul(5) = exp_var(5)
exp_eul(6) = (params(1)/params(3))*exp_var(6)/endogvar(10) - 0.5d0*endogvar(11)*invshk*(endogvar(17)-1.0d0)*(endogvar(17)-1.0d0)

polyappnew(1) = log(exp_eul(1))
polyappnew(2) = log(exp_eul(2))
polyappnew(3) = exp_eul(3)
polyappnew(4) = exp_eul(4)
polyappnew(5) = log(exp_eul(5))
polyappnew(6) = exp_eul(6)
polyappnew(7) = log( 1.0d0 + (1/params(18))*log(endogvar(15)/rkss) )

if ((zlbinfo .ne. 0) .and. (poly%zlbswitch .eqv. .true.)) then
   exp_eul(7) =  (params(1)/params(3))*liqshk*endogvarzlb(9)*exp_var(7)
   exp_eul(8) = (params(1)/params(3))*exp_var(8)/endogvarzlb(10)
   exp_eul(9) =  params(1)*exp_var(9)/(endogvarzlb(10)*endogvarzlb(7))+&
     (ep/params(7))*( endogvarzlb(14)-(ep-1.0d0)/ep )
   exp_eul(10) = params(1)*exp_var(10)+(params(10)/params(8))*endogvarzlb(10)*endogvarzlb(12)*&
           ( (params(4)*endogvarzlb(12)**params(6)/endogvarzlb(10))-((params(10)-1.0d0)/&
           params(10))*endogvarzlb(4) ) 
   exp_eul(11) = exp_var(11)
   exp_eul(12) = (params(1)/params(3))*exp_var(12)/endogvarzlb(10) - 0.5d0*endogvarzlb(11)*invshk*&
        (endogvarzlb(17)-1.0d0)*(endogvarzlb(17)-1.0d0)

   polyappnew(8) = log(exp_eul(7))
   polyappnew(9) = log(exp_eul(8))
   polyappnew(10) = exp_eul(9) 
   polyappnew(11) = exp_eul(10) 
   polyappnew(12) = log(exp_eul(11))
   polyappnew(13) = exp_eul(12)
   polyappnew(14) = log( 1.0d0 + (1/params(18))*log(endogvarzlb(15)/rkss) )
else
   exp_eul(7:12) = exp_eul(1:6)
   polyappnew(8:14) = polyappnew(1:7)
end if

errsum = 0.0d0
errmax = 0.0d0
imax = 0
do ifunc = 1,2*poly%nfunc
   abserror(ifunc) = abs(polyappnew(ifunc)-polyapp(ifunc))
   errsum = errsum + abserror(ifunc)
   if (abserror(ifunc) > errmax) then
      errmax = abserror(ifunc)
      imax = ifunc
   end if
end do
errsum = errsum/dble(2*poly%nfunc)

end subroutine decr_euler

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: get_shockdetails
!> @author Christopher Gust
!> @date 6-28-16

!> @brief Get the markov processes for the shocks.
!> @param[in] nparams Number of parameter values.
!> @param[in] nexog Number of exogenous shock processes.
!> @param[in] nexogshock Number of exogenous shock processes in finite-element part of approximation.
!> @param[in] nexogcont number of exogenous shocks used in poly approximation.
!> @param[in] ns Number of exogenous states.
!> @param[in] volatilityindex Indicator for shock that allows for two volatility regimes.
!> @param[in] number_shock_values Number of shock values.
!> @param[in] number_shock_values_sq Number of shock values squared.
!> @param[in] ngridshocks Array indexing number of realizations for each shock (length nexog). 
!> @param[in] params Model parameters.
!> @param[out] shockinfo Matrix containing the values of the shocks for each exogenous state.  
!> @param[out] pibig Vector of transition probabilities.
!> @param[out] trans_update Vector of cumulative transition probabilities for each shock.
!----------------------------------------------------------------------------------------------------------
subroutine get_shockdetails(nparams,nexog,nexogshock,nexogcont,ns,number_shock_values,ngridshocks,params,&
     exogvarinfo,exoggrid,shockbounds,shockdistance)

implicit none
integer, intent(in) :: nparams,nexog,nexogshock,nexogcont,ns
integer, intent(in) :: number_shock_values
integer, intent(in) :: ngridshocks(nexog-nexogcont)
integer, intent(in) :: exogvarinfo(nexog-nexogcont,ns)
double precision, intent(in)  :: params(nparams)
double precision, intent(out) :: exoggrid(nexog-nexogcont,ns)
double precision, intent(out) :: shockbounds(nexogshock,2),shockdistance(nexogshock)

integer :: i,nshocksum,ss,nall,shockpos
integer :: currentshockindex(nexogshock),nshocksum_vec(nexogshock)
double precision :: shockvalues(number_shock_values)
double precision :: rhog,rhoinv,rholiq,rhoint,rhoa
double precision :: sdevtech,sdevg,sdevinv,sdevliq,sdevint,sdeva
double precision :: rhovec(nexog),sigmavec(nexog)
double precision, allocatable :: xshock(:)

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

rhovec = (/rholiq,rhoinv,0.0d0,rhoint,rhog,rhoa/)
sigmavec = (/sdevliq,sdevinv,sdevtech,sdevint,sdevg,sdeva/)

nshocksum = 0
shockvalues = 0.0d0
do i = 1,nexogshock
   allocate(xshock(ngridshocks(i)))
   call finite_grid(ngridshocks(i),rhovec(i),sigmavec(i),shockdistance(i),xshock)
   shockbounds(i,1) = xshock(1)
   shockbounds(i,2) = xshock(ngridshocks(i))
   shockvalues(nshocksum+1:nshocksum+ngridshocks(i)) = xshock
   nshocksum = nshocksum + ngridshocks(i)
   deallocate(xshock)
end do

exoggrid = 0.0d0
do ss = 1,ns
   currentshockindex = exogvarinfo(1:nexogshock,ss)
   nshocksum_vec = 0
   nall = ngridshocks(1)
   shockpos = currentshockindex(1)
   exoggrid(1,ss) = shockvalues(currentshockindex(1))
   do i = 2,nexogshock
      nshocksum_vec(i) = nshocksum_vec(i-1) + ngridshocks(i-1)
      shockpos = shockpos + nall*(currentshockindex(i)-1)
      nall = nall*ngridshocks(i)
      exoggrid(i,ss) = shockvalues(nshocksum_vec(i) + currentshockindex(i))
   end do
end do

end subroutine get_shockdetails

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: calc_premium 
!> @author Christopher Gust
!> @date 12-18-15

!> @brief Calculates the capital financial premium in the model.  
!> @param[in] endogvar Endogenous variables and shocks (lagged) necessary for the particle filter.
!> @param[in] innovations Model innovations. 
!> @param[in] params Model parameters.  
!> @param[in] poly Polynomial details.
!> @param[in] alphacoeff Polynomial coefficients.
!> @param[out] premium Capital Financial Premium.
!----------------------------------------------------------------------------------------------------------
subroutine calc_premium(endogvar,params,poly,alphacoeff,premium)

implicit none
type(polydetails), intent(in)   :: poly
double precision, intent(in)    :: params(poly%nparams)
double precision, intent(in)    :: endogvar(poly%nvars+poly%nexog)
double precision, intent(in)    :: alphacoeff(poly%nfunc*poly%ngrid,2*poly%ns)
double precision, intent(out)   :: premium 

integer :: ss
double precision :: endogvarp(poly%nvars+poly%nexog)
double precision :: innovations(poly%nexog)
double precision :: techshkp,invshkp,ev,utilcostp,exp_var,rkss

rkss = params(3)/params(1)-1.0d0+params(16)

innovations = 0.0d0
exp_var = 0.0d0
do ss = 1,poly%nquad
   innovations(1:poly%nexogshock) = poly%ghnodes(:,ss) 
   call decr(endogvar,innovations,params,poly,alphacoeff,endogvarp)
   techshkp = exp(endogvarp(25))
   invshkp = exp(endogvarp(24))
   
   utilcostp = (rkss/params(18))*(exp(params(18)*(endogvarp(13)-1.0d0))-1.0d0)
   ev = ( (endogvarp(15)*endogvarp(13))-utilcostp+(1.0d0-params(16))*endogvarp(11) )*endogvarp(6)
   exp_var = exp_var + poly%ghweights(ss)*ev   
end do

premium = exp_var/(endogvar(11)*endogvar(9))

end subroutine calc_premium

end module model_details





