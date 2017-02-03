module class_Model

  use gensys
  use filter
  use prior, only: priordens

  implicit none

  !double precision, parameter :: M_PI = 3.14159265358979323846d0

  type Model 

     character(len=144) :: name = 'model_v5_2_linear'

     character(len=300) :: datafile = 'yy.txt'
     character(len=300) :: priorfile = 'prior.txt'
     character(len=300) :: transfile = 'trans.txt'

     double precision, allocatable :: yy(:,:) 
     integer :: nobs = 5, T = 125, ny = 5

     integer :: nexog = 8, neps = 8
     integer :: nvars = 42, neta = 8, neq = 42
     
     integer :: npara = 28
     double precision, allocatable :: params(:)
     double precision, allocatable :: TT(:, :), RR(:, :), QQ(:, :), DD(:), ZZ(:, :), HH(:,:)
     double precision, allocatable  :: GAM0(:, :), GAM1(:, :), C(:), PSI(:, :), PPI(:, :), CC(:)     
     double precision, allocatable :: trspec(:,:)
     double precision, allocatable :: pmean(:), pstdd(:), pfix(:)
     integer, allocatable :: pshape(:), pmask(:)

   contains
     procedure :: load_data
     procedure :: steadystate
     generic :: solve => solve_serial, solve_parallel
     procedure :: solve_serial
     procedure :: solve_parallel

     procedure :: g 
     procedure :: pdfy
     procedure :: logpdfy_kernel
     procedure :: describe
     procedure :: lik
     procedure :: pr
     procedure :: inbounds
     procedure :: lik_states
  end type Model

  interface Model
     module procedure new_Model
  end interface Model


contains

  double precision function lik(m, p0)

    class(Model) :: m
    double precision :: p0(m%npara)

    logical :: converged 
    converged = m%solve(p0, 0, 0)

    if (converged==.false.) then
       lik = -1000000000000.0d0
       return 
    end if

    lik = chand_recursion(m%yy, m%TT, m%RR, m%QQ, m%DD, m%ZZ, m%HH, m%ny, m%T, m%neps, m%neq, 0)

  end function lik

  subroutine  lik_states(m, p0, lik0, filter_states,filter_vars,filter_shocks, smooth_states,smooth_vars,smooth_shocks)

    class(Model), intent(inout) :: m
    double precision :: p0(m%npara)
    real(wp), intent(out) :: lik0, smooth_states(m%T,m%nvars), smooth_shocks(m%T,m%nexog), smooth_vars(m%nvars,m%nvars,m%T)
    real(wp), intent(out) :: filter_states(m%T,m%nvars), filter_shocks(m%T,m%nexog), filter_vars(m%nvars,m%nvars,m%T)

    logical :: converged 
    converged = m%solve(p0, 0, 0)

    if (converged==.false.) then
       lik0 = -1000000000000.0d0
       return 
    end if

    call kalman_filter_missing_with_states(m%yy, m%TT, m%RR, m%QQ, m%DD, m%ZZ, m%HH, m%ny, m%T, m%neps, m%neq, 0, lik0, filter_states, filter_vars, filter_shocks, smooth_states, smooth_vars, smooth_shocks)

  end subroutine lik_states

  double precision function pr(m, p0)

    class(Model) :: m
    real(wp), intent(inout) :: p0(m%npara)


    ! check whether the parameters are inbounds 
    if (m%inbounds(p0) .ne. .true.) then
       pr = -1000000000.d0 !REALLY_NEG
       return
    end if

    pr = priordens(p0,m%pmean,m%pstdd,m%pshape,m%pmask,m%pfix)

  end function pr

  logical function inbounds(m, p0)

    class(Model) :: m
    real(wp), intent(inout) :: p0(m%npara)


    real(wp) :: low, hi
    integer :: i

    do i = 1, m%npara
       low = m%trspec(2,i)
       hi = m%trspec(3,i)

       if (m%pmask(i)==1) then
          p0(i) = m%pfix(i)
       end if

       if ((p0(i) < low .or. p0(i) > hi) .and. (m%pmask(i)==0)) then
          !print*,'parameter ', i, 'out of bounds', p0(i), low, hi
          inbounds = .false.
          return
       end if


    end do

    inbounds = .true.

  end function inbounds

  subroutine load_data(m)
    class(Model) :: m
    integer :: i

    open(1, file=trim(m%datafile), action='read')

    do i = 1, m%T
       read(1,*) m%yy(:,i)
    end do

    close(1)

    open(1, file=m%priorfile, status='old', action='read')
    do i = 1, m%npara
       read(1, *) m%pshape(i), m%pmean(i), m%pstdd(i), m%pmask(i), m%pfix(i)
    end do
    close(1)


    open(1, file=m%transfile, status='old', action='read')
    do i = 1, m%npara
       read(1, *) m%trspec(:,i)
    end do
    close(1)



  end subroutine load_data


  double precision function steadystate(m, params) 

    class(Model) :: m
    double precision, intent(in)  :: params(m%npara) 

    steadystate = 0.0d0
  end function steadystate

  type(Model) function new_Model(blar) result(m)

    logical, intent(in), optional :: blar ! to have the same signature as nonlinear model

    allocate(m%params(m%npara))
    m%params = (/0.130169219986d0,100d0*(1.0025d0 - 1d0),exp(1.00586879d0/100d0),0.0d0,0.52323556d0,2.04724047d0,113.36316712d0,4029.72644029d0,1-0.26556192d0,1-0.68732933d0,0.32562593d0,3.48130086d0,5.64497737d0,0.69449115d0,1.9106d0,0.27180239d0,0.00420848d0*100d0,0.4374224d0,0.00177246d0*100d0,0.4415008d0,0.04946043d0*100d0,0.00137775d0*100d0,0.00098857d0,0.5d0,0.0005d0,0.5d0,0.0005d0,0.0005d0/)


    allocate(m%TT(m%neq, m%neq), m%RR(m%neq, m%neps), m%QQ(m%neps, m%neps), m%DD(m%ny), m%ZZ(m%ny, m%neq), m%HH(m%ny,m%ny))
    allocate(m%GAM0(m%neq, m%neq), m%GAM1(m%neq, m%neq), m%C(m%neq), m%PSI(m%neq, m%neps), m%PPI(m%neq, m%neta), m%CC(m%neq))




    allocate(m%yy(m%nobs, m%T))

    allocate(m%trspec(4,m%npara), m%pshape(m%npara), m%pmean(m%npara), m%pstdd(m%npara), m%pmask(m%npara), m%pfix(m%npara))



    call m%load_data

  end function new_Model

  function solve_serial(m, para) result(convergence)

    class(Model) :: m

    double precision, intent(in) :: para(m%npara)

    logical :: convergence

    convergence = m%solve_parallel(para,1,0)

  end function solve_serial

  function solve_parallel(m, para, nproc, rank) result(convergence)
    class(Model) :: m

    double precision, intent(in) :: para(m%npara)
    integer, intent(in) :: nproc, rank

    logical :: convergence

    double precision :: beta_tr,pibar_tr,gz_tr,gamxhp,gamma,sigmal,phi,phiw,ap_tr,aw_tr,alpha,phii,sigmaa,gam_rs,gam_dp,gamdy,sdevtech_tr,rhog,sdevg_tr,rhoinv,sdevinv_tr,sdevliq_tr,sdevint_tr,sdeva_tr
    double precision :: rbar,gz,rholiq,ep_tr,epw_tr,shrgy,delta,lamhp,psil,beta,pibar,ep,epw,ap,aw,bw,sdevtech,sdevg,sdevinv,sdevliq,sdevint,sdeva,gg,gamtil,mc,k2yrat,shriy,shrcy,labss,kappaw,kappap,kss,gdpss,invss,phii_jpt,css,rwss,mucss,lamss,rss,rkss,rhoint,rhoa

    double precision :: rhoelast,sdevelast_tr,rhoelastw,sdevelastw_tr,sdevelastw,sdevelast

    ! gensys
    double precision :: fmat, fwt, ywt, gev, loose, DIV
    integer :: eu(2)

    beta_tr = para(1)
    pibar_tr = para(2)
    gz_tr = para(3)
    gamxhp = para(4)
    gamma = para(5)
    sigmal = para(6)
    phi = para(7)
    phiw = para(8)
    ap_tr = para(9)
    aw_tr = para(10)
    alpha = para(11)
    phii = para(12)
    sigmaa = para(13)
    gam_rs = para(14)
    gam_dp = para(15)
    gamdy = para(16)
    sdevtech_tr = para(17)
    rhog = para(18)
    sdevg_tr = para(19)
    rhoinv = para(20)
    sdevinv_tr = para(21)
    sdevliq_tr = para(22)
    sdevint_tr = para(23)
    rhoelast = para(24)
    sdevelast_tr = para(25)
    rhoelastw = para(26)
    sdevelastw_tr = para(27)
    sdeva_tr = para(28)
    gz=exp(gz_tr/100.0)
    rholiq=0.85
    ep_tr=0.2
    epw_tr=0.2
    shrgy=0.2
    delta=0.025
    lamhp=1600.0
    psil=1.0
    beta=1.0/(beta_tr/100.0 + 1.0)
    pibar=1.0/100.0*pibar_tr + 1
    ep=(1.0+ep_tr)/ep_tr
    epw=(1.0+epw_tr)/epw_tr
    ap=1-ap_tr
    aw=1-aw_tr
    bw=aw
    sdevtech=sdevtech_tr / 100.0
    sdevg=sdevg_tr / 100.0
    sdevinv=sdevinv_tr / 100.0
    sdevliq=sdevliq_tr / 100.0
    sdevint=sdevint_tr / 100.0
    sdevelast=sdevelast_tr / 100.0
    sdevelastw=sdevelastw_tr / 100.0
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
    rhoint=0.0
    sdeva=sdeva_tr / 100.0



    m%GAM0 = 0.0d0
    m%GAM1 = 0.0d0
    m%PSI = 0.0d0
    m%PPI = 0.0d0
    m%C = 0.0d0

    m%GAM0(3, 1) = 1.0d0
    m%GAM0(2, 2) = -1.0d0*beta*gamtil**2.0d0 - 1.0d0
    m%GAM0(7, 2) = gg*shrcy
    m%GAM0(16, 2) = -1.0d0
    m%GAM0(21, 2) = gamtil*1.0/(-1.0d0*gamtil + 1.0d0)
    m%GAM0(37, 2) = -1.0d0
    m%GAM0(3, 3) = 1.0/gz*(-1.0d0*delta + 1.0d0) - 1.0d0
    m%GAM0(7, 3) = gg*shriy
    m%GAM0(11, 3) = phii_jpt*(beta + 1.0d0)
    m%GAM0(17, 3) = 1.0d0
    m%GAM0(22, 3) = -1.0d0
    m%GAM0(35, 3) = -1.0d0
    m%GAM0(4, 4) = -1.0d0
    m%GAM0(14, 4) = 1.0d0
    m%GAM0(15, 4) = 1.0d0
    m%GAM0(19, 4) = -1.0d0*kappaw
    m%GAM0(5, 5) = -1.0d0
    m%GAM0(9, 5) = 1.0d0
    m%GAM0(10, 5) = 1.0d0
    m%GAM0(4, 6) = -1.0d0
    m%GAM0(5, 6) = gam_dp*(-1.0d0*gam_rs + 1.0d0)
    m%GAM0(6, 6) = -1.0d0
    m%GAM0(18, 6) = 1.0d0
    m%GAM0(40, 6) = -1.0d0
    m%GAM0(5, 7) = gamdy*(-1.0d0*gam_rs + 1.0d0)
    m%GAM0(7, 7) = -1.0d0
    m%GAM0(12, 7) = 1.0d0
    m%GAM0(14, 7) = -1.0d0
    m%GAM0(5, 8) = gamxhp*(-1.0d0*gam_rs + 1.0d0)
    m%GAM0(8, 8) = -1.0d0
    m%GAM0(9, 9) = -1.0d0
    m%GAM0(1, 10) = -1.0d0
    m%GAM0(2, 10) = -1.0d0*(-1.0d0*gamtil + 1.0d0)*(-1.0d0*beta*gamtil + 1.0d0 &
         )
    m%GAM0(10, 10) = -1.0d0
    m%GAM0(19, 10) = -1.0d0*kappaw
    m%GAM0(42, 10) = -1.0d0
    m%GAM0(1, 11) = -1.0d0
    m%GAM0(11, 11) = -1.0d0
    m%GAM0(41, 11) = -1.0d0
    m%GAM0(8, 12) = -1.0d0*alpha + 1.0d0
    m%GAM0(12, 12) = alpha - 1.0d0
    m%GAM0(14, 12) = 1.0d0
    m%GAM0(15, 12) = 1.0d0
    m%GAM0(19, 12) = kappaw*sigmal
    m%GAM0(7, 13) = alpha*1.0/ep*gg*(ep - 1.0d0)
    m%GAM0(8, 13) = alpha
    m%GAM0(12, 13) = -1.0d0*alpha
    m%GAM0(13, 13) = -1.0d0
    m%GAM0(15, 13) = -1.0d0
    m%GAM0(6, 14) = kappap
    m%GAM0(14, 14) = -1.0d0
    m%GAM0(13, 15) = 1.0/sigmaa
    m%GAM0(15, 15) = -1.0d0
    m%GAM0(36, 15) = -1.0d0
    m%GAM0(16, 16) = gamtil - 1.0d0
    m%GAM0(17, 17) = -1.0d0
    m%GAM0(18, 18) = -1.0d0
    m%GAM0(19, 19) = -1.0d0
    m%GAM0(20, 19) = 1.0d0
    m%GAM0(38, 19) = -1.0d0
    m%GAM0(4, 20) = 1.0d0
    m%GAM0(20, 20) = -1.0d0
    m%GAM0(21, 21) = -1.0d0
    m%GAM0(22, 22) = -1.0d0
    m%GAM0(10, 23) = 1.0d0
    m%GAM0(27, 23) = -1.0d0
    m%GAM0(3, 24) = 1.0/gz*(-1.0d0*delta + 1.0d0) - 1.0d0
    m%GAM0(11, 24) = -1.0d0
    m%GAM0(28, 24) = -1.0d0
    m%GAM0(2, 25) = -1.0d0*gamtil
    m%GAM0(3, 25) = 1.0/gz*(-1.0d0*delta + 1.0d0)
    m%GAM0(4, 25) = -1.0d0
    m%GAM0(5, 25) = gamdy*(-1.0d0*gam_rs + 1.0d0)
    m%GAM0(11, 25) = phii_jpt
    m%GAM0(12, 25) = alpha
    m%GAM0(15, 25) = 1.0d0
    m%GAM0(16, 25) = -1.0d0*gamtil
    m%GAM0(17, 25) = 1.0d0
    m%GAM0(20, 25) = -1.0d0*bw + 1.0d0
    m%GAM0(29, 25) = -1.0d0
    m%GAM0(39, 25) = -1.0d0
    m%GAM0(5, 26) = 1.0d0
    m%GAM0(30, 26) = -1.0d0
    m%GAM0(7, 27) = 1.0d0
    m%GAM0(31, 27) = -1.0d0
    m%GAM0(32, 28) = -1.0d0
    m%GAM0(33, 29) = -1.0d0
    m%GAM0(12, 30) = alpha - 1.0d0
    m%GAM0(34, 30) = -1.0d0
    m%GAM0(23, 31) = -1.0d0
    m%GAM0(24, 32) = -1.0d0
    m%GAM0(25, 33) = -1.0d0
    m%GAM0(26, 34) = -1.0d0
    m%GAM0(11, 35) = -1.0d0*beta*phii_jpt
    m%GAM0(22, 35) = 1.0d0
    m%GAM0(1, 36) = -1.0d0*beta*1.0/gz*(-1.0d0*delta + 1.0d0) + 1.0d0
    m%GAM0(2, 37) = beta*gamtil
    m%GAM0(21, 37) = -1.0d0*1.0/(-1.0d0*gamtil + 1.0d0)
    m%GAM0(19, 38) = beta
    m%GAM0(1, 39) = -1.0d0
    m%GAM0(2, 39) = beta*gamtil
    m%GAM0(10, 39) = -1.0d0
    m%GAM0(11, 39) = -1.0d0*beta*phii_jpt
    m%GAM0(21, 39) = -1.0d0*gamtil*1.0/(-1.0d0*gamtil + 1.0d0)
    m%GAM0(22, 39) = 1.0d0
    m%GAM0(6, 40) = beta*1.0/(beta*(-1.0d0*ap + 1.0d0) + 1.0d0)
    m%GAM0(10, 40) = -1.0d0
    m%GAM0(1, 41) = beta*1.0/gz*(-1.0d0*delta + 1.0d0)
    m%GAM0(1, 42) = 1.0d0
    m%GAM0(10, 42) = 1.0d0

    m%GAM1(3, 1) = 1.0/gz*(-1.0d0*delta + 1.0d0)
    m%GAM1(12, 1) = alpha
    m%GAM1(15, 1) = 1.0d0
    m%GAM1(2, 2) = -1.0d0*gamtil
    m%GAM1(16, 2) = -1.0d0*gamtil
    m%GAM1(24, 2) = -1.0d0
    m%GAM1(11, 3) = phii_jpt
    m%GAM1(17, 3) = 1.0d0
    m%GAM1(25, 3) = -1.0d0
    m%GAM1(4, 4) = -1.0d0
    m%GAM1(26, 4) = -1.0d0
    m%GAM1(5, 5) = -1.0d0*gam_rs
    m%GAM1(6, 6) = -1.0d0*(-1.0d0*ap + 1.0d0)*1.0/(beta*(-1.0d0*ap + 1.0d0) + &
         1.0d0)
    m%GAM1(18, 6) = -1.0d0*ap + 1.0d0
    m%GAM1(20, 6) = aw - 1.0d0
    m%GAM1(5, 7) = gamdy*(-1.0d0*gam_rs + 1.0d0)
    m%GAM1(23, 7) = -1.0d0
    m%GAM1(27, 23) = -0.85d0
    m%GAM1(28, 24) = -1.0d0*rhoinv
    m%GAM1(30, 26) = -1.0d0*rhoint
    m%GAM1(31, 27) = -1.0d0*rhog
    m%GAM1(32, 28) = -1.0d0*rhoelast
    m%GAM1(33, 29) = -1.0d0*rhoelastw
    m%GAM1(35, 35) = -1.0d0
    m%GAM1(36, 36) = -1.0d0
    m%GAM1(37, 37) = -1.0d0
    m%GAM1(38, 38) = -1.0d0
    m%GAM1(39, 39) = -1.0d0
    m%GAM1(40, 40) = -1.0d0
    m%GAM1(41, 41) = -1.0d0
    m%GAM1(42, 42) = -1.0d0

    m%PSI(27, 1) = -1.0d0
    m%PSI(28, 2) = -1.0d0
    m%PSI(29, 3) = -1.0d0
    m%PSI(30, 4) = -1.0d0
    m%PSI(31, 5) = -1.0d0
    m%PSI(32, 6) = -1.0d0
    m%PSI(33, 7) = -1.0d0
    m%PSI(34, 8) = -1.0d0

    m%PPI(35, 1) = -1.0d0
    m%PPI(36, 2) = -1.0d0
    m%PPI(37, 3) = -1.0d0
    m%PPI(38, 4) = -1.0d0
m%PPI(39, 5) = -1.0d0
m%PPI(40, 6) = -1.0d0
m%PPI(41, 7) = -1.0d0
m%PPI(42, 8) = -1.0d0




call do_gensys(m%TT, m%CC, m%RR, fmat, fwt, ywt, gev, eu, loose, m%GAM0, m%GAM1, m%C, m%PSI, m%PPI, DIV)

    
    m%QQ = 0.0d0
    m%ZZ = 0.0d0
    m%DD = 0.0d0
    m%HH = 0.0d0



    !
    convergence = eu(1)*eu(2) == 1

    m%QQ(1, 1) = sdevliq**2.0d0
m%QQ(2, 2) = sdevinv**2.0d0
m%QQ(3, 3) = sdevtech**2.0d0
m%QQ(4, 4) = sdevint**2.0d0
m%QQ(5, 5) = sdevg**2.0d0
m%QQ(6, 6) = sdevelast**2.0d0
m%QQ(7, 7) = sdevelastw**2.0d0
m%QQ(8, 8) = sdeva**2.0d0

m%ZZ(2, 2) = 1.0d0
m%ZZ(3, 3) = 1.0d0
m%ZZ(4, 6) = 1.0d0
m%ZZ(1, 7) = 1.0d0
m%ZZ(5, 9) = 1.0d0
m%ZZ(1, 25) = 1.0d0
m%ZZ(2, 25) = 1.0d0
m%ZZ(3, 25) = 1.0d0
m%ZZ(1, 31) = -1.0d0
m%ZZ(2, 32) = -1.0d0
m%ZZ(3, 33) = -1.0d0


m%HH(1, 1) = 4.14900249066087d-6
m%HH(2, 2) = 2.61699653733386d-6
m%HH(3, 3) = 5.89253403646793d-5
m%HH(4, 4) = 5.94733694096693d-7
m%HH(5, 5) = 5.10160016213353d-6

m%DD(1) = log(gz)
m%DD(2) = log(gz)
m%DD(3) = log(gz)
m%DD(4) = log(pibar)
m%DD(5) = log(1.0/beta*gz*pibar)


  end function solve_parallel

  function g(m, states_old, shocks_curr, regime_curr) result(states_curr)
    class(Model) :: m

    double precision, intent(in) :: states_old(m%nvars), shocks_curr(m%nexog)
    integer, intent(in) :: regime_curr(2)
    double precision :: states_curr(m%nvars)

    double precision :: tmp(m%nexog)

    integer :: i

    states_curr = 0.0d0
    do i = 1, m%nexog
       tmp(i) = sqrt(m%QQ(i,i))*shocks_curr(i)
    end do


    call dgemv('n', m%neq, m%neps, 1.0d0, m%RR, m%neq, tmp, 1, 0.0d0, states_curr, 1)
    call dgemv('n', m%neq, m%neq, 1.0d0, m%TT, m%neq, states_old, 1, 1.0d0, states_curr, 1)

  end function g

  double precision function pdfy(m, t, states_new, states_old, para) result(pdf)
    class(Model) :: m 

    integer, intent(in) :: t

    double precision, intent(in) :: states_new(m%nvars), states_old(m%nvars)
    double precision, intent(in) :: para(m%npara)

    double precision :: z
    integer :: i

    pdf = 1.0d0
    do i = 1, m%nobs
       z = (m%yy(i,t) - m%DD(i) - dot_product(m%ZZ(i,:), states_new)) / sqrt(m%HH(i,i))
       pdf = pdf / sqrt(m%HH(i,i)) * exp(-0.5d0 * z**2)
    end do

    pdf = pdf * (1.0d0 / sqrt(2.0d0*M_PI))**m%nobs
  end function pdfy


  double precision function logpdfy_kernel(m, t, states_new, states_old, para) result(logpdf)
    class(Model) :: m 

    integer, intent(in) :: t

    double precision, intent(in) :: states_new(m%nvars), states_old(m%nvars)
    double precision, intent(in) :: para(m%npara)

    double precision :: z
    integer :: i

    logpdf = 0.0d0
    do i = 1, m%nobs
       z = (m%yy(i,t) - m%DD(i) - dot_product(m%ZZ(i,:), states_new)) / sqrt(m%HH(i,i))
       logpdf =  logpdf - 0.5d0 * z**2
    end do


  end function logpdfy_kernel



  subroutine describe(m)
    class(Model) :: m

    print*,m%name

  end subroutine describe

  subroutine describe_params(m)
    class(Model) :: m

    print*,m%name

  end subroutine describe_params

  subroutine cleanup(m)
    class(Model) :: m

    print*,m%name

  end subroutine cleanup




end module class_Model
