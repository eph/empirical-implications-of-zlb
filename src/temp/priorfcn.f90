subroutine priorfcn(nparams,paramvec,prior)
use pdfs
implicit none
integer, intent(in) :: nparams
real(8), intent(in) :: paramvec(nparams,1)
real(8), intent(out) :: prior
real(8) :: pdf
real(8) :: a1,b1
real(8) :: a2,b2
real(8) :: sf
prior = 0.0d0
pdf = 0.0d0

!structural parameters
!beta  = 0.9987d0
!pibar = 1.006d0
!gz = 1.0041d0
!psil = 1.0d0
!gamma = 0.858d0
!sigmal = 4.49d0
!phi = 95.0d0
!phiw = 8000.0d0
!ep = 6.0d0
!epw = 8.0d0
!ap = 0.87d0
!aw = 0.92d0
!bw = 0.92d0
!lamhp = 1600.0d0
!alpha = 0.167d0
!delta = 0.025d0
!phii = 0.5586d0
!sigmaa = 5.0d0
!gam_rs = 0.86d0
!gam_dp = 1.688d0
!gamxhp = 0.0d0
!gamdy = 0.21d0/(1.0d0-gam_rs)
!shrgy = 0.2d0
!sdevtech = 0.01d0
!rhog = 0.95d0
!sdevg = 0.01d0
!rhoinv = 0.77d0
!sdevinv = 0.05d0
!rholiq = 0.90d0
!sdevliq = 0.01d0
!rhoint = 0.0d0
!sdevint = 0.01d0
!rhoelast = 0.98d0
!sdevelast = 0.01d0
!rhoelastw = 0.98d0
!sdevelastw = 0.01d0
!me1 = 0.001d0
!me2 = 0.001d0
!me3 = 0.001d0
!me4 = 0.001d0
!me5 = 0.001d0
!me6 = 0.001d0
!me7 = 0.001d0

! Input prior densities here
call gamma_converter( 0.25d0, 0.1d0, a1, b1)
call gamma_pdf ( 100*(paramvec(1,1)**(-1)-1.0d0), a1, b1, pdf )    ! beta 
!print*, ' beta = ', log(pdf)
prior = prior + log(pdf)! + log(100.0d0*paramvec(1,1)**(-2))

!write(*,*) 'after 1 prior is ', prior
call normal_pdf ( 100.0d0*(paramvec(2,1)-1.0d0), 0.50d0, 0.10d0, pdf )  ! pibar
!print*, ' pibar= ', log(pdf)
prior = prior + log(pdf)! + log(100.0d0)

!write(*,*) 'after 2 prior is ', prior
call normal_pdf ( 100*log(paramvec(3,1)), 0.5d0, 0.03d0, pdf )                ! gz
!print*, ' gz= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 3 prior is ', prior
!call gamma_converter( 1.0d0, 0.10d0, a1, b1)
!call gamma_pdf ( paramvec(4,1), a1, b1, pdf )                      ! psil
!call normal_pdf ( paramvec(4,1), 1.0d0, 0.01d0, pdf )                ! psil
!call normal_pdf ( paramvec(4,1), 1.0d0, 0.1d0, pdf )                ! psil
!prior = prior + log(pdf)

!write(*,*) 'after 4 prior is ', prior
call beta_converter(0.6d0, 0.1d0, a1,b1)
call beta_pdf( paramvec(5,1), a1, b1, pdf )                         ! gamma
!print*, ' gamma= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 5 prior is ', prior
call gamma_converter( 2.0d0, 0.75d0, a1, b1)
call gamma_pdf( paramvec(6,1), a1, b1, pdf )                   ! sigmal
!print*, ' sigmal= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 6 prior is ', prior
call normal_pdf( paramvec(7,1), 100.0d0,  25.0d0, pdf )                  ! phi
!print*, ' phi= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 7 prior is ', prior
call normal_pdf( paramvec(8,1), 3000.d0, 5000.0d0, pdf )                   ! phiw
!print*, ' phiw= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 8 prior is ', prior
call normal_pdf( 1.0d0/ ( paramvec(9,1)-1.0d0) , 0.15d0,0.05d0, pdf )   ! ep
!print*, ' ep= ', log(pdf)
!prior = prior + log(pdf)! + log((paramvec(9 ,1)-1.0d0)**(-2))

! epw is fixed
!call normal_pdf( 1.0d0/ ( paramvec(10,1)-1.0d0), 0.15d0,0.05d0, pdf )   ! epw
!print*, ' epw= ', log(pdf)
!prior = prior + log(pdf)! + log((paramvec(10,1)-1.0d0)**(-2))
!write(*,*) 'after 10 prior is ', prior

call beta_converter(0.5d0, 0.15d0, a1, b1)
call beta_pdf( 1.0d0 - paramvec(11,1),   a1,   b1, pdf )         ! ap
!print*, ' ap= ', log(pdf)
prior = prior + log(pdf) 

call beta_pdf( 1.0d0 - paramvec(12,1), a1,b1, pdf )                       ! aw
!print*, ' aw= ', log(pdf)
prior = prior + log(pdf)

!bw = aw in the code - 5-21-15
!write(*,*) 'after 12 prior is ', prior
!call beta_pdf( paramvec(13,1), 2.0d0,2.0d0, pdf )  ! bw
!prior = prior + log(pdf)

!call normal_pdf( paramvec(14,1), 1600.0d0, 100.0d0, pdf )  ! lamhp
!call gamma_converter( 0.25d0, 0.1d0, a1, b1)
!call gamma_pdf ( 100.0d0*(paramvec(14,1)-1.0d0), a1, b1, pdf )    
!prior = prior + log(pdf)

!write(*,*) 'after 14 prior is ', prior
call normal_pdf( paramvec(15,1), 0.30d0,0.05d0, pdf )  ! alpha
!print*, ' alpha= ', log(pdf)
prior = prior + log(pdf)

!call normal_pdf( paramvec(16,1), 0.025d0,0.0005d0, pdf )  ! delta
!print*, ' delta= ', log(pdf)
!prior = prior + log(pdf)

!write(*,*) 'after 16 prior is ', prior

call gamma_converter( 4.0d0, 1.0d0, a1, b1)    !Updated a1,b1 -- 4-27-15
call gamma_pdf( paramvec(17,1), a1, b1, pdf )  ! phii
!print*, ' phii= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 17 prior is ', prior
call gamma_converter( 5.0d0, 1.0d0, a1, b1)   
call gamma_pdf( paramvec(18,1), a1, b1, pdf )  ! sigmaa
!print*, ' sigmaa= ', log(pdf)
prior = prior + log(pdf)


!write(*,*) 'after 18 prior is ', prior
call beta_converter(0.5d0, 0.10d0, a1, b1)
call beta_pdf( paramvec(19,1), a1, b1, pdf )      ! gam_rs
!call normal_pdf( paramvec(19,1), 0.0d0, 0.2d0, pdf )      ! gam_rs
!print*, ' gam_rs= ', log(pdf)
prior = prior + log(pdf)

call normal_pdf( paramvec(20,1), 1.70d0,0.3d0, pdf )      ! gam_dp
!print*, ' gam_dp= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 20 prior is ', prior
call normal_pdf( paramvec(21,1), 0.4d0, 0.3d0, pdf )      ! gamxhp
!!print*, ' gamxhp= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 21 prior is ', prior
call normal_pdf( paramvec(22,1), 0.4d0,0.3d0, pdf )      ! gamdy
!print*, ' gamdy= ', log(pdf)
prior = prior + log(pdf)
!write(*,*) 'after 22 prior is ', prior

!call normal_pdf( paramvec(23,1), 0.2d0, 0.005d0, pdf )      ! shrgy
!print*, ' shrgy= ', log(pdf)
!prior = prior + log(pdf)

!write(*,*) 'after 23 prior is ', prior
! This is mu = 0.1, std = 1
a1 = 0.006422351995801d0
b1 = 2.006358764352278d0

! This is mu = 0.5, std = 1
a2 = 0.19384964390576
b2 = 2.155079715124610   

!write(*,*) 'after 23 prior is ', prior

!call invgamma_pdf_1( 100.0d0 * paramvec(24,1), b2, a2, pdf )      ! sdevtech
pdf = logigpdf(100.0d0*paramvec(24,1),0.33d0,2.0d0)
prior = prior + pdf!+ log(100.0d0)
!print*, ' sdevtech= ', pdf

call beta_converter(0.6d0, 0.2d0, a1,b1)
call beta_pdf( paramvec(25,1), a1, b1, pdf )      ! rhog
!call normal_pdf( paramvec(25,1), 0.0d0, 0.1d0, pdf )      ! rhog
!print*, ' rhog= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 25 prior is ', prior
!call invgamma_pdf_1( 100.0d0 * paramvec(26,1), b2, a2, pdf )      ! sdevg
pdf = logigpdf(100.0d0*paramvec(26,1),0.33d0,2.0d0)
prior = prior + pdf! + log(100.0d0)
!print *,' sdevg= ', pdf
!write(*,*) 'after 26 prior is ', prior

!call beta_pdf( paramvec(27,1), 2.0d0, 2.0d0, pdf )      ! rhoinv
call beta_converter(0.6d0, 0.2d0, a1,b1)
call beta_pdf( paramvec(27,1), a1, b1, pdf )
!print*, ' rhoinv= ', log(pdf)
prior = prior + log(pdf)

!write(*,*) 'after 27 prior is ', prior
!call invgamma_pdf_1( 100.0d0 * paramvec(28,1), b2, a2, pdf )      ! sdevinv
pdf = logigpdf(100.0d0*paramvec(28,1),0.33d0,2.0d0)
prior = prior + pdf !+ log(100.0d0)
!print *,' sdevinv= ', pdf

!write(*,*) 'after 28 prior is ', prior
!call beta_pdf( paramvec(29,1), 2.0d0, 2.0d0, pdf )      ! rholiq
!call normal_pdf( paramvec(29,1), 0.0d0, 0.1d0, pdf )      ! rholiq
!print*, ' rholiq= ', log(pdf)
!prior = prior + log(pdf)

!call invgamma_pdf_1( 100.0d0 * paramvec(30,1), b1, a1, pdf )      ! sdevliq
pdf = logigpdf(100.0d0*paramvec(30,1),0.33d0,2.0d0)
prior = prior + pdf! + log(100.0d0)
!print *,' sdevliq= ', pdf

!write(*,*) 'after 30 prior is ', prior
!call beta_pdf( paramvec(31,1)/2.d0+0.5d0 , 2.0d0, 2.0d0, pdf )      ! rhoint
!prior = prior + log(pdf) + log(0.5d0)
!call beta_converter(0.5d0, 0.15d0, a1,b1)
!call beta_pdf( paramvec(31,1), a1, b1, pdf )
!call normal_pdf( paramvec(31,1), 0.0d0, 0.1d0, pdf )                 ! rhoint
!!print*, ' rhoint= ', log(pdf)
!prior = prior + log(pdf)
!write(*,*) 'after 31 prior is ', prior

!call invgamma_pdf_1( 100.0d0 * paramvec(32,1), b1, a1, pdf )      ! sdevint
pdf = logigpdf(100.0d0*paramvec(32,1),0.33d0,2.0d0)
!print*, ' sdevint= ', pdf
prior = prior + pdf! + log(100.0d0)

!write(*,*) 'after 32 prior is ', prior
!call beta_pdf( paramvec(33,1), 2.0d0, 2.0d0, pdf )      ! rhoelast
!prior = prior + log(pdf)

!call invgamma_pdf_1( 100.0d0 * paramvec(34,1), b1, a1, pdf )      ! sdevelast
!pdf = logigpdf(100.0d0*paramvec(34,1),0.33d0,2.0d0)
!prior = prior + pdf


!write(*,*) 'after 34 prior is ', prior
!call beta_pdf( paramvec(35,1), 2.0d0, 2.0d0, pdf )      ! rhoelastw
!prior = prior + log(pdf)

!call invgamma_pdf_1( .1d0 * paramvec(36,1), b1, a1, pdf )      ! stdevelastw
!prior = prior + log(pdf) + log(.1d0)

!write(*,*) 'after 36 prior is ', prior
! MEASUREMENT ERROR PRIORS


! this is mu = 0.1 and sigma = 0.005
!a1 =  0.548465619838185
!b1 =  56.709787515030897
! this is mu = 0.05 and sigma = 0.01
!a1 = 0.033012039987416
!b1 = 14.696938456698652


sf = 100.0d0
!sf = 1.0d0

!b1 = 56.7097875150313
!a1 = 0.229390411951865
a1 =  0.548465619838185
b1 =  56.709787515030897
!call invgamma_pdf_1(sf*paramvec(37,1),  b1 , a1, pdf )  ! me 1

!b1 = 0.000646715014129
!call rayleigh_pdf(sf*paramvec(37,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)
!write(*,*) 'after 37 prior is ', prior

!b1 = 56.7097875150313
!a1 = 0.144688758440556
a1 =  0.548465619838185
b1 =  56.709787515030897
!call invgamma_pdf_1(sf*paramvec(38,1),  b1 , a1, pdf )  ! me 1
!b1 = 0.000513620879143
!call rayleigh_pdf(sf*paramvec(38,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)
!write(*,*) 'after 38 prior is ', prior


b1 = 56.7097875150312944
a1 = 3.2578696357890338
!a1 =  0.548465619838185
!b1 =  56.709787515030897

!b1 = 56.7097875150312944
!a1 = 81.4467408947258633
!call invgamma_pdf_1(sf*paramvec(39,1),  b1 , a1, pdf )  ! me 1
!b1 = 0.002437205478900
!call rayleigh_pdf(sf*paramvec(39,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)
!write(*,*) 'after 39 prior is ', prior

b1 = 56.7097875150312944
a1 = 18.5119795955022681
!a1 =  0.548465619838185
!b1 =  56.709787515030897

!b1 = 56.7097875150312944
!a1 = 462.799489887556717
!call invgamma_pdf_1(sf*paramvec(40,1),  b1 , a1, pdf )  ! me 1
!b1 = 0.005809673510570
!call rayleigh_pdf(sf*paramvec(40,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)
!write(*,*) 'after 40 prior is ', prior

!b1 = 56.7097875150313
!a1 = 0.368668413111026
a1 =  0.548465619838185
b1 =  56.709787515030897
!call invgamma_pdf_1(sf*paramvec(41,1),  b1 , a1, pdf )  ! me 1

!b1 = 0.000819866718311
!call rayleigh_pdf(sf*paramvec(41,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)
!write(*,*) 'after 41 prior is ', prior

!b1= 56.7097875150313
!a1 = 0.0328816862042421
a1 =  0.548465619838185
b1 =  56.709787515030897
!call invgamma_pdf_1(sf*paramvec(42,1),  b1 , a1, pdf )  ! me 1
!b1 = 0.000244851273823
!call rayleigh_pdf(sf*paramvec(42,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)

!write(*,*) 'after 42 prior is ', prior

!b1 = 56.7097875150313
!a1 = 0.282057738101700
a1 =  0.548465619838185
b1 =  56.709787515030897
!call invgamma_pdf_1(sf*paramvec(43,1),  b1 , a1, pdf )  ! me 1
!b1 = 0.000717124078396
!call rayleigh_pdf(sf*paramvec(43,1),  b1 , pdf )  ! me 1
!prior = prior + log(pdf) + log(sf)
!write(*,*) 'after 43 prior is ', prior

end subroutine priorfcn
