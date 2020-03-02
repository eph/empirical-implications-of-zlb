!---------------------------------------------------------------------------------------------------
!
!PROGRAM: driver_irfs.f90
!> @author Christopher Gust
!> @version 1.0
!> @date 11-17-16
!
!DESCRIPTION:
!> Read in parameters and compute irfs to risk premium or MEI shock.
!_________________________________________________________________________________________________________
program driver_irfs
  use utils, only: read_matrix, write_matrix
  use class_model, only: model

  use flap

  implicit none
  include 'mpif.h'

  type(command_line_interface) :: cli
  type(model) :: m
  logical, parameter :: parallelswitch = .true.
  integer, parameter :: neulererrors = 18
  integer :: i,captirf,nsim,shockindex,dsetswitch
  logical :: zlbswitch,convergence
  double precision, allocatable :: endogirf(:,:)
  double precision, allocatable :: linirf(:,:)
  double precision, allocatable :: euler_errors(:,:)
  double precision :: shockscale
  integer :: rank
  integer :: nproc
  integer :: mpierror, error
  character(len=150) :: arg,param_type
  character (len=250) :: filename


  if (parallelswitch .eqv. .true.) then
     call MPI_init(mpierror)
     call MPI_Comm_size(MPI_COMM_WORLD,nproc,mpierror)
     call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierror)
     write(*,*) 'Hello from processor ', rank, 'I am 1 of ', nproc ,' processes running'
  else
     rank = 0
     write(*,*) 'You are running the serial version of the code.'
  end if


  call cli%init(progname = 'driver_irfs', &
       authors='Chris Gust & Ed Herbst', &
       description='Program for computing impulse response functions.')
  call cli%add(switch='--nsim',switch_ab='-n',required=.false.,def='5000',help='Number of draws for MC integration')
  call cli%add(switch='--capt',switch_ab='-c',required=.false.,def='20',help='Length of IRF')
  call cli%add(switch='--shockindex',switch_ab='-s',required=.false.,def='1',choices='1,2',help='Shock: (1) Risk Premium; (2) MEI')
  call cli%add(switch='--dsetswitch',switch_ab='-d',required=.false.,def='0',choices='0,1',help='Parameters: (0) Mean; (1) Medium')
  call cli%add(switch='--scaleshock', switch_ab='-a',required=.false.,def='2.5d0',help='size of shock in stds')
  call cli%parse(error=error)
  call cli%get(switch='-n',val=nsim)
  call cli%get(switch='-c',val=captirf)
  call cli%get(switch='-s',val=shockindex)
  call cli%get(switch='-d',val=dsetswitch)
  call cli%get(switch='-a',val=shockscale)
  !iniitialize solution details
  zlbswitch = .true.
  m = model(zlbswitch)  
  if (rank .eq. 0) then
     call m%describe()
     write(*,*) 'capt = ', captirf
     write(*,*) 'number of draws for MC integration = ', nsim
     if (shockindex .eq. 1) write(*,*) 'IRF to Risk Premium Shock'
     if (shockindex .eq. 2) write(*,*) 'IRF to MEI Shock'
     write(*,*) '----------------------------'
  end if

  !allocate matrices for IRFs
  allocate(endogirf(m%solution%poly%nvars+m%solution%poly%nexog+2,captirf))
  allocate(linirf(m%solution%poly%nvars+m%solution%poly%nexog+2,captirf))
  allocate(euler_errors(2*neulererrors,captirf))
  
  !get parameters from disk
  if (dsetswitch .eq. 0) then
     param_type = 'mean'
  else
     param_type = 'median'
  end if
  filename = 'input/' // trim(param_type) //  '.txt'
  call read_matrix(filename,m%solution%poly%nparams,1,m%params)

  !solve model
  if (parallelswitch .eqv. .true.) then
     convergence = m%solve_parallel(m%params,nproc,rank)
  else
     convergence = m%solve(m%params)
  end if

  if (convergence .eqv. .false.) then  !if no solution, report this back to disk
     if (rank .eq. 0) write(*,*) 'Failed to converge (driverirf)'     
  else  !if computed solution, simulate and send irf data to disk
     if (rank .eq. 0) then
        write(*,*) 'Successfully solved model (driverirf). Computing IRFs.' 
        call m%simulate_modelirfs(captirf,nsim,shockindex,endogirf,linirf,euler_errors,neulererrors, shockscale)
        !send irfs to disk
        if (zlbswitch .eqv. .true.) then
           write(filename,"(A,I0,A)") 'results/irf/nonlinearirf_' // trim(param_type) // '_', shockindex , '.txt'
        else
           write(filename,"(A,I0,A)") 'results/irf/nonlinearirf_unc_' // trim(param_type) // '_', shockindex , '.txt'
        end if
        call write_matrix(filename,m%solution%poly%nvars+m%solution%poly%nexog+2,captirf,endogirf)    
        write(filename,"(A,I0,A)") 'results/irf/linearirf_' // trim(param_type) // '_', shockindex , '.txt'
        call write_matrix(filename,m%solution%poly%nvars+m%solution%poly%nexog+2,captirf,linirf) 
        write(filename,"(A,I0,A)") 'results/irf/eulers_' // trim(param_type) // '_', shockindex, '.txt'
        call write_matrix(filename,2*neulererrors,captirf,euler_errors)    
     end if
  end if
 
  deallocate(endogirf,linirf,euler_errors)
  call m%cleanup()

  if (parallelswitch .eqv. .true.) then
     call MPI_finalize(mpierror)
  end if
   
end program driver_irfs
