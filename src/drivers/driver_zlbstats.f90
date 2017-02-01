!---------------------------------------------------------------------------------------------------
!
!PROGRAM: driver_zlbstats.f90
!> @author Christopher Gust
!> @version 1.0
!> @date 11-16-16
!
!DESCRIPTION:
!> Loop through parameter draws and calculate zlb probability and zlb duration data used to construct histograms.
!_________________________________________________________________________________________________________
program driver_zlbstats
  use utils, only: read_matrix, write_matrix
  use class_model, only: model
  use simulate_model, only: compute_zlbstats
  implicit none
  include 'mpif.h'

  type(model) :: m
  logical, parameter :: parallelswitch = .true.
  integer, parameter :: capt_sample = 125
  integer :: capt_simdata,ndsets
  integer :: i,id,number_samples,nspellmax,number_zlbspells,ndsets_used
  logical :: nonlinearswitch,zlbswitch,convergence
  double precision :: nspells_id(1)
  double precision, allocatable :: modeldata(:,:)
  double precision, allocatable :: nomrdata(:)
  double precision, allocatable :: zlbfrequency(:)
  integer, allocatable :: zlbduration(:)
  integer :: rank
  integer :: nproc
  integer :: mpierror
  character(len=150) :: arg
  character (len=250) :: filename

  if (parallelswitch .eq. .true.) then
     call MPI_init(mpierror)
     call MPI_Comm_size(MPI_COMM_WORLD,nproc,mpierror)
     call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierror)
     write(*,*) 'Hello from processor ', rank, 'I am 1 of ', nproc ,' processes running'
  else
     rank = 0
     write(*,*) 'You are running the serial version of the code.'
  end if

  capt_simdata = 100000
  ndsets = 1

  do i = 1, command_argument_count()
     call get_command_argument(i, arg)
     select case(arg)
     case ('--ndsets', '-n')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') ndsets
     case('--capt', '-c')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') capt_simdata
     end select
  end do

  !iniitialize solution details
  zlbswitch = .true.
  m = model(zlbswitch)  
  if (rank .eq. 0) then
     call m%describe()
     write(*,*) 'capt = ', capt_simdata
     write(*,*) 'number of dsets = ', ndsets
     write(*,*) '----------------------------'
  end if

  !allocate matrices for zlb frequency and duration data
  allocate(modeldata(m%solution%poly%nvars+2*m%solution%poly%nexog,capt_simdata),nomrdata(capt_simdata))
  number_samples = int(capt_simdata/capt_sample)
  allocate(zlbfrequency(number_samples))
  zlbfrequency = 0.0d0
  nspellmax = int(0.1*capt_simdata)
  allocate(zlbduration(nspellmax))
  zlbduration = 0
  
  ndsets_used = 0
  do id = 1,ndsets

     !get parameters from disk
     write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/final-final/by-draw/parasim', id-1 , '.txt'
     call read_matrix(filename,m%solution%poly%nparams,1,m%params)
 
     !solve model
     if (parallelswitch .eq. .true.) then
        convergence = m%solve_parallel(m%params,nproc,rank)
     else
        convergence = m%solve(m%params)
     end if

     if (convergence .eq. .false.) then  !if no solution, report this back to disk
        if (rank .eq. 0) then
           write(*,*) 'Failed to converge for parameters when id = ', id
           write(*,*) '----------------------------------'
           write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/final-final/zlbstat-results/nzlbspells', id-1 , '.txt'
           nspells_id(1) = 0.0d0
           call write_matrix(filename,1,1,nspells_id)
        end if
     else  !if computed solution, simulate and send zlb data to disk
        if (rank .eq. 0) then
           write(*,*) 'Successfully solved model when id = ', id 
           write(*,*) '---------------------------------'

           nonlinearswitch = .true.
           call m%simulate_modeldata(capt_simdata,nonlinearswitch,modeldata,seed=1221+id)
           nomrdata = 400.0d0*log(modeldata(5,:))
           call compute_zlbstats(nomrdata,capt_simdata,capt_sample,number_samples,nspellmax,number_zlbspells,zlbduration,&
                zlbfrequency)
           ndsets_used = ndsets_used + 1
           
           write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/zlbstat-results/nzlbspells' , id-1 , '.txt'
           nspells_id(1) = dble(number_zlbspells)
           call write_matrix(filename,1,1,nspells_id)           
           write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/zlbstat-results/zlbduration', id-1 , '.txt'
           call write_matrix(filename,number_zlbspells,1,dble(zlbduration))
           write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/zlbstat-results/zlbfrequency', id-1 , '.txt'
           call write_matrix(filename,number_samples,1,zlbfrequency)    
        end if
     end if
  end do

  if (rank .eq. 0) write(*,*) 'Number of dsets for zlb statistics = ', ndsets_used  
  
  deallocate(zlbduration,zlbfrequency,modeldata,nomrdata)
  call m%cleanup()

  if (parallelswitch .eq. .true.) then
     call MPI_finalize(mpierror)
  end if
   
end program driver_zlbstats
