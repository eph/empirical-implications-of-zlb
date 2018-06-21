!---------------------------------------------------------------------------------------------------
!
!PROGRAM: driver_simdata.f90
!> @author Christopher Gust
!> @version 1.0
!> @date 11-21-16
!
!DESCRIPTION:
!> Loop through parameter draws and calculate zlb probability and zlb duration data used to construct histograms.
!_________________________________________________________________________________________________________
program driver_simdata
  use utils, only: read_matrix, write_matrix
  use class_model, only: model
  use simulate_model, only: compute_zlbstats

  use flap 

  implicit none
  include 'mpif.h'

  type(command_line_interface) :: cli
  type(model) :: m
  logical, parameter :: parallelswitch = .false.
  integer :: capt_simdata,ndsets,modelindex
  integer :: i,id,ndsets_used
  logical :: nonlinearswitch,zlbswitch,convergence
  double precision :: dum(1)
  double precision, allocatable :: modeldata(:,:)
  integer :: rank
  integer :: nproc
  integer :: mpierror, error
  character(len=150) :: arg
  character (len=250) :: filename

  call cli%init(progname = 'driver_selectmoments', &
       authors='Chris Gust & Ed Herbst', &
       description='Program for computing ZLB probabilities and durations from simulated data.')
  call cli%add(switch='--ndsets', switch_ab='-n',required=.false.,def='1',help='Number of dsets')
  call cli%add(switch='--capt', switch_ab='-c',required=.false.,def='1000', help='Simulation Length')
  call cli%add(switch='--model',switch_ab='-m',required=.false.,def='0',choices='0,1',help='Model: (0) ZLB; (1) NOZLB')
  call cli%parse(error=error)
  call cli%get(switch='-n',val=ndsets)
  call cli%get(switch='-c',val=capt_simdata)
  call cli%get(switch='-m',val=modelindex)

  if (parallelswitch .eqv. .true.) then
     call MPI_init(mpierror)
     call MPI_Comm_size(MPI_COMM_WORLD,nproc,mpierror)
     call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierror)
     write(*,*) 'Hello from processor ', rank, 'I am 1 of ', nproc ,' processes running'
  else
     rank = 0
     write(*,*) 'You are running the serial version of the code.'
  end if

  !iniitialize solution details
  if (modelindex .eq. 0) then
     zlbswitch = .true. 
  else
     zlbswitch = .false.
  end if
  m = model(zlbswitch)  
  if (rank .eq. 0) then
     call m%describe()
     write(*,*) 'capt = ', capt_simdata
     write(*,*) 'number of dsets = ', ndsets
     write(*,*) '----------------------------'
  end if

  !allocate matrices for zlb frequency and duration data
  allocate(modeldata(m%solution%poly%nvars+2*m%solution%poly%nexog,capt_simdata))
  
  dum = 0.0d0
  ndsets_used = 0
  do id = 1,ndsets

     !get parameters from disk
     write(filename,"(A,I4.4,A)") 'results/thinned_posterior/parasim', id-1 , '.txt'
     call read_matrix(filename,m%solution%poly%nparams,1,m%params)
 
     !solve model
     if (parallelswitch .eqv. .true.) then
        convergence = m%solve_parallel(m%params,nproc,rank)
     else
        convergence = m%solve(m%params)
     end if

     if (convergence .eqv. .false.) then  !if no solution, report this back to disk
        if (rank .eq. 0) then
           write(*,*) 'Failed to converge for parameters when id = ', id
           write(*,*) '----------------------------------'
           if (modelindex .eq. 0) then
              write(filename,"(A,I4.4,A)") 'results/zlb/modeldata_file', id-1 , '.txt'
           else
              write(filename,"(A,I4.4,A)") 'results/zlb/modeldata_unc_file', id-1 , '.txt'
           end if
           call write_matrix(filename,1,1,dum)
           if (modelindex .eq. 0) then
              write(filename,"(A,I4.4,A)") 'reuslts/zlb/modeldata_linear_file', id-1 , '.txt'
              call write_matrix(filename,1,1,dum)
           end if


        end if
     else  !if computed solution, simulate and send zlb data to disk
        if (rank .eq. 0) then
           write(*,*) 'Successfully solved model when id = ', id 
           write(*,*) '---------------------------------'
           ndsets_used = ndsets_used + 1
           nonlinearswitch = .true.
           call m%simulate_modeldata(capt_simdata,nonlinearswitch,modeldata,seed=1221+id)
           write(filename,"(A,I4.4,A)") 'results/zlb/modeldata_file', id-1 , '.txt'
           call write_matrix(filename,m%solution%poly%nvars+2*m%solution%poly%nexog,capt_simdata,modeldata)  
           if (modelindex .eq. 0) then
              nonlinearswitch = .false.
              call m%simulate_modeldata(capt_simdata,nonlinearswitch,modeldata,seed=1221+id)
              write(filename,"(A,I4.4,A)") 'results/zlb/modeldata_linear_file', id-1 , '.txt'
              call write_matrix(filename,m%solution%poly%nvars+2*m%solution%poly%nexog,capt_simdata,modeldata) 
           end if
        end if
     end if
  end do

  if (rank .eq. 0) write(*,*) 'Number of dsets for selected moments table = ', ndsets_used  
  
  deallocate(modeldata)
  call m%cleanup()

  if (parallelswitch .eqv. .true.) then
     call MPI_finalize(mpierror)
  end if
   
end program driver_simdata
