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
  implicit none
  include 'mpif.h'

  type(model) :: m
  logical, parameter :: parallelswitch = .true.
  integer :: capt_simdata,ndsets
  integer :: i,id,ndsets_used,modelindex
  logical :: nonlinearswitch,zlbswitch,convergence
  double precision :: dum(1)
  double precision, allocatable :: modeldata(:,:)
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

  capt_simdata = 1000
  ndsets = 1
  modelindex = 0

  do i = 1, command_argument_count()
     call get_command_argument(i, arg)
     select case(arg)
     case ('--ndsets', '-n')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') ndsets
     case('--capt', '-c')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') capt_simdata
     case('--model','-m')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') modelindex
     end select
  end do

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
           if (modelindex .eq. 0) then
              write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/final-final/zlbstat-results/modeldata_file', id-1 , '.txt'
           else
              write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/final-final/zlbstat-results/modeldata_unc_file', id-1 , '.txt'
           end if
           call write_matrix(filename,1,1,dum)
           if (modelindex .eq. 0) then
              write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/final-final/zlbstat-results/modeldata_linear_file', id-1 , '.txt'
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
           if (modelindex .eq. 0) then
              write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/zlbstat-results/modeldata_file', id-1 , '.txt'
           else
              write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/zlbstat-results/modeldata_unc_file', id-1 , '.txt'
           end if
           call write_matrix(filename,m%solution%poly%nvars+2*m%solution%poly%nexog,capt_simdata,modeldata)  
           if (modelindex .eq. 0) then
              nonlinearswitch = .false.
              call m%simulate_modeldata(capt_simdata,nonlinearswitch,modeldata,seed=1221+id)
              write(filename,"(A,I4.4,A)") '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/zlbstat-results/modeldata_linear_file', id-1 , '.txt'
              call write_matrix(filename,m%solution%poly%nvars+2*m%solution%poly%nexog,capt_simdata,modeldata) 
           end if
        end if
     end if
  end do

  if (rank .eq. 0) write(*,*) 'Number of dsets for selected moments table = ', ndsets_used  
  
  deallocate(modeldata)
  call m%cleanup()

  if (parallelswitch .eq. .true.) then
     call MPI_finalize(mpierror)
  end if
   
end program driver_simdata
