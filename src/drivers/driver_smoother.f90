program ghlss_smoother_driver

  use class_model, only: model
  use class_ParticleSmoother, only: ParticleSmoother
  use class_RandomNumber, only: RandomNumber

  use json_module


  implicit none

  include 'mpif.h'

  type(Model) :: linear_dsge
  type(ParticleSmoother) :: ppf
  
  double precision, allocatable :: para(:)
  character(len=:), allocatable :: pmsv_file, out_file
  double precision, allocatable :: filtered_states(:,:,:), smoothed_states(:,:,:), filtered_shocks(:,:,:), smoothed_shocks(:,:,:)
  double precision, allocatable :: true_smooth_states(:,:), true_smooth_shocks(:,:), true_smooth_var(:,:,:)
  double precision, allocatable :: true_filter_states(:,:), true_filter_shocks(:,:), true_filter_var(:,:,:)

  double precision :: lik0 
  integer :: nsave, npart
  character(len=500) :: arg, charsimi, charvari
  integer :: i, n,j

  type(json_core) :: json 
  type(json_value), pointer :: p, inp, output, sim_i, sim_i_s

  integer :: rank, nproc, mpierror

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

  npart = 150000
  out_file = 'test.json'
  do i = 1, command_argument_count()
     call get_command_argument(i, arg)

     select case(arg)
     case('--pmsv','-p0')
        call get_command_argument(i+1,arg)
        pmsv_file = arg
     case('--outfile')
        call get_command_argument(i+1,arg)
        out_file = arg
     case('--npart')
        call get_command_argument(i+1,arg)
        read(arg, '(i)') npart
     end select
  end do

  linear_dsge = model(.true.)
  print*,'linear_dsge%T', linear_dsge%T
  allocate(para(linear_dsge%npara))

  open(49, file=pmsv_file, action='read')
  do i = 1,linear_dsge%npara
     read(49, *) para(i)
  end do
  close(49)
  
 ! allocate(true_smooth_states(linear_dsge%T, linear_dsge%nvars), & 
 !      true_smooth_shocks(linear_dsge%T, linear_dsge%neps), &
 !      true_smooth_var(linear_dsge%nvars, linear_dsge%nvars, linear_dsge%T))

 !1 allocate(true_filter_states(linear_dsge%T, linear_dsge%nvars), & 
 !1      true_filter_shocks(linear_dsge%T, linear_dsge%neps), &
 !1      true_filter_var(linear_dsge%nvars, linear_dsge%nvars, linear_dsge%T))

  ! print*,linear_dsge%lik(para)

  ! call linear_dsge%lik_states(para, lik0, true_filter_states,true_filter_var,true_filter_shocks,true_smooth_states,true_smooth_var,true_smooth_shocks)
  ! print*,linear_dsge%DD(5)
  ! print*,'initial likelihood = ', lik0, linear_dsge%lik(para)

  ! open(111, file='pmax_filtered_true.txt', action='write')
  ! do n = 1, linear_dsge%nvars
  !    write(111, '(1000f)') true_filter_states(:,n)
  ! end do

  ! close(111)


  ! open(111, file='pmax_smoothed_true.txt', action='write')
  ! do n = 1, linear_dsge%nvars
  !    write(111, '(1000f)') true_smooth_states(:,n)
  ! end do

  ! close(111)
  ! deallocate(true_smooth_states, true_smooth_shocks, true_smooth_var)
  ! deallocate(true_filter_states, true_filter_shocks, true_filter_var)

  ppf = ParticleSmoother(linear_dsge, npart=npart, seed=1848)
  ppf%adjusted_proposal(104) = .true.
  ppf%adjusted_proposal_std(104) = 1.2d0
  ppf%adjusted_proposal_mu(104,1) = 3.0d0

  nsave = 10
  allocate(filtered_states(0:linear_dsge%T, linear_dsge%nvars, nsave))
  allocate(smoothed_states(0:linear_dsge%T, linear_dsge%nvars, nsave))

  allocate(filtered_shocks(0:linear_dsge%T, linear_dsge%nexog, nsave))
  allocate(smoothed_shocks(0:linear_dsge%T, linear_dsge%nexog, nsave))


  call ppf%filter_and_smooth(para, nsave, filtered_states, smoothed_states, &
       filtered_shocks, smoothed_shocks)

  call json%create_object(p,'')
  call json%create_object(inp,'input')
  call json%add(p, inp)
  call json%add(inp, 'p0', para)
  call json%add(inp, 'npart', ppf%npart)
  nullify(inp)

  call json%create_object(output, 'output')
  call json%add(p, output)

  do i = 1, nsave
     write(charsimi, '(I3.3)') i 

     call json%create_object(sim_i, 'sim_'//trim(charsimi))
     call json%add(output, sim_i)

     call json%create_object(sim_i_s, 'filtered_states')
     call json%add(sim_i, sim_i_s) 
     do j = 1, linear_dsge%nvars
        write(charvari, '(I2.2)') j
        call json%add(sim_i_s,'endogvar_'//trim(charvari), filtered_states(:,j,i))
     end do
     nullify(sim_i_s)

     call json%create_object(sim_i_s, 'filtered_shocks')
     call json%add(sim_i, sim_i_s) 
     do j = 1, linear_dsge%nexog
        write(charvari, '(I2.2)') j 
        call json%add(sim_i_s,'exogvar_'//trim(charvari), filtered_shocks(:,j,i))
     end do
     nullify(sim_i_s)

     call json%create_object(sim_i_s, 'smoothed_shocks')
     call json%add(sim_i, sim_i_s) 
     do j = 1, linear_dsge%nexog
        write(charvari, '(I2.2)') j 
        call json%add(sim_i_s,'exogvar_'//trim(charvari), smoothed_shocks(:,j,i))
     end do
     nullify(sim_i_s)

     call json%create_object(sim_i_s, 'smoothed_states')
     call json%add(sim_i, sim_i_s) 
     do j = 1, linear_dsge%nvars
        write(charvari, '(I2.2)') j 
        call json%add(sim_i_s,'endogvar_'//trim(charvari), smoothed_states(:,j,i))
     end do
     nullify(sim_i_s)
     nullify(sim_i)
  end do
  nullify(output)

  call json%print(p, out_file)
  call json%destroy(p)

end program
