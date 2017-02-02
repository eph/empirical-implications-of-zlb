program ghlss_smoother_driver

  use class_model, only: model
  use class_ParticleSmoother, only: ParticleSmoother
  use class_RandomNumber, only: RandomNumber

  use json_module

  use flap, only: command_line_interface

  implicit none

  include 'mpif.h'

  type(command_line_interface) :: cli

  type(Model) :: dsge
  type(ParticleSmoother) :: ppf
  
  double precision, allocatable :: para(:)
  character(len=:), allocatable :: pmsv_file, out_file
  double precision, allocatable :: filtered_states(:,:,:), smoothed_states(:,:,:), filtered_shocks(:,:,:), smoothed_shocks(:,:,:)
  double precision, allocatable :: true_smooth_states(:,:), true_smooth_shocks(:,:), true_smooth_var(:,:,:)
  double precision, allocatable :: true_filter_states(:,:), true_filter_shocks(:,:), true_filter_var(:,:,:)

  double precision :: lik0 
  integer :: nsave, npart, seed, error
  character(len=500) :: arg, charsimi, charvari
  integer :: i, n,j

  type(json_core) :: json 
  type(json_value), pointer :: p, inp, output, sim_i, sim_i_s

  integer :: rank, nproc, mpierror

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)


  call cli%init(progname = 'driver_smoother', &
       version='0.0.1', &
       authors='Chris Gust & Ed Herbst', &
       description='Runs the particle filter and smoother for a given set of parameters')


  call cli%add(switch='--p0', switch_ab='-p0', required=.false., def='inputs/case-0/pmsv00.txt',help='Parameters values for model.')
  call cli%add(switch='--npart', required=.false., def='1500000',help='Number of Particles for PF')
  call cli%add(switch='--output-file', required=.false., def='output.json',help='Output file')
  call cli%add(switch='--nsave',required=.false.,def='10',help='Number of draws to save')
  call cli%add(switch='--seed', required=.false., def='1848',help='Seed for RNG')
  call cli%parse(error=error)
  call cli%get(switch='--npart', val=npart)
  call cli%get(switch='--nsave', val=nsave)
  call cli%get(switch='--p0', val=pmsv_file)
  call cli%get(switch='--seed', val=seed)
  call cli%get(switch='--output-file', val=out_file)


  ! npart = 150000
  ! out_file = 'test.json'
  ! do i = 1, command_argument_count()
  !    call get_command_argument(i, arg)

  !    select case(arg)
  !    case('--pmsv','-p0')
  !       call get_command_argument(i+1,arg)
  !       pmsv_file = arg
  !    case('--outfile')
  !       call get_command_argument(i+1,arg)
  !       out_file = arg
  !    case('--npart')
  !       call get_command_argument(i+1,arg)
  !       read(arg, '(i)') npart

  !    end select
  ! end do

  dsge = model(.true.)
  print*,'dsge%T', dsge%T
  allocate(para(dsge%npara))

  open(49, file=pmsv_file, action='read')
  do i = 1,dsge%npara
     read(49, *) para(i)
  end do
  close(49)

  ppf = ParticleSmoother(dsge, npart=npart, seed=seed)
  ppf%adjusted_proposal(104) = .true.
  ppf%adjusted_proposal_std(104) = 1.2d0
  ppf%adjusted_proposal_mu(104,1) = 3.0d0

  nsave = 10
  allocate(filtered_states(0:dsge%T, dsge%nvars, nsave))
  allocate(smoothed_states(0:dsge%T, dsge%nvars, nsave))

  allocate(filtered_shocks(0:dsge%T, dsge%nexog, nsave))
  allocate(smoothed_shocks(0:dsge%T, dsge%nexog, nsave))


  call ppf%filter_and_smooth(para, nsave, filtered_states, smoothed_states, &
       filtered_shocks, smoothed_shocks)

  call json%create_object(p,'')
  call json%create_object(inp,'input')
  call json%add(p, inp)
  call json%add(inp, 'p0', para)
  call json%add(inp, 'npart', ppf%npart)
  call json%add(inp, 'seed', seed)

  nullify(inp)

  call json%create_object(output, 'output')
  call json%add(p, output)

  do i = 1, nsave
     write(charsimi, '(I3.3)') i 

     call json%create_object(sim_i, 'sim_'//trim(charsimi))
     call json%add(output, sim_i)

     call json%create_object(sim_i_s, 'filtered_states')
     call json%add(sim_i, sim_i_s) 
     do j = 1, dsge%nvars
        write(charvari, '(I2.2)') j
        call json%add(sim_i_s,'endogvar_'//trim(charvari), filtered_states(:,j,i))
     end do
     nullify(sim_i_s)

     call json%create_object(sim_i_s, 'filtered_shocks')
     call json%add(sim_i, sim_i_s) 
     do j = 1, dsge%nexog
        write(charvari, '(I2.2)') j 
        call json%add(sim_i_s,'exogvar_'//trim(charvari), filtered_shocks(:,j,i))
     end do
     nullify(sim_i_s)

     call json%create_object(sim_i_s, 'smoothed_shocks')
     call json%add(sim_i, sim_i_s) 
     do j = 1, dsge%nexog
        write(charvari, '(I2.2)') j 
        call json%add(sim_i_s,'exogvar_'//trim(charvari), smoothed_shocks(:,j,i))
     end do
     nullify(sim_i_s)

     call json%create_object(sim_i_s, 'smoothed_states')
     call json%add(sim_i, sim_i_s) 
     do j = 1, dsge%nvars
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
