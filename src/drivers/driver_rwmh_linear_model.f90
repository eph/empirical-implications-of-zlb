program ghlss_rwmh_driver

  !use flap, only : command_line_interface
  use class_model, only: model
  use class_ParallelParticleFilter, only: ParallelParticleFilter
  use class_RandomNumber, only: RandomNumber

  use json_module

  implicit none

  include 'mpif.h'

  !type(command_line_interface) :: cli
  integer :: error

  character(len=400) :: arg     ! get_command_argument needs an allocated character array
  character(len=:), allocatable :: pmsv_file, cov_file, prior_file
  character(len=:), allocatable :: output_dir
  integer :: nsim, trial, npart

  double precision :: c, lik0, lik1, pr0, pr1, alp
  double precision, allocatable :: parasim(:,:), acptsim(:), postsim(:), M(:,:)
  double precision, allocatable :: eps(:,:), u(:,:), p0(:), p1(:), diag_M(:)

  type(model) :: dsge
  type(ParallelParticleFilter) :: ppf
  type(RandomNumber) :: random_number

  type(json_core) :: json 
  type(json_value), pointer :: p, inp, output

  logical :: convergence, zlb

  integer :: rank, nproc, mpierror, i , seed, acpt, j

  integer :: time0, time1, rate

  logical, external :: inbounds

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

  call system_clock(count_rate=rate)

  nsim = 50000
  pmsv_file = 'none'
  cov_file = '/msu/scratch3/m1cjg01/aer_revision_ed/input/proposal_variance_eph_mu.txt'
  output_dir = '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/other/pmcmc-linear-analysis/'
  prior_file = ''
  c = 0.1d0
  seed = 1847
  trial = 1
  npart = 1500000

  do i = 1, command_argument_count()
     call get_command_argument(i, arg)
     select case(arg)
     case ('--starting-value-file', '-p0')
        call get_command_argument(i+1, arg)
        pmsv_file = arg
     case ('--covariance-file', '-M')
        call get_command_argument(i+1, arg)
        cov_file = arg
     case ('--prior-file', '-pr')
        call get_command_argument(i+1, arg)
        prior_file = arg
     case ('--nsim', '-n')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') nsim
     case('--scaling', '-c')
        call get_command_argument(i+1, arg)
        read(arg, '(f)') c
     case('--seed', '-s')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') seed
     case('--trial','-i')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') trial
     end select
  end do

  dsge = model()

  allocate(parasim(dsge%npara, nsim), acptsim(nsim), postsim(nsim))
  allocate(p0(dsge%npara), p1(dsge%npara), M(dsge%npara,dsge%npara), diag_M(dsge%npara))

  if (pmsv_file=='none') then
     p0 = dsge%params
  else
     if (rank==0) print*,'pmsv_file', pmsv_file
     open(199, file=pmsv_file, action='read')
     do i = 1, dsge%npara
        read(199,*) p0(i)
     end do
     close(199)
  end if

  open(199, file=cov_file, action='read')
  do i = 1, dsge%npara
     read(199, *) M(i,:)
  end do
  close(199)

  ppf = ParallelParticleFilter(dsge, npart=npart, nproc=nproc, rank=rank, seed=rank)
  ppf%adjusted_proposal(104) = .true.
  ppf%adjusted_proposal_std(104) = 1.2d0
  ppf%adjusted_proposal_mu(104,1) = 3.0d0

  lik0 = ppf%lik(p0, rank=rank, nproc=nproc)
  pr0 = dsge%pr(p0)
  !call priorfcn(dsge%npara, p0, pr0)

  if (rank==0) then
     call json%create_object(p, '')
     call json%create_object(inp, 'input')
     call json%add(p, inp)

     call json%add(inp, 'name', dsge%name)
     call json%add(inp, 'p0', p0)
     call json%add(inp, 'lik0', lik0)
     call json%add(inp, 'pr0', pr0)
     call json%add(inp, 'c', c)
     !call json%add(inp, 'covariance', M)
     call json%add(inp, 'nsim', nsim)
     call json%add(inp, 'npart', npart)
     nullify(inp)

  end if


  random_number = RandomNumber(seed=seed)
  if (rank==0) print*,'Initial Log Likelihood:', lik0

  allocate(eps(dsge%npara, nsim), u(nsim,1))
  eps = random_number%norm_rvs(dsge%npara, nsim)
  u = random_number%uniform_rvs(nsim,1)
  acpt = 0

  do i = 2, nsim
     if (rank==0) then 
        call system_clock(time0)
        !p1 = p0 + c * diag_M * eps(:,i)
        p1 = p0
        call dgemv('n', dsge%npara, dsge%npara, c, M, dsge%npara, eps(:,i), 1, 1.0d0, p1, 1)
        
     end if

     call mpi_bcast(p1, dsge%npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
     call mpi_barrier(MPI_COMM_WORLD, mpierror)

     lik1 = ppf%lik(p1, rank=rank, nproc=nproc) 
     pr1 = dsge%pr(p1)

     alp = exp(lik1 + pr1 - lik0 - pr0)

     if (rank==0) then
        print*,''
        print*,'===================================================================='
        print*,'ITERATION ', i
        print*,'Current likelihood + prior:', lik0, pr0
        print*,'Proposed likelihood + prior:', lik1, pr1
     end if

     if (u(i,1) < alp) then 
        lik0 = lik1
        pr0 = pr1
        p0 = p1
        acpt = acpt + 1
        if (rank==0) print*,'ACCEPTED'           
     end if

     if (rank==0) then
        parasim(:,i) = p0
        postsim(i) = lik0 + pr0
        call system_clock(time1)
        print*, 'Current Acceptance Rate', acpt/real(i-1)
        print*,'iteration', i, ' took ', (time1-time0)/real(rate) 
        print*,'===================================================================='
     end if
  end do


  if (rank==0) then
     call json%create_object(output, 'output')
     call json%add(p, output)

     !call json%add(output, 'parasim', parasim)
     call json%add(output, 'postsim', postsim)
     call json%add(output, 'parasim_file', output_dir//'/parasim1.txt')
     nullify(output)
     call json%print(p, output_dir//'/output1.json')
     open(1, file=output_dir//'/parasim.txt', action='write')
     do i = 1, nsim
        write(1, '(100f)') parasim(:,i)
     end do
     close(1)
  call json%destroy(p)
  end if

  deallocate(eps, u, p0, p1, M, diag_M, parasim, postsim, acptsim)

  call mpi_finalize(mpierror)

  ! call cli%init(progname='ghlss_rwmh_driver', &
  !      version='0.0.3', &
  !      authors='Chris Gust, Matt Smith, and Ed Herbst', &
  !      description='Bayesian Estimation of a Nonlinear DSGE Model')

  ! call cli%add(switch='--nsim', switch_ab='-n', &
  !      help='Integer Length of chain generated by MCMC algorithm.', &
  !      act='store', required=.false., def='20000', error=error)

  ! call cli%add(switch='--scaling', switch_ab='-c', &
  !      help='Real number for scaling the variance of the proposal dist.', &
  !      act='store', required=.false., def='0.5')

  ! call cli%add(switch='--pmsv', switch_ab='-p0', &
  !      help='Path to file containing starting values for parameters.', &
  !      act='store', required=.false., def='p0.txt', error=error)

  ! call cli%add(switch='--covariance', switch_ab='-M', &
  !      help='Path to file containing the (unscaled)'// &
  !      ' covariance for the proposal.', &
  !      act='store', required=.false., def='covariance.txt', error=error)

  ! call cli%add(switch='--prior', switch_ab='-pr', &
  !      help='Path to file containing the prior.', &
  !      act='store', required=.false., def='prior.txt', error=error)

  ! call cli%add(switch='--output-dir', switch_ab='-od', &
  !      help='Path to desired output directory.', &
  !      act='store', required=.false., def='.')

  ! call cli%get(switch='-c', val=c, error=error)
  ! call cli%get(switch='-n', val=nsim, error=error)
  ! call cli%get(switch='-p0', val=pmsv_file, error=error)
  ! call cli%get(switch='-M', val=cov_file, error=error)

  ! model = Model()

  ! allocate(parasim(model%npara,nsim), acptsim(nsim), postsim(nsim))

  ! do i = 1, nsim

  ! end do

  ! deallocate(parasim, acptsim, postsim)

end program ghlss_rwmh_driver
