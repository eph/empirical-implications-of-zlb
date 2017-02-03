program ghlss_rwmh_driver

  use flap, only : command_line_interface
  use class_model, only: model
  use class_ParallelParticleFilter, only: ParallelParticleFilter
  use class_RandomNumber, only: RandomNumber
  !use class_TemperedParticleFilter, only: ParallelParticleFilter => TemperedParticleFilter 
  use json_module

  implicit none

  include 'mpif.h'

  type(command_line_interface) :: cli
  integer :: error

  character(len=400) :: pmsv_file, cov_file, prior_file
  character(len=400) :: output_dir, output_file, tmp_file
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

  integer :: rank, nproc, mpierror, i , seed, acpt, j, suffix

  integer :: time0, time1, rate

  logical, external :: inbounds

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

  call system_clock(count_rate=rate)


  call cli%init(progname = 'driver_prwmh', &
       version='0.0.1', &
       authors='Chris Gust & Ed Herbst', &
       description='Program for estimating model via particle MCMC')
  call cli%add(switch='--p0', switch_ab='-p0', required=.false., def='input/mean.txt',help='starting value')
  call cli%add(switch='--covariance', switch_ab='-H', required=.false., def='input/cholM.txt',help='Cholesky of covariance for proposal innvoations')
  call cli%add(switch='--prior-file', required=.false., def='inputs/case-0/prior.txt',help='Prior file')
  call cli%add(switch='--scaling',switch_ab='-c', required=.false., def='0.15',help='scaling of innovation')
  call cli%add(switch='--nsim', switch_ab='-n', required=.false., def='50000', help='Length of MCMC chain')
  call cli%add(switch='--npart', required=.false., def='1500000',help='Number of Particles for PF')
  call cli%add(switch='--zlb', required=.false., def='.true.',help='Impose ZLB')
  call cli%add(switch='--output-file', required=.false., def='output.json',help='Output file')
  call cli%add(switch='--seed', required=.false., def='1848',help='Seed for RNG')

  call cli%parse(error=error)
  call cli%get(switch='--zlb', val=zlb)
  call cli%get(switch='-H', val=cov_file)
  call cli%get(switch='--npart', val=npart)
  call cli%get(switch='--nsim', val=nsim)
  call cli%get(switch='--seed', val=seed)
  call cli%get(switch='--scaling', val=c)
  call cli%get(switch='--output-file', val=output_file)
  call cli%get(switch='--p0', val=pmsv_file)
  !call cli%get(switch='--prior-file', val=prior_file)


  dsge = model(zlb)

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

  ! zero out the fixed parameters
   !M(1,:) = 0.0d0
   !M(2,:) = 0.0d0
   M(4,:) = 0.0d0
   M(9,:) = 0.0d0
   M(10,:) = 0.0d0
   !M(11,:) = 0.0d0
   !M(12,:) = 0.0d0
   M(13,:) = 0.0d0   
   !M(14,:) = 0.001d0
   M(16,:) = 0.0d0
   M(23,:) = 0.0d0
   !M(21,:) = 0.0d0
   !M(22,:) = 0.0d0
   !M(28,:) = M(28,:)/10.0d0
   M(29,:) = 0.00d0
   M(31,:) = 0.00d0
   M(33,:) = 0.0d0
   M(34,:) = 0.0d0
   M(35,:) = 0.0d0
   M(36,:) = 0.0d0
   M(37,:) = 0.0d0
   M(38,:) = 0.0d0
   M(39,:) = 0.0d0
   M(40,:) = 0.0d0
   M(41,:) = 0.0d0
   M(42,:) = 0.0d0
   M(43,:) = 0.0d0

   ppf = ParallelParticleFilter(dsge, npart=npart, nproc=nproc, rank=rank, seed=seed+1+rank)
   ppf%adjusted_proposal(104) = .true.
   ppf%adjusted_proposal_std(104) = 1.2d0
   ppf%adjusted_proposal_mu(104,1) = 3.0d0

   lik0 = ppf%lik(p0, rank=rank, nproc=nproc)
   pr0 = dsge%pr(p0)

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


  suffix = index(output_file, '.json')
  tmp_file = output_file(1:suffix-1)//'_para_tmp.txt'




   random_number = RandomNumber(seed=seed)
   if (rank==0) print*,'Initial Log Likelihood:', lik0
   if (rank==0) open(1, file=tmp_file, action='write')

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
        print*,sum(parasim(:,13))/dble(i-1)
        print*,'===================================================================='

        write(1, '(100f)') parasim(:,i)

     end if
  end do


  if (rank==0) then
     close(1)
     call json%create_object(output, 'output')
     call json%add(p, output)


     call json%add(output, 'postsim', postsim)
     call json%add(output, inp)
     do i = 1, dsge%npara
        write(tmp_file, '(I2.2)') i
        call json%add(inp, 'para.'//trim(tmp_file), parasim(i,:))
     end do

     !call json%add(output, 'parasim_file', output_dir//'/parasim1.txt')
     nullify(output)
     call json%print(p, output_file)
  call json%destroy(p)
  end if

  deallocate(eps, u, p0, p1, M, diag_M, parasim, postsim, acptsim)

  call mpi_finalize(mpierror)

end program ghlss_rwmh_driver
