module class_ParallelParticleFilter

  use class_model, only: model
  use class_ParticleSystem, only: ParticleSystem
  use class_RandomNumber, only: RandomNumber

  implicit none

  include 'mpif.h'

  type ParallelParticleFilter

     integer :: npart = 1000
     integer :: nproc = 1
     integer :: nlocalpart, nforeignpart, naltpart

     integer :: initialization_T = 1

     logical, allocatable :: adjusted_proposal(:)
     double precision, allocatable :: adjusted_proposal_mu(:,:), adjusted_proposal_std(:)

     type(Model) :: m

     type(RandomNumber) :: rng



   contains

     !     procedure :: filter
     !     procedure :: smoother
     procedure :: lik 

  end type ParallelParticleFilter

  interface ParallelParticleFilter
     module procedure new_ParallelParticleFilter
  end interface ParallelParticleFilter

contains

  type(ParallelParticleFilter) function new_ParallelParticleFilter & 
       (m, npart, seed, nproc, rank) result(ppf)

    class(Model) :: m

    integer, optional, intent(in) :: npart, nproc, seed, rank

    integer :: rng_seed
    

    if (present(npart)) ppf%npart = npart
    if (present(npart)) ppf%nproc = nproc
    !if (present(nlocalpart)) ppf%nlocalpart = nlocalpart
    !if (present(nforeignpart)) ppf%nforeignpart = nforeignpart
    !if (present(naltpart)) ppf%naltpart = naltpart

    ppf%nlocalpart = ppf%npart / ppf%nproc
    ppf%nforeignpart = ppf%nlocalpart / 1.5d0
    ppf%naltpart = ppf%nlocalpart / 2.0d0

    if (rank==0) then
       print *, 'Initializing a parallel particle filter '
       print *, 'nlocalpart = ', ppf%nlocalpart, 'nproc = ', ppf%nproc
       print *, 'nforeignpart = ', ppf%nforeignpart, 'naltpart = ', ppf%naltpart
    end if

    if (present(seed)) then
       rng_seed = seed
    else
       rng_seed = 0
    end if

    ppf%rng = RandomNumber(seed=rng_seed)
    ppf%m = m

    allocate(ppf%adjusted_proposal(ppf%m%T), ppf%adjusted_proposal_mu(ppf%m%T, ppf%m%nexog), &
         ppf%adjusted_proposal_std(ppf%m%T))

    ppf%adjusted_proposal = .false.
    ppf%adjusted_proposal_mu = 0.0d0
    ppf%adjusted_proposal_std = 1.0d0

  end function new_ParallelParticleFilter


  double precision function lik(ppf, para, rank, nproc, save_states) 
    class(ParallelParticleFilter) :: ppf

    double precision, intent(in) :: para(ppf%m%npara)

    integer, intent(in) :: rank
    integer, intent(in) :: nproc
    logical, optional, intent(in) :: save_states

    logical :: converged 


    type(ParticleSystem) :: old_local, old_copy, old_foreign, new_local, alt_foreign, alt_copy
    double precision :: incwt, py, py_across_procs, relative_weight, py_vector(nproc)

    double precision :: ess, min_ess, randu(ppf%m%T,1), effective_procs
    integer :: foreignpartstart, ranked_procs(nproc), partner
    logical :: remaining_procs(nproc)
    double precision :: endogsteady(ppf%m%nvars), shocks(ppf%m%nexog, ppf%nlocalpart)
    double precision :: minwt, maxwt

    integer :: i, j, k, t

    double precision :: is_sampling_std, is_sampling_mu(ppf%m%nexog), kap

    integer :: explode_count

    integer :: regime(2)
    integer :: left, right
    integer :: mpierror
    integer :: mpi_recv_status1(MPI_Status_Size)
    integer :: mpi_recv_status2(MPI_Status_Size)
    integer :: mpi_recv_status3(MPI_Status_Size)
    integer :: mpi_recv_status4(MPI_Status_Size)


    integer :: mpi_send_status1(MPI_Status_Size)
    integer :: mpi_send_status2(MPI_Status_Size)
    integer :: mpi_send_status3(MPI_Status_Size)
    integer :: mpi_send_status4(MPI_Status_Size)

    integer :: mpi_req_recv1
    integer :: mpi_req_recv2
    integer :: mpi_req_recv3
    integer :: mpi_req_recv4

    integer :: mpi_req_send1
    integer :: mpi_req_send2
    integer :: mpi_req_send3
    integer :: mpi_req_send4

    character(len=5) :: rb, char_rank, char_t
    character(len=200) :: out_str
    logical :: save_states_opt

    save_states_opt = .false.
    if (present(save_states)) save_states_opt = save_states

    converged = ppf%m%solve(para, nproc, rank)
    !if (rank==0) print*,'entering particle filter'
    lik = 0.0d0
    call mpi_barrier(MPI_COMM_WORLD, mpierror)

    if (converged==.false.) then
       lik = -1000000000000.0d0
       return 
    end if

    randu = ppf%rng%uniform_rvs(ppf%m%T, 1)


    right = modulo(rank+1, nproc)
    left = modulo(rank-1, nproc) 
    foreignpartstart = ppf%nlocalpart - ppf%nforeignpart

    ! initialize part
    old_local = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%nlocalpart)
    new_local = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%nlocalpart)

    old_copy = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%nforeignpart)
    old_foreign = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%nforeignpart)

    alt_copy = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%naltpart)
    alt_foreign = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%naltpart)

    endogsteady = ppf%m%steadystate(para)
    old_local%particles = 0.0d0

    do j = 1, ppf%m%nvars
       old_local%particles(j,:) = endogsteady(j)
    end do


    !if (rank==0) call ppf%m%describe_params()
    !if (rank==0) print*, old_local%particles(:,1)
    call mpi_barrier(MPI_COMM_WORLD, mpierror)

    !if (rank==0) print*,'entering first initializaztion'
    ! first initialization
    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%nlocalpart)

       !new_local%particles(:,1) = ppf%m%g(old_local%particles(:,ppf%nlocalpart), shocks(:,1), [1,1])
       j = 2
       explode_count = 1
       do while (j <= ppf%nlocalpart)

          !do j = 2, ppf%nlocalpart
          old_local%particles(:,j) = ppf%m%g(old_local%particles(:,j-1), shocks(:,j), regime)
          ! check for explosions
          if ( isnan(sum(old_local%particles(:,j))) ) then
             j = max(j-3, 2)
             shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%nlocalpart)
             explode_count = explode_count + 1

             if (explode_count > 1) then
                print*,'particle',j,'explodes on rank',rank

                old_local%particles(:,j-1) = 0.0d0
                do k = 1, ppf%m%nvars
                   old_local%particles(k,j-1) = endogsteady(k)
                end do

             end if
          else
             j = j + 1
             explode_count = 1
          end if

       end do
    end do


    ! second initialization
    call mpi_barrier(MPI_COMM_WORLD, mpierror)
    !if (rank==0) print*,'entering second initialization'

    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%nlocalpart)

       do j = 1, ppf%nlocalpart
          new_local%particles(:,j) = ppf%m%g(old_local%particles(:,j), shocks(:,j), regime)
       end do
       old_local = new_local
    end do

    do j = 1, ppf%nlocalpart
       do k = 1, ppf%m%nvars
          if (isnan(old_local%particles(k,j))) then
             old_local%particles(:,j) = 0.0d0
             old_local%weights(j) = 0.0d0
          end if
       end do
    end do

    if (sum(old_local%weights) == 0.0d0) then
       print*,'divergence on rank ', rank
       lik = -1000000000.0
       stop
    end if
    old_copy = old_local 

    call mpi_barrier(MPI_COMM_WORLD, mpierror)

    old_copy%weights = 1.0d0 / (ppf%nlocalpart * nproc)
    old_local%weights = 1.0d0 / (ppf%nlocalpart * nproc)

    ! ! the particle filter
    !if (rank==0) print*,'entering loop'
    do t = 1, ppf%m%T

       call mpi_irecv(old_foreign%particles, ppf%m%nvars*ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
            left, 1, MPI_COMM_WORLD, mpi_req_recv1, mpierror)
       call mpi_irecv(old_foreign%weights, ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
            left, 1, MPI_COMM_WORLD, mpi_req_recv2, mpierror)
       call mpi_isend(old_copy%particles(:,1:ppf%nforeignpart), ppf%m%nvars*ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
            right, 1, MPI_COMM_WORLD, mpi_req_send1, mpierror)
       call mpi_isend(old_copy%weights(1:ppf%nforeignpart), ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
            right, 1, MPI_COMM_WORLD, mpi_req_send2, mpierror)

       ! draw shocks and regimes
       shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%nlocalpart)

       if (ppf%adjusted_proposal(t)) then
          is_sampling_std = ppf%adjusted_proposal_std(t)
          is_sampling_mu = ppf%adjusted_proposal_mu(t,:)

          ! should be speed up
          shocks = shocks * is_sampling_std
          do j = 1, ppf%m%nexog
             shocks(j,:) = shocks(j,:) + is_sampling_mu(j)
          end do
       end if

       do j = 1, ppf%nlocalpart - ppf%nforeignpart
          new_local%particles(:,j) = ppf%m%g(old_local%particles(:,ppf%nforeignpart+j), &
               shocks(:,j), [1,1])
       end do

       !calculate unnormalized weights 
       do j = 1, ppf%nlocalpart - ppf%nforeignpart
          incwt = ppf%m%pdfy(t, new_local%particles(:,j), old_local%particles(:,ppf%nforeignpart+j), para)
          new_local%weights(j) = old_local%weights(ppf%nforeignpart+j) * incwt 
       end do

       call mpi_wait(mpi_req_recv1, mpi_recv_status1, mpierror) 

       do j = 1, ppf%nforeignpart
          new_local%particles(:,foreignpartstart+j) = ppf%m%g(old_foreign%particles(:,j), &
               shocks(:,foreignpartstart+j), [1,1])
       end do

       call mpi_wait(mpi_req_recv2, mpi_recv_status2, mpierror)

       do j = 1, ppf%nforeignpart
          incwt = ppf%m%pdfy(t, new_local%particles(:, foreignpartstart+j), &
               old_foreign%particles(:,j), para)
          !if (incwt==0) print*,'***',j,rank,incwt,new_local%particles(:,foreignpartstart+j)
          new_local%weights(foreignpartstart+j) = old_foreign%weights(j) * incwt
       end do

       call mpi_wait(mpi_req_send1, mpi_send_status1,mpierror)
       call mpi_wait(mpi_req_send2, mpi_send_status2,mpierror)

       ! could be speed up
       if (ppf%adjusted_proposal(t)) then

          do j = 1, ppf%nlocalpart
             kap = is_sampling_std ** ppf%m%nexog &
                  * exp(-0.5d0 * dot_product(shocks(:,j),shocks(:,j))) &
                  / exp(-0.5d0 * dot_product(shocks(:,j)-is_sampling_mu, shocks(:,j)-is_sampling_mu) / is_sampling_std ** 2)
             new_local%weights(j) = new_local%weights(j) * kap
          end do
       end if


       minwt = minval(new_local%weights)
       maxwt = maxval(new_local%weights)

       do j = 1,ppf%nlocalpart
          if (isnan(new_local%weights(j))) new_local%weights(j) = 0.000000000001d0
       end do

       ! check the particles
       call new_local%normalize_weights(py)
       ess = new_local%ess()
       call mpi_allreduce(ess, min_ess, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierror)         
       !new_local%weights = new_local%weights * relative_weight 

       ! check lnpy

       ! if (sumproby ==0.0d0) then
       !    write(*,*) 'divergence of particle filter on rank', rank, 'at time ', t
       ! end if

       ! if (sumproby == nlocalparticles*minweight_param) then
       !    write(*,*) 'divergence of particle filter on rank', rank, 'avoided at time ', t
       ! end if

       ! if (isnan(sumproby)==.TRUE.) then
       !    write(*,*) 'prob(y) is nan on rank', rank
       ! end if
       ! totsumproby = 0.0d0
       py_across_procs = 0.0d0
       call mpi_allreduce(py, py_across_procs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       relative_weight = py / py_across_procs

       lik = lik + log(py_across_procs)

       if (min_ess < ppf%nlocalpart) then
          !if (rank==0) print*,'min_ess',min_ess, ppf%nlocalpart, t
          call new_local%systematic_resampling(randu(t,1))
          new_local%weights = relative_weight / ppf%nlocalpart
       end if
       old_local = new_local


       py = sum(old_local%weights)
       py_vector = 0.0d0
       call mpi_allgather(py, 1, MPI_DOUBLE_PRECISION, py_vector, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierror)
       effective_procs = 1.0d0 / sum(py_vector**2)

       !if (rank==0) print*,'before passing', t, effective_procs 
       if (effective_procs < (nproc/2.0d0)) then 
          ! sort processors by relative weight, partner 
          remaining_procs = .true.
          do i = 1, nproc
             ranked_procs(i) = minloc(py_vector, dim=1, mask=remaining_procs) - 1
             remaining_procs(ranked_procs(i)+1) = .false.
          end do
          partner = nproc - minloc(ranked_procs, dim=1, mask=ranked_procs==rank) + 1
          partner = ranked_procs(partner)
          !print*,'left', left, 'partner', partner

          ! pass the last nalt particlse
          call mpi_barrier(MPI_COMM_WORLD, mpierror)

          alt_copy%particles = old_local%particles(:,ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart)
          alt_copy%weights = old_local%weights(ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart)

          call mpi_irecv(alt_foreign%particles, ppf%m%nvars*ppf%naltpart, MPI_DOUBLE_PRECISION, partner, 2, MPI_COMM_WORLD, mpi_req_recv3, mpierror)
          call mpi_irecv(alt_foreign%weights, ppf%naltpart, MPI_DOUBLE_PRECISION, partner, 2, MPI_COMM_WORLD, mpi_req_recv4, mpierror)

          call mpi_isend(alt_copy%particles, ppf%m%nvars*ppf%naltpart, MPI_DOUBLE_PRECISION, partner, 2, MPI_COMM_WORLD, mpi_req_send3, mpierror)
          call mpi_isend(alt_copy%weights, ppf%naltpart, MPI_DOUBLE_PRECISION, partner, 2, MPI_COMM_WORLD, mpi_req_send4, mpierror)

          call mpi_wait(mpi_req_send3, mpi_send_status3, mpierror)
          call mpi_wait(mpi_req_send4, mpi_send_status4, mpierror)
          call mpi_wait(mpi_req_recv3, mpi_recv_status3, mpierror)
          call mpi_wait(mpi_req_recv4, mpi_recv_status4, mpierror)

          old_local%particles(:, ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart) = alt_foreign%particles
          old_local%weights(ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart) = alt_foreign%weights

          py = sum(old_local%weights)
          py_vector = 0.0d0
          call mpi_allgather(py, 1, MPI_DOUBLE_PRECISION, py_vector, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierror)
          effective_procs = 1.0d0 / sum(py_vector**2)
          !if (rank==0) print*, 'after passing', t, effective_procs 

       end if

       old_copy%particles = old_local%particles(:,1:ppf%nforeignpart)
       old_copy%weights = old_local%weights(1:ppf%nforeignpart)

       if (t>=105) then
          py = sum(old_copy%weights * 400.0d0*log(old_copy%particles(9,:)))
          py_vector = 0.0d0
          call mpi_allgather(py, 1, MPI_DOUBLE_PRECISION, py_vector, 1, & 
               MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierror)
          if (rank==0) print*, 'filtered nomr in ', t, sum(py_vector) 
       end if

       if (save_states_opt) then
          write(char_rank,'(I3.3)') rank
          write(char_t, '(I3.3)') t
          out_str = 'test-output/t_'//trim(adjustl(char_t))//'_rank_'//trim(adjustl(char_rank))//'_'
          call old_copy%write(out_str)
       end if

       call mpi_barrier(MPI_COMM_WORLD, mpierror)
       !print*,rank,py_vector
      end do

      ! deallocate
      call old_local%free()
      call old_copy%free()
      call old_foreign%free()
      call new_local%free()
      call alt_foreign%free()
      call alt_copy%free()

    end function lik
end module class_ParallelParticleFilter


