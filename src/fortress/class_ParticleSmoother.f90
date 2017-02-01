module class_ParticleSmoother

  use class_model, only: model
  use class_ParticleSystem, only: ParticleSystem
  use class_RandomNumber, only: RandomNumber

  implicit none

  include 'mpif.h'

  type ParticleSmoother

     integer :: npart = 1000
     integer :: nproc = 1

     integer :: initialization_T = 1

     logical, allocatable :: adjusted_proposal(:)
     double precision, allocatable :: adjusted_proposal_mu(:,:), adjusted_proposal_std(:)

     type(Model) :: m

     type(RandomNumber) :: rng



   contains

     !     procedure :: filter
     !     procedure :: smoother
     procedure :: filter_and_smooth

  end type ParticleSmoother

  interface ParticleSmoother
     module procedure new_ParticleSmoother
  end interface ParticleSmoother

contains

  type(ParticleSmoother) function new_ParticleSmoother & 
       (m, npart, seed) result(ppf)

    class(Model) :: m

    integer, optional, intent(in) :: npart, seed

    integer :: rng_seed

    if (present(npart)) ppf%npart = npart

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

  end function new_ParticleSmoother


  subroutine filter_and_smooth(ppf, para, nsave, filtered_states, smoothed_states, filtered_shocks, smoothed_shocks) 
    class(ParticleSmoother) :: ppf

    double precision, intent(in) :: para(ppf%m%npara)

    integer, intent(in) :: nsave

    double precision, intent(out) :: filtered_states(0:ppf%m%T, ppf%m%nvars, nsave) 
    double precision, intent(out) :: smoothed_states(0:ppf%m%T, ppf%m%nvars, nsave)
    double precision, intent(out) :: filtered_shocks(0:ppf%m%T, ppf%m%nexog, nsave) 
    double precision, intent(out) :: smoothed_shocks(0:ppf%m%T, ppf%m%nexog, nsave)

    logical :: converged 


    type(ParticleSystem) :: old_local, new_local

    double precision :: particle_family(0:ppf%m%T, ppf%m%nvars, ppf%npart)
    double precision :: shock_family(0:ppf%m%T, ppf%m%nexog, ppf%npart)

    integer :: resample_ind(ppf%npart)
    double precision :: incwt, py
    double precision :: ess, randu(ppf%m%T,1)


    double precision :: endogsteady(ppf%m%nvars), shocks(ppf%m%nexog, ppf%npart)
    double precision :: lik

    integer :: i, j, k, t, t_resample

    double precision :: is_sampling_std, is_sampling_mu(ppf%m%nexog), kap

    integer :: explode_count

    integer :: regime(2)

    character(len=5) :: rb, char_rank, char_t
    character(len=200) :: out_str

    converged = ppf%m%solve(para, 1, 0)
    !if (rank==0) print*,'entering particle filter'
    lik = 0.0d0


    if (converged==.false.) then
       lik = -1000000000000.0d0
       return 
    end if

    randu = ppf%rng%uniform_rvs(ppf%m%T, 1)

    ! initialize part
    old_local = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%npart)
    new_local = ParticleSystem(nvars=ppf%m%nvars, npart=ppf%npart)

    endogsteady = ppf%m%steadystate(para)
    old_local%particles = 0.0d0

    do j = 1, ppf%m%nvars
       old_local%particles(j,:) = endogsteady(j)
    end do



    !if (rank==0) print*,'entering first initializaztion'
    ! first initialization
    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%npart)

       !new_local%particles(:,1) = ppf%m%g(old_local%particles(:,ppf%npart), shocks(:,1), [1,1])
       j = 2
       explode_count = 1
       do while (j <= ppf%npart)

          !do j = 2, ppf%npart
          old_local%particles(:,j) = ppf%m%g(old_local%particles(:,j-1), shocks(:,j), regime)
          ! check for explosions
          if ( isnan(sum(old_local%particles(:,j))) ) then

             j = max(j-3, 2)
             shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%npart)
             explode_count = explode_count + 1

             if (explode_count > 1) then
                !print*,'particle',j,'explodes on rank',rank

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

    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%npart)

       do j = 1, ppf%npart
          new_local%particles(:,j) = ppf%m%g(old_local%particles(:,j), shocks(:,j), regime)
       end do
       old_local = new_local
    end do

    do j = 1, ppf%npart
       do k = 1, ppf%m%nvars
          if (isnan(old_local%particles(k,j))) then
             old_local%particles(:,j) = 0.0d0
             old_local%weights(j) = 0.0d0
          end if
       end do
    end do

    if (sum(old_local%weights) == 0.0d0) then
       print*,'divergence on rank 0'
       lik = -1000000000.0
       stop
    end if


    old_local%weights = 1.0d0 / (ppf%npart)

    particle_family(0,:,:) = old_local%particles
    shock_family(0,:,:) = 0.0d0

    filtered_states(0,:,:) = old_local%particles(:,1:nsave)


    ! ! the particle filter
    !if (rank==0) print*,'entering loop'
    do t = 1, ppf%m%T
       print*,t

       ! draw shocks and regimes
       shocks = ppf%rng%norm_rvs(ppf%m%nexog, ppf%npart)

       if (ppf%adjusted_proposal(t)) then
          is_sampling_std = ppf%adjusted_proposal_std(t)
          is_sampling_mu = ppf%adjusted_proposal_mu(t,:)

          ! should be speed up
          shocks = shocks * is_sampling_std
          do j = 1, ppf%m%nexog
             shocks(j,:) = shocks(j,:) + is_sampling_mu(j)
          end do
       end if

       do j = 1, ppf%npart
          new_local%particles(:,j) = ppf%m%g(old_local%particles(:,j), &
               shocks(:,j), [1,1])
          incwt = ppf%m%pdfy(t, new_local%particles(:,j), old_local%particles(:,j), para)
          new_local%weights(j) = old_local%weights(j) * incwt 
       end do

       
       ! could be speed up
       if (ppf%adjusted_proposal(t)) then

          do j = 1, ppf%npart
             kap = is_sampling_std ** ppf%m%nexog &
                  * exp(-0.5d0 * dot_product(shocks(:,j),shocks(:,j))) &
                  / exp(-0.5d0 * dot_product(shocks(:,j)-is_sampling_mu, shocks(:,j)-is_sampling_mu) / is_sampling_std ** 2)
             new_local%weights(j) = new_local%weights(j) * kap
          end do
       end if
       
       do j = 1,ppf%npart
          if (isnan(new_local%weights(j))) then
              new_local%weights(j) = 0.000000000001d0
           end if
       end do
       particle_family(t, :, :) = new_local%particles
       shock_family(t,:,:) = shocks
       ! check the particles
       call new_local%normalize_weights(py)
       ess = new_local%ess()
       print*,log(py),ess
       lik = lik + log(py)

       if (ess < ppf%npart) then
          call new_local%systematic_resampling(randu(t,1), resample_ind)
          new_local%weights = 1.0d0 / ppf%npart

          do t_resample = 0, t

             particle_family(t_resample, :, :) = particle_family(t_resample, :, resample_ind)
             shock_family(t_resample, :, :) = shock_family(t_resample, :, resample_ind)
             

          end do
       end if
       old_local = new_local




       if (t >= 105) then
          print*, t, sum(400.0d0*log(new_local%particles(5,:)))/ppf%npart
       end if


       filtered_states(t, :, :) = new_local%particles(:,1:nsave)
       filtered_shocks(t, :, :) = shock_family(t, :, 1:nsave)



       end do
       ! do smoothing
       do t = 0, ppf%m%T
          smoothed_states(t,:,:) = particle_family(t,:,1:nsave)
          smoothed_shocks(t,:,:) = shock_family(t,:,1:nsave)

          if (t >= 105) then
             print*, t, sum(400.0d0*log(particle_family(t,5,:)))/ppf%npart
          end if

       end do




      ! deallocate
      call old_local%free()
      call new_local%free()
      print*,'lik=',lik
      
    end subroutine filter_and_smooth
end module class_ParticleSmoother


 
