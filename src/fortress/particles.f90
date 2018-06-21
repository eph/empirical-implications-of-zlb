module class_ParticleSystem

  !use rng_parallel


  implicit none

  type ParticleSystem

     integer :: npart = 1000
     integer :: nvars = 1

     double precision, allocatable :: particles(:,:), weights(:)

   contains

     ! procedure :: mean, variance
     procedure :: mean
     procedure :: normalize_weights
     procedure :: ESS
     procedure :: describe 
     procedure :: systematic_resampling
     procedure :: free
     procedure :: write
     !final :: cleanup 

  end type ParticleSystem


  interface ParticleSystem
     module procedure new_ParticleSystem
  end interface ParticleSystem

contains

  type(ParticleSystem) function new_ParticleSystem(npart, nvars) result(p)

    integer, optional, intent(in) :: npart, nvars

    if (present(npart)) p%npart = npart
    if (present(nvars)) p%nvars = nvars

    allocate(p%particles(p%nvars, p%npart), p%weights(p%npart))

    p%particles = 0.0d0
    p%weights = 1.0d0 / p%npart

  end function new_ParticleSystem


  subroutine normalize_weights(p, z)
    class(ParticleSystem) :: p
    double precision, intent(out) :: z 

    z = sum(p%weights)

    p%weights = p%weights / z 

  end subroutine normalize_weights


  function ESS(p) result(effective_sample_size)
    class(ParticleSystem) :: p

    double precision :: effective_sample_size 

    effective_sample_size = 1.0d0 / sum(p%weights**2)

  end function ESS

  subroutine describe(p) 
    class(ParticleSystem) :: p

    double precision :: mu(p%nvars), std(p%nvars)

    integer :: i

    print*, 'Describing Particle Swarm' 
    print*, '# variables = ', p%nvars
    print*, '# particles = ', p%npart
    print*, 'ESS         = ', p%ESS()

    mu = p%mean()
    do i = 1, p%nvars
       print*, 'variable', i, mu(i)
    end do

  end subroutine describe

  function mean(p) result(mu)
    class(ParticleSystem) :: p
    double precision :: mu(p%nvars)

    integer :: i

    do i = 1,p%nvars
       mu(i) = dot_product(p%particles(i,:), p%weights)
    end do

  end function mean

  subroutine systematic_resampling(p, randu, resample_ind)
    class(ParticleSystem) :: p

    double precision, intent(in) :: randu
    integer, intent(out), optional :: resample_ind(p%npart)
    double precision :: cdf(p%npart), uu(p%npart)
    double precision :: part_rep(p%nvars, p%npart)

    integer :: i,j

    !wtold = wt/sum(wt)
    cdf(1) = p%weights(1)!wtold(1)
    do i=2,p%npart
       cdf(i) = cdf(i-1) + p%weights(i)
    end do

    uu = ( randu -1.0d0 + real( (/ (i, i=1,p%npart) /) ,8) ) / real(p%npart,8)


    j=1
    do i=1,p%npart
       ! move along the CDF
       do while (uu(i)>cdf(j))
          j=j+1
       end do
       ! shuffling
       part_rep(:,i) = p%particles(:,j)

       if (present(resample_ind)) resample_ind(i) = j
    end do

    p%particles = part_rep
    p%weights = 1.0d0 / p%npart

  end subroutine systematic_resampling

  subroutine write(self, file_prefix)
    class(ParticleSystem) :: self

    character(len=200), intent(in) :: file_prefix

    character(len=200) :: state_name, weight_name

    integer :: i

    state_name = trim(adjustl(trim(file_prefix)))//"states.txt"
    weight_name = trim(adjustl(trim(file_prefix)))//"weights.txt"

    open(345, file=state_name, action='write')
    open(346, file=weight_name, action='write')

    do i = 1, self%npart
       write(345, '(100f16.8)') self%particles(:, i)
       write(346, '(100f16.8)') self%weights(i)
    end do

    close(345)
    close(346)

  end subroutine write


  
  subroutine free(p) 
    class(ParticleSystem) :: p

    deallocate(p%particles, p%weights)

  end subroutine free

  subroutine cleanup(p) 
    type(ParticleSystem) :: p

    deallocate(p%particles, p%weights)

  end subroutine cleanup



end module
