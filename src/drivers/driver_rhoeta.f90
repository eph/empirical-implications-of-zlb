program driver_rhoeta

  !use flap, only : command_line_interface
  use class_model, only: model
  use class_ParallelParticleFilter, only: ParallelParticleFilter


  implicit none

  include 'mpif.h'

  type(model) :: dsge
  type(ParallelParticleFilter) :: ppf

  double precision :: rhoeta(30), lik0
  double precision, allocatable :: p0(:)

  character(len=:), allocatable :: pmsv_file
  integer :: rank, nproc, mpierror, i, j, seed, npart

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

  seed = 1848
  npart = 500000

  dsge = model(.true.)

  ppf = ParallelParticleFilter(dsge, npart=npart, nproc=nproc, rank=rank, seed=rank)
   
  ppf%adjusted_proposal(104) = .true.
  ppf%adjusted_proposal_std(104) = 1.2d0
  ppf%adjusted_proposal_mu(104,1) = 3.0d0
   

  allocate(p0(dsge%npara))
  
  rhoeta = [0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0,69]

  !pmsv_file = '/mq/home/m1eph00/tmp/aer_revision_ed/final_code/final-final/pf_stab_para.txt'
  pmsv_file = '/mq/home/m1eph00/tmp/aer_revision_ed/final_code/final-final/mean.txt'
  open(199, file=pmsv_file, action='read')
  do i = 1, dsge%npara
     read(199,*) p0(i)
  end do
  close(199)

  do i = 1, 30

     p0(29) = rhoeta(i)
     ppf%m%solution%startingguess = .false.
     do j = 1, 10
        

        lik0 = ppf%lik(p0, rank=rank, nproc=nproc)
        ppf%m%solution%startingguess = .true. 
        if (rank .eq. 0) then
           print*,'****',rhoeta(i),'lik = ',lik0
           print*,dsge%pr(p0)
        end if
     end do

  end do
  deallocate(p0)

  call MPI_finalize(mpierror)

end program
