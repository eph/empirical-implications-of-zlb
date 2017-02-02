program ghlss_smoother_driver

  use class_model, only: model
  use model_details, only: calc_premium

  use json_module

  use flap

  implicit none

  include 'mpif.h'

  type(command_line_interface) :: cli
  type(Model) :: dsge
  
  character(len=:), allocatable :: in_file, out_file, sim
  double precision, allocatable :: states(:,:), shocks(:,:), p0(:), equity_premium(:)
  double precision, allocatable :: tmp(:), states_new(:,:), shocks_t(:)


  integer :: nsave
  character(len=500) :: arg, charsimi, charvari, charstate
  integer :: i, n,j, t, i0, error

  type(json_core) :: json 
  type(json_value), pointer :: p, inp, output, sim_i, sim_i_s
  type(json_file) :: input_json

  integer :: rank, nproc, mpierror

  logical :: zlb, found, converged


  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)


  call cli%init(progname = 'driver_smoother', &
       version='0.0.1', &
       authors='Chris Gust & Ed Herbst & Matt Smith', &
       description='Computes the equity premium given a set of filtered estimates.')

  call cli%add(switch='--infile',required=.true.,help='The JSON file containing a set of filtered/smoothed estimates.')o
  call cli%add(switch='--sim',required=.true.,help='The alt simulation name')
  call cli%parse(error=error)
  call cli%get(switch='--infile',val=in_file)
  call cli%get(switch='--sim',val=sim)
  zlb = .true.

  i0 = index(in_file,'.json')


  dsge = model(zlb)

  allocate(p0(dsge%npara))
  
  allocate(states(0:dsge%T, dsge%nvars), states_new(0:dsge%T, dsge%nvars))
  allocate(shocks(0:dsge%T, dsge%nexog))
  allocate(tmp(0:dsge%T), shocks_t(dsge%nexog), equity_premium(0:dsge%T))

  call input_json%initialize()
  call input_json%load_file(in_file);
  !if (json_failed()) stop
  
  call input_json%get('input.p0', p0, found)

  converged = dsge%solve(p0, nproc, rank)
  call mpi_barrier(MPI_COMM_WORLD, mpierror)
  
  if (rank==0) then 
     print*,'Model solved'
     call json%create_object(p,'')
     call json%create_object(inp,'input')
     call json%add(p, inp)
     call json%add(inp, 'p0', p0)
     nullify(inp)

     call json%create_object(output,'output')
     call json%add(p, output)
     do i = 1, 10
        print*,'Running Simulation ', i, 'of 10'

        write(charsimi, '(I3.3)') i 
        call json%create_object(sim_i, 'sim_'//trim(charsimi))
        call json%add(output, sim_i)

        do j = 1,dsge%nvars
           write(charvari, '(I2.2)') j 

           charstate = 'output.sim_'//trim(charsimi)//'.smoothed_states.endogvar_'//trim(charvari)
           call input_json%get(charstate, tmp, found)
           states(:,j) = tmp

        end do

        do j = 1,dsge%nexog
           write(charvari, '(I2.2)') j 

           charstate = 'output.sim_'//trim(charsimi)//'.smoothed_shocks.exogvar_'//trim(charvari)
           call input_json%get(trim(charstate), tmp, found)
           shocks(:,j) = tmp

        end do


        equity_premium(0) = 0.0d0 !states_new(0,:) = states(0,:)
        do t = 0,dsge%T

           shocks_t = shocks(t,:)

           call calc_premium(states(t,:), p0, dsge%solution%poly, dsge%solution%alphacoeff, & 
                equity_premium(t)) 

        end do


        call json%create_object(sim_i_s, 'smoothed_states')
        call json%add(sim_i, sim_i_s) 
        call json%add(sim_i_s,'endogvar_29', equity_premium)
        nullify(sim_i_s)
        nullify(sim_i)

     end do


     nullify(output)
     call json%print(p, 'final-final/alt-sims/'//trim(sim)//trim(in_file(i0-10:)))
     call json%destroy(p)

  end if

  call dsge%cleanup()

  deallocate(p0, states, shocks,tmp, states_new, shocks_t, equity_premium)
  call MPI_finalize(mpierror)

end program

