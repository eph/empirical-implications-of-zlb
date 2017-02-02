program driver_altsim

  use class_model, only: model

  use json_module

  use flap

  implicit none

  include 'mpif.h'

  type(command_line_interface) :: cli
  type(Model) :: dsge
  
  character(len=:), allocatable :: in_file, out_file, sim
  double precision, allocatable :: states(:,:), shocks(:,:), p0(:)
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


  zlb = .true.

  out_file = 'test.json'

  call cli%init(progname = 'driver_altsim', &
       version='0.0.1', &
       authors='Chris Gust & Ed Herbst', &
       description='Program for running counterfactuals')

  call cli%add(switch='--infile',required=.true.,def='test.json',help='The JSON file containing a set of filtered/smoothed estimates.')
  call cli%add(switch='--sim',required=.true.,def='nozlb_2008', &
       choices='nozlb_2003,nozlb_2008,great_recession_liq,great_recession_inv,great_recession_tech,great_recession_demand,great_recession_mp', &
       help='The alternative simulation to run.')
  call cli%parse(error=error)
  call cli%get(switch='--infile',val=in_file)
  call cli%get(switch='--sim',val=sim)

  ! do i = 1, command_argument_count()
  !    call get_command_argument(i, arg)

  !    select case(arg)
  !    case('--infile')
  !       call get_command_argument(i+1,arg)
  !       in_file = arg
  !    case('--sim')
  !       call get_command_argument(i+1,arg)
  !       sim = arg
  !    end select
  ! end do

  i0 = index(in_file,'.json')


  if ( (sim=='nozlb_2008') .or. (sim=='nozlb_2003')) then
     zlb = .false.
  end if

  dsge = model(zlb)

  allocate(p0(dsge%npara))
  
  allocate(states(0:dsge%T, dsge%nvars), states_new(0:dsge%T, dsge%nvars))
  allocate(shocks(0:dsge%T, dsge%nexog))
  allocate(tmp(0:dsge%T), shocks_t(dsge%nexog))

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


        states_new(0,:) = states(0,:)
        do t = 1,dsge%T

           shocks_t = shocks(t,:)

           select case(sim)
           case ('great_recession_liq')
              if (t >= 101) then
                 shocks_t(2:dsge%nexog) = 0.0d0
              end if
           case ('great_recession_inv')
              if (t >= 101) then
                 shocks_t(1) = 0.0d0
                 shocks_t(3:dsge%nexog) = 0.0d0
              end if
           case ('great_recession_tech')
              if (t >= 101) then
                 shocks_t(1:2) = 0.0d0
                 shocks_t(4:dsge%nexog) = 0.0d0
              end if
           case ('great_recession_demand')
              if (t >= 101) then
                 shocks_t(1:3) = 0.0d0
                 shocks_t(5:dsge%nexog) = 0.0d0
              end if
           case ('great_recession_mp')
              if (t >= 101) then
                 shocks_t(1:4) = 0.0d0
                 shocks_t(6:dsge%nexog) = 0.0d0
              end if
           case('nozlb_2008')
              if (t==101) states_new(t-1,:) = states(t-1,:)
           case('nozlb_2003')
              if (t==81) states_new(t-1,:) = states(t-1,:)
           end select
           states_new(t,:) = dsge%g(states_new(t-1,:), shocks_t, [1,1])
           print*,states_new(t,1), states(t,1)
           stop

        end do


        call json%create_object(sim_i_s, 'smoothed_states')
        call json%add(sim_i, sim_i_s) 
        do j = 1, dsge%nvars
           write(charvari, '(I2.2)') j 
           call json%add(sim_i_s,'endogvar_'//trim(charvari), states_new(:,j))
        end do
        nullify(sim_i_s)
        nullify(sim_i)

     end do


     nullify(output)
     call json%print(p, 'final-final/alt-sims/'//trim(sim)//trim(in_file(i0-10:)))
     call json%destroy(p)

  end if

  call dsge%cleanup()

  deallocate(p0, states, shocks,tmp, states_new, shocks_t)
  call MPI_finalize(mpierror)

end program

