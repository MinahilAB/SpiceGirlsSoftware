program atmosphere_model
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq, nx, nz
  use physical_parameters, only : xlen, zlen
  use iodir, only : stdout
  use module_nvtx
  use mpi
  use parallel_parameters, only : rank, csize, left_rank, right_rank, cart_comm
  implicit none

  real(wp) :: etime
  real(wp) :: ptime
  real(wp) :: output_counter
  real(wp) :: pctime
  real(wp) :: mass0, te0
  real(wp) :: mass1, te1
  real(wp) :: mass_buf(2)
  integer(8) :: t1, t2, rate
  integer :: ierr
  integer :: n_args
  character(len=32) :: arg
  integer, dimension(1) :: dims
  logical, dimension(1) :: periods
  logical :: reorder



  n_args = command_argument_count()
  if (n_args == 2) then
    call get_command_argument(1, arg)
    read(arg, *) nx
    call get_command_argument(2, arg)
    read(arg, *) sim_time
  end if

  nz = int(nx * zlen/xlen)

  ! MPI initialisation
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, csize, ierr)

  dims(1) = csize
  periods(1) = .true.
  reorder = .true.
  call MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, cart_comm, ierr)

  call MPI_Cart_shift(cart_comm, 0, -1, ierr, left_rank, ierr)
  call MPI_Cart_shift(cart_comm, 0, 1, ierr, right_rank, ierr)
  
  if (rank == 0) write(stdout, '(/,A,/)') 'SIMPLE ATMOSPHERIC MODEL STARTING.'

  call init(etime,output_counter,dt)
  call total_mass_energy(mass0,te0)
  mass_buf(1) = mass0
  mass_buf(2) = te0
  call MPI_Allreduce(MPI_IN_PLACE, mass_buf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  mass0 = mass_buf(1)
  te0 = mass_buf(2)
  call create_output( )
  call write_record(oldstat,ref,etime)

  call system_clock(t1)

  ! Use NVTX to mark the main computational region for profiling
  call nvtx_push('main_loop')

#if defined(_OACC)
  !$acc data present(oldstat, newstat, flux, tend, ref)
#endif

  ptime = int(sim_time/10.0)
  do while (etime < sim_time)

    if (etime + dt > sim_time) dt = sim_time - etime
      call rungekutta(oldstat,newstat,flux,tend,dt)
    
    if ( (rank == 0) .and. (mod(etime,ptime) < dt) ) then
      pctime = (etime/sim_time)*100.0_wp
      write(stdout,'(1x,a,i2,a)') 'TIME PERCENT : ', int(pctime), '%'
    end if

    etime = etime + dt
    output_counter = output_counter + dt

    if (output_counter >= output_freq) then
      output_counter = output_counter - output_freq
      call write_record(oldstat,ref,etime)
    end if

  end do

  ! Exit the OpenACC data region.
#if defined(_OACC)
  !$acc end data
#endif

  ! Ensure all device operations are complete before diagnostics
#if defined(_OACC)
  !$acc wait
#endif

  call nvtx_pop()

  call total_mass_energy(mass1,te1)
  mass_buf(1) = mass1
  mass_buf(2) = te1
  call MPI_Allreduce(MPI_IN_PLACE, mass_buf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  mass1 = mass_buf(1)
  te1 = mass_buf(2)
  call close_output( )

  if (rank == 0) then
    write(stdout,'(/,A)') "----------------- Atmosphere check ----------------"
    write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
    write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
    write(stdout,'(A,/)') "---------------------------------------------------"
  end if

  call finalize()
  
  call system_clock(t2,rate)

  if (rank == 0) then
    write(stdout,'(A)') "SIMPLE ATMOSPHERIC MODEL RUN COMPLETED."
    write(stdout,'(A, F0.18, /)') "USED CPU TIME: ", dble(t2-t1)/dble(rate)
  end if

  call MPI_Finalize(ierr)
  
end program atmosphere_model
