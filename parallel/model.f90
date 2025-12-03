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
  use parallel_parameters, only : rank, csize
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

  ! --- MPI Initialization and Argument Parsing ---
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, csize, ierr)

  n_args = command_argument_count()
  if (n_args >= 2) then
    call get_command_argument(1, arg)
    read(arg, *) nx
    call get_command_argument(2, arg)
    read(arg, *) sim_time
  end if

  ! Set vertical resolution based on aspect ratio
  nz = int(nx * zlen/xlen)
  
  ! Initialize parallel parameters, data structures, and OpenACC data regions
  call init(etime, output_counter, dt)
  
  if (rank == 0) call create_output()

  ! Calculate Initial Mass and Energy
  call total_mass_energy(mass0,te0)
  mass_buf(1) = mass0
  mass_buf(2) = te0
  
  ! Global Sum
  call MPI_Allreduce(MPI_IN_PLACE, mass_buf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  mass0 = mass_buf(1)
  te0 = mass_buf(2)

  if (rank == 0) then
      write(stdout,*) 'Initial Mass: ', mass0
      write(stdout,*) 'Initial Energy: ', te0
      ! Set progress print time
      ptime = sim_time / 10.0_wp
  end if
  
  ! --- Main Time-Stepping Loop ---
  call nvtx_push("Main Loop") 

  do while (etime < sim_time)

    call nvtx_push("RungeKutta")
    call rungekutta(oldstat, newstat, flux, tend, ref, dt)
    call nvtx_pop()

    ! Output to NetCDF
    if ( output_counter >= output_freq ) then
       ! Move data from device to host for output calculation
       !$acc update self(oldstat%mem)
       if (rank==0) call write_record(oldstat,ref,etime)
       output_counter = output_counter - output_freq
    end if
    
    ! Print progress
    if (rank == 0 .and. (mod(etime,ptime) < dt) ) then
        pctime = (etime/sim_time)*100.0_wp
        write(stdout,'(1x,a,i2,a)') 'TIME PERCENT : ', int(pctime), '%'
    end if

    etime = etime + dt
  end do

  call nvtx_pop()

  ! Calculate Final Mass and Energy
  call total_mass_energy(mass1,te1)
  mass_buf(1) = mass1
  mass_buf(2) = te1
  call MPI_Allreduce(MPI_IN_PLACE, mass_buf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  mass1 = mass_buf(1)
  te1 = mass_buf(2)
  
  ! --- Finalization ---
  if (rank == 0) call close_output()
  call finalize()
  call MPI_Finalize(ierr)

  if (rank == 0) then
    write(stdout,'(/,A)') "----------------- Atmosphere check ----------------"
    write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
    write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
    write(stdout,'(A)') "---------------------------------------------------"
  end if

end program atmosphere_model
