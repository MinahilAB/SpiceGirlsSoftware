!> @file atmosphere_model.f90
!> @brief Main driver program for the simple atmospheric model.
!>
!> This program initializes the simulation environment (MPI/OpenACC, grid, state), 
!> executes the main time integration loop, handles I/O output, and calculates 
!> conservation checks (mass and energy).
program atmosphere_model
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq, nx, nz
  use physical_parameters, only : xlen, zlen
  use parallel_parameters, only : ierr, rank
  use iodir, only : stdout
  use module_nvtx
  use parallel_timer
  use mpi
#if defined(_OACC)
  use openacc
#endif
#if defined(_OMP)
  use omp_lib
#endif

  implicit none

  real(wp) :: etime
  real(wp) :: ptime
  real(wp) :: output_counter
  real(wp) :: pctime
  real(wp) :: mass0, te0
  real(wp) :: mass1, te1
  integer(8) :: t1, t2, rate
  integer :: n_args
  character(len=32) :: arg

  
  ! Argument parsing occuring here
  n_args = command_argument_count()
  if (n_args == 2) then
    call get_command_argument(1, arg)
    read(arg, *) nx
    call get_command_argument(2, arg)
    read(arg, *) sim_time
  end if

  nz = int(nx * zlen/xlen)

  ! Initialisation
  if (rank == 0) write(stdout, '(/,A,/)') 'SIMPLE ATMOSPHERIC MODEL STARTING.'

  ! Use NVTX to mark the main computational region for profiling
  call nvtx_push('tot')

  call init(etime,output_counter,dt)

  call total_mass_energy(mass0,te0)
  call create_output( )
  call write_record(oldstat,ref,etime)

  call system_clock(t1)

  ptime = int(sim_time/10.0)

  call nvtx_push('main_loop')
  
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

  call nvtx_pop()

  
  call total_mass_energy(mass1,te1)
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

  call nvtx_pop()

  call MPI_Finalize(ierr)
  
end program atmosphere_model