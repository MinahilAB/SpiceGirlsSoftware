program atmosphere_model
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq
  use iodir, only : stdout
  use mpi
  use parallel_timer
  implicit none

  real(wp) :: etime
  real(wp) :: ptime
  real(wp) :: output_counter
  real(wp) :: pctime
  real(wp) :: mass0, te0
  real(wp) :: mass1, te1
  integer(8) :: t1, t2, rate
  ! type(timer_type) :: t
  integer :: rank, size
  integer :: ierror

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

  write(stdout, *) 'SIMPLE ATMOSPHERIC MODEL STARTING.'
  call init(etime,output_counter,dt)
  call total_mass_energy(mass0,te0)
  call create_output( )
  call write_record(oldstat,ref,etime)

  call system_clock(t1)

#if defined(_OACC)
  !$acc data present(oldstat, newstat, flux, tend, ref)
#endif

 ! Use NVTX to mark the main computational region for profiling
  ! call nvtxRangeStartA('Main Time Loop')

  ptime = int(sim_time/10.0)
  do while (etime < sim_time)

    if (etime + dt > sim_time) dt = sim_time - etime
    
  
      call rungekutta(oldstat,newstat,flux,tend,dt)


    if ( mod(etime,ptime) < dt ) then
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

   ! call nvtxRangeEnd()

  ! Exit the OpenACC data region.
#if defined(_OACC)
  !$acc end data
#endif

  ! Ensure all device operations are complete before diagnostics
#if defined(_OACC)
  !$acc wait
#endif

  call total_mass_energy(mass1,te1)
  call close_output( )

  write(stdout,*) "----------------- Atmosphere check ----------------"
  write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
  write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
  write(stdout,*) "---------------------------------------------------"

  call finalize()
  call system_clock(t2,rate)

  write(stdout,*) "SIMPLE ATMOSPHERIC MODEL RUN COMPLETED."
  write(stdout,*) "USED CPU TIME: ", dble(t2-t1)/dble(rate)

  call MPI_FINALIZE(ierror)
end program atmosphere_model
