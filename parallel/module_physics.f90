module module_physics
  use calculation_types, only : wp
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use legendre_quadrature
  use dimensions
  use iodir
  use module_types
  use module_nvtx
#if defined(_OACC)
  use openacc
#endif

  implicit none

  private

  public :: init
  public :: finalize
  public :: rungekutta
  public :: total_mass_energy

  real(wp), public :: dt
  real(wp) :: dx, dz

  type(atmospheric_state), public :: oldstat
  type(atmospheric_state), public :: newstat
  type(atmospheric_tendency), public :: tend
  type(atmospheric_flux), public :: flux
  type(reference_state), public :: ref

  contains

  subroutine init(etime,output_counter,dt)
    implicit none
    real(wp), intent(out) :: etime, output_counter, dt
    integer :: i, k, ii, kk
    real(wp) :: x, z, r, u, w, t, hr, ht

    call nvtx_push('init')
    

    ! MPI initialisation here
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, csize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    dims(1) = csize
    periods(1) = .true.
    reorder = .true.
    call MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, cart_comm, ierr)
    call MPI_Cart_shift(cart_comm, 0, -1, rank_source_dummy, left_rank, ierr)
    call MPI_Cart_shift(cart_comm, 0, 1, rank_source_dummy, right_rank, ierr)
    
    
    ! Assigning GPU to rank here
#ifdef _OACC
    num_devices = acc_get_num_devices(acc_device_nvidia)
    
    if ( (num_devices == 0) .and. (rank == 0) ) then
      write(*, '(/,A,/)') 'FATAL: No NVIDIA GPUs detected by OpenACC. Aborting...'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    device_id = mod(rank, num_devices)
    call acc_set_device_num(device_id, acc_device_nvidia)
#endif

    ! Grid dimensions get initialised here
    rest = mod(nx, csize)
    nx_loc = nx / csize

    if (rank < rest) then
      nx_loc = nx_loc + 1
      i_beg = rank * nx_loc + 1
    else
      i_beg = rank * nx_loc + rest + 1
    end if
    i_end = i_beg + nx_loc - 1

  
    ! The rest of the initialisation happening now
    dx = xlen / nx
    dz = zlen / nz

    call oldstat%new_state( )
    call newstat%new_state( )
    call flux%new_flux( )
    call tend%new_tendency( )
    call ref%new_ref( )

    dt = min(dx,dz) / max_speed * cfl
    etime = 0.0_wp
    output_counter = 0.0_wp

!     write(stdout,'(/, A, I0, A)')  'R', rank, ': INITIALIZING MODEL STATUS.'
! #if defined(_OACC)
!     write(stdout,'(A, I0)')        'GPU ID     : ', device_id
! #endif
!     write(stdout,'(A, I0)')        'nx         : ', nx
!     write(stdout,'(A, I0)')        'nx_loc     : ', nx_loc
!     write(stdout,'(A, I0)')        'nz         : ', nz
!     write(stdout,'(A, F0.18)')     'dx         : ', dx
!     write(stdout,'(A, F0.18)')     'dz         : ', dz
!     write(stdout,'(A, F0.18)')     'dt         : ', dt
!     write(stdout,'(A, F0.18, /)')  'final time : ', sim_time

    call oldstat%set_state(0.0_wp)

#if defined(_OMP)
    !$omp parallel do collapse(2) private(i,k,ii,kk,x,z,r,u,w,t,hr,ht)
#endif
    do k = 1-hs, nz+hs
      do i = 1-hs, nx_loc+hs
        do kk = 1, nqpoints
          do ii = 1, nqpoints
            x = (i_beg-1 + i-0.5_wp) * dx + (qpoints(ii)-0.5_wp)*dx
            z = (k_beg-1 + k-0.5_wp) * dz + (qpoints(kk)-0.5_wp)*dz
            call thermal(x,z,r,u,w,t,hr,ht)
            oldstat%dens(i,k) = oldstat%dens(i,k) + &
                       r * qweights(ii)*qweights(kk)
            oldstat%umom(i,k) = oldstat%umom(i,k) + &
                      (r+hr)*u * qweights(ii)*qweights(kk)
            oldstat%wmom(i,k) = oldstat%wmom(i,k) + &
                      (r+hr)*w * qweights(ii)*qweights(kk)
            oldstat%rhot(i,k) = oldstat%rhot(i,k) + &
                    ( (r+hr)*(t+ht) - hr*ht ) * qweights(ii)*qweights(kk)
          end do
        end do
      end do
    end do
    newstat = oldstat
    ref%density(:) = 0.0_wp
    ref%denstheta(:) = 0.0_wp
#if defined(_OMP)
    !$omp end parallel do
#endif

#if defined(_OMP)
    !$omp parallel do collapse(2) private(k,kk,z,hr,ht)
#endif
    do k = 1-hs, nz+hs
      do kk = 1, nqpoints
        z = (k_beg-1 + k-0.5_wp) * dz + (qpoints(kk)-0.5_wp)*dz
        call thermal(0.0_wp,z,r,u,w,t,hr,ht)
        ref%density(k) = ref%density(k) + hr * qweights(kk)
        ref%denstheta(k) = ref%denstheta(k) + hr*ht * qweights(kk)
      end do
    end do

    do k = 1, nz+1
      z = (k_beg-1 + k-1) * dz
      call thermal(0.0_wp,z,r,u,w,t,hr,ht)
      ref%idens(k) = hr
      ref%idenstheta(k) = hr*ht
      ref%pressure(k) = c0*(hr*ht)**cdocv
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif

#if defined(_OACC)
    !$acc enter data create(oldstat, newstat, flux, tend, ref)
    !$acc enter data copyin(oldstat%mem, newstat%mem, flux%mem, tend%mem)
    !$acc enter data copyin(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc enter data attach(oldstat%mem, newstat%mem, flux%mem, tend%mem)
    !$acc enter data attach(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)

    
    allocate(send_left_d(hs*nz*NVARS), send_right_d(hs*nz*NVARS), &
      recv_left_d(hs*nz*NVARS), recv_right_d(hs*nz*NVARS))
    !$acc enter data create(send_left_d, send_right_d, recv_left_d, recv_right_d)
#else
    allocate(send_left(hs*nz*NVARS), send_right(hs*nz*NVARS), &
      recv_left(hs*nz*NVARS), recv_right(hs*nz*NVARS))
#endif

    call nvtx_pop()
  end subroutine init



  subroutine rungekutta(s0,s1,fl,tend,dt)
    implicit none

    type(atmospheric_state), intent(inout) :: s0
    type(atmospheric_state), intent(inout) :: s1
    type(atmospheric_flux), intent(inout) :: fl
    type(atmospheric_tendency), intent(inout) :: tend
    real(wp), intent(in) :: dt
    real(wp) :: dt1, dt2, dt3
    logical, save :: dimswitch = .true.

    call nvtx_push('rungekutta')

    dt1 = dt/1.0_wp
    dt2 = dt/2.0_wp
    dt3 = dt/3.0_wp
    if ( dimswitch ) then
      call step(s0, s0, s1, dt3, DIR_X, fl, tend)
      call step(s0, s1, s1, dt2, DIR_X, fl, tend)
      call step(s0, s1, s0, dt1, DIR_X, fl, tend)
      call step(s0, s0, s1, dt3, DIR_Z, fl, tend)
      call step(s0, s1, s1, dt2, DIR_Z, fl, tend)
      call step(s0, s1, s0, dt1, DIR_Z, fl, tend)
    else
      call step(s0, s0, s1, dt3, DIR_Z, fl, tend)
      call step(s0, s1, s1, dt2, DIR_Z, fl, tend)
      call step(s0, s1, s0, dt1, DIR_Z, fl, tend)
      call step(s0, s0, s1, dt3, DIR_X, fl, tend)
      call step(s0, s1, s1, dt2, DIR_X, fl, tend)
      call step(s0, s1, s0, dt1, DIR_X, fl, tend)
    end if
    dimswitch = .not. dimswitch

    call nvtx_pop()
  end subroutine rungekutta



  ! Semi-discretized step in time:
  ! s2 = s0 + dt * rhs(s1)
  subroutine step(s0, s1, s2, dt, dir, fl, tend)
    implicit none
    type(atmospheric_state), intent(in) :: s0
    type(atmospheric_state), intent(inout) :: s1
    type(atmospheric_state), intent(inout) :: s2
    type(atmospheric_flux), intent(inout) :: fl
    type(atmospheric_tendency), intent(inout) :: tend
    real(wp), intent(in) :: dt
    integer, intent(in) :: dir
    if (dir == DIR_X) then
      call tend%xtend(fl,ref,s1,dx,dt)
    else if (dir == DIR_Z) then
      call tend%ztend(fl,ref,s1,dz,dt)
    end if
    call s2%update(s0,tend,dt)
  end subroutine step



  subroutine thermal(x,z,r,u,w,t,hr,ht)
#if defined(_OACC)
    !$acc routine seq
#endif
    implicit none
    real(wp), intent(in) :: x, z
    real(wp), intent(out) :: r, u, w, t
    real(wp), intent(out) :: hr, ht
    call hydrostatic_const_theta(z,hr,ht)
    r = 0.0_wp
    t = 0.0_wp
    u = 0.0_wp
    w = 0.0_wp
    t = t + ellipse(x,z,3.0_wp,hxlen,p1,p1,p1)
  end subroutine thermal



  subroutine hydrostatic_const_theta(z,r,t)
#if defined(_OACC)
    !$acc routine seq
#endif
    implicit none
    real(wp), intent(in) :: z
    real(wp), intent(out) :: r, t
    real(wp) :: p,exner,rt
    t = theta0
    exner = exner0 - grav * z / (cp * theta0)
    p = p0 * exner**(cp/rd)
    rt = (p / c0)**cvocd
    r = rt / t
  end subroutine hydrostatic_const_theta



  elemental function ellipse(x,z,amp,x0,z0,x1,z1) result(val)
#if defined(_OACC)
    !$acc routine seq
#endif
    implicit none
    real(wp), intent(in) :: x, z
    real(wp), intent(in) :: amp
    real(wp), intent(in) :: x0, z0
    real(wp), intent(in) :: x1, z1
    real(wp) :: val
    real(wp) :: dist
    dist = sqrt( ((x-x0)/x1)**2 + ((z-z0)/z1)**2 ) * hpi
    if (dist <= hpi) then
      val = amp * cos(dist)**2
    else
      val = 0.0_wp
    end if
  end function ellipse



  subroutine finalize()
    implicit none

#if defined(_OACC)
    !$acc exit data delete(oldstat%mem, newstat%mem, flux%mem, tend%mem)
    !$acc exit data delete(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc exit data delete(oldstat, newstat, flux, tend, ref)
    !$acc exit data delete(send_left_d, send_right_d, recv_left_d, recv_right_d)

    deallocate(send_left_d, send_right_d, recv_left_d, recv_right_d)
#else
    deallocate(send_left, send_right, recv_left, recv_right)
#endif

    call oldstat%del_state( )
    call newstat%del_state( )
    call flux%del_flux( )
    call tend%del_tendency( )
    call ref%del_ref( )
  end subroutine finalize



  subroutine total_mass_energy(mass,te)
    implicit none
    real(wp), intent(out) :: mass, te
    integer :: i, k, ierr
    real(wp) :: r, u, w, th, p, t, ke, ie
    real(wp) :: mass_loc, te_loc

    mass_loc = 0.0_wp
    te_loc   = 0.0_wp

#if defined(_OACC)
    !$acc parallel loop collapse(2) present(oldstat%mem, ref%density, ref%denstheta) reduction(+:mass_loc, te_loc)
#elif defined(_OMP)
    !$omp parallel do collapse(2) private(i,k,r,u,w,th,p,t,ke,ie) reduction(+:mass_loc, te_loc)
#endif
    do k = 1, nz
      do i = 1, nx_loc
        r  = oldstat%mem(i, k, I_DENS) + ref%density(k)
        u  = oldstat%mem(i, k, I_UMOM) / r
        w  = oldstat%mem(i, k, I_WMOM) / r
        th = (oldstat%mem(i, k, I_RHOT) + ref%denstheta(k)) / r
        p  = c0*(r*th)**cdocv
        t  = th / (p0/p)**rdocp
        ke = r*(u*u + w*w)
        ie = r*cv*t
        mass_loc = mass_loc + r*dx*dz
        te_loc   = te_loc   + (ke + r*cv*t)*dx*dz
      end do
    end do
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif

    call MPI_Allreduce(mass_loc, mass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(te_loc, te, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  end subroutine total_mass_energy

end module module_physics