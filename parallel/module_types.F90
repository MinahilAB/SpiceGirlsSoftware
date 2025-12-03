module module_types
  use calculation_types
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use legendre_quadrature
  use dimensions
  use iodir
  use module_nvtx

  implicit none

  private

  public :: reference_state
  public :: atmospheric_state
  public :: atmospheric_flux
  public :: atmospheric_tendency

  public :: assignment(=)

  type reference_state
    real(wp), allocatable, dimension(:) :: density
    real(wp), allocatable, dimension(:) :: denstheta
    real(wp), allocatable, dimension(:) :: idens
    real(wp), allocatable, dimension(:) :: idenstheta
    real(wp), allocatable, dimension(:) :: pressure
    contains
    procedure, public :: new_ref
    procedure, public :: del_ref
  end type reference_state



  type atmospheric_state
    real(wp), pointer, dimension(:,:,:) :: mem => null( )
    real(wp), pointer, dimension(:,:) :: dens
    real(wp), pointer, dimension(:,:) :: umom
    real(wp), pointer, dimension(:,:) :: wmom
    real(wp), pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_state
    procedure, public :: set_state
    procedure, public :: del_state
    procedure, public :: update
    procedure, public :: exchange_halo_x
    procedure, public :: exchange_halo_z
  end type atmospheric_state



  type atmospheric_flux
    real(wp), pointer, dimension(:,:,:) :: mem => null( )
    real(wp), pointer, dimension(:,:) :: dens
    real(wp), pointer, dimension(:,:) :: umom
    real(wp), pointer, dimension(:,:) :: wmom
    real(wp), pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_flux
    procedure, public :: set_flux
    procedure, public :: del_flux
  end type atmospheric_flux



  type atmospheric_tendency
    real(wp), pointer, dimension(:,:,:) :: mem => null( )
    real(wp), pointer, dimension(:,:) :: dens
    real(wp), pointer, dimension(:,:) :: umom
    real(wp), pointer, dimension(:,:) :: wmom
    real(wp), pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_tendency
    procedure, public :: set_tendency
    procedure, public :: del_tendency
    procedure, public :: xtend
    procedure, public :: ztend
  end type atmospheric_tendency



  interface assignment(=)
    module procedure state_equal_to_state
  end interface assignment(=)

  contains

  subroutine new_state(atmo)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    allocate(atmo%mem(1-hs:nx_loc+hs, 1-hs:nz+hs, NVARS))
    atmo%dens(1-hs:,1-hs:) => atmo%mem(:,:,I_DENS)
    atmo%umom(1-hs:,1-hs:) => atmo%mem(:,:,I_UMOM)
    atmo%wmom(1-hs:,1-hs:) => atmo%mem(:,:,I_WMOM)
    atmo%rhot(1-hs:,1-hs:) => atmo%mem(:,:,I_RHOT)

#if defined(_OACC)
    !$acc enter data copyin(atmo)
    !$acc enter data create(atmo%mem)
    !$acc enter data attach(atmo%dens, atmo%umom, atmo%wmom, atmo%rhot)
#endif

  end subroutine new_state



  subroutine set_state(atmo, xval)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    real(wp), intent(in) :: xval
    if ( .not. associated(atmo%mem) ) then
      write(stderr,*) 'NOT ALLOCATED STATE ERROR AT LINE ', __LINE__
      stop
    end if
#if defined(_OACC)
    !$acc kernels present(atmo)
#endif
    atmo%mem(1-hs:nx_loc+hs, 1-hs:nz+hs, :) = xval
#if defined(_OACC)
    !$acc end kernels
#endif
  end subroutine set_state



  subroutine del_state(atmo)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    nullify(atmo%dens)
    nullify(atmo%umom)
    nullify(atmo%wmom)
    nullify(atmo%rhot)
  end subroutine del_state



  subroutine update(s2,s0,tend,dt)
    implicit none
    class(atmospheric_state), intent(inout) :: s2
    class(atmospheric_state), intent(in) :: s0
    class(atmospheric_tendency), intent(in) :: tend
    real(wp), intent(in) :: dt
    integer :: ll, k, i

#if defined(_OACC)
    !$acc parallel loop collapse(3) present(s2, s0, tend)
#elif defined(_OMP)
    !$omp parallel do collapse(3) default(none) shared(s2, s0, tend, nx_loc, nz, dt) private(ll, k, i)
#endif
    do ll = 1, NVARS
      do k = 1, nz
        do i = 1, nx_loc
          s2%mem(i,k,ll) = s0%mem(i,k,ll) + dt * tend%mem(i,k,ll)
        end do
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine update



  subroutine xtend(tendency,flux,ref,atmostat,dx,dt)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tendency
    class(atmospheric_flux), intent(inout) :: flux
    class(reference_state), intent(in) :: ref
    class(atmospheric_state), intent(inout) :: atmostat
    real(wp), intent(in) :: dx, dt
    integer :: i, k, ll, s
    real(wp) :: r, u, w, t, p, hv_coef
    real(wp), dimension(STEN_SIZE) :: stencil
    real(wp), dimension(NVARS) :: d3_vals, vals

    ! Halo exchange must complete before flux calculation
    call atmostat%exchange_halo_x( )

    hv_coef = -hv_beta * dx / (16.0_wp*dt)
#if defined(_OACC)
    !$acc parallel loop collapse(2) present(flux, ref, atmostat) &
    !$acc& private(i, k, ll, s, r, u, w, t, p, stencil, d3_vals, vals)
#elif defined(_OMP)
    !$omp parallel do collapse(2) default(none) shared(flux, ref, atmostat, dx, dt, hv_coef, nx_loc, nz) &
    !$omp& private(i, k, ll, s, r, u, w, t, p, stencil, d3_vals, vals)
#endif
    do k = 1, nz
      do i = 1, nx_loc+1
        do ll = 1, NVARS
          do s = 1, STEN_SIZE
            stencil(s) = atmostat%mem(i-hs-1+s,k,ll)
          end do
          vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &
                     + 7.0_wp * stencil(2)/12.0_wp &
                     + 7.0_wp * stencil(3)/12.0_wp &
                     - 1.0_wp * stencil(4)/12.0_wp
          d3_vals(ll) = - 1.0_wp * stencil(1) &
                        + 3.0_wp * stencil(2) &
                        - 3.0_wp * stencil(3) &
                        + 1.0_wp * stencil(4)
        end do
        r = vals(I_DENS) + ref%density(k)
        u = vals(I_UMOM) / r
        w = vals(I_WMOM) / r
        t = ( vals(I_RHOT) + ref%denstheta(k) ) / r
        p = c0*(r*t)**cdocv
        flux%dens(i,k) = r*u - hv_coef*d3_vals(I_DENS)
        flux%umom(i,k) = r*u*u+p - hv_coef*d3_vals(I_UMOM)
        flux%wmom(i,k) = r*u*w - hv_coef*d3_vals(I_WMOM)
        flux%rhot(i,k) = r*u*t - hv_coef*d3_vals(I_RHOT)
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif

    ! Second loop for tendency calculation
#if defined(_OACC)
    !$acc parallel loop collapse(3) present(tendency, flux)
#elif defined(_OMP)
    !$omp parallel do collapse(3) default(none) shared(tendency, flux, dx, nx_loc, nz) private(ll, k, i)
#endif
    do ll = 1, NVARS
      do k = 1, nz
        do i = 1, nx_loc
          tendency%mem(i,k,ll) = &
              -( flux%mem(i+1,k,ll) - flux%mem(i,k,ll) ) / dx
        end do
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine xtend



  subroutine ztend(tendency,flux,ref,atmostat,dz,dt)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tendency
    class(atmospheric_flux), intent(inout) :: flux
    class(reference_state), intent(in) :: ref
    class(atmospheric_state), intent(inout) :: atmostat
    real(wp), intent(in) :: dz, dt
    integer :: i, k, ll, s
    real(wp) :: r, u, w, t, p, hv_coef
    real(wp), dimension(STEN_SIZE) :: stencil
    real(wp), dimension(NVARS) :: d3_vals, vals

    ! Halo exchange must complete before flux calculation
    call atmostat%exchange_halo_z(ref)

    hv_coef = -hv_beta * dz / (16.0_wp*dt)
#if defined(_OACC)
    !$acc parallel loop collapse(2) present(flux, ref, atmostat) &
    !$acc& private(i, k, ll, s, r, u, w, t, p, stencil, d3_vals, vals)
#elif defined(_OMP)
    !$omp parallel do collapse(2) default(none) shared(flux, ref, atmostat, dz, dt, hv_coef, nx_loc, nz) &
    !$omp& private(i, k, ll, s, r, u, w, t, p, stencil, d3_vals, vals)
#endif
    do k = 1, nz+1
      do i = 1, nx_loc
        do ll = 1, NVARS
          do s = 1, STEN_SIZE
            stencil(s) = atmostat%mem(i,k-hs-1+s,ll)
          end do
          vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &
                     + 7.0_wp * stencil(2)/12.0_wp &
                     + 7.0_wp * stencil(3)/12.0_wp &
                     - 1.0_wp * stencil(4)/12.0_wp
          d3_vals(ll) = - 1.0_wp * stencil(1) &
                        + 3.0_wp * stencil(2) &
                        - 3.0_wp * stencil(3) &
                        + 1.0_wp * stencil(4)
        end do
        r = vals(I_DENS) + ref%idens(k)
        u = vals(I_UMOM) / r
        w = vals(I_WMOM) / r
        t = ( vals(I_RHOT) + ref%idenstheta(k) ) / r
        p = c0*(r*t)**cdocv - ref%pressure(k)
        if (k == 1 .or. k == nz+1) then
          w = 0.0_wp
          d3_vals(I_DENS) = 0.0_wp
        end if
        flux%dens(i,k) = r*w - hv_coef*d3_vals(I_DENS)
        flux%umom(i,k) = r*w*u - hv_coef*d3_vals(I_UMOM)
        flux%wmom(i,k) = r*w*w+p - hv_coef*d3_vals(I_WMOM)
        flux%rhot(i,k) = r*w*t - hv_coef*d3_vals(I_RHOT)
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif

    ! Second loop for tendency calculation
#if defined(_OACC)
    !$acc parallel loop collapse(3) present(tendency, flux, atmostat)
#elif defined(_OMP)
    !$omp parallel do collapse(3) default(none) shared(tendency, flux, atmostat, dz, nx_loc, nz) private(ll, k, i)
#endif
    do ll = 1, NVARS
      do k = 1, nz
        do i = 1, nx_loc
          tendency%mem(i,k,ll) = &
              -( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
          ! Source term application (gravity)
          if (ll == I_WMOM) then
            tendency%wmom(i,k) = tendency%wmom(i,k) - atmostat%dens(i,k)*grav
          end if
        end do
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine ztend



  subroutine exchange_halo_x(s)
    use mpi
    use dimensions, only : nz, nx_loc
    class(atmospheric_state), intent(inout) :: s
    integer :: ierr, ncount
    integer :: reqs(4)
    real(wp), allocatable :: send_left(:), send_right(:), recv_left(:), recv_right(:)


    call nvtx_push('exchange_halo')

#if defined(_OACC)
    !$acc update self(s)
#endif

    ncount = hs * nz * NVARS
    allocate(send_left(ncount), send_right(ncount), recv_left(ncount), recv_right(ncount))
  
    call pack_strip(s%mem, 1, 1, nx_loc, hs, send_left)
    call pack_strip(s%mem, nx_loc-hs+1, 1, nx_loc, hs, send_right)
    
    call MPI_Irecv(recv_left, ncount, MPI_DOUBLE_PRECISION, left_rank, 102, cart_comm, reqs(1), ierr)
    call MPI_Irecv(recv_right, ncount, MPI_DOUBLE_PRECISION, right_rank, 101, cart_comm, reqs(2), ierr)
    call MPI_Isend(send_right, ncount, MPI_DOUBLE_PRECISION, right_rank, 102, cart_comm, reqs(3), ierr)
    call MPI_Isend(send_left, ncount, MPI_DOUBLE_PRECISION, left_rank, 101, cart_comm, reqs(4), ierr)
    
    call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
    
    call unpack_strip(recv_left, s%mem, 1-hs, 1)
    call unpack_strip(recv_right, s%mem, nx_loc+1, 1)

#if defined(_OACC)
    !$acc update device(s)
#endif
    
    deallocate(send_left, send_right, recv_left, recv_right)

    call nvtx_pop()
  end subroutine exchange_halo_x



  subroutine pack_strip(mem, i_start, k_start, nx_loc, width, buffer)
    real(wp), intent(in) :: mem(1-hs:nx_loc+hs,1-hs:nz+hs,NVARS)
    integer, intent(in) :: i_start, k_start, nx_loc, width
    real(wp), intent(out) :: buffer(:)
    integer :: ll, k, i, pos
    pos = 0

    do ll = 1, NVARS
      do k = k_start, k_start+nz-1
        do i = i_start, i_start+width-1
          pos = pos + 1
          buffer(pos) = mem(i,k,ll)
        end do
      end do
    end do
  end subroutine pack_strip



  subroutine unpack_strip(buffer, mem, i_start, k_start)
    real(wp), intent(in) :: buffer(:)
    real(wp), intent(inout) :: mem(1-hs:nx_loc+hs,1-hs:nz+hs,NVARS)
    integer, intent(in) :: i_start, k_start
    integer :: ll, k, i, pos
    pos = 0

    do ll = 1, NVARS
      do k = k_start, k_start+nz-1
        do i = i_start, i_start+hs-1
          pos = pos + 1
          mem(i,k,ll) = buffer(pos)
        end do
      end do
    end do
  end subroutine unpack_strip




  subroutine exchange_halo_z(s,ref)
    implicit none
    class(atmospheric_state), intent(inout) :: s
    class(reference_state), intent(in) :: ref
    integer :: i, ll
#if defined(_OACC)
    !$acc parallel loop collapse(2) present(s, ref)
#elif defined(_OMP)
    !$omp parallel do private(ll, i) shared(s, ref)
#endif
    do ll = 1, NVARS
      do i = 1-hs,nx_loc+hs
        if (ll == I_WMOM) then
          s%mem(i,-1,ll) = 0.0_wp
          s%mem(i,0,ll) = 0.0_wp
          s%mem(i,nz+1,ll) = 0.0_wp
          s%mem(i,nz+2,ll) = 0.0_wp
        else if (ll == I_UMOM) then
          s%mem(i,-1,ll)   = s%mem(i,1,ll) /  &
              ref%density(1) * ref%density(-1)
          s%mem(i,0,ll)    = s%mem(i,1,ll) /  &
              ref%density(1) * ref%density(0)
          s%mem(i,nz+1,ll) = s%mem(i,nz,ll) / &
              ref%density(nz) * ref%density(nz+1)
          s%mem(i,nz+2,ll) = s%mem(i,nz,ll) / &
              ref%density(nz) * ref%density(nz+2)
        else
          s%mem(i,-1,ll) = s%mem(i,1,ll)
          s%mem(i,0,ll) = s%mem(i,1,ll)
          s%mem(i,nz+1,ll) = s%mem(i,nz,ll)
          s%mem(i,nz+2,ll) = s%mem(i,nz,ll)
        end if
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine exchange_halo_z



  subroutine new_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
    allocate(ref%density(1-hs:nz+hs))
    allocate(ref%denstheta(1-hs:nz+hs))
    allocate(ref%idens(nz+1))
    allocate(ref%idenstheta(nz+1))
    allocate(ref%pressure(nz+1))

#if defined(_OACC)
  !$acc data copyin(ref)
  !$acc enter data create(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
#endif
  end subroutine new_ref


  subroutine del_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
#if defined(_OACC)
    !$acc exit data delete(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc exit data delete(ref)
#endif
    deallocate(ref%density)
    deallocate(ref%denstheta)
    deallocate(ref%idens)
    deallocate(ref%idenstheta)
    deallocate(ref%pressure)
  end subroutine del_ref



  subroutine new_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) deallocate(flux%mem)
    allocate(flux%mem(1:nx_loc+1, 1:nz+1,NVARS))
    flux%dens => flux%mem(:,:,I_DENS)
    flux%umom => flux%mem(:,:,I_UMOM)
    flux%wmom => flux%mem(:,:,I_WMOM)
    flux%rhot => flux%mem(:,:,I_RHOT)

#if defined(_OACC)
    !$acc enter data copyin(flux)
    !$acc enter data create(flux%mem)
    !$acc enter data attach(flux%dens, flux%umom, flux%wmom, flux%rhot)
#endif
  end subroutine new_flux



  subroutine set_flux(flux, xval)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    real(wp), intent(in) :: xval
    if ( .not. associated(flux%mem) ) then
      write(stderr,*) 'NOT ALLOCATED FLUX ERROR AT LINE ', __LINE__
      stop
    end if
#if defined(_OACC)
    !$acc kernels present(flux)
! #elif defined(_OMP)
  ! !$omp parallel do collapse(3) default(none) shared(flux, xval, nx, nz, NVARS)
#endif
  flux%mem(:,:,:) = xval
#if defined(_OACC)
  !$acc end kernels
#endif
  end subroutine set_flux



  subroutine del_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
#if defined(_OACC)
    !$acc exit data delete(flux%mem)
    !$acc exit data delete(flux)
#endif
    if ( associated(flux%mem) ) deallocate(flux%mem)
    nullify(flux%dens)
    nullify(flux%umom)
    nullify(flux%wmom)
    nullify(flux%rhot)
  end subroutine del_flux



  subroutine new_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) deallocate(tend%mem)
    allocate(tend%mem(nx_loc, nz,NVARS))
    tend%dens => tend%mem(:,:,I_DENS)
    tend%umom => tend%mem(:,:,I_UMOM)
    tend%wmom => tend%mem(:,:,I_WMOM)
    tend%rhot => tend%mem(:,:,I_RHOT)
#if defined(_OACC)
    !$acc enter data copyin(tend)
    !$acc enter data create(tend%mem)
    !$acc enter data attach(tend%dens, tend%umom, tend%wmom, tend%rhot)
#endif
  end subroutine new_tendency



  subroutine set_tendency(tend, xval)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    real(wp), intent(in) :: xval
    if ( .not. associated(tend%mem) ) then
      write(stderr,*) 'NOT ALLOCATED FLUX ERROR AT LINE ', __LINE__
      stop
    end if
#if defined(_OACC)
    !$acc kernels present(tend)
! #elif defined(_OMP)
  ! !$omp parallel do collapse(3) default(none) shared(tend, xval, nx, nz, NVARS)
#endif
    tend%mem(:,:,:) = xval
#if defined(_OACC)
    !$acc end kernels
#endif
  end subroutine set_tendency



  subroutine del_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
#if defined(_OACC)
    !$acc exit data delete(tend%mem)
    !$acc exit data delete(tend)
#endif
    if ( associated(tend%mem) ) deallocate(tend%mem)
    nullify(tend%dens)
    nullify(tend%umom)
    nullify(tend%wmom)
    nullify(tend%rhot)
  end subroutine del_tendency



  subroutine state_equal_to_state(x,y)
    implicit none
    type(atmospheric_state), intent(inout) :: x
    type(atmospheric_state), intent(in) :: y
#if defined(_OACC)
    !$acc kernels present(x, y)
    ! #elif defined(_OMP)
    ! !$omp parallel do collapse(3) default(none) shared(x, y, nx, nz, NVARS, hs)
#endif
    x%mem(1-hs:nx_loc+hs, 1-hs:nz+hs, :) = y%mem(1-hs:nx_loc+hs, 1-hs:nz+hs, :)
#if defined(_OACC)
    !$acc end kernels
#endif
  end subroutine state_equal_to_state
end module module_types
