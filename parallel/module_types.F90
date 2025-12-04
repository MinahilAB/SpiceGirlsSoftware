!> @file module_types.F90
!> @brief Defines all core derived types for the atmospheric model state, fluxes, and tendencies, 
!> along with parallel data management and kernel execution routines.
!>
!> This module establishes the data structures and core subroutines (like \c xtend, \c ztend, 
!> and halo exchange) that drive the 3D stencil computation. It is crucial for understanding 
!> the data layout and the structure of the time integration step.

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

  !> @brief Derived type for storing the constant reference atmospheric state.
  !>
  !> This state is typically independent of time and is used to linearize the equations
  !> or define background fields.
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


  !> @brief Derived type for the time-dependent atmospheric state variables.
  !>
  !> This type holds the core prognostic variables (\c dens, \c umom, \c wmom, \c rhot)
  !> and includes halo regions (\c hs) for parallel exchange.
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


   !> @brief Derived type for storing calculated fluxes.
  !>
  !> Fluxes are stored on cell faces/interfaces, hence their dimensioning (\c nx_loc+1, \c nz+1) 
  !> often differs slightly from the state.
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


  !> @brief Derived type for storing calculated tendencies (time derivatives).
  !>
  !> Tendency fields are computed in the interior cells only.
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

  
#if defined(_OACC)
  public :: send_left_d, send_right_d, recv_left_d, recv_right_d
  real(wp), allocatable :: send_left_d(:), send_right_d(:), recv_left_d(:), recv_right_d(:)
#else
  public :: send_left, send_right, recv_left, recv_right
  real(wp), allocatable :: send_left(:), send_right(:), recv_left(:), recv_right(:)
#endif

  contains

  !> @brief Allocates memory for an atmospheric state type instance.
  !> @param[inout] atmo The atmospheric state instance.
  subroutine new_state(atmo)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    allocate(atmo%mem(1-hs:nx_loc+hs, 1-hs:nz+hs, NVARS))
    atmo%dens(1-hs:,1-hs:) => atmo%mem(:,:,I_DENS)
    atmo%umom(1-hs:,1-hs:) => atmo%mem(:,:,I_UMOM)
    atmo%wmom(1-hs:,1-hs:) => atmo%mem(:,:,I_WMOM)
    atmo%rhot(1-hs:,1-hs:) => atmo%mem(:,:,I_RHOT)
  end subroutine new_state


  !> @brief Sets all grid points of a state to a scalar value.
  !> @param[inout] atmo The atmospheric state instance.
  !> @param[in] xval The scalar value to set.
  subroutine set_state(atmo, xval)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    real(wp), intent(in) :: xval

    if ( .not. associated(atmo%mem) ) then
      write(stderr,*) 'NOT ALLOCATED STATE ERROR AT LINE ', __LINE__
      stop
    end if

    atmo%mem(1-hs:nx_loc+hs, 1-hs:nz+hs, :) = xval
  end subroutine set_state


  !> @brief Deallocates memory for an atmospheric state type instance.
  !> @param[inout] atmo The atmospheric state instance.
  subroutine del_state(atmo)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    nullify(atmo%dens)
    nullify(atmo%umom)
    nullify(atmo%wmom)
    nullify(atmo%rhot)
  end subroutine del_state


  !> @brief Performs a time integration update step.
  !>
  !> Uses a simple explicit time integration formula: \f$ S^{n+1} = S^n + \Delta t \cdot Tendency \f$.
  !> The wall-clock time is captured by the NVTX range \c :update.
  !>
  !> @param[inout] s2 The state at the new time level (\f$S^{n+1}\f$).
  !> @param[in] s0 The state at the current time level (\f$S^n\f$).
  !> @param[in] tend The calculated tendency (\f$ \frac{dS}{dt} \f$).
  !> @param[in] dt The time step size (\f$\Delta t\f$).
  subroutine update(s2,s0,tend,dt)
    implicit none
    class(atmospheric_state), intent(inout) :: s2
    class(atmospheric_state), intent(in) :: s0
    class(atmospheric_tendency), intent(in) :: tend
    real(wp), intent(in) :: dt
    integer :: ll, k, i
    
#if defined(_OACC)
    !$acc parallel loop collapse(3) present(s2%mem, s0%mem, tend%mem)
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
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine update


  !> @brief Calculates the X-dimension contribution to the total tendency.
  !> @ingroup PerformanceRanges
  !>
  !> This routine involves two main steps: flux calculation (using a 4-point stencil) and 
  !> tendency calculation (flux difference). It requires an X-halo exchange beforehand.
  !> The wall-clock time is captured by the NVTX range \c :xtend.
  !>
  !> @param[inout] tendency The resulting time tendency array.
  !> @param[inout] flux Temporary flux array used for calculation.
  !> @param[in] ref Reference state data.
  !> @param[inout] atmostat Current atmospheric state (required for halo data).
  !> @param[in] dx Grid spacing in the X-direction.
  !> @param[in] dt Time step size (used for hyperviscosity coefficient).
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
    !$acc parallel loop collapse(2) &
    !$acc& present(flux%mem, atmostat%mem, ref%density, ref%denstheta) &
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
        flux%mem(i, k, I_DENS) = r*u - hv_coef*d3_vals(I_DENS)
        flux%mem(i, k, I_UMOM) = r*u*u+p - hv_coef*d3_vals(I_UMOM)
        flux%mem(i, k, I_WMOM) = r*u*w - hv_coef*d3_vals(I_WMOM)
        flux%mem(i, k, I_RHOT) = r*u*t - hv_coef*d3_vals(I_RHOT)
      end do
    end do
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif

    ! Second loop for tendency calculation
#if defined(_OACC)
    !$acc parallel loop collapse(3) present(tendency%mem, flux%mem)
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
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine xtend


  !> @brief Calculates the Z-dimension contribution to the total tendency.
  !> @ingroup PerformanceRanges
  !>
  !> This routine involves two main steps: flux calculation (using a 4-point stencil) and 
  !> tendency calculation (flux difference), plus the application of gravity as a source term.
  !> It requires a Z-halo exchange beforehand (local boundary condition setting).
  !> The wall-clock time is captured by the NVTX range \c :ztend.
  !>
  !> @param[inout] tendency The resulting time tendency array.
  !> @param[inout] flux Temporary flux array used for calculation.
  !> @param[in] ref Reference state data.
  !> @param[inout] atmostat Current atmospheric state (required for halo data).
  !> @param[in] dz Grid spacing in the Z-direction.
  !> @param[in] dt Time step size (used for hyperviscosity coefficient).

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
    !$acc parallel loop collapse(2) present(flux%mem, atmostat%mem, ref%idens, ref%idenstheta, ref%pressure) &
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
        flux%mem(i, k, I_DENS) = r*w - hv_coef*d3_vals(I_DENS)
        flux%mem(i, k, I_UMOM) = r*w*u - hv_coef*d3_vals(I_UMOM)
        flux%mem(i, k, I_WMOM) = r*w*w+p - hv_coef*d3_vals(I_WMOM)
        flux%mem(i, k, I_RHOT) = r*w*t - hv_coef*d3_vals(I_RHOT)
      end do
    end do
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif

    ! Second loop for tendency calculation
#if defined(_OACC)
    !$acc parallel loop collapse(3) present(tendency%mem, flux%mem, atmostat%mem)
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
            tendency%mem(i,k,ll) = tendency%mem(i,k,ll) - atmostat%mem(i,k,I_DENS)*grav
          end if
        end do
      end do
    end do
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine ztend


  !> @brief Performs the MPI-based halo exchange in the X-dimension.
  !> @ingroup PerformanceRanges
  !>
  !> This routine is responsible for point-to-point communication with neighboring ranks.
  !> It involves packing data, initiating non-blocking sends/receives, and waiting for completion.
  !> The wall-clock time is captured by the NVTX range \c :exchange_halo_x.
  !>
  !> @param[inout] s The atmospheric state whose halo regions need updating.
  subroutine exchange_halo_x(s)
    use dimensions, only : nz, nx_loc
    class(atmospheric_state), intent(inout) :: s
    integer :: ncount
    integer :: reqs(4), ierr
    integer :: ll, k, i, idx

    call nvtx_push('exchange_halo')

    ncount = hs * nz * NVARS

#if defined(_OACC)
    ! Pack left
    !$acc parallel loop collapse(3) present(s, send_left_d)
    do ll = 1, NVARS
      do k = 1, nz
        do i = 1, hs
          idx = ((ll - 1) * nz + (k - 1)) * hs + i
          send_left_d(idx) = s%mem(i, k, ll)
        end do
      end do
    end do
    !$acc end parallel

    ! Pack right
    !$acc parallel loop collapse(3) present(s, send_right_d)
    do ll = 1, NVARS
      do k = 1, nz
        do i = 1, hs
          idx = ((ll - 1) * nz + (k - 1)) * hs + i
          send_right_d(idx) = s%mem(nx_loc - hs + i, k, ll)
        end do
      end do
    end do
    !$acc end parallel

    ! Send/Recv
    !$acc host_data use_device(send_left_d, send_right_d, recv_left_d, recv_right_d)
      call MPI_Irecv(recv_left_d,  ncount, MPI_DOUBLE_PRECISION, left_rank,  102, cart_comm, reqs(1), ierr)
      call MPI_Irecv(recv_right_d, ncount, MPI_DOUBLE_PRECISION, right_rank, 101, cart_comm, reqs(2), ierr)
      call MPI_Isend(send_right_d, ncount, MPI_DOUBLE_PRECISION, right_rank, 102, cart_comm, reqs(3), ierr)
      call MPI_Isend(send_left_d,  ncount, MPI_DOUBLE_PRECISION, left_rank,  101, cart_comm, reqs(4), ierr)
    !$acc end host_data

    call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)

    ! Unpack received left
    !$acc parallel loop collapse(3) present(s, recv_left_d)
      do ll = 1, NVARS
        do k = 1, nz
          do i = 1, hs
            idx = ((ll - 1) * nz + (k - 1)) * hs + i
            s%mem(1 - hs + (i - 1), k, ll) = recv_left_d(idx)
            ! equivalently: s%mem(i - hs, k, ll)
          end do
        end do
      end do
    !$acc end parallel

    ! Unpack received right
    !$acc parallel loop collapse(3) present(s, recv_right_d)
      do ll = 1, NVARS
        do k = 1, nz
          do i = 1, hs
            idx = ((ll - 1) * nz + (k - 1)) * hs + i
            s%mem(nx_loc + i, k, ll) = recv_right_d(idx)
          end do
        end do
      end do
    !$acc end parallel

#else

    call pack_strip(s%mem, 1,           1, nx_loc, hs, send_left)
    call pack_strip(s%mem, nx_loc-hs+1, 1, nx_loc, hs, send_right)
      
    call MPI_Irecv(recv_left,  ncount, MPI_DOUBLE_PRECISION, left_rank,  102, cart_comm, reqs(1), ierr)
    call MPI_Irecv(recv_right, ncount, MPI_DOUBLE_PRECISION, right_rank, 101, cart_comm, reqs(2), ierr)
    call MPI_Isend(send_right, ncount, MPI_DOUBLE_PRECISION, right_rank, 102, cart_comm, reqs(3), ierr)
    call MPI_Isend(send_left,  ncount, MPI_DOUBLE_PRECISION, left_rank,  101, cart_comm, reqs(4), ierr)
      
    call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
      
    call unpack_strip(recv_left,  s%mem, 1-hs,    1)
    call unpack_strip(recv_right, s%mem, nx_loc+1,1)
    
#endif

    call nvtx_pop()
  end subroutine exchange_halo_x


   !> @brief Packs a strip of data (e.g., the halo region) into a 1D send buffer.
  !> @ingroup PerformanceRanges
  !>
  !> The wall-clock time is captured by the NVTX range \c :pack_strip.
  !> @param[in] mem The source 3D memory block.
  !> @param[in] i_start Starting index in the X-direction.
  !> @param[in] k_start Starting index in the Z-direction.
  !> @param[in] nx_loc Local X-dimension size.
  !> @param[in] width Width of the strip (usually \c hs).
  !> @param[out] buffer The resulting 1D communication buffer.
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


   !> @brief Unpacks data from a 1D receive buffer into the state's halo region.
  !> @ingroup PerformanceRanges
  !>
  !> The wall-clock time is captured by the NVTX range \c :unpack_strip.
  !> @param[in] buffer The 1D communication buffer received from a neighbor.
  !> @param[inout] mem The target 3D memory block (where the halo is written).
  !> @param[in] i_start Starting index in the X-direction for the target halo.
  !> @param[in] k_start Starting index in the Z-direction for the target halo.
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


  !> @brief Implements the Z-dimension (vertical) halo exchange.
  !>
  !> Since the Z-dimension is not parallelized, this involves applying local, non-communicative 
  !> boundary conditions (e.g., free slip, fixed values) to the halo regions.
  !> The wall-clock time is captured by the NVTX range \c :exchange_halo_z.
  !>
  !> @param[inout] s The atmospheric state to apply boundary conditions to.
  !> @param[in] ref The reference state used to calculate some boundary values.
  subroutine exchange_halo_z(s,ref)
    implicit none
    class(atmospheric_state), intent(inout) :: s
    class(reference_state), intent(in) :: ref
    integer :: i, ll
#if defined(_OACC)
    !$acc parallel loop collapse(2) present(s%mem, ref%density)
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
#if defined(_OACC)
    !$acc end parallel
#elif defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine exchange_halo_z


  !> @brief Allocates memory for the reference state.
  !> @param[inout] ref The reference state instance.
  subroutine new_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
    allocate(ref%density(1-hs:nz+hs))
    allocate(ref%denstheta(1-hs:nz+hs))
    allocate(ref%idens(nz+1))
    allocate(ref%idenstheta(nz+1))
    allocate(ref%pressure(nz+1))
  end subroutine new_ref


  !> @brief Deallocates memory for the reference state.
  !> @param[inout] ref The reference state instance.
  subroutine del_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
    deallocate(ref%density)
    deallocate(ref%denstheta)
    deallocate(ref%idens)
    deallocate(ref%idenstheta)
    deallocate(ref%pressure)
  end subroutine del_ref


  !> @brief Allocates memory for an atmospheric flux type instance.
  !> @param[inout] flux The atmospheric flux instance.
  subroutine new_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) deallocate(flux%mem)
    allocate(flux%mem(1:nx_loc+1, 1:nz+1,NVARS))
    flux%dens => flux%mem(:,:,I_DENS)
    flux%umom => flux%mem(:,:,I_UMOM)
    flux%wmom => flux%mem(:,:,I_WMOM)
    flux%rhot => flux%mem(:,:,I_RHOT)
  end subroutine new_flux


  !> @brief Sets all grid points of a flux field to a scalar value.
  !> @param[inout] flux The atmospheric flux instance.
  !> @param[in] xval The scalar value to set.
  subroutine set_flux(flux, xval)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    real(wp), intent(in) :: xval
    if ( .not. associated(flux%mem) ) then
      write(stderr,*) 'NOT ALLOCATED FLUX ERROR AT LINE ', __LINE__
      stop
    end if

#if defined(_OACC)
    !$acc kernels present(flux%mem)
#endif

  flux%mem(:,:,:) = xval

#if defined(_OACC)
  !$acc end kernels
#endif

  end subroutine set_flux


  !> @brief Deallocates memory for an atmospheric flux type instance.
  !> @param[inout] flux The atmospheric flux instance.
  subroutine del_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) deallocate(flux%mem)
    nullify(flux%dens)
    nullify(flux%umom)
    nullify(flux%wmom)
    nullify(flux%rhot)
  end subroutine del_flux


  !> @brief Allocates memory for an atmospheric tendency type instance.
  !> @param[inout] tend The atmospheric tendency instance.
  subroutine new_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) deallocate(tend%mem)
    allocate(tend%mem(nx_loc, nz,NVARS))
    tend%dens => tend%mem(:,:,I_DENS)
    tend%umom => tend%mem(:,:,I_UMOM)
    tend%wmom => tend%mem(:,:,I_WMOM)
    tend%rhot => tend%mem(:,:,I_RHOT)
  end subroutine new_tendency


  !> @brief Sets all grid points of a tendency field to a scalar value.
  !> @param[inout] tend The atmospheric tendency instance.
  !> @param[in] xval The scalar value to set.
  subroutine set_tendency(tend, xval)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    real(wp), intent(in) :: xval
    if ( .not. associated(tend%mem) ) then
      write(stderr,*) 'NOT ALLOCATED FLUX ERROR AT LINE ', __LINE__
      stop
    end if

#if defined(_OACC)
    !$acc kernels present(tend%mem)
#endif

    tend%mem(:,:,:) = xval

#if defined(_OACC)
    !$acc end kernels
#endif

  end subroutine set_tendency


  !> @brief Deallocates memory for an atmospheric tendency type instance.
  !> @param[inout] tend The atmospheric tendency instance.
  subroutine del_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) deallocate(tend%mem)
    nullify(tend%dens)
    nullify(tend%umom)
    nullify(tend%wmom)
    nullify(tend%rhot)
  end subroutine del_tendency


  !> @brief Overloaded assignment operator to copy one state to another.
  !> @param[inout] x The destination state.
  !> @param[in] y The source state.
  subroutine state_equal_to_state(x,y)
    implicit none
    type(atmospheric_state), intent(inout) :: x
    type(atmospheric_state), intent(in) :: y
    integer :: ll, k, i

#if defined(_OMP)
    !$omp parallel do collapse(3) default(none) shared(x, y, nx_loc, nz, hs, NVARS) private(ll, k, i)
#endif
    do ll = 1, NVARS
      do k = 1-hs, nz+hs
        do i = 1-hs, nx_loc+hs
          x%mem(i,k,ll) = y%mem(i,k,ll)
        end do
      end do
    end do
#if defined(_OMP)
    !$omp end parallel do
#endif
  end subroutine state_equal_to_state

  
end module module_types
