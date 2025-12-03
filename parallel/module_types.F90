module module_types
  use calculation_types
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use legendre_quadrature
  use dimensions
  use iodir
  use mpi 

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
  end type atmospheric_state

  type atmospheric_flux
    real(wp), pointer, dimension(:,:,:) :: mem => null( )
    real(wp), pointer, dimension(:,:) :: dens
    real(wp), pointer, dimension(:,:) :: umom
    real(wp), pointer, dimension(:,:) :: wmom
    real(wp), pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_flux
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
    procedure, public :: del_tendency
    procedure, public :: set_tendency
  end type atmospheric_tendency

  interface assignment(=)
    module procedure state_equal_to_state
  end interface

contains

  ! ------------------------------------------------------------------
  ! Reference State
  ! ------------------------------------------------------------------
  subroutine new_ref(ref, nz)
    implicit none
    class(reference_state), intent(inout) :: ref
    integer, intent(in) :: nz
    allocate(ref%density(nz))
    allocate(ref%denstheta(nz))
    allocate(ref%idens(nz))
    allocate(ref%idenstheta(nz))
    allocate(ref%pressure(nz))
    !$acc enter data copyin(ref)
    !$acc enter data create(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
  end subroutine new_ref

  subroutine del_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
    !$acc exit data delete(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc exit data delete(ref)
    if (allocated(ref%density)) deallocate(ref%density)
    if (allocated(ref%denstheta)) deallocate(ref%denstheta)
    if (allocated(ref%idens)) deallocate(ref%idens)
    if (allocated(ref%idenstheta)) deallocate(ref%idenstheta)
    if (allocated(ref%pressure)) deallocate(ref%pressure)
  end subroutine del_ref

  ! ------------------------------------------------------------------
  ! Atmospheric State
  ! ------------------------------------------------------------------
  subroutine new_state(stat, nx, nz)
    implicit none
    class(atmospheric_state), intent(inout) :: stat
    integer, intent(in) :: nx, nz
    
    allocate(stat%mem(nx, nz, NVARS))
    stat%dens => stat%mem(:,:,I_DENS)
    stat%umom => stat%mem(:,:,I_UMOM)
    stat%wmom => stat%mem(:,:,I_WMOM)
    stat%rhot => stat%mem(:,:,I_RHOT)
    stat%mem = 0.0_wp

#if defined(_OACC)
    ! Split directives to avoid "mixed component and non-component accesses" error
    !$acc enter data copyin(stat)
    !$acc enter data create(stat%mem)
    !$acc enter data attach(stat%dens, stat%umom, stat%wmom, stat%rhot)
#endif
  end subroutine new_state

  subroutine del_state(stat)
    implicit none
    class(atmospheric_state), intent(inout) :: stat
#if defined(_OACC)
    ! Split directives to avoid mixed access errors
    !$acc exit data delete(stat%mem)
    !$acc exit data delete(stat)
#endif
    if ( associated(stat%mem) ) deallocate(stat%mem)
    nullify(stat%dens)
    nullify(stat%umom)
    nullify(stat%wmom)
    nullify(stat%rhot)
  end subroutine del_state

  subroutine set_state(stat, xval)
    implicit none
    class(atmospheric_state), intent(inout) :: stat
    real(wp), intent(in) :: xval
    !$acc kernels present(stat)
    stat%mem(:,:,:) = xval
    !$acc end kernels
  end subroutine set_state

  subroutine update(stat, tend, dt)
    implicit none
    class(atmospheric_state), intent(inout) :: stat
    type(atmospheric_tendency), intent(in) :: tend
    real(wp), intent(in) :: dt
    integer :: i, j, k ! <--- FIXED: Added loop variables declarations

    !$acc parallel loop collapse(3) present(stat, tend)
    do k=1,nz
       do j=1,nx
          do i=1,NVARS
             stat%mem(j,k,i) = stat%mem(j,k,i) + dt * tend%mem(j,k,i)
          end do
       end do
    end do
  end subroutine update

  ! ------------------------------------------------------------------
  ! Halo Exchange
  ! ------------------------------------------------------------------
  subroutine exchange_halo_x(stat)
    implicit none
    class(atmospheric_state), intent(inout) :: stat
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    
#if defined(_OACC)
    ! Update host memory before MPI
    !$acc update self(stat%mem)
#endif

    ! Send Right, Recv Left
    call MPI_Sendrecv(stat%mem(nx_loc-hs+1,:,:), hs*nz*NVARS, MPI_DOUBLE_PRECISION, right_rank, 0, &
                      stat%mem(1,:,:),            hs*nz*NVARS, MPI_DOUBLE_PRECISION, left_rank,  0, &
                      cart_comm, status, ierr)

    ! Send Left, Recv Right
    call MPI_Sendrecv(stat%mem(hs+1,:,:),         hs*nz*NVARS, MPI_DOUBLE_PRECISION, left_rank,  1, &
                      stat%mem(nx_loc+hs+1,:,:),  hs*nz*NVARS, MPI_DOUBLE_PRECISION, right_rank, 1, &
                      cart_comm, status, ierr)

#if defined(_OACC)
    ! Update device memory after MPI
    !$acc update device(stat%mem)
#endif
  end subroutine exchange_halo_x

  ! ------------------------------------------------------------------
  ! Atmospheric Flux
  ! ------------------------------------------------------------------
  subroutine new_flux(flux, nx, nz)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    integer, intent(in) :: nx, nz
    allocate(flux%mem(nx, nz, NVARS))
    flux%dens => flux%mem(:,:,I_DENS)
    flux%umom => flux%mem(:,:,I_UMOM)
    flux%wmom => flux%mem(:,:,I_WMOM)
    flux%rhot => flux%mem(:,:,I_RHOT)
    flux%mem = 0.0_wp
#if defined(_OACC)
    ! Split directives
    !$acc enter data copyin(flux)
    !$acc enter data create(flux%mem)
    !$acc enter data attach(flux%dens, flux%umom, flux%wmom, flux%rhot)
#endif
  end subroutine new_flux

  subroutine del_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
#if defined(_OACC)
    ! Split directives
    !$acc exit data delete(flux%mem)
    !$acc exit data delete(flux)
#endif
    if ( associated(flux%mem) ) deallocate(flux%mem)
    nullify(flux%dens, flux%umom, flux%wmom, flux%rhot)
  end subroutine del_flux

  ! ------------------------------------------------------------------
  ! Atmospheric Tendency
  ! ------------------------------------------------------------------
  subroutine new_tendency(tend, nx, nz)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    integer, intent(in) :: nx, nz
    allocate(tend%mem(nx, nz, NVARS))
    tend%dens => tend%mem(:,:,I_DENS)
    tend%umom => tend%mem(:,:,I_UMOM)
    tend%wmom => tend%mem(:,:,I_WMOM)
    tend%rhot => tend%mem(:,:,I_RHOT)
    tend%mem = 0.0_wp
#if defined(_OACC)
    ! Split directives
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
       stop 'Tendency not allocated'
    end if
    !$acc kernels present(tend)
    tend%mem(:,:,:) = xval
    !$acc end kernels
  end subroutine set_tendency

  subroutine del_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
#if defined(_OACC)
    ! Split directives
    !$acc exit data delete(tend%mem)
    !$acc exit data delete(tend)
#endif
    if ( associated(tend%mem) ) deallocate(tend%mem)
    nullify(tend%dens, tend%umom, tend%wmom, tend%rhot)
  end subroutine del_tendency

  subroutine state_equal_to_state(x,y)
    implicit none
    class(atmospheric_state), intent(out) :: x
    class(atmospheric_state), intent(in)  :: y
    !$acc kernels present(x,y)
    x%mem = y%mem
    !$acc end kernels
  end subroutine state_equal_to_state

end module module_types
