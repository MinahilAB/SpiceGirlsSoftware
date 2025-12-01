module calculation_types
  use netcdf, only : nf90_real, nf90_double
  use iso_fortran_env
  implicit none
  public
  integer, parameter :: wp = real64
  integer, parameter :: iowp = nf90_double
  integer, parameter :: ip = int32
end module calculation_types

module physical_constants
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: pi = 3.14159265358979323846264338327_wp
  real(wp), parameter :: hpi = 0.5_wp * pi
  real(wp), parameter :: grav = 9.80665_wp
  real(wp), parameter :: boltzk = 1.3806490e-23_wp
  real(wp), parameter :: navgdr = 6.02214076e23_wp
  real(wp), parameter :: rgasmol = navgdr*boltzk
  real(wp), parameter :: amd = 28.96454_wp
  real(wp), parameter :: amw = 18.01528_wp
  real(wp), parameter :: rgas = (rgasmol/amd)*1000.0_wp
  real(wp), parameter :: cp = 3.5_wp*rgas
  real(wp), parameter :: cv = 2.5_wp*rgas
  real(wp), parameter :: rd = rgas
  real(wp), parameter :: p0 = 1.e5_wp
  real(wp), parameter :: t0 = 298.0_wp
  real(wp), parameter :: cdocv = 3.5_wp/2.5_wp
  real(wp), parameter :: cvocd = 2.5_wp/3.5_wp
  real(wp), parameter :: rdocp = rd/cp
  real(wp), parameter :: rdocv = rd/cv
  real(wp), parameter :: c0 = (rd**cdocv)*p0**(-rdocv)
  real(wp), parameter :: theta0 = 300.0_wp
  real(wp), parameter :: exner0 = 1.0_wp
end module physical_constants

module physical_parameters
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: xlen = 20000.0_wp
  real(wp), parameter :: zlen = 10000.0_wp
  real(wp), parameter :: hxlen = 0.5_wp*xlen
  real(wp), parameter :: p1 = xlen/10.0_wp
  real(wp), parameter :: p2 = 4.0_wp*p1
  real(wp), parameter :: p3 = 2.0_wp*p1
  real(wp), parameter :: mnt_width = xlen/8.0_wp
  real(wp), parameter :: hv_beta = 0.05_wp
  real(wp), parameter :: cfl = 1.5_wp
  real(wp), parameter :: max_speed = 450.0_wp
  real(wp), parameter :: x0 = mnt_width
  real(wp), parameter :: z0 = 1000.0_wp
  real(wp), parameter :: xrad = 500.0_wp
  real(wp), parameter :: zrad = 500.0_wp
  real(wp), parameter :: amp = 0.01_wp
end module physical_parameters

module indexing
  implicit none
  public
  integer, parameter :: NVARS = 4
  integer, parameter :: STEN_SIZE = 4
  integer, parameter :: I_DENS  = 1
  integer, parameter :: I_UMOM  = 2
  integer, parameter :: I_WMOM  = 3
  integer, parameter :: I_RHOT  = 4
  integer, parameter :: DIR_X = 1
  integer, parameter :: DIR_Z = 2
end module indexing

module iodir
  use iso_fortran_env
  public
  integer, public, parameter :: stdin = input_unit
  integer, public, parameter :: stdout = output_unit
  integer, public, parameter :: stderr = error_unit
end module iodir

module parallel_parameters
  implicit none
  public
  integer, parameter :: i_beg = 1
  integer, parameter :: k_beg = 1
  integer, parameter :: hs = 2
end module parallel_parameters

module legendre_quadrature
  use calculation_types, only : wp
  implicit none
  public
  integer , parameter :: nqpoints = 3
  real(wp), parameter :: qpoints(nqpoints) =    &
      [ 0.112701665379258311482073460022E0_wp , &
        0.500000000000000000000000000000E0_wp , &
        0.887298334620741688517926539980E0_wp ]
  real(wp), parameter :: qweights(nqpoints) =   &
      [ 0.277777777777777777777777777779E0_wp , &
        0.444444444444444444444444444444E0_wp , &
        0.277777777777777777777777777779E0_wp ]
end module legendre_quadrature

module dimensions
  use calculation_types, only : wp
  use physical_parameters, only : zlen, xlen
  use indexing
  implicit none
  public
  integer , parameter :: nx = 100
  integer , parameter :: nz = int(nx * zlen/xlen)
  real(wp), parameter :: sim_time = 1000.0_wp
  real(wp), parameter :: output_freq = 10.0_wp
end module dimensions
