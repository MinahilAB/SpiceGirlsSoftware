!> @file module_output.f90
!> @brief Handles parallel NetCDF output operations, including file creation, data writing, and cleanup.
!>
!> This module manages the creation of the output file, defines NetCDF dimensions and variables, 
!> extracts diagnostic variables from the atmospheric state, and writes the data to disk using 
!> parallel NetCDF (PNETCDF) capabilities.
module module_output
  use calculation_types, only : wp, iowp
  use parallel_parameters, only : i_beg, k_beg, cart_comm, rank
  use dimensions, only : nx, nx_loc, nz
  use module_types, only : atmospheric_state, reference_state
  use iodir, only : stderr
  use indexing, only : I_DENS, I_UMOM, I_WMOM, I_RHOT
  use netcdf
  use module_nvtx
  use mpi
  
  implicit none

  private

  public :: create_output
  public :: write_record
  public :: close_output

  real(wp), allocatable, dimension(:,:) :: dens
  real(wp), allocatable, dimension(:,:) :: uwnd
  real(wp), allocatable, dimension(:,:) :: wwnd
  real(wp), allocatable, dimension(:,:) :: theta

  integer :: ncid
  integer :: dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid
  integer :: rec_out
  integer :: nc_comm

  contains

  !> @brief Creates the NetCDF output file and defines dimensions and variables.
  !>
  !> Sets up parallel NetCDF file creation and defines 4D variables (X, Z, Time) for 
  !> density, velocity components, and potential temperature.
  subroutine create_output
    implicit none
    integer :: t_dimid, x_dimid, z_dimid
    integer :: ierr

    call nvtx_push('create_output')

    allocate(dens(nx_loc,nz))
    allocate(uwnd(nx_loc,nz))
    allocate(wwnd(nx_loc,nz))
    allocate(theta(nx_loc,nz))
#if defined(_OACC)
    !$acc enter data create(dens, uwnd, wwnd, theta)
#endif

    call mpi_comm_dup(cart_comm, nc_comm, ierr)
    call ncwrap(nf90_create_par('output.nc', nf90_clobber, nc_comm, MPI_INFO_NULL, ncid), __LINE__)
    
    
    call ncwrap(nf90_def_dim(ncid,'time',nf90_unlimited,t_dimid), __LINE__)
    call ncwrap(nf90_def_dim(ncid,'x',nx,x_dimid), __LINE__)
    call ncwrap(nf90_def_dim(ncid,'z',nz,z_dimid), __LINE__)
    call ncwrap(nf90_def_var(ncid,'time',iowp,[t_dimid],t_varid), __LINE__)
    call ncwrap(nf90_def_var(ncid,'rho',iowp, &
        [x_dimid,z_dimid,t_dimid],dens_varid), __LINE__)
    call ncwrap(nf90_def_var(ncid,'u',iowp, &
        [x_dimid,z_dimid,t_dimid],uwnd_varid), __LINE__)
    call ncwrap(nf90_def_var(ncid,'w',iowp, &
        [x_dimid,z_dimid,t_dimid],wwnd_varid), __LINE__)
    call ncwrap(nf90_def_var(ncid,'theta',iowp, &
        [x_dimid,z_dimid,t_dimid],theta_varid), __LINE__)
    call ncwrap(nf90_enddef(ncid), __LINE__)
    rec_out = 1

    call nvtx_pop()
  end subroutine create_output


  !> @brief Writes a single record of the atmospheric state to the NetCDF file.
  !>
  !> Calculates diagnostic fields (velocity, potential temperature perturbation) 
  !> from prognostic variables and reference state, then writes the data using 
  !> parallel I/O with appropriate start (\c st3) and count (\c ct3) indices for the domain decomposition.
  !> @param[in] atmostat The current atmospheric state (prognostic variables).
  !> @param[in] ref The reference hydrostatic state.
  !> @param[in] etime The current elapsed simulation time.
  subroutine write_record(atmostat,ref,etime)
    implicit none
    type(atmospheric_state), intent(in) :: atmostat
    type(reference_state), intent(in) :: ref
    real(wp), intent(in) :: etime
    integer :: i, k
    integer, dimension(1) :: st1, ct1
    integer, dimension(3) :: st3, ct3
    real(wp), dimension(1) :: etimearr

    call nvtx_push('write_record')

#if defined(_OACC)
  !$acc parallel present(atmostat%mem, ref%density, ref%denstheta, dens, uwnd, wwnd, theta)
  !$acc loop gang vector collapse(2) private(i, k)
#elif defined(_OMP)
  !$omp parallel do collapse(2) private(i,k)
#endif
    do k = 1, nz
      do i = 1, nx_loc
        dens(i,k) = atmostat%mem(i, k, I_DENS)
        uwnd(i,k) = atmostat%mem(i, k, I_UMOM)/(ref%density(k)+dens(i,k))
        wwnd(i,k) = atmostat%mem(i, k, I_WMOM)/(ref%density(k)+dens(i,k))
        theta(i,k) = (atmostat%mem(i, k, I_RHOT) + ref%denstheta(k)) / &
            (ref%density(k) + dens(i,k)) - ref%denstheta(k)/ref%density(k)
      end do
    end do
#if defined(_OACC)
    !$acc end parallel
    !$acc update self(dens, uwnd, wwnd, theta)
#elif defined(_OMP)
    !$omp end parallel do
#endif

    st3 = [ i_beg, k_beg, rec_out ]
    ct3 = [ nx_loc, nz, 1 ]
    call ncwrap(nf90_put_var(ncid,dens_varid,dens,st3,ct3), __LINE__)
    call ncwrap(nf90_put_var(ncid,uwnd_varid,uwnd,st3,ct3), __LINE__)
    call ncwrap(nf90_put_var(ncid,wwnd_varid,wwnd,st3,ct3), __LINE__)
    call ncwrap(nf90_put_var(ncid,theta_varid,theta,st3,ct3), __LINE__)

    if (rank == 0) then
      st1 = [ rec_out ]
      ct1 = [ 1 ]
      etimearr(1) = etime
      call ncwrap(nf90_put_var(ncid,t_varid,etimearr,st1,ct1), __LINE__)
    end if

    rec_out = rec_out + 1

    call nvtx_pop()
  end subroutine write_record


  !> @brief Closes the NetCDF file, frees the MPI communicator, and deallocates arrays.
  subroutine close_output
    implicit none
    integer :: ierr

    call nvtx_push('close_output')

    if ( allocated(dens) ) then
      
#if defined(_OACC)
      !$acc exit data delete(dens, uwnd, wwnd, theta)
#endif

      deallocate(dens)
      deallocate(uwnd)
      deallocate(wwnd)
      deallocate(theta)
    end if

    call ncwrap(nf90_close(ncid), __LINE__)
    call mpi_comm_free(nc_comm, ierr)

    call nvtx_pop()
  end subroutine close_output


  !> @brief NetCDF error-wrapping subroutine.
  !>
  !> Prints the NetCDF error string and stops execution if an error is detected.
  !> @param[in] ierr The NetCDF status code.
  !> @param[in] line The line number where the NetCDF call was made.
  subroutine ncwrap(ierr,line)
    implicit none
    integer, intent(in) :: ierr
    integer, intent(in) :: line
    if (ierr /= nf90_noerr) then
      write(stderr,*) 'NetCDF Error at line: ', line
      write(stderr,*) nf90_strerror(ierr)
      stop
    end if
  end subroutine ncwrap

end module module_output
