module module_output
  use calculation_types, only : wp, iowp
  use parallel_parameters, only : i_beg, k_beg
  use dimensions, only : nx, nz
  use module_types, only : atmospheric_state, reference_state
  use iodir, only : stderr
  use netcdf

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

  contains

  subroutine create_output
    implicit none
    integer :: t_dimid, x_dimid, z_dimid
    allocate(dens(nx,nz))
    allocate(uwnd(nx,nz))
    allocate(wwnd(nx,nz))
    allocate(theta(nx,nz))
    call ncwrap(nf90_create('output.nc',nf90_clobber,ncid), __LINE__)
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
  end subroutine create_output

  subroutine write_record(atmostat,ref,etime)
    implicit none
    type(atmospheric_state), intent(in) :: atmostat
    type(reference_state), intent(in) :: ref
    real(wp), intent(in) :: etime
    integer :: i, k
    integer, dimension(1) :: st1, ct1
    integer, dimension(3) :: st3, ct3
    real(wp), dimension(1) :: etimearr

    do k = 1, nz
      do i = 1, nx
        dens(i,k) = atmostat%dens(i,k)
        uwnd(i,k) = atmostat%umom(i,k)/(ref%density(k)+dens(i,k))
        wwnd(i,k) = atmostat%wmom(i,k)/(ref%density(k)+dens(i,k))
        theta(i,k) = (atmostat%rhot(i,k) + ref%denstheta(k)) / &
            (ref%density(k) + dens(i,k)) - ref%denstheta(k)/ref%density(k)
      end do
    end do

    st3 = [ i_beg, k_beg, rec_out ]
    ct3 = [ nx, nz, 1 ]
    call ncwrap(nf90_put_var(ncid,dens_varid,dens,st3,ct3), __LINE__)
    call ncwrap(nf90_put_var(ncid,uwnd_varid,uwnd,st3,ct3), __LINE__)
    call ncwrap(nf90_put_var(ncid,wwnd_varid,wwnd,st3,ct3), __LINE__)
    call ncwrap(nf90_put_var(ncid,theta_varid,theta,st3,ct3), __LINE__)

    st1 = [ rec_out ]
    ct1 = [ 1 ]
    etimearr(1) = etime
    call ncwrap(nf90_put_var(ncid,t_varid,etimearr,st1,ct1), __LINE__)

    rec_out = rec_out + 1
  end subroutine write_record

  subroutine close_output
    implicit none
    if ( allocated(dens) ) then
      deallocate(dens)
      deallocate(uwnd)
      deallocate(wwnd)
      deallocate(theta)
    end if
    call ncwrap(nf90_close(ncid), __LINE__)
  end subroutine close_output

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
