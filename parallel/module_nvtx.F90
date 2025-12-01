module module_nvtx
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  implicit none

#ifdef USE_NVTX
  interface
     function nvtxRangePushA(name) bind(C, name="nvtxRangePushA")
       import :: c_int, c_char
       integer(c_int) :: nvtxRangePushA
       character(kind=c_char), dimension(*) :: name
     end function nvtxRangePushA

     function nvtxRangePop() bind(C, name="nvtxRangePop")
       import :: c_int
       integer(c_int) :: nvtxRangePop
     end function nvtxRangePop
  end interface
#endif

contains

  subroutine nvtx_push(label)
    character(len=*), intent(in) :: label
#ifdef USE_NVTX
    character(kind=c_char,len=:), allocatable :: c_label
    integer(c_int) :: istat

    ! Convert Fortran string to C string (null-terminated)
    c_label = trim(label)//c_null_char
    istat = nvtxRangePushA(c_label)
#else
    ! NVTX disabled: no-op
#endif
  end subroutine nvtx_push

  subroutine nvtx_pop()
#ifdef USE_NVTX
    integer(c_int) :: istat
    istat = nvtxRangePop()
#else
    ! NVTX disabled: no-op
#endif
  end subroutine nvtx_pop

end module module_nvtx