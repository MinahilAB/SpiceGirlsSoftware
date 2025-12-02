module parallel_timer
  use iso_c_binding
  implicit none
  private
  public :: timer_type, mytimer_create, mytimer_destroy, mytimer_print, mytimer_gather_stats

  type, bind(C) :: timer_type
     type(c_ptr) :: handle = c_null_ptr
  end type timer_type

  interface
     function mytimer_create_c(name) bind(C, name="mytimer_start")
        import :: c_ptr, c_char
        type(c_ptr) :: mytimer_create_c
        character(kind=c_char), intent(in) :: name(*)
     end function mytimer_create_c

     subroutine mytimer_destroy_c(handle) bind(C, name="mytimer_stop")
        import :: c_ptr
        type(c_ptr), value :: handle
     end subroutine mytimer_destroy_c

     subroutine mytimer_print_c() bind(C, name="mytimer_print")
     end subroutine mytimer_print_c

     subroutine mytimer_gather_c() bind(C, name="mytimer_gather_and_print")
     end subroutine mytimer_gather_c
  end interface

contains

  subroutine mytimer_create(timer, name)
    type(timer_type), intent(out) :: timer
    character(*), intent(in) :: name
    character(kind=c_char), dimension(:), allocatable :: cname
    integer :: n
    n = len(name)
    allocate(cname(n+1))
    cname(1:n) = transfer(name, cname(1:n))
    cname(n+1) = c_null_char
    timer%handle = mytimer_create_c(cname)
  end subroutine mytimer_create

  subroutine mytimer_destroy(timer)
    type(timer_type), intent(inout) :: timer
    if (c_associated(timer%handle)) then
       call mytimer_destroy_c(timer%handle)
       timer%handle = c_null_ptr
    end if
  end subroutine mytimer_destroy

  subroutine mytimer_print()
    call mytimer_print_c()
  end subroutine mytimer_print

  subroutine mytimer_gather_stats()
    call mytimer_gather_c()
  end subroutine mytimer_gather_stats

end module parallel_timer
