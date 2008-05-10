module command_line
  implicit none
!  interface
!     function iargc ()
!       integer :: iargc
!     end function iargc
!  end interface
contains
  subroutine cl_getarg (k, arg, len, ierr)
    implicit none
    integer,       intent (in)  :: k
    character (*), intent (out) :: arg
    integer,       intent (in)  :: len
    integer,       intent (out) :: ierr

!    interface
!       subroutine getarg (m, buf)
!         integer,       intent (in)  :: m
!         character (*), intent (out) :: buf
!       end subroutine getarg
!    end interface

    call getarg (k, arg)
    ierr = 0
  end subroutine cl_getarg
end module command_line
