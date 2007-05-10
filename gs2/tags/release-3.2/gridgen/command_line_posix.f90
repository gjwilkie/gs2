module command_line
  implicit none
  interface
     function iargc ()
       integer :: iargc
     end function iargc
  end interface
contains
  subroutine cl_getarg (k, arg, len, ierr)
    implicit none
    integer,       intent (in)  :: k
    character (*), intent (out) :: arg
    integer,       intent (in)  :: len
    integer,       intent (out) :: ierr

    interface
       subroutine pxfgetarg (m, buf, ilen, ierror)
         integer,       intent (in)  :: m
         character (*), intent (out) :: buf
         integer,       intent (in)  :: ilen
         integer,       intent (out) :: ierror
       end subroutine pxfgetarg
    end interface

    call pxfgetarg (k, arg, len, ierr)
  end subroutine cl_getarg
end module command_line
