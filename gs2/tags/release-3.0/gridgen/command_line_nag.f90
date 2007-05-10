module command_line
  use f90_unix, only: nag_iargc => iargc, getarg
  implicit none

contains

  function iargc ()
    integer :: iargc
    iargc = nag_iargc()
  end function iargc

  subroutine cl_getarg (k, arg, len, ierr)
    implicit none
    integer,       intent (in)  :: k
    character (*), intent (out) :: arg
    integer,       intent (in)  :: len
    integer,       intent (out) :: ierr

    call getarg (k, arg)
    ierr = 0
  end subroutine cl_getarg
end module command_line
