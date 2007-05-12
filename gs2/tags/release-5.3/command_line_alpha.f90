module command_line
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
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
       subroutine getarg (m, buf)
         integer,       intent (in)  :: m
         character (*), intent (out) :: buf
       end subroutine getarg
    end interface

    call getarg (k, arg)
    ierr = 0
  end subroutine cl_getarg

  subroutine ishell(string)
    ! execute the command in string in a system subshell

    character*(*), intent(in) :: string
    call system(string)
  end subroutine ishell

  real function second()
    implicit none
    real*4 etime,tarray(2)
    second = etime(tarray)
    return
  end function second

end module command_line
