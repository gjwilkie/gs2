# include "define.inc"

module command_line
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!
# if FCOMPILER == _NAG_ && ( ! defined POSIX )
  use f90_unix, only: nag_iargc => iargc, getarg
# endif
  implicit none

# if ( defined POSIX ) || ( FCOMPILER != _GFORTRAN_ && FCOMPILER != _NAG_ )
  interface
     function iargc ()
       integer :: iargc
     end function iargc
  end interface
# endif

contains

# if FCOMPILER == _NAG_ && ( ! defined POSIX )
  function iargc ()
    integer :: iargc
    iargc = nag_iargc()
  end function iargc
# endif

  subroutine cl_getarg (k, arg, len, ierr)
    implicit none
    integer,       intent (in)  :: k
    character (*), intent (out) :: arg
    integer,       intent (in)  :: len
    integer,       intent (out) :: ierr
# ifdef POSIX
    interface
       subroutine pxfgetarg (m, buf, ilen, ierror)
         integer,       intent (in)  :: m
         character (*), intent (out) :: buf
         integer,       intent (in)  :: ilen
         integer,       intent (out) :: ierror
       end subroutine pxfgetarg
    end interface
# elif FCOMPILER != _GFORTRAN_ && FCOMPILER != _NAG_
    interface
       subroutine getarg (m, buf)
         integer,       intent (in)  :: m
         character (*), intent (out) :: buf
       end subroutine getarg
    end interface
# endif
# ifdef POSIX
    call pxfgetarg (k, arg, len, ierr)
# else
    call getarg (k, arg)
    ierr = 0
# endif
  end subroutine cl_getarg

# if FCOMPILER == _ALPHA_ && ( ! defined POSIX )
  subroutine ishell(string)
    ! execute the command in string in a system subshell
    character (*), intent(in) :: string
    call system(string)
  end subroutine ishell

  real function second()
    implicit none
    real (kind=4) :: etime,tarray(2)
    second = etime(tarray)
    return
  end function second
# endif

end module command_line
