#include "define.inc"
!>  This module is intended to be used for runtime tests
!!  which interrogate what is functional/what compile time
!!  options were enabled/disabled. 
!!
module runtime_tests

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Tests for compilers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function compiler_pgi()
    logical :: compiler_pgi
    compiler_pgi = .false.
#if FCOMPILER == _PGI_
    compiler_pgi = .true.
#endif
  end function compiler_pgi

  function get_compiler_name()
    character(len=9) :: get_compiler_name
    get_compiler_name='unknown'
#if FCOMPILER == _PGI_
    get_compiler_name='pgi'
#elif FCOMPILER == _INTEL_
    get_compiler_name='intel'
#elif FCOMPILER == _GFORTRAN_
    get_compiler_name='gfortran'
#elif FCOMPILER == _XL_
    get_compiler_name='xl'
#elif FCOMPILER == _NAG_
    get_compiler_name='nag'
#elif FCOMPILER == _CRAY_
    get_compiler_name='cray'
#elif FCOMPILER == _G95_
    get_compiler_name='g95'
#elif FCOMPILER == _PATHSCALE_
    get_compiler_name='pathscale'
#elif FCOMPILER == _LAHEY_
    get_compiler_name='lahey'
#elif FCOMPILER == _ABSOFT_
    get_compiler_name='absoft'
#elif FCOMPILER == _ALPHA_
    get_compiler_name='alpha'
#elif FCOMPILER == _SUN_
    get_compiler_name='sun'
#endif
  end function get_compiler_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Tests for svn info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>This function returns the output of svnversion 
  !It would be nice if we could strip of trailing text
  !so that we're just left with the integer revision.
  function get_svn_rev()
    character(len=8) :: get_svn_rev
#ifndef SVN_REV
#define SVN_REV "unknown"
#endif
    get_svn_rev=SVN_REV
  end function get_svn_rev

  !>This function returns true if the source code has
  !been modified relative to repo
  function get_svn_modified()
    logical :: get_svn_modified
    integer :: indx
#ifndef SVN_REV
#define SVN_REV "unknown"
#endif
    indx=index(SVN_REV,"M")
    if(indx.eq.0)then
       get_svn_modified=.false.
    else
       get_svn_modified=.true.
    endif
  end function get_svn_modified
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module runtime_tests
