program gs2

! main program which calls the run_gs2 subroutine

! make_lib is a compiler flag used if running with 
! trinity (coupled flux tube code)

# ifndef MAKE_LIB 
  use gs2_main, only: run_gs2
# endif

  implicit none

  call run_gs2

end program gs2
